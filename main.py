from typing import List, Dict, Any, Tuple
import pandas as pd
import geopandas as gpd
import shapely.geometry.multilinestring
import shapely.ops
import datetime
import re
import Geometry as gm
from Geometry import Vector2
import pandas as pd
import numpy as np

# container for output of our algorithm
output: List[Dict[str, Any]] = []

###############################
# LOAD ROAD LINE FEATURES
df_all_roads: gpd.geodataframe.GeoDataFrame = gpd.read_file("data.gdb", layer="NTWK_IRIS_Road_Network_20210305")

###############################
# LOAD ROAD INTERSECTION FEATURES
df_intersections: gpd.geodataframe.GeoDataFrame = gpd.read_file("data.gdb", layer="NTWK_Intersections_20210305")

##################################################################################################
# SELECT SPECIFIC ROADS FOR WHICH TO PERFORM THE ANALYSIS SINCE THIS TAKES FLIPPING FOREVER TO RUN

# take the ["ROAD"] number column from the database (which returns a pd.Series object) and call the unique() function on it which returns a NumPy.array object of unique values.
# We can assume it is a normal python List
list_of_state_road_numbers: List[str] = list(df_all_roads.loc[df_all_roads["NETWORK_TYPE"]=="State Road", "ROAD"].unique())

# override above list of all roads with smaller list of roads.
# list_of_state_road_numbers = pd.read_csv("list_of_road_numbers_to_process.csv")
# list_of_state_road_numbers = list_of_state_road_numbers[list_of_state_road_numbers["enabled"]==1]
# list_of_state_road_numbers = list_of_state_road_numbers["road"].to_list()


def direction_of_node(row: gpd.GeoSeries, start_end: str = "START"):
	"""
	:param row: represents a multi-polyline 'network element'
	:param start_end: tells the function to find the direction at the start or end of the line
	:returns: radians"""
	a: Vector2 = Vector2()
	b: Vector2 = Vector2()
	c: Vector2 = Vector2()
	
	geom = row.geometry
	if type(geom) is shapely.geometry.multilinestring.MultiLineString:
		geom = shapely.ops.linemerge(geom)
	if len(geom.coords) < 2:
		return None
	elif len(geom.coords) == 2:
		a, b = map(Vector2, geom.coords)
		return (b - a).direction()
	
	if start_end == "START":
		a, b, c = map(Vector2, geom.coords[:3])  # first three vertices in line converted to vectors and assigned to a, b, c
	elif start_end == "END":
		a, b, c = map(Vector2, geom.coords[-3:])  # last three vertices in line converted to vectors and assigned to a, b, c
	else:
		raise Exception("start_end parameter must have a value of 'START' or 'END'")
	
	ab = b - a
	bc = c - b
	weight_ab = ab.magnitude()
	weight_bc = bc.magnitude()
	weight_sum = weight_ab + weight_bc
	
	return gm.interpolate_angles(ab.direction(), bc.direction(), weight_bc / weight_sum)


# subsequent operations may rely on this sort order:
# TODO: not sure if they still do rely on this sort order..
df_all_roads.sort_values(["ROAD", "START_TRUE_DIST"], ascending=True, inplace=True)

counter_total: int = len(list_of_state_road_numbers)
counter: int = 1
print(datetime.datetime.now())
for current_road_number in list_of_state_road_numbers:
	print(f"{counter} of {counter_total} roads: {current_road_number}    {datetime.datetime.now()}")
	counter += 1
	
	# obtain a list of road segments which EXCLUDES the current_road we are considering.
	# this is not used immediately but is important for later
	all_roads_except_current_road: gpd.GeoDataFrame = df_all_roads[df_all_roads["ROAD"] != current_road_number]
	
	# filter the state roads data frame to get all the segments/rows that make up the current_road that we are looking at
	current_road_segments_L: gpd.GeoDataFrame = df_all_roads[
		  (df_all_roads["ROAD"] == current_road_number)
		& (df_all_roads["CWY"]  == "Left")
		#& (df_all_roads["slk_from"] )
		#& (df_all_roads["slk_to"] )
	]
	current_road_segments_R: gpd.GeoDataFrame = df_all_roads[(df_all_roads["ROAD"] == current_road_number) & (df_all_roads["CWY"] == "Right")]
	current_road_segments_S: gpd.GeoDataFrame = df_all_roads[(df_all_roads["ROAD"] == current_road_number) & (df_all_roads["CWY"] == "Single")]
	
	# from the above segments, extract the list of all nodes on the current road
	# using set() instead of list() to automatically avoid duplicates
	current_road_nodes_L: set = set()
	current_road_nodes_R: set = set()
	current_road_nodes_S: set = set()
	for index, intersecting_segment in current_road_segments_L.iterrows():
		current_road_nodes_L.add(intersecting_segment["END_NODE_NO"])
		current_road_nodes_L.add(intersecting_segment["START_NODE_NO"])
	
	for index, intersecting_segment in current_road_segments_R.iterrows():
		current_road_nodes_R.add(intersecting_segment["END_NODE_NO"])
		current_road_nodes_R.add(intersecting_segment["START_NODE_NO"])
	
	for index, intersecting_segment in current_road_segments_S.iterrows():
		current_road_nodes_S.add(intersecting_segment["END_NODE_NO"])
		current_road_nodes_S.add(intersecting_segment["START_NODE_NO"])
	
	print(f"    {current_road_number} has the following number of nodes:           L {len(current_road_nodes_L)}, R {len(current_road_nodes_R)}, S {len(current_road_nodes_S)}")
	
	# not all of the nodes we just collected from the current_road correspond to intersections
	# obtain a list of intersections which have matching node numbers; to do this we use the Pandas.Series.isin() function which tests if each value in a column is a member of a set
	# for some reason the "NODE_NAME" of an intersection is actually its number
	intersections_on_road_L: gpd.GeoDataFrame = df_intersections[df_intersections["NODE_NAME"].astype(str).isin(current_road_nodes_L)]
	intersections_on_road_R: gpd.GeoDataFrame = df_intersections[df_intersections["NODE_NAME"].astype(str).isin(current_road_nodes_R)]
	intersections_on_road_S: gpd.GeoDataFrame = df_intersections[df_intersections["NODE_NAME"].astype(str).isin(current_road_nodes_S)]
	print(f"    Of these, the following number are intersections: L {len(intersections_on_road_L.index)}, R {len(intersections_on_road_R.index)}, S {len(intersections_on_road_S.index)}")
	
	
	# TODO: There are some sneaky nodes that are also secretly intersections with the current_road but do not share a node number with the current_road
	#  this is a massively expensive operation... but only has to be done once i guess :/
	
	# this function allows us to simplify the next for loops by executing three consecutive for loops in the background where each loop goes through a different dataframe
	def compound_iterator_generator_1() -> Tuple[str, gpd.GeoDataFrame, gpd.GeoSeries]:
		for index, l_current_intersection in intersections_on_road_L.iterrows():
			yield "L", current_road_segments_L, l_current_intersection
		for index, l_current_intersection in intersections_on_road_R.iterrows():
			yield "R", current_road_segments_R, l_current_intersection
		for index, l_current_intersection in intersections_on_road_S.iterrows():
			yield "S", current_road_segments_S, l_current_intersection
	
	
	# compile list of road numbers of all roads intersecting current_road... we will prevent these from being part of the "continuing_roads" output column
	intersecting_road_numbers = set()
	for current_road_cway, current_road_segments, current_intersection in compound_iterator_generator_1():
		current_intersection_node_number = str(current_intersection["NODE_NAME"])
		
		segments_starting_or_ending_at_current_node: pd.DataFrame = all_roads_except_current_road[
			(all_roads_except_current_road["START_NODE_NO"] == current_intersection_node_number) |
			(all_roads_except_current_road["END_NODE_NO"] == current_intersection_node_number)
		]
		intersecting_road_numbers.update(segments_starting_or_ending_at_current_node["ROAD"].to_list())
	
	print(f"    {current_road_number} intersects the following roads: {', '.join([str(x) for x in intersecting_road_numbers])}")
	print("    Processing +Intersections and -Intersecting Roads:", end='')
	orphaned_intersections = []
	# Loop through all intersections on the current_road again
	for current_road_cway, current_road_segments, current_intersection in compound_iterator_generator_1():
		print("+", end='')
		current_intersection_node_number = current_intersection["NODE_NAME"]
		node_lng, node_lat = current_intersection.geometry.coords[0]
		# obtain segments that are connected to this node, excluding segments that are part of the current_road
		# this is done by filtering df_all_roads_except_current_road which was prepared at the top of the main loop
		# in order to preserve directionality; we must extract two lists, one of segments starting from this node, and one of segments ending at this node.
		# This way we get to know the SLK direction of a road for free
		segments_starting_from_current_node: pd.DataFrame = all_roads_except_current_road[all_roads_except_current_road["START_NODE_NO"] == current_intersection_node_number]
		segments_ending_at_current_node: pd.DataFrame = all_roads_except_current_road[all_roads_except_current_road["END_NODE_NO"] == current_intersection_node_number]
		
		# get the segments of the current road (on the current_road_cway) that connect to the current node;
		# unless the data is very messed up we should get only zero or one rows on the following two queries.
		current_road_segments_starting_at_node = current_road_segments[current_road_segments["START_NODE_NO"] == current_intersection_node_number]
		current_road_segments_ending_at_node = current_road_segments[current_road_segments["END_NODE_NO"] == current_intersection_node_number]
		road_network_element_starting_at_node = None
		road_network_element_ending_at_node = None
		road_true_dist_at_element_starting_from_node = ""
		road_true_dist_at_element_ending_at_node = ""
		road_slk_at_element_starting_from_node = ""
		road_slk_at_element_ending_at_node = ""
		
		direction_at_starting = None
		direction_at_ending = None
		
		# TODO: i assume that we have guaranteed that current_road_segments_starting_at_node and current_road_segments_ending_at_node contain either 0 or 1 records
		#  since we are looping through each carriageway in turn. I have added the following assert statements to guarantee this fact
		assert len(current_road_segments_starting_at_node.index) in {0, 1}
		assert len(current_road_segments_ending_at_node.index) in {0, 1}
		if not current_road_segments_starting_at_node.empty:
			road_network_element_starting_at_node = current_road_segments_starting_at_node.iloc[0]["NETWORK_ELEMENT"]
			road_true_dist_at_element_starting_from_node = current_road_segments_starting_at_node.iloc[0]["START_TRUE_DIST"]
			road_slk_at_element_starting_from_node = current_road_segments_starting_at_node.iloc[0]["START_SLK"]
			direction_at_starting = direction_of_node(current_road_segments_starting_at_node.iloc[0], "START")
		if not current_road_segments_ending_at_node.empty:
			road_network_element_ending_at_node = current_road_segments_ending_at_node.iloc[0]["NETWORK_ELEMENT"]
			road_true_dist_at_element_ending_at_node = current_road_segments_ending_at_node.iloc[0]["END_TRUE_DIST"]
			road_slk_at_element_ending_at_node = current_road_segments_ending_at_node.iloc[0]["END_SLK"]
			# TODO: I have corrected what i believe to be an error by changing the line below; The old line is commented out with the direction found at the "START"
			#  the new line finds the direction at "END".
			# direction_at_ending = direction_of_node(current_road_segments_ending_at_node.iloc[0], "START")
			direction_at_ending = direction_of_node(current_road_segments_ending_at_node.iloc[0], "END")
		
		# sanitise road_true_dist_first
		road_true_dist_first = road_true_dist_at_element_starting_from_node
		if type(road_true_dist_at_element_starting_from_node) is not np.float64:
			road_true_dist_first = road_true_dist_at_element_ending_at_node
			if type(road_true_dist_at_element_ending_at_node) is not np.float64:
				road_true_dist_first = "None"
		
		# sanitise road_slk_first
		road_slk_first = road_slk_at_element_starting_from_node
		if type(road_slk_at_element_starting_from_node) is not np.float64:
			road_slk_first = road_slk_at_element_ending_at_node
			if type(road_slk_at_element_ending_at_node) is not np.float64:
				road_slk_first = "None"
		
		# determine the average direction of the current road;
		if direction_at_starting is not None and direction_at_ending is not None:
			current_road_direction = gm.interpolate_angles(direction_at_starting, direction_at_ending, 0.5)
		elif direction_at_starting is not None:
			current_road_direction = direction_at_starting
		elif direction_at_ending is not None:
			current_road_direction = direction_at_ending
		else:
			current_road_direction = None
		
		if len(segments_ending_at_current_node.index) + len(segments_starting_from_current_node.index) == 0:
			print("0", end='')
			orphaned_intersections.append(current_intersection.to_list())
		
		
		def compound_iterator_generator_3():
			for index, intersecting_segment in segments_starting_from_current_node.iterrows():
				yield intersecting_segment, "START"
			for index, intersecting_segment in segments_ending_at_current_node.iterrows():
				yield intersecting_segment, "END"
		
		
		for intersecting_segment, startend in compound_iterator_generator_3():
			print("-", end='')
			# find the node at the other end of the intersecting road;
			# then find the list of road segments connected to that continuing_node
			
			# TODO: This assumes there are no intersecting roads along the way; ie this is a simple freeway ramp.
			#  we skip all the way to the other end of the intersecting road and record whats there, missing all the intermediat intersections.
			#  we could apply further smarts to only get the continuing roads we really want to know about; but thats very hardddd i think
			
			# first get all the segments of the intersecting road;
			intersecting_road_segments = df_all_roads[df_all_roads["ROAD"] == intersecting_segment["ROAD"]]
			
			# find the node at the other side of the intersecting road be choosing the first or last node of the first or last segment respectively
			#  TODO: Note that we just grab the last row from the sorted list of 'intersecting road' segments, we dont pay attention to which carriageway of the intersecting road arrives at this node.
			#   We just assume both L and R of the intersecting road finish at the same node. This works fine for ramps which are Single carriageway anyway.
			#  also determine direction of intersecting road (increasing slk), and the direction away from the current road (from current_node outward);
			node_at_other_side_of_intersecting_road = None
			intersecting_road_direction = None
			intersecting_road_direction_away = None
			if startend == "START":
				# find last node of intersecting road;
				node_at_other_side_of_intersecting_road = intersecting_road_segments.iloc[-1]["END_NODE_NO"]
				intersecting_road_direction = direction_of_node(intersecting_road_segments.iloc[0], "START")
				intersecting_road_direction_away = intersecting_road_direction
			elif startend == "END":
				# find first node of intersecting road;
				node_at_other_side_of_intersecting_road = intersecting_road_segments.iloc[0]["START_NODE_NO"]
				intersecting_road_direction = direction_of_node(intersecting_road_segments.iloc[-1], "END")
				intersecting_road_direction_away = gm.opposite_angle(intersecting_road_direction)
			
			# Find the direction of the intersecting road, leading away from the current node, relative to the direction of the current road increasing SLK
			intersecting_road_relative_direction = None
			if intersecting_road_direction is not None and current_road_direction is not None:
				intersecting_road_relative_direction = gm.angle_difference(intersecting_road_direction_away, current_road_direction)
			
			# then find all roads connected to this node; (don't include roads that are also primary connections to current_road)
			other_road_segments = df_all_roads[df_all_roads["ROAD"] != intersecting_segment["ROAD"]]
			road_segments_at_the_other_side_of_intersecting_road = other_road_segments[
				((other_road_segments["START_NODE_NO"] == node_at_other_side_of_intersecting_road) |
				 (other_road_segments["END_NODE_NO"] == node_at_other_side_of_intersecting_road)) &
				~(other_road_segments["ROAD"].isin(intersecting_road_numbers))
			]
			
			# compile the names of these continuing roads into a list
			list_of_continuing_roads = ""
			if not road_segments_at_the_other_side_of_intersecting_road.empty:
				list_of_continuing_roads = set()
				for index, other_road_segment in road_segments_at_the_other_side_of_intersecting_road.iterrows():
					list_of_continuing_roads.add(f"{other_road_segment['ROAD']} -- {other_road_segment['ROAD_NAME']}")
				list_of_continuing_roads = " ; ".join(list_of_continuing_roads)
			
			# Make some measurements of the intersecting road
			# TODO: Note that we are measuring the entire road, both carriageways. Could be improved by starting from intersecting segment and paying attention to the carriageway. i think both these measures might a bit pointless.
			length_of_intersecting_road = float(intersecting_road_segments.iloc[-1]["END_TRUE_DIST"]) - float(intersecting_road_segments.iloc[0]["START_TRUE_DIST"])
			number_of_intersecting_road_segments = len(intersecting_road_segments.index)
			
			# Apply a regex to decide if we have a ramp/slip road type name
			intersecting_road_looks_like_ramp = "RAMP" if re.search(r"( off to )|( on to )|( slip)", intersecting_segment["ROAD_NAME"], flags=re.IGNORECASE) else ""
			if intersecting_road_looks_like_ramp == "":
				if number_of_intersecting_road_segments <= 3 and length_of_intersecting_road < 2:
					if re.search(r" to ", intersecting_segment["ROAD_NAME"], re.IGNORECASE) or re.search(r" ramp", intersecting_segment["ROAD_NAME"], re.IGNORECASE):
						intersecting_road_looks_like_ramp = "RAMP"
			
			if intersecting_road_looks_like_ramp == "RAMP":
				if startend == "START":
					intersecting_road_looks_like_ramp = "OFF RAMP"
				else:
					intersecting_road_looks_like_ramp = "ON RAMP"
			
			# guess
			
			# finally; populate our output table with everything we know;
			output.append({
				"node":                                         current_intersection_node_number,
				"road_no":                                      current_road_number,
				"road_cwy":                                     current_road_cway,
				"road_network_element_starting_at_node":        road_network_element_starting_at_node,
				"road_network_element_ending_at_node":          road_network_element_ending_at_node,

				"road_true_dist_at_element_starting_from_node": road_true_dist_at_element_starting_from_node,
				"road_true_dist_at_element_ending_at_node":     road_true_dist_at_element_ending_at_node,
				"road_true_dist_first":                         road_true_dist_first,
				
				"road_slk_at_element_starting_from_node":       road_slk_at_element_starting_from_node,
				"road_slk_at_element_ending_at_node":           road_slk_at_element_ending_at_node,
				"road_slk_first":                               road_slk_first,
				"road_direction_deg":                           gm.radians_to_degrees(current_road_direction),
				"intersecting_direction":                       gm.radians_to_degrees(intersecting_road_direction),
				"intersecting_road_away_direction":             gm.radians_to_degrees(intersecting_road_direction_away),
				"intersecting_road_relative_direction":         gm.radians_to_degrees(intersecting_road_relative_direction),
				"intersecting_road_relative_LR":                "L" if intersecting_road_relative_direction < 0 else "R",
				"intersecting_network_element":                 intersecting_segment["NETWORK_ELEMENT"],
				"intersecting_road_cwy":                        intersecting_segment["CWY"][0],
				"intersecting_no":                              intersecting_segment["ROAD"],
				"intersecting_name":                            intersecting_segment["ROAD_NAME"],
				"intersecting_road_length":                     length_of_intersecting_road,
				"intersecting_road_segments":                   number_of_intersecting_road_segments,
				"intersecting_road_is_ramp":                    intersecting_road_looks_like_ramp,
				"intersecting_road_extremity_at_node":          startend,
				"continuing_node":                              node_at_other_side_of_intersecting_road,
				"continuing_roads":                             list_of_continuing_roads,
				"lng":                                          node_lng,
				"lat":                                          node_lat
			})
		print(" ", end='')
	print(" ")
	
	if orphaned_intersections:
		print("    Nodes that are intersections but not connected to other roads: ")
		for index, item in enumerate(orphaned_intersections):
			print("        " + str(item))

print("\r\n======= PHASE 2 ===========")
print("converting output to dataframe")
df = pd.DataFrame(output)

df.sort_values(["road_no", "road_true_dist_at_element_starting_from_node"], ascending=True, inplace=True)


def make_agg_join(joinstr: str = ''):
	def inner_agg_join(series: pd.Series):
		s1 = series.fillna('')
		return joinstr.join(s1[s1 != ""].to_list())
	
	return inner_agg_join


def make_aggregator_with_default(groupby, agg_dict: dict, default='first'):
	def inner(df: pd.DataFrame):
		dfg = df.groupby(groupby)
		agg_dict_inner = dict.fromkeys(df, default)
		for k, v in agg_dict.items():
			agg_dict_inner[k] = v
		dfr = dfg.agg(agg_dict_inner)
		dfr = dfr.reset_index(drop=True)
		return dfr
	
	return inner


print("Grouping by node and consolidating by grouping by intersecting_network_element, selecting the first of each column but joining the values of road_cwy, road_network_element_starting_at_node and road_network_element_ending_at_node")
df = df.groupby("node").apply(
	make_aggregator_with_default(
		groupby="intersecting_network_element",
		agg_dict={
			"road_cwy":                              ''.join,
			"road_network_element_starting_at_node": make_agg_join(";"),
			"road_network_element_ending_at_node":   make_agg_join(";"),
		},
		default='first'
	)
).reset_index(drop=True)

print("Sorting by intersecting_road_cwy to make the next grouping happen in a nice order")
print("Grouping by node then by intersecting_no then by intersecting_road_extremity_at_node: Then selecting the first of each column but joining the values of intersecting_road_cwy")
df = df.sort_values("intersecting_road_cwy", ascending=True)
df = df.groupby("node").apply(
	make_aggregator_with_default(
		groupby=[
			"intersecting_no",
			"intersecting_road_extremity_at_node"
		],
		agg_dict={
			"intersecting_road_cwy": ''.join,
			# "intersecting_network_element": ';'.join,
			# "continuing_node": lambda item: ';'.join(item.astype(str))
		},
		default='first'
	)
).reset_index(drop=True)
print("final sort;")
df = df.sort_values(["road_no", "road_true_dist_first", "road_cwy"])

df.to_csv(f"output/gis extract - intersections with direction {datetime.datetime.now().strftime('%Y%m%d')}.csv", index_label="index")

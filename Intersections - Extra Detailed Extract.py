"""
Nicholas Archer 2020-04-16

This script extracts a list of intersections along each state road.
It is similar to the 'intersections' extract from IRIS Reporting Center, but extracts more information

It loops over every state road
	└┬─then over every node on that road
	 └───then over each row (network element) in the road database also touching that node (exclusing network elements forming part of the primary state road that we are busy looping over)
			Then it outputs a row
			sometimes if the primary road is dual carriageway, it will share 1 node for both carriageways
			sometimes there will be a node for each carriageway
			see diagram below which shows two of the many possible ways intersections can work, resulting in 3 output rows.


        >==== contining road A =====(CONTINUING NODE)===== continuing road A ===>                  (CONTINUING NODE)
                                      ||                                                            ||
                                      (output row 1)                                                (output row 3)
                                      ||                                                            ||
                                      (NODE) (* intermediate nodes not captured)                    ||
                                      ||                                                            ||
                                      ||                                                            ||
>===== primary road H001 (L) ========= ========================================================== (NODE) ======================================>
                                    (NODE)
>===== primary road H001 (R) ========= ========================================================== (NODE) ======================================>
                                      ||
                                      ||
                                      (output row2)
                                      ||
        <== contining road B =====(CONTINUING NODE)===== continuing road C ==>



Finally each row has information on the 'continuing node' or the node you would get to if you followed the intersecting road all the way to its far end.

	- primary road
		- road number
		- road slk
		- carriageway (L, R, LR, LRS, S) where LR means 'both' and 'LRS' occurs where the primary road changes from dual cway to single cway at the same node as the intersection
		- direction of primary road
		
	- intersecting road
		- intersecting road no
		- cardinal direction
		- intersecting carriageway (L, R, LR, S)
		- direction relative to the primary road in degrees
		- direction realative to the primary road (L, R) (for convienience; can be computed from the relative degrees)
		- length # TODO: note this is currently the total for the whole of the intersecting road, not as it should be; from the intersection with the primary road to the relevant start node / end node.
		- number of segments  # TODO: note this is currently the total for the whole of the intersecting road, not as it should be; from the intersection with the primary road to the relevant start node / end node
		- a guess if this is an ON RAMP or OFF RAMP
		  note that ramp direction is guessed based on the direction of increasing SLK of the intersecting road.
		  ie; the leach highway onramp, may be lablled as an offramp for output rows where the primary road is albany highway
		
	- continuing node
	  (node at the furthest extent of intersecting road, either the start node or end node of the intersecting road)
	  (this is useful for determining what road a freeway ramp connects to. It is a bit non-sensical to read this informaiton for most intersections on albany highway)
		- node number
		- continuing roads (other roads leading away from this node, not including the intersecting road reffered to by this row, and not including any other intersecting road which connects back to the primary road
		  listed in a single field with the format '<ROAD_NO> - <ROAD NAME> ; <ROAD_NO> - <ROAD NAME> ; ...'
		 




The extraction contains a lot of redundant rows because it outputs a row for each carriageway of each road, and each carriageway of each intersecting road

data.gdb includes informaiton extracted from the RIME Spatial Server

	NTWK_Intersections_20200424
		(coppied from RIME_Spatial.RIME.NTWK_Intersections)
		(one point feature per row per intersection)
		OBJECTID NODE_NAME	NODE_DESCR NODE_TYPE	NO_NODE_ID

	NTWK_IRIS_Road_Network_20200424
		(coppied from RIME_Spatial.RIME.NTWK_IRIS_Road_Network)
		(each road is in multiple rows, where each row is a multi-polyline feature. Although Multi-Polylines are )
		OBJECTID	ROAD	ROAD_NAME	COMMON_USAGE_NAME	START_SLK	END_SLK	CWY	START_TRUE_DIST	END_TRUE_DIST	NETWORK_TYPE	RA_NO	RA_NAME	LG_NO	LG_NAME	START_NODE_NO	START_NODE_NAME	END_NODE_NO	END_NODE_NAME	DATUM_NE_ID	NM_BEGIN_MP	NM_END_MP	NETWORK_ELEMENT	ROUTE_NE_ID	GEOLOC_STLength__

General Notes about the Code blow
	in the code below a "segment" refers to a piece of a road between two nodes; ie one row of the RIME_Spatial.RIME.NTWK_IRIS_Road_Network
	in formal graph math i suppose this should be called an "edge"
	'Nodes' may not be intersections in the graph; they can happen in the middle of a road, or where a road changes from single to dual carriageway
	for this reason we need to refer to the RIME_Spatial.RIME.NTWK_Intersections table; this tells us where the actual intersections are so that we can ignore the others
	
	In an earlier version this extracted OBJECTID... it turns out this is useless. ArcGIS changes that whenever it feels like it. Not usefull as a unique reference.
	This version uses the "NETWORK_ELEMENT" field

"""

from typing import List, Dict, Any
import pandas as pd
import geopandas as gpd
import shapely.geometry.multilinestring
import shapely.ops
import datetime
import re
import Geometry as gm
from Geometry import Vector2

import numpy as np

# container for output of our algoritim
output: List[Dict[str, Any]] = []

###############################
# LOAD ROAD LINE FEATURES
df_all_roads: gpd.geodataframe.GeoDataFrame = gpd.read_file("data.gdb", layer="NTWK_IRIS_Road_Network_20200424")
df_state_roads: gpd.geodataframe.GeoDataFrame = df_all_roads[df_all_roads["ROAD"].str.startswith(("H", "M"), na=False)]
# TODO: previously we would filter the intersectinos here to only show "State Road Node"s but i think that may have been the cause of missing so many.
#  Revmoving that filter has caused big drop in speed.

###############################
# LOAD ROAD INTERSECTION FEATURES
df_intersections: gpd.geodataframe.GeoDataFrame = gpd.read_file("data.gdb", layer="NTWK_Intersections_20200424")


###############################
# LOAD ROAD INTERSECTION FEATURES


# take the ["ROAD"] number column from the database (which returns a pd.Series object) and call the unique() function on it which returns a NumPy.array object of unique values.
# We can assume it is a normal python List
# list_of_state_road_numbers: List[str] = df_state_roads["ROAD"].unique()
# override above list of all roads with smaller list of roads.
list_of_state_road_numbers = [
	'H036',
	"H574",
	"H021",
	"H001",
	"H015",
	"H016",
	"H026",
	"H005",
	"H027",
	"H052",
	"H013",
	"H002",
	"H029"
]
list_of_state_road_numbers = ['H036', "H016"]


def direction_of_node(row: gpd.GeoSeries, start_end: str = "START"):
	""":return: radians"""
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
		a, b, c = map(Vector2, geom.coords[:3])  # first three verticies in line converted to vectors and assigned to a, b, c
	elif start_end == "END":
		a, b, c = map(Vector2, geom.coords[-3:])  # last three verticies in line converted to vectors and assigned to a, b, c
	else:
		raise Exception("start_end parameter must have a value of 'START' or 'END'")
	
	ab = b - a
	bc = c - b
	weight_ab = ab.magnitude()
	weight_bc = bc.magnitude()
	weight_sum = weight_ab + weight_bc
	
	return gm.interpolate_angles(ab.direction(), bc.direction(), weight_bc/weight_sum)


# subsequent opperations may rely on this sort order: # TODO: not sure if they still do rely on this sort order..
df_all_roads.sort_values(["ROAD", "START_TRUE_DIST"], ascending=True, inplace=True)

counter_total: int = len(list_of_state_road_numbers)
counter: int = 1
print(datetime.datetime.now())
for current_road_number in list_of_state_road_numbers:
	print(f"{counter} of {counter_total} roads: {current_road_number}    {datetime.datetime.now()}")
	counter += 1

	# obtain a list of road segments which DISCLUDES the current_road we are considering.
	# this is not used immediately but is important for later
	all_roads_except_current_road:gpd.GeoDataFrame = df_all_roads[df_all_roads["ROAD"] != current_road_number]
	
	# filter the state roads dataframe to get all the segments/rows that make up the current_road that we are looking at
	current_road_segments_L: gpd.GeoDataFrame = df_state_roads[(df_state_roads["ROAD"] == current_road_number) & (df_state_roads["CWY"] == "Left")]
	current_road_segments_R: gpd.GeoDataFrame = df_state_roads[(df_state_roads["ROAD"] == current_road_number) & (df_state_roads["CWY"] == "Right")]
	current_road_segments_S: gpd.GeoDataFrame = df_state_roads[(df_state_roads["ROAD"] == current_road_number) & (df_state_roads["CWY"] == "Single")]
	
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
	
	# TODO: obtain sneaky nodes that are also secret intersections but dont share a node number with the current_road
	#  this is a massively expensive opperation... but only has to be done once i guess :/
	
	# this function allows us to simplify the next for loops by executing three consecutive for loops in the background where each loop goes through a different dataframe
	def compound_iterator_generator_1() -> (str, gpd.GeoDataFrame, gpd.GeoSeries):
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
		
		segments_startfing_or_ending_at_current_node: pd.DataFrame = all_roads_except_current_road[
			(all_roads_except_current_road["START_NODE_NO"] == current_intersection_node_number) |
			(all_roads_except_current_road["END_NODE_NO"] == current_intersection_node_number)
		]
		intersecting_road_numbers.update(segments_startfing_or_ending_at_current_node["ROAD"].to_list())
	
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
		
		direction_at_starting = None
		direction_at_ending = None
		
		if not current_road_segments_starting_at_node.empty:
			road_network_element_starting_at_node = current_road_segments_starting_at_node.iloc[0]["NETWORK_ELEMENT"]
			road_true_dist_at_element_starting_from_node = current_road_segments_starting_at_node.iloc[0]["START_TRUE_DIST"]  # TODO: assumes true may have gap or overlap at node
			direction_at_starting = direction_of_node(current_road_segments_starting_at_node.iloc[0], "START")
		if not current_road_segments_ending_at_node.empty:
			road_network_element_ending_at_node = current_road_segments_ending_at_node.iloc[0]["NETWORK_ELEMENT"]
			road_true_dist_at_element_ending_at_node = current_road_segments_ending_at_node.iloc[0]["END_TRUE_DIST"]  # TODO: assumes true may have gap or overlap at node
			direction_at_ending = direction_of_node(current_road_segments_ending_at_node.iloc[0], "START")
		
		# sanitise road_true_dist_first
		road_true_dist_first = road_true_dist_at_element_starting_from_node
		if type(road_true_dist_at_element_starting_from_node) is not np.float64:
			road_true_dist_first = road_true_dist_at_element_ending_at_node
			if type(road_true_dist_at_element_ending_at_node) is not np.float64:
				road_true_dist_first = "None"
		
		
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
			
			# first get all the segements of the intersecting road;
			intersecting_road_segments = df_all_roads[df_all_roads["ROAD"] == intersecting_segment["ROAD"]]
			
			
			
			# find the node at the other side of the intersecing road be choosing the first or last node of the first or last segment respectively
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
			
			
			# then find all roads connected to this node; (dont include roads that are also primary connections to current_road
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
			# TODO: Note that we are measuring the entire road, both carriageways. Could be improeved by starting from intersecting segment and paying attention to the carriageway. i think both these measures might a bit pointless.
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
				"road_true_dist_first":							road_true_dist_first,
				"road_direction_deg":							gm.radians_to_degrees(current_road_direction),
				"intersecting_direction":						gm.radians_to_degrees(intersecting_road_direction),
				"intersecting_road_away_direction":				gm.radians_to_degrees(intersecting_road_direction_away),
				"intersecting_road_relative_direction":			gm.radians_to_degrees(intersecting_road_relative_direction),
				"intersecting_road_relative_LR":				"L" if intersecting_road_relative_direction < 0 else "R",
				"intersecting_network_element":                 intersecting_segment["NETWORK_ELEMENT"],
				"intersecting_road_cwy":                        intersecting_segment["CWY"][0],
				"intersecting_no":                              intersecting_segment["ROAD"],
				"intersecting_name":                            intersecting_segment["ROAD_NAME"],
				"intersecting_road_length":                     length_of_intersecting_road,
				"intersecting_road_segments":                   number_of_intersecting_road_segments,
				"intersecting_road_is_ramp":                     intersecting_road_looks_like_ramp,
				"intersecting_road_extremity_at_node":          startend,
				"continuing_node":                              node_at_other_side_of_intersecting_road,
				"continuing_roads":                             list_of_continuing_roads,
				"lng":											node_lng,
				"lat":											node_lat
			})
		print(" ", end='')
	print(" ")
	
	if orphaned_intersections:
		print("    Nodes that are intersections but not connected to other roads: ")
		for index, item in enumerate(orphaned_intersections):
			print("        "+str(item))


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


print("Grouping by node and consolidating by grouping by intersecting_network_element, selecting the first of each column but joining the valies of road_cwy, road_network_element_starting_at_node and road_network_element_ending_at_node")
df = df.groupby("node").apply(
	make_aggregator_with_default(
		groupby="intersecting_network_element",
		agg_dict={
			"road_cwy": ''.join,
			"road_network_element_starting_at_node": make_agg_join(";"),
			"road_network_element_ending_at_node": make_agg_join(";"),
		},
		default='first'
	)
).reset_index(drop=True)
print("Sorting by intersecting_road_cwy to make the next grouping happen in a nice order")
print("Grouping by node then by intersecting_no then by intersecting_road_extremity_at_node: Then selecting the first of each column but joining the valies of intersecting_road_cwy")
df = df.sort_values("intersecting_road_cwy", ascending=True)
df = df.groupby("node").apply(
	make_aggregator_with_default(
		groupby=[
			"intersecting_no",
			"intersecting_road_extremity_at_node"
		],
		agg_dict={
			"intersecting_road_cwy": ''.join,
			#"intersecting_network_element": ';'.join,
			#"continuing_node": lambda item: ';'.join(item.astype(str))
		},
		default='first'
	)
).reset_index(drop=True)
print("final sort;")
df = df.sort_values(["road_no", "road_true_dist_first", "road_cwy"])

df.to_csv(f"output/gis extract - intersections with direction {datetime.datetime.now().strftime('%Y%m%d')}.csv", index_label="index")

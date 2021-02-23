# Advanced Intersection Information Extractor

# Description
This script extracts a list of intersections along each state road, but with extra information.

It intended to supplement the dataset available here: https://portal-mainroads.opendata.arcgis.com/datasets/intersections

Additional information is extracted about the direction of intersection legs relative to each other, and the connectivity of intersecting roads with the left and right carriageway of the primary road.

The output data is intended to be filtered road-by-road.
If we filter the output data by `road_no="H001"` we can expect a list of all intersections on Albany Highway.
If we do the same for Leach Highway `road_no="H012"` then we would also get all intersections for Leach Highway.
Both extracts would include the common intersection between Albany Highway and Leach Highway but from a different point of view.

## Define some types before explaining the algorithm
```python
from geopandas import GeoDataFrame, GeoSeries
road_network: GeoDataFrame
intersections: GeoDataFrame
row: GeoSeries  
road:str # road number
node:str # node number

road_rows: GeoDataFrame
road_row: GeoSeries

```
Note that `row` represents a record from `road_network`

Each `row` represents a segment of a road, from some `'START_SLK'` to `'END_SLK'`
Each row has an `'START_NODE_NO'` and an `'END_NODE_NO'`

## Explaining the algorithm:
- For each `road` in `unique(road_network['ROAD'])`
  - let `road_rows = rows in road_network where ['ROAD'] == road` 
  - For each `node` in `road_rows`
    - For each `row` in `road_network` where `row` shares `node` and `row['ROAD'] is not road`
	  - Take data from `road_rows`, `node`, and `row` etc
      - output the measures to table called `intermediate_rows`
- Aggregate `intermediate_rows` such that there is only one output row per intersecting road at each `node`
- Output a `.csv` file with the columns described in the table below.

Note:Sometimes if the `road` is a dual carriageway, it will share one `node` for both carriageways.
Sometimes there will be two nodes; one for each carriageway. 
The diagram below shows two of the many possible ways intersections can work.
This example would result in three output rows.

```
        >==== continuing road A =====(CONTINUING NODE)===== continuing road A ===>                  (CONTINUING NODE)
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
        <== continuing road B =====(CONTINUING NODE)===== continuing road C ==>
```


## Output Columns

TODO: Currently SLK is not extracted. Only true distance is output to the table. It would be better to extract both.

|Column|Variable|Describes|Description|
|------|--------|---------|-----------|
|node | `current_intersection_node_number` |Node| Unique node ID number |
|road_no | `current_road_number` | Primary Road| Road Number|
|road_cwy | `current_road_cway` |  Primary Road| Carriageway (L, R, LR, LRS, S) where LR means 'both' and 'LRS' occurs where the primary road changes from dual cway to single cway at the same node as the intersection|
|road_network_element_starting_at_node | `road_network_element_starting_at_node` | Primary Road| |
|road_network_element_ending_at_node | `road_network_element_ending_at_node` | Primary Road| |
|road_true_dist_at_element_starting_from_node | `road_true_dist_at_element_starting_from_node` | Primary Road| |
|road_true_dist_at_element_ending_at_node | `road_true_dist_at_element_ending_at_node` | Primary Road | |
|road_true_dist_first | `road_true_dist_first` | Primary Road|The TRUE_DIST_START field for the segment of the primary road starting from the current node|
|road_direction_deg | `gm.radians_to_degrees(current_road_direction)` | Primary Road | Cardinal direction of the primary road in degrees|
|intersecting_direction | `gm.radians_to_degrees(intersecting_road_direction)`| Primary Road | Intersecting road cardinal direction in degrees |
|intersecting_road_away_direction | `gm.radians_to_degrees(intersecting_road_direction_away)`|Primary Road| |
|intersecting_road_relative_direction | `gm.radians_to_degrees(intersecting_road_relative_direction)`|Intersecting Road | Direction (in degrees) relative to the primary road (facing increasing SLK, negative values mean Left, positive values mean Right)|
|intersecting_road_relative_LR | `"L" if intersecting_road_relative_direction < 0 else "R"`|Intersecting Road  | Direction relative to the primary road ('L', 'R') (included for convenience since it can be computed from the relative degrees)|
|intersecting_network_element | `intersecting_segment["NETWORK_ELEMENT"]` |Intersecting Road ||
|intersecting_road_cwy | `intersecting_segment["CWY"][0]` |Intersecting Road | Carriageway (L, R, LR, S)|
|intersecting_no | `intersecting_segment["ROAD"]` | Intersecting Road | Road Number|
|intersecting_name | `intersecting_segment["ROAD_NAME"]` |Intersecting Road | Road name|
|intersecting_road_length | `length_of_intersecting_road`|Intersecting Road  | Length of entire intersecting road. (TODO: fix. It should only measure be from intersection to relevant start/end node, not the entire intersecting road)|
|intersecting_road_segments | `number_of_intersecting_road_segments`|Intersecting Road  |Number of segments in entire intersecting road. (TODO: fix. It should only count be from intersection to relevant start/end node, not the entire intersecting road)|
|intersecting_road_is_ramp | `intersecting_road_looks_like_ramp`|Intersecting Road  | A rough guess whether the intersecting road is a freeway ramp. Could be `"ON RAMP"` or `"OFF RAMP"`. Note that ramp direction is guessed based on the direction of increasing SLK of the intersecting road.  ie; the Leach Highway on-ramp, may be labeled as an off-ramp for output rows where the primary road is Albany Highway|
|intersecting_road_extremity_at_node | `startend`|Intersecting Road  |Does the intersecting road _segment_ start or end at the node... Could be `"START"` or `"END"`. ie. is the intersecting road's SLK increasing or decreasing when traveling away from the node|
|continuing_node | `node_at_other_side_of_intersecting_road`|Continuing Node  | The node number of the node at the far end of the intersecting road. Can be used to look up the name of the 'continuing roads' of a freeway ramp|
|continuing_roads | `list_of_continuing_roads`|Continuing Node |The names of the roads which continue from the 'continuing node'. Listed in a single field with the format `'<ROAD_NO> - <ROAD NAME> ; <ROAD_NO> - <ROAD NAME> ; ...'`|
|lng | `node_lng` |Node|Location of the intersection with the primary road|
|lat | `node_lat` |Node|Location of the intersection with the primary road|
	



# Source Data
Source data is not included in the repo due to size.
Similar data is available in GeoJSON and other formats at<br> https://portal-mainroads.opendata.arcgis.com/

The intersection dataset that is publicly available seems to be missing the NO_NODE_ID field which links it with the road network dataset.
With some effort it could be reconstructed from this dataset<br> https://portal-mainroads.opendata.arcgis.com/datasets/082e88d12c894956945ef5bcee0b39e2_17

`data.gdb` includes the following data extracts:

- `NTWK_Intersections_20200424`<br>
  (one point feature per row, per intersection)
    - `OBJECTID`
    - `NODE_NAME` 
    - `NODE_DESCR` 
    - `NODE_TYPE`
    - `NO_NODE_ID`


- `NTWK_IRIS_Road_Network_20200424`<br>
  (each road is in multiple rows, where each row is a MultiLineString feature.
  In this dataset all MultiLineString features are always composed of a single LineString object.)
    - `OBJECTID`
    - `ROAD`
    - `ROAD_NAME`
    - `COMMON_USAGE_NAME`
    - `START_SLK`
    - `END_SLK`
    - `CWY`
    - `START_TRUE_DIST`
    - `END_TRUE_DIST`
    - `NETWORK_TYPE`
    - `RA_NO`
    - `RA_NAME`
    - `LG_NO`
    - `LG_NAME`
    - `START_NODE_NO`
    - `START_NODE_NAME`
    - `END_NODE_NO`
    - `END_NODE_NAME`
    - `DATUM_NE_ID`
    - `NM_BEGIN_MP`
    - `NM_END_MP`
    - `NETWORK_ELEMENT`
    - `ROUTE_NE_ID`
    - `GEOLOC_STLength__`
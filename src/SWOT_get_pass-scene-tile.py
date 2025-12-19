
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 08:35:48 2025

@author: alchrist
"""



## Answers the question, what SWOT frame is my gauge located within?
import zipfile
from pathlib import Path
import fiona
import geopandas as gpd
fiona.drvsupport.supported_drivers['KML'] = 'rw'
import os
import tempfile
import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
import subprocess 
from shapely.geometry import box

base_dir = Path(os.getcwd())
# gauge_folder = base_dir / 'Validation/gauges/'
orbit_path = base_dir / 'orbit_files' 
Path(orbit_path).mkdir(parents=True, exist_ok=True)




print('Provide a CSV with either a list of points (name, latitude, longitude) or a list of aois with columns aoi, minx, miny, maxx, maxy')
print('Coordinates must be in WGS84 (EPSG 4326)')
input_file = Path(str(input('Path to CSV file: \n')))
csv = pd.read_csv(input_file)
print(csv.head())

# csv = csv.iloc[53:54]
try:
    points = gpd.GeoDataFrame({'name':csv['name'],'geometry': gpd.points_from_xy(csv['longitude'], csv['latitude'])}, crs="EPSG:4326")
    name_id = 'name'

except:
    points = gpd.GeoDataFrame({'aoi':csv['aoi'],'geometry': [box(minx=row["minx"], miny=row["miny"], maxx=row["maxx"], maxy=row["maxy"]) for i,row in csv.iterrows()]}, crs="EPSG:4326")
    name_id = 'aoi'

## Cal Val 1-day Repeat SWOT Orbit
file = orbit_path / 'swot_beta_preval_coverage_20231204.kmz'
if os.path.isfile(file)==False:
    url = 'https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/web-misc/swot_mission_docs/swot_beta_preval_coverage_20231204.kmz'
    subprocess.run(['curl','-L',url,'-o',file])


temp_dir = tempfile.mkdtemp()
with zipfile.ZipFile(file, 'r') as kmz:
    kml_filepath = next((f for f in kmz.namelist() if f.lower().endswith(".kml")), None)
    kmz.extract(kml_filepath, path=temp_dir)
    kml_path = os.path.join(temp_dir, kml_filepath)
    fiona.drvsupport.supported_drivers['KML'] = 'rw'

    passes = []
    scenes = []
    names = []
    tiles = []

    scenes_gdf = gpd.read_file(kml_path, driver='KML', layer='scenes')
    points_within_scenes = gpd.sjoin(scenes_gdf,points , how='right', predicate='intersects')
    points_within_scenes = points_within_scenes[~np.isnan(points_within_scenes['index_left'])].reset_index(drop=True)
    
    tiles_gdf = gpd.read_file(kml_path, driver='KML', layer='tiles')
    points_within_tiles = gpd.sjoin(tiles_gdf,points , how='right', predicate='intersects')
    points_within_tiles = points_within_tiles[~np.isnan(points_within_tiles['index_left'])].reset_index(drop=True)
    if len(points_within_tiles)>0:                
        for i in range(0,len(points_within_tiles)):
            names.append(points_within_tiles.iloc[i][name_id])


            passs = points_within_tiles.iloc[i]['pass_num'].zfill(3)
            print('Pass: %s' %(passs)) 
            passes.append(passs)
            tile = points_within_tiles.iloc[i]['tile_num'].zfill(3)+points_within_tiles.iloc[i]['tile_side']
            tiles.append(tile.zfill(3))
            print('Tile: %s' %(tile))
        
            scene = str(int(np.ceil(int(tile[:-1])/2))).zfill(3)+'F'

            print('Scene: %s' %(scene))
            scenes.append(scene)
    
    calval = pd.DataFrame({'pass':passes,'tiles':tiles,'scenes':scenes,'filename':names})
    # calval.to_csv(str(input_file).replace('.csv','_SWOT_CalValOrbit.csv'))




## Science 21-day Repeat SWOT Orbit
file = orbit_path / 'swot_science_coverage_20240319.kmz'
if os.path.isfile(file)==False:
    url = 'https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/web-misc/swot_mission_docs/swot_science_coverage_20240319.kmz'
    subprocess.run(['curl','-L',url,'-o',file])

temp_dir = tempfile.mkdtemp()
with zipfile.ZipFile(file, 'r') as kmz:
    kml_filepath = next((f for f in kmz.namelist() if f.lower().endswith(".kml")), None)
    kmz.extract(kml_filepath, path=temp_dir)
    kml_path = os.path.join(temp_dir, kml_filepath)
    fiona.drvsupport.supported_drivers['KML'] = 'rw'

    passes = []
    scenes = []
    names = []
    tiles = []
    for layer in fiona.listlayers(kml_path):
        if layer[0]=='P':
            print(layer)
            polygon_gdf = gpd.read_file(kml_path, driver='KML', layer=layer)
            points_within_polygons = gpd.sjoin(polygon_gdf,points , how='right', predicate='intersects')
            points_within_polygons = points_within_polygons[~np.isnan(points_within_polygons['index_left'])]
            if len(points_within_polygons)>0:
                print(layer)
                
                for i in range(0,len(points_within_polygons)):
                    soup = BeautifulSoup(points_within_polygons.iloc[i]['description'], "html.parser")
                    categories = soup('b', text=lambda text: text and text.endswith(":"))
                    if len(categories)>0:
                    
                        passes.append(categories[0].next_sibling.strip('" \n').split(", ")[0])
                        tiles.append(categories[1].next_sibling.strip('" \n').split(", ")[0])
                        scenes.append(categories[2].next_sibling.strip('" \n').split(", ")[0])
                        names.append(points_within_polygons.iloc[i][name_id])

    science = pd.DataFrame({'pass':passes,'tiles':tiles,'scenes':scenes,'filename':names})
    # science.to_csv(str(input_file).replace('.csv','_SWOT_ScienceOrbit.csv'))
  
    
  
scenes = []
tiles = []
passes = []
for row in csv[name_id]:
    
    this_science = science[science['filename']==row]
    this_tile = []
    this_scene = []
    this_pass = []
    for i in this_science.index:
        this_tile.append(this_science.loc[i]['tiles'])
        this_scene.append(this_science.loc[i]['scenes'])
        this_pass.append(this_science.loc[i]['pass'])

    
    this_calval = calval[calval['filename']==row]
    this_tile_c = []
    this_scene_c = []
    this_pass_c = []
    for i in this_calval.index:
        this_tile.append(this_calval.loc[i]['tiles'])
        this_scene.append(this_calval.loc[i]['scenes'])
        this_pass.append(this_calval.loc[i]['pass'])
    tiles.append(this_tile)
    scenes.append(this_scene)
    passes.append(this_pass)

csv['pass'] = passes
csv['scene'] = tiles
csv['tile'] = scenes
csv['firstdate'] = '20230101'
csv['lastdate'] = '20250515'


csv.to_csv(str(input_file).replace('.csv','_SWOTacquisitions.csv'))





#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 14:40:43 2026

@author: alchrist
"""

## Import Packages
import os
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import glob
import argparse

import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import ast
import sys
import fnmatch 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import earthaccess
from pyproj import CRS
import rioxarray as rio
from shapely.geometry import box
from collections import defaultdict

## Import custom packages 
from swot_utils import (haversine, 
                        add_WSE_HRL2,
                        calculate_delh_MeantTide_minus_TideFree_geoidheights, 
                        convert_geoid_fromMeanTide_toTideFree,
                        )   
from swot_plotting import plot1variable_HR_Raster                        

##### Get command inputs
parser = argparse.ArgumentParser(
    description='''Script to extract SWOT HR data within an AOI ''')

parser.add_argument('--interactive', dest='interactive', action="store_true", help='if interactive is chosen, you will be prompted to enter information need for search')
parser.add_argument('--aoi', dest='aoi', type=str, help='If not interactive, you must provide the AOI name')
parser.add_argument('--limited', dest='limited', action="store_true", help='Limited will only save a few attributes, other will save all attributes')
parser.add_argument('--crop', dest='crop', action="store_true", help='Crop option will loop for *crop.shp in the AOI directory')
parser.add_argument('--i', dest='index_array', type=int, help='if running job array, what index number')
args = parser.parse_args()
if len(sys.argv) < 2:
    args.interactive = True


##############################################################################
## Set Directory
base_dir = Path(os.path.dirname(os.getcwd()))
print(f"Base Directory: {base_dir}")

##############################################################################
## Get bounding information from CSV
try:
    bounding_LUT = pd.read_csv(base_dir / 'aoi_template.csv')
    aois = list(bounding_LUT['aoi'])
except Exception:
    bounding_LUT = pd.DataFrame([['default', np.nan, np.nan, np.nan, np.nan, 10, 15, -2, 2, 0, 1]], 
                                columns=['aoi', 'minx', 'miny', 'maxx', 'maxy', 'minssh', 'maxssh', 'minwse', 'maxwse', 'mindac', 'maxdac'])
    aois = []
    
if args.interactive:
    aoi = input('which AOI: %s ' % (aois))
else:
    aoi = args.aoi

## Get bounding areas of AOI
if aoi not in aois:
    LUT = bounding_LUT[bounding_LUT['aoi'] == 'default'].reset_index()
    print('what bounding box (lat/lon in decimal degrees)')
    LUT['minx'][0] = float(input('xmin: '))
    LUT['miny'][0] = float(input('ymin: '))
    LUT['maxx'][0] = float(input('xmax: '))
    LUT['maxy'][0] = float(input('ymax: '))
    print(LUT)
else:
    LUT = bounding_LUT[bounding_LUT['aoi'] == aoi].reset_index()

area = [LUT['minx'][0], LUT['miny'][0], LUT['maxx'][0], LUT['maxy'][0]]

## Correct any negative E-W coordinates
if area[0] < 0: 
    area2 = [360 + area[0], area[1], 360 + area[2], area[3]]
else:
    area2 = [area[0], area[1], area[2], area[3]]

print('Crop to Area: %s' % (area2))


##############################################################################
### SWOT Product Information
l2version = 'D'
if l2version == 'C':
    l2hr_prefix = 'SWOT_L2_HR_Raster'
    L2HR_folder = base_dir / 'Data' / (l2hr_prefix + '_2.0')
else:
    l2hr_prefix = 'SWOT_L2_HR_Raster' 
    L2HR_folder = base_dir / 'Data' / (l2hr_prefix + '_' + l2version)

l2_processing = 'P*%s*' % (l2version)
resolution = '100m'

### Create directory for L2 data
L2HR_fig_folder = L2HR_folder / aoi / 'plots'
Path(L2HR_fig_folder).mkdir(parents=True, exist_ok=True)

L2HR_output_folder = L2HR_folder / aoi / 'Extracted'
Path(L2HR_output_folder).mkdir(parents=True, exist_ok=True)

aoi_dir = base_dir / 'Examples' / aoi


##############################################################################
## Get tide-free geoid
geoid_name = 'geoid_egm2008_wgs84_tidefree_min1x1_20200109T140029_v100_fixedshift.nc'
geoid_file = glob.glob(str(base_dir / ('*/' + geoid_name)))

xvar = 'gauge_wse'

passes_tiles = ['%s_%s' % (i, j) for i, j in zip(ast.literal_eval(LUT['pass'][0]), ast.literal_eval(LUT['tile'][0]))] 

all_nc_files = glob.glob(str(L2HR_folder / '*_100m_*.nc'))
all_nc_files = [i for i in all_nc_files if '_'.join(i.split('/')[-1].split('_')[11:13]) in passes_tiles]

pass_tile_uni = np.unique(['_'.join(i.split('/')[-1].split('_')[11:13]) for i in all_nc_files])
all_nc_files = []
for pt in pass_tile_uni:
    this_pt = glob.glob(str(L2HR_folder / ('*100m_*%s*.nc' % (pt))))[0]
    all_nc_files.append(this_pt)
print('There are %s Raster 100m files' % (len(all_nc_files)))
all_nc_files.sort()

ds = xr.open_dataset(all_nc_files[0])
variables_from_l2hr = list(ds.keys())
cycles = np.unique(np.array([i.split('/')[-1].split('_')[10] for i in all_nc_files]))
dates = np.unique([all_nc_files[i].split('/')[-1].split('_')[13][:8] for i in range(len(all_nc_files))])
frames = np.unique(np.array([i.split('/')[-1].split('_')[12] for i in all_nc_files]))
passes = np.unique(np.array([i.split('/')[-1].split('_')[11] for i in all_nc_files]))

print('In %s, there are %s cycles with %s matching frames for %s' % (aoi, len(cycles), len(frames), l2hr_prefix))
print('Cycles: ', cycles)
print('Frames: ', frames)

##############################################################################
## Scan Granules: Print Pass - Scene - Tile - EPSG - Extent
extension = 10000  # Extend by 10 km buffer

pixel_votes = defaultdict(int)
grid_samples = {}
representative_files = {}
granule_metadata = []

print("\n--- Granule Spatial Reference & Extent Scan ---")
for f in all_nc_files:
    fname = Path(f).name
    parts = fname.split('_')
    
    # Parse Pass (index 11), Tile/Frame (index 12), Scene/Date (index 13 or 10)
    pass_id = parts[11] if len(parts) > 11 else "N/A"
    tile_id = parts[12] if len(parts) > 12 else "N/A"
    scene_id = parts[13] if len(parts) > 13 else (parts[10] if len(parts) > 10 else "N/A")

    with xr.open_dataset(f) as ds_curr:
        wkt = ds_curr.crs.attrs.get('crs_wkt')
        epsg = int(CRS.from_wkt(wkt).to_epsg())
        spacing = round(abs(ds_curr.x[1] - ds_curr.x[0]).item(), 2)
        num_pixels = ds_curr.x.size * ds_curr.y.size
        
        # Native extent [xmin, ymin, xmax, ymax]
        native_extent = [
            round(ds_curr.x.min().item(), 2),
            round(ds_curr.y.min().item(), 2),
            round(ds_curr.x.max().item(), 2),
            round(ds_curr.y.max().item(), 2)
        ]

        # Print Granule Summary: Pass - Scene - Tile - EPSG - Extent
        print(f"Pass: {pass_id} | Scene: {scene_id} | Tile: {tile_id} | EPSG:{epsg} | Extent: {native_extent}")

        grid_id = (epsg, spacing)
        pixel_votes[grid_id] += num_pixels
        
        if grid_id not in grid_samples:
            grid_samples[grid_id] = {'x': ds_curr.x.values, 'y': ds_curr.y.values}
            representative_files[grid_id] = f

        granule_metadata.append({
            'file': f,
            'pass': pass_id,
            'tile': tile_id,
            'scene': scene_id,
            'epsg': epsg,
            'spacing': spacing
        })

majority_grid_id = max(pixel_votes, key=pixel_votes.get)
majority_epsg, majority_spacing = majority_grid_id
print(f"\nElected Majority Grid: EPSG:{majority_epsg} with {majority_spacing}m spacing "
      f"({pixel_votes[majority_grid_id]} total pixels)")

anchor_x = grid_samples[majority_grid_id]['x'][0]
anchor_y = grid_samples[majority_grid_id]['y'][0]

##############################################################################
## Compute Unified Extent in Majority EPSG
print(f"Calculating unified global extent in EPSG:{majority_epsg}...")

all_min_x, all_max_x = np.inf, -np.inf
all_min_y, all_max_y = np.inf, -np.inf

for item in granule_metadata:
    f_path = item['file']
    f_epsg = item['epsg']
    
    with xr.open_dataset(f_path) as ds_curr:
        b_box = box(ds_curr.x.min().item(), ds_curr.y.min().item(), 
                    ds_curr.x.max().item(), ds_curr.y.max().item())
        gdf = gpd.GeoDataFrame({'geometry': [b_box]}, crs=f"EPSG:{f_epsg}")
        
        # Reproject bounds to majority EPSG if granule originates in a different UTM zone
        if f_epsg != majority_epsg:
            gdf = gdf.to_crs(f"EPSG:{majority_epsg}")
            
        bounds = gdf.total_bounds  # [minx, miny, maxx, maxy]
        
        all_min_x = min(all_min_x, bounds[0])
        all_min_y = min(all_min_y, bounds[1])
        all_max_x = max(all_max_x, bounds[2])
        all_max_y = max(all_max_y, bounds[3])

def snap(val, anchor, s, direction='min'):
    if direction == 'min':
        return anchor + np.floor((val - anchor) / s) * s
    else:
        return anchor + np.ceil((val - anchor) / s) * s

f_min_x = snap(all_min_x - extension, anchor_x, majority_spacing, 'min')
f_max_x = snap(all_max_x + extension, anchor_x, majority_spacing, 'max')
f_min_y = snap(all_min_y - extension, anchor_y, majority_spacing, 'min')
f_max_y = snap(all_max_y + extension, anchor_y, majority_spacing, 'max')

new_x = np.arange(f_min_x, f_max_x + majority_spacing, majority_spacing)
new_y = np.arange(f_max_y, f_min_y - majority_spacing, -majority_spacing)

print(f"Final Grid Size: {len(new_x)} columns x {len(new_y)} rows")

##############################################################################
## Create Master Reference Grid
master_grid = xr.Dataset(coords={"x": (["x"], new_x), "y": (["y"], new_y)})
master_grid.rio.write_crs(f"EPSG:{majority_epsg}", inplace=True)

x_mesh, y_mesh = np.meshgrid(master_grid.x.values, master_grid.y.values)
master_grid['x_gradient'] = (('y', 'x'), x_mesh.astype('f4'))

height, width = len(master_grid.y), len(master_grid.x)
rows = np.arange(height).reshape(-1, 1)
cols = np.arange(width).reshape(1, -1)
pixel_switch = (rows + cols) % 2

master_grid['pixel_test'] = (('y', 'x'), pixel_switch.astype('f4'))
master_grid['pixel_test'].attrs = {
    'long_name': 'Single Pixel Checkerboard',
    'grid_mapping': 'spatial_ref',
    'description': 'Alternates between 0 and 1 at every pixel'
}

master_grid['grid_extent'] = (('y', 'x'), np.ones((len(new_y), len(new_x)), dtype='f4'))
master_grid['grid_extent'].attrs['grid_mapping'] = 'spatial_ref'

master_grid.rio.write_transform(inplace=True)

##############################################################################
## Save Output to AOI_referencegrid.nc
output_filename = f"{aoi}_referencegrid.nc"
output_grid_path = aoi_dir / output_filename

master_grid.to_netcdf(output_grid_path)
print(f"Saved reference grid to: {output_grid_path}")
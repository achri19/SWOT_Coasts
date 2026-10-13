#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 09:33:03 2025

@author: alchrist
"""

## Extract SWOT HR L2 100m and 250m pixels around gauge locations
## Files are downloaded in previous step SWOT_download_files.py
## Set to process L2 (C version)
## Manually adjust if you want a different version
## You can process one Pass or multiple

## Import Packages
import os
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import glob
import argparse
from pyproj import Transformer

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
import rioxarray
## Import custom packages 
from swot_utils import add_WSE_HRL2
from swot_plotting import plot1variable_HR_Raster          
import rasterio              

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
##############################################################################
## Set Directory
base_dir = Path(os.path.realpath(__file__)).parent.parent
# base_dir = Path(os.path.dirname(os.getcwd()))
print(base_dir)

##############################################################################
##############################################################################
## Get bounding information from CSV. If you don't have this file or you're working on a new AOI, you can define boundaries below
try:
    bounding_LUT = pd.read_csv(base_dir / 'aoi_template.csv')
    aois = list(bounding_LUT['aoi'])
except:
    bounding_LUT = pd.DataFrame([['default', np.nan, np.nan, np.nan, np.nan, 10, 15, -2, 2, 0, 1]], columns=['aoi', 'minx', 'miny', 'maxx', 'maxy', 'minssh', 'maxssh', 'minwse', 'maxwse', 'mindac', 'maxdac'])
    aois = []
    
if args.interactive:
    aoi = input('which AOI: %s' % (aois))
else:
    aoi = args.aoi

## Get bounding areas of AOI
if aoi not in aois:
    LUT = bounding_LUT[bounding_LUT['aoi'] == 'default'].reset_index(drop=True)
    print('what bounding box (lat/lon in decimal degrees)')
    LUT.loc[0, 'minx'] = float(input('xmin: '))
    LUT.loc[0, 'miny'] = float(input('ymin: '))
    LUT.loc[0, 'maxx'] = float(input('xmax: '))
    LUT.loc[0, 'maxy'] = float(input('ymax: '))
    print(LUT)
else:
    LUT = bounding_LUT[bounding_LUT['aoi'] == aoi].reset_index(drop=True)

area = [LUT.loc[0, 'minx'], LUT.loc[0, 'miny'], LUT.loc[0, 'maxx'], LUT.loc[0, 'maxy']]

## Correct any negative E-W coordinates
if area[0] < 0: 
    area2 = [360 + area[0], area[1], 360 + area[2], area[3]]
else:
    area2 = [area[0], area[1], area[2], area[3]]

print('Crop to Area: %s' % (area2))


##############################################################################
##############################################################################
### SWOT Product Information
l2version = 'D' # 'D'
if l2version == 'C':
    l2hr_prefix = 'SWOT_L2_HR_Raster' ## Version C has shortname = 2.0, but filename = C
    L2HR_folder = base_dir / 'Data' / (l2hr_prefix + '_2.0')
else:
    l2hr_prefix = 'SWOT_L2_HR_Raster' 
    L2HR_folder = base_dir / 'Data' / (l2hr_prefix + '_' + l2version)

l2_processing = 'P*%s*' % (l2version) #Get all version C and use only the highest number
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
### Free2Mean conversion from Shailen
k2 = 0.3
hperm = -0.31460

## Search for all matching SWOT data files   
short_name = l2hr_prefix + '_' + l2version
try:
    auth = earthaccess.login() 
    results = earthaccess.search_data(short_name=short_name, 
                                      bounding_box=(area[0], area[1], area[2], area[3]),
                                      granule_name='*_100m_*%s*' % (l2_processing))
    search_list = [i['meta']['native-id'] + '.nc' for i in results]

    all_nc_files = glob.glob(str(L2HR_folder / '*.nc'))#('%s_%s_*%s*.nc' %(l2hr_prefix,resolution,l2_processing))))
    all_nc_files = [i for i in all_nc_files if i.split('/')[-1] in search_list]
except Exception as e:
    print(f"Earthaccess login or search skipped / offline mode: {e}")
    all_nc_files = glob.glob(str(L2HR_folder / '*_100m_*.nc'))

all_nc_files.sort()

with xr.open_dataset(all_nc_files[0]) as ds_init:
    variables_from_l2hr = list(ds_init.keys())

cycles = np.unique(np.array([i.split('/')[-1].split('_')[10] for i in all_nc_files]))
dates = np.unique([all_nc_files[i].split('/')[-1].split('_')[13][:8] for i in range(len(all_nc_files))])
frames = np.unique(np.array([i.split('/')[-1].split('_')[12] for i in all_nc_files]))
passes = np.unique(np.array([i.split('/')[-1].split('_')[11] for i in all_nc_files]))
print('In %s, there are %s cycles with %s matching frames for %s' % (aoi, len(cycles), len(frames), l2hr_prefix))
              
variables_to_save = variables_from_l2hr 
variables_to_add = ['egm08_meantide', 'egm08_tidefree', 'wse_egm08']
variables_to_save = variables_to_save + variables_to_add
variables_to_remove = ['crs']#,'wind_speed_rad','rad_surface_type_flag','rad_tmb_187','rad_tmb_238','rad_tmb_340','rad_water_vapor','rad_cloud_liquid_water']
variables_to_save = [x for x in variables_to_save if x not in variables_to_remove]


#############################################################################
#############################################################################
### Hardcoded Plotting Parameters
minsig0 = -10
maxsig0 = 100
minssh = -30
maxssh = 0
minwse = -2
maxwse = 2

## Determine which cycles to process (if not all)
if args.index_array is not None:
    cycles = [cycles[args.index_array]]
print('Process Cycles: %s' % (cycles))

ref_template = xr.open_dataset(aoi_dir / f"{aoi}_referencegrid.nc")
if ref_template.rio.crs is None:
    # Ensure reference has CRS (adjust if reference is not EPSG:4326)
    ref_epsg = CRS.from_wkt(ref_template.spatial_ref.attrs['crs_wkt']).to_epsg()
    ref_template.rio.write_crs(f"EPSG:{ref_epsg}", inplace=True)
else:
    ref_epsg = ref_template.rio.crs.to_epsg()

passes = np.unique(['%s' % (i) for i in ast.literal_eval(LUT['pass'][0])])
print('Passes: %s' % (passes)) 

## loop through chosen cycles and plot all frames
processed_granules = 0
for cycle in cycles[:]:
    print('\n[%s]' % (cycle)) 
    matching_cycles = fnmatch.filter(all_nc_files, '*%s_%s_*_%s_*_*_*_*_%s_*.nc' % (l2hr_prefix, resolution, cycle, l2_processing))
    matching_cycles.sort()

    disk_passes = np.unique([f.split('/')[-1].split('_')[11] for f in matching_cycles])
    active_passes = np.intersect1d(passes, disk_passes)

    for passs in active_passes:
        pass_tile = np.unique(['%s' % (j) for i, j in zip(ast.literal_eval(LUT['pass'][0]), ast.literal_eval(LUT['tile'][0])) if i == passs])

        matching_passes = fnmatch.filter(matching_cycles, '*%s_%s_*_%s_%s_*_*_*_%s_*.nc' % (l2hr_prefix, resolution, cycle, passs, l2_processing))
        datetimes = np.unique(np.array([i.split('/')[-1].split('_')[13][:9] for i in matching_passes]))
        print('[%s-%s]' % (cycle, passs)) 
        print('[%s-%s] Expected Tiles = %s' % (cycle, passs, pass_tile)) 

        print('[%s-%s] Dates: %s' % (cycle, passs, datetimes)) 
        for datefull in datetimes:
            print('[%s-%s-%s]' % (cycle, passs, datefull)) 
            
            matching_dates = fnmatch.filter(matching_passes, '*%s_%s_*_%s_%s_*_%s*_*_%s_*.nc' % (l2hr_prefix, resolution, cycle, passs, datefull, l2_processing))
            matching_dates.sort()
            date = pd.to_datetime(np.unique(np.array([i.split('/')[-1].split('_')[13][:15] for i in matching_dates]))).mean().strftime('%Y%m%d%H%M')
            year = int(date[:4])
            
            frames = np.unique(np.array([i.split('/')[-1].split('_')[12] for i in matching_dates]))            
            if np.isin(pass_tile, frames).all() == False:
                print('Not all tiles are available')

            prefixes = np.unique(['_'.join(matching_dates[i].split('_')[:-2]) for i in range(len(matching_dates))])
            print('[%s-%s-%s] There are %s matching frames: %s ' % (cycle, passs, datefull, len(matching_dates), frames))

            o_100 = '%s_%s_C%s_P%s_%s_%s_%s_%s' % (l2hr_prefix, resolution, cycle, passs, '-'.join(frames), date, l2_processing.replace('*', 'X'), aoi)
            if args.limited:
                o_100 += '_limited'
            if args.crop:
                o_100 += '_crop'

            output_filename = L2HR_output_folder / (o_100 + ".nc")
            if os.path.isfile(output_filename) == False:
                # lons,lats, wses, times, wse_qual, sig0, geoid, xover,n_wse,inc,water_frac  = [],[],[],[],[],[],[],[],[],[],[]
                mosaic = None
                for file in prefixes:
                    matching_files = fnmatch.filter(matching_dates, '*%s*' % (file.split('/')[-1]))
                    matching_files.sort()
                    file = matching_files[-1] ## choose the higher processing level (PGC0_02 over PGC0_01)
                    print('[%s-%s-%s] File: %s' % (cycle, passs, datefull, file.split('/')[-1]))                    
                    print('[%s-%s-%s] Pass: %s' % (cycle, passs, datefull, file.split('/')[-1].split('_')[11]))             
                    print('[%s-%s-%s] Frames: %s' % (cycle, passs, datefull, file.split('/')[-1].split('_')[12]))               
                    
                    ds = xr.open_dataset(file)#.drop_vars('crs', errors='ignore')
                    wkt_string = ds.crs.attrs['crs_wkt']
                    crs_obj = CRS.from_wkt(wkt_string)      
                    epsg_code = crs_obj.to_epsg()
                    
                    ## convert datetime types to float, needed for rioxarray to do reprojection
                    if 'illumination_time' in ds:
                        ds['illumination_time'] = ds['illumination_time'].astype('float64')
                    if 'illumination_time_tai' in ds:
                        ds['illumination_time_tai'] = ds['illumination_time_tai'].astype('float64')                                
                    
                    ds.rio.write_crs(f"EPSG:{epsg_code}", inplace=True)
                    ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
                    
                    # Separate variables by data type & assign NoData values before reprojection
                    float_vars = []
                    int_vars = []
                    
                    for v in ds.data_vars:
                        if 'x' not in ds[v].dims or 'y' not in ds[v].dims:
                            continue
                        
                        # Register NoData values if not present
                        if ds[v].rio.nodata is None:
                            if np.issubdtype(ds[v].dtype, np.floating):
                                ds[v] = ds[v].rio.write_nodata(np.nan)
                            else:
                                fill_val = ds[v].attrs.get('_FillValue', 255)
                                ds[v] = ds[v].rio.write_nodata(fill_val)

                        is_discrete = (
                            np.issubdtype(ds[v].dtype, np.integer) or 
                            any(k in v.lower() for k in ['qual', 'flag', 'mask'])
                        )
                        
                        if is_discrete:
                            int_vars.append(v)
                        else:
                            float_vars.append(v)
                    # Inspect grid alignment against reference template
                    src_x0 = float(ds.x.values[0])
                    src_y0 = float(ds.y.values[0])
                    ref_x0 = float(ref_template.x.values[0])
                    ref_y0 = float(ref_template.y.values[0])
                    res_x = float(abs(ref_template.x.diff('x').values[0]))
                    res_y = float(abs(ref_template.y.diff('y').values[0]))

                    # Calculate sub-pixel phase offset
                    x_offset = (src_x0 - ref_x0) % res_x
                    y_offset = (src_y0 - ref_y0) % res_y
                    
                    # Check if offset is an exact integer multiple of resolution (within 1 mm tolerance)
                    is_exact_lattice = (
                        (np.isclose(x_offset, 0.0, atol=1e-3) or np.isclose(x_offset, res_x, atol=1e-3)) and
                        (np.isclose(y_offset, 0.0, atol=1e-3) or np.isclose(y_offset, res_y, atol=1e-3))
                    )
                    crs_matches = (epsg_code == ref_epsg)

                    print(f"\t\t\t[Resample Diagnostic] CRS: EPSG:{epsg_code} -> EPSG:{ref_epsg} (Match: {crs_matches})")
                    print(f"\t\t\t[Resample Diagnostic] Phase Offset: dx = {x_offset:.2f} m, dy = {y_offset:.2f} m (Grid Size: {res_x:.1f} m)")
                    
                    if crs_matches and is_exact_lattice:
                        print("\t\t\t[Resample Diagnostic] -> Grid centers align exactly. Bilinear interpolation will yield exact, unblended values.")
                    elif crs_matches:
                        print("\t\t\t[Resample Diagnostic] -> Fractional grid offset detected. Bilinear interpolation will blend neighboring pixels.")
                    else:
                        print("\t\t\t[Resample Diagnostic] -> Non-matching CRS. Full reprojection and spatial warping active.")

                    print(f"\t\t\t[Resampling] Executing: {len(float_vars)} float vars (Bilinear) | {len(int_vars)} discrete vars (Nearest Neighbor)")

                    # Resample float variables using Bilinear interpolation
                    ds_floats = ds[float_vars].rio.reproject_match(
                        ref_template, 
                        resampling=rasterio.enums.Resampling.bilinear
                    )
                    
                    # Resample integer/flag variables using Nearest Neighbor interpolation
                    ds_ints = ds[int_vars].rio.reproject_match(
                        ref_template, 
                        resampling=rasterio.enums.Resampling.nearest
                    )
                    
                    # # Resample float variables using Bilinear interpolation
                    # ds_floats = ds[float_vars].rio.reproject_match(
                    #     ref_template, 
                    #     resampling=rasterio.enums.Resampling.bilinear
                    # )
                    
                    # # Resample integer/flag variables using Nearest Neighbor interpolation
                    # ds_ints = ds[int_vars].rio.reproject_match(
                    #     ref_template, 
                    #     resampling=rasterio.enums.Resampling.nearest
                    # )
                    
                    # Merge reprojected datasets and retain original global metadata
                    ds_aligned = xr.merge([ds_floats, ds_ints], combine_attrs='override')
                    ds_aligned.attrs = ds.attrs
                       
                    ## Exclude "bad" pixels, this is required before mosaicing files. wse_qual ==3 is used to filter pixels outside of the
                    ds_aligned = ds_aligned.where(~np.isnan(ds_aligned.sig0))
    
                    valid_pixels = ds_aligned['wse'].notnull().sum().item()
                    if valid_pixels == 0:
                        print('\t\t\tNo Valid Pixels')
                    else:
                        print(f"\t\t\tValid 'wse' pixels found: {valid_pixels}")
                        if mosaic is None:
                            mosaic = ds_aligned.copy()
                        else:
                            mosaic = mosaic.combine_first(ds_aligned)
                    
                    ds.close()
                            
                if mosaic is not None:
                    mosaic = mosaic.rename({'crs': 'spatial_ref'})
                    transformer = Transformer.from_crs("EPSG:4326", mosaic.rio.crs, always_xy=True)
                    utm_minx, utm_miny = transformer.transform(area[0], area[1])
                    utm_maxx, utm_maxy = transformer.transform(area[2], area[3])
                    mosaic_clip = mosaic.rio.clip_box(minx=utm_minx, miny=utm_miny, maxx=utm_maxx, maxy=utm_maxy)
                    del mosaic
                    
                    ## Return datetime format                
                    if mosaic_clip['illumination_time'].dtype != 'datetime64[ns]':
                        mosaic_clip['illumination_time'] = mosaic_clip['illumination_time'].where(mosaic_clip['illumination_time'] > 0).astype('datetime64[ns]')
                        mosaic_clip['illumination_time_tai'] = mosaic_clip['illumination_time_tai'].where(mosaic_clip['illumination_time_tai'] > 0).astype('datetime64[ns]')
    
                    ###########################################
                    ## Geoids
                    mosaic_clip['-egm08_meantide'] = -mosaic_clip.geoid.astype('float32')
    
                    ## Calculate del h, mean-tide geoid height - tide-free geoid height
                    ## SWOT User Handbook 11.3.1, Equation 11.4 and 11.5
                    # mosaic_clip = calculate_delh_MeantTide_minus_TideFree_geoidheights(mosaic_clip)
                    
                    ## Convert the geoid from mean-tide to tide-free
                    geoid = 'egm08_tidefree'
                    # mosaic_clip = convert_geoid_fromMeanTide_toTideFree(mosaic_clip,geoid_file) 
                    
                    mosaic_clip['free2mean'] = ((1.0 + k2) * hperm * np.sqrt(5 / (4 * np.pi)) * (-0.5 + 1.5 * np.sin(np.deg2rad(mosaic_clip['latitude'])) * np.sin(np.deg2rad(mosaic_clip['latitude'])))).astype('float32')
                    mosaic_clip['free2mean'].attrs = {'long_name': 'free_to_mean_conversion',
                                                       'units': 'm',
                                                       'comment': 'Conversion factor from mean tide geoid to tide free geoid using (1.0 + k2)*hperm*sqrt(5/(4*pi))*(-0.5 + 1.5*sin(lat)*sin(lat)) where k2 = 0.3 hperm = -0.31460'
                                                       }     
                    mosaic_clip['egm08_tidefree'] = mosaic_clip.geoid - mosaic_clip['free2mean'].astype('float32')
                    mosaic_clip['egm08_tidefree'].attrs = {'long_name': 'tidefree_EGM08_geoid_height_above_reference_ellipsoid',
                                                           'units': 'm',
                                                           'comment': 'Geoid height above the reference ellipsoid to refer the value to the tide-free system, i.e. excludes the permanent tide (zero frequency). Calculated as the original SWOT geoid (mean tide) - free2mean'
                                                           }     
        
                    ############################################
                    ## Calculate WSE as wse - geoid
                    final_corrections = [geoid, '-egm08_meantide']#,'delhpd']
                    mosaic_clip = add_WSE_HRL2(mosaic_clip, 'wse', final_corrections)    
                                    
                    time = pd.to_datetime(mosaic_clip['illumination_time'].values.flatten(), unit='ms')
                    time = time[(~pd.isnull(time)) & (pd.to_datetime(time).year == year)] 
                    
                    avgtime = pd.to_datetime(np.nanmean(time.values.astype(np.int64))).tz_localize('UTC')
                    if pd.notnull(avgtime):
                        time_label = avgtime.strftime('%Y%m%d %H%M')
                    else:
                        time_label = "00000000 0000"
                    mosaic_clip['avgtimestr'] = (('y', 'x'), np.full(mosaic_clip.wse.shape, time_label).astype('S13'))
                    mosaic_clip['avgtimestr'].attrs = {'long_name': 'granule_average_time_str',
                                                        'units': 'time string %Y%m%d %H%M',
                                                        'comment': 'Average time across the entire granule, formatted as %Y%m%d %H%M'
                                                        }                
                    ############################################
                    ## Filter pixels by cross track distance, wse_qual, and year
                    filters = xr.where(
                        (abs(mosaic_clip['cross_track']) >= 10000) & 
                        (abs(mosaic_clip['cross_track']) <= 60000) & 
                        (mosaic_clip['wse_qual'] == 0) & 
                        (mosaic_clip['illumination_time'].dt.year == int(year)), 1, np.nan)
                        
                    mosaic_clip['good'] = filters.astype('float32') 
                    mosaic_clip['good'].attrs = {'long_name': 'Good_pixels',
                                                 'comment': '1 if cross_track 10-60km, wse_qual = 0, and year matches',
                                                 }
                    
                    file_name = Path(file).name
                    mosaic_clip['file'] = (('y', 'x'), np.full(mosaic_clip.wse.shape, file_name, dtype=object))
                    mosaic_clip['file'].attrs = {'long_name': 'Name of the original file',
                                                 'units': 'string',
                                                 }      
                    
                    mosaic_clip = mosaic_clip.drop_vars('-egm08_meantide')
    
                    if args.crop:
                        try:
                            crop = gpd.read_file(aoi_dir / ('%s_crop.shp' % (aoi))).to_crs(ref_epsg)
                            mosaic_clip = mosaic_clip.rio.clip(crop, drop=True, all_touched=True)
                            print('\t\t\tcropped to %s_crop.shp' % (aoi))
                        except:
                            print('\t\t\tFull extent used')
                
                    if args.limited:
                        limited = ['x', 'y', 'longitude', 'latitude', 'wse', 'wse_qual', 'wse_uncert', 'water_frac', 'sig0', 'illumination_time', 'wse_egm08', 'geoid', 'egm08_tidefree', 'free2mean']
                        mosaic_clip = mosaic_clip.drop_vars([i for i in mosaic_clip.data_vars if i not in limited])
                    # Save the file
                    encoding = {}
                    compress = {'zlib': True, 'complevel': 5}
                    for var in mosaic_clip.variables:
                        if var in ['x', 'y']:
                            continue
                        if var == 'avgtimestr':
                            encoding[var] = {**compress, 'dtype': 'S13'}
                            continue
                        if var == 'file':
                            encoding[var] = {**compress, 'dtype': f'S{len(file_name)}'}
                            continue
                        if 'illumination' in var:
                            encoding[var] = {**compress, 'dtype': 'float64'}
                        else:
                            encoding[var] = {**compress, 'dtype': mosaic_clip[var].dtype}
                            
                    mosaic_clip.rio.write_crs(f"EPSG:{ref_epsg}", inplace=True)
                    mosaic_clip.rio.write_transform(inplace=True)
                    mosaic_clip = mosaic_clip.set_coords("spatial_ref")
                    for var in mosaic_clip.data_vars:
                        if 'grid_mapping' in mosaic_clip[var].attrs:
                            del mosaic_clip[var].attrs['grid_mapping']
                        mosaic_clip[var].attrs['grid_mapping'] = 'spatial_ref' 
    
                    has_data = mosaic_clip.wse.notnull().any().item()
                    if has_data:    
                        mosaic_clip.to_netcdf(output_filename, encoding=encoding, mode='w')
                    
                    del mosaic_clip
            else:
                print('\t\talready processed') 

            if os.path.exists(output_filename):
                print('Making plots for %s...' % (o_100))
                check1 = L2HR_fig_folder / (o_100 + '_%s_sig0.png' % (aoi))
                check3 = L2HR_fig_folder / (o_100 + '_%s_wse.png' % (aoi))
                aspect = (area[2] - area[0]) / (area[3] - area[1])

                with xr.open_dataset(output_filename) as proc_ds:
                    # Plot sig0
                    fig, ax = plt.subplots(1, 1, figsize=(10 * aspect, 10), subplot_kw=dict(projection=ccrs.PlateCarree()))
                    plot1variable_HR_Raster(proc_ds, 'L2 HR Raster 100m', minsig0, maxsig0, 'sig0', 'viridis', ax, area2)
                    fig.suptitle(datefull + ': sig0')
                    fig.tight_layout()
                    fig.savefig(check1, dpi=200)
                    plt.close(fig)

                    # Plot Water Surface Elevation
                    fig3, ax3 = plt.subplots(1, 1, figsize=(10 * aspect, 10), subplot_kw=dict(projection=ccrs.PlateCarree()))
                    plot1variable_HR_Raster(proc_ds, 'L2 HR Raster 100m', minwse, maxwse, 'wse_egm08', 'viridis', ax3, area2)
                    fig3.suptitle(datefull + ': Water Surface Elevation')
                    fig3.tight_layout()
                    fig3.savefig(check3, dpi=200)
                    plt.close(fig3)

                processed_granules += 1

print('\nProcessing all done (%s granules)' % (processed_granules))
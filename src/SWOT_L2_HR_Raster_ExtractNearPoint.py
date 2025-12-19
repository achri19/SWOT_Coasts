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

import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import ast
import sys

## Import custom packages 
from swot_utils import (haversine, 
                        add_WSE_HRL2,
                        calculate_delh_MeantTide_minus_TideFree_geoidheights, 
                        convert_geoid_fromMeanTide_toTideFree,
                        )                            

##### Get command inputs
parser=argparse.ArgumentParser(
    description='''Script to download SWOT data using earthaccess ''')

parser.add_argument('--interactive', dest='interactive',action="store_true",  help='if interactive is chosen, you will be prompted to enter information need for search')
parser.add_argument('--aoi', dest='aoi', type=str, help='If not interactive, you must provide the AOI name')
parser.add_argument('--i', dest='index_array',type = int,  help='if running job array, what index number')
parser.add_argument('--t', dest='types',type = str,  help='profiles or points')
parser.add_argument('--r', dest='search_radius_hr',type = float,  help='search radius in km')
args = parser.parse_args()
if len(sys.argv) < 2:
    args.interactive = True


##############################################################################
##############################################################################
## Set Directory
base_dir = Path(os.path.realpath(__file__)).parent.parent
# base_dir = Path(os.path.dirname(os.getcwd())).parent
print(base_dir)

##############################################################################
##############################################################################
## Get bounding information from CSV. If you don't have this file or you're working on a new AOI, you can define boundaries below
try:
    bounding_LUT = pd.read_csv(base_dir / 'point_template.csv')
    aois = list(np.unique(bounding_LUT['aoi']))
except:
    bounding_LUT = pd.DataFrame([['default',np.nan,np.nan,np.nan,np.nan,10,15,-2,2,0,1]],columns=['aoi','minx','miny','maxx','maxy','minssh','maxssh','minwse','maxwse','mindac','maxdac'])
    aois = []

    
if args.interactive:
    aoi = input('which AOI: %s'%(aois))
    types = input('profiles or points? ')
    search_radius_hr = float(input('search radius in km? '))
else:
    aoi = args.aoi
    types = args.types
    search_radius_hr = float(args.search_radius_hr)
 
# ## Get bounding areas of AOI
# if aoi not in aois:
#     LUT =  bounding_LUT[bounding_LUT['aoi']=='default'].reset_index()
#     print('what bounding box (lat/lon in decimal degrees)')
#     LUT['minx'][0] = input('xmin: ')
#     LUT['miny'][0] = input('ymin: ')
#     LUT['maxx'][0] = input('xmax: ')
#     LUT['maxy'][0] = input('ymax: ')
# else:
#     LUT =  bounding_LUT[bounding_LUT['aoi']==aoi].reset_index()

# area = [LUT['minx'][0],LUT['miny'][0],LUT['maxx'][0],LUT['maxy'][0]]

# ## Correct any negative E-W coordinates
# if area[0] <0: 
#     area2 = [360 + area[0],area[1],360+area[2],area[3]]
# else:
#     area2 = [  area[0],area[1],area[2],area[3]]
#     area = [  area[0]-360,area[1],area[2]-360,area[3]]

# print('Crop to Area: %s-%s' %(np.round(area[1],1),np.round(area[3],1)))
# print('Crop to Area: %s-%s' %(np.round(area2[1],1),np.round(area2[3],1)))

##############################################################################
##############################################################################
### SWOT Product Information
l2version = 'D' # 'C' ## Version C has shortname = 2.0, but filename = C
l2hr_prefix = 'SWOT_L2_HR_Raster' 
if l2version == 'C':
    L2HR_folder = base_dir / 'Data' / (l2hr_prefix + '_2.0') / aoi
else:
    L2HR_folder = base_dir / 'Data' / (l2hr_prefix + '_' + l2version) / aoi

l2_processing = 'P*%s*' %(l2version) #Get all version C and use only the highest number
resolution = '100m'

### Create directory for L2 data
L2HR_fig_folder = L2HR_folder / 'plots'
Path(L2HR_fig_folder).mkdir(parents=True, exist_ok=True)

L2HR_output_folder = L2HR_folder / 'final'
Path(L2HR_output_folder).mkdir(parents=True, exist_ok=True)

## Get L2 LR variables from reference dataset
ref_l2 = glob.glob(str(L2HR_folder / ('%s_%s_*%s*.nc' %(l2hr_prefix,resolution,l2_processing))))[0]
ds = xr.open_dataset(ref_l2)
variables_to_save = list(ds.keys())





















##### Get Point Information
output_dir = base_dir / 'Examples' / aoi
Path(output_dir).mkdir(parents=True, exist_ok=True)

coords = pd.read_csv(base_dir / 'point_template.csv')
coords = coords[coords['aoi']==aoi].reset_index()

g_names = coords['name']
gs = np.arange(0,len(g_names))
g_lons = coords['longitude'].astype(float)
g_lats = coords['latitude'].astype(float)


##############################################################################
## Get tide-free geoid
geoid_name = 'geoid_egm2008_wgs84_tidefree_min1x1_20200109T140029_v100_fixedshift.nc'
geoid_file = glob.glob(str(base_dir / ('*/' + geoid_name)))
### Free2Mean conversion from Shailen
# k2 = 0.3
# hperm = -0.31460
# def free2mean(lat,k2,hperm):
#     return (1.0 + k2)*hperm*np.sqrt(5/(4*np.pi))*(-0.5 + 1.5*np.sin(lat)*np.sin(lat))
# meantide = free2mean(meantide_egm08['lat'],k2,hperm)
# lat_2d = meantide.broadcast_like(meantide_egm08['geoid']) 
# search_radius_hr = 0.5# 1 #km

print(len(gs))
if args.index_array==None:
    gs = gs
else:
    gs = [gs[args.index_array]]

variables_to_add = ['egm08_meantide','egm08_tidefree','wse_egm08']
variables_to_save = variables_to_save + variables_to_add
variables_to_remove = ['crs']#,'wind_speed_rad','rad_surface_type_flag','rad_tmb_187','rad_tmb_238','rad_tmb_340','rad_water_vapor','rad_cloud_liquid_water']
variables_to_save = [x for x in variables_to_save if x not in variables_to_remove]
# swot_products = [l2hr_prefix+ '_' + i for i in ['100m','250m'] ]
swot_products = [l2hr_prefix+ '_' + i for i in ['100m'] ]

for g in gs:
    ## Gauge Name
    gauge_name = g_names[g]
    print('\n\n\n\n[%s-%s]' %(gauge_name,g))
    
    gauge_info = coords[coords['name']==gauge_name].reset_index()

    ## Gauge Coordinates
    g_lon = g_lons[g]
    g_lat = g_lats[g]
    if g_lon < 180:
        g_lon = g_lon + 360
    print('Lat: %s Lon: %s' %(g_lat,g_lon))
    
    ## Determine which tiles cover the gauge
    try:
        passes_tiles = ['%s_%s' %(i,j) for i,j in zip(ast.literal_eval(gauge_info['pass'][0]),ast.literal_eval(gauge_info['tile'][0]))] 
    except:
        print('points_template.csv does not contain pass, scene, and tile information. Go back and run SWOT_get_pass-scene-tile.py')
        break
    print(passes_tiles)
    
    for product in swot_products:
        print('\n\n\n\n[%s]' % product)
        
        
        
        

        



        output1 = output_dir / ('%s_%s_%s_%skm_%s_%s.geojson' %(aoi,product,l2_processing,search_radius_hr,gauge_name,types))
        if (os.path.isfile(output1)==False) :
            save_swot_df = {k: [] for k in variables_to_save + (['dist','unipix','avgtime','mode','cycle','pass_tile','file','avgtimestr','gauge','good'])}
    
            for pass_tile in passes_tiles:
                matching_passtiles = glob.glob(str(L2HR_folder/ ('%s_*_%s_*%s*.nc' %(product,pass_tile,l2_processing))))
                print('\n\n\n[%s %s %s]\tPass-Tile: %s ' %(aoi, product, gauge_name, pass_tile))
                pass_dates = np.unique([matching_passtiles[i].split('/')[-1].split('_')[13][:8] for i in range(len(matching_passtiles))])
                pass_dates.sort()
                print('\n\n\n[%s %s %s]\tPass-Tile: %s --> %s dates' %(aoi, product, gauge_name, pass_tile, len(pass_dates)))
                for date in pass_dates:
                    year = int(date[:4])
                    matches = sorted(glob.glob(str(L2HR_folder/ ('%s_*%s*%s*.nc' %(product,pass_tile,date)))))
                    if len(matches)>0:
                        print('\n[%s %s %s %s]\t\tDate: %s with %s matching files' %(aoi, product, gauge_name, pass_tile,date,len(matches)))
                        file = matches[-1]
                        fn = file.split('/')[-1]
                        cycle = fn.split('_')[10]
                        
                        
                        
                        ds = xr.open_dataset(file)
                        x = ds['longitude'].values.flatten()
                        y = ds['latitude'].values.flatten()
                        time = ds['illumination_time'].values.flatten()
                        
                        ############################################
                        ## Geoids
                        ds['egm08_meantide'] = ds.geoid
                        ds['-egm08_meantide'] = -ds.geoid
                        
                        ## Calculate del h, mean-tide geoid height - tide-free geoid height
                        ## SWOT User Handbook 11.3.1, Equation 11.4 and 11.5
                        ds = calculate_delh_MeantTide_minus_TideFree_geoidheights(ds)
                        
                        ## Convert the geoid from mean-tide to tide-free
                        geoid = 'egm08_tidefree'
                        ds = convert_geoid_fromMeanTide_toTideFree(ds,geoid_file) 
                        
                        ############################################
                        ## Calculate WSE as wse - geoid
                        final_corrections = [geoid,'-egm08_meantide']#,'delhpd']
                        ds = add_WSE_HRL2(ds,'wse',final_corrections)    
                        
                        ############################################
                        ## Find points within the search radius
                        msk = []
                        dists = []
                        for i in range(len(y)):
                            dist2gauge = haversine(g_lat,g_lon, y[[i]], x[[i]])
                            if dist2gauge <= search_radius_hr:
                                msk.append(i)        
                                dists.append(i)
                                
                        ############################################
                        ## Filter pixels by cross track distance, wse_qual, and year
                        filters = np.where((abs(ds.cross_track.values.flatten()[msk])>=10000) & (abs(ds.cross_track.values.flatten()[msk])<=60000) 
                                        & (ds['wse_qual'].values.flatten()[msk]==0)
                                        & (pd.to_datetime(time[msk]).year==year),True,False)
                        
                        
                                      
                                
                        
                        if len(msk)>0:                        
                            # swot_time = time[(~pd.isnull(time)) & (pd.to_datetime(time).year==year)] ## some swot files have outlier dates... remove them crudely by only including ones with correct year
                            # avg_swot_time = pd.to_datetime(np.nanmean(pd.to_datetime(swot_time,unit='ns').values.astype(np.int64)),unit='ns').tz_localize('UTC') #+ datetime.timedelta(minutes=26)
                            
                            time = pd.to_datetime(ds['illumination_time'].values.flatten()[msk],unit='ms')
                            time = time[(~pd.isnull(time)) & (pd.to_datetime(time).year==year)] 
                            avgtime = pd.to_datetime(np.nanmean(time.values.astype(np.int64))).tz_localize('UTC')
                            for var in variables_to_save:
                                save_swot_df[var].extend(ds[var].values.flatten()[msk].astype(float))
                                
                            save_swot_df['cycle'].extend([cycle for i in range(0,len(msk))])
                            save_swot_df['avgtime'].extend([avgtime for i in range(0,len(msk))])
                            save_swot_df['avgtimestr'].extend([avgtime.strftime('%Y%m%d %H%M') for i in range(0,len(msk))])
                            save_swot_df['good'].extend(filters)
    
                            save_swot_df['unipix'].extend([str(abs(a)) + str(abs(b)) for a,b in zip(x[msk],y[msk])])
                            save_swot_df['dist'].extend(dists)
                            save_swot_df['mode'].extend([product for i in range(0,len(msk))])
    
                            save_swot_df['file'].extend([file for i in msk])
                            save_swot_df['pass_tile'].extend([pass_tile for i in msk])
                            save_swot_df['gauge'].extend([gauge_name for i in msk])
                           
                            
                        else:
                            print('[%s %s P%s D%s]\t\t\tThere are no good SWOT data near the gauge' %(product, gauge_name, pass_tile,date))  
                        # print(len(save_swot_df['gauge']))
                        
                   
                                                
            if len(save_swot_df['file'])>0:
                ## Save full dataframe with SWOT variables
    
                tmp = pd.DataFrame(np.column_stack([np.column_stack(list(save_swot_df.values()))]),columns=variables_to_save+(['dist','unipix','avgtime','mode','cycle','pass_tile','file','avgtimestr','gauge','good']))
                save_swot_gdf = gpd.GeoDataFrame(tmp,geometry=gpd.points_from_xy(save_swot_df['longitude'],save_swot_df['latitude']))
                save_swot_gdf.rename(columns={'x': 'longitude', 'y': 'latitude'}, inplace=True)
                save_swot_gdf.crs = 'epsg:4326'
                
                bad_cycles = []
                save_swot_gdf = save_swot_gdf[~save_swot_gdf['cycle'].isin(bad_cycles)]
                
                try:
                    watermask = gpd.read_file(output_dir / ('%s_watermask.shp' %(aoi))).to_crs('4326').explode()
                    test = []
                    test2 = []
                    for i in range(0,len(watermask)): 
                        # save_swot_gdf = gpd.overlay(save_swot_gdf,watermask,how='intersection')
                        test.append(save_swot_gdf.within(watermask.iloc[i].geometry).values)
                    for i in range(0,len(save_swot_gdf)):
                        test3 = []
                        for j in range(0,len(watermask)):
                            test3.append(test[j][i])
                        test2.append(any(test3))
                    save_swot_gdf[str(output_dir / ('%s_water_connected_10.shp' %(aoi)))] = test2
                except:
                    print('Water mask shapefile not available so water flag is not included. If you want to manually flag water pixels, provide a water mask file in the Examples/%s subfolder' %(aoi))
                
                avgtimes = pd.to_datetime(save_swot_gdf['avgtime'])#.dt.tz_localize('UTC')
                save_swot_gdf['avgtime']=avgtimes
                save_swot_gdf.set_index('avgtimestr',inplace=True)
                if len(save_swot_gdf)>0: 
                    for var in save_swot_gdf.columns: 
                        if var not in ['avgtime','avgtimestr','gauge_time']:
                            try: 
                                save_swot_gdf[var] = pd.to_numeric(save_swot_gdf[var])
                            except: 
                                print(var)
                            
                                
                    save_swot_gdf.to_file(str(output1))
                    print('[%s] Output=' % product,output1)
                else:
                    print('[%s] There are no good SWOT pixels in any passes for this product' % product)
            print('[%s] Output=' % product,output1)        

     
            
     


       
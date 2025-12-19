#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 09:33:03 2025

@author: alchrist
"""

## Extract SWOT LR L2 and L3 2km and 250m pixels around gauge locations
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
import ast
import sys


## Import custom packages 
from swot_utils import (haversine,
                        )




##### Get command inputs
parser=argparse.ArgumentParser(
    description='''Script to download SWOT data using earthaccess ''')

parser.add_argument('--interactive', dest='interactive',action="store_true",  help='if interactive is chosen, you will be prompted to enter information need for search')
parser.add_argument('--aoi', dest='aoi', type=str, help='If not interactive, you must provide the AOI name')
parser.add_argument('--i', dest='index_array',type = int,  help='if running job array, what index number')
parser.add_argument('--t', dest='types',type = str,  help='profiles or points')
parser.add_argument('--r', dest='search_radius',type = float,  help='search radius in km')
args = parser.parse_args()
if len(sys.argv) < 2:
    args.interactive = True


##############################################################################
##############################################################################
## Set Directory
# base_dir = Path(os.path.realpath(__file__)).parent.parent.parent
base_dir = Path(os.path.dirname(os.getcwd())).parent


##############################################################################
##############################################################################
## Get bounding information from CSV. If you don't have this file or you're working on a new AOI, you can define boundaries below
try:
    bounding_LUT = pd.read_csv(base_dir / 'config'/ 'aoi_SWOT_all.csv')
    aois = list(bounding_LUT['aoi'])
except:
    bounding_LUT = pd.DataFrame([['default',np.nan,np.nan,np.nan,np.nan,10,15,-2,2,0,1]],columns=['aoi','minx','miny','maxx','maxy','minssh','maxssh','minwse','maxwse','mindac','maxdac'])
    aois = []

if args.interactive:
    aoi = input('which AOI: %s'%(aois))
    types = input('profiles or points? ')
    search_radius = float(input('search radius in km? '))
else:
    aoi = args.aoi
    types = args.types
    search_radius = float(args.search_radius)
 
## Get bounding areas of AOI
if aoi not in aois:
    LUT =  bounding_LUT[bounding_LUT['aoi']=='default'].reset_index()
    print('what bounding box (lat/lon in decimal degrees)')
    LUT['minx'][0] = input('xmin: ')
    LUT['miny'][0] = input('ymin: ')
    LUT['maxx'][0] = input('xmax: ')
    LUT['maxy'][0] = input('ymax: ')
else:
    LUT =  bounding_LUT[bounding_LUT['aoi']==aoi].reset_index()

area = [LUT['minx'][0],LUT['miny'][0],LUT['maxx'][0],LUT['maxy'][0]]

## Correct any negative E-W coordinates
if area[0] <0: 
    area2 = [360 + area[0],area[1],360+area[2],area[3]]
else:
    area2 = [  area[0],area[1],area[2],area[3]]
    area = [  area[0]-360,area[1],area[2]-360,area[3]]

print('Crop to Area: %s-%s' %(np.round(area[1],1),np.round(area[3],1)))
print('Crop to Area: %s-%s' %(np.round(area2[1],1),np.round(area2[3],1)))

##############################################################################
##############################################################################
### SWOT Product Information
l2version = 'D' # 'C' ## Version C has shortname = 2.0, but filename = C
l2lr_prefix = 'SWOT_L2_LR_SSH' 
if l2version == 'C':
    L2LR_folder = base_dir / 'Data' / (l2lr_prefix + '_2.0' )
else:
    L2LR_folder = base_dir / 'Data' / (l2lr_prefix + '_' + l2version)

l2_processing = 'P*%s*' %(l2version) #Get all version C and use only the highest number


### Create directory for L2 data
L2LR_fig_folder = L2LR_folder / 'plots'
Path(L2LR_fig_folder).mkdir(parents=True, exist_ok=True)

L2LR_output_folder = L2LR_folder / 'final'
Path(L2LR_output_folder).mkdir(parents=True, exist_ok=True)

## Get L2 LR variables from reference dataset
ref_l2 = glob.glob(str(L2LR_output_folder / ('*Expert_*_*_%s.geojson' %(aoi))))[0]
ds = gpd.read_file(ref_l2)
variables_to_savel2 = list(ds.columns)

l3version = 'v2.0.1' #'v1.0.2' ## updated to latest version 2025
l3lr_prefix = 'SWOT_L3_LR_SSH_' +l3version ## Change if you want to process a different version

### Create directory for L3 data
L3LR_folder = base_dir / 'Data' / l3lr_prefix
Path(L3LR_folder).mkdir(parents=True, exist_ok=True)

L3LR_fig_folder = L3LR_folder / 'plots'
Path(L3LR_fig_folder).mkdir(parents=True, exist_ok=True)

L3LR_output_folder = L3LR_folder / 'final'
Path(L3LR_output_folder).mkdir(parents=True, exist_ok=True)

## Get L2 LR variables from reference dataset
ref_l3 = glob.glob(str(L3LR_output_folder / ('*Expert_*_*_%s.geojson' %(aoi))))[0]
ds = gpd.read_file(ref_l3)
variables_to_savel3 = list(ds.columns)



##### Get Gauge data
gauge_folder = base_dir / "Validation" / 'gauges' / aoi #'stlawrence'
if aoi == 'louisiana':
    gauge_folder = '/projects/loac_hydro/alchrist/CRMS/'
if types == 'profiles':
    coords = pd.read_csv(base_dir / "Validation" / 'gauges' / 'profiles_SWOT_all.csv')
if types == 'points':
    coords = pd.read_csv(base_dir / "Validation" / 'gauges' / 'gauges_SWOT_all.csv')
coords = coords[coords['aoi']==aoi].reset_index()
passes = [ast.literal_eval(i) for i in coords['pass']]
scenes = [ast.literal_eval(i) for i in coords['scene']]
tiles = [ast.literal_eval(i) for i in coords['tile']]
allpasses = [x for xs in passes for x in xs]
allpasses = np.unique(allpasses)
alltiles = [x for xs in tiles for x in xs]
alltiles = np.unique(alltiles)

g_names = coords['gaugename']
gs = np.arange(0,len(g_names))
g_lons = coords['longitude'].astype(float)
g_lats = coords['latitude'].astype(float)
g_firstdates = coords['firstdate']
g_lastdates = coords['lastdate']
g_files = coords['filename']















print(len(gs))
if args.index_array==None:
    gs = gs
else:
    gs = [gs[args.index_array]]



variables_to_remove = ['crs']#,'wind_speed_rad','rad_surface_type_flag','rad_tmb_187','rad_tmb_238','rad_tmb_340','rad_water_vapor','rad_cloud_liquid_water']


swot_products = ['SWOT_L2_LR_SSH_Expert','SWOT_L2_LR_SSH_Unsmoothed','SWOT_L3_LR_SSH_Expert','SWOT_L3_LR_SSH_Unsmoothed']

for g in gs:
    ## Gauge Name
    gauge_name = g_names[g]
    print('\n\n\n\n[%s-%s]: %s' %(gauge_name,g, g_files[g]))
    
    gauge_info = coords[coords['gaugename']==gauge_name].reset_index()

    ## Gauge Coordinates
    g_lon = g_lons[g]
    g_lat = g_lats[g]
    if g_lon > 180:
        g_lon = g_lon - 360
    print('Lat: %s Lon: %s' %(g_lat,g_lon))
    g_wl = coords['wl_var'][g]
    g_date = coords['date_var'][g]
    
    ## Determine which passes cover the gauge
    gauge_passes = ast.literal_eval(gauge_info['pass'][0])
    print(gauge_passes)   
   
    for product in swot_products[:]:
        
        print('\n\n\n\n[%s]' % product)
        if 'L2' in product:
            processing = l2_processing
            variables_to_save = [x for x in variables_to_savel2 if x not in variables_to_remove]

        elif 'L3' in product:   
            processing = l3version
            variables_to_save = [x for x in variables_to_savel3 if x not in variables_to_remove]
        
        output1 = gauge_folder / ('%s_%s_%s_%skm_%s_%s.geojson' %(aoi,product,processing,search_radius,gauge_name,types))
        print('[%s] Output=' % product,output1)
        if (os.path.isfile(output1)==False) :
            save_swot_df = {k: [] for k in variables_to_save + (['dist','unipix','avgtime','mode','cycle','pass','file','avgtimestr','gauge'])}
                
            for passs in gauge_passes[:]:
                matching_passes = glob.glob(str(base_dir / 'Data' / '*' / 'final'/('%s_*_%s_*_%s_lat%s-%s_%s.geojson' %(product,passs,processing,np.round(area2[1],1),np.round(area2[3],1),aoi))))
                print('\n\n[%s %s]\tPass: %s = %s' %(product, gauge_name, passs,len(matching_passes)))           
                cycles = np.unique([matching_passes[i].split('/')[-1].split('_')[5] for i in range(len(matching_passes))])
                cycles.sort()
                for cycle in cycles:
                    
                    matches = sorted(glob.glob(str(base_dir / 'Data' / '*' / 'final'/('%s_%s_%s_*_%s_lat%s-%s_%s.geojson' %(product,cycle,passs,processing, np.round(area2[1],1),np.round(area2[3],1),aoi)))))
                    if len(matches)>0:

                        file = matches[-1]
                        fn = file.split('/')[-1]
                        date1 = file.split('/')[-1].split('_')[7]
                        date2 = file.split('/')[-1].split('_')[8]
                        year = int(date1[:4])
                        print('[%s %s P%s]\t\tCycle: %s = %s-%s' %(product, gauge_name, passs, cycle, date1,date2))  
                            
                        ds = gpd.read_file(file)
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        ############################################
                        ## Find points within the search radius
                        msk = []
                        dists = []
                        for i in range(len(ds)):
                            dist2gauge = haversine(g_lat,g_lon, ds['y'][i], ds['x'][i])
                            if dist2gauge <= search_radius:
                                msk.append(i)      
                                dists.append(i)
    
    
    
    
    
    
    
    
    
    
                        
                        if len(msk)>0: 
                            # time = pd.to_datetime(ds['time'].values.flatten()[msk],unit='ms')
                            # time = time[(~pd.isnull(time)) & (pd.to_datetime(time).year==year)] 
                            # avgtime = pd.to_datetime(np.nanmean(time.values.astype(np.int64))).tz_localize('UTC')  
                            for var in variables_to_save:
                                if var in ['time','time_tai']: 
                                    save_swot_df[var].extend(pd.to_datetime(ds[var][msk]).dt.tz_localize('UTC'))
                                elif var == 'geometry':
                                    save_swot_df[var].extend(ds[var][msk])
                                else:
                                    save_swot_df[var].extend(pd.to_numeric(ds[var][msk]))

                            save_swot_df['cycle'].extend([cycle for i in msk])
                            save_swot_df['avgtime'].extend([pd.to_datetime(date1) for i in range(0,len(msk))])
                            save_swot_df['avgtimestr'].extend([pd.to_datetime(date1).strftime('%Y%m%d %H%M') for i in range(0,len(msk))])
                            
                            save_swot_df['unipix'].extend([str(abs(a)) + str(abs(b)) for a,b in zip(ds['x'][msk],ds['y'][msk])])
                            save_swot_df['dist'].extend(dists)
                            save_swot_df['mode'].extend([product for i in msk])
                            
                            save_swot_df['file'].extend([file for i in msk])
                            save_swot_df['pass'].extend([passs for i in msk])
                            save_swot_df['gauge'].extend([gauge_name for i in msk])
                            
                        else:
                            print('[%s %s P%s]\t\t\tThere are no good SWOT data near the gauge' %(product, gauge_name, passs))  
                        print(len(save_swot_df['gauge']))
                     
               
                                                    
            if len(save_swot_df['file'])>0:
                ## Save full dataframe with SWOT variables
        
                tmp = pd.DataFrame(np.column_stack([np.column_stack(list(save_swot_df.values()))]),columns=variables_to_save+(['dist','unipix','avgtime','mode','cycle','pass','file','avgtimestr','gauge']))
                save_swot_gdf = gpd.GeoDataFrame(tmp,geometry=gpd.points_from_xy(save_swot_df['x2'],save_swot_df['y']))
                save_swot_gdf.rename(columns={'x2': 'longitude', 'y': 'latitude'}, inplace=True)
                save_swot_gdf.crs = 'epsg:4326'
              
                bad_cycles = []
                save_swot_gdf = save_swot_gdf[~save_swot_gdf['cycle'].isin(bad_cycles)]
                
                watermask = gpd.read_file(gauge_folder / ('%s_water_connected_10.shp' %(aoi))).to_crs('4326').explode()
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
                save_swot_gdf[str(gauge_folder / ('%s_water_connected_10.shp' %(aoi)))] = test2
                
                avgtimes = pd.to_datetime(save_swot_gdf['avgtime'])#.dt.tz_localize('UTC')
                save_swot_gdf['avgtime']=avgtimes
                save_swot_gdf.set_index('avgtimestr',inplace=True)
                if len(save_swot_gdf)>0: 
                    for var in save_swot_gdf.columns: 
                        if var not in ['avgtime','avgtimestr','gauge_time']:
                            try: 
                                save_swot_gdf[var] = pd.to_numeric(save_swot_gdf[var])
                            except: 
                                ''#print(var)
                            
                                
                    save_swot_gdf.to_file(str(output1))
                    print('[%s] Output=' % product,output1)
                else:
                    print('[%s] There are no good SWOT pixels in any passes for this product' % product)
                
     
                                







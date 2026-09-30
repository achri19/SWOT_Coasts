#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 08:04:53 2023

@author: alchrist
"""

## If you want to download files, this script will do that.
## Earth Data Login credentials are required
## If your AOI is not one already programmed, you'll need to define a bounding box
## Template file is stored in the repo root directory: aoi_template.csv

import os
from pathlib import Path
import argparse
import ast
import earthaccess
import pandas as pd
import sys
import paramiko
import getpass
import fnmatch
import re
import shutil
import tempfile
import zipfile
import geopandas as gpd
from shapely.geometry import box

#############################################################################
#############################################################################
#############################################################################
## Set Directory
base_dir = Path(os.path.realpath(__file__)).parent.parent


bounding_LUT = pd.read_csv(base_dir /'aoi_template.csv')
aois = list(bounding_LUT['aoi'])


def crop_lake_sp_archives(version_folder, selected_aoi, area):
    """Crop matching Lake SP shapefiles and remove processed source files."""
    aoi_rows = bounding_LUT[bounding_LUT['aoi'] == selected_aoi]
    if aoi_rows.empty:
        print('No AOI metadata found for %s; leaving Lake SP archives untouched.' % selected_aoi)
        return

    record = aoi_rows.iloc[0]
    lake_value = record.get('lake')
    lake_names = [] if pd.isna(lake_value) else [
        name.strip() for name in re.split(r'[;,]', str(lake_value)) if name.strip()
    ]
    if not lake_names:
        print('No lake names found for %s; leaving Lake SP archives untouched.' % selected_aoi)
        return
    tokens = set()
    for column in ['pass', 'scene', 'tile']:
        value = record.get(column)
        if isinstance(value, str):
            try:
                values = ast.literal_eval(value)
            except (ValueError, SyntaxError):
                values = []
            tokens.update(str(item).replace('PASS_', '') for item in values)

    if not tokens:
        print('No pass/scene/tile metadata for %s; leaving Lake SP archives untouched.' % selected_aoi)
        return

    def archive_matches(path):
        name = path.stem
        return any(re.search(r'(?<![A-Za-z0-9])%s(?![A-Za-z0-9])' % re.escape(token), name)
                   for token in tokens)

    archives = [path for path in version_folder.rglob('*.zip') if archive_matches(path)]
    if not archives:
        print('No matching Lake SP ZIP archives found for %s.' % selected_aoi)
        return

    crop_box = gpd.GeoDataFrame(
        geometry=[box(area[0], area[1], area[2], area[3])],
        crs='EPSG:4326',
    )
    sidecar_suffixes = ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qix', '.fix']

    for archive in archives:
        extracted = False
        with tempfile.TemporaryDirectory(prefix='lake_sp_') as temp_dir:
            with zipfile.ZipFile(archive) as zip_file:
                zip_file.extractall(temp_dir)
            shapefiles = list(Path(temp_dir).rglob('*.shp'))
            if not shapefiles:
                print('No shapefile in %s; leaving archive untouched.' % archive.name)
                continue

            for source_shp in shapefiles:
                relative_stem = source_shp.relative_to(temp_dir).with_suffix('')
                output_stem = version_folder / relative_stem
                output_stem.parent.mkdir(parents=True, exist_ok=True)
                source_gdf = gpd.read_file(source_shp)
                lake_columns = {
                    column.lower(): column for column in source_gdf.columns
                }
                lake_column = lake_columns.get('lake_name')
                if lake_column is None:
                    print('No lake_name field in %s; leaving archive untouched.' % archive.name)
                    extracted = False
                    break
                lake_pattern = '|'.join(re.escape(name) for name in lake_names)
                source_gdf = source_gdf[
                    source_gdf[lake_column].astype('string').str.contains(
                        lake_pattern, case=False, na=False, regex=True
                    )
                ]
                if source_gdf.crs is None:
                    source_gdf = source_gdf.set_crs('EPSG:4326')
                crop_geometry = crop_box.to_crs(source_gdf.crs)
                cropped = gpd.clip(source_gdf, crop_geometry)

                with tempfile.TemporaryDirectory(prefix='lake_sp_output_') as output_dir:
                    output_shp = Path(output_dir) / source_shp.name
                    schema = gpd.io.file.infer_schema(cropped)
                    for field, field_type in schema['properties'].items():
                        if field_type.startswith('float'):
                            schema['properties'][field] = 'float:24.6'
                        elif field_type.startswith('int'):
                            schema['properties'][field] = 'int:18'
                    cropped.to_file(
                        output_shp,
                        driver='ESRI Shapefile',
                        engine='fiona',
                        schema=schema,
                    )
                    for suffix in sidecar_suffixes:
                        original = output_stem.with_suffix(suffix)
                        if original.exists():
                            original.unlink()
                    for generated in Path(output_dir).glob('%s.*' % output_shp.stem):
                        shutil.move(str(generated), str(output_stem.with_suffix(generated.suffix)))
                extracted = True
                print('Cropped %s to %s features.' % (source_shp.name, len(cropped)))

        if extracted:
            archive.unlink()
            print('Deleted processed archive %s.' % archive.name)


#############################################################################
#############################################################################
#############################################################################
## Get search parameters 
parser=argparse.ArgumentParser(
    description='''Script to download SWOT data using earthaccess ''')

parser.add_argument('--interactive', dest='interactive',action='store_true',  help='if interactive is chosen, you will be prompted to enter information need for search')
parser.add_argument('--short_name', dest='short_name', type=str, help='If not interactive, you must provide the searchable short name formatted as SWOT_AA_BB_CC_DD where AA is Level, BB is Mode, CC is Product, and DD is version. See PODAAC search keywords for details. ')
parser.add_argument('--processing', dest='processing', type=str, help='If not interactive, you must provide the processing level such as PGC0_01')
parser.add_argument('--startdate', dest='startdate', type=str, help='If not interactive, you must provide the start date as %Y-%m-%d %H:%M:%S')
parser.add_argument('--enddate', dest='enddate', type=str, help='If not interactive, you must provide the end date as %Y-%m-%d %H:%M:%S')
parser.add_argument('--aoi', dest='aoi', type=str, help='If not interactive, you must provide the AOI name')
parser.add_argument('--area',dest='area',nargs='+',type=float,help='If your AOI is new, you must define the bounding box in lat/lon dec degree coordinates (xmin,ymin,xmax,ymax)')
args = parser.parse_args()
if len(sys.argv) < 2:
    args.interactive = True

if args.interactive:
    aoi = input('which AOI (if none, you need bounding box): %s '%(aois))
    LUT =  bounding_LUT[bounding_LUT['aoi']==aoi].reset_index()

    version = str(input('Which version C or D? '))
    processing = 'P*%s*' %(version) 

    if version == 'C':
        version = '2.0'
    products = {0:'SWOT_L2_HR_PIXC',1:'SWOT_L2_HR_Raster',2:'SWOT_L2_HR_RiverSP',3:'SWOT_L2_HR_LakeSP',5:'SWOT_L2_LR_SSH'}
    short_name = products[int(input(f"Which product:\n{products}"))] + '_' + version
    L3 = 'n'
    if 'SWOT_L2_LR_SSH' in short_name:
        mode = 'LR'
        print('For LR mode, only Expert and Unsmoothed products will be downloaded')
        # product = input('Which product? SSH_EXPERT SSH_BASIC SSH_UNSMOOTHED SSH_WINDWAVE for all, use SSH ')
        L3 = input('Do you want to download L3 v2.0.1 products from AVISO (username/password required)? y/n')
    else:
        mode = 'HR'
        if 'Raster' in short_name:
            processing = '100m*' + processing
    ## Get bounding areas and search dates
    if aoi not in aois:
        print('what bounding box (lat/lon in decimal degrees)')
        area = [input('xmin: '),input('ymin: '),input('xmax: '),input('ymax: ')]
        startdate = str(input('start date: %Y-%m-%d %H:%M:%S '))
        enddate = str(input('end date: %Y-%m-%d %H:%M:%S '))
    else:
        area = [LUT['minx'][0],LUT['miny'][0],LUT['maxx'][0],LUT['maxy'][0]]
        startdate = str(LUT['startdate'][0])
        enddate = str(LUT['enddate'][0])
        if (startdate=='nan') | (enddate=='nan') :
            startdate = str(input('start date: %Y-%m-%d %H:%M:%S '))
            enddate = str(input('end date: %Y-%m-%d %H:%M:%S '))
        print('Search from %s to %s' %(startdate,enddate))
    

    
else:
    aoi = args.aoi
    short_name = args.short_name
    processing = args.processing
    startdate = args.startdate
    enddate = args.enddate
    area = args.area
 
#############################################################################
#############################################################################
#############################################################################
### Create directory for data
product_folder = base_dir /'Data' / short_name
version_folder = product_folder 

print('Search bounding box: %s' %area)
Path(version_folder).mkdir(parents=True, exist_ok=True)
print('files will be saved here: ', version_folder)
## Search NASA Earth Data for matching data products
print('short_name: %s' %(short_name))
print('Search %s - %s' %(startdate,enddate))
auth = earthaccess.login() 
results = earthaccess.search_data(short_name = short_name, 
                                  temporal = (startdate, enddate), # can also specify by time
                                  bounding_box = (area[0],area[1],area[2],area[3]),
                                  granule_name= '*%s*' %(processing))
if mode == 'LR':
    results = [i for i in results if i['umm']['GranuleUR'].split('_')[4] in ['Expert','Unsmoothed']]
if 'SWOT_L2_HR_LakeSP' in short_name:
    results = [
        i for i in results
        if '_Obs_' in i['umm']['GranuleUR'].split('/')[-1]
    ]
    print('Lake SP Obs granules selected: ', len(results))
    
if len(results)>0:
    print('Number of Matching Results = ', len(results))
    earthaccess.download(results[:], version_folder,show_progress=True)
    if 'SWOT_L2_HR_LakeSP' in short_name:
        crop_lake_sp_archives(version_folder, aoi, area)
    
    
#############################################################################
#############################################################################
#############################################################################
## Code for downloading L3 products using FTP from AVISO. 
## This requires an AVISO account and password
if L3 != 'n' :
    
    l3version = input('Which L3 version? 2.0.1 is the latest ')
    l3_folder = base_dir /'Data' / ('SWOT_L3_LR_SSH_v%s' %(l3version))
    SFTP_HOST = "ftp-access.aviso.altimetry.fr"  # Example FTP server
    username = input(f"Enter username for {SFTP_HOST}: ")
    password = getpass.getpass(f"Enter password for {username}: ")
    
    print(f"\n[INFO] Connecting to {SFTP_HOST}...")
    SSH_Client= paramiko.SSHClient()
    SSH_Client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    SSH_Client.connect( hostname=SFTP_HOST, username=username, port = 2221,
                   password= password, look_for_keys= False
                 )
    sftp_client= SSH_Client.open_sftp() 
    
    for file in results[:]:
       print(file['umm']['GranuleUR'])
       l3_cycle = file['umm']['GranuleUR'].split('_')[5]
       l3_pass = file['umm']['GranuleUR'].split('_')[6]
       AVISO_path1 = '/swot_products/l3_karin_nadir/l3_lr_ssh/v%s/Expert/cycle_%s/' %( l3version.replace('.','_'),l3_cycle)
       AVISO_path2 = '/swot_products/l3_karin_nadir/l3_lr_ssh/v%s/Unsmoothed/cycle_%s/' %( l3version.replace('.','_'),l3_cycle)
       file_list = sftp_client.listdir(AVISO_path1)
       REMOTE_FILE_PATH = [f for f in file_list if fnmatch.fnmatch(f, ('*_%s_202*' %(l3_pass)))]
       
       if (len(REMOTE_FILE_PATH)>0):
           local_file = str(l3_folder / REMOTE_FILE_PATH[0].split('/')[-1])
           if(os.path.isfile(local_file)==False):
               print(f"[INFO] Attempting to download '{REMOTE_FILE_PATH}'...")
               sftp_client.get(AVISO_path1  + REMOTE_FILE_PATH[0], local_file)
       
       file_list = sftp_client.listdir(AVISO_path2)
       REMOTE_FILE_PATH = [f for f in file_list if fnmatch.fnmatch(f, ('*_%s_202*' %(l3_pass)))]

       if (len(REMOTE_FILE_PATH)>0):
           local_file = str(l3_folder / REMOTE_FILE_PATH[0].split('/')[-1])
           if (os.path.isfile(local_file)==False):
               print(f"[INFO] Attempting to download '{REMOTE_FILE_PATH}'...")
               sftp_client.get(AVISO_path2  + REMOTE_FILE_PATH[0], local_file)

   
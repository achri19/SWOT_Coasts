#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 08:04:53 2023

@author: alchrist
"""

## If you just want to download files, this script will do that.
## Earth Data Login credentials are required
## If your AOI is not one already programmed, you'll need to define a bounding box

import os
from pathlib import Path
import argparse
import earthaccess
import pandas as pd
import sys
import paramiko
import getpass
import fnmatch
##############################################################################
##############################################################################
## Set Directory
base_dir = Path(os.path.dirname(os.getcwd()))


bounding_LUT = pd.read_csv(base_dir / 'aoi_template.csv')
aois = list(bounding_LUT['aoi'])


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
    products = {0:'SWOT_L2_HR_PIXC',1:'SWOT_L2_HR_Raster',2:'SWOT_L2_HR_RiverSP',3:'SWOT_L2_HR_LakeSP',4:'SWOT_L2_HR_PIXCVec',5:'SWOT_L2_LR_SSH'}
    short_name = products[int(input(f"Which product:\n{products}"))] + '_' + version
    L3 = 'n'
    if 'SWOT_L2_LR_SSH' in short_name:
        mode = 'LR'
        # product = input('Which product? SSH_EXPERT SSH_BASIC SSH_UNSMOOTHED SSH_WINDWAVE for all, use SSH ')
        L3 = input('Do you want to download L3 v2.0.1 products from AVISO (username/password required)? y/n')
    else:
        mode = 'HR'
    if aoi not in aois:
        print('what bounding box (lat/lon in decimal degrees)')
        area = [input('xmin: '),input('ymin: '),input('xmax: '),input('ymax: ')]
        startdate = str(input('start date: %Y-%m-%d %H:%M:%S '))
        enddate = str(input('end date: %Y-%m-%d %H:%M:%S '))
    else:
        area = [LUT['minx'][0],LUT['miny'][0],LUT['maxx'][0],LUT['maxy'][0]]
        startdate = str(LUT['startdate'][0])
        enddate = str(LUT['enddate'][0])
    

    
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
## Get bounding areas of default AOIs


### Create directory for data
product_folder = base_dir /'Data' / short_name
# Path(product_folder).mkdir(parents=True, exist_ok=True)
if mode == 'LR': 
    version_folder = product_folder 
elif mode == 'HR':
    version_folder = product_folder / aoi
Path(version_folder).mkdir(parents=True, exist_ok=True)
print('files will be saved here: ', version_folder)
print('short_name: %s' %(short_name))
## Search NASA Earth Data for matching data products
auth = earthaccess.login() 
results = earthaccess.search_data(short_name = short_name, 
                                  temporal = (startdate, enddate), # can also specify by time
                                  bounding_box = (area[0],area[1],area[2],area[3]),
                                  granule_name= '*%s*' %(processing))
if mode == 'LR':
    results = [i for i in results if i['umm']['GranuleUR'].split('_')[4] in ['Expert','Unsmoothed']]
    
if len(results)>0:
    print('Number of Matching Results = ', len(results))
    earthaccess.download(results[:], version_folder)
    
    

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

   
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 09:06:15 2024

@author: alchrist
"""

## Process all available SWOT LR L2 (Expert and Unsmoothed) and optional L3 (Expert and Unsmoothed)
## Files were downloaded in using SWOT_download_files.py
## Set to process L2 (C version) and L3 (v2.0.1)
## Manually adjust if you want a different version
## You can process one Pass or multiple
## Outputs are 3 figures (wse, sig0, ssh)
## and geojsons of SWOT Pixels, one file for each track/frame contains multiple acquisitions

## Import Packages
import numpy as np
from pathlib import Path
import os
import sys
import glob
import xarray as xr
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import argparse
import ast
import shutil

## Import custom packages (modified the full utility toolbook to include functions specifically needed for LR)
from swot_utils import (interpolate_2km_to_250m, 
                            add_WSE_LRL2,add_wse_LRL3,
                            fix_latlon, savedf,
                            calculate_delh_MeantTide_minus_TideFree_geoidheights, 
                            convert_geoid_fromMeanTide_toTideFree,
                            other_filters_exclude,other_filters_range,
                            add_tide_corrections_toSSH, add_crossover_correction_toSSH,
                            add_permanent_tide,
                            crop
                            )
from swot_plotting import (plot1variable_LR_L2orL3_exp, plot1variable_LR_L2orL3_uns)

##############################################################################
##############################################################################
## Code below is for running command line only. Skipped if running in IDE
parser=argparse.ArgumentParser(
    description='''Script to process SWOT LR data. Code will save geojson of cropped LR pixels and simple plot''')
parser.add_argument('--interactive', dest='interactive',action="store_true",  help='if interactive is chosen, you will be prompted to enter information need for search')
parser.add_argument('--aoi', dest='aoi', type=str, help='If not interactive, you must provide the AOI name')
parser.add_argument('--pass', dest='passs', type=str, help='If not interactive, you must provide the pass number. Use --pass any for wildcard')
parser.add_argument('--i', dest='index_array',type = int,  help='if running job array, what index number')
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
    bounding_LUT = pd.read_csv(base_dir / 'aoi_template.csv')
    aois = list(bounding_LUT['aoi'])
except:
    bounding_LUT = pd.DataFrame([['default',np.nan,np.nan,np.nan,np.nan,10,15,-2,2,0,1]],columns=['aoi','minx','miny','maxx','maxy','minssh','maxssh','minwse','maxwse','mindac','maxdac'])
    aois = []


### Choose whether to process LR L3 products 
if args.interactive:
    aoi = input('which AOI: %s'%(aois)).lower()
    passs = str(input('which pass? or any = * ')) #'PIB0_01'
else:
    if (args.passs is None) | (args.aoi is None) :
        parser.error("The --pass and --aoi arguments cannot be None.")
    else:
        aoi = args.aoi
        passs = args.passs
        if passs == 'any':
            passs = '*'
print(aoi)

## Get bounding box of AOI
### Coordinates should be in Lat/Lon degrees
if aoi not in aois:
    LUT =  bounding_LUT[bounding_LUT['aoi']=='default'].reset_index()
    print('what bounding box (lat/lon in decimal degrees)')
    area = [float(input('xmin: ')),float(input('ymin:')),float(input('maxx:')),float(input('maxy:'))]
else:
    LUT =  bounding_LUT[bounding_LUT['aoi']==aoi].reset_index()
    area = [LUT['minx'][0],LUT['miny'][0],LUT['maxx'][0],LUT['maxy'][0]]
    allpasses = np.unique([ast.literal_eval(i) for i in LUT['pass']])
    print('All Passes: %s' %(allpasses))
        

## Correct any negative E-W coordinates
if area[0] <0: 
    area2 = [360 + area[0],area[1],360+area[2],area[3]]
else:
    area2 = [  area[0],area[1],area[2],area[3]]
    
print('Crop to Area: %s' %(area2))


##############################################################################
##############################################################################
### SWOT Product Information
## Level 2 Products
l2version = 'D' # 'D'
if l2version == 'C':
    l2_name = 'SWOT_L2_LR_SSH_2.0' ## Version C has shortname = 2.0, but filename = C
else:
    l2_name = 'SWOT_L2_LR_SSH_%s' %(l2version) 

l2_processing = 'P*%s*' %(l2version) #Get all version C and use only the highest number

## Level 3 Products
l3version = 'v2.0.1' #'v1.0.2' ## updated to latest version 2025
l3_name = 'SWOT_L3_LR_SSH_' +l3version ## Change if you want to process a different version


### Data directories
## Level 2
L2_folder = base_dir / 'Data' / l2_name
Path(L2_folder).mkdir(parents=True, exist_ok=True)

L2_fig_folder = L2_folder / 'plots'
Path(L2_fig_folder).mkdir(parents=True, exist_ok=True)

L2_output_folder = L2_folder / 'Extracted' / aoi
Path(L2_output_folder).mkdir(parents=True, exist_ok=True)

tmp_folder = L2_folder / 'tmp'
Path(tmp_folder ).mkdir(parents=True, exist_ok=True)

## Level 3
L3_folder = base_dir / 'Data' / l3_name
Path(L3_folder).mkdir(parents=True, exist_ok=True)

L3_output_folder = L3_folder / 'Extracted'
Path(L3_output_folder).mkdir(parents=True, exist_ok=True)

L3_fig_folder = L3_folder / 'plots'
Path(L3_fig_folder).mkdir(parents=True, exist_ok=True)

reference_list = glob.glob(str(L2_folder / ('SWOT_L2_*_Unsmoothed_*_%s_*.nc' %(passs))))
reference_list.sort()
reference_list = [i for i in reference_list if i.split('/')[-1].split('_')[6] in allpasses]
if len(reference_list)==0:
    print('No matching L2 LR SSH products, make sure you downloaded files with SWOT_download_files.py first')
    sys.exit(1)
else:
    print('There are %s matching L2 Unsmoothed cycles for pass %s' %(len(reference_list),passs))


##############################################################################
##############################################################################
## For batch processing on Gattaca2 only
if len(sys.argv) > 0:
    if args.index_array==None:
        reference_list = reference_list
    else:
        reference_list = [reference_list[args.index_array]]

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

#############################################################################
#############################################################################
## Filters Options
## Make adjustments here to increase/decrease filtering
filters_crosstrackl2 =  [['cross_track_distance',-10000.,10000.],
                            ['cross_track_distance',60000.,100000.],
                            ['cross_track_distance',-100000.,-60000.],
                            ]
filters_crosstrackl3 =  [['cross_track_distance',-10.,10.],
                            ['cross_track_distance',60.,100.],
                            ['cross_track_distance',-10.,-60.],
                            ]
filters_class = [['ancillary_surface_classification_flag',0.,0.]]

## L2 Quality Flag options:
    # 'height_cor_xover_qual'
    # 'orbit_qual'
    # 'sig0_karin_2_qual'
    # 'sig0_karin_qual'
    # 'ssh_karin_2_qual'
    # 'ssh_karin_qual'
    # 'ssha_karin_2_qual'
    # 'ssha_karin_qual'
    # 'swh_karin_qual'
    # 'wind_speed_karin_2_qual'
    # 'wind_speed_karin_qual'
    # qual_flag = 'ssh_karin_2_qual'        
l2_qual_flag = 'ssh_karin_2_qual'

## L3 Quality Flag options:
     ## Flag #102: No SSHa values available
     ## Flag #101: Pixels over land.
     ## Flag #100: Edges of swath. Only values between 10 to 60 km to the nadir are considered as valid data.
     ## Flag #70: Pixels impacted by spacecraft events.
     ## Flag #50: Abnormally high SSHA values.
     ## Flag #30: SSHA pixels out of the expected statistical distribution.
     ## Flag #20: Suspected sea-ice pixels.
     ## Flag #10: Suspected coastal pixels.
     ## Flag #5: SSHA pixels out of the local distribution.
     ## Flag #0: Valid data.
l3_qual_flag = 'quality_flag'

l3_exclusions = [[l3_qual_flag,[100,101,102]]]

#############################################################################
#############################################################################
### Hardcoded Plotting Parameters
minsig0 = -10
maxsig0 = 100
minssh = -30
maxssh = 0
minwse = -2
maxwse = 2


for i in range(len(reference_list)):
    ## L2 Unsmoothed File 
    l2_uns_file = reference_list[i]

    ## Get Pass and Date Information
    fn = l2_uns_file.split('/')[-1].split('.nc')[0]
    cycle = fn.split('_')[5]
    passs = fn.split('_')[6]
    date1 = fn.split('_')[7][:8]
    date2 = fn.split('_')[8][:8]
    year = int(date1[:4])
    print('\n\n\n\nPass: %s Date: %s-%s Cycle: %s' %(passs,date1,date2,cycle))
    
    ## Prefix and Suffix improve searchign for files
    prefix = fn.split('Unsmoothed')[0]
    suffix = fn.split('Unsmoothed')[1]
    print('### Prefix: %s Suffix:  %s' %(prefix,suffix))
    
    check1 = L2_fig_folder/('%s_%s_%s_%s-%s_sig0_karin.png' %(prefix,aoi,passs,date1,date2))
    check2 = L2_fig_folder/('%s_%s_%s_%s-%s_ssh_karin.png' %(prefix,aoi,passs,date1,date2))
    check3 = L2_fig_folder/('%s_%s_%s_%s-%s_wse.png' %(prefix,aoi,passs,date1,date2))
    if (os.path.isfile(check1)==False)|(os.path.isfile(check2)==False)|(os.path.isfile(check3)==False):
    
        l3_exp_list = glob.glob(str(L3_folder / ('SWOT_L3_*_Expert_%s_%s_*.nc' %(cycle,passs))))
        l3_exp_list.sort()
        l3_uns_list = glob.glob(str(L3_folder / ('SWOT_L3_*_Unsmoothed_%s_%s_*.nc' %(cycle,passs))))
        l3_uns_list.sort()
        plots = 2
        if len(l3_exp_list)>0: 
            plots = plots + 1
        if len(l3_uns_list)>0: 
            plots = plots + 1
            
        ## Plot normalized radar cross-section
        fig,ax = plt.subplots(1,2,figsize=(15, 5),subplot_kw=dict(projection=ccrs.PlateCarree()))
        ## Plot of Sea Surface Height
        fig2,ax2 = plt.subplots(1,2,figsize=(15, 5),subplot_kw=dict(projection=ccrs.PlateCarree()))
        ## Plot corrected Water Surface Elevation
        fig3,ax3 = plt.subplots(1,plots,figsize=(15, 5),subplot_kw=dict(projection=ccrs.PlateCarree()))
    
        ## Create name for saving files
        l2_exp    = '%sExpert%s_lat%s-%s_%s' %(prefix,suffix,np.round(area[1],1),np.round(area[3],1),aoi)
        l2_uns    = '%sUnsmoothed%s_lat%s-%s_%s' %(prefix,suffix,np.round(area[1],1),np.round(area[3],1),aoi)    
    
        
        ######################################################################                                                                                                                                                                                                                                                                                                                                                                                                             
        ######################################################################                                                                                                                                                                                                                                                                                                                                                                                                             
        ######### Process L2 Expert (2km)
        ######################################################################                                                                                                                                                                                                                                                                                                                                                                                                             
        ######################################################################  
        print('[LR L2 SSH EXPERT]')
        ## Search for the matching L2 Expert file, choose the last one in the list, which should be the highest re-processing number _01 or _02
        l2_exp_file = glob.glob(str(L2_folder / ('%sExpert*%s*%s*%s*%s*.nc' %(prefix,passs,date1,date2,l2_processing))))[-1]
        print('[LR L2 SSH EXPERT] File: ', l2_exp_file)
        l2_processing_full = '_'.join(l2_exp_file.split('.nc')[0].split('_')[-2:])
        ############################################
        ## Open L2 Expert File
        L2_expert = xr.open_dataset(l2_exp_file)
        ## Variables to Save
        # variables_to_save = [v for v in list(L2_expert.keys()) if L2_expert[v].ndim==2]
        variables_to_save = [v for v in list(L2_expert.keys()) if all(dim in L2_expert[v].dims for dim in ['num_lines','num_pixels'])]
        variables_to_add = ['egm08_meantide','egm08_tidefree','ssh_xover','ssh_tides','wse']
        variables_to_save = variables_to_save + variables_to_add
        variables_to_remove = [''] #['geoid','wind_speed_rad','rad_surface_type_flag','rad_tmb_187','rad_tmb_238','rad_tmb_340','rad_water_vapor','rad_cloud_liquid_water']
        variables_to_save = [x for x in variables_to_save if x not in variables_to_remove]
        
        ## Get Time information
        time = L2_expert.time.values.flatten()
        swot_time = time[(~pd.isnull(time)) & (pd.to_datetime(time).year==year)] ## some swot files have outlier dates... remove them crudely by only including ones with correct year
        avg_swot_time = pd.to_datetime(np.nanmean(pd.to_datetime(swot_time,unit='ns').values.astype(np.int64)),unit='ns').tz_localize('UTC') #+ datetime.timedelta(minutes=26)                                                                                                                                                                                                                                                                                                                                                                                                        
        
        ############################################
        ## Geoids
        L2_expert['egm08_meantide'] = L2_expert.geoid
        
        ## Calculate del h, mean-tide geoid height - tide-free geoid height
        ## SWOT User Handbook 11.3.1, Equation 11.4 and 11.5
        L2_expert = calculate_delh_MeantTide_minus_TideFree_geoidheights(L2_expert)
        
        ## Convert the geoid from mean-tide to tide-free
        geoid = 'egm08_tidefree'
        L2_expert = convert_geoid_fromMeanTide_toTideFree(L2_expert,geoid_file) 
        ############################################    
        
        
        ############################################
        ## Filter pixels
        ## Crop to new area
        L2_expert = crop(L2_expert,area2)
        
        
        try:
            L2_expert_filtered = other_filters_range(L2_expert,filters_crosstrackl2,False,False)
            L2_expert_filtered = other_filters_range(L2_expert_filtered,filters_class,True,False)
            # L2_expert_filtered =  get_qual_flag_mask(L2_expert_filtered,qual_flag,['bad']) 
            
            ############################################
            ## On-board processing already includes some instrument corrections when calculating SSH, including media delays (Handbook Section 11.1)
            ### Wet tropospheric range delays - either from radiometer or ECMWF model (deviation = 1cm). corrections range = 0-0.45m
            ### Dry tropospheric range delays - from ECMWF model sea surface pressure fields every 6 hours (error 1cm). corrections range = -2.1 to -2.4m
            ### Ionospheric range delays - corrections for Ka bands range = 0 to -0.03m. Ku and C- band NAlt gives direct measurement of TEC, and iono corrections. Also used to compute Global Ionosphere Maps as backup
            ### Sea State Bias - difference in signal reflectivity between peaks and troughs of ocean waves, from empirically derived models corrections range = 0 to -0.4m
            
            ## Check radiometer measurements  
            if np.count_nonzero(~np.isnan(L2_expert_filtered['ssh_karin'].data))==0:
                # This scene has no good radiometer data, using ssh_karin_2 instead of ssh_karin
                # Model based corrections for wet troposhperic range delays
                # Sea state bias from sigma0 atmospheric attenuation, wind speed, and SWH
                use_ssh = '_2'
                print('\n[LR L2 SSH EXPERT] Radiometer quality is bad, use ssh_karin_2', l2_exp_file)
            else:
                use_ssh = ''
                print(']n[LR L2 SSH EXPERT] Radiometer quality is good, use ssh_karin', l2_exp_file)

            
        
            # ### There are discrepancies in the calculation for SSHA for Version C and D.         
            # ### The 2025 PDD https://deotb6e7tfubr.cloudfront.net/s3-edaf5da92e0ce48fb61175c28b67e95d/podaac-ops-cumulus-docs.s3.us-west-2.amazonaws.com/web-misc/swot_mission_docs/pdd/D-56407_SWOT_Product_Description_L2_LR_SSH_20250224a_RevC_clean_sig.pdf?A-userid=None&Expires=1754431195&Signature=ZFxZ5qTr1lySfoUyjxA7lASIlG9dyhkcCuI3OayJNlYhoDRjNYX4eIkq238HcAKU4R0TZ6Aoi1k4YAJO3yRC1MXJjDjwz-0xtuNlI84oEVArwJ-EYBNpX3aLiOBSvrjGZDa09mJq87FiKl32ohcnH73rBkx58IVJQ~PKwZyLGmpzmsQVmaGVz~VaZ8nq7G0QGsITYxw6jj2HRdO-XquWiqEJfaSC~A4a4r8UINrp~am8L7s6onBFKhkjpiISmTINA88ADIqfiak8HW-B9IrxSbMsOjLgOZVg~Zulore2EoLziy9Bt8jvcNMwb01aVt5PQUPumab0hVOFFLFTpAxzSQ__&Key-Pair-Id=K2MYOFHSV8YFSM
            # ### Using their calculation, the SSHA in Version C is not correct when calculated from SSH. There is a +/- 0.003 difference
            # ### get SSH from SSHA by removing corrections for DAC, MSS, Ocean tide (fes, non eq, and eq), Internal tide (hret)
            # if version=='C':
            #     L2_expert_filtered = add_SSH_fromSSHA_VC(L2_expert_filtered,use_ssh)
            # elif version =='D':
            #     L2_expert_filtered = add_SSH_fromSSHA_VD(L2_expert_filtered,use_ssh)
                
            
            ############################################
            ## Apply crossover correction
            L2_expert_corrected,updated_ssh_var = add_crossover_correction_toSSH(L2_expert_filtered,use_ssh)

            ############################################
            ## Tide Corrections
            ### Options:
            ### pole_tide: model for sea surface height displacement from teh geocentric pole tide, sum of solid-earth pole tide height, modeled ocean pole tide and load pole tide heights. corrects for linear drift. Over land, value only includes solod earth pole tide and load pole tide
            ### solid_earth_tide: model for the solid earth (body) tide height from Cartwright/Taylor/Edden tide-generating potential coefficients and consists of 2nd and 3rd degree constituents. The permanent tide (zero frequency) is not included
            ### ocean_tide_fes: model for sea surface height displacement from the ocean tide from FES2022b, includes sum of short-period ocean tide, short-period load tide, and long-period equilibrium ocean tide
            ### ocean_tide_got: from GOT4.10c model, includes short-period ocean tie, short-period load tide, and long-period equilibrium ocean tide
            ### load_tide_fes: model for geocentric surface height displacement from load tide, from FES2022b, includes short-period and long-period components. The short-period component is already included in the ocean_tide_fes
            ### load_tide_got: from GOT4.10c, but only includes the short-period component, which is already included in ocean_tide_got
            ### ocean_tide_eq: model for sea surface height displacement from equilibrium long-period ocean tides, already included in the ocean_tide_fes and ocean_tide_got. Over land, this = 0
            ### ocean_tide_non_eq: model for sea surface height from non-equilibrium long-period ocean tides from FES2022b, includes the long-period ocean tide with respect to ocean_tide_eq and long-period load tide. Over land, only includes long-period load tide
            ### internal_tide_hret: model for sea surface height displacement from coherent internal tide, does not include the incoherent tide
            ### internal_tide_sol2: model for sea surface height displacement from coherent internal tide, alternative to internal_tide_hret
            ### mean_sea_surface_cnescls: model for the mean sea surface height above the reference ellipsoid from 2023 CNES/CLS/SIO/DTU hybrid model.
            ### mean_sea_surface_dtu: alternative model from DTU21
            
            ## Default corrections to include: Pole and Solid Earth Tides 
            tide_corrections = ['pole_tide','solid_earth_tide','load_tide_fes']
            L2_expert_corrected,updated_ssh_var = add_tide_corrections_toSSH(L2_expert_corrected,tide_corrections,updated_ssh_var)
            
            ## Add permanent (zero frequency), which is not included in the solid_earth_tide in SWOT
            ### If  the ground positioning software packages that are used to compute coordinates from in-situ surveys likely adopt 
            ### the solid Earth tide model from the IERS Conventions [38], which includes what they refer to as a “permanent deformation”. 
            ### If the positioning software package includes the permanent deformation in the background solid Earth tide model, 
            ### then the computed coordinates exclude that deformation and are referred to as “conventional tide free” values. 
            ### Calculate the permenant tide (deformation) del hpd that should be added to the tide-free (in situ) OR subtracted from the mean-tide (SWOT)
            # L2_expert_corrected = add_permanent_tide(L2_expert_corrected)
            
    
            ############################################
            ## Calculate WSE as SSH_corrected - geoid
            final_corrections = [geoid]#,'delhpd']
            L2_expert_corrected = add_WSE_LRL2(L2_expert_corrected,updated_ssh_var,final_corrections)    
    
            
        except:
            print('[LR L2 SSH EXPERT] Filters gave 0 results')
        else:
            if (L2_expert_corrected.num_lines.size>0) & (L2_expert_corrected.num_pixels.size>0): 
                ## Save to a new dataframe
                save_L2_expert = pd.DataFrame({'x':L2_expert_corrected.longitude.values.flatten(),'y':L2_expert_corrected.latitude.values.flatten()})        
                for variable in variables_to_save:
                    if len(L2_expert_corrected[variable].values.shape)==2:                
                        save_L2_expert[variable] = L2_expert_corrected[variable].values.flatten()
                ## Remove bad data
                # save_L2_expert = save_L2_expert[~np.isnan(save_L2_expert['ssh_karin_2'])]
                
                ## Correct the longitude coordinate to be -180 to 180 and save
                if len(save_L2_expert)>0:
                    save_L2_expert = fix_latlon(save_L2_expert)
                    savedf(save_L2_expert,L2_output_folder,l2_exp)             
                    L2_expert_corrected.to_netcdf(L2_output_folder / ('%s.nc' %( l2_exp)))
                    print('[LR L2 SSH EXPERT] Saving to Geojson and NetCDF') 
                                              
                ## Make Plots
                plot1variable_LR_L2orL3_exp(L2_expert_corrected,'L2 LR Expert',minsig0,maxsig0,'sig0_karin' + use_ssh,'viridis',ax[0],area2)
                plot1variable_LR_L2orL3_exp(L2_expert_corrected,'L2 LR Expert', minssh,maxssh,'ssh_karin' + use_ssh,'viridis',ax2[0],area2)
                plot1variable_LR_L2orL3_exp(L2_expert_corrected,'L2 LR Expert', minwse,maxwse,'wse','viridis',ax3[0],area2)

                
                ######################################################################                                                                                                                                                                                                                                                                                                                                                                                                             
                ######################################################################                                                                                                                                                                                                                                                                                                                                                                                                             
                ########## Process L2 Unsmoothed (250m)
                ######################################################################                                                                                                                                                                                                                                                                                                                                                                                                             
                ######################################################################                                                                                                                                                                                                                                                                                                                                                                                                             
                print('[LR L2 SSH UNSMOOTHED]')
                print('[LR L2 SSH UNSMOOTHED] File: ', l2_uns_file)
                L2_uns_right = xr.open_dataset(l2_uns_file, group='right',chunks={'num_lines': 10000})
                L2_uns_left = xr.open_dataset(l2_uns_file, group='left',chunks={'num_lines': 10000})
        
                try:
                    ## Filter pixels
                    ## Crop to area
                    L2_uns_right = crop(L2_uns_right,area2)
                    L2_uns_left = crop(L2_uns_left,area2)
            
                    filters = [['ancillary_surface_classification_flag',0.,0.]]
                    L2_uns_right = other_filters_range(L2_uns_right,filters_class,True,False)
                    L2_uns_left = other_filters_range(L2_uns_left,filters_class,True,False)

                    ## Interpolate missing variables using the matching Expert files
                    variables_to_interpret = [x for x in variables_to_save if (x not in list(L2_uns_right.keys())) & (x not in variables_to_add) ]
                    L2_uns_right_filtered = interpolate_2km_to_250m(L2_uns_right,L2_expert,variables_to_interpret + ['egm08_meantide'])
                    L2_uns_left_filtered = interpolate_2km_to_250m(L2_uns_left,L2_expert,variables_to_interpret + ['egm08_meantide'])

                    L2_uns_right_filtered = other_filters_range(L2_uns_right_filtered,filters_crosstrackl2,False,False)
                    L2_uns_left_filtered = other_filters_range(L2_uns_left_filtered,filters_crosstrackl2,False,False)
                    
                    ## Switch to a tide-free geoid
                    if len(geoid_file)>0:
                        geoid = 'egm08_tidefree'
                        L2_uns_right = convert_geoid_fromMeanTide_toTideFree(L2_uns_right_filtered,geoid_file) 
                        L2_uns_left = convert_geoid_fromMeanTide_toTideFree(L2_uns_left_filtered,geoid_file) 
                    
                    else:
                        print('[LR L2 SSH UNSMOOTHED] The tide-free geoid %s was not found, default geoid within the SWOT products will be used' %(geoid_name))
                        geoid = 'egm08_meantide'
            
            
                    # L2_uns_right_filtered =  get_qual_flag_mask(L2_uns_right_filtered,qual_flag,['bad']) 
                    # L2_uns_left_filtered =  get_qual_flag_mask(L2_uns_left_filtered,qual_flag,['bad']) 
                    
                    
                    # ## Tide Corrections
                    # ## There are discrepancies in the calculation for SSHA for Version C and D.         
                    # ## The 2025 PDD https://deotb6e7tfubr.cloudfront.net/s3-edaf5da92e0ce48fb61175c28b67e95d/podaac-ops-cumulus-docs.s3.us-west-2.amazonaws.com/web-misc/swot_mission_docs/pdd/D-56407_SWOT_Product_Description_L2_LR_SSH_20250224a_RevC_clean_sig.pdf?A-userid=None&Expires=1754431195&Signature=ZFxZ5qTr1lySfoUyjxA7lASIlG9dyhkcCuI3OayJNlYhoDRjNYX4eIkq238HcAKU4R0TZ6Aoi1k4YAJO3yRC1MXJjDjwz-0xtuNlI84oEVArwJ-EYBNpX3aLiOBSvrjGZDa09mJq87FiKl32ohcnH73rBkx58IVJQ~PKwZyLGmpzmsQVmaGVz~VaZ8nq7G0QGsITYxw6jj2HRdO-XquWiqEJfaSC~A4a4r8UINrp~am8L7s6onBFKhkjpiISmTINA88ADIqfiak8HW-B9IrxSbMsOjLgOZVg~Zulore2EoLziy9Bt8jvcNMwb01aVt5PQUPumab0hVOFFLFTpAxzSQ__&Key-Pair-Id=K2MYOFHSV8YFSM
                    # ## Using their calculation, the SSHA in Version C is not correct when calculated from SSH. There is a +/- 0.003 difference
                    # ## get SSH from SSHA by removing corrections for DAC, MSS, Ocean tide (fes, non eq, and eq), Internal tide (hret)
                    # if version=='C':
                    #     L2_uns_right_filtered = add_SSH_fromSSHA_VC(L2_uns_right_filtered,use_ssh)
                    #     L2_uns_left_filtered = add_SSH_fromSSHA_VC(L2_uns_left_filtered,use_ssh)
                    # elif version =='D':
                    #     L2_uns_right_filtered = add_SSH_fromSSHA_VD(L2_uns_right_filtered,use_ssh)
                    #     L2_uns_left_filtered = add_SSH_fromSSHA_VD(L2_uns_left_filtered,use_ssh)
                        
                    ############################################
                    ## Apply crossover correction
                    L2_uns_right_corrected,updated_ssh_varR = add_crossover_correction_toSSH(L2_uns_right_filtered,use_ssh)
                    L2_uns_left_corrected,updated_ssh_varL = add_crossover_correction_toSSH(L2_uns_left_filtered,use_ssh)

                    ############################################
                    ## Tide Corrections                    
                    L2_uns_right_corrected,updated_ssh_varR = add_tide_corrections_toSSH(L2_uns_right_corrected,tide_corrections,updated_ssh_varR)
                    L2_uns_left_corrected,updated_ssh_varL = add_tide_corrections_toSSH(L2_uns_left_corrected,tide_corrections,updated_ssh_varL)
                    
                    ## Add permanent (zero frequency), which is not included in the solid_earth_tide in SWOT
                    L2_uns_right_corrected = add_permanent_tide(L2_uns_right_corrected)
                    L2_uns_left_corrected = add_permanent_tide(L2_uns_left_corrected)
                    
                    
                    ############################################
                    ## Calculate WSE as SSH_corrected - geoid
                    L2_uns_right_corrected = add_WSE_LRL2(L2_uns_right_corrected,updated_ssh_var,final_corrections)    
                    L2_uns_left_corrected = add_WSE_LRL2(L2_uns_left_corrected,updated_ssh_var,final_corrections)    
            
                except:
                    print('[LR L2 SSH UNSMOOTHED] Filters gave 0 results')
                else:
                    if (L2_uns_right_corrected.num_lines.size>0) | (L2_uns_right_corrected.num_pixels.size>0) | (L2_uns_left_corrected.num_lines.size>0) & (L2_uns_left_corrected.num_pixels.size>0):
                        ## Save to GeoDataFrame that will be exported as geojson
                        save_L2_uns = pd.DataFrame({'x':np.append(L2_uns_right_corrected.longitude,L2_uns_left_corrected.longitude),
                                              'y':np.append(L2_uns_right_corrected.latitude,L2_uns_left_corrected.latitude)})        
                        for variable in variables_to_save:
                            save_L2_uns[variable] = np.append(L2_uns_right_corrected[variable].values,L2_uns_left_corrected[variable].values)
                        ## Remove bad dates
                        save_L2_uns = save_L2_uns[~np.isnan(save_L2_uns['ssh_karin_2'])]
                        
                        ## Correct the longitude coordinate to be -180 to 180 and save
                        if len(save_L2_uns)>0:
                            save_L2_uns = fix_latlon(save_L2_uns)
                            savedf(save_L2_uns,L2_output_folder,l2_uns)   
                            # L2_uns_right_corrected.to_netcdf(L2_output_folder / ('%s_right.nc' %( l2_uns_file)))
                            # L2_uns_left_corrected.to_netcdf(L2_output_folder / ('%s_left.nc' %( l2_uns_file)))
                            print('[LR L2 SSH UNSMOOTHED] Saving to Geojson only') 
                            
                        ## Make Plots
                        plot1variable_LR_L2orL3_uns(L2_uns_right_corrected,'L2 LR Unsmoothed', L2_uns_left_corrected,minsig0,maxsig0,'sig0_karin' + use_ssh,'viridis',ax[1],area2)
                        plot1variable_LR_L2orL3_uns(L2_uns_right_corrected,'L2 LR Unsmoothed', L2_uns_left_corrected,minssh,maxssh,'ssh_karin' + use_ssh,'viridis',ax2[1],area2)
                        plot1variable_LR_L2orL3_uns(L2_uns_right_corrected,'L2 LR Unsmoothed', L2_uns_left_corrected,minwse,maxwse,'wse','viridis',ax3[1],area2)
                
            
                        ############################################
                        ############################################
                        ############################################
                        ## L3 Expert File
                        print('[LR L3 SSH EXPERT]')
                        l3_exp = '%sExpert%s_lat%s-%s_%s' %(prefix.replace('L2','L3'),suffix.replace(l2_processing_full,l3version),np.round(area[1],1),np.round(area[3],1),aoi)
                        try:
                            l3_exp_file = glob.glob(str(L3_folder / ('%sExpert*%s_%s*%s*%s*.nc' %(prefix.replace('L2','L3'),passs,date1,date2,l3version))))[-1]
                        except:
                            print('[LR L3 SSH EXPERT] L3 Expert file not downloaded yet for pass: %s cycle: %s date: %s-%s' %(passs, cycle, date1,date2))
                        else:
                            print('[LR L3 SSH EXPERT] File: ', l3_exp_file)
                            ######################################################################                                                                                                                                                                                                                                                                                                                                                                                                             
                            ## Process L3 Expert (2km)
                            ######################################################################                                                                                                                                                                                                                                                                                                                                                                                                             
                            ######################################################################                                                                                                                                                                                                                                                                                                                                                                                                             
                            L3_expert = xr.open_dataset(l3_exp_file)
                
                            variables_to_save = list(L3_expert.keys()) 
                            variables_to_remove = ['i_num_line','i_num_pixel']
                            variables_to_add = ['wse','egm08_meantide','egm08_tidefree']
                            variables_to_save = variables_to_save + variables_to_add
                            variables_to_save = [x for x in variables_to_save if x not in variables_to_remove]
                        
                            
                            ## get name of ssha that is not filtered
                            if l3version == 'v2.0.1':
                                ssha_var = 'ssha_unedited'
                            else:
                                ssha_var = 'ssha_noiseless'
                    
                            ## Crop to area 
                            L3_expert = crop(L3_expert,area2)
                              
                                                 ## Add geoid variable from L2 product
                            L3_expert['egm08_meantide'] = L2_expert['egm08_meantide']
                            ## Switch to a tide-free geoid
                            if len(geoid_file)>0:
                                geoid = 'egm08_tidefree'
                                L3_expert = convert_geoid_fromMeanTide_toTideFree(L3_expert,geoid_file) 
                            
                            else:
                                print('[LR L3 SSH EXPERT] The tide-free geoid %s was not found, default geoid within the SWOT products will be used')
                                geoid = 'egm08_meantide'
                
                            try:
                                L3_expert_filtered = other_filters_exclude(L3_expert,l3_exclusions)                    
                                L3_expert_filtered = other_filters_range(L3_expert_filtered,filters_crosstrackl3,False,False)
    
                               
                                ############################################
                                ## Calculate WSE as SSH_corrected - geoid
                                corrections_to_remove = ['dac','mss','ocean_tide','internal_tide']
                                L3_expert_corrected = add_wse_LRL3(L3_expert_filtered,ssha_var,geoid,corrections_to_remove)    
                                
                                
                            except:
                                print('[LR L3 SSH EXPERT] Filter gave 0 results')
                            else:
                                if (L3_expert_corrected.num_lines.size>0) | (L3_expert_corrected.num_pixels.size>0) :
                                    ## Save to GeoDataFrame that will be exported as geojson
                                    save_L3_expert = pd.DataFrame({'x':L3_expert_corrected.longitude.values.flatten(),
                                                          'y':L3_expert_corrected.latitude.values.flatten()})
                            
                                    for variable in variables_to_save:
                                        save_L3_expert[variable] = L3_expert_corrected[variable].values.flatten()
                                    ## Remove bad data
                                    save_L3_expert = save_L3_expert[~np.isnan(save_L3_expert['ssha_unedited'])]
                                    ## Correct the longitude coordinate to be -180 to 180 and save
                                    if len(save_L3_expert)>0:
                                        save_L3_expert = fix_latlon(save_L3_expert)
                                        savedf(save_L3_expert,L3_output_folder,l3_exp)   
                                        L3_expert_corrected.to_netcdf(L3_output_folder / ('%s.nc' %( l3_exp)))
                                        print('[LR L3 SSH EXPERT] Saving to Geojson and NetCDF') 
                                        
                                    ## Make Plots
                                    plot1variable_LR_L2orL3_exp(L3_expert_corrected,'L3 LR Expert', minwse,maxwse,'wse','viridis',ax3[2],area2)
                                
                                 
                                
                                    ######################################################################                                                                                                                                                                                                                                                                                                                                                                                                             
                                    ######################################################################                                                                                                                                                                                                                                                                                                                                                                                                             
                                    ## Process L3 Unsmoothed (250m)
                                    ######################################################################                                                                                                                                                                                                                                                                                                                                                                                                             
                                    ######################################################################                                                                                                                                                                                                                                                                                                                                                                                                             
                                    ## L3 Unsmoothed File
                                    print('[LR L3 SSH UNSMOOTHED]')
                                    l3_uns    = '%sUnsmoothed%s_lat%s-%s_%s' %(prefix.replace('L2','L3'),suffix.replace(l2_processing_full,l3version),np.round(area[1],1),np.round(area[3],1),aoi)
                                    try:
                                        l3_uns_file = glob.glob(str(L3_folder / ('%sUnsmoothed*%s_%s*%s*%s*.nc' %(prefix.replace('L2','L3'),passs,date1,date2,l3version))))[-1]
                                    except:
                                        print('[LR L3 SSH UNSMOOTHED] L3 Unsmoothed file not downloaded yet for pass: %s cycle: %s date: %s-%s' %(passs, cycle, date1,date2))
                                    else:
                                        print('[LR L3 SSH UNSMOOTHED] File: ',l3_uns_file)
                                        
                                        L3_uns = xr.open_dataset(l3_uns_file)
                                    
                                        ## Crop to area 
                                        L3_uns = crop(L3_uns,area2)
                                        
                                        if l3version != 'v2.0.1': 
                                            variables_to_interpret = [x for x in variables_to_save if (x not in list(L3_uns.keys())) & (x not in variables_to_add) ]
                                            L3_uns = interpolate_2km_to_250m(L3_uns,L3_expert,variables_to_interpret)
                            
                                        
                                        ## Add geoid variable from L2 product
                                        L3_uns = interpolate_2km_to_250m(L3_uns,L3_expert,['egm08_meantide'])
                                        if len(geoid_file)>0:
                                            geoid = 'egm08_tidefree'
                                            L3_uns = convert_geoid_fromMeanTide_toTideFree(L3_uns,geoid_file) 
                                        
                                        else:
                                            print('[LR L3 SSH UNSMOOTHED] The tide-free geoid %s was not found, default geoid within the SWOT products will be used')
                                            geoid = 'egm08_meantide'
                            
                                     
                                        try:
                                            L3_uns_filtered = other_filters_exclude(L3_uns,l3_exclusions)
                                            L3_uns_filtered = other_filters_range(L3_uns_filtered,filters_crosstrackl3,False,False)
                                        
                                            ## Add WSE variable
                                            L3_uns_corrected = add_wse_LRL3(L3_uns_filtered,ssha_var,geoid,corrections_to_remove)    
                                
                                        except:
                                            print('[LR L3 SSH UNSMOOTHED] Filter gave 0 results')
                                        else:
                                            if (L3_uns_corrected.num_lines.size>0) | (L3_uns_corrected.num_pixels.size>0) :
                                                ## Save to GeoDataFrame that will be exported as geojson
                                                save_l3_uns = pd.DataFrame({'x':L3_uns_corrected.longitude.values.flatten(),
                                                                      'y':L3_uns_corrected.latitude.values.flatten()})    
                                                for variable in variables_to_save:
                                                    save_l3_uns[variable] = L3_uns_corrected[variable].values.flatten()
                                                ## Remove bad dates
                                                save_l3_uns = save_l3_uns[~np.isnan(save_l3_uns[ssha_var])]
                                                ## Correct the longitude coordinate to be -180 to 180 and save
                                                if len(save_l3_uns)>0:
                                                    save_l3_uns = fix_latlon(save_l3_uns)
                                                    savedf(save_l3_uns,L3_output_folder,l3_uns)   
                                                    L3_uns_corrected.to_netcdf(L3_output_folder / ('%s.nc' %( l3_uns)))
                                                    print('[LR L3 SSH UNSMOOTHED] Saving to Geojson and NetCDF') 
                                                ## Make Plots
                                                plot1variable_LR_L2orL3_exp(L3_uns_corrected,'L3 LR Unsmoothed', minwse,maxwse,'wse','viridis',ax3[3],area2)
                                        
                                
                        
                
                fig.suptitle(date1 + ': sig0_karin'+ use_ssh)
                fig.tight_layout()
                fig2.suptitle(date1 + ': ssh_karin' + use_ssh)
                fig2.tight_layout()
                fig3.suptitle(date1 + ': Water Surface Elevation')
                fig3.tight_layout()
                
                fig.savefig(check1,dpi=200)
                fig2.savefig(check2,dpi=200)
                fig3.savefig(check3,dpi=200)
    else:
        print('Already processed')
    
        
           
print('Done')
# shutil.rmtree(tmp_folder)





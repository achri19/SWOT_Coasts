#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Utilitiy Toolbox

import numpy as np
import scipy
import xarray as xr
import geopandas as gpd
import pandas as pd
from geopy import distance
import math
import netCDF4 as nc


def crop(ds,area):
    dim_0_name = ds['latitude'].dims[0]
    dim_1_name = ds['latitude'].dims[1]
    mask = (
        (ds['latitude'] >= area[1]) & (ds['latitude'] <= area[3]) &
        (ds['longitude'] >= area[0]) & (ds['longitude'] <= area[2])
        )
    mask_for_dim_0 = mask.any(dim=dim_1_name)
    mask_for_dim_1 = mask.any(dim=dim_0_name)
    mask_for_dim_0 = mask_for_dim_0.compute()
    mask_for_dim_1 = mask_for_dim_1.compute()
    mask_for_dim_0 = mask_for_dim_0.compute()
    mask_for_dim_1 = mask_for_dim_1.compute()
    ds_cropped = ds.sel({
        dim_0_name: mask_for_dim_0,
        dim_1_name: mask_for_dim_1
    })
    
    print(f"Original size: {dict(ds.dims)}")
    print(f"Cropped size:  {dict(ds_cropped.dims)}")
    
    return ds_cropped 


def other_filters_range(ds,filters,include=True,drop=True):
    print('\tNC dimensions: Lines = %s x Pixels = %s = %s' %(ds.dims['num_lines'],ds.dims['num_pixels'],ds.dims['num_lines'] * ds.dims['num_pixels']))
    for i in range(len(filters)): 
        var = filters[i][0]
        vmin = filters[i][1]
        vmax = filters[i][2]
        print('\n\t\tFilter by %s <= %s <= %s Inclusive = %s Drop = %s' %(vmin,var,vmax,include,drop)) 
        print('\t\tBefore: %s' %(np.count_nonzero(~np.isnan(ds[var]))))
        selection = (ds[var] >= vmin) & (ds[var] <= vmax)
        selection = selection.compute() 
                      
        try:
            if include:
                ds = ds.where(selection, drop=drop)  
            else:
                ds = ds.where(~selection,drop=drop)
            print('\t\tAfter: %s' %(np.count_nonzero(~np.isnan(ds[var]))))
        except:
            print('\t\tBad Filter')
            ds = xr.Dataset()
    print('\t\tNC dimensions: Lines (%s) x Pixels (%s) = %s' %(ds.dims['num_lines'],ds.dims['num_pixels'],ds.dims['num_lines'] * ds.dims['num_pixels']))
            
    return ds


def other_filters_exclude(ds,filters):
    for i in range(len(filters)): 
        
        var = filters[i][0]
        exclude = filters[i][1]
        print('\n\tFilter by %s != %s' %(var,exclude)) 
        print('\t\tBefore: %s' %(np.count_nonzero(~np.isnan(ds[var]))))
        selection = ~ds[var].isin(exclude)# 2684354816))
        selection = selection.compute()
        try:
            ds = ds.where(selection, drop=False)  
            print('\t\tAfter: %s' %(np.count_nonzero(~np.isnan(ds[var]))))
        except:
            print('\t\tBad Filter')
            ds = xr.Dataset()
            
    return ds 




def add_crossover_correction_toSSH(ds,use_ssh):
    print('\n\tAdd SSH Corrections')
    ssh = ds['ssh_karin' +use_ssh]
    ds['ssh_xover'] = (('num_lines','num_pixels'),ssh.values + ds.height_cor_xover.values)
    ds.ssh_xover.attrs = ssh.attrs
    ds.ssh_xover.attrs['long_name'] = 'SSH_xover'
    ds.ssh_xover.attrs['standard_name'] = 'SSH corrected with crossover calibration estimate'
    ds.ssh_xover.attrs['valid_min'] = ssh.attrs['valid_min'] + ds.height_cor_xover.attrs['valid_min']
    ds.ssh_xover.attrs['valid_max'] = ssh.attrs['valid_max'] + ds.height_cor_xover.attrs['valid_max'] 
    ds.ssh_xover.attrs['comment'] = 'ssh_xover is calculated as ssh_karin%s + height_cor_xover ' %(use_ssh)
    print('\t\t %s' %(ds.ssh_xover.attrs['comment']))
    return ds,'ssh_xover'

def add_tide_corrections_toPIXC(ds,tides,updated_ssh_var):
    print('\n\tAdd SSH Corrections')
    ssh = ds[updated_ssh_var]
    tmp = ssh.values
    txt = updated_ssh_var
    for tide in tides:
        tmp = tmp  - ds[tide].values
        txt = txt + ' - ' + tide

    ds['ht_tides'] = (('points'),tmp)
    ds.ht_tides.attrs = ssh.attrs
    ds.ht_tides.attrs['long_name'] = 'tide_corrections_applied'
    ds.ht_tides.attrs['standard_name'] = 'SSH corrected for tide corrections: %s' %(txt)
    ds.ht_tides.attrs['valid_min'] = ssh.attrs['valid_min'] 
    ds.ht_tides.attrs['valid_max'] = ssh.attrs['valid_max'] 
    ds.ht_tides.attrs['comment'] = 'ssh_tides is calculated as %s' %(txt)
    print('\t\t %s ' %(ds.ht_tides.attrs['comment']))
    return ds,'ht_tides'


def add_tide_corrections_toSSH(ds,tides,updated_ssh_var):
    print('\n\tAdd SSH Corrections')
    ssh = ds[updated_ssh_var]
    tmp = ssh.values
    txt = updated_ssh_var
    for tide in tides:
        tmp = tmp  - ds[tide].values
        txt = txt + ' minus ' + tide

    ds['ssh_tides'] = (('num_lines','num_pixels'),tmp)
    ds.ssh_tides.attrs = ssh.attrs
    ds.ssh_tides.attrs['long_name'] = 'SSH_tide_corrections_applied'
    ds.ssh_tides.attrs['standard_name'] = 'SSH corrected for tide corrections: %s' %(txt)
    ds.ssh_tides.attrs['valid_min'] = ssh.attrs['valid_min'] 
    ds.ssh_tides.attrs['valid_max'] = ssh.attrs['valid_max'] 
    ds.ssh_tides.attrs['comment'] = 'ssh_tides is calculated as %s' %(txt)
    print('\t\t %s' %(ds.ssh_tides.attrs['comment']))
    return ds,'ssh_tides'



def calculate_delh_MeantTide_minus_TideFree_geoidheights(ds):
    print('\n\t\tCalculate conversion from SWOT mean-tide geoid height to tide-free geoid heights based on equation 11.4 and 11.5 in User Handbook')
    ## From March 2025 handbook page 186
    ## SWOT geoid is provided in mean-tide system adn includes the zero-frequency (time-independent) permanent tide height. 
    ## The difference (del_h) between mean-tide (hmt) and tide-free (htf) geoid heights is a function of latitude
    ## ∆h = hmt − htf [eq 11.4 in handbook)= (1 + k2) * Hperm * sqrt(5/4pi) * (3/2 * sin2(lat deg) - 1/2) [eq 11.5 in handbook]
    hperm = -0.31460
    k2 = 0.3 #(Love number, change in the gravitational potential caused by tidal deformation)
    ds['geoid_delh'] = (1+k2)*hperm*np.sqrt(5/(4*np.pi))*((3/2)*np.sin(ds['latitude'])**2-(1/2))
    ds.geoid_delh.attrs['long_name'] = 'mean-tide - tide-free geoid height difference'
    ds.geoid_delh.attrs['standard_name'] = 'del h '
    ds.geoid_delh.attrs['comment'] = 'SWOT Handbook, Section 11.3.1, Equation 11.4 and 11.5: delh is calculated as hmt - htf, equal to (1+k2)Hperm * sqrt(5/4pi)(1.5sin2(latdeg)-0.5) where k2 = %s and hperm = %s' %(k2, hperm)
    # print('\t\t%s' %(ds.geoid_delh.attrs['comment']))
    return ds

def convert_geoid_fromMeanTide_toTideFreePIXC(ds,geoid_file):
    print('\n\tAdd geoid')
    
    if len(geoid_file)>0:
        ## An alternative geoid is provided
        ## The tide-free geoid heights were already calculated by Michael
        ## Pixel offset errors were also corrected in his geoid
        
        ## Use meantide_egm08 (must be a DataArray) to select nearest points
        # ds['egm08_tidefree'] = geoid['geoid'].sel(lat=ds['latitude'],lon=ds['longitude'],method='nearest')
        print('\n\t\tA tide-free geoid is provided as an ancillary nc file' %(geoid_file[0]))
        ds['egm08_tidefree'] = (('points'),interp_nc_geoid_at_coords(ds['longitude'],ds['latitude'],geoid_file[0]))
    
    else:
        ## Must calculated tide-free geoid using the geoid values within the SWOT product
        print('\n\t\tThe tide-free geoid was not found, tide-free geoid heights will be calculated with the geoid included in the SWOT product')
        ds['egm08_tidefree'] = ds['egm08_meantide'] - ds['geoid_delh']

    ds.egm08_tidefree.attrs = ds.egm08_meantide.attrs
    ds.egm08_tidefree.attrs['comment'] = ds.egm08_meantide.attrs['comment'].replace('mean tide','tide free')
    return ds

def convert_geoid_fromMeanTide_toTideFree(ds,geoid_file):
    print('\n\t\tAdd geoid')
    
    if len(geoid_file)>0:
        ## An alternative geoid is provided
        ## The tide-free geoid heights were already calculated by Michael
        ## Pixel offset errors were also corrected in his geoid
        
        ## Use meantide_egm08 (must be a DataArray) to select nearest points
        # ds['egm08_tidefree'] = geoid['geoid'].sel(lat=ds['latitude'],lon=ds['longitude'],method='nearest')
        print('\n\t\tA tide-free geoid, %s, is provided as an ancillary nc file' %(geoid_file[0].split('/')[-1]))
        ds['egm08_tidefree'] = (('num_lines','num_pixels'),interp_nc_geoid_at_coords(ds['longitude'],ds['latitude'],geoid_file[0]))
    
    else:
        ## Must calculated tide-free geoid using the geoid values within the SWOT product
        print('\n\t\tThe tide-free geoid %s was not found, tide-free geoid heights will be calculated with the geoid included in the SWOT product')
        ds['egm08_tidefree'] = ds['egm08_meantide'] - ds['geoid_delh']

    ds.egm08_tidefree.attrs = ds.egm08_meantide.attrs
    ds.egm08_tidefree.attrs['comment'] = ds.egm08_meantide.attrs['comment'].replace('mean tide','tide free')
    return ds

def add_permanent_tide(ds):
    ## From March 2025 handbook page 191
    ## As mentioned above, the solid Earth tide model in SWOT products excludes the permanent
    ## tide. In contrast, most ground positioning software packages that are used to compute coor-
    ## dinates from in-situ surveys likely adopt the solid Earth tide model from the IERS Conventions
    ## [38], which includes what they refer to as a “permanent deformation”. 
    ## If the positioning software package includes the permanent deformation in the background solid
    ## Earth tide model, then the computed coordinates exclude that deformation and are referred to
    ## as “conventional tide free” values. In this case, the permanent deformation should be added to
    ## the conventional tide free coordinates before comparing them to SWOT measurements. Specifically,
    ## the following representation of the permanent deformation ∆hpd should be added to the
    ## conventional tide free coordinates, where ∆hpd = h2 * Hperm * sqrt(5/4pi) * (3/2 * sin2(lat deg)-1/2)
    ## hperm = -0.31460 (permanent tide-potential amplitude)
    ## h2 = 0.6078 (Love number, vertical displacement of the body's surface due to tidal forces)
    print('\n\tCalculate permanent deformation from IERS conventions, which should be added to the conventional tide-free coordinates according to eq 11.7 in the User Handbook')
    h2 = 0.6078
    hperm = -0.31460
    ds['delhpd'] = h2*hperm*np.sqrt(5/(4*np.pi))*((3/2*np.sin(ds['latitude'])**2)-(1/2))
    ds.delhpd.attrs['long_name'] = 'permanent deformation'
    ds.delhpd.attrs['standard_name'] = 'permanent deformation'
    ds.delhpd.attrs['comment'] = 'delhpd is calculated as h2*hperm*sqrt(5/4pi)*(1.5sin2(latdeg)-0.5) where h2 = %s and hperm = %s' %(h2, hperm)
    print('\t\t%s' %(ds.delhpd.attrs['comment']))
    return ds
    
def add_WSE_HRPIXC(ds,wse_var,corrections):
    print('\n\tCalculate Water Surface Elevation')
    wse = ds[wse_var]
    tmp = wse.values
    txt = wse_var
    for cor in corrections:
        tmp = tmp  - ds[cor].values
        txt = txt + ' minus ' + cor

    ds['wse_egm08'] = (('points'),tmp)
    ds.wse_egm08.attrs = wse.attrs
    ds.wse_egm08.attrs['long_name'] = 'water surface elevation'
    ds.wse_egm08.attrs['standard_name'] = 'water surface elevation: %s' %(txt)
    ds.wse_egm08.attrs['valid_min'] = wse.attrs['valid_min'] 
    ds.wse_egm08.attrs['valid_max'] = wse.attrs['valid_max'] 
    ds.wse_egm08.attrs['comment'] = 'HR L2 PIXC is calculated as %s ' %(txt)
    print('\t\t%s' %(ds.wse_egm08.attrs['comment']))
    return ds


def add_WSE_HRL2(ds,wse_var,corrections):
    print('\n\t\tCalculate Water Surface Elevation')
    wse = ds[wse_var]
    tmp = wse.values
    txt = wse_var
    for cor in corrections:
        tmp = tmp  - ds[cor].values
        txt = txt + ' minus ' + cor

    ds['wse_egm08'] = (('num_lines','num_pixels'),tmp)
    ds.wse_egm08.attrs = wse.attrs
    ds.wse_egm08.attrs['long_name'] = 'water surface elevation'
    ds.wse_egm08.attrs['standard_name'] = 'water surface elevation: %s' %(txt)
    ds.wse_egm08.attrs['valid_min'] = wse.attrs['valid_min'] 
    ds.wse_egm08.attrs['valid_max'] = wse.attrs['valid_max'] 
    ds.wse_egm08.attrs['comment'] = 'HR L2 Raster WSE is calculated as %s ' %(txt)
    print('\t\t\t%s' %(ds.wse_egm08.attrs['comment']))
    return ds
    

def add_WSE_LRL2(ds,updated_ssh_var,corrections):
    print('\n\tCalculate Water Surface Elevation')
    ssh = ds[updated_ssh_var]
    tmp = ssh.values
    txt = updated_ssh_var
    for cor in corrections:
        tmp = tmp  - ds[cor].values
        txt = txt + ' minus ' + cor

    ds['wse'] = (('num_lines','num_pixels'),tmp)
    ds.wse.attrs = ssh.attrs
    ds.wse.attrs['long_name'] = 'water surface elevation'
    ds.wse.attrs['standard_name'] = 'water surface elevation: %s' %(txt)
    ds.wse.attrs['valid_min'] = ssh.attrs['valid_min'] 
    ds.wse.attrs['valid_max'] = ssh.attrs['valid_max'] 
    ds.wse.attrs['comment'] = 'LR L2 WSE is calculated as %s ' %(txt)
    print('\t\t%s' %(ds.wse.attrs['comment']))
    return ds
    
def add_wse_LRL3(ds,ssha_var,geoid,corrections_to_remove):
    print('\n\tCalculate Water Surface Elevation')
    ssha = ds[ssha_var]
    tmp = ssha.values 
    txt = ssha_var
    for cor in corrections_to_remove:
        tmp = tmp  + ds[cor].values
        txt = txt + ' - ' + cor
    tmp = tmp  - ds[geoid].values

    ds['wse'] = (('num_lines','num_pixels'),tmp)
    ds.wse.attrs = ds[ssha_var].attrs
    ds.wse.attrs['long_name'] = 'water surface elevation'
    ds.wse.attrs['standard_name'] = 'water surface elevation: %s' %(txt)
    ds.wse.attrs['comment'] = 'LR L3 WSE is calculated as %s ' %(txt)
    print('\t\t%s' %(ds.wse.attrs['comment']))
    return ds


def interpolate_2km_to_250m(ds_side,ds,variables_to_interpret):    
    print('\nInterpolate from EXPERT 2km to UNSMOOTHED 250m')
    print(variables_to_interpret)
    side_lon = np.ma.masked_array(ds_side.longitude.values)
    side_lat = np.ma.masked_array(ds_side.latitude.values)
    expert_lon = np.ma.masked_array(ds.longitude.values)
    expert_lat = np.ma.masked_array(ds.latitude.values)
    for variable in variables_to_interpret:
        print(variable, ds[variable].values.shape)
        expert = np.ma.masked_array(ds[variable].values)
        ds_side[variable] = (('num_lines','num_pixels'),interp_expert_to_unsmoothed_AlexFore(expert_lon, expert_lat, expert, side_lon, side_lat))
        ds_side[variable].attrs = ds[variable].attrs
        ds_side[variable].attrs['comment'] =  'INTERPOLATED from L2 EXPERT PRODUCT: ' + ds_side[variable].attrs['comment']
       
    return ds_side


def nearest_gauge(items,pivot):
    nearest_ts = min(items,key=lambda x: abs(x-pivot))
    nearest_loc = np.where(items==nearest_ts)[0][0]
    if abs(items.iloc[nearest_loc] - pivot).seconds/60>30:
        # print('No gauge data within 30 minutes')
        nearest_ts = ''
        nearest_loc = ''
    return  nearest_ts, nearest_loc

def func(x, a):
    gradient = 1 # fixed gradient, not optimized
    return gradient * x + a

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculates the distance in km between two points on Earth using the Haversine formula.

    Args:
        lat1: Latitude of the first point in degrees.
        lon1: Longitude of the first point in degrees.
        lat2: Latitude of the second point in degrees.
        lon2: Longitude of the second point in degrees.

    Returns:
        The distance between the two points in kilometers.
    """
    R = 6371  # Radius of the Earth in kilometers
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)
    
    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad
    a = math.sin(dlat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    dist = R * c
    return dist

def fix_latlon(df):
    if df['x'].iloc[0] > 180:
        df['x2'] =df['x']-360
    else:
        df['x2'] =df['x']
    return df

def savedf(df,folder,name):
    gdf = gpd.GeoDataFrame(df,geometry= gpd.points_from_xy(df.x2,df.y))
    gdf.to_file(folder / (name + '.geojson'))       
    
### Code from collaborators at JPL

def interp_expert_to_unsmoothed_AlexFore(expert_lon, expert_lat, expert_var, unsmoothed_lon, unsmoothed_lat):
    """
    Copyright (c) 2023, California Institute of Technology ("Caltech"). U.S.
    Government sponsorship acknowledged.
    All rights reserved.
     
    Author (s): Alex Fore
    """
    
    """
    Interpolates expert_var to unsmoothed grid
    """
    if not (np.ma.isMaskedArray(expert_var) and
            np.ma.isMaskedArray(expert_lon) and
            np.ma.isMaskedArray(expert_lat)):
        raise Exception((
            "Please pass in expert_lon, expert_lat, and expert_var as numpy "
            "mask arrays."))
 
    if not (len(expert_var.shape) == 2 and
            len(expert_lon.shape) == 2 and
            len(expert_lat.shape) == 2):
        raise Exception((
            "Please pass in expert_lon, expert_lat, and expert_var as full "
            "2D arrays as they are in the Expert file."))
 
    mask = ~np.ma.getmaskarray(expert_var)
    var_in = expert_var[mask].data
    lon_in = expert_lon[mask].data
    lat_in = expert_lat[mask].data
 
    unsmoothed_mask = np.logical_and(
        ~np.ma.getmaskarray(unsmoothed_lon),
        ~np.ma.getmaskarray(unsmoothed_lat))
 
    # Center longitude on median longitude
    center_lon = np.median(lon_in)
    lon_in -= center_lon
    unsmoothed_lon_in = unsmoothed_lon[unsmoothed_mask] - center_lon
 
    unsmoothed_var = np.ones(unsmoothed_lon.shape) * np.nan
    unsmoothed_var[unsmoothed_mask] = scipy.interpolate.griddata(
        (lon_in, lat_in), var_in,
        (unsmoothed_lon_in, unsmoothed_lat[unsmoothed_mask]),
        method='linear', fill_value=np.nan)
 
    return unsmoothed_var.astype(np.float32)



def bilinear_interpolate(data, x, y):
    ## Function by Michael Denbina, California Institute of Technology, Jet Propulsion Laboratory
    """Function to perform bilinear interpolation on the input array data, at
    the image coordinates given by input arguments x and y.

    Arguments
        data (array): 2D array containing raster data to interpolate.
        x (array): the X coordinate values at which to interpolate (in array
            indices, starting at zero).  Note that X refers to the second
            dimension of data (e.g., the columns).
        y (array): the Y coordinate values at which to interpolate (in array
            indices, starting at zero).  Note that Y refers to the first
            dimension of data (e.g., the rows).

    Returns:
        intdata (array): The 2D interpolated array, with same dimensions as
            x and y.

    """
    x = np.asarray(x)
    y = np.asarray(y)

    # Get lower and upper bounds for each pixel.
    x0 = np.floor(x).astype(int)
    x1 = x0 + 1
    y0 = np.floor(y).astype(int)
    y1 = y0 + 1

    # Clip the image coordinates to the size of the input data.
    x0 = np.clip(x0, 0, data.shape[1]-1)
    x1 = np.clip(x1, 0, data.shape[1]-1)
    y0 = np.clip(y0, 0, data.shape[0]-1)
    y1 = np.clip(y1, 0, data.shape[0]-1)

    data_ll = data[ y0, x0 ] # lower left corner image values
    data_ul = data[ y1, x0 ] # upper left corner image values
    data_lr = data[ y0, x1 ] # lower right corner image values
    data_ur = data[ y1, x1 ] # upper right corner image values

    w_ll = (x1-x) * (y1-y) # weight for lower left value
    w_ul = (x1-x) * (y-y0) # weight for upper left value
    w_lr = (x-x0) * (y1-y) # weight for lower right value
    w_ur = (x-x0) * (y-y0) # weight for upper right value

    # Where the x or y coordinates are outside of the image boundaries, set all
    # of the weights to nan, so that these values are nan in the output array.
    ind = (x < 0) | (x > data.shape[1]-1) | (y < 0) | (y > data.shape[0]-1)
    if np.any(ind):
        w_ll[ind] = np.nan
        w_ul[ind] = np.nan
        w_lr[ind] = np.nan
        w_ur[ind] = np.nan
        data_ll[ind] = np.nan
        data_ul[ind] = np.nan
        data_lr[ind] = np.nan
        data_ur[ind] = np.nan

    # Handle void data.
    ind_all_void = np.isnan(data_ll) & np.isnan(data_ul) & np.isnan(data_lr) & np.isnan(data_ur)
    ind_some_void = np.isnan(data_ll) & (~ind_all_void)
    if np.any(ind_some_void):
        w_ll[ind_some_void] = 0
        data_ll[ind_some_void] = 0
    ind_some_void = np.isnan(data_ul) & (~ind_all_void)
    if np.any(ind_some_void):
        w_ul[ind_some_void] = 0
        data_ul[ind_some_void] = 0
    ind_some_void = np.isnan(data_lr) & (~ind_all_void)
    if np.any(ind_some_void):
        w_lr[ind_some_void] = 0
        data_lr[ind_some_void] = 0
    ind_some_void = np.isnan(data_ur) & (~ind_all_void)
    if np.any(ind_some_void):
        w_ur[ind_some_void] = 0
        data_ur[ind_some_void] = 0

    intdata = w_ll*data_ll + w_ul*data_ul + w_lr*data_lr + w_ur*data_ur

    if np.any(ind_all_void):
        intdata[ind_all_void] = np.nan
    ind_low_weight = (w_ll + w_ul + w_lr + w_ur) < 0.99
    if np.any(ind_low_weight):
        intdata[ind_low_weight] = np.nan

    return intdata


def interp_nc_geoid_at_coords(lon, lat, geoid_file):
    
    ## Function by Michael Denbina, California Institute of Technology, Jet Propulsion Laboratory
    """ Interpolate NetCDF geoid at given coordinates. """
    geoid_ds = nc.Dataset(geoid_file)
    geoid_lat = geoid_ds['lat'][:]
    geoid_lon = geoid_ds['lon'][:]
    geoid_hgt = geoid_ds['geoid'][:]

    lon = np.where(lon<0,lon+360,lon) 
    # lon[lon < 0] += 360 # geoid uses 0-360 degree longitude

    # calculate raster coordinates x, y from lon, lat
    x = (lon - geoid_lon[0]) / (geoid_lon[1] - geoid_lon[0])
    y = (lat - geoid_lat[0]) / (geoid_lat[1] - geoid_lat[0])

    return bilinear_interpolate(geoid_hgt, x, y)

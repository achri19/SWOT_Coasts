#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 09:10:18 2024

@author: alchrist
"""

import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import sys

import geopandas as gpd
import fnmatch
import numpy as np
import pandas as pd
sys.path.append('../src/')
from swot_utils import nearest_gauge
import matplotlib.dates as mdates
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error






## Plot error vs. conditions for all pixels within search radius
def make_conditions_error_plots(df,variables_to_plot,gauge_folder,gauge,mode):
    rows = int(np.sqrt(len(variables_to_plot)))
    cols = int(np.ceil(len(variables_to_plot)/rows))
    fig, axes = plt.subplots(nrows=rows, ncols=cols, sharey=True,figsize=(2*cols,4*rows))    
    axes = axes.ravel()
    
    for i in range(len(variables_to_plot)):
        v = variables_to_plot[i]
    
        # axes[i].scatter(df[v],df['Insitu Gauge WSE']-df['SWOT WSE L2'],c='orange',s=1,label='L2')
        axes[i].scatter(df[v],df['gauge_wse']-df['wse'],c='orange',s=1,label='L2')
        # axes[i].scatter(df[v],df['Insitu Gauge WSE']-df['SWOT WSE L3'],c='b',s=1,label='L3')
        axes[i].set_xlabel(v,fontsize=5)
        axes[i].tick_params(axis='both', which='both', labelsize=5)
        axes[i].xaxis.offsetText.set_fontsize(5)
        # if (i)%cols==0:   
        #     axes[i].set_ylabel('Error (In situ - SWOT WSE) m')
    axes[0].legend()
    fig.text(0.005, 0.5, 'Error (In situ - SWOT WSE) m', fontsize=8,va='center', rotation='vertical')
    plt.suptitle(gauge)
    plt.tight_layout()
    plt.savefig('%s/%s_conditions_%s.png' %( gauge_folder , gauge,mode),dpi=300)



def add_colorbar(ax,cax,shrink,pad,aspect,long_name,area):
    cbar = plt.colorbar(cax, ax=ax,
            shrink=shrink, orientation='horizontal',
            pad=pad, aspect=aspect,label= long_name)
    cbar.ax.xaxis.label.set_wrap(True) 
    ax.add_feature(cfeature.COASTLINE)
    ax.set_extent([area[0],area[2],area[1],area[3]] ,crs=ccrs.PlateCarree())

def add_gauge(ax,gauge_lon,gauge_lat):
    ax.scatter(float(gauge_lon),float(gauge_lat),marker = '*', size=10, color = 'k')

def get_qual_flag_mask(qual,bad_list):    
    ### code modified from https://podaac.github.io/tutorials/notebooks/datasets/SWOT_quality_flag_demo.html
    atts = dict(qual.attrs.copy())
    flags = pd.DataFrame({'flag_meanings': atts.get("flag_meanings").split(" "), 
                      'flag_masks': atts.get("flag_masks"), }).set_index("flag_meanings")
    unq_ints_qual = sorted(list(set(qual.load().data.astype(np.int64).flatten().tolist())))
    ## Decompose bit flags
    def qual_bits_iter(n):
        while n:
            b = n & (~n+1)
            yield b
            n ^= b
    
    decomposed_qual = {}
    for i in sorted(unq_ints_qual):
        if i not in decomposed_qual:
            decomposed_qual[i] = qual_bits_iter(i)
    
    new_columns_qual = {}
    for i, j in enumerate(list(set(decomposed_qual))):
        new_columns_qual[j] = flags['flag_masks'] \
            .apply(lambda x: x in [b for b in qual_bits_iter(j)])
    
    lookup = flags.join(pd.DataFrame(new_columns_qual).sort_index(axis=1)).copy()
    flags_to_mask = sorted([i for i in lookup.index if i.startswith(tuple(bad_list))])
    tmp = lookup.copy()  # Make a copy of the table.
    
    del tmp['flag_masks']  # Remove the 'flag_masks' column.
    
    indices = tmp.loc[flags_to_mask].apply(lambda x: x.any()).copy()  # Get a boolean array of integers that comprise the mask.
    
    bad_integers = sorted(tmp.columns[indices].astype(np.int64).tolist())  # Get the columns/integers with one or more flag set to True.
    
    return bad_integers, flags_to_mask
    
def plot_LR_expert_qualityflags(bad_flags,flags_to_mask,qual_flag, fig_folder,o_2000,date,ds,area,variables):
    fig, axes = plt.subplots(2,int(np.ceil(len(variables)/2)), figsize=(15, 12),subplot_kw=dict(projection=ccrs.PlateCarree()))
    axes = axes.ravel()
    shrink = 0.9
    pad = 0.1
    aspect = 20

    selection = ~ds[qual_flag].isin(bad_flags)# 2684354816))
    ds = ds.where(selection, drop=True)  
    
    i=0
    for variable in variables['variable']:
        cax =  axes[i].pcolormesh(ds.longitude, ds.latitude, ds[variables.iloc[i]['variable']],cmap='plasma',vmin=variables['min'].iloc[i],vmax=variables['max'].iloc[i])
        add_colorbar(axes[i],cax,shrink,pad,aspect,ds[variables.iloc[i]['variable']].long_name,area)
        i=i+1
    
    plt.suptitle("%s: %s " %('_'.join(o_2000.split('_')[:4]), date),fontsize=20)
    plt.text(x=0.5, y=0.9, s= "Mask with %s: %s" %(qual_flag, flags_to_mask), fontsize=12, ha="center", transform=fig.transFigure)
    plt.subplots_adjust(left=0.05, right=0.98, top=0.88, bottom=0.05)

    plt.savefig(fig_folder / (o_2000 + '.png'))
    plt.close()


def plot_expert(fig_folder,o_2000,date,ds,area,minssh,maxssh,minwse,maxwse,avg_swot_time,gauge_folder,gauge_name,wl_var,date_var):
    fig, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2, 3, figsize=(15, 12),subplot_kw=dict(projection=ccrs.PlateCarree()))
    # gauge_name = 'Trois-Pistoles'
    gauge = [os.path.join(dirpath,f)
                for dirpath,dirnames, files in os.walk(gauge_folder)
                for f in fnmatch.filter(files,'*%s.csv' %(gauge_name))][0]
    df = gpd.read_file(gauge)
    df[wl_var].replace('', np.nan, inplace=True)
    df.dropna(subset=[wl_var], inplace=True)

    aoi = str(gauge_folder).split('/')[-1]
    gauge_locations = gpd.read_file(gauge_folder / (aoi + '_gauges.shp'))
    gauge_lat = gauge_locations[gauge_locations['name']==gauge_name].geometry.y.item()
    gauge_lon = gauge_locations[gauge_locations['name']==gauge_name].geometry.x.item()
    
    try:
        gauge_Time = pd.to_datetime(df[date_var]).dt.tz_localize('UTC') #- datetime.timedelta(hours=1)
    except:
        gauge_Time = pd.to_datetime(df[date_var]) #-  datetime.timedelta(hours=1)

    ## Make figures
    try:
        nearest_gauge_Time, loc = nearest_gauge(gauge_Time,avg_swot_time)    
        loc = loc[0].item() 
        interval = int(60/((gauge_Time.iloc[loc+1] - gauge_Time.iloc[loc]).seconds/60))
        ax1 = fig.add_subplot(2, 3, 1)
        ax1.plot(gauge_Time.iloc[loc-(interval*24*5):loc+(interval*24*5)],df[wl_var].astype(float).iloc[loc-(interval*24*5):loc+(interval*24*5)],c='k')
        ax1.scatter(gauge_Time.iloc[loc],df[wl_var].astype(float).iloc[loc],s=50,c='r')
        ax1.set_ylim((minwse,maxwse))
        ax1.set_ylabel('%s Water Surface Elevation (m geoid)' %(gauge_name))
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))
        ax1.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    except:
        print('no gauge')
    
    shrink = 0.9
    pad = 0.1
    aspect = 20

    ## Plot SWH
    cax2 =  ax2.pcolormesh(ds.longitude, ds.latitude, 
            ds.swh_karin,cmap='plasma',vmin=0,vmax=2)
    add_features(ax2,cax2,shrink,pad,aspect,gauge_lon,gauge_lat,area,ds.swh_karin.long_name)


    ## Plot SSH Quality flags
    cax3 =  ax3.pcolormesh(ds.longitude, ds.latitude, 
            ds.sig0_karin_2,vmin=-10,vmax=55,cmap='plasma')
    add_features(ax3,cax3,shrink,pad,aspect,gauge_lon,gauge_lat,area,ds.sig0_karin_2.long_name)

    ## Plot Sigma0 
    cax4 =  ax4.pcolormesh(ds.longitude, ds.latitude, 
            ds.ssh_karin_2,cmap='viridis',vmin=minssh,vmax=maxssh)
    add_features(ax4,cax4,shrink,pad,aspect,gauge_lon,gauge_lat,area,ds.ssh_karin_2.long_name)

    ## Plot Xover Height Corrections 
    cax5 =  ax5.pcolormesh(ds.longitude, ds.latitude, 
            ds.height_cor_xover,cmap='PiYG',vmin=-2,vmax=2)
    add_features(ax5,cax5,shrink,pad,aspect,gauge_lon,gauge_lat,area,ds.height_cor_xover.long_name)
    
    ## Plot Water Surface Elevation
    cax6 =  ax6.pcolormesh(ds.longitude, ds.latitude, 
            ds.wse,cmap='viridis',vmin=minwse,vmax=maxwse)
    add_features(ax6,cax6,shrink,pad,aspect,gauge_lon,gauge_lat,area,ds.wse.long_name)
    
    plt.suptitle("SWOT LR L2 Expert 2km: %s" %(date),fontsize=20)
    plt.tight_layout()
    plt.savefig(fig_folder / (o_2000 + '.png'))
    plt.close()
    
    
def plot_l3(fig_folder,o_l3,matching_l3,date,ds,area,mindac,maxdac,minssh,maxssh,minwse,maxwse,avg_swot_time,gauge_folder,gauge_name,wl_var,date_var):
    fig, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2, 3, figsize=(15, 12),subplot_kw=dict(projection=ccrs.PlateCarree()))
    gauge = [os.path.join(dirpath,f)
                for dirpath,dirnames, files in os.walk(gauge_folder)
                for f in fnmatch.filter(files,'*%s.csv' %(gauge_name))][0]
    df = gpd.read_file(gauge)
    df[wl_var].replace('', np.nan, inplace=True)
    df.dropna(subset=[wl_var], inplace=True)
    
    aoi = str(gauge_folder).split('/')[-1]
    gauge_locations = gpd.read_file(gauge_folder / (aoi + '_gauges.shp'))

    gauge_lat = gauge_locations[gauge_locations['name']==gauge_name].geometry.y.item()
    gauge_lon = gauge_locations[gauge_locations['name']==gauge_name].geometry.x.item()
    shrink = 0.9
    pad = 0.1
    aspect = 20
    
    try:
        gauge_Time = pd.to_datetime(df[date_var]).dt.tz_localize('UTC') #- datetime.timedelta(hours=1)
    except:
        gauge_Time = pd.to_datetime(df[date_var]) #-  datetime.timedelta(hours=1)
        
    ## Make figures
    try:
        nearest_gauge_Time, loc = nearest_gauge(gauge_Time,avg_swot_time)    
        loc = loc[0].item() 
        interval = int(60/((gauge_Time.iloc[loc+1] - gauge_Time.iloc[loc]).seconds/60))
        ax1 = fig.add_subplot(2, 3, 1)
        ax1.plot(gauge_Time.iloc[loc-(interval*24*5):loc+(interval*24*5)],df[wl_var].astype(float).iloc[loc-(interval*24*5):loc+(interval*24*5)],c='k')
        ax1.scatter(gauge_Time.iloc[loc],df[wl_var].astype(float).iloc[loc],s=50,c='r')
        ax1.set_ylim((minwse,maxwse))
        ax1.set_ylabel('%s Water Surface Elevation (m geoid)' %(gauge_name))
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))
        ax1.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    except:
        print('no gauge on this date')

    ## Plot SWH
    cax2 =  ax2.pcolormesh(ds.longitude, ds.latitude, 
            ds.ssha_unedited,cmap='plasma',vmin=minwse,vmax=maxwse)
    add_features(ax2,cax2,shrink,pad,aspect,gauge_lon,gauge_lat,area,ds.ssha_unedited.long_name)

    ## Plot SSH Quality flags
    cmap = colors.ListedColormap(['white','black','brown','red','orange','yellow','green','cyan','blue','purple'])
    bounds=[-2.5,2.5,7.5,15,25,40,60,99.5,100.5,101.5,102.5]
    # ticks=[0,5,10,20,30,50,70,100,101,102]
    norm = colors.BoundaryNorm(bounds, cmap.N)
   
    cax3 =  ax3.pcolormesh(ds.longitude, ds.latitude, 
            ds.quality_flag,cmap=cmap,norm=norm)
    add_features(ax3,cax3,shrink,pad,aspect,gauge_lon,gauge_lat,area,ds.quality_flag.long_name)

    ## Plot Sigma0 
    cax4 =  ax4.pcolormesh(ds.longitude, ds.latitude, 
            ds.sigma0,cmap='viridis',vmin=-10,vmax=55)
    add_features(ax4,cax4,shrink,pad,aspect,gauge_lon,gauge_lat,area,ds.sigma0.long_name)

    ## Plot Water Surface Elevation (modified with +MSS -Geoid)
    cax5 =  ax5.pcolormesh(ds.longitude, ds.latitude, 
            ds.wse_mss,cmap='viridis',vmin=minwse,vmax=maxwse)#1
    add_features(ax5,cax5,shrink,pad,aspect,gauge_lon,gauge_lat,area,ds.wse_mss.long_name)

    ## Plot Water Surface Elevation (original equation)
    cax6 =  ax6.pcolormesh(ds.longitude, ds.latitude, 
            ds.wse,cmap='viridis',vmin=minwse,vmax=maxwse)
    add_features(ax6,cax6,shrink,pad,aspect,gauge_lon,gauge_lat,area,ds.wse.long_name)

    plt.suptitle("%s: %s" %('_'.join(o_l3.split('_')[:7]),date),fontsize=20)
    plt.tight_layout()
    plt.savefig(fig_folder /(o_l3 + '.png'))
    plt.close()
    
    
    
def plot_raster(fig_folder, o_100, aoi,variables_to_plot,variables_to_save, full_df,variable_range,area,l2_name,resolution,time ):
     from matplotlib.font_manager import FontProperties
     print('\t\t\tplotting')
     aspect = (area[2]-area[0])/(area[3]-area[1])
     cols = 4
     rows = np.ceil(variables_to_plot/4).astype(int)
     fig2, ax = plt.subplots(rows,cols,figsize=(7*aspect,9),subplot_kw=dict(projection=ccrs.PlateCarree()))#subplot_kw={'projection': ccrs.UTM(zone=19,southern_hemisphere=True)})
     ax = ax.ravel()
     v = 0
     
     for variable in variables_to_plot:
         idx = variables_to_save.index(variable)
         im =  ax[v].scatter(full_df['x'],full_df['y'],c=full_df[variable], s=.02,cmap='viridis',vmin=variable_range[idx][0],vmax=variable_range[idx][1])
         cbar = plt.colorbar(im, ax=ax[v],
             shrink=0.6, orientation='horizontal',
             pad=0.1, aspect=20)
         cbar.set_label(label=variable,size=10,weight='bold')
         ax[v].set_ylim([area[1],area[3]])
         ax[v].set_xlim([area[0],area[2]])
         ax[v].add_feature(cfeature.COASTLINE)
         ax[v].add_feature(cfeature.RIVERS)
         # ax[v].set_extent([area[0],area[2],area[1],area[3]])
         v=v+1
     fp = FontProperties(family="sans-serif", size=25, weight="bold")
     plt.suptitle((l2_name+'_'+resolution) + ' @ ' + time).set_fontproperties(fp)
     plt.tight_layout()#rect=[0, 0.03, 1, 0.95])
     plt.savefig(fig_folder / ('%s_%s.png' %(aoi, o_100)),dpi=300)
     plt.close()
     


def plot1variable_LR_L2orL3_exp(ds,title, vmin,vmax,variable,cmap,ax,area):
    shrink = 0.9
    pad = 0.1
    aspect = 20
    cax =  ax.pcolormesh(ds.longitude, ds.latitude, ds[variable],cmap=cmap,vmin=vmin,vmax=vmax)
    add_colorbar(ax,cax,shrink,pad,aspect,ds[variable].long_name,area)
    ax.set_title("%s: %s " %(title, ds[variable].long_name),fontsize=10)
    print('\tPlotting %s %s' %(title, variable))


def plot1variable_LR_L2orL3_uns(ds_right,title,ds_left,vmin,vmax,variable,cmap,ax,area):
    shrink = 0.9
    pad = 0.1
    aspect = 20
    try:
        cax =  ax.pcolor(ds_right.longitude, ds_right.latitude, 
            ds_right[variable],cmap=cmap,vmin=vmin,vmax=vmax)
    except:''
    try:
        cax =  ax.pcolor(ds_left.longitude, ds_left.latitude, 
            ds_left[variable],cmap=cmap,vmin=vmin,vmax=vmax)
    except:''
    add_colorbar(ax,cax,shrink,pad,aspect,ds_right[variable].long_name,area)
    ax.set_title("%s: %s " %(title, ds_right[variable].long_name),fontsize=10)
    print('\tPlotting %s %s' %(title, variable))


from pyproj import Transformer

def safe_get_crs(ds):
    """Safely resolve the CRS from rioxarray, spatial_ref, or SWOT granule metadata."""
    # 1. Try rioxarray CRS
    try:
        if hasattr(ds, 'rio') and ds.rio.crs is not None and str(ds.rio.crs).strip() != "":
            return ds.rio.crs
    except Exception:
        pass

    # 2. Check for spatial_ref coordinate or variable
    if 'spatial_ref' in ds.coords or 'spatial_ref' in ds.data_vars:
        sref = ds['spatial_ref']
        if 'crs_wkt' in sref.attrs:
            return sref.attrs['crs_wkt']
        if 'spatial_ref' in sref.attrs:
            return sref.attrs['spatial_ref']

    # 3. Reconstruct EPSG from SWOT metadata attributes (utm_zone_num & mgrs_latitude_band)
    utm_zone = ds.attrs.get('utm_zone_num')
    if utm_zone is not None:
        try:
            zone = int(utm_zone)
            lat_band = str(ds.attrs.get('mgrs_latitude_band', 'N')).upper()
            is_north = lat_band >= 'N'
            epsg = 32600 + zone if is_north else 32700 + zone
            return f"EPSG:{epsg}"
        except (ValueError, TypeError):
            pass

    # 4. Default fallback for Arctic Alaska / Prudhoe Bay (Zone 6N)
    return "EPSG:32606"


def plot1variable_HR_Raster(ds, title, vmin, vmax, variable, cmap, ax, area):
    shrink = 0.6
    pad = 0.1 
    aspect = 20

    # Extract 2D lon/lat arrays
    lon = ds['longitude'].values if 'longitude' in ds else None
    lat = ds['latitude'].values if 'latitude' in ds else None

    # If lon/lat contain NaNs from masking/clipping, reconstruct a clean grid from x & y
    if lon is None or lat is None or np.isnan(lon).any() or np.isnan(lat).any():
        crs = safe_get_crs(ds)
        transformer = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
        x_mesh, y_mesh = np.meshgrid(ds.x.values, ds.y.values)
        lon, lat = transformer.transform(x_mesh, y_mesh)

    # Render continuous 2D raster with pcolormesh
    cax = ax.pcolormesh(
        lon,
        lat,
        ds[variable].values,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        shading='auto',
        transform=ccrs.PlateCarree()
    )
    
    add_colorbar(ax, cax, shrink, pad, aspect, variable, area)
    ax.set_title("%s: %s " % (title, variable), fontsize=10)
    print('\tPlotting %s %s' % (title, variable))
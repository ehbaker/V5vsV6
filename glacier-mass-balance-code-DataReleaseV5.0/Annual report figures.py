# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 19:14:42 2019

@author: cmcneil
"""

import numpy as np,pandas as pd, geopandas as gpd,matplotlib.pyplot as plt, datetime as dt
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
from shapely.geometry import Polygon
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

def append_mb_to_rgi(rgi,rgi_glacier_names,bmg_glacier_names,bmg_dir,year):
    "function to append USGS Benchmark glacier data to rgi shapefile"
    'rgi = geodataframe of rgi shapefile'
    'rgi_glacier_names = names of glaciers, as they appear in the rgi, that you want selected'
    'bmg_glacier_names = names of glaciers, as they appear in the BMG directory, you want selected'
    'bmg_dir = path to USGS bmg dataset'
    'year = year of USGS bmg dataset desired'
    
    
    bmgs=gpd.GeoDataFrame()
    for counter,glacier in enumerate(rgi_glacier_names):
        glacier_data=[]
        glacier_data=rgi[rgi.Name==glacier]
    
        fn='Output_'+bmg_glacier_names[counter]+'_Glacier_Wide_solutions_calibrated.csv'
        data_path=os.path.join(bmg_dir,bmg_glacier_names[counter],'Output',fn)
        mb_data=pd.read_csv(data_path)
        data=[]
        data=mb_data[mb_data.loc[:,'Year']==year]
        if data.empty:
            glacier_data.insert(1,"Bw",np.nan)
            glacier_data.insert(1,"Bs",np.nan)
            glacier_data.insert(1,"Ba",np.nan)
        else:
            glacier_data.insert(1,"Bw",np.float64(data['Bw_mwe']))
            glacier_data.insert(1,"Bs",np.float64(data['Bs_mwe']))
            glacier_data.insert(1,"Ba",np.float64(data['Ba_mwe']))
        bmgs=bmgs.append(glacier_data)
    return bmgs
    
def regional_mb_map(bmgs, glacier_shp, countries_shp,ax=None, field=None, title=None,ref_map=False, legend=False):
    ax=ax or plt.gca()
    buffer=500000
    buffer_2=1000000
    minx, miny, maxx, maxy = bmgs.geometry.total_bounds
    ax=plt.subplot(1,1,1)
    countries_shp.plot(ax=ax, facecolor='xkcd:sandstone', edgecolor='black')
    glacier_shp.plot(ax=ax, facecolor='white', edgecolor='white')
    norm = mpl.colors.Normalize(vmin=-5, vmax=5)
    if field is not None:
        bmgs.plot(ax=ax, facecolor=cm.ScalarMappable(norm=norm, cmap=cm.jet_r).to_rgba(np.float64(bmgs[field])), edgecolor='None',markersize=np.exp(np.abs(np.float64(bmgs[field])))*100,legend = True)
    ax.set_xlim(((minx-buffer),(maxx+buffer)))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylim(((miny-buffer),(maxy+buffer)))
    ax.set_facecolor(('xkcd:dusky blue'))
    if title is not None:
        ax.set_title(title, fontsize=18, fontweight='bold')
    if legend is True:
        for a in [5.0, 2.5 , 0.0, -2.5,  -5.0]:
            ax.scatter([], [], c=cm.ScalarMappable(norm=norm, cmap=cm.jet_r).to_rgba(np.float64(a)),label=str(a) + ' m w.e.')
    
        legend=ax.legend(fontsize=24)
        ax.legend(loc='upper right')
    if ref_map is True:
        axins = inset_axes(ax, width=3, height=3,loc=2)
        countries_shp.plot(ax=axins, facecolor='gray', edgecolor='black',alpha=.9)
        axins.set_xlim(((minx-buffer_2),(maxx+buffer_2)))
        axins.set_ylim(((miny-buffer_2),(maxy+buffer_2)))
        bbox=Polygon([(minx-buffer,maxy+buffer),(maxx+buffer,maxy+buffer),(maxx+buffer,miny-buffer),(minx-buffer,miny-buffer)])
        axins.plot(*bbox.exterior.xy, color='black')
        axins.set_xticks([])
        axins.set_yticks([])
    return 

def animate_mb(frame):
    year=1953+frame
    bmgs=gpd.GeoDataFrame()
    bmgs=append_mb_to_rgi(glacier_shp,glaciers,glacier_abbrevs,bmg_dir,year)  
    minx, miny, maxx, maxy = bmgs.geometry.total_bounds
    bmgs['Points']=bmgs['geometry'].centroid
    bmgs=bmgs.set_geometry('Points')
    regional_mb_map(bmgs, glacier_shp, countries_shp, ax=ax, field='Ba', title=str(int(year)),ref_map=True,legend=False)
    return

now=dt.datetime.now()

bmg_dir=r'C:\Users\cmcneil\Desktop\Code\BMG_GUI\glacier-mass-balance-code\data'

countries_fn=r'G:\GIS\World\countries\ne_50m_admin_0_countries.shp'
countries_shp=gpd.read_file(countries_fn)

rgi_fn=r'C:\Users\cmcneil\Desktop\Code\gnp_regional\data\rgi60\regions\rgi60_merge.shp'
glacier_shp=gpd.read_file(rgi_fn)

glaciers=["Gulkana Glacier", "Wolverine Glacier","Lemon Creek Glacier"]
glacier_abbrevs=["Gulkana", "Wolverine","LemonCreek"]

year=now.year
bmgs=append_mb_to_rgi(glacier_shp,glaciers,glacier_abbrevs,bmg_dir,year)
   

minx, miny, maxx, maxy = bmgs.geometry.total_bounds
bmgs['Points']=bmgs['geometry'].centroid
bmgs=bmgs.set_geometry('Points')
buffer=3
glacier_shp = glacier_shp.cx[(minx-buffer):(maxx+buffer), (miny-buffer):(maxy+buffer)]
countries_shp=countries_shp.cx[(minx-buffer):(maxx+buffer), (miny-buffer):(maxy+buffer)]

glacier_shp=glacier_shp.to_crs({'init': 'epsg:3338'})
countries_shp=countries_shp.to_crs({'init': 'epsg:3338'})   
bmgs=bmgs.to_crs({'init': 'epsg:3338'}) 

#generate figure containing subplots for Winter, Summer, and Annual mass balance solutions
fig, (ax1,ax2,ax3) = plt.subplots(3, figsize=(15, 15))
regional_mb_map(bmgs, glacier_shp, countries_shp, ax=ax1, field='Bw', title='2019 Bw',ref_map=False)
regional_mb_map(bmgs, glacier_shp, countries_shp, ax=ax2, field='Bs', title='2019 Bs',ref_map=False)
regional_mb_map(bmgs, glacier_shp, countries_shp, ax=ax3, field='Ba', title='2019 Ba',ref_map=True)
fig.savefig(bmg_dir+'/'+str(now.year)+'result map.jpeg', dpi=200)

#generate figure containing single season of mass balance solutions
regional_mb_map(bmgs, glacier_shp, countries_shp,field='Ba', title='2019 Ba',ref_map=True)

from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
#generate a movie of mass balance solutions through time
#Writer = animation.writers['ffmpeg'] #set animation writer
#writer = Writer(fps=2, metadata=dict(artist='Me'), bitrate=1800)
#year=1953
#bmgs=gpd.GeoDataFrame()
#bmgs=append_mb_to_rgi(glacier_shp,glaciers,glacier_abbrevs,bmg_dir,year)

#minx, miny, maxx, maxy = bmgs.geometry.total_bounds
#bmgs['Points']=bmgs['geometry'].centroid
#bmgs=bmgs.set_geometry('Points')
#fig, ax1 = plt.subplots(figsize=(15, 15)) 
#regional_mb_map(bmgs, glacier_shp, countries_shp, ax=ax1, field='Ba', title=str(int(year)),ref_map=True, legend=True)
#ani = animation.FuncAnimation(fig, animate_mb, frames=np.linspace(0,66, 67), repeat=True) #this calls the 
#ani.save(bmg_dir+"/glacierchangemove.mp4", writer=writer)


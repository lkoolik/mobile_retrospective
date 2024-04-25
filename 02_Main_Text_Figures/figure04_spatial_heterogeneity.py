#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PM2.5 exposure disparities persist despite strict vehicle emissions controls in California
Koolik et al. (2024)

Figure 4

"""

# Import libraries
import pandas as pd # v1.4.2
import geopandas as gpd # v0.10.2
import matplotlib.pyplot as plt # v3.5.1
from matplotlib import gridspec 
import seaborn as sns # v0.11.2
sns.set_theme(context=None, font_scale=1, rc=None, style='ticks')
from cmcrameri import cm # v1.7
import numpy as np # v1.22.3
from os import path

#%% Define Critical User Paths Up Front
# File path for the summarized PWM and Relative Disparity by Race-Ethnicity
# This file can be generated by running b_pre_process_by_policy.py
policy_pwm_fp = '.../disparities_by_policy.csv' # add path

# Paths to AB617 geography, California boundary, and California counties shapefiles
ab617_fp = '.../ab617_communities.shp' # add path
ca_fp = '.../ca_border.feather' # add path (can be downloaded from the ECHO-AIR materials)
counties_fp = '.../counties.feather' # add path (can be downloaded from the ECHO-AIR materials)

# Folder containing exposure shapefiles (e.g., all_cy2000_exposure_concentrations.shp)
main_data_path = '...' # add path

#%% Define Global Variables
# Vehicle Color Palette
vt_palette = cm.acton
vt_colors = {'OTH':(0.8,0.8,0.8), 
             'HDV':vt_palette(0.0), 
             'MDV':vt_palette(0.375), 
             'LDV':vt_palette(0.75)}
vt_color_list = [vt_colors['OTH'], vt_colors['HDV'], vt_colors['MDV'],
                 vt_colors['LDV']] # list format will be usefull too

# Define colors for each of the regions
region_palette = cm.imola
region_colors = [region_palette(0.25), #LA
                 region_palette(0.5), #CV
                 region_palette(0.75)] #SF

#%% Load the data and filter for just Hispanic Californians
# Load the data
policy_pwm = pd.read_csv(policy_pwm_fp)

# Filter for only 2010
policy_pwm = policy_pwm[policy_pwm['YEAR']==2010].copy()

# Grab the total statewide PWM data
statewide = policy_pwm[['YEAR','SOURCE','TOTAL_PWM']].drop_duplicates().copy().reset_index(drop=True)
statewide = statewide[['SOURCE','TOTAL_PWM']].set_index('SOURCE').transpose().reset_index()
statewide['CODE'], statewide['GROUP'] = 'TOTAL', 'TOTAL'
statewide = statewide[['CODE','GROUP','ALL','HDV','LDV','MDV','OTH']].copy()

# Remove unnecessary columns and pivot to get sources as columns
policy_pwm = policy_pwm[['YEAR','SOURCE','GROUP','CODE','PWM']].copy().reset_index(drop=True)
policy_pwm = policy_pwm.pivot(columns=['SOURCE'], index=['CODE','GROUP'], values='PWM').reset_index()

# Add statewide back in
policy_pwm = pd.concat([policy_pwm, statewide], ignore_index=True).reset_index(drop=True)

# Calculate the fraction from each source
for vt in vt_colors.keys():
    policy_pwm[vt] = policy_pwm[vt] / policy_pwm['ALL']

#%% Split the data up based on how it will be plotted for simplicity
# Grab just the AB617 communities
ab617_communities = policy_pwm[~policy_pwm['CODE'].isin(['TOTAL','AB617','SB535','CV','SF','LA'])].copy()

# Grab the summary data
summary_data = policy_pwm[policy_pwm['CODE'].isin(['TOTAL','AB617','SB535','CV','SF','LA'])].copy()

#%% Load the geography and get everything on the same CRS
# Load the AB617 data
ab617_geo = gpd.read_file(ab617_fp)

# Load the CA data
ca = gpd.read_feather(ca_fp)

# Load the county data
counties = gpd.read_feather(counties_fp)

# Load one concentration dataset to use that CRS
crs_to_use = gpd.read_file(path.join(main_data_path,
                                     'all_cy2000_exposure_concentrations.shp')).crs

# Project both AB617 and CA
ab617_geo = ab617_geo.to_crs(crs_to_use)
ca = ca.to_crs(crs_to_use)
counties = counties.to_crs(crs_to_use)

#%% Estimate locations for each of the AB617 communities to be plotted
# Estimate the centroids
ab617_geo['CENTROID'] = ab617_geo.centroid

# Add these to the AB617 communities data
ab617_plot = pd.merge(ab617_communities, ab617_geo[['code','CENTROID']], 
                      left_on='CODE', right_on='code')
ab617_plot = gpd.GeoDataFrame(ab617_plot, geometry=ab617_plot['CENTROID'],
                              crs=crs_to_use)

# We add offsets manually to help prevent overlap
ab617_offset_dict = {'East Los Angeles, Boyle Heights, West Commerce':[50000,75000],
                     'South Central Fresno':[0,0], 
                     'San Bernardino, Muscoy':[75000,0], 
                     'West Oakland':[-75000,10000],
                     'Portside Environmental Justice Neighborhoods':[-50000,25000],
                     'El Centro, Heber, Calexico':[-50000,-25000], 
                     'South East Los Angeles':[-100000,-75000],
                     'Wilmington, Carson, West Long Beach':[50000,-50000],
                     'South Sacramento - Florin':[50000,50000],
                     'South Los Angeles':[-100000,75000], 
                     'Shafter':[-40000,40000], 
                     'Richmond - San Pablo':[-50000,50000],
                     'Arvin, Lamont':[40000,40000], 
                     'Eastern Coachella Valley':[0,0], 
                     'Stockton':[50000,-50000],
                     'Bayview Hunters Point/Southeast San Francisco':[-75000,-50000], 
                     'East Oakland':[50000,-50000], 
                     'International Border Community':[-50000,-25000], 
                     'Northern Imperial Phase 1':[25000,25000]} 

#%% The other datapoints are provided locations manually
# Make a copy of the summary data dataframe
summary_plot = summary_data.copy()

# Start by providing all of them with the same longitude
summary_plot['LON'] = -2500000

# The latitude will be input manually via dictionary
summary_lats = {'TOTAL': 475000.,
                'AB617': -590000.,
                'SB535': -590000.,
                'CV': -125000.,
                'SF': 0.,
                'LA': -250000.}

# Add latitude to the summary_plot dataframe
summary_plot['LAT'] = summary_plot['CODE'].map(summary_lats)

# Add spacing in between AB617 and SB535
summary_plot.loc[summary_plot['CODE']=='AB617','LON'] = -2400000
summary_plot.loc[summary_plot['CODE']=='SB535','LON'] = -2600000

# Convert into a geodataframe
summary_plot = gpd.GeoDataFrame(summary_plot, 
                                geometry=gpd.points_from_xy(summary_plot['LON'], 
                                                            summary_plot['LAT']),
                                crs=crs_to_use)

#%% Map the relevant regions for plotting
# Create a dictionary to map counties to regions
county_mapper = {'ALAMEDA':'SF', 'ALPINE':'NA', 'AMADOR':'NA', 'BUTTE':'CV', 
                 'CALAVERAS':'NA', 'COLUSA':'CV', 'CONTRA COSTA':'SF', 'DEL NORTE':'NA', 
                 'EL DORADO':'NA', 'FRESNO':'CV', 'GLENN':'CV', 'HUMBOLDT':'NA', 'IMPERIAL':'NA', 
                 'INYO':'NA', 'KERN':'CV', 'KINGS':'CV', 'LAKE':'NA', 'LASSEN':'NA', 'LOS ANGELES':'LA', 
                 'MADERA':'CV', 'MARIN':'SF', 'MARIPOSA':'NA', 'MENDOCINO':'NA', 'MERCED':'CV', 
                 'MODOC':'NA', 'MONO':'NA', 'MONTEREY':'NA', 'NAPA':'SF', 'NEVADA':'NA', 'ORANGE':'LA',
                 'PLACER':'NA', 'PLUMAS':'NA', 'RIVERSIDE':'LA', 'SACRAMENTO':'CV', 'SAN BENITO':'NA',
                 'SAN BERNARDINO':'LA', 'SAN DIEGO':'NA', 'SAN FRANCISCO':'SF', 'SAN JOAQUIN':'CV',
                 'SAN LUIS OBISPO':'NA', 'SAN MATEO':'SF', 'SANTA BARBARA':'NA', 'SANTA CLARA':'SF',
                 'SANTA CRUZ':'NA', 'SHASTA':'CV', 'SIERRA':'NA', 'SISKIYOU':'NA', 
                 'SOLANO':'SF', 
                 'SONOMA':'SF','STANISLAUS':'CV', 'SUTTER':'CV', 'TEHAMA':'CV', 'TRINITY':'NA', 
                 'TULARE':'CV', 'TUOLUMNE':'NA', 'VENTURA':'LA', 'YOLO':'CV', 'YUBA':'CV'}

# Map this to the counties dataframe
counties['REGION'] = counties['NAME'].map(county_mapper)


#%% Define functions to help plot
# First, make a function that gets the values for the pie charts
def get_pie_vals(geo, data_plot, offset_dict):
    ''' Function for estimating the necessary input data for drawing the pie 
        chart map markers.
        INPUTS:
            - geo = the name of the geography
            - data_plot = dataframe containing the data to plot
            - offset_dict = the dictionary with the offsets to prevent point
              overlap'''
    
    # Clip data to individual geography
    data = data_plot[data_plot['GROUP']==geo].copy().reset_index(drop=True)
     
    # If relevant, grab the offset
    if geo in offset_dict.keys():
        offsets = (offset_dict[geo][0],
                   offset_dict[geo][1])
    else:
        offsets = (0,0)
    
    # Get xs and ys by adding the offset
    xs = data.geometry.x + offsets[0]
    ys = data.geometry.y + offsets[1]
        
    # The size should be correlated to the total PWM
    sizes = data.loc[0, 'ALL']*100
    
    # Put the ratios into a list for the pie chart function
    ratios = [data.loc[0,'OTH'], 
              data.loc[0,'HDV'],
              data.loc[0,'MDV'],
              data.loc[0,'LDV']]
    
    # There's a strange Python error that is causing these to be every so slightly above 1, fix this:
    ratios = np.array(ratios)
    ratios = ratios/ratios.sum()
    ratios = list(ratios)
    
    return xs, ys, ratios, sizes

def draw_pie_marker(geo, data_plot, offset_dict, colors, ax):
    ''' Function for drawing the pie chart marker for each unique geography
        INPUTS:
            - geo = the name of the geography
            - data_plot = dataframe containing the data to plot
            - offset_dict = the dictionary with the offsets to prevent point
              overlap
            - colors = the color palette for vehicle contributions
            - ax = the axis object to plot on '''
    
    # Get values from get_pie_vals
    xs, ys, ratios, sizes = get_pie_vals(geo, data_plot, offset_dict)
    
    # Store the wedges to be drawn
    wedges = []
    
    # Set the previous tracker to zero to initialize
    previous = 0
    
    # Calculate the points of the pie pieces individually
    for color, ratio in zip(colors, ratios): # Loop through each chunk
        
        # Calculate the new value wedge endpoint
        current = 2 * np.pi * ratio + previous
        
        # Get the x and y vectors for this new vector
        x  = [0] + np.sin(np.linspace(previous, current, 20)).tolist() + [0]
        y  = [0] + np.cos(np.linspace(previous, current, 20)).tolist() + [0]
        
        # Combine into a stacked array
        xy = np.column_stack([x, y])
        
        # Update the startpoint to the previous endpoint
        previous = current
        
        # Add to the list of markers
        wedges.append({'marker':xy, 's':np.abs(xy).max()**2*np.array(sizes), 'facecolor':color, 'edgecolor':'black', 'linewidth':0.25})

    # Add each wedge to the map
    for wedge in wedges:
        ax.scatter(xs, ys, **wedge)
        
    return #nothing

def draw_pie_legend(xs,ys, size, ax):
    ''' Function for drawing the legends for the pie chart sizes
        INPUTS:
            - xs = the x-coordinate of the center of the plot
            - ys = the y-coordinate of the center of the plot
            - size = the size that should be plotted as a legend
            - ax = the axis object to plot on 
            - lc = the axis object to plot on 
            
        In order to match the size exactly, this will plot with four equal 
        quadrants '''
        
    # Store the wedges to be drawn
    wedges = []

    # Set the previous tracker to zero to initialize
    previous = 0
    
    # Hardcode the ratios since we are just showing the size
    ratios=[0.25,0.25,0.25,0.25] 
    
    # Hardcode the colorscheme since it's a legend figure
    colors, lc = (['white','white','white','white'], 'black')
    
    # Calculate the points of the pie pieces individually
    for color, ratio in zip(colors, ratios): # Loop through each chunk
        
        # Calculate the new value wedge endpoint
        current = 2 * np.pi * ratio + previous
        
        # Get the x and y vectors for this new vector
        x  = [0] + np.sin(np.linspace(previous, current, 20)).tolist() + [0]
        y  = [0] + np.cos(np.linspace(previous, current, 20)).tolist() + [0]
        
        # Combine into a stacked array
        xy = np.column_stack([x, y])
        
        # Update the startpoint to the previous endpoint
        previous = current
        
        # Add to the list of markers
        wedges.append({'marker':xy, 's':np.abs(xy).max()**2*np.array(size), 'facecolor':color, 'edgecolor':lc, 'linewidth':0.25})

    # Add each wedge to the legend
    for wedge in wedges:
        ax.scatter(xs, ys, **wedge)
        
    # Add a label that says the size
    ax.text(xs+75000, ys, r'{} $\mu$g/m$^3$'.format(size/100), ha='left', va='center', fontsize=6)
    
    return #nothing

def add_pies_to_map(data_plot, ax, offset_dict=ab617_offset_dict, label=False):
    ''' Function for adding pie charts as markers for the map
        INPUTS:
            - data_plot = dataframe containing the data to plot
            - ax = the axis object to plot on
            - offset_dict = the dictionary with the offsets to prevent point
              overlap
            - label = Boolean to indicate if a label shoould be added '''
            
    # Iterate through sub-geographies
    for geo in data_plot['GROUP'].unique():
        draw_pie_marker(geo, data_plot, offset_dict, colors=vt_color_list, ax=ax)
        
        # For the summary ones, add a label
        if label:
            # Update the formatting of the label
            label_text = {'TOTAL':'Statewide\nAverage',
                          'San Francisco Bay Area':'San Francisco\nBay Area',
                          'Central Valley':'Central\nValley',
                          'Los Angeles Area':'Los Angeles\nArea',
                          'AB617':'AB617', 'SB535':'SB535'}[geo]
            
            # Add the text
            ax.text(x=data_plot[data_plot['GROUP']==geo].geometry.x, 
                    y=data_plot[data_plot['GROUP']==geo].geometry.y+40000,
                    s=label_text, fontsize=6, ha='center', va='bottom')
        
    return #nothing

#%% Draw the figure
# Initialize the figure
fig, ax = plt.subplots(figsize=(4,6))

# Add the California border
ca.plot(edgecolor='black', facecolor='none', linewidth=0.5, ax=ax)

# Add the region shading
counties[counties['REGION']=='LA'].plot(facecolor=region_colors[0],
                                        linewidth=0, alpha=0.5, ax=ax)
counties[counties['REGION']=='CV'].plot(facecolor=region_colors[1],
                                        linewidth=0, alpha=0.75, ax=ax)
counties[counties['REGION']=='SF'].plot(facecolor=region_colors[2],
                                        linewidth=0, alpha=0.75, ax=ax)

# Add the scaled pie charts for the summary data and the community-specific data
add_pies_to_map(summary_plot, ax, label=True)
add_pies_to_map(ab617_plot, ax)

# Add the legend pies
draw_pie_legend(-1700000, -175000 ,0.5*100, ax)
draw_pie_legend(-1700000, -110000 ,1*100, ax)
draw_pie_legend(-1700000, 0 ,5*100, ax)

# Format the full figure
ax.set_xlim([-2500000-150000,-1500000])
miny = ca.geometry.total_bounds[1] - np.abs(ca.geometry.total_bounds[1]*0.05)
maxy = ca.geometry.total_bounds[3] + np.abs(ca.geometry.total_bounds[3]*0.05)
ax.set_ylim([miny, maxy])
plt.axis('off')
fig.tight_layout()
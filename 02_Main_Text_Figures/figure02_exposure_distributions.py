#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PM2.5 exposure disparities persist despite strict vehicle emissions controls in California
Koolik et al. (2024)

Figure 2

"""

# Import libraries
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
sns.set_theme(context=None, font_scale=1, rc=None, style='ticks')
from cmcrameri import cm
import numpy as np

#%% Define Critical User Paths Up Front
# File path for the shapefile containing exposure to all vehicles for calendar year 2000
data_2000_fp = '/Users/libbykoolik/Documents/Research/OEHHA Project/Retrospective Analysis/Tool Outputs/out_all_cy2000/shapes/all_cy2000_exposure_concentrations.shp'

# File path for the shapefile containing exposure to all vehicles for calendar year 2019
data_2019_fp = '/Users/libbykoolik/Documents/Research/OEHHA Project/Retrospective Analysis/Tool Outputs/out_all_cy2019/shapes/all_cy2019_exposure_concentrations.shp'

#%% Define a few global variables
# Define racial-ethnic groups
groups = ['ASIAN', 'BLACK', 'HISLA', 'WHITE', 'OTHER']

# Define color palette for racial/ethnic groups
re_palette_func = cm.corkO
re_palette = {'ASIAN':re_palette_func(0.85), #Asian
              'BLACK':re_palette_func(0.5), #Black
              'HISLA':re_palette_func(0.05), #Hispanic/Latino
              'WHITE':re_palette_func(0.25), #White
              'OTHER':(0.8,0.8,0.8)} #Other

#%% Create a set of functions for streamlining code
def estimate_percentiles(data, group):
    ''' Function for estimating the percentile of exposure for the population
        INPUTS:
            - data = dataframe containing PM2.5 concentrations and population counts
            - group = group of interest
            
        OUTPUTS:
            - pctl_df = dataframe of PM2.5 concentration percentiles '''
            
    # Make a copy of the dataframe that is just the group and the exposure
    pctl_df = data[['PM25_UG_M3', group]].copy()
    
    # Sort by concentration 
    pctl_df.sort_values(by='PM25_UG_M3', inplace=True)
    
    # Reset index to clean up
    pctl_df.reset_index(drop=True, inplace=True)
    
    # Estimate the cumulative population
    pctl_df.loc[:, 'CUMULATIVE_POP'] = pctl_df.loc[:, group].cumsum()
    
    # Estimate the total population in that group, then divide the cumulative sum
    # to get the percentile
    total_pop_group = pctl_df[group].sum()
    pctl_df.loc[:, 'PCTL_'+group] = pctl_df['CUMULATIVE_POP']/total_pop_group
    
    return pctl_df

def get_pctl_val(pctl_df, group, pctl):
    ''' Looks up the nearest value to the exact percentile of interest within a group
        INPUTS:
            - pctl_df = dataframe of PM2.5 concentration percentiles
            - group = group of interest
            - pctl = percentile of interest
            
        OUTPUTS:
            - pctl_val = corresponding concentration at pctl '''
    
    # Find the concentration that corresponds to the percentile pctl
    idx = abs(pctl_df['PCTL_'+group]-pctl).idxmin()
    pctl_val = pctl_df.loc[idx, 'PM25_UG_M3']
    
    return pctl_val

def get_pop_counts(data, lower, upper, dec, groups=groups+['TOTAL']):
    ''' Estimates the number of people of each group in a dataset.
        INPUTS:
            - data = dataframe of PM2.5 concentrations and population counts
            - lower = lower bound of the concentration bin
            - upper = upper bound of the concentration bin
            - dec = the decile number 
            
        OUTPUTS:
            - pop_counts = dataframe with the population percent by group in the decile '''
            
    # Copy the dataframe
    data = data.copy().sort_values(by='PM25_UG_M3').reset_index(drop=True)
    
    # Trim the dataframe for just the cutoff
    data = data.loc[(data['PM25_UG_M3']<=upper)&(data['PM25_UG_M3']>lower)]
    
    # Get the total population
    total_population = data['TOTAL'].sum()
    
    # Return the percent of the population for each group in that decile
    pop_counts = pd.DataFrame(data[groups].sum()/total_population).transpose()
    
    # Label the decile
    pop_counts['DECILE'] = str(dec)
    
    return pop_counts

def process_distributions(data_fp, groups=groups):
    ''' Function for running the calculations necessary for each year of data.
        INPUTS:
            - data_fp = the filepath for the exposure concentration shapefile
            
        OUTPUTS:
            - pop_distrirbutions = dataframe of the distribution of population 
              within each exposure decile
            - cutoffs = array of decile bin limits
            - data = dataframe of PM2.5 concentrations and population counts '''
    
    # Import the data
    data = gpd.read_file(data_fp)
    
    # Create OTHER and remove negative population counts caused by rounding
    data['OTHER'] = data['TOTAL'] - data['ASIAN'] - data['BLACK'] - data['HISLA'] - data['WHITE']
    data.loc[data['OTHER']<0, 'OTHER'] = 0

    # Get rid of extra groups
    data = data[['ISRM_ID', 'PM25_UG_M3', 'TOTAL']+groups].copy()
    
    # Get the exposure decile bins
    pctl_df = estimate_percentiles(data, 'TOTAL')
    
    # Get the exposure cut offs for each decile
    deciles = np.linspace(0,1,11)
    cutoffs = [get_pctl_val(pctl_df, 'TOTAL', pctl) for pctl in deciles]
    
    # Calculate the distribution of population within each bin
    pop_distributions = [] # Placeholder list
    for i in range(len(cutoffs)-1): # Loop through each bin
        lower, upper = cutoffs[i], cutoffs[i+1]
        pop_distributions.append(get_pop_counts(data, lower, upper, i+1))

    # Make a dataframe from these
    pop_distributions = pd.concat(pop_distributions)

    # Remove the total
    pop_distributions = pop_distributions[groups+['DECILE']]
    
    return pop_distributions, cutoffs, data

# Run this for each year, saving the population data only once
pop_distribution_2000, cutoffs_2000, pop_data = process_distributions(data_2000_fp)
pop_distribution_2019, cutoffs_2019, _ = process_distributions(data_2019_fp)

#%% Estimate statewide population for comparison
# Estimate total population by simulating a bin that is the complete range of concentrations
total_pop = get_pop_counts(pop_data, 0, pop_data['PM25_UG_M3'].max(), 'TOTAL')

# Trim for simplicity
total_pop = total_pop[groups+['DECILE']]

#%% Automate the drawing of the panel
def draw_distributions(pop_distributions, cutoffs, ax, year, re_palette=re_palette, deciles=np.linspace(0,1,11)):
    ''' Helper function for plotting distribution of population by exposure
        decile.
        INPUTS:
            - pop_distrirbutions = dataframe of the distribution of population 
              within each exposure decile
            - cutoffs = array of decile bin limits
            - ax = axis object to plot on
            - year = year of data  '''
        
    
    # Plot the horizontal bar plot
    pop_distributions.plot(kind='barh', stacked=True, align='edge', width=1.0, 
                           legend=False, color=re_palette, ax=ax)
    
    # Format the X-Axis
    ax.set_xlabel('Racial-Ethnic Composition')
    ax.set_xlim([0,1])
    ax.set_xticks(ticks=[0, 0.25, 0.5, 0.75, 1], labels=['0%', '25%','50%','75%', '100%'])
    
    # Format the Y-Axis
    ax.set_ylim([0, len(deciles)-1])
    cutoff_labels = [str(np.round(c,1)) for c in cutoffs]
    ax.set_yticks(ticks=range(len(deciles-1)), labels=cutoff_labels)
    ax.set_ylabel('{} On-Road Mobile\n'.format(year)+r'PM$_{2.5}$ Deciles ($\mu$g/m$^3$)')

    return #nothing

#%% Plot this
# Initialize the figure
fig, ax = plt.subplots(figsize=(4, 6))

# Create three uneven panels for the plot
gs = gridspec.GridSpec(3, 1, height_ratios=[10, 10, 1]) 
ax0 = plt.subplot(gs[0]) # a) Distributions in 2000
ax1 = plt.subplot(gs[1]) # b) Distributions in 2019
ax2 = plt.subplot(gs[2]) # c) Statewide population distribution

# Plot the horizontal bar charts for 2000 and 2019
draw_distributions(pop_distribution_2000, cutoffs_2000, ax0, 2000)
draw_distributions(pop_distribution_2019, cutoffs_2019, ax1, 2019)

# Plot the statewide population distribution and format
total_pop.plot(kind='barh', stacked=True, align='edge', width=1.0, legend=False, 
               color=re_palette, ax=ax2)
ax2.set_xlabel('Racial-Ethnic Composition')
ax2.set_xlim([0,1])
ax2.set_xticks(ticks=[0, 0.25, 0.5, 0.75, 1], labels=['0%', '25%','50%','75%', '100%'])
ax2.set_ylim([0, 1])
ax2.set_yticks(ticks=[])
ax2.set_ylabel('Statewide\nPopulation', ha='right', va='center', rotation=0)

# Add a thicker frame for each panels
for axis in ['top','bottom','left','right']:
    ax0.spines[axis].set_linewidth(2)
    ax1.spines[axis].set_linewidth(2)
    ax2.spines[axis].set_linewidth(2)

# Final clean up
fig.tight_layout()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PM2.5 exposure disparities persist despite strict vehicle emissions controls in California
Koolik et al. (2024)

Figure 3

"""

# Import libraries
import pandas as pd # v1.4.2
import matplotlib.pyplot as plt # v3.5.1
from matplotlib import gridspec 
import seaborn as sns # v0.11.2
sns.set_theme(context=None, font_scale=1, rc=None, style='ticks')
from cmcrameri import cm # v1.7
import numpy as np # v1.22.3
from os import path

#%% Define Critical User Paths Up Front
# File path for the summarized PWM and Relative Disparity by Race-Ethnicity
# This file can be generated by running a_pre_process_by_race_ethnicity.py
re_pwm_fp = '/Users/libbykoolik/Documents/Research/OEHHA Project/Retrospective Analysis/Final_Scripts/processed_data/disparities_by_race_ethnicity.csv'

#%% Define Global Variables
# Vehicle Color Palette
vt_palette = cm.acton
vt_colors = {'OTH':(0.8,0.8,0.8), 
             'HDV':vt_palette(0.0), 
             'MDV':vt_palette(0.375), 
             'LDV':vt_palette(0.75)}
vt_color_list = [vt_colors['OTH'], vt_colors['HDV'], vt_colors['MDV'],
                 vt_colors['LDV']] # list format will be usefull too

# Define the years
years = range(2000,2020)

#%% Load the data and filter for just Hispanic Californians
# Load the data
re_pwm = pd.read_csv(re_pwm_fp)

# Filter for Hispanic Californians
hispanic_pwm = re_pwm[re_pwm['GROUP']=='HISPANIC'].copy()

# Remove unnecessary columns
hispanic_pwm = hispanic_pwm[['YEAR','SOURCE','ABSOLUTE_DISP','RELATIVE_DISP']].copy().reset_index(drop=True)

#%% Create a function for automating plotting the first two panels slightly
def make_stackplot(years, pwm_data, ax, normalize=False, vt_color_list=vt_color_list):
    ''' Function for plotting absolute disparity by vehicle type scatcked
        INPUTS:
            - years = array of relevant years of study
            - pwm_data = dataframe containing absolute disparity by group
              by source
            - ax = axis object to plot on
            - normalize = Boolean of whether or not to normalize by total 
              absolute disparity
            - vt_color_list = list of vehicle type colors '''
            
    # Create the vectors for each vehicle type
    oth = np.array(pwm_data.loc[pwm_data['SOURCE']=='OTH']['ABSOLUTE_DISP'])
    mdv = np.array(pwm_data.loc[pwm_data['SOURCE']=='MDV']['ABSOLUTE_DISP'])
    ldv = np.array(pwm_data.loc[pwm_data['SOURCE']=='LDV']['ABSOLUTE_DISP'])
    hdv = np.array(pwm_data.loc[pwm_data['SOURCE']=='HDV']['ABSOLUTE_DISP'])
    
    # Normalize, if relevant by dividing by the total disparity 
    if normalize:
        # Estimate the fractional contribution to disparity for each group
        othf = oth / (oth + mdv + ldv + hdv)
        mdvf = mdv / (oth + mdv + ldv + hdv)
        ldvf = ldv / (oth + mdv + ldv + hdv)
        hdvf = hdv / (oth + mdv + ldv + hdv)
    
        # Plot the data
        ax.stackplot(years, np.array([othf, hdvf, mdvf, ldvf]), 
                     colors=vt_color_list, linewidth=0)
    
    else:        
        # Plot the data
        ax.stackplot(years, np.array([oth, hdv, mdv, ldv]), 
                     colors=vt_color_list, linewidth=0)
    
    # Format the x-axis
    ax.set_xlabel('Year')
    ax.set_xticks(ticks=years, labels=['2000','','','','',
                                       '2005','','','','',
                                       '2010','','','','',
                                       '2015','','','',''])
    ax.set_xlim([2000,2019])
    
    # Format the y-axis
    if normalize:
        ax.set_ylim([0,1])
        ax.set_yticks(ticks=[0.,0.2, 0.4, 0.6, 0.8, 1.0], labels=['0','20','40','60','80','100'])
        ax.set_ylabel('Contribution to On-Road\nMobile Source Disparity (%)')
    else:
        ax.set_ylim([0,0.4])
        ax.set_yticks(ticks=np.arange(0,0.41,0.1))
        ax.set_ylabel('Absolute '+r'PM$_{2.5}$'+'\nDisparity '+r'($\mu$g/m$^3$)')
    
    return 

#%% Create the Figure
fig, ax = plt.subplots(figsize=(3.5,8))

# Create three uneven panels for the plot
gs = gridspec.GridSpec(3, 1, height_ratios=[4, 8, 8]) 
ax0 = plt.subplot(gs[0]) # a) Absolute disparity by source
ax1 = plt.subplot(gs[1]) # b) Fractional disparity by source
ax2 = plt.subplot(gs[2]) # c) Relative disparity by source

# # # # # # # # # # # # 
# Panels (a) and (b)  #
# # # # # # # # # # # #
# Plot panels A and B using the function defined above
make_stackplot(years, hispanic_pwm, ax0, normalize=False)
make_stackplot(years, hispanic_pwm, ax1, normalize=True)

# # # # # # # # # # # #
# Panel (c)           #
# # # # # # # # # # # #
# Plot the relative disparity as lines
sns.lineplot(data=hispanic_pwm[hispanic_pwm['SOURCE']!='ALL'], x='YEAR', 
              y='RELATIVE_DISP', hue='SOURCE', palette=vt_colors, legend=False,
              marker='o', ax=ax2)

# Format the x-axis
ax2.set_xlim([2000,2019])
ax2.set_xticks(ticks=years, labels=['2000','','','','',
                                    '2005','','','','',
                                    '2010','','','','',
                                    '2015','','','',''])
ax2.set_xlabel('Year')

# Format the y-axis
ax2.set_yticks(ax2.get_yticks(), 
               labels=['{:.0f}'.format(np.round(y*100,0)) for y in ax2.get_yticks()])
ax2.set_ylabel('Relative Disparity from\nOn-Road Mobile Source (%)')

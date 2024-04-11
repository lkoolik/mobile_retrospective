#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: libbykoolik
"""

# Import libraries
import pandas as pd
import matplotlib.pyplot as plt
from cmcrameri import cm
import seaborn as sns
import numpy as np

#%% Define Critical User Paths Up Front
# Where to output files:
out_dir = '/Users/libbykoolik/Documents/Research/OEHHA Project/Retrospective Analysis/Figures/final_manuscript_figures/figure_01'

# File path for the summarized PWM and Relative Disparity by Race-Ethnicity
re_pwm_fp = '/Users/libbykoolik/Documents/Research/OEHHA Project/Retrospective Analysis/Figures/emfac2021/all_statewide/disparities_all.csv'

# File path for the summarized PWM and Relative Disparity by Policy 
policy_pwm_fp = '/Users/libbykoolik/Documents/Research/OEHHA Project/Retrospective Analysis/Figures/emfac2021/all_statewide_by_policy/all/all_CES_PWM_Exposure_Over_Time.csv'

#%% Set up defaults for plot formatting
# Set seaborn plot style
sns.set_theme(context='notebook', style="ticks")

# Get years variable
years = np.array(range(2000,2020))

# Define color palette for racial/ethnic groups
re_palette_func = cm.corkO
palette = {'ASIAN':re_palette_func(0.85), #Asian
           'BLACK':re_palette_func(0.5), #Black
           'HISPANIC':re_palette_func(0.05), #Hispanic/Latino
           'WHITE':re_palette_func(0.25), #White
           'OTHER':(0.8,0.8,0.8), #Other
           'AB617':'lightsalmon',
           'SB535':'indianred'}

#%% Load the PWM and disparity data for race-ethnicity and by policy and combine
# Load Data by Race-Ethnicity
re_pwms = pd.read_csv(re_pwm_fp)

# Load Data by Policy Designation
policy_pwms = pd.read_csv(policy_pwm_fp)

# Cut to just relevant columns
policy_pwms = policy_pwms[['YEAR','Total','AB617', 'SB535']].copy()
re_pwms = re_pwms[['YEAR','GROUP','PWM']].copy()

# Update the data by race/ethnicity to match the formatting of the policy results
re_pwms = re_pwms.pivot_table(index='YEAR', values='PWM', columns='GROUP').reset_index()
re_pwms.columns = re_pwms.columns.str.upper()

# Update the data by policy for consistent formatting
policy_pwms.columns = policy_pwms.columns.str.upper()

# Combine both datasets
pwms = pd.merge(re_pwms, policy_pwms, on='YEAR')

#%% Calculate absolute and relative disparity per group per year
# Define the relevant groups
groups = ['ASIAN', 'BLACK', 'HISPANIC', 'OTHER', 'WHITE','AB617', 'SB535']

# Loop through to calculate disparity
for group in groups:
    pwms[group+'_ABS_DISP'] = pwms[group] - pwms['TOTAL']
    pwms[group+'_REL_DISP'] = pwms[group+'_ABS_DISP'] / pwms['TOTAL'] * 100.
    
#%% Finally, reformat the data slightly
# Split the dataframe into a PWM dataframe and a relative disparity dataframe
pwms_by_yr = pwms[['ASIAN', 'BLACK', 'HISPANIC', 'OTHER', 'WHITE', 'TOTAL',
                   'AB617', 'SB535', 'YEAR']].copy()
disp_by_yr = pwms[['ASIAN_REL_DISP', 'BLACK_REL_DISP', 'HISPANIC_REL_DISP',
                   'OTHER_REL_DISP', 'WHITE_REL_DISP', 'AB617_REL_DISP', 
                   'SB535_REL_DISP', 'YEAR']].copy()

# Melt both dataframes to streamline seaborn plotting
pwms_by_yr = pwms_by_yr.melt(id_vars='YEAR', var_name='GROUP', value_name='PWM')
disp_by_yr = disp_by_yr.melt(id_vars='YEAR', var_name='GROUP_RE', value_name='RELATIVE_DISPARITY')
disp_by_yr['GROUP'] = disp_by_yr['GROUP_RE'].str.split('_').str[0]

#%% Create Figure 1
# Initialize figure
fig, ax = plt.subplots(1,2, figsize=(9,4))

# # # # # # # #
# Panel (a)   #
# # # # # # # #

# Plot the data as lines per group
sns.lineplot(data=pwms_by_yr[pwms_by_yr['GROUP']!='TOTAL'], x='YEAR', y='PWM', 
             hue='GROUP', palette=palette, legend=False, ax=ax[0])

# Add the TOTAL line as a dashed line
ax[0].plot(pwms_by_yr.loc[pwms_by_yr['GROUP']=='TOTAL', 'YEAR'], 
           pwms_by_yr.loc[pwms_by_yr['GROUP']=='TOTAL', 'PWM'], 
           '--', color='black', linewidth=1.5)

# Add formatting for the x-axis
ax[0].set_xticks(ticks=years, labels=['2000','','','','','2005','','','','',
                                      '2010','','','','','2015','','','',''])
ax[0].set_xlim([years.min(), years.max()])
ax[0].set_xlabel(r'Year', fontsize=13)

# Add formatting for the y-axis
ax[0].set_ylabel('Population-Weighted Mean Exposure\nto On-Road Mobile Source'+r' PM$_{2.5}$ ($\mu$g/m$^3$)', fontsize=13)
ax[0].set_ylim([0,4.5])
ax[0].set_yticks(ticks=np.arange(0,5,0.5), labels=['0','','','1.5','','','3.0','','','4.5'])

# Add a thicker frame for both panels
for axis in ['top','bottom','left','right']:
    ax[0].spines[axis].set_linewidth(2)
    ax[1].spines[axis].set_linewidth(2)
    
# # # # # # # #
# Panel (b)   #
# # # # # # # #

# Plot the data as lines per group
sns.lineplot(data=disp_by_yr, x='YEAR', y='RELATIVE_DISPARITY', 
             hue='GROUP', palette=palette, legend=False, ax=ax[1])

# Add a zero line as a dashed line for reference
ax[1].plot(np.array([2000,2019]),np.array([0,0]), '--', color='black', linewidth=1.5)
# Add formatting for the x-axis
ax[1].set_xticks(ticks=years, labels=['2000','','','','','2005','','','','',
                                      '2010','','','','','2015','','','',''])
ax[1].set_xlim([years.min(), years.max()])
ax[1].set_xlabel(r'Year', fontsize=13)

# Add formatting for the y-axis
ax[1].set_ylabel('Relative Disparity in Exposure to\nOn-Road Mobile Source'+r' PM$_{2.5}$ (%)', fontsize=13)
ax[1].set_ylim([-50,50])
ax[1].set_yticks(ticks=np.arange(-50,51,12.5), labels=['-50','','-25','','0',
                                                       '','25','','50'])
# Add a thicker frame for both panels
for axis in ['top','bottom','left','right']:
    ax[0].spines[axis].set_linewidth(2)
    ax[1].spines[axis].set_linewidth(2)


# Final clean up
fig.tight_layout()
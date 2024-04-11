#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PM2.5 exposure disparities persist despite strict vehicle emissions controls in California
Koolik et al. (2024)

Figure 3

"""

# Import libraries
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(context=None, font_scale=1, rc=None, style='ticks')
from cmcrameri import cm
import numpy as np
from os import path

#%% Define Critical User Paths Up Front
# Folder containing exposure shapefiles 
main_data_path = '/Users/libbykoolik/Documents/Research/OEHHA Project/Retrospective Analysis/Final_Scripts/data_to_upload'

#%% Global definitions
# Define the vehicle types run through ECHO-AIR
vehicle_types = ['all','ldv','mdv','hdv']

# Define a list of demographic groups
demo_groups = ['TOTAL', 'ASIAN', 'BLACK', 'HISLA', 'WHITE', 'OTHER']

# Define the range of time periods
years = range(2000,2020)

#%% Create functions to expedite processing
def estimate_total_pwm(c20xx, group):
    ''' Returns the total population-weighted mean exposure for a group for an exposure year '''
    # Create a copy of the dataframe
    c20xx = c20xx.copy()
    
    # Estimate the population-weighted mean exposure
    pwm = (c20xx['PM25_UG_M3'] * c20xx[group]).sum() / (c20xx[group].sum())
    
    return pwm

def process_c20xx(year, source, demo_groups=demo_groups, main_data_path=main_data_path):
    ''' Runs functions for a year of exposure data '''
    
    # Load the data
    c20xx = gpd.read_file(path.join(main_data_path,
                                    '{}_cy{}_exposure_concentrations.shp'.format(source,year)))
    
    # Create OTHER and remove negative population counts caused by rounding
    c20xx['OTHER'] = c20xx['TOTAL'] - c20xx['ASIAN'] - c20xx['BLACK'] - c20xx['HISLA'] - c20xx['WHITE']
    c20xx.loc[c20xx['OTHER']<0, 'OTHER'] = 0
    
    # Create small geodataframe for groups
    demo_groups = ['TOTAL', 'ASIAN', 'BLACK', 'HISLA', 'WHITE', 'OTHER']
    pwms = {}
    
    # Get PWM for each demographic group
    for group in demo_groups:
        pwms[group] = [estimate_total_pwm(c20xx, group)]
        
    # Convert into dataframe
    pwm_dataframe = pd.DataFrame(pwms)
    
    # Add some formatting
    pwm_dataframe['YEAR'] = year
    pwm_dataframe['SOURCE'] = source.upper()
    
    return pwm_dataframe

#%% Load and compile the data
# Create a list for storing data
pwms = []

# Iterate through each vehicle type and year to get the PWMs
for vt in vehicle_types: # Loop through each vehicle type
    for year in years: # Loop through each year
        pwms.append(process_c20xx(year, vt))
        
#%% Process the population-weighted mean data
# Create one big dataframe
pwm_by_year = pd.concat(pwms, ignore_index=True).reset_index(drop=True)

# Melt this dataframe to calculate disparities
pwm_by_year = pd.melt(pwm_by_year, id_vars=['YEAR','SOURCE','TOTAL'], var_name='GROUP',
                      value_name='PWM') 

# Rename the TOTAL column for clarity
pwm_by_year.rename(columns={'TOTAL':'TOTAL_PWM'}, inplace=True)

#%% Need to calculate the "OTHER" source
# Create a copy of the dataframe by grabbing the LDV rows
other_pwms = pwm_by_year[pwm_by_year['SOURCE']=='LDV'][['YEAR','SOURCE','TOTAL_PWM','GROUP']].copy()

# Update the source name
other_pwms['SOURCE'] = 'OTH'

# Write a simple helper function
def calc_oth(year, group, pwm_by_year=pwm_by_year):
    ''' Helper function that estimates the impacts from all other vehicles '''
    # Trim the dataframe to the corresponding year and group
    tmp = pwm_by_year[(pwm_by_year['YEAR']==year)&(pwm_by_year['GROUP']==group)].copy()
    
    # Get the result for the full fleet
    all_pwm = tmp.loc[tmp['SOURCE']=='ALL','PWM'].sum()
    
    # Substract the result from LDV, MDV, and HDV from the full fleet PWM
    oth_pwm = all_pwm - tmp.loc[tmp['SOURCE']!='ALL','PWM'].sum()

    return oth_pwm

# Apply this function on the new oth_pwms dataframe
other_pwms['PWM'] = other_pwms.apply(lambda x: calc_oth(x['YEAR'], x['GROUP']), axis=1)

# Combine this dataframe with the pwm_by_year
pwm_by_year = pd.concat([pwm_by_year, other_pwms], ignore_index=True).reset_index(drop=True)

#%% Calculate the disparity
# Copy the dataframe for a new object
disparities = pwm_by_year.copy()

# Calculate absolute disparity
disparities['ABSOLUTE_DISP'] = disparities['TOTAL_PWM'] - disparities['PWM']

# Calculate relative disparity
disparities['RELATIVE_DISP'] = disparities['ABSOLUTE_DISP'] / disparities['TOTAL_PWM']
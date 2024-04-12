#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PM2.5 exposure disparities persist despite strict vehicle emissions controls in California
Koolik et al. (2024)

Data Processing Scripts

B: Process Data by Policy Designation

In this script, we summarize the population-weighted mean exposure concentration and
exposure disparity for AB617 and SB535 community residents for each vehicle type. 

"""

# Import libraries
import pandas as pd # v1.4.2
import geopandas as gpd # v0.10.2
import numpy as np  # v1.22.3
from os import path 

#%% Define Critical User Paths Up Front
# Folder containing exposure shapefiles 
main_data_path = '/Users/libbykoolik/Documents/Research/OEHHA Project/Retrospective Analysis/Final_Scripts/data_to_upload'

# Paths to necessary policy boundaries
ab617_fp = '/Users/libbykoolik/Documents/Research/OEHHA Project/data/CAPP Communities (dl 2023-09-01)/AB617_Shapes_Edited/ab617_communities.shp'
sb535_fp = '/Users/libbykoolik/Documents/Research/OEHHA Project/data/AB617 and SB535 from Lucas (dl 2023-06-14)/SB535_TRACTS/SB635_tracts_gdf.shp'

# Path to raw population data (can be downloaded from the ECHO-AIR materials)
pop_fp = '/Users/libbykoolik/Documents/Research/OEHHA Project/scripts/isrm_health_calculations/data/ca2010.feather'

# Folder to save the generated CSV files
out_path = '/Users/libbykoolik/Documents/Research/OEHHA Project/Retrospective Analysis/Final_Scripts/processed_data'

#%% Global definitions
# Define the vehicle types run through ECHO-AIR
vehicle_types = ['all','ldv','mdv','hdv']

# Define a list of policy groups
policy_groups = ['TOTAL', 'AB617', 'SB535']

# Define the range of time periods
years = range(2000,2020)

#%% Load the political boundary files and the population data and do some minor processing
# Read both boundary files
ab617 = gpd.read_file(ab617_fp)
sb535 = gpd.read_file(sb535_fp)

# Update the naming on the SB535 files for consistency with AB617
sb535 = sb535[['ApproxLoc','geometry']].copy()
sb535.rename(columns={'ApproxLoc':'Name'}, inplace=True)

# Read the population data
pop = gpd.read_feather(pop_fp)

# Reduce the total data dimensions
pop_geo = pop[['POP_ID','geometry']].drop_duplicates() # extract geo data
pop = pop.groupby(['POP_ID','YEAR'])[['TOTAL']].sum().reset_index() # sum across age bins
pop = gpd.GeoDataFrame(pd.merge(pop, pop_geo, on='POP_ID')) # merge for one smaller dataframe

#%% Project all three to match the CRS of the exposure data
# Load one file briefly to get the CRS
crs_to_use = gpd.read_file(path.join(main_data_path,
                                     'all_cy2000_exposure_concentrations.shp')).crs

# Project each geodataframe
pop = pop.to_crs(crs_to_use)
ab617 = ab617.to_crs(crs_to_use)
sb535 = sb535.to_crs(crs_to_use)

#%% There are a few problems with the underlying geometries that need to be resolved
# Write a helper function to remove multipolygons
def remove_multis(data, name, id_col, cols):
    ''' Removes the multipolygons from a geodataframe '''
    # Get rid of multipolygons
    data = data.explode(index_parts=True).reset_index(drop=True).reset_index()
    
    # Create a codename
    data[name+'_ID'] = data[id_col].astype(str)+'_'+data['index'].astype(str)

    # Clean up 
    data = data[[name+'_ID', 'geometry']+cols].copy()

    return data

# Run on each dataset
ab617 = remove_multis(ab617, 'AB617', 'code', ['code','Community'])
sb535 = remove_multis(sb535, 'SB535', 'Name', ['Name'])

# Buffer the population data to help with donut holes
pop['geometry'] = pop.buffer(0)

#%% Define a useful function for geoprocessing data
def intersect_data(data, boundary, intrinsic_vars, extrinsic_vars, name='', how='union'):
    ''' Function for intersecting data and boundary and allocating the variables 
         INPUTS:
             - data: the dataset that has variables that are being apportioned
             - boundary: the dataset that has boundaries but no variables 
               to be apportioned
             - intrinsic_vars: the names of the columns that are intrinsic 
               (e.g., concentration)
             - extrinsic_vars: the names of the columns that are extrinsic 
               (e.g., population)
             - name: the name to assign the filter column (e.g., AB617)'''
    
    # Create a copy of each dataset
    data_copy = data.copy(deep=True)
    boundary_copy = boundary.copy(deep=True)
    
    # If there are only intrinsic variables, make a copy
    if len(intrinsic_vars) > 0 and len(extrinsic_vars) == 0:
        data_copy = data_copy[intrinsic_vars + ['geometry']].copy()
        
    # Add an ID field to the dataset being apportioned
    data_copy['DATA_ID'] = np.arange(data_copy.shape[0])
    
    # Add an ID field to the boundary data
    if len(name) > 0:
        boundary_copy['FILTER'] = name.upper()

    ## Perform Intersect 
    # Create intersect object between data and boundary
    intersect = gpd.overlay(data_copy, boundary_copy, keep_geom_type=False, how=how)

    # Add the boundary ID where it was nan
    if len(name) > 0:
        intersect['FILTER'] = intersect['FILTER'].fillna('NOT_'+name.upper())
    
    # Some incorrect objects can be created, so remove these before we scale variables
    intersect = intersect[~intersect['POP_ID'].isna()].copy()
    intersect = intersect[intersect.geometry.area >= 10**-8]
    
    # If there are extrinsic values, apportion these using area fractions
    if len(extrinsic_vars) > 0:
        intersect['area_km2'] = intersect.geometry.area/(1000.0*1000.0)    
        data_total_area = intersect.groupby('DATA_ID').sum()['area_km2'].to_dict()
    
        # Add a total area and area fraction to the intersect object
        intersect['area_total'] = intersect['DATA_ID'].map(data_total_area)
        intersect['area_frac'] = intersect['area_km2'] / intersect['area_total']

        # Update the extrinsic values to scale by the area fraction
        for var in extrinsic_vars:
            intersect[var] = intersect['area_frac'] * intersect[var]  
            
            # Remove any null variables
            intersect[var] = intersect[var].fillna(0)
        
        ## Clean up
        intersect = intersect.drop(columns=['area_km2','area_total','area_frac'])
        
    # Check for intrinsic variables and get rid of nulls
    if len(intrinsic_vars) > 0:
        for var in intrinsic_vars:
            intersect[var] = intersect[var].fillna(0)
        
    # Clean up
    intersect = intersect.drop(columns=['DATA_ID'])
    intersect = intersect.reset_index(drop=True)
    
    return intersect

#%% Create intersections between population and EJ boundaries
# Run this function on the AB617 and SB535 boundaries
# These new objects are not aggregated at the community level, providing a high 
# level of granularity for future calculations
ab617_pop = intersect_data(pop, ab617, [], ['TOTAL'], name='AB617', how='union')
sb535_pop = intersect_data(pop, sb535, [], ['TOTAL'], name='SB535', how='union')

# Do a quick clean-up on boundaries to remove donut holes again
ab617_pop['geometry'] = ab617_pop.buffer(0)
sb535_pop['geometry'] = sb535_pop.buffer(0)

#%% Create a crosswalk for Concentration-Population-EJ Designation
# Note that the ISRM ID for each year is preserved, so we only need to intersect 
# the population-EJ boundary files with the ISRM grid once
def get_crosswalk(ej_pop, c20xx, filter_name):
    ''' '''
    # Perform an intersection between these objects
    ej_pop_isrm = intersect_data(ej_pop[ej_pop['FILTER']==filter_name].copy(), 
                                 c20xx, [], ['TOTAL'], how='intersection')
    
    return ej_pop_isrm

# Load one file and trim to only necessary columns
tmp = gpd.read_file(path.join(main_data_path,
                                     'all_cy2000_exposure_concentrations.shp'))
tmp = tmp[['ISRM_ID','geometry']].copy()

# Run on AB617-Population and SB535-Population objects
ab617_pop_isrm = get_crosswalk(ab617_pop, tmp, 'AB617')
sb535_pop_isrm = get_crosswalk(sb535_pop, tmp, 'SB535')

#%% Create a function for loading exposure data and combining with these crosswalks
#   Then, run for each year and source for both AB617 and SB535
def process_c20xx(year, source, ej_pop_isrm, main_data_path=main_data_path):
    ''' Runs functions for a year of exposure data from a source'''
    
    # Load the data
    c20xx = gpd.read_file(path.join(main_data_path,
                                    '{}_cy{}_exposure_concentrations.shp'.format(source,year)))
    
    # Estimate the statewide PWM
    total_pwm = (c20xx['TOTAL']*c20xx['PM25_UG_M3']).sum() / (c20xx['TOTAL'].sum())
    
    # Remove unnecessary columns
    c20xx = c20xx[['ISRM_ID','PM25_UG_M3']].copy()
    
    # Merge with the crosswalk
    ej_pop_c20xx = pd.merge(ej_pop_isrm, c20xx, on='ISRM_ID')
    
    # Calculate the PWM
    pwm = (ej_pop_c20xx['TOTAL']*ej_pop_c20xx['PM25_UG_M3']).sum()/(ej_pop_c20xx['TOTAL'].sum())
    
    # Create a list so it's easy to turn into a dataframe later
    pwm_list = [year, source, pwm, total_pwm]
    
    return pwm_list

# Create lists for storing data
ab617_pwms = []
sb535_pwms = []

# Iterate through each vehicle type and year to get the PWMs
for vt in vehicle_types: # Loop through each vehicle type
    for year in years: # Loop through each year
        ab617_pwms.append(process_c20xx(year, vt, ab617_pop_isrm))
        sb535_pwms.append(process_c20xx(year, vt, sb535_pop_isrm))
        
#%% Compile all of this data together
# Turn each list into a dataframe
ab617_pwms = pd.DataFrame(ab617_pwms, 
                          columns=['YEAR','SOURCE','PWM','TOTAL_PWM']).reset_index(drop=True)
sb535_pwms = pd.DataFrame(sb535_pwms, 
                          columns=['YEAR','SOURCE','PWM','TOTAL_PWM']).reset_index(drop=True)

# Merge these dataframes together
ej_pwms = pd.merge(ab617_pwms, sb535_pwms, on=['YEAR','SOURCE'],
                   suffixes = ("_AB617", "_SB535"))

# Confirm that the totals are equivalent
assert ((ej_pwms['TOTAL_PWM_AB617'] / ej_pwms['TOTAL_PWM_SB535']).sum() == ej_pwms.shape[0])

# Clean up this dataframe 
ej_pwms.rename(columns={'TOTAL_PWM_AB617':'TOTAL_PWM', 
                        'PWM_AB617':'AB617',
                        'PWM_SB535':'SB535'}, inplace=True)
ej_pwms = ej_pwms[['YEAR','SOURCE','TOTAL_PWM','AB617', 'SB535']].copy()
ej_pwms['SOURCE'] = ej_pwms['SOURCE'].str.upper()

#%% Need to calculate the "OTHER" source
# Create a copy of the dataframe by grabbing the LDV rows
other_pwms = ej_pwms[ej_pwms['SOURCE']=='LDV'][['YEAR','SOURCE','TOTAL_PWM','AB617','SB535']].copy()

# Update the source name
other_pwms['SOURCE'] = 'OTH'

# Write a simple helper function
def calc_oth(year, group, ej_pwms=ej_pwms):
    ''' Helper function that estimates the impacts from all other vehicles '''
    
    # Trim the dataframe to the corresponding year 
    tmp = ej_pwms[(ej_pwms['YEAR']==year)].copy()
    
    # Get the result for the full fleet
    all_pwm = tmp.loc[tmp['SOURCE']=='ALL',group].sum()
    
    # Substract the result from LDV, MDV, and HDV from the full fleet PWM
    oth_pwm = all_pwm - tmp.loc[tmp['SOURCE']!='ALL',group].sum()

    return oth_pwm

# Apply this function on the new oth_pwms dataframe
other_pwms['AB617'] = other_pwms.apply(lambda x: calc_oth(x['YEAR'], 'AB617'), axis=1)
other_pwms['SB535'] = other_pwms.apply(lambda x: calc_oth(x['YEAR'], 'SB535'), axis=1)

# Need to update the TOTAL PWM for the OTH group
other_total_pwm = other_pwms[['YEAR']].drop_duplicates().copy()
other_total_pwm['TOTAL_PWM'] = other_total_pwm.apply(lambda x: calc_oth(x['YEAR'], 'TOTAL_PWM'), axis=1)

# Merge together
other_pwms = pd.merge(other_pwms[['YEAR','SOURCE','AB617','SB535']], other_total_pwm[['YEAR','TOTAL_PWM']], on='YEAR')
other_pwms = other_pwms[['YEAR','SOURCE','TOTAL_PWM','AB617','SB535']].copy()

# Combine this dataframe with the pwm_by_year
ej_pwms = pd.concat([ej_pwms, other_pwms], ignore_index=True).reset_index(drop=True)

#%% Calculate the disparity
# Copy the dataframe for a new object
disparities = ej_pwms.copy()

# Calculate absolute disparity for each EJ group
disparities['AB617_ABSOLUTE_DISP'] = disparities['AB617'] - disparities['TOTAL_PWM']
disparities['SB535_ABSOLUTE_DISP'] = disparities['SB535'] - disparities['TOTAL_PWM']

# Calculate relative disparity
disparities['AB617_RELATIVE_DISP'] = disparities['AB617_ABSOLUTE_DISP'] / disparities['TOTAL_PWM']
disparities['SB535_RELATIVE_DISP'] = disparities['SB535_ABSOLUTE_DISP'] / disparities['TOTAL_PWM']

#%% Output the disparities dataframe to a CSV
disparities.to_csv(path.join(out_path, 'disparities_by_policy.csv'), index=False)
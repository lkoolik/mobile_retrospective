#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PM2.5 exposure disparities persist despite strict vehicle emissions controls in California
Koolik et al. (2024)

Data Processing Scripts

B: Process Data by Policy Designation

In this script, we summarize the population-weighted mean exposure concentration and
exposure disparity from each vehicle type for each year for 
    - AB617 community residents in aggregate
    - SB535 community residents in aggregate
    - Ab617 community residents by community
    - Regional boundaries defined in Koolik et al. (2024)

"""

# Import libraries
import pandas as pd # v1.4.2
import geopandas as gpd # v0.10.2
import numpy as np  # v1.22.3
from os import path 

#%% Define Critical User Paths Up Front
# Folder containing exposure shapefiles (e.g., all_cy2000_exposure_concentrations.shp)
main_data_path = '...' # add path

# Paths to necessary political boundaries
ab617_fp = '.../ab617_communities.shp' # add path
sb535_fp = '.../SB635_tracts_gdf.shp' # add path
counties_fp = '.../counties.feather' # add path (can be downloaded from the ECHO-AIR materials)

# Path to raw population data (can be downloaded from the ECHO-AIR materials)
pop_fp = '.../ca2010.feather' # add path

# Folder to save the generated CSV files
out_path = '...' # add path

#%% Global definitions
# Define the vehicle types run through ECHO-AIR
vehicle_types = ['all','ldv','mdv','hdv']

# Define a list of policy groups
policy_groups = ['TOTAL', 'AB617', 'SB535']

# Define the range of time periods
years = range(2000,2020)

#%% Load the political boundary files and the population data and do some minor processing
# Read the boundary files
ab617 = gpd.read_file(ab617_fp)
sb535 = gpd.read_file(sb535_fp)
counties = gpd.read_feather(counties_fp)

# Make a few cosmetic updates to these dataframes
sb535 = sb535[['ApproxLoc','geometry']].copy()
sb535.rename(columns={'ApproxLoc':'Name'}, inplace=True)

# Read the population data
pop = gpd.read_feather(pop_fp)

# Reduce the total data dimensions
pop_geo = pop[['POP_ID','geometry']].drop_duplicates() # extract geo data
pop = pop.groupby(['POP_ID','YEAR'])[['TOTAL']].sum().reset_index() # sum across age bins
pop = gpd.GeoDataFrame(pd.merge(pop, pop_geo, on='POP_ID')) # merge for one smaller dataframe

#%% We are more interested in regional results than county-specific results
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

# Remove the 'NA' values
# counties = counties[counties['REGION'].isin(['SF','CV','LA'])].copy()

# Dissolve counties to get a regional dataframe
regions = counties[['REGION','geometry']].dissolve(by='REGION').reset_index()

#%% Project all boundary files to match the CRS of the exposure data
# Load one file briefly to get the CRS
crs_to_use = gpd.read_file(path.join(main_data_path,
                                     'all_cy2000_exposure_concentrations.shp')).crs

# Project each geodataframe
pop = pop.to_crs(crs_to_use)
ab617 = ab617.to_crs(crs_to_use)
sb535 = sb535.to_crs(crs_to_use)
regions = regions.to_crs(crs_to_use)

#%% There are a few problems with the underlying geometries that need to be resolved
# Write a helper function to remove multipolygons
def remove_multis(data, name, id_col, cols):
    ''' Removes the multipolygons from a geodataframe
        INPUTS:
            - data = geodataframe containing multipolygons
            - name = string used to assign a unique ID to new polygons
            - id_col = column containing former ID 
            - cols = columns to keep in the new geodataframe
            
        OUTPUTS:
            - data = updated, simplified geodataframe  '''
            
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
regions = remove_multis(regions, 'REGION', 'REGION', ['REGION'])

# Buffer the population data to help with donut holes
pop['geometry'] = pop.buffer(0)

#%% Define a useful function for geoprocessing data
def intersect_data(data, boundary, intrinsic_vars, extrinsic_vars, name='', how='union'):
    ''' Function for intersecting data and boundary and allocating the variables 
         INPUTS:
             - data = the dataset that has variables that are being apportioned
             - boundary = the dataset that has boundaries but no variables 
               to be apportioned
             - intrinsic_vars = the names of the columns that are intrinsic 
               (e.g., concentration)
             - extrinsic_vars = the names of the columns that are extrinsic 
               (e.g., population)
             - name = the name to assign the filter column (e.g., AB617)
             - how = method for combining geospatial information
         
            OUTPUTS:
             - intersect: a geodataframe that contains smaller polygons representing 
               the intersection of the data and the boundary '''
    
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
# Run this function on the boundary files
# These new objects are not aggregated at the community level, providing a high 
# level of granularity for future calculations
ab617_pop = intersect_data(pop, ab617, [], ['TOTAL'], name='AB617', how='union')
sb535_pop = intersect_data(pop, sb535, [], ['TOTAL'], name='SB535', how='union')
region_pop = intersect_data(pop, regions, [], ['TOTAL'], name='REGION', how='intersection')

# Do a quick clean-up on boundaries to remove donut holes again
ab617_pop['geometry'] = ab617_pop.buffer(0)
sb535_pop['geometry'] = sb535_pop.buffer(0)
region_pop['geometry'] = region_pop.buffer(0)

#%% Create a crosswalk for Concentration-Population-EJ Designation
# Note that the ISRM ID for each year is preserved, so we only need to intersect 
# the population-EJ boundary files with the ISRM grid once
def get_crosswalk(ej_pop, c20xx, filter_name):
    ''' Creates a geospatial intersection between the EJ boundary, population data,
        and ISRM geography
        INPUTS:
            - ej_pop = geodataframe containing the intersection between the EJ
              boundary and the population data
            - c20xx = geodataframe containing concentration at the ISRM grid cell
            
        OUTPUTS:
            - ej_pop_isrm = geodataframe containing the intersection between the 
              EJ boundary, the population data, and the ISRM ID '''
            
    # Perform an intersection between these objects
    ej_pop_isrm = intersect_data(ej_pop[ej_pop['FILTER']==filter_name].copy(), 
                                 c20xx, [], ['TOTAL'], how='intersection')
    
    return ej_pop_isrm

# Load one file and trim to only necessary columns
tmp = gpd.read_file(path.join(main_data_path,
                                     'all_cy2000_exposure_concentrations.shp'))
tmp = tmp[['ISRM_ID','geometry']].copy()

# Run on AB617-Population, SB535-Population, and County-Population objects
ab617_pop_isrm = get_crosswalk(ab617_pop, tmp, 'AB617')
sb535_pop_isrm = get_crosswalk(sb535_pop, tmp, 'SB535')
region_pop_isrm = get_crosswalk(region_pop, tmp, 'REGION')

#%% Make a few quick updates to columns for simplifying later functions
# Rename the AB617 columns
ab617_pop_isrm.rename(columns={'code':'CODE', 'Community':'NAME'}, inplace=True)

# Make the region dataframe match the AB617 dataframe structure
region_pop_isrm['NAME'] = region_pop_isrm['REGION'].map({'CV':'Central Valley',
                                                         'LA':'Los Angeles Area',
                                                         'SF':'San Francisco Bay Area'})
region_pop_isrm.rename(columns={'REGION':'CODE'}, inplace=True)

# Save these as a crosswalks dictionary
crosswalks = {'AB617':ab617_pop_isrm,
              'SB535':sb535_pop_isrm,
              'REGION':region_pop_isrm}

#%% Create a functions for loading exposure data, combining with these crosswalks,
#   and calculating population-weighted mean exposure concentrations
def get_pwm(data, pop_col):
    '''Estimates the population-weighted mean exposure from a data set
        INPUTS:
            - data = dataframe containing population counts and concentration
            - pop_col = column name for population counts
            
        OUTPUTS:
            - pwm = the population-weighted mean concentration '''
            
    # Estimate the PWM
    pwm = (data['PM25_UG_M3']*data[pop_col]).sum() / (data[pop_col].sum())
    
    return pwm

# As part of this pre-processing, we should also get community-specific PWMs 
# for each year and source for AB617
def get_detailed_pwms(year, source, ej_pop_c20xx):
    ''' Estimates the population-weighted mean exposure for each year for each source
        for each AB617 community or each region
        INPUTS:
            - year = year corresponding to modeled data
            - source = vehicle type group name corresponding to modeled data
            - ej_pop_c20xx = geodataframe containing the intersection between the 
              boundary, the population data, and the concentration
            
        OUTPUTS:
            - detailed_pwms = a dataframe containing the year, the source, 
              and the population-weighted mean concentration for each AB617 
              community or each region''' 
              
    # Create a dataframe of community names
    detailed_pwms = ej_pop_c20xx[['NAME','CODE']].drop_duplicates().copy().reset_index(drop=True)
    
    # Add details to the dataframe
    detailed_pwms['YEAR'] = year
    detailed_pwms['SOURCE'] = source.upper()
    
    # Call the PWM function for each community
    detailed_pwm_tmp = {}
    for code in detailed_pwms['CODE'].unique():
        tmp = ej_pop_c20xx[ej_pop_c20xx['CODE']==code].copy()
        detailed_pwm_tmp[code] = get_pwm(tmp, 'TOTAL')
        
    # Add this dictionary as a column
    detailed_pwms['PWM'] = detailed_pwms['CODE'].map(detailed_pwm_tmp)
    
    # Clean up
    detailed_pwms = detailed_pwms[['YEAR','SOURCE','NAME','CODE','PWM']]
    
    return detailed_pwms

def process_c20xx(year, source, crosswalks=crosswalks, main_data_path=main_data_path):
    ''' Estimates the population-weighted mean exposure for each year for each source
        INPUTS:
            - year = year corresponding to modeled data
            - source = vehicle type group name corresponding to modeled data
            - crosswalks = a dictionary containing each of the crosswalks between
              political boundaries, population, and ISRM grid cell
            - main_data_path = filepath where all ECHO-AIR outputs are stored
            
        OUTPUTS:
            - pwms = a list or dataframe containing the year, the source, the population-
              weighted mean concentration for the EJ group, and the total 
              statewide population-weighted mean concentration '''
    
    # Load the data
    c20xx = gpd.read_file(path.join(main_data_path,
                                    '{}_cy{}_exposure_concentrations.shp'.format(source,year)))
    
    # Estimate the statewide PWM
    total_pwm = get_pwm(c20xx, 'TOTAL')
    
    # Remove unnecessary columns
    c20xx = c20xx[['ISRM_ID','PM25_UG_M3']].copy()
    
    # Merge with each crosswalk
    ab617_pop_c20xx = pd.merge(crosswalks['AB617'], c20xx, on='ISRM_ID')
    sb535_pop_c20xx = pd.merge(crosswalks['SB535'], c20xx, on='ISRM_ID')
    region_pop_c20xx = pd.merge(crosswalks['REGION'], c20xx, on='ISRM_ID')
    
    # Calculate the PWM for AB617 and SB535
    ab617_pwm = get_pwm(ab617_pop_c20xx, 'TOTAL')
    sb535_pwm = get_pwm(sb535_pop_c20xx, 'TOTAL')
    
    # Combine these into a dataframe
    ej_pwms = pd.DataFrame([[year, source.upper(), 'TOTAL', 'TOTAL', total_pwm],
                            [year, source.upper(), 'AB617', 'AB617', ab617_pwm],
                            [year, source.upper(), 'SB535', 'SB535', sb535_pwm]],
                           columns = ['YEAR','SOURCE','NAME','CODE','PWM'])
    
    # Calculate the PWM for each AB617 community
    community_pwms = get_detailed_pwms(year, source, ab617_pop_c20xx)
    
    # Calculate the PWM for each region
    regional_pwms = get_detailed_pwms(year, source, region_pop_c20xx)
    
    # Combine all three dataframes
    pwms = pd.concat([ej_pwms, community_pwms, regional_pwms], 
                     ignore_index=True).reset_index(drop=True)
    
    return pwms

#%% Process all of the data
# Create list for storing data
pwms = []

# Iterate through each vehicle type and year to get the PWMs
for vt in vehicle_types: # Loop through each vehicle type
    for year in years: # Loop through each year
        pwms.append(process_c20xx(year, vt))
        
#%% Compile all of this data together
# Turn this list into a dataframe
pwms_by_year = pd.concat(pwms).reset_index(drop=True)

# Save a lookup for the NAME-CODE combinations
lookup = pwms_by_year[['NAME','CODE']].drop_duplicates().copy()
lookup = dict(zip(lookup['CODE'], lookup['NAME']))

# Pivot this data to estimate the contribution from all other sources
pwms_pivot = pwms_by_year[['YEAR','SOURCE','CODE','PWM']].pivot(columns=['SOURCE'],
                                                                index=['YEAR','CODE'], 
                                                                values='PWM').reset_index()
pwms_pivot['OTH'] = pwms_pivot['ALL'] - pwms_pivot['LDV'] - pwms_pivot['MDV'] - pwms_pivot['HDV']

#%% Restructure for a useful output file
# Melt the data to the longest format possible
pwms_df = pwms_pivot.melt(id_vars=['YEAR','CODE'], var_name='SOURCE', value_name='PWM').reset_index(drop=True)

# Extract the totals column
pwms_total = pwms_df[pwms_df['CODE']=='TOTAL'].copy()
pwms_df = pwms_df[pwms_df['CODE']!='TOTAL'].copy()

# Merge the totals back in
pwms_df = pd.merge(pwms_df, pwms_total[['YEAR','SOURCE','PWM']], 
                   on=['YEAR','SOURCE'], suffixes=('','_TOTAL'))

# Add the community names back in
pwms_df['NAME'] = pwms_df['CODE'].map(lookup)

# Remove the extra region and clean up slightly
pwms_df = pwms_df[pwms_df['CODE'] != 'NA'].copy().reset_index(drop=True)
pwms_df = pwms_df[['YEAR','SOURCE','CODE','NAME','PWM_TOTAL','PWM']].copy()

#%% Calculate the disparity
# Copy the dataframe for a new object
disparities = pwms_df.copy()

# Calculate absolute disparity for each EJ group
disparities['ABSOLUTE_DISP'] = disparities['PWM'] - disparities['PWM_TOTAL']

# Calculate relative disparity
disparities['RELATIVE_DISP'] = disparities['ABSOLUTE_DISP'] / disparities['PWM_TOTAL']

#%% Output the disparities dataframe to a CSV
# Update to match the format of the racial-ethnic dataset
disparities.rename(columns={'PWM_TOTAL':'TOTAL_PWM',
                            'NAME':'GROUP'}, inplace=True)
disparities = disparities[['YEAR','SOURCE','TOTAL_PWM','CODE','GROUP','PWM',
                           'ABSOLUTE_DISP','RELATIVE_DISP']].copy()

# Export
disparities.to_csv(path.join(out_path, 'disparities_by_policy.csv'), index=False)
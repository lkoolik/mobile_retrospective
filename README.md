# California Mobile Retrospective

<i>Libby H. Koolik, Álvaro Alvarado, Amy Budahn, Laurel Plummer, Julian D. Marshall, and Joshua S. Apte</i>

This Github repository contains the Python code for the paper titled "PM<sub>2.5</sub> exposure disparities persist despite strict vehicle emissions controls in California", currently in review at Science Advances. For more information on the analysis, including methodological details, please see the paper: [link](https://chemrxiv.org/engage/chemrxiv/article-details/6584780166c1381729bcf0b0).

Below are links to the key tools and datasets used in this paper:

* ECHO-AIR: [https://echo-air-model.github.io/](https://echo-air-model.github.io/)
* California InMAP Source-Receptor Matrix: [https://zenodo.org/record/2589760](https://zenodo.org/records/7548607). (Accessed Mar. 28, 2024)
* Year-2010 race-ethnicity data at Census block level: [https://www.nhgis.org/](https://www.nhgis.org/).

This repository contains scripts for reproducing the main text figures. The scripts are organized into two folders, as described below.

----

## Code Description

### ECHO-AIR Output Processing

Our estimates of PM<sub>2.5</sub> concentration are modeled using [ECHO-AIR](https://echo-air-model.github.io/), a novel open source pipeline based on the InMAP Source-Receptor Matrix. All of the relevant data outputs from ECHO-AIR are included in our data repository. 

To simplify the scripts used to reproduce the main text figures, we have first included [two pre-processing scripts](https://github.com/lkoolik/mobile_retrospective/tree/main/01_ECHO_AIR_Output_Processing) for summarizing these gridded concentration estimates.

1. [`a_pre_process_by_race_ethnicity.py`](https://github.com/lkoolik/mobile_retrospective/blob/main/01_ECHO_AIR_Output_Processing/a_pre_process_by_race_ethnicity.py): summarizes the population-weighted mean exposure concentration and exposure disparity for each racial-ethnic group for each vehicle type.
2. [`b_pre_process_by_policy.py`](https://github.com/lkoolik/mobile_retrospective/blob/main/01_ECHO_AIR_Output_Processing/b_pre_process_by_policy.py): summarizes the population-weighted mean exposure concentration and exposure disparity from each vehicle type for each year for AB617 community residents in aggregate, SB535 community residents in aggregate, AB617 community residents by community, and regional boundaries defined in Koolik et al. (2024)

### Main Text Figures

The four figures from the main text of Koolik et al. (2024) are fully reproducible via [four scripts](https://github.com/lkoolik/mobile_retrospective/tree/main/02_Main_Text_Figures) included in this repository. 

1. [`figure01_pwm_and_re.py`](https://github.com/lkoolik/mobile_retrospective/blob/main/02_Main_Text_Figures/figure01_pwm_and_re.py): creates the line plots of population-weighted mean exposure and relative disparity in exposure by each demographic group 
2. [`figure02_exposure_distributions.py`](https://github.com/lkoolik/mobile_retrospective/blob/main/02_Main_Text_Figures/figure02_exposure_distributions.py): creates the distribution of population by race-ethnicity at each decile of exposure to the full vehicle fleet in 2000 and 2019
3. [`figure03_disparity_by_source.py`](https://github.com/lkoolik/mobile_retrospective/blob/main/02_Main_Text_Figures/figure03_disparity_by_source.py): creates the three-panel figure that breaks down absolute and relative disparity by vehicle fleet for Hispanic Californians
4. [`figure04_spatial_heterogeneity.py`](https://github.com/lkoolik/mobile_retrospective/blob/main/02_Main_Text_Figures/figure04_spatial_heterogeneity.py): creates the map of contributions to disparity by vehicle type for the full state, three geographic regions, the two aggregate categories of environmental justice communities, and each of California's AB617 communities

----

### Requirements

All six scripts included in this repository are written in Python 3.10.2. For reproducibility, we have tried to minimize the number of libraries used to generate these figures. Each library (and corresponding version) is listed below.

* cmcrameri (v1.7)
* geopandas (v0.10.2)
* matplotlib (v3.5.1)
* numpy (v1.22.3)
* pandas (v1.4.2)
* seaborn (v0.11.2)

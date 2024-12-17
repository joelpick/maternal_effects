# Code README

**Project title:** Simple maternal effect animal models provide biased estimates of additive genetic and maternal variation
**Code Author:** Joel Pick
**Code License:** MIT



## R Code Files

### 00_functions.R
- Contains the functions that are commonly used in other code files. This does not need to be run - it is called within the other code files when needed

### 01_parameters.R
- Defines the parameters used in the simulations
- Produces *Data/Intermediate/parameters.Rdata* which is used in *02_simulations.R* and *03_small_ped_sims.R*

### 02_simulations.R
- Simulates pedigrees and datasets, and analyses them with a simple maternal effects model, and saves results
- Imports simulation parameters from *Data/Intermediate/parameters.Rdata*
- Produces *Data/Intermediate/mge_sims3.Rdata*
- NOTE: require ASREML licence to run the models

### 03_small_ped_sims.R
- Simulates small pedigrees and datasets, and analyses them with 4 different models, and saves results
- Imports simulation parameters from *Data/Intermediate/parameters.Rdata*
- Produces *mge_sims_small_ped.Rdata*
- NOTE: require ASREML licence to run the models

### 04_synthesis_metrics.R
- Extracts various metrics from data from different synthesis used in this study.
- Imports 
    - *Data/Raw/young_postma_2023_data.csv*
    - *Data/Raw/postma2014_SM.csv*
    - *Data/Raw/moore_2019_datawithSE.csv*
    - *Data/Raw/Bonnet/*

### Fig_2_literature.R
- Generates results from literature survey
- Imports *Data/Raw/MGE - all_est.csv*
- Produces figure 2
    - *Figures/Fig2_lit.pdf*
    - *Figures/FigS2_lit_excluded.pdf*
    - Latex code for Table S1

### Fig3_mat_links.R
- Generates figure 3
- Imports 
    - *Data/Intermediate/mge_sims3.Rdata*
    - *Data/Raw/BT_Pick.csv*
    - *Data/Raw/RD_Gauzere.txt*
- Produces:
    - *Figures/Figures/Fig3_mat_links.pdf*

### Fig4_bias.R
- Generates figures 4 and 5
- Imports 
    - *Data/Intermediate/mge_sims3.Rdata*
- Produces
    - *Figures/Fig4_bias.pdf*
    - *Figures/Fig5_Va_Vm.pdf*
    - *Figures/FigSM_extra_scenarios.pdf*

### Fig6_total_va.R
- Generates figure 6
- Imports 
    - *Data/Intermediate/mge_sims3.Rdata*
- Produces
    - *Figures/Fig6_totalVa.pdf*

### Fig7_small_ped_sims.R
- Generates figure 7
- Imports 
    - *Data/Intermediate/parameters.Rdata*
    - *Data/Intermediate/mge_sims_small_ped.Rdata*
- Produces
    - *Figures/Fig7_small_ped.pdf*

### FigSM_all_sim_bv.R
- Generates figure S
- Imports 
    - 
- Produces
    - 

### FigSM_all_sim.R
- Generates all FigSM_all_sim... figures
- Imports 
    - 
- Produces
    - 


### FigSM_covs.R
- Generates
- Imports 
    - 
- Produces
    - 


### FigSM_links.R
- Generates
- Imports 
    - 
- Produces
    - 

### FigSM_ped_depth.R
- Generates
- Imports 
    - 
- Produces
    - 

### FigSM_ped_links.R
- Generates
- Imports 
    - 
- Produces
    - 

### FigSM_samp_cov.R
- Generates
- Imports 
    - 
- Produces
    - 



## Packages

asreml version 4.1.0

pedAgree version 0.0.1 (installed from github)

squidSim version 0.2.3 (installed from github)

parallel 

beeswarm
scales

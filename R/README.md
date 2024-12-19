# Code README

**Project title:** Simple maternal effect animal models provide biased estimates of additive genetic and maternal variation

**Code Author:** Joel Pick

**Code License:** MIT (see R/LICENSE)



## R Code Files

### 00_functions.R
- Contains the functions that are commonly used in other code files. This does not need to be run - it is called within the other code files when needed

### 01_parameters.R
- Defines the parameters used in the simulations
- Produces *Data/Intermediate/parameters.Rdata* which is used for the simulations in in *02_simulations.R* and *03_small_ped_sims.R*

### 02_simulations.R
- Code for the main set of simulations, looking a the bias in Va in simple maternal effects models. Simulates pedigrees and datasets, and analyses them with a simple maternal effects model, and saves results
- Imports simulation parameters from *Data/Intermediate/parameters.Rdata*
- Produces *Data/Intermediate/mge_sims3.Rdata*
- NOTE: require ASREML licence to run the models. Model output is saved for use in other code files

### 03_small_ped_sims.R
- Code for the main set of simulations, assessing the accuracy of more complex models with small pedigrees. Simulates small pedigrees and datasets, and analyses them with 4 different models, and saves results
- Imports simulation parameters from *Data/Intermediate/parameters.Rdata*
- Produces *mge_sims_small_ped.Rdata*
- NOTE: require ASREML licence to run the models. Model output is saved for use in other code files

### 04_synthesis_metrics.R
- Extracts various metrics from data from different synthesis used in this study.
- Imports 
    - *Data/Raw/young_postma_2023_data.csv*
    - *Data/Raw/postma2014_SM.csv*
    - *Data/Raw/moore_2019_datawithSE.csv*
    - *Data/Raw/Bonnet/*

### Fig_2_literature.R
- Generates results from literature survey; Figure 2 and S3, and Table S1
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
- Generates figures 4 and 5, and S7
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

### FigSM_all_sim_small_ped.R
- Generates Table S2-3, Figures S11-14 and S27-S38
- Imports 
    - mge_sims3.Rdata
    - parameters.Rdata
    - mge_sims_small_ped.Rdata
- Produces
    - all Figures/FigSM_all_sim_bv... figures
    - Figures/FigSM_small_ped_Va.pdf
    - Figures/FigSM_small_ped_tVa.pdf
    - Figures/FigSM_small_ped_Vm.pdf
    - Figures/FigSM_small_ped_convergence.pdf
    - Table S2
    - Table S3

### FigSM_all_sim.R
- Generates Figures S15-S26
- Imports 
    - mge_sims3.Rdata
    - parameters.Rdata
- Produces
    - all Figures/FigSM_all_sim

### FigSM_covs.R
- Generates figures S1-2 showing how the covariance between different relatives is generated.
- Imports 
    - Data/Raw/covariances.csv
- Produces
    - Figures/FigS1_relative_cov.pdf
    - Figures/FigS1_relative_cov_eg.pdf

### FigSM_links.R
- Generates Figure S5
- Imports 
    - mge_sims3.Rdata
- Produces
    - Figures/FigS3_links_cor.pdf

### FigSM_ped_depth.R
- Runs pedigree depth simulations and generates figure S6
- Produces
    - Figures/FigSM_ped_depth.pdf

### FigSM_ped_links.R
- Generates Figure S4
- Imports 
    - mge_sims3.Rdata
- Produces
    - Figures/FigS4_links_ped.pdf

### FigSM_samp_cov.R
- Generates S8-10
- Imports 
    - mge_sims3.Rdata
- Produces
    - Figures/FigSM_samp_cov_mat_links.pdf
    - Figures/FigSM_samp_cov_ped.pdf
    - Figures/FigSM_samp_cov_scenario.pdf



## Packages

asreml version 4.1.0 - Note this is propriety software and requires a license

pedAgree version 0.0.1 (installed from github)

squidSim version 0.2.3 (installed from github)

parallel version xxx

beeswarm version xxx

scales version xxx

viridis version xxx
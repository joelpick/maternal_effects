rm(list=ls())

library(asreml)

wd <- "/Users/joelpick/github/maternal_effects/"
data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))

## loads fs_peds,hs_peds,fs_data,hs_data
load(paste0(data_wd,"gaussian_data.Rdata"))

model1_fs <- model_func(m1_func,fs_peds,fs_data)
model1_hs <- model_func(m1_func,hs_peds,hs_data)
model2_fs <- model_func(m2_func,fs_peds,fs_data)
model2_hs <- model_func(m2_func,hs_peds,hs_data)
model3_fs <- model_func(m3_func,fs_peds,fs_data)
model3_hs <- model_func(m3_func,hs_peds,hs_data)
model4_fs <- model_func(m4_func,fs_peds,fs_data)
model4_hs <- model_func(m4_func,hs_peds,hs_data)
model5_fs <- model_func(m5_func,fs_peds,fs_data) # many had this error: Likelihood evaluation failed with fault 1009 ; trying with reduced updates - This means it hasn't worked
model5_hs <- model_func(m5_func,hs_peds,hs_data)# Likelihood evaluation failed with fault 1009 ; trying with reduced updates
model6_fs <- model_func(m6_func,fs_peds,fs_data)
model6_hs <- model_func(m6_func,hs_peds,hs_data)
model7_fs <- model_func(m7_func,fs_peds,fs_data)
model7_hs <- model_func(m7_func,hs_peds,hs_data)
model8_fs <- model_func(m8_func,fs_peds,fs_data)
model8_hs <- model_func(m8_func,hs_peds,hs_data)

save(
	model1_hs,model1_fs,
	model2_hs,model2_fs,
	model3_hs,model3_fs,
	model4_hs,model4_fs, 
	model5_hs, model5_fs,
	model6_hs, model6_fs,
	model7_hs, model7_fs,
	model8_hs, model8_fs,  
	file=paste0(data_wd,"gaussian_sims.Rdata")
)

ped_bt <- read.csv(paste0(wd,"Data/Raw/ped_BT.csv"))[]
ped_rd <- read.csv(paste0(wd,"Data/Raw/ped_RD.csv"))
ped_ss <- read.csv(paste0(wd,"Data/Raw/ped_SSH.csv"))

lapply()
head(ped_bt)
model1_fs <- model_func(m1_func,ped_bt,bt_data)




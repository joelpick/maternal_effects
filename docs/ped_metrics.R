rm(list=ls())
wd <- "/Users/joelpick/github/maternal_effects/"

# devtools::install_github("cran/pedantics")
library(pedantics)

# Difference between blue tit and red deer pedigree

ped_bt <-  MasterBayes::insertPed(read.csv(paste0(wd,"Data/Raw/ped_BT.csv"))[,1:3])

bt<-pedigreeStats(ped_bt, includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE)

pedStatSummary(bt)
# maternities is mother-offspring links, but not unique mothers

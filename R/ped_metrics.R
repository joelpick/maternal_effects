rm(list=ls())
wd <- "/Users/joelpick/github/maternal_effects/"

source(paste0(wd,"R/extract_cousins.R"))
# devtools::install_github("cran/pedantics")
library(pedantics)

# Difference between blue tit and red deer pedigree

ped_bt <-  MasterBayes::insertPed(read.csv(paste0(wd,"Data/Raw/ped_BT.csv"))[,1:3])
ped_RD<-read.csv(paste0(wd,"Data/Raw/ped_RD.csv"))
ped_SSH<-read.csv(paste0(wd,"Data/Raw/ped_SSH.csv"))
ped_SFW<-read.csv(paste0(wd,"Data/Raw/ped_SFW.csv"))
ped_SV<-read.csv(paste0(wd,"Data/Raw/ped_SV.csv"))



bt_prop <- prop_stat(ped_bt)
rd_prop <- prop_stat(ped_RD)
ssh_prop <- prop_stat(ped_SSH)
sfw_prop <- prop_stat(ped_SFW)
sv_prop <- prop_stat(ped_SV)

sum(bt_prop)
sum(rd_prop)
sum(ssh_prop)
sum(sfw_prop)
sum(sv_prop)


cbind(cousins(ped_bt),
cousins(ped_RD),
cousins(ped_SFW),
cousins(ped_SV))

round(cbind(bt_prop,
rd_prop,
sfw_prop,
sv_prop),4)

barplot(rbind(bt_prop,
rd_prop,
sfw_prop,
sv_prop), beside=TRUE)


bt<-pedigreeStats(ped_bt, includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE)


pedStatSummary(bt)
# maternities is mother-offspring links, but not unique mothers

## total links
N<-pedStatSummary(bt)[1] * (pedStatSummary(bt)[1] - 1) / 2

round(pedStatSummary(bt)[2:12]/N,6)*100


source("/Users/joelpick/github/squidPed/R/simulate_pedigree.R")

generations=5
n_females=100
fecundity=4

fs_ped <-	simulate_pedigree(
		years = generations,
		n_females = n_females,
		fecundity = fecundity,
		p_sire = 1, 				# mating system (0-1, 1= one male per female, 0=complete random mating)
		juv_surv = 2/fecundity, # insures no population growth
		adult_surv = 0,					# discrete generations
		immigration = 0, 				# closed population
		constant_pop = TRUE     # constant population size
		)$pedigree

hs_ped <-	simulate_pedigree(
		years = generations,
		n_females = n_females,
		fecundity = fecundity,
		p_sire = 0, 				# mating system (0-1, 1= one male per female, 0=complete random mating)
		juv_surv = 2/fecundity, # insures no population growth
		adult_surv = 0,					# discrete generations
		immigration = 0, 				# closed population
		constant_pop = TRUE     # constant population size
		)$pedigree

fhs_ped <-	simulate_pedigree(
		years = generations,
		n_females = n_females,
		fecundity = fecundity,
		p_sire = 0.75, 				# mating system (0-1, 1= one male per female, 0=complete random mating)
		juv_surv = 2/fecundity, # insures no population growth
		adult_surv = 0,					# discrete generations
		immigration = 0, 				# closed population
		constant_pop = TRUE     # constant population size
		)$pedigree

fhs10_ped <- simulate_pedigree(
		years = generations,
		n_females = n_females,
		fecundity = 10,
		p_sire = 0.75, 				# mating system (0-1, 1= one male per female, 0=complete random mating)
		juv_surv = 0.2, # insures no population growth
		adult_surv = 0,					# discrete generations
		immigration = 0, 				# closed population
		constant_pop = TRUE     # constant population size
		)$pedigree



fs_stat<-pedigreeStats(fs_ped[,1:3], includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE)
hs_stat<-pedigreeStats(hs_ped[,1:3], includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE)
fhs_stat<-pedigreeStats(fhs_ped[,1:3], includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE)
fhs10_stat<-pedigreeStats(fhs10_ped[,1:3], includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE)


# pedStatSummary(fs_stat)

## total links
total_links<-function(ped) nrow(ped) * (nrow(ped) - 1) / 2

cbind(
	round(pedStatSummary(fs_stat)[2:12]/total_links(fs_ped),6)*100,
	round(pedStatSummary(hs_stat)[2:12]/total_links(hs_ped),6)*100,
	round(pedStatSummary(fhs_stat)[2:12]/total_links(fhs_ped),6)*100,
	round(pedStatSummary(fhs10_stat)[2:12]/total_links(fhs10_ped),6)*100
)

cbind(
	pedStatSummary(fs_stat)[2:12],
	pedStatSummary(hs_stat)[2:12],
	pedStatSummary(fhs_stat)[2:12],
	pedStatSummary(fhs10_stat)[2:12]
)

# make fhs10 samller so total sample size is smaller 
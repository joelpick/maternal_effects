rm(list=ls())
wd <- "/Users/joelpick/github/maternal_effects/"

source(paste0(wd,"R/extract_cousins.R"))
# devtools::install_github("cran/pedantics")
library(pedantics)




# Difference between blue tit and red deer pedigree

ped_bt <-  MasterBayes::insertPed(read.csv(paste0(wd,"Data/Raw/ped_BT.csv"))[,1:3])

head(ped_bt)

ped_RD<-read.csv(paste0(wd,"Data/Raw/ped_RD.csv"))
ped_SSH<-read.csv(paste0(wd,"Data/Raw/ped_SSH.csv"))
ped_SFW<-read.csv(paste0(wd,"Data/Raw/ped_SFW.csv"))
ped_SV<-read.csv(paste0(wd,"Data/Raw/ped_SV.csv"))

po(ped_RD)
gp(ped_RD)
pedSum(ped_RD)

pp<- sample (ped_RD[,1],1500)

cbind(ped_stat(ped_RD),ped_stat(ped_RD,pp))


cbind(
	ped_stat(ped_RD),
	ped_stat(ped_bt),
	ped_stat(ped_SSH),
	ped_stat(ped_SFW),
	ped_stat(ped_SV)
	)

ml <- cbind(
	RD=mat_links(ped_RD),
	SS=mat_links(ped_SSH),
	SV=mat_links(ped_SV),
	SFW=mat_links(ped_SFW),
	BT=mat_links(ped_bt)
	)

cov <- read.csv(paste0(wd,"Data/Raw/covariances.csv"))

rd_stat<-ped_stat(ped_RD)

biases <- function(ped,cov){
	rownames(cov)<-cov$relationship
	non_me <- c("dam","sire","PHS","MG","PG","au_D_FS","au_S_FS","au_D_MHS","au_S_MHS","au_D_PHS","au_S_PHS","cousin_D_FS","cousin_DS_FS","cousin_S_FS","cousin_D_HS","cousin_DS_HS","cousin_S_HS")
	
	stat<-ped_stat(ped)

	c(Vmg = sum(stat[non_me] * cov[non_me,"Vmg"]) / sum(stat[non_me]),
		COVam = sum(stat[non_me] * cov[non_me,"COVam"]) / sum(stat[non_me]))

}
## or is it to do with how much Va is as well - so need to multiply by 1/Va cov

biases(ped_bt,cov)
biases(ped_RD,cov)
biases(ped_SSH,cov)
biases(ped_SV,cov)
biases(ped_SFW,cov)


sum(rd_stat[match(cov$relationship,names(rd_stat))] * cov[,"Vmg"])/ sum(rd_stat[match(cov$relationship,names(rd_stat))])

sum(rd_stat[match(cov$relationship,names(rd_stat))] * cov[,"COVam"])/ sum(rd_stat[match(cov$relationship,names(rd_stat))])



barplot(apply(ml,2,function(x) x/sum(x)))


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
		n_females = n_females/2,
		fecundity = fecundity*2,
		p_sire = 0.75, 				# mating system (0-1, 1= one male per female, 0=complete random mating)
		juv_surv = 0.25, # insures no population growth
		adult_surv = 0,					# discrete generations
		immigration = 0, 				# closed population
		constant_pop = TRUE     # constant population size
		)$pedigree

fs10_ped <-	simulate_pedigree(
		years = generations,
		n_females = n_females/2,
		fecundity = fecundity*2,
		p_sire = 1, 				# mating system (0-1, 1= one male per female, 0=complete random mating)
		juv_surv = 0.25, # insures no population growth
		adult_surv = 0,					# discrete generations
		immigration = 0, 				# closed population
		constant_pop = TRUE     # constant population size
		)$pedigree

fs_ped_I <-	simulate_pedigree(
		years = generations,
		n_females = n_females,
		fecundity = fecundity,
		p_sire = 1, 				# mating system (0-1, 1= one male per female, 0=complete random mating)
		juv_surv = 1.5/fecundity, # insures no population growth
		adult_surv = 0,					# discrete generations
		immigration = 0.25, 				# closed population
		constant_pop = TRUE     # constant population size
		)$pedigree

fs_ped_I2 <-	simulate_pedigree(
		years = generations,
		n_females = n_females,
		fecundity = fecundity,
		p_sire = 1, 				# mating system (0-1, 1= one male per female, 0=complete random mating)
		juv_surv = c(1.5/fecundity,2/fecundity), # insures no population growth
		adult_surv = 0,					# discrete generations
		immigration = c(0.25,0), 				# closed population
		constant_pop = TRUE     # constant population size
		)$pedigree

cbind(
	mat_links(fs_ped),
	mat_links(fs10_ped),
	mat_links(hs_ped),
	mat_links(fhs_ped),
	mat_links(fhs10_ped),
	mat_links(fs_ped_I),
	mat_links(fs_ped_I2)
)


nrow(fs_ped)
nrow(hs_ped)
nrow(fhs_ped)
## difference in total size of pedigree due to starting number of females - shouldn't be a difference in phenotyped individuals
nrow(fhs10_ped)
nrow(fs10_ped)
nrow(fs_ped_I) ## different in total size of pedigree because there is immigrants and the offspring that leave - shouldn't be a difference in phenotyped individuals
nrow(fs_ped_I2)

cbind(cousins(fs_ped),
cousins(hs_ped),
cousins(fhs_ped),
cousins(fhs10_ped))

barplot(rbind(prop_stat(fs_ped[,1:3]),
prop_stat(hs_ped[,1:3]),
prop_stat(fhs_ped[,1:3]),
prop_stat(fhs10_ped[,1:3])),beside=TRUE)

##why no cousin_FS_FM in FS design?? there are a_u_FS_FM

library(pedantics)
fs_stat<-pedigreeStats(fs_ped[,1:3], includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE)
pedStatSummary(fs_stat)[2:12]
sum(!is.na(fs_ped[,2]))
sum(!is.na(fs_ped[,3]))
gp(fs_ped)


hs_stat<-pedigreeStats(hs_ped[,1:3], includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE)
fhs_stat<-pedigreeStats(fhs_ped[,1:3], includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE)
fhs10_stat<-pedigreeStats(fhs10_ped[,1:3], includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE)
fs10_stat<-pedigreeStats(fs10_ped[,1:3], includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE)
fs_stat_I<-pedigreeStats(fs_ped_I[,1:3], includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE)
fs_stat_I2<-pedigreeStats(fs_ped_I2[,1:3], includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE)


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
	pedStatSummary(fs10_stat)[2:12],
	pedStatSummary(hs_stat)[2:12],
	pedStatSummary(fhs_stat)[2:12],
	pedStatSummary(fhs10_stat)[2:12],
	pedStatSummary(fs_stat_I)[2:12],
	pedStatSummary(fs_stat_I2)[2:12]
)

cbind(cousins(fs_ped),
	cousins(fs10_ped),
cousins(hs_ped),
cousins(fhs_ped),
cousins(fhs10_ped),
cousins(fs_ped_I),cousins(fs_ped_I2)
)
1

# make fhs10 samller so total sample size is smaller 
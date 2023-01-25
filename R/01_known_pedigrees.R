rm(list=ls())
wd <- "/Users/joelpick/github/maternal_effects/"

source(paste0(wd,"R/extract_cousins.R"))
source(paste0(wd,"R/00_functions.R"))

# ped_bt <-  MasterBayes::insertPed(read.csv(paste0(wd,"Data/Raw/ped_BT.csv"))[,1:3])

dat_bt<-subset(read.csv("/Users/joelpick/Dropbox/0_blue_tits/skew/Data/Raw/tPED.csv"), ped_type=="DSCM")[,1:4]
ped_bt<-pedantics::fixPedigree(dat_bt[,1:3])
ped_bt[,"cohort"] <- dat_bt[match(ped_bt[,1],dat_bt[,1]),4]
ped_bt[,4] <- as.numeric(substr(ped_bt[,4],1,2))
names(ped_bt)[1] <- "animal"

head(ped_bt)

# ped_rd<-read.csv(paste0(wd,"Data/Raw/ped_RD.csv"))
ped_rd<-read.table("/Users/joelpick/github/pedigree_simulations/Data/red_deer/Pedigree_File.txt", header=TRUE)
names(ped_rd) <- c("animal","sire","dam")
dat_rd<-read.table("/Users/joelpick/github/pedigree_simulations/Data/red_deer/Individual_Data.txt", header=TRUE)

ped_rd$cohort <- dat_rd[match(ped_rd[,"ID"],dat_rd[,"ID"]),"BirthYear"]



par(mfrow=c(2,1))
plot(table(table(ped_bt$dam)))
mean(table(ped_bt$dam))

plot(table(table(ped_rd$dam)))
mean(table(ped_rd$dam))

ped<-ped_bt
immigration <- function(ped, sex_specific=TRUE){
	founders <- ped$animal[is.na(ped$sire) & is.na(ped$sire)]
	ped$dam_founder <- ped$dam%in%founders
	ped$sire_founder <- ped$sire%in%founders
	year_d <- aggregate(cohort~dam,subset(ped,!dam%in%founders),min)$cohort
	year_d_F <- aggregate(cohort~dam,subset(ped,dam%in%founders),min)$cohort
	year_s <- aggregate(cohort~sire,subset(ped,!sire%in%founders),min)$cohort
	year_s_F <- aggregate(cohort~sire,subset(ped,sire%in%founders),min)$cohort

	year_d <- aggregate(cohort~dam+dam_founder,ped,min)
	year_s <- aggregate(cohort~sire+sire_founder,ped,min)

	all_d <- table(year_d$cohort,year_d$dam_founder)
	all_s <- table(year_s$cohort,year_s$sire_founder)

	all_s <- table(year_s)
	all_d_F <- table(year_d_F)
	all_s_F <- table(year_s_F)


	print(all_d)
	print(all_s)

	c(
		f=mean(all_d[-1]/all_d[1]),
		m=mean(all_s[-1]/all_s[1])
		)
}

immigration(ped_bt)

plot(as.numeric(table(ped_rd$cohort)))

plot(aggregate(dam~cohort,ped_rd,function(x) length(unique(x)))$dam)
plot(aggregate(sire~cohort,ped_rd,function(x) length(unique(x)))$sire)






rm(list=ls())


library(asreml)
library(parallel)

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/extract_cousins.R"))
source(paste0(wd,"R/00_functions.R"))
source("/Users/joelpick/github/squidPed/R/simulate_pedigree.R")
# devtools::load_all("~/github/squidSim/R")

run=FALSE

n_sims <-100
generations=5
n_females=100
fecundity=10

n_females=n_females
fecundity=fecundity
p_sire = 1
p_polyandry=0.5
juv_surv = rep(2/fecundity,2)
immigration = c(0,0)


ped <- simulate_pedigree(
	years = generations,
	n_females = n_females,
	fecundity = fecundity,
	p_sire = p_sire,
	p_polyandry=1,
	juv_surv = juv_surv,
	adult_surv = 0,					# discrete generations
	immigration = immigration, 				# closed population
	constant_pop = TRUE     # constant population size
)$pedigree



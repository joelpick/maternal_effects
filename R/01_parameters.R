
######
#  This script defines the parameters used in the simulations
######

rm(list=ls())

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

##------------------------
## Pedigree parameters
##------------------------

## total offspring per generation
n_offspring=600

#Fecundity values
## Chosen to be comparable to range from bonnet et al 
m_fecundity = 6 # median is 6.5
h_fecundity = 12 #highest is 13.3
l_fecundity = 3 #lowest is 1.5
## these numbers are easier to maintain total number of phenotypes individuals

## first number is female, then male immigration rate
no_immigration = c(0,0)
f_immigration=c(0.4,0.1)   # female bias
m_immigration=c(0.1,0.4)   # male bias
b_immigration=c(0.25,0.25) # same overall immigration but no sex bias

# mating system parameters (full versus half sib)
fhs = 0.75
fs = 1
hs = 0

## put these together
fecundity <- c(m_fecundity,h_fecundity,l_fecundity)
immigration <- rbind(no_immigration,f_immigration,m_immigration,b_immigration)
p_sire <- c(fhs,hs,fs)

f_names <- c("mF","hF","lF")
i_names <- c("nI","fI","mI","bI")
ms_names <- c("fhs","hs","fs")

## make all combos 
combos<-expand.grid(p_sire=1:3,immigration=1:4,fecundity=1:3)

## make parameters for pedigree simulations
peds_param <- cbind(
	n_females=n_offspring / fecundity[combos[,"fecundity"]],
	fecundity=fecundity[combos[,"fecundity"]],
	p_sire = p_sire[combos[,"p_sire"]],
	immigration_f = immigration[combos[,"immigration"],1],
	immigration_m = immigration[combos[,"immigration"],2])

## work out juvenile survival, so that populations are stationary
peds_param <- cbind(peds_param, 
	juv_surv_f=2*(1 - peds_param[,"immigration_f"])/peds_param[,"fecundity"],
	juv_surv_m=2*(1 - peds_param[,"immigration_m"])/peds_param[,"fecundity"]
)

## generate pedigree names
ped_names <- rownames(peds_param) <- paste(
	ms_names[combos[,"p_sire"]],
	f_names[combos[,"fecundity"]],
	i_names[combos[,"immigration"]], sep="_")




##------------------------
## Reduced set of pedigree parameters
##------------------------

# Postma 2014, the median sample size of studies using animal models is 361
# in young and postma 2023 its 420
    #  0%     10%     25%     50%     75%     90%    100% 
    # 6.0   105.4   174.5   420.0  1106.0  1895.0 38024.0 

peds_param_reduced <- rbind(
	cbind(generations=2,peds_param[grep("fhs_lF",ped_names),]),
	cbind(generations=4,peds_param[grep("fhs_lF",ped_names),])
	)
peds_param_reduced[,"n_females"] <- rep(c(20,30),each=4)

ped_names_reduced <- rownames(peds_param_reduced) <- c(
	paste0(ped_names[grep("fhs_lF",ped_names)],"_small"),
	paste0(ped_names[grep("fhs_lF",ped_names)],"_medium")
	)



##-------------------------------
## Different simulation scenarios
##-------------------------------

scenarios <- rbind(	
	# A) High Maternal genetic only
	a=c(Va=0, Vmg=0.5, Vme=0, r_amg=0),
	# B) Direct genetic and maternal environment
	b=c(Va=0, Vmg=0.25, Vme=0.25, r_amg=0),
	# C) High Maternal environment only
	c=c(Va=0, Vmg=0, Vme=0.5, r_amg=0),

	# D) Maternal genetic only
	d=c(Va=0, Vmg=0.25, Vme=0, r_amg=0),
	# E) Direct, maternal genetic and maternal environment
	e=c(Va=0.25, Vmg=0.25, Vme=0.25, r_amg=0),

	# F) Direct and maternal genetic, no covariance
	f=c(Va=0.25, Vmg=0.25, Vme=0, r_amg=0),
	# G+H) Direct and maternal genetic, + covariance
	g=c(Va=0.25, Vmg=0.25, Vme=0, r_amg=0.3),
	h=c(Va=0.25, Vmg=0.25, Vme=0, r_amg=0.6),
	# I+J) Direct and maternal genetic, - covariance
	i=c(Va=0.25, Vmg=0.25, Vme=0, r_amg=-0.3),
	j=c(Va=0.25, Vmg=0.25, Vme=0, r_amg=-0.6),

	# K+L) just Va
	k= c(Va=0.25, Vmg=0, Vme=0, r_amg=0),
	l= c(Va=0.25, Vmg=0, Vme=0.25, r_amg=0)	
)




save(ped_names,peds_param,ped_names_reduced,peds_param_reduced,scenarios,file=paste0(data_wd,"parameters.Rdata"))







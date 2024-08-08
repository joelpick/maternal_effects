
rm(list=ls())

library(asreml)
library(parallel)
library(squidPed)

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

# source(paste0(wd,"R/extract_cousins.R"))
source(paste0(wd,"R/00_functions.R"))
# source("/Users/joelpick/github/squidPed/R/simulate_pedigree.R")
# devtools::load_all("~/github/squidSim/R")

run=TRUE

n_sims <-100

##------------------------
## Pedigree parameters
##------------------------

generations=5
n_offspring=600

## Range from bonnet et al 
m_fecundity = 6 # mean is 6.9
h_fecundity = 12 #highest is 13.3
l_fecundity = 3 #lowest is 1.5
## these numbers are easier to maintain total number of phenotypes individuals

## first number is female, then male immigration rate
no_immigration = c(0,0)
f_immigration=c(0.4,0.1)
m_immigration=c(0.1,0.4)
b_immigration=c(0.25,0.25) # same overall immigration but no sex bias

fhs = 0.75
fs = 1
hs = 0

fecundity <- c(m_fecundity,h_fecundity,l_fecundity)
immigration <- rbind(no_immigration,f_immigration,m_immigration,b_immigration)
p_sire <- c(fhs,hs,fs)

f_names <- c("mF","hF","lF")
i_names <- c("nI","fI","mI","bI")
ms_names <- c("fhs","hs","fs")

## make all combos 
combos<-expand.grid(p_sire=1:3,immigration=1:4,fecundity=1:3)

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
	j=c(Va=0.25, Vmg=0.25, Vme=0, r_amg=-0.6)		
)


###------------------------
### Pedigree and Phenotype Simulations, and Running Models
###------------------------

if(run){
	cores<-6
	set.seed(20230126)

	ped_str <- vector("list",length=nrow(peds_param))
	names(ped_str) <- ped_names
	# ped_str_mat <- vector("list",length=nrow(peds_param))
	# names(ped_str_mat) <- ped_names

# k="fs_lF_mI"
	## make pedigrees
	
	for(k in ped_names){
		cat(k, "\n")
		cat("Simulating Pedigrees\n")
		peds <- mclapply(1:n_sims,	function(i){
			simulate_pedigree(
				years = generations,
				n_females = peds_param[k,"n_females"],
				fecundity = peds_param[k,"fecundity"],
				p_sire = peds_param[k,"p_sire"],
				p_polyandry=1,
				juv_surv = c(peds_param[k,"juv_surv_f"],peds_param[k,"juv_surv_m"]),
				adult_surv = 0,					# discrete generations
				immigration = c(peds_param[k,"immigration_f"],peds_param[k,"immigration_m"]), 				# closed population
				constant_pop = TRUE     # constant population size
				)$pedigree
		}, mc.cores=cores)
		# assign(paste0(k ,"_peds"),peds)	

		cat("Generating Pedigree Metrics\n")
		ped_str[[k]]<- do.call(rbind,mclapply(peds,FUN= function(x){
			ped_stat(x,phenotyped=x[!is.na(x[,"dam"]),"animal"])}, mc.cores=cores))
		# ped_str_mat[[k]]<- do.call(rbind,mclapply(peds,ped_stat2, mc.cores=cores))


		cat("Simulating Data\n")
		## simulate data
		dat<-mclapply(peds, function(i){
			x<-vector("list", nrow(scenarios))
			for(j in 1:nrow(scenarios)){
				x[[j]]<- mge_sim(i[,1:3], param=scenarios[j,])
			}
			x
		}, mc.cores=cores)
		# assign(paste0(k,"_data"),dat)

	## run models
		cat("Running models: \n")
	
		# cat("Model 1: ")
		# model1 <- model_func(m1a_func,peds,dat,mc.cores=cores)
		# assign(paste0("model1_",k),model1)
		cat("\nModel 2: ")
		model2 <- model_func(m2_func,peds,dat,mc.cores=cores)
		assign(paste0("model2_",k),model2)
		cat("\n")
		rm(peds,dat)
	}

#paste0("model1_",ped_names),
	save(list=(c("ped_names","peds_param","scenarios","ped_str",paste0("model2_",ped_names))),file=paste0(data_wd,"mge_sims3.Rdata"))#"ped_str_mat",
}else{
	load(paste0(data_wd,"mge_sims3.Rdata"))	
}



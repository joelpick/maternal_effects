
rm(list=ls())

library(asreml)
library(parallel)
library(squidPed)
library(squidSim)
# devtools::install_github("squid-group/squidSim")

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))
load(paste0(data_wd,"parameters.Rdata"))

run <- TRUE
cores <-6
n_sims <- 100

generations <- 5


###------------------------
### Pedigree and Phenotype Simulations, and Running Models
###------------------------

if(run){
	
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
				immigration = c(peds_param[k,"immigration_f"],peds_param[k,"immigration_m"]), 	
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
}



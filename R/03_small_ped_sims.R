

rm(list=ls())

# devtools::install_github("squidgroup/squidSim")
# devtools::install_github("squidgroup/pedAgree")

library(asreml)
library(parallel)
library(squidSim)
library(pedAgree)

wd <- "/Users/joelpick/github/maternal_effects/"
data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))
load(paste0(data_wd,"parameters.Rdata"))


cores<-6
n_sims <-100


###------------------------
### Pedigree and Phenotype Simulations, and Running Models
###------------------------
	
set.seed(20240822)

	

for(k in ped_names_reduced){
	cat(k, "\n")
	cat("Simulating Pedigrees\n")
	## make pedigrees
	peds <- lapply(1:n_sims,	function(i){
		ped<-simulate_pedigree(
			years = peds_param_reduced[k,"generations"],
			n_females = peds_param_reduced[k,"n_females"],
			fecundity = peds_param_reduced[k,"fecundity"],
			fixed_fecundity = TRUE,
			p_sire = peds_param_reduced[k,"p_sire"],
			p_polyandry=1,
			p_breed=1,
			juv_surv = c(peds_param_reduced[k,"juv_surv_f"],peds_param_reduced[k,"juv_surv_m"]),
			adult_surv = 0,					# discrete generations
			immigration = c(peds_param_reduced[k,"immigration_f"],peds_param_reduced[k,"immigration_m"]),
			constant_pop = TRUE   # constant population size
			)$pedigree
	})

	cat("Simulating Data\n")
	## simulate data
	dat<-mclapply(peds, function(i){
		x<-vector("list", nrow(scenarios))
		for(j in 1:nrow(scenarios)){
			x[[j]]<- mge_sim(i[,1:3], param=scenarios[j,])
		}
		x
	}, mc.cores=cores)

## run models
	cat("Running models: \n")

	cat("\nModel 1: ")
	model1 <- model_func(m1_func,peds,dat,mc.cores=cores)

	cat("\nModel 2: ")
	model2 <- model_func(m2_func,peds,dat,mc.cores=cores)

	cat("\nModel 4: ")
	model4 <- model_func(m4_func,peds,dat,mc.cores=cores)

	cat("\nModel 5: ")
	model5 <- model_func(m5_func,peds,dat,mc.cores=cores)

	assign(paste0("model1_",k),model1)
	assign(paste0("model2_",k),model2)
	assign(paste0("model4_",k),model4)
	assign(paste0("model5_",k),model5)
	
	cat("\n")
	rm(peds,dat)
}


save(list=c(
	paste0("model1_",ped_names_reduced),
	paste0("model2_",ped_names_reduced),
	paste0("model4_",ped_names_reduced),
	paste0("model5_",ped_names_reduced)
	),
	file=paste0(data_wd,"mge_sims_small_ped.Rdata"))

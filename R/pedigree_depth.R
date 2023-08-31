
rm(list=ls())

library(parallel)

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source("/Users/joelpick/github/squidPed/R/simulate_pedigree.R")

max_gen <- 10
n_sims <- 100
cores<-6
	ped_str <- vector("list",length=max_gen)
	
	for(k in 1:max_gen){
		cat(k, "\n")
		ped_str[[k]] <- do.call(rbind,mclapply(1:n_sims,	function(i){
			ped <- simulate_pedigree(
				years = k,
				n_females = 100,
				fecundity = 8,
				p_sire = 0.75,
				p_polyandry=1,
				juv_surv = 0.25,
				adult_surv = 0,					# discrete generations
				immigration = 0, 				# closed population
				constant_pop = TRUE     # constant population size
				)$pedigree
			ped_stat(ped)
		}, mc.cores=cores))		
	}
x<-ped_str[[1]]
	mat_ratio_all<-sapply(ped_str,function(x){
	rowSums(x[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(x[,-(1:2)]) 
})
	mat_ratio <- colMeans(mat_ratio_all)

plot(mat_ratio, pch=19, ylim=c(0,0.3), xlab="Pedigree Depth (Generations)", ylab="Proportion of non-sibling maternal links")

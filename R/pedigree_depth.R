
rm(list=ls())

library(parallel)

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source("/Users/joelpick/github/squidPed/R/simulate_pedigree.R")

max_gen <- 10
n_sims <- 50
cores<-6
fecundity=8

immigration <- rbind(
	no_immigration = c(0,0),
	f_immigration=c(0.4,0.1),
	m_immigration=c(0.1,0.4),
	b_immigration=c(0.25,0.25))

combos<-expand.grid(gens=1:10,immigration=1:4)

peds_param <- cbind(
	years= combos[,"gens"],
	immigration_f = immigration[combos[,"immigration"],1],
	immigration_m = immigration[combos[,"immigration"],2])

peds_param <- cbind(peds_param, 
	juv_surv_f=2*(1 - peds_param[,"immigration_f"])/fecundity,
	juv_surv_m=2*(1 - peds_param[,"immigration_m"])/fecundity
)

# ped_names <- rownames(peds_param) <- paste(
# 	ms_names[combos[,"p_sire"]],
# 	f_names[combos[,"fecundity"]],
# 	i_names[combos[,"immigration"]], sep="_")


	ped_str <- vector("list",length=max_gen)
	
	for(k in 1:length(peds_param)){
		cat(k, "\n")
		ped_str[[k]] <- do.call(rbind,mclapply(1:n_sims,	function(i){
			ped <- simulate_pedigree(
				years = k,
				n_females = 100,
				fecundity = fecundity,
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

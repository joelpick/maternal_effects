
rm(list=ls())

source("/Users/joelpick/github/pedigree_simulations/R/simulate_pedigree.R")
# devtools::load_all("~/github/squidSim/R")

library(asreml)
library(parallel)

wd <- 

source("/Users/joelpick/github/pedigree_simulations/R/simulate_pedigree.R")

 
n_sims <-100
generations=5
n_females=100
fecundity=4


#### are number of maternal links the same in the two pedigrees?
#### number of maternal half sibs will differ presumably
#### moore et al meta-analysis mentions taking into account pedigree links in the supp matt

## create 100 full sib pedigrees 

fs_peds <- lapply(1:n_sims,	function(i){
	simulate_pedigree(
		years = generations,
		n_females = n_females,
		fecundity = fecundity,
		p_sire = 1, 				# mating system (0-1, 1= one male per female, 0=complete random mating)
		juv_surv = 2/fecundity, # insures no population growth
		adult_surv = 0,					# discrete generations
		immigration = 0, 				# closed population
		constant_pop = TRUE     # constant population size
		)$pedigree
})

hs_peds <- lapply(1:n_sims,	function(i){
	simulate_pedigree(
		years = generations,
		n_females = n_females,
		fecundity = fecundity,
		p_sire = 0, 				# mating system (0-1, 1= one male per female, 0=complete random mating)
		juv_surv = 2/fecundity, # insures no population growth
		adult_surv = 0,					# discrete generations
		immigration = 0, 				# closed population
		constant_pop = TRUE     # constant population size
		)$pedigree
})

scenarios <- rbind(
	a=c(Va=0.2, Vmg=0, r_amg=0, Vme=0),
		# B) Maternal environment only
		b=c(Va=0, Vmg=0, r_amg=0, Vme=0.4),
		# C) Maternal genetic only
		c=c(Va=0, Vmg=0.4, r_amg=0, Vme=0),
		# D) Direct genetic and maternal environment
		d=c(Va=0.2, Vmg=0, r_amg=0, Vme=0.4),
		# E) Maternal genetic and maternal environment
		e=c(Va=0, Vmg=0.4, r_amg=0, Vme=0.2),#####
		# F) Direct and maternal genetic, no covariance
		f=c(Va=0.1, Vmg=0.2, r_amg=0, Vme=0),
		# G) Direct and maternal genetic, positive covariance
		g=c(Va=0.1, Vmg=0.2, r_amg=0.5, Vme=0),
		# H) Direct and maternal genetic, negative covariance
		h=c(Va=0.1, Vmg=0.2, r_amg=-0.5, Vme=0),
		# I) Direct and maternal genetic, no covariance and maternal environment
		i=c(Va=0.1, Vmg=0.2, r_amg=0, Vme=0.2),
		# J) Direct and maternal genetic, positive covariance and maternal environment
		j=c(Va=0.1, Vmg=0.2, r_amg=0.5, Vme=0.2),
		# K) Direct and maternal genetic, negative covariance and maternal environment
		k=c(Va=0.1, Vmg=0.2, r_amg=-0.5, Vme=0.2)
	)

fs_data <- lapply(fs_peds, function(i){
	list(
		# A) Direct genetic effects only
		a=mge_sim(i[,1:3], Va=0.2, Vmg=0, r_amg=0, Vme=0),
		# B) Maternal environment only
		b=mge_sim(i[,1:3], Va=0, Vmg=0, r_amg=0, Vme=0.4),
		# C) Maternal genetic only
		c=mge_sim(i[,1:3], Va=0, Vmg=0.4, r_amg=0, Vme=0),
		# D) Direct genetic and maternal environment
		d=mge_sim(i[,1:3], Va=0.2, Vmg=0, r_amg=0, Vme=0.4),
		# E) Maternal genetic and maternal environment
		e=mge_sim(i[,1:3], Va=0, Vmg=0.4, r_amg=0, Vme=0.2),#####
		# F) Direct and maternal genetic, no covariance
		f=mge_sim(i[,1:3], Va=0.1, Vmg=0.2, r_amg=0, Vme=0),
		# G) Direct and maternal genetic, positive covariance
		g=mge_sim(i[,1:3], Va=0.1, Vmg=0.2, r_amg=0.5, Vme=0),
		# H) Direct and maternal genetic, negative covariance
		h=mge_sim(i[,1:3], Va=0.1, Vmg=0.2, r_amg=-0.5, Vme=0),
		# I) Direct and maternal genetic, no covariance and maternal environment
		i=mge_sim(i[,1:3], Va=0.1, Vmg=0.2, r_amg=0, Vme=0.2),
		# J) Direct and maternal genetic, positive covariance and maternal environment
		j=mge_sim(i[,1:3], Va=0.1, Vmg=0.2, r_amg=0.5, Vme=0.2),
		# K) Direct and maternal genetic, negative covariance and maternal environment
		k=mge_sim(i[,1:3], Va=0.1, Vmg=0.2, r_amg=-0.5, Vme=0.2)
	)
})

hs_data <- lapply(hs_peds, function(i){
	list(
		# A) Direct genetic effects only
		a=mge_sim(i[,1:3], Va=0.2, Vmg=0, r_amg=0, Vme=0),
		# B) Maternal environment only
		b=mge_sim(i[,1:3], Va=0, Vmg=0, r_amg=0, Vme=0.4),
		# C) Maternal genetic only
		c=mge_sim(i[,1:3], Va=0, Vmg=0.4, r_amg=0, Vme=0),
		# D) Direct genetic and maternal environment
		d=mge_sim(i[,1:3], Va=0.2, Vmg=0, r_amg=0, Vme=0.4),
		# E) Maternal genetic and maternal environment
		e=mge_sim(i[,1:3], Va=0, Vmg=0.4, r_amg=0, Vme=0.2),#####
		# F) Direct and maternal genetic, no covariance
		f=mge_sim(i[,1:3], Va=0.1, Vmg=0.2, r_amg=0, Vme=0),
		# G) Direct and maternal genetic, positive covariance
		g=mge_sim(i[,1:3], Va=0.1, Vmg=0.2, r_amg=0.5, Vme=0),
		# H) Direct and maternal genetic, negative covariance
		h=mge_sim(i[,1:3], Va=0.1, Vmg=0.2, r_amg=-0.5, Vme=0),
		# I) Direct and maternal genetic, no covariance and maternal environment
		i=mge_sim(i[,1:3], Va=0.1, Vmg=0.2, r_amg=0, Vme=0.2),
		# J) Direct and maternal genetic, positive covariance and maternal environment
		j=mge_sim(i[,1:3], Va=0.1, Vmg=0.2, r_amg=0.5, Vme=0.2),
		# K) Direct and maternal genetic, negative covariance and maternal environment
		k=mge_sim(i[,1:3], Va=0.1, Vmg=0.2, r_amg=-0.5, Vme=0.2)
	)
})

save(scenarios,fs_peds,hs_peds,fs_data,hs_data, file="/Users/joelpick/Dropbox/0_fitness/fitness_timing/Data/Intermediate/gaussian_data.Rdata")


rm(list=ls())

source("/Users/joelpick/github/pedigree_simulations/R/simulate_pedigree.R")
# devtools::load_all("~/github/squidSim/R")

library(asreml)
library(parallel)

rbv0 <- function(pedigree, G){
	X <- matrix(0, nrow=nrow(pedigree), ncol=nrow(G))
	index <- which(diag(G)!=0)
	if(any(diag(G)==0)) G <- G[index,index]
	X2 <- MCMCglmm::rbv(pedigree=pedigree, G=G)
	X[,index] <- X2
	X
}
 
 generations=5
 n_females=100
 fecundity=4
 p_sire=0 
 models=1:8
 Va=0.5
 Vmg=0
 r_amg=0
 Vme=0
 Ve=0.5


 ## function to simulate pedigree, phenotypes from the pedigree and analyse with different models
mge_sim <- function(generations=5, n_females=100, fecundity=4, p_sire=0, Va, Vmg,r_amg, Vme, Vp=1, models=1:8){

## create pedigree according to specified 
	ped <- simulate_pedigree(
		years = generations,
		n_females = n_females,
		fecundity = fecundity,
		p_sire = p_sire, 				# mating system (0-1, 1= one male per female, 0=complete random mating)
		juv_surv = 2/fecundity, # insures no population growth
		adult_surv = 0,					# discrete generations
		immigration = 0, 				# closed population
		constant_pop = TRUE     # constant population size
		)$pedigree

	## make genetic covariance matrix between direct and maternal genetic effects
	D <- diag(sqrt(c(Va,Vmg)))
	R <- matrix(c(1,r_amg,r_amg,1),nrow=2)
	G <- D%*%R%*%D

	Ve <- Vp - sum(G) - Vme

	## simulate direct, maternal genetic and environmental effects, add together to make phenotype
	g <- rbv0(ped[,1:3],G)
	a <- g[,1]
	mg <- g[match(ped[,2],ped[,1]),2]
	me <- rnorm(nrow(ped),0,sqrt(Vme))[match(ped[,2],ped[,1])]
	e <- rnorm(nrow(ped),0,sqrt(Ve))
	p <- a + mg + me + e

	## prepare data for models
	data <- data.frame(cbind(p=p,ped))
	data$mother <- as.factor(data$dam)
	data$mother_PE <- as.factor(data$dam)
	data$animal <- as.factor(data$animal)
	data <- subset(data, !is.na(mother)) 
	#individual without a know mother are excluded, as they cant be used  for estimating Va assuming phenotype of mother 

	## inverse A matrix for models
	assign("ped.ainv", asreml::ainverse(ped), envir = .GlobalEnv) 

## asreml.options(workspace="256mb")

  ## Model 1 - additive genetic effects, assuming phenotype is trait of offspring
	if(1 %in% models){
		mod_1 <- asreml(fixed= p~1
		                   , random= ~vm(animal,ped.ainv)
		                   , data= data
		                   , trace=FALSE
		                 	)
		m1 <- summary(mod_1)$varcomp
		m1_sum <- c(
			A= m1["vm(animal, ped.ainv)",1], 
			Me=NA, 
			Mg=NA,
			cov_AMg =NA,
			E = m1["units!R",1])
	}


  ## Model 2 - additive genetic effects with maternal environment effects, assuming phenotype is trait of offspring
	if(2 %in% models){
		mod_2 <- asreml(fixed= p~1
		                   , random= ~vm(animal,ped.ainv) + mother_PE
		                   , data= data
		                   , trace=FALSE
		                 	)
		m2<-summary(mod_2)$varcomp
		m2_sum <- c(
			A= m2["vm(animal, ped.ainv)",1], 
			Me=m2["mother_PE",1], 
			Mg=NA,
			cov_AMg =NA,
			E = m2["units!R",1])
	}

  ## Model 3 - additive and maternal genetic effects, assuming phenotype is trait of offspring
	if(3 %in% models){
		mod_3 <- asreml(fixed= p~1
	                   , random= ~vm(animal,ped.ainv) + vm(mother,ped.ainv) 
	                   , data= data
	                   , trace=FALSE
	                 	)
		m3<-summary(mod_3)$varcomp
		m3_sum <- c(
			A= m3["vm(animal, ped.ainv)",1], 
			Me=NA, 
			Mg=m3["vm(mother, ped.ainv)",1],
			cov_AMg =NA,
			E = m3["units!R",1])

	}

	## Model 4 - 
	if(4 %in% models){
		mod_4 <- asreml(fixed= p~1
		                   , random= ~vm(animal,ped.ainv) + vm(mother,ped.ainv) + mother_PE
		                   , data= data
											 , trace=FALSE
		                 	)
		m4<-summary(mod_4)$varcomp
		m4_sum <- c(
			A= m4["vm(animal, ped.ainv)",1], 
			Me=m4["mother_PE",1], 
			Mg=m4["vm(mother, ped.ainv)",1],
			cov_AMg =NA,
			E = m4["units!R",1])
	}

  ## Model 5 - additive and maternal genetic effects and covariance, assuming phenotype is trait of offspring
	if(5 %in% models){
		mod_5 <- asreml(fixed= p~1
	                   , random= ~str(~vm(animal,ped.ainv) +vm(mother,ped.ainv) ,~us(2):vm(animal,ped.ainv)) + mother_PE
	                   , data= data
	                   , trace=FALSE
	                   ,maxit=50
	                 	)
		m5<-summary(mod_5)$varcomp
		m5_sum <- c(
			A= m5[1,1], 
			Me=m5["mother_PE",1], 
			Mg=m5[3,1],
			cov_AMg =m5[2,1],
			E = m5["units!R",1])
	}

asreml.options(maxit=20 )


  ## Model 6 - additive genetic effects, assuming phenotype is trait of mother
	if(6 %in% models){
		mod_6 <- asreml(fixed= p~1
	                   , random= ~vm(mother,ped.ainv)
	                   , data= data
	                   , trace=FALSE
	                 	)
		m6<-summary(mod_6)$varcomp
		m6_sum <- c(
			A= NA, 
			Me=NA, 
			Mg=m6["vm(mother, ped.ainv)",1],
			cov_AMg =NA,
			E = m6["units!R",1])
	}

	## Model 7 - additive genetic effects assuming phenotype is trait of mother, with permanent (maternal) environment effects, 
	if(7 %in% models){
		mod_7 <- asreml(fixed= p~1
		                   , random= ~vm(mother,ped.ainv) + mother_PE
		                   , data= data
											 , trace=FALSE
		                 	)
		m7<-summary(mod_7)$varcomp
		m7_sum <- c(
			A= NA, 
			Me=m7["mother_PE",1], 
			Mg=m7["vm(mother, ped.ainv)",1],
			cov_AMg =NA,
			E = m7["units!R",1])
	}

 
  ## Model 8 - additive genetic effects, assuming phenotype is trait of mother, whilst taking the mean of the offspring phenotypes
	if(8 %in% models){
		data_means <- aggregate(p~mother,data,mean)
		mod_8 <- asreml(fixed= p~1
	                   , random= ~vm(mother,ped.ainv)
	                   , data= data_means
	                   , trace=FALSE
	                 	)
		m8<-summary(mod_8)$varcomp
		m8_sum <- c(
			A= NA, 
			Me=NA, 
			Mg=m8["vm(mother, ped.ainv)",1],
			cov_AMg =NA,
			E = m8["units!R",1])
	}

	rm("ped.ainv", envir = .GlobalEnv)

	rbind(
		m1=if(1 %in% models) {m1_sum}else{NA},
		m2=if(2 %in% models) {m2_sum}else{NA},
		m3=if(3 %in% models) {m3_sum}else{NA},
		m4=if(4 %in% models) {m4_sum}else{NA},
		m5=if(5 %in% models) {m5_sum}else{NA},
		m6=if(6 %in% models) {m6_sum}else{NA},
		m7=if(7 %in% models) {m7_sum}else{NA},
		m8=if(8 %in% models) {m8_sum}else{NA}
	)
}



# Scenario 1 - Direct genetic effects only
s1_hs <- mclapply(1:100,function(i){
	out<-mge_sim(Va=0.25, Vmg=0, r_amg=0, Vme=0, p_sire=0)
	cat(i," ")
	out
}, mc.cores=8)

s1_fs <- mclapply(1:100,function(i){
	out<-mge_sim(Va=0.25, Vmg=0, r_amg=0, Vme=0, p_sire=1)
	cat(i," ")
	out
}, mc.cores=8)


# Scenario 2 - Maternal environmental effects only
s2_hs <- mclapply(1:100,function(i){
	out<-mge_sim(Va=0, Vmg=0, r_amg=0, Vme=0.3, p_sire=0)
	cat(i," ")
	out
}, mc.cores=8)

s2_fs <- mclapply(1:100,function(i){
	out<-mge_sim(Va=0, Vmg=0, r_amg=0, Vme=0.3, p_sire=1)
	cat(i," ")
	out
}, mc.cores=8)



# Scenario 2 - Maternal genetic effects only
s2_hs <- mclapply(1:100,function(i){
	out<-mge_sim(Va=0, Vmg=0.5, r_amg=0, Vme=0, p_sire=0)
	cat(i," ")
	out
}, mc.cores=8)

s2_fs <- mclapply(1:100,function(i){
	out<-mge_sim(Va=0, Vmg=0.5, r_amg=0, Vme=0, p_sire=1)
	cat(i," ")
	out
}, mc.cores=8)


# Scenario 3 - Both Uncorrelated

s3_hs <- mclapply(1:100,function(i){
	out<-mge_sim(Va=0.25, Vmg=0.25, r_amg=0, Vme=0, p_sire=0)
	cat(i," ")
	out
}, mc.cores=8)

s3_fs <- mclapply(1:100,function(i){
	out<-mge_sim(Va=0.25, Vmg=0.25, r_amg=0, Vme=0, p_sire=1)
	cat(i," ")
	out
}, mc.cores=8)



# Scenario 4 - Both Correlated


s4_hs <- mclapply(1:100,function(i){
	out<-mge_sim(Va=0.2, Vmg=0.2, r_amg=0.25, Vme=0, p_sire=0)
	cat(i," ")
	out
}, mc.cores=8)

s4_fs <- mclapply(1:100,function(i){
	out<-mge_sim(Va=0.2, Vmg=0.2, r_amg=0.25, Vme=0, p_sire=1)
	cat(i," ")
	out
}, mc.cores=8)


save(s1_hs,s1_fs,s2_hs,s2_fs,s3_hs,s3_fs,s4_hs,s4_fs, file="/Users/joelpick/Dropbox/0_fitness/fitness_timing/Data/Intermediate/gaussian_sims.Rdata")








out<-mge_sim(generations=5, n_females=100, Va=0.5, Vmg=0, r_amg=0, Vme=0, p_sire=0)



## repeated measured, need maternal id as well as 
## how does this translate to fitness where there is one measure 

# Scenario 1 - Direct genetic effects only
s1_hs <- mclapply(1:100,function(i){
	out<-mge_sim(generations=20, n_females=100, Va=0.25, Vmg=0, r_amg=0, Vme=0,Ve=0.75, p_sire=0, models=c(3,4,7))
	cat(i," ")
	out
}, mc.cores=8)

s1_hs2 <- mclapply(1:100,function(i){
	out<-mge_sim(generations=2, n_females=1000, Va=0.5, Vmg=0, r_amg=0, Ve=0.5, p_sire=0, models=c(3,4,7))
	cat(i," ")
	out
}, mc.cores=8)

s1_fs <- mclapply(1:100,function(i){
	out<-mge_sim(generations=20, n_females=100, Va=0.5, Vmg=0, r_amg=0, Ve=0.5, p_sire=1, models=c(3,4,7))
	cat(i," ")
	out
}, mc.cores=8)
s1_fs2 <- mclapply(1:100,function(i){
	out<-mge_sim(generations=2, n_females=1000, Va=0.5, Vmg=0, r_amg=0, Ve=0.5, p_sire=1, models=c(3,4,7))
	cat(i," ")
	out
}, mc.cores=8)

s1_fs3 <- mclapply(1:100,function(i){
	out<-mge_sim(generations=5, n_females=100, Va=0.5, Vmg=0, r_amg=0, Ve=0.5, p_sire=1, models=c(3,4,7))
	cat(i," ")
	out
}, mc.cores=8)


rbind(
	colMeans(do.call(rbind,lapply(s1_hs,function(x) x[1,]))),
	colMeans(do.call(rbind,lapply(s1_hs,function(x) x[2,]))),
	colMeans(do.call(rbind,lapply(s1_hs,function(x) x[3,])))
	)[,c(2,3,5)]
rbind(
	colMeans(do.call(rbind,lapply(s1_hs2,function(x) x[1,]))),
	colMeans(do.call(rbind,lapply(s1_hs2,function(x) x[2,]))),
	colMeans(do.call(rbind,lapply(s1_hs2,function(x) x[3,])))
	)[,c(2,3,5)]

rbind(
	colMeans(do.call(rbind,lapply(s1_fs,function(x) x[1,]))),
	colMeans(do.call(rbind,lapply(s1_fs,function(x) x[2,]))),
	colMeans(do.call(rbind,lapply(s1_fs,function(x) x[3,])))
	)[,c(2,3,5)]
rbind(
	colMeans(do.call(rbind,lapply(s1_fs2,function(x) x[1,]))),
	colMeans(do.call(rbind,lapply(s1_fs2,function(x) x[2,]))),
	colMeans(do.call(rbind,lapply(s1_fs2,function(x) x[3,])))
	)[,c(2,3,5)]

rbind(
	colMeans(do.call(rbind,lapply(s1_fs,function(x) x[1,]))),
	colMeans(do.call(rbind,lapply(s1_fs,function(x) x[2,]))),
	colMeans(do.call(rbind,lapply(s1_fs,function(x) x[3,])))
	)[,c(2,3,5)]

# Scenario 2 - Maternal genetic effects only
s2_hs <- mclapply(1:100,function(i){
	out<-mge_sim(Va=0, Vmg=0.5, r_amg=0, Ve=0.5, p_sire=0, models=c(1,2,3,5))
	cat(i," ")
	out
}, mc.cores=8)

s2_fs <- mclapply(1:100,function(i){
	out<-mge_sim(Va=0, Vmg=0.5, r_amg=0, Ve=0.5, p_sire=1, models=c(1,2,3,5))
	cat(i," ")
	out
}, mc.cores=8)


# Scenario 3 - Both Uncorrelated

s3_hs <- mclapply(1:100,function(i){
	out<-mge_sim(Va=0.25, Vmg=0.25, r_amg=0, Ve=0.5, p_sire=0, models=c(1,2,3,5))
	cat(i," ")
	out
}, mc.cores=8)

s3_fs <- mclapply(1:100,function(i){
	out<-mge_sim(Va=0.25, Vmg=0.25, r_amg=0, Ve=0.5, p_sire=1, models=c(1,2,3,5))
	cat(i," ")
	out
}, mc.cores=8)



# Scenario 4 - Both Correlated


s4_hs <- mclapply(1:100,function(i){
	out<-mge_sim(Va=0.2, Vmg=0.2, r_amg=0.25, Ve=0.5, p_sire=0, models=c(1,2,3,5,6))
	cat(i," ")
	out
}, mc.cores=8)

s4_fs <- mclapply(1:100,function(i){
	out<-mge_sim(Va=0.2, Vmg=0.2, r_amg=0.25, Ve=0.5, p_sire=1, models=c(1,2,3,5,6))
	cat(i," ")
	out
}, mc.cores=8)


save(s1_hs,s1_fs,s2_hs,s2_fs,s3_hs,s3_fs,s4_hs,s4_fs, file="/Users/joelpick/Dropbox/0_fitness/fitness_timing/Data/Intermediate/gaussian_sims.Rdata")


# maternal_relatedness_ped <- function(ped){ 
# 	ped<-na.omit(ped) ## might depend on whether unknown males paternity likely means unknown male - in which case excluding will increase full-sib
# 	c(sapply(split(ped,ped$dam), function(x){
# 		if(nrow(x)==1){
# 		 NULL
# 		}else{
# 			fs <- NULL
# 			for(i in 1:(nrow(x)-1)){
# 				fs <- c(fs, x$sire[i] == x$sire[(i+1):nrow(x)])
# 			}
# 			mean(0.25 + fs*0.25)
# 		}
# 	}), recursive=TRUE)
# }

# maternal_relatedness_A <- function(A,dat){ 
# 	c(parallel::mclapply(split(dat,dat$dam), function(x){
# 		if(nrow(x)==1){
# 		 NULL
# 		}else{
# 			sub_A <-A[as.character(x$animal),as.character(x$animal)]
# 			mean(sub_A[which(lower.tri(sub_A))])
# 		}
# 	},mc.cores=8), recursive=TRUE)
# }

# 	ped <- simulate_pedigree(
# 		years = generations,
# 		n_females = n_females,
# 		fecundity = fecundity,
# 		p_sire = 0.5, 				# mating system (0-1, 1= one male per female, 0=complete random mating)
# 		juv_surv = 2/fecundity, # insures no population growth
# 		adult_surv = 0,					# discrete generations
# 		immigration = 0, 				# closed population
# 		constant_pop = TRUE     # constant population size
# 		)$pedigree

# hist(maternal_relatedness_ped(ped))

# A<-nadiv::makeA(ped[,1:3])
# hist(maternal_relatedness_A(A,ped))
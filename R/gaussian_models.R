rm(list=ls())

library(asreml)
load( file="/Users/joelpick/Dropbox/0_fitness/fitness_timing/Data/Intermediate/gaussian_data.Rdata")
## loads fs_peds,hs_peds,fs_data,hs_data


  ## Model 1 - additive genetic effects, assuming phenotype is trait of offspring
m1_func <- function(data){
	mod_1 <- asreml(
		fixed= p~1
		, random= ~vm(animal,ped.ainv)
	  , data= data, trace=FALSE)
	m1 <- summary(mod_1)$varcomp
	c(A= m1["vm(animal, ped.ainv)",1], 
		Me=NA, 
		Mg=NA,
		cov_AMg =NA,
		E = m1["units!R",1])
}

  ## Model 2 - additive genetic effects with maternal environment effects, assuming phenotype is trait of offspring
m2_func <- function(data){
	mod_2 <- asreml(
		fixed= p~1
    , random= ~vm(animal,ped.ainv) + mother_PE
    , data= data, trace=FALSE)
	m2<-summary(mod_2)$varcomp
	m2_sum <- c(
		A= m2["vm(animal, ped.ainv)",1], 
		Me=m2["mother_PE",1], 
		Mg=NA,
		cov_AMg =NA,
		E = m2["units!R",1])
}

m3_func <- function(data){
	mod_3 <- asreml(
		fixed= p~1
		, random= ~vm(animal,ped.ainv) + vm(mother,ped.ainv) 
		, data= data, trace=FALSE)
	m3<-summary(mod_3)$varcomp
	m3_sum <- c(
		A= m3["vm(animal, ped.ainv)",1], 
		Me=NA, 
		Mg=m3["vm(mother, ped.ainv)",1],
		cov_AMg =NA,
		E = m3["units!R",1])
}

m4_func <- function(data){
  mod_4 <- asreml(
		fixed= p~1
	  , random= ~vm(animal,ped.ainv) + vm(mother,ped.ainv) + mother_PE
	  , data= data, trace=FALSE)
	m4<-summary(mod_4)$varcomp
	m4_sum <- c(
		A= m4["vm(animal, ped.ainv)",1], 
		Me=m4["mother_PE",1], 
		Mg=m4["vm(mother, ped.ainv)",1],
		cov_AMg =NA,
		E = m4["units!R",1])
}

m5_func <- function(data){
	mod_5 <- asreml(
		fixed= p~1
    , random= ~str(~vm(animal,ped.ainv) +vm(mother,ped.ainv) ,~us(2):vm(animal,ped.ainv)) + mother_PE
    , data= data, trace=FALSE,maxit=50)
	m5<-summary(mod_5)$varcomp
	m5_sum <- c(
		A= m5[1,1], 
		Me=m5["mother_PE",1], 
		Mg=m5[3,1],
		cov_AMg =m5[2,1],
		E = m5["units!R",1])
}	

		assign("ped.ainv", asreml::ainverse(hs_peds[[2]]), envir = .GlobalEnv) 
		out <- do.call(rbind,parallel::mclapply(hs_data[[2]], m5_func, mc.cores=6))
		rm("ped.ainv", envir = .GlobalEnv)
		out
 
	mod_5 <- asreml(
		fixed= p~1
    , random= ~str(~vm(animal,ped.ainv) +vm(mother,ped.ainv) ,~us(2):vm(animal,ped.ainv)) + mother_PE
    , data= hs_data[[2]][[3]], trace=FALSE,maxit=50)


m6_func <- function(data){
	mod_6 <- asreml(
		fixed= p~1
    , random= ~vm(mother,ped.ainv)
    , data= data, trace=FALSE)
	m6<-summary(mod_6)$varcomp
	m6_sum <- c(
		A= NA, 
		Me=NA, 
		Mg=m6["vm(mother, ped.ainv)",1],
		cov_AMg =NA,
		E = m6["units!R",1])
}	


m7_func <- function(data){
	mod_7 <- asreml(
		fixed= p~1
	  , random= ~vm(mother,ped.ainv) + mother_PE
	  , data= data, trace=FALSE)
	m7<-summary(mod_7)$varcomp
	m7_sum <- c(
		A= NA, 
		Me=m7["mother_PE",1], 
		Mg=m7["vm(mother, ped.ainv)",1],
		cov_AMg =NA,
		E = m7["units!R",1])
}	


m8_func <- function(data){
	data_means <- aggregate(p~mother,data,mean)
	mod_8 <- asreml(
		fixed= p~1
    , random= ~vm(mother,ped.ainv)
    , data= data_means, trace=FALSE)
	m8<-summary(mod_8)$varcomp
	m8_sum <- c(
		A= NA, 
		Me=NA, 
		Mg=m8["vm(mother, ped.ainv)",1],
		cov_AMg =NA,
		E = m8["units!R",1])
}	


model_func <- function(FUN,peds,data){
	lapply(1:n_sims, function(i){
		ped <- peds[[i]]
		assign("ped.ainv", asreml::ainverse(ped), envir = .GlobalEnv) 
		out <- do.call(rbind,parallel::mclapply(data[[i]], FUN, mc.cores=6))
		rm("ped.ainv", envir = .GlobalEnv)
		cat(i," ")
		out
	})
}

model1_fs <- model_func(m1_func,fs_peds,fs_data)
model1_hs <- model_func(m1_func,hs_peds,hs_data)
model2_fs <- model_func(m2_func,fs_peds,fs_data)
model2_hs <- model_func(m2_func,hs_peds,hs_data)
model3_fs <- model_func(m3_func,fs_peds,fs_data)
model3_hs <- model_func(m3_func,hs_peds,hs_data)
model4_fs <- model_func(m4_func,fs_peds,fs_data)
model4_hs <- model_func(m4_func,hs_peds,hs_data)
model5_fs <- model_func(m5_func,fs_peds,fs_data) # many had this error: Likelihood evaluation failed with fault 1009 ; trying with reduced updates - This means it hasn't worked
model5_hs <- model_func(m5_func,hs_peds,hs_data)# Likelihood evaluation failed with fault 1009 ; trying with reduced updates
model6_fs <- model_func(m6_func,fs_peds,fs_data)
model6_hs <- model_func(m6_func,hs_peds,hs_data)
model7_fs <- model_func(m7_func,fs_peds,fs_data)
model7_hs <- model_func(m7_func,hs_peds,hs_data)
model8_fs <- model_func(m8_func,fs_peds,fs_data)
model8_hs <- model_func(m8_func,hs_peds,hs_data)



save(model1_hs,model1_fs,
	model2_hs,model2_fs,
	model3_hs,model3_fs,
	model4_hs,model4_fs, 
	model5_hs, model5_fs,
	model6_hs, model6_fs,
	model7_hs, model7_fs,
	model8_hs, model8_fs,  
	file="/Users/joelpick/Dropbox/0_fitness/fitness_timing/Data/Intermediate/gaussian_sims2.Rdata")

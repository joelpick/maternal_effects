
rm(list=ls())


library(asreml)
library(parallel)

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/extract_cousins.R"))
source(paste0(wd,"R/00_functions.R"))
source("/Users/joelpick/github/squidPed/R/simulate_pedigree.R")
# devtools::load_all("~/github/squidSim/R")

run=TRUE

n_sims <-100

##------------------------
## Pedigree parameters
##------------------------

generations=5

	load(paste0(data_wd,"mge_sims3.Rdata"))	

	scenarios[1,]

if(run){
	cores<-4
	set.seed(20230126)

	ped_str <- vector("list",length=nrow(peds_param))
	names(ped_str) <- ped_names
	
	# k="fs_lF_mI"
	
	for(k in ped_names){
		cat(k, "\n")
		cat("Simulating Pedigrees\n")
	## make pedigrees
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
		ped_str[[k]]<- do.call(rbind,mclapply(peds,ped_stat, mc.cores=cores))
		
		cat("Simulating Data\n")
		## simulate data
	
		dat<-mclapply(peds, function(i){	
			list(mge_sim(i[,1:3], param=scenarios[1,]))
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

		cat("\nModel 4: ")
		model4 <- model_func(m4_func,peds,dat,mc.cores=cores)
		assign(paste0("model4_",k),model4)

		cat("\n")
		rm(peds,dat)
	}

#paste0("model1_",ped_names),
	save(list=(c("ped_names","ped_str",paste0("model2_",ped_names),paste0("model4_",ped_names))),file=paste0(data_wd,"mge_sims_bv.Rdata"))
}


mat_ratio_all<-sapply(ped_str,function(x){
	rowSums(x[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(x[,-(1:2)]) 
})
mat_ratio <- colMeans(mat_ratio_all)


mod2<-do.call(rbind,lapply(ped_names,function(k) {
	mod2 <- do.call(rbind,lapply(get(paste0("model2_",k)), function(x) {
			data.frame(
				r=k,
				scenario=1,
				Va_est = x[1,"A"],
				Vm_est = x[1,"Me"],
				Va_sim=scenarios[1,"Va"],
				Vm_sim =sum(scenarios[1,c("Vmg","Vme")]),
				Vmg_sim=scenarios[1,"Vmg"])
	}))
	# assign(paste0("mod2_",k),mod2)
}))
mod2$Va_bias <- mod2$Va_est - mod2$Va_sim
mod2$Vm_bias <- mod2$Vm_est - mod2$Vm_sim


mod4<-do.call(rbind,lapply(ped_names,function(k) {
	mod4 <- do.call(rbind,lapply(get(paste0("model4_",k)), function(x) {
			data.frame(
				r=k,
				scenario=1,
				Va_est = x[1,"A"],
				Vme_est = x[1,"Me"],
				Vmg_est = x[1,"Mg"],
				Vm_est = x[1,"Me"] + x[1,"Mg"],
				Va_sim=scenarios[1,"Va"],
				Vm_sim =sum(scenarios[1,c("Vmg","Vme")]),
				Vmg_sim=scenarios[1,"Vmg"])
	}))
	# assign(paste0("mod2_",k),mod2)
}))
mod4$Va_bias <- mod4$Va_est - mod4$Va_sim
mod4$Vm_bias <- mod4$Vm_est - mod4$Vm_sim


# mod2$ln_Va_bias <- log(mod2$Va_bias)
va2<-aggregate(cbind(Va_bias,Vmg_sim,Vm_sim,Vm_bias)~ scenario+r, mod2,mean)
va2_prec<-aggregate(cbind(Va_est,Vm_est)~ scenario+r, mod2,function(x) sd(x))

va4<-aggregate(cbind(Va_bias,Vmg_sim,Vm_sim,Vm_bias)~ scenario+r, mod4,mean)
va4_prec<-aggregate(cbind(Va_est,Vm_est)~ scenario+r, mod4,function(x) sd(x))


va2_se<-aggregate(cbind(Va_bias,Vm_bias)~ scenario+r, mod2,se)

r_order<- sapply(va2$r, function(x) which(ped_names==x))
va4_prec$mat_ratio <- va2_prec$mat_ratio <- va4$mat_ratio <- va2$mat_ratio <- mat_ratio[r_order]
{
par(mfrow=c(2,2))
plot(Va_bias~mat_ratio,va2)
points(Va_bias~mat_ratio,va4, pch=19)

plot(Va_est~mat_ratio,va2_prec)
points(Va_est~mat_ratio,va4_prec, pch=19)

plot(Vm_bias~mat_ratio,va2)
points(Vm_bias~mat_ratio,va4, pch=19)

plot(Vm_est~mat_ratio,va2_prec)
points(Vm_est~mat_ratio,va4_prec, pch=19)
}


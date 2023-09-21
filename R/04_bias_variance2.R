
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
cores<-6

##------------------------
## Pedigree parameters
##------------------------

generations=3

fecundity=4

immigration <- c(f=0.5,m=0.5)

juv_surv <- c(
	f=2*(1 - immigration["f"])/fecundity,
	m=2*(1 - immigration["m"])/fecundity
)


scenario <-   c(Va = 0.1, Vmg= 0.2, Vme = 0.05, r_amg =0)
 

if(run){
	set.seed(20230920)
	
	cat("Simulating Pedigrees\n")
	## make pedigrees
		peds <- mclapply(1:n_sims,	function(i){
			simulate_pedigree(
				years = generations,
				n_females = 30,
				fecundity = fecundity,
				p_sire = 0.75,
				p_polyandry=1,
				juv_surv = juv_surv,
				adult_surv = 0,					# discrete generations
				immigration = immigration, 				# closed population
				constant_pop = TRUE     # constant population size
				)$pedigree
		}, mc.cores=cores)

		cat("Generating Pedigree Metrics\n")
		ped_str<- do.call(rbind,mclapply(peds,ped_stat, mc.cores=cores))
		
		cat("Simulating Data\n")
		## simulate data
	
		dat<-mclapply(peds, function(i){	
			list(mge_sim(i[,1:3], param=scenario))
		}, mc.cores=cores)
		# assign(paste0(k,"_data"),dat)

	## run models
		cat("Running models: \n")
	
		cat("Model 1: ")
		model1 <- model_func(m1_func,peds,dat,mc.cores=cores)

		cat("\nModel 2: ")
		model2 <- model_func(m1a_func,peds,dat,mc.cores=cores)
		
		cat("\nModel 3: ")
		model3 <- model_func(m2_func,peds,dat,mc.cores=cores)

		cat("\nModel 4: ")
		model4 <- model_func(m4_func,peds,dat,mc.cores=cores)
		
		# cat("\nModel 5: ")
		# model5 <- model_func(m5_func,peds,dat,mc.cores=cores)
		cat("\n")
		rm(peds,dat)
	}

#paste0("model1_",ped_names),
	save(list=(c("ped_names","ped_str",paste0("model2_",ped_names),paste0("model4_",ped_names))),file=paste0(data_wd,"mge_sims_bv.Rdata"))
}

	load(file=paste0(data_wd,"mge_sims_bv.Rdata"))

mat_ratio_all<-rowSums(ped_str[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(ped_str[,-(1:2)]) 

mat_ratio <- mean(mat_ratio_all)

mods <- do.call(rbind,lapply(1:n_sims, function(x) {
	data.frame(
		model=c(1:4),
		Va_est = c(model1[[x]][1,"A"],model2[[x]][1,"A"],model3[[x]][1,"A"],model4[[x]][1,"A"]),
		Vme_est = c(NA,model2[[x]][1,"Me"],model3[[x]][1,"Me"],model4[[x]][1,"Me"]),
		Vmg_est = c(NA,NA,NA,model4[[x]][1,"Mg"]),
		Vm_est = c(NA,model2[[x]][1,"Me"],model3[[x]][1,"Me"],sum(model4[[x]][1,c("Me","Mg")])),
		Va_sim=as.numeric(scenario["Va"]),
		Vm_sim =as.numeric(sum(scenario[c("Vmg","Vme")])),
		Vme_sim=as.numeric(scenario["Vme"]),
		Vmg_sim=as.numeric(scenario["Vmg"]))

}))

mods$Va_bias <- mods$Va_est - mods$Va_sim
mods$Vm_bias <- mods$Vm_est - mods$Vm_sim

# mod2$ln_Va_bias <- log(mod2$Va_bias)
va<-aggregate(cbind(Va_bias,Vmg_sim,Vm_sim,Vm_bias)~ model, mods,mean)
va_prec<-aggregate(cbind(Va_est,Vm_est)~ model, mods,function(x) sd(x))

va2_se<-aggregate(cbind(Va_bias,Vm_bias)~ model, mods,se)

par(mfrow=c(2,2))
boxplot(Va_est~model,mods); abline(h=scenario["Va"])
boxplot(Vm_est~model,mods); abline(h=scenario["Vme"] + scenario["Vmg"])

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


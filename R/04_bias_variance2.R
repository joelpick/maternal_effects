

rm(list=ls())


library(asreml)
library(parallel)
library(squidSim)
library(squidPed)

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/extract_cousins.R"))
source(paste0(wd,"R/00_functions.R"))
# source("/Users/joelpick/github/squidPed/R/simulate_pedigree.R")
# devtools::load_all("~/github/squidSim/R")





run=TRUE

n_sims <-100
cores<-6

##------------------------
## Pedigree parameters
##------------------------

generations=3

fecundity=3

immigration <- c(f=0.5,m=0.5)

juv_surv <- c(
	f=2*(1 - immigration["f"])/fecundity,
	m=2*(1 - immigration["m"])/fecundity
)


# Postma 2014, the median sample size of studies using animal models is 361
# in young and postma 2023 its 420
scenario <-   c(Va = 0.1, Vmg= 0.2, Vme = 0.05, r_amg =0)
 
ped<-simulate_pedigree(
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
nrow(ped)
sum(!is.na(ped$dam))

table()

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

		cat("\nModel 2: ")
		model2 <- model_func(m2_func,peds,dat,mc.cores=cores)

		cat("\nModel 4: ")
		model4 <- model_func(m4_func,peds,dat,mc.cores=cores)

		cat("\nModel 9: ")
		model9 <- model_func(m9_func,peds,dat,mc.cores=cores)

		cat("\nModel 5: ")
		model5 <- model_func(m5_func,peds,dat,mc.cores=cores)
		cat("\n")
		rm(peds,dat)
	}
			assign("ped.ainv", asreml::ainverse(peds[[2]]), envir = .GlobalEnv) 

	
m5_func(dat[[2]][[1]])


#paste0("model1_",ped_names),
	save(ped_names,ped_str,model2,
		model4,
		model9
		)),file=paste0(data_wd,"mge_sims_bv.Rdata"))
}

	load(file=paste0(data_wd,"mge_sims_bv.Rdata"))

mat_ratio_all<-rowSums(ped_str[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(ped_str[,-(1:2)]) 

mat_ratio <- mean(mat_ratio_all)

mods <- do.call(rbind,lapply(1:n_sims, function(x) {
	data.frame(
		model=c(2,4,9),
		Va_est = c(model2[[x]]$ml[1,"A"],model4[[x]]$ml[1,"A"],model9[[x]]$ml[1,"A"]),
		Vme_est = c(model2[[x]]$ml[1,"Me"],model4[[x]]$ml[1,"Me"],model9[[x]]$ml[1,"Me"]),
		Vmg_est = c(NA,model4[[x]]$ml[1,"Mg"],NA),
		Vml_est = c(NA,NA,model9[[x]]$ml[1,"Ml"]),
		Vm_est = c(model2[[x]]$ml[1,"Me"],sum(model4[[x]]$ml[1,c("Me","Mg")]),sum(model9[[x]]$ml[1,c("Me","Ml")])),
		Va_sim=as.numeric(scenario["Va"]),
		Vm_sim =as.numeric(sum(scenario[c("Vmg","Vme")])),
		Vme_sim=as.numeric(scenario["Vme"]),
		Vmg_sim=as.numeric(scenario["Vmg"]),

		Va_se=c(model2[[x]]$samp_cov["animal","animal",1],model4[[x]]$samp_cov["animal","animal",1],model9[[x]]$samp_cov["animal","animal",1]),
		Vmg_se=c(model2[[x]]$samp_cov["animal","animal",1],model4[[x]]$samp_cov["animal","animal",1],model9[[x]]$samp_cov["animal","animal",1])

		)

}))

mods$Va_bias <- mods$Va_est - mods$Va_sim
mods$Vm_bias <- mods$Vm_est - mods$Vm_sim

# mod2$ln_Va_bias <- log(mod2$Va_bias)
va<-aggregate(cbind(Va_bias,Vmg_sim,Vm_sim,Vm_bias)~ model, mods,mean)
aggregate(cbind(Va_est,Vm_est)~ model, mods,mean)

aggregate(cbind(Va_est,Vm_est, Va_se)~ model, mods,mean)

va_prec<-aggregate(cbind(Va_est,Vm_est)~ model, mods,function(x) sd(x))

va2_se<-aggregate(cbind(Va_bias,Vm_bias)~ model, mods,se)

par(mfrow=c(2,2))
boxplot(Va_est~model,mods); abline(h=scenario["Va"])
boxplot(Vm_est~model,mods); abline(h=scenario["Vme"] + scenario["Vmg"])
boxplot(Vme_est~model,mods); abline(h=scenario["Vme"])

boxplot(sqrt(Va_se)~model,mods); abline(h=scenario["Vme"] + scenario["Vmg"])

boxplot(Va_se~model,mods); abline(h=scenario["Vme"] + scenario["Vmg"])

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


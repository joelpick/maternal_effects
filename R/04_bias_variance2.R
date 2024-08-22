

rm(list=ls())

# devtools::install_github("squidgroup/squidSim")
# devtools::install_github("squidgroup/squidPed")

library(asreml)
library(parallel)
library(squidSim)
library(squidPed)

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))
load(paste0(data_wd,"parameters.Rdata"))

### WOULD A DAM SIRE MODEL WORK?

run=TRUE

cores<-6
n_sims <-100


# sapply( (peds_param_2[,"juv_surv_f"] * peds_param_2[,"fecundity"])/2 + peds_param_2[,"immigration_f"],all.equal,1)
# scenario <-   c(Va = 0.3, Vmg= 0.1, Vme = 0.1, r_amg =0)


###------------------------
### Pedigree and Phenotype Simulations, and Running Models
###------------------------

if(run){
	
	set.seed(20240822)

# k="fhs_lF_nI_small"
	## make pedigrees
	
	for(k in ped_names_reduced){
		cat(k, "\n")
		cat("Simulating Pedigrees\n")
		peds <- lapply(1:n_sims,	function(i){
			# print(i)
			# set.seed(85)
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
	
		# cat("Model 1: ")
		# model1 <- model_func(m1a_func,peds,dat,mc.cores=cores)
		# assign(paste0("model1_",k),model1)
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
}else{
	load(paste0(data_wd,"mge_sims_small_ped.Rdata"))
}







ped<-simulate_pedigree(
				years = peds_param_2[k,"generations"],
				n_females = peds_param_2[k,"n_females"],
				fecundity = peds_param_2[k,"fecundity"],
				fixed_fecundity = TRUE,
				p_sire = peds_param_2[k,"p_sire"],
				p_polyandry=1,
				p_breed=1,
				juv_surv = c(peds_param_2[k,"juv_surv_f"],peds_param_2[k,"juv_surv_m"]),
				adult_surv = 0,					# discrete generations
				immigration = c(peds_param_2[k,"immigration_f"],peds_param_2[k,"immigration_m"]),
				constant_pop = TRUE     # constant population size
				)$pedigree
nrow(ped)
sum(!is.na(ped$dam))

table(ped$cohort)

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

		# cat("\nModel 5: ")
		# model5 <- model_func(m5_func,peds,dat,mc.cores=cores)
		cat("\n")
		rm(peds,dat)
	}


	



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

		Va_se=sqrt(c(model2[[x]]$samp_cov["animal","animal",1],model4[[x]]$samp_cov["animal","animal",1],model9[[x]]$samp_cov["animal","animal",1])),
		Vme_se=sqrt(c(model2[[x]]$samp_cov["mother_PE","mother_PE",1],model4[[x]]$samp_cov["mother_PE","mother_PE",1],model9[[x]]$samp_cov["mother_PE","mother_PE",1])),
		Vmg_se=sqrt(c(NA,model4[[x]]$samp_cov["mother","mother",1],model9[[x]]$samp_cov["matriline","matriline",1]))
		)
}))

mods$Va_bias <- mods$Va_est - mods$Va_sim
mods$Vm_bias <- mods$Vm_est - mods$Vm_sim
mods$Va_z <- mods$Va_se/mods$Va_est
mods$Vmg_z <- mods$Vmg_se/mods$Vmg_est
mods$Vme_z <- mods$Vme_se/mods$Vme_est


# mod2$ln_Va_bias <- log(mod2$Va_bias)

va<-aggregate(cbind(Va_bias,Vm_bias,Va_est,Vm_est, Va_se)~ model, mods,mean,)


aggregate(cbind(Va_est,Vm_est, Va_se)~ model, mods,mean)

va_prec<-aggregate(cbind(Va_est,Vm_est)~ model, mods,function(x) sd(x))

aggregate(cbind(Va_bias)~ model, mods, function(x) sqrt(mean(x^2)))

par(mfrow=c(2,2))
boxplot(Va_est~model,mods); abline(h=scenario["Va"])
points(Va_est~c(1,2,3),va, pch=19, col="red")

boxplot(Vm_est~model,mods); abline(h=scenario["Vme"] + scenario["Vmg"])
points(Vm_est~c(1,2,3),va, pch=19, col="red")

# boxplot(Vme_est~model,mods); abline(h=scenario["Vme"])

boxplot(Va_se~model,mods)
points(Va_se~c(1,2,3),va, pch=19, col="red")

boxplot(Vme_se~model,mods)
boxplot(Vmg_se~model,mods)



library(beeswarm)
par(mfrow=c(3,2))
	beeswarm(Va_est~model,mods, pch=19, cex=0.4, col=scales::alpha(1,0.3),method = "compactswarm",corral="wrap")
points(Va_est~c(1,2,3),va, pch=19, col="red")

	beeswarm(Va_se~model,mods, pch=19, cex=0.4, col=scales::alpha(1,0.3),method = "compactswarm",corral="wrap")
	points(Va_se~c(1,2,3),va, pch=19, col="red")

# 	beeswarm(Vm_est~model,mods, pch=19, cex=0.4, col=scales::alpha(1,0.3),method = "compactswarm",corral="wrap")
# points(Vm_est~c(1,2,3),va, pch=19, col="red")

	beeswarm(Vmg_est~model,mods, pch=19, cex=0.4, col=scales::alpha(1,0.3),method = "compactswarm",corral="wrap")
points(Vmg_est~c(1,2,3),va, pch=19, col="red")


	beeswarm(Vmg_se~model,mods, pch=19, cex=0.4, col=scales::alpha(1,0.3),method = "compactswarm",corral="wrap")

		beeswarm(Vme_est~model,mods, pch=19, cex=0.4, col=scales::alpha(1,0.3),method = "compactswarm",corral="wrap")

	beeswarm(Vme_se~model,mods, pch=19, cex=0.4, col=scales::alpha(1,0.3),method = "compactswarm",corral="wrap")


plot(Vmg_est~Vme_est,mods)

par(mfrow=c(1,3))

	beeswarm(Va_z~model,mods, pch=19, cex=0.4, col=scales::alpha(1,0.3),method = "compactswarm",corral="wrap", ylim=c(0,10))
	beeswarm(Vmg_z~model,mods, pch=19, cex=0.4, col=scales::alpha(1,0.3),method = "compactswarm",corral="wrap", ylim=c(0,10))
	beeswarm(Vme_z~model,mods, pch=19, cex=0.4, col=scales::alpha(1,0.3),method = "compactswarm",corral="wrap", ylim=c(0,10))



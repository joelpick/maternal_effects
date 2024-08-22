
rm(list=ls())

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))
load(paste0(data_wd,"parameters.Rdata"))
	load(paste0(data_wd,"mge_sims_small_ped.Rdata"))

library(beeswarm)
library(scales)

mod1<-do.call(rbind,lapply(ped_names_reduced,function(k) {
	mod1 <- do.call(rbind,lapply(get(paste0("model1_",k)), function(x) {
			data.frame(
				r=k,
				scenario=1:nrow(scenarios),
				Va_est = x[["ml"]][,"A"],
				Va_sim=scenarios[,"Va"],
				Vm_sim = rowSums(scenarios[,c("Vmg","Vme")]),
				Vmg_sim=scenarios[,"Vmg"])
	}))
}))

mod2<-do.call(rbind,lapply(ped_names_reduced,function(k) {
	mod2 <- do.call(rbind,lapply(get(paste0("model2_",k)), function(x) {
			data.frame(
				r=k,
				scenario=1:nrow(scenarios),
				Va_est = x[["ml"]][,"A"],
				Vm_est = x[["ml"]][,"Me"],
				Va_sim=scenarios[,"Va"],
				Vm_sim = rowSums(scenarios[,c("Vmg","Vme")]),
				Vmg_sim=scenarios[,"Vmg"])
	}))
}))

mod4<-do.call(rbind,lapply(ped_names_reduced,function(k) {
	mod4 <- do.call(rbind,lapply(get(paste0("model4_",k)), function(x) {
			data.frame(
				r=k,
				scenario=1:nrow(scenarios),
				Va_est = x[["ml"]][,"A"],
				Vmg_est = x[["ml"]][,"Mg"],
				Vme_est = x[["ml"]][,"Me"],
				Va_sim=scenarios[,"Va"],
				Vm_sim = rowSums(scenarios[,c("Vmg","Vme")]),
				Vmg_sim=scenarios[,"Vmg"])
	}))
}))

order_exp <- expand.grid(imm=c("mI","nI","bI","fI"),fec=c("lF","mF","hF"),ms=c("fs","fhs","hs"))
order<-apply(order_exp,1,function(x) paste(x[3:1],collapse="_"))

mod2$order <- match(mod2$r,order)


par(mfrow=c(3,1))
beeswarm(Va_est~ r, subset(mod1,scenario==11),pch=19, cex=0.2, col=alpha(palette.colors()[1:4],0.5),method = "compactswarm",corral="wrap", ylab="Estimated Va", xaxt="n")
beeswarm(Va_est~ r, subset(mod2,scenario==11),pch=19, cex=0.2, col=alpha(palette.colors()[1:4],0.5),method = "compactswarm",corral="wrap", ylab="Estimated Va", xaxt="n")
beeswarm(Va_est~ r, subset(mod4,scenario==11),pch=19, cex=0.2, col=alpha(palette.colors()[1:4],0.5),method = "compactswarm",corral="wrap", ylab="Estimated Va", xaxt="n")



va<- rbind(
cbind(model=1,aggregate(Va_est ~ scenario+r, mod1,sd)),
cbind(model=2,aggregate(Va_est ~ scenario+r, mod2,sd)),
cbind(model=4,aggregate(Va_est ~ scenario+r, mod4,sd)))


beeswarm(Va_est~ model+r, va,pch=19, cex=0.4, col=alpha(palette.colors()[1:4],0.5),method = "compactswarm",corral="wrap", ylab="Estimated Va", xaxt="n")




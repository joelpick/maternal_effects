
rm(list=ls())

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))
load(paste0(data_wd,"parameters.Rdata"))
	load(paste0(data_wd,"mge_sims_small_ped.Rdata"))

library(beeswarm)
library(scales)


order_exp <- expand.grid(fec=c("lF"),ms=c("fhs"), model=c(1,2,4),imm=c("mI","nI","bI","fI"),size=c("small","medium"))
order<-apply(order_exp,1,function(x) paste(x[c(2,1,4,5,3)],collapse="_"))

mod1<-do.call(rbind,lapply(ped_names_reduced,function(k) {
	mod1 <- do.call(rbind,lapply(get(paste0("model1_",k)), function(x) {
			data.frame(
				r=k,
				scenario=1:nrow(scenarios),
				Va_est = x[["ml"]][,"A"],
				Va_sim=scenarios[,"Va"])
	}))
}))

mod2<-do.call(rbind,lapply(ped_names_reduced,function(k) {
	mod2 <- do.call(rbind,lapply(get(paste0("model2_",k)), function(x) {
			data.frame(
				r=k,
				scenario=1:nrow(scenarios),
				Va_est = x[["ml"]][,"A"],
				Va_sim=scenarios[,"Va"])
	}))
}))

mod4<-do.call(rbind,lapply(ped_names_reduced,function(k) {
	mod4 <- do.call(rbind,lapply(get(paste0("model4_",k)), function(x) {
			data.frame(
				r=k,
				scenario=1:nrow(scenarios),
				Va_est = x[["ml"]][,"A"],
				Va_sim=scenarios[,"Va"])
	}))
}))

all_mod<-rbind(
	cbind(model=1,mod1),
	cbind(model=2,mod2),
	cbind(model=4,mod4)
)
all_mod$r_model <- paste0(all_mod$r,"_",all_mod$model)
all_mod$order <- match(all_mod$r_model,order)





par(mfrow=c(1,1))
beeswarm(Va_est~ order, subset(all_mod,scenario==11),pch=19, cex=0.2, col=rep(alpha(palette.colors()[1:4],0.5), each=3),method = "compactswarm",corral="wrap", ylab="Estimated Va",las=2, label=order)



va<- rbind(
	cbind(model=1,aggregate(Va_est ~ scenario+r, mod1,mean)),
	cbind(model=2,aggregate(Va_est ~ scenario+r, mod2,mean)),
	cbind(model=4,aggregate(Va_est ~ scenario+r, mod4,mean)))

va$Va_sd<- rbind(
	cbind(model=1,aggregate(Va_est ~ scenario+r, mod1,sd)),
	cbind(model=2,aggregate(Va_est ~ scenario+r, mod2,sd)),
	cbind(model=4,aggregate(Va_est ~ scenario+r, mod4,sd)))$Va_est

va$r_model <- paste0(va$r,"_",va$model)
va$order <- match(va$r_model,order)


va$Va_z<-va$Va_sd/va$Va_est
# va$Va_z2<-va$Va_sd/va$Va_sim

par(mar=c(10,5,1,1))
beeswarm(Va_est~ order, subset(va,scenario==11),pch=19, cex=0.4, col=rep(alpha(palette.colors()[1:4],0.5), each=3),method = "compactswarm",corral="wrap", ylab="Estimated Va", las=2, label=order)
abline(h=0.25)

beeswarm(Va_sd~ order, subset(va,scenario==11),pch=19, cex=0.4, col=alpha(palette.colors()[1:4],0.5),method = "compactswarm",corral="wrap", ylab="Estimated Va", las=2)

beeswarm(Va_z~ order, subset(va,scenario==11),pch=19, cex=0.4, col=alpha(palette.colors()[1:4],0.5),method = "compactswarm",corral="wrap", ylab="Estimated Va", las=2)






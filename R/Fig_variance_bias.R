
rm(list=ls())

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))
load(paste0(data_wd,"parameters.Rdata"))
load(paste0(data_wd,"mge_sims_small_ped.Rdata"))

library(beeswarm)
library(scales)


order_exp <- expand.grid(fec=c("lF"),ms=c("fhs"), model=c(1,2,4,5),imm=c("mI","nI","bI","fI"),size=c("small","medium"))
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

mod5<-do.call(rbind,lapply(ped_names_reduced,function(k) {
	mod5 <- do.call(rbind,lapply(get(paste0("model5_",k)), function(x) {
			data.frame(
				r=k,
				scenario=1:nrow(scenarios),
				Va_est = x[["ml"]][,"A"],
				Va_sim=scenarios[,"Va"])
	}))
}))
## are these 0 when one of the variances is 0?

head(mod5)
mean(is.na(mod5$Va_est))
table(is.na(mod5$Va_est),mod5$scenario)
table(is.na(mod5$Va_est),mod5$r)

all_mod<-rbind(
	cbind(model=1,mod1),
	cbind(model=2,mod2),
	cbind(model=4,mod4),	
	cbind(model=5,mod5)
)
all_mod$r_model <- paste0(all_mod$r,"_",all_mod$model)
all_mod$order <- match(all_mod$r_model,order)

all_mod$Va_bias <- all_mod$Va_est - all_mod$Va_sim




par(mfrow=c(1,1),mar=c(10,5,1,1))
beeswarm(Va_bias~ order, subset(all_mod,scenario==11),pch=19, cex=0.2, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Bias in Va",las=2, label=order)


par(mfrow=c(1,1),mar=c(10,5,1,1))
beeswarm(Va_bias~ order, subset(all_mod,scenario==1),pch=19, cex=0.2, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Bias in Va",las=2, label=order)


par(mfrow=c(1,1),mar=c(10,5,1,1))
beeswarm(Va_bias~ order, subset(all_mod,scenario==6),pch=19, cex=0.2, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Bias in Va",las=2, label=order)




va<- aggregate(Va_est ~ scenario+r_model+r+model, all_mod,mean)

va$Va_sd<- aggregate(Va_est ~ scenario+r_model+r+model, all_mod,sd)$Va_est
va$Va_precision<- aggregate(Va_est ~ scenario+r_model+r+model, all_mod,function(x)1/sd(x))$Va_est
va$Va_rel_prec<-va$Va_est/va$Va_sd

va$Va_bias <- aggregate(Va_bias ~ scenario+r_model+r+model, all_mod, function(x) mean(x))$Va_bias

va$Va_abs_bias <- aggregate(Va_bias ~ scenario+r_model+r+model, all_mod, function(x) mean(abs(x)))$Va_bias

va$order <- match(va$r_model,order)


va$Va_z<-va$Va_sd/va$Va_est
# va$Va_z2<-va$Va_sd/va$Va_sim

ps<-11
pchs <- c(15,16,17)
{
	par(mfrow=c(5,1),mar=c(1,5,1,1))
beeswarm(Va_est~ order, subset(va,scenario==ps),pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Estimated Va", las=2)#, label=order)
abline(h=0.25)

# par(mar=c(10,5,1,1))
beeswarm(Va_bias~ order, subset(va,scenario==ps),pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Bias in Va", las=2)#, label=order)
abline(h=0)

# par(mar=c(10,5,1,1))
beeswarm(Va_sd~ order, subset(va,scenario==ps),pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="SE Va", las=2)#, label=order)

# par(mar=c(10,5,1,1))
beeswarm(Va_z~ order, subset(va,scenario==ps),pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Z Va", las=2)#, label=order)
abline(h=0)

# par(mar=c(10,5,1,1))
beeswarm(Va_abs_bias~ order, subset(va,scenario==ps),pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Absolute Error Va", las=2)#, label=order)
abline(h=0)
}


### 
ps<-11
pchs <- c(15,16,17)
{
	par(mfrow=c(4,1),mar=c(1,5,1,1))
beeswarm(Va_bias~ order, va,pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Bias in Va", las=2)#, 
abline(h=0)
beeswarm(Va_sd~ order, va,pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Precision Va", las=2)#, 
beeswarm(Va_z~ order, va,pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Realtive precision Va", las=2)#, 
beeswarm(Va_abs_bias~ order, va,pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Absolute Error Va", las=2)#, label=order)
}


pchs <- c(15,16,17)
cols <- alpha(palette.colors()[1:4],0.5)

plot_va <- subset(va, r %in% ped_names_reduced[grep("bI",ped_names_reduced)])
{
	par(mfrow=c(4,1),mar=c(1,5,1,1))
beeswarm(Va_bias~ order, plot_va,pch=pchs, cex=1, col=cols,method = "compactswarm",corral="wrap", ylab="Bias in Va", las=2)#, 
abline(h=0)
beeswarm(Va_precision~ order, plot_va,pch=pchs, cex=1, col=cols,method = "compactswarm",corral="wrap", ylab="Precision Va", las=2)#, 
beeswarm(Va_rel_prec~ order, plot_va,pch=pchs, cex=1, col=cols,method = "compactswarm",corral="wrap", ylab="Realtive precision Va", las=2)#, 
beeswarm(Va_abs_bias~ order, plot_va,pch=pchs, cex=1, col=cols,method = "compactswarm",corral="wrap", ylab="Absolute Error Va", las=2)#, label=order)
}


plot(Va_abs_bias~ order, va,pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=12), ylab="Mean Absolute Error in Va", las=2)#, 
for(i in 1:nrow(scenarios)){
	for(j in 1:length(ped_names_reduced)){
	lines(Va_abs_bias~order,subset(va,r==ped_names_reduced[j] & scenario==i), col=alpha(1,0.2))}}

{
	breaks=seq(-0.2,0.2,0.02)
par(mfrow=c(3,2))
par(mar=c(4,4,1,1))
## mod 2 better than mod 1

# hist(subset(va, model==1)$Va_abs_bias-subset(va, model==2)$Va_abs_bias, breaks=seq(-0.5,0.5,0.02))
# hist(subset(va, model==1 & scenario%in%c(7:10))$Va_abs_bias-subset(va, model==2& scenario%in%c(7:10))$Va_abs_bias, breaks=seq(-0.5,0.5,0.02), add=TRUE, col="blue")

## mod 4 better than mod 3
hist(subset(va, model==2)$Va_abs_bias-subset(va, model==4)$Va_abs_bias, breaks=breaks)
hist(subset(va, model==2 & scenario%in%c(3,7:12))$Va_abs_bias-subset(va, model==4& scenario%in%c(3,7:12))$Va_abs_bias, breaks=breaks, add=TRUE, col="red")
hist(subset(va, model==2 & scenario%in%c(3,11,12))$Va_abs_bias-subset(va, model==4& scenario%in%c(3,11,12))$Va_abs_bias, breaks=breaks, add=TRUE, col="blue")


## mod 5 better than mod 4
hist(subset(va, model==4)$Va_abs_bias-subset(va, model==5)$Va_abs_bias, breaks=breaks)

hist(subset(va, model==4 & scenario%in%c(7:10))$Va_abs_bias-subset(va, model==5& scenario%in%c(7:10))$Va_abs_bias, breaks=breaks, add=TRUE, col="blue")

hist(subset(va, model==4 & scenario%in%c(1:4))$Va_abs_bias-subset(va, model==5& scenario%in%c(1:4))$Va_abs_bias, breaks=breaks, add=TRUE, col="red")}
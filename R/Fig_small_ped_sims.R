
rm(list=ls())

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))
load(paste0(data_wd,"parameters.Rdata"))
load(paste0(data_wd,"mge_sims_small_ped.Rdata"))

library(beeswarm)
library(scales)


# order_exp <- expand.grid(fec=c("lF"),ms=c("fhs"), model=c(1,2,4,5),imm=c("mI","nI","bI","fI"),size=c("small","medium"))
# order<-apply(order_exp,1,function(x) paste(x[c(2,1,4,5,3)],collapse="_"))
order_exp <- expand.grid(ms=c("fhs"),fec=c("lF"),imm=c("mI","nI","bI","fI"), model=c(1,2,4,5),size=c("small","medium"))
order<-apply(order_exp,1,function(x) paste(x[c("ms","fec","imm","size","model")],collapse="_"))


scenarios <- cbind(scenarios, cov_amg= scenarios[,"r_amg"]*sqrt(scenarios[,"Va"] * scenarios[,"Vmg"])) 

extract_func <- function(x,k) {
	data.frame(
		r=k,
		scenario=1:nrow(scenarios),
		Va_est = x[["ml"]][,"A"],
		Vmg_est = x[["ml"]][,"Mg"],
		Vm_est = rowSums(x[["ml"]][,c("Mg","Me")], na.rm=TRUE),
		cov_amg_est = x[["ml"]][,"cov_AMg"],
		Va_sim=scenarios[,"Va"],
		Vmg_sim=scenarios[,"Vmg"],
		Vm_sim = scenarios[,"Vmg"]+scenarios[,"Vme"],
		cov_amg_sim = scenarios[,"cov_amg"]
	)
}


all_mod<-do.call(rbind,lapply(ped_names_reduced,function(k) {
	
	mod1 <- do.call(rbind,lapply(get(paste0("model1_",k)), extract_func,k=k))

	mod2 <- do.call(rbind,lapply(get(paste0("model2_",k)), extract_func,k=k))

	mod4 <- do.call(rbind,lapply(get(paste0("model4_",k)),
			extract_func,k=k))

	mod5 <- do.call(rbind,lapply(get(paste0("model5_",k)), extract_func,k=k))

	rbind(
		cbind(model=1,mod1),
		cbind(model=2,mod2),
		cbind(model=4,mod4),	
		cbind(model=5,mod5)
	)
}))

all_mod$r_model <- paste0(all_mod$r,"_",all_mod$model)
all_mod$order <- match(all_mod$r_model,order)

all_mod$Va_bias <- all_mod$Va_est - all_mod$Va_sim

all_mod$Vm_bias <- all_mod$Vm_est - all_mod$Vm_sim



## calculation of total Va

all_mod$Vmg_est <- ifelse(is.na(all_mod$Vmg_est),0,all_mod$Vmg_est)
all_mod$cov_amg_est <- ifelse(is.na(all_mod$cov_amg_est),0,all_mod$cov_amg_est)

all_mod$tVa_sim <- all_mod$Va_sim + 0.5*all_mod$Vmg_sim+ 1.5*all_mod$cov_amg_sim

all_mod$tVa_est <- all_mod$Va_est + 0.5*all_mod$Vmg_est+ 1.5*all_mod$cov_amg_est

all_mod$tVa_bias <- all_mod$tVa_est - all_mod$tVa_sim






########
## MAKE SUMMARIES
########

va<- aggregate(cbind(Va_bias,Vm_bias,tVa_bias) ~ scenario+r_model+r+model, all_mod,mean)


va$Va_precision<- aggregate(Va_est ~ scenario+r_model+r+model, all_mod,function(x)1/sd(x))$Va_est

va$Va_rel_prec<- aggregate(Va_est ~ scenario+r_model+r+model, all_mod,function(x)mean(x)/sd(x))$Va_est

va$Va_abs_bias <- aggregate(Va_bias ~ scenario+r_model+r+model, all_mod, function(x) mean(abs(x)))$Va_bias



va$tVa_bias<- aggregate(tVa_bias ~ scenario+r_model+r+model, all_mod,mean)$tVa_bias
va$tVa_abs_bias<- aggregate(tVa_bias ~ scenario+r_model+r+model, all_mod, function(x) mean(abs(x)))$tVa_bias

va$order <- match(va$r_model,order)






### 
pchs <- c(15,16,17)
{
	par(mfrow=c(4,1),mar=c(1,5,1,1))
beeswarm(Va_bias~ order, va,pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Bias in Va", las=2)#, 
abline(h=0)
beeswarm(Va_precision~ order, va,pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Precision Va", las=2)#, 
beeswarm(Va_rel_prec~ order, va,pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Relative precision Va", las=2)#, 
beeswarm(Va_abs_bias~ order, va,pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Absolute Error Va", las=2)#, label=order)
}


beeswarm(tVa_bias~ order, va,pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="bias in tVa", las=2)#, label=order)
abline(h=0)
beeswarm(tVa_abs_bias~ order, va,pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="bias in tVa", las=2)#, label=order)
abline(h=0)


pchs <- c(21:24)
cols <- alpha(palette.colors()[1:4],0.5)

plot_va <- subset(va, r %in% ped_names_reduced[grep("bI",ped_names_reduced)])
{
	par(mfrow=c(4,1),mar=c(4,5,1,1))
beeswarm(Va_bias~ order, plot_va,pch=pchs, cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", xlab='Model',ylab="Bias in Va", label=c(1:4,1:4))#, 
abline(h=0)
beeswarm(Va_precision~ order, plot_va,pch=pchs, cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", xlab='Model',ylab="Precision Va", label=c(1:4,1:4))#, 
beeswarm(Va_rel_prec~ order, plot_va,pch=pchs, cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", xlab='Model',ylab="Realtive precision Va", label=c(1:4,1:4))#, 
beeswarm(Va_abs_bias~ order, plot_va,pch=pchs, cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", xlab='Model',ylab="Absolute Error Va", label=c(1:4,1:4))#, label=order)
}

	par(mfrow=c(1,1),mar=c(4,5,1,1))

s_col <- rep(alpha(palette.colors()[1:4],0.5), 8)[va$order]
plot(Va_abs_bias~ order, va,pch=pchs, cex=1, col=s_col, bg=s_col ylab="Mean Absolute Error in Va", las=2)#, 
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



par(mfrow=c(2,2))
hist(subset(va, model==2)$Va_abs_bias-subset(va, model==4)$Va_abs_bias, breaks=breaks)


hist(subset(va, model==2)$tVa_abs_bias-subset(va, model==4)$tVa_abs_bias, breaks=breaks)

hist(subset(va, model==4)$Va_abs_bias-subset(va, model==5)$Va_abs_bias, breaks=breaks)


hist(subset(va, model==4)$tVa_abs_bias-subset(va, model==5)$tVa_abs_bias, breaks=breaks)


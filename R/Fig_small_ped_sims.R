
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

# code above calcaultes a vm also for modelt aht didnt converge - undo this
all_mod$Vm_est <- ifelse(is.na(all_mod$Va_est),NA,all_mod$Vm_est)

all_mod$Vm_bias <- all_mod$Vm_est - all_mod$Vm_sim



## calculation of total Va

all_mod$Vmg_est <- ifelse(all_mod$model %in% c(1,2) & is.na(all_mod$Vmg_est),0,all_mod$Vmg_est)
all_mod$cov_amg_est <- ifelse(all_mod$model %in% c(1,2,4) & is.na(all_mod$cov_amg_est),0,all_mod$cov_amg_est)

all_mod$tVa_sim <- all_mod$Va_sim + 0.5*all_mod$Vmg_sim+ 1.5*all_mod$cov_amg_sim

all_mod$tVa_est <- all_mod$Va_est + 0.5*all_mod$Vmg_est+ 1.5*all_mod$cov_amg_est

all_mod$tVa_bias <- all_mod$tVa_est - all_mod$tVa_sim



########
## MAKE SUMMARIES
########

va<- aggregate(cbind(Va_bias,Vm_bias,tVa_bias) ~ scenario+r_model+r+model, all_mod,mean)

va[,c("Va_precision","Vm_precision","tVa_precision")] <- aggregate(cbind(Va_est,Vm_est,tVa_est) ~ scenario+r_model+r+model, all_mod, function(x)1/sd(x))[,c("Va_est","Vm_est","tVa_est")]

va[,c("Va_rel_prec","Vm_rel_prec","tVa_rel_prec")] <- aggregate(cbind(Va_est,Vm_est,tVa_est) ~ scenario+r_model+r+model, all_mod, function(x)mean(x)/sd(x))[,c("Va_est","Vm_est","tVa_est")]

va[,c("Va_abs_bias","Vm_abs_bias","tVa_abs_bias")] <- aggregate(cbind(Va_bias,Vm_bias,tVa_bias) ~ scenario+r_model+r+model, all_mod, function(x) mean(abs(x)))[,c("Va_bias","Vm_bias","tVa_bias")]

va$order <- match(va$r_model,order)



####
#  mean difference in absolute error across all simulations  
####
mod2 <- subset(all_mod,model==2)
mod4 <- subset(all_mod,model==4)
mod5 <- subset(all_mod,model==5)

## difference 
va_4_5 <- tapply(abs(mod4$Va_bias) - abs(mod5$Va_bias), list(mod2$r_model,mod2$scenario),mean, na.rm=TRUE)

vm_4_5 <- tapply(abs(mod4$Vm_bias) - abs(mod5$Vm_bias), list(mod2$r_model,mod2$scenario),mean, na.rm=TRUE)

tva_4_5 <- tapply(abs(mod4$tVa_bias) - abs(mod5$tVa_bias), list(mod2$r_model,mod2$scenario),mean, na.rm=TRUE)


va_2_4 <- tapply(abs(mod2$Va_bias) - abs(mod4$Va_bias), list(mod2$r_model,mod2$scenario),mean, na.rm=TRUE)

vm_2_4 <- tapply(abs(mod2$Vm_bias) - abs(mod4$Vm_bias), list(mod2$r_model,mod2$scenario),mean, na.rm=TRUE)

tva_2_4 <- tapply(abs(mod2$tVa_bias) - abs(mod4$tVa_bias), list(mod2$r_model,mod2$scenario),mean, na.rm=TRUE)


### 

hl_hist <- function(x, breaks, main="", ylim=NULL, col1, col2,...){
	if(is.null(ylim)) ylim<-range(hist(x,plot=FALSE,breaks=breaks)$counts)
	hist(x[x>=0], breaks=breaks, main=main, col=scales::alpha(col1,0.8),ylim=ylim,...)
	hist(x[x<0], breaks=breaks, main=main, col=scales::alpha(col2,0.8), add=TRUE)
}


pchs <- c(21:24)
cols <- alpha(palette.colors()[1:4],0.5)

plot_va <- subset(va, r %in% ped_names_reduced[grep("bI_small",ped_names_reduced)])

plot_va$Vm_abs_bias <- ifelse(plot_va$model==1,NA,plot_va$Vm_abs_bias)
plot_va$model<-ifelse(plot_va$model==4, 3,
	ifelse(plot_va$model==5, 4,
	plot_va$model))

breaks=seq(-0.2,0.2,0.02)

setEPS()
pdf(paste0(wd,"Figures/Fig7_small_ped.pdf"), height=6, width=10)
{
layout(matrix(c(1:15),ncol=3,byrow=TRUE), width=c(1.5,1,1), height=c(1,4,4,4,1))

	par(mar=c(0,5,0,0))#mfrow=c(3,1),
	plot(NA,xlim=c(0.5,4.5),ylim=c(0,1),bty="n",xaxt="n",yaxt="n", ylab="")
	arrows(1,0.2,4,0.2,code=2, xpd=TRUE, length=0.15, lwd=3)
	text(2.5,0.6,"Increasing model complexity", cex=1.5)

	plot(NA,xlim=c(0.5,4.5),ylim=c(0,1),bty="n",xaxt="n",yaxt="n", ylab="")
	text(2.5,0.6,"Model 2 - Model 3", cex=1.5)

	plot(NA,xlim=c(0.5,4.5),ylim=c(0,1),bty="n",xaxt="n",yaxt="n", ylab="")
	text(2.5,0.6,"Model 3 - Model 4", cex=1.5)

par(mar=c(1,5,0,1), las=1, cex.lab=1.5)
beeswarm(Va_abs_bias~ order, plot_va,pch=pchs, cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap",ylab=expression(V[A]~Absolute~Error), xaxt="n", ylim=c(0,0.6))
# for(i in 1:nrow(scenarios)){
# 	lines(Va_abs_bias~model,subset(plot_va, scenario==i), col=alpha(1,0.1))
# }


par(mar=c(2,5,0.5,0.5))
hl_hist(as.vector(va_2_4), breaks=breaks, col1=palette.colors()[3],col2=palette.colors()[2], ylim=c(0,60), ylab="")
hl_hist(va_4_5, breaks=breaks, ylim=c(0,60), col1=palette.colors()[4],col2=palette.colors()[3], ylab="")

par(mar=c(1,5,0,1))
beeswarm(Vm_abs_bias~ order, plot_va,pch=pchs, cex=1, col=cols[2:4],bg=cols[2:4],method = "compactswarm",corral="wrap", xaxt="n",ylab=
expression(V[M]~Absolute~Error), xlim=c(-0.5,3.5), ylim=c(0,0.6))
legend (-0.5,0.5,c("1","2","3","4"), pch=pchs, col=cols, pt.bg=cols, title="Model",bty="n", cex=1.25)
# for(i in 1:nrow(scenarios)){
# 	lines(Vm_abs_bias~I(model-1),subset(plot_va, scenario==i), col=alpha(1,0.1))
# }

par(mar=c(2,5,0.5,0.5))
hl_hist(vm_2_4, breaks=breaks, col1=palette.colors()[3],col2=palette.colors()[2], ylim=c(0,60))
hl_hist(vm_4_5, breaks=breaks, ylim=c(0,60), col1=palette.colors()[4],col2=palette.colors()[3], ylab="")

par(mar=c(1,5,0,1))
beeswarm(tVa_abs_bias~ order, plot_va,pch=pchs, cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap",ylab=expression(V[At]~Absolute~Error), ylim=c(0,0.6), xaxt="n")
axis(1,1:4,1:4, cex.axis=1.25)
mtext("Model",1,outer=TRUE,adj=0.23, line=-1.5, cex=1.5)

# for(i in 1:nrow(scenarios)){
# 	lines(tVa_abs_bias~model,subset(plot_va, scenario==i), col=alpha(1,0.1))
# }

par(mar=c(2,5,0.5,0.5))
hl_hist(tva_2_4, breaks=breaks, col1=palette.colors()[3],col2=palette.colors()[2], ylim=c(0,60), ylab="")
hl_hist(tva_4_5, breaks=breaks, ylim=c(0,60), col1=palette.colors()[4],col2=palette.colors()[3], ylab="")

mtext("Mean Difference in Absolute Error",1,outer=TRUE,adj=0.8, line=-1.5, cex=1.5)

}
dev.off()


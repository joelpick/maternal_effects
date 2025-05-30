
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
head(model4_fhs_lF_nI_small)

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

# code above calculates a vm also for model that didnt converge - undo this
all_mod$Vm_est <- ifelse(is.na(all_mod$Va_est),NA,all_mod$Vm_est)

all_mod$Vm_bias <- all_mod$Vm_est - all_mod$Vm_sim




## calculation of total Va

all_mod$Vmg_est <- ifelse(all_mod$model %in% c(1,2) & is.na(all_mod$Vmg_est),0,all_mod$Vmg_est)
all_mod$cov_amg_est <- ifelse(all_mod$model %in% c(1,2,4) & is.na(all_mod$cov_amg_est),0,all_mod$cov_amg_est)

all_mod$tVa_sim <- all_mod$Va_sim + 0.5*all_mod$Vmg_sim+ 1.5*all_mod$cov_amg_sim

all_mod$tVa_est <- all_mod$Va_est + 0.5*all_mod$Vmg_est+ 1.5*all_mod$cov_amg_est

all_mod$tVa_bias <- all_mod$tVa_est - all_mod$tVa_sim



all_mod$Vmg_bias <- all_mod$Vmg_est - all_mod$Vmg_sim
all_mod$cov_amg_bias <- all_mod$cov_amg_est - all_mod$cov_amg_sim

all_mod$cor_amg_est <- all_mod$cov_amg_est/(sqrt(all_mod$Vmg_est)*sqrt(all_mod$Va_est))
all_mod$cor_amg_est[all_mod$cor_amg_est>1|all_mod$cor_amg_est< -1]<-NA


hist(all_mod$cor_amg_est[all_mod$model==5], breaks=100)
all_mod$cor_amg_sim <- all_mod$cov_amg_sim/(sqrt(all_mod$Vmg_sim)*sqrt(all_mod$Va_sim))
all_mod$cor_amg_bias <- all_mod$cov_amg_est - all_mod$cov_amg_sim

subset(all_mod,cor_amg_est>1|cor_amg_est< -1)
$cor_amg_est

########
## MAKE SUMMARIES
########

va<- aggregate(cbind(Va_bias,Vm_bias,tVa_bias,Vmg_bias,cov_amg_bias,cor_amg_bias) ~ scenario+r_model+r+model, all_mod,mean)

va[,c("Va_precision","Vm_precision","tVa_precision")] <- aggregate(cbind(Va_est,Vm_est,tVa_est) ~ scenario+r_model+r+model, all_mod, function(x)1/sd(x))[,c("Va_est","Vm_est","tVa_est")]

va[,c("Va_rel_prec","Vm_rel_prec","tVa_rel_prec")] <- aggregate(cbind(Va_est,Vm_est,tVa_est) ~ scenario+r_model+r+model, all_mod, function(x)mean(x)/sd(x))[,c("Va_est","Vm_est","tVa_est")]

va[,c("Va_abs_bias","Vm_abs_bias","tVa_abs_bias")] <- aggregate(cbind(Va_bias,Vm_bias,tVa_bias) ~ scenario+r_model+r+model, all_mod, function(x) mean(abs(x)))[,c("Va_bias","Vm_bias","tVa_bias")]

va[,c("Va_rmse","Vm_rmse","tVa_rmse")] <- aggregate(cbind(Va_bias,Vm_bias,tVa_bias) ~ scenario+r_model+r+model, all_mod, function(x) sqrt(mean(x^2)))[,c("Va_bias","Vm_bias","tVa_bias")]

va$order <- match(va$r_model,order)

## remove models where Vm wasnt estimated
va$Vm_rmse <- ifelse(va$model==1,NA,va$Vm_rmse)
va$Vm_abs_bias <- ifelse(va$model==1,NA,va$Vm_abs_bias)


plot_va <- subset(va, r %in% ped_names_reduced[grep("bI_small",ped_names_reduced)])

# plot_va$Vm_abs_bias <- ifelse(plot_va$model==1,NA,plot_va$Vm_abs_bias)


plot_va$model<-ifelse(plot_va$model==4, 3,
	ifelse(plot_va$model==5, 4,
	plot_va$model))


va$ped_size<- substring(va$r,11)
va$scenario2<- ifelse(va$scenario%in%c(5:10),"both",
	ifelse(va$scenario%in%c(1:4), "no_va"
		,"no_vmg"
		))
va$scenario3<- ifelse(va$scenario%in%c(7:8),"cov_positive",
	ifelse(va$scenario%in%c(9:10), "cov_negative"
		,va$scenario2
		))

vmg <- subset(va, model %in% c(4,5))





plot_bias <- function(x, scenarios, cols,ylab=x,xlab="",pchs=19,bgs=cols,ylim=c(0,0.6),...){
	form <- as.formula(paste(x,"~ model")) 
	plot(form,subset(plot_va, scenario==scenarios[1]), type="b",ylim=ylim, col=cols[1], xaxt="n", pch=pchs[1],ylab=ylab,xlab=xlab,bg=bgs[1],...)#, ylim=range(plot_va[,x],na.rm=TRUE)
	axis(1,1:4,1:4)
	for(i in 2:length(scenarios)){
	lines(form,subset(plot_va, scenario==scenarios[i]), type="b", col=cols[i], pch=pchs[i],bg=bgs[i])	
	}
}



cols <- c(palette.colors(),palette.colors()[1:3])
bgs <- c(palette.colors(),rep("white",3))
pchs <- rep(21:25,3)
{
par(mfrow=c(2,3))
plot_bias("Vmg_bias",5:10, pch=pchs[5:10], cols=cols[5:10], bgs=bgs[5:10], ylab=expression(italic(hat(V)[M])~RMSE),ylim=c(-0.6,0.6))

plot_bias("Vmg_bias",1:4, pch=pchs[1:4], cols=cols[1:4], bgs=bgs[1:4], ylab=expression(italic(hat(V)[M])~RMSE),ylim=c(-0.6,0.6))

plot_bias("Vmg_bias",11:12, pch=pchs[11:12], cols=cols[11:12], bgs=bgs[11:12], ylab="",ylim=c(-0.6,0.6))

plot_bias("cov_amg_bias",5:10, pch=pchs[5:10], cols=cols[5:10], bgs=bgs[5:10], ylab=expression(italic(hat(V)[M])~RMSE),ylim=c(-0.6,0.6))

plot_bias("cov_amg_bias",1:4, pch=pchs[1:4], cols=cols[1:4], bgs=bgs[1:4], ylab=expression(italic(hat(V)[M])~RMSE),ylim=c(-0.6,0.6))

plot_bias("cov_amg_bias",11:12, pch=pchs[11:12], cols=cols[11:12], bgs=bgs[11:12], ylab="",ylim=c(-0.6,0.6))
}

cols <- alpha(palette.colors()[1:5],0.5)

vmg$cov_amg_bias <- ifelse(vmg$model == 4,0,vmg$cov_amg_bias)

{par(mfrow=c(3,1), mar=c(0,5,4,1))
	beeswarm(Va_bias~ scenario2 + ped_size + model, vmg,pch=c(15:17), cex=1, pwcol=cols[as.factor(vmg$scenario3)],bg=cols,method = "compactswarm",corral="wrap", ylab=expression(Bias~"in"~italic(hat(V)[A])), labels=c("both","no_va","no_vmg"))
		abline(h=0)
	abline(v=(1:7)*3+0.5, col=alpha(c("grey","black"),0.5))
		axis(3,c(2,5,8,11),c("Medium","Small","Medium","Small"), lwd.ticks=0, padj=1, cex.axis=1.25)

		axis(3,c(2,3.5,5),c(" ","model 3 (no cov)"," "), lwd.ticks=0, padj=1, cex.axis=1.25, line=2)
		axis(3,c(8,9.5,11),c(" ","model 4 (cov)"," "), lwd.ticks=0, padj=1, cex.axis=1.25, line=2)

par(mar=c(2,5,2,1))

	beeswarm(Vmg_bias~ scenario2 + ped_size + model, vmg,pch=c(15:17), cex=1, pwcol=cols[as.factor(vmg$scenario3)],bg=cols,method = "compactswarm",corral="wrap", ylab=expression(Bias~"in"~italic(hat(V)[Mg])))
		abline(h=0)
	abline(v=(1:7)*3+0.5, col=alpha(c("grey","black"),0.5))

par(mar=c(3,5,1,1))

	beeswarm(cor_amg_bias~ scenario2 + ped_size + model, vmg,pch=c(15:17), cex=1, pwcol=cols[as.factor(vmg$scenario3)],bg=cols,method = "compactswarm",corral="wrap", ylab=expression(Bias~"in"~italic(hat(COV)[A,Mg])))
		abline(h=0)
	abline(v=(1:7)*3+0.5, col=alpha(c("grey","black"),0.5))
}

## make into correlation


{par(mfrow=c(3,1), mar=c(2,5,1,1))
pchs <- c(21:24)
cols <- alpha(palette.colors()[1:4],0.5)


	beeswarm(Va_bias~ order, vmg,pwpch=c(15:17)[as.factor(vmg$scenario2)], cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", ylab=expression(Bias~"in"~italic(hat(V)[A])), las=2,xaxt="n", xlim=c(0,16))#, 
	abline(h=0)
	abline(v=(1:7)*4+0.5, col=alpha(c("grey","black"),0.5))
legend("topleft", c("both","no_va","no_vmg"),pch=c(15:17), )

	# par(mar=c(0.5,7,0.5,0.5), cex.lab=1.5)#mfrow=c(3,1),
	beeswarm(Vmg_bias~ order, vmg,pwpch=c(15:17)[as.factor(vmg$scenario2)], cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", ylab=expression(Bias~"in"~italic(hat(V)[Mg])), las=2,xaxt="n", xlim=c(0,16))#, 
	abline(h=0)
	abline(v=(1:7)*4+0.5, col=alpha(c("grey","black"),0.5))
legend("topleft", c("both","no_va","no_vmg"),pch=c(15:17), )

pch=c(15:17)[as.factor(va$scenario2)]

	# par(mar=c(0.5,7,0.5,0.5), cex.lab=1.5)#mfrow=c(3,1),
	beeswarm(cov_amg_bias~ order, vmg,pwpch=c(15:17)[as.factor(vmg$scenario2)], cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", ylab=expression(Bias~"in"~italic(hat(COV)[A,Mg])), las=2,xaxt="n")#, 
	abline(h=0)
	abline(v=4.5, col=alpha("grey",0.5))
legend("bottomright", c("both","no_va","no_vmg"),pch=c(15:17), )

va$ped_size<- substring(va$r,11)
va$scenario2<- ifelse(va$scenario%in%c(5:10),"both",
	ifelse(va$scenario%in%c(1:4), "no_va","no_vmg"
		))
}

par(mfrow=c(1,3))
plot(Vmg_bias~ cov_amg_bias, subset(va,model==5), col=c(1:3)[as.factor(va$scenario2)], pch=c(15:16)[as.factor(va$ped_size)])
abline(h=0,v=0)

plot(Va_bias~ cov_amg_bias, subset(va,model==5), col=c(1:3)[as.factor(va$scenario2)], pch=c(15:16)[as.factor(va$ped_size)])
abline(h=0,v=0)

plot(Vmg_bias~ Va_bias, subset(va,model==5), col=c(1:3)[as.factor(va$scenario2)], pch=c(15:16)[as.factor(va$ped_size)])
abline(h=0,v=0)


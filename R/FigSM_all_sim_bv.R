
rm(list=ls())

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))

load(paste0(data_wd,"mge_sims3.Rdata"))
load(paste0(data_wd,"parameters.Rdata"))
load(paste0(data_wd,"mge_sims_small_ped.Rdata"))

library(beeswarm)
library(scales)

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
		cov_amg_sim = scenarios[,"cov_amg"])
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


order_exp <- expand.grid(ms=c("fhs"),fec=c("lF"),imm=c("mI","nI","bI","fI"),size=c("small","medium"), model=c(1,2,4,5))
order<-apply(order_exp,1,function(x) paste(x[c("ms","fec","imm","size","model")],collapse="_"))

all_mod$r_model <- paste0(all_mod$r,"_",all_mod$model)
all_mod$order <- match(all_mod$r_model,order)


for(j in 1:12){
	setEPS()
	pdf(paste0(wd,"Figures/FigSM_all_sim_bv",j,".pdf"), height=8, width=13)
{	
	
	layout(matrix(c(1,1,1,2,2,2,3,3,4),nrow=3, byrow=TRUE))

	par(mar=c(1,7,5,0.5), cex.lab=1.5)#mfrow=c(3,1),
	
	beeswarm(Va_est~ order, subset(all_mod,scenario==j),pch=19, cex=0.3, col=alpha(palette.colors()[1:4],0.5),method = "compactswarm",corral="wrap", ylab=expression(Estimated~V[A]),las=2, xaxt="n", xlim=c(1.5,31.5))
	abline(h=scenarios[j,"Va"], col=alpha(palette.colors()[8],0.5))
	abline(v=(1:7)*4+0.5, col=alpha(c("grey","black"),0.5))

	for(i in c(0,8,16,24)){
		axis(3,c(1,4.5,8)+i,c("",paste("Model",i/8+1),""), lwd.ticks=0, line=2.5, padj=1, cex.axis=1.25)
		axis(3,c(1,2.5,4)+i,c("","Small",""), lwd.ticks=0, line=0.5, padj=1, cex.axis=1.25)
		axis(3,c(1,2.5,4)+4+i,c("","Medium",""), lwd.ticks=0, line=0.5, padj=1, cex.axis=1.25)
	}

	mtext("Model", side=3, line=-2.5, outer=TRUE, adj=0)
	mtext("Pedigree Size", side=3, line=-4.5, outer=TRUE, adj=0)

####

	par(mar=c(6,7,0,0.5))

	beeswarm(Vm_est~ order, subset(all_mod,scenario==j),pch=19, cex=0.2, col=alpha(palette.colors()[1:4],0.5),method = "compactswarm",corral="wrap", ylab=expression(Estimated~V[M]), labels=c("M","N","U","F"), xlab="", xlim=c(1.5,31.5))
	abline(v=(1:7)*4+0.5, col=alpha(c("grey","black"),0.5))
	abline(h=scenarios[j,"Vmg"]+scenarios[j,"Vme"], col=alpha(palette.colors()[8],0.5))
	
	mtext("Immigration", side=2, line=1, adj=0.9, las=1, padj=10)


####
		par(mar=c(3,7,3,0.5))

	beeswarm(Vmg_est~ order, subset(all_mod,scenario==j & model %in% c(4,5)),pch=19, cex=0.2, col=alpha(palette.colors()[1:4],0.5),method = "compactswarm",corral="wrap", ylab=expression(Estimated~V[Mg]), labels=c("M","N","U","F"), xlab="", xlim=c(1,16))
	abline(h=scenarios[j,"Vmg"], col=alpha(palette.colors()[8],0.5))
	abline(v=(1:3)*4+0.5, col=alpha(c("grey","black"),0.5))

	for(i in c(0,8)){
		axis(3,c(1,4.5,8)+i,c("",paste("Model",i/8+3),""), lwd.ticks=0, line=2.5, padj=1, cex.axis=1.25)
		axis(3,c(1,2.5,4)+i,c("","Small",""), lwd.ticks=0, line=0.5, padj=1, cex.axis=1.25)
		axis(3,c(1,2.5,4)+4+i,c("","Medium",""), lwd.ticks=0, line=0.5, padj=1, cex.axis=1.25)
	}

	mtext("Immigration", side=2, line=1, adj=0.9, las=1, padj=10)


	beeswarm(cov_amg_est~ order, subset(all_mod,scenario==j & model ==5),pch=19, cex=0.2, col=alpha(palette.colors()[1:4],0.5),method = "compactswarm",corral="wrap", ylab=expression(Estimated~COV["A,Mg"]), labels=c("M","N","U","F"), xlab="", xlim=c(0.5,8.5))
	abline(h=scenarios[j,"cov_amg"], col=alpha(palette.colors()[8],0.5))
	abline(v=4.5, col=alpha(c("grey"),0.5))
	axis(3,c(1,4.5,8),c("","Model 4",""), lwd.ticks=0, line=2.5, padj=1, cex.axis=1.25)
	axis(3,c(1,2.5,4),c("","Small",""), lwd.ticks=0, line=0.5, padj=1, cex.axis=1.25)
	axis(3,c(1,2.5,4)+4,c("","Medium",""), lwd.ticks=0, line=0.5, padj=1, cex.axis=1.25)

}	

dev.off()
}


#############
## CONVERGENCE
#############

head(mod5)
mean(is.na(mod5$Va_est))
par(mfrow=c(2,1))
barplot(table(is.na(mod5$Va_est),mod5$scenario))
barplot(table(is.na(mod5$Va_est),mod5$r))




########
## MAKE SUMMARIES
########

va<- aggregate(cbind(Va_bias,Vm_bias,tVa_bias) ~ scenario+r_model+r+model, all_mod,mean)

va[,c("Va_precision","Vm_precision","tVa_precision")] <- aggregate(cbind(Va_est,Vm_est,tVa_est) ~ scenario+r_model+r+model, all_mod, function(x)1/sd(x))[,c("Va_est","Vm_est","tVa_est")]

va[,c("Va_rel_prec","Vm_rel_prec","tVa_rel_prec")] <- aggregate(cbind(Va_est,Vm_est,tVa_est) ~ scenario+r_model+r+model, all_mod, function(x)mean(x)/sd(x))[,c("Va_est","Vm_est","tVa_est")]

va[,c("Va_abs_bias","Vm_abs_bias","tVa_abs_bias")] <- aggregate(cbind(Va_bias,Vm_bias,tVa_bias) ~ scenario+r_model+r+model, all_mod, function(x) mean(abs(x)))[,c("Va_bias","Vm_bias","tVa_bias")]

va$order <- match(va$r_model,order)


pchs <- c(15,16,17)

#direct comparison with meyer
{
va$Va_sd<- aggregate(Va_est ~ scenario+r_model+r+model, all_mod,sd)$Va_est

par(mar=c(10,5,1,1))
beeswarm(Va_sd~ order, subset(va,scenario==11),pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="SE Va", las=2, label=order)
}



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

{
		vm <- subset(va, model%in% c(2,4,5))

	par(mfrow=c(4,1),mar=c(1,5,1,1))
beeswarm(Vm_bias~ order, vm,pch=pchs, cex=1, col=rep(alpha(palette.colors()[2:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Bias in Vm", las=2)#, 
abline(h=0)
beeswarm(Vm_precision~ order, vm,pch=pchs, cex=1, col=rep(alpha(palette.colors()[2:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Precision Vm", las=2)#, 
beeswarm(Vm_rel_prec~ order, vm,pch=pchs, cex=1, col=rep(alpha(palette.colors()[2:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Relative precision Vm", las=2)#, 
beeswarm(Vm_abs_bias~ order, vm,pch=pchs, cex=1, col=rep(alpha(palette.colors()[2:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Absolute Error Vm", las=2)#, label=order)
}

{
	par(mfrow=c(4,1),mar=c(1,5,1,1))
beeswarm(tVa_bias~ order, va,pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Bias in total Va", las=2)#, 
abline(h=0)
beeswarm(tVa_precision~ order, va,pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Precision total Va", las=2)#, 
beeswarm(tVa_rel_prec~ order, va,pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Relative precision total Va", las=2)#, 
beeswarm(tVa_abs_bias~ order, va,pch=pchs, cex=1, col=rep(alpha(palette.colors()[1:4],0.5), each=4),method = "compactswarm",corral="wrap", ylab="Absolute Error total Va", las=2)#, label=order)
}


plot(subset(va,r=="fhs_lF_nI_small")$Va_bias,subset(va,r=="fhs_lF_mI_small")$Va_bias)

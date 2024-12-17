
rm(list=ls())

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))

load(paste0(data_wd,"mge_sims3.Rdata"))
load(paste0(data_wd,"parameters.Rdata"))
load(paste0(data_wd,"mge_sims_small_ped.Rdata"))

library(xtable)
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


order_exp <- expand.grid(ms=c("fhs"),fec=c("lF"),imm=c("mI","nI","bI","fI"),size=c("small","medium"), model=c(1,2,4,5))
order<-apply(order_exp,1,function(x) paste(x[c("ms","fec","imm","size","model")],collapse="_"))

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




pchs <- c(21:24)
cols <- alpha(palette.colors()[1:4],0.5)




for(j in 1:12){
	setEPS()
	pdf(paste0(wd,"Figures/FigSM_all_sim_bv",j,".pdf"), height=8, width=13)
{	
	
	layout(matrix(c(1,1,1,2,2,2,3,3,4),nrow=3, byrow=TRUE))

	par(mar=c(1,7,5,0.5), cex.lab=1.5)#mfrow=c(3,1),
	
	beeswarm(Va_est~ order, subset(all_mod,scenario==j),pch=19, cex=0.3, col=cols,method = "compactswarm",corral="wrap", ylab=expression(Estimated~V[A]),las=2, xaxt="n", xlim=c(1.5,31.5))
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

	beeswarm(Vm_est~ order, subset(all_mod,scenario==j),pch=19, cex=0.2, col=cols,method = "compactswarm",corral="wrap", ylab=expression(Estimated~V[M]), labels=c("M","N","U","F"), xlab="", xlim=c(1.5,31.5))
	abline(v=(1:7)*4+0.5, col=alpha(c("grey","black"),0.5))
	abline(h=scenarios[j,"Vmg"]+scenarios[j,"Vme"], col=alpha(palette.colors()[8],0.5))
	
	mtext("Immigration", side=2, line=1, adj=0.9, las=1, padj=10)


####
		par(mar=c(3,7,3,0.5))

	beeswarm(Vmg_est~ order, subset(all_mod,scenario==j & model %in% c(4,5)),pch=19, cex=0.2, col=cols,method = "compactswarm",corral="wrap", ylab=expression(Estimated~V[Mg]), labels=c("M","N","U","F"), xlab="", xlim=c(1,16))
	abline(h=scenarios[j,"Vmg"], col=alpha(palette.colors()[8],0.5))
	abline(v=(1:3)*4+0.5, col=alpha(c("grey","black"),0.5))

	for(i in c(0,8)){
		axis(3,c(1,4.5,8)+i,c("",paste("Model",i/8+3),""), lwd.ticks=0, line=2.5, padj=1, cex.axis=1.25)
		axis(3,c(1,2.5,4)+i,c("","Small",""), lwd.ticks=0, line=0.5, padj=1, cex.axis=1.25)
		axis(3,c(1,2.5,4)+4+i,c("","Medium",""), lwd.ticks=0, line=0.5, padj=1, cex.axis=1.25)
	}

	mtext("Immigration", side=2, line=1, adj=0.9, las=1, padj=10)


	beeswarm(cov_amg_est~ order, subset(all_mod,scenario==j & model ==5),pch=19, cex=0.2, col=cols,method = "compactswarm",corral="wrap", ylab=expression(Estimated~COV["A,Mg"]), labels=c("M","N","U","F"), xlab="", xlim=c(0.5,8.5))
	abline(h=scenarios[j,"cov_amg"], col=alpha(palette.colors()[8],0.5))
	abline(v=4.5, col=alpha(c("grey"),0.5))
	axis(3,c(1,4.5,8),c("","Model 4",""), lwd.ticks=0, line=2.5, padj=1, cex.axis=1.25)
	axis(3,c(1,2.5,4),c("","Small",""), lwd.ticks=0, line=0.5, padj=1, cex.axis=1.25)
	axis(3,c(1,2.5,4)+4,c("","Medium",""), lwd.ticks=0, line=0.5, padj=1, cex.axis=1.25)

}	

dev.off()
}





########
## MAKE SUMMARIES
########

va<- aggregate(cbind(Va_bias,Vm_bias,tVa_bias) ~ scenario+r_model+r+model, all_mod,mean)

va[,c("Va_precision","Vm_precision","tVa_precision")] <- aggregate(cbind(Va_est,Vm_est,tVa_est) ~ scenario+r_model+r+model, all_mod, function(x)1/sd(x))[,c("Va_est","Vm_est","tVa_est")]

va[,c("Va_rel_prec","Vm_rel_prec","tVa_rel_prec")] <- aggregate(cbind(Va_est,Vm_est,tVa_est) ~ scenario+r_model+r+model, all_mod, function(x)mean(x)/sd(x))[,c("Va_est","Vm_est","tVa_est")]

va[,c("Va_abs_bias","Vm_abs_bias","tVa_abs_bias")] <- aggregate(cbind(Va_bias,Vm_bias,tVa_bias) ~ scenario+r_model+r+model, all_mod, function(x) mean(abs(x)))[,c("Va_bias","Vm_bias","tVa_bias")]

va$order <- match(va$r_model,order)





#############
### SUMMARY PLOTS
#############

# similarity of results across pedigrees

abs_bias_cor<-round(cor(sapply(ped_names_reduced,function(x)subset(va,r==x)$Va_abs_bias)),3)
abc_names <- sub("fhs_lF_","",rownames(abs_bias_cor))
abc_names <- sub("_"," ",abc_names)
abc_names <- sub("bI","uI",abc_names)

colnames(abs_bias_cor)<-rownames(abs_bias_cor)<-abc_names
abs_bias_cor[upper.tri(abs_bias_cor)]<-""

xtable(abs_bias_cor)





#direct comparison with meyer
# order_exp2 <- expand.grid(ms=c("fhs"),fec=c("lF"), model=c(1,2,4,5),imm=c("mI","nI","bI","fI"),size=c("small","medium"))
# order2<-apply(order_exp2,1,function(x) paste(x[c("ms","fec","imm","size","model")],collapse="_"))
# va$order2 <- match(va$r_model,order2)


# setEPS()
# pdf(paste0(wd,"Figures/FigSM_small_ped_meyer.pdf"), height=8, width=13)

# {
# va$Va_sd<- aggregate(Va_est ~ scenario+r_model+r+model, all_mod,sd)$Va_est

# ## split into two plots by ped size
# par(mfrow=c(2,1),mar=c(10,5,1,1))
# beeswarm(Va_sd~ order2, subset(va,scenario==11),cex=1, pch=pchs,col=cols,bg=cols,method = "compactswarm",corral="wrap", ylab="SE Va", las=2, label=order)
# beeswarm(Va_sd~ order2, subset(va,scenario==11),cex=1, pch=pchs,col=cols,bg=cols,method = "compactswarm",corral="wrap", ylab="SE Va", las=2, label=order)
# }

# dev.off()





### 


pchs <- c(21:24)
cols <- alpha(palette.colors()[1:4],0.5)


setEPS()
pdf(paste0(wd,"Figures/FigSM_small_ped_Va.pdf"), height=8, width=13)
{
	layout(matrix(1:6,ncol=1), height=c(1,2,2,2,2,1))
	# par(mfrow=c(4,1))
	# par(mar=c(1,7,5,0.5), cex.lab=1.5)#mfrow=c(3,1),
	par(mar=c(0,0,0,0))#mfrow=c(3,1),
	plot(NA,xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")


	par(mar=c(0.5,7,0.5,0.5), cex.lab=1.5)#mfrow=c(3,1),
	beeswarm(Va_bias~ order, va,pch=pchs, cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", ylab=expression(Bias~"in"~V[A]), las=2,xaxt="n", xlim=c(1.5,31.5))#, 
	abline(h=0)
	abline(v=(1:7)*4+0.5, col=alpha(c("grey","black"),0.5))

	for(i in c(0,8,16,24)){
		axis(3,c(1,4.5,8)+i,c("",paste("Model",i/8+1),""), lwd.ticks=0, line=2.5, padj=1, cex.axis=1.25)
		axis(3,c(1,2.5,4)+i,c("","Small",""), lwd.ticks=0, line=0.5, padj=1, cex.axis=1.25)
		axis(3,c(1,2.5,4)+4+i,c("","Medium",""), lwd.ticks=0, line=0.5, padj=1, cex.axis=1.25)
	}
	
	beeswarm(Va_precision~ order, va,pch=pchs, cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", ylab="Precision", las=2,xaxt="n", xlim=c(1.5,31.5))#, 
		abline(v=(1:7)*4+0.5, col=alpha(c("grey","black"),0.5))

	beeswarm(Va_rel_prec~ order, va,pch=pchs, cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", ylab="Relative Precision", las=2,xaxt="n", xlim=c(1.5,31.5))#, 
		abline(v=(1:7)*4+0.5, col=alpha(c("grey","black"),0.5))

	beeswarm(Va_abs_bias~ order, va,pch=pchs, cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", ylab="Mean Absolute Error", las=2,xaxt="n", xlim=c(1.5,31.5))#, label=order)
		abline(v=(1:7)*4+0.5, col=alpha(c("grey","black"),0.5))
		axis(1,1:32,rep(c("M","N","U","F"),8))
}
dev.off()


setEPS()
pdf(paste0(wd,"Figures/FigSM_small_ped_tVa.pdf"), height=8, width=13)
{
	layout(matrix(1:6,ncol=1), height=c(1,2,2,2,2,1))
	# par(mfrow=c(4,1))
	# par(mar=c(1,7,5,0.5), cex.lab=1.5)#mfrow=c(3,1),
	par(mar=c(0,0,0,0))#mfrow=c(3,1),
	plot(NA,xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")


	par(mar=c(0.5,7,0.5,0.5), cex.lab=1.5)#mfrow=c(3,1),
	beeswarm(tVa_bias~ order, va,pch=pchs, cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", ylab=expression(Bias~"in"~V[At]), las=2,xaxt="n", xlim=c(1.5,31.5))#, 
	abline(h=0)
	abline(v=(1:7)*4+0.5, col=alpha(c("grey","black"),0.5))

	for(i in c(0,8,16,24)){
		axis(3,c(1,4.5,8)+i,c("",paste("Model",i/8+1),""), lwd.ticks=0, line=2.5, padj=1, cex.axis=1.25)
		axis(3,c(1,2.5,4)+i,c("","Small",""), lwd.ticks=0, line=0.5, padj=1, cex.axis=1.25)
		axis(3,c(1,2.5,4)+4+i,c("","Medium",""), lwd.ticks=0, line=0.5, padj=1, cex.axis=1.25)
	}
	
	beeswarm(tVa_precision~ order, va,pch=pchs, cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", ylab="Precision", las=2,xaxt="n", xlim=c(1.5,31.5))#, 
		abline(v=(1:7)*4+0.5, col=alpha(c("grey","black"),0.5))

	beeswarm(tVa_rel_prec~ order, va,pch=pchs, cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", ylab="Relative Precision", las=2,xaxt="n", xlim=c(1.5,31.5))#, 
		abline(v=(1:7)*4+0.5, col=alpha(c("grey","black"),0.5))

	beeswarm(tVa_abs_bias~ order, va,pch=pchs, cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", ylab="Mean Absolute Error", las=2,xaxt="n", xlim=c(1.5,31.5))#, label=order)
		abline(v=(1:7)*4+0.5, col=alpha(c("grey","black"),0.5))
		axis(1,1:32,rep(c("M","N","U","F"),8))
}
dev.off()


## only model 2-4 have maternal effects estimated
vm <- subset(va, model%in% c(2,4,5))

setEPS()
pdf(paste0(wd,"Figures/FigSM_small_ped_Vm.pdf"), height=8, width=13)
{

	layout(matrix(1:6,ncol=1), height=c(1,2,2,2,2,1))
	# par(mfrow=c(4,1))
	# par(mar=c(1,7,5,0.5), cex.lab=1.5)#mfrow=c(3,1),
	par(mar=c(0,0,0,0))#mfrow=c(3,1),
	plot(NA,xlim=c(0,1),ylim=c(0,1),bty="n",xaxt="n",yaxt="n")


	par(mar=c(0.5,7,0.5,0.5), cex.lab=1.5)#mfrow=c(3,1),
	beeswarm(Vm_bias~ order, vm,pch=pchs, cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", ylab=expression(Bias~"in"~V[M]), las=2,xaxt="n", xlim=c(1.5,23.5))#, 
	abline(h=0)
	abline(v=(1:5)*4+0.5, col=alpha(c("grey","black"),0.5))

	for(i in c(0,8,16)){
		axis(3,c(1,4.5,8)+i,c("",paste("Model",i/8+2),""), lwd.ticks=0, line=2.5, padj=1, cex.axis=1.25)
		axis(3,c(1,2.5,4)+i,c("","Small",""), lwd.ticks=0, line=0.5, padj=1, cex.axis=1.25)
		axis(3,c(1,2.5,4)+4+i,c("","Medium",""), lwd.ticks=0, line=0.5, padj=1, cex.axis=1.25)
	}
	
	beeswarm(Vm_precision~ order, vm,pch=pchs, cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", ylab="Precision", las=2,xaxt="n", xlim=c(1.5,23.5))#, 
		abline(v=(1:5)*4+0.5, col=alpha(c("grey","black"),0.5))

	beeswarm(Vm_rel_prec~ order, vm,pch=pchs, cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", ylab="Relative Precision", las=2,xaxt="n", xlim=c(1.5,23.5))#, 
		abline(v=(1:5)*4+0.5, col=alpha(c("grey","black"),0.5))

	beeswarm(Vm_abs_bias~ order, vm,pch=pchs, cex=1, col=cols,bg=cols,method = "compactswarm",corral="wrap", ylab="Mean Absolute Error", las=2,xaxt="n", xlim=c(1.5,23.5))#, label=order)
		abline(v=(1:5)*4+0.5, col=alpha(c("grey","black"),0.5))
		axis(1,1:24,rep(c("M","N","U","F"),6))
}
dev.off()



#############
## CONVERGENCE
#############
order_cov <- expand.grid(fec=c("lF"),ms=c("fhs"),imm=c("mI","nI","bI","fI"),size=c("small","medium"))
order_cov2<-apply(order_cov,1,function(x) paste(x[c(2,1,3,4)],collapse="_"))

mod5 <- subset(all_mod,model==5)

mean(is.na(mod5$Va_est))

setEPS()
pdf(paste0(wd,"Figures/FigSM_small_ped_convergence.pdf"), height=8, width=7)

{
par(mfrow=c(2,1), mar=c(6,5,2,1))
converge_scenario<-table(is.na(mod5$Va_est),mod5$scenario)

barplot(converge_scenario, xlab="Scenario", names=LETTERS[1:12], ylab="Number of models")
barplot(converge_scenario[,1:4], xlab="Scenario",col=scales::alpha("red",c(0.9,0.3)), add=TRUE,xaxt="n")
mtext("A)",side=3,line=-2,adj=0.025, outer=TRUE, cex=1.25)

d<-barplot(table(is.na(mod5$Va_est),mod5$r)[,order_cov2], xlab="Pedigree", names=rep(c("M","N","U","F"),2), col=grey.colors(2), ylab="Number of models")
axis(3,c(d[1],mean(d[1:4]),d[4]),c("","Small",""), lwd.ticks=0, line=0.5, padj=1, cex.axis=1.25)
axis(3,c(d[5],mean(d[5:8]),d[8]),c("","Medium",""), lwd.ticks=0, line=0.5, padj=1, cex.axis=1.25)
mtext("B)",side=3,line=-22,adj=0.025, outer=TRUE, cex=1.25)
}
dev.off()


### comparing results from model 3 in datasets where model did and did not converge

mod5_sub <- subset(all_mod,model==5 & scenario %in% 1:4)
mod4_sub <- subset(all_mod,model==4 & scenario %in% 1:4)
## va estimated in model 3 when model 4 doesnt converge
mod4_sub$na_mod5 <- !is.na(mod5_sub$Va_est)
mod5_sub$na_mod5 <- !is.na(mod5_sub$Va_est)
mod45_sub <- rbind(mod4_sub,mod5_sub)


conv_sum<- cbind(tapply(mod4_sub$Va_bias,list(mod4_sub$r,mod4_sub$na_mod5), function(x) mean(abs(x))),tapply(mod45_sub$Va_bias,list(mod45_sub$r,mod45_sub$model), function(x) mean(abs(x),na.rm=TRUE)))

rownames(conv_sum) <- sub("fhs_lF_","",rownames(conv_sum))

colnames(conv_sum) <- c("m4 not converged","m4 converged","all","Model 4")

xtable(conv_sum,digits=4)

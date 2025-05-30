
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

# code above calculates a vm also for model that didnt converge - undo this
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






plot_bias <- function(x, scenarios, cols,ylab=x,xlab="",pchs=19,bgs=cols,...){
	form <- as.formula(paste(x,"~ model")) 
	plot(form,subset(plot_va, scenario==scenarios[1]), type="b",ylim=c(0,0.6), col=cols[1], xaxt="n", pch=pchs[1],ylab=ylab,xlab=xlab,bg=bgs[1],...)#, ylim=range(plot_va[,x],na.rm=TRUE)
	axis(1,1:4,1:4)
	for(i in 2:length(scenarios)){
	lines(form,subset(plot_va, scenario==scenarios[i]), type="b", col=cols[i], pch=pchs[i],bg=bgs[i])	
	}
}



cols <- c(palette.colors(),palette.colors()[1:3])
bgs <- c(palette.colors(),rep("white",3))
pchs <- rep(21:25,3)

setEPS()
pdf(paste0(wd,"Figures/Fig7_small_ped.pdf"), height=5.5, width=8)
{
# par(mfrow=c(3,3))
	layout(matrix(c(1,2,3,11,4,5,6,10,7,8,9,12),nrow=3, byrow=TRUE), width=c(5,5,5,1))


par(mar=c(1,5,5,1))
plot_bias("Va_abs_bias",5:10, ylab=expression(italic(hat(V)[A])~Absolute~Error), main=expression(V[A]~and~V[Mg]), pch=pchs[5:10], cols=cols[5:10], bgs=bgs[5:10])
mtext("A)",side=3, at=0, cex=1, las=1, line=0.5)

plot_bias("Va_abs_bias",1:4, ylab="", main=expression(No~V[A]), pch=pchs[1:4], cols=cols[1:4], bgs=bgs[1:4])
mtext("B)",side=3, at=0, cex=1, las=1, line=0.5)

plot_bias("Va_abs_bias",11:12, ylab="", main=expression(No~V[Mg]), pch=pchs[11:12], cols=cols[11:12], bgs=bgs[11:12])
# plot_bias("Va_abs_bias",5:6, cols=cols[5:6] )
mtext("C)",side=3, at=0, cex=1, las=1, line=0.5)

par(mar=c(3,5,3,1))
plot_bias("Vm_abs_bias",5:10, pch=pchs[5:10], cols=cols[5:10], bgs=bgs[5:10], ylab=expression(italic(hat(V)[M])~Absolute~Error))
mtext("D)",side=3, at=0, cex=1, las=1, line=0.5)

plot_bias("Vm_abs_bias",1:4, pch=pchs[1:4], cols=cols[1:4], bgs=bgs[1:4], ylab="")
mtext("E)",side=3, at=0, cex=1, las=1, line=0.5)

plot_bias("Vm_abs_bias",11:12, pch=pchs[11:12], cols=cols[11:12], bgs=bgs[11:12], ylab="")
mtext("F)",side=3, at=0, cex=1, las=1, line=0.5)
# plot_bias("Vm_abs_bias",5:6, cols=cols[5:6] )

par(mar=c(5,5,1,1))
plot_bias("tVa_abs_bias",5:10, pch=pchs[5:10], cols=cols[5:10], bgs=bgs[5:10], ylab=expression(italic(hat(V)[At])~Absolute~Error))
mtext("G)",side=3, at=0, cex=1, las=1, line=0.5)

plot_bias("tVa_abs_bias",1:4, pch=pchs[1:4], cols=cols[1:4], bgs=bgs[1:4], ylab="", xlab="Model")
mtext("H)",side=3, at=0, cex=1, las=1, line=0.5)

plot_bias("tVa_abs_bias",11:12, pch=pchs[11:12], cols=cols[11:12], bgs=bgs[11:12], ylab="")
mtext("I)",side=3, at=0, cex=1, las=1, line=0.5)
# plot_bias("tVa_abs_bias",5:6, cols=cols[5:6] )


	par(mar=c(0,0,0,0))
	plot(NA, xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), xlab="",ylab="",bty="n")
	legend("center",letters[1:12], pch=pchs, pt.bg=bgs, col=cols, bty="n", cex=1,title="Scenario")

}
dev.off()



setEPS()
pdf(paste0(wd,"Figures/FigSM_small_ped_rmse.pdf"), height=5.5, width=8)
{
# par(mfrow=c(3,3))
	layout(matrix(c(1,2,3,11,4,5,6,10,7,8,9,12),nrow=3, byrow=TRUE), width=c(5,5,5,1))


par(mar=c(1,5,5,1))
plot_bias("Va_rmse",5:10, ylab=expression(italic(hat(V)[A])~RMSE), main=expression(V[A]~and~V[Mg]), pch=pchs[5:10], cols=cols[5:10], bgs=bgs[5:10])
mtext("A)",side=3, at=0, cex=1, las=1, line=0.5)

plot_bias("Va_rmse",1:4, ylab="", main=expression(No~V[A]), pch=pchs[1:4], cols=cols[1:4], bgs=bgs[1:4])
mtext("B)",side=3, at=0, cex=1, las=1, line=0.5)

plot_bias("Va_rmse",11:12, ylab="", main=expression(No~V[Mg]), pch=pchs[11:12], cols=cols[11:12], bgs=bgs[11:12])
# plot_bias("Va_rmse",5:6, cols=cols[5:6] )
mtext("C)",side=3, at=0, cex=1, las=1, line=0.5)

par(mar=c(3,5,3,1))
plot_bias("Vm_rmse",5:10, pch=pchs[5:10], cols=cols[5:10], bgs=bgs[5:10], ylab=expression(italic(hat(V)[M])~RMSE))
mtext("D)",side=3, at=0, cex=1, las=1, line=0.5)

plot_bias("Vm_rmse",1:4, pch=pchs[1:4], cols=cols[1:4], bgs=bgs[1:4], ylab="")
mtext("E)",side=3, at=0, cex=1, las=1, line=0.5)

plot_bias("Vm_rmse",11:12, pch=pchs[11:12], cols=cols[11:12], bgs=bgs[11:12], ylab="")
mtext("F)",side=3, at=0, cex=1, las=1, line=0.5)
# plot_bias("Vm_rmse",5:6, cols=cols[5:6] )

par(mar=c(5,5,1,1))
plot_bias("tVa_rmse",5:10, pch=pchs[5:10], cols=cols[5:10], bgs=bgs[5:10], ylab=expression(italic(hat(V)[At])~RMSE))
mtext("G)",side=3, at=0, cex=1, las=1, line=0.5)

plot_bias("tVa_rmse",1:4, pch=pchs[1:4], cols=cols[1:4], bgs=bgs[1:4], ylab="", xlab="Model")
mtext("H)",side=3, at=0, cex=1, las=1, line=0.5)

plot_bias("tVa_rmse",11:12, pch=pchs[11:12], cols=cols[11:12], bgs=bgs[11:12], ylab="")
mtext("I)",side=3, at=0, cex=1, las=1, line=0.5)
# plot_bias("tVa_abs_bias",5:6, cols=cols[5:6] )


	par(mar=c(0,0,0,0))
	plot(NA, xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), xlab="",ylab="",bty="n")
	legend("center",letters[1:12], pch=pchs, pt.bg=bgs, col=cols, bty="n", cex=1,title="Scenario")

}
dev.off()









setEPS()
pdf(paste0(wd,"Figures/FigSM_mae_rmse.pdf"), height=10, width=10)
{
par(mar=c(5,5,1,1))

par(mfrow=c(2,2))
plot(Va_rmse~Va_abs_bias,va, pch=19, col=scales::alpha(1,0.3), cex=0.75,ylab=expression(italic(hat(V)[A])~RMSE), xlab=expression(italic(hat(V)[A])~MAE)); abline(0,1)
text(0.8*max(va$Va_abs_bias),0.2*max(va$Va_rmse),paste("r = ",round(cor(va$Va_rmse,va$Va_abs_bias),3)))

plot(Vm_rmse~Vm_abs_bias,va, pch=19, col=scales::alpha(1,0.3), cex=0.75,ylab=expression(italic(hat(V)[M])~RMSE), xlab=expression(italic(hat(V)[M])~MAE)); abline(0,1)
text(0.8*max(va$Vm_abs_bias,na.rm=TRUE),0.2*max(va$Vm_rmse,na.rm=TRUE),paste("r = ",round(cor(va$Vm_rmse,va$Vm_abs_bias, use="complete.obs"),3)))

plot(tVa_rmse~tVa_abs_bias,va, pch=19, col=scales::alpha(1,0.3), cex=0.75, ylab=expression(italic(hat(V)[At])~RMSE), xlab=expression(italic(hat(V)[At])~MAE)); abline(0,1)
text(0.8*max(va$tVa_abs_bias),0.2*max(va$tVa_rmse),paste("r = ",round(cor(va$tVa_rmse,va$tVa_abs_bias),3)))


}
dev.off()

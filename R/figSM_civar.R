
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
		cov_amg_sim = scenarios[,"cov_amg"],
		r_amg_sim = scenarios[,"r_amg"]

# Va_Vm_cov = x[["samp_cov"]]["mother_PE","animal",],
# Va_Vm_cor = apply(x[["samp_cov"]], 3, function(y) cov2cor(y)["mother_PE","animal"]),

	)
}

all_mod<-do.call(rbind,lapply(ped_names_reduced,function(k) {

	mod4 <- do.call(rbind,lapply(get(paste0("model4_",k)),
			extract_func,k=k))

	mod5 <- do.call(rbind,lapply(get(paste0("model5_",k)), extract_func,k=k))

	rbind(
		cbind(model=4,mod4),	
		cbind(model=5,mod5)
	)
}))

all_mod$r_model <- paste0(all_mod$r,"_",all_mod$model)
all_mod$order <- match(all_mod$r_model,order)

all_mod$Va_bias <- all_mod$Va_est - all_mod$Va_sim

all_mod$Vmg_bias <- all_mod$Vmg_est - all_mod$Vmg_sim

all_mod$cov_amg_est <- ifelse(all_mod$model ==4 & is.na(all_mod$cov_amg_est),0,all_mod$cov_amg_est)

all_mod$cov_amg_bias <- all_mod$cov_amg_est - all_mod$cov_amg_sim

all_mod$r_amg_est <- all_mod$cov_amg_est/(sqrt(all_mod$Vmg_est)*sqrt(all_mod$Va_est))
all_mod$r_amg_est[all_mod$r_amg_est>1|all_mod$r_amg_est< -1]<-NA
# all_mod$r_amg_est <- ifelse(all_mod$model ==4 & is.na(all_mod$r_amg_est),0,all_mod$cor_amg_est)


# hist(all_mod$cor_amg_est[all_mod$model==5], breaks=100)
all_mod$r_amg_bias <- all_mod$r_amg_est - all_mod$r_amg_sim

# subset(all_mod,cor_amg_est>1|cor_amg_est< -1)
# $cor_amg_est

########
## MAKE SUMMARIES
########

va<- aggregate(cbind(Va_bias,Vmg_bias,cov_amg_bias,r_amg_bias) ~ scenario+r_model+r+model, all_mod,mean,na.rm=TRUE)


va$order <- match(va$r_model,order)


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

vmg$cov_amg_bias <- ifelse(vmg$model == 4,0,vmg$cov_amg_bias)
vmg$r_amg_bias <- ifelse(vmg$model == 4,0,vmg$r_amg_bias)



cols <- c(alpha(palette.colors()[1:4],0.5),palette.colors()[5])

setEPS()
pdf(paste0(wd,"Figures/FigSM_small_ped_COV.pdf"), height=10, width=8)
{par(mfrow=c(4,1), mar=c(0,5,5,1))
	beeswarm(Va_bias~ scenario2 + ped_size + model, vmg,pch=c(15:17), cex=1, pwcol=cols[as.factor(vmg$scenario3)],bg=cols,method = "compactswarm",corral="wrap", ylab=expression(Bias~"in"~italic(hat(V)[A])), labels=c("both","no_va","no_vmg"))
		abline(h=0)
	abline(v=c(3.5,6.5,9.5), col=alpha(c("grey","black","grey"),0.5))
		axis(3,c(2,5,8,11),c("Medium","Small","Medium","Small"), lwd.ticks=0, padj=1, cex.axis=1.25)

		axis(3,c(2,3.5,5),c(" ","model 3 (no cov)"," "), lwd.ticks=0, padj=1, cex.axis=1.25, line=2)
		axis(3,c(8,9.5,11),c(" ","model 4 (cov)"," "), lwd.ticks=0, padj=1, cex.axis=1.25, line=2)

par(mar=c(1,5,4,1))

	beeswarm(Vmg_bias~ scenario2 + ped_size + model, vmg,pch=c(15:17), cex=1, pwcol=cols[as.factor(vmg$scenario3)],bg=cols,method = "compactswarm",corral="wrap", ylab=expression(Bias~"in"~italic(hat(V)[Mg])), labels=c("both","no_va","no_vmg"))
		abline(h=0)
	abline(v=c(3.5,6.5,9.5), col=alpha(c("grey","black","grey"),0.5))

par(mar=c(2,5,3,1))

	beeswarm(cov_amg_bias~ scenario2 + ped_size + model, vmg,pch=c(15:17), cex=1, pwcol=cols[as.factor(vmg$scenario3)],bg=cols,method = "compactswarm",corral="wrap", ylab=expression(Bias~"in"~italic(hat(COV)['A,Mg'])), labels=c("both","no_va","no_vmg"))
		abline(h=0)
	abline(v=c(3.5,6.5,9.5), col=alpha(c("grey","black","grey"),0.5))

par(mar=c(4,5,1,1))

	beeswarm(r_amg_bias~ scenario2 + ped_size + model, vmg,pch=c(15:17), cex=1, pwcol=cols[as.factor(vmg$scenario3)],bg=cols,method = "compactswarm",corral="wrap", ylab=expression(Bias~"in"~italic(hat(r)['A,Mg'])), labels=c("both","no_va","no_vmg"))
		abline(h=0)
	abline(v=c(3.5,6.5,9.5), col=alpha(c("grey","black","grey"),0.5))

	legend("bottomleft",c("both (no COV)","both (- COV)","both (+ COV)","no_va","no_vmg"), pch=19,col=cols)
}
dev.off()

setEPS()
pdf(paste0(wd,"Figures/FigSM_small_ped_COV2.pdf"), height=5, width=15)

{
par(mfrow=c(1,3), mar=c(6,6,1,1), cex.lab=2, cex.axis=1.5)
plot(Vmg_bias~ r_amg_bias, subset(va,model==5), col=cols[as.factor(va$scenario3)], pch=c(15:16)[as.factor(va$ped_size)],cex=2)
abline(h=0,v=0)
	legend("bottomleft",c("both (no COV)","both (- COV)","both (+ COV)","no_va","no_vmg"), pch=19,col=cols)


plot(Va_bias~ r_amg_bias, subset(va,model==5), col=cols[as.factor(va$scenario3)], pch=c(15:16)[as.factor(va$ped_size)],cex=2)
abline(h=0,v=0)

plot(Vmg_bias~ Va_bias, subset(va,model==5), col=cols[as.factor(va$scenario3)], pch=c(15:16)[as.factor(va$ped_size)],cex=2)
abline(h=0,v=0)

}
dev.off()


### work out sampling covariances?

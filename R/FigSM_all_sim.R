
rm(list=ls())

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))

load(paste0(data_wd,"mge_sims3.Rdata"))
load(paste0(data_wd,"parameters.Rdata"))

library(beeswarm)
library(scales)

mat_ratio_all<-sapply(ped_str,function(x){
	rowSums(x[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(x[,-(1:2)]) 
})

mat_ratio <- colMeans(mat_ratio_all)


mod2<-do.call(rbind,lapply(ped_names,function(k) {
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
	# assign(paste0("mod2_",k),mod2)
}))
order_exp <- expand.grid(imm=c("mI","nI","bI","fI"),fec=c("lF","mF","hF"),ms=c("fs","fhs","hs"))
order<-apply(order_exp,1,function(x) paste(x[3:1],collapse="_"))

mod2$order <- match(mod2$r,order)


for(j in 1:12){
	setEPS()
	pdf(paste0(wd,"Figures/FigSM_all_sim",j,".pdf"), height=8, width=13)
{	
	par(mfrow=c(2,1), mar=c(0.5,5,5,1))
	beeswarm(Va_est~ order, subset(mod2,scenario==j),pch=19, cex=0.2, col=alpha(palette.colors()[1:4],0.5),method = "compactswarm",corral="wrap", ylab="Estimated Va", xaxt="n")
	abline(v=(1:2)*12+0.5, col=alpha("grey",0.5))
	abline(h=scenarios[j,"Va"], col=alpha(palette.colors()[8],0.5))

	for(i in c(0,12,24)){
		axis(3,c(1,2.5,4)+i,c("","Low",""), lwd.ticks=0, line=1, padj=1, cex.axis=1)
		axis(3,c(1,2.5,4)+4+i,c("","Medium",""), lwd.ticks=0, line=1, padj=1, cex.axis=1)
		axis(3,c(1,2.5,4)+8+i,c("","High",""), lwd.ticks=0, line=1, padj=1, cex.axis=1)
	}
	axis(3,c(1,6.5,12),c("","Full-Sib",""), lwd.ticks=0, line=3, padj=1, cex.axis=1)
	axis(3,c(1,6.5,12) + 12,c("","Mixed",""), lwd.ticks=0, line=3, padj=1, cex.axis=1)
	axis(3,c(1,6.5,12) + 24,c("","Half-Sib",""), lwd.ticks=0, line=3, padj=1, cex.axis=1)

	mtext("Mating System", side=3, line=-1.75, outer=TRUE, adj=0)
	mtext("Fecundity", side=3, line=-3.75, outer=TRUE, adj=0)


	par(mar=c(5,5,1,0.5))

	beeswarm(Vm_est~ order, subset(mod2,scenario==j),pch=19, cex=0.2, col=alpha(palette.colors()[1:4],0.5),method = "compactswarm",corral="wrap", ylab="Estimated Vm", labels=c("M","N","U","F"), xlab="Immigration")
	abline(v=(1:2)*12+0.5, col=alpha("grey",0.5))
	abline(h=scenarios[j,"Vmg"]+scenarios[j,"Vme"], col=alpha(palette.colors()[8],0.5))
}
	dev.off()
}


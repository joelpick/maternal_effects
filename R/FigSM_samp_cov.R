
rm(list=ls())

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))

load(paste0(data_wd,"mge_sims3.Rdata"))


# ped_names2 <- ped_names[!grepl("fI",ped_names)]

# ped_str<-ped_str[ped_names2]
# ped_str_mat<-ped_str_mat[ped_names2]

# ped_sum<-sapply(ped_str,colMeans)
# ped_sum_mat<-sapply(ped_str_mat,colMeans)

mat_ratio_all<-sapply(ped_str,function(x){
	rowSums(x[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(x[,-(1:2)]) 
})

matsib_ratio_all<-sapply(ped_str,function(x){
	rowSums(x[,c("FS","MHS")])/rowSums(x[,-(1:2)]) 
})

cov_ratio_all<-sapply(ped_str,function(x){
	rowSums(x[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS",
		"sire","PG","au_S_FS","au_S_MHS","au_D_PHS","cousin_DS_FS","cousin_DS_HS"
		)])/rowSums(x[,-(1:2)]) 
})

mat_ratio <- colMeans(mat_ratio_all)
mat_ratio_se <- apply(mat_ratio_all,2,se)

matsib_ratio <- colMeans(matsib_ratio_all)

model2_fs_lF_bI [["samp_cov"]]

samp_cov<-do.call(rbind,lapply(ped_names,function(k) {
	do.call(rbind,lapply(get(paste0("model2_",k)), function(x) {
			data.frame(
				r=k,
				scenario=1:nrow(scenarios),
				Va_Vm_cov = x[["samp_cov"]]["mother_PE","animal",],
				Va_Vm_cor = apply(x[["samp_cov"]], 3, function(y) cov2cor(y)["mother_PE","animal"]),
				Vmg_sim=scenarios[,"Vmg"]
			)
	}))
	# assign(paste0("samp_cov_",k),samp_cov)
}))

order_exp <- expand.grid(imm=c("mI","nI","bI","fI"),fec=c("lF","mF","hF"),ms=c("fs","fhs","hs"))
order<-apply(order_exp,1,function(x) paste(x[3:1],collapse="_"))

samp_cov$order <- match(samp_cov$r,order)



nrow(samp_cov)
head(samp_cov,20)
tail(samp_cov,30)
hist(samp_cov$Va_Vm_cor,breaks=50)

samp_cov[1000:1020,]



r_order2<- sapply(samp_cov$r, function(x) which(ped_names==x))

samp_cov$mat_ratio <- mat_ratio[r_order2]
samp_cov$matsib_ratio <- matsib_ratio[r_order2]

summary(lme4::lmer(Va_Vm_cor ~ mat_ratio +matsib_ratio + (1|r),samp_cov, subset=scenario=="1"))

# samp_cov$ln_Va_bias <- log(samp_cov$Va_bias)
Va_Vm_cor_mean<-aggregate(cbind(Va_Vm_cor,Vmg_sim,mat_ratio,matsib_ratio)~ scenario+r, samp_cov,mean)
Va_Vm_cor_se<-aggregate(cbind(Va_Vm_cor)~ scenario+r, samp_cov,se)
# r_order<- sapply(Va_Vm_cor_mean$r, function(x) which(ped_names==x))

# which(va1$r)

# Va_Vm_cor_mean$mat_ratio <- mat_ratio[r_order]
# Va_Vm_cor_se$mat_ratio <- mat_ratio_se[r_order]

head(Va_Vm_cor_mean)
Va_Vm_cor_mean$Vmg_n0 <- as.numeric(Va_Vm_cor_mean$Vmg_sim>0) +1

Va_Vm_cor_mean[,c("ms","fec","imm")] <- do.call(rbind,strsplit(Va_Vm_cor_mean$r,"_"))
# va2$matM_ratio2<- matM_ratio2[r_order]


# plot(matsib_ratio ~ mat_ratio,Va_Vm_cor_mean, col=Va_Vm_cor_mean$scenario, pch=19)


	
	setEPS()
	pdf(paste0(wd,"Figures/FigSM_samp_cov_mat_links.pdf"), height=6, width=13)


{	
		layout(matrix(c(1,2),nrow=1, byrow=TRUE), width=c(10,1))

par(mar=c(5,5,1,1),cex.lab=1.5, cex.axis=1.1 )

cols<-c(palette.colors(),1)#viridis::viridis(10)
	bgs= c(palette.colors(),0)
	pch=rep(21:25,2)

plot(Va_Vm_cor ~ mat_ratio,Va_Vm_cor_mean, pch= pch[Va_Vm_cor_mean$scenario],col= cols[Va_Vm_cor_mean$scenario],bg= bgs[Va_Vm_cor_mean$scenario], xlab="Proportion non-sibling maternal links", ylab="Estimated Sampling covariance")
	abline(h=0, col=alpha("grey",0.5))

par(mar=c(0,0,0,0))
	plot(NA, xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), xlab="",ylab="",bty="n")
	legend("center",LETTERS[1:10], pch=pch, pt.bg=bgs, col=cols, bty="n", cex=1.1,title="Scenario")

}
	dev.off()

# plot(Va_Vm_cor ~ matsib_ratio,Va_Vm_cor_mean, col=Va_Vm_cor_mean$scenario, pch=19)

# boxplot(Va_Vm_cor ~ scenario,Va_Vm_cor_mean)
# boxplot(Va_Vm_cor ~ r,Va_Vm_cor_mean)

	setEPS()
	pdf(paste0(wd,"Figures/FigSM_samp_cov_ped.pdf"), height=6, width=13)


{	par(mar=c(5,5,5,1),cex.lab=1.5, cex.axis=1.1)
	beeswarm(Va_Vm_cor~ order, samp_cov,pch=19, cex=0.1, col=alpha(palette.colors()[1:4],0.5),method = "compactswarm",corral="wrap", ylab="Estimated Sampling covariance", labels=c("M","N","U","F"), xlab="Immigration")
	abline(h=0, col=alpha("grey",0.5))


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

}
	dev.off()


	setEPS()
	pdf(paste0(wd,"Figures/FigSM_samp_cov_scenario.pdf"), height=6, width=13)

{
	par(mar=c(5,5,1,1),cex.lab=1.5, cex.axis=1.1)
	beeswarm::beeswarm(Va_Vm_cor~ scenario, samp_cov,pch=19, cex=0.1, col=scales::alpha(c(1,1,2,rep(1,7)),0.3),method = "compactswarm",corral="wrap",labels=LETTERS[1:10],ylab="Estimated Sampling covariance", xlab="Scenario")
	abline(h=0, col=alpha("grey",0.5))
	}	
	dev.off()
# beeswarm::beeswarm(Va_Vm_cov~ scenario, samp_cov,pch=19, cex=0.1, col=scales::alpha(c(1,1,2,rep(1,7)),0.5),method = "compactswarm",corral="wrap")

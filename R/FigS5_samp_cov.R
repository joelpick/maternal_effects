
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

samp_cov<-do.call(rbind,lapply(ped_names,function(k) {
	do.call(rbind,lapply(get(paste0("model2_",k)), function(x) {
			data.frame(
				r=k,
				scenario=1:nrow(scenarios),
				Va_Vm = apply(x[["samp_cov"]], 3, function(y) cov2cor(y)["mother_PE","animal"]),
				Vmg_sim=scenarios[,"Vmg"]
			)
	}))
	# assign(paste0("samp_cov_",k),samp_cov)
}))

nrow(samp_cov)
head(samp_cov,20)
hist(samp_cov$Va_Vm,breaks=50)

r_order2<- sapply(samp_cov$r, function(x) which(ped_names==x))

samp_cov$mat_ratio <- mat_ratio[r_order2]
samp_cov$matsib_ratio <- matsib_ratio[r_order2]

summary(lme4::lmer(Va_Vm ~ mat_ratio +matsib_ratio + (1|r),samp_cov, subset=scenario=="1"))

# samp_cov$ln_Va_bias <- log(samp_cov$Va_bias)
Va_Vm_mean<-aggregate(cbind(Va_Vm,Vmg_sim,mat_ratio,matsib_ratio)~ scenario+r, samp_cov,mean)
Va_Vm_se<-aggregate(cbind(Va_Vm)~ scenario+r, samp_cov,se)
# r_order<- sapply(Va_Vm_mean$r, function(x) which(ped_names==x))

# which(va1$r)

# Va_Vm_mean$mat_ratio <- mat_ratio[r_order]
# Va_Vm_se$mat_ratio <- mat_ratio_se[r_order]

head(Va_Vm_mean)
Va_Vm_mean$Vmg_n0 <- as.numeric(Va_Vm_mean$Vmg_sim>0) +1

Va_Vm_mean[,c("ms","fec","imm")] <- do.call(rbind,strsplit(Va_Vm_mean$r,"_"))
# va2$matM_ratio2<- matM_ratio2[r_order]



plot(Va_Vm ~ mat_ratio,Va_Vm_mean, col=Va_Vm_mean$scenario, pch=19)
plot(Va_Vm ~ matsib_ratio,Va_Vm_mean, col=Va_Vm_mean$scenario, pch=19)

boxplot(Va_Vm ~ scenario,Va_Vm_mean)

boxplot(Va_Vm ~ r,Va_Vm_mean)


beeswarm::beeswarm(Va_Vm~ scenario, samp_cov,pch=19, cex=0.1, col=scales::alpha(c(1,1,2,rep(1,7)),0.5),method = "compactswarm",corral="wrap")
abline(h=0)
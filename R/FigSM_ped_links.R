
rm(list=ls())

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))
load(paste0(data_wd,"mge_sims3.Rdata"))	


mat_ratio_all<-sapply(ped_str,function(x){
	rowSums(x[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(x[,-(1:2)])
})
mat_ratio <- colMeans(mat_ratio_all)


ratio_dat <- as.data.frame(do.call(rbind,strsplit(names(mat_ratio),"_")))
names(ratio_dat) <- c("ms","fec","imm")
ratio_dat$nsml <- mat_ratio
ratio_dat$total <- rowSums(t(sapply(ped_str,colMeans))[,-(1:2)])

ratio_dat <- cbind(ratio_dat,t(sapply(ped_str,colMeans)))


setEPS()
pdf(paste0(wd,"Figures/FigSM_links_ped.pdf"), height=10, width=5)
{
par(mfrow=c(3,1), mar=c(4,4,1,1), cex.lab=1.4,mgp=c(2,0.5,0))

boxplot(nsml ~ ms,ratio_dat, at=c(2,3,1), names=c("Mixed", "Full-sibs", "Half-sibs"), ylab="Proportion non-sibling maternal links", xlab="Mating system")
boxplot(nsml ~ fec,ratio_dat, at=c(1,3,2), names=c("Low","High", "Mid"), ylab="Proportion non-sibling maternal links", xlab="Fecundity")
boxplot(nsml ~ imm,ratio_dat, at=c(3,2,4,1), names=c("Both sexes", "Female", "Male", "None"), ylab="Proportion non-sibling maternal links", xlab="Dispersal")
}
dev.off()





### dam-offspring are fixed in all sims

{
par(mfrow=c(3,2), mar=c(4,4,1,1), cex.lab=1.4,mgp=c(2,0.5,0))
boxplot(nsml ~ fec,ratio_dat, at=c(1,3,2), names=c("Low","High", "Mid"), ylab="Prop. relationships", xlab="Fecundity")
boxplot(FS + MHS ~ fec,ratio_dat, at=c(1,3,2), names=c("Low","High", "Mid"), ylab="Maternal siblings", xlab="Fecundity")
boxplot(total ~ fec,ratio_dat, at=c(1,3,2), names=c("Low","High", "Mid"), ylab="total links", xlab="Fecundity")
boxplot(MG ~ fec,ratio_dat, at=c(1,3,2), names=c("Low","High", "Mid"), ylab="grandmaternal-offspring", xlab="Fecundity")
boxplot(au_D_FS + au_D_MHS ~ fec,ratio_dat, at=c(1,3,2), names=c("Low","High", "Mid"), ylab="Maternal aunt/uncle", xlab="Fecundity")
boxplot(cousin_D_FS+cousin_D_HS ~ fec,ratio_dat, at=c(1,3,2), names=c("Low","High", "Mid"), ylab="maternal cousins", xlab="Fecundity")
}

{
par(mfrow=c(3,2), mar=c(4,4,1,1), cex.lab=1.4,mgp=c(2,0.5,0))
boxplot(nsml ~ ms,ratio_dat, at=c(2,3,1), names=c("Mixed", "Full-sibs", "Half-sibs"), ylab="Prop. relationships", xlab="Mating system")
boxplot(FS + MHS ~ ms,ratio_dat, at=c(2,3,1), names=c("Mixed", "Full-sibs", "Half-sibs"), ylab="Maternal siblings", xlab="Mating system")
boxplot(total ~ ms,ratio_dat, at=c(2,3,1), names=c("Mixed", "Full-sibs", "Half-sibs"), ylab="total links", xlab="Mating system")
boxplot(MG ~ ms,ratio_dat, at=c(2,3,1), names=c("Mixed", "Full-sibs", "Half-sibs"), ylab="grandmaternal-offspring", xlab="Mating system")
boxplot(au_D_FS + au_D_MHS ~ ms,ratio_dat, at=c(2,3,1), names=c("Mixed", "Full-sibs", "Half-sibs"), ylab="Maternal aunt/uncle", xlab="Mating system")
boxplot(cousin_D_FS+cousin_D_HS ~ ms,ratio_dat, at=c(2,3,1), names=c("Mixed", "Full-sibs", "Half-sibs"), ylab="maternal cousins", xlab="Mating system")
}

# boxplot(cousin_D_FS ~ ms,ratio_dat, at=c(2,3,1), names=c("Mixed", "Full-sibs", "Half-sibs"), ylab="maternal cousins", xlab="Mating system")
# boxplot(cousin_D_HS ~ ms,ratio_dat, at=c(2,3,1), names=c("Mixed", "Full-sibs", "Half-sibs"), ylab="maternal cousins", xlab="Mating system")


{
par(mfrow=c(3,2), mar=c(4,4,1,1), cex.lab=1.4,mgp=c(2,0.5,0))
boxplot(nsml ~ imm,ratio_dat, at=c(3,2,4,1), names=c("Both sexes", "Female", "Male", "None"), ylab="Prop. relationships", xlab="Dispersal")
boxplot(FS + MHS ~ imm,ratio_dat, at=c(3,2,4,1), names=c("Both sexes", "Female", "Male", "None"), ylab="Maternal siblings", xlab="Dispersal")
boxplot(total ~ imm,ratio_dat, at=c(3,2,4,1), names=c("Both sexes", "Female", "Male", "None"), ylab="total links", xlab="Dispersal")
boxplot(MG ~ imm,ratio_dat, at=c(3,2,4,1), names=c("Both sexes", "Female", "Male", "None"), ylab="grandmaternal-offspring", xlab="Dispersal")
boxplot(au_D_FS + au_D_MHS ~ imm,ratio_dat, at=c(3,2,4,1), names=c("Both sexes", "Female", "Male", "None"), ylab="Maternal aunt/uncle", xlab="Dispersal")
boxplot(cousin_D_FS+cousin_D_HS ~ imm,ratio_dat, at=c(3,2,4,1), names=c("Both sexes", "Female", "Male", "None"), ylab="maternal cousins", xlab="Dispersal")
}
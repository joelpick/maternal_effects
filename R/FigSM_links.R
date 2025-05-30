
rm(list=ls())

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

load(paste0(data_wd,"mge_sims3.Rdata"))


mat_ratio_all<-sapply(ped_str,function(x){
	rowSums(x[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(x[,-(1:2)]) 
})

cov_ratio_all<-sapply(ped_str,function(x){
	rowSums(x[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS",
		"sire","PG","au_S_FS","au_S_MHS","au_D_PHS","cousin_DS_FS","cousin_DS_HS"
		)])/rowSums(x[,-(1:2)]) 
})

mat_ratio <- colMeans(mat_ratio_all)

cov_ratio <- colMeans(cov_ratio_all)



setEPS()
pdf(paste0(wd,"Figures/FigSM_links_cor.pdf"), height=5, width=5)
par(mar=c(6,6,1,1))
plot(cov_ratio~mat_ratio, pch=19, xlab=expression(Proportion~non-sibling~maternal~links~(italic(V[Mg]))), ylab=expression('Proportion non-sibling maternal links \nthrough single mother'~(italic(COV["A,Mg"]))))
dev.off()



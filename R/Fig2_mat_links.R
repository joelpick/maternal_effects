

rm(list=ls())

library(parallel)

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))



load(paste0(data_wd,"mge_sims3.Rdata"))
mat_ratio_all<-sapply(ped_str,function(x){
	rowSums(x[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(x[,-(1:2)])
})
mat_ratio <- colMeans(mat_ratio_all)


ped_BT<-nadiv::prepPed(read.csv(paste0(wd,"Data/Raw/ped_BT.csv"))[,1:3])
ped_RD<-read.csv(paste0(wd,"Data/Raw/ped_RD.csv"))[,1:3]
ped_SFW<-read.csv(paste0(wd,"Data/Raw/ped_SFW.csv"))[,1:3]
ped_SSH<-read.csv(paste0(wd,"Data/Raw/ped_SSH.csv"))[,1:3]
ped_SV<-read.csv(paste0(wd,"Data/Raw/ped_SV.csv"))[,1:3]

chick_stat <- function(ped) ped_stat(ped, ped[!is.na(ped$dam),1])
adult_stat <- function(ped) ped_stat(ped, unique(c(ped$dam,ped$sire)))
adultM_stat <- function(ped) ped_stat(ped, unique(ped$sire))
adultF_stat <- function(ped) ped_stat(ped, unique(ped$dam))



# stat1<-rbind(
# 	BT_chick=chick_stat(ped_BT),
# 	BT_adult=adult_stat(ped_BT),
# 	RD_chick=chick_stat(ped_RD),
# 	RD_adult=adult_stat(ped_RD),
# 	SFW_chick=chick_stat(ped_SFW),
# 	SFW_adult=adult_stat(ped_SFW),
# 	SSH_chick=chick_stat(ped_SSH),
# 	SSH_adult=adult_stat(ped_SSH),
# 	SV_chick=chick_stat(ped_SV),
# 	SV_adult=adult_stat(ped_SV)
# )


stat1<-rbind(
	BT_chick=chick_stat(ped_BT),
	BT_adult=adult_stat(ped_BT),
	# BT_adultM=adultM_stat(ped_BT),
	# BT_adultF=adult_stat(ped_BT),
	RD_chick=chick_stat(ped_RD),
	RD_adult=adult_stat(ped_RD)
	# RD_adultM=adultM_stat(ped_RD),
	# RD_adultF=adult_stat(ped_RD)
)

stat_mr<-rowSums(stat1[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(stat1[,-(1:2)])



setEPS()
pdf(paste0(wd,"Figures/fig2_mat_links.pdf"), height=6, width=8)
# order <- c(17,13,16,14,15)-5
order <- c(10,11)
{
par(mar=c(5,1,1,1), cex.lab=1.25)
hist(mat_ratio, xlim=c(0,0.5), breaks=15, ylim=c(0,14), yaxt="n", xlab="Proportion non-sibling maternal links", ylab=
	"",main="")
arrows(stat_mr[c(1,3,5,7,9)],order,stat_mr[c(2,4,6,8,10)],order,code=0)
points(stat_mr, rep(order,each=2), pch=19, col=c("blue","red"))
# text(0.025,order,c("Blue tit","Red deer","Superb fairy wren","Soay sheep","Snow vole"))
text(0.025,order,c("Blue tit","Red deer"))
# text(stat_mr[1:2],c(15,15),c("Juvenile","Adult"),col=c("blue","red"))

legend("top",c("Juvenile","Adult"),pch=19,col=c("blue","red"),bty="n")
}
dev.off()


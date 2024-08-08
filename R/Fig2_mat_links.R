

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


ped_BT<-read.csv(paste0(wd,"Data/Raw/BT_Pick.csv"))[,1:3]
ped_RD<-read.table(paste0(wd,"Data/Raw/RD_Gauzere.txt"),header=TRUE)[,c(1,3,2)]
names(ped_RD) <- c("animal","dam","sire")
nrow(ped_BT)
nrow(ped_RD)

juv_stat <- function(ped) ped_stat(ped, ped[!is.na(ped$dam),1])
adult_stat <- function(ped) ped_stat(ped, unique(c(ped$dam,ped$sire)))


stat1<-rbind(
	BT_juv=juv_stat(ped_BT),
	BT_adult=adult_stat(ped_BT),
	RD_juv=juv_stat(ped_RD),
	RD_adult=adult_stat(ped_RD)
)

stat_mr<-rowSums(stat1[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(stat1[,-(1:2)])


setEPS()
pdf(paste0(wd,"Figures/fig2_mat_links.pdf"), height=6, width=8)
order <- c(7,8)
{
par(mar=c(5,1,1,1), cex.lab=1.5, cex.axis=1.25)
hist(mat_ratio, xlim=c(0,0.5), breaks=15, ylim=c(0,11), yaxt="n", xlab="Proportion non-sibling maternal links", ylab=
	"",main="")
arrows(stat_mr[c(1,3)],order,stat_mr[c(2,4,6,8,10)],order,code=0)
points(stat_mr, rep(order,each=2), pch=19, col=c("blue","red"), cex=1.5)
# text(0.025,order,c("Blue tit","Red deer","Superb fairy wren","Soay sheep","Snow vole"))
text(0.025,order,c("Blue tit","Red deer"), cex=1.5)
# text(stat_mr[1:2],c(15,15),c("Juvenile","Adult"),col=c("blue","red"))

legend("top",c("Juvenile","Adult"),pch=19,col=c("blue","red"),bty="n", cex=1.5)
}
dev.off()



rm(list=ls())

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
pdf(paste0(wd,"Figures/Fig3_mat_links.pdf"), height=4, width=7)
{
par(mfrow=c(1,1),mar=c(5,1,1,1), cex.lab=1.5, cex.axis=1.25)
hist(mat_ratio, xlim=c(0,0.5), breaks=15, ylim=c(0,8), yaxt="n", xlab="Proportion non-sibling maternal links", ylab=
	"",main="")

arrows(stat_mr,rep(6,4),stat_mr,rep(0,4),code=2,lwd=4,col=c("blue","blue","red","red"), length=0.15, lty=c(1,6))

legend("topleft",c("Juvenile","Adult","Blue tit","Red deer"),lwd=2,lty=c(1,6,0,0),pch=c(NA,NA,19,19),col=c(1,1,"blue","red"),bty="n", cex=1.1)

}
dev.off()


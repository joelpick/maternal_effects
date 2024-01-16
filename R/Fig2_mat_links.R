

rm(list=ls())

library(parallel)

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/extract_cousins.R"))
source(paste0(wd,"R/00_functions.R"))
source("/Users/joelpick/github/squidPed/R/simulate_pedigree.R")

max_gen <- 10
n_sims <- 50
cores<-6
fecundity=8

immigration <- rbind(
	no_immigration = c(0,0),
	f_immigration=c(0.4,0.1),
	m_immigration=c(0.1,0.4),
	b_immigration=c(0.25,0.25))

combos<-expand.grid(gens=1:10,immigration=1:4)

peds_param <- cbind(
	years= combos[,"gens"],
	immigration_f = immigration[combos[,"immigration"],1],
	immigration_m = immigration[combos[,"immigration"],2])

peds_param <- cbind(peds_param, 
	juv_surv_f=2*(1 - peds_param[,"immigration_f"])/fecundity,
	juv_surv_m=2*(1 - peds_param[,"immigration_m"])/fecundity
)

# ped_names <- rownames(peds_param) <- paste(
# 	ms_names[combos[,"p_sire"]],
# 	f_names[combos[,"fecundity"]],
# 	i_names[combos[,"immigration"]], sep="_")


ped_str_pd <- vector("list",length=nrow(peds_param))

for(k in 1:nrow(peds_param)){
	cat(k, " ")
	ped_str_pd[[k]] <- do.call(rbind,mclapply(1:n_sims,	function(i){
		ped <- simulate_pedigree(
			years = peds_param[k,"years"],
			n_females = 100,
			fecundity = fecundity,
			p_sire = 0.75,
			p_polyandry=1,
			juv_surv = c(peds_param[k,"juv_surv_f"],peds_param[k,"juv_surv_m"]),
			adult_surv = 0,					# discrete generations
			immigration = c(peds_param[k,"immigration_f"],peds_param[k,"immigration_m"]), 				# closed population
			constant_pop = TRUE     # constant population size
			)$pedigree
		ped_stat(ped)
	}, mc.cores=cores))		
}

mat_ratio_all_pd<-sapply(ped_str_pd,function(x){
	rowSums(x[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(x[,-(1:2)])
})

mat_ratio_pd <- colMeans(mat_ratio_all_pd)




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

stat1<-rbind(
	BT_chick=chick_stat(ped_BT),
	BT_adult=adult_stat(ped_BT),
	RD_chick=chick_stat(ped_RD),
	RD_adult=adult_stat(ped_RD),
	SFW_chick=chick_stat(ped_SFW),
	SFW_adult=adult_stat(ped_SFW),
	SSH_chick=chick_stat(ped_SSH),
	SSH_adult=adult_stat(ped_SSH),
	SV_chick=chick_stat(ped_SV),
	SV_adult=adult_stat(ped_SV)
)

stat_mr<-rowSums(stat1[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(stat1[,-(1:2)])



setEPS()
pdf(paste0(wd,"Figures/fig2_mat_links.pdf"), height=6, width=15)
order <- c(17,13,16,14,15)-5
{
par(mfrow=c(1,2), mar=c(5,1,1,1), cex.lab=1.25)
hist(mat_ratio, xlim=c(0,0.5), breaks=15, ylim=c(0,14), yaxt="n", xlab="Proportion non-sibling maternal links", ylab=
	"",main="")
arrows(stat_mr[c(1,3,5,7,9)],order,stat_mr[c(2,4,6,8,10)],order,code=0)
points(stat_mr, rep(order,each=2), pch=19, col=c("blue","red"))
text(0.025,order,c("Blue tit","Red deer","Superb fairy wren","Soay sheep","Snow vole"))
# text(stat_mr[1:2],c(15,15),c("Juvenile","Adult"),col=c("blue","red"))

legend("top",c("Juvenile","Adult"),pch=19,col=c("blue","red"),bty="n")


	par(mar=c(5,5,1,1))
plot(mat_ratio_pd~rep(1:10,4), pch=19, col=rep(viridis::viridis(4),each=10),ylim=c(0,0.4), xlab="Pedigree Depth (Generations)", ylab="Proportion of non-sibling maternal links")
for(i in 1:4) lines(mat_ratio_pd[(1:10)+10*(i-1)]~c(1:10),col=viridis::viridis(4)[i])

legend("bottomright",c("None","Female biased","Male biased","No bias"), col=viridis::viridis(4), pch=19, title="Immigration", bty="n")
}
dev.off()


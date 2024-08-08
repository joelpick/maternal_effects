

rm(list=ls())

library(parallel)
library(squidPed)

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))


max_gen <- 10
n_sims <- 50
cores<-6
fecundity=6

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
		ped_stat(ped,phenotyped=ped[!is.na(ped[,"dam"]),"animal"])
	}, mc.cores=cores))		
}

mat_ratio_all_pd<-sapply(ped_str_pd,function(x){
	rowSums(x[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(x[,-(1:2)])
})

mat_ratio_pd <- colMeans(mat_ratio_all_pd)
# mat_ratio_pd_sd <- apply(mat_ratio_all_pd,2,sd)




setEPS()
pdf(paste0(wd,"Figures/FigSM_ped_depth.pdf"), height=5, width=7)
{
	par(mar=c(5,5,1,1),cex.lab=1.3, cex.axis=1.1)
plot(mat_ratio_pd~rep(1:10,4), pch=19, col=rep(viridis::viridis(4),each=10),ylim=c(0,0.4), xlab="Pedigree Depth (Generations)", ylab="Proportion of non-sibling maternal links")
for(i in 1:4) lines(mat_ratio_pd[(1:10)+10*(i-1)]~c(1:10),col=viridis::viridis(4)[i])

legend("bottomright",c("None","Female biased","Male biased","No bias"), col=viridis::viridis(4), pch=19, title="Immigration", bty="n")
}
dev.off()
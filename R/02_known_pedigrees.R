rm(list=ls())
wd <- "/Users/joelpick/github/maternal_effects/"
data_wd <- paste0(wd,"Data/Intermediate/")
source(paste0(wd,"R/extract_cousins.R"))
source(paste0(wd,"R/00_functions.R"))

# ped_bt <-  MasterBayes::insertPed(read.csv(paste0(wd,"Data/Raw/ped_BT.csv"))[,1:3])



## extract data from a couple of additional pedigrees
ped_BT <-  MasterBayes::insertPed(read.csv(paste0(wd,"Data/Raw/ped_BT.csv"))[,1:3])
#
 ped_BT_social <- MasterBayes::insertPed(read.csv(paste0(wd,"Data/Raw/ped_BT.csv"))[,c(1,2,4)])
head(ped_BT_social)
 ped_BT_socialD <- MasterBayes::insertPed(read.csv(paste0(wd,"Data/Raw/ped_BT.csv"))[,c(1,2,5)])
 colnames(ped_BT_social)[3] <- colnames(ped_BT_socialD)[3] <- "sire"



dat_bt<-subset(read.csv("/Users/joelpick/Dropbox/0_blue_tits/skew/Data/Raw/tPED.csv"), ped_type=="DSCM")[,1:4]
ped_bt<-pedantics::fixPedigree(dat_bt[,1:3])

ped_bt <- ped_BT_social

ped_bt[,"cohort"] <- dat_bt[match(ped_bt[,1],dat_bt[,1]),4]
ped_bt[,4] <- as.numeric(substr(ped_bt[,4],1,2))
names(ped_bt)[1] <- "animal"

head(ped_bt)

# ped_rd<-read.csv(paste0(wd,"Data/Raw/ped_RD.csv"))
ped_rd<-pedantics::fixPedigree(read.table("/Users/joelpick/github/pedigree_simulations/Data/Raw/red_deer/Pedigree_File.txt", header=TRUE))
names(ped_rd) <- c("animal","sire","dam")
dat_rd<-read.table("/Users/joelpick/github/pedigree_simulations/Data/Raw/red_deer/Individual_Data.txt", header=TRUE)

ped_rd$cohort <- dat_rd[match(ped_rd[,"animal"],dat_rd[,"ID"]),"BirthYear"]
table(ped_rd$cohort)

ped_rd_short <- pedantics::fixPedigree(subset(ped_rd,cohort<2005))




par(mfrow=c(2,1))
plot(table(table(ped_bt$dam)))
mean(table(ped_bt$dam))

plot(table(table(ped_rd$dam)))
mean(table(ped_rd$dam))




mat_ratio2<-
ped_sum2[2,]/colSums(ped_sum2)

ped_sum2[2,]/ped_sum2[3,]

## blue tit pedigree - there will be a lot of unknown females that abandon on eggs, but will be in ped through x-fostering



ped<-ped_bt
# immigration <- function(ped, sex_specific=TRUE){
	founders <- ped$animal[is.na(ped$sire) & is.na(ped$sire)]
	ped$dam_founder <- ped$dam%in%founders
	ped$sire_founder <- ped$sire%in%founders
	# year_d <- aggregate(cohort~dam,subset(ped,!dam%in%founders),min)$cohort
	# year_d_F <- aggregate(cohort~dam,subset(ped,dam%in%founders),min)$cohort
	# year_s <- aggregate(cohort~sire,subset(ped,!sire%in%founders),min)$cohort
	# year_s_F <- aggregate(cohort~sire,subset(ped,sire%in%founders),min)$cohort

	year_d <- aggregate(cohort~dam+dam_founder,ped,min)
	year_s <- aggregate(cohort~sire+sire_founder,ped,min)

	all_d <- table(year_d$cohort,year_d$dam_founder)
	all_s <- table(year_s$cohort,year_s$sire_founder)
colMeans(	all_d)
colMeans(	all_s)
mean(all_d[,2]/rowSums(all_d))
mean(all_s[,2]/rowSums(all_s))

plot(	all_d[,2]/rowSums(all_d)~rownames(all_d));abline(h=0.05)
plot(	all_s[,2]/rowSums(all_s)~rownames(all_s));abline(h=0.20)

plot(	all_d[,2]/rowSums(all_d)~rownames(all_d), col="red", ylim=c(0,1))
points(	all_s[,2]/rowSums(all_s)~rownames(all_s), col="blue")

plot(	all_d[,2]~rownames(all_d), col="red", ylim=range(c(all_d,all_s)))
points(	all_s[,2]~rownames(all_s), col="blue")
points(	all_d[,1]~rownames(all_s), col="red", pch=19)
points(	all_s[,1]~rownames(all_s), col="blue", pch=19)


# 	all_s <- table(year_s)
# 	all_d_F <- table(year_d_F)
# 	all_s_F <- table(year_s_F)


# 	print(all_d)
# 	print(all_s)

# 	c(
# 		f=mean(all_d[-1]/all_d[1]),
# 		m=mean(all_s[-1]/all_s[1])
# 		)
# }






immigration(ped_bt)

plot(as.numeric(table(ped_rd$cohort)))

plot(aggregate(dam~cohort,ped_rd,function(x) length(unique(x)))$dam)
plot(aggregate(sire~cohort,ped_rd,function(x) length(unique(x)))$sire)





## range offf offspring per mother from bonnet et al 2022

main_dir <- "/Users/joelpick/github/maternal_effects/Data/Raw/Bonnet/"
mean_offspring<-	sapply(dir(main_dir), function(i){
	dat<-read.csv(paste0(main_dir,i,"/data_",i,".csv"))
	#table(table(dat$dam))
	mean(table(dat$dam))
})
mean(mean_offspring)













rm(list=ls())


library(asreml)
library(parallel)

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/extract_cousins.R"))
source(paste0(wd,"R/00_functions.R"))
source("/Users/joelpick/github/squidPed/R/simulate_pedigree.R")
# devtools::load_all("~/github/squidSim/R")

run=FALSE

n_sims <-100
generations=5
n_females=100
fecundity=10

n_females=n_females
fecundity=fecundity
p_sire = 1
p_polyandry=0.5
juv_surv = rep(2/fecundity,2)
immigration = c(0,0)


ped <- simulate_pedigree(
	years = generations,
	n_females = n_females,
	fecundity = fecundity,
	p_sire = p_sire,
	p_polyandry=1,
	juv_surv = juv_surv,
	adult_surv = 0,					# discrete generations
	immigration = immigration, 				# closed population
	constant_pop = TRUE     # constant population size
)$pedigree







############

scenarios <- rbind(	
	# # C) Maternal genetic only
	# c=c(Va=0, Vmg=0.25, r_amg=0, Vme=0),
	# # D) Direct genetic and maternal environment
	# e=c(Va=0, Vmg=0.25, r_amg=0, Vme=0.25),#####
	# # F) Direct and maternal genetic, no covariance
	# f=c(Va=0.25, Vmg=0.25, r_amg=0, Vme=0),
	# # I) Direct and maternal genetic, no covariance and maternal environment
	# i=c(Va=0.25, Vmg=0.25, r_amg=0, Vme=0.25),
	rd=c(Va=0.05, Vmg=0.4, r_amg=0, Vme=0.05)
)
ped_names <- c("rd","rd_short","bt")

ped_str_known <- cbind(
	rd=ped_stat(ped_rd),
	rd_short=ped_stat(ped_rd_short),
	bt=ped_stat(ped_bt))

ped_sum2 <- rbind(mat_sibs=colSums(ped_str_known[c("FS","MHS"),]), mat_links=colSums(ped_str_known[c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS"),]), other=colSums(ped_str_known[!rownames(ped_str_known)%in%c( "individuals" ,"links" ,"FS","MHS","dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS"),]) )
ped_sum2[2,]/colSums(ped_sum2)


	cat("\n\nSimulating Data: \n")
	## simulate data
	for(k in ped_names){
		dat<-mclapply(1:100, function(i){
			x<-vector("list", nrow(scenarios))
			for(j in 1:nrow(scenarios)){
				x[[j]]<- mge_sim(get(paste0("ped_",k))[,1:3], param=scenarios[j,])
			}
			x
		}, mc.cores=8)
		assign(paste0(k,"_data"),dat)
		cat(k, " ")
	}

	## run models
	cat("\n\nRunning models: \n")
	for(k in ped_names){
		cat(k, "\n")
		cat("Model 1\n")
		model1 <- model_func(m1_func,get(paste0("ped_",k)),get(paste0(k,"_data")),mc.cores=8)
		assign(paste0("model1_",k),model1)
		cat("\nModel 2\n")
		model2 <- model_func(m2_func,get(paste0("ped_",k)),get(paste0(k,"_data")),mc.cores=8)
		assign(paste0("model2_",k),model2)
		cat("\n")
	}

	save(list=(c("ped_str_known",paste0("model1_",ped_names),paste0("model2_",ped_names))),file=paste0(data_wd,"mge_sims_known_rd.Rdata"))






mod2<-do.call(rbind,lapply(ped_names,function(k) {
	mod2 <- do.call(rbind,lapply(get(paste0("model2_",k)), function(x) {
		do.call(rbind,lapply(1:nrow(scenarios), function(i) data.frame(r=k,scenario=i,Va_est = x[i,"A"],Vm_est = x[i,"Me"],Va_sim=scenarios[i,"Va"],Vm_sim =sum(scenarios[i,c("Vmg","Vme")]),Vmg_sim=scenarios[i,"Vmg"])))
	}))
	# assign(paste0("mod2_",k),mod2)
}))
#,sum(x[i,c("Mg","Me")])
head(mod2,20)
mod2$Va_bias <- mod2$Va_est - mod2$Va_sim
# mod2$ln_Va_bias <- log(mod2$Va_bias)
va2<-aggregate(cbind(Va_bias,Vmg_sim,Vm_sim)~ scenario+r, mod2,mean)







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
SV_adult=adult_stat(ped_SV))

stat_mr<-rowSums(stat1[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(stat1[,-(1:2)])

order <- c(17,15,16,14,13)
{
par(mfrow=c(1,1), mar=c(5,1,1,1))
hist(mat_ratio, xlim=c(0,0.5), breaks=10, ylim=c(0,20), yaxt="n", ylab=
	"",main="")
points(stat_mr, rep(order,each=2), pch=19, col=c("blue","red"))
arrows(stat_mr[c(1,3,5,7,9)],order,stat_mr[c(2,4,6,8,10)],order)
text(0.05,order,c("Blue tit","Red deer","Superb Fairy Wren","Soay Sheep","Snow Voles"))
}



stat2<-rbind(RD=ped_stat2(ped_RD),
SFW=ped_stat2(ped_SFW),
SSH=ped_stat2(ped_SSH),
SV=ped_stat2(ped_SV))

(stat2[,2]-stat2[,1])/stat2[,3]


stat1<-rbind(BT=ped_stat(ped_BT),
	RD=ped_stat(ped_RD),
SFW=ped_stat(ped_SFW),
SSH=ped_stat(ped_SSH),
SV=ped_stat(ped_SV))






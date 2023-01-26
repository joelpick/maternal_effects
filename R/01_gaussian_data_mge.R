
rm(list=ls())


library(asreml)
library(parallel)

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/extract_cousins.R"))
source(paste0(wd,"R/00_functions.R"))
source("/Users/joelpick/github/squidPed/R/simulate_pedigree.R")
# devtools::load_all("~/github/squidSim/R")

run=TRUE

n_sims <-100
generations=5
n_females=100

## Range from bonnet et al 
m_fecundity = 6 # mean is 6.9
h_fecundity = 12 #highest is 13.3
l_fecundity = 3 #lowest is 1.5
## these numbers are easier to maintain total number of phenotuypes indidivuals

## first number is female, then male immigration rate
no_immigration = c(0,0)
f_immigration=c(0.5,0)
m_immigration=c(0,0.5)
b_immigration=c(0.25,0.25) # same overall immigration but no sex bias

fhs = 0.75
fs = 1
hs = 0

#### are number of maternal links the same in the two pedigrees?
#### number of maternal half sibs will differ presumably
#### moore et al meta-analysis mentions taking into account pedigree links in the supp matt

## create 100 full sib pedigrees 

peds_param <- rbind(
	##
baseline = c(n_females=n_females,fecundity=m_fecundity,p_sire = fhs,immigration = no_immigration),

fs=c(n_females=n_females,fecundity=m_fecundity,p_sire = fs,immigration = no_immigration),
hs=c(n_females=n_females,fecundity=m_fecundity,p_sire = hs,immigration = no_immigration),

hF=c(n_females = n_females/2,fecundity = h_fecundity,p_sire = fhs,immigration = no_immigration),
lF=c(n_females = n_females*2,fecundity = l_fecundity,p_sire = fhs,immigration = no_immigration),

IF=c(n_females = n_females, fecundity = m_fecundity, p_sire = fhs, immigration = f_immigration),
IM=c(n_females = n_females, fecundity = m_fecundity, p_sire = fhs, immigration = m_immigration),
IB=c(n_females = n_females, fecundity = m_fecundity, p_sire = fhs, immigration = b_immigration)
)

# 1=(juv_surv_f * fecundity)/2 + immigration_f
# juv_surv_f = 2*(1 - immigration_f)/fecundity

peds_param <- cbind(peds_param, 
	juv_surv1=2*(1 - peds_param[,"immigration1"])/peds_param[,"fecundity"],
	juv_surv2=2*(1 - peds_param[,"immigration2"])/peds_param[,"fecundity"]
)

# peds_param[,"n_females"]*peds_param[,"juv_surv1"]*peds_param[,"fecundity"]
# peds_param[,"n_females"]*peds_param[,"juv_surv2"]*peds_param[,"fecundity"]

ped_names <- rownames(peds_param)


scenarios <- rbind(	
	# C) Maternal genetic only
	c=c(Va=0, Vmg=0.25, r_amg=0, Vme=0),
	# D) Direct genetic and maternal environment
	e=c(Va=0, Vmg=0.25, r_amg=0, Vme=0.25),#####
	# F) Direct and maternal genetic, no covariance
	f=c(Va=0.25, Vmg=0.25, r_amg=0, Vme=0),
	# I) Direct and maternal genetic, no covariance and maternal environment
	i=c(Va=0.25, Vmg=0.25, r_amg=0, Vme=0.25)
)

if(run){

	set.seed(20230126)

	ped_str <- vector("list",length=nrow(peds_param))

	## make pedigrees
	cat("Simulating Pedigrees:\n")
	for(j in 1:nrow(peds_param)){
		peds <- mclapply(1:n_sims,	function(i){
			simulate_pedigree(
				years = generations,
				n_females = peds_param[j,"n_females"],
				fecundity = peds_param[j,"fecundity"],
				p_sire = peds_param[j,"p_sire"],
				p_polyandry=1,
				juv_surv = c(peds_param[j,"juv_surv1"],peds_param[j,"juv_surv2"]),
				adult_surv = 0,					# discrete generations
				immigration = c(peds_param[j,"immigration1"],peds_param[j,"immigration2"]), 				# closed population
				constant_pop = TRUE     # constant population size
				)$pedigree
		}, mc.cores=8)
		assign(paste0(ped_names[j] ,"_peds"),peds)	
		ped_str[[ped_names[j]]]<- do.call(rbind,mclapply(peds,ped_stat, mc.cores=8))
		cat(j, " ")
	}


lapply(ped_str,colMeans)

	cat("\nSimulating Data: \n")
	## simulate data
	for(k in ped_names){
		dat<-mclapply(get(paste0(k,"_peds")), function(i){
			x<-vector("list", nrow(scenarios))
			for(j in 1:nrow(scenarios)){
				x[[j]]<- mge_sim(i[,1:3], param=scenarios[j,])
			}
			x
		}, mc.cores=8)
		assign(paste0(k,"_data"),dat)
		cat(k, " ")
	}

	## run models
	cat("\nRunning models: \n")
	for(k in ped_names){
		model1 <- model_func(m1_func,get(paste0(k,"_peds")),get(paste0(k,"_data")),mc.cores=8)
		assign(paste0("model1_",k),model1)
		model2 <- model_func(m2_func,get(paste0(k,"_peds")),get(paste0(k,"_data")),mc.cores=8)
		assign(paste0("model2_",k),model2)
		cat(k, " ")
	}

	save(list=(c(ped_str,paste0("model1_",ped_names),paste0("model2_",ped_names))),file=paste0(data_wd,"mge_sims.Rdata"))

}

load(paste0(data_wd,"mge_sims.Rdata"))

# for(k in ped_names){
# 	mod2 <- do.call(rbind,lapply(get(paste0("model2_",k)), function(x) {
# 		do.call(rbind,lapply(1:nrow(scenarios), function(i) data.frame(r=k,scenario=i,comp=c("A","Me"),estimate=x[i,c("A","Me")])))
# 	}))
# 	assign(paste0("mod2_",k),mod2)
# }

total_Va<-do.call(rbind,lapply(ped_names,function(k) {
	mod2 <- do.call(rbind,lapply(get(paste0("model2_",k)), function(x) {
		do.call(rbind,lapply(1:nrow(scenarios), function(i) data.frame(r=k,scenario=i,max=x[i,"A"] + 0.5*x[i,"Me"], min=x[i,"A"] )
		))
	}))}
	))
total_Va$coverage <- total_Va$max>0.2 & total_Va$min<0.2
aggregate(coverage~ r+scenario, total_Va, mean)

mod2<-do.call(rbind,lapply(ped_names,function(k) {
	mod2 <- do.call(rbind,lapply(get(paste0("model2_",k)), function(x) {
		do.call(rbind,lapply(1:nrow(scenarios), function(i) data.frame(r=k,scenario=i,Va_est = x[i,"A"],Vm_est = x[i,"Me"],Va_sim=scenarios[i,"Va"],Vm_sim =sum(scenarios[i,c("Vmg","Vme")]),Vmg_sim=scenarios[i,"Vmg"])))
	}))
	# assign(paste0("mod2_",k),mod2)
}))
#,sum(x[i,c("Mg","Me")])
head(mod2,20)
mod2$Va_bias <- mod2$Va_est - mod2$Va_sim



 # mod2<-rbind(mod2_fs,mod2_hs,mod2_fhs,mod2_fhs_highF,mod2_fs_highF,mod2_fs_I,mod2_fs_IF,mod2_fhs_IF)

cols <- viridis::inferno(6)[c(2,3,5,4,6)]


scenarios2 <- scenarios[,c(1,4,2,3)]
scenarios2[,"r_amg"] <- scenarios2[,"r_amg"]*sqrt(scenarios2[,"Va"])*sqrt(scenarios2[,"Vmg"])
scenarios2 <- cbind(scenarios2,1 - rowSums(scenarios2) - scenarios2[,"r_amg"])
colnames(scenarios2) <- c("A","Me","Mg","cov_AMg","E")

s3 <- cbind(scenarios2, 0)[,c(4,6,1,2,3,5)]
s4<-t(apply(s3,1,function(x) c(x[1], cumsum(x[2:6]))))

ped_n <- length(ped_names)

library(beeswarm)
par(mfrow=c(2,2),mar=c(4,4,1,1))
# layout(matrix(c(1:3,0,4:11),nrow=3, byrow=TRUE))
for(i in 1:4){#nrow(scenarios)
	beeswarm(estimate~ r, mod2, subset=scenario==i&comp=="A",pch=19, cex=0.2, col=scales::alpha(1,0.3),method = "compactswarm",corral="wrap",  ylim=c(0,0.8),xlim=c(-1,12))#, xlim=c(-1,8)
	for(j in 1:5){
		polygon(x=c(-1,0,0,-1),y=c(s4[i,j],s4[i,j],s4[i,j+1],s4[i,j+1]), col=cols[j])	
	}
	means<-aggregate(estimate~ r, mod2,mean, subset=scenario==i&comp=="A")$estimate

	arrows((1:(ped_n*2))-0.25,means,(1:(ped_n*2))+0.25, code=0, col="blue")
	arrows((1:(ped_n*2))-0.25,rep(c(scenarios2[i,1],sum(scenarios2[i,2:3])),each=ped_n),(1:(ped_n*2))+0.25, code=0, col="red")

	# text(1,0.4,paste(paste(colnames(scenarios),"=",scenarios[i,]),collapse=", "),pos=4)
	
}

par(mfrow=c(2,2),mar=c(4,4,1,1))
for(i in 1:4){#nrow(scenarios)
	beeswarm(Va_bias~ r, mod2, subset=scenario==i,pch=19, cex=0.2, col=scales::alpha(1,0.3),method = "compactswarm",corral="wrap",  xlim=c(1,12), ylim=c(-0.1,0.35))

	means<-aggregate(Va_bias~ r, mod2,mean, subset=scenario==i)$Va_bias

	arrows((1:(ped_n*2))-0.25,means,(1:(ped_n*2))+0.25, code=0, col="blue")

	# text(1,0.4,paste(paste(colnames(scenarios),"=",scenarios[i,]),collapse=", "),pos=4)	
}

va<-aggregate(cbind(Va_bias,Vmg_sim)~ scenario+r, mod2,mean)
va$bias_prop <- va$Va_bias/va$Vmg_sim

beeswarm(bias_prop~ r, va,pch=19, cex=1, col= va$scenario,method = "compactswarm",corral="wrap")
plot(bias_prop~ as.numeric(as.factor(r)), va,pch=19, cex=1,col= va$scenario)
plot(Va_bias~ as.numeric(as.factor(r)), va,pch=19, cex=1,col= va$scenario)

#scales::alpha(1,0.3)

## correlate means with ratio of maternal link 



ped_str<-mclapply(ped_names,function(k){
	ped <- get(paste0(k,"_peds"))[[1]]
	ped_stat(ped, phenotyped=ped[!is.na(ped[,2]),1])
},mc.cores=8)
ped_str<-do.call(cbind,ped_str)

dam





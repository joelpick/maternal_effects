
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

##------------------------
## Pedigree parameters
##------------------------

generations=5
n_offspring=600

## Range from bonnet et al 
m_fecundity = 6 # mean is 6.9
h_fecundity = 12 #highest is 13.3
l_fecundity = 3 #lowest is 1.5
## these numbers are easier to maintain total number of phenotypes individuals

## first number is female, then male immigration rate
no_immigration = c(0,0)
f_immigration=c(0.4,0.1)
m_immigration=c(0.1,0.4)
b_immigration=c(0.25,0.25) # same overall immigration but no sex bias

fhs = 0.75
fs = 1
hs = 0

fecundity <- c(m_fecundity,h_fecundity,l_fecundity)
immigration <- rbind(no_immigration,f_immigration,m_immigration,b_immigration)
p_sire <- c(fhs,hs,fs)

f_names <- c("mF","hF","lF")
i_names <- c("nI","fI","mI","bI")
ms_names <- c("fhs","hs","fs")

## make all combos 
combos<-expand.grid(p_sire=1:3,immigration=1:4,fecundity=1:3)

peds_param <- cbind(
	n_females=n_offspring / fecundity[combos[,"fecundity"]],
	fecundity=fecundity[combos[,"fecundity"]],
	p_sire = p_sire[combos[,"p_sire"]],
	immigration_f = immigration[combos[,"immigration"],1],
	immigration_m = immigration[combos[,"immigration"],2])

## work out juvenile survival, so that populations are stationary
peds_param <- cbind(peds_param, 
	juv_surv_f=2*(1 - peds_param[,"immigration_f"])/peds_param[,"fecundity"],
	juv_surv_m=2*(1 - peds_param[,"immigration_m"])/peds_param[,"fecundity"]
)

## generate pedigree names
ped_names <- rownames(peds_param) <- paste(
	ms_names[combos[,"p_sire"]],
	f_names[combos[,"fecundity"]],
	i_names[combos[,"immigration"]], sep="_")


##

scenarios <- rbind(	
	# A) Maternal genetic only
	a=c(Va=0, Vmg=0.25, Vme=0, r_amg=0),
	# B) Direct genetic and maternal environment
	b=c(Va=0, Vmg=0.25, Vme=0.25, r_amg=0),
	# C) high Maternal genetic only
	c=c(Va=0, Vmg=0.5, Vme=0, r_amg=0),
	# D) high Maternal environment only
	d=c(Va=0, Vmg=0, Vme=0.5, r_amg=0),
	# E) Direct and maternal genetic, no covariance
	e=c(Va=0.25, Vmg=0.25, Vme=0, r_amg=0),
	# F) Direct and maternal environment, no covariance
	f=c(Va=0.25, Vmg=0, Vme=0.25, r_amg=0),
	# G) Direct and maternal genetic, no covariance and maternal environment
	g=c(Va=0.25, Vmg=0.25, Vme=0.25, r_amg=0),
	# H) Direct and maternal genetic, + covariance and maternal environment
	h=c(Va=0.25, Vmg=0.25, Vme=0, r_amg=0.5),
	# I) Direct and maternal genetic, - covariance and maternal environment
	i=c(Va=0.25, Vmg=0.25, Vme=0, r_amg=-0.5)
		
)

if(run){

	set.seed(20230126)

	ped_str <- vector("list",length=nrow(peds_param))
	names(ped_str) <- ped_names
	ped_str_mat <- vector("list",length=nrow(peds_param))
	names(ped_str_mat) <- ped_names

# k=ped_names[6]
	## make pedigrees
	
	for(k in ped_names){
		cat(k, "\n")
		cat("Simulating Pedigrees\n")
		peds <- mclapply(1:n_sims,	function(i){
			simulate_pedigree(
				years = generations,
				n_females = peds_param[k,"n_females"],
				fecundity = peds_param[k,"fecundity"],
				p_sire = peds_param[k,"p_sire"],
				p_polyandry=1,
				juv_surv = c(peds_param[k,"juv_surv_f"],peds_param[k,"juv_surv_m"]),
				adult_surv = 0,					# discrete generations
				immigration = c(peds_param[k,"immigration_f"],peds_param[k,"immigration_m"]), 				# closed population
				constant_pop = TRUE     # constant population size
				)$pedigree
		}, mc.cores=8)
		# assign(paste0(k ,"_peds"),peds)	

		cat("Generating Pedigree Metrics\n")
		ped_str[[k]]<- do.call(rbind,mclapply(peds,ped_stat, mc.cores=8))
		ped_str_mat[[k]]<- do.call(rbind,mclapply(peds,ped_stat2, mc.cores=8))
		
	

		cat("Simulating Data\n")
		## simulate data
	
		dat<-mclapply(peds, function(i){
			x<-vector("list", nrow(scenarios))
			for(j in 1:nrow(scenarios)){
				x[[j]]<- mge_sim(i[,1:3], param=scenarios[j,])
			}
			x
		}, mc.cores=8)
		# assign(paste0(k,"_data"),dat)

	## run models
		cat("Running models: \n")
	
		cat("Model 1: ")
		model1 <- model_func(m1a_func,peds,dat,mc.cores=8)
		assign(paste0("model1_",k),model1)
		cat("\nModel 2: ")
		model2 <- model_func(m2_func,peds,dat,mc.cores=8)
		assign(paste0("model2_",k),model2)
		cat("\n")
		rm(peds,dat)
	}

	save(list=(c("ped_names","scenarios","ped_str","ped_str_mat",paste0("model1_",ped_names),paste0("model2_",ped_names))),file=paste0(data_wd,"mge_sims3.Rdata"))
}

load(paste0(data_wd,"mge_sims2.Rdata"))



ped_sum<-sapply(ped_str,colMeans)
ped_sum_mat<-sapply(ped_str_mat,colMeans)

ped_sum_se<-sapply(ped_str,function(x) apply(x,2,se))
ped_sum_mat_se<-sapply(ped_str_mat,function(x) apply(x,2,se))


pedC_sum<-sapply(ped_cors,colMeans)
rownames(pedC_sum) <- c("A-Mg","A-Me","Mg-Me")

ped_sum2 <- rbind(mat_sibs=colSums(ped_sum[c("FS","MHS"),]), mat_links=colSums(ped_sum[c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS"),]), other=colSums(ped_sum[!rownames(ped_sum)%in%c( "individuals" ,"links" ,"FS","MHS","dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS"),]) )

plot(ped_sum_mat["mat_links",],colSums(ped_sum2[c("mat_sibs","mat_links"),]))

plot(ped_sum_mat["total_links",],ped_sum2[c("other"),])

plot((ped_sum_mat["mat_links",] - ped_sum_mat["mat_sib",])/ped_sum_mat["total_links",],ped_sum2[2,]/colSums(ped_sum2)); abline(0,1)

ps <- t(t(ped_sum2)/colSums(ped_sum2))

plot(ps[2,],pedC_sum[1,])
plot(ps[2,],pedC_sum[3,])

# ps <- t(t(ped_sum2[-3,])/colSums(ped_sum2[-3,]))
barplot(ps)
par(mfrow=c(3,1), mar=c(4,4,1,1), cex.lab=1.4,mgp=c(2,0.5,0))
barplot(ps[,c(5,1,4)], names=c("Low", "Mid", "High"), xlab="Fecundity", col=c("lightblue","orange","white"))
barplot(ps[,c(1,8,6,7)], names=c("None", "Both sexes", "Female", "Male"), xlab="Dispersal", ylab="Prop. relationships",col=c("lightblue","orange","white"))
barplot(ps[,c(3,1,2)], names=c("Half-sibs", "Mixed", "Full- sibs"), xlab="Mating system",col=c("lightblue","orange","white"))

ps

mat_ratio<-ped_sum2[2,]/ped_sum2[1,]
mat_ratio2<-ped_sum2[2,]/colSums(ped_sum2)
mat_ratio3<-ped_sum2[1,]/colSums(ped_sum2)
mat_ratio4<-ped_sum2[1,]/ped_sum2[3,]
mat_ratio5<-ped_sum2[2,]/ped_sum2[3,]
mat_ratio6<-ped_sum2[3,]/colSums(ped_sum2)

mat_ratio2<-ped_sum2[2,]/colSums(ped_sum2)

matM_ratio <- (ped_sum_mat["mat_links",] - ped_sum_mat["mat_sib",])/ped_sum_mat["total_links",]

matM_ratio2 <- ped_sum_mat["maternal",]/ped_sum_mat["total",]
matM_ratio3 <- ped_sum2["mat_sibs",]/ped_sum_mat["maternal",]
plot(matM_ratio,mat_ratio2)
plot(matM_ratio2,mat_ratio3)

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
# total_Va$coverage <- total_Va$max>0.2 & total_Va$min<0.2
# aggregate(coverage~ r+scenario, total_Va, mean)


mod1<-do.call(rbind,lapply(ped_names,function(k) {
	mod1 <- do.call(rbind,lapply(get(paste0("model1_",k)), function(x) {
		do.call(rbind,lapply(1:nrow(scenarios), function(i) data.frame(r=k,scenario=i,Va_est = x[i,"A"],Va_sim=scenarios[i,"Va"],Vm_sim =sum(scenarios[i,c("Vmg","Vme")]),Vmg_sim=scenarios[i,"Vmg"])))
	}))
	# assign(paste0("mod2_",k),mod2)
}))
mod1$Va_bias <- mod1$Va_est - mod1$Va_sim
mod1$ln_Va_bias <- log(mod1$Va_bias)
# va1<-aggregate(cbind(Va_bias,Vmg_sim,Vm_sim)~ scenario+r, mod1,median)
va1<-aggregate(cbind(Va_bias,Vmg_sim,Vm_sim,ln_Va_bias)~ scenario+r, mod1,mean)
# which(va1$r)

r_order<- sapply(va1$r, function(x) which(rownames(peds_param)==x))

plot(Va_bias~ log(mat_ratio[r_order]), va1,pch=19, cex=1,col= va1$scenario)
# plot(Va_bias~ log(mat_ratio2[r_order]), va1,pch=19, cex=1,col= va1$scenario)

plot(ln_Va_bias~ log(mat_ratio[r_order]), va1,pch=19, cex=1,col= va1$scenario)

plot(Va_bias~ r_order, va1,pch=19, cex=1,col= va1$scenario, xaxt="n")
axis(1,1:8,rownames(peds_param))



cov <- read.csv(paste0(wd,"Data/Raw/covariances.csv"))
sum((ped_sum[cov$relationship,1]*cov[,"Va"]*scenarios[1,"Va"] + ped_sum[cov$relationship,1]*cov[,"Vmg"]*scenarios[1,"Vmg"] + ped_sum[cov$relationship,1]*cov[,"Vme"]*scenarios[1,"Vme"])/sum(ped_sum[cov$relationship,1])*1/cov[,"Va"])


mod2<-do.call(rbind,lapply(ped_names,function(k) {
	mod2 <- do.call(rbind,lapply(get(paste0("model2_",k)), function(x) {
		do.call(rbind,lapply(1:nrow(scenarios), function(i) data.frame(r=k,scenario=i,Va_est = x[i,"A"],Vm_est = x[i,"Me"],Va_sim=scenarios[i,"Va"],Vm_sim =sum(scenarios[i,c("Vmg","Vme")]),Vmg_sim=scenarios[i,"Vmg"])))
	}))
	# assign(paste0("mod2_",k),mod2)
}))
r_order<- sapply(va2$r, function(x) which(rownames(peds_param)==x))


#,sum(x[i,c("Mg","Me")])
head(mod2,20)
mod2$Va_bias <- mod2$Va_est - mod2$Va_sim
# mod2$ln_Va_bias <- log(mod2$Va_bias)
va2<-aggregate(cbind(Va_bias,Vmg_sim,Vm_sim)~ scenario+r, mod2,mean)
va2_se<-aggregate(Va_bias~ scenario+r, mod2,se)


# which(va1$r)

va2$mat_ratio <- mat_ratio2[r_order]
va2$matM_ratio<- matM_ratio[r_order]
va2[,c("ms","fec","imm")] <- do.call(rbind,strsplit(va2$r,"_"))
# va2$matM_ratio2<- matM_ratio2[r_order]
par(mfrow=c(1,2), mar=c(5,5,1,1), cex.lab=1.75, cex.axis=1.25 )
plot(Va_bias~ mat_ratio, va2, subset=scenario==1, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~h^2), col=(1:4)[as.factor(va2$imm)], pch=c(15:17)[as.factor(va2$ms)])
arrows(va2$mat_ratio,va2$Va_bias+va2_se$Va_bias,va2$mat_ratio,va2$Va_bias-va2_se$Va_bias,code=3,angle=90,length=0.1)
plot(Va_bias~ matM_ratio, va2, subset=scenario==1, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~h^2), col=(1:4)[as.factor(va2$imm)], pch=c(15:17)[as.factor(va2$ms)])
arrows(va2$matM_ratio,va2$Va_bias+va2_se$Va_bias,va2$matM_ratio,va2$Va_bias-va2_se$Va_bias,code=3,angle=90,length=0.1)
# plot(Va_bias~ matM_ratio2, va2, subset=scenario==1,  pch=19, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~h^2))


mod2_1 <- subset(mod2,scenario==1)

mat_propM<- do.call(rbind,ped_str_mat)
mod2_1$mat_propM <- mat_propM[,1]/mat_propM[,2]

mat_prop<- do.call(rbind,ped_str)
mod2_1$mat_sibs <- rowSums(mat_prop[,c("FS","MHS")])

mat_prop2 <- cbind(mat_sibs=rowSums(mat_prop[,c("FS","MHS")]), mat_links=rowSums(mat_prop[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")]), other=rowSums(mat_prop[,!colnames(mat_prop)%in%c( "individuals" ,"links" ,"FS","MHS","dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")]) )


plot(Va_bias~ mat_propM,mod2_1, pch=19, col=c(1:7)[as.factor(mod2_1$r)])
cor(mod2_1$Va_bias, mod2_1$mat_propM)

sapply(ped_names,function(i){
	x<-subset(mod2_1, r==i)
	cor(x$Va_bias, x$mat_propM)
})


plot(Va_bias~ mat_ratio, va2,pch=19, cex=1,col= va2$scenario, xlab="Proportion non-sibling maternal links", ylab="Bias in Va")
legend("topleft",apply(scenarios[,c(2,1,4)],1,function(x) paste(colnames(scenarios[,c(2,1,4)]),x,collapse=", ", sep="=") ), pch=19, col=1:4)

# plot(Va_bias~ pedC_sum[3,][r_order], va2,pch=19, cex=1,col= va2$scenario, xlab="Proportion non-sibling maternal links", ylab="Bias in Va")


# mod2<-sapply(ped_names,function(k) {
# 	mod2 <- do.call(rbind,lapply(get(paste0("model2_",k)), function(x) {
# 		do.call(rbind,lapply(1:nrow(scenarios), function(i) 

# 			c(r=k,scenario=i,Va_est = x[i,"A"],Vm_est = x[i,"Me"],Va_sim=scenarios[i,"Va"],Vm_sim =sum(scenarios[i,c("Vmg","Vme")]),Vmg_sim=scenarios[i,"Vmg"])))
# 	}))
# 	# assign(paste0("mod2_",k),mod2)
# })


# plot(Va_bias~ log(mat_ratio2[r_order]), va1,pch=19, cex=1,col= va1$scenario)

plot(ln_Va_bias~ log(mat_ratio[r_order]), va1,pch=19, cex=1,col= va1$scenario)

plot(Va_bias~ r_order, va1,pch=19, cex=1,col= va1$scenario, xaxt="n")
axis(1,1:8,rownames(peds_param))





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
	beeswarm(Va_est~ r, mod2, subset=scenario==i,pch=19, cex=0.2, col=scales::alpha(1,0.3),method = "compactswarm",corral="wrap",  ylim=c(0,0.8),xlim=c(-1,12))#, xlim=c(-1,8)
	for(j in 1:5){
		polygon(x=c(-1,0,0,-1),y=c(s4[i,j],s4[i,j],s4[i,j+1],s4[i,j+1]), col=cols[j])	
	}
	means<-aggregate(Va_est~ r, mod2,mean, subset=scenario==i)$Va_est

	arrows((1:(ped_n*2))-0.25,means,(1:(ped_n*2))+0.25, code=0, col="blue")
	arrows((1:(ped_n*2))-0.25,rep(c(scenarios2[i,1],sum(scenarios2[i,2:3])),each=ped_n),(1:(ped_n*2))+0.25, code=0, col="red")

	# text(1,0.4,paste(paste(colnames(scenarios),"=",scenarios[i,]),collapse=", "),pos=4)
	
}

par(mfrow=c(2,2),mar=c(4,4,1,1))
for(i in 1:4){#nrow(scenarios)
	beeswarm(Va_bias~ r, mod2, subset=scenario==i,pch=19, cex=0.2, col=scales::alpha(1,0.3),method = "compactswarm",corral="wrap",  xlim=c(1,8), ylim=c(-0.1,0.35))

	means<-aggregate(Va_bias~ r, mod2,mean, subset=scenario==i)$Va_bias

	arrows((1:(ped_n))-0.25,means,(1:(ped_n))+0.25, code=0, col="blue")

	# text(1,0.4,paste(paste(colnames(scenarios),"=",scenarios[i,]),collapse=", "),pos=4)	
}


par(mfrow=c(2,2),mar=c(4,4,1,1))
for(i in 1:4){#nrow(scenarios)
	beeswarm(Vm_est~ r, mod2, subset=scenario==i,pch=19, cex=0.2, col=scales::alpha(1,0.3),method = "compactswarm",corral="wrap",  xlim=c(1,8), ylim=c(0,0.7))

	means<-aggregate(Vm_sim~ r, mod2,mean, subset=scenario==i)$Vm_sim
	means<-aggregate(Vm_sim~ r, mod2,mean, subset=scenario==i)$Vm_sim

	arrows((1:(ped_n))-0.25,means,(1:(ped_n))+0.25, code=0, col="blue")

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





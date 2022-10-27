rm(list=ls())

library(viridis)
library(beeswarm)

wd <- "/Users/joelpick/github/maternal_effects/"
data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))
# load(paste0(data_wd,"gaussian_data.Rdata"))
load(paste0(data_wd,"gaussian_sims.Rdata"))

comp_names <- c("A","Me","Mg","cov_AMg","E")
model_names <- c(
	"A",
	"A + Me",
	"A + Mg",
  "A + Mg + Me",
  "A + Mg + cov + Me",
  "Mg",
  "Mg + Me",
  "Mg average"
	)

scenarios2 <- scenarios[,c(1,4,2,3)]
scenarios2[,"r_amg"] <- scenarios2[,"r_amg"]*sqrt(scenarios2[,"Va"])*sqrt(scenarios2[,"Vmg"])
scenarios2 <- cbind(scenarios2,1 - rowSums(scenarios2) - scenarios2[,"r_amg"])
colnames(scenarios2) <- c("A","Me","Mg","cov_AMg","E")

s3 <- cbind(scenarios2, 0)[,c(4,6,1,2,3,5)]
s4<-t(apply(s3,1,function(x) c(x[1], cumsum(x[2:6]))))


# list2array(model2_fs)[1,"A",]

s <- letters[c(1,2,4,3,5,6,9,7,8,10,11)]

mod2_fs<-do.call(rbind,lapply(model2_fs, function(x) {
	do.call(rbind,lapply(s, function(i) data.frame(r="fs",scenario=i,comp=c("A","Me"),estimate=x[i,c("A","Me")])))
}))
mod2_hs<-do.call(rbind,lapply(model2_hs, function(x) {
	do.call(rbind,lapply(s, function(i) data.frame(r="hs",scenario=i,comp=c("A","Me"),estimate=x[i,c("A","Me")])))
}))
mod2_fhs<-do.call(rbind,lapply(model2_fhs, function(x) {
	do.call(rbind,lapply(s, function(i) data.frame(r="fhs",scenario=i,comp=c("A","Me"),estimate=x[i,c("A","Me")])))
}))
mod2_fhs10<-do.call(rbind,lapply(model2_fhs10, function(x) {
	do.call(rbind,lapply(s, function(i) data.frame(r="fhs10",scenario=i,comp=c("A","Me"),estimate=x[i,c("A","Me")])))
}))

mod2<-rbind(mod2_fs,mod2_hs,mod2_fhs,mod2_fhs10)
nrow(mod2)

cols <- inferno(6)[c(2,3,5,4,6)]

par(mfrow=c(3,4))
layout(matrix(c(1:3,0,4:11),nrow=3, byrow=TRUE))
for(i in s){
	beeswarm(estimate~ r+comp, mod2, subset=scenario==i,pch=19, cex=0.2, col=scales::alpha(1,0.3),method = "compactswarm",corral="wrap",  ylim=c(0,1), xlim=c(-1,8))
	for(j in 1:5){
		polygon(x=c(-1,0,0,-1),y=c(s4[i,j],s4[i,j],s4[i,j+1],s4[i,j+1]), col=cols[j])	
	}
	arrows((1:4)-0.25,rep(c(scenarios2[i,1],sum(scenarios2[i,2:3])),each=2),(1:4)+0.25, code=0, col="red")
	
}



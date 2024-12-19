
rm(list=ls())

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Raw/")

############
## young and postma 2023
############

YP <- read.csv(paste0(data_wd,"young_postma_2023_data.csv"))

YP2 <- subset(YP,YP$method %in% c("animal model (MCMC)","animal model (REML)"))
head(YP2)

# pedigree sizes
quantile(YP2$n,c(0,0.1,0.25,0.5,0.75,0.90,0.95,1),na.rm=TRUE)
hist((YP2$n), xlim=c(0,4000), breaks=seq(0,40000,100))

## mean SE
mean(YP2$SE.h2,na.rm=TRUE)


############
## postma 2014
############

postma<-read.csv(paste0(data_wd,"postma2014_SM.csv"))

postma2 <- subset(postma,postma$method %in% c("animal model"))
head(postma2)

## mean SE
mean(postma2$SE.of.h2,na.rm=TRUE)


############
## moore et al.
############

moore<-read.csv(paste0(data_wd,"moore_2019_datawithSE.csv"))
head(moore)

## mean SE
mean(moore$h2.samp.err)


############
## Fecundities from Bonnet et al. 2022
############

pops<- list.files(paste0(data_wd,"Bonnet/"))

fecundies<-sapply(pops, function(i){
	dat <- read.csv( paste0(data_wd,"Bonnet/",i))
	mean(table(dat$dam))
})
range(fecundies)
median(fecundies)

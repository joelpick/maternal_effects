rm(list=ls())
wd <- "/Users/joelpick/github/maternal_effects/"

source(paste0(wd,"R/extract_cousins.R"))
source("/Users/joelpick/github/squidPed/R/simulate_pedigree.R")

generations=10
n_females=100
fecundity=4
juv_surv = 2/fecundity # 1/fecundity
immigration = 0
adult_surv = 0#0.5

ped <-	simulate_pedigree(
		years = generations,
		n_females = n_females,
		fecundity = fecundity,
		p_sire = 0.5, 				# mating system (0-1, 1= one male per female, 0=complete random mating)
		juv_surv = juv_surv, # insures no population growth
		adult_surv = adult_surv,					# discrete generations
		immigration = immigration, 				# closed population
		constant_pop = TRUE     # constant population size
		)$pedigree
 ped_stat(ped)

 barplot(mat_links(ped))
 mat_links(ped)/sum(mat_links(ped))
pedSum(ped[,1:3])

ped[sample(1:nrow(ped),200),3]<-NA
cousins(ped)- cousins2(ped)
au(ped)- au2(ped)

out<-parallel::mclapply(1:100,function(x){
	ped <-	simulate_pedigree(
		years = generations,
		n_females = n_females,
		fecundity = fecundity,
		p_sire = 1, 				# mating system (0-1, 1= one male per female, 0=complete random mating)
		juv_surv = juv_surv, # insures no population growth
		adult_surv = adult_surv,					# discrete generations
		immigration = immigration, 				# closed population
		constant_pop = TRUE     # constant population size
		)$pedigree
	ped_stat(ped)	
}, mc.cores=8)

 
rowMeans(do.call(cbind,out))



pedSum(ped[,1:3])





#expected relationships

# number of mothers
# n_females + generations * fecundity * juv_surv_f/2 + immigration_f*n_females
# starting females + new mothers produces each year * years
n_mothers <- n_females + (generations-1) * (n_females * fecundity * juv_surv/2 + immigration*n_females)

length(unique(na.omit(ped[,2])))



# average number of offspring 
# fecundity * 1/(1-adult_surv)
LRS <- fecundity * 1/(1-adult_surv)


#dam-offspring
# number of mothers * n_offspring
## but this assumes their females will all live out their life
## maybe something with geometric probs for each year? 
# dgeom(1:10,0.5)*n_females
# n_mothers * LRS
generations * fecundity * n_females # if constant pop - just how many offspring are produced over the observation period 


## sibs
# number of mothers  * (n_offspring*(n_offspring-1)/2)
n_mothers * LRS*(LRS-1)/2

## full sibs
prop_FS<-0.5
n_mothers * (LRS*prop_FS) *((LRS*prop_FS)-1)/2



## half sibs
## paternal half sibs == number of sibs??

n_recruits <- fecundity * juv_surv
## grandoffspring
# number of mothers * juv_surv * n_offspring^2

## number of maternal granddaughters produced per year =  recruitment * fecundity - last generation wont have any grandoffspring
(generations-1) * juv_surv * fecundity^2 * n_females / 2 ## for each PGM and MGM - depending on sex ratio

# and then immigration
# + immigration*n_females*(generations-2) * juv_surv * fecundity^2 / 2 ?
# MGO of immigrant females - immigrants only appear in generation 2 
# maybe this isnt needed if immigration counteracts lower JS to produce same pop size each year

# aunts/uncles

n_recruits * fecundity * (fecundity-1) * (generations-1) * n_females /2
## divide by two for realted through mother or fathers siblings

## yes for FS designs
## also AU through phs?


# cousins

(utils::combn(rep(fecundity,n_recruits), m =2)[1, ] * utils::combn(rep(fecundity,n_recruits), m =2)[2, ] )* (generations-1) * n_females /2

biases(ped,cov)


#################################################
rbind_notAnnoying <- function(..., names=NULL){
  x <- list(...)
  y <- lapply(x, function(y){
    names(y) <- if(is.null(names)) names(x[[1]]) else names
    return(y)  
  } )
  do.call(rbind,y)
}

## function to make products of all combinations and sum them 
combo_prod <- function(x) if(length(x)>1){sum( combn(x, m =2)[1, ] * combn(x, m =2)[2, ])}else{0}


combo_au <- function(x) if(length(x)>1){sum(x * (length(x)-1))}else{0}


n_cousin <- function(p,gp,ped){
  f1 <- formula(paste("animal~",p,"+",gp))
  f2 <- formula(paste("animal~",gp))
  d1 <- aggregate(f1, ped,length)
	sum(aggregate(f2,d1,combo_prod)$animal)
}

n_au <- function(link,ped){
	
	p <- if(substr(link,1,1)=="M"){"dam"}else{"sire"}
	gp <- if(substr(link,3,3)=="M"){"dam"}else if(substr(link,3,3)=="F"){"sire"} else{"pair"}
  f1 <- formula(paste("animal~",p,"+",link))
  f2 <- formula(paste("animal~",gp))
  d1 <- aggregate(f1, ped,length)
  d2 <- aggregate(f2, ped,length)
	sum(d1$animal * (d2[match(d1[,link],d2[,gp]),"animal"]-1))
}
# cousins(ped)






n_au(p="dam", gp="MGP", link=, ped=ped)


ped[2089:2096,]

nrow(subset(ped,is.na(MGF) & !is.na(MGM)))
nrow(subset(ped,is.na(PGF) & !is.na(MGM)))
cbind(cousins(ped),cousins2(ped))

au2 <- function(ped){	

 	dummy_dam <- which(is.na(ped$dam) & !is.na(ped$sire))
	dummy_sire <- which(is.na(ped$sire) & !is.na(ped$dam))

	ped[dummy_sire,"sire"] <- paste0("ds_",seq_along(dummy_sire))
	ped[dummy_dam,"dam"] <- paste0("dd_",seq_along(dummy_dam))


	ped$MGM <- ped[match(ped[,2],ped[,1]),2]
	ped$MGF <- ped[match(ped[,2],ped[,1]),3]
	ped$PGM <- ped[match(ped[,3],ped[,1]),2]
	ped$PGF <- ped[match(ped[,3],ped[,1]),3]

	ped$MGP <- paste(ped$MGM,ped$MGF)
	ped$MGP<-ifelse(ped$MGP== "NA NA", NA, ped$MGP)
	ped$PGP <- paste(ped$PGM,ped$PGF)
	ped$PGP<-ifelse(ped$PGP== "NA NA", NA, ped$PGP)
# ped$MGP<-ifelse(grepl("NA",MGP), NA, ped$MGP)
# ped$PGP<-ifelse(grepl("NA",PGP), NA, ped$PGP)

	ped$pair<-paste(ped$dam,ped$sire)
	ped$pair<-ifelse(ped$pair== "NA NA", NA, ped$pair)


	GM <- rbind_notAnnoying(ped[,c("animal","dam","MGM","MGF","MGP")], ped[,c("animal","sire","PGM","PGF","PGP")])

	# related through a maternal grandmother (au_FS_F + au_MHS_F ?)
	au_MS_F <- n_au("MGM",ped)

	# related through a maternal grandfather (au_FS_F + au_PHS_F ?)
	au_PS_F <- n_au("MGF",ped)

	# related through a maternal grandparents (au_FS_F ?)
	au_FS_F <- n_au("MGP",ped)

	au_PHS_F <- au_PS_F - au_FS_F
	au_MHS_F <- au_MS_F - au_FS_F

	# related through a paternal grandmother (au_FS_M + au_MHS_M ?)
	au_MS_M <- n_au("PGM",ped)

	# related through a paternal grandfather (au_FS_M + au_PHS_F ?)
	au_PS_M <- n_au("PGF",ped)

	# related through a paternal grandparents (au_FS_M ?)
	d1 <- aggregate(animal~sire+PGP,ped,length)
	d2 <- aggregate(animal~pair,ped,length)
	sum(d1$animal * (d2[match(d1$PGP,d2$pair),"animal"]-1))
	au_FS_M <- n_au("PGP",ped)

	au_PHS_M <- au_PS_M - au_FS_M
	au_MHS_M <- au_MS_M - au_FS_M




	c(au_FS_F=au_FS_F,au_FS_M=au_FS_M,au_MHS_F=au_MHS_F,au_MHS_M=au_MHS_M,au_PHS_F=au_PHS_F,au_PHS_M=au_PHS_M)

# c("au_FS_F","au_FS_M","au_MHS_F","au_MHS_M","au_PHS_F","au_PHS_M")

}
au2(ped)
system.time(au(ped))
system.time(au2(ped))


# system.time({
# sum(apply(table(ped$dam,ped$MGF),2, function(x){
#  y <- x[x>0]
# if(length(y)>1){sum( combn(y, m =2)[1, ] * combn(y, m =2)[2, ])}else{0}
# })
# )})

# system.time(n_cousin(p="dam", gp="MGF", ped=ped))

    2028     1679     2772     2695     2688     2373 
    1996     1651     2804     2723     2720     2401 

dim(table(ped$dam,ped$MGF))

system.time({	
	mat <- aggregate(animal~dam+MGM+MGF+MGP, ped,length, na.action=na.pass)
	cousin_FS_FF <- sum(aggregate(animal~MGP,mat,combo_prod)$animal)
	cousin_MS_FF <- sum(aggregate(animal~MGM,mat,combo_prod)$animal)
	cousin_PS_FF <- sum(aggregate(animal~MGF,mat,combo_prod)$animal)
  cousin_MHS_FF <- cousin_MS_FF - cousin_FS_FF
  cousin_PHS_FF <- cousin_PS_FF - cousin_FS_FF

  pat <- aggregate(animal~sire+PGM+PGF+PGP, ped,length)
	cousin_FS_MM <- sum(aggregate(animal~PGP,pat,combo_prod)$animal)
	cousin_MS_MM <- sum(aggregate(animal~PGM,pat,combo_prod)$animal)
	cousin_PS_MM <- sum(aggregate(animal~PGF,pat,combo_prod)$animal)
  cousin_MHS_MM <- cousin_MS_MM - cousin_FS_MM
  cousin_PHS_MM <- cousin_PS_MM - cousin_FS_MM

  mp <- aggregate(animal~dam+MGM+MGF+MGP, GM,length)
	cousin_FS <- sum(aggregate(animal~MGP,mp,combo_prod)$animal)
	cousin_MS <- sum(aggregate(animal~MGM,mp,combo_prod)$animal)
	cousin_PS <- sum(aggregate(animal~MGF,mp,combo_prod)$animal)

	cousin_FS_FM <- cousin_FS-cousin_FS_FF-cousin_FS_MM
	cousin_MHS_FM <- cousin_MS-cousin_FS-cousin_MHS_FF-cousin_MHS_MM
	cousin_PHS_FM <- cousin_PS-cousin_FS-cousin_PHS_FF-cousin_PHS_MM

})
c(cousin_FS_FF=cousin_FS_FF,cousin_FS_FM=cousin_FS_FM,cousin_FS_MM=cousin_FS_MM,cousin_MHS_FF=cousin_MHS_FF,cousin_MHS_FM=cousin_MHS_FM,cousin_MHS_MM=cousin_MHS_MM,cousin_PHS_FF=cousin_PHS_FF,cousin_PHS_FM=cousin_PHS_FM,cousin_PHS_MM=cousin_PHS_MM)


system.time({	
	# mat <- aggregate(animal~dam+MGM+MGF+MGP, ped,length)
mat1 <- aggregate(animal~dam+MGM, ped,length)
mat2 <- aggregate(animal~dam+MGF, ped,length)
mat3 <- aggregate(animal~dam+MGP, ped,length)
mat<-merge(merge(mat1,mat2,all=TRUE),mat3,all=TRUE)


	# mat <- aggregate(animal~dam+MGM+MGF+MGP, ped,length, na.action=na.pass)

	cousin_FS_FF <- sum(aggregate(animal~MGP,mat,combo_prod)$animal)
	cousin_MS_FF <- sum(aggregate(animal~MGM,mat,combo_prod)$animal)
	cousin_PS_FF <- sum(aggregate(animal~MGF,mat,combo_prod)$animal)
  cousin_MHS_FF <- cousin_MS_FF - cousin_FS_FF
  cousin_PHS_FF <- cousin_PS_FF - cousin_FS_FF


	pat1 <- aggregate(animal~sire+PGM, ped,length)
	pat2 <- aggregate(animal~sire+PGF, ped,length)
	pat3 <- aggregate(animal~sire+PGP, ped,length)
	pat<-merge(merge(pat1,pat2,all=TRUE),pat3,all=TRUE)


	cousin_FS_MM <- sum(aggregate(animal~PGP,pat,combo_prod)$animal)
	cousin_MS_MM <- sum(aggregate(animal~PGM,pat,combo_prod)$animal)
	cousin_PS_MM <- sum(aggregate(animal~PGF,pat,combo_prod)$animal)
  cousin_MHS_MM <- cousin_MS_MM - cousin_FS_MM
  cousin_PHS_MM <- cousin_PS_MM - cousin_FS_MM

mp1 <- aggregate(animal~dam+MGM, GM,length)
mp2 <- aggregate(animal~dam+MGF, GM,length)
mp3 <- aggregate(animal~dam+MGP, GM,length)
mp<-merge(merge(mp1,mp2,all=TRUE),mp3,all=TRUE)

	cousin_FS <- sum(aggregate(animal~MGP,mp,combo_prod)$animal)
	cousin_MS <- sum(aggregate(animal~MGM,mp,combo_prod)$animal)
	cousin_PS <- sum(aggregate(animal~MGF,mp,combo_prod)$animal)

	cousin_FS_FM <- cousin_FS-cousin_FS_FF-cousin_FS_MM
	cousin_MHS_FM <- cousin_MS-cousin_FS-cousin_MHS_FF-cousin_MHS_MM
	cousin_PHS_FM <- cousin_PS-cousin_FS-cousin_PHS_FF-cousin_PHS_MM

})
c(cousin_FS_FF=cousin_FS_FF,cousin_FS_FM=cousin_FS_FM,cousin_FS_MM=cousin_FS_MM,cousin_MHS_FF=cousin_MHS_FF,cousin_MHS_FM=cousin_MHS_FM,cousin_MHS_MM=cousin_MHS_MM,cousin_PHS_FF=cousin_PHS_FF,cousin_PHS_FM=cousin_PHS_FM,cousin_PHS_MM=cousin_PHS_MM)




system.time({	
cousin_FS_FF <- n_cousin(p="dam", gp="MGP", ped=ped)
cousin_MHS_FF <- n_cousin(p="dam", gp="MGM", ped=ped) - cousin_FS_FF
cousin_PHS_FF <- n_cousin(p="dam", gp="MGF", ped=ped) - cousin_FS_FF


cousin_FS_MM <- n_cousin(p="sire", gp="PGP", ped=ped)
cousin_MHS_MM <- n_cousin(p="sire", gp="PGM", ped=ped) - cousin_FS_MM
cousin_PHS_MM <- n_cousin(p="sire", gp="PGF", ped=ped) - cousin_FS_MM


cousin_FS <- n_cousin(p="dam", gp="MGP", ped=GM)
cousin_FS_FM <- cousin_FS-cousin_FS_FF-cousin_FS_MM

cousin_MS <- n_cousin(p="dam", gp="MGM", ped=GM)
cousin_MHS_FM <- cousin_MS-cousin_FS-cousin_MHS_FF-cousin_MHS_MM

cousin_PS <- n_cousin(p="dam", gp="MGF", ped=GM)
cousin_PHS_FM <- cousin_PS-cousin_FS-cousin_PHS_FF-cousin_PHS_MM


})
c(cousin_FS_FF=cousin_FS_FF,cousin_FS_FM=cousin_FS_FM,cousin_FS_MM=cousin_FS_MM,cousin_MHS_FF=cousin_MHS_FF,cousin_MHS_FM=cousin_MHS_FM,cousin_MHS_MM=cousin_MHS_MM,cousin_PHS_FF=cousin_PHS_FF,cousin_PHS_FM=cousin_PHS_FM,cousin_PHS_MM=cousin_PHS_MM)







system.time(cousins(ped))
system.time(cousins2(ped))

d1<-aggregate(animal~dam + MGF, GM,length) 
sum(aggregate(animal~MGF,d1,combo_prod)$animal)




	
	MGO_FS<-table(ped$MGP)
	PGO_FS<-table(ped$PGP)
	sum(MGO_FS*(MGO_FS-1)/2)+		sum(PGO_FS*(PGO_FS-1)/2)

subset(ped,MGP=="3_305 3_325")

table(ped$dam)
table(ped$sire)


sum(table(ped$MGM)*(table(ped$MGM)-1)/2)
sum(table(ped$MGM)*(table(ped$MGM)-1)/2)


sum(MGO_FS*(MGO_FS-1)/2)







# x<-subset(d1,MGP=="0_13 0_154")$animal
# sum( combn(x, m =2)[1, ] * combn(x, m =2)[2, ])



cousins(ped)




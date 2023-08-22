
rm(list=ls())

extract<-FALSE
wd <- "/Users/joelpick/github/maternal_effects/"
main_dir <- paste0(wd,"Data/Raw/Bonnet/")
data_dir <- paste0(wd,"Data/Intermediate/")

covariances <- read.csv(paste0(wd,"Data/Raw/covariances.csv"))

source(paste0(wd,"R/extract_cousins.R"))
source(paste0(wd,"R/00_functions.R"))
source("/Users/joelpick/github/squidPed/R/simulate_pedigree.R")

pops<- list.files(main_dir)

# i=pops[11]
	## run export R_MAX_VSIZE=32000000000 in terminal before running this 
all<-do.call(rbind,parallel::mclapply(pops[c(1:3,9:12,15:19)], function(i){
	pop_dir <- paste0(main_dir,i)
	load(paste0(pop_dir,"/Ainv_",i))


# x<-Matrix::solve(ainv,sparse=TRUE)
# A<-Matrix::forceSymmetric(round(Matrix::solve(ainv,sparse=TRUE),10))
	# A<-as(A,"symmetricMatrix")


	A <- Matrix::forceSymmetric(round(Matrix::solve(ainv,sparse=FALSE),10),uplo="L")
	A<-as(A,"dsCMatrix")
	colnames(A) <- rownames(A) <- rownames(ainv)
	# rm(ainv)
	dat<-read.csv(paste0(pop_dir,"/data_",i,".csv"))
	A_dam<-A[as.character(dat$dam),as.character(dat$dam)]
	A_id<-A[as.character(dat$id),as.character(dat$id)]
A_dam[1:10,1:10]
A_id[1:10,1:10]

	sDam<-Matrix::summary(A_dam)
	names(sDam)[3] <- "Vmg"
	sID<-Matrix::summary(A_id)
	names(sID)[3] <- "Va"

nrow(sDam)
nrow(sID)

	both<-as.data.frame(dplyr::left_join(sID,sDam,by=c("i","j")))
	both[is.na(both[,"Vmg"]),"Vmg"] <- 0
	both[,c("Va","Vmg")]<-round(both[,c("Va","Vmg")],6)
	mat_sibs<- table(dat$dam)
	c(
	all= (nrow(dat)^2-nrow(dat))/2,
	mat_sibs = sum(mat_sibs * (mat_sibs-1)/2),
	mat_links = nrow(sDam)-nrow(dat),
	total_links = nrow(sID)-nrow(dat),
	linkSub=sum(both$Va>=0.0625),
	matlinkSub=sum(both$Va>=0.0625 & both$Vmg>=0.25),
	dam = sum(both$Va==0.5 & both$Vmg==0.5),
	FS = sum(both$Va==0.5 & both$Vmg==1),
	MHS =	sum(both$Va==0.25 & both$Vmg==1),
	MG=	sum(both$Va==0.25 & both$Vmg==0.25),
	au_D_FS	=sum(both$Va==0.25 & both$Vmg==0.5),
	au_D_MHS =sum(both$Va==0.125 & both$Vmg==0.25),
	cousin_D_FS	=sum(both$Va==0.125 & both$Vmg==0.5),
	cousin_D_HS	=sum(both$Va==0.0625 & both$Vmg==0.25)
	)
}, mc.cores=6))
all
rownames(all)<-pops[c(1:3,9:12,15:19)]
round(all/rowSums(all),3)
rowSums(all[,2:3])/rowSums(all[,-(2:3)])

rowSums(all[,-c(1:3,5:6)])/all[,3]


plot(all[,"matlinkSub"],all[,"mat_links"],)

all[,"mat_links"]/all[,"total_links"]


A_dam<-A[as.character(dat$dam),as.character(dat$dam)]
A_id<-A[as.character(dat$id),as.character(dat$id)]
A_id_dam<-A[as.character(dat$id),as.character(dat$dam)]

		# ainv <- methods::as(ainv, "dgCMatrix")

# cor(as.numeric(A_id),as.numeric(A_dam))

sum(round(A_dam,6)[lower.tri(round(A_dam,6))]>0)

lt_id<-round(A_id,6)[lower.tri(round(A_id,6))]
lt2_id<-lt_id[lt_id>0]
hist(lt2_id, breaks=1000)
abline(v=1/2^c(1:9), col="red",lwd=0.5)


lt_dam<-round(A_dam,6)[lower.tri(round(A_dam,6))]
lt2_dam<-lt_dam[lt_dam>0]
hist(lt2_dam, breaks=1000)
abline(v=1/2^c(1:8), col="red",lwd=0.5)

lt<-data.frame(Va=lt_id,Vmg=lt_dam)
# lt_combos<-aggregate(n~lt_id+lt_dam,lt, sum)
# nrow(lt_combos[lt_combos$n>100,])



# nrow(lt[lt$Va==0.5 & lt$Vmg==0.5,])
# nrow(lt[lt$Va==0.5 & lt$Vmg==0,])
# nrow(lt[lt$Va==0.5 & lt$Vmg==1,])
# nrow(lt[lt$Va==0.25 & lt$Vmg==1,])
# nrow(lt[lt$Va==0.125 & lt$Vmg==0.25,])

rels<-sapply(1:nrow(covariances), function(i){
	nrow(lt[lt$Va==covariances$Va[i] & lt$Vmg==covariances$Vmg[i],])
})


## simulate pedigree
## work out number of cousins etc
## work it out from A matrix 
## work it out from inverse A matrix 

ped <- simulate_pedigree(
	years = 5,
	n_females = 1000,
	fecundity = 4,
	p_sire = 0,
	p_polyandry=0.2,
	juv_surv = c(0.5),
	adult_surv = 0,					# discrete generations
	immigration = c(0), 				# closed population
	constant_pop = TRUE     # constant population size
	)$pedigree

## function to match individuals to certain bounds of relatedness - looks up given values of va and vmg, and then also accepts upto (but not including) 2* those relatedness' - this is because inbreeding will lead to slightly higher relatedness' then expected. If the relatedness moves into the next category through high inbreeding then the individuals will not be caught in this category
 rel<-function(x,va,vmg) sum(x$Va>=va & x$Va<2*va & x$Vmg>=vmg & x$Vmg<2*vmg)

 rel2<-function(x,va,vmg,cov) sum(x$Va>=va & x$Va<(2*va) & x$Vmg>=vmg & x$Vmg<(2*vmg) &x$COVamg==cov)


ped_statM <- function(ped){
	dat <- subset(ped,!is.na(dam))
	A <- nadiv::makeA(ped[,1:3])
	A_dam<-A[as.character(dat$dam),as.character(dat$dam)]
	A_id<-A[as.character(dat$animal),as.character(dat$animal)]
	A_id_dam<-A[as.character(dat$animal),as.character(dat$dam)]
	A_id_dam2<-t(A_id_dam) + A_id_dam
		sID_Dam<-Matrix::summary(dd)

	sDam<-Matrix::summary(A_dam)
	names(sDam)[3] <- "Vmg"
	sID<-Matrix::summary(A_id)
	names(sID)[3] <- "Va"
	sID_Dam<-Matrix::summary(A_id_dam2)
	names(sID_Dam)[3] <- "COVamg"
	df1<-as.data.frame(dplyr::left_join(sID,sDam,by=c("i","j")))
	df2<-as.data.frame(dplyr::left_join(df1,sID_Dam,by=c("i","j")))
	df2[is.na(df2[,"Vmg"]),"Vmg"] <- 0
	df2[is.na(df2[,"COVamg"]),"COVamg"] <- 0

	df2[,c("Va","Vmg","COVamg")]<-round(df2[,c("Va","Vmg","COVamg")],6)
	c(
		dam = rel(df2,0.5,0.5),
		FS = rel(df2,0.5, 1),
		MHS =	rel(df2,0.25, 1),
		MG=	rel(df2,0.25, 0.25),
		au_D_FS	=rel(df2,0.25, 0.5),
		au_D_MHS =rel(df2,0.125, 0.25),
		cousin_D_FS	=rel(df2,0.125, 0.5),
		cousin_D_HS	=rel(df2,0.0625, 0.25)
	)

	c(
		dam = rel2(df2,0.5,0.5,1.25),
		FS = rel2(df2,0.5, 1,1),
		MHS =	rel2(df2,0.25, 1,1),
		MG=	rel2(df2,0.25, 0.25,0.625),
		au_D_FS	=rel2(df2,0.25, 0.5,0.75),
		au_D_MHS =rel2(df2,0.125, 0.25,0.25),
		cousin_D_FS	=rel2(df2,0.125, 0.5,0.5),
		cousin_D_HS	=rel2(df2,0.0625, 0.25,0.25)
	)
	# A_id-A_dam
	# lt_id<-round(A_id,6)[lower.tri(round(A_id,6))]
	# lt2_id<-lt_id[lt_id>0]
	# lt_dam<-round(A_dam,6)[lower.tri(round(A_dam,6))]
	# lt2_dam<-lt_dam[lt_dam>0]
	# lt<-data.frame(Va=lt_id,Vmg=lt_dam)
	# lt2<-subset(lt, Va>0)
	# rels<-sapply(1:nrow(covariances), function(i){
	# 	nrow(lt2[lt2$Va==covariances$Va[i] & lt2$Vmg==covariances$Vmg[i],])
	# })
	# names(rels) <- covariances$relationship
	# rels
}

cor(both$Va,both$Vmg)

A_id_dam2[1:10,1:10]
Matrix::summary(A_id_dam)

mat_links(ped)
ped_stat(ped, phenotyped=ped[!is.na(ped$dam),"animal"])
ped_statM(ped)



dat <- subset(ped,!is.na(dam))

A <- nadiv::makeA(ped[,1:3])

#### have I made A_dam right??? - it seems to not be working great?
### maybe its OK?

A_dam<-A[as.character(dat$dam),as.character(dat$dam)]
A_damE<-matrix(as.numeric(A_dam>=1),nrow(A_dam))

damDM<-model.matrix(~dam-1,dat)
colnames(damDM) <- gsub("dam","",colnames(damDM))
rownames(damDM) <- as.character(dat$dam)
A_damE <- damDM[,as.character(dat$dam)]

A_id<-A[as.character(dat$animal),as.character(dat$animal)]
# E <- diag(nrow(dat))
A_cov<-A[as.character(dat$animal),as.character(dat$dam)]

### look with very small pedigree to work out discrepancies

m_rel<-sum(A_dam[lower.tri(A_dam)]>=0.25)
ped_sum <- ped_stat(ped,phenotyped=dat$animal)
c(mat_sibs=sum(ped_sum[c("FS","MHS")]), mat_links=sum(ped_sum[c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")]), other=sum(ped_sum[!names(ped_sum)%in%c( "individuals" ,"links" ,"FS","MHS","dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")]) )

sum(ped_sum[c("dam","MG","FS","MHS","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])
sum(A_dam[lower.tri(A_dam)]>0.125)
nrow(Matrix::summary(A_dam))

sum(A_dam[lower.tri(A_dam)]>0.5)
sum(ped_sum[c("dam","FS","MHS","au_D_FS","cousin_D_FS")])


ped$matriline <- matriline(ped)

table(ped$matriline)
ped_sub <- subset(ped,matriline=="0_33")[-1,1:3]
ped_sub_full <- nadiv::prepPed(ped_sub)
A2 <- nadiv::makeA(ped_sub_full)

A_dam2<-A2[as.character(ped_sub$dam),as.character(ped_sub$dam)]

A_id2<-A2[as.character(ped_sub$animal),as.character(ped_sub$animal)]

ped_sum2 <- ped_stat(ped_sub_full,phenotyped=ped_sub$animal)
dam_sum<-Matrix::summary(A_dam2)
dam_sum2 <- dam_sum[dam_sum[,"i"]!=dam_sum[,"j"],]
nrow(dam_sum2)
sum(round(dam_sum2[,"x"],6)>=0.25)

sum(ped_sum2[c("dam","FS","MHS","au_D_FS","cousin_D_FS")])

sum(ped_sum2[c("dam","MG","FS","MHS","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])

ped_stat2 <- function(ped){
	dat <- ped[!is.na(ped$dam),]
	A <- nadiv::makeA(ped[,1:3])
	A_dam<-A[as.character(dat$dam),as.character(dat$dam)]
	A_id<-A[as.character(dat$animal),as.character(dat$animal)]
	dam_sum<-Matrix::summary(A_dam)
	dam_sum2 <- dam_sum[dam_sum[,"i"]!=dam_sum[,"j"],]
	id_sum<-Matrix::summary(A_id)
	id_sum2 <- id_sum[id_sum[,"i"]!=id_sum[,"j"],]
	c(maternal=nrow(dam_sum2),total=nrow(id_sum2))
}
ped_stat2(ped)

## matriline doesn't completely cover it, because half-sib mothers related through their father would have different matrilines, but their offspring (cousins) could share maternal genetic effects 

A_id[1:10,1:10]
A_dam[1:10,1:10]
A_cov[1:20,1:20]
A_damE[1:10,1:10]
A_dam[3990:4000,3990:4000]
A_damE[3990:4000,3990:4000]

E[1:10,1:10]

cor(as.numeric(A_id),as.numeric(A_dam))
cor(as.numeric(A_id),as.numeric(A_damE))
cor(as.numeric(A_dam),as.numeric(A_damE))
cor(as.numeric(A_dam),as.numeric(E))
cor(as.numeric(A_damE),as.numeric(E))
cor(as.numeric(A_id),as.numeric(E))

betas <- sqrt(c(0.5,0.2,0.3))
mat_cor <- cor(cbind(A_id[lower.tri(A_id)],A_dam[lower.tri(A_dam)],A_damE[lower.tri(A_damE)]))
# mat_cor <- cor(cbind(as.numeric(A_id),as.numeric(A_dam),as.numeric(A_damE)))
# mat_cor <- diag(3)

covs<- betas %*% mat_cor

covs %*% solve(mat_cor)

betas^2
(covs[c(1,3)] %*% solve(mat_cor[c(1,3),c(1,3)]))^2


cA<-chol(0.5 * A_id + sqrt(0) * A_dam + 0.3 * A_damE)

preds <- cA %*% rnorm(nrow(dat)) 
var(preds)
p <- preds + rnorm(nrow(dat),0,sqrt(0.5))

G<-diag(c(0.5,0.5))
Vme <- 0
Ve<-0.5
	g <- rbv0(ped[,1:3],G)
	a <- g[,1]
	mg <- g[match(ped[,2],ped[,1]),2]
	me <- rnorm(nrow(ped),0,sqrt(Vme))[match(ped[,2],ped[,1])]
	e <- rnorm(nrow(ped),0,sqrt(Ve))
	p <- a + mg + me + e
	data <- data.frame(cbind(p=p,ped))
	data$matriline <- as.factor(matriline(ped))

cor(mg,me, use="complete.obs")

# data <- data.frame(cbind(p=p,dat))
	data$mother <- as.factor(data$dam)
	data$mother_PE <- as.factor(data$dam)
	data$animal <- as.factor(data$animal)
	data <- subset(data, !is.na(mother)) 


			assign("ped.ainv", asreml::ainverse(ped), envir = .GlobalEnv) 

mod_1 <- asreml::asreml(
		fixed= p~1
	  , random= ~vm(animal,ped.ainv)
	  , data= data, trace=FALSE)
summary(mod_1)$varcomp

mod_2 <- asreml::asreml(
		fixed= p~1
	  , random= ~vm(animal,ped.ainv) + mother_PE
	  , data= data, trace=FALSE)
summary(mod_2)$varcomp

mod_2b <- asreml::asreml(
		fixed= p~1
	  , random= ~vm(animal,ped.ainv) + matriline
	  , data= data, trace=FALSE)
summary(mod_2b)$varcomp


mod_3 <- asreml::asreml(
		fixed= p~1
	  , random= ~vm(animal,ped.ainv) + mother_PE + matriline
	  , data= data, trace=FALSE)
summary(mod_3)$varcomp

mod_4 <- asreml::asreml(
		fixed= p~1
	  , random= ~vm(animal,ped.ainv) + vm(mother,ped.ainv) + mother_PE
	  , data= data, trace=FALSE)

summary(mod_4)$varcomp

## try out matriline


cM<-Matrix::chol(A_damE)
predsM <- cM %*% rnorm(nrow(dat)) 

diag(A_damE)
isSymmetric(A_damE)
eigen(A_damE, symmetric=TRUE,only.values=TRUE)$values



rm(list=ls())

extract<-FALSE
wd <- "/Users/joelpick/github/maternal_effects/"
main_dir <- paste0(wd,"Data/Raw/Bonnet/")
data_dir <- paste0(wd,"Data/Intermediate/")


source(paste0(wd,"R/extract_cousins.R"))
source(paste0(wd,"R/00_functions.R"))
source("/Users/joelpick/github/squidPed/R/simulate_pedigree.R")

pops<- list.files(main_dir)

i=pops[1]

		pop_dir <- paste0(main_dir,i)
		load(paste0(pop_dir,"/Ainv_",i))
		A <- Matrix::solve(ainv)
		colnames(A) <- rownames(A) <- rownames(ainv)
		rm(ainv)
		dat<-read.csv(paste0(pop_dir,"/data_",i,".csv"))


A_dam<-A[as.character(dat$dam),as.character(dat$dam)]
A_id<-A[as.character(dat$id),as.character(dat$id)]

cor(as.numeric(A_id),as.numeric(A_dam))

rownames(A_id)

sum(round(A_dam,6)[lower.tri(round(A_dam,6))]>0)




ped <- simulate_pedigree(
	years = 5,
	n_females = 200,
	fecundity = 4,
	p_sire = 0,
	p_polyandry=1,
	juv_surv = 0.5,
	adult_surv = 0,					# discrete generations
	immigration = 0, 				# closed population
	constant_pop = TRUE     # constant population size
	)$pedigree

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

m_rel<-sum(A_dam[lower.tri(A_dam)]>=0.25)
ped_sum <- ped_stat(ped,phenotyped=dat$animal)
c(mat_sibs=sum(ped_sum[c("FS","MHS")]), mat_links=sum(ped_sum[c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")]), other=sum(ped_sum[!rownames(ped_sum)%in%c( "individuals" ,"links" ,"FS","MHS","dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")]) )

sum(ped_sum[c("dam","MG","FS","MHS","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])
sum(A_dam[lower.tri(A_dam)]>0)

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

G<-diag(c(0.5,0))
Vme <- 0.3
Ve<-0.5
	g <- rbv0(ped[,1:3],G)
	a <- g[,1]
	mg <- g[match(ped[,2],ped[,1]),2]
	me <- rnorm(nrow(ped),0,sqrt(Vme))[match(ped[,2],ped[,1])]
	e <- rnorm(nrow(ped),0,sqrt(Ve))
	p <- a + mg + me + e
	data <- data.frame(cbind(p=p,ped))

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

mod_3 <- asreml::asreml(
		fixed= p~1
	  , random= ~vm(animal,ped.ainv) + mother_PE
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


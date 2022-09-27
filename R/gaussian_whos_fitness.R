
rm(list=ls())

source("/Users/joelpick/github/pedigree_simulations/R/simulate_pedigree.R")
devtools::load_all("~/github/squidSim/R")


full_sib <- simulate_pedigree(
	years = 5,
	n_females = 100,
	p_breed = 1, 
	fecundity = 4,
	p_sire = 1, 
	p_retain = 0.5,
	juv_surv = 0.5,
	adult_surv = 0,
	immigration = 0,
	constant_pop = TRUE ,
	known_age_structure = FALSE)



half_sib<- simulate_pedigree(
	years = 5,
	n_females = 100,
	p_breed = 1, 
	fecundity = 4,
	p_sire = 0, 
	p_retain = 0.5,
	juv_surv = 0.5,
	adult_surv = 0,
	immigration = 0,
	constant_pop = TRUE ,
	known_age_structure = FALSE)

head(half_sib$pedigree)




# a_hs <- MCMCglmm::rbv(half_sib$pedigree[,1:3],0.5)

a_hs2 <- MCMCglmm::rbv(half_sib$pedigree[,1:3],matrix(c(0.5,0.3,0.3,0.5),nrow=2))

a_hs <- a_hs2[,1]
mg_hs <- a_hs2[match(half_sib$pedigree[,2],half_sib$pedigree[,1]),2]
e_hs <- rnorm(nrow(half_sib$pedigree),0,sqrt(0.5))
p <- a_hs + mg_hs + e_hs

hs_data <- data.frame(cbind(p=p,half_sib$pedigree[,1:3]))

hs_data$mother <- as.factor(hs_data$dam)
hs_data$animal <- as.factor(hs_data$animal)
hs_data<-subset(hs_data, !is.na(mother))
head(hs_data)


library(asreml)
Ainv_hs<-ainverse(half_sib$pedigree[,1:3])


mod_hs1 <- asreml(fixed= p~1
                   , random= ~vm(animal,Ainv_hs) + mother
                   , data= hs_data
                 	)
summary(mod_hs1)$varcomp


mod_hs2 <- asreml(fixed= p~1
                   , random= ~vm(mother,Ainv_hs)
                   , data= hs_data
                 	)
summary(mod_hs2)$varcomp

mod_hs3 <- asreml(fixed= p~1
                   , random= ~vm(mother,Ainv_hs) + vm(animal,Ainv_hs)
                   , data= hs_data
                 	)
summary(mod_hs3)$varcomp

mod_hs3 <- asreml(fixed= p~1
                   , random= ~str(~vm(mother,Ainv_hs) + vm(animal,Ainv_hs),~us(2):vm(animal,Ainv_hs))
                   , data= hs_data
                 	)
summary(mod_hs3)$varcomp

mge_analysis <- function(data,ped){
	ped.ainv1<-ainverse(ped)
	assign("ped.ainv", ped.ainv1, envir = .GlobalEnv) 

	mod_hs1 <- asreml(fixed= p~1
	                   , random= ~vm(animal,ped.ainv)
	                   , data= data
	                 	)

	mod_hs2 <- asreml(fixed= p~1
	                   , random= ~vm(animal,ped.ainv) + mother
	                   , data= data
	                 	)

	mod_hs3 <- asreml(fixed= p~1
	                   , random= ~vm(mother,ped.ainv)
	                   , data= data
	                 	)

	mod_hs4 <- asreml(fixed= p~1
	                   , random= ~vm(mother,ped.ainv) + mother
	                   , data= data
	                 	)

	mod_hs5 <- asreml(fixed= p~1
	                   , random= ~vm(mother,ped.ainv) + vm(animal,ped.ainv)
	                   , data= data
	                 	)

	mod_hs6 <- asreml(fixed= p~1
	                   , random= ~str(~vm(animal,ped.ainv) +vm(mother,ped.ainv) ,~us(2):vm(animal,ped.ainv))
	                   , data= data
	                 	)

rm("ped.ainv")

	m1<-summary(mod_hs1)$varcomp
	m2<-summary(mod_hs2)$varcomp
	m3<-summary(mod_hs3)$varcomp
	m4<-summary(mod_hs4)$varcomp
	m5<-summary(mod_hs5)$varcomp
	m6<-summary(mod_hs6)$varcomp

	cbind(
		A = c(m1["vm(animal, ped.ainv)",1], m2["vm(animal, ped.ainv)",1], NA, NA, m5["vm(animal, ped.ainv)",1], m6[1,1]),
		Me = c(NA, m2["mother",1],NA,m4["mother",1],NA,NA),
		Mg = c(NA, NA, m3["vm(mother, ped.ainv)",1], m4["vm(mother, ped.ainv)",1], m5["vm(mother, ped.ainv)",1], m6[3,1]),
		cov_AMg = c(NA,NA,NA,NA,NA, m6[2,1]),
		E = c(m1["units!R",1],m2["units!R",1],m3["units!R",1],m4["units!R",1],m5["units!R",1],m6["units!R",1])
	)
}

out <- mge_analysis(hs_data,half_sib$pedigree[,1:3])

plot(NA,ylim=c(0,1),xlim=c(1,6))
sapply(1:6,function(x) points(rep(x,5)+(1:5)/10,out[x,], pch=19, col=1:5))
abline(h=0.5)


0.5 + 3/2*0.3 + 1/4*0.5

0.5 + 1/4*0.5 +  sqrt(0.25) * 0.3*2

x<-matrix(c(0.5,0.3,0.3,0.5),nrow=2)
y<-diag(sqrt(c(1,0.25)))

sum(y%*%x%*%y)




Ainv_hs<-MCMCglmm::inverseA(half_sib$pedigree[,1:3])$Ainv

mod_hs1 <- MCMCglmm::MCMCglmm(p~1, random=~ animal,data=hs_data,ginverse=list(animal=Ainv_hs))
summary(mod_hs1)

mod_hs2 <- MCMCglmm::MCMCglmm(p~1, random=~ dam,data=hs_data,ginverse=list(dam=Ainv_hs))
summary(mod_hs2)



a_fs <- MCMCglmm::rbv(full_sib$pedigree[,1:3],0.5)
e_fs <- rnorm(nrow(full_sib$pedigree),0,sqrt(0.5))
p<-a_fs+e_fs
fs_data <- data.frame(cbind(p=p,full_sib$pedigree[,1:3]))
fs_data$mother <- fs_data$dam
head(fs_data)

Ainv_fs<-MCMCglmm::inverseA(full_sib$pedigree[,1:3])$Ainv

mod_fs1 <- MCMCglmm::MCMCglmm(p~1, random=~ animal,data=fs_data,ginverse=list(animal=Ainv_fs))
summary(mod_fs1)

mod_fs2 <- MCMCglmm::MCMCglmm(p~1, random=~ dam,data=fs_data,ginverse=list(dam=Ainv_fs))
summary(mod_fs2)


se <- function(x) sd(x)/sqrt(length(x))


#####
#--- relatedness functions
####

maternal_relatedness_ped <- function(ped){ 
	ped<-na.omit(ped) ## might depend on whether unknown males paternity likely means unknown male - in which case excluding will increase full-sib
	c(sapply(split(ped,ped$dam), function(x){
		if(nrow(x)==1){
		 NULL
		}else{
			fs <- NULL
			for(i in 1:(nrow(x)-1)){
				fs <- c(fs, x$sire[i] == x$sire[(i+1):nrow(x)])
			}
			mean(0.25 + fs*0.25)
		}
	}), recursive=TRUE)
}

## if pedigree is big, will take a long time/crash
maternal_relatedness_ainv <- function(ainv,dat){ 
	# ped<-na.omit(ped) ## might depend on whether unknown males paternity likely means unknown male - in which case excluding will increase full-sib
	A <- Matrix::solve(ainv)
	# A <- chol2inv(chol(ainv))
	colnames(A) <- rownames(A) <- rownames(ainv)
	c(sapply(split(dat,dat$dam), function(x){
		if(nrow(x)==1){
		 NULL
		}else{
			mean(A[as.character(x$id),as.character(x$id)][lower.tri(A[x$id,x$id])])
		}
	}), recursive=TRUE)
}

maternal_relatedness_A <- function(A,dat){ 
	c(parallel::mclapply(split(dat,dat$dam), function(x){
		if(nrow(x)==1){
		 NULL
		}else{
			sub_A <-A[as.character(x$id),as.character(x$id)]
			mean(sub_A[which(lower.tri(sub_A))])
		}
	},mc.cores=8), recursive=TRUE)
}


matriline <- function(ped){
	matriline <- vector(mode="character",length=nrow(ped))
	for(i in 1:nrow(ped)){
		matriline[i] <- ifelse(is.na(ped[i,"dam"]), ped[i,"animal"],matriline[match(ped[i,"dam"],ped[,"animal"])])
	}
	matriline
}


pedCor <- function(ped){
	dat <- subset(ped,!is.na(dam))
	A <- nadiv::makeA(ped[,1:3])
	A_dam<-A[as.character(dat$dam),as.character(dat$dam)]
	# A_damE<-matrix(as.numeric(A_dam>=1),nrow(A_dam))
	damDM<-	Matrix::fac2sparse(dat$dam)
	# colnames(damDM) <- gsub("dam","",colnames(damDM))
	colnames(damDM) <- as.character(dat$dam)
	A_damE <- as(damDM[as.character(dat$dam),], "symmetricMatrix")
	A_id<-A[as.character(dat$animal),as.character(dat$animal)]
	E <- diag(nrow(dat))
	A_cov<-A[as.character(dat$animal),as.character(dat$dam)]

	cor(cbind(A_id[lower.tri(A_id)],A_dam[lower.tri(A_dam)],A_damE[lower.tri(A_damE)]))
}

	# nrow(Matrix::summary(A_damE))-nrow(dat)

ped_stat2 <- function(ped){
	colnames(ped)[1]<-"animal"
	dat <- ped[!is.na(ped$dam),]
	A <- nadiv::makeA(ped[,1:3])
	A_dam<-A[as.character(dat$dam),as.character(dat$dam)]
	A_id<-A[as.character(dat$animal),as.character(dat$animal)]
	dam_sum<-Matrix::summary(A_dam)
	id_sum<-Matrix::summary(A_id)
	mat_sibs<- table(ped$dam)
	c(mat_sib = sum(mat_sibs * (mat_sibs-1)/2),
		mat_links = nrow(dam_sum)-nrow(dat),
		total_links = nrow(id_sum)-nrow(dat))
}


rbv0 <- function(pedigree, G){
	X <- matrix(0, nrow=nrow(pedigree), ncol=nrow(G))
	index <- which(diag(G)!=0)
	if(any(diag(G)==0)) G <- G[index,index]
	if(length(index)>0){
		X2 <- MCMCglmm::rbv(pedigree=pedigree, G=G)
		X[,index] <- X2	
	}
	X
}

 mge_sim <- function(ped,param, Vp=1){
 	
 	colnames(ped) <- c("animal","dam","sire")

 	Va <- param["Va"]
 	Vmg <- param["Vmg"]
 	r_amg <- param["r_amg"]
 	Vme <- param["Vme"]

 	D <- diag(sqrt(c(Va,Vmg)))
	R <- matrix(c(1,r_amg,r_amg,1),nrow=2)
	G <- D%*%R%*%D

	Ve <- Vp - sum(G) - Vme

	## simulate direct, maternal genetic and environmental effects, add together to make phenotype
	
	g <- rbv0(ped[,1:3],G)
	a <- g[,1]
	mg <- g[match(ped[,2],ped[,1]),2]
	me <- rnorm(nrow(ped),0,sqrt(Vme))[match(ped[,2],ped[,1])]
	e <- rnorm(nrow(ped),0,sqrt(Ve))
	p <- a + mg + me + e
	data <- data.frame(cbind(p=p,ped))
	data$mother <- as.factor(data$dam)
	data$mother_PE <- as.factor(data$dam)
	data$animal <- as.factor(data$animal)
	data <- subset(data, !is.na(mother)) 
	data
 }


## paternal genetic effects 
pge_sim <- function(ped, Va, Vmg, Vpg, R, Vme, Vpe, Vp=1){
 	D <- diag(sqrt(c(Va,Vmg,Vpg)))
	# R <- matrix(c(1,r_amg,r_amg,1),nrow=2)
	G <- D%*%R%*%D

	Ve <- Vp - sum(G) - Vme

	## simulate direct, maternal genetic and environmental effects, add together to make phenotype
	
	g <- rbv0(ped[,1:3],G)
	a <- g[,1]
	mg <- g[match(ped[,2],ped[,1]),2]
	## need to be linked to social father
	pg <- g[match(ped[,"social_sire"],ped[,1]),3]

	me <- rnorm(nrow(ped),0,sqrt(Vme))[match(ped[,2],ped[,1])]
	e <- rnorm(nrow(ped),0,sqrt(Ve))
	p <- a + mg + me + e
	data <- data.frame(cbind(p=p,ped))
	data$mother <- as.factor(data$dam)
	data$mother_PE <- as.factor(data$dam)
	data$animal <- as.factor(data$animal)
	data <- subset(data, !is.na(mother)) 
	data
 }


  ## Model 1 - additive genetic effects, assuming phenotype is trait of offspring
m1_func <- function(data){
	mod_1 <- asreml(
		fixed= p~1
		, random= ~vm(animal,ped.ainv)
	  , data= data, trace=FALSE)
	m1 <- summary(mod_1)$varcomp
	c(A= m1["vm(animal, ped.ainv)",1], 
		Me=NA, 
		Mg=NA,
		cov_AMg =NA,
		E = m1["units!R",1])
}

m1a_func <- function(data){
	mod_1 <- asreml(
		fixed= p~1
		, random= ~mother_PE
	  , data= data, trace=FALSE)
	m1 <- summary(mod_1)$varcomp
	c(A= NA, 
		Me=m1["mother_PE",1], 
		Mg=NA,
		cov_AMg =NA,
		E = m1["units!R",1])
}

  ## Model 2 - additive genetic effects with maternal environment effects, assuming phenotype is trait of offspring
m2_func <- function(data){
	mod_2 <- asreml(
		fixed= p~1
    , random= ~vm(animal,ped.ainv) + mother_PE
    , data= data, trace=FALSE)
	m2<-summary(mod_2)$varcomp
	m2_sum <- c(
		A= m2["vm(animal, ped.ainv)",1], 
		Me=m2["mother_PE",1], 
		Mg=NA,
		cov_AMg =NA,
		E = m2["units!R",1])
}

m3_func <- function(data){
	mod_3 <- asreml(
		fixed= p~1
		, random= ~vm(animal,ped.ainv) + vm(mother,ped.ainv) 
		, data= data, trace=FALSE)
	m3<-summary(mod_3)$varcomp
	m3_sum <- c(
		A= m3["vm(animal, ped.ainv)",1], 
		Me=NA, 
		Mg=m3["vm(mother, ped.ainv)",1],
		cov_AMg =NA,
		E = m3["units!R",1])
}

m4_func <- function(data){
  mod_4 <- asreml(
		fixed= p~1
	  , random= ~vm(animal,ped.ainv) + vm(mother,ped.ainv) + mother_PE
	  , data= data, trace=FALSE)
	m4<-summary(mod_4)$varcomp
	m4_sum <- c(
		A= m4["vm(animal, ped.ainv)",1], 
		Me=m4["mother_PE",1], 
		Mg=m4["vm(mother, ped.ainv)",1],
		cov_AMg =NA,
		E = m4["units!R",1])
}

m5_func <- function(data){
	mod_5 <- asreml(
		fixed= p~1
    , random= ~str(~vm(animal,ped.ainv) +vm(mother,ped.ainv) ,~us(2):vm(animal,ped.ainv)) + mother_PE
    , data= data, trace=FALSE,maxit=50)
	m5<-summary(mod_5)$varcomp
	m5_sum <- c(
		A= m5[1,1], 
		Me=m5["mother_PE",1], 
		Mg=m5[3,1],
		cov_AMg =m5[2,1],
		E = m5["units!R",1])
}	

	# 	assign("ped.ainv", asreml::ainverse(hs_peds[[2]]), envir = .GlobalEnv) 
	# 	out <- do.call(rbind,parallel::mclapply(hs_data[[2]], m5_func, mc.cores=6))
	# 	rm("ped.ainv", envir = .GlobalEnv)
	# 	out
 
	# mod_5 <- asreml(
	# 	fixed= p~1
 #    , random= ~str(~vm(animal,ped.ainv) +vm(mother,ped.ainv) ,~us(2):vm(animal,ped.ainv)) + mother_PE
 #    , data= hs_data[[2]][[3]], trace=FALSE,maxit=50)


m6_func <- function(data){
	mod_6 <- asreml(
		fixed= p~1
    , random= ~vm(mother,ped.ainv)
    , data= data, trace=FALSE)
	m6<-summary(mod_6)$varcomp
	m6_sum <- c(
		A= NA, 
		Me=NA, 
		Mg=m6["vm(mother, ped.ainv)",1],
		cov_AMg =NA,
		E = m6["units!R",1])
}	


m7_func <- function(data){
	mod_7 <- asreml(
		fixed= p~1
	  , random= ~vm(mother,ped.ainv) + mother_PE
	  , data= data, trace=FALSE)
	m7<-summary(mod_7)$varcomp
	m7_sum <- c(
		A= NA, 
		Me=m7["mother_PE",1], 
		Mg=m7["vm(mother, ped.ainv)",1],
		cov_AMg =NA,
		E = m7["units!R",1])
}	


m8_func <- function(data){
	data_means <- aggregate(p~mother,data,mean)
	mod_8 <- asreml(
		fixed= p~1
    , random= ~vm(mother,ped.ainv)
    , data= data_means, trace=FALSE)
	m8<-summary(mod_8)$varcomp
	m8_sum <- c(
		A= NA, 
		Me=NA, 
		Mg=m8["vm(mother, ped.ainv)",1],
		cov_AMg =NA,
		E = m8["units!R",1])
}	

model_func <- function(FUN,peds,data,mc.cores=8){
	mclapply(1:length(data), function(i){
		if(is.list(peds)&!is.data.frame(peds)) { 
			ped <- peds[[i]] 
		}else{
			ped <- peds
		}
		colnames(ped) <- c("animal","dam","sire")
		assign("ped.ainv", asreml::ainverse(ped), envir = .GlobalEnv) 
		out <- do.call(rbind,lapply(data[[i]], FUN))
		rm("ped.ainv", envir = .GlobalEnv)
		cat(i," ")
		out
	}, mc.cores=mc.cores)
}



list2array <- function(x) array(unlist(x), dim = c(nrow(x[[1]]), ncol(x[[1]]), length(x)), dimnames=list(NULL,c("A","Me","Mg","cov_AMg","E"),NULL))

change2zero <- function(x) ifelse(is.na(x), 0, x)

total_va <- function(x, j, m){
	out<-rep(NA, nrow(x))
	y <- change2zero(x)
	out[j]<- y[j,"A"] + 3/2 * y[j,"cov_AMg"] + 1/2 * y[j,"Mg"]
	out[m]<-  y[m,"Mg"]
	out
}

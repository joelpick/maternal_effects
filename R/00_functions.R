
########################
# generic functions
######################## 


####
#--- function for standard error
####

se <- function(x) sd(x)/sqrt(length(x))


####
#--- function to change NAs to 0s
####

change2zero <- function(x) ifelse(is.na(x), 0, x)


####
#--- rbind that doesnt mind abut names not matching, and takes the first set of names, or a given set of names
####
rbind_notAnnoying <- function(..., names=NULL){
  x <- list(...)
  y <- lapply(x, function(y){
    names(y) <- if(is.null(names)) names(x[[1]]) else names
    return(y)  
  } )
  do.call(rbind,y)
}


########################
# RELATEDNESS FUNCTIONS
######################## 

#################################################

## function to make products of all combinations and sum them 
combo_prod <- function(x) if(length(x)>1){sum( utils::combn(x, m =2)[1, ] * utils::combn(x, m =2)[2, ])}else{0}


combo_au <- function(x) if(length(x)>1){sum(x * (length(x)-1))}else{0}

n_cousin <- function(p,gp,ped){
	if(all(is.na(ped[,gp]))){ ## stops working when there no links through a certain grandparent type
		0
	}else{
		f1 <- formula(paste("animal~",p,"+",gp))
	  f2 <- formula(paste("animal~",gp))
	  d1 <- aggregate(f1, ped,length)
		sum(aggregate(f2,d1,combo_prod)$animal)	
	}
}


n_au <- function(link,ped){
	if(all(is.na(ped[,link]))){
		0
	}else{	
		p <- if(substr(link,1,1)=="M"){"dam"}else{"sire"}
		gp <- if(substr(link,3,3)=="M"){"dam"}else if(substr(link,3,3)=="F"){"sire"} else{"pair"}
	  f1 <- formula(paste("animal~",p,"+",link))
	  f2 <- formula(paste("animal~",gp))
	  d1 <- aggregate(f1, ped,length)
	  d2 <- aggregate(f2, ped,length)
		sum(d1$animal * (d2[match(d1[,link],d2[,gp]),"animal"]-1), na.rm=TRUE)
	}
}

total_links <- function(ped) nrow(ped) * (nrow(ped) - 1) / 2

non_zero_links <- function(ped){
	# pedA<-nadiv::makeA(ped[,1:3])
	# sum(pedA[lower.tri(pedA)]>0)
nrow(Matrix::summary(nadiv::makeA(ped[,1:3])))
}


# ped<-ped_sub_full

ped_stat <- function(ped, phenotyped=NULL){	
	colnames(ped) <- c("animal","dam","sire")

	# dummy_dam <- which(is.na(ped$dam) & !is.na(ped$sire))
	# dummy_sire <- which(is.na(ped$sire) & !is.na(ped$dam))

	# ped[dummy_sire,"sire"] <- paste0("ds_",seq_along(dummy_sire))
	# ped[dummy_dam,"dam"] <- paste0("dd_",seq_along(dummy_dam))


	## parent pairs
	ped$pair<-paste(ped$dam,ped$sire)
	ped$pair<-ifelse(ped$pair== "NA NA", NA, ped$pair)

  ## grandparents
	ped$MGM <- ped[match(ped[,2],ped[,1]),2]
	ped$MGF <- ped[match(ped[,2],ped[,1]),3]
	ped$PGM <- ped[match(ped[,3],ped[,1]),2]
	ped$PGF <- ped[match(ped[,3],ped[,1]),3]

	## grandparent pairs
	ped$MGP <- paste(ped$MGM,ped$MGF)
	ped$PGP <- paste(ped$PGM,ped$PGF)
	# ped$MGP<-ifelse(ped$MGP== "NA NA", NA, ped$MGP)
	# ped$PGP<-ifelse(ped$PGP== "NA NA", NA, ped$PGP)
	ped$MGP<-ifelse(grepl("NA",ped$MGP), NA, ped$MGP)
	ped$PGP<-ifelse(grepl("NA",ped$PGP), NA, ped$PGP)

	## if dont specify phenotyped, then assume all phenotyped
	if(is.null(phenotyped)) phenotyped <- ped[,1]
	## subset to phenotyped individuals
	ped2 <- subset(ped,animal %in% phenotyped)	


	####################
	# ---- parents and grandparents
	####################

	gp <- apply(ped[ped[,1] %in% phenotyped,c("dam","sire","MGM","MGF","PGM","PGF")],2,function(x) sum(x%in% phenotyped) )


	####################
	# ---- sibs
	####################

	mat_sib <- sum(table(ped2$dam) * (table(ped2$dam)-1)/2)
	pat_sib <- sum(table(ped2$sire) * (table(ped2$sire)-1)/2)
	ped2$pair <- paste(ped2$dam,ped2$sire)
	n_pair <-table(ped2$pair[!grepl("NA",ped2$pair)])
	FS <- sum(n_pair * (n_pair-1)/2)
	sibs <- c(FS=FS, MHS=mat_sib-FS, PHS=pat_sib-FS)



	####################
	# ---- Cousins
	####################
		## stack maternal and paternal grandparents, to get allows for all cousins relationships, not just through mothers or father
	GM <- rbind_notAnnoying(ped2[,c("animal","dam","MGM","MGF","MGP")], ped2[,c("animal","sire","PGM","PGF","PGP")])

	cousin_D_FS <- n_cousin(p="dam", gp="MGP", ped=ped2)
	cousin_D_MHS <- n_cousin(p="dam", gp="MGM", ped=ped2) - cousin_D_FS
	cousin_D_PHS <- n_cousin(p="dam", gp="MGF", ped=ped2) - cousin_D_FS

	cousin_S_FS <- n_cousin(p="sire", gp="PGP", ped=ped2)
	cousin_S_MHS <- n_cousin(p="sire", gp="PGM", ped=ped2) - cousin_S_FS
	cousin_S_PHS <- n_cousin(p="sire", gp="PGF", ped=ped2) - cousin_S_FS


	cousin_FS <- n_cousin(p="dam", gp="MGP", ped=GM)
	cousin_DS_FS <- cousin_FS-cousin_D_FS-cousin_S_FS

	cousin_MS <- n_cousin(p="dam", gp="MGM", ped=GM)
	cousin_DS_MHS <- cousin_MS-cousin_FS-cousin_D_MHS-cousin_S_MHS

	cousin_PS <- n_cousin(p="dam", gp="MGF", ped=GM)
	cousin_DS_PHS <- cousin_PS-cousin_FS-cousin_D_PHS-cousin_S_PHS

	# cousins <- c(cousin_D_FS=cousin_D_FS,
	# 	cousin_DS_FS=cousin_DS_FS,
	# 	cousin_S_FS=cousin_S_FS,
	# 	cousin_D_MHS=cousin_D_MHS,
	# 	cousin_DS_MHS=cousin_DS_MHS,
	# 	cousin_S_MHS=cousin_S_MHS,
	# 	cousin_D_PHS=cousin_D_PHS,
	# 	cousin_DS_PHS=cousin_DS_PHS,
	# 	cousin_S_PHS=cousin_S_PHS)

	cousins <- c(cousin_D_FS=cousin_D_FS,
		cousin_DS_FS=cousin_DS_FS,
		cousin_S_FS=cousin_S_FS,
		cousin_D_HS=cousin_D_MHS + cousin_D_PHS,
		cousin_DS_HS=cousin_DS_MHS + cousin_DS_PHS,
		cousin_S_HS=cousin_S_MHS + cousin_S_PHS)


####################
# ---- aunts and uncles
####################

	# related through a maternal grandmother (au_D_FS + au_D_MHS ?)
	au_D_MS <- n_au("MGM",ped2)

	# related through a maternal grandfather (au_D_FS + au_D_PHS ?)
	au_D_PS <- n_au("MGF",ped2)

	# related through a maternal grandparents (au_D_FS ?)
	au_D_FS <- n_au("MGP",ped2)

	au_D_PHS <- au_D_PS - au_D_FS
	au_D_MHS <- au_D_MS - au_D_FS

	# related through a paternal grandmother (au_S_FS + au_S_MHS ?)
	au_S_MS <- n_au("PGM",ped2)

	# related through a paternal grandfather (au_S_FS + au_D_PHS ?)
	au_S_PS <- n_au("PGF",ped2)

	# related through a paternal grandparents (au_S_FS ?)
	au_S_FS <- n_au("PGP",ped2)

	au_S_PHS <- au_S_PS - au_S_FS
	au_S_MHS <- au_S_MS - au_S_FS

	au <- c(au_D_FS=au_D_FS,au_S_FS=au_S_FS,au_D_MHS=au_D_MHS,au_S_MHS=au_S_MHS,au_D_PHS=au_D_PHS,au_S_PHS=au_S_PHS)

		####################
		# ---- put together
		####################


	stat <- c(
		individuals = length(phenotyped),
		links = length(phenotyped) * (length(phenotyped) - 1) / 2,
		gp[c("dam","sire")],
		MG=as.vector(gp["MGM"] + gp["MGF"]), PG=as.vector(gp["PGM"] + gp["PGF"]),
		sibs,au,cousins
	)
	all_rel <- c("individuals","links","dam","sire","FS","MHS","PHS","MG","PG","au_D_FS","au_S_FS","au_D_MHS","au_S_MHS","au_D_PHS","au_S_PHS","cousin_D_FS","cousin_DS_FS","cousin_S_FS","cousin_D_HS","cousin_DS_HS","cousin_S_HS")#,"cousin_D_FS","cousin_DS_FS","cousin_S_FS","cousin_D_MHS","cousin_DS_MHS","cousin_S_MHS","cousin_D_PHS","cousin_DS_PHS","cousin_S_PHS"

	out<-rep(0,length(all_rel))
	names(out)<-all_rel

	out[names(stat)]<-stat
	out
}





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



########################
# SIMULATION FUNCTIONS
######################## 


####
#--- modification of MCMCglmms rbv function for simulating breeding values, that accepts 0 variance 
####
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

####
#--- simulate maternal genetic effects
####
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


########################
# MODELLING FUNCTIONS
######################## 

####
#--- gets the sampling covariance matrix for the variance components from an asreml model. Note that it needs residual = ~idv(units) to have been run in the model to work
####
inv_hessian_varcomp <- function(mod){
	inv_hess <- mod$ai
	inv_hess <- as.matrix(inv_hess[rownames(inv_hess)!="units!R",colnames(inv_hess)!="units!R"])
	ih_names <- rownames(inv_hess)
	ih_names <- sub("!.+$", "", ih_names)
	ih_names <- sub(".*\\((\\w*)\\,.+\\).*$", "\\1", ih_names)
	colnames(inv_hess)<-rownames(inv_hess)<-ih_names
	return(inv_hess)
}


## Model 1 - additive genetic effects, assuming phenotype is trait of offspring
m1_func <- function(data){
	mod <- asreml(
		fixed= p~1
		, random= ~vm(animal,ped.ainv)
	  , residual = ~idv(units)
    , data= data, trace=FALSE)
	m1 <- summary(mod)$varcomp
	list(
		samp_cov = inv_hessian_varcomp(mod),
		ml = c(
			A= m1["vm(animal, ped.ainv)",1], 
			Me=NA, 
			Mg=NA,
			cov_AMg =NA,
			E = m1["units!R",1])
	)
	
}

m1a_func <- function(data){
	mod <- asreml(
		fixed= p~1
		, random= ~mother_PE
	  , residual = ~idv(units)
    , data= data, trace=FALSE)
	m1 <- summary(mod)$varcomp
	list(
		samp_cov = inv_hessian_varcomp(mod),
		ml= c(
			A= NA, 
			Me=m1["mother_PE",1], 
			Mg=NA,
			cov_AMg =NA,
			E = m1["units!R",1])
		)
}

  ## Model 2 - additive genetic effects with maternal environment effects, assuming phenotype is trait of offspring
m2_func <- function(data){
	mod <- asreml(
		fixed= p~1
    , random= ~vm(animal,ped.ainv) + mother_PE
    , residual = ~idv(units)
    , data= data, trace=FALSE)
	m2<-summary(mod)$varcomp
	list(
		samp_cov = inv_hessian_varcomp(mod),
		ml= c(
			A= m2["vm(animal, ped.ainv)",1], 
			Me=m2["mother_PE",1], 
			Mg=NA,
			cov_AMg =NA,
			E = m2["units!R",1])
		)
}

m3_func <- function(data){
	mod <- asreml(
		fixed= p~1
		, random= ~vm(animal,ped.ainv) + vm(mother,ped.ainv) 
		, residual = ~idv(units)
    , data= data, trace=FALSE)
	m3<-summary(mod)$varcomp
	list(
		samp_cov = inv_hessian_varcomp(mod),
		ml= c(
			A= m3["vm(animal, ped.ainv)",1], 
			Me=NA, 
			Mg=m3["vm(mother, ped.ainv)",1],
			cov_AMg =NA,
			E = m3["units!R",1])
		)
}

m4_func <- function(data){
  mod <- asreml(
		fixed= p~1
	  , random= ~vm(animal,ped.ainv) + vm(mother,ped.ainv) + mother_PE
	  , residual = ~idv(units)
    , data= data, trace=FALSE)
	m4<-summary(mod)$varcomp
	list(
		samp_cov = inv_hessian_varcomp(mod),
		ml= c(
			A= m4["vm(animal, ped.ainv)",1], 
			Me=m4["mother_PE",1], 
			Mg=m4["vm(mother, ped.ainv)",1],
			cov_AMg =NA,
			E = m4["units!R",1])
		)
}

m5_func <- function(data){
	mod <- asreml(
		fixed= p~1
    , random= ~str(~vm(animal,ped.ainv) +vm(mother,ped.ainv) ,~us(2):vm(animal,ped.ainv)) + mother_PE
    , residual = ~idv(units)
    , data= data, trace=FALSE,maxit=50)
	m5<-summary(mod)$varcomp
	list(
		samp_cov = inv_hessian_varcomp(mod),
		ml= c(
			A= m5[1,1], 
			Me=m5["mother_PE",1], 
			Mg=m5[3,1],
			cov_AMg =m5[2,1],
			E = m5["units!R",1])
		)
}	


m6_func <- function(data){
	mod <- asreml(
		fixed= p~1
    , random= ~vm(mother,ped.ainv)
    , residual = ~idv(units)
    , data= data, trace=FALSE)
	m6<-summary(mod)$varcomp
	list(
		samp_cov = inv_hessian_varcomp(mod),
		ml= c(
			A= NA, 
			Me=NA, 
			Mg=m6["vm(mother, ped.ainv)",1],
			cov_AMg =NA,
			E = m6["units!R",1])
		)
}	


m7_func <- function(data){
	mod <- asreml(
		fixed= p~1
	  , random= ~vm(mother,ped.ainv) + mother_PE
	  , residual = ~idv(units)
    , data= data, trace=FALSE)
	m7<-summary(mod)$varcomp
	list(
		samp_cov = inv_hessian_varcomp(mod),
		ml= c(
			A= NA, 
			Me=m7["mother_PE",1], 
			Mg=m7["vm(mother, ped.ainv)",1],
			cov_AMg =NA,
			E = m7["units!R",1])
		)
}	


m8_func <- function(data){
	data_means <- aggregate(p~mother,data,mean)
	mod <- asreml(
		fixed= p~1
    , random= ~vm(mother,ped.ainv)
    , residual = ~idv(units)
    , data= data_means, trace=FALSE)
	m8<-summary(mod)$varcomp
	list(
		samp_cov = inv_hessian_varcomp(mod),
		ml= c(
			A= NA, 
			Me=NA, 
			Mg=m8["vm(mother, ped.ainv)",1],
			cov_AMg =NA,
			E = m8["units!R",1])
		)
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
		out <- lapply(data[[i]], FUN)
		out2 <- list(
			ml=do.call(rbind,lapply(out,function(x) x[["ml"]])),
			samp_cov=list2array(lapply(out,function(x) x[["samp_cov"]]))
		)
		rm("ped.ainv", envir = .GlobalEnv)
		cat(i," ")
		out2
	}, mc.cores=mc.cores)
}


list2array <- function(x) array(unlist(x), dim = c(nrow(x[[1]]), ncol(x[[1]]), length(x)), dimnames=list(rownames(x[[1]]),colnames(x[[1]]),NULL))

# list2array <- function(x) array(unlist(x), dim = c(nrow(x[[1]]), ncol(x[[1]]), length(x)), dimnames=list(NULL,c("A","Me","Mg","cov_AMg","E"),NULL))


total_va <- function(x, j, m){
	out<-rep(NA, nrow(x))
	y <- change2zero(x)
	out[j]<- y[j,"A"] + 3/2 * y[j,"cov_AMg"] + 1/2 * y[j,"Mg"]
	out[m]<-  y[m,"Mg"]
	out
}

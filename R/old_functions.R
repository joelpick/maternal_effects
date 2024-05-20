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


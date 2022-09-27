
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

 mge_sim <- function(ped,Va, Vmg,r_amg, Vme, Vp=1){
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


list2array <- function(x) array(unlist(x), dim = c(nrow(x[[1]]), ncol(x[[1]]), length(x)), dimnames=list(NULL,c("A","Me","Mg","cov_AMg","E"),NULL))

change2zero <- function(x) ifelse(is.na(x), 0, x)

total_va <- function(x, j, m){
	out<-rep(NA, nrow(x))
	y <- change2zero(x)
	out[j]<- y[j,"A"] + 3/2 * y[j,"cov_AMg"] + 1/2 * y[j,"Mg"]
	out[m]<-  y[m,"Mg"]
	out
}

mat.summ   <- Matrix::summary(A)

lower.summ <- subset(mat.summ, i != j)
sum(mat.summ$i != mat.summ$j)
head(mat.summ[order(mat.summ$i),])

full_ind <- which(upper.tri(matrix(0,nrow=nrow(dat),ncol=nrow(dat))), arr.ind=TRUE)

x<-lower.summ[,"x"][match(paste(full_ind[,"row"],full_ind[,"col"]),paste(lower.summ[,"i"],lower.summ[,"j"]))]

### check whether damA are all within A, because then could just do with with in the sparse summary of A

x<-pedCor(ped)
x[lower.tri(x)]

pedCor <- function(ped){
	dat <- subset(ped,!is.na(dam))
	A <- nadiv::makeA(ped[,1:3])
	A_dam<-A[as.character(dat$dam),as.character(dat$dam)]
	A_damE<-matrix(as.numeric(A_dam>=1),nrow(A_dam))
	damDM<-model.matrix(~dam-1,dat)
	colnames(damDM) <- gsub("dam","",colnames(damDM))
	rownames(damDM) <- as.character(dat$dam)
	A_damE <- damDM[,as.character(dat$dam)]
	A_id<-A[as.character(dat$animal),as.character(dat$animal)]
	E <- diag(nrow(dat))
	A_cov<-A[as.character(dat$animal),as.character(dat$dam)]

	cor(cbind(A_id[lower.tri(A_id)],A_dam[lower.tri(A_dam)],A_damE[lower.tri(A_damE)]))
}

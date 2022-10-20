rm(list=ls())

library(viridis)

wd <- "/Users/joelpick/github/maternal_effects/"
data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))
# load(paste0(data_wd,"gaussian_data.Rdata"))
load(paste0(data_wd,"gaussian_sims.Rdata"))



comp_names <- c("A","Me","Mg","cov_AMg","E")
model_names <- c(
	"A",
	"A + Me",
	"A + Mg",
  "A + Mg + Me",
  "A + Mg + cov + Me",
  "Mg",
  "Mg + Me",
  "Mg average"
	)

scenarios2 <- scenarios[,c(1,4,2,3)]
scenarios2[,"r_amg"] <- scenarios2[,"r_amg"]*sqrt(scenarios2[,"Va"])*sqrt(scenarios2[,"Vmg"])
scenarios2 <- cbind(scenarios2,1 - rowSums(scenarios2) - scenarios2[,"r_amg"])
colnames(scenarios2) <- c("A","Me","Mg","cov_AMg","E")


model5_fs <- lapply(model5_fs, function(x){
  x[x=="Error in asreml(fixed = p ~ 1, random = ~str(~vm(animal, ped.ainv) + vm(mother,  : \n  Variance structure is not positive definite\n"]<-NA
	matrix(as.numeric(x),nrow=nrow(x),ncol=ncol(x))
})
model5_hs <- lapply(model5_hs, function(x){
  x[x=="Error in asreml(fixed = p ~ 1, random = ~str(~vm(animal, ped.ainv) + vm(mother,  : \n  Variance structure is not positive definite\n"]<-NA
	matrix(as.numeric(x),nrow=nrow(x),ncol=ncol(x))
})
## some the covariance models are systematically not working well for only some scenarios ?!
table(c(sapply(model5_fs, function(x) which(apply(x,1,function(i) all(is.na(i))))),recursive=TRUE))
table(c(sapply(model5_hs, function(x) which(apply(x,1,function(i) all(is.na(i))))),recursive=TRUE))

## summaries for all the models
for(i in 1:8){
	for(j in c("fs","hs")){
		dat <- get(paste0("model",i,"_",j))
		assign(paste0("m",i,"_",j),apply(list2array(dat), c(1,2), mean, na.rm=TRUE))
	}
}

for(i in 1:2){
	for(j in c("fhs","fhs10")){
		dat <- get(paste0("model",i,"_",j))
		assign(paste0("m",i,"_",j),apply(list2array(dat), c(1,2), mean, na.rm=TRUE))
	}
}

load(paste0(data_wd,"real_sims.Rdata"))

for(i in 1:2){
	for(j in c("bt_socialD","bt_social","bt","rd","ss")){
		dat <- get(paste0("model",i,"_",j))
		assign(paste0("m",i,"_",j),apply(list2array(dat), c(1,2), mean, na.rm=TRUE))
	}
}

{
par(mfrow=c(2,1), mar=c(0,4,4,0))


	nr <-11
	nc <- ncol(scenarios2)

	plot(NA, ylim=c(0,0.9), xlim=c(1,nr+1), xaxt="n", ylab="Variance estimate")
	axis(3,1:nr+0.5,LETTERS[1:nr], cex.axis=2)
	# axis(3,)
	legend("topleft",c("full-sib","half=sib"),pch=c(17,19))
	
	for(i in 1:nc){
			arrows((1:nr)+i/nc,scenarios2[1:nr,i],(1:nr)+(i-1)/nc,scenarios2[1:nr,i], col=i, code=0)
			points((1:nr)+(i-0.3)/nc,m1_fs[1:nr,i], col=i, pch=17)
			points((1:nr)+(i-0.7)/nc,m1_hs[1:nr,i], col=i, pch=19)
	}
	abline(v=2:nr, col="grey")

	# nr <-6
	# nc <- ncol(scenarios2)
par( mar=c(4,4,0,0))
	plot(NA, ylim=c(0,0.9), xlim=c(1,nr+1), xaxt="n", ylab="Variance estimate")
	axis(1,rep(seq(0,1,length.out=7)[2:6],6) + rep(1:6,each=5),rep(c("A","Me","Mg","cov_AMg","E"),6))
	# axis(3,)

	
	for(i in 1:nc){
			arrows((1:nr)+i/nc,scenarios2[1:nr,i],(1:nr)+(i-1)/nc,scenarios2[1:nr,i], col=i, code=0)
			points((1:nr)+(i-0.3)/nc,m2_fs[1:nr,i], col=i, pch=17)
			points((1:nr)+(i-0.7)/nc,m2_hs[1:nr,i], col=i, pch=19)
	}
	abline(v=2:nr, col="grey")
}


### stacked barplots real

{

cols <- inferno(6)[2:6]
order <- c(4,1,3,2,5)
	all_mod<-lapply(c(1,2,4,3,5,6,9,7,8,10,11), function(i){
		rbind(
			scenarios2[i,order],
			m1_fs[i,order],m1_hs[i,order],m1_bt_socialD[i,order],m1_bt_social[i,order],m1_bt[i,order],m1_rd[i,order],m1_ss[i,order],
			m2_fs[i,order],m2_hs[i,order],m2_bt_socialD[i,order],m2_bt_social[i,order],m2_bt[i,order],m2_rd[i,order],m2_ss[i,order]
		)	
	})

	layout(matrix(c(1,2,3,12,4:11),byrow=TRUE,ncol=4))
	par(mar=c(3,4,1,1))
	for(i in 1:11){
		bp<-barplot(t(change2zero(all_mod[[i]])),space=c(0,1,0,0.3,0,0,0,0,0.5,0,0.3,0,0,0,0), col=cols, ylim=c(-0.1,1.2))
		axis(1,bp[c(1,1:8*2)] +c(0,rep(0.5,8)),c("sim",1:8))
	}
	

	plot(NULL,xaxt="n",yaxt="n",ylim=c(-1,1),xlim=c(-1,1), bty="n", ylab="", xlab="")
	legend("left", comp_names[order], pch=19, col=cols, bty="n")
	legend("right", model_names, pch=as.character(1:8), bty="n")
}




{

cols <- inferno(6)[2:6]
order <- c(4,1,3,2,5)
	all_mod<-lapply(c(1,2,4,3,5,6,9,7,8,10,11), function(i){
		rbind(
			scenarios2[i,order],
			m1_fs[i,order],m1_hs[i,order],m1_fhs[i,order],m1_fhs10[i,order],
			m2_fs[i,order],m2_hs[i,order],m2_fhs[i,order],m2_fhs10[i,order]
		)	
	})

	layout(matrix(c(1,2,3,12,4:11),byrow=TRUE,ncol=4))
	par(mar=c(3,4,1,1))
	for(i in 1:11){
		bp<-barplot(t(change2zero(all_mod[[i]])),space=c(0,1,0,0,0,0.3,0,0,0), col=cols, ylim=c(-0.1,1.2))
		axis(1,bp,c("sim",rep(c())))
	}
	

	plot(NULL,xaxt="n",yaxt="n",ylim=c(-1,1),xlim=c(-1,1), bty="n", ylab="", xlab="")
	legend("left", comp_names[order], pch=19, col=cols, bty="n")
	legend("right", model_names, pch=as.character(1:8), bty="n")
}





{

cols <- inferno(6)[2:6]
order <- c(4,1,3,2,5)
	all_mod<-lapply(c(1,2,4,3,5,6,9,7,8,10,11), function(i){
		rbind(
			scenarios2[i,order],
			m1_fs[i,order],m1_hs[i,order],
			m2_fs[i,order],m2_hs[i,order]
		)	
	})

	layout(matrix(c(1,2,3,12,4:11),byrow=TRUE,ncol=4))
	par(mar=c(3,4,1,1))
	for(i in 1:11){
		bp<-barplot(t(change2zero(all_mod[[i]])),space=c(0,1,0,0.3,0), col=cols, ylim=c(-0.1,1.2))
		axis(1,bp,c("sim",rep(c("fs","hs"),2)))
	}
	

	plot(NULL,xaxt="n",yaxt="n",ylim=c(-1,1),xlim=c(-1,1), bty="n", ylab="", xlab="")
	legend("left", comp_names[order], pch=19, col=cols, bty="n")
	legend("right", model_names, pch=as.character(1:8), bty="n")
}







### stacked barplots

{

cols <- inferno(6)[2:6]
order <- c(4,1,3,2,5)
	all_mod<-lapply(1:11, function(i){
		rbind(
			scenarios2[i,order],
			m1_fs[i,order],m1_hs[i,order],
			m2_fs[i,order],m2_hs[i,order],
			m3_fs[i,order],m3_hs[i,order],
			m4_fs[i,order],m4_hs[i,order],
			m5_fs[i,order],m5_hs[i,order],
			m6_fs[i,order],m6_hs[i,order],
			m7_fs[i,order],m7_hs[i,order],
			m8_fs[i,order],m8_hs[i,order]
		)	
	})

	par(mfrow=c(3,4),mar=c(3,4,1,1))
	for(i in 1:11){
		bp<-barplot(t(change2zero(all_mod[[i]])),space=c(0,1,0,rep(c(0.3,0),(nrow(all_mod[[i]])-3)/2)), col=cols, ylim=c(-0.1,1.2))
		axis(1,bp[c(1,1:8*2)] +c(0,rep(0.5,8)),c("sim",1:8))
	}
	

	plot(NULL,xaxt="n",yaxt="n",ylim=c(-1,1),xlim=c(-1,1), bty="n", ylab="", xlab="")
	legend("left", comp_names[order], pch=19, col=cols, bty="n")
	legend("right", model_names, pch=as.character(1:8), bty="n")
}






### total Va
{
simulated_tVa <- total_va(scenarios2,1:11,NULL)

fs_tVa <- NULL
hs_tVa <- NULL

for(i in 1:8){
		fs_tVa <- cbind(fs_tVa,total_va(get(paste0("m",i,"_fs")),1:11,NULL))
		hs_tVa <- cbind(hs_tVa,total_va(get(paste0("m",i,"_hs")),1:11,NULL))
}

par(mfrow=c(1,1))
plot(NA, xlim=c(1,12), ylim=c(0,0.6))
abline(v=1:11, col="grey")
nc<-ncol(fs_tVa)
arrows(1:11,simulated_tVa, 2:12,simulated_tVa, code=0)
for(i in 1:nc){
	points(fs_tVa[,i]~I(1:11 +i/(nc+1)), pch=19, col=i)
	points(hs_tVa[,i]~I(1:11 +i/(nc+1)), pch=17, col=i)	
}

}

{
nr<-ncol(fs_tVa)
nc<-nrow(fs_tVa)
par(mfrow=c(1,1))
plot(NA, xlim=c(1,nr+1), ylim=c(0,0.6))
abline(v=1:nr, col="grey")
arrows(rep(1:nr,each=nc) + rep((1:nc+0.5)/(nc+1),nr),rep(simulated_tVa,nr), rep(1:nr,each=nc) + rep((1:nc-0.5)/(nc+1),nr),rep(simulated_tVa,nr), code=0)

for(i in 1:nc){
	points(fs_tVa[i,]~I(1:nr +i/(nc+1)), pch=19, col=i)
	points(hs_tVa[i,]~I(1:nr +i/(nc+1)), pch=17, col=i)	
}

text(c(1:8+0.5),rep(0.6,8),model_names)
}
# model with just Va overestimate Va, unless there is no maternal (environment or genetic effects)

## although Va is often biased total Va well estimated with maternal ID, in pedigrees with half siblings

## maternal genetic but no maternal environment, typically overestimates total Va

## if there is maternal environment effects, but no maternal genetic and you model genetic but not maternal E then you estimate maternal genetic
# - i.e. unmodelled Me will appear as Mg - interestingly the total variance in model is also overestimated


bias_fs <- apply(fs_tVa,2,function(x)x-simulated_tVa)
bias_hs <- apply(hs_tVa,2,function(x)x-simulated_tVa)
{
nr<-ncol(bias_fs)
nc<-nrow(bias_fs)
par(mfrow=c(1,1), mar=c(1,4,1,1))
plot(NA, xlim=c(1,nr+1), ylim=c(-0.25,0.45), xaxt="n", ylab="Bias in total Va")
# axis(3,1:8 + 0.5, LETTERS[1:8])
abline(v=1:nr, col="grey")
abline(h=0, col=1)

for(i in 1:nc){
	points(bias_fs[i,]~I(1:nr +i/(nc+1)), pch=19, col=i)#i)
	points(bias_hs[i,]~I(1:nr +i/(nc+1)), pch=17, col=i)	
}


text(c(1:8+0.5),rep(0.4,8),sub("cov","\ncov", model_names))
}




mgS <- c(1,1,2,1,rep(2,7))

{
nr<-ncol(bias_fs)
nc<-nrow(bias_fs)
par(mfrow=c(1,1), mar=c(1,4,1,1))
plot(NA, xlim=c(1,nr+1), ylim=c(-0.25,0.45), xaxt="n", ylab="Bias in total Va")
# axis(3,1:8 + 0.5, LETTERS[1:8])
abline(v=1:nr, col="grey")
abline(h=0, col=1)

for(i in 1:nc){
	points(bias_fs[i,]~I(1:nr +i/(nc+1)), pch=19, col=mgS[i])#i)
	points(bias_hs[i,]~I(1:nr +i/(nc+1)), pch=17, col=mgS[i])	
}


text(c(1:8+0.5),rep(0.4,8),sub("cov","\ncov", model_names))
}






plot(bias_hs[,2],pch=17, col=scales::alpha(mgS,0.7), ylim=c(-0.2,0.05))
points(bias_fs[,2],pch=19, col=mgS)
abline(h=0, col=1)
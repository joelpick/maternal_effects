


rm(list=ls())

extract<-FALSE

main_dir <- "/Users/joelpick/Dropbox/0_fitness/fitness_timing/Data/Raw/Bonnet/"
data_dir <- "/Users/joelpick/Dropbox/0_fitness/fitness_timing/Data/Intermediate/"


### when using A - what happens to unknown fathers - what is relatedness between siblings?

#####
#--- functions
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


#####
# - Extract data from Ainv in Bonnet et al supp mat
####
if(extract){
	## run export R_MAX_VSIZE=32000000000 in terminal before running this 
	pops <- dir(main_dir)

	# 5:7
	for(i in pops[c(1:5,8:19)]){
		pop_dir <- paste0(main_dir,i)
		load(paste0(pop_dir,"/Ainv_",i))
		A <- Matrix::solve(ainv)
		colnames(A) <- rownames(A) <- rownames(ainv)
		rm(ainv)
		dat<-read.csv(paste0(pop_dir,"/data_",i,".csv"))
		mr<-maternal_relatedness_A(A,dat)
		save(mr, file=paste0(data_dir,"/mr_",i,".Rdat"))
		rm(A,dat,mr)
		cat(i," ")
	}

## extract data from a couple of additional pedigrees
	ped_BT<-read.csv("/Users/joelpick/Dropbox/0_fitness/fitness_timing/Data/Raw/ped_BT.csv")
	mr<-maternal_relatedness_ped(ped_BT)
	save(mr,file="/Users/joelpick/Dropbox/0_fitness/fitness_timing/Data/Intermediate/mr_btD.Rdat")

	ped_HS<-read.table("/Users/joelpick/Dropbox/0_fitness/fitness_timing/Data/Raw/ped_HS.txt", header=TRUE)
	mr<-maternal_relatedness_ped(ped_HS)
	save(mr,file="/Users/joelpick/Dropbox/0_fitness/fitness_timing/Data/Intermediate/mr_hsL.Rdat")

}


par(mfrow=c(4,5),mar=c(3,3,3,0))
	for(i in dir("/Users/joelpick/Dropbox/0_fitness/fitness_timing/Data/Raw/Bonnet/")){
		dat<-read.csv(paste0(main_dir,i,"/data_",i,".csv"))
		hist(dat$inbreeding, xlim=c(0,0.35), main=i)
	}

#####
### plotting 
####
breaks<-seq(0.25,0.8,length.out=56)


meta<-read.csv("/Users/joelpick/Dropbox/0_fitness/fitness_timing/Data/Raw/pops.csv")

par(mfrow=c(4,5),mar=c(3,3,3,0))
for(i in 1:19){
	load(paste0(data_dir,"/mr_",meta$code[i],".Rdat"))
	hist(mr, main=meta$population[i], breaks=breaks, xlab="Average relatedness within mothers", xlim=c(0.25,0.8), col=ifelse(meta$taxa[i]=="bird","yellow","green"))
	abline(v=mean(mr), col="red", lwd=3)
	abline(v=c(0.25,0.5), col="blue", lwd=1)
}

par(mfrow=c(2,3),mar=c(3,3,3,0))
for(i in 1:6){
	load(paste0(data_dir,"/mr_",meta$code[i],".Rdat"))
	hist(mr, main=meta$population[i], breaks=breaks, xlab="Average relatedness within mothers", xlim=c(0.25,0.8), col=ifelse(meta$taxa[i]=="bird","yellow","green"))
	abline(v=mean(mr), col="red", lwd=3)
	abline(v=c(0.25,0.5), col="blue", lwd=1)
}

par(mfrow=c(2,2),mar=c(3,3,3,0))
for(i in 7:10){
	load(paste0(data_dir,"/mr_",meta$code[i],".Rdat"))
	hist(mr, main=meta$population[i], breaks=breaks, xlab="Average relatedness within mothers", xlim=c(0.25,0.8), col=ifelse(meta$taxa[i]=="bird","yellow","green"))
	abline(v=mean(mr), col="red", lwd=3)
	abline(v=c(0.25,0.5), col="blue", lwd=1)
}

par(mfrow=c(3,3),mar=c(3,3,3,0))
for(i in 11:19){
	load(paste0(data_dir,"/mr_",meta$code[i],".Rdat"))
	hist(mr, main=meta$population[i], breaks=breaks, xlab="Average relatedness within mothers", xlim=c(0.25,0.8), col=ifelse(meta$taxa[i]=="bird","yellow","green"))
	abline(v=mean(mr), col="red", lwd=3)
	abline(v=c(0.25,0.5), col="blue", lwd=1)
}


sapply(1:19, function(i){
	load(paste0(data_dir,"/mr_",meta$code[i],".Rdat"))
	mean(mr)
})



## compare maternal relatedness pedigree to ainv 
ped_RD<-read.csv("/Users/joelpick/Dropbox/0_fitness/fitness_timing/Data/Raw/ped_RD.csv")
ped_SSH<-read.csv("/Users/joelpick/Dropbox/0_fitness/fitness_timing/Data/Raw/ped_SSH.csv")
ped_SFW<-read.csv("/Users/joelpick/Dropbox/0_fitness/fitness_timing/Data/Raw/ped_SFW.csv")
ped_SV<-read.csv("/Users/joelpick/Dropbox/0_fitness/fitness_timing/Data/Raw/ped_SV.csv")


{
par(mfrow=c(2,4),mar=c(3,3,1,0))
hist(maternal_relatedness_ped(ped_RD), main="Red Deer", breaks=breaks, xlab="Average relatedness within mothers", xlim=c(0.25,0.8)); abline(v=mean(maternal_relatedness_ped(ped_RD)), col="red", lwd=3)
hist(maternal_relatedness_ped(ped_SSH), main="Soay Sheep", breaks=breaks, xlab="Average relatedness within mothers", xlim=c(0.25,0.8)); abline(v=mean(maternal_relatedness_ped(ped_SSH)), col="red", lwd=3)
hist(maternal_relatedness_ped(ped_SV), main="Snow Voles", breaks=breaks, xlab="Average relatedness within mothers", xlim=c(0.25,0.8)); abline(v=mean(maternal_relatedness_ped(ped_SV)), col="red", lwd=3)
hist(maternal_relatedness_ped(ped_SFW), main="Superb Fairy Wren", breaks=breaks, xlab="Average relatedness within mothers", xlim=c(0.25,0.8)); abline(v=mean(maternal_relatedness_ped(ped_SFW)), col="red", lwd=3)

for(i in c(c("rdR","ssS","svG","sfC"))){
	load(paste0(data_dir,"/mr_",i,".Rdat"))
	hist(mr, main="", breaks=breaks, xlab="Average relatedness within mothers", xlim=c(0.25,0.8))
	abline(v=mean(mr), col="red", lwd=3)
	abline(v=c(0.25,0.5), col="blue", lwd=1)
}
}

par(mfrow=c(3,1),mar=c(3,3,1,0))

hist(maternal_relatedness_ped(ped_SFW), main="Superb Fairy Wren", breaks=breaks, xlab="Average relatedness within mothers", xlim=c(0.25,0.8)); abline(v=mean(maternal_relatedness_ped(ped_SFW)), col="red", lwd=3)

load(paste0(data_dir,"/mr_sfC.Rdat"))
hist(mr, main="", breaks=breaks, xlab="Average relatedness within mothers", xlim=c(0.25,0.8))

SV_A<-nadiv::makeA(ped_SV)
mr_SV<-maternal_relatedness_A(SV_A,ped_SV)
hist(mr_SV, breaks=breaks, xlab="Average relatedness within mothers", xlim=c(0.25,0.8)); abline(v=mean(mr_SV), col="red", lwd=3)



# mean(maternal_relatedness_ped(ped_RD))
# mean(maternal_relatedness_ped(ped_SSH))
# mean(maternal_relatedness_ped(ped_SV))
# mean(maternal_relatedness_ped(ped_SFW))
# mean(maternal_relatedness_ped(ped_BT))
# mean(maternal_relatedness_ped(ped_HS))

# par(mfrow=c(3,2))
# hist(maternal_relatedness_ped(ped_RD), main="Red Deer", breaks=breaks, xlab="Average relatedness within mothers"); abline(v=mean(maternal_relatedness_ped(ped_RD)), col="red", lwd=3)
# hist(maternal_relatedness_ped(ped_SSH), main="Soay Sheep", breaks=breaks, xlab="Average relatedness within mothers"); abline(v=mean(maternal_relatedness_ped(ped_SSH)), col="red", lwd=3)
# hist(maternal_relatedness_ped(ped_SV), main="Snow Voles", breaks=breaks, xlab="Average relatedness within mothers"); abline(v=mean(maternal_relatedness_ped(ped_SV)), col="red", lwd=3)
# hist(maternal_relatedness_ped(ped_SFW), main="Superb Fairy Wren", breaks=breaks, xlab="Average relatedness within mothers"); abline(v=mean(maternal_relatedness_ped(ped_SFW)), col="red", lwd=3)
# hist(maternal_relatedness_ped(ped_BT), main="Blue Tits", breaks=breaks, xlab="Average relatedness within mothers"); abline(v=mean(maternal_relatedness_ped(ped_BT)), col="red", lwd=3)
# hist(maternal_relatedness_ped(ped_HS), main="House Sparrows", breaks=breaks, xlab="Average relatedness within mothers"); abline(v=mean(maternal_relatedness_ped(ped_HS)), col="red", lwd=3)


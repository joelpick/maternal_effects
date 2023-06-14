
rm(list=ls())

extract<-FALSE
wd <- "/Users/joelpick/github/maternal_effects/"
main_dir <- paste0(wd,"Data/Raw/Bonnet/")
data_dir <- paste0(wd,"Data/Intermediate/")

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

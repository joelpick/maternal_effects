## make this faster- all in one function, but could take the cousin/au specific parts out and then make into a function - would also cut down on repetitive dam and sire code. Would mean that splitting was happening much less.

#################################################
rbind_notAnnoying <- function(..., names=NULL){
  x <- list(...)
  y <- lapply(x, function(y){
    names(y) <- if(is.null(names)) names(x[[1]]) else names
    return(y)  
  } )
  do.call(rbind,y)
}

## function to make products of all combinations and sum them 
combo_prod <- function(x) if(length(x)>1){sum( combn(x, m =2)[1, ] * combn(x, m =2)[2, ])}else{0}


combo_au <- function(x) if(length(x)>1){sum(x * (length(x)-1))}else{0}


n_cousin <- function(p,gp,ped){
  f1 <- formula(paste("animal~",p,"+",gp))
  f2 <- formula(paste("animal~",gp))
  d1 <- aggregate(f1, ped,length)
	sum(aggregate(f2,d1,combo_prod)$animal)
}

n_au <- function(link,ped){
	
	p <- if(substr(link,1,1)=="M"){"dam"}else{"sire"}
	gp <- if(substr(link,3,3)=="M"){"dam"}else if(substr(link,3,3)=="F"){"sire"} else{"pair"}
  f1 <- formula(paste("animal~",p,"+",link))
  f2 <- formula(paste("animal~",gp))
  d1 <- aggregate(f1, ped,length)
  d2 <- aggregate(f2, ped,length)
	sum(d1$animal * (d2[match(d1[,link],d2[,gp]),"animal"]-1), na.rm=TRUE)
}
total_links <- function(ped) nrow(ped) * (nrow(ped) - 1) / 2

non_zero_links <- function(ped){
	pedA<-nadiv::makeA(ped[,1:3])
	sum(pedA[lower.tri(pedA)]>0)
}





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




cousins <- function(ped, phenotyped=NULL){	
	colnames(ped) <- c("animal","dam","sire")

	# dummy_dam <- which(is.na(ped$dam) & !is.na(ped$sire))
	# dummy_sire <- which(is.na(ped$sire) & !is.na(ped$dam))

	# ped[dummy_sire,"sire"] <- paste0("ds_",seq_along(dummy_sire))
	# ped[dummy_dam,"dam"] <- paste0("dd_",seq_along(dummy_dam))

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

	ped2 <- subset(ped,animal %in% phenotyped)	

	## stack maternal and paternal grandparents, to get allows for all cousins relationships, not just through mothers or father
	GM <- rbind_notAnnoying(ped2[,c("animal","dam","MGM","MGF","MGP")], ped2[,c("animal","sire","PGM","PGF","PGP")])

	cousin_FS_FF <- n_cousin(p="dam", gp="MGP", ped=ped2)
	cousin_MHS_FF <- n_cousin(p="dam", gp="MGM", ped=ped2) - cousin_FS_FF
	cousin_PHS_FF <- n_cousin(p="dam", gp="MGF", ped=ped2) - cousin_FS_FF

	cousin_FS_MM <- n_cousin(p="sire", gp="PGP", ped=ped2)
	cousin_MHS_MM <- n_cousin(p="sire", gp="PGM", ped=ped2) - cousin_FS_MM
	cousin_PHS_MM <- n_cousin(p="sire", gp="PGF", ped=ped2) - cousin_FS_MM

	cousin_FS <- n_cousin(p="dam", gp="MGP", ped=GM)
	cousin_FS_FM <- cousin_FS-cousin_FS_FF-cousin_FS_MM

	cousin_MS <- n_cousin(p="dam", gp="MGM", ped=GM)
	cousin_MHS_FM <- cousin_MS-cousin_FS-cousin_MHS_FF-cousin_MHS_MM

	cousin_PS <- n_cousin(p="dam", gp="MGF", ped=GM)
	cousin_PHS_FM <- cousin_PS-cousin_FS-cousin_PHS_FF-cousin_PHS_MM


	c(cousin_FS_FF=cousin_FS_FF,
		cousin_FS_FM=cousin_FS_FM,
		cousin_FS_MM=cousin_FS_MM,
		cousin_MHS_FF=cousin_MHS_FF,
		cousin_MHS_FM=cousin_MHS_FM,
		cousin_MHS_MM=cousin_MHS_MM,
		cousin_PHS_FF=cousin_PHS_FF,
		cousin_PHS_FM=cousin_PHS_FM,
		cousin_PHS_MM=cousin_PHS_MM)
}


au <- function(ped, phenotyped=NULL){	
	colnames(ped) <- c("animal","dam","sire")

 	# dummy_dam <- which(is.na(ped$dam) & !is.na(ped$sire))
	# dummy_sire <- which(is.na(ped$sire) & !is.na(ped$dam))

	# ped[dummy_sire,"sire"] <- paste0("ds_",seq_along(dummy_sire))
	# ped[dummy_dam,"dam"] <- paste0("dd_",seq_along(dummy_dam))


	ped$MGM <- ped[match(ped[,2],ped[,1]),2]
	ped$MGF <- ped[match(ped[,2],ped[,1]),3]
	ped$PGM <- ped[match(ped[,3],ped[,1]),2]
	ped$PGF <- ped[match(ped[,3],ped[,1]),3]

	ped$MGP <- paste(ped$MGM,ped$MGF)
	ped$MGP<-ifelse(ped$MGP== "NA NA", NA, ped$MGP)
	ped$PGP <- paste(ped$PGM,ped$PGF)
	ped$PGP<-ifelse(ped$PGP== "NA NA", NA, ped$PGP)
# ped$MGP<-ifelse(grepl("NA",MGP), NA, ped$MGP)
# ped$PGP<-ifelse(grepl("NA",PGP), NA, ped$PGP)

	ped$pair<-paste(ped$dam,ped$sire)
	ped$pair<-ifelse(ped$pair== "NA NA", NA, ped$pair)
	

	## if dont specify phenotyped, then assume all phenotyped
	if(is.null(phenotyped)) phenotyped <- ped[,1]

	ped2 <- subset(ped,animal %in% phenotyped)	


	GM <- rbind_notAnnoying(ped2[,c("animal","dam","MGM","MGF","MGP")], ped2[,c("animal","sire","PGM","PGF","PGP")])

	# related through a maternal grandmother (au_FS_F + au_MHS_F ?)
	au_MS_F <- n_au("MGM",ped2)

	# related through a maternal grandfather (au_FS_F + au_PHS_F ?)
	au_PS_F <- n_au("MGF",ped2)

	# related through a maternal grandparents (au_FS_F ?)
	au_FS_F <- n_au("MGP",ped2)

	au_PHS_F <- au_PS_F - au_FS_F
	au_MHS_F <- au_MS_F - au_FS_F

	# related through a paternal grandmother (au_FS_M + au_MHS_M ?)
	au_MS_M <- n_au("PGM",ped2)

	# related through a paternal grandfather (au_FS_M + au_PHS_F ?)
	au_PS_M <- n_au("PGF",ped2)

	# related through a paternal grandparents (au_FS_M ?)
	au_FS_M <- n_au("PGP",ped2)

	au_PHS_M <- au_PS_M - au_FS_M
	au_MHS_M <- au_MS_M - au_FS_M

	c(au_FS_F=au_FS_F,au_FS_M=au_FS_M,au_MHS_F=au_MHS_F,au_MHS_M=au_MHS_M,au_PHS_F=au_PHS_F,au_PHS_M=au_PHS_M)


}



gp <- function(ped, phenotyped=NULL){
	#ped=ped_bt
	colnames(ped) <- c("animal","dam","sire")

## more efficient to just work out everyones maternal and paternal grand parents and then sum the non-NAs?

	if(is.null(phenotyped)) phenotyped <- ped[,1]

	ped$MGM <- ped[match(ped[,2],ped[,1]),2]
	ped$MGF <- ped[match(ped[,2],ped[,1]),3]
	ped$PGM <- ped[match(ped[,3],ped[,1]),2]
	ped$PGF <- ped[match(ped[,3],ped[,1]),3]

	apply(ped[ped[,1] %in% phenotyped,c("dam","sire","MGM","MGF","PGM","PGF")],2,function(x) sum(x%in% phenotyped) )
}




sibs <- function(ped, phenotyped=NULL){
	colnames(ped) <- c("animal","dam","sire")

	## if dont specify phenotyped, then assume all phenotyped
	if(is.null(phenotyped)) phenotyped <- ped[,1]
	
	ped2 <- subset(ped,animal %in% phenotyped)	

	mat_sib <- sum(table(ped2$dam) * (table(ped2$dam)-1)/2)
	pat_sib <- sum(table(ped2$sire) * (table(ped2$sire)-1)/2)
	ped2$pair <- paste(ped2$dam,ped2$sire)
	n_pair <-table(ped2$pair[!grepl("NA",ped2$pair)])
	FS <- sum(n_pair * (n_pair-1)/2)
	c(FS=FS, MHS=mat_sib-FS, PHS=pat_sib-FS)

}






# c("dam","MGF","MGM","au_FS_F","cousin_FS_FF","cousin_MHS_FF")


# # ped_stat(fs_peds[[1]], phenotyped=fs_peds[[1]][!is.na(fs_peds[[1]][,2]),1])
# dd<-ped_stat(ped)
# sum(dd[c("dam","MGF","MGM","au_FS_F","cousin_FS_FF","cousin_MHS_FF")])
# sum(dd[c("FS","MHS","PHS")])

mat_links <- function(ped, phenotyped=NULL){
	# ped_stat(ped, phenotyped=ped[!is.na(ped[,2]),1])
	dd<-ped_stat(ped)
	c(mat_sibs=sum(dd[c("FS","MHS")]), mat_links=sum(dd[c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")]), other=sum(dd[!names(dd)%in%c( "individuals" ,"links" ,"FS","MHS","dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")]) )
	# sum(dd[c("FS","MHS","PHS")])

}


pedSum <- function(ped) pedantics::pedStatSummary(pedantics::pedigreeStats(ped, includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE))[2:12]

prop_stat <- function(ped) c(pedantics::pedStatSummary(pedantics::pedigreeStats(ped[,1:3], includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE))[2:12],cousins(ped))/non_zero_links(ped)


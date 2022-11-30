## make this faster- all in one function, but could take the cousin/au specific parts out and then make into a function - would also cut down on repetitive dam and sire code. Would mean that splitting was happening much less.




cousins2 <- function(ped, phenotyped=NULL){	
	colnames(ped) <- c("animal","dam","sire")

	# dummy_dam <- which(is.na(ped$dam) & !is.na(ped$sire))
	# dummy_sire <- which(is.na(ped$sire) & !is.na(ped$dam))

	# ped[dummy_sire,"sire"] <- paste0("ds_",seq_along(dummy_sire))
	# ped[dummy_dam,"dam"] <- paste0("dd_",seq_along(dummy_dam))


	## if dont specify phenotyped, then assume all phenotyped
	if(is.null(phenotyped)) phenotyped <- ped[,1]

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

	ped2 <- subset(ped,animal %in% phenotyped)	

	## stack maternal and paternal grandparents, to get allows for all cousins relationships, not just through mothers or father
	GM <- rbind_notAnnoying(ped2[,c("animal","dam","MGM","MGF","MGP")], ped2[,c("animal","sire","PGM","PGF","PGP")])

	## group by maternal grandparents
	# mat <- aggregate(animal~dam+MGM+MGF+MGP, ped2,length)
	# cousin_FS_FF <- sum(aggregate(animal~MGP,mat,combo_prod)$animal)
	# cousin_MS_FF <- sum(aggregate(animal~MGM,mat,combo_prod)$animal)
	# cousin_PS_FF <- sum(aggregate(animal~MGF,mat,combo_prod)$animal)
  # cousin_MHS_FF <- cousin_MS_FF - cousin_FS_FF
  # cousin_PHS_FF <- cousin_PS_FF - cousin_FS_FF

	# ## group by paternal grandparents
  # pat <- aggregate(animal~sire+PGM+PGF+PGP, ped2,length)
	# cousin_FS_MM <- sum(aggregate(animal~PGP,pat,combo_prod)$animal)
	# cousin_MS_MM <- sum(aggregate(animal~PGM,pat,combo_prod)$animal)
	# cousin_PS_MM <- sum(aggregate(animal~PGF,pat,combo_prod)$animal)
  # cousin_MHS_MM <- cousin_MS_MM - cousin_FS_MM
  # cousin_PHS_MM <- cousin_PS_MM - cousin_FS_MM

  # ## group by all grandparents
  # mp <- aggregate(animal~dam+MGM+MGF+MGP, GM,length)
	# cousin_FS <- sum(aggregate(animal~MGP,mp,combo_prod)$animal)
	# cousin_MS <- sum(aggregate(animal~MGM,mp,combo_prod)$animal)
	# cousin_PS <- sum(aggregate(animal~MGF,mp,combo_prod)$animal)

	# cousin_FS_FM <- cousin_FS-cousin_FS_FF-cousin_FS_MM
	# cousin_MHS_FM <- cousin_MS-cousin_FS-cousin_MHS_FF-cousin_MHS_MM
	# cousin_PHS_FM <- cousin_PS-cousin_FS-cousin_PHS_FF-cousin_PHS_MM

	cousin_FS_FF <- n_cousin(p="dam", gp="MGP", ped=ped)
	cousin_MHS_FF <- n_cousin(p="dam", gp="MGM", ped=ped) - cousin_FS_FF
	cousin_PHS_FF <- n_cousin(p="dam", gp="MGF", ped=ped) - cousin_FS_FF

	cousin_FS_MM <- n_cousin(p="sire", gp="PGP", ped=ped)
	cousin_MHS_MM <- n_cousin(p="sire", gp="PGM", ped=ped) - cousin_FS_MM
	cousin_PHS_MM <- n_cousin(p="sire", gp="PGF", ped=ped) - cousin_FS_MM


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


au2 <- function(ped){	

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


	GM <- rbind_notAnnoying(ped[,c("animal","dam","MGM","MGF","MGP")], ped[,c("animal","sire","PGM","PGF","PGP")])

	# related through a maternal grandmother (au_FS_F + au_MHS_F ?)
	au_MS_F <- n_au("MGM",ped)

	# related through a maternal grandfather (au_FS_F + au_PHS_F ?)
	au_PS_F <- n_au("MGF",ped)

	# related through a maternal grandparents (au_FS_F ?)
	au_FS_F <- n_au("MGP",ped)

	au_PHS_F <- au_PS_F - au_FS_F
	au_MHS_F <- au_MS_F - au_FS_F

	# related through a paternal grandmother (au_FS_M + au_MHS_M ?)
	au_MS_M <- n_au("PGM",ped)

	# related through a paternal grandfather (au_FS_M + au_PHS_F ?)
	au_PS_M <- n_au("PGF",ped)

	# related through a paternal grandparents (au_FS_M ?)
	d1 <- aggregate(animal~sire+PGP,ped,length)
	d2 <- aggregate(animal~pair,ped,length)
	sum(d1$animal * (d2[match(d1$PGP,d2$pair),"animal"]-1))
	au_FS_M <- n_au("PGP",ped)

	au_PHS_M <- au_PS_M - au_FS_M
	au_MHS_M <- au_MS_M - au_FS_M




	c(au_FS_F=au_FS_F,au_FS_M=au_FS_M,au_MHS_F=au_MHS_F,au_MHS_M=au_MHS_M,au_PHS_F=au_PHS_F,au_PHS_M=au_PHS_M)

# c("au_FS_F","au_FS_M","au_MHS_F","au_MHS_M","au_PHS_F","au_PHS_M")

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




# au(ped)
# cousins(ped)
# gp(ped) 

total_links <- function(ped) nrow(ped) * (nrow(ped) - 1) / 2

non_zero_links <- function(ped){
	pedA<-nadiv::makeA(ped[,1:3])
	sum(pedA[lower.tri(pedA)]>0)
}

ped_stat <- function(ped,phenotyped=NULL){
	if(is.null(phenotyped)) phenotyped <- ped[,1]
	stat <- c(
		individuals = length(phenotyped),
		links = length(phenotyped) * (length(phenotyped) - 1) / 2,
		# po(ped,phenotyped),
		sibs(ped,phenotyped),
		gp(ped,phenotyped),
		au(ped,phenotyped),
		cousins(ped,phenotyped)
	)

	all_rel <- c("individuals","links","dam","sire","FS","MHS","PHS","MGF","MGM","PGF","PGM","au_FS_F","au_FS_M","au_MHS_F","au_MHS_M","au_PHS_F","au_PHS_M","cousin_FS_FF","cousin_FS_FM","cousin_FS_MM","cousin_MHS_FF","cousin_MHS_FM","cousin_MHS_MM","cousin_PHS_FF","cousin_PHS_FM","cousin_PHS_MM")

	out<-rep(0,length(all_rel))
	names(out)<-all_rel

	out[names(stat)]<-stat
	out
}


# c("dam","MGF","MGM","au_FS_F","cousin_FS_FF","cousin_MHS_FF")


# # ped_stat(fs_peds[[1]], phenotyped=fs_peds[[1]][!is.na(fs_peds[[1]][,2]),1])
# dd<-ped_stat(ped)
# sum(dd[c("dam","MGF","MGM","au_FS_F","cousin_FS_FF","cousin_MHS_FF")])
# sum(dd[c("FS","MHS","PHS")])

mat_links <- function(ped, phenotyped=NULL){
	# ped_stat(ped, phenotyped=ped[!is.na(ped[,2]),1])
	dd<-ped_stat(ped)
	c(dd[c("FS","MHS","PHS")], mat_links=sum(dd[c("dam","MGF","MGM","au_FS_F","cousin_FS_FF","cousin_MHS_FF")]), other=sum(dd[!names(dd)%in%c("FS","MHS","PHS","dam","MGF","MGM","au_FS_F","cousin_FS_FF","cousin_MHS_FF")]) )
	# sum(dd[c("FS","MHS","PHS")])

}


pedSum <- function(ped) pedantics::pedStatSummary(pedantics::pedigreeStats(ped, includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE))[2:12]

prop_stat <- function(ped) c(pedantics::pedStatSummary(pedantics::pedigreeStats(ped[,1:3], includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE))[2:12],cousins(ped))/non_zero_links(ped)


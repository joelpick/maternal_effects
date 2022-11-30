## make this faster- all in one function, but could take the cousin/au specific parts out and then make into a function - would also cut down on repetitive dam and sire code. Would mean that splitting was happening much less.

cousins <- function(ped, phenotyped=NULL){
	
	## number sharing grandparents?



	#ped=ped_bt
	colnames(ped) <- c("animal","dam","sire")

	## if dont specify phenotyped, then assume all phenotyped
	if(is.null(phenotyped)) phenotyped <- ped[,1]

	d_s<-split(ped,ped$dam)
	s_s<-split(ped,ped$sire)

	# x<-d_s[[1]]
	out_d<-do.call(rbind,lapply(d_s, function(x){
			## work out which offspring are themselves parents
		sisters <- x$animal[x$animal %in% ped$dam]
		brothers <- x$animal[x$animal %in% ped$sire]
		sibs <- c(sisters,brothers)
		
		# all_sibs <- x$animal
		# all_sires <- x$sire

		## sires of those offspring that are parents
		sires <- 	x$sire[match(sibs,x$animal)]

		## sex of those that are parents
		sexes <- c(rep("F",length(sisters)),rep("M",length(brothers)))
		
		## number of *phenotyped* offspring they have
		# off <- c(sapply(sisters, function(i) nrow(d_s[[i]])),	sapply(brothers, function(i) nrow(s_s[[i]])), recursive=TRUE)
		off <- c(
			sapply(sisters, function(i) sum(d_s[[i]]$animal %in% phenotyped)),#nrow(d_s[[i]])),
			sapply(brothers, function(i) sum(s_s[[i]]$animal %in% phenotyped)),#nrow(s_s[[i]])),
			recursive=TRUE)

	## this will fail with full sib matings, because will double count

		if(length(sibs)>1){

			ex_fam <- do.call(rbind,lapply(1: (length(sibs)-1), function(i){
				## loop thorugh all the comibiations of sibling parents
				rbind(
					## number of cousin relationship will be the number of offspring of one sibling * n offspring of the other
					data.frame(
						relationship="cousin",
						n = off[i]*off[(i+1):length(sibs)],
						parent_r = ifelse(sires[(i+1):length(sibs)] %in% sires[i], "FS","MHS"),
						sex = paste0(sexes[i],sexes[(i+1):length(sibs)])
					)		
					# # Number of aunt/uncle relationships from siblings to focal individuals offspring
					#, data.frame(
					# 	relationship="a_u",
					# 	n = off[i],
					# 	parent_r = ifelse(sires[(i+1):length(sibs)] %in% sires[i], "FS","MHS"),
					# 	sex = paste0(sexes[(i+1):length(sibs)],sexes[i])
					# ),
					# ## Number of aunt/uncle relationships from focal individual to siblings offspring
					# data.frame(
					# 	relationship="a_u",
					# 	n = off[(i+1):length(sibs)],
					# 	parent_r = ifelse(sires[(i+1):length(sibs)] %in% sires[i], "FS","MHS"),
					# 	sex = paste0(sexes[i],sexes[(i+1):length(sibs)])
					# )
				)
			}))
	## have to do the same for paternal half sibs - because linked through dad doesn't matter exact relation?
	}else{
		NULL
	}


	}))

	out_s<-do.call(rbind,lapply(s_s, function(x){
			## work out which offspring are themselves parents
		sisters <- x$animal[x$animal %in% ped$dam]
		brothers <- x$animal[x$animal %in% ped$sire]
		sibs <- c(sisters,brothers)
		## aunt and uncles don't need to have offspring themselves


		## sires of those offspring that are parents
		dams <- 	x$dam[match(sibs,x$animal)]

		## sex of those that are parents
		sexes <- c(rep("F",length(sisters)),rep("M",length(brothers)))
		
		## number of *phenotyped* offspring they have
		# off <- c(sapply(sisters, function(i) nrow(d_s[[i]])),	sapply(brothers, function(i) nrow(s_s[[i]])), recursive=TRUE)
		off <- c(
			sapply(sisters, function(i) sum(d_s[[i]]$animal %in% phenotyped)),#nrow(d_s[[i]])),
			sapply(brothers, function(i) sum(s_s[[i]]$animal %in% phenotyped)),#nrow(s_s[[i]])),
			recursive=TRUE)
	## this will fail with full sib matings, because will double count

		if(length(sibs)>1){

			ex_fam <- do.call(rbind,lapply(1: (length(sibs)-1), function(i){
				## loop through all the combinations of sibling parents
				rbind(
					## number of cousin relationship will be the number of offspring of one sibling * n offspring of the other
					data.frame(
						relationship="cousin",
						n = off[i]*off[(i+1):length(sibs)],
						parent_r = ifelse(dams[(i+1):length(sibs)] %in% dams[i], "FS","PHS"),
						sex = paste0(sexes[i],sexes[(i+1):length(sibs)])
					)
					# ## Number of aunt/uncle relationships from siblings to focal individuals offspring
					#, data.frame(
					# 	relationship="a_u",
					# 	n = off[i],
					# 	parent_r = ifelse(dams[(i+1):length(sibs)] %in% dams[i], "FS","PHS"),
					# 	sex = paste0(sexes[(i+1):length(sibs)],sexes[i])
					# ),
					# ## Number of aunt/uncle relationships from focal individual to siblings offspring
					# data.frame(
					# 	relationship="a_u",
					# 	n = off[(i+1):length(sibs)],
					# 	parent_r = ifelse(dams[(i+1):length(sibs)] %in% dams[i], "FS","PHS"),
					# 	sex = paste0(sexes[i],sexes[(i+1):length(sibs)])
					# )
				)
			}))
	## have to do the same for paternal half sibs - because linked through dad doesn't matter exact relation?
	}else{
		NULL
	}


	}))

	out <- rbind(out_d,subset(out_s,parent_r!="FS"))
	out$sex <- ifelse(out$sex=="MF","FM",out$sex)
	c_sum <- aggregate(n~sex+parent_r+relationship,out,sum)
	c_sum2 <- c_sum[,4]
	names(c_sum2) <- paste(c_sum[,3],c_sum[,2],c_sum[,1],sep="_")
	c_sum2
}



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

n_cousin <- function(p,gp,ped){
  f1 <- formula(paste("animal~",p,"+",gp))
  f2 <- formula(paste("animal~",gp))
  d1 <- aggregate(f1, ped,length)
	sum(aggregate(f2,d1,combo_prod)$animal)
}

cousins2 <- function(ped, phenotyped){	
	colnames(ped) <- c("animal","dam","sire")

	## if dont specify phenotyped, then assume all phenotyped
	if(is.null(phenotyped)) phenotyped <- ped[,1]

  ## grandparents
	ped$MGM <- ped[match(ped[,2],ped[,1]),2]
	ped$MGF <- ped[match(ped[,2],ped[,1]),3]
	ped$PGM <- ped[match(ped[,3],ped[,1]),2]
	ped$PGF <- ped[match(ped[,3],ped[,1]),3]

	## grandparent pairs
	ped$MGP <- paste(ped$MGM,ped$MGF)
	ped$MGP<-ifelse(ped$MGP== "NA NA", NA, ped$MGP)
	ped$PGP <- paste(ped$PGM,ped$PGF)
	ped$PGP<-ifelse(ped$PGP== "NA NA", NA, ped$PGP)

	ped2 <- subset(ped,animal %in% phenotyped)	

	## stack maternal and paternal grandparents, to get allows for all cousins relationships, not just through mothers or father
	GM <- rbind_notAnnoying(ped2[,c("animal","dam","MGM","MGF","MGP")], ped2[,c("animal","sire","PGM","PGF","PGP")])

	## group by maternal grandparents
	mat <- aggregate(animal~dam+MGM+MGF+MGP, ped2,length)
	cousin_FS_FF <- sum(aggregate(animal~MGP,mat,combo_prod)$animal)
	cousin_MS_FF <- sum(aggregate(animal~MGM,mat,combo_prod)$animal)
	cousin_PS_FF <- sum(aggregate(animal~MGF,mat,combo_prod)$animal)
  cousin_MHS_FF <- cousin_MS_FF - cousin_FS_FF
  cousin_PHS_FF <- cousin_PS_FF - cousin_FS_FF

	## group by paternal grandparents
  pat <- aggregate(animal~sire+PGM+PGF+PGP, ped2,length)
	cousin_FS_MM <- sum(aggregate(animal~PGP,pat,combo_prod)$animal)
	cousin_MS_MM <- sum(aggregate(animal~PGM,pat,combo_prod)$animal)
	cousin_PS_MM <- sum(aggregate(animal~PGF,pat,combo_prod)$animal)
  cousin_MHS_MM <- cousin_MS_MM - cousin_FS_MM
  cousin_PHS_MM <- cousin_PS_MM - cousin_FS_MM

  ## group by all grandparents
  mp <- aggregate(animal~dam+MGM+MGF+MGP, GM,length)
	cousin_FS <- sum(aggregate(animal~MGP,mp,combo_prod)$animal)
	cousin_MS <- sum(aggregate(animal~MGM,mp,combo_prod)$animal)
	cousin_PS <- sum(aggregate(animal~MGF,mp,combo_prod)$animal)

	cousin_FS_FM <- cousin_FS-cousin_FS_FF-cousin_FS_MM
	cousin_MHS_FM <- cousin_MS-cousin_FS-cousin_MHS_FF-cousin_MHS_MM
	cousin_PHS_FM <- cousin_PS-cousin_FS-cousin_PHS_FF-cousin_PHS_MM

	c(cousin_FS_FF=cousin_FS_FF,cousin_FS_FM=cousin_FS_FM,cousin_FS_MM=cousin_FS_MM,cousin_MHS_FF=cousin_MHS_FF,cousin_MHS_FM=cousin_MHS_FM,cousin_MHS_MM=cousin_MHS_MM,cousin_PHS_FF=cousin_PHS_FF,cousin_PHS_FM=cousin_PHS_FM,cousin_PHS_MM=cousin_PHS_MM)
}



au <- function(ped, phenotyped=NULL){
	#ped=ped_bt
	colnames(ped) <- c("animal","dam","sire")
	
	## if dont specify phenotyped, then assume all phenotyped
	if(is.null(phenotyped)) phenotyped <- ped[,1]

	# ped2 <- subset(ped,animal %in% phenotyped)

	## all dam and sire	families
	d_s<-split(ped,ped$dam)
	s_s<-split(ped,ped$sire)

	# d_s_P<- d_s[names(d_s) %in% phenotyped] 
	# s_s_P<- s_s[names(s_s) %in% phenotyped] 

	# x<-d_s[[2]]
	out_d<-do.call(rbind,lapply(d_s, function(x){
			
			## work out which offspring are themselves parents
		sisters <- x$animal[x$animal %in% ped$dam]
		brothers <- x$animal[x$animal %in% ped$sire]
		sibs <- c(sisters,brothers)
		
		## phenotyped sibs
		# all_sibs <- x$animal
		# all_sires <- x$sire
		all_sibs <- x$animal[x$animal %in% phenotyped]
		all_sires <- x$sire[x$animal %in% phenotyped]

		## sires of those offspring that are parents
		sires <- 	x$sire[match(sibs,x$animal)]

		## sex of those that are parents
		sexes <- c(rep("F",length(sisters)),rep("M",length(brothers)))
		
		## number of *phenotyped* offspring they have
		off <- c(
			sapply(sisters, function(i) sum(d_s[[i]]$animal %in% phenotyped)),#nrow(d_s[[i]])),
			sapply(brothers, function(i) sum(s_s[[i]]$animal %in% phenotyped)),#nrow(s_s[[i]])),
			recursive=TRUE)

	## this will fail with full sib matings, becasue will double count

		if(length(all_sibs)>1 & length(sibs)>0){

			ex_fam <- do.call(rbind,lapply(1: length(sibs), function(i){
				data.frame(
					relationship="au",
					n = off[i],
					parent_r = ifelse(sires[i] == all_sires[which(all_sibs!=sibs[i])], "FS","MHS"),
					parent_sex = sexes[i]
					)
			}))

	## have to do the same for paternal half sibs - because linked through dad doesn't matter exact relation?
		}else{
			NULL
		}


	}))
	# x<-s_s[17]
	out_s<-do.call(rbind,lapply(s_s, function(x){
			## work out which offspring are themselves parents
		# for(j in 1:length(s_s)){
			# x<-s_s[[j]]

		sisters <- x$animal[x$animal %in% ped$dam]
		brothers <- x$animal[x$animal %in% ped$sire]
		sibs <- c(sisters,brothers)
		## aunt and uncles don't need to have offspring themselves

		# all_sibs <- x$animal
		# all_dams <- x$dam
		all_sibs <- x$animal[x$animal %in% phenotyped]
		all_dams <- x$dam[x$animal %in% phenotyped]

		## dams of those offspring that are parents
		dams <- 	x$dam[match(sibs,x$animal)]

		## sex of those that are parents
		sexes <- c(rep("F",length(sisters)),rep("M",length(brothers)))
		
		## number of *phenotyped* offspring they have
		off <- c(
			sapply(sisters, function(i) sum(d_s[[i]]$animal %in% phenotyped)),#nrow(d_s[[i]])),
			sapply(brothers, function(i) sum(s_s[[i]]$animal %in% phenotyped)),#nrow(s_s[[i]])),
			recursive=TRUE)
	## this will fail with full sib matings, because will double count

		if(length(all_sibs)>1 & length(sibs)>0){

			ex_fam <- do.call(rbind,lapply(1: length(sibs), function(i){
				data.frame(
					relationship="au",
					n = off[i],
					parent_r = ifelse(dams[i] == all_dams[which(all_sibs!=sibs[i])], "FS","PHS"),
					parent_sex = sexes[i]
					)
			}))

	## have to do the same for paternal half sibs - because linked through dad doesn't matter exact relation?
		}else{
			NULL
		}
# }

	}))

	out <- rbind(out_d,subset(out_s,parent_r!="FS"))
	c_sum <- aggregate(n~parent_sex+parent_r+relationship,out,sum)
	c_sum2 <- c_sum[,4]
	names(c_sum2) <- paste(c_sum[,3],c_sum[,2],c_sum[,1],sep="_")
	c_sum2
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


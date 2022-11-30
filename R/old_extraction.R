cousins <- function(ped, phenotyped=NULL){
	
	## number sharing grandparents?



	#ped=ped_bt
	colnames(ped) <- c("animal","dam","sire")

	## if dont specify phenotyped, then assume all phenotyped
	if(is.null(phenotyped)) phenotyped <- ped[,1]

	d_s<-split(ped,ped$dam)
	s_s<-split(ped,ped$sire)

	# x<-d_s[['3_27']]

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


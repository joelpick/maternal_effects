source("/Users/joelpick/github/squidPed/R/simulate_pedigree.R")

generations=5
n_females=100
fecundity=10

hs_ped <-	simulate_pedigree(
		years = generations,
		n_females = n_females,
		fecundity = fecundity,
		p_sire = 0.5, 				# mating system (0-1, 1= one male per female, 0=complete random mating)
		juv_surv = 2/fecundity, # insures no population growth
		adult_surv = 0,					# discrete generations
		immigration = 0, 				# closed population
		constant_pop = TRUE     # constant population size
		)$pedigree

cousins <- function(ped){
	#ped=ped_bt
	colnames(ped) <- c("animal","dam","sire")
	d_s<-split(ped,ped$dam)
	s_s<-split(ped,ped$sire)

	# x<-d_s[[1]]
	out_d<-do.call(rbind,lapply(d_s, function(x){
			## work out which offspring are themselves parents
		sisters <- x$animal[x$animal %in% ped$dam]
		brothers <- x$animal[x$animal %in% ped$sire]
		sibs <- c(sisters,brothers)
		
		## sires of those offspring that are parents
		sires <- 	x$sire[match(sibs,x$animal)]

		## sex of those that are parents
		sexes <- c(rep("F",length(sisters)),rep("M",length(brothers)))
		
		## number of offspring they have
		off <- c(sapply(sisters, function(i) nrow(d_s[[i]])),	sapply(brothers, function(i) nrow(s_s[[i]])), recursive=TRUE)

	## this will fail with full sib matings, becasue will double count

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
					),
					## Number of aunt/uncle relationships from siblings to focal individuals offspring
					data.frame(
						relationship="a_u",
						n = off[i],
						parent_r = ifelse(sires[(i+1):length(sibs)] %in% sires[i], "FS","MHS"),
						sex = paste0(sexes[(i+1):length(sibs)],sexes[i])
					),
					## Number of aunt/uncle relationships from focal individual to siblings offspring
					data.frame(
						relationship="a_u",
						n = off[(i+1):length(sibs)],
						parent_r = ifelse(sires[(i+1):length(sibs)] %in% sires[i], "FS","MHS"),
						sex = paste0(sexes[i],sexes[(i+1):length(sibs)])
					)
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
		
		## sires of those offspring that are parents
		dams <- 	x$dam[match(sibs,x$animal)]

		## sex of those that are parents
		sexes <- c(rep("F",length(sisters)),rep("M",length(brothers)))
		
		## number of offspring they have
		off <- c(sapply(sisters, function(i) nrow(d_s[[i]])),	sapply(brothers, function(i) nrow(s_s[[i]])), recursive=TRUE)
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
					),
					## Number of aunt/uncle relationships from siblings to focal individuals offspring
					data.frame(
						relationship="a_u",
						n = off[i],
						parent_r = ifelse(dams[(i+1):length(sibs)] %in% dams[i], "FS","PHS"),
						sex = paste0(sexes[(i+1):length(sibs)],sexes[i])
					),
					## Number of aunt/uncle relationships from focal individual to siblings offspring
					data.frame(
						relationship="a_u",
						n = off[(i+1):length(sibs)],
						parent_r = ifelse(dams[(i+1):length(sibs)] %in% dams[i], "FS","PHS"),
						sex = paste0(sexes[i],sexes[(i+1):length(sibs)])
					)
				)
			}))
	## have to do the same for paternal half sibs - because linked through dad doesn't matter exact relation?
	}else{
		NULL
	}


	}))

	out <- rbind(out_d,subset(out_s,parent_r!="FS"))
	out$sex <- ifelse(out$sex=="MF","FM",out$sex)
	aggregate(n~sex+parent_r+relationship,out,sum)
}

cousins(hs_ped)

head(out_s,10)
sum(subset(out_d, relationship=="cousin" & parent_r=="FS")$n)
sum(subset(out_s, relationship=="cousin" & parent_r=="FS")$n)
sum(subset(out_d, relationship=="a_u" & parent_r=="FS")$n)
sum(subset(out_s, relationship=="a_u" & parent_r=="FS")$n)



sum(subset(out_d, relationship=="cousin" & parent_r=="MHS")$n)
sum(subset(out_s, relationship=="cousin" & parent_r=="PHS")$n)


pedA<-nadiv::makeA(hs_ped[,1:3])
all_link<-sum(pedA[lower.tri(pedA)]>0)
all_link/total_links(hs_ped)

sum(subset(out_s, relationship=="cousin" & parent_r=="FS")$n)/all_link
sum(subset(out_d, relationship=="cousin" & parent_r=="MHS")$n)/all_link
sum(subset(out_s, relationship=="cousin" & parent_r=="PHS")$n)/all_link


pedantics::pedStatSummary(pedantics::pedigreeStats(hs_ped[,1:3], includeA=FALSE,lowMem=TRUE,graphicalReport=FALSE))

## total links
total_links<-function(ped) nrow(ped) * (nrow(ped) - 1) / 2

fs_m, fs_mf, fs_f, hs_m, hs_mf, hs_f
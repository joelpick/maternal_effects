rm(list=ls())

# that not accounting for maternal effects inflates Va has permeated through the literature,
# but that unmodelled maternal genetic effects does the same has not. For example, in a large meta-analysis of h2 and m2 (modelled almost exlusively 97.8% as Vme), no mention is made of this bias, even though is vmg is as ubiquitous as current studies show, it will apply in most cases, meanng that we are almost universally overestimating va and underestimating Vmg
# This may be because the best known examples in WAM are empirical 
# Simulation results in Clements has not gained traction in WAM (check citations)

# correlations - where does variance go (vmg) and so how does that affect total Va
# is cov picked up when not there

# total Va equation gives response when selection is maximal
# total Va equation assumes no selection on parental traits, which intuitively is unlikely, but has been empircally difficult to show

# if quail work biased?
# what happens to trait based models without vmg?


## social peds might et a double whammy
## inflation of Va due to full sib structure AND inflation of Va due to misassignment of sires

## bt ped seems better with me - more offspring per mother?



library(asreml)
library(parallel)

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

load(paste0(data_wd,"gaussian_data.Rdata"))

rm(fs_peds,hs_peds,fs_data,hs_data)
ls()
source(paste0(wd,"R/00_functions.R"))
 
n_sims <-100

ped_bt <-  MasterBayes::insertPed(read.csv(paste0(wd,"Data/Raw/ped_BT.csv"))[,1:3])
ped_bt_social <- MasterBayes::insertPed(read.csv(paste0(wd,"Data/Raw/ped_BT.csv"))[,c(1,2,4)])
ped_bt_socialD <- MasterBayes::insertPed(read.csv(paste0(wd,"Data/Raw/ped_BT.csv"))[,c(1,2,5)])
head(ped_bt)
head(ped_bt_social)
head(ped_bt_socialD)


# ped<-read.csv(paste0(wd,"Data/Raw/ped_BT.csv"))
# head(ped)
# mean(ped[,3]==ped[,4], na.rm=TRUE)
# apply(ped,2,function(x)sum(is.na(x)))

# sum(is.na(ped[,2])&is.na(ped[,3]))

# data_wd <- "~/Dropbox/0_blue_tits/skew/Data/Raw/"
# tMORPH <- read.csv(paste0(data_wd,"tMORPH.csv"))
# adult_males <- subset(tMORPH, agecode%in%c(5,6) & !is.na(nest) & morph_sex=="M")[,c("bird_id","nest")]
# tPED <- subset(read.csv(paste0(data_wd,"tPED.csv")), ped_type=="DSCM")[,1:4]
# tPED$social_male <- adult_males$bird_id[match(tPED$nest_orig,adult_males$nest)]
# ## fill in dummy males
# tPED$social_male_dummy <- ifelse(is.na(tPED$social_male), paste0(tPED$nest_orig,"M"), tPED$social_male)
# head(tPED)
# mean(tPED[,3]==tPED[5], na.rm=TRUE)
# mean(tPED[,3]==tPED[6], na.rm=TRUE)

# write.csv(tPED[,c(1:3,5,6)],file=paste0(wd,"Data/Raw/ped_BT.csv"),row.names=FALSE)
plot(table(table(ped_bt$dam)))
plot(table(table(ped_rd$dam)))
plot(table(table(ped_ss$dam)))

table(table(ped_bt$dam))/sum(length(unique(ped_bt$dam)))
table(table(ped_rd$dam))/sum(length(unique(ped_rd$dam)))
table(table(ped_ss$dam))/sum(length(unique(ped_ss$dam)))

maternal_fs_ped <- function(ped){ 
	ped<-na.omit(ped) ## might depend on whether unknown males paternity likely means unknown male - in which case excluding will increase full-sib
	t(sapply(split(ped,ped$dam), function(x) c(length(unique(x$sire)),nrow(x)) ))
}

## possibly number of offspring per female
## possibly need both half and full sibs?
## but 


# x<-subset(ped_rd,dam=="NOB76")
# Tri2M(fs, diag=FALSE)

rd1<-maternal_fs_ped(ped_rd)
rd1a <- rd1[,2]/rd1[,1]
rd2<-maternal_relatedness_ped(ped_rd)
plot(rd1a[names(rd1a) %in% names(rd2)]~rd2)

plot(rd1[rownames(rd1) %in% names(rd2),1]~rd2)


bt1<-maternal_fs_ped(ped_bt)
bt1a <- bt1[,1]/bt1[,2]
bt2<-maternal_relatedness_ped(ped_bt)
plot(bt1a[names(bt1a) %in% names(bt2)]~bt2)

plot(bt1[rownames(bt1) %in% names(bt2),1]~bt2)




ped_rd <- read.csv(paste0(wd,"Data/Raw/ped_RD.csv"))
ped_ss <- read.csv(paste0(wd,"Data/Raw/ped_SSH.csv"))


rd_data <- lapply(1:n_sims, function(i){
	list(
		a=mge_sim(ped_rd, param=scenarios["a",]),
		b=mge_sim(ped_rd, param=scenarios["b",]),
		c=mge_sim(ped_rd, param=scenarios["c",]),
		d=mge_sim(ped_rd, param=scenarios["d",]),
		e=mge_sim(ped_rd, param=scenarios["e",]),
		f=mge_sim(ped_rd, param=scenarios["f",]),
		g=mge_sim(ped_rd, param=scenarios["g",]),
		h=mge_sim(ped_rd, param=scenarios["h",]),
		i=mge_sim(ped_rd, param=scenarios["i",]),
		j=mge_sim(ped_rd, param=scenarios["j",]),
		k=mge_sim(ped_rd, param=scenarios["k",])
	)
})

model1_rd <- model_func(m1_func,ped_rd,rd_data)
model2_rd <- model_func(m2_func,ped_rd,rd_data)

ss_data <- lapply(1:n_sims, function(i){
	list(
		a=mge_sim(ped_ss, param=scenarios["a",]),
		b=mge_sim(ped_ss, param=scenarios["b",]),
		c=mge_sim(ped_ss, param=scenarios["c",]),
		d=mge_sim(ped_ss, param=scenarios["d",]),
		e=mge_sim(ped_ss, param=scenarios["e",]),
		f=mge_sim(ped_ss, param=scenarios["f",]),
		g=mge_sim(ped_ss, param=scenarios["g",]),
		h=mge_sim(ped_ss, param=scenarios["h",]),
		i=mge_sim(ped_ss, param=scenarios["i",]),
		j=mge_sim(ped_ss, param=scenarios["j",]),
		k=mge_sim(ped_ss, param=scenarios["k",])
	)
})

model1_ss <- model_func(m1_func,ped_ss,ss_data)
model2_ss <- model_func(m2_func,ped_ss,ss_data)

bt_data <- mclapply(1:n_sims, function(i){
	cat(i," ")
	list(
		a=mge_sim(ped_bt, param=scenarios["a",]),
		b=mge_sim(ped_bt, param=scenarios["b",]),
		c=mge_sim(ped_bt, param=scenarios["c",]),
		d=mge_sim(ped_bt, param=scenarios["d",]),
		e=mge_sim(ped_bt, param=scenarios["e",]),
		f=mge_sim(ped_bt, param=scenarios["f",]),
		g=mge_sim(ped_bt, param=scenarios["g",]),
		h=mge_sim(ped_bt, param=scenarios["h",]),
		i=mge_sim(ped_bt, param=scenarios["i",]),
		j=mge_sim(ped_bt, param=scenarios["j",]),
		k=mge_sim(ped_bt, param=scenarios["k",])
	)
}, mc.cores=6)

model1_bt <- model_func(m1_func,ped_bt,bt_data)
model2_bt <- model_func(m2_func,ped_bt,bt_data)


bt_social_data <- mclapply(1:n_sims, function(i){
	cat(i," ")
	list(
		a=mge_sim(ped_bt_social, param=scenarios["a",]),
		b=mge_sim(ped_bt_social, param=scenarios["b",]),
		c=mge_sim(ped_bt_social, param=scenarios["c",]),
		d=mge_sim(ped_bt_social, param=scenarios["d",]),
		e=mge_sim(ped_bt_social, param=scenarios["e",]),
		f=mge_sim(ped_bt_social, param=scenarios["f",]),
		g=mge_sim(ped_bt_social, param=scenarios["g",]),
		h=mge_sim(ped_bt_social, param=scenarios["h",]),
		i=mge_sim(ped_bt_social, param=scenarios["i",]),
		j=mge_sim(ped_bt_social, param=scenarios["j",]),
		k=mge_sim(ped_bt_social, param=scenarios["k",])
	)
}, mc.cores=6)

model1_bt_social <- model_func(m1_func,ped_bt_social,bt_social_data)
model2_bt_social <- model_func(m2_func,ped_bt_social,bt_social_data)

bt_socialD_data <- mclapply(1:n_sims, function(i){
	cat(i," ")
	list(
		a=mge_sim(ped_bt_socialD, param=scenarios["a",]),
		b=mge_sim(ped_bt_socialD, param=scenarios["b",]),
		c=mge_sim(ped_bt_socialD, param=scenarios["c",]),
		d=mge_sim(ped_bt_socialD, param=scenarios["d",]),
		e=mge_sim(ped_bt_socialD, param=scenarios["e",]),
		f=mge_sim(ped_bt_socialD, param=scenarios["f",]),
		g=mge_sim(ped_bt_socialD, param=scenarios["g",]),
		h=mge_sim(ped_bt_socialD, param=scenarios["h",]),
		i=mge_sim(ped_bt_socialD, param=scenarios["i",]),
		j=mge_sim(ped_bt_socialD, param=scenarios["j",]),
		k=mge_sim(ped_bt_socialD, param=scenarios["k",])
	)
}, mc.cores=6)

model1_bt_socialD <- model_func(m1_func,ped_bt_socialD,bt_socialD_data)
model2_bt_socialD <- model_func(m2_func,ped_bt_socialD,bt_socialD_data)




save(scenarios, model1_ss,model2_ss,model1_rd,model2_rd,model1_bt,model2_bt, model1_bt_social,model2_bt_social, model1_bt_socialD,model2_bt_socialD, file=paste0(data_wd,"real_sims.Rdata"))

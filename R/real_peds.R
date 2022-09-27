rm(list=ls())


library(asreml)
library(parallel)

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

load(paste0(data_wd,"gaussian_data.Rdata"))

rm(fs_peds,hs_peds,fs_data,hs_data)
ls()
source(paste0(wd,"R/00_functions.R"))
 
n_sims <-100

ped_bt <- read.csv(paste0(wd,"Data/Raw/ped_BT.csv"))[]
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

for()

model1_fs <- model_func(m1_func,ped_bt,bt_data)

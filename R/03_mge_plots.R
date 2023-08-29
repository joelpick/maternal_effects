
rm(list=ls())


scenarios[,"Va"] + 0.5*scenarios[,"Vmg"]+ 1.5*scenarios[,"r_amg"]*sqrt(scenarios[,"Va"] * scenarios[,"Vmg"])

library(asreml)
library(parallel)

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/extract_cousins.R"))
source(paste0(wd,"R/00_functions.R"))
source("/Users/joelpick/github/squidPed/R/simulate_pedigree.R")
# devtools::load_all("~/github/squidSim/R")


load(paste0(data_wd,"mge_sims3.Rdata"))


# ped_names2 <- ped_names[!grepl("fI",ped_names)]

# ped_str<-ped_str[ped_names2]
# ped_str_mat<-ped_str_mat[ped_names2]

# ped_sum<-sapply(ped_str,colMeans)
# ped_sum_mat<-sapply(ped_str_mat,colMeans)

mat_ratio_all<-sapply(ped_str,function(x){
	rowSums(x[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(x[,-(1:2)]) 
})

mat_ratio <- colMeans(mat_ratio_all)
mat_ratio_se <- apply(mat_ratio_all,2,se)


matM_ratio_all<-sapply(ped_str_mat,function(x){
	(x[,"mat_links"] - x[,"mat_sib"])/x[,"total_links"]
})

matM_ratio <- colMeans(matM_ratio_all)
matM_ratio_se <- apply(matM_ratio_all,2,se)


# mod1<-do.call(rbind,lapply(ped_names,function(k) {
# 	mod1 <- do.call(rbind,lapply(get(paste0("model1_",k)), function(x) {
# 		do.call(rbind,lapply(1:nrow(scenarios), function(i) 
# 			data.frame(
# 				r=k,
# 				scenario=i,
# 				Va_est = x[i,"A"],
# 				Vm_est = x[i,"Me"],
# 				Va_sim=scenarios[i,"Va"],
# 				Vm_sim =sum(scenarios[i,c("Vmg","Vme")]),
# 				Vmg_sim=scenarios[i,"Vmg"])))
# 	}))
# 	# assign(paste0("mod2_",k),mod2)
# }))

# va1<-aggregate(cbind(Vm_sim,Vm_est,Vm_bias)~ scenario+r, mod2,mean)

# plot(va1$Vm_bias,va2$Vm_bias,)


mod2<-do.call(rbind,lapply(ped_names,function(k) {
	mod2 <- do.call(rbind,lapply(get(paste0("model2_",k)), function(x) {
		do.call(rbind,lapply(1:nrow(scenarios), function(i) 
			data.frame(
				r=k,
				scenario=i,
				Va_est = x[i,"A"],
				Vm_est = x[i,"Me"],
				Va_sim=scenarios[i,"Va"],
				Vm_sim =sum(scenarios[i,c("Vmg","Vme")]),
				Vmg_sim=scenarios[i,"Vmg"])))
	}))
	# assign(paste0("mod2_",k),mod2)
}))

sapply(ped_names,function(k) {
	sum(sapply(get(paste0("model2_",k)),is.matrix))
	
})


#,sum(x[i,c("Mg","Me")])
head(mod2,20)
mod2$Va_bias <- mod2$Va_est - mod2$Va_sim
mod2$Vm_bias <- mod2$Vm_est - mod2$Vm_sim

# mod2$ln_Va_bias <- log(mod2$Va_bias)
va2<-aggregate(cbind(Va_bias,Vmg_sim,Vm_sim,Vm_bias)~ scenario+r, mod2,mean)
va2_se<-aggregate(cbind(Va_bias,Vm_bias)~ scenario+r, mod2,se)
r_order<- sapply(va2$r, function(x) which(ped_names==x))

# which(va1$r)

va2$mat_ratio <- mat_ratio[r_order]
va2$matM_ratio<- matM_ratio[r_order]

va2_se$mat_ratio <- mat_ratio_se[r_order]
va2_se$matM_ratio <- matM_ratio_se[r_order]


va2[,c("ms","fec","imm")] <- do.call(rbind,strsplit(va2$r,"_"))
# va2$matM_ratio2<- matM_ratio2[r_order]
scenarios
s <- c(1,3)
s <- c(1,2,5,7)
s <- c(2:4)
s <- c(4,6)
s <- c(5,8,9)
s<-1:9


plot_func <-function(s){
	dd<-subset(va2, scenario %in% s)
dd_se<-subset(va2_se, scenario %in% s)
cols<-viridis::viridis(length(s))
par(mfrow=c(2,2), mar=c(5,5,1,1), cex.lab=1.75, cex.axis=1.25 )
plot(Va_bias~ mat_ratio, dd, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~h^2), col=(cols)[as.factor(dd$scenario)], pch=c(15:17)[as.factor(dd$ms)])
arrows(dd$mat_ratio,dd$Va_bias+dd_se$Va_bias,dd$mat_ratio,dd$Va_bias-dd_se$Va_bias,code=3,angle=90,length=0.01, col=(cols)[as.factor(dd$scenario)])
# arrows(dd$mat_ratio+dd_se$mat_ratio,dd$Va_bias,dd$mat_ratio-dd_se$mat_ratio,dd$Va_bias,code=3,angle=90,length=0.01)
abline(h=0)


plot(Vm_bias~ mat_ratio, dd, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~m^2), col=(cols)[as.factor(dd$scenario)], pch=c(15:17)[as.factor(dd$fec)])
# arrows(dd$mat_ratio,dd$Va_bias+dd_se$Va_bias,dd$mat_ratio,dd$Va_bias-dd_se$Va_bias,code=3,angle=90,length=0.1)
arrows(dd$mat_ratio,dd$Vm_bias+dd_se$Vm_bias,dd$mat_ratio,dd$Vm_bias-dd_se$Vm_bias,code=3,angle=90,length=0.01, col=(cols)[as.factor(dd$scenario)])
abline(h=0)


plot(Vm_bias~ Va_bias, dd, cex=1, xlab=expression(Bias~"in"~h^2), ylab=expression(Bias~"in"~m^2), col=(cols)[as.factor(dd$scenario)], pch=c(15:18)[as.factor(dd$imm)])
abline(0,-0.5)
arrows(dd$Va_bias,dd$Vm_bias+dd_se$Vm_bias,dd$Va_bias,dd$Vm_bias-dd_se$Vm_bias,code=3,angle=90,length=0.01, col=(cols)[as.factor(dd$scenario)])
arrows(dd$Va_bias+dd_se$Va_bias,dd$Vm_bias,dd$Va_bias-dd_se$Va_bias,dd$Vm_bias,code=3,angle=90,length=0.01, col=(cols)[as.factor(dd$scenario)])

plot(NA, xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), xlab="",ylab="",bty="n")
legend("center",apply(scenarios[s,],1, function(x) paste(colnames(scenarios[s,]),"=",formatC(x,digits=2,format="f"), collapse=", ")), pch=19, col=cols, bty="n")


}

scenarios
setEPS()
pdf(paste0(wd,"Figures/mge_fig1.pdf"), height=10, width=10)
plot_func(1:9)
dev.off()

setEPS()
pdf(paste0(wd,"Figures/mge_fig2.pdf"), height=10, width=10)
plot_func(c(1,2,5:7))
dev.off()

setEPS()
pdf(paste0(wd,"Figures/mge_fig3.pdf"), height=10, width=10)
plot_func(c(2:4))
dev.off()

setEPS()
pdf(paste0(wd,"Figures/mge_fig4.pdf"), height=10, width=10)
plot_func(c(5,8,9))
dev.off()


plot_func(c(1,3))
plot_func(c(1,2,5:7))
plot_func(c(4,6))





plot(Va_bias~ matM_ratio, va2, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~h^2), col=(cols)[as.factor(va2$scenario)], pch=c(15:17)[as.factor(va2$ms)])
arrows(va2$matM_ratio,va2$Va_bias+va2_se$Va_bias,va2$matM_ratio,va2$Va_bias-va2_se$Va_bias,code=3,angle=90,length=0.1)
# plot(Va_bias~ matM_ratio2, va2, subset=scenario==1,  pch=19, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~h^2))








# par(mfrow=c(1,2), mar=c(5,5,1,1), cex.lab=1.75, cex.axis=1.25 )
# plot(Va_bias~ mat_ratio, va2, subset=scenario==1, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~h^2), col=(1:4)[as.factor(va2$imm)], pch=c(15:17)[as.factor(va2$ms)])
# arrows(va2$mat_ratio,va2$Va_bias+va2_se$Va_bias,va2$mat_ratio,va2$Va_bias-va2_se$Va_bias,code=3,angle=90,length=0.1)
# plot(Va_bias~ matM_ratio, va2, subset=scenario==1, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~h^2), col=(1:4)[as.factor(va2$imm)], pch=c(15:17)[as.factor(va2$ms)])
# arrows(va2$matM_ratio,va2$Va_bias+va2_se$Va_bias,va2$matM_ratio,va2$Va_bias-va2_se$Va_bias,code=3,angle=90,length=0.1)













total_Va<-do.call(rbind,lapply(ped_names,function(k) {
	mod2 <- do.call(rbind,lapply(get(paste0("model2_",k)), function(x) {
		do.call(rbind,lapply(1:nrow(scenarios), function(i) data.frame(
			r=k,
			scenario=i,
			max=x[i,"A"] + 0.5*x[i,"Me"], 
			min=x[i,"A"], 
			sim=scenarios[i,"Va"] + 0.5*scenarios[i,"Vmg"]+ 1.5*scenarios[i,"r_amg"]*sqrt(scenarios[i,"Va"] * scenarios[i,"Vmg"])
			)
		))
	}))}
	))
total_Va$coverage <- total_Va$max>total_Va$sim & total_Va$min<total_Va$sim
total_Va$bias <- total_Va$min - total_Va$sim

# hist(total_Va$bias)

# hist(as.numeric(total_Va$coverage))



tVa_means <- aggregate(cbind(coverage,bias)~ r+scenario, total_Va, mean)
r_order<- sapply(tVa_means$r, function(x) which(ped_names==x))
tVa_means$mat_ratio <- mat_ratio[r_order]
tVa_means$matM_ratio<- matM_ratio[r_order]
tVa_means[,c("ms","fec","imm")] <- do.call(rbind,strsplit(tVa_means$r,"_"))

setEPS()
pdf(paste0(wd,"Figures/mge_fig5.pdf"), height=6, width=6)
	par(mfrow=c(1,1), mar=c(5,5,1,1), cex.lab=1.75, cex.axis=1.25 )

	cols<-viridis::viridis(9)
	plot(bias~mat_ratio,tVa_means, pch=19, col=(cols)[as.factor(tVa_means$scenario)], ylab="Bias in Total Va", xlab="Proportion non-sibling maternal links");abline(h=0)

dev.off()


{
	cols<-viridis::viridis(9)

	par(mfrow=c(1,2), mar=c(5,5,1,1), cex.lab=1.75, cex.axis=1.25 )

	plot(bias~mat_ratio,tVa_means, pch=19, col=(cols)[as.factor(tVa_means$scenario)]);abline(h=0)
	# legend("topleft",apply(scenarios[,c(1,2,4)],1, function(x) paste(colnames(scenarios[,c(1,2,4)]),"=",x, collapse=",")), pch=19, col=cols)

	plot(coverage~mat_ratio,tVa_means, pch=19, col=(cols)[as.factor(tVa_means$scenario)])
}
## take home here is that an inference based on just Va where maternal genetic effects are likely is not very informative, as we are likely underestmating, but be overesitmaint g  


hist(tVa_means$mat_ratio)

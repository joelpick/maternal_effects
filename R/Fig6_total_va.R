# This script generates figure 6

rm(list=ls())


wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))

load(paste0(data_wd,"mge_sims3.Rdata"))



mat_ratio_all<-sapply(ped_str,function(x){
	rowSums(x[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(x[,-(1:2)]) 
})

mat_ratio <- colMeans(mat_ratio_all)
mat_ratio_se <- apply(mat_ratio_all,2,se)


total_Va<-do.call(rbind,lapply(ped_names,function(k) {
	mod2 <- do.call(rbind,lapply(get(paste0("model2_",k)), function(x) {
		do.call(rbind,lapply(1:nrow(scenarios), function(i) data.frame(
			r=k,
			scenario=i,
			max=x[["ml"]][i,"A"] + 0.5*x[["ml"]][i,"Me"], 
			min=x[["ml"]][i,"A"], 
			sim=scenarios[i,"Va"] + 0.5*scenarios[i,"Vmg"]+ 1.5*scenarios[i,"r_amg"]*sqrt(scenarios[i,"Va"] * scenarios[i,"Vmg"])
			)
		))
	}))}
	))
total_Va$bias <- total_Va$min - total_Va$sim


tVa_means <- aggregate(bias~ r+scenario, total_Va, mean)
r_order<- sapply(tVa_means$r, function(x) which(ped_names==x))
tVa_means$mat_ratio <- mat_ratio[r_order]
tVa_means[,c("ms","fec","imm")] <- do.call(rbind,strsplit(tVa_means$r,"_"))

setEPS()
pdf(paste0(wd,"Figures/Fig6_totalVa.pdf"), height=6, width=8)
{	

	layout(matrix(c(1,2),nrow=1, byrow=TRUE), width=c(10,1))

	par( mar=c(5,5,1,1), cex.lab=1.5, cex.axis=1.1 )

	cols<-c(palette.colors(),1)
	bgs= c(palette.colors(),0)
	pch=rep(21:25,2)
	plot(bias~mat_ratio,tVa_means, pch=pch[as.factor(tVa_means$scenario)], col=(cols)[as.factor(tVa_means$scenario)],bg=(bgs)[as.factor(tVa_means$scenario)], ylab=expression(Bias~"in"~italic(hat(V)[At])), xlab="Proportion non-sibling maternal links");abline(h=0)

	par(mar=c(0,0,0,0))

	plot(NA, xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), xlab="",ylab="",bty="n")
	legend("center",letters[1:10], pch=pch, pt.bg=bgs, col=cols, bty="n", cex=1,title="Scenario")



}
dev.off()
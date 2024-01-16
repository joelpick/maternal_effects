
rm(list=ls())

library(asreml)
library(parallel)

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/extract_cousins.R"))
source(paste0(wd,"R/00_functions.R"))
source("/Users/joelpick/github/squidPed/R/simulate_pedigree.R")
# devtools::load_all("~/github/squidSim/R")


load(paste0(data_wd,"mge_sims3.Rdata"))



mat_ratio_all<-sapply(ped_str,function(x){
	rowSums(x[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(x[,-(1:2)]) 
})

mat_ratio <- colMeans(mat_ratio_all)
mat_ratio_se <- apply(mat_ratio_all,2,se)

# matM_ratio_all<-sapply(ped_str_mat,function(x){
# 	(x[,"mat_links"] - x[,"mat_sib"])/x[,"total_links"]
# })

# matM_ratio <- colMeans(matM_ratio_all)
# matM_ratio_se <- apply(matM_ratio_all,2,se)


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
# tVa_means$matM_ratio<- matM_ratio[r_order]
tVa_means[,c("ms","fec","imm")] <- do.call(rbind,strsplit(tVa_means$r,"_"))

setEPS()
pdf(paste0(wd,"Figures/fig4_totalVa.pdf"), height=6, width=12)
{	

		layout(matrix(c(1,2),nrow=1, byrow=TRUE), width=c(3,2))

	par( mar=c(5,5,1,1), cex.lab=1.75, cex.axis=1.25 )

	cols<-viridis::viridis(10)
	pch=rep(21:25,2)
	plot(bias~mat_ratio,tVa_means, pch=pch[as.factor(tVa_means$scenario)], col=(cols)[as.factor(tVa_means$scenario)],bg=(cols)[as.factor(tVa_means$scenario)], ylab="Bias in Total Va", xlab="Proportion non-sibling maternal links");abline(h=0)

	par(mar=c(0,0,0,0))
	# scenarios2 <- formatC(scenarios,digits=2,format="f")
	# scenarios2[scenarios2!="0.00"] <- paste0("bold(",scenarios2[scenarios2!="0.00"],")")
	# s=1:2

	# x<-parse(text="a")
# expression(bquote(.(x)^2))

	legend_text<-apply(scenarios,1, function(x) paste(colnames(scenarios),"=",formatC(x,digits=2,format="f"), collapse=", "))
	# legend_text<-(c(apply(scenarios2[s,,drop=FALSE],1, function(x) paste(colnames(scenarios2[s,,drop=FALSE]),"=",x, collapse=", ")),recursive=TRUE))

	plot(NA, xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), xlab="",ylab="",bty="n")
	legend("center",legend_text, pch=pch, pt.bg=cols, col=cols, bty="n", cex=1)



}
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


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
			max=x[["ml"]][i,"A"] + 0.5*x[["ml"]][i,"Me"], 
			min=x[["ml"]][i,"A"], 
			sim=scenarios[i,"Va"] + 0.5*scenarios[i,"Vmg"]+ 1.5*scenarios[i,"r_amg"]*sqrt(scenarios[i,"Va"] * scenarios[i,"Vmg"])
			)
		))
	}))}
	))
total_Va$bias <- total_Va$min - total_Va$sim

# hist(total_Va$bias)

# hist(as.numeric(total_Va$coverage))



tVa_means <- aggregate(cbind(coverage,bias)~ r+scenario, total_Va, mean)
r_order<- sapply(tVa_means$r, function(x) which(ped_names==x))
tVa_means$mat_ratio <- mat_ratio[r_order]
# tVa_means$matM_ratio<- matM_ratio[r_order]
tVa_means[,c("ms","fec","imm")] <- do.call(rbind,strsplit(tVa_means$r,"_"))

setEPS()
pdf(paste0(wd,"Figures/Fig6_totalVa.pdf"), height=6, width=8)
{	

	layout(matrix(c(1,2),nrow=1, byrow=TRUE), width=c(10,1))

	par( mar=c(5,5,1,1), cex.lab=1.5, cex.axis=1.1 )

	cols<-c(palette.colors(),1)
	bgs= c(palette.colors(),0)
	pch=rep(21:25,2)
	plot(bias~mat_ratio,tVa_means, pch=pch[as.factor(tVa_means$scenario)], col=(cols)[as.factor(tVa_means$scenario)],bg=(bgs)[as.factor(tVa_means$scenario)], ylab=expression(Bias~in~Total~V[A]), xlab="Proportion non-sibling maternal links");abline(h=0)

	par(mar=c(0,0,0,0))
	# scenarios2 <- formatC(scenarios,digits=2,format="f")
	# scenarios2[scenarios2!="0.00"] <- paste0("bold(",scenarios2[scenarios2!="0.00"],")")
	# s=1:2

	# x<-parse(text="a")
# expression(bquote(.(x)^2))

	# legend_text<-apply(scenarios,1, function(x) paste(colnames(scenarios),"=",formatC(x,digits=2,format="f"), collapse=", "))
	# legend_text<-(c(apply(scenarios2[s,,drop=FALSE],1, function(x) paste(colnames(scenarios2[s,,drop=FALSE]),"=",x, collapse=", ")),recursive=TRUE))

	plot(NA, xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), xlab="",ylab="",bty="n")
	legend("center",LETTERS[1:10], pch=pch, pt.bg=bgs, col=cols, bty="n", cex=1,title="Scenario")



}
dev.off()
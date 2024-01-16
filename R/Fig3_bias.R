
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


# ped_names2 <- ped_names[!grepl("fI",ped_names)]

# ped_str<-ped_str[ped_names2]
# ped_str_mat<-ped_str_mat[ped_names2]

# ped_sum<-sapply(ped_str,colMeans)
# ped_sum_mat<-sapply(ped_str_mat,colMeans)

mat_ratio_all<-sapply(ped_str,function(x){
	rowSums(x[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS")])/rowSums(x[,-(1:2)]) 
})

matsib_ratio_all<-sapply(ped_str,function(x){
	rowSums(x[,c("FS","MHS")])/rowSums(x[,-(1:2)]) 
})

cov_ratio_all<-sapply(ped_str,function(x){
	rowSums(x[,c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS",
		"sire","PG","au_S_FS","au_S_MHS","au_D_PHS","cousin_DS_FS","cousin_DS_HS"
		)])/rowSums(x[,-(1:2)]) 
})

mat_ratio <- colMeans(mat_ratio_all)
mat_ratio_se <- apply(mat_ratio_all,2,se)

matsib_ratio <- colMeans(matsib_ratio_all)


matM_ratio_all<-sapply(ped_str_mat,function(x){
	(x[,"mat_links"] - x[,"mat_sib"])/x[,"total_links"]
})

matM_ratio <- colMeans(matM_ratio_all)
matM_ratio_se <- apply(matM_ratio_all,2,se)



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


head(mod2,20)
mod2$Va_bias <- mod2$Va_est - mod2$Va_sim
mod2$Vm_bias <- mod2$Vm_est - mod2$Vm_sim
r_order2<- sapply(mod2$r, function(x) which(ped_names==x))

mod2$mat_ratio <- mat_ratio[r_order2]
mod2$matsib_ratio <- matsib_ratio[r_order2]

summary(lme4::lmer(Va_bias ~ mat_ratio +matsib_ratio + (1|r),mod2, subset=scenario=="1"))

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

s=1:2


plot_func <-function(s, legend_parts=1:4, legend_order=1:length(s), lines=TRUE, cols=viridis::viridis(10), pchs=rep(21:25,2), Va_lim=c(-0.11,0.3), Vm_lim=c(-0.17,0.06)){
	dd<-subset(va2, scenario %in% s)
	dd_se<-subset(va2_se, scenario %in% s)
	cols<-viridis::viridis(length(s))

	# layout(matrix(c(1,2,3,4,4,4),nrow=2, byrow=TRUE), height=c(5,2))
	par(mar=c(5,5,1,1), cex.lab=1.75, cex.axis=1.25 )
	plot(Va_bias~ mat_ratio, dd, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~V[A]), 
		pch=pchs[as.factor(dd$scenario)], 
		col=cols[as.factor(dd$scenario)],
		bg=cols[as.factor(dd$scenario)],
		ylim=Va_lim)
	arrows(dd$mat_ratio,dd$Va_bias+dd_se$Va_bias,dd$mat_ratio,dd$Va_bias-dd_se$Va_bias,code=3,angle=90,length=0.01, col=(cols)[
		as.factor(dd$scenario)])
	# arrows(dd$mat_ratio+dd_se$mat_ratio,dd$Va_bias,dd$mat_ratio-dd_se$mat_ratio,dd$Va_bias,code=3,angle=90,length=0.01)
	abline(h=0)

	if(lines){
		coefsA<-sapply(s,function(i)(coef(lm(Va_bias~mat_ratio,dd,subset=scenario==i))))
		sapply(1:length(s),function(x) abline(coefsA[1,x],coefsA[2,x],col=cols[x], lty=2))
	}


	plot(Vm_bias~ mat_ratio, dd, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~V[M]), 
		pch=pchs[as.factor(dd$scenario)], 
		col=cols[as.factor(dd$scenario)],
		bg=cols[as.factor(dd$scenario)],
		ylim=Vm_lim)
	# arrows(dd$mat_ratio,dd$Va_bias+dd_se$Va_bias,dd$mat_ratio,dd$Va_bias-dd_se$Va_bias,code=3,angle=90,length=0.1)
	arrows(dd$mat_ratio,dd$Vm_bias+dd_se$Vm_bias,dd$mat_ratio,dd$Vm_bias-dd_se$Vm_bias,code=3,angle=90,length=0.01, col=(cols)[as.factor(dd$scenario)])
	abline(h=0)

	if(lines){
		coefsM<-sapply(s,function(i)(coef(lm(Vm_bias~mat_ratio,dd,subset=scenario==i))))
		sapply(1:length(s),function(x) abline(coefsM[1,x],coefsM[2,x],col=cols[x], lty=2))
	}


	plot(Vm_bias~ Va_bias, dd, cex=1, xlab=expression(Bias~"in"~V[A]), ylab=expression(Bias~"in"~V[M]), 
		pch=pchs[as.factor(dd$scenario)], 
		col=cols[as.factor(dd$scenario)],
		bg=cols[as.factor(dd$scenario)],
		ylim=Vm_lim,
		xlim=Va_lim)
	abline(0,-0.5)
	arrows(dd$Va_bias,dd$Vm_bias+dd_se$Vm_bias,dd$Va_bias,dd$Vm_bias-dd_se$Vm_bias,code=3,angle=90,length=0.01, col=(cols)[as.factor(dd$scenario)])
	arrows(dd$Va_bias+dd_se$Va_bias,dd$Vm_bias,dd$Va_bias-dd_se$Va_bias,dd$Vm_bias,code=3,angle=90,length=0.01, col=(cols)[as.factor(dd$scenario)])
	legend("topright",expression(m^2~"="~"-"*0.5*h^2),lty=1,bty="n")

	# par(mar=c(0,0,0,0))
	# scenarios2 <- formatC(scenarios,digits=2,format="f")
	# scenarios2[scenarios2!="0.00"] <- paste0("bold(",scenarios2[scenarios2!="0.00"],")")
	# s=1:2

	# x<-parse(text="a")
# expression(bquote(.(x)^2))

	# legend_text<-apply(scenarios[s,legend_parts,drop=FALSE],1, function(x) paste(colnames(scenarios[s,legend_parts,drop=FALSE]),"=",formatC(x,digits=2,format="f"), collapse=", "))
	# # legend_text<-(c(apply(scenarios2[s,,drop=FALSE],1, function(x) paste(colnames(scenarios2[s,,drop=FALSE]),"=",x, collapse=", ")),recursive=TRUE))

	# plot(NA, xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), xlab="",ylab="",bty="n")
	# legend("center",
	# legend_text[legend_order]
	# 	, pch=19, col=cols[legend_order], bty="n", cex=2)
}

scenarios

par(mfrow=c(4,3))

plot_func(1, 2)
plot_func(c(2,5,6),2:3,c(2,1,3))
plot_func(c(1:4), 1:3)
plot_func(c(3,7:10),4,c(2:3,1,4:5))


setEPS()
pdf(paste0(wd,"Figures/mge_fig1.pdf"), height=6, width=15)

dev.off()

setEPS()
pdf(paste0(wd,"Figures/mge_fig2.pdf"), height=6, width=15)
plot_func(c(1,4), 1:3)
dev.off()

setEPS()
pdf(paste0(wd,"Figures/mge_fig3.pdf"), height=6, width=15)
plot_func(c(2,5,6),2:3,c(2,1,3))
dev.off()

setEPS()
pdf(paste0(wd,"Figures/mge_fig4.pdf"), height=6, width=15)
plot_func(c(3,7:10),4,c(2:3,1,4:5))
dev.off()


setEPS()
pdf(paste0(wd,"Figures/cov_mge.pdf"), height=5, width=5)
par(mar=c(6,6,1,1))
plot(cov_ratio~mat_ratio, pch=19, xlab="Proportion non-sibling maternal links (Vmg)", ylab="Proportion non-sibling links \nthrough single mother (COVa,mg)")
dev.off()


plot(Va_bias~ matM_ratio, va2, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~h^2), col=(cols)[as.factor(va2$scenario)], pch=c(15:17)[as.factor(va2$ms)])
arrows(va2$matM_ratio,va2$Va_bias+va2_se$Va_bias,va2$matM_ratio,va2$Va_bias-va2_se$Va_bias,code=3,angle=90,length=0.1)
# plot(Va_bias~ matM_ratio2, va2, subset=scenario==1,  pch=19, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~h^2))








# par(mfrow=c(1,2), mar=c(5,5,1,1), cex.lab=1.75, cex.axis=1.25 )
# plot(Va_bias~ mat_ratio, va2, subset=scenario==1, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~h^2), col=(1:4)[as.factor(va2$imm)], pch=c(15:17)[as.factor(va2$ms)])
# arrows(va2$mat_ratio,va2$Va_bias+va2_se$Va_bias,va2$mat_ratio,va2$Va_bias-va2_se$Va_bias,code=3,angle=90,length=0.1)
# plot(Va_bias~ matM_ratio, va2, subset=scenario==1, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~h^2), col=(1:4)[as.factor(va2$imm)], pch=c(15:17)[as.factor(va2$ms)])
# arrows(va2$matM_ratio,va2$Va_bias+va2_se$Va_bias,va2$matM_ratio,va2$Va_bias-va2_se$Va_bias,code=3,angle=90,length=0.1)




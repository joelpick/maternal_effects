
rm(list=ls())

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))

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
				Va_est = x[["ml"]][i,"A"],
				Vm_est = x[["ml"]][i,"Me"],
				Va_sim=scenarios[i,"Va"],
				Vm_sim =sum(scenarios[i,c("Vmg","Vme")]),
				Vmg_sim=scenarios[i,"Vmg"])))
	}))
	# assign(paste0("mod2_",k),mod2)
}))

nrow(mod2)
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



plot_func <-function(s, legend_parts=1:4, legend_order=1:length(s), lines=TRUE, cols=viridis::viridis(10), pchs=rep(21:25,2), Va_lim=c(-0.11,0.3), Vm_lim=c(-0.17,0.06)){
	dd<-subset(va2, scenario %in% s)
	dd_se<-subset(va2_se, scenario %in% s)
	# cols<-viridis::viridis(length(s))

 	pch<- pchs[dd$scenario] 
	col<-cols[dd$scenario]
		

	# layout(matrix(c(1,2,3,4,4,4),nrow=2, byrow=TRUE), height=c(5,2))
	par(mar=c(5,5,1,1), cex.lab=1.75, cex.axis=1.25 )
	plot(Va_bias~ mat_ratio, dd, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~V[A]), 
		pch=pch, 
		col=col,
		bg=col,
		ylim=Va_lim)
	arrows(dd$mat_ratio,dd$Va_bias+dd_se$Va_bias,dd$mat_ratio,dd$Va_bias-dd_se$Va_bias,code=3,angle=90,length=0.01, col=col)
	# arrows(dd$mat_ratio+dd_se$mat_ratio,dd$Va_bias,dd$mat_ratio-dd_se$mat_ratio,dd$Va_bias,code=3,angle=90,length=0.01)
	abline(h=0)

	if(lines){
		coefsA<-sapply(s,function(i)(coef(lm(Va_bias~mat_ratio,dd,subset=scenario==i))))
		sapply(1:length(s),function(x) abline(coefsA[1,x],coefsA[2,x],col=cols[s[x]], lty=2))
	}


	plot(Vm_bias~ mat_ratio, dd, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~V[M]), 
		pch=pch, 
		col=col,
		bg=col,
		ylim=Vm_lim)
	# arrows(dd$mat_ratio,dd$Va_bias+dd_se$Va_bias,dd$mat_ratio,dd$Va_bias-dd_se$Va_bias,code=3,angle=90,length=0.1)
	arrows(dd$mat_ratio,dd$Vm_bias+dd_se$Vm_bias,dd$mat_ratio,dd$Vm_bias-dd_se$Vm_bias,code=3,angle=90,length=0.01, col=col)
	abline(h=0)

	if(lines){
		coefsM<-sapply(s,function(i)(coef(lm(Vm_bias~mat_ratio,dd,subset=scenario==i))))
		sapply(1:length(s),function(x) abline(coefsM[1,x],coefsM[2,x],col=cols[s[x]], lty=2))
	}


	plot(Vm_bias~ Va_bias, dd, cex=1, xlab=expression(Bias~"in"~V[A]), ylab=expression(Bias~"in"~V[M]), 
		pch=pch, 
		col=col,
		bg=col,
		ylim=Vm_lim,
		xlim=Va_lim)
	abline(0,-0.5)
	arrows(dd$Va_bias,dd$Vm_bias+dd_se$Vm_bias,dd$Va_bias,dd$Vm_bias-dd_se$Vm_bias,code=3,angle=90,length=0.01, col=col)
	arrows(dd$Va_bias+dd_se$Va_bias,dd$Vm_bias,dd$Va_bias-dd_se$Va_bias,dd$Vm_bias,code=3,angle=90,length=0.01, col=col)
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




setEPS()
pdf(paste0(wd,"Figures/fig3_bias.pdf"), height=20, width=15)

{
	layout(matrix(1:15, byrow=TRUE, ncol=3,nrow=5),height=c(1,1,1,1,0.5))
	# par(mfrow=c(4,3))
	cols <- viridis::viridis(11)[c(10,5,6, 1 ,2,9, 1,4,8,11)]
	pchs <- rep(21:25,2)

	legend_text<-apply(scenarios,1, function(x) paste(colnames(scenarios),"=",formatC(x,digits=2,format="f"), collapse=", "))
	# legend_text<-(c(apply(scenarios2[s,,drop=FALSE],1, function(x) paste(colnames(scenarios2[s,,drop=FALSE]),"=",x, collapse=", ")),recursive=TRUE))

	legend_text2<-paste0(LETTERS[1:9], ": ", legend_text[c(5,2,6,1,3,7:10)])
	col2 <- cols[c(5,2,6,1,3,7:10)]
	pch2 <- pchs[c(5,2,6,1,3,7:10)]

	plot_func(c(2,5,6), cols=cols)
	plot_func(c(1,2), cols=cols)
	plot_func(c(1,3), cols=cols)
	plot_func(c(3,7:10), cols=cols)
par(mar=c(0,0,0,0))
	plot(NA, xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), xlab="",ylab="",bty="n")
	legend("center", legend_text2[1:3], pch=pch2[1:3], col=col2[1:3],pt.bg=col2[1:3], bty="n", cex=1.5)

	plot(NA, xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), xlab="",ylab="",bty="n")
	legend("center", legend_text2[4:6], pch=pch2[4:6], col=col2[4:6],pt.bg=col2[4:6], bty="n", cex=1.5)

		plot(NA, xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), xlab="",ylab="",bty="n")
	legend("center", legend_text2[7:9], pch=pch2[7:9], col=col2[7:9],pt.bg=col2[7:9], bty="n", cex=1.5)
}
dev.off()



	col2 <- c(palette.colors(10),palette.colors()[5])#cols[c(5,2,6,1,3,7:10)]
	pchs <- rep(21:25,2)
	legend_text<-apply(scenarios[,1:3],1, function(x) paste(colnames(scenarios[,1:3]),"=",formatC(x,digits=2,format="f"), collapse=", "))


	par(mfrow=c(3,3))
	plot_func(c(1:3), cols=col2)
	legend("bottomleft", legend_text[c(2,5,6)], pch=pchs[c(2,5,6)], col=col2[c(2,5,6)],pt.bg=col2[c(2,5,6)], bty="n", cex=1.25)

	plot_func(c(4:5), cols=col2)
		legend("bottomleft", legend_text[c(1,4)], pch=pchs[c(1,4)], col=col2[c(1,4)],pt.bg=col2[c(1,4)], bty="n", cex=1.25)

	legend_text<-apply(scenarios[,c(1:2,4)],1, function(x) paste(colnames(scenarios[,c(1:2,4)]),"=",formatC(x,digits=2,format="f"), collapse=", "))
	plot_func(c(6:10), cols=col2)
	legend("bottomright", legend_text[c(3,7:10)], pch=pchs[c(3,7:10)], col=col2[c(3,7:10)],pt.bg=col2[c(3,7:10)], box.col="white", cex=1.25,bg="white")




plot_func2 <-function(s, lines=TRUE, Va_lim=c(-0.11,0.3), Vm_lim=c(-0.17,0.06)){
	dd<-subset(va2, scenario %in% s)
	dd_se<-subset(va2_se, scenario %in% s)
	# cols<-viridis::viridis(length(s))

	pch= rep(21:25,2)[dd$scenario]
	col= c(palette.colors(),1)[dd$scenario]
	bg= c(palette.colors(),"white")[dd$scenario]

	# layout(matrix(c(1,2,3,4,4,4),nrow=2, byrow=TRUE), height=c(5,2))
	plot(Va_bias~ mat_ratio, dd, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~V[A]), 
		pch=pch, 
		col=col,
		bg=bg,
		ylim=Va_lim)
	arrows(dd$mat_ratio,dd$Va_bias+dd_se$Va_bias,dd$mat_ratio,dd$Va_bias-dd_se$Va_bias,code=3,angle=90,length=0.01, col=col)
	# arrows(dd$mat_ratio+dd_se$mat_ratio,dd$Va_bias,dd$mat_ratio-dd_se$mat_ratio,dd$Va_bias,code=3,angle=90,length=0.01)
	abline(h=0)

	legend("topleft",letters[s],
		pch= rep(21:25,2)[s], 
		col= c(palette.colors(),1)[s],
		pt.bg= c(palette.colors(),"white")[s],bty="n",title="Scenario")

	if(lines){
		coefsA<-sapply(s,function(i)(coef(lm(Va_bias~mat_ratio,dd,subset=scenario==i))))
		sapply(1:length(s),function(x) abline(coefsA[1,x],coefsA[2,x],col=cols[s[x]], lty=2))
	}


	plot(Vm_bias~ mat_ratio, dd, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~V[M]), 
		pch=pch, 
		col=col,
		bg=col,
		ylim=Vm_lim)
	# arrows(dd$mat_ratio,dd$Va_bias+dd_se$Va_bias,dd$mat_ratio,dd$Va_bias-dd_se$Va_bias,code=3,angle=90,length=0.1)
	arrows(dd$mat_ratio,dd$Vm_bias+dd_se$Vm_bias,dd$mat_ratio,dd$Vm_bias-dd_se$Vm_bias,code=3,angle=90,length=0.01, col=col)
	abline(h=0)

	if(lines){
		coefsM<-sapply(s,function(i)(coef(lm(Vm_bias~mat_ratio,dd,subset=scenario==i))))
		sapply(1:length(s),function(x) abline(coefsM[1,x],coefsM[2,x],col=cols[s[x]], lty=2))
	}

}


plot_func3 <-function(s, Va_lim=c(-0.11,0.3), Vm_lim=c(-0.17,0.06), legend_pos="bottomleft"){
	dd<-subset(va2, scenario %in% s)
	dd_se<-subset(va2_se, scenario %in% s)
	
cols <- c(palette.colors(),1)[dd$scenario]

	plot(Vm_bias~ Va_bias, dd, cex=1, xlab=expression(Bias~"in"~V[A]), ylab=expression(Bias~"in"~V[M]), 
		pch= rep(21:25,2)[dd$scenario], 
		col= cols,
		bg= c(palette.colors(),"white")[dd$scenario],
		ylim=Vm_lim,
		xlim=Va_lim)
	abline(0,-0.5)
	arrows(dd$Va_bias,dd$Vm_bias+dd_se$Vm_bias,dd$Va_bias,dd$Vm_bias-dd_se$Vm_bias,code=3,angle=90,length=0.01, col=cols)
	arrows(dd$Va_bias+dd_se$Va_bias,dd$Vm_bias,dd$Va_bias-dd_se$Va_bias,dd$Vm_bias,code=3,angle=90,length=0.01, col=cols)
	legend(legend_pos,,letters[s],
		pch= rep(21:25,2)[s], 
		col= c(palette.colors(),1)[s],
		pt.bg= c(palette.colors(),"white")[s],bg="white",box.col=0,title="Scenario")


}



### add in one that is open symbols


setEPS()
pdf(paste0(wd,"Figures/Fig3_bias.pdf"), height=10, width=8)

{	
	par(mfrow=c(3,2),	mar=c(5,5,1,1), cex.lab=1.5, cex.axis=1.1)

	plot_func2(c(1:3))
	plot_func2(c(4:5))
	plot_func2(c(6:10))
}
dev.off()

setEPS()
pdf(paste0(wd,"Figures/Fig4_Va_Vm.pdf"), height=5, width=10)

{
	par(mfrow=c(1,2),	mar=c(5,5,1,1), cex.lab=1.5, cex.axis=1.1)
	plot_func3(c(1:5))
	legend("topright",expression(m^2~"="~"-"*0.5*h^2),lty=1,bty="n")

	plot_func3(c(6:10), legend_pos="bottomright")
}
dev.off()
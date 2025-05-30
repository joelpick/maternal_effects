
rm(list=ls())

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

source(paste0(wd,"R/00_functions.R"))

load(paste0(data_wd,"mge_sims3.Rdata"))

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


mod2<-do.call(rbind,lapply(ped_names,function(k) {
	mod2 <- do.call(rbind,lapply(get(paste0("model2_",k)), function(x) {
			data.frame(
				r=k,
				scenario=1:nrow(scenarios),
				Va_est = x[["ml"]][,"A"],
				Vm_est = x[["ml"]][,"Me"],
				Va_sim=scenarios[,"Va"],
				Vm_sim = rowSums(scenarios[,c("Vmg","Vme")]),
				Vmg_sim=scenarios[,"Vmg"])
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

summary(lme4::lmer(Va_bias ~ mat_ratio  + (1|r),mod2, subset=scenario=="1"))

summary(lme4::lmer(Va_bias ~ 1+ (1|r),mod2, subset=scenario=="1"))


# mod2$ln_Va_bias <- log(mod2$Va_bias)
va2<-aggregate(cbind(Va_bias,Vmg_sim,Vm_sim,Vm_bias)~ scenario+r, mod2,mean)
va2_se<-aggregate(cbind(Va_bias,Vm_bias)~ scenario+r, mod2,se)
r_order<- sapply(va2$r, function(x) which(ped_names==x))

# which(va1$r)

va2$mat_ratio <- mat_ratio[r_order]
va2_se$mat_ratio <- mat_ratio_se[r_order]

va2$matsib_ratio <- matsib_ratio[r_order]


va2[,c("ms","fec","imm")] <- do.call(rbind,strsplit(va2$r,"_"))
# va2$matM_ratio2<- matM_ratio2[r_order]
plot(Va_bias ~ matsib_ratio,va2, col=va2$scenario, pch=19)


plot(matsib_ratio)


plot_fig3 <-function(s, lines=TRUE, Va_lim=c(-0.11,0.35), Vm_lim=c(-0.17,0.06),cols=c(palette.colors(),1),pchs=rep(21:25,2), bgs=c(palette.colors(),"white"), labels=c("","") , label_at=0.09){
	dd<-subset(va2, scenario %in% s)
	dd_se<-subset(va2_se, scenario %in% s)
	# cols<-viridis::viridis(length(s))
	pch= pchs[dd$scenario]
	col= cols[dd$scenario]
	bg= bgs[dd$scenario]

	# layout(matrix(c(1,2,3,4,4,4),nrow=2, byrow=TRUE), height=c(5,2))
	plot(Va_bias~ mat_ratio, dd, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~italic(hat(V)[A])), 
		pch=pch, 
		col=col,
		bg=bg,
		ylim=Va_lim)
	arrows(dd$mat_ratio,dd$Va_bias+dd_se$Va_bias,dd$mat_ratio,dd$Va_bias-dd_se$Va_bias,code=3,angle=90,length=0.01, col=col)
	# arrows(dd$mat_ratio+dd_se$mat_ratio,dd$Va_bias,dd$mat_ratio-dd_se$mat_ratio,dd$Va_bias,code=3,angle=90,length=0.01)
	abline(h=0)

	legend("topleft",letters[s],
		pch= pchs[s], 
		col= cols[s],
		pt.bg= bgs[s],bty="n",title="Scenario")

	if(lines){
		coefsA <- sapply(s,function(i){
			coefs <- coef(lm(Va_bias~mat_ratio,dd,subset=scenario==i))
			preds <- data.frame(x=range(dd$mat_ratio))
			preds$y <- coefs[1] + coefs[2]*preds$x
			print(preds)
			lines(y~x,preds,col=cols[i], lty=2)
		})
	}
	mtext(labels[1],side=3, at=label_at, cex=1.4, las=1, line=-0.6)


	plot(Vm_bias~ mat_ratio, dd, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~italic(hat(V)[M])), 
		pch=pch, 
		col=col,
		bg=bg,
		ylim=Vm_lim)
	# arrows(dd$mat_ratio,dd$Va_bias+dd_se$Va_bias,dd$mat_ratio,dd$Va_bias-dd_se$Va_bias,code=3,angle=90,length=0.1)
	arrows(dd$mat_ratio,dd$Vm_bias+dd_se$Vm_bias,dd$mat_ratio,dd$Vm_bias-dd_se$Vm_bias,code=3,angle=90,length=0.01, col=col)
	abline(h=0)

	if(lines){
		coefsM <- sapply(s,function(i){
			coefs <- coef(lm(Vm_bias~mat_ratio,dd,subset=scenario==i))
			preds <- data.frame(x=range(dd$mat_ratio))
			preds$y <- coefs[1] + coefs[2]*preds$x
			print(preds)
			lines(y~x,preds,col=cols[i], lty=2)
		})
	}
	mtext(labels[2],side=3, at=label_at, cex=1.4, las=1, line=-0.6)

}

plot_fig4 <-function(s, Va_lim=c(-0.11,0.35), Vm_lim=c(-0.17,0.06), legend_pos="bottomleft"){
	dd<-subset(va2, scenario %in% s)
	dd_se<-subset(va2_se, scenario %in% s)
	
cols <- c(palette.colors(),1)[dd$scenario]

	plot(Vm_bias~ Va_bias, dd, cex=1, xlab=expression(Bias~"in"~italic(hat(V)[A])), ylab=expression(Bias~"in"~italic(hat(V)[M])), 
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
pdf(paste0(wd,"Figures/Fig4_bias.pdf"), height=10, width=8)
{	
	par(mfrow=c(3,2),	mar=c(5,5,1,1), cex.lab=1.5, cex.axis=1.1)

	plot_fig3(c(1:3), labels=c("A)","B)"))
	plot_fig3(c(4:5), labels=c("C)","D)"))
	plot_fig3(c(6:10), labels=c("E)","F)"))
}
dev.off()

setEPS()
pdf(paste0(wd,"Figures/Fig5_Va_Vm.pdf"), height=5, width=10)
{
	par(mfrow=c(1,2),	mar=c(5,5,1,1), cex.lab=1.5, cex.axis=1.1)
	plot_fig4(c(1:6))
	legend("topright",expression(m^2~"="~"-"*0.5*h^2),lty=1,bty="n")
	mtext("A)",side=3, at=-0.225, cex=1.4, las=1, line=-0.5)

	plot_fig4(c(6:10), legend_pos="bottomright")

	mtext("B)",side=3, at=-0.225, cex=1.4, las=1, line=-0.5)

}
dev.off()



setEPS()
pdf(paste0(wd,"Figures/FigSM_extra_scenarios.pdf"), height=4, width=10)
{
	par(mfrow=c(1,2),	mar=c(5,5,1,1), cex.lab=1.5, cex.axis=0.9)
	plot_fig3(c(11:12), 
		cols=c(palette.colors(),palette.colors()),
		bgs=c(palette.colors(),palette.colors()), 
		pchs=rep(21:25,3),
		Va_lim=c(-0.1,0.1), Vm_lim=c(-0.05,0.05),
		labels=c("A)","B)"), label_at=0.075
		)
}
dev.off()



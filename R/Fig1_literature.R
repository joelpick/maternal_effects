
rm(list=ls())

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

dd<-read.csv(paste0(data_wd,"MGE - all_est.csv"))

setEPS()
pdf("/Users/joelpick/github/maternal_effects/Figures/h2_diff.pdf", height=5, width=5)

{
	par(mar=c(4,4,1,1))
hist(round(dd$h2_1-dd$h2_2,3)[dd$h2_1-dd$h2_2>=0], breaks=15, xlim=c(-0.1,0.2), main="",xlab=expression(Difference~between~h^2~from~simple~and~full~models))
hist(round(dd$h2_1-dd$h2_2,3)[dd$h2_1-dd$h2_2<0], breaks=5,col="red", add=TRUE)
}
dev.off()

dd$h2_diff <- (dd$h2_1-dd$h2_2)/dd$h2_2

str(dd)

## genetic correlations
hist(dd$r_2, breaks=20)

hist(dd$h2_diff, breaks=5000, xlim=c(-1,2))
mean(dd$h2_diff, na.rm=TRUE)



dd$h2_1/dd$h2_2 -1 
## massive one (dd[32,]) might be a typo? Otherwise two on the <0 ones have a negative covariance, and many of the 0s have no maternal variance

plot(dd$h2_1-dd$h2_2 ~ dd$m2_2)


setEPS()
pdf("/Users/joelpick/github/maternal_effects/Figures/lit_Vmg.pdf", height=5, width=10)

{
par(mfrow=c(1,2),mar=c(5,5,1,1))
# hist(dd$m2_2, breaks=20, xlab="Mg^2", main="") 
# hist(dd$m2_2[dd$juv_trait==1], breaks=20, xlab="Mg^2", main="") 


hist(dd$m2_2, breaks=20, xlab=expression(m[g]^2), main="", col="blue") 
legend("topright",c("1st year traits","Older traits"),col=c("lightblue","blue"), pch=15, bty="n", cex=1)
hist(dd$m2_2[dd$juv_trait==1], breaks=20, main="", add=TRUE, col="lightblue") 
text(0.4,16,paste("1st year mean =",round(mean(dd$m2_2[dd$juv_trait==1],na.rm=TRUE),3)), adj=0)
text(0.4,15,paste("All mean =",round(mean(dd$m2_2,na.rm=TRUE),3)), adj=0)

# x<-"a"
# expression(bquote(.(x)^2))

## proportion of variance due to Mge
hist((dd$m2_2/(dd$m2_2+dd$c2_2))[(dd$m2_2+dd$c2_2)>0.05], breaks=20, xlab="", main="", col="blue")
mtext(expression(frac(V[mg] , V[mg]+V[me])), side=1, line=4, cex=1)
mtext(expression( (m^2>0.05)), side=1, line=3.5, cex=1, adj=1)

# hist((dd$m2_2/(dd$m2_2+dd$c2_2)), breaks=20, xlab="Vmg / Vmg+Vme", main="") 
hist((dd$m2_2/(dd$m2_2+dd$c2_2))[dd$juv_trait==1 & (dd$m2_2+dd$c2_2)>0.05], breaks=20,add=TRUE, col="lightblue") 
text(0.1,4.7,paste("1st year mean =",round(mean((dd$m2_2/(dd$m2_2+dd$c2_2))[dd$juv_trait==1 &(dd$m2_2+dd$c2_2)>0.05],na.rm=TRUE),3)), adj=0)
text(0.1,4.4,paste("All mean =",round(mean((dd$m2_2/(dd$m2_2+dd$c2_2))[(dd$m2_2+dd$c2_2)>0.05],na.rm=TRUE),3)), adj=0)




}
dev.off()





cbind(dd$m2_2,dd$c2_2,dd$m2_2+dd$c2_2)
plot(dd$m2_2/(dd$m2_2+dd$c2_2),dd$m2_2+dd$c2_2)

setEPS()
pdf("/Users/joelpick/github/maternal_effects/Figures/fig1_lit.pdf", height=5, width=15)

{
par(mfrow=c(1,3),mar=c(5,5,1,1))
# hist(dd$m2_2, breaks=20, xlab="Mg^2", main="") 
# hist(dd$m2_2[dd$juv_trait==1], breaks=20, xlab="Mg^2", main="") 


hist(dd$m2_2, breaks=20, xlab=expression(m[g]^2), main="", col="blue") 
legend("topright",c("1st year traits","Older traits"),col=c("lightblue","blue"), pch=15, bty="n", cex=1)
hist(dd$m2_2[dd$juv_trait==1], breaks=20, main="", add=TRUE, col="lightblue") 
text(0.4,16,paste("1st year mean =",round(mean(dd$m2_2[dd$juv_trait==1],na.rm=TRUE),3)), adj=0)
text(0.4,15,paste("All mean =",round(mean(dd$m2_2,na.rm=TRUE),3)), adj=0)

# x<-"a"
# expression(bquote(.(x)^2))

## proportion of variance due to Mge
hist((dd$m2_2/(dd$m2_2+dd$c2_2))[(dd$m2_2+dd$c2_2)>0.05], breaks=20, xlab="", main="", col="blue")
mtext(expression(frac(V[mg] , V[mg]+V[me])), side=1, line=4, cex=1)
mtext(expression( (m^2>0.05)), side=1, line=3.5, cex=1, adj=1)

# hist((dd$m2_2/(dd$m2_2+dd$c2_2)), breaks=20, xlab="Vmg / Vmg+Vme", main="") 
hist((dd$m2_2/(dd$m2_2+dd$c2_2))[dd$juv_trait==1 & (dd$m2_2+dd$c2_2)>0.05], breaks=20,add=TRUE, col="lightblue") 
text(0.1,4.7,paste("1st year mean =",round(mean((dd$m2_2/(dd$m2_2+dd$c2_2))[dd$juv_trait==1 &(dd$m2_2+dd$c2_2)>0.05],na.rm=TRUE),3)), adj=0)
text(0.1,4.4,paste("All mean =",round(mean((dd$m2_2/(dd$m2_2+dd$c2_2))[(dd$m2_2+dd$c2_2)>0.05],na.rm=TRUE),3)), adj=0)


hist(round(dd$h2_1-dd$h2_2,3)[dd$h2_1-dd$h2_2>=0], breaks=15, xlim=c(-0.1,0.2), main="",xlab=expression(Difference~between~h^2~from~simple~and~full~models))
hist(round(dd$h2_1-dd$h2_2,3)[dd$h2_1-dd$h2_2<0], breaks=5,col="red", add=TRUE)

}
dev.off()

plot((dd$h2_1-dd$h2_2),(dd$c2_1-(dd$m2_2+dd$c2_2)),pch=19,col=c(1,2)[as.factor(dd$h2_1-dd$h2_2>=0)])
abline(h=0,v=0)
abline(0,-0.5)

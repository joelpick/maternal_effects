dd<-read.csv("~/github/maternal_effects/Data/Intermediate/MGE - all_est.csv")

setEPS()
pdf("/Users/joelpick/github/maternal_effects/Figures/h2_diff.pdf", height=5, width=5)

{
	par(mar=c(4,4,1,1))
hist(round(dd$h2_1-dd$h2_2,3)[dd$h2_1-dd$h2_2>=0], breaks=15, xlim=c(-0.1,0.2), main="",xlab="Difference between h^2 in simple and full model")
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


hist(dd$m2_2, breaks=20, xlab="Mg^2", main="", col="blue") 
legend("topright",c("1st year traits","Older traits"),col=c("lightblue","blue"), pch=15, bty="n", cex=1)
hist(dd$m2_2[dd$juv_trait==1], breaks=20, xlab="Mg^2", main="", add=TRUE, col="lightblue") 
mean(dd$m2_2,na.rm=TRUE)
mean(dd$m2_2[dd$juv_trait==1],na.rm=TRUE)

#3 proportion of variance due to Mge
hist((dd$m2_2/(dd$m2_2+dd$c2_2))[(dd$m2_2+dd$c2_2)>0.05], breaks=20, xlab="Vmg / Vmg+Vme (m^2>0.05)", main="", col="blue") 
# hist((dd$m2_2/(dd$m2_2+dd$c2_2)), breaks=20, xlab="Vmg / Vmg+Vme", main="") 
hist((dd$m2_2/(dd$m2_2+dd$c2_2))[dd$juv_trait==1 & (dd$m2_2+dd$c2_2)>0.05], breaks=20,add=TRUE, col="lightblue") 


}
dev.off()

rm(list=ls())

library(xtable)

wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

dd<-read.csv(paste0(data_wd,"MGE - all_est.csv"))



########
## Figure 2
########

setEPS()
pdf("/Users/joelpick/github/maternal_effects/Figures/fig1_lit.pdf", height=5, width=15)

{
par(mfrow=c(1,3),mar=c(6,6,1,1), cex.lab=2, cex.axis=1.6, mgp=c(4,1,0))
# hist(dd$m2_2, breaks=20, xlab="Mg^2", main="") 
# hist(dd$m2_2[dd$juv_trait==1], breaks=20, xlab="Mg^2", main="") 


hist(dd$m2_2, breaks=20, xlab=expression(m[g]^2), main="", col="blue") 
legend("topright",c("1st year traits","Older traits"),col=c("lightblue","blue"), pch=15, bty="n", cex=2)
hist(dd$m2_2[dd$juv_trait==1], breaks=20, main="", add=TRUE, col="lightblue") 
text(0.35,17,paste("1st year mean =",round(mean(dd$m2_2[dd$juv_trait==1],na.rm=TRUE),3)), adj=0, cex=1.75)
text(0.35,16,paste("All mean =",round(mean(dd$m2_2,na.rm=TRUE),3)), adj=0, cex=1.75)

# x<-"a"
# expression(bquote(.(x)^2))

## proportion of variance due to Mge
hist((dd$m2_2/(dd$m2_2+dd$c2_2))[(dd$m2_2+dd$c2_2)>0.05], breaks=20, xlab="", main="", col="blue")
mtext(expression(frac(V[mg] , V[mg]+V[me])), side=1, line=5, cex=1)
mtext(expression( (m^2>0.05)), side=1, line=3.5, cex=1, adj=1)

# hist((dd$m2_2/(dd$m2_2+dd$c2_2)), breaks=20, xlab="Vmg / Vmg+Vme", main="") 
hist((dd$m2_2/(dd$m2_2+dd$c2_2))[dd$juv_trait==1 & (dd$m2_2+dd$c2_2)>0.05], breaks=20,add=TRUE, col="lightblue") 
text(0.05,4.7,paste("1st year mean =",round(mean((dd$m2_2/(dd$m2_2+dd$c2_2))[dd$juv_trait==1 &(dd$m2_2+dd$c2_2)>0.05],na.rm=TRUE),3)), adj=0, cex=1.75)
text(0.05,4.4,paste("All mean =",round(mean((dd$m2_2/(dd$m2_2+dd$c2_2))[(dd$m2_2+dd$c2_2)>0.05],na.rm=TRUE),3)), adj=0, cex=1.75)


hist(round(dd$h2_1-dd$h2_2,3)[dd$h2_1-dd$h2_2>=0], breaks=15, xlim=c(-0.1,0.2), main="",xlab=expression(h^2~difference))
hist(round(dd$h2_1-dd$h2_2,3)[dd$h2_1-dd$h2_2<0], breaks=5,col="red", add=TRUE)

}
dev.off()



########
## Figure S3
########


setEPS()
pdf("/Users/joelpick/github/maternal_effects/Figures/FigS2_lit_excluded.pdf", height=6, width=6)
par(mar=c(6,5,1,1))
plot(dd$m2_2/(dd$m2_2+dd$c2_2),dd$m2_2+dd$c2_2, pch=19, col=scales::alpha(1,0.3), xlab="", ylab=expression(m^2), cex.lab=1.5)
mtext(expression(frac(V[mg] , V[mg]+V[me])), side=1, line=5, cex=1.5)

abline(h=0.05, col=scales::alpha(2,0.75))
dev.off()



########
## TABLE S1
########

dd$Age <- ifelse(dd$juv_trait==1,"Juvenile","Adult")
xt <- xtable(dd[,c(names(dd)[1:4],"Age",names(dd)[6:11])],digits=3)
names(xt)[6:10] <- c('$h^2_1$','$c^2_1$','$h^2_2$','$m^2_2$','$c^2_2$' )
print(xt, include.rownames=FALSE,sanitize.text.function=function(x){x},only.contents=TRUE)







subset(dd,m2_2>0.6)


# prop m2 of total m2
dd$h2m_2 <-dd$m2_2/(dd$m2_2+dd$c2_2)

(dd$m2_2/(dd$m2_2+dd$c2_2))[(dd$m2_2+dd$c2_2)>0.05]


cbind(dd$m2_2,dd$c2_2,dd$m2_2+dd$c2_2)

subset(dd,m2_2+c2_2<0.05)
subset(dd,h2m_2<0.05)
subset(dd,m2_2==0)


plot((dd$h2_1-dd$h2_2),(dd$c2_1-(dd$m2_2+dd$c2_2)),pch=19,col=c(1,2)[as.factor(dd$h2_1-dd$h2_2>=0)])
abline(h=0,v=0)
abline(0,-0.5)

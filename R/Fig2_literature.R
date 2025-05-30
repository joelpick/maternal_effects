
rm(list=ls())

library(xtable)

wd <- "/Users/joelpick/github/maternal_effects/"

dd<-read.csv(paste0(wd,"Data/Raw/MGE - all_est.csv"))

## n studies
length(unique(dd$Study))

##n estimates
nrow(dd)

##adult and juvenile estimates
table(dd$juv_trait)


########
## Figure 2
########

setEPS()
pdf(paste0(wd,"Figures/Fig2_lit.pdf"), height=10, width=12)

{
par(mfrow=c(2,2),mar=c(6,6,1,1), cex.lab=2, cex.axis=1.6, mgp=c(4,1,0))

hist(dd$m2_2, breaks=20, xlab=expression(italic(hat(m)[g]^2)), main="", col="blue") 
legend("topright",c("1st year traits","Older traits"),col=c("lightblue","blue"), pch=15, bty="n", cex=1.6)
hist(dd$m2_2[dd$juv_trait==1], breaks=20, main="", add=TRUE, col="lightblue") 
text(0.65,16.5,paste("1st year mean =",round(mean(dd$m2_2[dd$juv_trait==1],na.rm=TRUE),3)), adj=0, cex=1.75, pos=2)
text(0.65,15,paste("All mean =",round(mean(dd$m2_2,na.rm=TRUE),3)), adj=0, cex=1.75, pos=2)
mtext("A)",side=3, adj = -0.175, cex=1.7, las=1, line=-0.7)

# x<-"a"
# expression(bquote(.(x)^2))

## proportion of variance due to Mge
hist((dd$m2_2/(dd$m2_2+dd$c2_2))[(dd$m2_2+dd$c2_2)>0.05], breaks=20, xlab="", main="", col="blue")
mtext(expression(italic(frac(hat(V)[Mg] , hat(V)[Mg]+hat(V)[Me]))), side=1, line=5, cex=1)
mtext(expression( italic(hat(m)^2>0.05)), side=1, line=3.5, cex=1, adj=1)

# hist((dd$m2_2/(dd$m2_2+dd$c2_2)), breaks=20, xlab="Vmg / Vmg+Vme", main="") 
hist((dd$m2_2/(dd$m2_2+dd$c2_2))[dd$juv_trait==1 & (dd$m2_2+dd$c2_2)>0.05], breaks=20,add=TRUE, col="lightblue") 
text(0.5,4.65,paste("1st year mean =",round(mean((dd$m2_2/(dd$m2_2+dd$c2_2))[dd$juv_trait==1 &(dd$m2_2+dd$c2_2)>0.05],na.rm=TRUE),3)), adj=0, cex=1.75, pos=2)
text(0.5,4.2,paste("All mean =",round(mean((dd$m2_2/(dd$m2_2+dd$c2_2))[(dd$m2_2+dd$c2_2)>0.05],na.rm=TRUE),3)), adj=0, cex=1.75, pos=2)
mtext("B)",side=3, adj = -0.175, cex=1.7, las=1, line=-0.7)


# hist(round(dd$h2_1-dd$h2_2,3)[dd$h2_1-dd$h2_2>=0], breaks=15, xlim=c(-0.1,0.2), main="",xlab=expression(h^2~difference))
# hist(round(dd$h2_1-dd$h2_2,3)[dd$h2_1-dd$h2_2<0], breaks=5,col="red", add=TRUE)

h2_diff <- dd$h2_1-dd$h2_2
hist(h2_diff[h2_diff>=0], breaks=10, xlim=c(-0.1,0.2), main="",xlab=expression(italic(hat(h)^2)~difference),col=scales::alpha("red",0.8))
hist(h2_diff[h2_diff<0], breaks=2, add=TRUE)
mtext("C)",side=3, adj = -0.175, cex=1.7, las=1, line=-0.7)


m2_diff <- dd$c2_1 - (dd$m2_2+dd$c2_2)
hist(m2_diff[m2_diff<0], breaks=15, col=scales::alpha("red",0.8), xlim=c(-0.15,0.05),xlab=expression(italic(hat(m)^2)~difference), main="")
hist(m2_diff[m2_diff>0], breaks=2, add=TRUE)
mtext("D)",side=3, adj = -0.175, cex=1.7, las=1, line=-0.7)

}
dev.off()

## number with simple and complex estimates:
sum!is.na(h2_diff))

########
## Figure S3
########


setEPS()
pdf(paste0(wd,"Figures/FigSM_lit_excluded.pdf"), height=6, width=6)
par(mar=c(6,5,1,1))
plot(dd$m2_2/(dd$m2_2+dd$c2_2),dd$m2_2+dd$c2_2, pch=19, col=scales::alpha(1,0.3), xlab="", ylab=expression(italic(hat(m)^2)), cex.lab=1.5)
mtext(expression(frac(italic(hat(V)[Mg]) , italic(hat(V)[Mg])+italic(hat(V)[Me]))), side=1, line=5, cex=1.5)

abline(h=0.05, col=scales::alpha(2,0.75))
dev.off()



########
## TABLE S1
########

dd$Age <- ifelse(dd$juv_trait==1,"Juvenile","Adult")
xt <- xtable(dd[,c(names(dd)[1:4],"Age",names(dd)[6:11])],digits=3)
names(xt)[6:10] <- c('$h^2_1$','$c^2_1$','$h^2_2$','$m^2_2$','$c^2_2$' )
print(xt, include.rownames=FALSE,sanitize.text.function=function(x){x},only.contents=TRUE)



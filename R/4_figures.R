rm(list=ls())


wd <- "/Users/joelpick/github/maternal_effects/"

data_wd <- paste0(wd,"Data/Intermediate/")

load(paste0(data_wd,"mge_sims_rd.Rdata"))
load(paste0(data_wd,"mge_sims_known_rd.Rdata"))


ped_sum<-cbind(sapply(ped_str,colMeans), ped_str_known)

ped_sum2 <- rbind(mat_sibs=colSums(ped_sum[c("FS","MHS"),]), mat_links=colSums(ped_sum[c("dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS"),]), other=colSums(ped_sum[!rownames(ped_sum)%in%c( "individuals" ,"links" ,"FS","MHS","dam","MG","au_D_FS","au_D_MHS","cousin_D_FS","cousin_D_HS"),]) )



ps <- t(t(ped_sum2)/colSums(ped_sum2))
# ps <- t(t(ped_sum2[-3,])/colSums(ped_sum2[-3,]))
# barplot(ps)
par(mfrow=c(3,1), mar=c(4,4,1,1), cex.lab=1.4,mgp=c(2,0.5,0))
barplot(ps[,c(5,1,4)], names=c("Low", "Mid", "High"), xlab="Fecundity", col=c("lightblue","orange","white"))
barplot(ps[,c(1,8,6,7)], names=c("None", "Both sexes", "Female", "Male"), xlab="Dispersal", ylab="Prop. relationships",col=c("lightblue","orange","white"))
barplot(ps[,c(3,1,2)], names=c("Half-sibs", "Mixed", "Full- sibs"), xlab="Mating system",col=c("lightblue","orange","white"))

ps

mat_ratio<-ped_sum2[2,]/ped_sum2[1,]
mat_ratio2<-ped_sum2[2,]/colSums(ped_sum2)
mat_ratio3<-ped_sum2[1,]/colSums(ped_sum2)
mat_ratio4<-ped_sum2[1,]/ped_sum2[3,]
mat_ratio5<-ped_sum2[2,]/ped_sum2[3,]
mat_ratio6<-ped_sum2[3,]/colSums(ped_sum2)


ped_names <- colnames(ped_sum2)

mod2<-do.call(rbind,lapply(ped_names,function(k) {
	mod2 <- do.call(rbind,lapply(get(paste0("model2_",k)), function(x) {
		do.call(rbind,lapply(1:nrow(scenarios), function(i) data.frame(r=k,scenario=i,Va_est = x[i,"A"],Vm_est = x[i,"Me"],Va_sim=scenarios[i,"Va"],Vm_sim =sum(scenarios[i,c("Vmg","Vme")]),Vmg_sim=scenarios[i,"Vmg"])))
	}))
	# assign(paste0("mod2_",k),mod2)
}))
#,sum(x[i,c("Mg","Me")])
head(mod2,20)
mod2$Va_bias <- (mod2$Va_est - mod2$Va_sim) 
# mod2$ln_Va_bias <- log(mod2$Va_bias)
va2<-aggregate(cbind(Va_bias,Vmg_sim,Vm_sim)~ scenario+r, mod2,mean)
# which(va1$r)

r_order<- sapply(va2$r, function(x) which(ped_names==x))

va2$mat_ratio <- mat_ratio2[r_order]
va2$known_ped <- ifelse(va2$r %in% c("bt","rd"), "2","1")
par(mfrow=c(1,1), mar=c(5,5,1,1), cex.lab=1.75, cex.axis=1.25 )
plot(Va_bias~ mat_ratio, va2, subset=scenario==1,  pch=19, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression("%"~Bias~"in"~h^2), col=va2$known_ped)

# plot(Va_bias~ mat_ratio, va2,  pch=19, cex=1, xlab="Proportion non-sibling maternal links", ylab=expression(Bias~"in"~h^2))



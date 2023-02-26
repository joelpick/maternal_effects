dd<-read.csv("~/github/maternal_effects/Data/Intermediate/MGE - all_est.csv")


dd$h2_diff <- (dd$h2_1-dd$h2_2)/dd$h2_2

str(dd)

## genetic correlations
hist(dd$r_2, breaks=20)

hist(dd$h2_diff, breaks=5000, xlim=c(-1,2))
mean(dd$h2_diff, na.rm=TRUE)

dd$h2_1/dd$h2_2 -1 
## massive one (dd[32,]) might be a typo? Otherwise two on the <0 ones have a negative covariance, and many of the 0s have no maternal variance

plot(dd$h2_1-dd$h2_2 ~ dd$m2_2)

#3 proportion of variance due to Mge
hist(dd$m2_2/(dd$m2_2+dd$c2_2), breaks=20) 


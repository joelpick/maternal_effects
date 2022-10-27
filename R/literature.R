dd<-read.csv("~/Downloads/MGE - all_est.csv")

str(dd)

hist(dd$r_2, breaks=10)

hist(dd$h2_1-dd$h2_2, breaks=10)

plot(dd$h2_1-dd$h2_2 ~ dd$m2_2)

hist(dd$m2_2/(dd$m2_2+dd$c2_2), breaks=10) 
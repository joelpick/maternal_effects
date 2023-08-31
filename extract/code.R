x<-metaDigitise::getExtracted("/Users/joelpick/github/maternal_effects/extract", summary=FALSE)[[1]][[1]]

round(diff(x[x$group==1,"y"]),3)
round(diff(x[x$group==2,"y"]),3)

0.134 0.370
0.041 0.429 0.035



focal <- c(x=1,y=1)

m_focal <- c(x=1,y=2)

relative <- c(x=2,y=1)

m_relative <- c(x=2,y=2)

{

par(mar=c(0,0,0,0))
plot(NA,xlim=c(0.5,2.5),ylim=c(0.5,2.5),xaxt="n",yaxt="n",xlab="",ylab="",bty="n")

text(focal["x"],focal["y"],"Focal \n Individual")
text(m_focal["x"],m_focal["y"],"Focal's \n Mother")
text(relative["x"],relative["y"],"Relative")
text(m_relative["x"],m_relative["y"],"Relative's \n Mother")

arrows(focal["x"]+0.1,focal["y"],relative["x"]-0.1,relative["y"],code=3)
text(mean(c(focal["x"],relative["x"])),mean(c(focal["y"],relative["y"])),"r1",adj=c(0.5,1))

arrows(m_focal["x"]+0.1,m_focal["y"],m_relative["x"]-0.1,m_relative["y"],code=3)
text(mean(c(m_focal["x"],m_relative["x"])),mean(c(m_focal["y"],m_relative["y"])),"r2",adj=c(0.5,0))

arrows(m_focal["x"]+0.1,m_focal["y"]-0.1,relative["x"]-0.1,relative["y"]+0.1,code=3)
text(mean(c(m_focal["x"],relative["x"])),mean(c(m_focal["y"],relative["y"])),"r3",adj=c(0,0))


arrows(focal["x"]+0.1,focal["y"]+0.1,m_relative["x"]-0.1,m_relative["y"]-0.1,code=3)
text(mean(c(focal["x"],m_relative["x"])),mean(c(focal["y"],m_relative["y"])),"r4",adj=c(1,1))

}


arrow_colour <- function(x) if(!is.numeric(x) || x>0) 1 else "grey"
arrow_lwd <- function(x) if(!is.numeric(x) || x==0) 1 else x*2+1

rel_cov_plot <- function(relative_name="Relative",relative_mother_name="Relative's \n Mother", r1="r1", r2="r2",r3="r3" ,r4="r4", arrow_length=0.1, offset=0.1,offset_diag=0.15){
	focal <- c(x=1,y=1)
	m_focal <- c(x=1,y=2)
	relative <- c(x=2,y=1)
	m_relative <- c(x=2,y=2)

	par(mar=c(0,0,0,0))
	plot(NA,xlim=c(0.75,2.25),ylim=c(0.4,2.4),xaxt="n",yaxt="n",xlab="",ylab="",bty="n")

	text(focal["x"],focal["y"],"Focal")
	text(m_focal["x"],m_focal["y"],"Mother")
	text(relative["x"],relative["y"],relative_name)
	text(m_relative["x"],m_relative["y"],relative_mother_name)

	arrows(focal["x"]+offset,focal["y"],relative["x"]-offset,relative["y"],code=3, length=arrow_length, col=arrow_colour(r1), lwd=arrow_lwd(r1))
	text(mean(c(focal["x"],relative["x"])),mean(c(focal["y"],relative["y"]))-offset,r1)


	arrows(m_focal["x"]+offset,m_focal["y"],m_relative["x"]-offset,m_relative["y"],code=3, length=arrow_length, col=arrow_colour(r2), lwd=arrow_lwd(r2))
	text(mean(c(m_focal["x"],m_relative["x"])),mean(c(m_focal["y"],m_relative["y"]))+offset,r2)

	arrows(focal["x"]+offset_diag,focal["y"]+offset_diag,m_relative["x"]-offset_diag,m_relative["y"]-offset_diag,code=3, length=arrow_length,, col=arrow_colour(r3), lwd=arrow_lwd(r3))
	text(1.6+offset,1.575,r3)

	arrows(m_focal["x"]+offset_diag,m_focal["y"]-offset_diag,relative["x"]-offset_diag,relative["y"]+offset_diag,code=3, length=arrow_length, col=arrow_colour(r4), lwd=arrow_lwd(r4))
	text(1.4-offset,1.575,r4)
	
	if(is.numeric(r1)&is.numeric(r2)&is.numeric(r3)&is.numeric(r4)){
		print_r2 <- if(r2>0 & r2<1) paste0(" + ",r2,"~V[Mg]") else if(r2==1) " + ~V[Mg]" else ""
		print_r_34 <- if((r3+r4)==1) " + ~COV['A,Mg']" else if((r3+r4)>0) paste0(" + ",(r3+r4),"~COV['A,Mg']") else ""
		print_r <- paste0(r1,"~V[A]",print_r2,print_r_34)
		text(1.5,0.5, parse(text=print_r), col="blue")
	}else{
		text(1.5,0.5, expression(r[1]*V[A] + r[2]*V[Mg] + (r[3]+r[4])*COV["A,Mg"]), col="blue")
	}

}


dat <- read.csv("/Users/joelpick/github/maternal_effects/Data/Raw/covariances\ -\ Sheet5.csv")
head(dat)

relative2 <- sub(" "," \n ",dat$relative)
relative_mother2 <- sub(" ","\n",dat$relative_mother)



setEPS()
pdf(paste0(wd,"Figures/FigS1_relative_cov.pdf"), height=10, width=7)

{	
par(mfrow=c(7,3))

for(i in 1:nrow(dat)) {
	rel_cov_plot(
		relative_name=relative2[i],#dat[i,"relative"],
		relative_mother_name=relative_mother2[i],#dat[i,"relative_mother"], 
		r1=dat[i,"r1"], r2=dat[i,"r2"],r3=dat[i,"r3.f.rm"] ,r4=dat[i,"r4.r.m"],
		offset=0.2,offset_diag=0.175)
	text(0.75,2.2, paste0(LETTERS[i],")"), cex=1.4 )
}
}
dev.off()


dd <- paste(expression(frac(V[mg] , V[mg]+V[me])),1)
expression(eval(dd))

setEPS()
pdf(paste0(wd,"Figures/FigS1_relative_cov_eg.pdf"), height=4, width=4)

{	
par(mfrow=c(1,1))
rel_cov_plot(r1=expression(r[1]),r2=expression(r[2]),r3=expression(r[3]),r4=expression(r[4]),offset=0.15,offset_diag=0.125)
}

dev.off()
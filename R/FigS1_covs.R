


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



rel_cov_plot <- function(relative_name="Relative",relative_mother_name="Relative's \n Mother", r1="r1", r2="r2",r3="r3" ,r4="r4", arrow_length=0.1, offset=0.1,offset_diag=0.1){
	focal <- c(x=1,y=1)
	m_focal <- c(x=1,y=2)
	relative <- c(x=2,y=1)
	m_relative <- c(x=2,y=2)

	par(mar=c(0,0,0,0))
	plot(NA,xlim=c(0.5,2.5),ylim=c(0.5,2.5),xaxt="n",yaxt="n",xlab="",ylab="",bty="n")

	text(focal["x"],focal["y"],"Focal")
	text(m_focal["x"],m_focal["y"],"Mother")
	text(relative["x"],relative["y"],relative_name)
	text(m_relative["x"],m_relative["y"],relative_mother_name)

	arrows(focal["x"]+offset,focal["y"],relative["x"]-offset,relative["y"],code=3, length=arrow_length)
	text(mean(c(focal["x"],relative["x"])),mean(c(focal["y"],relative["y"]))-offset,r1)

	arrows(m_focal["x"]+offset,m_focal["y"],m_relative["x"]-offset,m_relative["y"],code=3, length=arrow_length)
	text(mean(c(m_focal["x"],m_relative["x"])),mean(c(m_focal["y"],m_relative["y"]))+offset,r2)

	arrows(focal["x"]+offset_diag,focal["y"]+offset_diag,m_relative["x"]-offset_diag,m_relative["y"]-offset_diag,code=3, length=arrow_length)
	text(1.65+offset,1.65,r3)

	arrows(m_focal["x"]+offset_diag,m_focal["y"]-offset_diag,relative["x"]-offset_diag,relative["y"]+offset_diag,code=3, length=arrow_length)
	text(1.35-offset,1.65,r4)
	
	print_r2 <- if(r2>0 & r2<1) paste0(" + ",r2,"Vmg") else if(r2==1) " + Vmg" else ""
	print_r_34 <- if((r3+r4)==1) " + COVa,mg" else if((r3+r4)>0) paste0(" + ",(r3+r4),"COVa,mg") else ""
	text(1.5,0.5, paste0(r1,"Va",print_r2,print_r_34))

}




dat <- read.csv("/Users/joelpick/github/maternal_effects/Data/Raw/covariances\ -\ Sheet5.csv")
head(dat)
par(mfrow=c(7,3))

for(i in 1:nrow(dat)) rel_cov_plot(
	relative_name=dat[i,"relative"],
	relative_mother_name=dat[i,"relative_mother"], 
	r1=dat[i,"r1"], r2=dat[i,"r2"],r3=dat[i,"r3.f.rm"] ,r4=dat[i,"r4.r.m"],
	offset=0.2,offset_diag=0.15)

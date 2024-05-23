


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
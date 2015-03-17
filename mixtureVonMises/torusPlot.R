install.packages('NPCirc');
library(NPCirc)

### circular-circular
data(wind)
wind6 <- circular(wind$wind.dir[seq(7,1752,by=24)])
wind12 <- circular(wind$wind.dir[seq(13,1752,by=24)])
estNW <- kern.reg.circ.circ(wind6,wind12,t=NULL,bw=6.1,method="NW")
estLL <- kern.reg.circ.circ(wind6,wind12,t=NULL,bw=2.25,method="LL")
# Torus representation
plot(estNW, plot.type="circle", points.plot=TRUE, line.col=2, lwd=2, points.col=2,
units="degrees")
lines(estLL, plot.type="circle", line.col=3, lwd=2)
# Linear representation
plot(estNW, plot.type="line", points.plot=TRUE, xlab="Wind direction at 6 a.m.",
ylab="Wind direction at noon")
lines(estLL, plot.type="line", line.col=2)
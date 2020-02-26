library(tidyverse)

irs4oe <- tribble(
  ~time,~v1,~v2,~v3,~oe1,~oe2,~oe3,~y1,~y2,~y3,
  0,5,5,5,5,5,5,5,5,5,
  24,6.57,6.31,6.27,7.44,7.6,7.21,7.77,7.61,7.98,
  48,9.48,11.1,10.6,14.3,16,15.4,16,15.2,17.6,
  72,15.4,16.5,17.3,29.5,30.2,29.4,33.1,32.3,33.5,
  96,23,24.3,23.1,73.3,75.9,75.6,76.8,80.9,83.9
)

dr <- tribble(
  ~conc,~v1,~v2,~v3,~oe1,~oe2,~oe3,~y1,~y2,~y3,
  -12,	1.009149,	1.011366,	0.9794852,	1.004869,	0.9992209,	0.9959098,	0.9452502,	0.9562835,	1.098466,
  -10,	1.003481,	1.035063,	0.9929107,	1.041213,	1.16236,	0.8745683,	0.8260403,	1.056208,	1.184701,
  -9,	0.7350116,	1.06575,	1.08331,	1.1266,	1.016984,	0.8929155,	0.7212707,	1.07043,	1.184701,
  -8,	0.6863812,	0.9867733,	1.008126,	1.032488,	1.074246,	0.8593373,	0.8956783,	1.021671,	1.17198,
  -7,	0.5745013,	0.7294709,	0.7378672,	1.008063,	0.8430156,	0.7652634,	0.7442124,	0.8931779,	0.8417932,
  -6,	0.3198414,	0.4887907,	0.4705063,	0.8080349,	0.7812346,	0.6096424,	0.6291596,	0.7572774,	0.648632,
  -5,	0.2779451,	0.360587,	0.4054242,	0.6588413,	0.6138105,	0.5053237,	0.5420808,	0.6138755,	0.5688046,
  -4,	0.04178269,	0.04958231,	0.05537876,	0.05178279,	0.04391409,	0.03701924,	0.03479819,	0.04333104,	0.03529828
)

cairo_pdf("~kjyi/Projects/thymus_single_cell/final2/figures/irs4_overexpression.pdf",width=1150/2/254,height=1050/254,pointsize = 12*0.7)

par(mfrow=c(2,1),mar=c(3,4,.5,.5))

plot(200,xlim=c(1,5),ylim=c(0,80),xlab="",ylab="",xaxt="n",las=2)
mtext("Time (h)",side=1,line = 2)
mtext(expression("Cell number ("*"x10"^4*")"),side=2,line=2)
axis(1,at=1:5,labels=irs4oe$time)
lines(rowMeans(irs4oe[,2:4]))
points(rowMeans(irs4oe[,2:4]),pch=16)
lines(rowMeans(irs4oe[,5:7]),col="red")
points(rowMeans(irs4oe[,5:7]),pch=16,col="red")
arrows(2:5,apply((irs4oe[,2:4]),1,min)[2:5],2:5,apply((irs4oe[,2:4]),1,max)[2:5],angle=90,code = 3,length = 0.05,col="black")
arrows(2:5,apply((irs4oe[,5:7]),1,min)[2:5],2:5,apply((irs4oe[,5:7]),1,max)[2:5],angle=90,code = 3,length = 0.05,col="red")
legend("topleft",lty=1,pch=16,legend = c('Vector','IRS4 overexpression'),col=c("black","red"),bty='n')





plot(200,xlim=c(1,9),ylim=c(0,1.2),xlab="",ylab="",xaxt="n",las=2)
mtext(expression(log[10]*"[OSI-906] (M)"),side=1,line = 2)
mtext("Cell viability",side=2,line=2.5)
axis(1,at=c(1,3:9),labels=dr$conc)

lines(c(1,3:9),rowMeans(dr[,2:4]))
points(c(1,3:9),rowMeans(dr[,2:4]),pch=16)
lines(c(1,3:9),rowMeans(dr[,5:7]),col="red")
points(c(1,3:9),rowMeans(dr[,5:7]),pch=16,col="red")
arrows(3:9,apply((dr[-1,2:4]),1,min),3:9,apply((dr[-1,2:4]),1,max),angle=90,code = 3,length = 0.05,col="black")
arrows(3:9,apply((dr[-1,5:7]),1,min),3:9,apply((dr[-1,5:7]),1,max),angle=90,code = 3,length = 0.05,col="red")
# legend("topleft",lty=1,pch=16,legend = c('Vector','IRS4 overexpression'),col=c("black","red"),bty='n')
dev.off()


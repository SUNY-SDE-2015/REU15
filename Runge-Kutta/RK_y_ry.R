setwd("~/Documentos/REU/RungeKutta")
data<-read.csv('RK_y_ry.csv')
steps=data$steps[1]

y<-seq(from=0, to=1, by=0.01)
x<-seq(from=0, to=1, by=1/steps)

plot(x,data$Approximation,main="Runge Kutta",ylab="y")
lines(y,exp(y*data$r[1]),'l')


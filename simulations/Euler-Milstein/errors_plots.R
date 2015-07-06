approximation<-read.csv("eulermaruyama_milstein.csv")
alfa<-2;
beta<-1;
euler_maruyama<-approximation$ErrorEuler[(approximation$alpha==alfa) & (approximation$beta==beta)]
milstein<-approximation$ErrorMilstein[(approximation$alpha==alfa) & (approximation$beta==beta)]
time<-approximation$dt[(approximation$alpha==alfa) & (approximation$beta==beta)]
plot(time,euler_maruyama,'l',col=9,ylab="Error",xlab="dt")
lines(time,milstein,'l',col=3)
legend("topleft",c("Euler","Milstein"),col = c(9,3),lty = c(2, 2))
title("alfa=2 beta=1")


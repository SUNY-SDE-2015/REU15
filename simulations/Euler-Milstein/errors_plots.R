approximation<-read.csv("eulermaruyama_milstein.csv")
alfa<-2;
beta<-0.5;

euler_maruyama <-approximation$ErrorEuler[
                     (approximation$alpha==alfa) & (approximation$beta==beta)]
milstein      <- approximation$ErrorMilstein[
                     (approximation$alpha==alfa) & (approximation$beta==beta)]
time          <- approximation$dt[
                     (approximation$alpha==alfa) & (approximation$beta==beta)]

plot(time,euler_maruyama,'l',col=4,ylab="Error",xlab="dt",lty=2)
lines(time,milstein,'l',col=2,lty=1)

legend("topleft",c("Euler","Milstein"),col = c(4,2),lty = c(2, 1))
title(substitute(paste("Approximation Errors for ",
                       alpha,'=',a,', ',
                       beta,'=',b),list(a=alfa,b=beta)))


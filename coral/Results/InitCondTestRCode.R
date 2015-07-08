data<-read.csv("trialsA.csv",header=TRUE)

a<-0
m<-seq(0,40,by=1)
n<-seq(1,16400,by=1)
for (i in m)
	{
		a[i+1]<-0
	}
for (i in m)
	{
		for (j in n)
			{
				if (data$coral[j]>0.55 && data$theta[j]>pi/80*i-0.01 && data$theta[j]<pi/80*i+0.01)
					{
						a[i+1]<-a[i+1]+1
					}
			}
	}



phi<-seq(0,pi/2,pi/80)

b<-a/400

fit1<-glm(cbind(a,c)~phi,family=binomial(link=logit))
summary(fit1)

fit2<-glm(a~phi,family=poisson(link=log))
summary(fit2)

fit3<-lm(b~phi)
summary(fit3)

e<-seq(0,pi/2,by=0.001)
d<-exp(0.33078+3.67635*e)/(1+exp(0.33078+3.67635*e))
f<-exp(5.69056+0.24592 *e)/400
g<-0.72963+0.22166*e

plot(phi,b,main="Beta=0.5,g=0.3,Tau=0.5",xlab="Theta",ylab="Proportion of Trials to Coral Point")
lines(e,d,col='blue')
lines(e,f,col='red')
lines(e,g,col='black')
legend (0.7,0.4,c("Binomial Logit","Poisson Log","OLS"),col=c(4,2,1), lty=c(1,1,1))

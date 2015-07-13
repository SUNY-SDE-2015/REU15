data<-read.csv("trialsA.csv",header=TRUE)
a<-length(unique(data$beta))
b<-length(unique(data$theta))
c<-seq(1,a,by=1)
d<-seq(1,b,by=1)
l<-seq(1,a*b,by=1)
succ<-0					#number of sucess
total<-0					#total number of trials
prop<-0					#output values
beta<-0					#beta values
theta<-0					#theta values
for (i in l)
	{
		succ[i]<-0
		total[i]<-0
		prop[i]<-0
		beta[i]<-0
		theta[i]<-0
	}
k<-1
for(i in c)
	{
		for(j in d)
			{
				succ[k]<-length(data$coral[data$beta==unique(data$beta)[i] & data$theta==unique(data$theta)[j] & data$coral>0.55 & (!is.na(data$coral))])
				total[k]<-length(data$coral[data$beta==unique(data$beta)[i] & data$theta==unique(data$theta)[j] & (!is.na(data$coral))])
				##prop[k]<-length(data$coral[data$beta==unique(data$beta)[i] & data$theta==unique(data$theta)[j] & data$coral>0.55 & (!is.na(data$coral))])/length(data$coral[data$beta==unique(data$beta)[i] & data$theta==unique(data$theta)[j] & (!is.na(data$coral))])
				beta[k]<-unique(data$beta)[i]
				theta[k]<-unique(data$theta)[j]
				k<-k+1
			}
	}
prop<-succ/total
fail<-total-succ

theta.logit<-glm(cbind(succ[beta==0.5],total[beta==0.5]-succ[beta==0.5])~theta[beta==0.5],family=binomial(link=logit))
summary(theta.logit)

theta.poisson<-glm(succ[beta==0.5]~theta[beta==0.5],family=poisson(link=log))
summary(theta.poisson)

#library(MASS)
#theta.nb<-glm.nb(succ[beta==0.5]~theta[beta==0.5])
#summary(theta.nb)

#theta.ols<-lm(prop[beta==0.5]~theta[beta==0.5])
#summary(theta.ols)

theta.ols<-glm(prop[beta==0.5]~theta[beta==0.5],family=gaussian(link=identity))
summary(theta.ols)

#Plot of Initial Condition for beta 0.5
plot(theta[beta==0.50],prop[beta==0.50],main=expression(paste("Initial Conditions for ", beta, "=0.5")),xlab=expression(theta),ylab='Proportion of Coral Stable Point')

#Initializes Plot for the Subsequent Models
plot(theta[beta==0.50],prop[beta==0.50],main=expression(paste("Three Models for Initial Conditions for ", beta, "=0.5")),xlab=expression(theta),ylab='Proportion of Coral Stable Point')

#Plot for Logit Model
lines(seq(0,pi/2,by=pi/10000),
exp(summary(theta.logit)$coefficients[[2]]*seq(0,pi/2,by=pi/10000)+summary(theta.logit)$coefficients[[1]])
/(1+exp(summary(theta.logit)$coefficients[[2]]*seq(0,pi/2,by=pi/10000)+summary(theta.logit)$coefficients[[1]])),
col='blue')

#Plot for Poisson Model
lines(seq(0,pi/2,by=pi/10000),
exp(summary(theta.poisson)$coefficients[[2]]*seq(0,pi/2,by=pi/10000)+summary(theta.poisson)$coefficients[[1]])/400,col='green')

#Plot for OLS
lines(seq(0,pi/2,by=pi/10000),(summary(theta.ols)$coefficient[[2]]*(seq(0,pi/2,by=pi/10000))+summary(theta.ols)$coefficient[[1]]))

#Legend
legend (0.87,0.5,c("Binomial Logit","Poisson Log","OLS"),col=c(4,11,1), lty=c(1,1), lwd=5)



#TESTS FOR G VALUES
value.beta<-0.4
value.theta<-unique(theta)[(length(unique(dataVals$theta))+1)/2]
value.tau<-0.2

lgg<-gval[theta==value.theta & beta==value.beta & tau==value.tau]
lgscount<-lgsucc[theta==value.theta & beta==value.beta & tau==value.tau]
lgfcount<-lgfail[theta==value.theta & beta==value.beta & tau==value.tau]
lgsprop<-lgscount/(lgfcount+lgscount)

g.logit<-glm(cbind(lgscount,lgfcount)~lgg,family=binomial(link=logit))
g.poisson<-glm(lgscount~lgg,family=poisson(link=log))
g.ols<-glm(lgsprop~lgg)
summary(g.logit)
summary(g.poisson)
summary(g.ols)
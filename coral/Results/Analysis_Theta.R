#TESTS FOR THETA VALUES
value.beta<-0.65
value.g<-0.15
value.tau<-0.5

lgscount<-lgsucc[beta==value.beta & gval==value.g & tau==value.tau]
lgfcount<-lgfail[beta==value.beta & gval==value.g & tau==value.tau]
lgsprop<-lgscount/(lgfcount+lgscount)
lgtheta<-theta[beta==value.beta & gval==value.g & tau==value.tau]

#plot(lgtheta,lgsprop)


theta.logit<-glm(cbind(lgscount,lgfcount)~lgtheta,family=binomial(link=logit))
theta.poisson<-glm(lgscount~lgtheta,family=poisson(link=log))
theta.nb<-glm.nb(lgfcount~lgtheta)
theta.ols<-glm(lgsprop~lgtheta,family=gaussian(link=identity))
#theta.ols.log<-glm(lgsprop~log(lgtheta),family=gaussian(link=identity))
summary(theta.logit)
summary(theta.poisson)
summary(theta.ols)
#summary(theta.ols.log)
#TESTS FOR BETA VALUES
value.g<-0.15
value.theta<-unique(theta)[3]
value.tau<-0.5

lgbeta<-beta[theta==value.theta & gval==value.g & tau==value.tau]
lgscount<-lgsucc[theta==value.theta & gval==value.g & tau==value.tau]
lgfcount<-lgfail[theta==value.theta & gval==value.g & tau==value.tau]
lgsprop<-lgscount/(lgfcount+lgscount)

#plot(lgbeta,lgsprop)

beta.logit<-glm(cbind(lgscount,lgfcount)~lgbeta,family=binomial(link=logit))
beta.poisson<-glm(lgscount~lgbeta,family=poisson(link=log))
beta.ols<-glm(lgsprop~lgbeta)
summary(beta.logit)
summary(beta.poisson)
summary(beta.ols)
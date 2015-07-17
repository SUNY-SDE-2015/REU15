
#This is a code with functions in R, I changed some things because I'm not sure 
#I understood what Michael did. I get some warnings with the glm fit though. 
#So, you have to run all this code once, and then just call the function you want to use.

if(!exists('dataVals')){
  dataVals<-read.csv("trials.csv",header=TRUE)
}
  theta<-sort(unique(dataVals$theta))
  beta<-sort(unique(dataVals$beta))
  gvals<-sort(unique(dataVals$g))
  tau<-sort(unique(dataVals$tau))


test_theta<-function(b,ta,g){
  
  proportion<-rep(0,length(theta))
  succ<-rep(0,length(theta))
  total<-rep(0,length(theta))
  fail<-rep(0,length(theta))
  
  for (i in length(theta)){
    succ[i]<-length(dataVals$lgCoral[dataVals$beta==b & dataVals$theta==theta[i] & dataVals$g==g & dataVals$lgCoral>0.50 & (!is.na(dataVals$lgCoral))  & dataVals$tau==ta])
    total[i]<-length(dataVals$lgCoral[dataVals$beta==b & dataVals$theta==theta[i] & dataVals$g==g & (!is.na(dataVals$lgCoral))  & dataVals$tau==ta])
    fail[i]<-total[i]-succ[i]
    proportion[i]<-succ[i]/total[i]
    }
  
  plot(theta,proportion)
  
  
  theta.logit<-glm(cbind(succ,fail)~theta,family=binomial(link=logit))
  theta.poisson<-glm(succ~theta,family=poisson(link=log))
  theta.nb<-glm.nb(fail~theta)
  theta.ols<-glm(proportion~theta,family=gaussian(link=identity))
  #theta.ols.log<-glm(lgsprop~log(lgtheta),family=gaussian(link=identity))
  summary(theta.logit)
  summary(theta.poisson)
  summary(theta.ols)
  #summary(theta.ols.log)
}


test_g<-function(b,ta){
  
  th<-unique(theta)[(length(unique(dataVals$theta))+1)/2]
  proportion<-rep(0,length(gvals))
  succ<-rep(0,length(gvals))
  total<-rep(0,length(gvals))
  fail<-rep(0,length(gvals))
  
  for (i in length(gvals)){
    succ[i]<-length(dataVals$lgCoral[dataVals$beta==b & dataVals$theta==th & dataVals$g==gvals[i] & dataVals$lgCoral>0.50 & (!is.na(dataVals$lgCoral))  & dataVals$tau==ta])
    total[i]<-length(dataVals$lgCoral[dataVals$beta==b & dataVals$theta==th & dataVals$g==gvals[i] & (!is.na(dataVals$lgCoral))  & dataVals$tau==ta])
    fail[i]<-total[i]-succ[i]
    proportion[i]<-succ[i]/total[i]
  }
  
  plot(gvals,proportion)
  
  
  g.logit<-glm(cbind(succ,fail)~gvals,family=binomial(link=logit))
  g.poisson<-glm(succ~gvals,family=poisson(link=log))
  g.ols<-glm(proportion~gvals)
  summary(g.logit)
  summary(g.poisson)
  summary(g.ols)
}


test_beta<-function(ta,th,g){
  
  proportion<-rep(0,length(beta))
  succ<-rep(0,length(beta))
  total<-rep(0,length(beta))
  fail<-rep(0,length(beta))
  
  for (i in length(beta)){
    succ[i]<-length(dataVals$lgCoral[dataVals$beta==beta[i] & dataVals$theta==th & dataVals$g==g & dataVals$lgCoral>0.50 & (!is.na(dataVals$lgCoral))  & dataVals$tau==ta])
    total[i]<-length(dataVals$lgCoral[dataVals$beta==beta[i] & dataVals$theta==th & dataVals$g==g & (!is.na(dataVals$lgCoral))  & dataVals$tau==ta])
    fail[i]<-total[i]-succ[i]
    proportion[i]<-succ[i]/total[i]
  }
  
  plot(beta,proportion)
  
  beta.logit<-glm(cbind(succ,fail)~beta,family=binomial(link=logit))
  beta.poisson<-glm(succ~beta,family=poisson(link=log))
  beta.ols<-glm(proportion~beta)
  summary(beta.logit)
  summary(beta.poisson)
  summary(beta.ols)
}

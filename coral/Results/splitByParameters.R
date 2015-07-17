if(!exists('dataVals'))
    {
        dataVals<-read.csv("trialsB.csv",header=TRUE)
    }
a<-length(unique(dataVals$beta))
b<-length(unique(dataVals$theta))
c<-seq(1,a,by=1)
d<-seq(1,b,by=1)
e<-length(unique(dataVals$g))
f<-seq(1,e,by=1)
v<-length(unique(dataVals$tau))
w<-seq(1,v,by=1)
l<-seq(1,a*b*e*v,by=1)
#succ<-0  				#number of sucess
#total<-0					#total number of trials
lgsucc<-0
lgtotal<-0
prop<-0					#output values
beta<-0					#beta values
theta<-0					#theta values
gval<-0            #g values
tau<-0
for (i in l)
{
  #succ[i]<-0
  #total[i]<-0
  lgsucc[i]<-0
  lgtotal[i]<-0  
  prop[i]<-0
  beta[i]<-0
  theta[i]<-0
  gval[i]<-0
  tau[i]<-0
}
k<-1
for(m in w)
{
for(h in f)
{
for(i in c)
{
  for(j in d)
  {
    #succ[k]<-length(dataVals$coral[dataVals$beta==unique(dataVals$beta)[i] & dataVals$theta==unique(dataVals$theta)[j] & dataVals$g==unique(dataVals$g)[h] & dataVals$coral>0.50 & (!is.na(dataVals$coral)) & dataVals$tau==unique(dataVals$tau)[m]])
    #total[k]<-length(dataVals$coral[dataVals$beta==unique(dataVals$beta)[i] & dataVals$theta==unique(dataVals$theta)[j] & dataVals$g==unique(dataVals$g)[h] & (!is.na(dataVals$coral))  & dataVals$tau==unique(dataVals$tau)[m]])
    lgsucc[k]<-length(dataVals$lgCoral[dataVals$beta==unique(dataVals$beta)[i] & dataVals$theta==unique(dataVals$theta)[j] & dataVals$g==unique(dataVals$g)[h] & dataVals$lgCoral>0.50 & (!is.na(dataVals$lgCoral))  & dataVals$tau==unique(dataVals$tau)[m]])
    lgtotal[k]<-length(dataVals$lgCoral[dataVals$beta==unique(dataVals$beta)[i] & dataVals$theta==unique(dataVals$theta)[j] & dataVals$g==unique(dataVals$g)[h] & (!is.na(dataVals$lgCoral))  & dataVals$tau==unique(dataVals$tau)[m]])
    ##prop[k]<-length(dataVals$coral[dataVals$beta==unique(dataVals$beta)[i] & dataVals$theta==unique(dataVals$theta)[j] & dataVals$coral>0.55 & (!is.na(dataVals$coral))])/length(dataVals$coral[dataVals$beta==unique(dataVals$beta)[i] & dataVals$theta==unique(dataVals$theta)[j] & (!is.na(dataVals$coral))])
    beta[k]<-unique(dataVals$beta)[i]
    theta[k]<-unique(dataVals$theta)[j]
    gval[k]<-unique(dataVals$g)[h]
    tau[k]<-unique(dataVals$tau)[m]
    k<-k+1
  }
}
}
}
#prop<-succ/total
lgprop<-lgsucc/lgtotal
#fail<-total-succ
lgfail<-lgtotal-lgsucc
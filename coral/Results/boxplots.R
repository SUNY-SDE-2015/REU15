#if(!exists('dataVals'))
#{
  dataVals<-read.csv("trials_lucy.csv",header=TRUE)
#}

thetaVals   <- unique(dataVals$theta)
numTheta    <- length(thetaVals)
thetaCuts   <- c(0,0.5*(thetaVals[2:numTheta]+thetaVals[1:(numTheta-1)]),max(thetaVals)+1.0)
thetaLabels <- format(c(thetaVals),digits=2)

dataVals$thetaVals <- cut(dataVals$theta,thetaCuts,includelowest=TRUE,labels=thetaLabels,right=FALSE)


for(beta in sort(unique(dataVals$beta)))
{
  for(g in sort(unique(dataVals$g)))
  {
    for(theta in sort(unique(dataVals$theta)))
    {
      which <- (dataVals$beta==beta)& (dataVals$g==g) &
        !is.na(dataVals$lgMacro)&!is.na(dataVals$lgCoral)&
        dataVals$theta==theta
      algae <- dataVals$lgMacro[which]
      coral <- dataVals$lgCoral[which]
      
      t    <- format(pi/theta,digits=2)
      b    <- format(beta,digits=2)
      gval <- format(g,digits=2)
      plot(algae,coral,
           xlim=c(0,1),ylim=c(0,1),
           main=substitute(paste("Coral/Macro Algae Endpoints for ",
                                 beta,"=",b,' and ',
                                 g,"=",gval,' and ',
                                 theta,"=",pi,'/(',t,')'),list(b=b,t=t,gval=g)),
           xlab="Macroalgae",ylab="Coral")
      
      #smoothScatter(algae,coral,
      #     xlim=c(0,1),ylim=c(0,1),
      #     main=substitute(paste("Coral/Macro Algae Endpoints for ",
      #         beta,"=",t,' and ',
      #         g,"=",gval,' and ',
      #         theta,"=",b),list(b=b,t=t)),
      #     xlab="Macroalgae",ylab="Coral")
      
      readline(prompt = "Press <Enter> for next plot")
    }
  }
  
  
  print(paste("Examination of impact of theta for beta = ",beta))
  which <- (dataVals$beta==beta)&
    !is.na(dataVals$macroalgae)&!is.na(dataVals$coral)
  algae <- dataVals$macroalgae[which]
  phi <- dataVals$thetaVals[which]
  results <- cut(algae,c(0.0,0.4,1.1),labels=c("C","M"),include.lowest=TRUE)
  fit1 <- glm(results~phi,family=binomial(link=logit))
  print(summary(fit1))
  
  
}
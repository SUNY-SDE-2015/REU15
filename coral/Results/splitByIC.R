data<-read.csv("trialsA.csv",header=TRUE)

thetaVals   <- unique(data$theta)
numTheta    <- length(thetaVals)
thetaCuts   <- c(0,0.5*(thetaVals[2:numTheta]+thetaVals[1:(numTheta-1)]),max(thetaCuts)+1.0)
thetaLabels <- format(thetaVals,digits=2)

data$thetaVals <- cut(data$theta,thetaCuts,includelowest=TRUE,labels=thetaLabels,right=FALSE)


for(beta in sort(unique(data$beta)))
    {
        for(theta in sort(unique(data$theta)))
            {
                which <- (data$beta==beta)&
                    !is.na(data$macroalgae)&!is.na(data$coral)&
                        data$theta==theta
                algae <- data$macroalgae[which]
                coral <- data$coral[which]

                t <- format(pi/theta,digits=2)
                b <- format(beta,digits=2)
                plot(algae,coral,
                     xlim=c(0,1),ylim=c(0,1),
                     main=substitute(paste("Coral/Macro Algae Endpoints for ",
                         beta,"=",b,' and ',
                         theta,"=",pi,'/(',t,')'),list(b=b,t=t)),
                     xlab="Macroalgae",ylab="Coral")

                #smoothScatter(algae,coral,
                #     xlim=c(0,1),ylim=c(0,1),
                #     main=substitute(paste("Coral/Macro Algae Endpoints for ",
                #         beta,"=",t,' and ',
                #         theta,"=",b),list(b=b,t=t)),
                #     xlab="Macroalgae",ylab="Coral")

                         
            }


        print(paste("Examination of impact of theta for beta = ",beta))
        which <- (data$beta==beta)&
            !is.na(data$macroalgae)&!is.na(data$coral)
        algae <- data$macroalgae[which]
        phi <- data$thetaVals[which]
        results <- cut(algae,c(0.0,0.4,1.1),labels=c("C","M"),include.lowest=TRUE)
        fit1 <- glm(results~phi,family=binomial(link=logit))
        print(summary(fit1))


    }

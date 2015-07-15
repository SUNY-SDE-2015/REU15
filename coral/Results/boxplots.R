data<-read.csv("trials_lucy.csv")

for(taus in sort(unique(data$tau))){
 for(gs in sort(unique(data$g))){
  for(betas in sort(unique(data$beta))){
    coral<-data$coral[(data$tau==taus)&(data$g==gs)&(data$beta==betas)]
    thetas<-data$theta[(data$tau==taus)&(data$g==gs)&(data$beta==betas)]
    boxplot(coral~thetas,
        main=substitute(paste("tau=",tau,", beta =",b,"and g =",gval),list(b=betas,gval=gs,tau=taus)),
        ylab="Coral surface",xlab="Theta")
    readline(prompt = "Press <Enter> for the macroalgae plot")
    
    #I tried to make the title look nicer, as Kelly did, but for some reason I couldn't ):
    
    macroalgae<-data$macroalgae[(data$tau==taus)&(data$g==gs)&(data$beta==betas)]
    thetas<-data$theta[(data$tau==taus)&(data$g==gs)&(data$beta==betas)]
    boxplot(macroalgae~thetas,
        main=substitute(paste("tau=",tau,", beta =",b,"and g =",gval),list(b=betas,gval=gs, tau=taus)),
        ylab="Macroalgae surface",xlab="Theta")

    readline(prompt = "Press <Enter> for next plot")
    }
  }
}

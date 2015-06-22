# Euler/Milstein

r <- 2; K <- 1; beta <- 0.25; Xzero <- 0.5;   # problem parameters
T <- 1; N <- 2^(11); dt <- T/N;               #             
M <- 500;                                     # number of paths sampled
R <- c(1, 16, 32, 64, 128);                   # Milstein stepsizes are R*dt

dW <- sqrt(dt)* matrix(runif(M*N),M);
Xmil <- matrix(M, 5);
for(p in 1:5) {
    Dt = R[p]*dt;
    L = N/R[p];
    Xtemp = Xzero*rep(M,1);
    for(j in 1:L) {
        Winc = sum(dW[, seq(R[p]*(j-1)+1, R[p]*(j))]);
        Xtemp = Xzero + Dt*r*Xtemp*(K-Xtemp) + beta*Xtemp*Winc + .5*beta^2*Xtemp*(Winc*Winc - Dt);
    }
    Xmil[,p] <- Xtemp; #store Milstein solution at t=1
}

Xref = Xmil[,1];                             # Reference solution
Xerr = abs(Xmil[,2:5] - repmat(Xref,1,4));   # Error in each path
mean(Xerr);                                  # Mean pathwise errors
Dtvals = dt*R(2:5);                          # Milstein timesteps used

#subplot(224)
plot(Dtvals, 
     mean(Xerr), 
     log = "xy", 
     type = "l");
lines(Dtvals, Dtvals, log = "xy",
#     axis([1e-3 1e-1 1e-4 1],
      xlab = 'dt',
      ylab = 'Sample average of | X(T) - X_L|',
      main = 'milstein')
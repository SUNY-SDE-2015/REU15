# RungeKutta Approximation and Heun's Method
r     <- 5;
a     <- 1/2;
b     <- 1/2;
alpha <- 1;
beta  <- 1;


final <- 5;
steps <- 100000;
dt    <- final/steps;


n    <- 0;
y    <- 1;
i    <- c(0, steps);
j    <- c(0, steps);
i[1] <- 0;
j[1] <- 1;


while(n < steps) {
    k1     <- dt*r*y;
    k2     <- dt*r*(y + beta*k1);
    y      <- y + a*k1 + b*k2;
    i[n+1] <- dt*(n+1);
    j[n+1] <- abs(y - exp(r*dt*(n + 1)));
    n      <- n + 1;
}

 plot(i, 
      j,
      ylim = c(0,30),
      #log = 'y',
      ylab = 'Error', 
      xlab = 'Time', 
      main = 'r = 5 and dt = 10^{-2}')
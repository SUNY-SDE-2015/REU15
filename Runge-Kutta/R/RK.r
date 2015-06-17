# RungeKutta Approximation and Heun's Method
r     <- -5;
a     <- 1/2;
b     <- 1/2;
alpha <- 1;
beta  <- 1;


final <- 25;
steps <- 1000;
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

plot( i, 
      j,
      type = 'l',
      xlim = c(0,25),
      ylim = c(0,5),
      # log = 'xy',
      ylab = 'Error', 
      xlab = 'Time', 
      main = 'r = -5 and dt = ')
      
file.create("/Users/kylemcgrath/Documents/SDE-REU-2015/SUNY-SDE-2015/REU15/Runge-Kutta/R/error-dat.csv");
data <- file("/Users/kylemcgrath/Documents/SDE-REU-2015/SUNY-SDE-2015/REU15/Runge-Kutta/R/error-dat.csv");
write(c("dt","step", "error"), file = data, ncolumns = 3, sep = ", ");
write(c(dt,i,j), append = TRUE, file = data, ncolumns = 3, sep = ", ")

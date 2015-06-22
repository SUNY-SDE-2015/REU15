
final <- 10;
steps <- 50;
dt    <- final / steps;

n <- 0;
y <- 1;
i <- c(0, steps);
j <- c(0, steps);

while( n < steps ) {
    y =  y + dt * -y;
    i[n+1] <- n *dt;
    j[n+1] <- y;
    n = n + 1;
}

plot(type = 'l', i, j, col = 4, ylim=c(-2,2));
lines(i,exp(-i), col = 2);
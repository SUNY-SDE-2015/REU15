huen <- function(steps,final)
{
	# RungeKutta Approximation and Heun's Method
	r     <- -5;
	a     <- 1/2;
	b     <- 1/2;
	alpha <- 1;
	beta  <- 1;
	
	#final <- 25;
	#steps <- 1000;
	dt    <- final/steps;
	
	n    <- 0;
	y    <- 1;
	#i    <- double(steps);
	#j    <- double(steps);
	
	
	while(n < steps) {
	    k1     <- dt*r*y;
	    k2     <- dt*r*(y + beta*k1);
	    y      <- y + a*k1 + b*k2;
	    #i[n+1] <- dt*(n+1);
	    #j[n+1] <- abs(y - exp(r*dt*(n + 1)));
	    n      <- n + 1;
	}
	
	#plot( 
    #  i, 
    #  j,
    #  type = 'l',
    #  xlim = c(0,25),
    #  ylim = c(0,5),
    #  log = 'xy',
    #  ylab = 'Error', 
    #  xlab = 'Time', 
    #  main = 'r = -5 and dt = ')

	
	return(y)
}


#x <- data.frame(dt=i,x=j);
#write.csv(x,file="blah blah.csv")


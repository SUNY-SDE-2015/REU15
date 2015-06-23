a <- 1
g <- 1
gamma <- 1
r <- 1
d <- 1

cMin <- -0.1
cMax <-  1.1
mMin <- -0.1
mMax <-  1.1
par(xpd=FALSE)

x1 <- seq(cMin, cMax, by=.01)
y1 <- (x1*(1-x1)-gamma*x1/(x1+D))/alpha     ### x null-cline

plot(x1,y1, xlim=c(cMin, cMax),ylim=c(mMin, mMax),col=2,type="l",
     main="Phase Plane for Coral Without Delay",
     xlab="Species x",ylab="Species y")
points (c(0,0),c(0,1.1),col=2,type="l")

y2 <- seq(mMin, mMax, by=.01)
x2 <- (rho*y2*(1-y2)-delta*y2/(y2+R))/beta  ### y null-cline
points (x2,y2,col=3,type="l")
points (c(0,1.1),c(0,0),col=3,type="l")

par(xpd=TRUE)

legend (.4,1.4,c("x Nullcline","y Nullcline"),col=c(2,3), lty=c(1,1))

vectorField <- function(x,y,alpha,beta,delta,gamma,rho,R,D)
    {
        fx <- x*x*(1-x)-alpha*x*y-gamma*x*x/(D+x)
        fy <- rho*y*y*(1-y)-beta*x*y-delta*y*y/(R+y)
        return(c(fx,fy))
    }

x1 <- seq(0, 1.1, by=.033)
y2 <- x1
vx <- numeric(length(x1))
vy <- numeric(length(y2))
dt <- 0.015;
for (x in x1)
{
    lupe <- 0
    for (y in y2)
        {
            lupe <- lupe + 1
            v <- vectorField(x,y,alpha,beta,delta,gamma,rho,R,D)
            vx[lupe] <- v[1]
            vy[lupe] <- v[2]
        }
    vecLength = sqrt(vx*vx+vy*vy)
    endx <- x + vx*dt/vecLength
    endy <- y2 + vy*dt/vecLength
    arrows(x,y2, endx, endy,length=dt*3, col=1,angle=20);
}

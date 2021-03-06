a <- 0.1
g <- 0.3
gamma <- 0.8
r <- 1
d <- 0.44

cMin <- -0.1
cMax <-  1.1
mMin <- -0.1
mMax <-  1.1
par(xpd=FALSE)
par(cex.main=1.0)
par(cex.lab=1.0)
par(cex.axis=1.0)

c1 <- seq(cMin, cMax, by=.01)
m1 <- (a/gamma*c1-(g/gamma)/(1.0-c1)+1.0-c1)    ### x null-cline

pdf("nullclines.pdf",width=10,height=10,pointsize=18)
plot(m1,c1, xlim=c(mMin, mMax),ylim=c(cMin, cMax),col=2,type="l",lwd=5,
     main="",xlab="",ylab="",cex.axis=1.5)
points (c(0.0,0.0),c(cMin,cMax),col=2,type="l",lwd=5)

c2 <- seq(cMin, cMax, by=.01)
m2 <- (r-r*c2-d)/(r+a)   ### y null-cline
points (m2,c2,col=3,type="l",lwd=5)
points (c(mMin,mMax),c(0.0,0.0),col=3,type="l",lwd=5)

par(xpd=TRUE)
par(adj=0.2)
title(main="Phase Plane for Coral Without Delay",
      xlab="Macroalgae Area",ylab="Coral Area")
legend (0.87,1.33,c("M Nullcline","C Nullcline"),col=c(2,3), lty=c(1,1), lwd=5)

vectorField <- function(mat,coral,a,g,gamma,r,d)
	    {
		        turf  <- 1.0 - mat - coral
		        fx    <- a*mat*coral-g*mat/(1.0-coral)+gamma*mat*turf
		        fy    <- r*turf*coral-d*coral-a*mat*coral
		        return(c(fx,fy))
		    }

m1 <- seq(0, 1.1, by=.045)
c1 <- m1
vx <- numeric(length(m1))
vy <- numeric(length(c1))
dt <- 0.025;
for (m in m1)
	{
		    lupe <- 0
		    for (c in c1)
			        {
				            lupe <- lupe + 1
				            v <- vectorField(m,c,a,g,gamma,r,d)
				            vx[lupe] <- v[1]
				            vy[lupe] <- v[2]
				        }
		    vecLength = sqrt(vx*vx+vy*vy)
		    endx <- m + vx*dt/vecLength
		    endy <- c1 + vy*dt/vecLength
		    arrows(m,c1, endx, endy,length=dt*3, col=1,angle=20,lwd=3);
		}

dev.off()
#dev.print(device=png,"nullclines.png",width=2048,height=2048)

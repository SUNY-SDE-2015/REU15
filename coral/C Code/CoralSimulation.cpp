#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double linear(long steps, double a,double gamma,double r,double d,double g,double *x,double *y, double dt)
	{
		double *z;
		for (double k=0;k<steps;k++)
			{
				x[0]=y[0]+(y[0]*(gamma-gamma*y[0]+(a-gamma)*y[1]-(g/(1-y[1]))))*dt;
				x[1]=y[1]+(y[1]*(r-d-(a+r)*y[0]-r*y[1]))*dt;
				//printf("%f\t%f\n",x[0],x[1]);
				z=x;
				x=y;
				y=x;
			}
		printf("%f\t%f\n",x[0],x[1]);
		return 0;
	}

int main(void)
	{
		long steps;
		double *x,*y,a,g,gamma,r,d,dt,final;
		int trials;
		final=1;
		trials=1;
		a=0.1;
		gamma=0.8;
		r=1;
		d=0.44;
		g=0.6;
		x=(double *) calloc(2,sizeof(double));
		y=(double *) calloc(2,sizeof(double));
		y[0]=0.8;
		y[1]=0.1;		
		dt=0.0000001;
		steps=long(final/dt);
		for (int k=0;k<trials;k++)
			{
				linear(steps,a,gamma,r,d,g,x,y,dt);
			}
	}

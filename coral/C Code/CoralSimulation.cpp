#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/* ****************************************************************
	 linear
	 
	 Function to evaluate the nonlinear operator.
**************************************************************** */
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

/* ****************************************************************
	 main
	 
	 The main function. Called at the start of the program.
**************************************************************** */
int main(int argc,char **argv)
	{
		long steps;
		double *x,*y;

		double a     = 0.1;
		double g     = 0.6;
		double gamma = 0.8;
		double r     = 1.0;
		double d     = 0.44;

		double dt,final;
		int trials;
		
		final=1;
		trials=1;

		x=(double *) calloc(2,sizeof(double));
		y=(double *) calloc(2,sizeof(double));
		y[0]=0.8;
		y[1]=0.1;		
		dt=0.0000001;
		steps=(long)(final/dt);
		for (int k=0;k<trials;k++)
			{
				linear(steps,a,gamma,r,d,g,x,y,dt);
			}

		free(x);
		free(y);
	}

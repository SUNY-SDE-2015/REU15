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
		long steps;    // The number of steps to take in a single simulation.
		double *x,*y;  // The variables used for the state of the system.

		/*
			Define the constants.

			These are the parameters that are used in the operator for the
			differential equations.
		 */
		double a     = 0.1;
		double g     = 0.6;
		double gamma = 0.8;
		double r     = 1.0;
		double d     = 0.44;

		double dt,final;    // The time step and the final time.
		int trials;         // The number of simulations to make.
		
		final=1;  // Set the final time.
		trials=1; // Set the number of trials to perform.

		// Allocate the space for the state of the system and define the
		// initial condition.
		x=(double *) calloc(2,sizeof(double));
		y=(double *) calloc(2,sizeof(double));

		// Set the time step and iterate through for the number of trials.
		dt=0.0000001;
		steps=(long)(final/dt);
		for (int k=0;k<trials;k++)
			{
				y[0]=0.8; // Initialize the value of x
				y[1]=0.1;	// Initialize the value of y

				// Perform a single simulation.
				linear(steps,a,gamma,r,d,g,x,y,dt);
			}

		free(x);
		free(y);
	}

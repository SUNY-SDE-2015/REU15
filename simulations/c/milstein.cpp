/* ***************************************************************************
   Approximation of the SDE
   dx = alpha x dt + beta x dW
   
   Makes use of the Euler-Maruyama approximation scheme.

   Date:
   Author:

***************************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


/* 
   Routine to generate two normally distributed random numbers.  This
   routine uses the Box-Muller transformation and generates two
   numbers.
*/
void normalDistRand(double stdDev,double* randomNumbers)
{
	/* Generate a random number in polar coordinates */
	double radius = sqrt(-2.0*log(drand48()));
	double angle  = 2.0*M_PI*drand48();

	/* transform the number back into cartesian coordinates */
	randomNumbers[0] = stdDev*radius*sin(angle);
	randomNumbers[1] = stdDev*radius*cos(angle);
}



/*
  Description goes here.
 */
double theoretical(FILE *fp,
				  int trial,
				  double alpha, double beta,
				  double dt,long steps) //does all the heavy lifting
	{

		double xzero=1;
		double w=0;
		double B[2];
		int calcRandom = 0;
		long m=0;
		
		/*int R=4;
		double Dt=R*dt;
		double L=1/(R*dt);*/
		double xtemptrue=1;
		double xtemp=1;
		double xmil=1;
		while (m++<steps)	
			{
				if(calcRandom==0)
					normalDistRand(sqrt(dt),B); //noise in most recent time step
				w=w+B[calcRandom];				//keeps a running tab on the noise

				//computes the theoretical values
				xtemptrue=xzero*exp((alpha-0.5*beta*beta)*(m*dt)+beta*w);

				//computes the Euler-Maruyama values
				xtemp += dt*alpha*xtemp+beta*xtemp*B[calcRandom];

				//computes the Milstein values
				xmil  += alpha*xmil*dt+
					beta*xmil*B[calcRandom]+0.5*beta*xmil*beta*(B[calcRandom]*B[calcRandom]-dt);

				
				calcRandom = (calcRandom+1)%2; // update which random number to use.
			}
		
		fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f\n",
				dt,alpha,beta,xtemptrue,xtemp,xmil,
				fabs(xtemptrue-xtemp),fabs(xtemptrue-xmil));
		return w;
	}


int main(void)
	{
		srand48(time(NULL));
		
		FILE*fp;
		fp=fopen("eulermaruyama_milstein.csv","w");	

        double dt,alpha,beta;
		long steps;
		long trials = 1;
	   
		fprintf(fp,"dt,alpha,beta,xTrue,xEuler,xMilstein,ErrorEuler,ErrorMilstein\n");

        for(dt= .001; dt <= .01; dt += .001)
	    	for(alpha=-2.0;alpha<=2.0;alpha+=0.5)
			    for(beta=0.25;beta<=4.0;beta+=0.25)
			    	{
			     		steps=(long)(1.0/dt);
			     		printf("%f\t%f\t%f\n",alpha,beta,dt);  //optional line.
					                                       // shows current progress
					                                       // but slows down program		

			    		for (long stepsTwo=0;stepsTwo<trials;++stepsTwo)
							    theoretical(fp,stepsTwo,alpha,beta,dt,steps);
			     	}

		fclose(fp);
		return 0;
	}

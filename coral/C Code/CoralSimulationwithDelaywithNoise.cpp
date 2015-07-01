#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void normalDistRand(double stdDev,double* randomNumbers)
	{
	    /* Generate a random number in polar coordinates */
		double radius = sqrt(-2.0*log(drand48()));
		double angle = 2.0*M_PI*drand48();
		/* transform the number back into cartesian coordinates */
		randomNumbers[0] = stdDev*radius*sin(angle);
		randomNumbers[1] = stdDev*radius*cos(angle);
	}

double linear(long steps,
			  double a,double gamma,double r,double d,double g,
			  double *x,double *y,double *z, double dt,
			  int n,
			  double beta,double tau,
			  FILE *fp,
			  int q,
			  double h,double s,double *v,double *w)
	{
		int m=0;
		int p;
		p=(m+n-1)%(n-1);
		double B[2];
		int calcRandom=0;
		for (long k=0;k<steps;k++)
		{
			if(calcRandom==0)
				normalDistRand(sqrt(dt),B); //noise in most recent time step

			//x[0]=y[m]+(y[m]*(gamma-gamma*y[m]+(a-gamma)*z[m])-(g*y[p]/(1-z[p])))*dt
			//x[0]=x[0]+beta*y[m]*(1-y[m])+0.5*beta*y[m]*(1-y[m])*beta*(1-2*y[m])*(B[calcRandom]*B[calcRandom]-dt);
			
			//Computes the deterministic component for Macroalgae
			x[0]=y[m]+(y[m]*(gamma-gamma*y[m]+(a-gamma)*z[m])-(g*y[p]/(1-z[p])))*dt;
			
			//Adds in the noise	
			//x[0]=x[0]+beta*y[m]*B[calcRandom]+0.5*beta*y[m]*beta*(B[calcRandom]*B[calcRandom]-dt);
			x[0]=x[0]+beta*y[m]*B[calcRandom]+0.5*beta*beta*y[m]*(B[calcRandom]*B[calcRandom]-dt);
			
			//Computes the deterministic component for Coral
			x[1]=z[m]+(z[m]*(r-d-(a+r)*y[m]-r*z[m]))*dt;
			
			
			x[2]=v[m]+(v[m]*(gamma-gamma*v[m]+(a-gamma)*w[m])-(g*v[p]/(1-w[p])))*dt;
			x[2]=x[2]+beta*v[m]*(1-v[m])*B[calcRandom]+0.5*beta*(1-2*v[m])*beta*v[m]*(1-v[m])*(B[calcRandom]*B[calcRandom]-dt);
			x[3]=w[m]+(w[m]*(r-d-(a+r)*v[m]-r*w[m]))*dt;
			
	        /****************************************************************
	            Account for extinction and overgrowing!!
			****************************************************************/
			for(int i=0;i<4;++i)
				{
					if(x[i]<0.0)
						x[i] = 0.0;
					else if (x[i]>1.0)
						x[i] = 1.0;
				}

			//Updates delay and cell index
			m=(m+1)%(n-1);
			p=(m+1)%(n-1);
			y[m]=x[0];
			z[m]=x[1];
			v[m]=x[2];
			w[m]=x[3];
			calcRandom = (calcRandom+1)%2; // update which random number to use.
			fprintf(fp,"%f,%f,%f,%f,%f,%f\n",x[0],x[1],1-x[0]-x[1],x[2],x[3],1-x[2]-x[3]);
		}

		//printf("%f\t%f\t%f\t%f\t%f\n",dt,beta,tau,x[0],x[1]);
		//fprintf(fp,"%f,%f,%f,%i,%f,%f,%f,%f,%f,%f\n",dt,beta,tau,q+1,h,s,1-h-s,x[0],x[1],1-x[0]-x[1]);
		return 0;
	}


/* ****************************************************************
	 main
	 
	 The main function. Called at the start of the program.
**************************************************************** */
int main(int argc,char **argv)
	{
		
		long steps;    // The number of steps to take in a single simulation.
		double *v,*w,*x,*y,*z;  // The variables used for the state of the system.

		/*
			Define the constants.

			These are the parameters that are used in the operator for the
			differential equations.
		 */
		double a     = 0.1;
		double g     = 0.3;
		double gamma = 0.8;
		double r     = 1.0;
		double d     = 0.44;
		double tau	 = 0.5;
		double beta  = 0.1;

		double dt,final;    // The time step and the final time.
		int trials;         // The number of simulations to make.
		
		final=25;  // Set the final time.
		trials=1; // Set the number of trials to perform.
		
		// Set the smallest time step
		dt=0.0001;
		
		// Sets the seed for the random numbers
		srand48(time(NULL));
		
		// The number of cells needed for the delay (changes with dt)
		int n;				
		n=(int) round(tau/dt);
		
		// Allocate the space for the state of the system
		x=(double *) calloc(4,sizeof(double));
		y=(double *) calloc(n,sizeof(double));		//macroalgae for multiplicative noise
		z=(double *) calloc(n,sizeof(double));		//coral for multiplicative noise
		v=(double *) calloc(n,sizeof(double));		//macroalgae for logistic noise
		w=(double *) calloc(n,sizeof(double));		//coral for logistic noise
		
		// Create a CSV File
		FILE*fp;
		fp=fopen("trials.csv","w");
		//fprintf(fp,"dt,beta,tau,trial,initMacro,initCoral,initTurf,macroalgae,coral,turf\n");
		fprintf(fp,"macroalgae,coral,turf,lgmacroalgae,lgcoral,lgturf\n");
		//printf("dt\t\tbeta\t\ttau\t\tMacroalgae\tCoral\n");
		
		
		for (double h=0.24;h<=0.24;h=h+0.1)
			for (double s=0.27;s<=0.27;s=s+0.1)
			{
				dt=0.0001;
				while (dt<=0.0001)
				{
					//printf("%f\t%f\n",dt,fmod(tau,dt));
					if ((int)round(10000*tau)%(int)round(dt*10000)==0)
					{
						//index = tau/dt;
						n=(int) round(tau/dt);
						//printf("%i\n",n);
						steps=(long)(final/dt);
						for (int k=0;k<trials;k++)
						{
							y[0]=h; //initial Macroalgae level
							z[0]=s; //initial Coral level
							v[0]=h;
							w[0]=s;
							for (int l=1;l<n;l++) //fills in "negative" times for both y and z
							{
								y[l]=y[0];
								z[l]=z[0];
								v[l]=v[0];
								w[l]=w[0];
							}
							fprintf(fp,"%f,%f,%f,%f,%f,%f\n",y[0],z[0],1-y[0]-z[0],v[0],w[0],1-v[0]-w[0]);
							linear(steps,a,gamma,r,d,g,x,y,z,dt,n,beta,tau,fp,k,h,s,v,w);
						}
					}
					dt=dt+0.0001;
				}
			}
		free(x);
		free(y);
		free(z);
		free(v);
		free(w);
		fclose(fp);
	    return 0;
	}

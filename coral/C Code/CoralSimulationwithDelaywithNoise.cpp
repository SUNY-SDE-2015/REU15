#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void normalDistRand(double stdDev,double* randomNumbers)
{
	/* Generate a random number in polar coordinates */
	double radius = sqrt(-2.0*log(drand48()));
	double angle  = 2.0*M_PI*drand48();

	/* transform the number back into cartesian coordinates */
	randomNumbers[0] = stdDev*radius*sin(angle);
	randomNumbers[1] = stdDev*radius*cos(angle);
}

double linear(long steps, double a,double gamma,double r,double d,double g,double *x,double *y,double *z, double dt,int n,double beta,double tau)
	{
		int m=0;
		int p;
		p=(m+n-1)%(n-1);
		double B[2];
		int calcRandom=0;
		for (double k=0;k<steps;k++)
			{
/*				
				x[0]=y[0]+(y[0]*(gamma-gamma*y[0]+(a-gamma)*y[1]-(g/(1-y[1]))))*dt;
				x[1]=y[1]+(y[1]*(r-d-(a+r)*y[0]-r*y[1]))*dt;
				//printf("%f\t%f\n",x[0],x[1]);
				z=x;
				x=y;
				y=x;
*/
				if(calcRandom==0)
					normalDistRand(sqrt(dt),B); //noise in most recent time step
				
				x[0]=y[m]+(y[m]*(gamma-gamma*y[m]+(a-gamma)*z[m])-(g*y[p]/(1-z[p])))*dt+beta*y[m]*(1-y[m])+0.5*beta*y[m]*(1-y[m])*beta*(1-2*y[m])*(B[calcRandom]*B[calcRandom]-dt);
				x[1]=z[m]+(z[m]*(r-d-(a+r)*y[m]-r*z[m]))*dt;
				if (x[0]<0)
					{
						x[0]=0;
					}
				if (x[1]<0)
					{
						x[1]=0;
					}
				m=(m+1)%(n-1);
				p=(m+1)%(n-1);									
				y[m]=x[0];
				z[m]=x[1];


				//printf("%i\t%i\n",m,p);
				
				calcRandom = (calcRandom+1)%2; // update which random number to use.
			}
		printf("%f\t%f\t%f\t%f\t%f\n",dt,beta,tau,x[0],x[1]);
		return 0;
	}

int main(void)
	{
		long steps;
		double *x,*y,*z,a,g,gamma,r,d,dt,final,tau,index,beta;
		int trials,n;
		final=50;		//change time to value
		trials=1;		//change number of trials
		a=0.1;
		gamma=0.8;
		r=1;
		d=0.44;
		g=0.8;			//change g value
		tau=1;		//change tau value
		beta=.01;		//change beta value
		//dt=0.0000001;	//change dt value
		//printf("%f\n",index);
		n=1;
		x=(double *) calloc(2,sizeof(double));
		dt=0.0001;
        printf("dt\t\tbeta\t\ttau\t\tMacroalgae\tCoral\n");
		while (dt<=0.01)
			{
				//printf("%f\t%f\n",dt,fmod(tau,dt));
				if (fmod(tau,dt)-dt>-0.0000001 || fmod(tau,dt)<0.0000001)
					{
						index=tau/dt;
						n=1;
						while (n<index)
							{
								n++;
							}
						if (n-0.25>index)
							{
								n--;
							}
						n++;
						//printf("%i\n",n);
						if (dt==0.0001)
							{
								y=(double *) calloc(n,sizeof(double));
								z=(double *) calloc(n,sizeof(double));
							}
						steps=long(final/dt);						
						for (int k=0;k<trials;k++)
							{
								y[0]=0.38692413985;		//initial Macroalgae level
								z[0]=0.6;		//initial Coral level
								for (int l=1;l<n;l++)		//fills in "negative" times for both y and z
									{
										y[l]=y[0];
										z[l]=z[0];
									}			
								linear(steps,a,gamma,r,d,g,x,y,z,dt,n,beta,tau);
							}
					}
				dt=dt+0.0001;
			}
		return 0;
	}

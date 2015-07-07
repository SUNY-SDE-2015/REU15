#include <QCoreApplication>

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265359
#endif

#define SHOW_PROGRESS

double drand48()
{
    return((double)(rand())/((double)RAND_MAX));
}

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
        long m=n-1;
        long p=0;
        double B[2];
        int calcRandom=0;
        for (long k=0;k<steps;k++)
        {
            if(calcRandom==0)
                normalDistRand(sqrt(dt),B); //noise in most recent time step

            //x[0]=y[m]+(y[m]*(gamma-gamma*y[m]+(a-gamma)*z[m])-(g*y[p]/(1-z[p])))*dt
            //x[0]=x[0]+beta*y[m]*(1-y[m])+0.5*beta*y[m]*(1-y[m])*beta*(1-2*y[m])*(B[calcRandom]*B[calcRandom]-dt);

            //Computes the deterministic component for Macroalgae
            x[0] = y[m]+(y[m]*(gamma-gamma*y[m]+(a-gamma)*z[m])-(g*y[p]/(1-z[p])))*dt;

            //Adds in the noise
            //x[0]=x[0]+beta*y[m]*B[calcRandom]+0.5*beta*y[m]*beta*(B[calcRandom]*B[calcRandom]-dt);
            x[0] += beta*y[m]*B[calcRandom]+0.5*beta*beta*y[m]*(B[calcRandom]*B[calcRandom]-dt);

            //Computes the deterministic component for Coral
            x[1] = z[m]+(z[m]*(r-d-(a+r)*y[m]-r*z[m]))*dt;


            x[2] =  v[m]+(v[m]*(gamma-gamma*v[m]+(a-gamma)*w[m])-(g*v[p]/(1-w[p])))*dt;
            x[2] += beta*v[m]*(1-v[m])*B[calcRandom]+0.5*beta*(1-2*v[m])*beta*v[m]*(1-v[m])*(B[calcRandom]*B[calcRandom]-dt);
            x[3] =  w[m]+(w[m]*(r-d-(a+r)*v[m]-r*w[m]))*dt;

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
            m=(m+1)%n;
            p=(p+1)%n;
            y[m]=x[0];
            z[m]=x[1];
            v[m]=x[2];
            w[m]=x[3];
            calcRandom = (calcRandom+1)%2; // update which random number to use.
            //fprintf(fp,"%f,%f,%f,%f,%f,%f\n",x[0],x[1],1-x[0]-x[1],x[2],x[3],1-x[2]-x[3]);
        }

        //printf("%f\t%f\t%f\t%f\t%f\n",dt,beta,tau,x[0],x[1]);
        fprintf(fp,"%f,%f,%f,%f,%i,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",dt,beta,g,tau,q+1,h,s,1-h-s,x[0],x[1],1-x[0]-x[1],x[2],x[3],1-x[2]-x[3]);
        return 0;
    }

void Trials(double a,double gamma,double r,double d,double g,double dt,double macroSaddle, double coralSaddle,
            double beta,double tau,double final,int trials,long steps, bool arc){

    double *v,*w,*x,*y,*z;
    FILE*fp;
    if (arc==false){
        fp=fopen("trials.csv","w");
        fprintf(fp,"dt,beta,g,tau,trial,initMacro,initCoral,initTurf,macroalgae,coral,turf,lgMacro,lgCoral,lgTurf\n");
    }
    if (arc==true)
    fp=fopen("arc.csv","a+");


    int n;
    n=(int)abs(tau/dt+0.5);

    x=(double *) calloc(4,sizeof(double));
    y=(double *) calloc(n,sizeof(double));		//macroalgae for multiplicative noise
    z=(double *) calloc(n,sizeof(double));		//coral for multiplicative noise
    v=(double *) calloc(n,sizeof(double));		//macroalgae for logistic noise
    w=(double *) calloc(n,sizeof(double));		//coral for logistic noise
        while (dt<=0.0001)
        {
            //printf("%f\t%f\n",dt,fmod(tau,dt));
#ifdef SHOW_PROGRESS
            std::cout << "dt = " << dt << std::endl;
#endif
            if ((int)(10000*tau+.5)%(int)(dt*10000+.5)==0)
            {
                //index = tau/dt;
                n=(int)(tau/dt+.5);
                //printf("%i\n",n);
                steps=(long)(final/dt);
                for (int k=0;k<trials;k++)
                {
                    y[0]=macroSaddle; //initial Macroalgae level
                    z[0]=coralSaddle; //initial Coral level
                    v[0]=macroSaddle;
                    w[0]=coralSaddle;
                    for (int l=1;l<n;l++) //fills in "negative" times for both y and z
                    {
                        y[l]=y[0];
                        z[l]=z[0];
                        v[l]=v[0];
                        w[l]=w[0];
                    }
                    //fprintf(fp,"%f,%f,%f,%f,%f,%f\n",y[0],z[0],1-y[0]-z[0],v[0],w[0],1-v[0]-w[0]);
                    linear(steps,a,gamma,r,d,g,x,y,z,dt,n,beta,tau,fp,k,macroSaddle,coralSaddle,v,w);
#ifdef SHOW_PROGRESS
                    if(k%20 == 0)
                        std::cout << "  Simulation number " << k << std::endl;
#endif
                }
            }
            dt=dt+0.0001;
        }
    free(x);
    free(y);
    free(z);
    free(v);
    free(w);
    //}
    //}

fclose(fp);


}

void Arc(double a,double gamma,double r,double d,double g,double dt,double beta,double tau,
         double final,int trials,long steps, double radius, int number_points){

    int n=0;
    double angle,x,y;
    FILE*fp;
        fp=fopen("arc.csv","w");
        fprintf(fp,"dt,beta,g,tau,trial,initMacro,initCoral,initTurf,macroalgae,coral,turf,lgMacro,lgCoral,lgTurf\n");
        fclose(fp);
    while(n<number_points){
        angle=2.0*M_PI*drand48();
        x=fabs(radius*cos(angle));
        y=fabs(radius*sin(angle));
        std::cout << "Point number " << n << std::endl;
        Trials(a,gamma,r,d,g,dt,x,y,beta,tau,final,trials,steps,true);
        n++;
    }


}

/* ****************************************************************
     main

     The main function. Called at the start of the program.
**************************************************************** */
int main(int argc, char *argv[])
{
    QCoreApplication b(argc, argv);
        /*
            Define the constants.

            These are the parameters that are used in the operator for the
            differential equations.
         */
        double a    	 = 0.1;
        double g    	 = 0.3;
        double gamma	 = 0.8;
        double r    	 = 1.0;
        double d    	 = 0.44;
        double beta 	 = .5;
        double chi		 = r*gamma/(r+a)-gamma+a;					//Intermediate Step
        double xi		 = -(d*gamma/(r+a)+a);						//Intermediate Step
        double cbar		 = (-xi-sqrt(xi*xi-4*chi*g))/(2*chi);		//Intermediate Step
        double coralSaddle		 = 1-cbar;									//Saddle point value for coral
        double macroSaddle		 = (r-r*coralSaddle-d)/(r+a);						//Saddle point value for macroalgae
        //double gZero	 = ((d*a*r+d*d)*(gamma-a))/(r*r);
        //double gOne		 = (gamma*(a+d))/(a+r);
        //double omega	 = sqrt((r*r*(g*g-gZero*gZero))/(d*d));
        //double tauZero	 = (1/omega)*acos(gZero/g);

        double final=10;    // The time step and the final time.
        int trials=100;         // The number of simulations to make.

        // Set the smallest time step
        double dt=0.0001;

        // Set tau
        double tau = 0.5;

        // The number of steps to take in a single simulation.
        long steps=final/dt;

        // Sets the seed for the random numbers
        srand(time(NULL));


        //Trials(a,gamma,r,d,g,dt,macroSaddle,coralSaddle,beta,tau,final,trials,steps,false);

        //The equilibrium points
        double macroalgae=0.1*(1-g/gamma);
        double coral=0.1*(1-d/r);
        double radius=sqrt(pow(macroalgae,2)+pow(coral,2));
        int number_points=1000;

        final=100;
        trials=1;
        Arc(a,gamma,r,d,g,dt,beta,tau,final,trials,steps,radius,number_points);
        std::cout << " Done ";
        return b.exec();
}

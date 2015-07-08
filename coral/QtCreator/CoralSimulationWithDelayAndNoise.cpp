#include <QCoreApplication>
#include <QtDebug>

#include <iostream>
#include <fstream>
//#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cerrno>
#include <cmath>
#include <thread>
#include <mutex>


#define OUTPUT_FILE "./trials.csv"
#define NUMBER_THREADS 7
#define USE_MULTIPLE_THREADS
// create a mutex that is used to protect the writing of the data to
// the file.
std::mutex writeToFile;



#ifndef M_PI
#define M_PI 3.14159265359
#endif

#define SHOW_PROGRESS
//#define SHOW_INTERMEDIATE
#define BASE_DT 0.0001

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


void printToCSVFile(double dt, double beta, double g,double tau,
                    double q, double theta, double h,double s,
                    double *x,std::ofstream *fp)
{

    std::lock_guard<std::mutex> guard(writeToFile);  // Make sure that
                                                       // this routine
                                                       // can only
                                                       // access the file once
                                                       // at any one time.

    *fp << dt << ","
        << beta << ","
        << g << ","
        << tau << ","
        << q+1 << ","
        << theta << ","
        << h << ","
        << s << ","
        << 1-h-s << ","
        << x[0] << ","
        << x[1] << ","
        << 1-x[0]-x[1] << ","
        << x[2] << ","
        << x[3] << ","
        << 1-x[2]-x[3]
        << std::endl;

    (*fp).flush();
}


void linear(long steps,
              double a,double gamma,double r,double d,double g,
              double *x,double *y,double *z, double dt,
              int n,
              double beta,double tau,
              std::ofstream *fp,
              int q,
              double h,double s,double *v,double *w, double theta)
    {

        long m=n-1;
        long p=0;
        double B[2];
        int calcRandom=0;

        // Step through every time step.
        for (long k=0;k<steps;k++)
        {

            // Calculate a new set of random numbers if necessary.
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

        }

        //printf("%f\t%f\t%f\t%f\t%f\n",dt,beta,tau,x[0],x[1]);
        printToCSVFile(dt,beta,g,tau,q,theta,h,s,x,fp);

 #ifdef SHOW_INTERMEDIATE
        qDebug() << dt << ","
            << beta << ","
            << g << ","
            << tau << ","
            << q+1 << ","
            << h << ","
            << s << ","
            << 1-h-s << ","
            << x[0] << ","
            << x[1] << ","
            << 1-x[0]-x[1] << ","
            << x[2] << ","
            << x[3] << ","
            << 1-x[2]-x[3];
        qDebug() << errno << EDOM << ERANGE;
#endif

    }


/* ****************************************************************
     main

     The main function. Called at the start of the program.
**************************************************************** */
int main(int argc, char *argv[])
{
    QCoreApplication b(argc, argv);

        long steps;    // The number of steps to take in a single simulation.
        double *v,*w,*x,*y,*z;  // The variables used for the state of the system.

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
    //	double tau		 = 0.5;
        double beta 	 = .5;
        double chi		 = r*gamma/(r+a)-gamma+a;					//Intermediate Step
        double xi		 = -(d*gamma/(r+a)+a);						//Intermediate Step
        double cbar		 = (-xi-sqrt(xi*xi-4*chi*g))/(2*chi);		//Intermediate Step
        double coralSaddle		 = 1-cbar;						    //Saddle point value for coral
        double macroSaddle		 = (r-r*coralSaddle-d)/(r+a);		//Saddle point value for macroalgae
        double gZero	 = ((d*a*r+d*d)*(gamma-a))/(r*r);
        double gOne		 = (gamma*(a+d))/(a+r);
        double omega	 = sqrt((r*r*(g*g-gZero*gZero))/(d*d));
        double tauZero	 = (1/omega)*acos(gZero/g);

        double dt;          // Set the initial time step
        long   numberDT;    // The current iteration for the number assocated with the value of dt.
        double final;       // The final time for each simulation.
        long   trials;      // The number of simulations to make.

<<<<<<< HEAD
        final=50;  // Set the final time.
        trials=400; // Set the number of trials to perform.
=======
        final=50.0;  // Set the final time.
        trials=50;   // Set the number of trials to perform.
>>>>>>> afa3c18ff2f48dc4ee88e7e8158f5b78dae8cb48


        // Set the time delay, tau
        double tau = 0.5;

        // set up the variables for using different approximations on different threads.
        std::thread simulation[NUMBER_THREADS];
        int numberThreads = 0;

        // Sets the seed for the random numbers
        srand(time(NULL));



        // Create a CSV File
        std::ofstream fp;
                //String fileName = "trials-g" + std::toString(g) + "-tau" + std::toString(tau);
        fp.open(OUTPUT_FILE,std::ios::out | std::ios::trunc);

        fp << "dt,beta,g,tau,trial,theta,initMacro,initCoral,initTurf,macroalgae,coral,turf,lgMacro,lgCoral,lgTurf" << std::endl;


/*		for (g=0.1;g<=0.8;g=g+0.02)
            {
                //Redefine initial conditions and critical points for varying g
                cbar		 = (-xi-sqrt(xi*xi-4*chi*g))/(2*chi);
                coralSaddle		 = 1-cbar;
                macroSaddle		 = (r-r*coralSaddle-d)/(r+a);
                omega	 = sqrt((r*r*(g*g-gZero*gZero))/(d*d));
                tauZero	 = (1/omega)*acos(gZero/g);
                if (coralSaddle>0 && coralSaddle<1 && macroSaddle>0 && macroSaddle<1 && tauZero>0)
            for (double aleph=0;aleph<=5;aleph=aleph+1)
            {

                tau=.3*(1+aleph)*tauZero;
                dt=0.0001;
*/

           // Determine the number of time steps required to move back to the delay in time.
           // The number of cells needed for the delay (changes with dt)
            int n;
            n=(int)(tau/BASE_DT+0.5);

            // Allocate the space for the states of the system
            x=(double *) calloc(4,sizeof(double));
            y=(double *) calloc(n,sizeof(double));		//macroalgae for multiplicative noise
            z=(double *) calloc(n,sizeof(double));		//coral for multiplicative noise
            v=(double *) calloc(n,sizeof(double));		//macroalgae for logistic noise
            w=(double *) calloc(n,sizeof(double));		//coral for logistic noise
<<<<<<< HEAD
            //for(numberDT=1;numberDT<5;++numberDT)
                    dt = BASE_DT; //(double)numberDT;
                    //printf("%f\t%f\n",dt,fmod(tau,dt));
=======


            // Make different approximations for different values for the time steps.
            for(numberDT=1;numberDT<5;++numberDT)
                {
                    dt = BASE_DT*(double)numberDT;
>>>>>>> afa3c18ff2f48dc4ee88e7e8158f5b78dae8cb48
#ifdef SHOW_PROGRESS
                    std::cout << "dt = " << dt << std::endl;
#endif
                    if ((int)(10000.0*tau+.5)%(int)(dt*10000.0+.5)==0)
                    {
                        //index = tau/dt;
                        n=(int)(tau/dt+.5);
                        //printf("%i\n",n);
                        steps=(long)(final/dt);
<<<<<<< HEAD
                        for (double theta=0;theta<=M_PI/2;theta+=(M_PI*0.5)*0.025)
=======

                        // Make an approximation for different initial conditions.
                        // Make an arc through 0 to pi/2 radians from the origin.
                        for (double theta=0;theta<=M_PI/2;theta+=(M_PI*0.5)*0.05)
>>>>>>> afa3c18ff2f48dc4ee88e7e8158f5b78dae8cb48
                        {
                            for (int k=0;k<trials;k++)
                            {
                                y[0]=0.06*cos(theta); //initial Macroalgae level
                                z[0]=0.06*sin(theta); //initial Coral level
                                v[0]=0.06*cos(theta);
                                w[0]=0.06*sin(theta);
                                for (int l=1;l<n;l++) //fills in the past times for y, z, v, and w
                                {
                                    y[l]=y[0];
                                    z[l]=z[0];
                                    v[l]=v[0];
                                    w[l]=w[0];
                                }
                                //fprintf(fp,"%f,%f,%f,%f,%f,%f\n",y[0],z[0],1-y[0]-z[0],v[0],w[0],1-v[0]-w[0]);
#ifdef USE_MULTIPLE_THREADS


                                // Make a simulation for each of the available threads.
                                // First check to see how many threads are running.
                                if(numberThreads >= NUMBER_THREADS)
                                {
                                    // There are too many threads. Wait for each run to end.
                                    while(numberThreads>0)
                                    {
#ifdef THREAD_DEBUG
                                        std::cout << "Waiting on thread "
                                                  << simulation[numberThreads-1].get_id()
                                                  << std::endl;
#endif
                                        simulation[--numberThreads].join();
                                    }
                                }


                                // Make a run in a separate thread.
                                simulation[numberThreads++] = std::thread(linear,
                                                                          steps,a,gamma,r,d,g,x,y,z,dt,n,beta,tau,&fp,k,y[0],z[0],v,w,theta);

#else

                                // Ignore the different threads. Just make one approximation in this one thread.
                                linear(steps,a,gamma,r,d,g,x,y,z,dt,n,beta,tau,&fp,k,y[0],z[0],v,w,theta);

#endif


#ifdef SHOW_PROGRESS
                                if(k%20 == 0)
                                    std::cout << "  Simulation number " << k << std::endl;
#endif

                            }
                        }
            }

            // Free up the allocated memory.
            free(x);
            free(y);
            free(z);
            free(v);
            free(w);
            //}
            //}

        fp.close();

#ifdef SHOW_PROGRESS
            std::cout << "all done" << std::endl;
#endif

        return b.exec();
}

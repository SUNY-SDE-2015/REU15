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
#define NUMBER_THREADS 1
#define USE_MULTIPLE_THREADS
#define LOGISTIC_NOISE
// create a mutex that is used to protect the writing of the data to
// the file.
std::mutex writeToFile;



#ifndef M_PI
#define M_PI 3.14159265359
#endif

#define SHOW_PROGRESS
//#define SHOW_INTERMEDIATE
//#define THREAD_DEBUG
#define BASE_DT 0.00001
#define RADIUS 0.06


#ifndef Q_OS_UNIX
double drand48()
{
    return((double)(rand())/((double)RAND_MAX));
}
#endif


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
                    double *reefState,std::ofstream *fp)
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
        << reefState[0] << ","
        << reefState[1] << ","
        << 1-reefState[0]-reefState[1]
        << std::endl;

    (*fp).flush();
}


void linear(long steps,
              double a,double gamma,double r,double d,double g,
              double *reefState,double *y,double *z, double dt,
              int n,
              double beta,double tau,
              std::ofstream *fp,
              int q,
              double h,double s,double theta)
    {

        long m=n-1;
        long p=0;
        double B[2];
        int calcRandom=0;

#ifdef THREAD_DEBUG
        std::cout << "My thread id: " << std::this_thread::get_id() << std::endl;
#endif

        // Step through every time step.
        for (long k=0;k<steps;k++)
        {

            // Calculate a new set of random numbers if necessary.
            if(calcRandom==0)
                normalDistRand(sqrt(dt),B); //noise in most recent time step

#ifdef LOGISTIC_NOISE
            // Update the macroalgae term
            reefState[0] =  y[m]+(y[m]*(gamma-gamma*y[m]+(a-gamma)*z[m])-(g*y[p]/(1.0-z[p])))*dt
                    +  beta*y[m]*(1.0-y[m])*B[calcRandom]+0.5*beta*(1.0-2.0*y[m])*beta*y[m]*(1-y[m])*(B[calcRandom]*B[calcRandom]-dt);
            //Computes the deterministic component for Coral
            reefState[1] =  z[m]+(z[m]*(r-d-(a+r)*y[m]-r*z[m]))*dt;
#else
            // Update the macroalgae term
            reefState[0] += beta*y[m]*B[calcRandom]+0.5*beta*beta*y[m]*(B[calcRandom]*B[calcRandom]-dt);
            //Computes the deterministic component for Coral
            reefState[1] = z[m]+(z[m]*(r-d-(a+r)*y[m]-r*z[m]))*dt;
#endif

            /****************************************************************
                Account for extinction and overgrowing!!
            ****************************************************************/
            for(int i=0;i<2;++i)
                {
                    if(reefState[i]<0.0)
                        reefState[i] = 0.0;
                    else if (reefState[i]>1.0)
                        reefState[i] = 1.0;
                }

            //Updates delay and cell index
            m=(m+1)%n;
            p=(p+1)%n;
            y[m]=reefState[0];
            z[m]=reefState[1];
            calcRandom = (calcRandom+1)%2; // update which random number to use.

        }

        //printf("%f\t%f\t%f\t%f\t%f\n",dt,beta,tau,reefState[0],reefState[1]);
        printToCSVFile(dt,beta,g,tau,q,theta,h,s,reefState,fp);

 #ifdef SHOW_INTERMEDIATE
        qDebug() << dt << ","
            << beta << ","
            << g << ","
            << tau << ","
            << q+1 << ","
            << h << ","
            << s << ","
            << 1-h-s << ","
            << reefState[0] << ","
            << reefState[1] << ","
            << 1-reefState[0]-reefState[1];
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


        long steps;           // The number of steps to take in a single simulation.
        double reefState[2];  // The state of the system, coral and macroalgae
        double *y,*z;         // The variables used for the state of the system.

        /*
            Define the constants.

            These are the parameters that are used in the operator for the
            differential equations.
         */
        double a    	 = 0.1;
        double g    ;//	 = 0.3;
        double gamma	 = 0.8;
        double r    	 = 1.0;
        double d    	 = 0.44;

        double tau;
        double beta ;//	 = .5;
//        double gZero	 = ((d*a*r+d*d)*(gamma-a))/(r*r);
//        double gOne		 = (gamma*(a+d))/(a+r);
//        double chi		 = r*gamma/(r+a)-gamma+a;					//Intermediate Step
//        double xi		 = -(d*gamma/(r+a)+a);						//Intermediate Step
//        double cbar		 = (-xi-sqrt(xi*xi-4*chi*g))/(2*chi);		//Intermediate Step
//        double coralSaddle		 = 1-cbar;						    //Saddle point value for coral
//        double macroSaddle		 = (r-r*coralSaddle-d)/(r+a);		//Saddle point value for macroalgae


//        long   numberDT;    // The current iteration for the number assocated with the value of dt.
//        double omega	 = sqrt((r*r*(g*g-gZero*gZero))/(d*d));
//        double tauZero	 = (1/omega)*acos(gZero/g);

        double dt        = BASE_DT;          // Set the initial time step
        double final;       // The final time for each simulation.
        long   trials;      // The number of simulations to make.

        final=50.0;  // Set the final time.
        trials=400;   // Set the number of trials to perform.


        // Set the time delay, tau
        //double tau = 0.5;

        // set up the variables for using different approximations on different threads.
        std::thread simulation[NUMBER_THREADS];
        int numberThreads = 0;

        // Sets the seed for the random numbers
#ifdef Q_OS_UNIX
        srand48(time(NULL));
        //qDebug() << "This is a POSIX system!" ;
#else
        srand(time(NULL));
#endif



        // Create a CSV File
        std::ofstream fp;
                //String fileName = "trials-g" + std::toString(g) + "-tau" + std::toString(tau);
        fp.open(OUTPUT_FILE,std::ios::out | std::ios::trunc);

#ifdef LOGISTIC_NOISE
        fp << "dt,beta,g,tau,trial,theta,initMacro,initCoral,initTurf,lgMacro,lgCoral,lgTurf" << std::endl;
#else
        fp << "dt,beta,g,tau,trial,theta,initMacro,initCoral,initTurf,macroalgae,coral,turf" << std::endl;
#endif

        for(tau = .2; tau <= .4; tau += .2 )
        {
           // Determine the number of time steps required to move back to the delay in time.
           // The number of cells needed for the delay (changes with dt)
            int n;
            if(tau > 0.0)
                n=(int)(tau/BASE_DT+0.5);
            else
                n = 1;
            // Allocate the space for the states of the system
            y=(double *) calloc(n,sizeof(double));		//macroalgae density in the past
            z=(double *) calloc(n,sizeof(double));		//coral density in the past

            if((y==NULL) || (z==NULL))
            {
                std::cout << "Error - unable to allocate necessary memory." << std::endl;
                free(y);
                free(z);
                return b.exec();
            }

            for(g=.2; g<.8; g += .2)
            {
                //double omega	 = sqrt((r*r*(g*g-gZero*gZero))/(d*d));
                // double tauZero	 = (1/omega)*acos(gZero/g);

                // Make different approximations for different values of beta.
                for(beta=.2;beta<=1.0; beta += .2)
                {


#ifdef SHOW_PROGRESS
                    std::cout << "dt = " << dt << std::endl;
#endif
                        //index = tau/dt;

                        //printf("%i\n",n);
                        steps=(long)(final/dt);

                        for (double theta=0;theta<=M_PI/2;theta+=(M_PI*0.5)*0.025)
                        {

                        // Make an approximation for different initial conditions.
                        // Make an arc through 0 to pi/2 radians from the origin.

                            for (int k=0;k<trials;k++)
                            {
                                y[0] = RADIUS*cos(theta); //initial Macroalgae level
                                z[0] = RADIUS*sin(theta); //initial Coral level
                                for (int l=1;l<n;l++) //fills in the past times for y, z, v, and w
                                {
                                    y[l]=y[0];
                                    z[l]=z[0];
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
                                } // if(numberThreads)


                                // Make a run in a separate thread.
                                simulation[numberThreads++] = std::thread(linear,
                                                                          steps,a,gamma,r,d,g,reefState,y,z,dt,n,beta,tau,&fp,k,y[0],z[0],theta);

#else

                                // Ignore the different threads. Just make one approximation in this one thread.
                                linear(steps,a,gamma,r,d,g,reefState,y,z,dt,n,beta,tau,&fp,k,y[0],z[0],theta);

#endif


#ifdef SHOW_PROGRESS
                                if(k%20 == 0)
                                    std::cout << "  Simulation number " << k << std::endl;
#endif

                            } // for(k<trials)
                        } // for(theta)
                   } // for(beta)
            } // for(g)

            // Free up the allocated memory.
            free(y);
            free(z);

        } // for(tau)


        fp.close();

#ifdef SHOW_PROGRESS
        std::cout << "all done" << std::endl;
#endif

        return b.exec();
}

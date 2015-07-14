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


/* ********************************************
 *  Define the values used for which parameters
 *  are sampled.
 ********************************************** */

#define TAU_START    0.2
#define TAU_END      0.6
#define NUMBER_TAU   3

#define G_START      0.2
#define G_END        0.6
#define NUMBER_G     3

#define BETA_START   0.2
#define BETA_END     1.2
#define NUMBER_BETA  6

#define THETA_START  0.0
#define THETA_END    M_PI*0.5
#define NUMBER_THETA 20

#define NUMBER_TRIALS 3
#define FINAL_TIME    35.0

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
        << reefState[1]
        << std::endl;

    (*fp).flush();
}


void linear(long steps,
              double a,double gamma,double r,double d,double g,
              double *reefState,double *macroalgaePast,double *coralPast, double dt,
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
            reefState[0] =  macroalgaePast[m] +
                    (macroalgaePast[m]*(gamma-gamma*macroalgaePast[m]+(a-gamma)*coralPast[m])-(g*macroalgaePast[p]/(1.0-coralPast[p])))*dt
                    + beta*macroalgaePast[m]*(1.0-macroalgaePast[m])*B[calcRandom]
                    + 0.5*beta*(1.0-2.0*macroalgaePast[m])*beta*macroalgaePast[m]*(1-macroalgaePast[m])*(B[calcRandom]*B[calcRandom]-dt);
            //Computes the deterministic component for Coral
            reefState[1] =  coralPast[m]+(coralPast[m]*(r-d-(a+r)*macroalgaePast[m]-r*coralPast[m]))*dt;
#else
            // Update the macroalgae term
            reefState[0] += beta*macroalgaePast[m]*B[calcRandom]+0.5*beta*beta*macroalgaePast[m]*(B[calcRandom]*B[calcRandom]-dt);
            //Computes the deterministic component for Coral
            reefState[1] = coralPast[m]+(coralPast[m]*(r-d-(a+r)*macroalgaePast[m]-r*coralPast[m]))*dt;
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
            macroalgaePast[m]=reefState[0];
            coralPast[m]=reefState[1];
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


        long steps;             // The number of steps to take in a single simulation.
        double reefState[2];    // The state of the system, coral and macroalgae
        double *macroalgaePast; // Past values for the macroalgae density
        double *coralPast;      // Past values for the coral density

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
        double theta;
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
        fp << "dt,beta,g,tau,trial,theta,initMacro,initCoral,initTurf,lgMacro,lgCoral" << std::endl;
#else
        fp << "dt,beta,g,tau,trial,theta,initMacro,initCoral,initTurf,macroalgae,coral" << std::endl;
#endif

        for(int tauStep=0;tauStep<NUMBER_TAU;++tauStep)
        {
            tau = TAU_START+(TAU_END-TAU_START)/((double)(NUMBER_TAU-1))*((double)tauStep);

           // Determine the number of time steps required to move back to the delay in time.
           // The number of cells needed for the delay (changes with dt)
            int n;
            if(tau > 0.0)
                n=(int)(tau/BASE_DT+0.5);
            else
                n = 1;
            // Allocate the space for the states of the system
            macroalgaePast = (double *) calloc(n,sizeof(double));		//macroalgae density in the past
            coralPast      = (double *) calloc(n,sizeof(double));		//coral density in the past

            if((macroalgaePast==NULL) || (coralPast==NULL))
            {
                std::cout << "Error - unable to allocate necessary memory." << std::endl;
                free(macroalgaePast);
                free(coralPast);
                return b.exec();
            }

            for(int gStep=0;gStep<NUMBER_G;++gStep)
            {
                g = G_START+(G_END-G_START)/((double)(NUMBER_G-1))*((double)gStep);

                //double omega	 = sqrt((r*r*(g*g-gZero*gZero))/(d*d));
                // double tauZero	 = (1/omega)*acos(gZero/g);

                // Make different approximations for different values of beta.
                for(int betaStep=0;betaStep<NUMBER_BETA;++betaStep)
                {
                    beta = BETA_START+(BETA_END-BETA_START)/((double)(NUMBER_BETA-1))*((double)betaStep);


#ifdef SHOW_PROGRESS
                    std::cout << "tau = " << tau << " and g = " << g << std::endl;
#endif
                        //index = tau/dt;

                        //printf("%i\n",n);
                        steps=(long)(FINAL_TIME/dt);

                        for(int thetaStep=1;thetaStep<=NUMBER_THETA;++thetaStep)
                        {

                            theta = THETA_START+(THETA_END-THETA_START)/((double)(NUMBER_THETA+1))*((double)thetaStep);

                            // Make an approximation for different initial conditions.
                            // Make an arc through 0 to pi/2 radians from the origin.

                            for (int k=0;k<NUMBER_TRIALS;k++)
                            {
                                macroalgaePast[0] = RADIUS*cos(theta); //initial Macroalgae level
                                coralPast[0]      = RADIUS*sin(theta); //initial Coral level
                                for (int l=1;l<n;l++) //fills in the past times for y, z, v, and w
                                {
                                    macroalgaePast[l]=macroalgaePast[0];
                                    coralPast[l]=coralPast[0];
                                }
                                //fprintf(fp,"%f,%f,%f,%f,%f,%f\n",macroalgaePast[0],coralPast[0],1-macroalgaePast[0]-coralPast[0],v[0],w[0],1-v[0]-w[0]);
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
                                                                          steps,a,gamma,r,d,g,reefState,macroalgaePast,coralPast,dt,n,beta,tau,&fp,k,macroalgaePast[0],coralPast[0],theta);

#else

                                // Ignore the different threads. Just make one approximation in this one thread.
                                linear(steps,a,gamma,r,d,g,reefState,macroalgaePast,coralPast,dt,n,beta,tau,&fp,k,macroalgaePast[0],coralPast[0],theta);

#endif


#ifdef SHOW_PROGRESS
                                if(k%20 == 0)
                                    std::cout << "  Simulation number " << k << std::endl;
#endif

                            } // for(k<trials)
                        } // for(thetaStep)
                   } // for(betaStep)
            } // for(gStep)

            // Free up the allocated memory.
            free(macroalgaePast);
            free(coralPast);

        } // for(tauStep)


        fp.close();

#ifdef SHOW_PROGRESS
        std::cout << "all done" << std::endl;
#endif

        return b.exec();
}

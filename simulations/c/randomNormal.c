/* ************************************************************************
 randomNormal.c
 
 Program to provide an example of how to generate normally distributed
 random numbers with mean zero and a given standard deviation.
 
 ************************************************************************ */




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

int main(int argc,char **argv)
{

	/* Set the initial seed for the random number generator. */
	srand48(time(NULL));

	/* Now set up an example of how to use the routine. */
	double x[2];  /* allocate an array of length two to get the numbers. */

	/* Call the routine, and the results are in the vector x. Print out the results.*/
	normalDistRand(0.1,x);
	printf("%f,%f\n",x[0],x[1]);

	return(0);
}


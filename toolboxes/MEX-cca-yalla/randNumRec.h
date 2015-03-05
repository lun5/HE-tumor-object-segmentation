/* File randNumRec.h */

#ifndef TRUE
# define TRUE 1
#endif
#ifndef FALSE
# define FALSE 0
#endif

/* Seed the generator with a user provided number. */
void RanSeed( long number );

/* Seed the generator using the system clock. */
long RanTimeSeed();

/* Get the next pseudo-random number.
   The result is in the interval (0,1). */
float RanUniform();

/* The RanExponential() function returns the waiting time to the
   next event of a Poisson distribution which has
   Lambda events per second. */
float RanExponential( float Lambda );

/* Returns a normally distributed random number with zero
   mean and unit variance */
float RanNormal();

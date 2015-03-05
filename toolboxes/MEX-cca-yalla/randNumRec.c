/******************************************

 File:
  randNumRec.c	

 Description:
	Random number generator.

 Modified from TIPS package.
 (c) 1997, Thomas F. El-Maraghi

*******************************************/


#include <math.h>
#include <time.h>
#include "randNumRec.h"

/* The number of entries in the mSeedTable shuffle table */
#define SEED_TABLE_SIZE 32

static long mSeed;
static long mSeedTable[SEED_TABLE_SIZE];
static long mNextSeed;
static int   fExtraNormal;
static float ExtraNormal;


/* The following definitions are identical to those
   in the original implementation as found in
   "Numerical Recipies"  */
#define IA 16807
#define IM 2147483647
#define IQ 127773
#define IR 2836
#define AM (1.0/IM)
#define NDIV (1+(IM-1)/SEED_TABLE_SIZE)
#define EPS 1.2e-7
#define RandomMax (1.0-EPS)

/* Seed the generator with a user provided number. */
void RanSeed( long number )
{
        int i;
        long k;
	/* The Seed() function accepts a seed from the user
	   and uses it to fill in the mSeedTable shuffle table.  Values
	   in the mSeedTable shuffle table are then used by the
	   get function to generate the pseudo-random numbers. */
	mSeed = (number<0) ? number : -number;
	if( -mSeed < 1)
		mSeed = 1;
	else
		mSeed = -mSeed;
	for( i = SEED_TABLE_SIZE + 7 ; i >= 0 ; i-- ) {
		k = mSeed / IQ;
		mSeed = IA*(mSeed - k * IQ) - IR*k;
		if( mSeed < 0 )
			mSeed += IM;
		if( i < SEED_TABLE_SIZE )
			mSeedTable[i] = mSeed;
	}
	mNextSeed = mSeedTable[0];
}

/* Seed the generator using the system clock. */
long RanTimeSeed()
{
	/* A useful variant of the Seed() function which uses the system
	   clock to generate the initial seed. */
        long seed;
        seed = (long)time( 0 );
	RanSeed( seed );
	return seed;
}

/* Get the next pseudo-random number. 
   The result is in the interval (0,1). */
float RanUniform()
{
        long k;
        int Pick;
        float RetVal;
	if( mSeed <= 0 || !mNextSeed )
		RanTimeSeed();
	k = mSeed / IQ;
	mSeed = IA*(mSeed - k * IQ) - IR*k;
	if( mSeed < 0 )
		mSeed += IM;
	Pick = (int)( mNextSeed / NDIV );
	mNextSeed = mSeedTable[Pick];
	mSeedTable[Pick] = mSeed;
	RetVal = (float)AM * mNextSeed;
	if( RetVal > RandomMax )
		return (float)RandomMax;
	else
		return RetVal;
}

/* The RanExponential() function returns the waiting time to the
   next event of a Poisson distribution which has
   Lambda events per second. */
float RanExponential( float Lambda )
{
	float Dummy;
	do
	   Dummy = RanUniform();
	while( Dummy == 0.0 );
	return (float)(-log(Dummy)/Lambda);
}

/* Returns a normally distributed random number with zero
   mean and unit variance */
float RanNormal()
{
	float uniform1, uniform2, radius_square, factor;

	if( fExtraNormal )
	{
		/* if we have an extra deviate available, then use it */
		fExtraNormal = FALSE;
		return ExtraNormal;
	}
	
	/* Get two uniform deviates inside the unit circle */
	do {
		uniform1 = (float)(2.0*RanUniform()-1.0);
		uniform2 = (float)(2.0*RanUniform()-1.0);
		radius_square = uniform1*uniform1 + uniform2*uniform2;
	} while( radius_square >= 1.0 || radius_square == 0.0 );

	/* Compute 2 gaussian deviates, saving one for later */
	factor = (float)sqrt(-2.0*log(radius_square)/radius_square);
	ExtraNormal = uniform1*factor;
	fExtraNormal = TRUE;
	return uniform2*factor;
}


/********************************************************
 *
 * mextools.c
 *
 * C utilities for MEX programs.
 * 
 * Thomas El-Maraghi
 * September 2000
 *
 ********************************************************/


#include "mex.h"


/********************************************************
 * Convert a 2D mxArray to a 2D C array.  The converted
 * array is returned in 'out'.  The size of the array
 * is returned in 'nx' and 'ny' (i.e., columns, rows).
 ********************************************************/
void mex2float( const mxArray *in, float **out, int *nx, int *ny )
{
  int r, c;
  double *pin;
  float *pout;

  if( !in )
    mexErrMsgTxt( "Invalid mxArray" );
  pin = mxGetPr(in); /* get pointer to real part */
  if( !pin )
    mexErrMsgTxt( "Could not get pointer to mxArray data" );
  *ny = mxGetM(in);
  *nx = mxGetN(in);
  grabFloatMemory( &pout, *nx * *ny, "mxArray2Float" );
  for( c = 0; c < *nx; c++ )
    for( r = 0; r < *ny; r++, pin++ )
      *(pout + r * *nx + c) = (float)*pin;
  *out = pout;
}


/********************************************************
 * Convert a 2D mxArray to a 1D C array.  The converted
 * array is returned in 'out'.  The size of the array
 * is returned in 'n' (i.e., columns*rows).
 ********************************************************/
void mex2floatv( const mxArray *in, float **out, int *n )
{
  int i;
  double *pin;
  float *pout;

  if( !in )
    mexErrMsgTxt( "Invalid mxArray" );
  *n = mxGetM(in) * mxGetN(in);
  pin = mxGetPr(in);
  if( !pin )
    mexErrMsgTxt( "Could not get pointer to mxArray data" );
  grabFloatMemory( &pout, *n, "mxArray2FloatVector" );
  for( i = 0; i < *n; i++ )
    pout[i] = (float)pin[i];
  *out = pout;
}
 

/********************************************************
 * Convert a 2D C array to a 2D mxArray.  The converted
 * array is returned in 'out'.  The size of the array
 * to convert is 'nx' by 'ny' (i.e., columns, rows).
 ********************************************************/
void float2mex( const float *in, mxArray **out, int nx, int ny )
{
  int r, c;
  double *pout;
 
  if( !in )
    mexErrMsgTxt( "Invalid C array" );
  *out=mxCreateDoubleMatrix(ny,nx,mxREAL);
  if( !(*out) )
    mexErrMsgTxt( "Could not create mxArray" );
  pout=mxGetPr(*out);
  if( !pout )
    mexErrMsgTxt( "Could not get pointer to mxArray data" );
  for( c = 0; c < nx; c++ )
    for( r = 0; r < ny; r++, pout++ )
      *pout=(double)*(in + r * nx + c);  
}

 

/********************************************************
 *
 * mextools.h
 *
 * C utilities for MEX programs.
 * 
 * Thomas El-Maraghi
 * September 2000
 *
 ********************************************************/

#ifndef _mextools_h_
#define _mextools_h_


#include "mex.h"


/********************************************************
 * Convert a 2D mxArray to a 2D C array.  The converted
 * array is returned in 'out'.  The size of the array
 * is returned in 'nx' and 'ny' (i.e., columns, rows).
 ********************************************************/
void mex2float( const mxArray *in, float **out, int *nx, int *ny );


/********************************************************
 * Convert a 2D mxArray to a 1D C array.  The converted
 * array is returned in 'out'.  The size of the array
 * is returned in 'n' (i.e., columns*rows).
 ********************************************************/
void mex2floatv( const mxArray *in, float **out, int *n );
 

/********************************************************
 * Convert a 2D C array to a 2D mxArray.  The converted
 * array is returned in 'out'.  The size of the array
 * to convert is 'nx' by 'ny' (i.e., columns, rows).
 ********************************************************/
void float2mex( const float *in, mxArray **out, int nx, int ny );

 
#endif

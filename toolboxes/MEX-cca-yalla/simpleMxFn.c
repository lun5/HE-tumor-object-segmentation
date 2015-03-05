#include "mextools.h"
#include "utils.h"

void simple (float *image, int nx, int ny)
{
  int k;

  /* multiply every element by 2 */
  for (k = 0; k < nx*ny; k++) {
    image[k] *= 2.0;
  }
}

/* MEX interface function */
/* you CANNOT change the name of this routine */
/* this is the entry point */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  int nx, ny;

  int *rowPtr, *colInfo;
  double *data;
  int rows, cols, elements;

  /* check for the required input arguments */
  if( nrhs != 1 )
    mexErrMsgTxt( "One Input expected " );

  /* Check data type of input argument  */
  if (!(mxIsSparse(prhs[0]))) {
        mexErrMsgTxt("Input argument must be a sparse matrix.");
    }
 
  /* Get the size and pointers to input data */
  rows  = mxGetM(prhs[0]); /* number of rows */
  cols  = mxGetN(prhs[0]); /* number of cols */
 
  data  = mxGetPr(prhs[0]); /* pointer to real data */
  rowPtr = mxGetIr(prhs[0]); /* pointer to row info */
  colInfo = mxGetJc(prhs[0]); /* pointer to col info */

  /* number of non-zero elements */
  elements=mxGetNumberOfElements(prhs[0]);
   
  mexPrintf("#of rows cols: (%d %d)\n",rows,cols);
  mexPrintf("#of non-zeros: (%d)\n",colInfo[cols]);

  for (nx=0; nx < cols; nx++) {
    int noElemCol;
    noElemCol = colInfo[nx+1] - colInfo[nx];
    if (noElemCol > 0) {
      for (ny=colInfo[nx]; (ny < colInfo[nx+1])&&(rowPtr[ny] < nx); ny++) {
	/*for (ny=colInfo[nx]; (ny < colInfo[nx+1]); ny++) {*/
	mexPrintf("%d %d %f \n",rowPtr[ny]+1,nx+1,data[ny]);
      }
    }
  }

}

void mexFunction_old( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  float *image, *tmp;
  int nx, ny;

  /* check for the required input arguments */
  if( nrhs < 1 )
    mexErrMsgTxt( "Input expected: <floating-point-Array>" );

  /* get the image to be filtered */
  mex2float( prhs[0], &image, &nx, &ny);
  /* perform a simple operation */
  simple(image,nx,ny);

  /* return the image to matlab */
  float2mex(image, &plhs[0], nx,  ny ); 

  /* free up memory */
  utilFree( (void **)&image );

}

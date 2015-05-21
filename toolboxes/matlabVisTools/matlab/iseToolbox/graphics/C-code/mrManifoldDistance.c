/*
   mrManifoldDistance

   AUTHOR:  Engel, Wandell
   DATE:    Nov., 1994
   PURPOSE:
   This code is used to create a mex-file for matlab for the routine
   mrManifoldDistance().
   
   The input is an array of sample point coordinates that should form
   a connected manifold in three-space.

   The point of the routine is to compute the distances between a point 
   in three-space and a set of other points in three-space.  The distance
   is measured through the connected space of points. 

   DESCRIPTION:

    dist = mrManifoldDistance(grayM,iSize,numSlices,startPt,[noVal],[radius])

   ARGUMENTS:
    grayM:  A volume of binary data indicating where the gray matter is.
    iSize:  The size of an image slice through the gray matter
    numSlices: The number of image slices
    startPt:   3d coordinates defining where to start the flood fill

   OPTIONAL ARGUMENTS:
    dimdist:Array of y,x and z separations between points.
    noVal:  The value returned for unreached locations (default 0)
    radius: The max distance to flood out
     (default 0 == flood as far as can)

   RETURNS:
    dist:  distances to each point from startPt -- same size as grayM
    nPntsReached:  The number of points reached from the start point.
   
*/

#include "mex.h"
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>
#include "cityheap.h"

#define POINT_UNREACHED -1
#define XCOORD 0
#define YCOORD 1
#define ZCOORD 2
#define DIAG  -1
#define HEAP_SIZE 	200000
#define POINT_REACHED 	0

#define XSEPARATION .9375;
#define YSEPARATION .9375;
#define ZSEPARATION .70

/* Routines used here */
int floodfill(double *dist, int *grayM, int *iSize, int numSlices, 
	       int *startPt, double *dimdist, int radius);
long mr3d21d(int *samp, int *iSize);
active_elem *make_active_elem(int *samp, float dist, int prev);


/* The utility routine mr3d21d converts a 3d coordinate, samp, into an
   integer index.  This integer can be used to index into the
   volume array that contains the data
*/
long mr3d21d(int *samp, int *iSize)
{      

  return( (long) samp[1] + 
	         (samp[0] - 1)*iSize[0] +
	         (samp[2] - 1)*iSize[0]*iSize[1] - 1); 

}

/* Allocate memory and set the values of an active element along
   the search path
*/
active_elem *make_active_elem(int *samp, float dist, int pdir)
{
  active_elem *elt;

  elt = (active_elem *)mxCalloc(1,sizeof(active_elem));
  if(!elt) {
    mexErrMsgTxt("make_active_elem(): unable to allocate memory");
  }

  (elt->coord)[0] = samp[0];
  (elt->coord)[1] = samp[1];
  (elt->coord)[2] = samp[2];

  elt->prevdir = pdir;
  elt->dist = dist;

  return(elt);
}

/* Compute the distance from the startPoint to each of the other
   gray-matter points in the volume.  The distances are returned in
   the vector double *dist.
*/

int floodfill(double *dist, int *grayM, int *iSize, int numSlices, 
	       int *startPt, double *dimdist, int radius)
{
  Heap *heep;
  active_elem *the_path;
  float nudist;
  long idx;
  int tmp[3];
  int i,j,k;
  int thecoord, incr, pdir, count;
  long nPixels = 0;

  nPixels = numSlices*iSize[0]*iSize[1];

  /* What is the heap used for?  Should this size be hardcoded?  */
  /* Initialize the heap */
  heep = make_empty_heap(HEAP_SIZE);

  /* Initialize distances and stuff at the starting point */
  /* Mark the gray matter index at the start point so we don't return here */
  idx = mr3d21d(startPt, iSize);
  dist[idx] = 0;
  grayM[idx] = POINT_REACHED;
  heap_insert(heep,make_active_elem(startPt,0,DIAG));

  /* We are going to create a wavefront that emanates from the start point. */
  /* Each point in the gray matter will be assigned a "the_path" structure */
  /* that has associated with it a distance, coordinate, and previous */
  /* direction.  I am a little worried whether this is the proper algorithm */
  /* when the dimdists are not equal.  -- BW */

  /* Each element of the heap records: its coordinates, the shortest */
  /* distance to that coordinate from the start point, and the previous */
  /* direction along the path that got there from the start point */
  /* You visit each voxel once.  Once you have gotten there, you know the */
  /* the distance there */

  /* Let's count the number of points we reached */
  count = 0;

  /* The while conditional sets the_path to the maximal element in the heap */
  /* and removes that element from the heap. When there are no more */
  /* terms on the heap, the_path is NULL and the while() loop exits */
  while (( the_path = heap_extract_max(heep) ) != NULL) {

    /* If we are within the radius, or radius is 0 search */
    if((radius == 0) || (the_path->dist < radius)) {

      /* Check all of the unvisited neighbors that are in the gray matter */
      /* We calculate distances from the start point, put them on the heap */
      /* and we mark them as visited.   */
      /* We check in six directions, + and - in x,y,z.  But, the distance */
      /* is calculated using a (clever?) hack to get the corners right */
      for(thecoord = XCOORD; thecoord <= ZCOORD; thecoord++)
	for(incr = -1; incr < 2; incr+=2) {

	  /* Put the coordinate into tmp */
	  tmp[XCOORD] = ((the_path->coord)[XCOORD]);
	  tmp[YCOORD] = ((the_path->coord)[YCOORD]);
	  tmp[ZCOORD] = ((the_path->coord)[ZCOORD]);

	  /* Check an adjacent coordinate */
	  tmp[thecoord] += incr;

	  /* If the adjacent coordinate is in range, do this stuff */
	  if((tmp[ZCOORD] >= 1) && (tmp[ZCOORD] <= numSlices) &&
	     (tmp[YCOORD] >= 1) && (tmp[YCOORD] <= iSize[0]) &&
	     (tmp[XCOORD] >= 1) && (tmp[XCOORD] <= iSize[1])) {

	    /* Find the index for this adjacent coordinate */
	    idx = mr3d21d(tmp,iSize);
	    if (idx > nPixels) {
	      mexPrintf("idx seems too big: %d\n",idx);
	    }	    

	    /* If this coordinate is in the gray-matter */
	    /* fill in its the_path parameters properly */
	    if(grayM[idx] > 0) {
	      count++;

	      /* Don't come here again */
	      grayM[idx] = POINT_REACHED;

	      /* Copy the previous direction to this shorter name  */
	      /* so the code looks nicer */
	      pdir = the_path->prevdir;

	      /* Now, assign a distance according to whether we are on */
	      /* a diagonal or not.  We compute a distance */
	      /* to the current location from the start point. */

	      /* If we are moving along the same direction, or */
	      /* the previous step was diagonal, increment distance this way */
	      if((thecoord == pdir) || (pdir == DIAG)) {
		nudist = (the_path->dist)+dimdist[thecoord];
		pdir = thecoord;
	      }
	      /* If we are moving along a new diagonal step, then */
	      /* increment the distance this way */
	      else  {
		nudist = (the_path->dist)-dimdist[pdir];
		nudist = nudist +
		  sqrt(dimdist[pdir]*dimdist[pdir] +
		       dimdist[thecoord]*dimdist[thecoord]);
		pdir = DIAG;
	      }

	      /* Assign the new distance to the_path and to dist */
	      /* dist[idx] = (double) (the_path->dist); */
	      dist[idx] = nudist;

	      /* Add the results from this voxel to the heap */
	      /* so we can check its neighbors later */
	      heap_insert(heep,make_active_elem(tmp,nudist,pdir));
	    }
	  }
	}
    }
    mxFree(the_path);

  }
  free_heap(heep);

  return(count);

}

void mexFunction(int nlhs,	/* number of arguments on lhs */
		 Matrix	*plhs[], /* Matrices on lhs      */
		 int nrhs,	/* no. of mat on rhs    */
		 Matrix	*prhs[]	/* Matrices on rhs      */
		 )
{
  double *dist, *tmp, *nPointsReached;
  int iSize[2],  numSlices, startPt[3], *grayM;
  long i,nPixels;
  int nPnts, n, radius;
  double nuldist, dimdist[3];
  int count;

  /* Check for proper number of arguments */

  if (nrhs == 0) {
    printf("dist = mrManifoldDistance(grayM,iSize,numSlices,startPt,[noVal],[radius]) \n");
  }
 
  else {

    /* Create space for return arguments on the lhs */

    /* The size of dist is equal to the size of grayM */
    plhs[0]=mxCreateFull(mxGetM(prhs[0]),mxGetN(prhs[0]),REAL);
    dist=mxGetPr(plhs[0]);

    /* Interpret the input arguments on the rhs */

    if (nrhs < 4) {
      mexErrMsgTxt("mrManifoldDistance: At least four arguments are needed.");
    }
    
    /* 1.  Extract the data in grayM and converted to ints */
    tmp=mxGetPr(prhs[0]);
    grayM = (int *) mxCalloc(mxGetN(prhs[0])*mxGetM(prhs[0]),sizeof(int));

    /* I don't understand what is going on here.  If grayM is binary, */
    /* why are we doing this work? */
    for (i = 0; i < mxGetN(prhs[0])*mxGetM(prhs[0]); i++) {
      grayM[i] = floor( 0.5 + tmp[i]);
    }
    nPnts = mxGetN(prhs[0])*mxGetM(prhs[0]);
    fflush(stderr);

    /* 2.  Extract the image size */
    tmp=mxGetPr(prhs[1]);
    iSize[0] = floor( 0.5 + tmp[0]);
    iSize[1] = floor( 0.5 + tmp[1]);

    /* 3.  Extract the number of slices */
    tmp=mxGetPr(prhs[2]);
    numSlices = floor( 0.5 + tmp[0]);

    nPixels = iSize[0]*iSize[1]*numSlices;
    if (nPixels != mxGetM(prhs[0])*mxGetN(prhs[0])) {
      mexErrMsgTxt("mrManifoldDistance: grayM not the right size.");
    }

    /* 4.  Extract the start point coordinates */
    tmp=mxGetPr(prhs[3]);
    startPt[0] = floor( 0.5 + tmp[0]);
    startPt[1] = floor( 0.5 + tmp[1]);
    startPt[2] = floor( 0.5 + tmp[2]);

    /* Analyze the optional arguments */

    /* 5.  Figure out the distances in and between planes */
    if (nrhs < 5) {
      dimdist[0] = YSEPARATION;
      dimdist[1] = XSEPARATION;
      dimdist[2] = ZSEPARATION;
    }
    else {
      tmp = mxGetPr(prhs[4]);
      dimdist[0] = tmp[0];/* row, col, and plane separations */
      dimdist[1] = tmp[1];/* which are y,x and z */
      dimdist[2] = tmp[2];
    }

    /* 6.  Choose the default value when points are unreached */
    if (nrhs < 6) {
      nuldist = POINT_UNREACHED;
    }
    else {
      tmp=mxGetPr(prhs[5]);
      nuldist = tmp[0];
    }

    /* Initialize the distance values in the matrix to the default */
    for (i = 0; i < nPnts; i++) {
      dist[i] = nuldist;
    }

    /*7.   Read the radius */
    if (nrhs >= 7) {
      tmp=mxGetPr(prhs[6]);
      radius = floor( 0.5 + tmp[0]);
    }
    else {
      radius = 0;
    }
    /* mexPrintf("radius:  %d\n",radius); */

    /* Call the floodfill to measure the distances */
    count = floodfill(dist, grayM, iSize, numSlices, startPt, dimdist, radius);
    /* mexPrintf("count:  %d\n",count); */

    if(nlhs > 1) {
      plhs[1]=mxCreateFull(1,1,REAL);
      nPointsReached = mxGetPr(plhs[1]);
      *nPointsReached = (double) count;
     /* mexPrintf("count:  %f\n", *nPointsReached);*/
    }

    /* free the large gray-matter vector */
    mxFree(grayM);
  }      

}



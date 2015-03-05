/* generate binary cca.mexglx */
/*  The routines are from ADJ's connected components toolbox. */
/*  CC made modifications to handle affinity matrices. */

#include "mex.h"
#include "utils.h"
#include "macros.h"
#include "imageFile-io.h"
#include "randNumRec.h"
#include "mextools.h"
#include "connCompLinks2.h"

/* Link directions in E, SW, S, SE order */
static int dx[NUMLINKS] = {1, -1, 0, 1};
static int dy[NUMLINKS] = {0, 1, 1, 1};

/* just so i remember */
/* # define  index(IX, IY, NX) ((IY)*(NX)+(IX)) */

extern void connectMarkov(ccLabelStruct *label[], float *target, 
			       float minTarget,int nX, int nY) ;
extern void connectSymmetric(ccLabelStruct *label[], float *target, 
			     float minTarget,int nX, int nY) ;
extern void connectSymmetricDbl(ccLabelStruct *label[], double *target, 
			     double minTarget,int nX, int nY) ;

void connectSymmetricSparse(ccLabelStruct *label[], 
			    mwIndex *rowPtr, mwIndex *colInfo, double *data, 
			    double minTarget, int nX, int nY) ;

extern void processMarkov();
extern void processLines();

extern void emptyLabels(ccLabelStruct *label[], int nx, int ny);
extern void printMarkovLabels(ccLabelStruct **label,int nX);

extern void analyze_sparse(const mxArray *array_ptr);


/*============================================================*/

void labels2mex(ccLabelStruct **label, mxArray **out, int nx)
{
  int r, c;
  double *pout;
  
  *out=mxCreateDoubleMatrix(nx,1,mxREAL);
  if( !(*out) )
    mexErrMsgTxt( "Could not create mxArray" );
  pout=mxGetPr(*out);
  if( !pout )
    mexErrMsgTxt( "Could not get pointer to mxArray data" );
  for( c = 0; c < nx; c++, pout++ )
    {
      if (label[c] != NULL)
	*pout=(double) label[c]->count;
      else
	*pout=(double) -1;
    }
}

/*============================================================*/
/* symmetric affinity matrix is double */
/*============================================================*/
void connectSymmetricDbl(ccLabelStruct *label[], double *target, 
			 double minTarget, int nX, int nY)
{
  int ix, iy, ixy, ixn, iyn, ixyn, kDir;
  ccLabelStruct *currentRoot, *nhbrLabel, *nhbrRoot;

  for (iy = 0; iy < nY; iy++) {
    ixy = index(iy,iy,nX);
    label[iy] = NULL;
    currentRoot =  NULL;
    /*mexPrintf("[%d %2.3f] \n",iy,target[ixy]);*/
    if (target[ixy] >= minTarget) {
      if (iy > 0) {
	/* loop over neighbors */
	for (ix = iy-1; ix >= 0; ix--) {
	  /* symmetric matrix. it does not matter how the */
	  /* following is written */
	  ixyn = index(ix,iy,nX);
	  /*ixyn = index(iy,ix,nX);*/
	  nhbrLabel = label[ix];
	  /*mexPrintf("\n (%d %d %2.3f) \n",iy,ix,target[ixyn]);*/
	  if ((nhbrLabel != NULL) && (target[ixyn] >= minTarget)) {
	    /* Connect this neighbour to current pixel */
	    nhbrRoot = findLabelRoot(nhbrLabel);
	    
	    if (currentRoot == NULL) { /* Current pixel unlabelled */
	      label[iy] = nhbrRoot;  /* Use label of connected neighbour */
	      currentRoot = nhbrRoot;
	      currentRoot->size += 1;
	      currentRoot->weight = max(currentRoot->weight, target[ixy]);
	    } else if (nhbrRoot != currentRoot)
	      /* Merge label trees */
	      currentRoot = mergeLabelRoots(nhbrRoot, currentRoot);
	  }
	} /* End: loop over neighbours */
      }
      if (currentRoot == NULL) { /* No labelled neighbours found */
	/* Generate new label for this pixel */
	/*mexPrintf("\n no neighbors %d \n",iy);*/
	label[iy] = newLabel();
	label[iy]->size = 1;
	label[iy]->weight = target[ixy];
      }
    }/* target[ixy] >= minTarget */
  } /* End: go over the diagonal */
  /* Completed labelling. */

  /* Reset all the labels to refer to only the root node */
  flattenLabels(label, nX, 1);
 
} /* End: connectSymmetricDbl */


/*============================================================*/
/* symmetric affinity matrix is sparse and double */
/*============================================================*/
void connectSymmetricSparse(ccLabelStruct *label[], 
			    mwIndex *rowPtr, mwIndex *colInfo, double *data, 
			    double minTarget, int cols, int rows)
{/* begin connectSymmetricSparse */
  ccLabelStruct *currentRoot, *nhbrLabel, *nhbrRoot;
  int nx, ny;
  int debugId;

  label[0] = newLabel();
  label[0]->size = 1;
  label[0]->weight = data[0];

  debugId = 4063;
  /*
  for (nx=0; nx < cols; nx++) {
    mexPrintf("rowBou %d %d \n",nx,colInfo[nx+1]);    
  }
  */

  for (nx=0; nx < cols; nx++) {
    int noElemCol, ok;
    double tmpdata = 0;

    noElemCol = colInfo[nx+1] - colInfo[nx];
    /*mexPrintf("noElemCol %d %d %d \n",noElemCol,colInfo[nx+1],colInfo[nx]);*/
    label[nx] = NULL;
    currentRoot =  NULL;

    /* why check the diagonal elements? */
    ok = 1;

    /*
    for (ny=colInfo[nx]; (ny < colInfo[nx+1])&&(!ok); ny++) {
      if (rowPtr[ny] == nx) {
	ok = (data[ny] >= minTarget) ? 1 : 0;
	tmpdata = data[ny];
      }
    }
    mexPrintf("%d [%d %2.3f] \n",ok,nx,tmpdata);
    */

    if ((nx+1)==debugId) {
      for (ny=colInfo[nx]; (ny < colInfo[nx+1]); ny++) {
	mexPrintf("%d [%d %2.3f] \n",rowPtr[ny]+1,ny,data[ny]);
      }
    }


    if (ok) { 
      /* loop over neighbors */
      if (noElemCol) {

	/* making use of symmetry by reading only the upper-triangular matrix */
	/*for (ny=colInfo[nx]; (ny < colInfo[nx+1])&&(rowPtr[ny] < nx); ny++) {*/
	for (ny=colInfo[nx]; (ny < colInfo[nx+1]); ny++) {

	  nhbrLabel = label[rowPtr[ny]];
	  

	  if ((nx+1)==debugId) {
	    /*mexPrintf("label: %f (%d %d %f) \n",label[rowPtr[ny]]->weight, 
	      rowPtr[ny]+1,nx+1,data[ny]);*/
	    mexPrintf("label: (%d %d %f) \n",rowPtr[ny]+1,nx+1,data[ny]);
	  }


	  if ((nhbrLabel != NULL) && (data[ny] >= minTarget)) {

	    if ((nx+1)==debugId) {
	      /*mexPrintf("label: %f (%d %d %f) \n",label[rowPtr[ny]]->weight, 
		rowPtr[ny]+1,nx+1,data[ny]);*/
	      mexPrintf("nhbrlabel is not NULL: (%d %d %f) \n",rowPtr[ny]+1,nx+1,data[ny]);
	    }

	    /* Connect this neighbour to current pixel */
	    nhbrRoot = findLabelRoot(nhbrLabel);
	    
	    if (currentRoot == NULL) { /* Current pixel unlabelled */
	      label[nx] = nhbrRoot;  /* Use label of connected neighbour */
	      currentRoot = nhbrRoot;
	      currentRoot->size += 1;
	      currentRoot->weight = max(currentRoot->weight, data[ny]);
	    } else if (nhbrRoot != currentRoot)
	      /* Merge label trees */
	      currentRoot = mergeLabelRoots(nhbrRoot, currentRoot);
	  }
	} /* End: loop over neighbours */
      }
      if (currentRoot == NULL) { /* No labelled neighbours found */
	/* Generate new label for this pixel */
	/*mexPrintf("\n no neighbors %d \n",iy);*/

	if ((nx+1)==debugId) {
	  /*mexPrintf("label: %f (%d %d %f) \n",label[rowPtr[ny]]->weight, 
	    rowPtr[ny]+1,nx+1,data[ny]);*/
	  mexPrintf("generating new label: (%d %d %f) \n",rowPtr[ny]+1,nx+1,data[ny]);
	}

	label[nx] = newLabel();
	label[nx]->size = 1;
	label[nx]->weight = tmpdata;
      }
    }/* target[ixy] >= minTarget */
  } /* End: go over the diagonal */
  /* Completed labelling. */

  /* Reset all the labels to refer to only the root node */
  flattenLabels(label, cols, 1);
 
} /* End: connectSymmetricSparse */


void connectSymmetricSparse_old(ccLabelStruct *label[], 
			    int *rowPtr, int *colInfo, double *data, 
			    double minTarget, int cols, int rows)
{/* begin connectSymmetricSparse */
  ccLabelStruct *currentRoot, *nhbrLabel, *nhbrRoot;
  int nx, ny;


  label[0] = newLabel();
  label[0]->size = 1;
  label[0]->weight = data[0];

  for (nx=1; nx < cols; nx++) {
    int noElemCol, ok;
    double tmpdata;

    noElemCol = colInfo[nx+1] - colInfo[nx];
    mexPrintf("noElemCol %d %d %d \n",noElemCol,colInfo[nx+1],colInfo[nx]);
    label[nx] = NULL;
    currentRoot =  NULL;

    ok = 0;
    for (ny=colInfo[nx]; (ny < colInfo[nx+1])&&(!ok); ny++) {
      if (rowPtr[ny] == nx) {
	ok = (data[ny] >= minTarget) ? 1 : 0;
	tmpdata = data[ny];
      }
    }
    mexPrintf("%d [%d %2.3f] \n",ok,nx,tmpdata);
    if (ok) { 
      /* loop over neighbors */
      if (noElemCol > 1) {
	for (ny=colInfo[nx]; (ny < colInfo[nx+1])&&(rowPtr[ny] < nx); ny++) {
	  nhbrLabel = label[rowPtr[ny]];
	  
	  /*mexPrintf("label: %f (%d %d %f) \n",label[rowPtr[ny]]->weight, 
	    rowPtr[ny]+1,nx+1,data[ny]);*/

	  if ((nhbrLabel != NULL) && (data[ny] >= minTarget)) {
	    /* Connect this neighbour to current pixel */
	    nhbrRoot = findLabelRoot(nhbrLabel);
	    
	    if (currentRoot == NULL) { /* Current pixel unlabelled */
	      label[nx] = nhbrRoot;  /* Use label of connected neighbour */
	      currentRoot = nhbrRoot;
	      currentRoot->size += 1;
	      currentRoot->weight = max(currentRoot->weight, data[ny]);
	    } else if (nhbrRoot != currentRoot)
	      /* Merge label trees */
	      currentRoot = mergeLabelRoots(nhbrRoot, currentRoot);
	  }
	} /* End: loop over neighbours */
      }
      if (currentRoot == NULL) { /* No labelled neighbours found */
	/* Generate new label for this pixel */
	/*mexPrintf("\n no neighbors %d \n",iy);*/
	label[nx] = newLabel();
	label[nx]->size = 1;
	label[nx]->weight = tmpdata;
      }
    }/* target[ixy] >= minTarget */
  } /* End: go over the diagonal */
  /* Completed labelling. */

  /* Reset all the labels to refer to only the root node */
  flattenLabels(label, cols, 1);
 
} /* End: connectSymmetricSparse */


/*============================================================*/
/* ability to handle both sparse and full affinity matrix     */ 
/*============================================================*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  float *image, *tmp;
  int nX, nY;
  ccLabelStruct **label;
  int ixy;
  int numComp, t;
  double minTarget; /* 0.001 = boxes02-add.pgm */
  double *target;
  float debug = 0;

  /* for sparse array */
  mwIndex *rowPtr, *colInfo;
  double *data;
  int rows, cols;

  int sparseFlag = 0;

  RanTimeSeed;

  /*mexPrintf("nlhs:%d nrhs: %d\n",nlhs,nrhs);*/

  /* check for the required input arguments */
  if( nrhs < 2 )
    mexErrMsgTxt( "Input expected: <symmetric-positive-array-for-cca> <affty-threshold> <1-for-debug>" );

  /* get the symmetric matrix */
  /*mex2float( prhs[0], &target, &nX, &nY);*/
  nX = mxGetN(prhs[0]); /* cols */ 
  cols = nX;

  nY = mxGetM(prhs[0]); /* rows */
  rows = nY;
  
  mexPrintf("size: %d %d\n",nX,nY);

  if (mxIsSparse(prhs[0])) {
    sparseFlag = 1;
    mexPrintf("Affy matrix Sparse? : %d\n",sparseFlag);
  }

  if (sparseFlag) {
    mwIndex *jc;

    data    = mxGetPr(prhs[0]); /* pointer to real data */
    rowPtr  = mxGetIr(prhs[0]); /* pointer to row info */
    colInfo = mxGetJc(prhs[0]); /* pointer to col info */
    jc = mxGetJc(prhs[0]);

    /*
    for (ixy=0; ixy < mxGetN(prhs[0]); ixy++) {
      mexPrintf("rowBouMain %d (%d %d) \n",ixy,jc[ixy+1],colInfo[ixy+1]);    
    }

    analyze_sparse(prhs[0]);
    */
  } else {
    target = (double *)mxGetPr(prhs[0]);
  }


  /* get the threshold */
  minTarget = (double) mxGetScalar(prhs[1]);

  /* check for debug flag */
  if (nrhs == 3) {
    debug = mxGetScalar(prhs[2]);
  }

  if (debug)
    mexPrintf("Affty threshold: %2.3e\n",minTarget);

  /* set up empty labels */
  grabByteMemory((char **) &label,sizeof(ccLabelStruct *) * nX,"label");
  for (ixy=0;ixy < nX; ixy++)
    label[ixy] = NULL;

  if (debug) {
    mexPrintf("Computing connected components....");
  }

  if (sparseFlag) {
    connectSymmetricSparse(label,rowPtr,colInfo,data, minTarget, nX, nY);
  } else {
    connectSymmetricDbl(label, target, minTarget, nX, nY);
  }

  if (debug) {
    mexPrintf("done\n");
  }

  /* yalla edition */
  /*
  if (debug) {
    mexPrintf("Pruning small or weak connected components....");
  }
  pruneLabelComp(label, nX, 1, 1, 0);
  */

  if (debug) {
    mexPrintf("done\n");
  }

  /* component count */
  numComp = countLabelComp(label, nX, 1);
  if (debug) {
    mexPrintf("Found %d components.\n\n", numComp);
  }

  /*
  if (numComp > 0)
    printMarkovLabels(label,nX);
  */

  /* return the labels to matlab */
  if (nlhs > 0)
    labels2mex(label, &plhs[0], nX ); 
  if (nlhs > 1)
    {
      double *pout;
      plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
      pout = mxGetPr(plhs[1]);
      *pout = (double) numComp;
    }
  /* free up memory. */
  /* target is a pointer to prhs[0]. so no need to free. */
  /*utilFree( (void **)&target );*/
  utilFree((void **)label);
  /*
  for (ixy = 0; ixy < nX; ixy++)
    free(label[ixy]);
  */

}

/*============================================================*/
/* MEX interface function: replica of processLines()          */
/* do not declare new arrays but instead work with            */
/* arrays that have been already declared in MATLAB           */
/* also the tau threshold is now read as double.              */ 
/*============================================================*/
void mexFunction_double(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  float *image, *tmp;
  int nX, nY;
  ccLabelStruct **label;
  int ixy;
  int numComp, t;
  double minTarget; /* 0.001 = boxes02-add.pgm */
  double *target;
  float debug = 0;

  RanTimeSeed;

  /*mexPrintf("nlhs:%d nrhs: %d\n",nlhs,nrhs);*/

  /* check for the required input arguments */
  if( nrhs < 2 )
    mexErrMsgTxt( "Input expected: <symmetric-positive-array-for-cca> <affty-threshold> <1-for-debug>" );

  /* get the symmetric matrix */
  /*mex2float( prhs[0], &target, &nX, &nY);*/
  nX = mxGetN(prhs[0]); /* cols */ 
  nY = mxGetM(prhs[0]); /* rows */
  target = (double *)mxGetPr( prhs[0] );

  /*mexPrintf("size: %d %d\n",nX,nY);*/

  /* get the threshold */
  minTarget = (double) mxGetScalar(prhs[1]);
  if (nrhs == 3) {
    debug = mxGetScalar(prhs[2]);
  }

  if (debug)
    mexPrintf("Affty threshold: %2.3e\n",minTarget);

  /* set up empty labels */
  grabByteMemory((char **) &label,sizeof(ccLabelStruct *) * nX,"label");
  for (ixy=0;ixy < nX; ixy++)
    label[ixy] = NULL;

  if (debug) {
    mexPrintf("Computing connected components....");
  }
  connectSymmetricDbl(label, target, minTarget, nX, nY);
  if (debug) {
    mexPrintf("done\n");
  }

  if (debug) {
    mexPrintf("Pruning small or weak connected components....");
  }
  pruneLabelComp(label, nX, 1, 1, 0);
  if (debug) {
    mexPrintf("done\n");
  }

  /* component count */
  numComp = countLabelComp(label, nX, 1);
  if (debug) {
    mexPrintf("Found %d components.\n\n", numComp);
  }

  /*
  if (numComp > 0)
    printMarkovLabels(label,nX);
  */

  /* return the labels to matlab */
  if (nlhs > 0)
    labels2mex(label, &plhs[0], nX ); 
  if (nlhs > 1)
    {
      double *pout;
      plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
      pout = mxGetPr(plhs[1]);
      *pout = (double) numComp;
    }
  /* free up memory. */
  /* target is a pointer to prhs[0]. so no need to free. */
  /*utilFree( (void **)&target );*/
  utilFree((void **)label);
  /*
  for (ixy = 0; ixy < nX; ixy++)
    free(label[ixy]);
  */

}

/*============================================================*/
/* MEX interface function */
/* replica of processLines */
/*============================================================*/
void mexFunction_float( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  float *image, *tmp;
  int nX, nY;
  ccLabelStruct **label;
  int ixy;
  int numComp, t;
  float minTarget; /* 0.001 = boxes02-add.pgm */
  float *target;
  float debug = 0;

  RanTimeSeed;

  /*mexPrintf("nlhs:%d nrhs: %d\n",nlhs,nrhs);*/

  /* check for the required input arguments */
  if( nrhs < 2 )
    mexErrMsgTxt( "Input expected: <symmetric-positive-array-for-cca> <affty-threshold> <1-for-debug>" );

  /* get the symmetric matrix */
  mex2float( prhs[0], &target, &nX, &nY);
  /*mexPrintf("size: %d %d\n",nX,nY);*/

  /* get the threshold */
  minTarget = mxGetScalar(prhs[1]);
  if (nrhs == 3) {
    debug = mxGetScalar(prhs[2]);
  }

  if (debug)
    mexPrintf("Affty threshold: %2.3e\n",minTarget);

  /* set up empty labels */
  grabByteMemory((char **) &label,sizeof(ccLabelStruct *) * nX,"label");
  for (ixy=0;ixy < nX; ixy++)
    label[ixy] = NULL;

  if (debug) {
    mexPrintf("Computing connected components....");
  }
  connectSymmetric(label, target, minTarget, nX, nY);
  if (debug) {
    mexPrintf("done\n");
  }

  if (debug) {
    mexPrintf("Pruning small or weak connected components....");
  }
  pruneLabelComp(label, nX, 1, 1, 0);
  if (debug) {
    mexPrintf("done\n");
  }

  /* component count */
  numComp = countLabelComp(label, nX, 1);
  if (debug) {
    mexPrintf("Found %d components.\n\n", numComp);
  }

  /*
  if (numComp > 0)
    printMarkovLabels(label,nX);
  */

  /* return the labels to matlab */
  if (nlhs > 0)
    labels2mex(label, &plhs[0], nX ); 
  if (nlhs > 1)
    {
      double *pout;
      plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
      pout = mxGetPr(plhs[1]);
      *pout = (double) numComp;
    }
  /* free up memory */
  utilFree( (void **)&target );
  utilFree((void **)label);
  /*
  for (ixy = 0; ixy < nX; ixy++)
    free(label[ixy]);
  */

}

/*============================================================*/
/* following works for symmetric matrix. targe array is float */
/* changed it to double in the code above: connectSymmetric() */
/*============================================================*/
void connectSymmetric(ccLabelStruct *label[], float *target, 
			  float minTarget,int nX, int nY)
{
  int ix, iy, ixy, ixn, iyn, ixyn, kDir;
  ccLabelStruct *currentRoot, *nhbrLabel, *nhbrRoot;

  for (iy = 0; iy < nY; iy++) {
    ixy = index(iy,iy,nX);
    currentRoot = label[iy] = NULL;
    /*fprintf(stderr,"[%d %2.3f] \n",iy,target[ixy]);*/
    if (target[ixy] >= minTarget) {
      if (iy > 0) {
	/* loop over neighbors */
	for (ix = iy-1; ix >= 0; ix--) {
	  /* symmetric matrix. it does not matter how the */
	  /* following is written */
	  ixyn = index(ix,iy,nX);
	  /*ixyn = index(iy,ix,nX);*/
	  nhbrLabel = label[ix];
	  /*fprintf(stderr,"\n (%d %d %2.3f) \n",ix,iy,target[ixyn]);*/
	  if ((nhbrLabel != NULL) && (target[ixyn] >= minTarget)) {
	    /* Connect this neighbour to current pixel */
	    nhbrRoot = findLabelRoot(nhbrLabel);
	    
	    if (currentRoot == NULL) { /* Current pixel unlabelled */
	      label[iy] = nhbrRoot;  /* Use label of connected neighbour */
	      currentRoot = nhbrRoot;
	      currentRoot->size += 1;
	      currentRoot->weight = max(currentRoot->weight, target[ixy]);
	    } else if (nhbrRoot != currentRoot)
	      /* Merge label trees */
	      currentRoot = mergeLabelRoots(nhbrRoot, currentRoot);
	  }
	} /* End: loop over neighbours */
      }
      if (currentRoot == NULL) { /* No labelled neighbours found */
	/* Generate new label for this pixel */
	/*fprintf(stderr,"\n no neighbors %d \n",iy);*/
	label[iy] = newLabel();
	label[iy]->size = 1;
	label[iy]->weight = target[ixy];
      }
    }/* target[ixy] >= minTarget */
  } /* End: go over the diagonal */
  /* Completed labelling. */

  /* Reset all the labels to refer to only the root node */
  flattenLabels(label, nX, 1);
 
} /* End: connectSymmetric */

void analyze_sparse(const mxArray *array_ptr)
{
  double  *pr, *pi;
  mwIndex  *ir, *jc;
  mwSize      col, total=0;
  mwIndex   starting_row_index, stopping_row_index, current_row_index;
  mwSize      n, ixy;
  
  /* Get the starting positions of all four data arrays. */ 
  pr = mxGetPr(array_ptr); /* real elements */
  pi = mxGetPi(array_ptr); /* imaginary elements */
  ir = mxGetIr(array_ptr); /* row pointer */
  jc = mxGetJc(array_ptr); /* col pointer */
  /*
  for (ixy=0; ixy < mxGetN(array_ptr); ixy++) {
    mexPrintf("rowBou %d %d \n",ixy,jc[ixy+1]);    
  }
  */

  /* Display the nonzero elements of the sparse array. */ 
  n = mxGetN(array_ptr);
  for (col=0; col<n; col++)  { 
    starting_row_index = jc[col]; 
    stopping_row_index = jc[col+1]; 
    mexPrintf("%d %d\n",starting_row_index,stopping_row_index);
    if (starting_row_index == stopping_row_index)
      continue;
    else {
      for (current_row_index = starting_row_index; 
	   current_row_index < stopping_row_index; 
	   current_row_index++)  {
	if (mxIsComplex(array_ptr))  {
	  mexPrintf("\t(%"FMT_SIZE_T"u,%"FMT_SIZE_T"u) = %g+%g i\n", 
                    ir[current_row_index]+1, 
                    col+1, pr[total], pi[total]);
	  total++;
	} else {
	  mexPrintf("\t(%"FMT_SIZE_T"u,%"FMT_SIZE_T"u) = %g\n", 
                    ir[current_row_index]+1, 
		    col+1, pr[total++]);
        }
      }
    }
  }
}

/*============================================================*/
/* REST OF THE CODE BELOW WAS REALLY FOR main_old()           */
/* IN CASE I WANT TO REMOVE THE MEX INTERFACE                 */
/*============================================================*/
void makeRandIm(float *target, int nX, int nY)
{
  int ixy;

  for (ixy=0;ixy < nX*nY; ixy++) {
    target[ixy] = RanUniform();
  }
}

void printRandIm(float *target,int nX, int nY, float thresh)
{
  int ix, iy, ixy;
  /*
  for (ixy=0;ixy < nX*nY; ixy++) {
    fprintf(stderr,"%2.3f ",target[ixy]);
  }
  fprintf(stderr,"\n");
  */

  for (iy=0; iy<nY; iy++) {
    for (ix=0; ix<nX; ix++) {
      float tmp;

      ixy = index(ix,iy,nX);
      tmp = target[ixy];
      if (tmp < thresh)
	fprintf(stderr,"      ");
      else
	fprintf(stderr,"%2.3f ",tmp);
    }
    fprintf(stdout,"\n");
  }
  fprintf(stdout,"\n");
}

void setLinksRandIm(unsigned char *links[NUMLINKS], int nX, int nY )
{
  int ixy, i, ix, iy, kDir, ixs, iys;

  for (ixy=0; ixy<nX*nY; ixy++)
    for (i=0; i<NUMLINKS; i++)
      links[i][ixy] = 0;
  
  for (iy=0; iy<nY; iy++)
    for (ix=0; ix<nX; ix++) 
      {
	ixy = index(ix,iy,nX);
	for(kDir=0; kDir<NUMLINKS; kDir++) {
	  ixs = ix + dx[kDir];
	  iys = iy + dy[kDir];
	  if ((ixs>=0) && (ixs < nX) && (iys>=0) && (iys<nY)) {
	    links[kDir][ixy] = 1;
	  }
	}
      }
}

void printLabels(ccLabelStruct **label,int nX, int nY)
{
  int ix, iy, ixy;

  for (iy=0; iy<nY; iy++) {
    float tmp = 0;
    for (ix=0; ix<nX; ix++) 
      {
	ixy = index(ix,iy,nX);
	if (label[ixy] != NULL)
	  fprintf(stderr,"%4d",label[ixy]->count);
	else
	  fprintf(stderr,"    ");	  
      }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");
}

void processRandom()
{
  float *target;
  int nX = 10;
  int nY = 10;
  ccLabelStruct **label;
  unsigned char *links[NUMLINKS];
  int numLinks, i, kDir;
  int ix,iy, ixy, ixs, iys, ixys;
  float ampThresh;
  int numComp;
  float minTarget = 0.5;

  /* generate a random image as target */
  grabFloatMemory((float **)&target, nX*nY, "target");
  RanTimeSeed;

  /* generate a markov matrix */
  makeRandIm(target,nX,nY);
  printRandIm(target,nX,nY,0.5);

  /* set up empty labels */
  grabByteMemory((char **) &label,sizeof(ccLabelStruct *) * nX*nY,"label");
  for (ixy=0;ixy < nX*nY; ixy++) {
    label[ixy] = NULL;
  }

  /* Mark neighbouring pairs of edges for cocircularity */
  for (i=0; i<NUMLINKS; i++) 
    grabByteMemory((char **) links+i, nX*nY, "coCircularLinks");
  setLinksRandIm(links,nX,nY);

  fprintf(stderr, "Computing connected components....");
  connectLinks(label, target, minTarget, links,  nX, nY);
  fprintf(stderr, "done\n");

  /*
  fprintf(stderr, "Pruning small or weak connected components....");
  pruneLabelComp(label, nX, nY, 1, 0);
  fprintf(stderr, "done\n");
  */

  /* component count */
  numComp = countLabelComp(label, nX, nY);
  fprintf(stderr, "Found %d components.\n\n", numComp);

  printLabels(label,nX,nY);
  /* free stuff */
  free(target);
  free(label);
  for(i=0; i<NUMLINKS; i++) {
    free(links[i]);
  }

}

/*============================================================*/
float makeMarkovIm(float *target, int nX, int nY)
{
  int ixy, ix, iy;
  float tmp;
  int d;
  float minTarget;

  for (ixy=0;ixy < nX*nY; ixy++) {
    target[ixy] = RanUniform();
  }

  for (ix = 0; ix < nX; ix++)
    {
      d = index(ix,ix,nX);
      tmp = target[d];
      target[d] = 0.8 > tmp ? 0.8 : tmp;
    }

  for (iy=0; iy<nY; iy++) {
    tmp = 0;
    for (ix=0; ix<nX; ix++) {
      ixy = index(ix,iy,nX);
      tmp += target[ixy];
    }
    for (ix=0; ix<nX; ix++) {
      ixy = index(ix,iy,nX);
      target[ixy] /= tmp;
    }
  }

  minTarget = 100.0;
  for (ix=0;ix < nX; ix++) {
    d = index(ix,ix,nX);
    minTarget = minTarget < target[d] ? minTarget : target[d];
  }
  /*fprintf(stderr,"minTarget: %2.3f \n",minTarget);*/
  return minTarget;

}

void printMarkovIm(float *target,int nX, int nY, float thresh)
{
  int ix, iy, ixy;

  fprintf(stdout,"M = [ \n");
  for (iy=0; iy<nY; iy++) {
    for (ix=0; ix<nX; ix++) {
      float tmp;
      ixy = index(ix,iy,nX);
      tmp = target[ixy];
      if (tmp < thresh)
	fprintf(stderr,"      ");
      else
	fprintf(stderr,"%2.3f ",tmp);
      /*fprintf(stdout,"%2.3f ",tmp);*/
    }
    fprintf(stdout,"\n");
  }
  fprintf(stdout,"] ;\n");
  fflush(stdout);
}

void printMarkovLabels(ccLabelStruct **label,int nX)
{
  int ix ;

  /*fprintf(stdout,"lbl = [ \n");*/
  for (ix = 0; ix < nX; ix++)
    {
      if (label[ix] != NULL)
	{
	  /*fprintf(stdout,"%3d \n",label[ix]->count);*/
	  mexPrintf("%3d \n",label[ix]->count);
	}
      else
	{
	  /*fprintf(stdout,"-1 \n");*/
	  mexPrintf("-1 \n");
	}
    }
  /*fprintf(stdout,"] ;\n");*/
  fflush(stdout);
}

float makeSymmetricIm(float *target, int nX, int nY)
{
  int ixy1, ixy2, ix, iy;
  float minTarget, tmp;
  int d;

  for (ixy1 = 0; ixy1 < nX*nY; ixy1++)
    target[ixy1] = 0;

  /* set up all the off-diagonal elements */
  for (iy=0; iy<nY-1; iy++) {
    for (ix=iy+1; ix<nX; ix++) {
      ixy1 = index(ix,iy,nX);
      target[ixy1] = RanUniform();
      ixy2 = index(iy,ix,nX);
      target[ixy2] = target[ixy1];
    }
  }

  /* load the diagonal */
  for (ix = 0; ix < nX; ix++)
    {
      d = index(ix,ix,nX);
      tmp = RanUniform();
      /*target[d] = tmp;*/
      target[index(ix,ix,nX)] = 0.51 > tmp ? 0.51 : tmp;
      fprintf(stderr,"(%d %2.3f) ",d,target[d]);
    }
  fprintf(stderr,"\n");

  minTarget = 100.0;
  for (ix=0;ix < nX; ix++) {
    d = index(ix,ix,nX);
    minTarget = minTarget < target[d] ? minTarget : target[d];
  }
  /*fprintf(stderr,"minTarget: %2.3f \n",minTarget);*/
  return minTarget;

}



void processLines()
{
  /*
  int nX = 354;
  int nY = 354;
  */
  int nX = 342;
  int nY = 342;
  ccLabelStruct **label;
  int ix,iy, ixy, ixs, iys, ixys;
  float ampThresh;
  int numComp, t;
  float minTarget; /* 0.001 = boxes02-add.pgm */
  float *target;
  FILE *fp;
  char fnInput[MAXLEN];

  RanTimeSeed;

  fprintf(stderr, "filename for lines data > ");
  scanf("%s", fnInput);

  /* generate a random image as target */
  /*
    grabFloatMemory((float **)&target, nX*nY, "target");
    if ((fp = fopen("lines.LL","rb")) == NULL)
    {
    fprintf(stderr,"Error: Cannot Open LinesFile %s\n");
    exit(2);
    } 
    t = fread(target,sizeof(float),nX*nY,fp);
    if (t != nX*nY)
    {
    fprintf(stderr,"Incorrect data length:%d \n",t);
    exit(2);
    }
  */

  if ((fp = fopen(fnInput,"rb")) == NULL)
    {
      fprintf(stderr,"Error: Cannot Open Lines File %s\n",fnInput);
      exit(2);
    }
  t = read_pfm_image(&target, fp, fnInput, &nX, &nY);
  if (t < 0)
    {
      exit(2);
    }

  fprintf(stderr, "affty threshold  > ");
  scanf("%f", &minTarget);

  /*
  printMarkovIm(target,nX,nY,0);
  printMarkovIm(target,nX,nY,minTarget);
  */

  /* set up empty labels */
  grabByteMemory((char **) &label,sizeof(ccLabelStruct *) * nX,"label");
  for (ixy=0;ixy < nX; ixy++)
    label[ixy] = NULL;

  fprintf(stderr, "Computing connected components....");
  connectSymmetric(label, target, minTarget, nX, nY);
  fprintf(stderr, "done\n");

  fprintf(stderr, "Pruning small or weak connected components....");
  pruneLabelComp(label, nX, 1, 1, 0);
  fprintf(stderr, "done\n");

  /* component count */
  numComp = countLabelComp(label, nX, 1);
  fprintf(stderr, "Found %d components.\n\n", numComp);

  if (numComp > 0)
    printMarkovLabels(label,nX);

  /* free stuff */
  free(target);
  free(label);
}

/*============================================================*/

int main_old(int argc, char *argv[])
{
  int MARKOV = 0;
  int LINES = 0;
  int RANDOM = 0;
  char method[MAXLEN];

  fprintf(stderr, "which method > ");
  scanf("%s", method);

  if (strncmp(method,"mar",3) == 0)
    MARKOV = 1;
  if (strncmp(method,"lin",3) == 0)
    LINES = 1;
  if (strncmp(method,"ran",3) == 0)
    RANDOM = 1;
  
  if (MARKOV)
    processMarkov(); /* symmetric/markov matrix */
  if (LINES)
    processLines(); /* information from lines */
  if (RANDOM)
    processRandom(); /* 4 neigbors */
}

/*============================================================*/
void processMarkov()
{
  int nX = 5;
  int nY = 5;
  ccLabelStruct **label;
  unsigned char *links[NUMLINKS];
  int numLinks, i, kDir;
  int ix,iy, ixy, ixs, iys, ixys;
  float ampThresh;
  int numComp;
  float minTarget ;
  int MARKOV = 0;
  float *target;
  float target1[] = {0.510, 0.258, 0.211, 0.356, 0.391,
		    0.258, 0.818, 0.704, 0.858, 0.117, 
		    0.211, 0.704, 0.510, 0.460, 0.674, 
		    0.356, 0.858, 0.460, 0.510, 0.060, 
		    0.391, 0.117, 0.674, 0.060, 0.510}; /* 2 components */

  RanTimeSeed;

  /* generate a random image as target */
  grabFloatMemory((float **)&target, nX*nY, "target");
  minTarget = 0.5;

  /* generate a markov matrix */
  /*minTarget = makeMarkovIm(target,nX,nY);*/

  minTarget = makeSymmetricIm(target,nX,nY); 
  fprintf(stderr,"minTarget: %2.3f \n",minTarget);

  printMarkovIm(target,nX,nY,0);
  printMarkovIm(target,nX,nY,minTarget);

  /* set up empty labels */
  grabByteMemory((char **) &label,sizeof(ccLabelStruct *) * nX,"label");
  for (ixy=0;ixy < nX; ixy++)
    label[ixy] = NULL;

  fprintf(stderr, "Computing connected components....");
  /*connectMarkov(label, target, minTarget, nX, nY);*/
  connectSymmetric(label, target, minTarget, nX, nY);
  fprintf(stderr, "done\n");

  fprintf(stderr, "Pruning small or weak connected components....");
  pruneLabelComp(label, nX, 1, 1, 0);
  fprintf(stderr, "done\n");

  /* component count */
  numComp = countLabelComp(label, nX, 1);
  fprintf(stderr, "Found %d components.\n\n", numComp);

  if (numComp > 0)
    printMarkovLabels(label,nX);

  /* free stuff */
  /*free(target);*/
  free(label);
}

/*============================================================*/
/* following works for markov asymmetric matrix */
void connectMarkov(ccLabelStruct *label[], float *target, float minTarget,
		   int nX, int nY)
{
  int ix, iy, ixy, ixn, iyn, ixyn, kDir;
  ccLabelStruct *currentRoot, *nhbrLabel, *nhbrRoot;

  for (iy = 0; iy < nY; iy++) {
    ixy = index(iy,iy,nX);
    currentRoot = label[iy] = NULL;
    if (target[ixy] >= minTarget) {
      if (iy > 0) {
	/* loop over neighbors */
	for (ix = 0; ix < nX; ix++) {

	  if (ix != iy) {
	    /*ixyn = index(ix,iy,nX);*/
	    ixyn = index(iy,ix,nX);
	    nhbrLabel = label[ix];
	    fprintf(stderr,"\n (%d %d %2.3f) \n",ix,iy,target[ixyn]);
	    if ((nhbrLabel != NULL) && (target[ixyn] >= minTarget)) {
	      /* Connect this neighbour to current pixel */
	      nhbrRoot = findLabelRoot(nhbrLabel);
	      
	      if (currentRoot == NULL) { /* Current pixel unlabelled */
		label[iy] = nhbrRoot;  /* Use label of connected neighbour */
		currentRoot = nhbrRoot;
		currentRoot->size += 1;
		currentRoot->weight = max(currentRoot->weight, target[ixy]);
	      } else if (nhbrRoot != currentRoot)
		/* Merge label trees */
		currentRoot = mergeLabelRoots(nhbrRoot, currentRoot);
	    }
	  } /* ignore the diagonal: ix == iy */

	} /* End: loop over neighbours */
      }
      if (currentRoot == NULL) { /* No labelled neighbours found */
	/* Generate new label for this pixel */
	fprintf(stderr,"\n no neighbors %d \n",iy);
	label[iy] = newLabel();
	label[iy]->size = 1;
	label[iy]->weight = target[ixy];
      }
    }/* target[ixy] >= minTarget */
  } /* End: go over the diagonal */
  /* Completed labelling. */

  /* Reset all the labels to refer to only the root node */
  flattenLabels(label, nX, 1);
 
} /* End: connectMarkov */


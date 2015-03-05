#include "macros.h"
#include "utils.h"
#include "imageFile-io.h"
/*
#include "phase.h"
*/
#include "connCompLinks2.h"
#include "randNumRec.h"

/* Link directions in E, SW, S, SE order */
static int dx[NUMLINKS] = {1, -1, 0, 1};
static int dy[NUMLINKS] = {0, 1, 1, 1};

/* just so i remember */
/* # define  index(IX, IY, NX) ((IY)*(NX)+(IX)) */

extern void connectMarkov(ccLabelStruct *label[], float *target, 
			       float minTarget,int nX, int nY) ;
extern void connectSymmetric(ccLabelStruct *label[], float *target, 
			     float minTarget,int nX, int nY) ;

extern void processMarkov();
extern void processLines();

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
	fprintf(stdout,"%3d \n",label[ix]->count);
      else
	fprintf(stdout,"-1 \n");
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

  for (iy=0; iy<nY-1; iy++) {
    for (ix=iy+1; ix<nX; ix++) {
      ixy1 = index(ix,iy,nX);
      target[ixy1] = RanUniform();
      ixy2 = index(iy,ix,nX);
      target[ixy2] = target[ixy1];
    }
  }

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

int main(int argc, char *argv[])
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
/* following works for symmetric matrix */
void connectSymmetric(ccLabelStruct *label[], float *target, float minTarget,
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



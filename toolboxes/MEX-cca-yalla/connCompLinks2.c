/***********************************************************
File: connCompLinks2.c   Version 0.1

Performs connected components using links set using
the standard 8-neigbourhood of image pixels.  See diagram below.

Algorithm: Tree implementation of Merge-Find Sets using
path compression (see Aho, Hopcroft, Ullman, Data Structures
and Algorithms, p. 185.)

The labels are organized in a forest, with each tree corresponding
to one connected component.  The root of each tree provides 
the label for the component.  Merging components is done by merging trees.
Path compression is used whenever a tree is ascended to determine
the root.

ADJ  Sept, 2000
***********************************************************/
#include "macros.h"
#include "utils.h"
#include "connCompLinks2.h"

/* Here the links are from the current pixel to 
   the neighbour pixel, denoted by * :

            current -E- * 
            /  |  \
           /   |   \
         SW    S    SE
        /      |      \
       /       |       \
      *        *        *  
   Where:
     pixel current = (ix, iy),
     pixel nghbr * = (ix+dx[k], iy+dy[k]), for k = 0, ..., 3,
   with the nghbr's listed in the  E, SW, S, SE order.
******************************************************/
static int dx[NUMLINKS] = {1, -1, 0, 1};
static int dy[NUMLINKS] = {0, 1, 1, 1};

/** Return the shifts from the current pixel (ix, iy), to the
    neighbour pixel (ix+cdx[k], iy+cdy[k]) for each of the
    directions k = 0,..., NUMLINKS-1.
    Precondition: cdx and cdy must be arrays of length >= NUMLINKS.
**/
void getNghbrShifts(int cdx[], int cdy[])
{
  int i;
  for(i=0; i<NUMLINKS; i++) {
    cdx[i] = dx[i];
    cdy[i] = dy[i];
  }
} /* End: getNghbrShifts */

/** Generate a new empty label struct **/
ccLabelStruct * newLabel() 
{
  ccLabelStruct * n;
  grabByteMemory((char **)&n, sizeof(ccLabelStruct), "ccLabel");
  n->parent = NULL;
  n->size = 0;
  n->weight = 0.0;
  n->count = 0;
  return(n);
}

/** Merge trees to represent merging two components. 
    Return the root of the new tree.  **/
ccLabelStruct * mergeLabelRoots(ccLabelStruct * root1, ccLabelStruct * root2)
{
  /* Merge the smaller tree into the larger tree */
  if (root1->size > root2->size) {
    mergeLabelHelper(root2, root1);
    return(root1);
  }
  mergeLabelHelper(root1, root2);
  return(root2);
} /* End: mergeLabelRoots */

/** Merge smaller tree into larger tree **/
void mergeLabelHelper(ccLabelStruct * smallRoot, ccLabelStruct * largeRoot)
{
  smallRoot->parent = largeRoot;
  largeRoot->size += smallRoot->size;
  largeRoot->weight = max(largeRoot->weight, smallRoot->weight); 
} /* End: mergeLabelHelper */
  
/** Determine the root label for a node.  
    Precondition: 
      1) label != NULL
      2) for every node in the tree, current->size >
         sum of the sizes of all of the children of the current node.
         That is, every node must be directly referred to by at least
         one pixel.
    Postcondition:
      1) Performs path compression, so the label tree will be modified,
         if necessary, to ensure either label == root or 
         label->parent == root 
**/
ccLabelStruct * findLabelRoot(ccLabelStruct * label)
{
  ccLabelStruct *current, *next, *root;

  /* Ascend the tree to find the root */
  current = label;
  while((next = current->parent) != NULL)
   current = next;
  root = current;
   
  /* Do path compression.  That is, make all non-root members
     along the search path ascended in the loop above
     point directly to the root instead of pointing at some
     intermediate ancestor.  When relinking these nodes, modify 
     the sizes of the parents of relinked nodes to reflect the
     change in the number of pixels now contributed  to
     the root component.
  */
  if (label != root) {
    current = label;
    next = current->parent;  /* Cannot be null, since current != root. */
    while( next != root) {
      /* Link the current node to the root instead of it's parent */
      current->parent = root;
      
      next->size -= current->size; /* these are now pointing at the root,
                                      not at the next node.  */
      /* Note next->size remains > 0 by precondition 2 above. */

      /* Ascend the tree one step */
      current = next;
      next = current->parent;
    }
  }
  
  return(root);
} /* End: findLabelRoot */

/* 
Loop through array target, labelling any pixel with a 
target value > minTarget. 
Precondition:
 1) The arrays links[k] must be set > 0 for connected neighbours.
    The neighbourhood geometry (dx[k], dy[k]) is specified according to
    the diagram at the beginning of this file.
Postcondition:
 1) The returned labels have been flattened.  That is all
    merge trees consist only of the root nodes.  Therefore, either
        label[ixy] == NULL (unlabelled pixel)
    or
        label[ixy]->parent == NULL (root label).
 2) The size of each component, in terms of the number of pixels,
    and the maximum target weight occurring anywhere on a particular
    component, are stored in the fields
      label[ixy]->size, and  
      label[ixy]->weight,
     respectively.
 3) The labels are NOT counted (this can be done after pruning).
    The counts are cleared, that is,
      label[ixy]->count = 0
    for all labelled pixels.  See countLabelComp() below.
*/
void connectLinks(ccLabelStruct *label[], float *target, float minTarget,
		  unsigned char *links[NUMLINKS], int nx, int ny )
{
  int ix, iy, ixy, ixn, iyn, ixyn, kDir;
  ccLabelStruct *currentRoot, *nhbrLabel, *nhbrRoot;

  /* 
   To use links between neighbours in the directions (dx,dy) specified above,
   it is most convenient to work backwards through the image, starting
   in the bottom right corner. That way, the current node
   gets linked to previously labelled nodes (i.e. the *'d nodes).
  */
  for(iy=ny-1; iy>=0; iy--)
    for(ix=nx-1; ix>=0; ix--) { /* loop backwards over image */
      ixy = index(ix,iy, nx);
      currentRoot = label[ixy] = NULL;
      if (target[ixy] >= minTarget) {
        for(kDir = 0; kDir < NUMLINKS; kDir++) {  /* loop over neighbours */
          /* Neighbour index */
          ixn = ix + dx[kDir];
          iyn = iy + dy[kDir];
	  if ((ixn>=0) && (ixn < nx) && (iyn>=0) && (iyn<ny)) {
	    ixyn = index(ixn, iyn, nx);
            nhbrLabel = label[ixyn];
	    if (nhbrLabel != NULL && links[kDir][ixy]) {
              /* Connect this neighbour to current pixel */
              nhbrRoot = findLabelRoot(nhbrLabel);
              if (currentRoot == NULL) { /* Current pixel unlabelled */
                label[ixy] = nhbrRoot;  /* Use label of connected neighbour */
                currentRoot = nhbrRoot;
                currentRoot->size += 1;
                currentRoot->weight = max(currentRoot->weight, target[ixy]);
              } else if (nhbrRoot != currentRoot)
                /* Merge label trees */
                currentRoot = mergeLabelRoots(nhbrRoot, currentRoot);
            } /* connect the two components */
          } /* neighbour within bounds */
        } /* End: loop over neighbours */
        if (currentRoot == NULL) { /* No labelled neighbours found */
          /* Generate new label for this pixel */
          label[ixy] = newLabel();
          label[ixy]->size = 1;
          label[ixy]->weight = target[ixy];
        }
      } /* target pixel */
    } /* End: loop over image */
  /* Completed labelling. */

  /* Reset all the labels to refer to only the root node */
  flattenLabels(label, nx, ny);
 
} /* End: connectLinks */

/** Reset all labels to refer to the corresponding root node. 
    This flattens all label trees to consist of only the root node.
    Labels which become empty (unused) in this process are freed. 
**/
void flattenLabels(ccLabelStruct *label[], int nx, int ny)
{
  int ixy;
  ccLabelStruct *current, *root;
  for(ixy=0; ixy<nx*ny; ixy++) { /* loop over image */
    current = label[ixy];
    if (current != NULL && current->parent != NULL) {  
      /* Non-root label found */
      root = findLabelRoot(current);  /* This does path-compression,
                                         so current->parent now equals root */
      label[ixy] = root;  /* Switch label to refer to root node. */
      /* Maintain size of component labelled by old label by
         decrementing by one, to account for the label switch above.
         Delete this old label if it becomes unused. */
      if (current->size == 1)
	{
	  /*free(current);*/
	  utilFree((void **)&current);
	}
      else
        current->size -= 1;
    } /* Found non-root label */
  } /* End: loop over image */
} /* End: flattenLabels */

/** Free all malloc'd labels. 
    Precondition:
    1) the label trees have been flattened (i.e. either
       label[ixy] == NULL or label[ixy]->parent == NULL).
    2) the size field of the labels has been properly maintained to
       provide the number of pixels in each component. **/
void emptyLabels(ccLabelStruct *label[], int nx, int ny)
{
  int ixy;
  ccLabelStruct *current;

  for(ixy=0; ixy<nx*ny; ixy++) /* loop over image */
    if (label[ixy] != NULL) {
      current = label[ixy];
      label[ixy] = NULL;
      /* We are assuming the labels have been flattened, so
         current->parent == NULL, see precondition. */
      if (current->size == 1) {
        /*free(current);*/
	utilFree((void **)&current);
      } else {
        current->size -= 1;
      }
    } /* Non-null label */
} /* End: emptyLabels */

/** Prune components which are either small (i.e size < minSize)
    or weak (i.e. weight < minWeight). 
    Precondition:
    1) the label trees have been flattened (i.e. either
       label[ixy] == NULL or label[ixy]->parent == NULL).
    2) the size field of the labels has been properly maintained to
       provide the number of pixels in each component. **/
void pruneLabelComp(ccLabelStruct *label[], int nx, int ny, 
               int pruneSize, float minWeight)
{
  int nPrune, ixy;
  ccLabelStruct *current;

  nPrune=0;
  for(ixy=0; ixy<nx*ny; ixy++) {
    current = label[ixy];
    if ((current != NULL) && 
        ((current->size < pruneSize) || (current->weight < minWeight))) {
      label[ixy] = NULL;
      if (current->size == 1) {
        nPrune++;
        /*free(current);*/
	utilFree((void **)&current);
      } else
        current->size -= 1;
    }
  }
  /*fprintf(stderr," pruned %d comp...",nPrune);*/
}

/** Count distinct components. 
    Precondition:
    1) the label trees have been flattened (i.e. either
       label[ixy] == NULL or label[ixy]->parent == NULL,
       see flattenLabels()).
    2) the label counts have all been initialized to 0 (see
       newLabel() and clearLabelCounts()).
**/
long countLabelComp(ccLabelStruct *label[], int nx, int ny)
{
  int ixy;
  ccLabelStruct *current;
  long count = 0;
  /* Count all the distinct labels. */
  for(ixy=0; ixy<nx*ny; ixy++)
    if ((current = label[ixy]) != NULL) {
      if (current->count == 0) { /* New component found */
        current->count = ++count;
        /** Debug*/
  	/*fprintf(stderr,"ixy: %d, label: s %d, w %f, n %d\n",
	  current->size, current->weight, current->count);*/
      }
      /*
      else {
	fprintf(stderr,"ixy: %d, label: s %d, w %f, n %d\n",
		ixy, current->size, current->weight, current->count);
      }
      */
    }
  return(count);
} /* End: countLabelComp */

/** Clear counts stored in component labels. **/
void clearLabelCounts(ccLabelStruct *label[], int nx, int ny)
{
  int ixy;
  ccLabelStruct *current;
  for(ixy=0; ixy<nx*ny; ixy++)
    if ((current = label[ixy]) != NULL) 
      current->count = 0; /* old code: current->count == 0; */
} /* End: clearLabelCounts */


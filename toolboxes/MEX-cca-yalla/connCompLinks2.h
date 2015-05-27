/***********************************************************
File: connCompLinks2.c
***********************************************************/
#define NUMLINKS 4  
/* Links from each pixel to 4 (E, SW, S, SE) neighbours.
   Here the links are from the current pixel to 
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

#ifndef _CONN_COMP_H
#define _CONN_COMP_H

/* Unique labels for each component.  
   The labels are organized in a forest, with each tree corresponding
   to one connected component.  The root of each tree provides 
   the label for the component.  Merging components is done by merging trees.
   Path compression is used whenever a tree is ascended to determine
   the root (see Aho, Hopcroft and Ullman, Data structures and algorithms.)
*/
typedef struct cclabelStruct {
  struct cclabelStruct * parent;  /* Parent label structure (null if the
                                     current node is a root) */
  long size;    /* Number of edgels in current connected component */
  float weight; /* Maximum edgel weight in component */
  long count;   /* Unique integer label for component, labels
                   are not necessarily consecutive. */
} ccLabelStruct;

/** Function prototypes **/

/** Return the shifts from the current pixel (ix, iy), to the
    neighbour pixel (ix+cdx[k], iy+cdy[k]) for each of the
    directions k = 0,..., NUMLINKS-1.
    Precondition: cdx and cdy must be arrays of length >= NUMLINKS.
**/
void getNghbrShifts(int cdx[], int cdy[]);

/** Generate a new empty label struct **/
ccLabelStruct * newLabel();

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
		  unsigned char *links[NUMLINKS], int nx, int ny );

/** Prune components which are either small (i.e size < minSize)
    or weak (i.e. weight < minWeight). 
    Precondition:
    1) the label trees have been flattened (i.e. either
       label[ixy] == NULL or label[ixy]->parent == NULL).
    2) the size field of the labels has been properly maintained to
       provide the number of pixels in each component. **/
void pruneLabelComp(ccLabelStruct *label[], int nx, int ny, 
               int pruneSize, float minWeight);

/** Count distinct components. 
    Precondition:
    1) the label trees have been flattened (i.e. either
       label[ixy] == NULL or label[ixy]->parent == NULL,
       see flattenLabels()).
    2) the label counts have been initialized to 0 (see
       newLabel() and clearLabelCounts()).
**/
long countLabelComp(ccLabelStruct *label[], int nx, int ny);

/** Clear counts stored in component labels. **/
void clearLabelCounts(ccLabelStruct *label[], int nx, int ny);

/** Merge trees to represent merging two components. 
    Return the root of the new tree.  **/
ccLabelStruct * mergeLabelRoots(ccLabelStruct * root1, ccLabelStruct * root2);

/** Merge smaller tree into larger tree **/
void mergeLabelHelper(ccLabelStruct * smallRoot, ccLabelStruct * largeRoot);

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
ccLabelStruct * findLabelRoot(ccLabelStruct * label);

/** Reset all labels to refer to the corresponding root node. 
    This flattens all label trees to consist of only the root node.
    Labels which become empty (unused) in this process are freed. 
**/
void flattenLabels(ccLabelStruct *label[], int nx, int ny);

/** Free all malloc'd labels. 
    Precondition:
    1) the label trees have been flattened (i.e. either
       label[ixy] == NULL or label[ixy]->parent == NULL).
    2) the size field of the labels has been properly maintained to
       provide the number of pixels in each component. **/
void emptyLabels(ccLabelStruct *label[], int nx, int ny);

#endif

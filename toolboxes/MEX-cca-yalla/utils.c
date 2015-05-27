#include "string.h"
#include "stdio.h"
#include "stdlib.h"

#include "mex.h"

#include "utils.h"


static int memDebug = FALSE;
static FILE * fpMemDebug = NULL; /* Output file to dump allocate/free info */

void error(const char *s1, const char *s2)
 {
   char buffer[256];
   sprintf( buffer, "Error: %s %s\n", s1, s2 );
   mexErrMsgTxt( buffer );
 }

 /** Returns 0 when fname is Null, NULL, null, etc. **/
int checkForNull(char *fname) 
 {
  if (strcasecmp(fname, "NULL") == 0)
   return((int) 0);
  else
   return((int) 1);
 }

void swapPointers(void **a, void **b) 
{
  void *tmp;
  tmp = *a;
  *a = *b;
  *b = tmp;
}
 
void grabFloatMemory(float **x, int n, char * blurb)
 {
   if ((*x = (float *) mxMalloc((size_t) n * sizeof(float)) ) == NULL)
     error(" Unable to allocate memory for ", blurb);
   if (memDebug)
     echoMemAllocate((void *) *x, blurb);
 }
 
void grabByteMemory(char **x, int n, char *blurb)
 {
   if ((*x = (char *) mxMalloc((size_t) n * sizeof(char)) ) == NULL)
     error(" Unable to allocate memory for ", blurb);
   if (memDebug)
     echoMemAllocate((void *) *x, blurb);
 }
 
void grabIntMemory(int **x, int n, char *blurb)
 {
   if ((*x = (int *) mxMalloc((size_t) n * sizeof(int)) ) == NULL)
     error(" Unable to allocate memory for ", blurb);
   if (memDebug)
     echoMemAllocate((void *) *x, blurb);
 }

void utilFree(void ** p) {
   if (*p == NULL) {
     // fprintf("Warning: freeing null pointer.");
     return;
   }
   if (memDebug)
     echoMemFree(*p);
   mxFree(*p);
   *p = NULL;
}

/*** Memory Allocation/Free Debug Aids ***/
void setMemDebug(int setTo, char *fname) 
{
  closeEchoMem();  // Close any previous memDebug output file
  if (setTo) {  // Turn memory debugging on, use file fname 
    openEchoMem(fname);
  }
  memDebug = setTo;
}
 
void openEchoMem(char *fname) 
{
  if (fpMemDebug != NULL) 
    closeEchoMem();
  fpMemDebug = fopen(fname, "w");
  if (fpMemDebug == NULL)
    error("echoAllocate: cannot open memory debug output file", fname);
}

void closeEchoMem() 
{
  if (fpMemDebug != NULL) {
    fclose(fpMemDebug);
    fpMemDebug == NULL;
  }
}

void echoMemAllocate(void * p, char * blurb) 
{
  if (fpMemDebug == NULL) 
    error("echoMemAllocate:", "memory debug file not open");
  fprintf(fpMemDebug, "alloc:\t%x\t%s\n", p, blurb);
}

void echoMemFree(void * p) 
{
  if (fpMemDebug == NULL) 
    error("echoMemAllocate:", "memory debug file not open");
  fprintf(fpMemDebug, "free:\t%x\n", p);
}


/* input a single line from a file (up to the character '\n' or
   maxSize-1 characters, whichever comes first).  
   getLine returns:
     1    if maxSize-1 is reached (line continues beyond buffer)
     0    if the carriage return is read (end of line read)
    -1    if eof is bumped into or could not read one character from file.
************************************************************************/
int getLine(FILE *fpIn, char *line, int maxSize)
 {
  int k, endLine;
  char *p;
  endLine = 0;
  for(p=line, k=0; k<maxSize-1; k++, p++) {
   if (!fread(p, (int) 1, 1, fpIn)) {
    *p = '\0';
    return(-1);
   }
   if (*p == '\n') {
    endLine = 1;
    *p = '\0';
    break;
   }
  }
  if (endLine == 0) {
   line[maxSize-1] = '\0';
   return(1);
  } else
   return(0);
 }

/* input a single line from a file (up to the character '\n' or
   maxSize-1 characters, whichever comes first).  
   Strip leading blanks and tabs.
   getStripLine returns:
     1    if maxSize-1 is reached (line continues beyond buffer)
     0    if the carriage return is read (end of line read)
    -1    if eof is bumped into or could not read one character from file.
************************************************************************/
int getStripLine(FILE *fpIn, char *line, int maxSize)
 {
  int k, endLine, foundNonWht;
  char *p;
  endLine = 0;
  foundNonWht = 0;
  
  for(p=line, k=0; k<maxSize-1; k+=foundNonWht, p+=foundNonWht) {
   if (!fread(p, (int) 1, 1, fpIn)) {
    *p = '\0';
    return(-1);
   }
   if (*p == '\n') {
    endLine = 1;
    *p = '\0';
    break;
   }
   if (!foundNonWht) {
    /* Check for first appearance of non blank or tab */
    if (!(*p == ' ' || *p == '\t'))
     foundNonWht = 1;
   }
  }
  if (endLine == 0) {
   line[maxSize-1] = '\0';
   return(1);
  } else
   return(0);
 }

/* Flush a line from an input file (up to the character '\n' or EOF)
   flushLine returns:
    0 if successful (character '\n' found)
   -1 if eof bumped into or cannot read character.
*************************************************************/
int flushLine(FILE *fpIn)
 {
  char p;
  do {
   if (!fread(&p, (int) 1, 1, fpIn))
    return(-1);
  } while (p != '\n');
  return(0);
 }

/* Read next integer in input file, ignoring comments
   which are indicated by # and continue to the next newline.
   Returns 1 (success), 0 (failure),  else EOF */
int readOneInt(FILE *fp, int *pn)
{
  char comment[MAXLEN];
  int numRead, more;
 
  do {
    numRead = fscanf(fp, "%d", pn);

    if (numRead != 1) { /* Strip expected comment */
      if ((more=getStripLine(fp, comment, MAXLEN)) == EOF) {
        return(EOF);
      }
      if (comment[0] == '#') { /* Comment beginning detected */
        /* Strip comment */
        numRead = 0;
        if (more == 1)  /* Comment longer than MAXLEN */
          if (flushLine(fp) == EOF)
            return(EOF);
      } else { /* Garbled input */
        return(0);
      }
    }
  } while( numRead != 1);

  return(numRead);
}

/* Read next float in input file, ignoring comments
   which are indicated by # and continue to the next newline 
   Returns 1 (success), 0 (failure),  else EOF */
int readOneFloat(FILE *fp, float *px)
{
  char comment[MAXLEN];
  int numRead, more;
 
  do {
    numRead = fscanf(fp, "%e", px);
  
    if (numRead != 1) { /* Strip expected comment */
      if ((more=getStripLine(fp, comment, MAXLEN)) == EOF) {
        return(EOF);
      }
      if (comment[0] == '#') { /* Comment beginning detected */
        numRead = 0;
        if (more == 1)  /* Comment longer than MAXLEN, flush it */
          if (flushLine(fp) == EOF)
            return(EOF);
      } else { /* Garbled input */
        return(0);
      }
    }
  } while( numRead != 1);

  return(numRead);
}


  /* Read next string in input file, ignoring comments
   which are indicated by # and continue to the next newline.
   It is assumed the string is long enough to contain the input. 
   Returns 1 (success), else EOF */
int readOneString(FILE *fp, char *s)
{
  char comment[MAXLEN];
  int numRead;
 
  do {
    numRead = fscanf(fp, "%s", s);
    if (numRead != 1) 
      return(EOF);
    if (s[0] == '#') { /* Strip comment */
      numRead = 0;
      if (flushLine(fp) == EOF)
        return(EOF);
    }
  } while( numRead != 1);
  return(numRead);
}

/**********************************************
Sort float array w into decreasing order, using insertion sort
Return permutation array iPerm *... ith sorted value is w[iPerm[i]] 
***********************************************/
void sortFloatArray( int *iPerm, float *w, int numSamp)
{
  float tmp;
  int i, j, iTmp;
   
  for(i=0; i<numSamp; i++)
    iPerm[i] = i;

  for(i=0; i<numSamp-1; i++) {
    /* find max in remaining set i..numSamp */
    tmp = w[iPerm[i]];
    iTmp = i;
    for(j=i+1; j<numSamp; j++)
      if (w[iPerm[j]] > tmp) {
	tmp = w[iPerm[j]];
	iTmp = j;
      }
    /* Switch */
    j = iPerm[i];
    iPerm[i] = iPerm[iTmp];
    iPerm[iTmp] = j;
  }
}


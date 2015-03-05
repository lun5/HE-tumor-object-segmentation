
# include "macros.h"

void error(const char *s1, const char *s2);

/** Returns 0 when fname is Null, NULL, null, etc. **/
int checkForNull(char *fname);
 
void swapPointers(void **a, void **b);

void grabFloatMemory(float **x, int n,char *blurb);

void grabByteMemory(char **x, int n, char *blurb);

void grabIntMemory(int **x, int n, char *blurb);
 
void utilFree(void ** p);

/*** Memory Allocation/Free Debug Aids ***/
void setMemDebug(int setTo, char *fname) ;
void openEchoMem(char *fname);
void closeEchoMem();
void echoMemAllocate(void * p, char * blurb);
void echoMemFree(void * p);

/* input a single line from a file (up to the character '\n' or
   maxSize-1 characters, whichever comes first).  
   getLine returns:
     1    if maxSize-1 is reached (line continues beyond buffer)
     0    if the carriage return is read (end of line read)
    -1    if eof is bumped into or could read one character from file.
************************************************************************/
int getLine(FILE *fpIn, char *line, int maxSize);

/* input a single line from a file (up to the character '\n' or
   maxSize-1 characters, whichever comes first).
   Strip leading blanks and tabs.
   getStripLine returns:
     1    if maxSize-1 is reached (line continues beyond buffer)
     0    if the carriage return is read (end of line read)
    -1    if eof is bumped into or could not read one character from file.
************************************************************************/
int getStripLine(FILE *fpIn, char *line, int maxSize);

/* Flush a line from a file (up to the character '\n' or EOF)
   flushLine returns:
    0 if successful (character '\n' found)
   -1 if eof bumped into or cannot read character.
*************************************************************/
int flushLine(FILE *fpIn);

/* Read next integer in input file, ignoring comments
   which are indicated by # and continue to the next newline. 
   Returns 1 (success), 0 (failure), else EOF */
int readOneInt(FILE *fp, int *pn);

/* Read next float in input file, ignoring comments
   which are indicated by # and continue to the next newline 
   Returns 1 (success), 0 (failure),  else EOF */
int readOneFloat(FILE *fp, float *px);

/* Read next string in input file, ignoring comments
   which are indicated by # and continue to the next newline.
   It is assumed the string is long enough to contain the input. 
   Returns 1 (success), else EOF */
int readOneString(FILE *fp, char *s);

/**********************************************
Sort float array w into decreasing order, using insertion sort
Return permutation array iPerm *... ith sorted value is w[iPerm[i]] 
***********************************************/
void sortFloatArray( int *iPerm, float *w, int numSamp);



/************************* Version 2.0.1 *********************
 Changes from version 2.0:
  Added range_byte_image
  Changed all scanf's to annotated input form, i.e. used readOne*'s
*****/
#include "utils.h"
#include "mex.h"

# include "endianness.h"
#ifdef _LITTLE_ENDIAN
#define BYTE_ORDER_TYPE 'P'   /* for PC style byte order */
#else
#ifdef _BIG_ENDIAN
# define BYTE_ORDER_TYPE 'U'   /* for byte order in Unix boxes */
#else
# define BYTE_ORDER_TYPE '?'  /* unknown byte order */
#endif
#endif

#include "imageFile-io.h"

# define   	MAXPAD		6 /* Max length of file sequence number */
# define	whtspc(A)	((A)==' ' || (A)=='\t' || (A)=='\n')

static char *fltFileTail=".flt";
static char *byteFileTail=".bin";
static char *headerFileTail=".hdr";
static char *pgmFileTail=".pgm";
static char *ppmFileTail=".ppm";
static char *pfmFileTail=".pfm";

static char *blank=" ";
static int minHeadSize=15; /* Number of char's in smallest pgm header, useful
                            when the header is to be stripped based on size. */

long strtol();

char getByteOrderType()
{
  return(BYTE_ORDER_TYPE);
}

int read_image_params(char *fnInput, int *pPgmInput,
  int *pNxImage, int *pNyImage, int *pNHead)
 /* Read image filename, image type, and necessary size information */
 {
  char dumStr[255];

  /*** Read image file type and filename ****/
  fprintf(stderr, "Input image file type (pgm or raw) and  filename\n    ");

  if ((readOneString(stdin, dumStr) == EOF) ||
      (readOneString(stdin, fnInput) == EOF))
   return(EOF);
  
  if (strcmp(dumStr, "pgm") == 0) {
   *pPgmInput = 1;
  } else if (strcmp(fnInput, "raw") == 0) {
    *pPgmInput = 0;
    fprintf(stderr, "image size: nHeader, nx, ny :");
    readOneInt(stdin, pNHead);
    readOneInt(stdin, pNxImage);
    readOneInt(stdin, pNyImage);
  } else {
    error("Unrecognized file type",fnInput);
  }
  return(0);
 }

int read_file_path_and_root(char *pathName, char *fileRoot)
 {
  int nullFlag;

  if ((readOneString(stdin, pathName) == EOF) ||
      (readOneString(stdin, fileRoot) == EOF)) 
   return(EOF);

  nullFlag = 0;
  /* Ignore filenames or path "NULL" */
  if (checkForNull(pathName) == 0) {
   pathName[0] = 0;
   nullFlag = NULL_NAME_FLAG;
  }
  if (checkForNull(fileRoot) == 0) {
   fileRoot[0] = 0;
   nullFlag = NULL_NAME_FLAG;
  }

  return(nullFlag);
 }
  
/******************
fopenInputFile:
   Attempt to open input image file with the name fnIn
   Returns fpInput if successful, NULL otherwise
************/
FILE *fopenInputFile(char *fnIn)
{
    FILE *infile;
    if ((infile = fopen(fnIn, "rb")) == NULL) {
      return((FILE *) NULL);
    }
    return(infile);
}

/******************
freopenInputFile:
   Attempt to reopen input image file with the name fnIn
   on an existing stream fp or, if fp==NULL, then it
   attempts to open a new stream. 
   Returns fpInput if successful, NULL otherwise
******************/
FILE * freopenInputFile(char* fnIn, FILE* fp)
{
    FILE *infile;
    if (fp == NULL) 
     infile = fopen(fnIn,"rb");
    else 
     infile = freopen(fnIn,"rb",fp);
    return(infile);
}

/******************
freopenOutputFile:
   Attempt to reopen output image file with the name fnIn
   on an existing stream fp or, if fp==NULL, then it
   attempts to open a new stream. 
   Returns fpOutput if successful, NULL otherwise
******************/
FILE * freopenOutputFile(char* fnOut, FILE* fp)
{
    FILE *outfile;
    if (fp == NULL) 
     outfile = fopen(fnOut,"wb");
    else 
     outfile = freopen(fnOut,"wb",fp);
    return(outfile);
}

/******************
fopenInputFileNum:
   Attempt to open input image file with the name of the form:
    filename = headIn0..0###TailIn
   Where the image number 0...0### is possibly padded with leading
   zeroes to have a total length of at most MAXPAD.
   
   Routine openInputFileNum returns FILE* infile if successful, NULL otherwise
************/
FILE* fopenInputFileNum(char* fnIn, char* headIn, 
                       int* ppad, int kImage, char* tailIn)
{
    char fileNum[MAXLEN];
    int numLen;
    FILE *infile;

    sprintf(fileNum, "%d", kImage);
    numLen = strlen( fileNum );

    *ppad = numLen;
    while((*ppad)<=MAXPAD) {
      build_filename(fnIn, headIn, *ppad, kImage, tailIn);
      /***
      fprintf(stderr,"DEBUG: Attempting to open %s\n", fnIn);
      ***/
      if ( (infile = fopenInputFile(fnIn)) == NULL ) {
       (*ppad)++;
      } else {
       break;
      }
    }
  
    if (*ppad <= MAXPAD) {
      /*
     fprintf(stderr," Opened input file %s\n",fnIn);
      */
     return(infile); /* Success */
    } else { /* Return unpadded name and error code of NULL */
     *ppad = 0;
     build_filename(fnIn, headIn, *ppad, kImage, tailIn);
/***
     fprintf(stderr," Unsuccesful opening input file %s\n", fnIn);
****/
     return(NULL);
    }
}

/******************
build_filename:
   Construct filename of the form:
    filename = head###tail
   Where there are at least `pad' digits in the file number, padded with
   leading zeroes if needed.

************/
void build_filename(char *fname, char *head, int pad, int numImage, char *tail)
{
    char fileNum[MAXLEN];
    int numLen, k;

    sprintf(fileNum, "%d", numImage);
    numLen = strlen( fileNum );

    /* copy header to filename */
    strcpy(fname, head);

    /* put leading zeroes onto filename */
    for(k=numLen; k<pad; k++) 
     strcat(fname, "0");

    /* put file number onto end of filename */
    strcat(fname, fileNum);

    /* cat tail onto end of filename */
    strcat(fname, tail);
  
}

/* return an integer as a string, padded with leading
   zeroes to make the string at least of length pad. */
void pad_integer(char *padInt, int n, int pad)
{
    char tailInt[MAXLEN];
    int numLen, k;

    sprintf(tailInt, "%d", n);
    numLen = strlen( tailInt );

    /* put leading zeroes onto filename */
    strcpy(padInt,"");
    for(k=numLen; k<pad; k++) 
     strcat(padInt, "0");

    /* put file number onto end of filename */
    strcat(padInt, tailInt);
}

/*************************************************************
next_image_name:
   Find the next image name in a sequence, adjust padding
   zeroes and skip missing file numbers.
   Return: (-1) on failure
           (0) success, with:
   On successful return,
               imageName is the desired filename
               *pFileNum is the desired sequence number
*************************************************************/
int next_image_name(char *imageName, char *imageFileRoot, char *imageFileTail,
 int *pFileNum, int maxNum) 
 {
  int pad;
  FILE* infile;

  /* Try to open image file */
  while(*pFileNum <= maxNum) {
   if ((infile = fopenInputFileNum(imageName, imageFileRoot,
                           &pad, *pFileNum, imageFileTail)) != NULL) {
     fclose(infile); /* Just want the name and number, close the file for now*/
     break; /* Opened file successfully */
   }
   (*pFileNum)++;
  }
  if (*pFileNum > maxNum) {
       return(-1);  /* Failure return code */
  }
  return(0);
 }

void read_byte_image_strip_bytes(unsigned char* image, FILE* fpIn, char* fnIn,
                            int nx, int ny, int bytes)
{
  char *header;
  int sizeInput;

  if (bytes > 0) {

   grabByteMemory((char **) &header, bytes, "strip-header");

   sizeInput = bytes * sizeof( char );
   if (!fread((char *) header, sizeInput,1, fpIn))
    error("unexpected end of file on", fnIn);

   utilFree((void **) &header);
  }

  sizeInput = nx * ny * sizeof( char );
  if (!fread((char *) image, sizeInput, 1, fpIn)) {
    fprintf(stderr,"error: unexpected end of file on %s\n", fnIn);
    exit(1);
  }
 
}

void read_byte_image(unsigned char* image, FILE* fpIn, char* fnIn, 
                     int nx, int ny)
{
  int sizeInput;

  sizeInput = nx * ny * sizeof( char );
  if ( !fread((char *) image, sizeInput,1,fpIn)) {
    fprintf(stderr,"error: unexpected end of file on %s\n", fnIn);
    exit(1);
  }
 
}

void read_float_image_strip_bytes(float* imageFlt, FILE* fpIn, char* fnIn,
				  int nx, int ny, int bytes)
{
  char *header;
  int sizeInput;

  if (bytes > 0) {

   grabByteMemory((char **) &header, bytes, "strip-header");

   sizeInput = bytes * sizeof( char );
   if ( !fread((char *) header, sizeInput, 1, fpIn)) 
    error("unexpected end of file on", fnIn);

   utilFree((void **) &header);
  }

  sizeInput = nx * ny * sizeof( float );
  if ( !fread((char *) imageFlt, sizeInput, 1, fpIn)) {
    fprintf(stderr,"error: unexpected end of file on %s\n", fnIn);
    exit(1);
  }
 
}

void read_float_image(float* imageFlt, FILE* fpIn, char* fnIn, int nx, int ny)
{
  int sizeInput;

  sizeInput = nx * ny * sizeof( float );
  if ( !fread((char *) imageFlt, sizeInput,1,fpIn)) {
    fprintf(stderr,"error: unexpected end of file on %s\n", fnIn);
    exit(1);
  }
 
}

void read_int_image_strip_bytes(unsigned int * image, FILE * fpIn, char * fnIn,
				int nx, int ny, int bytes)
{
  char *header;
  int sizeInput;

  if (bytes > 0) {

   grabByteMemory((char **) &header, bytes, "strip-header");

   sizeInput = bytes * sizeof( char );
   if ( !fread((char *) header, sizeInput, 1, fpIn)) 
    error("unexpected end of file on", fnIn);

   utilFree((void **) &header);
  }

  sizeInput = nx * ny * sizeof( unsigned int );
  if ( !fread((char *) image, sizeInput, 1, fpIn)) {
    fprintf(stderr,"error: unexpected end of file on %s\n", fnIn);
    exit(1);
  }
 
}

void read_int_image(unsigned int *image, FILE *fpIn, char *fnIn, int nx, int ny)
{
  int sizeInput;

  sizeInput = nx * ny * sizeof( unsigned int );
  if ( !fread((char *) image, sizeInput,1,fpIn)) {
    fprintf(stderr,"error: unexpected end of file on %s\n", fnIn);
    exit(1);
  }
 
}

void int_to_float_image(float **pimage, unsigned int *intImage, int nx, int ny)
{
  int i;

  grabFloatMemory( pimage, nx*ny, "int_to_float_image");
  for(i=0; i<nx*ny; i++)
   (*pimage)[i] = (float) intImage[i];
 
}

/*
Read in image from path and store in image.
This function reads:
  1) images saved by save_image where image name is path
  2) straight byte images named path. */

void load_image(char *path, float **pImage, int *pnx, int *pny)
{
  unsigned char *temp;
  float *image;
  char dataIn[MAXLEN], headIn[MAXLEN];
  int ixy, i, j, sizeInput;
  int nx, ny;
  float min, max, scale;
  char *datatail=".bin", *headtail=".hdr";
  FILE *fpIn;

  /* Read header file */
  strcpy(headIn, "");
  strcpy(dataIn, "");
  strcat(strcat(headIn, path), headtail);
  fprintf(stderr, "headIn: %s\n", headIn);
  if ( (fpIn = fopenInputFile(headIn)) == NULL ) {
    error(" cannot open header input file ", headIn);
  }
  else {
    readOneInt(fpIn, &nx);
    readOneInt(fpIn, &ny);
    readOneFloat(fpIn, &min);
    readOneFloat(fpIn, &max);
    readOneFloat(fpIn, &scale);
    fprintf(stderr, "nx: %d\n", nx);
    fprintf(stderr, "ny: %d\n", ny);
  }
  fclose(fpIn);

  /* Read image file */
  strcat(strcat(dataIn, path), datatail);
  fprintf(stderr, "dataIn: %s,  path: %s  \n\n", dataIn, path);
  if ((fpIn = freopenInputFile(dataIn, fpIn)) == NULL) 
      error("Can't open input files, load_image %s", dataIn);
  grabByteMemory((char **)&temp, nx*ny, "loadImage: temp");
  grabFloatMemory(pImage, nx*ny, "loadImage: pImage");
  image = *pImage;

  sizeInput = nx*ny;
  if (!fread((char *) temp, sizeInput, 1, fpIn)) { 
    utilFree((void **) &temp);
    utilFree((void **) pImage);
    *pnx = *pny = 0;
    error("unexpected end of file on", dataIn);
  }
  *pnx = nx;
  *pny = ny;
  for(i=0;i<ny;i++){
    for(j=0;j<nx;j++){
      ixy = (i*nx) + j;
      image[ixy] = (((float) temp[ixy]) / scale) + (float) min;
    }
  }
  utilFree((void **) &temp);
  fclose(fpIn);
}

int read_pgm_image(unsigned char **pimage, FILE *fpIn, char *fnIn, 
		   int *pnx, int *pny)
{
  int sizeInput;

  /* Read image size off of pgm header, if possible */
  if (stripPgmHeader(fpIn, pnx, pny) < 0)
   return(-1);

  sizeInput = (*pnx) * (*pny) * sizeof( char );
  grabByteMemory((char **)pimage, sizeInput, fnIn);

  if (!fread((char *) (*pimage), sizeInput, 1, fpIn)) { 
    utilFree((void **) pimage);
    error("unexpected end of file on", fnIn);
  }

  return(0);
 
}

int read_ppm_image(unsigned char *pimage[3], FILE *fpIn, char *fnIn, 
		   int *pnx, int *pny)
{
  unsigned char *raw;
  int i0, ixy, sizeInput;

  /* Read image size off of pgm header, if possible */
  if (stripPpmHeader(fpIn, pnx, pny) < 0)
   return(-1);

  sizeInput = (*pnx) * (*pny) * sizeof( char );
  grabByteMemory((char **)pimage, sizeInput, fnIn);
  grabByteMemory((char **)pimage+1, sizeInput, fnIn);
  grabByteMemory((char **)pimage+2, sizeInput, fnIn);

  sizeInput = 3 * (*pnx) * (*pny) * sizeof( char );
  grabByteMemory((char **)&raw, sizeInput, fnIn);
  if ( !fread((char *) raw, sizeInput,1,fpIn)) 
    error("unexpected end of file on", fnIn);

  for(ixy=0; ixy< (*pnx) * (*pny); ixy++) {
   i0 = 3 * ixy;
   pimage[0][ixy] = raw[i0];
   pimage[1][ixy] = raw[i0+1];
   pimage[2][ixy] = raw[i0+2];
  }
  utilFree((void **) &raw);
  return(0);
 
}

int read_pfm_image(float **pimage, FILE *fpIn, char *fnIn, 
		   int *pnx, int *pny)
{
  int sizeInput;
  char byteOrderType;

  /* Read image size off of pfm header, if possible */
  if (stripPfmHeader(fpIn, pnx, pny, &byteOrderType) < 0) {
   fprintf(stderr, "cannot read header of file %s\n", fnIn);
   return(-1);
  }

  if (byteOrderType == '?') 
    error("unrecognized byteOrderType type in file", fnIn);

  sizeInput = (*pnx) * (*pny);
  grabFloatMemory(pimage, sizeInput, fnIn);

  if (!fread((char *) (*pimage), sizeInput * sizeof(float), 1, fpIn))
    error("unexpected end of file on", fnIn);

  /* Swap byte order if file written on a different machine (PC/Unix box). */
  if (byteOrderType != BYTE_ORDER_TYPE) {
    fprintf(stderr, " Swapping bytes in float image from file %s\n", fnIn);
    floatImageSwapBytes(*pimage, *pnx, *pny);
  }

  return(0);
 
}

void readFilterOutputs(FILE *fpIn, char* fname, char header[][HEADER], 
		       int nHeader, int nByte, int nOrient, int nDim, 
		       int sizeImage, 
		       char **image, int *pnxImage, int *pnyImage)
 /****************************************************************
  Read image data from file "fname", assumes open with descriptor fd.
    Format:  Assumes data has nOrient different channels, and each
             channel is of dimension nDim  (eg. a complex filter response
             could have 4 orientations each of which is of dimension 2)
             Data has nByte 's per pixel.
    If nHead > 0 the procedure reads the size of the image from the
    header, else
    if nHead == 0, then the procedure assumes a file of size nByte.
    
    Order:
             global header  ( size nHeader )
             header orientation 1 ( size nHeader )
             orientation 1, dimension 1
             orientation 1, dimension 2
             ... 
             orientation 1, dimension nDim
             header orientation 2 ( size nHeader )
             orientation 2, dimension 1
             ...
    Returns:
             headers for each orientation (not the global header)
         and pointers to memory blocks (in above order)
  *********************************************************************/      
      
 {
    int nx0, ny0, iMap, i, j, checkParam;
    unsigned int sizeInput, sizeJunk;
    char name[HEADER], scr[MAXLEN];

    /* Read global header into local variable char * name */
    if (nHeader > 0) {
      /* Read first HEADER characters of header */
      sizeInput = min(HEADER, nHeader);
      if ( ! fread(name, sizeInput,1, fpIn) ) 
        mexErrMsgTxt("unexpected end of file on input file");
      sscanf(name,"%s %d %d %d %d", scr, &nx0, &ny0, &i, &j);
   
      *pnxImage = nx0;
      *pnyImage = ny0;

      /* Check header data with calling parameters */
      checkParam = 1;
      if ((nByte == 4) && (strcmp(scr,"float") != 0))
        checkParam = 0;
      if ((nByte == 1) && (strcmp(scr,"byte") != 0))
        checkParam = 0;
      if (i < nOrient || j != nDim)
        checkParam = 0;

      if (checkParam == 0) {
        mexPrintf("DEBUG: nByte %d scr %s i %d nOrient %d j %d nDim %d\n",
              nByte, scr, i, nOrient, j, nDim);
        mexErrMsgTxt("unexpected format of input file");
      }
  
      j = nHeader - sizeInput;
      while (j>0) {
       /* Strip remaining characters from header */
       sizeJunk = min(HEADER, j);
       if ( ! fread(name, (unsigned int) sizeJunk, 1, fpIn) ) 
        mexErrMsgTxt("unexpected end of file on input file");
       j -= sizeJunk;
      }
    } 

    iMap = 0;
    for(i=0; i<nOrient; i++) {
     if (nHeader > 0) {
      /* Read first HEADER characters of header */
      sizeInput = min(HEADER, nHeader);
      if ( ! fread( header[i], sizeInput, 1, fpIn) ) 
        mexErrMsgTxt("unexpected end of file on input file");

      j = nHeader - sizeInput;
      while (j>0) {
       /* Strip remaining characters from header */
       sizeJunk = min(HEADER, j);
       if ( ! fread(name, (unsigned int) sizeJunk, 1, fpIn) ) 
        mexErrMsgTxt("unexpected end of file on input file");
       j -= sizeJunk;
      }
     } 

     /* Set image size:  Use header info if parameter sizeImage <=0
                         Otherwise use parameter sizeImage */
     if (sizeImage <= 0 && nHeader > 0)
      sizeInput = nx0 * ny0 * nByte;
     else 
      sizeInput = sizeImage;
     
     /* Read each dimension */
     for(j=0; j< nDim; j++) {
      
        /* Get enough space for image */
        grabByteMemory((char**)&(image[iMap]), sizeInput, "readFilterOutputs");
        if( image[iMap] == NULL )
           mexErrMsgTxt(" Unable to allocate memory while reading");
   
        /* Read filter amplitude and phase data */
        if ( ! fread((char *)image[iMap], sizeInput, 1, fpIn) )
           mexErrMsgTxt("unexpected end of file on input image");
           
        iMap++;
        
     }
    }        
    fclose(fpIn);
 }

/* Strip a pgm form header from input file, return size of
   image in *pnx, *pny, and return the file type string
   in *fType.
   The recognized types are:
     - "P5" - pgm image
     - "P6" - ppm image
     - "FP" - pfm image, image of floats in a PC format
     - "FU","F5" - pfm image, image of floats in a UNIX machine format
   stripPnmHeader returns -1 if file type not recognized.
*********/
int stripPnmHeader(FILE *fpIn, int *pnx, int *pny, char *fType)
{
  return(stripPgmFormHeader(fpIn, pnx, pny, fType));
}

/* Strip a pgm form header from input file, return size of
   image in *pnx, *pny, and return the file type string
   in *fType.
   The recognized types are:
     - "P5" - pgm image
     - "P6" - ppm image
     - "FP" - pfm image, image of floats in a PC format
     - "FU","F5" - pfm image, image of floats in a UNIX machine format
   stripPgmFormHeader returns -1 if file type not recognized.
*********/
int stripPgmFormHeader(FILE *fpIn, int *pnx, int *pny, char *fType)
 {
   char line[MAXLINE_PGM], *cPtr;
   int i, paramCount;
   long param[3];

   /* Strip file type in first two characters */
   fread(line, (int) 2, 1, fpIn);
   line[2] = '\0';
   strcpy(fType, line);
   /* Check if file type is one of the recognized types */
   if (!((strcmp(line, "P5") == 0)
         || (strcmp(line, "P6") == 0)
         || (strcmp(line, "FP") == 0)
         || (strcmp(line, "F5") == 0)
         || (strcmp(line, "FU") == 0)))
     return(-1);   /* unrecognized file type */

   fread( line, (int) 1, 1, fpIn);
   if (!whtspc(line[0]))
        return(-1);

   /* Read parameters pnx, pny, ngray */
   paramCount = 0;
   while(paramCount < 3) {  /* Seek next parameter */

        /* Strip white space */
        while(whtspc(line[0])) {
         if (! fread(line, (int) 1, 1, fpIn)) 
          return(-1);
        }

        if (line[0] == '#') {

         /* Strip comment */
         while(line[0] != '\n') { 
          if (! fread(line, (int) 1, 1, fpIn)) 
           return(-1);
         } 

        } else {

         /* Read parameter */
         i = 0;
         while( !whtspc(line[i]) && i < MAXLINE_PGM) { 
          i++;
          if (!fread( &(line[i]), (int) 1, 1, fpIn)) 
           return(-1);
         }
         param[paramCount] = strtol(line, &cPtr, 10);
         if (cPtr == NULL) 
           return(-1);
         paramCount++;
         line[0] = line[i]; /* re-prime whitespace stream */
        }
   }

   *pnx = param[0];
   *pny = param[1];

   if (*pnx < 0 || *pnx > MAXSIZE_PGM || *pny < 0 || *pny > MAXSIZE_PGM )
    /* ASSUMES an error reading the header if pnx or pny out of range */
    return(-1); 
   else
    return(0);
 } 

int stripPgmHeader(FILE *fpIn, int *pnx, int *pny) 
{
  char fType[3];
  if (stripPgmFormHeader(fpIn, pnx, pny, fType) < 0)
    return(-1);
  if (strcmp(fType, "P5") != 0)
    return(-1);
  return(0);
}

int stripPpmHeader(FILE *fpIn, int *pnx, int *pny)
{
  char fType[3];
  if (stripPgmFormHeader(fpIn, pnx, pny, fType) < 0)
    return(-1);
  if (strcmp(fType, "P6") != 0)
    return(-1);
  return(0);
}

/* Strip a pgm form header from a float file (i.e. a pfm format
   file).
   Return the image size, and a single character indicating
   the machine type:
   byteOrderType = 'P' for a PC
                  = 'U' for a Unix box (Sun, SGI)
   Function returns 0 if successful, otherwise -1.
*********/
int stripPfmHeader(FILE *fpIn, int *pnx, int *pny, char* pByteOrderType)
{
  char fType[3];
  *pByteOrderType = '?';
  if (stripPgmFormHeader(fpIn, pnx, pny, fType) < 0)
    return(-1);
  if (fType[0] != 'F')
    return(-1);
  if ((fType[1] == 'P') || (fType[1] == 'U'))
    *pByteOrderType = fType[1];
  else if (fType[1] == '5')  /* Previously used only on an SGI */
    *pByteOrderType = 'U'; 
  else
    return(-1);
  return(0);
}

/* Read a pgm form header from input file, return the first comment line
   found (when *headerState = 0).  *headerState is modified to indicate
   the state of the pgm header read so far.  If this value of *headerState
   is used in a second call to readPnmCommentLine, then a second comment
   line is read.  And so on.  
   The recognized types are:
     - "P5" - pgm image
     - "P6" - ppm image
     - "FP" - pfm image, image of floats in a PC format
     - "FU","F5" - pfm image, image of floats in a UNIX machine format
   readPnmCommentLine returns:
      -1 if error found
       0 if no comment found
       1 found a comment line
*********/
int readPnmCommentLine(FILE *fpIn, char *commentLine, int *headerState)
{
  char line[MAXLINE_PGM], fileType[MAXLEN], *cPtr;
  int i, k, paramCount, foundComment;

  /* Strip file type in first two characters */
  if (*headerState == 0) {
    fread(line, (int) 2, 1, fpIn);
    line[2] = '\0';
    /* Check if file type is one of the recognized types */
    if (!((strcmp(line, "P5") == 0)
	  || (strcmp(line, "P6") == 0)
	  || (strcmp(line, "FP") == 0)
	  || (strcmp(line, "F5") == 0)
	  || (strcmp(line, "FU") == 0)))
      return(-1);   /* unrecognized file type */
    fread( line, (int) 1, 1, fpIn);
    if (!whtspc(line[0]))
      return(-1);
    *headerState = 1;
  }
  paramCount = (*headerState)-1;

  /* Read parameters nx, ny, ngray */
  foundComment = 0;
  while(paramCount < 3 && foundComment == 0) {  
    /* Seek next parameter or comment */

    /* Strip white space */
    while(whtspc(line[0])) {
      if (! fread(line, (int) 1, 1, fpIn)) 
	return(-1);
    }

    if (line[0] == '#') {

      /* Read comment line */
      k = getStripLine(fpIn, commentLine, MAXLINE_PGM);
      if (k==1) { /* Could not find a '\n' before MAXLINE_PGM,
		     return first portion of comment line, and flush
		     the rest. */
	if (flushLine(fpIn) < 0)
	  return(-1); /* Eof bumped into unexpectedly */
      } else if (k < 0) 
	return(-1); /* Eof bumped into unexpectedly */
      foundComment = 1;  
    } else {

      /* Read parameter */
      i = 0;
      while( !whtspc(line[i]) && i < MAXLINE_PGM) { 
	i++;
	if (!fread( &(line[i]), (int) 1, 1, fpIn)) 
	  return(-1);
      }
      strtol(line, &cPtr, 10);
      if (cPtr == NULL) 
	return(-1);
      paramCount++;
      line[0] = line[i]; /* re-prime whitespace stream */
    }
  }
  *headerState = paramCount+1;
  return(foundComment);
} 

void save_float_image_from_float(char *path, float *image, int nx, int ny)
{
  char dataOut[MAXLEN];
  int sizeOutput;
  FILE *outfile;

  /* Save float image */
  strcpy(dataOut, "");
  strcat(strcat(dataOut, path), fltFileTail);
  if ( (outfile = fopen(dataOut, "wb")) == NULL)
      error("Can't open output file", dataOut);
  sizeOutput = nx * ny * sizeof(float);
  fwrite((char *) image, sizeOutput, 1, outfile);
  fclose(outfile);

}

int range_crop_float_image(float *imageFlt, int nx, int ny, 
			   float crop[2], float range[2])

 /* computes range of values in nx by ny float image, with
    values outside the cropping interval [crop[0], crop[1]] ignored

    Returns -1 if no values found within this range,
            0 for success
 */

 {
     float val;
     int cut, ix, iy;

     cut = 0;
     range[0] = crop[1];  /* Starting `min' value ... as large as possible */
     range[1] = crop[0];  /* Starting `max' value ... as small as possible */

     for( iy=0; iy< ny; iy++) {
      for( ix=0; ix< nx; ix++) {
       
       val = imageFlt[index(ix,iy,nx)];

       if ((val >= crop[0]) && (val <= crop[1]) ) {
        range[0] = min( range[0], val);
        range[1] = max( range[1], val);
       }
       else cut++;

      }  
     }

     if (cut < nx * ny)
      return(0);
     else
      return(-1);  /* No image values within crop interval */

 }

void range_float_image(float *imageFlt, int nx, int ny, float range[2])

 /* computes range of values in nx by ny float image */

 {
     float val;
     int ix, iy;

     range[0] = range[1] = imageFlt[0];

     for( iy=0; iy< ny; iy++) {
      for( ix=0; ix< nx; ix++) {
       
       val = imageFlt[index(ix,iy,nx)];
       range[0] = min( range[0], val);
       range[1] = max( range[1], val);

      }  
     }

 }


void range_byte_image(unsigned char *image, int nx, int ny, 
                      unsigned char range[2])

 /* computes range of values in nx by ny byte image */

 {
     unsigned char val;
     int ix, iy;

     range[0] = range[1] = image[0];

     for( iy=0; iy< ny; iy++) {
      for( ix=0; ix< nx; ix++) {
       
       val = image[index(ix,iy,nx)];
       range[0] = min( range[0], val);
       range[1] = max( range[1], val);

      }  
     }

 }

void quantize_crop_float_image(unsigned char **pimage, float* imageFlt, 
			       int nx, int ny, 
			       float crop[2], float out_of_range[2], float par[2])
 /***
  Params:
    out_of_range[2]:  -0.05 , 1.10  means values below crop[0]
                      are assigned a float value of less than the
		      minimum value by 5% of the range.  While values
		      above crop[1] are assigned a float value of more
		      than the maximum value by 10% of the range.
    par[2]:  return pedestal, scl respectively
 ****/
 {
   float range[2], extent, cropVal[2], pedestal, scl;
   unsigned char *image;
   int ix, iy, ixy;

   if (crop[1] < crop[0])
    error(" Crop interval is null","");
   if ((out_of_range[0] > 0.0 && out_of_range[0] < 1.0)
       || (out_of_range[1] > 0.0 && out_of_range[1] < 1.0) )
    fprintf(stderr,"Warning: cropped image values are aliased\n");

   grabByteMemory((char **) pimage, nx*ny, "quantized image");
   image = *pimage;

   if (range_crop_float_image(imageFlt, nx, ny, crop, range) < 0) {
    fprintf(stderr,"Warning: no image values in crop range %e, %e\n", crop[0], crop[1]);
    range[0] = crop[0]; range[1] = crop[1];
   }
/***
  fprintf(stderr," Range of cropped image values %e, %e\n", range[0], range[1]);
***/
   extent = range[1] - range[0];
   if (extent <= 0.0) 
    extent = 1.0;
   
   /* Compute outlier values */
   cropVal[0] = range[0] + out_of_range[0] * extent;
   cropVal[1] = range[0] + out_of_range[1] * extent;

   /* update range values */
   range[0] = min(cropVal[0], range[0]);
   range[0] = min(cropVal[1], range[0]);
   range[1] = max(cropVal[0], range[1]);
   range[1] = max(cropVal[1], range[1]);
   extent = range[1] - range[0];
   if (extent <= 0.0) 
    extent = 1.0;
   
   par[0] = pedestal = range[0];
   par[1] = scl = 255.0/extent;
   for(iy=0; iy<ny; iy++)
    for(ix=0; ix<nx; ix++) {
     ixy = index(ix,iy,nx);
     image[ixy] = quantVal( imageFlt[ixy], pedestal, scl, crop, cropVal);
    }

 }

void quantize_crop_range_float_image(unsigned char **pimage, float *imageFlt, 
				     int nx, int ny, 
				     float crop[2], float out_of_range[2],
				     float par[2])
 /******
 Params:
   crop[2]  The scaling is forced to be crop[0] to crop[1],
            with, perhaps values designating out of range 
   out_of_range[2]  -0.05 , 1.10  means values below crop[0]
                    are assigned a float value of less than the
		    minimum value by 5% of the range.  While values
		    above crop[1] are assigned a float value of more
		    than the maximum value by 10% of the range. 
   par[2]  return pedestal, scl respectively 
 ******/
 {
   float range[2], extent, cropVal[2], pedestal, scl;
   float middle, val, res;
   unsigned char *image;
   int ix, iy, ixy;

   if (crop[1] < crop[0])
    error(" Crop interval is null","");
   if ((out_of_range[0] > 0.0 && out_of_range[0] < 1.0)
       || (out_of_range[1] > 0.0 && out_of_range[1] < 1.0) )
    fprintf(stderr,"Warning: cropped image values are aliased\n");


   grabByteMemory((char **) pimage, nx*ny, "quantized image");
   image = *pimage;

   range[0] = crop[0];
   range[1] = crop[1];
   extent = range[1] - range[0];
   
   if (extent <= 0.0) 
    extent = 1.0;
   
   /* Compute outlier values */
   cropVal[0] = range[0] + out_of_range[0] * extent;
   cropVal[1] = range[0] + out_of_range[1] * extent;

   /* update range values */
   range[0] = min(cropVal[0], range[0]);
   range[0] = min(cropVal[1], range[0]);
   range[1] = max(cropVal[0], range[1]);
   range[1] = max(cropVal[1], range[1]);
   extent = range[1] - range[0];
   if (extent <= 0.0) 
    extent = 1.0;
   
   par[1] = scl = 256.0/extent;  
   pedestal = range[0];
   /* Make middle map to the center of a bin */
   middle = (crop[0] + crop[1])/2.0;
   val = (middle-pedestal)*scl;
   res = val - ROUND(val);
   pedestal += res/scl;
   par[0] = pedestal;
   /*******
   if (0.0>=crop[0] && 0.0<crop[1])
     fprintf(stderr," quantization sets 0.0 to %f -> %d\n", 
	      (0.0 - pedestal)*scl,
       (int) quantVal((float) 0.0, pedestal, scl, crop, cropVal));
   **********/

   for(iy=0; iy<ny; iy++)
    for(ix=0; ix<nx; ix++) {
     ixy = index(ix,iy,nx);
     image[ixy] = quantVal( imageFlt[ixy], pedestal, scl, crop, cropVal);
    }

 }

void quantize_crop_range_byte_image(unsigned char **pimage, 
   unsigned char *imageIn,
   int nx, int ny, float crop[2], float out_of_range[2], float par[2])
 /****
 Params:
   crop[2] The scaling is forced to be crop[0] to crop[1],
           with, perhaps values designating out of range
   out_of_range[2]  -0.05 , 1.10  means values below crop[0]
                    are assigned a float value of less than the
                    minimum value by 5% of the range.  While values
                    above crop[1] are assigned a float value of more
                    than the maximum value by 10% of the range.
   par[2]  return pedestal, scl respectively
 ******/
 {
   float range[2], extent, cropVal[2], pedestal, scl;
   unsigned char *image;
   int ix, iy, ixy;

   if (crop[1] < crop[0])
    error(" Crop interval is null","");
   if ((out_of_range[0] > 0.0 && out_of_range[0] < 1.0)
       || (out_of_range[1] > 0.0 && out_of_range[1] < 1.0) )
    fprintf(stderr,"Warning: cropped image values are aliased\n");


   grabByteMemory((char **) pimage, nx*ny, "quantized image");
   image = *pimage;

   range[0] = crop[0];
   range[1] = crop[1];
   extent = range[1] - range[0];
   if (extent <= 0.0) 
    extent = 1.0;
   
   /* Compute outlier values */
   cropVal[0] = range[0] + out_of_range[0] * extent;
   cropVal[1] = range[0] + out_of_range[1] * extent;

   /* update range values */
   range[0] = min(cropVal[0], range[0]);
   range[0] = min(cropVal[1], range[0]);
   range[1] = max(cropVal[0], range[1]);
   range[1] = max(cropVal[1], range[1]);
   extent = range[1] - range[0];
   if (extent <= 0.0) 
    extent = 1.0;
   
   par[0] = pedestal = range[0];
   par[1] = scl = 255.0/extent;

   for(iy=0; iy<ny; iy++)
    for(ix=0; ix<nx; ix++) {
     ixy = index(ix,iy,nx);
     image[ixy] = quantVal( (float)imageIn[ixy], pedestal, scl, crop, cropVal);
    }

 }

void quantize_float_image(unsigned char **pimage, float *imageFlt, 
			  int nx, int ny)
 /* Useful for quantizing interpolated images.
    Assumes float values range between 0, 255, so
    pedestal = 0, scl = 1.0.  Values outside this range are
    saturated at 0 or 255 */
 {
   float par[2], crop[2], pedestal, scl;
   unsigned char *image;
   int ix, iy, ixy;


   crop[0] = 0.0;
   crop[1] = 255.0;
   par[0] = pedestal = 0.0;
   par[1] = scl = 1.0;

   grabByteMemory((char **) pimage, nx*ny, "quantized image");
   image = *pimage;

   for(iy=0; iy<ny; iy++)
    for(ix=0; ix<nx; ix++) {
     ixy = index(ix,iy,nx);
     image[ixy] = quantVal( imageFlt[ixy], pedestal, scl, crop, crop);
    }
 }

void reverseContrast(unsigned char *image, int nx, int ny)
 {
  int ixy;

  for(ixy=0; ixy<nx*ny; ixy++)
   image[ixy] = 255 - image[ixy];
 }

unsigned char quantVal( float val, float pedestal, float scl, 
			float crop[2], float cropVal[2])
 {
  float sval;
  int qVal;

  /* Crop */
  if (val < crop[0])
   sval = cropVal[0];
  else if (val > crop[1])
   sval = cropVal[1];
  else 
   sval = val;

  /* Map to standard interval */
  sval = (sval-pedestal)*scl;
  qVal = ROUND(sval);

  /* Crop to [0, 255] */
  qVal = min( 255, qVal);
  qVal = max(0, qVal);
  return((unsigned char) qVal);
 }

void save_byte_image_from_float_scaled(char *path, float *imageFlt, 
       int nx, int ny, float crop[2], float out_of_range[2])
{
  char dataOut[MAXLEN], headOut[MAXLEN];
  int sizeOutput;
  float par[2];
  unsigned char *image;
  FILE *outfile;
 
  if (crop[1] < crop[0]) {
   /* Do not perform cropping */
   range_float_image(imageFlt, nx, ny, crop);
   out_of_range[0] = 0.0;
   out_of_range[1] = 1.0;
  }
  quantize_crop_float_image(&image, imageFlt, nx, ny, crop, out_of_range, par);
      
  strcpy(headOut, "");
  strcat(strcat(headOut, path), headerFileTail);
  
  if ((outfile = fopen(headOut, "wb")) == NULL)
    error("Unable to open file", headOut);

  fprintf(stderr, "%d %d %f %f %f %f %f\n", nx,ny, par[0], par[0]+255.0/par[1],
        par[1], crop[0], crop[1]);
  fprintf(outfile, "%d %d %f %f %f %f %f\n", nx,ny, par[0], par[0]+255.0/par[1],
        par[1], crop[0], crop[1]);
  fclose(outfile);

  /* Save byte image */
  strcpy(dataOut,"");
  strcat(strcat(dataOut, path), byteFileTail);
  if ( (outfile = fopen(dataOut, "wb")) == NULL) 
      error("Can't open output file", dataOut);

  sizeOutput = nx * ny * sizeof( char );
  if (! fwrite((char *) image, sizeOutput, 1, outfile))
   error("Unable to write", dataOut);
  fclose(outfile);

  utilFree((void **) &image);
}

void save_byte_image_from_float_unscaled(char *path, float *imageFlt, 
					 int nx, int ny)
{
  char dataOut[MAXLEN], headOut[MAXLEN];
  int sizeOutput;
  float crop[2], par[2];
  unsigned char *image;
  FILE *outfile;
 
  quantize_float_image(&image, imageFlt, nx, ny);
  par[0] = 0.0;
  par[1] = 1.0;
  crop[0] = 0.0;
  crop[1] = 255.0;
      
  strcpy(headOut, "");
  strcat(strcat(headOut, path), headerFileTail);
  
  if ((outfile = fopen(headOut, "wb")) == NULL)
    error("Unable to write file", headOut);

  fprintf(stderr, "%d %d %f %f %f %f %f\n", nx,ny, par[0], par[0]+255.0/par[1],
        par[1], crop[0], crop[1]);
  fprintf(outfile, "%d %d %f %f %f %f %f\n", nx,ny, par[0], par[0]+255.0/par[1],
        par[1], crop[0], crop[1]);
  fclose(outfile);

  /* Save byte image */
  strcpy(dataOut,"");
  strcat(strcat(dataOut, path), byteFileTail);
  if ( (outfile = fopen(dataOut, "wb")) == NULL ) 
      error("Can't open output file", dataOut);

  sizeOutput = nx * ny * sizeof( char );
  if (! fwrite((char *) image, sizeOutput, 1, outfile))
   error("Unable to write", dataOut);
  fclose(outfile);

  utilFree((void **) &image);
}

FILE * write_pgm_image(FILE *outfile, char *path,
                     unsigned char *image, int nx, int ny)
{
  char dataOut[MAXLEN], header[MAXLEN];
  int  sizeOutput, nHead;
      
  strcpy(dataOut, "");
  sprintf(dataOut,"%d %d\n%d\n", nx, ny, 255);
  strcpy(header, "P5\n");
  nHead = strlen(header) + strlen(dataOut);
  while(nHead < minHeadSize) {
   strcat(header, blank);
   nHead++;
  }
  strcat(header, dataOut);
  
  /* Open image */
  strcpy(dataOut,"");
  strcat(strcat(dataOut, path), pgmFileTail);
  if ( (outfile = freopenOutputFile(dataOut, outfile)) == NULL ) 
      error("Can't open output file", dataOut);

  /* Save header */
  sizeOutput = strlen(header) * sizeof( char );
  if (!fwrite(header, sizeOutput, 1, outfile))
   error("Unable to write", dataOut);

  /* concatenate byte image in output image file */
  sizeOutput = nx * ny * sizeof( char );
  if (!fwrite((char *) image, sizeOutput, 1, outfile))
   error("Unable to write", dataOut);

  fclose(outfile);
  return(outfile);
}

void save_pgm_image(char *path, unsigned char *image, int nx, int ny)
{
  char dataOut[MAXLEN], header[MAXLEN];
  int  sizeOutput, nHead;
  FILE *outfile; 
  
  write_pgm_image((FILE *) NULL, path, image, nx, ny);
}

void save_pgm_image_wcomment(char *path, char *comment,
                             unsigned char *image, int nx, int ny)
{
  char dataOut[MAXLEN], header[MAXLEN];
  int  sizeOutput, nHead;
  FILE *fpOut;
 
  sprintf(dataOut,"%d %d\n%d\n", nx, ny, 255);
  strcpy(header, "P5\n#");
  strcat(header, comment);

  nHead = strlen(header) + strlen(dataOut);
  while(nHead < minHeadSize) {
   strcat(header, blank);
   nHead++;
  }
  strcat(header, dataOut);
  
  /* Open image */
  strcpy(dataOut,"");
  strcat(strcat(dataOut, path), pgmFileTail);
  if ( (fpOut = fopen(dataOut,"wb")) == NULL  ) 
      error("Can't open output file", dataOut);

  /* Save header */
  sizeOutput = strlen(header) * sizeof( char );
  if ( ! fwrite(header, sizeOutput, 1, fpOut))
   error("Unable to write", dataOut);

  /* concatenate uchar image in output image file */
  sizeOutput = nx * ny * sizeof( char );
  if ( ! fwrite((char *) image, sizeOutput, 1, fpOut))
   error("Unable to write", dataOut);

  fclose(fpOut);
}

void save_ppm_image(char *path, unsigned char *cImage[3], int nx, int ny)
{
  unsigned char *raw;
  char dataOut[MAXLEN], header[MAXLEN];
  int sizeOutput, nHead, ixy, k;
  FILE *outfile;
 
  grabByteMemory((char **) &raw, 3 * nx * ny, "rawColourImage");   

  strcpy(dataOut, "");
  sprintf(dataOut,"%d %d\n%d\n", nx, ny, 255);
  strcpy(header, "P6\n");
  nHead = strlen(header) + strlen(dataOut);
  while(nHead < minHeadSize) {
   strcat(header, blank);
   nHead++;
  }
  strcat(header, dataOut);
  
  /* Open image */
  strcpy(dataOut,"");
  strcat(strcat(dataOut, path), ppmFileTail);
  if ( (outfile = fopen(dataOut, "wb")) == NULL  ) 
      error("Can't open output file", dataOut);

  /* Save header */
  sizeOutput = strlen(header) * sizeof( char );
  if ( ! fwrite(header, sizeOutput, 1, outfile))
   error("Unable to write", dataOut);

  /* concatenate RBG byte images */
  for(ixy=0; ixy<nx * ny; ixy++) {
   k = 3*ixy;
   raw[k] = cImage[0][ixy];
   raw[k+1] = cImage[1][ixy];
   raw[k+2] = cImage[2][ixy];
  }
  /* write raw byte image in output image file */
  sizeOutput = 3 * nx * ny * sizeof( char );
  if ( ! fwrite((char *) raw, sizeOutput, 1, outfile))
   error("Unable to write", dataOut);

  fclose(outfile);
  utilFree((void **) &raw);
}

void save_pfm_image(char *path, float *image, int nx, int ny)
{
  char dataOut[MAXLEN], header[MAXLEN];
  int sizeOutput, nHead;
  FILE *fpOut;

  if (BYTE_ORDER_TYPE == '?')
    error("save_pfm_image: Unknown byte order type.","");
      
  strcpy(dataOut, "");
  sprintf(dataOut,"%d %d\n%d\n", nx, ny, 0);
  sprintf(header, "F%c\n", BYTE_ORDER_TYPE);
  nHead = strlen(header) + strlen(dataOut);
  while(nHead < minHeadSize) {
   strcat(header, blank);
   nHead++;
  }
  strcat(header, dataOut);
  
  /* Open image */
  strcpy(dataOut,"");
  strcat(strcat(dataOut, path), pfmFileTail);
  if ( (fpOut = fopen(dataOut, "wb")) == NULL  ) 
      error("Can't open output file", dataOut);

  /* Save header */
  sizeOutput = strlen(header) * sizeof( char );
  if (!fwrite(header, sizeOutput, 1, fpOut))
   error("Unable to write", dataOut);

  /* concatenate float image in output image file */
  sizeOutput = nx * ny * sizeof( float );
  if (!fwrite((char *) image, sizeOutput, 1, fpOut))
   error("Unable to write", dataOut);

  fclose(fpOut);

}

void save_pfm_image_wcomment(char *path, char *comment, 
                             float *image, int nx, int ny)
{
  char dataOut[MAXLEN], header[MAXLEN];
  int sizeOutput, nHead;
  FILE *fpOut;
      
  if (BYTE_ORDER_TYPE == '?')
    error("save_pfm_image: Unknown byte order type.","");
      
  strcpy(dataOut, "");
  sprintf(dataOut,"%d %d\n%d\n", nx, ny, 0);
  sprintf(header, "F%c\n#", BYTE_ORDER_TYPE);
  strcat(header, comment);

  nHead = strlen(header) + strlen(dataOut);
  while(nHead < minHeadSize) {
   strcat(header, blank);
   nHead++;
  }
  strcat(header, dataOut);
  
  /* Open image */
  strcpy(dataOut,"");
  strcat(strcat(dataOut, path), pfmFileTail);
  if ( (fpOut = fopen(dataOut,"wb")) == NULL  ) 
      error("Can't open output file", dataOut);

  /* Save header */
  sizeOutput = strlen(header) * sizeof( char );
  if (!fwrite(header, sizeOutput, 1, fpOut))
   error("Unable to write", dataOut);

  /* concatenate float image in output image file */
  sizeOutput = nx * ny * sizeof( float );
  if (!fwrite((char *) image, sizeOutput, 1, fpOut))
   error("Unable to write", dataOut);

  fclose(fpOut);

}

void save_byte_image(char *path, unsigned char *image, int nx, int ny)
{
  char dataOut[MAXLEN];
  int sizeOutput;
  FILE *outfile;
      
  strcpy(dataOut, "");
  strcat(strcat(dataOut, path), byteFileTail);
  if ( (outfile = fopen(dataOut, "wb")) == NULL ) 
      error("Can't open output file", dataOut);

  sizeOutput = nx * ny * sizeof( char );
  if ( ! fwrite((char *) image, sizeOutput, 1, outfile))
   error("Unable to write", dataOut);

  fclose(outfile);

}

void writeFilterOutputs(char *fname, char *fidentifier, 
 int nByte, int nOrient, int nDim, int nx, int ny, int numPar, 
 float *outPar, float **resp)
 /****************************************************************
  Writes filter response data onto file "fname"
    Format:  Assumes data has nOrient different channels, and each
             channel is of dimension nDim  (eg. a complex filter response
             could have 4 orientations each of which is of dimension 2)
             Data to be stored with nByte's per pixel.
             (so far only float output implemented, nByte==4)
    
    Order:
             global header  ( size nHeader )
             header orientation 1 ( size nHeader )
             orientation 1, dimension 1
             orientation 1, dimension 2
             ... 
             orientation 1, dimension nDim
             header orientation 2 (size nHeader)
             orientation 2, dimension 1
             ...
  *********************************************************************/      
      
 {
    float pedestal, scl;
    int iOrient, iDim, k, sizeOutput, iChannel;
    char header[HEADER], scrPar[32];
    FILE *outfile;

    if ( (outfile = fopen(fname, "wb")) == NULL ) 
      mexErrMsgTxt("Can't open output file");

    pedestal = 0.0; scl = 1.0;
    sprintf(header,"%s %d %d %d %d", "float", nx, ny, nOrient, nDim); 
    if ( ! fwrite((char *) header, HEADER, 1, outfile))
      mexErrMsgTxt("Unable to write to output file");

    sizeOutput = nx * ny * sizeof(float);
    for(iOrient=0; iOrient<nOrient; iOrient++) {
     sprintf(header,"%s %e %e", fidentifier,
       pedestal, scl);
     for(k=0; k<numPar; k++) {
      sprintf(scrPar, " %e", outPar[k+iOrient*numPar]);
      strcat(header, scrPar);
     }
     if ( ! fwrite((char *) header, HEADER, 1, outfile))
       mexErrMsgTxt("Unable to write to output file");
     for(iDim=0; iDim<nDim; iDim++) {
      iChannel = index(iDim, iOrient, nDim);
      if (!fwrite((char *) resp[iChannel], sizeOutput, 1, outfile))
	mexErrMsgTxt("Unable to write to output file");
     }
    }
    fclose(outfile);
 }

unsigned char colour(float v, float minV, float scl)
 {
  float val;
  int iVal;
  val = (v-minV)*scl;
  iVal = ROUND(val);
  iVal = min( 255, iVal);
  iVal = max(0, iVal);
  return((unsigned char) iVal);
 }

void swapWordByteOrder( void *p)
 {
  char bytes[4];
  char *a = (char *)p;
  int k;

  for(k=0; k<4; k++)
    bytes[3-k] = a[k];
  for(k=0; k<4; k++)
    a[k] = bytes[k];
 }

void floatImageSwapBytes(float *image, int nx, int ny)
{
 int i;
 float *p = image;
 for(i=0; i<nx*ny; i++, p++)
  swapWordByteOrder((void *)p);
}

/* Expand an image by duplicating each pixel to 
   a patch of size xpnd x xpnd  */
void expandImage(unsigned char **pBigImage, unsigned char *image,
		 int nx, int ny, int xpnd)
{
 unsigned char val;
 int ix, iy, bx, by, i,j;
 int bnx, bny;
 bnx = xpnd*nx;
 bny = xpnd*ny;
 grabByteMemory((char **)pBigImage, bnx*bny, "big image");
 for(iy=0; iy<ny; iy++) {
   by = xpnd*iy;
   for(ix=0; ix<nx; ix++) {
     bx = xpnd*ix;
     val = image[index(ix,iy,nx)];
     for(j=0; j<xpnd; j++)
       for(i=0; i<xpnd; i++) 
	 (*pBigImage)[index(bx+i, by+j, bnx)] = val;
   }
 }
}


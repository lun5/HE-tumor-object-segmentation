
/************************* Version 2.0 **************************/
# include <stdio.h>

#ifndef _IMAGEFILE_IO_H
# define _IMAGEFILE_IO_H

# define        HEADER      	512
# define        NULL_NAME_FLAG   -2
# define	MAXLINE_PGM	256  /* specification (71) not always obeyed */
# define	MAXSIZE_PGM	4096

#endif

/* Return the byte order type 'P', 'U', '?' for the current system */
char getByteOrderType();

/* Read image filename, image type, and necessary size information
   Returns *pPgmInput = 1 if file is a pgm file, otherwise
                      = 0, file is assumed to be raw bytes */
int read_image_params(char *fnInput, int *pPgmInput,
  int *pNxImage, int *pNyImage, int *pNHead);

int read_file_path_and_root(char *pathName, char *fileRoot);
  
/******************
fopenInputFile:
   Attempt to open input image file with the name fnIn
   Returns fpInput if successful, NULL otherwise
************/
FILE* fopenInputFile(char *fnIn);

/******************
freopenInputFile:
   Attempt to reopen input image file with the name fnIn
   on an existing stream fp or, if fp==NULL, then it
   attempts to open a new stream.
   Returns fpInput if successful, NULL otherwise
******************/
FILE* freopenInputFile(char *fnIn, FILE * fp);

/******************
freopenOutputFile:
   Attempt to reopen output image file with the name fnIn
   on an existing stream fp or, if fp==NULL, then it
   attempts to open a new stream. 
   Returns fpOutput if successful, NULL otherwise
******************/
FILE * freopenOutputFile(char* fnOut, FILE* fp);

/******************
fopenInputFileNum:
   Attempt to open input image file with the name of the form:
    filename = headIn0..0###TailIn
   Where the image number 0...0### is possibly padded with leading
   zeroes to have a total length of at most MAXPAD.
   
   Routine openInputFileNum returns FILE* infile if successful, NULL otherwise
************/
FILE* fopenInputFileNum(char* fnIn, char* headIn, 
                       int* ppad, int kImage, char* tailIn);

/******************
build_filename:
   Construct filename of the form:
    filename = head###tail
   Where there are at least `pad' digits in the file number, padded with
   leading zeroes if needed.

************/
void build_filename(char *fname, char *head, 
		    int pad, int numImage, char *tail);

/* return an integer as a string, padded with leading
   zeroes to make the string at least of length pad. */
void pad_integer(char *padInt, int n, int pad);

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
 int *pFileNum, int maxNum); 

void read_byte_image_strip_bytes(unsigned char* image, FILE* fpIn, char* fnIn,
                            int nx, int ny, int bytes);

void read_byte_image(unsigned char* image, FILE* fpIn, char* fnIn, 
                     int nx, int ny);

void read_float_image_strip_bytes(float* imageFlt, FILE* fpIn, char* fnIn,
				  int nx, int ny, int bytes);

void read_float_image(float* imageFlt, FILE* fpIn, char* fnIn, int nx, int ny);

void read_int_image_strip_bytes(unsigned int * image, FILE* fdIn, char * fnIn,
				int nx, int ny, int bytes);

void read_int_image(unsigned int *image, FILE *fpIn, char *fnIn, int nx, int ny);

void int_to_float_image(float **pimage, unsigned int *intImage, int nx, int ny);

void load_image(char *path, float **pimage, int *nx, int *ny);

int read_pgm_image(unsigned char **pimage, FILE *fpIn, char *fnIn, 
		   int *pnx, int *pny);

int read_ppm_image(unsigned char *pimage[3], FILE *fpIn, char *fnIn, 
		   int *pnx, int *pny);

int read_pfm_image(float **pimage, FILE *fpIn, char *fnIn, 
		   int *pnx, int *pny);

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
             headers and pointers to memory blocks (in above order)
  *********************************************************************/      
void readFilterOutputs(FILE *fpIn, char* fname, char header[][HEADER], 
		       int nHeader, int nByte, int nOrient, int nDim, 
		       int sizeImage, 
		       char **image, int *pnxImage, int *pnyImage);

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
int stripPgmFormHeader(FILE *fpIn, int *pnx, int *pny, char fType[3]);
/* Same function */
int stripPnmHeader(FILE *fpIn, int *pnx, int *pny, char fType[3]);

int stripPgmHeader(FILE *fpIn, int *pnx, int *pny);
 
int stripPpmHeader(FILE *fpIn, int *pnx, int *pny);

/* Reads the header from a float pgm form image (pfm)
   Returns nx, ny, and 
     SysType = 'U' if the image was written on a Unix box
             = 'P' if the image was written on a PC
             = '?' otherwise
**********/
int stripPfmHeader(FILE *fpIn, int *pnx, int *pny, char *pSysType);

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
int readPnmCommentLine(FILE *fpIn, char *commentLine, int *headerState);

void save_float_image_from_float(char *path, float *image, int nx, int ny);

 /* computes range of values in nx by ny float image, with
    values outside the cropping interval [crop[0], crop[1]] ignored

    Returns -1 if no values found within this range,
            0 for success
 */
int range_crop_float_image(float *imageFlt, int nx, int ny, 
			   float crop[2], float range[2]);


 /* computes range of values in nx by ny float image */
void range_float_image(float *imageFlt, int nx, int ny, float range[2]);

 /* computes range of values in nx by ny byte image */
void range_byte_image(unsigned char *image, int nx, int ny, 
                      unsigned char range[2]);

 /***
  Params:
    out_of_range[2]:  -0.05 , 1.10  means values below crop[0]
                      are assigned a float value of less than the
		      minimum value by 5% of the range.  While values
		      above crop[1] are assigned a float value of more
		      than the maximum value by 10% of the range.
    par[2]:  return pedestal, scl respectively
 ****/
void quantize_crop_float_image(unsigned char **pimage, float* imageFlt, 
		    int nx, int ny, 
		    float crop[2], float out_of_range[2], float par[2]);

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
void quantize_crop_range_float_image(unsigned char **pimage, float *imageFlt, 
				     int nx, int ny, 
				     float crop[2], float out_of_range[2],
				     float par[2]);

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
void quantize_crop_range_byte_image(unsigned char **pimage, 
   unsigned char *imageIn,
   int nx, int ny, float crop[2], float out_of_range[2], float par[2]);

 /* Useful for quantizing interpolated images.
    Assumes float values range between 0, 255, so
    pedestal = 0, scl = 1.0.  Values outside this range are
    saturated at 0 or 255 */
void quantize_float_image(unsigned char **pimage, float *imageFlt, 
			  int nx, int ny);

void reverseContrast(unsigned char *image, int nx, int ny);

unsigned char quantVal( float val, float pedestal, float scl, 
			float crop[2], float cropVal[2]);

void save_byte_image_from_float_scaled(char *path, float *imageFlt, 
       int nx, int ny, float crop[2], float out_of_range[2]);

void save_byte_image_from_float_unscaled(char *path, float *imageFlt, 
					 int nx, int ny);

FILE * write_pgm_image(FILE *outfile, char *path,
                     unsigned char *image, int nx, int ny);

void save_pgm_image(char *path, unsigned char *image, int nx, int ny);

void save_pgm_image_wcomment(char *path, char *comment,
                             unsigned char *image, int nx, int ny);

void save_ppm_image(char *path, unsigned char *cImage[3], int nx, int ny);

void save_pfm_image(char *path, float *image, int nx, int ny);

void save_pfm_image_wcomment(char *path, char *comment,
                             float *image, int nx, int ny);

void save_byte_image(char *path, unsigned char *image, int nx, int ny);

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
 void writeFilterOutputs(char *fname, char *fidentifier, 
 int nByte, int nOrient, int nDim, int nx, int ny, int numPar, 
 float *outPar, float **resp);

unsigned char colour(float v, float minV, float scl);

void swapWordByteOrder( void *p);

void floatImageSwapBytes(float *image, int nx, int ny);

/* Expand an image by duplicating each pixel to 
   a patch of size xpnd x xpnd  */
void expandImage(unsigned char **pBigImage, unsigned char *image,
		 int nx, int ny, int xpnd);


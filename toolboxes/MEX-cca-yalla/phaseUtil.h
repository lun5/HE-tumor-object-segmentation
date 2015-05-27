/* FILE: phaseUtil.h 
   For phaseUtil.V4.0.c
   May 1999 */

#ifndef _PHASE_H
# include "phase.h"  /* Provides constants, NUMSCALE ... */
#endif

#ifndef _FREEMAN_PYRAMID_H
# include "freemanPyramid.h"  /* Provides freemanPyramid structure typedefs */
#endif

/**********************************************************
  Function prototypes for phaseUtil.V4.*.c
**********************************************************/

/***
   Reset the global tunable parameters according to values
   in the FreemanPyramidCntrlStruct. 
   DO NOT use this unless you know what you are doing. */
void resetFreemanPyramidParams(FreemanPyramidCntrlStruct * pc);
                         
/***
 recoverFastFilterOutputs:
   The input param bytesPerPixel must be 1 or 4, indicating that 
   the input image is of the type (unsigned char *) or (float *), respectively.
   When bytesPerPixel == 4, and *localPyrBase is returned to be 0 (false),
   then the input float image is being used for the base level of
   the pyramid.  That is, the float * variable pyramid[0] is aliased
   with the user's input image.  So be careful not to free pyramid[0] when
   *localPyrBase is returned to be 0.
  See filterDemo.c and freemanPyramid.c for an example of the use of this
  routine.
  General Idea:
    For each scale, this routine checks what information
    is desired, according to the storeMap array.  It 
    figures out what needs to be computed (computeMap[]).
    It checks the cache for what is needed, if the data isn't
    found in the cashe, then it computes whatever is needed.
    New computations are cached when indicated by the global array
    cacheMap[] (see phaseUtil.V4.0.c:defaultFreemanCacheCntrl()).
 ***/
void recoverFastFilterOutputs(void * image, int imageBytesPerPixel,
	int nxImage, int nyImage,
	char *pathNameCache, char *rootNameCache, int *cacheFiles,
	float *lambda, int nScale, 
	float *thetaDeg, int nOrient,
	int *storeMap,
	float **pyrMap, int *nPyr, int *localPyrBase,
	float **basisMap, float **filterMap, 
	float **pMap, float **aMap, float **gMap,
	float epsAmp, float epsC, 
	int *nxPyr, int *nyPyr, 
	int *nxFilt, int *nyFilt);

/* Return the Freeman pyramid maps specifed by the array of flags, storeMap,
   computed for the given image.  
   The input param bytesPerPixel must be 1 or 4, indicating that 
   the input image is of the type (unsigned char *) or (float *), respectively.
   When bytesPerPixel == 4, and *localPyrBase is returned to be 0 (false),
   then the input float image is being used for the base level of
   the pyramid.  That is, the float * variable pyramid[0] is aliased
   with the user's input image.  Be careful not to free pyramid[0] when
   *localPyrBase == 0.
   Intermediate maps may be computed to form the ones that are returned, 
   but these are freed by this routine. */
void computeFastFilterOutputs(void * image, int bytesPerPixel,
        int nxImage, int nyImage,
        float *lambda, int nScale,
        float *thetaDeg, int nOrient,
        int *storeMap,
        float **pyrMap, int *nPyr, int *localPyrBase,
        float **basisMap, float **filterMap,
        float **pMap, float **aMap, float **gMap,
        float epsAmp, float epsC,
        int *nxPyr, int *nyPyr,
        int *nxFilt, int *nyFilt);


 /*******************************************************
  recoverPyramid:
   The input param bytesPerPixel must be 1 or 4, indicating that 
   the input image is of the type (unsigned char *) or (float *), respectively.
   When bytesPerPixel == 4, and *localPyrBase is returned to be 0 (false),
   then the input float image is being used for the base level of
   the pyramid.  That is, the float * variable pyramid[0] is aliased
   with the user's input image.  Be careful not to free pyramid[0] when
   *localPyrBase == 0.
  General Idea:
    This routine computes an image pyramid with nPyr levels,
    subsampling the original image at 1, 2, 4, ... 2^(nPyr-1)
    pixels, respectively.
    For each level, this routine checks the cache. It then
    figures out what needs to be computed.
    New computations are cached when the flag cachePyr > 0
 *******************************************************/
void recoverPyramid(void * image, int imageBytesPerPixel,
		    int nxImage, int nyImage,
		    char *pathNameCache, char *rootNameCache, 
		    int nPyr, int *localPyrBase,
		    float **pyramid, int *nxPyr, int *nyPyr, 
		    int cachePyr);

/* Free allocated images for the current set of filter responses.
   Here storeMap[isPyr] !=0 indicates the pyramid images are to be freed,
   and similarly for the basis and various filter responses.
   The scaleRange array gives the range of scales to be freed,
   to free all scales, set scaleRange = (0, nScale).  Here scaleRange[1]
   is NOT deleted...just iScale = scaleRange[0]..(scaleRange[1]-1).
   Finally he lowest level in the pyramid (pyrMap[0]) is treated differently,
   since this may refer to a user's float image (indicated by localPyrBase == 0).
   In that case we take care NOT to free it here.  The user needs to do that. */
void freeFastFilterOutputs(int scaleRange[2], int nOrient,
			   int *storeMap, 
			   float **pyrMap, int nPyr, int localPyrBase,
			   float **basisMap, float **filterMap, 
			   float **pMap, float **aMap, float **gMap);

/* computePyramid: same as recoverPyramid, only does not use cache. 
   The input param bytesPerPixel must be 1 or 4, indicating that 
   the input image is of the type (unsigned char *) or (float *), respectively.
   When bytesPerPixel == 4, and *localPyrBase is returned to be 0 (false),
   then the input float image is being used for the base level of
   the pyramid.  That is, the float * variable pyramid[0] is aliased
   with the user's input image.  Be careful not to free pyramid[0] when
   *localPyrBase == 0.
*/
void computePyramid(void *image, int bytesPerPixel, 
	       int nxImage, int nyImage, 
               int nPyr, int *localPyrBase, 
	       float **pyramid, int *nxPyr, int *nyPyr);


/** Return just one level of the pyramid. Provides simpler
    interface to computePyramid.  (Not called within phaseUtil.c) **/
void samplePyramid(void *image, int bytesPerPixel, 
		   int nxImage, int nyImage, 
		   int pyrLevel, int *localPyrBase,
		   float **pDwnImageFlt, int *pnxDwn, int *pnyDwn);

/* getGradMap:  compute phase gradient images, PhiX and PhiY, 
   given filter responses for a given scale and orientation. */
void getGradMap(float *filtResp[2], int nxFilt, int nyFilt, 
	   float lambda, float thetaDeg, 
           float **pPhiX, float **pPhiY, float *pMaxAmp);

/* getGradMap:  compute phase image, *pMap, 
   given filter responses for a given scale and orientation. */
void getPhaseMap(float **pMap, float **filtMap, 
		 int nx, int ny);

/* getAmpMap:  compute amplitude image, *aMap, 
   given filter responses for a given scale and orientation.
   Small amplitudes (below rockBottom) are censored (set
   to BIG_NEG. */
void getAmpMap(float **aMap, float **filtMap, int nx, int ny, float *pMax);

/* getFastFilterMaps:  compute G2,H2 filter responses for a given scale
   and orientation given the seven basis images. */
void getFastFilterMaps(float **filtResp, float thetaDeg, 
		  float **respG2, float **respH2, 
		  int nxFilt, int nyFilt, float *pMaxAmp);

/* steerFilterMaps:  given a steering orientation image, thetaDeg[],
   compute G2,H2 filter responses for a given scale
   at the steering orientation for each corresponding pixel
   given the seven G2,H2 basis images. */
void steerFilterMaps(float **filtResp, float *thetaMap, 
		  float **respG2, float **respH2, 
		  int nxFilt, int nyFilt, float *pMaxAmp);

/* ampFilterMaps:  compute amplitude image, *pAmp, 
   given filter response images filtResp[2] for a given scale and
   orientation.  Small amplitudes are NOT censored. */
void ampFilterMaps(float **pAmp, float **filtResp,
		   int nxFilt, int nyFilt, float *pMaxAmp);

/* Get orientation energy coefficience C1, C2 and C3 (see p.905 Freeman and
   Adelson, Table XI, with sign of theta reversed) */
void getOrientEng(float **basisG2, float **basisH2, int nx, int ny, 
		  float epsContrast, float epsAmp, 
		  float **pUnifMap, float **pAmpMap, float **pThetaMap);

/* Get G2/H2 basis maps using 1D convolution. */
void getFastBasisMaps(void *image, int bytesPerPixel, int dwnSmplImage,
		      int nxImage, int nx0, int ny0, float lambda, 
		      float **respG2, float **respH2,
		      int *pnxFilt, int *pnyFilt);

/* Downsample an image using a Gaussian blur prefilter */
void downSample(void *image, int bytesPerPixel, 
		int nxImage, int nx0, int ny0, int dwnSmpl,
                float **resp, int *pnxFilt, int *pnyFilt);

/* censor an image by doing an OR of a given image with a previously
   censored mask image.  That is, given the mask function maskFlt, and
   the image imageFlt, set any pixel ixy of imageFlt to BIG_NEG for which
   maskFlt[ixy] <= BIG_NEG. */
void censor(float * imageFlt, float *maskFlt, int nx, int ny);

/* Keep only the responses for which the phase gradient is closest to
   the filter tuning.   That is,set any pixel ixy of imageFlt to BIG_NEG
   for which  gradPhase is BIG_NEG or closer to another channel. */
void keepNearestFilter(float *ampMap[], float *gradMap[], 
                       int nxFilt[], int nyFilt[],
                       float thetaDeg[], int nOrient,
                       float lambda[], int nScale,
                       int subSample[]);

/* Perform a gaussian blur on an input image. */
void gaussianBlur(void *image, int bytesPerPixel, 
		  int nxImage, int nx0, int ny0,
		  float sigma, int dwnSmpl, 
                  float **resp, int *pnxFilt, int *pnyFilt);


/* Get G2/H2 filter maps for a given scale and orientation using
   2D correlation.  This is slower than getFastFilterMaps, and should
   only be used for a sanity check. */
void getFilterMaps(void *image, int bytesPerPixel, 
		   int nx0, int ny0, 
		   float lambda, float thetaDeg,
                   int *pnxFiltered, int *pnyFiltered, int *pSubSample,
                   float **filtResp, float *pMaxAmp);

/* Correlate input image with a 2D filter mask, and downsample the result */
void corrByte(unsigned char *input, int dwnSmpl, int nxIn,
	      int nx, int ny, int sample,
	      float *filter, int filter_size,
	      float *output, int nxOut,
	      float *bound);
 
/* Correlate a float image with a 2D filter mask, and downsample the result */
void corrFlt(float *input, int dwnSmpl, int nxIn,
	      int nx, int ny, int sample,
	      float *filter, int filter_size,
	      float *output, int nxOut,
	      float *bound);

/* Separable correlations on an image, and downsample the result */
void corrSepByte(unsigned char *input, int dwnSmpl,int nxIn,
	    int nx0, int ny0, int sample,
	    float *filter_X, float *filter_Y, int filter_size,
	    float *temp, float *output,int  nxOut);

/* Separable correlations on a float image, and downsample the result */
void corrSepFlt(float *input, int dwnSmpl,int nxIn,
	    int nx0, int ny0, int sample,
	    float *filter_X, float *filter_Y, int filter_size,
	    float *temp, float *output,int  nxOut);

/* Do non-max suppression on amp Map */
void thinAmpMap(unsigned char **indMap, float radius, 
		float **aMap, float *thetaDeg, int nOrient, 
		int nx, int ny);

/* Compute the maximum height of the pyramid needed for
   a given set of wavelengths */
int pyramidHeight(float *lambda, int nScale);

/* Compute the subsampling and downsampling rates for an
   array of scales, and the height of the minimum pyramid. 
   Here:
      downSample[iScale]: is the downSampling rate of the low-pass pyramid
                          used for the filtering wrt the original image.
      subSample[iScale]: is the subsampling rate wrt to the downsampled
                         low-pass pyramid level used when filtering for the
                         G2/H2 responses.
      downSample[iScale]*subSample[iScale] :
                        is the subsampling rate for
                        the G2/H2 filters at scale iScale wrt to
                        the original image.
      lowPassSample[iPyr] : subsampling rate (wrt the original image)
                        for this level of the low-pass pyramid.
*/
void pyramidSubSample(int *downSample, int *subSample, 
                     int *lowPassSample, int *pnPyr,
                     float *lambda, int nScale);

/* Use bilinear interpolation to interpolate image map at x,y */
float interpMap(float *fMap, float x, float y, int nx, int ny);

/* Check the cache for existing freeman pyramid images */
void checkPhaseCache(int computeMap[NUMCODE][NUMSCALE], int *storeMap,
		     float *lambda, int nScale, 
		     float *thetaDeg, int nOrient,
		     float *allOrientMaxAmp,
		     char *pathNameCache, char *rootNameCache, 
		     float *respPnt[NUMCODE][MAXLEN],
		     int *nxFilt, int *nyFilt);

/* Check the cache for existing Gaussian pyramid images */
void checkPyramidCache(float **pyramid, int *localPyrBase, 
		       int *computePyr, 
		       int nPyr, int *nxPyr, int *nyPyr, 
		       char *pathNameCache, char *rootNameCache);

/* Compute params of Freeman filter given the peak wavelength */
void freemanConstants(float lambda, int *pdwnSmpl, int *psubSample, 
		      float *psigma, float *prho);

/* Correlate (slow 2D method) a Freeman filter with an image */
void applyFreeman(void *image, int bytesPerPixel, 
		  int nxImage, int nx0, int ny0, 
		  float lambda, float thetaDeg, 
		  float *respR, float *respI, int nxResp, 
		  float *pMaxAmp);

/* Allocate and generate separable Freeman filters */
void build_sepFreeman_filter(float **G2x, float **G2y, 
			     float **H2x, float **H2y, 
			     int *pfilt_size, 
			     float lambda, float sigma, float rho);

/* Allocate and generate a 2D Freeman filter */
void build_freeman_filter(float **pFreemanR, float **pFreemanI, 
			  int *pfilter_size, 
                          float thetaDeg, float lambda,
			  float sigma, float rho);

/* Sample freeman pair G2/H2 at a tap position x[2] */ 
void freeman(float *x, float sigma, float theta, float *pair);

/* Compute the 5-point discrete derivative operator for
   a band-pass filter channel tuned for peak wavelength lambda
   and orientation thetaDeg.  Also estimate the covariance matrix C
   for the frequency response. */
void derivMask(float dPhi[2], float C[2][2], 
	       float maskX[5][2], float maskY[5][2],
	       float lambda, float thetaDeg,
	       float sigma, float rho, int subSample); 
   
/***************************************
 imageLocalFreqFlt:  Input
   nx, ny, gCos, gSin -
          float arrays of size nx by ny (raster order)
          gCos and gSin, which contain real and imaginary part of a
          Hilbert transform pair of filters applied to a real image.
   subSample - subSampling rate of bandpass convolution filters forming
               gCos and gSin
   maskX, maskY - derivative mask coefficients (demodulate differentiate)
                  see derivMask routine above.
   C - 2x2 covariance matrix for filter envelope
   lambda - wavelength for basic filter (in terms of pixels in original
           image)
   thetaDeg - spatial orientation of filter (degrees)
   epsAmp - "rock bottom" amplitude, below which phase etc not computed.
   tau - singularity threshold
   border - size of border to leave out around outside of image
   
 Output:  
   phiX, phiY - float images of the X and Y derivatives of phase
                Gradient in units radians/pixel divided by PI 
                Phase singularities marked by phiX=phiY= BIG_NEG
   maxAmp - float, maximum amplitude of image 
*****************************************************/       
void imageLocalFreqFlt(int nx, int ny, 
		       float *gCos, float *gSin, 
		       int subSample, 
		       float maskX[5][2], float maskY[5][2], float C[2][2],
		       float lambda, float thetaDeg, 
		       float epsAmp, float tau, 
		       int border, 
		       float *phiX, float *phiY, float *maxAmp );

/* Set cacheName to be the concatenation:
      pathName + prefix + root + postfix
   It is assumed that cachName has been declared to be sufficiently long
   to hold this concatenated string. */
void makeCacheName(char *cacheName, 
		   char *pathName, char *prefix, char *root, char *postfix);

/****  Version 4.0, May 1999. *******************************************

  Here is pixel coord x is to the left, y down.
         Theta = 0  -- i.e. theta == 0 corresponds to horizontal line
                 45 /
                 90 |
                 135 \

   Filter tuning:
     Peak Frequency:
       omega = 2*PI/lambda , lambda specified by user (typically >=3.5)
       theta  spatial orientation.
       
       dPhi[*] = (Expected phase difference)/PI in x,y directions
                 for steps of one (original) pixel.
       
       dPhi[0] = omega*sin(theta)/PI;
       dPhi[1] = omega*cos(theta)/PI;

Modified: May 2001 to use utilFree instead of free

Modifications from Version 3.*
 1) A simpler interface is provided in freemanPyramid.c and freemanPyramid.h 
    All you need to know is probably there.
 2) The recover* routines no longer open and read the original image files.
    Instead either the original image needs to be passed to these routines,
    or, if all the available information can be recovered from the cache,
    then a null image pointer can be provided.  Routines raise error if
    they cannot determine all the required information.
 3) Functions which compute the pyramid and/or filter responses
    also return an integer flag (localPyrBase) indicating whether the base
    pyramid image points to new storage (which can safely be freed), or
    existing storage allocated elsewhere by the user (which cannot be freed
    here).
 4) freeFastFilterOutputs() has been changed to use this flag.
 5) The global parameters can be set by invoking:
       resetFreemanPyramidParams(FreemanPyramidCntrlStruct * pc)
    See freemanPyramid.h and the function defaultFreemanPyramidCntrlStruct() in
    freemanPyramid.c
 6) The tuning constants lamToSig, rhoThetaScale, filterDiam have
    been changed.  See ./design/readme.m 

Contents: Version 4.0
=====================

Tuning params
--------------

resetFreemanPyramidParams() -set the global tunable parameters according
                           to values in the FreemanPyramidCntrlStruct argument.

Get filter responses given an image, using the cache:
--------------------
 recoverFastFilterOutputs()  - build pyramid, G2,H2 basis, filtered, amp/phase
                               and grad files for input image file. Uses cache.
                               Does multiple orientations and scales.

 recoverPyramid()  - build image pyramid for input image.  Uses cache.

Get filter responses given an image array (bytes or floats), no caching
of intermediate results:
-----------------------------------
 computeFastFilterOutputs()  - same as recover, but does not use cache.

 freeFastFilterOutputs() - frees images allocated by computeFastFilterOutputs

 computePyramid()  - build image pyramid for input image.  No cache.

 samplePyramid() - computes a particular level of the image pyramid.

Separate steps (not for public consumption):
---------------------------------------------
 getGradMap() - given freeman G2, H2 responses, computes gradient of
                    phase map for one orientation.

 getPhaseMap() - given complex filter responses, computes phase map

 getAmpMap() - given complex filter responses, computes amplitude map

 getFastFilterMaps()  - compute G2,H2 filter (one orientation) response
                        from precomputed basis maps.

 steerFilterMaps() - compute G2,H2 images given an input image of
                     steering directions

 ampFilterMaps() - given G2,H2 images, compute amplitude image,
                       amp = sqrt(G2^2+H2^2)

 getOrientEng()  - get oriented energy according to Freeman and Adelson paper.

 getFastBasisMaps()   - compute G2,H2 basis maps from image array.

 downSample()	  - downSample an image using a Gaussian filter.
 
 censor()	- mark low amplitude and phase singularities in image map.

 gaussianBlur() - blur an image using a Gaussian filter.

Using 2D convolution:
---------------------
 getFilterMaps()      - compute G2,H2 filter responses from image array
                        by 2D convolution. One orientation. 
                       (Slower version of previous two routines combined.)

Correlation routines:
---------------------
 corrByte()
 corrFlt()

 corrSepByte()
 corrSepFlt()

Misc. routines:
---------------
 thinAmpMap()
 int pyramidHeight()
 pyramidSubSample()

Local routines (primarily for internal use):
--------------------------------------------
 checkPhaseCache()
 checkPyramidCache()

 freemanConstants()
 applyFreeman()       

 build_sepFreeman_filter()
 build_freeman_filter()
 freeman()

 derivMask()
 imageLocalFreqFlt()
***************************************************************************/
#include        "mex.h"
#include        "macros.h"
#include        "utils.h"
#include        "imageFile-io.h"
#include        "phase.h"
#include        "freemanPyramid.h"
#include        "phaseUtil.h"

/* Constants for the definition of G2/H2 filters */
static float lamToSig = 4.1572; /* conversion factor for sigma(lambda) */
static float filterDiam=7.0; /* filter diameter is roughly 7.0*sigma*/

static char *pyrPrefix =   "_Pyr_";
static char *basisPrefix = "_Basis_";
static char *filterPrefix ="_Filt_";
static char *ampPrefix =   "_Amp_";
static char *phasePrefix = "_Phase_";
static char *gradPrefix =  "_Grad_";

/* Global tunable parameters.  These can be left as is. */
static int border = 0; /* border around image in which phase, amp and grad 
                          maps are not computed */
static float preSigma = 0.0;  /* First filter raw images by Gaussian
				 with std-dev = preSigma (0->no prefilter) */
static float rockBottom = ROCK_BOTTOM; /* Minimum amplitude of filter response
				          to consider non-zero */
static float tauThres=1.3; /* tau treshold for phase singularity constraint */
static float rhoThetaScale = 1.323; /* ratio of std devs for the theta
		and scale directions for the amplitude spectrum of G2+iH2 */
static float thetaTol = 30.0; /* maximum orientation difference for G2/H2
                                 wrt center filter tuning */

/* Reset the global tunable parameters according to values
   in the FreemanPyramidCntrlStruct.
   DO NOT use this unless you know what you are doing. */
void resetFreemanPyramidParams(FreemanPyramidCntrlStruct * pc)
{
  {
    /* Compute the subsampling rate for each scale, record in
       the pyramid control structure */
    int iScale;
    float sigma, rho;
    int dwnSample, subSampleDwn;
    for(iScale=0; iScale<pc->nScale; iScale++) {
      freemanConstants(pc->lambda[iScale], &dwnSample, &subSampleDwn, 
		       &sigma, &rho);
      pc->subSample[iScale] = dwnSample * subSampleDwn;
    }
  }

  /* Tunable constants */
  border = pc->border;             /* No border around filtered image */
  preSigma = pc->preSigma;         /* Don't prefilter input image */

  /* Don't change the following unless you know what you are doing */
  rockBottom = pc->rockBottom;
  tauThres = pc->tauThres; /* tau treshold for phase singularity constraint */
  filterDiam = pc->filterDiam; /* filter radius is roughly 7*sigma */
  thetaTol = pc->thetaTol; /* Orientation tolerance wrt center filter tuning */
}
                         
void recoverFastFilterOutputs(void *vimage, int imageBytesPerPixel,
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
	int *nxFilt, int *nyFilt)
 /*******************************************************
    General Idea:
    For each scale, this routine checks what information
    is desired, according to the storeMap array.  It 
    figures out what needs to be computed (computeMap[]).
    It checks the cache for what is needed, if the data isn't
    found it computes it.  New computations are cached
    when indicated by the global array cacheMap[]
 *******************************************************/
 {
  float *respPnt[NUMCODE][MAXLEN];
  float sigma, rho;
  float epsA, allOrientMaxAmp[NUMSCALE], maxAmp[NUMCHANNEL];
  float tmp;
  float outPar[MAXLEN];
 
  int pyrLevel, bytesPerPixel;
  int computeMap[NUMCODE][NUMSCALE];
  int nxPyrLocal[NUMSCALE], nyPyrLocal[NUMSCALE];
  int nxFiltLocal[NUMSCALE], nyFiltLocal[NUMSCALE];
  int numPar, dwnSmpl[NUMSCALE];
  int iChannelIn, iChannelOut;
  int iScale, iOrient, ix, iy, ixy, k, subSampleDwn;

  char filterFileTail[10], fname[MAXLEN];
 
  /* Check what is cached and what needs to be recomputed. */
  checkPhaseCache(computeMap, storeMap, lambda, nScale, thetaDeg, nOrient,
             allOrientMaxAmp,
             pathNameCache, rootNameCache, 
             respPnt, nxFiltLocal, nyFiltLocal);
  
  /* Decide what is needed from image pyramid */
  for(iScale=0; iScale<NUMSCALE; iScale++)  {
   computeMap[isPyr][iScale] = 0;
   nxPyrLocal[iScale] = nyPyrLocal[iScale] = 0;
  }
  
  *nPyr = max(0, *nPyr);
  for(iScale=0; iScale<nScale; iScale++) {
   if (storeMap[isPyr] || computeMap[isBasis][iScale] >= 2) {
    freemanConstants(lambda[iScale], &(dwnSmpl[iScale]),
               &subSampleDwn, &sigma, &rho);
    tmp = log((double)dwnSmpl[iScale])/lnTwo;
    pyrLevel = ROUND( tmp );
    computeMap[isPyr][pyrLevel] = 2;
    *nPyr = max(pyrLevel+1, *nPyr);
   }
  }
  if (*nPyr > NUMSCALE) {
    /* mexPrintf(" Requesting nPyr = %d : ", *nPyr); */
    error(" NUMSCALE (phase.h) too small","");
  }

  /* Recover what is needed in the image pyramid. */
  if (*nPyr > 0)  {
   recoverPyramid(vimage, imageBytesPerPixel, nxImage, nyImage,
                   pathNameCache, rootNameCache, *nPyr, localPyrBase,
                   respPnt[isPyr], nxPyrLocal, nyPyrLocal, cacheFiles[isPyr]);

   for(k=0; k< *nPyr; k++)
    computeMap[isPyr][k] = 1;
   /* Save pointers to pyramid, if desired */
   if (storeMap[isPyr] == 1)
    for(k=0; k< *nPyr; k++) {
     pyrMap[k] = respPnt[isPyr][k];
     nxPyr[k] = nxPyrLocal[k];
     nyPyr[k] = nyPyrLocal[k];
    }
  }

  for(iScale=0; iScale<nScale; iScale++) {
   /* Compute basis results */
   if (computeMap[isBasis][iScale] >= 2) {
    /* Find what down sampling rate is needed */
    freemanConstants(lambda[iScale], &(dwnSmpl[iScale]),
               &subSampleDwn, &sigma, &rho);
    tmp = log((double)dwnSmpl[iScale])/lnTwo;
    pyrLevel = ROUND( tmp );
    /*mexPrintf( " Computing basis maps for lambda = %4.1f, using pyramid level %d\n", lambda[iScale], pyrLevel);*/

    bytesPerPixel = 4;
    iChannelOut = index(0, iScale, NUMBASIS);

    getFastBasisMaps((void *)respPnt[isPyr][pyrLevel], bytesPerPixel, dwnSmpl[iScale],
         nxPyrLocal[pyrLevel], nxPyrLocal[pyrLevel], nyPyrLocal[pyrLevel], lambda[iScale],
         respPnt[isBasis]+iChannelOut, respPnt[isBasis]+iChannelOut+3,
         &(nxFiltLocal[iScale]), &(nyFiltLocal[iScale]));

    computeMap[isBasis][iScale] = 1;
    
    /*mexPrintf( " Basis map size %d %d\n", nxFiltLocal[iScale], nyFiltLocal[iScale]);*/

    /* Cache basis results */
    if (cacheFiles[isBasis] > 0) {
     sprintf(filterFileTail,".w%d",ROUND(lambda[iScale]));
     makeCacheName(fname, pathNameCache, basisPrefix, rootNameCache, filterFileTail);
     numPar = 1;
     outPar[0] = lambda[iScale];
     iChannelOut = index(0, iScale, NUMBASIS);
     writeFilterOutputs(fname, "freemanBasisFlt",
                         sizeof(float), 1, NUMBASIS, 
                         nxFiltLocal[iScale], nyFiltLocal[iScale],
                         numPar, outPar, respPnt[isBasis]+iChannelOut);
    }

   } /* End: Compute basis image */

   /* Save pointers to basis results */
   if (storeMap[isBasis] > 0) {
    iChannelOut = index(0, iScale, NUMBASIS);
    for(k=0; k<NUMBASIS; k++)
     basisMap[iChannelOut+k] = respPnt[isBasis][iChannelOut+k];
   } 
  } /* Over scales */

  /* Free pyramid stuff, if not to be returned to caller */
  if (*nPyr > 0 && storeMap[isPyr] == 0) {
   if ((*localPyrBase == 1) && (computeMap[isPyr][0]==1))
     utilFree((void **) &(respPnt[isPyr][0]));
   for(k=1; k< *nPyr; k++) {
    if (computeMap[isPyr][k] == 1)
     utilFree((void **) &(respPnt[isPyr][k]));
   }
  }

  /* Begin: Compute filter responses */
  for(iScale=0; iScale<nScale; iScale++) {
   if (computeMap[isFilter][iScale] >= 2){
 
     /*mexPrintf( " Computing filter maps for lambda = %4.1f:\n",
       lambda[iScale]);*/
      iChannelIn  = index(0, iScale, NUMBASIS);
      for(iOrient=0; iOrient<nOrient; iOrient++) {
       iChannelOut = 2 * index(iOrient, iScale, nOrient);
       getFastFilterMaps(respPnt[isFilter]+iChannelOut, thetaDeg[iOrient],
          respPnt[isBasis]+iChannelIn, respPnt[isBasis]+iChannelIn+3,
          nxFiltLocal[iScale], nyFiltLocal[iScale], &(maxAmp[iOrient]));
       /*mexPrintf( " %5.1f: %5.1f %d %d\n", thetaDeg[iOrient], 
	 maxAmp[iOrient], nxFiltLocal[iScale], nyFiltLocal[iScale]);*/
      }
      computeMap[isFilter][iScale] = 1;

      /* Cache filter results */
      if (cacheFiles[isFilter] >= 1) {
       sprintf(filterFileTail,".w%d",ROUND(lambda[iScale]));
       makeCacheName(fname, pathNameCache, filterPrefix, rootNameCache, 
                     filterFileTail);
       numPar = 3;
       for(k=0; k<nOrient; k++) {
        outPar[k*numPar] = lambda[iScale];
        outPar[k*numPar+1] = thetaDeg[k];
        outPar[k*numPar+2] = maxAmp[iOrient];
       }
       iChannelOut = 2 * index(0, iScale, nOrient);
       writeFilterOutputs(fname, "freemanFilterFlt",
                          sizeof(float), nOrient, 2, 
                          nxFiltLocal[iScale], nyFiltLocal[iScale],
                          numPar, outPar, respPnt[isFilter]+iChannelOut);
      } /* End cache filter results */

      /* Free basis results at this scale */
      if (storeMap[isBasis] == 0) {
       iChannelIn  = index(0, iScale, NUMBASIS);
       if (computeMap[isBasis][iScale] == 1) {
        for(k=0; k<NUMBASIS; k++) 
         utilFree((void **) &(respPnt[isBasis][k+iChannelIn]));
       }
      }

   } /* End: Compute filter results */
 
   /* Save pointers to filter results */
   if (storeMap[isFilter] > 0) {
    iChannelOut = index(0, iScale, 2*nOrient);
    for(k=0; k<nOrient*2;k++)
     filterMap[iChannelOut+k] = respPnt[isFilter][iChannelOut+k];
   }

   /* Compute amplitude results */
   if (computeMap[isAmp][iScale] >= 2)  { 

    allOrientMaxAmp[iScale] = 0;
    for(iOrient=0; iOrient<nOrient; iOrient++) {
     iChannelOut = index(iOrient, iScale, nOrient);
     iChannelIn  = 2 * index(iOrient, iScale, nOrient);
     getAmpMap(respPnt[isAmp]+iChannelOut, respPnt[isFilter]+iChannelIn,
          nxFiltLocal[iScale], nyFiltLocal[iScale], &(maxAmp[iOrient]));
     allOrientMaxAmp[iScale] = max(maxAmp[iOrient], allOrientMaxAmp[iScale]);
    }
    computeMap[isAmp][iScale] = 1;
 
    /* Cache amplitude results */
    if (cacheFiles[isAmp] >= 1) {
       sprintf(filterFileTail,".w%d",ROUND(lambda[iScale]));
       makeCacheName(fname, pathNameCache, ampPrefix, rootNameCache, filterFileTail);
       numPar = 3;
       for(k=0; k<nOrient; k++) {
        outPar[k*numPar] = lambda[iScale];
        outPar[k*numPar+1] = thetaDeg[k];
        outPar[k*numPar+2] = maxAmp[iOrient];
       }
       iChannelOut = index(0, iScale, nOrient);
       writeFilterOutputs(fname, "freemanAmpFlt",
                          sizeof(float), nOrient, 1, 
                          nxFiltLocal[iScale], nyFiltLocal[iScale],
                          numPar, outPar, respPnt[isAmp]+iChannelOut);
    } /* End cache amplitude results */
   } /* End: Compute amplitude results */

   /* Save pointers to amplitude results */
   if (storeMap[isAmp] > 0) {
    iChannelOut = index(0, iScale, nOrient);
    for(k=0; k<nOrient;k++)
     aMap[iChannelOut+k] = respPnt[isAmp][iChannelOut+k];
   }

   /* Compute phase results */
   if (computeMap[isPhase][iScale] >= 2)  { 

    for(iOrient=0; iOrient<nOrient; iOrient++) {
     iChannelOut = index(iOrient, iScale, nOrient);
     iChannelIn  = 2 * index(iOrient, iScale, nOrient);
     getPhaseMap(respPnt[isPhase]+iChannelOut, respPnt[isFilter]+iChannelIn,
          nxFiltLocal[iScale], nyFiltLocal[iScale]);
    }
    computeMap[isPhase][iScale] = 1;
 
    /* Cache phase results */
    if (cacheFiles[isPhase] >= 1) {
       sprintf(filterFileTail,".w%d",ROUND(lambda[iScale]));
       makeCacheName(fname, pathNameCache, phasePrefix, rootNameCache, filterFileTail);
       numPar = 2;
       for(k=0; k<nOrient; k++) {
        outPar[k*numPar] = lambda[iScale];
        outPar[k*numPar+1] = thetaDeg[k];
       }
       iChannelOut = index(0, iScale, nOrient);
       writeFilterOutputs(fname, "freemanPhaseFlt",
                          sizeof(float), nOrient, 1, 
                          nxFiltLocal[iScale], nyFiltLocal[iScale],
                          numPar, outPar, respPnt[isPhase]+iChannelOut);
    } /* End cache phase results */
   } /* End: Compute phase results */

   /* Save pointers to phase results */
   if (storeMap[isPhase] > 0) {
    iChannelOut = index(0, iScale, nOrient);
    for(k=0; k<nOrient;k++)
     pMap[iChannelOut+k] = respPnt[isPhase][iChannelOut+k];
   }

   /* Compute grad results */
   if (computeMap[isGrad][iScale] >= 2)  { 

    allOrientMaxAmp[iScale] = 0.0;
    /*mexPrintf( " Computing gradient maps for lambda = %4.1f:\n",
      lambda[iScale]);*/
    for(iOrient=0; iOrient<nOrient; iOrient++) {
     iChannelIn  = 2 * index(iOrient, iScale, nOrient);
     iChannelOut = iChannelIn;
     getGradMap(respPnt[isFilter]+iChannelIn,
                nxFiltLocal[iScale], nyFiltLocal[iScale],
                lambda[iScale], thetaDeg[iOrient],
                respPnt[isGrad]+iChannelOut, respPnt[isGrad]+iChannelOut+1,
                maxAmp+iOrient);
     allOrientMaxAmp[iScale] = max(maxAmp[iOrient], allOrientMaxAmp[iScale]);
    }
    computeMap[isGrad][iScale] = 1;
 
    /* Cache grad results */
    if (cacheFiles[isGrad] >= 1) {
       sprintf(filterFileTail,".w%d",ROUND(lambda[iScale]));
       makeCacheName(fname, pathNameCache, gradPrefix, rootNameCache, filterFileTail);
       numPar = 3;
       for(k=0; k<nOrient; k++) {
        outPar[k*numPar] = lambda[iScale];
        outPar[k*numPar+1] = thetaDeg[k];
        outPar[k*numPar+2] = maxAmp[k];
       }
       iChannelOut = 2*index(0, iScale, nOrient);
       writeFilterOutputs(fname, "freemanGradFlt",
                          sizeof(float), nOrient, 2, 
                          nxFiltLocal[iScale], nyFiltLocal[iScale],
                          numPar, outPar, respPnt[isGrad]+iChannelOut);
    } /* End cache grad results */
   } /* End: Compute grad results */

   /* Save pointers to grad results */
   if (storeMap[isGrad] > 0) {
    iChannelOut = 2 * index(0, iScale, nOrient);
    for(k=0; k<2*nOrient;k++)
     gMap[iChannelOut+k] = respPnt[isGrad][iChannelOut+k];
   }

   /* Now have either reconstructed and cached gradient data,
     or read it from a cache */
   if (storeMap[isGrad] == 1) {
    /* Threshold low amplitude responses */
    epsA = max( epsAmp, allOrientMaxAmp[iScale]*epsC);
    /* Throw out low amplitude responses */
    for(iOrient=0; iOrient<nOrient; iOrient++) {
      iChannelIn = index(iOrient, iScale, nOrient);
      for(iy=0; iy<nyFiltLocal[iScale]; iy++) {
       for(ix=0; ix<nxFiltLocal[iScale]; ix++) {
          ixy = index(ix, iy, nxFiltLocal[iScale]);
          if (respPnt[isAmp][iChannelIn][ixy] <= epsA) {
            /* phase singularities marked with grad_x,y = BIG_NEG */
            gMap[2*iChannelIn][ixy] = BIG_NEG;
            gMap[2*iChannelIn+1][ixy] = BIG_NEG;
          }
       }
      }
    }  /* End thresholding filtered images */
   }

   /* Free filter results */
   if (storeMap[isFilter] == 0 && computeMap[isFilter][iScale] == 1) {
    iChannelOut = index(0, iScale, 2*nOrient);
    for(k=0; k<nOrient*2;k++)
      utilFree((void **) &(respPnt[isFilter][iChannelOut+k]));
   }

   if (storeMap[isBasis] == 1 || storeMap[isFilter] == 1 || storeMap[isAmp] == 1 || storeMap[isPhase] == 1 || storeMap[isGrad] == 1) {
    nxFilt[iScale] = nxFiltLocal[iScale];
    nyFilt[iScale] = nyFiltLocal[iScale];
   }


  } /* over scales */
 }

void recoverPyramid(void *vimage, int imageBytesPerPixel,
		    int nxImage, int nyImage,
		    char *pathNameCache, char *rootNameCache, 
		    int nPyr, int *localPyrBase,
		    float **pyramid, int *nxPyr, int *nyPyr, 
		    int cachePyr)
 /*******************************************************
    General Idea:
    This routine computes an image pyramid with nPyr levels,
    subsampling the original image at 1, 2, 4, ... 2^(nPyr-1)
    pixels, respectively.
    For each level, this routine checks the cache. It then
    figures out what needs to be computed.
    New computations are cached when the flag cachePyr > 0
 *******************************************************/
 {
  float outPar[MAXLEN];
  unsigned char * image;
  int numPar, dwnSmpl, incrementalDwnSmpl, bytesPerPixel;
  int k, ixy, iPyr;
  int computePyr[NUMSCALE];
  char fname[MAXLEN], pyrFileTail[10];
  
  *localPyrBase = 0;
  for(k=0; k<nPyr; k++)  /* Compute the first nPyr levels */
   computePyr[k] = 2;

  /* Check what is cached and what needs to be recomputed. */
  checkPyramidCache( pyramid, localPyrBase, computePyr, nPyr, nxPyr, nyPyr, 
             pathNameCache, rootNameCache);
  dwnSmpl = 1;
  for(iPyr=0; iPyr<nPyr; iPyr++) {
    
   if (computePyr[iPyr] >= 2) { /* Compute pyramid image */
  
    if (iPyr == 0) { /* Convert input image to floats, if necessary */

     if (preSigma > 0.0) {
      gaussianBlur(vimage, imageBytesPerPixel, nxImage, nxImage, nyImage, 
              preSigma, (int) 1, &(pyramid[0]), nxPyr, nyPyr);
      *localPyrBase = 1;
     } else {
      if (imageBytesPerPixel == 4) {
       pyramid[0] = (float *) vimage;
       *localPyrBase = 0;  /* Reuse float image file as pyramid base level */
      } else if (imageBytesPerPixel == 1) {
       image = (unsigned char *) vimage;
       grabFloatMemory((float **) &(pyramid[0]), nxImage*nyImage,"pyramid[0]");
       *localPyrBase = 1;
       for(ixy=0; ixy<nxImage * nyImage; ixy++) 
        pyramid[0][ixy] = (float) image[ixy];
      } else {
        error("invalid number of bytesPerPixel for input image","");
      }
      nxPyr[0] = nxImage;
      nyPyr[0] = nyImage;
     }
     computePyr[0] = 1;
    } else { /* Downsample image */
     bytesPerPixel = 4;
     incrementalDwnSmpl = 2;
     downSample((void *)pyramid[iPyr-1], bytesPerPixel,
		nxPyr[iPyr-1], nxPyr[iPyr-1], nyPyr[iPyr-1],
		incrementalDwnSmpl, pyramid+iPyr, nxPyr+iPyr, nyPyr+iPyr);
     computePyr[iPyr] = 1;
    } /* Down sample > 1 */

    /* Cache the result */
    if (cachePyr > 0) {
      sprintf(pyrFileTail,".s%d",dwnSmpl);
      makeCacheName(fname, pathNameCache, pyrPrefix, rootNameCache, pyrFileTail);
      numPar = 1;
      outPar[0] = (float) dwnSmpl;
      writeFilterOutputs(fname, "pyramidFlt",
                         sizeof(float), (int) 1, (int) 1, 
                         nxPyr[iPyr], nyPyr[iPyr],
                         numPar, outPar, pyramid+iPyr);
    }

   } /* End: Compute pyramid image */
   dwnSmpl *= 2;
  } /* Over levels, iPyr */

 }

/* Return the Freeman pyramid maps specifed by the array of flags, storeMap,
   computed for the given image.  (If bytesPerPixel == 4, the input image
   is assumed to be a float image.)  Intermediate maps may be computed to
   form the ones that are returned, but these are freed by this routine. */
void computeFastFilterOutputs(void * vimage, int bytesPerPixel,
        int nxImage, int nyImage,
        float *lambda, int nScale,
        float *thetaDeg, int nOrient,
        int *storeMap,
        float **pyrMap, int *nPyr, int *localPyrBase,
        float **basisMap, float **filterMap,
        float **pMap, float **aMap, float **gMap,
        float epsAmp, float epsC,
        int *nxPyr, int *nyPyr,
        int *nxFilt, int *nyFilt)
 {
  float *respPnt[NUMCODE][MAXLEN];
  float sigma, rho, tmp;
  float epsA, allOrientMaxAmp[NUMSCALE], maxAmp[NUMCHANNEL];
 
  int computeMap[NUMCODE];
  int nxPyrLocal[NUMSCALE], nyPyrLocal[NUMSCALE];
  int nxFiltLocal[NUMSCALE], nyFiltLocal[NUMSCALE];
  int pyrLevel, dwnSmpl[NUMSCALE];
  int iChannelIn, iChannelOut;
  int iScale, iOrient, ix, iy, ixy, k, subSampleDwn;
  
  /* Set compute map to indicate the most dependent map, and
     all the maps required to generate this map */
  for(k=0; k<NUMCODE; k++) {
   computeMap[k] = 0;
   if (storeMap[k] > 0)
    computeMap[k] = 2;
  }
  if (computeMap[isGrad]>=2)
   computeMap[isAmp] = 2;
  if ( computeMap[isGrad] || computeMap[isPhase] || computeMap[isAmp])
   computeMap[isFilter] = 2;
  if (computeMap[isFilter]>=2)
   computeMap[isBasis] = 2;
  if (computeMap[isBasis]>=2)
   computeMap[isPyr] = 2;

  /* Initialize an empty pyramid */
  for(iScale=0; iScale<NUMSCALE; iScale++) {
   nxPyrLocal[iScale] = nyPyrLocal[iScale] = 0;
  }
  
  *nPyr = max(0, *nPyr);
  if (computeMap[isPyr] >= 2) {
   for(iScale=0; iScale<nScale; iScale++) {
    freemanConstants(lambda[iScale], &(dwnSmpl[iScale]),
               &subSampleDwn, &sigma, &rho);
    tmp = log((double)dwnSmpl[iScale])/lnTwo;
    pyrLevel = ROUND( tmp );
    *nPyr = max(pyrLevel+1, *nPyr);
   }
   if (*nPyr > NUMSCALE) {
     /*mexPrintf(" Requesting nPyr = %d : ", *nPyr);*/
    error(" NUMSCALE (phase.h) too small","");
   }

   /* Recover what is needed in the image pyramid. */
   computePyramid(vimage, bytesPerPixel, nxImage, nyImage, 
               *nPyr, localPyrBase, respPnt[isPyr], nxPyrLocal, nyPyrLocal);

   if (storeMap[isPyr] == 1)
    for(k=0; k<*nPyr; k++) {
     pyrMap[k] = respPnt[isPyr][k];
     nxPyr[k] = nxPyrLocal[k];
     nyPyr[k] = nyPyrLocal[k];
    }
  }

  if (computeMap[isBasis] >= 2) {
   for(iScale=0; iScale<nScale; iScale++) {
    /* Compute basis results */
    /* Find what down sampling rate is needed */
    freemanConstants(lambda[iScale], &(dwnSmpl[iScale]),
               &subSampleDwn, &sigma, &rho);
    tmp = log((double)dwnSmpl[iScale])/lnTwo;
    pyrLevel = ROUND( tmp );
    /*mexPrintf( " Computing basis maps for lambda = %4.1f, using pyramid level %d\n", lambda[iScale], pyrLevel);*/

    bytesPerPixel = 4;
    iChannelOut = index(0, iScale, NUMBASIS);

    getFastBasisMaps((void *)respPnt[isPyr][pyrLevel], bytesPerPixel,
	 dwnSmpl[iScale], nxPyrLocal[pyrLevel], 
	 nxPyrLocal[pyrLevel], nyPyrLocal[pyrLevel], lambda[iScale],
         respPnt[isBasis]+iChannelOut, respPnt[isBasis]+iChannelOut+3,
         &(nxFiltLocal[iScale]), &(nyFiltLocal[iScale]));

    /*mexPrintf( " Basis map size %d %d\n", nxFiltLocal[iScale], nyFiltLocal[iScale]);*/

    /* Save pointers to basis results */
    if (storeMap[isBasis] > 0) {
     iChannelOut = index(0, iScale, NUMBASIS);
     for(k=0; k<NUMBASIS; k++)
      basisMap[iChannelOut+k] = respPnt[isBasis][iChannelOut+k];
    } 
    
   } /* Over scales */

   /* Free pyramid stuff, if not to be returned to caller */
   if (storeMap[isPyr] == 0) {
    if (*localPyrBase == 1)
     utilFree((void **) &(respPnt[isPyr][0]));
    for(k=1; k<*nPyr; k++)
     utilFree((void **) &(respPnt[isPyr][k]));
   }
  }

  if (computeMap[isFilter] >= 2) {
   for(iScale=0; iScale<nScale; iScale++) {
 
    /* Begin: Compute filter responses */
     /*mexPrintf( " Computing filter maps for lambda = %4.1f:\n",
       lambda[iScale]);*/
    iChannelIn  = index(0, iScale, NUMBASIS);
    for(iOrient=0; iOrient<nOrient; iOrient++) {
     iChannelOut = 2 * index(iOrient, iScale, nOrient);
     getFastFilterMaps(respPnt[isFilter]+iChannelOut, thetaDeg[iOrient],
          respPnt[isBasis]+iChannelIn, respPnt[isBasis]+iChannelIn+3,
          nxFiltLocal[iScale], nyFiltLocal[iScale], &(maxAmp[iOrient]));
     /*mexPrintf( " %5.1f: %5.1f %d %d\n", thetaDeg[iOrient], maxAmp[iOrient], nxFiltLocal[iScale], nyFiltLocal[iScale]);*/
    }

    /* Free basis results at this scale */
    if (storeMap[isBasis] == 0) {
       iChannelIn  = index(0, iScale, NUMBASIS);
       for(k=0; k<NUMBASIS; k++) 
        utilFree((void **) &(respPnt[isBasis][k+iChannelIn]));
    }

    /* Save pointers to filter results */
    if (storeMap[isFilter] > 0) {
     iChannelOut = index(0, iScale, 2*nOrient);
     for(k=0; k<nOrient*2;k++)
      filterMap[iChannelOut+k] = respPnt[isFilter][iChannelOut+k];
    }

    /* Compute amplitude results */
    if (computeMap[isAmp] >= 2)  { 

     allOrientMaxAmp[iScale] = 0;
     for(iOrient=0; iOrient<nOrient; iOrient++) {
      iChannelOut = index(iOrient, iScale, nOrient);
      iChannelIn  = 2 * index(iOrient, iScale, nOrient);
      getAmpMap(respPnt[isAmp]+iChannelOut, respPnt[isFilter]+iChannelIn,
          nxFiltLocal[iScale], nyFiltLocal[iScale], &(maxAmp[iOrient]));
      allOrientMaxAmp[iScale] = max(maxAmp[iOrient], allOrientMaxAmp[iScale]);
     }
 
     /* Save pointers to amplitude results */
     if (storeMap[isAmp] > 0) {
      iChannelOut = index(0, iScale, nOrient);
      for(k=0; k<nOrient;k++)
       aMap[iChannelOut+k] = respPnt[isAmp][iChannelOut+k];
     }
    } /* End: Compute amplitude results */

    /* Compute phase results */
    if (computeMap[isPhase] >= 2)  { 

     for(iOrient=0; iOrient<nOrient; iOrient++) {
      iChannelOut = index(iOrient, iScale, nOrient);
      iChannelIn  = 2 * index(iOrient, iScale, nOrient);
      getPhaseMap(respPnt[isPhase]+iChannelOut, respPnt[isFilter]+iChannelIn,
          nxFiltLocal[iScale], nyFiltLocal[iScale]);
     }

     /* Save pointers to phase results */
     if (storeMap[isPhase] > 0) {
      iChannelOut = index(0, iScale, nOrient);
      for(k=0; k<nOrient;k++)
       pMap[iChannelOut+k] = respPnt[isPhase][iChannelOut+k];
     }

    } /* End: Compute phase results */

    /* Compute grad results */
    if (computeMap[isGrad] >= 2)  { 

     allOrientMaxAmp[iScale] = 0.0;
     /*mexPrintf( " Computing gradient maps for lambda = %4.1f:\n",
       lambda[iScale]);*/
     for(iOrient=0; iOrient<nOrient; iOrient++) {
      iChannelIn  = 2 * index(iOrient, iScale, nOrient);
      iChannelOut = iChannelIn;
      getGradMap(respPnt[isFilter]+iChannelIn,
                 nxFiltLocal[iScale], nyFiltLocal[iScale],
                 lambda[iScale], thetaDeg[iOrient],
                 respPnt[isGrad]+iChannelOut, respPnt[isGrad]+iChannelOut+1,
                 maxAmp+iOrient);
      allOrientMaxAmp[iScale] = max(maxAmp[iOrient], allOrientMaxAmp[iScale]);
     }

     /* Save pointers to grad results */
     if (storeMap[isGrad] > 0) {
      iChannelOut = 2 * index(0, iScale, nOrient);
      for(k=0; k<2*nOrient;k++)
       gMap[iChannelOut+k] = respPnt[isGrad][iChannelOut+k];
     }

     if (storeMap[isGrad] == 1) {
      /* Threshold low amplitude responses */
      epsA = max( epsAmp, allOrientMaxAmp[iScale]*epsC);
      /* Throw out low amplitude responses */
      for(iOrient=0; iOrient<nOrient; iOrient++) {
       iChannelIn = index(iOrient, iScale, nOrient);
       for(iy=0; iy<nyFiltLocal[iScale]; iy++) {
        for(ix=0; ix<nxFiltLocal[iScale]; ix++) {
         ixy = index(ix, iy, nxFiltLocal[iScale]);
         if (respPnt[isAmp][iChannelIn][ixy] <= epsA)
            /* phase singularities marked with grad_x = BIG_NEG */
            gMap[2*iChannelIn][ixy] = BIG_NEG;
        }
       }
      }  /* End thresholding filtered images */
     }
 
    } /* End: Compute grad results */

    /* free filter results, if not wanted */
    if (computeMap[isFilter] >= 1 && storeMap[isFilter] == 0) {
     iChannelOut = 2 * index(0, iScale, nOrient);
     for(k=0; k<nOrient*2;k++)      
       utilFree((void **) &(respPnt[isFilter][iChannelOut+k]));
    }

    if (storeMap[isBasis] == 1 || storeMap[isFilter] == 1 || storeMap[isAmp] == 1 || storeMap[isPhase] == 1 || storeMap[isGrad] == 1) {
     nxFilt[iScale] = nxFiltLocal[iScale];
     nyFilt[iScale] = nyFiltLocal[iScale];
    }

   } /* over scales */
  } /* End: Compute filter */
 }

/* Free allocated images for the current set of filter responses.
   Here storeMap[isPyr] !=0 indicates the pyramid images are to be freed,
   for scales 0, ..., nPyr-1.
   The lowest level in the pyramid (pyrMap[0]) is treated differently,
   since this may refer to a user's float image (indicated by
   localPyrBase == 0 (FALSE)).  In that case we take care NOT to
   free it here, rather the user needs to do that.
   Similarly the basis and various filter responses are to be freed over
   the range of scales specified by scaleRange.
   To free all scales, set scaleRange = (0, nScale).  Here scaleRange[1]
   is NOT deleted...just iScale = scaleRange[0]..(scaleRange[1]-1). */
void freeFastFilterOutputs(int scaleRange[2], int nOrient,
			   int *storeMap, 
			   float **pyrMap, int nPyr, int localPyrBase,
			   float **basisMap, float **filterMap, 
			   float **pMap, float **aMap, float **gMap)
 {
  int k, iScale, iShift;

  /* Free pyramid images */
  if (storeMap[isPyr]) {
    if (localPyrBase) /* Careful not to free user's image */
      utilFree((void **) &(pyrMap[0]));
    for(iScale=1; iScale<nPyr; iScale++)
      utilFree((void **) &(pyrMap[iScale]));
  } 

  /* Free G2/H2 images over range of scales */
  for(iScale=scaleRange[0]; iScale<scaleRange[1]; iScale++) {
   if (storeMap[isBasis]) {
    iShift = iScale*NUMBASIS;
    for(k=0; k<NUMBASIS; k++) 
     utilFree((void **) &(basisMap[k+iShift]));
   }
   if (storeMap[isFilter]) {
    iShift = iScale*nOrient*2;
    for(k=0; k<2*nOrient; k++) 
     utilFree((void **) &(filterMap[k+iShift]));
   }
   if (storeMap[isPhase]) {
    iShift = iScale*nOrient;
    for(k=0; k<nOrient; k++) 
     utilFree((void **) &(pMap[k+iShift]));
   }
   if (storeMap[isAmp]) {
    iShift = iScale*nOrient;
    for(k=0; k<nOrient; k++) 
     utilFree((void **) &(aMap[k+iShift]));
   }
   if (storeMap[isGrad]) {
    iShift = iScale*nOrient*2;
    for(k=0; k<2*nOrient; k++) 
     utilFree((void **) &(gMap[k+iShift]));
   }
  }
 }

/* computePyramid: same as recoverPyramid, only does not use cache. */
void computePyramid(void *vimage, int bytesPerPixel, 
	       int nxImage, int nyImage, 
               int nPyr, int *localPyrBase, 
	       float **pyramid, int *nxPyr, int *nyPyr)
 {
  unsigned char *image;
  int dwnSmpl, incrementalDwnSmpl;
  int ixy, iPyr;

  if (0<nPyr) {
   nxPyr[0] = nxImage;
   nyPyr[0] = nyImage;
   if (bytesPerPixel == 1 ) { /* Convert image to floats */
    image = (unsigned char *) vimage;
    grabFloatMemory((float **) &(pyramid[0]), nxImage * nyImage,
		    "pyramid[0]");
    *localPyrBase = TRUE;
    for(ixy=0; ixy<nxImage * nyImage; ixy++) {
      pyramid[0][ixy] = (float) image[ixy];
    }
   } else {  /* Set pyramid[0] to point to existing float image */
    pyramid[0] = (float *) vimage;
    *localPyrBase = FALSE;
   }
  }

  dwnSmpl = 1;
  for(iPyr=1; iPyr<nPyr; iPyr++) {
   dwnSmpl *= 2;
   bytesPerPixel = 4;
   incrementalDwnSmpl = 2;
   downSample((void *)pyramid[iPyr-1], bytesPerPixel, nxPyr[iPyr-1], 
              nxPyr[iPyr-1], nyPyr[iPyr-1],
              incrementalDwnSmpl, pyramid+iPyr, nxPyr+iPyr, nyPyr+iPyr);
  } /* Over levels, iPyr */
 }

/** Return just one level of the pyramid. Provides simpler
    interface to computePyramid.  (Not called within phaseUtil.c) **/
void samplePyramid(void *vimage, int bytesPerPixel, 
		   int nxImage, int nyImage, 
		   int pyrLevel, int *localPyrBase,
		   float **pDwnImageFlt, int *pnxDwn, int *pnyDwn)
 {
  float *scrPyramid[NUMSCALE];
  int l, nxPyrLocal[NUMSCALE], nyPyrLocal[NUMSCALE];

  if (pyrLevel < 0) 
   error("samplPyramid requested to expand an image","");

  if (pyrLevel >= NUMSCALE)
   error("samplPyramid: NUMSCALE too small for requested pryLevel","");

  computePyramid(vimage, bytesPerPixel, nxImage, nyImage, 
               pyrLevel+1, localPyrBase, scrPyramid, nxPyrLocal, nyPyrLocal);
  
  for(l=0; l< pyrLevel; l++)
   utilFree((void **) &(scrPyramid[l]));
  if (pyrLevel != 0)
   *localPyrBase = 0;
 
  *pDwnImageFlt = scrPyramid[pyrLevel];
  *pnxDwn = nxPyrLocal[pyrLevel];
  *pnyDwn = nyPyrLocal[pyrLevel];
  
 }

/* getGradMap:  compute phase gradient images, PhiX and PhiY, 
   given filter responses for a given scale and orientation. */
void getGradMap(float *filtResp[2], int nxFilt, int nyFilt, 
	   float lambda, float thetaDeg, 
           float **pPhiX, float **pPhiY, float *pMaxAmp)
 {
    float sigma, rho;
    float maskX[5][2], maskY[5][2], dPhi[2], C[2][2];
    int subSample, dwnSmpl, subSampleDwn, sizeFiltered;

    freemanConstants(lambda, &dwnSmpl, &subSampleDwn, &sigma, &rho);
    subSample = dwnSmpl * subSampleDwn;

    /* Build derivative mask for filter responses */
    derivMask(dPhi, C, maskX, maskY,
               lambda, thetaDeg, sigma, rho, subSample);
  
    /* Malloc space for phase, amp, gradient responses */
    sizeFiltered = nxFilt * nyFilt;
    grabFloatMemory(pPhiX, sizeFiltered, "grad_phase_x map");
    grabFloatMemory(pPhiY, sizeFiltered, "grad_phase_y map");

    /* Compute phi and grad phi maps from filter responses */
    imageLocalFreqFlt(nxFilt, nyFilt, filtResp[0], filtResp[1], subSample,
          maskX, maskY, C, lambda, thetaDeg, rockBottom, tauThres, border,
          *pPhiX, *pPhiY, pMaxAmp );
 }

/* getGradMap:  compute phase image, *pMap, 
   given filter responses for a given scale and orientation. */
void getPhaseMap(float **pMap, float **filtMap, 
		 int nx, int ny)
 {
  float respCos, respSin, amp, amp2, tmp;
  int ixy;

  grabFloatMemory( pMap, nx * ny, "phaseMap");

  for(ixy=0; ixy<nx*ny; ixy++) {
      respCos = filtMap[0][ixy];
      respSin = filtMap[1][ixy];
      amp2 = respCos*respCos + respSin*respSin;
      amp = sqrt(amp2) ;

      if (amp < rockBottom) {
       (*pMap)[ixy] = BIG_NEG;  /* Signifies no valid value */
      } else {
       tmp = atan2( respSin, respCos );
       /* Compute principal phase, scaled by PI */
       tmp = tmp/PI; 
       while(tmp<= -1.0) tmp += 2.0;
       while(tmp> 1.0) tmp -= 2.0;
       (*pMap)[ixy] = tmp;  /* Phase is divided by PI !!!!!! */
      } 
  }
 }

/* getAmpMap:  compute amplitude image, *aMap, 
   given filter responses for a given scale and orientation.
   Small amplitudes (below rockBottom) are censored (set
   to BIG_NEG. */
void getAmpMap(float **aMap, float **filtMap, int nx, int ny, float *pMax)
 {
  float respCos, respSin, amp, amp2;
  int ixy;

  grabFloatMemory( aMap, nx * ny, "ampMap");
  *pMax = 0.0;
  for(ixy=0; ixy<nx*ny; ixy++) {
      respCos = (filtMap[0])[ixy];
      respSin = (filtMap[1])[ixy];
      amp2 = respCos*respCos + respSin*respSin;
      amp = sqrt(amp2) ;
      *pMax = max( *pMax, amp);

      if (amp < rockBottom) {
       (*aMap)[ixy] = BIG_NEG;  /* Signifies no valid value */
      } else {
       (*aMap)[ixy] = amp;
      } 
  }
 }

/* getFastFilterMaps:  compute G2,H2 filter responses for a given scale
   and orientation given the seven basis images. */
void getFastFilterMaps(float **filtResp, float thetaDeg, 
		  float **respG2, float **respH2, 
		  int nxFilt, int nyFilt, float *pMaxAmp)
 {
    float theta, ctheta, stheta, G2k[3], H2k[4];
    float bound[2];
    int sizeFiltered, i,j,k, ixy;

    /* Malloc space for filter responses */
    sizeFiltered = nxFilt * nyFilt;
    grabFloatMemory(filtResp, sizeFiltered, "respR");
    grabFloatMemory(filtResp+1, sizeFiltered, "respI");

    /* Compute coefficients for basis maps */
    theta = PI * thetaDeg / 180.0;
    ctheta = cos(theta); stheta =sin(theta);
    G2k[0]=ctheta * ctheta;
    G2k[1]=2.0*ctheta*stheta;  /* In Freeman's paper he has - ve signs on */
    G2k[2]=stheta * stheta ;   /* some of these coefficients to make theta */
    H2k[0]=ctheta*ctheta*ctheta;/* +ve ccw - however, since Allan's code    */
    H2k[1]=3.0*ctheta*ctheta*stheta;/* reverses X and Y when generating the*/
    H2k[2]=3.0*ctheta*stheta*stheta;/* filter kernels, this becomes redundant*/
    H2k[3]=stheta * stheta * stheta ;

    /* Compute filter responses from basis maps */
    bound[0] = bound[1] = 0.0;
    for (j=0;j<nyFilt;j++)
     for (i=0; i<nxFilt; i++) {

      ixy = index(i,j,nxFilt);

      for (k=0, filtResp[0][ixy]=0.0; k<3; k++)
       filtResp[0][ixy] += G2k[k] * respG2[k][ixy];

      for (k=0, filtResp[1][ixy]=0.0; k<4; k++)
       filtResp[1][ixy] += H2k[k] * respH2[k][ixy];

      if (i>=border && i < nxFilt-border && j>= border && j< nyFilt-border){
       bound[0] = min(bound[0], filtResp[0][ixy]);
       bound[1] = max(bound[1], filtResp[0][ixy]);
       bound[0] = min(bound[0], filtResp[1][ixy]);
       bound[1] = max(bound[1], filtResp[1][ixy]);
      }

     }
  *pMaxAmp = abs(bound[0]);
  *pMaxAmp = max( (*pMaxAmp), abs(bound[1]));
 }

/* steerFilterMaps:  given a steering orientation image, thetaDeg[],
   compute G2,H2 filter responses for a given scale
   at the steering orientation for each corresponding pixel
   given the seven G2,H2 basis images. */
void steerFilterMaps(float **filtResp, float *thetaMap, 
		  float **respG2, float **respH2, 
		  int nxFilt, int nyFilt, float *pMaxAmp)
 {
    float theta, ctheta, stheta, G2k[3], H2k[4];
    float bound[2];
    int sizeFiltered, i,j,k, ixy;

    /* Malloc space for filter responses */
    sizeFiltered = nxFilt * nyFilt;
    grabFloatMemory(filtResp, sizeFiltered, "respR");
    grabFloatMemory(filtResp+1, sizeFiltered, "respI");

    /* Compute filter responses from basis maps */
    bound[0] = bound[1] = 0.0;
    for (j=0;j<nyFilt;j++)
     for (i=0; i<nxFilt; i++) {

      ixy = index(i,j,nxFilt);
      filtResp[0][ixy] = filtResp[1][ixy] = 0.0;

      if ((theta = thetaMap[ixy]) > -10.0) {
       /* Compute coefficients for basis maps */
       ctheta = cos(theta); stheta =sin(theta);
       G2k[0]=ctheta * ctheta;
       G2k[1]=2.0*ctheta*stheta;  /* In Freeman's paper he has - ve signs on */
       G2k[2]=stheta * stheta ;   /* some of these coefficients to make theta */
       H2k[0]=ctheta*ctheta*ctheta;/* +ve ccw - however, since Allan's code    */
       H2k[1]=3.0*ctheta*ctheta*stheta;/* reverses X and Y when generating the*/
       H2k[2]=3.0*ctheta*stheta*stheta;/* filter kernels, this becomes redundant*/
       H2k[3]=stheta * stheta * stheta ;

 
       for (k=0; k<3; k++)
        filtResp[0][ixy] += G2k[k] * respG2[k][ixy];

       for (k=0; k<4; k++)
        filtResp[1][ixy] += H2k[k] * respH2[k][ixy];

       if (i>=border && i < nxFilt-border && j>= border && j< nyFilt-border){
        bound[0] = min(bound[0], filtResp[0][ixy]);
        bound[1] = max(bound[1], filtResp[0][ixy]);
        bound[0] = min(bound[0], filtResp[1][ixy]);
        bound[1] = max(bound[1], filtResp[1][ixy]);
       }
      }
     }
  *pMaxAmp = abs(bound[0]);
  *pMaxAmp = max( (*pMaxAmp), abs(bound[1]));
 }

/* ampFilterMaps:  compute amplitude image, *pAmp, 
   given filter response images filtResp[2] for a given scale and
   orientation.  Small amplitudes are NOT censored. */
void ampFilterMaps(float **pAmp, float *filtResp[2],
		   int nxFilt, int nyFilt, float *pMaxAmp)
 {
    float amp, bound[2];
    int sizeFiltered, i,j,ixy;

    /* Malloc space for filter responses */
    sizeFiltered = nxFilt * nyFilt;
    grabFloatMemory(pAmp, sizeFiltered, "engMap");

    /* Compute energy from filter maps */
    bound[0] = bound[1] = 0.0;
    for (j=0;j<nyFilt;j++)
     for (i=0; i<nxFilt; i++) {

      ixy = index(i,j,nxFilt);
      amp = filtResp[0][ixy] * filtResp[0][ixy];
      amp += filtResp[1][ixy] * filtResp[1][ixy];
      (*pAmp)[ixy] = sqrt(amp);

      if (i>=border && i < nxFilt-border && j>= border && j< nyFilt-border){
        bound[0] = min(bound[0], amp);
        bound[1] = max(bound[1], amp);
      }
     }
  *pMaxAmp = abs(bound[0]);
  *pMaxAmp = max( (*pMaxAmp), abs(bound[1]));
 }

/* Get orientation energy coefficience C1, C2 and C3 (see p.905 Freeman and
   Adelson, Table XI, with sign of theta reversed) */
void getOrientEng(float **basisG2, float **basisH2, int nx, int ny, 
		  float epsContrast, float epsAmp, 
		  float **pUnifMap, float **pAmpMap, float **pThetaMap)
 {
  float amp, C1, C2, C3;
  int ix, iy, ixy;

  grabFloatMemory(pUnifMap, nx*ny, "orientUnifMap");
  grabFloatMemory(pAmpMap, nx*ny, "orientAmpMap");
  grabFloatMemory(pThetaMap, nx*ny, "orientThetaMap");
  
  /* From p.905 Freeman and Adelson, Table XI, with sign of theta reversed */
  
  for(iy=0; iy<ny; iy++)
   for(ix=0; ix<nx; ix++) {
    ixy = index(ix,iy,nx);
    C1 = 0.5 * basisG2[1][ixy] * basisG2[1][ixy];
    C1 += 0.25 * basisG2[0][ixy] * basisG2[2][ixy];
    C1 += 0.375 * (basisG2[0][ixy] * basisG2[0][ixy] + basisG2[2][ixy] * basisG2[2][ixy]);
    C1 += 0.3125 * (basisH2[0][ixy] * basisH2[0][ixy] + basisH2[3][ixy] * basisH2[3][ixy]);
    C1 += 0.5625 * (basisH2[1][ixy] * basisH2[1][ixy] + basisH2[2][ixy] * basisH2[2][ixy]);
    C1 += 0.375 * (basisH2[0][ixy] * basisH2[2][ixy] + basisH2[1][ixy] * basisH2[3][ixy]);

    C2 = 0.5 * ( basisG2[0][ixy] * basisG2[0][ixy] - basisG2[2][ixy] * basisG2[2][ixy]);
    C2 += 0.46875 * (basisH2[0][ixy] * basisH2[0][ixy] - basisH2[3][ixy] * basisH2[3][ixy]);
    C2 += 0.28125 * (basisH2[1][ixy] * basisH2[1][ixy] - basisH2[2][ixy] * basisH2[2][ixy]);
    C2 += 0.1875 * (basisH2[0][ixy] * basisH2[2][ixy] - basisH2[1][ixy] * basisH2[3][ixy]);
     
    C3 = basisG2[0][ixy] * basisG2[1][ixy] + basisG2[1][ixy] * basisG2[2][ixy];
    C3 += 0.9375 * ( basisH2[2][ixy] * basisH2[3][ixy] + basisH2[0][ixy] * basisH2[1][ixy]);
    C3 += 1.6875 * basisH2[1][ixy] * basisH2[2][ixy];
    C3 += 0.1875 * basisH2[0][ixy] * basisH2[3][ixy];

    (*pUnifMap)[ixy] = C1;
    amp =sqrt(C2 * C2 + C3 * C3);
    (*pAmpMap)[ixy] = amp;

    if ((*pAmpMap)[ixy] > epsAmp && C1 > 0.0 && amp/C1 > epsContrast) 
     (*pThetaMap)[ixy] = atan2(C3, C2)/2.0;
    else
     (*pThetaMap)[ixy] = -10.0; /* large negative response when zero amp */
   }

 }

/* Get G2/H2 basis maps using 1D convolution. */
void getFastBasisMaps(void *vimage, int bytesPerPixel,
		      int dwnSmplImage,
		      int nxImage, int nx0, int ny0, float lambda, 
		      float **respG2, float **respH2,
		      int *pnxFilt, int *pnyFilt)
 {
  float *G2x[3], *G2y[3], *H2x[4], *H2y[4], *temp;
  float lambdaDwn, sigma, sigmaDwn, rho;
  int subSample, dwnSmpl, nx, ny;
  int filtSize, sizeResp;
  int i, nxTemp;
  
  freemanConstants(lambda, &dwnSmpl, &subSample, &sigma, &rho);
  if (dwnSmpl != dwnSmplImage)
   error("Wrong downsampling rate in basis computation","");
  lambdaDwn = lambda/dwnSmpl;
  sigmaDwn = sigma/dwnSmpl;

  *pnxFilt = nx = nx0/subSample;
  *pnyFilt = ny = ny0/subSample;

  build_sepFreeman_filter(G2x, G2y, H2x, H2y, &filtSize, lambdaDwn, sigmaDwn, rho);

  /**********************************************************/
  /*** allocate memory for basis responses     **************/
  /**********************************************************/

  sizeResp = nx * ny * sizeof(float);
  for(i=0; i<3; i++)
   grabFloatMemory(respG2+i, sizeResp, "respG2");
  for(i=0; i<4; i++)
   grabFloatMemory(respH2+i, sizeResp, "respH2");
  nxTemp = nx0;
  grabFloatMemory(&temp, nxTemp*ny*sizeof(float), "temp");

  /**********************************************************/
  /*** Perform Basic Filter Correlations ********************/
  /**********************************************************/
  /*mexPrintf("Performing G2 correlations ...");*/
  if (bytesPerPixel==1) {
   for (i=0; i<3; i++) {
    corrSepByte((unsigned char *)vimage, (int) 1, nxImage, nx0, ny0, subSample,
         G2x[i], G2y[i], filtSize, 
         temp, respG2[i], nx);
   }
  } else {
   for (i=0; i<3; i++) {
    corrSepFlt((float *)vimage, (int) 1, nxImage, nx0, ny0, subSample,
         G2x[i], G2y[i], filtSize, 
         temp, respG2[i], nx);
   }
  }
  /*mexPrintf(" done\nPerforming H2 correlations ...");*/
  if (bytesPerPixel==1) {
   for (i=0; i<4; i++) {
    corrSepByte((unsigned char *)vimage, (int) 1, nx0, nx0, ny0, subSample,
         H2x[i], H2y[i], filtSize, 
         temp, respH2[i], nx);
   }
  } else {
   for (i=0; i<4; i++) {
    corrSepFlt((float *)vimage, (int) 1, nx0, nx0, ny0, subSample,
         H2x[i], H2y[i], filtSize, 
         temp, respH2[i], nx);
   }
  }
  /*mexPrintf(" done\n");*/

  utilFree((void **) &(temp));
  for(i=0; i<3; i++) {
   utilFree((void **) &(G2x[i]));
   utilFree((void **) &(G2y[i]));
  }
  for(i=0; i<4; i++) {
   utilFree((void **) &(H2x[i]));
   utilFree((void **) &(H2y[i]));
  }
 }

/* Downsample an image using a Gaussian blur prefilter */
void downSample(void *vimage, int bytesPerPixel, 
		int nxImage, int nx0, int ny0, int dwnSmpl,
                float **resp, int *pnxFilt, int *pnyFilt)
 {
  float *G_filt, *temp;
  float sigma, sum;
  int nx, ny;
  int filtSize, filtSize2, sizeResp;
  int i;
  
  if (dwnSmpl < 1) {
   mexPrintf("\nWarning: dwnSmpl = %d, ignored\n", dwnSmpl);
   return;
  }

  /* Build Gaussian blur for dwnSmpl */
  sigma = dwnSmpl/2.208;  /* Chosen to make contributions at sample
                             point, and half way between sample points equal */

  filtSize = 2 * ROUND(filterDiam/2.0 * sigma) + 1; 
  filtSize2 = filtSize/2;
  grabFloatMemory( &G_filt, filtSize, "Gaussian filter");
  sum = G_filt[filtSize2] = 1.0;
  for(i=1; i<=filtSize2; i++) {
   G_filt[filtSize2 -i] = G_filt[filtSize2 +i] = exp( -0.5 * (i*i)/sigma/sigma);
   sum += 2.0 * G_filt[filtSize2 +i];
  }
  for(i=0; i<filtSize; i++) 
   G_filt[i] /= sum;

  *pnxFilt = nx = nx0/dwnSmpl;
  *pnyFilt = ny = ny0/dwnSmpl;

  /**********************************************************/
  /*** allocate memory for basis responses     **************/
  /**********************************************************/

  sizeResp = nx * ny * sizeof(float);
  grabFloatMemory(&temp, nx0*ny*sizeof(float), "temp");
  grabFloatMemory( resp, sizeResp, "respGaussianBlur");

  /**********************************************************/
  /*** Perform Basic Filter Correlations ********************/
  /**********************************************************/
  /*mexPrintf("Performing Gaussian blur  correlations ...");*/
  if (bytesPerPixel==1) {
    corrSepByte((unsigned char *)vimage, 1, nxImage, nx0, ny0, dwnSmpl,
         G_filt, G_filt, filtSize, temp, *resp, nx);
  } else {
    corrSepFlt((float *)vimage, 1, nxImage, nx0, ny0, dwnSmpl,
         G_filt, G_filt, filtSize, temp, *resp, nx);
  }
  /*mexPrintf(" done\n");*/

  utilFree((void **) &(temp));
  utilFree((void **) &(G_filt));
 }

/* censor an image by doing an OR of a given image with a previously
   censored mask image.  That is, given the mask function maskFlt, and
   the image imageFlt, set any pixel ixy of imageFlt to BIG_NEG for which
   maskFlt[ixy] <= BIG_NEG. */
void censor(float * imageFlt, float *maskFlt, int nx, int ny)
 {
  int ixy;
  for(ixy=0; ixy < nx * ny; ixy++)
   if (maskFlt[ixy] <= BIG_NEG)
    imageFlt[ixy] = BIG_NEG;
 }

/* Perform a gaussian blur on an input image. */
void gaussianBlur(void *vimage, int bytesPerPixel, 
		  int nxImage, int nx0, int ny0,
		  float sigma, int dwnSmpl, 
                  float **resp, int *pnxFilt, int *pnyFilt)
 {
  float *G_filt, *temp;
  float sum;
  int nx, ny;
  int filtSize, filtSize2, sizeResp;
  int i;
  
  if (dwnSmpl < 1) {
   mexPrintf("\nWarning: dwnSmpl = %d, ignored\n", dwnSmpl);
   return;
  }

  /* Build Gaussian blur for dwnSmpl */
  filtSize = 2 * ROUND(2.5 * sigma) + 1; 
  filtSize2 = filtSize/2;
  grabFloatMemory( &G_filt, filtSize, "Gaussian filter");
  sum = G_filt[filtSize2] = 1.0;
  for(i=1; i<=filtSize2; i++) {
   G_filt[filtSize2 -i] = G_filt[filtSize2 +i] = exp( -0.5 * (i*i)/sigma/sigma);
   sum += 2.0 * G_filt[filtSize2 +i];
  }
  for(i=0; i<filtSize; i++) 
   G_filt[i] /= sum;

  *pnxFilt = nx = nx0/dwnSmpl;
  *pnyFilt = ny = ny0/dwnSmpl;

  /****************************************************/
  /*** allocate memory for responses     **************/
  /****************************************************/

  sizeResp = nx * ny * sizeof(float);
  grabFloatMemory(&temp, nx0*ny*sizeof(float), "temp");
  grabFloatMemory( resp, sizeResp, "respGaussianBlur");

  /**********************************************************/
  /*** Perform Basic Filter Correlations ********************/
  /**********************************************************/
  /*mexPrintf("Performing Gaussian blur  correlations ...");*/
  if (bytesPerPixel==1) {
    corrSepByte((unsigned char *)vimage, 1, nxImage, nx0, ny0, dwnSmpl,
         G_filt, G_filt, filtSize, temp, *resp, nx);
  } else {
    corrSepFlt((float *)vimage, 1, nxImage, nx0, ny0, dwnSmpl,
         G_filt, G_filt, filtSize, temp, *resp, nx);
  }
  /*mexPrintf(" done\n");*/

  utilFree((void **) &(temp));
  utilFree((void **) &(G_filt));
 }

/* Get G2/H2 filter maps for a given scale and orientation using
   2D correlation.  This is slower than getFastFilterMaps, and should
   only be used for a sanity check. */
void getFilterMaps(void *vimage, int bytesPerPixel, 
		   int nx0, int ny0, 
		   float lambda, float thetaDeg,
                   int *pnxFiltered, int *pnyFiltered, int *pSubSample,
                   float **filtResp, float *pMaxAmp)
 {
    float sigma, rho;
    int dwnSmpl, subSampleDwn, sizeFiltered;

    freemanConstants(lambda, &dwnSmpl, &subSampleDwn, &sigma, &rho);
    *pSubSample = dwnSmpl * subSampleDwn;

    *pnxFiltered = nx0/(*pSubSample);
    *pnyFiltered = ny0/(*pSubSample);

    /* Malloc space for correlation response */
    sizeFiltered = (*pnxFiltered) * (*pnyFiltered);
    grabFloatMemory(filtResp, sizeFiltered, "respR");
    grabFloatMemory(filtResp+1, sizeFiltered, "respI");

    /* Compute convolution with complex-valued filter */
    applyFreeman(vimage, bytesPerPixel, nx0, nx0, ny0, lambda, thetaDeg,
                     filtResp[0], filtResp[1], *pnxFiltered, pMaxAmp) ;
 }

/* Correlate input image with a 2D filter mask, and downsample the result */
void corrByte(unsigned char *input, int dwnSmpl, int nxIn,
	      int nx, int ny, int sample,
	      float *filter, int filter_size,
	      float *output, int nxOut,
	      float *bound)
  /* Input image size nx, ny,  Storage like: input[* by nxIn];
   * Downsample image at rate dwnSmpl
   * Correlate downsampled image, with filter
   * Sample result every sample pixels in DOWNSAMPLED image;
   * Store result in unsigned char image output[* by nxOut]
   * return upper and lower bounds of result
   *********************************************/
  {
    float val;
    int ds;
    int i,ii, ir, j, jj, jr, i0, j0;

    ds = dwnSmpl * sample;
    bound[0] = bound[1] = 0;
    for(i=0;i<ny/ds;i++){
     for(j=0;j<nx/ds;j++){
        val = 0.0;
        i0 = i*ds;
        j0 = j*ds;
        for(ii= 0;ii<filter_size;ii++){
         /* With wrap around 
         ir = (i0 + (ii - filter_size/2)*dwnSmpl + ny) % ny;
         *****/
         /* With nearest pixel */
         ir = i0 + (ii - filter_size/2)*dwnSmpl;
         ir = max(ir, 0);
         ir = min(ir, ny-1);

         for(jj= 0;jj<filter_size;jj++){
           /* With wrap around 
           jr = (j0 + (jj - filter_size/2)*dwnSmpl + nx) % nx;
           */
           /* With nearest pixel */
           jr = j0 + (jj - filter_size/2)*dwnSmpl;
           jr = max(jr, 0);
           jr = min(jr, nx-1);

           val += filter[index(jj, ii, filter_size)] *
                  ((int) input[index(jr,ir,nxIn)]);
         }
        }
        output[index(j,i,nxOut)] = val;
        if ( i>=border && j>=border && i<ny/ds-border && j<nx/ds-border) {
         bound[0] = min( bound[0] , val );
         bound[1] = max( bound[1] , val );
        }         
     }
    }
  }

/* Correlate a float image with a 2D filter mask, and downsample the result */
void corrFlt(float *input, int dwnSmpl, int nxIn,
	      int nx, int ny, int sample,
	      float *filter, int filter_size,
	      float *output, int nxOut,
	      float *bound)
  /* Input image size nx, ny,  Storage like: input[* by nxIn];
   * Downsample image at rate dwnSmpl
   * Correlate downsampled image, with filter
   * Sample result every sample pixels in DOWNSAMPLED image;
   * Store result in unsigned char image output[* by nxOut]
   * return upper and lower bounds of result
   *********************************************/
  {
    float val;
    int ds;
    int i,ii, ir, j, jj, jr, i0, j0;

    ds = dwnSmpl * sample;
    for(i=0;i<ny/ds;i++){
     for(j=0;j<nx/ds;j++){
        val = 0.0;
        i0 = i*ds;
        j0 = j*ds;
        for(ii= 0;ii<filter_size;ii++){
         /* With wrap around 
         ir = (i0 + (ii - filter_size/2)*dwnSmpl + ny) % ny;
         *****/
         /* With nearest pixel */
         ir = i0 + (ii - filter_size/2)*dwnSmpl;
         ir = max(ir, 0);
         ir = min(ir, ny-1);

         for(jj= 0;jj<filter_size;jj++){
           /* With wrap around 
           jr = (j0 + (jj - filter_size/2)*dwnSmpl + nx) % nx;
           */
           /* With nearest pixel */
           jr = j0 + (jj - filter_size/2)*dwnSmpl;
           jr = max(jr, 0);
           jr = min(jr, nx-1);

           val += filter[index(jj, ii, filter_size)] *
                   input[index(jr,ir,nxIn)];
         }
        }
        output[index(j,i,nxOut)] = val;
        if ( i==0 && j==0) 
         bound[0] = bound[1] = val;
        else {
         bound[0] = min( bound[0] , val );
         bound[1] = max( bound[1] , val );
        }         
     }
    }
  }

/* Separable correlations on an image, and downsample the result */
void corrSepByte(unsigned char *input, int dwnSmpl,int nxIn,
	    int nx0, int ny0, int sample,
	    float *filter_X, float *filter_Y, int filter_size,
	    float *temp, float *output,int  nxOut)
  {
    float val;
    int i,j, i0,j0, ii, jj, ir, jr, fs2;
    int ds, nxTemp;

    fs2 = filter_size / 2 ;
    ds = dwnSmpl * sample;

    /* Filter in Y direction */
    /* Result: temp is nx0/dwnSmpl by ny0/(dwnSmpl * sample) */
    nxTemp = nx0/dwnSmpl;
    for (i=0; i<ny0/ds; i++) {
     i0 = i*ds;
     for (j=0; j<nx0/dwnSmpl; j++){
      val = 0.0 ;
      j0 = j *dwnSmpl;
      for (ii=0; ii<filter_size; ii++) {
       ir = i0 + (ii-fs2)*dwnSmpl;
       /* Use nearest pixel within image */
       ir = max(ir, 0);
       ir = min(ir, ny0-1);
       val += filter_Y[ii] * (int) input[index(j0,ir,nxIn)];
      }
      temp[index(j,i,nxTemp)] = val ;
     }
    }

    /* Filter in X direction */
    for (j=0; j<nx0/ds; j++){
     j0 = j * sample;
     for (i=0; i<ny0/ds; i++) { 
      val = 0.0 ;
      for (jj=0; jj<filter_size; jj++) {
       jr = j0 + (jj-fs2);
       /* Use nearest pixel within image */
       jr = max(jr, 0);
       jr = min(jr, nxTemp-1);
       val += filter_X[jj] * temp[index(jr,i,nxTemp)] ;
      }
      output[index(j,i,nxOut)] = val ;
     }
    }

  } /* End: corrSepByte */

/* Separable correlations on a float image, and downsample the result */
void corrSepFlt(float *input, int dwnSmpl,int nxIn,
	    int nx0, int ny0, int sample,
	    float *filter_X, float *filter_Y, int filter_size,
	    float *temp, float *output,int  nxOut)
  {
    float val;
    int i,j, i0,j0, ii, jj, ir, jr, fs2;
    int ds, nxTemp;

    fs2 = filter_size / 2 ;
    ds = dwnSmpl * sample;

    /* Filter in Y direction */
    /* Result: temp is nx0/dwnSmpl by ny0/(dwnSmpl * sample) */
    nxTemp = nx0/dwnSmpl;
    for (i=0; i<ny0/ds; i++) {
     i0 = i*ds;
     for (j=0; j<nx0/dwnSmpl; j++){
      val = 0.0 ;
      j0 = j *dwnSmpl;
      for (ii=0; ii<filter_size; ii++) {
       ir = i0 + (ii-fs2)*dwnSmpl;
       /* Use nearest pixel within image */
       ir = max(ir, 0);
       ir = min(ir, ny0-1);
       val += filter_Y[ii] * input[index(j0,ir,nxIn)];
      }
      temp[index(j,i,nxTemp)] = val ;
     }
    }

    /* Filter in X direction */
    for (j=0; j<nx0/ds; j++){
     j0 = j * sample;
     for (i=0; i<ny0/ds; i++) { 
      val = 0.0 ;
      for (jj=0; jj<filter_size; jj++) {
       jr = j0 + (jj-fs2);
       /* Use nearest pixel within image */
       jr = max(jr, 0);
       jr = min(jr, nxTemp-1);
       val += filter_X[jj] * temp[index(jr,i,nxTemp)] ;
      }
      output[index(j,i,nxOut)] = val ;
     }
    }

  } /* End: corrSepFlt */

/* Do non-max suppression on amp Map */
void thinAmpMap(unsigned char **indMap, float radius, 
		float **aMap, float *thetaDeg, int nOrient, 
		int nx, int ny)
 {
  float val, theta, pVec[2], stepLength, nVal;
  int ix, iy, ixy, k, j;

  grabByteMemory((char **)indMap, nx*ny, "indicator map");
  for(iy=0; iy<ny; iy++)
   for(ix=0; ix<nx; ix++) {
    ixy = index(ix,iy, nx);

    /* Find local max amplitude */
    k = 0;
    val = aMap[0][ixy];
    for(j=1; j<nOrient; j++)
     if (aMap[j][ixy] > val) {
      val = aMap[j][ixy];
      k = j;
     }
    if (val <= BIG_NEG)
     (*indMap)[ixy] = outlierIndex;
    else {
     (*indMap)[ixy] = (unsigned char) k;

     if (radius >= 1.0) {
      /* Look along perp direction, up to a distance `radius' */
      theta = thetaDeg[k] * PI/180.0;
      pVec[0] = sin(theta);
      pVec[1] = cos(theta);
      for(stepLength=1.0; stepLength <= radius+0.5; stepLength+=1.0){
       nVal = interpMap(aMap[k], ix+stepLength*pVec[0], iy+stepLength*pVec[1], nx, ny);
       if (nVal > val) {
        (*indMap)[ixy] = outlierIndex;
        break;
       }
       nVal = interpMap(aMap[k], ix-stepLength*pVec[0], iy-stepLength*pVec[1], nx, ny);
       if (nVal > val) {
        (*indMap)[ixy] = outlierIndex;
        break;
       }
      }
     }
    }
   }
 }

/* Compute the maximum height of the pyramid needed for
   a given set of wavelengths */
int pyramidHeight(float *lambda, int nScale)
 {
  float sigma, rho, tmp;
  int nPyr, iScale, pyrLevel, downSample, subSample;
  
  nPyr = 0;
  for(iScale=0; iScale<nScale; iScale++) {
   freemanConstants(lambda[iScale], &downSample,
               &subSample, &sigma, &rho);
   tmp = log((double)downSample)/lnTwo;
   pyrLevel = ROUND( tmp );
   nPyr = max(pyrLevel+1, nPyr);
  }
  
  return(nPyr);
 }

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
                     float *lambda, int nScale)
 {
  float sigma, rho, tmp;
  int iScale, pyrLevel;
  
  *pnPyr = max(*pnPyr, 0);
  for(iScale=0; iScale<nScale; iScale++) {
   freemanConstants(lambda[iScale], downSample + iScale,
               subSample+iScale, &sigma, &rho);
   tmp = log((double)downSample[iScale])/lnTwo;
   pyrLevel = ROUND( tmp );
   *pnPyr = max(pyrLevel+1, *pnPyr);
  }
  if (*pnPyr > NUMSCALE) {
    mexPrintf(" Requesting nPyr = %d : ", *pnPyr);
    error(" NUMSCALE (phase.h) too small","");
  }
  lowPassSample[0] = 1;
  for(pyrLevel = 1; pyrLevel < *pnPyr; pyrLevel++) {
    lowPassSample[pyrLevel] = 2 * lowPassSample[pyrLevel-1];
  }
 }

/* Use bilinear interpolation to interpolate image map at x,y */
float interpMap(float *fMap, float x, float y, int nx, int ny)
 {
  float rx, ry, val0, val1;
  int floor_x, ceil_x;
  int floor_y, ceil_y;

  floor_x = (int) x;
  ceil_x = ceil(x);
  floor_y = (int) y;
  ceil_y = ceil(y);
  if (floor_x>=0 && ceil_x<nx && floor_y>=0 && ceil_y<ny) {
   rx = x - floor_x;
   ry = y - floor_y;
   val0 = linearInterp( fMap[index(floor_x, floor_y, nx)], fMap[index(ceil_x, floor_y, nx)], rx);
   val1 = linearInterp( fMap[index(floor_x, ceil_y, nx)], fMap[index(ceil_x, ceil_y, nx)], rx);
   return((float)linearInterp(val0, val1, ry));
  } else {
   return(BIG_NEG);
  }
 }

/* Check the cache for existing freeman pyramid images */
void checkPhaseCache(int computeMap[NUMCODE][NUMSCALE], int *storeMap,
		     float *lambda, int nScale, 
		     float *thetaDeg, int nOrient,
		     float *allOrientMaxAmp,
		     char *pathNameCache, char *rootNameCache, 
		     float *respPnt[NUMCODE][MAXLEN],
		     int *nxFilt, int *nyFilt)
 {
  float maxAmp[NUMORIENT];
  float pedestal, scl, lambdaRead, thetaDegRead;
  
  int iScale, iChannel, iOrient, checkParam, nxTmp, nyTmp;
  int k;
  char fname[MAXLEN], scrName[MAXLEN], header[NUMORIENT][HEADER];
  char filterFileTail[10];
  FILE *fpIn;

  /* Look over scales, check what is stored in cache */
  for(iScale=0; iScale<nScale; iScale++) {
   sprintf(filterFileTail,".w%d",ROUND(lambda[iScale]));
   allOrientMaxAmp[iScale] = -1.0;
  
   computeMap[isGrad][iScale] = 0;
   if (storeMap[isGrad] > 0) { /* Need to return gradient maps */
    /* Begin: Decide whether gradient response needs to be computed */
    computeMap[isGrad][iScale] = 2;
 
    /* Try opening and reading grad_ file */
    makeCacheName(fname, pathNameCache, gradPrefix, rootNameCache, filterFileTail);

    fpIn = fopenInputFile(fname);

    if ((fpIn=fopenInputFile(fname))!= NULL) {
     /* Begin: Read gradient file, check it */

     computeMap[isGrad][iScale] = 1;  /* found grad results */
     mexPrintf(" Reading gradient data from file %s\n", fname);

     /* Read phase gradient image */
     iChannel = 2 * index(0, iScale, nOrient);
     readFilterOutputs(fpIn, fname, header, HEADER, 
        sizeof(float), nOrient, (int) 2, (int) 0, 
        (char **) ((respPnt[isGrad])+ iChannel),
        &nxTmp, &nyTmp);
     /* readFilterOutputs closes file */

     /* Check phase gradient image */
     checkParam = 1;
     allOrientMaxAmp[iScale] = 0;
     for(iOrient=0; iOrient<nOrient; iOrient++) {
      sscanf(header[iOrient],"%s %e %e %e %e %e", 
             scrName, &pedestal, &scl, &lambdaRead,
             &thetaDegRead, &(maxAmp[iOrient]) );
      if (strcmp(scrName, "freemanGradFlt") != 0)
       checkParam = 0;
      if (ROUND(lambdaRead*100.0) != ROUND(lambda[iScale]*100.0) )
       checkParam = 0;
      if (ROUND(thetaDegRead*100.0) != ROUND(thetaDeg[iOrient]*100.0))
       checkParam = 0;
          
      allOrientMaxAmp[iScale] = max(maxAmp[iOrient], allOrientMaxAmp[iScale]);
     }

     if (checkParam == 0) { /* Gradient files stale */
      mexPrintf(" Gradient files found to be out of date.  Recomputing.\n");
      computeMap[isGrad][iScale] = 2; /* Need to recompute grad results */
      for(iOrient=0; iOrient<nOrient; iOrient++) {
       utilFree((void **) &(respPnt[isGrad][2*iOrient+iChannel]));
       utilFree((void **) &(respPnt[isGrad][2*iOrient+1+iChannel]));
      } 
     } else { /* Gradient files ok, remember size */
      nxFilt[iScale] = nxTmp;
      nyFilt[iScale] = nyTmp;
     }
    } /* End: Read gradient file, check it */
    /* At this point we have pointers to the cached gradient response,
       now stored in the array respPnt, if the file checked out ok. */
   } /* End: Gradient resp needed? */

   computeMap[isPhase][iScale] = 0;
   if (storeMap[isPhase] > 0) { /* Need to return phase maps */
    /* Check out phase results at this scale */
    computeMap[isPhase][iScale] = 2;
 
    /* Try opening and reading phase_ file */
    makeCacheName(fname, pathNameCache, phasePrefix, rootNameCache, filterFileTail);

    if ((fpIn=fopenInputFile(fname)) != NULL) {
     /* Begin: Read phase file, check it */

     computeMap[isPhase][iScale] = 1;  /* found phase results */
     mexPrintf(" Reading phase data from file %s\n", fname);

     /* Read phase image */
     iChannel = index(0, iScale, nOrient);
     readFilterOutputs(fpIn, fname, header, HEADER, 
        sizeof(float), nOrient, (int) 1, (int) 0, 
        (char **) ((respPnt[isPhase])+ iChannel),
        &nxTmp, &nyTmp);
     /* readFilterOutputs closes file */

     /* Check phase image */
     checkParam = 1;
     for(iOrient=0; iOrient<nOrient; iOrient++) {
      sscanf(header[iOrient],"%s %e %e %e %e", 
             scrName, &pedestal, &scl, &lambdaRead,
             &thetaDegRead );
      if (strcmp(scrName, "freemanPhaseFlt") != 0)
       checkParam = 0;
      if (ROUND(lambdaRead*100.0) != ROUND(lambda[iScale]*100.0) )
       checkParam = 0;
      if (ROUND(thetaDegRead*100.0) != ROUND(thetaDeg[iOrient]*100.0))
       checkParam = 0;
     }

     if (checkParam == 0) { /* Phase file stale */
      mexPrintf(" Phase file found to be out of date.  Recomputing.\n");
      computeMap[isPhase][iScale] = 2; /* Need to recompute phase results */
      for(iOrient=0; iOrient<nOrient; iOrient++) {
       utilFree((void **) &(respPnt[isPhase][iOrient+iChannel]));
      } 
     } else { /* Phase files ok, remember size */
      nxFilt[iScale] = nxTmp;
      nyFilt[iScale] = nyTmp;
     }
    } /* End: Read phase file, check it */
    /* At this point we have pointers to the cached phase response,
       now stored in the array respPnt, if the file checked out ok. */
   } /* End: Phase resp needed? */

   computeMap[isAmp][iScale] = 0;
   if (storeMap[isAmp] > 0 || storeMap[isGrad] > 0) {
             /* Need to compute amplitude maps */
    /* Begin: Decide whether amplitude response needs to be computed */
    computeMap[isAmp][iScale] = 2;
 
    /* Try opening and reading _amp_ file */
    makeCacheName(fname, pathNameCache, ampPrefix, rootNameCache, filterFileTail);

    if ((fpIn = fopenInputFile(fname)) != NULL) {
     /* Begin: Read gradient file, check it */

     computeMap[isAmp][iScale] = 1;  /* found amp results */
     mexPrintf(" Reading amplitude data from file %s\n", fname);

     /* Read amplitude image */
     iChannel = index(0, iScale, nOrient);
     readFilterOutputs(fpIn, fname, header, HEADER, 
        sizeof(float), nOrient, (int) 1, (int) 0, 
        (char **) (respPnt[isAmp]+ iChannel),
        &nxTmp, &nyTmp);
     /* readFilterOutputs closes file */

     /* Check amp image */
     checkParam = 1;
     allOrientMaxAmp[iScale] = 0;
     for(iOrient=0; iOrient<nOrient; iOrient++) {
      sscanf(header[iOrient],"%s %e %e %e %e %e", 
             scrName, &pedestal, &scl, &lambdaRead,
             &thetaDegRead, &(maxAmp[iOrient]) );
      if (strcmp(scrName, "freemanAmpFlt") != 0)
       checkParam = 0;
      if (ROUND(lambdaRead*100.0) != ROUND(lambda[iScale]*100.0) )
       checkParam = 0;
      if (ROUND(thetaDegRead*100.0) != ROUND(thetaDeg[iOrient]*100.0))
       checkParam = 0;
          
      allOrientMaxAmp[iScale] = max(maxAmp[iOrient], allOrientMaxAmp[iScale]);
     }

     if (checkParam == 0) { /* Amplitude files stale */
      mexPrintf(" Amplitude files found to be out of date.  Recomputing.\n");
      computeMap[isAmp][iScale] = 2; /* Need to recompute amp results */
      for(iOrient=0; iOrient<nOrient; iOrient++) {
       utilFree((void **) &(respPnt[isAmp][iOrient+iChannel]));
      } 
     } else { /* Amplitude files ok, remember size */
      nxFilt[iScale] = nxTmp;
      nyFilt[iScale] = nyTmp;
     }
    } /* End: Read amplitude file, check it */
    /* At this point we have pointers to the cached amplitude response,
       now stored in the array respPnt, if the file checked out ok. */
   } /* End: Amplitude resp needed? */

   /* Decide whether to compute filter responses */
   computeMap[isFilter][iScale] = 0; 
   if (storeMap[isFilter] > 0)  /* Need to return filter maps */
    computeMap[isFilter][iScale] = 2;
   else if ((computeMap[isAmp][iScale]>=2) ||
            (computeMap[isPhase][iScale]>=2) ||
            (computeMap[isGrad][iScale]>=2)) {
    computeMap[isFilter][iScale] = 2; /* Need to use filter map */
   }
 
   if (computeMap[isFilter][iScale] > 1) { /* Need to recover filter map */
    /* Begin: Decide whether filter response needs to be computed */
    /* Try opening and reading filt_ file */
    makeCacheName(fname, pathNameCache, filterPrefix, rootNameCache, filterFileTail);

    if ((fpIn=fopenInputFile(fname)) != NULL) {
     /* Begin: Read filter file, check it */

     computeMap[isFilter][iScale] = 1;  /* found filter results */
     mexPrintf(" Reading filter data from file %s\n", fname);

     /* Read filter image */
     iChannel = 2 * index(0, iScale, nOrient);
     readFilterOutputs(fpIn, fname, header, HEADER, 
        sizeof(float), nOrient, (int) 2, (int) 0, 
        (char **) (respPnt[isFilter] + iChannel),
        &nxTmp, &nyTmp);
     /* readFilterOutputs closes file */

     /* Check filter image */
     checkParam = 1;
     for(iOrient=0; iOrient<nOrient; iOrient++) {
      sscanf(header[iOrient],"%s %e %e %e %e", 
             scrName, &pedestal, &scl, &lambdaRead, &thetaDegRead );
      if (strcmp(scrName, "freemanFilterFlt") != 0)
       checkParam = 0;
      if (ROUND(lambdaRead*100.0) != ROUND(lambda[iScale]*100.0) )
       checkParam = 0;
      if (ROUND(thetaDegRead*100.0) != ROUND(thetaDeg[iOrient]*100.0))
       checkParam = 0;

     }

     if (checkParam == 0) { /* Filter files stale */
      mexPrintf(" Filter files found to be out of date.  Recomputing.\n");
      computeMap[isFilter][iScale] = 2; /* Need to recompute filter results */
      iChannel = 2 * index(0, iScale, nOrient);
      for(iOrient=0; iOrient<nOrient; iOrient++) {
       utilFree((void **) &(respPnt[isFilter][2*iOrient+iChannel]));
       utilFree((void **) &(respPnt[isFilter][2*iOrient+1+iChannel]));
      } 
     } else { /* Filter files ok, remember size */
      nxFilt[iScale] = nxTmp;
      nyFilt[iScale] = nyTmp;
     }
    } /* End: Read filter file, check it */
    /* At this point we have pointers to the cached filter response,
       now stored in the array respPnt, if the file checked out ok. */
   } /* End: Filter resp needed? */

   /* Decide whether to compute basis responses */
   computeMap[isBasis][iScale] = 0; 
   if (storeMap[isBasis] > 0)  /* Need to return basis maps */
    computeMap[isBasis][iScale] = 2;
   else if (computeMap[isFilter][iScale]>=2) {
    computeMap[isBasis][iScale] = 2; /* Need to use basis map */
   }

   if (computeMap[isBasis][iScale] > 1) { /* Need to recover basis map */
    /* Begin: Decide whether basis response needs to be computed */
    /* Try opening and reading basis_ file */
    makeCacheName(fname, pathNameCache, basisPrefix, rootNameCache, filterFileTail);

    if ((fpIn=fopenInputFile(fname)) != NULL) {
     /* Begin: Read basis file, check it */

     computeMap[isBasis][iScale] = 1;  /* found basis results */
     mexPrintf(" Reading basis data from file %s\n", fname);

     /* Read basis image */
     iChannel = index(0, iScale, NUMBASIS);
     readFilterOutputs(fpIn, fname, header, HEADER, 
        sizeof(float), (int) 1, NUMBASIS, (int) 0, 
        (char **) (respPnt[isBasis]+ iChannel),
        &nxTmp, &nyTmp);
     /* readFilterOutputs closes file */

     /* Check basis image */
     checkParam = 1;
     sscanf(header[0],"%s %e %e %e", 
             scrName, &pedestal, &scl, &lambdaRead );
     if (strcmp(scrName, "freemanBasisFlt") != 0)
       checkParam = 0;
     if (ROUND(lambdaRead*100.0) != ROUND(lambda[iScale]*100.0) )
       checkParam = 0;
     if (checkParam == 0) { /* Basis files stale */
      mexPrintf(" Basis files found to be out of date.  Recomputing.\n");
      computeMap[isBasis][iScale] = 2; /* Need to recompute basis results */
      for(k=0; k<NUMBASIS; k++)
       utilFree((void **) &(respPnt[isBasis][k+iChannel]));
     } else { /* Basis files ok, remember size */
      nxFilt[iScale] = nxTmp;
      nyFilt[iScale] = nyTmp;
     }
    } /* End: Read basis file, check it */
    /* At this point we have pointers to the cached basis response,
       now stored in the array respPnt, if the file checked out ok. */
   } /* End: Basis resp needed? */

  } /* Over scales */

 } /* End: checkPhaseCache */

/* Check the cache for existing Gaussian pyramid images */
void checkPyramidCache(float **pyramid, int *localPyrBase, 
		       int *computePyr, 
		       int nPyr, int *nxPyr, int *nyPyr, 
		       char *pathNameCache, char *rootNameCache)
 {
  float pedestal, scl, tmp;
  char fname[MAXLEN], scrName[MAXLEN];
  char pyrFileTail[10], header[NUMORIENT][HEADER];
  int k, sampleRate, checkParam;
  FILE *fpIn;

  *localPyrBase = 0;
  /* Propagate dependencies backwards in scale */
  for(k=nPyr-2; k>=0; k--) {
   if (computePyr[k+1] >= 2) computePyr[k] = 2;
  }

  sampleRate = 1;
  for(k=0; k<nPyr; k++) {
   /* Decide whether to compute each pyramid response */
   if (computePyr[k] > 1) {  /* Need to return pyramid maps */

    sprintf(pyrFileTail,".s%d",sampleRate);
   
    /* Try opening and reading pyramid_ file */
    makeCacheName(fname, pathNameCache, pyrPrefix, rootNameCache, pyrFileTail);

    if ((fpIn=fopenInputFile(fname)) != NULL) {
     /* Begin: Read pyramid file, check it */

     computePyr[k] = 1;  /* found pyramid results */
     mexPrintf(" Reading pyramid data from file %s\n", fname);

     /* Read pyramid image */
     readFilterOutputs(fpIn, fname, header, HEADER, 
        sizeof(float), (int) 1, (int) 1, (int) 0, 
        (char **) (pyramid + k), nxPyr+k, nyPyr+k);
     /* readFilterOutputs closes file */

     /* Check pyramid image */
     checkParam = 1;
     sscanf(header[0],"%s %e %e %e", 
             scrName, &pedestal, &scl, &tmp );
     if (strcmp(scrName, "pyramidFlt") != 0)
       checkParam = 0;
     if (ROUND(tmp) != sampleRate)
       checkParam = 0;

     if (checkParam == 0) { /* Pyramid files stale */
      mexPrintf(" Pyramid file %s found to be out of date. Recomputing.\n",
          fname);
      computePyr[k] = 2; /* Need to recompute pyramid image */
      utilFree((void **) &(pyramid[k]));
     } else {
      if (k==0) *localPyrBase = 1;
     }

    } /* End: Read pyramid file, check it */
    /* At this point we have pointers to the cached pyramid response,
       now stored in the array pyramid, if the file checked out ok. */
   } /* End: Pyramid resp needed? */
   sampleRate *= 2;
  } /* Over scales */

} /* End: checkPyramidCache() */

/* Compute params of Freeman filter given the peak wavelength */
void freemanConstants(float lambda, int *pdwnSmpl, int *psubSample, 
		      float *psigma, float *prho)
 {
  float band, tmp;
  int subSample0, pyrLevel;
  /* Determine sampling rate for filters, namely lambda/4 */
  subSample0 = (int) (lambda/4.0 + 0.0001);
  subSample0 = max(1, subSample0);
  
  /* Using an image pyramid, we need the subsampling rate to
     be subSample0 >= 2^(dwnSmpl+1) */
  if (subSample0 >= 4) {

    tmp = log((double)(subSample0/2.0))/lnTwo;
    pyrLevel = (int) ( tmp + 5.0 * MACHINE_EPS );
    tmp = pow( (double)2.0, (double) pyrLevel);
    *pdwnSmpl = ROUND(tmp);
  } else {
    *pdwnSmpl = 1;
  }
  *psubSample = subSample0/(*pdwnSmpl); 
                /* subSample0 should be 2^k for efficiency */
  *psigma = lambda / lamToSig;  /* Constant chosen to give mean frequency
                                  having wavelength lambda (see g2h2_1D.m) */
  *prho = rhoThetaScale; 
  band = 1.00; /* Empirically, based on 1 sigma contours
                  of the Fourier transform for wavelengths 8 or larger */

 }

/* Correlate (slow 2D method) a Freeman filter with an image */
void applyFreeman(void *vimage, int bytesPerPixel, 
		  int nxImage, int nx0, int ny0, 
		  float lambda, float thetaDeg, 
		  float *respR, float *respI, int nxResp, 
		  float *pMaxAmp)
 {
  float *filtR, *filtI; 
  float lambdaDwn, sigma, sigmaDwn, rho;
  float bound[2];
  int subSample, ds, dwnSmpl, nx, ny;
  int filtSize;

  freemanConstants(lambda, &dwnSmpl, &subSample, &sigma, &rho);
  lambdaDwn = lambda/dwnSmpl;
  sigmaDwn = sigma/dwnSmpl;

  ds = subSample * dwnSmpl;
  nx = nx0/ds;
  ny = ny0/ds;


  build_freeman_filter(&filtR, &filtI, &filtSize, 
                        thetaDeg, lambdaDwn, sigmaDwn, rho);

  *pMaxAmp = 0.0; 
  if (bytesPerPixel == 1) {
   corrByte((unsigned char *)vimage,dwnSmpl, nxImage, nx0,ny0,subSample,filtR,filtSize,
         respR, nxResp, bound);
  } else {
   corrFlt((float *) vimage, dwnSmpl, nxImage, nx0,ny0,subSample,filtR,filtSize,
         respR, nxResp, bound);
  }
  *pMaxAmp = max( (*pMaxAmp), abs(bound[0]));
  *pMaxAmp = max( (*pMaxAmp), abs(bound[1]));
  if (bytesPerPixel == 1) {
   corrByte((unsigned char *)vimage,dwnSmpl, nxImage, nx0,ny0,subSample,filtI,filtSize,
         respI, nxResp, bound);
  } else {
   corrFlt((float *) vimage, dwnSmpl, nxImage, nx0,ny0,subSample,filtR,filtSize,
         respI, nxResp, bound);
  }
  *pMaxAmp = max( (*pMaxAmp), abs(bound[0]));
  *pMaxAmp = max( (*pMaxAmp), abs(bound[1]));
    
 }

/* Allocate and generate separable Freeman filters */
void build_sepFreeman_filter(float **G2x, float **G2y, 
			     float **H2x, float **H2y, 
			     int *pfilt_size, 
			     float lambda, float sigma, float rho)
{
  int   i,j, filtSize, filtSize2 ;
  int   sizeMem;
  float x, x2, ex2, omega;

  omega = 2.0*PI/lambda;
  filtSize = 2 * ROUND(filterDiam/2.0 * sigma) + 1;
  filtSize2 = filtSize/2;
  *pfilt_size = filtSize;

  /* Allocate space for filter coefficients */
  sizeMem = filtSize * sizeof(float);
  for (i=0; i<3; i++) {
   grabFloatMemory(G2x+i, sizeMem, "G2x_filt");
   grabFloatMemory(G2y+i, sizeMem, "G2y_filt");
  }
  for (i=0; i<4; i++) {
   grabFloatMemory(H2x+i, sizeMem, "H2x_filt");
   grabFloatMemory(H2y+i, sizeMem, "H2y_filt");
  }

  for (i=0; i<filtSize; i++) {

    x   = (i - filtSize2)/sigma ;
    x2  = x * x  ;
    ex2 = exp( - x2 / 2.0 );

    /************************************************************************/
    /* It may appear that I've generated the filter masks in reverse order  */
    /* as compared to Freeman's paper, Tables III & IV, but Allan's code    */
    /* reverses the order of X and Y (together with ommitting the -ve signs */
    /* from the weights from those tables), so I've chosen to follow suit.  */
    /************************************************************************/
     
    G2x[2][i] = sqrt(0.9213) * (x2 - 1.0) * ex2 ;
    G2x[1][i] = sqrt(0.9213) * x          * ex2 ;
    G2x[0][i] = sqrt(0.9213)              * ex2 ;

    G2y[0][i] = G2x[2][i] ;
    G2y[1][i] = G2x[1][i] ;
    G2y[2][i] = G2x[0][i] ;

    H2x[3][i] = (-0.3458 * x2 + 1.559) * x * ex2 ;
    H2x[2][i] = (-0.3458 * x2 + 1.559/3.0) * ex2 ;
    H2x[1][i] = x * ex2  ;
    H2x[0][i] = ex2      ;

    H2y[0][i] = H2x[3][i];
    H2y[1][i] = H2x[2][i];
    H2y[2][i] = H2x[1][i];
    H2y[3][i] = H2x[0][i];
  }

  /* Normalize according to plane wave, wave number omega */
  {
      float tmpG2 = 0.0, tmpH2 = 0.0, sum = 0.0, tmp, y ;
      int   k = 0 ;

      for (i=0; i<filtSize; i++) 
        for (j=0; j<filtSize; j++) {
            y = j - filtSize2 ;
            /* apply plane wave for scaling along dominant */
            /* direction of filter                         */
            tmp = sin( omega * y);
            tmpG2 += G2x[k][i] * G2y[k][j] * tmp ;
            sum   += G2x[k][i] * G2y[k][j];
          }

      for (i=0; i<filtSize; i++) 
        for (j=0; j<filtSize; j++) {
            y = j - filtSize2 ;
            tmp = sin( omega *y );
            tmpH2 += H2x[k][i] * H2y[k][j] * tmp ;
          }

      tmp = 1.0 / sqrt(tmpG2 * tmpG2 + tmpH2 * tmpH2);
      tmp = sqrt(tmp); /* take sqrt() since the masks multiply each other */

      for (k=0; k<3; k++) {
        for (i=0; i<filtSize; i++) {
            G2x[k][i] *= tmp ;
            G2y[k][i] *= tmp ;
          }
        }

      for (k=0; k<4; k++) {
        for (i=0; i<filtSize; i++) {
            H2x[k][i] *= tmp ;
            H2y[k][i] *= tmp ;
          }
        }

  }

} /* End: build_sepFreeman_filter() **/

/* Allocate and generate a 2D Freeman filter */
void build_freeman_filter(float **pFreemanR, float **pFreemanI, 
			  int *pfilter_size, 
                          float thetaDeg, float lambda,
			  float sigma, float rho)
{
   float omega, x[2], pair[2];
   float theta, ctheta, stheta, sum, tmp, tmpI, tmpR, val;
   int ix, iy, ixy, filtSize, filtSize2, sizeOutput;

   omega = 2.0*PI/lambda;
   filtSize = 2 * ROUND(filterDiam/2 * sigma) + 1; 
   filtSize = (filtSize & ~0x1) + 1;
   filtSize2 = filtSize/2;
   *pfilter_size = filtSize;

   sizeOutput = sizeof(float) * filtSize * filtSize;
   /* Allocate space for filter coefficients */
   grabFloatMemory(pFreemanR, sizeOutput, "freeman_R_filt");
   grabFloatMemory(pFreemanI, sizeOutput, "freeman_I_filt");
   theta = PI * thetaDeg / 180.0;
     
   ctheta = cos(theta); stheta =sin(theta);
   sum=tmpI=tmpR=0.0;
   for(iy=0; iy<filtSize; iy++)
    for(ix=0; ix<filtSize; ix++) {
        x[0] = iy - filtSize2; x[1] = ix-filtSize2;
        freeman(x, sigma, theta, pair);
        ixy = index(ix, iy, filtSize);
        (*pFreemanR)[ixy] = pair[0];
        (*pFreemanI)[ixy] = -pair[1];  /* Use complex congugate filter, then correlate */
        sum += (*pFreemanR)[ixy];
        /* Use plane wave of wavelength lambda, amp 1 for scaling purposes */
        val = sin(omega * (ctheta * x[0] + stheta * x[1]));
        tmpR += (*pFreemanR)[ixy] * val;
        tmpI += (*pFreemanI)[ixy] * val;
    }
   tmp = sqrt(tmpR * tmpR + tmpI * tmpI);
   tmp = 1.0/tmp;
   for(iy=0; iy<filtSize; iy++) {
    for(ix=0; ix<filtSize; ix++) {
        /*  Rescale to have gain of 1.0 for wavelength lambda, and orientation
           theta   */
        ixy = index(ix, iy, filtSize);
        (*pFreemanR)[ixy] *= tmp;
        (*pFreemanI)[ixy] *= tmp;
    }
  }
}
 
/* Sample freeman pair G2/H2 at a tap position x[2] */ 
void freeman(float *x, float sigma, float theta, float *pair)
{
   /* Freeman pair G2, H2.  Use envelop of the form x^2/(2.0*sigma^2)
      instead of (x/sigma)^2.  This changes some coefficients.  Also,
      in Freeman's thesis there is an extra factor of two in the real
      part */
   float windowFun, rot_x[2], x2;

   rot_x[0] = (x[0] * cos( theta ) + x[1] * sin( theta ))/sigma;
   rot_x[1] = (x[0] * -sin( theta ) + x[1] * cos( theta ))/sigma;
   windowFun = -( x[0] * x[0] + x[1] * x[1] )/ (2.0 * sigma * sigma);
   if ( windowFun < -20.0 )      
    windowFun = 0.0;
   else {
    windowFun = exp(windowFun);
   }
   x2 = rot_x[0] * rot_x[0];
   pair[0] = 0.9213 * ( x2 - 1.0) * windowFun;
   pair[1] = rot_x[0] * (0.3458 * x2 - 1.559) * windowFun;
}

/* Compute the 5-point discrete derivative operator for
   a band-pass filter channel tuned for peak wavelength lambda
   and orientation thetaDeg.  Also estimate the inverse covariance matrix C
   for the frequency response. */
void derivMask(float dPhi[2], float C[2][2], 
	       float maskX[5][2], float maskY[5][2],
	       float lambda, float thetaDeg,
	       float sigma, float rho, int subSample)
 {
      float dxRot[2], dyRot[2], tmp;
      float omega, theta, ctheta, stheta;
      float sigmaTheta, sigmaPerp, varTheta, varPerp;
      int i,j;
     
       omega = 2.0*PI/lambda;
       theta = PI * thetaDeg / 180.0;
       ctheta = cos(theta); stheta =sin(theta);

       sigmaTheta = sigma;
       sigmaPerp = sigma * rho;
       varTheta =  sigmaTheta * sigmaTheta;
       varPerp =  sigmaPerp * sigmaPerp;

       /* Compute singularity neighbourhood filter parameters.
          C = inverse covariance matrix for bandpass filters. 
          C = sigma^2 * [ I + (rho^2-1) (s c)^T (s c) ] */
       C[0][0] = ctheta * ctheta * varTheta + stheta * stheta * varPerp ;
       C[0][1] = C[1][0] =  stheta * ctheta * ( varPerp - varTheta );
       C[1][1] =  stheta * stheta * varTheta + ctheta * ctheta * varPerp ;
 
       /*  Set up mean frequency dPhi, relative to subsampled grid */
       /*  Here dPhi[0] is phase difference in x direction, dPhi[1]
           is frequency in y direction.
         Theta = 0  - i.e. corresponds to horizontal line
                 45 /
                 90 |
                 135 \
         dPhi[*] = Expected phase difference in x,y directions
                   on one step of subsampled grid.
       **************/
       dPhi[0] = omega*stheta*subSample;
       dPhi[1] = omega*ctheta*subSample;
 
       /* Demodulate by counter rotating by amount dPhi[*] */
       dxRot[0] = cos( -dPhi[0]); dxRot[1] = sin( -dPhi[0] );
       dyRot[0] = cos( -dPhi[1] ); dyRot[1] = sin( -dPhi[1] );
     
       for(i=0; i<2; i++) {
        maskX[2][i] = 0.0;
        maskX[3][i] = 8.0 * dxRot[i];
        maskY[2][i] = 0.0;
        maskY[3][i] = 8.0 * dyRot[i];
       }
       maskX[4][0] = -(dxRot[0] * dxRot[0] - dxRot[1] * dxRot[1]);
       maskX[4][1] = -(2.0 * dxRot[0] * dxRot[1]);
       maskY[4][0] = -(dyRot[0] * dyRot[0] - dyRot[1] * dyRot[1]);
       maskY[4][1] = -(2.0 * dyRot[0] * dyRot[1]);
       for(i=-2; i<=-1; i++) {
        maskX[i+2][0] = -maskX[2-i][0];
        maskX[i+2][1] = maskX[2-i][1];
        maskY[i+2][0] = -maskY[2-i][0];
        maskY[i+2][1] = maskY[2-i][1];
       }
       tmp = 12.0 * subSample;
       for(i=0; i<=4; i++) 
        for(j=0; j<2; j++) {
         maskX[i][j] /= tmp;
         maskY[i][j] /= tmp;
        }
 }
   
/***************************************
 Input: 
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
		       float *phiX, float *phiY, float *maxAmp )
 {
 
  float thres, amp, amp2, tmp, respCos, respSin;
  float dPhi[2], dAmp[2], omega, ctheta, stheta, theta, thetaOrient;
  float dFx[2], dFy[2];
  int sx[5], sy[5], i, j, ix, iy, ixy;
  
  omega = 2.0 * PI / lambda;
  theta = thetaDeg * PI / 180.0;
  ctheta = cos(theta); stheta = sin(theta);
  *maxAmp = 0.0;
  
  /* Loop through image */
  for(iy=0; iy<ny; iy++) 
   for(ix=0; ix<nx; ix++) {
    ixy = index(ix, iy, nx);
    phiX[ixy] = phiY[ixy] = BIG_NEG;
    
    /* Check border condition */
    if (ix < border || ix >= nx-border ||
        iy < border || iy >= ny-border) {
      phiX[ixy] = phiY[ixy] = BIG_NEG; /* BIG_NEG Signifies no valid response */
    } else { /* pixel ix, iy is sufficiently far from image border */
    
      respCos = gCos[ixy];
      respSin = gSin[ixy];
      amp2 = respCos*respCos + respSin*respSin;
      amp = sqrt(amp2) ;

      if (amp <= epsAmp) {
       phiX[ixy] = phiY[ixy] = BIG_NEG; /* Signifies no valid value */
      } else {
       
       /* Compute phase gradient */

       /* Build indirect references for neigbouring pixels */
       for(i=0; i<=4; i++) {
        sx[i] = (ix + i - 2);
        sx[i] = max(sx[i],0);
        sx[i] = min(sx[i],nx-1);
        sy[i] = (iy + i - 2);
        sy[i] = max(sy[i],0);
        sy[i] = min(sy[i],ny-1);
       }
       for(i=0; i<2; i++)
        dFx[i] = dFy[i] = 0.0;
       for(i=0; i<=4; i++) {
        j = index(sx[i], iy, nx);
        dFx[0] += maskX[i][0] * gCos[j] - maskX[i][1] * gSin[j];
        dFx[1] += maskX[i][1] * gCos[j] + maskX[i][0] * gSin[j];
        j = index(ix, sy[i], nx);
        dFy[0] += maskY[i][0] * gCos[j] - maskY[i][1] * gSin[j];
        dFy[1] += maskY[i][1] * gCos[j] + maskY[i][0] * gSin[j];
       }       

       dPhi[0] = (dFx[1] * respCos - dFx[0] * respSin )/amp2 ;
       dPhi[1] = (dFy[1] * respCos - dFy[0] * respSin )/amp2 ;

       dAmp[0] = (dFx[0] * respCos + dFx[1] * respSin )/amp2;
       dAmp[1] = (dFy[0] * respCos + dFy[1] * respSin )/amp2;

       /****   Amplitude component of threshold ***/
       tmp = dAmp[0] * ( C[0][0]*dAmp[0] + C[0][1]*dAmp[1] );
       tmp += dAmp[1] * ( C[1][0]*dAmp[0] + C[1][1]*dAmp[1] );
       thres = tmp;
       /****   Phase component of threshold ***/
       tmp = dPhi[0] * ( C[0][0]*dPhi[0] + C[0][1]*dPhi[1] );
       tmp += dPhi[1] * ( C[1][0]*dPhi[0] + C[1][1]*dPhi[1] );
       thres += tmp;

       thres = sqrt(thres);

       /*  Gradient in units radians/pixel divided by PI */
       phiX[ixy] = (omega*stheta + dPhi[0])/PI;
       phiY[ixy] = (omega*ctheta + dPhi[1])/PI;
       
       if (thres > tau) /* Fails singularity test */
        phiX[ixy] = phiY[ixy] = BIG_NEG;
       else {
        /* Check orientation is within range */
        thetaOrient = atan2( phiX[ixy] , phiY[ixy])/PI*180.0;
        /** get direction modulo 180 **/
        while(thetaOrient - thetaDeg < -90.0) thetaOrient += 180.0;
        while(thetaOrient - thetaDeg > 90.0) thetaOrient -= 180.0;
        if ( abs(thetaOrient - thetaDeg) > thetaTol ) { /* in Degrees */
          phiX[ixy] = phiY[ixy] = BIG_NEG; /* Remove measurement */
        } else
         *maxAmp = max( *maxAmp, amp );
         
       } /* passes singularity test */
      } /* passes amplitude test */
     } /* passes image border test */
   } /* loop over ix, iy */
 }

/* Set cacheName to be the concatenation:
      pathName + prefix + root + postfix
   It is assumed that cachName has been declared to be sufficiently long
   to hold this concatenated string. */
void makeCacheName(char *cacheName, 
		   char *pathName, char *prefix, char *root, char *postfix)
 {
  strcpy(cacheName, pathName);
  strcat(cacheName, prefix);
  strcat(cacheName, root);
  strcat(cacheName, postfix);
 }

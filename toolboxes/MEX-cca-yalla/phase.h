#ifndef PI
# define        PI      	M_PI
#endif

#ifndef _PHASE_H
# define        _PHASE_H

# define        halfPI      	1.5707963
# define        twoPI      	6.2831853
# define        lnTwo      	0.6931472
# define        DIM     	512 /***** Mac 128  SGI 512  *****/
# define        HEADER      	512
# define	NUMORIENT	4
# define	NUMSCALE  	6
# define	NUMCHANNEL  	24  /* NUMSCALE * NUMORIENT */
# define	NUMBASIS	7  /* number of basis elements in sep freeman */
# define	NUMTEMPLATE  	10

# define        NO_VALUE        0
# define        NULL_NAME_FLAG  -2
# define        MAX_NHBR        90 /* number of elements in previous nhbr-hood */

# define	ROCK_BOTTOM	0.1  /* minimum non-neglibible filter amp */

# define	isGrad		0  /* Codes for storeMap, cachMap */
# define	isPhase		1
# define	isAmp		2
# define	isFilter	3
# define	isBasis		4
# define	isPyr		5
# define	NUMCODE		6
# define        outlierIndex    255
#endif

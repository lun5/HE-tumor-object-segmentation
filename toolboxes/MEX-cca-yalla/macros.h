
#include        "stdio.h"
#include        <stdlib.h>
#include        <math.h>
#include        <string.h>

#ifndef PI
#ifdef M_PI
 # define        PI              M_PI
#else
 # define        PI              3.14159265358979323846
#endif
#endif

#ifndef TRUE
# define TRUE 1
#endif
#ifndef FALSE
# define FALSE 0
#endif

#ifndef _MY_MACROS_H
#define _MY_MACROS_H

# define        ROUND(A)        ((A)<0.0?-((int)(-(A)+0.5)):((int)((A)+0.5)))
# define        max(A,B)        ((A)<(B)?(B):(A))
# define        min(A,B)        ((A)<(B)?(A):(B))
# define        sign(A)		((A)<0?(-1):(A)>0?(1):0)
# define        giveSign(A,B)	((A)<0?-(B):(B))
# define        abs(A)		((A)<0?-(A):(A))
# define	index(IX, IY, NX)	((IY)*(NX)+(IX))
# define	index3d(i,j,k, n1, n2) ((i)+(n1)*((j)+(k)*n2))
# define	linearInterp(x0, x1, dx) ((x0)*(1-(dx))+(x1)*(dx))

# define 	MACHINE_EPS 	0.5e-6
# define 	MAXLEN 	        256
# define        BIG_NEG         (float)-1.0e+38

#endif

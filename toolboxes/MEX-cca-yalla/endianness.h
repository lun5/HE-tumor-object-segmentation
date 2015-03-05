/*
  FILE: endianness.h   AJ Dec.10, 1999
 */
#ifndef	ENDIANNESS_H

#define	ENDIANNESS_H

/* This file checks to see if _BIG_ENDIAN or _LITTLE_ENDIAN 
   is defined.  If not, it tries to determine the appropriate
   one to define.  It causes a compile time error if it gets
   confused, or if you are running on a PDP_ENDIAN machine.
*/

#if !(defined(_BIG_ENDIAN) || defined(_LITTLE_ENDIAN) || defined(_PDP_ENDIAN))

#ifndef BYTE_ORDER
#if (BSD >= 199103)
# include <machine/endian.h>
#else
#ifdef linux
# include <endian.h>
#else

#ifndef LITTLE_ENDIAN
#define LITTLE_ENDIAN   1234    /* least-significant byte first (vax, pc) */
#endif

#ifndef BIG_ENDIAN
#define BIG_ENDIAN      4321    /* most-significant byte first (IBM, net) */
#endif

#ifndef PDP_ENDIAN
#define PDP_ENDIAN      3412    /* LSB first in word, MSW first in long (pdp)*/
#endif

#if defined(vax) || defined(ns32000) || defined(sun386) || defined(i386) || \
    defined(MIPSEL) || defined(_MIPSEL) || defined(BIT_ZERO_ON_RIGHT) || \
    defined(__alpha__) || defined(__alpha)
#define BYTE_ORDER      LITTLE_ENDIAN
#endif

#if defined(sel) || defined(pyr) || defined(mc68000) || defined(sparc) || \
    defined(is68k) || defined(tahoe) || defined(ibm032) || defined(ibm370) || \
    defined(MIPSEB) || defined(_MIPSEB) || defined(_IBMR2) || defined(DGUX) ||\
    defined(apollo) || defined(__convex__) || defined(_CRAY) || \
    defined(__hppa) || defined(__hp9000) || \
    defined(__hp9000s300) || defined(__hp9000s700) || \
    defined (BIT_ZERO_ON_LEFT) || defined(m68k)
#define BYTE_ORDER      BIG_ENDIAN
#endif
#endif /* linux */
#endif /* BSD */
#endif /* BYTE_ORDER */

#if !defined(BYTE_ORDER) || \
    (BYTE_ORDER != BIG_ENDIAN && BYTE_ORDER != LITTLE_ENDIAN && \
    BYTE_ORDER != PDP_ENDIAN)
        /* you must determine what the correct bit order is for
         * your compiler - the next line is an intentional error
         * which will force your compiles to bomb until you fix
         * the above macros.
         */
  error "Undefined or invalid BYTE_ORDER";
#endif
/* BYTE_ORDER has now been defined */

#if (BYTE_ORDER == BIG_ENDIAN)
#define _BIG_ENDIAN
#elif (BYTE_ORDER == LITTLE_ENDIAN)
#define _LITTLE_ENDIAN
#elif (BYTE_ORDER == PDP_ENDIAN)
#define _PDP_ENDIAN
#endif

#endif  /* _BIG/_LITTLE/_PDP ENDIAN not defined */

#ifdef _PDP_ENDIAN
        /* The next line is an intentional error
         * which will force your compiles to bomb until you fix
         * imageFile-io.*.c to deal with PDP_ENDIAN float pfm images.
         */
  error "Byte swapping not currently defined for PDP_ENDIAN";
#endif

#endif  /* ENDIANNESS_H */

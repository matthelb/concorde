/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

#ifndef  __BIGGUY_H
#define  __BIGGUY_H

#include "util.h"

#undef CC_BIGGUY_LONG
#undef CC_BIGGUY_LONGLONG

#ifdef  CC_BIGGUY_LONGLONG
typedef long long CCbigguy;
#define CC_BIGGUY_BUILTIN
#endif

#ifdef CC_BIGGUY_LONG
typedef long CCbigguy;
#define CC_BIGGUY_BUILTIN
#endif

#ifdef CC_BIGGUY_BUILTIN

#define CCbigguy_FRACBITS 32
#define CCbigguy_DUALSCALE (((CCbigguy) 1) << CCbigguy_FRACBITS)
#define CCbigguy_FRACPART(x) ((x) & (CCbigguy_DUALSCALE-1))
#define CCbigguy_MAXBIGGUY (((((CCbigguy) 1) << 62) - 1) + \
                            (((CCbigguy) 1) << 62))
#define CCbigguy_MINBIGGUY (-CCbigguy_MAXBIGGUY)
#define CCbigguy_bigguytod(x) (((double) (x)) / ((double) CCbigguy_DUALSCALE))
#define CCbigguy_itobigguy(d) ((CCbigguy) ((d) * (double) CCbigguy_DUALSCALE))
#define CCbigguy_ceil(x) (CCbigguy_FRACPART(x) ? \
        ((x) + (CCbigguy_DUALSCALE - CCbigguy_FRACPART(x))) : (x))
#define CCbigguy_cmp(x,y) (((x) < (y)) ? -1 : ((x) > (y)) ? 1 : 0)
#define CCbigguy_ZERO ((CCbigguy) 0)
#define CCbigguy_ONE ((CCbigguy) CCbigguy_DUALSCALE)
#define CCbigguy_addmult(x,y,m) ((*x) += (y)*(m))
#define CCbigguy_dtobigguy(d) ((CCbigguy) ((d) * (double) CCbigguy_DUALSCALE))

#else /* CC_BIGGUY_BUILTIN */

typedef struct CCbigguy {
    unsigned short ihi;
    unsigned short ilo;
    unsigned short fhi;
    unsigned short flo;
} CCbigguy;

extern const CCbigguy CCbigguy_MINBIGGUY;
extern const CCbigguy CCbigguy_MAXBIGGUY;
extern const CCbigguy CCbigguy_ZERO;
extern const CCbigguy CCbigguy_ONE;


    void
        CCbigguy_addmult (CCbigguy *x, CCbigguy y, int m);

    int
        CCbigguy_cmp (CCbigguy x, CCbigguy y);

    double
        CCbigguy_bigguytod (CCbigguy x);

    CCbigguy
        CCbigguy_itobigguy (int d),
        CCbigguy_dtobigguy (double d),
        CCbigguy_ceil (CCbigguy x);


#endif /* CC_BIGGUY_BUILTIN */

#define CCbigguy_add(x,y) (CCbigguy_addmult(x,y,1))
#define CCbigguy_sub(x,y) (CCbigguy_addmult(x,y,-1))


int
    CCbigguy_swrite (CC_SFILE *f, CCbigguy x),
    CCbigguy_sread (CC_SFILE *f, CCbigguy *x);


#endif /* __BIGGUY_H */

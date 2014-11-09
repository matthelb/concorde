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

#ifndef __FMATCH_H
#define __FMATCH_H

#include "util.h"


int
    CCfmatch_fractional_2match (int ncount, int ecount, int *elist, int *elen,
        CCdatagroup *dat, double *val, int *thematching, int *thedual,
        int *thebasis, int wantbasic, int silent, CCrandstate *rstate);


#endif  /* __FMATCH_H */

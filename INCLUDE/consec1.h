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

#ifndef  __CONSEC1_H
#define  __CONSEC1_H

#include "tsp.h"

int
    CCpq_consecutiveones (CCtsp_lpcut_in **cuts, int *cutcount,
        CCtsp_cuttree *ctree, CCtsp_lpcuts *pool, int ecount, int *elist,
        double *x);


#endif /* __CONSEC1_H */

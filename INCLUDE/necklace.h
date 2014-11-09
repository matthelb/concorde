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

#ifndef  __NECKLACE_H
#define  __NECKLACE_H

#include "tsp.h"

int
    CCpq_necklaces (CCtsp_lpcut_in **cuts, int *cutcount, CCtsp_cuttree *ctree,
        int ecount, int *elist, double *x, CCrandstate *rstate);


#endif /* __NECKLACE_H */

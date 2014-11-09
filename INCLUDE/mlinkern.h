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

#ifndef __MLINKERN_H
#define __MLINKERN_H

#include "util.h"
#include "kdtree.h"

int
    CCedgegen_mlinkern (int ncount, CCdatagroup *dat, int wantlist,
        int *ecount, int **elist, CCkdtree *kt, int iterations,
        CCrandstate *rstate);

#endif

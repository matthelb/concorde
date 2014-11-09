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

#ifndef __HELDKARP_H
#define __HELDKARP_H

#define HELDKARP_ERROR               -1
#define HELDKARP_SEARCHLIMITEXCEEDED  1

#include "util.h"



int
    CCheldkarp_small (int ncount, CCdatagroup *dat, double *upbound,
             double *optval, int *foundtour, int anytour, int *tour_elist,
             int nodelimit, int silent),
    CCheldkarp_small_elist (int ncount, int ecount, int *elist, int *elen,
             double *upbound, double *optval, int *foundtour, int anytour,
             int *tour_elist, int nodelimit, int silent);


#endif  /* __HELDKARP_H */

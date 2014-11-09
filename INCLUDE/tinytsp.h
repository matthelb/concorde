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

#ifndef __TINYTSP_H
#define __TINYTSP_H

#include "util.h"

#define CC_TINYTSP_ERROR               -1
#define CC_TINYTSP_SEARCHLIMITEXCEEDED  1
#define CC_TINYTSP_INFEASIBLE           2

#define CC_TINYTSP_MAXIMIZE (-1)
#define CC_TINYTSP_MINIMIZE (1)


int
    CCtiny_bnc_tsp (int ncount, CCdatagroup *dat, double *upbound,
        double *optval, int nodelimit),
    CCtiny_bnc_msp (int ncount, int ecount, int *elist, int *elen, int depot,
        int *lower, int *upper, double *upperbound, int objsense,
        double *optval, int *xsol, int checkresult, int searchlimit),
    CCtiny_bnb_tsp (int nnodes, int nedges, int *elist, int *weight,
        int *lbound, int *ubound, double *objlimit, int objdir,
        double *objval, int *xsol, int searchlimit),
    CCtiny_bnb_msp (int nnodes, int nedges, int *elist, int *weight, int depot,
        int *lbound, int *ubound, double *objlimit, int objdir,
        double *objval, int *xsol, int searchlimit),
    CCtiny_benttsp_elist (int ncount, int ecount, int *elist, int *elen,
        double *upbound, double *optval, int *foundtour, int anytour,
        int searchlimit, int silent);


#endif  /* __TINYTSP_H */

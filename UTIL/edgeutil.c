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

/****************************************************************************/
/*                                                                          */
/*                    EDGE LIST UTILITY ROUTINES                            */
/*                                                                          */
/*                              TSP CODE                                    */
/*                                                                          */
/*  Written by: Applegate, Bixby, Chvatal, and Cook                         */
/*  Date: February 8, 1995                                                  */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_edge_to_cycle (int ncount, int *elist, int *yesno,           */
/*      int *cyc)                                                           */
/*    CONVERTS an edgelist to a cycle.                                      */
/*     -ncount is the number of nodes.                                      */
/*     -elist is an edgelist in end1 end2 format.                           */
/*     -yesno returns 1 if elist describes a tour and 0 otherwise.          */
/*     -cyc returns the cycle in permutation format if it is not NULL       */
/*      (if cyc is not NULL, then it should point to an array of            */
/*      length at least ncount).                                            */
/*     Returns a nonzero value if there was an error.                       */
/*                                                                          */
/*  void CCelist_init (CCelist *elist)                                      */
/*     initialize a CCelist                                                 */
/*                                                                          */
/*  void CCelistl_init (CCelistl *elist)                                    */
/*     initialize a CCelistl                                                */
/*                                                                          */
/*  void CCelistw_init (CCelistw *elist)                                    */
/*     initialize a CCelistw                                                */
/*                                                                          */
/*  void CCelistlw_init (CCelistlw *elist)                                  */
/*     initialize a CCelistlw                                               */
/*                                                                          */
/*  void CCelist_free (CCelist *elist)                                      */
/*     free a CCelist                                                       */
/*                                                                          */
/*  void CCelistl_free (CCelistl *elist)                                    */
/*     free a CCelistl                                                      */
/*                                                                          */
/*  void CCelistw_free (CCelistw *elist)                                    */
/*     free a CCelistw                                                      */
/*                                                                          */
/*  void CCelistlw_free (CCelistlw *elist)                                  */
/*     free a CCelistlw                                                     */
/*                                                                          */
/*  int CCelist_alloc (CCelist *elist, int ecount)                          */
/*     allocate space for a CCelist                                         */
/*                                                                          */
/*  int CCelistl_alloc (CCelistl *elist, int ecount)                        */
/*     allocate space for a CCelistl                                        */
/*                                                                          */
/*  int CCelistw_alloc (CCelistw *elist, int ecount)                        */
/*     allocate space for a CCelistw                                        */
/*                                                                          */
/*  int CCelistlw_alloc (CCelistlw *elist, int ecount)                      */
/*     CCelistlw                                                            */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

void CCelist_init (CCelist *elist)
{
    elist->ecount = 0;
    elist->ends = (int *) NULL;
}

void CCelistl_init (CCelistl *elist)
{
    elist->ecount = 0;
    elist->ends = (int *) NULL;
    elist->len = (int *) NULL;
}

void CCelistw_init (CCelistw *elist)
{
    elist->ecount = 0;
    elist->ends = (int *) NULL;
    elist->weight = (double *) NULL;
}

void CCelistlw_init (CCelistlw *elist)
{
    elist->ecount = 0;
    elist->ends = (int *) NULL;
    elist->len = (int *) NULL;
    elist->weight = (double *) NULL;
}

void CCelist_free (CCelist *elist)
{
    elist->ecount = 0;
    CC_IFFREE (elist->ends, int);
}

void CCelistl_free (CCelistl *elist)
{
    elist->ecount = 0;
    CC_IFFREE (elist->ends, int);
    CC_IFFREE (elist->len, int);
}

void CCelistw_free (CCelistw *elist)
{
    elist->ecount = 0;
    CC_IFFREE (elist->ends, int);
    CC_IFFREE (elist->weight, double);
}

void CCelistlw_free (CCelistlw *elist)
{
    elist->ecount = 0;
    CC_IFFREE (elist->ends, int);
    CC_IFFREE (elist->len, int);
    CC_IFFREE (elist->weight, double);
}

int CCelist_alloc (CCelist *elist, int ecount)
{
    elist->ends = CC_SAFE_MALLOC (ecount*2, int);
    if (elist->ends == (int *) NULL) {
        CCelist_free (elist);
        return 1;
    }
    elist->ecount = ecount;
    return 0;
}

int CCelistl_alloc (CCelistl *elist, int ecount)
{
    elist->ends = CC_SAFE_MALLOC (ecount*2, int);
    elist->len = CC_SAFE_MALLOC (ecount, int);
    if (elist->ends == (int *) NULL ||
        elist->len == (int *) NULL) {
        CCelistl_free (elist);
        return 1;
    }
    elist->ecount = ecount;
    return 0;
}

int CCelistw_alloc (CCelistw *elist, int ecount)
{
    elist->ends = CC_SAFE_MALLOC (ecount*2, int);
    elist->weight = CC_SAFE_MALLOC (ecount, double);
    if (elist->ends == (int *) NULL ||
        elist->weight == (double *) NULL) {
        CCelistw_free (elist);
        return 1;
    }
    elist->ecount = ecount;
    return 0;
}

int CCelistlw_alloc (CCelistlw *elist, int ecount)
{
    elist->ends = CC_SAFE_MALLOC (ecount*2, int);
    elist->len = CC_SAFE_MALLOC (ecount, int);
    elist->weight = CC_SAFE_MALLOC (ecount, double);
    if (elist->ends == (int *) NULL ||
        elist->len == (int *) NULL ||
        elist->weight == (double *) NULL) {
        CCelistlw_free (elist);
        return 1;
    }
    elist->ecount = ecount;
    return 0;
}

int CCutil_edge_to_cycle (int ncount, int *elist, int *yesno, int *cyc)
{
    int *Lside, *Rside;
    int i, k, end1, end2, prev, this, next, start, okfirst, first = 0;
    int rval = 0;

    *yesno = 0;

    Lside = CC_SAFE_MALLOC (ncount, int);
    Rside = CC_SAFE_MALLOC (ncount, int);
    if (!Lside || !Rside) {
        fprintf (stderr, "out of memory in CCutil_edge_to_cycle\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) {
        Lside[i] = Rside[i] = -1;
    }

    for (i = 0, k = 0; i < ncount; i++) {
        end1 = elist[k++];
        end2 = elist[k++];
        if (Lside[end1] == -1)
            Lside[end1] = end2;
        else
            Rside[end1] = end2;
        if (Lside[end2] == -1)
            Lside[end2] = end1;
        else
            Rside[end2] = end1;
    }

    for (i = 0, k = 0; i < ncount; i++) {
        end1 = elist[k++];
        end2 = elist[k++];
        if (Lside[end1] == -1 || Rside[end1] == -1 ||
            Lside[end2] == -1 || Rside[end2] == -1) {
            *yesno = 0;  goto CLEANUP;
        }
    }
    start = elist[0];
    prev = -2;
    this = start;
    k = 0;
    okfirst = 0;
    do {
        if (this == first)
           okfirst = 1;
        if (Lside[this] != prev)
            next = Lside[this];
        else
            next = Rside[this];
        prev = this;
        this = next;
        k++;
    } while (next != start && k < ncount);

    if (k != ncount || !okfirst) {
        *yesno = 0;  goto CLEANUP;
    }

    *yesno = 1;

    if (cyc) {
        start = first;
        prev = -2;
        this = start;
        k = 0;
        do {
            cyc[k++] = this;
            if (Lside[this] != prev)
                next = Lside[this];
            else
                next = Rside[this];
            prev = this;
            this = next;
        } while (next != start && k < ncount);
    }


CLEANUP:

    CC_IFFREE (Lside, int);
    CC_IFFREE (Rside, int);

    return rval;
}

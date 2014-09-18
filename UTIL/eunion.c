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
/*                   Find the Union of a Set of Edge Set                    */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*  Written by: Applegate, Bixby, Chvatal, and Cook                         */
/*  Date: May 3, 1999                                                       */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_edge_file_union (int ncount, int nfiles, char **flist,       */
/*      int *ecount, int **elist, int **elen, int *foundtour,               */
/*      int *besttourlen)                                                   */
/*    MERGES a list of edge sets.                                           */
/*     -ncount is the number of nodes.                                      */
/*     -nfiles is the number of files to be read.                           */
/*     -flist is the list of the files.                                     */
/*     -ecount, elist, elen returns the merged edge set.                    */
/*     -foundtour will return a 1 if at least one of the files is the       */
/*      edgeset of a tour (it can be NULL).                                 */
/*     -besttour will return the length of the best tour amongst the        */
/*      edge sets (it can be NULL)                                          */
/*     Returns a nonzero value if there was and error.                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

int CCutil_edge_file_union (int ncount, int nfiles, char **flist, int *ecount,
        int **elist, int **elen, int *foundtour, double *besttourlen) 
{
    int i, j, rval = 0;
    CCutil_edgehash h;
    int tcount;
    int *tlist = (int *) NULL;
    int *telen = (int *) NULL;

    *ecount = 0;
    *elist = (int *) NULL;
    *elen  = (int *) NULL;
    if (foundtour) *foundtour = 0;
    if (besttourlen) *besttourlen = CCutil_MAXDOUBLE;

    rval = CCutil_edgehash_init (&h, 2*ncount);
    if (rval) {
        fprintf (stderr, "CCutil_edgehash_init failed\n"); goto CLEANUP;
    }

    for (i = 0; i < nfiles; i++) {
        rval = CCutil_getedgelist (ncount, flist[i], &tcount, &tlist,
                                   &telen, 0);
        if (rval) {
            fprintf (stderr, "CCutil_getedgelist failed\n"); goto CLEANUP;
        }

        for (j = 0; j < tcount; j++) {
            rval = CCutil_edgehash_set (&h, tlist[2*j], tlist[2*j+1],
                                        telen[j]);
            if (rval) {
                fprintf (stderr, "CCutil_edgehash_set failed\n"); 
                goto CLEANUP;
            }
        }

        if (foundtour && (tcount == ncount)) {
            int yesno;
            rval = CCutil_edge_to_cycle (ncount, tlist, &yesno, (int *) NULL);
            if (rval) {
                fprintf (stderr, "CCutil_edge_to_cycle failed\n"); goto CLEANUP;
            }
            if (yesno) {
                *foundtour = 1;
                if (besttourlen) {
                    double len = 0.0;
                    for (j = 0; j < tcount; j++) {
                        len += (double) telen[j];
                    }
                    if (len < *besttourlen) {
                        *besttourlen = len;
                    }
                }
            }

        }

        CC_IFFREE (tlist, int);
        CC_IFFREE (telen, int);
    }

    rval = CCutil_edgehash_getall (&h, ecount, elist, elen);
    if (rval) {
        fprintf (stderr, "CCutil_edgehash_getall failed\n");
        goto CLEANUP;
    }


CLEANUP:

    if (rval) {
        *ecount = 0;
        CC_IFFREE (*elist, int);
        CC_IFFREE (*elen, int);
        if (foundtour) *foundtour = 0;
        if (besttourlen) *besttourlen = CCutil_MAXDOUBLE;
    }
    CC_IFFREE (tlist, int);
    CC_IFFREE (telen, int);
    CCutil_edgehash_free (&h);

    return rval;
}

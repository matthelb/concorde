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
/*              NEAREST NEIGHBORS FOR X-NORMS AND JUNK-NORMS                */
/*                                                                          */
/*                               TSP CODE                                   */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: March 2, 1995                                                     */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCedgegen_x_k_nearest (int ncount, int k, CCdatagroup *dat,         */
/*      double *wcoord, int wantlist, int *ecount, int **elist,             */
/*      int silent)                                                         */
/*    RETURNS the k_nearest neighbor graph (for X-Norms)                    */
/*      -ncount is the number of nodes                                      */
/*      -k is the number of nearest neighbors wanted                        */
/*      -dat contains the info to generate edge lengths                     */
/*      -wcoord are nodeweights for Held-Karp style edge lengths, using     */
/*       len[i,j] + wcoord[i] + wcoord[j] (wcoord can be NULL)              */
/*      -wantlist should be set to 0 if you don't want the edges            */
/*      -ecount returns the number of edges if wantlist is 1                */
/*      -elist returns the edges in end1 end2 format if wantlist is 1       */
/*                                                                          */
/*  int CCedgegen_x_quadrant_k_nearest (int ncount, int k,                  */
/*      CCdatagroup *dat, double *wcoord, int wantlist, int *ecount,        */
/*      int **elist, int silent)                                            */
/*    RETURNS the quadrant k_nearest_graph (for X-Norms)                    */
/*                                                                          */
/*  int CCedgegen_x_node_k_nearest (CCxnear *xn, int n, int k,              */
/*      int ncount, int *list)                                              */
/*    RETURNS the k nearest neighbors from node n (for X-Norms              */
/*      -xn is a structure built by a call to CCedgegen_xnear_build ()      */
/*      -list returns the neighbors of n. The calling routine should        */
/*       be sure that list points to an array of length at least num.       */
/*                                                                          */
/*  int CCedgegen_x_node_quadrant_k_nearest (CCxnear *xn, int n, int k,     */
/*      int ncount, int *list)                                              */
/*    RETURNS the quadrant k nearest to node n (for X-Norms)                */
/*      -xn is a structure built by a call to CCedgegen_xnear_build ()      */
/*      -list returns the neighbors of n. The calling routine should        */
/*       be sure that list points to a sufficiently large array (4*num      */
/*       for D2_SIZE norms and 8*num for D3_SIZE norms)                     */
/*                                                                          */
/*  int CCedgegen_x_node_nearest (CCxnear *xn, int ncount, int ni,          */
/*      char *marks)                                                        */
/*    RETURNS the nearest unmarked node to node n (as the return value)     */
/*      -marks is an array. The entries that are nonzero correspond to      */
/*       nodes that will not be looked at in the search.                    */
/*                                                                          */
/*  int CCedgegen_x_nearest_neighbor_tour (int ncount, int start,           */
/*      CCdatagroup *dat, int *outcycle, double *val)                       */
/*    RETURNS a nearest neighbor tour, starting at node start.              */
/*      -outcycle will contain the tour if it is not NULL (the calling      */
/*       routine should be sure it points to an array of length at          */
/*       least ncount if it is not set to NULL)                             */
/*      -val will return the length of the tour.                            */
/*                                                                          */
/*  int CCedgegen_junk_k_nearest (int ncount, int k, CCdatagroup *dat,      */
/*      double *wcoord, int wantlist, int *ecount, int **elist,             */
/*      int silent)                                                         */
/*    RETURNS the k-nearest graph (for JUNK-Norms)                          */
/*      -see CCedgegen_x_k_nearest (above) for the variables                */
/*                                                                          */
/*  int CCedgegen_junk_node_k_nearest (CCdatagroup *dat, double *wcoord,    */
/*      int n, int k, int ncount, int *list)                                */
/*    RETURNS the k nearest neighbors to node n (for JUNK-Norms)            */
/*      -list returns the neighbors of n. The calling routine should        */
/*       be sure that list points to an array of length at least num.       */
/*                                                                          */
/*  int CCedgegen_junk_node_nearest (CCdatagroup *dat, double *wcoord,      */
/*      int ncount, int n, char *marks)                                     */
/*    RETURNS the nearest unmarked node to node n (as the return value)     */
/*      -marks is an array, the nodes with marks[i] nonzero are ignored.    */
/*                                                                          */
/*  int CCedgegen_junk_nearest_neighbor_tour (int ncount, int start,        */
/*      CCdatagroup *dat, int *outcycle, double *val, int silent)           */
/*    RETURNS a nearest neighbor tour starting at node start. Note that     */
/*      this will be slow for large problems (it is a quadratic routine)    */
/*      -see the describtion of CCedgegen_x_nearest_neighbor_tour above     */
/*                                                                          */
/*  int CCedgegen_junk_greedy_tour (int ncount, CCdatagroup *dat,           */
/*      int *outcycle, double *val, int ecount, int *elist, int silent)     */
/*     RETURNS a greedy tour using edges from elist.  Disjoint segments     */
/*       from elist will be connected by nearest neighbors (which are slow  */
/*       for large problems).                                               */
/*                                                                          */
/*  int CCedgegen_junk_qboruvka_tour (int ncount, CCdatagroup *dat,         */
/*      int *outcycle, double *val, int ecount, int *elist, int silent)     */
/*     RETURNS a quick-boruvka tour using edges from elist.  Disjoint       */
/*       segments from elist will be connected by nearest neighbors (which  */
/*       are slow for large problems).  See the description of              */
/*       CCkdtree_qboruvka_tour for a description of quick-boruvka.         */
/*                                                                          */
/*  int CCedgegen_xnear_build (int ncount, CCdatagroup *dat,                */
/*      double *wcoord, CCxnear *xn)                                        */
/*    RETURNS the structure needed for calls to                             */
/*      CCedgegen_x_node_k_nearest and                                      */
/*      CCedgegen_x_node_quadrant_k_nearest (the calling routine should     */
/*      be sure that xn points to such a structure). All this routine       */
/*      does is permute the data so that the x coordinates are in           */
/*      nonincreasing order.                                                */
/*                                                                          */
/*  void CCedgegen_xnear_free (CCxnear *xn)                                 */
/*        FREES the CCxnear structure pointed to by xn                      */
/*                                                                          */
/*    NOTES:                                                                */
/*         All routines other than CCedgegen_xnear_free and                 */
/*       CCedgegen_x_node_nearest return 0 on success and 1 on failure      */
/*       (normally due to running out of memory).                           */
/*         The X-Norm functions will also work for KD-Norms, but they are   */
/*       much slower than the KD-Norm functions.                            */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "edgegen.h"
#include "util.h"
#include "macrorus.h"

#define BIGDOUBLE (1e30)
#define NEAR_HEAP_CUTOFF 100  /* When to switch from list to heap       */

#define dtrunc(x) (((x)>0.0)?floor(x):ceil(x))

typedef struct shortedge {
    double length;
    int end;
} shortedge;

typedef struct intptr {
    int this;
    struct intptr *next;
} intptr;

typedef struct tabledat {
    intptr **table;
    CCptrworld intptr_world;
} tabledat;

static void
    add_to_list_and_reset (int *list, int *lcount, shortedge *nearlist,
        int nearnum, int *nodenames),
    insert (int n, int m, shortedge *nearlist, CCdatagroup *dat,
        double *wcoord),
    x_quicksort (int *list, double *x, int l, int u);
static int
    run_x_k_nearest (int ncount, int num, CCdatagroup *dat, double *wcoord,
        int wantlist, int *ecount, int **elist, int doquad, int silent),
    put_in_table (tabledat *td, int i, int j, int *added);



CC_PTRWORLD_LIST_ROUTINES (intptr, int, intptralloc, intptr_bulk_alloc,
        intptrfree, intptr_listadd, intptr_listfree)
CC_PTRWORLD_LEAKS_ROUTINE (intptr, intptr_check_leaks, this, int)


int CCedgegen_x_k_nearest (int ncount, int num, CCdatagroup *dat,
        double *wcoord, int wantlist, int *ecount, int **elist, int silent)
{
    return run_x_k_nearest (ncount, num, dat, wcoord, wantlist, ecount,
                            elist, 0, silent);
}

int CCedgegen_x_quadrant_k_nearest (int ncount, int num, CCdatagroup *dat,
        double *wcoord, int wantlist, int *ecount, int **elist, int silent)
{
    return run_x_k_nearest (ncount, num, dat, wcoord, wantlist, ecount,
                            elist, 1, silent);
}

int CCedgegen_junk_k_nearest (int ncount, int num, CCdatagroup *dat,
        double *wcoord, int wantlist, int *ecount, int **elist, int silent)
{
    return run_x_k_nearest (ncount, num, dat, wcoord, wantlist, ecount,
                            elist, 0, silent);
}


static int run_x_k_nearest (int ncount, int num, CCdatagroup *dat,
        double *wcoord, int wantlist, int *ecount, int **elist, int doquad,
        int silent)
{
    int rval = 0;
    int i, n;
    intptr *ip, *ipnext;
    int total, onlist;
    int added, ntotal = 0;
    int *list = (int *) NULL;
    int goal, usex;
    CCxnear xn;
    int norm;
    tabledat td;

    td.table = (intptr **) NULL;
    CCptrworld_init (&td.intptr_world);
    
    xn.nodenames = (int *) NULL;
    xn.invnames = (int *) NULL;
    xn.w = (double *) NULL;
    CCutil_init_datagroup (&xn.dat);

    if (wantlist) {
        *ecount = 0;
        *elist = (int *) NULL;
    }

    CCutil_dat_getnorm (dat, &norm);
    if ((norm & CC_NORM_BITS) == CC_X_NORM_TYPE ||
        (norm & CC_NORM_BITS) == CC_KD_NORM_TYPE) {
        usex = 1;
        if (!silent) {
            printf ("Using x-norm nearest code\n"); fflush (stdout);
        }
    } else {
        usex = 0;
        if (!silent) {
            printf ("Using junk-norm nearest code\n"); fflush (stdout);
        }
    }

    if (wcoord != (double *) NULL) {
        for (i = 0; i < ncount; i++) {
            if (wcoord[i] < -0.00000001) {
                fprintf (stderr, "Cannot use CCxnear with negative weights\n");
                return 1;
            }
        }
    }

    if (doquad && dat->z != (double *) NULL)
        goal = 8 * num;
    else if (doquad)
        goal = 4 * num;
    else
        goal = num;

    if (usex) {
        if (CCedgegen_xnear_build (ncount, dat, wcoord, &xn)) {
            fprintf (stderr, "build_nodes failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }

    td.table = CC_SAFE_MALLOC (ncount, intptr *);
    if (!td.table) {
        rval = 1;
        goto CLEANUP;
    }
    for (i = 0; i < ncount; i++)
        td.table[i] = (intptr *) NULL;
    list = CC_SAFE_MALLOC (goal, int);
    if (!list) {
        rval = 1;
        goto CLEANUP;
    }

    if (!usex && doquad) {
        printf ("NOTE: Cannot run quadrant nearest with a JUNK norm.\n");
        printf ("      Running nearest instead.\n");
        fflush (stdout);
    }

    for (n = 0; n < ncount; n++) {
        if (usex) {
            if (doquad) {
                if (CCedgegen_x_node_quadrant_k_nearest (&xn, n, num, ncount,
                                                         list)) {
                    fprintf (stderr, "CCedgegen_x_node_quadrant_k_nearest failed\n");
                    rval = 1;
                    goto CLEANUP;
                }
            } else {
                if (CCedgegen_x_node_k_nearest (&xn, n, num, ncount, list)) {
                    fprintf (stderr, "CCedgegen_x_node_k_nearest failed\n");
                    rval = 1;
                    goto CLEANUP;
                }
            }
        } else {
            if (CCedgegen_junk_node_k_nearest (dat, wcoord, n, num, ncount,
                                               list)) {
                fprintf (stderr, "junk_node_k_nearest_failed\n");
                rval = 1;
                goto CLEANUP;
            }
        }
        for (i = 0; i < goal; i++) {
            if (list[i] != -1) {
                if (put_in_table (&td, n, list[i], &added))  {
                    fprintf (stderr, "put_in_table failed\n");
                    rval = 1;
                    goto CLEANUP;
                } else {
                    ntotal += added;
                }
            }
        }
        if (!silent) {
            if (n % 1000 == 999) {
                printf (".");
                fflush (stdout);
            }
        }
    }
  
    if (!silent) {
        printf (" %d edges\n", ntotal); fflush (stdout);
    }

    if (wantlist) {
        int j = 0;
        *elist = CC_SAFE_MALLOC (2 * ntotal, int);
        if (!(*elist)) {
            rval = 1;
            goto CLEANUP;
        }
        *ecount = ntotal;
        for (i = 0; i < ncount; i++) {
            for (ip = td.table[i]; ip; ip = ipnext) {
                ipnext =  ip->next;
                (*elist)[j++] = i;
                (*elist)[j++] = ip->this;
                intptrfree (&td.intptr_world, ip);
            }
            td.table[i] = (intptr *) NULL;
        }
    } else {
        for (i = 0; i < ncount; i++) {
            intptr_listfree (&td.intptr_world, td.table[i]);
            td.table[i] = (intptr *) NULL;
        }
    }

    if (intptr_check_leaks (&td.intptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding intptrs in kdnear\n",
                 total - onlist);
    }

CLEANUP:

    CCptrworld_delete (&td.intptr_world);
    CC_IFFREE (list, int);
    CC_IFFREE (td.table, intptr *);
    if (usex)
        CCedgegen_xnear_free (&xn);

    return rval;
}

int CCedgegen_x_node_quadrant_k_nearest (CCxnear *xn, int ni, int nearnum,
        int ncount, int *list)
{
    int i, j, ntotal = 0;
    shortedge *nearlist = (shortedge *) NULL;
    double scale;
    int goal = (xn->dat.z == (double *) NULL ? 4 * nearnum : 8 * nearnum);
    int n = xn->invnames[ni];
    int norm;

    nearlist = CC_SAFE_MALLOC (nearnum + 1, shortedge);
    if (!nearlist)
        return 1;
    for (i = 0; i < nearnum; i++)
        nearlist[i].length = BIGDOUBLE;
    nearlist[nearnum].length = -BIGDOUBLE;

    CCutil_dat_getnorm (&xn->dat, &norm);
    if (norm == CC_GEOGRAPHIC) scale = CC_GEOGRAPHIC_SCALE;
    else if (norm == CC_GEOM)  scale = CC_GEOM_SCALE;
    else if (norm == CC_ATT)   scale = CC_ATT_SCALE;
    else                       scale = 1.0;
    if ((norm & CC_NORM_SIZE_BITS) == CC_D3_NORM_SIZE) {
        double ny = xn->dat.y[n];
        double nz = xn->dat.z[n];
        for (j = n - 1; j >= 0 && dtrunc((xn->dat.x[n] - xn->dat.x[j]) * scale)
                    < nearlist[0].length; --j)  {
            if (xn->dat.y[j] <= ny && xn->dat.z[j] <= nz)
                insert (n, j, nearlist, &(xn->dat), xn->w);
        }
        add_to_list_and_reset (list, &ntotal, nearlist, nearnum,
                               xn->nodenames);
        for (j = n - 1; j >= 0 && dtrunc((xn->dat.x[n] - xn->dat.x[j]) * scale)
                    < nearlist[0].length; --j)  {
            if (xn->dat.y[j] <= ny && xn->dat.z[j] >= nz)
                insert (n, j, nearlist, &(xn->dat), xn->w);
        }
        add_to_list_and_reset (list, &ntotal, nearlist, nearnum,
                               xn->nodenames);
        for (j = n - 1; j >= 0 && dtrunc((xn->dat.x[n] - xn->dat.x[j]) * scale)
                    < nearlist[0].length; --j)  {
            if (xn->dat.y[j] >= ny && xn->dat.z[j] <= nz)
                insert (n, j, nearlist, &(xn->dat), xn->w);
        }
        add_to_list_and_reset (list, &ntotal, nearlist, nearnum,
                               xn->nodenames);
        for (j = n - 1; j >= 0 && dtrunc((xn->dat.x[n] - xn->dat.x[j]) * scale)
                    < nearlist[0].length; --j)  {
            if (xn->dat.y[j] >= ny && xn->dat.z[j] >= nz)
                insert (n, j, nearlist, &(xn->dat), xn->w);
        }
        add_to_list_and_reset (list, &ntotal, nearlist, nearnum,
                               xn->nodenames);
        for (j = n + 1; j < ncount &&
                        dtrunc((xn->dat.x[j] - xn->dat.x[n]) * scale)
                        < nearlist[0].length; j++) {
            if (xn->dat.y[j] <= ny && xn->dat.z[j] <= nz)
                insert (n, j, nearlist, &(xn->dat), xn->w);
        }
        add_to_list_and_reset (list, &ntotal, nearlist, nearnum,
                               xn->nodenames);
        for (j = n + 1; j < ncount &&
                    dtrunc((xn->dat.x[j] - xn->dat.x[n]) * scale)
                    < nearlist[0].length; j++) {
            if (xn->dat.y[j] <= ny && xn->dat.z[j] >= nz)
                insert (n, j, nearlist, &(xn->dat), xn->w);
        }
        add_to_list_and_reset (list, &ntotal, nearlist, nearnum,
                               xn->nodenames);
        for (j = n + 1; j < ncount &&
                    dtrunc((xn->dat.x[j] - xn->dat.x[n]) * scale)
                    < nearlist[0].length; j++) {
            if (xn->dat.y[j] >= ny && xn->dat.z[j] <= nz)
                insert (n, j, nearlist, &(xn->dat), xn->w);
        }
        add_to_list_and_reset (list, &ntotal, nearlist, nearnum,
                               xn->nodenames);
        for (j = n + 1; j < ncount &&
                    dtrunc((xn->dat.x[j] - xn->dat.x[n]) * scale)
                    < nearlist[0].length; j++) {
            if (xn->dat.y[j] >= ny && xn->dat.z[j] >= nz)
                insert (n, j, nearlist, &(xn->dat), xn->w);
        }
        add_to_list_and_reset (list, &ntotal, nearlist, nearnum,
                               xn->nodenames);
    } else {
        double ny = xn->dat.y[n];
        for (j = n - 1; j >= 0 &&
                    dtrunc((xn->dat.x[n] - xn->dat.x[j]) * scale)
                    < nearlist[0].length; --j)  {
            if (xn->dat.y[j] <= ny)
                insert (n, j, nearlist, &(xn->dat), xn->w);
        }
        add_to_list_and_reset (list, &ntotal, nearlist, nearnum,
                               xn->nodenames);
            for (j = n - 1; j >= 0 &&
                    dtrunc((xn->dat.x[n] - xn->dat.x[j]) * scale)
                    < nearlist[0].length; --j)  {
            if (xn->dat.y[j] >= ny)
                insert (n, j, nearlist, &(xn->dat), xn->w);
        }
        add_to_list_and_reset (list, &ntotal, nearlist, nearnum,
                               xn->nodenames);
        for (j = n + 1; j < ncount &&
                    dtrunc((xn->dat.x[j] - xn->dat.x[n]) * scale)
                    < nearlist[0].length; j++) {
            if (xn->dat.y[j] <= ny)
                insert (n, j, nearlist, &(xn->dat), xn->w);
        }
        add_to_list_and_reset (list, &ntotal, nearlist, nearnum,
                               xn->nodenames);
        for (j = n + 1; j < ncount &&
                    dtrunc((xn->dat.x[j] - xn->dat.x[n]) * scale)
                    < nearlist[0].length; j++) {
            if (xn->dat.y[j] >= ny)
                insert (n, j, nearlist, &(xn->dat), xn->w);
        }
        add_to_list_and_reset (list, &ntotal, nearlist, nearnum,
                               xn->nodenames);
    }

    if (ntotal < goal) {
        for (i = ntotal; i < goal; i++)
            list[i] = -1;
    }
    CC_IFFREE (nearlist, shortedge);
    return 0;
}

static void add_to_list_and_reset (int *list, int *lcount, shortedge *nearlist,
                                   int nearnum, int *nodenames)
{
    int i;

    for (i = 0; i < nearnum; i++) {
        if (nearlist[i].length < BIGDOUBLE) {
            list[(*lcount)++] = nodenames[nearlist[i].end];
            nearlist[i].length = BIGDOUBLE;
        }
    }
}

int CCedgegen_x_node_k_nearest (CCxnear *xn, int ni, int nearnum, int ncount,
        int *list)
{
    int i, j, ntotal;
    shortedge *nearlist = (shortedge *) NULL;
    double scale;
    int n = xn->invnames[ni];
    int norm;

    nearlist = CC_SAFE_MALLOC (nearnum + 1, shortedge);
    if (!nearlist)
        return 1;
    for (i = 0; i < nearnum; i++)
        nearlist[i].length = BIGDOUBLE;
    nearlist[nearnum].length = -BIGDOUBLE;

    CCutil_dat_getnorm (&xn->dat, &norm);
    if (norm == CC_GEOGRAPHIC) scale = CC_GEOGRAPHIC_SCALE;
    else if (norm == CC_GEOM)  scale = CC_GEOM;
    else if (norm == CC_ATT)   scale = CC_ATT_SCALE;
    else                       scale = 1.0;
    for (j = n - 1; j >= 0 &&
               dtrunc((xn->dat.x[n] - xn->dat.x[j]) * scale)
               < nearlist[0].length; --j)
        insert (n, j, nearlist, &(xn->dat), xn->w);
    for (j = n + 1; j < ncount &&
               dtrunc((xn->dat.x[j] - xn->dat.x[n]) * scale)
               < nearlist[0].length; j++)
        insert (n, j, nearlist, &(xn->dat), xn->w);

    ntotal = 0;
    for (i = 0; i < nearnum; i++) {
        if (nearlist[i].length < BIGDOUBLE)
            list[ntotal++] = xn->nodenames[nearlist[i].end];
    }
    if (ntotal < nearnum) {
        fprintf (stderr, "WARNING: There do not exist %d neighbors\n",
                 nearnum);
        for (i = ntotal; i < nearnum; i++)
            list[i] = -1;
        return 1;
    }

    CC_IFFREE (nearlist, shortedge);
    return 0;
}

int CCedgegen_x_node_nearest (CCxnear *xn, int ncount, int ni, char *marks)
{
    int n = xn->invnames[ni];
    int j, bestnode = 0;
    double scale, thisdist, bestdist = BIGDOUBLE;
    int norm;

    CCutil_dat_getnorm (&xn->dat, &norm);
    if (norm == CC_GEOGRAPHIC) scale = CC_GEOGRAPHIC_SCALE;
    else if (norm == CC_GEOM)  scale = CC_GEOM;
    else if (norm == CC_ATT)   scale = CC_ATT_SCALE;
    else                       scale = 1.0;
    for (j = n - 1; j >= 0 &&
               dtrunc((xn->dat.x[n] - xn->dat.x[j]) * scale)
               < bestdist; --j) {
        if (!marks[xn->nodenames[j]]) {
            thisdist = CCutil_dat_edgelen (n, j, &(xn->dat));
            if (xn->w)
                thisdist += (xn->w[n] + xn->w[j]);
            if (thisdist < bestdist) {
                bestdist = thisdist;
                bestnode = j;
            }
        }
    }
    for (j = n + 1; j < ncount &&
               dtrunc((xn->dat.x[j] - xn->dat.x[n]) * scale)
               < bestdist; j++) {
        if (!marks[xn->nodenames[j]]) {
            thisdist = CCutil_dat_edgelen (n, j, &(xn->dat));
            if (xn->w)
                thisdist += (xn->w[n] + xn->w[j]);
            if (thisdist < bestdist) {
                bestdist = thisdist;
                bestnode = j;
            }
        }
    }
    return xn->nodenames[bestnode];
}

static void insert (int n, int m, shortedge *nearlist, CCdatagroup *dat,
                    double *wcoord)
{
    int i;
    int thisdist;

    thisdist = CCutil_dat_edgelen (n, m, dat);

    if (wcoord != (double *) NULL)
        thisdist += (wcoord[n] + wcoord[m]);

    if (thisdist < nearlist[0].length) {
        for (i = 0; nearlist[i+1].length > thisdist; i++) {
            nearlist[i].end = nearlist[i + 1].end;
            nearlist[i].length = nearlist[i + 1].length;
        }
        nearlist[i].length = thisdist;
        nearlist[i].end = m;
    }
}

int CCedgegen_junk_node_k_nearest (CCdatagroup *dat, double *wcoord, int n,
        int nearnum, int ncount, int *list)
{
    int i, j, ntotal;
    shortedge *nearlist = (shortedge *) NULL;

    nearlist = CC_SAFE_MALLOC (nearnum + 1, shortedge);
    if (!nearlist)
        return 1;
    for (i = 0; i < nearnum; i++)
        nearlist[i].length = BIGDOUBLE;
    nearlist[nearnum].length = -BIGDOUBLE;

    for (j = n - 1; j >= 0; j--) {
        insert (n, j, nearlist, dat, wcoord);
    }
    for (j = n + 1; j < ncount; j++) {
        insert (n, j, nearlist, dat, wcoord);
    }

    ntotal = 0;
    for (i = 0; i < nearnum; i++) {
        if (nearlist[i].length < BIGDOUBLE)
            list[ntotal++] = nearlist[i].end;
    }
    if (ntotal < nearnum) {
        fprintf (stderr, "WARNING: There do not exist %d neighbors\n",
                 nearnum);
        for (i = ntotal; i < nearnum; i++)
            list[i] = -1;
        return 1;
    }

    CC_IFFREE (nearlist, shortedge);
    return 0;
}

int CCedgegen_junk_node_nearest (CCdatagroup *dat, double *wcoord, int ncount,
        int n, char *marks)
{
    int j, bestnode = 0;
    double thisdist, bestdist = BIGDOUBLE;


    if (wcoord) {
        for (j = n - 1; j >= 0; j--) {
            if (!marks[j]) {
                thisdist = CCutil_dat_edgelen (n, j, dat) +
                           (wcoord[n] + wcoord[j]);
                if (thisdist < bestdist) {
                    bestdist = thisdist;
                    bestnode = j;
                }
            }
        }
        for (j = n + 1; j < ncount; j++) {
            if (!marks[j]) {
                thisdist = CCutil_dat_edgelen (n, j, dat) +
                           (wcoord[n] + wcoord[j]);
                if (thisdist < bestdist) {
                    bestdist = thisdist;
                    bestnode = j;
                }
            }
        }
    } else {
        for (j = n - 1; j >= 0; j--) {
            if (!marks[j]) {
                thisdist = CCutil_dat_edgelen (n, j, dat);
                if (thisdist < bestdist) {
                    bestdist = thisdist;
                    bestnode = j;
                }
            }
        }
        for (j = n + 1; j < ncount; j++) {
            if (!marks[j]) {
                thisdist = CCutil_dat_edgelen (n, j, dat);
                if (thisdist < bestdist) {
                    bestdist = thisdist;
                    bestnode = j;
                }
            }
        }
    }

    return bestnode;
}

int CCedgegen_x_nearest_neighbor_tour (int ncount, int start, CCdatagroup *dat,
                             int *outcycle, double *val)
{
    double len;
    int i, current, next;
    CCxnear xn;
    char *marks;
    int norm;

    /*
        printf ("Grow nearest neighbor tour from node %d\n", start);
        fflush (stdout);
    */

    if (ncount < 3) {
        fprintf (stderr, "Cannot find tour in an %d node graph\n", ncount);
        return 1;
    }
    CCutil_dat_getnorm (dat, &norm);
    if ((norm & CC_NORM_BITS) != CC_X_NORM_TYPE &&
        (norm & CC_NORM_BITS) != CC_KD_NORM_TYPE) {
        fprintf (stderr, "Cannot run x_nearest with norm %d\n", norm);
        return 1;
    }

    if (CCedgegen_xnear_build (ncount, dat, (double *) NULL, &xn)) {
        fprintf (stderr, "Unable to build CCxnear\n");
        return 1;
    }

    marks = CC_SAFE_MALLOC (ncount, char );
    if (!marks) {
        CCedgegen_xnear_free (&xn);
        return 1;
    }

    for (i = 0; i < ncount; i++)
        marks[i] = 0;


    len = 0.0;
    current = start;
    if (outcycle != (int *) NULL)
        outcycle[0] = start;

    for (i = 1; i < ncount; i++) {
        marks[current] = 1;
        next = CCedgegen_x_node_nearest (&xn, ncount, current, marks);
        if (outcycle != (int *) NULL)
            outcycle [i] = next;
        len += (double) CCutil_dat_edgelen (current, next, dat);
        current = next;
    }
    len += (double) CCutil_dat_edgelen (current, start, dat);
    *val = len;
    CCedgegen_xnear_free (&xn);
    CC_IFFREE (marks, char);
    return 0;
}

int CCedgegen_x_greedy_tour (int ncount, CCdatagroup *dat, int *outcycle,
        double *val, int ecount, int *elist, int silent)
{
    int rval;
    int *perm = (int *) NULL;
    int *elen = (int *) NULL;
    int *tcyc = (int *) NULL;
    int *degree = (int *) NULL;
    int *tail = (int *) NULL;
    char *marks = (char *) NULL;
    int norm;
    CCxnear xn;
    int tcount = 0;
    int x;
    int y;
    int z;
    int i;
    double len;

    if (!silent) {
        printf ("Grow a greedy tour \n");
        fflush (stdout);
    }

    CCutil_dat_getnorm (dat, &norm);
    if ((norm & CC_NORM_BITS) != CC_X_NORM_TYPE &&
        (norm & CC_NORM_BITS) != CC_KD_NORM_TYPE) {
        fprintf (stderr, "Cannot run x_nearest with norm %d\n", norm);
        return 1;
    }

    if (CCedgegen_xnear_build (ncount, dat, (double *) NULL, &xn)) {
        fprintf (stderr, "Unable to build CCxnear\n");
        return 1;
    }

    perm = CC_SAFE_MALLOC (ecount, int);
    elen = CC_SAFE_MALLOC (ecount, int);
    if (perm == (int *) NULL ||
        elen == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCedgegen_x_greedy_tour\n");
        rval = 1; goto CLEANUP;
    }
    for (i=0; i<ecount; i++) {
        perm[i] = i;
        elen[i] = CCutil_dat_edgelen (elist[2*i], elist[2*i+1], dat);
    }
    CCutil_int_perm_quicksort (perm, elen, ecount);

    if (outcycle) {
        tcyc = CC_SAFE_MALLOC (2 * ncount, int);
        if (!tcyc) {
            rval = 1;
            goto CLEANUP;
        }
    }

    degree = CC_SAFE_MALLOC (ncount, int);
    if (!degree) {
        rval = 1;
        goto CLEANUP;
    }
    tail = CC_SAFE_MALLOC (ncount, int);
    if (!tail) {
        rval = 1;
        goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) {
        degree[i] = 0;
        tail[i] = -1;
    }

    len = 0.0;
    for (i = 0; i < ecount; i++) {
        x = elist[2*perm[i]];
        y = elist[2*perm[i]+1];
        if (degree[x] < 2 && degree[y] < 2 && y != tail[x]) {
            if (tcyc) {
                tcyc[2*tcount] = x;
                tcyc[2*tcount+1] = y;
            }
            tcount++;
            len += elen[perm[i]];
            (degree[x])++;
            (degree[y])++;
            if (tail[x] == -1) {
                if (tail[y] == -1) {
                    tail[x] = y;
                    tail[y] = x;
                } else {
                    tail[x] = tail[y];
                    tail[tail[y]] = x;
                }
            } else if (tail[y] == -1) {
                tail[tail[x]] = y;
                tail[y] = tail[x];
            } else {
                tail[tail[x]] = tail[y];
                tail[tail[y]] = tail[x];
            }
        }
    }
    CC_FREE (perm, int);
    CC_FREE (elen, int);
    marks = CC_SAFE_MALLOC (ncount, char);
    if (!marks) {
        rval = 1; goto CLEANUP;
    }
    for (i=0; i<ncount; i++) {
        if (degree[i] == 2) marks[i] = 1;
        else marks[i] = 0;
    }
    for (x = 0; marks[x] == 1; x++);
    if (tail[x] != -1) marks[tail[x]] = 1;
    while (tcount < ncount-1) {
        marks[x] = 1;
        y = CCedgegen_x_node_nearest (&xn, ncount, x, marks);
        if (tcyc) {
            tcyc[2*tcount] = x;
            tcyc[2*tcount+1] = y;
        }
        tcount++;
        len += (double) CCutil_dat_edgelen (x, y, dat);
        if (tail[y] == -1) z = y;
        else z = tail[y];
        if (degree[x]) marks[x] = 1;
        (degree[x])++;
        if (degree[y]) marks[y] = 1;
        (degree[y])++;
        if (tail[x] == -1) {
            if (tail[y] == -1) {
                tail[x] = y;
                tail[y] = x;
            } else {
                tail[x] = tail[y];
                tail[tail[y]] = x;
            }
        } else if (tail[y] == -1) {
            tail[tail[x]] = y;
            tail[y] = tail[x];
        } else {
            tail[tail[x]] = tail[y];
            tail[tail[y]] = tail[x];
        }
        x = z;
    }
    for (x = 0; degree[x] != 1; x++);
    for (y = x + 1; degree[y] != 1; y++);
    if (tcyc) {
        tcyc[2*tcount] = x;
        tcyc[2*tcount+1] = y;
    }
    tcount++;
    len += (double) CCutil_dat_edgelen (x, y, dat);
    *val = len;
    if (!silent) {
        printf ("Length of Greedy Tour: %.2f\n", len);
        fflush (stdout);
    }

    if (tcyc) {
        int istour;
        rval = CCutil_edge_to_cycle (ncount, tcyc, &istour, outcycle);
        if (rval) {
            fprintf (stderr, "CCutil_edge_to_cycle failed\n"); goto CLEANUP;
        }
        if (istour == 0) {
            fprintf (stderr, "ERROR: greedy tour is not a tour\n");
            rval = 1; goto CLEANUP;
        }
    }
    rval = 0;

CLEANUP:
    CCedgegen_xnear_free (&xn);
    CC_IFFREE (tcyc, int);
    CC_IFFREE (degree, int);
    CC_IFFREE (tail, int);
    CC_IFFREE (marks, char);
    CC_IFFREE (perm, int);
    CC_IFFREE (elen, int);
    return rval;
}

int CCedgegen_x_qboruvka_tour (int ncount, CCdatagroup *dat, int *outcycle,
        double *val, int ecount, int *elist, int silent)
{
    int rval;
    int *perm = (int *) NULL;
    int *elen = (int *) NULL;
    int *deg  = (int *) NULL;
    int **adj = (int **) NULL;
    int *adjspace = (int *) NULL;
    int *tcyc = (int *) NULL;
    int *degree = (int *) NULL;
    int *tail = (int *) NULL;
    char *marks = (char *) NULL;
    int tcount = 0;
    double len;
    int norm;
    CCxnear xn;
    int x;
    int y = 0;
    int i;
    int j;
    int k;
    int count;
    int found;

    if (!silent) {
        printf ("Grow a Quick-Boruvka tour \n");
        fflush (stdout);
    }

    CCutil_dat_getnorm (dat, &norm);
    if ((norm & CC_NORM_BITS) != CC_X_NORM_TYPE &&
        (norm & CC_NORM_BITS) != CC_KD_NORM_TYPE) {
        fprintf (stderr, "Cannot run x_nearest with norm %d\n", norm);
        return 1;
    }

    if (CCedgegen_xnear_build (ncount, dat, (double *) NULL, &xn)) {
        fprintf (stderr, "Unable to build CCxnear\n");
        return 1;
    }

    perm = CC_SAFE_MALLOC (ecount, int);
    elen = CC_SAFE_MALLOC (ecount, int);
    if (perm == (int *) NULL ||
        elen == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCedgegen_x_greedy_tour\n");
        rval = 1; goto CLEANUP;
    }
    for (i=0; i<ecount; i++) {
        perm[i] = i;
        elen[i] = CCutil_dat_edgelen (elist[2*i], elist[2*i+1], dat);
    }
    CCutil_int_perm_quicksort (perm, elen, ecount);

    CC_FREE (elen, int);
    
    deg  = CC_SAFE_MALLOC (ncount, int);
    adj  = CC_SAFE_MALLOC (ncount, int *);
    adjspace = CC_SAFE_MALLOC (ecount*2, int);
    if (deg  == (int *) NULL ||
        adj  == (int **) NULL ||
        adjspace == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCedgegen_x_greedy_tour\n");
        rval = 1; goto CLEANUP;
    }
    for (i=0; i<ncount; i++) {
        deg[i] = 0;
    }
    
    for (i=0; i<ecount; i++) {
        deg[elist[2*i]]++;
        deg[elist[2*i+1]]++;
    }
    j = 0;
    for (i=0; i<ncount; i++) {
        adj[i] = &(adjspace[j]);
        j += deg[i];
        deg[i] = 0;
    }

    for (i=0; i<ecount; i++) {
        j = elist[2*perm[i]];
        k = elist[2*perm[i]+1];
        adj[j][deg[j]++] = k;
        adj[k][deg[k]++] = j;
    }
    CC_FREE (perm, int);

    if (outcycle) {
        tcyc = CC_SAFE_MALLOC (2 * ncount, int);
        if (!tcyc) {
            rval = 1;
            goto CLEANUP;
        }
    }

    degree = CC_SAFE_MALLOC (ncount, int);
    if (!degree) {
        rval = 1;
        goto CLEANUP;
    }
    tail = CC_SAFE_MALLOC (ncount, int);
    if (!tail) {
        rval = 1;
        goto CLEANUP;
    }
    marks = CC_SAFE_MALLOC (ncount, char);
    if (!marks) {
        rval = 1;
        goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) {
        degree[i] = 0;
        tail[i] = -1;
        marks[i] = 0;
    }

    len = 0.0;
    count = 1;
    while (count < ncount) {
        for (x = 0; x < ncount && count < ncount; x++) {
            if (degree[x] != 2) {
                found = 0;
                for (j=0; j<deg[x] && !found; j++) {
                    y = adj[x][j];
                    if (degree[y] != 2 && y != tail[x]) {
                        found = 1;
                    }
                }
                if (!found) {
                    if (tail[x] != -1) {
                        marks[tail[x]] = 1;
                        y = CCedgegen_x_node_nearest (&xn, ncount, x, marks);
                        marks[tail[x]] = 0;
                    } else {
                        y = CCedgegen_x_node_nearest (&xn, ncount, x, marks);
                    }
                }

                /* add (x, y) to the tour */
                if (degree[x] != 0)
                    marks[x] = 1;
                if (degree[y] != 0)
                    marks[y] = 1;
                len += (double) CCutil_dat_edgelen (x, y, dat);
                degree[x]++;
                degree[y]++;
                if (tcyc) {
                    tcyc[tcount++] = x;
                    tcyc[tcount++] = y;
                }
                if (tail[x] == -1) {
                    if (tail[y] == -1) {
                        tail[x] = y;
                        tail[y] = x;
                    } else {
                        tail[x] = tail[y];
                        tail[tail[y]] = x;
                    }
                } else if (tail[y] == -1) {
                    tail[tail[x]] = y;
                    tail[y] = tail[x];
                } else {
                    tail[tail[x]] = tail[y];
                    tail[tail[y]] = tail[x];
                }
                if (count % 10000 == 9999) {
                    printf (".");
                    fflush (stdout);
                }
                count++;
            }
        }
    }
    for (x = 0; degree[x] != 1; x++);
    for (y = x + 1; degree[y] != 1; y++);
    if (tcyc) {
        tcyc[tcount++] = x;
        tcyc[tcount++] = y;
    }
    len += (double) CCutil_dat_edgelen (x, y, dat);
    *val = len;
    if (ncount >= 10000)
        printf ("\n");
    printf ("Length of Quick-Boruvka Tour: %.2f\n", len);

    if (tcyc) {
        int istour;
        rval = CCutil_edge_to_cycle (ncount, tcyc, &istour, outcycle);
        if (rval) {
            fprintf (stderr, "CCutil_edge_to_cycle failed\n"); goto CLEANUP;
        }
        if (istour == 0) {
            fprintf (stderr, "ERROR: greedy tour is not a tour\n");
            rval = 1; goto CLEANUP;
        }
    }
    rval = 0;

CLEANUP:
    CCedgegen_xnear_free (&xn);
    CC_IFFREE (tcyc, int);
    CC_IFFREE (degree, int);
    CC_IFFREE (tail, int);
    CC_IFFREE (marks, char);
    CC_IFFREE (deg, int);
    CC_IFFREE (adj, int *);
    CC_IFFREE (adjspace, int);
    CC_IFFREE (perm, int);
    CC_IFFREE (elen, int);
    return rval;
}

int CCedgegen_junk_nearest_neighbor_tour (int ncount, int start,
        CCdatagroup *dat, int *outcycle, double *val, int silent)
{
    double len;
    int i, current, next;
    char *marks;

/*
    printf ("Grow nearest neighbor tour from node %d\n", start);
*/
    if (!silent) {
        printf ("This is a JUNK norm, so expect a quadratic running time\n");
        fflush (stdout);
    }

    if (ncount < 3) {
        fprintf (stderr, "Cannot find tour in an %d node graph\n", ncount);
        return 1;
    }

    marks = CC_SAFE_MALLOC (ncount, char );
    if (!marks) {
        return 1;
    }

    for (i = 0; i < ncount; i++)
        marks[i] = 0;

    len = 0.0;
    current = start;
    if (outcycle != (int *) NULL)
        outcycle[0] = start;

    for (i = 1; i < ncount; i++) {
        marks[current] = 1;
        next = CCedgegen_junk_node_nearest (dat, (double *) NULL, ncount,
                                            current, marks);
        if (outcycle != (int *) NULL)
            outcycle [i] = next;
        len += (double) CCutil_dat_edgelen (current, next, dat);
        current = next;
    }
    len += (double) CCutil_dat_edgelen (current, start, dat);
    *val = len;
    CC_IFFREE (marks, char);
    return 0;
}

int CCedgegen_junk_greedy_tour (int ncount, CCdatagroup *dat, int *outcycle,
        double *val, int ecount, int *elist, int silent)
{
    int rval;
    int *perm = (int *) NULL;
    int *elen = (int *) NULL;
    int *tcyc = (int *) NULL;
    int *degree = (int *) NULL;
    int *tail = (int *) NULL;
    char *marks = (char *) NULL;
    int tcount = 0;
    int x;
    int y;
    int z;
    int i;
    double len;

    if (!silent) {
        printf ("Grow a greedy tour \n");
        fflush (stdout);
    }

    perm = CC_SAFE_MALLOC (ecount, int);
    elen = CC_SAFE_MALLOC (ecount, int);
    if (perm == (int *) NULL ||
        elen == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCedgegen_junk_greedy_tour\n");
        rval = 1; goto CLEANUP;
    }
    for (i=0; i<ecount; i++) {
        perm[i] = i;
        elen[i] = CCutil_dat_edgelen (elist[2*i], elist[2*i+1], dat);
    }
    CCutil_int_perm_quicksort (perm, elen, ecount);

    if (outcycle) {
        tcyc = CC_SAFE_MALLOC (2 * ncount, int);
        if (!tcyc) {
            rval = 1;
            goto CLEANUP;
        }
    }

    degree = CC_SAFE_MALLOC (ncount, int);
    if (!degree) {
        rval = 1;
        goto CLEANUP;
    }
    tail = CC_SAFE_MALLOC (ncount, int);
    if (!tail) {
        rval = 1;
        goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) {
        degree[i] = 0;
        tail[i] = -1;
    }

    len = 0.0;
    for (i = 0; i < ecount; i++) {
        x = elist[2*perm[i]];
        y = elist[2*perm[i]+1];
        if (degree[x] < 2 && degree[y] < 2 && y != tail[x]) {
            if (tcyc) {
                tcyc[2*tcount] = x;
                tcyc[2*tcount+1] = y;
            }
            tcount++;
            len += elen[perm[i]];
            (degree[x])++;
            (degree[y])++;
            if (tail[x] == -1) {
                if (tail[y] == -1) {
                    tail[x] = y;
                    tail[y] = x;
                } else {
                    tail[x] = tail[y];
                    tail[tail[y]] = x;
                }
            } else if (tail[y] == -1) {
                tail[tail[x]] = y;
                tail[y] = tail[x];
            } else {
                tail[tail[x]] = tail[y];
                tail[tail[y]] = tail[x];
            }
        }
    }
    CC_FREE (perm, int);
    CC_FREE (elen, int);
    marks = CC_SAFE_MALLOC (ncount, char);
    if (!marks) {
        rval = 1; goto CLEANUP;
    }
    for (i=0; i<ncount; i++) {
        if (degree[i] == 2) marks[i] = 1;
        else marks[i] = 0;
    }
    for (x = 0; marks[x] == 1; x++);
    if (tail[x] != -1) marks[tail[x]] = 1;
    while (tcount < ncount-1) {
        marks[x] = 1;
        y = CCedgegen_junk_node_nearest (dat, (double *) NULL, ncount, x,
                                         marks);
        if (tcyc) {
            tcyc[2*tcount] = x;
            tcyc[2*tcount+1] = y;
        }
        tcount++;
        len += (double) CCutil_dat_edgelen (x, y, dat);
        if (tail[y] == -1) z = y;
        else z = tail[y];
        if (degree[x]) marks[x] = 1;
        (degree[x])++;
        if (degree[y]) marks[y] = 1;
        (degree[y])++;
        if (tail[x] == -1) {
            if (tail[y] == -1) {
                tail[x] = y;
                tail[y] = x;
            } else {
                tail[x] = tail[y];
                tail[tail[y]] = x;
            }
        } else if (tail[y] == -1) {
            tail[tail[x]] = y;
            tail[y] = tail[x];
        } else {
            tail[tail[x]] = tail[y];
            tail[tail[y]] = tail[x];
        }
        x = z;
    }
    for (x = 0; degree[x] != 1; x++);
    for (y = x + 1; degree[y] != 1; y++);
    if (tcyc) {
        tcyc[2*tcount] = x;
        tcyc[2*tcount+1] = y;
    }
    tcount++;
    len += (double) CCutil_dat_edgelen (x, y, dat);
    *val = len;
    if (!silent) {
        printf ("Length of Greedy Tour: %.2f\n", len);
        fflush (stdout);
    }

    if (tcyc) {
        int istour;
        rval = CCutil_edge_to_cycle (ncount, tcyc, &istour, outcycle);
        if (rval) {
            fprintf (stderr, "CCutil_edge_to_cycle failed\n"); goto CLEANUP;
        }
        if (istour == 0) {
            fprintf (stderr, "ERROR: greedy tour is not a tour\n");
            rval = 1; goto CLEANUP;
        }
    }
    rval = 0;

CLEANUP:
    CC_IFFREE (tcyc, int);
    CC_IFFREE (degree, int);
    CC_IFFREE (tail, int);
    CC_IFFREE (marks, char);
    CC_IFFREE (perm, int);
    CC_IFFREE (elen, int);
    return rval;
}

int CCedgegen_junk_qboruvka_tour (int ncount, CCdatagroup *dat, int *outcycle,
        double *val, int ecount, int *elist, int silent)
{
    int rval;
    int *perm = (int *) NULL;
    int *elen = (int *) NULL;
    int *deg  = (int *) NULL;
    int **adj = (int **) NULL;
    int *adjspace = (int *) NULL;
    int *tcyc = (int *) NULL;
    int *degree = (int *) NULL;
    int *tail = (int *) NULL;
    char *marks = (char *) NULL;
    int tcount = 0;
    double len;
    int x;
    int y = 0;
    int i;
    int j;
    int k;
    int count;
    int found;

    if (!silent) {
        printf ("Grow a Quick-Boruvka tour \n");
        fflush (stdout);
    }

    perm = CC_SAFE_MALLOC (ecount, int);
    elen = CC_SAFE_MALLOC (ecount, int);
    if (perm == (int *) NULL ||
        elen == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCedgegen_junk_greedy_tour\n");
        rval = 1; goto CLEANUP;
    }
    for (i=0; i<ecount; i++) {
        perm[i] = i;
        elen[i] = CCutil_dat_edgelen (elist[2*i], elist[2*i+1], dat);
    }
    CCutil_int_perm_quicksort (perm, elen, ecount);

    CC_FREE (elen, int);
    
    deg  = CC_SAFE_MALLOC (ncount, int);
    adj  = CC_SAFE_MALLOC (ncount, int *);
    adjspace = CC_SAFE_MALLOC (ecount*2, int);
    if (deg  == (int *) NULL ||
        adj  == (int **) NULL ||
        adjspace == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCedgegen_junk_greedy_tour\n");
        rval = 1; goto CLEANUP;
    }
    for (i=0; i<ncount; i++) {
        deg[i] = 0;
    }
    
    for (i=0; i<ecount; i++) {
        deg[elist[2*i]]++;
        deg[elist[2*i+1]]++;
    }
    j = 0;
    for (i=0; i<ncount; i++) {
        adj[i] = &(adjspace[j]);
        j += deg[i];
        deg[i] = 0;
    }

    for (i=0; i<ecount; i++) {
        j = elist[2*perm[i]];
        k = elist[2*perm[i]+1];
        adj[j][deg[j]++] = k;
        adj[k][deg[k]++] = j;
    }
    CC_FREE (perm, int);

    if (outcycle) {
        tcyc = CC_SAFE_MALLOC (2 * ncount, int);
        if (!tcyc) {
            rval = 1;
            goto CLEANUP;
        }
    }

    degree = CC_SAFE_MALLOC (ncount, int);
    if (!degree) {
        rval = 1;
        goto CLEANUP;
    }
    tail = CC_SAFE_MALLOC (ncount, int);
    if (!tail) {
        rval = 1;
        goto CLEANUP;
    }
    marks = CC_SAFE_MALLOC (ncount, char);
    if (!marks) {
        rval = 1;
        goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) {
        degree[i] = 0;
        tail[i] = -1;
        marks[i] = 0;
    }

    len = 0.0;
    count = 1;
    while (count < ncount) {
        for (x = 0; x < ncount && count < ncount; x++) {
            if (degree[x] != 2) {
                found = 0;
                for (j=0; j<deg[x] && !found; j++) {
                    y = adj[x][j];
                    if (degree[y] != 2 && y != tail[x]) {
                        found = 1;
                    }
                }
                if (!found) {
                    if (tail[x] != -1) {
                        marks[tail[x]] = 1;
                        y = CCedgegen_junk_node_nearest (dat, (double *) NULL,
                                ncount, x, marks);
                        marks[tail[x]] = 0;
                    } else {
                        y = CCedgegen_junk_node_nearest (dat, (double *) NULL,
                                ncount, x, marks);
                    }
                }

                /* add (x, y) to the tour */
                if (degree[x] != 0)
                    marks[x] = 1;
                if (degree[y] != 0)
                    marks[y] = 1;
                len += (double) CCutil_dat_edgelen (x, y, dat);
                degree[x]++;
                degree[y]++;
                if (tcyc) {
                    tcyc[tcount++] = x;
                    tcyc[tcount++] = y;
                }
                if (tail[x] == -1) {
                    if (tail[y] == -1) {
                        tail[x] = y;
                        tail[y] = x;
                    } else {
                        tail[x] = tail[y];
                        tail[tail[y]] = x;
                    }
                } else if (tail[y] == -1) {
                    tail[tail[x]] = y;
                    tail[y] = tail[x];
                } else {
                    tail[tail[x]] = tail[y];
                    tail[tail[y]] = tail[x];
                }
                if (count % 10000 == 9999) {
                    printf (".");
                    fflush (stdout);
                }
                count++;
            }
        }
    }
    for (x = 0; degree[x] != 1; x++);
    for (y = x + 1; degree[y] != 1; y++);
    if (tcyc) {
        tcyc[tcount++] = x;
        tcyc[tcount++] = y;
    }
    len += (double) CCutil_dat_edgelen (x, y, dat);
    *val = len;
    if (ncount >= 10000)
        printf ("\n");
    printf ("Length of Quick-Boruvka Tour: %.2f\n", len);

    if (tcyc) {
        int istour;
        rval = CCutil_edge_to_cycle (ncount, tcyc, &istour, outcycle);
        if (rval) {
            fprintf (stderr, "CCutil_edge_to_cycle failed\n"); goto CLEANUP;
        }
        if (istour == 0) {
            fprintf (stderr, "ERROR: greedy tour is not a tour\n");
            rval = 1; goto CLEANUP;
        }
    }
    rval = 0;

CLEANUP:
    CC_IFFREE (tcyc, int);
    CC_IFFREE (degree, int);
    CC_IFFREE (tail, int);
    CC_IFFREE (marks, char);
    CC_IFFREE (deg, int);
    CC_IFFREE (adj, int *);
    CC_IFFREE (adjspace, int);
    CC_IFFREE (perm, int);
    CC_IFFREE (elen, int);
    return rval;
}

static int put_in_table (tabledat *td, int i, int j, int *added)
{
    intptr *ip;

    if (j < i) {
        int temp;
        CC_SWAP(i, j, temp);
    }

    for (ip = td->table[i]; ip; ip = ip->next) {
        if (ip->this == j) {
            *added = 0;
            return 0;
        }
    }
    if (intptr_listadd (&td->table[i], j, &td->intptr_world)) {
        *added = 0;
        return 1;
    }
    *added = 1;
    return 0;
}

int CCedgegen_xnear_build (int ncount, CCdatagroup *dat, double *wcoord,
        CCxnear *xn)
{
    int i;
    int norm;

    xn->nodenames = (int *) NULL;
    xn->invnames = (int *) NULL;
    xn->w = (double *) NULL;
    CCutil_init_datagroup (&xn->dat);

    CCutil_dat_getnorm (dat, &norm);
    if (CCutil_dat_setnorm (&xn->dat, norm))
        return 1;

    xn->nodenames = CC_SAFE_MALLOC (ncount, int);
    if (!xn->nodenames)
        return 1;
    for (i = 0; i < ncount; i++)
        xn->nodenames[i] = i;
    xn->dat.x = CC_SAFE_MALLOC (ncount, double);
    if (!xn->dat.x) {
        CC_FREE (xn->nodenames, int);
        return 1;
    }
    for (i = 0; i < ncount; i++)
        xn->dat.x[i] = dat->x[i];

    for (i = 1; i < ncount && dat->x[i] >= dat->x[i - 1]; i++);
    if (i < ncount) {
        x_quicksort (xn->nodenames, xn->dat.x, 0, ncount - 1);
    }

    xn->invnames = CC_SAFE_MALLOC (ncount, int);
    if (!xn->invnames) {
        CC_FREE (xn->nodenames, int);
        CCutil_freedatagroup (&(xn->dat));
        return 1;
    }
    for (i = 0; i < ncount; i++)
        xn->invnames[xn->nodenames[i]] = i;

    xn->dat.y = CC_SAFE_MALLOC (ncount, double);
    if (!xn->dat.y) {
        CC_FREE (xn->nodenames, int);
        CC_FREE (xn->invnames, int);
        CCutil_freedatagroup (&(xn->dat));
        return 1;
    }
    for (i = 0; i < ncount; i++)
        xn->dat.y[i] = dat->y[xn->nodenames[i]];
    if (dat->z != (double *) NULL) {
        xn->dat.z = CC_SAFE_MALLOC (ncount, double);
        if (!xn->dat.z) {
            CC_FREE (xn->nodenames, int);
            CC_FREE (xn->invnames, int);
            CCutil_freedatagroup (&(xn->dat));
            return 1;
        }
        for (i = 0; i < ncount; i++)
            xn->dat.z[i] = dat->z[xn->nodenames[i]];
    }
    if (wcoord != (double *) NULL) {
        xn->w = CC_SAFE_MALLOC (ncount, double);
        if (!xn->w) {
            CC_FREE (xn->nodenames, int);
            CC_FREE (xn->invnames, int);
            CCutil_freedatagroup (&(xn->dat));
            return 1;
        }
        for (i = 0; i < ncount; i++)
            xn->w[i] = wcoord[xn->nodenames[i]];
    }
    return 0;
}

void CCedgegen_xnear_free (CCxnear *xn)
{
    CC_IFFREE (xn->nodenames, int);
    CC_IFFREE (xn->invnames, int);
    CC_IFFREE (xn->w, double);
    CCutil_freedatagroup (&(xn->dat));
}

static void x_quicksort (int *list, double *x, int l, int u)
{
    int i, j, itemp;
    double t, dtemp;

    if (l >= u)
        return;

    CC_SWAP (x[l], x[(l+u)/2], dtemp);
    CC_SWAP (list[l], list[(l+u)/2], itemp);

    i = l;
    j = u + 1;
    t = x[l];

    while (1) {
        do i++; while (i <= u && x[i] < t);
        do j--; while (x[j] > t);
        if (j < i) break;
        CC_SWAP (x[i], x[j], dtemp);
        CC_SWAP (list[i], list[j], itemp);
    }
    CC_SWAP (x[l], x[j], dtemp);
    CC_SWAP (list[l], list[j], itemp);
    x_quicksort (list, x, l, j - 1);
    x_quicksort (list, x, i, u);
}

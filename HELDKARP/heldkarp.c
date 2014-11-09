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
/*                 HELD-KARP FOR SMALL TSP INSTANCES                        */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*               Sanjeeb Dash (Prim's Algorithm)                            */
/*  Date: February 17, 1998                                                 */
/*        October 23, 2003 (bico, add code to return the tour)              */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCheldkarp_small (int ncount, CCdatagroup *dat, double *upbound,    */
/*      double *optval, int *foundtour, int anytour, int *tour_elist,       */
/*      int nodelimit, int silent)                                          */
/*    -ncount is the number of nodes in the graph.                          */
/*    -dat specifies the information needed to compute the edge lengths.    */
/*    -upbound is an upperbound on the optimal tour length (it can be       */
/*     NULL)                                                                */
/*    -optval returns the length of an optimal tour (or the upperbound      */
/*     if no better tour is found, or the value of the first tour found     */
/*     that it better than upbound if anytour = 1 is specified)             */
/*    -foundtour will be set to 1 if a tour better than upbound is found    */
/*    -anytour should be set to 1 to cut off the search after any tour      */
/*     better than upbound is found (it may not be an optimal tour)         */
/*    -tour_elist will return an edgelist for the tour in end0 end1         */
/*     format; if not NULL then it should point to an array of length       */
/*     at least 2*ncount.  Can be NULL.                                     */
/*    -nodelimit specifies a limit on the number of search nodes (use -1    */
/*     to impose no limit)                                                  */
/*    -silent should be set to 1 to restrict the output and 2 to            */
/*     disable all normal output                                            */
/*                                                                          */
/*  int CCheldkarp_small_elist (int ncount, int ecount, int *elist,         */
/*      int *elen, int *upbound, int *optval, int *foundtour,               */
/*      int anytour, int *tour_elist, int nodelimit, int silent)            */
/*     USES edgelist rather than datagroup.                                 */
/*      -ecount is the number of edges in the graph.                        */
/*      -elist is the list of edges in end0 end1 format.                    */
/*      -elen is a list of the edge lengths.                                */
/*                                                                          */
/*    NOTES: The upperbound will be converted to an int.                    */
/*           Graph can have at most MAX_NODES with edge lengths no greater  */
/*           than  WEIGHT_MAX_EDGE                                          */
/*           The code was designed for problems in the range of 25 to 35    */
/*           nodes.                                                         */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "heldkarp.h"
#include "util.h"
#include "macrorus.h"

#define LINE_LEN (75)

#define MAX_NODES  (100)
#define WEIGHT_ADJUST (5)
#define WEIGHT_MULT (1 << WEIGHT_ADJUST)
#define WEIGHT_MAX_EDGE (1 << (20 - WEIGHT_ADJUST))  /* no overflow */ 
#define WEIGHT_MAX_NODE (1 << 21)

typedef struct treenode {
    int  deg;
    int  parent;
    int  mark;
    int  *adj;
    int  *eadj;
    int  parentedge;
    int  parentlen;
} treenode;


static void
    initial_y (int ncount, int ecount, int *elist, int *len, int *y),
    hk_work (int ncount, int *elist, int *elen, int *len, int **adjlist,
        int *zadjlist, int *y, int *deg, int *upperbound, int *tree,
        int *foundtour, int *besttour, int *efix, int *degfix, int depth,
        int *bbcount, int just_verify, int silent, int nodelimit),
    held_karp_bound (int ncount, int *elist, int *elen, int *len,
        int **adjlist, int *zadjlist, int *y, int *deg, int upperbound,
        int *tree, int *val, int *newtour, int *besttour, int maxiter,
        double beta, int silent),
    one_tree (int ncount, int *elist, int *len, int **adjlist, int *zadjlist,
        int *y, int *tree, int *notree),
    span_tree(int nnodes, int **adjlist, int elen[], int y[], int sptree[],
        int *notree),
    edge_select (int ncount, int *elist, int *len, int *y, int *tree,
        int *efix, int *ebranch),
    set_adjlist (int n0, int n1, int **adjlist, int *zadjlist, int val);


int CCheldkarp_small (int ncount, CCdatagroup *dat, double *upbound,
        double *optval, int *foundtour, int anytour, int *tour_elist,
        int nodelimit, int silent)
{
    int rval = 0;
    int i, j, k, ecount;
    int *elist = (int *) NULL;
    int *elen  = (int *) NULL;

    ecount = ncount * (ncount-1) / 2;
    elist = CC_SAFE_MALLOC (ecount*2, int);
    elen  = CC_SAFE_MALLOC (ecount, int);
    if (elist == (int *) NULL || elen == (int *) NULL) {
        fprintf (stderr, "out of memory in CCheldkarp_small\n");
        rval = HELDKARP_ERROR; goto CLEANUP;
    }
    for (i = 0, k = 0; i < ncount; i++) {
        for (j = 0; j < i; j++) {
            elist[2*k] = i;
            elist[2*k+1] = j;
            elen[k] = CCutil_dat_edgelen (i, j, dat);
            k++;
        }
    }

    rval = CCheldkarp_small_elist (ncount, ecount, elist, elen, upbound,
                                   optval, foundtour, anytour, tour_elist,
                                   nodelimit, silent);

CLEANUP:

    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);

    return rval;
}

/* In adjacency list
 *   if (i,j) = k'th edge (starting from 0), then adj(i,j) = k+1
 *   adj(i,j) = 0 => undefined edge.
 *   adj(i,j) < 0 => -adj(i,j) = k+1 (means edge is fixed to 1)
 */

int CCheldkarp_small_elist (int ncount, int ecount, int *elist, int *elen,
        double *upbound, double *optval, int *foundtour, int anytour,
        int *tour_elist, int nodelimit, int silent)
{
    int rval    = 0;
    int bbcount = 0;
    int init_ub = ncount * WEIGHT_MAX_EDGE + 1;
    int n1, n2, i, upperbound, val;
    int *p;
    int **adjlist = (int **) NULL;
    int *padjlist = (int *) NULL;
    int *zadjlist = (int *) NULL;
    int *degfix   = (int *) NULL;
    int *tree     = (int *) NULL;
    int *efix     = (int *) NULL;
    int *len      = (int *) NULL;
    int *deg      = (int *) NULL;
    int *y        = (int *) NULL;
    int *besttour = (int *) NULL;

    *foundtour = 0;

    if (upbound) upperbound = (int) (*upbound);
    else         upperbound = init_ub;

    val = upperbound;

    if (ncount > MAX_NODES) {
        fprintf (stderr, "too many nodes\n");
        rval = HELDKARP_ERROR; goto CLEANUP;
    }
    for (i = 0; i < ecount; i++) {
        if (elen[i] > WEIGHT_MAX_EDGE || -elen[i] > WEIGHT_MAX_EDGE) {
            fprintf (stderr, "edge too long\n"); 
            rval = HELDKARP_ERROR; goto CLEANUP;
        }
    }

    /* build adjlist for graph with node 0 deleted */

    adjlist  = CC_SAFE_MALLOC (ncount-1, int *);
    padjlist = CC_SAFE_MALLOC ((ncount-1)*(ncount-1), int);
    zadjlist = CC_SAFE_MALLOC (ncount, int);
    len      = CC_SAFE_MALLOC (ecount, int);
    if (adjlist  == (int **) NULL || padjlist == (int *) NULL ||
        zadjlist == (int *) NULL  ||      len == (int *) NULL) {
       fprintf (stderr, "out of memory in tiny_heldkarp\n");
       rval = HELDKARP_ERROR; goto CLEANUP;
    }
    for (i=0, p = padjlist; i<ncount-1; i++, p += (ncount-1)) {
       adjlist[i] = p;
    }
    for (i = 0; i < (ncount-1)*(ncount-1); i++) padjlist[i] = 0;
    for (i = 0; i < ncount; i++) zadjlist[i] = 0;

    /* fill in edge # in adj list; 0 stands for no edge; i+1 <-> edge i */

    for (i=0; i<ecount; i++){
       len[i] = (elen[i] << WEIGHT_ADJUST);
       n1 = elist[2*i]; n2 = elist[2*i+1];
       if (n1 == 0) {
           zadjlist[n2] = i+1;
       } else if (n2 == 0) {
           zadjlist[n1] = i+1;
       } else {
           adjlist[n1-1][n2-1] = adjlist[n2-1][n1-1] = i+1;
       }
    }

    tree = CC_SAFE_MALLOC (ncount, int);
    y    = CC_SAFE_MALLOC (ncount, int);
    deg  = CC_SAFE_MALLOC (ncount, int);
    besttour = CC_SAFE_MALLOC (ncount, int);
    if (tree == (int *) NULL || y == (int *) NULL || deg == (int *) NULL ||
        besttour == (int *) NULL) {
        fprintf (stderr, "out of memory in tiny_heldkarp\n"); 
        rval = HELDKARP_ERROR; goto CLEANUP;
    }
    initial_y (ncount, ecount, elist, len, y);

    efix   = CC_SAFE_MALLOC (ecount, int);
    degfix = CC_SAFE_MALLOC (ncount, int);
    if (efix == (int *) NULL || degfix == (int *) NULL) {
        fprintf (stderr, "out of memory in tiny_heldkarp\n"); 
        rval = HELDKARP_ERROR; goto CLEANUP;
    }
    for (i = 0; i < ecount; i++) efix[i] = 0;
    for (i = 0; i < ncount; i++) degfix[i] = 0;

    hk_work (ncount, elist, elen, len, adjlist, zadjlist, y, deg, &val, tree,
             foundtour, besttour, efix, degfix, 0, &bbcount, anytour, silent,
             nodelimit);
    if (silent<2) { printf ("BBnodes: %d\n", bbcount); fflush (stdout); }

    if (nodelimit != -1 && bbcount > nodelimit) {
        rval = HELDKARP_SEARCHLIMITEXCEEDED;
    } else {
        *optval = (double) val;
    }

    if (*foundtour && tour_elist) {
        for (i = 0; i < ncount; i++) {
            tour_elist[2*i] = elist[2*besttour[i]];
            tour_elist[2*i+1] = elist[2*besttour[i]+1];
        }
    }

CLEANUP:

    CC_IFFREE (adjlist, int *);
    CC_IFFREE (padjlist, int);
    CC_IFFREE (zadjlist, int);
    CC_IFFREE (degfix, int);
    CC_IFFREE (tree, int);
    CC_IFFREE (efix, int);
    CC_IFFREE (len, int);
    CC_IFFREE (deg, int);
    CC_IFFREE (y, int);
    CC_IFFREE (besttour, int);
    return rval;
}

static void initial_y (int ncount, int ecount, int *elist, int *len, int *y)
{
    int i;

    for (i = 0; i < ncount; i++) y[i] = INT_MAX;
    for (i=0; i < ecount; i++) {
        if (len[i] < y[elist[2*i]]) y[elist[2*i]] = len[i];
        if (len[i] < y[elist[2*i+1]]) y[elist[2*i+1]] = len[i];
    }
    for (i = 0; i < ncount; i++) {
        y[i] /= 2;
    }
}

static void hk_work (int ncount, int *elist, int *elen, int *len,
        int **adjlist, int *zadjlist, int *y, int *deg, int *upperbound,
        int *tree, int *foundtour, int *besttour, int *efix, int *degfix,
        int depth, int *bbcount, int just_verify, int silent, int nodelimit)
{
    int ebranch, n0, n1, maxiter, val, newtour;
    double beta;

    (*bbcount)++;
    if (nodelimit != -1 && *bbcount > nodelimit) return;
    maxiter = (depth > 0 ? 10 : 1000); 
    beta    = (depth > 0 ? 0.9 : 0.99);
    held_karp_bound (ncount, elist, elen, len, adjlist, zadjlist, y, deg,
                     *upperbound, tree, &val, &newtour, besttour, maxiter,
                     beta, silent);
    if (newtour == 1) { *foundtour = 1; *upperbound = val; return; }
    if (val >= *upperbound) return;

    edge_select (ncount, elist, len, y, tree, efix, &ebranch);
    if (ebranch == -1) return;
    n0 = elist[2*ebranch];
    n1 = elist[2*ebranch+1];
    set_adjlist (n0, n1, adjlist, zadjlist, 0);

    if (!silent && depth < LINE_LEN) { printf ("0"); fflush (stdout); }
    hk_work (ncount, elist, elen, len, adjlist, zadjlist, y, deg, upperbound,
             tree, foundtour, besttour, efix, degfix, depth+1, bbcount,
             just_verify, silent, nodelimit);
    if (!silent && depth < LINE_LEN) { printf ("\b \b"); fflush (stdout); }
    if (*foundtour == 1 && just_verify == 1) {
        set_adjlist (n0, n1, adjlist, zadjlist, ebranch+1); return;
    }

    if (degfix[n0] < 2 && degfix[n1] < 2) {
        efix[ebranch] = 1;
        degfix[n0]++;
        degfix[n1]++;
        set_adjlist (n0, n1, adjlist, zadjlist, -(ebranch+1));

        if (!silent && depth < LINE_LEN) { printf ("1"); fflush (stdout); }
        hk_work (ncount, elist, elen, len, adjlist, zadjlist, y, deg,
                 upperbound, tree, foundtour, besttour, efix, degfix, depth+1,
                 bbcount, just_verify, silent, nodelimit);
        if (!silent && depth < LINE_LEN) { printf ("\b \b"); fflush (stdout); }

        efix[ebranch] = 0;
        degfix[n0]--;
        degfix[n1]--;
    }
    set_adjlist (n0, n1, adjlist, zadjlist, ebranch+1);
}

static void held_karp_bound (int ncount, int *elist, int *elen, int *len,
        int **adjlist, int *zadjlist, int *y, int *deg, int upperbound,
        int *tree, int *val, int *newtour, int *besttour, int maxiter,
        double beta, int silent)
{
    int i, k, t, tlen, ysum, square, notree, newsum;
    int abound    = (upperbound << WEIGHT_ADJUST);
    int goal      = ((upperbound-1) << WEIGHT_ADJUST);
    int bestbound = -INT_MAX;
    int iter      = 0;
    double alpha  = 2.0;

    *newtour = 0;
    for (i = 0, ysum = 0; i < ncount; i++) {
       ysum += y[i];
    }

    do {
        one_tree (ncount, elist, len, adjlist, zadjlist, y, tree, &notree);
        if (notree == 1) {
            *val = INT_MAX; return;
        }
        for (i = 0, tlen = 2*ysum; i < ncount; i++) {
            k = tree[i];
            tlen += (len[k] - y[elist[2*k]] - y[elist[2*k+1]]);
        }
        if (tlen > bestbound) bestbound = tlen;
        if (tlen > goal) break;
  
        for (i = 0; i < ncount; i++) deg[i] = 2;
        for (i = 0; i < ncount; i++) {
            deg[elist[2*tree[i]]]--;
            deg[elist[2*tree[i]+1]]--;
        }
        for (i = 1, square = 0; i < ncount; i++) {
            square += ((deg[i])*(deg[i]));
        }
        if (square == 0) {
            *newtour = 1;
            for (i = 0; i < ncount; i++) {
                besttour[i] = tree[i];
            }
            for (i = 0, *val = 0; i < ncount; i++) *val += elen[tree[i]];
            if (silent < 2) {
                printf ("Tour found: %d\n", *val); fflush (stdout);
            }
            return;
        }
        if (++iter >= maxiter) break;

        t =  (int) (alpha * (double) ((abound-tlen)) / (double) square);
        if (t < 2) break;
        alpha *= beta;

        newsum = 0;
        for (i = 1; i < ncount; i++) {
            y[i] += (t*(deg[i]));   /* ysum does not change */
            if (y[i] > WEIGHT_MAX_NODE) {  /*  prevent overflow */
                y[i] = WEIGHT_MAX_NODE; newsum = 1; 
            }
        }
        if (newsum) {
            for (i = 0, ysum = 0; i < ncount; i++) ysum += y[i];
        }
    } while (1);

    *val = (bestbound >> WEIGHT_ADJUST);
    if (bestbound % WEIGHT_MULT) (*val)++;
}

static void one_tree (int ncount, int *elist, int *len, int **adjlist,
        int *zadjlist, int *y, int *tree, int *notree)
{
    int min1, min2, emin1, emin2, i, w, e;

    *notree = 0;
    span_tree (ncount-1, adjlist, len, y+1, tree, notree);
    if (*notree) return;

    min1 = INT_MAX;
    min2 = INT_MAX;
    emin1 = -1;
    emin2 = -1;

    for (i = 1; i < ncount; i++) {
        e = zadjlist[i];
        if (e > 0) {
            e--;
            w = len[e] - y[elist[2*e]] - y[elist[2*e+1]];
            if (w < min1) {
                min2 = min1;
                min1 = w;
                emin2 = emin1;
                emin1 = e;
            } else if (w < min2) {
                min2 = w;
                emin2 = e;
            }
        } else if (e < 0) {
            if (min2 == -INT_MAX) { *notree = 1; return; }
            min2 = min1;
            min1 = -INT_MAX;
            emin2 = emin1;
            emin1 = (-e)-1;
        }
    }
    if (min2 == INT_MAX) {
        *notree = 1;
    } else {
        tree[ncount-2] = emin1;
        tree[ncount-1] = emin2;
    }
}

static void span_tree (int nnodes, int **adjlist, int elen[], int y[],
        int sptree[], int *notree)
{
    int nbnd,nadd,nrem;
    int cur,minnode,i;
    int we,e,ycur;
    int nremain[MAX_NODES];
    int nedge[MAX_NODES];
    int nlen[MAX_NODES];
    int minlen;
    int *tptr;
 
    nbnd = nnodes-1;
    for (i=0; i<nbnd; i++){
        nremain[i] = i+1;
        nedge[i] = 0;
        nlen[i] = INT_MAX;
    }
    cur = 0;
    nbnd = nnodes-1;
    nrem = nbnd;
    for (nadd = 0; nadd<nbnd; nadd++){
        minlen = INT_MAX;
        minnode = -1;
        ycur = y[cur];
        tptr = adjlist[cur];

        for (i=nrem-1; i>=0; i--){
            e = tptr[nremain[i]];
            if (e > 0 ) {
                we = elen[e-1] - ycur - y[nremain[i]];
                if (we < nlen[i]){
                    nedge[i] = e;
                    nlen[i] = we;
                }
            } else if (e < 0) {      /*  => edge is fixed */
                if (nlen[i] == -INT_MAX){
                    *notree = 1;
                    return;
                } else {
                    nedge[i] = -e;
                    nlen[i] = -INT_MAX;
                }
            }
            if (nlen[i] < minlen){
                minlen = nlen[i];
                minnode = i;
            }
        }
        if (minnode == -1) {   /* Graph not connected */
            *notree = 1;  
            return;
        } else {
            e = nedge[minnode]-1;
            sptree[nadd] = e;
            cur = nremain[minnode];
            nrem--;
            nremain[minnode] = nremain[nrem];
            nedge[minnode] = nedge[nrem];
            nlen[minnode] = nlen[nrem];
        }
    }
}

static void edge_select (int ncount, int *elist, int *len, int *y, int *tree,
        int *efix, int *ebranch)
{
    int i, e, w;
    int min = INT_MAX;
    int emin = -1;

    for (i = 0; i < ncount; i++) {
        e = tree[i];
        if (efix[e] == 0) {
            w = len[e] - y[elist[2*e]] - y[elist[2*e+1]];
            if (w < min) { min = w; emin = e; }
        }
    }
    *ebranch = emin;
}

static void set_adjlist (int n0, int n1, int **adjlist, int *zadjlist,
        int val)
{
    if      (n0 == 0) zadjlist[n1] = val;
    else if (n1 == 0) zadjlist[n0] = val; 
    else              adjlist[n0-1][n1-1] = adjlist[n1-1][n0-1] = val;

}

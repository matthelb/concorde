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

#include "machdefs.h"

#define MAXNODE 16
#define MAXEDGE (MAXNODE * (MAXNODE-1)/2)
#define MAXNUM 20000
#define MAXCUTSIZE 1000

#define EDGENUM(i,j) ((i<j) ? (j*(j-1)/2+i) : (i*(i-1)/2+j))

#define CC_SWAP(a,b,t) ((t)=(a),(a)=(b),(b)=(t))

#define ABS(a) ((a<0)?(-a):a)

static int lhs[MAXEDGE];
static int rhs;
static int minrhs;
static int nnodes, nedges;
static int ntours = 0;
static int ntight = 0;
static int nindep = 0;
static double ind_tri[MAXEDGE][MAXEDGE];
static int elim_coord[MAXEDGE];
static char minval[MAXEDGE];
static char maxval[MAXEDGE];


typedef struct cut {
    int node[MAXCUTSIZE];
    int size;
    struct cut *next;
} cut;

static int nodetype[MAXNUM];
static int typesize[MAXNUM];
static int ntypes = 0;

static void
    check_indep (char *tour, char *touredges),
    checktour (char *tour, char *touredges),
    enumtours (int d, char *tour, char *touredges),
    split_type (cut *c, int t);


static void check_indep (char *tour, char *touredges)
{
    int i, j;
    double mul;

    for (i=0; i<nedges; i++) {
        ind_tri[nindep][i] = touredges[i];
        if (touredges[i] < minval[i]) minval[i] = touredges[i];
        if (touredges[i] > maxval[i]) maxval[i] = touredges[i];
    }
    for (i=0; i<nindep; i++) {
        if (ind_tri[nindep][elim_coord[i]]) {
            mul = ind_tri[nindep][elim_coord[i]] / ind_tri[i][elim_coord[i]];
            for (j=0; j<nedges; j++) {
                ind_tri[nindep][j] -= mul * ind_tri[i][j];
            }
        }
    }
    j = 0;
    mul = ABS(ind_tri[nindep][0]);
    for (i=1; i<nedges; i++) {
        if (ABS(ind_tri[nindep][i]) > mul) {
            mul = ABS(ind_tri[nindep][i]);
            j = i;
        }
    }
    if (mul > 1e-6) {
        elim_coord[nindep] = j;
        printf ("%d:",nindep+1);
        for (i=0; i<nnodes; i++) printf (" %d",tour[i]);
        printf (" (%g)\n",mul);
/*
        for (i=0; i<nedges; i++) {
            printf ("%g ",ind_tri[nindep][i]);
        }
        printf ("\n");
*/
        nindep++;
    }
}

static void checktour (char *tour, char *touredges)
{
    int v;
    int i;

    ntours++;

    for (i=0, v=0; i<nedges; i++) {
        v += touredges[i] * lhs[i];
    }

    if (v < rhs) {
        fprintf (stderr, "Whoa, tour ");
        for (i=0; i<nnodes; i++) fprintf (stderr, "%d ",tour[i]);
        fprintf (stderr, "has rhs %d\n",v);
        minrhs = v;
    }
    if (v == rhs) {
        ntight++;
        check_indep (tour, touredges);
    }
}

static void enumtours (int d, char *tour, char *touredges)
{
    int t;
    int i;

    if (d == nnodes) {
        if (tour[1] < tour[nnodes-1]) {
            touredges[EDGENUM(tour[0],tour[nnodes-1])] = 1;
            checktour (tour, touredges);
            touredges[EDGENUM(tour[0],tour[nnodes-1])] = 0;
        }
    } else {
        for (i=d; i<nnodes; i++) {
            CC_SWAP (tour[d], tour[i], t);
            touredges[EDGENUM(tour[d-1],tour[d])] = 1;
            enumtours (d+1, tour, touredges);
            touredges[EDGENUM(tour[d-1],tour[d])] = 0;
            CC_SWAP (tour[d], tour[i], t);
        }
    }
}

static void split_type (cut *c, int t)
{
    int cnt;
    int i;

    for (i=0, cnt=0; i<c->size; i++) {
        if (nodetype[c->node[i]] == t) cnt++;
    }
    if (cnt > 0 && cnt < typesize[t]) {
        typesize[ntypes] = 0;
        for (i=0; i<c->size; i++) {
            if (nodetype[c->node[i]] == t) {
                nodetype[c->node[i]] = ntypes;
                typesize[ntypes]++;
                typesize[t]--;
            }
        }
        ntypes++;
    }
}

int main (CC_UNUSED int ac, CC_UNUSED char **av)
{
    char touredges[MAXEDGE];
    char tour[MAXNODE];
    int i, j, k;
    int mark[MAXNODE];
    cut *cutlist = (cut *) NULL;
    cut *c;

    for (i=0; i<MAXNUM; i++) nodetype[i] = 0;
    typesize[0] = MAXNUM;
    ntypes = 1;

    printf ("cut: "); fflush (stdout);
    scanf ("%d", &i);
    while (i != -1) {
        c = (cut *) malloc (sizeof (cut));
        if (!c) {
            fprintf (stderr, "Out of memory\n");
            exit (1);
        }
        c->size = 0;
        while (i != -1) {
            if (i < 0 || i >= MAXNUM) {
                fprintf (stderr, "Can only handle nodes >= 0 and < %d\n",
                         MAXNUM);
                exit (1);
            }
            c->node[c->size++] = i;
            scanf ("%d", &i);
        }
        c->next = cutlist;
        cutlist = c;
        printf ("cut: "); fflush (stdout);
        scanf ("%d", &i);
    }

    printf ("rhs: "); fflush (stdout);
    scanf ("%d", &rhs);

    for (c = cutlist; c; c = c->next) {
        for (i=0; i<ntypes; i++) {
            split_type (c, i);
        }
    }

    printf ("%d node types\n", ntypes);
    if (ntypes > MAXNODE) {
        fprintf (stderr, "Too many types\n");
        exit (1);
    }

    printf ("nodemap:");
    for (i=0; i<MAXNUM; i++) {
        if (nodetype[i] != 0) printf (" %5d", i);
    }
    printf ("\n        ");
    for (i=0; i<MAXNUM; i++) {
        if (nodetype[i] != 0) printf (" %5d", nodetype[i]);
    }
    printf ("\n");

    nnodes = ntypes;
    nedges = nnodes * (nnodes-1) / 2;
    for (i=0; i<nedges; i++) {
        lhs[i] = 0;
    }

    for (c = cutlist; c; c = c->next) {
        for (i=0; i<nnodes; i++) mark[i] = 0;
        for (i=0; i<c->size; i++) mark[nodetype[c->node[i]]] = 1;
        for (i=0; i<nnodes; i++) {
            for (j=i+1; j<nnodes; j++) {
                if (mark[i] != mark[j]) {
                    lhs[EDGENUM(i,j)]++;
                }
            }
        }
    }

    printf ("ineq:");
    for (i=0; i<nnodes; i++) {
        for (j=i+1; j<nnodes; j++) {
            if (lhs[EDGENUM(i,j)]) {
                printf (" %d*(%d-%d)", lhs[EDGENUM(i,j)], i, j);
            }
        }
    }
    printf (" >= %d\n", rhs);
    fflush (stdout);

    for (i=0; i<nedges; i++) {
        touredges[i] = 0;
        minval[i] = 1;
        maxval[i] = 0;
    }

    for (i=0; i<nnodes; i++) tour[i] = i;

    minrhs = 1<<30;

    enumtours (1, tour, touredges);

    printf ("%d tours, %d tight, %d independent\n", ntours, ntight, nindep);
    if (minrhs < rhs) {
        printf ("VIOLATED - minimum %d\n", minrhs);
    } else if (nindep == nedges - nnodes) {
        printf ("FACET!!!\n");
    } else {
        printf ("VALID, but not FACET\n");
    }
    for (i=0; i<nnodes; i++) {
        for (j=i+1; j<nnodes; j++) {
            k = EDGENUM(i,j);
            if (minval[k] == maxval[k]) {
                printf ("Edge %d,%d always %d\n",i,j,minval[k]);
            }
        }
    }
    return 0;
}

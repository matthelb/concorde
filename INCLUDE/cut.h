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
/****************************************************************************/
/*                                                                          */
/*                      PROTOTYPES FOR FILES IN CUT                         */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/


#ifndef  __CUT_H
#define  __CUT_H

#include "util.h"

#define CC_LINSUB_NO_END 0
#define CC_LINSUB_LEFT_END 1
#define CC_LINSUB_RIGHT_END 2
#define CC_LINSUB_BOTH_END 3

#define CC_MINCUT_BIGDOUBLE   (1e30)
#define CC_MINCUT_ONE_EPSILON (0.000001)


int
    CCcut_mincut (int ncount, int ecount, int *elist, double *dlen,
        double *cutval, int **cut, int *cutcount),
    CCcut_violated_cuts (int ncount, int ecount, int *elist, double *dlen,
        double cutoff, int (*doit_fn) (double, int, int *, void *),
        void *pass_param),
    CCcut_shrink_cuts (int ncount, int ecount, int *elist, double *dlen,
        double cutoff, int (*doit_fn) (double, int, int *, void *),
        void *pass_param),
    CCcut_mincut_containing_set (int ncount, int ecount, int *elist,
        double *dlen, int scount, int *slist, double *cutval, int **cut,
        int *cutcount, int quickshrink, CCrandstate *rstate),
    CCcut_mincut_st (int ncount, int ecount, int *elist, double *ecap,
        int s, int t, double *value, int **cut, int *cutcount),
    CCcut_linsub (int ncount, int ecount, int *endmark, int *elist, double *x,
        double maxval, void *u_data, int (*cut_callback) (double cut_val,
        int cut_start, int cut_end, void *u_data)),
    CCcut_linsub_allcuts (int ncount, int ecount, int *perm, int *endmark,
        int *elist, double *x, double maxval, void *u_data,
        int (*cut_callback) (double cut_val, int cut_start, int cut_end,
        void *u_data)),
    CCcut_connect_components (int ncount, int ecount, int *elist, double *x,
        int *ncomp, int **compscount, int **comps);



/****************************************************************************/
/*                                                                          */
/*                           gomoryhu.c                                     */
/*                                                                          */
/****************************************************************************/

typedef struct CC_GHnode {
    struct CC_GHnode *parent;
    struct CC_GHnode *sibling;
    struct CC_GHnode *child;
    double            cutval;
    int               ndescendants;
    int               special;
    int              *nlist;
    int               listcount;
    int               num;
} CC_GHnode;

typedef struct CC_GHtree {
    struct CC_GHnode *root;
    struct CC_GHnode *supply;
    int              *listspace;
} CC_GHtree;



int
    CCcut_gomory_hu (CC_GHtree *T, int ncount, int ecount, int *elist,
        double *ecap, int markcount, int *marks, CCrandstate *rstate);

void
    CCcut_GHtreeinit (CC_GHtree *T),
    CCcut_GHtreefree (CC_GHtree *T),
    CCcut_GHtreeprint (CC_GHtree *T);




/****************************************************************************/
/*                                                                          */
/*                             shrink.c                                     */
/*                                                                          */
/****************************************************************************/

typedef struct CC_SRKnode {
    struct CC_SRKedge  *adj;
    struct CC_SRKnode  *next;
    struct CC_SRKnode  *prev;
    struct CC_SRKnode  *members;
    struct CC_SRKnode  *parent;
    struct CC_SRKnode  *qnext;
    double           prweight;
    double           weight;
    int              num;
    int              newnum;
    int              onecnt;
    int              onqueue;
    int              mark;
} CC_SRKnode;

typedef struct CC_SRKedge {
    struct CC_SRKnode  *end;
    struct CC_SRKedge  *other;
    struct CC_SRKedge  *next;
    struct CC_SRKedge  *prev;
    double           weight;
} CC_SRKedge;

typedef struct CC_SRKgraph {
    struct CC_SRKnode  *nodespace;
    struct CC_SRKedge  *edgespace;
    struct CC_SRKnode  *head;
    struct CC_SRKedge **hit;
    int              original_ncount;
    int              original_ecount;
    int              marker;
} CC_SRKgraph;

typedef struct CC_SRKexpinfo {
    int ncount;
    int *members;
    int *memindex;
} CC_SRKexpinfo;

typedef struct CC_SRKcallback {
    double cutoff;
    void *pass_param;
    int (*doit_fn) (double, int, int *, void *);
} CC_SRKcallback;


void
    CCcut_SRK_identify_paths (CC_SRKgraph *G, int *newcount, int onecnt_okay),
    CCcut_SRK_identify_paths_to_edges (CC_SRKgraph *G, int *newcount,
        int onecnt_okay),
    CCcut_SRK_identify_ones (CC_SRKgraph *G, int *count, double epsilon),
    CCcut_SRK_identify_one_triangles (CC_SRKgraph *G, int *count,
        CC_SRKnode *qstart, double epsilon, double cutoff, int unmarked),
    CCcut_SRK_identify_tight_triangles (CC_SRKgraph *G, int *count,
        double cutoff, int unmarked),
    CCcut_SRK_identify_tight_squares (CC_SRKgraph *G, int *count,
        double cutoff, int unmarked),
    CCcut_SRK_identify_triangle_square (CC_SRKgraph *G, int *count,
        double epsilon, int unmarked),
    CCcut_SRK_identify_one_square (CC_SRKgraph *G, int *count,
        double epsilon, double cutoff, int unmarked),
    CCcut_SRK_identify_nodes (CC_SRKgraph *G, CC_SRKnode *n, CC_SRKnode *m),
    CCcut_SRK_init_graph (CC_SRKgraph *G),
    CCcut_SRK_free_graph (CC_SRKgraph *G),
    CCcut_SRK_init_expinfo (CC_SRKexpinfo *expand),
    CCcut_SRK_free_expinfo (CC_SRKexpinfo *expand),
    CCcut_SRK_init_callback (CC_SRKcallback *cb),
    CCcut_SRK_increment_marker (CC_SRKgraph *G),
    CCcut_SRK_set_mark (CC_SRKgraph *G, int marker);

int
    CCcut_SRK_buildgraph (CC_SRKgraph *G, int ncount, int ecount, int *elist,
        double *dlen),
    CCcut_SRK_subtour_shrink (CC_SRKgraph *G, double *minval, double epsilon,
        CC_SRKcallback *cb, int **cut, int *cutcount),
    CCcut_SRK_crowder_padberg (CC_SRKgraph *G, double epsilon,
        CC_SRKcallback *cb),
    CCcut_SRK_identify_pr_edges (CC_SRKgraph *G, double *minval, int *count,
        CC_SRKnode *qstart, double epsilon, CC_SRKcallback *cb, int **cut,
        int *cutcount),
    CCcut_SRK_identify_set (CC_SRKgraph *G, int scount, int *slist),
    CCcut_SRK_defluff (CC_SRKgraph *G),
    CCcut_SRK_grab_edges (CC_SRKgraph *G, int *oncount, int *oecount,
        int **olist, double **olen, CC_SRKexpinfo *expand),
    CCcut_SRK_grab_nodes (CC_SRKgraph *G, CC_SRKexpinfo *expand),
    CCcut_SRK_trivial (int ncount, CC_SRKexpinfo *expand),
    CCcut_SRK_expand (CC_SRKexpinfo *expand, int *arr, int size, int **pnewarr,
        int *pnewsize),
    CCcut_SRK_original_ncount (CC_SRKexpinfo *expand);


#endif  /* __CUT_H */

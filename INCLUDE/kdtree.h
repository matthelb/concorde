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

#ifndef __KDTREE_H
#define __KDTREE_H

#include "util.h"

typedef struct CCkdnode {
    double cutval;
    struct CCkdnode *loson;
    struct CCkdnode *hison;
    struct CCkdnode *father;
    struct CCkdnode *next;
    struct CCkdbnds *bnds;
    int              lopt;
    int              hipt;
    char             bucket;
    char             empty;
    char             cutdim;
} CCkdnode;

typedef struct CCkdtree {
    CCkdnode        *root;
    CCkdnode       **bucketptr;
    int             *perm;
    CCptrworld       kdnode_world;
    CCptrworld       kdbnds_world;
} CCkdtree;

typedef struct CCkdbnds {
    double           x[2];
    double           y[2];
    struct CCkdbnds *next;
} CCkdbnds;


void
    CCkdtree_free (CCkdtree *kt),
    CCkdtree_delete (CCkdtree *kt, int k),
    CCkdtree_delete_all (CCkdtree *kt, int ncount),
    CCkdtree_undelete (CCkdtree *kt, int k),
    CCkdtree_undelete_all (CCkdtree *kt, int ncount);

int
    CCkdtree_build (CCkdtree *kt, int ncount, CCdatagroup *dat,
        double *wcoord, CCrandstate *rstate),
    CCkdtree_k_nearest (CCkdtree *kt, int ncount, int k, CCdatagroup *dat,
        double *wcoord, int wantlist, int *ocount, int **olist,
        int silent, CCrandstate *rstate),
    CCkdtree_quadrant_k_nearest (CCkdtree *kt, int ncount, int k,
        CCdatagroup *dat, double *wcoord, int wantlist, int *ocount,
        int **olist, int silent, CCrandstate *rstate),
    CCkdtree_node_k_nearest (CCkdtree *kt, int ncount, int n, int k,
        CCdatagroup *dat, double *wcoord, int *list, CCrandstate *rstate),
    CCkdtree_node_quadrant_k_nearest (CCkdtree *kt, int ncount, int n, int k,
        CCdatagroup *dat, double *wcoord, int *list, CCrandstate *rstate),
    CCkdtree_node_nearest (CCkdtree *kt, int n, CCdatagroup *dat,
        double *wcoord),
    CCkdtree_fixed_radius_nearest (CCkdtree *kt, CCdatagroup *dat,
        double *wcoord, int n, double rad, int (*doit_fn) (int, int, void *),
        void *pass_param),
    CCkdtree_nearest_neighbor_tour (CCkdtree *kt, int ncount, int start,
        CCdatagroup *dat, int *outcycle, double *val, CCrandstate *rstate),
    CCkdtree_nearest_neighbor_2match (CCkdtree *kt, int ncount, int start,
        CCdatagroup *dat, int *outmatch, double *val, CCrandstate *rstate),
    CCkdtree_prim_spanningtree (CCkdtree *kt, int ncount, CCdatagroup *dat,
        double *wcoord, int *outtree, double *val, CCrandstate *rstate),
    CCkdtree_greedy_tour (CCkdtree *kt, int ncount, CCdatagroup *dat,
        int *outcycle, double *val, int silent, CCrandstate *rstate),
    CCkdtree_far_add_tour (CCkdtree *kt, int ncount, int start,
        CCdatagroup *dat, int *outcycle, double *val, CCrandstate *rstate),
    CCkdtree_qboruvka_tour (CCkdtree *kt, int ncount, CCdatagroup *dat,
        int *outcycle, double *val, CCrandstate *rstate),
    CCkdtree_boruvka_tour (CCkdtree *kt, int ncount, CCdatagroup *dat,
        int *outcycle, double *val, CCrandstate *rstate),
    CCkdtree_twoopt_tour (CCkdtree *kt, int ncount, CCdatagroup *dat,
        int *incycle, int *outcycle, double *val, int run_two_and_a_half_opt,
        int silent, CCrandstate *rstate),
    CCkdtree_3opt_tour (CCkdtree *kt, int ncount, CCdatagroup *dat,
        int *incycle, int *outcycle, double *val, int silent,
        CCrandstate *rstate);


#endif  /* __KDTREE_H */


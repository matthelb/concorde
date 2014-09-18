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

#ifndef  __COMBS_H
#define  __COMBS_H

#define CCcombs_BLOCK_ZERO_EPSILON (1e-10)

typedef struct CC_GCedge {
    int         to;
    double      weight;
} CC_GCedge;

typedef struct CC_GCnode {
    int         deg;
    CC_GCedge  *adj;
    int         mark;
    int         qhandle;
    double      flow;
    int         setloc;
    int         setdeg;
    int         status;
} CC_GCnode;

typedef struct CC_GCgraph {
    int        ncount;
    int        ecount;
    CC_GCnode *nodelist;
    CC_GCedge *edgespace;
} CC_GCgraph;



int
    CCcombs_find_blocks (int ncount, int ecount, int *elist, double *x,
        int *nblocks, int **blockcnts, int ***blocks, int *ncutnodes,
        int **cutnodes),
    CCcombs_greedy_cut (CC_GCgraph *g, int *setsize, int *set, int mark_fixed,
        int forced_moves, int bad_moves, int fixed_moves, int *moves_done,
        double *cut_val),
    CCcombs_GC_build_graph (CC_GCgraph *G, int ncount, int ecount, int *elist,
        double *x);

void
    CCcombs_GC_init_graph (CC_GCgraph *G),
    CCcombs_GC_free_graph (CC_GCgraph *G);


#endif /* __COMBS_H */

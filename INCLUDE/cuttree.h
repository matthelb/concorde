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

#ifndef  __CUTTREE_H
#define  __CUTTREE_H

#include "tsp.h"
#include "pq.h"

void
    CCpq_cuttree_init (CCtsp_cuttree *t),
    CCpq_cuttree_freetree (CCtsp_cuttree *t),
    CCpq_check_clique (CCpq_tree *pqt, CCtsp_lpclique *c, int *status),
    CCpq_cuttree_display (CCtsp_cuttree *t),
    CCpq_cuttree_describe (CCtsp_cuttree *t);

int
    CCpq_cuttree_trivial (CCtsp_cuttree *t, int nodecount, int extern_node),
    CCpq_cuttree_update_clean (CCtsp_cuttree *t, int edgecount, int *elist,
        double *x),
    CCpq_cuttree_improve_quick (CCtsp_cuttree *t, CCtsp_lpcuts *pool,
        int edgecount, int *elist, double *x),
    CCpq_apply_clique (CCpq_tree *T, CCtsp_lpclique *c, int *status),
    CCpq_cuttree_gen_cliques (CCtsp_cuttree *t, void *u_data,
        int (*cut_callback) (int *arr, int cnt, int *stop, void *u_data)),
    CCpq_cuttree_build_necklaces (CCtsp_cuttree *t, int ecount, int *elist,
        double *x, int *p_neckcount, CCtsp_cutnode ***p_necklist,
        int *necknum);


#endif /* __CUTTREE_H */

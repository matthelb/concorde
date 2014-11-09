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

#ifndef __PQ_H
#define __PQ_H

#include "pqsets.h"
#include "util.h"

struct CCtsp_cuttree;

typedef struct CCpq_node {
    int number;
    struct CCpq_node *next;

    CCpq_elem queue_elem;

/* the size of the children_set will not necessarily be correct for Q nodes */
    CCpq_set children_set;
    CCpq_elem children_elem;

    CCpq_set full_children_set;
    CCpq_elem full_children_elem;

    CCpq_set partial_children_set;
    CCpq_elem partial_children_elem;

    CCpq_elem blocked_elem;

    CCpq_elem leaves_elem;

    struct CCpq_node *parent;

    int pertinent_child_count;
    int pertinent_leaf_count;

    int mark;
#define IS_UNINITIALIZED(x,T) ((x)->mark < (T)->markbase)
#define UNMARKED(T) ((T)->markbase+0)
#define QUEUED(T) ((T)->markbase+1)
#define BLOCKED(T) ((T)->markbase+2)
#define UNBLOCKED(T) ((T)->markbase+3)

    int type;
    int parenttype;
#define PQ_LEAF 0
#define PQ_PNODE 1
#define PQ_QNODE 2
#define PQ_EXTERN 3
#define PQ_ROOT 4

    int label;
#define IS_EMPTY(x,T) ((x)->label <= (T)->markbase)
#define EMPTY(T) ((T)->markbase + 0)
#define PARTIAL(T) ((T)->markbase + 1)
#define FULL(T) ((T)->markbase + 2)

} CCpq_node;

typedef struct CCpq_tree {
    int nodecount;
    int extern_node;
    CCpq_node *elems;
    CCpq_node *leaflist;
    int markbase;
    CCpq_node pseudo_root;
    int node_counter;
    int nontrivial;
    CCptrworld pqnode_world;
} CCpq_tree;

#define CCpq_STATUS_NOSOL 1
#define CCpq_STATUS_TRIVIAL 2
#define CCpq_STATUS_NONTRIVIAL 3
#define CCpq_STATUS_BUBBLEOK 4

#define CCpq_clear_leaflist(T) ((T)->leaflist = (CCpq_node *) NULL)
#define CCpq_add_leaflist(T,i) ((T)->elems[(i)].next = (T)->leaflist, \
                              (T)->leaflist = &((T)->elems[(i)]))
#define CCpq_set_leaflist(T,l) ((T)->leaflist = l)


void
    CCpq_tree_init (CCpq_tree *T),
    CCpq_tree_free (CCpq_tree *T),
    CCpq_describe_solution (CCpq_tree *T),
    CCpq_dump_solution (CCpq_tree *T);

int
    CCpq_check (CCpq_tree *T, int *status),
    CCpq_apply (CCpq_tree *T, int *status),
    CCpq_tree_trivial (CCpq_tree *T, int nodecount, int extern_node),
    CCpq_cuttree_to_pq (struct CCtsp_cuttree *ct, CCpq_tree *pqT);

CCpq_node
   *CCpq_find_root (CCpq_tree *T);


#endif  /* __PQ_H */

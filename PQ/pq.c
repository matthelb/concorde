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
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  void CCpq_tree_init (CCpq_tree *T)                                      */
/*    NONE                                                                  */
/*                                                                          */
/*  void CCpq_tree_free (CCpq_tree *T)                                      */
/*    NONE                                                                  */
/*                                                                          */
/*  void CCpq_describe_solution (CCpq_tree *T)                              */
/*    NONE                                                                  */
/*                                                                          */
/*  void CCpq_dump_solution (CCpq_tree *T)                                  */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCpq_check (CCpq_tree *T, int *status)                              */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCpq_apply (CCpq_tree *T, int *status)                              */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCpq_tree_trivial (CCpq_tree *T, int nodecount, int extern_node)    */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCpq_cuttree_to_pq (struct CCtsp_cuttree *ct, CCpq_tree *pqT)       */
/*    NONE                                                                  */
/*                                                                          */
/*  CCpq_node *CCpq_find_root (CCpq_tree *T)                                */
/*    NONE                                                                  */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "pq.h"
#include "tsp.h"
#include "cuttree.h"
#include "macrorus.h"

#define PQSTATUS_UNMATCHED 5
#define PQSTATUS_MATCHED 6


static void
    bubble (CCpq_tree *T, CCpq_node *l, int *status),
    template_l1 (CCpq_tree *T, CCpq_node *x, int *status),
    template_p1 (CCpq_tree *T, CCpq_node *x, int *status),
    template_q1 (CCpq_tree *T, CCpq_node *x, int *status),
    template_q2 (CCpq_tree *T, CCpq_node *x, int real, int *status),
    template_q3_root (CCpq_tree *T, CCpq_node *x, int real, int *status),
    label_full (CCpq_tree *T, CCpq_node *x),
    label_partial (CCpq_tree *T, CCpq_node *x),
    neighbor_replace (CCpq_node *x, CCpq_set *s, CCpq_node *q, CCpq_node *r),
    merge_qnode (CCpq_tree *T, CCpq_node *x),
    replace_node (CCpq_node *old, CCpq_node *new),
    node_init (CCpq_tree *T, CCpq_node *x),
    subtree_free (CCpq_node *x, CCptrworld *pqnode_world),
    describe_subtree (CCpq_node *x),
    dump_subtree (CCpq_node *x);

static int
    apply_checked (CCpq_tree *T, int *status),
    reduce (CCpq_tree *T, CCpq_node *l, int real, int *status),
    template_p3_notroot (CCpq_tree *T, CCpq_node *x, int real, int *status),
    template_p5_notroot (CCpq_tree *T, CCpq_node *x, int real, int *status),
    template_p2_root (CCpq_tree *T, CCpq_node *x, int real, int *status),
    template_p4_root (CCpq_tree *T, CCpq_node *x, int real, int *status),
    template_p6_root (CCpq_tree *T, CCpq_node *x, int real, int *status),
    collect_full_children (CCpq_tree *T, CCpq_node *x, int t),
    collect_empty_children (CCpq_tree *T, CCpq_node **p_x, int t),
    left_end_is_full (CCpq_tree *T, CCpq_node *x);

static CCpq_node
   *check_invert (CCpq_tree *T),
   *invert_leaflist (CCpq_tree *T),
   *cuttree_to_pqtree_work (CCtsp_cutnode *x, CCtsp_cutnode *nodelist,
        CCpq_tree *T);

#ifndef NDEBUG

static int
    check_tree (CCpq_tree *T),
    check_subtree (CCpq_tree *T, CCpq_node *x),
    check_node (CCpq_tree *T, CCpq_node *x),
    check_node_pert (CCpq_tree *T, CCpq_node *x),
    check_node_pert_work (CCpq_tree *T, CCpq_node *x);

#endif /* NDEBUG */


CC_PTRWORLD_ROUTINES (CCpq_node, PQ_node_alloc, PQ_node_bulkalloc,
                      PQ_node_free)
#ifndef NDEBUG
CC_PTRWORLD_LEAKS_ROUTINE (CCpq_node, PQ_node_leaks, mark, int)
#endif

/* This is an implementation of the PQ-tree algorithms described in "Testing
 * for the Consecutive Ones Property, Interval Graphs, and Graph Planarity
 * Using PQ-Tree Algorithms" by Kellogg Booth and George Leuker, Journal of
 * Computer and Systems Sciences 13 (1976) pp 335-379.  For details, proofs
 * of correctness, etc. see the paper */

void CCpq_tree_init (CCpq_tree *T)
{
    CCptrworld_init (&T->pqnode_world);
    T->elems = (CCpq_node *) NULL;
    T->leaflist = (CCpq_node *) NULL;
    T->nodecount = -1;
    T->extern_node = -1;
    T->markbase = 0;
    T->node_counter = 0;
    T->nontrivial = 0;
    node_init (T, &(T->pseudo_root));
}

void CCpq_tree_free (CCpq_tree *T)
{
    if (T->elems != (CCpq_node *) NULL) {
        subtree_free (CCpq_find_root (T), &T->pqnode_world);

        CCptrworld_delete (&T->pqnode_world);

        CC_FREE (T->elems, CCpq_node);
    }
    CCpq_tree_init (T);
}

int CCpq_check (CCpq_tree *T, int *status)
{
    int rval;
    CCpq_node *l;

    l = check_invert (T);
    bubble (T, l, status);
    if (*status == CCpq_STATUS_NOSOL) {
        return 0;
    }
    T->nontrivial = 0;
    rval = reduce (T, l, 0, status);
    if (rval) {
        fprintf (stderr, "reduce failed\n");
        return 1;
    }
    return 0;
} /* END PQ_CHECK */

static int apply_checked (CCpq_tree *T, int *status)
{
    int rval;
    CCpq_node *l;

    l = check_invert (T);
    bubble (T, l, status);
    if (*status == CCpq_STATUS_NOSOL) {
        return 0;
    }
    T->nontrivial = 0;
    rval = reduce (T, l, 1, status);
    if (rval) {
        fprintf (stderr, "reduce failed\n");
        return 1;
    }
    return 0;
} /* END apply_checked */

int CCpq_apply (CCpq_tree *T, int *status)
{
    int rval;

    rval = CCpq_check (T, status);
    if (rval) {
        fprintf (stderr, "CCpq_check failed\n");
        return 1;
    }
    if (*status == CCpq_STATUS_NONTRIVIAL) {
        rval = apply_checked (T, status);
        if (rval) {
            fprintf (stderr, "apply_checked failed\n");
            return 1;
        }
        if (*status != CCpq_STATUS_NONTRIVIAL) {
            fprintf (stderr, "ERROR: apply_checked status != CCpq_check status\n");
            return -1;
        }
    }
    return 0;
}

static CCpq_node *check_invert (CCpq_tree *T)
{
    CCpq_node *l;

    for (l=T->leaflist; l; l = l->next) {
        if (l->type == PQ_EXTERN) {
            return invert_leaflist (T);
        }
    }
    return T->leaflist;
}

static CCpq_node *invert_leaflist (CCpq_tree *T)
{
    CCpq_node *l;
    int i;
    int nodecount = T->nodecount;
    CCpq_node *elems = T->elems;

    for (i=0; i<nodecount; i++) {
        elems[i].label = 0;
    }
    for (l=T->leaflist; l; l = l->next) {
        l->label = 1;
    }
    l = (CCpq_node *) NULL;
    for (i=0; i<nodecount; i++) {
        if (elems[i].label == 0) {
            elems[i].next = l;
            l = &elems[i];
        } else {
            elems[i].label = 0;
        }
    }

    return l;
}

int CCpq_tree_trivial (CCpq_tree *T, int nodecount, int extern_node)
{
    int i;
    CCpq_node *root = (CCpq_node *) NULL;
    CCpq_node *elems = (CCpq_node *) NULL;
    int rval;

    CCpq_tree_free (T);

    if (nodecount < 3) {
        fprintf (stderr, "Can't build PQ tree with %d nodes\n", nodecount);
        rval = 1; goto CLEANUP;
    }

    elems = CC_SAFE_MALLOC (nodecount, CCpq_node);

    if (elems == (CCpq_node *) NULL) {
        fprintf (stderr, "Out of memory in CCpq_tree_trivial\n");
        rval = 1; goto CLEANUP;
    }

    root = PQ_node_alloc (&T->pqnode_world);
    if (root == (CCpq_node *) NULL) {
        fprintf (stderr, "PQ_node_alloc failed\n");
        rval = 1; goto CLEANUP;
    }

    T->nodecount = nodecount;
    T->extern_node = extern_node;
    T->elems = elems;
    T->leaflist = (CCpq_node *) NULL;
    T->markbase = 0;
    T->node_counter = -1;
    T->nontrivial = 0;

    node_init (T, root);
    root->type = PQ_PNODE;
    root->parent = (CCpq_node *) NULL;
    root->parenttype = PQ_PNODE;

    for (i=0; i<nodecount; i++) {
        elems[i].number = i;
        if (i == extern_node) {
            elems[i].type = PQ_EXTERN;
        } else {
            elems[i].parent = root;
            elems[i].parenttype = PQ_PNODE;
            CCpq_set_ADD (&elems[i], root->children_set, children_elem);
            elems[i].label = 0;
            elems[i].mark = 0;
            elems[i].type = PQ_LEAF;
            CCpq_set_INIT (elems[i].children_set);
            CCpq_set_INIT (elems[i].full_children_set);
            CCpq_set_INIT (elems[i].partial_children_set);
            elems[i].pertinent_child_count = 0;
            elems[i].pertinent_leaf_count = 0;
        }
    }

    rval = 0;
  CLEANUP:
    if (rval) {
        if (root != (CCpq_node *) NULL) PQ_node_free (&T->pqnode_world, root);
        CC_IFFREE (elems, CCpq_node);
    }
    return rval;
}

/* bubble forms the first half of the PQ-tree implementation.
 *
 * It takes as input a list of CCpq_node(s) linked by (next).
 *
 * It fills in the parent field of all the nodes in the PQ-tree, notices
 * certain unsatisfiable requests, fills in the pertinent_child_counts, ...
 */

static void bubble (CCpq_tree *T, CCpq_node *l, int *status)
{
    CCpq_set q;
    int block_count;
    CCpq_node *x;
    CCpq_node *zprev, *z, *znext;
    CCpq_node *leftsib, *rightsib;
    int i;
    int nbs;
    int off_the_top;
    CCpq_set blocked_set;

    assert (check_tree (T));

    T->markbase += 4;
    CCpq_set_INIT (q);
    block_count = 0;
    CCpq_set_INIT (blocked_set);
    off_the_top = 0;

    for (x = l; x; x = x->next) {
        CCpq_set_ADD_LEFT (x, q, queue_elem);
    }

    while (CCpq_set_SIZE (q) + block_count + off_the_top > 1) {
        if (CCpq_set_ISEMPTY (q)) {
            *status = CCpq_STATUS_NOSOL;
            return;
        }
        x = CCpq_set_RIGHT_ELEM (q);
        CCpq_set_DELETE (x, q, queue_elem);
        x->mark = BLOCKED(T);

        x->pertinent_leaf_count = 0;
        CCpq_set_INIT (x->full_children_set);
        CCpq_set_INIT (x->partial_children_set);

        nbs = 0;

        leftsib = x->children_elem.ptr1;
        rightsib = x->children_elem.ptr2;

        if (leftsib && leftsib->mark == BLOCKED(T)) {
            nbs++;
        }
        if (rightsib && rightsib->mark == BLOCKED(T)) {
            nbs++;
        }
        if (x->parenttype == PQ_PNODE || leftsib == NULL
            || rightsib == NULL) {
            x->mark = UNBLOCKED(T);
        } else if (leftsib->mark == UNBLOCKED(T)) {
            z = x;
            zprev = x->children_elem.ptr2;
            CCpq_set_FOREACH_FROM (z, children_elem, zprev, znext) {
                if (z == leftsib) {
                    break;
                }
                z->parent = leftsib->parent;
                z->mark = UNBLOCKED(T);
            }
            assert (z == leftsib);
        } else if (rightsib->mark == UNBLOCKED(T)) {
            z = x;
            zprev = x->children_elem.ptr1;
            CCpq_set_FOREACH_FROM (z, children_elem, zprev, znext) {
                if (z == rightsib) {
                    break;
                }
                z->parent = rightsib->parent;
                z->mark = UNBLOCKED(T);
            }
            assert (z == rightsib);
        }
        if (x->mark == UNBLOCKED(T)) {
            if (x->parent && IS_UNINITIALIZED (x->parent,T)) {
                x->parent->mark = UNMARKED(T);
                x->parent->pertinent_child_count = 0;
            }
            CCpq_set_FOREACH_ADJ (x, children_elem, z, i) {
                zprev = x;
                CCpq_set_FOREACH_FROM (z, children_elem, zprev, znext) {
                    if (z->mark == BLOCKED(T)) {
                        z->parent = x->parent;
                        z->mark = UNBLOCKED(T);
                        if (x->parent) {
                            x->parent->pertinent_child_count++;
                        }
                        CCpq_set_DELETE (z, blocked_set, blocked_elem);
                    } else {
                        break;
                    }
                }
            }

            if (x->parent) {
                x->parent->pertinent_child_count++;
                if (x->parent->mark == UNMARKED(T)) {
                    CCpq_set_ADD_LEFT (x->parent, q, queue_elem);
                    x->parent->mark = QUEUED(T);
                }
            } else {
                off_the_top = 1;
            }
            block_count -= nbs;
        } else {
            block_count = block_count + 1 - nbs;
            CCpq_set_ADD (x, blocked_set, blocked_elem);
        }
    }
    if (!CCpq_set_ISEMPTY (q)) {
        x = CCpq_set_LEFT_ELEM (q);
        x->pertinent_leaf_count = 0;
        if (x->pertinent_child_count == 1) {
            x->pertinent_child_count = 0;
        }
        CCpq_set_INIT (x->full_children_set);
        CCpq_set_INIT (x->partial_children_set);
    }
    if (block_count) {
        T->pseudo_root.type = PQ_QNODE;
        T->pseudo_root.parent = (CCpq_node *) NULL;
        T->pseudo_root.pertinent_child_count = 0;
        CCpq_set_FOREACH (blocked_set, x, blocked_elem, zprev, znext) {
            x->parent = &T->pseudo_root;
            T->pseudo_root.pertinent_child_count++;
        }
        T->pseudo_root.pertinent_leaf_count = 0;
        CCpq_set_INIT (T->pseudo_root.full_children_set);
        CCpq_set_INIT (T->pseudo_root.partial_children_set);
        T->pseudo_root.label = 0;
    }
    assert (check_tree (T));
    *status = CCpq_STATUS_BUBBLEOK;
    return;
} /* END bubble */

static int reduce (CCpq_tree *T, CCpq_node *l, int real, int *status)
{
    CCpq_set q;
    CCpq_node *x;
    int ssize;
    CCpq_node *y;
    int rval;

    assert (check_tree (T));

    CCpq_set_INIT (q);

    for (x = l, ssize = 0; x; x = x->next, ssize++) {
        CCpq_set_ADD_LEFT (x, q, queue_elem);
        x->pertinent_leaf_count = 1;
    }
    while (!CCpq_set_ISEMPTY (q)) {
        x = CCpq_set_RIGHT_ELEM (q);
        CCpq_set_DELETE (x, q, queue_elem);

        assert (check_node_pert (T, x));

        if (x->pertinent_leaf_count < ssize) {
            /* x is not the root of the pertinent subtree */
            y = x->parent;
            y->pertinent_leaf_count += x->pertinent_leaf_count;
            y->pertinent_child_count--;
            if (y->pertinent_child_count == 0) {
                CCpq_set_ADD_LEFT (y, q, queue_elem);
            }

            /* now, the templates */
            *status = PQSTATUS_UNMATCHED;

            template_l1 (T, x, status);
            if (*status == PQSTATUS_MATCHED) continue;

            template_p1 (T, x, status);
            if (*status == PQSTATUS_MATCHED) continue;

            rval = template_p3_notroot (T, x, real, status);
            if (rval) {
                fprintf (stderr, "template_p3_notroot failed\n");
                return 1;
            }
            if (*status == PQSTATUS_MATCHED) continue;

            rval = template_p5_notroot (T, x, real, status);
            if (rval) {
                fprintf (stderr, "template_p5_notroot failed\n");
                return 1;
            }
            if (*status == PQSTATUS_MATCHED) continue;

            template_q1 (T, x, status);
            if (*status == PQSTATUS_MATCHED) continue;

            template_q2 (T, x, real, status);
            if (*status == PQSTATUS_MATCHED) continue;
            *status = CCpq_STATUS_NOSOL;
            return 0;
        } else {
            /* x is the root of the pertinent subtree */

            *status = PQSTATUS_UNMATCHED;
            assert (CCpq_set_ISEMPTY (q));

            template_l1 (T, x, status);
            if (*status == PQSTATUS_MATCHED) {
                *status = (T->nontrivial) ? CCpq_STATUS_NONTRIVIAL
                                          : CCpq_STATUS_TRIVIAL;
                return 0;
            }

            template_p1 (T, x, status);
            if (*status == PQSTATUS_MATCHED) {
                *status = (T->nontrivial) ? CCpq_STATUS_NONTRIVIAL
                                          : CCpq_STATUS_TRIVIAL;
                return 0;
            }

            rval = template_p2_root (T, x, real, status);
            if (rval) {
                fprintf (stderr, "template_p2_root failed\n");
                return 1;
            }
            if (*status == PQSTATUS_MATCHED) {
                *status = (T->nontrivial) ? CCpq_STATUS_NONTRIVIAL
                                          : CCpq_STATUS_TRIVIAL;
                return 0;
            }

            rval = template_p4_root (T, x, real, status);
            if (rval) {
                fprintf (stderr, "template_p4_root failed\n");
                return 1;
            }
            if (*status == PQSTATUS_MATCHED) {
                *status = (T->nontrivial) ? CCpq_STATUS_NONTRIVIAL
                                          : CCpq_STATUS_TRIVIAL;
                return 0;
            }

            rval = template_p6_root (T, x, real, status);
            if (rval) {
                fprintf (stderr, "template_p6_root failed\n");
                return 1;
            }
            if (*status == PQSTATUS_MATCHED) {
                *status = (T->nontrivial) ? CCpq_STATUS_NONTRIVIAL
                                          : CCpq_STATUS_TRIVIAL;
                return 0;
            }

            template_q1 (T, x, status);
            if (*status == PQSTATUS_MATCHED) {
                *status = (T->nontrivial) ? CCpq_STATUS_NONTRIVIAL
                                          : CCpq_STATUS_TRIVIAL;
                return 0;
            }

            template_q2 (T, x, real, status);
            if (*status == PQSTATUS_MATCHED) {
                *status = (T->nontrivial) ? CCpq_STATUS_NONTRIVIAL
                                          : CCpq_STATUS_TRIVIAL;
                return 0;
            }

            template_q3_root (T, x, real, status);
            if (*status == PQSTATUS_MATCHED) {
                *status = (T->nontrivial) ? CCpq_STATUS_NONTRIVIAL
                                          : CCpq_STATUS_TRIVIAL;
                return 0;
            }

            *status = CCpq_STATUS_NOSOL;
            return 0;
        }
    }
    *status = CCpq_STATUS_NOSOL;
    return 0;
} /* END reduce */

static void template_l1 (CCpq_tree *T, CCpq_node *x, int *status)
{
    if (x->type != PQ_LEAF) {
        return;
    }
    label_full (T, x);
    *status = PQSTATUS_MATCHED;
    return;
} /* END template_l1 */

static void template_p1 (CCpq_tree *T, CCpq_node *x, int *status)
{
    if (x->type != PQ_PNODE) {
        return;
    }
    if (CCpq_set_SIZE (x->full_children_set) !=
        CCpq_set_SIZE (x->children_set)) {
        return;
    }
    label_full (T, x);
    *status = PQSTATUS_MATCHED;
    return;
} /* END template_p1 */

static int template_p3_notroot (CCpq_tree *T, CCpq_node *x, int real,
        int *status)
{
    int rval;

    if (x->type != PQ_PNODE) {
        return 0;
    }
    if (!CCpq_set_ISEMPTY (x->partial_children_set)) {
        return 0;
    }
    assert (!CCpq_set_ISEMPTY (x->full_children_set));
    assert (CCpq_set_SIZE (x->full_children_set) !=
            CCpq_set_SIZE (x->children_set));

    T->nontrivial = 1;

    if (real) {
        rval = collect_full_children (T, x, PQ_QNODE);
        if (rval) {
            fprintf (stderr, "collect_full_children failed\n");
            return 1;
        }
        rval = collect_empty_children (T, &x, PQ_QNODE);
        if (rval) {
            fprintf (stderr, "collect_empty_children failed\n");
            return 1;
        }
    }

    label_partial (T, x);
    *status = PQSTATUS_MATCHED;

    return 0;
} /* END template_p3_NOTROOT */

static int template_p5_notroot (CCpq_tree *T, CCpq_node *x, int real,
        int *status)
{
    CCpq_node *partial_child;
    CCpq_node *full_child;
    CCpq_node *z;
    int rval;

    if (x->type != PQ_PNODE) {
        return 0;
    }
    if (CCpq_set_SIZE (x->partial_children_set) != 1) {
        return 0;
    }

    T->nontrivial = 1;

    if (real) {
        rval = collect_full_children (T, x, PQ_QNODE);
        if (rval) {
            fprintf (stderr, "collect_full_children failed\n");
            return 1;
        }

        partial_child = CCpq_set_LEFT_ELEM (x->partial_children_set);
        CCpq_set_DELETE (partial_child, x->children_set, children_elem);
        CCpq_set_DELETE (partial_child, x->partial_children_set,
                       partial_children_elem);

        if (!CCpq_set_ISEMPTY (x->full_children_set)) {
            full_child = CCpq_set_LEFT_ELEM (x->full_children_set);
            CCpq_set_DELETE (full_child, x->children_set, children_elem);
            CCpq_set_DELETE (full_child, x->full_children_set,
                           full_children_elem);
            full_child->parent = partial_child;
            full_child->parenttype = PQ_QNODE;
            if (left_end_is_full (T, partial_child)) {
                CCpq_set_ADD_LEFT (full_child, partial_child->children_set,
                                 children_elem);
                CCpq_set_ADD_LEFT (full_child,
                        partial_child->full_children_set, full_children_elem);
            } else {
                CCpq_set_ADD_RIGHT (full_child, partial_child->children_set,
                                  children_elem);
                CCpq_set_ADD_RIGHT (full_child,
                        partial_child->full_children_set, full_children_elem);
            }
        }
        replace_node (x, partial_child);

        if (CCpq_set_ISEMPTY (x->children_set)) {
            PQ_node_free (&T->pqnode_world, x);
        } else {
            if (CCpq_set_SIZE (x->children_set) == 1) {
                z = CCpq_set_LEFT_ELEM (x->children_set);
                PQ_node_free (&T->pqnode_world, x);
                x = z;
            }
            x->parent = partial_child;
            x->parenttype = PQ_QNODE;
            if (left_end_is_full (T, partial_child)) {
                CCpq_set_ADD_RIGHT (x, partial_child->children_set,
                                  children_elem);
            } else {
                CCpq_set_ADD_LEFT (x, partial_child->children_set,
                                 children_elem);
            }
        }
        label_partial (T, partial_child);
    } else {
        label_partial (T, x);
    }
    *status = PQSTATUS_MATCHED;
    return 0;
} /* END template_p5_NOTROOT */

static void template_q1 (CCpq_tree *T, CCpq_node *x, int *status)
{
    CCpq_node *z, *zprev, *znext;

    if (x->type != PQ_QNODE) {
        return;
    }
    if (x == &T->pseudo_root) {
        return;
    }
    CCpq_set_FOREACH (x->children_set, z, children_elem, zprev, znext) {
        if (z->label != FULL(T)) {
            return;
        }
    }
    label_full (T, x);
    *status = PQSTATUS_MATCHED;
    return;
} /* END template_q1 */

static void template_q2 (CCpq_tree *T, CCpq_node *x, int real, int *status)
{
    CCpq_node *partial_child;
    int full_count;
    CCpq_node *z, *zprev, *znext;
    CCpq_node *full_end;

    if (x->type != PQ_QNODE) {
        return;
    }
    if (x == &T->pseudo_root) {
        return;
    }
    if (CCpq_set_SIZE (x->partial_children_set) > 1) {
        return;
    }
    if (left_end_is_full (T, x)) {
        full_end = CCpq_set_LEFT_ELEM (x->children_set);
    } else {
        full_end = CCpq_set_RIGHT_ELEM (x->children_set);
    }

    full_count = 0;
    z = full_end;
    zprev = (CCpq_node *) NULL;

    CCpq_set_FOREACH_FROM (z, children_elem, zprev, znext) {
        if (z->label == FULL(T)) {
            full_count++;
        } else {
            break;
        }
    }

    if (full_count != CCpq_set_SIZE (x->full_children_set)) {
        return;
    }
    if (!CCpq_set_ISEMPTY (x->partial_children_set) &&
        z != CCpq_set_LEFT_ELEM (x->partial_children_set)) {
        return;
    }
    assert (!CCpq_set_ISEMPTY (x->full_children_set) ||
            !CCpq_set_ISEMPTY (x->partial_children_set));

    if (!CCpq_set_ISEMPTY (x->partial_children_set)) {
        T->nontrivial = 1;
        if (real) {
            partial_child = CCpq_set_LEFT_ELEM (x->partial_children_set);
            CCpq_set_DELETE (partial_child, x->partial_children_set,
                           partial_children_elem);
            merge_qnode (T, partial_child);
        }
    }
    label_partial (T, x);
    *status = PQSTATUS_MATCHED;
    return;
} /* END template_q2 */

static int template_p2_root (CCpq_tree *T, CCpq_node *x, int real, int *status)
{
    int rval;

    if (x->type != PQ_PNODE) {
        return 0;
    }
    if (!CCpq_set_ISEMPTY (x->partial_children_set)) {
        return 0;
    }
    assert (!CCpq_set_ISEMPTY (x->full_children_set));
    assert (CCpq_set_SIZE (x->full_children_set) !=
            CCpq_set_SIZE (x->children_set));

    if (CCpq_set_SIZE (x->full_children_set) > 1) {
        T->nontrivial = 1;
    }

    if (real) {
        rval = collect_full_children (T, x, PQ_PNODE);
        if (rval) {
            fprintf (stderr, "collect_full_children failed\n");
            return 1;
        }
    }
    *status = PQSTATUS_MATCHED;
    return 0;
} /* END template_p2_ROOT */

static int template_p4_root (CCpq_tree *T, CCpq_node *x, int real, int *status)
{
    CCpq_node *partial_child;
    CCpq_node *full_child;
    int rval;

    if (x->type != PQ_PNODE) {
        return 0;
    }
    if (CCpq_set_SIZE (x->partial_children_set) != 1) {
        return 0;
    }

    partial_child = CCpq_set_LEFT_ELEM (x->partial_children_set);

    if (!CCpq_set_ISEMPTY (x->full_children_set)) {
        T->nontrivial = 1;
        if (real) {
            rval = collect_full_children (T, x, PQ_QNODE);
            if (rval) {
                fprintf (stderr, "collect_full_children failed\n");
                return 1;
            }
            full_child = CCpq_set_LEFT_ELEM (x->full_children_set);
            CCpq_set_DELETE (full_child, x->children_set, children_elem);
            CCpq_set_DELETE (full_child, x->full_children_set,
                           full_children_elem);
            full_child->parent = partial_child;
            full_child->parenttype = PQ_QNODE;
            if (left_end_is_full (T, partial_child)) {
                CCpq_set_ADD_LEFT (full_child, partial_child->children_set,
                                 children_elem);
                CCpq_set_ADD_LEFT (full_child,
                        partial_child->full_children_set, full_children_elem);
            } else {
                CCpq_set_ADD_RIGHT (full_child, partial_child->children_set,
                                  children_elem);
                CCpq_set_ADD_RIGHT (full_child,
                        partial_child->full_children_set, full_children_elem);
            }
            if (CCpq_set_SIZE (x->children_set) == 1) {
                full_child = CCpq_set_LEFT_ELEM (x->children_set);
                replace_node (x, full_child);
                PQ_node_free (&T->pqnode_world, x);
            }
        }
    }
    *status = PQSTATUS_MATCHED;
    return 0;
} /* END template_p4_ROOT */

static int template_p6_root (CCpq_tree *T, CCpq_node *x, int real, int *status)
{
    CCpq_node *partial_child1, *partial_child2;
    CCpq_node *full_end1, *full_end2, *empty_end2;
    CCpq_node *z, *zprev, *znext;
    CCpq_node *full_child;
    int rval;

    if (x->type != PQ_PNODE) {
        return 0;
    }
    if (CCpq_set_SIZE (x->partial_children_set) > 2) {
        return 0;
    }
    assert (CCpq_set_SIZE (x->partial_children_set) == 2);

    T->nontrivial = 1;

    if (real) {
        partial_child1 = CCpq_set_LEFT_ELEM (x->partial_children_set);
        partial_child2 = CCpq_set_RIGHT_ELEM (x->partial_children_set);

        rval = collect_full_children (T, x, PQ_QNODE);
        if (rval) {
            fprintf (stderr, "collect_full_children failed\n");
            return 1;
        }

        if (!CCpq_set_ISEMPTY (x->full_children_set)) {
            full_child = CCpq_set_LEFT_ELEM (x->full_children_set);
            CCpq_set_DELETE (full_child, x->children_set, children_elem);
            CCpq_set_DELETE (full_child, x->full_children_set,
                           full_children_elem);

            full_child->parent = partial_child1;
            full_child->parenttype = PQ_QNODE;
            if (left_end_is_full (T, partial_child1)) {
                CCpq_set_ADD_LEFT (full_child, partial_child1->children_set,
                                 children_elem);
                CCpq_set_ADD_LEFT (full_child,
                        partial_child1->full_children_set, full_children_elem);
            } else {
                CCpq_set_ADD_RIGHT (full_child, partial_child1->children_set,
                                  children_elem);
                CCpq_set_ADD_RIGHT (full_child,
                                  partial_child1->full_children_set,
                                  full_children_elem);
            }
        }
        if (left_end_is_full (T, partial_child1)) {
            full_end1 = CCpq_set_LEFT_ELEM (partial_child1->children_set);
        } else {
            full_end1 = CCpq_set_RIGHT_ELEM (partial_child1->children_set);
        }
        if (left_end_is_full (T, partial_child2)) {
            full_end2 = CCpq_set_LEFT_ELEM (partial_child2->children_set);
            empty_end2 = CCpq_set_RIGHT_ELEM (partial_child2->children_set);
        } else {
            full_end2 = CCpq_set_RIGHT_ELEM (partial_child2->children_set);
            empty_end2 = CCpq_set_LEFT_ELEM (partial_child2->children_set);
        }

        CCpq_set_PTR_REPLACE (full_end1->children_elem, 0, full_end2);
        CCpq_set_PTR_REPLACE (full_end2->children_elem, 0, full_end1);

        if (partial_child1->children_set.left == full_end1) {
            partial_child1->children_set.left = empty_end2;
        } else {
            partial_child1->children_set.right = empty_end2;
        }
        partial_child1->children_set.size += partial_child2->children_set.size;
        empty_end2->parent = partial_child1;

        CCpq_set_FOREACH_DEL (partial_child2->full_children_set, z,
                            full_children_elem, zprev, znext) {
            CCpq_set_DELETE (z, partial_child2->full_children_set,
                           full_children_elem);
            CCpq_set_ADD (z, partial_child1->full_children_set,
                        full_children_elem);
            z->parent = partial_child1;
        }

        CCpq_set_DELETE (partial_child2, x->children_set, children_elem);
        CCpq_set_DELETE (partial_child2, x->partial_children_set,
                       partial_children_elem);

        PQ_node_free (&T->pqnode_world, partial_child2);

        if (CCpq_set_SIZE (x->children_set) == 1) {
            full_child = CCpq_set_LEFT_ELEM (x->children_set);
            replace_node (x, full_child);
            PQ_node_free (&T->pqnode_world, x);
        }
    }
    *status = PQSTATUS_MATCHED;
    return 0;
} /* END template_p6_ROOT */

static void template_q3_root (CCpq_tree *T, CCpq_node *x, int real,
        int *status)
{
    int full_count;
    CCpq_node *z, *zprev, *znext, *zstart;

    /* note that we cannot talk about x->children_set, since x might be the
     * pseudo_root */

    if (x->type != PQ_QNODE) {
        return;
    }
    if (CCpq_set_ISEMPTY (x->partial_children_set)) {
        zstart = CCpq_set_LEFT_ELEM (x->full_children_set);
        assert (zstart->label == FULL(T));
        full_count = 1;
        zprev = zstart;
        z = zstart->children_elem.ptr1;

        CCpq_set_FOREACH_FROM (z, children_elem, zprev, znext) {
            if (z->label == FULL(T)) {
                full_count++;
            } else {
                break;
            }
        }

        zprev = zstart;
        z = zstart->children_elem.ptr2;

        CCpq_set_FOREACH_FROM (z, children_elem, zprev, znext) {
            if (z->label == FULL(T)) {
                full_count++;
            } else {
                break;
            }
        }
        if (full_count != CCpq_set_SIZE (x->full_children_set)) {
            return;
        }
    } else {
        zstart = CCpq_set_LEFT_ELEM (x->partial_children_set);
        z = zstart->children_elem.ptr1;
        if (!z || IS_EMPTY (z,T)) {
            z = zstart->children_elem.ptr2;
        }

        zprev = zstart;
        full_count = 0;
        CCpq_set_FOREACH_FROM (z, children_elem, zprev, znext) {
            if (z->label == FULL(T)) {
                full_count++;
            } else {
                break;
            }
        }
        if (full_count != CCpq_set_SIZE (x->full_children_set)) {
            return;
        }
        if (z && z->label == PARTIAL(T) &&
            CCpq_set_SIZE (x->partial_children_set) != 2) {
            return;
        }
        if ((!z || z->label != PARTIAL(T)) &&
            CCpq_set_SIZE (x->partial_children_set) != 1) {
            return;
        }
    }

    if (CCpq_set_SIZE (x->partial_children_set) > 0) {
        T->nontrivial = 1;

        if (real) {
            CCpq_set_FOREACH_DEL (x->partial_children_set, z,
                                partial_children_elem, zprev, znext) {
                CCpq_set_DELETE (z, x->partial_children_set,
                               partial_children_elem);
                merge_qnode (T, z);
            }
        }
    }
    *status = PQSTATUS_MATCHED;
    return;
} /* END template_q3_root */

/* label_full marks a node as full, and if it has a parent, adds it to the
 * parents list of full children */

static void label_full (CCpq_tree *T, CCpq_node *x)
{
    x->label = FULL(T);
    if (x->parent) {
        CCpq_set_ADD (x, x->parent->full_children_set, full_children_elem);
    }
} /* label_full */

/* label_partial marks a node as partial, and if it has a parent, adds it
 * to the parents list of partial children */

static void label_partial (CCpq_tree *T, CCpq_node *x)
{
    x->label = PARTIAL(T);
    if (x->parent) {
        CCpq_set_ADD (x, x->parent->partial_children_set,
                    partial_children_elem);
    }
} /* END label_partial */

/* collect_full_children collects all the full children of x, and places
 * them under a new P-node below x.  It sets the parenttype of the new node
 * to t */

static int collect_full_children (CCpq_tree *T, CCpq_node *x, int t)
{
    CCpq_node *new_node;
    CCpq_node *z, *zprev, *znext;

    assert (x->type == PQ_PNODE);
    if (CCpq_set_SIZE (x->full_children_set) > 1) {
        new_node = PQ_node_alloc (&T->pqnode_world);
        if (new_node == (CCpq_node *) NULL) {
            fprintf (stderr, "PQ_node_alloc failed\n");
            return 1;
        }
        node_init (T, new_node);
        new_node->label = FULL(T);
        new_node->type = PQ_PNODE;
        CCpq_set_FOREACH_DEL (x->full_children_set, z, full_children_elem,
                            zprev, znext) {
            CCpq_set_DELETE (z, x->children_set, children_elem);
            CCpq_set_DELETE (z, x->full_children_set, full_children_elem);
            z->parent = new_node;
            z->parenttype = PQ_PNODE;
            CCpq_set_ADD (z, new_node->children_set, children_elem);
            CCpq_set_ADD (z, new_node->full_children_set, full_children_elem);
        }
        assert (CCpq_set_ISEMPTY (x->full_children_set));
        new_node->parent = x;
        new_node->parenttype = t;
        CCpq_set_ADD (new_node, x->children_set, children_elem);
        CCpq_set_ADD (new_node, x->full_children_set, full_children_elem);
    } else if (CCpq_set_SIZE (x->full_children_set) == 1) {
        z = CCpq_set_LEFT_ELEM (x->full_children_set);
        z->parenttype = t;
    }
    return 0;
} /* collect_full_children */

/* collect_empty_children collects all the emtpy children of x, and places
 * them under a new P-node below x.  It sets the parenttype of the new node
 * to t.  If there is more than 1 empty child, it actually does this by
 * creating a new node to replace x, and moves all the full children of x to
 * the new node, and moves x down.  It returns a pointer to the node that is
 * where x used to be. */

static int collect_empty_children (CCpq_tree *T, CCpq_node **p_x, int t)
{
    CCpq_node *x = *p_x;
    CCpq_node *y;
    CCpq_node *z, *zprev, *znext;
    CCpq_node *new_node;

    assert (x->type == PQ_PNODE);
    assert (CCpq_set_ISEMPTY (x->partial_children_set));
    assert (CCpq_set_SIZE (x->children_set) >
            CCpq_set_SIZE (x->full_children_set));

    if (CCpq_set_SIZE (x->children_set) -
        CCpq_set_SIZE (x->full_children_set) > 1) {
        new_node = PQ_node_alloc (&T->pqnode_world);
        if (new_node == (CCpq_node *) NULL) {
            fprintf (stderr, "PQ_node_alloc failed\n");
            return 1;
        }

        node_init (T, new_node);
        new_node->label = x->label;
        new_node->type = t;

        CCpq_set_FOREACH_DEL (x->full_children_set, z, full_children_elem,
                            zprev, znext) {
            CCpq_set_DELETE (z, x->children_set, children_elem);
            CCpq_set_DELETE (z, x->full_children_set, full_children_elem);
            z->parent = new_node;
            z->parenttype = t;
            CCpq_set_ADD (z, new_node->children_set, children_elem);
            CCpq_set_ADD (z, new_node->full_children_set, full_children_elem);
        }

        y = x->parent;
        new_node->parent = y;
        new_node->parenttype = x->parenttype;
        new_node->children_elem.ptr1 = x->children_elem.ptr1;
        new_node->children_elem.ptr2 = x->children_elem.ptr2;
        if (y) {
            neighbor_replace (new_node->children_elem.ptr1,
                              &y->children_set, x, new_node);
            neighbor_replace (new_node->children_elem.ptr2,
                              &y->children_set, x, new_node);
        }
        x->parent = new_node;
        x->parenttype = t;

        x->label = EMPTY(T);
        CCpq_set_ADD (x, new_node->children_set, children_elem);

        x = new_node;
    } else if (CCpq_set_SIZE (x->children_set) -
               CCpq_set_SIZE (x->full_children_set) == 1) {
        CCpq_set_FOREACH (x->children_set, z, children_elem, zprev, znext) {
            if (IS_EMPTY (z,T)) {
                z->parenttype = t;
                break;
            }
        }
        x->type = t;
    }
    *p_x = x;
    return 0;
} /* END collect_empty_children */

/* neighbor_replace replaces the pointer to q in the children list with a
 * pointer to r */

static void neighbor_replace (CCpq_node *x, CCpq_set *s, CCpq_node *q,
        CCpq_node *r)
{
    if (x) {
        CCpq_set_PTR_REPLACE (x->children_elem, q, r);
    } else {
        if (s->left == q) {
            s->left = r;
        } else {
            assert (s->right == q);
            s->right = r;
        }
    }
} /* END neighbor_replace */

/* merge_qnode merges x into its parent */

static void merge_qnode (CCpq_tree *T, CCpq_node *x)
{
    CCpq_node *parent = x->parent;
    CCpq_node *full_end, *empty_end;
    CCpq_node *z, *zprev, *znext;
    CCpq_node *empty_neighbor, *full_neighbor;

    assert (x->type == PQ_QNODE);
    if (left_end_is_full (T, x)) {
        full_end = CCpq_set_LEFT_ELEM (x->children_set);
        empty_end = CCpq_set_RIGHT_ELEM (x->children_set);
    } else {
        full_end = CCpq_set_RIGHT_ELEM (x->children_set);
        empty_end = CCpq_set_LEFT_ELEM (x->children_set);
    }

    CCpq_set_FOREACH_DEL (x->full_children_set, z, full_children_elem,
                        zprev, znext) {
        CCpq_set_DELETE (z, x->full_children_set, full_children_elem);
        CCpq_set_ADD (z, parent->full_children_set, full_children_elem);
        z->parent = parent;
    }

    full_end->parent = parent;
    empty_end->parent = parent;
    empty_neighbor = x->children_elem.ptr1;
    full_neighbor = x->children_elem.ptr2;
    if ((empty_neighbor && !IS_EMPTY (empty_neighbor,T)) ||
        (full_neighbor && IS_EMPTY (full_neighbor,T))) {
        CC_SWAP (empty_neighbor, full_neighbor, z);
    }
    assert (full_end != empty_end);
    assert (!empty_neighbor || IS_EMPTY (empty_neighbor,T));
    assert (!full_neighbor || full_neighbor->label == FULL(T) ||
            full_neighbor->label == PARTIAL(T));
    assert (empty_neighbor || full_neighbor);

    CCpq_set_PTR_REPLACE (full_end->children_elem, 0, full_neighbor);
    CCpq_set_PTR_REPLACE (empty_end->children_elem, 0, empty_neighbor);
    neighbor_replace (full_neighbor, &parent->children_set, x, full_end);
    neighbor_replace (empty_neighbor, &parent->children_set, x, empty_end);

    parent->children_set.size += x->children_set.size - 1;

    PQ_node_free (&T->pqnode_world, x);
} /* END merge_qnode */

/* replace_node replaces old with new */

static void replace_node (CCpq_node *old, CCpq_node *new)
{
    new->parent = old->parent;
    new->parenttype = old->parenttype;

    new->children_elem.ptr1 = old->children_elem.ptr1;
    new->children_elem.ptr2 = old->children_elem.ptr2;

    if (new->parent) {
        neighbor_replace (new->children_elem.ptr1,
                          &new->parent->children_set, old, new);
        neighbor_replace (new->children_elem.ptr2,
                          &new->parent->children_set, old, new);
    } else {
        if (new->children_elem.ptr1) {
            CCpq_set_PTR_REPLACE (new->children_elem.ptr1->children_elem,
                                old, new);
        }
        if (new->children_elem.ptr2) {
            CCpq_set_PTR_REPLACE (new->children_elem.ptr2->children_elem,
                                old, new);
        }
    }
} /* END replace_node */

/* left_end_is_full is a boolean which indicates whether the left end of a
 * partial PQ_QNODE is full or not */

static int left_end_is_full (CCpq_tree *T, CCpq_node *x)
{
    CCpq_node *left_end;
    CCpq_node *right_end;

    assert (x->type == PQ_QNODE);

    left_end = CCpq_set_LEFT_ELEM (x->children_set);
    right_end = CCpq_set_RIGHT_ELEM (x->children_set);
    return IS_EMPTY (right_end,T) || left_end->label == FULL(T);
} /* END left_end_is_full */

static void node_init (CCpq_tree *T, CCpq_node *x)
{
    x->next = (CCpq_node *) NULL;
    x->number = T->node_counter;
    T->node_counter--;
    x->queue_elem.ptr1 = (CCpq_node *) NULL;
    x->queue_elem.ptr2 = (CCpq_node *) NULL;
    CCpq_set_INIT (x->children_set);
    CCpq_set_INIT (x->partial_children_set);
    CCpq_set_INIT (x->full_children_set);
    x->children_elem.ptr1 = (CCpq_node *) NULL;
    x->children_elem.ptr2 = (CCpq_node *) NULL;
    x->partial_children_elem.ptr1 = (CCpq_node *) NULL;
    x->partial_children_elem.ptr2 = (CCpq_node *) NULL;
    x->full_children_elem.ptr1 = (CCpq_node *) NULL;
    x->full_children_elem.ptr2 = (CCpq_node *) NULL;
    x->blocked_elem.ptr1 = (CCpq_node *) NULL;
    x->blocked_elem.ptr2 = (CCpq_node *) NULL;
    x->leaves_elem.ptr1 = (CCpq_node *) NULL;
    x->leaves_elem.ptr2 = (CCpq_node *) NULL;
    x->parent = (CCpq_node *) NULL;
    x->pertinent_child_count = 0;
    x->pertinent_leaf_count = 0;
    x->mark = 0;
    x->type = 0;
    x->parenttype = 0;
    x->label = 0;
} /* END node_init */

/* subtree_free recursively frees the tree rooted at x */

static void subtree_free (CCpq_node *x, CCptrworld *pqnode_world)
{
    CCpq_node *z, *zprev, *znext;

    CCpq_set_FOREACH_DEL (x->children_set, z, children_elem, zprev, znext) {
        CCpq_set_DELETE (z, x->children_set, children_elem);
        subtree_free (z, pqnode_world);
    }
    if (x->type != PQ_LEAF) {
        PQ_node_free (pqnode_world, x);
    }
} /* END subtree_free */

/* CCpq_find_root is only used by the dump and check routines, PQ_free_all,
   and pqtree_to_cuttree */

CCpq_node *CCpq_find_root (CCpq_tree *T)
{
    CCpq_node *zprev, *znext;
    CCpq_node *z;

    if (T->extern_node == 0) z = T->elems + 1;
    else z = T->elems;

    for (;;) {
        if (z->parenttype == PQ_PNODE) {
            if (z->parent)
                z = z->parent;
            else
                return z;
        } else {
            zprev = z->children_elem.ptr1;
            CCpq_set_FOREACH_FROM (z, children_elem, zprev, znext);
            if (zprev->parent)
                z = zprev->parent;
            else
                return zprev;
        }
    }
} /* END CCpq_find_root */

void CCpq_describe_solution (CCpq_tree *T)
{
    describe_subtree (CCpq_find_root (T));
    printf ("\n");
    fflush (stdout);
} /* END DESCRIBE_SOLUTION */

static void describe_subtree (CCpq_node *x)
{
    int lcnt;
    int ccnt;
    CCpq_node *z, *zprev, *znext;

    if (x->type == PQ_LEAF) return;
    if (x->type == PQ_PNODE) {
        printf ("(");
        lcnt = 0;
        ccnt = 0;
        CCpq_set_FOREACH (x->children_set, z, children_elem, zprev, znext) {
            if (z->type == PQ_LEAF) {
                lcnt++;
            } else {
                describe_subtree (z);
            }
            ccnt++;
        }
        if (lcnt) {
            printf ("%d-L ", lcnt);
        }
        printf (")<%d/%d> ", x->number, ccnt);
    } else {
        printf ("[");
        lcnt = 0;
        ccnt = 0;
        CCpq_set_FOREACH (x->children_set, z, children_elem, zprev, znext) {
            if (z->type == PQ_LEAF) {
                lcnt++;
            } else {
                if (lcnt) {
                    printf ("%d-L ", lcnt);
                    lcnt=0;
                }
                describe_subtree (z);
            }
            ccnt++;
        }
        if (lcnt) {
            printf ("%d-L ", lcnt);
            lcnt=0;
        }
        printf ("]<%d/%d> ", x->number, ccnt);
    }
}

void CCpq_dump_solution (CCpq_tree *T)
{
    dump_subtree (CCpq_find_root (T));
    printf ("\n");
    fflush (stdout);
} /* END DUMP_SOLUTION */

static void dump_subtree (CCpq_node *x)
{
    CCpq_node *z, *zprev, *znext;

    if (x->type == PQ_LEAF) {
        printf ("%d ", x->number);
    } else if (x->type == PQ_PNODE) {
        printf ("(");
        CCpq_set_FOREACH (x->children_set, z, children_elem, zprev, znext) {
            dump_subtree (z);
        }
        printf (")<%d> ", x->number);
    } else {
        printf ("[");
        CCpq_set_FOREACH (x->children_set, z, children_elem, zprev, znext) {
            dump_subtree (z);
        }
        printf ("]<%d> ", x->number);
    }
} /* END dump_subtree */

#ifndef NDEBUG

/* the various check routines are for debugging use only, to verify that the
 * state of the CCpq_tree is valid. */

static int check_tree (CCpq_tree *T)
{
    if (!check_subtree (T, CCpq_find_root (T))) {
        CCpq_dump_solution (T);
        return 0;
    } else {
        return 1;
    }
} /* END CHECK_TREE */

static int check_subtree (CCpq_tree *T, CCpq_node *x)
{
    CCpq_node *z, *zprev, *znext;

    if (!check_node (T, x)) {
        return 0;
    }
    CCpq_set_FOREACH (x->children_set, z, children_elem, zprev, znext) {
        if (!check_subtree (T, z)) {
            return 0;
        }
    }
    return 1;
} /* END CHECK_SUBTREE */

static int check_node (CCpq_tree *T, CCpq_node *x)
{
    CCpq_node *z, *zprev, *znext;
    int cnt;

    cnt = 0;

    if (x == &T->pseudo_root) {
        return 1;
    }
    CCpq_set_FOREACH (x->children_set, z, children_elem, zprev, znext) {
        cnt++;
        if (z->parenttype != x->type) {
            printf ("Node %d has parenttype %d, parent %d has type %d\n",
                    z->number, z->parenttype, x->number, x->type);
            return 0;
        }
        if (x->type == PQ_PNODE && z->parent != x) {
            printf ("Node %d has parent %d, should be P-node %d\n",
                    z->number, z->parent->number,x->number);
            return 0;
        }
        if (!((z->children_elem.ptr1 == zprev &&
               z->children_elem.ptr2 != zprev) ||
              (z->children_elem.ptr2 == zprev &&
               z->children_elem.ptr1 != zprev))) {
            printf ("Node %d has children ptr1 %d ptr2 %d, should be =, != %d\n",
                    z->number, z->children_elem.ptr1->number,
                    z->children_elem.ptr2->number, zprev->number);
            return 0;
        }
    }
    if (x->type == PQ_PNODE && cnt != CCpq_set_SIZE (x->children_set)) {
        printf ("P-node %d has %d children, set size %d\n",
                x->number, cnt, CCpq_set_SIZE (x->children_set));
        return 0;
    }
    if (x->type == PQ_QNODE &&
        CCpq_set_LEFT_ELEM (x->children_set)->parent != x) {
        printf ("Q_node %d left child %d has parent %d\n",
                x->number, CCpq_set_LEFT_ELEM (x->children_set)->number,
                CCpq_set_LEFT_ELEM (x->children_set)->parent->number);
        return 0;
    }
    if (x->type == PQ_QNODE &&
        CCpq_set_RIGHT_ELEM (x->children_set)->parent != x) {
        printf ("Q_node %d right child %d has parent %d\n",
                x->number, CCpq_set_RIGHT_ELEM (x->children_set)->number,
                CCpq_set_RIGHT_ELEM (x->children_set)->parent->number);
        return 0;
    }
    if (x->type != PQ_LEAF && CCpq_set_ISEMPTY (x->children_set)) {
        printf ("node %d is type %d, but has no children\n",
                x->number, x->type);
        return 0;
    }
    if (x->type == PQ_LEAF && !CCpq_set_ISEMPTY (x->children_set)) {
        printf ("leaf node %d has children\n", x->number);
        return 0;
    }
    if (x->type == PQ_PNODE && cnt < 2) {
        printf ("P-node %d only has %d children\n", x->number, cnt);
        return 0;
    }
    if (x->type == PQ_QNODE && cnt < ((x->label == PARTIAL(T)) ? 2 : 3)) {
        printf ("Q-node %d only has %d children\n", x->number, cnt);
        return 0;
    }
    return 1;
} /* END CHECK_NODE */

static int check_node_pert (CCpq_tree *T, CCpq_node *x)
{
    if (!check_node_pert_work (T, x)) {
        CCpq_dump_solution (T);
        return 0;
    } else {
        return 1;
    }
} /* END CHECK_NODE_PERT */

static int check_node_pert_work (CCpq_tree *T, CCpq_node *x)
{
    CCpq_node *z, *zprev, *znext;
    int cnt, full_cnt, partial_cnt, empty_cnt;

    if (!check_node (T, x)) {
        return 0;
    }
    if (x != &T->pseudo_root) {
        empty_cnt = 0;
        full_cnt = 0;
        partial_cnt = 0;
        CCpq_set_FOREACH (x->children_set, z, children_elem, zprev, znext) {
            if (IS_EMPTY (z,T))
                empty_cnt++;
            else if (z->label == PARTIAL(T))
                partial_cnt++;
            else if (z->label == FULL(T))
                full_cnt++;
            else {
                printf ("Node %d has label %d\n", z->number,
                        z->label);
                return 0;
            }
            if (!IS_EMPTY (z,T) && x != &T->pseudo_root && z->parent != x) {
                printf ("Node %d has parent %d, should be %d\n",
                        z->number, z->parent->number,
                        x->number);
                return 0;
            }
        }

        if (full_cnt != CCpq_set_SIZE (x->full_children_set)) {
            printf ("node %d has full children size %d, but has %d children full\n",
                    x->number, CCpq_set_SIZE (x->full_children_set),
                    full_cnt);
            return 0;
        }
        if (partial_cnt != CCpq_set_SIZE (x->partial_children_set)) {
            printf ("node %d has partial children size %d, but has %d children partial\n",
                    x->number,
                    CCpq_set_SIZE (x->partial_children_set),
                    partial_cnt);
            return 0;
        }
    }
    cnt = 0;
    CCpq_set_FOREACH (x->full_children_set, z, full_children_elem,
                    zprev, znext) {
        cnt++;
        if (!((z->full_children_elem.ptr1 == zprev &&
               z->full_children_elem.ptr2 != zprev) ||
              (z->full_children_elem.ptr2 == zprev &&
               z->full_children_elem.ptr1 != zprev) ||
              (z->full_children_elem.ptr1 == zprev &&
               z->full_children_elem.ptr2 == zprev &&
               zprev == (CCpq_node *) NULL))) {
            printf ("Node %d has full_children ptr1 %d ptr2 %d, should be =, != %d\n",
                    z->number,
                    z->full_children_elem.ptr1->number,
                    z->full_children_elem.ptr2->number,
                    zprev->number);
            return 0;
        }
        if (z->label != FULL(T)) {
            printf ("Node %d has label %d, but is in full set\n",
                    z->number, z->label);
            return 0;
        }
    }
    if (cnt != CCpq_set_SIZE (x->full_children_set)) {
        printf ("node %d has %d full children, set size %d\n",
                x->number, cnt,
                CCpq_set_SIZE (x->full_children_set));
        return 0;
    }
    cnt = 0;
    CCpq_set_FOREACH (x->partial_children_set, z, partial_children_elem,
                    zprev, znext) {
        cnt++;
        if (!((z->partial_children_elem.ptr1 == zprev &&
               z->partial_children_elem.ptr2 != zprev) ||
              (z->partial_children_elem.ptr2 == zprev &&
               z->partial_children_elem.ptr1 != zprev) ||
              (z->partial_children_elem.ptr1 == zprev &&
               z->partial_children_elem.ptr2 == zprev &&
               zprev == (CCpq_node *) NULL))) {
            printf ("Node %d has partial_children ptr1 %d ptr2 %d, should be =, != %d\n",
                    z->number,
                    z->partial_children_elem.ptr1->number,
                    z->partial_children_elem.ptr2->number,
                    zprev->number);
            return 0;
        }
        if (z->label != PARTIAL(T)) {
            printf ("Node %d has label %d, but is in partial set\n",
                    z->number, z->label);
            return 0;
        }
    }
    if (cnt != CCpq_set_SIZE (x->partial_children_set)) {
        printf ("node %d has %d partial children, set size %d\n",
                x->number, cnt,
                CCpq_set_SIZE (x->partial_children_set));
        return 0;
    }
    return 1;
} /* END CHECK_NODE_PERT_WORK */

#endif /* NDEBUG */

int CCpq_cuttree_to_pq (CCtsp_cuttree *ct, CCpq_tree *pqT)
{
    int i;
    int nodecount = ct->nodecount;
    int extern_node = ct->extern_node;
    CCpq_node *elems = (CCpq_node *) NULL;
    CCpq_node *root;
    CCtsp_cutnode *r;
    int rval;

    CCpq_tree_free (pqT);

    if (nodecount < 3) {
        fprintf (stderr, "Can't build PQ tree with %d nodes\n",
                 ct->nodecount);
        rval = 1; goto CLEANUP;
    }

    elems = CC_SAFE_MALLOC (nodecount, CCpq_node);

    if (elems == (CCpq_node *) NULL) {
        fprintf (stderr, "Out of memory in CCpq_tree_trivial\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<nodecount; i++) {
        elems[i].number = i;
        if (i == extern_node) {
            elems[i].type = PQ_EXTERN;
        } else {
            elems[i].parent = (CCpq_node *) NULL;
            elems[i].label = 0;
            elems[i].mark = 0;
            elems[i].type = PQ_LEAF;
            CCpq_set_INIT (elems[i].children_set);
            CCpq_set_INIT (elems[i].full_children_set);
            CCpq_set_INIT (elems[i].partial_children_set);
            elems[i].pertinent_child_count = 0;
            elems[i].pertinent_leaf_count = 0;
        }
    }

    r = ct->root->child;
    if (r->type == CCtsp_CUT_EXTERN) {
        r = r->sibling;
    }

    pqT->nodecount = nodecount;
    pqT->extern_node = extern_node;
    pqT->elems = elems;
    pqT->leaflist = (CCpq_node *) NULL;
    pqT->markbase = 0;
    pqT->node_counter = -1;
    pqT->nontrivial = 0;

    root = cuttree_to_pqtree_work (r, ct->nodelist, pqT);
    if (root == (CCpq_node *) NULL) {
        fprintf (stderr, "cuttree_to_pqtree_work failed\n");
        rval = 1; goto CLEANUP;
    }

    root->parenttype = PQ_PNODE;

    assert (check_tree (pqT));

    rval = 0;

  CLEANUP:
    if (rval) {
        CC_IFFREE (elems, CCpq_node);
        CCpq_tree_init (pqT);
    }
    return rval;
}

static CCpq_node *cuttree_to_pqtree_work (CCtsp_cutnode *x,
        CCtsp_cutnode *nodelist, CCpq_tree *T)
{
    CCpq_node *n = (CCpq_node *) NULL;
    CCpq_node *c;
    CCtsp_cutnode *p;

    if (x->type == CCtsp_CUT_LEAF) {
        n = &T->elems[x - nodelist];
        n->type = PQ_LEAF;
        return n;
    }

    n = PQ_node_alloc (&T->pqnode_world);
    if (n == (CCpq_node *) NULL) {
        fprintf (stderr, "Out of memory in cuttree_to_pqtree_work\n");
        goto FAILURE;
    }

    node_init (T, n);

    if (x->type == CCtsp_CUT_PNODE) {
        n->type = PQ_PNODE;
    } else if (x->type == CCtsp_CUT_QNODE) {
        n->type = PQ_QNODE;
    } else {
        fprintf (stderr, "Unknown node type %d\n", x->type);
        goto FAILURE;
    }

    for (p = x->child; p; p = p->sibling) {
        c = cuttree_to_pqtree_work (p, nodelist, T);
        if (c == (CCpq_node *) NULL) goto FAILURE;
        c->parenttype = n->type;
        c->parent = n;
        CCpq_set_ADD (c, n->children_set, children_elem);
    }

    return n;

  FAILURE:
    if (n) {
        subtree_free (n, &T->pqnode_world);
    }
    return (CCpq_node *) NULL;
}

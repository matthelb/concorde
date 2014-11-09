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
/*                             SHRINK ROUTINES                              */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook (TSP Code)             */
/*  Date: July 19, 1996                                                     */
/*        October 21, 1996 (bico)                                           */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  void CCcut_SRK_init_graph (CC_SRKgraph *G)                              */
/*    INITIALIZES the fields of the CC_SRKgraph.                            */
/*                                                                          */
/*  int CCcut_SRK_buildgraph (CC_SRKgraph *G, int ncount, int ecount,       */
/*      int *elist, double *dlen)                                           */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCcut_SRK_free_graph (CC_SRKgraph *G)                              */
/*    FREES the CC_SRKgraph.                                                */
/*                                                                          */
/*  void CCcut_SRK_init_expinfo (CC_SRKexpinfo *expand)                     */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCcut_SRK_free_expinfo (CC_SRKexpinfo *expand)                     */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCcut_SRK_init_callback (CC_SRKcallback *cb)                       */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCcut_SRK_grab_edges (CC_SRKgraph *G, int *oncount, int *oecount,   */
/*      int **olist, double **olen, CC_SRKexpinfo *expand)                  */
/*          int **omembers)                                                 */
/*    RETURNS the edges and member lists for the shrunk graph.              */
/*     -G is a pointer to a shrunk graph                                    */
/*     -oncount returns the number of nodes in the shrunk graph             */
/*     -oecount returns the number of edges in the shrunk graph             */
/*     -olist returns the edges in node node format                         */
/*     -olen returns the edge lengths                                       */
/*     -expand will be filled in with a memindex and members array;         */
/*      memindex returns pointers into the members array, the               */
/*      members of node i will be stored in from memindex[i] to             */
/*      memindex[i+1] - 1, so memindex is ncount + 1 long; members          */
/*      returns the nodes lists corresponding to each node in the           */
/*      shrunk graph. (expand can be NULL)                                  */
/*                                                                          */
/*  int CCcut_SRK_grab_nodes (CC_SRKgraph *G, CC_SRKexpinfo *expand)        */
/*    RETURNS the member lists for the shrunk graph (see above)             */
/*                                                                          */
/*  int CCcut_SRK_trivial (int ncount, CC_SRKexpinfo *expand)               */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCcut_SRK_expand (CC_SRKexpinfo *expand, int *arr, int size,        */
/*      int **pnewarr, int *pnewsize)                                       */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCcut_SRK_identify_paths (CC_SRKgraph *G, int *newcount,           */
/*      int onecnt_okay)                                                    */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCcut_SRK_identify_paths_to_edges (CC_SRKgraph *G, int *newcount,  */
/*      int onecnt_okay)                                                    */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCcut_SRK_identify_ones (CC_SRKgraph *G, int *count,               */
/*      double epsilon)                                                     */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCcut_SRK_identify_one_triangles (CC_SRKgraph *G, int *count,      */
/*      CC_SRKnode *qstart, double epsilon, double cutoff, int unmarked)    */
/*    SHRINKS any one edge that sits in a tight triangle.                   */
/*     -G is the current shrunk graph                                       */
/*     -count returns the number of shrunk triangles (can be NULL)          */
/*     -qstart can point to the start of a queue (linked by qnext)          */
/*     -epsilon is used to determine one edges (at least 1.0 - epsilon)     */
/*     -cutoff is used to determine tight triangles (weight cutoff)         */
/*     -unmarked should be nonzero if only unmarked nodes (determined       */
/*      by G->marker) should be involved in shrinks                         */
/*                                                                          */
/*  void CCcut_SRK_identify_tight_triangles (CC_SRKgraph *G, int *count,    */
/*      double cutoff, int unmarked)                                        */
/*    SHRINKS any tight triangle into a single node.                        */
/*     -G is the current shrunk graph                                       */
/*     -count returns the number of shrunk triangles (can be NULL)          */
/*     -cutoff is used to decide if a triangle is tight (shrunk any T       */
/*      with x(T) >= cutoff)                                                */
/*     -unmarked should be nonzero if only unmarked nodes (determined       */
/*      by G->marker) should be involved in shrinks                         */
/*    NOTES: All new shrunk nodes will be marked.                           */
/*                                                                          */
/*  void CCcut_SRK_identify_tight_squares (CC_SRKgraph *G, int *count,      */
/*      double cutoff, int unmarked)                                        */
/*    SHRINKS tight squares into a single nodes.                            */
/*      -see above.                                                         */
/*                                                                          */
/*  void CCcut_SRK_identify_triangle_square (CC_SRKgraph *G, int *count,    */
/*      double epsilon, int unmarked)                                       */
/*    SHRINKS tight triangles within tight squares.                         */
/*      -epsilon is used to determine the tight triangle and square.        */
/*                                                                          */
/*  void CCcut_SRK_identify_one_square (CC_SRKgraph *G, int *count,         */
/*      double epsilon, double cutoff, int unmarked)                        */
/*    SHRINKS two opposite 1-edge in a tight 4-square                       */
/*                                                                          */
/*  void CCcut_SRK_identify_nodes (CC_SRKgraph *G, CC_SRKnode *n,           */
/*      CC_SRKnode *m)                                                      */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCcut_SRK_identify_set (CC_SRKgraph *G, int scount, int *slist)     */
/*    SHRINKS the nodes int slist to a single node.  NOTE this works        */
/*    by matching the indices in slist to the order in G's linked list of   */
/*    nodes.  (So it will work immediately after a call to grab_edges       */
/*    if the set is based on the ends of the edges.)                        */
/*                                                                          */
/*  void CCcut_SRK_increment_marker (CC_SRKgraph *G)                        */
/*    INCREASES the field used to mark nodes by 1.                          */
/*                                                                          */
/*  int CCcut_SRK_subtour_shrink (CC_SRKgraph *G, double *minval,           */
/*      double epsilon, CC_SRKcallback *cb, int **cut, int *cutcount)       */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCcut_SRK_identify_pr_edges (CC_SRKgraph *G, double *minval,        */
/*      int *count, CC_SRKnode *qstart, double epsilon,                     */
/*      CC_SRKcallback *cb, int **cut, int *cutcount)                       */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCcut_SRK_defluff (CC_SRKGRAPH *G)                                  */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCcut_SRK_set_mark (CC_SRKgraph *G, int marker)                    */
/*    SETS the mark field of all active nodes to marker.                    */
/*                                                                          */
/*  int CCcut_SRK_original_ncount (CC_SRKexpinfo *expand)                   */
/*      RETURNS the number of nodes in the original (unshrunk) graph.       */
/*                                                                          */
/*    NOTES: Cyles of 1's will be shrunk into single nodes.                 */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "cut.h"

#define ADD_TO_PR_QUEUE(n) {                                             \
    if (!(n)->onqueue) {                                                 \
        (n)->qnext = (CC_SRKnode *) NULL;                                \
        if (qtail)                                                       \
            qtail->qnext = (n);                                          \
        else                                                             \
            qhead = (n);                                                 \
        qtail = (n);                                                     \
        (n)->onqueue = 1;                                                \
    }                                                                    \
}

#undef  PR_USE3   /* PR_USE3 and PR_USE4 cannot be defined in current code */
#undef  PR_USE4
#undef  DEBUG_SHRINK
#define SRK_ZERO_EPSILON (1e-10)


static void
#ifdef DEBUG_SHRINK
    printgraph (CC_SRKgraph *G),
#endif
    count_ones (CC_SRKgraph *G),
    merge_adj (CC_SRKgraph *G, CC_SRKnode *n, CC_SRKnode *m),
    initial_queue (CC_SRKgraph *G, CC_SRKnode **outqhead,
            CC_SRKnode **outqtail, int unmarked);

static int
    test_node (CC_SRKnode *n, double *minval, CC_SRKcallback *cb, int **cut,
            int *cutcount),
    expand_the_node (CC_SRKnode *n, int *cutcount, int **cut),
    expand_and_pass (CC_SRKnode *n, int (*doit_fn) (double, int, int *,
            void *), void *pass_param);


int CCcut_SRK_subtour_shrink (CC_SRKgraph *G, double *minval, double epsilon,
        CC_SRKcallback *cb, int **cut, int *cutcount)
{
    int rval = 0;
    int k;
    int cnt = 0;
    CC_SRKnode *n;

    for (n = G->head; n; n = n->next) {
        cnt++;
    }

    /* In the tsp, paths of 1's can be shrunk via a call to the function  */
    /* CCcut_SRK_identify_paths.                                          */

    /* Could call a version of CCcut_SRK_identify_ones */

    /* printf ("Identify PR edges ....\n"); fflush (stdout); */
    rval = CCcut_SRK_identify_pr_edges (G, minval, &k, (CC_SRKnode *) NULL,
                   epsilon, cb, cut, cutcount);
    if (rval) {
        fprintf (stderr, "CCcut_SRK_identify_pr_edges failed\n"); goto CLEANUP;
    }

    cnt -= k;
    /* printf ("Graph shrunk to %d nodes\n", cnt); fflush (stdout); */

CLEANUP:

    return rval;
}

static void count_ones (CC_SRKgraph *G)
{
    int k;
    CC_SRKnode *n;
    CC_SRKedge *e;

    for (n = G->head; n; n = n->next) {
        for (k = 0, e = n->adj; e; e = e->next) {
            if (e->weight == 1.0)
                k++;
        }
        n->onecnt = k;
    }
}

void CCcut_SRK_identify_paths (CC_SRKgraph *G, int *newcount, int onecnt_okay)
{
    CC_SRKnode *n, *m, *last;
    CC_SRKedge *e, *f;
    int dropcnt = 0;
    double dropweight = 0.0;
    int k;

    /* printf ("Identify paths ...\n"); fflush (stdout); */

    if (!onecnt_okay)
        count_ones (G);
    for (n = G->head; n; n = n->next) {
        if (n->onecnt == 1) {
            e = n->adj;
            while (e->weight != 1.0)
                e = e->next;
            last = n;
            m = e->end;
            while (m->onecnt != 1) {
                m->parent = n;
                m->members = n->members;
                n->members = m;
                e = m->adj;
                while (e->weight != 1.0 || e->end == last)
                    e = e->next;
                last = m;
                m = e->end;
            }
            m->parent = n;
            m->members = n->members;
            n->members = m;
            m->onecnt = 3;
        }
    }

    for (n = G->head; n->parent != n; n = n->next);
    G->head = n;
    G->head->prev = (CC_SRKnode *) NULL;
    for (n = G->head->next; n; n = n->next) {
        if (n->parent != n) {
            n->prev->next = n->next;
            if (n->next)
               n->next->prev = n->prev;
        }
    }
    for (k = 0, n = G->head; n; n = n->next) {
        k++;
        if (n->members) {
            for (e = n->members->adj; e; e = e->next) {
                e->other->end = n;
            }
            for (m = n->members->members; m; m = m->members) {
                for (e = m->adj; e; e = e->next) {
                    if (e->weight == 1.0) {
                        e->other->end = n;
                    } else {
                        /* drop fluff edge */
                        dropcnt++;
                        dropweight += e->weight;
                        f = e->other;
                        if (f->prev) {
                            f->prev->next = f->next;
                        } else {
                            e->end->adj = f->next;
                        }
                        if (f->next) {
                            f->next->prev = f->prev;
                        }
                    }
                }
            }
            n->weight = n->weight + n->members->weight;
            merge_adj (G, n, n->members);
        }
    }

    if (dropcnt > 0) {
        printf ("dropped %d edges of total weight %f\n", dropcnt, dropweight);
        fflush (stdout);
    }

    *newcount = k;
}

int CCcut_SRK_defluff (CC_SRKgraph *G)
{
    CC_SRKnode *n;
    CC_SRKedge *e, *enext;
    int k;
    int ndel = 0;
    double delweight = 0.0;

    for (n = G->head; n; n = n->next) {
        for (k = 0, e = n->adj; e; e = e->next) {
            if (e->weight >= 1.0 - SRK_ZERO_EPSILON) {
                e->weight = 1.0;
                k++;
            }
        }
        n->onecnt = k;
    }

    for (n = G->head; n; n = n->next) {
        for (e = n->adj, n->adj = (CC_SRKedge *) NULL; e; e = enext) {
            enext = e->next;
            if (e->weight == 1.0 ||
                (n->onecnt < 2 && e->end->onecnt < 2
                 && e->weight > SRK_ZERO_EPSILON)) {
                if (n->adj) n->adj->prev = e;
                e->next = n->adj;
                n->adj = e;
                e->prev = (CC_SRKedge *) NULL;
            } else {
                delweight += e->weight;
                ndel++;
            }
        }
    }

    if (ndel & 1) {
        fprintf (stderr, "Whoa, deleted %d (odd) endpoints in CCcut_SRK_defluff\n",
                 ndel);
        return -1;
    }
/*
    printf ("CCcut_SRK_defluff deleted %d endpoints (weight %.6f)\n", ndel,
            delweight);
    fflush (stdout);
*/
    return 0;
}

void CCcut_SRK_identify_paths_to_edges (CC_SRKgraph *G, int *newcount,
        int onecnt_okay)
{
    CC_SRKnode *n, *p, *m, *last;
    CC_SRKedge *e;
    int k;

    /* NOTE: We should add in code to drop fluff edges */

    if (!onecnt_okay)
        count_ones (G);
    for (n = G->head; n; n = n->next) {
        if (n->onecnt == 1) {
            e = n->adj;
            while (e->weight != 1.0)
                e = e->next;
            m = e->end;
            if (m->onecnt != 1) {
                e = m->adj;
                while (e->weight != 1.0 || e->end == n)
                    e = e->next;
                last = m;
                p = e->end;
                while (p->onecnt != 1) {
                    p->parent = m;
                    p->members = m->members;
                    m->members = p;
                    e = p->adj;
                    while (e->weight != 1.0 || e->end == last)
                         e = e->next;
                    last = p;
                    p = e->end;
                }
                p->parent = m;
                p->members = m->members;
                m->members = p;
                p->onecnt = 3;
                m->mark = G->marker;  /* indicate node was in a contraction */
            }
        }
    }


    for (n = G->head; n->parent != n; n = n->next);
    G->head = n;
    G->head->prev = (CC_SRKnode *) NULL;
    for (n = G->head->next; n; n = n->next) {
        if (n->parent != n) {
            n->prev->next = n->next;
            if (n->next)
               n->next->prev = n->prev;
        }
    }
    for (k = 0, n = G->head; n; n = n->next) {
        k++;
        if (n->members) {
            for (m = n->members; m; m = m->members) {
                for (e = m->adj; e; e = e->next)
                    e->other->end = n;
            }
            merge_adj (G, n, n->members);
        }
    }
    *newcount = k;
}

void CCcut_SRK_identify_ones (CC_SRKgraph *G, int *count, double epsilon)
{
    CC_SRKnode *n;
    CC_SRKedge *e;
    double tol = 1.0 - epsilon;

    /*  printf ("Identify ones ....\n"); fflush (stdout); */

    *count = 0;

    for (n = G->head; n; n = n->next) {
        do {
            for (e = n->adj; e; e = e->next) {
                if (e->weight >= tol) {
                    CCcut_SRK_identify_nodes (G, n, e->end);
                    (*count)++;
                    break;
                }
            }
        } while (e);
    }
}

int CCcut_SRK_crowder_padberg (CC_SRKgraph *G, double epsilon,
        CC_SRKcallback *cb)
{
    int rval = 0;
    CC_SRKnode *n;
    CC_SRKedge *e, *h;
    double tol, minval = 0.0;
    CC_SRKnode *qhead, *qtail;

    for (n = G->head; n->next; n = n->next) {
        n->qnext = n->next;
        n->onqueue = 1;
    }
    qhead = G->head;
    qtail = n;
    qtail->onqueue = 1;
    qtail->qnext = (CC_SRKnode *) NULL;

    while (qhead) {
        n = qhead;
        qhead = qhead->qnext;
        if (!qhead)
            qtail = (CC_SRKnode *) NULL;
        if (n->parent != n)
            continue;
        n->onqueue = 0;
        tol = 1.0 - epsilon;

        for (e = n->adj; e && e->weight < tol; e = e->next);
        if (e) {
            CCcut_SRK_identify_nodes (G, n, e->end);
            ADD_TO_PR_QUEUE(n);
            for (h = n->adj; h; h = h->next) {
                ADD_TO_PR_QUEUE(h->end);
            }
            rval = test_node (n, &minval, cb, (int **) NULL, (int *) NULL);
            if (rval) { fprintf (stderr, "test_node failed\n"); goto CLEANUP; }
        }
    }

CLEANUP:

    return rval;
}

int CCcut_SRK_identify_pr_edges (CC_SRKgraph *G, double *minval, int *count,
        CC_SRKnode *qstart, double epsilon, CC_SRKcallback *cb, int **cut,
        int *cutcount)
{
    int rval = 0;
    CC_SRKnode *n;
    CC_SRKedge *e, *f, *h;
    double tol, tol1, tol2;
    CC_SRKnode *qhead, *qtail;

    /* Trivial PR Test: If w(uv) >= c(u)/2, then we can shrink edge uv.   */

    /* First PR Test: If we have a triangle uv, uw, vw with               */
    /* w(uv) + w(vw) >= c(v)/2 and w(uw) + w(vw) >= c(w)/2, then we can   */
    /* shrink edge vw.                                                    */

    /* Second PR Test: Compute a lower bound on any cut that separates    */
    /* the ends of an edge by summing the minimum values of the edges to  */
    /* common neighbors of the ends. If the bound is >= "2", then we can  */
    /* shrink the edge.                                                   */

    /* Third PR Test: Edge uv with common neighbors xy. If                */
    /* w(ux) + w(uy) + w(uv) >= "1" and w(vx) + w(vy) + w(uv) >= "1" and  */
    /* at least one of w(uy) + w(uv) and w(vx) + w(uv) is >= "1" and      */
    /* at least one of w(ux) + w(uv) and w(vy) + w(uv) is >= "1" then we  */
    /* can shrink the edge uv.                                            */

    *count = 0;

    if (cut && !cutcount) {
        fprintf (stderr, "cut defined, but not cutcount\n");
        rval = 1; goto CLEANUP;
    }

    if (qstart) {
        qhead = qstart;
        for (n = qstart; n->next; n = n->next)
            n->onqueue = 1;
        qtail = n;
        qtail->onqueue = 1;
    } else {
        for (n = G->head; n->next; n = n->next) {
            n->qnext = n->next;
            n->onqueue = 1;
        }
        qhead = G->head;
        qtail = n;
        qtail->onqueue = 1;
        qtail->qnext = (CC_SRKnode *) NULL;
    }

    while (qhead) {
        n = qhead;
        qhead = qhead->qnext;
        if (!qhead)
            qtail = (CC_SRKnode *) NULL;
        if (n->parent != n)
            continue;
        n->onqueue = 0;
        tol = n->weight/2.0 - epsilon;

        for (e = n->adj; e && e->weight < tol; e = e->next);
        if (e) {
            rval = test_node (n, minval, cb, cut, cutcount);
            if (rval) { fprintf (stderr, "test_node failed\n"); goto CLEANUP; }
            rval = test_node (e->end, minval, cb, cut, cutcount);
            if (rval) { fprintf (stderr, "test_node failed\n"); goto CLEANUP; }
            CCcut_SRK_identify_nodes (G, n, e->end);
            (*count)++;
            ADD_TO_PR_QUEUE(n);
            for (h = n->adj; h; h = h->next) {
                ADD_TO_PR_QUEUE(h->end);
            }
        } else {
            for (e = n->adj; e; e = e->next)
                e->end->prweight = e->weight;
            for (e = n->adj; e; e = e->next) {
                tol1 = (n->weight/2.0) - e->weight - epsilon;
                tol2 = (e->end->weight/2.0) - e->weight - epsilon;
                for (f = e->end->adj; f; f = f->next) {
                    if (f->weight >= tol2 && f->end->prweight >= tol1) {
                        rval = test_node (n, minval, cb, cut, cutcount);
                        if (rval) {
                            fprintf (stderr, "test_node failed\n");
                            goto CLEANUP;
                        }
                        rval = test_node (e->end, minval, cb, cut, cutcount);
                        if (rval) {
                            fprintf (stderr, "test_node failed\n");
                            goto CLEANUP;
                        }
                        CCcut_SRK_identify_nodes (G, n, e->end);
                        (*count)++;
                        ADD_TO_PR_QUEUE(n);
                        for (h = n->adj; h; h = h->next) {
                            ADD_TO_PR_QUEUE(h->end);
                        }
                        goto GET_OUT;
                    }
                }
            }

#ifdef PR_USE3
    -Must modify to use node current min cut value.
            for (e = n->adj; e; e = e->next) {
                tol = e->weight;
                for (f = e->end->adj; f; f = f->next) {
                    if (f->end->prweight >= f->weight)
                        tol += f->weight;
                    else if (f->end->prweight > -CC_MINCUT_BIGDOUBLE)
                        tol += f->end->prweight;
                }
                if (tol >= 1.0 + onetol) {
                    printf ("X"); fflush (stdout);
                    CCcut_SRK_identify_nodes (G, n, e->end);
                    (*count)++;
                    ADD_TO_PR_QUEUE(n);
                    for (h = n->adj; h; h = h->next) {
                        ADD_TO_PR_QUEUE(h->end);
                    }
                    goto GET_OUT;
                }
            }
#endif

#ifdef PR_USE4
    -Must modify to use node weights.
            for (e = n->adj; e; e = e->next) {
                tol = onetol - e->weight;
                for (f = e->end->adj; f; f = f->next) {
                    if (f->end->prweight > -CC_MINCUT_BIGDOUBLE) {
                        for (h = f->next; h; h = h->next) {
                            if (h->end->prweight > -CC_MINCUT_BIGDOUBLE) {
                                if (f->weight + h->weight >= tol
                       && f->end->prweight + h->end->prweight >= tol
                       && (f->weight >= tol || h->end->prweight >= tol)
                       && (h->weight >= tol || f->end->prweight >= tol)) {
                                    CCcut_SRK_identify_nodes (G, n, e->end);
                                    (*count)++;
                                        ADD_TO_PR_QUEUE(n);
                                    for (h = n->adj; h; h = h->next) {
                                        ADD_TO_PR_QUEUE(h->end);
                                    }
                                    goto GET_OUT;
                                }
                            }
                        }
                    }
                }
            }
#endif

GET_OUT:
            for (e = n->adj; e; e = e->next)
                e->end->prweight = -CC_MINCUT_BIGDOUBLE;
        }
    }

CLEANUP:

    return rval;
}

static int test_node (CC_SRKnode *n, double *minval, CC_SRKcallback *cb,
        int **cut, int *cutcount)
{
    int rval = 0;

    if (n->weight < *minval) {
        *minval = n->weight;
        /* printf ("New minimum: %f\n", *minval); */
        if (cut) {
            CC_IFFREE (*cut, int);
            rval = expand_the_node (n, cutcount, cut);
            if (rval) {
                fprintf (stderr, "expand_the_node failed\n"); goto CLEANUP;
            }
        }
    }
    if (cb) {
        if (n->weight <= cb->cutoff) {
            rval = expand_and_pass (n, cb->doit_fn, cb->pass_param);
            if (rval) {
                fprintf (stderr,"expand_and_pass failed\n"); goto CLEANUP;
            }
        }
    }

CLEANUP:

    return rval;
}

static int expand_and_pass (CC_SRKnode *n, int (*doit_fn) (double, int, int *,
            void *), void *pass_param)
{
    int rval = 0;
    int cutcount;
    int *cut = (int *) NULL;

    if (!doit_fn) goto CLEANUP;

    rval = expand_the_node (n, &cutcount, &cut);
    if (rval) {
        fprintf (stderr, "expand_the_node failed\n"); fflush (stdout);
    }

    rval = doit_fn (n->weight, cutcount, cut, pass_param);
    if (rval) {
        fprintf (stderr, "doit_fn failed\n"); goto CLEANUP;
    }

CLEANUP:

    CC_IFFREE (cut, int);
    return rval;
}

static int expand_the_node (CC_SRKnode *n, int *cutcount, int **cut)
{
    int rval = 0;
    int cnt;
    int *tcut = (int *) NULL;
    CC_SRKnode *v;

    *cutcount = 0;
    *cut = (int *) NULL;

    cnt = 1;
    for (v = n->members; v; v = v->members) {
        cnt++;
    }
    tcut = CC_SAFE_MALLOC (cnt, int);
    if (!tcut) {
        fprintf (stderr, "out of memory in expand_the_node\n");
        rval = 1; goto CLEANUP;
    }

    tcut[0] = n->num;
    cnt = 1;
    for (v = n->members; v; v = v->members) {
        tcut[cnt++] = v->num;
    }

    *cutcount = cnt;
    *cut = tcut;

CLEANUP:

    return rval;
}

void CCcut_SRK_identify_one_triangles (CC_SRKgraph *G, int *count,
        CC_SRKnode *qstart, double epsilon, double cutoff, int unmarked)
{
    CC_SRKnode *n;
    CC_SRKedge *e, *f, *h;
    CC_SRKnode *qhead, *qtail;
    double tol = 1.0 - epsilon;
    double ctol = cutoff - 1.0 - epsilon;

    /* Identify any edge of weight one that is contained in a triangle of */
    /* weight >= cutoff */

    if (count) *count = 0;
    if (qstart) {
        qhead = qstart;
        for (n = qstart; n->next; n = n->next)
            n->onqueue = 1;
        qtail = n;
        qtail->onqueue = 1;
    } else {
        for (n = G->head; n->next; n = n->next) {
            n->qnext = n->next;
            n->onqueue = 1;
        }
        qhead = G->head;
        qtail = n;
        qtail->onqueue = 1;
        qtail->qnext = (CC_SRKnode *) NULL;
    }

    while (qhead) {
        n = qhead;
        qhead = qhead->qnext;
        if (!qhead)
            qtail = (CC_SRKnode *) NULL;
        if (n->parent != n || (unmarked && n->mark == G->marker))
            continue;
        n->onqueue = 0;

        for (e = n->adj; e && e->weight < tol; e = e->next);
        if (e) {
            for (f = n->adj; f; f = f->next)
                f->end->prweight = f->weight;
            for (f = e->end->adj; f; f = f->next) {
                if (unmarked && f->end->mark == G->marker)
                    continue;
                if (f->weight + f->end->prweight >= ctol ||
                   (f->weight >= ctol && f->end != n)) {
                    CCcut_SRK_identify_nodes (G, n, e->end);
                    if (count) (*count)++;
                    ADD_TO_PR_QUEUE(n);
                    for (h = n->adj; h; h = h->next) {
                        ADD_TO_PR_QUEUE(h->end);
                    }
                    n->mark = G->marker;
                    goto GET_OUT;
                }
            }
GET_OUT:
            for (e = n->adj; e; e = e->next)
                e->end->prweight = -CC_MINCUT_BIGDOUBLE;
        }
    }
}

void CCcut_SRK_identify_tight_triangles (CC_SRKgraph *G, int *count,
        double cutoff, int unmarked)
{
    CC_SRKnode *qhead, *qtail, *n, *m, *o;
    CC_SRKedge *e, *f, *h;
    double ctol;

    if (count) *count = 0;
    initial_queue (G, &qhead, &qtail, unmarked);

    while (qhead) {
        n = qhead;
        qhead = qhead->qnext;
        if (!qhead)
            qtail = (CC_SRKnode *) NULL;
        if (n->parent != n || (unmarked && n->mark == G->marker))
            continue;
        n->onqueue = 0;

        for (e = n->adj; e; e = e->next) {
            e->end->prweight = e->weight;
        }

        for (e = n->adj; e; e = e->next) {
            m = e->end;
            if (unmarked && m->mark == G->marker)
                continue;
            ctol = cutoff - e->weight;
            for (f = m->adj; f; f = f->next) {
                o = f->end;
                if (unmarked && o->mark == G->marker)
                    continue;
                if (f->weight + o->prweight >= ctol ||
                   (f->weight >= ctol && o != n)) {
                    CCcut_SRK_identify_nodes (G, n, m);
                    CCcut_SRK_identify_nodes (G, n, o);
                    if (count) (*count)++;
                    if (!unmarked) {
                        ADD_TO_PR_QUEUE(n);
                        for (h = n->adj; h; h = h->next) {
                            ADD_TO_PR_QUEUE(h->end);
                        }
                    }
                    n->mark = G->marker;
                    goto GET_OUT;
                }
            }
        }

GET_OUT:
        for (e = n->adj; e; e = e->next) {
            e->end->prweight = -CC_MINCUT_BIGDOUBLE;
        }
    }
}


void CCcut_SRK_identify_tight_squares (CC_SRKgraph *G, int *count,
        double cutoff, int unmarked)
{
    CC_SRKnode *qhead, *qtail, *n, *m, *o, *p;
    CC_SRKedge *e, *f, *h, *d;
    double ctol;

    if (count) *count = 0;
    initial_queue (G, &qhead, &qtail, unmarked);

    for (n = G->head; n; n = n->next) {
        n->prweight = 0.0;
    }

    while (qhead) {
        n = qhead;
        qhead = qhead->qnext;
        if (!qhead)
            qtail = (CC_SRKnode *) NULL;
        if (n->parent != n || (unmarked && n->mark == G->marker))
            continue;
        n->onqueue = 0;

        for (e = n->adj; e; e = e->next) {
            e->end->prweight = e->weight;
        }

        for (e = n->adj; e; e = e->next) {
            m = e->end;
            if (unmarked && m->mark == G->marker)
                continue;
            ctol = cutoff - e->weight;
            for (f = m->adj; f; f = f->next) {
                f->end->prweight += f->weight;
            }
            for (f = m->adj; f; f = f->next) {
                o = f->end;
                if (o == n || (unmarked && o->mark == G->marker))
                    continue;
                for (h = o->adj; h; h = h->next) {
                    p = h->end;
                    if (p == n || p == m || (unmarked && p->mark == G->marker))
                        continue;
                    if (p->prweight + o->prweight + h->weight >= ctol) {
                        if (count) (*count)++;
                        for (d = m->adj; d; d = d->next) {
                            d->end->prweight -= d->weight;
                        }
                        CCcut_SRK_identify_nodes (G, n, m);
                        CCcut_SRK_identify_nodes (G, n, o);
                        CCcut_SRK_identify_nodes (G, n, p);
                        if (!unmarked) {
                            ADD_TO_PR_QUEUE(n);
                            for (d = n->adj; d; d = d->next) {
                                ADD_TO_PR_QUEUE(d->end);
                            }
                        }
                        n->mark = G->marker;
                        goto GET_OUT;
                    }
                }
            }
            for (f = m->adj; f; f = f->next) {
                f->end->prweight -= f->weight;
            }
        }
GET_OUT:
        for (e = n->adj; e; e = e->next) {
            e->end->prweight = 0.0;
        }

    }

    for (n = G->head; n; n = n->next) {
        n->prweight = -CC_MINCUT_BIGDOUBLE;
    }
}

void CCcut_SRK_identify_triangle_square (CC_SRKgraph *G, int *count,
        double epsilon, int unmarked)
{
    CC_SRKnode *qhead, *qtail, *n, *m, *o, *p;
    CC_SRKedge *e, *f, *h;
    int shrinkit = 0;

    /* shrink tight triangles within tight squares */

    if (count) *count = 0;
    initial_queue (G, &qhead, &qtail, unmarked);

    for (n = G->head; n; n = n->next) {
        n->prweight = 0.0;
    }

    while (qhead) {
        n = qhead;
        qhead = qhead->qnext;
        if (!qhead)
            qtail = (CC_SRKnode *) NULL;
        if (n->parent != n || (unmarked && n->mark == G->marker))
            continue;
        n->onqueue = 0;

        for (e = n->adj; e; e = e->next) {
            e->end->prweight = e->weight;
        }

        for (e = n->adj; e; e = e->next) {
            m = e->end;
            if (unmarked && m->mark == G->marker)
                continue;
            for (f = m->adj; f; f = f->next) {
                f->end->prweight += f->weight;
            }
            for (f = m->adj; f; f = f->next) {
                o = f->end;
                if (o == n || (unmarked && o->mark == G->marker))
                    continue;
                if (e->weight + o->prweight >= 2.0 - epsilon) {
                    for (h = o->adj; h && !shrinkit; h = h->next) {
                        p = h->end;
                        if (p != n && p != m &&
                           (p->prweight + h->weight >= 1.0 - epsilon)) {
                            shrinkit = 1;
                        }
                    }
                    if (!shrinkit) {
                        for (h = m->adj; h && !shrinkit; h = h->next) {
                            p = h->end;
                            if (p != n && p != m &&
                               (p->prweight >= 1 - epsilon)) {
                                shrinkit = 1;
                            }
                        }
                    }
                    if (shrinkit) {
                        shrinkit = 0;
                        for (h = m->adj; h; h = h->next) {
                            h->end->prweight -= h->weight;
                        }
                        CCcut_SRK_identify_nodes (G, n, m);
                        CCcut_SRK_identify_nodes (G, n, o);
                        if (!unmarked) {
                            ADD_TO_PR_QUEUE(n);
                            for (h = n->adj; h; h = h->next) {
                                ADD_TO_PR_QUEUE(h->end);
                            }
                        }
                        n->mark = G->marker;
                        goto GET_OUT;
                    }
                }
            }
            for (f = m->adj; f; f = f->next) {
                f->end->prweight -= f->weight;
            }
        }
GET_OUT:
        for (e = n->adj; e; e = e->next) {
            e->end->prweight = 0.0;
        }

    }

    for (n = G->head; n; n = n->next) {
        n->prweight = -CC_MINCUT_BIGDOUBLE;
    }
}

void CCcut_SRK_identify_one_square (CC_SRKgraph *G, int *count,
        double epsilon, double cutoff, int unmarked)
{
    CC_SRKnode *qhead, *qtail, *n, *m, *o, *p;
    CC_SRKedge *e, *f, *h, *d;

    /* shrink the opposing 1's in a square where the edges between the 1's */
    /* have weight at least cutoff                                         */

    if (count) *count = 0;
    initial_queue (G, &qhead, &qtail, unmarked);

    for (n = G->head; n; n = n->next) {
        n->prweight = 0.0;
    }

    while (qhead) {
        n = qhead;
        qhead = qhead->qnext;
        if (!qhead)
            qtail = (CC_SRKnode *) NULL;
        if (n->parent != n || (unmarked && n->mark == G->marker))
            continue;
        n->onqueue = 0;

        for (e = n->adj; e && e->weight < 1.0 - epsilon; e = e->next);
        if (e) {
            m = e->end;
            if (!unmarked || m->mark != G->marker) {
                for (f = n->adj; f; f = f->next)
                    f->end->prweight = f->weight;
                for (f = m->adj; f; f = f->next)
                    f->end->prweight += f->weight;
                for (f = m->adj; f; f = f->next) {
                    o = f->end;
                    if (o == n || (unmarked && o->mark == G->marker))
                        continue;
                    for (h = o->adj; h && h->weight < 1.0 - epsilon;
                                                        h = h->next);
                    if (h) {
                        p = h->end;
                        if (p != n && p != m && 
                           (!unmarked || p->mark != G->marker) &&
                           (p->prweight + o->prweight >= cutoff)) {
                            CCcut_SRK_identify_nodes (G, n, m);
                            CCcut_SRK_identify_nodes (G, o, p);
                            if (count) (*count)++;
                            if (!unmarked) {
                                ADD_TO_PR_QUEUE(n);
                                for (d = n->adj; d; d = d->next) {
                                    ADD_TO_PR_QUEUE(d->end);
                                }
                                ADD_TO_PR_QUEUE(o);
                                for (d = o->adj; d; d = d->next) {
                                    ADD_TO_PR_QUEUE(d->end);
                                }
                            }
                            n->mark = G->marker;
                            o->mark = G->marker;
                            goto GET_OUT;
                        }
                    }
                }
GET_OUT:
                for (e = n->adj; e; e = e->next)
                    e->end->prweight = 0.0;
                for (e = m->adj; e; e = e->next)
                    e->end->prweight = 0.0;
            }
        }
    }

    for (n = G->head; n; n = n->next) {
        n->prweight = -CC_MINCUT_BIGDOUBLE;
    }
}

static void initial_queue (CC_SRKgraph *G, CC_SRKnode **outqhead,
       CC_SRKnode **outqtail, int unmarked)
{
    CC_SRKnode *n, *qhead, *qtail;

    if (unmarked) {
        qhead = (CC_SRKnode *) NULL;
        qtail = (CC_SRKnode *) NULL;
        for (n = G->head; n; n = n->next) {
            n->onqueue = 0;
            if (n->mark != G->marker) {
                ADD_TO_PR_QUEUE (n);
            }
        } 
    } else {
        for (n = G->head; n->next; n = n->next) {
            n->qnext = n->next;
            n->onqueue = 1;
        }
        qhead = G->head;
        qtail = n;
        qtail->onqueue = 1;
        qtail->qnext = (CC_SRKnode *) NULL;
    }

    *outqhead = qhead;
    *outqtail = qtail;
}

int CCcut_SRK_identify_set (CC_SRKgraph *G, int scount, int *slist) 
{
    int rval = 0;
    CC_SRKnode **nlist = (CC_SRKnode **) NULL;
    CC_SRKnode *n;
    int i, k, cnt;
    int *perm = (int *) NULL;

    nlist = CC_SAFE_MALLOC (scount, CC_SRKnode *);
    CCcheck_NULL (nlist, "out of memory for nlist");

    perm = CC_SAFE_MALLOC (scount, int);
    CCcheck_NULL (perm, "out of memory for perm");

    for (i = 0; i < scount; i++) perm[i] = i;
    CCutil_int_perm_quicksort (perm, slist, scount);

    k = 0;
    for (n = G->head, cnt = 0; n && k < scount; n = n->next, cnt++) {
        if (cnt == slist[perm[k]]) {
            nlist[perm[k]] = n;
            k++;
        }
    }

    if (k != scount) {
        fprintf (stderr, "Error - did not find all nodes in set\n");
        rval = 1;  goto CLEANUP;
    }

    for (i = 1; i < scount; i++) {
        CCcut_SRK_identify_nodes (G, nlist[0], nlist[i]);
    }

CLEANUP:

    CC_IFFREE (nlist, CC_SRKnode *);
    CC_IFFREE (perm, int);
    return rval;
}

void CCcut_SRK_identify_nodes (CC_SRKgraph *G, CC_SRKnode *n, CC_SRKnode *m)
{
    CC_SRKedge *e;

    m->parent = n;
    n->weight += m->weight;

    if (!n->members) {
        n->members = m;
    } else if (!m->members) {
        m->members = n->members;
        n->members = m;
    } else {
        CC_SRKnode *t;
        for (t = n->members; t->members; t = t->members);
        t->members = m;
    }

    for (e = m->adj; e; e = e->next) {
        e->other->end = n;
    }

    merge_adj (G, n, m);

    if (m->prev)
        m->prev->next = m->next;
    else
        G->head = m->next;
    if (m->next)
        m->next->prev = m->prev;
}

static void merge_adj (CC_SRKgraph *G, CC_SRKnode *n, CC_SRKnode *m)
{
    CC_SRKedge *e, *f, *last, *work;
    CC_SRKedge **hit = G->hit;

    /* String together the two lists */

    if (n->adj) {
        for (last = n->adj; last->next; last = last->next);
        last->next = m->adj;
        if (m->adj)
            m->adj->prev = last;
        work = n->adj;
    } else {
        work = m->adj;
    }

    /* Remove any edges that become loops */

    while (work && work->end == n) {
        n->weight -= work->weight;
        work = work->next;
    }

    if (work) {
        work->prev = (CC_SRKedge *) NULL;
        for (e = work->next; e; e = e->next) {
            if (e->end == n) {
                n->weight -= e->weight;
                e->prev->next = e->next;
                if (e->next)
                    e->next->prev = e->prev;
            }
        }
    } else {
        n->adj = (CC_SRKedge *) NULL;
        return;
    }

    /* Remove the duplicate edges in the work list */

    hit[work->end->num] = work;
    for (e = work->next; e; e = e->next) {
        if (hit[e->end->num]) {
            /* A duplicate edge */

            hit[e->end->num]->weight += e->weight;
            hit[e->end->num]->other->weight = hit[e->end->num]->weight;
            e->prev->next = e->next;
            if (e->next)
                e->next->prev = e->prev;

            /* Fix the adj of the other end of the duplicate */

            f = e->other;
            if (f->prev) {
                f->prev->next = f->next;
            } else {
                e->end->adj = f->next;
            }
            if (f->next) {
                f->next->prev = f->prev;
            }
        } else {
            hit[e->end->num] = e;
        }
    }

    for (e = work; e; e = e->next)
        hit[e->end->num] = (CC_SRKedge *) NULL;
    n->adj = work;
}

int CCcut_SRK_buildgraph (CC_SRKgraph *G, int ncount, int ecount, int *elist,
                    double *dlen)
{
    int i;
    int *degree = (int *) NULL;
    int newecount = 0;
    CC_SRKnode *nodespace, *n, *n1, *n2;
    CC_SRKedge *e, *adj1, *adj2;
    CC_SRKedge **hit;

    G->nodespace = CC_SAFE_MALLOC(ncount, CC_SRKnode);
    G->hit = CC_SAFE_MALLOC(ncount, CC_SRKedge *);
    if (!G->nodespace || !G->hit) {
        fprintf (stderr, "out of memory in CCcut_SRK_buildgraph\n");
        CC_IFFREE(G->nodespace, CC_SRKnode);
        CC_IFFREE(G->hit, CC_SRKedge *);
        return 1;
    }
    nodespace = G->nodespace;
    hit = G->hit;
    G->head = nodespace;
    G->original_ncount = ncount;
    G->original_ecount = ecount;
    G->marker          = 0;

    degree = CC_SAFE_MALLOC(ncount, int);
    if (!degree) {
        fprintf (stderr, "out of memory in CCcut_SRK_buildgraph\n");
        CC_IFFREE(G->nodespace, CC_SRKnode);
        CC_IFFREE(G->hit, CC_SRKedge *);
        return 1;
    }

    for (i = 0, n = nodespace; i < ncount; i++, n++) {
        n->prev = n - 1;
        n->next = n + 1;
        n->num = i;
        n->members = (CC_SRKnode *) NULL;
        n->parent = n;
        n->prweight = -CC_MINCUT_BIGDOUBLE;
        n->weight = 0.0;
        hit[i] = (CC_SRKedge *) NULL;
        degree[i] = 0;
        n->onecnt = 0;
        n->mark   = 0;
    }
    nodespace[0].prev = (CC_SRKnode *) NULL;
    nodespace[ncount - 1].next = (CC_SRKnode *) NULL;

    for (i = 0; i < ecount; i++) {
        if (dlen[i] > SRK_ZERO_EPSILON) {
            newecount++;
            degree[elist[2*i]]++;
            degree[elist[2*i+1]]++;
        }
    }
    G->edgespace = CC_SAFE_MALLOC(2*newecount, CC_SRKedge);
    if (!G->edgespace) {
        fprintf (stderr, "out of memory in CCcut_SRK_buildgraph\n");
        CC_IFFREE(G->nodespace, CC_SRKnode);
        CC_IFFREE(G->hit, CC_SRKedge *);
        return 1;
    }

    for (e = G->edgespace, i = 0; i < ncount; i++) {
        nodespace[i].adj = e;
        e += degree[i];
    }
    for (i = 0; i < ecount; i++) {
        if (dlen[i] > SRK_ZERO_EPSILON) {
            n1 = nodespace + elist[2 * i];
            n2 = nodespace + elist[2 * i + 1];
            adj1 = n1->adj;
            adj2 = n2->adj;
            adj1->end = n2;
            adj1->weight = dlen[i];
            adj1->next = adj1 + 1;
            adj1->prev = adj1 - 1;
            adj1->other = adj2;
            adj2->end = n1;
            adj2->weight = dlen[i];
            adj2->next = adj2 + 1;
            adj2->prev = adj2 - 1;
            adj2->other = adj1;
            n1->adj++;
            n2->adj++;
            if (dlen[i] == 1.0) {
                n1->onecnt++;
                n2->onecnt++;
            }
        }
    }
    for (e = G->edgespace, i = 0; i < ncount; i++) {
        if (degree[i]) {
            (nodespace[i].adj - 1)->next = (CC_SRKedge *) NULL;
            nodespace[i].adj = e;
            nodespace[i].adj->prev = (CC_SRKedge *) NULL;
            e += degree[i];
        } else {
            nodespace[i].adj = (CC_SRKedge *) NULL;
        }
    }

    for (i = 0, n = nodespace; i < ncount; i++, n++) {
        for (e = n->adj; e; e = e->next) {
            n->weight += e->weight;
        }
    }

    CC_IFFREE(degree, int);
    return 0;
}

int CCcut_SRK_grab_edges (CC_SRKgraph *G, int *oncount, int *oecount,
        int **olist, double **olen, CC_SRKexpinfo *expand)
{
    int rval = 0;
    int k, num;
    int ncount = 0, ecount = 0;
    CC_SRKnode *n;
    CC_SRKedge *e;

    *oncount = 0;
    *oecount = 0;
    *olist = (int *) NULL;
    *olen = (double *) NULL;
    if (expand) {
        CCcut_SRK_init_expinfo (expand);
    }

    for (n = G->head; n; n = n->next) {
        n->newnum = ncount;
        for (e = n->adj; e; e = e->next)
            ecount++;
        ncount++;
    }

    if (ecount % 2) {
        fprintf (stderr, "Error in grab_edges\n");
        rval = 1; goto CLEANUP;
    } else {
        ecount /= 2;
    }

    if (ecount == 0) {
        rval = 0; goto CLEANUP;
    }

    *olist = CC_SAFE_MALLOC (ecount * 2, int);
    *olen  = CC_SAFE_MALLOC (ecount, double);
    if (!(*olist) || !(*olen)) {
        fprintf (stderr, "out of memory in grab_edges\n");
        rval = 1; goto CLEANUP;
    }

    k = 0;
    for (n = G->head; n; n = n->next) {
        num = n->newnum;
        for (e = n->adj; e; e = e->next) {
            if (num < e->end->newnum) {
                (*olist)[2 * k] = num;
                (*olist)[2 * k + 1] = e->end->newnum;
                (*olen)[k++] = e->weight;
            }
        }
    }
    if (k != ecount) {
        fprintf (stderr, "Error in grab_edges\n");
        rval = 1; goto CLEANUP;
    }

    *oncount = ncount;
    *oecount = ecount;

    if (expand) {
        rval = CCcut_SRK_grab_nodes (G, expand);
        if (rval) {
            fprintf (stderr, "CCcut_SRK_grab_nodes failed\n"); goto CLEANUP;
        }
    }

CLEANUP:

    if (rval) {
        CC_IFFREE (*olist, int);
        CC_IFFREE (*olen, double);
        if (expand) {
            CCcut_SRK_free_expinfo (expand);
        }
    }

    return rval;
}

int CCcut_SRK_grab_nodes (CC_SRKgraph *G, CC_SRKexpinfo *expand)
{
    int rval = 0;
    int i;
    CC_SRKnode *n, *m;
    int total  = 0;
    int ncount = 0;

    if (!expand) {
        fprintf (stderr, "CCcut_SRK_grab_nodes called without an expand struct\n");
        rval = 1; goto CLEANUP;
    }

    for (n = G->head; n; n = n->next) {
        ncount++;
    }

    CCcut_SRK_init_expinfo (expand);
    expand->ncount   = ncount;
    expand->members  = CC_SAFE_MALLOC (G->original_ncount, int);
    expand->memindex = CC_SAFE_MALLOC (ncount + 1, int);
    if (!(expand->members) || !(expand->memindex)) {
        fprintf (stderr, "out of memory in grab_nodes\n");
        rval = 1; goto CLEANUP;
    }
    for (n = G->head, i = 0; n; n = n->next, i++) {
        expand->memindex[i] = total;
        expand->members[total++] = n->num;
        for (m = n->members; m; m = m->members)
            expand->members[total++] = m->num;
    }
    expand->memindex[i] = total;

CLEANUP:

    if (rval) CCcut_SRK_free_expinfo (expand);
    return rval;
}

void CCcut_SRK_init_graph (CC_SRKgraph *G)
{
    if (G) {
        G->nodespace = (CC_SRKnode *) NULL;
        G->edgespace = (CC_SRKedge *) NULL;
        G->head      = (CC_SRKnode *) NULL;
        G->hit       = (CC_SRKedge **) NULL;
        G->original_ncount = 0;
        G->original_ecount = 0;
        G->marker          = 0;
    }
}

void CCcut_SRK_free_graph (CC_SRKgraph *G)
{
    if (G) {
        CC_IFFREE(G->nodespace, CC_SRKnode);
        CC_IFFREE(G->edgespace, CC_SRKedge);
        CC_IFFREE(G->hit, CC_SRKedge *);
    }
}

void CCcut_SRK_increment_marker (CC_SRKgraph *G)
{
    if (G) {
        G->marker++;
    }
}

void CCcut_SRK_set_mark (CC_SRKgraph *G, int marker)
{
    if (G) {
        CC_SRKnode *n;

        G->marker = marker;
        for (n = G->head; n; n = n->next) {
            n->mark = marker;
        }
    }
}

void CCcut_SRK_init_callback (CC_SRKcallback *cb)
{
    if (cb) {
        cb->cutoff     = 0.0;
        cb->pass_param = (void *) NULL;
        cb->doit_fn    = (int (*) (double, int, int *, void *)) NULL;
    }
}

int CCcut_SRK_trivial (int ncount, CC_SRKexpinfo *expand)
{
    int i;

    CCcut_SRK_init_expinfo (expand);
    expand->memindex = CC_SAFE_MALLOC (ncount+1, int);
    if (!expand->memindex) {
        fprintf (stderr, "Out of memory in CCcut_SRK_trivial\n");
        return -1;
    }
    expand->members = CC_SAFE_MALLOC (ncount, int);
    if (!expand->members) {
        fprintf (stderr, "Out of memory in CCcut_SRK_trivial\n");
        CC_FREE (expand->memindex, int);
        return -1;
    }
    for (i=0; i<ncount; i++) {
        expand->members[i] = i;
        expand->memindex[i] = i;
    }
    expand->memindex[ncount] = ncount;
    expand->ncount = ncount;
    return 0;
}

void CCcut_SRK_init_expinfo (CC_SRKexpinfo *expand)
{
    expand->ncount   = 0;
    expand->memindex = (int *) NULL;
    expand->members  = (int *) NULL;
}

void CCcut_SRK_free_expinfo (CC_SRKexpinfo *expand)
{
    expand->ncount = 0;
    CC_IFFREE (expand->memindex, int);
    CC_IFFREE (expand->members, int);
}

int CCcut_SRK_expand (CC_SRKexpinfo *expand, int *arr, int size, int **pnewarr,
                int *pnewsize)
{
    int newsize = 0;
    int *newarr = (int *) NULL;
    int i, j, jend;

    *pnewsize = 0;
    *pnewarr = (int *) NULL;
    for (i=0; i<size; i++) {
        newsize += expand->memindex[arr[i]+1] - expand->memindex[arr[i]];
    }
    newarr = CC_SAFE_MALLOC (newsize, int);
    if (!newarr) {
        fprintf (stderr, "Out of memory in CCcut_SRK_expand\n");
        return -1;
    }
    newsize = 0;
    for (i=0; i<size; i++) {
        for (j=expand->memindex[arr[i]], jend = expand->memindex[arr[i]+1];
             j < jend; j++) {
            newarr[newsize++] = expand->members[j];
        }
    }
    *pnewarr = newarr;
    *pnewsize = newsize;
    return 0;
}

int CCcut_SRK_original_ncount (CC_SRKexpinfo *expand)
{
    return expand->memindex[expand->ncount];
}

#ifdef DEBUG_SHRINK
static void printgraph (CC_SRKgraph *G)
{
    CC_SRKnode *n;
    CC_SRKedge *e;

    for (n = G->head; n; n = n->next) {
        printf ("Node %d: ", n->num);
        fflush (stdout);
        for (e = n->adj; e; e = e->next) {
            printf ("%d [%.2f] ", e->end->num, e->weight);
            fflush (stdout);
            if (e->other->end != n || e->other->weight != e->weight) {
                printf ("(Whoops) ");
                fflush (stdout);
            }
        }
        printf ("\n");
    }
}
#endif

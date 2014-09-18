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
/*                 MIN WEIGHT FRACTIONAL 2-MATCHINGS                        */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 24, 1995 - Modified 10/4/95 (Bico)                       */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCfmatch_fractional_2match (int ncount, int ecount, int *elist,     */
/*      int *elen, CCdatagroup *dat, double *val, int *thematching,         */
/*      int *thedual, int *thebasis, int wantbasic, int silent,             */
/*      CCrandstate *rstate)                                                */
/*       int ncount (the number of nodes in the graph)                      */
/*       int ecount (the number of edges)                                   */
/*       int *elist (the edgelist in end1 end2 format)                      */
/*       int *elen (the weights on the edges)                               */
/*       CCdatagroup *dat (the info to price edges - NULL if no pricing)    */
/*       double *xcoord (the x-coordinates for geometric problems - this    */
/*                       field can be NULL)                                 */
/*       double *ycoord (the y-coordinates)                                 */
/*       int innorm (the NORM for pricing the complete edgeset)             */
/*       double *val (returns the optimal weight)                           */
/*       int *thematching (if non-NULL, then returns the optimal matching   */
/*                         in end1 end2 value format, where value is 1 if   */
/*                         edge gets assigned 0.5 and value is 2 if edge    */
/*                         gets 1.0 - note that the array should be         */
/*                         allocated by the calling routine, and should     */
/*                         be 6 * ncount + 1 long - it is terminated by a   */
/*                         -1)                                              */
/*       int *thedual (if non-NULL, then returns the optimal dual solution  */
/*                     with values twice their actual value (so they will   */
/*                     be integers - the array should be alloced by the     */
/*                     calling routine, and should be ncount long)          */
/*       int *thebasis (if non-NULL, then returns the edges in the optimal  */
/*                      basis in end1 end2 format)                          */
/*       int wantbasis (if nonzero, then the optimal basic solution will    */
/*                      be found)                                           */
/*       int silent (if nonzero, will suppress print messages)              */
/*                                                                          */
/*    NOTES:                                                                */
/*       Use to find an optimal basis for the initial tsp LP. By changing   */
/*    MATCHDEGREE to 1, it should find min-wieght fractional matchings.     */
/*       The nodes should be numbered from 0 up to ncount - 1. If dat       */
/*    is specified, then the code will use a price-repair loop to solve     */
/*    the problem over the complete graph.                                  */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "fmatch.h"
#include "kdtree.h"
#include "util.h"
#include "macrorus.h"

#define MATCHDEGREE 2    /* Only set up for 1 or 2? */
#define MAXWEIGHT 1000000000

#define ZERO ((unsigned char) 0)
#define ONE ((unsigned char) 1)
#define TWO ((unsigned char) 2)
#ifdef FALSE
#undef FALSE
#endif
#define FALSE ((unsigned char) 0)
#ifdef TRUE
#undef TRUE
#endif
#define TRUE ((unsigned char) 1)

#define OTHEREND(e,n) (e->ends[0] == n ? e->ends[1] : e->ends[0])

typedef struct edgeptr {
    struct edge *this;
    struct node *other;
    struct edgeptr *next;
} edgeptr;

typedef struct edge {
    struct edge *next;
    struct edge *pnext;
    struct node *ends[2];
    int weight;
    int z;
    unsigned char x;
    unsigned char basic;
} edge;

typedef struct shortedge {
    struct node *ends[2];
    struct shortedge *next;
} shortedge;

typedef struct node {
    struct node *next;
    struct node *pnext;
    edgeptr *adj;
    edge *parentedge;
    int name;
    int y;
    int label;
    struct {
        int order;
        struct node *next;
        struct node **prev;
    } sort;
    unsigned char flag;
    unsigned char matchcnt;
} node;

typedef struct graph {
    int ncount;
    node *nodelist;
    edge *edgelist;
    node **nodenames;
    int PLUS;
    int MINUS;
    CCptrworld node_world;
    CCptrworld edge_world;
    CCptrworld edgeptr_world;
    CCptrworld shortedge_world;
} graph;



static void
    init_graph (graph *G),
    free_graph (graph *G),
    basic_grab_basic (node *n, int parity, int PLUS, int MINUS),
    basic_mark_component_as_done (node *n),
    basic_expand (node *n, int *hit_odd_circuit, int PLUS, int MINUS),
    basic_minalpha (node *n, node **new, int *alpha, int flip_plus_and_minus,
        int PLUS, int MINUS),
    basic_subalpha (node *n, int alpha, int flip_plus_and_minus, int PLUS,
        int MINUS),
    initmat (graph *G),
    initlist (graph *G, node *head, node *tail, node *head2, node *tail2,
       CCdatagroup *dat),
    expand (node *n, int *found_one, int PLUS, int MINUS),
    minalpha (node *n, node **new, int *alpha, int PLUS, int MINUS),
    subalpha (node *n, int alpha, int PLUS, int MINUS),
    augmentpath (node *n),
    augmentpath1 (node *n),
    augmentplushalf (node *n, edge *e),
    augmentminushalf (node *n, edge *e),
    flipcycle (node *n, edge *e, unsigned char v),
    augmentpair (node *n, node *n1, edge *e),
    setflag (node *n, unsigned char v),
    flipseg (node *n, node *n2),
    stringup (node *n, node *n1, node *n2, edge *e),
    y_quicksort (node **list, int *y, int l, int u),
    free_linked_world (graph *G);

static int
    buildgraph (graph *G, int ncount, int ecount, int *elist, int *elen),
    basicrun (graph *G),
    basicscan (graph *G, node *n),
    basic_check_scan (graph *G, node *n),
    basicgrow (graph *G, node *n),
    basic_grab_ones (node *n, int parity, edge **odd_circuit, int PLUS,
        int MINUS),
    basic_checkout_basic (node *n, int parity, edge **odd_circuit, int PLUS,
        int MINUS),
    twomatch (graph *G),
    chkmat (graph *G, double *val),
    fixmatch (graph *G, int *radded, CCdatagroup *dat, CCrandstate *rstate),
    kd_fixmatch (graph *G, int *radded, CCdatagroup *dat, CCrandstate *rstate),
    x_fixmatch (graph *G, int *radded, CCdatagroup *dat),
    junk_fixmatch (graph *G, int *radded, CCdatagroup *dat),
    checkoutedge (graph *G, node *n1, node *n2, int *hit, CCdatagroup *dat),
    precheckoutedge (node *n1, node *n2, shortedge **list, CCdatagroup *dat,
        CCptrworld *shortedge_world),
    addbadedge (graph *G, node *n1, node *n2, int w),
    augment (graph *G, node *n);

static node
    *basic_dualchange (node *n, int PLUS, int MINUS),
    *dualchange (node *n, int PLUS, int MINUS),
    *findflag (node *n),
    *findhole (node *n, node *n2);

static edge
    *newedge (graph *G, node *n1, node *n2),
    *findedge (node *n1, node *n2);




/********** Allocation routines **********/


CC_PTRWORLD_ROUTINES (node, nodealloc, node_bulk_alloc, nodefree)
CC_PTRWORLD_LEAKS_ROUTINE (node, node_check_leaks, flag, unsigned char)

CC_PTRWORLD_ROUTINES (edge, edgealloc, edge_bulk_alloc, edgefree)
CC_PTRWORLD_LISTFREE_ROUTINE (edge, edge_listfree, edgefree)
CC_PTRWORLD_LEAKS_ROUTINE (edge, edge_check_leaks, basic, unsigned char)

CC_PTRWORLD_ROUTINES (edgeptr, edgeptralloc, edgeptr_bulk_alloc, edgeptrfree)
CC_PTRWORLD_LISTFREE_ROUTINE (edgeptr, edgeptr_listfree, edgeptrfree)
CC_PTRWORLD_LEAKS_ROUTINE (edgeptr, edgeptr_check_leaks, this, edge *)

CC_PTRWORLD_ROUTINES (shortedge, shortedgealloc, shortedge_bulk_alloc,
        shortedgefree)
CC_PTRWORLD_LEAKS_ROUTINE (shortedge, shortedge_check_leaks, ends[0], node *)


static void free_linked_world (graph *G)
{
    int ntotal;
    int nreserve;
    if (node_check_leaks (&G->node_world, &ntotal, &nreserve)) {
        fprintf (stderr, "WARNING: Outstanding nodes %d (total %d)\n",
                 ntotal - nreserve, ntotal);
    }
    CCptrworld_delete (&G->node_world);

    if (edge_check_leaks (&G->edge_world, &ntotal, &nreserve)) {
        fprintf (stderr, "WARNING: Outstanding edges %d (total %d)\n",
                 ntotal - nreserve, ntotal);
    }
    CCptrworld_delete (&G->edge_world);

    if (edgeptr_check_leaks (&G->edgeptr_world, &ntotal, &nreserve)) {
        fprintf (stderr, "WARNING: Outstanding edgeptrs %d (total %d)\n",
                 ntotal - nreserve, ntotal);
    }
    CCptrworld_delete (&G->edgeptr_world);

    if (shortedge_check_leaks (&G->shortedge_world, &ntotal, &nreserve)) {
        fprintf (stderr, "WARNING: Outstanding shortedges %d (total %d)\n",
                 ntotal - nreserve, ntotal);
    }
    CCptrworld_delete (&G->shortedge_world);
}


/********** main routines **********/


static int buildgraph (graph *G, int ncount, int ecount, int *elist, int *elen)
{
    int i;
    node *n;
    edge *e;

    init_graph (G);

    G->ncount = ncount;
    G->nodenames = CC_SAFE_MALLOC (ncount, node *);
    if (G->nodenames == (node **) NULL) {
        fprintf (stderr, "out of memory in buildgraph\n"); return 1;
    }

    for (i = 0; i < ncount; i++) {
        n = nodealloc (&G->node_world);
        if (n == (node *) NULL) {
            fprintf (stderr, "out of memory in buildgraph\n"); return 1;
        }
        n->name = i;
        n->adj = (edgeptr *) NULL;
        n->next = G->nodelist;
        G->nodelist = n;
        G->nodenames[i] = n;
    }
    for (i = 0; i < ecount; i++) {
        e = newedge (G, G->nodenames[elist[2 * i]],
                     G->nodenames[elist[(2 * i) + 1]]);
        if (e == (edge *) NULL) {
            fprintf (stderr, "out of memory in buildgraph\n"); return 1;
        }
        e->weight = elen[i] + elen[i];
    }
    return 0;
}

static void init_graph (graph *G)
{
    G->ncount    = 0;
    G->nodenames = (node **) NULL;
    G->nodelist  = (node *) NULL;
    G->edgelist  = (edge *) NULL;
    G->PLUS      = 1;
    G->MINUS     = 2;
    CCptrworld_init (&G->node_world);
    CCptrworld_init (&G->edge_world);
    CCptrworld_init (&G->edgeptr_world);
    CCptrworld_init (&G->shortedge_world);
}

static void free_graph (graph *G)
{
    node *n, *nnext;

    for (n = G->nodelist; n; n = nnext) {
        nnext = n->next;
        edgeptr_listfree (&G->edgeptr_world, n->adj);
        nodefree (&G->node_world, n);
    }
    G->nodelist = (node *) NULL;

    edge_listfree (&G->edge_world, G->edgelist);
    G->edgelist = (edge *) NULL;

    CC_IFFREE (G->nodenames, node *);
    G->ncount = 0;
    G->PLUS   = 0;
    G->MINUS  = 0;
}


static edge *newedge (graph *G, node *n1, node *n2)
{
    edge *e;
    edgeptr *ep;

    e = edgealloc (&G->edge_world);
    if (e == (edge *) NULL)
        return (edge *) NULL;
    e->ends[0] = n1;
    e->ends[1] = n2;
    e->next = G->edgelist;
    G->edgelist = e;

    ep = edgeptralloc (&G->edgeptr_world);
    if (ep == (edgeptr *) NULL) {
        edgefree (&G->edge_world, e);
        return (edge *) NULL;
    }
    ep->this = e;
    ep->other = n2;
    ep->next = n1->adj;
    n1->adj = ep;

    ep = edgeptralloc (&G->edgeptr_world);
    if (ep == (edgeptr *) NULL) {
        edgefree (&G->edge_world, e);
        return (edge *) NULL;
    }
    ep->this = e;
    ep->other = n1;
    ep->next = n2->adj;
    n2->adj = ep;

    return e;
}

int CCfmatch_fractional_2match (int ncount, int ecount, int *elist, int *elen,
        CCdatagroup *dat, double *val, int *thematching, int *thedual,
        int *thebasis, int wantbasic, int silent, CCrandstate *rstate)
{
    double v, vbasic, szeit, tzeit;
    int added;
    int rval = 0;
    graph G;

    tzeit = CCutil_zeit ();

    rval = buildgraph (&G, ncount, ecount, elist, elen);
    if (rval) {
        fprintf (stderr, "buildgraph failed\n"); goto CLEANUP;
    }

    initmat (&G);
    rval = twomatch (&G);
    if (rval) {
        fprintf (stderr, "twomatch failed\n"); goto CLEANUP;
    }
    rval = chkmat (&G, &v);
    if (rval) {
        fprintf (stderr, "Chkmat found error in matching\n"); goto CLEANUP;
    }

    if (!silent) {
        printf ("Fractional Matching: %.1f\n", v);
        printf ("Initial Running Time: %.2f (seconds)\n",
                CCutil_zeit () - tzeit);
        fflush (stdout);
    }

    if (wantbasic) {
        szeit = CCutil_zeit ();
        if (basicrun (&G)) {
            fprintf (stderr, "Did not find a basic optimal solution\n");
            rval = 1;
            goto CLEANUP;
        }
        if (chkmat (&G, &vbasic)) {
            fprintf (stderr, "Chkmat found error in matching\n");
            rval = 1;
            goto CLEANUP;
        }
        if (vbasic != v) {
            fprintf (stderr, "ERROR: Basis routine altered objective\n");
            rval = 1;
            goto CLEANUP;
        }
        if (!silent) {
            printf ("Basis Running Time: %.2f (seconds)\n",
                    CCutil_zeit () - szeit);
            fflush (stdout);
        }
    }

    if (dat != (CCdatagroup *) NULL) {
        if (!silent) {
            printf ("Price-Repair ...\n"); fflush (stdout);
        }
        szeit = CCutil_zeit ();
        rval = fixmatch (&G, &added, dat, rstate);
        if (rval) {
            fprintf (stderr, "fixmatch failed\n"); goto CLEANUP;
        }
        if (chkmat (&G, &v)) {
            fprintf (stderr, "Chkmat found error in matching\n");
            rval = 1;
            goto CLEANUP;
        }
        if (wantbasic) {
            do  {
                if (!silent) {
                    printf ("Find basis ...\n"); fflush (stdout);
                }
                if (basicrun (&G)) {
                    fprintf (stderr, "Did not find a basic solution\n");
                    rval = 1;
                    goto CLEANUP;
                }
                if (chkmat (&G, &vbasic)) {
                    fprintf (stderr, "Chkmat found error in matching\n");
                    rval = 1;
                    goto CLEANUP;
                }
                if (vbasic != v) {
                    fprintf (stderr, "ERROR: Basis routine altered obj\n");
                    rval = 1;
                    goto CLEANUP;
                }

                if (!silent) {
                    printf ("Price-repair basic solution ...\n");
                    fflush (stdout);
                }
                rval = fixmatch (&G, &added, dat, rstate);
                if (rval) {
                    fprintf (stderr, "fixmatch failed\n"); goto CLEANUP;
                }
                if (chkmat (&G, &v)) {
                    fprintf (stderr, "Chkmat found error in matching\n");
                    rval = 1;
                    goto CLEANUP;
                }
            } while (added);
        }
        if (!silent) {
            printf ("Running Time for Price-Repair: %.2f\n",
                    CCutil_zeit () - szeit);
            printf ("Fractional Matching on Complete Graph: %.1f\n",v);
            fflush (stdout);
        }
    }

    *val = v;

    if (thematching != (int *) NULL) {
        int k = 0;
        edge *e;
        for (e = G.edgelist; e; e = e->next) {
            if ((unsigned int) e->x != (unsigned int) ZERO) {
                thematching[k++] = e->ends[0]->name;
                thematching[k++] = e->ends[1]->name;
                thematching[k++] = ((unsigned int) e->x == (unsigned int) ONE ? 1 : 2);
            }
        }
        thematching[k] = -1;
    }
    if (thedual != (int *) NULL) {
        int i;
        for (i = 0; i < G.ncount; i++)
            thedual[i] = G.nodenames[i]->y;
    }
    if (wantbasic && thebasis != (int *) NULL) {
        int k = 0;
        edge *e;
        for (e = G.edgelist; e; e = e->next) {
            if ((unsigned int) e->basic) {
                thebasis[k++] = e->ends[0]->name;
                thebasis[k++] = e->ends[1]->name;
            }
        }
    }

CLEANUP:

    free_graph (&G);
    free_linked_world (&G);

    if (!silent) {
        printf ("Total fractional matching time: %.2f (seconds)\n",
                 CCutil_zeit () - tzeit);
        fflush (stdout);
    }

    return rval;
}

static void initmat (graph *G)
{
    node *n;
    edge *e;

    for (n = G->nodelist; n; n = n->next) {
        n->y = 0;
        n->matchcnt = 2 - MATCHDEGREE;
        n->parentedge = (edge *) NULL;
        n->label = 0;
    }
    for (e = G->edgelist; e; e = e->next) {
        e->z = 0;
        e->x = ZERO;
        e->pnext = (edge *) NULL;
    }
    G->PLUS = 1;
    G->MINUS = 2;
}

static int chkmat (graph *G, double *val)   /* val may be 1/2 integer and big */
{
    int k;
    double v = 0.0, dualval = 0.0;
    node *n;
    edge *e;
    edgeptr *ep;

    for (n = G->nodelist; n; n = n->next) {
        k = 0;
        for (ep = n->adj; ep; ep = ep->next)
            k += (unsigned int) ep->this->x;
        if (k != 2 * MATCHDEGREE) {
            fprintf (stderr, "Not a matching, node %d has 2-degree %d\n",
                                                               n->name, k);
            return 1;
        }
        dualval += (double) n->y;
    }
    dualval *= (double) MATCHDEGREE;
    for (e = G->edgelist; e; e = e->next) {
        switch ((unsigned int) e->x) {
        case (unsigned int) TWO:
            if (e->z < 0 ||
                    e->z != e->ends[0]->y + e->ends[1]->y - e->weight) {
                fprintf (stderr, "Error in dual solution - 2\n");
                return 1;
            }
            v += (double) e->weight;
            v += (double) e->weight;
            dualval -= (double) e->z;
            break;
        case (unsigned int) ONE:
            if (e->z != 0 ||
                    e->ends[0]->y + e->ends[1]->y != e->weight) {
                fprintf (stderr, "Error in dual solution - 1\n");
                return 1;
            }
            v += (double) e->weight;
            break;
        case (unsigned int) ZERO:
            if (e->z != 0 || e->ends[0]->y + e->ends[1]->y > e->weight) {
                fprintf (stderr, "Error in dual solution - 0\n");
                return 1;
            }
            break;
        default:
            fprintf (stderr, "Error in matching values\n");
            return 1;
        }
    }
    v /= 4.0;
    dualval /= 2.0;

    if (v != dualval) {
        fprintf (stderr, "The primal and dual objective values differ.\n");
        return 1;
    }
    *val = v;
    return 0;
}


/**********  Core fractional matching routines **********/


static int twomatch (graph *G)
{
    node *n;
    int rval = 0;

    for (n = G->nodelist; n; n = n->next) {
        while ((unsigned int) n->matchcnt != (unsigned int) TWO) {
            rval = augment (G, n);
            if (rval) {
                fprintf (stderr, "augment failed - probably no fmatching\n");
                return rval;
            }
        }
    }
    return 0;
}

static int augment (graph *G, node *n)
{
    node *auglist;
    int found_augmenting_path = 0;

    G->PLUS += 2;
    G->MINUS += 2;
    n->label = G->PLUS;
    n->parentedge = (edge *) NULL;
    n->matchcnt = (unsigned int) n->matchcnt + 1;
    expand (n, &found_augmenting_path, G->PLUS, G->MINUS);
    if (found_augmenting_path)
        return 0;
    while ((auglist = dualchange (n, G->PLUS, G->MINUS)) != (node *) NULL) {
        while (auglist) {
            expand (auglist, &found_augmenting_path, G->PLUS, G->MINUS);
            if (found_augmenting_path)
                return 0;
            auglist = auglist->pnext;
        }
    }
    fprintf (stderr, "Error - dual change did not create new edges\n");
    return 1;
}

static void expand (node *n, int *found_one, int PLUS, int MINUS)
{
    node *n1;
    edgeptr *ep;
    edge *e;

    if (n->label == PLUS) {
        for (ep = n->adj; ep; ep = ep->next) {
            if ((unsigned int) ep->this->x == (unsigned int) ONE) {
                augmentplushalf (n, ep->this);
                *found_one = 1;
                return;
            }
        }
        for (ep = n->adj; ep; ep = ep->next) {
            if ((unsigned int) ep->this->x == (unsigned int) ZERO &&
                        ep->other->y + n->y == ep->this->weight) {
                e = ep->this;
                n1 = ep->other;
                if (n1->label < PLUS) {         /* n1 has no label */
                    n1->label = MINUS;
                    n1->parentedge = e;
                    if ((unsigned int) n1->matchcnt != (unsigned int) TWO) {
                        augmentpath (n1);
                        *found_one = 1;
                        return;
                    }
                    expand (n1, found_one, PLUS, MINUS);
                    if (*found_one)
                        return;
                } else if (n1->label == PLUS) {
                    augmentpair (n, n1, e);
                    *found_one = 1;
                    return;
                }
            }
        }
    } else {                    /* MINUS */
        for (ep = n->adj; ep; ep = ep->next) {
            if ((unsigned int) ep->this->x == (unsigned int) ONE) {
                augmentminushalf (n, ep->this);
                *found_one = 1;
                return;
            }
        }
        for (ep = n->adj; ep; ep = ep->next) {
            if ((unsigned int) ep->this->x == (unsigned int) TWO &&
                ep->this->z == 0) {
                e = ep->this;
                n1 = ep->other;
                if (n1->label < PLUS) {         /* n1 has no label */
                    n1->label = PLUS;
                    n1->parentedge = e;
                    expand (n1, found_one, PLUS, MINUS);
                    if (*found_one)
                        return;
                } else if (n1->label == MINUS) {
                    augmentpair (n, n1, e);
                    *found_one = 1;
                    return;
                }
            }
        }
    }
    return;
}

static node *dualchange (node *n, int PLUS, int MINUS)
{
    node *new = (node *) NULL;
    int alpha = MAXWEIGHT;

    minalpha (n, &new, &alpha, PLUS, MINUS);
    if (alpha == MAXWEIGHT) {
        fprintf (stderr, "Dual change required, but no candidate edges\n");
        return (node *) NULL;
    }
    if (alpha & 0x1) {
        fprintf (stderr, "Whoops, 2 * alpha = %d, not even\n", alpha);
        return (node *) NULL;
    }
    alpha /= 2;
    subalpha (n, alpha, PLUS, MINUS);
    return new;
}

static void minalpha (node *n, node **new, int *alpha, int PLUS, int MINUS)
{
    int minv = MAXWEIGHT;
    int thisv;
    node *n1;
    edgeptr *ep;
    edge *e;

    if (n->label == PLUS) {
        for (ep = n->adj; ep; ep = ep->next) {
            e = ep->this;
            if ((unsigned int) e->x == (unsigned int) ZERO) {
                n1 = ep->other;
                if (n1->label < PLUS) {         /* n1 is unlabeled */
                    thisv = e->weight - n->y - n1->y;
                    thisv += thisv;
                    if (thisv < minv)
                        minv = thisv;
                } else if (n1->label == PLUS) {
                    thisv = e->weight - n->y - n1->y;
                    if (thisv < minv)
                        minv = thisv;
                } else {        /* n1 has a minus label */
                    if (n1->parentedge == e)
                        minalpha (n1, new, alpha, PLUS, MINUS);
                }
            }
        }
    } else {                    /* MINUS case */
        for (ep = n->adj; ep; ep = ep->next) {
            e = ep->this;
            if ((unsigned int) e->x == (unsigned int) TWO) {
                n1 = ep->other;
                if (n1->label < PLUS) {
                    thisv = e->z + e->z;
                    if (thisv < minv)
                        minv = thisv;
                } else if (n1->label == PLUS) {
                    if (n1->parentedge == e)
                        minalpha (n1, new, alpha, PLUS, MINUS);
                } else {        /* n1 has a MINUS label */
                    thisv = e->z;
                    if (thisv < minv)
                        minv = thisv;
                }
            }
        }
    }
    if (minv < *alpha) {
        *new = n;
        n->pnext = (node *) NULL;
        *alpha = minv;
    } else if (minv == *alpha) {
        n->pnext = *new;
        *new = n;
    }
}

static void subalpha (node *n, int alpha, int PLUS, int MINUS)
{
    edgeptr *ep;
    edge *e;
    node *n1;

    if (n->label == PLUS) {
        n->y += alpha;
        for (ep = n->adj; ep; ep = ep->next) {
            e = ep->this;
            if ((unsigned int) e->x == (unsigned int) TWO)
                e->z += alpha;
            else if ((unsigned int) e->x == (unsigned int) ZERO) {
                n1 = ep->other;
                if (n1->parentedge == e && n1->label >= PLUS)
                    subalpha (n1, alpha, PLUS, MINUS);
            }
        }
    } else {                    /* MINUS */
        n->y -= alpha;
        for (ep = n->adj; ep; ep = ep->next) {
            e = ep->this;
            if ((unsigned int) e->x == (unsigned int) TWO) {
                e->z -= alpha;
                n1 = ep->other;
                if (n1->parentedge == e && n1->label >= PLUS)
                    subalpha (n1, alpha, PLUS, MINUS);
            }
        }
    }
}

static void augmentpath (node *n)
{
    n->matchcnt = (unsigned int) n->matchcnt + 1;
    augmentpath1 (n);
}

static void augmentpath1 (node *n)
{
    while (n->parentedge != (edge *) NULL) {
        n->parentedge->x = (unsigned int) TWO - (unsigned int) n->parentedge->x;
        n = OTHEREND(n->parentedge, n);
    }
}

static void augmentplushalf (node *n, edge *e)
{
    flipcycle (n, e, TWO);
    augmentpath1 (n);
}

static void augmentminushalf (node *n, edge *e)
{
    flipcycle (n, e, ZERO);
    augmentpath1 (n);
}

static void flipcycle (node *n, edge *e, unsigned char v)
{
    edge *e1;
    edge *e2;

    e1 = e->pnext;
    if (e1->ends[0] == n || e1->ends[1] == n)
        e = e1;
    e1 = e;
    do {
        e1->x = v;
        v = (unsigned int) TWO - (unsigned int) v;
        e2 = e1->pnext;
        e1->pnext = (edge *) NULL;
        e1 = e2;
    } while (e1 != e);
}

static void augmentpair (node *n, node *n1, edge *e)
{
    node *n2, *n3;

    setflag (n, FALSE);
    setflag (n1, TRUE);
    n2 = findflag (n);
    if ((n3 = findhole (n, n2)) != (node *) NULL) {
        n3->matchcnt = (unsigned int) n3->matchcnt + 1;
        flipseg (n, n3);
        e->x = (unsigned int) TWO - (unsigned int) e->x;
        augmentpath1 (n1);
        return;
    }
    if ((n3 = findhole (n1, n2)) != (node *) NULL) {
        n3->matchcnt = (unsigned int) n3->matchcnt + 1;
        flipseg (n1, n3);
        e->x = (unsigned int) TWO - (unsigned int) e->x;
        augmentpath1 (n);
        return;
    }
    stringup (n, n1, n2, e);
    augmentpath1 (n2);
}

static void setflag (node *n, unsigned char v)
{
    n->flag = v;
    while (n->parentedge != (edge *) NULL) {
        n = OTHEREND(n->parentedge, n);
        n->flag = v;
    }
}

static node *findflag (node *n)
{
    while ((unsigned int) n->flag == 0) {
        n = OTHEREND(n->parentedge, n);
    }
    return n;
}

static node *findhole (node *n, node *n2)
{
    while (n != n2) {
        if ((unsigned int) n->matchcnt != (unsigned int) TWO)
            return n;
        n = OTHEREND(n->parentedge, n);
    }
    return (unsigned int) n2->matchcnt != (unsigned int) TWO ? n2
        : (node *) NULL;
}

static void flipseg (node *n, node *n2)
{
    while (n != n2) {
        n->parentedge->x = (unsigned int) TWO - (unsigned int) n->parentedge->x;
        n = OTHEREND(n->parentedge, n);
    }
}

static void stringup (node *n, node *n1, node *n2, edge *e)
{
    edge *preve, *savee;

    preve = e;
    while (n != n2) {
        n->parentedge->x = ONE;
        preve->pnext = n->parentedge;
        preve = n->parentedge;
        n = OTHEREND(n->parentedge, n);
    }
    savee = preve;
    preve = e;
    while (n1 != n2) {
        n1->parentedge->x = ONE;
        n1->parentedge->pnext = preve;
        preve = n1->parentedge;
        n1 = OTHEREND(n1->parentedge, n1);
    }
    savee->pnext = preve;
    e->x = ONE;
}

static edge *findedge (node *n1, node *n2)
{
    edgeptr *ep;

    for (ep = n1->adj; ep; ep = ep->next) {
        if (ep->other == n2)
            return ep->this;
    }
    return  (edge *) NULL;
}


/********** Basis finding routines **********/


static int basicrun (graph *G)
{
    node *n;
    edge *e;

    G->PLUS = 1;
    G->MINUS = 2;

    for (n = G->nodelist; n; n = n->next) {
        n->label = 0;
        n->flag = 0;
    }
    for (e = G->edgelist; e; e = e->next)
        e->basic = 0;

    for (n = G->nodelist; n; n = n->next) {
        if (n->label == 0) {
            if (basicscan (G, n))
                return 1;
        }
    }
    for (n = G->nodelist; n; n = n->next) {
        if ((unsigned int) n->flag == 0) {
            if (basicgrow (G, n))
                return 1;
        }
    }

    for (n = G->nodelist; n; n = n->next)
        n->label = 0;
    for (n = G->nodelist; n; n = n->next) {
        if (n->label == 0) {
            if (basic_check_scan (G, n))
                return 1;
        }
    }

    return 0;
}

static int basicscan (graph *G, node *n)
{
    edge *odd_circuit = (edge *) NULL;

    G->PLUS += 2;
    G->MINUS += 2;
    n->parentedge = (edge *) NULL;
    if (basic_grab_ones (n, 0, &odd_circuit, G->PLUS, G->MINUS))
        return 1;
    if (odd_circuit != (edge *) NULL) {
        basic_mark_component_as_done (n);
    }
    return 0;
}

static int basic_check_scan (graph *G, node *n)
{
    edge *odd_circuit = (edge *) NULL;

    G->PLUS += 2;
    G->MINUS += 2;
    n->parentedge = (edge *) NULL;
    if (basic_checkout_basic (n, 0, &odd_circuit, G->PLUS, G->MINUS))
        return 1;
    if (odd_circuit == (edge *) NULL) {
        printf ("No odd circuit\n");
        return 1;
    }
    return 0;
}

static int basicgrow (graph *G, node *n)
{
    int hit_odd_circuit = 0;
    node *expandlist;

    G->PLUS += 2;
    G->MINUS += 2;
    basic_grab_basic (n, 0, G->PLUS, G->MINUS);

    n->parentedge = (edge *) NULL;
    basic_expand (n, &hit_odd_circuit, G->PLUS, G->MINUS);
    if (hit_odd_circuit)
        return 0;
    else {
        while ((expandlist = basic_dualchange (n, G->PLUS, G->MINUS))
                      != (node *) NULL) {
            while (expandlist) {
                basic_expand (expandlist, &hit_odd_circuit, G->PLUS, G->MINUS);
                if (hit_odd_circuit)
                    return 0;
                expandlist = expandlist->pnext;
            }
        }
        fprintf (stderr, "ERROR: No dual change in basis finding code\n");
        return 1;
    }
}

static int basic_grab_ones (node *n, int parity, edge **odd_circuit,
        int PLUS, int MINUS)
{
    edge *e;
    edgeptr *ep;
    node *n1;

    n->label = (parity == 0 ? PLUS : MINUS);
    for (ep = n->adj; ep; ep = ep->next) {
        e = ep->this;
        if ((unsigned int) e->x == (unsigned int) ONE && e != n->parentedge) {
            n1 = ep->other;
            if (n1->label == 0) {
                n1->parentedge = e;
                e->basic = 1;
                if (basic_grab_ones (n1, 1 - parity, odd_circuit, PLUS, MINUS))
                    return 1;
            } else if (n1->label == n->label) {
                if (*odd_circuit == (edge *) NULL) {
                    *odd_circuit = e;
                    e->basic = 1;
                } else if (*odd_circuit != e) {
                    fprintf (stderr, "ERROR: Two odd circuits in 1-graph\n");
                    printf ("Circuit forming edges: %d-%d  %d-%d\n",
                         (*odd_circuit)->ends[0]->name,
                         (*odd_circuit)->ends[1]->name,
                         e->ends[0]->name,
                         e->ends[1]->name);
                    return 1;
                }
            } else {
                fprintf (stderr, "ERROR: Even circuit in 1-graph\n");
                printf ("Circuit forming edge: %d-%d\n",
                              e->ends[0]->name,
                              e->ends[1]->name);
                return 1;
            }
        }
    }
    return 0;
}

static int basic_checkout_basic (node *n, int parity, edge **odd_circuit,
        int PLUS, int MINUS)
{
    edge *e;
    edgeptr *ep;
    node *n1;

    n->label = (parity == 0 ? PLUS : MINUS);
    for (ep = n->adj; ep; ep = ep->next) {
        e = ep->this;
        if ((unsigned int) e->basic && e != n->parentedge) {
            n1 = ep->other;
            if (n1->label == 0) {
                n1->parentedge = e;
                if (basic_checkout_basic (n1, 1 - parity, odd_circuit,
                                          PLUS, MINUS))
                    return 1;
            } else if (n1->label == n->label) {
                if (*odd_circuit == (edge *) NULL) {
                    *odd_circuit = e;
                } else if (*odd_circuit != e) {
                    fprintf (stderr, "ERROR: Two odd circuits in basish\n");
                    printf ("Circuit forming edges: %d-%d  %d-%d\n",
                         (*odd_circuit)->ends[0]->name,
                         (*odd_circuit)->ends[1]->name,
                         e->ends[0]->name,
                         e->ends[1]->name);
                    return 1;
                }
            } else {
                fprintf (stderr, "ERROR: Even circuit in basis\n");
                printf ("Circuit forming edge: %d-%d\n",
                              e->ends[0]->name,
                              e->ends[1]->name);
                return 1;
            }
        }
    }
    return 0;
}

static void basic_grab_basic (node *n, int parity, int PLUS, int MINUS)
{
    edgeptr *ep;

    n->label = (parity == 0 ? PLUS : MINUS);
    n->flag = 1;
    for (ep = n->adj; ep; ep = ep->next) {
        if ((unsigned int) ep->this->basic) {
            if ((unsigned int) ep->other->flag == 0)
                basic_grab_basic (ep->other, 1 - parity, PLUS, MINUS);
        }
    }
}

static void basic_mark_component_as_done (node *n)
{
    edgeptr *ep;

    n->flag = 1;
    for (ep = n->adj; ep; ep = ep->next) {
        if ((unsigned int) ep->this->basic) {
            if ((unsigned int) ep->other->flag == 0)
                basic_mark_component_as_done (ep->other);
        }
    }
}

static void basic_expand (node *n, int *hit_odd_circuit, int PLUS, int MINUS)
{
    edge *e;
    edgeptr *ep;
    node *n1;

    for (ep = n->adj; ep; ep = ep->next) {
        e = ep->this;
        if (e != n->parentedge) {
            n1 = ep->other;
            if ((unsigned int) e->basic) {
                n1->parentedge = e;
                basic_expand (n1, hit_odd_circuit, PLUS, MINUS);
                if (*hit_odd_circuit)
                    return;
            } else if (n->y + n1->y == e->weight) {
                if (n1->label < PLUS) {
                    e->basic = 1;
                    if ((unsigned int) n1->flag) {
                        *hit_odd_circuit = 1;
                        return;
                    } else {
                        n1->parentedge = e;
                        if (n->label == PLUS)
                            basic_grab_basic (n1, 1, PLUS, MINUS);
                        else
                            basic_grab_basic (n1, 0, PLUS, MINUS);
                        basic_expand (n1, hit_odd_circuit, PLUS, MINUS);
                        if (*hit_odd_circuit)
                            return;
                    }
                } else if (n1->label == n->label) {
                    e->basic = 1;
                    *hit_odd_circuit = 1;
                    return;
                }
            }
        }
    }
}

static node *basic_dualchange (node *n, int PLUS, int MINUS)
{
    node *new = (node *) NULL;
    int alpha = MAXWEIGHT;

    basic_minalpha (n, &new, &alpha, 0, PLUS, MINUS);
    if (alpha != MAXWEIGHT) {
        alpha /= 2;
        basic_subalpha (n, alpha, 0, PLUS, MINUS);
    } else {
        /* reverse sense of PLUS and MINUS */
        basic_minalpha (n, &new, &alpha, 1, PLUS, MINUS);
        if (alpha == MAXWEIGHT) {
            printf ("Basic dual change required, but no candidate edges\n");
            return (node *) NULL;
        }
        alpha /= 2;
        basic_subalpha (n, alpha, 1, PLUS, MINUS);
    }
    return new;
}

static void basic_minalpha (node *n, node **new, int *alpha,
        int flip_plus_and_minus, int PLUS, int MINUS)
{
    int minv = MAXWEIGHT;
    int thisv;
    node *n1;
    edge *e;
    edgeptr *ep;

    if ((n->label == PLUS && !flip_plus_and_minus) ||
        (n->label == MINUS && flip_plus_and_minus)) {
        for (ep = n->adj; ep; ep = ep->next) {
            e = ep->this;
            n1 = ep->other;
            if ((unsigned int) e->x != (unsigned int) TWO) {
                if (n1->label < PLUS) {         /* n1 is unlabeled */
                    thisv = e->weight - n->y - n1->y;
                    thisv += thisv;
                    if (thisv < minv)
                        minv = thisv;
                } else if ((n1->label == PLUS && !flip_plus_and_minus) ||
                           (n1->label == MINUS && flip_plus_and_minus)) {
                    thisv = e->weight - n->y - n1->y;
                    if (thisv < minv)
                        minv = thisv;
                } else {        /* n1 has a minus label */
                    if (n1->parentedge == e)
                        basic_minalpha (n1, new, alpha, flip_plus_and_minus,
                                        PLUS, MINUS);
                }
            } else if ((unsigned int) e->basic && n1->parentedge == e) {
                basic_minalpha (n1, new, alpha, flip_plus_and_minus,
                                PLUS, MINUS);
            }
        }
    } else {                    /* MINUS case */
        for (ep = n->adj; ep; ep = ep->next) {
            e = ep->this;
            n1 = ep->other;
            if ((unsigned int) e->x == (unsigned int) TWO) {
                if (n1->label < PLUS) {
                    thisv = e->z + e->z;
                    if (thisv < minv)
                        minv = thisv;
                } else if ((n1->label == PLUS && !flip_plus_and_minus) ||
                           (n1->label == MINUS && flip_plus_and_minus)) {
                    if (n1->parentedge == e)
                        basic_minalpha (n1, new, alpha, flip_plus_and_minus,
                                        PLUS, MINUS);
                } else {        /* n1 has a MINUS label */
                    thisv = e->z;
                    if (thisv < minv)
                        minv = thisv;
                }
            } else if ((unsigned int) e->basic && n1->parentedge == e) {
                basic_minalpha (n1, new, alpha, flip_plus_and_minus,
                                PLUS, MINUS);
            }
        }
    }

    if (minv < *alpha) {
        *new = n;
        n->pnext = (node *) NULL;
        *alpha = minv;
    } else if (minv == *alpha) {
        n->pnext = *new;
        *new = n;
    }
}

static void basic_subalpha (node *n, int alpha, int flip_plus_and_minus,
        int PLUS, int MINUS)
{
    edge *e;
    edgeptr *ep;

    if ((n->label == PLUS && !flip_plus_and_minus) ||
        (n->label == MINUS && flip_plus_and_minus)) {
        n->y += alpha;
        for (ep = n->adj; ep; ep = ep->next) {
            e = ep->this;
            if ((unsigned int) e->x == (unsigned int) TWO)
                e->z += alpha;
            if ((unsigned int) e->basic) {
                if (ep->other->parentedge == e)
                    basic_subalpha (ep->other, alpha, flip_plus_and_minus,
                                    PLUS, MINUS);
            }
        }
    } else {                    /* MINUS */
        n->y -= alpha;
        for (ep = n->adj; ep; ep = ep->next) {
            e = ep->this;
            if ((unsigned int) e->x == (unsigned int) TWO)
                e->z -= alpha;
            if ((unsigned int) e->basic) {
                if (ep->other->parentedge == e)
                    basic_subalpha (ep->other, alpha, flip_plus_and_minus,
                                    PLUS, MINUS);
            }
        }
    }
}


/********** Price - Repair routines **********/


static int checkoutedge (graph *G, node *n1, node *n2, int *hit,
        CCdatagroup *dat)
{
    int w, wbar;
    edge *e;

    *hit = 0;
    w = CCutil_dat_edgelen (n1->name, n2->name, dat);
    w += w;
    wbar = w - n1->y - n2->y;
    if (wbar < 0) {
        if ((e = findedge (n1, n2)) != (edge *) NULL) {
            if (e->z != -wbar) {
                printf ("Hmmm.  edge (%d-%d) has z %d, wbar %d\n",
                e->ends[0]->name, e->ends[1]->name, e->z, wbar);
            }
        } else {
            if (addbadedge (G, n1, n2, w)) {
                fprintf (stderr, "addbadedge failed\n");
                return 1;
            }
            *hit = 1;
        }
    }
    return 0;
}

static int precheckoutedge (node *n1, node *n2, shortedge **list,
        CCdatagroup *dat, CCptrworld *shortedge_world)
{
    int w, wbar;
    edge *e;
    shortedge *s;

    w = CCutil_dat_edgelen (n1->name, n2->name, dat);
    w += w;
    wbar = w - n1->y - n2->y;
    if (wbar < 0) {
        if ((e = findedge (n1, n2)) != (edge *) NULL) {
            if (e->z != -wbar) {
                printf ("Hmmm.  edge (%d-%d) has z %d, wbar %d\n",
                e->ends[0]->name, e->ends[1]->name, e->z, wbar);
            }
        } else {
            s = shortedgealloc (shortedge_world);
            s->ends[0] = n1;
            s->ends[1] = n2;
            s->next = *list;
            *list = s;
            return 1;
        }
    }
    return 0;
}

static int fixmatch (graph *G, int *radded, CCdatagroup *dat,
        CCrandstate *rstate)
{
    int datnorm;

    CCutil_dat_getnorm (dat, &datnorm);
    if ((datnorm & CC_NORM_BITS) == CC_KD_NORM_TYPE)
        return kd_fixmatch (G, radded, dat, rstate);
    else if ((datnorm & CC_NORM_BITS) == CC_X_NORM_TYPE)
        return x_fixmatch (G, radded, dat);
    else
        return junk_fixmatch (G, radded, dat);
}

#define NEAR_TRY_NUM 1   /* The number of nearest (wbar) neighbors         */

#define PULL_DIVISOR 100 /* Do not pull more than ncount/PULL_DIVISOR.     */
#define PULL_UNIT    100 /* Pull if PULL_UNIT nodes will cut the spread.   */
#define PULL_CUT       2 /* A unit must cut at least spread/PULL_CUT.      */

static int kd_fixmatch (graph *G, int *radded, CCdatagroup *dat,
        CCrandstate *rstate)
{
    int rval = 0;
    int i, j, added, totaladded = 0;
    int hit, passcount = 0;
    int maxy = -MAXWEIGHT;
    int miny =  MAXWEIGHT;
    double *xcoord = (double *) NULL;
    double *ycoord = (double *) NULL;
    double *wcoord = (double *) NULL;
    CCdatagroup ldat;
    node *n;
/*
    NEEDED WHEN RADIX SORT IS WORKING
    node *ysorted;
*/
    node **heavy, **light, **order = (node **) NULL;
    int top, bottom, nlight, nheavy = 0;  /* nheavy should be set to 0 */
    int *invnames = (int *) NULL;
    double *lbnds = (double *) NULL, *gbnds = (double *) NULL;
    int datnorm;

    *radded = 0;

    CCutil_init_datagroup (&ldat);

    CCutil_dat_getnorm (dat, &datnorm);
    if (CCutil_dat_setnorm (&ldat, datnorm)) {
        rval = 1;
        goto CLEANUP;
    }

    for (n = G->nodelist; n; n = n->next) {
        if (n->y < miny)
            miny = n->y;
        if (n->y > maxy)
            maxy = n->y;
    }
    printf ("Node weight spread: (%d, %d)\n", miny, maxy);
    fflush (stdout);

/*
    THIS CODE CANNOT BE USED UNDER OS2 WITH CURRENT RADIX
    for (n = G->nodelist; n; n = n->next)
        n->pnext = n->next;
    ysorted = (node *) CCutil_linked_radixsort ((char *) G->nodelist,
        (char *) &(G->nodelist->pnext),
        (char *) &(G->nodelist->y), sizeof (int));

    order = CC_SAFE_MALLOC (G->ncount, node *);
    if (!order) {
        rval = 1;
        goto CLEANUP;
    }

    THIS IS THE CODE WHEN RADIXSORT WORKS WITH NEGATIVES
    for (i = 0, n = ysorted; n; i++, n = n->pnext) {
        order[i] = n;
    }

    INSTEAD, THIS CODE WORKS WITH CURRENT RADIX ON THE RS6000
    for (n = ysorted; n; n = n->pnext)
        if (n->y < 0)
            break;
    for (i = 0; n; n = n->pnext, i++)
        order[i] = n;
    for (n = ysorted; n && n->y >= 0; n = n->pnext, i++)
        order[i] = n;
*/

    {
        /* ONLY HERE UNTIL RADIX WORKS */
        int *y;

        order = CC_SAFE_MALLOC (G->ncount, node *);
        if (!order) {
            rval = 1;
            goto CLEANUP;
        }
        y = CC_SAFE_MALLOC (G->ncount, int);
        if (!y) {
            rval = 1;
            goto CLEANUP;
        }
        for (i = 0; i < G->ncount; i++) {
            order[i] = G->nodenames[i];
            y[i] = G->nodenames[i]->y;
        }
        y_quicksort (order, y, 0, G->ncount - 1);

        CC_FREE (y, int);
    }


    {
        int new, newspread, newtop, newbottom, spread;

        newtop = top = -1;
        newbottom = bottom = G->ncount;
        nheavy = 0;
        spread = maxy - miny;

        do {
            new = 0;
            while (new < PULL_UNIT && nheavy + new < G->ncount/ PULL_DIVISOR) {
                if (order[newbottom - 1]->y > -2 * order[newtop + 1]->y)
                    newbottom--;
                else
                    newtop++;
                new++;
            }
            newspread = order[newbottom - 1]->y - order[newtop + 1]->y;
            if (spread - newspread > spread / PULL_CUT) {
                spread = newspread;
                bottom = newbottom;
                top = newtop;
                nheavy += new;
            }
        } while (spread == newspread && nheavy < G->ncount/PULL_DIVISOR);
    }


    printf ("Truncated %d nodes to get spread: (%d, %d)\n",
        nheavy, order[top + 1]->y, order[bottom - 1]->y);
    fflush (stdout);


    if (nheavy) {
        heavy = order;
        light = order + nheavy;
        nlight = G->ncount - nheavy;
/*
        THIS IS THE CODE WHEN RADIXSORT WORKS WITH NEGATIVES
        for (i = 0, n = ysorted; i <= top; i++, n = n->pnext) {
            heavy[i] = n;
        }
        for (i = 0; i < nlight; i++, n = n->pnext) {
            light[i] = n;
        }
        for (i = top + 1; i < nheavy; i++, n = n->pnext) {
            heavy[i] = n;
        }
*/
        {
            node **temporder = (node **) NULL;
            int k;

            temporder = CC_SAFE_MALLOC (G->ncount, node *);
            if (!temporder) {
                rval = 1;
                goto CLEANUP;
            }
            for (i = 0; i < G->ncount; i++)
                temporder[i] = order[i];

            for (i = 0, k = 0; i <= top; i++)
                heavy[i] = temporder[k++];
            for (i = 0; i < nlight; i++)
                light[i] = temporder[k++];
            for (i = top + 1; i < nheavy; i++)
                heavy[i] = temporder[k++];
                CC_FREE (temporder, node *);
        }

        lbnds = CC_SAFE_MALLOC (nheavy, double);
        if (!lbnds) {
            rval = 1;
            goto CLEANUP;
        }
        gbnds = CC_SAFE_MALLOC (nheavy, double);
        if (!gbnds) {
            rval = 1;
            goto CLEANUP;
        }
        xcoord = CC_SAFE_MALLOC (nlight, double);
        if (!xcoord) {
            rval = 1;
            goto CLEANUP;
        }
        ldat.x = xcoord;
        ycoord = CC_SAFE_MALLOC (nlight, double);
        if (!ycoord) {
            rval = 1;
            goto CLEANUP;
        }
        ldat.y = ycoord;
        for (i = 0; i < nlight; i++) {
            xcoord[i] = dat->x[light[i]->name];
            ycoord[i] = dat->y[light[i]->name];
        }
    } else {
        nlight = G->ncount;
        light = G->nodenames;
        heavy = (node **) NULL;
        xcoord = dat->x;
        ycoord = dat->y;
        ldat.x = xcoord;
        ldat.y = ycoord;
        CC_FREE (order, node *);
    }

    wcoord = CC_SAFE_MALLOC (nlight, double);
    if (!wcoord) {
        rval = 1;
        goto CLEANUP;
    }
    invnames = CC_SAFE_MALLOC (G->ncount, int);
    if (!invnames) {
        rval = 1;
        goto CLEANUP;
    }
    for (i = 0; i < nlight; i++)
        invnames[light[i]->name] = i;
    for (i = 0; i < nheavy; i++)
        invnames[heavy[i]->name] = -i;


    do {
        int nodeschecked = 0;
        int saver = 0;
        int list[NEAR_TRY_NUM];
        shortedge *s, *snext, *slist = (shortedge *) NULL;
        CCkdtree localkt;
        added = 0;
        maxy = -MAXWEIGHT;
        for (i = 0; i < nlight; i++) {
            if (light[i]->y > maxy)
                maxy = light[i]->y;
        }
        for (i = 0; i < nlight; i++)
            wcoord[i] = ((double) (maxy - light[i]->y)) * 0.5;
        if (CCkdtree_build (&localkt, nlight, &ldat, wcoord, rstate)) {
            fprintf (stderr, "Unable to build CCkdtree\n");
            rval = 1;
            goto CLEANUP;
        }
        for (i = 0; i < nheavy; i++) {                    /* 1.0 for floats */
            lbnds[i] = dat->x[heavy[i]->name] -
                          (((double) (heavy[i]->y)) * 0.5) - 1.0;
            gbnds[i] = dat->x[heavy[i]->name] +
                          (((double) (heavy[i]->y)) * 0.5) + 1.0;
        }

        for (i = 0; i < nlight; i++) {
            if (light[i]->label != -1) {
                edgeptr *ep;
                hit = 0;
                for (ep = light[i]->adj; ep; ep = ep->next) {
                    if ((unsigned int) ep->this->x != (unsigned int) ZERO &&
                        invnames[ep->other->name] >= 0)
                        CCkdtree_delete (&localkt, invnames[ep->other->name]);
                }
                nodeschecked++;
                if (CCkdtree_node_k_nearest (&localkt, nlight, i, NEAR_TRY_NUM,
                                           &ldat, wcoord, list, rstate)) {
                    fprintf (stderr, "node nearest failed\n");
                    CCkdtree_free (&localkt);
                    rval = 1;
                    goto CLEANUP;
                }
                for (j = NEAR_TRY_NUM - 1; j >= 0; j--) {
                    if (list[j] != -1)
                        hit += precheckoutedge (light[i], light[list[j]],
                                        &slist, dat, &G->shortedge_world);
                }
                for (j = 0; j < nheavy; j++) {
                    if (heavy[j]->label == -1) {
                        if (xcoord[i] +
                               (((double) light[i]->y) * 0.5) > lbnds[j] &&
                            xcoord[i] -
                               (((double) light[i]->y) * 0.5) < gbnds[j]) {
                            hit += precheckoutedge (light[i], heavy[j],
                                            &slist, dat, &G->shortedge_world);
                        } else {
                            saver++;
                        }
                    }
                }
                added += hit;
                if (hit == 0)
                    light[i]->label = -1;
                for (ep = light[i]->adj; ep; ep = ep->next) {
                    if ((unsigned int) ep->this->x != (unsigned int) ZERO &&
                        invnames[ep->other->name] >= 0)
                        CCkdtree_undelete (&localkt,
                                           invnames[ep->other->name]);
                }
            }
        }
        for (j = 0; j < nheavy; j++) {
            if (heavy[j]->label != -1) {
                hit = 0;
                nodeschecked++;
                for (i = 0; i < nlight; i++) {
                    if (xcoord[i]+(((double) light[i]->y) * 0.5) > lbnds[j] &&
                        xcoord[i]-(((double) light[i]->y) * 0.5) < gbnds[j]) {
                        hit += precheckoutedge (light[i], heavy[j], &slist,
                                                dat, &G->shortedge_world);
                    } else {
                        saver++;
                    }
                }
                for (i = 0; i < j; i++) {
                    if (dat->x[heavy[i]->name] +
                          (((double) heavy[i]->y) * 0.5) > lbnds[j] &&
                        dat->x[heavy[i]->name] -
                          (((double) heavy[i]->y) * 0.5) < gbnds[j]) {
                        hit += precheckoutedge (heavy[i], heavy[j], &slist,
                                                dat, &G->shortedge_world);
                    }
                }
                for (i = j + 1; i < nheavy; i++) {
                    if (heavy[i]->label == -1) {
                        if (dat->x[heavy[i]->name] +
                              (((double) heavy[i]->y) * 0.5) > lbnds[j] &&
                            dat->x[heavy[i]->name] -
                              (((double) heavy[i]->y) * 0.5) < gbnds[j]) {
                            hit += precheckoutedge (heavy[i], heavy[j],
                                             &slist, dat, &G->shortedge_world);
                        }
                    }
                }
                added += hit;
                if (hit == 0)
                    heavy[j]->label = -1;
            }
        }

        printf ("Need to check %d edges (saved %d checks)\n", added, saver);
        fflush (stdout);
        CCkdtree_free (&localkt);

        added = 0;
        for (s = slist; s; s = snext) {
            snext = s->next;
            if (checkoutedge (G, s->ends[0], s->ends[1], &hit, dat)) {
                fprintf (stderr, "checkoutedge failed\n");
                rval = 1;
                goto CLEANUP;
            }
            added += hit;
            shortedgefree (&G->shortedge_world, s);
        }
        totaladded += added;
        printf ("Pass %d: %d edges added (%d total), %d nodes checked\n",
                              passcount++, added, totaladded, nodeschecked);
        fflush (stdout);
    } while (added);
    *radded = totaladded;

CLEANUP:

    CC_IFFREE (invnames, int);
    CC_IFFREE (wcoord, double);
    CC_IFFREE (order, node *);
    if (nheavy) {
        CC_IFFREE (xcoord, double);
        CC_IFFREE (ycoord, double);
        CC_IFFREE (lbnds, double);
        CC_IFFREE (gbnds, double);
    }
    return rval;
}

static void initlist (graph *G, node *head, node *tail, node *head2,
        node *tail2, CCdatagroup *dat)
{
    node *n, *p;
    int bound;
    int datnorm;
    double *xcoord = dat->x;
    double scale;

    CCutil_dat_getnorm (dat, &datnorm);
    if (datnorm == CC_GEOGRAPHIC) scale = CC_GEOGRAPHIC_SCALE;
    else if (datnorm == CC_GEOM)  scale = CC_GEOM;
    else if (datnorm == CC_ATT)   scale = CC_ATT_SCALE;
    else                          scale = 1.0;

    head->sort.order = -MAXWEIGHT;
    tail->sort.order = MAXWEIGHT;
    head->sort.next = tail;
    head->sort.prev = (node **) NULL;
    tail->sort.next = (node *) NULL;
    tail->sort.prev = &(head->sort.next);
    head2->sort.order = MAXWEIGHT;
    tail2->sort.order = -MAXWEIGHT;
    head2->sort.next = tail2;
    head2->sort.prev = (node **) NULL;
    tail2->sort.next = (node *) NULL;
    tail2->sort.prev = &(head2->sort.next);

    for (n = G->nodelist; n; n = n->next) {
        bound = (2 * ((int) (scale * xcoord[n->name]))) - n->y;
        for (p = head->sort.next; p->sort.order < bound; p = p->sort.next);
        n->sort.order = bound;
        n->sort.next = p;
        n->sort.prev = p->sort.prev;
        *(n->sort.prev) = n;
        p->sort.prev = &(n->sort.next);
    }
}

static int x_fixmatch (graph *G, int *radded, CCdatagroup *dat)
{
    node *n1, *n2;
    int i;
    int added, hit;
    int nodeschecked;
    int edgeschecked;
    int totaladded = 0;
    node high_fakehead, high_faketail, low_fakehead, low_faketail;
    int bound;
    double *xcoord = dat->x;
    double scale;
    int datnorm;

    CCutil_dat_getnorm (dat, &datnorm);
    if ((datnorm & CC_NORM_BITS) != CC_X_NORM_TYPE &&
        (datnorm & CC_NORM_BITS) != CC_KD_NORM_TYPE) {
        fprintf (stderr, "Cannot run x_fixmatch with norm %d\n", datnorm);
        return 1;
    }

    if (datnorm == CC_GEOGRAPHIC) scale = CC_GEOGRAPHIC_SCALE;
    else if (datnorm == CC_GEOM)  scale = CC_GEOM;
    else if (datnorm == CC_ATT)   scale = CC_ATT_SCALE;
    else                          scale = 1.0;

    initlist (G, &high_fakehead, &high_faketail, &low_fakehead, &low_faketail,
              dat);

    do {
        added = 0;
        nodeschecked = 0;
        edgeschecked = 0;
        for (i = 0; i < G->ncount; i++) {
            n1 = G->nodenames[i];
            *(n1->sort.prev) = n1->sort.next;
            n1->sort.next->sort.prev = n1->sort.prev;
            if (n1->label != -1) {
                nodeschecked++;
                n1->label = -1;
                bound = (2 * ((int) (scale * xcoord[n1->name]))) + n1->y + 3;
                 /* Need the +3 to handle floating point data */
                for (n2 = high_fakehead.sort.next; n2->sort.order < bound;
                     n2 = n2->sort.next) {
                    edgeschecked++;
                    if (checkoutedge (G, n1, n2, &hit, dat)) {
                        fprintf (stderr, "checkoutedge failed\n");
                        return 1;
                    }
                    added += hit;
                }
                bound = (2 * ((int) (scale * xcoord[n1->name]))) - n1->y - 3;
                for (n2 = low_fakehead.sort.next; n2->sort.order > bound;
                     n2 = n2->sort.next) {
                    edgeschecked++;
                    if (checkoutedge (G, n1, n2, &hit, dat)) {
                        fprintf (stderr, "checkoutedge failed\n");
                        return 1;
                    }
                    added += hit;
                }
            }
            bound = (2 * ((int) (scale * xcoord[n1->name]))) + n1->y;
            for (n2 = low_fakehead.sort.next; n2->sort.order > bound;
                 n2 = n2->sort.next);
            n1->sort.order = bound;
            n1->sort.next = n2;
            n1->sort.prev = n2->sort.prev;
            *(n1->sort.prev) = n1;
            n2->sort.prev = &n1->sort.next;
        }
        totaladded += added;
        printf ("Forward pass completed, %d nodes checked, %d edges checked\n",
                nodeschecked, edgeschecked);
        printf ("    %d edges added, total %d edges added\n",
                added, totaladded);
        if (added == 0)
            break;
        added = 0;
        nodeschecked = 0;
        edgeschecked = 0;
        for (i = G->ncount - 1; i >= 0; i--) {
            n1 = G->nodenames[i];
            *(n1->sort.prev) = n1->sort.next;
            n1->sort.next->sort.prev = n1->sort.prev;
            if (n1->label != -1) {
                nodeschecked++;
                n1->label = -1;
                bound = (2 * ((int) (scale * xcoord[n1->name]))) + n1->y + 3;
                for (n2 = high_fakehead.sort.next; n2->sort.order < bound;
                     n2 = n2->sort.next) {
                    edgeschecked++;
                    if (checkoutedge (G, n1, n2, &hit, dat)) {
                        fprintf (stderr, "checkoutedge failed\n");
                        return 1;
                    }
                    added += hit;
                }
                bound = (2 * ((int) (scale * xcoord[n1->name]))) - n1->y - 3;
                for (n2 = low_fakehead.sort.next; n2->sort.order > bound;
                     n2 = n2->sort.next) {
                    edgeschecked++;
                    if (checkoutedge (G, n1, n2, &hit, dat)) {
                        fprintf (stderr, "checkoutedge failed\n");
                        return 1;
                    }
                    added += hit;
                }
            }
            bound = (2 * ((int) (scale * xcoord[n1->name]))) - n1->y;
            for (n2 = high_fakehead.sort.next; n2->sort.order < bound;
                 n2 = n2->sort.next);
            n1->sort.order = bound;
            n1->sort.next = n2;
            n1->sort.prev = n2->sort.prev;
            *(n1->sort.prev) = n1;
            n2->sort.prev = &n1->sort.next;
        }
        totaladded += added;
        printf ("Backward pass completed, %d nodes checked, %d edges checked\n",
                nodeschecked, edgeschecked);
        printf ("    %d edges added, total %d edges added\n",
                added, totaladded);
    } while (added);
    *radded = totaladded;
    return 0;
}

static int junk_fixmatch (graph *G, int *radded, CCdatagroup *dat)
{
    int added, hit, totaladded = 0;
    node *n1, *n2;

    *radded = 0;

    do {
        added = 0;
        for (n1 = G->nodelist; n1; n1 = n1->next) {
            for (n2 = n1->next; n2; n2 = n2->next) {
                if (checkoutedge (G, n1, n2, &hit, dat)) {
                    fprintf (stderr, "checkoutedge failed\n");
                    return 1;
                }
                added += hit;
            }
        }
        totaladded += added;
        printf ("Pass completed: %d edges added, total %d edges added\n",
                 added, totaladded);
        fflush (stdout);
    } while (added);

    *radded = totaladded;
    return 0;
}

static int addbadedge (graph *G, node *n1, node *n2, int w)
{
    int wbar = -(w - n1->y - n2->y);
    edgeptr *ep;
    edge *newe;
    node *other1 = 0, *other2 = 0;

    for (ep = n1->adj; ep; ep = ep->next) {
        switch ((unsigned int) ep->this->x) {
        case (unsigned int) ONE:
            flipcycle (n1, ep->this, ZERO);
            n1->matchcnt = (unsigned int) n1->matchcnt - 1;
            break;
        case (unsigned int) TWO:
            if ((ep->this->z -= wbar) < 0) {
                ep->this->z = 0;
                ep->this->x = ZERO;
                n1->matchcnt = (unsigned int) n1->matchcnt - 1;
                ep->other->matchcnt = (unsigned int) ep->other->matchcnt - 1;
                other2 = other1;
                other1 = ep->other;
            }
            break;
        }
    }
    n1->y -= wbar;
    newe = newedge (G, n1, n2);
    if (!newe)
        return 1;
    newe->weight = w;
    newe->z = 0;
    newe->x = ZERO;
    newe->pnext = (edge *) NULL;
    while ((unsigned int) n1->matchcnt != (unsigned int) TWO)
        augment (G, n1);
    if (other1) {
        while ((unsigned int) other1->matchcnt != (unsigned int) TWO)
            augment (G, other1);
    }
    if (other2) {
        while ((unsigned int) other2->matchcnt != (unsigned int) TWO)
            augment (G, other2);
    }
    return 0;
}

static void y_quicksort (node **list, int *y, int l, int u)
{
    int i, j, t;
    int itemp;
    node *ntemp;

    if (l >= u)
        return;

    CC_SWAP (y[l], y[(l+u)/2], itemp);
    CC_SWAP (list[l], list[(l+u)/2], ntemp);

    i = l;
    j = u + 1;
    t = y[l];

    while (1) {
        do i++; while (i <= u && y[i] < t);
        do j--; while (y[j] > t);
        if (j < i) break;
        CC_SWAP (y[i], y[j], itemp);
        CC_SWAP (list[i], list[j], ntemp);
    }
    CC_SWAP (y[l], y[j], itemp);
    CC_SWAP (list[l], list[j], ntemp);
    y_quicksort (list, y, l, j - 1);
    y_quicksort (list, y, i, u);
}

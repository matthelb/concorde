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
/*  int CCtsp_tighten_lpcut_in (CCtsp_lpgraph *g, CCtsp_lpcut_in *c,        */
/*      double *x, CCtsp_lpcut_in *d, CCtsp_tighten_info *stats,            */
/*      double *pimprove)                                                   */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_tighten_lpcut (CCtsp_lpgraph *g, CCtsp_lpclique *cliques,     */
/*      CCtsp_lpcut *c, double *x, CCtsp_lpcut_in *d,                       */
/*      CCtsp_tighten_info *stats, double *pimprove)                        */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCtsp_init_tighten_info (CCtsp_tighten_info *stats)                */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCtsp_print_tighten_info (CCtsp_tighten_info *stats)               */
/*    MISSING                                                               */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"

#define USE_INLINES
#undef  TIGHTEN_NOTIEADD
#undef  NOTIGHTEN
#undef  DEBUG_SKELETONS

#define MAX_ALLOWABLE_CLIQUE_COUNT (100)    /* Due to memory requirements */
#define EPS (1e-6)

typedef struct tighten_node {
    struct tighten_node *next;              /* to maintain atomlist */
    struct tighten_node *prev;
    struct tighten_node **atomlist;         /* alpha(v) */
    struct qu_elem *moves;                    /* size ncliques */
    struct tighten_node *next_initialized;  /* to maintain initialized_nodes */
} tighten_node;

typedef struct qu_elem {
    double delta;
    tighten_node *v;                        /* node to move */
    int i;                                  /* clique to adjust */
    int chi;                                /* 0==outside, last del tied */
                                            /* 1==inside, 2 otherwise */
    int queue_handle;
} qu_elem;

typedef union atomfind {
    union atomfind *child[2];
    struct {
        tighten_node *nodelist;                 /* A(I) in writeup */
        int           atomnumber;
    } leaf;
    union atomfind *next;                   /* for freelist */
} atomfind;

typedef struct tighten_graph {
    CCtsp_lpgraph *g;
    struct tighten_node *nodes;
    double *x;
    int ncliques;
    int rhs;

    atomfind *atomtree;
    CCptrworld atomfind_world;
    int atomcount;
    tighten_node *initialized_nodes;       /* Vstar in writeup */
    CCpriority queue;
} tighten_graph;

#define NODE_IN_SINGLETON(v) ((v)->atomlist && (*((v)->atomlist))       \
                              && (!(*((v)->atomlist))->next))

static void
    add_to_atom (tighten_node *v),
    delete_from_atom (tighten_node *v),
    collect_atomfinder (atomfind *tree, int depth, CCtsp_skeleton *skel,
        tighten_node *nodelist),
    cleanup_atomfinder (CCptrworld *atomfind_world, atomfind *tree,
        int depth),
    cleanup_graph (tighten_graph *tg);

static int
    update_queue (CCpriority *q, qu_elem *qe),
    initialize_graph (CCtsp_lpgraph *g, double *x, tighten_graph *tg),
    initialize_lpcut_in (tighten_graph *tg, CCtsp_lpcut_in *c),
    initialize_lpcut (tighten_graph *tg, CCtsp_lpclique *cliques,
        CCtsp_lpcut *c),
    initialize_node (tighten_graph *tg, tighten_node *v, int add_atom),
    process_qu_elem (tighten_graph *tg, qu_elem *q),
    process_qu_elem_atom (tighten_graph *tg, tighten_node *v),
    process_qu_elem_delta (tighten_graph *tg, qu_elem *q),
    collect_new_cut (tighten_graph *tg, CCtsp_lpcut_in *cout),
    collect_skeleton (tighten_graph *tg, CCtsp_skeleton *skel),
    tighten_cut (tighten_graph *tg, CCtsp_tighten_info *stats,
        double *pimprove);

static tighten_node
  **find_atomlist (tighten_graph *tg, tighten_node *v, int add);

#ifndef USE_INLINES
static double
    qu_elem_key (qu_elem *q);
static int
    qu_elem_is_active (qu_elem *q);
#endif

CC_PTRWORLD_ROUTINES (atomfind, atomfind_alloc, atomfind_bulkalloc,
        atomfind_free)
CC_PTRWORLD_LEAKS_ROUTINE (atomfind, atomfind_leaks, child[1], atomfind *)

#ifdef USE_INLINES

#define qu_elem_key(q) (((q)->delta > EPS) ? (-(3.0 + (q)->delta))        \
                                         : (-((q)->chi + (q)->delta)))
#ifdef TIGHTEN_NOTIEADD
#define qu_elem_is_active(q) ((!NODE_IN_SINGLETON((q)->v))                \
                            && ((q)->delta > EPS                        \
                                || ((q)->delta >= -EPS && (q)->chi == 1)))
#else
#define qu_elem_is_active(q) ((!NODE_IN_SINGLETON((q)->v))                \
                            && ((q)->delta > EPS                        \
                                || ((q)->delta >= -EPS && (q)->chi > 0)))
#endif

#else /* USE_INLINES */

static double qu_elem_key (qu_elem *q)
{
    double key;

    if (q->delta > EPS) {
        key = 3.0 + q->delta;
    } else {
        key = q->chi + q->delta;
    }
    return -key;
}

static int qu_elem_is_active (qu_elem *q)
{
#ifdef TIGHTEN_NOTIEADD
    return ((!NODE_IN_SINGLETON(q->v))
            && ((q->delta > EPS)
                || ((q)->delta >= -EPS && (q)->chi == 1)));
#else
    return ((!NODE_IN_SINGLETON(q->v))
            && ((q->delta > EPS)
                || ((q)->delta >= -EPS && (q)->chi > 0)));
#endif
}

#endif

static void add_to_atom (tighten_node *v)
{
    if (v->atomlist) {
        v->next = *(v->atomlist);
        if (v->next) v->next->prev = v;
        v->prev = (tighten_node *) NULL;
        *(v->atomlist) = v;
    } else {
        v->next = (tighten_node *) NULL;
        v->prev = (tighten_node *) NULL;
    }
}

static void delete_from_atom (tighten_node *v)
{
    if (v->next) v->next->prev = v->prev;
    if (v->prev) v->prev->next = v->next;
    else if (v->atomlist) *(v->atomlist) = v->next;
}

static int update_queue (CCpriority *q, qu_elem *qe)
{
    double newkey;
    int handle;

    if (qu_elem_is_active (qe)) {
        newkey = qu_elem_key (qe);
        if (qe->queue_handle >= 0) {
            CCutil_priority_changekey (q, qe->queue_handle, newkey);
        } else {
            handle = CCutil_priority_insert (q, (void *) qe, newkey);
            if (handle < 0) {
                return handle;
            }
            qe->queue_handle = handle;
        }
    } else {
        if (qe->queue_handle >= 0) {
            CCutil_priority_delete (q, qe->queue_handle);
            qe->queue_handle = -1;
        }
    }
    return 0;
}

static tighten_node **find_atomlist (tighten_graph *tg, tighten_node *v,
                                     int add)
{
    atomfind *f = tg->atomtree;
    atomfind **pf = &(tg->atomtree);
    int i;
    int dir;

    for (i=0; i<tg->ncliques; i++) {
        if (!f) {
            if (add) {
                (*pf) = f = atomfind_alloc(&tg->atomfind_world);
                if (!f) return (tighten_node **) NULL;
                f->child[0] = (atomfind *) NULL;
                f->child[1] = (atomfind *) NULL;
            } else {
                return (tighten_node **) NULL;
            }
        }

        dir = v->moves[i].chi & 1;
        pf = &(f->child[dir]);
        f = f->child[dir];
    }
    if (!f) {
        if (add) {
            (*pf) = f = atomfind_alloc(&tg->atomfind_world);
            if (!f) return (tighten_node **) NULL;
            f->leaf.nodelist = (tighten_node *) NULL;
            f->leaf.atomnumber = tg->atomcount++;
        } else {
            return (tighten_node **) NULL;
        }
    }

    return &(f->leaf.nodelist);
}

static void collect_atomfinder (atomfind *tree, int depth,
        CCtsp_skeleton *skel, tighten_node *nodelist)
{
    tighten_node *p;
    tighten_node *pmin;
    
    if (depth == 0) {
        pmin = tree->leaf.nodelist;
        for (p = tree->leaf.nodelist; p; p = p->next) {
            if (p < pmin) pmin = p;
        }
        skel->atoms[tree->leaf.atomnumber] = (int) (pmin - nodelist);
    } else {
        if (tree->child[0]) {
            collect_atomfinder (tree->child[0], depth-1, skel, nodelist);
        }
        if (tree->child[1]) {
            collect_atomfinder (tree->child[1], depth-1, skel, nodelist);
        }
    }

}

static void cleanup_atomfinder (CCptrworld *atomfind_world, atomfind *tree,
        int depth)
{
    if (depth > 0) {
        if (tree->child[0]) {
            cleanup_atomfinder (atomfind_world, tree->child[0], depth-1);
        }
        if (tree->child[1]) {
            cleanup_atomfinder (atomfind_world, tree->child[1], depth-1);
        }
    }
    atomfind_free (atomfind_world, tree);
}

static int initialize_graph (CCtsp_lpgraph *g, double *x, tighten_graph *tg)
{
    int rval;

    tg->g = g;
    tg->x = x;
    CCptrworld_init (&tg->atomfind_world);

    tg->nodes = (tighten_node *) CC_SAFE_MALLOC (g->ncount, tighten_node);
    if (!tg->nodes) return -1;

    rval = CCutil_priority_init (&tg->queue, 1000);
    if (rval) {
        CC_FREE (tg->nodes, tighten_node);
        return rval;
    }

    tg->initialized_nodes = (tighten_node *) NULL;
    tg->atomtree = (atomfind *) NULL;
    tg->atomcount = 0;
    tg->ncliques = 0;
    g->nodemarker++;

    return 0;
}

static int initialize_lpcut_in (tighten_graph *tg, CCtsp_lpcut_in *c)
{
    int i,tmp,k,l,m;
    CCtsp_lpgraph *g = tg->g;
    CCtsp_lpclique *cl;
    tighten_node *v;
    tighten_node **vval;
    int rval;

    tg->ncliques = c->cliquecount;
    tg->rhs = c->rhs;
    for (i=0; i<tg->ncliques; i++) {
        cl = &(c->cliques[i]);
        CC_FOREACH_NODE_IN_CLIQUE (k, *cl, tmp) {
            if (g->nodes[k].mark < g->nodemarker) {
                rval = initialize_node (tg, &(tg->nodes[k]), 0);
                if (rval) return rval;
            }
            tg->nodes[k].moves[i].chi = 1;
            for (l=0; l<g->nodes[k].deg; l++) {
                m = g->nodes[k].adj[l].to;
                if (g->nodes[m].mark < g->nodemarker) {
                    rval = initialize_node (tg, &(tg->nodes[m]), 0);
                    if (rval) return rval;
                }
            }
        }
    }

    if (c->skel.atomcount == 0) {
        fprintf (stderr, "error, cut in tighten has empty skeleton\n");
        return 1;
    }
    
    /* force atoms in skeleton */
    for (i=0; i<c->skel.atomcount; i++) {
        m = c->skel.atoms[i];
        if (g->nodes[m].mark < g->nodemarker) {
            rval = initialize_node (tg, &tg->nodes[m], 0);
            if (rval) return rval;
        }
        vval = find_atomlist (tg, &(tg->nodes[c->skel.atoms[i]]), 1);
        if (vval == (tighten_node **) NULL) {
            return -1;
        }
    }

    for (v = tg->initialized_nodes; v; v = v->next_initialized) {
        v->atomlist = find_atomlist (tg, v, 0);
        if (v->atomlist != (tighten_node **) NULL) {
            add_to_atom (v);
        }
    }

    for (v = tg->initialized_nodes; v; v = v->next_initialized) {
        k = (int) (v - tg->nodes);
        for (l=0; l<g->nodes[k].deg; l++) {
            m = g->nodes[k].adj[l].to;
            if (g->nodes[m].mark == g->nodemarker) {
                for (i=0; i<tg->ncliques; i++) {
                    if (v->moves[i].chi != tg->nodes[m].moves[i].chi) {
                        v->moves[i].delta += tg->x[g->nodes[k].adj[l].edge];
                    }
                }
            }
        }
        for (i=0; i<tg->ncliques; i++) {
            rval = update_queue (&tg->queue, &v->moves[i]);
            if (rval) return rval;
        }
    }

    return 0;
}

static int initialize_lpcut (tighten_graph *tg, CCtsp_lpclique *cliques,
        CCtsp_lpcut *c)
{
    int i,tmp,k,l,m;
    CCtsp_lpgraph *g = tg->g;
    CCtsp_lpclique *cl;
    tighten_node *v;
    tighten_node **vval;
    int rval;

    tg->ncliques = c->cliquecount;
    tg->rhs = c->rhs;
    for (i=0; i<tg->ncliques; i++) {
        cl = &cliques[c->cliques[i]];
        CC_FOREACH_NODE_IN_CLIQUE (k, *cl, tmp) {
            if (g->nodes[k].mark < g->nodemarker) {
                rval = initialize_node (tg, &(tg->nodes[k]), 0);
                if (rval) return rval;
            }
            tg->nodes[k].moves[i].chi = 1;
            for (l=0; l<g->nodes[k].deg; l++) {
                m = g->nodes[k].adj[l].to;
                if (g->nodes[m].mark < g->nodemarker) {
                    rval = initialize_node (tg, &(tg->nodes[m]), 0);
                    if (rval) return rval;
                }
            }
        }
    }

    /* force atoms in skeleton */
    for (i=0; i<c->skel.atomcount; i++) {
        vval = find_atomlist (tg, &(tg->nodes[c->skel.atoms[i]]), 1);
        if (vval == (tighten_node **) NULL) {
            return -1;
        }
    }

    for (v = tg->initialized_nodes; v; v = v->next_initialized) {
        v->atomlist = find_atomlist (tg, v, 0);
        if (v->atomlist != (tighten_node **) NULL) {
            add_to_atom (v);
        }
    }

    for (v = tg->initialized_nodes; v; v = v->next_initialized) {
        k = (int) (v - tg->nodes);
        for (l=0; l<g->nodes[k].deg; l++) {
            m = g->nodes[k].adj[l].to;
            if (g->nodes[m].mark == g->nodemarker) {
                for (i=0; i<tg->ncliques; i++) {
                    if (v->moves[i].chi != tg->nodes[m].moves[i].chi) {
                        v->moves[i].delta += tg->x[g->nodes[k].adj[l].edge];
                    }
                }
            }
        }
        for (i=0; i<tg->ncliques; i++) {
            rval = update_queue (&tg->queue, &v->moves[i]);
            if (rval) return rval;
        }
    }

    return 0;
}

static int initialize_node (tighten_graph *tg, tighten_node *v, int add_atom)
{
    int i;

    v->moves = CC_SAFE_MALLOC (tg->ncliques, qu_elem);
    if (!v->moves) return -1;

    for (i=0; i<tg->ncliques; i++) {
        v->moves[i].v = v;
        v->moves[i].i = i;
        v->moves[i].delta = -1.0;
        v->moves[i].chi = 2;
        v->moves[i].queue_handle = -1;
    }
    if (add_atom) {
        v->atomlist = find_atomlist (tg, v, 0);
        if (v->atomlist) {
            add_to_atom (v);
        }
    } else {
        v->atomlist = (tighten_node **) NULL;
    }

    v->next_initialized = tg->initialized_nodes;
    tg->initialized_nodes = v;

    tg->g->nodes[v - tg->nodes].mark = tg->g->nodemarker;

    return 0;
}

static int process_qu_elem (tighten_graph *tg, qu_elem *q)
{
    int rval;

    if (q->chi == 1) {
        if (q->delta <= EPS) {
            q->chi = 0;
        } else {
            q->chi = 2;
        }
    } else {
        q->chi = 1;
    }
    rval = process_qu_elem_atom (tg, q->v);
    if (rval) return rval;
    rval = process_qu_elem_delta (tg, q);
    if (rval) return rval;
    return 0;
}

static int process_qu_elem_atom (tighten_graph *tg, tighten_node *v)
{
    tighten_node *w;
    int j;
    int rval;

    if (v->atomlist) {
        delete_from_atom (v);

        if (NODE_IN_SINGLETON(v)) {
            w = *(v->atomlist);
            for (j=0; j<tg->ncliques; j++) {
                rval = update_queue (&tg->queue, &w->moves[j]);
                if (rval) return rval;
            }
        }
    }

    /* chi was updated before the call */
    v->atomlist = find_atomlist (tg, v, 0);

    if (v->atomlist) {
        w = *(v->atomlist);

        add_to_atom (v);

        if (w->next == (tighten_node *) NULL) {
            /* v->atomlist was the singleton w */
            for (j=0; j<tg->ncliques; j++) {
                rval = update_queue (&tg->queue, &w->moves[j]);
                if (rval) return rval;
            }
        }
    }
    return 0;
}

static int process_qu_elem_delta (tighten_graph *tg, qu_elem *q)
{
    tighten_node *v = q->v;
    int i = q->i;
    tighten_node *w;
    int j;
    CCtsp_lpgraph *g = tg->g;
    int nv = (int) (v - tg->nodes);
    int nw;
    int e;
    int rval;

    v->moves[i].delta = -v->moves[i].delta;
    rval = update_queue (&tg->queue, &v->moves[i]);
    if (rval) return rval;

    for (j=0; j<g->nodes[nv].deg; j++) {
        nw = g->nodes[nv].adj[j].to;
        e = g->nodes[nv].adj[j].edge;
        w = &tg->nodes[nw];
        if (g->nodes[nw].mark < g->nodemarker) {
            rval = initialize_node (tg, w, 1);
            if (rval) return rval;
        }

        if ((v->moves[i].chi & 1) == (w->moves[i].chi & 1)) {
            w->moves[i].delta -= tg->x[e];
        } else {
            w->moves[i].delta += tg->x[e];
        }
        rval = update_queue (&tg->queue, &w->moves[i]);
        if (rval) return rval;
    }

    return 0;
}

static void cleanup_graph (tighten_graph *tg)
{
    tighten_node *v, *vnext;

    for (v = tg->initialized_nodes; v; v = vnext) {
        vnext = v->next_initialized;
        CC_FREE (v->moves, qu_elem);
    }
    cleanup_atomfinder (&tg->atomfind_world, tg->atomtree, tg->ncliques);

#if 0
    {
        int total, onlist, leak;
        if ((leak = atomfind_leaks(&tg->atomfind_world, &total, &onlist))) {
            fprintf (stderr, "TIGHTEN leaked %d atomfind's (total %d onlist %d)\n",
                     leak, total, onlist);
        }
    }
#endif

    CCptrworld_delete (&tg->atomfind_world);

    CCutil_priority_free (&tg->queue);

    CC_FREE (tg->nodes, tighten_node);
}

static int collect_lpclique (tighten_graph *tg, int cnum, CCtsp_lpclique *c)
{
    tighten_node *v;
    int cnt;
    int *arr;
    int rval;

    cnt = 0;
    for (v = tg->initialized_nodes; v; v = v->next_initialized) {
        if (v->moves[cnum].chi == 1) cnt++;
    }
    arr = CC_SAFE_MALLOC (cnt, int);
    if (!arr) return -1;

    cnt = 0;
    for (v = tg->initialized_nodes; v; v = v->next_initialized) {
        if (v->moves[cnum].chi == 1) {
            arr[cnt++] = (int) (v - tg->nodes);
        }
    }

    rval = CCtsp_array_to_lpclique (arr, cnt, c);
    CC_FREE (arr, int);
    return rval;
}

static int collect_new_cut (tighten_graph *tg, CCtsp_lpcut_in *cout)
{
    int i, j;
    int rval;

    CCtsp_init_lpcut_in (cout);
    
    cout->cliquecount = tg->ncliques;
    cout->rhs = tg->rhs;
    cout->cliques = CC_SAFE_MALLOC (tg->ncliques, CCtsp_lpclique);
    if (!cout->cliques) {
        return -1;
    }

    for (i=0; i<tg->ncliques; i++) {
        rval = collect_lpclique (tg, i, &(cout->cliques[i]));
        if (rval) {
            for (j=0; j<i; j++) {
                CC_FREE (cout->cliques[i].nodes, CCtsp_segment);
            }
            CC_FREE (cout->cliques, CCtsp_lpclique);
            return rval;
        }
    }

    return 0;
}

static int collect_skeleton (tighten_graph *tg, CCtsp_skeleton *skel)
{
    int i;
    int rval;

    CCtsp_init_skeleton (skel);
    
    skel->atoms = CC_SAFE_MALLOC (tg->atomcount, int);
    if (skel->atoms == (int *) NULL) {
        fprintf (stderr, "Out of memory in collect_skeleton\n");
        rval = 1; goto CLEANUP;
    }
    skel->atomcount = tg->atomcount;

    for (i=0; i<skel->atomcount; i++) {
        skel->atoms[i] = -1;
    }
    
    collect_atomfinder (tg->atomtree, tg->ncliques, skel, tg->nodes);

    for (i=0; i<skel->atomcount; i++) {
        if (skel->atoms[i] == -1) {
            fprintf (stderr, "collect_atomfinder didn't find atom %d\n", i);
            rval = 1; goto CLEANUP;
        }
    }

    CCutil_int_array_quicksort (skel->atoms, skel->atomcount);

    rval = 0;
    
 CLEANUP:
    if (rval) {
        CCtsp_free_skeleton (skel);
    }
    return rval;
}

static int tighten_cut (tighten_graph *tg, CCtsp_tighten_info *stats,
                        double *pimprove)
{
    qu_elem *q;
    int rval;
    double improve = 0.0;

    while ((q = (qu_elem *) CCutil_priority_deletemin (&tg->queue,
                                                       (double *) NULL))) {
        if (q->v->moves[q->i].chi == 1) {
            stats->ndel++;
            stats->del_delta += q->delta;
            if (q->delta <= EPS) stats->ndel_tied++;
        } else {
            stats->nadd++;
            stats->add_delta += q->delta;
            if (q->delta <= EPS) stats->nadd_tied++;
        }
        improve += q->delta;
        q->queue_handle = -1;
        rval = process_qu_elem (tg, q);
        if (rval) return rval;
    }
    if (pimprove) *pimprove = improve;
    return 0;
}

int CCtsp_tighten_lpcut_in (CCtsp_lpgraph *g, CCtsp_lpcut_in *c, double *x,
        CCtsp_lpcut_in *cout, CCtsp_tighten_info *stats, double *pimprove)
{
    tighten_graph tg;
    int rval = 0;
    double szeit = CCutil_zeit();

    if (c->branch != 0) {
        fprintf (stderr, "try to tighten a branch cut\n"); return 1;
    }
    if (c->sense != 'G') {
        fprintf (stderr, "try to tighten a <= cut\n"); return 1;
    }
    if (c->cliquecount > MAX_ALLOWABLE_CLIQUE_COUNT) {
        if (pimprove) *pimprove = 0.0;
        rval = CCtsp_copy_lpcut_in (c, cout);
        if (rval) {
            fprintf (stderr, "CCtsp_copy_lpcut_in failed\n"); return rval;
        }
        return 0;
    }

#ifdef DEBUG_SKELETONS
    {
        int i;
        
        printf ("lpcut_in in  skeleton:");
        for (i=0; i<c->skel.atomcount; i++) {
            printf (" %d", c->skel.atoms[i]);
        }
        printf ("\n");
        CCtsp_print_lpcut_in (c);
        fflush (stdout);
    }
#endif
    
    rval = initialize_graph (g, x, &tg);
    if (rval) return rval;

    rval = initialize_lpcut_in (&tg, c);
    if (rval) goto CLEANUP;

#ifndef NOTIGHTEN
    rval = tighten_cut (&tg, stats, pimprove);
    if (rval) goto CLEANUP;
#endif

    rval = collect_new_cut (&tg, cout);
    if (rval) goto CLEANUP;

    rval = collect_skeleton (&tg, &cout->skel);
    if (rval) goto CLEANUP;
    
    cout->branch = c->branch;
    cout->sense  = c->sense;

#ifdef DEBUG_SKELETONS
    {
        int i;
        
        printf ("lpcut_in out skeleton:");
        for (i=0; i<cout->skel.atomcount; i++) {
            printf (" %d", cout->skel.atoms[i]);
        }
        printf ("\n");
        CCtsp_print_lpcut_in (cout);
        fflush (stdout);
    }
#endif
    
    rval = 0;

CLEANUP:

    stats->ncall++;
    if (rval) stats->nfail++;
    stats->time += CCutil_zeit() - szeit;
    cleanup_graph (&tg);
    return rval;
}

int CCtsp_tighten_lpcut (CCtsp_lpgraph *g, CCtsp_lpclique *cliques,
        CCtsp_lpcut *c, double *x, CCtsp_lpcut_in *cout,
        CCtsp_tighten_info *stats, double *pimprove)
{
    tighten_graph tg;
    int rval = 0;
    double szeit = CCutil_zeit();

    if (c->branch != 0) {
        fprintf (stderr, "try to tighten a branch cut\n"); return 1;
    }
    if (c->sense != 'G') {
        fprintf (stderr, "try to tighten a <= cut\n"); return 1;
    }

    if (c->cliquecount > MAX_ALLOWABLE_CLIQUE_COUNT) {
        CCtsp_lpcuts dummy;

        dummy.cliques = cliques;
        if (pimprove) *pimprove = 0.0;
        rval = CCtsp_lpcut_to_lpcut_in (&dummy, c, cout);
        if (rval) {
            fprintf (stderr, "CCtsp_lpcut_to_lpcut_in failed\n"); return rval;
        }
        return 0;
    }

    rval = initialize_graph (g, x, &tg);
    if (rval) return rval;

    rval = initialize_lpcut (&tg, cliques, c);
    if (rval) goto CLEANUP;

    rval = tighten_cut (&tg, stats, pimprove);
    if (rval) goto CLEANUP;

    rval = collect_new_cut (&tg, cout);
    if (rval) goto CLEANUP;

    cout->branch = c->branch;
    cout->sense  = c->sense;

CLEANUP:

    stats->ncall++;
    if (rval) stats->nfail++;
    stats->time += CCutil_zeit() - szeit;
    cleanup_graph (&tg);
    return rval;
}

void CCtsp_init_tighten_info (CCtsp_tighten_info *stats)
{
    stats->ncall = 0;
    stats->nfail = 0;
    stats->nadd = 0;
    stats->nadd_tied = 0;
    stats->ndel = 0;
    stats->ndel_tied = 0;
    stats->add_delta = 0.0;
    stats->del_delta = 0.0;
    stats->time = 0.0;
}

void CCtsp_print_tighten_info (CCtsp_tighten_info *stats)
{
    printf ("TIGHTEN STATS: %d calls (%d failed), %.2f improvement, %.2f seconds\n",
            stats->ncall, stats->nfail, stats->add_delta + stats->del_delta,
            stats->time);
    printf ("               %d adds, %d tied, %.2f improvement\n",
            stats->nadd, stats->nadd_tied, stats->add_delta);
    printf ("               %d dels, %d tied, %.2f improvement\n",
            stats->ndel, stats->ndel_tied, stats->del_delta);
    fflush (stdout);
}

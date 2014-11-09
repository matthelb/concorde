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
/*  int CCchunk_finder (int ncount, int ecount, int *elist, double *elen,   */
/*      double eps, CCchunk_flag flags, CCchunk_find_timer *timer,          */
/*      CCchunk_chunk_callback *callback, CCrandstate *rstate),             */
/*    MISSING                                                               */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"
#include "localcut.h"
#include "macrorus.h"

#define MINCHUNKSIZE   5        /* do not try chunks with fewer vertices */

#undef  DUMPBADCHUNKS

#define EPS2 0.001
/* disregard chunks with edges of weight > 1+EPS2 */

#define EPS3 0.0001
/* disregard edges of weight <= EPS3 */

#define EPS4 0.0001
/* edges of weight >= 1-EPS4 are considered tight */

#define HEADER 1
#define NONHEADER 0

#define OUT -1

#define OTHEREND(f,v) (f->end1 == v ? f->end2 : f->end1)

typedef struct vertex {
    struct item        *adj;
    struct edge        *into;
    int                 name;
    int                 deg;
    double              val;
} vertex;

typedef struct edge {
    struct vertex      *end1;
    struct vertex      *end2;
    struct edge        *link;
    double              weight;
    int                 rank;
    char                status;
} edge;

typedef struct item {
    struct edge        *edgeptr;
    struct item        *next;
} item;

typedef struct graph {
    struct vertex      *vlist;
    struct edge        *elist;
    struct item        *supply;
    int                 vcount;
    int                 ecount;
} graph;

typedef struct graph_chunk_old {
    int                 ncount;
    int                 ecount;
    int                *end0;
    int                *end1;
    int                *upper;
    int                *lower;
    double             *weight;
    int               **members;
}   graph_chunk_old;


static int
    buildgraph (graph *G, int vcount, int ecount, int *elist, double *elen),
    equiv_chunkfinder (graph *G, CCchunk_flag flags, 
        CCutil_timer *timer, CCchunk_chunk_callback *callback),
    grab_all_chunks (graph *G, CCchunk_flag flags, CCutil_timer *timer,
        CCchunk_chunk_callback *callback),
   *member_dup (int *omem),
    checkout_chunk (graph_chunk_old *c, CCutil_timer *timer,
        CCchunk_chunk_callback *callback),
    chunk_ecount (vertex **v, int count),
    sphere_chunkfinder (graph *G, double eps, CCchunk_flag flags, 
        CCutil_timer *timer, CCchunk_chunk_callback *callback),
    dummy_chunkfinder (graph *G, double eps, CCchunk_flag flags,
        CCutil_timer *timer, CCchunk_chunk_callback *callback),
    permute_chunkfinder (graph *G, double eps, CCchunk_flag flags,
        CCutil_timer *timer, CCchunk_chunk_callback *callback,
        CCrandstate *rstate),
    weighted_chunkfinder (graph *G, double eps, CCchunk_flag flags,
        CCutil_timer *timer, CCchunk_chunk_callback *callback),
    dummy_grab_chunk (graph *G, int count, int *list,
        graph_chunk_old **c);

static void
    graph_chunk_old_free (graph_chunk_old *c),
    freegraph (graph *G),
    permute_edgelist (graph *G, CCrandstate *rstate),
    make_adjlists (graph *G),
    enlarge_chunk (vertex **chunkvertex, int *count, int maxchunksize),
    civilize_chunk (graph_chunk_old *c),
    merge_cycles (vertex *v),
    merge_classes (edge *e, edge *f),
    convert_lists (graph *G),
    CC_UNUSED dump_graph (graph *G);

static edge
    *cclass (edge *e);

static graph_chunk_old
    *graph_chunk_old_alloc (int ncount, int ecount),
    *grab_chunk (graph *G, vertex **chunkvertex, int count),
    *get_sphere (graph *G, vertex *v, double eps, CCchunk_flag flags),
    *get_weighted_sphere (graph *G, vertex *v, double eps,
        CCchunk_flag flags),
    *get_card_sphere (graph *G, vertex *v, double eps,
        CCchunk_flag flags);

static CCchunk_graph
    *strip_outside (graph_chunk_old *c);


int CCchunk_finder (int vcount, int ecount, int *elist, double *elen,
        double eps, CCchunk_flag flags, CCchunk_find_timer *timer,
        CCchunk_chunk_callback *callback, CCrandstate *rstate)
{
    graph G;

    CCutil_start_timer (&timer->all);
    
    if (buildgraph (&G, vcount, ecount, elist, elen)) {
        printf ("buildgraph failed\n");
        CCutil_stop_timer (&timer->all, 0);
        return 1;
    }

    if (flags.dummy) {
        if (dummy_chunkfinder (&G, eps, flags, &timer->locate, callback)) {
            fprintf (stderr, "dummy_chunkfinder failed\n");
            freegraph (&G);
            CCutil_stop_timer (&timer->all, 0);
            return 1;
        }
    } else if (flags.permute) {
        if (permute_chunkfinder (&G, eps, flags, &timer->locate, callback,
                                 rstate)) {
            fprintf (stderr, "permute_chunkfinder failed\n");
            freegraph (&G);
            CCutil_stop_timer (&timer->all, 0);
            return 1;
        }
    } else if (flags.weighted) {
        if (weighted_chunkfinder (&G, eps, flags, &timer->locate, callback)) {
            fprintf (stderr, "weighted_chunkfinder failed\n");
            freegraph (&G);
            CCutil_stop_timer (&timer->all, 0);
            return 1;
        }
    } else if (flags.spheres) {
        if (sphere_chunkfinder (&G, eps, flags, &timer->locate, callback)) {
            fprintf (stderr, "sphere_chunkfinder failed\n");
            freegraph (&G);
            CCutil_stop_timer (&timer->all, 0);
            return 1;
        }
    } else {
        if (equiv_chunkfinder (&G, flags, &timer->locate, callback)) {
            fprintf (stderr, "equiv_chunkfinder failed\n");
            freegraph (&G);
            CCutil_stop_timer (&timer->all, 0);
            return 1;
        }
    }

    freegraph (&G);
    CCutil_stop_timer (&timer->all, 0);
    return 0;
}

static graph_chunk_old *graph_chunk_old_alloc (int ncount, int ecount)
{
    graph_chunk_old *c;
    int i;

    c = CC_SAFE_MALLOC (1, graph_chunk_old);
    if (!c) {
        printf ("out of memory in graph_chunk_old_alloc\n");
        return (graph_chunk_old *) NULL;
    }

    c->end0 = CC_SAFE_MALLOC (ecount, int);
    c->end1 = CC_SAFE_MALLOC (ecount, int);
    c->upper = CC_SAFE_MALLOC (ecount, int);
    c->lower = CC_SAFE_MALLOC (ecount, int);
    c->weight = CC_SAFE_MALLOC (ecount, double);
    c->members = CC_SAFE_MALLOC (ncount, int *);

    if (!c->end0 || !c->end1 || !c->upper || !c->lower
                 || !c->weight || !c->members) {
        printf ("out of memory in graph_chunk_old_alloc\n");
        return (graph_chunk_old *) NULL;
    }

    for (i = 0; i < ncount; i++)
        c->members[i] = (int *) NULL;

    c->ncount = ncount;
    c->ecount = ecount;

    return c;
}

static void graph_chunk_old_free (graph_chunk_old *c)
{
    int i;

    if (c) {
        CC_IFFREE (c->end0, int);
        CC_IFFREE (c->end1, int);
        CC_IFFREE (c->upper, int);
        CC_IFFREE (c->lower, int);
        CC_IFFREE (c->weight, double);
        if (c->members) {
            for (i = 0; i < c->ncount; i++) {
                CC_IFFREE (c->members[i], int);
            }
            CC_FREE (c->members, int *);
        }
        CC_FREE (c, graph_chunk_old);
    }
}

static int buildgraph (graph *G, int vcount, int ecount, int *elist,
                       double *elen)
{
    int k, i;

    G->elist = (edge *) NULL;
    G->vlist = (vertex *) NULL;
    G->ecount = 0;
    G->vcount = 0;


    G->vlist = CC_SAFE_MALLOC (vcount, vertex);
    G->elist = CC_SAFE_MALLOC (ecount, edge);
    G->supply = CC_SAFE_MALLOC (2*ecount, item);
    if (!G->vlist || !G->elist || !G->supply) {
        printf ("out of memory in getedges\n");
        G->vlist = (vertex *) NULL;
        G->elist = (edge *) NULL;
        G->supply = (item *) NULL;
        G->vcount = 0;
        G->ecount = 0;
        return 1;
    }

    G->vcount = vcount;

    for (i = 0, k = 0; i < ecount; i++)
        if (elen[i] > EPS3)
            k++;
    G->ecount = k;

    for (i = 0, k = 0; i < ecount; i++) {
        if (elen[i] > EPS3) {
            G->elist[k].end1 = G->vlist + elist[2*i];
            G->elist[k].end2 = G->vlist + elist[2*i+1];
            G->elist[k].weight = elen[i];
            k++;
        }
    }
    make_adjlists (G);

    return 0;
}

static void freegraph (graph *G)
{
    CC_IFFREE (G->vlist, vertex);
    CC_IFFREE (G->elist, edge);
    CC_IFFREE (G->supply, item);
    G->vcount = 0;
    G->ecount = 0;
}

static void permute_edgelist (graph *G, CCrandstate *rstate)
{
    int i, k;
    edge temp;

    for (i = G->ecount; i > 1; i--) {
        k = CCutil_lprand (rstate) % i;
        CC_SWAP (G->elist[i-1], G->elist[k], temp);
    }
    make_adjlists (G);
}

static void make_adjlists (graph *G)
{
    edge *e, *lastedge = G->elist + G->ecount;
    vertex *v, *lastvertex = G->vlist + G->vcount;
    item *s = G->supply;

    for (v = G->vlist; v != lastvertex; v++)
        v->adj = (item *) NULL;

    for (e = G->elist; e != lastedge; e++) {
        s->edgeptr = e;
        s->next = e->end1->adj;
        e->end1->adj = s;
        s++;

        s->edgeptr = e;
        s->next = e->end2->adj;
        e->end2->adj = s;
        s++;
    }
}

static int equiv_chunkfinder (graph *G, CCchunk_flag flags, 
        CCutil_timer *timer, CCchunk_chunk_callback *callback)
{
    edge *e, *lastedge = G->elist + G->ecount;
    vertex *v, *lastvertex = G->vlist + G->vcount;

    CCutil_start_timer (timer);
    
    for (e = G->elist; e != lastedge; e++) {
        e->link = e;
        e->rank = 0;
    }

    for (v = G->vlist; v != lastvertex; v++)
        v->into = (edge *) NULL;

    for (v = G->vlist; v != lastvertex; v++)
        merge_cycles (v);

    convert_lists (G);

    if (grab_all_chunks (G, flags, timer, callback)) {
        printf ("grab_all_chunks failed\n");
        CCutil_stop_timer (timer, 0);
        return 1;
    }
    CCutil_stop_timer (timer, 0);
    return 0;
}

static void merge_cycles (vertex *v)
{
    vertex *a, *b, *w;
    edge *e, *f, *g;
    item *pi, *pj, *pk;

    /* run through all choices of distinct neighbors a,b of v; notation e=va,
     * f=vb; whenever avb or some avbw is a cycle, put its edges in the same
     * equivalence class */

    /* fix a */
    for (pi = v->adj; pi; pi = pi->next) {
        e = pi->edgeptr;
        a = OTHEREND(e, v);

        /* mark all neighbors of a */
        for (pk = a->adj; pk; pk = pk->next) {
            g = pk->edgeptr;
            w = OTHEREND(g, a);
            w->into = g;
        }

        /* fix b */
        for (pj = pi->next; pj; pj = pj->next) {
            f = pj->edgeptr;
            if (cclass (e) == cclass (f))
                continue;
            b = OTHEREND(f, v);

            /* scan all neighbors of b */
            for (pk = b->adj; pk; pk = pk->next) {
                g = pk->edgeptr;
                if (g == f)
                    continue;
                w = OTHEREND(g, b);
                if (w == a) {
                    /* avb is a cycle */
                    merge_classes (e, f);
                    merge_classes (e, g);
                }
                if (w->into) {
                    /* avbw is a cycle */
                    merge_classes (e, f);
                    merge_classes (e, w->into);
                    merge_classes (e, g);
                }
            }
        }

        /* unmark neighbors of a */;

        for (pk = a->adj; pk; pk = pk->next) {
            pk->edgeptr->end1->into = (edge *) NULL;
            pk->edgeptr->end2->into = (edge *) NULL;
        }
    }
}

static edge *cclass (edge *e)
{
    edge *r, *s, *t;

    for (r = e; r->link != r; r = r->link);
    for (s = e; s->link != r; s = t) {
        t = s->link;
        s->link = r;
    }
    return r;
}

static void merge_classes (edge *e, edge *f)
{
    e = cclass (e);
    f = cclass (f);
    if (e == f)
        return;

    if (e->rank == f->rank)
        (e->rank)++;
    if (e->rank > f->rank)
        f->link = e;
    else
        e->link = f;
}

static void convert_lists (graph *G)
{
    edge *e, *f;
    edge *lastedge = G->elist + G->ecount;

    for (e = G->elist; e != lastedge; e++)
        e->status = (cclass (e) == e ? HEADER : NONHEADER);

    for (e = G->elist; e != lastedge; e++)
        if (e->status == HEADER) {
            e->link =  (edge *) NULL;
        }

    for (e = G->elist; e != lastedge; e++)
        if (e->status == NONHEADER) {
            f = e->link;
            e->link = f->link;
            f->link = e;
        }
}

static int grab_all_chunks (graph *G, CCchunk_flag flags,
        CCutil_timer *timer, CCchunk_chunk_callback *callback)
{
    edge *f;
    int i, count;
    vertex **chunkvertex = (vertex **) NULL;
    vertex *v, *lastvertex = G->vlist + G->vcount;
    edge *e, *lastedge = G->elist + G->ecount;
    graph_chunk_old *c = (graph_chunk_old *) NULL;
    int rval = 0;

    chunkvertex = CC_SAFE_MALLOC (flags.maxchunksize, vertex *);
    if (!chunkvertex) {
        fprintf (stderr, "Out of memory in grab_all_chunks\n");
        rval = -1;
        goto CLEANUP;
    }

    for (v = G->vlist; v != lastvertex; v++)
        v->name = OUT;

    for (e = G->elist; e != lastedge; e++) {
        if (e->status != NONHEADER) {
            count = 1;
            for (f = e;
                 f && ((unsigned) count) < ((unsigned) flags.maxchunksize);
                 f = f->link) {
                v = f->end1;
                if (v->name == OUT) {
                    v->name = count++;
                    chunkvertex[v->name] = v;
                }
                v = f->end2;
                if (v->name == OUT &&
                    ((unsigned) count) < ((unsigned) flags.maxchunksize)) {
                    v->name = count++;
                    chunkvertex[v->name] = v;
                }
            }

            if (count > MINCHUNKSIZE &&
                ((unsigned) count) < ((unsigned) flags.maxchunksize)) {
                if (!(flags.uncivilized)) {
                    enlarge_chunk (chunkvertex, &count, flags.maxchunksize);
                }
                c = grab_chunk (G, chunkvertex, count);
                if (!c) {
                    fprintf (stderr, "grab_chunk failed\n");
                    rval = -1;
                    goto CLEANUP;
                }

                if (c->ncount <= G->vcount) { /* there is an outside */
                    rval = checkout_chunk (c, timer, callback);
                    if (rval) {
                        fprintf (stderr, "checkout_chunk failed\n");
/* FAILURE could just be determinant
                        graph_chunk_old_free (c);
                        goto CLEANUP;
*/
                    }
                }
                graph_chunk_old_free (c);
                c = (graph_chunk_old *) NULL;
            }
            for (i = 1; i < count; i++) {
                chunkvertex[i]->name = OUT;
            }
        }
    }
    rval = 0;

  CLEANUP:
    CC_IFFREE (chunkvertex, vertex *);
    return rval;
}

static int dummy_grab_chunk (graph *G, int count, int *list,
        graph_chunk_old **c)
{
    vertex **chunkvertex = (vertex **) NULL;
    int i, rval = 0;

    chunkvertex = CC_SAFE_MALLOC (count + 1, vertex *);
    if (!chunkvertex) {
        fprintf (stderr, "out of memory in dummy_grab_chunk\n");
        rval = 1; goto CLEANUP;
    }
   
    for (i = 0; i < count; i++) {
        chunkvertex[i+1] = G->vlist + list[i];
        G->vlist[list[i]].name = i+1;
    }

    *c = grab_chunk (G, chunkvertex, count+1);

    for (i = 1; i <= count; i++) {
        chunkvertex[i]->name = OUT;
    }

CLEANUP:

    CC_IFFREE (chunkvertex, vertex *);
    return rval;
}

static graph_chunk_old *grab_chunk (graph *G, vertex **chunkvertex, int count)
{
    vertex *v;
    edge *e;
    item *pi;
    int i;
    double weight;
    int outside, lower, upper, cecount;
    int k = 0;
    graph_chunk_old *c = (graph_chunk_old *) NULL;

    if (count >= G->vcount) count = G->vcount-1;

    cecount = chunk_ecount (chunkvertex, count);
    c = graph_chunk_old_alloc (count, cecount);
    if (!c) {
        printf ("graph_chunk_old_alloc failed\n");
        return (graph_chunk_old *) NULL;
    }
    for (i = 1; i < count; i++) {
        v = chunkvertex[i];
        weight = 0.0;
        lower = outside = 0;
        upper = 2;

        for (pi = v->adj; pi; pi = pi->next) {
            e = pi->edgeptr;
            if ((e->end1->name != OUT) && (e->end2->name != OUT) &&
                (e->end1->name < count) && (e->end2->name < count)) {
                if (e->weight >= 1 - EPS4)
                    upper = 1;
                if (e->end2 == v) {
                    c->end0[k] = e->end1->name;
                    c->end1[k] = e->end2->name;
                    if (e->weight >= 1 - EPS4 && e->weight <= 1 + EPS4) {
                        c->lower[k] = 1;
                        c->upper[k] = 1;
                        c->weight[k] = 1.0;
                    } else {
                        c->lower[k] = 0;
                        c->upper[k] = 1;
                        c->weight[k] = e->weight;
                    }
                    k++;
                }
            } else {
                outside++;
                weight += e->weight;
                if (e->weight >= 1 - EPS4)
                    lower = 1;
            }
        }

        if (outside) {
            if (outside == 1)
                upper = 1;
            c->end0[k] = 0;
            c->end1[k] = i;
            c->weight[k] = weight;
            c->lower[k] = lower;
            c->upper[k] = upper;
            k++;
        }
    }

    for (i = 0; i < count; i++) {
        c->members[i] = CC_SAFE_MALLOC (2, int);
        if (!c->members[i]) {
            printf ("out of memory in grab_chunks\n");
            graph_chunk_old_free (c);
            return (graph_chunk_old *) NULL;
        }
        if (i == 0) {
            c->members[i][0] = -1;
        } else {
            c->members[i][0] = (int) (chunkvertex[i] - G->vlist);
            c->members[i][1] = -1;
        }
    }

    return c;
}

static int validate_chunk (graph_chunk_old *c)
{
    int i;

    for (i=0; i<c->ecount; i++) {
        if (c->weight[i] > 1.0 + EPS2 &&
            c->end0[i] != 0 && c->end1[i] != 0) {
            return 1;
        }
    }
    return 0;
}

static CCchunk_graph *strip_outside (graph_chunk_old *c)
{
    CCchunk_graph *cnew = (CCchunk_graph *) NULL;
    int i;
    int ecount;

    ecount = 0;
    for (i=0; i<c->ecount; i++) {
        if (c->end0[i] != 0 && c->end1[i] != 0) {
            ecount++;
        }
    }

    if (c->ncount <= 4 || ecount <= 1) {
        goto CLEANUP;
    }

    cnew = CCchunk_graph_alloc (c->ncount-1, ecount);
    if (!cnew) {
        fprintf (stderr, "CCchunk_graph_alloc failed\n");
        goto CLEANUP;
    }

    for (i=0; i<cnew->ncount; i++) {
        cnew->equality[i] = 1;
        cnew->members[i] = member_dup (c->members[i+1]);
        if (!cnew->members[i]) {
            fprintf (stderr, "member_dup failed\n");
            goto CLEANUP;
        }
    }

    ecount = 0;
    for (i=0; i<c->ecount; i++) {
        if (c->end0[i] == 0) {
            if (c->lower[i] != 0 || c->upper[i] != 2) {
                fprintf (stderr, "strip_outside, outside edge bounds %d %d\n",
                         c->lower[i], c->upper[i]);
                goto CLEANUP;
            }
            cnew->equality[c->end1[i]-1] = 0;
        } else if (c->end1[i] == 0) {
            if (c->lower[i] != 0 || c->upper[i] != 2) {
                fprintf (stderr, "strip_outside, outside edge bounds %d %d\n",
                         c->lower[i], c->upper[i]);
                goto CLEANUP;
            }
            cnew->equality[c->end0[i]-1] = 0;
        } else {
            cnew->end0[ecount] = c->end0[i]-1;
            cnew->end1[ecount] = c->end1[i]-1;
            cnew->weight[ecount] = c->weight[i];
            if (c->lower[i] == 0 && c->upper[i] == 1) {
                cnew->fixed[ecount] = -1;
            } else if (c->lower[i] == 0 && c->upper[i] == 0) {
                cnew->fixed[ecount] = 0;
            } else if (c->lower[i] == 1 && c->upper[i] == 1) {
                cnew->fixed[ecount] = 1;
            } else {
                fprintf (stderr, "strip_outside, non-outside edge bounds %d %d\n",
                         c->lower[i], c->upper[i]);
                goto CLEANUP;
            }
            ecount++;
        }
    }

    return cnew;

  CLEANUP:
    if (cnew) CCchunk_graph_free (cnew);
    return (CCchunk_graph *) NULL;
}

static int *member_dup (int *omem)
{
    int i;
    int *nmem = (int *) NULL;
    int cnt;

    for (i=0, cnt=0; omem[i] != -1; i++) cnt++;
    cnt++;

    nmem = CC_SAFE_MALLOC (cnt, int);
    if (!nmem) {
        fprintf (stderr, "Out of memory in member_dup\n");
        return (int *) NULL;
    }

    for (i=0; omem[i] != -1; i++) nmem[i] = omem[i];
    nmem[i] = -1;
    return nmem;
}


static int checkout_chunk (graph_chunk_old *c, CCutil_timer *timer,
        CCchunk_chunk_callback *callback)
{
    int rval;
    int faulty = 0;
    CCchunk_graph *cnew = (CCchunk_graph *) NULL;

/*
    printf ("Q"); fflush (stdout);
*/

    civilize_chunk (c);
    rval = validate_chunk (c);
    if (rval) return 0;

    cnew = strip_outside (c);
    if (cnew == (CCchunk_graph *) NULL) {
        rval = 0;
        goto CLEANUP;
    }

    CCutil_suspend_timer (timer);
    rval = (*callback->func) (cnew, &faulty, callback->u_data);
    CCutil_resume_timer (timer);
#ifdef DUMPBADCHUNKS
    if (rval) {
        int i;
        printf ("failed chunk:\n");
        printf ("%d %d\n", c->ncount, c->ecount);
        for (i=0; i<c->ecount; i++) {
            printf ("%d %d %.16f %d %d\n",c->end0[i], c->end1[i], c->weight[i],
                    c->lower[i], c->upper[i]);
        }
        fflush (stdout);
    }
#endif
 CLEANUP:
    if (cnew) CCchunk_graph_free (cnew);
    return rval;
}

static void enlarge_chunk (vertex **chunkvertex, int *count, int maxchunksize)
{
    vertex *v;
    edge *e;
    item *pi;
    int i;

    for (i = 1; i < *count && *count < maxchunksize; i++) {
        v = chunkvertex[i];
        for (pi = v->adj; pi && *count < maxchunksize; pi = pi->next) {
            e = pi->edgeptr;
            if ( ((e->end1->name == OUT) || (e->end2->name == OUT))
                  && (e->weight >= 1 - EPS4)){
                if (e->end1->name == OUT) {
                    e->end1->name = (*count)++;
                    chunkvertex[e->end1->name] = e->end1;
                } else {
                    e->end2->name = (*count)++;
                    chunkvertex[e->end2->name] = e->end2;
                }
            }

        }
    }
}

static void civilize_chunk (graph_chunk_old *c)
{
    int i;

    for (i=0; i<c->ecount; i++) {
        if (c->end0[i] == 0 || c->end1[i] == 0) {
            if (c->lower[i] == 1) c->lower[i] = 0;
            if (c->upper[i] == 1) c->upper[i] = 2;
        }
    }
}

static int chunk_ecount (vertex **chunkvertex, int count)
{
    edge *e;
    vertex *v;
    item *pi;
    int ecount = 0;
    int outside, i;

    for (i = 1; i < count; i++) {
        v = chunkvertex[i];
        outside = 0;
        for (pi = v->adj; pi; pi = pi->next) {
            e = pi->edgeptr;
            if ((e->end1->name != OUT) && (e->end2->name != OUT) &&
                (e->end1->name < count) && (e->end2->name < count)) {
                if (e->end2 == v)
                    ecount++;
            } else {
                outside++;
            }
        }
        if (outside)
            ecount++;
    }

    return ecount;
}

static int sphere_chunkfinder (graph *G, double eps, CCchunk_flag flags,
        CCutil_timer *timer, CCchunk_chunk_callback *callback)
{
    vertex *v, *lastvertex = G->vlist + G->vcount;
    graph_chunk_old *c;
    int rval;

    CCutil_start_timer (timer);
    
    for (v = G->vlist; v != lastvertex; v++)
        v->name = OUT;

    for (v = G->vlist; v != lastvertex; v++) {
        if (v->adj) {
            c = get_sphere (G, v, eps, flags);
            if (!c) {
/* failure could just be size */
                CCutil_stop_timer (timer, 0);
                return 0;
            }

            if (c->ncount <= G->vcount) { /* there is an outside */
                rval = checkout_chunk (c, timer, callback);
                if (rval) {
                    fprintf (stderr, "checkout_chunk failed\n");
/* failure could just be determinant
                    graph_chunk_old_free (c);
                    CCutil_stop_timer (timer, 0);
                    return rval;
*/
                }
            }
            graph_chunk_old_free (c);
        }
    }
    CCutil_stop_timer (timer, 0);
    return 0;
}

static int dummy_chunkfinder (graph *G, double eps, CCchunk_flag flags,
        CCutil_timer *timer, CCchunk_chunk_callback *callback)
{
    vertex *v, *lastvertex = G->vlist + G->vcount;
    graph_chunk_old *c;
    int i, rval = 0;

    CCutil_start_timer (timer);

    for (v = G->vlist; v != lastvertex; v++) v->name = OUT;
    for (v = G->vlist; v != lastvertex; v++) v->val = 0.0;
    for (v = G->vlist; v != lastvertex; v++) v->deg = 0;

    for (i = 0; i < G->ecount; i++) {
        if (G->elist[i].weight > eps) {
            G->elist[i].end1->deg++;
            G->elist[i].end2->deg++;
        }
    }

    for (v = G->vlist; v != lastvertex; v++) {
        if (v->deg) {
            c = get_card_sphere (G, v, eps, flags);
            if (!c) {
                CCutil_stop_timer (timer, 0);
                return 0;
            }

            if (c->ncount <= G->vcount) { /* there is an outside */
                rval = checkout_chunk (c, timer, callback);
                if (rval) {
                    fprintf (stderr, "checkout_chunk failed\n");
                }
            }
            graph_chunk_old_free (c);
        }
    }
    

    CCutil_stop_timer (timer, 0);
    return 0;
}

static int permute_chunkfinder (graph *G, double eps, CCchunk_flag flags,
        CCutil_timer *timer, CCchunk_chunk_callback *callback,
        CCrandstate *rstate)
{
    vertex *v, *lastvertex = G->vlist + G->vcount;
    graph_chunk_old *c;
    int rval;

    CCutil_start_timer (timer);

    for (v = G->vlist; v != lastvertex; v++) v->name = OUT;

    permute_edgelist (G, rstate);

    for (v = G->vlist; v != lastvertex; v++) {
        if (v->adj) {
            c = get_sphere (G, v, eps, flags);
            if (!c) {
                CCutil_stop_timer (timer, 0);
                return 0;
            }

            if (c->ncount <= G->vcount) { /* there is an outside */
                rval = checkout_chunk (c, timer, callback);
                if (rval) {
                    fprintf (stderr, "checkout_chunk failed\n");
                }
            }
            graph_chunk_old_free (c);
        }
    }

    CCutil_stop_timer (timer, 0);
    return 0;
}

static int weighted_chunkfinder (graph *G, double eps, CCchunk_flag flags,
        CCutil_timer *timer, CCchunk_chunk_callback *callback)
{
    vertex *v, *lastvertex = G->vlist + G->vcount;
    graph_chunk_old *c;
    int rval;

    CCutil_start_timer (timer);

    for (v = G->vlist; v != lastvertex; v++) v->name = OUT;
    for (v = G->vlist; v != lastvertex; v++) v->val = 0.0;

    for (v = G->vlist; v != lastvertex; v++) {
        if (v->adj) {
            c = get_weighted_sphere (G, v, eps, flags);
            if (!c) {
                CCutil_stop_timer (timer, 0);
                return 0;
            }

            if (c->ncount <= G->vcount) { /* there is an outside */
                rval = checkout_chunk (c, timer, callback);
                if (rval) {
                    fprintf (stderr, "checkout_chunk failed\n");
                }
            }
            graph_chunk_old_free (c);
        }
    }
    
    CCutil_stop_timer (timer, 0);
    return 0;
}

static graph_chunk_old *get_sphere (graph *G, vertex *v, double eps,
        CCchunk_flag flags)
{
    vertex *u, *w;
    edge *e;
    item *pi;
    int i, count;
    vertex **chunkvertex = (vertex **) NULL;
    graph_chunk_old *c = (graph_chunk_old *) NULL;

    chunkvertex = CC_SAFE_MALLOC (flags.maxchunksize, vertex *);
    if (!chunkvertex) {
        fprintf (stderr, "Out of memory in get_sphere\n");
        goto CLEANUP;
    }

    v->name = 1;
    chunkvertex[1] = v;
    count = 2;

    for (i = 1;
         i < count && ((unsigned) count) < ((unsigned) flags.spheresize);
         i++) {
        u = chunkvertex[i];
        for (pi = u->adj;
             pi && ((unsigned) count) < ((unsigned) flags.spheresize);
             pi = pi->next) {
            e = pi->edgeptr;
            if (e->weight > eps) {
                w = OTHEREND(e, u);
                if (w->name == OUT) {
                    w->name = count++;
                    chunkvertex[w->name] = w;
                }
            }
        }
    }

    if (!flags.uncivilized) {
        enlarge_chunk (chunkvertex, &count, flags.maxchunksize);
    }

    if (count > MINCHUNKSIZE) {
        c = grab_chunk (G, chunkvertex, count);
        if (!c) {
            fprintf (stderr, "grab_chunk failed\n");
            goto CLEANUP;
        }

    }
    for (i = 1; i < count; i++)
        chunkvertex[i]->name = OUT;

  CLEANUP:
    CC_IFFREE (chunkvertex, vertex *);
    return c;
}

static graph_chunk_old *get_weighted_sphere (graph *G, vertex *v,
        double eps, CCchunk_flag flags)
{
    vertex *u, *w;
    vertex *maxn = (vertex *) NULL;
    edge *e;
    item *pi;
    int i, j, count;
    double max;
    vertex **chunkvertex = (vertex **) NULL;
    graph_chunk_old *c = (graph_chunk_old *) NULL;

    chunkvertex = CC_SAFE_MALLOC (flags.maxchunksize, vertex *);
    if (!chunkvertex) {
        fprintf (stderr, "Out of memory in get_sphere\n");
        goto CLEANUP;
    }

    v->name = 1;
    chunkvertex[1] = v;
    count = 2;

    for (i = 1;
         i < count && ((unsigned) count) < ((unsigned) flags.spheresize);
         i++) {
        u = chunkvertex[i];
        for (pi = u->adj; pi; pi = pi->next) {
            e = pi->edgeptr;
            if (e->weight > eps) {
                w = OTHEREND (e, u);
                w->val += e->weight;
            }
        }
        max = -1.0;
        for (j = 1; j <= i; j++) {
            u = chunkvertex[j];
            for (pi = u->adj; pi; pi = pi->next) {
                e = pi->edgeptr;
                if (e->weight > eps) {
                    w = OTHEREND (e, u);
                    if (w->name == OUT && w->val > max)  {
                        max = w->val;
                        maxn = w;
                    }
                }
            }
        }
        if (max != -1.0) {
            maxn->name = count++;
            chunkvertex[maxn->name] = maxn;
        }
    }

    if (!flags.uncivilized) {
        enlarge_chunk (chunkvertex, &count, flags.maxchunksize);
    }

    if (count > MINCHUNKSIZE) {
        c = grab_chunk (G, chunkvertex, count);
        if (!c) {
            fprintf (stderr, "grab_chunk failed\n");
            goto CLEANUP;
        }
    }

    for (i = 1; i < count; i++) {
        u = chunkvertex[i];
        u->name = OUT;
        u->val = 0.0;
        for (pi = u->adj; pi; pi = pi->next) {
            e = pi->edgeptr;
            w = OTHEREND (e, u);
            w->val = 0.0;
        }
    }

CLEANUP:

    CC_IFFREE (chunkvertex, vertex *);
    return c;
}

static graph_chunk_old *get_card_sphere (graph *G, vertex *v, double eps,
        CCchunk_flag flags)
{
    vertex *u, *w, *maxn;
    edge *e;
    item *pi;
    int i, j, count;
    double max;
    vertex **chunkvertex = (vertex **) NULL;
    graph_chunk_old *c = (graph_chunk_old *) NULL;

    chunkvertex = CC_SAFE_MALLOC (flags.maxchunksize, vertex *);
    if (!chunkvertex) {
        fprintf (stderr, "Out of memory in get_sphere\n");
        goto CLEANUP;
    }

    v->name = 1;
    chunkvertex[1] = v;
    count = 2;

    for (i = 1;
         i < count && ((unsigned) count) < ((unsigned) flags.spheresize);
         i++) {
        u = chunkvertex[i];
        for (pi = u->adj; pi; pi = pi->next) {
            e = pi->edgeptr;
            if (e->weight > eps) {
                w = OTHEREND (e, u);
                w->val += 1.0;
            }
        }
        max = (double) -(G->vcount);
        maxn = (vertex *) NULL;
        for (j = 1; j <= i; j++) {
            u = chunkvertex[j];
            for (pi = u->adj; pi; pi = pi->next) {
                e = pi->edgeptr;
                if (e->weight > eps) {
                    w = OTHEREND (e, u);
                    if (w->name == OUT && (2 * w->val) - (double) w->deg > max)  {
                        max = 2 * w->val - (double) w->deg;
                        maxn = w;
                    }
                }
            }
        }
        if (maxn != (vertex *) NULL) {
            maxn->name = count++;
            chunkvertex[maxn->name] = maxn;
        }
    }

    if (!flags.uncivilized) {
        enlarge_chunk (chunkvertex, &count, flags.maxchunksize);
    }

    if (count > MINCHUNKSIZE) {
        c = grab_chunk (G, chunkvertex, count);
        if (!c) {
            fprintf (stderr, "grab_chunk failed\n");
            goto CLEANUP;
        }
    }

    for (i = 1; i < count; i++) {
        u = chunkvertex[i];
        u->name = OUT;
        u->val = 0.0;
        for (pi = u->adj; pi; pi = pi->next) {
            e = pi->edgeptr;
            w = OTHEREND (e, u);
            w->val = 0.0;
        }
    }

CLEANUP:

    CC_IFFREE (chunkvertex, vertex *);
    return c;
}

static void CC_UNUSED dump_graph (graph *G)
{
    int i;
    item *a;

    printf ("%d %d\n", G->vcount, G->ecount);
    for (i=0; i<G->ecount; i++) {
        printf ("%d %d %.16f\n", (int) (G->elist[i].end1-G->vlist),
                (int) (G->elist[i].end2-G->vlist),G->elist[i].weight);
    }
    printf ("VERTEX ADJ:\n");
    for (i=0; i<G->vcount; i++) {
        printf ("%d (%d):", i, G->vlist[i].name);
        for (a=G->vlist[i].adj; a; a = a->next) {
            printf (" %d %.16f", a->edgeptr->end1 == &G->vlist[i] ?
                    (int) (a->edgeptr->end2-G->vlist) :
                    (int) (a->edgeptr->end1-G->vlist),
                    a->edgeptr->weight);
        }
        printf ("\n");
    }
}

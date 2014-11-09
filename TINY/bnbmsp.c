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
/*  CCtiny_bnb_msp (int nnodes, int nedges, int *elist, int *weight,        */
/*      int depot, int *lbound, int *ubound, double *objlimit,              */
/*      int objdir, double *objval, int *xsol, int searchlimit)             */
/*    MISSING                                                               */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tinytsp.h"

/* the node which can have degree > 2 */
#define OUTSIDE (0)
#define OUTSIDEEDGE(n0,n1) ((n0) == OUTSIDE || (n1) == OUTSIDE)

#define BIGDOUBLE (1e30)
#undef  USE_DUMBO

#ifdef USE_DUMBO
#define TINYBRANCH_WEIGHT 10.0
#define TINYBRANCH_VAL(x,y) (((x)<(y)) ? TINYBRANCH_VAL2((y),(x)) \
                                       : TINYBRANCH_VAL2((x),(y)))
#define TINYBRANCH_VAL2(x,y) (((x)+(y)*TINYBRANCH_WEIGHT) / \
                              (TINYBRANCH_WEIGHT+1.0))
#endif

#define SELECT_INFEASIBLE 1
#define SELECT_TOUR 2

/* #define DEBUG 1 */
#undef TIMINGS

#define SEARCHLIMIT (1<<30)

typedef struct node {
    int availdeg;
    int deg;
/* OUTSIDE's pathend == OUTSIDE */
    int pathend;
    struct edge *adjvec;
    int adjlist;
} node;

typedef struct edge {
    int next;
    int prev;
    int weight;
    int lo;
    int hi;
    int value;
    int bestvalue;
} edge;

typedef struct histent {
    int from;
    int to;
    int pfrom;
    int pto;
} histent;

typedef struct tspsearch {
    int nnodes;
    node *nodes;
    histent *searchhist;
    int histtop;
    double currentval;
    double bestval;
    int searchcount;
    int searchlimit;
} tspsearch;

static int
    checkout_node (tspsearch *s, int i),
    build_paths (tspsearch *s),
    set_edge(tspsearch *s, int n0, int n1, int v),
    select_edge (tspsearch *s, int *n0, int *n1);

static void
    build_adj (tspsearch *s),
    search (tspsearch *s),
    unset_edge(tspsearch *s, histent *h),
    savetour (tspsearch *s),
    unroll_stack (tspsearch *s, int smark);

static double
    lowerbound (tspsearch *s);

int CCtiny_bnb_msp (int nnodes, int nedges, int *elist, int *weight, int depot,
        int *lbound, int *ubound, double *objlimit, int objdir,
        double *objval, int *xsol, int searchlimit)
{
    tspsearch s;
    int i, j;
    int rval;
    int n0, n1;
    double initlimit;
#ifdef TIMINGS
    double sz = CCutil_zeit();
#endif

    if (depot != OUTSIDE) {
        fprintf (stderr, "bnbtsp doesn't know how to handle depot %d\n",
                 depot);
        return CC_TINYTSP_ERROR;
    }

    rval = CC_TINYTSP_ERROR; /* Failures in this part are from out of memory */
    s.nnodes = nnodes;
    s.nodes = (node *) NULL;
    s.searchhist = (histent *) NULL;
    s.histtop = 0;
    s.currentval = 0.0;
    if (objlimit) {
        initlimit = objdir * (*objlimit);
    } else {
        initlimit = BIGDOUBLE;
    }
    s.bestval = initlimit;

    if (searchlimit < 0) {
        s.searchlimit = SEARCHLIMIT;
    } else {
        s.searchlimit = searchlimit;
    }

    s.nodes = CC_SAFE_MALLOC (nnodes, node);
    if (!s.nodes) goto CLEANUP;
    for (i=0; i<nnodes; i++) {
        s.nodes[i].adjvec = (edge *) NULL;
    }
    s.searchhist = CC_SAFE_MALLOC (nnodes * nnodes, histent);
    if (!s.searchhist) goto CLEANUP;
    for (i=0; i<nnodes; i++) {
        s.nodes[i].adjvec = CC_SAFE_MALLOC (nnodes, edge);
        if (!s.nodes[i].adjvec) goto CLEANUP;
    }

    for (i=0; i<nnodes; i++) {
        s.nodes[i].pathend = i;
        for (j=0; j<nnodes; j++) {
            s.nodes[i].adjvec[j].weight = 0;
            s.nodes[i].adjvec[j].lo = 0;
            s.nodes[i].adjvec[j].hi = 0;
            s.nodes[i].adjvec[j].value = -1;
            s.nodes[i].adjvec[j].bestvalue = -1;
        }
    }

    for (i=0; i<nedges; i++) {
        n0 = elist[2*i];
        n1 = elist[2*i+1];
        s.nodes[n0].adjvec[n1].weight = objdir * weight[i];
        s.nodes[n1].adjvec[n0].weight = objdir * weight[i];
        s.nodes[n0].adjvec[n1].lo = lbound[i];
        s.nodes[n1].adjvec[n0].lo = lbound[i];
        s.nodes[n0].adjvec[n1].hi = ubound[i];
        s.nodes[n1].adjvec[n0].hi = ubound[i];
    }

    build_adj (&s);

#ifdef DEBUG
    printf ("Performing initial checks\n");
#endif
    if (build_paths (&s)) {
        rval = CC_TINYTSP_INFEASIBLE;
        goto CLEANUP;
    }

    for (i=0; i<nnodes; i++) {
        if (checkout_node (&s, i)) {
            rval = CC_TINYTSP_INFEASIBLE;
            goto CLEANUP;
        }
        if (s.nodes[i].deg == 1 && i != OUTSIDE &&
            s.nodes[i].pathend != OUTSIDE &&
            s.nodes[i].adjvec[s.nodes[i].pathend].value == -1) {
#ifdef DEBUG
#if DEBUG>1
            printf ("node %d pathend, forcing %d to 0\n",
                    i,s.nodes[i].pathend);
#endif
#endif
            if (set_edge (&s, i, s.nodes[i].pathend, 0)) {
                rval = CC_TINYTSP_INFEASIBLE;
                goto CLEANUP;
            }
        }
    }

#ifdef DEBUG
    printf ("Beginning search\n");
#endif
    s.searchcount = 0;
    search (&s);
#ifdef DEBUG
    printf ("Search finished\n");
#endif

    if (s.searchcount >= s.searchlimit) {
#ifdef DEBUG
        fprintf (stderr, "search node limit %d exceeded\n", s.searchlimit);
        fprintf (stderr, "SNLE obj:\n");
        fprintf (stderr, "%d %d\n", nnodes, nedges);
        for (i=0; i<nedges; i++) {
            fprintf (stderr, "%d %d %d %d %d\n", elist[2*i], elist[2*i+1],
                     lbound[i], ubound[i], objdir * weight[i]);
        }
#endif
#ifdef TIMINGS
        printf ("BNB Search Limit Nodes: %d   Time %.3f\n",
                s.searchcount, CCutil_zeit() - sz); fflush (stdout);
#endif
        rval = CC_TINYTSP_SEARCHLIMITEXCEEDED;
        goto CLEANUP;
    }

#ifdef TIMINGS
    printf ("BNB Search in Nodes: %d   Time %.3f\n",
            s.searchcount, CCutil_zeit() - sz); fflush (stdout);
#endif
#ifdef DEBUG
    printf ("bnbtsp search finished in %d nodes\n", s.searchcount);
    fflush (stdout);
#endif

    if (objval) *objval = objdir * s.bestval;
    if (s.bestval >= initlimit) {
        rval = CC_TINYTSP_INFEASIBLE;
        goto CLEANUP;
    }

    if (xsol) {
        for (i=0; i<nedges; i++) {
            xsol[i] = s.nodes[elist[2*i]].adjvec[elist[2*i+1]].bestvalue;
        }
    }
    rval = 0;

  CLEANUP:
    if (s.nodes) {
        for (i=0; i<nnodes; i++) {
            CC_IFFREE (s.nodes[i].adjvec, edge);
        }
        CC_FREE (s.nodes, node);
    }
    CC_IFFREE (s.searchhist, histent);
    return rval;
}

static int checkout_node (tspsearch *s, int i)
{
    node *nodes = s->nodes;
    int j;

    if (nodes[i].availdeg < 2 || (i != OUTSIDE && nodes[i].deg > 2)) return 1;
    if (nodes[i].deg < 2) {
        if (nodes[i].availdeg == 2) {
            while ((j = nodes[i].adjlist) != -1) {
#ifdef DEBUG
#if DEBUG>1
                printf ("node %d avail 2, forcing %d to %d\n",
                        i,j,nodes[i].adjvec[j].hi);
#endif
#endif
                if (set_edge (s, i, j, nodes[i].adjvec[j].hi)) return 1;
            }
        }
    }
    j = nodes[i].adjlist;
    if (j != -1) {
        if (nodes[i].adjvec[j].next == -1) {
            if (i == OUTSIDE) {
                if (nodes[i].deg & 1) {
#ifdef DEBUG
#if DEBUG>1
                    printf ("node %d single adj + odd deg, forcing %d to %d\n",
                             i,j,nodes[i].adjvec[j].lo + 1);
#endif
#endif
                    if (set_edge (s, i, j, nodes[i].adjvec[j].lo + 1))
                        return 1;
                } else if (nodes[i].adjvec[j].hi - nodes[i].adjvec[j].lo
                           == 1) {
#ifdef DEBUG
#if DEBUG>1
                    printf ("node %d single adj + even deg, forcing %d to %d\n",
                             i,j,nodes[i].adjvec[j].lo);
#endif
#endif
                    if (set_edge (s, i, j, nodes[i].adjvec[j].lo)) return 1;
                }
            } else {
#ifdef DEBUG
#if DEBUG>1
                printf ("node %d single adj, forcing %d to %d\n",
                        i,j,2-nodes[i].deg);
#endif
#endif

                if (set_edge (s, i, j, 2 - nodes[i].deg)) return 1;
            }
        }
    }
    if (i != OUTSIDE && nodes[i].deg == 2) {
        while ((j = nodes[i].adjlist) != -1) {
#ifdef DEBUG
#if DEBUG>1
            printf ("node %d deg 2, forcing %d to %d\n",
                    i,j,nodes[i].adjvec[j].lo);
#endif
#endif
            if (set_edge (s, i, j, nodes[i].adjvec[j].lo)) return 1;
        }
    }

    return 0;
}

static int build_paths (tspsearch *s)
{
    int i, j;
    int pi, pj;
    int nnodes = s->nnodes;
    node *nodes = s->nodes;

    for (i=0; i<nnodes; i++) {
        nodes[i].pathend = i;
    }
    for (i=0; i<nnodes; i++) {
        for (j=i+1; j<nnodes; j++) {
            if (s->nodes[i].adjvec[j].lo > 0) {
                if (nodes[i].pathend == j) return 1;
                pi = nodes[i].pathend;
                pj = nodes[j].pathend;
                if (i != OUTSIDE) {
                    nodes[i].pathend = pj;
                    if (pi != OUTSIDE && pi != i) nodes[pi].pathend = pj;
                }
                if (j != OUTSIDE) {
                    nodes[j].pathend = pi;
                    if (pj != OUTSIDE && pj != j) nodes[pj].pathend = pi;
                }
            }
        }
    }
    return 0;
}

static void build_adj (tspsearch *s)
{
    int i, j;
    int p;
    int nnodes = s->nnodes;
    node *nodes = s->nodes;

    for (i=0; i<nnodes; i++) {
        p = -1;
        nodes[i].adjlist = -1;
        nodes[i].deg = 0;
        nodes[i].availdeg = 0;
        for (j=0; j<nnodes; j++) {
            nodes[i].deg += nodes[i].adjvec[j].lo;
            nodes[i].availdeg += nodes[i].adjvec[j].hi;
            if (nodes[i].adjvec[j].lo < nodes[i].adjvec[j].hi) {
                nodes[i].adjvec[j].prev = p;
                nodes[i].adjvec[j].next = -1;
                if (p == -1) {
                    nodes[i].adjlist = j;
                } else {
                    nodes[i].adjvec[p].next = j;
                }
                p = j;
            } else {
                nodes[i].adjvec[j].value = nodes[i].adjvec[j].lo;
                if (i <= j) { /* only add in each edge once */
                    s->currentval += nodes[i].adjvec[j].value *
                        nodes[i].adjvec[j].weight;
                }
            }
        }
    }
}

/* minimization */
static void search (tspsearch *s)
{
    int n0, n1;
    int smark;
    int i, lo, hi;
    int val;

    s->searchcount++;
    if (s->searchcount >= s->searchlimit) return;

    if (s->currentval + lowerbound(s) >= s->bestval) {
#ifdef DEBUG
        printf ("search, currentval %.0f lowerbound %.0f, bestval %.0f, returning\n",
                s->currentval, lowerbound(s), s->bestval);
#endif
        return;
    }
    val = select_edge (s, &n0, &n1);
    if (val == SELECT_TOUR) {
        if (s->currentval + lowerbound(s) < s->bestval) {
            savetour(s);
        }
        return;
    } else if (val == SELECT_INFEASIBLE) {
        return;
    }

    smark = s->histtop;
    lo = s->nodes[n0].adjvec[n1].lo;
    hi = s->nodes[n0].adjvec[n1].hi;
    if (s->nodes[n0].adjvec[n1].weight > 0) {
        for (i = lo; i<=hi; i++) {
#ifdef DEBUG
            printf ("branching %d-%d to %d (lo %d hi %d)\n",
                    n0,n1,i,lo,hi);
#endif
            if (!set_edge (s, n0, n1, i)) {
                search (s);
            }
            unroll_stack (s, smark);
        }
    } else {
        for (i = hi; i>=lo; i--) {
#ifdef DEBUG
            printf ("branching %d-%d to %d (lo %d hi %d)\n",
                    n0,n1,i,lo,hi);
#endif
            if (!set_edge (s, n0, n1, i)) {
                search (s);
            }
            unroll_stack (s, smark);
        }
    }
#ifdef DEBUG
    printf ("done branching on %d-%d (lo %d hi %d)\n",
            n0,n1,lo,hi);
#endif
    return;
}

/* 1 means contradiction found */
static int set_edge(tspsearch *s, int n0, int n1, int v)
{
    int p0, p1;
    edge *e;
    node *nodes = s->nodes;
    int lo, hi;

#ifdef DEBUG
#if DEBUG>2
    printf ("set_edge %d-%d = %d\n",n0,n1,v);
#endif
#endif
    if (n1 == OUTSIDE) {
        n1 = n0;
        n0 = OUTSIDE;
    }
    if (nodes[n0].adjvec[n1].value != -1) {
        if (nodes[n0].adjvec[n1].value == v) return 0;
        else return 1;
    }
    lo = nodes[n0].adjvec[n1].lo;
    hi = nodes[n0].adjvec[n1].hi;
    if (n0 != OUTSIDE && nodes[n0].deg - lo + v > 2) return 1;
    if (nodes[n1].deg - lo + v > 2) return 1;
    if (nodes[n0].availdeg - hi + v < 2) return 1;
    if (nodes[n1].availdeg - hi + v < 2) return 1;
    p0 = nodes[n0].pathend;
    p1 = nodes[n1].pathend;
    if (v >= 1 && n0 != OUTSIDE && n1 != OUTSIDE && p0 == n1) return 1;

    nodes[n0].deg += v - lo;
    nodes[n1].deg += v - lo;
    nodes[n0].availdeg += v - hi;
    nodes[n1].availdeg += v - hi;

    if (v >= 1) { /* The only situation in which pathends matter */
        if (p0 != OUTSIDE) nodes[p0].pathend = p1;
        if (p1 != OUTSIDE) nodes[p1].pathend = p0;
    }

    e = &nodes[n0].adjvec[n1];
    if (e->next != -1) nodes[n0].adjvec[e->next].prev = e->prev;
    if (e->prev == -1) {
        nodes[n0].adjlist = e->next;
    } else {
        nodes[n0].adjvec[e->prev].next = e->next;
    }
    e->value = v;

    e = &nodes[n1].adjvec[n0];
    if (e->next != -1) nodes[n1].adjvec[e->next].prev = e->prev;
    if (e->prev == -1) {
        nodes[n1].adjlist = e->next;
    } else {
        nodes[n1].adjvec[e->prev].next = e->next;
    }
    e->value = v;

    s->currentval += v * e->weight;

    s->searchhist[s->histtop].from = n0;
    s->searchhist[s->histtop].to = n1;
    s->searchhist[s->histtop].pfrom = p0;
    s->searchhist[s->histtop].pto = p1;
    s->histtop++;

    if (v >= 1) {
        if (p0 != OUTSIDE && p1 != OUTSIDE &&
            (p0 != n0 || p1 != n1)) {
#ifdef DEBUG
#if DEBUG>1
            printf ("Avoiding subtour, forcing %d-%d to 0\n",
                    p0,p1);
#endif
#endif
            if (set_edge (s, p0, p1, 0)) return 1;
        }
    }
    checkout_node (s, n0);
    checkout_node (s, n1);

    return 0;
}

static void unset_edge(tspsearch *s, histent *h)
{
    int n0 = h->from;
    int n1 = h->to;
    edge *e;
    node *nodes = s->nodes;
    int v = nodes[n0].adjvec[n1].value;
    int lo, hi;

    if (v == -1) return;

#ifdef DEBUG
#if DEBUG>2
    printf ("unset_edge %d-%d from %d\n",n0,n1,v);
#endif
#endif
    lo = nodes[n0].adjvec[n1].lo;
    hi = nodes[n0].adjvec[n1].hi;
    nodes[n0].deg += lo - v;
    nodes[n1].deg += lo - v;
    nodes[n0].availdeg += hi - v;
    nodes[n1].availdeg += hi - v;

    if (v >= 1) { /* The only situation in which pathends matter */
        if (n0 != OUTSIDE) nodes[n0].pathend = h->pfrom;
        if (n1 != OUTSIDE) nodes[n1].pathend = h->pto;
        if (h->pfrom != OUTSIDE && h->pfrom != n0)
            nodes[h->pfrom].pathend = n0;
        if (h->pto != OUTSIDE && h->pto != n1) nodes[h->pto].pathend = n1;
    }

    e = &nodes[n0].adjvec[n1];
    if (e->next != -1) nodes[n0].adjvec[e->next].prev = n1;
    if (e->prev == -1) {
        nodes[n0].adjlist = n1;
    } else {
        nodes[n0].adjvec[e->prev].next = n1;
    }
    e->value = -1;

    e = &nodes[n1].adjvec[n0];
    if (e->next != -1) nodes[n1].adjvec[e->next].prev = n0;
    if (e->prev == -1) {
        nodes[n1].adjlist = n0;
    } else {
        nodes[n1].adjvec[e->prev].next = n0;
    }
    e->value = -1;

    s->currentval -= v * e->weight;
}

#ifdef USE_DUMBO

static int select_edge (tspsearch *s, int *n0, int *n1)
{
    int i, j;
    int bestn0, bestn1;
    double bestval;
    node *nodes = s->nodes;
    int nnodes = s->nnodes;
    edge *adjvec;
    int nforced;
    double val;
    double val0;
    double val1;
    int infeas0;
    int infeas1;
    int smark;

    do {
        nforced = 0;
        bestn0 = bestn1 = -1;
        bestval = 0.0;
        for (i=1; i<nnodes; i++) {
            adjvec = nodes[i].adjvec;
            for (j=nodes[i].adjlist; j != -1; j = adjvec[j].next) {
                if (j != OUTSIDE && adjvec[j].value == -1) {
/* adjvec[j].value is tested since the loop can change the list */
/* since set_edge and unset_edge preserve list order, and don't */
/* change the next/prev of the deleted items, this is safe */
                    smark = s->histtop;
                    infeas0 = set_edge (s, i, j, 0);
                    if (infeas0 == 0) {
                        val0 = s->currentval + lowerbound(s);
                    }
                    unroll_stack (s, smark);
                    if (infeas0 == 0 && val0 < s->bestval) {
                        infeas1 = set_edge (s, i, j, 1);
                        if (infeas1 == 0) {
                            val1 = s->currentval + lowerbound(s);
                        }
                        unroll_stack (s, smark);
                    }

#ifdef DEBUG
#if DEBUG>1
                    printf ("eval edge %d-%d =  %.0f %.0f = %.0f\n",
                            i,j,infeas0?-999.0:val0,infeas1?-999.0:val1,
                            TINYBRANCH_VAL (val0, val1));
#endif
#endif
                    if ((infeas0 || val0 >= s->bestval) &&
                        (infeas1 || val1 >= s->bestval)) {
                        return SELECT_INFEASIBLE;
                    } else if (infeas0 || val0 >= s->bestval) {
                        if (set_edge (s, i, j, 1)) {
                            return SELECT_INFEASIBLE;
                        }
                        nforced++;
                    } else if (infeas1 || val1 >= s->bestval) {
                        if (set_edge (s, i, j, 0)) {
                            return SELECT_INFEASIBLE;
                        }
                        nforced++;
                    } else {
                        val = TINYBRANCH_VAL(val0, val1);
                        if (bestn0 == -1 || val > bestval) {
                            bestn0 = i;
                            bestn1 = j;
                            bestval = val;
                        }
                    }
                }
            }
        }
#ifdef DEBUG
        printf ("%d edges forced (bound %.0f val %.0f)\n",
                nforced, s->currentval + lowerbound(s), bestval);
        fflush (stdout);
#endif
    } while (nforced);

    if (bestn0 == -1) {
        return SELECT_TOUR;
    }
    *n0 = bestn0;
    *n1 = bestn1;
    return 0;
}

#else /* USE_DUMBO */

static int select_edge (tspsearch *s, int *n0, int *n1)
{
    int i, j;
    int bestn0 = -1, bestn1 = -1;
    int bestval;
    edge *adjvec;
    node *nodes = s->nodes;
    int nnodes = s->nnodes;

    bestval = -1;

    for (i=0; i<nnodes; i++) {
        if (i != OUTSIDE) {
            adjvec = nodes[i].adjvec;
            for (j=nodes[i].adjlist; j != -1; j = adjvec[j].next) {
                if (j != OUTSIDE) {
                    if (adjvec[j].weight > bestval) {
                        bestn0 = i;
                        bestn1 = j;
                        bestval = adjvec[j].weight;
                    } else if (-adjvec[j].weight > bestval) {
                        bestn0 = i;
                        bestn1 = j;
                        bestval = -adjvec[j].weight;
                    }
                }
            }
        }
    }
    if (bestval >= 0) {
        *n0 = bestn0;
        *n1 = bestn1;
        return 0;
    }

    i = OUTSIDE;
    adjvec = nodes[i].adjvec;
    for (j=nodes[i].adjlist; j != -1; j = adjvec[j].next) {
        if (adjvec[j].weight > bestval) {
            bestn0 = i;
            bestn1 = j;
            bestval = adjvec[j].weight;
        } else if (-adjvec[j].weight > bestval) {
            bestn0 = i;
            bestn1 = j;
            bestval = -adjvec[j].weight;
        }
    }

    if (bestval == -1) {
        return SELECT_TOUR;
    } else {
        *n0 = bestn0;
        *n1 = bestn1;
        return 0;
    }
}

#endif /* USE_DUMBO */

static void savetour (tspsearch *s)
{
    int i, j;

#ifdef DEBUG
    printf ("Solution found, saving tour, value %.0f\n", s->currentval);
#endif

    s->bestval = s->currentval;
    for (i=0; i<s->nnodes; i++) {
        for (j=0; j<s->nnodes; j++) {
            s->nodes[i].adjvec[j].bestvalue = s->nodes[i].adjvec[j].value;
        }
    }
}

static void unroll_stack (tspsearch *s, int smark)
{
    int i = s->histtop;

    while (i > smark) {
        i--;
        unset_edge (s, &s->searchhist[i]);
    }
    s->histtop = smark;
}

static double lowerbound (tspsearch *s)
{
    int i, j;
    node *nodes = s->nodes;
    int nnodes = s->nnodes;
    double lb = 0.0;
    int min1, min2;
    int cnt;

    for (i=0; i<nnodes; i++) {
        if (i == OUTSIDE) {
            min1 = INT_MAX;
            cnt = 0;
            for (j=nodes[i].adjlist; j != -1; j = nodes[i].adjvec[j].next) {
                if (nodes[i].adjvec[j].weight < 0) {
                    lb += nodes[i].adjvec[j].weight;
                    cnt++;
                    if (-nodes[i].adjvec[j].weight < min1) {
                        min1 = -nodes[i].adjvec[j].weight;
                    }
                } else if (nodes[i].adjvec[j].weight < min1) {
                    min1 = nodes[i].adjvec[j].weight;
                }
            }
            if ((nodes[i].deg + cnt) & 1) {
                lb += min1;
            }
        } else if (nodes[i].deg == 0) {
            min1 = INT_MAX;
            min2 = INT_MAX;
            for (j=nodes[i].adjlist; j != -1; j = nodes[i].adjvec[j].next) {
                if (nodes[i].adjvec[j].weight < min2) {
                    if (nodes[i].adjvec[j].weight < min1) {
                        if (nodes[i].adjvec[j].hi - nodes[i].adjvec[j].lo
                            <= 1) {
                            min2 = min1;
                            min1 = nodes[i].adjvec[j].weight;
                        } else {
                            min1 = nodes[i].adjvec[j].weight;
                            min2 = min1;
                        }
                    } else {
                        min2 = nodes[i].adjvec[j].weight;
                    }
                }
            }
            lb += min1 + min2;
        } else if (nodes[i].deg == 1) {
            min1 = INT_MAX;
            for (j=nodes[i].adjlist; j != -1; j = nodes[i].adjvec[j].next) {
                if (nodes[i].adjvec[j].weight < min1) {
                    min1 = nodes[i].adjvec[j].weight;
                }
            }
            lb += min1;
        }
    }
    return ceil(lb/2.0);
}

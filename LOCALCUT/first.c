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
/*  void CCchunk_init_separate_timer (CCchunk_separate_timer *timer)        */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCchunk_init_find_timer (CCchunk_find_timer *timer)                */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCchunk_init_lift_timer (CCchunk_lift_timer *timer)                */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCchunk_init_oracle_timer (CCchunk_oracle_timer *timer)            */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCchunk_init_localcut_timer (CCchunk_localcut_timer *timer)        */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCchunk_print_separate_timer (CCchunk_separate_timer *timer)       */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCchunk_print_find_timer (CCchunk_find_timer *timer)               */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCchunk_print_lift_timer (CCchunk_lift_timer *timer)               */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCchunk_print_oracle_timer (CCchunk_oracle_timer *timer)           */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCchunk_print_localcut_timer (CCchunk_localcut_timer *timer)       */
/*    MISSING                                                               */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "localcut.h"
#include "util.h"
#include "macrorus.h"

#undef  STUDY5          /* to print objective functions */
#undef  STUDY6          /* to print new signatures */
#undef  STUDY12         /* to print all integerizes */
#define STUDY13         /* to print integerizes with high errors */
#undef  STUDY14         /* to print lp solver & oracle timings */

#define MAXEFFORT 200   /* in oracle */
#undef  INIT_SIGNATURES /* build an initial set of signatures */
#undef  INIT_EXTRASIGS

/* only one of these three representations should be used */
#define REP_NONBASIS    /* represent using a nonbasis */
#undef  REP_EHALF       /* represent using E 1/2 */
#undef  REP_ENONZERO    /* represent using E 1/2 union E 1 */

#define MAXITERATIONS 200

#define CFRAC_LIMIT INT_MAX   /* avoid numbers > this in cfrac_approx */

#undef  DUMPCHUNKS

#define INTEGERIZE_MUL (16*9*5*7*11*13*17)

typedef struct adjinfo {
    int n;
    int e;
} adjinfo;

typedef struct node {
    int deg;
    adjinfo *adj;
    int equality;
#ifdef REP_NONBASIS
#define NM_UNSEEN 0
#define NM_NEROOT 1
#define NM_EVEN 2
#define NM_ODD 3
    int mark;
#endif
} node;

typedef struct edge {
    int end0;
    int end1;
    int fixed;
    int in_rep;
} edge;

typedef struct graph {
    int ncount;
    int ecount;
    int bcount;
    node *nodelist;
    edge *edgelist;
    adjinfo *adjspace;
    double *target;
    CCchunklp *lp;
    CCchunk_graph *chunk;
} graph;


static int
    separate (graph *g, CCchunk_graph *chunk, CCchunk_separate_timer *timer,
        CCchunk_fault_callback *callback),
    collect_solutions (graph *g, int *sig_list, int nsigs,
                              int *p_nsols, int **p_sols),
    graph_build (graph *g, CCchunk_graph *chunk),
    get_signature (graph *g, int *obj, int rhs, int *x, int *found,
            CCchunk_oracle_timer *timer);

static void
    build_ineq (int ecount, int *c, int rhs, CCchunk_ineq *f),
    graph_init (graph *g),
    graph_free (graph *g),
    build_representation (graph *g),
#ifdef REP_NONBASIS
    get_component (graph *g, int n, int *oddedge),
#endif
    integerize_check (int count, double *dvec, int *ivec),
    integerize_vector (int count, double *dvec, int *ivec),
    cfrac_approx (double d, int *p_num, int *p_den),
    integerize_vector2 (int count, double *dvec, int *ivec),
    report_objective (int ecount, int *obj, int rhs),
    report_new_signature (int ecount, int *x, int objslack);

static double
    get_cutslack (int ecount, int *cut, int rhs, double *target);


int CCchunk_separate (CCchunk_graph *chunk, CCchunk_separate_timer *timer,
        CCchunk_fault_callback *callback)
{
    int rval = 0;
    graph g;
#ifdef DUMPCHUNKS
    int i;
#endif

    CCutil_start_timer (&timer->all);
    
    graph_init (&g);

#ifdef DUMPCHUNKS
    printf ("CCchunk_separate\n");
    printf ("%d %d\n", chunk->ncount, chunk->ecount);
    for (i=0; i<chunk->ecount; i++) {
        printf ("%d %d %.16f %d\n",chunk->end0[i], chunk->end1[i], chunk->weight[i],
                chunk->fixed[i]);
    }
    for (i=0; i<chunk->ncount; i++) {
        if (i>0) printf (" ");
        printf ("%d",chunk->equality[i]);
    }
    printf ("\n");
#endif /* DUMPCHUNKS */

    rval = graph_build (&g, chunk);
    if (rval) {
        fprintf (stderr, "graph_build failed\n");
        goto  CLEANUP;
    }

    rval = separate (&g, chunk, timer, callback);
    if (rval) {
        fprintf (stderr, "separate failed\n");
        goto CLEANUP;
    }

  CLEANUP:
    CCutil_stop_timer (&timer->all, 0);
    graph_free (&g);
    return rval;
}

static int separate (graph *g, CCchunk_graph *chunk, CCchunk_separate_timer *timer,
        CCchunk_fault_callback *callback)
{
    int ecount = g->ecount;
    int bcount;
    edge *edgelist = g->edgelist;
    int rval = 0;
    int lpstat;
    double *c = (double *) NULL; /* alpha is stored in c[0] */
    int    *c_int = (int *) NULL;
    int    *c2_int = (int *) NULL;
    double *x = (double *) NULL;
    int    *x_int = (int *) NULL;
    int     maxsignatures = 0;
    int    *sig_list = (int *) NULL;
    int iter;
    int nsigs = 0;
    int i;
    int j;
    int finished = 0;
    int found = 0;
    CCchunk_fault fault;
#ifdef STUDY14
    double s;
    double szeit = CCutil_zeit();
    double szeit2;
    double oracle_zeit = 0.0;
    double lpsolver_zeit = 0.0;
#endif

    fault.a.coef = (int *) NULL;
    fault.nsols = 0;
    fault.sols = (int *) NULL;

    build_representation (g);

    bcount = g->bcount;

    if (bcount == 0) {
        /*
        fprintf (stderr, "Representation has 0 edges\n");
        */
        rval = 0; goto CLEANUP;
    }

    maxsignatures = MAXITERATIONS;
#ifdef INIT_SIGNATURES
    maxsignatures += bcount;
#ifdef INIT_EXTRASIGS
    maxsignatures += bcount;
#endif
#endif

    c            = CC_SAFE_MALLOC (bcount+1, double);
    c_int        = CC_SAFE_MALLOC (bcount+1, int);
    c2_int       = CC_SAFE_MALLOC (ecount+1, int);
    x            = CC_SAFE_MALLOC (bcount, double);
    x_int        = CC_SAFE_MALLOC (ecount, int);
    fault.a.coef = CC_SAFE_MALLOC (ecount, int);
    sig_list     = CC_SAFE_MALLOC (ecount * maxsignatures, int);

    if (c              == (double *) NULL ||
        c_int          == (int *) NULL ||
        c2_int         == (int *) NULL ||
        x              == (double *) NULL ||
        x_int          == (int *) NULL ||
        fault.a.coef   == (int *) NULL ||
        sig_list       == (int *) NULL) {
        fprintf (stderr, "Out of memory in separate\n");
        rval = -1;
        goto CLEANUP;
    }


    for (i=0, j=0; i<ecount; i++) {
        if (edgelist[i].in_rep) {
            x[j++] = g->target[i];
        }
    }
    
#ifdef STUDY14
    szeit2 = CCutil_zeit();
#endif

    CCutil_start_timer (&timer->lpsolver);
    rval = CCchunk_lpinit (&g->lp, "firstlp", bcount, x);
    CCutil_stop_timer (&timer->lpsolver, 0);
    if (rval) {
        fprintf (stderr, "CCchunk_lpinit failed\n");
        goto CLEANUP;
    }
#ifdef STUDY14
    lpsolver_zeit += CCutil_zeit() - szeit2;
#endif

    for (iter = 0; iter < MAXITERATIONS; iter++) {
#ifdef STUDY14
        szeit2 = CCutil_zeit();
#endif
        CCutil_start_timer (&timer->lpsolver);
        rval = CCchunk_lpsolve (g->lp, &lpstat, c, &c[bcount]);
        CCutil_stop_timer (&timer->lpsolver, 0);
        if (rval) {
            fprintf (stderr, "CCchunk_lpsolve failed\n");
            goto CLEANUP;
        }
#ifdef STUDY14
        lpsolver_zeit += CCutil_zeit() - szeit2;
#endif

        if (lpstat == CC_CHUNK_LPFEASIBLE) {
            /* not faulty */
            rval = 0;
#ifdef STUDY14
            printf ("faultless in %d iterations time %.2f lpsolver %.2f oracle %.2f misc\n",
                    iter, lpsolver_zeit, oracle_zeit,
                    CCutil_zeit() - szeit - lpsolver_zeit - oracle_zeit);
#endif
            goto CLEANUP;
        } else if (lpstat != CC_CHUNK_LPINFEASIBLE) {
            fprintf (stderr, "CCchunk_lpsolve lpstat %d\n", lpstat);
            rval = -1;
            goto CLEANUP;
        }

        integerize_vector (bcount+1, c, c_int);

#ifdef STUDY14
        szeit2 = CCutil_zeit();
#endif

        for (i=0, j=0; i<ecount; i++) {
            if (edgelist[i].in_rep) {
                c2_int[i] = c_int[j++];
            } else {
                c2_int[i] = 0;
            }
        }
        c2_int[ecount] = c_int[bcount];
        
        rval = get_signature (g, c2_int, c2_int[ecount], x_int, &found,
                              &timer->oracle);
        if (rval) {
            fprintf (stderr, "get_signature failed\n");
            goto CLEANUP;
        }
        if (found) {
#ifdef STUDY14
            oracle_zeit += CCutil_zeit() - szeit2;
#endif
            for (i=0, j=0; i<ecount; i++) {
                if (edgelist[i].in_rep) {
                    x[j++] = (double) x_int[i];
                }
            }
            rval = CCchunk_lpaddcol (g->lp, x);
            if (rval) {
                fprintf (stderr, "CCchunk_lpaddcol failed\n");
                goto CLEANUP;
            }
            if (nsigs >= maxsignatures) {
                fprintf (stderr, "TOO MANY SIGNATURES\n");
            } else {
                for (i=0; i<ecount; i++) sig_list[nsigs*ecount+i] = x_int[i];
                nsigs++;
            }
        } else {
#ifdef STUDY14
            oracle_zeit += CCutil_zeit() - szeit2;
            s = get_cutslack (ecount, c2_int, c2_int[ecount], g->target);
            printf ("   faulty (%.2f) in %d iterations time %.2f lpsolver %.2f oracle %.2f misc\n",
                    s, iter, lpsolver_zeit, oracle_zeit,
                    CCutil_zeit() - szeit - lpsolver_zeit - oracle_zeit);
#endif
            build_ineq (ecount, c2_int, c2_int[ecount], &fault.a);
            rval = collect_solutions (g, sig_list, nsigs, &fault.nsols,
                                      &fault.sols);
            if (rval) {
                fprintf (stderr, "collect_solutions failed\n");
                goto CLEANUP;
            }

            CCutil_suspend_timer (&timer->all);
            rval = (*callback->func) (chunk, &fault, &finished,
                                      callback->u_data);
            CCutil_resume_timer (&timer->all);
            if (rval) {
                fprintf (stderr, "fault callback failed\n");
                goto CLEANUP;
            }

            CC_IFFREE (fault.sols, int);

            if (finished) {
                rval = 0;
                goto CLEANUP;
            }

            for (i=0; i<bcount; i++) {
                if (c_int[i] != 0) {
                    break;
                }
            }
            if (i >= bcount) {
                /* all cuts found */
                rval = 0;
                goto CLEANUP;
            }

#if 0
            printf ("relaxing row %d\n", i);
            fflush (stdout);
#endif

            rval = CCchunk_lprelaxrow (g->lp, i);
            if (rval) {
                fprintf (stderr, "CCchunk_lprelaxrow failed\n");
                goto CLEANUP;
            }
        }
    }

    printf ("ITERATION LIMIT %d exceeded in separate\n",
            MAXITERATIONS);

  CLEANUP:
    CC_IFFREE (fault.a.coef, int);
    CC_IFFREE (fault.sols, int);
    CC_IFFREE (x_int, int);
    CC_IFFREE (x, double);
    CC_IFFREE (c_int, int);
    CC_IFFREE (c2_int, int);
    CC_IFFREE (c, double);
    CC_IFFREE (sig_list, int);

    return rval;
}

static void build_ineq (int ecount, int *c, int rhs, CCchunk_ineq *f)
{
    int i;

    for (i=0; i<ecount; i++) f->coef[i] = c[i];
    f->rhs = rhs;
}

static int collect_solutions (graph *g, int *sig_list, int nsigs,
                               int *p_nsols, int **p_sols)
{
    int rval;
    int *bas = (int *) NULL;
    int *sols = (int *) NULL;
    int nsols;
    int ecount = g->ecount;
    int i;
    int j;

    *p_sols = (int *) NULL;
    *p_nsols = 0;

    if (nsigs == 0) return 0;

    bas = CC_SAFE_MALLOC (nsigs, int);
    if (bas == (int *) NULL) {
        fprintf (stderr, "Out of memory in collect_solutions\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCchunk_lpbasis (g->lp, nsigs, bas);
    if (rval) {
        fprintf (stderr, "CCchunk_lpbasis failed\n");
        goto CLEANUP;
    }

    for (i=0, nsols=0; i<nsigs; i++) {
        if (bas[i] == 1) {
            nsols++;
        }
    }

    if (nsols == 0) {
        *p_nsols = 0;
        *p_sols = (int *) NULL;
        rval = 0;
        goto CLEANUP;
    }

    sols = CC_SAFE_MALLOC (ecount * nsols, int);
    if (sols == (int *) NULL) {
        fprintf (stderr, "Out of memory in collect_solutions\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0, nsols=0; i<nsigs; i++) {
        if (bas[i] == 1) {
            for (j=0; j<ecount; j++) {
                sols[ecount * nsols + j] = sig_list[ecount * i + j];
            }
            nsols++;
        }
    }

    *p_nsols = nsols;
    *p_sols = sols;

    rval = 0;

  CLEANUP:
    CC_IFFREE (bas, int);
    if (rval) {
        CC_IFFREE (sols, int);
    }
    return rval;
}

static void graph_init (graph *g)
{
    g->ncount = 0;
    g->ecount = 0;
    g->bcount = 0;
    g->nodelist      = (node *) NULL;
    g->edgelist      = (edge *) NULL;
    g->adjspace      = (adjinfo *) NULL;
    g->target        = (double *) NULL;
    g->lp            = (CCchunklp *) NULL;
}

static int graph_build (graph *g, CCchunk_graph *chunk)
{
    int rval = 0;
    int end0, end1;
    adjinfo *p;
    int i, j;

    g->nodelist      = (node *) NULL;
    g->edgelist      = (edge *) NULL;
    g->adjspace      = (adjinfo *) NULL;
    g->target        = (double *) NULL;
    g->chunk         = chunk;

    g->ncount = chunk->ncount;
    g->ecount = chunk->ecount;
    g->bcount = 0;

    g->nodelist = CC_SAFE_MALLOC (g->ncount, node);
    g->edgelist = CC_SAFE_MALLOC (g->ecount, edge);
    g->adjspace = CC_SAFE_MALLOC (g->ecount*2, adjinfo);
    g->target   = CC_SAFE_MALLOC (g->ecount, double);
    if (g->nodelist == (node *) NULL ||
        g->edgelist == (edge *) NULL ||
        g->adjspace == (adjinfo *) NULL ||
        g->target   == (double *) NULL) {
        fprintf (stderr, "Out of memory in graph_build\n");
        rval = -1;
        goto CLEANUP;
    }

    for (i=0; i<chunk->ncount; i++) {
        g->nodelist[i].equality = chunk->equality[i];
        g->nodelist[i].deg = 0;
        g->nodelist[i].adj = (adjinfo *) NULL;
    }

    for (i=0; i<chunk->ecount; i++) {
        g->edgelist[i].end0 = chunk->end0[i];
        g->edgelist[i].end1 = chunk->end1[i];
        g->edgelist[i].fixed = chunk->fixed[i];
        g->target[i] = chunk->weight[i];
        g->nodelist[chunk->end0[i]].deg++;
        g->nodelist[chunk->end1[i]].deg++;
    }

    for (i=0, p=g->adjspace; i<chunk->ncount; i++) {
        g->nodelist[i].adj = p;
        p += g->nodelist[i].deg;
        g->nodelist[i].deg = 0;
    }

    for (i=0; i<chunk->ecount; i++) {
        end0 = g->edgelist[i].end0;
        end1 = g->edgelist[i].end1;
        j = g->nodelist[end0].deg++;
        g->nodelist[end0].adj[j].n = end1;
        g->nodelist[end0].adj[j].e = i;
        j = g->nodelist[end1].deg++;
        g->nodelist[end1].adj[j].n = end0;
        g->nodelist[end1].adj[j].e = i;
    }

    return 0;

  CLEANUP:
    CC_IFFREE (g->nodelist, node);
    CC_IFFREE (g->edgelist, edge);
    CC_IFFREE (g->adjspace, adjinfo);
    CC_IFFREE (g->target,   double);
    g->ncount = 0;
    g->ecount = 0;
    g->bcount = 0;
    return rval;
}

static void graph_free (graph *g)
{
    CC_IFFREE (g->nodelist, node);
    CC_IFFREE (g->edgelist, edge);
    CC_IFFREE (g->adjspace, adjinfo);
    CC_IFFREE (g->target, double);
    if (g->lp) {
        CCchunk_lpfree (&g->lp);
    }
    g->ncount = 0;
    g->ecount = 0;
    g->bcount = 0;
}

static double get_cutslack (int ecount, int *cut, int rhs, double *target)
{
    double sum = 0.0;
    int i;

    for (i=0; i<ecount; i++) {
        sum += cut[i] * target[i];
    }
    return ((double) rhs) - sum;
}

static void build_representation (graph *g)
{
#ifdef REP_ENONZERO
    int i;
    int ecount = g->ecount;
    edge *edgelist = g->edgelist;

    for (i=0; i<ecount; i++) {
        edgelist[i].in_rep = 1;
    }
    g->bcount = ecount;
#endif /* REP_ENONZERO */
#ifdef REP_EHALF
    int i;
    int j;
    int ecount = g->ecount;
    int bcount = 0;
    edge *edgelist = g->edgelist;

    for (i=0; i<ecount; i++) {
        if (edgelist[i].fixed == -1) {
            edgelist[i].in_rep = 1;
            bcount++;
        } else {
            edgelist[i].in_rep = 0;
        }
    }
    g->bcount = bcount;
#endif /* REP_EHALF */
#ifdef REP_NONBASIS
    int i;
    int ncount = g->ncount;
    int ecount = g->ecount;
    int bcount;
    node *nodelist = g->nodelist;
    edge *edgelist = g->edgelist;
    int oddedge;

    for (i=0; i<ncount; i++) {
        if (nodelist[i].equality) {
            nodelist[i].mark = NM_UNSEEN;
        } else {
            nodelist[i].mark = NM_NEROOT;
        }
    }

    for (i=0; i<ecount; i++) {
        edgelist[i].in_rep = -1;
    }
    
    for (i=0; i<ncount; i++) {
        if (nodelist[i].mark == NM_NEROOT) {
            oddedge = 0;
            nodelist[i].mark = NM_ODD;
            get_component (g, i, &oddedge);
        }
    }
    for (i=0; i<ncount; i++) {
        if (nodelist[i].mark == NM_UNSEEN) {
            oddedge = 1;
            nodelist[i].mark = NM_ODD;
            get_component (g, i, &oddedge);
        }
    }

    bcount = 0;
    for (i=0; i<ecount; i++) {
        if (edgelist[i].in_rep) {
            bcount++;
        }
    }
    g->bcount = bcount;
#endif /* REP_NONBASIS */
}

#ifdef REP_NONBASIS
static void get_component (graph *g, int n, int *oddedge)
{
    node *nodelist = g->nodelist;
    edge *edgelist = g->edgelist;
    int newn;
    int newe;
    int i;

    for (i=0; i<nodelist[n].deg; i++) {
        newe = nodelist[n].adj[i].e;
        newn = nodelist[n].adj[i].n;
        if (edgelist[newe].fixed == -1) {
            if (nodelist[newn].mark == NM_UNSEEN) {
                nodelist[newn].mark = (NM_ODD + NM_EVEN) - nodelist[n].mark;
                edgelist[newe].in_rep = 0;
                get_component (g, newn, oddedge);
            } else if (*oddedge && nodelist[n].mark == nodelist[newn].mark) {
                edgelist[newe].in_rep = 0;
                (*oddedge)--;
            } else if (edgelist[newe].in_rep == -1) {
                edgelist[newe].in_rep = 1;
            }
        } else {
            edgelist[newe].in_rep = 0;
        }
    }
}
#endif /* REP_NONBASIS */

static void integerize_check (int count, double *dvec, int *ivec)
{
    double error;
    double scale;
    int i;

    scale = 1.0;
    for (i=0; i<count; i++) {
        if (ivec[i]) {
            scale = CC_OURABS(dvec[i]/ivec[i]);
            break;
        }
    }
    error = 0.0;
    for (i=0; i<count; i++) {
        error += CC_OURABS(dvec[i]/scale - ivec[i]);
    }
    if (error > 0.001) {
        printf ("bad integerize error %f\n", error);
        printf ("from:");
        for (i=0; i<count; i++) printf (" %f", dvec[i]);
        printf ("\n");
        printf ("to:");
        for (i=0; i<count; i++) printf (" %d", ivec[i]);
        printf ("\n");
        fflush (stdout);
    }
}

static void cfrac_approx (double d, int *p_num, int *p_den)
{
    int sgn = 1;
    int a;
    int g2 = 0;
    int g1 = 1;
    int g = 0;
    int h2 = 1;
    int h1 = 0;
    int h = 1;
    int bestg = 0;
    int besth = 1;
    double besterr = 1.0;
    int lim;

    if (d < 0.0) {
        d = -d;
        sgn = -1;
    }

    for (;;) {
        a = (int) d;
        if (a) {
            lim = CFRAC_LIMIT/a;
        } else {
            lim = CFRAC_LIMIT;
        }
        if (g1 > lim || h1 > lim ||
            g2 > CFRAC_LIMIT - g1*a || h2 > CFRAC_LIMIT - h1*a)
            break;
        g = a*g1 + g2;
        h = a*h1 + h2;
        g2 = g1;
        h2 = h1;
        g1 = g;
        h1 = h;
        if (d - (double) a < besterr) {
            bestg = g;
            besth = h;
            besterr = d - (double) a;
            if (besterr < 0.0001) break;
        }
        d = 1 / (d - a);
    }
    *p_num = bestg * sgn;
    *p_den = besth;
}

static void integerize_vector (int count, double *dvec, int *ivec)
{
    int i;
    int num;
    int den;
    int g;
    int scale = 1;
    double dmax;

#ifdef STUDY12
    printf ("integerize:");
    for (i=0; i<count; i++) printf (" %f", dvec[i]);
    printf ("\n");
#endif

    dmax = 0.0;
    for (i=0; i<count; i++) {
        if (CC_OURABS(dvec[i]) > dmax) dmax = CC_OURABS(dvec[i]);
    }
    dmax = 1.0/dmax;

    for (i=0; i<count; i++) {
        cfrac_approx (dvec[i]*dmax, &num, &den);
        g = CCutil_our_gcd (den, scale);
        scale /= g;
        if (INT_MAX/CC_OURABS(scale) < CC_OURABS(den)) goto FAILURE;
        scale *= den;
    }

    for (i=0; i<count; i++) {
        cfrac_approx (dvec[i]*dmax, &num, &den);
        den = scale / den;
        if (INTEGERIZE_MUL/CC_OURABS(den) < CC_OURABS(num)) goto FAILURE;
        ivec[i] = num * den;
    }

#ifdef STUDY12
    printf ("integerized");
    for (i=0; i<count; i++) printf (" %d", ivec[i]);
    printf ("\n");
#endif

#ifdef STUDY13
    integerize_check (count, dvec, ivec);
#endif

    return;

  FAILURE:
    printf ("Overflow in integerize_vector.  Calling integerize_vector2\n");
    integerize_vector2 (count, dvec, ivec);
}

static void integerize_vector2 (int count, double *dvec, int *ivec)
{
    int i;
    double scale;
    int divv;

#ifdef STUDY12
    printf ("integerize(2):");
    for (i=0; i<count; i++) printf (" %f", dvec[i]);
    printf ("\n");
#endif

    scale = CC_OURABS(dvec[0]);
    for (i=1; i<count; i++) {
        if (CC_OURABS(dvec[i]) > scale) scale = CC_OURABS(dvec[i]);
    }
    if (scale == 0.0) {
        for (i=0; i<count; i++) {
            ivec[i] = 0;
        }
        return;
    }

    for (i=0; i<count; i++) {
        if (dvec[i] < 0.0) {
            ivec[i] = -((int) (-dvec[i] * INTEGERIZE_MUL / scale + 0.5));
        } else {
            ivec[i] = ((int) (dvec[i] * INTEGERIZE_MUL / scale + 0.5));
        }
    }
    divv = ivec[0];
    for (i=1; i<count; i++) {
        divv = CCutil_our_gcd (divv, ivec[i]);
    }

    if (divv > 1) {
        for (i=0; i<count; i++) {
            ivec[i] /= divv;
        }
    }

#ifdef STUDY12
    printf ("integerized");
    for (i=0; i<count; i++) printf (" %d", ivec[i]);
    printf ("\n");
#endif

#ifdef STUDY13
    integerize_check (count, dvec, ivec);
#endif
}

static int get_signature (graph *g, int *obj, int rhs, int *x, int *found,
    CCchunk_oracle_timer *timer)
{
    CCchunk_graph *ch = g->chunk;
    int ecount = ch->ecount;
    CCchunk_ineq a;
    int objval = 0;
    int rval;

    a.coef = (int *) NULL;
    *found = 0;

#ifdef STUDY5
    report_objective (ecount, obj, rhs);
#endif

    a.coef = CC_SAFE_MALLOC (ecount, int);
    if (a.coef == (int *) NULL) {
        fprintf (stderr, "Out of memory in get_signature\n");
        rval = 1; goto CLEANUP;
    }

    build_ineq (ecount, obj, rhs, &a);

    rval = CCchunk_oracle (ch, &a, x, &objval, 1, MAXEFFORT, timer);
    if (rval == CC_CHUNK_ORACLE_INFEASIBLE) {
        *found = 0;
        rval = 0; goto CLEANUP;
    } else if (rval == CC_CHUNK_ORACLE_SEARCHLIMITEXCEEDED) {
        fprintf (stderr, "CCchunk_oracle node limit exceeded\n");
        goto CLEANUP;
    } else if (rval) {
        fprintf (stderr, "CCchunk_oracle failed\n");
        goto CLEANUP;
    }

#ifdef STUDY6
    report_new_signature (ecount, x, objval - rhs);
#endif
    *found = 1;
    rval = 0;

 CLEANUP:
    CC_IFFREE (a.coef, int);
    return rval;
}

static void report_objective (int ecount, int *obj, int rhs)
{
    int i;

    printf ("objective ");
    for (i=0; i < ecount; i++) {
        printf (" %d", obj[i]);
    }
    printf ("> %d\n", rhs);
}

static void report_new_signature (int ecount, int *x, int objslack)
{
    int i;

    printf ("got new signature (by %d)", objslack);
    for (i=0; i < ecount; i++) {
        printf (" %d", x[i]);
    }
    printf ("\n");
}


CCchunk_graph *CCchunk_graph_alloc (int ncount, int ecount)
{
    CCchunk_graph *c;
    int i;

    c = CC_SAFE_MALLOC (1, CCchunk_graph);
    if (!c) {
        printf ("out of memory in CCchunk_graph_alloc\n");
        return (CCchunk_graph *) NULL;
    }

    c->end0 = CC_SAFE_MALLOC (ecount, int);
    c->end1 = CC_SAFE_MALLOC (ecount, int);
    c->fixed = CC_SAFE_MALLOC (ecount, int);
    c->weight = CC_SAFE_MALLOC (ecount, double);
    c->equality = CC_SAFE_MALLOC (ncount, int);
    c->members = CC_SAFE_MALLOC (ncount, int *);

    if (!c->end0 || !c->end1 || !c->fixed || !c->weight ||
        !c->equality || !c->members) {
        printf ("out of memory in CCchunk_graph_alloc\n");
        return (CCchunk_graph *) NULL;
    }

    for (i = 0; i < ncount; i++)
        c->members[i] = (int *) NULL;

    c->ncount = ncount;
    c->ecount = ecount;

    return c;
}

void CCchunk_graph_free (CCchunk_graph *c)
{
    int i;

    if (c) {
        CC_IFFREE (c->end0, int);
        CC_IFFREE (c->end1, int);
        CC_IFFREE (c->fixed, int);
        CC_IFFREE (c->weight, double);
        CC_IFFREE (c->equality, int);
        if (c->members) {
            for (i = 0; i < c->ncount; i++) {
                CC_IFFREE (c->members[i], int);
            }
            CC_FREE (c->members, int *);
        }
        CC_FREE (c, CCchunk_graph);
    }
}

void CCchunk_init_separate_timer (CCchunk_separate_timer *timer)
{
    CCutil_init_timer (&timer->all,      "separate");
    CCutil_init_timer (&timer->lpsolver, "  separate lpsolver");
    CCchunk_init_oracle_timer (&timer->oracle);
}

void CCchunk_init_find_timer (CCchunk_find_timer *timer)
{
    CCutil_init_timer (&timer->all,      "find chunks");
    CCutil_init_timer (&timer->shrink,   "  find shrink");
    CCutil_init_timer (&timer->locate,   "  find locate");
}

void CCchunk_init_lift_timer (CCchunk_lift_timer *timer)
{
    CCutil_init_timer (&timer->all,                 "lift");
    CCutil_init_timer (&timer->liberate_equality,   "  lift liberate equality");
    CCutil_init_timer (&timer->liberate_fixed,      "  lift liberate fixed");
    CCutil_init_timer (&timer->strengthen_edges,    "  lift strengthen edges");
    CCutil_init_timer (&timer->strengthen_equality, "  lift strengthen equality");
    CCutil_init_timer (&timer->strengthen_work,     "  lift strengthen work");
    CCutil_init_timer (&timer->decompose,           "  lift decompose");
    CCutil_init_timer (&timer->tilt_oracle,         "  lift tilt oracle");
    CCutil_init_timer (&timer->liberate_oracle,     "  lift liberate oracle");
    CCutil_init_timer (&timer->verify_oracle,       "  lift verify oracle");

    CCchunk_init_oracle_timer (&timer->oracle);
}

void CCchunk_init_oracle_timer (CCchunk_oracle_timer *timer)
{
    CCutil_init_timer (&timer->all,                 "oracle");
    CCutil_init_timer (&timer->bnbtsp,              "  oracle bnbtsp");
    CCutil_init_timer (&timer->tinytsp,             "  oracle tinytsp");
}

void CCchunk_init_localcut_timer (CCchunk_localcut_timer *timer)
{
    CCchunk_init_find_timer (&timer->find);
    CCchunk_init_separate_timer (&timer->separate);
    CCchunk_init_lift_timer (&timer->lift);
    CCutil_init_timer (&timer->all,                 "localcuts");
}

void CCchunk_print_separate_timer (CCchunk_separate_timer *timer)
{
    CCutil_total_timer (&timer->all,                 3);
    CCutil_total_timer (&timer->lpsolver,            3);
    CCchunk_print_oracle_timer (&timer->oracle);
}

void CCchunk_print_find_timer (CCchunk_find_timer *timer)
{
    CCutil_total_timer (&timer->all,                 3);
    CCutil_total_timer (&timer->shrink,              3);
    CCutil_total_timer (&timer->locate,              3);
}

void CCchunk_print_lift_timer (CCchunk_lift_timer *timer)
{
    CCutil_total_timer (&timer->all,                 3);
    CCutil_total_timer (&timer->liberate_equality,   3);
    CCutil_total_timer (&timer->liberate_fixed,      3);
    CCutil_total_timer (&timer->strengthen_edges,    3);
    CCutil_total_timer (&timer->strengthen_equality, 3);
    CCutil_total_timer (&timer->strengthen_work,     3);
    CCutil_total_timer (&timer->decompose,           3);
    CCutil_total_timer (&timer->tilt_oracle,         3);
    CCutil_total_timer (&timer->liberate_oracle,     3);
    CCutil_total_timer (&timer->verify_oracle,       3);

    CCchunk_print_oracle_timer (&timer->oracle);
}

void CCchunk_print_oracle_timer (CCchunk_oracle_timer *timer)
{
    CCutil_total_timer (&timer->all,               3);
    CCutil_total_timer (&timer->bnbtsp,            3);
    CCutil_total_timer (&timer->tinytsp,           3);
}

void CCchunk_print_localcut_timer (CCchunk_localcut_timer *timer)
{
    CCchunk_print_find_timer (&timer->find);
    CCchunk_print_separate_timer (&timer->separate);
    CCchunk_print_lift_timer (&timer->lift);
    CCutil_total_timer (&timer->all, 3);
}

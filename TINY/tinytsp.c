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
/*               LP-based Branch-and-Bound for Tiny TSPs                    */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: April 1, 1997                                                     */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCtiny_bnc_tsp (int ncount, CCdatagroup *dat, double *upbound,      */
/*      double *optval, int nodelimit)                                      */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtiny_bnc_msp (int ncount, int ecount, int *elist, int *elen,      */
/*      int depot, int *lower, int *upper, double *upperbound,              */
/*      int objsense, double *optval, int *xsol, int checkresult,           */
/*      int searchlimit)                                                    */
/*    MISSING                                                               */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tinytsp.h"
#include "lp.h"
#include "bigguy.h"
#include "cut.h"
#include "tsp.h"

#define CC_TINYTSP_MAXDOUBLE   (1e15)
#define CC_TINYTSP_EPSILON     (1e-10)
#define CC_TINYTSP_INTTOL      (1e-10)
#define CC_TINYTSP_SEARCHLIMIT (1 << 20)

#undef  TINYDEBUG
#undef  TINYNOISY
#undef  TIMINGS

#undef  TINY_FIRST_BRANCH
#define TINY_BEST_BRANCH

typedef struct tinyadj {
    int              edge;
    int              to;
} tinyadj;

typedef struct tinynode {
    tinyadj         *adj;
    int              deg;
    int              mark;
    int              magic;
    int              parity;          /* blossom algorithm */
} tinynode;

typedef struct tinyedge {
    struct tinyedge *next;
    int              ends[2];
    int              len;
    int              coef;
    CCbigguy         rc;
} tinyedge;

typedef struct tinygraph {
    tinynode        *nodelist;
    tinyadj         *adjspace;
    tinyedge        *edgelist;
    int              ncount;
    int              ecount;
    int              magiclabel;
} tinygraph;

typedef struct tinytooth {
    int              ends[2];
} tinytooth;

typedef struct tinycut {
    int             *nodes;
    tinytooth       *teeth;
    int              count;
    int              tcount;
    struct tinycut *next;
} tinycut;

typedef struct tinycomp {
    tinycut         *complist;
    int             *grab;
    int             *stack;
    int              ncomp;
} tinycomp;

typedef struct tiny_lp {
    tinygraph        graph;
    tinycomp         comp;
    CClp            *lp;
    double          *upper;
    double          *lower;
    double          *x;
    double          *xsol;
    double           val;
    double           upperbound;
    double           lowerbound;
    tinycut         *cuts;
    int              cutsize;
    int              ncuts;
    int              nbbnodes;
    int              searchlimit;
    int              status;
    int              foundtour;
    CCbigguy        *node_pi;
    CCbigguy        *cut_pi;
    int              cutpisize;
    double          *pi_double;
    int              pisize_double;
    int              checkresult;
    int              depot;
    int              objsense;
} tiny_lp;


static int
    permute_edges (int ecount, int *elist, int *elen,
        int *lower, int *upper, int *newecount, int **newelist,
        int **newelen, int **newlower, int **newupper, int **perm),
    init_tinylp (tiny_lp *lp, int ncount, int ecount, int *elist, int *elen,
        int *lower, int *upper, double *upperbound, int checkresult,
        int depot, int searchlimit, int objsense),
    load_tinylp (tiny_lp *lp),
    add_tinycut_list (tiny_lp *lp, tinycut **list, int *nadded),
    add_tinycut (tiny_lp *lp, tinycut *c, int *added),
    tiny_setbounds (tiny_lp *lp, int col, int lower, int upper),
    optimize_tinylp (tiny_lp *lp),
    tiny_pi (tiny_lp *lp),
    build_graph (tinygraph *g, int ncount, int ecount, int *elist, int *elen,
        int objsense),
    init_tinycomp (tinycomp *t, int ncount),
    tiny_connectcut (tiny_lp *lp, tinycomp *t, int *nadded),
    tiny_connect (tinygraph *g, tinycomp *tc, double *x, int depot),
    tiny_blossom (tiny_lp *lp, int *nadded),
    blossom_grab (tinygraph *g, int depot, double *x, int *handle,
        int count, tinycut **list, int *foundcut),
    tiny_brancher (tiny_lp *lp, int depth),
    tiny_find_branch (tinygraph *g, int depot, double *x),
    tiny_checkbound (tiny_lp *lp, int *cutoff),
    tiny_checktour (tiny_lp *lp, double *tourval),
    tiny_checkdual (tiny_lp *lp, int *cutoff),
    exact_subtours (tiny_lp *lp),
    add_exact (double val, int count, int *cutarray, void *pass_param),
    exact_blossoms (tiny_lp *lp),
    grab_nonzero_x (tiny_lp *lp, int *ecount, int **elist, double **x,
        double tol),
    lpcut_in_to_tinycut (CCtsp_lpcut_in *in, tinycut **out);

static void
    tiny_loadxsol (tiny_lp *lp, int newecount, int ecount, int *xsol,
        int *perm),
    tiny_loadbounds (tiny_lp *lp, int *lower, int *upper),
    init_tiny_lp_struct (tiny_lp *lp),
    free_tinylp (tiny_lp *lp),
    free_tinygraph (tinygraph *g),
    free_tinycomp (tinycomp *t),
    connect_search (tinygraph *g, int n, int marker, int *dstack, double *x,
        int *comp, int *count),
    blossom_search (tinygraph *g, int n, int marker, int *dstack,
        double *x, int *comp, int *count, int *parity),
#ifdef TINYDEBUG
    print_tinycut (tinycut *c),
#endif
    free_tinycut (tinycut *c),
    tiny_price (tiny_lp *lp, int phase1);


int CCtiny_bnc_tsp (int ncount, CCdatagroup *dat, double *upbound,
        double *optval, int nodelimit)
{
    int rval;
    int *lower = (int *) NULL;
    int *upper = (int *) NULL;
    int *elist = (int *) NULL;
    int *elen = (int *) NULL;
    int edgecount = ncount * (ncount-1) / 2;
    int i;
    int j;
    int k;

    lower = CC_SAFE_MALLOC (edgecount, int);
    upper = CC_SAFE_MALLOC (edgecount, int);
    elist = CC_SAFE_MALLOC (edgecount*2, int);
    elen  = CC_SAFE_MALLOC (edgecount, int);
    if (lower == (int *) NULL ||
        upper == (int *) NULL ||
        elist == (int *) NULL ||
        elen  == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCtiny_bnc_tsp\n");
        rval = -1;
        goto CLEANUP;
    }
    for (i=0, k=0; i<ncount; i++) {
        for (j=0; j<i; j++) {
            lower[k] = 0;
            upper[k] = 1;
            elist[2*k] = i;
            elist[2*k+1] = j;
            elen[k] = CCutil_dat_edgelen (i, j, dat);
            k++;
        }
    }

    rval = CCtiny_bnc_msp (ncount, edgecount, elist, elen, -1,
                    lower, upper, upbound, CC_TINYTSP_MINIMIZE,
                    optval, (int *) NULL, 1, nodelimit);
  CLEANUP:
    CC_IFFREE (lower, int);
    CC_IFFREE (upper, int);
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    return rval;
}

int CCtiny_bnc_msp (int ncount, int ecount, int *elist, int *elen, int depot,
        int *lower, int *upper, double *upperbound, int objsense,
        double *optval, int *xsol, int checkresult, int searchlimit)
{
    int rval = 0;
    int cutoff;
    tiny_lp lp;
#ifdef TIMINGS
    double sz = CCutil_zeit();
#endif
    int *perm     = (int *) NULL;
    int *newelist = (int *) NULL;
    int *newelen  = (int *) NULL;
    int *newupper = (int *) NULL;
    int *newlower = (int *) NULL;
    int newecount;

    /* printf ("t"); fflush (stdout); */

    if (optval) {
        *optval = 0.0;
    }
    init_tiny_lp_struct (&lp);

    rval = permute_edges (ecount, elist, elen, lower, upper,
             &newecount, &newelist, &newelen, &newlower, &newupper, &perm); 
    if (rval) {
        fprintf (stderr, "permute_edges failed\n");
        rval = CC_TINYTSP_ERROR;  goto CLEANUP;
    }

    rval = init_tinylp (&lp, ncount, newecount, newelist, newelen, newlower,
                        newupper, upperbound, checkresult, depot, searchlimit,
                        objsense);
    if (rval) {
        fprintf (stderr, "init_tinylp failed\n");
        rval = CC_TINYTSP_ERROR;  goto CLEANUP;
    }

    rval = optimize_tinylp (&lp);
    if (rval && rval != 2) {
        fprintf (stderr, "optimize_tinylp failed\n"); goto CLEANUP;
    }
#ifdef TINYNOISY
    printf ("Initial LP Value: %.4f\n", lp.val); fflush (stdout);
#endif

    rval = tiny_checkbound (&lp, &cutoff);
    if (rval) {
        fprintf (stderr, "tiny_checkbound failed\n");
        rval = CC_TINYTSP_ERROR; goto CLEANUP;
    }
    if (!cutoff) {
        rval = tiny_brancher (&lp, 0);
        if (rval) {
            fprintf (stderr, "tiny_brancher failed\n");
            rval = CC_TINYTSP_ERROR; goto CLEANUP;
        }
        if (lp.status == CC_TINYTSP_SEARCHLIMITEXCEEDED) {
#ifdef TIMINGS
            printf ("Hit Search Limit Nodes: %d   Cuts: %d   Time %.3f\n",
                    lp.nbbnodes, lp.ncuts, CCutil_zeit() - sz);
            fflush (stdout);
#endif
            rval = CC_TINYTSP_SEARCHLIMITEXCEEDED; goto CLEANUP;
        }
    }

#ifdef TIMINGS
    printf ("Number of Search Nodes: %d   Cuts: %d   Time %.3f\n",
            lp.nbbnodes, lp.ncuts, CCutil_zeit() - sz); fflush (stdout);
#endif

    if (!lp.foundtour) {
#ifdef TINYNOISY
        printf ("No tour found\n"); fflush (stdout);
#endif
        rval = CC_TINYTSP_INFEASIBLE;
    } else {
        if (lp.objsense == CC_TINYTSP_MAXIMIZE) {
            lp.upperbound *= -1;
        }
#ifdef TINYNOISY
        printf ("Optimal Value: %.4f\n", lp.upperbound); fflush (stdout);
#endif
        if (optval) {
            *optval = lp.upperbound;
        }
        tiny_loadxsol (&lp, newecount, ecount, xsol, perm);
        rval = 0;
    }

CLEANUP:

    CC_IFFREE (perm, int);
    CC_IFFREE (newelist, int);
    CC_IFFREE (newelen, int);
    CC_IFFREE (newlower, int);
    CC_IFFREE (newupper, int);

    free_tinylp (&lp);

    /* printf ("T"); fflush (stdout); */

    return rval;
}

static void tiny_loadxsol (tiny_lp *lp, int newecount, int ecount, int *xsol,
        int *perm)
{
    int i;

    if (xsol) {
        for (i = 0; i < newecount; i++) {
            xsol[perm[i]] = (int) (lp->xsol[i] + CC_TINYTSP_INTTOL);
        }
        for (; i < ecount; i++) {
            xsol[perm[i]] = 0;
        }
    }
}

static int permute_edges (int ecount, int *elist, int *elen,
        int *lower, int *upper, int *newecount, int **newelist,
        int **newelen, int **newlower, int **newupper, int **perm)
{
    int i, k, head, tail;

    *newecount = ecount;
    *newlower = (int *) NULL;
    *newupper = (int *) NULL;

    *newelist = CC_SAFE_MALLOC (2 * ecount, int);
    *newelen = CC_SAFE_MALLOC (ecount, int);
    *perm = CC_SAFE_MALLOC (ecount, int);
    if (!(*newelist) || !(*newelen) || !(*perm)) goto CLEANUP;

    if (lower) {
        *newlower = CC_SAFE_MALLOC (ecount, int);
        if (!(*newlower)) goto CLEANUP;
    }
    if (upper) {
        *newupper = CC_SAFE_MALLOC (ecount, int);
        if (!(*newupper)) goto CLEANUP;
    }

    head = 0;
    tail = ecount;

    for (i = 0; i < ecount; i++) {
        if (!upper || upper[i] != 0) {
            k = head++;
        } else {
            k = --tail;
        }
        (*newelist)[2*k] = elist[2*i];
        (*newelist)[2*k+1] = elist[2*i+1];
        (*newelen)[k] = elen[i];
        if (lower) {
            (*newlower)[k] = lower[i];
        }
        if (upper) {
            (*newupper)[k] = upper[i];
        }
        (*perm)[k] = i;
    }

    *newecount = head;
    return 0;

CLEANUP:

    fprintf (stderr, "out of memory in permute_depot_edges\n");
    CC_IFFREE (*newelist, int);
    CC_IFFREE (*newelen, int);
    CC_IFFREE (*perm, int);
    CC_IFFREE (*newlower, int);
    CC_IFFREE (*newupper, int);

    return 1;
}

static int init_tinylp (tiny_lp *lp, int ncount, int ecount, int *elist,
        int *elen, int *lower, int *upper, double *upperbound,
        int checkresult, int depot, int searchlimit, int objsense)
{
    int rval = 0;
    int i;

    init_tiny_lp_struct (lp);

    lp->objsense = objsense;
    if (upperbound) {
        if (lp->objsense == CC_TINYTSP_MAXIMIZE) {
            lp->upperbound = -(*upperbound);
        } else {
            lp->upperbound =   *upperbound;
        }
    }

    rval = build_graph (&(lp->graph), ncount, ecount, elist, elen,
                        objsense);
    if (rval) {
        fprintf (stderr, "build_graph failed\n"); goto CLEANUP;
    }
    rval = init_tinycomp (&(lp->comp), ncount);
    if (rval) {
        fprintf (stderr, "init_tinycomp failed\n"); goto CLEANUP;
    }

    lp->lower = CC_SAFE_MALLOC (ecount, double);
    lp->upper = CC_SAFE_MALLOC (ecount, double);
    if (!lp->upper || !lp->lower) {
        fprintf (stderr, "out of memory in init_tinylp\n");
        rval = 1; goto CLEANUP;
    }
    tiny_loadbounds (lp, lower, upper);

    lp->x    = CC_SAFE_MALLOC (ecount, double);
    lp->xsol = CC_SAFE_MALLOC (ecount, double);
    if (!lp->x || !lp->xsol) {
        fprintf (stderr, "out of memory in init_tinylp\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ecount; i++) {
        lp->x[i] = 0.0;
        lp->xsol[i] = 0.0;
    }

    if (checkresult) {
        lp->checkresult = 1;
        lp->node_pi   = CC_SAFE_MALLOC (ncount, CCbigguy);
        lp->pi_double = CC_SAFE_MALLOC (ncount, double);
        if (!lp->node_pi || !lp->pi_double) {
            fprintf (stderr, "out of memory in init_tinylp\n");
            rval = 1; goto CLEANUP;
        }
        lp->pisize_double = ncount;
    }

    lp->depot = depot;
    if (searchlimit < 0) {
        lp->searchlimit = CC_TINYTSP_SEARCHLIMIT;
    } else {
        lp->searchlimit = searchlimit;
    }

    rval = CClp_init (&(lp->lp));
    if (rval) {
        fprintf (stderr, "CClp_init failed\n"); goto CLEANUP;
    }
    rval = load_tinylp (lp);
    if (rval) {
        fprintf (stderr, "load_tinylp failed\n"); goto CLEANUP;
    }
    rval = CClp_tune_small (lp->lp);
    if (rval) {
        fprintf (stderr, "CClp_tune_small failed\n"); goto CLEANUP;
    }

    return rval;

CLEANUP:

    free_tinylp (lp);
    return rval;
}

static void tiny_loadbounds (tiny_lp *lp, int *lower, int *upper)
{
    int i;
    int ecount = lp->graph.ecount;

    if (lower && upper) {
        for (i = 0; i < ecount; i++) {
            lp->lower[i] = (double) lower[i];
            lp->upper[i] = (double) upper[i];
        }
    } else  {
        for (i = 0; i < ecount; i++) {
            lp->lower[i] = 0.0;
            lp->upper[i] = 1.0;
        }
    }
}

static int load_tinylp (tiny_lp *lp)
{
    int rval = 0;
    int     i, j;
    double  obj[1], lb[1], ub[1];
    int     cmatbeg[1], cmatind[2];
    double  cmatval[2];

    rval = CClp_create (lp->lp, "tinylp");
    if (rval) {
        fprintf (stderr, "CClp_create failed\n");
        goto CLEANUP;
    }

    /* Create the rows */

    for (i = 0; i < lp->graph.ncount; i++) {
        rval = CClp_new_row (lp->lp, 'E', 2.0);
        if (rval) {
            fprintf (stderr, "CClp_new_row failed\n");
            goto CLEANUP;
        }
    }
    if (lp->depot != -1) {
        rval = CClp_change_sense (lp->lp, lp->depot, 'G');
        if (rval) {
            fprintf (stderr, "CClp_change_sense failed\n"); goto CLEANUP;
        }
    }

    /* Fill in the columns */

    cmatbeg[0] = 0;
    cmatval[0] = 1.0;
    cmatval[1] = 1.0;
    for (j = 0; j < lp->graph.ecount; j++) {
        obj[0]     = (double) lp->graph.edgelist[j].len;
        lb[0]      = lp->lower[j];
        ub[0]      = lp->upper[j];
        cmatind[0] = lp->graph.edgelist[j].ends[0];
        cmatind[1] = lp->graph.edgelist[j].ends[1];
        rval = CClp_addcols (lp->lp, 1 /* no. new variables */,
                             2 /* no. new nonzeros */,
                           obj, cmatbeg, cmatind, cmatval,
                           lb, ub);
        if (rval) {
            fprintf (stderr, "CClp_addcols failed\n"); goto CLEANUP;
        }
    }

CLEANUP:
    return rval;
}

static int tiny_setbounds (tiny_lp *lp, int col, int lower, int upper)
{
    int rval = 0;

    lp->lower[col] = (double) lower;
    lp->upper[col] = (double) upper;

    rval = CClp_setbnd (lp->lp, col, 'L', (double) lower);
    if (rval) {
        fprintf (stderr, "CClp_setbnd failed\n"); return rval;
    }
    rval = CClp_setbnd (lp->lp, col, 'U', (double) upper);
    if (rval) {
        fprintf (stderr, "CClp_setbnd failed\n"); return rval;
    }
    rval = optimize_tinylp (lp);
    if (rval && rval != 2) {
        fprintf (stderr, "optimize_tinylp failed\n"); return rval;
    }
    return 0;
}

static int add_tinycut_list (tiny_lp *lp, tinycut **list, int *nadded)
{
    int rval = 0;
    int added;
    tinycut *c, *cnext;

    if (nadded) *nadded = 0;
    for (c = *list; c; c = cnext) {
        cnext = c->next;
        if (lp->depot != -1 || c->count <= lp->graph.ncount / 2) {
            rval = add_tinycut (lp, c, &added);
            if (rval) {
                fprintf (stderr, "add_tinycut failed\n"); goto CLEANUP;
            }
            if (added) {
                if (nadded) (*nadded)++;
                if (lp->ncuts >= lp->cutsize) {
                    rval = CCutil_reallocrus_scale ((void **) &lp->cuts,
                       &lp->cutsize, lp->ncuts + 1, 1.3,
                       sizeof (tinycut));
                    if (rval) {
                        fprintf (stderr, "CCutil_reallocrus_scale failed\n");
                        rval = 1; goto CLEANUP;
                    }
                }
                lp->cuts[lp->ncuts++] = *c;
            } else {
                free_tinycut (c);
            }
        } else {
            free_tinycut (c);
        }
        CC_FREE (c, tinycut);
    }
    *list = (tinycut *) NULL;

CLEANUP:

    return rval;
}

static int add_tinycut (tiny_lp *lp, tinycut *c, int *added)
{
    int rval = 0;
    tinygraph *g = &lp->graph;
    tinynode *nodelist = lp->graph.nodelist;
    tinyedge *edgelist = lp->graph.edgelist;
    tinyedge *nzlist = (tinyedge *) NULL;
    tinyadj  *a;
    int i, j, k, deg, nzcount = 0;
    double rhs[1];
    char   sense[1];
    int    rmatbeg[1];
    int    *rmatind = (int *) NULL;
    double *rmatval = (double *) NULL;

    if (added) *added = 0;

    /* All cuts added to LP go through this call, so we can check here */
    /* that cuts are legal.                                            */

    /*****************   Checking the cut  *****************/

    if (c->count <= 2) {
        fprintf (stderr, "TT Warning: Handle with less than 3 nodes\n");
        fprintf (stderr, "TT HANDLE: ");
        for (i = 0; i < c->count; i++) {
            printf ("%d ",  c->nodes[i]);
        }
        printf ("  Depot = %d  ", lp->depot);
        if (c->teeth) {
            printf (" Nteeth = %d\n", c->tcount);
        } else {
            printf ("\n");
        }
        fflush (stdout);
        return 0;
    }

    if (c->teeth) {
        if (c->tcount % 2 == 0) {
            fprintf (stderr, "TT Warning: Even number of teeth\n");
            return 0;
        }
        if (c->tcount <= 2) {
            /*
            fprintf (stderr, "TT Warning: Blossom with less than 3 teeth\n");
            */
            return 0;
        }
        g->magiclabel++;
        for (i = 0; i < c->tcount; i++) {
            /*
            if (c->teeth[i].ends[0] == lp->depot ||
                c->teeth[i].ends[1] == lp->depot) {
                fprintf (stderr, "TT Warning: Tooth contains the depot\n");
                return 0;
            }
            */
            if (nodelist[c->teeth[i].ends[0]].magic == g->magiclabel &&
                nodelist[c->teeth[i].ends[1]].magic == g->magiclabel) {
                /*
                fprintf (stderr, "TT Warning: Possible duplicate tooth\n");
                */
                return 0;
            } else {
                nodelist[c->teeth[i].ends[0]].magic = g->magiclabel;
                nodelist[c->teeth[i].ends[1]].magic = g->magiclabel;
            }
        }
        g->magiclabel++;
        for (i = 0; i < c->count; i++) {
            nodelist[c->nodes[i]].magic = g->magiclabel;
        }
        for (i = 0; i < c->tcount; i++) {
            if (nodelist[c->teeth[i].ends[0]].magic != g->magiclabel &&
                nodelist[c->teeth[i].ends[1]].magic != g->magiclabel) {
                fprintf (stderr, "TT Warning: Tooth outside of Handle\n");
                return 0;
            }
            if (nodelist[c->teeth[i].ends[0]].magic == g->magiclabel &&
                nodelist[c->teeth[i].ends[1]].magic == g->magiclabel) {
                fprintf (stderr, "TT Warning: Tooth inside of Handle\n");
                return 0;
            }
        }
    }

    /*****************   Checking complete  *****************/

    g->magiclabel++;
    for (i = 0; i < c->count; i++) {
        nodelist[c->nodes[i]].magic = g->magiclabel;
    }
    for (i = 0; i < c->count; i++) {
        a = nodelist[c->nodes[i]].adj;
        deg = nodelist[c->nodes[i]].deg;
        for (j = 0; j < deg; j++) {
            if (nodelist[a[j].to].magic != g->magiclabel) {
                if (!g->edgelist[a[j].edge].coef) {
                    g->edgelist[a[j].edge].next = nzlist;
                    nzlist = &(g->edgelist[a[j].edge]);
                    nzcount++;
                }
                g->edgelist[a[j].edge].coef++;
            }
        }
    }
    if (c->teeth) {
        for (k = 0; k < c->tcount; k++) {
            g->magiclabel++;
            for (i = 0; i < 2; i++) {
                nodelist[c->teeth[k].ends[i]].magic = g->magiclabel;
            }
            for (i = 0; i < 2; i++) {
                a = nodelist[c->teeth[k].ends[i]].adj;
                deg = nodelist[c->teeth[k].ends[i]].deg;
                for (j = 0; j < deg; j++) {
                    if (nodelist[a[j].to].magic != g->magiclabel) {
                        if (!g->edgelist[a[j].edge].coef) {
                            g->edgelist[a[j].edge].next = nzlist;
                            nzlist = &(g->edgelist[a[j].edge]);
                            nzcount++;
                        }
                        g->edgelist[a[j].edge].coef++;
                    }
                }
            }
        }
    }

    rmatind = CC_SAFE_MALLOC (nzcount, int);
    rmatval = CC_SAFE_MALLOC (nzcount, double);
    if (!rmatind || !rmatval) {
        fprintf (stderr, "out of memory in add_tinycut\n");
        for (; nzlist; nzlist = nzlist->next) {
            nzlist->coef = 0;
        }
        goto CLEANUP;
    }
    for (i = 0; nzlist; nzlist = nzlist->next, i++) {
        rmatind[i] = (int) (nzlist - edgelist);
        rmatval[i] = (double) nzlist->coef;
        nzlist->coef = 0;
    }
    if (c->teeth) {
        rhs[0] = 3.0 * ((double) c->tcount) + 1.0;
    } else {
        rhs[0] = 2.0;
    }
    sense[0] = 'G';
    rmatbeg[0] = 0;

    rval = CClp_addrows (lp->lp, 1, nzcount, rhs, sense, rmatbeg, rmatind,
                         rmatval);
    if (rval) {
        fprintf (stderr, "CClp_addrows failed\n"); goto CLEANUP;
    }
    if (added) *added = 1;

CLEANUP:

    CC_IFFREE (rmatind, int);
    CC_IFFREE (rmatval, double);

    return rval;
}

static int optimize_tinylp (tiny_lp *lp)
{
    int rval;

    rval = CClp_opt (lp->lp, CClp_METHOD_DUAL);
    if (rval == 2) {
        lp->val = CC_TINYTSP_MAXDOUBLE;
        return 2;
    }
    if (rval) {
        fprintf (stderr, "CClp_opt failed\n");
        return 1;
    }

    rval = CClp_objval (lp->lp, &lp->val);
    if (rval) {
        fprintf (stderr, "CClp_objval failed\n");
        return rval;
    }

    rval = CClp_x (lp->lp, lp->x);
    if (rval) {
        fprintf (stderr, "CClp_x failed\n");
        return rval;
    }
    return 0;
}

static int tiny_pi (tiny_lp *lp)
{
    int rval = 0;
    int i, ncount = lp->graph.ncount;

    if (lp->ncuts > lp->cutpisize) {
        rval = CCutil_reallocrus_scale ((void **) &lp->cut_pi,
           &lp->cutpisize, lp->ncuts, 1.3, sizeof (CCbigguy));
        if (rval) {
            fprintf (stderr, "CCutil_reallocrus_scale failed\n");
            rval = 1; goto CLEANUP;
        }
    }
    if (ncount + lp->ncuts > lp->pisize_double) {
        rval = CCutil_reallocrus_scale ((void **) &lp->pi_double,
           &lp->pisize_double, ncount + lp->ncuts, 1.3, sizeof (double));
        if (rval) {
            fprintf (stderr, "CCutil_reallocrus_scale failed\n");
            rval = 1; goto CLEANUP;
        }
    }

    rval = CClp_pi (lp->lp, lp->pi_double);
    if (rval) {
        fprintf (stderr, "CClp_pi failed\n"); goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) {
        lp->node_pi[i] = CCbigguy_dtobigguy (lp->pi_double[i]);
    }
    for (i = 0; i < lp->ncuts; i++) {
        lp->cut_pi[i] = CCbigguy_dtobigguy (lp->pi_double[ncount + i]);
    }


CLEANUP:

    return rval;
}

static void init_tiny_lp_struct (tiny_lp *lp)
{
    lp->graph.nodelist = (tinynode *) NULL;
    lp->graph.adjspace = (tinyadj *) NULL;
    lp->graph.edgelist = (tinyedge *) NULL;
    lp->graph.ncount = 0;
    lp->graph.ecount = 0;

    lp->comp.complist = (tinycut *) NULL;
    lp->comp.grab = (int *) NULL;
    lp->comp.stack = (int *) NULL;
    lp->comp.ncomp = 0;

    lp->lp = (CClp *) NULL;
    lp->upper = (double *) NULL;
    lp->lower = (double *) NULL;
    lp->x = (double *) NULL;
    lp->xsol = (double *) NULL;
    lp->val = -CC_TINYTSP_MAXDOUBLE;
    lp->upperbound =  CC_TINYTSP_MAXDOUBLE;
    lp->lowerbound = -CC_TINYTSP_MAXDOUBLE;
    lp->cuts = (tinycut *) NULL;
    lp->cutsize = 0;
    lp->ncuts = 0;
    lp->nbbnodes = 0;
    lp->searchlimit = 0;
    lp->status = 0;
    lp->foundtour = 0;

    lp->node_pi = (CCbigguy *) NULL;
    lp->cut_pi = (CCbigguy *) NULL;
    lp->cutpisize = 0;
    lp->pi_double = (double *) NULL;
    lp->pisize_double = 0;
    lp->checkresult = 0;
    lp->depot = -1;
}

static void free_tinylp (tiny_lp *lp)
{
    int i;

    if (lp) {
        free_tinygraph (&(lp->graph));
        free_tinycomp (&(lp->comp));
        CClp_free (&(lp->lp));
        CC_IFFREE (lp->upper, double);
        CC_IFFREE (lp->lower, double);
        CC_IFFREE (lp->x, double);
        CC_IFFREE (lp->xsol, double);
        lp->upperbound =  CC_TINYTSP_MAXDOUBLE;
        lp->lowerbound = -CC_TINYTSP_MAXDOUBLE;
        if (lp->cuts) {
            for (i = 0; i < lp->ncuts; i++) {
                free_tinycut (&(lp->cuts[i]));
            }
            CC_FREE (lp->cuts, tinycut);
        }
        lp->cutsize = 0;
        lp->ncuts = 0;
        CC_IFFREE (lp->node_pi, CCbigguy);
        CC_IFFREE (lp->cut_pi, CCbigguy);
        lp->cutpisize = 0;
        CC_IFFREE (lp->pi_double, double);
        lp->pisize_double = 0;
    }
}

static void free_tinycut (tinycut *c)
{
    if (c) {
        CC_IFFREE (c->nodes, int);
        CC_IFFREE (c->teeth, tinytooth);
    }
}

static int build_graph (tinygraph *g, int ncount, int ecount, int *elist,
        int *elen, int objsense)
{
    int i;
    tinyadj *p;
    tinynode *nlist;

    g->ncount = ncount;
    g->ecount = ecount;

    g->nodelist = CC_SAFE_MALLOC (g->ncount, tinynode);
    g->adjspace = CC_SAFE_MALLOC (2*g->ecount, tinyadj);
    g->edgelist = CC_SAFE_MALLOC (g->ecount, tinyedge);
    if (!g->nodelist || !g->adjspace || !g->edgelist) {
        fprintf (stderr, "out of memory in build_graph\n");
        CC_IFFREE (g->nodelist, tinynode);
        CC_IFFREE (g->adjspace, tinyadj);
        CC_IFFREE (g->edgelist, tinyedge);
        return 1;
    }
    nlist = g->nodelist;

    for (i = 0; i < ncount; i++) {
        nlist[i].deg = 0;
        nlist[i].mark = 0;
        nlist[i].magic = 0;
    }
    for (i = 0; i < ecount; i++) {
        nlist[elist[2*i]].deg++;
        nlist[elist[2*i+1]].deg++;
    }

    p = g->adjspace;
    for (i = 0; i < ncount; i++) {
        nlist[i].adj = p;
        p += nlist[i].deg;
        nlist[i].deg = 0;
    }

    for (i = 0; i < ecount; i++) {
        nlist[elist[2*i]].adj[nlist[elist[2*i]].deg].edge = i;
        nlist[elist[2*i]].adj[nlist[elist[2*i]].deg++].to = elist[2*i+1];
        nlist[elist[2*i+1]].adj[nlist[elist[2*i+1]].deg].edge = i;
        nlist[elist[2*i+1]].adj[nlist[elist[2*i+1]].deg++].to = elist[2*i];
    }

    if (objsense == CC_TINYTSP_MAXIMIZE) {
        for (i = 0; i < ecount; i++) {
            g->edgelist[i].ends[0] = elist[2*i];
            g->edgelist[i].ends[1] = elist[2*i+1];
            g->edgelist[i].len = -elen[i];
            g->edgelist[i].coef = 0;
        }
    } else {
        for (i = 0; i < ecount; i++) {
            g->edgelist[i].ends[0] = elist[2*i];
            g->edgelist[i].ends[1] = elist[2*i+1];
            g->edgelist[i].len =  elen[i];
            g->edgelist[i].coef = 0;
        }
    }

    g->magiclabel = 0;

    return 0;
}

static void free_tinygraph (tinygraph *g)
{
    if (g) {
        CC_IFFREE (g->nodelist, tinynode);
        CC_IFFREE (g->adjspace, tinyadj);
        CC_IFFREE (g->edgelist, tinyedge);
    }
    g->ncount = 0;
    g->ecount = 0;
}

static int tiny_connectcut (tiny_lp *lp, tinycomp *t, int *nadded)
{
    int rval = 0;
    int added;

    if (nadded) *nadded = 0;
    do {
        rval = tiny_connect (&lp->graph, t, lp->x, lp->depot);
        if (rval) {
            fprintf (stderr, "tiny_connect failed\n"); goto CLEANUP;
        }
#ifdef TINYNOISY
        printf ("Number of Components: %d\n", t->ncomp); fflush (stdout);
#endif

        if (t->ncomp > 1) {
            rval = add_tinycut_list (lp, &t->complist, &added);
            if (rval) {
                fprintf (stderr, "add_tinycut_list failed\n"); goto CLEANUP;
            }
            if (added) {
                if (nadded) (*nadded) += added;
                rval = optimize_tinylp (lp);
                if (rval == 2) {
#ifdef TINYNOISY
                    printf ("LP is infeasible\n"); fflush (stdout);
#endif
                    rval = 0; goto CLEANUP;
                } else if (rval) {
                    fprintf (stderr, "optimize_tinylp failed\n"); goto CLEANUP;
                } else {
#ifdef TINYNOISY
                    printf ("LP Value: %.4f\n", lp->val); fflush (stdout);
#endif
                }
            } else {
                fprintf (stderr, "error in connect loop\n");
                rval = 1; goto CLEANUP;
            }
        }
    } while (t->ncomp > 1);

CLEANUP:
    return rval;
}

static int init_tinycomp (tinycomp *t, int ncount)
{
    t->ncomp = 0;
    t->complist = (tinycut *) NULL;
    t->grab = CC_SAFE_MALLOC (ncount, int);
    t->stack = CC_SAFE_MALLOC (ncount, int);
    if (!t->grab || !t->stack) {
        fprintf (stderr, "out of memory in init_tinycomp\n");
        CC_IFFREE (t->grab, int);
        CC_IFFREE (t->stack, int);
        return 1;
    }
    return 0;
}

static void free_tinycomp (tinycomp *t)
{
    if (t) {
        t->ncomp = 0;
        t->complist = (tinycut *) NULL;
        CC_IFFREE (t->grab, int);
        CC_IFFREE (t->stack, int);
    }
}

static int tiny_connect (tinygraph *g, tinycomp *tc, double *x, int depot)
{
    int i, j, count;
    int ncount = g->ncount;
    tinycut *c;

    tc->ncomp = 0;
    tc->complist = (tinycut *) NULL;
    for (i = 0; i < ncount; i++) {
        g->nodelist[i].mark = 0;
    }

    if (depot != -1) {
        tc->ncomp++;
        connect_search (g, depot, tc->ncomp, tc->stack, x, tc->grab, &count);
    }

    for (i = 0; i < ncount; i++) {
        if (!g->nodelist[i].mark) {
            tc->ncomp++;
            connect_search (g, i, tc->ncomp, tc->stack, x, tc->grab, &count);
            if (count < ncount && count > 2 /* for biconnect case */) {
                c = CC_SAFE_MALLOC (1, tinycut);
                if (!c) {
                    fprintf (stderr, "out of memory in tiny_connect\n");
                    return 1;
                }
                c->teeth = (tinytooth *) NULL;
                c->tcount = 0;
                c->nodes = CC_SAFE_MALLOC (count, int);
                if (!c->nodes) {
                    fprintf (stderr, "out of memory in tiny_connect\n");
                    CC_FREE (c, tinycut);
                    return 1;
                }
                for (j = 0; j < count; j++) {
                    c->nodes[j] = tc->grab[j];
                }
                c->count = count;
                c->next = tc->complist;
                tc->complist = c;
            }
        }
    }
    return 0;
}

static void connect_search (tinygraph *g, int n, int marker, int *dstack,
        double *x, int *comp, int *count)
{
    int i, k, cnt = 0, head = 0;
    tinynode *nodelist = g->nodelist;


    nodelist[n].mark = marker;
    comp[cnt++] = n;
    dstack[head++] = n;

    while (head > 0) {
        n = dstack[--head];
        for (i = 0; i < nodelist[n].deg; i++) {
            if (x[nodelist[n].adj[i].edge] > CC_TINYTSP_EPSILON) {
                k = nodelist[n].adj[i].to;
                if (!nodelist[k].mark) {
                    nodelist[k].mark = marker;
                    comp[cnt++] = k;
                    dstack[head++] =  k;
                }
            }
        }
    }
    *count = cnt;
}

typedef struct exactsub_param {
    int             nodecount;
    int             cutcount;
    int             depot;
    int             ncount;
    tinycut        *cuts;
} exactsub_param;

static int exact_subtours (tiny_lp *lp)
{
    int rval = 0;
    int i, added;
    int *elist = (int *) NULL;
    tinygraph *g = &lp->graph;
    exactsub_param p;

    elist = CC_SAFE_MALLOC (2*g->ecount, int);
    if (!elist) {
        fprintf (stderr, "out of memory in exact_subtours\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < g->ecount; i++) {
        elist[2*i]   = g->edgelist[i].ends[0];
        elist[2*i+1] = g->edgelist[i].ends[1];
    }

    p.nodecount = g->ncount;
    p.cutcount  = 0;
    p.cuts      = (tinycut *) NULL;
    p.depot     = lp->depot;
    p.ncount    = lp->graph.ncount;

    rval = CCcut_violated_cuts (g->ncount, g->ecount, elist, lp->x,
                   2.0 - 0.0001, add_exact, (void *) &p);
    if (rval) {
        fprintf (stderr, "CCcut_violated_cuts failed\n"); goto CLEANUP;
    }

    if (p.cutcount > 0) {
        rval = add_tinycut_list (lp, &p.cuts, &added);
        if (rval) {
            fprintf (stderr, "add_tinycut_list failed\n"); goto CLEANUP;
        }
        if (added) {
            rval = optimize_tinylp (lp);
            if (rval == 2) {
#ifdef TINYNOISY
                printf ("LP is infeasible\n"); fflush (stdout);
#endif
                rval = 0; goto CLEANUP;
            } else if (rval) {
                fprintf (stderr, "optimize_tinylp failed\n"); goto CLEANUP;
            } else {
#ifdef TINYNOISY
                printf ("LP Value: %.4f\n", lp->val); fflush (stdout);
#endif
            }
        }
    }

CLEANUP:

    CC_IFFREE (elist, int);
    return rval;
}

static int add_exact (double val, int count, int *cutarray, void *pass_param)
{
    int rval = 0;
    int j;
    exactsub_param *p = (exactsub_param *) pass_param;
    tinycut *c = (tinycut *) NULL;

    if (val > 2.0) {
        fprintf (stderr, "TT Warning: Cut of value %f in add_exact\n", val);
        goto CLEANUP;
    }

    if (count > 2 && count < p->ncount - 2) {
        for (j = 0; j < count; j++) {
            if (cutarray[j] == p->depot) {
                goto CLEANUP; 
            }
        }

        c = CC_SAFE_MALLOC (1, tinycut);
        if (!c) {
            fprintf (stderr, "out of memory in add_exact\n");
            rval = 1; goto CLEANUP;
        }
        c->teeth = (tinytooth *) NULL;
        c->nodes = (int *) NULL;
        c->tcount = 0;

        c->nodes = CC_SAFE_MALLOC (count, int);
        if (!c->nodes) {
            fprintf (stderr, "out of memory in add_exact\n");
            rval = 1; goto CLEANUP;
        }
        for (j = 0; j < count; j++) {
            c->nodes[j] = cutarray[j];
        }
        c->count = count;

        c->next = p->cuts;
        p->cuts = c;
        p->cutcount++;
    }

CLEANUP:

    if (rval) {
        free_tinycut (c);
        CC_IFFREE (c, tinycut);
    }

    return rval;
}

#define X_FLUFF (1e-10)

static int exact_blossoms (tiny_lp *lp)
{
    int ecount, cutcount;
    int *elist = (int *) NULL;
    double *x = (double *) NULL;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    CCrandstate rstate;
    int seed = 99;
    int rval = 0;

    CCutil_sprand (seed, &rstate);

    rval = grab_nonzero_x (lp, &ecount, &elist, &x, X_FLUFF);
    if (rval) {
        fprintf (stderr, "grab_nonzero_x failed\n"); goto CLEANUP;
    }

    rval =  CCtsp_exactblossom (&cuts, &cutcount, lp->graph.ncount, ecount,
                                elist, x, &rstate);
    if (rval) {
        fprintf (stderr, "CCtsp_exactblossom failed\n"); goto CLEANUP;
    }

    if (cutcount) {
        CCtsp_lpcut_in *cut, *cnext;
        tinycut *tc, *tcuts = (tinycut *) NULL;
        int added;

        for (cut = cuts; cut; cut = cnext) {
            cnext = cut->next;
            rval = lpcut_in_to_tinycut (cut, &tc);
            if (rval) {
                fprintf (stderr, "lpcut_in_to_tinycut failed\n"); goto CLEANUP;
            }
            if (tc) {
                tc->next = tcuts;
                tcuts = tc;
            }
            CCtsp_free_lpcut_in (cut);
            CC_FREE (cut, CCtsp_lpcut_in);
        }
        rval = add_tinycut_list (lp, &tcuts, &added);
        if (rval) {
            fprintf (stderr, "add_tinycut_list failed\n"); goto CLEANUP;
        }
        if (added) {
            rval = optimize_tinylp (lp);
            if (rval == 2) {
#ifdef TINYNOISY
                printf ("LP is infeasible\n"); fflush (stdout);
#endif
                rval = 0; goto CLEANUP;
            } else if (rval) {
                fprintf (stderr, "optimize_tinylp failed\n"); goto CLEANUP;
            } else {
#ifdef TINYNOISY
                printf ("LP Value: %.4f\n", lp->val); fflush (stdout);
#endif
            }
        }
    }

CLEANUP:

    CC_IFFREE (x, double);
    CC_IFFREE (elist, int);

    return rval;
}

static int grab_nonzero_x (tiny_lp *lp, int *ecount, int **elist, double **x,
       double tol)
{
    int i, count;
    int gecount = lp->graph.ecount;
    int rval = 0;

    *ecount = 0;
    *elist = (int *) NULL;
    *x     = (double *) NULL;

    for (i = 0, count = 0; i < gecount; i++) {
        if (lp->x[i] > tol) count++;
    }

    *elist = CC_SAFE_MALLOC (2*count, int);
    *x = CC_SAFE_MALLOC (count, double);
    if (!(*elist) || !(*x)) {
        fprintf (stderr, "out of memory in grab_nonzero_x\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0, count = 0; i < gecount; i++) {
        if (lp->x[i] > tol) {
            (*elist)[2*count]   = lp->graph.edgelist[i].ends[0];
            (*elist)[2*count+1] = lp->graph.edgelist[i].ends[1];
            (*x)[count] = lp->x[i];
            count++;
        }
    }
    *ecount = count;

CLEANUP:

    if (rval) {
        CC_IFFREE (*elist, int);
        CC_IFFREE (*x, double);
        *ecount = 0;
    }
    return rval;
}

static int lpcut_in_to_tinycut (CCtsp_lpcut_in *in, tinycut **out)
{
    tinycut *c = (tinycut *) NULL;
    int rval = 0;
    int i, icount;
    int *ar = (int *) NULL;

    *out = (tinycut *) NULL;

    c = CC_SAFE_MALLOC (1, tinycut);
    if (!c) {
        fprintf (stderr, "out of memory in lpcut_in_to_tinycut\n");
        rval = 1; goto CLEANUP;
    }

    c->teeth = (tinytooth *) NULL;
    c->nodes = (int *) NULL;
    c->tcount = 0;
    c->count = 0;
    c->next  = (tinycut *) NULL;

    rval = CCtsp_clique_to_array (&in->cliques[0], &c->nodes, &c->count);
    if (rval) {
        fprintf (stderr, "CCtsp_clique_to_array failed\n"); goto CLEANUP;
    }

    if (in->cliquecount > 1) {
        c->teeth = CC_SAFE_MALLOC (in->cliquecount - 1, tinytooth);
        if (!c->teeth) {
            fprintf (stderr, "out of memory in lpcut_in_to_tinycut\n");
            rval = 1; goto CLEANUP;
        }
        for (i = 1; i < in->cliquecount; i++) {
            rval = CCtsp_clique_to_array (&in->cliques[i], &ar, &icount);
            if (rval) {
                fprintf (stderr, "CCtsp_clique_to_array failed\n");
                goto CLEANUP;
            }
            if (icount != 2) {   /* Not a blossom */
                printf ("Not a blossom\n"); fflush (stdout);
                free_tinycut (c);
                CC_IFFREE (c, tinycut);
                goto CLEANUP;
            }
            c->teeth[c->tcount].ends[0] = ar[0];
            c->teeth[c->tcount].ends[1] = ar[1];
            c->tcount++;
            CC_FREE (ar, int);
        }
    }

    *out = c;

CLEANUP:
    if (rval) {
        free_tinycut (c);
        CC_IFFREE (c, tinycut);
    }
    CC_IFFREE (ar, int);

    return rval;
}

static int tiny_blossom (tiny_lp *lp, int *nadded)
{
    tinygraph *g = &lp->graph;
    tinycomp *tc = &lp->comp;
    int ncount = g->ncount;
    int ecount = g->ecount;
    int i, count, parity;
    int rval = 0;

    if (nadded) *nadded = 0;

    tc->ncomp = 0;
    tc->complist = (tinycut *) NULL;
    for (i = 0; i < ncount; i++) {
        g->nodelist[i].mark = 0;
        g->nodelist[i].parity = 0;
    }
    for (i = 0; i < ecount; i++) {
        if (lp->x[i] >= 1 - CC_TINYTSP_EPSILON) {
            g->nodelist[g->edgelist[i].ends[0]].parity ^= 1;
            g->nodelist[g->edgelist[i].ends[1]].parity ^= 1;
        }
    }

    if (lp->depot != -1) {
        tc->ncomp++;
        blossom_search (g, lp->depot, tc->ncomp, tc->stack, lp->x, tc->grab,
                        &count, &parity);
    }

    for (i = 0; i < ncount; i++) {
        if (!g->nodelist[i].mark) {
            tc->ncomp++;
            blossom_search (g, i, tc->ncomp, tc->stack, lp->x, tc->grab,
                            &count, &parity);
            if (count < ncount && count > 2 && parity == 1) {
                rval = blossom_grab (g, lp->depot, lp->x, tc->grab, count,
                                     &tc->complist, (int *) NULL);
                if (rval) {
                    fprintf (stderr, "blossom_grab failed\n");
                    return rval;
                }
            }
        }
    }

    if (tc->complist) {
        rval = add_tinycut_list (lp, &tc->complist, nadded);
        if (rval) {
            fprintf (stderr, "add_tinycut_list failed\n");
            return 0;
        }
    }
    if (*nadded) {
        rval = optimize_tinylp (lp);
        if (rval == 2) {
#ifdef TINYNOISY
            printf ("LP is infeasible\n"); fflush (stdout);
#endif
            return 0;
        } else if (rval) {
            fprintf (stderr, "optimize_tinylp failed\n");
            return rval;
        } else {
#ifdef TINYNOISY
            printf ("Blossom LP Value: %.4f\n", lp->val);
            fflush (stdout);
#endif
        }
    }

    return 0;
}

static void blossom_search (tinygraph *g, int n, int marker, int *dstack,
        double *x, int *comp, int *count, int *parity)
{
    int i, k, cnt = 0, head = 0, par = 0;
    tinynode *nodelist = g->nodelist;


    nodelist[n].mark = marker;
    comp[cnt++] = n;
    par += nodelist[n].parity;
    dstack[head++] = n;

    while (head > 0) {
        n = dstack[--head];
        for (i = 0; i < nodelist[n].deg; i++) {
            if (x[nodelist[n].adj[i].edge] > CC_TINYTSP_EPSILON &&
                x[nodelist[n].adj[i].edge] < 1 - CC_TINYTSP_EPSILON) {
                k = nodelist[n].adj[i].to;
                if (!nodelist[k].mark) {
                    nodelist[k].mark = marker;
                    comp[cnt++] = k;
                    par += nodelist[k].parity;
                    dstack[head++] =  k;
                }
            }
        }
    }
    *count = cnt;
    *parity = par % 2;
}

static int blossom_grab (tinygraph *g, int depot, double *x, int *handle,
        int count, tinycut **list, int *foundcut)
{
    int i, j, tcount = 0;
    tinynode *n;
    tinycut *c;

    if (foundcut) *foundcut = 0;

    for (i = 0; i < count; i++) {
        if (handle[i] == depot) {
            return 0;
        }
        if (g->nodelist[handle[i]].parity) {
            n = &(g->nodelist[handle[i]]);
            for (j = 0; j < n->deg; j++) {
                if (x[n->adj[j].edge] >= 1 - CC_TINYTSP_EPSILON) {
                    if (n->adj[j].to == depot) {
                        return 0;
                    } else if (g->nodelist[n->adj[j].to].mark != n->mark) {
                        tcount++;
                    }
                }
            }
        }
    }

    if (tcount % 2 == 0) {
        fprintf (stderr, "TT Warning: Blossom with even number of teeth\n");
        return 0;
    }

    if (tcount < 3) {
        /* Only one tooth, no need to add the blossom */
        return 0;
    }

    /* depot is not in handle or teeth */

    c = CC_SAFE_MALLOC (1, tinycut);
    if (!c) {
        fprintf (stderr, "out of memory in blossom_grab\n");
        return 1;
    }
    c->nodes = CC_SAFE_MALLOC (count, int);
    if (!c->nodes) {
        fprintf (stderr, "out of memory in blossom_grab\n");
        CC_FREE (c, tinycut);
        return 1;
    }
    for (i = 0; i < count; i++) {
        c->nodes[i] = handle[i];
    }
    c->count = count;

    c->teeth = CC_SAFE_MALLOC (tcount, tinytooth);
    if (!c->teeth) {
        fprintf (stderr, "out of memory in blossom_grab\n");
        CC_FREE (c->nodes, int);
        CC_FREE (c, tinycut);
        return 1;
    }
    c->tcount = tcount;
    tcount = 0;
    for (i = 0; i < count; i++) {
        if (g->nodelist[handle[i]].parity) {
            n = &(g->nodelist[handle[i]]);
            for (j = 0; j < n->deg; j++) {
                if (x[n->adj[j].edge] >= 1 - CC_TINYTSP_EPSILON &&
                     g->nodelist[n->adj[j].to].mark != n->mark) {
                    c->teeth[tcount].ends[0] =
                                  g->edgelist[n->adj[j].edge].ends[0];
                    c->teeth[tcount].ends[1] =
                                  g->edgelist[n->adj[j].edge].ends[1];
                    tcount++;
                }
            }
        }
    }
    c->next = *list;
    *list = c;
    if (foundcut) *foundcut = 1;

    return 0;
}

static int tiny_brancher (tiny_lp *lp, int depth)
{
    int rval = 0;
    int i, next, cutoff, nadded;
    double val;

    if (lp->nbbnodes >= lp->searchlimit) {
        lp->status = CC_TINYTSP_SEARCHLIMITEXCEEDED;
        return 0;
    }

    lp->nbbnodes++;

    rval = tiny_connectcut (lp, &lp->comp, (int *) NULL);
    if (rval) {
        fprintf (stderr, "tiny_connectcut failed\n"); return rval;
    }
    rval = tiny_checkbound (lp, &cutoff);
    if (rval) {
        fprintf (stderr, "tiny_checkbound failed\n"); return 1;
    }
    if (cutoff) return 0;

    do {
        rval = exact_subtours (lp);
        if (rval) {
            fprintf (stderr, "exact_subtour failed\n"); return rval;
        }
        rval = tiny_checkbound (lp, &cutoff);
        if (rval) {
            fprintf (stderr, "tiny_checkbound failed\n"); return 1;
        }
        if (cutoff) return 0;

        rval = tiny_blossom (lp, &nadded);
        if (rval) {
            fprintf (stderr, "tiny_blossom failed\n"); return rval;
        }
        rval = tiny_checkbound (lp, &cutoff);
        if (rval) {
            fprintf (stderr, "tiny_checkbound failed\n"); return 1;
        }
        if (cutoff) return 0;
    } while (nadded);

    if (depth == 0) {
        rval = exact_blossoms (lp);
        if (rval) {
            fprintf (stderr, "tiny_tighten failed\n"); return rval;
        }
        rval = tiny_checkbound (lp, &cutoff);
        if (rval) {
            fprintf (stderr, "tiny_checkbound failed\n"); return 1;
        }
        if (cutoff) return 0;

        rval = exact_subtours (lp);
        if (rval) {
            fprintf (stderr, "exact_subtour failed\n"); return rval;
        }
        rval = tiny_checkbound (lp, &cutoff);
        if (rval) {
            fprintf (stderr, "tiny_checkbound failed\n"); return 1;
        }
        if (cutoff) return 0;
    }

    rval = tiny_connectcut (lp, &lp->comp, (int *) NULL);
    if (rval) {
        fprintf (stderr, "tiny_connectcut failed\n"); return rval;
    }
    rval = tiny_checkbound (lp, &cutoff);
    if (rval) {
        fprintf (stderr, "tiny_checkbound failed\n"); return 1;
    }
    if (cutoff) return 0;

#ifdef TINYNOISY
    printf ("Cut LP Value: %.4f\n", lp->val); fflush (stdout);
#endif

/*
    for (i = 0; i < lp->graph.ecount; i++) {
        if (lp->x[i] > 0.0) {
            printf ("%d %d %f\n", lp->graph.edgelist[i].ends[0],
                                  lp->graph.edgelist[i].ends[1],
                                  lp->x[i]);
        }
    }
    exit (1);
*/

    next = tiny_find_branch (&(lp->graph), lp->depot, lp->x);
    if (next == -1) {
        rval = tiny_checktour (lp, &val);
        if (rval) {
            fprintf (stderr, "tiny_checktour failed\n"); return 1;
        } else {
#ifdef TINYNOISY
            printf ("INTEGRAL SOLUTION AT %.2f\n", val);
#endif
            lp->upperbound = val;
            lp->foundtour = 1;
            for (i = 0; i < lp->graph.ecount; i++) {
                lp->xsol[i] = lp->x[i];
            }
            return 0;
        }
    }

#ifdef TINYNOISY
    printf ("Branch on %d\n", next); fflush (stdout);
#endif

    rval = tiny_setbounds (lp, next, 0, 0);
    if (rval) {
        fprintf (stderr, "tiny_setbounds failed\n"); return rval;
    }
    rval = tiny_checkbound (lp, &cutoff);
    if (rval) {
        fprintf (stderr, "tiny_checkbound failed\n"); return 1;
    }
    if (cutoff) {
#ifdef TINYNOISY
        printf ("0-Side Does not Need to be Searched: %f\n", lp->val);
        fflush (stdout);
#endif
    } else {
#ifdef TINYNOISY
        printf ("Evaluate 0-Side of %d (Depth %d, Best %.2f)\n",
                 next, depth, lp->upperbound);
        fflush (stdout);
#endif
        rval = tiny_brancher (lp, depth + 1);
        if (rval) return rval;
        if (lp->status == CC_TINYTSP_SEARCHLIMITEXCEEDED) return 0;
    }

    rval = tiny_setbounds (lp, next, 1, 1);
    if (rval) {
        fprintf (stderr, "tiny_setbounds failed\n"); return rval;
    }

    rval = tiny_checkbound (lp, &cutoff);
    if (rval) {
        fprintf (stderr, "tiny_checkbound failed\n"); return 1;
    }
    if (cutoff) {
#ifdef TINYNOISY
        printf ("1-Side of %d Does not Need to be Searched: %f\n",
                  next, lp->val);
        fflush (stdout);
#endif
    } else {
#ifdef TINYNOISY
        printf ("Evaluate 1-Side of %d (Depth %d, Best %.2f)\n",
                 next, depth, lp->upperbound);
        fflush (stdout);
#endif
        rval = tiny_brancher (lp, depth + 1);
        if (rval) return rval;
        if (lp->status == CC_TINYTSP_SEARCHLIMITEXCEEDED) return 0;
    }

    rval = tiny_setbounds (lp, next, 0, 1);
    if (rval) {
        fprintf (stderr, "tiny_setbounds failed\n"); return rval;
    }

    return 0;
}


#ifdef TINY_BEST_BRANCH
static int tiny_find_branch (tinygraph *g, int depot, double *x)
{
    int i;
    int besti = -1;
    int ecount = g->ecount;
    tinyedge *elist = g->edgelist;
    double maxdiff = 0.0;


    for (i = 0; i < ecount; i++) {
        if (x[i] > 0.5) {
            if (1.0 - x[i] > maxdiff && elist[i].ends[0] != depot &&
                                        elist[i].ends[1] != depot) {
                maxdiff = 1.0 - x[i];
                besti = i;
            }
        } else {
            if (x[i] > maxdiff && elist[i].ends[0] != depot &&
                                  elist[i].ends[1] != depot ) {
                maxdiff = x[i];
                besti = i;
            }
        }
    }

    if (maxdiff < CC_TINYTSP_INTTOL)
        return -1;
    else
        return  besti;
}
#endif

#ifdef TINY_FIRST_BRANCH
static int tiny_find_branch (tinygraph *g, int depot, double *x)
{
    int i;
    int ecount = g->ecount;
    tinyedge *elist = g->edgelist;

    for (i = 0; i < ecount; i++) {
        if (x[i] > CC_TINYTSP_INTTOL && x[i] < 1.0 - CC_TINYTSP_INTTOL &&
            elist[i].ends[0] != depot && elist[i].ends[1] != depot) {
            return i;
        }
    }
    return -1;
}
#endif

static int tiny_checkbound (tiny_lp *lp, int *cutoff)
{
    int rval = 0;


    if (lp->val == CC_TINYTSP_MAXDOUBLE ||
                   lp->val > lp->upperbound - 1.0 + CC_TINYTSP_EPSILON) {
        if (lp->checkresult) {
            rval = tiny_checkdual (lp, cutoff);
            if (rval) {
                fprintf (stderr, "tiny_checkdual failed\n");
                return rval;
            }
        } else {
            *cutoff = 1;
        }
    } else {
        *cutoff = 0;
    }

    return 0;
}

static int tiny_checkdual (tiny_lp *lp, int *cutoff)
{
    int rval = 0;
    int i, k, phase1;
    tinyedge *elist = lp->graph.edgelist;
    int ecount = lp->graph.ecount;
    CCbigguy bnd, rhs_sum;

#ifdef TINYNOISY
    printf ("tiny_checkdual ...\n"); fflush (stdout);
#endif

    rval = tiny_pi (lp);
    if (rval) {
        fprintf (stderr, "tiny_pi failed\n");
        return rval;
    }

    rhs_sum = CCbigguy_ZERO;
    for (i = 0; i < lp->graph.ncount; i++) {
        CCbigguy_addmult (&rhs_sum, lp->node_pi[i], 2);
    }
    for (i = 0; i < lp->ncuts; i++) {
        if (lp->cuts[i].teeth) {
            k = 3 * lp->cuts[i].tcount + 1;
        } else {
            k = 2;
        }
        CCbigguy_addmult (&rhs_sum, lp->cut_pi[i], k);
    }

    phase1 = (lp->val == CC_TINYTSP_MAXDOUBLE);
    tiny_price (lp, phase1);
#ifdef TINYNOISY
    printf ("rhs_sum = %f\n", CCbigguy_bigguytod (rhs_sum)); fflush (stdout);
#endif
    for (i = 0; i < ecount; i++) {
        k = CCbigguy_cmp (elist[i].rc, CCbigguy_ZERO);
        if (k < 0) {
            if (lp->upper[i] > 0.0) {
                CCbigguy_addmult (&rhs_sum, elist[i].rc, (int) lp->upper[i]);
            }
        } else if (k > 0) {
            if (lp->lower[i] > 0.0) {
                CCbigguy_addmult (&rhs_sum, elist[i].rc, (int) lp->lower[i]);
            }
        }
    }

    if (phase1) {
        if (CCbigguy_cmp (rhs_sum, CCbigguy_ZERO) <= 0) {
            CClp_dump_lp (lp->lp, "dump.sav");
            fprintf (stderr, "Infeasible LP with Farkas RHS %f\n",
                   CCbigguy_bigguytod (rhs_sum));
            return 1;
        }
        *cutoff = 1;
    } else {
        if (lp->val > CCbigguy_bigguytod (rhs_sum) + 0.5) {
           fprintf (stderr, "Val: %f   Exact Val: %f\n",
                   lp->val, CCbigguy_bigguytod (rhs_sum));
           return 1;
        }
        bnd = CCbigguy_dtobigguy (lp->upperbound - 1.0);
        if (CCbigguy_cmp (rhs_sum, bnd) > 0) {
            *cutoff = 1;
        } else {
            *cutoff = 0;
        }
    }
    return 0;
}

static void tiny_price (tiny_lp *lp, int phase1)
{
    tinyedge *elist = lp->graph.edgelist;
    tinynode *nlist = lp->graph.nodelist;
    int ecount = lp->graph.ecount;
    tinygraph *g = &(lp->graph);
    int i, j, k, t, n;
    CCbigguy x;

    if (phase1) {
        for (i = 0; i < ecount; i++) {
            elist[i].rc = CCbigguy_ZERO;
        }
    } else {
        for (i = 0; i < ecount; i++) {
            elist[i].rc = CCbigguy_itobigguy (elist[i].len);
        }
    }

    for (i = 0; i < ecount; i++) {
        CCbigguy_sub (&elist[i].rc, lp->node_pi[elist[i].ends[0]]);
        CCbigguy_sub (&elist[i].rc, lp->node_pi[elist[i].ends[1]]);
    }
    for (i = 0; i < lp->ncuts; i++) {
        x = lp->cut_pi[i];
        g->magiclabel++;
        for (j = 0;  j < lp->cuts[i].count; j++) {
            nlist[lp->cuts[i].nodes[j]].magic = g->magiclabel;
        }
        for (j = 0;  j < lp->cuts[i].count; j++) {
            n = lp->cuts[i].nodes[j];
            for (k = 0; k < nlist[n].deg; k++)  {
                if (nlist[nlist[n].adj[k].to].magic != g->magiclabel) {
                    CCbigguy_sub (&elist[nlist[n].adj[k].edge].rc, x);
                }
            }
        }
        if (lp->cuts[i].teeth) {
            for (t = 0; t < lp->cuts[i].tcount; t++) {
                g->magiclabel++;
                for (j = 0; j < 2; j++) {
                    nlist[lp->cuts[i].teeth[t].ends[j]].magic = g->magiclabel;
                }
                for (j = 0; j < 2; j++) {
                    n = lp->cuts[i].teeth[t].ends[j];
                    for (k = 0; k < nlist[n].deg; k++)  {
                        if (nlist[nlist[n].adj[k].to].magic != g->magiclabel) {
                            CCbigguy_sub (&elist[nlist[n].adj[k].edge].rc, x);
                        }
                    }
                }
            }
        }
    }
}

static int tiny_checktour (tiny_lp *lp, double *tourval)
{
    int rval = 0;
    int i;
    int val = 0;
    tinygraph *g = &lp->graph;

    for (i = 0; i < g->ncount; i++) {
        g->nodelist[i].mark = 0;
    }
    for (i = 0; i < g->ecount; i++) {
        if (lp->x[i] < lp->lower[i] - CC_TINYTSP_INTTOL ||
            lp->x[i] > lp->upper[i] + CC_TINYTSP_INTTOL) {
            fprintf (stderr, "variable not between bounds\n");
            return 1;
        }
        if (lp->x[i] > CC_TINYTSP_INTTOL) {
            val += g->edgelist[i].len;
            g->nodelist[g->edgelist[i].ends[0]].mark++;
            g->nodelist[g->edgelist[i].ends[1]].mark++;
            if (lp->x[i] > 1 + CC_TINYTSP_INTTOL) {
                val += g->edgelist[i].len;
                g->nodelist[g->edgelist[i].ends[0]].mark++;
                g->nodelist[g->edgelist[i].ends[1]].mark++;
            }
        }
    }
    if (lp->val < ((double) val) - 0.5)  {
        fprintf (stderr, "LP val and computed tour length do not agree\n");
        return 1;
    }
    for (i = 0; i < g->ncount; i++) {
        if (g->nodelist[i].mark != 2 && i != lp->depot) {
            fprintf (stderr, "node in tour does not have degree 2\n");
            return 1;
        }
    }

    rval = tiny_connect (&lp->graph, &lp->comp, lp->x, lp->depot);
    if (rval) {
        fprintf (stderr, "tiny_connect failed\n"); return rval;
    }
    if (lp->comp.ncomp > 1) {
        tinycut *c, *cnext;
        for (c = lp->comp.complist; c; c = cnext) {
            cnext = c->next;
            free_tinycut (c);
            CC_FREE (c, tinycut);
        }
        fprintf (stderr, "tour is not connected\n");
        return 1;
    }

    *tourval = (double) val;
    return 0;
}


#ifdef TINYDEBUG
static void print_tinycut (tinycut *c)
{
    int i;

    printf ("CUT: "); fflush (stdout);
    for (i = 0; i < c->count; i++) {
        printf ("%d ", c->nodes[i]);
    }
    printf ("\n"); fflush (stdout);
}
#endif /* TINYDEBUG */


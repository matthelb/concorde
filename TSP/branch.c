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
/*                  ROUTINES TO EXECUTE BRANCHING                           */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: May 22, 1997                                                      */
/*  Modified: June 17, 1997 (bix)                                           */
/*            June 27, 1997 (bico)                                          */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCtsp_find_branch (CCtsp_lp *lp, int nwant, int *ngot,              */
/*      CCtsp_branchobj **bobj, double *val, int **cyc,                     */
/*      int usecliques, int longedge_branching, int silent)                 */
/*    FINDS a set of branching edges and cliques.                           */
/*     -usecliques should be set to 1 to allow branching on cliques         */
/*     -val returns the length of a tour if one is detected.                */
/*     -cyc returns the tour (it can be NULL)                               */
/*     -longedge_branching if nonzero selects the alternative longedge      */
/*      branching rule.                                                     */
/*                                                                          */
/*  int CCtsp_find_branch_edge (CCtsp_lp *lp, int *n0, int *n1,             */
/*      double *val, int **cyc, int branchtype, int silent)                 */
/*    FINDS a branching edge or detects that solution is integral.          */
/*     -lp points to an optimized lp.                                       */
/*     -n0, n1 return the edges of the branching edge; n0 is set to -1      */
/*          if the current lp solution is a tour                            */
/*     -val returns the value the tour if n0 is set to -1                   */
/*     -branchtype determines the strategy for choosing the branching       */
/*          edge; choices for branchtype are given in tsp.h                 */
/*                                                                          */
/*  int CCtsp_check_integral (CCtsp_lp *lp, double *val, int **cyc,         */
/*      int *yesno, int silent)                                             */
/*    TESTS if the current x-vector is a tour.                              */
/*     -yesno is set to 1 if it is a tour and 0 otherwise.                  */
/*                                                                          */
/*  int CCtsp_find_branch_cliques (CCtsp_lp *lp, int nwant,                 */
/*      int longedge_branching, int *ngot, CCtsp_lpclique **bcliques,       */
/*      double *bval, int silent)                                           */
/*    FINDS branching cliques (it may return ngot == 0)                     */
/*    -bval will return the stongbranching function evaluation for          */
/*     each clique (it can be NULL)                                         */
/*                                                                          */
/*  void CCtsp_init_branchobj (CCtsp_branchobj *b)                          */
/*    INITITALIZES the fields in the CCtsp_branchobj pointed to by b.       */
/*                                                                          */
/*  void CCtsp_free_branchobj (CCtsp_branchobj *b)                          */
/*    FREES the fields in the CCtsp_branchobj pointed to by b.              */
/*                                                                          */
/*  void CCtsp_print_branchhistory (CCtsp_lp *lp)                           */
/*    PRINT to stdout the branch history of the lp                          */
/*                                                                          */
/*  int CCtsp_execute_branch (CCtsp_lp *lp, CCtsp_branchobj *b,             */
/*      int silent, CCrandstate *rstate)                                    */
/*    SETS the lp to realize the branch described in b                      */
/*    NOTES: returns 2 if the LP becomes infeasible.                        */
/*                                                                          */
/*  int CCtsp_execute_unbranch (CCtsp_lp *lp, CClp_warmstart *warmstart,    */
/*      int silent, CCrandstate *rstate)                                    */
/*    UNDOS the changes to the lp caused by the most recent branch that     */
/*     has not yet been unbranched (used in dfs)                            */
/*     -warmstart can specify a warmstart to help resolve the LP            */
/*      (it can be NULL)                                                    */
/*    NOTES: returns 2 if the LP is infeasible.                             */
/*                                                                          */
/*  int CCtsp_add_branchhistory_to_lp (CCtsp_lp *lp)                        */
/*    SETS the lp to realize the branches in branch history                 */
/*                                                                          */
/*  int CCtsp_bb_find_branch (char *probname, int probnum, int ncount,      */
/*      CCdatagroup *dat, int *ptour, double *upperbound,                   */
/*      CCtsp_lpcuts *pool, int nwant, int *ngot, CCtsp_branchobj **b,      */
/*      int usecliques, int longedge_branching, int *foundtour,             */
/*      int *besttour, int silent, CCrandstate *rstate),                    */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_splitprob (CCtsp_lp *lp, CCtsp_branchobj *b, int child0,      */
/*      int child1, int silent, CCrandstate *rstate)                        */
/*    EXECUTES a branch on the lp and writes to two child lps               */
/*     -b contains the branching information (the rhs side value is set     */
/*      by this function to give 0 and 1 for edge branches and 2 & 4 for    */
/*      clique branches; the sense is set by this function for clique       */
/*      branches to realize <= 2 and >= 4)                                  */
/*     -child0 and child1 are the ids of the children                       */
/*                                                                          */
/*  int CCtsp_bb_splitprob (char *probname, int probnum, int ncount,        */
/*      CCdatagroup *dat, int *ptour, double initial_ub,                    */
/*      CCtsp_lpcuts *pool, CCtsp_branchobj *b, int child0,                 */
/*      int child1, double *val0, double *val1, int *prune0,                */
/*      int *prune1, int silent, CCranstate *rstate)                        */
/*    CALLS splitprob after reading the lp file and building an lp; this    */
/*     function will also price the lp and attempt to verify infeasible     */
/*     lps.                                                                 */
/*     -val0 and val1 return the (priced) lp-bounds for the children; if    */
/*      an lp is infeasible then the val is set to CCtsp_LP_MAXDOUBLE and   */
/*      the lp is not written.                                              */
/*     -prune0 and prune1 will be set to 1 if the child can be pruned       */
/*      (in which case the lp is not written)                               */
/*                                                                          */
/*  int CCtsp_dumptour (int ncount, CCdatagroup *dat, int *perm,            */
/*      char *probname, int *tour, char *fname, int writeedges, int silent) */
/*    WRITES the tour file to fname or to probname.sol.                     */
/*     -dat is used to compute the length (it can be NULL)                  */
/*     -perm is a permutation tour                                          */
/*     -tour gives the tour (perm[tour[i]] with be printed)                 */
/*     -fname specifies the output name (can be NULL)                       */
/*     -writeedges should be set to nonzero to write an edge file (default  */
/*      is to write a tour file)                                            */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "macrorus.h"
#include "tsp.h"
#include "lp.h"
#include "cut.h"

#undef  DEBUG
#undef  LONGEDGE_TRIANGLES

#define TSP_BRANCH_STRONG_ALL_CHOICES      (-1)
#define TSP_BRANCH_STRONG_ITERATIONS       (100)  /* 100 */
#define TSP_BRANCH_STRONG_EXTRA_ITERATIONS (500)  /* 500 */
#define TSP_BRANCH_STRONG_CHOICES_MULT     (5)    /* 5 */
#define TSP_BRANCH_STRONG_LP_CHOICES_MULT  (2)    /* 2 */

#define TSP_STRONG_CUT_CHOICES_MULT        (25)   /* 50 */
#define TSP_STRONG_CUT_LP_CHOICES_MULT     (5)    /* 5 */
#define TSP_STRONG_CUT_CANDIDATES          (1000) /* 1000 */
#define TSP_STRONG_CUT_CANDIDATES_MULT     (100)  /* 100 */

#define TSP_BRANCH_STRONG_FIRST_WEIGHT (10.0)
#define TSP_BRANCH_STRONG_WEIGHT (100.0)  
#define TSP_BRANCH_STRONG_CUT_NORM_WEIGHT (100.0)

#define TSP_BRANCH_STRONG_FIRST_VAL(v0,v1)                                 \
    (((v0) < (v1) ? (TSP_BRANCH_STRONG_FIRST_WEIGHT * (v0) + (v1))         \
                  : (TSP_BRANCH_STRONG_FIRST_WEIGHT * (v1) + (v0)))        \
                    / (TSP_BRANCH_STRONG_FIRST_WEIGHT + 1.0))

#define TSP_BRANCH_STRONG_VAL(v0,v1)                                       \
    (((v0) < (v1) ? (TSP_BRANCH_STRONG_WEIGHT * (v0) + (v1))               \
                  : (TSP_BRANCH_STRONG_WEIGHT * (v1) + (v0)))              \
                    / (TSP_BRANCH_STRONG_WEIGHT + 1.0))

#define TSP_BRANCH_STRONG_CUT_NORM_VAL(v0,v1)                              \
    (((v0) < (v1) ? (TSP_BRANCH_STRONG_CUT_NORM_WEIGHT * (v0) + (v1))      \
                  : (TSP_BRANCH_STRONG_CUT_NORM_WEIGHT * (v1) + (v0)))     \
                    / (TSP_BRANCH_STRONG_CUT_NORM_WEIGHT + 1.0))

/* a clique must be this far from 2.0 and 4.0 to be used in a branch */

#define TSP_BRANCH_CLIQUE_TOL 0.01

typedef struct sbitem {
    int    name;
    int    name2;
    int    name3;
    int    name4;
    double val;
} sbitem;

typedef struct ds_node {
    struct ds_node *parent;
    int rank;
} ds_node;

#define LE_OTHEREND(e,n) (((e)->ends[0] == (n)) ? ((e)->ends[1]) : ((e)->ends[0]))

typedef struct le_adj {
    struct le_edge *this;
    struct le_adj *next;
} le_adj;

typedef struct le_node {
    int number;
    int mark;
    struct le_node *next_member;
    struct le_edge *edgein;
    le_adj *adj;
} le_node;

typedef struct le_edge {
    int len;
    int live;
    double x;
    struct le_node *ends[2];
} le_edge;

typedef struct le_graph {
    int nnodes;
    int nedges;
    le_node *nodes;
    le_edge *edges;
    le_adj *adjs;
} le_graph;

typedef struct le_penalty {
    double delta;
    double delta_bad;
    double uppen_bad;
    double downpen_bad;
    double uppen;
    double downpen;
    double uppen2;
    double downpen2;
    double uppen3;
    double downpen3;
    double uppen4;
    double downpen4;
    int inlen;
    int outlen;
    int ncnt;
    int nodes[4];
} le_penalty;
    
    
typedef struct le_info {
    double maxup;
    double maxdown;
    int nwant;
    CCtsp_lpclique *clist;
    double *cval;
    le_penalty *penalties;
    int *workarr;
} le_info;

static int
    checkclique (CCtsp_lp *lp, CCtsp_lpclique *cliq, double delta, int *mark,
        double *downpen, double *uppen),
    CC_UNUSED checkclique2 (CCtsp_lp *lp, CCtsp_lpclique *cliq, double *x,
        double delta, int *mark, double *downpen, double *uppen),
    newclique (int nwant, CCtsp_lpclique *clist, double *cval,
        CCtsp_lpclique *c, double down, double up),
    newclique2 (int nwant, CCtsp_lpclique *clist, double *cval,
        int ncnt, le_node **nodes, double down, double up, int *workarr,
        le_penalty *pen, le_penalty *penalties),
    build_le_graph (le_graph *g, int nnodes, int nedges, int *elist,
        int *elen, double *x),
    le_findmax (int ncnt, le_node **nodes, le_info *le_dat),
    le_getcliques (int ncnt, le_node **nodes, le_info *le_dat),
    le_foreach_clique (le_graph *g, int (*func)(int ncnt,
        le_node **nodes, le_info *le_dat), le_info *le_dat),
    CC_UNUSED le_foreach_clique2 (le_graph *g, CCdatagroup *dat,
        int (*func)(int ncnt, le_node **nodes, le_info *le_dat),
        le_info *le_dat),
    merge_edge_clique (CCtsp_lp *lp, int nwant, int *ngot,
        CCtsp_branchobj **bobj, int ecount, int *elist, double *eval,
        int ccount, CCtsp_lpclique *clist, double *cval),
    find_strong_branch (CCtsp_lp *lp, int *n0, int *n1, int silent),
    find_strongbranch_edges (CCtsp_lp *lp, int nwant, int *ngot, int **elist,
        double **eval, int silent),
    find_candidate_edges (CCtsp_lp *lp, int nwant, int *ngot, int **list,
        int silent),
    find_all_candidate_edges (CCtsp_lp *lp, int *ngot, int **list),
    find_candidate_cliques (CCtsp_lp *lp, int nwant, int *ngot,
        CCtsp_lpclique **list, int use_getweight, int silent),
    find_branched_clique (CCtsp_lp *lp, CCtsp_lpclique *c, char sense, int rhs,
        int *cutnum),
    CC_UNUSED find_longedge_cliques (CCtsp_lp *lp, int nwant, int *ngot,
        CCtsp_lpclique **list, int silent),
    find_longedge_cliques2 (CCtsp_lp *lp, int nwant, int *ngot,
        CCtsp_lpclique **list, int silent),
    branch_side (CCtsp_lp *lp, CCtsp_branchobj *b, int side, int child,
        double *val, int *prune, int silent, CCrandstate *rstate),
    test_cut_branch (CCtsp_lp *lp, CCtsp_lpclique *c, double *down,
        double *up, int iter, int silent, CClp_warmstart *warmstart);

static void
    ds_makeset (ds_node *v),
    print_branchobj (CCtsp_branchobj *b),
    init_le_graph (le_graph *g),
    free_le_graph (le_graph *g),
    le_contract_work (le_node *m, le_node *newn, le_adj **p_newadj),
    le_contract_finish (le_node *newn, le_adj *adj),
    le_cliquevals (int ncnt, le_node **nodes, double *p_delta,
        int *p_inlen, int *p_outlen, le_penalty *pen),
    init_sblist (sbitem *list, int count),
    insert_sblist (sbitem *list, double val, int name),
    insert_sblist4 (sbitem *list, double val, int name, int name2, int name3,
                    int name4);

static le_node
   *le_contract (le_node *end0, le_node *end1);

static ds_node
   *ds_find (ds_node *v),
   *ds_link (ds_node *x, ds_node *y);


void CCtsp_init_branchobj (CCtsp_branchobj *b)
{
    b->depth     = 0;
    b->rhs       = 0;
    b->ends[0]   = -1;
    b->ends[1]   = -1;
    b->sense     = 'X';
    b->clique    = (CCtsp_lpclique *) NULL;
}

void CCtsp_free_branchobj (CCtsp_branchobj *b)
{
    if (!b) return;

    b->depth     = 0;
    b->rhs       = 0;
    b->ends[0]   = -1;
    b->ends[1]   = -1;
    b->sense     = 'X';
    if (b->clique) {
        CCtsp_free_lpclique (b->clique);
        CC_FREE (b->clique, CCtsp_lpclique);
    }
}

void CCtsp_print_branchhistory (CCtsp_lp *lp)
{
    int j;
    printf ("Branch History\n"); fflush (stdout);
    if (lp->branchdepth == 0) {
        printf ("    Root Node\n");
    } else {
        for (j = 0; j < lp->branchdepth; j++) {
            printf ("    ");
            print_branchobj (&lp->branchhistory[j]);
        }
    }
    fflush (stdout);
}

static void print_branchobj (CCtsp_branchobj *b)
{
    int i;

    printf ("Depth %d:  ", b->depth);
    if (b->ends[0] != -1) {
        printf ("Edge (%d,%d) set to %d\n", b->ends[0], b->ends[1], b->rhs);
    } else {
        printf ("Clique ");
        for (i = 0; i < b->clique->segcount; i++) {
            printf ("%d->%d ", b->clique->nodes[i].lo, b->clique->nodes[i].hi);
        }
        if (b->sense == 'L') {
            printf ("at most %d\n", b->rhs);
        } else {
            printf ("at least %d\n", b->rhs);
        }
    }
    fflush (stdout);
}

int CCtsp_find_branch (CCtsp_lp *lp, int nwant, int *ngot,
        CCtsp_branchobj **bobj, double *val, int **cyc, int usecliques,
        int longedge_branching, int silent)
{
    int rval = 0;
    int egot = 0;
    int cgot = 0;
    int      *elist = (int *) NULL;
    CCtsp_lpclique *clist = (CCtsp_lpclique *) NULL;
    double   *eval  = (double *) NULL;
    double   *cval  = (double *) NULL;
    int i, n0, n1;

    CCutil_start_timer (&lp->stats.strongbranch);

    *ngot = 0;
    *bobj = (CCtsp_branchobj *) NULL;
    if (cyc) *cyc = (int *) NULL;

    if (nwant <= 0) {
        fprintf (stderr, "CCtsp_find_branch called with no nwant\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_find_branch_edge (lp, &n0, &n1, val, cyc,
                                   CCtsp_BRANCH_MIDDLE, silent);
    if (rval) {
        fprintf (stderr, "CCtsp_find_branch_edge failed\n"); goto CLEANUP;
    }

    if (n0 == -1 && n1 == -1) {
        if (!silent) {
            printf ("Integral solution: %f\n", *val); fflush (stdout);
        }
        goto CLEANUP;
    }

    rval = find_strongbranch_edges (lp, nwant, &egot, &elist, &eval, silent);
    if (rval) {
        fprintf (stderr, "find_strongbranch_edges failed\n"); goto CLEANUP;
    }

    if (longedge_branching || usecliques) {
        rval = CCtsp_find_branch_cliques (lp, nwant, longedge_branching,
                                          &cgot, &clist, &cval, silent);
        if (rval) {
            fprintf (stderr, "CCtsp_find_branch_cliques failed\n");
            goto CLEANUP;
        }
        if (!silent) {
            printf ("Cliques found:\n"); fflush (stdout);
            for (i = 0; i < cgot; i++) {
                CCtsp_print_lpclique (&clist[i]);
            }
            fflush (stdout);
        }
    }

    if (egot + cgot > 0) {
        rval = merge_edge_clique (lp, nwant, ngot, bobj, egot, elist, eval,
                                                         cgot, clist, cval);
        if (rval) {
            fprintf (stderr, "merge_edge_clique failed\n"); goto CLEANUP;
        }
    } else {
        CCtsp_branchobj *b;
        if (!silent) {
            printf ("found no edges or cliques, use the middle branch edge\n");
            fflush (stdout);
        }

        b = CC_SAFE_MALLOC (1, CCtsp_branchobj);
        if (!b) {
            fprintf (stderr, "out of memory in CCtsp_find_branch\n");
            rval = 1; goto CLEANUP;
        }
        CCtsp_init_branchobj (b);
        b->ends[0] = n0;
        b->ends[1] = n1;
        *bobj = b;
        *ngot = 1;
    }

CLEANUP:

    if (!silent) {
        CCutil_stop_timer (&lp->stats.strongbranch, 1);
    } else {
        CCutil_stop_timer (&lp->stats.strongbranch, 0);
    }

    CC_IFFREE (elist, int);
    for (i = 0; i < cgot; i++) {
        CCtsp_free_lpclique (&clist[i]);
    }
    CC_IFFREE (clist, CCtsp_lpclique);
    CC_IFFREE (cval, double);
    CC_IFFREE (eval, double);
    return rval;
}

int CCtsp_find_fast_branch (CCtsp_lp *lp, int *ngot, CCtsp_branchobj **bobj,
        double *val, int **cyc, int usecliques, int longedge_branching,
        int silent)
{
    int rval;
    int n0, n1, i, besti;
    int ccount = 0;
    double sval, bestval, up, down;
    CCtsp_lpclique *cliques = (CCtsp_lpclique *) NULL;
    CCtsp_branchobj *b;
    CCutil_timer timer;
    CClp_warmstart *warmstart = (CClp_warmstart *) NULL;

    CCutil_init_timer (&timer, "Strong fast branch");

    *ngot = 0;

    rval = CCtsp_find_branch_edge (lp, &n0, &n1, val, cyc,
                                   CCtsp_BRANCH_MIDDLE, silent);
    if (rval) {
        fprintf (stderr, "CCtsp_find_branch_edge failed\n"); goto CLEANUP;
    }

    if (n0 == -1 && n1 == -1) {
        if (!silent) {
            printf ("Integral solution: %f\n", *val); fflush (stdout);
        }
        goto CLEANUP;
    }

    if (longedge_branching || usecliques) {
        if (longedge_branching) {
            rval = find_longedge_cliques2 (lp, 10, &ccount, &cliques, silent);
            if (rval) {
                fprintf (stderr, "find_longedge_cliques failed\n");
                goto CLEANUP;
            }
        } else {
            rval = find_candidate_cliques (lp, 10, &ccount, &cliques, 1, silent);
            if (rval) {
                fprintf (stderr, "find_candidate_cliques failed\n"); goto CLEANUP;
            }
        }
        if (!silent) {
            printf ("Found %d candidate cliques\n", ccount); fflush (stdout);
        }
    }

    b = CC_SAFE_MALLOC (1, CCtsp_branchobj);
    if (!b) {
        fprintf (stderr, "out of memory in CCtsp_find_branch\n");
        rval = 1; goto CLEANUP;
    }
    CCtsp_init_branchobj (b);

    besti = -1;
    bestval = -CCutil_MAXDOUBLE;

    if (ccount) {
        rval = CClp_get_warmstart (lp->lp, &warmstart);
        if (rval) {
            fprintf (stderr, "CClp_get_warmstart failed\n"); goto CLEANUP;
        }
        for (i = 0; i < ccount; i++) {
            double st;

            CCutil_start_timer (&timer);
            rval = test_cut_branch (lp, &cliques[i], &down, &up, 25, silent,
                                    warmstart);
            if (rval) {
                fprintf (stderr, "test_cut_branch failed\n");
                goto CLEANUP;
            }
            st = CCutil_stop_timer (&timer, 0);
            if (!silent) {
                printf ("SB CLIQUE %d:  %f  %f  (%.2f seconds)\n",
                           i, down, up, st);
                fflush (stdout);
            }
            sval = TSP_BRANCH_STRONG_VAL (down, up);
            if (sval > bestval) {
                bestval = sval;
                besti = i;
            }
        }

        rval = CClp_opt (lp->lp, CClp_METHOD_DUAL);
        if (rval) {
            fprintf (stderr, "CClp_opt failed\n");
        }

        if (besti == -1) {
            fprintf (stderr, "ERROR in CCtsp_find_fast_branch\n");
            rval = 1; goto CLEANUP;
        } else {
            b->clique = CC_SAFE_MALLOC (1, CCtsp_lpclique);
            if (!b->clique) {
                fprintf (stderr, "out of memory in merge_edge_clique\n");
                rval = 1; goto CLEANUP;
            }
            rval = CCtsp_copy_lpclique (&(cliques[besti]), b->clique);
            if (rval) {
                fprintf (stderr, "CCtsp_copy_clique failed\n");
                CC_FREE (b, CCtsp_branchobj);
                goto CLEANUP;
            }
        }
    } else {
        b->ends[0] = n0;
        b->ends[1] = n1;
    }

    *bobj = b;
    *ngot = 1;

CLEANUP:

    if (cliques) {
        for (i = 0; i < ccount; i++) {
            CCtsp_free_lpclique (&(cliques[i]));
        }
        CC_FREE (cliques, CCtsp_lpclique);
    }
    CClp_free_warmstart (&warmstart);

    return rval;
}

/* disjoint sets ala Tarjan (from Data Structures and Network Algorithms
   by Robert Tarjan) */

static void ds_makeset (ds_node *v)
{
    v->parent = v;
    v->rank = 0;
}

static ds_node *ds_find (ds_node *v)
{
    ds_node *p = v->parent;

    return v == p ? v : (v->parent = ds_find (p));
}

static ds_node *ds_link (ds_node *x, ds_node *y)
{
    ds_node *t;

    x = ds_find (x);
    y = ds_find (y);

    if (x != y) {
        if (x->rank > y->rank) {
            CC_SWAP (x,y,t);
        } else if (x->rank == y->rank) {
            y->rank++;
        }
        x->parent = y;
    }
    return y;
}

static int checkclique (CCtsp_lp *lp, CCtsp_lpclique *cliq, double delta,
        int *mark, double *downpen, double *uppen)
{
    int i;
    int j;
    int k;
    int m;
    int e = 0;
    CCtsp_lpnode *nodes = lp->graph.nodes;
    CCtsp_lpedge *edges = lp->graph.edges;
    int noutside = 0;
    int ninside = 0;
    int ncomp = 0;
    int *inside = (int *) NULL;
    int *outside = (int *) NULL;
    int *inlen = (int *) NULL;
    int *outlen = (int *) NULL;
    int *inperm = (int *) NULL;
    int *outperm = (int *) NULL;
    ds_node *innodes = (ds_node *) NULL;
    int rval;

    CC_FOREACH_NODE_IN_CLIQUE (i, *cliq, j) {
        mark[i] = 1;
    }
    CC_FOREACH_NODE_IN_CLIQUE (i, *cliq, j) {
        for (k=0; k<nodes[i].deg; k++) {
            m = nodes[i].adj[k].to;
            if (mark[m] == 0) {
                noutside++;
            } else if (i < m) {
                ninside++;
            }
        }
    }
    inside = CC_SAFE_MALLOC (ninside, int);
    outside = CC_SAFE_MALLOC (noutside, int);
    inlen = CC_SAFE_MALLOC (ninside, int);
    outlen = CC_SAFE_MALLOC (noutside, int);
    inperm = CC_SAFE_MALLOC (ninside, int);
    outperm = CC_SAFE_MALLOC (noutside, int);
    innodes = CC_SAFE_MALLOC (lp->graph.ncount, ds_node);
    if (inside == (int *) NULL || outside == (int *) NULL ||
        inlen == (int *) NULL || outlen == (int *) NULL ||
        inperm == (int *) NULL || outperm == (int *) NULL ||
        innodes == (ds_node *) NULL) {
        fprintf (stderr, "Out of memory in checkclique\n");
        rval = 1; goto CLEANUP;
    }

    ninside = 0;
    noutside = 0;
    CC_FOREACH_NODE_IN_CLIQUE (i, *cliq, j) {
        ds_makeset (&innodes[i]);
        ncomp++;
        for (k=0; k<nodes[i].deg; k++) {
            m = nodes[i].adj[k].to;
            if (mark[m] == 0) {
                e = nodes[i].adj[k].edge;
                outside[noutside] = e;
                outlen[noutside] = edges[e].len;
                outperm[noutside] = noutside;
                noutside++;
            } else if (i < m) {
                e = nodes[i].adj[k].edge;
                inside[ninside] = e;
                inlen[ninside] = edges[e].len;
                inperm[ninside] = ninside;
                ninside++;
            }
        }
    }

    CC_FOREACH_NODE_IN_CLIQUE (i, *cliq, j) {
        mark[i] = 0;
    }

    CCutil_int_perm_quicksort (inperm, inlen, ninside);
    CCutil_int_perm_quicksort (outperm, outlen, noutside);
    
    if (uppen != (double *) NULL) {
        *uppen = ((4.0 - delta) * outlen[outperm[0]]);
    }

    for (i=0; i<ninside && ncomp > 1; i++) {
        e = inside[inperm[i]];
        if (ds_find (&innodes[edges[e].ends[0]]) !=
            ds_find (&innodes[edges[e].ends[1]])) {
            ncomp--;
            ds_link (&innodes[edges[e].ends[0]], &innodes[edges[e].ends[1]]);
        }
    }

    if (downpen != (double *) NULL) {
        *downpen = (((delta - 2.0) / 2.0) * edges[e].len);
    }
#ifdef DEBUG
    printf ("clique ");
    CCtsp_print_lpclique (cliq);
    printf ("delta %.6f outlen %d inlen %d uppen %.6f downpen %.6f\n",
            delta, outlen[outperm[0]], edges[e].len, *uppen, *downpen);
#endif /* DEBUG */
    
    rval = 0;

 CLEANUP:
    CC_IFFREE (inside, int);
    CC_IFFREE (outside, int);
    CC_IFFREE (inlen, int);
    CC_IFFREE (outlen, int);
    CC_IFFREE (inperm, int);
    CC_IFFREE (outperm, int);
    CC_IFFREE (innodes, ds_node);

    return rval;
}

static int CC_UNUSED checkclique2 (CCtsp_lp *lp, CCtsp_lpclique *cliq,
        double *x, double delta, int *mark, double *downpen, double *uppen)
{
    int i;
    int j;
    int k;
    int m;
    int e = 0;
    CCtsp_lpnode *nodes = lp->graph.nodes;
    CCtsp_lpedge *edges = lp->graph.edges;
    int noutside = 0;
    int ninside = 0;
    int ncomp = 0;
    int *inside = (int *) NULL;
    int *outside = (int *) NULL;
    double *inlen = (double *) NULL;
    double *outlen = (double *) NULL;
    int *inperm = (int *) NULL;
    int *outperm = (int *) NULL;
    ds_node *innodes = (ds_node *) NULL;
    double sum;
    double cost;
    int rval;

    CC_FOREACH_NODE_IN_CLIQUE (i, *cliq, j) {
        mark[i] = 1;
    }
    CC_FOREACH_NODE_IN_CLIQUE (i, *cliq, j) {
        for (k=0; k<nodes[i].deg; k++) {
            m = nodes[i].adj[k].to;
            if (mark[m] == 0) {
                noutside++;
            } else if (i < m) {
                ninside++;
            }
        }
    }
    inside = CC_SAFE_MALLOC (ninside, int);
    outside = CC_SAFE_MALLOC (noutside, int);
    inlen = CC_SAFE_MALLOC (ninside, double);
    outlen = CC_SAFE_MALLOC (noutside, double);
    inperm = CC_SAFE_MALLOC (ninside, int);
    outperm = CC_SAFE_MALLOC (noutside, int);
    innodes = CC_SAFE_MALLOC (lp->graph.ncount, ds_node);
    if (inside == (int *) NULL || outside == (int *) NULL ||
        inlen == (double *) NULL || outlen == (double *) NULL ||
        inperm == (int *) NULL || outperm == (int *) NULL ||
        innodes == (ds_node *) NULL) {
        fprintf (stderr, "Out of memory in checkclique2\n");
        rval = 1; goto CLEANUP;
    }

    ninside = 0;
    noutside = 0;
    CC_FOREACH_NODE_IN_CLIQUE (i, *cliq, j) {
        ds_makeset (&innodes[i]);
        ncomp++;
        for (k=0; k<nodes[i].deg; k++) {
            m = nodes[i].adj[k].to;
            if (mark[m] == 0) {
                e = nodes[i].adj[k].edge;
                outside[noutside] = e;
                outlen[noutside] = edges[e].len * (1.0 - x[e]);
                outperm[noutside] = noutside;
                noutside++;
            } else if (i < m) {
                e = nodes[i].adj[k].edge;
                inside[ninside] = e;
                inlen[ninside] = edges[e].len * (1.0 - x[e]);
                inperm[ninside] = ninside;
                ninside++;
            }
        }
    }

    CC_FOREACH_NODE_IN_CLIQUE (i, *cliq, j) {
        mark[i] = 0;
    }

    CCutil_double_perm_quicksort (inperm, inlen, ninside);
    CCutil_double_perm_quicksort (outperm, outlen, noutside);

    sum = delta;

    cost = 0.0;
    for (i=0; i<noutside && sum < 4.0; i++) {
        e = outside[outperm[i]];
        if (sum + (1.0 - x[e]) <= 4.0) {
            cost += (1.0 - x[e]) * edges[e].len;
            sum += 1.0 - x[e];
        } else {
            cost += (4.0 - sum) * edges[e].len;
            sum += 4.0 - sum;
        }
    }

    if (uppen != (double *) NULL) {
        *uppen = cost;
    }

    for (i=0; i<ninside && ncomp > 1; i++) {
        e = inside[inperm[i]];
        if (ds_find (&innodes[edges[e].ends[0]]) !=
            ds_find (&innodes[edges[e].ends[1]])) {
            ncomp--;
            ds_link (&innodes[edges[e].ends[0]], &innodes[edges[e].ends[1]]);
        }
    }

    if (downpen != (double *) NULL) {
        *downpen = ((1.0 - x[e]) * edges[e].len);
    }
    
    rval = 0;

 CLEANUP:
    CC_IFFREE (inside, int);
    CC_IFFREE (outside, int);
    CC_IFFREE (inlen, double);
    CC_IFFREE (outlen, double);
    CC_IFFREE (inperm, int);
    CC_IFFREE (outperm, int);
    CC_IFFREE (innodes, ds_node);

    return rval;
}

static int newclique (int nwant, CCtsp_lpclique *clist, double *cval,
        CCtsp_lpclique *c, double down, double up)
{
    double angle = atan2(down, up);
    int bin = (int) (nwant * angle * 2.0 / M_PI);
    double val = down * down + up * up;
    int rval;

    if (bin == nwant) bin = nwant-1;

    if (val > cval[bin]) {
        cval[bin] = val;
        CCtsp_free_lpclique (&clist[bin]);
        rval = CCtsp_copy_lpclique (c, &clist[bin]);
        if (rval) {
            fprintf (stderr, "CCtsp_copy_lpclique failed\n");
            return rval;
        }
    }
    return 0;
}

static int CC_UNUSED find_longedge_cliques (CCtsp_lp *lp, int nwant, int *ngot,
        CCtsp_lpclique **list, int silent)
{
    int ncount = lp->graph.ncount;
    CCdatagroup *dat = lp->dat;
    int rval = 0;
    int i, j, k, l;
    int nbig;
    int *big = (int *) NULL;
    double    delta = 0.0;
    int seglist[4];
    CCtsp_lpclique cliq;
    CCtsp_lpclique *clist = (CCtsp_lpclique *) NULL;
    double *cval = (double *) NULL;
    double *x = (double *) NULL;
    int *mark = (int *) NULL;
    double *len = (double *) NULL;
    int *perm = (int *) NULL;
    int xcount;
    double downpen;
    double uppen;
    double maxdown;
    double maxup;
    double st;
    CCutil_timer timer;

    *ngot = 0;
    CCtsp_init_lpclique (&cliq);
    CCutil_init_timer (&timer, "Long-edge search");
    CCutil_start_timer (&timer);

    if (!silent) {
        printf ("Finding long-edge cliques\n"); fflush (stdout);
    }
    
    if (nwant <= 0) {
        fprintf (stderr, "find_longedge_cliques called with no nwant\n");
        rval = 1; goto CLEANUP;
    }

    if (dat == (CCdatagroup *) NULL) {
        fprintf (stderr, "find_longedge_cliques called with no dat\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL, &xcount,
                                (int **) NULL, &x, (double **) NULL,
                                (double **) NULL, (double **) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_get_lp_result failed\n");
        goto CLEANUP;
    }

    if (xcount != lp->graph.ecount) {
        fprintf (stderr, "Whoops, xcount %d ecount %d\n", xcount,
                 lp->graph.ecount);
        rval = 1; goto CLEANUP;
    }

    mark = CC_SAFE_MALLOC (ncount, int);
    if (mark == (int *) NULL) {
        fprintf (stderr, "Out of memory in find_longedge_cliques\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<ncount; i++) {
        mark[i] = 0;
    }

    len = CC_SAFE_MALLOC (ncount, double);
    perm = CC_SAFE_MALLOC (ncount, int);
    if (len == (double *) NULL ||
        perm == (int *) NULL) {
        fprintf (stderr, "Out of memory in find_longedge_cliques\n");
        rval = 1; goto CLEANUP;
    }
    for (i=0; i<ncount-1; i++) {
        len[i] = CCutil_dat_edgelen (i, i+1, dat);
        perm[i] = i;
    }
    len[ncount-1] = CCutil_dat_edgelen (ncount-1, 0, dat);
    perm[ncount-1] = ncount-1;

    nbig = (int) (sqrt((double) ncount)+0.5);

    CCutil_double_perm_quicksort (perm, len, ncount);
    CCutil_int_array_quicksort (perm + ncount - nbig, nbig);
    
    big = CC_SAFE_MALLOC (nbig, int);
    if (big == (int *) NULL) {
        fprintf (stderr, "Out of memory in find_longedge_cliques\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<nbig; i++) {
        big[i] = perm[ncount - nbig + i];
    }

    CC_IFFREE (perm, int);
    CC_IFFREE (len, double);

    maxup = 0.0;
    maxdown = 0.0;
    for (i=0; i<nbig; i++) {
        for (j=i+1; j<nbig; j++) {
            seglist[0] = big[i]+1;
            seglist[1] = big[j];
            rval = CCtsp_seglist_to_lpclique (1, seglist, &cliq);
            if (rval) {
                fprintf (stderr, "CCtsp_seglist_to_lpclique failed\n");
                goto CLEANUP;
            }
            rval = CCtsp_clique_delta (&lp->graph, x, &cliq, &delta);
            if (rval) {
                fprintf (stderr, "CCtsp_clique_delta failed\n");
                goto CLEANUP;
            }
            if (delta > (2.0 + TSP_BRANCH_CLIQUE_TOL) &&
                delta < (4.0 - TSP_BRANCH_CLIQUE_TOL)) {
                rval = checkclique (lp, &cliq, delta, mark, &downpen, &uppen);
                if (rval) {
                    fprintf (stderr, "checkclique failed\n");
                    goto CLEANUP;
                }
                if (downpen > maxdown) maxdown = downpen;
                if (uppen > maxup) maxup = uppen;
            }
            CCtsp_free_lpclique (&cliq);
            for (k=j+1; k<nbig; k++) {
                for (l=k+1; l<nbig; l++) {
                    seglist[2] = big[k]+1;
                    seglist[3] = big[l];
                    rval = CCtsp_seglist_to_lpclique (2, seglist, &cliq);
                    if (rval) {
                        fprintf (stderr, "CCtsp_seglist_to_lpclique failed\n");
                        goto CLEANUP;
                    }
                    rval = CCtsp_clique_delta (&lp->graph, x, &cliq, &delta);
                    if (rval) {
                        fprintf (stderr, "CCtsp_clique_delta failed\n");
                        goto CLEANUP;
                    }
                    if (delta > (2.0 + TSP_BRANCH_CLIQUE_TOL) &&
                        delta < (4.0 - TSP_BRANCH_CLIQUE_TOL)) {
                        rval = checkclique (lp, &cliq, delta, mark, &downpen,
                                            &uppen);
                        if (rval) {
                            fprintf (stderr, "checkclique failed\n");
                            goto CLEANUP;
                        }
                        if (downpen > maxdown) maxdown = downpen;
                        if (uppen > maxup) maxup = uppen;
                    }
                    CCtsp_free_lpclique (&cliq);
                }
            }
        }
    }

    if (maxdown == 0.0) maxdown = 1.0;
    if (maxup == 0.0) maxup = 1.0;
/****/
    maxup = 3.0;
    maxdown = 5.0;

    clist = CC_SAFE_MALLOC (nwant, CCtsp_lpclique);
    cval = CC_SAFE_MALLOC (nwant, double);
    if (clist == (CCtsp_lpclique *) NULL ||
        cval == (double *) NULL) {
        fprintf (stderr, "Out of memory in find_longedge_cliques\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<nwant; i++) {
        cval[i] = -1.0;
        CCtsp_init_lpclique (&clist[i]);
    }
    
    for (i=0; i<nbig; i++) {
        for (j=i+1; j<nbig; j++) {
            seglist[0] = big[i]+1;
            seglist[1] = big[j];
            rval = CCtsp_seglist_to_lpclique (1, seglist, &cliq);
            if (rval) {
                fprintf (stderr, "CCtsp_seglist_to_lpclique failed\n");
                goto CLEANUP;
            }
            rval = CCtsp_clique_delta (&lp->graph, x, &cliq, &delta);
            if (rval) {
                fprintf (stderr, "CCtsp_clique_delta failed\n");
                goto CLEANUP;
            }
            if (delta > (2.0 + TSP_BRANCH_CLIQUE_TOL) &&
                delta < (4.0 - TSP_BRANCH_CLIQUE_TOL)) {
                rval = checkclique (lp, &cliq, delta, mark, &downpen, &uppen);
                if (rval) {
                    fprintf (stderr, "checkclique failed\n");
                    goto CLEANUP;
                }
                rval = newclique (nwant, clist, cval, &cliq,
                                  downpen / maxdown, uppen / maxup);
                if (rval) {
                    fprintf (stderr, "newclique failed\n");
                    goto CLEANUP;
                }
            }
            CCtsp_free_lpclique (&cliq);
            for (k=j+1; k<nbig; k++) {
                for (l=k+1; l<nbig; l++) {
                    seglist[2] = big[k]+1;
                    seglist[3] = big[l];
                    rval = CCtsp_seglist_to_lpclique (2, seglist, &cliq);
                    if (rval) {
                        fprintf (stderr, "CCtsp_seglist_to_lpclique failed\n");
                        goto CLEANUP;
                    }
                    rval = CCtsp_clique_delta (&lp->graph, x, &cliq, &delta);
                    if (rval) {
                        fprintf (stderr, "CCtsp_clique_delta failed\n");
                        goto CLEANUP;
                    }
                    if (delta > (2.0 + TSP_BRANCH_CLIQUE_TOL) &&
                        delta < (4.0 - TSP_BRANCH_CLIQUE_TOL)) {
                        rval = checkclique (lp, &cliq, delta, mark, &downpen,
                                            &uppen);
                        if (rval) {
                            fprintf (stderr, "checkclique failed\n");
                            goto CLEANUP;
                        }
                        rval = newclique (nwant, clist, cval, &cliq,
                                          downpen / maxdown, uppen / maxup);
                        if (rval) {
                            fprintf (stderr, "newclique failed\n");
                            goto CLEANUP;
                        }
                    }
                    CCtsp_free_lpclique (&cliq);
                }
            }
        }
    }

    k = 0;
    for (i = 0, k=0; i<nwant; i++) {
        if (cval[i] >= 0.0) {
            k++;
        }
    }

    if (k == 0) {
        rval = 0;
        goto CLEANUP;
    }

    *list = CC_SAFE_MALLOC (k, CCtsp_lpclique);
    if (!(*list)) {
        fprintf (stderr, "Out of memory in find_longedge_cliques\n");
        rval = 1; goto CLEANUP;
    }
    *ngot = k;

    for (i=0, k=0; i<nwant; i++) {
        if (cval[i] >= 0.0) {
            if (!silent) {
                printf ("%d: (%.6f) ", i, cval[i]);
                CCtsp_print_lpclique (&clist[i]);
            }
            (*list)[k] = clist[i];
            k++;
        }
    }

    st = CCutil_stop_timer (&timer, 0);
    if (!silent) {
        printf ("%d long-edge cliques found in %.2f seconds\n", *ngot, st);
        fflush (stdout);
    }
    
    CC_IFFREE (clist, CCtsp_lpclique);
    rval = 0;

CLEANUP:
    if (clist) {
        for (i=0; i<nwant; i++) {
            CCtsp_free_lpclique (&clist[i]);
        }
    }
    CCtsp_free_lpclique (&cliq);
    CC_IFFREE (clist, CCtsp_lpclique);
    CC_IFFREE (cval, double);
    CC_IFFREE (big, int);
    CC_IFFREE (x, double);
    CC_IFFREE (mark, int);
    CC_IFFREE (len, double);
    CC_IFFREE (perm, int);

    return rval;
}

static void init_le_graph (le_graph *g)
{
    g->nnodes = 0;
    g->nedges = 0;
    g->nodes = (le_node *) NULL;
    g->edges = (le_edge *) NULL;
    g->adjs = (le_adj *) NULL;
}

static void free_le_graph (le_graph *g)
{
    CC_IFFREE (g->nodes, le_node);
    CC_IFFREE (g->edges, le_edge);
    CC_IFFREE (g->adjs, le_adj);
}

static int build_le_graph (le_graph *g, int nnodes, int nedges, int *elist,
        int *elen, double *x)
{
    int rval;
    int i;
    int j;
    le_node *n;
    le_node *nodes;
    le_edge *edges;
    int *perm = (int *) NULL;
    
    init_le_graph (g);
    g->nnodes = nnodes;
    g->nedges = nedges;
    g->nodes = CC_SAFE_MALLOC (nnodes, le_node);
    g->edges = CC_SAFE_MALLOC (nedges, le_edge);
    g->adjs  = CC_SAFE_MALLOC (nedges*2, le_adj);
    perm     = CC_SAFE_MALLOC (nedges, int);
    if (g->nodes == (le_node *) NULL ||
        g->edges == (le_edge *) NULL ||
        g->adjs  == (le_adj *)  NULL ||
        perm     == (int *)     NULL) {
        fprintf (stderr, "Out of memory in build_le_graph\n");
        rval = 1; goto CLEANUP;
    }

    nodes = g->nodes;
    edges = g->edges;
    for (i=0; i<nnodes; i++) {
        nodes[i].number = i;
        nodes[i].next_member = &(nodes[i]);
        nodes[i].mark = 0;
        nodes[i].adj = (le_adj *) NULL;
        nodes[i].edgein = (le_edge *) NULL;
    }

    for (i=0; i<nedges; i++) {
        perm[i] = i;
    }
    
    CCutil_int_perm_quicksort (perm, elen, nedges);
    
    for (i=0; i<nedges; i++) {
        j = perm[i];
        edges[i].len     = elen[j];
        edges[i].x       = x[j];
        edges[i].live    = 1;
        edges[i].ends[0] = &(nodes[elist[2*j]]);
        edges[i].ends[1] = &(nodes[elist[2*j+1]]);
        edges[i].ends[0]->mark++;
        edges[i].ends[1]->mark++;
    }
    j=0;
    for (i=0; i<nnodes; i++) {
        nodes[i].adj = g->adjs + j;
        j += nodes[i].mark;
        nodes[i].mark = 0;
    }
    for (i=0; i<nedges; i++) {
        n = edges[i].ends[0];
        n->adj[n->mark].this = &edges[i];
        n->adj[n->mark].next = &(n->adj[n->mark+1]);
        n->mark++;
        n = edges[i].ends[1];
        n->adj[n->mark].this = &edges[i];
        n->adj[n->mark].next = &(n->adj[n->mark+1]);
        n->mark++;
    }
    for (i=0; i<nnodes; i++) {
        nodes[i].adj[nodes[i].mark-1].next = (le_adj *) NULL;
        nodes[i].mark = 0;
    }

#if 0
    {
        static int n = 0;
        int i;
        char fname[50];
        FILE *f;
        sprintf (fname, "legraph.%d.out", n);
        n++;
        f = fopen (fname, "w");
        fprintf (f, "%d %d\n", g->nnodes, g->nedges);
        for (i=0; i<g->nedges; i++) {
            fprintf (f, "%d %d %.6f %d\n", g->edges[i].ends[0]->number,
                     g->edges[i].ends[1]->number, g->edges[i].x,
                     g->edges[i].len);
        }
        fclose (f);
    }
#endif
    
    rval = 0;
    
 CLEANUP:
    if (rval) {
        free_le_graph (g);
    }
    CC_IFFREE (perm, int);
    return rval;
}

static void le_contract_work (le_node *m, le_node *newn, le_adj **p_newadj)
{
    le_adj *a;
    le_adj *anext;
    le_edge *e;
    le_node *n;
    
    le_adj *newadj = *p_newadj;
    
    for (a = m->adj; a; a = anext) {
        e = a->this;
        anext = a->next;
        if (e->live) {
            if (e->ends[0] == m) {
                e->ends[0] = newn;
                n = e->ends[1];
            } else {
                e->ends[1] = newn;
                n = e->ends[0];
            }
            if (n->edgein == (le_edge *) NULL) {
                a->next = newadj;
                newadj = a;
                n->edgein = e;
            } else {
                if (e->len < n->edgein->len) {
                    n->edgein->live = 0;
                    e->x += n->edgein->x;
                    n->edgein = e;
                } else {
                    e->live = 0;
                    n->edgein->x += e->x;
                }
            }
        }
    }

    *p_newadj = newadj;
}

static void le_contract_finish (le_node *newn, le_adj *adj)
{
    le_adj *a;
    le_adj *anext;
    le_node *n;
    le_edge *e;
    le_adj *newadj = (le_adj *) NULL;

    for (a = adj; a; a = anext) {
        e = a->this;
        anext = a->next;
        a->next = newadj;
        newadj = a;
        n = LE_OTHEREND(e, newn);
        e = n->edgein;
        a->this = e;
        n->edgein = (le_edge *) NULL;
    }
    newn->adj = newadj;
}
        
static le_node *le_contract (le_node *end0, le_node *end1)
{
    le_adj *newadj;
    le_node *n;
    le_edge dummy;

    dummy.x = 0.0;
    dummy.len = 0;
    end0->edgein = &dummy;
    end1->edgein = &dummy;

    newadj = (le_adj *) NULL;

    le_contract_work (end0, end0, &newadj);
    le_contract_work (end1, end0, &newadj);
    le_contract_finish (end0, newadj);

    n = end0->next_member;
    end0->next_member = end1->next_member;
    end1->next_member = n;

    end0->edgein = (le_edge *) NULL;

    return end0;
}

static int newclique2 (int nwant, CCtsp_lpclique *clist, double *cval,
        int ncnt, le_node **nodes, double down, double up, int *workarr,
        le_penalty *pen, le_penalty *penalties)
{
    double angle = atan2(down, up);
    int bin = (int) (nwant * angle * 2.0 / M_PI);
    double val = down * down + up * up;
    int rval;
    le_node *n;
    int cnt;
    int i;

    if (bin == nwant) bin = nwant-1;

    if (val > cval[bin]) {
        cval[bin] = val;
        CCtsp_free_lpclique (&clist[bin]);
        cnt = 0;
        for (i=0; i<ncnt; i++) {
            n = nodes[i];
            do {
                workarr[cnt++] = n->number;
                n = n->next_member;
            } while (n != nodes[i]);
        }

        rval = CCtsp_array_to_lpclique (workarr, cnt, &clist[bin]);
        if (rval) {
            fprintf (stderr, "CCtsp_array_to_lpclique failed\n");
            return rval;
        }
        penalties[bin] = *pen;
    }
    return 0;
}

static void le_cliquevals (int ncnt, le_node **nodes, double *p_delta,
        int *p_inlen, int *p_outlen, le_penalty *pen)
{
    double delta = 0.0;
    double delta_bad = 0.0;
    double indelta;
    double outdelta;
    int inlen = 0;
    int outlen = 0;
    int inlen_bad = 0;
    int outlen_bad = 0;
    int minlen_bad = 0;
    int minlen = 0;
    int moutlen = 0;
    int gotin = 0;
    int gotout = 0;
    int gotmin;
    int gotin_bad = 0;
    int gotout_bad = 0;
    int gotmin_bad;
    int gotmout;
    int i;
    le_node *n;
    le_adj *a;
    le_edge *e;
    le_node *m;
    double bestdiff = 0.0;
    double bestdiff2 = 0.0;
    double sumtwo = 0.0;
    double summove = 0.0;
    double bestmove = 0.0;
    double bestmove2 = 0.0;
    double onemove;
    double twomove;
    double diff;
    double diffmove;
    int gotuc = 0;
    int gotuc2 = 0;
    double bestucost2 = 0.0;
    double bestucost = 0.0;
    double bestumove2 = 0.0;
    double bestumove = 0.0;
    double ucost;
    double umove;

    for (i=0; i<ncnt; i++) {
        n = nodes[i];
        n->mark = 1;
    }
    
    for (i=0; i<ncnt; i++) {
        n = nodes[i];
        gotmin_bad = 0;
        for (a = n->adj; a; a = a->next) {
            e = a->this;
            if (e->live) {
                m = LE_OTHEREND (e, n);
                if (m->mark) {
                    if (!gotmin_bad || e->len < minlen_bad) {
                        minlen_bad = e->len;
                        gotmin_bad = 1;
                    }
                } else {
                    delta_bad += e->x;
                    if (!gotout_bad || e->len < outlen_bad) {
                        outlen_bad = e->len;
                        gotout_bad = 1;
                    }
                }
            }
        }
        if (!gotin_bad || minlen_bad > inlen_bad) {
            inlen_bad = minlen_bad;
            gotin_bad = 1;
        }
    }

    for (i=0; i<ncnt; i++) {
        n = nodes[i];
        gotmin = 0;
        gotmout = 0;
        indelta = 0.0;
        outdelta = 0.0;
        for (a = n->adj; a; a = a->next) {
            e = a->this;
            if (e->live) {
                m = LE_OTHEREND (e, n);
                if (m->mark) {
                    indelta += e->x;
                    if (!gotmin || e->len < minlen) {
                        minlen = e->len;
                        gotmin = 1;
                    }
                } else {
                    delta += e->x;
                    outdelta += e->x;
                    if (!gotout || e->len < outlen) {
                        outlen = e->len;
                        gotout = 1;
                    }
                    if (!gotmout || e->len < moutlen) {
                        moutlen = e->len;
                        gotmout = 1;
                    }
                }
            }
        }
        if (!gotin || minlen > inlen) {
            inlen = minlen;
            gotin = 1;
        }
        if (indelta >= 2.0) twomove = 0.0;
        else                twomove = (2.0 - indelta);
        if (indelta >= 1.0) onemove = 0.0;
        else                onemove = (1.0 - indelta);
        sumtwo += twomove * minlen;
        summove += twomove;
        diffmove = (twomove - onemove);
        diff = diffmove * minlen;
        if (diff > bestdiff2) {
            if (diff > bestdiff) {
                bestdiff2 = bestdiff;
                bestmove2 = bestmove;
                bestdiff = diff;
                bestmove = diffmove;
            } else {
                bestdiff2 = diff;
                bestmove2 = diffmove;
            }
        }
        if (outdelta >= 2.0) umove = 0.0;
        else                 umove = (2.0 - outdelta);
        ucost = umove * moutlen;
        if (!gotuc2 || ucost < bestucost2) {
            if (!gotuc || ucost < bestucost) {
                bestucost2 = bestucost;
                bestumove2 = bestumove;
                bestucost = ucost;
                bestumove = umove;
                gotuc = 1;
            } else {
                bestucost2 = ucost;
                bestumove2 = umove;
                gotuc2 = 1;
            }
        }

    }
    sumtwo -= bestdiff + bestdiff2;
    summove -= bestmove + bestmove2;

    ucost = bestucost + bestucost2;
    umove = bestumove + bestumove2;

    for (i=0; i<ncnt; i++) {
        n = nodes[i];
        n->mark = 0;
    }
    if (pen) {
        pen->delta_bad = delta_bad;
        pen->delta = delta;
        pen->inlen = inlen;
        pen->outlen = outlen;
        pen->uppen_bad = (4.0 - delta_bad) * outlen_bad;
        pen->downpen_bad = (delta_bad - 2.0) * inlen_bad;
        pen->uppen = (4.0 - delta) * outlen;
        pen->downpen = (delta - 2.0) * inlen;
        pen->uppen2 = (4.0 - delta) * (outlen - inlen * 0.5);
        pen->downpen2 = (delta - 2.0) * (inlen * 0.5 - outlen);
        pen->downpen3 = sumtwo;
        if (summove == 0.0) pen->downpen4 = sumtwo;
        else                pen->downpen4 = sumtwo * (delta - 2.0) / summove;
        pen->uppen3 = ucost;
        if (umove == 0.0) pen->uppen4 = ucost;
        else              pen->uppen4 = ucost * (4.0 - delta) / umove;
        pen->ncnt = ncnt;
        for (i=0; i<ncnt && i < 4; i++) {
            pen->nodes[i] = nodes[i]->number;
        }
    }

/*
    *p_delta = delta_bad;
    *p_inlen = inlen_bad;
    *p_outlen = outlen_bad;
    */


    *p_delta = delta;
    *p_inlen = inlen;
    *p_outlen = outlen;

}

static int le_findmax (int ncnt, le_node **nodes, le_info *le_dat)
{
    double delta = 0.0;
    int inlen = 0;
    int outlen = 0;
    double downpen;
    double uppen;

    le_cliquevals (ncnt, nodes, &delta, &inlen, &outlen, (le_penalty *) NULL);

    if (delta <= (2.0 + TSP_BRANCH_CLIQUE_TOL) ||
        delta >= (4.0 - TSP_BRANCH_CLIQUE_TOL)) return 0;

    uppen = (4.0 - delta) * outlen;
    downpen = (delta - 2.0) * inlen;

    if (uppen > le_dat->maxup) le_dat->maxup = uppen;
    if (downpen > le_dat->maxdown) le_dat->maxdown = downpen;

    return 0;
}

static int le_getcliques (int ncnt, le_node **nodes, le_info *le_dat)
{
    double delta = 0.0;
    int inlen = 0;
    int outlen = 0;
    double downpen;
    double uppen;
    int rval = 0;
    le_penalty pen;

    le_cliquevals (ncnt, nodes, &delta, &inlen, &outlen, &pen);
    
    if (delta <= (2.0 + TSP_BRANCH_CLIQUE_TOL) ||
        delta >= (4.0 - TSP_BRANCH_CLIQUE_TOL)) return 0;
    
    uppen = (4.0 - delta) * outlen;
    downpen = (delta - 2.0) * inlen;

#ifdef DEBUG
    {
        le_node *n;
        CCtsp_lpclique cliq;
        int *arr = (int *) NULL;
        int cnt = 0;
        int i;

        CCtsp_init_lpclique (&cliq);

        for (i=0; i<ncnt; i++) {
            n = nodes[i];
            do {
                cnt++;
                n = n->next_member;
            } while (n != nodes[i]);
        }

        arr = CC_SAFE_MALLOC (cnt, int);
        if (arr == (int *) NULL) {
            return 1;
        }
        cnt = 0;

        for (i=0; i<ncnt; i++) {
            n = nodes[i];
            do {
                arr[cnt++] = n->number;
                n = n->next_member;
            } while (n != nodes[i]);
        }

        rval = CCtsp_array_to_lpclique (arr, cnt, &cliq);
        if (rval) {
            CC_IFFREE (arr, int);
            CCtsp_free_lpclique (&cliq);
            return rval;
        }
        
        printf ("clique ");
        CCtsp_print_lpclique (&cliq);
        printf ("delta %.6f outlen %d inlen %d uppen %.6f downpen %.6f\n",
                delta, outlen, inlen, uppen, downpen);
        fflush (stdout);
        CCtsp_free_lpclique (&cliq);
        CC_IFFREE (arr, int);
    }
#endif /* DEBUG */

    uppen /= le_dat->maxup;
    downpen /= le_dat->maxdown;

    rval = newclique2 (le_dat->nwant, le_dat->clist, le_dat->cval, ncnt, nodes,
                       downpen, uppen, le_dat->workarr, &pen,
                       le_dat->penalties);
    if (rval) {
        fprintf (stderr, "newclique2 failed\n");
        return rval;
    }

    return rval;
}

static int le_foreach_clique (le_graph *g, int (*func)(int ncnt,
        le_node **nodes, le_info *le_dat), le_info *le_dat)
{
    int i;
    int nedges = g->nedges;
    int nnodes = g->nnodes;
    le_edge *edges = g->edges;
    le_node *n;
    le_adj *a;
#ifdef LONGEDGE_TRIANGLES
    le_adj *b;
    le_node *m;
#endif
    int rval;
    le_node *cnodes[3];

    for (i=0; i<nedges; i++) {
        cnodes[0] = edges[i].ends[0];
        cnodes[1] = edges[i].ends[1];
        rval = (*func)(2, cnodes, le_dat);
        if (rval) {
            fprintf (stderr, "le_foreach_clique callback failed\n");
            return rval;
        }
    }

#ifdef LONGEDGE_TRIANGLES
    for (i=0; i<nnodes; i++) {
        n = &g->nodes[i];
        cnodes[0] = n;
        for (a = n->adj; a; a = a->next) {
            if (a->this->live) {
                cnodes[1] = LE_OTHEREND (a->this, n);
                for (b = a->next; b; b = b->next) {
                    if (b->this->live) {
                        cnodes[2] = LE_OTHEREND (b->this, n);
                    }
                    rval = (*func)(3, cnodes, le_dat);
                    if (rval) {
                        fprintf (stderr, "le_foreach_clique callback failed\n");
                        return rval;
                    }
                }
            }
        }
    }
#endif
                
    for (i=0; i<nedges && nnodes >= 4; i++) {
        if (edges[i].live) {
            if (edges[i].ends[0] == edges[i].ends[1]) {
                fprintf (stderr, "Whoops, trying to contract self-loop\n");
                continue;
            }
#ifdef DEBUG
            printf ("contracting %d: %d - %d\n",
                    i, edges[i].ends[0]->number, edges[i].ends[1]->number);
#endif /* DEBUG */
            n = le_contract (edges[i].ends[0], edges[i].ends[1]);
            for (a = n->adj; a; a = a->next) {
                if (a->this->live) {
                    cnodes[0] = a->this->ends[0];
                    cnodes[1] = a->this->ends[1];
                    rval = (*func)(2, cnodes, le_dat);
                    if (rval) {
                        fprintf (stderr, "le_foreach_clique callback failed\n");
                        return rval;
                    }
#ifdef LONGEDGE_TRIANGLES
                    if (nnodes >= 5) {
                        for (b = a->next; b; b = b->next) {
                            if (b->this->live) {
                                cnodes[2] = LE_OTHEREND (b->this, n);
                                rval = (*func) (3, cnodes, le_dat);
                                if (rval) {
                                    fprintf (stderr, "le_foreach_clique callback failed\n");
                                    return rval;
                                }
                            }
                        }
                        m = LE_OTHEREND (a->this, n);
                        for (b = m->adj; b; b = b->next) {
                            if (b->this->live) {
                                cnodes[2] = LE_OTHEREND (b->this, m);
                                if (cnodes[2] != n) {
                                    rval = (*func) (3, cnodes, le_dat);
                                    if (rval) {
                                        fprintf (stderr, "le_foreach_clique callback failed\n");
                                        return rval;
                                    }
                                }
                            }
                        }
                    }
#endif

                }
            }
            nnodes--;
        }
    }
        
    return 0;
}

static int CC_UNUSED le_foreach_clique2 (le_graph *g, CCdatagroup *dat,
        int (*func)(int ncnt, le_node **nodes, le_info *le_dat),
        le_info *le_dat)
{
    int i;
    int nedges = g->nedges;
    int nnodes = g->nnodes;
    le_edge *edges = g->edges;
    le_node *n;
    le_adj *a;
    int *len = (int *) NULL;
    int *perm = (int *) NULL;
    le_node *cnodes[2];
    int rval = 0;

    len = CC_SAFE_MALLOC (nnodes, int);
    perm = CC_SAFE_MALLOC (nnodes, int);
    if (len == (int *) NULL ||
        perm == (int *) NULL) {
        fprintf (stderr, "Out of memory in le_foreach_clique2\n");
        rval = 1; goto CLEANUP;
    }
    len[0] = CCutil_dat_edgelen (0, nnodes-1, dat);
    for (i=1; i<nnodes; i++) {
        len[i] = CCutil_dat_edgelen (i-1, i, dat);
        perm[i] = i;
    }
    
    for (i=0; i<nedges; i++) {
        cnodes[0] = edges[i].ends[0];
        cnodes[1] = edges[i].ends[1];
        rval = (*func)(2, cnodes, le_dat);
        if (rval) {
            fprintf (stderr, "le_foreach_clique callback failed\n");
            goto CLEANUP;
        }
    }

    for (i=0; i<g->nnodes && nnodes >= 4; i++) {
        if (edges[i].live && edges[i].ends[0] != edges[i].ends[1]) {
#ifdef DEBUG
            printf ("contracting %d: %d - %d\n",
                    i, edges[i].ends[0]->number, edges[i].ends[1]->number);
#endif /* DEBUG */
            n = le_contract (edges[i].ends[0], edges[i].ends[1]);
            for (a = n->adj; a; a = a->next) {
                cnodes[0] = a->this->ends[0];
                cnodes[1] = a->this->ends[1];
                rval = (*func)(2, cnodes, le_dat);
                if (rval) {
                    fprintf (stderr, "le_foreach_clique callback failed\n");
                    goto CLEANUP;
                }
            }
            nnodes--;
        }
    }

 CLEANUP:
    CC_IFFREE (len, int);
    CC_IFFREE (perm, int);
    
    return rval;
}


static int find_longedge_cliques2 (CCtsp_lp *lp, int nwant, int *ngot,
        CCtsp_lpclique **list, int silent)
{
    int ncount = lp->graph.ncount;
    le_graph g;
    int rval = 0;
    int *len = (int *) NULL;
    int xcount;
    int *xlist = (int *) NULL;
    double *x = (double *) NULL;
    CCtsp_lpclique *clist = (CCtsp_lpclique *) NULL;
    double *cval = (double *) NULL;
    double st;
    int *workarr = (int *) NULL;
    int i, k;
    CCutil_timer timer;
    le_penalty *penalties = (le_penalty *) NULL;
    le_info le_dat;

    *ngot = 0;
    init_le_graph (&g);
    CCutil_init_timer (&timer, "Long-edge search");
    CCutil_start_timer (&timer);

    if (!silent) {
        printf ("Finding long-edge cliques\n"); fflush (stdout);
    }
    
    if (nwant <= 0) {
        fprintf (stderr, "find_longedge_cliques called with no nwant\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL, &xcount,
                                &xlist, &x, (double **) NULL,
                                (double **) NULL, (double **) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_get_lp_result failed\n");
        goto CLEANUP;
    }

    if (xcount != lp->graph.ecount) {
        fprintf (stderr, "Whoops, xcount %d ecount %d\n", xcount,
                 lp->graph.ecount);
        rval = 1; goto CLEANUP;
    }

    len = CC_SAFE_MALLOC (xcount, int);
    if (len == (int *) NULL) {
        fprintf (stderr, "Out of memory in find_longedge_cliques\n");
        rval = 1; goto CLEANUP;
    }
    for (i=0; i<xcount; i++) {
        len[i] = lp->graph.edges[i].len;
    }

#if 0
    rval = build_le_graph (&g, ncount, xcount, xlist, len, x);
    if (rval) {
        fprintf (stderr, "build_le_graph failed\n");
        goto CLEANUP;
    }

    le_dat.maxup = 0.0;
    le_dat.maxdown = 0.0;
    
    rval = le_foreach_clique (&g, le_findmax, &le_dat);
    if (rval) {
        fprintf (stderr, "le_foreach_clique failed\n");
        goto CLEANUP;
    }

/*
    rval = le_foreach_clique2 (&g, lp->dat, le_findmax, &le_dat);
    if (rval) {
        fprintf (stderr, "le_foreach_clique2 failed\n");
        goto CLEANUP;
    }
*/

    free_le_graph (&g);
    
    if (le_dat.maxdown == 0.0) le_dat.maxdown = 1.0;
    if (le_dat.maxup == 0.0) le_dat.maxup = 1.0;

#endif /* 0 */

    /*****/
    le_dat.maxdown = 10.0;
    le_dat.maxup = 3.0;
    
    clist = CC_SAFE_MALLOC (nwant, CCtsp_lpclique);
    cval = CC_SAFE_MALLOC (nwant, double);
    penalties = CC_SAFE_MALLOC (nwant, le_penalty);
    workarr = CC_SAFE_MALLOC (ncount, int);
    if (clist == (CCtsp_lpclique *) NULL ||
        cval == (double *) NULL ||
        penalties == (le_penalty *) NULL ||
        workarr == (int *) NULL) {
        fprintf (stderr, "Out of memory in find_longedge_cliques\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<nwant; i++) {
        cval[i] = -1.0;
        CCtsp_init_lpclique (&clist[i]);
    }

    le_dat.nwant = nwant;
    le_dat.clist = clist;
    le_dat.penalties = penalties;
    le_dat.cval = cval;
    le_dat.workarr = workarr;
    
    rval = build_le_graph (&g, ncount, xcount, xlist, len, x);
    if (rval) {
        fprintf (stderr, "build_le_graph failed\n");
        goto CLEANUP;
    }

    rval = le_foreach_clique (&g, le_getcliques, &le_dat);
    if (rval) {
        fprintf (stderr, "le_foreach_clique failed\n");
        goto CLEANUP;
    }

    free_le_graph (&g);

    k = 0;
    for (i = 0, k=0; i<nwant; i++) {
        if (cval[i] >= 0.0) {
            k++;
        }
    }

    if (k == 0) {
        rval = 0;
        goto CLEANUP;
    }

    *list = CC_SAFE_MALLOC (k, CCtsp_lpclique);
    if (!(*list)) {
        fprintf (stderr, "Out of memory in find_longedge_cliques\n");
        rval = 1; goto CLEANUP;
    }
    *ngot = k;

    for (i=0, k=0; i<nwant; i++) {
        if (cval[i] >= 0.0) {
            if (!silent) {
                int j;
                printf ("%d: (%.6f) ", i, cval[i]);
                CCtsp_print_lpclique (&clist[i]);
                printf ("%d: delta %.6f inlen %d outlen %d uppen_bad %.6f downpen_bad %.6f uppen %.6f downpen %.6f uppen2 %.6f downpen2 %.6f uppen3 %.6f downpen3 %.6f uppen4 %.6f downpen4 %.6f\n",
                    i, penalties[i].delta, penalties[i].inlen, penalties[i].outlen,
                    penalties[i].uppen_bad, penalties[i].downpen_bad, 
                    penalties[i].uppen, penalties[i].downpen, 
                    penalties[i].uppen2, penalties[i].downpen2, 
                    penalties[i].uppen3, penalties[i].downpen3, 
                    penalties[i].uppen4, penalties[i].downpen4);
                printf ("%d: nodes", i);
                for (j=0; j<penalties[i].ncnt && j < 4; j++) {
                    printf (" %d", penalties[i].nodes[j]);
                }
                printf ("\n");
            }
            (*list)[k] = clist[i];
            k++;
        }
    }

    st = CCutil_stop_timer (&timer, 0);
    if (!silent) {
        printf ("%d long-edge cliques found in %.2f seconds\n", *ngot, st);
        fflush (stdout);
    }

    CC_IFFREE (clist, CCtsp_lpclique);
    rval = 0;

CLEANUP:
    if (clist) {
        for (i=0; i<nwant; i++) {
            CCtsp_free_lpclique (&clist[i]);
        }
    }
    free_le_graph (&g);
    CC_IFFREE (clist, CCtsp_lpclique);
    CC_IFFREE (penalties, le_penalty);
    CC_IFFREE (cval, double);
    CC_IFFREE (x, double);
    CC_IFFREE (xlist, int);
    CC_IFFREE (len, int);
    CC_IFFREE (workarr, int);

    return rval;
}

static int merge_edge_clique (CCtsp_lp *lp, int nwant, int *ngot,
        CCtsp_branchobj **bobj, int ecount, int *elist, double *eval,
        int ccount, CCtsp_lpclique *clist, double *cval)
{
    int rval = 0;
    int i, k;
    sbitem *slist = (sbitem *) NULL;
    CCtsp_branchobj *b;

    *ngot = 0;
    *bobj = (CCtsp_branchobj *) NULL;

    if (ecount + ccount == 0) {
        fprintf (stderr, "no elements in merge_edge_clique\n");
        rval = 1; goto CLEANUP;
    }

    slist = CC_SAFE_MALLOC (nwant + 1, sbitem);
    if (!slist) {
        fprintf (stderr, "out of memory in merge_edge_clique\n");
        rval = 1; goto CLEANUP;
    }
    init_sblist (slist, nwant);

    for (i = 0; i < ecount; i++) {
        insert_sblist (slist, eval[i], i);
    }
    for (i = 0; i < ccount; i++) {
        insert_sblist (slist, cval[i], i + ecount);
    }

    for (i = 0, k = 0; i < nwant; i++) {
        if (slist[i].name != -1) {
            k++;
        }
    }
    if (k == 0) {
        fprintf (stderr, "nothing appeares in merge_edge_clique\n");
        rval = 1; goto CLEANUP;
    }

    *bobj = CC_SAFE_MALLOC (k, CCtsp_branchobj);
    if (!(*bobj)) {
        fprintf (stderr, "out of memory in merge_edge_clique\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0, k = 0; i < nwant; i++) {
        if (slist[i].name != -1) {
            b = &((*bobj)[k]);
            CCtsp_init_branchobj (b);
            if (slist[i].name < ecount) {
                b->ends[0] = lp->graph.edges[elist[slist[i].name]].ends[0];
                b->ends[1] = lp->graph.edges[elist[slist[i].name]].ends[1];
            } else {
                b->clique = CC_SAFE_MALLOC (1, CCtsp_lpclique);
                if (!b->clique) {
                    fprintf (stderr, "out of memory in merge_edge_clique\n");
                    rval = 1; goto CLEANUP;
                } else {
                    rval = CCtsp_copy_lpclique (&clist[slist[i].name - ecount],
                                                b->clique);
                }
                if (!b->clique || rval) {
                    fprintf (stderr, "CCtsp_copy_clique failed\n");
                    for (i = 0; i < k; i++) {
                        if ((*bobj)[i].clique) {
                            CCtsp_free_lpclique ((*bobj)[i].clique);
                            CC_IFFREE ((*bobj)[i].clique, CCtsp_lpclique);
                        }
                    }
                    CC_IFFREE (b->clique, CCtsp_lpclique);
                    CC_FREE (*bobj, CCtsp_branchobj);
                    goto CLEANUP;
                }
            }
            k++;
        }
    }
    *ngot = k;

CLEANUP:

    CC_IFFREE (slist, sbitem);
    return rval;
}

int CCtsp_check_integral (CCtsp_lp *lp, double *val, int **cyc, int *yesno,
        int silent)
{
    int rval = 0;
    double *x = (double *) NULL;
    double eval = 0.0;
    int *xlist = (int *) NULL;
    int xcount;
    int *comps = (int *) NULL;
    int *compscount = (int *) NULL;
    int ncomp;
    int *elist = (int *) NULL;
    int ncount = lp->graph.ncount;
    int i, j, ecount;

    *yesno = 0;
    *val = 0.0;
    if (cyc) *cyc = (int *) NULL;

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL, &xcount,
                          &xlist, &x, (double **) NULL, (double **) NULL,
                          (double **) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_get_lp_result failed\n"); goto CLEANUP;
    }

    for (i = 0; i < xcount; i++) {
        if (x[i] > 0.5) {
            if (1.0 - x[i] > CCtsp_INTTOL) goto CLEANUP;
        } else {
            if (x[i] > CCtsp_INTTOL) goto CLEANUP;
        }
    }

    elist = CC_SAFE_MALLOC (2*ncount, int);
    if (!elist) {
        fprintf (stderr, "out of memory in CCtsp_check_integral\n");
    }
    ecount = 0;

    for (i = 0; i < xcount; i++) {
        if (x[i] > CCtsp_INTTOL) {
            j = CCtsp_find_edge (&lp->graph, xlist[2*i], xlist[2*i+1]);
            if (j < 0) {
                fprintf (stderr, "x edge not in graph\n");
                rval = 1; goto CLEANUP;
            }
            eval += ((double) lp->graph.edges[j].len);
            elist[2*ecount]   = lp->graph.edges[j].ends[0];
            elist[2*ecount+1] = lp->graph.edges[j].ends[1];
            ecount++;
        }
    }
    rval = CCcut_connect_components (ncount, ecount, elist, (double *) NULL,
                                     &ncomp, &compscount, &comps);
    if (rval) {
        fprintf (stderr, "CCcut_connect_components failed\n"); goto CLEANUP;
    }
    if (ncomp > 1) {
        if (!silent) {
            printf ("integral solution not connected\n"); fflush (stdout);
        }
        goto CLEANUP;
    }
    if (!silent) {
        printf ("Integral Solution of Value %.2f\n", eval); fflush (stdout);
    }

    if (cyc) {
        int istour;
        *cyc = CC_SAFE_MALLOC (ncount, int);
        if (!(*cyc)) {
            fprintf (stderr, "out of memory in CCtsp_check_integral\n");
            rval = 1; goto CLEANUP;
        }
        rval = CCutil_edge_to_cycle (ncount, elist, &istour, *cyc);
        if (rval) {
            fprintf (stderr, "CCutil_edge_to_cycle failed\n");
            CC_FREE (*cyc, int);
            goto CLEANUP;
        }
        if (istour == 0) {
            fprintf (stderr, "ERROR: integral, connected solution not tour\n");
            rval = 1; goto CLEANUP;
        }
    }
    *yesno = 1;
    *val = eval;

CLEANUP:

    CC_IFFREE (x, double);
    CC_IFFREE (xlist, int);
    CC_IFFREE (comps, int);
    CC_IFFREE (compscount, int);
    CC_IFFREE (elist, int);

    return rval;
}

int CCtsp_find_branch_edge (CCtsp_lp *lp, int *n0, int *n1, double *val,
        int **cyc, int branchtype, int silent)
{
    int rval = 0;
    double *x = (double *) NULL;
    int *xlist = (int *) NULL;
    int xcount;
    double maxdiff;
    int i, test, besti;
    int j;

    *n0 = -2;
    *n1 = -2;
    *val = 0.0;
    if (cyc) *cyc = (int *) NULL;

    rval = CCtsp_check_integral (lp, val, cyc, &test, silent);
    if (rval) {
        fprintf (stderr, "CCtsp_check_integral failed\n");
        goto CLEANUP;
    }
    if (test) {
        if (!silent) {
            printf ("Integral solution detected in CCtsp_find_branch_edge\n");
            fflush (stdout);
        }
        *n0 = -1;
        *n1 = -1;
        goto CLEANUP;
    }

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL, &xcount,
                          &xlist, &x, (double **) NULL, (double **) NULL,
                          (double **) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_get_lp_result failed\n"); goto CLEANUP;
    }

    besti = -1;
    maxdiff = -1.0;

    for (i = 0; i < xcount; i++) {
        j = CCtsp_find_edge (&lp->graph, xlist[2*i], xlist[2*i+1]);
        if (j < 0) {
            fprintf (stderr, "edge should be in LP\n");
            rval = 1; goto CLEANUP;
        }
        if (!lp->graph.edges[j].fixed && !lp->graph.edges[j].branch) {
            if (x[i] > 0.5) {
                if (1.0 - x[i] > maxdiff) {
                    maxdiff = 1.0 - x[i];
                    besti = i;
                }
            } else {
                if (x[i] > maxdiff) {
                    maxdiff = x[i];
                    besti = i;
                }
            }
        }
    }

    if (besti == -1) {
        fprintf (stderr, "All edges are either branched or fixed\n");
        rval = 1; goto CLEANUP;
    }
    
    switch (branchtype) {
    case CCtsp_BRANCH_MIDDLE:
        *n0 = xlist[2*besti];
        *n1 = xlist[2*besti+1];
        break;
    case CCtsp_BRANCH_STRONG:
        rval = find_strong_branch (lp, n0, n1, silent);
        if (rval) {
            fprintf (stderr, "find_strong_branch failed\n");
            goto CLEANUP;
        }
        if (*n0 == -1) {
            *n0 = xlist[2*besti];
            *n1 = xlist[2*besti+1];
        }
        break;
    default:
        fprintf (stderr, "unknown branchtype\n");
        rval = 1; goto CLEANUP;
    }

CLEANUP:

    CC_IFFREE (x, double);
    CC_IFFREE (xlist, int);

    return rval;
}

static int find_strong_branch (CCtsp_lp *lp, int *n0, int *n1, int silent)
{
    int rval = 0;
    int *elist = (int *) NULL;
    int ngot;

    *n0 = -1;
    *n1 = -1;

    rval = find_strongbranch_edges (lp, 1, &ngot, &elist, (double **) NULL,
                                    silent);
    if (rval) {
        fprintf (stderr, "find_strongbranch_edges failed\n");
        goto CLEANUP;
    }

    if (ngot == 0) {
        if (!silent) {
            printf ("WARNING: nothing from find_strongbranch_edges\n");
            fflush (stdout);
        }
        goto CLEANUP;
    }

    *n0 = lp->graph.edges[elist[0]].ends[0];
    *n1 = lp->graph.edges[elist[0]].ends[1];

    if (!silent) {
        printf ("STRONG branch edge: %d %d\n", *n0, *n1);
        fflush (stdout);
    }

CLEANUP:

    CC_IFFREE (elist, int);
    return rval;
}


static int find_strongbranch_edges (CCtsp_lp *lp, int nwant, int *ngot,
        int **elist, double **eval, int silent)
{
    int rval = 0;
    int *candlist = (int *) NULL;
    int ncand = 0;
    int lpcand = 0;
    int i, k, lpwant;
    double *downpen = (double *) NULL;
    double *uppen   = (double *) NULL;
    double sval;
    double meanval = 0.0;
    sbitem *slist = (sbitem *) NULL;

    *ngot = 0;
    *elist = (int *) NULL;
    if (eval) {
        *eval = (double *) NULL;
    }

    rval = find_candidate_edges (lp, nwant * TSP_BRANCH_STRONG_CHOICES_MULT,
                             &ncand, &candlist, silent);
    if (rval) {
        fprintf (stderr, "find_candidate_edges failed\n");
        goto CLEANUP;
    }

    if (ncand == 0) {
        if (!silent) {
            printf ("WARNING: find_candidate edges did not find anything\n");
            fflush (stdout);
        }
        goto CLEANUP;
    }

    if (!silent) {
        printf ("Run strongbranch with %d candidate edges\n", ncand);
        fflush (stdout);
    }

    downpen = CC_SAFE_MALLOC (ncand, double);
    uppen   = CC_SAFE_MALLOC (ncand, double);
    if (!downpen || !uppen) {
        fprintf (stderr, "out of memory in find_strongbranch_edges\n");
        rval = 1; goto CLEANUP;
    }

    CCutil_start_timer (&lp->stats.strongbranch_opt);
    rval = CClp_strongbranch (lp->lp, candlist, ncand, downpen, uppen,
                TSP_BRANCH_STRONG_ITERATIONS, lp->upperbound);
    if (rval) {
        fprintf (stderr, "CClp_strongbranch failed\n"); goto CLEANUP;
    }
    if (!silent) {
        CCutil_stop_timer (&lp->stats.strongbranch_opt, 1);
    } else {
        CCutil_stop_timer (&lp->stats.strongbranch_opt, 0);
    }

    lpwant = TSP_BRANCH_STRONG_LP_CHOICES_MULT * nwant;
    slist = CC_SAFE_MALLOC (lpwant + 1, sbitem);
    if (!slist) {
        fprintf (stderr, "out of memory in find_strongbranch_edges\n");
        rval = 1; goto CLEANUP;
    }
    init_sblist (slist, lpwant);

    for (i = 0; i < ncand; i++) {
        if (!silent) {
            printf ("SB Edge (%d, %d):  %.4f  %.4f\n",
                    lp->graph.edges[candlist[i]].ends[0],
                    lp->graph.edges[candlist[i]].ends[1],
                    downpen[i], uppen[i]);
            fflush (stdout);
        }
        sval = TSP_BRANCH_STRONG_VAL (downpen[i], uppen[i]);
        insert_sblist (slist, sval, candlist[i]);
        meanval += sval;
    }

    if (ncand && !silent) {
        printf ("Average Edge Value: %f\n", meanval / ((double) ncand));
        fflush (stdout);
    }
    
    for (i = lpwant - 1, lpcand = 0; i >= 0; i--) {
        if (slist[i].name != -1) {
            if (!silent) {
                printf ("First Stage Top Edge: (%d, %d)  (%f)\n",
                        lp->graph.edges[slist[i].name].ends[0],
                        lp->graph.edges[slist[i].name].ends[1],
                        slist[i].val);
                fflush (stdout);
            }
            candlist[lpcand++] = slist[i].name;
        }
    }
 
    if (lpcand == 0 && !silent) {
        printf ("WARNING: no edges appeared in strongbranch\n");
        goto CLEANUP;
    }

    CCutil_start_timer (&lp->stats.strongbranch_opt);
    rval = CClp_strongbranch (lp->lp, candlist, lpcand, downpen, uppen,
                TSP_BRANCH_STRONG_EXTRA_ITERATIONS, lp->upperbound);
    if (rval) {
        fprintf (stderr, "CClp_strongbranch failed\n"); goto CLEANUP;
    }
    if (!silent) {
        CCutil_stop_timer (&lp->stats.strongbranch_opt, 1);
    } else {
        CCutil_stop_timer (&lp->stats.strongbranch_opt, 0);
    }

    meanval = 0.0;
    init_sblist (slist, nwant);
    for (i = 0; i < lpcand; i++) {
        if (!silent) {
            printf ("SB2 Edge (%d, %d):  %.4f  %.4f\n",
                        lp->graph.edges[candlist[i]].ends[0],
                        lp->graph.edges[candlist[i]].ends[1],
                        downpen[i], uppen[i]);
            fflush (stdout);
        }
        sval = TSP_BRANCH_STRONG_VAL (downpen[i], uppen[i]);
        insert_sblist (slist, sval, candlist[i]);
        meanval += sval;
    }

    if (lpcand && !silent) {
        printf ("Average Edge Value: %f\n", meanval / ((double) lpcand));
        fflush (stdout);
    }
    
    for (i = nwant - 1, k = 0; i >= 0; i--) {
        if (slist[i].name != -1) {
            k++;
            if (!silent) {
                printf ("Top Edge: (%d, %d)  (%f)\n",
                    lp->graph.edges[slist[i].name].ends[0],
                    lp->graph.edges[slist[i].name].ends[1],
                    slist[i].val);
                fflush (stdout);
            }
        }
    }

    if (k == 0) {
        if (!silent) {
            printf ("WARNING: no edges appeared in strongbranch\n");
            fflush (stdout);
        }
        goto CLEANUP;
    }

    *elist = CC_SAFE_MALLOC (k, int);
    if (!(*elist)) {
        fprintf (stderr, "out of memory in find_strongbranch_edges\n");
        rval = 1; goto CLEANUP;
    }
    if (eval) {
        *eval = CC_SAFE_MALLOC (k, double);
        if (!(*eval)) {
            fprintf (stderr, "out of memory in find_strongbranch_edges\n");
            CC_IFFREE (*elist, int);
            rval = 1; goto CLEANUP;
        }
    }
    *ngot = k;
    for (i = 0, k = 0; i < nwant; i++) {
        if (slist[i].name != -1) {
            (*elist)[k] = slist[i].name;
            if (eval) {
                (*eval)[k] = slist[i].val;
            }
            k++;
        }
    }

CLEANUP:

    CC_IFFREE (slist, sbitem);
    CC_IFFREE (candlist, int);
    CC_IFFREE (downpen, double);
    CC_IFFREE (uppen, double);

    return rval;
}

static int find_candidate_edges (CCtsp_lp *lp, int nwant, int *ngot,
        int **list, int silent)
{
    int *goodlist     = (int *) NULL;
    double *downpen = (double *) NULL;
    double *uppen   = (double *) NULL;
    double sval;
    int nrows, ngood;
    int i, k;
    sbitem *slist = (sbitem *) NULL;
    int rval = 0;

    *ngot = 0;
    *list = (int *) NULL;

    if (nwant == TSP_BRANCH_STRONG_ALL_CHOICES) {
        rval = find_all_candidate_edges (lp, ngot, list);
        if (rval) {
            fprintf (stderr, "find_all_candidate_edges failed\n");
            goto CLEANUP;
        }
        goto CLEANUP;
    }

    slist = CC_SAFE_MALLOC (nwant + 1, sbitem);
    if (!slist) {
        fprintf (stderr, "out of memory in find_strongbranch\n");
        rval = 1; goto CLEANUP;
    }
    init_sblist (slist, nwant);

    nrows     = lp->graph.ncount + lp->cuts.cutcount;
    goodlist  = CC_SAFE_MALLOC (nrows, int);
    downpen = CC_SAFE_MALLOC (nrows, double);
    uppen   = CC_SAFE_MALLOC (nrows, double);
    if (!goodlist || !downpen || !uppen) {
        fprintf (stderr, "out of memory in find_strongbranch\n");
        rval = 1; goto CLEANUP;
    }
    rval = CClp_getgoodlist (lp->lp, goodlist, &ngood, downpen, uppen);
    if (rval) {
        fprintf (stderr, "CClp_getgoodlist failed\n"); goto CLEANUP;
    }
    if (!silent) {
        printf ("Found %d good edges\n", ngood); fflush (stdout);
    }

    for (i = 0; i < ngood; i++) {
        if (downpen[i] > lp->upperbound) downpen[i] = lp->upperbound;
        if (uppen[i] > lp->upperbound) uppen[i] = lp->upperbound;
        sval = TSP_BRANCH_STRONG_FIRST_VAL (downpen[i], uppen[i]);
        insert_sblist (slist, sval, goodlist[i]);
    }

    for (i = 0, k = 0; i < nwant; i++) {
        if (slist[i].name != -1) {
            k++;
        }
    }
    if (k == 0) {
        if (!silent) {
            printf ("WARNING: CClp_getgoodlist returned no edges\n");
            fflush (stdout);
        }
        goto CLEANUP;
    }

    *list = CC_SAFE_MALLOC (k, int);
    if (!(*list)) {
        fprintf (stderr, "out of memory in find_candidate list\n");
        rval = 1; goto CLEANUP;
    }
    *ngot = k;
    for (i = 0, k = 0; i < nwant; i++) {
        if (slist[i].name != -1) {
            (*list)[k++] = slist[i].name;
        }
    }

CLEANUP:

    CC_IFFREE (goodlist, int);
    CC_IFFREE (downpen, double);
    CC_IFFREE (uppen, double);
    CC_IFFREE (slist, sbitem);

    return rval;
}

static int find_all_candidate_edges (CCtsp_lp *lp, int *ngot, int **list)
{
    int rval = 0;
    double *x  = (double *) NULL;
    int *xlist = (int *) NULL;
    int i, j, xcount, count;

    *ngot = 0;
    *list = (int *) NULL;

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL, &xcount,
                          &xlist, &x, (double **) NULL, (double **) NULL,
                          (double **) NULL);
    if (rval) {
        fprintf (stderr, "get_lp_result failed\n"); goto CLEANUP;
    }

    count = 0;
    for (i = 0; i < xcount; i++) {
        if (x[i] >= CCtsp_INTTOL && x[i] <= 1.0 - CCtsp_INTTOL) {
            count++;
        }
    }

    if (!count) {
        fprintf (stderr, "WARNING: The solution is integral\n");
        goto CLEANUP;
    }

    *list = CC_SAFE_MALLOC (count, int);
    if (!(*list)) {
        fprintf (stderr, "out of memory in find_all_candidate_edges\n");
        rval = 1; goto CLEANUP;
    }

    count = 0;
    for (i = 0; i < xcount; i++) {
        if (x[i] >= CCtsp_INTTOL && x[i] <= 1.0 - CCtsp_INTTOL) {
            j = CCtsp_find_edge (&lp->graph, xlist[2 * i], xlist[2 * i + 1]);
            if (j < 0) {
                fprintf (stderr, "edge not in lp in find_all_candiate_edges\n");
                CC_IFFREE (*list, int);
                rval = 1; goto CLEANUP;
            }
            (*list)[count++] = j;
        }
    }
    *ngot = count;

CLEANUP:

    CC_IFFREE (x, double);
    CC_IFFREE (xlist, int);

    return rval;
}

static void init_sblist (sbitem *list, int count)
{
    int i;

    for (i = 0; i < count; i++) {
        list[i].val   = -1.0;
        list[i].name  = -1;
        list[i].name2 = -1;
        list[i].name3 = -1;
        list[i].name4 = -1;
    }
    list[count].val = CCtsp_LP_MAXDOUBLE;
}

static void insert_sblist (sbitem *list, double val, int name)
{
    insert_sblist4 (list, val, name, -1, -1, -1);
}

static void insert_sblist4 (sbitem *list, double val, int name, int name2,
        int name3, int name4)
{
    int k;

    if (list[0].val < val) {
        for (k = 0; list[k+1].val < val; k++) {
            list[k] = list[k+1];
        }
        list[k].val   = val;
        list[k].name  = name;
        list[k].name2 = name2;
        list[k].name3 = name3;
        list[k].name4 = name4;
    }
}

int CCtsp_find_branch_cliques (CCtsp_lp *lp, int nwant, int longedge_branching,
        int *ngot, CCtsp_lpclique **bcliques, double **bval, int silent)
{
    int      ccount   = 0;
    CCtsp_lpclique *cliques = (CCtsp_lpclique *) NULL;
    sbitem   *slist   = (sbitem *) NULL;
    int lpcount = 0;
    int i, k, lpwant;
    int rval = 0;
    int *lpcand = (int *) NULL;
    double up, down, sval;
    double meanval = 0.0;
    CCutil_timer timer;
    CClp_warmstart *warmstart = (CClp_warmstart *) NULL;

    CCutil_init_timer (&timer, "Strong cut branch");
    
    *ngot     = 0;
    *bcliques = (CCtsp_lpclique *) NULL;
    if (bval) {
        *bval = (double *) NULL;
    }

    if (longedge_branching) {
        rval = find_longedge_cliques2 (lp, nwant * TSP_STRONG_CUT_CHOICES_MULT,
                                       &ccount, &cliques, silent);
        if (rval) {
            fprintf (stderr, "find_longedge_cliques failed\n");
            goto CLEANUP;
        }
    } else {
        rval = find_candidate_cliques (lp, nwant * TSP_STRONG_CUT_CHOICES_MULT,
                                       &ccount, &cliques, 1, silent);
        if (rval) {
            fprintf (stderr, "find_candidate_cliques failed\n"); goto CLEANUP;
        }
    }
    if (!silent) {
        printf ("Found %d candidate cliques\n", ccount); fflush (stdout);
    }

    lpwant = TSP_STRONG_CUT_LP_CHOICES_MULT * nwant;
    slist = CC_SAFE_MALLOC (lpwant + 1, sbitem);
    if (!slist) {
        fprintf (stderr, "out of memory in find_branch_cliques\n");
        rval = 1; goto CLEANUP;
    }
    init_sblist (slist, lpwant);

    rval = CClp_get_warmstart (lp->lp, &warmstart);
    if (rval) {
        fprintf (stderr, "CClp_get_warmstart failed\n"); goto CLEANUP;
    }
    for (i = 0; i < ccount; i++) {
        double st;
        CCutil_start_timer (&timer);
        rval = test_cut_branch (lp, &cliques[i], &down, &up,
                              TSP_BRANCH_STRONG_ITERATIONS, silent, warmstart);
        if (rval) {
            fprintf (stderr, "test_cut_branch failed\n");
            goto CLEANUP;
        }
        st = CCutil_stop_timer (&timer, 0);
        if (!silent) {
            printf ("SB CLIQUE %d:  %f  %f  (%.2f seconds)\n",
                       i, down, up, st);
            fflush (stdout);
        }
        sval = TSP_BRANCH_STRONG_VAL (down, up);
        insert_sblist (slist, sval, i);
        meanval += sval;
    }

    rval = CClp_opt (lp->lp, CClp_METHOD_DUAL);
    if (rval) {
        fprintf (stderr, "CClp_opt failed\n");
    }

    if (ccount && !silent) {
        printf ("Average Clique Value: %f\n", meanval / ((double) ccount));
        fflush (stdout);
    }

    lpcand = CC_SAFE_MALLOC (lpwant, int);
    if (lpcand == (int *) NULL) {
        fprintf (stderr, "out of memory in find_branch_cliques\n");
        rval = 1; goto CLEANUP;
    }
    for (i = lpwant - 1, lpcount = 0; i >= 0; i--) {
        if (slist[i].name != -1) {
            if (!silent) {
                printf ("First Stage Top Clique: %d  (%f)\n", slist[i].name,
                                                              slist[i].val);
                fflush (stdout);
            }
            lpcand[lpcount++] = slist[i].name;
        }
    }
    if (lpcount == 0) {
        if (!silent) {
            printf ("WARNING: no branching cliques were found\n");
            fflush (stdout);
        }
        goto CLEANUP;
    }

    meanval = 0.0;
    init_sblist (slist, nwant);
    for (i = 0; i < lpcount; i++) {
        double st;
        CCutil_start_timer (&timer);
        rval = test_cut_branch (lp, &cliques[lpcand[i]], &down, &up,
                      TSP_BRANCH_STRONG_EXTRA_ITERATIONS, silent, warmstart);
        if (rval) {
            fprintf (stderr, "test_cut_branch failed\n");
            goto CLEANUP;
        }
        st = CCutil_stop_timer (&timer, 0);
        if (!silent) {
            printf ("SB2 CLIQUE %d:  %f  %f  (%.2f seconds)\n",
                       i, down, up, st);
            fflush (stdout);
        }
        sval = TSP_BRANCH_STRONG_VAL (down, up);
        insert_sblist (slist, sval, lpcand[i]);
        meanval += sval;
    }

    rval = CClp_opt (lp->lp, CClp_METHOD_DUAL);
    if (rval) {
        fprintf (stderr, "CClp_opt failed\n");
    }

    if (lpcount && !silent) {
        printf ("Average Clique Value: %f\n", meanval / ((double) lpcount));
        fflush (stdout);
    }

    for (i = nwant - 1, k = 0; i >= 0; i--) {
        if (slist[i].name != -1) {
            k++;
            if (!silent) {
                printf ("Top Clique: %d  (%f)\n", slist[i].name, slist[i].val);
                fflush (stdout);
            }
        }
    }

    if (k == 0) {
        if (!silent) {
            printf ("WARNING: no branching cliques were found\n");
            fflush (stdout);
        }
        goto CLEANUP;
    }

    *bcliques = CC_SAFE_MALLOC (k, CCtsp_lpclique);
    if (!(*bcliques)) {
        fprintf (stderr, "out of memory in CCtsp_find_branch_cliques\n");
        rval = 1; goto CLEANUP;
    }
    if (bval) {
        *bval = CC_SAFE_MALLOC (k, double);
        if (!(*bval)) {
            fprintf (stderr, "out of memory in CCtsp_find_branch_cliques\n");
            CC_IFFREE (*bcliques, CCtsp_lpclique);
            rval = 1; goto CLEANUP;
        }
    }
    *ngot = k;
    for (i = 0, k = 0; i < nwant; i++) {
        if (slist[i].name != -1) {
            rval = CCtsp_copy_lpclique (&cliques[slist[i].name],
                                        &(*bcliques)[k]);
            if (rval) {
                fprintf (stderr, "CCtsp_copy_clique failed\n");
                for (i = 0; i < k; i++) {
                    CC_IFFREE ((*bcliques)[i].nodes, CCtsp_segment);
                }
                CC_FREE (*bcliques, CCtsp_lpclique);
                goto CLEANUP;
            }
            if (bval) {
                (*bval)[k] = slist[i].val;
            }
            k++;
        }
    }
    *ngot = k;


CLEANUP:

    CC_IFFREE (slist, sbitem);
    for (i = 0; i < ccount; i++) {
        CC_IFFREE (cliques[i].nodes, CCtsp_segment);
    }
    CC_IFFREE (cliques, CCtsp_lpclique);
    CC_IFFREE (lpcand, int);
    CClp_free_warmstart (&warmstart);

    return 0;
}

#define TSP_STRONG_GETWEIGHT_BATCH 250

static int find_candidate_cliques (CCtsp_lp *lp, int nwant, int *ngot,
        CCtsp_lpclique **list, int use_getweight, int silent)
{
    int      xcount   = 0;
    double   *x       = (double *) NULL;
    int      *xlist   = (int *) NULL;
    int      ccount   = 0;
    double   *cval    = (double *) NULL;
    CCtsp_lpclique *cliques = (CCtsp_lpclique *) NULL;
    double   *weights = (double *) NULL;
    sbitem   *slist   = (sbitem *) NULL;
    CCtsp_lprow cr;
    CCtsp_lpcut_in cu;
    int i, k, j, batch, ntry, nzlist;
    double down, up, t, sval;
    int rval = 0;

    CCtsp_init_lprow (&cr);

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL, &xcount,
                     &xlist, &x, (double **) NULL, (double **) NULL,
                     (double **) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_get_lp_result failed\n"); goto CLEANUP;
    }

    ntry = (use_getweight ? (TSP_STRONG_CUT_CANDIDATES +
                            (TSP_STRONG_CUT_CANDIDATES_MULT * nwant)) : nwant);

    rval = CCtsp_branch_cutpool_cliques (lp->pool, &cliques, &ccount,
              lp->graph.ncount, xcount, xlist, x, ntry, &cval, silent);
    if (rval) {
        fprintf (stderr, "CCtsp_branch_cutpool_cliques failed\n");
        goto CLEANUP;
    }
    if (ccount == 0) {
        fprintf (stderr, "WARNING: no cutpool cliques were found\n");
        goto CLEANUP;
    }

    if (use_getweight) {
        CCtsp_init_lpcut_in (&cu);
        cu.cliquecount = 1;
        cu.rhs  = 2;

        weights = CC_SAFE_MALLOC (TSP_STRONG_GETWEIGHT_BATCH, double);
        if (!weights) {
            fprintf (stderr, "out of memory in find_candidate_cliques\n");
            rval = 1; goto CLEANUP;
        }

        slist = CC_SAFE_MALLOC (nwant + 1, sbitem);
        if (!slist) {
            fprintf (stderr, "out of memory in find_candidate_cliques\n");
            rval = 1; goto CLEANUP;
        }
        init_sblist (slist, nwant);

        for (i = 0, batch = 0; i < ccount; i++) {
            cu.cliques = &(cliques[i]);
            nzlist = CCtsp_lpcut_in_nzlist (&lp->graph, &cu);
            rval = CCtsp_add_nzlist_to_lp (lp, nzlist, 2, 'G', &cr);
            if (rval) {
                fprintf (stderr, "CCtsp_add_nzlist_to_lp failed\n");
                goto CLEANUP;
            }
            if ((i+1) % TSP_STRONG_GETWEIGHT_BATCH == 0 || i == ccount - 1) {
                rval = CCutil_reallocrus_count ((void **) &(cr.begin),
                            cr.rowcnt + 1, sizeof (int));
                if (rval) {
                    fprintf (stderr, "out of memory\n");
                    goto CLEANUP;
                }
                cr.begin[cr.rowcnt] = cr.nzcnt;
                rval = CClp_getweight (lp->lp, cr.rowcnt, cr.begin,
                                       cr.indices, cr.entries, weights);
                if (rval) {
                    fprintf (stderr, "CClp_getweight failed\n"); goto CLEANUP;
                }

                for (j = 0; j < cr.rowcnt; j++) {
                    t = cval[j] - 2.0;
                    down = (t * t) /  weights[j];
                    t = 4.0 - cval[j];
                    up = (t * t) / weights[j];
                    sval = TSP_BRANCH_STRONG_CUT_NORM_VAL (down, up);
                    insert_sblist (slist, sval, batch+j);
                }
                CCtsp_free_lprow (&cr);
                batch = i+1;
            }
        }

        for (i = 0, k = 0; i < nwant; i++) {
            if (slist[i].name != -1) {
                k++;
            }
        }
        if (k == 0) {
            if (!silent) {
                printf ("WARNING: no candidate branching cliques were found\n");
                fflush (stdout);
            }
            goto CLEANUP;
        }

        *list = CC_SAFE_MALLOC (k, CCtsp_lpclique);
        if (!(*list)) {
            fprintf (stderr, "out of memory in CCtsp_find_branch_cliques\n");
            rval = 1; goto CLEANUP;
        }
        for (i = 0, k = 0; i < nwant; i++) {
            if (slist[i].name != -1) {
                rval = CCtsp_copy_lpclique (&cliques[slist[i].name],
                                            &(*list)[k++]);
                if (rval) {
                    fprintf (stderr, "CCtsp_copy_clique failed\n");
                    for (i = 0; i < k; i++) {
                        CC_IFFREE ((*list)[i].nodes, CCtsp_segment);
                    }
                    CC_FREE (*list, CCtsp_lpclique);
                    goto CLEANUP;
                }
            }
         }
        *ngot = k;
    } else {
        *list = cliques;
        *ngot = ccount;
    }

CLEANUP:

    CCtsp_free_lprow (&cr);
    CC_IFFREE (slist, sbitem);
    CC_IFFREE (x, double);
    CC_IFFREE (xlist, int);
    CC_IFFREE (cval, double);
    CC_IFFREE (weights, double);
    if (use_getweight) {
        for (i = 0; i < ccount; i++) {
            CC_IFFREE (cliques[i].nodes, CCtsp_segment);
        }
        CC_IFFREE (cliques, CCtsp_lpclique);
    }

    return rval;
}

static int test_cut_branch (CCtsp_lp *lp, CCtsp_lpclique *c, double *down,
        double *up, int iter, int silent, CClp_warmstart *warmstart)
{
    CCtsp_lprow cr;
    CCtsp_lpcut_in cu;
    int nzlist, status;
    int rval = 0;

    /* Note: this function does not reoptimize the last LP - so after a */
    /* series of calls, the calling program should call optimize.       */

    *down = -CCtsp_LP_MAXDOUBLE;
    *up   = -CCtsp_LP_MAXDOUBLE;
    CCtsp_init_lprow (&cr);

    CCtsp_init_lpcut_in (&cu);
    cu.cliquecount = 1;
    cu.cliques = c;
    cu.rhs  = 2;

    nzlist = CCtsp_lpcut_in_nzlist (&lp->graph, &cu);
    rval = CCtsp_add_nzlist_to_lp (lp, nzlist, 2, 'L', &cr);
    if (rval) {
        fprintf (stderr, "CCtsp_add_nzlist_to_lp failed\n"); goto CLEANUP;
    }
    rval = CCtsp_add_multiple_rows (lp, &cr);
    if (rval) {
        fprintf (stderr, "CCtsp_add_multiple_rows failed\n"); goto CLEANUP;
    }
    CCutil_start_timer (&lp->stats.strongbranch_opt);
    
    rval = CClp_limited_dualopt (lp->lp, iter, &status, &lp->upperbound);
    if (rval) {
        fprintf (stderr, "CClp_limited_dualopt failed\n"); goto CLEANUP;
    }
    CCutil_stop_timer (&lp->stats.strongbranch_opt, 0);
    if (status == CClp_INFEASIBLE) {
        if (!silent) {
            printf ("Down side of cut branch is infeasible\n");
            fflush (stdout);
        }
        *down = lp->upperbound;
    } else if (status == CClp_UNKNOWN) {
        if (!silent) {
            printf ("Down side information is not available\n");
            fflush (stdout);
        }
        *down = lp->lowerbound;
    } else {
        rval = CClp_objval (lp->lp, down);
        if (rval) {
            fprintf (stderr, "CClp_objval failed\n"); goto CLEANUP;
        }
        if (*down > lp->upperbound) *down = lp->upperbound;
    }

    rval = CCtsp_delete_cut (lp, lp->cuts.cutcount);
    if (rval) {
        fprintf (stderr, "CCtsp_delete_cut failed\n"); goto CLEANUP;
    }
    rval = CClp_load_warmstart (lp->lp, warmstart);
    if (rval) {
        fprintf (stderr, "CClp_load_warmstart failed\n"); goto CLEANUP;
    }

    cr.sense[0] = 'G';
    cr.rhs[0] = 4.0;
    rval = CCtsp_add_multiple_rows (lp, &cr);
    if (rval) {
        fprintf (stderr, "CCtsp_add_multiple_rows failed\n"); goto CLEANUP;
    }
    CCutil_start_timer (&lp->stats.strongbranch_opt);
    rval = CClp_limited_dualopt (lp->lp, iter, &status, &lp->upperbound);
    if (rval) {
        fprintf (stderr, "CClp_limited_dualopt failed\n"); goto CLEANUP;
    }
    CCutil_stop_timer (&lp->stats.strongbranch_opt, 0);
    if (status == CClp_INFEASIBLE) {
        if (!silent) {
            printf ("Up side of cut branch is infeasible\n"); fflush (stdout);
        }
        *up = lp->upperbound;
    } else if (status == CClp_UNKNOWN) {
        if (!silent) {
            printf ("Up side information is not available\n"); fflush (stdout);
        }
        *up = lp->lowerbound;
    } else {
        rval = CClp_objval (lp->lp, up);
        if (rval) {
            fprintf (stderr, "CClp_objval failed\n"); goto CLEANUP;
        }
        if (*up > lp->upperbound) *up = lp->upperbound;
    }

    rval = CCtsp_delete_cut (lp, lp->cuts.cutcount);
    if (rval) {
        fprintf (stderr, "CCtsp_delete_cut failed\n"); goto CLEANUP;
    }
    rval = CClp_load_warmstart (lp->lp, warmstart);
    if (rval) {
        fprintf (stderr, "CClp_load_warmstart failed\n"); goto CLEANUP;
    }

CLEANUP:

    CCtsp_free_lprow (&cr);
    return rval;
}

int CCtsp_execute_branch (CCtsp_lp *lp, CCtsp_branchobj *b, int silent,
        CCrandstate *rstate)
{
    int n0      = -1;
    int n1      = -1;
    CCtsp_lpclique *c = (CCtsp_lpclique *) NULL;
    int rval    = 0;
    int i, j;

    if (!b) {
        fprintf (stderr, "CCtsp_execute_branch called without a CCtsp_branchobj\n");
        rval = 1; goto CLEANUP;
    }

    if (b->ends[0] != -1) {
        n0 = b->ends[0];
        n1 = b->ends[1];
        if (!silent) {
            printf ("Branch Edge (%d,%d), to value %d\n", n0, n1, b->rhs);
            fflush (stdout);
        }

        if (n0 >= lp->graph.ncount || n0 < 0 ||
            n1 >= lp->graph.ncount || n1 < 0) {
            fprintf (stderr, "CCtsp_execute_branch has invalid nodes\n");
            rval = 1; goto CLEANUP;
        }

        if (n0 > n1) {
            CC_SWAP (n0, n1, j);
        }

        j = CCtsp_find_edge (&lp->graph, n0, n1);
        if  (j < 0) {
            fprintf (stderr, "branching edge is not in the LP edgeset\n");
            rval = 1; goto CLEANUP;
        }
        if (lp->graph.edges[j].fixed) {
            fprintf (stderr, "branching edge is fixed to 1 in the LP\n");
            rval = 1; goto CLEANUP;
        }
        if (lp->graph.edges[j].branch) {
            fprintf (stderr, "branching edge has already been branched\n");
            rval = 1; goto CLEANUP;
        }

        if (b->rhs) {
            rval = CClp_setbnd (lp->lp, j, 'L', 1.0);
            if (rval) {
                fprintf (stderr, "CClp_setbnd failed\n");
                rval = 1; goto CLEANUP;
            }
            lp->graph.edges[j].branch = lp->branchdepth + 1;
        } else {
            rval = CClp_setbnd (lp->lp, j, 'U', 0.0);
            if (rval) {
                fprintf (stderr, "CClp_setbnd failed\n");
                rval = 1; goto CLEANUP;
            }
            lp->graph.edges[j].branch = -(lp->branchdepth + 1);
        }
    } else {
        CCtsp_lprow cr;
        CCtsp_lpcut_in d;

        if (!b->clique) {
            fprintf (stderr, "CCtsp_branchobj has no edge or clique\n");
            rval = 1; goto CLEANUP;
        }

        if (!silent) {
            printf ("Branch Clique "); fflush (stdout);
            for (i = 0; i < b->clique->segcount; i++) {
                printf ("%d->%d ", b->clique->nodes[i].lo,
                                   b->clique->nodes[i].hi);
                fflush (stdout);
            }
            if (b->sense == 'G') {
                printf ("to at least %d\n", b->rhs);
            } else {
                printf ("to at most %d\n", b->rhs);
            }
            fflush (stdout);
        }

        c = CC_SAFE_MALLOC (1, CCtsp_lpclique);
        if (!c) {
            fprintf (stderr, "out of memory in CCtsp_execute_branch\n");
            rval = 1; goto CLEANUP;
        }
        rval = CCtsp_copy_lpclique (b->clique, c);
        if (rval) {
            fprintf (stderr, "CCtsp_copy_lpclique failed\n");
            rval = 1; goto CLEANUP;
        }

        CCtsp_init_lpcut_in (&d);
        d.cliquecount = 1;
        d.rhs = b->rhs;
        d.sense = b->sense;
        d.branch = 1;
        d.cliques = c;

        rval = CCtsp_construct_skeleton (&d, lp->graph.ncount);
        if (rval) {
            fprintf (stderr, "CCtsp_construct_skeleton failed\n");
            rval = 1; goto CLEANUP;
        }
        
        CCtsp_init_lprow (&cr);
        rval = CCtsp_add_cut (lp, &d, &cr);
        if (rval) {
            fprintf (stderr, "CCtsp_add_cut failed\n");
            rval = 1; goto CLEANUP;
        }
        rval = CCtsp_add_multiple_rows (lp, &cr);
        if (rval) {
            fprintf (stderr, "CCtsp_add_multiple_rows failed\n");
            rval = 1; goto CLEANUP;
        }
        CCtsp_free_lprow (&cr);
        CCtsp_free_lpcut_in (&d);
    }

    rval = CClp_opt (lp->lp, CClp_METHOD_DUAL);
    if (rval == 2) {
        rval = CCtsp_infeas_recover (lp, silent, rstate);
        if (rval == 2) {
            int tval;
            if (!silent) {
                printf ("Problem is really infeasible (CCtsp_execute_branch)\n");
                fflush (stdout);
            }

/* Bico Added on February 3, 1998 */
            tval = CCtsp_update_result (lp);
            if (tval) {
                fprintf (stderr, "CCtsp_update_result failed - ignoring\n");
                /* rval = 1; goto CLEANUP;  Bico Deleted on October 10, 2003 */
            }
            CCtsp_free_bigdual (&lp->exact_dual);
/* End Bico */

            goto CLEANUP;
        } else if (rval) {
            fprintf (stderr, "CCtsp_infeas_recover failed\n");
            rval = 1; goto CLEANUP;
        }
    } else if (rval) {
        fprintf (stderr, "CClp_opt failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_update_result (lp);
    if (rval) {
        fprintf (stderr, "CCtsp_update_result failed\n");
        rval = 1; goto CLEANUP;
    }
    CCtsp_free_bigdual (&lp->exact_dual);

CLEANUP:

    if (rval == 0 || rval == 2) {
        int sval = 0;
        sval = CCutil_reallocrus_count ((void **) &(lp->branchhistory),
                   lp->branchdepth + 1, sizeof (CCtsp_branchobj));
        if (sval) {
            fprintf (stderr, "CCutil_reallocrus_count failed\n"); return 1;
        }
        CCtsp_init_branchobj (&lp->branchhistory[lp->branchdepth]);
        lp->branchhistory[lp->branchdepth].depth   = lp->branchdepth + 1;
        lp->branchhistory[lp->branchdepth].ends[0] = n0;
        lp->branchhistory[lp->branchdepth].ends[1] = n1;
        lp->branchhistory[lp->branchdepth].rhs     = b->rhs;
        if (b->clique) {
            c = CC_SAFE_MALLOC (1, CCtsp_lpclique);
            if (!c) {
                fprintf (stderr, "out of memory in CCtsp_execute_branch\n");
                return 1;
            }
            sval = CCtsp_copy_lpclique (b->clique, c);
            if (sval) {
                fprintf (stderr, "CCtsp_copy_lpclique failed\n"); return 1;
            }
            lp->branchhistory[lp->branchdepth].clique = c;
        } else {
            lp->branchhistory[lp->branchdepth].clique =
                (CCtsp_lpclique *) NULL;
        }
        lp->branchhistory[lp->branchdepth].sense = b->sense;
        lp->branchdepth++;
    }
    return rval;
}

int CCtsp_execute_unbranch (CCtsp_lp *lp, CClp_warmstart *warmstart,
        int silent, CCrandstate *rstate)
{
    int rval = 0;
    int n0, n1;
    int num;
    int depth = lp->branchdepth;
    CCtsp_branchobj *b;
    int j;

    if (depth <= 0) {
        fprintf (stderr, "CCtsp_execute_unbranch called at depth 0\n");
        rval = 1; goto CLEANUP;
    }

    if (lp->branchhistory[depth - 1].depth != depth) {
        fprintf (stderr, "branchhistory is corrupted\n");
        rval = 1; goto CLEANUP;
    }
    b = &lp->branchhistory[depth - 1];

    if (lp->branchhistory[depth - 1].ends[0] != -1) {
        n0    = b->ends[0];
        n1    = b->ends[1];
        if (!silent) {
            printf ("Unbranch Edge (%d,%d), from value %d\n", n0, n1, b->rhs);
            fflush (stdout);
        }

        if (n0 > n1) {
            CC_SWAP (n0, n1, j);
        }

        j = CCtsp_find_edge (&lp->graph, n0, n1);
        if  (j < 0) {
            fprintf (stderr, "ERROR: unbranching 1-edge is not in LP\n");
            rval = 1; goto CLEANUP;
        }
        if (b->rhs) {
            if (lp->graph.edges[j].branch <= 0) {
                fprintf (stderr, "unbranching 1-edge not branched to 1\n");
                rval = 1; goto CLEANUP;
            }
            rval = CClp_setbnd (lp->lp, j, 'L', 0.0);
            if (rval) {
                fprintf (stderr, "CClp_setbnd failed\n");
                rval = 1; goto CLEANUP;
            }
        } else {
            if (lp->graph.edges[j].branch >= 0) {
                fprintf (stderr, "unbranching 0-edge not branched to 0\n");
                rval = 1; goto CLEANUP;
            }

            rval = CClp_setbnd (lp->lp, j, 'U', 1.0);
            if (rval) {
                fprintf (stderr, "CClp_setbnd failed\n");
                rval = 1; goto CLEANUP;
            }
        }
        lp->graph.edges[j].branch = 0;
    } else {
        if (!b->clique) {
            fprintf (stderr, "branchhistory has no edge or clique\n");
            rval = 1; goto CLEANUP;
        }
        rval = find_branched_clique (lp, b->clique, b->sense, b->rhs, &num);
        if (rval) {
            fprintf (stderr, "find_branched_clique failed\n");
            rval = 1; goto CLEANUP;
        }
        if (!silent) {
            printf ("The unbranching clique is cut %d\n", num); fflush (stdout);
         }
        if (lp->cuts.cuts[num].branch == 0) {
            fprintf (stderr, "the unbranching clique is not set to branch\n");
            rval = 1; goto CLEANUP;
        }

        if (!silent) {
            int q;
            CCtsp_lpcut *cu = &lp->cuts.cuts[num];
            CCtsp_lpclique *t;

            printf ("Sense: %c  RHS: %d  Cliques: %d  Branch: %d\n",
                 cu->sense, cu->rhs, cu->cliquecount, cu->branch);
            t = &lp->cuts.cliques[cu->cliques[0]];
            printf ("Clique: ");
            for (q = 0; q < t->segcount; q++) {
                printf ("%d->%d ", t->nodes[q].lo, t->nodes[q].hi);
            }
            printf ("\n"); fflush (stdout);
        }

        rval = CCtsp_delete_cut (lp, num);
        if (rval) {
            fprintf (stderr, "CCtsp_delete_cut failed\n");
            rval = 1; goto CLEANUP;
        }
        CCtsp_delete_cut_from_cutlist (&lp->cuts, num);
    }

    if (warmstart) {
        rval = CClp_load_warmstart (lp->lp, warmstart);
        if (rval) {
            fprintf (stderr, "CClp_load_warmstart failed\n");
            rval = 1; goto CLEANUP;
        }
    }

    rval = CClp_opt (lp->lp, CClp_METHOD_DUAL);

/* Bico Added on February 3, 1998 */
    if (rval == 2) {
        fprintf (stderr, "infeasible lp in CCtsp_execute_unbranch\n");
        rval = CCtsp_infeas_recover (lp, silent, rstate);
        if (rval == 2) {
            int tval;
            if (!silent) {
                printf ("Problem is really infeasible (CCtsp_execute_unbranch)\n");
                fflush (stdout);
            }
            tval = CCtsp_update_result (lp);
            if (tval) {
                fprintf (stderr, "CCtsp_update_result failed\n");
                rval = 1; goto CLEANUP;
            }
            CCtsp_free_bigdual (&lp->exact_dual);
            goto CLEANUP;
        } else if (rval) {
            fprintf (stderr, "CCtsp_infeas_recover failed\n");
            rval = 1; goto CLEANUP;
        }
/* End Bico */
    } else if (rval) {
        fprintf (stderr, "CClp_opt failed\n"); goto CLEANUP;
    }

    rval = CCtsp_update_result (lp);
    if (rval) {
        fprintf (stderr, "CCtsp_update_result failed\n");
        rval = 1;  goto CLEANUP;
    }
    CCtsp_free_bigdual (&lp->exact_dual);

CLEANUP:

    if (!rval || rval == 2) {
        CCtsp_free_branchobj (&lp->branchhistory[lp->branchdepth - 1]);
        lp->branchdepth--;
    }
    return rval;
}

int CCtsp_add_branchhistory_to_lp (CCtsp_lp *lp)
{
    int i, k, num;
    int rval = 0;
    CCtsp_branchobj *b;

    for (i = 0; i < lp->branchdepth; i++) {
        b = &lp->branchhistory[i];
        if (b->ends[0] != -1) {
            if (lp->graph.ecount == 0) {
                printf ("No graph - can't add edge %d,%d = %d from branch history\n",
                        b->ends[0], b->ends[1], b->rhs);
            } else {
                k = CCtsp_find_edge (&lp->graph, b->ends[0], b->ends[1]);
                if (k == -1) {
                    fprintf (stderr, "edge in branch history is not in LP\n");
                    rval = 1; goto CLEANUP;
                }
                if (lp->graph.edges[k].fixed || lp->graph.edges[k].branch) {
                    fprintf (stderr, "edge in branch history is fixed/branched\n");
                    rval = 1; goto CLEANUP;
                }
                if (b->rhs) {
                    rval = CClp_setbnd (lp->lp, k, 'L', 1.0);
                    if (rval) {
                        fprintf (stderr, "CClp_setbnd failed\n"); goto CLEANUP;
                    }
                    lp->graph.edges[k].branch = b->depth;
                } else {
                    rval = CClp_setbnd (lp->lp, k, 'U', 0.0);
                    if (rval) {
                        fprintf (stderr, "CClp_setbnd failed\n"); goto CLEANUP;
                    }
                    lp->graph.edges[k].branch = -(b->depth);
                }
            }
        } else {
            rval = find_branched_clique (lp, b->clique, b->sense,
                                             b->rhs, &num);
            if (rval) {
                fprintf (stderr, "find_branch_clique failed\n");
                goto CLEANUP;
            }
            lp->cuts.cuts[num].branch = 1;
        }
    }

CLEANUP:

    return rval;
}

static int find_branched_clique (CCtsp_lp *lp, CCtsp_lpclique *c, char sense,
                                 int rhs, int *cutnum)
{
    int i;
    CCtsp_lpcut *cu;
    CCtsp_lpcut *cuts       = lp->cuts.cuts;
    CCtsp_lpclique *cliques = lp->cuts.cliques;
    int cutcount      = lp->cuts.cutcount;
    int diff = 0;

    *cutnum = -1;

    for (i = 0; i < cutcount; i++) {
        cu = &cuts[i];
        if (cu->cliquecount == 1 &&
            cu->sense == sense && cu->rhs == rhs) {
            CCtsp_lpclique_compare (&cliques[cu->cliques[0]], c, &diff);
            if (!diff) {
                if (*cutnum == -1) {
                    *cutnum = i;
                } else {
                    fprintf (stderr, "two copies of branched clique\n");
                    return 1;
                }
            }
        }
    }

    if (*cutnum == -1) {
        fprintf (stderr, "did not find branched clique\n");
        return 1;
    } else {
        return 0;
    }
}

int CCtsp_bb_find_branch (char *probname, int probnum, int ncount,
        CCdatagroup *dat, int *ptour, double *upperbound, CCtsp_lpcuts *pool,
        int nwant, int *ngot, CCtsp_branchobj **b, int usecliques,
        int longedge_branching, int *prune, int *foundtour, int *besttour,
        int silent, CCrandstate *rstate)
{
    int rval = 0;
    double tval;
    CCtsp_lp *lp = (CCtsp_lp *) NULL;
    int *cyc = (int *) NULL;
    int i, test;

    *foundtour = 0;
    *prune = 0;
    *b = (CCtsp_branchobj *) NULL;

    rval = CCtsp_bb_init_lp (&lp, probname, probnum, ncount, dat, ptour,
               *upperbound, pool, silent, rstate);
    if (rval) {
        fprintf (stderr, "CCtsp_bb_init_lp failed\n"); goto CLEANUP;
    }

    if (lp->lowerbound >= *upperbound - 0.9) {
        printf ("Do not branch, the lp is within 1.0 of the upperbound\n");
        fflush (stdout);

        rval = CCtsp_verify_lp_prune (lp, &test, silent);
        if (rval) {
            fprintf (stderr, "CCtsp_verify_lp_prune failed\n");
            goto CLEANUP;
        }
        if (test) {
            if (!silent) {
                printf ("verified that LP can be pruned\n");
                fflush (stdout);
            }
            *prune = 1;
            goto CLEANUP;
        } else {
            printf ("exact pricing could not prune search - need to branch\n");
            fflush (stdout);
        }
    }

    rval = CCtsp_find_branch (lp, nwant, ngot, b, &tval, &cyc, usecliques,
                              longedge_branching, silent);
    if (rval) {
        fprintf (stderr, "CCtsp_find_branch failed\n");
        goto CLEANUP;
    }

    if (*ngot == 0) {
        if (!silent) {
            printf ("No branch, found tour of value %.2f\n", tval);
            fflush (stdout);
        }
        if (tval < lp->upperbound) lp->upperbound = tval;
        rval = CCtsp_verify_lp_prune (lp, &test, silent);
        if (rval) {
            fprintf (stderr, "CCtsp_verify_lp_prune failed\n");
            goto CLEANUP;
        }
        if (test) {
            if (!silent) {
                printf ("verified that LP can now be pruned\n");
                fflush (stdout);
            }
            *foundtour = 1;
            if (tval < *upperbound) {
                *upperbound = tval;
                if (besttour) {
                    for (i = 0; i < ncount; i++) {
                        besttour[i] = cyc[i];
                    }
                }
            }
            goto CLEANUP;
        } else {
            fprintf (stderr, "new tour did not permit exact pruning\n");
            rval = 1; goto CLEANUP;
        }
    } else {
        if (!silent) {
            printf ("found branch\n"); fflush (stdout);
        }
    }


CLEANUP:

    if (lp) CCtsp_free_tsp_lp_struct (&lp);
    CC_IFFREE (cyc, int);
    return rval;
}

int CCtsp_bb_splitprob (char *probname, int probnum, int ncount,
        CCdatagroup *dat, int *ptour, double initial_ub, CCtsp_lpcuts *pool,
        CCtsp_branchobj *b, int child0, int child1, double *val0, double *val1,
        int *prune0, int *prune1, int silent, CCrandstate *rstate)
{
    int rval = 0;
    CCtsp_lp *lp = (CCtsp_lp *) NULL;

    *val0 = 0.0;
    *val1 = 0.0;
    *prune0 = 0;
    *prune1 = 0;

    rval = CCtsp_bb_init_lp (&lp, probname, probnum, ncount, dat, ptour,
               initial_ub, pool, silent, rstate);
    if (rval) {
        fprintf (stderr, "CCtsp_bb_init_lp failed\n"); goto CLEANUP;
    }

    rval = branch_side (lp, b, 0, child0, val0, prune0, silent, rstate);
    if (rval) {
        fprintf (stderr, "branch_side failed\n"); goto CLEANUP;
    }

    CCtsp_free_tsp_lp_struct (&lp);
    rval = CCtsp_bb_init_lp (&lp, probname, probnum, ncount, dat, ptour,
               initial_ub, pool, silent, rstate);
    if (rval) {
        fprintf (stderr, "CCtsp_bb_init_lp failed\n"); goto CLEANUP;
    }

    rval = branch_side (lp, b, 1, child1, val1, prune1, silent, rstate);
    if (rval) {
        fprintf (stderr, "branch_side failed\n"); goto CLEANUP;
    }


CLEANUP:

    if (lp) CCtsp_free_tsp_lp_struct (&lp);
    return rval;
}

static int branch_side (CCtsp_lp *lp, CCtsp_branchobj *b, int side, int child,
        double *val, int *prune, int silent, CCrandstate *rstate)
{
    int rval = 0;
    int oldid       =  lp->id;
    int oldparent   =  lp->parent_id;
    double oldbound = lp->lowerbound;
    double newbound;
    int test;

    *val = 0.0;
    *prune = 0;

    if (b->ends[0] != -1) {
        if (!silent) {
            printf ("Creating child %d of LP %d: Set Edge (%d, %d) to %d\n",
                         side, lp->id, b->ends[0], b->ends[1], side);
            fflush (stdout);
        }
        if (side == 0) {
            b->rhs = 0;
        } else {
            b->rhs = 1;
        }
    } else {
        if (side == 0) {
            if (!silent) {
                printf ("Creating child 0 of LP %d: Set Clique <= 2\n", lp->id);
                fflush (stdout);
            }
            b->rhs = 2; b->sense = 'L';
        } else {
            if (!silent) {
                printf ("Creating child 1 of LP %d: Set Clique >= 4\n", lp->id);
                fflush (stdout);
            }
            b->rhs = 4; b->sense = 'G';
        }
    }
    fflush (stdout);

    rval = CCtsp_execute_branch (lp, b, silent, rstate);
    if (rval && rval != 2) {
        fprintf (stderr, "CCtsp_execute_branch failed\n"); goto CLEANUP;
    } else if (rval == 2) {
        printf ("Branched-LP is infeasible\n"); fflush (stdout);
        rval = CCtsp_verify_infeasible_lp (lp, &test, silent);
        if (rval) {
            fprintf (stderr, "CCtsp_verify_infeasible_lp failed\n");
            goto CLEANUP;
        }
        if (test) {
            if (!silent) {
                printf ("Creating child leafnode - infeasible\n");
                fflush (stdout);
            }
            *val = CCtsp_LP_MAXDOUBLE;
            *prune = 1;
            lp->parent_id = oldid;
            lp->id = child;
            rval = CCtsp_write_probleaf_id (lp);
            if (rval) {
                fprintf (stderr, "CCtsp_write_probleaf_id failed\n");
                goto CLEANUP;
            }
            rval = 0;
        } else {
            fprintf (stderr, "did not verify an infeasible LP\n");
            rval = 1; goto CLEANUP;
        }
    } else {
        rval = CCtsp_pricing_loop (lp, &newbound, silent, rstate);
        CCcheck_rval (rval, "CCtsp_pricing_loop failed");

        *val = newbound;
        lp->lowerbound = newbound;

        if (lp->lowerbound >= lp->upperbound - 0.9) {
            rval = CCtsp_verify_lp_prune (lp, &test, silent);
            CCcheck_rval (rval, "CCtsp_verify_lp_prune failed");
            if (test) {
                if (!silent) {
                    printf ("verified that child can be pruned\n");
                    fflush (stdout);
                }
                *prune = 1;
                lp->parent_id = oldid;
                lp->id = child;
                rval = CCtsp_write_probleaf_id (lp);
                CCcheck_rval (rval,"CCtsp_write_probleaf_id failed");
            } else {
                printf ("exact pricing could not prune child\n");
                fflush (stdout);
            }
        }

        if (*prune == 0) {
            lp->parent_id = oldid;
            lp->id = child;
            rval = CCtsp_write_probfile_id (lp);
            CCcheck_rval (rval, "CCtsp_write_probfile_id failed");
            lp->parent_id = oldparent;
            lp->id = oldid;
        }

        lp->lowerbound = oldbound;
    }

CLEANUP:

    return rval;
}

int CCtsp_splitprob (CCtsp_lp *lp, CCtsp_branchobj *b, int child0, int child1,
        int silent, CCrandstate *rstate)
{
    int oldid     =  lp->id;
    int oldparent =  lp->parent_id;
    CClp_warmstart *warmstart = (CClp_warmstart *) NULL;
    int rval = 0;

    rval = CClp_get_warmstart (lp->lp, &warmstart);
    if (rval) {
        fprintf (stderr, "CClp_get_warmstart failed\n"); goto CLEANUP;
    }

    lp->parent_id = lp->id;

    if (b->ends[0] != -1) {
        b->rhs = 0;
    } else {
        b->rhs = 2;
        b->sense = 'L';
    }

    lp->id = child0;
    rval = CCtsp_execute_branch (lp, b, silent, rstate);
    if (rval == 2) {
        rval = 0;
        printf ("The down side of the branch was infeasible\n");
        fflush (stdout);
    } else if (rval) {
        fprintf (stderr, "CCtsp_execute_branch failed\n"); goto CLEANUP;
    }
    rval = CCtsp_write_probfile_id (lp);
    if (rval) {
        fprintf (stderr, "CCtsp_write_probfile_id failed\n"); goto CLEANUP;
    }
    rval = CCtsp_execute_unbranch (lp, warmstart, silent, rstate);
    if (rval) {
        fprintf (stderr, "CCtsp_execute_unbranch failed\n"); goto CLEANUP;
    }

    if (b->ends[0] != -1) {
        b->rhs = 1;
    } else {
        b->rhs = 4;
        b->sense = 'G';
    }

    lp->id = child1;
    rval = CCtsp_execute_branch (lp, b, silent, rstate);
    if (rval == 2) {
        rval = 0;
        printf ("The up side of the branch was infeasible\n");
        fflush (stdout);
    } else if (rval) {
        fprintf (stderr, "CCtsp_execute_branch failed\n"); goto CLEANUP;
    }
    rval = CCtsp_write_probfile_id (lp);
    if (rval) {
        fprintf (stderr, "CCtsp_write_probfile_id failed\n"); goto CLEANUP;
    }
    rval = CCtsp_execute_unbranch (lp, warmstart, silent, rstate);
    if (rval) {
        fprintf (stderr, "CCtsp_execute_unbranch failed\n"); goto CLEANUP;
    }


CLEANUP:

    CClp_free_warmstart (&warmstart);
    lp->parent_id = oldparent;
    lp->id = oldid;

    return rval;
}

int CCtsp_dumptour (int ncount, CCdatagroup *dat, int *perm, char *probname,
        int *tour, char *fname, int writeedges, int silent)
{
    int rval = 0;
    int *cyc = (int *) NULL;
    int i;
    double len = 0.0;
    FILE *fout = (FILE *) NULL;
    char buf[1024];

    if (!perm || !tour) {
        fprintf (stderr, "bad input for CCtsp_dumptour\n");
        rval = 1; goto CLEANUP;
    }

    if (!fname) {
        sprintf (buf, "%s.sol", probname);
    } else {
        sprintf (buf, "%s", fname);
    }

    cyc = CC_SAFE_MALLOC (ncount, int);
    if (!cyc) {
        fprintf (stderr, "out of memory in CCtsp_dumptour\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) {
        cyc[i] = 0;
    }
    for (i = 0; i < ncount; i++) {
        cyc[tour[i]] = 1;
    }
    for (i = 0; i < ncount; i++) {
        if (cyc[i] == 0) {
            fprintf (stderr, "array is not a tour in CCtsp_dumptour\n");
            rval = 1; goto CLEANUP;
        }
    }
    for (i = 0; i < ncount; i++) {
        cyc[i] = perm[tour[i]];
    }

    if (dat) {
        for (i = 1; i < ncount; i++) {
            len += (double) CCutil_dat_edgelen (tour[i-1], tour[i], dat);
        }
        len += (double) CCutil_dat_edgelen (tour[ncount-1], tour[0], dat);
        if (!silent) {
            printf ("Write tour of length %.2f to %s\n", len, buf);
            fflush (stdout);
        }
    } else {
        if (!silent) {
            printf ("Write tour to %s\n", buf); fflush (stdout);
        }
    }

    if (!writeedges) {
        rval = CCutil_writecycle (ncount, buf, cyc, 0);
        if (rval) {
            fprintf (stderr, "CCutil_writecycle failed\n"); goto CLEANUP;
        }
    } else {
        if (!dat) {
            fprintf (stderr, "need datagroup to write edge file\n");
            rval = 1; goto CLEANUP;
        } else {
            fout = fopen (buf, "w");
            if (fout == (FILE *) NULL) {
                perror (buf);
                fprintf (stderr, "Unable to open %s for output\n", buf);
                rval = 1; goto CLEANUP;
            }
            fprintf (fout, "%d %d\n", ncount, ncount);
            for (i = 1; i < ncount; i++) {
                fprintf (fout, "%d %d %d\n", cyc[i-1], cyc[i],
                         CCutil_dat_edgelen (tour[i-1], tour[i], dat));
            }
            fprintf (fout, "%d %d %d\n", cyc[ncount - 1], cyc[0],
                     CCutil_dat_edgelen (tour[ncount - 1], tour[0], dat));
        }
    }

CLEANUP:

    CC_IFFREE (cyc, int);
    return rval;
}

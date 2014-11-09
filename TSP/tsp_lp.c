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
/*            ROUTINES TO BUILD LPS AND CALL THE LP SOLVER                  */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: September 26, 1995                                                */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  void CCtsp_init_tsp_lpcuts_struct (CCtsp_lpcuts *c)                     */
/*    INITIALIZES the CCtsp_lpcuts struture with NULL values                */
/*                                                                          */
/*  void CCtsp_init_tsp_lp_struct (CCtsp_lp *lp)                            */
/*    INITIALIZES the CCtsp_lp struture with NULL values                    */
/*                                                                          */
/*  void CCtsp_free_tsp_lp_struct (CCtsp_lp **lp)                           */
/*     free a CCtsp_lp                                                      */
/*                                                                          */
/*  int CCtsp_init_lp (CCtsp_lp **lp, char *probloc, int probnum,           */
/*      char *probfilename, int ncount, CCdatagroup *dat, int ecount,       */
/*      int *elist, int *elen, int excount, int *exlist, int *exlen,        */
/*      int exvalid, int *ptour, double initial_ub,                         */
/*      CCtsp_lpcuts *pool, int silent, CCrandstate *rstate)                */
/*    BUILDS/READS the problem, and loads it into the LP solver. If         */
/*     probnum < 0, init_lp will build an initial problem according to      */
/*     elist and elen; otherwise it will read the problem from disk. If     */
/*     the problem is read from disk, then the elist is ignored.            */
/*     -lp is a handle to the tsp lp (filled in by init_lp)                 */
/*     -probname is the name for the problem                                */
/*     -probnum is the number for the problem                               */
/*     -ncount is the number of nodes                                       */
/*     -dat is a handle on the complete graph                               */
/*     -ecount, elist, elen specify an initial edge set; if a prob is       */
/*      read from a file, then this list is ignored                         */
/*     -excount, exlist, exlen specify an full edge set (they can be        */
/*      0, NULL, NULL); if the probfile already has an full edge set,       */
/*      then this values are ignored.                                       */
/*     -exvalid indicates whether or not the edges specified in exlist      */
/*      are a valid complete set of edges (0 no, 1 yes)                     */
/*     -pool is a pointer to a cutpool (can be NULL)                        */
/*     -silent will suppress print messages if set to a nonzero value       */
/*    NOTES: If init_lp returns 2, then the LP is infeasible (even after    */
/*     considering the full edge set).                                      */
/*                                                                          */
/*  int CCtsp_bb_init_lp (CCtsp_lp **lp, char *probname, int probnum,       */
/*      int ncount, CCdatagroup *dat, int *ptour, double initial_ub,        */
/*      CCtsp_lpcuts *pool, int silent, CCrandstate *rstate)                */
/*    SHORT form of CCtsp_init_lp for use in the branch and bound.          */
/*                                                                          */
/*  int CCtsp_get_lp_result (CCtsp_lp *lp, double *lb, double *ub,          */
/*      int *ecount, int **elist, double **x, double **rc,                  */
/*      double **node_pi, double **cut_pi)                                  */
/*    RETURNS a copy of the values cached in lp->result. However, it        */
/*     allows a single point of locking for the threaded version. Any       */
/*     return argument can be NULL.                                         */
/*     -lp is a pointer to the tsp lp                                       */
/*     -obj returns the location for the current objective value            */
/*     -ecount returns the location for the number of nonzero edges         */
/*     -elist returns the  location for the nonzero edges in end1 end2      */
/*      format                                                              */
/*     -x returns location for the edge values                              */
/*     -rc returns location for the edge values                             */
/*     -node_pi returns the values on the degree constraints                */
/*     -cut_pi returns the dual values on the cuts                          */
/*    NOTES: node_pi and cut_pi go to the LP to fetch the results.          */
/*                                                                          */
/*  int CCtsp_lpcut_in_nzlist (CCtsp_lpgraph *g, CCtsp_lpcut_in *c)         */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_add_cut (CCtsp_lp *lp, CCtsp_lpcut_in *d, CCtsp_lprow *cr)    */
/*    ADDS cut d to the lp structure and to cr (a call to                   */
/*     CCtsp_add_multiple will put the cut into the lp solver)              */
/*                                                                          */
/*  int CCtsp_add_nzlist_to_lp (CCtsp_lp *lp, int nzlist, int nrhs,         */
/*      char sense, CCtsp_lprow *cr)                                        */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_add_multiple_rows (CCtsp_lp *lp, CCtsp_lprow *cr)             */
/*    HANDS the cuts in cr to the lp solver.                                */
/*                                                                          */
/*  int CCtsp_delete_cut (CCtsp_lp *lp, int i)                              */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_write_probfile_sav (CCtsp_lp *lp)                             */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_write_probfile_id (CCtsp_lp *lp)                              */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_write_probroot_id (char *probloc, CCtsp_lp *lp)               */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCtsp_add_cuts_to_queue (CCtsp_lp *lp, CCtsp_lpcut_in **clist)     */
/*    ADDS clist to the queue of cuts to be processed by the lp solver;     */
/*     clist will be set to NULL                                            */
/*     -lp is a pointer to the tsp lp                                       */
/*     -clist is the head of a NULL terminated linked list of cuts          */
/*                                                                          */
/*  int CCtsp_process_cuts (CCtsp_lp *lp, int *pnadded, int tighten,        */
/*      int silent, CCrandstate *rstate)                                    */
/*     -lp is a pointer to the tsp lp                                       */
/*     -pnadded returns the location for the number of cuts added           */
/*     -tighten is a flag to indicate whether or not the tighten routine    */
/*      should be called for each cut before it is added to the LP          */
/*    NOTES: process_cuts runs through all the cuts in the queue;           */
/*      process_cuts also calls add_to_cutpool().  If process_cuts          */
/*      returns 2, then the LP is infeasible, even after considering        */
/*      the full edge set.                                                  */
/*                                                                          */
/*  int CCtsp_addbad_variables (CCtsp_lp *lp, CCtsp_edgegenerator *eg,      */
/*      double *ppenalty, int *pnadded, double rcthresh,                    */
/*      double maxpenalty, int phase1, int *feasible, int silent,           */
/*      CCrandstate *rstate)                                                */
/*    ADDS negative reduced cost edges to the LP; if phase1 is nonzero      */
/*     then the added edges attempt to make a feasible LP (in this          */
/*     case the eg variable is ignored and the edges are taked either       */
/*     from fulladj (if they are valid) or from dat)                        */
/*     -lp is a pointer to the tsp lp                                       */
/*     -eg is a generator for the edges to check                            */
/*     -ppenalty is the penalty from the last pass of pricing               */
/*     -pnadded is the number of negative reduced cost edges added          */
/*     -rcthresh is the threshold on the reduced cost of edges to be        */
/*      added (it should be something <= 0.0)                               */
/*     -maxpenalty is the maximum sum of penalties that is permitted        */
/*      before the rounds of pricing stop                                   */
/*     -phase1 should be 0 for normal column generation, and nonzero        */
/*      to try to fix an infeasible LP                                      */
/*     -feasible can be NULL, otherwise it is set to 1 if phase 1           */
/*      gets to a feasible LP and 0 if the LP really is infeasible          */
/*                                                                          */
/*  int CCtsp_eliminate_variables (CCtsp_lp *lp, int eliminate_sparse,      */
/*      int silent)                                                         */
/*    SETS edges to 0 or 1 if possible, based on reduced costs              */
/*     -lp is a pointer to the tsp lp                                       */
/*     -if eliminate_sparse is nonzero and the full edge list in lp is      */
/*      present and valid, eliminate variables from that list.  Otherwise   */
/*      eliminate from the complete graph                                   */
/*                                                                          */
/*  double CCtsp_cutprice (CCtsp_lpgraph *g, CCtsp_lpcut_in *c,             */
/*      double *x)                                                          */
/*    RETURNS the slack of cut c                                            */
/*     -g is a pointer to an CCtsp_lpgraph that matches the vector x        */
/*                                                                          */
/*  int CCtsp_add_vars_to_lp (CCtsp_lp *lp, CCtsp_predge *prlist, int n)    */
/*    ADDS the columns to the lp.                                           */
/*     -n the number of edges listed in prlist                              */
/*                                                                          */
/*  int CCtsp_update_result (CCtsp_lp *lp)                                  */
/*    UPDATES the solution information in the lp structure                  */
/*                                                                          */
/*  int CCtsp_infeas_recover (CCtsp_lp *lp, int silent,                     */
/*      CCrandstate *rstate)                                                */
/*    TRIES to add columns to lp to regain feasibiblity                     */
/*    NOTES: Returns 2 if the full lp is infeasible                         */
/*                                                                          */
/*  int CCtsp_build_lpgraph (CCtsp_lpgraph *g, int ncount, int ecount,      */
/*      int *elist, int *elen)                                              */
/*    BUILDS the node and edge lists for the CCtsp_lpgraph pointed to by    */
/*      g.                                                                  */
/*     -elen contains the edge lengths (it can be NULL, in which case       */
/*      the lengths are set to 0).                                          */
/*                                                                          */
/*  int CCtsp_build_lpadj (CCtsp_lpgraph *g, int estart, int eend)          */
/*    BUILDS the incidence list for the graph *g                            */
/*     -estart is the index of the first edge to include in the list        */
/*     -eend is the index of the last edge + 1                              */
/*                                                                          */
/*  int CCtsp_find_edge (CCtsp_lpgraph *g, int from, int to)                */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCtsp_init_lpgraph_struct (CCtsp_lpgraph *g)                       */
/*    INITIALIZES the CCtsp_lpgraph struct pointed to by g.                 */
/*                                                                          */
/*  void CCtsp_free_lpgraph (CCtsp_lpgraph *g)                              */
/*    FREES the fields in the CCtsp_lpgraph pointed to by g.                */
/*                                                                          */
/*  void CCtsp_init_lprow (CCtsp_lprow *cr)                                 */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCtsp_free_lprow (CCtsp_lprow *cr)                                 */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_inspect_full_edges (CCtsp_lp *lp)                             */
/*    CHECKS that full edge set contains the current LP edge set; it        */
/*     returns 0 if it is okay and 1 if some edge is not present            */
/*     -lp is the CCtsp_lp                                                  */
/*                                                                          */
/*  int CCtsp_lpcut_nzlist (CCtsp_lpgraph *g, CCtsp_lpcut *c,               */
/*      CCtsp_lpclique *cliques, CCtsp_lpdomino *dominos, int do_mods)      */
/*                                                                          */
/*  int CCtsp_resparsify_lp (CCtsp_lp *lp, int silent)                      */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_read_probfile (CCtsp_lp *lp, char *fname, char *probloc,      */
/*      int *ncount, int silent)                                            */
/*    READS a tsp file and loads the results into lp                        */
/*     -lp is an initialized lp (via a call to init_tsp_lp_struct; the      */
/*      results are returned in this struct                                 */
/*     -fname is the tsp file                                               */
/*     -probloc is the problem location.  If probloc == NULL, then the      */
/*      this is derived from fname.                                         */
/*     -ncount is a pointer to the number of nodes; if it is nonnull and    */
/*      *ncount != 0, it is used as a check to see if the tsp file          */
/*      matches.  If it is nonnull and *ncount == 0, it returns the         */
/*      from the tsp file.                                                  */
/*                                                                          */
/*  int CCtsp_read_probfile_id (CCtsp_lp *lp, char *fname, int id,          */
/*      int *ncount, int silent)                                            */
/*    READS a tsp file and loads the results into lp, where the filename    */
/*     is obtained by using the id.                                         */
/*                                                                          */
/*  int CCtsp_reduced_cost_nearest (CCtsp_lp *lp, int k, int *ecount,       */
/*      int **elist, double **elen, int sparse)                             */
/*    RETURNS the smallest k reduced cost edge graph in ecount, elist,      */
/*     elen.                                                                */
/*                                                                          */
/*  int CCtsp_dump_rc_nearest (CCtsp_lp *lp, int k, char *fname,            */
/*      int sparse)                                                         */
/*    WRITES the smallest k reduced cost edge graph to fname.               */
/*                                                                          */
/*  int CCtsp_dump_x (CCtsp_lp *lp, char *fname)                            */
/*    WRITES the lp solution to fname.                                      */
/*    NOTES: The vector contains the original node names.                   */
/*                                                                          */
/****************************************************************************/

/*
Remark:
  Sparified cut is viewed as cut + CCtsp_sparser*stars (ie, a positive value
  in the CCtsp_sparser means that the cut in the lp is the original cut + the
  star for that node
*/

#include "machdefs.h"
#include "util.h"
#include "macrorus.h"
#include "fmatch.h"
#include "edgegen.h"
#include "linkern.h"
#include "tsp.h"
#include "lp.h"
#include "bigguy.h"
#include "cut.h"
#include "pq.h"
#include "cuttree.h"
#include "verify.h"

/* PERTURB_LPS is helpful for fl3795 */
#define PERTURB_LPS

/*
   parameters for addbad_variables.  Addbad generates edges (x,y) with
   len(x,y) - pi_ub[x] - pi_ub[y] < 0 in chunks of PRICE_GEN, computes
   rc for that chunk, and adds those with rc(x,y) < CCtsp_PRICE_RCTHRESH to
   the pool.  This continues until the pool contains more than PRICE_POOL
   edges, at which point the worst PRICE_ADD edges are added to the lp, new
   reduced are computed for the pool, and only edges with rc(x,y) <
   CCtsp_PRICE_RCTHRESH are retained.  This overall process continues until
   the total penalty for a full pass through the edges without adding any is
   < CCtsp_PRICE_MAXPENALTY.
*/
/* Some of these should probably depend on the problem size */

#define REALLOC_FACTOR    1.25
#define PRICE_GEN         20000
#define PRICE_GEN_FACTOR  3
#define PRICE_POOL        1000     /* 10000 */
#define PRICE_ADD         100

#define PROBTYPE_NORMAL 0
#define PROBTYPE_ROOT   1
#define PROBTYPE_LEAF   2
#define PROBTYPE_SAVE   3


static int
    first_lp (CCtsp_lp *lp, char *probloc, int ncount, int ecount, int *elist,
        int *elen, int silent, CCrandstate *rstate),
    find_edge_full (CCtsp_lp *lp, int from, int to),
    load_lp (CCtsp_lp *lp, int silent),
    build_lp_cols (CCtsp_lpgraph *g, CCtsp_lpcuts *cuts, int estart, int eend,
        int *pnzcnt, double **pobj, int **pmatbeg, int **pmatcnt,
        int **pmatind, double **pmatval, double **plb, double **pub),
    checkout_cut (CCtsp_lp *lp, CCtsp_lpcut_in *c, double *x, CCtsp_lprow *cr,
        int tighten),
    update_newcuts (CCtsp_lp *lp, int silent, CCrandstate *rstate),
    phase1_test_edge (int end1, int end2, double *node_piest),
    phase1_generate_edges (CCtsp_lp *lp, double *node_piest, int nwant,
           int *ngot, int *genlist, int *genlen, int start, int *ni, int *nj,
           int *finished),
    pricing_duals (CCtsp_lp *lp, double *node_pi, double *node_piest,
           double *cut_pi, double *clique_pi, double *domino_pi),
    price_list (CCtsp_lp *lp, int ecount, CCtsp_predge *elist, double *node_pi,
           double *cut_pi, double *clique_pi, double *domino_pi, int phase1),
    age_cuts (CCtsp_lp *lp, int *ndeleted),
    age_edges (CCtsp_lp *lp, int *ndeleted),
    get_pi (CCtsp_lp *lp, double *node_pi, double *cut_pi),
    addrow_to_list (int nzcnt, double drhs, char sense,
           int *rmatind, double *rmatval, CCtsp_lprow *cr),
    lp_addcols (CCtsp_lp *lp, int ncols, int nzcnt,
           double *obj, int *matbeg, int *matind, double *matval,
           double *lb, double *ub),
    lp_delete_cut_set (CCtsp_lp *lp, int *del),
    lp_delete_var_set (CCtsp_lp *lp, int *del),
    sparse_nearest (CCtsp_lp *lp, int k, double *node_pi, double *cut_pi,
           double *clique_pi, double *domino_pi, CCutil_edgehash *eh),
    write_probfile (CCtsp_lp *lp, int probtype, char *fname),
    read_probfile (CCtsp_lp *lp, CCtsp_PROB_FILE *p, int *ncount, int silent);

static void
    lpcut_nonzero_work (CCtsp_lpgraph *g, CCtsp_lpclique *c, int *pnzlist),
    lpcut_nonzero_semicut (CCtsp_lpgraph *g, CCtsp_lpdomino *c, int *pnzlist),
    lpcut_nonzero_parity_handle (CCtsp_lpgraph *g, CCtsp_lpclique *c,
        int *pnzlist),
    lpcut_nonzero_domino (CCtsp_lpgraph *g, CCtsp_lpdomino *c, int *pnzlist),
    lpcut_nonzero_modify (CCtsp_lpgraph *g, int modcount, CCtsp_sparser *mods,
           int *pnzlist),
    clear_nzlist (CCtsp_lpgraph *g, int nzlist),
    pr_select (int nsel, int n, CCtsp_predge *list),
    add_cut_to_queue (CCtsp_lp *lp, CCtsp_lpcut_in *c),
    pushsparse (int i, int n, double d, double e, double **blist,
        double **llist, int **ilist);



int CCtsp_bb_init_lp (CCtsp_lp **lp, char *probname, int probnum, int ncount,
        CCdatagroup *dat, int *ptour, double initial_ub, CCtsp_lpcuts *pool,
        int silent, CCrandstate *rstate)
{
    return CCtsp_init_lp (lp, probname, probnum, (char *) NULL,
                  ncount, dat,
                  0, (int *) NULL, (int *) NULL,       /* elist  */
                  0, (int *) NULL, (int *) NULL, 0,    /* exlist */
                  ptour, initial_ub, pool, (CCtsp_lpcuts *) NULL, silent,
                  rstate);
}

int CCtsp_init_lp (CCtsp_lp **lp, char *probloc, int probnum,
        char *probfilename, int ncount, CCdatagroup *dat, int ecount,
        int *elist, int *elen, int excount, int *exlist, int *exlen,
        int exvalid, int *ptour, double initial_ub, CCtsp_lpcuts *pool,
        CCtsp_lpcuts *dominopool, int silent, CCrandstate *rstate)
{
    int rval = 0;
    double st, val;

    if (!ptour) {
        fprintf (stderr, "must have a permutation tour in CCtsp_init_lp\n");
        rval = 1; goto CLEANUP;
    }

    *lp = CC_SAFE_MALLOC (1, CCtsp_lp);
    if ( !(*lp) ) {
       rval = 1; goto CLEANUP;
    }

    CCtsp_init_tsp_lp_struct (*lp);
    (*lp)->perm = ptour;
    (*lp)->dat = dat;
    (*lp)->pool = pool;
    (*lp)->dominopool = dominopool;

    if (probfilename) {
        rval = CCtsp_read_probfile (*lp, probfilename, probloc, &ncount,
                                    silent);
        CCtsp_free_bigdual (&((*lp)->exact_dual));
    } else if (probnum != -1) {
        rval = CCtsp_read_probfile_id (*lp, probloc, probnum, &ncount, silent);
        CCtsp_free_bigdual (&((*lp)->exact_dual));
    } else {
        rval = first_lp (*lp, probloc, ncount, ecount, elist, elen, silent,
                         rstate);
    }
    if (rval) {
        fprintf (stderr, "CCtsp_read_probfile or first_lp failed\n");
        CCtsp_free_tsp_lp_struct (lp);
        rval = 1; goto CLEANUP;
    }

    if ((*lp)->fullcount == 0) {
        if (excount) {
            rval = CCtsp_edgelist_to_genadj (ncount, excount, exlist, exlen,
                           &((*lp)->fulladj), &((*lp)->fulladjspace));
            if (rval) {
                fprintf (stderr, "CCtsp_edgelist_to_genadj failed\n");
                CCtsp_free_tsp_lp_struct (lp);
                rval = 1; goto CLEANUP;
            }
            (*lp)->fullcount = excount;
            if (exvalid)
                (*lp)->full_edges_valid = 1;
        }
    }

    if (initial_ub < (*lp)->upperbound) {
        if (!silent) {
            printf ("Setting upperbound to the initial bound: %.2f\n",
                     initial_ub);
            fflush (stdout);
        }
        (*lp)->upperbound = initial_ub;
    }

    rval = CClp_init (&((*lp)->lp));
    if (rval) {
        fprintf (stderr, "CClp_init failed\n");
        CCtsp_free_tsp_lp_struct (lp);
        rval = 1; goto CLEANUP;
    }

#ifdef PERTURB_LPS
    if ( CClp_force_perturb ((*lp)->lp) ) {
        fprintf (stdout, "CClp_force_perturb failed, continuing anyway\n");
    }
#endif

    rval = load_lp (*lp, silent);
    if (rval) {
        fprintf (stderr, "load_lp failed\n");
        CCtsp_free_tsp_lp_struct (lp);
        rval = 1; goto CLEANUP;
    }

    CClp_free_warmstart (&((*lp)->warmstart));

    (*lp)->edge_life = CCtsp_EDGE_LIFE;
    (*lp)->cut_life  = CCtsp_CUT_LIFE;

    rval = CCtsp_add_branchhistory_to_lp (*lp);
    if (rval) {
        fprintf (stderr, "CCtsp_add_branchhistory_to_lp failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCpq_cuttree_trivial (&(*lp)->tightcuts, ncount, 0);
    if (rval) {
        fprintf (stderr, "CCpq_cuttree_trivial failed\n");
        rval = 1; goto CLEANUP;
    }

    CCutil_start_timer (&(*lp)->stats.misc_opt);
    rval = CClp_opt ((*lp)->lp, CClp_METHOD_DUAL);

    st = CCutil_stop_timer (&(*lp)->stats.misc_opt, 0);
    if (!silent) {
        printf ("Dual opt returned after %.2f seconds\n", st); fflush (stdout);
    }

    if (rval == 2) {
        fprintf (stderr, "Initial lp infeasible\n");
        goto CLEANUP;
    } else if (rval) {
        fprintf (stderr, "CClp_opt failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_update_result (*lp);
    if (rval) {
        fprintf (stderr, "CCtsp_update_result failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_get_lp_result (*lp, &val, (double *) NULL, (int *) NULL,
                        (int **) NULL, (double **) NULL, (double **) NULL,
                        (double **) NULL, (double **) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_get_lp_result failed\n");
        rval = 1; goto CLEANUP;
    }

    if (!silent) {
        printf ("Initial LP value: %.6f\n", val);
        fflush (stdout);
    }
    
CLEANUP:

    return rval;
}

static int first_lp (CCtsp_lp *lp, char *probloc, int ncount,
        int ecount, int *elist, int *elen, int silent, CCrandstate *rstate)
{
    double v;
    int *basis = (int *) NULL;
    int *match = (int *) NULL;
    int rval = 0;
    int i;
    int e;
    double st;
    CClp_info *info = (CClp_info *) NULL;

    CCutil_start_timer (&lp->stats.misc);

    rval = CCtsp_init_cliquehash (&lp->cuts, 2*ncount);
    CCcheck_rval (rval, "CCtsp_init_cliquehash failed");
    rval = CCtsp_init_dominohash (&lp->cuts, 2*ncount);
    CCcheck_rval (rval, "CCtsp_init_dominohash failed");
    rval = CCtsp_build_lpgraph (&lp->graph, ncount, ecount, elist, elen);
    CCcheck_rval (rval, "CCtsp_build_lpgraph failed");
    rval = CCtsp_build_lpadj (&lp->graph, 0, ecount);
    CCcheck_rval (rval, "CCtsp_build_lpadj failed");

    lp->probloc = CCutil_strdup (probloc);
    lp->problabel = CCtsp_problabel (probloc);
    if (lp->probloc == (char *) NULL ||
        lp->problabel == (char *) NULL) {
        rval = 1;
        goto CLEANUP;
    }

    lp->id = 0;
    lp->parent_id = -1;
    lp->full_edges_valid = 0;

    basis = CC_SAFE_MALLOC (2 * ncount, int);
    if (!basis) {
        rval = 1;
        goto CLEANUP;
    }
    match = CC_SAFE_MALLOC ((6 * ncount) + 1, int);
    if (!match) {
        rval = 1;
        goto CLEANUP;
    }

    rval = CCfmatch_fractional_2match (ncount, ecount, elist, elen,
            (CCdatagroup *) NULL, &v, match, (int *) NULL, basis, 1, 
            silent, rstate);
             
    if (rval) {
        fprintf (stderr, "Fractional matching routine failed\n");
        CC_FREE (basis, int);
        CC_FREE (match, int);
        return 0;
    }

    rval = CClp_create_info (&info, ncount, ecount);
    if (rval) {
        fprintf (stderr, "CClp_create_info failed\n");
        goto CLEANUP;
    }

    for (i=0; i<ncount; i++) {
        CClp_set_row_active (info, i);
    }
    for (i=0; i<ecount; i++) {
        CClp_set_col_inactive (info, i);
    }

    for (i = 0; i < ncount; i++) {
        e = CCtsp_find_edge (&lp->graph, basis[2*i], basis[2*i+1]);
        if (e < 0) {
            fprintf (stderr, "Basis contains edge %d,%d not in edgelist\n",
                     basis[2*i],basis[2*i+1]);
        } else {
            CClp_set_col_active (info, e);
        }
    }

    for (i=0; i <= 6*ncount && match[i] > -1; i += 3) {
        e = CCtsp_find_edge (&lp->graph, match[i], match[i+1]);
        if (e < 0) {
            fprintf (stderr, "Matching contains edge %d,%d not in edgelist\n",
                     match[i], match[i+1]);
        } else {
            if (!CClp_is_col_active (info, e)) {
                if (match[i+2] == 1) {
                    fprintf (stderr, "Edge at 0.5 not in basis\n");
                }
                CClp_set_col_upper (info, e);
            }
        }
    }

    CC_FREE (basis, int);
    CC_FREE (match, int);

    rval = CClp_build_warmstart (&(lp->warmstart), info);
    if (rval) {
        fprintf (stderr, "CClp_build_warmstart failed\n");
        goto CLEANUP;
    }
    
    st = CCutil_stop_timer (&lp->stats.misc, 0);
    if (!silent) {
        printf ("Total Time for first_lp: %.2f (seconds)\n", st);
        fflush (stdout);
    }

    rval = 0;

CLEANUP:

    CC_IFFREE (basis, int);
    CC_IFFREE (match, int);
    CClp_free_info (&info);

    return rval;
}

int CCtsp_build_lpgraph (CCtsp_lpgraph *g, int ncount, int ecount,
        int *elist, int *elen)
{
    int i;
    CCtsp_lpnode *n;
    CCtsp_lpedge *e;

    g->ncount = ncount;
    g->ecount = ecount;
    g->nodes = CC_SAFE_MALLOC (ncount, CCtsp_lpnode);
    if (!g->nodes) {
        return 1;
    }
    g->edges = CC_SAFE_MALLOC (ecount, CCtsp_lpedge);
    if (!g->edges) {
        CC_FREE (g->nodes, CCtsp_lpnode);
        return 1;
    }
    g->espace = ecount;
    n = g->nodes;
    e = g->edges;

    for (i = 0; i < ncount; i++) {
        n[i].mark = 0;
    }
    for (i=0; i<ecount; i++) {
        if (elist[2*i] < elist[2*i+1]) {
            e[i].ends[0] = elist[2*i];
            e[i].ends[1] = elist[2*i+1];
        } else {
            e[i].ends[0] = elist[2*i+1];
            e[i].ends[1] = elist[2*i];
        }
        e[i].fixed = 0;
        e[i].branch = 0;
        e[i].age = 0;
        if (elen) {
            e[i].len = elen[i];
        } else {
            e[i].len = 0;
        }
        e[i].coefnext = -2;
        e[i].coef = 0;
    }
    return 0;
}

int CCtsp_build_lpadj (CCtsp_lpgraph *g, int estart, int eend)
{
    CCtsp_lpadj *a;
    CCtsp_lpnode *n = g->nodes;
    CCtsp_lpedge *e = g->edges;
    int i, j;

    if (g->adjspace) {
        if (g->adjstart == estart && g->adjend == eend) {
            return 0;
        } else {
            CC_FREE (g->adjspace, CCtsp_lpadj);
        }
    }

    if (estart >= eend) {
        g->adjstart = estart;
        g->adjend = eend;
        for (i=0; i<g->ncount; i++) {
            n[i].deg = 0;
            n[i].adj = (CCtsp_lpadj *) NULL;
        }
        return 0;
    }

    g->adjspace = CC_SAFE_MALLOC ((eend - estart)*2, CCtsp_lpadj);
    if (!g->adjspace) {
        return 1;
    }
    a = g->adjspace;
    for (i=0; i<g->ncount; i++) {
        n[i].deg = 0;
    }
    for (i=estart; i<eend; i++) {
        n[e[i].ends[0]].deg++;
        n[e[i].ends[1]].deg++;
    }
    for (i=0; i<g->ncount; i++) {
        n[i].adj = a;
        a += n[i].deg;
        n[i].deg = 0;
    }
    for (i=estart; i<eend; i++) {
        j = e[i].ends[0];
        a = &n[j].adj[n[j].deg];
        a->to = e[i].ends[1];
        a->edge = i;
        n[j].deg++;
        j = e[i].ends[1];
        a = &n[j].adj[n[j].deg];
        a->to = e[i].ends[0];
        a->edge = i;
        n[j].deg++;
    }
    g->adjstart = estart;
    g->adjend = eend;

    return 0;
}

int CCtsp_find_edge (CCtsp_lpgraph *g, int from, int to)
{
    int t;
    int i;
    CCtsp_lpnode *f;

    if (from > to) {
        CC_SWAP (from, to, t);
    }

    f = &g->nodes[from];
    for (i=0; i<f->deg; i++) {
        if (f->adj[i].to == to) {
            return f->adj[i].edge;
        }
    }
    return -1;
}

static int find_edge_full (CCtsp_lp *lp, int from, int to)
{
    int i;
    CCtsp_genadjobj *a;

    /* returns 1 if it is there */

    if (!lp->fulladj)
        return 0;

    if (from > to) {
        CC_SWAP (from, to, i);
    }

    a = lp->fulladj[from].list;
    for (i = lp->fulladj[from].deg-1; i >= 0; i--) {
        if (a[i].end == to) {
            return 1;
        }
    }
    return 0;
}

int CCtsp_inspect_full_edges (CCtsp_lp *lp)
{
    int i;
    int ecount = lp->graph.ecount;
    CCtsp_lpedge *edges = lp->graph.edges;

    for (i = 0; i < ecount; i++) {
        if (find_edge_full (lp, edges[i].ends[0], edges[i].ends[1]) == 0) {
            printf ("edge (%d,%d) not in full list\n",
                                 edges[i].ends[0], edges[i].ends[1]);
            fflush (stdout);
            return 1;
        }
    }
    return 0;
}

void CCtsp_init_tsp_lpcuts_struct (CCtsp_lpcuts *c)
{
    if (c) {
        c->cutcount        = 0;
        c->savecount       = 0;
        c->cliqueend       = 0;
        c->cutspace        = 0;
        c->cliquespace     = 0;
        c->cliquehashsize  = 0;
        c->cliquefree      = 0;
        c->cliquehash      = (int *) NULL;
        c->cuts            = (CCtsp_lpcut *) NULL;
        c->cliques         = (CCtsp_lpclique *) NULL;
        c->cuthash         = (CCgenhash *) NULL;
        c->tempcuthash     = (char *) NULL;
        c->tempcuthashsize = 0;
        c->dominoend       = 0;
        c->dominospace     = 0;
        c->dominohashsize  = 0;
        c->dominofree      = 0;
        c->dominohash      = (int *) NULL;
        c->dominos         = (CCtsp_lpdomino *) NULL;
        c->workloads       = (double *) NULL;
    }
}

void CCtsp_init_tsp_lp_struct (CCtsp_lp *lp)
{
    CCtsp_init_lpgraph_struct (&(lp->graph));
    CCtsp_init_tsp_lpcuts_struct (&(lp->cuts));

    lp->sparsifier = (CCtsp_qsparsegroup *) NULL;
    lp->perm       = (int *) NULL;
    lp->dat        = (CCdatagroup *) NULL;
    lp->lp         = (CClp *) NULL;


    lp->pool             = (CCtsp_lpcuts *) NULL;
    lp->dominopool       = (CCtsp_lpcuts *) NULL;
    lp->fullcount        = 0;
    lp->fulladj          = (CCtsp_genadj *) NULL;
    lp->fulladjspace     = (CCtsp_genadjobj *) NULL;
    lp->fixededges       = (int *) NULL;
    lp->nfixededges      = 0;
    lp->problabel        = (char *) NULL;
    lp->probloc          = (char *) NULL;
    lp->id               = -1;
    lp->parent_id        = -1;
    lp->root             = 0;
    lp->upperbound       =  CCtsp_LP_MAXDOUBLE;
    lp->lowerbound       = -CCtsp_LP_MAXDOUBLE;
    lp->exact_lowerbound = CCbigguy_MINBIGGUY;
    lp->exact_dual       = (CCtsp_bigdual *) NULL;
    lp->infeasible       = 0;
    lp->full_edges_valid = 0;
    lp->warmstart        = (CClp_warmstart *) NULL;

    CCtsp_init_lpcut_in (&lp->cutqueue);
    lp->cutqueue.next = &lp->cutqueue;
    lp->cutqueue.prev = &lp->cutqueue;

    lp->result.ub     = 0.0;
    lp->result.lb     = 0.0;
    lp->result.ecount = 0;
    lp->result.elist  = (int *) NULL;
    lp->result.x      = (double *) NULL;
    lp->result.rc     = (double *) NULL;

    lp->branchhistory = (CCtsp_branchobj *) NULL;
    lp->branchdepth   = 0;

    CCtsp_init_statistics (&lp->stats);

    CCpq_cuttree_init (&lp->tightcuts);
}

void CCtsp_free_tsp_lp_struct (CCtsp_lp **lp)
{
    int i;

    if (!(*lp))  return;

    CCtsp_free_lpgraph (&((*lp)->graph));

    (*lp)->perm = (int *) NULL;      /* perm is owned by calling program */
    if ((*lp)->sparsifier) {
        CCtsp_free_qsparsify (&((*lp)->sparsifier));
    }
    (*lp)->dat = (CCdatagroup *) NULL; /* dat is owned by calling program */
    if ((*lp)->fulladjspace) {
        CC_FREE ((*lp)->fulladjspace, CCtsp_genadjobj);
        CC_IFFREE ((*lp)->fulladj, CCtsp_genadj);
    }
    (*lp)->fullcount = 0;
    CC_IFFREE ((*lp)->fixededges, int);
    (*lp)->nfixededges = 0;

    CClp_free (&((*lp)->lp));

    if ((*lp)->cuts.cuts) {
        for (i=0; i<(*lp)->cuts.cutcount; i++) {
            CC_IFFREE ((*lp)->cuts.cuts[i].cliques, int);
            CC_IFFREE ((*lp)->cuts.cuts[i].dominos, int);
            CC_IFFREE ((*lp)->cuts.cuts[i].mods, CCtsp_sparser);
            CCtsp_free_skeleton (&(*lp)->cuts.cuts[i].skel);
        }
        CC_FREE ((*lp)->cuts.cuts, CCtsp_lpcut);
    }
    if ((*lp)->cuts.cliques) {
        for (i=0; i<(*lp)->cuts.cliqueend; i++) {
            CC_IFFREE ((*lp)->cuts.cliques[i].nodes, CCtsp_segment);
        }
        CC_FREE ((*lp)->cuts.cliques, CCtsp_lpclique);
    }
    if ((*lp)->cuts.dominos) {
        for (i=0; i<(*lp)->cuts.dominoend; i++) {
            CCtsp_free_lpdomino (&((*lp)->cuts.dominos[i]));
        }
        CC_FREE ((*lp)->cuts.dominos, CCtsp_lpdomino);
    }

    (*lp)->pool = (CCtsp_lpcuts *) NULL;  /* owned by calling routine */

    if ((*lp)->exact_dual) {
        CC_IFFREE ((*lp)->exact_dual->node_pi, CCbigguy);
        CC_IFFREE ((*lp)->exact_dual->cut_pi, CCbigguy);
        CC_FREE ((*lp)->exact_dual, CCtsp_bigdual);
    }
    CC_IFFREE ((*lp)->cuts.cliquehash, int);
    CC_IFFREE ((*lp)->cuts.dominohash, int);
    CC_IFFREE ((*lp)->problabel, char);
    CC_IFFREE ((*lp)->probloc, char);
    CClp_free_warmstart (&((*lp)->warmstart));

    CC_IFFREE ((*lp)->result.elist, int);
    CC_IFFREE ((*lp)->result.x, double);
    CC_IFFREE ((*lp)->result.rc, double);

    CCpq_cuttree_freetree (&(*lp)->tightcuts);

    if ((*lp)->branchhistory) {
        for (i = 0; i < (*lp)->branchdepth; i++) {
            CCtsp_free_branchobj (&((*lp)->branchhistory[i]));
        }
        CC_FREE ((*lp)->branchhistory, CCtsp_branchobj);
        (*lp)->branchdepth = 0;
    }

    CC_IFFREE (*lp, CCtsp_lp);
}

void CCtsp_init_lpgraph_struct (CCtsp_lpgraph *g)
{
    g->ncount = 0;
    g->ecount = 0;
    g->nodes = (CCtsp_lpnode *) NULL;
    g->edges = (CCtsp_lpedge *) NULL;
    g->adjspace = (CCtsp_lpadj *) NULL;
    g->adjstart = 0;
    g->adjend = 0;
    g->nodemarker = 0;
    g->espace = 0;
}

void CCtsp_free_lpgraph (CCtsp_lpgraph *g)
{
    CC_IFFREE (g->nodes, CCtsp_lpnode);
    CC_IFFREE (g->edges, CCtsp_lpedge);
    CC_IFFREE (g->adjspace, CCtsp_lpadj);
    g->espace = 0;
}

void CCtsp_init_statistics (CCtsp_statistics *stats)
{
    CCutil_init_timer (&stats->cutting_loop,
                       "Cutting Loop");
    CCutil_init_timer (&stats->cutting_inner_loop,
                       "Cutting inside loop");
    CCutil_init_timer (&stats->cuts_filecut,
                       "File cuts");
    CCutil_init_timer (&stats->cuts_filecut_opt,
                       "File cuts optimize");
    CCutil_init_timer (&stats->cuts_cutpool,
                       "Cutpool cuts");
    CCutil_init_timer (&stats->cuts_cutpool_opt,
                       "Cutpool cuts optimize");
    CCutil_init_timer (&stats->cuts_connect,
                       "Connect cuts");
    CCutil_init_timer (&stats->cuts_connect_opt,
                       "Connect cuts optimize");
    CCutil_init_timer (&stats->cuts_segment,
                       "Segment cuts");
    CCutil_init_timer (&stats->cuts_segment_opt,
                       "Segment cuts optimize");
    CCutil_init_timer (&stats->cuts_remotepool,
                       "Remote pool cuts");
    CCutil_init_timer (&stats->cuts_remotepool_opt,
                       "Remote pool cuts optimize");
    CCutil_init_timer (&stats->cuts_blockcomb,
                       "Block comb cuts");
    CCutil_init_timer (&stats->cuts_blockcomb_opt,
                       "Block comb cuts optimize");
    CCutil_init_timer (&stats->cuts_growcomb,
                       "Grow comb cuts");
    CCutil_init_timer (&stats->cuts_growcomb_opt,
                       "Grow comb cuts optimize");
    CCutil_init_timer (&stats->cuts_exactsubtour,
                       "Exact subtour cuts");
    CCutil_init_timer (&stats->cuts_exactsubtour_opt,
                       "Exact subtour cuts optimize");
    CCutil_init_timer (&stats->cuts_tighten_lp,
                       "Tighten lp cuts");
    CCutil_init_timer (&stats->cuts_tighten_lp_opt,
                       "Tighten lp cuts optimize");
    CCutil_init_timer (&stats->cuts_tighten_lp_close,
                       "Tighten lp (close) cuts");
    CCutil_init_timer (&stats->cuts_tighten_lp_close_opt,
                       "Tighten lp (close) cuts optimize");
    CCutil_init_timer (&stats->cuts_decker_lp,
                       "Ddecker lp cuts");
    CCutil_init_timer (&stats->cuts_decker_lp_opt,
                       "Ddecker lp cuts optimize");
    CCutil_init_timer (&stats->cuts_decker_lp_close,
                       "Ddecker lp (close) cuts");
    CCutil_init_timer (&stats->cuts_decker_lp_close_opt,
                       "Ddecker lp (close) cuts optimize");
    CCutil_init_timer (&stats->cuts_star_lp,
                       "Star lp cuts");
    CCutil_init_timer (&stats->cuts_star_lp_opt,
                       "Star lp cuts optimize");
    CCutil_init_timer (&stats->cuts_handling_lp,
                       "Handling lp cuts");
    CCutil_init_timer (&stats->cuts_handling_lp_opt,
                       "Handling lp cuts optimize");
    CCutil_init_timer (&stats->cuts_cliquetree_lp,
                       "Cliquetree lp cuts");
    CCutil_init_timer (&stats->cuts_cliquetree_lp_opt,
                       "Cliquetree lp cuts optimize");
    CCutil_init_timer (&stats->cuts_teething_lp,
                       "Teething lp cuts");
    CCutil_init_timer (&stats->cuts_teething_lp_opt,
                       "Teething lp cuts optimize");

    CCutil_init_timer (&stats->cuts_fastblossom,
                       "Fast blossom cuts");
    CCutil_init_timer (&stats->cuts_fastblossom_opt,
                       "Fast blossom cuts optimize");

    CCutil_init_timer (&stats->cuts_ghfastblossom,
                       "Fast GH-blossom cuts");
    CCutil_init_timer (&stats->cuts_ghfastblossom_opt,
                       "Fast GH-blossom cuts optimize");

    CCutil_init_timer (&stats->cuts_exactblossom,
                       "Exact blossom cuts");
    CCutil_init_timer (&stats->cuts_exactblossom_opt,
                       "Exact blossom cuts optimize");

    CCutil_init_timer (&stats->cuts_tighten_pool,
                       "Tighten pool cuts");
    CCutil_init_timer (&stats->cuts_tighten_pool_opt,
                       "Tighten pool cuts optimize");
    CCutil_init_timer (&stats->cuts_decker_pool,
                       "Ddecker pool cuts");
    CCutil_init_timer (&stats->cuts_decker_pool_opt,
                       "Ddecker pool cuts optimize");

    CCutil_init_timer (&stats->cuts_star_pool,
                       "Star pool cuts");
    CCutil_init_timer (&stats->cuts_star_pool_opt,
                       "Star pool cuts optimize");
    CCutil_init_timer (&stats->cuts_handling_pool,
                       "Handling lp cuts");
    CCutil_init_timer (&stats->cuts_handling_pool_opt,
                       "Handling lp cuts optimize");

    CCutil_init_timer (&stats->cuts_teething_pool,
                       "Teething pool cuts");
    CCutil_init_timer (&stats->cuts_teething_pool_opt,
                       "Teething pool cuts optimize");
    CCutil_init_timer (&stats->cuts_consecutiveones,
                       "Consecutiveones cuts");
    CCutil_init_timer (&stats->cuts_consecutiveones_opt,
                       "Consecutiveones cuts optimize");
    CCutil_init_timer (&stats->cuts_necklace,
                       "Necklace cuts");
    CCutil_init_timer (&stats->cuts_necklace_opt,
                       "Necklace cuts optimize");
    CCutil_init_timer (&stats->cuts_localcut,
                       "Localcut cuts");
    CCutil_init_timer (&stats->cuts_localcut_opt,
                       "Localcut cuts optimize");
    CCutil_init_timer (&stats->cuts_extraconnect,
                       "Extra connect cuts");
    CCutil_init_timer (&stats->cuts_extraconnect_opt,
                       "Extra connect cuts optimize");
    CCutil_init_timer (&stats->sparse_edge_check,
                       "Sparse edge check");
    CCutil_init_timer (&stats->full_edge_check,
                       "Full edge check");
    CCutil_init_timer (&stats->addcuts,
                       "Add cuts");
    CCutil_init_timer (&stats->addcuts_opt,
                       "Add cuts optimize");
    CCutil_init_timer (&stats->agecuts,
                       "Age cuts");
    CCutil_init_timer (&stats->agecuts_opt,
                       "Age cuts optimize");
    CCutil_init_timer (&stats->ageedges,
                       "Age edges");
    CCutil_init_timer (&stats->ageedges_opt,
                       "Age edges optimize");
    CCutil_init_timer (&stats->addbad,
                       "Add edges");
    CCutil_init_timer (&stats->addbad_opt,
                       "Add edges optimize");
    CCutil_init_timer (&stats->strongbranch,
                       "Strong branching");
    CCutil_init_timer (&stats->strongbranch_opt,
                       "Strong branching optimize");
    CCutil_init_timer (&stats->linkern,
                       "Lin-Kernighan");
    CCutil_init_timer (&stats->misc,
                       "Miscellaneous");
    CCutil_init_timer (&stats->misc_opt,
                       "Miscellaneous optimize");
    CCutil_init_timer (&stats->total,
                       "Total");
    stats->problem_cnt = 0;

    CCtsp_init_tighten_info (&stats->tighten_stats);
    CCtsp_init_tighten_info (&stats->extra_tighten_stats);

    fflush (stdout);
}

void CCtsp_output_statistics (CCtsp_statistics *stats)
{
    double y = 0.0;
    double z = 0.0;
    
    printf ("    Cutting times:\n");
    y += CCutil_total_timer (&stats->cuts_filecut, 4);
    z += CCutil_total_timer (&stats->cuts_filecut_opt, 0);
    y += CCutil_total_timer (&stats->cuts_cutpool, 4);
    z += CCutil_total_timer (&stats->cuts_cutpool_opt, 0);
    y += CCutil_total_timer (&stats->cuts_connect, 4);
    z += CCutil_total_timer (&stats->cuts_connect_opt, 0);
    y += CCutil_total_timer (&stats->cuts_segment, 4);
    z += CCutil_total_timer (&stats->cuts_segment_opt, 0);
    y += CCutil_total_timer (&stats->cuts_remotepool, 4);
    z += CCutil_total_timer (&stats->cuts_remotepool_opt, 0);
    y += CCutil_total_timer (&stats->cuts_blockcomb, 4);
    z += CCutil_total_timer (&stats->cuts_blockcomb_opt, 0);
    y += CCutil_total_timer (&stats->cuts_growcomb, 4);
    z += CCutil_total_timer (&stats->cuts_growcomb_opt, 0);
    y += CCutil_total_timer (&stats->cuts_exactsubtour, 4);
    z += CCutil_total_timer (&stats->cuts_exactsubtour_opt, 0);
    y += CCutil_total_timer (&stats->cuts_tighten_lp, 4);
    z += CCutil_total_timer (&stats->cuts_tighten_lp_opt, 0);
    y += CCutil_total_timer (&stats->cuts_tighten_lp_close, 4);
    z += CCutil_total_timer (&stats->cuts_tighten_lp_close_opt, 0);
    y += CCutil_total_timer (&stats->cuts_decker_lp, 4);
    z += CCutil_total_timer (&stats->cuts_decker_lp_opt, 0);
    y += CCutil_total_timer (&stats->cuts_decker_lp_close, 4);
    z += CCutil_total_timer (&stats->cuts_decker_lp_close_opt, 0);

    y += CCutil_total_timer (&stats->cuts_star_lp, 4);
    z += CCutil_total_timer (&stats->cuts_star_lp_opt, 0);
    y += CCutil_total_timer (&stats->cuts_handling_lp, 4);
    z += CCutil_total_timer (&stats->cuts_handling_lp_opt, 0);

    y += CCutil_total_timer (&stats->cuts_cliquetree_lp, 4);
    z += CCutil_total_timer (&stats->cuts_cliquetree_lp_opt, 0);
    y += CCutil_total_timer (&stats->cuts_teething_lp, 4);
    z += CCutil_total_timer (&stats->cuts_teething_lp_opt, 0);

    y += CCutil_total_timer (&stats->cuts_ghfastblossom, 4);
    z += CCutil_total_timer (&stats->cuts_ghfastblossom_opt, 0);

    y += CCutil_total_timer (&stats->cuts_fastblossom, 4);
    z += CCutil_total_timer (&stats->cuts_fastblossom_opt, 0);

    y += CCutil_total_timer (&stats->cuts_exactblossom, 4);
    z += CCutil_total_timer (&stats->cuts_exactblossom_opt, 0);

    y += CCutil_total_timer (&stats->cuts_tighten_pool, 4);
    z += CCutil_total_timer (&stats->cuts_tighten_pool_opt, 0);
    y += CCutil_total_timer (&stats->cuts_decker_pool, 4);
    z += CCutil_total_timer (&stats->cuts_decker_pool_opt, 0);

    y += CCutil_total_timer (&stats->cuts_star_pool, 4);
    z += CCutil_total_timer (&stats->cuts_star_pool_opt, 0);
    y += CCutil_total_timer (&stats->cuts_handling_pool, 4);
    z += CCutil_total_timer (&stats->cuts_handling_pool_opt, 0);

    y += CCutil_total_timer (&stats->cuts_teething_pool, 4);
    z += CCutil_total_timer (&stats->cuts_teething_pool_opt, 0);
    y += CCutil_total_timer (&stats->cuts_consecutiveones, 4);
    z += CCutil_total_timer (&stats->cuts_consecutiveones_opt, 0);
    y += CCutil_total_timer (&stats->cuts_necklace, 4);
    z += CCutil_total_timer (&stats->cuts_necklace_opt, 0);
    y += CCutil_total_timer (&stats->cuts_localcut, 4);
    z += CCutil_total_timer (&stats->cuts_localcut_opt, 0);
    y += CCutil_total_timer (&stats->cuts_extraconnect, 4);
    z += CCutil_total_timer (&stats->cuts_extraconnect_opt, 0);
    
    printf ("    Cutting totals:\n");
    printf ("    Cutting plane heuristics     %9.2f seconds\n", y);
    printf ("    Cutting plane optimization   %9.2f seconds\n", z);
    CCutil_total_timer (&stats->cutting_inner_loop, 4);
    CCutil_total_timer (&stats->cutting_loop, 4);
    CCutil_total_timer (&stats->sparse_edge_check, 4);
    CCutil_total_timer (&stats->full_edge_check, 4);

    printf ("    LP processing, included in optimization and edge check totals\n");
    CCutil_total_timer (&stats->addcuts, 4);
    CCutil_total_timer (&stats->addcuts_opt, 4);
    CCutil_total_timer (&stats->agecuts, 4);
    CCutil_total_timer (&stats->agecuts_opt, 4);
    CCutil_total_timer (&stats->addbad, 4);
    CCutil_total_timer (&stats->addbad_opt, 4);
    CCutil_total_timer (&stats->ageedges, 4);
    CCutil_total_timer (&stats->ageedges_opt, 4);

    printf ("    Other times\n");
    CCutil_total_timer (&stats->strongbranch, 4);
    CCutil_total_timer (&stats->strongbranch_opt, 4);
    CCutil_total_timer (&stats->linkern, 4);
    CCutil_total_timer (&stats->misc, 4);
    CCutil_total_timer (&stats->misc_opt, 4);
    CCutil_total_timer (&stats->total, 4);

    printf ("General cut tightening\n");
    CCtsp_print_tighten_info (&stats->tighten_stats);
    printf ("Heuristic lp and pool tightening\n");
    CCtsp_print_tighten_info (&stats->extra_tighten_stats);
}

static int load_lp (CCtsp_lp *lp, int silent)
{
    int rval = 0;
    int nzcnt;
    int i, j;
    int xrhs;
    double st;

    double *obj    = (double *) NULL;
    double *rhs    = (double *) NULL;
    char   *sense  = (char *)   NULL;
    int    *matbeg = (int *)    NULL;
    int    *matcnt = (int *)    NULL;
    int    *matind = (int *)    NULL;
    double *matval = (double *) NULL;
    double *lb     = (double *) NULL;
    double *ub     = (double *) NULL;

    CCutil_start_timer (&lp->stats.misc);
    
    if (!silent) {
        printf ("Loading lp..."); fflush (stdout);
    }

    rhs      = CC_SAFE_MALLOC (lp->graph.ncount + lp->cuts.cutcount, double);
    sense    = CC_SAFE_MALLOC (lp->graph.ncount + lp->cuts.cutcount, char);
    if (!rhs || !sense ) {
       fprintf (stderr, "not enough memory to load problem\n");
       rval = 1;
       goto CLEANUP;
    }

    rval = build_lp_cols (&lp->graph, &lp->cuts, 0, lp->graph.ecount,
                          &nzcnt, &obj, &matbeg, &matcnt, &matind, &matval,
                          &lb, &ub);
    if (rval)  goto CLEANUP;

    for (i = 0; i < lp->graph.ncount; i++) {
        rhs[i] = 2.0;
        sense[i] = 'E';
    }
    for (i = 0; i < lp->cuts.cutcount; i++) {
        xrhs = lp->cuts.cuts[i].rhs;
        for (j = 0; j < lp->cuts.cuts[i].modcount; j++) {
            xrhs += 2*(((int) lp->cuts.cuts[i].mods[j].mult)-128);
        }
        rhs[lp->graph.ncount + i] = xrhs;
        sense[lp->graph.ncount + i] = lp->cuts.cuts[i].sense;
    }

    rval = CClp_loadlp (lp->lp, lp->problabel, lp->graph.ecount,
            lp->graph.ncount + lp->cuts.cutcount, 1, obj, rhs, sense,
            matbeg, matcnt, matind, matval, lb, ub);

    if (rval) {
        fprintf (stderr, "couldn't load problem\n");
        goto CLEANUP;
    }

    if (lp->warmstart) {
        rval = CClp_load_warmstart (lp->lp, lp->warmstart);
        if (rval) {
            fprintf (stderr, "CClp_load_warmstart failed\n");
        }
    } else {
        fprintf (stderr, "No warmstart, stumbling on anyway\n");
    }

    st = CCutil_stop_timer (&lp->stats.misc, 0);
    if (!silent) {
        printf ("done in %.2f seconds\n", st);
        fflush (stdout);
    }

    if (!silent) {
        printf ("LP has:  %d rows  %d columns  %d nonzeros\n",
                CClp_nrows (lp->lp), CClp_ncols (lp->lp),
                CClp_nnonzeros (lp->lp)); fflush (stdout);
        fflush (stdout);
    }

CLEANUP:

    CC_IFFREE (obj, double);
    CC_IFFREE (rhs, double);
    CC_IFFREE (sense, char);
    CC_IFFREE (matbeg, int);
    CC_IFFREE (matcnt, int);
    CC_IFFREE (matind, int);
    CC_IFFREE (matval, double);
    CC_IFFREE (lb, double);
    CC_IFFREE (ub, double);

    return rval;
}

static int build_lp_cols (CCtsp_lpgraph *g, CCtsp_lpcuts *cuts, int estart,
        int eend, int *pnzcnt, double **pobj, int **pmatbeg, int **pmatcnt,
        int **pmatind, double **pmatval, double **plb, double **pub)
{
    int rval;
    int nzcnt;
    int ncols = eend - estart;
    double *obj = (double *) NULL;
    int *matbeg = (int *) NULL;
    int *matcnt = (int *) NULL;
    int *matind = (int *) NULL;
    double *matval = (double *) NULL;
    double *lb = (double *) NULL;
    double *ub = (double *) NULL;
    int i;
    int nzlist, nznext;

    if (estart >= eend) {
        fprintf (stderr, "No columns for build_lp_cols to build\n");
        return 1;
    }

    rval = CCtsp_build_lpadj (g, estart, eend);
    if (rval) goto CLEANUP;

    obj = CC_SAFE_MALLOC (ncols, double);
    lb = CC_SAFE_MALLOC (ncols, double);
    ub = CC_SAFE_MALLOC (ncols, double);
    matbeg = CC_SAFE_MALLOC (ncols, int);
    matcnt = CC_SAFE_MALLOC (ncols, int);
    if (!obj || !lb || !ub || !matbeg || !matcnt) {
        rval = 1;
        goto CLEANUP;
    }
    for (i = 0; i < ncols; i++) {
        obj[i] = g->edges[i+estart].len;
        if (g->edges[i+estart].fixed || g->edges[i+estart].branch > 0) {
            lb[i] = 1.0;
        } else {
            lb[i] = 0.0;
        }
        if (g->edges[i+estart].branch < 0) {
            ub[i] = 0.0;
        } else {
            ub[i] = 1.0;
        }
        matcnt[i] = 2;
    }

    for (i=0; i<cuts->cutcount; i++) {
        nzlist = CCtsp_lpcut_nzlist (g, &cuts->cuts[i], cuts->cliques,
                               cuts->dominos, 1);
        while (nzlist != -1) {
            nznext = g->edges[nzlist].coefnext;
            g->edges[nzlist].coefnext = -2;
            if (g->edges[nzlist].coef) {
                g->edges[nzlist].coef = 0;
                matcnt[nzlist - estart]++;
            }
            nzlist = nznext;
        }
    }

    nzcnt = 0;
    for (i=0; i<ncols; i++) {
        matbeg[i] = nzcnt;
        nzcnt += matcnt[i];
        matcnt[i] = 0;
    }
    matind = CC_SAFE_MALLOC (nzcnt, int);
    matval = CC_SAFE_MALLOC (nzcnt, double);
    if (!matind || !matval) {
        rval = 1;
        goto CLEANUP;
    }

    for (i=0; i<ncols; i++) {
        matval[matbeg[i] + matcnt[i]] = 1.0;
        matind[matbeg[i] + matcnt[i]] = g->edges[estart+i].ends[0];
        matcnt[i]++;
        matval[matbeg[i] + matcnt[i]] = 1.0;
        matind[matbeg[i] + matcnt[i]] = g->edges[estart+i].ends[1];
        matcnt[i]++;
    }


    for (i=0; i<cuts->cutcount; i++) {
        nzlist = CCtsp_lpcut_nzlist (g, &cuts->cuts[i], cuts->cliques,
                                    cuts->dominos, 1);
        while (nzlist != -1) {
            nznext = g->edges[nzlist].coefnext;
            g->edges[nzlist].coefnext = -2;
            if (g->edges[nzlist].coef) {
                matval[matbeg[nzlist-estart] + matcnt[nzlist-estart]] =
                        g->edges[nzlist].coef;
                matind[matbeg[nzlist-estart] + matcnt[nzlist-estart]] =
                        g->ncount + i;
                matcnt[nzlist-estart]++;
                g->edges[nzlist].coef = 0;
            }
            nzlist = nznext;
        }
    }

    if (pnzcnt) *pnzcnt = nzcnt;
    if (pobj) *pobj = obj;
    else CC_FREE (obj, double);
    if (pmatbeg) *pmatbeg = matbeg;
    else CC_FREE (matbeg, int);
    if (pmatcnt) *pmatcnt = matcnt;
    else CC_FREE (matcnt, int);
    if (pmatind) *pmatind = matind;
    else CC_FREE (matind, int);
    if (pmatval) *pmatval = matval;
    else CC_FREE (matval, double);
    if (plb) *plb = lb;
    else CC_FREE (lb, double);
    if (pub) *pub = ub;
    else CC_FREE (ub, double);

    return 0;

CLEANUP:

    CC_IFFREE (obj, double);
    CC_IFFREE (matbeg, int);
    CC_IFFREE (matcnt, int);
    CC_IFFREE (matind, int);
    CC_IFFREE (matval, double);
    CC_IFFREE (lb, double);
    CC_IFFREE (ub, double);

    return rval;
}

int CCtsp_resparsify_lp (CCtsp_lp *lp, int silent)
{
    CCtsp_lpcuts *lpcuts = &lp->cuts;
    int cutcount = lpcuts->cutcount;
    CCtsp_lpcut *cuts = lpcuts->cuts;
    CCtsp_lpgraph *g = &lp->graph;
    CCtsp_sparser *newmods = (CCtsp_sparser *) NULL;
    CClp_info *binfo = (CClp_info *) NULL;
    int newmodcount = 0;
    int saved;
    int nzlist;
    int rval;
    int i;

    CCutil_start_timer (&lp->stats.misc);

    if (!silent) {
        printf ("LP has:  %d rows  %d columns  %d nonzeros\n",
                CClp_nrows (lp->lp), CClp_ncols (lp->lp),
                CClp_nnonzeros (lp->lp));
        fflush (stdout);
    }
   

    for (i=0; i< cutcount; i++) {
        CC_IFFREE (cuts[i].mods, CCtsp_sparser);
        cuts[i].modcount = 0;
        nzlist = CCtsp_lpcut_nzlist (g, &cuts[i], lpcuts->cliques,
                                      lpcuts->dominos, 1);
        rval = CCtsp_qsparsify (&lp->sparsifier, g, &nzlist, &newmodcount,
                                &newmods, &saved);
        if (rval) {
            fprintf (stderr, "CCtsp_qsparsify failed\n");
            clear_nzlist (g, nzlist);
            goto CLEANUP;
        }
        cuts[i].mods = newmods;
        cuts[i].modcount = newmodcount;
        newmods = (CCtsp_sparser *) NULL;
        newmodcount = 0;
        clear_nzlist (g, nzlist);
    }

    CCutil_stop_timer (&lp->stats.misc, 0);
    
    /* Bico, 29.5.00 - since we alter the LP, the norms will no longer be */
    /*                 valid, so build just basis via a call to get_info. */

    rval = CClp_get_info (lp->lp, &binfo);
    if (rval) {
        fprintf (stderr, "CClp_get_info failed\n");
        goto CLEANUP;
    }

    rval = CClp_build_warmstart (&lp->warmstart, binfo);
    if (rval) {
        fprintf (stderr, "CClp_build_warmstart failed\n");
        goto CLEANUP;
    }
    CClp_free_info (&binfo);

    CClp_freelp (&lp->lp);

    rval = load_lp (lp, silent);
    if (rval) {
        fprintf (stderr, "load_lp failed\n");
        goto CLEANUP;
    }

    rval = 0;

    CCutil_start_timer (&lp->stats.misc_opt);
    rval = CClp_opt (lp->lp, CClp_METHOD_DUAL);
    if (rval) {
        fprintf (stderr, "reoptimization after sparsify failed\n");
        goto CLEANUP;
    }
    CCutil_stop_timer (&lp->stats.misc_opt, 0);
    
CLEANUP:
    CClp_free_info (&binfo);
    CClp_free_warmstart (&lp->warmstart);
    CC_IFFREE (newmods, CCtsp_sparser);
    return rval;
}

int CCtsp_update_result (CCtsp_lp *lp)
{
    CCtsp_lp_result new;
    int i;

    if (CClp_objval (lp->lp, &new.lb)) {
        return 1;
    }
    new.ub = lp->upperbound;
    new.elist = CC_SAFE_MALLOC (lp->graph.ecount*2, int);
    if (!new.elist) return 1;
    new.x = CC_SAFE_MALLOC (lp->graph.ecount, double);
    if (!new.x) {
        CC_FREE (new.elist, int);
        return 1;
    }
    new.rc = CC_SAFE_MALLOC (lp->graph.ecount, double);
    if (!new.rc) {
        CC_FREE (new.x, double);
        CC_FREE (new.elist, int);
        return 1;
    }

    if (CClp_x (lp->lp, new.x)) {
        CC_FREE (new.rc, double);
        CC_FREE (new.x, double);
        CC_FREE (new.elist, int);
        return 1;
    }

    if (CClp_rc (lp->lp, new.rc)) {
        CC_FREE (new.rc, double);
        CC_FREE (new.x, double);
        CC_FREE (new.elist, int);
        return 1;
    }

    new.ecount = lp->graph.ecount;
    for (i=0; i<new.ecount; i++) {
        new.elist[2*i] = lp->graph.edges[i].ends[0];
        new.elist[2*i+1] = lp->graph.edges[i].ends[1];
    }

    CC_IFFREE (lp->result.elist, int);
    CC_IFFREE (lp->result.x, double);
    CC_IFFREE (lp->result.rc, double);

    lp->result = new;

    return 0;
}

int CCtsp_get_lp_result (CCtsp_lp *lp, double *lb, double *ub, int *ecount,
        int **elist, double **x, double **rc, double **node_pi,
        double **cut_pi)
{
    int *myelist = (int *) NULL;
    double *myx = (double *) NULL;
    double *myrc = (double *) NULL;
    double *mynode_pi = (double *) NULL;
    double *mycut_pi = (double *) NULL;
    int i;
    int rval = 0;

    if ((elist || x || rc) && lp->result.ecount == 0){
       fprintf (stderr, "lp->result is not initialized\n");
       return 1;
    }

    if (elist) {
        myelist = CC_SAFE_MALLOC (2*lp->result.ecount, int);
        if (!myelist) {
            fprintf (stderr, "out of memory in CCtsp_get_lp_result\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    if (x) {
        myx = CC_SAFE_MALLOC (lp->result.ecount, double);
        if (!myx) {
            fprintf (stderr, "out of memory in CCtsp_get_lp_result\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    if (rc) {
        myrc = CC_SAFE_MALLOC (lp->result.ecount, double);
        if (!myrc) {
            fprintf (stderr, "out of memory in CCtsp_get_lp_result\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    if (node_pi) {
        mynode_pi = CC_SAFE_MALLOC (lp->graph.ncount, double);
        if (!mynode_pi) {
            fprintf (stderr, "out of memory in CCtsp_get_lp_result\n");
            rval = 1;
            goto CLEANUP;
        }
    }
    if (cut_pi && lp->cuts.cutcount) {
        mycut_pi = CC_SAFE_MALLOC (lp->cuts.cutcount, double);
        if (!mycut_pi) {
            fprintf (stderr, "out of memory in CCtsp_get_lp_result\n");
            rval = 1;
            goto CLEANUP;
        }
    }

    if (elist) {
        for (i=0; i<2*lp->result.ecount; i++) {
            myelist[i] = lp->result.elist[i];
        }
        *elist = myelist;
    }
    if (x) {
        for (i=0; i<lp->result.ecount; i++) {
            myx[i] = lp->result.x[i];
        }
        *x = myx;
    }
    if (rc) {
        for (i=0; i<lp->result.ecount; i++) {
            myrc[i] = lp->result.rc[i];
        }
        *rc = myrc;
    }
    if (node_pi || cut_pi) {
        rval = get_pi (lp, mynode_pi, mycut_pi);
        if (rval) {
            fprintf (stderr, "get_pi failed\n");
            goto CLEANUP;
        }
        *node_pi = mynode_pi;
        *cut_pi = mycut_pi;
    }
    if (lb) *lb = lp->result.lb;
    if (ub) *ub = lp->result.ub;
    if (ecount) *ecount = lp->result.ecount;

    return 0;

CLEANUP:

    CC_IFFREE (myelist, int);
    CC_IFFREE (myx, double);
    CC_IFFREE (myrc, double);
    CC_IFFREE (mynode_pi, double);
    CC_IFFREE (mycut_pi, double);

    return rval;
}

void CCtsp_add_cuts_to_queue (CCtsp_lp *lp, CCtsp_lpcut_in **clist)
{
    CCtsp_lpcut_in *c, *cnext;

    for (c = *clist; c; c = cnext) {
        cnext = c->next;
        add_cut_to_queue (lp, c);
    }
    *clist = (CCtsp_lpcut_in *) NULL;
}

static void add_cut_to_queue (CCtsp_lp *lp, CCtsp_lpcut_in *c)
{
    assert (c->sense == 'G' || c->branch != 0);
    c->next = &lp->cutqueue;
    c->prev = lp->cutqueue.prev;
    c->next->prev = c;
    c->prev->next = c;
}

double CCtsp_cutprice (CCtsp_lpgraph *g, CCtsp_lpcut_in *c, double *x)
{
    double slack;
    int nzlist, nznext;

    slack =  (double) -(c->rhs);

    nzlist = CCtsp_lpcut_in_nzlist (g, c);

    while (nzlist != -1) {
        nznext = g->edges[nzlist].coefnext;
        g->edges[nzlist].coefnext = -2;
        slack += g->edges[nzlist].coef * x[nzlist];
        g->edges[nzlist].coef = 0;
        nzlist = nznext;
    }
    return slack;
}

/* Domino cuts are handled by first adding all semi-cuts, then */
/* running through nzlist and incrementing the coef of every   */
/* odd edge that is not in delta(H), then running through      */
/* delta(H) and incrementing every even edge.                  */

int CCtsp_lpcut_in_nzlist (CCtsp_lpgraph *g, CCtsp_lpcut_in *c)
{
    int nzlist = -1;
    int i;

    if (c->dominocount == 0) {
        for (i = 0; i < c->cliquecount; i++) {
            lpcut_nonzero_work (g, &c->cliques[i], &nzlist);
        }
    } else {
        if (c->cliquecount != 1) {
            fprintf (stderr, "Yipes! Domino without a handle\n");
            exit (1);
        }
        for (i = 0; i < c->dominocount; i++) {
            lpcut_nonzero_semicut (g, &c->dominos[i], &nzlist);
        }
        lpcut_nonzero_parity_handle (g, &c->cliques[0], &nzlist);
        for (i = 0; i < c->dominocount; i++) {
            lpcut_nonzero_domino (g, &c->dominos[i], &nzlist);
        }
    }
    return nzlist;
}

int CCtsp_lpcut_nzlist (CCtsp_lpgraph *g, CCtsp_lpcut *c,
        CCtsp_lpclique *cliques, CCtsp_lpdomino *dominos, int do_mods)
{
    int nzlist = -1;
    int i;

    if (c->dominocount == 0) {
        for (i = 0; i < c->cliquecount; i++) {
            lpcut_nonzero_work (g, &cliques[c->cliques[i]], &nzlist);
        }
    } else {
        if (c->cliquecount != 1) {
            fprintf (stderr, "Yipes! Domino without a handle\n");
            exit (1);
        }
        if (!c->dominos) {
            fprintf (stderr, "Yipes! Dominocount with no dominos\n");
            exit (1);
        }
        for (i = 0; i < c->dominocount; i++) {
            lpcut_nonzero_semicut (g, &dominos[c->dominos[i]], &nzlist);
        }
        lpcut_nonzero_parity_handle (g, &cliques[c->cliques[0]], &nzlist);
        for (i = 0; i < c->dominocount; i++) {
            lpcut_nonzero_domino (g, &dominos[c->dominos[i]], &nzlist);
        }
    }
    if (do_mods) {
        lpcut_nonzero_modify (g, c->modcount, c->mods, &nzlist);
    }
    return nzlist;
}

static void lpcut_nonzero_work (CCtsp_lpgraph *g, CCtsp_lpclique *c,
        int *pnzlist)
{
    int nzlist = *pnzlist;
    int nodemarker;
    CCtsp_lpadj *a;
    int e;
    int tmp, k, l;

    g->nodemarker++;
    nodemarker = g->nodemarker;

    CC_FOREACH_NODE_IN_CLIQUE (k, *c, tmp) {
        g->nodes[k].mark = nodemarker;
    }

    CC_FOREACH_NODE_IN_CLIQUE (k, *c, tmp) {
        a = g->nodes[k].adj;
        for (l=0; l<g->nodes[k].deg; l++) {
            if (g->nodes[a[l].to].mark != nodemarker) {
                e = a[l].edge;
                if (g->edges[e].coefnext == -2) {
                    g->edges[e].coefnext = nzlist;
                    nzlist = e;
                }
                g->edges[e].coef++;
            }
        }
    }
    *pnzlist = nzlist;
}

static void lpcut_nonzero_semicut (CCtsp_lpgraph *g, CCtsp_lpdomino *d,
        int *pnzlist)
{
    /* Add in the nonzeros for the edges in the semi-cut (A:B) */

    int nzlist = *pnzlist;
    int nodemarker;
    CCtsp_lpadj *a;
    int e;
    int tmp, k, l;

    g->nodemarker++;
    nodemarker = g->nodemarker;

    CC_FOREACH_NODE_IN_CLIQUE (k, d->sets[0], tmp) {
        g->nodes[k].mark = nodemarker;
    }

    CC_FOREACH_NODE_IN_CLIQUE (k, d->sets[1], tmp) {
        a = g->nodes[k].adj;
        for (l=0; l<g->nodes[k].deg; l++) {
            if (g->nodes[a[l].to].mark == nodemarker) {
                e = a[l].edge;
                if (g->edges[e].coefnext == -2) {
                    g->edges[e].coefnext = nzlist;
                    nzlist = e;
                }
                g->edges[e].coef++;
            }
        }
    }
    *pnzlist = nzlist;
}

static void lpcut_nonzero_parity_handle (CCtsp_lpgraph *g, CCtsp_lpclique *H,
        int *pnzlist)
{
    /* Add in the nonzeros for the edges F, defined as the edges in an */
    /* odd number of semi-cuts (so an odd value in nzlist) and not in  */
    /* delta(H), or in an even number of semi-cuts and in delta(H).    */

    int nzlist = *pnzlist;
    int nodemarker;
    CCtsp_lpadj *a;
    int e;
    int tmp, k, l;

    g->nodemarker++;
    nodemarker = g->nodemarker;

    CC_FOREACH_NODE_IN_CLIQUE (k, *H, tmp) {
        g->nodes[k].mark = nodemarker;
    }

    for (e = nzlist; e != -1; e = g->edges[e].coefnext) {
        if (g->edges[e].coef % 2 == 1) {
            if ((g->nodes[g->edges[e].ends[0]].mark != nodemarker &&
                 g->nodes[g->edges[e].ends[1]].mark != nodemarker) || 
                (g->nodes[g->edges[e].ends[0]].mark == nodemarker &&
                 g->nodes[g->edges[e].ends[1]].mark == nodemarker)) {
                    g->edges[e].coef++;
            }
        }
    }

    CC_FOREACH_NODE_IN_CLIQUE (k, *H, tmp) {
        a = g->nodes[k].adj;
        for (l=0; l<g->nodes[k].deg; l++) {
            if (g->nodes[a[l].to].mark != nodemarker) {
                e = a[l].edge;
                if (g->edges[e].coefnext == -2) {
                    g->edges[e].coefnext = nzlist;
                    nzlist = e;
                    g->edges[e].coef++;
                } else if (g->edges[e].coef % 2 == 0) {
                    g->edges[e].coef++;
                }
            }
        }
    }
    *pnzlist = nzlist;
}

static void lpcut_nonzero_domino (CCtsp_lpgraph *g, CCtsp_lpdomino *d,
        int *pnzlist)
{
    /* Add in the nonzeros for the edges in the cut delta(A union B) */

    int nzlist = *pnzlist;
    int nodemarker;
    CCtsp_lpadj *a;
    int e;
    int tmp, i, k, l;

    g->nodemarker++;
    nodemarker = g->nodemarker;

    for (i = 0; i < 2; i++) {
        CC_FOREACH_NODE_IN_CLIQUE (k, d->sets[i], tmp) {
            g->nodes[k].mark = nodemarker;
        }
    }

    for (i = 0; i < 2; i++) {
        CC_FOREACH_NODE_IN_CLIQUE (k, d->sets[i], tmp) {
            a = g->nodes[k].adj;
            for (l=0; l<g->nodes[k].deg; l++) {
                if (g->nodes[a[l].to].mark != nodemarker) {
                    e = a[l].edge;
                    if (g->edges[e].coefnext == -2) {
                        g->edges[e].coefnext = nzlist;
                        nzlist = e;
                    }
                    g->edges[e].coef++;
                }
            }
        }
    }
    *pnzlist = nzlist;
}

static void lpcut_nonzero_modify (CCtsp_lpgraph *g, int modcount,
        CCtsp_sparser *mods, int *pnzlist)
{
    int nzlist = *pnzlist;
    int i,j,k,l;
    CCtsp_lpadj *a;
    int e;

    for (i=0; i<modcount; i++) {
        k = mods[i].node;
        j = ((int) mods[i].mult) - 128;
        a = g->nodes[k].adj;
        for (l=0; l<g->nodes[k].deg; l++) {
            e = a[l].edge;
            if (g->edges[e].coefnext == -2) {
                g->edges[e].coefnext = nzlist;
                nzlist = e;
            }
            g->edges[e].coef += j;
        }
    }
    *pnzlist = nzlist;
}

int CCtsp_process_cuts (CCtsp_lp *lp, int *pnadded, int tighten,
        int silent, CCrandstate *rstate)
{
    int nadded = 0;
    int addpiece = 0;
    CCtsp_lpcut_in *c = (CCtsp_lpcut_in *) NULL;
    double *x = (double *) NULL;
    int ecount_save;
    int rval;
    CCtsp_lprow cr;

/*
    printf ("CCtsp_process_cuts (%d, %d)\n", tighten, silent); fflush (stdout);
*/

    CCutil_start_timer (&lp->stats.addcuts);
    
    *pnadded = 0;
    CCtsp_init_lprow (&cr);
    ecount_save = lp->graph.ecount;
    x = CC_SAFE_MALLOC (ecount_save, double);
    if (!x) {
        return 1;
    }

    rval = CClp_x (lp->lp, x);
    if (rval) {
        CC_FREE (x, double);
        return rval;
    }

    nadded = 0;
    addpiece = 0;
    while (lp->cutqueue.next != &lp->cutqueue) {
        c = lp->cutqueue.next;
        c->next->prev = c->prev;
        c->prev->next = c->next;
        rval = checkout_cut (lp, c, x, &cr, tighten);
        if (rval > 1) {
            rval = 1;
            goto CLEANUP;
        } else if (rval == 1) {
            nadded++;
            addpiece++;
            if (cr.rowcnt >= CCtsp_STORE_BATCH) {
                if (cr.rowcnt > 0) {
                    rval = CCtsp_add_multiple_rows (lp, &cr);
                    if (rval) {
                        fprintf (stderr, "CCtsp_add_multiple_rows failed\n");
                        goto CLEANUP;
                    }
                    CCtsp_free_lprow (&cr);
                }
            }
            if (addpiece >= CCtsp_CUT_BATCH) {
                if (cr.rowcnt > 0) {
                    rval = CCtsp_add_multiple_rows (lp, &cr);
                    if (rval) {
                        fprintf (stderr, "CCtsp_add_multiple_rows failed\n");
                        goto CLEANUP;
                    }
                    CCtsp_free_lprow (&cr);
                }
                rval = update_newcuts (lp, silent, rstate);
                if (rval == 2) {
                    printf ("LP is really infeasible (processs_cuts)\n");
                    fflush (stdout);
                    goto CLEANUP;
                } else if (rval) {
                    goto CLEANUP;
                }
                if (lp->graph.ecount != ecount_save) {
                    CC_FREE (x, double);
                    ecount_save = lp->graph.ecount;
                    x = CC_SAFE_MALLOC (ecount_save, double);
                    if (!x) {
                        rval = 1;
                        goto CLEANUP;
                    }
                }
                rval = CClp_x (lp->lp, x);
                if (rval) {
                    goto CLEANUP;
                }
                addpiece = 0;
            }
        }
        CCtsp_free_lpcut_in (c);
        CC_FREE (c, CCtsp_lpcut_in);
    }
    if (addpiece > 0) {
        if (cr.rowcnt > 0) {
            rval = CCtsp_add_multiple_rows (lp, &cr);
            if (rval) {
                fprintf (stderr, "CCtsp_add_multiple_rows failed\n");
                goto CLEANUP;
            }
            CCtsp_free_lprow (&cr);
        }
        rval = update_newcuts (lp, silent, rstate);
        if (rval == 2) {
            printf ("LP is really infeasible (processs_cuts)\n");
            fflush (stdout);
            goto CLEANUP;
        } else if (rval) {
            goto CLEANUP;
        }
    }
    CC_FREE (x, double);
    *pnadded = nadded;

    rval = CCtsp_update_result (lp);
    if (rval) goto CLEANUP;

    CCutil_stop_timer (&lp->stats.addcuts, 0);
    return 0;

CLEANUP:

    CCtsp_free_lprow (&cr);
    if (c) {
        CCtsp_free_lpcut_in (c);
        CC_FREE (c, CCtsp_lpcut_in);
    }
    CC_IFFREE (x, double);
    *pnadded = nadded;
    return rval;
}

void CCtsp_init_lprow (CCtsp_lprow *cr)
{
    cr->rowcnt = 0;
    cr->nzcnt = 0;
    cr->sense = (char *) NULL;
    cr->rhs = (double *) NULL;
    cr->begin = (int *) NULL;
    cr->indexspace = 0;
    cr->indices = (int *) NULL;
    cr->entryspace = 0;
    cr->entries = (double *) NULL;
}

void CCtsp_free_lprow (CCtsp_lprow *cr)
{
    if (cr) {
        cr->rowcnt = 0;
        cr->nzcnt = 0;
        CC_IFFREE (cr->sense, char);
        CC_IFFREE (cr->rhs, double);
        CC_IFFREE (cr->begin, int);
        CC_IFFREE (cr->indices, int);
        cr->indexspace = 0;
        CC_IFFREE (cr->entries, double);
        cr->entryspace = 0;
    }
}

static int checkout_cut (CCtsp_lp *lp, CCtsp_lpcut_in *c, double *x,
        CCtsp_lprow *cr, int tighten)
{
    int rval;
    CCtsp_lpcut_in d;
    double slack;
    CCtsp_lpgraph *g = &lp->graph;

    if (tighten && c->dominocount == 0) {
        rval = CCtsp_tighten_lpcut_in (g, c, x, &d, &lp->stats.tighten_stats,
                                       (double *) NULL);
        if (rval) {
            fprintf (stderr, "CCtsp_tighten_lpcut_in failed\n");
            return 2;
        }
    } else {
        rval = CCtsp_copy_lpcut_in (c, &d);
        if (rval) {
            fprintf (stderr, "CCtsp_copy_lpcut_in failed\n");
            return 2;
        }
    }

    slack = CCtsp_cutprice (g, &d, x);
    if (slack >= -CCtsp_MIN_VIOL) {
/*
        printf ("Slack = %f  ", slack); fflush (stdout); 
        CCtsp_print_lpcut_in (&d);
*/
        CCtsp_free_lpcut_in (&d);
        return 0;
    }

#if 0
    rval = CCverify_cut (&d, CC_TYPE_ALL, (int *) NULL);
    if (rval) {
        fprintf (stderr, "Discarding invalid cut\n");
        CCtsp_print_lpcut_in (&d);
        CCtsp_free_lpcut_in (&d);
        return 2;
    }
#endif

    rval = CCtsp_add_cut (lp, &d, cr);
    if (rval) {
        CCtsp_free_lpcut_in (&d);
        return 2;
    }
    CCtsp_free_lpcut_in (&d);
    return 1;
}

static int update_newcuts (CCtsp_lp *lp, int silent, CCrandstate *rstate)
{
    int rval;
    int ndeleted;
    double objval;

/*
    printf ("update_newcuts (%d)\n", silent); fflush (stdout);
*/

    CCutil_start_timer (&lp->stats.addcuts_opt);
    rval = CClp_opt (lp->lp, CClp_METHOD_DUAL);
    if (rval == 2) {
        rval = CCtsp_infeas_recover (lp, silent, rstate);
        if (rval == 2) {
            if (!silent) {
                printf ("Problem is really infeasible (update_newcuts)\n");
                fflush (stdout);
                CCutil_stop_timer (&lp->stats.addcuts_opt, 1);
            } else {
                CCutil_stop_timer (&lp->stats.addcuts_opt, 0);
            }
            return 2;
        } else if (rval) {
            return 1;
        }
    } else if (rval) {
        fprintf (stderr, "CClp_opt failed\n");
        return 1;
    }
    CCutil_stop_timer (&lp->stats.addcuts_opt, 0);

    rval = CClp_objval (lp->lp, &objval);
    if (rval) {
        fprintf (stderr, "CClp_objval failed\n");
        return 1;
    }

    if (objval > lp->result.lb + 0.000001) {
        CCutil_start_timer (&lp->stats.agecuts);
        ndeleted = 0;
        rval = age_cuts (lp, &ndeleted);
        if (rval) {
            fprintf (stderr, "age_cuts failed\n");
            return 1;
        }
        CCutil_stop_timer (&lp->stats.agecuts, 0);

        CCutil_start_timer (&lp->stats.agecuts_opt);
        if ( ndeleted ) {
            rval = CClp_opt (lp->lp, CClp_METHOD_DUAL);
        }
        if (rval) {
            fprintf (stderr, "CClp_opt failed\n");
            return 1;
        }
        CCutil_stop_timer (&lp->stats.agecuts_opt, 0);

        CCutil_start_timer (&lp->stats.ageedges);
        ndeleted = 0;
        rval = age_edges (lp, &ndeleted);
        if (rval) {
            fprintf (stderr, "age_edges failed\n");
            return 1;
        }
        CCutil_stop_timer (&lp->stats.ageedges, 0);

        CCutil_start_timer (&lp->stats.ageedges_opt);
        if ( ndeleted ) {
            rval = CClp_opt (lp->lp, CClp_METHOD_DUAL);
        }
        if (rval) {
            fprintf (stderr, "CClp_opt failed\n");
            return 1;
        }
        CCutil_stop_timer (&lp->stats.ageedges_opt, 0);
    }

    return 0;
}

int CCtsp_infeas_recover (CCtsp_lp *lp, int silent, CCrandstate *rstate)
{
    double penalty;
    int nadded, feasible;
    int rval;

    printf ("infeas_recover ...\n"); fflush (stdout);

    rval = CCtsp_addbad_variables (lp, (CCtsp_edgegenerator *) NULL, &penalty,
            &nadded, CCtsp_PHASE1_RCTHRESH, CCtsp_PHASE1_MAXPENALTY, 1,
            &feasible, silent, rstate);
    if (rval) {
        fprintf (stderr, "CCtsp_addbad_variables failed\n");
        return 1;
    }

    if (feasible) {
        printf ("Recovered a feasible LP\n");
        fflush (stdout);
        return 0;
    } else {
        printf ("Could not recover a feasible LP\n");
        fflush (stdout);
        return 2;
    }
}

int CCtsp_add_cut (CCtsp_lp *lp, CCtsp_lpcut_in *d, CCtsp_lprow *cr)
{
    int nzlist;
    CCtsp_lpgraph *g = &lp->graph;
    int saved;
    CCtsp_lpcut new;
    int rval = 0;
    int newloc;
    int rhs;
    int i;

    assert (d->branch == 'G' || d->branch != 0);

    CCtsp_init_lpcut (&new);

    new.rhs         = d->rhs;
    new.sense       = d->sense;
    new.branch      = d->branch;
    rval = CCtsp_register_cliques (&lp->cuts, d, &new);
    CCcheck_rval (rval, "CCtsp_register_cliques failed");
    rval = CCtsp_register_dominos (&lp->cuts, d, &new);
    CCcheck_rval (rval, "CCtsp_register_dominos failed");

    nzlist = CCtsp_lpcut_in_nzlist (g, d);
    rval = CCtsp_qsparsify (&lp->sparsifier, g, &nzlist, &new.modcount,
            &new.mods, &saved);
    if (rval) {
        fprintf (stderr, "CCtsp_qsparsify failed\n");
        CCtsp_unregister_cliques (&lp->cuts, &new);
        CCtsp_unregister_dominos (&lp->cuts, &new);
        CC_IFFREE (new.mods, CCtsp_sparser);
        clear_nzlist (g, nzlist);
        goto CLEANUP;
    }
    new.age = CCtsp_NEWCUT_AGE;
    rval = CCtsp_copy_skeleton (&d->skel, &new.skel);
    if (rval) {
        fprintf (stderr, "CCtsp_copy_skeleton failed\n");
        CCtsp_unregister_cliques (&lp->cuts, &new);
        CCtsp_unregister_dominos (&lp->cuts, &new);
        CC_IFFREE (new.mods, CCtsp_sparser);
        clear_nzlist (g, nzlist);
        goto CLEANUP;
    }
        

    newloc = CCtsp_add_cut_to_cutlist (&lp->cuts, &new);
    if (newloc == -1) {
        fprintf (stderr, "CCtsp_add_cut_to_cutlist failed\n");
        CCtsp_unregister_cliques (&lp->cuts, &new);
        CCtsp_unregister_dominos (&lp->cuts, &new);
        CC_IFFREE (new.mods, CCtsp_sparser);
        CCtsp_free_skeleton (&new.skel);
        clear_nzlist (g, nzlist);
        rval = 1;  goto CLEANUP;
    }
    rhs = new.rhs;
    for (i=0; i<new.modcount; i++) {
        rhs += 2*(((int) new.mods[i].mult) - 128);
    }
    rval = CCtsp_add_nzlist_to_lp (lp, nzlist, rhs, new.sense, cr);
    if (rval) {
        fprintf (stderr, "CCtsp_add_nzlist_to_lp failed\n");
        CCtsp_delete_cut_from_cutlist (&lp->cuts, newloc);
        goto CLEANUP;
    }

CLEANUP:

    return rval;
}

int CCtsp_add_nzlist_to_lp (CCtsp_lp *lp, int nzlist, int rhs, char sense,
                            CCtsp_lprow *cr)
{
    int nzcnt;
    int i;
    CCtsp_lpgraph *g = &lp->graph;
    int *rmatind = (int *) NULL;
    double *rmatval = (double *) NULL;
    double drhs = (double) rhs;
    int rval;

    nzcnt = 0;
    for (i=nzlist; i != -1; i = g->edges[i].coefnext) {
        if (g->edges[i].coef) nzcnt++;
    }

    if (nzcnt != 0) {
        rmatind = CC_SAFE_MALLOC (nzcnt, int);
        if (!rmatind) {
            clear_nzlist (g, nzlist);
            return 1;
        }
        rmatval = CC_SAFE_MALLOC (nzcnt, double);
        if (!rmatval) {
            CC_FREE (rmatind, int);
            clear_nzlist (g, nzlist);
            return 1;
        }
        for (nzcnt = 0; nzlist != -1; nzlist = i) {
            i = g->edges[nzlist].coefnext;
            g->edges[nzlist].coefnext = -2;
            if (g->edges[nzlist].coef) {
                rmatind[nzcnt] = nzlist;
                rmatval[nzcnt] = g->edges[nzlist].coef;
                g->edges[nzlist].coef = 0;
                nzcnt++;
            }
        }
    } else {
        printf ("WARNING: Adding an empty cut to the LP\n");
        fflush (stdout);
    }

    rval = addrow_to_list (nzcnt, drhs, sense, rmatind, rmatval, cr);

    CC_FREE (rmatind, int);
    CC_FREE (rmatval, double);

    return rval;
}

static void clear_nzlist (CCtsp_lpgraph *g, int nzlist)
{
    int nznext;

    while (nzlist != -1) {
        nznext = g->edges[nzlist].coefnext;
        g->edges[nzlist].coefnext = -2;
        g->edges[nzlist].coef = 0;
        nzlist = nznext;
    }
}

int CCtsp_addbad_variables (CCtsp_lp *lp, struct CCtsp_edgegenerator *eg,
        double *ppenalty, int *pnadded, double rcthresh, double maxpenalty,
        int phase1, int *feasible, int silent, CCrandstate *rstate)
{
    int nadded;
    double penalty;
    int finished;
    CCutil_edgehash eh;
    int rval;
    int *genlist = (int *) NULL;
    int *genlen = (int *) NULL;
    int gencount;
    CCtsp_predge *inlist = (CCtsp_predge *) NULL;
    int incount;
    CCtsp_predge *prlist = (CCtsp_predge *) NULL;
    int prcount;
    double *node_pi = (double *) NULL;
    double *node_piest = (double *) NULL;
    double *clique_pi = (double *) NULL;
    double *domino_pi = (double *) NULL;
    double *cut_pi = (double *) NULL;
    int i, iend, n;
    int start = 0, ni, nj;
    int ngen = PRICE_GEN + PRICE_GEN_FACTOR * lp->graph.ncount;

    CCutil_start_timer (&lp->stats.addbad);

    if (phase1) {
        printf ("phase 1 addbad_variables\n"); fflush (stdout);
    }

    if (feasible)
        *feasible = 0;

    rval = CCutil_edgehash_init (&eh, (int) (ngen * 1.5));
    CCcheck_rval (rval, "CCutil_edgehash_init failed");

    gencount = 0;
    incount = 0;
    prcount = 0;

    genlist    = CC_SAFE_MALLOC (ngen*2, int);
    genlen     = CC_SAFE_MALLOC (ngen, int);
    inlist     = CC_SAFE_MALLOC (ngen, CCtsp_predge);
    prlist     = CC_SAFE_MALLOC (PRICE_POOL + ngen, CCtsp_predge);
    node_pi    = CC_SAFE_MALLOC (lp->graph.ncount, double);
    node_piest = CC_SAFE_MALLOC (lp->graph.ncount, double);
    if (!genlist || !genlen || !inlist || !prlist || !node_pi || !node_piest) {
        fprintf (stderr, "out of memory in CCtsp_addbad_variables\n");
        rval = 1; goto CLEANUP;
    }

    if (lp->cuts.cliqueend) {
        clique_pi = CC_SAFE_MALLOC (lp->cuts.cliqueend, double);
        CCcheck_NULL (clique_pi, "out of memory in CCtsp_addbad_variables");
    }
    if (lp->cuts.dominoend) {
        domino_pi = CC_SAFE_MALLOC (lp->cuts.dominoend, double);
        CCcheck_NULL (clique_pi, "out of memory in CCtsp_addbad_variables");
    }
    if (lp->cuts.cutcount) {
        cut_pi = CC_SAFE_MALLOC (lp->cuts.cutcount, double);
        CCcheck_NULL (cut_pi, "out of memory in CCtsp_addbad_variables");
    }

    *pnadded = 0;
    rval = pricing_duals (lp, node_pi, node_piest, cut_pi, clique_pi,
                          domino_pi);
    CCcheck_rval (rval, "pricing_duals failed");

    if (phase1) {
        start = ni = 0;
        nj = (lp->full_edges_valid ? 0 : start + 1);
    } else {
        rval = CCtsp_reset_edgegenerator (eg, node_piest, silent);
        CCcheck_rval (rval, "CCtsp_reset_edgegenerator failed");
    }
    finished = 0;
    nadded = 0;
    penalty = 0.0;
    while (!finished) {
        if (phase1) {
            rval = phase1_generate_edges (lp, node_piest, ngen,
                 &gencount, genlist, genlen, start, &ni, &nj, &finished);
            CCcheck_rval (rval, "phase1_generate_edges failed");
        } else {
            rval = CCtsp_generate_edges (eg, ngen, &gencount, genlist, genlen,
                                   &finished, silent, rstate);
            CCcheck_rval (rval, "CCtsp_generate_edges failed");
        }

        for (i = 0, incount = 0; i < gencount; i++) {
            int hval;
            if (CCutil_edgehash_find (&eh, genlist[2*i], genlist[2*i+1], &hval)
                != 0
                && CCtsp_find_edge (&lp->graph, genlist[2*i], genlist[2*i+1])
                == -1) {
                rval = CCutil_edgehash_add (&eh, genlist[2*i], genlist[2*i+1],
                                            1);
                if (rval) {
                    fprintf (stderr, "CCutil_edgehash_add failed\n");
                    goto CLEANUP;
                }
                if (genlist[2*i] < genlist[2*i+1]) {
                    inlist[incount].ends[0] = genlist[2*i];
                    inlist[incount].ends[1] = genlist[2*i+1];
                } else {
                    inlist[incount].ends[0] = genlist[2*i+1];
                    inlist[incount].ends[1] = genlist[2*i];
                }
                inlist[incount].len = genlen[i];
                incount++;
            }
        }
        rval = price_list (lp, incount, inlist, node_pi, cut_pi, clique_pi,
                           domino_pi, phase1);
        CCcheck_rval (rval, "price_list failed");

        for (i=0; i<incount; i++) {
            if (inlist[i].rc < 0.0) penalty += inlist[i].rc;
            if (inlist[i].rc < rcthresh) {
                prlist[prcount++] = inlist[i];
            } else {
                rval = CCutil_edgehash_del (&eh, inlist[i].ends[0],
                                           inlist[i].ends[1]);
                if (rval) goto CLEANUP;
            }
        }
        nadded = 0;
        while ((!finished && prcount >= PRICE_POOL) ||
                (finished && penalty < -maxpenalty && prcount > 0)) {
            n = PRICE_ADD;
            if (n >= prcount) {
                n = prcount;
            } else {
                pr_select (prcount - n, prcount, prlist);
            }
            rval = CCtsp_add_vars_to_lp (lp, prlist + prcount - n, n);
            if (rval) goto CLEANUP;
            nadded += n;
            *pnadded += n;
            prcount -= n;
            CCutil_start_timer (&lp->stats.addbad_opt);
            rval = CClp_opt (lp->lp, CClp_METHOD_DUAL);
            if (!silent) {
                CCutil_stop_timer (&lp->stats.addbad_opt, 1);
            } else {
                CCutil_stop_timer (&lp->stats.addbad_opt, 0);
            }
            if (phase1) {
                if (rval == 0) {
                    printf ("LP is now feasible\n"); fflush (stdout);
                    if (feasible) *feasible = 1;
                    goto DONE;
                } else if (rval != 2) {
                    fprintf (stderr, "CClp_opt failed\n");
                    goto CLEANUP;
                }
            } else {
                if (rval == 2) {
                    fprintf (stderr, "Adding variables made LP infeasible!\n");
                    rval = 1;
                    goto CLEANUP;
                } else if (rval) {
                    fprintf (stderr, "CClp_opt failed\n");
                    goto CLEANUP;
                }
            }
            rval = pricing_duals (lp, node_pi, node_piest, cut_pi,
                                  clique_pi, domino_pi);
            CCcheck_rval (rval, "pricing_duals failed");

            rval = price_list (lp, prcount, prlist, node_pi, cut_pi, clique_pi,
                               domino_pi, phase1);
            CCcheck_rval (rval, "price_list failed");

            penalty = 0.0;
            for (i=0, iend = prcount, prcount = 0; i<iend; i++) {
                if (prlist[i].rc < 0.0) penalty += prlist[i].rc;
                if (prlist[i].rc < rcthresh &&
                    CCtsp_find_edge (&lp->graph, prlist[i].ends[0],
                                               prlist[i].ends[1]) == -1) {
                    prlist[prcount++] = prlist[i];
                }
            }
        }
        if (nadded > 0) {
            if (phase1) {
                start = ni;
                nj = (lp->full_edges_valid ? 0 : start + 1);
            } else {
                rval = CCtsp_reset_edgegenerator (eg, node_piest, silent);
                if (rval) goto CLEANUP;
            }
            finished = 0;
            CCutil_edgehash_delall (&eh);
            for (i = 0; i < prcount; i++) {
                CCutil_edgehash_add (&eh, prlist[i].ends[0],
                                     prlist[i].ends[1], 1);
            }
        }
    }

DONE:

    if (!silent) {
        CCutil_stop_timer (&lp->stats.addbad, 1);
    } else {
        CCutil_stop_timer (&lp->stats.addbad, 0);
    }
    *ppenalty = penalty;
    if (!phase1 || (feasible && *feasible)) {
        rval = CCtsp_update_result (lp);
    }

CLEANUP:

    CC_IFFREE (cut_pi, double);
    CC_IFFREE (domino_pi, double);
    CC_IFFREE (clique_pi, double);
    CC_IFFREE (node_piest, double);
    CC_IFFREE (node_pi, double);
    CC_IFFREE (prlist, CCtsp_predge);
    CC_IFFREE (inlist, CCtsp_predge);
    CC_IFFREE (genlen, int);
    CC_IFFREE (genlist, int);
    CCutil_edgehash_free (&eh);
    return rval;
}

static int phase1_generate_edges (CCtsp_lp *lp, double *node_piest, int nwant,
       int *ngot, int *genlist, int *genlen, int start, int *ni, int *nj,
       int *finished)
{
    int i = *ni;
    int j = *nj;
    CCtsp_genadj *adj;
    CCdatagroup *dat;
    int cnt = 0;
    int ncount = lp->graph.ncount;

    *ngot = 0;
    *finished = 0;

    if (!lp->dat && !lp->full_edges_valid) {
        fprintf (stderr, "no source of edges to generate\n");
        return 1;
    }

    if (i >= ncount) {
        i = 0;
        j = (lp->full_edges_valid ? 0 : i + 1);
    }

    if (lp->full_edges_valid) {
        adj = lp->fulladj;
        for (; j < adj[i].deg; j++) {
            if (phase1_test_edge (i, adj[i].list[j].end, node_piest)) {
                genlist[2*cnt] = i;
                genlist[2*cnt+1] = adj[i].list[j].end;
                genlen[cnt] = adj[i].list[j].len;
                cnt++;
                if (cnt == nwant) {
                    goto NOT_FINISHED;
                }
            }
        }
        while ((i = (i+1) % ncount) != start) {
            for (j = 0; j < adj[i].deg; j++) {
                if (phase1_test_edge (i, adj[i].list[j].end, node_piest)) {
                    genlist[2*cnt] = i;
                    genlist[2*cnt+1] = adj[i].list[j].end;
                    genlen[cnt] = adj[i].list[j].len;
                    cnt++;
                    if (cnt == nwant) {
                        goto NOT_FINISHED;
                    }
                }
            }
        }
    } else {
        dat = lp->dat;
        for (; j < ncount; j++) {
            if (phase1_test_edge (i, j, node_piest)) {
                genlist[2*cnt] = i;
                genlist[2*cnt+1] = j;
                genlen[cnt] = CCutil_dat_edgelen (i, j, dat);
                cnt++;
                if (cnt == nwant) {
                    goto NOT_FINISHED;
                }
            }
        }
        while ((i = (i+1) % ncount) != start) {
            for  (j = i + 1; j < ncount; j++) {
                if (phase1_test_edge (i, j, node_piest)) {
                    genlist[2*cnt] = i;
                    genlist[2*cnt+1] = j;
                    genlen[cnt] = CCutil_dat_edgelen (i, j, dat);
                    cnt++;
                    if (cnt == nwant) {
                        goto NOT_FINISHED;
                    }
                }
            }
        }
    }

    *finished = 1;
    *ngot = cnt;
    return 0;

NOT_FINISHED:

    *finished = 0;
    *ngot = cnt;
    *ni = i; *nj = j + 1;
    return 0;
}

static int phase1_test_edge (int end1, int end2, double *node_piest)
{
    return (node_piest[end1] + node_piest[end2] > 0.0);
}

int CCtsp_eliminate_variables (CCtsp_lp *lp, int eliminate_sparse, int silent)
{
    int i, j, k;
    int rval;
    CCbigguy ub;

    /* This function should not create an infeasible LP */

    if (lp->upperbound == CCtsp_LP_MAXDOUBLE ||
            CCbigguy_cmp (lp->exact_lowerbound, CCbigguy_MINBIGGUY) == 0) {
        printf ("Can't elmininate without upper and lower bounds\n");
        fflush (stdout);
        return 0;
    }

    ub = CCbigguy_dtobigguy (lp->upperbound - 1.0);
    if (CCbigguy_cmp (lp->exact_lowerbound, ub) > 0) {
        printf ("No need for elimination, bounds are optimal\n");
        fflush (stdout);
        return 0;
    }

    assert (lp->infeasible == 0);

    rval = CCtsp_edge_elimination (lp, eliminate_sparse, silent);
    if (rval) {
        fprintf (stderr, "tsp_edge_elimination failed\n");
        return rval;
    }

    assert (lp->infeasible == 0);

    for (i = 0; i < lp->nfixededges; i++) {
        k = CCtsp_find_edge (&(lp->graph), lp->fixededges[2*i],
                                           lp->fixededges[2*i+1]);
        if (k != -1) {
            rval = CClp_setbnd (lp->lp, k, 'L', 1.0);
            lp->graph.edges[k].fixed = 1;
        } else {
            printf ("WARNING: Fixed edge (%d, %d) is not in LP\n",
                     lp->fixededges[2*i], lp->fixededges[2*i+1]);
            fflush (stdout);
        }
    }

    CC_IFFREE (lp->graph.adjspace, CCtsp_lpadj);
    for (i = lp->graph.ecount - 1; i >= 0; i--) {
        if (!find_edge_full (lp, lp->graph.edges[i].ends[0],
                                 lp->graph.edges[i].ends[1])) {
            if (!lp->graph.edges[i].fixed && !lp->graph.edges[i].branch) {
                rval = CClp_delete_column (lp->lp, i);
                if (rval) {
                    fprintf (stderr, "CClp_delete_column failed\n");
                    return rval;
                }
                lp->graph.edges[i].ends[0] = 0;
                lp->graph.edges[i].ends[1] = 0;
            } else {
                printf ("WARNING: Tried to eliminate a fixed/branch edge\n");
                fflush (stdout);
            }
        }
    }
    for (i = 0, j = 0; i < lp->graph.ecount; i++) {
        if (lp->graph.edges[i].ends[1] != 0 ||
            lp->graph.edges[i].ends[0] != 0) {
            lp->graph.edges[j] = lp->graph.edges[i];
            j++;
        }
    }
    if (!silent) {
        printf ("Eliminated %d LP edges\n", lp->graph.ecount - j);
        fflush (stdout);
    }

    assert (lp->infeasible == 0);

    lp->graph.ecount = j;
    rval = CCtsp_build_lpadj (&lp->graph, 0, lp->graph.ecount);
    if (rval) {
        fprintf (stderr, "CCtsp_build_lpadj failed\n");
        return rval;
    }

    rval = CClp_opt (lp->lp, CClp_METHOD_DUAL);
    if (rval == 2) {
        fprintf (stderr, "ERROR: edge_elimination created an infeasible LP\n");
        return 1;
    } else if (rval) {
        fprintf (stderr, "CClp_opt failed\n");
        return rval;
    }
    rval = CCtsp_update_result (lp);
    if (rval) {
        fprintf (stderr, "CCtsp_update_result failed\n");
        return rval;
    }

    assert (lp->infeasible == 0);

    return 0;
}

static int pricing_duals (CCtsp_lp *lp, double *node_pi, double *node_piest,
        double *cut_pi, double *clique_pi, double *domino_pi)
{
    double x;
    int i, j, tmp, k;
    int rval = 0;

    rval = get_pi (lp, node_pi, cut_pi);
    CCcheck_rval (rval, "get_pi failed");

    for (i = 0; i < lp->cuts.cliqueend; i++) {
        clique_pi[i] = 0.0;
    }

    for (i = 0; i < lp->cuts.dominoend; i++) {
        domino_pi[i] = 0.0;
    }

    for (i = 0; i < lp->cuts.cutcount; i++) {
        x = cut_pi[i];
        for (k = 0; k < lp->cuts.cuts[i].modcount; k++) {
            node_pi[lp->cuts.cuts[i].mods[k].node] += x *
                (((int) lp->cuts.cuts[i].mods[k].mult) - 128);
        }
        for (k = 0; k < lp->cuts.cuts[i].cliquecount; k++) {
            clique_pi[lp->cuts.cuts[i].cliques[k]] += x;
        }
        for (k = 0; k < lp->cuts.cuts[i].dominocount; k++) {
            domino_pi[lp->cuts.cuts[i].dominos[k]] += x;
        }
    }

    for (i = 0; i < lp->graph.ncount; i++) {
        node_piest[i] = node_pi[i];
    }

    for (i = 0; i < lp->cuts.cliqueend; i++) {
        x = clique_pi[i];
        if (x > 0.0) {
            CC_FOREACH_NODE_IN_CLIQUE (k, lp->cuts.cliques[i], tmp) {
                node_pi[k] += x;
                node_piest[k] += x;
            }
        } else if (x < 0.0) {
            CC_FOREACH_NODE_IN_CLIQUE (k, lp->cuts.cliques[i], tmp) {
                node_pi[k] += x;
            }
        }
    }

    for (i = 0; i < lp->cuts.dominoend; i++) {
        x = 2.0 * domino_pi[i];
        if (x > 0.0) {
            for (j = 0; j < 2; j++) {
                CC_FOREACH_NODE_IN_CLIQUE (k,lp->cuts.dominos[i].sets[j],tmp) {
                    node_pi[k] += x;
                    node_piest[k] += x;
                }
            }
        } else if (x < -0.001) {
            fprintf (stderr, "YIPES: We have a negative domino: %f\n", x);
            rval = 1;  goto CLEANUP;
        }
    }

    /* For dominos, we want to add pi for each node in the handle and */
    /* 2*pi for each node in each of the dominos.                     */

    /* For now, no munching of edges.  To do in future, just add
       loop here for each edge clique, (if find_edge,)? zero the edge
       clique's pi */

CLEANUP:

    return rval;
}

static int price_list (CCtsp_lp *lp, int ecount, CCtsp_predge *elist,
        double *node_pi, double *cut_pi, double *clique_pi, double *domino_pi,
        int phase1)
{
    CCtsp_lpadj *adjspace = (CCtsp_lpadj *) NULL;
    CCtsp_lpnode *n = (CCtsp_lpnode *) NULL;
    int i, j, tmp, l, k, ci, nzlist, nznext;
    CCtsp_lpadj *a;
    int marker = 0;
    double x;
    int ncount = lp->graph.ncount;
    int ccount = lp->cuts.cliqueend;
    CCtsp_lpclique *c = lp->cuts.cliques;
    double *nodom_pi = (double *) NULL;
    double *nodom_clique_pi = (double *) NULL;
    CCtsp_lpgraph g;
    int *temp_elist = (int *) NULL;
    int rval = 0;

    CCtsp_init_lpgraph_struct (&g);

    if (ecount == 0) goto CLEANUP;

    /* remove domino contributions from node_pi to get nodom_pi */

    nodom_pi = CC_SAFE_MALLOC (ncount, double);
    CCcheck_NULL (nodom_pi, "out of memory in price_list");
    for (i = 0; i < ncount; i++) {
        nodom_pi[i] = node_pi[i];
    }
    for (i = 0; i < lp->cuts.dominoend; i++) {
        x = 2.0 * domino_pi[i];
        if (x > 0.0) {
            for (j = 0; j < 2; j++) {
                CC_FOREACH_NODE_IN_CLIQUE (k,lp->cuts.dominos[i].sets[j],tmp) {
                    nodom_pi[k] -= x;
                }
            }
        } else if (x < -0.001) {
            fprintf (stderr, "YIPES: We have a negative domino: %f\n", x);
            rval = 1; goto CLEANUP;
        }
    }

    /* remove domino handles from clique_pi to get nodom_clique_pi */
    /* and remove the handle contributions from nodom_pi           */

    if (ccount > 0) {
        nodom_clique_pi = CC_SAFE_MALLOC (ccount, double);
        CCcheck_NULL (nodom_clique_pi, "out of memory in price_list");
    }
    for (i = 0; i < ccount; i++) {
        nodom_clique_pi[i] = clique_pi[i];
    }
    if (lp->cuts.dominoend > 0) {
        for (i = 0; i < lp->cuts.cutcount; i++) {
            if (lp->cuts.cuts[i].dominocount > 0) {
                if (lp->cuts.cuts[i].cliquecount != 1) {
                    fprintf (stderr, "domoino with no handle\n");
                    rval = 1; goto CLEANUP;
                }
                ci = lp->cuts.cuts[i].cliques[0];
                x = cut_pi[i];
                if (x > 0.0) {
                    nodom_clique_pi[ci] -= x;
                    CC_FOREACH_NODE_IN_CLIQUE (k, lp->cuts.cliques[ci], tmp) {
                        nodom_pi[k] -= x;
                    }
                } else if (x < -0.001) {
                    fprintf (stderr, "YIPES: negative domino %f\n", x);
                    rval = 1; goto CLEANUP;
                }
            }
        }
    }

    /* build adjlist to make it easy to run through the cliques and */
    /* adjust the reduced costs; if dominos look promising, this    */
    /* should be changed to work with an lpadj, since we need to    */
    /* build this in any case to call nzlist for the domino cuts    */

    n = CC_SAFE_MALLOC (ncount, CCtsp_lpnode);
    CCcheck_NULL (n, "out of memory in price_list");

    adjspace = CC_SAFE_MALLOC (2*ecount, CCtsp_lpadj);
    CCcheck_NULL (adjspace, "out of memory in price_list");


    for (i = 0; i < ncount; i++) {
        n[i].deg = 0;
        n[i].mark = 0;
    }
    for (i = 0; i < ecount; i++) {
        if (phase1) {
            elist[i].rc = - nodom_pi[elist[i].ends[0]]
                          - nodom_pi[elist[i].ends[1]];

        } else {
            elist[i].rc = elist[i].len
                          - nodom_pi[elist[i].ends[0]]
                          - nodom_pi[elist[i].ends[1]];
        }

        n[elist[i].ends[0]].deg++;
        n[elist[i].ends[1]].deg++;
    }
    a = adjspace;
    for (i = 0; i < ncount; i++) {
        n[i].adj = a;
        a += n[i].deg;
        n[i].deg = 0;
    }
    for (i=0; i<ecount; i++) {
        j = elist[i].ends[0];
        n[j].adj[n[j].deg].to = elist[i].ends[1];
        n[j].adj[n[j].deg].edge = i;
        n[j].deg++;
        j = elist[i].ends[1];
        n[j].adj[n[j].deg].to = elist[i].ends[0];
        n[j].adj[n[j].deg].edge = i;
        n[j].deg++;
    }

    for (i = 0; i < ccount; i++) {
        if (nodom_clique_pi[i]) {
            x = nodom_clique_pi[i] * 2;
            marker++;
            CC_FOREACH_NODE_IN_CLIQUE (j, c[i], tmp) {
                a = n[j].adj;
                for (l = 0; l < n[j].deg; l++) {
                    if (n[a[l].to].mark == marker) {
                        elist[a[l].edge].rc += x;
                        /* We could test if rc>0, and delete e from adj*/
                    }
                }
                n[j].mark = marker;
            }
        }
    }

    /*  Run through domino-parity cuts, subtracting dual values from rc */

    if (lp->cuts.dominoend > 0) {
        temp_elist = CC_SAFE_MALLOC (2*ecount, int);
        CCcheck_NULL (temp_elist, "out of memory in price_list");
        for (i = 0; i < ecount; i++) {
            temp_elist[2*i]   = elist[i].ends[0];
            temp_elist[2*i+1] = elist[i].ends[1];
        }
        rval = CCtsp_build_lpgraph (&g, ncount, ecount, temp_elist,
                                   (int *) NULL);
        CCcheck_rval (rval, "CCtsp_build_lpgraph failed");
        rval = CCtsp_build_lpadj (&g, 0, ecount);
        CCcheck_rval (rval, "CCtsp_build_lpadj failed");

        for (i = 0; i < lp->cuts.cutcount; i++) {
            x = cut_pi[i];
            if (lp->cuts.cuts[i].dominocount > 0 && x > 0.0) {
                nzlist = CCtsp_lpcut_nzlist (&g, &(lp->cuts.cuts[i]),
                                       lp->cuts.cliques, lp->cuts.dominos, 0);
                while (nzlist != -1) {
                    nznext = g.edges[nzlist].coefnext;
                    g.edges[nzlist].coefnext = -2;
                    if (g.edges[nzlist].coef) {
                        elist[nzlist].rc -=
                              (((double) g.edges[nzlist].coef) * x);
                        g.edges[nzlist].coef = 0;
                    }
                    nzlist = nznext;
                }
            }
        }
    }

CLEANUP:

    CC_IFFREE (n, CCtsp_lpnode);
    CC_IFFREE (adjspace, CCtsp_lpadj);
    CC_IFFREE (nodom_pi, double);
    CC_IFFREE (nodom_clique_pi, double);
    CC_IFFREE (temp_elist, int);
    CCtsp_free_lpgraph (&g);
    return rval;
}

static void pr_select (int nsel, int n, CCtsp_predge *list)
{
    double s1, s2, s3, s;
    int i, j, k, l;
    CCtsp_predge t;

    s1 = list[n/4].rc;
    s2 = list[n/2].rc;
    s3 = list[3*n/4].rc;
    if (s1 < s2) {
        if (s2 < s3) s = s2;
        else if (s1 < s3) s = s3;
        else s = s1;
    } else {
        if (s1 < s3) s = s1;
        else if (s2 < s3) s = s3;
        else s = s2;
    }

    /* 0 <= x < i ==> r[i] == s
     * i <= x < j ==> r[i] > s
     * j <= x <= k ==> r[i] unknown
     * k < x <= l ==> r[i] < s
     * l < x < n ==> r[i] == s
     */

    i = j = 0;
    k = l = n-1;
    while (j <= k) {
        while (j <= k && list[j].rc >= s) {
            if (list[j].rc == s) {
                CC_SWAP (list[i], list[j], t);
                i++;
            }
            j++;
        }
        while (j <= k && list[k].rc <= s) {
            if (list[k].rc == s) {
                CC_SWAP (list[l], list[k], t);
                l--;
            }
            k--;
        }
        if (j <= k) {
            CC_SWAP (list[j], list[k], t);
            j++;
            k--;
        }
    }
    while (i > 0) {
        i--;
        j--;
        CC_SWAP (list[i], list[j], t);
    }
    while (l < n-1) {
        k++;
        l++;
        CC_SWAP (list[k], list[l], t);
    }
    /* 0 <= x < j ==> r[i] > s
     * j <= x <= k ==> r[i] == s
     * k < x < n ==> r[i] < s
     */
    if (nsel < j) {
        pr_select (nsel, j, list);
    } else if (nsel > k) {
        pr_select (nsel - (k+1), n - (k+1), list + (k+1));
    }
}

int CCtsp_add_vars_to_lp (CCtsp_lp *lp, CCtsp_predge *prlist, int n)
{
    int i;
    CCtsp_lpedge *e;
    int rval;
    int nzcnt;
    double *obj;
    int *matbeg;
    int *matind;
    double *matval;
    double *lb;
    double *ub;

    while (lp->graph.ecount + n > lp->graph.espace) {
        if (CCutil_reallocrus_scale ((void **) &lp->graph.edges,
                &lp->graph.espace, lp->graph.ecount + n, 1.3,
                sizeof (CCtsp_lpedge))) {
            return 1;
        }
    }
    e = lp->graph.edges + lp->graph.ecount;
    for (i = 0; i < n; i++) {
        e[i].ends[0] = prlist[i].ends[0];
        e[i].ends[1] = prlist[i].ends[1];
        e[i].fixed = 0;
        e[i].branch = 0;
        e[i].age = 0;
        e[i].len = prlist[i].len;
        e[i].coefnext = -2;
        e[i].coef = 0;
    }

    rval = build_lp_cols (&lp->graph, &lp->cuts, lp->graph.ecount,
                          lp->graph.ecount + n, &nzcnt, &obj, &matbeg,
                          (int **) NULL, &matind, &matval, &lb, &ub);
    if (rval) return rval;

    rval = lp_addcols (lp, n, nzcnt, obj, matbeg, matind, matval,
                          lb, ub);
    if (rval) goto CLEANUP;

    lp->graph.ecount += n;
    rval = CCtsp_build_lpadj (&lp->graph, 0, lp->graph.ecount);
    if (rval) goto CLEANUP;
    rval = 0;

CLEANUP:

    CC_IFFREE (obj, double);
    CC_IFFREE (matbeg, int);
    CC_IFFREE (matind, int);
    CC_IFFREE (matval, double);
    CC_IFFREE (lb, double);
    CC_IFFREE (ub, double);
    return rval;
}

static int age_cuts (CCtsp_lp *lp, int *ndeleted)
{
    int *del = (int *) NULL;
    double *cut_pi = (double *) NULL;
    int ncount = lp->graph.ncount;
    int rval;
    int i, j;
    CClp_info *b = (CClp_info *) NULL;

/*
    printf ("age_cuts () ...\n"); fflush (stdout); 
*/

    *ndeleted = 0;

    rval = CClp_get_info (lp->lp, &b);
    if (rval) {
        fprintf (stderr, "CClp_get_info failed\n");
        goto CLEANUP;
    }

    if (lp->cuts.cutcount) {
        cut_pi = CC_SAFE_MALLOC (lp->cuts.cutcount, double);
        if (!cut_pi) {
            fprintf (stderr, "Out of memory in age_cuts\n");
            rval = 1; goto CLEANUP;
        }
    }
    rval = get_pi (lp, (double *) NULL, cut_pi);
    if (rval) {
        fprintf (stderr, "get_pi failed\n");
        goto CLEANUP;
    }

    for (i = lp->cuts.cutcount-1; i >= 0  /* 5000 DFS HACK */; i--) {
        if (lp->cuts.cuts[i].branch == 0) {
            if ((cut_pi[i] > CCtsp_DUAL_DUST ||
                 lp->cuts.cuts[i].age == CCtsp_NEWCUT_AGE) &&
                CClp_is_row_active (b, ncount+i)) {
                if (lp->cuts.cuts[i].age == CCtsp_NEWCUT_AGE) {
                    if (lp->cuts.cuts[i].dominocount == 0) {
                        rval = CCtsp_add_to_cutpool (lp->pool, &lp->cuts,
                                                &lp->cuts.cuts[i]);
                        CCcheck_rval (rval, "CCtsp_add_to_cutpool failed");
                    } else {
                        rval = CCtsp_add_to_dominopool (lp->dominopool,
                                                &lp->cuts, &lp->cuts.cuts[i]);
                        CCcheck_rval (rval, "CCtsp_add_to_dominopool failed");
                    }
                }
                lp->cuts.cuts[i].age = 0;
            } else {
                if (lp->cuts.cuts[i].age != CCtsp_NEWCUT_AGE
                 /* && !CClp_is_row_active (b, ncount+i)  no-pivotin HACK */) {
                    lp->cuts.cuts[i].age++;
                }
                if (lp->cuts.cuts[i].age == CCtsp_NEWCUT_AGE || 
                    lp->cuts.cuts[i].age >= lp->cut_life) {
                    CCtsp_unregister_cliques (&lp->cuts, &lp->cuts.cuts[i]);
                    CCtsp_unregister_dominos (&lp->cuts, &lp->cuts.cuts[i]);
                    CC_IFFREE (lp->cuts.cuts[i].mods, CCtsp_sparser);
                    CCtsp_free_skeleton (&lp->cuts.cuts[i].skel);
                    lp->cuts.cuts[i].cliquecount = 0;
                    lp->cuts.cuts[i].dominocount = 0;
                    lp->cuts.cuts[i].modcount = 0;
                }
            }
        }
    }
    CC_FREE (cut_pi, double);
    CClp_free_info (&b);

    del = CC_SAFE_MALLOC (lp->cuts.cutcount, int);
    if (!del) {
        fprintf (stderr, "out of memory in age_cuts\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0, j = 0; i < lp->cuts.cutcount; i++) {
        if (lp->cuts.cuts[i].cliquecount) {
            lp->cuts.cuts[j] = lp->cuts.cuts[i];
            j++;
            del[i] = 0;
        } else {
            del[i] = 1;
        }
    }
    if (j < lp->cuts.cutcount) {
        rval = lp_delete_cut_set (lp, del);
        if (rval) {
            fprintf (stderr, "lp_delete_cut_set failed\n");
            goto CLEANUP;
        }
    }

    *ndeleted = lp->cuts.cutcount - j;
    lp->cuts.cutcount = j;
    rval = 0;

  CLEANUP:
    CC_IFFREE (cut_pi, double);
    CC_IFFREE (del, int);
    CClp_free_info (&b);
    return rval;
}

static int age_edges (CCtsp_lp *lp, int *ndeleted)
{
    int *del = (int *) NULL;
    int rval;
    int i, j;
    CClp_info *b = (CClp_info *) NULL;
    double *x = (double *) NULL;

    *ndeleted = 0;

    rval = CClp_get_info (lp->lp, &b);
    if (rval) {
        fprintf (stderr, "CClp_get_info failed\n");
        goto CLEANUP;
    }

    CC_IFFREE (lp->graph.adjspace, CCtsp_lpadj);

    x = CC_SAFE_MALLOC (lp->graph.ecount, double);
    if (x == (double *) NULL) {
        fprintf (stderr, "Out of memory in age_edges\n");
        rval = 1; goto CLEANUP;
    }

    rval = CClp_x (lp->lp, x);
    if (rval) {
        fprintf (stderr, "CClp_x failed\n");
        goto CLEANUP;
    }

    for (i = lp->graph.ecount-1; i >= 0; i--) {
        if (lp->graph.edges[i].fixed == 0 &&
            lp->graph.edges[i].branch == 0) {
            if (x[i] > CCtsp_EDGE_DUST && CClp_is_col_active (b, i)) {
                lp->graph.edges[i].age = 0;
            } else {
                lp->graph.edges[i].age++;
/*
                if (lp->graph.edges[i].age >= lp->edge_life) {
                Bico:  29.5.00 - don't delete basis columns  HACK
*/
                if (lp->graph.edges[i].age >= lp->edge_life &&
                    CClp_is_col_active (b, i) == 0) {
                    lp->graph.edges[i].ends[0] = 0;
                    lp->graph.edges[i].ends[1] = 0;
                }
            }
        }
    }

    CClp_free_info (&b);
    CC_FREE (x, double);

    del = CC_SAFE_MALLOC (lp->graph.ecount, int);
    if (!del) {
        fprintf (stderr, "out of memory in age_edges\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0, j = 0; i < lp->graph.ecount; i++) {
        if (lp->graph.edges[i].ends[1] != 0 ||
            lp->graph.edges[i].ends[0] != 0) {
            lp->graph.edges[j] = lp->graph.edges[i];
            j++;
            del[i] = 0;
        } else {
            del[i] = 1;
        }
    }
    if (j < lp->graph.ecount) {
        rval = lp_delete_var_set (lp, del);
        if (rval) {
            fprintf (stderr, "lp_delete_var_set failed\n");
            rval = 1; goto CLEANUP;
        }
    }

    *ndeleted = lp->graph.ecount - j;
    lp->graph.ecount = j;
    rval = CCtsp_build_lpadj (&lp->graph, 0, lp->graph.ecount);
    if (rval) {
        fprintf (stderr, "CCtsp_build_lpadj failed\n");
        goto CLEANUP;
    }

    rval = 0;

  CLEANUP:
    CC_IFFREE (del, int);
    CC_IFFREE (x, double);
    CClp_free_info (&b);
    return rval;
}

static int get_pi (CCtsp_lp *lp, double *node_pi, double *cut_pi)
{
    int rval;
    double *pi = (double *) NULL;
    int ncount = lp->graph.ncount;
    int ncuts = lp->cuts.cutcount;
    int nrows = ncount + ncuts;
    int i;

    pi = CC_SAFE_MALLOC (nrows, double);
    if (pi == (double *) NULL) {
        fprintf (stderr, "Out of memory in get_pi\n");
        rval = 1; goto CLEANUP;
    }

    rval = CClp_pi (lp->lp, pi);
    if (rval) {
        fprintf (stderr, "CClp_pi failed\n");
        goto CLEANUP;
    }


    if (node_pi) {
        for (i = 0; i < ncount; i++) {
            node_pi[i] = pi[i];
        }
    }
    if (cut_pi) {
        for (i = 0; i < ncuts; i++) {
            cut_pi[i] = pi[ncount + i];
        }
    }

    rval = 0;
  CLEANUP:
    CC_IFFREE (pi, double);
    return rval;
}

int CCtsp_add_multiple_rows (CCtsp_lp *lp, CCtsp_lprow *cr)
{
    int rval = 0;

    if (!cr->rowcnt)
        return 0;

    rval = CClp_addrows (lp->lp, cr->rowcnt, cr->nzcnt, cr->rhs, cr->sense,
                         cr->begin, cr->indices, cr->entries);
    if (rval) {
        fprintf (stderr, "Couldn't add rows into LP\n");
    }
    return rval;
}

static int addrow_to_list (int nzcnt, double drhs, char sense,
            int *rmatind, double *rmatval, CCtsp_lprow *cr)
{
    int i, rval;
    int *ip;
    double *dp;
    double rhs = drhs; /* copy to work around problem with sun cc -dalign */

    rval = CCutil_reallocrus_count ((void **) &(cr->sense), cr->rowcnt + 1,
                             sizeof (char));
    if (rval) goto CLEANUP;
    rval = CCutil_reallocrus_count ((void **) &(cr->rhs), cr->rowcnt + 1,
                             sizeof (double));
    if (rval) goto CLEANUP;
    rval = CCutil_reallocrus_count ((void **) &(cr->begin), cr->rowcnt + 1,
                             sizeof (int));
    if (rval) goto CLEANUP;

    if (nzcnt + cr->nzcnt > cr->indexspace) {
        rval = CCutil_reallocrus_scale ((void **) &(cr->indices),
                &(cr->indexspace), nzcnt + cr->nzcnt, 1.3, sizeof (int));
        if (rval) goto CLEANUP;
    }
    if (nzcnt + cr->nzcnt > cr->entryspace) {
        rval = CCutil_reallocrus_scale ((void **) &(cr->entries),
                &(cr->entryspace), nzcnt + cr->nzcnt, 1.3, sizeof (double));
        if (rval) goto CLEANUP;
    }

    cr->sense[cr->rowcnt] = sense;
    cr->rhs[cr->rowcnt] = rhs;
    cr->begin[cr->rowcnt] = cr->nzcnt;
    cr->rowcnt++;

    ip = cr->indices + cr->nzcnt;
    dp = cr->entries + cr->nzcnt;
    for (i = 0; i < nzcnt; i++) {
        ip[i] = rmatind[i];
        dp[i] = rmatval[i];
    }
    cr->nzcnt += nzcnt;

    return 0;

CLEANUP:

    fprintf (stderr, "out of memory in addrow_to_list\n");
    return rval;
}

static int lp_addcols (CCtsp_lp *lp, int ncols, int nzcnt,
                  double *obj, int *matbeg, int *matind, double *matval,
                  double *lb, double *ub)
{
    int rval = 0;

    rval = CClp_addcols (lp->lp, ncols, nzcnt, obj,
                         matbeg, matind, matval, lb, ub);
    if (rval) {
        fprintf (stderr, "Couldn't add columns into LP\n");
    }
    return rval;
}

int CCtsp_delete_cut (CCtsp_lp *lp, int i)
{
    int rval = 0;
    rval = CClp_delete_row (lp->lp, lp->graph.ncount + i);
    if (rval) fprintf (stderr, "CClp_delete_row failed\n");
    return rval;
}

static int lp_delete_cut_set (CCtsp_lp *lp, int *del)
{
    int i, rval = 0;
    int *delstat = (int *) NULL;
    int ncount = lp->graph.ncount;
    int cutcount = lp->cuts.cutcount;

    delstat = CC_SAFE_MALLOC (ncount + cutcount, int);
    if (!delstat) {
        fprintf (stderr, "out of memory in lp_delete_cut_set\n");
        return 1;
    }
    for (i = 0; i < ncount; i++)
        delstat[i] = 0;
    for (i = 0; i < cutcount; i++)
        delstat[i + ncount] = del[i];

    rval = CClp_delete_set_of_rows (lp->lp, delstat);
    if (rval) {
        fprintf (stderr, "CClp_delete_set_of_rows failed\n");
    }

    CC_FREE (delstat, int);
    return rval;
}

int CCtsp_reduced_cost_nearest (CCtsp_lp *lp, int k, int *ecount, int **elist,
        double **elen, int sparse)
{
    CCdatagroup *dat = lp->dat;
    double *node_pi = (double *) NULL;
    double *node_piest = (double *) NULL;
    double *clique_pi = (double *) NULL;
    double *domino_pi = (double *) NULL;
    double *cut_pi = (double *) NULL;
    double *bestk  = (double *) NULL;
    double *bestkl = (double *) NULL;
    int    *bestki = (int *) NULL;
    int    *tlen = (int *) NULL;
    CCtsp_predge *inlist = (CCtsp_predge *) NULL;
    CCutil_edgehash eh;
    int incount;
    int ncount = lp->graph.ncount;
    int i, j, l;
    double d, e;
    int rval = 0;

    rval = CCutil_edgehash_init (&eh, (int) (k * ncount * 1.5));
    if (rval) return rval;

    if (sparse && !lp->fulladj) {
        printf ("Need fulladj to run sparse nearest\n");
        rval = 1;  goto CLEANUP;
    }

    node_pi    = CC_SAFE_MALLOC (lp->graph.ncount, double);
    node_piest = CC_SAFE_MALLOC (lp->graph.ncount, double);
    if (!node_pi || !node_piest) {
        fprintf (stderr, "out of memory in CCtsp_reduced_cost_nearest, pi\n");
        rval = 1; goto CLEANUP;
    }

    if (lp->cuts.cliqueend) {
        clique_pi = CC_SAFE_MALLOC (lp->cuts.cliqueend, double);
        CCcheck_NULL (clique_pi, "out of memory in CCtsp_reduced_cost_nearest");
    }
    if (lp->cuts.dominoend) {
        domino_pi = CC_SAFE_MALLOC (lp->cuts.dominoend, double);
        CCcheck_NULL (domino_pi, "out of memory in CCtsp_reduced_cost_nearest");
    }
    if (lp->cuts.cutcount) {
        cut_pi = CC_SAFE_MALLOC (lp->cuts.cutcount, double);
        CCcheck_NULL (cut_pi, "out of memory in CCtsp_reduced_cost_nearest");
    }

    rval = pricing_duals (lp, node_pi, node_piest, cut_pi, clique_pi,
                          domino_pi);
    CCcheck_rval (rval, "pricing_duals failed");
  
    if (sparse) {
        rval = sparse_nearest (lp, k, node_pi, cut_pi, clique_pi, domino_pi,
                               &eh);
        CCcheck_rval (rval, "sparse_nearest failed");
    } else {
        inlist     = CC_SAFE_MALLOC (ncount, CCtsp_predge);
        bestk      = CC_SAFE_MALLOC (k+1, double);
        bestkl     = CC_SAFE_MALLOC (k+1, double);
        bestki     = CC_SAFE_MALLOC (k+1, int);
        if (!inlist || !bestk || !bestkl || !bestki) {
            fprintf (stderr, "out of memory in CCtsp_reduced_cost_nearest\n");
            rval = 1; goto CLEANUP;
        }

        for (i = 0; i < ncount; i++) {
            if (i % 100 == 99) {
                printf ("."); fflush (stdout);
                if (i % 5000 == 4999) {
                    printf ("\n"); fflush (stdout);
                }
            }
            incount = 0;
            for (j = 0; j < ncount; j++) {
                if (j != i) {
                    inlist[incount].ends[0] = i;
                    inlist[incount].ends[1] = j;
                    inlist[incount++].len = CCutil_dat_edgelen (i, j, dat);
                }
            }
            rval = price_list (lp, incount, inlist, node_pi, cut_pi, clique_pi, 
                               domino_pi, 0);
            CCcheck_rval (rval, "price_list failed");

            for (l = 0; l < k; l++) {
                bestk[l] = CCtsp_LP_MAXDOUBLE;
                bestkl[l] = CCtsp_LP_MAXDOUBLE;
                bestki[l] = -1;
            }
            bestk[k]  = -CCtsp_LP_MAXDOUBLE;
            bestkl[k] = -CCtsp_LP_MAXDOUBLE;

            for (j = 0; j < incount; j++) {
                d = inlist[j].rc;
                e = inlist[j].len;
                if (d < bestk[0] || (d == bestk[0] && e < bestkl[0])) {
                    for (l = 0; d < bestk[l+1] ||
                          (d == bestk[l+1] && e < bestkl[l+1]); l++) {
                        bestk[l] = bestk[l+1];
                        bestkl[l] = bestkl[l+1];
                        bestki[l] = bestki[l+1];
                    }
                    bestk[l] = d;
                    bestkl[l] = e;
                    bestki[l] = j;
                }
            }        
            for (l = k-1; l >= 0; l--) {
                int hval;
                if (bestk[l] != CCtsp_LP_MAXDOUBLE) {
                    if (CCutil_edgehash_find (&eh, inlist[bestki[l]].ends[0],
                            inlist[bestki[l]].ends[1], &hval) != 0) {
                        rval = CCutil_edgehash_add (&eh, inlist[bestki[l]].ends[0],
                                                         inlist[bestki[l]].ends[1],
                                                         inlist[bestki[l]].len);
                        if (rval) {
                            fprintf (stderr, "CCutil_edgehash_add failed\n");
                            goto CLEANUP;
                        }
                    }
                }
            }
        }
    }

    if (sparse) {  /* Add the lp edges */
        int e1, e2, hval, len;

        for (i = 0; i < lp->graph.ecount; i++) {
            e1 = lp->graph.edges[i].ends[0];
            e2 = lp->graph.edges[i].ends[1];
            len = CCutil_dat_edgelen (e1, e2, dat);
            if (CCutil_edgehash_find (&eh, e1, e2, &hval) != 0) {
                printf ("+"); fflush (stdout);
                rval = CCutil_edgehash_add (&eh, e1, e2, len);
                if (rval) {
                    fprintf (stderr, "CCutil_edgehash_add failed\n");
                    goto CLEANUP;
                }
            }
        }
    }

    rval = CCutil_edgehash_getall (&eh, ecount, elist, &tlen);
    if (rval) {
        fprintf (stderr, "CCutil_edgehash_getall failed\n"); goto CLEANUP;
    }

    CC_IFFREE (inlist, CCtsp_predge);
    inlist = CC_SAFE_MALLOC (*ecount, CCtsp_predge);
    CCcheck_NULL (inlist, "out of memory in CCtsp_reduced_cost_nearest");

    for (i = 0; i < *ecount; i++) {
        inlist[i].ends[0] = (*elist)[2*i];
        inlist[i].ends[1] = (*elist)[2*i+1];
        inlist[i].len     = tlen[i];
    }

    rval = price_list (lp, *ecount, inlist, node_pi, cut_pi, clique_pi,
                       domino_pi, 0);
    CCcheck_rval (rval, "price_list failed");

    *elen = CC_SAFE_MALLOC (*ecount, double);
    CCcheck_NULL (*elen, "out of memory in CCtsp_reduced_cost_nearest");

    for (i = 0; i < *ecount; i++) {
        (*elen)[i] = inlist[i].rc;
    }

CLEANUP:

    CC_IFFREE (cut_pi, double);
    CC_IFFREE (domino_pi, double);
    CC_IFFREE (clique_pi, double);
    CC_IFFREE (node_piest, double);
    CC_IFFREE (node_pi, double);
    CC_IFFREE (bestk, double);
    CC_IFFREE (bestkl, double);
    CC_IFFREE (bestki, int);
    CC_IFFREE (tlen, int);
    CC_IFFREE (inlist, CCtsp_predge);
    CCutil_edgehash_free (&eh);

    return rval;
}


#define N_MAX_COUNT (10000)

static int sparse_nearest (CCtsp_lp *lp, int k, double *node_pi,
        double *cut_pi, double *clique_pi, double *domino_pi,
        CCutil_edgehash *eh)
{
    int rval = 0;
    int i, j, n0, n1, l, cnt, incount, mdeg;
    double d, e;
    int ncount = lp->graph.ncount;
    CCtsp_genadj *adj = lp->fulladj;
    CCdatagroup *dat = lp->dat;
    CCtsp_predge *inlist = (CCtsp_predge *) NULL;
    int **ilist = (int **) NULL;
    double **blist = (double **) NULL;
    double **llist = (double **) NULL;

    mdeg = 0;
    for (i = 0; i < ncount; i++) {
        if (adj[i].deg > mdeg) mdeg = adj[i].deg;
    }
    if (mdeg == 0) {
        fprintf (stderr, "no edges in fulladj\n");
        rval = 1;  goto CLEANUP;
    }

    inlist = CC_SAFE_MALLOC (mdeg * N_MAX_COUNT, CCtsp_predge);
    CCcheck_NULL (inlist, "out of memory in sparse_nearest");
    
    blist = CC_SAFE_MALLOC (ncount, double *);
    CCcheck_NULL (blist, "out of memory in sparse_nearest");
    for (i = 0; i < ncount; i++) {
        blist[i] = (double *) NULL;
    }

    llist = CC_SAFE_MALLOC (ncount, double *);
    CCcheck_NULL (llist, "out of memory in sparse_nearest");
    for (i = 0; i < ncount; i++) {
        llist[i] = (double *) NULL;
    }

    ilist = CC_SAFE_MALLOC (ncount, int *);
    CCcheck_NULL (ilist, "out of memory in sparse_nearest");
    for (i = 0; i < ncount; i++) {
        ilist[i] = (int *) NULL;
    }

    for (i = 0; i < ncount; i++) {
        blist[i] = CC_SAFE_MALLOC (k+1, double);
        llist[i] = CC_SAFE_MALLOC (k+1, double);
        ilist[i] = CC_SAFE_MALLOC (k+1, int);
        if (!blist[i] || !llist[i] || !ilist[i]) {
            fprintf (stderr, "out of memory in CCtsp_reduced_cost_nearest\n");
            rval = 1; goto CLEANUP;
        }
        for (l = 0; l < k; l++) {
            blist[i][l] = CCtsp_LP_MAXDOUBLE;
            llist[i][l] = CCtsp_LP_MAXDOUBLE;
            ilist[i][l] = -1;
        }
        blist[i][k]  = -CCtsp_LP_MAXDOUBLE;
        llist[i][k] = -CCtsp_LP_MAXDOUBLE;
    }

    i = 0;
    while (i < ncount) {
        printf ("Process %d\n", i); fflush (stdout);

        incount = 0;
        for (cnt = 0; cnt < N_MAX_COUNT && i < ncount; cnt++, i++) {
            for (j = 0; j < adj[i].deg; j++) {
                inlist[incount].ends[0] = i;
                inlist[incount].ends[1] = adj[i].list[j].end;
                inlist[incount++].len =
                          CCutil_dat_edgelen (i, adj[i].list[j].end, dat);
            }
        }

        rval = price_list (lp, incount, inlist, node_pi, cut_pi, clique_pi,
                           domino_pi, 0);
        CCcheck_rval (rval, "price_list failed");

        for (j = 0; j < incount; j++) {
            d = inlist[j].rc;
            e = inlist[j].len;
            n0 = inlist[j].ends[0];
            n1 = inlist[j].ends[1];
            pushsparse (n0, n1, d, e, blist, llist, ilist);
            pushsparse (n1, n0, d, e, blist, llist, ilist);
        }        
    }
    for (i = 0; i < ncount; i++) {
        for (l = k-1; l >= 0; l--) {
            int hval;
            if (blist[i][l] != CCtsp_LP_MAXDOUBLE) {
                if (CCutil_edgehash_find (eh, i, ilist[i][l], &hval) != 0) {
                    rval = CCutil_edgehash_add (eh, i, ilist[i][l], 
                                CCutil_dat_edgelen (i, ilist[i][l], dat));
                    if (rval) {
                        fprintf (stderr, "CCutil_edgehash_add failed\n");
                        goto CLEANUP;
                    }
                }
            }
        }
    }

CLEANUP:

    CC_IFFREE (inlist, CCtsp_predge);
    if (ilist) {
        for (i = 0; i < ncount; i++) {
            CC_IFFREE (ilist[i], int);
        }
        CC_FREE (ilist, int *);
    }
    if (llist) {
        for (i = 0; i < ncount; i++) {
            CC_IFFREE (llist[i], double);
        }
        CC_FREE (llist, double *);
    }
    if (blist) {
        for (i = 0; i < ncount; i++) {
            CC_IFFREE (blist[i], double);
        }
        CC_FREE (blist, double *);
    }

    return rval;
}

static void pushsparse (int i, int n, double d, double e, double **blist,
        double **llist, int **ilist)
{
    int l;

    if (d < blist[i][0] || (d == blist[i][0] && e < llist[i][0])) {
        for (l = 0; d < blist[i][l+1] ||
              (d == blist[i][l+1] && e < llist[i][l+1]); l++) {
            blist[i][l] = blist[i][l+1];
            llist[i][l] = llist[i][l+1];
            ilist[i][l] = ilist[i][l+1];
        }
        blist[i][l] = d;
        llist[i][l] = e;
        ilist[i][l] = n;
    }
}        

static int lp_delete_var_set (CCtsp_lp *lp, int *del)
{
    int rval = 0;

    rval = CClp_delete_set_of_columns (lp->lp, del);
    if (rval) {
        fprintf (stderr, "CClp_delete_set_of_columns failed\n");
    }
    return rval;
}

int CCtsp_write_probfile_sav (CCtsp_lp *lp)
{
    return write_probfile (lp, PROBTYPE_SAVE, (char *) NULL);
}


int CCtsp_write_probroot_id (char *problname, CCtsp_lp *lp)
{
    return write_probfile (lp, PROBTYPE_ROOT, problname);
}

int CCtsp_write_probleaf_id (CCtsp_lp *lp)
{
    return write_probfile (lp, PROBTYPE_LEAF, (char *) NULL);
}

int CCtsp_write_probfile_id (CCtsp_lp *lp)
{
    return write_probfile (lp, PROBTYPE_NORMAL, (char *) NULL);
}

static int write_probfile (CCtsp_lp *lp, int probtype, char *fname)
{
    CCtsp_PROB_FILE *p = (CCtsp_PROB_FILE *) NULL;
    char nambuf[1024];

    if (!lp) {
        fprintf (stderr, "write_probfile called without an lp\n");
        return 1;
    }
    if (!lp->graph.ecount) {
        fprintf (stderr, "write_probfile called with an edgeset\n");
        return 1;
    }

    if (probtype == PROBTYPE_SAVE) {
        if (strlen (lp->problabel) + 40 > sizeof (nambuf)) {
            sprintf (nambuf, "probfile.sav");
        } else {
            sprintf (nambuf, "%s.sav", lp->problabel);
        }

        p = CCtsp_prob_write_name (nambuf);
        if (!p) {
            fprintf (stderr, "CCtsp_prob_write_name failed\n");
            return 1;
        }
    } else {
        if (fname == (char *) NULL) {
            fname = lp->probloc;
        }
        p = CCtsp_prob_write (fname, lp->id);
        if (!p) {
            fprintf (stderr, "CCtsp_prob_write failed\n");
            return 1;
        }
    }

    if (CCtsp_prob_putname (p, lp->problabel)) {
        fprintf (stderr, "CCtsp_prob_putname failed\n"); goto CLEANUP;
    }
    if (CCtsp_prob_putid (p, lp->id)) {
        fprintf (stderr, "CCtsp_prob_putid failed\n"); goto CLEANUP;
    }
    if (CCtsp_prob_putparent (p, lp->parent_id)) {
        fprintf (stderr, "CCtsp_prob_putparent failed\n"); goto CLEANUP;
    }
    if (CCtsp_prob_putnnodes (p, lp->graph.ncount)) {
        fprintf (stderr, "CCtsp_prob_putnodes failed\n"); goto CLEANUP;
    }
    if (CCtsp_prob_putub (p, lp->upperbound))  {
        fprintf (stderr, "CCtsp_prob_putub failed\n"); goto CLEANUP;
    }
    if (CCtsp_prob_putlb (p, lp->lowerbound)) {
        fprintf (stderr, "CCtsp_prob_putlb failed\n"); goto CLEANUP;
    }
    if (CCtsp_prob_putexactlb (p, lp->exact_lowerbound)) {
        fprintf (stderr, "CCtsp_prob_puteactlb failed\n"); goto CLEANUP;
    }
    if (CCtsp_prob_putinfeasible (p, lp->infeasible))  {
        fprintf (stderr, "CCtsp_prob_putinfeasible failed\n"); goto CLEANUP;
    }

    if (probtype != PROBTYPE_LEAF) {
        int *elist = (int *) NULL;
        int *elen = (int *) NULL;
        int ecount = lp->graph.ecount;
        int i;

        elist = CC_SAFE_MALLOC (2 * ecount, int);
        elen = CC_SAFE_MALLOC (ecount, int);

        if (!elist || !elen) {
            fprintf (stderr, "out of memory in write_probfile\n");
            CC_IFFREE (elist, int);
            CC_IFFREE (elen, int);
            goto CLEANUP;
        }

        for (i = 0; i < ecount; i++) {
            elist[2*i]     = lp->graph.edges[i].ends[0];
            elist[2*i + 1] = lp->graph.edges[i].ends[1];
            elen[i]        = lp->graph.edges[i].len;
        }

        if (CCtsp_prob_putedges (p, lp->graph.ncount, ecount, elist, elen)) {
            fprintf (stderr, "CCtsp_prob_putedges failed\n");
            CC_FREE (elist, int);
            CC_FREE (elen, int);
            goto CLEANUP;
        }
        CC_FREE (elist, int);
        CC_FREE (elen, int);
    }

    if (CCtsp_prob_putcuts (p, lp->graph.ncount, &(lp->cuts))) {
        fprintf (stderr, "CCtsp_prob_putcuts failed\n");
        goto CLEANUP;
    }

    if (probtype != PROBTYPE_LEAF) {
        CClp_warmstart *w = (CClp_warmstart *) NULL;

        if (CClp_get_warmstart (lp->lp, &w)) {
            printf ("No warmstart to add to probfile\n");
            fflush (stdout);
        } else {
            if (CCtsp_prob_putwarmstart (p, w)) {
                fprintf (stderr, "CCtsp_prob_putwarmstart failed\n");
                CClp_free_warmstart (&w);
                goto CLEANUP;
            }
        }
        CClp_free_warmstart (&w);
    }

    if (probtype != PROBTYPE_LEAF && lp->nfixededges > 0) {
        if (CCtsp_prob_putfixed (p, lp->graph.ncount, lp->nfixededges,
                lp->fixededges)) {
            fprintf (stderr, "CCtsp_prob_putfixed failed\n");
            goto CLEANUP;
        }
    }

    if (probtype != PROBTYPE_LEAF && lp->fullcount > 0 &&
        lp->full_edges_valid) {
        if (CCtsp_prob_putfulladj (p, lp->graph.ncount, lp->fullcount,
                              lp->fulladj)) {
            fprintf (stderr, "CCtsp_prob_putfulladj failed\n");
            goto CLEANUP;
        }
    }

    if ((probtype == PROBTYPE_LEAF || probtype == PROBTYPE_ROOT) &&
        lp->exact_dual) {
        if (CCtsp_prob_putexactdual (p, lp->exact_dual, lp->graph.ncount)) {
            fprintf (stderr, "CCtsp_prob_putexact_dual failed\n");
            goto CLEANUP;
        }
    }

    if (CCtsp_prob_puthistory (p, lp->branchdepth, lp->branchhistory)) {
        fprintf (stderr, "CCtsp_prob_puthistory failed\n");
        goto CLEANUP;
    }

    if (CCtsp_prob_wclose (p)) {
        fprintf (stderr, "CCtsp_prob_wclose failed\n");
        return 1;
    }

    return 0;

CLEANUP:

    fprintf (stderr, "write_probfile failed\n");

    if (p)
        CCtsp_prob_wclose (p);
    return 1;
}

int CCtsp_read_probfile (CCtsp_lp *lp, char *fname, char *probloc, int *ncount,
        int silent)
{
    CCtsp_PROB_FILE *p = (CCtsp_PROB_FILE *) NULL;

    p = CCtsp_prob_read_name (fname);
    if (!p) {
        fprintf (stderr, "could not open %s for reading\n", fname);
        return 1;
    }

    lp->problabel = CCtsp_problabel (fname);
    if (lp->problabel == (char *) NULL) {
        fprintf (stderr, "CCtsp_problabel failed\n");
        CCtsp_prob_rclose (p);
        return 1;
    }

    if (probloc == (char *) NULL) {
        lp->probloc = CCutil_strdup (lp->problabel);
    } else {
        lp->probloc = CCutil_strdup (probloc);
    }
    
    return read_probfile (lp, p, ncount, silent);
}

int CCtsp_read_probfile_id (CCtsp_lp *lp, char *name, int id, int *ncount,
        int silent)
{
    CCtsp_PROB_FILE *p = (CCtsp_PROB_FILE *) NULL;

    p = CCtsp_prob_read (name, id);
    if (!p) {
        fprintf (stderr, "could not open %s for reading\n", name);
        return 1;
    }

    lp->probloc = CCutil_strdup (name);
    if (lp->probloc == (char *) NULL) {
        fprintf (stderr, "CCutil_strdup failed\n");
        CCtsp_prob_rclose (p);
        return 1;
    }
    
    return read_probfile (lp, p, ncount, silent);
}

static int read_probfile (CCtsp_lp *lp, CCtsp_PROB_FILE *p, int *ncount,
        int silent)
{
    int tncount, ecount;
    int *elist = (int *) NULL;
    int *elen = (int *) NULL;
    int i, k;
    int rval = 0;

    rval = CCtsp_prob_getnnodes (p, &tncount);
    if (rval == -1) goto CLEANUP;
    if (rval == 0) {
        if (ncount != (int *) NULL && *ncount != 0 && tncount != *ncount) {
            fprintf (stderr, "node counts differ in probfile and input\n");
            rval = 1; goto CLEANUP;
        }
        if (ncount != (int *) NULL && *ncount == 0) {
            *ncount = tncount;
        }
    } else {
        if (ncount == (int *) NULL || *ncount == 0) {
            fprintf (stderr, "node count not present in probfile or input\n");
            rval = 1; goto CLEANUP;
        } else {
            tncount = *ncount;
        }
    }

    rval = CCtsp_init_cliquehash (&lp->cuts, 2*tncount);
    if (rval) return rval;

    rval = CCtsp_init_dominohash (&lp->cuts, 2*tncount);
    if (rval) return rval;

    lp->problabel = CC_SAFE_MALLOC (CCtsp_PROB_FILE_NAME_LEN, char);
    if (!lp->problabel) {
        fprintf (stderr, "out of memory in read_probfile\n");
        rval = 1;
        goto CLEANUP;
    }
    rval = CCtsp_prob_getname (p, lp->problabel);
    if (rval == -1) goto CLEANUP;
    if (!silent) {
        printf ("Prob Name: %s\n", lp->problabel); fflush (stdout);
    }

    rval = CCtsp_prob_getid (p, &(lp->id));
    if (rval == -1) goto CLEANUP;
    if (!silent) {
        printf ("Prob ID: %d\n", lp->id); fflush (stdout);
    }

    rval = CCtsp_prob_getparent (p, &(lp->parent_id));
    if (rval == -1) goto CLEANUP;
    if (!silent) {
        printf ("Prob Parent ID: %d\n", lp->parent_id); fflush (stdout);
    }

    rval = CCtsp_prob_getub (p, &(lp->upperbound));
    if (rval == -1) goto CLEANUP;
    rval = CCtsp_prob_getlb (p, &(lp->lowerbound));
    if (rval == -1) goto CLEANUP;
    if (!silent) {
        printf ("Prob Bounds: (%f, %f)\n", lp->lowerbound, lp->upperbound);
        fflush (stdout);
    }

    rval = CCtsp_prob_getexactlb (p, &(lp->exact_lowerbound));
    if (rval == -1) goto CLEANUP;
    if (CCbigguy_cmp (lp->exact_lowerbound, CCbigguy_MINBIGGUY) != 0) {
        if (!silent) {
            printf ("Prob Exact Lowerbound: %f\n",
                             CCbigguy_bigguytod (lp->exact_lowerbound));
            fflush (stdout);
        }
    }

    rval = CCtsp_prob_getinfeasible (p, &(lp->infeasible));
    if (rval == -1) goto CLEANUP;
    if (lp->infeasible) {
        if (!silent) {
            printf ("Prob stored is tagged as infeasible\n"); fflush (stdout);
        }
    }

    rval = CCtsp_prob_getcuts (p, &tncount, &(lp->cuts), silent);
    if (rval == -1) goto CLEANUP;

    rval = CCtsp_prob_getedges (p, tncount, &ecount, &elist, &elen, silent);
    if (rval == -1) goto CLEANUP;
    if (!rval) {
        rval = CCtsp_build_lpgraph (&lp->graph, tncount, ecount, elist, elen);
        if (rval) goto CLEANUP;
        rval = CCtsp_build_lpadj (&lp->graph, 0, ecount);
        if (rval) goto CLEANUP;
        CC_FREE (elist, int);
        CC_FREE (elen, int);
    }

    rval = CCtsp_prob_getfixed (p, tncount, &(lp->nfixededges),
            &(lp->fixededges), silent);
    if (rval == -1) goto CLEANUP;
    if (!rval) {
        if (!silent) {
            printf ("Read %d LP fixed edges\n", lp->nfixededges);
            fflush (stdout);
        }
        for (i = 0; i < lp->nfixededges; i++) {
            k = CCtsp_find_edge (&(lp->graph), lp->fixededges[2*i],
                                               lp->fixededges[2*i+1]);
            if (k != -1) {
                lp->graph.edges[k].fixed = 1;
            } else {
                printf ("WARNING: File want's to fix a non-lp edge\n");
                fflush (stdout);
            }
        }
    }

    rval = CCtsp_prob_getwarmstart (p, &(lp->warmstart), silent);
    if (rval == -1) {
        fprintf (stderr, "CCtsp_prob_getwarmstart failed\n");
        goto CLEANUP;
    }

    rval = CCtsp_prob_getfulladj (p, tncount, &(lp->fullcount),
                             &(lp->fulladj), &(lp->fulladjspace), silent);
    if (rval == -1) {
        fprintf (stderr, "CCtsp_prob_getfulladj failed\n");
        goto CLEANUP;
    }
    if (!rval) {
        if (!silent) {
            printf ("Read LP full adj\n"); fflush (stdout);
        }
        if (lp->fullcount) {
            lp->full_edges_valid = 1;
        }
    }

    rval = CCtsp_prob_getexactdual (p, tncount, &(lp->exact_dual), silent);
    if (rval == -1) {
        fprintf (stderr, "CCtsp_prob_getexactdual failed\n");
        goto CLEANUP;
    }
    if (!rval) {
        if (!silent) {
            printf ("Read LP exact dual values\n"); fflush (stdout);
        }
    }

    rval = CCtsp_prob_gethistory (p, &lp->branchdepth, &lp->branchhistory,
                                  silent);
    if (rval == -1) {
        fprintf (stderr, "CCtsp_prob_gethistory failed\n");
        goto CLEANUP;
    }
    if (!rval) {
        if (!silent) {
            CCtsp_print_branchhistory (lp);
        }
    }

    rval = 0;


CLEANUP:

    if (CCtsp_prob_rclose (p)) {
        fprintf (stderr, "CCtsp_prob_rclose failed\n");
        return 1;
    }

    if (!silent) {
        printf ("Done with read_probfile\n"); fflush (stdout);
    }
    return rval;
}

int CCtsp_dump_rc_nearest (CCtsp_lp *lp, int k, char *fname, int sparse)
{
    int rc_count;
    int *rc_list = (int *) NULL;
    double *rc_len = (double *) NULL;
    int ncount = lp->graph.ncount;
    int i;
    int rval = 0;

    printf ("Dumping the %d-nearest rc vector to %s\n", k, fname);
    fflush (stdout);

    rval = CCtsp_reduced_cost_nearest (lp, k, &rc_count, &rc_list, &rc_len,
                                       sparse);
    if (rval) {
        fprintf (stderr, "CCtsp_reduced_cost_nearest failed\n");
        goto CLEANUP;
    }

    for (i = 0; i < rc_count; i++) {
        rc_list[2*i] = lp->perm[rc_list[2*i]];
        rc_list[2*i+1] = lp->perm[rc_list[2*i+1]];
    }

    rval = CCutil_writeedges_double (ncount, fname, rc_count, rc_list,
                                     rc_len, 0);
    if (rval) {
        fprintf (stderr, "CCutil_writeedges_int failed\n");
        fflush (stdout);
    }

CLEANUP:


    CC_IFFREE (rc_list, int);
    CC_IFFREE (rc_len, double);
    return rval;
}

int CCtsp_dump_x (CCtsp_lp *lp, char *fname)
{
    int i, xcount;
    int nonzero = 0;
    int *xlist = (int *) NULL;
    double *x  = (double *) NULL;
    FILE *out;
    int rval = 0;

    printf ("Dumping the x vector to %s ... ", fname); fflush (stdout);

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL, &xcount,
             &xlist, &x, (double **) NULL, (double **) NULL, (double **) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_get_lp_result failed\n"); return rval;
    }

    for (i = 0; i < xcount; i++) {
        if (x[i] > CCtsp_INTTOL) {
            nonzero++;
        }
    }

    printf ("%d edges, %d nonzero ", xcount, nonzero); fflush (stdout);
    
    out = fopen (fname, "w");
    if (out == (FILE *) NULL) {
        fprintf (stderr, "could not open %s for writing\n", fname);
        rval = 1; goto CLEANUP;
    }

    fprintf (out, "%d %d\n", lp->graph.ncount, nonzero);
    for (i = 0; i < xcount; i++) {
        if (x[i] > CCtsp_INTTOL) {
            fprintf (out, "%d %d %f\n", lp->perm[xlist[2*i]],
                                        lp->perm[ xlist[2*i+1]],
                                        x[i]);
/*
            fprintf (out, "%d %d %f\n", xlist[2*i],xlist[2*i+1], x[i]);
*/
        }
    }
    fclose (out);
    printf ("DONE\n"); fflush (stdout);

CLEANUP:

    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);
    return rval;
}

int CCtsp_depot_valid (CCtsp_lp *lp, int ndepot, int *yesno)
{
    int i, xcount, norig;
    int *xlist = (int *) NULL;
    double *x  = (double *) NULL;
    double total = 0.0;
    int rval = 0;

    if (yesno) *yesno = 0;

    norig = lp->graph.ncount - ndepot;

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL, &xcount,
             &xlist, &x, (double **) NULL, (double **) NULL, (double **) NULL);
    CCcheck_rval (rval, "CCtsp_get_lp_result failed");

    for (i = 0; i < xcount; i++) {
        if (xlist[2*i] >= norig && xlist[2*i+1] >= norig) {
            total += x[i];
        }
    }

    printf ("Depot Total: %lf\n", total); fflush (stdout);
    if (total >= 0.5) {
        printf ("Depot nodes are valid\n"); fflush (stdout);
        if (yesno) *yesno = 1;
    } else {
        printf ("SUB: WARNING - not enough depot nodes\n");
        if (yesno) *yesno = 0;
    }
    fflush (stdout);

CLEANUP:

    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);
    return rval;
}

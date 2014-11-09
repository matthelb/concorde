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
/*                  THE CONTROLLER FOR CUTTING PLANES                       */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: June 27, 1997                                                     */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  void CCtsp_init_cutselect (CCtsp_cutselect *s)                          */
/*    INITIALIZES the cut selections                                        */
/*                                                                          */
/*  void CCtsp_cutselect_dominos (CCtsp_cutselect *s, int domsel)           */
/*    SETS the domino field according to value of domsel                    */
/*                                                                          */
/*  void CCtsp_cutselect_tighten (CCtsp_cutselect *s, int tighten)          */
/*    SETS the usetighten field according to value of tighten               */
/*                                                                          */
/*  void CCtsp_cutselect_chunksize (CCtsp_cutselect *s, int chunksize)      */
/*    SETS the maxchunksize filed according to value of chunksize           */
/*                                                                          */
/*  void CCtsp_cutselect_filecuts (CCtsp_cutselect *s, char *fname)         */
/*    SETS the cutselector to read cuts from file fname                     */
/*                                                                          */
/*  void CCtsp_cutselect_remotepool (CCtsp_cutselect *s, char *cutbossname) */
/*    SETS the cutselector to use the specified cutboss                     */
/*                                                                          */
/*  void CCtsp_cutselect_domboss (CCtsp_cutselect *s, char *dombossname)    */
/*    SETS the cutselector to use the specified domino boss                 */
/*                                                                          */
/*  void CCtsp_init_tentative_cutselect (CCtsp_cutselect *s)                */
/*    INITIALIZES the cut selections for tenative branching                 */
/*                                                                          */
/*  int CCtsp_cutselect_set_tols (CCtsp_cutselect *s, CCtsp_lp *lp,         */
/*      int level, int silent)                                              */
/*    SETS the tolerances for the cut selections                            */
/*     -level should be set to 0 for tentative cutting                      */
/*    NOTES: The lp should be solved before this call.                      */
/*                                                                          */
/*  int CCtsp_cutting_multiple_loop (CCtsp_lp *lp, CCtsp_cutselect *sel,    */
/*      int savelp, int maxlocal, int update_tol, int silent,               */
/*      CCrandstate *rstate)                                                */
/*    CALLS CCtsp_cutting_loop multiple times, incrementing maxchunksize    */
/*      from 16 up to maxlocal, going up by 4 each time.                    */
/*                                                                          */
/*  int CCtsp_cutting_loop (CCtsp_lp *lp, CCtsp_cutselect *sel,             */
/*      int savelp, int silent, CCrandstate *rstate)                        */
/*    CALLS the cutting plane and pricing routines.                         */
/*     -sel should be set with the desired cut selection.                   */
/*     -savelp should be set to a nonzero value to write the lps to after   */
/*      rounds of cuts                                                      */
/*     -silent turns off most output if set to a nonzero value              */
/*    NOTES: It returns a 2 if the lp becomes infeasible                    */
/*                                                                          */
/*  int CCtsp_subtour_loop (CCtsp_lp *lp, int silent,                       */
/*      CCrandstate *rstate)                                                */
/*    CALLS the cutting and pricing to optimize over the subtour polytope.  */
/*                                                                          */
/*  int CCtsp_blossom_loop (CCtsp_lp *lp, int silent,                       */
/*      CCrandstate *rstate)                                                */
/*    CALLS the cutting and pricing to optimize over the blossom polytope.  */
/*                                                                          */
/*  int CCtsp_subtour_and_blossom_loop (CCtsp_lp *lp, int silent,           */
/*      CCrandstate *rstate)                                                */
/*    CALLS the cutting and princing to optimize over subtours and          */
/*     trivial blossoms.                                                    */
/*                                                                          */
/*  int CCtsp_pricing_loop (CCtsp_lp *lp, double *bnd, int silent,          */
/*      CCrandstate *rstate)                                                */
/*    ADDS negative reduced costs edges to lp and returns the current       */
/*     lowerbound.                                                          */
/*     -bnd can be NULL                                                     */
/*    NOTES: The LP must have full_edges_valid.                             */
/*                                                                          */
/*  int CCtsp_call_x_heuristic (CCtsp_lp *lp, double *val, int *outcyc,     */
/*      int silent, CCrandstate *rstate)                                    */
/*    CALLS the x-greedy LK heuristic with the current LP solution.         */
/*     -val returns the length of the tour.                                 */
/*     -outcyc will return the tour in node-node-node format if the         */
/*      length of the tour is less than lp->upperbound; the array should    */
/*      at least of length ncount (it can be NULL)                          */
/*                                                                          */
/*  int CCtsp_bb_cutting (char *probname, int probnum, int prob_newnum,     */
/*      int ncount, CCdatagroup *dat, int *ptour, double *upbound,          */
/*      CCtsp_lpcuts *pool, CCtsp_cutselect *sel, double *val,              */
/*      int *prune, int *foundtour, int *besttour, int level,               */
/*      int silent, CCrandstate *rstate)                                    */
/*    CALLS the cutting loop after reading the lp; writes the result to     */
/*     prob file prob_newnum; using exact price to verify pruned runs       */
/*     -upbound should be passed in as the current bound; if a better       */
/*      tour is found then upbound will be updated                          */
/*     -val returns the lp bound; it is CCtsp_LP_MAXDOUBLE if infeasible    */
/*     -prune is set to 1 if bbnode can be pruned                           */
/*     -foundtour is set to 1 if a better tour is found.                    */
/*     -besttour (if not NULL) will return a better tour if one is found.   */
/*     -level should be set to 0 for tentative cutting                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"
#include "bigguy.h"
#include "cut.h"
#include "pq.h"
#include "cuttree.h"
#include "consec1.h"
#include "necklace.h"
#include "localcut.h"
#if 0
#include "verify.h"
#endif

#define BIX_CHUNKS
#undef  DAVE_CHUNKS
#undef  OLD_CHUNKS
#undef  CUTTREE_DISPLAY
#undef  POLISHED_CHUNKS

static int
    call_add_cuts (CCtsp_lp *lp, CCtsp_lpcut_in **cuts, int *cut_added,
        int *xcount, int **xlist, double **x, double *val, int tighten,
        int *istour, int silent, CCrandstate *rstate),
    lp_value (CCtsp_lp *lp, double *val),
    lp_x (CCtsp_lp *lp, int *xcount, int **xlist, double **x),
    full_edge_check (CCtsp_lp *lp, int *nadded, int silent,
        CCrandstate *rstate),
    sparse_edge_check (CCtsp_lp *lp, CCtsp_edgegenerator *eg, int *nadded,
        double *bnd, int silent, CCrandstate *rstate),
    bb_cutting_work (CCtsp_lp **lp, char *probname, int probnum,
        int ncount, CCdatagroup *dat, int *ptour, double upbound,
        CCtsp_lpcuts *pool, CCtsp_cutselect *sel, double *val, int level,
        int silent, CCrandstate *rstate),
    grab_close_x (int ncount, int xcount, int *xlist, double *x,
        int *newcount, int **newlist, double **newx, double mult),
    no_tighten (int ncount, int xcount, int *xlist, double *x, int *test,
        double tol),
    CC_UNUSED grab_polished_x (CCtsp_lp *lp, double dust_val, int *newcount,
        int **newlist, double **newx),
    CC_UNUSED grab_polished2_x (CCtsp_lp *lp, double dust_val, int *newcount,
        int **newlist, double **newx);



void CCtsp_init_cutselect (CCtsp_cutselect *s)
{
    s->cutpool          = 1;  /* 1 */
    s->remotepool       = 0;
    s->remotehost       = (char *) NULL;
    s->remoteport       = 0;
    s->domboss          = 0;
    s->dombosshost      = (char *) NULL;
    s->connect          = 1;  /* 1 */
    s->segments         = 1;  /* 1 */
    s->blockcombs       = 1;  /* 1 */
    s->growcombs        = 0;  /* 0 */
    s->prclique         = 0;  /* 0 */
    s->exactsubtour     = 1;  /* 1 */
    s->exactblossom     = 0;  /* 0 */
    s->tighten_lp       = 1;  /* 1 */
    s->teething_lp      = 1;  /* 1 */
    s->tighten_pool     = 0;  /* 0 */   /* Slow when pool is large */
    s->teething_pool    = 0;  /* 0 */
    s->cliquetree_lp    = 0;  /* 0 */
    s->fastblossom      = 1;  /* 1 */
    s->ghfastblossom    = 1;  /* 1 */
    s->consecutiveones  = 1;  /* 1 */   /* Slow on large instances */
    s->necklace         = 0;  /* 0 */   /* Uses a big chunk of memory */
    s->usetighten       = 0;  /* 0 */
    s->extra_connect    = 0;  /* 0 */
    s->decker_lp        = 1;  /* 1 */
    s->decker_pool      = 0;  /* 0 */
    s->star_lp          = 0;  /* 0 */   /* Slow on very large instances */
    s->star_pool        = 0;  /* 0 */   /* Slow on most instances */
    s->handling_lp      = 1;  /* 1 */
    s->handling_pool    = 0;  /* 0 */
    s->maxchunksize     = 16; /* 16 */
    s->filecuts         = 0;  /* 0 */
    s->filecutname      = (char *) NULL;
    s->nexttol          = 0.0;
    s->roundtol         = 0.0;
    s->dominos          = 0;
    s->shrunk_dominos   = 0;
    s->fastcuts         = 0;           /* Keep this at 0 (changes tols) */
}

void CCtsp_cutselect_dominos (CCtsp_cutselect *s, int domsel)
{
    if (domsel == 1 || domsel == 3) s->dominos = 1;
    if (domsel == 2 || domsel == 3) s->shrunk_dominos = 1;
}

void CCtsp_cutselect_tighten (CCtsp_cutselect *s, int tighten)
{
    s->usetighten = tighten;
}

void CCtsp_cutselect_chunksize (CCtsp_cutselect *s, int chunksize)
{
    s->maxchunksize = chunksize;
}

void CCtsp_cutselect_filecuts (CCtsp_cutselect *s, char *fname)
{
    s->filecuts = 1;
    s->filecutname = fname;
}

void CCtsp_cutselect_remotepool (CCtsp_cutselect *s, char *cutbossname)
{
    s->remotepool = 1;
    s->remotehost = cutbossname;
    s->remoteport = CCtsp_CUT_PORT;
}

void CCtsp_cutselect_domboss (CCtsp_cutselect *s, char *dombossname)
{
    s->domboss = 1;
    s->dombosshost = dombossname;
}

void CCtsp_init_tentative_cutselect (CCtsp_cutselect *s)
{
    s->cutpool          = 1;
    s->remotepool       = 0;
    s->remotehost       = (char *) NULL;
    s->remoteport       = 0;
    s->domboss          = 0;
    s->dombosshost      = (char *) NULL;
    s->connect          = 1;
    s->segments         = 1;
    s->blockcombs       = 1;
    s->growcombs        = 0;
    s->prclique         = 0;
    s->exactsubtour     = 1;
    s->exactblossom     = 0;
    s->tighten_lp       = 1;
    s->teething_lp      = 1;
    s->tighten_pool     = 0;
    s->teething_pool    = 0;
    s->cliquetree_lp    = 0;
    s->fastblossom      = 1;
    s->ghfastblossom    = 1;
    s->consecutiveones  = 0;
    s->necklace         = 0; /* 1 */
    s->usetighten       = 1;
    s->extra_connect    = 0;
    s->decker_lp        = 1;
    s->decker_pool      = 0;
    s->star_lp          = 1;
    s->star_pool        = 0;
    s->handling_lp      = 1;
    s->handling_pool    = 0;
    s->maxchunksize     = 0;
    s->filecuts         = 0;
    s->filecutname      = (char *) NULL;
    s->nexttol          = 0.0;
    s->roundtol         = 0.0;
    s->dominos          = 0;
    s->shrunk_dominos   = 0;
    s->fastcuts         = 0;
}

void CCtsp_init_simple_cutselect (CCtsp_cutselect *s)
{
    s->cutpool          = 0;
    s->remotepool       = 0;
    s->remotehost       = (char *) NULL;
    s->remoteport       = 0;
    s->domboss          = 0;
    s->dombosshost      = (char *) NULL;
    s->connect          = 1;
    s->segments         = 1;
    s->blockcombs       = 0;
    s->growcombs        = 0;
    s->prclique         = 0;
    s->exactsubtour     = 0;
    s->exactblossom     = 0;
    s->tighten_lp       = 1;
    s->teething_lp      = 0;
    s->tighten_pool     = 0;
    s->teething_pool    = 0;
    s->cliquetree_lp    = 0;
    s->fastblossom      = 0;
    s->ghfastblossom    = 0;
    s->consecutiveones  = 0;
    s->necklace         = 0;
    s->usetighten       = 0;
    s->extra_connect    = 0;
    s->decker_lp        = 1;
    s->decker_pool      = 0;
    s->star_lp          = 0;
    s->star_pool        = 0;
    s->handling_lp      = 0;
    s->handling_pool    = 0;
    s->maxchunksize     = 0; 
    s->filecuts         = 0;
    s->filecutname      = (char *) NULL;
    s->nexttol          = 0.0;
    s->roundtol         = 0.0;
    s->dominos          = 0;
    s->shrunk_dominos   = 0;
    s->fastcuts         = 0;
}

void CCtsp_init_fast_cutselect (CCtsp_cutselect *s)
{
    s->cutpool          = 0;
    s->remotepool       = 0;
    s->remotehost       = (char *) NULL;
    s->remoteport       = 0;
    s->domboss          = 0;
    s->dombosshost      = (char *) NULL;
    s->connect          = 1;
    s->segments         = 1;
    s->blockcombs       = 1;
    s->growcombs        = 0;
    s->prclique         = 0;
    s->exactsubtour     = 1;
    s->exactblossom     = 0;
    s->tighten_lp       = 1;
    s->teething_lp      = 0;
    s->tighten_pool     = 0;
    s->teething_pool    = 0;
    s->cliquetree_lp    = 0;
    s->fastblossom      = 1;
    s->ghfastblossom    = 1;
    s->consecutiveones  = 0;
    s->necklace         = 0;
    s->usetighten       = 0;
    s->extra_connect    = 0;
    s->decker_lp        = 1;
    s->decker_pool      = 0;
    s->star_lp          = 0;
    s->star_pool        = 0;
    s->handling_lp      = 0;
    s->handling_pool    = 0;
    s->maxchunksize     = 0; 
    s->filecuts         = 0;
    s->filecutname      = (char *) NULL;
    s->nexttol          = 0.0;
    s->roundtol         = 0.0;
    s->dominos          = 0;
    s->shrunk_dominos   = 0;
    s->fastcuts         = 1;
}

int CCtsp_cutselect_set_tols (CCtsp_cutselect *s, CCtsp_lp *lp, int level,
        int silent)
{
    double ub, lb, beta;
    int rval;

    rval = CCtsp_get_lp_result (lp, &lb, &ub, (int *) NULL,
              (int **) NULL, (double **) NULL, (double **) NULL,
              (double **) NULL, (double **) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_get_lp_result failed\n"); return rval;
    }

#ifdef CCtsp_CUTS_DELTA
    if (ub - lb < 1.0) {
        beta = 1.0;
    } else {
        if (ub > 2*lb) beta = lb;
        else           beta = ub - lb;
    }
#else
    if (ub > 2*lb) beta = lb;
    else           beta = ub;
#endif

    if (level == -1) {
        s->nexttol  = 10.0 * CCtsp_TENTATIVE_CUTS_NEXT_TOL * beta;
        s->roundtol = 10.0 * CCtsp_TENTATIVE_CUTS_NEXT_ROUND * beta;
    } else if (level == 0) {
        s->nexttol  = CCtsp_TENTATIVE_CUTS_NEXT_TOL * beta;
        s->roundtol = CCtsp_TENTATIVE_CUTS_NEXT_ROUND * beta;
    } else {
        s->nexttol  = CCtsp_CUTS_NEXT_TOL * beta; 
        s->roundtol = CCtsp_CUTS_NEXT_ROUND * beta;
    }

    /* Simple Branching 
    s->nexttol =  2 * CCtsp_CUTS_NEXT_TOL * beta;
    s->roundtol = 2 * CCtsp_CUTS_NEXT_ROUND * beta;
    */


    if (!silent) {
        printf ("Setting tolerances: next cuts %.4f next round %.4f\n",
                s->nexttol, s->roundtol);
        fflush (stdout);
    }

    return 0;
}

int CCtsp_cutting_multiple_loop (CCtsp_lp *lp, CCtsp_cutselect *sel,
        int savelp, int maxlocal, int update_tol, int silent,
        CCrandstate *rstate)
{
    int rval = 0;
    int k;

    if (maxlocal < 16) {
        rval = CCtsp_cutting_loop (lp, sel, savelp, silent, rstate);
    } else {
        for (k = 16; k <= maxlocal; k += 4) {
            sel->maxchunksize = k;
            printf ("SETTING MAXCHUNKSIZE = %d\n", k); fflush (stdout);
            if (update_tol) {
                rval = CCtsp_cutselect_set_tols (sel, lp, 1, 0);
                if (rval) {
                    fprintf (stderr, "CCtsp_cutselect_set_tols failed\n");
                    /* Reset rval, since 2 has a special meaning */
                    rval = 1; goto CLEANUP; 
                }
            }
            rval = CCtsp_cutting_loop (lp, sel, savelp, silent, rstate);
            if (rval) {
                fprintf (stderr, "CCtsp_cutting_loop returned %d\n", rval);
                goto CLEANUP;
            }
        }
        if (maxlocal % 4 != 0) {
            sel->maxchunksize = maxlocal;
            printf ("SETTING MAXCHUNKSIZE = %d\n", maxlocal); fflush (stdout);
            if (update_tol) {
                rval = CCtsp_cutselect_set_tols (sel, lp, 1, 0);
                if (rval) {
                    fprintf (stderr, "CCtsp_cutselect_set_tols failed\n");
                    /* Reset rval, since 2 has a special meaning */
                    rval = 1; goto CLEANUP; 
                }
            }
            rval = CCtsp_cutting_loop (lp, sel, savelp, silent, rstate);
        }
    }

CLEANUP:

    return rval;
}


#define LOOP_FULL (25)      /* to force a full price after 25 inner loops */
#define CC_NO_NEAREST (50)  /* the initial sparse graph for pricing       */

int CCtsp_cutting_loop (CCtsp_lp *lp, CCtsp_cutselect *sel, int savelp,
        int silent, CCrandstate *rstate)
{
    int xcount, cutcount, cutcount_connect, cut_added, edge_added, closecount;
    int ttest = 0;
    int *xlist = (int *) NULL;
    int *closelist = (int *) NULL;
    int outside = 0;
    double newval, oldval, ub;
    double priceval;
    double *x = (double *) NULL;
    double *closex = (double *) NULL;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    CCtsp_edgegenerator eginside;
    int rval = 0;
    int loopcount = 0;
    int istour = 0;
    double maxviol;
    double save_nexttol  = sel->nexttol;
    double save_roundtol = sel->roundtol;
    double lcimprove = CCtsp_LP_MAXDOUBLE;
#if defined(BIX_CHUNKS) || defined(DAVE_CHUNKS) || defined(POLISHED_CHUNKS)
    double otherimprove = 0.0;
#endif
#if defined(POLISHED_CHUNKS)
    int i;
#endif
    double z;
    double dval;
    double szeit = CCutil_zeit ();
    CCchunk_localcut_timer lc_timer;

    CCutil_start_timer (&lp->stats.cutting_loop);
    CCchunk_init_localcut_timer (&lc_timer);
    
    eginside.ncount = 0;
    if (lp->fulladj) {
        rval = CCtsp_init_edgegenerator (&eginside, lp->graph.ncount, lp->dat,
                                   lp->fulladj, 0, silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_init_edgegenerator (sparse) failed\n");
            rval = 1; goto CLEANUP;
        }
    } else if (lp->dat) {
        rval = CCtsp_init_edgegenerator (&eginside, lp->graph.ncount, lp->dat,
                 (CCtsp_genadj *) NULL, CC_NO_NEAREST, silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_init_edgegenerator (sparse) failed\n");
            rval = 1; goto CLEANUP;
        }
    }

    rval = CCtsp_get_lp_result (lp, (double *) NULL, &ub, (int *) NULL,
              (int **) NULL, (double **) NULL, (double **) NULL,
              (double **) NULL, (double **) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_get_lp_result failed\n");
        rval = 1; goto CLEANUP;
    }

    do {
        loopcount = 0;
        cutcount_connect = 0;
        do {
            CCutil_start_timer (&lp->stats.cutting_inner_loop);

            rval = CCtsp_resparsify_lp (lp, silent);
            if (rval) {
                fprintf (stderr, "CCtsp_resparsify_lp failed\n");
                rval = 1; goto CLEANUP;
            }

            cut_added = 0;
            rval = lp_value (lp, &oldval);
            if (rval) {rval = 1; goto CLEANUP;}

            newval = oldval;

            rval = lp_x (lp, &xcount, &xlist, &x);
            if (rval) {rval = 1; goto CLEANUP;}

            rval = CCtsp_check_integral (lp, &dval, (int **) NULL, &istour,
                                         silent);
            if (rval) {
                fprintf (stderr, "CCtsp_check_integral failed\n"); goto CLEANUP;
            }
            if (istour) goto OUT_LOOP;

            if (sel->filecuts) {
                CCutil_start_timer (&lp->stats.cuts_filecut);
                rval = CCtsp_file_cuts (sel->filecutname, &cuts, &cutcount,
                           lp->graph.ncount, lp->perm);
                if (rval) {
                    fprintf (stderr, "CCtsp_file_cuts failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_filecut, 0);
                if (!silent) {
                    printf ("Found %2d file cuts in %.2f seconds\n",
                             cutcount, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_filecut_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_filecut_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->cutpool) {
                int qqq = 0;
                do {
                    CCutil_start_timer (&lp->stats.cuts_cutpool);
                    rval = CCtsp_search_cutpool (lp->pool, &cuts, &cutcount,
                          &maxviol, lp->graph.ncount, xcount, xlist, x, 0,
                          rstate);
                    if (rval) {
                        fprintf (stderr, "CCtsp_search_cutpool failed\n");
                        rval = 1; goto CLEANUP;
                    }
                    z = CCutil_stop_timer (&lp->stats.cuts_cutpool, 0);
                    if (cutcount)  {
                        if (!silent) {
                            printf ("Found %2d pool cuts (viol %.4f) in %.2f seconds\n",
                                     cutcount, maxviol, z);
                            fflush (stdout);
                        }
                        CCutil_start_timer (&lp->stats.cuts_cutpool_opt);
                        rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                        if (rval) {
                            fprintf (stderr, "call_add_cuts failed\n");
                            goto CLEANUP;
                        }
                        CCutil_stop_timer (&lp->stats.cuts_cutpool_opt, 0);
                        if (istour) goto OUT_LOOP;
                    }
                    qqq++;
                } while (cutcount && qqq < 1 /* 10 */);
            }

            if (sel->connect) {
                CCutil_start_timer (&lp->stats.cuts_connect);
                rval = CCtsp_connect_cuts (&cuts, &cutcount, lp->graph.ncount,
                                           xcount, xlist, x);
                if (rval) {
                    fprintf (stderr, "CCtsp_connect_cuts failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_connect, 0);
                if (!silent) {
                    printf ("Found %2d connect cuts in %.2f seconds\n",
                             cutcount, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_connect_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_connect_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->segments) {
                CCutil_start_timer (&lp->stats.cuts_segment);
                rval = CCtsp_segment_cuts (&cuts, &cutcount, lp->graph.ncount,
                                          xcount, xlist, x);
                if (rval) {
                    fprintf (stderr,  "CCtsp_segment_cuts failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_segment, 0);
                if (!silent) {
                    printf ("Found %2d segment cuts in %.2f seconds\n",
                             cutcount, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_segment_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_segment_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->fastblossom) {
                CCutil_start_timer (&lp->stats.cuts_fastblossom);
                rval = CCtsp_fastblossom (&cuts, &cutcount, lp->graph.ncount,
                                  xcount, xlist, x);
                if (rval) {
                    fprintf (stderr, "CCtsp_fastblossom failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_fastblossom, 0);
                if (!silent) {
                    printf ("Found %2d Fast Blossoms in %.2f seconds\n",
                             cutcount, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_fastblossom_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_fastblossom_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->ghfastblossom) {
                CCutil_start_timer (&lp->stats.cuts_ghfastblossom);
                rval = CCtsp_ghfastblossom (&cuts, &cutcount, lp->graph.ncount,
                                  xcount, xlist, x);
                if (rval) {
                    fprintf (stderr, "CCtsp_ghfastblossom failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_ghfastblossom, 0);
                if (!silent) {
                    printf ("Found %2d Groetschel-Holland Blossoms in %.2f seconds\n",
                             cutcount, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_ghfastblossom_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_ghfastblossom_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->remotepool) {
                int qqq = 0;
                do {
                    CCutil_start_timer (&lp->stats.cuts_remotepool);
                    rval = CCtsp_search_remotepool (sel->remotehost,
                            sel->remoteport, &cuts, &cutcount, &maxviol,
                            lp->graph.ncount, xcount, xlist, x);
                    if (rval) {
                        fprintf (stderr, "CCtsp_search_remotepool failed\n");
                        rval = 1; goto CLEANUP;
                    }
                    z = CCutil_stop_timer (&lp->stats.cuts_remotepool, 0);
                    if (cutcount)  {
                        if (!silent) {
                            printf ("%d remote pool cuts found (viol %.3f) in %.2f seconds\n",
                                     cutcount, maxviol, z);
                            fflush (stdout);
                        }
                        CCutil_start_timer (&lp->stats.cuts_remotepool_opt);
                        rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                      &xlist, &x, &newval, sel->usetighten,
                                      &istour, silent, rstate);
                        if (rval) {
                            fprintf (stderr, "call_add_cuts failed\n");
                            goto CLEANUP;
                        }
                        CCutil_stop_timer (&lp->stats.cuts_remotepool_opt, 0);
                        if (istour) goto OUT_LOOP;
                    }
                    qqq++;
                } while (cutcount && qqq < 1 /* 1 */);
            }

            if (sel->blockcombs) {
                CCutil_start_timer (&lp->stats.cuts_blockcomb);
                rval = CCtsp_block_combs (&cuts, &cutcount, lp->graph.ncount,
                                            xcount, xlist, x, silent);
                if (rval) {
                    fprintf (stderr, "CCtsp_block_combs failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_blockcomb, 0);
                if (!silent) { 
                    printf ("Found %2d block_combs in %.2f seconds\n",
                             cutcount, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_blockcomb_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_blockcomb_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->growcombs) {
                CCutil_start_timer (&lp->stats.cuts_growcomb);
                rval = CCtsp_edge_comb_grower (&cuts, &cutcount,
                        lp->graph.ncount, xcount, xlist, x,
                        &lp->stats.extra_tighten_stats);
                if (rval) {
                    fprintf (stderr, "CCtsp_block_combs failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_growcomb, 0);
                if (!silent) {
                    printf ("Found %2d grown combs in %.2f seconds\n",
                             cutcount, z);
                             fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_growcomb_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_growcomb_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->prclique) {
                rval = CCtsp_pr_cliquetree (&cuts, &cutcount,
                               lp->graph.ncount, xcount, xlist, x,
                               &lp->stats.extra_tighten_stats);
                if (rval) {
                    fprintf (stderr, "CCtsp_pr_cliquetree failed\n");
                    rval = 1; goto CLEANUP;
                }
                if (!silent) {
                    printf ("Found %2d PR cliquetrees\n", cutcount);
                    fflush (stdout);
                }
                if (cutcount) {
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->exactsubtour) {
                CCutil_start_timer (&lp->stats.cuts_exactsubtour);
                rval = CCtsp_exact_subtours (&cuts, &cutcount,
                            lp->graph.ncount, xcount, xlist, x);
                if (rval) {
                    fprintf (stderr, "CCtsp_exact_subtours failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_exactsubtour, 0);
                if (!silent) {
                    printf ("Found %2d exact subtours in %.2f seconds\n",
                            cutcount, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_exactsubtour_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_exactsubtour_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

#ifdef CCtsp_USE_DOMINO_CUTS
            if (sel->shrunk_dominos) {
                do {
                    cut_added = 0;
                    CCutil_start_timer (&lp->stats.cuts_exactsubtour);
                    rval = CCtsp_exact_subtours (&cuts, &cutcount,
                                lp->graph.ncount, xcount, xlist, x);
                    if (rval) {
                        fprintf (stderr, "CCtsp_exact_subtours failed\n");
                        rval = 1; goto CLEANUP;
                    }
                    z = CCutil_stop_timer (&lp->stats.cuts_exactsubtour, 0);
                    if (!silent) {
                        printf ("Found %2d exact subtours in %.2f seconds\n",
                                cutcount, z);
                        fflush (stdout);
                    }
                    if (cutcount) {
                        CCutil_start_timer (&lp->stats.cuts_exactsubtour_opt);
                        rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                        if (rval) {
                            fprintf (stderr, "call_add_cuts failed\n");
                            goto CLEANUP;
                        }
                        CCutil_stop_timer (&lp->stats.cuts_exactsubtour_opt, 0);
                        if (istour) goto OUT_LOOP;
                    }
                } while (cut_added > 0);

                CCutil_start_timer (&lp->stats.cuts_exactsubtour);
                rval = CCtsp_shrink_domino (&cuts, &cutcount, lp->graph.ncount,
                        xcount, xlist, x, 1, 5, rstate, sel->dombosshost);
                CCcheck_rval (rval, "CCtsp_new_domino failed");

                z = CCutil_stop_timer (&lp->stats.cuts_exactsubtour, 0);
                if (0 || !silent) {
                    printf ("Found %2d domino cuts in %.2f seconds\n",
                            cutcount, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_exactsubtour_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, 0,
                                          &istour, 0 /* silent */, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_exactsubtour_opt, 0);
                    if (istour) goto OUT_LOOP;
                    printf ("Added %d shrunk domino cuts\n", cut_added);
                    fflush (stdout);
                }
            }

            if (sel->dominos) {
                do {
                    cut_added = 0;
                    CCutil_start_timer (&lp->stats.cuts_exactsubtour);
                    rval = CCtsp_exact_subtours (&cuts, &cutcount,
                                lp->graph.ncount, xcount, xlist, x);
                    if (rval) {
                        fprintf (stderr, "CCtsp_exact_subtours failed\n");
                        rval = 1; goto CLEANUP;
                    }
                    z = CCutil_stop_timer (&lp->stats.cuts_exactsubtour, 0);
                    if (!silent) {
                        printf ("Found %2d exact subtours in %.2f seconds\n",
                                cutcount, z);
                        fflush (stdout);
                    }
                    if (cutcount) {
                        CCutil_start_timer (&lp->stats.cuts_exactsubtour_opt);
                        rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                        if (rval) {
                            fprintf (stderr, "call_add_cuts failed\n");
                            goto CLEANUP;
                        }
                        CCutil_stop_timer (&lp->stats.cuts_exactsubtour_opt, 0);
                        if (istour) goto OUT_LOOP;
                    }
                } while (cut_added > 0);

                CCutil_start_timer (&lp->stats.cuts_exactsubtour);
                rval = CCtsp_new_domino (&cuts, &cutcount, lp->graph.ncount,
                        xcount, xlist, x, sel->dombosshost);
                CCcheck_rval (rval, "CCtsp_new_domino failed");

                z = CCutil_stop_timer (&lp->stats.cuts_exactsubtour, 0);
                if (0 || !silent) {
                    printf ("Found %2d domino cuts in %.2f seconds\n",
                            cutcount, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_exactsubtour_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, 0,
                                          &istour, 0 /* silent */, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_exactsubtour_opt, 0);
                    if (istour) goto OUT_LOOP;
                    printf ("Added %d domino cuts\n", cut_added);
                    fflush (stdout);
                }
            }
#endif

            if (sel->exactblossom) {
                CCutil_start_timer (&lp->stats.cuts_exactblossom);
                rval = CCtsp_exactblossom (&cuts, &cutcount, lp->graph.ncount,
                                  xcount, xlist, x, rstate);
                if (rval) {
                    fprintf (stderr, "CCtsp_exactblossom failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_exactblossom, 0);
                if (!silent) {
                    printf ("Found %2d Exact Blossoms in %.2f seconds\n",
                            cutcount, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_exactblossom_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_exactblossom_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->tighten_lp) {
                CCutil_start_timer (&lp->stats.cuts_tighten_lp);
                rval = CCtsp_tighten_lp (&lp->cuts,
                        &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                        lp->graph.ncount, xcount, xlist, x, 0.5, 500, &maxviol,
                        rstate);
                if (rval) {
                    fprintf (stderr, "CCtsp_tighten_lp failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_tighten_lp, 0);
                if (!silent) {
                    printf ("Found %2d tighten_lp cuts (viol %.4f) in %.2f seconds\n",
                            cutcount, maxviol, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_tighten_lp_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_tighten_lp_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
                CCutil_start_timer (&lp->stats.cuts_tighten_lp_close);
                rval = grab_close_x (lp->graph.ncount, xcount, xlist, x,
                        &closecount, &closelist, &closex, 0.5);
                if (rval) {
                    fprintf (stderr, "grab_close_x failed\n");
                    goto CLEANUP;
                }
                rval = CCtsp_tighten_lp (&lp->cuts,
                        &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                        lp->graph.ncount, closecount, closelist, closex,
                        0.5, 500, &maxviol, rstate);
                if (rval) {
                    fprintf (stderr, "CCtsp_tighten_lp failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_tighten_lp_close, 0);
                if (!silent) {
                    printf ("Found %2d CLOSE tighten_lp cuts (viol %.4f) in %.2f seconds\n",
                            cutcount, maxviol, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_tighten_lp_close_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_tighten_lp_close_opt,
                                       0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->decker_lp) {
                CCutil_start_timer (&lp->stats.cuts_decker_lp);
                rval = CCtsp_double_decker_lp (&lp->cuts,
                        &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                        lp->graph.ncount, xcount, xlist, x, 10.0, 500,
                        &maxviol, rstate);
                if (rval) {
                    fprintf (stderr, "CCtsp_double_decker_lp failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_decker_lp, 0);
                if (!silent) {
                    printf ("Found %2d double deckers (viol %.4f) in %.2f seconds\n",
                            cutcount, maxviol, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_decker_lp_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_decker_lp_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
                CCutil_start_timer (&lp->stats.cuts_decker_lp_close);
                rval = grab_close_x (lp->graph.ncount, xcount, xlist, x,
                        &closecount, &closelist, &closex, 0.5);
                if (rval) {
                    fprintf (stderr, "grab_close_x failed\n");
                    goto CLEANUP;
                }
                rval = CCtsp_double_decker_lp (&lp->cuts,
                        &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                        lp->graph.ncount, closecount, closelist, closex,
                        10.0, 500, &maxviol, rstate);
                if (rval) {
                    fprintf (stderr, "CCtsp_double_decker_lp failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_decker_lp_close, 0);
                if (!silent) {
                    printf ("Found %2d CLOSE double deckers (viol %.4f) in %.2f seconds\n",
                            cutcount, maxviol, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_decker_lp_close_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_decker_lp_close_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }


            if (sel->star_lp) {
                CCutil_start_timer (&lp->stats.cuts_star_lp);
                rval = CCtsp_star_lp (&lp->cuts,
                        &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                        lp->graph.ncount, xcount, xlist, x, 10.0, 500,
                        &maxviol, rstate);
                if (rval) {
                    fprintf (stderr, "CCtsp_star_lp failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_star_lp, 0);
                if (!silent) {
                    printf ("Found %2d Star Inequalities (viol %.4f) in %.2f seconds\n",
                            cutcount, maxviol, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_star_lp_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                &xlist, &x, &newval, sel->usetighten, &istour,
                                silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_star_lp_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->handling_lp) {
                rval = no_tighten (lp->graph.ncount, xcount, xlist, x, &ttest,
                                   0.20);
                if (rval) {
                    fprintf (stderr, "no_tighten failed\n");  goto CLEANUP;
                }
            }

            if (sel->handling_lp && ttest == 0) {
                CCutil_start_timer (&lp->stats.cuts_handling_lp);
                rval = CCtsp_handling_lp (&lp->cuts,
                        &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                        lp->graph.ncount, xcount, xlist, x, 10.0, 500,
                        &maxviol, rstate);
                if (rval) {
                    fprintf (stderr, "CCtsp_handling_lp failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_handling_lp, 0);
                if (!silent) {
                    printf ("Found %2d Handling Inequalities (viol %.4f) in %.2f seconds\n",
                            cutcount, maxviol, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_handling_lp_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                &xlist, &x, &newval, sel->usetighten, &istour,
                                silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_handling_lp_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->cliquetree_lp) {
                CCutil_start_timer (&lp->stats.cuts_cliquetree_lp);
                rval = CCtsp_cliquetree_lp (&lp->cuts,
                        &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                        lp->graph.ncount, xcount, xlist, x, 10.0, 500,
                        &maxviol, rstate);
                if (rval) {
                    fprintf (stderr, "CCtsp_cliqutree_lp failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_cliquetree_lp, 0);
                if (!silent) {
                    printf ("Found %2d comb cliquetrees (viol %.4f) in %.2f seconds\n",
                            cutcount, maxviol, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_cliquetree_lp_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_cliquetree_lp_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->teething_lp) {
                rval = no_tighten (lp->graph.ncount, xcount, xlist, x, &ttest,
                                   0.20);
                if (rval) {
                    fprintf (stderr, "no_tighten failed\n");  goto CLEANUP;
                }
            }

            if (sel->teething_lp && ttest == 0) {
                CCutil_start_timer (&lp->stats.cuts_teething_lp);
                rval = CCtsp_teething_lp (&lp->cuts,
                        &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                        lp->graph.ncount, xcount, xlist, x, 10.0, 500,
                        &maxviol, rstate);
                if (rval) {
                    fprintf (stderr, "CCtsp_teething_lp failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_teething_lp, 0);
                if (!silent) {
                    printf ("Found %2d teethed combs (viol %.4f) in %.2f seconds\n",
                            cutcount, maxviol, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_teething_lp_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_teething_lp_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->tighten_pool) {
                rval = no_tighten (lp->graph.ncount, xcount, xlist, x, &ttest,
                                   0.2);
                if (rval) {
                    fprintf (stderr, "no_tighten failed\n");  goto CLEANUP;
                }
            }

            if (sel->tighten_pool && newval < oldval + sel->nexttol &&
                ttest == 0) {
                CCutil_start_timer (&lp->stats.cuts_tighten_pool);
                rval = CCtsp_tighten_lp (lp->pool,
                        &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                        lp->graph.ncount, xcount, xlist, x, 0.1, 1000,
                        &maxviol, rstate);
                if (rval) {
                    fprintf (stderr, "CCtsp_tighten_lp failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_tighten_pool, 0);
                if (!silent) {
                    printf ("Found %2d tighten pool cuts (viol %.4f) in %.2f seconds\n",
                            cutcount, maxviol, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_tighten_pool_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_tighten_pool_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->decker_pool && newval < oldval + sel->nexttol) {
                CCutil_start_timer (&lp->stats.cuts_decker_pool);
                rval = CCtsp_double_decker_lp (lp->pool,
                        &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                        lp->graph.ncount, xcount, xlist, x, 2.0, 1000,
                        &maxviol, rstate);
                if (rval) {
                    fprintf (stderr, "CCtsp_double_decker_lp failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_decker_pool, 0);
                if (!silent) {
                    printf ("Found %2d pool double deckers (viol %.4f) in %.2f seconds\n",
                            cutcount, maxviol, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_decker_pool_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_decker_pool_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->star_pool && newval < oldval + sel->nexttol) {
                CCutil_start_timer (&lp->stats.cuts_star_pool);
                rval = CCtsp_star_lp (lp->pool,
                        &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                        lp->graph.ncount, xcount, xlist, x, 2.0, 1000,
                        &maxviol, rstate);
                if (rval) {
                    fprintf (stderr, "CCtsp_star_lp failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_star_pool, 0);
                if (!silent) {
                    printf ("Found %2d Pool Star Inequalities (viol %.4f) in %.2f seconds\n",
                                cutcount, maxviol, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_star_pool_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                &xlist, &x, &newval, sel->usetighten, &istour,
                                silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_star_pool_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->handling_pool && newval < oldval + sel->nexttol) {
                CCutil_start_timer (&lp->stats.cuts_handling_pool);
                rval = CCtsp_handling_lp (lp->pool,
                        &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                        lp->graph.ncount, xcount, xlist, x, 2.0, 1000,
                        &maxviol, rstate);
                if (rval) {
                    fprintf (stderr, "CCtsp_handling_lp failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_handling_pool, 0);
                if (!silent) {
                    printf ("Found %2d Pool Handling Inequalities (viol %.4f) in %.2f seconds\n",
                             cutcount, maxviol, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_handling_pool_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                &xlist, &x, &newval, sel->usetighten, &istour,
                                silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_handling_pool_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->teething_pool && newval < oldval + sel->nexttol) {
                CCutil_start_timer (&lp->stats.cuts_teething_pool);
                rval = CCtsp_teething_lp (lp->pool,
                        &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                        lp->graph.ncount, xcount, xlist, x, 0.5, 1000,
                        &maxviol, rstate);
                if (rval) {
                    fprintf (stderr, "CCtsp_teething_lp failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_teething_pool, 0);
                if (!silent) {
                    printf ("Found %2d pool teething combs (viol %.4f) in %.2f seconds\n",
                             cutcount, maxviol, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_teething_pool_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_teething_pool_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->consecutiveones && newval < oldval + sel->nexttol) {
                CCutil_start_timer (&lp->stats.cuts_consecutiveones);
                rval = CCpq_cuttree_improve_quick (&lp->tightcuts, lp->pool,
                            xcount, xlist, x);
                if (rval) {
                    fprintf (stderr, "CCpq_cuttree_improve_quick failed\n");
                    rval = 1; goto CLEANUP;
                }

                rval = CCpq_consecutiveones (&cuts, &cutcount, &lp->tightcuts,
                            lp->pool, xcount, xlist, x);
                if (rval) {
                    fprintf (stderr, "CCpq_consecutiveones failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_consecutiveones, 0);
                if (!silent) {
                    printf ("Found %2d consecutiveones cuts in %.2f seconds\n",
                            cutcount, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_consecutiveones_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_consecutiveones_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->necklace && newval < oldval + sel->nexttol) {
                CCutil_start_timer (&lp->stats.cuts_necklace);
                rval = CCpq_cuttree_improve_quick (&lp->tightcuts, lp->pool,
                            xcount, xlist, x);
                if (rval) {
                    fprintf (stderr, "CCpq_cuttree_improve_quick failed\n");
                    rval = 1; goto CLEANUP;
                }

                rval = CCpq_necklaces (&cuts, &cutcount, &lp->tightcuts,
                            xcount, xlist, x, rstate);
                if (rval) {
                    fprintf (stderr, "CCpq_necklaces failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_necklace, 0);

                if (!silent) {
                    printf ("Found %2d necklace cuts in %.2f seconds\n",
                            cutcount, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_necklace_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_necklace_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

#ifdef  DAVE_CHUNKS
            otherimprove = newval - oldval;
            if (sel->maxchunksize > 0 && newval < oldval + sel->nexttol &&
                otherimprove <  0.5 * lcimprove) {
                CCchunk_flag flags;

                flags.uncivilized = 0;
                flags.noshrink = 0;
                flags.nolift = 0;
                flags.maxchunksize = sel->maxchunksize;

                CCutil_start_timer (&lp->stats.cuts_localcut);
                rval = CCchunk_localcuts (&cuts, &cutcount, lp->graph.ncount,
                              xcount, xlist, x, 0.0, flags, &lc_timer, silent,
                              rstate);
                if (rval) {
                    fprintf (stderr, "LocalCuts failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_localcut, 0);
                if (cutcount) {
                    if (!silent) {
                        printf ("Found %2d LocalCuts in %.2f seconds\n",
                                 cutcount, z);
                        fflush (stdout);
                    }
                    CCutil_start_timer (&lp->stats.cuts_localcut_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_localcut_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }
#endif /* DAVE_CHUNKS */

#ifdef  BIX_CHUNKS
            otherimprove = newval - oldval;
            if (sel->maxchunksize > 0 && newval < oldval + sel->nexttol &&
                otherimprove <  0.5 * lcimprove) {
                int  maxchunksize, firstsize;
                CCchunk_flag flags;

                flags.dummy = 0;
                flags.permute = 0;
                flags.weighted = 0;
                flags.spheres = 0;
                flags.uncivilized = 0;
                flags.noshrink = 0;
                flags.nolift = 0;

                if (sel->maxchunksize < 8)  firstsize = sel->maxchunksize;
                else                        firstsize = 8;
                for (maxchunksize = firstsize;
                     maxchunksize <= sel->maxchunksize;
                     maxchunksize++) {
                    flags.maxchunksize = maxchunksize;
                    flags.spheresize   = maxchunksize - 2;

                    CCutil_start_timer (&lp->stats.cuts_localcut);
                    rval = CCchunk_localcuts (&cuts, &cutcount,
                         lp->graph.ncount, xcount, xlist, x, 0.0, flags,
                         &lc_timer, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "LocalCuts failed\n");
                        rval = 1; goto CLEANUP;
                    }
                    z = CCutil_stop_timer (&lp->stats.cuts_localcut, 0);
                    if (!silent) {
                        printf ("Found %2d LocalCuts in %.2f seconds\n",
                                cutcount, z);
                        fflush (stdout);
                    }
                    if (cutcount) {
                        CCutil_start_timer (&lp->stats.cuts_localcut_opt);
                        rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                         &xlist, &x, &newval, sel->usetighten,
                                         &istour, silent, rstate);
                        if (rval) {
                            fprintf (stderr, "call_add_cuts failed\n");
                            goto CLEANUP;
                        }
                        CCutil_stop_timer (&lp->stats.cuts_localcut_opt, 0);
                        if (istour) goto OUT_LOOP;
                    }
                    if (newval >= oldval + sel->nexttol)  break;
                }
            }

            if (sel->maxchunksize > 0 && newval < oldval + sel->nexttol) {
                int  maxchunksize, firstsize;
                CCchunk_flag flags;
                double beforeval = newval;

                flags.dummy = 0;
                flags.permute = 0;
                flags.weighted = 0;
                flags.spheres = 1;
                flags.uncivilized = 0;
                flags.noshrink = 0;
                flags.nolift = 0;

                if (sel->maxchunksize < 8)  firstsize = sel->maxchunksize;
                else                        firstsize = 8;
                for (maxchunksize = firstsize;
                     maxchunksize <= sel->maxchunksize;
                     maxchunksize++) {
                    flags.maxchunksize = maxchunksize;
                    flags.spheresize = maxchunksize - 2;

                    CCutil_start_timer (&lp->stats.cuts_localcut);
                    rval = CCchunk_localcuts (&cuts, &cutcount,
                             lp->graph.ncount, xcount, xlist, x, 0.0, flags,
                             &lc_timer, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "LocalCuts failed\n");
                        rval = 1; goto CLEANUP;
                    }
                    z = CCutil_stop_timer (&lp->stats.cuts_localcut, 0);
                    if (!silent) {
                        printf ("Found %2d LocalCuts in %.2f seconds\n",
                                 cutcount, z);
                        fflush (stdout);
                    }
                    if (cutcount) {
                        CCutil_start_timer (&lp->stats.cuts_localcut_opt);
                        rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                         &xlist, &x, &newval, sel->usetighten,
                                         &istour, silent, rstate);
                        if (rval) {
                            fprintf (stderr, "call_add_cuts failed\n");
                            goto CLEANUP;
                        }
                        CCutil_stop_timer (&lp->stats.cuts_localcut_opt, 0);
                        if (istour) goto OUT_LOOP;
                    }
                    if (newval >= oldval + sel->nexttol)  break;
                }
                lcimprove = newval - beforeval;
            }
#if 0
            if (sel->maxchunksize > 0 && newval < oldval + sel->nexttol) {
                int i;
                for (i = 0; i < 3; i++) {
                    int  maxchunksize, firstsize;
                    CCchunk_flag flags;
                    double beforeval = newval;

                    flags.dummy = 0;
                    flags.permute = 0;
                    flags.weighted = 0;
                    flags.spheres = 0;
                    flags.uncivilized = 0;
                    flags.noshrink = 0;
                    flags.nolift = 0;

                    if (i == 0)      flags.dummy = 1;
                    else if (i == 1) flags.permute = 1;
                    else if (i == 2) flags.weighted = 1;
                   

                    if (sel->maxchunksize < 8)  firstsize = sel->maxchunksize;
                    else                        firstsize = 8;
                    for (maxchunksize = firstsize;
                         maxchunksize <= sel->maxchunksize;
                         maxchunksize++) {
                        flags.maxchunksize = maxchunksize;
                        flags.spheresize = maxchunksize - 2;

                        CCutil_start_timer (&lp->stats.cuts_localcut);
                        rval = CCchunk_localcuts (&cuts, &cutcount,
                                  lp->graph.ncount, xcount, xlist, x, 0.0,
                                  flags, &lc_timer, silent, rstate);
                        if (rval) {
                            fprintf (stderr, "LocalCuts failed\n");
                            rval = 1; goto CLEANUP;
                        }
                        z = CCutil_stop_timer (&lp->stats.cuts_localcut, 0);
                        if (!silent) {
                            printf ("Found %2d LocalCuts in %.2f seconds\n",
                                     cutcount, z);
                            fflush (stdout);
                        }
                        if (cutcount) {
                            CCutil_start_timer (&lp->stats.cuts_localcut_opt);
                            rval = call_add_cuts (lp, &cuts, &cut_added,
                                    &xcount, &xlist, &x, &newval,
                                    sel->usetighten, &istour, silent,
                                    rstate);
                            if (rval) {
                                fprintf (stderr, "call_add_cuts failed\n");
                                goto CLEANUP;
                            }
                            CCutil_stop_timer (&lp->stats.cuts_localcut_opt, 0);
                            if (istour) goto OUT_LOOP;
                        }
                        if (newval >= oldval + sel->nexttol)  break;
                    }
                    lcimprove = newval - beforeval;
                }
            }
#endif
#endif /* BIX_CHUNKS */
#ifdef OLD_CHUNKS
            if (sel->maxchunksize > 0 && newval < oldval + sel->nexttol) {
                CCchunk_flag flags;

                flags.dummy = 0;
                flags.permute = 0;
                flags.weighted = 0;
                flags.spheres = 0;
                flags.uncivilized = 0;
                flags.noshrink = 0;
                flags.nolift = 0;

                flags.maxchunksize = sel->maxchunksize;
                flags.spheresize = flags.maxchunksize - 2;
                CCutil_start_timer (&lp->stats.cuts_localcut);
                rval = CCchunk_localcuts (&cuts, &cutcount, lp->graph.ncount,
                           xcount, xlist, x, 0.0, flags, &lc_timer, silent,
                           rstate);
                if (rval) {
                    fprintf (stderr, "LocalCuts failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_localcut, 0);
                if (!silent) {
                    printf ("Found %2d LocalCuts in %.2f seconds\n",
                             cutcount, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_localcut_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                        &xlist, &x, &newval, sel->usetighten,
                                        &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_localcut_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }

            if (sel->maxchunksize > 0 && newval < oldval + sel->nexttol) {
                CCchunk_flag flags;

                flags.dummy = 0;
                flags.permute = 0;
                flags.weighted = 0;
                flags.spheres = 1;
                flags.uncivilized = 0;
                flags.noshrink = 0;
                flags.nolift = 0;

                flags.maxchunksize = sel->maxchunksize;
                flags.spheresize = flags.maxchunksize - 2;
                CCutil_start_timer (&lp->stats.cuts_localcut);
                rval = CCchunk_localcuts (&cuts, &cutcount, lp->graph.ncount,
                         xcount, xlist, x, 0.0, flags, &lc_timer, silent,
                         rstate);
                if (rval) {
                    fprintf (stderr, "LocalCuts failed\n");
                    rval = 1; goto CLEANUP;
                }
                z = CCutil_stop_timer (&lp->stats.cuts_localcut, 0);
                if (!silent) {
                    printf ("Found %2d LocalCuts in %.2f seconds\n",
                             cutcount, z);
                    fflush (stdout);
                }
                if (cutcount) {
                    CCutil_start_timer (&lp->stats.cuts_localcut_opt);
                    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                          &xlist, &x, &newval, sel->usetighten,
                                          &istour, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "call_add_cuts failed\n");
                        goto CLEANUP;
                    }
                    CCutil_stop_timer (&lp->stats.cuts_localcut_opt, 0);
                    if (istour) goto OUT_LOOP;
                }
            }
#endif /* OLD_CHUNKS */

#ifdef  POLISHED_CHUNKS
/* polished localcuts */
            rval = grab_polished2_x (lp, 0.001,
                                    &closecount, &closelist, &closex);
            if (rval) {
                fprintf (stderr, "grab_polished_x failed\n");
                goto CLEANUP;
            }
            
            {
                static int zzy = 1;
                char nbuf[1024];
                FILE *fout;
                int zzz;
                sprintf (nbuf, "x.polished.%d", zzy);
                fout =  fopen (nbuf,"w");
                fprintf (fout, "%d %d\n", lp->graph.ncount, closecount);
                for (zzz = 0; zzz < closecount; zzz++) {
                    fprintf (fout, "%d %d %.6f\n", closelist[2*zzz],
                             closelist[2*zzz+1], closex[zzz]);
                }
                fclose (fout);
                sprintf (nbuf, "x.dusty.%d", zzy);
                fout = fopen (nbuf, "w");
                fprintf (fout, "%d %d\n", lp->graph.ncount, xcount);
                for (zzz = 0; zzz < xcount; zzz++) {
                    fprintf (fout, "%d %d %.6f\n", xlist[2*zzz],
                             xlist[2*zzz+1], x[zzz]);
                }
                fclose (fout);
                zzy++;
            }
                
            otherimprove = newval - oldval;
            if (sel->maxchunksize > 0 && newval < oldval + sel->nexttol &&
                otherimprove <  0.5 * lcimprove) {
                int  maxchunksize, firstsize;
                CCchunk_flag flags;

                flags.dummy = 0;
                flags.permute = 0;
                flags.weighted = 0;
                flags.spheres = 0;
                flags.uncivilized = 0;
                flags.noshrink = 0;
                flags.nolift = 0;

                if (sel->maxchunksize < 8)  firstsize = sel->maxchunksize;
                else                        firstsize = 8;
                for (maxchunksize = firstsize;
                     maxchunksize <= sel->maxchunksize;
                     maxchunksize++) {
                    flags.maxchunksize = maxchunksize;
                    flags.spheresize   = maxchunksize - 2;

                    CCutil_start_timer (&lp->stats.cuts_localcut);
                    rval = CCchunk_localcuts (&cuts, &cutcount,
                             lp->graph.ncount, closecount, closelist, closex,
                             0.0, flags, &lc_timer, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "LocalCuts failed\n");
                        rval = 1; goto CLEANUP;
                    }
                    z = CCutil_stop_timer (&lp->stats.cuts_localcut, 0);
                    if (!silent) {
                        printf ("Found %2d POLISHED LocalCuts in %.2f seconds\n",
                                 cutcount, z);
                        fflush (stdout);
                    }
                    if (cutcount) {
                        CCutil_start_timer (&lp->stats.cuts_localcut_opt);
                        rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                         &xlist, &x, &newval, sel->usetighten,
                                         &istour, silent, rstate);
                        if (rval) {
                            fprintf (stderr, "call_add_cuts failed\n");
                            goto CLEANUP;
                        }
                        CCutil_stop_timer (&lp->stats.cuts_localcut_opt, 0);
                        if (istour) goto OUT_LOOP;
                    }
                    if (newval >= oldval + sel->nexttol)  break;
                }
            }

            if (sel->maxchunksize > 0 && newval < oldval + sel->nexttol) {
                int  maxchunksize, firstsize;
                CCchunk_flag flags;
                double beforeval = newval;

                flags.dummy = 0;
                flags.permute = 0;
                flags.weighted = 0;
                flags.spheres = 1;
                flags.uncivilized = 0;
                flags.noshrink = 0;
                flags.nolift = 0;

                if (sel->maxchunksize < 8)  firstsize = sel->maxchunksize;
                else                        firstsize = 8;
                for (maxchunksize = firstsize;
                     maxchunksize <= sel->maxchunksize;
                     maxchunksize++) {
                    flags.maxchunksize = maxchunksize;
                    flags.spheresize = maxchunksize - 2;

                    CCutil_start_timer (&lp->stats.cuts_localcut);
                    rval = CCchunk_localcuts (&cuts, &cutcount,
                             lp->graph.ncount, closecount, closelist, closex,
                             0.0, flags, &lc_timer, silent, rstate);
                    if (rval) {
                        fprintf (stderr, "LocalCuts failed\n");
                        rval = 1; goto CLEANUP;
                    }
                    z = CCutil_stop_timer (&lp->stats.cuts_localcut, 0);
                    if (!silent) {
                        printf ("Found %2d POLISHED LocalCuts in %.2f seconds\n",
                                 cutcount, z);
                        fflush (stdout);
                    }
                    if (cutcount) {
                        CCutil_start_timer (&lp->stats.cuts_localcut_opt);
                        rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                         &xlist, &x, &newval, sel->usetighten,
                                         &istour, silent, rstate);
                        if (rval) {
                            fprintf (stderr, "call_add_cuts failed\n");
                            goto CLEANUP;
                        }
                        CCutil_stop_timer (&lp->stats.cuts_localcut_opt, 0);
                        if (istour) goto OUT_LOOP;
                    }
                    if (newval >= oldval + sel->nexttol)  break;
                }
                lcimprove = newval - beforeval;
            }
            if (sel->maxchunksize > 0 && newval < oldval + sel->nexttol) {
                for (i = 0; i < 3; i++) {
                    int  maxchunksize, firstsize;
                    CCchunk_flag flags;
                    double beforeval = newval;

                    flags.dummy = 0;
                    flags.permute = 0;
                    flags.weighted = 0;
                    flags.spheres = 0;
                    flags.uncivilized = 0;
                    flags.noshrink = 0;
                    flags.nolift = 0;

                    if (i == 0)      flags.dummy = 1;
                    else if (i == 1) flags.permute = 1;
                    else if (i == 2) flags.weighted = 1;
                   

                    if (sel->maxchunksize < 8)  firstsize = sel->maxchunksize;
                    else                        firstsize = 8;
                    for (maxchunksize = firstsize;
                         maxchunksize <= sel->maxchunksize;
                         maxchunksize++) {
                        flags.maxchunksize = maxchunksize;
                        flags.spheresize = maxchunksize - 2;

                        CCutil_start_timer (&lp->stats.cuts_localcut);
                        rval = CCchunk_localcuts (&cuts, &cutcount,
                                 lp->graph.ncount, closecount, closelist,
                                 closex, 0.0, flags, &lc_timer, silent,
                                 rstate);
                        if (rval) {
                            fprintf (stderr, "LocalCuts failed\n");
                            rval = 1; goto CLEANUP;
                        }
                        z = CCutil_stop_timer (&lp->stats.cuts_localcut, 0);
                        if (!silent) {
                            printf ("Found %2d POLISHED LocalCuts in %.2f seconds\n",
                                     cutcount, z);
                            fflush (stdout);
                        }
                        if (cutcount) {
                            CCutil_start_timer (&lp->stats.cuts_localcut_opt);
                            rval = call_add_cuts (lp, &cuts, &cut_added,
                                    &xcount, &xlist, &x, &newval,
                                    sel->usetighten, &istour, silent,
                                    rstate);
                            if (rval) {
                                fprintf (stderr, "call_add_cuts failed\n");
                                goto CLEANUP;
                            }
                            CCutil_stop_timer (&lp->stats.cuts_localcut_opt, 0);
                            if (istour) goto OUT_LOOP;
                        }
                        if (newval >= oldval + sel->nexttol)  break;
                    }
                    lcimprove = newval - beforeval;
                }
            }
#endif /* POLISHED_CHUNKS */

OUT_LOOP:
            
            CC_IFFREE (xlist, int);
            CC_IFFREE (x, double);

            CCutil_start_timer (&lp->stats.sparse_edge_check);
            rval = sparse_edge_check (lp, &eginside, &edge_added,
                                      (double *) NULL, silent, rstate);
            if (rval) {
                fprintf (stderr, "sparse_edge_check failed\n");
                rval = 1; goto CLEANUP;
            }
            if (!silent) {
                CCutil_stop_timer (&lp->stats.sparse_edge_check, 1);
            } else {
                CCutil_stop_timer (&lp->stats.sparse_edge_check, 0);
            }
            
            if (savelp) {
                rval = CCtsp_write_probfile_sav (lp);
                if (rval) {
                    fprintf (stderr, "CCtsp_write_probfile_sav failed\n");
                    rval = 1; goto CLEANUP;
                }
            }
            if (lp->pool && savelp) {
                char buf[1024];
                if (!silent) {
                    printf ("Write Pool: %d cuts\n", lp->pool->cutcount);
                    fflush (stdout);
                }
                sprintf (buf, "%s.pul", lp->problabel);
                rval = CCtsp_write_cutpool (lp->graph.ncount, buf, lp->pool);
                if (rval) {
                    fprintf (stderr, "CCtsp_write_cutpool failed\n");
                    rval = 1; goto CLEANUP;
                }
            }
            if (lp->pool && sel->remotepool &&
                lp->pool->cutcount > lp->pool->savecount) {
                rval = CCtsp_send_newcuts (lp->graph.ncount, lp->pool,
                        sel->remotehost, sel->remoteport);
                if (rval) {
                    fprintf (stderr, "CCtsp_send_newcuts failed\n");
                    rval = 0;
                }
            }
            if (!silent) {
                CCutil_stop_timer (&lp->stats.cutting_inner_loop, 1);
            } else {
                CCutil_stop_timer (&lp->stats.cutting_inner_loop, 0);
            } 

            rval = lp_value (lp, &priceval);
            if (rval) {rval = 1; goto CLEANUP;}

            if (lp->lowerbound >= lp->upperbound - 0.9) {
                if (!silent) {
                    printf ("Stop cutting, lp bound is within 0.9 of upperbound\n");
                    fflush (stdout);
                }
                goto CLEANUP;
            }
            loopcount++;
            if (silent && !lp->full_edges_valid) {
                printf ("  LP Value %2d: %f  (%.2f seconds)\n", loopcount,
                     priceval, CCutil_zeit () - szeit);
                fflush (stdout);
            }
        } while ((newval > oldval + sel->roundtol ||
                  priceval < newval - sel->roundtol) &&
                 loopcount < LOOP_FULL &&
                 (lp->full_edges_valid || priceval < lp->upperbound));

        CCutil_start_timer (&lp->stats.full_edge_check);
        rval = full_edge_check (lp, &edge_added, silent, rstate);
        if (rval) {
            fprintf (stderr, "full_edge_check failed\n");
            rval = 1; goto CLEANUP;
        }
        if (!silent) {
            CCutil_stop_timer (&lp->stats.full_edge_check, 1);
        } else {
            CCutil_stop_timer (&lp->stats.full_edge_check, 0);
        }
        if (savelp) {
            rval = CCtsp_write_probfile_sav (lp);
            if (rval) {
                fprintf (stderr, "CCtsp_write_probfile_sav failed\n");
                rval = 1; goto CLEANUP;
            }
        }
        rval = lp_value (lp, &priceval);
        if (rval) {rval = 1; goto CLEANUP;}

        if (sel->extra_connect && priceval >= newval - sel->roundtol &&
            loopcount != LOOP_FULL) {
            if (!silent) {
                printf ("Check connectivity before exiting cutting_loop\n");
                fflush (stdout);
            }

            CCutil_start_timer (&lp->stats.cuts_extraconnect);
            
            rval = lp_x (lp, &xcount, &xlist, &x);
            if (rval) {rval = 1; goto CLEANUP;}

            rval = CCtsp_connect_cuts (&cuts, &cutcount_connect,
                        lp->graph.ncount, xcount, xlist, x);
            if (rval) {
                fprintf (stderr, "CCtsp_connect_cuts failed\n");
                rval = 1; goto CLEANUP;
            }
            z = CCutil_stop_timer (&lp->stats.cuts_extraconnect, 0);
            if (!silent) {
                printf ("Found %2d extra connect cuts in %.2f seconds\n",
                         cutcount_connect, z);
                fflush (stdout);
            }
            if (cutcount_connect) {
                CCutil_start_timer (&lp->stats.cuts_extraconnect_opt);
                rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                                      &xlist, &x, &newval, sel->usetighten,
                                      (int *) NULL, silent, rstate);
                if (rval) {
                    fprintf (stderr, "call_add_cuts failed\n");
                    goto CLEANUP;
                }
                CCutil_stop_timer (&lp->stats.cuts_extraconnect_opt, 0);
            }

            CC_FREE (xlist, int);
            CC_FREE (x, double);
        }
        outside++;

        /* This last bit will cause a second pass for large probs, with */
        /* an updated tolerance.                                        */

        if (!sel->fastcuts && lp->full_edges_valid == 0 && outside == 1 &&
            lp->graph.ncount >= 400 &&
            lp->lowerbound < lp->upperbound - 0.9) {
            rval = CCtsp_cutselect_set_tols (sel, lp, 1, silent);
            if (rval) {
                fprintf (stderr, "CCtsp_cutselect_set_tols failed\n");
                rval = 1;  goto CLEANUP;
            }
            loopcount = LOOP_FULL;  /* to run again */
        }
    } while (priceval < newval - sel->roundtol || loopcount == LOOP_FULL ||
             cutcount_connect);

CLEANUP:

    if (rval == 2) {
        printf ("LP is infeasible in cutting_loop\n");
        fflush (stdout);
    } else if (rval) {
        fprintf (stderr, "failure in cutting_loop\n");
    }
    if (!silent) {
        CCutil_stop_timer (&lp->stats.cutting_loop, 1);
        printf ("Number of outside rounds: %d\n", outside);
        fflush (stdout);
    } else {
        CCutil_stop_timer (&lp->stats.cutting_loop, 0);
    }

    if (eginside.ncount)
        CCtsp_free_edgegenerator (&eginside);
    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);
    CC_IFFREE (closelist, int);
    CC_IFFREE (closex, double);

    sel->nexttol = save_nexttol;
    sel->roundtol = save_roundtol;

    return rval;
}

#define CC_NO_NEAREST_SUBTOUR 50
#define CC_SUBTOUR_ROUNDS     5

int CCtsp_subtour_loop (CCtsp_lp *lp, int silent, CCrandstate *rstate)
{
    int rval = 0;
    int xcount, cutcount, cut_added, edge_added;
    int outside = 0;
    int inside = 0;
    int tighten = 0;
    double newval, priceval;
    int *xlist = (int *) NULL;
    double *x = (double *) NULL;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    CCtsp_edgegenerator eginside;
    double z;

    CCutil_start_timer (&lp->stats.cutting_loop);
    eginside.ncount = 0;
    if (lp->fulladj) {
        rval = CCtsp_init_edgegenerator (&eginside, lp->graph.ncount, lp->dat,
                                         lp->fulladj, 0, silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_init_edgegenerator (sparse) failed\n");
            rval = 1; goto CLEANUP;
        }
    } else if (lp->dat) {
        rval = CCtsp_init_edgegenerator (&eginside, lp->graph.ncount, lp->dat,
                (CCtsp_genadj *) NULL, CC_NO_NEAREST_SUBTOUR, silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_init_edgegenerator (sparse) failed\n");
            rval = 1; goto CLEANUP;
        }
    }

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL,
              (int *) NULL, (int **) NULL, (double **) NULL, (double **) NULL,
              (double **) NULL, (double **) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_get_lp_result failed\n");
        rval = 1; goto CLEANUP;
    }

    do {
        do {
            cut_added = 0;

            CCutil_start_timer (&lp->stats.cutting_inner_loop);

            rval = lp_x (lp, &xcount, &xlist, &x);
            if (rval) {rval = 1; goto CLEANUP;}

            /**** Connect Cuts ****/

            CCutil_start_timer (&lp->stats.cuts_connect);
            rval = CCtsp_connect_cuts (&cuts, &cutcount, lp->graph.ncount,
                                       xcount, xlist, x);
            if (rval) {
                fprintf (stderr, "CCtsp_connect_cuts failed\n");
                rval = 1; goto CLEANUP;
            }
            z = CCutil_stop_timer (&lp->stats.cuts_connect, 0);
            if (!silent) {
                printf ("Found %2d connect cuts in %.2f seconds\n",
                         cutcount, z);
                fflush (stdout);
            }
            if (cutcount) {
                CCutil_start_timer (&lp->stats.cuts_connect_opt);
                rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                         &xlist, &x, &newval, tighten, (int *) NULL, 
                         silent, rstate);
                if (rval) {
                    fprintf (stderr, "call_add_cuts failed\n"); goto CLEANUP;
                }
                CCutil_stop_timer (&lp->stats.cuts_connect_opt, 0);
            }

            /**** Shrink Cuts ****/

            CCutil_start_timer (&lp->stats.cuts_exactsubtour);
            rval = CCtsp_shrink_subtours (&cuts, &cutcount, lp->graph.ncount,
                                         xcount, xlist, x);
            if (rval) {
                fprintf (stderr, "CCtsp_shrink_subtours failed\n");
                rval = 1; goto CLEANUP;
            }
            z = CCutil_stop_timer (&lp->stats.cuts_exactsubtour, 0);
            if (!silent) {
                printf ("Found %2d shrink subtours in %.2f seconds\n",
                         cutcount, z);
                fflush (stdout);
            }
            if (cutcount) {
                CCutil_start_timer (&lp->stats.cuts_exactsubtour_opt);
                rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                          &xlist, &x, &newval, tighten, (int *) NULL,
                          silent, rstate);
                if (rval) {
                    fprintf (stderr, "call_add_cuts failed\n"); goto CLEANUP;
                }
                CCutil_stop_timer (&lp->stats.cuts_exactsubtour_opt, 0);
            }


            /**** Linear Cuts ****/

            CCutil_start_timer (&lp->stats.cuts_segment);
            rval = CCtsp_segment_cuts (&cuts, &cutcount, lp->graph.ncount,
                                      xcount, xlist, x);
            if (rval) {
                fprintf (stderr,  "CCtsp_segment_cuts failed\n");
                rval = 1; goto CLEANUP;
            }
            z = CCutil_stop_timer (&lp->stats.cuts_segment, 0);
            if (!silent) {
                printf ("Found %2d segment cuts in %.2f seconds\n",
                         cutcount, z);
                fflush (stdout);
            }
            if (cutcount) {
                CCutil_start_timer (&lp->stats.cuts_segment_opt);
                rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                          &xlist, &x, &newval, tighten, (int *) NULL, 
                          silent, rstate);
                if (rval) {
                    fprintf (stderr, "call_add_cuts failed\n"); goto CLEANUP;
                }
                CCutil_stop_timer (&lp->stats.cuts_segment_opt, 0);
            }


            /**** Exact Cuts ****/

            CCutil_start_timer (&lp->stats.cuts_exactsubtour);
            rval = CCtsp_exact_subtours (&cuts, &cutcount, lp->graph.ncount,
                                         xcount, xlist, x);
            if (rval) {
                fprintf (stderr, "CCtsp_exact_subtours failed\n");
                rval = 1; goto CLEANUP;
            }
            z = CCutil_stop_timer (&lp->stats.cuts_exactsubtour, 0);
            if (!silent) {
                printf ("Found %2d exact subtours in %.2f seconds\n",
                         cutcount, z);
                fflush (stdout);
            }
            if (cutcount) {
                CCutil_start_timer (&lp->stats.cuts_exactsubtour_opt);
                rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                          &xlist, &x, &newval, tighten, (int *) NULL,
                          silent, rstate);
                if (rval) {
                    fprintf (stderr, "call_add_cuts failed\n"); goto CLEANUP;
                }
                CCutil_stop_timer (&lp->stats.cuts_exactsubtour_opt, 0);
            }

            CC_FREE (xlist, int);
            CC_FREE (x, double);

            if (!cut_added || (inside % CC_SUBTOUR_ROUNDS) == 0) {
                CCutil_start_timer (&lp->stats.sparse_edge_check);
                rval = sparse_edge_check (lp, &eginside, &edge_added,
                                          (double *) NULL, silent, rstate);
                if (rval) {
                    fprintf (stderr, "sparse_edge_check failed\n");
                    rval = 1; goto CLEANUP;
                }
                if (!silent) {
                    CCutil_stop_timer (&lp->stats.sparse_edge_check, 1);
                } else {
                    CCutil_stop_timer (&lp->stats.sparse_edge_check, 0);
                }
            }
            if (!silent) {
                CCutil_stop_timer (&lp->stats.cutting_inner_loop, 1);
            } else {
                CCutil_stop_timer (&lp->stats.cutting_inner_loop, 0);
            }
            inside++;
            if (silent) {
                rval = lp_value (lp, &priceval);
                if (rval) {rval = 1; goto CLEANUP;}
                printf ("  LP Value %2d: %f\n", inside, priceval);
                fflush (stdout);
            }
        } while (edge_added || cut_added);

        CCutil_start_timer (&lp->stats.full_edge_check);
        rval = full_edge_check (lp, &edge_added, silent, rstate);
        if (rval) {
            fprintf (stderr, "full_edge_check failed\n");
            rval = 1; goto CLEANUP;
        }
        if (!silent) {
            CCutil_stop_timer (&lp->stats.full_edge_check, 1);
        } else {
            CCutil_stop_timer (&lp->stats.full_edge_check, 0);
        }
        outside++;
    } while (edge_added);

CLEANUP:

    if (rval == 2) {
        printf ("LP is infeasible in subtour_loop\n");
        fflush (stdout);
    } else if (rval) {
        fprintf (stderr, "failure in subtour_loop\n");
    }
    z = CCutil_stop_timer (&lp->stats.cutting_loop, 1);
    printf ("Time in cutting routine: %.2f\n", z);
    CCutil_total_timer (&lp->stats.cuts_connect, 1);
    CCutil_total_timer (&lp->stats.cuts_segment, 1);
    CCutil_total_timer (&lp->stats.cuts_exactsubtour, 1);

    printf ("Number of outside rounds: %d (%d inside)\n", outside, inside);
    fflush (stdout);

    if (eginside.ncount)
        CCtsp_free_edgegenerator (&eginside);
    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);

    return rval;
}

int CCtsp_blossom_loop (CCtsp_lp *lp, int silent, CCrandstate *rstate)
{
    int rval = 0;
    int xcount, cutcount, cut_added, edge_added;
    int outside = 0;
    int inside = 0;
    int tighten = 0;
    double newval, priceval;
    int *xlist = (int *) NULL;
    double *x = (double *) NULL;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    CCtsp_edgegenerator eginside;
    double z;

    CCutil_start_timer (&lp->stats.cutting_loop);
    eginside.ncount = 0;
    if (lp->fulladj) {
        rval = CCtsp_init_edgegenerator (&eginside, lp->graph.ncount, lp->dat,
                                         lp->fulladj, 0, silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_init_edgegenerator (sparse) failed\n");
            rval = 1; goto CLEANUP;
        }
    } else if (lp->dat) {
        rval = CCtsp_init_edgegenerator (&eginside, lp->graph.ncount, lp->dat,
                (CCtsp_genadj *) NULL, CC_NO_NEAREST_SUBTOUR, silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_init_edgegenerator (sparse) failed\n");
            rval = 1; goto CLEANUP;
        }
    }

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL,
              (int *) NULL, (int **) NULL, (double **) NULL, (double **) NULL,
              (double **) NULL, (double **) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_get_lp_result failed\n");
        rval = 1; goto CLEANUP;
    }

    do {
        do {
            cut_added = 0;

            CCutil_start_timer (&lp->stats.cutting_inner_loop);

            rval = lp_x (lp, &xcount, &xlist, &x);
            if (rval) {rval = 1; goto CLEANUP;}


            /****  Fast Blossoms ****/

            CCutil_start_timer (&lp->stats.cuts_fastblossom);
            rval = CCtsp_fastblossom (&cuts, &cutcount, lp->graph.ncount,
                                  xcount, xlist, x);
            if (rval) {
                fprintf (stderr, "CCtsp_fastblossom failed\n");
                rval = 1; goto CLEANUP;
            }
            z = CCutil_stop_timer (&lp->stats.cuts_fastblossom, 0);
            if (!silent) {
                printf ("Found %2d Fast Blossoms in %.2f seconds\n",
                         cutcount, z);
                fflush (stdout);
            }
            if (cutcount) {
                CCutil_start_timer (&lp->stats.cuts_fastblossom_opt);
                rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                          &xlist, &x, &newval, tighten, (int *) NULL,
                          silent, rstate);
                if (rval) {
                    fprintf (stderr, "call_add_cuts failed\n"); goto CLEANUP;
                }
                CCutil_stop_timer (&lp->stats.cuts_fastblossom_opt, 0);
            }


            /****  Groetschel-Holland Fast Blossoms ****/
 
            CCutil_start_timer (&lp->stats.cuts_ghfastblossom);
            rval = CCtsp_ghfastblossom (&cuts, &cutcount, lp->graph.ncount,
                                  xcount, xlist, x);
            if (rval) {
                fprintf (stderr, "CCtsp_ghfastblossom failed\n");
                rval = 1; goto CLEANUP;
            }
            z = CCutil_stop_timer (&lp->stats.cuts_ghfastblossom, 0);
            if (!silent) {
                printf ("Found %2d Groetschel-Holland Blossoms in %.2f seconds\n",
                    cutcount, z);
                fflush (stdout);
            }
            if (cutcount) {
                CCutil_start_timer (&lp->stats.cuts_ghfastblossom_opt);
                rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                          &xlist, &x, &newval, tighten, (int *) NULL,
                          silent, rstate);
                if (rval) {
                    fprintf (stderr, "call_add_cuts failed\n"); goto CLEANUP;
                }
                CCutil_stop_timer (&lp->stats.cuts_ghfastblossom_opt, 0);
            }


            /**** Exact Blossoms ****/

            CCutil_start_timer (&lp->stats.cuts_exactblossom);
            rval = CCtsp_exactblossom (&cuts, &cutcount, lp->graph.ncount,
                                  xcount, xlist, x, rstate);
            if (rval) {
                fprintf (stderr, "CCtsp_exactblossom failed\n");
                rval = 1; goto CLEANUP;
            }
            z = CCutil_stop_timer (&lp->stats.cuts_exactblossom, 0);
            if (!silent) {
                printf ("Found %2d Exact Blossoms in %.2f seconds\n",
                         cutcount, z);
                fflush (stdout);
            }
            if (cutcount) {
                CCutil_start_timer (&lp->stats.cuts_exactblossom_opt);
                rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                          &xlist, &x, &newval, tighten, (int *) NULL,
                          silent, rstate);
                if (rval) {
                    fprintf (stderr, "call_add_cuts failed\n"); goto CLEANUP;
                }
                CCutil_stop_timer (&lp->stats.cuts_exactblossom_opt, 0);
            }

            CC_FREE (xlist, int);
            CC_FREE (x, double);

            if (!cut_added || (inside % CC_SUBTOUR_ROUNDS) == 0) {
                CCutil_start_timer (&lp->stats.sparse_edge_check);
                rval = sparse_edge_check (lp, &eginside, &edge_added,
                                          (double *) NULL, silent, rstate);
                if (rval) {
                    fprintf (stderr, "sparse_edge_check failed\n");
                    rval = 1; goto CLEANUP;
                }
                if (!silent) {
                    CCutil_stop_timer (&lp->stats.sparse_edge_check, 1);
                } else {
                    CCutil_stop_timer (&lp->stats.sparse_edge_check, 0);
                }
            }
            if (!silent) {
                CCutil_stop_timer (&lp->stats.cutting_inner_loop, 1);
            } else {
                CCutil_stop_timer (&lp->stats.cutting_inner_loop, 0);
            }
            inside++;
            if (silent) {
                rval = lp_value (lp, &priceval);
                if (rval) {rval = 1; goto CLEANUP;}
                printf ("  LP Value %2d: %f\n", inside, priceval);
                fflush (stdout);
            }
        } while (edge_added || cut_added);

        CCutil_start_timer (&lp->stats.full_edge_check);
        rval = full_edge_check (lp, &edge_added, silent, rstate);
        if (rval) {
            fprintf (stderr, "full_edge_check failed\n");
            rval = 1; goto CLEANUP;
        }
        if (!silent) {
            CCutil_stop_timer (&lp->stats.full_edge_check, 1);
        } else {
            CCutil_stop_timer (&lp->stats.full_edge_check, 0);
        }
        outside++;
    } while (edge_added);

CLEANUP:

    if (rval == 2) {
        printf ("LP is infeasible in blossom_loop\n");
        fflush (stdout);
    } else if (rval) {
        fprintf (stderr, "failure in blossom_loop\n");
    }
    z = CCutil_stop_timer (&lp->stats.cutting_loop, 1);
    printf ("Time in cutting routine: %.2f\n", z);
    CCutil_total_timer (&lp->stats.cuts_fastblossom, 1);
    CCutil_total_timer (&lp->stats.cuts_exactblossom, 1);

    printf ("Number of outside rounds: %d (%d inside)\n", outside, inside);
    fflush (stdout);

    if (eginside.ncount)
        CCtsp_free_edgegenerator (&eginside);
    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);

    return rval;
}

int CCtsp_subtour_and_blossom_loop (CCtsp_lp *lp, int silent,
        CCrandstate *rstate)
{
    int rval = 0;
    int xcount, cutcount, cut_added, edge_added, blossom_added;
    int outside = 0;
    int inside = 0;
    int tighten = 0;
    double newval, priceval;
    int *xlist = (int *) NULL;
    double *x = (double *) NULL;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    CCtsp_edgegenerator eginside;
    double z;

    CCutil_start_timer (&lp->stats.cutting_loop);
    eginside.ncount = 0;
    if (lp->fulladj) {
        rval = CCtsp_init_edgegenerator (&eginside, lp->graph.ncount, lp->dat,
                                         lp->fulladj, 0, silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_init_edgegenerator (sparse) failed\n");
            rval = 1; goto CLEANUP;
        }
    } else if (lp->dat) {
        rval = CCtsp_init_edgegenerator (&eginside, lp->graph.ncount, lp->dat,
                (CCtsp_genadj *) NULL, CC_NO_NEAREST_SUBTOUR, silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_init_edgegenerator (sparse) failed\n");
            rval = 1; goto CLEANUP;
        }
    }

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL,
              (int *) NULL, (int **) NULL, (double **) NULL, (double **) NULL,
              (double **) NULL, (double **) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_get_lp_result failed\n");
        rval = 1; goto CLEANUP;
    }

    do {
        do {
            cut_added = blossom_added = 0;

            CCutil_start_timer (&lp->stats.cutting_inner_loop);

            rval = lp_x (lp, &xcount, &xlist, &x);
            if (rval) {rval = 1; goto CLEANUP;}


            /**** Connect Cuts ****/

            CCutil_start_timer (&lp->stats.cuts_connect);
            rval = CCtsp_connect_cuts (&cuts, &cutcount, lp->graph.ncount,
                                       xcount, xlist, x);
            if (rval) {
                fprintf (stderr, "CCtsp_connect_cuts failed\n");
                rval = 1; goto CLEANUP;
            }
            z = CCutil_stop_timer (&lp->stats.cuts_connect, 0);
            if (!silent) {
                printf ("Found %2d connect cuts in %.2f seconds\n",
                         cutcount, z);
                fflush (stdout);
            }
            if (cutcount) {
                CCutil_start_timer (&lp->stats.cuts_connect_opt);
                rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                         &xlist, &x, &newval, tighten, (int *) NULL,
                          silent, rstate);
                if (rval) {
                    fprintf (stderr, "call_add_cuts failed\n"); goto CLEANUP;
                }
                CCutil_stop_timer (&lp->stats.cuts_connect_opt, 0);
            }


            /**** Linear Cuts ****/

            CCutil_start_timer (&lp->stats.cuts_segment);
            rval = CCtsp_segment_cuts (&cuts, &cutcount, lp->graph.ncount,
                                      xcount, xlist, x);
            if (rval) {
                fprintf (stderr,  "CCtsp_segment_cuts failed\n");
                rval = 1; goto CLEANUP;
            }
            z = CCutil_stop_timer (&lp->stats.cuts_segment, 0);
            if (!silent) {
                printf ("Found %2d segment cuts in %.2f seconds\n",
                         cutcount, z);
                fflush (stdout);
            }
            if (cutcount) {
                CCutil_start_timer (&lp->stats.cuts_segment_opt);
                rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                          &xlist, &x, &newval, tighten, (int *) NULL, 
                          silent,rstate);
                if (rval) {
                    fprintf (stderr, "call_add_cuts failed\n"); goto CLEANUP;
                }
                CCutil_stop_timer (&lp->stats.cuts_segment_opt, 0);
            }


            /**** Exact Cuts ****/

            CCutil_start_timer (&lp->stats.cuts_exactsubtour);
            rval = CCtsp_exact_subtours (&cuts, &cutcount, lp->graph.ncount,
                                         xcount, xlist, x);
            if (rval) {
                fprintf (stderr, "CCtsp_exact_subtours failed\n");
                rval = 1; goto CLEANUP;
            }
            z = CCutil_stop_timer (&lp->stats.cuts_exactsubtour, 0);
            if (!silent) {
                printf ("Found %2d exact subtours in %.2f seconds\n",
                         cutcount, z);
                fflush (stdout);
            }
            if (cutcount) {
                CCutil_start_timer (&lp->stats.cuts_exactsubtour_opt);
                rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                          &xlist, &x, &newval, tighten, (int *) NULL,
                          silent, rstate);
                if (rval) {
                    fprintf (stderr, "call_add_cuts failed\n"); goto CLEANUP;
                }
                CCutil_stop_timer (&lp->stats.cuts_exactsubtour_opt, 0);
            }

            /****  Fast Blossoms ****/

            CCutil_start_timer (&lp->stats.cuts_fastblossom);
            rval = CCtsp_fastblossom (&cuts, &cutcount, lp->graph.ncount,
                                  xcount, xlist, x);
            if (rval) {
                fprintf (stderr, "CCtsp_fastblossom failed\n");
                rval = 1; goto CLEANUP;
            }
            z = CCutil_stop_timer (&lp->stats.cuts_fastblossom, 0);
            if (!silent) {
                printf ("Found %2d Fast Blossoms in %.2f seconds\n",
                         cutcount, z);
                fflush (stdout);
            }
            if (cutcount) {
                CCutil_start_timer (&lp->stats.cuts_fastblossom_opt);
                rval = call_add_cuts (lp, &cuts, &blossom_added, &xcount,
                          &xlist, &x, &newval, tighten, (int *) NULL,
                          silent, rstate);
                if (rval) {
                    fprintf (stderr, "call_add_cuts failed\n"); goto CLEANUP;
                }
                CCutil_stop_timer (&lp->stats.cuts_fastblossom_opt, 0);
            }


            CC_FREE (xlist, int);
            CC_FREE (x, double);

            if (!cut_added || (inside % CC_SUBTOUR_ROUNDS) == 0) {
                CCutil_start_timer (&lp->stats.sparse_edge_check);
                rval = sparse_edge_check (lp, &eginside, &edge_added,
                                          (double *) NULL, silent, rstate);
                if (rval) {
                    fprintf (stderr, "sparse_edge_check failed\n");
                    rval = 1; goto CLEANUP;
                }
                if (!silent) {
                    CCutil_stop_timer (&lp->stats.sparse_edge_check, 1);
                } else {
                    CCutil_stop_timer (&lp->stats.sparse_edge_check, 0);
                }
            }
            if (!silent) {
                CCutil_stop_timer (&lp->stats.cutting_inner_loop, 1);
            } else {
                CCutil_stop_timer (&lp->stats.cutting_inner_loop, 0);
            }
            inside++;
            if (silent) {
                rval = lp_value (lp, &priceval);
                if (rval) {rval = 1; goto CLEANUP;}
                printf ("  LP Value %2d: %f\n", inside, priceval);
                fflush (stdout);
            }
        } while (edge_added || cut_added || blossom_added);

        CCutil_start_timer (&lp->stats.full_edge_check);
        rval = full_edge_check (lp, &edge_added, silent, rstate);
        if (rval) {
            fprintf (stderr, "full_edge_check failed\n");
            rval = 1; goto CLEANUP;
        }
        if (!silent) {
            CCutil_stop_timer (&lp->stats.full_edge_check, 1);
        } else {
            CCutil_stop_timer (&lp->stats.full_edge_check, 0);
        }
        outside++;

    } while (edge_added);

CLEANUP:

    if (rval == 2) {
        printf ("LP is infeasible in subtour_and_blossom_loop\n");
        fflush (stdout);
    } else if (rval) {
        fprintf (stderr, "failure in subtour_and_blossom_loop\n");
    }
    z = CCutil_stop_timer (&lp->stats.cutting_loop, 1);
    printf ("Time in cutting routine: %.2f\n", z);
    CCutil_total_timer (&lp->stats.cuts_connect, 1);
    CCutil_total_timer (&lp->stats.cuts_segment, 1);
    CCutil_total_timer (&lp->stats.cuts_exactsubtour, 1);
    CCutil_total_timer (&lp->stats.cuts_fastblossom, 1);
    

    printf ("Number of outside rounds: %d (%d inside)\n", outside, inside);
    fflush (stdout);

    if (eginside.ncount)
        CCtsp_free_edgegenerator (&eginside);
    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);

    return rval;
}

static int call_add_cuts (CCtsp_lp *lp, CCtsp_lpcut_in **cuts, int *cut_added,
        int *xcount, int **xlist, double **x, double *val, int tighten,
        int *istour, int silent, CCrandstate *rstate)
{
    int rval = 0;
    double dval;
    double szeit = CCutil_zeit ();

    if (istour) *istour = 0;

    CC_IFFREE (*xlist, int);
    CC_IFFREE (*x, double);

    CCtsp_add_cuts_to_queue (lp, cuts);
    rval = CCtsp_process_cuts (lp, cut_added, tighten, silent, rstate);
    if (rval) {
        fprintf (stderr, "process_cuts failed\n"); goto CLEANUP;
    }

    rval = lp_value (lp, val);
    if (rval) {
        fprintf (stderr, "lp_value failed\n"); rval = 1; goto CLEANUP;
    }
    if (!silent) {
        printf ("  Add %2d cuts (Total %d), LP: %f (%.2f seconds)\n",
                       *cut_added, lp->cuts.cutcount, *val,
                        CCutil_zeit () - szeit);
        fflush (stdout);
    }

    rval = lp_x (lp, xcount, xlist, x);
    if (rval) {
        fprintf (stderr, "lp_x failed\n"); rval = 1; goto CLEANUP;
    }

    if (istour) {
        rval = CCtsp_check_integral (lp, &dval, (int **) NULL, istour,
                                     silent);
        if (rval) {
            fprintf (stderr, "CCtsp_check_integral failed\n"); goto CLEANUP;
        }
    } 

CLEANUP:

    return rval;
}

static int lp_value (CCtsp_lp *lp, double *val)
{
    int rval;

    rval = CCtsp_get_lp_result (lp, val, (double *) NULL, (int *) NULL,
                 (int **) NULL, (double **) NULL, (double **) NULL,
                 (double **) NULL, (double **) NULL);
    if (rval) fprintf (stderr, "CCtsp_get_lp_result failed\n");
    return rval;
}

static int lp_x (CCtsp_lp *lp, int *xcount, int **xlist, double **x)
{
    int rval;

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL, xcount,
                     xlist, x, (double **) NULL, (double **) NULL,
                     (double **) NULL);
    if (rval) fprintf (stderr, "CCtsp_get_lp_result failed\n");
    return rval;
}

int CCtsp_pricing_loop (CCtsp_lp *lp, double *bnd, int silent,
        CCrandstate *rstate)
{
    CCtsp_edgegenerator eg;
    int nadded;
    int rval = 0;

    eg.ncount = 0;
    if (!lp->full_edges_valid) {
        fprintf (stderr, "CCtsp_pricing_loop called without valid edges\n");
        rval = 1; goto CLEANUP;
    }


    rval = CCtsp_init_edgegenerator (&eg, lp->graph.ncount, lp->dat,
                                     lp->fulladj, 0, silent, rstate);
    if (rval) {
        fprintf (stderr, "CCtsp_init_edgegenerator failed\n"); goto CLEANUP;
    }
    rval = sparse_edge_check (lp, &eg, &nadded, bnd, silent, rstate);
    if (rval) {
        fprintf (stderr, "sparse_edge_check failed\n"); goto CLEANUP;
    }

CLEANUP:

    if (eg.ncount) {
        CCtsp_free_edgegenerator (&eg);
    }
    return rval;
}

static int full_edge_check (CCtsp_lp *lp, int *nadded, int silent,
        CCrandstate *rstate)
{
    int rval;
    CCtsp_edgegenerator eg;
    double val, penalty;

    if (lp->dat && (!lp->full_edges_valid)) {
        rval = CCtsp_init_edgegenerator (&eg, lp->graph.ncount, lp->dat,
                    (CCtsp_genadj *) NULL, CCtsp_PRICE_COMPLETE_GRAPH, silent,
                    rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_init_edgegenerator failed\n"); return rval;
        }

        rval = CCtsp_addbad_variables (lp, &eg, &penalty, nadded,
                      CCtsp_PRICE_RCTHRESH, CCtsp_PRICE_MAXPENALTY, 0,
                      (int *) NULL, silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_addbad_variables failed\n");
            CCtsp_free_edgegenerator (&eg);
            return rval;
        }
        CCtsp_free_edgegenerator (&eg);
        if (!silent) {
            printf ("%d edges added, penalty %f\n", *nadded, penalty);
            fflush (stdout);
        }

        rval = lp_value (lp, &val);
        if (rval) return rval;

        if (val + penalty > lp->lowerbound) {
            printf ("New lower bound: %f\n", val+ penalty);
            fflush (stdout);
            lp->lowerbound = val + penalty;
        }
    } else {
        *nadded = 0;
    }
    return 0;
}

static int sparse_edge_check (CCtsp_lp *lp, CCtsp_edgegenerator *eg,
        int *nadded, double *bnd, int silent, CCrandstate *rstate)
{
    double val, penalty;
    int rval;

    if (bnd) *bnd = -CCtsp_LP_MAXDOUBLE;

    if (eg->ncount > 0) {
        rval = CCtsp_addbad_variables (lp, eg, &penalty, nadded,
                  CCtsp_PRICE_RCTHRESH, CCtsp_PRICE_MAXPENALTY, 0,
                  (int *) NULL, silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_addbad_variables failed\n"); return rval;
        }

        rval = lp_value (lp, &val);
        if (rval) { fprintf (stderr, "lp_value failed\n"); return rval; }

        if (!silent) {
            printf ("(SPARSE) %d edges added, penalty %f, val %f\n",
                      *nadded, penalty, val);
            fflush (stdout);
        }

        if (lp->full_edges_valid) {
            if (val + penalty > lp->lowerbound) {
                if (!silent) {
                    printf ("New (node) lower bound: %f\n", val + penalty);
                    fflush (stdout);
                }
                lp->lowerbound = val + penalty;
            }
            if (bnd) *bnd = val + penalty;
        }
    } else {
        *nadded = 0;
    }
    return 0;
}

int CCtsp_bb_cutting (char *probname, int probnum, int prob_newnum, int ncount,
        CCdatagroup *dat, int *ptour, double *upbound, CCtsp_lpcuts *pool,
        CCtsp_cutselect *sel, double *val, int *prune, int *foundtour,
        int *besttour, int level, int silent, CCrandstate *rstate)
{
    int rval = 0;
    CCtsp_lp *lp = (CCtsp_lp *) NULL;
    double cval, tourval;
    int test;

    *val = 0.0;
    *prune = 0;
    *foundtour = 0;

    rval = bb_cutting_work (&lp, probname, probnum, ncount, dat, ptour,
                  *upbound, pool, sel, &cval, level, silent, rstate);
    if (rval) {
        fprintf (stderr, "bb_cutting_work failed\n"); fflush (stdout);
        goto CLEANUP;
    }

    if (lp != (CCtsp_lp *) NULL) {
        lp->id = prob_newnum;
    }

    if (cval == CCtsp_LP_MAXDOUBLE) {
        rval = CCtsp_verify_infeasible_lp (lp, &test, silent);
        if (rval) {
            fprintf (stderr, "CCtsp_verify_infeasible_lp failed\n");
            goto CLEANUP;
        }
        if (test) {
            printf ("verified infeasible LP\n"); fflush (stdout);
            *val = CCtsp_LP_MAXDOUBLE;
            *prune = 1;
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
        rval = CCtsp_pricing_loop (lp, val, silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_pricing_loop failed\n");
            rval = 1; goto CLEANUP;
        }
        lp->lowerbound = *val;
        if (lp->upperbound < *upbound) *upbound = lp->upperbound;

        if (lp->lowerbound < lp->upperbound - 0.9) {
            CCutil_start_timer (&lp->stats.linkern);
            rval = CCtsp_call_x_heuristic (lp, &tourval, besttour, silent,
                                           rstate);
            if (rval) {
                fprintf (stderr, "CCtsp_call_x_heuristic failed\n");
                goto CLEANUP;
            }
            if (!silent) {
                CCutil_stop_timer (&lp->stats.linkern, 1);
            } else {
                CCutil_stop_timer (&lp->stats.linkern, 0);
            }
            if (tourval < lp->upperbound) {
                printf ("New upperbound from x-heuristic: %.2f\n", tourval);
                lp->upperbound = tourval;
                *upbound = tourval;
                *foundtour = 1;
            }
        }

        if (lp->lowerbound >= lp->upperbound - 0.9) {
            rval = CCtsp_verify_lp_prune (lp, &test,  silent);
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
                rval = CCtsp_write_probleaf_id (lp);
                if (rval) {
                    fprintf (stderr, "CCtsp_write_probleaf_id failed\n");
                    goto CLEANUP;
                }
            } else {
                printf ("exact pricing could not prune the search\n");
                fflush (stdout);
                rval = CCtsp_write_probfile_id (lp);
                CCcheck_rval (rval, "CCtsp_write_probfile_id failed");
            }
        } else {
            rval = CCtsp_write_probfile_id (lp);
            CCcheck_rval (rval, "CCtsp_write_probfile_id failed");
        }
    }

CLEANUP:

    if (lp) CCtsp_free_tsp_lp_struct (&lp);
    return rval;
}

int CCtsp_call_x_heuristic (CCtsp_lp *lp, double *val, int *outcyc,
        int silent, CCrandstate *rstate)
{
    int rval = 0;
    int *cyc   = (int *) NULL;
    int *xlist = (int *) NULL;
    double *x   = (double *) NULL;
    int ncount = lp->graph.ncount;
    int xcount, i;

    *val = CCtsp_LP_MAXDOUBLE;

    if (!lp->dat) goto CLEANUP;

    cyc = CC_SAFE_MALLOC (ncount, int);
    if (!cyc) {
        fprintf (stderr, "out of memory for cycle\n");
        rval = 1; goto CLEANUP;
    }
    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL,
         &xcount, &xlist, &x, (double **) NULL, (double **) NULL,
         (double **) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_get_lp_result failed\n");
        goto CLEANUP;
    }

    rval = CCtsp_x_greedy_tour_lk (lp->dat, ncount, xcount, xlist, x,
                   cyc, val, silent, rstate);
    if (rval) {
        fprintf (stderr, "CCtsp_x_greedy_tour_lk failed\n"); goto CLEANUP;
    }
    if (!silent) {
        printf ("x-heuristic lk  gives: %.2f\n", *val); fflush (stdout);
    }
    if (*val < lp->upperbound) {
        if (outcyc) {
            for (i = 0; i < ncount; i++) {
                outcyc[i] = cyc[i];
            }
        }
    }

CLEANUP:

    CC_IFFREE (cyc, int);
    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);
    return rval;
}

static int bb_cutting_work (CCtsp_lp **lp, char *probname, int probnum,
        int ncount, CCdatagroup *dat, int *ptour, double initial_ub,
        CCtsp_lpcuts *pool, CCtsp_cutselect *sel, double *val, int level,
        int silent, CCrandstate *rstate)
{
    int rval = 0;

    *lp = (CCtsp_lp *) NULL;
    *val = 0.0;

    rval = CCtsp_bb_init_lp (lp, probname, probnum, ncount, dat, ptour,
               initial_ub, pool, silent, rstate);
    if (rval == 2) {
        printf ("LP is reported to be infeasible\n"); fflush (stdout);
        *val = CCtsp_LP_MAXDOUBLE;
        rval = 0; goto CLEANUP;
    } else if (rval) {
        fprintf (stderr, "CCtsp_bb_init_lp failed\n"); goto CLEANUP;
    }
    CCutil_start_timer (&(*lp)->stats.total);

    if ((*lp)->lowerbound >= (*lp)->upperbound - 0.9) {
        printf ("Do not cut, the lp is within 1.0 of the upperbound\n");
        fflush (stdout);
        *val = (*lp)->lowerbound;
        goto CLEANUP;
    } else {
        rval = CCtsp_cutselect_set_tols (sel, *lp, level, silent);
        if (rval) {
            fprintf (stderr, "CCtsp_cutselect_set_tols failed\n");
            goto CLEANUP;
        }
        
        rval = CCtsp_cutting_loop (*lp, sel, 0, silent, rstate);
        if (rval == 2) {
            printf ("Cut LP is reported to be infeasible\n"); fflush (stdout);
            *val = CCtsp_LP_MAXDOUBLE;
            rval = 0;
        } else if (rval) {
            fprintf (stderr, "CCtsp_cutting_loop failed\n"); goto CLEANUP;
        } else {
            *val = (*lp)->lowerbound;
        }
    }

CLEANUP:

    if (!silent) {
        CCutil_stop_timer (&(*lp)->stats.total, 1);
    } else {
        CCutil_stop_timer (&(*lp)->stats.total, 0);
    }
    /* CCtsp_output_statistics (&(*lp)->stats); */

    if (!silent) {
        printf ("Final LP has %d rows, %d columns, %d nonzeros\n",
                CClp_nrows ((*lp)->lp), CClp_ncols ((*lp)->lp),
                CClp_nnonzeros ((*lp)->lp));
        fflush (stdout);
    }
    
    return rval;
}

static int grab_close_x (int ncount, int xcount, int *xlist, double *x,
        int *newcount, int **newlist, double **newx, double mult)
{
    int rval = 0;
    char *marks = (char *) NULL;
    int i, k, tmp, n1, n2;

    CC_IFFREE (*newlist, int);
    CC_IFFREE (*newx, double);

    *newx    = CC_SAFE_MALLOC (xcount + ncount, double);
    *newlist = CC_SAFE_MALLOC (2 * (xcount + ncount), int);
    marks    = CC_SAFE_MALLOC (ncount, char);
    if (!(*newx) || !(*newlist) || !marks) {
        fprintf (stderr, "out of memory in grab_close_x\n");
        CC_IFFREE (*newx, double);
        CC_IFFREE (*newlist, int);
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ncount; i++) {
        marks[i] = 0;
    }
    k = 0;
    for (i = 0; i < xcount; i++) {
        n1 = xlist[2*i];
        n2 = xlist[2*i+1];
        if (n2 < n1) {
            tmp = n1; n1 = n2; n2 = tmp;
        }
        if (n1 == (n2 - 1)) {
            (*newx)[i] = mult * x[i] + (1.0 - mult);
            marks[n1] = 1;
        } else if (n1 == 0 && n2 == ncount - 1) {
            (*newx)[i] = mult * x[i] + (1.0 - mult);
            marks[n2] = 1;
        } else {
            (*newx)[i] = mult * x[i];
        }
        (*newlist)[k++] = n1;
        (*newlist)[k++] = n2;
    }
    *newcount = xcount;
    for (i = 0; i < ncount; i++) {
        if (marks[i] == 0) {
            (*newx)[(*newcount)++] = 1.0 - mult;
            (*newlist)[k++] = i;
            (*newlist)[k++] = (i + 1) % ncount;
        }
    }

CLEANUP:

    CC_IFFREE (marks, char);
    return rval;
}

static int CC_UNUSED grab_polished_x (CCtsp_lp *lp, double dust_val,
        int *newcount, int **newlist, double **newx)
{
    int rval = 0;
    CClp_warmstart *warmstart = (CClp_warmstart *) NULL;
    double *x = (double *) NULL;
    char *marks = (char *) NULL;
    int i;
    int nset;
    int ncols;
    double szeit;
    double objval;
    double origval;

    CC_IFFREE (*newlist, int);
    CC_IFFREE (*newx, double);

    ncols = CClp_ncols (lp->lp);
    
    x = CC_SAFE_MALLOC (ncols, double);
    marks = CC_SAFE_MALLOC (ncols, char);
    if (x == (double *) NULL || marks == (char *) NULL) {
        fprintf (stderr, "Out of memory in grab_polished_x\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<ncols; i++) {
        if (lp->graph.edges[i].fixed || lp->graph.edges[i].branch) {
            marks[i] = -1;
        } else {
            marks[i] = 0;
        }
    }
    
    rval = CClp_get_warmstart (lp->lp, &warmstart);
    if (rval) {
        fprintf (stderr, "CClp_get_warmstart failed\n"); goto CLEANUP;
    }

    rval = CClp_objval (lp->lp, &origval);
    if (rval) {
        fprintf (stderr, "CClp_objval failed\n"); goto CLEANUP;
    }

    printf ("Polishing, objval %.6f\n", origval);
    fflush (stdout);
    do {
        szeit = CCutil_zeit();
        
        rval = CClp_x (lp->lp, x);
        if (rval) {
            fprintf (stderr, "CClp_x failed\n"); goto CLEANUP;
        }

        nset = 0;
        for (i=0; i<ncols; i++) {
            if (x[i] != 0.0 && x[i] < dust_val && marks[i] == 0) {
                marks[i] = 1;
                rval = CClp_setbnd (lp->lp, i, 'U', 0.0);
                if (rval) {
                    fprintf (stderr, "CClp_setbnd failed\n"); goto CLEANUP;
                }
                nset++;
            } else if (x[i] != 1.0 && x[i] > 1.0-dust_val && marks[i] == 0) {
                marks[i] = 2;
                rval = CClp_setbnd (lp->lp, i, 'L', 1.0);
                if (rval) {
                    fprintf (stderr, "CClp_setbnd failed\n"); goto CLEANUP;
                }
                nset++;
            }
        }

        if (nset > 0) {
            rval = CClp_opt (lp->lp, CClp_METHOD_DUAL);
            if (rval == 2) {
                fprintf (stderr, "Polished LP infeasible\n"); goto CLEANUP;
            } else if (rval) {
                fprintf (stderr, "CClp_opt failed\n"); goto CLEANUP;
            }
        }

        rval = CClp_objval (lp->lp, &objval);
        if (rval) {
            fprintf (stderr, "CClp_objval failed\n"); goto CLEANUP;
        }
        printf ("polished away %d dusty edges in %.2f seconds, objval %.6f\n",
                nset, CCutil_zeit() - szeit, objval);
        fflush (stdout);
    } while (nset > 0);

    for (i=0; i<ncols; i++) {
        if (marks[i] == 1) {
            rval = CClp_setbnd (lp->lp, i, 'U', 1.0);
            if (rval) {
                fprintf (stderr, "CClp_setbnd failed\n"); goto CLEANUP;
            }
        } else if (marks[i] == 2) {
            rval = CClp_setbnd (lp->lp, i, 'L', 0.0);
            if (rval) {
                fprintf (stderr, "CClp_setbnd failed\n"); goto CLEANUP;
            }
        }
    }
            
    rval = CClp_load_warmstart (lp->lp, warmstart);
    if (rval) {
        fprintf (stderr, "CClp_load_warmstart failed\n"); goto CLEANUP;
    }
    
    rval = CClp_opt (lp->lp, CClp_METHOD_DUAL);
    if (rval == 2) {
        fprintf (stderr, "restored LP infeasible\n"); goto CLEANUP;
    } else if (rval) {
        fprintf (stderr, "CClp_opt failed\n"); goto CLEANUP;
    }
    
    *newx    = x;
    x = (double *) NULL;

    *newlist = CC_SAFE_MALLOC (2 * ncols, int);
    if ((*newlist) == (int *) NULL) {
        fprintf (stderr, "out of memory in grab_polished_x\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ncols; i++) {
        (*newlist)[2*i] = lp->graph.edges[i].ends[0];
        (*newlist)[2*i+1] = lp->graph.edges[i].ends[1];
    }

    *newcount = ncols;
    rval = 0;

 CLEANUP:
    CClp_free_warmstart (&warmstart);
    
    CC_IFFREE (x, double);
    CC_IFFREE (marks, char);

    if (rval) {
        CC_IFFREE (*newlist, int);
        CC_IFFREE (*newx, double);
    }
    return rval;
}

static int CC_UNUSED grab_polished2_x (CCtsp_lp *lp, double dust_val,
        int *newcount, int **newlist, double **newx)
{
    int rval = 0;
    CClp_warmstart *warmstart = (CClp_warmstart *) NULL;
    double *x = (double *) NULL;
    char *marks = (char *) NULL;
    int i;
    int ncols;
    double szeit;
    double objval;
    double origval;
    double setval;
    int ninfeas = 0;
    int nfixed = 0;

    CC_IFFREE (*newlist, int);
    CC_IFFREE (*newx, double);

    ncols = CClp_ncols (lp->lp);
    
    x = CC_SAFE_MALLOC (ncols, double);
    marks = CC_SAFE_MALLOC (ncols, char);
    if (x == (double *) NULL || marks == (char *) NULL) {
        fprintf (stderr, "Out of memory in grab_polished_x\n");
        rval = 1; goto CLEANUP;
    }

    rval = CClp_get_warmstart (lp->lp, &warmstart);
    if (rval) {
        fprintf (stderr, "CClp_get_warmstart failed\n"); goto CLEANUP;
    }

    rval = CClp_objval (lp->lp, &origval);
    if (rval) {
        fprintf (stderr, "CClp_objval failed\n"); goto CLEANUP;
    }

    printf ("Polishing, objval %.6f\n", origval);
    fflush (stdout);

    rval = CClp_x (lp->lp, x);
    if (rval) {
        fprintf (stderr, "CClp_x failed\n"); goto CLEANUP;
    }

    for (i=0; i<ncols; i++) {
        if (lp->graph.edges[i].fixed || lp->graph.edges[i].branch) {
            marks[i] = -1;
        } else {
            if (x[i] - dust_val > 0.0) {
                rval = CClp_setbnd (lp->lp, i, 'L', x[i] - dust_val);
                if (rval) {
                    fprintf (stderr, "CClp_setbnd failed\n"); goto CLEANUP;
                }
            }
            if (x[i] + dust_val < 1.0) {
                rval = CClp_setbnd (lp->lp, i, 'U', x[i] + dust_val);
                if (rval) {
                    fprintf (stderr, "CClp_setbnd failed\n"); goto CLEANUP;
                }
            }
        }
    }

    szeit = CCutil_zeit();
    for (i=0; i<ncols; i++) {
      if (marks[i] == 0 && (x[i] < dust_val || x[i] > 1.0 - dust_val)) {
        printf ("edge %d oldval %.6f", i, x[i]); fflush (stdout);
        if (x[i] < dust_val) setval = 0.0;
        else                 setval = 1.0;
        rval = CClp_setbnd (lp->lp, i, 'L', setval);
        if (rval) {
          fprintf (stderr, "CClp_setbnd failed\n"); goto CLEANUP;
        }
        rval = CClp_setbnd (lp->lp, i, 'U', setval);
        if (rval) {
          fprintf (stderr, "CClp_setbnd failed\n"); goto CLEANUP;
        }
        rval = CClp_opt (lp->lp, CClp_METHOD_DUAL);
        if (rval != 0 && rval != 2) {
          fprintf (stderr, "CClp_opt failed\n"); goto CLEANUP;
        } else if (rval == 2) {
          rval = CClp_setbnd (lp->lp, i, 'L', 0.0);
          if (rval) {
            fprintf (stderr, "CClp_setbnd failed\n"); goto CLEANUP;
          }
          rval = CClp_setbnd (lp->lp, i, 'U', 1.0);
          if (rval) {
            fprintf (stderr, "CClp_setbnd failed\n"); goto CLEANUP;
          }
          ninfeas++;
        } else {
          nfixed++;
        }
        rval = CClp_objval (lp->lp, &objval);
        if (rval) {
          fprintf (stderr, "CClp_objval failed\n"); goto CLEANUP;
        }
        printf (" dusted, cum time %.2f, objval %.6f\n",
                CCutil_zeit() - szeit, objval);
        fflush (stdout);
      }
    }

    rval = CClp_opt (lp->lp, CClp_METHOD_DUAL);
    if (rval) {
      fprintf (stderr, "CClp_opt failed\n"); goto CLEANUP;
    }
    rval = CClp_objval (lp->lp, &objval);
    if (rval) {
      fprintf (stderr, "CClp_objval failed\n"); goto CLEANUP;
    }
    printf ("polished away %d dusty edges (%d stubborn) in %.2f seconds, objval %.6f\n",
             nfixed, ninfeas, CCutil_zeit() - szeit, objval);
    fflush (stdout);

    rval = CClp_x (lp->lp, x);
    if (rval) {
        fprintf (stderr, "CClp_x failed\n"); goto CLEANUP;
    }

    for (i=0; i<ncols; i++) {
        if (marks[i] == 0) {
            rval = CClp_setbnd (lp->lp, i, 'L', 0.0);
            if (rval) {
                fprintf (stderr, "CClp_setbnd failed\n"); goto CLEANUP;
            }
            rval = CClp_setbnd (lp->lp, i, 'U', 1.0);
            if (rval) {
                fprintf (stderr, "CClp_setbnd failed\n"); goto CLEANUP;
            }
        }
    }
            
    rval = CClp_load_warmstart (lp->lp, warmstart);
    if (rval) {
        fprintf (stderr, "CClp_load_warmstart failed\n"); goto CLEANUP;
    }
    
    rval = CClp_opt (lp->lp, CClp_METHOD_DUAL);
    if (rval == 2) {
        fprintf (stderr, "restored LP infeasible\n"); goto CLEANUP;
    } else if (rval) {
        fprintf (stderr, "CClp_opt failed\n"); goto CLEANUP;
    }
    
    *newx    = x;
    x = (double *) NULL;

    *newlist = CC_SAFE_MALLOC (2 * ncols, int);
    if ((*newlist) == (int *) NULL) {
        fprintf (stderr, "out of memory in grab_polished_x\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ncols; i++) {
        (*newlist)[2*i] = lp->graph.edges[i].ends[0];
        (*newlist)[2*i+1] = lp->graph.edges[i].ends[1];
    }

    *newcount = ncols;
    rval = 0;

 CLEANUP:
    CClp_free_warmstart (&warmstart);
    
    CC_IFFREE (x, double);
    CC_IFFREE (marks, char);

    if (rval) {
        CC_IFFREE (*newlist, int);
        CC_IFFREE (*newx, double);
    }
    return rval;
}

static int no_tighten (int ncount, int xcount, int *xlist, double *x, int *test,
        double tol)
{
    CC_SRKgraph G;
    int rval;
    int k;

    *test = 0;
    CCcut_SRK_init_graph (&G);

    rval = CCcut_SRK_buildgraph (&G, ncount, xcount, xlist, x);
    if (rval) {
        fprintf (stderr, "CCcut_SRK_buildgraph failed\n");
        goto CLEANUP;
    }
    CCcut_SRK_increment_marker (&G);

    rval = CCcut_SRK_defluff (&G);
    if (rval) {
        fprintf (stderr, "CCcut_SRK_defluff failed\n");
        goto CLEANUP;
    }

    CCcut_SRK_identify_paths_to_edges (&G, &k, 0);

    if (k < (tol * ncount)) {
        *test = 1;
    }


CLEANUP:

    CCcut_SRK_free_graph (&G);
    return rval;
}


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
/*            INTERFACE ROUTINES TO THE EXACT TSP SOLVER                    */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: April 23, 1996                                                    */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCtsp_solve_sparse (int ncount, int ecount, int *elist,             */
/*      int *elen, int *in_tour, int *out_tour, double *in_val,             */
/*      double *optval, int *success, int *foundtour, char *name,           */
/*      double *timebound, int *hit_timebound, int silent,                  */
/*      CCrandstate *rstate)                                                */
/*    SOLVES the TSP over the graph specfied in the edgelist.               */
/*     -elist is an array giving the ends of the edges (in pairs)           */
/*     -elen is an array giving the weights of the edges.                   */
/*     -in_tour gives a starting tour in node node node format (it can      */
/*      be NULL)                                                            */
/*     -out_tour will return the optimal tour (it can be NULL, if it is     */
/*      not NULL then it should point to an array of length at least        */
/*      ncount.                                                             */
/*     -in_val can be used to specify an initial upperbound (it can be      */
/*      NULL)                                                               */
/*     -optval will return the value of the optimal tour.                   */
/*     -success will be set to 1 if the run finished normally, and set to   */
/*      if the search was terminated early (by hitting some predefined      */
/*      limit)                                                              */
/*     -foundtour will be set to 1 if a tour has been found (if success     */
/*      is 0, then it may not be the optimal tour)                          */
/*     -name specifes a char string that will be used to name various       */
/*      files that are written during the branch and bound search (if it    */
/*      is NULL, then "noname" will be used - this will cause problems      */
/*      in a multithreaded program, so specify a distinct name in that      */
/*      case).                                                              */
/*     -silent will suppress most output if set to a nonzero value.         */
/*                                                                          */
/*  int CCtsp_solve_dat (int ncount, CCdatagroup *indat, int *in_tour,      */
/*      int *out_tour, double *in_val, double *optval, int *success,        */
/*      int *foundtour, char *name, double *timebound, int *hit_timebound,  */
/*      int silent, CCrandstate *rstate)                                    */
/*    SOLVES the TSP over the graph specified in the datagroup.             */
/*    LIKE CCtsp_solve_sparse.                                              */
/*                                                                          */
/****************************************************************************/

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


static void
     perm_bound (double *bound, int ncount, CCdatagroup *dat);

static int
     tsp_solve_lp (CCtsp_lp *lp, CCtsp_cutselect *sel, int *out_tour,
         double *optval, int *success, int *foundtour, double *timebound,
         int *hit_timebound, int silent, CCrandstate *rstate),
     find_good_tour (int ncount, CCdatagroup *dat, int *tour,
         double *val, int trials, CCrandstate *rstate),
     build_edges (int ncount, CCdatagroup *dat, int *ecount, int **elist,
         int **elen, CCrandstate *rstate),
     build_extra_edges (int ncount, CCdatagroup *dat, int *ecount,
         int **elist, int **elen, int *valid),
     grab_plan_edges (int ncount, CCdatagroup *dat, CCedgegengroup *plan,
         int *ecount, int **elist, int **elen, CCrandstate *rstate);


int CCtsp_solve_sparse (int ncount, int ecount, int *elist, int *elen,
        int *in_tour, int *out_tour, double *in_val, double *optval,
        int *success, int *foundtour, char *name, double *timebound,
        int *hit_timebound, int silent, CCrandstate *rstate)
{
    int rval = 0;
    CCdatagroup dat;

    *success = 0;
    *foundtour = 0;

    CCutil_init_datagroup (&dat);

    rval = CCutil_graph2dat_sparse (ncount, ecount, elist, elen, 0, &dat);
    if (rval) {
        fprintf (stderr, "CCutil_graph2dat_sparse failed\n"); goto CLEANUP;
    }

    rval = CCtsp_solve_dat (ncount, &dat, in_tour, out_tour, in_val, optval,
                            success, foundtour, name, timebound, hit_timebound,
                            silent, rstate);
    if (rval) {
        fprintf (stderr, "CCtsp_solve_dat failed\n"); goto CLEANUP;
    }

CLEANUP:

    CCutil_freedatagroup (&dat);
    return rval;
}

int CCtsp_solve_dat (int ncount, CCdatagroup *indat, int *in_tour,
        int *out_tour, double *in_val, double *optval, int *success,
        int *foundtour, char *name, double *timebound, int *hit_timebound,
        int silent, CCrandstate *rstate)
{
    int i, norm, newtour, rval = 0;
    CCdatagroup dat;
    CCtsp_cutselect sel;
    int *itour   = (int *) NULL;
    int *otour   = (int *) NULL;
    int *mytour;
    int iecount;
    int *ielist  = (int *) NULL;
    int *ielen   = (int *) NULL;
    int iexcount = 0;
    int iexvalid = 0;
    int *iexlist = (int *) NULL;
    int *iexlen  = (int *) NULL;
    double val, upperbound, tbound, *mytbound;
    double szeit = CCutil_zeit ();
    CCtsp_lpcuts *pool = (CCtsp_lpcuts *) NULL;
    CCtsp_lpcuts *dominopool = (CCtsp_lpcuts *) NULL;
    CCtsp_lp *lp = (CCtsp_lp *) NULL;
    char pname[1024];
  
    printf ("CCtsp_solve_dat ...\n"); fflush (stdout);

    CCutil_init_datagroup (&dat);
    CCtsp_init_cutselect (&sel);

    if (name) {
        sprintf (pname, "%s", name);
    } else {
        sprintf (pname, "noname");
    }

    rval = CCutil_copy_datagroup (ncount, indat, &dat);
    if (rval) {
        fprintf (stderr, "CCutil_copy_datagroup failed\n");
        goto CLEANUP;
    }

    if (!in_tour) {
        itour = CC_SAFE_MALLOC (ncount, int);
        if (!itour) {
            fprintf (stderr, "out of memory in CCtsp_solve_dat\n");
            rval = 1; goto CLEANUP;
        }
        mytour = itour;

        if (!in_val) {
            rval = find_good_tour (ncount, &dat, mytour, &val, 1, rstate);
        } else {
            printf ("Initial bnd %f given - use short tour run\n", *in_val);
            fflush (stdout);
            rval = find_good_tour (ncount, &dat, mytour, &val, 0, rstate);
        }
        if (rval) {
            fprintf (stderr, "find_good_tour failed\n"); goto CLEANUP;
        }
    } else {
        mytour = in_tour;
    }

    rval = CCutil_datagroup_perm (ncount, &dat, mytour);
    if (rval) {
        fprintf (stderr, "CCutil_datagroup_perm failed\n"); goto CLEANUP;
    }

    otour = CC_SAFE_MALLOC (ncount, int);
    if (!otour) {
        fprintf (stderr, "out of memory in CCtsp_solve_dat\n");
        rval = 1; goto CLEANUP;
    }

    if (in_val) upperbound = *in_val;
    else        upperbound = CCtsp_LP_MAXDOUBLE;

    perm_bound (&val, ncount, &dat);
    if (val <= upperbound) {
        printf ("Set initial upperbound to %.0f (from tour)\n", val);
        fflush (stdout);
        upperbound = val;
        *foundtour = 1;
        for (i = 0; i < ncount; i++) {
            otour[i] = i;
        }
    }

    rval = build_edges (ncount, &dat, &iecount, &ielist, &ielen, rstate);
    if (rval) {
        fprintf (stderr, "build_edges failed\n"); goto CLEANUP;
    }

    rval = build_extra_edges (ncount, &dat, &iexcount, &iexlist, &iexlen,
                              &iexvalid);
    if (rval) {
        fprintf (stderr, "build_extra_edges failed\n"); goto CLEANUP;
    }

    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &pool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");
    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &dominopool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed for dominos");

    rval = CCtsp_init_lp (&lp, pname, -1, (char *) NULL, ncount, &dat,
               iecount, ielist, ielen, iexcount, iexlist, iexlen, iexvalid,
               mytour, upperbound, pool, dominopool, silent, rstate);

    if (rval == 2) {
        int is_infeasible;

        printf ("CCtsp_init_lp reports an infeasible LP\n");
        fflush (stdout);
        rval = CCtsp_verify_infeasible_lp (lp, &is_infeasible, silent);
        if (rval) {
            fprintf (stderr, "CCtsp_verify_infeasible_lp failed\n");
            goto CLEANUP;
        }
        if (!is_infeasible) {
            printf ("Couldn't verify infeasible LP\n"); fflush (stdout);
            rval = 1; goto CLEANUP;
        }
        *optval = CCtsp_LP_MAXDOUBLE;
        *success = 1;
        goto CLEANUP;
    } else if (rval) {
        fprintf (stderr, "CCtsp_init_lp failed\n"); goto CLEANUP;
    }

    if (timebound) {
        tbound = (*timebound) - (CCutil_zeit () - szeit);
        if (tbound <= 0.0) {
            printf ("Hit time bound\n"); fflush (stdout);
            if (hit_timebound) *hit_timebound = 1;
            goto DONE;
        }
        mytbound = &tbound;
    } else {
        mytbound = (double *) NULL;
    }

    rval = tsp_solve_lp (lp, &sel, otour, optval, success, &newtour,
                         mytbound, hit_timebound, silent, rstate);
    if (rval) {
        fprintf (stderr, "tsp_solve_lp failed\n"); goto CLEANUP;
    }
    if (newtour == 1) *foundtour = 1;

    if (*success) {
        int istour;

        CCutil_dat_getnorm (&dat, &norm);
        if ((norm == CC_SPARSE) && (*foundtour)) {
            rval = CCutil_sparse_real_tour (ncount, &dat, otour, &istour);
            if (rval) {
                fprintf (stderr, "CCutil_sparse_real_tour failed\n");
                goto CLEANUP;
            }
            if (istour == 0) {
                printf ("Tour uses artificial edges\n"); fflush (stdout);
                *optval = CCtsp_LP_MAXDOUBLE;
                *foundtour = 0;
            } else {
                printf ("Optimal tour: %.0f\n", *optval); fflush (stdout);
            }
        } else if (*foundtour) {
            printf ("Optimal tour: %.0f\n", *optval); fflush (stdout);
        } else {
            printf ("Did not find a tour\n"); fflush (stdout);
        }
    } else {
        printf ("Did not succeed in finding optimal tour\n"); fflush (stdout);
    }

DONE:

    if ((*foundtour == 1) && (out_tour != (int *) NULL)) {
        for (i = 0; i < ncount; i++) {
            out_tour[i] = mytour[otour[i]];
        }
    }

    printf ("Total Time to solve TSP: %.2f\n", CCutil_zeit () - szeit);
    fflush (stdout);

CLEANUP:

    if (rval == 0) {
        char buf[1024];
        int sval;

        sprintf (buf, "%s.pul", pname);
        sval = CCutil_sdelete_file (buf);
        sval = CCutil_sdelete_file_backup (buf);
        sprintf (buf, "%s.sav", pname);
        sval = CCutil_sdelete_file (buf);
        sval = CCutil_sdelete_file_backup (buf);
    }

    CCtsp_free_tsp_lp_struct (&lp);
    if (pool) { CCtsp_free_cutpool (&pool); }
    if (dominopool) { CCtsp_free_cutpool (&dominopool); }
    CCutil_freedatagroup (&dat);
    CC_IFFREE (itour, int);
    CC_IFFREE (otour, int);
    CC_IFFREE (ielist, int);
    CC_IFFREE (ielen, int);

    return rval;
}

static int tsp_solve_lp (CCtsp_lp *lp, CCtsp_cutselect *sel, int *out_tour,
        double *optval, int *success, int *foundtour, double *timebound,
        int *hit_timebound, int silent, CCrandstate *rstate)
{
    int i, rval = 0;
    int ncount = lp->graph.ncount;
    int *tour    = (int *) NULL;
    CCbigguy bupper, exactbound;
    double tbound, *mytbound;
    double szeit = CCutil_zeit ();

    *success = 0;
    *optval = CCtsp_LP_MAXDOUBLE;
    *foundtour = 0;

    if (!lp->dat) {
        fprintf (stderr, "tsp_solve_lp called without a datagroup\n");
        rval = 1; goto CLEANUP;
    }

    tour = CC_SAFE_MALLOC (ncount, int);
    if (!tour) {
        fprintf (stderr, "out of memory in tsp_solve_lp failed\n");
        rval = 1;  goto CLEANUP;
    }

    rval = CCtsp_cutselect_set_tols (sel, lp, 1, silent);
    if (rval) {
        fprintf (stderr, "CCtsp_cutselect_set_tols failed\n");
        goto CLEANUP;
    }

    rval = CCtsp_cutting_loop (lp, sel, 1, silent, rstate);
    if (rval == 2) {
        int is_infeasible;

        printf ("CCtsp_cutting_loop reports an infeasible LP\n");
        fflush (stdout);
        rval = CCtsp_verify_infeasible_lp (lp, &is_infeasible, silent);
        if (rval) {
            fprintf (stderr, "CCtsp_verify_infeasible_lp failed\n");
            goto CLEANUP;
        }
        if (!is_infeasible) {
            printf ("Couldn't verify infeasibile LP\n"); fflush (stdout);
            rval = 1; goto CLEANUP;
        }
        *optval = CCtsp_LP_MAXDOUBLE;
        *success = 1;
        printf ("Final LP has %d rows, %d columns, %d nonzeros\n",
                CClp_nrows (lp->lp), CClp_ncols (lp->lp),
                CClp_nnonzeros (lp->lp));
        fflush (stdout);
        goto CLEANUP;
    } else if (rval) {
        fprintf (stderr, "Ctsp_cutting_loop failed\n"); goto CLEANUP;
    }

    {
        double tourval;
        rval = CCtsp_call_x_heuristic (lp, &tourval, tour, silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_call_x_heuristic failed\n");
            goto CLEANUP;
        }
        if (tourval < lp->upperbound) {
            printf ("Upperbound from x-heuristic: %.2f\n", tourval);
            lp->upperbound = tourval;
            *foundtour = 1;
            if (out_tour) {
                for (i = 0; i < ncount; i++) {
                    out_tour[i] = tour[i];
                }
            }
        }
    }

    if (ncount >= 100000) {
        fprintf (stderr, "Not running branching on problems of this size\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_exact_price (lp, &exactbound, 0, 0, silent);
    if (rval) {
        fprintf (stderr, "CCtsp_exact_price failed\n"); goto CLEANUP;
    }
    lp->exact_lowerbound = exactbound;
    printf ("Exact lower bound: %.6f\n", CCbigguy_bigguytod (exactbound));
    printf ("DIFF: %f\n", lp->lowerbound - CCbigguy_bigguytod (exactbound));
    fflush (stdout);

    bupper = CCbigguy_dtobigguy (lp->upperbound);
    CCbigguy_sub (&bupper, CCbigguy_ONE);

    if (CCbigguy_cmp (lp->exact_lowerbound, bupper) > 0) {
        printf ("Established Bound: %.0f\n", lp->upperbound); fflush (stdout);
        *optval = lp->upperbound;
        *success = 1;
        goto CLEANUP;
    }

    rval = CCtsp_eliminate_variables (lp, 1, silent);
    if (rval) {
        fprintf (stderr, "CCtsp_eliminate_variables failed\n"); goto CLEANUP;
    }

    if (timebound) {
        tbound = (*timebound) - (CCutil_zeit () - szeit);
        if (tbound <= 0.0) {
            printf ("Hit time bound\n"); fflush (stdout);
            if (hit_timebound) *hit_timebound = 1;
            goto CLEANUP;
        }
        mytbound = &tbound;
    } else {
        mytbound = (double *) NULL;
    }

    {
        int saveproof            = 0;
        int bbcount              = 0;
        int usebranchcliques     = 1;
        int tentative_branch_num = 0;
        unsigned short hostport  = (unsigned short) 0;
        double branchzeit        = 0.0;
        double upbound           = lp->upperbound;

        bbcount = 0;

        rval = CCtsp_write_probroot_id (lp->probloc, lp);
        if (rval) {
            fprintf (stderr, "CCtsp_write_probroot_id failed\n");
            goto CLEANUP;
        }

        rval = CCtsp_bfs_brancher (lp->probloc, lp->id, lp->lowerbound, sel,
                sel, &upbound, &bbcount, usebranchcliques,  lp->dat,
                lp->perm, lp->pool, ncount, tour, hostport, &branchzeit,
                saveproof, tentative_branch_num, 0, mytbound, hit_timebound,
                silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_bfs_brancher failed\n"); goto CLEANUP;
        }

        if (upbound < lp->upperbound && out_tour) {
            *foundtour = 1;
            for (i = 0; i < ncount; i++) {
                out_tour[i] = tour[i];
            }
        }
        printf ("Total number of nodes in search tree: %d\n", bbcount);
        fflush (stdout);
        *optval = upbound;
        *success = 1;
    }
 
CLEANUP:

    CC_IFFREE (tour, int);
    return rval;
}

static int find_good_tour (int ncount, CCdatagroup *dat, int *tour,
        double *tval, int trials, CCrandstate *rstate)
{
    int rval = 0;
    CCedgegengroup plan;
    int ecount;
    int *elist = (int *) NULL;
    int tcount;
    int *tlist = (int *) NULL;
    int *bestcyc = (int *) NULL;
    int *cyc     = (int *) NULL;
    int *tmp;
    double val, bestval, szeit;
    int i, kicks, istour;

    szeit = CCutil_zeit ();
    bestval = CCtsp_LP_MAXDOUBLE;

    if (trials == -1) {
        kicks = (ncount > 400 ? 100 : ncount/4);
    } else {
        kicks = (ncount > 1000 ? 500 : ncount/2);
    }

    printf ("Finding a good tour for compression ...\n"); fflush (stdout);

    cyc     = CC_SAFE_MALLOC (ncount, int);
    bestcyc = CC_SAFE_MALLOC (ncount, int);
    if (!cyc || !bestcyc) {
        fprintf (stderr, "out of memory in find_good_tour\n");
        rval = 1; goto CLEANUP;
    }

    CCedgegen_init_edgegengroup (&plan);
    plan.quadnearest = 2;
    rval = CCedgegen_edges (&plan, ncount, dat, (double *) NULL, &ecount,
                            &elist, 1, rstate);
    if (rval) {
        fprintf (stderr, "CCedgegen_edges failed\n"); goto CLEANUP;
    }
    plan.quadnearest = 0;


    plan.tour.greedy = 1;
    rval = CCedgegen_edges (&plan, ncount, dat, (double *) NULL, &tcount,
                            &tlist, 1, rstate);
    if (rval) {
        fprintf (stderr, "CCedgegen_edges failed\n"); goto CLEANUP;
    }

    if (tcount != ncount) {
        fprintf (stderr, "wrong edgeset from CCedgegen_edges\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_edge_to_cycle (ncount, tlist, &istour, cyc);
    if (rval) {
        fprintf (stderr, "CCutil_edge_to_cycle failed\n");
        rval = 1; goto CLEANUP;
    }
    if (!istour) {
        fprintf (stderr, "Starting tour has an error\n");
        rval = 1; goto CLEANUP;
    }
    CC_FREE (tlist, int);

    rval = CClinkern_tour (ncount, dat, ecount, elist, ncount, kicks,
                    cyc, bestcyc, &bestval, 0, 0.0, 0.0, (char *) NULL,
                    CC_LK_GEOMETRIC_KICK, rstate);
    if (rval) {
        fprintf (stderr, "CClinkern_tour failed\n"); goto CLEANUP;
    }
    printf ("LK Initial Run: %.1f\n", bestval); fflush (stdout);

    for (i = 0; i < trials; i++) {
        rval = CClinkern_tour (ncount, dat, ecount, elist, ncount, kicks,
                        (int *) NULL, cyc, &val, 1, 0.0, 0.0, (char *) NULL,
                        CC_LK_GEOMETRIC_KICK, rstate);
        if (rval) {
            fprintf (stderr, "CClinkern_tour failed\n"); goto CLEANUP;
        }
        printf ("LK Run %d: %.1f\n", i, val); fflush (stdout);
        if (val < bestval) {
            CC_SWAP (cyc, bestcyc, tmp);
            bestval = val;
        }
    }

    if (trials > 0) {
        rval = CClinkern_tour (ncount, dat, ecount, elist, ncount, 2 * kicks,
                        bestcyc, tour, tval, 1, 0.0, 0.0, (char *) NULL,
                        CC_LK_GEOMETRIC_KICK, rstate);
        if (rval) {
            fprintf (stderr, "CClinkern_tour failed\n"); goto CLEANUP;
        }
        printf ("LK Run from best tour: %.1f\n", *tval); fflush (stdout);
    } else {
        for (i = 0; i < ncount; i++) {
            tour[i] = bestcyc[i];
        }
        *tval = bestval;
    }

    printf ("Time to find compression tour: %.2f (seconds)\n",
            CCutil_zeit() - szeit);
    fflush (stdout);

CLEANUP:

    CC_IFFREE (cyc, int);
    CC_IFFREE (bestcyc, int);
    CC_IFFREE (elist, int);
    CC_IFFREE (tlist, int);
    return rval;
}

static void perm_bound (double *bound, int ncount, CCdatagroup *dat)
{
    double bnd;
    int i;

    bnd = CCutil_dat_edgelen (ncount - 1, 0, dat);
    for (i = 1; i < ncount; i++) {
        bnd += CCutil_dat_edgelen (i-1, i, dat);
    }
    *bound = bnd;
}

static int build_edges (int ncount, CCdatagroup *dat, int *ecount, int **elist,
        int **elen, CCrandstate *rstate)
{
    int rval = 0;
    CCedgegengroup plan;
    int norm;

    *ecount = 0;
    *elist  = (int *) NULL;
    *elen   = (int *) NULL;

    CCutil_dat_getnorm (dat, &norm);

    if (norm == CC_SPARSE) {
        if (dat->sparse_ecount <=  2 * ncount) {
            printf ("Use entire sparse graph as initial edge set\n");
            fflush (stdout);
            rval = CCutil_get_sparse_dat_edges (ncount, dat, ecount, elist,
                                                elen);
            if (rval) {
                fprintf (stderr, "CCutil_get_sparse_dat_edges failed\n");
                goto CLEANUP;
            }
        } else {
            int tecount;
            int *telist = (int *) NULL;
            int *telen  = (int *) NULL;

            CCedgegen_init_edgegengroup (&plan);
            plan.nearest = 4;
            rval = grab_plan_edges (ncount, dat, &plan, &tecount, &telist,
                                    &telen, rstate);
            if (rval) {
                fprintf (stderr, "grab_plan_edges failed\n"); goto CLEANUP;
            }

            rval = CCutil_sparse_strip_edges (dat, tecount, telist, telen,
                                              ecount, elist, elen);
            if (rval) {
                fprintf (stderr, "CCutil_sparse_strip_edges failed\n");
                CC_IFFREE (telist, int);
                CC_IFFREE (telen, int);
                goto CLEANUP;
            }

            CC_IFFREE (telist, int);
            CC_IFFREE (telen, int);
        }
    } else {
        CCedgegen_init_edgegengroup (&plan);
        plan.linkern.count = 10;
        plan.linkern.quadnearest = 2;
        plan.linkern.greedy_start = 0;
        plan.linkern.nkicks = (ncount / 100) + 1;

        rval = grab_plan_edges (ncount, dat, &plan, ecount, elist, elen,
                                rstate);
        if (rval) {
            fprintf (stderr, "grab_plan_edges failed\n"); goto CLEANUP;
        }
    }

CLEANUP:

    if (rval) {
        CC_IFFREE (*elist, int);
        CC_IFFREE (*elen, int);
    }

    return rval;
}

static int  build_extra_edges (int ncount, CCdatagroup *dat, int *ecount,
        int **elist, int **elen, int *valid)
{
    int norm, rval = 0;

    *ecount = 0;
    *elist  = (int *) NULL;
    *elen   = (int *) NULL;
    *valid  = 0;

    CCutil_dat_getnorm (dat, &norm);

    if (norm == CC_SPARSE) {
        rval = CCutil_get_sparse_dat_edges (ncount, dat, ecount, elist, elen);
        if (rval) {
            fprintf (stderr, "CCutil_get_sparse_dat_edges failed\n");
            goto CLEANUP;
        }
        *valid = 1;
    }

CLEANUP:

    if (rval) {
        CC_IFFREE (*elist, int);
        CC_IFFREE (*elen, int);
    }

    return rval;

}

static int grab_plan_edges (int ncount, CCdatagroup *dat,
    CCedgegengroup *plan, int *ecount, int **elist, int **elen,
    CCrandstate *rstate)
{
    int i, rval = 0;

    rval = CCedgegen_edges (plan, ncount, dat, (double *) NULL, ecount,
                            elist, 1, rstate);
    if (rval) {
        fprintf (stderr, "CCedgegen_edges failed\n"); goto CLEANUP;
    }

    *elen = CC_SAFE_MALLOC (*ecount, int);
    if (!(*elen)) {
        fprintf (stderr, "out of memory in grab_plan_edges\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < *ecount; i++) {
        (*elen)[i] = CCutil_dat_edgelen ((*elist)[2*i],
                                         (*elist)[(2*i) + 1], dat);
    }

CLEANUP:

    return rval;
} 



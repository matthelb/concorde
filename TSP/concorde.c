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
/*                  THE MAIN PROGRAM FOR CONCORDE                           */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: September 25, 1995                                                */
/*        November 28, 2003 (bico)                                          */
/*  SEE short decsribtion in usage ().                                      */
/*                                                                          */
/*  NOTE:  When using CC_SPARSE edge sets, it is important to specify a     */
/*   a tour or an upperbound if you have one, since our heuristics are      */
/*   not designed for finding Hamilton circuits in sparse graphs.           */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "edgegen.h"
#include "tsp.h"
#include "linkern.h"
#include "heldkarp.h"
#include "bigguy.h"
#include "macrorus.h"

#define CC_JUST_SUBTOUR (1)
#define CC_JUST_BLOSSOM (2)
#define CC_JUST_SUBTOUR_AND_BLOSSOM (3)
#define CC_JUST_FAST_CUTS (4)

static int norm = CC_EUCLIDEAN;
static char *datfname        = (char *) NULL;
static char *edgegenfname    = (char *) NULL;
static char *problname       = (char *) NULL;
static char *probfname       = (char *) NULL;
static char *edgefname       = (char *) NULL;
static char *fullfname       = (char *) NULL;
static char *tourfname       = (char *) NULL;
static char *masterfname     = (char *) NULL;
static char *poolfname       = (char *) NULL;
#ifdef CCtsp_USE_DOMINO_CUTS
static char *dominopoolfname = (char *) NULL;
#endif
static char *restartfname    = (char *) NULL;
static char *xfname          = (char *) NULL;
static char *outfname        = (char *) NULL;
static char *filecutname     = (char *) NULL;
static int seed                 = 0;
static int nnodes_want          = 0;
static int binary_in            = 0;
static int tsplib_in            = 1;
static int gridsize             = 0;
static int just_cuts            = 0;
static int dontcutroot          = 0;
static int usetighten           = 0;
static int usedominos           = 0;
static int maxchunksize         = 16;
static int multiple_chunker     = 0;
static int valid_edges          = 0;
static int dfs_branching        = 0;
static int bfs_branching        = 1;
static int simple_branching     = 0;
static int usebranchcliques     = 1;  
static int tentative_branch_num = 0;
static int complete_price       = 0;
static int want_rcnearest       = 0;
static int output_tour_as_edges = 0;
static int run_silently         = 1;
static int be_nethost           = 0;
static int unlink_files         = 0;
static double initial_ub = CCtsp_LP_MAXDOUBLE;
static unsigned short hostport = CCtsp_HOST_PORT;
static char *grunthostname = (char *) NULL;
static char *cutbossname   = (char *) NULL;
static char *dombossname   = (char *) NULL;
static int eliminate_edges = -1;   /* Set to 1 to force elim, 0 to not elim */
static int eliminate_sparse = 0;   /* Set to 1 to elim from full edge list  */
                                   /* if the full edge list is valid        */
static int longedge_branching = 1; /* Set to 0 to turn off           */
static int save_proof = 0;         /* Set to 1 to save the proof     */
static int standalone_branch = 0;  /* Set to 1 to do a manual branch */


static void
    adjust_upbound (double *bound, int ncount, CCdatagroup *dat),
    usage (char *f);

static int
    handle_just_cuts (CCtsp_lp *lp, int the_cuts, CCrandstate *rstate,
       int silent),
    run_hk (int ncount, CCdatagroup *dat, int *hk_tour),
    build_edges (int *p_ecount, int **p_elist, int **p_elen,
        int ncount, int *ptour, CCdatagroup *dat, char *in_edgefname,
        char *in_edgegenfname, int in_just_cuts, int silent,
        CCrandstate *rstate),
    build_fulledges (int *p_excount, int **p_exlist, int **p_exlen,
        int ncount, int *ptour, char *in_fullfname),
    parseargs (int ac, char **av),
    find_tour (int ncount, CCdatagroup *dat, int *perm, double *ub,
            int trials, int silent, CCrandstate *rstate),
    getedges (CCdatagroup *dat, CCedgegengroup *plan, int ncount, int *ecount,
            int **elist, int **elen, int silent, CCrandstate *rstate),
    dump_rc (CCtsp_lp *lp, int count, char *pname, int usesparse);



int main (int ac, char **av)
{
    int rval = 0;
    char *probname = (char *) NULL;
    int i, ncount, silent, allow_dups, use_gridsize;
    CCdatagroup dat;
    CCtsp_lp *lp = (CCtsp_lp *) NULL;
    CCtsp_cutselect sel, tentativesel;
    CCtsp_lpcuts *pool = (CCtsp_lpcuts *) NULL;
    CCtsp_lpcuts *dominopool = (CCtsp_lpcuts *) NULL;
    int *ptour = (int *) NULL;
    int ecount = 0;
    int *elist = (int *) NULL;
    int *elen = (int *) NULL;
    int excount = 0;
    int *exlist = (int *) NULL;
    int *exlen = (int *) NULL;
    int *besttour = (int *) NULL;
    int is_infeasible = 0;
    int bbcount = 0;
    double szeit;
    double upbound = 0.0;
    double branchzeit = 0.0;
    CCrandstate rstate;
    char buf[1024];

    CCutil_init_datagroup (&dat);

    rval = CCutil_print_command (ac, av);
    CCcheck_rval (rval, "CCutil_print_command failed");
    
    szeit = CCutil_zeit ();
    seed = (int) CCutil_real_zeit ();

    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    CCutil_printlabel ();
    CCutil_signal_init ();
    CCutil_sprand (seed, &rstate);
    printf ("Using random seed %d\n", seed); fflush (stdout);

    silent = run_silently;

    if (grunthostname) {
#ifdef CC_NETREADY
        rval = CCtsp_grunt (grunthostname, hostport, poolfname, cutbossname,
                            problname, silent, &rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_grunt failed\n");
        }
        goto CLEANUP;
#else  /* CC_NETREADY */
        fprintf (stderr, "Networking not enabled\n");
        rval = 1; goto CLEANUP;
#endif /* CC_NETREADY */
    }

    if (!be_nethost) {
        hostport = 0;
    }

    CCtsp_init_cutselect (&sel);
    CCtsp_init_tentative_cutselect (&tentativesel);
    CCtsp_cutselect_tighten (&sel, usetighten);
    CCtsp_cutselect_tighten (&tentativesel, usetighten);
    CCtsp_cutselect_chunksize (&sel, maxchunksize);
    CCtsp_cutselect_dominos (&sel, usedominos);
    if (filecutname) CCtsp_cutselect_filecuts (&sel, filecutname);
#ifdef CC_NETREADY
    if (cutbossname != (char *) NULL) {
        CCtsp_cutselect_remotepool (&sel, cutbossname);
        CCtsp_cutselect_remotepool (&tentativesel, cutbossname);
    }
    if (dombossname != (char *) NULL) {
        CCtsp_cutselect_domboss (&sel, dombossname);
        CCtsp_cutselect_domboss (&tentativesel, dombossname);
    }
#endif /* CC_NETREADY */

    if (problname)        probname = CCtsp_problabel (problname);
    else if (datfname)    probname = CCtsp_problabel (datfname);
    else if (masterfname) probname = CCtsp_problabel (masterfname);
    else                  probname = CCtsp_problabel ("unnamed");
    if (probname == (char *) NULL) {
        fprintf (stderr, "CCtsp_problabel failed\n");
        rval = 1; goto CLEANUP;
    }
    if (problname == (char *) NULL) {
        problname = probname;
    }

    if (masterfname) {
        rval = CCutil_getmaster (masterfname, &ncount, &dat, &ptour);
        CCcheck_rval (rval, "CCutil_getmaster failed");
        if (ncount < 10) {
            fprintf (stderr, "Master file has less than 10 nodes - Abort.\n");
            rval = 1; goto CLEANUP;
        }
    }  else {
        CCutil_init_datagroup (&dat);
        if (tsplib_in && datfname != (char *) NULL) {
            rval = CCutil_gettsplib (datfname, &ncount, &dat);
            CCcheck_rval (rval, "CCutil_gettsplib failed");
        } else {
            ncount = nnodes_want;
            if (gridsize < 0) {
                use_gridsize = -gridsize;
                allow_dups = 0;
            } else if (gridsize > 0) {
                use_gridsize = gridsize;
                allow_dups = 1;
            } else {
                use_gridsize = nnodes_want;
                allow_dups = 0;
            }
            rval = CCutil_getdata (datfname, binary_in, norm, &ncount, &dat,
                                   use_gridsize, allow_dups, &rstate);
            CCcheck_rval (rval, "CCutil_getdata failed");
        }

        /* Handle small instances */

        if (ncount < 3) {
            printf ("Only %d nodes -- must have at least 3 nodes in a TSP\n",
                     ncount);
            fflush (stdout);
            rval = 1;  goto CLEANUP;
        } else if (ncount < 10) {
            besttour = CC_SAFE_MALLOC (ncount, int);
            CCcheck_NULL (besttour, "out of memory for besttour");
            if (ncount == 3) {
                for (i = 0; i < ncount; i++) besttour[i] = i;
            } else {
                rval = run_hk (ncount, &dat, besttour);
                CCcheck_rval (rval, "run_hk failed");
            }
            ptour = CC_SAFE_MALLOC (ncount, int);
            CCcheck_NULL (ptour, "out of memory for ptour");
            for (i = 0; i < ncount; i++) ptour[i] = i;
            rval = CCtsp_dumptour (ncount, &dat, ptour, probname, besttour,
                                   outfname, output_tour_as_edges, silent);
            CCcheck_rval (rval, "CCtsp_dumptour failed");

            printf ("Total Running Time: %.2f (seconds)\n",
                     CCutil_zeit () - szeit);
            fflush (stdout);
            goto CLEANUP;
        }

        /***** Get the permutation tour and permute the data  *****/

        ptour = CC_SAFE_MALLOC (ncount, int);
        CCcheck_NULL (ptour, "out of memory for ptour");

        if (tourfname) {
            rval = CCutil_getcycle (ncount, tourfname, ptour, 0);
            CCcheck_rval (rval, "CCutil_getcycle failed");
        } else {
            double bnd;
            if (just_cuts > 0) {
                rval = find_tour (ncount, &dat, ptour, &bnd, -1, silent,
                                  &rstate);
            } else if (initial_ub == CCtsp_LP_MAXDOUBLE) {
                rval = find_tour (ncount, &dat, ptour, &bnd, 1, silent,
                                  &rstate);
            } else {
                if (!silent) {
                    printf ("Initial bnd %f - use short LK\n", initial_ub);
                    fflush (stdout);
                }
                rval = find_tour (ncount, &dat, ptour, &bnd, 0, silent,
                                 &rstate);
            }
            CCcheck_rval (rval, "find_tour failed");
        }

        rval = CCutil_datagroup_perm (ncount, &dat, ptour);
        CCcheck_rval (rval, "CCutil_datagroup_perm failed");

        sprintf (buf, "%s.mas", probname);
        rval = CCutil_putmaster (buf, ncount, &dat, ptour);
        CCcheck_rval (rval, "CCutil_putmaster failed");
    }

    adjust_upbound (&initial_ub, ncount, &dat);

    if (!probfname && !restartfname) {
        rval = build_edges (&ecount, &elist, &elen, ncount, ptour,
                            &dat, edgefname, edgegenfname, just_cuts,
                            silent, &rstate);
        CCcheck_rval (rval, "build_edges failed");
    }

    rval = build_fulledges (&excount, &exlist, &exlen, ncount, ptour,
                            fullfname);
    CCcheck_rval (rval, "build_fulledges failed");
    
    rval = CCtsp_init_cutpool (&ncount, poolfname, &pool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");
#ifdef CCtsp_USE_DOMINO_CUTS
    rval = CCtsp_init_cutpool (&ncount, dominopoolfname, &dominopool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed for dominos");
#endif

    /***** Initialize besttour to be the permutation tour  ****/

    besttour = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (besttour, "out of memory for besttour");
    for (i = 0; i < ncount; i++) {
        besttour[i] = i;
    }

    if (restartfname) {
        upbound  = initial_ub;
        bbcount = 0;

        rval = CCtsp_bfs_restart (problname, restartfname, &sel,
                &tentativesel, &upbound, &bbcount, usebranchcliques, &dat,
                ptour, pool, ncount, besttour, hostport, &branchzeit,
                save_proof, tentative_branch_num, longedge_branching,
                (double *) NULL, (int *) NULL, silent, &rstate);
        CCcheck_rval (rval, "CCtsp_bfs_restart failed");
        goto DONE;
    }

    rval = CCtsp_dumptour (ncount, &dat, ptour, probname, besttour,
                           (char *) NULL, 0, silent);
    CCcheck_rval (rval, "CCtsp_dumptour failded");

    rval = CCtsp_init_lp (&lp, problname, -1, probfname, ncount, &dat,
                    ecount, elist, elen, excount, exlist, exlen, valid_edges,
                    ptour, initial_ub, pool, dominopool, silent, &rstate);
    if (rval == 2) {
        printf ("CCtsp_init_lp reports an infeasible LP\n");
        rval = CCtsp_verify_infeasible_lp (lp, &is_infeasible, silent);
        CCcheck_rval (rval, "CCtsp_verify_infeasible_lp failed");
        if (!is_infeasible) {
            printf ("Couldn't verify infeasible LP\n"); fflush (stdout);
            rval = 1; goto CLEANUP;
        }
        upbound = CCtsp_LP_MAXDOUBLE;
        bbcount = 1;
        goto DONE;
    } else if (rval) {
        fprintf (stderr, "CCtsp_init_lp failed\n"); goto CLEANUP;
    }

    CCutil_start_timer (&lp->stats.total);
    
    ecount = 0;
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    excount = 0;
    CC_IFFREE (exlist, int);
    CC_IFFREE (exlen, int);

    if (0 && lp->full_edges_valid) {
        if (CCtsp_inspect_full_edges (lp)) {
            fprintf (stderr, "full edge set does not contain all LP edges\n");
            rval = 1; goto CLEANUP;
        }
    }

    if (standalone_branch) {
        rval = CCtsp_do_interactive_branch (lp, silent, &rstate);
        CCcheck_rval (rval, "CCtsp_do_interactive_branch failed");
        printf ("Total Running Time: %.2f (seconds)\n", CCutil_zeit () - szeit);
        goto CLEANUP;
    }

    if (just_cuts > 0) {
        rval = handle_just_cuts (lp, just_cuts, &rstate, silent);
        CCcheck_rval (rval, "handle_just_cuts failed");
        if (want_rcnearest) {
            rval = dump_rc (lp, want_rcnearest, probname, 0);
            CCcheck_rval (rval, "dump_rc failed");
        }
        if (xfname) {
            rval = CCtsp_dump_x (lp, xfname);
            CCcheck_rval (rval, "CCtsp_dump_x failed");
        }
        goto DONE;
    }

    rval = CCtsp_cutselect_set_tols (&sel, lp, 1, silent);
    CCcheck_rval (rval, "CCtsp_cutselect_set_tols failed");

    if (dontcutroot == 0) {
        if (multiple_chunker) {
            rval = CCtsp_cutting_multiple_loop (lp, &sel, 1, maxchunksize,
                                    1, silent, &rstate);
        } else {
            rval = CCtsp_cutting_loop (lp, &sel, 1, silent, &rstate);
        }
        if (rval == 2) {
            printf ("CCtsp_cutting_loop reports an infeasible LP\n");
            rval = CCtsp_verify_infeasible_lp (lp, &is_infeasible, silent);
            CCcheck_rval (rval, "CCtsp_verify_infeasible_lp failed");
            if (!is_infeasible) {
                printf ("Couldn't verify infeasibile LP\n");
                fflush (stdout);
                rval = 1; goto CLEANUP;
            }
            upbound = CCtsp_LP_MAXDOUBLE;
            bbcount = 1;
            CCutil_stop_timer (&lp->stats.total, 1);
            printf ("Final LP has %d rows, %d columns, %d nonzeros\n",
                    CClp_nrows (lp->lp), CClp_ncols (lp->lp),
                    CClp_nnonzeros (lp->lp));

            goto DONE;
        } else if (rval) {
            fprintf (stderr, "cutting_loop failed\n");
            goto CLEANUP;
        }
    }

    {
        double tourval;
        CCutil_start_timer (&lp->stats.linkern);
        rval = CCtsp_call_x_heuristic (lp, &tourval, besttour, silent, &rstate);
        CCcheck_rval (rval, "CCtsp_call_x_heuristic failed");

        if (!silent) CCutil_stop_timer (&lp->stats.linkern, 1);
        else         CCutil_stop_timer (&lp->stats.linkern, 0);

        if (tourval < lp->upperbound) {
            printf ("New upperbound from x-heuristic: %.2f\n", tourval);
            lp->upperbound = tourval;
            rval = CCtsp_dumptour (ncount, &dat, ptour, probname, besttour,
                                   (char *) NULL, 0, silent);
            CCcheck_rval (rval, "CCtsp_dumptour failed");
        }
        printf ("Final lower bound %f, upper bound %f\n", lp->lowerbound,
                                                          lp->upperbound);
        fflush (stdout);
    }

    if (xfname) {
        rval = CCtsp_dump_x (lp, xfname);
        CCcheck_rval (rval, "CCtsp_dump_x failed");
    }
    if (want_rcnearest) {
        rval = dump_rc (lp, want_rcnearest, probname, 0);
        CCcheck_rval (rval, "dump_rc failed");
    }

    if (lp->graph.ncount < 100000 || complete_price) {
        CCbigguy bound;
        CCbigguy bupper;
        rval = CCtsp_exact_price (lp, &bound, complete_price, 0, silent);
        if (rval) {
            fprintf (stderr, "CCtsp_exact_price failed\n");
            goto CLEANUP;
        }
        lp->exact_lowerbound = bound;
        printf ("Exact lower bound: %.6f\n", CCbigguy_bigguytod (bound));
        if (1 || !silent) {
            printf ("DIFF: %f\n", lp->lowerbound - CCbigguy_bigguytod (bound));
            fflush (stdout);
        }

        bupper = CCbigguy_dtobigguy (lp->upperbound);
        CCbigguy_sub (&bupper, CCbigguy_ONE);

        if (CCbigguy_cmp (lp->exact_lowerbound, bupper) > 0) {
            upbound = lp->upperbound;
            bbcount = 1;
            if (!dfs_branching && !bfs_branching) {
                printf ("Optimal Solution: %.2f\n", upbound);
                printf ("Number of bbnodes: %d\n", bbcount);
                fflush (stdout);
            }
            if (!silent) {
                CCutil_stop_timer (&lp->stats.total, 1);
            } else {
                CCutil_stop_timer (&lp->stats.total, 0);
            }
            printf ("Final LP has %d rows, %d columns, %d nonzeros\n",
                    CClp_nrows (lp->lp), CClp_ncols (lp->lp),
                    CClp_nnonzeros (lp->lp));

            if (dat.ndepot > 0) {
                rval = CCtsp_depot_valid (lp, dat.ndepot, (int *) NULL);
                CCcheck_rval (rval, "CCtsp_depot_valid failed");
            }
            goto DONE;
        }

        if (dat.ndepot == 0 && eliminate_edges) {
            rval = CCtsp_eliminate_variables (lp, eliminate_sparse, silent);
            CCcheck_rval (rval, "CCtsp_eliminate_variables failed");
        }
    } else {
        printf ("During testing, do not exact price large problems\n");
        fflush (stdout);
        CCutil_stop_timer (&lp->stats.total, 1);
        printf ("Final LP has %d rows, %d columns, %d nonzeros\n",
                CClp_nrows (lp->lp), CClp_ncols (lp->lp),
                CClp_nnonzeros (lp->lp));

        goto DONE;
    }

    CCutil_stop_timer (&lp->stats.total, 1);
    printf ("Final LP has %d rows, %d columns, %d nonzeros\n",
            CClp_nrows (lp->lp), CClp_ncols (lp->lp),
            CClp_nnonzeros (lp->lp));
    fflush (stdout);

    if (dat.ndepot > 0) {
        rval = CCtsp_depot_valid (lp, dat.ndepot, (int *) NULL);
        CCcheck_rval (rval, "CCtsp_depot_valid failed");
        goto DONE;
    }
    
    if (dfs_branching) {
        upbound = lp->upperbound;
        bbcount = 0;

        if (simple_branching) CCtsp_init_simple_cutselect (&sel);
        rval = CCtsp_easy_dfs_brancher (lp, &sel, 0, &upbound, &bbcount,
                     usebranchcliques, besttour, longedge_branching,
                     simple_branching, silent, &rstate);
        CCcheck_rval (rval, "CCtsp_easy_dfs_brancher failed");
    } else if (bfs_branching) {
        double lowbound = lp->lowerbound;
        int id          = lp->id;

        upbound  = lp->upperbound;
        bbcount = 0;

        rval = CCtsp_write_probroot_id (problname, lp);
        CCcheck_rval (rval, "CCtsp_write_probroot_id failed");
        CCtsp_free_tsp_lp_struct (&lp);

        rval = CCtsp_bfs_brancher (problname, id, lowbound, &sel, 
                &tentativesel, &upbound, &bbcount, usebranchcliques, &dat,
                ptour, pool, ncount, besttour, hostport, &branchzeit,
                save_proof, tentative_branch_num, longedge_branching,
                (double *) NULL, (int *) NULL, silent, &rstate);
        CCcheck_rval (rval, "CCtsp_bfs_brancher failed");
    }

DONE:

    if (dfs_branching || bfs_branching || restartfname) {
        printf ("Optimal Solution: %.2f\n", upbound);
        printf ("Number of bbnodes: %d\n", bbcount);
        fflush (stdout);
        rval = CCtsp_dumptour (ncount, &dat, ptour, probname, besttour,
                               outfname, output_tour_as_edges, silent);
        CCcheck_rval (rval, "CCtsp_dumptour failed");
    } else {
        rval = CCtsp_write_probfile_sav (lp);
        CCcheck_rval (rval, "CCtsp_write_probfile_sav failed");
    }

    printf ("Total Running Time: %.2f (seconds)", CCutil_zeit () - szeit);
    if (branchzeit != 0.0) {
        printf ("  Branching Time: %.2f (seconds)", branchzeit);
    }
    printf ("\n"); fflush (stdout);

    /*  CCtsp_output_statistics (&lp->stats);  */

    if (pool && pool->cutcount) {
        if (!silent) {
            printf ("Final Pool: %d cuts\n", pool->cutcount); fflush (stdout);
        }
        sprintf (buf, "%s.pul", probname);
        rval = CCtsp_write_cutpool (ncount, buf, pool);
        CCcheck_rval (rval, "CCtsp_write_cutpool failed");
    }

#ifdef CCtsp_USE_DOMINO_CUTS
    if (dominopool && dominopool->cutcount) {
        if (1 || !silent) {
            printf ("Final Domino Pool: %d cuts\n", dominopool->cutcount);
            fflush (stdout);
        }
        sprintf (buf, "%s.dominopul", probname);
        rval = CCtsp_write_cutpool (ncount, buf, dominopool);
        CCcheck_rval (rval, "CCtsp_write_cutpool failed");
    }
#endif

    if (sel.remotepool && pool && pool->cutcount > pool->savecount) {
        rval = CCtsp_send_newcuts (ncount, pool, sel.remotehost,
                sel.remoteport);
        if (rval) {
            fprintf (stderr, "CCtsp_send_newcuts failed\n");
            rval = 0;
        }
    }
        
    rval = 0;

CLEANUP:

    if (unlink_files) {
        if (!run_silently) {
            printf ("Delete the temporary files: pul sav mas\n");
            fflush (stdout);
        }

        sprintf (buf, "%s.pul", probname);
        rval = unlink (buf);
        if (rval && !run_silently) {
            printf ("CCutil_sdelete_file failed for %s\n", buf);
        }

        sprintf (buf, "O%s.pul", probname);
        rval = unlink (buf);
        if (rval && !run_silently) {
            printf ("CCutil_sdelete_file failed for %s\n", buf);
        }

        sprintf (buf, "%s.sav", probname);
        rval = unlink (buf);
        if (rval && !run_silently) {
            printf ("CCutil_sdelete_file failed for %s\n", buf);
        }

        sprintf (buf, "O%s.sav", probname);
        rval = unlink (buf);
        if (rval && !run_silently) {
            printf ("CCutil_sdelete_file failed for %s\n", buf);
        }

        sprintf (buf, "%s.mas", probname);
        rval = unlink (buf);
        if (rval && !run_silently) {
            printf ("CCutil_sdelete_file failed for %s\n", buf);
        }

        sprintf (buf, "O%s.mas", probname);
        rval = unlink (buf);
        if (rval && !run_silently) {
            printf ("CCutil_sdelete_file failed for %s\n", buf);
        }
    }

    if (lp) CCtsp_free_tsp_lp_struct (&lp);
    if (pool) { CCtsp_free_cutpool (&pool); }
    if (dominopool) { CCtsp_free_cutpool (&dominopool); }

    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    CC_IFFREE (exlist, int);
    CC_IFFREE (exlen, int);
    CC_IFFREE (ptour, int);
    CC_IFFREE (besttour, int);
    CC_IFFREE (probname, char);
    CCutil_freedatagroup (&dat);

    return rval;
}

static int handle_just_cuts (CCtsp_lp *lp, int the_cuts, CCrandstate *rstate,
       int silent)
{
    int rval = 0;
    CCtsp_cutselect sel;

    if (the_cuts == CC_JUST_FAST_CUTS) {
        CCtsp_init_fast_cutselect (&sel);
        rval = CCtsp_cutselect_set_tols (&sel, lp, -1, silent);
        CCcheck_rval (rval, "CCtsp_cutselect_set_tols failed");
        rval = CCtsp_cutting_loop (lp, &sel, 1, silent, rstate);
        CCcheck_rval (rval, "CCtsp_cutting_loop failed");
    } else if (the_cuts == CC_JUST_SUBTOUR) {
        rval = CCtsp_subtour_loop (lp, silent, rstate);
        CCcheck_rval (rval, "CCtsp_subtour_loop failed");
    } else if (just_cuts == CC_JUST_BLOSSOM) {
        rval = CCtsp_blossom_loop (lp, silent, rstate);
        CCcheck_rval (rval, "CCtsp_blossom_loop failed");
    } else if (just_cuts == CC_JUST_SUBTOUR_AND_BLOSSOM) {
        rval = CCtsp_subtour_and_blossom_loop (lp, silent, rstate);
        CCcheck_rval (rval, "CCtsp_subtour_and_blossom_loop failed");
    }

    printf ("Bound: %f\n", lp->lowerbound); fflush (stdout);
    CCutil_stop_timer (&lp->stats.total, 1);
    printf ("Final Root LP has %d rows, %d columns, %d nonzeros\n",
            CClp_nrows (lp->lp), CClp_ncols (lp->lp), CClp_nnonzeros (lp->lp));

CLEANUP:

    return rval;
}

static int run_hk (int ncount, CCdatagroup *dat, int *hk_tour)
{
    double hk_val;
    int hk_found, hk_yesno;
    int *hk_tlist = (int *) NULL;
    int rval = 0;

    hk_tlist = CC_SAFE_MALLOC (2*ncount, int);
    CCcheck_NULL (hk_tlist, "out of memory for hk_tlist");

    rval = CCheldkarp_small (ncount, dat, (double *) NULL, &hk_val,
                             &hk_found, 0, hk_tlist, 1000000, 2);
    CCcheck_rval (rval, "CCheldkarp_small failed");
    printf ("Optimal Solution: %.2f\n", hk_val); fflush (stdout);

    rval = CCutil_edge_to_cycle (ncount, hk_tlist, &hk_yesno, hk_tour);
    CCcheck_rval (rval, "CCutil_edge_to_cycle failed");

    if (hk_yesno == 0) {
        fprintf (stderr, "Held-Karp returned list that is not a tour\n");
        rval = 1;  goto CLEANUP;
    }

CLEANUP:

     CC_IFFREE (hk_tlist, int);
     return rval;
}

static void adjust_upbound (double *bound, int ncount, CCdatagroup *dat)
{
    double bnd;
    int i;

    bnd = CCutil_dat_edgelen (ncount - 1, 0, dat);
    for (i = 1; i < ncount; i++) {
        bnd += CCutil_dat_edgelen (i-1, i, dat);
    }
    if (bnd < *bound) {
        printf ("Set initial upperbound to %.0f (from tour)\n", bnd);
        fflush (stdout);
        *bound = bnd;
    }
}

static int build_edges (int *p_ecount, int **p_elist, int **p_elen,
        int ncount, int *ptour, CCdatagroup *dat, char *in_edgefname,
        char *in_edgegenfname, int in_just_cuts, int silent,
        CCrandstate *rstate)
{
    int rval = 0;
    int *elist = (int *) NULL;
    int ecount;
    int i;
    
    if (in_edgefname) {
        int *invperm = (int *) NULL;

        printf ("Read initial edge set\n"); fflush (stdout);
        
        rval = CCutil_getedgelist (ncount, in_edgefname, p_ecount, p_elist,
                                   p_elen, 0);
        CCcheck_rval (rval, "CCutil_getedgelist failed");
        ecount = *p_ecount;
        elist = *p_elist;
        printf ("Initial edgeset: %d edges (%d nodes)\n", ecount, ncount);
        printf ("Rearrange the edges to match the tour order\n");
        fflush (stdout);

        invperm = CC_SAFE_MALLOC (ncount, int);
        CCcheck_NULL (invperm, "out of memory for invperm");
        for (i = 0; i < ncount; i++) invperm[ptour[i]] = i;
        for (i = 0; i < 2*ecount; i++) elist[i] = invperm[elist[i]];
        CC_FREE (invperm, int);
    } else if (dat) {
        CCedgegengroup plan;
        
        if (in_edgegenfname) {
            rval = CCedgegen_read (in_edgegenfname, &plan);
            CCcheck_rval (rval, "CCedgegen_read failed");
        } else {
            CCedgegen_init_edgegengroup (&plan);
            if (in_just_cuts == CC_JUST_SUBTOUR ||
                in_just_cuts == CC_JUST_BLOSSOM ||
                in_just_cuts == CC_JUST_SUBTOUR_AND_BLOSSOM) {
                plan.tour.greedy = 1;
                plan.f2match_nearest.number = 4;
            } else {
                plan.linkern.count = 10;
                plan.linkern.quadnearest = 2;
                plan.linkern.greedy_start = 0;
                plan.linkern.nkicks = (ncount / 100) + 1;
            }
        }

        rval = getedges (dat, &plan, ncount, p_ecount, p_elist, p_elen,
                         silent, rstate);
        CCcheck_rval (rval, "getedges failed");
    }

CLEANUP:

    return rval;
}

static int build_fulledges (int *p_excount, int **p_exlist, int **p_exlen,
        int ncount, int *ptour, char *in_fullfname)
{
    int i;
    int rval = 0;
    int *exlist;
    int excount;
    
    if (in_fullfname) {
        int *invperm = (int *) NULL;

        rval = CCutil_getedgelist (ncount, in_fullfname, p_excount, p_exlist,
                                   p_exlen, 0);
        CCcheck_rval (rval, "CCutil_getedgelist failed");

        invperm = CC_SAFE_MALLOC (ncount, int);
        CCcheck_NULL (invperm, "out of memory for invperm");
        for (i = 0; i < ncount; i++) invperm[ptour[i]] = i;
        excount = *p_excount;
        exlist = *p_exlist;
        for (i = 0; i < 2*excount; i++) exlist[i] = invperm[exlist[i]];
        CC_FREE (invperm, int);
    } else {
        *p_excount = 0;
    }

CLEANUP:

    return rval;
}

static int find_tour (int ncount, CCdatagroup *dat, int *perm, double *ub,
        int trials, int silent, CCrandstate *rstate)
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
    int kicks, i, istour;

    szeit = CCutil_zeit ();
    bestval = CCtsp_LP_MAXDOUBLE;

    if (trials == -1) {
        kicks = (ncount > 400 ? 100 : ncount/4);
    } else {
        kicks = (ncount > 1000 ? 500 : ncount/2);
    }

    if (!silent) {
        printf ("Finding a good tour for compression: %d\n", trials);
        fflush (stdout);
    }

    cyc = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (cyc, "out of memory for cyc");
    bestcyc = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (bestcyc, "out of memory for bestcyc");

    CCedgegen_init_edgegengroup (&plan);
    plan.quadnearest = 2;
    rval = CCedgegen_edges (&plan, ncount, dat, (double *) NULL, &ecount,
                            &elist, silent, rstate);
    CCcheck_rval (rval, "CCedgegen_edges failed");
    plan.quadnearest = 0;

    plan.tour.greedy = 1;
    rval = CCedgegen_edges (&plan, ncount, dat, (double *) NULL, &tcount,
                            &tlist, silent, rstate);
    CCcheck_rval (rval, "CCedgegen_edges failed");

    if (tcount != ncount) {
        fprintf (stderr, "wrong edgeset from CCedgegen_edges\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_edge_to_cycle (ncount, tlist, &istour, cyc);
    CCcheck_rval (rval, "CCutil_edge_to_cycle failed");
    if (istour == 0) {
        fprintf (stderr, "Starting tour has an error\n");
        rval = 1; goto CLEANUP;
    }
    CC_FREE (tlist, int);

    rval = CClinkern_tour (ncount, dat, ecount, elist, ncount, kicks,
                    cyc, bestcyc, &bestval, silent, 0.0, 0.0,
                    (char *) NULL,
                    CC_LK_GEOMETRIC_KICK, rstate);
    CCcheck_rval (rval, "CClinkern_tour failed");

    for (i = 0; i < trials; i++) {
        rval = CClinkern_tour (ncount, dat, ecount, elist, ncount, kicks,
                        (int *) NULL, cyc, &val, silent, 0.0, 0.0,
                        (char *) NULL, CC_LK_GEOMETRIC_KICK, rstate);
        CCcheck_rval (rval, "CClinkern_tour failed");
        if (val < bestval) {
            CC_SWAP (cyc, bestcyc, tmp);
            bestval = val;
        }
    }

    if (trials > 0) {
        rval = CClinkern_tour (ncount, dat, ecount, elist, ncount, 2 * kicks,
                        bestcyc, perm, ub, silent, 0.0, 0.0,
                        (char *) NULL, CC_LK_GEOMETRIC_KICK, rstate);
        CCcheck_rval (rval, "CClinkern_tour failed");
    } else {
        for (i = 0; i < ncount; i++) {
            perm[i] = bestcyc[i];
        }
    }

    if (!silent) {
        printf ("Time to find compression tour: %.2f (seconds)\n",
                CCutil_zeit() - szeit);
        fflush (stdout);
    }

CLEANUP:

    CC_IFFREE (cyc, int);
    CC_IFFREE (bestcyc, int);
    CC_IFFREE (elist, int);
    CC_IFFREE (tlist, int);
    return rval;
}

static int getedges (CCdatagroup *dat, CCedgegengroup *plan, int ncount,
        int *ecount, int **elist, int **elen, int silent,
        CCrandstate *rstate)
{
    int i;
    int rval = 0;

    *elist = (int *) NULL;
    *elen = (int *) NULL;

    if (dat == (CCdatagroup *) NULL || plan == (CCedgegengroup *) NULL) {
        fprintf (stderr, "getedges needs CCdatagroup and CCedgegengroup\n");
        rval = 1;  goto CLEANUP;
    }

    rval = CCedgegen_edges (plan, ncount, dat, (double *) NULL, ecount, elist,
                            silent, rstate);
    CCcheck_rval (rval, "CCedgegen_edges failed");

    *elen = CC_SAFE_MALLOC(*ecount, int);
    CCcheck_NULL (*elen, "out of memory for elen");

    for (i = 0; i < *ecount; i++) {
        (*elen)[i] = CCutil_dat_edgelen ((*elist)[2*i], (*elist)[(2*i)+1], dat);
    }

CLEANUP:

    if (rval) {
        CC_IFFREE (*elist, int);
        CC_IFFREE (*elen, int);
    }
    return rval;
}

static int dump_rc (CCtsp_lp *lp, int count, char *pname, int usesparse)
{
    int rval = 0;
    char rcnname[1024];

    sprintf (rcnname, "%s.rcn", pname);
    rval = CCtsp_dump_rc_nearest (lp, count, rcnname, usesparse);
    CCcheck_rval (rval, "CCtsp_dump_rc failed");

CLEANUP:

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c, inorm;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    /* remaining: aAbcGHjlLpQWY */

    while ((c = CCutil_bix_getopt (ac, av, "BC:dD:e:E:fF:g:hiIJ:k:K:mM:n:N:o:P:qr:R:s:S:t:T:u:UvVwX:xyz:Z:", &boptind, &boptarg)) != EOF)
        switch (c) {
        case 'B':
            bfs_branching = 0;
            break;
        case 'C':
            maxchunksize = atoi(boptarg);
            break;
        case 'd':
            dfs_branching = 1;
	    bfs_branching = 0;
            break;
        case 'D':
            edgegenfname = boptarg;
            break;
        case 'e':
            edgefname = boptarg;
            break;
        case 'E':
            fullfname = boptarg;
	    valid_edges = 1;
            break;
        case 'f':
            output_tour_as_edges = 1;
            break;
        case 'F':
            filecutname = boptarg;
            break;
        case 'g':
            grunthostname = boptarg;
            break;
        case 'h':
            be_nethost = 1;
            break;
        case 'i':
            just_cuts = CC_JUST_BLOSSOM;
            break;
        case 'I':
            just_cuts = CC_JUST_SUBTOUR;
            break;
        case 'J':
            tentative_branch_num = atoi (boptarg);
            break;
        case 'k':
            nnodes_want = atoi (boptarg);
            break;
        case 'K':
            cutbossname = boptarg;
            break;
        case 'm':
            multiple_chunker = 1;
            break;
        case 'M':
            masterfname = boptarg;
            break;
        case 'n':
            problname = boptarg;
            break;
        case 'o':
            outfname = boptarg;
            break;
        case 'P':
            poolfname = boptarg;
            break;
        case 'q':
            dontcutroot = 1;
            break;
        case 'r':
            gridsize = atoi(boptarg);
            break;
        case 'R':
            restartfname = boptarg;
            break;
        case 's':
            seed = atoi (boptarg);
            break;
        case 'S':
            probfname = boptarg;
            break;
        case 't':
            tourfname =  boptarg;
            break;
#ifdef CCtsp_USE_DOMINO_CUTS
        case 'Z':
            usedominos = atoi(boptarg);
            break;
        case 'T':
            dombossname = boptarg;
            if (!usedominos) usedominos = 1;
            break;
#endif
        case 'u':
            initial_ub = atof (boptarg);
            break;
        case 'U':
            usebranchcliques = 0;
            break;
        case 'v':
            run_silently = 0;
            break;
        case 'V':
            just_cuts = CC_JUST_FAST_CUTS;
            maxchunksize = 0;
            break;
        case 'w':
            just_cuts = CC_JUST_SUBTOUR_AND_BLOSSOM;
            break;
        case 'X':
            xfname = boptarg;
            break;
        case 'x':
            unlink_files = 1;
            break;
        case 'y':
            simple_branching = 1;
            break;
        case 'z':
            want_rcnearest = atoi (boptarg);
            break;
        case 'N':
            inorm = atoi (boptarg);
            switch (inorm) {
            case 0: norm = CC_MAXNORM; break;
            case 1: norm = CC_MANNORM; break;
            case 2: norm = CC_EUCLIDEAN; break;
            case 3: norm = CC_EUCLIDEAN_3D; break;
            case 4: norm = CC_USER; break;
            case 5: norm = CC_ATT; break;
            case 6: norm = CC_GEOGRAPHIC; break;
            case 7: norm = CC_MATRIXNORM; break;
            case 8: norm = CC_DSJRANDNORM; break;
            case 9: norm = CC_CRYSTAL; break;
            case 10: norm = CC_SPARSE; break;
            case 11: norm = CC_RHMAP1; break;
            case 12: norm = CC_RHMAP2; break;
            case 13: norm = CC_RHMAP3; break;
            case 14: norm = CC_RHMAP4; break;
            case 15: norm = CC_RHMAP5; break;
            case 16: norm = CC_EUCTOROIDAL; break;
            case 17: norm = CC_GEOM; break;
            case 18: norm = CC_EUCLIDEAN_CEIL; break;
            default:
                usage (av[0]);
                return 1;
            }
            tsplib_in = 0;
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    if (boptind < ac) {
        datfname = av[boptind++];
    }

    if (boptind != ac) {
        usage (av[0]);
        return 1;
    }

    if (datfname == (char *) NULL && nnodes_want == 0 && 
        probfname == (char *) NULL && edgefname == (char *) NULL &&
        masterfname == (char *) NULL && grunthostname == (char *) NULL) {
        usage (av[0]);
        return 1;
    }

    if (datfname == (char *) NULL && masterfname == (char *) NULL &&
        edgefname != (char *) NULL) {
        fprintf (stderr, "cannot give edgefile without a dat or master file\n");
        return 1;
    }

    if (eliminate_edges < 0) {
        if (bfs_branching || dfs_branching) {
            eliminate_edges = 1;
        } else {
            eliminate_edges = 0;
        }
    }
    
    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s [-see below-] [dat_file]\n", f);
    fprintf (stderr, "   -B    do not branch\n");
    fprintf (stderr, "   -C #  maximum chunk size in localcuts (default 16)\n");
    fprintf (stderr, "   -d    use dfs branching instead of bfs\n");
    fprintf (stderr, "   -D f  edgegen file for initial edge set\n");
    fprintf (stderr, "   -e f  initial edge file\n");
    fprintf (stderr, "   -E f  full edge file (must contain initial edge set)\n");
    fprintf (stderr, "   -f    write optimal tour as edge file (default is tour file)\n");
    fprintf (stderr, "   -F f  read extra cuts from file\n");
    fprintf (stderr, "   -g h  be a grunt for boss h\n");
    fprintf (stderr, "   -h    be a boss for the branching\n");
    fprintf (stderr, "   -i    just solve the blossom polytope\n");
    fprintf (stderr, "   -I    just solve the subtour polytope\n");
    fprintf (stderr, "   -J #  number of tentative branches\n");
    fprintf (stderr, "   -k #  number of nodes for random problem\n");
    fprintf (stderr, "   -K h  use cut server h\n");
    fprintf (stderr, "   -M f  master file\n");
    fprintf (stderr, "   -m    use multiple passes of cutting loop\n");
    fprintf (stderr, "   -n s  problem location (just a name or host:name, not a file name)\n");
    fprintf (stderr, "   -o f  output file name (for optimal tour)\n");
    fprintf (stderr, "   -P f  cutpool file\n");
    fprintf (stderr, "   -q    do not cut the root lp\n");
    fprintf (stderr, "   -r #  use #x# grid for random points, no dups if #<0\n");
    fprintf (stderr, "   -R f  restart file\n");
    fprintf (stderr, "   -s #  random seed\n");
    fprintf (stderr, "   -S f  problem file\n");
    fprintf (stderr, "   -t f  tour file (in node node node format)\n");
#ifdef CCtsp_USE_DOMINO_CUTS
    fprintf (stderr, "   -Z #  dp-cuts (#=1 normal, #=2 shrunk, #=3 both\n");
    fprintf (stderr, "   -T h  use domino server h\n");
#endif
    fprintf (stderr, "   -u v  initial upperbound\n");
    fprintf (stderr, "   -U    do not permit branching on subtour inequalities\n");
    fprintf (stderr, "   -v    verbose (turn on lots of messages)\n");
    fprintf (stderr, "   -V    just run fast cuts\n");
    fprintf (stderr, "   -w    just subtours and trivial blossoms\n");
    fprintf (stderr, "   -x    delete files on completion (sav pul mas)\n");
    fprintf (stderr, "   -X f  write the last root fractional solution to f\n");
    fprintf (stderr, "   -y    use simple cutting and branching in DFS\n");
    fprintf (stderr, "   -z #  dump the #-lowest reduced cost edges to file xxx.rcn\n");
    fprintf (stderr, "   -N #  norm (must specify if dat file is not a TSPLIB file)\n");
    fprintf (stderr, "         0=MAX, 1=L1, 2=L2, 3=3D, 4=USER, 5=ATT, 6=GEO, 7=MATRIX,\n");
    fprintf (stderr, "         8=DSJRAND, 9=CRYSTAL, 10=SPARSE, 11-15=RH-norm 1-5, 16=TOROIDAL\n");
    fprintf (stderr, "         17=GEOM, 18=JOHNSON\n");
}

#include <stdio.h>
#include "machdefs.h"
#include "util.h"
#include "tsp.h"


static int process_subproblem (char *probname, int id, int ncount,
        CCdatagroup *dat, int *ptour, double *lbound, CCrandstate *rstate);
static int build_edges (CCdatagroup *dat, int ncount, int *ecount, int **elist,
         int **elen, int silent, CCrandstate *rstate);

static FILE *jokeout = (FILE *) NULL;
static int jokecount = 0;
static int jokeneg = 0;
static double *jokepi = (double *) NULL;

int main (int ac, char **av)
{
    char *bosshost = (char *) NULL;
    double lbound = -1.0, rtime = 0.0;
    char probname[CCutil_FILE_NAME_LEN];
    int id = -1;
    int rval = 0;
    CC_SFILE *s = (CC_SFILE *) NULL;
    double szeit;
    CCrandstate rstate;
    int seed = 99;
    int ncount;
    int *perm = (int *) NULL;
    CCdatagroup dat;

    if (ac != 2) {
        fprintf (stderr, "Usage: %s boss\n", av[0]);
        rval = 1; goto CLEANUP;
    }

    CCutil_printlabel ();

    CCutil_sprand (seed, &rstate);


    /* Hacky for subtours */
    jokeout = fopen ("all.cuts", "w");
    if (!jokeout) {
        fprintf (stderr, "unable to open all.cuts for writing\n");
        rval = 1;  goto CLEANUP;
    }
    jokepi = CC_SAFE_MALLOC (100000, double);
    CCcheck_NULL (jokepi, "out of memory for jokepi");
    {
        int i;
        for (i = 0; i < 100000; i++) jokepi[i] = -1000.0;
    }

    /* end Hacky for subtours */

    bosshost = av[1];

    while (1) {
        s = CCutil_snet_open (bosshost, CC_SUBDIV_PORT);
        if (!s) {
            fprintf (stderr, "CCutil_snet_open failed\n");
            rval = 1;  goto CLEANUP;
        }

        rval = CCutil_swrite_int (s, id);
        CCcheck_rval (rval, "CCutil_swrite_int failed (id)");
        rval = CCutil_swrite_double (s, rtime);
        CCcheck_rval (rval, "CCutil_swrite_double failed (rtime)");
        rval = CCutil_swrite_double (s, lbound);
        CCcheck_rval (rval, "CCutil_swrite_double failed (rtime)");

        rval = CCutil_sread_int (s, &id);
        CCcheck_rval (rval, "CCutil_sread_int failed (id)");
  
        if (id == -1) {
            CCutil_sclose (s);
            s = (CC_SFILE *) NULL;
            goto DONE;
        }

        rval = CCutil_sread_string (s, probname, CCutil_FILE_NAME_LEN);
        CCcheck_rval (rval, "CCutil_sread_string failed (probname)");

        CCutil_init_datagroup (&dat);
        perm = (int *) NULL;

        rval = CCutil_readmaster (s, &ncount, &dat, &perm);
        CCcheck_rval (rval, "CCutil_readmaster failed");
        
        CCutil_sclose (s);
        s = (CC_SFILE *) NULL;

        printf ("PROCESSING %s subproblem %d\n", probname, id);
        fflush (stdout);

        szeit = CCutil_zeit ();

        rval = process_subproblem (probname, id, ncount, &dat, perm, &lbound,
                                   &rstate);
        CCcheck_rval (rval, "process_subproblem failed");

        rtime = CCutil_zeit () - szeit;

        CC_IFFREE (perm, int);
        CCutil_freedatagroup (&dat);
    }

DONE:
    
    printf ("No work available.  Shutting down.\n"); fflush (stdout);
    printf ("Found %d cuts in total\n", jokecount);
    printf ("Found %d negative pi values\n", jokeneg);
    {
        FILE *jokepiout = fopen ("all.pi", "w");
        int i = 0;

        if (!jokepiout) {
            fprintf (stderr, "could not open all.pi for writing\n");
            rval = 1;  goto CLEANUP;
        }

        while (jokepi[i] != -1000.0) {
            fprintf (jokepiout, "%f\n", jokepi[i]);
            i++;
        }
        printf ("Found %d pi values\n", i);
        fclose (jokepiout);
    }

CLEANUP:

    if (s != (CC_SFILE *) NULL) {
        CCutil_sclose (s);
    }

    if (jokeout) fclose (jokeout);
    CC_IFFREE (jokepi, double);

    return rval;
}

static int process_subproblem (char *probname, int id, int ncount,
        CCdatagroup *dat, int *ptour, double *lbound, CCrandstate *rstate) 
{
    char name[CCutil_FILE_NAME_LEN];
    int i, ecount, yesno, rval = 0, silent = 1;
    int *elist = (int *) NULL;
    int *elen = (int *) NULL;
    int *besttour = (int *) NULL;
    CCtsp_lp *lp = (CCtsp_lp *) NULL;
    CCtsp_lpcuts *pool = (CCtsp_lpcuts *) NULL;
    double ptour_len, lb;
    CCtsp_cutselect sel;
    CCbigguy exbound;


    printf ("ncount = %d, Depots = %d\n", ncount, dat->ndepot);
    fflush (stdout);

    CCutil_cycle_len (ncount, dat, ptour, &ptour_len);
    printf ("initial tour: %.2f\n", ptour_len); fflush (stdout);

    rval = build_edges (dat, ncount, &ecount, &elist, &elen, silent, rstate);
    CCcheck_rval (rval, "build_edges failed");

    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &pool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");

    besttour = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (besttour, "out or memory in process_subproblem");
    for (i = 0; i < ncount; i++) besttour[i] = i;

    sprintf (name, "%s_%d", probname, id);

    rval = CCtsp_init_lp (&lp, name,  -1, (char *) NULL, ncount, dat,
                    ecount, elist, elen, 0, (int *) NULL, (int *) NULL,
                    0, ptour, ptour_len, pool, 0, rstate);
    CCcheck_rval (rval, "CCtsp_init_lp failed");
   
    CCtsp_init_cutselect (&sel);
    rval = CCtsp_get_lp_result (lp, &lb, (double *) NULL, (int *) NULL,
              (int **) NULL, (double **) NULL, (double **) NULL,
              (double **) NULL, (double **) NULL);
    CCcheck_rval (rval, "CCtsp_get_lp_result failed");
    sel.roundtol = 0.0001 * lb;
    sel.nexttol  = 0.001 * lb;
    printf ("Setting tolerances: next cuts %.4f next round %.4f\n",
                sel.nexttol, sel.roundtol);
    fflush (stdout);
/*
    rval = CCtsp_cutselect_set_tols (&sel, lp, 1, silent);
    CCcheck_rval (rval, "CCtsp_cutselect_set_tols failed");
*/

/*
    rval = CCtsp_cutting_loop (lp, &sel, 1, silent, rstate);
    CCcheck_rval (rval, "CCtsp_cutting_loop failed");
*/

/*
    rval = CCtsp_cutting_multiple_loop (lp, &sel, 1, 24, 0, silent, rstate);
    CCcheck_rval (rval, "CCtsp_cutting_multiple_loop failed");
*/

    /* Hacky code for subtour picture */

    rval = CCtsp_subtour_loop (lp, silent, rstate);
    CCcheck_rval (rval, "CCtsp_subtour_loop failed");

    {
        int j, k, acount, ocount;
        CCtsp_lpcut_in c;
        int *ar = (int *) NULL;
        int orig = lp->dat->orig_ncount;
        int ncuts = lp->cuts.cutcount;
        int nrows = ncount + ncuts;
        double *pi = (double *) NULL;
        double *cut_pi;
        double cx;

/*
        for (i = 0; i < ncount; i++) {
            rval = CClp_change_sense (lp->lp, i, 'G');
            CCcheck_rval (rval, "CClp_change_sense failed");
        }
    
        rval = CClp_opt (lp->lp, CClp_METHOD_DUAL);
        CCcheck_rval (rval, "CClp_opt failed");
*/

        pi = CC_SAFE_MALLOC (nrows, double);
        CCcheck_NULL (pi, "out of memory in print_the_subs");

        rval = CClp_pi (lp->lp, pi);
        CCcheck_rval (rval, "CClp_pi failed");
        cut_pi = pi + ncount;

        for (i = 0; i < ncuts; i++) {
            cx = cut_pi[i];
            for (k = 0; k < lp->cuts.cuts[i].modcount; k++) {
                pi[lp->cuts.cuts[i].mods[k].node] += cx *
                    (((int) lp->cuts.cuts[i].mods[k].mult) - 128);
            }
        }

        for (i = 0; i < ncount; i++) {
            if (lp->perm[i] < orig) {
                if (pi[i] < -0.001) {
                    fprintf (stderr, "Oh no, negative %f\n", pi[i]);
                    jokeneg++;
                } 
            }
        }

        for (i = 0; i < ncount; i++) {
            if (lp->perm[i] < orig) {
              jokepi[lp->dat->orig_names[lp->perm[i]]] = pi[i];
            }
        }

        for (i = 0; i < ncuts; i++) {
            if (cut_pi[i] < 0.001) continue;
            rval = CCtsp_lpcut_to_lpcut_in (&lp->cuts, &lp->cuts.cuts[i], &c);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");

            if (c.cliquecount != 1) {
                fprintf (stderr, "HEY!  This is not a subtour\n"); 
                exit (1);
            }
 
            rval = CCtsp_clique_to_array (&c.cliques[0], &ar, &acount);
            if (rval) {
                fprintf (stderr, "CCtsp_clique_to_array failed");
                rval = 1;  goto CLEANUP;
            }

            ocount = 0;
            for (j = 0; j < acount; j++) {
                if (lp->perm[ar[j]] >= orig) {
                    ocount++;
                    ar[j] = ar[acount-1];
                    acount--;
                    j--;
                }
            }

            if (ocount == 0 && acount >= 3) {
                for (j = 0; j < acount; j++) {
                    ar[j] = lp->dat->orig_names[lp->perm[ar[j]]];
                }
                CCutil_int_array_quicksort (ar, acount);

                fprintf (jokeout, "%f %d ", cut_pi[i], acount);
                for (j = 0; j < acount; j++) {
                    fprintf (jokeout, "%d ", ar[j]);
                }
                fprintf (jokeout, "\n");
                jokecount++;
            }

            CCtsp_free_lpcut_in (&c);
            CC_IFFREE (ar, int);
        }
        CC_IFFREE (pi, double);
    }

    /* end subtour picture hack */

/*
    {
        double tourval;
        CCutil_start_timer (&lp->stats.linkern);
        rval = CCtsp_call_x_heuristic (lp, &tourval, besttour, silent,
                                       rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_call_x_heuristic failed\n");
            goto CLEANUP;
        }
        CCutil_stop_timer (&lp->stats.linkern, 0);
        if (tourval < lp->upperbound) {
            printf ("New upperbound from x-heuristic: %.2f\n", tourval);
            lp->upperbound = tourval;
        }
    }
*/


    printf ("Final LP has %d rows, %d columns, %d nonzeros\n",
            CClp_nrows (lp->lp), CClp_ncols (lp->lp),
            CClp_nnonzeros (lp->lp));
    printf ("Final lower bound %f, upper bound %f\n", lp->lowerbound,
                                                      lp->upperbound);
    fflush (stdout);

    rval = CCtsp_exact_price (lp, &exbound, 1, 0, silent);
    CCcheck_rval (rval, "CCtsp_exact_price failed");

    lp->exact_lowerbound = exbound;
    printf ("Exact lower bound: %.6f\n", CCbigguy_bigguytod (exbound));
    printf ("DIFF: %f\n", lp->lowerbound - CCbigguy_bigguytod (exbound));
    fflush (stdout);

    rval = CCtsp_depot_valid (lp, dat->ndepot, &yesno);
    CCcheck_rval (rval, "CCtsp_depot_valid failed");

    if (yesno == 0) *lbound = -1.0;
    else            *lbound = CCbigguy_bigguytod (exbound);

CLEANUP:

    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    if (pool) CCtsp_free_cutpool (&pool);
    CCtsp_free_tsp_lp_struct (&lp);
    return rval;
}

static int build_edges (CCdatagroup *dat, int ncount, int *ecount, int **elist,
         int **elen, int silent, CCrandstate *rstate)
{
    int i, rval = 0;
    CCedgegengroup plan;

    CCedgegen_init_edgegengroup (&plan);
    plan.linkern.count = 10;
    plan.linkern.quadnearest = 2;
    plan.linkern.greedy_start = 0;
    plan.linkern.nkicks = (ncount / 100) + 1;

    rval = CCedgegen_edges (&plan, ncount, dat, (double *) NULL, ecount,
                            elist, silent, rstate);
    CCcheck_rval (rval, "CCedgegen_edges failed");
    
    *elen = CC_SAFE_MALLOC (*ecount, int);
    CCcheck_NULL (*elen, "out of memory in build_edges");

    for (i = 0; i < *ecount; i++) {
        (*elen)[i] = CCutil_dat_edgelen ((*elist)[2*i], (*elist)[(2*i) + 1],
                                         dat);
    }

CLEANUP:

    return rval;
}

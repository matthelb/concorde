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

#include "machdefs.h"
#include "util.h"
#include "tsp.h"

static char *edgefname = (char *) NULL;
static char *tourfname = (char *) NULL;
static char *outfname  = (char *) NULL;
static int tsplib = 0;
static int seed = 0;
static int run_silently = 0;
static double init_ub = CCtsp_LP_MAXDOUBLE;
static double in_timebound = 0.0;


int
    main (int ac, char **av);

static void
    usage (char *f);

static int
    parseargs (int ac, char **av);



int main (int ac, char **av)
{
    int rval = 0;
    double szeit, optval, *mybnd, *mytimebound;
    int ncount, ecount, success, foundtour, hit_timebound = 0;
    int *elist = (int *) NULL;
    int *elen = (int *) NULL;
    int *in_tour = (int *) NULL;
    int *out_tour = (int *) NULL;
    CCrandstate rstate;
    char *probname = (char *) NULL;
    CCdatagroup dat;

    CCutil_init_datagroup (&dat);

    seed = (int) CCutil_real_zeit ();
    if (parseargs (ac, av)) goto CLEANUP;

    CCutil_sprand (seed, &rstate);
    printf ("Using random seed %d\n", seed); fflush (stdout);

    mybnd = (init_ub == CCtsp_LP_MAXDOUBLE ? (double *) NULL : &init_ub);
    mytimebound = (in_timebound == 0.0 ? (double *) NULL : &in_timebound);

    if (tsplib == 0) {
        rval = CCutil_getedgelist_n (&ncount, edgefname, &ecount, &elist,
                                     &elen, 0);
        if (rval) {
            fprintf (stderr, "CCutil_getedgelist_n failed\n"); goto CLEANUP;
        }
    } else {
        rval = CCutil_gettsplib (edgefname, &ncount, &dat);
        if (rval) {
            fprintf (stderr, "CCutil_gettsplib failed\n"); goto CLEANUP;
        }
    }

    out_tour  = CC_SAFE_MALLOC (ncount, int);
    if (!out_tour) {
        fprintf (stderr, "out of memory in main\n");
        rval = 1; goto CLEANUP;
    }

    probname = CCtsp_problabel (edgefname);
    if (!probname) {
        fprintf (stderr, "CCtsp_problabel failed\n");
        rval = 1; goto CLEANUP;
    }

    if (tourfname) {
        in_tour = CC_SAFE_MALLOC (ncount, int);
        if (!in_tour) {
            fprintf (stderr, "out of memory in main\n");
            rval = 1; goto CLEANUP;
        }
        rval = CCutil_getcycle (ncount, tourfname, in_tour, 0);
        if (rval) {
            fprintf (stderr, "CCutil_getcycle failed\n"); goto CLEANUP;
        }
    }

    szeit = CCutil_zeit ();

    if (tsplib == 0) {
        rval = CCtsp_solve_sparse (ncount, ecount, elist, elen, in_tour,
                  out_tour, mybnd, &optval, &success, &foundtour, 
                  probname, mytimebound, &hit_timebound, run_silently, &rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_solve_sparse failed\n"); goto CLEANUP;
        }
    } else {
        rval = CCtsp_solve_dat (ncount, &dat, in_tour, out_tour, mybnd,
                 &optval, &success, &foundtour, probname, mytimebound,
                 &hit_timebound, run_silently, &rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_solve_dat failed\n"); goto CLEANUP;
        }
    }

    if (hit_timebound) {
        printf ("Code hit time bound\n");
        if (foundtour) {
            printf ("A possibly non-optimal tour has been found\n");
        }
        fflush (stdout);
    } else if (success) {
        if (optval == CCtsp_LP_MAXDOUBLE) {
            printf ("Graph does not contain a tour\n"); fflush (stdout);
        } else {
            if (foundtour == 0) {
                printf ("Lowerbound %.0f established\n", optval);
            } else {
                printf ("Found the optimal tour: %.0f\n", optval);
            }
            fflush (stdout);
        }
    } else {
        printf ("Code hit limit, did not succeed\n");
        if (foundtour) {
            printf ("A possibly non-optimal tour has been found\n");
        }
        fflush (stdout);
    }

    if (foundtour && outfname) {
        rval = CCutil_writecycle (ncount, outfname, out_tour, 0);
        if (rval) {
            fprintf (stderr, "CCutil_writecycle failed\n"); fflush (stdout);
        }
    }

    printf ("Total Running Time: %.2f (seconds)\n", CCutil_zeit () - szeit);
    fflush (stdout);

CLEANUP:

    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    CC_IFFREE (probname, char);
    CCutil_freedatagroup (&dat);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "o:Qs:t:Tu:z:", &boptind, &boptarg)) != EOF)
        switch (c) {
        case 'o':
            outfname = boptarg;
            break;
        case 'Q':
            run_silently = 1;
            break;
        case 's':
            seed = atoi (boptarg);
            break;
        case 't':
            tourfname = boptarg;
            break;
        case 'T':
            tsplib = 1;
            break;
        case 'u':
            init_ub = atof (boptarg);
            break;
        case 'z':
            in_timebound = atof (boptarg);
            break;
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    if (boptind >= ac) {
        usage (av[0]);
        return 1;
    }
    edgefname = av[boptind++];
    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "usage: %s [- below -] file\n", f);
    fprintf (stderr, "   -o f  write tour (node node node format)\n");
    fprintf (stderr, "   -Q    run quietly (don't generate so much output)\n");
    fprintf (stderr, "   -s #  random seed\n");
    fprintf (stderr, "   -t f  starting tour file (node node node format)\n");
    fprintf (stderr, "   -T    file is a TSPLIB file (default is edge file)\n");
    fprintf (stderr, "   -u #  bound on tour length\n");
    fprintf (stderr, "   -z #  bound on running time\n");
}





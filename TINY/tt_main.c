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
/*                        Demo of TINYTSP Solver                            */
/*                                                                          */
/*  Written by:  David Applegate and Bill Cook                              */
/*  Date: April 1, 1997                                                     */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tinytsp.h"

static char *tspfile = (char *) NULL;
static int verify_result = 0;
static int depot = -1;
static int searchlimit = -1;
static int mspfile = 0;
static int haveupper = 0;
static int objsense = CC_TINYTSP_MINIMIZE;
static double upperbound = 0.0;
static int use_bnbtsp = 0;


int
    main (int ac, char **av);

static void
    usage (char *f);

static int
    getmsp (char *fname, int *ncount, int *ecount, int **elist, int **elen,
            int **lower, int **upper),
    parseargs (int ac, char **av);


int main (int ac, char **av)
{
    double optval = 0.0;
    double szeit = 0.0;
    double *up = (double *) NULL;
    int ncount = 0, ecount = 0, i;
    int rval = 0;
    int *elist = (int *) NULL;
    int *elen = (int *) NULL;
    int *lower = (int *) NULL;
    int *upper = (int *) NULL;
    int *xsol = (int *) NULL;
    const char *algname = "unknown";
    int foundtour = 0;

    if (parseargs (ac, av))
        return 0;

    if (use_bnbtsp && depot == -1) {
        fprintf (stderr, "Can only use bnb with msp format\n");
        return 0;
    }

    if (haveupper)
        up = &upperbound;


    if (mspfile) {
        depot = 0;
        rval = getmsp (tspfile, &ncount, &ecount, &elist, &elen, &lower,
                       &upper);
        if (rval) {
            fprintf (stderr, "getmsp failed\n"); goto CLEANUP;
        }
    } else {
        rval = CCutil_getedgelist_n (&ncount, tspfile, &ecount, &elist,
                                     &elen, 0);
        if (rval) {
            fprintf (stderr, "CCutil_getedgelist_n failed\n"); goto CLEANUP;
        }
        lower = CC_SAFE_MALLOC (ecount, int);
        upper = CC_SAFE_MALLOC (ecount, int);
        if (lower == (int *) NULL ||
            upper == (int *) NULL) {
            fprintf (stderr, "Out of memory in main\n"); goto CLEANUP;
        }
        for (i=0; i<ecount; i++) {
            lower[i] = 0;
            upper[i] = 1;
        }
    }

    if (use_bnbtsp) {
        algname = "CCtiny_bnb_msp";
    } else {
        algname = "CCtiny_bnc_msp";
    }
    printf ("Problem: %d nodes, %d edges (%s)\n", ncount, ecount, algname);
    fflush (stdout);

    szeit = CCutil_zeit ();

    xsol = CC_SAFE_MALLOC (ecount, int);
    if (!xsol) {
        fprintf (stderr, "out of memory in main\n"); goto CLEANUP;
    }

    if (use_bnbtsp) {
        rval = CCtiny_bnb_msp (ncount, ecount, elist, elen, depot, lower,
                        upper, up, objsense, &optval, xsol, searchlimit);
    } else {
        rval = CCtiny_bnc_msp (ncount, ecount, elist, elen, depot, lower,
               upper, up, objsense, &optval, xsol, verify_result, searchlimit);
    }

    if (rval == CC_TINYTSP_INFEASIBLE) {
        printf ("There is No Tour\n"); fflush (stdout);
        rval = 0; goto CLEANUP;
    } else if (rval == CC_TINYTSP_SEARCHLIMITEXCEEDED) {
        printf ("Exeeded the search limit\n"); fflush (stdout);
        rval = 0; goto CLEANUP;
    } else if (rval) {
        fprintf (stderr, "CCtiny_bnc_msp failed\n"); goto CLEANUP;
    }

    printf ("Tour Length: %.2f\n", optval);
    printf ("Tour: ");
    for (i = 0; i < ecount; i++) {
        if (xsol[i]) {
            printf ("(%d, %d)", elist[2*i], elist[2*i+1]);
            if (xsol[i] == 2.0) {
                printf ("* ");
            } else {
                printf ("  ");
            }
            fflush (stdout);
        }
    }
    printf ("\n"); fflush (stdout);

CLEANUP:
    if (rval == 0) {
        printf ("Running Time: %.2f (seconds)\n", CCutil_zeit () - szeit);
        fflush (stdout);
    }

    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    CC_IFFREE (xsol, int);
    CC_IFFREE (lower, int);
    CC_IFFREE (upper, int);
    return rval;
}

static int getmsp (char *fname, int *ncount, int *ecount, int **elist,
                   int **elen, int **lower, int **upper)
{
    FILE *in;
    int i;

    *elist = (int *) NULL;
    *elen = (int *) NULL;
    *lower = (int *) NULL;
    *upper = (int *) NULL;

    if ((in = fopen (fname, "r")) == (FILE *) NULL) {
        perror (fname);
        fprintf (stderr, "Unable to open %s for input\n", fname);
        return 1;
    }

    *ncount = CCutil_readint (in);
    *ecount = CCutil_readint (in);

    *elist = CC_SAFE_MALLOC(2 * (*ecount), int);
    *elen = CC_SAFE_MALLOC(*ecount, int);
    *lower = CC_SAFE_MALLOC (*ecount, int);
    *upper = CC_SAFE_MALLOC (*ecount, int);
    if (!(*elist) || !(*elen) || !(*lower) || !(*upper)) {
        fprintf (stderr, "out of memory in getmsp\n");
        fclose (in);
        CC_IFFREE (*elist, int);
        CC_IFFREE (*elen, int);
        CC_IFFREE (*lower, int);
        CC_IFFREE (*upper, int);
        return 1;
    }

    for (i = 0; i < *ecount; i++) {
        (*elist)[2*i] = CCutil_readint (in);
        (*elist)[2*i+1] = CCutil_readint (in);
        (*lower)[i] = CCutil_readint (in);
        (*upper)[i] = CCutil_readint (in);
        (*elen)[i] = CCutil_readint (in);
    }

    fclose (in);
    return 0;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "bd:mMS:u:v", &boptind, &boptarg)) != EOF)
        switch (c) {
        case 'b':
            use_bnbtsp = 1; break;
        case 'd':
            depot = atoi (boptarg); break;
        case 'm':
            mspfile = 1; break;
        case 'M':
            objsense  = CC_TINYTSP_MAXIMIZE; break;
        case 'S':
            searchlimit = atoi (boptarg); break;
        case 'u':
            haveupper = 1;
            upperbound = atof (boptarg); break;
        case 'v':
            verify_result = 1; break;
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    if (boptind >= ac) {
        usage (av[0]);
        return 1;
    }

    tspfile = av[boptind++];
    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "usage: %s [- below -] edge_file\n", f);
    fprintf (stderr, "    b     use bnb tsp solver\n");
    fprintf (stderr, "    d #   use node # as a depot\n");
    fprintf (stderr, "    m     edge_file is in msp format\n");
    fprintf (stderr, "    M     maximize the tour length\n");
    fprintf (stderr, "    S #   limit on number of search nodes\n");
    fprintf (stderr, "    u #   upperbound on tour length\n");
    fprintf (stderr, "    v     verify results\n");
}

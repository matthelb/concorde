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
/*          A PROGRAM TO COMPUTE THE LENGTH OF A TOUR IN A FILE             */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: May 3, 1999                                                       */
/*                                                                          */
/*  SEE short decsribtion in usage ().                                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

static char *cycfname  = (char *) NULL;
static char *tspfname  = (char *) NULL;
static char *edgefname = (char *) NULL;
static int eformat = 0;
static int simpletour = 0;
static int seed = 0;


int
    main (int ac, char **av);

static int
    parseargs (int ac, char **av);

static void
    usage (char *fname);


int main (int ac, char **av)
{
    int ncount, rval = 0;
    int *tour = (int *) NULL;
    double val;
    CCdatagroup dat;
    CCrandstate rstate;

    CCutil_init_datagroup (&dat);

    seed = (int) CCutil_real_zeit ();

    rval = parseargs (ac, av);
    if (rval) return 1;

    if ((!edgefname && !tspfname) || (edgefname && tspfname)) {
        usage (av[0]);
        return 1;
    }

    CCutil_sprand (seed, &rstate);

    if (tspfname) {
        rval = CCutil_gettsplib (tspfname, &ncount, &dat);
        if (rval) {
            fprintf (stderr, "CCutil_gettsplib failed\n"); goto CLEANUP;
        }
    } else {
        rval = CCutil_getdata (edgefname, 0, CC_SPARSE, &ncount, &dat,
                               0, 0, &rstate);
        if (rval) {
            fprintf (stderr, "CCutil_getdata failed\n"); goto CLEANUP;
        }
    }

    tour = CC_SAFE_MALLOC (ncount, int);
    if (!tour) {
        fprintf (stderr, "out of memory in main\n");
        rval = 1; goto CLEANUP;
    }

    if (eformat == 0) {
        if (simpletour) {
            rval = CCutil_getcycle (ncount, cycfname, tour, 0);
            if (rval) {
                fprintf (stderr, "CCutil_getcycle failed\n"); goto CLEANUP;
            }
        } else {
            rval = CCutil_getcycle_tsplib (ncount, cycfname, tour);
            CCcheck_rval (rval, "CCutil_getcycle_tsplib failed");
        }
    } else {
        rval = CCutil_getcycle_edgelist (ncount, cycfname, tour, 0);
        if (rval) {
            fprintf (stderr, "CCutil_getcycle_edgelist failed\n");
            goto CLEANUP;
        }
    }

    {
        int *chk = (int *) NULL;
        int i;

        chk = CC_SAFE_MALLOC (ncount, int);
        CCcheck_NULL (chk, "out of memory in main");

        for (i = 0; i < ncount; i++) chk[i] = 0;
        for (i = 0; i < ncount; i++) {
            if (chk[tour[i]] == 1) {
                fprintf (stderr, "duplicate node in tour: %d, position %d\n",
                                  tour[i], i);
                rval = 1;  goto CLEANUP;
            }
            chk[tour[i]] = 1;
        }
        CC_IFFREE (chk, int);
    }

    if (edgefname) {
        int istour;
        rval = CCutil_sparse_real_tour (ncount, &dat, tour, &istour);
        if (rval) {
            fprintf (stderr, "CCutil_sparse_real_tour failed\n");
            goto CLEANUP;
        }
        if (istour == 0) {
            printf ("Tour is not contained in the sparse edge set\n");
            fflush (stdout);
            goto CLEANUP;
        }
    }

    CCutil_cycle_len (ncount, &dat, tour, &val);
    printf ("Tour Length: %.0f\n", val); fflush (stdout);


CLEANUP:

    CC_IFFREE (tour, int);
    CCutil_freedatagroup (&dat);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "eE:tT:", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'e':
            eformat = 1;
            break;
        case 'E':
            edgefname = boptarg;
            break;
        case 't':
            simpletour = 1;
            break;
        case 'T':
            tspfname = boptarg;
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac) {
        cycfname = av[boptind++];
    } else {
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] tour_file\n", fname);
    fprintf (stderr, "   -e    tour file is an edge list (default is TSPLIB)\n");
    fprintf (stderr, "   -t    tour file in concorde permutation format (default TSPLIB)\n"); 
    fprintf (stderr, "   -E f  edge file to specify lengths\n");
    fprintf (stderr, "   -T f  TSPLIB file to specify lengths\n");
}

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
/*            UTIL PROGRAM TO CHECK A CYCLE FILE IN EDGE FORMAT             */
/*                                                                          */
/*                              TSP CODE                                    */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: May 22, 1995                                                      */
/*                                                                          */
/*  SEE short decsribtion in usage ().                                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

static char *cyclefilename = (char *) NULL;
static char *outfname = (char *) NULL;
static int dummy;


int
    main (int ac, char **av);
static void
    usage (char *f);
static int
    parseargs (int ac, char **av);



int main (int ac, char **av)
{
    double szeit = CCutil_zeit ();
    int ncount = 0, rval = 0;
    int *tour = (int *) NULL;
    double len = 0.0;
    FILE *out = (FILE *) NULL;

    if (parseargs (ac, av))
        return 1;

    {
        int i, dum, w;
        FILE *in = fopen (cyclefilename, "r");

        if (in == (FILE *) NULL) {
            perror (cyclefilename);
            fprintf (stderr, "Unable to open %s for input\n", cyclefilename);
            rval = 1; goto CLEANUP;
        }
        ncount = CCutil_readint (in);
        i = CCutil_readint (in);
        printf ("Number of nodes: %d\n", ncount);
        fflush (stdout);
        if (i != ncount) {
            fprintf (stderr, "Not an edge-cycle file\n");
            fclose (in);
            rval = 1; goto CLEANUP;
        }
        for (i = 0; i < ncount; i++) {
            dum = CCutil_readint (in);
            dum = CCutil_readint (in);
            w = CCutil_readint (in);
            len += (double) w;
        }
        printf ("Tour Length: %.0f\n", len);
        fflush (stdout);
        fclose (in);
    }

    tour = CC_SAFE_MALLOC (ncount, int);
    if (!tour) {
        fprintf (stderr, "out of memory in main\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_getcycle_edgelist (ncount, cyclefilename, tour, 0);
    if (rval) {
        fprintf (stderr, "CCutil_getcycle_edgelist failed\n");
        goto CLEANUP;
    }


    if (outfname) {
        int i;

        out = fopen (outfname, "w");
        if (!out) {
            perror (cyclefilename);
            fprintf (stderr, "Unable to open %s for input\n", cyclefilename);
            rval = 1; goto CLEANUP;
        }

        fprintf (out, "%d\n", ncount);
        for (i = 0; i < ncount; i++) {
            fprintf (out, "%d\n", tour[i]);
        }
    }

    printf ("Total Running Time: %.2f (seconds)\n", CCutil_zeit () - szeit);
    fflush (stdout);

CLEANUP:

    CC_IFFREE (tour, int);
    if (out) fclose (out);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "do:", &boptind, &boptarg)) != EOF)
        switch (c) {
        case 'd':
            dummy = 0;
            break;
        case 'o':
            outfname = boptarg;
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
    cyclefilename = av[boptind++];

    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s edge_cycle_file\n", f);
    fprintf (stderr,  "    -o f output file for the tour\n"); 
}

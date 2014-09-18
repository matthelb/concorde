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
/*                Demo of Comb Separation Algorithms                        */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal,  and Cook                       */
/*  Date: September 29, 1997                                                */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "combs.h"

static int binary_in = 0;
static char *fname = (char *) NULL;
static int blocktest = 0;


int
    main (int ac, char **av);

static void
    usage (char *f);

static int
    test_blocks (int ncount, int ecount, int *elist, double *ecap),
    parseargs (int ac, char **av);


int main (int ac, char **av)
{
    int rval = 0;
    double szeit;
    int ncount, ecount;
    int *elist = (int *) NULL;
    double *ecap = (double *) NULL;

    if (parseargs (ac, av)) goto CLEANUP;

    rval = CCutil_getedges_double (&ncount, fname, &ecount, &elist, &ecap,
                                   binary_in);
    if (rval) {
        fprintf (stderr, "CCutil_getedges_double failed\n"); goto CLEANUP;
    }

    szeit = CCutil_zeit ();

    if (blocktest) {
        rval = test_blocks (ncount, ecount, elist, ecap);
        if (rval) {
            fprintf (stderr, "test_blocks failed\n");
            goto CLEANUP;
        }
    } 

    printf ("Total Running Time: %.2f (seconds)\n", CCutil_zeit () - szeit);
    fflush (stdout);

CLEANUP:

    CC_IFFREE (elist, int);
    CC_IFFREE (ecap, double);
    return rval;
}

static int test_blocks (int ncount, int ecount, int *elist, double *ecap)
{
    int rval = 0;
    int *cutnodes = (int *) NULL;
    int ncutnodes = 0;
    int **blocks = (int **) NULL;
    int *blockcnt = (int *) NULL;
    int nblocks = 0;
    int i, j;

    rval = CCcombs_find_blocks (ncount, ecount, elist, ecap,
             &nblocks, &blockcnt, &blocks, &ncutnodes, &cutnodes);
    if (rval) {
        fprintf (stderr, "CCcombs_find_blocks failed\n"); goto CLEANUP;
    }

    printf ("Cutnodes: ");
    for (i = 0; i < ncutnodes; i++) {
        printf ("%d ", cutnodes[i]);
    }
    printf ("\n");
    fflush (stdout);

    for (i = 0; i < nblocks; i++) {
        printf ("Block %d: ", i); fflush (stdout);
        for (j = 0; j < blockcnt[i]; j++) {
            printf ("%d ", blocks[i][j]);
        }
        printf ("\n"); fflush (stdout);
    }

CLEANUP:

    CC_IFFREE (cutnodes, int);
    for (i = 0; i < nblocks; i++) {
        CC_IFFREE (blocks[i], int);
    }
    CC_IFFREE (blocks, int *);
    CC_IFFREE (blockcnt, int);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "bB", &boptind, &boptarg)) != EOF)
        switch (c) {
        case 'b':
            binary_in = 1;
            break;
        case 'B':
            blocktest = 1;
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

    fname = av[boptind++];
    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "usage: %s [- below -] edge_file\n", f);
    fprintf (stderr, "    b:   binary input file\n");
    fprintf (stderr, "    B:   find block decomposition\n");
}

#if 0
static int display_all_cuts (double val, int cnt, int *cut, void *pass_param)
{
    if (pass_param) {
        fprintf (stderr, "don't know about pass_param in display_all_cuts\n");
        return 1;
    }

    if (cut && cnt) {
        printf ("Found cut of value %f\n", val); fflush (stdout);
    }
    return 0;
}
#endif

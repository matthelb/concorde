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
/*           A PROGRAM TO CONVERT EDGE/DAT FILES TO TSPLIB FILES            */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: January 28, 1998                                                  */
/*                                                                          */
/*  SEE short decsribtion in usage ().                                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

static char *fname = (char *) NULL;
static int norm = CC_EUCLIDEAN;
static int seed = 0;
static int binary_in = 0;
static int isdat = 0;
static int nnodes_want = 0;


int
    main (int ac, char **av);

static int
    parseargs (int ac, char **av);

static void
    usage (char *f);


int main (int ac, char **av)
{
    int rval = 0;
    int ncount, ecount;
    int *elist = (int *) NULL;
    int *elen = (int *) NULL;
    CCdatagroup dat;
    CCrandstate rstate;

    CCutil_init_datagroup (&dat);

    seed = (int) CCutil_real_zeit ();

    rval = parseargs (ac, av);
    if (rval) return 1;

    if (!nnodes_want && !fname) {
        usage (av[0]);
        return 1;
    }

    CCutil_sprand (seed, &rstate);
    printf ("Using random seed %d\n", seed); fflush (stdout);

    if (isdat || !fname) {
        ncount = nnodes_want;
        rval = CCutil_getdata (fname, binary_in, norm, &ncount, &dat, 0, 1, &rstate);
        if (rval) {
            fprintf (stderr, "CCutil_getdata failed\n"); goto CLEANUP;
        }
    } else {
        rval = CCutil_getedgelist_n (&ncount, fname, &ecount, &elist,
                                     &elen, 0);
        if (rval) {
            fprintf (stderr, "CCutil_getedgelist failed\n"); goto CLEANUP;
        }
        rval = CCutil_graph2dat_matrix (ncount, ecount, elist, elen, 1000000,
                                        &dat);
        if (rval) {
            fprintf (stderr, "CCutil_graph2dat_matrix failed\n"); goto CLEANUP;
        }
    }

    rval = CCutil_writetsplib ("out.tsp", ncount, &dat);
    if (rval) {
        fprintf (stderr, "CCutil_writetsplib failed\n"); goto CLEANUP;
    }

CLEANUP:

    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    CCutil_freedatagroup (&dat);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int inorm, c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "bdk:N:s:", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'b':
            binary_in = 1;
            break;
        case 'd':
            isdat = 1;
            break;
        case 'k':
            nnodes_want = atoi (boptarg);
            break;
        case 's':
            seed = atoi (boptarg);
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
                printf ("unknown norm %d\n", inorm); fflush (stdout);
                usage (av[0]);
                return 1;
            }
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac) {
        fname = av[boptind++];
    }
    if (boptind < ac) {
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s [-flags below] [filename]\n", f);
    fprintf (stderr, "   -b    dat file is in binary\n");
    fprintf (stderr, "   -d    dat file (default is edge file)\n");
    fprintf (stderr, "   -k #  number of nodes for random problem\n");
    fprintf (stderr, "   -s #  random seed\n");
    fprintf (stderr, "   -N #  norm (default is L2)\n");
    fprintf (stderr, "         0=MAX, 1=L1, 2=L2, 3=3D, 4=USER, 5=ATT, 6=GEO, 7=MATRIX,\n");
    fprintf (stderr, "         8=DSJRAND, 9=CRYSTAL, 10=SPARSE, 11-15=RH-norm 1-5, 16=TOROIDAL\n");
    fprintf (stderr, "         17=GEOM, 18=JOHNSON\n");
}

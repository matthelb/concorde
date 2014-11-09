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
/*               CODE FOR TESTING linkern_path                              */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: May 6, 2003                                                       */
/*                                                                          */
/*  For a short describtion see usage ()                                    */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "linkern.h"
#include "util.h"
#include "kdtree.h"
#include "edgegen.h"
#include "macrorus.h"

static int norm = CC_EUCLIDEAN;
static int seed = 0;
static int binary_in = 0;
static int tsplib_in = 1;
static int nnodes_want = 0;
static int gridsize = 0;
static int quadtry = 2;
static int run_silently = 0;

static int in_repeater = -1;

static char *nodefile = (char *) NULL;
static char *goodfname = (char *) NULL;
static char *cycfname = (char *) NULL;
static char *outfname = (char *) NULL;


int
   main (int, char **);

static int
   print_command (int ac, char **av),
   parseargs (int, char **);

static void
   usage (char *f);


int main (int ac, char **av)
{
    int k, ncount;
    double val, best;
    double startzeit;
    int tempcount, *templist;
    int *incycle = (int *) NULL, *outcycle = (int *) NULL;
    CCdatagroup dat;
    int rval = 0;
    CCrandstate rstate;
    int allow_dups;
    int use_gridsize;

    CCutil_printlabel ();
    CCutil_init_datagroup (&dat);

    rval = print_command (ac, av);
    CCcheck_rval (rval, "print_command failed");

    seed = (int) CCutil_real_zeit ();
    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    CCutil_sprand (seed, &rstate);

    printf ("Chained Lin-Kernighan with seed %d\n", seed);
    fflush (stdout);

    if ((!nnodes_want && !nodefile) || (tsplib_in && !nodefile)) {
        usage (av[0]);
        goto CLEANUP;
    }

    startzeit = CCutil_zeit ();

    if (tsplib_in) {
        if (CCutil_gettsplib (nodefile, &ncount, &dat)) {
            fprintf (stderr, "could not read the TSPLIB file\n");
            rval = 1;
            goto CLEANUP;
        }
        CCutil_dat_getnorm (&dat, &norm);
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
        if (CCutil_getdata (nodefile, binary_in, norm, &ncount, &dat,
                            use_gridsize, allow_dups, &rstate)) {
            rval = 1;
            goto CLEANUP;
        }
    }

    if (in_repeater == -1) in_repeater = ncount / 10;


    if (cycfname) {
        incycle = CC_SAFE_MALLOC (ncount, int);
        CCcheck_NULL (incycle, "out of memory in main");

        rval =  CCutil_getcycle (ncount, cycfname, incycle, 0);
        CCcheck_rval (rval, "CCutil_getcycle failed");
    }

    if (goodfname) {
        int *templen = (int *) NULL;
        rval = CCutil_getedgelist (ncount, goodfname, &tempcount, &templist,
                                   &templen, 0);
        CCcheck_rval (rval, "CCutil_getedgelist failed");

        CC_IFFREE (templen, int);
        printf ("Read good-edge file: %d edges\n", tempcount);
        fflush (stdout);
    } else {
        CCedgegengroup plan;

        CCedgegen_init_edgegengroup (&plan);
        plan.quadnearest = quadtry;

        rval = CCedgegen_edges (&plan, ncount, &dat, (double *) NULL,
                                &tempcount, &templist, 0, &rstate);
        CCcheck_rval (rval, "CCedgegen_edges failed");
    } 


    outcycle = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (outcycle, "out of memory in main");

    rval = CClinkern_path (ncount, &dat, tempcount, templist, 
               in_repeater, incycle, outcycle, &val, run_silently, &rstate);
    CCcheck_rval (rval, "CClinkern_path failed");

    printf ("First: %d   Last: %d\n", outcycle[0], outcycle[ncount-1]);
    printf ("Total Running Time: %.2f\n", CCutil_zeit () - startzeit);
    fflush (stdout);

    if (outfname) {
        rval = CCutil_writecycle (ncount, outfname, outcycle, 0);
        CCcheck_rval (rval, "CCutil_writecycle failed");
    }

CLEANUP:

    CC_IFFREE (templist, int);
    CC_IFFREE (incycle, int);
    CC_IFFREE (outcycle, int);
    CCutil_freedatagroup (&dat);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c, k;
    int inorm;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "bBg:G:k:lN:o:q:QR:s:y:", &boptind, &boptarg)) != EOF)
        switch (c) {
        case 'B':
            binary_in = 1;
            break;
        case 'b':
            binary_in = 2;
            break;
        case 'g':
            goodfname = boptarg;
            break;
        case 'G':
            gridsize = atoi(boptarg);
            break;
        case 'k':
            nnodes_want = atoi (boptarg);
            tsplib_in = 0;
            break;
        case 'N':
            inorm = atoi(boptarg);
            switch (inorm) {
            case 0: norm = CC_MAXNORM; break;
            case 1: norm = CC_EUCLIDEAN_CEIL; break;
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
            case 18: norm = CC_MANNORM; break;
            default:
                usage (av[0]);
                return 1;
            }
            tsplib_in = 0;
            break;
        case 'o':
            outfname = boptarg;
            break;
        case 'q':
            quadtry = atoi (boptarg);
            break;
        case 'Q':
            run_silently++;
            break;
        case 'R':
            in_repeater = atoi (boptarg);
            break;
        case 's':
            seed = atoi (boptarg);
            break;
        case 'y':
            cycfname = boptarg;
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    if (boptind < ac)
        nodefile = av[boptind++];

    if (boptind > ac) {
        usage (av[0]);
        return 1;
    }
    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "usage: %s [- see below -] [tsplib_file or dat_file]\n", f);
    fprintf (stderr, "   -s #  random number seed\n");
    fprintf (stderr, "   -k #  number of nodes for random problem\n");
    fprintf (stderr, "   -G #  use #x# grid for random points, no dups if #<0\n");
    fprintf (stderr, "   -o f  save final tour\n");
    fprintf (stderr, "   -q #  use quad #-nearest as the sparse set (default is 2)\n");
    fprintf (stderr, "   -g f  use the edges in file f as the sparse edge set\n");
    fprintf (stderr, "   -R #  number of kicks in iterated Lin-Kernighan (default is #nodes)\n");
    fprintf (stderr, "   -y f  starting cycle\n");
    fprintf (stderr, "   -Q    run silently\n");
    fprintf (stderr, "   -b    dat file in binary doubles\n");
    fprintf (stderr, "   -B    dat file in binary ints\n");
    fprintf (stderr, "   -N #  norm (must specify if dat file is not a TSPLIB file)\n");
    fprintf (stderr, "         0=MAX, 1=JOHNSON, 2=L2, 3=3D, 4=USER, 5=ATT, 6=GEO, 7=MATRIX,\n");
    fprintf (stderr, "         8=DSJRAND, 9=CRYSTAL, 10=SPARSE, 11-15=RH-norm 1-5, 16=TOROIDAL\n");
    fprintf (stderr, "         17=GEOM, 18=L1\n");
}

static int print_command (int ac, char **av)
{
    int rval = 0;
    int i, cmdlen = 0;
    char *cmdout = (char *) NULL;

    for (i=0; i<ac; i++) {
        cmdlen += strlen(av[i]) + 1;
    }
    cmdout = CC_SAFE_MALLOC (cmdlen, char);
    CCcheck_NULL (cmdout, "out of memory in print_command");

    cmdlen = 0;
    for (i=0; i<ac; i++) {
        strcpy (cmdout + cmdlen, av[i]);
        cmdlen += strlen(av[i]);
        cmdout[cmdlen] = ' ';
        cmdlen++;
    }
    cmdout[cmdlen-1] = '\0';
    printf ("%s\n", cmdout); fflush (stdout);

CLEANUP:

    CC_IFFREE (cmdout, char);
    return rval;
}

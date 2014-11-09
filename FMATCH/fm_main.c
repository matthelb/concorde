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
/*                TEST PROGRAM FOR FRACTIONAL MATCHINGS                     */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 24, 1995  (cofeb24)                                      */
/*                                                                          */
/*  SEE short decsription in usage ().                                      */
/*                                                                          */
/*  Link with:                                                              */
/*    SEE fmmake.grs                                                        */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "fmatch.h"
#include "util.h"
#include "kdtree.h"
#include "edgegen.h"

static int norm = CC_EUCLIDEAN;
static int wantbasic = 0;
static char *edgefilename = (char *) NULL;
static char *datfilename = (char *) NULL;
static char *edgegenfname = (char *) NULL;
static int dumpmatch = 0;
static int dumpdual = 0;
static int dumpbasis = 0;
static int binary_in = 0;
static int tsplib_in = 1;
static int nnodes_want = 0;
static int gridsize = 0;
static int seed = 0;
static int usenn2match = 0;
static int quadtry = 2;
static int run_silently = 0;



int
    main (int ac, char **av);

static void
    usage (char *f),
    dump_match (int *thematching),
    dump_dual (int *thedual, int ncount),
    dump_basis (int *thebasis, int ncount);

static int
    parseargs (int ac, char **av),
    getgraph (char *edgefile, CCdatagroup *dat, int *ncount, int *ecount,
        int **elist, int **elen, int silent, CCrandstate *rstate);



int main (int ac, char **av)
{
    double v, szeit;
    int ncount, ecount;
    int *elist = (int *) NULL, *elen = (int *) NULL;
    int *thematching = (int *) NULL, *thedual = (int *) NULL;
    int *thebasis = (int *) NULL;
    CCdatagroup dat, *mydat = (CCdatagroup *) NULL;
    int rval = 0;
    CCrandstate rstate;
    int allow_dups;
    int use_gridsize;

    CCutil_init_datagroup (&dat);
    
    seed = (int) CCutil_real_zeit ();
    if (parseargs (ac, av))
        return 0;
    CCutil_sprand (seed, &rstate);

    if (edgefilename == (char *) NULL && datfilename == (char *) NULL &&
                                         nnodes_want == 0) {
        usage (av[0]);
        return 0;
    }
    ncount = nnodes_want;

    szeit = CCutil_zeit ();
    if (tsplib_in && datfilename != (char *) NULL) {
        rval = CCutil_gettsplib (datfilename, &ncount, &dat);
        if (rval) {
            fprintf (stderr, "could not read the TSPLIB file\n"); goto CLEANUP;
        }
        CCutil_dat_getnorm (&dat, &norm);
        mydat = &dat;
    } else if (edgefilename == (char *) NULL || datfilename != (char *) NULL) {
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
        rval = CCutil_getdata (datfilename, binary_in, norm, &ncount, &dat,
                               use_gridsize, allow_dups, &rstate);
        if (rval) {
            fprintf (stderr, "Could not create data set\n"); goto CLEANUP;
        }
        mydat = &dat;
    } else {
        mydat = (CCdatagroup *) NULL;
    }

    rval = getgraph (edgefilename, &dat, &ncount, &ecount, &elist, &elen,
                     run_silently, &rstate);
    if (rval) {
        fprintf (stderr, "getgraph failed\n"); goto CLEANUP;
    }
    printf ("Initial edgeset: %d edges (%d nodes)\n", ecount, ncount);
    printf ("Time to generate graph: %.2f (seconds)\n",
            CCutil_zeit () - szeit);
    fflush (stdout);

    if (dumpmatch) {
        thematching = CC_SAFE_MALLOC((6 * ncount) + 1, int);
        if (!thematching) {
            fprintf (stderr, "out of memory in main\n");
            rval = 1; goto CLEANUP;
        }
    }
    if (dumpdual) {
        thedual = CC_SAFE_MALLOC(ncount, int);
        if (!thedual) {
            fprintf (stderr, "out of memory in main\n");
            rval = 1; goto CLEANUP;
        }
    }
    if (dumpbasis) {
        thebasis = CC_SAFE_MALLOC(2 * ncount, int);
        if (!thebasis) {
            fprintf (stderr, "out of memory in main\n");
            rval = 1; goto CLEANUP;
        }
    }


    szeit = CCutil_zeit ();
    rval = CCfmatch_fractional_2match (ncount, ecount, elist, elen, mydat,
                   &v, thematching, thedual, thebasis, wantbasic, 
                   run_silently, &rstate);
    if (rval) {
        fprintf (stderr, "Fractional matching routine failed\n");
        goto CLEANUP;
    }

    printf ("Total Running Time: %.2f (seconds)\n", CCutil_zeit () - szeit);
    printf ("Final matching weight: %.1f\n", v);
    fflush (stdout);

    if (dumpmatch)
        dump_match (thematching);
    if (dumpdual)
        dump_dual (thedual, ncount);
    if (dumpbasis)
        dump_basis (thebasis, ncount);


CLEANUP:

    CC_IFFREE(thematching, int);
    CC_IFFREE(thedual, int);
    CC_IFFREE(thebasis, int);
    CC_IFFREE(elist, int);
    CC_IFFREE(elen, int);

    if (mydat != (CCdatagroup *) NULL) {
        CCutil_freedatagroup (mydat);
    }

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c, inorm;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "BbD:e:k:mn:N:q:Qr:s:xyz", &boptind, &boptarg)) != EOF)
        switch (c) {
        case 'B':
            wantbasic = 1;
            break;
        case 'b':
            binary_in = 1;
            break;
        case 'D':
            edgegenfname = boptarg;
            break;
        case 'e':
            edgefilename = boptarg;
            break;
        case 'k':
            nnodes_want = atoi (boptarg);
            break;
        case 'm':
            usenn2match = 1;
             break;
        case 'n':
            datfilename = boptarg;
            break;
        case 'q':
            quadtry = atoi (boptarg);
            break;
        case 'Q':
            run_silently = 1;
            break;
        case 'r':
            gridsize = atoi(boptarg);
            break;
        case 's':
            seed = atoi (boptarg);
            break;
        case 'x':
            dumpmatch = 1;
            break;
        case 'y':
            dumpdual = 1;
            break;
        case 'z':
            dumpbasis = 1;
            wantbasic = 1;
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
    if (boptind != ac) {
        usage (av[0]);
        return 1;
    }
    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s [-see below-]\n", f);
    fprintf (stderr, "   -B    find basic optimal solution\n");
    fprintf (stderr, "   -b    datfile in integer binary format\n");
    fprintf (stderr, "   -D f  edgegen file for initial edge set\n");
    fprintf (stderr, "   -e f  edge file - initial edge set\n");
    fprintf (stderr, "   -k #  number of nodes for random problem\n");
    fprintf (stderr, "   -m    NN 2-matching as initial edge set\n");
    fprintf (stderr, "   -n f  dat file - for fmatch on complete graph\n");
    fprintf (stderr, "   -q #  quad-nearest # as initial edge set (default 2)\n");
    fprintf (stderr, "   -Q    run quietly (don't generate so much output)\n");
    fprintf (stderr, "   -r #  use #x# grid for random points, no dups if #<0\n");
    fprintf (stderr, "   -s #  random seed\n");
    fprintf (stderr, "   -x    dump matching to match.out\n");
    fprintf (stderr, "   -y    dump dual solution to dual.out\n");
    fprintf (stderr, "   -z    dump basic edges to basis.out\n");
    fprintf (stderr, "   -N #  norm for pricing (must specify if dat file is not a TSPLIB file)\n");
    fprintf (stderr, "         0=MAX, 1=L1, 2=L2, 3=3D, 4=USER, 5=ATT, 6=GEO, 7=MATRIX,\n");
    fprintf (stderr, "         8=DSJRAND, 9=CRYSTAL, 11-15=RH-norm 1-5, 16=TOROIDAL\n");
    fprintf (stderr, "         17=GEOM, 18=JOHNSON\n");
}

#define QUAD_TRY 2

static int getgraph (char *edgefile, CCdatagroup *dat, int *ncount,
        int *ecount, int **elist, int **elen, int silent, CCrandstate *rstate)
{
    FILE *in;
    int i, k;
    int datnorm;

    *elist = (int *) NULL;
    *elen = (int *) NULL;

    CCutil_dat_getnorm (dat, &datnorm);
    if (edgefile != (char *) NULL) {
        if ((in = fopen (edgefile, "r")) == (FILE *) NULL) {
            perror (edgefile);
            fprintf (stderr, "Unable to open %s for input\n", edgefile);
            return 1;
        }

        k = CCutil_readint (in);
        if (*ncount != 0 && k != *ncount) {
            fprintf (stderr, "Edge file does not match dat file\n");
            fclose (in);
            return 1;
        }
        *ncount = k;
        *ecount = CCutil_readint (in);

        *elist = CC_SAFE_MALLOC(2 * (*ecount), int);
        if (!(*elist)) {
            fclose (in);
            return 1;
        }
        *elen = CC_SAFE_MALLOC(*ecount, int);
        if (!(*elen)) {
            fclose (in);
            CC_FREE (*elist, int);
            return 1;
        }

        for (i = 0, k = 0; i < *ecount; i++) {
            (*elist)[k++] = CCutil_readint (in);
            (*elist)[k++] = CCutil_readint (in);
            (*elen)[i] = CCutil_readint (in);
        }
        fclose (in);
    } else if (edgegenfname) {
        CCedgegengroup plan;
        if (CCedgegen_read (edgegenfname, &plan)) {
            fprintf (stderr, "CCedgegen_read failed\n");
            return 1;
        }
        if (CCedgegen_edges (&plan, *ncount, dat, (double *) NULL, ecount,
                             elist, silent, rstate)) {
            fprintf (stderr, "CCedgegen_edges failed\n");
            return 1;
        }
        *elen = CC_SAFE_MALLOC(*ecount, int);
        if (!(*elen)) {
            CC_FREE (*elist, int);
            return 1;
        }
        for (i = 0; i < *ecount; i++) {
            (*elen)[i] = CCutil_dat_edgelen ((*elist)[2*i],
                    (*elist)[(2*i) + 1], dat);
        }
    } else {
        if (usenn2match) {
            double val;
            *ecount = *ncount;
            *elist = CC_SAFE_MALLOC(2 * (*ecount), int);
            if (!(*elist))
                return 1;
            if ((datnorm & CC_NORM_BITS) == CC_KD_NORM_TYPE) {
                printf ("Using nearest neighbor 2-matching graph\n");
                fflush (stdout);
                if (CCkdtree_nearest_neighbor_2match ((CCkdtree *) NULL,
                        *ncount, CCutil_lprand (rstate) % (*ncount), dat,
                        *elist, &val, rstate)) {
                   fprintf (stderr, "nearest 2-matching code failed\n");
                   CC_FREE (*elist, int);
                   return 1;
                }
            } else {
                int *cyc = (int *) NULL;
                printf ("Not setup for nearest 2-match with x or junk norms\n");
                printf ("Using nearest neighbour tour graph\n");
                fflush (stdout);
                cyc = CC_SAFE_MALLOC (*ncount, int);
                if (!cyc) {
                    CC_FREE (*elist, int);
                    return 1;
                }
                if ((datnorm & CC_NORM_BITS) == CC_X_NORM_TYPE) {
                    if (CCedgegen_x_nearest_neighbor_tour (*ncount,
                        CCutil_lprand (rstate) % (*ncount), dat, cyc, &val)) {
                        CC_FREE (*elist, int);
                        CC_FREE (cyc, int);
                        return 1;
                    }
                } else {
                    if (CCedgegen_junk_nearest_neighbor_tour (*ncount,
                        CCutil_lprand (rstate) % (*ncount), dat, cyc,
                        &val, silent)) {
                        CC_FREE (*elist, int);
                        CC_FREE (cyc, int);
                        return 1;
                    }
                }
                for (i = 0; i < *ncount - 1; i++) {
                    (*elist)[2 * i] = cyc[i];
                    (*elist)[(2 * i) + 1] = cyc[i + 1];
                }
                (*elist)[2 * i] = cyc[i];
                (*elist)[(2 * i) + 1] = cyc[0];
                CC_FREE (cyc, int);
            }
        } else {
            if ((datnorm & CC_NORM_BITS) == CC_KD_NORM_TYPE) {
                printf ("Using quadrant nearest %d graph\n", quadtry);
                fflush (stdout);
                if (CCkdtree_quadrant_k_nearest ((CCkdtree *) NULL, *ncount,
                                      quadtry, dat, (double *) NULL, 1, 
                                      ecount, elist, silent, rstate)) {
                    fprintf (stderr, "CCkdtree-quad nearest code failed\n");
                    CC_FREE (*elist, int);
                    return 1;
                }
            } else if ((datnorm & CC_NORM_BITS) == CC_X_NORM_TYPE) {
                printf ("Using quadrant nearest %d graph\n", quadtry);
                fflush (stdout);
                if (CCedgegen_x_quadrant_k_nearest (*ncount, quadtry, dat,
                       (double *) NULL, 1, ecount, elist, run_silently)) {
                    fprintf (stderr, "CCedgegen_x_quadrant_k_nearest code failed\n");
                    CC_FREE (*elist, int);
                    return 1;
                }
            } else {
                printf ("No junk quad nearest, using %d nearest graph\n",
                         4 * quadtry);
                fflush (stdout);
                if (CCedgegen_junk_k_nearest (*ncount, 4 * quadtry, dat,
                       (double *) NULL, 1, ecount, elist, run_silently)) {
                    fprintf (stderr, "CCedgegen_junk_k_nearest code failed\n");
                    CC_FREE (*elist, int);
                    return 1;
                }
            }
        }

        *elen = CC_SAFE_MALLOC(*ecount, int);
        if (!(*elen)) {
            CC_FREE (*elist, int);
            return 1;
        }
        for (i = 0; i < *ecount; i++) {
            (*elen)[i] = CCutil_dat_edgelen ((*elist)[2*i],
                    (*elist)[(2*i) + 1], dat);
        }
    }

    return 0;
}

static void dump_match (int *thematching)
{
    FILE *out = fopen ("match.out", "w");
    int k;

    if (out == (FILE *) NULL) {
        perror ("match.out");
        fprintf (stderr, "Unable to open match.out for output\n");
        return;
    }

    k = 0;
    while (thematching[k] != -1) {
        fprintf (out, "%d %d %d\n", thematching[k], thematching[k + 1],
                       thematching[k + 2]);
        k += 3;
    }
}
static void dump_dual (int *thedual, int ncount)
{
    FILE *out = fopen ("dual.out", "w");
    int k;

    if (out == (FILE *) NULL) {
        perror ("dual.out");
        fprintf (stderr, "Unable to open dual.out for output\n");
        return;
    }

    for (k = 0; k < ncount; k++)
        fprintf (out, "%d\n", thedual[k]);
}


static void dump_basis (int *thebasis, int ncount)
{
    FILE *out = fopen ("basis.out", "w");
    int k;

    if (out == (FILE *) NULL) {
        perror ("basis.out");
        fprintf (stderr, "Unable to open basis.out for output\n");
        return;
    }

    for (k = 0; k < ncount; k++)
        fprintf (out, "%d %d\n", thebasis[2 * k], thebasis[(2 * k) + 1]);
}

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
/*                  CODE FOR TESTING KD-TREE ROUTINES                       */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 24, 1995  (cofeb24)                                      */
/*                                                                          */
/*  For a short describtion see usage ()                                    */
/*                                                                          */
/****************************************************************************/


#include "machdefs.h"
#include "util.h"
#include "kdtree.h"

static int norm = CC_EUCLIDEAN;
static int usenodeweights = 0;
static int random_weight_limit = 0;
static int seed = 0;
static int nnodes_want = 0;
static int gridsize = 0;
static int nearnum = 0;
static int quadnearnum = 0;
static int find_nearest_tour = 0;
static int find_nearest_2match = 0;
static int find_greedy_tour = 0;
static int find_fa_tour = 0;
static int find_qboruvka_tour = 0;
static int find_boruvka_tour = 0;
static int find_twoopt_tour = 0;
static int find_3opt_tour = 0;
static int find_spanning_tree = 0;
static int run_two_and_a_half_opt = 0;
static int binary_in = 0;
static int tsplib_in = 1;
static int run_silently = 0;

static char *nodefile = (char *) NULL;
static char *weightfile = (char *) NULL;
static char *cycle_for_twoopt = (char *) NULL;
static char *outfile = (char *) NULL;



int
    main (int ac, char **av);
static void
    usage (char *f);
static int
    parseargs (int ac, char **av);


int main (int ac, char **av)
{
    double val, szeit;
    CCkdtree kt;
    CCdatagroup dat;
    double *wcoord = (double *) NULL;
    int ncount;
    int *ttour = (int *) NULL, *tour2 = (int *) NULL;
    int rval = 0;
    int ecount;
    int *elist = (int *) NULL;
    CCrandstate rstate;
    int use_gridsize, allow_dups;

    CCutil_init_datagroup (&dat);
    
    seed = (int) CCutil_real_zeit ();
    if (parseargs (ac, av))
        return 1;
    CCutil_sprand (seed, &rstate);

    if ((!nnodes_want && !nodefile) || (tsplib_in && !nodefile)) {
        usage (av[0]);
        return 1;
    }

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
    if ((norm & CC_NORM_BITS) != CC_KD_NORM_TYPE) {
        fprintf (stderr, "Cannot run CCkdtree with norm %d\n", norm);
        rval = 1;
        goto CLEANUP;
    }

    if (usenodeweights) {
        if (CCutil_getnodeweights (weightfile, ncount, random_weight_limit,
                                   &wcoord, &rstate)) {
            fprintf (stderr, "could not read the nodeweight file\n");
            rval = 1;
            goto CLEANUP;
        }
    }

    if (find_nearest_tour || find_greedy_tour || find_twoopt_tour ||
        find_fa_tour || find_qboruvka_tour || find_boruvka_tour ||
        find_3opt_tour) {
        ttour = CC_SAFE_MALLOC (ncount, int);
        if (!ttour) {
            rval = 1;
            goto CLEANUP;
        }
    }
    if (find_twoopt_tour || find_3opt_tour) {
        tour2 = CC_SAFE_MALLOC (ncount, int);
        if (!tour2) {
            rval = 1;
            goto CLEANUP;
        }
    }

    if (cycle_for_twoopt) {
        if (CCutil_getcycle (ncount, cycle_for_twoopt, ttour, 0)) {
            fprintf (stderr, "Getcycle failed\n");
            rval = 1;
            goto CLEANUP;
        }
    } else if (find_nearest_tour) {
        szeit = CCutil_zeit ();
        if (CCkdtree_nearest_neighbor_tour ((CCkdtree *) NULL, ncount,
               CCutil_lprand (&rstate) % ncount, &dat, ttour, &val, &rstate)) {
            fprintf (stderr, "Nearest neighbor failed\n");
            rval = 1;
            goto CLEANUP;
        }
        printf ("NN tour: %.2f (%.2f seconds)\n", val, CCutil_zeit () - szeit);
        fflush (stdout);
    } else if (find_fa_tour) {
        szeit = CCutil_zeit ();
        if (CCkdtree_far_add_tour ((CCkdtree *) NULL, ncount,
               CCutil_lprand (&rstate) % ncount, &dat, ttour, &val, &rstate)) {
            fprintf (stderr, "Farthest Addition Tour failed\n");
            rval = 1;
            goto CLEANUP;
        }
        printf ("FA tour: %.2f (%.2f seconds)\n", val, CCutil_zeit () - szeit);
        fflush (stdout);
        {
            char *marks;
            int i;

            marks = CC_SAFE_MALLOC (ncount, char);
            if (!marks) {
                rval = 1;
                goto CLEANUP;
            }
            for (i = 0; i < ncount; i++)
                marks[i] = 0;
            for (i = 0; i < ncount; i++) {
                if (ttour[i] < 0 || ttour[i] >= ncount) {
                    fprintf (stderr, "MADE NODE IN FA TOUR: %d\n", ttour[i]);
                    rval = 1;
                    goto CLEANUP;
                }
                if (marks[ttour[i]]) {
                    fprintf (stderr, "REPEAT NODE IN FA-TOUR: %d\n", ttour[i]);
                    fprintf (stderr, "BAD INDEX: %d\n", i);
                    rval = 1;
                    goto CLEANUP;
                } else {
                    marks[ttour[i]] = 1;
                }
            }
            CC_FREE (marks, char);
        }
    } else if (find_boruvka_tour) {
        szeit = CCutil_zeit ();
        if (CCkdtree_boruvka_tour ((CCkdtree *) NULL, ncount, &dat, ttour,
                                   &val, &rstate)) {
            fprintf (stderr, "Boruvka tour failed\n");
            rval = 1;
            goto CLEANUP;
        }
        printf ("Boruvka tour: %.2f (%.2f seconds)\n",
                val, CCutil_zeit () - szeit);
        fflush (stdout);
    } else if (find_qboruvka_tour) {
        szeit = CCutil_zeit ();
        if (CCkdtree_qboruvka_tour ((CCkdtree *) NULL, ncount, &dat, ttour,
                                    &val, &rstate)) {
            fprintf (stderr, "Quick-Boruvka tour failed\n");
            rval = 1;
            goto CLEANUP;
        }
        printf ("Quick-Boruvka tour: %.2f (%.2f seconds)\n",
                         val, CCutil_zeit () - szeit);
        fflush (stdout);
    } else if (find_greedy_tour) {
        szeit = CCutil_zeit ();
        if (CCkdtree_greedy_tour ((CCkdtree *) NULL, ncount, &dat, ttour,
                                  &val, run_silently, &rstate)) {
            fprintf (stderr, "Greedy tour failed\n");
            rval = 1;
            goto CLEANUP;
        }
        printf ("Greedy tour: %.2f (%.2f seconds)\n",
                val, CCutil_zeit () - szeit);
        fflush (stdout);
    }

    if (find_twoopt_tour) {
        szeit = CCutil_zeit ();
        if (CCkdtree_twoopt_tour ((CCkdtree *) NULL, ncount, &dat, ttour,
                       tour2, &val, run_two_and_a_half_opt, 0, &rstate)) {
            fprintf (stderr, "Two-opt failed\n");
            rval = 1;
            goto CLEANUP;
        }
        printf ("2-opt tour: %.2f (%.2f seconds))\n",
                val, CCutil_zeit () - szeit);
        fflush (stdout);
    } else if (find_3opt_tour) {
        szeit = CCutil_zeit ();
        if (CCkdtree_3opt_tour ((CCkdtree *) NULL, ncount, &dat, ttour, tour2,
                                &val, 0, &rstate)) {
            fprintf (stderr, "3-opt failed\n");
            rval = 1;
            goto CLEANUP;
        }
        printf ("3-opt tour: %.2f (%.2f seconds))\n",
                val, CCutil_zeit () - szeit);
        fflush (stdout);
        {
            char *marks;
            int i;

            marks = CC_SAFE_MALLOC (ncount, char);
            if (!marks) {
                rval = 1;
                goto CLEANUP;
            }
            for (i = 0; i < ncount; i++)
                marks[i] = 0;
            for (i = 0; i < ncount; i++) {
                if (tour2[i] < 0 || tour2[i] >= ncount) {
                    fprintf (stderr, "MADE NODE IN TOUR2: %d\n", tour2[i]);
                    rval = 1;
                    goto CLEANUP;
                }
                if (marks[tour2[i]]) {
                    fprintf (stderr, "REPEATED NODE IN TOUR2: %d\n",
                             tour2[i]);
                    rval = 1;
                    goto CLEANUP;
                } else {
                    marks[tour2[i]] = 1;
                }
            }
            CC_FREE (marks, char);
        }
    }

    if (find_spanning_tree) {
        if (outfile) {
            ecount = ncount - 1;
            elist = CC_SAFE_MALLOC (2 * ecount, int);
            if (!elist) {
                rval = 1;
                goto CLEANUP;
            }
        }
        szeit = CCutil_zeit ();
        if (CCkdtree_prim_spanningtree ((CCkdtree *) NULL, ncount, &dat,
                wcoord, elist, &val, &rstate)) {
            fprintf (stderr, "Prim_spanningtree failed\n");
            rval = 1;
            goto CLEANUP;
        }
        printf ("Min spanning tree: %.2f (%.2f seconds)\n", val,
                 CCutil_zeit () - szeit);
        fflush (stdout);
    } else if (nearnum) {
        int wantlist = (outfile ? 1 : 0);
        szeit = CCutil_zeit ();
        if (CCkdtree_k_nearest ((CCkdtree *) NULL, ncount, nearnum, &dat,
              wcoord, wantlist, &ecount, &elist, run_silently, &rstate)) {
            fprintf (stderr, "k-nearest failed\n");
            rval = 1;
            goto CLEANUP;
        }
        printf ("Nearest %d: %.2f (seconds)\n",
                nearnum, CCutil_zeit () - szeit);
        fflush (stdout);
    } else if (quadnearnum) {
        int wantlist = (outfile ? 1 : 0);
        szeit = CCutil_zeit ();
        if (CCkdtree_quadrant_k_nearest ((CCkdtree *) NULL, ncount,
              quadnearnum, &dat, wcoord, wantlist, &ecount, &elist,
              run_silently, &rstate)) {
            fprintf (stderr, "k-nearest failed\n");
            rval = 1;
            goto CLEANUP;
        }
        printf ("Quadrant-Nearest %d: %.2f (seconds)\n", quadnearnum,
                 CCutil_zeit () - szeit);
        fflush (stdout);
    } else if (find_nearest_2match) {
        if (outfile) {
            ecount = ncount;
            elist = CC_SAFE_MALLOC (2 * ecount, int);
            if (!elist) {
                rval = 1;
                goto CLEANUP;
            }
        }
        szeit = CCutil_zeit ();
        if (CCkdtree_nearest_neighbor_2match ((CCkdtree *) NULL, ncount,
               CCutil_lprand (&rstate) % ncount, &dat, elist, &val, &rstate)) {
            fprintf (stderr, "Nearest neighbor 2-matching failed\n");
            rval = 1;
            goto CLEANUP;
        }
        printf ("Nearest 2-matching: %.2f (%.2f seconds)\n", val,
                 CCutil_zeit () - szeit);
        fflush (stdout);
    } else if (!find_nearest_tour && !find_greedy_tour && !find_twoopt_tour &&
               !find_qboruvka_tour && !find_boruvka_tour && !find_fa_tour &&
               !find_3opt_tour) {
        szeit = CCutil_zeit ();
        if (CCkdtree_build (&kt, ncount, &dat, wcoord, &rstate)) {
            fprintf (stderr, "CCkdtree_build failed\n");
            rval = 1;
            goto CLEANUP;
        }
        printf ("Built CCkdtree: %.2f (seconds)\n", CCutil_zeit () - szeit);
        fflush (stdout);
        CCkdtree_free (&kt);
    }

    if (outfile) {
        if (find_twoopt_tour || find_3opt_tour) {
            if (CCutil_writecycle (ncount, outfile, tour2, 0)) {
                fprintf (stderr, "Could not write tour\n");
                rval = 1;
                goto CLEANUP;
            }
        } else if (find_nearest_tour || find_greedy_tour || find_fa_tour ||
                   find_qboruvka_tour || find_boruvka_tour) {
            if (CCutil_writecycle (ncount, outfile, ttour, 0)) {
                fprintf (stderr, "Could not write tour\n");
                rval = 1;
                goto CLEANUP;
            }
        } else if (find_spanning_tree ||  find_nearest_2match || nearnum ||
                                                             quadnearnum) {
            if (CCutil_writeedges (ncount, outfile, ecount, elist, &dat, 0)) {
                fprintf (stderr, "Could not write the edge set\n");
                rval = 1;
                goto CLEANUP;
            }
        }
    }


CLEANUP:

    CCutil_freedatagroup (&dat);
    CC_IFFREE (wcoord, double);
    CC_IFFREE (ttour, int);
    CC_IFFREE (tour2, int);
    CC_IFFREE (elist, int);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int inorm;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "bfgG:hjk:mn:N:o:pq:s:tvw:W:x:Xz:Z?", &boptind, &boptarg)) != EOF)
        switch (c) {
        case 'b':
            binary_in = 1;
            break;
        case 'f':
            find_fa_tour = 1;
            break;
        case 'j':
            find_qboruvka_tour = 1;
            break;
        case 'w':
            usenodeweights = 1;
            weightfile = boptarg;
            break;
        case 'W':
            usenodeweights = 1;
            random_weight_limit = atoi (boptarg);
            break;
        case 'k':
            nnodes_want = atoi (boptarg);
            tsplib_in = 0;
            break;
        case 'G':
            gridsize = atoi(boptarg);
            break;
        case 's':
            seed = atoi (boptarg);
            break;
        case 'n':
            nearnum = atoi (boptarg);
            break;
        case 'm':
            find_nearest_2match++;
            break;
        case 'q':
            quadnearnum = atoi (boptarg);
            break;
        case 't':
            find_nearest_tour++;
            break;
        case 'g':
            find_greedy_tour++;
            break;
        case 'v':
            find_boruvka_tour++;
            break;
        case 'Z':
            find_twoopt_tour++;
            find_greedy_tour++;
            break;
        case 'z':
            cycle_for_twoopt = boptarg;
            find_twoopt_tour++;
            break;
        case 'X':
            find_3opt_tour++;
            find_greedy_tour++;
            break;
        case 'x':
            cycle_for_twoopt = boptarg;
            find_3opt_tour++;
            break;
        case 'h':
            run_two_and_a_half_opt++;
            break;
        case 'o':
            outfile = boptarg;
            break;
        case 'p':
            find_spanning_tree++;
            break;
        case 'N':
            inorm = atoi(boptarg);
            switch (inorm) {
            case 0: norm = CC_MAXNORM; break;
            case 1: norm = CC_MANNORM; break;
            case 2: norm = CC_EUCLIDEAN; break;
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
    fprintf (stderr, "Usage: %s [- see below -] [tsplib_file or dat_file]\n", f);
    fprintf (stderr, "   -b:   dat file in binary-ints\n");
    fprintf (stderr, "   -w f  use node weights from file\n");
    fprintf (stderr, "   -W #  use random node weights (0, #)\n");
    fprintf (stderr, "   -k #  number of nodes for random problem\n");
    fprintf (stderr, "   -G #  use #x# grid for random points, no dups if #<0\n"
);

    fprintf (stderr, "   -s #  random seed\n");
    fprintf (stderr, "   -n #  find # nearest graph\n");
    fprintf (stderr, "   -q #  find quadrant # nearest graph\n");
    fprintf (stderr, "   -t    nearest neighbor tour\n");
    fprintf (stderr, "   -g    greedy tour\n");
    fprintf (stderr, "   -j    quick-boruvka tour\n");
    fprintf (stderr, "   -v    boruvka tour\n");
    fprintf (stderr, "   -f    farthest addition tour\n");
    fprintf (stderr, "   -z f  two_opt the given cycle\n");
    fprintf (stderr, "   -Z    run two_opt (default: on greedy)\n");
    fprintf (stderr, "   -x f  3_opt the given cycle\n");
    fprintf (stderr, "   -X    run 3_opt (default: on greedy)\n");
    fprintf (stderr, "   -h    use limited 3-swaps in two_opt\n");
    fprintf (stderr, "   -m    nearest neighbor 2-matcing\n");
    fprintf (stderr, "   -p    min spanning tree (prim)\n");
    fprintf (stderr, "   -o f  write the cycle or edge set to f\n");
    fprintf (stderr, "   -N #  norm (must specify if dat file is not a TSPLIB file)\n");
    fprintf (stderr, "         0=MAX, 1=L1, 2=L2, 18=JOHNSON\n");
}

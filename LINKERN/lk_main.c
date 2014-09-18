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
/*               CODE FOR TESTING ITERATED LIN-KERNIGHAN                    */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: March 17, 1995                                                    */
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

#define BIGDOUBLE (1e30)

#define LK_RANDOM   (0)
#define LK_NEIGHBOR (1)
#define LK_GREEDY   (2)
#define LK_BORUVKA  (3)
#define LK_QBORUVKA (4)


static int norm = CC_EUCLIDEAN;
static int seed = 0;
static int binary_in = 0;
static int binary_edges = 0;
static int tsplib_in = 1;
static int nnodes_want = 0;
static int gridsize = 0;
static int nearnum = 0;
static int quadtry = 2;
static int run_silently = 0;
static int kick_type = CC_LK_WALK_KICK;
static int tour_type = LK_QBORUVKA;

static int in_repeater = -1;
static int number_runs = 0;
static double time_bound = -1.0;
static double length_bound = -1.0;

static char *nodefile = (char *) NULL;
static char *goodfname = (char *) NULL;
static char *cycfname = (char *) NULL;
static char *edgegenfname = (char *) NULL;
static char *edgecycfname = (char *) NULL;
static char *saveit_final = (char *) NULL;
static char *saveit_name = (char *) NULL;


int
   main (int, char **);
static int
   print_command (int ac, char **av),
   parseargs (int, char **);
static void
   randcycle (int ncount, int *cyc, CCrandstate *rstate),
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
    if (parseargs (ac, av))
        return 1;
    CCutil_sprand (seed, &rstate);

    printf ("Chained Lin-Kernighan with seed %d\n", seed);
    fflush (stdout);

    if ((!nnodes_want && !nodefile) || (tsplib_in && !nodefile)) {
        usage (av[0]);
        return 1;
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

    if (in_repeater == -1) in_repeater = ncount;

    incycle = CC_SAFE_MALLOC (ncount, int);
    if (!incycle) {
        rval = 1;
        goto CLEANUP;
    }
    if (cycfname) {
        if (CCutil_getcycle (ncount, cycfname, incycle, binary_edges)) {
            fprintf (stderr, "CCutil_getcycle failed\n");
            rval = 1;
            goto CLEANUP;
        }
    } else if (edgecycfname) {
        if (CCutil_getcycle_edgelist (ncount, edgecycfname, incycle,
                                      binary_edges)) {
            fprintf (stderr, "CCutil_getcycle_edgelist failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }

    if (goodfname) {
        int *templen = (int *) NULL;
        if (CCutil_getedgelist (ncount, goodfname, &tempcount, &templist,
                                &templen, binary_edges)) {
            rval = 1;
            goto CLEANUP;
        }
        if (templen)
            CC_FREE (templen, int);
        printf ("Read good-edge file: %d edges\n", tempcount);
        fflush (stdout);
    } else if (edgegenfname) {
        CCedgegengroup plan;
        if (CCedgegen_read (edgegenfname, &plan)) {
            fprintf (stderr, "CCedgegen_read failed\n");
            rval = 1;
            goto CLEANUP;
        }
        if (CCedgegen_edges (&plan, ncount, &dat, (double *) NULL, &tempcount,
                     &templist, 0, &rstate)) {
            fprintf (stderr, "CCedgegen_edges failed\n");
            rval = 1;
            goto CLEANUP;
        }
    }

    if ((norm & CC_NORM_BITS) == CC_KD_NORM_TYPE) {
        CCkdtree localkt;
        double kzeit = CCutil_zeit ();

        if ((!goodfname && !edgegenfname) || (!cycfname && !edgecycfname)) {
            if (CCkdtree_build (&localkt, ncount, &dat, (double *) NULL,
                                &rstate)) {
                fprintf (stderr, "CCkdtree_build failed\n");
                rval = 1;
                goto CLEANUP;
            }
            printf ("Time to build kdtree: %.2f\n", CCutil_zeit () - kzeit);
            fflush (stdout);

            if (!goodfname && !edgegenfname) {
                kzeit = CCutil_zeit ();
                if (nearnum) {
                    if (CCkdtree_k_nearest (&localkt, ncount, nearnum, &dat,
                         (double *) NULL, 1, &tempcount, &templist,
                         run_silently, &rstate)) {
                        fprintf (stderr, "CCkdtree_k_nearest failed\n");
                        rval = 1;
                        goto CLEANUP;
                    }
                    if (!run_silently) {
                        printf ("Time to find %d-nearest: %.2f\n", nearnum,
                                                     CCutil_zeit () - kzeit);
                        fflush (stdout);
                    }
                } else {
                    if (CCkdtree_quadrant_k_nearest (&localkt, ncount, quadtry,
                           &dat, (double *) NULL, 1, &tempcount, &templist,
                           run_silently, &rstate)) {
                        fprintf (stderr, "CCkdtree-quad nearest code failed\n");
                        rval = 1;
                        goto CLEANUP;
                    }
                    if (!run_silently) {
                        printf ("Time to find quad %d-nearest: %.2f\n",
                                quadtry, CCutil_zeit () - kzeit);
                        fflush (stdout);
                    }
                }
            }
            if (!cycfname && !edgecycfname) {
                kzeit = CCutil_zeit ();
                if (tour_type == LK_GREEDY) {
                    if (CCkdtree_greedy_tour (&localkt, ncount,
                              &dat, incycle, &val, run_silently, &rstate)) {
                        fprintf (stderr, "CCkdtree greedy-tour failed\n");
                        rval = 1;
                        goto CLEANUP;
                    }
                } else if (tour_type == LK_QBORUVKA) {
                    if (CCkdtree_qboruvka_tour (&localkt, ncount,
                              &dat, incycle, &val, &rstate)) {
                        fprintf (stderr, "CCkdtree qboruvka-tour failed\n");
                        rval = 1;
                        goto CLEANUP;
                    }
                } else if (tour_type == LK_BORUVKA) {
                    if (CCkdtree_boruvka_tour (&localkt, ncount,
                              &dat, incycle, &val, &rstate)) {
                        fprintf (stderr, "CCkdtree boruvka-tour failed\n");
                        rval = 1;
                        goto CLEANUP;
                    }
                } else if (tour_type == LK_RANDOM) {
                    randcycle (ncount, incycle, &rstate);
                } else {
                    if (CCkdtree_nearest_neighbor_tour (&localkt, ncount,
                               CCutil_lprand (&rstate) % ncount, &dat,
                               incycle, &val, &rstate)) {
                        fprintf (stderr, "CCkdtree NN-tour failed\n");
                        rval = 1;
                        goto CLEANUP;
                    }
                }
                if (!run_silently) {
                    printf ("Time to grow tour: %.2f\n",
                            CCutil_zeit () - kzeit);
                    fflush (stdout);
                }
            }
            CCkdtree_free (&localkt);
        }
    } else if ((norm & CC_NORM_BITS) == CC_X_NORM_TYPE) {
        double xzeit = CCutil_zeit ();
        if (!goodfname && !edgegenfname) {
            if (nearnum) {
                if (CCedgegen_x_k_nearest (ncount, nearnum, &dat,
                        (double *) NULL, 1, &tempcount, &templist,
                        run_silently)) {
                    fprintf (stderr, "CCedgegen_x_k_nearest failed\n");
                    rval = 1;
                    goto CLEANUP;
                }
                if (!run_silently) {
                    printf ("Time to find %d-nearest: %.2f\n", nearnum,
                                                 CCutil_zeit () - xzeit);
                    fflush (stdout);
                }
            } else {
                if (CCedgegen_x_quadrant_k_nearest (ncount, quadtry, &dat,
                                 (double *) NULL, 1, &tempcount, &templist,
                                 run_silently)) {
                    fprintf (stderr, "x-quad nearest code failed\n");
                    rval = 1;
                    goto CLEANUP;
                }
                if (!run_silently) {
                    printf ("Time to find quad %d-nearest: %.2f\n", quadtry,
                                                 CCutil_zeit () - xzeit);
                    fflush (stdout);
                }
            }
        }
        if (!cycfname && !edgecycfname) {
            xzeit = CCutil_zeit ();
            if (tour_type == LK_GREEDY) {
                if (CCedgegen_x_greedy_tour (ncount, &dat, incycle, &val,
                        tempcount, templist, run_silently)) {
                    fprintf (stderr, "CCedgegen_x_greedy_tour failed\n");
                    rval = 1; goto CLEANUP;
                }
            } else if (tour_type == LK_QBORUVKA) {
                if (CCedgegen_x_qboruvka_tour (ncount, &dat, incycle, &val,
                        tempcount, templist, run_silently)) {
                    fprintf (stderr, "CCedgegen_x_qboruvka_tour failed\n");
                    rval = 1; goto CLEANUP;
                }
            } else if (tour_type == LK_RANDOM) {
                randcycle (ncount, incycle, &rstate);
            } else {
                if (CCedgegen_x_nearest_neighbor_tour (ncount,
                      CCutil_lprand (&rstate) % ncount, &dat, incycle, &val)) {
                    fprintf (stderr, "CCedgegen_x_nearest_neighbor_tour failed\n");
                    rval = 1;
                    goto CLEANUP;
                }
            }
            if (!run_silently) {
                printf ("Time to grow tour: %.2f\n", CCutil_zeit () - xzeit);
                fflush (stdout);
            }
        }
    } else {
        double jzeit = CCutil_zeit ();
        if (!goodfname && !edgegenfname) {
            if (!nearnum)
                nearnum = 4 * quadtry;
            if (CCedgegen_junk_k_nearest (ncount, nearnum, &dat,
                    (double *) NULL, 1, &tempcount, &templist, run_silently)) {
                fprintf (stderr, "CCedgegen_junk_k_nearest failed\n");
                rval = 1;
                goto CLEANUP;
            }
            if (!run_silently) {
                printf ("Time to find %d nearest: %.2f\n",
                         nearnum, CCutil_zeit () - jzeit);
                fflush (stdout);
            }
        }
        if (!cycfname && !edgecycfname) {
            jzeit = CCutil_zeit();
            if (tour_type == LK_GREEDY) {
                if (CCedgegen_junk_greedy_tour (ncount, &dat, incycle, &val,
                        tempcount, templist, run_silently)) {
                    fprintf (stderr, "CCedgegen_junk_greedy_tour failed\n");
                    rval = 1; goto CLEANUP;
                }
            } else if (tour_type == LK_QBORUVKA) {
                if (CCedgegen_junk_qboruvka_tour (ncount, &dat, incycle, &val,
                        tempcount, templist, run_silently)) {
                    fprintf (stderr, "CCedgegen_junk_qboruvka_tour failed\n");
                    rval = 1; goto CLEANUP;
                }
            } else if (tour_type == LK_RANDOM) {
                randcycle (ncount, incycle, &rstate);
            } else {
                if (CCedgegen_junk_nearest_neighbor_tour (ncount,
                       CCutil_lprand (&rstate) % ncount, &dat, incycle,
                       &val, run_silently)) {
                    fprintf (stderr, "CCedgegen_junk_nearest_neighbor_tour failed\n");
                    rval = 1;
                    goto CLEANUP;
                }
            }
            if (!run_silently) {
                printf ("Time to grow tour: %.2f\n", CCutil_zeit () - jzeit);
                fflush (stdout);
            }
        }
    }

    outcycle = CC_SAFE_MALLOC (ncount, int);
    if (!outcycle) {
        rval = 1;
        goto CLEANUP;
    }

    if (number_runs) {
        k = 0;
        best = BIGDOUBLE;
        do {
            printf ("\nStarting Run %d\n", k);
            if (CClinkern_tour (ncount, &dat, tempcount, templist, 100000000,
                   in_repeater, incycle, outcycle, &val, run_silently,
                   time_bound, length_bound, (char *) NULL, kick_type,
                   &rstate)) {
                fprintf (stderr, "CClinkern_tour failed\n");
                rval = 1;
                goto CLEANUP;
            }
            if (val < best) {
                best = val;
                if (saveit_final) {
                    if (CCutil_writecycle_edgelist (ncount, saveit_final, 
                            outcycle, &dat, binary_edges)) {
                        fprintf (stderr, "could not write the cycle\n");
                        rval = 1;
                        goto CLEANUP;
                    }
                }
            }
        } while (++k < number_runs);
        printf ("Overall Best Cycle: %.0f\n", val);
        fflush (stdout);
    } else {
        double lkzeit = CCutil_zeit ();
        int attempt = 1;
        do {
            if (CClinkern_tour (ncount, &dat, tempcount, templist, 10000000,
                   in_repeater, incycle, outcycle, &val, run_silently,
                   time_bound, length_bound, saveit_name, kick_type,
                   &rstate)) {
                fprintf (stderr, "CClinkern_tour failed\n");
                rval = 1;
                goto CLEANUP;
            }
            if (length_bound != -1 && val > length_bound) {
                printf ("Cycle of value %.0f  -  did not reach %.0f\n",
                    val, length_bound);
                printf ("Try again. Number of attempts: %d\n", ++attempt);
            }
        } while (length_bound != -1 && val > length_bound);
        if (saveit_final) {
            if (CCutil_writecycle_edgelist (ncount, saveit_final,
                        outcycle, &dat, binary_edges)) {
                fprintf (stderr, "could not write the cycle\n");
                rval = 1;
                goto CLEANUP;
            }
        }
        if (run_silently)
            printf ("Lin-Kernighan Running Time: %.2f\n",
                    CCutil_zeit () - lkzeit);
        printf ("Final Cycle: %.0f\n", val);
        fflush (stdout);
    }
    printf ("Total Running Time: %.2f\n", CCutil_zeit () - startzeit);
    fflush (stdout);

CLEANUP:

#ifndef BIG_PROBLEM
    CC_IFFREE (templist, int);
#endif
    CC_IFFREE (incycle, int);
    CC_IFFREE (outcycle, int);
    CCutil_freedatagroup (&dat);

    return rval;
}

static void randcycle (int ncount, int *cyc, CCrandstate *rstate)
{
    int i, k, temp;

    for (i = 0; i < ncount; i++)
        cyc[i] = i;

    for (i = ncount; i > 1; i--) {
        k = CCutil_lprand (rstate) % i;
        CC_SWAP (cyc[i - 1], cyc[k], temp);
    }
}


static int parseargs (int ac, char **av)
{
    int c, k;
    int inorm;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "a:bBD:E:g:G:h:k:lI:K:N:o:q:Qr:R:s:S:t:y:Y:", &boptind, &boptarg)) != EOF)
        switch (c) {
        case 'a':
            nearnum = atoi (boptarg);
            break;
        case 'B':
            binary_in = 1;
            break;
        case 'b':
            binary_in = 2;
            break;
        case 'E':
            binary_edges = 1;
            break;
        case 'D':
            edgegenfname = boptarg;
            break;
        case 'g':
            goodfname = boptarg;
            break;
        case 'G':
            gridsize = atoi(boptarg);
            break;
        case 'h':
            length_bound = atof (boptarg);
            break;
        case 'k':
            nnodes_want = atoi (boptarg);
            tsplib_in = 0;
            break;
        case 'I':
            k = atoi (boptarg);
            if (k == LK_RANDOM)        tour_type = LK_RANDOM;
            else if (k == LK_NEIGHBOR) tour_type = LK_NEIGHBOR;
            else if (k == LK_GREEDY)   tour_type = LK_GREEDY;
            else if (k == LK_BORUVKA)  tour_type = LK_BORUVKA;
            else if (k == LK_QBORUVKA) tour_type = LK_QBORUVKA;
            else fprintf (stderr, "unknown tour type, using default\n");
            break;
        case 'K':
            k = atoi (boptarg);
            if (k == CC_LK_RANDOM_KICK)         kick_type = CC_LK_RANDOM_KICK;
            else if (k == CC_LK_GEOMETRIC_KICK) kick_type = CC_LK_GEOMETRIC_KICK;
            else if (k == CC_LK_CLOSE_KICK)     kick_type = CC_LK_CLOSE_KICK;
            else if (k == CC_LK_WALK_KICK)      kick_type = CC_LK_WALK_KICK;
            else fprintf (stderr, "unknown kick type, using default\n");
            break;
        case 'N':
            inorm = atoi(boptarg);
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
                usage (av[0]);
                return 1;
            }
            tsplib_in = 0;
            break;
        case 'o':
            saveit_final = boptarg;
            break;
        case 'q':
            quadtry = atoi (boptarg);
            break;
        case 'Q':
            run_silently++;
            break;
        case 'r':
            number_runs = atoi (boptarg);
            break;
        case 'R':
            in_repeater = atoi (boptarg);
            break;
        case 's':
            seed = atoi (boptarg);
            break;
        case 'S':
            saveit_name = boptarg;
            break;
        case 't':
            time_bound = atof (boptarg);
            break;
        case 'y':
            cycfname = boptarg;
            break;
        case 'Y':
            edgecycfname = boptarg;
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
    fprintf (stderr, "   -K #  kick (%d-Random, %d-Geometric, %d-Close, %d-Random_Walk [default])\n",
           CC_LK_RANDOM_KICK, CC_LK_GEOMETRIC_KICK,
           CC_LK_CLOSE_KICK, CC_LK_WALK_KICK);
    fprintf (stderr, "   -o f  save final tour\n");
    fprintf (stderr, "   -S f  save tour in f after every 10000 kicks\n");
    fprintf (stderr, "   -D f  edgegen description file for the sparse edge set\n");
    fprintf (stderr, "   -q #  use quad #-nearest as the sparse set (default is 3)\n");
    fprintf (stderr, "   -a #  use #-nearest as the sparse edge set\n");
    fprintf (stderr, "   -g f  use the edges in file f as the sparse edge set\n");
    fprintf (stderr, "   -r #  number of runs\n");
    fprintf (stderr, "   -R #  number of kicks in iterated Lin-Kernighan (default is #nodes)\n");
    fprintf (stderr, "   -I #  generate starting cycle\n");
    fprintf (stderr, "           (%d-Rand, %d-NNeigh, %d-Greedy, %d-Boruvka, %d-QBoruvka[default])\n",
               LK_RANDOM, LK_NEIGHBOR, LK_GREEDY, LK_BORUVKA, LK_QBORUVKA);
    fprintf (stderr, "   -y f  starting cycle\n");
    fprintf (stderr, "   -Y f  starting cycle (as an edgelist)\n");
    fprintf (stderr, "   -t d  running time bound in seconds\n");
    fprintf (stderr, "   -h d  tour length bound (stop when we hit d)\n");
    fprintf (stderr, "   -Q    run silently\n");
    fprintf (stderr, "   -b    dat file in binary doubles\n");
    fprintf (stderr, "   -B    dat file in binary ints\n");
    fprintf (stderr, "   -E    edge files in binary\n");
    fprintf (stderr, "   -N #  norm (must specify if dat file is not a TSPLIB file)\n");
    fprintf (stderr, "         0=MAX, 1=L1, 2=L2, 3=3D, 4=USER, 5=ATT, 6=GEO, 7=MATRIX,\n");
    fprintf (stderr, "         8=DSJRAND, 9=CRYSTAL, 10=SPARSE, 11-15=RH-norm 1-5, 16=TOROIDAL\n");
    fprintf (stderr, "         17=GEOM, 18=JOHNSON\n");
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

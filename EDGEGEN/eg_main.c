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
/*               CODE FOR TESTING EDGE GENERATION ROUTINES                  */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 25, 1995                                                 */
/*                                                                          */
/*  For a short describtion see usage ()                                    */
/*                                                                          */
/****************************************************************************/


#include "machdefs.h"
#include "util.h"
#include "edgegen.h"

#define DEFAULT_EDGE_RANGE (5)   /* for random edge sets */

static int norm = CC_EUCLIDEAN;
static int seed = 0;
static int nnodes_want = 0;
static int gridsize = 0;
static int nearnum = 0;
static int quadnearnum = 0;
static int f2match_nearnum = 0;
static int usenodeweights = 0;
static int random_weight_limit = 0;
static int random_edge_count = 0;
static int random_edge_range = -1;
static int random_tour_count = 0;
static int nearest_tour_count = 0;
static int linkern_tour_count = 0;
static int linkern_kicks = 100;
static int use_greedy_in_linkern = 0;
static int use_boruvka_in_linkern = 0;
static int use_qboruvka_in_linkern = 0;
static int use_random_in_linkern = 0;
static int linkern_nearnum = 0;
static int linkern_quadnum = 0;
static int twoopt_tour_count = 0;
static int twoopt5_tour_count = 0;
static int threeopt_tour_count = 0;
static int nearest_twomatch_count = 0;
static int find_greedy_tour = 0;
static int find_boruvka_tour = 0;
static int find_qboruvka_tour = 0;
static int find_fractional_2match = 0;
static int find_spanning_tree = 0;
static int find_delaunay_edges = 0;
static int find_mlinkern_edges = 0;
static int binary_in = 0;
static int tsplib_in = 1;
static int binary_out = 0;

static char *nodefile = (char *) NULL;
static char *weightfile = (char *) NULL;
static char *outfile = (char *) NULL;
static char *pointfile = (char *) NULL;
static char *describefile = (char *) NULL;



int
    main (int ac, char **av);
static void
    usage (char *f);
static int
    parseargs (int ac, char **av);


int main (int ac, char **av)
{
    double szeit;
    double *wcoord = (double *) NULL;
    int ncount, maxlen;
    int rval = 0;
    int ecount = 0;
    int *elist = (int *) NULL;
    int *elen  = (int *) NULL;
    CCdatagroup dat;
    CCedgegengroup plan;
    CCrandstate rstate;

    CCutil_init_datagroup (&dat);

    seed = (int) CCutil_real_zeit ();
    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    CCutil_sprand (seed, &rstate);

    if ((!nnodes_want && !nodefile)) {
        usage (av[0]);
        rval = 1; goto CLEANUP;
    }
    if (nnodes_want) tsplib_in = 0;

    if (!nodefile && norm == CC_SPARSE) {
        if (random_edge_count == 0) {
            fprintf (stderr, "Must specify the number of edges\n");
            rval = 1; goto CLEANUP;
        }
        if (nearnum != 0 || quadnearnum != 0 || random_tour_count != 0 ||
            nearest_tour_count != 0 || linkern_tour_count != 0 ||
            twoopt_tour_count != 0 || twoopt5_tour_count != 0 ||
            threeopt_tour_count != 0 || nearest_twomatch_count != 0 ||
            find_greedy_tour != 0 || find_boruvka_tour != 0 ||
            find_qboruvka_tour != 0 || find_fractional_2match != 0 ||
            find_spanning_tree != 0 || find_delaunay_edges != 0 ||
            find_mlinkern_edges != 0) {
            fprintf (stderr, "Only permitted operation with SPARSE norm is a random edge set\n");
            rval = 1; goto CLEANUP;
        }
        ncount = nnodes_want;
        ecount = random_edge_count;
        if (random_edge_range == -1) {
            maxlen = DEFAULT_EDGE_RANGE * ncount;
        } else {
            maxlen = random_edge_range;
        }
        rval = CCutil_genedgelist (ncount, ecount, &elist, &elen,
                     (CCdatagroup *) NULL, maxlen, &rstate);
        if (rval) {
            fprintf (stderr, "CCutil_genedgelist failed\n"); goto CLEANUP;
        }
        if (outfile) {
            rval = CCutil_writeedges_int (ncount, outfile, ecount, elist,
                                          elen, binary_out);
            if (rval) {
                fprintf (stderr, "CCutil_writeedges_int failed\n");
                goto CLEANUP;
            }
        }
        goto CLEANUP;
    }

    if (tsplib_in) {
        if (CCutil_gettsplib (nodefile, &ncount, &dat)) {
            fprintf (stderr, "could not read the TSPLIB file\n");
            rval = 1;
            goto CLEANUP;
        }
        CCutil_dat_getnorm (&dat, &norm);
    } else {
        int allow_dups;
        int use_gridsize;
        
        ncount = nnodes_want;
        if (gridsize < 0) {
            use_gridsize = -gridsize;
            allow_dups = 0;
        } else if (gridsize > 0) {
            use_gridsize = gridsize;
            allow_dups = 1;
        } else {
            use_gridsize = ncount;
            allow_dups = 0;
        }
        if (CCutil_getdata (nodefile, binary_in, norm, &ncount, &dat,
                            use_gridsize, allow_dups, &rstate)) {
            rval = 1;
            goto CLEANUP;
        }
    }

    if (usenodeweights) {
        if (CCutil_getnodeweights (weightfile, ncount, random_weight_limit,
                                   &wcoord, &rstate)) {
            fprintf (stderr, "could not read the nodeweight file\n");
            rval = 1;
            goto CLEANUP;
        }
    }

    if (describefile) {
        if (CCedgegen_read (describefile, &plan)) {
            fprintf (stderr, "CCedgegen_read failed\n");
            rval = 1;
            goto CLEANUP;
        }
    } else {
        CCedgegen_init_edgegengroup (&plan);

        plan.random = random_edge_count;
        plan.nearest = nearnum;
        plan.quadnearest = quadnearnum;
        plan.delaunay = find_delaunay_edges;
        plan.mlinkern = find_mlinkern_edges;
        plan.tour.random_count = random_tour_count;
        plan.tour.nearest_count = nearest_tour_count;
        plan.tour.greedy = find_greedy_tour;
        plan.tour.boruvka = find_boruvka_tour;
        plan.tour.qboruvka = find_qboruvka_tour;
        plan.want_tree = find_spanning_tree;
        plan.tour.twoopt_count = twoopt_tour_count;
        plan.tour.twoopt5_count = twoopt5_tour_count;
        plan.tour.threeopt_count = threeopt_tour_count;
        if (linkern_tour_count) {
            plan.linkern.count = linkern_tour_count;
            if (!linkern_quadnum && !linkern_nearnum) {
                linkern_quadnum = 3;
            } else {
                plan.linkern.quadnearest = linkern_quadnum;
                plan.linkern.nearest = linkern_nearnum;
            }
            if (!use_greedy_in_linkern && !use_random_in_linkern &&
                !use_boruvka_in_linkern && !use_qboruvka_in_linkern) {
                plan.linkern.nearest_start = 1;
            } else {
                plan.linkern.greedy_start = use_greedy_in_linkern;
                plan.linkern.random_start = use_random_in_linkern;
                plan.linkern.boruvka_start = use_boruvka_in_linkern;
                plan.linkern.qboruvka_start = use_qboruvka_in_linkern;
            }
            plan.linkern.nkicks = linkern_kicks;
        }
        plan.nearest_twomatch_count = nearest_twomatch_count;
        plan.f2match.wantit = find_fractional_2match;
        if (f2match_nearnum) {
            plan.f2match_nearest.number = f2match_nearnum;
            plan.f2match_nearest.basic = 0;
            plan.f2match_nearest.priced = 1;
        }
    }

    szeit = CCutil_zeit ();

    if (CCedgegen_edges (&plan, ncount, &dat, wcoord, &ecount, &elist,
                         0, &rstate)) {
        fprintf (stderr, "CCedgegen_edges failed\n");
        rval = 1;
        goto CLEANUP;
    }

    printf ("Edgegen running time: %.2f (seconds)\n", CCutil_zeit () - szeit);
    fflush (stdout);

    if (outfile && ecount) {
        if (CCutil_writeedges (ncount, outfile, ecount, elist, &dat,
                               binary_out)) {
            fprintf (stderr, "Could not write the edge set\n");
            rval = 1;
            goto CLEANUP;
        }
    }

    if (pointfile && ncount) {
        if (CCutil_writedata (pointfile, binary_out, ncount, &dat)) {
            fprintf (stderr, "Could not write the point set\n");
            rval = 1;
            goto CLEANUP;
        }
    }

CLEANUP:

    CCutil_freedatagroup (&dat);
    CC_IFFREE (wcoord, double);
    CC_IFFREE (elist, int);
    CC_IFFREE (elist, int);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c, inorm;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "A:bB:C:de:D:f:FGhHiIk:K:L:s:m:M:n:N:o:Op:q:r:R:ST:uU:vw:W:x:y:?", &boptind, &boptarg))
            != EOF)
        switch (c) {
        case 'A':
            twoopt_tour_count = atoi (boptarg);
            break;
        case 'b':
            binary_in = 1;
            break;
        case 'B':
            twoopt5_tour_count = atoi (boptarg);
            break;
        case 'C':
            threeopt_tour_count = atoi (boptarg);
            break;
        case 'd':
            find_delaunay_edges = 1;
            break;
        case 'e':
            random_edge_count = atoi (boptarg);
            break;
        case 'D':
            describefile = boptarg;
            break;
        case 'f':
            f2match_nearnum = atoi (boptarg);
            break;
        case 'F':
            find_fractional_2match = 1;
            break;
        case 'G':
            find_greedy_tour = 1;
            break;
        case 'h':
            find_boruvka_tour = 1;
            break;
        case 'H':
            use_boruvka_in_linkern = 1;
            break;
        case 'i':
            find_qboruvka_tour = 1;
            break;
        case 'I':
            use_qboruvka_in_linkern = 1;
            break;
        case 'k':
            nnodes_want = atoi (boptarg);
            break;
        case 'K':
            random_edge_range = atoi (boptarg);
            break;
        case 'L':
            linkern_tour_count = atoi (boptarg);
            break;
        case 'm':
            find_mlinkern_edges = atoi (boptarg);
            break;
        case 'M':
            nearest_twomatch_count = atoi (boptarg);
            break;
        case 'n':
            nearnum = atoi (boptarg);
            break;
        case 'p':
            pointfile = boptarg;
            break;
        case 'T':
            nearest_tour_count = atoi (boptarg);
            break;
        case 'o':
            outfile = boptarg;
            break;
        case 'O':
            binary_out = 1;
            break;
        case 'q':
            quadnearnum = atoi (boptarg);
            break;
        case 'r':
            gridsize = atoi(boptarg);
            break;
        case 'R':
            linkern_kicks = atoi (boptarg);
            break;
        case 's':
            seed = atoi (boptarg);
            break;
        case 'S':
            find_spanning_tree = 1;
            break;
        case 'u':
            use_greedy_in_linkern = 1;
            use_random_in_linkern = 0;
            break;
        case 'U':
            random_tour_count = atoi (boptarg);
            break;
        case 'v':
            use_random_in_linkern = 1;
            use_greedy_in_linkern = 0;
            break;
        case 'w':
            usenodeweights = 1;
            weightfile = boptarg;
            break;
        case 'W':
            usenodeweights = 1;
            random_weight_limit = atoi (boptarg);
            break;
        case 'x':
            linkern_quadnum = atoi (boptarg);
            break;
        case 'y':
            linkern_nearnum = atoi (boptarg);
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
    fprintf (stderr, "Usage: %s [- see below -] [dat file]\n", f);
    fprintf (stderr, "   -b    dat file in binary-ints\n");
    fprintf (stderr, "   -w f  node weight file\n");
    fprintf (stderr, "   -W #  use random node weights, from 0 to # - 1\n");
    fprintf (stderr, "   -k #  number of nodes for random problem\n");
    fprintf (stderr, "   -K #  in a SPARSE-norm problem, use edge weights from 0 to # - 1\n");
    fprintf (stderr, "   -s #  random seed\n");
    fprintf (stderr, "   -D f  description file\n");
    fprintf (stderr, "   -e #  find # random edges\n");
    fprintf (stderr, "   -n #  find # nearest graph\n");
    fprintf (stderr, "   -q #  find quadrant # nearest graph\n");
    fprintf (stderr, "   -d    find Delaunay triangulation\n");
    fprintf (stderr, "   -f #  find # nearest using f2match reduced costs\n");
    fprintf (stderr, "   -U #  find # random tours\n");
    fprintf (stderr, "   -T #  find # nearest neighbor tours\n");
    fprintf (stderr, "   -G    find greedy tour\n");
    fprintf (stderr, "   -h    find boruvka tour\n");
    fprintf (stderr, "   -i    find quick boruvka tour\n");
    fprintf (stderr, "   -(A B C) #  find # (2opt, 2.5opt, 3opt) tours\n");
    fprintf (stderr, "   -L #  find # linkern tours\n");
    fprintf (stderr, "   -m #  find # linkern matchings\n");
    fprintf (stderr, "   -r n  use nXn grid for random points, no dups if n<0\n");
    fprintf (stderr, "   -R #  use # kicks in linkern (default: 100)\n");
    fprintf (stderr, "   -u    use greedy starting tour for linkern\n");
    fprintf (stderr, "   -H    use boruvka starting tour for linkern\n");
    fprintf (stderr, "   -I    use quick boruvka starting tour for linkern\n");
    fprintf (stderr, "   -v    use random starting tours for linkern\n");
    fprintf (stderr, "   -x #  use # quadnearest in linkern\n");
    fprintf (stderr, "   -y #  use # nearest in linkern (can use x & y)\n");
    fprintf (stderr, "   -M #  find # nearest neighbor 2-matchings\n");
    fprintf (stderr, "   -S    find min spanning tree\n");
    fprintf (stderr, "   -F    find fractional twomatch (not priced)\n");
    fprintf (stderr, "   -o f  write the cycle or edge set to f\n");
    fprintf (stderr, "   -p f  write the point set to f\n");
    fprintf (stderr, "   -O    use binary output\n");
    fprintf (stderr, "   -N #  norm (must specify if dat file is not a TSPLIB file)\n");
    fprintf (stderr, "         0=MAX, 1=L1, 2=L2, 3=3D, 4=USER, 5=ATT, 6=GEO, 7=MATRIX,\n");
    fprintf (stderr, "         8=DSJRAND, 9=CRYSTAL, 10=SPARSE, 11-15=RH-norm 1-5, 16=TOROIDAL\n");
    fprintf (stderr, "         17=GEOM, 18=JOHNSON\n");
}

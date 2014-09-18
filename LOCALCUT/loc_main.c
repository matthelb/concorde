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

#include "machdefs.h"
#include "util.h"
#include "cut.h"
#include "tsp.h"
#include "localcut.h"

static char *edgefile = (char *) NULL;
static int seed = 0;
static int spheres = 1;
static int binary_in = 0;
static int noshrink = 0;
static int uncivilized = 0;
static int maxchunksize = 15;
static int maxspheresize = 0;
static int dummy = 0;
static int permute = 0;
static int weighted = 0;
static int dumpcuts = 0;
static int nolift = 1;
static double eps = 0.0;
static int run_silently = 0;

static int
    parseargs (int ac, char **av);

static void
    usage (char *f);



int main (int ac, char **av)
{
    int ncount;
    int ecount;
    int *elist = (int *) NULL;
    double *x = (double *) NULL;
    double szeit;
    int rval = 0;
    CCchunk_flag flags;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    int cutcount = 0;
    CCtsp_lpcut_in *c = (CCtsp_lpcut_in *) NULL;
    CCtsp_lpcut_in *cnext = (CCtsp_lpcut_in *) NULL;
    CCtsp_lpgraph g;
    double cutval;
    double maxval;
    CCchunk_localcut_timer timer;
    int i;
    int cuthist[14];
    static double cutthresh[14] = {0.0, 0.001, 0.002, 0.005, 0.01, 0.02,
                                   0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 1000.0};
    CCrandstate rstate;

    seed = (int) CCutil_real_zeit ();
    CCtsp_init_lpgraph_struct (&g);

    rval = parseargs (ac, av);
    if (rval) return rval;

    CCutil_sprand (seed, &rstate);
    printf ("Using random seed %d\n", seed); fflush (stdout);

    CCchunk_init_localcut_timer (&timer);

    rval = CCutil_getedges_double (&ncount, edgefile, &ecount, &elist,
                                   &x, binary_in);
    if (rval) {
        fprintf (stderr, "getedges failed\n");
        goto CLEANUP;
    }

    szeit = CCutil_zeit ();

    if (maxspheresize == 0) maxspheresize = maxchunksize - 2;

    flags.dummy        = dummy;
    flags.permute      = permute;
    flags.weighted     = weighted;
    flags.spheres      = spheres;
    flags.uncivilized  = uncivilized;
    flags.noshrink     = noshrink;
    flags.maxchunksize = maxchunksize;
    flags.spheresize   = maxspheresize;
    flags.nolift       = nolift;

    rval = CCchunk_localcuts (&cuts, &cutcount, ncount, ecount, elist, x,
                              eps, flags, &timer, run_silently, &rstate);
    if (rval) {
        fprintf (stderr, "CCchunk_localcuts failed\n");
        goto CLEANUP;
    }


    szeit = CCutil_zeit () - szeit;
    
    CCchunk_print_localcut_timer (&timer);

    rval = CCtsp_build_lpgraph (&g, ncount, ecount, elist, (int *) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_build_lpgraph failed\n");
        goto CLEANUP;
    }

    rval = CCtsp_build_lpadj (&g, 0, g.ecount);
    if (rval) {
        fprintf (stderr, "CCtsp_build_lpadj failed\n");
        goto CLEANUP;
    }

    for (i=0; i<14; i++) cuthist[i] = 0;
    maxval = 0.0;
    for (c = cuts; c; c = c->next) {
        cutval = -CCtsp_cutprice (&g, c, x);
        for (i=0; i<13; i++) {
            if (cutval <= cutthresh[i]) {
                cuthist[i]++;
                break;
            }
        }
        if (i == 13) cuthist[i]++;
        
        if (cutval > maxval) maxval = cutval;
        if (dumpcuts) {
            printf ("Cut violation %f:\n", cutval);
            CCtsp_print_lpcut_in (c);
        }
    }
    if (dumpcuts) {
        printf ("\n");
    }

    printf ("cut slack distribution:\n");
    if (cuthist[0]) {
        printf ("        s <= %.3f: %4d\n", cutthresh[0], cuthist[0]);
    }
    for (i=1; i<13; i++) {
        if (cuthist[i]) {
            printf ("%.3f < s <= %.3f: %4d\n", cutthresh[i-1], cutthresh[i],
                    cuthist[i]);
        }
    }
    if (cuthist[13]) {
        printf ("%.3f < s         : %4d\n", cutthresh[13], cuthist[13]);
    }
    printf ("%d cuts (max viol %f) found in %.2f seconds\n",
            cutcount, maxval, szeit);
    fflush (stdout);

    rval = 0;

CLEANUP:
    for (c = cuts; c; c = cnext) {
        cnext = c->next;
        CCtsp_free_lpcut_in (c);
    }

    CCtsp_free_lpgraph (&g);
    CC_IFFREE(elist, int);
    CC_IFFREE(x, double);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt(ac, av, "bcC:dDeE:fFLpr:sS:w", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'b':
            binary_in = 1;
            break;
        case 'c':
            uncivilized = 1;
            break;
        case 'C':
            maxchunksize = atoi(boptarg);
            break;
        case 'd':
            dummy = 1;
            break;
        case 'D':
            dumpcuts = 1;
            break;
        case 'e':
            spheres = 0;
            break;
        case 'E':
            eps = atof(boptarg);
            break;
        case 'f':
            nolift = 2;
            break;
        case 'F':
            nolift = 3;
            break;
        case 'L':
            nolift = 0;
            break;
        case 'p':
            permute = 1;
            break;
        case 'r':
            seed = atoi(boptarg);
            break;
        case 's':
            noshrink = 1;
            break;
        case 'S':
            maxspheresize = atoi(boptarg);
            break;
        case 'w':
            weighted = 1;
            break;
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind >= ac) {
        usage (av[0]);
        return 1;
    }
    edgefile = av[boptind++];

    if (boptind != ac) {
        usage (av[0]);
        return 1;
    }
    
    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s [-see below-] edge_file\n", f);
    fprintf (stderr, "   -b    edge_file in binary format\n");
    fprintf (stderr, "   -c    don't civilize chunks\n");
    fprintf (stderr, "   -C n  max chunk size n\n");
    fprintf (stderr, "   -d    use dummy spheres\n");
    fprintf (stderr, "   -D    output the cuts found\n");
    fprintf (stderr, "   -e    use equivalence classes instead of spheres\n");
    fprintf (stderr, "   -E f  use epsilon f (default 0.0)\n");
    fprintf (stderr, "   -f    dump the faulty chunks found\n");
    fprintf (stderr, "   -F    dump the faulty chunks and tight tours found\n");
    fprintf (stderr, "   -L    lift the faults found (not with -f or -F)\n");
    fprintf (stderr, "   -p    use permuted spheres\n");
    fprintf (stderr, "   -r n  use random seed n\n");
    fprintf (stderr, "   -s    do not shrink the vector\n");
    fprintf (stderr, "   -S n  max sphere size n\n");
    fprintf (stderr, "   -w    use weighted spheres\n");
}

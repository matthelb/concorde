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
/*        Demo of Min-Cuts for Directed and Undirected Graphs               */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal,  and Cook                       */
/*  Date: September, 1994 (Bonn)                                            */
/*        July 28, 1997 (bico)                                              */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "cut.h"

#undef  USE_DIRECTED_GRAPH

static int source = -1;
static int sink = -1;
static int showcut = 0;
static int violatedcuts = 0;
static int use_shrink = 1;
static int find_global = 0;
static int find_gomoryhu = 0;
static int binary_in = 0;
static char *fname = (char *) NULL;
static char *setfname = (char *) NULL;
static int seed = 0;

int
    main (int ac, char **av);

static void
    usage (char *f),
    display_cut (int *cut, int count);

static int
    parseargs (int ac, char **av),
    shrink_ones (int ncount, int ecount, int *elist, double *dlen,
            int *oncount, int *oecount, int **olist, double **olen,
            double *minval),
    display_all_cuts (double val, int cnt, int *cut, void *pass_param);


#ifdef USE_DIRECTED_GRAPH
static int
    duplicate_edges (int ncount, int ecount, int *elist, double *ecap,
                     int **tlist, double **tcap);
#endif

int main (int ac, char **av)
{
    int rval = 0;
    double val;
    double szeit;
    int i, cutcount, scount = 0;
    int ncount, ecount;
    int *elist = (int *) NULL;
    int *cut   = (int *) NULL;
    int **mycut = (int **) NULL;
    int *slist = (int *) NULL;
    double *ecap = (double *) NULL;
    CCrandstate rstate;
    FILE *setin = (FILE *) NULL;

    seed = (int) CCutil_real_zeit ();

    if (parseargs (ac, av)) goto CLEANUP;

    CCutil_sprand (seed, &rstate);
    
    if (find_global && violatedcuts) {
        fprintf (stderr, "use at most one of -p and -c arguments\n");
        goto CLEANUP;
    }

#ifdef USE_DIRECTED
    if (find_global || violatedcuts) {
        fprintf (stderr, "not set up for global cut in directed graphs\n");
        goto CLEANUP;
    }
#endif


    rval = CCutil_getedges_double (&ncount, fname, &ecount, &elist, &ecap,
                                   binary_in);
    if (rval) {
        fprintf (stderr, "CCutil_getedges_double failed\n"); goto CLEANUP;
    }

    if (showcut) {
        mycut = &cut;
    } else {
        mycut = (int **) NULL;
    }

    if (setfname) {
        setin = fopen (setfname, "r");
        if (!setin) {
            fprintf (stderr, "could not open %s for reading set\n", setfname);
            rval = 1;  goto CLEANUP;
        }
        fscanf (setin, "%d", &scount);
        slist = CC_SAFE_MALLOC (scount, int);
        CCcheck_NULL (slist, "out of memory for slist");

        for (i = 0; i < scount; i++) {
            fscanf (setin, "%d", &slist[i]);
        }
        printf ("Find min-cut containing nodes:");
        for (i = 0; i < scount; i++) printf (" %d", slist[i]);
        printf ("\n"); fflush (stdout);

        rval = CCcut_mincut_containing_set (ncount, ecount, elist, ecap,
                          scount, slist, &val, &cut, &cutcount, 0, &rstate);
        CCcheck_rval (rval, "CCcut_mincut_containing_set failed");

        if (showcut) display_cut (cut, cutcount);
        goto CLEANUP;
    } else if (find_gomoryhu) {
        CC_GHtree T;

        szeit = CCutil_zeit ();
        CCcut_GHtreeinit (&T);

        rval = CCcut_gomory_hu (&T, ncount, ecount, elist, ecap, 0,
                                (int *) NULL, &rstate);
        if (rval) {
            fprintf (stderr, "CCcut_gomory_hu failed\n"); goto CLEANUP;
        }
        printf ("Gomory-Hu Tree: %.2f seconds)\n", CCutil_zeit() - szeit);
        fflush (stdout);
        if (showcut) {
            CCcut_GHtreeprint (&T);
        }
        CCcut_GHtreefree (&T);
        goto CLEANUP;
    } else if (find_global) {
        szeit = CCutil_zeit ();
        rval = CCcut_mincut (ncount, ecount, elist, ecap, &val, mycut,
                            &cutcount);
        if (rval) {
            fprintf (stderr, "CCcut_mincut failed\n"); goto CLEANUP;
        }
        printf ("Minimum Cut Value: %f (%.2f seconds)\n", val,
                CCutil_zeit() - szeit);
        fflush (stdout);
        if (showcut) display_cut (cut, cutcount);
        goto CLEANUP;
    } else if (violatedcuts) {
        szeit = CCutil_zeit ();
        rval = CCcut_violated_cuts (ncount, ecount, elist, ecap,
                2.0 - CC_MINCUT_ONE_EPSILON, display_all_cuts, (void *) NULL);
        if (rval) {
            fprintf (stderr, "CCcut_violated_cuts failed\n");
            goto CLEANUP;
        }
        printf ("Running time: %.2f seconds\n", CCutil_zeit() - szeit);
        goto CLEANUP;
    }

    szeit = CCutil_zeit ();
    if (use_shrink) {
        int tncount, tecount;
        int *telist = (int *) NULL;
        double *tecap = (double *) NULL;
        double minval = CC_MINCUT_BIGDOUBLE;
        double sszeit = CCutil_zeit ();

        rval = shrink_ones (ncount, ecount, elist, ecap, &tncount, &tecount,
                            &telist, &tecap, &minval);
        if (rval) {
            fprintf (stderr, "shrink_ones failed\n"); goto CLEANUP;
        }
        printf ("Shrunk graph has %d nodes and %d edges\n", tncount, tecount);
        if (minval != CC_MINCUT_BIGDOUBLE) {
            printf ("Shrinking found cut of value %f\n", minval);
        }
        fflush (stdout);
        printf ("Time for shrinking: %.2f seconds\n", CCutil_zeit () - sszeit);
        fflush (stdout);
        CC_IFFREE(elist, int);
        CC_IFFREE(ecap, double);
        ncount = tncount;
        ecount = tecount;
        elist = telist;
        ecap = tecap;
    }

#ifdef USE_DIRECTED_GRAPH
    {
        int tncount, tecount;
        int *telist = (int *) NULL;
        double *tecap = (double *) NULL;

        rval = duplicate_edges (ncount, ecount, elist, ecap, &telist, &tecap);
        if (rval) {
            fprintf (stderr, "duplicated_edges failed\n"); goto CLEANUP;
        }
        CC_IFFREE(elist, int);
        CC_IFFREE(ecap, double);
        ecount *= 2;
        elist = telist;
        ecap = tecap;
    }
#endif


    if (source == -1)
        source = 0;

    if (sink != -1) {
        if (source < 0 || sink < 0 || source >= ncount ||
            sink >= ncount || source == sink) {
            printf ("Bad source sink pair\n");
            fflush (stdout);
            goto CLEANUP;
        }

        szeit = CCutil_zeit ();
        rval = CCcut_mincut_st (ncount, ecount, elist, ecap, source, sink,
                            &val, mycut, &cutcount);
        if (rval) {
           fprintf (stderr, "mincut_st failed\n"); goto CLEANUP;
        }
        printf ("Cut value: %f\n", val);
        printf ("Running Time: %.2f (seconds)\n", CCutil_zeit () - szeit);
        fflush (stdout);

        if (showcut) display_cut (cut, cutcount);
    } else {
        double minval = CC_MINCUT_BIGDOUBLE;
        double fzeit = CCutil_zeit ();

        printf ("compute all cuts from source node %d\n", source);
        fflush (stdout);
        for (i = 0; i < ncount; i++) {
            if (i != source) {
                rval = CCcut_mincut_st (ncount, ecount, elist, ecap, source, i,
                                        &val, (int **) NULL, (int *) NULL);
                if (rval) {
                    fprintf (stderr, "mincut_digraph failed\n"); goto CLEANUP;
                }
                if (val < minval)
                    minval = val;
                printf ("."); fflush (stdout);
                if (i % 75 == 0)
                    printf ("%d\n", i);
            }
        }
        printf ("\nMinimum Cut Value: %f\n", minval);
        printf ("Running Time: %.2f (seconds)\n", CCutil_zeit () - fzeit);
        fflush (stdout);
    }

    printf ("Total Running Time: %.2f (seconds)\n", CCutil_zeit () - szeit);
    fflush (stdout);

CLEANUP:

    CC_IFFREE (elist, int);
    CC_IFFREE (ecap, double);
    CC_IFFREE (cut, int);
    CC_IFFREE (slist, int);
    if (setin) fclose (setin);
    return 0;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "bcCGpr:s:St:T:", &boptind, &boptarg)) != EOF)
        switch (c) {
        case 'b':
            binary_in = 1;
            break;
        case 'c':
            violatedcuts = 1;
            break;
        case 'C':
            showcut++;
            break;
        case 'G':
            find_gomoryhu = 1;
            break;
        case 'p':
            find_global = 1;
            break;
        case 'r':
            seed = atoi (boptarg);
            break;
        case 's':
            source = atoi (boptarg);
            break;
        case 'S':
            use_shrink = 0;
            break;
        case 't':
            sink = atoi (boptarg);
            break;
        case 'T':
            setfname = boptarg;
            use_shrink = 0;
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
    fprintf (stderr, "    C:   display the min cut (or G-H tree)\n");
    fprintf (stderr, "    c:   display all cuts < 2.0\n");
    fprintf (stderr, "    G:   find Gomory-Hu tree (no shrinking)\n");
    fprintf (stderr, "    p:   use Padberg-Rinaldi style shrinking\n");
    fprintf (stderr, "    r #: random seed\n");
    fprintf (stderr, "    S:   do not use the TSP shrink routines\n");
    fprintf (stderr, "    s #: source\n");
    fprintf (stderr, "    t #: sink\n");
    fprintf (stderr, "    T f: find min cut containing node set in file f\n");
}

static void display_cut (int *cut, int count)
{
    int i;

    printf ("MIN CUT:\n");
    for (i = 0; i < count; i++) {
        printf ("%3d ", cut[i]);
        if (i % 10 == 9)
            printf ("\n");
    }
    if (i % 10) printf ("\n");
    fflush (stdout);
}

#ifdef USE_DIRECTED_GRAPH
static int duplicate_edges (int ncount, int ecount, int *elist, double *ecap,
                            int **tlist, double **tcap)
{
    int i;

    *tlist = (int *) NULL;
    *tcap = (double *) NULL;


    /* Convert the graph to a digraph */

    *tlist = CC_SAFE_MALLOC (4 * ecount, int);
    *tcap  = CC_SAFE_MALLOC (2 * ecount, double);
    if (!*tlist || !*tcap) {
        fprintf (stderr, "Out of memory in duplicate_edges\n");
        CC_IFFREE (*tlist, int);
        CC_IFFREE (*tcap, double);
        return 1;
    }

    for (i = 0; i < ecount; i++) {
        (*tlist)[4 * i] = (*tlist)[4 * i + 3] = elist[2 * i];
        (*tlist)[4 * i + 1] = (*tlist)[4 * i + 2] = elist[2 * i + 1];
        (*tcap)[2 * i] = (*tcap)[2 * i + 1] = ecap[i];
    }

    return 0;
}
#endif


static int shrink_ones (int ncount, int ecount, int *elist, double *dlen,
        int *oncount, int *oecount, int **olist, double **olen, double *minval)
{
    int rval = 0;
    CC_SRKgraph G;

    CCcut_SRK_init_graph (&G);
    rval = CCcut_SRK_buildgraph (&G, ncount, ecount, elist, dlen);
    if (rval) {
        fprintf (stderr, "buildgraph failed in shrink_ones\n"); goto CLEANUP;
    }
    CCcut_SRK_subtour_shrink (&G, minval, CC_MINCUT_ONE_EPSILON,
            (CC_SRKcallback *) NULL, (int **) NULL, (int *) NULL);

    rval = CCcut_SRK_grab_edges (&G, oncount, oecount, olist, olen,
                          (CC_SRKexpinfo *) NULL);
    if (rval) {
        fprintf (stderr, "grab edges failed in shrink_ones\n"); goto CLEANUP;
    }


CLEANUP:

    CCcut_SRK_free_graph (&G);
    return rval;
}

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

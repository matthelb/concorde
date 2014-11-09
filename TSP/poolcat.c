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
/*                         MERGE CUTPOOLS                                   */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: May 2, 1997                                                       */
/*                                                                          */
/*  SEE short decsribtion in usage ().                                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"

static int seed = 0;


static void
    usage (char *f);



int main (int ac, char **av)
{
    double szeit;
    int i, k, count;
    CCtsp_lpcut_in *c, *cnext;
    CCtsp_lpcut_in *cuts   = (CCtsp_lpcut_in *) NULL;
    CCtsp_lpcuts *pool     = (CCtsp_lpcuts *) NULL;
    CCtsp_lpcuts *nextpool = (CCtsp_lpcuts *) NULL;
    int ncount;
    int textout, textin;
    int *ptour = (int *) NULL;
    int *qtour = (int *) NULL;
    int *perm  = (int *) NULL;
    int rval = 0;
    CCrandstate rstate;

    if (ac < 2) {
        usage (*av);
        return 0;
    }

    szeit = CCutil_zeit ();
    seed = (int) CCutil_real_zeit ();
    CCutil_sprand (seed, &rstate);

    k = 1;
    if (av[k][0] == '-' && av[k][1] == 't') {
        CCdatagroup dat;
        printf ("Read master ...\n");
        if (CCutil_getmaster (av[k+1], &ncount, &dat, &ptour)) {
            fprintf (stderr, "CCutil_getmaster failed\n");
            return 1;
        }
        CCutil_freedatagroup (&dat);
        textout = 1;
        k += 2;
    } else {
        textout = 0;
    }

    if (av[k][0] == '-' && av[k][1] == 's') {
        CCdatagroup dat;
        printf ("Read master ...\n");
        if (CCutil_getmaster (av[k+1], &ncount, &dat, &qtour)) {
            fprintf (stderr, "CCutil_getmaster failed\n");
            return 1;
        }
        CCutil_freedatagroup (&dat);
        textin = 1;
        k += 2;
    } else {
        textin = 0;
    }

    if (av[k][0] == '-' && av[k][1] == 'p') {
        FILE *tin = fopen (av[k+1], "r");

        printf ("Read permutation tour ...\n");
        if (!tin) {
            fprintf (stderr, "could not open %s for reading\n", av[k+1]);
            rval = 1; goto CLEANUP;
        }

        if (fscanf (tin, "%d", &ncount) != 1) {
            fprintf (stderr, "perm file in wrong format\n");
            rval = 1; fclose (tin); goto CLEANUP;
        }
        perm = CC_SAFE_MALLOC (ncount, int);
        if (!perm) {
            fprintf (stderr, "out of memory in main\n");
            rval = 1; fclose (tin); goto CLEANUP;
        }
        for (i = 0; i < ncount; i++) {
            if (fscanf (tin, "%d", &(perm[i])) != 1) {
                fprintf (stderr, "perm file in wrong format\n");
                rval = 1; fclose (tin); goto CLEANUP;
            }
        }
        fclose (tin);
        k += 2;
    }

    if (textin) {
        printf ("Number of Nodes: %d\n", ncount); fflush (stdout);

        rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &pool);
        if (rval) {
            fprintf (stderr, "CCtsp_init_cutpool failed\n"); goto CLEANUP;
        }

        cuts = (CCtsp_lpcut_in *) NULL;
        rval = CCtsp_file_cuts (av[k], &cuts, &count, ncount, qtour);
        if (rval) {
            fprintf (stderr, "CCtsp_file_cuts failed\n"); goto CLEANUP;
        }
        printf ("File has %d cuts\n", count); fflush (stdout);
        for (c = cuts; c; c = cnext) {
            cnext = c->next;
            rval = CCtsp_add_to_cutpool_lpcut_in (pool, c);
            if (rval) {
                fprintf (stderr, "CCtsp_add_to_cutpool_lpcut_in failed\n");
                goto CLEANUP;
            }
            CCtsp_free_lpcut_in (c);
            CC_FREE (c, CCtsp_lpcut_in);
        }
        k++;
    } else {
        ncount = 0;
        rval = CCtsp_init_cutpool (&ncount, av[k], &pool);
        if (rval) {
            fprintf (stderr, "CCtsp_init_cutpool failed\n"); goto CLEANUP;
        }
        printf ("Initial Pool: %d nodes %d cuts\n", ncount, pool->cutcount);
        fflush (stdout);
        k++;
    }

    for (; k < ac; k++) {
        printf ("Adding Pool %s ... ", av[k]);
        fflush (stdout);

        if (textin) {
            cuts = (CCtsp_lpcut_in *) NULL;
            rval = CCtsp_file_cuts (av[k], &cuts, &count, ncount, qtour);
            if (rval) {
                fprintf (stderr, "CCtsp_file_cuts failed\n"); goto CLEANUP;
            }
            for (c = cuts; c; c = cnext) {
                cnext = c->next;
                rval = CCtsp_add_to_cutpool_lpcut_in (pool, c);
                if (rval) {
                    fprintf (stderr, "CCtsp_add_to_cutpool_lpcut_in failed\n");
                    goto CLEANUP;
                }
                CCtsp_free_lpcut_in (c);
                CC_FREE (c, CCtsp_lpcut_in);
            }
        } else {
            rval = CCtsp_init_cutpool (&ncount, av[k], &nextpool);
            if (rval) {
                fprintf (stderr, "CCtsp_init_cutpool failed\n"); goto CLEANUP;
            }

            for (i = 0; i < nextpool->cutcount; i++) {
                rval = CCtsp_add_to_cutpool (pool, nextpool,
                                             &(nextpool->cuts[i]));
                if (rval) {
                    fprintf (stderr, "CCtsp_add_to_cutpool failed\n");
                    goto CLEANUP;
                }
            }
            CCtsp_free_cutpool (&nextpool);
        }
        printf ("%d\n", pool->cutcount); fflush (stdout);
    }


    printf ("Final Pool: %d cuts\n", pool->cutcount);
    fflush (stdout);

    if (textout) {
        printf ("Write text file ...\n"); fflush (stdout);
        if (perm) {
            int *pperm = (int *) NULL;

            pperm = CC_SAFE_MALLOC (ncount, int);
            if (!pperm) {
                fprintf (stderr, "out of memory in main\n");
                rval = 1; goto CLEANUP;
            }
            for (i = 0; i < ncount; i++) {
                pperm[i] = perm[ptour[i]];
            }
            for (i = 0; i < ncount; i++) {
                ptour[i] = pperm[i];
            }
            CC_FREE (pperm, int);
        }
        rval = CCtsp_file_cuts_write ("merge.txt", pool, ptour);
        if (rval) {
            fprintf (stderr, "CCtsp_file_cuts_write failed\n");
            goto CLEANUP;
        }
    } else {
        rval = CCtsp_write_cutpool (ncount, "merge.pul", pool);
        if (rval) {
            fprintf (stderr, "CCtsp_write_cutpool failed\n");
            goto CLEANUP;
        }
    }

    printf ("Running Time: %.2f seconds\n", CCutil_zeit () - szeit);
    fflush (stdout);

CLEANUP:

    if (pool) {
        CCtsp_free_cutpool (&pool);
    }
    if (nextpool) {
        CCtsp_free_cutpool (&nextpool);
    }
    CC_IFFREE (ptour, int);
    CC_IFFREE (qtour, int);

    return rval;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s [-t] pool1 pool2 ... \n", f);
    fprintf (stderr, "       -t f:  to write a text file specify a master\n");
    fprintf (stderr, "       -p f:  to permute nodes when writing text\n");
    fprintf (stderr, "       -s f:  to read text files specify a master\n");
    fprintf (stderr, "Note: merged pool will be written merge.pul (.txt)\n");
    fprintf (stderr, "      the master files are used to map the nodes\n");
    fprintf (stderr, "      the permuation given by p can be used to\n");
    fprintf (stderr, "      handle the different node orders in dat and tsp\n");
}

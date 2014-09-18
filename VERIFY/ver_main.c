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
#include "tsp.h"
#include "util.h"
#include "verify.h"

static int textin = 0;
static int probin = 0;
static int poolin = 0;
static char *cutname = (char *) NULL;
static int in_ncount = 0;
static int run_silently = 0;


static int
    verify_text_pool (char *poolname, int nodecount),
    verify_binary_pool (char *poolname),
    verify_problem (char *probname, int silent),
    parseargs (int ac, char **av);

static void
    usage (char *fname);


int main (int ac, char **av)
{
    int rval;
    CCutil_timer z;

    CCutil_init_timer (&z, "Verification");
    CCutil_start_timer (&z);
    
    rval = parseargs (ac, av);
    if (rval) return -1;

    if (textin) {
        rval = verify_text_pool (cutname, in_ncount);
        if (rval) {
            fprintf (stderr, "verify_text_pool failed\n");
            goto CLEANUP;
        }
    } else if (probin) {
        rval = verify_problem (cutname, run_silently);
        if (rval) {
            fprintf (stderr, "verify_problem failed\n");
            goto CLEANUP;
        }
    } else if (poolin) {
        rval = verify_binary_pool (cutname);
        if (rval) {
            fprintf (stderr, "verify_binary_pool failed\n");
            goto CLEANUP;
        }
    }

  CLEANUP:
    printf ("Verification completed (%s) in %.2f seconds\n",
            rval ? "unsuccessful" : "successful",
            CCutil_stop_timer (&z, 0));
    fflush (stdout);
    return rval;
}

static int verify_text_pool (char *poolname, int nodecount)
{
    int rval;
    int *tour = (int *) NULL;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    CCtsp_lpcut_in *cnext;
    int cutcount;
    int nfail = 0;
    int nsuccess = 0;
    int i;

    tour = CC_SAFE_MALLOC (nodecount, int);
    if (tour == (int *) NULL) {
        fprintf (stderr, "Out of memory in verify_text_pool\n");
        rval = 1; goto CLEANUP;
    }
    for (i=0; i<nodecount; i++) tour[i] = i;
    
    rval = CCtsp_file_cuts (poolname, &cuts, &cutcount, nodecount, tour);
    if (rval) {
        fprintf (stderr, "CCtsp_file_cuts failed\n");
        goto CLEANUP;
    }

    while (cuts) {
        cnext = cuts->next;
        rval = CCverify_cut (cuts, CC_TYPE_ALL, (int *) NULL);
        if (rval) {
            fprintf (stderr, "CCverify_cut failed\n");
            nfail++;
        } else {
            nsuccess++;
        }
        CCtsp_free_lpcut_in (cuts);
        CC_FREE (cuts, CCtsp_lpcut_in);
        cuts = cnext;
    }

    printf ("%d of %d cuts failed verification\n", nfail, nfail + nsuccess);

    if (nfail == 0) rval = 0;
    else rval = -1;

  CLEANUP:
    CC_IFFREE (tour, int);
    while (cuts) {
        cnext = cuts->next;
        CCtsp_free_lpcut_in (cuts);
        CC_FREE (cuts, CCtsp_lpcut_in);
        cuts = cnext;
    }
    return rval;
}

static int verify_binary_pool (char *poolname)
{
    int rval;
    int ncount = 0;
    CCtsp_lpcuts *pool = (CCtsp_lpcuts *) NULL;
    int i;
    CCtsp_lpcut_in cut;
    int nfail = 0;
    int nsuccess = 0;

    rval = CCtsp_init_cutpool (&ncount, poolname, &pool);
    if (rval) {
        fprintf (stderr, "CCtsp_init_cutpool failed\n");
        goto CLEANUP;
    }

    for (i=0; i<pool->cutcount; i++) {
        rval = CCtsp_lpcut_to_lpcut_in (pool, &(pool->cuts[i]), &cut);
        if (rval) {
            fprintf (stderr, "CCtsp_lpcut_to_lpcut_in failed\n");
            goto CLEANUP;
        }
        rval = CCverify_cut (&cut, CC_TYPE_ALL, (int *) NULL);
        if (rval) {
            fprintf (stderr, "CCverify_cut failed\n");
            nfail++;
        } else {
            nsuccess++;
        }
        CCtsp_free_lpcut_in (&cut);
    }

    printf ("%d of %d cuts failed verification\n", nfail, nfail + nsuccess);

    if (nfail == 0) rval = 0;
    else rval = -1;

  CLEANUP:
    if (pool) {
        CCtsp_free_cutpool (&pool);
    }
    return rval;
}

static int verify_problem (char *probname, int silent)
{
    int i;
    int rval;
    int nfail = 0;
    int nsuccess = 0;
    CCtsp_lp *lp = (CCtsp_lp *) NULL;
    CCtsp_lpcuts *cuts;
    CCtsp_lpcut_in cut;
    int ncount = 0;

    lp = CC_SAFE_MALLOC (1, CCtsp_lp);
    if (lp == (CCtsp_lp *) NULL) {
        fprintf (stderr, "Out of memory in verify_problem\n");
        rval = 1; goto CLEANUP;
    }
    CCtsp_init_tsp_lp_struct (lp);
    rval = CCtsp_read_probfile (lp, probname, (char *) NULL, &ncount, silent);
    if (rval) {
        fprintf (stderr, "CCtsp_read_probfile failed\n");
        goto CLEANUP;
    }

    rval = CCtsp_build_lpadj (&lp->graph, 0, lp->graph.ecount);
    if (rval) {
        fprintf (stderr, "CCtsp_build_lpadj failed\n");
        goto CLEANUP;
    }
    
    rval = CCtsp_add_branchhistory_to_lp (lp);
    if (rval) {
        fprintf (stderr, "CCtsp_add_branchhistory_to_lp failed\n");
        goto CLEANUP;
    }

    cuts = &lp->cuts;
    for (i=0; i<cuts->cutcount; i++) {
        if (cuts->cuts[i].branch == 0) {
            rval = CCtsp_lpcut_to_lpcut_in (cuts, &(cuts->cuts[i]), &cut);
            if (rval) {
                fprintf (stderr, "CCtsp_lpcut_to_lpcut_in failed\n");
                goto CLEANUP;
            }
            rval = CCverify_cut (&cut, CC_TYPE_ALL, (int *) NULL);
            if (rval) {
                fprintf (stderr, "CCverify_cut failed\n");
                nfail++;
            } else {
                nsuccess++;
            }
            CCtsp_free_lpcut_in (&cut);
        }
    }

    printf ("%d cuts failed verification, %d passed\n", nfail, nsuccess);

    if (nfail == 0) rval = 0;
    else rval = -1;

  CLEANUP:
    CCtsp_free_tsp_lp_struct (&lp);
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "NPt:", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'N': probin = 1; break;
        case 'P': poolin = 1; break;
        case 't': textin = 1; in_ncount = atoi(boptarg); break;
        default: usage (av[0]); return 1;
        }
    }
    if (textin == 0 && probin == 0 && poolin == 0) {
        usage (av[0]);
        return 1;
    }
    if (boptind >= ac) {
        usage (av[0]);
        return 1;
    }
    cutname = av[boptind++];
    if (boptind != ac) {
        usage (av[0]);
        return 1;
    }
    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [NPt:] cutfile\n", fname);
    fprintf (stderr, "   -N        cutfile is problem file\n");
    fprintf (stderr, "   -P        cutfile is a cut pool\n");
    fprintf (stderr, "   -t n      cutfile is text file, n is the nodecount\n");
    fprintf (stderr, "   One of -N, -P, or -t must be specified\n");
}

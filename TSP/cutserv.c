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

#ifdef CC_NETREADY

#include "util.h"
#include "tsp.h"

static char *poolfname = (char *) NULL;
static char *inpoolfname = (char *) NULL;
static unsigned short cutport = CCtsp_CUT_PORT;
static int nodecount = 0;
static int numthreads = 0;
static int seed = 0;

#define SAVE_INTERVAL 500

#define STATUS_NORMAL 0
#define STATUS_SAVE   1
#define STATUS_EXIT   2


static int 
    serve_request (CC_SFILE *f, int ncount, CCtsp_lpcuts *pool, int *status,
        int nthreads, CCrandstate *rstate),
    find_cuts (CC_SFILE *f, int ncount, CCtsp_lpcuts *pool, int nthreads,
        CCrandstate *rstate),
    receive_cuts (CC_SFILE *f, int ncount, CCtsp_lpcuts *pool),
    parseargs (int ac, char **av);

static void
    usage (char *f);


int main (int ac, char **av)
{
    CC_SPORT *p = (CC_SPORT *) NULL;
    CC_SFILE *f = (CC_SFILE *) NULL;
    CCtsp_lpcuts *pool = (CCtsp_lpcuts *) NULL;
    int rval = 0;
    int savecount;
    int status = STATUS_NORMAL;
    CCrandstate rstate;

    seed = (int) CCutil_real_zeit ();
    
    if (parseargs (ac, av))
        return 0;

    CCutil_sprand (seed, &rstate);
    
    CCutil_signal_init ();

    p = CCutil_snet_listen (cutport);
    if (p == (CC_SPORT *) NULL) {
        fprintf (stderr, "CCutil_snet_listen failed\n");
        rval = 1;
        goto CLEANUP;
    }

    rval = CCtsp_init_cutpool (&nodecount, inpoolfname, &pool);
    if (rval) {
        fprintf (stderr, "CCtsp_init_cutpool failed\n");
        rval = 1; goto CLEANUP;
    }

    savecount = pool->cutcount;
    
    do {
        f = CCutil_snet_receive (p);
        if (f == (CC_SFILE *) NULL) {
            fprintf (stderr, "CCutil_snet_receive failed\n");
            continue;
        }
        rval = serve_request (f, nodecount, pool, &status, numthreads,
                              &rstate);
        if (rval) {
            fprintf (stderr, "serve_request failed\n");
            if (CCutil_sclose (f)) {
                fprintf (stderr, "CCutil_sclose failed\n");
            }
            continue;
        }
        rval = CCutil_sclose (f);
        if (rval) {
            fprintf (stderr, "CCutil_sclose failed\n");
        }
        f = (CC_SFILE *) NULL;

        if (pool->cutcount > 0 && (status == STATUS_SAVE ||
                status == STATUS_EXIT || 
                (pool->cutcount - savecount > SAVE_INTERVAL &&
                 pool->cutcount - savecount > pool->cutcount/100))) {
            printf ("Writing pool: %d cuts\n", pool->cutcount);
            fflush (stdout);
            rval = CCtsp_write_cutpool (nodecount, poolfname, pool);
            if (rval) {
                fprintf (stderr, "CCtsp_write_cutpool failed\n");
                goto CLEANUP;
            }
            savecount = pool->cutcount;
        }
    } while (status != STATUS_EXIT);

  CLEANUP:
    if (f != (CC_SFILE *) NULL) CCutil_sclose (f);
    if (p != (CC_SPORT *) NULL) CCutil_snet_unlisten (p);
    return rval;
}

static int serve_request (CC_SFILE *f, int ncount, CCtsp_lpcuts *pool,
        int *status, int nthreads, CCrandstate *rstate)
{
    char request;
    int rval;

    *status = STATUS_NORMAL;

    rval = CCutil_sread_char (f, &request);
    if (rval) {
        fprintf (stderr, "CCutil_sread_char failed\n");
        return rval;
    }
    
    switch (request) {
      case CCtsp_POOL_GETCUTS:
        rval = find_cuts (f, ncount, pool, nthreads, rstate);
        if (rval) {
            fprintf (stderr, "find_cuts failed\n");
            return rval;
        }
        break;
      case CCtsp_POOL_PUTCUTS:
        rval = receive_cuts (f, ncount, pool);
        if (rval) {
            fprintf (stderr, "receive_cuts failed\n");
            return rval;
        }
        break;
      case CCtsp_POOL_SAVECUTS:
        *status = STATUS_SAVE;
        break;
      case CCtsp_POOL_EXIT:
        *status = STATUS_EXIT;
        break;
      default:
        fprintf (stderr, "Invalid request %c\n", request);
        return 1;
    }

    return 0;
}

static int find_cuts (CC_SFILE *f, int ncount, CCtsp_lpcuts *pool,
        int nthreads, CCrandstate *rstate)
{
    int rval;
    int ecount;
    int ncount2;
    int *elist = (int *) NULL;
    double *x = (double *) NULL;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    CCtsp_lpcut_in *cuts_next;
    int cutcount = 0;
    double maxviol;
    double szeit = CCutil_zeit();
    double real_szeit = CCutil_real_zeit();
    int i;

    printf ("Price "); fflush (stdout);
    rval = CCutil_sread_int (f, &ncount2);
    if (rval) {
        fprintf (stderr, "CCutil_sread_int failed\n");
        goto CLEANUP;
    }
    if (ncount != ncount2) {
        fprintf (stderr, "Request for wrong number of nodes (%d != %d)\n",
                 ncount, ncount2);
        rval = 1; goto CLEANUP;
    }
    rval = CCutil_sread_int (f, &ecount);
    if (rval) {
        fprintf (stderr, "CCutil_sread_int failed\n");
        goto CLEANUP;
    }

    elist = CC_SAFE_MALLOC (ecount*2, int);
    x = CC_SAFE_MALLOC (ecount, double);
    if (elist == (int *)    NULL ||
        x     == (double *) NULL) {
        fprintf (stderr, "Out of memory in find_cuts\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<ecount; i++) {
        rval = CCutil_sread_int (f, &elist[2*i]);
        if (rval) {
            fprintf (stderr, "CCutil_sread_int failed\n");
            goto CLEANUP;
        }
        rval = CCutil_sread_int (f, &elist[2*i+1]);
        if (rval) {
            fprintf (stderr, "CCutil_sread_int failed\n");
            goto CLEANUP;
        }
        rval = CCutil_sread_double (f, &x[i]);
        if (rval) {
            fprintf (stderr, "CCutil_sread_double failed\n");
            goto CLEANUP;
        }
    }

    printf ("%d edges...", ecount); fflush (stdout);

    rval = CCtsp_search_cutpool (pool, &cuts, &cutcount, &maxviol, ncount,
            ecount, elist, x, nthreads, rstate);
    if (rval) {
        fprintf (stderr, "CCtsp_search_cutpool failed\n");
        goto CLEANUP;
    }

    printf ("%3d cuts (max %.3f)...", cutcount, maxviol); fflush (stdout);
    
    rval = CCutil_swrite_int (f, cutcount);
    if (rval) {
        fprintf (stderr, "CCutil_swrite_int failed\n");
        goto CLEANUP;
    }

    rval = CCutil_swrite_double (f, maxviol);
    if (rval) {
        fprintf (stderr, "CCutil_swrite_double failed\n");
        goto CLEANUP;
    }

    for (i=0; i<cutcount; i++) {
        rval = CCtsp_write_lpcut_in (f, cuts, ncount);
        if (rval) {
            fprintf (stderr, "CCtsp_write_lpcut_in failed\n");
            goto CLEANUP;
        }
        cuts_next = cuts->next;
        CCtsp_free_lpcut_in (cuts);
        CC_IFFREE (cuts, CCtsp_lpcut_in);
        cuts = cuts_next;
    }

    rval = CCutil_sflush (f);
    f = (CC_SFILE *) NULL;
    if (rval) {
        fprintf (stderr, "CCutil_sflush failed\n");
        goto CLEANUP;
    }

    printf ("sent in %.2f (%.2f real)\n", CCutil_zeit() - szeit,
            CCutil_real_zeit() - real_szeit);
    fflush (stdout);

    rval = 0;

  CLEANUP:
    CC_IFFREE (elist, int);
    CC_IFFREE (x, double);
    while (cuts != (CCtsp_lpcut_in *) NULL) {
        cuts_next = cuts->next;
        CCtsp_free_lpcut_in (cuts);
        CC_IFFREE (cuts, CCtsp_lpcut_in);
        cuts = cuts_next;
    }
    return rval;
}

static int receive_cuts (CC_SFILE *f, int ncount, CCtsp_lpcuts *pool)
{
    int cutcount;
    int rval;
    int ncount2;
    CCtsp_lpcut_in newc;
    double szeit = CCutil_zeit();
    double real_szeit = CCutil_real_zeit();
    int savecount;
    int i;

    newc.cliques = (CCtsp_lpclique *) NULL;
    
    printf ("Cuts  "); fflush (stdout);
    rval = CCutil_sread_int (f, &ncount2);
    if (rval) {
        fprintf (stderr, "CCutil_sread_int failed\n");
        goto CLEANUP;
    }
    if (ncount != ncount2) {
        fprintf (stderr, "Request for wrong number of nodes (%d != %d)\n",
                 ncount, ncount2);
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_sread_int (f, &cutcount);
    if (rval) {
        fprintf (stderr, "CCutil_sread_int failed\n");
        goto CLEANUP;
    }
    printf ("%3d...", cutcount); fflush (stdout);

    savecount = pool->cutcount;
    for (i=0; i<cutcount; i++) {
        rval = CCtsp_read_lpcut_in (f, &newc, ncount);
        if (rval) {
            fprintf (stderr, "CCtsp_read_lpcut_in failed\n");
            goto CLEANUP;
        }
        rval = CCtsp_add_to_cutpool_lpcut_in (pool, &newc);
        if (rval) {
            fprintf (stderr, "CCtsp_add_to_cutpool_lpcut_in failed\n");
            goto CLEANUP;
        }
        CCtsp_free_lpcut_in (&newc);
    }
    printf ("%3d new in %.2f (%.2f real)  TOTAL %d\n", 
            pool->cutcount - savecount, CCutil_zeit() - szeit,
            CCutil_real_zeit() - real_szeit, pool->cutcount);
    fflush (stdout);
    rval = 0;

  CLEANUP:
    CCtsp_free_lpcut_in (&newc);
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "n:p:P:s:T:", &boptind, &boptarg)) != EOF) switch (c) {
    case 'n':
        nodecount = atoi(boptarg);
        break;
    case 'p':
        cutport = atoi(boptarg);
        break;
    case 'P':
        inpoolfname = boptarg;
        break;
    case 's':
        seed = atoi(boptarg);
        break;
    case 'T':
        numthreads = atoi(boptarg);
        break;
    default:
        usage (av[0]);
        return 1;
    }

    if (boptind >= ac) {
        usage (av[0]);
        return 1;
    }

    poolfname = av[boptind++];

    if (boptind != ac) {
        usage (av[0]);
        return 1;
    }

    if (nodecount == 0 && inpoolfname == (char *) NULL) {
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s [-see below-] poolname\n", f);
    fprintf (stderr, "  -n n # of nodes\n");
    fprintf (stderr, "  -p n use port n\n");
    fprintf (stderr, "  -P f initial pool file\n");
    fprintf (stderr, "  -s n use random seed n\n");
    fprintf (stderr, "  -T n use n threads\n");
    fprintf (stderr, "   At least one of -n and -P must be provided\n");
}

#else /* CC_NETREADY */

int main (int ac, char **av)
{
    fprintf (stderr, "Networking code not enabled - not able to serve cuts\n");
    return -1;
}

#endif /* CC_NETREADY */

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

static int debug = 0;
static char *probname = (char *) NULL;
static unsigned short probport = CCtsp_PROB_PORT;
static int run_silently = 0;


static int 
    serve_file (CC_SFILE *f, char *probfname, int silent),
    serve_read (CC_SFILE *f, char *probfname, int id, int silent),
    serve_write (CC_SFILE *f, char *probfname, int id, int silent),
    serve_delete (char *probfname, int id),
    parseargs (int ac, char **av);

static void
    usage (char *f);


int main (int ac, char **av)
{
    CC_SPORT *p = (CC_SPORT *) NULL;
    CC_SFILE *f = (CC_SFILE *) NULL;
    int rval = 0;

    if (parseargs (ac, av))
        return 0;

    CCutil_signal_init ();
    
    if (debug) {
        printf ("Serving files for %s\n", probname); fflush (stdout);
    }
    p = CCutil_snet_listen (probport);
    if (p == (CC_SPORT *) NULL) {
        fprintf (stderr, "CCutil_snet_listen failed\n");
        rval = 1;
        goto CLEANUP;
    }

    for (;;) {
        if (debug) {
            printf ("Waiting for connection\n"); fflush (stdout);
        }
        f = CCutil_snet_receive (p);
        if (f == (CC_SFILE *) NULL) {
            fprintf (stderr, "CCutil_snet_receive failed\n");
            continue;
        }
        if (debug) {
            printf ("Received connection\n"); fflush (stdout);
        }
        rval = serve_file (f, probname, run_silently);
        if (rval) {
            fprintf (stderr, "serve_file failed\n");
            if (CCutil_sclose (f)) {
                fprintf (stderr, "CCutil_sclose failed\n");
            }
            continue;
        }
        if (debug) {
            printf ("Closing connection\n"); fflush (stdout);
        }
        rval = CCutil_sclose (f);
        if (rval) {
            fprintf (stderr, "CCutil_sclose failed\n");
        }
        f = (CC_SFILE *) NULL;
    }

  CLEANUP:
    if (f != (CC_SFILE *) NULL) CCutil_sclose (f);
    if (p != (CC_SPORT *) NULL) CCutil_snet_unlisten (p);
    return rval;
}

static int serve_file (CC_SFILE *f, char *probfname, int silent)
{
    char request;
    char probbuf[1024];
    int id;
    int rval;

    rval = CCutil_sread_char (f, &request);
    if (rval) {
        fprintf (stderr, "CCutil_sread_char failed\n");
        return rval;
    }
    rval = CCutil_sread_string (f, probbuf, sizeof (probbuf));
    if (rval) {
        fprintf (stderr, "CCutil_sread_string failed\n");
        return rval;
    }
    rval = CCutil_sread_int (f, &id);
    if (rval) {
        fprintf (stderr, "CCutil_sread_int failed\n");
        return rval;
    }

    if (strcmp (probfname, probbuf)) {
        fprintf (stderr, "ERROR - serving %s, request %s\n", probfname,
                 probbuf);
        return rval;
    }

    switch (request) {
      case CCtsp_Pread:
        rval = serve_read (f, probfname, id, silent);
        if (rval) {
            fprintf (stderr, "serve_read failed\n");
            return rval;
        }
        break;
      case CCtsp_Pwrite:
        rval = serve_write (f, probfname, id, silent);
        if (rval) {
            fprintf (stderr, "serve_write failed\n");
            return rval;
        }
        break;
      case CCtsp_Pdelete:
        rval = serve_delete (probname, id);
        if (rval) {
            fprintf (stderr, "serve_delete failed\n");
            return rval;
        }
        break;
      default:
        fprintf (stderr, "Invalid request %c\n", request);
        return 1;
    }

    return 0;
}

static int serve_read (CC_SFILE *f, char *probfname, int id, int silent)
{
    CCtsp_PROB_FILE *local = (CCtsp_PROB_FILE *) NULL;
    CCtsp_PROB_FILE *remote = (CCtsp_PROB_FILE *) NULL;
    char request;
    int rval;

    if (debug) {
        printf ("serving read %s %d\n", probfname, id);
        fflush (stdout);
    }
    
    local = CCtsp_prob_read (probfname, id);
    if (local == (CCtsp_PROB_FILE *) NULL) {
        fprintf (stderr, "CCtsp_prob_read failed\n");
        rval = 1; goto CLEANUP;
    }

    remote = CCtsp_prob_server (f);
    if (remote == (CCtsp_PROB_FILE *) NULL) {
        fprintf (stderr, "CCtsp_prob_server failed\n");
        rval = 1; goto CLEANUP;
    }

    for (;;) {
        rval = CCutil_sread_char (f, &request);
        if (rval) {
            fprintf (stderr, "CCutil_sread_char failed\n"); goto CLEANUP;
        }
        if (debug) {
            printf ("Request %c...", request);
            fflush (stdout);
        }
        if (request == CCtsp_Pexit) goto CLEANUP;
        
        rval = CCtsp_prob_copy_section (local, remote, request, silent);
        if (rval) {
            fprintf (stderr, "CCtsp_prob_copy_section failed\n"); goto CLEANUP;
        }
        if (debug) {
            printf ("done\n");
            fflush (stdout);
        }
    }

  CLEANUP:
    if (debug) {
        printf ("exit\n");
        fflush (stdout);
    }
    
    rval |= CCtsp_prob_rclose (local);
    if (rval) {
        fprintf (stderr, "CCtsp_prob_rclose failed\n");
    }
    rval |= CCtsp_prob_wclose (remote);
    if (rval) {
        fprintf (stderr, "CCtsp_prob_wclose failed\n");
    }
    return rval;
}
    
static int serve_write (CC_SFILE *f, char *probfname, int id, int silent)
{
    CCtsp_PROB_FILE *local = (CCtsp_PROB_FILE *) NULL;
    CCtsp_PROB_FILE *remote = (CCtsp_PROB_FILE *) NULL;
    char request;
    int rval;

    if (debug) {
        printf ("serving write %s %d\n", probfname, id);
        fflush (stdout);
    }
    
    local = CCtsp_prob_write (probfname, id);
    if (local == (CCtsp_PROB_FILE *) NULL) {
        fprintf (stderr, "CCtsp_prob_write failed\n");
        rval = 1; goto CLEANUP;
    }

    remote = CCtsp_prob_server (f);
    if (remote == (CCtsp_PROB_FILE *) NULL) {
        fprintf (stderr, "CCtsp_prob_server failed\n");
        rval = 1; goto CLEANUP;
    }

    for (;;) {
        rval = CCutil_sread_char (f, &request);
        if (rval) {
            fprintf (stderr, "CCutil_sread_char failed\n"); goto CLEANUP;
        }
        if (debug) {
            printf ("Request %c...", request);
            fflush (stdout);
        }
        if (request == CCtsp_Pexit) goto CLEANUP;
        
        rval = CCtsp_prob_copy_section (remote, local, request, silent);
        if (rval) {
            fprintf (stderr, "CCtsp_prob_copy_section failed\n"); goto CLEANUP;
        }
        if (debug) {
            printf ("done\n");
            fflush (stdout);
        }
    }

  CLEANUP:
    if (debug) {
        printf ("exit\n");
        fflush (stdout);
    }
    
    rval |= CCtsp_prob_wclose (local);
    if (rval) {
        fprintf (stderr, "CCtsp_prob_wclose failed\n");
    }
    rval |= CCtsp_prob_rclose (remote);
    if (rval) {
        fprintf (stderr, "CCtsp_prob_rclose failed\n");
    }
    return rval;
}

static int serve_delete (char *probfname, int id)
{
    int rval;

    rval = CCtsp_prob_file_delete (probfname, id);
    if (rval) {
        fprintf (stderr, "CCtsp_prob_file_delete failed\n");
        return rval;
    }
    return 0;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "dp:", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'd':
            debug++;
            break;
        case 'p':
            probport = atoi(boptarg);
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

    probname = av[boptind++];

    if (boptind != ac) {
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s [-see below-] probname\n", f);
    fprintf (stderr, "  -d   turn on debugging\n");
    fprintf (stderr, "  -p n use port n\n");
}

#else /* CC_NETREADY */

int main (int ac, char **av)
{
    fprintf (stderr, "Networking code not enabled - not able to serve problems\n");
    return -1;
}

#endif /* CC_NETREADY */

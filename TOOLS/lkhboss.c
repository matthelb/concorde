#include <stdio.h>
#include "util.h"

#define SUB_STAT_OPEN 0
#define SUB_STAT_DONE 1
#define SUB_STAT_WORK 2

static char *indexfname = (char *) NULL;
static char *tourfname = (char *) NULL;
static char *tspfname = (char *) NULL;
static int simpletour = 0;
static int binary_in = 0;
static int tsplib_in = 1;
static int seed = 0;
static int norm = CC_EUCLIDEAN;

static int
    parseargs (int ac, char **av);

static void
    usage (char *fname);

int main (int ac, char **av)
{
    int i, id, ncount, nremain, scount, start, nfail = 0;
    int rval = 0;
    int *tour = (int *) NULL;
    CC_SPORT *lport = (CC_SPORT *) NULL;
    CC_SFILE *s;
    double delta, rtime, newlen, tourlen,  cumtime = 0.0, cumimprove = 0.0;
    double val;
    CCsubdiv_lkh *slist = (CCsubdiv_lkh *) NULL;
    CCsubdiv_lkh *p, *curloc;
    char *problabel = (char *) NULL;
    CCrandstate rstate;
    CCdatagroup dat;

    CCutil_init_datagroup (&dat);
    seed = (int) CCutil_real_zeit ();

    rval = parseargs (ac, av);
    if (rval) return 1;

    CCutil_sprand (seed, &rstate);

    if (tsplib_in) {
        rval = CCutil_gettsplib (tspfname, &ncount, &dat);
        CCcheck_rval (rval, "CCutil_gettsplib failed");
        CCutil_dat_getnorm (&dat, &norm);
    } else {
        rval = CCutil_getdata (tspfname, binary_in, norm, &ncount, &dat,
                               0, 0, &rstate);
        CCcheck_rval (rval, "CCutil_getdata failed");
    }

    if ((norm & CC_NORM_SIZE_BITS) != CC_D2_NORM_SIZE) {
        fprintf (stderr, "Only set up for 2D norms\n");
        rval = 1;  goto CLEANUP;
    }

    if ((norm & CC_NORM_BITS) != CC_KD_NORM_TYPE) {
        fprintf (stderr, "Only set up for KD-tree norms\n");
        rval = 1;  goto CLEANUP;
    }

    printf ("Reading index ..."); fflush (stdout);
    rval = CCutil_read_subdivision_lkh_index (indexfname, &problabel, &i,
                                          &scount, &slist, &tourlen);
    CCcheck_rval (rval, "CCutil_read_subdivision_index failed");
    printf ("DONE\n"); fflush (stdout);

    if (i != ncount) {
        fprintf (stderr, "index file does not match tsp file\n");
        rval = 1;  goto CLEANUP;
    }

    tour = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (tour, "out of memory in main");

    if (simpletour) {
        rval = CCutil_getcycle (ncount, tourfname, tour, 0);
        if (rval) {
            fprintf (stderr, "CCutil_getcycle_tsplib failed\n");
            goto CLEANUP;
        }
    } else {
        rval = CCutil_getcycle_tsplib (ncount, tourfname, tour);
        if (rval) {
            fprintf (stderr, "CCutil_getcycle_tsplib failed\n");
            goto CLEANUP;
        }
    }

    CCutil_cycle_len (ncount, &dat, tour, &val);
    if (val != tourlen) {
        fprintf (stderr, "Cycle length does not match index file\n");
        rval = 1;  goto CLEANUP;
    }

    nremain = 0;

    for (i = 0; i < scount; i++) {
        if (slist[i].newlen <= 0.0) {
            slist[i].status = SUB_STAT_OPEN;
            nremain++;
        } else {
            slist[i].status = SUB_STAT_DONE;
            delta = slist[i].origlen - slist[i].newlen;
            if (delta > 0.0) cumimprove += delta;
        }
    }

    printf ("BEGINNING SUBDIV NET PROCESSING: %s\n\n", problabel);
    fflush (stdout);

    lport = CCutil_snet_listen (CC_SUBDIV_PORT);
    if (lport == (CC_SPORT *) NULL) {                                           
        fprintf (stderr, "CCutil_snet_listen failed\n");
        rval = 1; goto CLEANUP;
    }


    while (nremain) {
        s = CCutil_snet_receive (lport);
        if (!s) {
            fprintf (stderr, "CCutil_snet_receive failed, ignoring\n");
            continue;
        }
        rval = CCutil_sread_int (s, &id);
        if (rval) {
            fprintf (stderr, "CCutil_sread_int failed, abort connection\n");
            rval = 0;
            goto CLOSE_CONN;
        }
        rval = CCutil_sread_double (s, &rtime);
        if (rval) {
            rval = 0;
            fprintf (stderr, "CCutil_sread_double failed, abort connection\n");
            goto CLOSE_CONN;
        }
        rval = CCutil_sread_double (s, &newlen);
        if (rval) {
            rval = 0;
            fprintf (stderr, "CCutil_sread_double failed, abort connection\n");
            goto CLOSE_CONN;
        }
        p = (CCsubdiv_lkh *) NULL;
        if (id == -1) goto GIVE_WORK;
        
        if (id < 0 || id >= scount) {
            fprintf (stderr, "Finished unknown id %d, ignoring\n", id);
            goto GIVE_WORK;
        }
        p = &slist[id];

        if (p->status != SUB_STAT_OPEN && p->status != SUB_STAT_WORK) {
            fprintf (stderr, "Finished completed node %d, ignoring\n", id);
            goto GIVE_WORK;
        }
        if (newlen >= 0.0)  {
            p->newlen = newlen;
            delta = p->origlen - newlen;
            if (delta > 0.0) cumimprove += delta;
        } else {
            printf ("SUBPROBLEM failed\n"); fflush (stdout);
            nfail++;
        }
        p->status = SUB_STAT_DONE;

        rval = CCutil_write_subdivision_lkh_index (problabel, ncount, scount,
                                                   slist, tourlen);
        CCcheck_rval (rval, "CCutil_write_subdivision_index failed");

        cumtime += rtime;
        nremain--;

    GIVE_WORK:

        if (p) {
            printf ("DONE %3d %7.2f sec %4.0f cum %7.0f delta %.2f ", p->id,
                    p->newlen, rtime, cumtime, cumimprove);
        }

        curloc = (CCsubdiv_lkh *) NULL;
        for (i = 0; i < scount; i++) {
            if (slist[i].status == SUB_STAT_OPEN) {
                curloc = &slist[i];
                break;
            }
        }

        if (curloc != (CCsubdiv_lkh *) NULL) {
             rval = CCutil_swrite_int (s, curloc->id);
             CCcheck_rval (rval, "CCutil_swrite_int failed");
        
             rval = CCutil_swrite_string (s, problabel);
             CCcheck_rval (rval, "CCutil_swrite_string failed");

             rval = CCutil_swrite_int (s, norm);
             CCcheck_rval (rval, "CCutil_swrite_int failed");

             rval = CCutil_swrite_int (s, curloc->cnt);
             CCcheck_rval (rval, "CCutil_swrite_int failed");

             start = curloc->start;

             for (i = 0; i < curloc->cnt; i++) {
                 rval = CCutil_swrite_double (s, dat.x[tour[start+i]]); 
                 CCcheck_rval (rval, "CCutil_swrite_double failed");
                 rval = CCutil_swrite_double (s, dat.y[tour[start+i]]); 
                 CCcheck_rval (rval, "CCutil_swrite_double failed");
             }

             printf ("%-4s %2d rem %2d\n",
                     curloc->status == SUB_STAT_OPEN ? "WORK" : "REWK",
                     curloc->id, nremain);
             fflush (stdout);
             curloc->status = SUB_STAT_WORK;
        } else {
             printf ("\n"); fflush (stdout);
             CCutil_swrite_int (s, -1);
        }

    CLOSE_CONN:

        CCutil_sclose (s);
    }
    printf ("\nFINISHED (%d failures): %.2f seconds\n", nfail, cumtime);
    printf ("Total Improvement: %lf\n", cumimprove);
    printf ("New Tour: %lf\n", tourlen - cumimprove);
    fflush (stdout);

CLEANUP:

    CCutil_snet_unlisten (lport);
    CC_IFFREE (slist, CCsubdiv_lkh);
    CC_IFFREE (problabel, char);
    CC_IFFREE (tour, int);
    CCutil_freedatagroup (&dat);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c, inorm;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "btN:?", &boptind, &boptarg)) != EOF) { 
        switch (c) {
        case 'b':
            binary_in = 1;
            break;
        case 't':
            simpletour = 1;
            break;
        case 'N':
            inorm = atoi (boptarg);
            switch (inorm) {
            case 0: norm = CC_MAXNORM; break;
            case 1: norm = CC_MANNORM; break;
            case 2: norm = CC_EUCLIDEAN; break;
            case 3: norm = CC_EUCLIDEAN_3D; break;
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
    }

    if (boptind < ac) {
        indexfname = av[boptind++];
        if (boptind < ac) {
            tourfname = av[boptind++];
            if (boptind < ac) {
                tspfname = av[boptind++];
            } else {
                fprintf (stderr, "Missing one file\n");
                usage (av[0]);
                return 1;
            }
        } else {
            fprintf (stderr, "Missing two files\n");
            usage (av[0]);
            return 1;
        }
    } else {
        fprintf (stderr, "Missing indexfile, tourfile, and tspfile\n");
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] index_file tour_file tsp_or_dat_file\n", fname);
    fprintf (stderr, "   -b    datfile in integer binary format\n");
    fprintf (stderr, "   -t    tour file in concorde format (default TSPLIB)\n");
    fprintf (stderr, "   -N #  norm (must specify if dat file is not a TSPLIB file)\n");
    fprintf (stderr, "         0=MAX, 1=L1, 2=L2, 17=GEOM, 18=JOHNSON\n");
}

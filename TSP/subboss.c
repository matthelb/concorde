#include <stdio.h>
#include "util.h"

#define SUB_STAT_OPEN 0
#define SUB_STAT_DONE 1
#define SUB_STAT_WORK 2

static int idcmp (void *v1, void *v2, CC_UNUSED void *u_data)
{
    return ((long) v1) - ((long) v2);
}

static unsigned int idhash (void *v1, CC_UNUSED void *u_data)
{
    return (unsigned int) (unsigned long) v1;
}


int main (int ac, char **av)
{
    int i, id, ncount, nremain, scount, nfail = 0;
    int rval = 0;
    CCgenhash idmap;
    int genhash_created = 0;
    CC_SPORT *lport = (CC_SPORT *) NULL;
    CC_SFILE *s;
    double rtime, rbound, cumtime = 0.0, cumbound = 0.0;
    CCsubdiv *slist = (CCsubdiv *) NULL;
    CCsubdiv *p, *curloc;
    char *problabel = (char *) NULL;

    if (ac != 2) {
        fprintf (stderr, "Usage: %s index_file\n", av[0]);
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_read_subdivision_index (av[1], &problabel, &ncount,
                                          &scount, &slist);
    CCcheck_rval (rval, "CCutil_read_subdivision_index failed");

    for (i = 0; i < scount; i++) {
        if (slist[i].bound <= 0.0) slist[i].status = SUB_STAT_OPEN;
        else                       slist[i].status = SUB_STAT_DONE;
    }

    printf ("BEGINNING SUBDIV NET PROCESSING: %s\n\n", problabel);
    fflush (stdout);

    rval = CCutil_genhash_init (&idmap, scount, idcmp, idhash, NULL, 0.8,
                                0.4);
    CCcheck_rval (rval, "CCutil_genhash_init failed");
    genhash_created = 1;

    for (i=0; i < scount; i++) {
        rval = CCutil_genhash_insert (&idmap, (void *) slist[i].id,
                                      (void *) &slist[i]);
        CCcheck_rval (rval, "CCutil_genhash_insert failed");
    }
    nremain = scount;

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
        rval = CCutil_sread_double (s, &rbound);
        if (rval) {
            rval = 0;
            fprintf (stderr, "CCutil_sread_double failed, abort connection\n");
            goto CLOSE_CONN;
        }
        p = (CCsubdiv *) NULL;
        if (id == -1) goto GIVE_WORK;
        p = (CCsubdiv *) CCutil_genhash_lookup (&idmap, (void *) id);
        if (!p) {
            fprintf (stderr, "Finished unknown id %d, ignoring\n", id);
            goto GIVE_WORK;
        }
        if (p->status != SUB_STAT_OPEN && p->status != SUB_STAT_WORK) {
            fprintf (stderr, "Finished completed node %d, ignoring\n", id);
            goto GIVE_WORK;
        }
        if (rbound >= 0.0)  {
            p->bound = rbound;
            cumbound += rbound;
        } else {
            printf ("SUBPROBLEM failed\n"); fflush (stdout);
            nfail++;
        }
        p->status = SUB_STAT_DONE;

        rval = CCutil_write_subdivision_index (problabel, ncount, scount,
                                               slist);
        CCcheck_rval (rval, "CCutil_write_subdivision_index failed");

        cumtime += rtime;
        nremain--;

    GIVE_WORK:

        if (p) {
            printf ("DONE %3d %7.2f sec %4.0f cum %7.0f bnd %.2f ", p->id,
                    p->bound, rtime, cumtime, cumbound);
        }

        curloc = (CCsubdiv *) NULL;
        for (i = 0; i < scount; i++) {
            if (slist[i].status == SUB_STAT_OPEN) {
                curloc = &slist[i];
                break;
            }
        }

        if (curloc != (CCsubdiv *) NULL) {
             int t_ncount;
             CCdatagroup t_dat;
             int *t_perm = (int *) NULL;
             char buf[1024];

             CCutil_init_datagroup (&t_dat);

             rval = CCutil_swrite_int (s, curloc->id);
             CCcheck_rval (rval, "CCutil_swrite_int failed");
        
             rval = CCutil_swrite_string (s, problabel);
             CCcheck_rval (rval, "CCutil_swrite_string failed");

             sprintf (buf, "%s_%d.mas", problabel, curloc->id);

             rval = CCutil_getmaster (buf, &t_ncount, &t_dat, &t_perm);
             CCcheck_rval (rval, "CCutil_getmaster failed");

             rval = CCutil_writemaster (s, t_ncount, &t_dat, t_perm);
             CCcheck_rval (rval, "CCutil_writemaster failed");

             printf ("%-4s %2d rem %2d\n",
                     curloc->status == SUB_STAT_OPEN ? "WORK" : "REWK",
                     curloc->id, nremain);
             fflush (stdout);
             curloc->status = SUB_STAT_WORK;


             CC_IFFREE (t_perm, int);
             CCutil_freedatagroup (&t_dat);
        } else {
             printf ("\n"); fflush (stdout);
             CCutil_swrite_int (s, -1);
        }

    CLOSE_CONN:

        CCutil_sclose (s);
    }
    printf ("\nFINISHED (%d failures): %.2f seconds\n", nfail, cumtime);
    printf ("Lower bound: %lf\n", cumbound);
    fflush (stdout);

CLEANUP:

    CCutil_snet_unlisten (lport);
    if (genhash_created) {
      CCutil_genhash_free (&idmap, NULL);
      genhash_created = 0;
    }
    CC_IFFREE (slist, CCsubdiv);
    CC_IFFREE (problabel, char);

    return rval;
}

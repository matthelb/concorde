#include <stdio.h>
#include "util.h"

#define SUB_WIDE      /* SUB_LOCAL or SUB_WIDE */

#define SUB_STAT_OPEN 0
#define SUB_STAT_DONE 1
#define SUB_STAT_WORK 2

#ifdef SUB_LOCAL
#define SUBBOSS_PORT   CC_SUBDIV_PORT
#define SUBGRUNT_PORT  CC_SUBGATE_PORT
#endif

#ifdef SUB_WIDE
#define SUBBOSS_PORT   CC_SUBGATE_PORT
#define SUBGRUNT_PORT  CC_SUBDIV_PORT
#endif

int main (int ac, char **av)
{
    int new_id, id, nprocessed = 0, rval = 0;
    char new_probname[CCutil_FILE_NAME_LEN];
    CC_SPORT *lport = (CC_SPORT *) NULL;
    CC_SFILE *s     = (CC_SFILE *) NULL;
    CC_SFILE *sgate = (CC_SFILE *) NULL;
    double rtime, rbound, cumtime = 0.0;
    char *bosshost = (char *) NULL;

    if (ac != 2) {
        fprintf (stderr, "Usage: %s home_gateway\n", av[0]);
        rval = 1; goto CLEANUP;
    }

    bosshost = av[1];


    printf ("BEGINNING Gateway NET PROCESSING\n\n"); fflush (stdout);


    lport = CCutil_snet_listen (SUBGRUNT_PORT);
    if (lport == (CC_SPORT *) NULL) {                                           
        fprintf (stderr, "CCutil_snet_listen failed\n");
        rval = 1; goto CLEANUP;
    }


    while (1) {
        s = CCutil_snet_receive (lport);
        if (!s) {
            fprintf (stderr, "CCutil_snet_receive failed, ignoring\n");
            continue;
        }
        rval = CCutil_sread_int (s, &id);
        CCcheck_rval (rval, "CCutil_sread_int failed");

        rval = CCutil_sread_double (s, &rtime);
        CCcheck_rval (rval, "CCutil_sread_double failed");

        rval = CCutil_sread_double (s, &rbound);
        CCcheck_rval (rval, "CCutil_sread_double failed");

        sgate = CCutil_snet_open (bosshost, SUBBOSS_PORT);
        if (!sgate) {
            fprintf (stderr, "CCutil_snet_open failed\n");
            rval = 1;  goto CLEANUP;
        }

        rval = CCutil_swrite_int (sgate, id);
        CCcheck_rval (rval, "CCutil_swrite_int failed (id)");

        rval = CCutil_swrite_double (sgate, rtime);
        CCcheck_rval (rval, "CCutil_swrite_double failed (rtime)");

        rval = CCutil_swrite_double (sgate, rbound);
        CCcheck_rval (rval, "CCutil_swrite_double failed (rtime)");

        
        if (id != -1) {
            nprocessed++;
            cumtime += rtime;
            printf ("DONE %3d %7.2f sec %4.0f cum %7.0f avg %.2f ", id,
                    rbound, rtime, cumtime, cumtime / (double) nprocessed);
        }

        rval = CCutil_sread_int (sgate, &new_id);
        CCcheck_rval (rval, "CCutil_sread_int failed (new_id)");

        if (new_id != -1) {
             int t_ncount;
             CCdatagroup t_dat;
             int *t_perm = (int *) NULL;

             CCutil_init_datagroup (&t_dat);

             rval = CCutil_swrite_int (s, new_id);
             CCcheck_rval (rval, "CCutil_swrite_int failed");

             rval = CCutil_sread_string (sgate, new_probname,
                                         CCutil_FILE_NAME_LEN);
             CCcheck_rval (rval, "CCutil_sread_string failed (new_probname)");
        
             rval = CCutil_swrite_string (s, new_probname);
             CCcheck_rval (rval, "CCutil_swrite_string failed");

             rval = CCutil_readmaster (sgate, &t_ncount, &t_dat, &t_perm);
             CCcheck_rval (rval, "CCutil_readmaster failed");

             rval = CCutil_writemaster (s, t_ncount, &t_dat, t_perm);
             CCcheck_rval (rval, "CCutil_writemaster failed");

             printf ("WORK %2d\n", new_id);
             fflush (stdout);

             CC_IFFREE (t_perm, int);
             CCutil_freedatagroup (&t_dat);
        } else {
             CCutil_swrite_int (s, -1);
             printf ("STOP\n"); fflush (stdout);
        }

        if (sgate) CCutil_sclose (sgate);
        if (s) CCutil_sclose (s);
    }

CLEANUP:

    CCutil_snet_unlisten (lport);
    return rval;
}

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
/*                           KILL A GRUNT                                   */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: October 29, 1999                                                  */
/*                                                                          */
/*  This program is used to kill a grunt, and inform the boss of a          */
/*  parallel branching run that the grunt is dead.  It reads the grunt's    */
/*  log from stdin, and from that determines the boss, the current task,    */
/*  and the grunt's pid.                                                    */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"

static int kill_process = 1;

static void
    parse_log (char **p_boss, int *p_task, int *p_process),
    usage (char *f);

static int
    parseargs (int ac, char **av);

extern int
    kill (pid_t pid, int sig);

int main (int ac, char **av)
{
    int rval;
    char *boss = (char *) NULL;
    int task;
    int process;
    CC_SFILE *f = (CC_SFILE *) NULL;

    rval = parseargs (ac, av);
    if (rval) {
        goto CLEANUP;
    }
    
    parse_log (&boss, &task, &process);

    if (kill_process) {
        printf ("killing process %d, task %d, and reporting to %s\n",
                process, task, boss);
        fflush (stdout);
    } else {
        printf ("reporting dead task %d to %s\n", task, boss);
        fflush (stdout);
    }

    if (kill_process) {
        rval = kill ((pid_t) process, SIGTERM);
        if (rval) {
            perror ("kill");
            fprintf (stderr, "Unable to kill process %d\n", process);
/*
            if (errno == ESRCH) {   /* Does not work on Red Hat 8 */
                fprintf (stderr, "Process does not exist, telling boss anyway\n");
            } else {
                goto CLEANUP;
            }
*/
            fprintf (stderr, "Process does not exist, telling boss anyway\n");
        }
    }

#ifdef CC_NETREADY
    if (boss != (char *) NULL && task >= 0) {
        f = CCutil_snet_open (boss, CCtsp_HOST_PORT);
        if (f == (CC_SFILE *) NULL) {
            fprintf (stderr, "Could not open connection to host %s\n", boss);
            rval = 1; goto CLEANUP;
        }
        rval = CCutil_swrite_char (f, CCtsp_BBREQ_DEADNODE);
        if (rval) {
            fprintf (stderr, "CCutil_swrite_char failed\n");
            goto CLEANUP;
        }
        rval = CCutil_swrite_int (f, task);
        if (rval) {
            fprintf (stderr, "CCutil_swrite_int failed\n");
            goto CLEANUP;
        }

        rval = CCutil_sclose (f);
        f = (CC_SFILE *) NULL;
        if (rval) {
            fprintf (stderr, "CCutil_sclose failed\n");
            goto CLEANUP;
        }
    }
    rval = 0;

#else /* CC_NETREADY */

    fprintf (stderr, "Networking disabled, unable to tell boss\n");
    rval = 1;

#endif
  CLEANUP:
    if (f != (CC_SFILE *) NULL) {
        CCutil_sclose (f);
    }
    CC_IFFREE (boss, char);
    return rval;
}

static void parse_log (char **p_boss, int *p_task, int *p_process)
{
    char logline[1024];
    char logline2[1024];
    char *curline = logline;
    char *prevline = (char *) NULL;
    char *altline = logline2;
    char *p;

    *p_boss = (char *) NULL;
    *p_task = -1;
    *p_process = -1;
    
    while (fgets (curline, sizeof (logline), stdin) != (char *) NULL) {
        if (!strncmp (curline, "Host: ", 6)) {
            if (prevline) {
                for (p=prevline; *p; p++) {
                    if (!strncmp (p, " -g", 3)) {
                        p = p+3;
                        while (*p == ' ' || *p == '\t' || *p == '\n' ||
                               *p == '\r') p++;
                        CC_IFFREE (*p_boss, char);
                        *p_boss = CCutil_strdup2 (p);
                        break;
                    }
                }
            }
            for (p = curline; *p; p++) {
                if (!strncmp (p, " id: ", 5)) {
                    *p_process = atoi (p+5);
                    break;
                }
            }
        } else if (!strncmp (curline, "Task ", 4)) {
            for (p = curline; *p; p++) {
                if (!strncmp (p, " node ", 6)) {
                    *p_task = atoi (p+6);
                    break;
                }
            }
        }
        prevline = curline;
        curline = altline;
        altline = prevline;
    }
}

static int parseargs (int ac, char **av)
{
    int boptind = 1;
    int c;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "k", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'k':
            kill_process = 0;
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }
    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s [- see below -] < log_file\n", f);
    fprintf (stderr, "  -k   do not kill process\n");
}

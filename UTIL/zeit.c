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
/*                        TIMING FUNCTIONS                                  */
/*                                                                          */
/*                            TSP CODE                                      */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: Summer 1994  (cofeb16)                                            */
/*        December 1997 (dla)                                               */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  double CCutil_zeit (void)                                               */
/*        - To measure cpu time.                                            */
/*    To use this, set double t = CCutil_zeit (), run the function you      */
/*    want to time, then compute CCutil_zeit () - t.                        */
/*                                                                          */
/*  double CCutil_real_zeit (void)                                          */
/*    - To measure wall clock time.                                         */
/*                                                                          */
/*    To use this, set double t = CCutil_real_zeit (), run the function     */
/*    you want to time, then compute CCutil_real_zeit () - t.               */
/*                                                                          */
/*  void CCutil_init_timer (CCutil_timer *t, const char *name)              */
/*    - Initializes a CCutil_timer, and gives it a name.                    */
/*    - The name is silently truncated if it is too long.                   */
/*                                                                          */
/*  void CCutil_start_timer (CCutil_timer *t)                               */
/*    - Starts the timer.                                                   */
/*                                                                          */
/*  void CCutil_suspend_timer (CCutil_timer *t)                             */
/*    - Suspends the timer.  Similar to CCutil_stop_timer, but doesn't      */
/*      count a call, and doesn't output.                                   */
/*                                                                          */
/*  void CCutil_resume_timer (CCutil_timer *t)                              */
/*    - Resumes the timer after a suspend.                                  */
/*                                                                          */
/*  double CCutil_stop_timer (CCutil_timer *t, int printit)                 */
/*    - Stops the timer, and returns the time since the last start.         */
/*    - if printit == 1, outputs the time spent.                            */
/*    - if printit == 2, outputs the time spent only if nonzero             */
/*    - if printit == 3,4, like 1,2, except brief, table-form output        */
/*                                                                          */
/*  double CCutil_total_timer (CCutil_timer *t, int printit)                */
/*    - Returns the cumulative time for this timer.                         */
/*    - if printit == 1, outputs the cumulative time.                       */
/*    - if printit == 2, outputs the cumulative time only if nonzero        */
/*    - if printit == 3,4, like 1,2, except brief, table-form output        */
/*                                                                          */
/****************************************************************************/


#include "machdefs.h"
#include "util.h"

#ifdef HAVE_GETRUSAGE

#ifdef HAVE_SYS_RESOURCE_H
# include <sys/resource.h>
#endif

double CCutil_zeit (void)
{
    struct rusage ru;

    getrusage (RUSAGE_SELF, &ru);

    return ((double) ru.ru_utime.tv_sec) +
            ((double) ru.ru_utime.tv_usec) / 1000000.0;
}
#else /* HAVE_GETRUSAGE */

#ifdef HAVE_TIMES

#ifdef HAVE_SYS_PARAM_H
# include <sys/param.h>
#endif
#ifdef HAVE_SYS_TIMES_H
# include <sys/times.h>
#endif

#ifdef CLK_TCK
#define MACHINE_FREQ CLK_TCK
#else
#define MACHINE_FREQ HZ
#endif

double CCutil_zeit (void)
{
    struct tms now;

    times (&now);
    return ((double) now.tms_utime) / ((double) MACHINE_FREQ);
}
#else /* HAVE_TIMES */

#ifdef HAVE_CLOCK

#ifndef CLOCKS_PER_SEC
#ifdef CLK_TCK
#define CLOCKS_PER_SEC CLK_TCK
#else
#define CLOCKS_PER_SEC 60
#endif
#endif

double CCutil_zeit (void)
{
    return ((double) clock()) / ((double) CLOCKS_PER_SEC);
}

#else /* HAVE_CLOCK */

double CCutil_zeit (void)
{
    return 0.0;
}
#endif /* HAVE_CLOCK */
#endif /* HAVE_TIMES */
#endif /* HAVE_GETRUSAGE */

double CCutil_real_zeit (void)
{
    return (double) time (0);
}

void CCutil_init_timer (CCutil_timer *t, const char *name)
{
    t->szeit    = -1.0;
    t->cum_zeit = 0.0;
    t->count    = 0;
    if (name == (char *) NULL || name [0] == '\0') {
        strncpy (t->name, "ANONYMOUS", sizeof (t->name)-1);
    } else {
        strncpy (t->name, name, sizeof (t->name)-1);
    }
    t->name[sizeof (t->name) - 1] = '\0';
}

void CCutil_start_timer (CCutil_timer *t)
{
    if (t->szeit != -1.0) {
        fprintf (stderr, "Warning: restarting running timer %s\n", t->name);
    }
    t->szeit = CCutil_zeit ();
}

void CCutil_suspend_timer (CCutil_timer *t)
{
    if (t->szeit == -1.0) {
        fprintf (stderr, "Warning: suspended non-running timer %s\n", t->name);
        return;
    }
    
    t->cum_zeit += CCutil_zeit() - t->szeit;
    t->szeit = -1.0;
}

void CCutil_resume_timer (CCutil_timer *t)
{
    if (t->szeit != -1.0) {
        fprintf (stderr, "Warning: resuming running timer %s\n", t->name);
        return;
    }
    t->szeit = CCutil_zeit ();
}

double CCutil_stop_timer (CCutil_timer *t, int printit)
{
    double z;
    
    if (t->szeit == -1.0) {
        fprintf (stderr, "Warning: stopping non-running timer %s\n", t->name);
        return 0.0;
    }
    z = CCutil_zeit() - t->szeit;
    t->szeit = -1.0;
    t->cum_zeit += z;
    t->count++;
    if (printit == 1 || (printit == 2 && z > 0.0)) {
        printf ("Time for %s: %.2f seconds (%.2f total in %d calls)\n",
                t->name, z, t->cum_zeit, t->count);
        fflush (stdout);
    } else if (printit == 3 || (printit == 4 && z > 0.0)) {
        printf ("T %-34.34s %9.2f %9.2f %d\n",
                t->name, z, t->cum_zeit, t->count);
        fflush (stdout);
    }
    return z;
}

double CCutil_total_timer (CCutil_timer *t, int printit)
{
    double z = t->cum_zeit;

    if (t->szeit != -1.0) z += CCutil_zeit() - t->szeit;
    if (printit == 1 || (printit == 2 && z > 0.0)) {
        printf ("Total time for %-34.34s %.2f seconds in %d%s calls\n",
                t->name, z, t->count, t->szeit == -1.0 ? "" : "+1");
        fflush (stdout);
    } else if (printit == 3 || (printit == 4 && z > 0.0)) {
        printf ("CT %-34.34s %9.2f %6d%s\n",
                t->name, z, t->count, t->szeit == -1.0 ? "" : "+1");
        fflush (stdout);
    }
    return z;
}

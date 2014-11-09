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
/*                  SIGNAL HANDLING ROUTINES                                */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/* Written by:  Applegate, Bixby, Chvatal, and Cook                         */
/* Date: October 11, 1997                                                   */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_signal_handler (int ccsignum, CCutil_handler handler)        */
/*     -ccsignum is one of the CCutil_SIG* values defined in util.h         */
/*     -handler is the new signal handler.                                  */
/*    INSTALLS handler as the handler for ccsignum signals.                 */
/*                                                                          */
/*    CCutil_handler is just a typedef for a pointer to a signal            */
/*    handling function.  It is a void (*)(int).  The signal handler        */
/*    is called whenever the signal is raised.  The integer argument        */
/*    passed to the signal handler is the signal number, in the             */
/*    operating system's numbering.  CCutil_sig_to_ccsig can be used        */
/*    to convert this signal number into a CCutil_SIG* value.               */
/*    CCutil_handler_fatal and CCutil_handler_warn are provided as          */
/*    signal handlers for two common cases.                                 */
/*                                                                          */
/*  int CCutil_signal_default (int ccsignum)                                */
/*     -ccsignum is one of the CCutil_SIG* values defined in util.h         */
/*    RESTORES the default handling for ccsignum signals.                   */
/*                                                                          */
/*  int CCutil_signal_ignore (int ccsignum)                                 */
/*     -ccsignum is one of the CCutil_SIG* values defined in util.h         */
/*    IGNORES ccsignum signals.                                             */
/*                                                                          */
/*  int CCutil_sig_to_ccsig (int signum)                                    */
/*    CONVERTS a signal number from the operating system's terms to         */
/*    a CCutil_SIG* value.  If there is no corresponding value, returns     */
/*    -1.                                                                   */
/*                                                                          */
/*  void CCutil_signal_init (void)                                          */
/*    INITIALIZES signal handlers to CCutil_handler_fatal,                  */
/*    CCutil_handler_warn or CCutil_handler_exit as appropriate.            */
/*                                                                          */
/*  void CCutil_handler_fatal (int signum)                                  */
/*    HANDLES a signal by printing an error message, waiting a long         */
/*    time to provide an opportunity to examine the program state with      */
/*    a debugger, and then exiting.  Designed to be used as a handler       */
/*    argument to CCutil_signal_handler.                                    */
/*                                                                          */
/*  void CCutil_handler_warn (int signum)                                   */
/*    HANDLES a signal by printing a warning message and continuing.        */
/*    Designed to be used as a handler argument to CCutil_signal_handler.   */
/*                                                                          */
/*  void CCutil_handler_exit (int signum)                                   */
/*    HANDLES a signal by printing a warning message and exiting.           */
/*    Designed to be used as a handler argument to CCutil_signal_handler.   */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* NOTES: One and only one of CCSIGNAL_SIGACTION, CCSIGNAL_SIGNAL,          */
/*  or CCSIGNAL_NONE should be defined to select the signal                 */
/*  interface used by the operating system.  CCSIGNAL_NONE does nothing     */
/*  to handle signals.                                                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"


static int
    ccsig_to_sig (int signum);

static const char
   *ccsignal_name (int ccsignum);


#ifdef CCSIGNAL_SIGACTION

int CCutil_signal_handler (int ccsignum, CCutil_handler handler)
{
    int signum = ccsig_to_sig (ccsignum);
    struct sigaction new;
    int rval;

    if (signum == -1) {
#if 0
        fprintf (stderr, "Signal %s doesn't exist on this system\n",
                 ccsignal_name (ccsignum));
#endif
        return -1;
    }

    rval = sigemptyset (&new.sa_mask);
    if (rval) {
        perror ("sigemptyset");
        fprintf (stderr, "sigemptyset failed for signal %s\n",
                 ccsignal_name (ccsignum));
        return -1;
    }
    new.sa_handler = handler;
    new.sa_flags = 0;
#ifdef SA_RESTART
    new.sa_flags |= SA_RESTART;
#endif
    rval = sigaction (signum, &new, (struct sigaction *) NULL);
    if (rval) {
        perror ("sigaction");
        fprintf (stderr, "Sigaction for signal %s failed\n",
                 ccsignal_name (ccsignum));
        return -1;
    }

    return 0;
}

#endif /* CCSIGNAL_SIGACTION */

#ifdef CCSIGNAL_SIGNAL

int CCutil_signal_handler (int ccsignum, CCutil_handler handler)
{
    int signum = ccsig_to_sig (ccsignum);
    CCutil_handler sigval;

    if (signum == -1) {
#if 0
        fprintf (stderr, "Signal %s doesn't exist on this system\n",
                 ccsignal_name (ccsignum));
#endif
        return -1;
    }

    sigval = signal (signum, handler);
    if (sigval == SIG_ERR) {
        perror ("signal");
        fprintf (stderr, "Signal() for signal %s failed\n",
                 ccsignal_name (ccsignum));
        return -1;
    }

    return 0;
}

#endif /* CCSIGNAL_SIGNAL */

#ifdef CCSIGNAL_NONE

int CCutil_signal_handler (int ccsignum, CCutil_handler handler)
{
    if (ccsignum || handler) {
        fprintf (stderr, "No signal handling enabled\n");
    }
    return -1;
}

#endif /* CCSIGNAL_NONE */

int CCutil_signal_default (int ccsignum)
{
    int rval;

    rval = CCutil_signal_handler (ccsignum, SIG_DFL);
    if (rval) {
        fprintf (stderr, "CCutil_signal_handler failed\n");
        return rval;
    }
    return 0;
}

int CCutil_signal_ignore (int ccsignum)
{
    int rval;

    rval = CCutil_signal_handler (ccsignum, SIG_IGN);
    if (rval) {
        fprintf (stderr, "CCutil_signal_handler failed\n");
        return rval;
    }
    return 0;
}

void CCutil_signal_init (void)
{
    (void) CCutil_signal_handler (CCutil_SIGHUP,     CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGINT,     CCutil_handler_exit);
    (void) CCutil_signal_handler (CCutil_SIGQUIT,    CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGILL,     CCutil_handler_fatal);
    (void) CCutil_signal_handler (CCutil_SIGTRAP,    CCutil_handler_fatal);
    (void) CCutil_signal_handler (CCutil_SIGABRT,    CCutil_handler_fatal);
    (void) CCutil_signal_handler (CCutil_SIGEMT,     CCutil_handler_fatal);
    (void) CCutil_signal_handler (CCutil_SIGFPE,     CCutil_handler_fatal);
    (void) CCutil_signal_handler (CCutil_SIGBUS,     CCutil_handler_fatal);
    (void) CCutil_signal_handler (CCutil_SIGSEGV,    CCutil_handler_fatal);
    (void) CCutil_signal_handler (CCutil_SIGSYS,     CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGPIPE,    CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGALRM,    CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGTERM,    CCutil_handler_exit);
    (void) CCutil_signal_handler (CCutil_SIGUSR1,    CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGUSR2,    CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGCHLD,    CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGPWR,     CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGWINCH,   CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGURG,     CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGIO,      CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGTSTP,    CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGCONT,    CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGTTIN,    CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGTTOU,    CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGVTALRM,  CCutil_handler_warn);
    /* If profiling under some systems, you may need to comment out the */
    /* following line (concerning SIGPROF)                              */
/*
    (void) CCutil_signal_handler (CCutil_SIGPROF,    CCutil_handler_warn);
*/
    (void) CCutil_signal_handler (CCutil_SIGXCPU,    CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGXFSZ,    CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGSTKFLT,  CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGIOT,     CCutil_handler_fatal);
    (void) CCutil_signal_handler (CCutil_SIGPOLL,    CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGMSG,     CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGDANGER,  CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGMIGRATE, CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGPRE,     CCutil_handler_warn);
    (void) CCutil_signal_handler (CCutil_SIGVIRT,    CCutil_handler_warn);
}


static int ccsig_to_sig (int ccsignum)
{
    switch (ccsignum) {
      case CCutil_SIGHUP:
#ifdef SIGHUP
        return SIGHUP;
#else
        return -1;
#endif
      case CCutil_SIGINT:
#ifdef SIGINT
        return SIGINT;
#else
        return -1;
#endif
      case CCutil_SIGQUIT:
#ifdef SIGQUIT
        return SIGQUIT;
#else
        return -1;
#endif
      case CCutil_SIGILL:
#ifdef SIGILL
        return SIGILL;
#else
        return -1;
#endif
      case CCutil_SIGTRAP:
#ifdef SIGTRAP
        return SIGTRAP;
#else
        return -1;
#endif
      case CCutil_SIGABRT:
#ifdef SIGABRT
        return SIGABRT;
#else
        return -1;
#endif
      case CCutil_SIGEMT:
#ifdef SIGEMT
        return SIGEMT;
#else
        return -1;
#endif
      case CCutil_SIGFPE:
#ifdef SIGFPE
        return SIGFPE;
#else
        return -1;
#endif
      case CCutil_SIGKILL:
#ifdef SIGKILL
        return SIGKILL;
#else
        return -1;
#endif
      case CCutil_SIGBUS:
#ifdef SIGBUS
        return SIGBUS;
#else
        return -1;
#endif
      case CCutil_SIGSEGV:
#ifdef SIGSEGV
        return SIGSEGV;
#else
        return -1;
#endif
      case CCutil_SIGSYS:
#ifdef SIGSYS
        return SIGSYS;
#else
        return -1;
#endif
      case CCutil_SIGPIPE:
#ifdef SIGPIPE
        return SIGPIPE;
#else
        return -1;
#endif
      case CCutil_SIGALRM:
#ifdef SIGALRM
        return SIGALRM;
#else
        return -1;
#endif
      case CCutil_SIGTERM:
#ifdef SIGTERM
        return SIGTERM;
#else
        return -1;
#endif
      case CCutil_SIGUSR1:
#ifdef SIGUSR1
        return SIGUSR1;
#else
        return -1;
#endif
      case CCutil_SIGUSR2:
#ifdef SIGUSR2
        return SIGUSR2;
#else
        return -1;
#endif
      case CCutil_SIGCHLD:
#ifdef SIGCHLD
        return SIGCHLD;
#else
#ifdef SIGCLD
        return SIGCLD;
#else
        return -1;
#endif
#endif
      case CCutil_SIGPWR:
#ifdef SIGPWR
        return SIGPWR;
#else
        return -1;
#endif
      case CCutil_SIGWINCH:
#ifdef SIGWINCH
        return SIGWINCH;
#else
        return -1;
#endif
      case CCutil_SIGURG:
#ifdef SIGURG
        return SIGURG;
#else
        return -1;
#endif
      case CCutil_SIGIO:
#ifdef SIGIO
        return SIGIO;
#else
        return -1;
#endif
      case CCutil_SIGSTOP:
#ifdef SIGSTOP
        return SIGSTOP;
#else
        return -1;
#endif
      case CCutil_SIGTSTP:
#ifdef SIGTSTP
        return SIGTSTP;
#else
        return -1;
#endif
      case CCutil_SIGCONT:
#ifdef SIGCONT
        return SIGCONT;
#else
        return -1;
#endif
      case CCutil_SIGTTIN:
#ifdef SIGTTIN
        return SIGTTIN;
#else
        return -1;
#endif
      case CCutil_SIGTTOU:
#ifdef SIGTTOU
        return SIGTTOU;
#else
        return -1;
#endif
      case CCutil_SIGVTALRM:
#ifdef SIGVTALRM
        return SIGVTALRM;
#else
        return -1;
#endif
      case CCutil_SIGPROF:
#ifdef SIGPROF
        return SIGPROF;
#else
        return -1;
#endif
      case CCutil_SIGXCPU:
#ifdef SIGXCPU
        return SIGXCPU;
#else
        return -1;
#endif
      case CCutil_SIGXFSZ:
#ifdef SIGXFSZ
        return SIGXFSZ;
#else
        return -1;
#endif
      case CCutil_SIGSTKFLT:
#ifdef SIGSTKFLT
        return SIGSTKFLT;
#else
        return -1;
#endif
      case CCutil_SIGIOT:
#ifdef SIGIOT
        return SIGIOT;
#else
        return -1;
#endif
      case CCutil_SIGPOLL:
#ifdef SIGPOLL
        return SIGPOLL;
#else
        return -1;
#endif
      case CCutil_SIGMSG:
#ifdef SIGMSG
        return SIGMSG;
#else
        return -1;
#endif
      case CCutil_SIGDANGER:
#ifdef SIGDANGER
        return SIGDANGER;
#else
        return -1;
#endif
      case CCutil_SIGMIGRATE:
#ifdef SIGMIGRATE
        return SIGMIGRATE;
#else
        return -1;
#endif
      case CCutil_SIGPRE:
#ifdef SIGPRE
        return SIGPRE;
#else
        return -1;
#endif
      case CCutil_SIGVIRT:
#ifdef SIGVIRT
        return SIGVIRT;
#else
        return -1;
#endif
      default:
        fprintf (stderr, "Invalid signal number %d in ccsig_to_sig\n",
                 ccsignum);
        return -1;
    }
}

int CCutil_sig_to_ccsig (int signum)
{
    int i;

    for (i=1; i<CCutil_MAXSIG; i++) {
        if (ccsig_to_sig (i) == signum) {
            return i;
        }
    }
    return -1;
}

static const char *ccsignal_name (int ccsignum)
{
    switch (ccsignum) {
      case CCutil_SIGHUP:
        return "SIGHUP";
      case CCutil_SIGINT:
        return "SIGINT";
      case CCutil_SIGQUIT:
        return "SIGQUIT";
      case CCutil_SIGILL:
        return "SIGILL";
      case CCutil_SIGTRAP:
        return "SIGTRAP";
      case CCutil_SIGABRT:
        return "SIGABRT";
      case CCutil_SIGEMT:
        return "SIGEMT";
      case CCutil_SIGFPE:
        return "SIGFPE";
      case CCutil_SIGKILL:
        return "SIGKILL";
      case CCutil_SIGBUS:
        return "SIGBUS";
      case CCutil_SIGSEGV:
        return "SIGSEGV";
      case CCutil_SIGSYS:
        return "SIGSYS";
      case CCutil_SIGPIPE:
        return "SIGPIPE";
      case CCutil_SIGALRM:
        return "SIGALRM";
      case CCutil_SIGTERM:
        return "SIGTERM";
      case CCutil_SIGUSR1:
        return "SIGUSR1";
      case CCutil_SIGUSR2:
        return "SIGUSR2";
      case CCutil_SIGCHLD:
        return "SIGCHLD";
      case CCutil_SIGPWR:
        return "SIGPWR";
      case CCutil_SIGWINCH:
        return "SIGWINCH";
      case CCutil_SIGURG:
        return "SIGURG";
      case CCutil_SIGIO:
        return "SIGIO";
      case CCutil_SIGSTOP:
        return "SIGSTOP";
      case CCutil_SIGTSTP:
        return "SIGTSTP";
      case CCutil_SIGCONT:
        return "SIGCONT";
      case CCutil_SIGTTIN:
        return "SIGTTIN";
      case CCutil_SIGTTOU:
        return "SIGTTOU";
      case CCutil_SIGVTALRM:
        return "SIGVTALRM";
      case CCutil_SIGPROF:
        return "SIGPROF";
      case CCutil_SIGXCPU:
        return "SIGXCPU";
      case CCutil_SIGXFSZ:
        return "SIGXFSZ";
      case CCutil_SIGSTKFLT:
        return "SIGSTKFLT";
      case CCutil_SIGIOT:
        return "SIGIOT";
      case CCutil_SIGPOLL:
        return "SIGPOLL";
      case CCutil_SIGMSG:
        return "SIGMSG";
      case CCutil_SIGDANGER:
        return "SIGDANGER";
      case CCutil_SIGMIGRATE:
        return "SIGMIGRATE";
      case CCutil_SIGPRE:
        return "SIGPRE";
      case CCutil_SIGVIRT:
        return "SIGVIRT";
      default:
        return "UNKNOWN SIGNAL";
    }
}

void CCutil_handler_fatal (int signum)
{
    int ccsignum = CCutil_sig_to_ccsig (signum);
    int i;

#ifdef CCSIGNAL_SYSV
    CCutil_signal_handler (ccsignum, CCutil_handler_fatal);
#endif
    fprintf (stderr, "FATAL ERROR - received signal %s (%d/%d)\n",
             ccsignal_name (ccsignum), ccsignum, signum);
    fflush (stdout);
#ifdef HAVE_SLEEP
    for (i=1; i>0; i--) {
        fprintf (stderr, "sleeping %d more hours to permit debugger access\n",
                 i);
        sleep (3600);
    }
#endif /* HAVE_SLEEP */
    fprintf (stderr, "FATAL ERROR - exiting\n");
    exit (-1);
}

void CCutil_handler_warn (int signum)
{
    int ccsignum = CCutil_sig_to_ccsig (signum);
    
#ifdef CCSIGNAL_SIGNAL
    CCutil_signal_handler (ccsignum, CCutil_handler_warn);
#endif
    fprintf (stderr, "WARNING - received signal %s (%d/%d)\n",
             ccsignal_name (ccsignum), ccsignum, signum);
}

void CCutil_handler_exit (int signum)
{
    int ccsignum = CCutil_sig_to_ccsig (signum);
    
    fprintf (stderr, "EXITING - received signal %s (%d/%d)\n",
             ccsignal_name (ccsignum), ccsignum, signum);
    exit (1);
}

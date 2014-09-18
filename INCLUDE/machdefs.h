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

#ifndef __MACHDEFS_H
#define __MACHDEFS_H

#define NDEBUG

#include "config.h"

#ifdef CC_POSIXTHREADS
#ifdef CC_SIGNAL_BEFORE_PTHREAD
#include <signal.h>
#endif
#include <pthread.h>
#endif

#include <stdio.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif
#ifdef HAVE_MATH_H
# include <math.h>
#endif
#ifdef HAVE_STRING_H
# include <string.h>
#else
# ifdef HAVE_STRINGS_H
#  include <strings.h>
# endif
#endif
#ifdef HAVE_ERRNO_H
# include <errno.h>
#endif
#ifdef HAVE_ASSERT_H
# include <assert.h>
#endif
#ifdef TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# ifdef HAVE_SYS_TIME_H
#  include <sys/time.h>
# else
#  ifdef HAVE_TIME_H
#   include <time.h>
#  endif
# endif
#endif
#ifdef HAVE_STDDEF_H
# include <stddef.h>
#endif
#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif
#ifdef HAVE_MALLOC_H
# include <malloc.h>
#endif
#ifdef HAVE_SYS_TYPES_H
# include <sys/types.h>
#endif
#ifdef HAVE_SYS_STAT_H
# include <sys/stat.h>
#endif
#ifdef HAVE_SYS_RESOURCE_H
# include <sys/resource.h>
#endif
#ifdef HAVE_FCNTL_H
# include <fcntl.h>
#endif
#ifdef HAVE_SIGNAL_H
# include <signal.h>
#endif
#ifdef HAVE_SYS_SOCKET_H
# include <sys/socket.h>
#endif
#ifdef HAVE_NETDB_H
# include <netdb.h>
#endif
#ifdef HAVE_NETINET_IN_H
# include <netinet/in.h>
#endif

#ifdef HAVE_SOCKET
#define CC_NETREADY
#endif

#ifdef CC_ATTRIBUTE
#define CC_UNUSED __attribute__ ((unused))
#else
#define CC_UNUSED
#endif

#ifdef CC_PROTO_PRINTF
/* assume that if you're missing printf, you're missing a bunch */
extern int
    printf (const char *, ...),
    fprintf (FILE *, const char *, ...),
    fflush (FILE *),
    scanf (const char *, ...),
    sscanf (const char *, const char *, ...),
    fscanf (FILE *, const char *, ...),
    fclose (FILE *),
    ungetc (int, FILE *),
    _filbuf (FILE *),
    time (int *);
#ifdef CC_NETREADY
extern int
    socket (int, int, int),
    connect (int, const struct sockaddr *, int),
    accept (int, struct sockaddr *, int *),
    bind (int, const struct sockaddr *, int),
    listen (int, int);
#endif
extern void
   *memset (void *, int, size_t),
    perror (const char *);
#endif

#ifdef CC_PROTO_RENAME
extern int
    rename (const char *, const char *);
#endif

#ifdef CC_PROTO_GETHOSTNAME
extern int
    gethostname (char *, int);
#endif

#ifdef CC_PROTO_GETRUSAGE
extern int
    getrusage (int, struct rusage *);
#endif

#ifdef CC_PROTO___VFORK
extern pid_t
    __vfork (void);
#endif

#ifndef NULL
#define NULL (0)
#endif

#ifndef INT_MAX
#define INT_MAX ((int) (~(((unsigned) 1) << ((8*sizeof(int))-1))))
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef SEEK_SET
#ifdef L_SET
#define SEEK_SET L_SET
#else
#define SEEK_SET 0
#endif
#endif

#ifdef CC_BADSIGDEF_CAST

#undef SIG_ERR
#undef SIG_DFL
#undef SIG_IGN
#define SIG_ERR ((void(*)(int))-1)
#define SIG_DFL ((void(*)(int))0)
#define SIG_IGN ((void(*)(int))1)

#endif

#endif  /* __MACHDEFS_H */

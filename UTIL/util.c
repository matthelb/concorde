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
/*               MISCELLANEOUS UTILITY ROUTINES                             */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: October 12, 1995                                                  */
/*  Date: September 28, 1997                                                */
/*        April 7, 2003 (bico)                                              */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  unsigned int CCutil_nextprime (unsigned int x)                          */
/*    FINDS the smallest positive prime >= x                                */
/*                                                                          */
/*  int CCutil_our_gcd (int a, int b)                                       */
/*    COMPUTES gcd(a,b)                                                     */
/*    -gcd(a,b) is always >= 0                                              */
/*    -a and b can be negative, positive, or zero                           */
/*                                                                          */
/*  int CCutil_our_lcm (int a, int b)                                       */
/*    COMPUTES lcm(a,b)                                                     */
/*    -lcm(a,b) is always >= 0                                              */
/*    -a and b can be negative, positive, or zero                           */
/*                                                                          */
/*  char *CCutil_strchr (char *s, int c)                                    */
/*    RETURNS a pointer to the first occurrence of c in s, or NULL if c     */
/*    does not occur in s                                                   */
/*                                                                          */
/*  char *CCutil_strrchr (char *s, int c)                                   */
/*    RETURNS a pointer to the last occurrence of c in s, or NULL if c      */
/*    does not occur in s                                                   */
/*                                                                          */
/*  const char *CCutil_strchr_c (const char *s, int c)                      */
/*    RETURNS a pointer to the first occurrence of c in s, or NULL if c     */
/*    does not occur in s.  A variant for const strings.                    */
/*                                                                          */
/*  const char *CCutil_strrchr_c (const char *s, int c)                     */
/*    RETURNS a pointer to the last occurrence of c in s, or NULL if c      */
/*    does not occur in s.  A variant for const strings.                    */
/*                                                                          */
/*  char *CCutil_strdup (const char *s)                                     */
/*    RETURNS a pointer to a copy of s, allocated with CC_SAFE_MALLOC,      */
/*    or NULL if unable to allocate space for the string                    */
/*                                                                          */
/*  char *CCutil_strdup2 (const char *s)                                    */
/*    RETURNS a pointer to a copy of s up until the first whitespace,       */
/*    allocated with CC_SAFE_MALLOC, or NULL if unable to allocate space    */
/*    for the string.                                                       */
/*                                                                          */
/*  void CCutil_readstr (FILE *f, char *s, int len)                         */
/*    READS a string from f into s.  The string is terminated by a          */
/*    whitespace character (space, tab, newline) or EOF.  The entire        */
/*    string including the terminating character is read, but only at       */
/*    most len characters are stored in s, including the terminating        */
/*    NULL.                                                                 */
/*                                                                          */
/*  void CCutil_printlabel (void)                                           */
/*    PRINTS information identifying a machine                              */
/*                                                                          */
/*  int CCutil_print_command (int ac, char **av)                            */
/*    PRINTS the command line                                               */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "macrorus.h"
#include "util.h"


static int
    isprime (unsigned int x);


unsigned int CCutil_nextprime (unsigned int x)
{
    if (x < 3) return 3;
    x |= 1;
    while (!isprime (x)) x += 2;
    return x;
}

static int isprime (unsigned int p)
{
    unsigned int i;

    if ((p&1) == 0) return 0;
    for (i=3; i*i<=p; i+=2) {
        if (p%i == 0) return 0;
    }
    return 1;
}

int CCutil_our_gcd (int a, int b)
{
    int c;

    if (a < 0) a = -a;
    if (b < 0) b = -b;
    if (a > b) CC_SWAP (a, b, c);

    while (a) {
        c = b % a;
        b = a;
        a = c;
    }
    return b;
}

int CCutil_our_lcm (int a, int b)
{
    int c;

    if (a < 0) a = -a;
    if (b < 0) b = -b;

    c = CCutil_our_gcd (a, b);

    return (a/c) * b;
}

char *CCutil_strchr (char *s, int c)
{
    while (*s) {
        if (*s == c) return s;
        s++;
    }
    return (char *) NULL;
}

const char *CCutil_strchr_c (const char *s, int c)
{
    while (*s) {
        if (*s == c) return s;
        s++;
    }
    return (char *) NULL;
}

char *CCutil_strrchr (char *s, int c)
{
    char *l = (char *) NULL;

    while (*s) {
        if (*s == c) l = s;
        s++;
    }
    return l;
}

const char *CCutil_strrchr_c (const char *s, int c)
{
    const char *l = (char *) NULL;

    while (*s) {
        if (*s == c) l = s;
        s++;
    }
    return l;
}

char *CCutil_strdup (const char *s)
{
    char *p = CC_SAFE_MALLOC (strlen(s)+1, char);

    if (p == (char *) NULL) {
        fprintf (stderr, "Out of memory in CCutil_strdup\n");
        return (char *) NULL;
    }
    strcpy (p, s);
    return p;
}

char *CCutil_strdup2 (const char *s)
{
    char *p;
    int len;

    for (len = 0; s[len]; len++) {
        if (s[len] == ' ' || s[len] == '\r' || s[len] == '\t' ||
            s[len] == '\n') {
            break;
        }
    }

    p = CC_SAFE_MALLOC (len+1, char);

    if (p == (char *) NULL) {
        fprintf (stderr, "Out of memory in CCutil_strdup2\n");
        return (char *) NULL;
    }
    strncpy (p, s, len);
    p[len] = '\0';

    return p;
}

void CCutil_readstr (FILE *f, char *s, int len)
{
    int c;

    while ((c = getc (f)) != EOF && c != ' ' && c != '\t' && c != '\r' &&
           c != '\n') {
        if (len > 1) {
            *s++ = c;
            len--;
        }
    }
    *s = '\0';
}

void CCutil_printlabel (void)
{
#ifdef CC_NETREADY
    char buf[1024];

    gethostname (buf, 1024);
    printf ("Host: %s  Current process id: %d\n", buf, (int) getpid());
    fflush (stdout);
#else
    printf ("No label - need function to non-NETREADY machines\n");
    fflush (stdout);
#endif
}

int CCutil_print_command (int ac, char **av)
{
    int rval = 0;
    int i, cmdlen = 0;
    char *cmdout = (char *) NULL;

    for (i=0; i<ac; i++) {
        cmdlen += strlen(av[i]) + 1;
    }
    cmdout = CC_SAFE_MALLOC (cmdlen, char);
    CCcheck_NULL (cmdout, "out of memory in print_command");

    cmdlen = 0;
    for (i=0; i<ac; i++) {
        strcpy (cmdout + cmdlen, av[i]);
        cmdlen += strlen(av[i]);
        cmdout[cmdlen] = ' ';
        cmdlen++;
    }
    cmdout[cmdlen-1] = '\0';
    printf ("%s\n", cmdout); fflush (stdout);

CLEANUP:

    CC_IFFREE (cmdout, char);
    return rval;
}


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
/*              MACHINE INDEPENDENT RANDOM NUMBER GENERATOR                 */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  DIMACS  (modified for TSP)                                 */
/*  Date: February 7, 1995  (cofeb16)                                       */
/*        September 18, 2001  (billenium fix)                               */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  void CCutil_sprand (int seed, CCrandstate *r)                           */
/*    - Call once to initialize the generator.                              */
/*                                                                          */
/*  int CCutil_lprand (CCrandstate *r)                                      */
/*    - Returns an integer in the range 0 to CC_PRANDMAX - 1.               */
/*                                                                          */
/*  double CCutil_normrand (CCrandstate *r)                                 */
/*    - Returns a normally-distributed random value with mean 0 and         */
/*      deviation 1.                                                        */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*    NOTES (from DIMACS):                                                  */
/*        This file contains a set of c-language functions for generating   */
/*    uniform integers.   This is a COMPLETELY PORTABLE generator. It will  */
/*    give IDENTICAL sequences of random numbers for any architecture with  */
/*    at least 30-bit integers, regardless of the integer representation,   */
/*    INT_MAX value, or roundoff/truncation method, etc.                    */
/*        This Truly Remarkable RNG is described more fully in              */
/*    J. Bentley's column, ``The Software Exploratorium ''. It is based on  */
/*    one in Knuth, Vol 2, Section 3.2.2 (Algorithm A).                     */
/*                                                                          */
/*  CCutil_normrand is not from DIMACS or Bentley, but rather just uses     */
/*  the Box-Muller transformation to generate a normally-distributed        */
/*  random variable from two uniform ones.                                  */
/*                                                                          */
/****************************************************************************/


#include "machdefs.h"
#include "util.h"


void CCutil_sprand (int seed, CCrandstate *r)
{
    int i, ii;
    int last, next;
    int *arr = r->arr;

    seed %= CC_PRANDMAX;
    if (seed < 0) seed += CC_PRANDMAX;

    arr[0] = last = seed;
    next = 1;
    for (i = 1; i < 55; i++) {
        ii = (21 * i) % 55;
        arr[ii] = next;
        next = last - next;
        if (next < 0)
            next += CC_PRANDMAX;
        last = arr[ii];
    }
    r->a = 0;
    r->b = 24;
    for (i = 0; i < 165; i++)
        last = CCutil_lprand (r);
}


int CCutil_lprand (CCrandstate *r)
{
    int t;

    if (r->a-- == 0)
        r->a = 54;
    if (r->b-- == 0)
        r->b = 54;

    t = r->arr[r->a] - r->arr[r->b];

    if (t < 0)
        t += CC_PRANDMAX;

    r->arr[r->a] = t;

    return t;
}


#ifdef      TRY_CODE

/*-----------------------------------------------*/
/* This is a little driver program so you can    */
/* test the code.                                */
/* Typing: a.out 0 3 1                           */
/* should produce                                */
/*     921674862                                 */
/*     250065336                                 */
/*     377506581                                 */
/*  Typing: a.out 1000000 1 2                    */
/*  should produce                               */
/*     57265995                                  */
/*-----------------------------------------------*/

int main (int ac, char **av)
{
    int i;
    int j;
    int n;
    int m;
    int seed;
    CCrandstate rstate;

    if (ac < 4) {
        fprintf (stderr, "Usage: #discard #print #seed\n");
        return 0;
    }
    m = atoi (av[1]);           /* Number to discard initially */
    n = atoi (av[2]);           /* Number to print */
    seed = atoi (av[3]);        /* Seed */

    CCutil_sprand (seed, &rstate);

    for (i = 0; i < m; i++)
        j = CCutil_lprand (&rstate);
    for (i = 0; i < n; i++)
        printf ("%ld\n", CCutil_lprand (&rstate));
    return 0;
}

#endif  /* TRY_CODE */


double CCutil_normrand (CCrandstate *r)
{
    double x1 = ((double) CCutil_lprand(r)) / ((double) CC_PRANDMAX);
    double x2 = ((double) CCutil_lprand(r)) / ((double) CC_PRANDMAX);

    return sqrt (-2*log(x1)) * cos(2*M_PI*x2);

}

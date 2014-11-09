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
/*                LINEAR SUBTOUR SEPARATION ROUTINES                        */
/*                           preliminary                                    */
/*                                                                          */
/*                              TSP CODE                                    */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: July 25, 1995                                                     */
/*  Thanks: Phil Gibbons, for helpful ideas                                 */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCcut_linsub (int ncount, int ecount, int *endmark, int *elist,     */
/*      double *x, double maxval, void *u_data,                             */
/*      int (*cut_callback) (double cut_val, int cut_start, int cut_end,    */
/*      void *u_data))                                                      */
/*    -ncount is the number of nodes                                        */
/*    -ecount is the number of edges                                        */
/*    -endmark indicates which segments are of interest, by indicating      */
/*     whether a node can be a right or left end of a segment               */
/*    -elist contains the LP edges in node node format                      */
/*    -x is an LP solution                                                  */
/*    -maxval is the maximum cut value desired                              */
/*    -u_data is user data to be passed to cut_callback                     */
/*    -cut_callback is a function to be called for segments which           */
/*     define a cut of value < cutlim.  The cut is cut_start,               */
/*     cut_start+1, ..., cut_end, and has value cut_val.  cut_callback      */
/*     will be called for the minimum segment cut starting at each          */
/*     endpoint marked as a right end, provided that cut has value <        */
/*     cutlim.                                                              */
/*                                                                          */
/*  int CCcut_linsub_allcuts (int ncount, int ecount, int *perm,            */
/*      int *endmark, int *elist, double *x, double maxval,                 */
/*      void *u_data, int (*cut_callback) (double cut_val,                  */
/*      int cut_start, int cut_end, void *u_data))                          */
/*        -ncount is the number of nodes                                    */
/*        -ecount is the number of edges                                    */
/*        -perm is a permutation of the nodes (if perm == (int *) NULL,     */
/*         the identity permutation will be used)                           */
/*        -elist contains the LP edges in node node format                  */
/*        -endmark indicates which segments are of interest, by indicating  */
/*         whether a node can be a right or left end of a segment           */
/*        -x is an LP solution                                              */
/*        -maxval is the maximum cut value desired                          */
/*        -u_data is data to be passed to the callback                      */
/*        -cut_callback is a function to be called for every segment which  */
/*         defines a cut of value <= cutlim.  The cut is perm[cut_start],   */
/*         perm[cut_start+1], ..., perm[cut_end], and has value cut_val.    */
/*         if cut_callback returns a nonzero value, CCcut_linsub_allcuts    */
/*         will terminate.                                                  */
/*                                                                          */
/*    NOTES:                                                                */
/*        CCcut_linsub runs in time O(m log n).                             */
/*        CCcut_linsub_allcuts runs in time O(m log n + |C| log n) where    */
/*        |C| is the number of cuts found.                                  */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "cut.h"

#define LINSUB_INF 1e20

typedef struct psh {
    int base;
    int size;
    double *sum;
    double *minpre;
} psh;


static void
    psh_free (psh *p),
    psh_add (psh *p, int i, double v);

static int
    psh_init (psh *p, int k, int *endmark),
    psh_minloc (psh *p),
    psh_enum (psh *p, int cut_start, double maxval, void *u_data,
        int (*cut_callback)(double cut_val, int cut_start, int cut_end,
        void *u_data)),
    psh_enum_work (psh *p, int n, int mul, double pre_sum, int cut_start,
        double maxslack, void *u_data, int (*cut_callback)(double cut_val,
        int cut_start, int cut_end, void *u_data));

static double
    psh_minval (psh *p);


int CCcut_linsub (int ncount, int ecount, int *endmark, int *elist, double *x,
        double maxval, void *u_data, int (*cut_callback) (double cut_val,
        int cut_start, int cut_end, void *u_data))
{
    psh p;
    int i, j;
    double v;
    int rval = 0;
    int *perm = (int *) NULL;
    int *eperm = (int *) NULL;
    int *ends = (int *) NULL;
    double *xends = (double *) NULL;

    if (psh_init (&p, ncount, endmark)) {
        return -1;
    }

    /*  arrange elist into the array ends (with xends being the x values   */
    /*  so that ends[2*i] is less than ends[2*i+1], and ends is sorted in  */
    /*  increasing  order of ends[2*i]                                     */

    perm = CC_SAFE_MALLOC (ecount, int);
    eperm = CC_SAFE_MALLOC (ecount, int);
    if (!perm || !eperm) {
        fprintf (stderr, "out of memory in CCcut_linsub\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ecount; i++) {
        eperm[i] = (elist[2*i] < elist[2*i+1] ? elist[2*i] : elist[2*i+1]);
        perm[i] = i;
    }
    CCutil_int_perm_quicksort (perm, eperm, ecount);

    ends = CC_SAFE_MALLOC (2*ecount, int);
    xends = CC_SAFE_MALLOC (ecount, double);
    if (!ends || !xends) {
        fprintf (stderr, "out of memory in CCcut_linsub\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ecount; i++) {
        j = perm[i];
        if (elist[2*j] < elist[2*j+1]) {
            ends[2*i] = elist[2*j];
            ends[2*i+1] = elist[2*j+1];
        } else {
            ends[2*i] = elist[2*j+1];
            ends[2*i+1] = elist[2*j];
        }
        xends[i] = x[j];
    }
    CC_FREE (perm, int);
    CC_FREE (eperm, int);

    for (i = ncount - 1, j = ecount - 1; i > 0; i--) {
        while (j >= 0 && ends[2*j] == i) {
            psh_add (&p, ends[2*j+1], -xends[j]);

            j--;
        }
        if (endmark[i] & CC_LINSUB_LEFT_END) {
            v = 2.0 + 2.0 * psh_minval (&p);
            if (v < maxval) {
                rval = (*cut_callback) (v, i, psh_minloc (&p), u_data);
                if (rval) {
                    fprintf (stderr, "cut_callback failed\n"); goto CLEANUP;
                }
            }
        }
        psh_add (&p, i, 1.0);
    }

CLEANUP:

    psh_free (&p);
    CC_IFFREE (ends, int);
    CC_IFFREE (xends, double);
    CC_IFFREE (perm, int);
    CC_IFFREE (eperm, int);

    return rval;
}

int CCcut_linsub_allcuts (int ncount, int ecount, int *perm, int *endmark,
        int *elist, double *x, double maxval, void *u_data,
        int (*cut_callback) (double cut_val, int cut_start, int cut_end,
        void *u_data))
{
    psh p;
    int i, j;
    int rval = 0;
    int *perm_inv = (int *) NULL;
    int *esort = (int *) NULL;
    int *eperm = (int *) NULL;
    int *ends = (int *) NULL;
    int *pendmark = (int *) NULL;
    double *xends = (double *) NULL;

    perm_inv = CC_SAFE_MALLOC (ncount, int);
    pendmark = CC_SAFE_MALLOC (ncount, int);
    if (perm_inv == (int *) NULL ||
        pendmark == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCcut_linsub_allcuts\n");
        CC_IFFREE (perm_inv, int);
        CC_IFFREE (pendmark, int);
        return -1;
    }

    if (perm == (int *) NULL) {
        for (i=0; i<ncount; i++) {
            perm_inv[i] = i;
            pendmark[i] = endmark[i];
        }
    } else {
        for (i=0; i<ncount; i++) {
            perm_inv[perm[i]] = i;
            pendmark[i] = endmark[perm[i]];
        }
    }

    if (psh_init (&p, ncount, pendmark)) {
        fprintf (stderr, "psh_init failed\n");
        CC_IFFREE (perm_inv, int);
        CC_IFFREE (pendmark, int);
        return -1;
    }

    /*  arrange elist into the array ends (with xends being the x values   */
    /*  so that ends[2*i] is less than ends[2*i+1], and ends is sorted in  */
    /*  increasing  order of ends[2*i]                                     */

    eperm = CC_SAFE_MALLOC (ecount, int);
    esort = CC_SAFE_MALLOC (ecount, int);
    if (eperm == (int *) NULL || esort == (int *) NULL) {
        fprintf (stderr, "out of memory in CCcut_linsub_allcuts\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ecount; i++) {
        esort[i] = (perm_inv[elist[2*i]] < perm_inv[elist[2*i+1]]
                    ? perm_inv[elist[2*i]]
                    : perm_inv[elist[2*i+1]]);
        eperm[i] = i;
    }
    CCutil_int_perm_quicksort (eperm, esort, ecount);

    ends = CC_SAFE_MALLOC (2*ecount, int);
    xends = CC_SAFE_MALLOC (ecount, double);
    if (ends == (int *) NULL || xends == (double *) NULL) {
        fprintf (stderr, "out of memory in CCcut_linsub_allcuts\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ecount; i++) {
        j = eperm[i];
        if (perm_inv[elist[2*j]] < perm_inv[elist[2*j+1]]) {
            ends[2*i] = perm_inv[elist[2*j]];
            ends[2*i+1] = perm_inv[elist[2*j+1]];
        } else {
            ends[2*i] = perm_inv[elist[2*j+1]];
            ends[2*i+1] = perm_inv[elist[2*j]];
        }
        xends[i] = x[j];
    }
    CC_FREE (eperm, int);
    CC_FREE (esort, int);

    for (i = ncount - 1, j = ecount - 1; i > 0; i--) {
        while (j >= 0 && ends[2*j] == i) {
            psh_add (&p, ends[2*j+1], -xends[j]);
            j--;
        }
        if (pendmark[i] & CC_LINSUB_LEFT_END) {
            rval = psh_enum (&p, i, maxval, u_data, cut_callback);
            if (rval) {
                fprintf (stderr, "psh_enum failed\n");
                goto CLEANUP;
            }
        }
        psh_add (&p, i, 1.0);
    }

    rval = 0;

CLEANUP:

    psh_free (&p);
    CC_IFFREE (ends, int);
    CC_IFFREE (xends, double);
    CC_IFFREE (eperm, int);
    CC_IFFREE (esort, int);
    CC_IFFREE (pendmark, int);
    CC_IFFREE (perm_inv, int);

    return rval;
}

/****************************************************************************/
/*                                                                          */
/*                      PREFIX SUM HEAP ROUTINES                            */
/*                                                                          */
/*                              TSP CODE                                    */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: July 25, 1995                                                     */
/*  Thanks: Phil Gibbons, for the algorithm                                 */
/*                                                                          */
/*  These routines implement a heap of prefix sums.  The routines are       */
/*  self-contained, but unlikely to be useful outside of the linear         */
/*  subtour separation code.                                                */
/*                                                                          */
/*  PREFIX SUM HEAP FUNCTIONS:                                              */
/*    int psh_init (psh *p, int k, int *endmark)                            */
/*      -p should point to a psh struct.                                    */
/*      -k is the number of elements in the psh.                            */
/*      -endmark marks which nodes can have prefix sums of interest         */
/*      -this initializes the value of all k elements to 0.0.               */
/*      -return value 0 means success, 1 means failure.                     */
/*    void psh_free (psh *p)                                                */
/*      -frees the spaces allocated by psh_init.                            */
/*    void psh_add (psh *p, int i, double v)                                */
/*      -adds v to the value of element i.                                  */
/*    double psh_minval (psh *p)                                            */
/*      -returns min_j sum_{i=0}^j value[i].                                */
/*       over j with endmark[j] & CC_LINSUB_RIGHT_END                       */
/*    int psh_minloc (psh *p)                                               */
/*      -returns the smallest j which achieves psh_minval(p).               */
/*       over j with endmark[j] & CC_LINSUB_RIGHT_END                       */
/*    int psh_enum (psh *p, int cut_start, double maxval, void *u_data,     */
/*            int (*cut_callback)(double cut_val, int cut_start,            */
/*            int cut_end, void *u_data))                                   */
/*      -calls cut_callback for each interval starting at cut_start         */
/*       defining a cut <= maxval                                           */
/*       with endmark[cut_start] & CC_LINSUB_LEFT_END, and                  */
/*       with endmark[cut_end] & CC_LINSUB_RIGHT_END.                       */
/*                                                                          */
/*  NOTES:                                                                  */
/*      A k-element heap will malloc 32k bytes of memory. If memory is      */
/*  tight, using integer values (instead of doubles) brings it down to      */
/*  16k bytes.                                                              */
/*      psh_init takes O(k) time.  psh_add and psh_minloc take O(log k)     */
/*  time.  psh_free and psh_minval take O(1) time.                          */
/*                                                                          */
/*      It is likely that using a ternary tree instead of binary would      */
/*  improve performance.  Also, psh_add could take advantage of knowing     */
/*  which child has changed.                                                */
/*      psh_minloc could be changed to return all elements which achieve    */
/*  the min, in time O(log k) per element returned.                         */
/*                                                                          */
/****************************************************************************/

static int psh_init (psh *p, int k, int *endmark)
{
    int i;
    int space;

    p->size = k;
    p->base = 1;
    while (p->base < k) p->base *= 2;

    space = p->base*2;

    p->sum = CC_SAFE_MALLOC (space, double);
    if (!p->sum)
        return 1;
    p->minpre = CC_SAFE_MALLOC (space, double);
    if (!p->minpre) {
        CC_FREE (p->sum, double);
        return 1;
    }

    for (i=0; i<k; i++) {
        p->sum[p->base + i] = 0.0;
        if (endmark[i] & CC_LINSUB_RIGHT_END) {
            p->minpre[p->base + i] = 0.0;
        } else {
            p->minpre[p->base + i] = LINSUB_INF;
        }
    }
    for (i=k; i<p->base; i++) {
        p->sum[p->base + i] = 0.0;
        p->minpre[p->base + i] = LINSUB_INF;
    }
    for (i=p->base - 1; i>=1; i--) {
        p->sum[i] = p->sum[2*i] + p->sum[2*i+1];
        if (p->minpre[2*i] < p->sum[2*i] + p->minpre[2*i+1]) {
            p->minpre[i] = p->minpre[2*i];
        } else {
            p->minpre[i] = p->sum[2*i] + p->minpre[2*i+1];
        }
    }

    return 0;
}

static void psh_free (psh *p)
{
    CC_FREE (p->minpre, double);
    CC_FREE (p->sum, double);
    p->size = 0;
    p->base = 0;
}

static void psh_add (psh *p, int i, double v)
{
    double *s = p->sum;
    double *m = p->minpre;

    i += p->base;
    s[i] += v;
    m[i] += v;
    i /= 2;
    while (i >= 1) {
        s[i] += v;
        if (m[2*i] < s[2*i] + m[2*i+1]) {
            m[i] = m[2*i];
        } else {
            m[i] = s[2*i] + m[2*i+1];
        }
        i /= 2;
    }
}

static int psh_minloc (psh *p)
{
    int i = 1;
    double *m = p->minpre;

    while (i < p->base) {
        if (m[i] == m[2*i]) {
            i = 2*i;
        } else {
            i = 2*i+1;
        }
    }
    if (i - p->base < p->size) {
        return i - p->base;
    } else {
        /* we're lost */
        return p->size - 1;
    }
}

static double psh_minval (psh *p)
{
    return p->minpre[1];
}

static int psh_enum (psh *p, int cut_start, double maxval, void *u_data,
        int (*cut_callback)(double cut_val, int cut_start, int cut_end,
        void *u_data))
{
    int rval;

    rval = psh_enum_work (p, 1, p->base, 0.0, cut_start,
                          (maxval - 2.0) / 2.0, u_data, cut_callback);
    if (rval) {
        fprintf (stderr, "psh_enum_work failed\n");
        return rval;
    }
    return 0;
}

static int psh_enum_work (psh *p, int n, int mul, double pre_sum,
        int cut_start, double maxslack, void *u_data,
        int (*cut_callback)(double cut_val, int cut_start, int cut_end,
        void *u_data))
{
    int rval;

    if (n >= p->base) {
        rval = (*cut_callback) (2.0 + 2.0*(pre_sum + p->sum[n]),
                cut_start, n - p->base, u_data);
        if (rval) {
            fprintf (stderr, "cut_callback failed\n");
            return rval;
        }
        return 0;
    }

    mul /= 2;
    n *= 2;
    if (n*mul < p->base + p->size &&
        n*mul + mul > p->base + cut_start &&
        pre_sum + p->minpre[n] <= maxslack) {
        rval = psh_enum_work (p, n, mul, pre_sum, cut_start,
                              maxslack, u_data, cut_callback);
        if (rval) return rval;
    }
    if (n*mul + mul < p->base + p->size &&
        pre_sum + p->sum[n] + p->minpre[n+1] <= maxslack) {
        rval = psh_enum_work (p, n+1, mul, pre_sum + p->sum[n], cut_start,
                              maxslack, u_data, cut_callback);
        if (rval) return rval;
    }
    return 0;
}


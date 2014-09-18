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
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCchunk_intmat_build (CCchunk_intmat *mat_p, int ncols)             */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCchunk_intmat_addrow (CCchunk_intmat *mat_p, int *row)             */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCchunk_intmat_ortho (CCchunk_intmat *mat_p, int *ortho,            */
/*      int *pcol_p, int *taboo)                                            */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCchunk_intmat_init (CCchunk_intmat *mat_p)                        */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCchunk_intmat_free (CCchunk_intmat *mat_p)                        */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCchunk_intmat_dellastrows (CCchunk_intmat *mat_p, int ndel)       */
/*    MISSING                                                               */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "localcut.h"
#include "macrorus.h"

/* Macros */

#define MATRIX(i,j)  (matrix[(n*(i))+(j)])
#define FACTOR(i,j)  (factor[(n*(i))+(j)])
#define ROW(i)       (factor+(n*(i)))
#define MULTLIM(a)   (CC_MATVAL_MAX/CC_OURABS(a))

/* Number of rows added to rowspace in reallocation */

#define ROWGROWTH  10

#define SOLVE_TEST
#undef  DUMP_BASIS

/* Static functions */


static int
    intmat_realloc (CCchunk_intmat *mat_p),
    intmat_solve (int m, int n, int *rperm, int *cperm, CCmatval *factor,
        CCmatval *x, int j_rhs),
    vector_multiple (int ncols, CCmatval *v1, int *v2);

static CCmatval
    matval_gcd (CCmatval a, CCmatval b);


static void
    intmat_reducevec (int beg, int end, int *perm, CCmatval *vec),
    print_matrix (CCchunk_intmat *mat_p);



void CCchunk_intmat_init (CCchunk_intmat *mat_p)
{
    mat_p->matrix = (int *) NULL;
    mat_p->factor = (CCmatval *) NULL;
    mat_p->csize  = (int *) NULL;
    mat_p->rperm  = (int *) NULL;
    mat_p->cperm  = (int *) NULL;
    mat_p->x      = (CCmatval *) NULL;
    mat_p->best_x = (CCmatval *) NULL;

    mat_p->nrows    = 0;
    mat_p->ncols    = 0;
    mat_p->rowspace = 0;
} /* END CCchunk_intmat_init */


void CCchunk_intmat_free (CCchunk_intmat *mat_p)
{
    CC_IFFREE (mat_p->matrix, int);
    CC_IFFREE (mat_p->factor, CCmatval);
    CC_IFFREE (mat_p->csize, int);
    CC_IFFREE (mat_p->rperm, int);
    CC_IFFREE (mat_p->cperm, int);
    CC_IFFREE (mat_p->x, CCmatval);
    CC_IFFREE (mat_p->best_x, CCmatval);

    mat_p->nrows    = 0;
    mat_p->ncols    = 0;
    mat_p->rowspace = 0;
} /* END CCchunk_intmat_free */


int CCchunk_intmat_build (CCchunk_intmat *mat_p, int ncols)
{
    int  rval = 0;
    int  m, n, mn;

    m  = ncols + 1;
    n  = ncols + 1;
    mn = m * n;

    mat_p->matrix = (int *) CC_SAFE_MALLOC (mn, int);
    mat_p->factor = (CCmatval *) CC_SAFE_MALLOC (mn, CCmatval);
    mat_p->csize  = (int *) CC_SAFE_MALLOC (n, int);
    mat_p->rperm  = (int *) CC_SAFE_MALLOC (m, int);
    mat_p->cperm  = (int *) CC_SAFE_MALLOC (n, int);
    mat_p->x      = (CCmatval *) CC_SAFE_MALLOC (n, CCmatval);
    mat_p->best_x = (CCmatval *) CC_SAFE_MALLOC (n, CCmatval);

    if ( mat_p->matrix == (int *) NULL ||
         mat_p->factor == (CCmatval *) NULL ||
         mat_p->csize  == (int *) NULL ||
         mat_p->rperm  == (int *) NULL ||
         mat_p->cperm  == (int *) NULL ||
         mat_p->x      == (CCmatval *) NULL ||
         mat_p->best_x == (CCmatval *) NULL   ) {
        rval = CC_CHUNK_INTMAT_MEMORY;
        fprintf (stderr, "Out of memory in CCchunk_intmat_build\n");
        goto CLEANUP;
    }

    mat_p->ncols    = ncols;
    mat_p->nrows    = 0;
    mat_p->rowspace = m;

CLEANUP:

    if ( rval )  CCchunk_intmat_free (mat_p);
    return rval;

} /* END CCchunk_intmat_build */


static int intmat_realloc (CCchunk_intmat *mat_p)
{
    int  rval = 0;

    int  *matrix = (int *) NULL;
    CCmatval  *factor = (CCmatval *) NULL;
    int  *rperm  = (int *) NULL;
    int  m, n, mn;

    m  = mat_p->rowspace + ROWGROWTH;
    n  = mat_p->ncols + 1;
    mn = m * n;

    matrix = (int *) CCutil_reallocrus (mat_p->matrix, mn*sizeof(int));
    factor = (CCmatval *) CCutil_reallocrus (mat_p->factor, mn*sizeof(CCmatval));
    rperm  = (int *) CCutil_reallocrus (mat_p->rperm, m*sizeof(int));

    if ( matrix == (int *) NULL ||
         factor == (CCmatval *) NULL ||
         rperm  == (int *) NULL   ) {
        if ( matrix != (int *) NULL )  mat_p->matrix = matrix;
        if ( factor != (CCmatval *) NULL )  mat_p->factor = factor;
        if ( rperm  != (int *) NULL )  mat_p->rperm  = rperm;
        rval = CC_CHUNK_INTMAT_MEMORY;
        fprintf (stderr, "Out of memory in intmat_realloc\n");
        return rval;
    }
    mat_p->matrix   = matrix;
    mat_p->factor   = factor;
    mat_p->rperm    = rperm;
    mat_p->rowspace = m;
    return rval;

} /* END intmat_realloc */


int CCchunk_intmat_addrow (CCchunk_intmat *mat_p, int *row)
{
    int  rval = 0;
    int  j, n;
    int  *matrix;

    if ( mat_p->rowspace <= mat_p->nrows ) {
        rval = intmat_realloc (mat_p);
        if (rval) return rval;
    }

    matrix = mat_p->matrix;
    n      = mat_p->ncols + 1;

    for (j = 0; j < mat_p->ncols; j++) {
        MATRIX (mat_p->nrows,j) = row[j];
    }
    MATRIX (mat_p->nrows,mat_p->ncols) = -1;
    mat_p->nrows++;
    return rval;

} /* END CCchunk_intmat_addrow */


void CCchunk_intmat_dellastrows (CCchunk_intmat *mat_p, int ndel)
{
    mat_p->nrows -= ndel;
    if ( mat_p->nrows < 0 )  mat_p->nrows = 0;
} /* END CCchunk_intmat_delllastrows */


int CCchunk_intmat_ortho (CCchunk_intmat *mat_p, int *ortho, int *pcol_p,
        int *taboo)
{
    int  rval = 0;

    int  m, n, mn, i, j, prestage, stage;
    int  itemp, jtemp, ipiv, jpiv;
    CCmatval a;
    CCmatval pval;
    CCmatval eval;
    CCmatval t;
    CCmatval *prow;
    CCmatval *erow;
    CCmatval lim1;
    CCmatval lim2;
    CCmatval t1;
    CCmatval t2;
    int  best_x_size;
    int  x_size;
    int  scan_rval;

    int  ncols   = mat_p->ncols;
    int  nrows   = mat_p->nrows;
    int  *matrix = mat_p->matrix;
    CCmatval  *factor = mat_p->factor;
    int  *csize  = mat_p->csize;
    int  *rperm  = mat_p->rperm;
    int  *cperm  = mat_p->cperm;
    CCmatval *x      = mat_p->x;
    CCmatval *best_x = mat_p->best_x;

    m  = nrows;
    n  = ncols + 1;
    mn = m * n;

    /* Initialize solution */

    for (j = 0; j < ncols; j++)  ortho[j] = 0;
    *pcol_p = 0;

    /* Initialize row and column permutations */

    for (i = 0; i < m; i++)  rperm[i] = i;
    for (j = 0; j < n; j++)  cperm[j] = j;

    /* Copy matrix into space where factorization will be done */

    for (i = 0; i < mn; i++)  factor[i] = (CCmatval) matrix[i];

    /* Initialize the column counts */

    for (j=0; j<n; j++) csize[j] = 0;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            if ( FACTOR(i,j) )  csize[j]++;
        }
    }

    /* Peel off the 1-element columns */

    for (prestage = 0; prestage < m; prestage++) {
        ipiv = -1;
        jpiv = -1;

        for (j = prestage; j < n; j++) {
            if ( csize[cperm[j]] == 1 ) {
                for (i = prestage; i < m; i++) {
                    if ( FACTOR(rperm[i],cperm[j]) ) {
                        ipiv = i;
                        jpiv = j;
                        break;
                    }
                }
                break;
            }
        }

        if ( ipiv == -1  ||  jpiv == -1 )  break;

        /* Switch the indices */

        itemp           = rperm[prestage];
        jtemp           = cperm[prestage];
        rperm[prestage] = rperm[ipiv];
        cperm[prestage] = cperm[jpiv];
        rperm[ipiv]     = itemp;
        cperm[jpiv]     = jtemp;

        ipiv = rperm[prestage];
        jpiv = cperm[prestage];

        prow = ROW (ipiv);
        for (j = prestage+1; j < n; j++) {
            if ( prow[cperm[j]] )  csize[cperm[j]]--;
        }
    }

    /* Complete factorization */

    for (stage = prestage; stage < m; stage++) {
        ipiv = -1;
        jpiv = -1;
        pval = CC_MATVAL_MAX;

        /* Find a smallest nonzero element in the remaining matrix */

        for (i = stage; i < m; i++) {
            prow = ROW(rperm[i]);
            for (j = stage; j < n; j++) {
                t = CC_OURABS (prow[cperm[j]]);
                if ( t  &&  t < pval ) {
                    pval = t;
                    ipiv = i;
                    jpiv = j;
                    if ( pval == 1 )  break;
                }
            }
            if ( pval == 1 )  break;
        }

        if ( ipiv == -1  ||  jpiv == -1 )  break;

        /* Switch the indices */

        itemp        = rperm[stage];
        jtemp        = cperm[stage];
        rperm[stage] = rperm[ipiv];
        cperm[stage] = cperm[jpiv];
        rperm[ipiv]  = itemp;
        cperm[jpiv]  = jtemp;

        ipiv = rperm[stage];
        jpiv = cperm[stage];

        pval = FACTOR(ipiv,jpiv);
        prow = ROW(ipiv);
        for (i = stage+1; i < m; i++) {
            if ( (eval = FACTOR (rperm[i],jpiv)) ) {
                t = pval;
                /* Find the gcd of (eval,pval) */
                if ( CC_OURABS (eval) != 1  &&  CC_OURABS (t) != 1 ) {
                    a = matval_gcd (eval, pval);
                    t    /= a;
                    eval /= a;
                }
                erow = ROW(rperm[i]);
                erow[cperm[stage]] = 0;
                lim1 = MULTLIM (t);
                lim2 = MULTLIM (eval);
                for (j = stage+1; j < n; j++) {
                    if ( CC_OURABS (erow[cperm[j]]) > lim1 ||
                         CC_OURABS (prow[cperm[j]]) > lim2   ) {
                        rval = CC_CHUNK_INTMAT_OVERFLOW_M;
                        goto CLEANUP;
                    }
                    t1 = t    * erow[cperm[j]];
                    t2 = eval * prow[cperm[j]];

                    /* - CC_MATVAL_MAX + t2 <= t1 <= CC_MATVAL_MAX + t2 */

                    if ( (t2 > 0 && t1 < -CC_MATVAL_MAX + t2) ||
                         (t2 < 0 && t1 >  CC_MATVAL_MAX + t2)   ) {
                        rval = CC_CHUNK_INTMAT_OVERFLOW_A;
                        goto CLEANUP;
                    }
                    erow[cperm[j]] = t1 - t2;
                }
                /* Reduce erow by its gcd */
                intmat_reducevec (stage+1, n, cperm, erow);
            }
        }
    }

    if ( stage < m ) {
        fprintf (stderr, "Rows dependent\n");
        rval = CC_CHUNK_INTMAT_ERROR;
        goto CLEANUP;
    }

    if ( stage < n-1 ) {
        int  j_rhs;

        best_x_size = -1;
        scan_rval = CC_CHUNK_INTMAT_NOORTHO;

        for (j_rhs = n-1; j_rhs >= stage; j_rhs--) {
            rval = intmat_solve (stage, n, rperm, cperm, factor, x, j_rhs);
            if ( rval )  {
                printf ("col %d intmat_solve overflow, continuing\n", j_rhs);
                scan_rval = rval;
                continue;
            }

            intmat_reducevec (0, ncols+1, (int *) NULL, x);

#ifdef DUMP_BASIS
            printf ("BASIS %d:", j_rhs);
            for (j=0; j<=ncols; j++) {
                printf (" %d", x[j]);
            }
            printf ("\n");
#endif

            x_size = 0;
            for (j = 0; j < ncols; j++) {
                if (CC_OURABS(x[j]) > x_size) x_size = CC_OURABS(x[j]);
            }

            if (best_x_size != -1 && x_size > best_x_size) {
#ifdef DUMP_BASIS
                printf ("col %d size %d bestsize %d; not as good, continuing\n",
                         j_rhs, x_size, best_x_size);
#endif
                continue;
            }

            if ( ! vector_multiple ( ncols, x, taboo ) ) {
#ifdef DUMP_BASIS
                printf ("col %d size %d bestsize %d new best\n", j_rhs,
                        x_size, best_x_size);
#endif
                for (j=0; j<=ncols; j++) {
                    best_x[j] = x[j];
                }
                best_x_size = x_size;
            }

        }

        if ( best_x_size == -1) {
            rval = scan_rval;
            goto CLEANUP;
        } else {
            for (j = 0; j < ncols; j++) {
                if (CC_OURABS (best_x[j]) <= INT_MAX) {
                    ortho[j] = (int) best_x[j];
                } else {
                    rval = CC_CHUNK_INTMAT_OVERFLOW_S;
                    goto CLEANUP;
                }
            }
            if (CC_OURABS (best_x[ncols]) <= INT_MAX) {
                *pcol_p = (int) best_x[ncols];
            } else {
                rval = CC_CHUNK_INTMAT_OVERFLOW_S;
                goto CLEANUP;
            }
            rval = 0;
        }

#ifdef  SOLVE_TEST
        for (i = 0; i < m; i++) {
            t = 0;
            for (j = 0; j < n; j++) {
                t += MATRIX(i,j) * best_x[j];
            }
            if ( t )  {
                fprintf (stderr, "WHOA, non-orthogonal orthogonal vector\n");
                rval = CC_CHUNK_INTMAT_ERROR;
                goto CLEANUP;
            }
        }
#endif
    }
    else {
        rval = CC_CHUNK_INTMAT_NOORTHO;
        goto CLEANUP;
    }

CLEANUP:

    if ( rval == CC_CHUNK_INTMAT_OVERFLOW_M ) {
        fprintf (stderr, "Multiply overflow in CCchunk_intmat_ortho\n");
    }
    if ( rval == CC_CHUNK_INTMAT_OVERFLOW_A ) {
        fprintf (stderr, "Addition overflow in CCchunk_intmat_ortho\n");
    }

    return rval;

} /* END CCchunk_intmat_ortho */

static int vector_multiple (int ncols, CCmatval *v1, int *v2)
{
    /* assumes that v1 (at least) is reduced */
    int j;
    CCmatval a;
    CCmatval lim1;

    for (j = 0; j < ncols; j++) {
        if ( v1[j]  ||  v2[j] )  break;
    }

    if ( v1[j] == 0  ||  v2[j] == 0 )  return 0;

    a = v2[j] / v1[j];
    if ( a == 0 ) return 0;

    lim1 = MULTLIM (a);
    for ( ; j < ncols; j++) {
        if ( CC_OURABS (v1[j]) > lim1  ||  v2[j] != a*v1[j] ) {
            return 0;
        }
    }
    return 1;
}

static CCmatval matval_gcd (CCmatval a, CCmatval b)
{
    CCmatval  c;

    if ( a < 0 )  a = - a;
    if ( b < 0 )  b = - b;
    if ( a > b ) {
        c = b;
        b = a;
        a = c;
    }
    while ( a ) {
        c = b % a;
        b = a;
        a = c;
    }
    return b;

} /* END matval_gcd */


static void intmat_reducevec (int beg, int end, int *perm, CCmatval *vec)
{
    int  j;
    CCmatval a = 0;

    if ( perm != NULL ) {
        for (j = beg; j < end; j++) {
            if ( vec[perm[j]] )  a = matval_gcd (a, vec[perm[j]]);
        }
        if ( a > 1 ) {
            for (j = beg; j < end; j++) {
                vec[perm[j]] /= a;
            }
        }
    } else {
        for (j = beg; j < end; j++) {
            if ( vec[j] )  a = matval_gcd (a, vec[j]);
        }
        if ( a > 1 ) {
            for (j = beg; j < end; j++) {
                vec[j] /= a;
            }
        }
    }
} /* END intmat_reducevec */


static int intmat_solve (int m, int n, int *rperm, int *cperm, CCmatval *factor,
                         CCmatval *x, int j_rhs)
{
    int  rval = 0;
    int  i, j, k;
    CCmatval b;
    CCmatval a;
    CCmatval c;
    CCmatval lim;

    for (j = 0; j < n; j++)  x[j] = 0;
    x[cperm[j_rhs]] = -1;

    for (k = m-1; k >= 0; k--) {
        b = FACTOR(rperm[k],cperm[j_rhs]);
        if ( b ) {
            a = FACTOR(rperm[k],cperm[k]);
            c = matval_gcd (a, b);
            a /= c;
            b /= c;
            x[cperm[k]] = b;
            if ( a != 1  &&  a != -1 ) {
                lim = MULTLIM (a);
                if ( CC_OURABS (x[cperm[j_rhs]]) > lim ) {
                    rval = CC_CHUNK_INTMAT_OVERFLOW_M;  goto CLEANUP;
                }
                for (i = m-1; i > k; i--) {
                    if ( CC_OURABS (x[cperm[i]]) > lim ) {
                        rval = CC_CHUNK_INTMAT_OVERFLOW_M;  goto CLEANUP;
                    }
                }
                for (i = 0; i < k; i++) {
                    if ( CC_OURABS (FACTOR(rperm[i],cperm[j_rhs])) > lim ) {
                        rval = CC_CHUNK_INTMAT_OVERFLOW_M;  goto CLEANUP;
                    }
                }
            }
            if ( a != 1 ) {
                x[cperm[j_rhs]] *= a;
                for (i = m-1; i > k; i--) {
                    x[cperm[i]] *= a;
                }
                for (i = 0; i < k; i++) {
                    FACTOR(rperm[i],cperm[j_rhs]) *= a;
                }
            }
            if ( b != 1  &&  b != -1 ) {
                lim = MULTLIM (b);
                for (i = 0; i < k; i++) {
                    if ( CC_OURABS (FACTOR(rperm[i],cperm[k])) > lim ) {
                        rval = CC_CHUNK_INTMAT_OVERFLOW_M;  goto CLEANUP;
                    }
                }
            }
            /* Update rhs */
            for (i = 0; i < k; i++) {
                a = b * FACTOR (rperm[i],cperm[k]);
                if ( (a > 0 && FACTOR(rperm[i],cperm[j_rhs]) < -CC_MATVAL_MAX+a) ||
                     (a < 0 && FACTOR(rperm[i],cperm[j_rhs]) >  CC_MATVAL_MAX+a)  ) {
                    rval = CC_CHUNK_INTMAT_OVERFLOW_A;  goto CLEANUP;
                }
                FACTOR(rperm[i],cperm[j_rhs]) -= a;
            }
        }
    }

CLEANUP:

    return rval;

} /* END intmat_solve */

static void print_matrix (CCchunk_intmat *mat_p)
{
    int i, j;
    int m, n;

    int ncols = mat_p->ncols;
    int nrows = mat_p->nrows;
    int *matrix = mat_p->matrix;

    m = nrows;
    n = ncols + 1;

    for (i=0; i<m; i++) {
        for (j=0; j<n; j++) {
            printf ("%d ",MATRIX(i, j));
        }
        printf ("\n");
    }
    printf ("\n");
}

static void print_factor (CCchunk_intmat *mat_p)
{
    int i, j;
    int m, n;

    int ncols = mat_p->ncols;
    int nrows = mat_p->nrows;
    CCmatval *factor = mat_p->factor;
    int *cperm = mat_p->cperm;
    int *rperm = mat_p->rperm;

    m = nrows;
    n = ncols + 1;

    for (i=0; i<m; i++) {
        for (j=0; j<n; j++) {
            printf ("%4ld ",FACTOR(rperm[i], cperm[j]));
        }
        printf ("\n");
    }
    printf ("\n");
}

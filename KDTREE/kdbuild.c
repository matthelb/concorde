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
/*                   ROUTINE FOR BUILDING KDTREES                           */
/*                                                                          */
/*  (Based on Jon Bentley's paper "K-d trees for semidynamic point sets")   */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 24, 1995 (cofeb24)                                       */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCkdtree_build (CCkdtree *kt, int ncount, CCdatagroup *dat,         */
/*      double *wcoord, CCrandstate *rstate)                                */
/*    -When called, intree should point to a CCkdtree struct that the       */
/*     funtion will load with the tree it builds. The wcoord array          */
/*     is used for node weights (like in Held-Karp), it can be NULL.        */
/*     The node weights must be nonegative (for cutoffs).                   */
/*                                                                          */
/*  void CCkdtree_free (CCkdtree *kt)                                       */
/*    -Frees the space (including the ptrs) used by kt.                     */
/*                                                                          */
/*  void CCkdtree_delete (CCkdtree *kt, int k)                              */
/*    -Deletes the point k from the CCkdtree kt.                            */
/*                                                                          */
/*  void CCkdtree_undelete (CCkdtree *kt, int k)                            */
/*    -Puts the previously deleted point k back into kt.                    */
/*                                                                          */
/*  void CCkdtree_delete_all (CCkdtree *kt, int ncount)                     */
/*    -Deletes all points in kt.                                            */
/*                                                                          */
/*  void CCkdtree_undelete_all (CCkdtree *kt, int ncount)                   */
/*         -Puts all deleted points back in kt. Used to cleanup trees.      */
/*                                                                          */
/*    NOTES:                                                                */
/*       On a 32 bit machine, a CCkdtree on n nodes needs about 52n         */
/*     bytes of memory. CCkdtree_build will return 1 if an error            */
/*     occurs (most likely running out of memory).                          */
/*       CCutil_sprand () should be called before calling                   */
/*     CCkdtree_build ().                                                   */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "kdtree.h"
#include "macrorus.h"

#define CUTOFF 5
#define BNDS_DEPTH 5   /* When bnds info is recorded */
#define BIGDOUBLE (1e30)


static void
    kdtree_free_work (CCkdnode *p, CCptrworld *kdnode_world,
        CCptrworld *kdbnds_world),
    kdtree_free_world (CCptrworld *kdnode_world, CCptrworld *kdbnds_world);
static unsigned char
    findmaxspread (int l, int u, CCkdtree *thetree, double *datx,
           double *daty, double *datw);
static CCkdnode
    *build (int l, int u, int *depth, double *current_bnds_x,
           double *current_bnds_y, CCkdtree *thetree, double *datx,
           double *daty, double *datw, CCrandstate *rstate);


CC_PTRWORLD_ROUTINES (CCkdnode, kdnodealloc, kdnode_bulk_alloc, kdnodefree)
CC_PTRWORLD_LEAKS_ROUTINE (CCkdnode, kdnode_check_leaks, empty, char)

CC_PTRWORLD_ROUTINES (CCkdbnds, kdbndsalloc, kdbnds_bulk_alloc, kdbndsfree)
CC_PTRWORLD_LEAKS_ROUTINE (CCkdbnds, kdbnds_check_leaks, x[0], double)

int CCkdtree_build (CCkdtree *intree, int ncount, CCdatagroup *dat,
        double *wcoord, CCrandstate *rstate)
{
    int i;
    int depth;
    double current_bnds_x[2];
    double current_bnds_y[2];
    CCkdtree *thetree;

    CCptrworld_init (&intree->kdnode_world);
    CCptrworld_init (&intree->kdbnds_world);
    
    if (wcoord != (double *) NULL) {
        for (i = 0; i < ncount; i++) {
            if (wcoord[i] < -0.00000001) {
                fprintf (stderr, "Cannot build with negative node weights\n");
                return 1;
            }
        }
    }

    thetree = intree;
    thetree->perm = CC_SAFE_MALLOC (ncount, int);
    if (!thetree->perm)
        return 1;
    for (i = 0; i < ncount; i++)
        thetree->perm[i] = i;

    thetree->bucketptr = CC_SAFE_MALLOC (ncount, CCkdnode *);
    if (!thetree->bucketptr) {
        CC_FREE (thetree->perm, int);
        return 1;
    }

    depth = 0;
    current_bnds_x[0] = -BIGDOUBLE;
    current_bnds_x[1] =  BIGDOUBLE;
    current_bnds_y[0] = -BIGDOUBLE;
    current_bnds_y[1] =  BIGDOUBLE;

    thetree->root = build (0, ncount - 1, &depth, current_bnds_x,
                     current_bnds_y, thetree, dat->x, dat->y, wcoord, rstate);
    if (!(thetree->root)) {
        fprintf (stderr, "Unable to build CCkdtree\n");
        CC_FREE (thetree->perm, int);
        CC_FREE (thetree->bucketptr, CCkdnode *);
        return 1;
    } else {
        thetree->root->father = (CCkdnode *) NULL;
        return 0;
    }
}

void CCkdtree_free (CCkdtree *kt)
{
    if (kt->perm)
        CC_FREE (kt->perm, int);
    if (kt->bucketptr)
        CC_FREE (kt->bucketptr, CCkdnode *);
    kdtree_free_work (kt->root, &kt->kdnode_world, &kt->kdbnds_world);
    kt->root = (CCkdnode *) NULL;
    
    kdtree_free_world (&kt->kdnode_world, &kt->kdbnds_world);
}

static void kdtree_free_world (CCptrworld *kdnode_world,
        CCptrworld *kdbnds_world)
{
    int total, onlist;

    if (kdnode_check_leaks (kdnode_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding kdnodes\n",
                 total - onlist);
    }
    if (kdbnds_check_leaks (kdbnds_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding kdbnds\n", total - onlist);
    }
    CCptrworld_delete (kdnode_world);
    CCptrworld_delete (kdbnds_world);
}

static void kdtree_free_work (CCkdnode *p, CCptrworld *kdnode_world,
        CCptrworld *kdbnds_world)
{
    if (p->bucket) {
        if (p->bnds)
            kdbndsfree (kdbnds_world, p->bnds);
        kdnodefree (kdnode_world, p);
    } else {
        kdtree_free_work (p->loson, kdnode_world, kdbnds_world);
        kdtree_free_work (p->hison, kdnode_world, kdbnds_world);
        if (p->bnds)
            kdbndsfree (kdbnds_world, p->bnds);
        kdnodefree (kdnode_world, p);
    }
}

static CCkdnode *build (int l, int u, int *depth, double *current_bnds_x,
        double *current_bnds_y, CCkdtree *thetree, double *datx, double *daty,
        double *datw, CCrandstate *rstate)
{
    CCkdnode *p;
    int i, m;
    double savebnd;

    (*depth)++;
    p = kdnodealloc (&thetree->kdnode_world);
    if (!p) {
        (*depth)--;
        return (CCkdnode *) NULL;
    }
    p->empty = 0;

    if (u - l + 1 < CUTOFF) {
        p->bucket = 1;
        p->lopt = l;
        p->hipt = u;
        for (i = l; i <= u; i++)
            thetree->bucketptr[thetree->perm[i]] = p;
        p->bnds = (CCkdbnds *) NULL;
    } else {
        p->bucket = 0;
        if (!((*depth) % BNDS_DEPTH)) {
            p->bnds = kdbndsalloc (&thetree->kdbnds_world);
            if (!p->bnds) {
                (*depth)--;
                kdnodefree (&thetree->kdbnds_world, p);
                return (CCkdnode *) NULL;
            }
            p->bnds->x[0] = current_bnds_x[0];
            p->bnds->x[1] = current_bnds_x[1];
            p->bnds->y[0] = current_bnds_y[0];
            p->bnds->y[1] = current_bnds_y[1];
        } else {
            p->bnds = (CCkdbnds *) NULL;
        }

        p->cutdim = findmaxspread (l, u, thetree, datx, daty, datw);
        m = (l + u) / 2;
        switch (p->cutdim) {
        case 0:
            CCutil_rselect (thetree->perm, l, u, m, datx, rstate);
            p->cutval = datx[thetree->perm[m]];

            savebnd = current_bnds_x[1];
            current_bnds_x[1] = p->cutval;
            p->loson = build (l, m, depth, current_bnds_x, current_bnds_y,
                              thetree, datx, daty, datw, rstate);
            if (!p->loson) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }
            current_bnds_x[1] = savebnd;

            savebnd = current_bnds_x[0];
            current_bnds_x[0] = p->cutval;
            p->hison = build (m + 1, u, depth, current_bnds_x, current_bnds_y,
                              thetree, datx, daty, datw, rstate);
            if (!p->hison) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }
            current_bnds_x[0] = savebnd;

            break;
        case 1:
            CCutil_rselect (thetree->perm, l, u, m, daty, rstate);
            p->cutval = daty[thetree->perm[m]];

            savebnd = current_bnds_y[1];
            current_bnds_y[1] = p->cutval;
            p->loson = build (l, m, depth, current_bnds_x, current_bnds_y,
                              thetree, datx, daty, datw, rstate);
            if (!p->loson) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }
            current_bnds_y[1] = savebnd;

            savebnd = current_bnds_y[0];
            current_bnds_y[0] = p->cutval;
            p->hison = build (m + 1, u, depth, current_bnds_x, current_bnds_y,
                              thetree, datx, daty, datw, rstate);
            if (!p->hison) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }
            current_bnds_y[0] = savebnd;

            break;
        case 2:
            CCutil_rselect (thetree->perm, l, u, m, datw, rstate);
            p->cutval = datw[thetree->perm[m]];

            p->loson = build (l, m, depth, current_bnds_x, current_bnds_y,
                              thetree, datx, daty, datw, rstate);
            if (!p->loson) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }

            p->hison = build (m + 1, u, depth, current_bnds_x, current_bnds_y,
                              thetree, datx, daty, datw, rstate);
            if (!p->hison) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }

            break;
        }
        p->loson->father = p;
        p->hison->father = p;
    }
    (*depth)--;
    return p;
}

static unsigned char findmaxspread (int l, int u, CCkdtree *thetree,
        double *datx, double *daty, double *datw)
{
    int i;
    double xmax, xmin, xval, xspread;
    double ymax, ymin, yval, yspread;
    double wmax, wmin, wval, wspread;

    wmax = (double) 0.0;
    wmin = (double) 0.0;

    if (datw != (double *) NULL) {
        wmin = datw[thetree->perm[l]];
        wmax = wmin;
    }
    xmin = datx[thetree->perm[l]];
    xmax = xmin;
    ymin = daty[thetree->perm[l]];
    ymax = ymin;
    for (i = l + 1; i <= u; i++) {
        xval = datx[thetree->perm[i]];
        if (xval < xmin)
            xmin = xval;
        else if (xval > xmax)
            xmax = xval;
        yval = daty[thetree->perm[i]];
        if (yval < ymin)
            ymin = yval;
        else if (yval > ymax)
            ymax = yval;
        if (datw != (double *) NULL) {
            wval = datw[thetree->perm[i]];
            if (wval < wmin)
                wmin = wval;
            else if (wval > wmax)
                wmax = wval;
        }
    }

    xspread = xmax - xmin;
    yspread = ymax - ymin;

    if (datw != (double *) NULL) {
        wspread = (wmax - wmin);
        if (xspread >= yspread && xspread >= wspread)
            return (unsigned char) 0;
        else if (yspread >= xspread && yspread >= wspread)
            return (unsigned char) 1;
        else {
            return (unsigned char) 2;
        }
    } else {
        if (xspread >= yspread)
            return (unsigned char) 0;
        else
            return (unsigned char) 1;
    }
}

void CCkdtree_delete (CCkdtree *kt, int k)
{
    int j, temp;
    CCkdnode *p;

    p = kt->bucketptr[k];
    j = p->lopt;
    while (kt->perm[j] != k)
        j++;
    CC_SWAP (kt->perm[j], kt->perm[p->hipt], temp);
    (p->hipt)--;
    if (p->lopt > p->hipt) {
        p->empty = 1;
        while ((p = p->father) != (CCkdnode *) NULL &&
                p->loson->empty && p->hison->empty)
            p->empty = 1;
    }
}

void CCkdtree_delete_all (CCkdtree *kt, int ncount)
{
    int k;

    for (k = 0; k < ncount; k++)
        CCkdtree_delete (kt, k);
}


void CCkdtree_undelete (CCkdtree *kt, int k)
{
    int j, temp;
    CCkdnode *p;

    p = kt->bucketptr[k];
    j = p->lopt;
    while (kt->perm[j] != k)
        j++;
    if (j > p->hipt) {
        (p->hipt)++;
        CC_SWAP (kt->perm[j], kt->perm[p->hipt], temp);
        if (p->empty) {
            p->empty = 0;
            while ((p = p->father) != (CCkdnode *) NULL && p->empty)
                p->empty = 0;
        }
    }
}

void CCkdtree_undelete_all (CCkdtree *kt, int ncount)
{
    int k;

    for (k = 0; k < ncount; k++)
        CCkdtree_undelete (kt, k);
}

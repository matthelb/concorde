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
/*                   UTILITY ROUTINES FOR SUBDIVISION                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: March 17, 2003                                                    */
/*  Changes:                                                                */
/*                                                                          */
/*                                                                          */
/*  EXPORTED FUNCTIONS:                                                     */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "macrorus.h"

typedef struct kpnode {
    struct kpnode *loson;
    struct kpnode *hison;
    struct kpnode *next;
    int lopt;
    int hipt;
    double x[2];
    double y[2];
} kpnode;

typedef struct kptree {
    kpnode  *root;
    int     *perm;
    CCptrworld kpnode_world;
} kptree;

CC_PTRWORLD_ROUTINES (kpnode, kpnodealloc, kpnode_bulk_alloc, kpnodefree)

static int
    findmaxspread (int l, int u, kptree *thetree, double *datx,
        double *daty),
    grab_parts (kptree *thetree, kpnode *n, int *ind, CCsubdiv *slist,
        int **plist);

static kpnode *build (int l, int u, double xlo, double xhi, double ylo,
        double yhi, kptree *thetree, double *datx, double *daty,
        int bucketsize, int *cnt, CCrandstate *rstate);


int CCutil_karp_partition (int ncount, CCdatagroup *dat, int partsize,
        int *p_scount, CCsubdiv **p_slist, int ***partlist,
        CCrandstate *rstate)
{
    int i, m, mult, target = 0, rval = 0;
    int norm = dat->norm;
    double xlo, xhi, ylo, yhi;
    int scount = 0;
    kptree thetree;
    kpnode *p, *q;
    double pcutval, qcutval;

    printf ("Create a Karp partition with bucketsize %d\n", partsize);
    fflush (stdout);

    *p_slist = (CCsubdiv *) NULL;
    *partlist = (int **) NULL;

    CCptrworld_init (&thetree.kpnode_world);
    thetree.perm = (int *) NULL;
    thetree.root = (kpnode *) NULL;

    if ((norm & CC_NORM_SIZE_BITS) != CC_D2_NORM_SIZE) {
        fprintf (stderr, "Only set up for 2D norms\n");
        rval = 1;  goto CLEANUP;
    }

    if (ncount < 2*partsize) {
        fprintf (stderr, "two few nodes for partition size\n");
        rval = 1;  goto CLEANUP;
    }

    if (norm == CC_GEOM) {
        mult = 1;
        do {
            mult *= 2;
        } while (ncount/(mult+2) > partsize);
        target = ncount/(mult+2);
        printf ("With polar caps, using bucketsize %d\n", target);
    }

    thetree.perm = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (thetree.perm, "out of memory in CCutil_karp_partitiion");

    for (i = 0; i < ncount; i++) {
        thetree.perm[i] = i;
    }

    xlo = xhi = dat->x[0];
    ylo = yhi = dat->y[0];
    for (i = 0; i < ncount; i++) {
        if (dat->x[i] < xlo) xlo = dat->x[i];
        if (dat->x[i] > xhi) xhi = dat->x[i];
        if (dat->y[i] < ylo) ylo = dat->y[i];
        if (dat->y[i] > yhi) yhi = dat->y[i];
    }

    if (norm == CC_GEOM) {
        CCutil_rselect (thetree.perm, 0, ncount-1, target, dat->x, rstate);
        pcutval = dat->x[thetree.perm[target]];

        p = kpnodealloc (&thetree.kpnode_world);
        CCcheck_NULL (p, "pnodealloc failed");

        thetree.root = p;

        p->loson = build (0, target, xlo, pcutval, ylo, yhi, &thetree, dat->x,
                          dat->y, partsize, &scount, rstate);
        if (!p->loson) {
            fprintf (stderr, "initial loson build failed\n");
            rval = 1;  goto CLEANUP;
        } 

        q = kpnodealloc (&thetree.kpnode_world);
        CCcheck_NULL (q, "pnodealloc failed");

        p->hison = q;

        m = ncount-1-target;
        CCutil_rselect (thetree.perm, target+1, ncount-1, m, dat->x, rstate);
        qcutval = dat->x[thetree.perm[m]];

        q->hison = build (m+1, ncount-1, qcutval, xhi, ylo, yhi, &thetree,
                          dat->x, dat->y, partsize, &scount, rstate);
        if (!q->hison) {
            fprintf (stderr, "q hison build failed\n");
            rval = 1;  goto CLEANUP;
        } 

        q->loson = build (target+1, m, pcutval, qcutval, ylo, yhi, &thetree,
                          dat->x, dat->y, partsize, &scount, rstate);
        if (!q->loson) {
            fprintf (stderr, "q loson build failed\n");
            rval = 1;  goto CLEANUP;
        } 
    } else {
        thetree.root = build (0, ncount-1, xlo, xhi, ylo, yhi, &thetree,
                              dat->x, dat->y, partsize, &scount, rstate);
        if (!thetree.root) {
            fprintf (stderr, "unable to build partition tree\n");
            rval = 1;  goto CLEANUP;
        }
    }

    *partlist = CC_SAFE_MALLOC (scount, int *);
    CCcheck_NULL (*partlist, "out of memory in CCutil_karp_partition");
    for (i = 0; i < scount; i++) {
        (*partlist)[i] = (int *) NULL;
    }

    *p_slist = CC_SAFE_MALLOC (scount, CCsubdiv);
    CCcheck_NULL (*p_slist, "out of memory in CCutil_karp_partition");
    
    i = 0;
    rval = grab_parts (&thetree, thetree.root, &i, *p_slist, *partlist);
    CCcheck_rval (rval, "grab_parts failed");

    *p_scount = scount;

CLEANUP:

    CCptrworld_delete (&thetree.kpnode_world);
    CC_IFFREE (thetree.perm, int);
    if (rval) {
        if (*partlist) {
            for (i = 0; i < scount; i++) {
                CC_IFFREE ((*partlist)[i], int);
            }
            CC_FREE (*partlist, int *);
        }
        CC_IFFREE (*p_slist, CCsubdiv);
    }

    return rval;
}

static kpnode *build (int l, int u, double xlo, double xhi, double ylo,
        double yhi, kptree *thetree, double *datx, double *daty,
        int bucketsize, int *cnt, CCrandstate *rstate)
{
    int rval = 0;
    int m, cutdim;
    double cutval;
    kpnode *p;

    p = kpnodealloc (&thetree->kpnode_world);
    CCcheck_NULL (p, "pnodealloc failed");

    if (u - l + 1 <= bucketsize) {
        p->hison = (kpnode *) NULL;
        p->loson = (kpnode *) NULL;
        p->x[0] = xlo;
        p->x[1] = xhi;
        p->y[0] = ylo;
        p->y[1] = yhi;
        p->lopt = l;
        p->hipt = u;
        printf ("Part %d [%.2f, %.2f] [%.2f, %.2f]: %d points\n", 
                    *cnt, xlo, xhi, ylo, yhi, u - l + 1);
        fflush (stdout);
        (*cnt)++;
    } else {
        cutdim = findmaxspread (l, u, thetree, datx, daty);
        m = (l + u) / 2;
        if (cutdim == 0) {
            CCutil_rselect (thetree->perm, l, u, m, datx, rstate);
            cutval = datx[thetree->perm[m]];
            p->loson = build (l, m, xlo, cutval, ylo, yhi, thetree, datx,
                              daty, bucketsize, cnt, rstate);
            if (!p->loson) {
                rval = 1;  goto CLEANUP;
            }
            p->hison = build (m+1, u, cutval, xhi, ylo, yhi, thetree, datx,
                              daty, bucketsize, cnt, rstate);
            if (!p->hison) {
                rval = 1;  goto CLEANUP;
            }
        } else {
            CCutil_rselect (thetree->perm, l, u, m, daty, rstate);
            cutval = daty[thetree->perm[m]];
            p->loson = build (l, m, xlo, xhi, ylo, cutval, thetree, datx,
                              daty, bucketsize, cnt, rstate);
            if (!p->loson) {
                rval = 1;  goto CLEANUP;
            }
            p->hison = build (m+1, u, xlo, xhi, cutval, yhi, thetree, datx,
                              daty, bucketsize, cnt, rstate);
            if (!p->hison) {
                rval = 1;  goto CLEANUP;
            }
        }
    }

CLEANUP:

    if (rval) return (kpnode *) NULL;
    else     return p;
}

static int findmaxspread (int l, int u, kptree *thetree, double *datx,
        double *daty)
{
    int i;
    double xmax, xmin, xval, xspread;
    double ymax, ymin, yval, yspread;

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
    }

    xspread = xmax - xmin;
    yspread = ymax - ymin;

    if (xspread >= yspread) return 0;
    else                    return 1;
}

static int grab_parts (kptree *thetree, kpnode *n, int *ind, CCsubdiv *slist,
        int **plist)
{
    int rval = 0;
    int i, k, cnt;

    if (n->loson) {
        rval = grab_parts (thetree, n->loson, ind, slist, plist);
        if (rval) goto CLEANUP;

        rval = grab_parts (thetree, n->hison, ind, slist, plist);
        if (rval) goto CLEANUP;
    } else {
        cnt = n->hipt - n->lopt + 1;

        k = *ind;
        plist[k] = CC_SAFE_MALLOC (cnt, int);
        for (i = 0; i < cnt; i++) {
            plist[k][i] = thetree->perm[n->lopt + i];
        }
        slist[k].id = k;
        slist[k].xrange[0] = n->x[0];
        slist[k].xrange[1] = n->x[1];
        slist[k].yrange[0] = n->y[0];
        slist[k].yrange[1] = n->y[1];
        slist[k].cnt = cnt;
        slist[k].bound = -1.0;

        (*ind)++;
    }

CLEANUP:

    return rval;
}

int CCutil_write_subdivision_index (char *problabel, int ncount, int scount,
        CCsubdiv *slist)
{
    FILE *f = (FILE *) NULL;
    char *index_name = (char *) NULL;
    char *new_name = (char *) NULL;
    char *back_name = (char *) NULL;
    int tval, rval = 0;
    int i;
    size_t len;

    len = strlen(problabel) + 10;
    index_name = CC_SAFE_MALLOC (len, char);
    new_name = CC_SAFE_MALLOC (len, char);
    back_name = CC_SAFE_MALLOC (len, char);
    if (index_name == (char *) NULL ||
        new_name == (char *) NULL ||
        back_name == (char *) NULL) {
        fprintf (stderr, "Out of memory in CCutil_write_subdivision_index\n");
        rval = 1; goto CLEANUP;
    }

    strcpy (index_name, problabel);
    strcat (index_name, ".index");
    strcpy (new_name, "N");
    strcat (new_name, index_name);
    strcpy (back_name, "O");
    strcat (back_name, index_name);
    
    f = fopen (new_name, "w");
    if (f == (FILE*) NULL) {
        perror (new_name);
        fprintf (stderr,
           "Unable to open %s for output in CCutil_write_subdivision_index\n",
            new_name);
        rval = 1; goto CLEANUP;
    }

    tval = fprintf (f, "%s %d\n", problabel, ncount);
    if (tval <= 0) {
        perror (new_name);
        fprintf (stderr, "fprintf to %s failed\n", new_name);
        rval = 1; goto CLEANUP;
    }

    fprintf (f, "%d\n", scount);
    for (i = 0; i < scount; i++) {
        fprintf (f, "%d %lf %lf %lf %lf %d %lf\n",
                      slist[i].id, slist[i].xrange[0], slist[i].xrange[1],
                      slist[i].yrange[0], slist[i].yrange[1],
                      slist[i].cnt, slist[i].bound);
    }
    
    tval = fclose (f);
    if (tval) {
        perror (new_name);
        fprintf (stderr, "fclose %s failed\n", new_name);
        rval = 1; goto CLEANUP;
    }
    f = (FILE *) NULL;

    rename (index_name, back_name);
    tval = rename (new_name, index_name);
    if (tval) {
        perror (index_name);
        fprintf (stderr, "rename %s to %s failed\n", new_name, index_name);
        rval = 1; goto CLEANUP;
    }
    
CLEANUP:

    if (f != (FILE *) NULL) {
        fclose (f);
    }
    CC_IFFREE (new_name, char);
    CC_IFFREE (index_name, char);
    CC_IFFREE (back_name, char);
    return rval;
}

int CCutil_read_subdivision_index (char *index_name, char **p_problabel,
        int *p_ncount, int *p_scount, CCsubdiv **p_slist)
{
    FILE *f = (FILE *) NULL;
    char *problabel = (char *) NULL;
    CCsubdiv *slist = (CCsubdiv *) NULL;
    int i, scount, rval = 0, tval;

    f = fopen (index_name, "r");
    if (f == (FILE*) NULL) {
        perror (index_name);
        fprintf (stderr,
             "Unable to open %s for input in CCutil_read_subdivision_index\n",
                 index_name);
        rval = 1; goto CLEANUP;
    }

    problabel = CC_SAFE_MALLOC (1024, char);
    if (problabel == (char *) NULL) {
        fprintf (stderr, "Out of memory in read_index\n");
        rval = 1; goto CLEANUP;
    }

    CCutil_readstr (f, problabel, 1024);
    
    tval = fscanf (f, "%d\n", p_ncount);
    if (tval <= 0) {
        perror (index_name);
        fprintf (stderr, "fscanf from %s failed\n", index_name);
        rval = 1; goto CLEANUP;
    }

    fscanf (f, "%d", &scount);
    if (scount <= 0) {
        fprintf (f, "Bad count in index file: %d\n", scount);
        rval = 1; goto CLEANUP;
    }
    slist = CC_SAFE_MALLOC (scount, CCsubdiv); 
    CCcheck_NULL (slist, "out of memory in CCutil_read_subdivision_index");

    for (i = 0; i < scount; i++) {
        fscanf (f, "%d %lf %lf %lf %lf %d %lf\n",
                      &slist[i].id, &slist[i].xrange[0], &slist[i].xrange[1],
                      &slist[i].yrange[0], &slist[i].yrange[1],
                      &slist[i].cnt, &slist[i].bound);
    }
    
    tval = fclose (f);
    if (tval) {
        perror (index_name);
        fprintf (stderr, "fclose %s failed\n", index_name);
        rval = 1; goto CLEANUP;
    }
    f = (FILE *) NULL;

    *p_problabel = problabel;
    *p_scount = scount;
    *p_slist = slist;
    
CLEANUP:

    if (f != (FILE *) NULL) {
        fclose (f);
    }
    if (rval) {
        CC_IFFREE (problabel, char);
        CC_IFFREE (slist, CCsubdiv);
    }
    return rval;
}

int CCutil_write_subdivision_lkh_index (char *problabel, int ncount,
        int scount, CCsubdiv_lkh *slist, double tourlen)
{
    FILE *f = (FILE *) NULL;
    char *index_name = (char *) NULL;
    char *new_name = (char *) NULL;
    char *back_name = (char *) NULL;
    int tval, rval = 0;
    int i;
    size_t len;

    len = strlen(problabel) + 10;
    index_name = CC_SAFE_MALLOC (len, char);
    new_name = CC_SAFE_MALLOC (len, char);
    back_name = CC_SAFE_MALLOC (len, char);
    if (index_name == (char *) NULL ||
        new_name == (char *) NULL ||
        back_name == (char *) NULL) {
        fprintf (stderr, "Out of memory in CCutil_write_subdivision_index\n");
        rval = 1; goto CLEANUP;
    }

    strcpy (index_name, problabel);
    strcat (index_name, ".lindex");
    strcpy (new_name, "N");
    strcat (new_name, index_name);
    strcpy (back_name, "O");
    strcat (back_name, index_name);

    f = fopen (new_name, "w");
    if (f == (FILE*) NULL) {
        perror (new_name);
        fprintf (stderr,
           "Unable to open %s for output in CCutil_write_subdivision_index\n",
            new_name);
        rval = 1; goto CLEANUP;
    }

    tval = fprintf (f, "%s %d %0.0f\n", problabel, ncount, tourlen);
    if (tval <= 0) {
        perror (new_name);
        fprintf (stderr, "fprintf to %s failed\n", new_name);
        rval = 1; goto CLEANUP;
    }

    fprintf (f, "%d\n", scount);
    for (i = 0; i < scount; i++) {
        fprintf (f, "%d %d %d %.0f %.0f\n",
                      slist[i].id, slist[i].cnt, slist[i].start,
                      slist[i].origlen, slist[i].newlen);
    }

    tval = fclose (f);  
    if (tval) {
        perror (new_name);
        fprintf (stderr, "fclose %s failed\n", new_name);
        rval = 1; goto CLEANUP;
    }
    f = (FILE *) NULL;

    rename (index_name, back_name);
    tval = rename (new_name, index_name);
    if (tval) {
        perror (index_name);
        fprintf (stderr, "rename %s to %s failed\n", new_name, index_name);
        rval = 1; goto CLEANUP;
    }
    
CLEANUP:

    if (f != (FILE *) NULL) {
        fclose (f);
    }
    CC_IFFREE (new_name, char);
    CC_IFFREE (index_name, char);
    CC_IFFREE (back_name, char);
    return rval;
}

int CCutil_read_subdivision_lkh_index (char *index_name, char **p_problabel,
        int *p_ncount, int *p_scount, CCsubdiv_lkh **p_slist,
        double *p_tourlen)
{
    FILE *f = (FILE *) NULL;
    char *problabel = (char *) NULL;
    CCsubdiv_lkh *slist = (CCsubdiv_lkh *) NULL;
    int i, scount, rval = 0, tval;

    f = fopen (index_name, "r");
    if (f == (FILE*) NULL) {
        perror (index_name);
        fprintf (stderr,
             "Unable to open %s for input in CCutil_read_subdivision_index\n",
                 index_name);
        rval = 1; goto CLEANUP;
    }

    problabel = CC_SAFE_MALLOC (1024, char);
    if (problabel == (char *) NULL) {
        fprintf (stderr, "Out of memory in read_index\n");
        rval = 1; goto CLEANUP;
    }

    CCutil_readstr (f, problabel, 1024);
    
    tval = fscanf (f, "%d %lf\n", p_ncount, p_tourlen);
    if (tval <= 0) {
        perror (index_name);
        fprintf (stderr, "fscanf from %s failed\n", index_name);
        rval = 1; goto CLEANUP;
    }

    fscanf (f, "%d", &scount);
    if (scount <= 0) {
        fprintf (f, "Bad count in index file: %d\n", scount);
        rval = 1; goto CLEANUP;
    }
    slist = CC_SAFE_MALLOC (scount, CCsubdiv_lkh); 
    CCcheck_NULL (slist, "out of memory in CCutil_read_subdivision_index");

    for (i = 0; i < scount; i++) {
        fscanf (f, "%d %d %d %lf %lf\n",
                      &slist[i].id, &slist[i].cnt, &slist[i].start,
                      &slist[i].origlen, &slist[i].newlen);
    }
    
    tval = fclose (f);
    if (tval) {
        perror (index_name);
        fprintf (stderr, "fclose %s failed\n", index_name);
        rval = 1; goto CLEANUP;
    }
    f = (FILE *) NULL;

    *p_problabel = problabel;
    *p_scount = scount;
    *p_slist = slist;
    
CLEANUP:

    if (f != (FILE *) NULL) {
        fclose (f);
    }
    if (rval) {
        CC_IFFREE (problabel, char);
        CC_IFFREE (slist, CCsubdiv_lkh);
    }
    return rval;
}

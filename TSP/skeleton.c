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
/*                  ROUTINES TO MAINTAIN SKELETONS                          */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: October 19, 1997                                                  */
/*        January 28, 2003 (bico)                                           */
/*                                                                          */
/*  These routines maintain skeletons for CCtsp_lpcut and CCtsp_lpcut_in.   */
/*  A skeleton is a set of atoms of the cut, and is represented by a        */
/*  sorted array of the lowest numbered node in each atom.  The only        */
/*  function that cares about which representative is chosen for each       */
/*  atom, or the order of the atoms is CCtsp_compare_skeleton.              */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  void CCtsp_init_skeleton (CCtsp_skeleton *skel)                         */
/*    INITIALIZES a skeleton (to the NULL skeleton)                         */
/*                                                                          */
/*  void CCtsp_free_skeleton (CCtsp_skeleton *skel)                         */
/*    FREES the memory used by skel                                         */
/*                                                                          */
/*  int CCtsp_copy_skeleton (CCtsp_skeleton *old, CCtsp_skeleton *new)      */
/*    COPIES a skeleton                                                     */
/*                                                                          */
/*  int CCtsp_read_skeleton (CC_SFILE *f, CCtsp_skeleton *skel,             */
/*      int ncount)                                                         */
/*    READS a skeleton from f                                               */
/*                                                                          */
/*  int CCtsp_write_skeleton (CC_SFILE *f, CCtsp_skeleton *skel,            */
/*      int ncount)                                                         */
/*    WRITES a skeleton to f                                                */
/*                                                                          */
/*  int CCtsp_construct_skeleton (CCtsp_lpcut_in *c, int nodecount)         */
/*    CONSTRUCTS a skeleton for c, representing all atoms in c              */
/*                                                                          */
/*  void CCtsp_compare_skeletons (CCtsp_skeleton *a, CCtsp_skeleton *b,     */
/*      int *diff)                                                          */
/*    COMPARES two skeletons, setting *diff=0 if they are the same, and     */
/*    *diff=1 if they are not.                                              */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"

#define SKELETON_WILD 0

#undef  DEBUG_CONSTRUCT

void CCtsp_init_skeleton (CCtsp_skeleton *skel)
{
    skel->atomcount = 0;
    skel->atoms = (int *) NULL;
}

void CCtsp_free_skeleton (CCtsp_skeleton *skel)
{
    if (skel != (CCtsp_skeleton *) NULL) {
        skel->atomcount = 0;
        CC_IFFREE (skel->atoms, int);
    }
}

int CCtsp_copy_skeleton (CCtsp_skeleton *old, CCtsp_skeleton *new)
{
    int i;

    CCtsp_init_skeleton (new);

    if (old->atomcount == 0) return 0;
    new->atoms = CC_SAFE_MALLOC (old->atomcount, int);
    if (new->atoms == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCtsp_copy_skeleton\n");
        return 1;
    }
    for (i=0; i<old->atomcount; i++) {
        new->atoms[i] = old->atoms[i];
    }
    new->atomcount = old->atomcount;

    return 0;
}

int CCtsp_read_skeleton (CC_SFILE *f, CCtsp_skeleton *skel, int ncount)
{
    int rval;
    int atomcount;
    int i;
    int nbits = CCutil_sbits (ncount);
    char type;
    
    CCtsp_init_skeleton (skel);

    rval = CCutil_sread_char (f, &type);
    if (rval) {
        fprintf (stderr, "CCutil_sread_char failed\n");
        goto CLEANUP;
    }

    switch (type) {
    case SKELETON_WILD:
        rval = CCutil_sread_bits (f, &atomcount, nbits);
        if (rval) {
            fprintf (stderr, "CCutil_sread_bits failed\n");
            goto CLEANUP;
        }
        skel->atomcount = atomcount;
        
        if (atomcount == 0) {
            skel->atoms = (int *) NULL;
            rval = 0;
            goto CLEANUP;
        }
        
        skel->atoms = CC_SAFE_MALLOC (atomcount, int);
        if (skel->atoms == (int *) NULL) {
            fprintf (stderr, "Out of memory in CCtsp_read_skeleton\n");
            rval = 1; goto CLEANUP;
        }

        for (i=0; i<atomcount; i++) {
            rval = CCutil_sread_bits (f, &skel->atoms[i], nbits);
            if (rval) {
                fprintf (stderr, "CCutil_sread_bits failed\n");
                goto CLEANUP;
            }
        }
        break;
    default:
        fprintf (stderr, "Unknown skeleton type %ud\n", (unsigned) type);
        rval = 1; goto CLEANUP;
    }
    rval = 0;

 CLEANUP:
    if (rval) {
        CCtsp_free_skeleton (skel);
    }
    return rval;
}

int CCtsp_write_skeleton (CC_SFILE *f, CCtsp_skeleton *skel, int ncount)
{
    int rval;
    int atomcount = skel->atomcount;
    int i;
    int nbits = CCutil_sbits (ncount);

    rval = CCutil_swrite_char (f, SKELETON_WILD);
    if (rval) {
        fprintf (stderr, "CCutil_swrite_char failed\n");
        goto CLEANUP;
    }
    
    rval = CCutil_swrite_bits (f, atomcount, nbits);
    if (rval) {
        fprintf (stderr, "CCutil_swrite_bits failed\n");
        goto CLEANUP;
    }

    for (i=0; i<atomcount; i++) {
        rval = CCutil_swrite_bits (f, skel->atoms[i], nbits);
        if (rval) {
            fprintf (stderr, "CCutil_swrite_bits failed\n");
            goto CLEANUP;
        }
    }

    rval = 0;

 CLEANUP:
    return rval;
}

int CCtsp_construct_skeleton (CCtsp_lpcut_in *c, int nodecount)
{
    int cliquecount = c->cliquecount;
    CCtsp_lpclique *cliques = c->cliques;
    int *label = (int *) NULL;
    int ccount;
    int *cnodes = (int *) NULL;
    int atomcount;
    int atomcount_save;
    int *atomsize = (int *) NULL;
    int *atomnew = (int *) NULL;
    int *atomwork = (int *) NULL;
    int *atoms = (int *) NULL;
    int i;
    int j;
    int tmp;
    int rval = 0;

    if (c->dominocount != 0) {
        printf ("Skeleton Yipes %d\n", c->dominocount);
        fflush (stdout);
        exit (1);
    }
    
    CCtsp_init_skeleton (&c->skel);
    if (c->dominocount > 0) goto CLEANUP;   /* don't build for dominos */

    label = CC_SAFE_MALLOC (nodecount, int);
    if (label == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCtsp_construct_skeleton\n");
        rval = 1; goto CLEANUP;
    }

    /* initialize labels */
    for (i=0; i<cliquecount; i++) {
        CC_FOREACH_NODE_IN_CLIQUE (j, cliques[i], tmp) {
            label[j] = 0;
        }
    }

    /* count nodes */
    ccount = 0;
    for (i=0; i<cliquecount; i++) {
        CC_FOREACH_NODE_IN_CLIQUE (j, cliques[i], tmp) {
            if (label[j] == 0) {
                label[j] = 1;
                ccount++;
            }
        }
    }

    cnodes   = CC_SAFE_MALLOC (ccount, int);
    atomsize = CC_SAFE_MALLOC (ccount+1, int);
    atomnew = CC_SAFE_MALLOC (ccount+1, int);
    atomwork = CC_SAFE_MALLOC (ccount+1, int);
    if (cnodes   == (int *) NULL ||
        atomsize == (int *) NULL ||
        atomnew == (int *) NULL ||
        atomwork == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCtsp_construct_skeleton\n");
        rval = 1; goto CLEANUP;
    }

    /* collect nodes */
    ccount = 0;
    for (i=0; i<cliquecount; i++) {
        CC_FOREACH_NODE_IN_CLIQUE (j, cliques[i], tmp) {
            if (label[j] == 1) {
                label[j] = 0;
                cnodes[ccount++] = j;
            }
        }
    }

    CCutil_int_array_quicksort (cnodes, ccount);

    /* refine atoms */
    atomsize[0] = ccount;
    atomcount = 1;
    for (i=0; i<cliquecount; i++) {
        for (j=0; j<atomcount; j++) {
            atomwork[j] = 0;
        }
        CC_FOREACH_NODE_IN_CLIQUE (j, cliques[i], tmp) {
            atomwork[label[j]]++;
        }
        atomcount_save = atomcount;
        for (j=0; j<atomcount_save; j++) {
            if (atomwork[j] == 0) {
                atomnew[j] = -1;
            } else if (atomwork[j] == atomsize[j]) {
                atomnew[j] = j;
            } else {
                atomsize[atomcount] = atomwork[j];
                atomsize[j] -= atomwork[j];
                atomnew[j] = atomcount;
                atomcount++;
            }
        }
        CC_FOREACH_NODE_IN_CLIQUE (j, cliques[i], tmp) {
            label[j] = atomnew[label[j]];
        }
    }

    atomcount_save = atomcount;
    if (ccount < nodecount) {
        atomcount += 1;
    }
    atoms = CC_SAFE_MALLOC (atomcount, int);
    if (atoms == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCtsp_construct_skeleton\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<atomcount; i++) {
        atoms[i] = -1;
    }

    /* find representatives */
    for (i=0; i<ccount; i++) {
        if (atoms[label[cnodes[i]]] == -1) {
            atoms[label[cnodes[i]]] = cnodes[i];
        }
    }

    if (ccount < nodecount) {
        if (cnodes[ccount-1] == ccount-1) {
            atoms[atomcount_save] = ccount;
        } else {
            for (i=0; i<ccount; i++) {
                if (cnodes[i] != i) {
                    atoms[atomcount_save] = i;
                    break;
                }
            }
        }
    }

    CCutil_int_array_quicksort (atoms, atomcount);
    
    c->skel.atoms = atoms;
    c->skel.atomcount = atomcount;
    atoms = (int *) NULL;

#ifdef DEBUG_CONSTRUCT
    printf ("skeleton:");
    for (i=0; i<c->skel.atomcount; i++) {
        printf (" %d", c->skel.atoms[i]);
    }
    printf ("\n");
    fflush (stdout);
#endif

    rval = 0;

 CLEANUP:
    CC_IFFREE (atoms, int);
    CC_IFFREE (atomwork, int);
    CC_IFFREE (atomnew, int);
    CC_IFFREE (atomsize, int);
    CC_IFFREE (cnodes, int);
    CC_IFFREE (label, int);
    if (rval) {
        CCtsp_free_skeleton (&c->skel);
    }
    return rval;
}

void CCtsp_compare_skeletons (CCtsp_skeleton *a, CCtsp_skeleton *b, int *diff)
{
    int i;
    
    *diff = 1;
    if (a->atomcount != b->atomcount) {
        return;
    }
    for (i=0; i<a->atomcount; i++) {
        if (a->atoms[i] != b->atoms[i]) {
            return;
        }
    }
    *diff = 0;
}

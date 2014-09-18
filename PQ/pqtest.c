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

#include "machdefs.h"
#include "util.h"
#include "pq.h"

int main (int ac, char **av)
{
    int nnodes = 20;
    CCpq_tree pqt;
    CCpq_node *n;
    int i;
    int rval;
    int status;

    printf ("%d\n", (int) (sizeof (CCpq_node)));
    if (ac > 1) {
        nnodes = atoi(av[1]);
    }
    CCpq_tree_init (&pqt);

    rval = CCpq_tree_trivial (&pqt, nnodes, 0);
    if (rval) {
        fprintf (stderr, "CCpq_tree_trivial failed\n");
        goto CLEANUP;
    }

    for (;;) {
        CCpq_dump_solution (&pqt);
        CCpq_clear_leaflist (&pqt);
        while (scanf("%d",&i) == 1 && i != -1) {
            CCpq_add_leaflist (&pqt, i);
        }
        if (pqt.leaflist == (CCpq_node *) NULL) break;
        printf ("adding:");
        for (n = pqt.leaflist; n; n = n->next) {
            printf (" %d", (int) (n - pqt.elems));
        }
        fflush (stdout);
        rval = CCpq_apply (&pqt, &status);
        if (rval) {
            fprintf (stderr, "CCpq_apply failed\n");
            goto CLEANUP;
        }
        if (status == CCpq_STATUS_TRIVIAL) {
            printf ("  (trivial)\n");
        } else if (status == CCpq_STATUS_NONTRIVIAL) {
            printf (" (nontrivial)\n");
        } else if (status == CCpq_STATUS_NOSOL) {
            printf (" (failed)\n");
        } else {
            fprintf (stderr, "Unknown PQ status %d\n", status);
            rval = 1; goto CLEANUP;
        }
    }
    rval = 0;
  CLEANUP:
    CCpq_tree_free (&pqt);
    return rval;
}

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
/*                        Teething Combs                                    */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: March 13, 1997                                                    */
/*        July 15, 1997                                                     */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCtsp_teething (CCtsp_lpgraph *g, double *x,                        */
/*      CCtsp_lpcut_in *cut, CCtsp_lpcut_in **newcut)                       */
/*       CALLS teething with the handle and bigteeth of the comb given by   */
/*        cut (the code will test that it really is a comb or subtour).     */
/*        -g is the graph of the active edges in the LP                     */
/*        -x is an LP solution                                              */
/*        -cut is the base comb                                             */
/*        -newcut returns the comb found after teething                     */
/*    NOTES:                                                                */
/*        The new cut may be just a subtour inequality. This code is based  */
/*    on the "Teething" section of the tsp notes.                           */
/*                                                                          */
/*  int CCtsp_teething_list (CCtsp_lpgraph *g, double *x,                   */
/*      CCtsp_lpclique *handle, int nbig, CCtsp_lpclique **bigteeth,        */
/*      CCtsp_lpcut_in **newcut)                                            */
/*    CALLS teething with the given handle and bigteeth.                    */
/*    NOTES: bigteeth should be filled in from 1 to nbig (not 0 to          */
/*           nbig-1).                                                       */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"

typedef struct intptr {
    int this;
    struct intptr *next;
} intptr;

typedef struct Rrecord {
    int add0;
    int add1;
    struct Rrecord *next;
} Rrecord;



static void
    teething_free_world (CCptrworld *intptr_world, CCptrworld *Rrecord_world),
    identify_big_teeth (CCtsp_lpcut_in *cut, int handle, int *nbig,
            CCtsp_lpclique **bigteeth);

static int
    optimal_pseudocomb (CCtsp_lpgraph *g, double *x, CCtsp_lpclique *handle,
        int nbig, CCtsp_lpclique **bigteeth, CCtsp_lpcut_in *d,
        CCptrworld *intptr_world, CCptrworld *Rrecord_world),
    grab_record (Rrecord *r, int parity, intptr **list,
        CCptrworld *intptr_world),
    clean_pseudocomb (CCtsp_lpgraph *g, double *x, CCtsp_lpcut_in *c,
        CCtsp_lpcut_in *d, CCptrworld *intptr_world),
    add_to_record (CCptrworld *Rrecord_world, Rrecord **R, int add0, int add1);


CC_PTRWORLD_LIST_ROUTINES (intptr, int, intptralloc, intptr_bulkalloc,
        intptrfree, intptr_listadd, intptr_listfree)
CC_PTRWORLD_LEAKS_ROUTINE (intptr, intptr_check_leaks, this, int)

CC_PTRWORLD_ROUTINES (Rrecord, Rrecordalloc, Rrecord_bulkalloc, Rrecordfree)
CC_PTRWORLD_LISTFREE_ROUTINE (Rrecord, Rrecord_listfree, Rrecordfree)
CC_PTRWORLD_LEAKS_ROUTINE (Rrecord, Rrecord_check_leaks, add0, int)


int CCtsp_teething (CCtsp_lpgraph *g, double *x, CCtsp_lpcut_in *cut,
        CCtsp_lpcut_in **newcut)
{
    int rval = 0;
    int test, ihandle, nbig;
    CCtsp_lpclique **bigteeth = (CCtsp_lpclique **) NULL;
    CCtsp_lpclique *handle = (CCtsp_lpclique *) NULL;

    *newcut = (CCtsp_lpcut_in *) NULL;

    if (cut->cliquecount > 1) {
        rval = CCtsp_test_pure_comb (g->ncount, cut, &test, &ihandle);
        if (rval) {
            fprintf (stderr, "CCtsp_test_pure_comb failed\n"); goto CLEANUP;
        }
        if (test != 1) goto CLEANUP;
    } else {
        ihandle = 0;
    }

    handle = &(cut->cliques[ihandle]);
    bigteeth = CC_SAFE_MALLOC (cut->cliquecount, CCtsp_lpclique *);
    if (!bigteeth) {
        fprintf (stderr, "out of memory in CCtsp_teething\n");
        rval = 1; goto CLEANUP;
    }
    identify_big_teeth (cut, ihandle, &nbig, bigteeth);

    rval = CCtsp_teething_list (g, x, handle, nbig, bigteeth, newcut);
    if (rval) {
        fprintf (stderr, "CCtsp_teething_list failed\n"); goto CLEANUP;
    }

CLEANUP:

    CC_IFFREE (bigteeth, CCtsp_lpclique *);
    return rval;
}

int CCtsp_teething_list (CCtsp_lpgraph *g, double *x, CCtsp_lpclique *handle,
        int nbig, CCtsp_lpclique **bigteeth, CCtsp_lpcut_in **newcut)
{
    int rval = 0;
    CCtsp_lpcut_in general;
    CCtsp_lpcut_in pseudo;
    CCptrworld intptr_world;
    CCptrworld Rrecord_world;
    int test;

    CCptrworld_init (&intptr_world);
    CCptrworld_init (&Rrecord_world);
    
    *newcut = (CCtsp_lpcut_in *) NULL;

    rval = optimal_pseudocomb (g, x, handle, nbig, bigteeth, &pseudo,
                               &intptr_world, &Rrecord_world);
    if (rval) {
        fprintf (stderr, "optimal_pseudocomb failed\n"); goto CLEANUP;
    }

    if (pseudo.cliquecount == 2 || pseudo.cliquecount == 1) {
        CCtsp_free_lpcut_in (&pseudo);
        goto CLEANUP;
    }
    rval = CCtsp_test_pseudocomb (g->ncount, &pseudo, 0, &test);
    if (rval) {
        fprintf (stderr, "CCtsp_test_pseudocomb failed\n"); goto CLEANUP;
    }
    if (test != 1) {
        fprintf (stderr, "Not a pseudocomb\n");
        CCtsp_print_lpcut_in (&pseudo);
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_test_teeth_disjoint (g->ncount, &pseudo, 0, &test);
    if (rval) {
        fprintf (stderr, "CCtsp_test_teeth_disjoint failed\n"); goto CLEANUP;
    }

    if (!test) {
        rval = clean_pseudocomb (g, x, &pseudo, &general, &intptr_world);
        CCtsp_free_lpcut_in (&pseudo);
        if (rval) {
            fprintf (stderr, "clean_pseudocomb failed\n");
            goto CLEANUP;
        }
        if (general.cliquecount <= 2) {
            CCtsp_free_lpcut_in (&general);
            goto CLEANUP;
        }
    }

    *newcut = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
    CCcheck_NULL (*newcut, "out of memory in CCtsp_teething");

    if (test) {
        **newcut = pseudo;
    } else {
        **newcut = general;
    }

    if ((*newcut)->dominocount != 0) {
        printf ("TEETHING Yipes %d\n", (*newcut)->dominocount);
        fflush (stdout);
        exit (1);
    }

    rval = CCtsp_construct_skeleton (*newcut, g->ncount);
    if (rval) {
        fprintf (stderr, "CCtsp_construct_skeleton failed\n");
        CCtsp_free_lpcut_in (*newcut);
        CC_IFFREE (*newcut, CCtsp_lpcut_in);
        goto CLEANUP;
    }

    rval = CCtsp_test_pure_comb (g->ncount, *newcut, &test, (int *) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_test_pure_comb failed\n");
        CCtsp_print_lpcut_in (&general);
        CC_FREE (*newcut, CCtsp_lpcut_in);
        goto CLEANUP;
    }

CLEANUP:

    teething_free_world (&intptr_world, &Rrecord_world);
    return rval;
}

static int optimal_pseudocomb (CCtsp_lpgraph *g, double *x,
        CCtsp_lpclique *handle, int nbig, CCtsp_lpclique **bigteeth,
        CCtsp_lpcut_in *d, CCptrworld *intptr_world, CCptrworld *Rrecord_world)
{
    int rval = 0;
    int *hhit           = (int *) NULL;
    int *thit           = (int *) NULL;
    int *grab0          = (int *) NULL;
    int *grab1          = (int *) NULL;
    int *usebig         = (int *) NULL;
    double *ro0         = (double *) NULL;
    double *ro1         = (double *) NULL;
    intptr *newteeth    = (intptr *) NULL;
    Rrecord **R         = (Rrecord **) NULL;
    int i, j, tmp, sj, u, v, e;
    int ncount = g->ncount;
    int add0, add1, parity;
    double nu0, nu1, delta;
    intptr *ip;
    int ends[2];

    /* Note: bigteeth should be filled in from 1 up to nbig */

    CCtsp_init_lpcut_in (d);

    hhit = CC_SAFE_MALLOC (ncount, int);
    thit = CC_SAFE_MALLOC (ncount, int);
    if (!hhit || !thit) {
        fprintf (stderr, "out of memory in optimal_pseudocombn");
        rval = 1; goto CLEANUP;
    }

    ro0 = CC_SAFE_MALLOC (nbig + 1, double);
    ro1 = CC_SAFE_MALLOC (nbig + 1, double);
    grab0 = CC_SAFE_MALLOC (nbig + 1, int);
    grab1 = CC_SAFE_MALLOC (nbig + 1, int);
    usebig = CC_SAFE_MALLOC (nbig + 1, int);
    if (!ro0 || !ro1 || !grab0 || !grab1 || !usebig) {
        fprintf (stderr, "out of memory in optimal_pseudocomb\n");
        rval = 1; goto CLEANUP;
    }
    R = CC_SAFE_MALLOC (nbig + 1, Rrecord *);
    if (!R) {
        fprintf (stderr, "out of memory in optimal_pseudocomb\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i <= nbig; i++) {
        R[i]  = (Rrecord *) NULL;
        ro0[i] = 0.0;
        ro1[i] = CCtsp_LP_MAXDOUBLE;
    }

    CCtsp_mark_clique_and_neighbors (g, handle, hhit, 0);
    CCtsp_mark_clique_and_neighbors (g, handle, thit, 0);
    for (i = 1; i <= nbig; i++) {
        CCtsp_mark_clique_and_neighbors (g, bigteeth[i], hhit, 0);
        CCtsp_mark_clique_and_neighbors (g, bigteeth[i], thit, 0);
    }

    CCtsp_mark_clique (handle, hhit, 1);
    for (i = 1; i <= nbig; i++) {
        CCtsp_mark_clique (bigteeth[i], thit, i);
    }

    CC_FOREACH_NODE_IN_CLIQUE (u, *handle, tmp) {
        for (sj = 0; sj < g->nodes[u].deg; sj++) {
            v = g->nodes[u].adj[sj].to;
            e = g->nodes[u].adj[sj].edge;
            if (hhit[v] != 1 && x[e] > 0.0) {
                if (thit[v] != 0 && thit[v] == thit[u]) {
                    i = thit[v];
                } else {
                    i = 0;
                }
                nu0 = ro0[i];
                nu1 = ro1[i];
                if (nu1 + (1.0 - 2.0*x[e]) < nu0) {
                    ro0[i] = nu1 + (1.0 - 2.0*x[e]);
                    add0 = e;
                } else {
                    add0 = -1;
                }
                if (nu0 + (1.0 - 2.0*x[e]) < nu1) {
                    ro1[i] = nu0 + (1.0 - 2.0*x[e]);
                    add1 = e;
                } else {
                    add1 = -1;
                }
                if (add0 != -1 || add1 != -1) {
                    rval = add_to_record (Rrecord_world, &R[i], add0, add1);
                    if (rval) {
                        fprintf (stderr, "add_to_record failed\n");
                        goto CLEANUP;
                    }
                }
            }
        }
    }

    for (i = 1; i <= nbig; i++) {
        rval = CCtsp_clique_delta (g, x, bigteeth[i], &delta);
        if (rval) {
            fprintf (stderr, "CCtsp_clique_delta failed\n"); goto CLEANUP;
        }
        if (delta - 3.0 < ro1[i]) {
            ro1[i] = delta - 3.0;
            usebig[i] = 1;
        } else {
            usebig[i] = 0;
        }
        nu0 = ro0[0];
        nu1 = ro1[0];
        if (nu1 + ro1[i] < nu0 + ro0[i]) {
            ro0[0] = nu1 + ro1[i];
            grab0[i] = 1;
        } else {
            ro0[0] = nu0 + ro0[i];
            grab0[i] = 0;
        }
        if (nu0 + ro1[i] < nu1 + ro0[i]) {
            ro1[0] = nu0 + ro1[i];
            grab1[i] = 1;
        } else {
            ro1[0] = nu1 + ro0[i];
            grab1[i] = 0;
        }
    }

    parity = 1;
    for (i = nbig; i > 0; i--) {
        if (parity) {
            if (grab1[i]) {
                if (usebig[i]) {
                    rval = intptr_listadd (&newteeth, -i, intptr_world);
                    if (rval) goto CLEANUP;
                } else {
                    rval = grab_record (R[i], 1, &newteeth, intptr_world);
                    if (rval) goto CLEANUP;
                }
                parity = 0;
            } else {
                rval = grab_record (R[i], 0, &newteeth, intptr_world);
                if (rval) goto CLEANUP;
                parity = 1;
            }
        } else {
            if (grab0[i]) {
                if (usebig) {
                    rval = intptr_listadd (&newteeth, -i, intptr_world);
                    if (rval) goto CLEANUP;
                } else {
                    rval = grab_record (R[i], 1, &newteeth, intptr_world);
                    if (rval) goto CLEANUP;
                }
                parity = 1;
            } else {
                rval = grab_record (R[i], 0, &newteeth, intptr_world);
                if (rval) goto CLEANUP;
                parity = 0;
            }
        }
    }
    rval = grab_record (R[0], parity, &newteeth, intptr_world);
    if (rval) goto CLEANUP;

    d->cliquecount = 1;
    for (ip = newteeth; ip; ip = ip->next) {
        d->cliquecount++;
    }
    d->cliques = CC_SAFE_MALLOC (d->cliquecount, CCtsp_lpclique);
    if (!d->cliques) {
        fprintf (stderr, "out of memory in optimal_pseudocombs\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_copy_lpclique (handle, &d->cliques[0]);
    if (rval) {
        fprintf (stderr, "CCtsp_copy_lpclique failed\n");
        CC_FREE (d->cliques, CCtsp_lpclique);
        goto CLEANUP;
    }
    for (i = 1, ip = newteeth; ip; ip = ip->next, i++) {
        if (ip->this < 0) {
            rval = CCtsp_copy_lpclique (bigteeth[-ip->this], &d->cliques[i]);
            if (rval) {
                fprintf (stderr, "CCtsp_copy_lpclique failed\n");
                for (j = 0; j < i; j++) {
                    CCtsp_free_lpclique (&d->cliques[j]);
                }
                CC_FREE (d->cliques, CCtsp_lpclique);
                goto CLEANUP;
            }
        } else {
            ends[0] = g->edges[ip->this].ends[0];
            ends[1] = g->edges[ip->this].ends[1];
            rval = CCtsp_array_to_lpclique (ends, 2, &d->cliques[i]);
            if (rval) {
                fprintf (stderr, "CCtsp_array_to_lpclique failed\n");
                for (j = 0; j < i; j++) {
                    CCtsp_free_lpclique (&d->cliques[j]);
                }
                CC_FREE (d->cliques, CCtsp_lpclique);
                goto CLEANUP;
            }
        }
    }

    d->sense = 'G';
    d->rhs = 3 * d->cliquecount - 2;

CLEANUP:

    CC_IFFREE (hhit, int);
    CC_IFFREE (thit, int);
    if (R) {
        for (i = 0; i <= nbig; i++) {
            Rrecord_listfree (Rrecord_world, R[i]);
        }
        CC_FREE (R, Rrecord *);
    }
    intptr_listfree (intptr_world, newteeth);
    CC_IFFREE (ro0, double);
    CC_IFFREE (ro1, double);
    CC_IFFREE (grab0, int);
    CC_IFFREE (grab1, int);
    CC_IFFREE (usebig, int);

    return rval;
}

static int grab_record (Rrecord *r, int parity, intptr **list,
        CCptrworld *intptr_world)
{
    int rval = 0;

    while (r) {
        if (parity) {
            if (r->add1 != -1) {
                rval = intptr_listadd (list, r->add1, intptr_world);
                if (rval) goto CLEANUP;
                parity = 0;
            }
        } else {
            if (r->add0 != -1) {
                rval = intptr_listadd (list, r->add0, intptr_world);
                if (rval) goto CLEANUP;
                parity = 1;
            }
        }
        r = r->next;
    }

CLEANUP:

    return rval;
}

static int clean_pseudocomb (CCtsp_lpgraph *g, double *x, CCtsp_lpcut_in *c,
        CCtsp_lpcut_in *d, CCptrworld *intptr_world)
{
    int rval = 0;
    int *hmarks      = (int *) NULL;
    int *activeteeth = (int *) NULL;
    int *cardteeth   = (int *) NULL;
    int *harray      = (int *) NULL;
    intptr **inteeth = (intptr **) NULL;
    intptr *hlist    = (intptr *) NULL;
    int ccount = c->cliquecount;
    CCtsp_lpclique *handle = &c->cliques[0];
    CCtsp_lpclique *cl;
    int i, tmp, k, si, ismarked, cnt, big, ibig, ismall, hcnt;
    double delta, smalldelta;
    intptr *ip;

    CCtsp_init_lpcut_in (d);

    hmarks       = CC_SAFE_MALLOC (g->ncount, int);
    activeteeth  = CC_SAFE_MALLOC (ccount, int);
    cardteeth    = CC_SAFE_MALLOC (ccount, int);
    if (!hmarks || !activeteeth || !cardteeth) {
        fprintf (stderr, "out of memory in clean_pseudocomb\n");
        rval = 1; goto CLEANUP;
    }

    CCtsp_mark_cut (c, hmarks, 0);
    CCtsp_mark_clique (handle, hmarks, 1);
    for (i = 1; i < ccount; i++) {
        activeteeth[i] = 1;
        CCtsp_clique_count (&c->cliques[i], &cardteeth[i]);
    }

    inteeth = CC_SAFE_MALLOC (g->ncount, intptr *);
    if (!inteeth) {
        fprintf (stderr, "out of memory in clean_pseudocomb\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 1; i < ccount; i++) {
        cl = &c->cliques[i];
        CC_FOREACH_NODE_IN_CLIQUE (k, *cl, tmp) {
            inteeth[k] = (intptr *) NULL;
        }
    }
    CC_FOREACH_NODE_IN_CLIQUE (k, *handle, tmp) {
        inteeth[k] = (intptr *) NULL;
    }

    for (i = 1; i < ccount; i++) {
        cl = &c->cliques[i];
        CC_FOREACH_NODE_IN_CLIQUE (k, *cl, tmp) {
            intptr_listadd (&inteeth[k], i, intptr_world);
        }
    }
    for (i = 1; i < ccount; i++) {
        if (!activeteeth[i]) continue;
        cl = &c->cliques[i];
        CC_FOREACH_NODE_IN_CLIQUE (k, *cl, tmp) {
            if (!hmarks[k]) {
                cnt = 0;
                big = 0;
                ibig = -1;
                for (ip = inteeth[k]; ip; ip = ip->next) {
                    if (activeteeth[ip->this]) {
                        cnt++;
                        if (cardteeth[ip->this] > big) {
                            big = cardteeth[ip->this];
                            ibig = ip->this;
                        }
                    }
                }
                if (cnt > 1) {
                    CCtsp_mark_clique (&c->cliques[ibig], hmarks, 1);
                    for (si = 1; si < ccount; si++) {
                        if (activeteeth[si]) {
                            CCtsp_is_clique_marked (&c->cliques[si],
                                                    hmarks, 0, &ismarked);
                            activeteeth[si] = ismarked;
                        }
                    }
                }
            }
        }
    }
    CC_FOREACH_NODE_IN_CLIQUE (k, *handle, tmp) {
        if (hmarks[k]) {
            cnt = 0;
            big = 0;
            ibig = -1;
            for (ip = inteeth[k]; ip; ip = ip->next) {
                if (activeteeth[ip->this]) {
                    cnt++;
                    if (cardteeth[ip->this] > big) {
                        big = cardteeth[ip->this];
                        ibig = ip->this;
                    }
                }
            }
            if (cnt > 1) {
                CCtsp_mark_clique (&c->cliques[ibig], hmarks, 0);
                for (si = 1; si < ccount; si++) {
                    if (activeteeth[si]) {
                        CCtsp_is_clique_marked (&c->cliques[si],
                                                hmarks, 1, &ismarked);
                        activeteeth[si] = ismarked;
                    }
                }
            }
        }
    }

    for (i = 0, hcnt = 0; i < ccount; i++) {
        cl = &c->cliques[i];
        CC_FOREACH_NODE_IN_CLIQUE (k, *cl, tmp) {
            if (hmarks[k] == 1) {
                rval = intptr_listadd (&hlist, k, intptr_world);
                if (rval) goto CLEANUP;
                hmarks[k] = 2;
                hcnt++;
            }
        }
    }

    if (!hcnt) {
        /*
        printf ("WARNING: generalized comb gets an empty handle\n");
        fflush (stdout);
        */
        goto CLEANUP;
    }
    harray = CC_SAFE_MALLOC (hcnt, int);
    if (!harray) {
        fprintf (stderr, "out of memory in clean_pseudocomb\n");
        rval = 1; goto CLEANUP;
    }
    for (ip = hlist, hcnt = 0; ip; ip = ip->next) {
        harray[hcnt++] = ip->this;
    }

    cnt = 0;
    for (i = 1; i < ccount; i++) {
        if (activeteeth[i]) cnt++;
    }
    if (!cnt) {
/*
        printf ("WARNING: generalized comb gets no teeth\n");
        fflush (stdout);
*/
        goto CLEANUP;
    }

    if (cnt % 2 == 0) {
        smalldelta = CCtsp_LP_MAXDOUBLE;
        ismall = -1;
        for (i = 1; i < ccount; i++) {
            if (activeteeth[i]) {
                rval = CCtsp_clique_delta (g, x, &c->cliques[i], &delta);
                if (rval) {
                    fprintf (stderr, "CCtsp_clique_delta failed\n");
                     goto CLEANUP;
                }
                if (delta < smalldelta) {
                    smalldelta = delta;
                    ismall = i;
                }
            }
        }
        activeteeth[ismall] = 0;
        cnt--;
    }
    cnt++;

    d->cliques = CC_SAFE_MALLOC (cnt, CCtsp_lpclique);
    if (!d->cliques) {
        fprintf (stderr, "out of memory in clean_pseudocomb\n");
        rval = 1; goto CLEANUP;
    }
    rval = CCtsp_array_to_lpclique (harray, hcnt, &d->cliques[0]);
    if (rval) {
        fprintf (stderr, "CCtsp_array_to_lpclique failed\n");
        CC_FREE (d->cliques, CCtsp_lpclique);
        goto CLEANUP;
    }
    for (i = 1, cnt = 1; i < ccount; i++) {
        if (activeteeth[i]) {
            rval = CCtsp_copy_lpclique (&c->cliques[i], &d->cliques[cnt++]);
            if (rval) {
                fprintf (stderr, "CCtsp_copy_lpclique failed\n");
                for (k = 0; k < cnt; k++) {
                    CCtsp_free_lpclique (&d->cliques[k]);
                }
                CC_FREE (d->cliques, CCtsp_lpclique);
                goto CLEANUP;
            }
        }
    }
    d->cliquecount = cnt;
    d->rhs = 3 * d->cliquecount - 2;
    d->sense = 'G';

CLEANUP:

    CC_IFFREE (hmarks, int);
    CC_IFFREE (activeteeth, int);
    CC_IFFREE (cardteeth, int);
    CC_IFFREE (harray, int);
    if (inteeth) {
        for (i = 1; i < ccount; i++) {
            cl = &c->cliques[i];
            CC_FOREACH_NODE_IN_CLIQUE (k, *cl, tmp) {
                if (inteeth[k]) {
                    intptr_listfree (intptr_world, inteeth[k]);
                    inteeth[k] = (intptr *) NULL;
                }
            }
        }
        CC_FREE (inteeth, intptr *);
    }
    intptr_listfree (intptr_world, hlist);

    return rval;
}

static void teething_free_world (CCptrworld *intptr_world,
        CCptrworld *Rrecord_world)
{
/*
    int total, onlist;

    if (intptr_check_leaks (intptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding intptrs in teething\n",
                 total - onlist);
    }
*/
    CCptrworld_delete (intptr_world);

/*
    if (Rrecord_check_leaks (Rrecord_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding Rrecords in teething\n",
                 total - onlist);
    }
*/
    CCptrworld_delete (Rrecord_world);
}

static void identify_big_teeth (CCtsp_lpcut_in *cut, int handle, int *nbig,
        CCtsp_lpclique **bigteeth)
{
    int i, k;

    /* teeth filled into positions 1 through nbig */

    *nbig = 0;
    for (i = 0; i < cut->cliquecount; i++) {
        if (i != handle) {
            CCtsp_clique_count (&(cut->cliques[i]), &k);
            if (k >= 3) {
                bigteeth[++(*nbig)] = &(cut->cliques[i]);
            }
        }
    }
}

static int add_to_record (CCptrworld *Rrecord_world, Rrecord **R, int add0,
        int add1)
{
    Rrecord *p;

    p = Rrecordalloc (Rrecord_world);
    if (!p) { fprintf (stderr, "Rrecordalloc failed\n"); return 1; }
    p->add0 = add0;
    p->add1 = add1;
    p->next = *R;
    *R = p;
    return 0;
}



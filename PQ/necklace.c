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
/*  CCpq_necklaces (CCtsp_lpcut_in **cuts, int *cutcount,                   */
/*      CCtsp_cuttree *ctree, int ecount, int *elist, double *x,            */
/*      CCrandstate *rstate)                                                */
/*    NONE                                                                  */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"
#include "pq.h"
#include "cuttree.h"
#include "necklace.h"
#include "verify.h"
#include "macrorus.h"

#undef  DEBUG
#define DEBUGLVL 2
#undef  DUMP_UNREALIZABLE
#undef  ONLY_EXACT

#define NECK_ENUM_CUTOFF 5
#define NECK_ENUM_NTRIES 50
#define NECK_NEXTTRY(x) (((x)*2)/3)

#ifdef ONLY_EXACT
#undef NECK_ENUM_CUTOFF
#define NECK_ENUM_CUTOFF -1
#endif

#define BINSYS_INFEAS 1
#define BINSYS_TRIVIAL 2
#define BINSYS_NONTRIV 3

#define NECK_DUST 0.00001

typedef struct intptr {
    int this;
    struct intptr *next;
} intptr;

typedef struct intptrptr {
    intptr *this;
    struct intptrptr *next;
} intptrptr;

typedef struct neckedge {
    struct necknode *ends[2];
    double x;
    int insystem;
    int inspanning;
    int necklabel;
    int label;
    struct neckedge *next;
} neckedge;

typedef struct neckedgeptr {
    neckedge *this;
    struct necknode *to;
    struct neckedgeptr *next;
} neckedgeptr;

typedef struct necknode {
    neckedgeptr *adj;
    neckedge *entered;
    struct {
        struct necknode *parent;
        int rank;
    } setinfo;
    struct necknode *next;
    intptr *toroot;
    int type;
    struct necknode *child;
    struct necknode *sibling;
    int magiclabel;
    int labeled_children_count;
} necknode;

typedef struct neckgraph {
    int ncount;
    int ecount;
    int neckcount;
    int magicnum;
    necknode *nodelist;
    neckedge *edgelist;
    necknode **necklist;
    necknode *spanroot;
    necknode *cutroot;
    CCptrworld intptr_world;
    CCptrworld intptrptr_world;
    CCptrworld neckedgeptr_world;
    CCptrworld necknode_world;
    CCptrworld eqn_world;
} neckgraph;

typedef struct eqn {
    intptr *lhs;
    int rhs;
    struct eqn *next;
    int pivot;
    int hitdense;
} eqn;

typedef struct bin_var {
#define VALUE_UNKNOWN (-1)
    int value;
    int fixed;
    struct eqn *elim;
} bin_var;

typedef struct bin_system {
    int nvars;
    int nfreevars;
    struct bin_var *vars;
    struct eqn *sparselist;
    struct eqn *denseeqn;
    CCptrworld *intptr_world;
    CCptrworld *eqn_world;
    CCrandstate *rstate;
} bin_system;

CC_PTRWORLD_LIST_ROUTINES (intptr, int, intptr_alloc, intptr_bulk_alloc,
        intptr_free, intptr_listadd, intptr_listfree)
CC_PTRWORLD_LEAKS_ROUTINE (intptr, intptr_check_leaks, this, int)

CC_PTRWORLD_LIST_ROUTINES (intptrptr, intptr *, intptrptr_alloc,
        intptrptr_bulk_alloc, intptrptr_free, intptrptr_listadd,
        intptrptr_listfree)
CC_PTRWORLD_LEAKS_ROUTINE (intptrptr, intptrptr_check_leaks, this, intptr *)

CC_PTRWORLD_ROUTINES (neckedgeptr, neckedgeptr_alloc, neckedgeptr_bulk_alloc,
        neckedgeptr_free)
CC_PTRWORLD_LISTFREE_ROUTINE (neckedgeptr, neckedgeptr_listfree,
        neckedgeptr_free)
CC_PTRWORLD_LEAKS_ROUTINE (neckedgeptr, neckedgeptr_check_leaks, this,
        neckedge *)

CC_PTRWORLD_ROUTINES (necknode, necknode_alloc, necknode_bulk_alloc,
        necknode_free)
CC_PTRWORLD_LISTFREE_ROUTINE (necknode, necknode_listfree, necknode_free)
CC_PTRWORLD_LEAKS_ROUTINE (necknode, necknode_check_leaks, type, int)

CC_PTRWORLD_ROUTINES (eqn, eqn_alloc, eqn_bulk_alloc, eqn_free)
CC_PTRWORLD_LISTFREE_ROUTINE (eqn, eqn_listfree, eqn_free)
CC_PTRWORLD_LEAKS_ROUTINE (eqn, eqn_check_leaks, rhs, int)

static void
    neckgraph_init (neckgraph *g),
    neckgraph_free (neckgraph *g),
    subnecktree_free (necknode *n, CCptrworld *necknode_world),
    free_equation (eqn *sys, bin_system *s),
    intptr_add_destruc (intptr *a, intptr *b, intptr **p_c,
        CCptrworld *intptr_world),
    collect_neck_tooth_leaf (necknode *n, necknode *nodelist,
        int mode, int label, int *arr, int *cnt),
    collect_neck_tooth_work (neckgraph *g, necknode *n, int mode,
        int label, int *arr, int *cnt, necknode *avoid),
    binsys_random_solution (bin_system *s),
    binsys_eval_pivot (bin_system *s, eqn *e),
    binsys_free_system (bin_system *s),
    intptrptr_list_freeall (intptrptr *p, CCptrworld *intptr_world,
        CCptrworld *intptrptr_world),
    ds_makeset (necknode *v);

static int
    necklace_build_graph (neckgraph *g, int ncount, int ecount,
        int *elist, double *x, int *necknum),
    neckgraph_build_adj (neckgraph *g),
    cuttree_to_necktree (CCtsp_cuttree *ctree, neckgraph *g, int neckcount,
        CCtsp_cutnode **necklist),
    necklace_build_spantree (neckgraph *g),
    necklace_crunch_cuts (CCtsp_lpcut_in **cuts, int *cutcount, neckgraph *g,
        CCrandstate *rstate),
    necklace_edge_to_eqn (neckedge *e, eqn **p_eq, bin_system *s),
    necklace_add_edge_to_sys (neckedge *e, bin_system *s, int *status),
    necklace_try_solutions (CCtsp_lpcut_in **cuts, int *cutcount,
        neckgraph *g, bin_system *necksys, intptrptr **p_found),
    compute_toroots (necknode *n, CCptrworld *intptr_world),
    intptr_add (intptr *a, intptr *b, intptr **p_c, CCptrworld *intptr_world),
    intptr_addto (intptr *a, intptr *b, intptr **p_c,
        CCptrworld *intptr_world),
    intptr_copy (intptr *a, intptr **p_b, CCptrworld *intptr_world),
    eqn_addto (eqn *a, eqn *b, bin_system *s),
    intptr_list_size (intptr *p),
    intptrlist_equal (intptr *a, intptr *b),
    find_solution (intptr *sollst, intptrptr *found),
    necklace_checkout_solution (intptr *sollst, neckgraph *g,
        CCtsp_lpcut_in **cuts, int *cutcount),
    count_labeled_children (necknode *n, int label),
    check_realization (necknode *n, int fullcnt, necknode **lefthalf,
        necknode **righthalf),
    collect_necklace_tooth (neckgraph *g, necknode *lefthalf,
        necknode *righthalf, necknode *n, CCtsp_lpclique *c),
    collect_necklace_handle (neckgraph *g, CCtsp_lpclique *c),
    collect_necklace_label (necknode *x, int label),
    binsys_init (bin_system *s, int nvars),
    binsys_add_dense (bin_system *s, eqn *e, int *status),
    binsys_add_sparse (bin_system *s, eqn *e, int *status),
    binsys_elim (bin_system *s, eqn *e),
    binsys_random_minimal_solution (bin_system *s, intptr **sollst),
    binsys_force_zero (bin_system *s, int v, int *status),
    binsys_pop_sparse (bin_system *s),
    binsys_list_solution (bin_system *s, intptr **p_sollst),
    find_label (intptr *e, int label);

static necknode
   *cuttree_to_necktree_work (CCtsp_cutnode *n, CCtsp_cutnode *nodelist,
        necknode *newnodelist, CCgenhash *neckmap, CCptrworld *necknode_world),
   *ds_find (necknode *v),
   *ds_link (necknode *x, necknode *y);

#ifdef DEBUG

static void
    dump_necklaces (neckgraph *g),
    dump_necklace_work (neckgraph *g, necknode *n),
    dump_neckgraph (neckgraph *g);

static int
    find_neck_label (necknode *n, necknode **necklist, int neckcount);

#endif /* DEBUG */


int CCpq_necklaces (CCtsp_lpcut_in **cuts, int *cutcount,
        CCtsp_cuttree *ctree, int ecount, int *elist, double *x,
        CCrandstate *rstate)
{
    int *necknum = (int *) NULL;
    int rval;
    CCtsp_cutnode **necklist = (CCtsp_cutnode **) NULL;
    int neckcount = 0;
    neckgraph ng;
    int ncount = ctree->nodecount;
/*    double szeit = CCutil_zeit ();*/

    neckgraph_init (&ng);

    *cutcount = 0;

    necknum = CC_SAFE_MALLOC (ecount, int);
    if (necknum == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCpq_necklaces\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCpq_cuttree_build_necklaces (ctree, ecount, elist, x, &neckcount,
                                    &necklist, necknum);
    if (rval) {
        fprintf (stderr, "CCpq_cuttree_build_necklaces failed\n");
        goto CLEANUP;
    }

    rval = necklace_build_graph (&ng, ncount, ecount, elist, x, necknum);
    if (rval) {
        fprintf (stderr, "necklace_build_graph failed\n");
        goto CLEANUP;
    }

    rval = cuttree_to_necktree (ctree, &ng, neckcount, necklist);
    if (rval) {
        fprintf (stderr, "cuttree_to_necktree failed\n");
        goto CLEANUP;
    }

#ifdef DEBUG
#if (DEBUGLVL>1)
    printf ("Necklaces:\n");
    dump_necklaces (&ng);
#endif
#endif

    CC_IFFREE (necklist, CCtsp_cutnode *);

    rval = necklace_build_spantree (&ng);
    if (rval) {
        fprintf (stderr, "necklace_build_spantree failed\n");
        goto CLEANUP;
    }

    if (ng.spanroot == (necknode *) NULL) {
        fprintf (stderr, "necklace_build_spantree couldn't build spantree\n");
        rval = 0;
        goto CLEANUP;
    }

    ng.spanroot->toroot = (intptr *) NULL;
    rval = compute_toroots (ng.spanroot, &ng.intptr_world);
    if (rval) {
        fprintf (stderr, "compute_toroots failed\n");
        goto CLEANUP;
    }

#ifdef DEBUG
#if (DEBUGLVL>1)
    dump_neckgraph (&ng);
#endif
#endif

    rval = necklace_crunch_cuts (cuts, cutcount, &ng, rstate);
    if (rval) {
        fprintf (stderr, "necklace_crunch_cuts failed\n");
        goto CLEANUP;
    }

/*
    printf ("CCpq_necklaces finished in %.2f seconds\n",
            CCutil_zeit() - szeit);
    fflush (stdout);
*/

    rval = 0;

  CLEANUP:
    neckgraph_free (&ng);
    CC_IFFREE (necknum, int);
    CC_IFFREE (necklist, CCtsp_cutnode *);

    return rval;
}

static void neckgraph_init (neckgraph *g)
{
    g->ncount = 0;
    g->ecount = 0;
    g->neckcount = 0;
    g->magicnum = 2; /* so we can use 0 and 1 with impunity */
    g->nodelist = (necknode *) NULL;
    g->edgelist = (neckedge *) NULL;
    g->spanroot = (necknode *) NULL;
    g->necklist = (necknode **) NULL;
    g->cutroot = (necknode *) NULL;
    CCptrworld_init (&g->intptr_world);
    CCptrworld_init (&g->intptrptr_world);
    CCptrworld_init (&g->neckedgeptr_world);
    CCptrworld_init (&g->necknode_world);
    CCptrworld_init (&g->eqn_world);
}

static void neckgraph_free (neckgraph *g)
{
    int i;
    necknode *nodelist = g->nodelist;
    int ncount = g->ncount;
    int total;
    int onlist;

    if (g->cutroot) {
        subnecktree_free (g->cutroot, &g->necknode_world);
    }

    if (nodelist) {
        for (i=0; i<ncount; i++) {
            neckedgeptr_listfree (&g->neckedgeptr_world, nodelist[i].adj);
            intptr_listfree (&g->intptr_world, nodelist[i].toroot);
        }
        CC_FREE (g->nodelist, necknode);
    }
    CC_IFFREE (g->edgelist, neckedge);
    CC_IFFREE (g->necklist, necknode *);
    g->ncount = 0;
    g->ecount = 0;
    g->neckcount = 0;

    if (intptr_check_leaks (&g->intptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding NECKLACE-intptr's\n",
                 total - onlist);
    }
    CCptrworld_delete (&g->intptr_world);

    if (intptrptr_check_leaks (&g->intptrptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding NECKLACE-intptrptr's\n",
                 total - onlist);
    }
    CCptrworld_delete (&g->intptrptr_world);
    
    if (neckedgeptr_check_leaks (&g->neckedgeptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding NECKLACE-neckedgeptr's\n",
                 total - onlist);
    }
    CCptrworld_delete (&g->neckedgeptr_world);
    
    if (necknode_check_leaks (&g->necknode_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding NECKLACE-necknode's\n",
                 total - onlist);
    }
    CCptrworld_delete (&g->necknode_world);
    
    if (eqn_check_leaks (&g->eqn_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding NECKLACE-eqn's\n",
                 total - onlist);
    }
    CCptrworld_delete (&g->eqn_world);
}

static int necklace_build_graph (neckgraph *g, int ncount, int ecount,
        int *elist, double *x, int *necknum)
{
    int *perm = (int *) NULL;
    necknode *nodelist = (necknode *) NULL;
    neckedge *edgelist = (neckedge *) NULL;
    int ecount2;
    int i;
    int j;
    int rval;

    neckgraph_free (g);

    perm = CC_SAFE_MALLOC (ecount, int);
    if (perm == (int *) NULL) {
        fprintf (stderr, "Out of memory in necklace_build_graph\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<ecount; i++) perm[i] = i;

    CCutil_double_perm_quicksort (perm, x, ecount);

    nodelist = CC_SAFE_MALLOC (ncount, necknode);
    edgelist = CC_SAFE_MALLOC (ecount, neckedge);
    if (nodelist == (necknode *) NULL ||
        edgelist == (neckedge *) NULL) {
        fprintf (stderr, "Out of memory in necklace_build_graph\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<ncount; i++) {
        nodelist[i].adj = (neckedgeptr *) NULL;
        nodelist[i].entered = (neckedge *) NULL;
        nodelist[i].toroot = (intptr *) NULL;
    }

    ecount2 = 0;
    for (i=0; i<ecount; i++) {
        j = perm[ecount-i-1];
        if (x[j] > NECK_DUST) {
            edgelist[ecount2].ends[0] = &nodelist[elist[2*j]];
            edgelist[ecount2].ends[1] = &nodelist[elist[2*j+1]];
            edgelist[ecount2].x = x[j];
            edgelist[ecount2].insystem = 0;
            edgelist[ecount2].inspanning= 0;
            edgelist[ecount2].necklabel = necknum[j];
            edgelist[ecount2].label = 0;
            edgelist[ecount2].next = (neckedge *) NULL;
            ecount2++;
        }
    }

    g->ncount = ncount;
    g->ecount = ecount2;
    g->nodelist = nodelist;
    g->edgelist = edgelist;

    rval = neckgraph_build_adj (g);
    if (rval) {
        fprintf (stderr, "neckgraph_build_adj failed\n");
        goto CLEANUP;
    }

    rval = 0;

  CLEANUP:
    if (rval) {
        CC_IFFREE (nodelist, necknode);
        CC_IFFREE (edgelist, neckedge);
        g->nodelist = (necknode *) NULL;
        g->edgelist = (neckedge *) NULL;
        g->ncount = 0;
        g->ecount = 0;
    }
    CC_IFFREE (perm, int);
    return rval;
}

static int neckgraph_build_adj (neckgraph *g)
{
    neckedgeptr *p;
    int ecount = g->ecount;
    int ncount = g->ncount;
    neckedge *edgelist = g->edgelist;
    necknode *nodelist = g->nodelist;
    necknode *n0;
    necknode *n1;
    int i;
    int rval;

    for (i=0; i<ncount; i++) {
        nodelist[i].adj = (neckedgeptr *) NULL;
    }

    for (i=0; i<ecount; i++) {
        n0 = edgelist[i].ends[0];
        n1 = edgelist[i].ends[1];
        p = neckedgeptr_alloc (&g->neckedgeptr_world);
        if (p == (neckedgeptr *) NULL) {
            fprintf (stderr, "neckedgeptr_alloc failed\n");
            rval = 1; goto CLEANUP;
        }
        p->this = &edgelist[i];
        p->to = n1;
        p->next = n0->adj;
        n0->adj = p;

        p = neckedgeptr_alloc (&g->neckedgeptr_world);
        if (p == (neckedgeptr *) NULL) {
            fprintf (stderr, "neckedgeptr_alloc failed\n");
            rval = 1; goto CLEANUP;
        }
        p->this = &edgelist[i];
        p->to = n0;
        p->next = n1->adj;
        n1->adj = p;
    }

    rval = 0;

  CLEANUP:
    if (rval) {
        for (i=0; i<ncount; i++) {
            neckedgeptr_listfree (&g->neckedgeptr_world, nodelist[i].adj);
            nodelist[i].adj = (neckedgeptr *) NULL;
        }
    }
    return rval;
}

static int neckcmp (void *key1, void *key2, CC_UNUSED void *u_data)
{
    if (key1 == key2) return 0;
    else return 1;
}

static unsigned int neckhash (void *key1, CC_UNUSED void *u_data)
{
    unsigned long v1 = (unsigned long) key1;
    return (unsigned int) v1;
}

static int cuttree_to_necktree (CCtsp_cuttree *ctree, neckgraph *g,
        int neckcount, CCtsp_cutnode **necklist)
{
    int rval = 0;
    necknode **newnecklist = (necknode **) NULL;
    int ncount = ctree->nodecount;
    int i;
    CCgenhash neckmap;

    rval = CCutil_genhash_init (&neckmap, 2*ncount, neckcmp, neckhash,
                         (void *) NULL, 1.0, 0.5);
    if (rval) {
        fprintf (stderr, "CCutil_genhash_init failed\n");
        return 1;
    }

    g->cutroot = cuttree_to_necktree_work (ctree->root, ctree->nodelist,
                                    g->nodelist, &neckmap, &g->necknode_world);
    if (g->cutroot == (necknode *) NULL) {
        fprintf (stderr, "cuttree_to_necktree failed\n");
        rval = 1; goto CLEANUP;
    }

    newnecklist = CC_SAFE_MALLOC (neckcount, necknode *);
    if (newnecklist == (necknode **) NULL) {
        fprintf (stderr, "Out of memory in cuttree_to_necktree\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<neckcount; i++) {
        newnecklist[i] = (necknode *) CCutil_genhash_lookup (&neckmap,
                                                             necklist[i]);
        if (newnecklist[i] == (necknode *) NULL) {
            fprintf (stderr, "Couldn't map necklist\n");
            rval = 1; goto CLEANUP;
        }
    }

    g->necklist = newnecklist;
    g->neckcount = neckcount;
    rval = 0;

  CLEANUP:
    if (rval) {
        CC_IFFREE (newnecklist,  necknode *);
    }
    CCutil_genhash_free (&neckmap, (void (*) (void *, void *, void *)) NULL);

    return rval;
}

static necknode *cuttree_to_necktree_work (CCtsp_cutnode *n,
        CCtsp_cutnode *nodelist, necknode *newnodelist, CCgenhash *neckmap,
        CCptrworld *necknode_world)
{
    necknode *nn;
    necknode *cn;
    CCtsp_cutnode *c;
    int rval;

    if (n->type == CCtsp_CUT_LEAF || n->type == CCtsp_CUT_EXTERN) {
        nn = &newnodelist[n - nodelist];
        nn->child = (necknode *) NULL;
        nn->sibling = (necknode *) NULL;
        nn->type = n->type;
        return nn;
    }

    nn = necknode_alloc (necknode_world);
    if (nn == (necknode *) NULL) {
        fprintf (stderr, "out of memory in cuttree_to_necktree_work\n");
        goto FAILURE;
    }
    nn->child = (necknode *) NULL;
    nn->sibling = (necknode *) NULL;
    nn->magiclabel = 0;
    nn->type = n->type;

    rval = CCutil_genhash_insert (neckmap, (void *) n, (void *) nn);
    if (rval) {
        fprintf (stderr, "CCutil_genhash_insert failed\n");
        goto FAILURE;
    }

    for (c = n->child; c; c = c->sibling) {
        cn = cuttree_to_necktree_work (c, nodelist, newnodelist, neckmap,
                                       necknode_world);
        if (cn == (necknode *) NULL) goto FAILURE;
        cn->sibling = nn->child;
        nn->child = cn;
    }

    return nn;

  FAILURE:
    if (nn) {
        subnecktree_free (nn, necknode_world);
    }
    return (necknode *) NULL;
}

static void subnecktree_free (necknode *n, CCptrworld *necknode_world)
{
    necknode *c, *cnext;

    for (c = n->child; c; c = cnext) {
        cnext = c->sibling;
        subnecktree_free (c, necknode_world);
    }

    if (n->type != CCtsp_CUT_LEAF && n->type != CCtsp_CUT_EXTERN) {
        necknode_free (necknode_world, n);
    }
}

static int necklace_build_spantree (neckgraph *g)
{
    int ncount = g->ncount;
    int ecount = g->ecount;
    neckedge *edgelist = g->edgelist;
    necknode *nodelist = g->nodelist;
    int i;
    int k;
    neckedge *e;
    neckedgeptr *it;
    necknode *root;
    necknode *front;
    necknode *back;
    necknode *new;
    int cnt;

    for (i=0; i<ncount; i++) {
        ds_makeset (&nodelist[i]);
        nodelist[i].magiclabel = 0;
    }
    for (i=0; i<ecount; i++) {
        edgelist[i].inspanning = 0;
    }

    for (i=0, k=ncount-1; k && i < ecount; i++) {
        if (ds_find (edgelist[i].ends[0]) != ds_find (edgelist[i].ends[1])) {
            ds_link (edgelist[i].ends[0], edgelist[i].ends[1]);
            edgelist[i].inspanning = 1;
            k--;
        }
    }
    if (k) {
        fprintf (stderr, "necklace graph is not connected\n");
        g->spanroot = (necknode *) NULL;
        return 0;
    }

    /* this is a bit clumsy; first we build the spanning tree, then we do
       a breadth-first search through it to root it */

    root = &nodelist[0];
    cnt = 0;

    front = root;
    front->magiclabel = 1;
    cnt++;
    front->entered = (neckedge *) NULL;
    front->next = (necknode *) NULL;
    back = front;

    while (front) {
        for (it = front->adj; it; it = it->next) {
            e = it->this;
            if (e->inspanning) {
                new = it->to;
                if (new->magiclabel == 0) {
                    new->magiclabel = 1;
                    cnt++;
                    new->entered = e;
                    back->next = new;
                    back = new;
                    back->next = (necknode *) NULL;
                }
            }
        }
        front = front->next;
    }
    if (cnt < ncount) {
        fprintf (stderr, "lost the spanning tree\n");
        g->spanroot = (necknode *) NULL;
        return 1;
    }

    g->spanroot = root;
    return 0;
}

static int necklace_crunch_cuts (CCtsp_lpcut_in **cuts, int *cutcount,
        neckgraph *g, CCrandstate *rstate)
{
    bin_system necksys;
    int neckcount = g->neckcount;
    int ecount = g->ecount;
    neckedge *edgelist = g->edgelist;
    eqn *oneseqn = (eqn *) NULL;
    int i;
    int tried;
    intptrptr *found;
    int trynext;
    int status;
    int rval;

    rval = binsys_init (&necksys, neckcount);
    if (rval) {
        fprintf (stderr, "binsys_init failed\n");
        goto CLEANUP;
    }

    necksys.intptr_world = &g->intptr_world;
    necksys.eqn_world = &g->eqn_world;
    necksys.rstate = rstate;
    
    for (i=0; i<ecount; i++) {
        edgelist[i].insystem = edgelist[i].inspanning;
    }

    oneseqn = eqn_alloc (&g->eqn_world);
    if (oneseqn == (eqn *) NULL) {
        fprintf (stderr, "eqn_alloc failed\n");
        rval = 1; goto CLEANUP;
    }

    oneseqn->lhs = (intptr *) NULL;
    for (i = neckcount-1; i>=0; i--) {
        rval = intptr_listadd (&oneseqn->lhs, i, &g->intptr_world);
        if (rval) goto CLEANUP;
    }
    oneseqn->rhs = 1;

    rval = binsys_add_dense (&necksys, oneseqn, &status);
    if (rval) {
        fprintf (stderr, "binsys_add_dense failed\n");
        goto CLEANUP;
    }
    oneseqn = (eqn *) NULL;
    if (status == BINSYS_INFEAS) {
        fprintf (stderr, "ZZZ ODDNESS CONSTRAINT FAILED\n");
        rval = 1; goto CLEANUP;
    }
#ifdef DEBUG
    printf ("1"); fflush (stdout);
#endif

    found = (intptrptr *) NULL;
#ifdef ONLY_EXACT
    trynext = -1;
#else
    trynext = NECK_NEXTTRY (neckcount);
#endif

    tried = 0;
    for (i=0; i<ecount && necksys.nfreevars > NECK_ENUM_CUTOFF; i++) {
        if (!edgelist[i].inspanning) {
            rval = necklace_add_edge_to_sys (&edgelist[i], &necksys, &status);
            if (rval) {
                fprintf (stderr, "necklace_add_edge_to_sys failed\n");
                goto CLEANUP;
            }
            if (status == BINSYS_NONTRIV) {
                edgelist[i].insystem = 1;
#ifdef DEBUG
                printf ("+"); fflush (stdout);
#endif
                tried = 0;
                if (necksys.nfreevars <= trynext) {
#ifdef DEBUG
                    printf (" (%.2f:%d)",edgelist[i].x,necksys.nfreevars);
                    fflush (stdout);
#endif
                    rval = necklace_try_solutions (cuts, cutcount, g,
                                                   &necksys, &found);
                    if (rval) {
                        fprintf (stderr, "necklace_try_solutions failed\n");
                        goto CLEANUP;
                    }
#ifdef DEBUG
                    printf ("\n"); fflush (stdout);
#endif
                    tried = 1;
                    trynext = NECK_NEXTTRY (trynext);
                }
            } else if (status == BINSYS_TRIVIAL) {
                edgelist[i].insystem = 1;
#ifdef DEBUG
                printf ("."); fflush (stdout);
#endif
            } else {
#ifdef DEBUG
                printf ("-"); fflush (stdout);
#endif
#ifdef ONLY_EXACT
                rval = 0; goto CLEANUP;
#endif
            }
        }
    }
    if (!tried) {
#ifdef DEBUG
        printf (" (%.2f:%d)",edgelist[i-1].x,necksys.nfreevars);

        fflush (stdout);
#endif
        rval = necklace_try_solutions (cuts, cutcount, g, &necksys, &found);
        if (rval) {
            fprintf (stderr, "necklace_try_solutions failed\n");
            goto CLEANUP;
        }
#ifdef DEBUG
        printf ("\n"); fflush (stdout);
#endif
    }

    rval = 0;
  CLEANUP:
    free_equation (oneseqn, &necksys);
    binsys_free_system (&necksys);
    intptrptr_list_freeall (found, &g->intptr_world, &g->intptrptr_world);
    return rval;
}

static int necklace_edge_to_eqn (neckedge *e, eqn **p_eq, bin_system *s)
{
    eqn *neweqn;
    intptr *tmp;
    int rval;

    *p_eq = (eqn *) NULL;

    neweqn = eqn_alloc (s->eqn_world);
    if (neweqn == (eqn *) NULL) {
        fprintf (stderr, "eqn_alloc failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = intptr_add (e->ends[0]->toroot, e->ends[1]->toroot, &neweqn->lhs,
                       s->intptr_world);
    if (rval) {
        fprintf (stderr, "intptr_add failed\n");
        goto CLEANUP;
    }

    if (e->necklabel != -1) {
        tmp = intptr_alloc (s->intptr_world);
        if (tmp == (intptr *) NULL) {
            fprintf (stderr, "intptr_alloc failed\n");
            rval = 1; goto CLEANUP;
        }
        tmp->this = e->necklabel;
        tmp->next = (intptr *) NULL;
        intptr_add_destruc (neweqn->lhs, tmp, &neweqn->lhs, s->intptr_world);
    }
    neweqn->rhs = 0;
    rval = 0;
  CLEANUP:
    if (rval) {
        free_equation (neweqn, s);
    } else {
        *p_eq = neweqn;
    }
    return rval;
}

static int necklace_add_edge_to_sys (neckedge *e, bin_system *s, int *status)
{
    eqn *neweqn = (eqn *) NULL;
    int rval;

    rval = necklace_edge_to_eqn (e, &neweqn, s);
    if (rval) {
        fprintf (stderr, "necklace_edge_to_eqn failed\n");
        goto CLEANUP;
    }

    rval = binsys_add_sparse (s, neweqn, status);
    if (rval) {
        fprintf (stderr, "binsys_add_sparse failed\n");
        goto CLEANUP;
    }
    rval = 0;

  CLEANUP:
    return rval;
}

static int necklace_try_solutions (CCtsp_lpcut_in **cuts, int *cutcount,
        neckgraph *g, bin_system *necksys, intptrptr **p_found)
{
    int i;
    intptr *sollst = (intptr *) NULL;
    intptrptr *new = (intptrptr *) NULL;
    int rval;

    for (i=0; i<NECK_ENUM_NTRIES; i++) {
        rval = binsys_random_minimal_solution (necksys, &sollst);
        if (rval) {
            fprintf (stderr, "binsys_random_minimal_solution failed\n");
            goto CLEANUP;
        }
        if (intptr_list_size (sollst) >= 3 &&
            !find_solution (sollst, *p_found)) {
            new = intptrptr_alloc (&g->intptrptr_world);
            if (new == (intptrptr *) NULL) {
                fprintf (stderr, "intptrptr_alloc failed\n");
                rval = 1; goto CLEANUP;
            }
            new->this = sollst;
            sollst = (intptr *) NULL;
            new->next = *p_found;
            *p_found = new;
            rval = necklace_checkout_solution (new->this, g, cuts, cutcount);
            if (rval) {
                fprintf (stderr, "necklace_checkout_solution failed\n");
                goto CLEANUP;
            }
        } else {
            intptr_listfree (&g->intptr_world, sollst);
            sollst = (intptr *) NULL;
        }
    }
    rval = 0;

  CLEANUP:
    intptr_listfree (&g->intptr_world, sollst);
    return rval;
}

static void free_equation (eqn *sys, bin_system *s)
{
    if (sys) {
        intptr_listfree (s->intptr_world, sys->lhs);
        eqn_free (s->eqn_world, sys);
    }
}

static int compute_toroots (necknode *n, CCptrworld *intptr_world)
{
    neckedgeptr *ep;
    neckedge *e;
    necknode *m;
    intptr tmp;
    int rval;

    for (ep = n->adj; ep; ep=ep->next) {
        if (ep->this->inspanning) {
            e = ep->this;
            m = ep->to;
            if (m->entered == e) {
                if (e->necklabel != -1) {
                    tmp.this = e->necklabel;
                    tmp.next = (intptr *) NULL;
                    rval = intptr_add (n->toroot, &tmp, &m->toroot,
                                       intptr_world);
                    if (rval) {
                        fprintf (stderr, "intptr_add failed\n");
                        return rval;
                    }
                } else {
                    rval = intptr_copy (n->toroot, &m->toroot, intptr_world);
                    if (rval) {
                        fprintf (stderr, "intptr_copy failed\n");
                        return rval;
                    }
                }
                rval = compute_toroots (m, intptr_world);
                if (rval) {
                    return rval;
                }
            }
        }
    }
    return 0;
}

static int intptr_add (intptr *a, intptr *b, intptr **p_c,
        CCptrworld *intptr_world)
{
    intptr *sum;
    intptr **sumend;
    int rval;

    sum = (intptr *) NULL;
    sumend = &sum;

    while (a && b) {
        if (a->this == b->this) {
            a = a->next;
            b = b->next;
        } else if (a->this < b->this) {
            *sumend = intptr_alloc (intptr_world);
            if (*sumend == (intptr *) NULL) {
                fprintf (stderr, "intptr_alloc failed\n");
                rval = 1; goto CLEANUP;
            }
            (*sumend)->this = a->this;
            sumend = & (*sumend)->next;
            a = a->next;
        } else {
            assert (a->this > b->this);
            *sumend = intptr_alloc (intptr_world);
            if (*sumend == (intptr *) NULL) {
                fprintf (stderr, "intptr_alloc failed\n");
                rval = 1; goto CLEANUP;
            }
            (*sumend)->this = b->this;
            sumend = & (*sumend)->next;
            b = b->next;
        }
    }
    if (a) {
        rval = intptr_copy (a, sumend, intptr_world);
        if (rval) {
            fprintf (stderr, "intptr_copy failed\n");
            goto CLEANUP;
        }
    } else if (b) {
        rval = intptr_copy (b, sumend, intptr_world);
        if (rval) {
            fprintf (stderr, "intptr_copy failed\n");
            goto CLEANUP;
        }
    } else {
        *sumend = (intptr *) NULL;
    }
    rval = 0;

  CLEANUP:
    if (rval) {
        intptr_listfree (intptr_world, sum);
        *p_c = (intptr *) NULL;
    } else {
        *p_c = sum;
    }
    return rval;
}

static void intptr_add_destruc (intptr *a, intptr *b, intptr **p_c,
        CCptrworld *intptr_world)
{
    intptr **sumend;
    intptr *tnext;

    *p_c = (intptr *) NULL;
    sumend = p_c;

    while (a && b) {
        if (a->this == b->this) {
            tnext = a->next;
            intptr_free (intptr_world, a);
            a = tnext;
            tnext = b->next;
            intptr_free (intptr_world, b);
            b = tnext;
        } else if (a->this < b->this) {
            *sumend = a;
            sumend = &a->next;
            a = a->next;
        } else {
            assert (a->this > b->this);
            *sumend = b;
            sumend = &b->next;
            b = b->next;
        }
    }
    if (a) {
        *sumend = a;
    } else if (b) {
        *sumend = b;
    } else {
        *sumend = (intptr *) NULL;
    }
}

/* destroys a, but leaves b alone */
static int intptr_addto (intptr *a, intptr *b, intptr **p_c,
        CCptrworld *intptr_world)
{
    intptr *sum;
    intptr **sumend;
    intptr *tnext;
    int rval;

    sum = (intptr *) NULL;
    sumend = &sum;

    while (a && b) {
        if (a->this == b->this) {
            tnext = a->next;
            intptr_free (intptr_world, a);
            a = tnext;
            b = b->next;
        } else if (a->this < b->this) {
            *sumend = a;
            sumend = &a->next;
            a = a->next;
        } else {
            assert (a->this > b->this);
            *sumend = intptr_alloc (intptr_world);
            if (*sumend == (intptr *) NULL) {
                fprintf (stderr, "intptr_alloc failed\n");
                rval = 1; goto CLEANUP;
            }
            (*sumend)->this = b->this;
            sumend = & (*sumend)->next;
            b = b->next;
        }
    }
    if (a) {
        *sumend = a;
    } else if (b) {
        rval = intptr_copy (b, sumend, intptr_world);
        if (rval) {
            fprintf (stderr, "intptr_copy failed\n");
            goto CLEANUP;
        }
    } else {
        *sumend = (intptr *) NULL;
    }
    rval = 0;

  CLEANUP:
    if (rval) {
        intptr_listfree (intptr_world, a);
        intptr_listfree (intptr_world, sum);
        *p_c = (intptr *) NULL;
    } else {
        *p_c = sum;
    }
    return rval;
}

static int intptr_copy (intptr *a, intptr **p_b, CCptrworld *intptr_world)
{
    intptr *sum;
    intptr **sumend;
    int rval;

    sum = (intptr *) NULL;
    sumend = &sum;

    while (a) {
        *sumend = intptr_alloc (intptr_world);
        if (*sumend == (intptr *) NULL) {
            fprintf (stderr, "intptr_alloc failed\n");
            rval = 1; goto CLEANUP;
        }
        (*sumend)->this = a->this;
        sumend = & (*sumend)->next;
        a = a->next;
    }
    *sumend = (intptr *) NULL;
    rval = 0;

  CLEANUP:
    if (rval) {
        intptr_listfree (intptr_world, sum);
        *p_b = (intptr *) NULL;
    } else {
        *p_b = sum;
    }
    return rval;
}

static int eqn_addto (eqn *a, eqn *b, bin_system *s)
{
    int rval;

    rval = intptr_addto (a->lhs, b->lhs, &a->lhs, s->intptr_world);
    if (rval) {
        fprintf (stderr, "intptr_addto failed\n");
        eqn_free (s->eqn_world, a);
        return rval;
    }
    a->rhs ^= b->rhs;
    return 0;
}

static int intptr_list_size (intptr *p)
{
    int i;

    for (i=0; p; p = p->next) {
        i++;
    }
    return i;
}

static int intptrlist_equal (intptr *a, intptr *b)
{
    while (a && b) {
        if (a->this != b->this) {
            return 0;
        }
        a = a->next;
        b = b->next;
    }
    return (a == b);
}

static int find_solution (intptr *sollst, intptrptr *found)
{
    while (found) {
        if (intptrlist_equal (sollst, found->this)) {
            return 1;
        }
        found = found->next;
    }
    return 0;
}

static int necklace_checkout_solution (intptr *sollst, neckgraph *g,
                                       CCtsp_lpcut_in **cuts, int *cutcount)
{
    int ncount = g->ncount;
    int ecount = g->ecount;
    neckedge *edgelist = g->edgelist;
    CCtsp_lpcut_in *cut = (CCtsp_lpcut_in *) NULL;
    necknode **necklist = g->necklist;
    int cnt;
    int cnt2;
    intptr *p;
    necknode *lefthalf, *righthalf;
    int i;
    int rval;

#ifdef DEBUG
    printf ("s");
    fflush (stdout);
#endif

    g->magicnum++;

    for (p = sollst, cnt = 0; p; p = p->next) {
        cnt++;
        necklist[p->this]->magiclabel = g->magicnum;
    }

    if (cnt % 2 == 0) {
#ifdef DEBUG
        printf ("E");
        fflush (stdout);
#endif
        rval = 0; goto CLEANUP;
    }

    cnt2 = count_labeled_children (g->cutroot, g->magicnum);

    assert (cnt2 == cnt);

    for (p = sollst; p; p = p->next) {
        if (!check_realization (necklist[p->this], cnt, &lefthalf,
                                &righthalf)) {
#ifdef DEBUG
            printf ("U");
            fflush (stdout);
#endif
#ifdef DUMP_UNREALIZABLE
            printf ("Unrealizable cut, %d total children\n", cnt);
            for (p = sollst; p; p = p->next) {
                printf ("%d(%d): ", p->this,
                        necklist[p->this]->labeled_children_count);
                dump_necklace_work (g, necklist[p->this]);
                printf ("\n");
            }
#endif /* DUMP_UNREALIZABLE */
            rval = 0; goto CLEANUP;
        }
    }
#ifdef DEBUG
    printf ("c");
    fflush (stdout);
#endif

    cut = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
    CCcheck_NULL (cut, "out of memory in necklace_checkout_solution");
    CCtsp_init_lpcut_in (cut);

    cut->cliquecount = cnt + 1;
    cut->cliques = CC_SAFE_MALLOC (cut->cliquecount, CCtsp_lpclique);
    if (cut->cliques == (CCtsp_lpclique *) NULL) {
        fprintf (stderr, "Out of memory in necklace_checkout_solution\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<cut->cliquecount; i++) {
        cut->cliques[i].nodes = (CCtsp_segment *) NULL;
    }

    for (i=0; i<ecount; i++) {
        edgelist[i].label = 0;
    }

    for (p = sollst, cnt2=1; p; p = p->next) {
        if (!check_realization (necklist[p->this], cnt, &lefthalf,
                                &righthalf)) {
            fprintf (stderr, "ZZZ Whoops, necklace broke\n");
            rval = 1; goto CLEANUP;
        }
        rval = collect_necklace_tooth (g, lefthalf, righthalf,
                                       necklist[p->this], &cut->cliques[cnt2]);
        if (rval) {
            fprintf (stderr, "collect_necklace_tooth failed\n");
            goto CLEANUP;
        }
        cnt2++;
    }
    rval = collect_necklace_handle (g, &cut->cliques[0]);
    if (rval) {
        fprintf (stderr, "collect_necklace_handle failed\n");
        goto CLEANUP;
    }
    cut->rhs = CCtsp_COMBRHS(cut);
    cut->next = *cuts;
    cut->sense = 'G';
    cut->branch = 0;

    rval = CCtsp_construct_skeleton (cut, ncount);
    if (rval) {
        fprintf (stderr, "CCtsp_construct_skeleton failed\n");
        goto CLEANUP;
    }

    rval = CCverify_cut (cut, CC_TYPE_COMB, (int *) NULL);
    if (rval == 0) {
        *cuts = cut;
        (*cutcount)++;
        cut = (CCtsp_lpcut_in *) NULL;
    } else {
/*#ifdef DEBUG*/
        printf ("V");
        fflush (stdout);
/*#endif*/
    }

    rval = 0;
  CLEANUP:
    if (cut) {
        CCtsp_free_lpcut_in (cut);
        CC_FREE (cut, CCtsp_lpcut_in);
    }
    return rval;
}

static int count_labeled_children (necknode *n, int label)
{
    int cnt = 0;
    necknode *c;

    for (c = n->child; c; c = c->sibling) {
        cnt += count_labeled_children (c, label);
    }
    if (n->magiclabel == label) {
        cnt++;
    }
    n->labeled_children_count = cnt;

    return cnt;
}

static int check_realization (necknode *n, int fullcnt, necknode **lefthalf,
                              necknode **righthalf)
{
    necknode *c;
    necknode *cnext;

    c = n->child;
    if (c == (necknode *) NULL) return 0;
    cnext = c->sibling;
    while (cnext) {
        if (c->labeled_children_count == 0 &&
            cnext->labeled_children_count == 0) {
            *lefthalf = c;
            *righthalf = cnext;
            return 1;
        }
        c = cnext;
        cnext = cnext->sibling;
    }
    if (n->labeled_children_count == fullcnt) {
        if (c->labeled_children_count == 0) {
            *lefthalf = c;
            *righthalf = (necknode *) NULL;
            return 1;
        }
        if (n->child->labeled_children_count == 0) {
            *lefthalf = n->child;
            *righthalf = (necknode *) NULL;
            return 1;
        }
    }
    return 0;
}

static int collect_necklace_tooth (neckgraph *g, necknode *lefthalf,
        necknode *righthalf, necknode *n, CCtsp_lpclique *c)
{
    int *arr = (int *) NULL;
    int cnt = 0;
    int rval;

    if (lefthalf == (necknode *) NULL) {
        lefthalf = righthalf;
        righthalf = (necknode *) NULL;
    }

    arr = CC_SAFE_MALLOC (g->ncount, int);
    if (arr == (int *) NULL) {
        fprintf (stderr, "out of memory in collect_necklace_tooth\n");
        rval = 1; goto CLEANUP;
    }

    g->magicnum++;

    collect_neck_tooth_work (g, lefthalf, 0, g->magicnum, arr, &cnt,
                             (necknode *) NULL);
    if (righthalf) {
        collect_neck_tooth_work (g, righthalf, 1, g->magicnum, arr, &cnt,
                                 (necknode *) NULL);
    } else {
        collect_neck_tooth_work (g, g->cutroot, 1, g->magicnum, arr, &cnt, n);
    }

    rval = CCtsp_array_to_lpclique (arr, cnt, c);
    if (rval) {
        fprintf (stderr, "CCtsp_array_to_lpclique failed\n");
        goto CLEANUP;
    }
    rval = 0;

  CLEANUP:
    CC_IFFREE (arr, int);
    return rval;
}

static void collect_neck_tooth_leaf (necknode *n, necknode *nodelist,
        int mode, int label, int *arr, int *cnt)
{
    neckedgeptr *ep;
    necknode *n2;

    arr[(*cnt)++] = (int) (n - nodelist);

    if (mode == 0) {
        n->magiclabel = label;
    } else {
        for (ep = n->adj; ep; ep = ep->next) {
            n2 = ep->to;
            if (n2->magiclabel == label) {
                ep->this->label = 1;
            }
        }
    }
}

static void collect_neck_tooth_work (neckgraph *g, necknode *n, int mode,
        int label, int *arr, int *cnt, necknode *avoid)
{
    necknode *c;

    if (n == avoid) return;
    if (n->type == CCtsp_CUT_LEAF || n->type == CCtsp_CUT_EXTERN) {
        collect_neck_tooth_leaf (n, g->nodelist, mode, label, arr, cnt);
    } else {
        for (c = n->child; c; c = c->sibling) {
            collect_neck_tooth_work (g, c, mode, label, arr, cnt, avoid);
        }
    }
}

static int collect_necklace_handle (neckgraph *g, CCtsp_lpclique *c)
{
    int *arr = (int *) NULL;
    int ncount = g->ncount;
    necknode *nodelist = g->nodelist;
    int cnt;
    int i;
    int cnt2;
    int rval;

    arr = CC_SAFE_MALLOC (ncount, int);
    if (arr == (int *) NULL) {
        fprintf (stderr, "out of memory in collect_necklace_handle\n");
        rval = 1; goto CLEANUP;
    }

    g->magicnum++;
    cnt = collect_necklace_label (nodelist, g->magicnum);

/* an ugly fix - don't invert a handle which is the entire graph */
    if (cnt*2 <= ncount || cnt == ncount) {
        cnt2 = 0;
        for (i=0; i<ncount; i++) {
            if (nodelist[i].magiclabel == g->magicnum) {
                arr[cnt2++] = i;
            }
        }
    } else {
        cnt2 = 0;
        for (i=0; i<ncount; i++) {
            if (nodelist[i].magiclabel != g->magicnum) {
                arr[cnt2++] = i;
            }
        }
    }
    rval = CCtsp_array_to_lpclique (arr, cnt2, c);
    if (rval) {
        fprintf (stderr, "CCtsp_array_to_lpclique failed\n");
        goto CLEANUP;
    }
    rval = 0;

  CLEANUP:
    CC_IFFREE (arr, int);
    return rval;
}

static int collect_necklace_label (necknode *x, int label)
{
    int cnt;
    neckedgeptr *ep;
    necknode *y;

    x->magiclabel = label;
    cnt = 1;

    for (ep = x->adj; ep; ep = ep->next) {
        if (ep->this->label == 0 && ep->this->insystem) {
            y = ep->to;
            if (y->magiclabel != label) {
                cnt += collect_necklace_label (y, label);
            }
        }
    }
    return cnt;
}

static int binsys_init (bin_system *s, int nvars)
{
    int i;
    int rval;

    s->vars = (bin_var *) NULL;
    s->sparselist = (eqn *) NULL;
    s->denseeqn = (eqn *) NULL;

    s->vars = CC_SAFE_MALLOC (nvars, bin_var);
    if (s->vars == (bin_var *) NULL) {
        fprintf (stderr, "out of memory in binsys_init\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<nvars; i++) {
        s->vars[i].elim = (eqn *) NULL;
        s->vars[i].value = VALUE_UNKNOWN;
        s->vars[i].fixed = 0;
    }

    s->nvars = nvars;
    s->nfreevars = nvars;
    rval = 0;

  CLEANUP:
    if (rval) {
        CC_IFFREE (s->vars, bin_var);
    }
    return rval;
}

static int binsys_add_dense (bin_system *s, eqn *e, int *status)
{
    int rval;

    assert (s->denseeqn == (eqn *) NULL);

    rval = binsys_elim (s, e);
    if (rval) {
        fprintf (stderr, "binsys_elim failed\n");
        return rval;
    }

    if (e->lhs == (intptr *) NULL) {
        if (e->rhs == 0) {
            eqn_free (s->eqn_world, e);
            *status = BINSYS_TRIVIAL;
            return 0;
        } else {
            eqn_free (s->eqn_world, e);
            *status = BINSYS_INFEAS;
            return 0;
        }
    }
    s->denseeqn = e;
    e->next = (eqn *) NULL;
    s->nfreevars--;
    *status = BINSYS_NONTRIV;
    return 0;
}

static int binsys_add_sparse (bin_system *s, eqn *e, int *status)
{
    int rval;

    rval = binsys_elim (s, e);
    if (rval) {
        fprintf (stderr, "binsys_elim failed\n");
        return rval;
    }

    if (e->lhs == (intptr *) NULL) {
        if (e->rhs == 0) {
            eqn_free (s->eqn_world, e);
            *status = BINSYS_TRIVIAL;
            return 0;
        } else {
            eqn_free (s->eqn_world, e);
            *status = BINSYS_INFEAS;
            return 0;
        }
    }
    if (s->denseeqn && intptrlist_equal (e->lhs, s->denseeqn->lhs)) {
        if (e->rhs == s->denseeqn->rhs) {
            free_equation (e, s);
            *status = BINSYS_TRIVIAL;
            return 0;
        } else {
            free_equation (e, s);
            *status = BINSYS_INFEAS;
            return 0;
        }
    }
    e->pivot = e->lhs->this;
    s->vars[e->pivot].elim = e;
    if (s->denseeqn && find_label (s->denseeqn->lhs, e->pivot)) {
        rval = eqn_addto (s->denseeqn, e, s);
        if (rval) {
            fprintf (stderr, "eqn_addto failed\n");
            s->denseeqn = (eqn *) NULL;
            free_equation (e, s);
            return rval;
        }
        assert (s->denseeqn->lhs);
        e->hitdense = 1;
    } else {
        e->hitdense = 0;
    }
    e->next = s->sparselist;
    s->sparselist = e;
    s->nfreevars--;

    *status = BINSYS_NONTRIV;
    return 0;
}

static int binsys_elim (bin_system *s, eqn *e)
{
    intptr **p;
    intptr *q;
    eqn *f;
    int rval;

    p = &e->lhs;
    while ((q = *p) != (intptr *) NULL) {
        if ((f = s->vars[q->this].elim) != (eqn *) NULL) {
            if (f->lhs && f->lhs->this == q->this) {
                rval = eqn_addto (e, f, s);
                if (rval) {
                    fprintf (stderr, "eqn_addto failed\n");
                    return rval;
                }
            } else {
                rval = eqn_addto (e, f, s);
                if (rval) {
                    fprintf (stderr, "eqn_addto failed\n");
                    return rval;
                }
                p = &e->lhs;
            }
        } else {
            p = & (*p)->next;
        }
    }
    return 0;
}

static int binsys_random_minimal_solution (bin_system *s, intptr **sollst)
{
    int nadded;
    int rval;
    int status;
    int i;

    *sollst = (intptr *) NULL;

    for (i=0; i<s->nvars; i++) {
        s->vars[i].fixed = 0;
    }

    nadded = 0;
    while (s->nfreevars) {
        binsys_random_solution (s);
        for (i=0; i<s->nvars; i++) {
            if (s->vars[i].value == 0 && !s->vars[i].fixed) {
                s->vars[i].fixed = 1;
                rval = binsys_force_zero (s, i, &status);
                if (rval) {
                    fprintf (stderr, "binsys_force_zero failed\n");
                    goto CLEANUP;
                }
                if (status == BINSYS_NONTRIV) {
                    nadded++;
                }
            }
        }
    }
    binsys_random_solution (s);
    rval = binsys_list_solution (s, sollst);
    if (rval) {
        fprintf (stderr, "binsys_list_solution failed\n");
        goto CLEANUP;
    }
    rval = 0;

  CLEANUP:
    for (i=0; i<nadded; i++) {
        if (binsys_pop_sparse (s)) {
            fprintf (stderr, "binsys_pop_sparse failed\n");
            rval = 1;
        }
    }
    return rval;
}

static void binsys_random_solution (bin_system *s)
{
    eqn *e;
    int i;

    for (i=0; i<s->nvars; i++) {
        if (s->vars[i].elim) {
            s->vars[i].value = VALUE_UNKNOWN;
        } else {
            s->vars[i].value = (CCutil_lprand (s->rstate) & 1);
        }
    }
    if (s->denseeqn && s->denseeqn->lhs) {
        s->denseeqn->pivot = s->denseeqn->lhs->this;
        binsys_eval_pivot (s, s->denseeqn);
    }
    for (e = s->sparselist; e; e = e->next) {
        binsys_eval_pivot (s, e);
    }
}

static void binsys_eval_pivot (bin_system *s, eqn *e)
{
    int val;
    intptr *p;

    val = e->rhs;
    for (p = e->lhs; p; p = p->next) {
        if (p->this != e->pivot) {
            assert (s->vars[p->this].value != VALUE_UNKNOWN);
            val ^= s->vars[p->this].value;
        }
    }
    s->vars[e->pivot].value = val;
}

static int binsys_force_zero (bin_system *s, int v, int *status)
{
    eqn *e = (eqn *) NULL;
    intptr *p = (intptr *) NULL;
    int rval;

    e = eqn_alloc (s->eqn_world);
    if (e == (eqn *) NULL) {
        fprintf (stderr, "eqn_alloc failed\n");
        rval = 1; goto CLEANUP;
    }

    e->rhs = 0;
    p = intptr_alloc (s->intptr_world);
    if (p == (intptr *) NULL) {
        fprintf (stderr, "intptr_alloc failed\n");
        rval = 1; goto CLEANUP;
    }

    p->this = v;
    p->next = (intptr *) NULL;
    e->lhs = p;

    rval = binsys_add_sparse (s, e, status);
    return rval;
  CLEANUP:
    if (e) eqn_free (s->eqn_world, e);
    if (p) intptr_free (s->intptr_world, p);
    return rval;
}

static int binsys_pop_sparse (bin_system *s)
{
    eqn *e;
    int rval;

    e = s->sparselist;
    s->sparselist = e->next;
    s->vars[e->pivot].elim = (eqn *) NULL;
    if (e->hitdense && s->denseeqn) {
        rval = eqn_addto (s->denseeqn, e, s);
        if (rval) {
            fprintf (stderr, "eqn_addto failed\n");
            return rval;
        }
    }
    s->nfreevars++;
    free_equation (e, s);
    return 0;
}

static int binsys_list_solution (bin_system *s, intptr **p_sollst)
{
    intptr **lstend;
    intptr *new;
    int rval;
    int i;

    *p_sollst = (intptr *) NULL;
    lstend = p_sollst;
    for (i=0; i<s->nvars; i++) {
        if (s->vars[i].value == 1) {
            new = intptr_alloc (s->intptr_world);
            if (new == (intptr *) NULL) {
                fprintf (stderr, "intptr_alloc failed\n");
                rval = 1; goto CLEANUP;
            }
            new->this = i;
            *lstend = new;
            lstend = &new->next;
        }
    }
    rval = 0;
  CLEANUP:
    *lstend = (intptr *) NULL;
    if (rval) {
        intptr_listfree (s->intptr_world, *p_sollst);
        *p_sollst = (intptr *) NULL;
    }
    return rval;
}

static void binsys_free_system (bin_system *s)
{
    eqn *e, *enext;

    CC_FREE (s->vars, bin_var);

    for (e = s->sparselist; e; e = enext) {
        enext = e->next;
        free_equation (e, s);
    }
    if (s->denseeqn) {
        free_equation (s->denseeqn, s);
    }

    s->vars = (bin_var *) NULL;
    s->sparselist = (eqn *) NULL;
    s->denseeqn = (eqn *) NULL;
    s->nvars = 0;
    s->nfreevars = 0;
}

#ifdef DEBUG
static void dump_necklaces (neckgraph *g)
{
    int i;

    for (i=0; i<g->neckcount; i++) {
        printf ("%d: ", i);
        dump_necklace_work (g, g->necklist[i]);
        printf ("\n");
    }
}

static void dump_necklace_work (neckgraph *g, necknode *n)
{
    necknode *c;
    int j;

    for (c = n->child; c; c = c->sibling) {
        if (c->type == CCtsp_CUT_QNODE || (c->type == CCtsp_CUT_PNODE &&
                c->child && c->child->sibling &&
                !c->child->sibling->sibling)) {
            j = find_neck_label (c, g->necklist, g->neckcount);
            printf ("%d ", j);
        } else if (c->type == CCtsp_CUT_LEAF) {
            printf ("n%d ", c - g->nodelist);
        } else {
            printf ("(");
            dump_necklace_work (g, c);
            printf (") ");
        }
    }
}

static void dump_neckgraph (neckgraph *g)
{
    int i;
    int ncount = g->ncount;
    necknode *nodelist = g->nodelist;
    int ecount = g->ecount;
    neckedge *edgelist = g->edgelist;
    intptr *p;

    printf ("Neck nodes:\n");
    for (i=0; i<ncount; i++) {
        printf ("%d:", i);
        for (p = nodelist[i].toroot; p; p = p->next) {
            printf (" %d", p->this);
        }
        printf ("\n");
    }
    printf ("Neck edges:\n");
    for (i=0; i<ecount; i++) {
        printf ("%d (%d %d) %.2f:", i, edgelist[i].ends[0] - nodelist,
                edgelist[i].ends[1] - nodelist, edgelist[i].x);
        if (edgelist[i].insystem) printf (" sys");
        if (edgelist[i].inspanning) printf (" span");
        if (edgelist[i].necklabel != -1)
            printf (" neck %d", edgelist[i].necklabel);
        printf ("\n");
    }
    fflush (stdout);
}

static int find_neck_label (necknode *n, necknode **necklist, int neckcount)
{
    int i;

    for (i=0; i<neckcount; i++) {
        if (necklist[i] == n) return i;
    }
    return -1;
}
#endif /* DEBUG */

static int find_label (intptr *e, int label)
{
    while (e) {
        if (e->this == label) return 1;
        if (e->this > label) return 0;
        e = e->next;
    }
    return 0;
}

static void intptrptr_list_freeall (intptrptr *p, CCptrworld *intptr_world,
        CCptrworld *intptrptr_world)
{
    intptrptr *pnext;

    while (p) {
        pnext = p->next;
        intptr_listfree (intptr_world, p->this);
        intptrptr_free (intptrptr_world, p);
        p = pnext;
    }
}

/* disjoint sets ala Tarjan (from Data Structures and Network Algorithms
   by Robert Tarjan) */

static void ds_makeset (necknode *v)
{
    v->setinfo.parent = v;
    v->setinfo.rank = 0;
}

static necknode *ds_find (necknode *v)
{
    necknode *p = v->setinfo.parent;

    return v == p ? v : (v->setinfo.parent = ds_find (p));
}

static necknode *ds_link (necknode *x, necknode *y)
{
    necknode *t;

    x = ds_find (x);
    y = ds_find (y);

    if (x != y) {
        if (x->setinfo.rank > y->setinfo.rank) {
            CC_SWAP (x,y,t);
        } else if (x->setinfo.rank == y->setinfo.rank) {
            y->setinfo.rank++;
        }
        x->setinfo.parent = y;
    }
    return y;
}

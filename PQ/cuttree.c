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
/*  void CCpq_cuttree_init (CCtsp_cuttree *t)                               */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCpq_cuttree_freetree (CCtsp_cuttree *t)                           */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCpq_check_clique (CCpq_tree *pqt, CCtsp_lpclique *c,              */
/*      int *status)                                                        */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCpq_cuttree_display (CCtsp_cuttree *t)                            */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCpq_cuttree_describe (CCtsp_cuttree *t)                           */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCpq_cuttree_trivial (CCtsp_cuttree *t, int nodecount,              */
/*      int extern_node)                                                    */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCpq_cuttree_update_clean (CCtsp_cuttree *t, int edgecount,         */
/*      int *elist, double *x)                                              */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCpq_cuttree_improve_quick (CCtsp_cuttree *t, CCtsp_lpcuts *pool,   */
/*      int edgecount, int *elist, double *x)                               */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCpq_apply_clique (CCpq_tree *T, CCtsp_lpclique *c, int *status)    */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCpq_cuttree_gen_cliques (CCtsp_cuttree *t, void *u_data,           */
/*      int (*cut_callback) (int *arr, int cnt, int *stop, void *u_data))   */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCpq_cuttree_build_necklaces (CCtsp_cuttree *t, int ecount,         */
/*      int *elist, double *x, int *p_neckcount,                            */
/*      CCtsp_cutnode ***p_necklist, int *necknum)                          */
/*    MISSING                                                               */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"
#include "pq.h"
#include "cuttree.h"
#include "cut.h"

#undef  DEBUG
#undef  VERBOSE

#define CUT_TOL 0.001

typedef struct cutgradj {
    double weight;
    struct cutgrnode *to;
    struct cutgrnode *from;
    struct cutgradj *next;
    int num;
} cutgradj;

typedef struct cutgrnode {
    int type;
    struct cutgrnode *parent;
    struct cutgrnode *child;
    struct cutgrnode *sibling;
    double int_weight;
    double ext_weight;
    int subtree_size;
    cutgradj *adj;
    cutgradj *adjinto;
    struct {
        int rank;
        struct cutgrnode *parent;
    } disj_set;
    struct cutgrnode *ancestor;
    struct cutgrnode *next;
    CCpq_node *leaf_elems;
    CCpq_node **leaf_elems_end;
    CCpq_node *save_leaf_next;
    struct CCtsp_cutnode *cn_base;
    int mark;
} cutgrnode;

typedef struct cutgrtree {
    int nodecount;
    int extern_node;
    cutgrnode *nodelist;
    cutgrnode *root;
    CCptrworld cutgrnode_world;
    CCptrworld cutgradj_world;
} cutgrtree;

struct pq_callback_struct {
    CCpq_tree *pqt;
    int nontrivial;
    int cutcount;
};


static void
    subcuttree_free (CCtsp_cutnode *n, CCptrworld *cutnode_world),
    cutgrnode_init (cutgrnode *n),
    cutgrtree_init (cutgrtree *t),
    cutgrtree_freetree (cutgrtree *t),
    subcutgrtree_free (cutgrnode *n, CCptrworld *cutgrnode_world,
        CCptrworld *cutgradj_world),
    cutgrtree_free_leafgraph (cutgrtree *t),
    cutgrtree_leafgraph_init (cutgrnode *n),
    cutgrtree_cleanadj (cutgrnode *n, CCptrworld *cutgradj_world),
    build_leaf_elems (cutgrnode *n, CCpq_node *elems, cutgrnode *nodelist),
    make_set (cutgrnode *x),
    label_necklaces (cutgrnode *n, CCtsp_cutnode **necklist, int *cnt,
        int *necknum),
    cuttree_display_work (CCtsp_cutnode *n, CCtsp_cutnode *nodelist),
    cuttree_describe_work (CCtsp_cutnode *n);

static int
    PQ_tree_to_cuttree (CCpq_tree *pqt, CCtsp_cuttree *ct),
    cuttree_to_cutgrtree (CCtsp_cuttree *t, cutgrtree *tgr),
    PQ_tree_load_segments (CCpq_tree *pqt, double cut_tol, int edgecount,
                           int *elist, double *x),
    load_seg_callback (double cut_val, int cut_start, int cut_end,
                       void *callbackdat),
    PQ_tree_load_pool (CCpq_tree *pqt,  double cut_tol, CCtsp_lpcuts *pool,
                       int edgecount, int *elist, double *x),
    cutgrtree_build_leafgraph (cutgrtree *t, int edgecount, int *elist,
                               double *x),
    cutgrtree_addadj (cutgrnode *n, cutgrnode *from, cutgrnode *to,
        double weight, int num, CCptrworld *cutgradj_world),
    cutgrtree_loadx (cutgrtree *t, int edgecount, int *elist, double *x,
                     int cleanup),
    cutgrtree_ancestor (cutgrnode *n, CCptrworld *cutgradj_world),
    cutgrtree_penultimate (cutgrnode *n, CCptrworld *cutgradj_world),
    PQ_add_tight_nodes (double tol, CCpq_tree *pqt, cutgrtree *t),
    add_tight_nodes_work (cutgrnode *n, double tol, CCpq_tree *pqt),
    add_node (cutgrnode *n, CCpq_tree *pqt),
    add_node_pair (cutgrnode *n1, cutgrnode *n2, CCpq_tree *pqt),
    add_node_complement (cutgrnode *p, cutgrnode *n, CCpq_tree *pqt),
    PQ_add_tight_edges (double tol, CCpq_tree *pqt, cutgrtree *t),
    add_tight_edges_work (cutgrnode *n, double tol, CCpq_tree *pqt),
    check_tight_path (cutgrnode *n, double tol, CCpq_tree *pqt),
    add_tight_path (cutgrnode *end1, cutgrnode *end2,
                    cutgrnode *avoid, double tol, CCpq_tree *pqt),
    check_tight_ext (cutgrnode *n, double tol, CCpq_tree *pqt),
    mark_tree (cutgrnode *n, int m),
    cuttree_gen_work (CCtsp_cutnode *n, CCtsp_cutnode *nodelist, int *nodenums,
        int *size, int *stop, void *u_data,
        int (*cut_callback) (int *arr, int cnt, int *stop, void *u_data)),
    count_necklaces (cutgrnode *n);

static CCtsp_cutnode
   *pqtree_to_cuttree_work (CCpq_node *x, CCpq_node *elems,
        CCtsp_cutnode *nodelist, CCptrworld *cutnode_world);

static cutgrnode
   *cuttree_to_cutgrtree_work (CCtsp_cutnode *n, CCtsp_cutnode *nodelist,
        cutgrnode *grnodelist, CCptrworld *cutgrnode_world,
        CCptrworld *cutgradj_world),
   *find_ext_node (cutgrnode *n, cutgrnode *avoid, double tol),
   *find_set (cutgrnode *x),
   *link_set (cutgrnode *x, cutgrnode *y);

static cutgradj
   *find_edge (cutgrnode *end1, cutgrnode *end2),
   *find_path_edge (cutgrnode *n, cutgrnode *avoid, double tol);


#ifdef DEBUG
static void
    verify_cutgrtree (cutgrtree *t, cutgrnode *n, int edgecount,
                      int *elist, double *x);
#endif


CC_PTRWORLD_ROUTINES (CCtsp_cutnode, cutnode_alloc, cutnode_bulk_alloc,
        cutnode_free)
CC_PTRWORLD_LEAKS_ROUTINE (CCtsp_cutnode, cutnode_check_leaks, type, int)

CC_PTRWORLD_ROUTINES (cutgradj, cutgradj_alloc, cutgradj_bulk_alloc,
        cutgradj_free)
CC_PTRWORLD_LISTFREE_ROUTINE (cutgradj, cutgradj_listfree, cutgradj_free)
CC_PTRWORLD_LEAKS_ROUTINE (cutgradj, cutgradj_check_leaks, num, int)

CC_PTRWORLD_ROUTINES (cutgrnode, cutgrnode_alloc, cutgrnode_bulk_alloc,
        cutgrnode_free)
CC_PTRWORLD_LEAKS_ROUTINE (cutgrnode, cutgrnode_check_leaks, type, int)

void CCpq_cuttree_init (CCtsp_cuttree *t)
{
    t->root = (CCtsp_cutnode *) NULL;
    t->nodelist = (CCtsp_cutnode *) NULL;
    t->nodecount = -1;
    t->extern_node = -1;
    CCptrworld_init (&t->cutnode_world);
}

void CCpq_cuttree_freetree (CCtsp_cuttree *t)
{
    int total;
    int onlist;

    if (t->root) {
        subcuttree_free (t->root, &t->cutnode_world);
        t->root = (CCtsp_cutnode *) NULL;
    }
    CC_IFFREE (t->nodelist, CCtsp_cutnode);

    if (cutnode_check_leaks (&t->cutnode_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding cutnodes\n",
                 total - onlist);
    }
    CCptrworld_delete (&t->cutnode_world);
}

static void subcuttree_free (CCtsp_cutnode *n, CCptrworld *cutnode_world)
{
    CCtsp_cutnode *p;
    CCtsp_cutnode *pnext;

    if (CCtsp_CUT_INNODELIST (n->type)) {
        return;
    }

    for (p = n->child; p; p = pnext) {
        pnext = p->sibling;
        subcuttree_free (p, cutnode_world);
    }
    cutnode_free (cutnode_world, n);
}

int CCpq_cuttree_trivial (CCtsp_cuttree *t, int nodecount, int extern_node)
{
    CCtsp_cutnode *intern = (CCtsp_cutnode *) NULL;
    CCtsp_cutnode *root = (CCtsp_cutnode *) NULL;
    CCtsp_cutnode *nodelist = (CCtsp_cutnode *) NULL;
    int i;

    CCpq_cuttree_freetree (t);
    CCpq_cuttree_init (t);

    if (nodecount < 3) {
        fprintf (stderr, "CCpq_cuttree_trivial can only handle graphs with >= 3 nodes\n");
        goto FAILURE;
    }

    nodelist = CC_SAFE_MALLOC (nodecount, CCtsp_cutnode);
    if (nodelist == (CCtsp_cutnode *) NULL) {
        fprintf (stderr, "Out of memory in CCpq_cuttree_trivial\n");
        goto FAILURE;
    }

    root = cutnode_alloc (&t->cutnode_world);
    intern = cutnode_alloc (&t->cutnode_world);
    if (root == (CCtsp_cutnode *) NULL || intern == (CCtsp_cutnode *) NULL) {
        fprintf (stderr, "cutnode_alloc failed\n");
        goto FAILURE;
    }

    root->type = CCtsp_CUT_ROOT;
    root->sibling = (CCtsp_cutnode *) NULL;
    root->child = &nodelist[extern_node];

    nodelist[extern_node].type = CCtsp_CUT_EXTERN;
    nodelist[extern_node].sibling = intern;
    nodelist[extern_node].child = (CCtsp_cutnode *) NULL;

    intern->type = CCtsp_CUT_PNODE;
    intern->sibling = (CCtsp_cutnode *) NULL;
    intern->child = (CCtsp_cutnode *) NULL;

    for (i=0; i<nodecount; i++) {
        if (i != extern_node) {
            nodelist[i].type = CCtsp_CUT_LEAF;
            nodelist[i].sibling = intern->child;
            intern->child = &nodelist[i];
            nodelist[i].child = (CCtsp_cutnode *) NULL;
        }
    }

    t->nodelist = nodelist;
    t->root = root;
    t->nodecount = nodecount;
    t->extern_node = extern_node;
    return 0;

  FAILURE:
    if (intern) cutnode_free (&t->cutnode_world, intern);
    if (root) cutnode_free (&t->cutnode_world, root);
    CC_IFFREE (nodelist, CCtsp_cutnode);
    CCpq_cuttree_freetree (t);
    return -1;
}

int CCpq_cuttree_update_clean (CCtsp_cuttree *t, int edgecount, int *elist,
        double *x)
{
    CCtsp_cuttree tnew;
    cutgrtree tgr;
    CCpq_tree pqt;
    int rval;

    cutgrtree_init (&tgr);
    CCpq_tree_init (&pqt);
    CCpq_cuttree_init (&tnew);

    rval = CCpq_tree_trivial (&pqt, t->nodecount, t->extern_node);
    if (rval) {
        fprintf (stderr, "CCpq_tree_trivial failed\n");
        goto CLEANUP;
    }

    rval = cuttree_to_cutgrtree (t, &tgr);
    if (rval) {
        fprintf (stderr, "cuttree_to_cutgrtree failed\n");
        goto CLEANUP;
    }

    rval = cutgrtree_loadx (&tgr, edgecount, elist, x, 1);
    if (rval) {
        fprintf (stderr, "cutgrtree_loadx failed\n");
        goto CLEANUP;
    }

    rval = PQ_add_tight_nodes (CUT_TOL, &pqt, &tgr);
    if (rval) {
        fprintf (stderr, "PQ_add_tight_nodes failed\n");
        goto CLEANUP;
    }

    cutgrtree_freetree (&tgr);

    rval = PQ_tree_to_cuttree (&pqt, &tnew);
    if (rval) {
        fprintf (stderr, "PQ_tree_to_cuttree failed\n");
        goto CLEANUP;
    }

    rval = cuttree_to_cutgrtree (&tnew, &tgr);
    if (rval) {
        fprintf (stderr, "cuttree_to_cutgrtree failed\n");
        goto CLEANUP;
    }

    CCpq_cuttree_freetree (&tnew);

    rval = cutgrtree_loadx (&tgr, edgecount, elist, x, 1);
    if (rval) {
        fprintf (stderr, "cutgrtree_loadx failed\n");
        goto CLEANUP;
    }

    rval = PQ_add_tight_edges (CUT_TOL, &pqt, &tgr);
    if (rval) {
        fprintf (stderr, "PQ_add_tight_edges failed\n");
        goto CLEANUP;
    }

    cutgrtree_freetree (&tgr);

    rval = PQ_tree_to_cuttree (&pqt, &tnew);
    if (rval) {
        fprintf (stderr, "PQ_tree_to_cuttree failed\n");
        goto CLEANUP;
    }

    CCpq_cuttree_freetree (t);
    *t = tnew;
    CCpq_cuttree_init (&tnew);

    rval = 0;
  CLEANUP:
    CCpq_cuttree_freetree (&tnew);
    cutgrtree_freetree (&tgr);
    CCpq_tree_free (&pqt);
    return rval;
}

int CCpq_cuttree_improve_quick (CCtsp_cuttree *t, CCtsp_lpcuts *pool,
        int edgecount, int *elist, double *x)
{
    CCtsp_cuttree tnew;
    CCpq_tree pqt;
    int rval;
/*
    double szeit = CCutil_zeit();
*/

    CCpq_cuttree_init (&tnew);
    CCpq_tree_init (&pqt);

    rval = CCpq_cuttree_update_clean (t, edgecount, elist, x);
    if (rval) {
        fprintf (stderr, "CCpq_cuttree_update_clean failed\n");
        goto CLEANUP;
    }

    rval = CCpq_cuttree_to_pq (t, &pqt);
    if (rval) {
        fprintf (stderr, "CCpq_cuttree_to_pq failed\n");
        goto CLEANUP;
    }

    rval = PQ_tree_load_segments (&pqt, CUT_TOL, edgecount, elist, x);
    if (rval) {
        fprintf (stderr, "PQ_tree_load_segments failed\n");
        goto CLEANUP;
    }

    rval = PQ_tree_load_pool (&pqt, CUT_TOL, pool, edgecount, elist, x);
    if (rval) {
        fprintf (stderr, "PQ_tree_load_pool failed\n");
        goto CLEANUP;
    }

    rval = PQ_tree_to_cuttree (&pqt, &tnew);
    if (rval) {
        fprintf (stderr, "PQ_tree_to_cuttree failed\n");
        goto CLEANUP;
    }

    rval = CCpq_cuttree_update_clean (&tnew, edgecount, elist, x);
    if (rval) {
        fprintf (stderr, "CCpq_cuttree_update_clean failed\n");
        goto CLEANUP;
    }

    CCpq_cuttree_freetree (t);
    *t = tnew;
    CCpq_cuttree_init (&tnew);

/*
    printf ("CCtsp_cuttree improved (quick) in %.2f seconds\n",
            CCutil_zeit() - szeit);
    fflush (stdout);
*/

    rval = 0;
  CLEANUP:
    CCpq_cuttree_freetree (&tnew);
    CCpq_tree_free (&pqt);
    return rval;
}

static int PQ_tree_to_cuttree (CCpq_tree *pqt, CCtsp_cuttree *ct)
{
    int i;
    CCtsp_cutnode *ext = (CCtsp_cutnode *) NULL;
    CCtsp_cutnode *newroot = (CCtsp_cutnode *) NULL;
    int nodecount = pqt->nodecount;
    int rval;

    ct->nodelist = CC_SAFE_MALLOC (nodecount, CCtsp_cutnode);
    if (ct->nodelist == (CCtsp_cutnode *) NULL) {
        fprintf (stderr, "Out of memory in PQ_tree_to_cuttree\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<nodecount; i++) {
        ct->nodelist[i].type = CCtsp_CUT_EXTERN;
    }

    ct->root = pqtree_to_cuttree_work (CCpq_find_root (pqt), pqt->elems,
                                       ct->nodelist, &ct->cutnode_world);
    if (ct->root == (CCtsp_cutnode *) NULL) {
        fprintf (stderr, "pqtree_to_cuttree_work failed\n");
        rval = 1; goto CLEANUP;
    }

    ct->extern_node = pqt->extern_node;

    ext=&ct->nodelist[ct->extern_node];

    newroot = cutnode_alloc (&ct->cutnode_world);
    if (newroot == (CCtsp_cutnode *) NULL) {
        fprintf (stderr, "Out of memory in PQ_tree_to_cuttree\n");
        rval = 1; goto CLEANUP;
    }

    newroot->type = CCtsp_CUT_ROOT;
    newroot->sibling = (CCtsp_cutnode *) NULL;
    newroot->child = ext;

    ext->type = CCtsp_CUT_EXTERN;
    ext->sibling = ct->root;
    ext->child = (CCtsp_cutnode *) NULL;

    ct->root = newroot;
    ct->nodecount = nodecount;
    newroot = (CCtsp_cutnode *) NULL;

    rval = 0;

  CLEANUP:
    if (rval) {
        if (newroot) cutnode_free (&ct->cutnode_world, newroot);
        CCpq_cuttree_freetree (ct);
    }
    return 0;
}

static CCtsp_cutnode *pqtree_to_cuttree_work (CCpq_node *x, CCpq_node *elems,
        CCtsp_cutnode *nodelist, CCptrworld *cutnode_world)
{
    CCtsp_cutnode *n = (CCtsp_cutnode *) NULL;
    CCtsp_cutnode *c;
    CCpq_node *p;
    CCpq_node *pnext;
    CCpq_node *pprev;

    if (x->type == PQ_LEAF) {
        n = &nodelist[x - elems];
        n->type = CCtsp_CUT_LEAF;
        n->child = (CCtsp_cutnode *) NULL;
        n->sibling = (CCtsp_cutnode *) NULL;
        return n;
    }

    n = cutnode_alloc (cutnode_world);
    if (n == (CCtsp_cutnode *) NULL) {
        fprintf (stderr, "Out of memory in pqtree_to_cuttree_work\n");
        goto FAILURE;
    }

    n->child = (CCtsp_cutnode *) NULL;
    n->sibling = (CCtsp_cutnode *) NULL;

    if (x->type == PQ_PNODE) {
        n->type = CCtsp_CUT_PNODE;
    } else if (x->type == PQ_QNODE) {
        n->type = CCtsp_CUT_QNODE;
    } else {
        fprintf (stderr, "Unknown node type %d\n", x->type);
        goto FAILURE;
    }

    CCpq_set_FOREACH (x->children_set, p, children_elem, pprev, pnext) {
        c = pqtree_to_cuttree_work (p, elems, nodelist, cutnode_world);
        if (c == (CCtsp_cutnode *) NULL) goto FAILURE;
        c->sibling = n->child;
        n->child = c;
    }

    return n;

  FAILURE:
    if (n) {
        subcuttree_free (n, cutnode_world);
    }
    return (CCtsp_cutnode *) NULL;
}

static int cuttree_to_cutgrtree (CCtsp_cuttree *t, cutgrtree *tgr)
{
    int nodecount = t->nodecount;

    tgr->nodelist = CC_SAFE_MALLOC (nodecount, cutgrnode);
    if (tgr->nodelist == (cutgrnode *) NULL) {
        fprintf (stderr, "Out of memory in cuttree_to_cutgrtree\n");
        return -1;
    }

    tgr->root = cuttree_to_cutgrtree_work (t->root, t->nodelist,
                  tgr->nodelist, &tgr->cutgrnode_world, &tgr->cutgradj_world);
    if (tgr->root == (cutgrnode *) NULL) {
        fprintf (stderr, "cuttree_to_cutgrtree_work failed\n");
        cutgrtree_freetree (tgr);
        return -1;
    }
    tgr->nodecount = t->nodecount;
    tgr->extern_node = t->extern_node;

    return 0;
}

static cutgrnode *cuttree_to_cutgrtree_work (CCtsp_cutnode *n,
        CCtsp_cutnode *nodelist, cutgrnode *grnodelist,
        CCptrworld *cutgrnode_world, CCptrworld *cutgradj_world)
{
    cutgrnode *ngr;
    cutgrnode *cgr;
    CCtsp_cutnode *c;

    if (n->type == CCtsp_CUT_LEAF || n->type == CCtsp_CUT_EXTERN) {
        ngr = &grnodelist[n - nodelist];
        cutgrnode_init (ngr);
        ngr->type = n->type;
        ngr->cn_base = n;
        return ngr;
    }

    ngr = cutgrnode_alloc (cutgrnode_world);
    if (ngr == (cutgrnode *) NULL) {
        fprintf (stderr, "Out of memory in cuttree_cutgrtree_work\n");
        goto FAILURE;
    }
    cutgrnode_init (ngr);

    ngr->type = n->type;
    ngr->cn_base = n;

    for (c = n->child; c; c = c->sibling) {
        cgr = cuttree_to_cutgrtree_work (c, nodelist, grnodelist,
                                         cutgrnode_world, cutgradj_world);
        if (cgr == (cutgrnode *) NULL) goto FAILURE;
        cgr->parent = ngr;
        cgr->sibling = ngr->child;
        ngr->child = cgr;
    }

    return ngr;

  FAILURE:
    if (ngr) {
        subcutgrtree_free (ngr, cutgrnode_world, cutgradj_world);
    }
    return (cutgrnode *) NULL;
}

static void cutgrnode_init (cutgrnode *n)
{
    n->parent = (cutgrnode *) NULL;
    n->sibling = (cutgrnode *) NULL;
    n->child = (cutgrnode *) NULL;
    n->adj = (cutgradj *) NULL;
    n->type = CCtsp_CUT_EXTERN;
}

static int PQ_tree_load_segments (CCpq_tree *pqt, double cut_tol,
                                  int edgecount, int *elist, double *x)

{
    int rval;
    struct pq_callback_struct pcs;
    int i;
    int nodecount = pqt->nodecount;
    int *endmark = (int *) NULL;
    int status;
    int nontriv;

    endmark = CC_SAFE_MALLOC (nodecount, int);
    if (endmark == (int *) NULL) {
        fprintf (stderr, "Out of memory in PQ_tree_load_segments\n");
        rval = -1; goto CLEANUP;
    }

    for (i=0; i<nodecount; i++) {
        endmark[i] = 0;
    }
    nontriv = 0;
    for (i=0; i<edgecount; i++) {
        if (x[i] >= 1.0 - cut_tol/2.0) {
            endmark[elist[2*i]]++;
            endmark[elist[2*i+1]]++;
            CCpq_clear_leaflist (pqt);
            CCpq_add_leaflist (pqt, elist[2*i]);
            CCpq_add_leaflist (pqt, elist[2*i+1]);
            rval = CCpq_apply (pqt, &status);
            if (rval) {
                fprintf (stderr, "CCpq_apply failed\n");
                goto CLEANUP;
            }
            if (status == CCpq_STATUS_NONTRIVIAL) {
                nontriv++;
            }
        }
    }

#ifdef VERBOSE
    if (nontriv) {
        printf ("%d nontrivial edge cuts added\n", nontriv);
    }
#endif

    for (i=0; i<nodecount; i++) {
        if (endmark[i] == 2) {
            endmark[i] = CC_LINSUB_NO_END;
        } else {
            endmark[i] = CC_LINSUB_BOTH_END;
        }
    }

    pcs.pqt = pqt;
    pcs.nontrivial = 0;
    pcs.cutcount = 0;

    rval = CCcut_linsub_allcuts (nodecount, edgecount, (int *) NULL, endmark,
                           elist, x, 2.0 + cut_tol, (void *) &pcs,
                           load_seg_callback);
    if (rval) {
        fprintf (stderr, "libsub_allcuts failed\n");
        goto CLEANUP;
    }

#ifdef VERBOSE
    printf ("CCcut_linsub_allcuts found %d (%d nontriv) cuts in %.2f seconds\n",
            pcs.cutcount, pcs.nontrivial, CCutil_zeit() - szeit);
    fflush (stdout);
#endif

    rval = 0;

  CLEANUP:

    CC_IFFREE (endmark, int);
    return rval;
}

static int load_seg_callback (CC_UNUSED double cut_val, int cut_start,
        int cut_end, void *callbackdat)
{
    struct pq_callback_struct *pcs = (struct pq_callback_struct *)
        callbackdat;
    CCpq_tree *pqt = pcs->pqt;
    int i;
    int rval;
    int status;

    if (cut_start == cut_end) return 0;

    CCpq_clear_leaflist (pqt);
    for (i=cut_start; i<=cut_end; i++) {
        CCpq_add_leaflist (pqt, i);
    }

    rval = CCpq_apply (pqt, &status);
    if (rval) {
        fprintf (stderr, "CCpq_apply failed\n");
        return rval;
    }
    if (status == CCpq_STATUS_NONTRIVIAL) pcs->nontrivial++;
    pcs->cutcount++;

    return 0;
}

static int PQ_tree_load_pool (CCpq_tree *pqt, double cut_tol,
        CCtsp_lpcuts *pool, int edgecount, int *elist, double *x)
{
    int *cliquenums = (int *) NULL;
    int cliquecount;
    int i;
    int rval;
    CCtsp_lpclique *c;
    int nontrivial = 0;
    int status;

    rval = CCtsp_get_clique_prices (pool, &cliquenums, (double **) NULL,
                                    2.0 - cut_tol, 2.0 + cut_tol,
                                    &cliquecount, pqt->nodecount, edgecount,
                                    elist, x);
    if (rval) {
        fprintf (stderr, "CCtsp_get_clique_prices failed\n");
        goto CLEANUP;
    }

    for (i=0; i<cliquecount; i++) {
        rval = CCtsp_get_clique (pool, cliquenums[i], &c);
        if (rval) {
            fprintf (stderr, "CCtsp_get_clique failed\n");
            goto CLEANUP;
        }
        rval = CCpq_apply_clique (pqt, c, &status);
        if (rval) {
            fprintf (stderr, "CCpq_apply_clique failed\n");
            goto CLEANUP;
        }
        if (status == CCpq_STATUS_NONTRIVIAL) nontrivial++;
    }

#ifdef VERBOSE
    printf ("pool found %d (%d nontriv) cuts in %.2f seconds\n",
            cliquecount, nontrivial, CCutil_zeit() - szeit);
    fflush (stdout);
#endif

    rval = 0;

  CLEANUP:
    CC_IFFREE (cliquenums, int);
    return rval;
}

int CCpq_apply_clique (CCpq_tree *pqt, CCtsp_lpclique *c, int *status)
{
    int i;
    int j;
    int rval;

    CCpq_clear_leaflist (pqt);
    CC_FOREACH_NODE_IN_CLIQUE (i, *c, j) {
        CCpq_add_leaflist (pqt, i);
    }

    rval = CCpq_apply (pqt, status);
    if (rval) {
        fprintf (stderr, "CCpq_apply failed\n");
        goto CLEANUP;
    }

  CLEANUP:
    return rval;
}

void CCpq_check_clique (CCpq_tree *pqt, CCtsp_lpclique *c, int *status)
{
    int i;
    int j;

    CCpq_clear_leaflist (pqt);
    CC_FOREACH_NODE_IN_CLIQUE (i, *c, j) {
        CCpq_add_leaflist (pqt, i);
    }

    CCpq_check (pqt, status);
}

static void cutgrtree_init (cutgrtree *t)
{
    t->root = (cutgrnode *) NULL;
    t->nodelist = (cutgrnode *) NULL;
    t->nodecount = -1;
    t->extern_node = -1;
    CCptrworld_init (&t->cutgrnode_world);
    CCptrworld_init (&t->cutgradj_world);
}

static void cutgrtree_freetree (cutgrtree *t)
{
    int total;
    int onlist;

    if (t->root) {
        subcutgrtree_free (t->root, &t->cutgrnode_world, &t->cutgradj_world);
        t->root = (cutgrnode *) NULL;
    }
    CC_IFFREE (t->nodelist, cutgrnode);

    if (cutgradj_check_leaks (&t->cutgradj_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding cutgradjs\n",
                 total - onlist);
    }
    CCptrworld_delete (&t->cutgradj_world);

    if (cutgrnode_check_leaks (&t->cutgrnode_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding cutgrnodes\n",
                 total - onlist);
    }
    CCptrworld_delete (&t->cutgrnode_world);
}

static void subcutgrtree_free (cutgrnode *n, CCptrworld *cutgrnode_world,
        CCptrworld *cutgradj_world)
{
    cutgrnode *p;
    cutgrnode *pnext;

    cutgradj_listfree (cutgradj_world, n->adj);
    n->adj = (cutgradj *) NULL;

    if (n->type == CCtsp_CUT_LEAF || n->type == CCtsp_CUT_EXTERN) return;

    for (p = n->child; p; p = pnext) {
        pnext = p->sibling;
        subcutgrtree_free (p, cutgrnode_world, cutgradj_world);
    }
    cutgrnode_free (cutgrnode_world, n);
}

static int cutgrtree_build_leafgraph (cutgrtree *t, int edgecount,
                                      int *elist, double *x)
{
    int i;
    cutgrnode *nodelist = t->nodelist;
    cutgrnode *e1;
    cutgrnode *e2;
    int rval;

    for (i=0; i<edgecount; i++) {
        if (x[i] > 0.0) {
            e1 = nodelist + elist[2*i];
            e2 = nodelist + elist[2*i+1];
            rval = cutgrtree_addadj (e1, e1, e2, x[i], i, &t->cutgradj_world);
            if (rval) {
                fprintf (stderr, "cutgrtree_addadj failed\n");
                return rval;
            }
            rval = cutgrtree_addadj (e2, e2, e1, x[i], i, &t->cutgradj_world);
            if (rval) {
                fprintf (stderr, "cutgrtree_addadj failed\n");
                return rval;
            }
        }
    }
    return 0;
}

static int cutgrtree_addadj (cutgrnode *n, cutgrnode *from, cutgrnode *to,
        double weight, int num, CCptrworld *cutgradj_world)
{
    cutgradj *a = cutgradj_alloc (cutgradj_world);

    if (a == (cutgradj *) NULL) {
        fprintf (stderr, "Out of memory in cutgrtree_addadj\n");
        return -1;
    }
    a->to = to;
    a->from = from;
    a->weight = weight;
    a->num = num;
    a->next = n->adj;
    n->adj = a;

    return 0;
}

static void cutgrtree_free_leafgraph (cutgrtree *t)
{
    cutgrnode *nodelist = t->nodelist;
    int nodecount = t->nodecount;
    int i;

    for (i=0; i<nodecount; i++) {
        cutgradj_listfree (&t->cutgradj_world, nodelist[i].adj);
        nodelist[i].adj = (cutgradj *) NULL;
    }
}

static int cutgrtree_loadx (cutgrtree *t, int edgecount, int *elist,
                            double *x, int cleanup)
{
    int rval;

    cutgrtree_leafgraph_init (t->root);

    rval = cutgrtree_build_leafgraph (t, edgecount, elist, x);
    if (rval) {
        fprintf (stderr, "cutgrtree_build_leafgraph failed\n");
        cutgrtree_free_leafgraph (t);
        return rval;
    }
    rval = cutgrtree_ancestor (t->root, &t->cutgradj_world);
    if (rval) {
        fprintf (stderr, "cutgrtree_ancestor failed\n");
        cutgrtree_free_leafgraph (t);
        return rval;
    }
    rval = cutgrtree_penultimate (t->root, &t->cutgradj_world);
    if (rval) {
        fprintf (stderr, "cutgrtree_penultimate failed\n");
        cutgrtree_free_leafgraph (t);
        return rval;
    }

    if (cleanup) {
        cutgrtree_cleanadj (t->root, &t->cutgradj_world);
    }

#ifdef VERBOSE
    printf ("cutgrtree graph loaded in %.2f seconds\n", CCutil_zeit() - szeit);
    fflush (stdout);
#endif

#ifdef DEBUG
    szeit = CCutil_zeit();
    verify_cutgrtree (t, t->root, edgecount, elist, x);
    printf ("tree verified in %.2f seconds\n", CCutil_zeit() - szeit);
    fflush (stdout);
#endif

    return 0;
}

static void cutgrtree_leafgraph_init (cutgrnode *n)
{
    cutgrnode *c;

    n->mark = 0;
    n->int_weight = 0.0;
    n->ext_weight = 0.0;
    n->subtree_size = 0;
    n->adjinto = (cutgradj *) NULL;
    n->adj = (cutgradj *) NULL;
    for (c = n->child; c; c = c->sibling) {
        cutgrtree_leafgraph_init (c);
    }
}

/* note: build_leafgraph places each edge on 2 leaf adjacency lists.  lca
   frees one copy and places the other on the lca's adjacency list.
*/

static int cutgrtree_ancestor (cutgrnode *n, CCptrworld *cutgradj_world)
{
    cutgrnode *c;
    cutgradj *a;
    cutgradj *anext;
    cutgrnode *anc;
    int rval;

    make_set (n);
    find_set (n)->ancestor = n;
    for (c = n->child; c; c = c->sibling) {
        rval = cutgrtree_ancestor (c, cutgradj_world);
        if (rval) return rval;
        link_set (n, c)->ancestor = n;
    }
    n->mark = 1;
    if (n->type == CCtsp_CUT_LEAF || n->type == CCtsp_CUT_EXTERN) {
        for (a = n->adj, n->adj = (cutgradj *) NULL; a; a = anext) {
            anext = a->next;
            if (a->to->mark == 1) {
                anc = find_set (a->to)->ancestor;
                rval = cutgrtree_addadj (anc, a->from, a->to, a->weight,
                                         a->num, cutgradj_world);
                if (rval) {
                    fprintf (stderr, "cutgrtree_addadj failed\n");
                    return rval;
                }
                anc->int_weight += a->weight;
            }
            cutgradj_free (cutgradj_world, a);
        }
    }
    return 0;
}

static int cutgrtree_penultimate (cutgrnode *n, CCptrworld *cutgradj_world)
{
    cutgrnode *c;
    cutgradj *a;
    cutgradj *anext;
    cutgrnode *f;
    cutgrnode *t;
    int rval;

    for (c = n->child; c; c = c->sibling) {
        rval = cutgrtree_penultimate (c, cutgradj_world);
        if (rval) return rval;
    }

    for (a = n->adj, n->adj = (cutgradj *) NULL; a; a = anext) {
        anext = a->next;
        f = find_set (a->from)->ancestor;
        t = find_set (a->to)->ancestor;

#ifdef DEBUG
        if (f == t) {
            fprintf (stderr, "Whoops, f = t\n");
        }
        if (f->parent != n) {
            fprintf (stderr, "Whoops, f->parent != n\n");
        }
        if (t->parent != n) {
            fprintf (stderr, "Whoops, t->parent != n\n");
        }
#endif

        rval = cutgrtree_addadj (f, f, t, a->weight, a->num, cutgradj_world);
        if (rval) {
            fprintf (stderr, "cutgrtree_addadj failed\n");
            return rval;
        }
        rval = cutgrtree_addadj (t, t, f, a->weight, a->num, cutgradj_world);
        if (rval) {
            fprintf (stderr, "cutgrtree_addadj failed\n");
            return rval;
        }
        cutgradj_free (cutgradj_world, a);
    }

    make_set (n);
    for (c = n->child; c; c = c->sibling) {
        link_set (n, c);
    }
    find_set (n)->ancestor = n;
    return 0;
}

static void cutgrtree_cleanadj (cutgrnode *n, CCptrworld *cutgradj_world)
{
    cutgrnode *c;
    double sum = n->int_weight;
    int size;
    cutgradj *a;
    cutgradj *anext;
    cutgrnode *v;
    double ext_weight;

    if (n->type == CCtsp_CUT_LEAF || n->type == CCtsp_CUT_EXTERN) {
        size = 1;
    } else {
        size = 0;
    }

    for (c = n->child; c; c = c->sibling) {
        cutgrtree_cleanadj (c, cutgradj_world);
        sum += c->int_weight;
        size += c->subtree_size;
    }
    n->int_weight = sum;
    n->subtree_size = size;

    ext_weight = 2.0 * (size - sum);

    for (a=n->adj, n->adj=(cutgradj *) NULL; a; a = anext) {
        anext = a->next;
        v = a->to;
        ext_weight -= a->weight;
        if (v->adjinto == (cutgradj *) NULL) {
            v->adjinto = a;
            a->next = n->adj;
            n->adj = a;
        } else {
            v->adjinto->weight += a->weight;
            cutgradj_free (cutgradj_world, a);
        }
    }
    for (a=n->adj; a; a = a->next) {
        a->to->adjinto = (cutgradj *) NULL;
    }
    n->ext_weight = ext_weight;
}

static int PQ_add_tight_nodes (double tol, CCpq_tree *pqt, cutgrtree *t)
{
    int rval;

    build_leaf_elems (t->root, pqt->elems, t->nodelist);

    rval = add_tight_nodes_work (t->root, tol, pqt);
    if (rval) {
        fprintf (stderr, "add_tight_nodes_work failed\n");
    }
    return rval;
}

static void build_leaf_elems (cutgrnode *n, CCpq_node *elems,
                              cutgrnode *nodelist)
{
    cutgrnode *c;
    CCpq_node **lasttail;
    CCpq_node *pqn;

    if (n->type == CCtsp_CUT_LEAF || n->type == CCtsp_CUT_EXTERN) {
        pqn = &elems[n - nodelist];
        pqn->next = (CCpq_node *) NULL;
        n->leaf_elems = pqn;
        n->leaf_elems_end = &(pqn->next);
    } else {
        lasttail = &(n->leaf_elems);
        for (c = n->child; c; c = c->sibling) {
            build_leaf_elems (c, elems, nodelist);
            *lasttail = c->leaf_elems;
            lasttail = c->leaf_elems_end;
        }
        n->leaf_elems_end = lasttail;
        *lasttail = (CCpq_node *) NULL;
    }
}

static int add_tight_nodes_work (cutgrnode *n, double tol, CCpq_tree *pqt)
{
    cutgrnode *c;
    int rval;
    double cut_weight = 2.0 * (n->subtree_size - n->int_weight);

    if ((n->type == CCtsp_CUT_PNODE || n->type == CCtsp_CUT_QNODE) &&
        (cut_weight >= 2.0-tol && cut_weight <= 2.0+tol)) {
        rval = add_node (n, pqt);
        if (rval) {
            fprintf (stderr, "add_node failed\n");
            return rval;
        }
    }

    for (c = n->child; c; c = c->sibling) {
        rval = add_tight_nodes_work (c, tol, pqt);
        if (rval) return rval;
    }

    return 0;
}

static int add_node (cutgrnode *n, CCpq_tree *pqt)
{
    int rval;
    int status;
    CCpq_node *save_next;

    save_next = *(n->leaf_elems_end);
    *(n->leaf_elems_end) = (CCpq_node *) NULL;
    CCpq_set_leaflist (pqt, n->leaf_elems);
    rval = CCpq_apply (pqt, &status);
    *(n->leaf_elems_end) = save_next;
    if (rval) {
        fprintf (stderr, "CCpq_apply failed\n");
        return rval;
    }
    if (status == CCpq_STATUS_NOSOL) {
        printf ("WARNING: add_node tight cut wouldn't fit\n");
    }
    return 0;
}

static int add_node_pair (cutgrnode *n1, cutgrnode *n2, CCpq_tree *pqt)
{
    int rval;
    int status;
    CCpq_node *save_next1;
    CCpq_node *save_next2;

    save_next1 = *(n1->leaf_elems_end);
    save_next2 = *(n2->leaf_elems_end);
    *(n1->leaf_elems_end) = n2->leaf_elems;
    *(n2->leaf_elems_end) = (CCpq_node *) NULL;
    CCpq_set_leaflist (pqt, n1->leaf_elems);
    rval = CCpq_apply (pqt, &status);
    *(n1->leaf_elems_end) = save_next1;
    *(n2->leaf_elems_end) = save_next2;
    if (rval) {
        fprintf (stderr, "CCpq_apply failed\n");
        return rval;
    }
    if (status == CCpq_STATUS_NOSOL) {
        printf ("WARNING: add_node_pair tight cut wouldn't fit\n");
    }
    return 0;
}

static int add_node_complement (cutgrnode *p, cutgrnode *n, CCpq_tree *pqt)
{
    int rval;
    int status;
    cutgrnode *c;
    CCpq_node *list = (CCpq_node *) NULL;

    for (c = p->child; c; c = c->sibling) {
        if (c != n) {
            c->save_leaf_next = *(c->leaf_elems_end);
            *(c->leaf_elems_end) = list;
            list = c->leaf_elems;
        }
    }
    CCpq_set_leaflist (pqt, list);
    rval = CCpq_apply (pqt, &status);

    for (c = p->child; c; c = c->sibling) {
        if (c != n) {
            *(c->leaf_elems_end) = c->save_leaf_next;
        }
    }

    if (rval) {
        fprintf (stderr, "CCpq_apply failed\n");
        return rval;
    }
    if (status == CCpq_STATUS_NOSOL) {
        printf ("WARNING: add_complement tight cut wouldn't fit\n");
    }
    return 0;
}

static int PQ_add_tight_edges (double tol, CCpq_tree *pqt, cutgrtree *t)
{
    int rval;

    tol /= 2.0;

    build_leaf_elems (t->root, pqt->elems, t->nodelist);
    mark_tree (t->root, 0);

    rval = add_tight_edges_work (t->root, tol, pqt);
    if (rval) {
        fprintf (stderr, "add_tight_edges_work failed\n");
    }
    return rval;
}

static int add_tight_edges_work (cutgrnode *n, double tol, CCpq_tree *pqt)
{
    cutgrnode *c;
    int rval;

    if (n->mark == 0) {
        rval = check_tight_path (n, tol, pqt);
        if (rval) {
            fprintf (stderr, "check_tight_path failed\n");
            return rval;
        }
        n->mark = 1;
    }

    rval = check_tight_ext (n, tol, pqt);
    if (rval) {
        fprintf (stderr, "check_tight_ext failed\n");
        return rval;
    }

    for (c = n->child; c; c = c->sibling) {
        rval = add_tight_edges_work (c, tol, pqt);
        if (rval) return rval;
    }
    return 0;
}

static int check_tight_path (cutgrnode *n, double tol, CCpq_tree *pqt)
{
    cutgradj *a;
    cutgrnode *end1;
    cutgrnode *end2;
    cutgrnode *end2alt;
    int pathlength = 0;
    int rval;

    a = find_path_edge (n, (cutgrnode *) NULL, tol);
    if (a == (cutgradj *) NULL) {
        return 0;
    }

    end2alt = a->to;
    end1 = n;
    while (a && a->to != n) {
        end1 = a->to;
        pathlength++;
        a = find_path_edge (end1, a->from, tol);
    }

    if (a) {
        rval = add_tight_path (end2alt, a->from, n, tol, pqt);
        if (rval) {
            fprintf (stderr, "add_tight_path failed\n");
        }
        return rval;
    }

    end2 = n;
    a = find_path_edge (n, end2alt, tol);
    while (a) {
        end2alt = end2;
        end2 = a->to;
        pathlength++;
        a = find_path_edge (end2, end2alt, tol);
    }

    if ((end1->ext_weight >= 1.0-tol && end1->ext_weight <= 1.0+tol) ||
        (end2->ext_weight >= 1.0-tol && end2->ext_weight <= 1.0+tol)) {
        fprintf (stderr, "Whoa, unexpected external path\n");
        return -1;
    }
    if (pathlength > 1) {
        a = find_edge (end1, end2);
        if (a && a->weight > tol) {
            end2->mark = 1;
            end2 = end2alt;
        }
    }

    rval = add_tight_path (end1, end2, (cutgrnode *) NULL, tol, pqt);
    if (rval) {
        fprintf (stderr, "add_tight_path failed\n");
    }
    return rval;
}

static cutgradj *find_edge (cutgrnode *end1, cutgrnode *end2)
{
    cutgradj *a;

    for (a = end1->adj; a; a = a->next) {
        if (a->to == end2) return a;
    }
    return (cutgradj *) NULL;
}

static cutgradj *find_path_edge (cutgrnode *n, cutgrnode *avoid, double tol)
{
    cutgradj *a;

    for (a = n->adj; a; a = a->next) {
        if (a->weight >= 1.0 - tol && a->weight <= 1.0 + tol &&
            a->to != avoid) {
            return a;
        }
    }
    return (cutgradj *) NULL;
}

static cutgrnode *find_ext_node (cutgrnode *n, cutgrnode *avoid, double tol)
{
    cutgrnode *c;

    for (c = n->child; c; c = c->sibling) {
        if (c->ext_weight >= 1.0 - tol && c->ext_weight <= 1.0 + tol &&
            c != avoid) {
            return c;
        }
    }
    return (cutgrnode *) NULL;
}

static int add_tight_path (cutgrnode *end1, cutgrnode *end2,
                           cutgrnode *avoid, double tol, CCpq_tree *pqt)
{
    cutgradj *a;
    int rval;

    end1->mark = 1;
    while (end1 != end2) {
        a = find_path_edge (end1, avoid, tol);
        if (a == (cutgradj *) NULL) {
            fprintf (stderr, "Whoa, path vanished in add_tight_path\n");
            return -1;
        }
        rval = add_node_pair (a->from, a->to, pqt);
        if (rval) {
            fprintf (stderr, "add_node_pair failed\n");
            return rval;
        }
        avoid = end1;
        end1 = a->to;
        end1->mark = 1;
    }
    return 0;
}

static int check_tight_ext (cutgrnode *n, double tol, CCpq_tree *pqt)
{
    cutgrnode *ext_end1;
    cutgrnode *int_end1;
    cutgrnode *int_end1_alt;
    cutgrnode *ext_end2;
    cutgrnode *int_end2;
    cutgrnode *int_end2_alt;
    cutgradj *a;
    int rval;

    ext_end1 = find_ext_node (n, (cutgrnode *) NULL, tol);
    if (ext_end1 == (cutgrnode *) NULL) {
        return 0;
    }

    int_end1_alt = (cutgrnode *) NULL;
    int_end1 = ext_end1;
    a = find_path_edge (int_end1, (cutgrnode *) NULL, tol);
    while (a) {
        int_end1_alt = int_end1;
        int_end1 = a->to;
        a = find_path_edge (int_end1, int_end1_alt, tol);
    }

    ext_end2 = find_ext_node (n, ext_end1, tol);
    if (ext_end2 == (cutgrnode *) NULL) {
        int_end2 = (cutgrnode *) NULL;
        int_end2_alt = (cutgrnode *) NULL;
    } else {
        int_end2_alt = (cutgrnode *) NULL;
        int_end2 = ext_end2;
        a = find_path_edge (int_end2, (cutgrnode *) NULL, tol);
        while (a) {
            int_end2_alt = int_end2;
            int_end2 = a->to;
            a = find_path_edge (int_end2, int_end2_alt, tol);
        }
    }

    rval = add_tight_path (int_end1, ext_end1, (cutgrnode *) NULL,
                           tol, pqt);
    if (rval) {
        fprintf (stderr, "add_tight_path failed\n");
        return rval;
    }

    if (ext_end2 == (cutgrnode *) NULL) {
        if (int_end1 == ext_end1 || int_end1->ext_weight <= tol) {
            rval = add_node_complement (n, ext_end1, pqt);
            if (rval) {
                fprintf (stderr, "add_node_complement failed\n");
                return rval;
            }
        }
        return 0;
    }
    if (ext_end2 == int_end1) {
        /* external cycle, done */
        return 0;
    }
    if (ext_end1 == int_end2) {
        fprintf (stderr, "Whoa, assymetric path\n");
        return -1;
    }
    a = find_edge (int_end1, int_end2);
    if (a && a->weight > tol) {
        if (int_end2_alt == (cutgrnode *) NULL) {
            rval = add_node_complement (n, ext_end1, pqt);
            if (rval) {
                fprintf (stderr, "add_node_complement failed\n");
            }
            return rval;
        }
        int_end2->mark = 1;
        int_end2 = int_end2_alt;
    }
    rval = add_tight_path (ext_end2, int_end2, (cutgrnode *) NULL,
                           tol, pqt);
    if (rval) {
        fprintf (stderr, "add_tight_path failed\n");
        return rval;
    }
    rval = add_node_complement (n, ext_end1, pqt);
    if (rval) {
        fprintf (stderr, "add_node_complement failed\n");
        return rval;
    }
    rval = add_node_complement (n, ext_end2, pqt);
    if (rval) {
        fprintf (stderr, "add_node_complement failed\n");
        return rval;
    }
    return 0;
}

static void make_set (cutgrnode *x)
{
    x->disj_set.parent = x;
    x->disj_set.rank = 0;
}

static cutgrnode *find_set (cutgrnode *x)
{
    cutgrnode *p;
    cutgrnode *xnext;

    for (p = x; p != p->disj_set.parent; p = p->disj_set.parent)
        ;

    while (x != p) {
        xnext = x->disj_set.parent;
        x->disj_set.parent = p;
        x = xnext;
    }
    return p;
}

static cutgrnode *link_set (cutgrnode *x, cutgrnode *y)
{
    x = find_set (x);
    y = find_set (y);
    if (x == y) return x;

    if (x->disj_set.rank < y->disj_set.rank) {
        x->disj_set.parent = y;
        return y;
    } else if (x->disj_set.rank > y->disj_set.rank) {
        y->disj_set.parent = x;
        return x;
    } else {
        x->disj_set.parent = y;
        y->disj_set.rank++;
        return y;
    }
}

static int mark_tree (cutgrnode *n, int m)
{
    cutgrnode *c;
    int cnt;

    n->mark = m;
    if (n->type == CCtsp_CUT_LEAF || n->type == CCtsp_CUT_EXTERN) {
        cnt = 1;
    } else {
        cnt = 0;
    }
    for (c = n->child; c; c = c->sibling) {
        cnt += mark_tree (c, m);
    }
    return cnt;
}

int CCpq_cuttree_gen_cliques (CCtsp_cuttree *t, void *u_data,
        int (*cut_callback) (int *arr, int cnt, int *stop, void *u_data))
{
    int rval;
    int *nodenums = (int *) NULL;
    int size = 0;
    int stop = 0;

    nodenums = CC_SAFE_MALLOC (t->nodecount, int);
    if (nodenums == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCpq_cuttree_gen_cliques\n");
        rval = 1; goto CLEANUP;
    }

    rval = cuttree_gen_work (t->root, t->nodelist, nodenums, &size,
                             &stop, u_data, cut_callback);
    if (rval) {
        fprintf (stderr, "cuttree_gen_work failed\n");
        goto CLEANUP;
    }

    rval = 0;

  CLEANUP:
    CC_IFFREE (nodenums, int);
    return rval;
}

static int cuttree_gen_work (CCtsp_cutnode *n, CCtsp_cutnode *nodelist,
        int *nodenums, int *size, int *stop, void *u_data,
        int (*cut_callback) (int *arr, int cnt, int *stop, void *u_data))
{
    CCtsp_cutnode *c;
    int loc1;
    int loc2;
    int rval;

    if (n->type == CCtsp_CUT_LEAF || n->type == CCtsp_CUT_EXTERN) {
        nodenums[*size] = (int) (n - nodelist);
        (*size)++;
    } else if (n->type == CCtsp_CUT_PNODE) {
        loc1 = *size;
        for (c = n->child; c; c = c->sibling) {
            rval = cuttree_gen_work (c, nodelist, nodenums, size, stop,
                                     u_data, cut_callback);
            if (rval) return rval;
            if (*stop) return 0;
        }
        rval = (*cut_callback) (nodenums + loc1, (*size) - loc1,
                                stop, u_data);
        if (rval) {
            fprintf (stderr, "cut_callback failed\n");
            return rval;
        }
        if (*stop) return 0;
    } else if (n->type == CCtsp_CUT_QNODE) {
        loc1 = -1;
        loc2 = *size;
        for (c = n->child; c; c = c->sibling) {
            rval = cuttree_gen_work (c, nodelist, nodenums, size, stop,
                                     u_data, cut_callback);
            if (rval) return rval;
            if (*stop) return 0;
            if (loc1 != -1) {
                rval = (*cut_callback) (nodenums + loc1, (*size) - loc1,
                                        stop, u_data);
                if (rval) {
                    fprintf (stderr, "cut_callback failed\n");
                    return rval;
                }
                if (*stop) return 0;
            }
            loc1 = loc2;
            loc2 = *size;
        }
    } else {
        for (c = n->child; c; c = c->sibling) {
            rval = cuttree_gen_work (c, nodelist, nodenums, size, stop,
                                     u_data, cut_callback);
            if (rval) return rval;
            if (*stop) return 0;
        }
    }
    return 0;
}

int CCpq_cuttree_build_necklaces (CCtsp_cuttree *t, int ecount, int *elist,
        double *x, int *p_neckcount, CCtsp_cutnode ***p_necklist, int *necknum)
{
    cutgrtree cgt;
    int neckcount;
    CCtsp_cutnode **necklist = (CCtsp_cutnode **) NULL;
    int i;
    int rval;

    cutgrtree_init (&cgt);
    *p_neckcount = 0;
    *p_necklist = (CCtsp_cutnode **) NULL;

    rval = cuttree_to_cutgrtree (t, &cgt);
    if (rval) {
        fprintf (stderr, "cuttree_to_cutgrtree failed\n");
        goto CLEANUP;
    }

    rval = cutgrtree_loadx (&cgt, ecount, elist, x, 0);
    if (rval) {
        fprintf (stderr, "cutgrtree_loadx failed\n");
        goto CLEANUP;
    }

    neckcount = count_necklaces (cgt.root);

    necklist = CC_SAFE_MALLOC (neckcount, CCtsp_cutnode *);
    if (necklist == (CCtsp_cutnode **) NULL) {
        fprintf (stderr, "Out of memory in CCpq_cuttree_build_necklaces\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<ecount; i++) {
        necknum[i] = -1;
    }

    i = 0;
    label_necklaces (cgt.root, necklist, &i, necknum);
    if (i != neckcount) {
        fprintf (stderr, "label_necklaces != count_necklaces\n");
        rval = 1; goto CLEANUP;
    }
    *p_neckcount = neckcount;
    *p_necklist = necklist;

    rval = 0;

  CLEANUP:
    if (rval) {
        CC_IFFREE (necklist, CCtsp_cutnode *);
    }
    cutgrtree_freetree (&cgt);
    return rval;
}

static int count_necklaces (cutgrnode *n)
{
    int cnt = 0;
    int childcount = 0;
    cutgrnode *c;

    for (c = n->child; c; c = c->sibling) {
        cnt += count_necklaces (c);
        childcount++;
    }

    if (n->type == CCtsp_CUT_QNODE || (n->type == CCtsp_CUT_PNODE &&
                                       childcount == 2)) {
        cnt++;
    }
    return cnt;
}

static void label_necklaces (cutgrnode *n, CCtsp_cutnode **necklist,
        int *cnt, int *necknum)
{
    int childcount = 0;
    cutgrnode *c;
    cutgrnode *c2;
    cutgradj *a;

    for (c = n->child; c; c = c->sibling) {
        label_necklaces (c, necklist, cnt, necknum);
        childcount++;
    }

    if (n->type == CCtsp_CUT_QNODE || (n->type == CCtsp_CUT_PNODE &&
                                       childcount == 2)) {
        c = n->child;
        c2 = c->sibling;
        for (a = c->adj; a; a = a->next) {
            if (a->to == c2) {
                necknum[a->num] = *cnt;
            }
        }
        necklist[(*cnt)++] = n->cn_base;
    }
}

void CCpq_cuttree_display (CCtsp_cuttree *t)
{
    cuttree_display_work (t->root, t->nodelist);
    printf ("\n");
    fflush (stdout);
}

void CCpq_cuttree_describe (CCtsp_cuttree *t)
{
    cuttree_describe_work (t->root);
    printf ("\n");
    fflush (stdout);
}

static void cuttree_display_work (CCtsp_cutnode *n, CCtsp_cutnode *nodelist)
{
    CCtsp_cutnode *c;

    if (n->type == CCtsp_CUT_EXTERN || n->type == CCtsp_CUT_LEAF) {
        printf ("%d", (int) (n - nodelist));
    } else if (n->type == CCtsp_CUT_ROOT) {
        printf ("{");
        for (c = n->child; c; c = c->sibling) {
            cuttree_display_work (c, nodelist);
            if (c->sibling) printf (" ");
        }
        printf ("}");
    } else if (n->type == CCtsp_CUT_PNODE) {
        printf ("(");
        for (c = n->child; c; c = c->sibling) {
            cuttree_display_work (c, nodelist);
            if (c->sibling) printf (" ");
        }
        printf (")");
    } else if (n->type == CCtsp_CUT_QNODE) {
        printf ("[");
        for (c = n->child; c; c = c->sibling) {
            cuttree_display_work (c, nodelist);
            if (c->sibling) printf (" ");
        }
        printf ("]");
    } else {
        printf ("?{");
        for (c = n->child; c; c = c->sibling) {
            cuttree_display_work (c, nodelist);
            if (c->sibling) printf (" ");
        }
        printf ("}?");
    }
}

static void cuttree_describe_work (CCtsp_cutnode *n)
{
    CCtsp_cutnode *c;
    int lcnt;
    int ccnt;

    if (n->type == CCtsp_CUT_EXTERN || n->type == CCtsp_CUT_LEAF) {
        return;
    } else if (n->type == CCtsp_CUT_ROOT || n->type == CCtsp_CUT_PNODE) {
        printf (n->type == CCtsp_CUT_ROOT ? "{" : "(");
        lcnt = 0;
        ccnt = 0;
        for (c = n->child; c; c = c->sibling) {
            if (c->type == CCtsp_CUT_EXTERN || c->type == CCtsp_CUT_LEAF) {
                lcnt++;
            }
            ccnt++;
        }
        if (lcnt) {
            printf ("%d-L", lcnt);
        }
        if (ccnt > lcnt) {
            printf (" ");
            for (c = n->child; c; c = c->sibling) {
                if (c->type != CCtsp_CUT_EXTERN && c->type != CCtsp_CUT_LEAF) {
                    cuttree_describe_work (c);
                    lcnt++;
                    if (lcnt < ccnt) printf (" ");
                }
            }
        }
        printf (n->type == CCtsp_CUT_ROOT ? "}" : ")");
        printf ("<%d>", ccnt);
    } else {
        printf (n->type == CCtsp_CUT_QNODE ? "[" : "?{");
        lcnt = 0;
        ccnt = 0;
        for (c = n->child; c; c = c->sibling) {
            if (c->type == CCtsp_CUT_LEAF || c->type == CCtsp_CUT_EXTERN) {
                lcnt++;
            } else {
                if (lcnt) {
                    printf ("%d-L ", lcnt);
                    lcnt = 0;
                }
                cuttree_describe_work (c);
                if (c->sibling) printf (" ");
            }
            ccnt++;
        }
        if (lcnt) {
            printf ("%d-L", lcnt);
        }
        printf (n->type == CCtsp_CUT_QNODE ? "]" : "?}");
        printf ("<%d>", ccnt);
    }
}

#ifdef DEBUG

static void verify_cutgrtree (cutgrtree *t, cutgrnode *n, int edgecount,
                         int *elist, double *x)
{
    double int_sum = 0.0;
    double int_sum2 = 0.0;
    double ext_sum = 0.0;
    cutgrnode *nodelist = t->nodelist;
    int cnt;
    int i;
    cutgrnode *c;
    cutgradj *a;

    cnt = mark_tree (n, 1);
    for (i=0; i<edgecount; i++) {
        if (nodelist[elist[2*i]].mark == 1 &&
            nodelist[elist[2*i+1]].mark == 1) {
            int_sum += x[i];
        }
    }
    for (c = n->child; c; c = c->sibling) {
        for (a=c->adj; a; a = a->next) {
            int_sum2 += a->weight;
        }
        int_sum2 += 2*c->int_weight;
    }
    int_sum2 /= 2.0;

    ext_sum = 0.0;
    if (n->parent) {
        mark_tree (n->parent, 2);
        mark_tree (n, 1);
        for (i=0; i<edgecount; i++) {
            if ((nodelist[elist[2*i]].mark == 0 &&
                 nodelist[elist[2*i+1]].mark == 1) ||
                (nodelist[elist[2*i]].mark == 1 &&
                 nodelist[elist[2*i+1]].mark == 0)) {
                ext_sum += x[i];
            }
        }
        mark_tree (n->parent, 0);
    } else {
        ext_sum = 0.0;
    }


    if (int_sum - n->int_weight > 0.0000001 ||
        int_sum - n->int_weight < -0.0000001 ||
        cnt != n->subtree_size ||
        ext_sum - n->ext_weight > 0.0000001 ||
        ext_sum - n->ext_weight < -0.0000001 ||
        int_sum2 - n->int_weight > 0.0000001 ||
        int_sum2 - n->int_weight < -0.0000001) {
        printf ("verify problem, node %d: cnt %d size %d\n",
                n - nodelist, cnt, n->subtree_size);
        printf ("   int_sum %.7f int_sum2 %.7f int_weight %.7f\n",
                int_sum, int_sum2, n->int_weight);
        printf ("   ext_sum %.7f ext_weight %.7f\n",
                ext_sum, n->ext_weight);
        fflush (stdout);
    } else {
        putchar ('.'); fflush (stdout);
    }
    mark_tree (n, 0);

    for (c = n->child; c; c = c->sibling) {
        verify_cutgrtree (t, c, edgecount, elist, x);
    }
}
#endif /* DEBUG */

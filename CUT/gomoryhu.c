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
/*                         GOMORY-HU TREES                                  */
/*                                                                          */
/*                            TSP CODE                                      */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: Spring 1989 (Dave)                                                */
/*        January 8, 1998 (Modified by Bico)                                */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCcut_gomory_hu (CC_GHtree *T, int ncount, int ecount,              */
/*      int *elist, double *ecap, int markcount, int *marks,                */
/*      CCrandstate *rstate)                                                */
/*    COMPUTES the Gomory-Hu tree of the marked nodes in G.                 */
/*     -T returns the tree (a description is given in the code below)       */
/*     -ncount, ecount, elist specify the input graph                       */
/*     -ecap lists the capacities of the edges                              */
/*     -markcount is the length of the array marks (if markcount is 0,      */
/*      then every node is a terminal)                                      */
/*     -marks lists the special nodes (the terminals)                       */
/*                                                                          */
/*  void CCcut_GHtreefree (CC_GHtree *T)                                    */
/*    FREES the tree pointed by T.                                          */
/*                                                                          */
/*  void CCcut_GHtreeinit (CC_GHtree *T)                                    */
/*    INITIALIZES the fields of T to NULL.                                  */
/*                                                                          */
/*  void CCcut_GHtreeprint (CC_GHtree *T)                                   */
/*       PRINTS the Gomory-Hu tree to stdout.                               */
/*                                                                          */
/*    NOTES:                                                                */
/*                                                                          */
/*      This code has only been tested on the instances that arise in       */
/*      exact blossom seperation.                                           */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "cut.h"

#define GH_MAXDOUBLE  (1e30)
#define OTHEREND(e,n) ((e)->ends[0] == (n) ? (e)->ends[1] : (e)->ends[0])

typedef struct edge {
    struct node          *ends[2];
    double                x;
    int                   magiclabel; 
    struct edge          *next;
} edge;

typedef struct edgeptr {
    struct edge          *this;
    struct edgeptr       *next;
} edgeptr;

typedef struct edgeset {
    struct edgeptr      *head;
    struct edgeptr      *tail;
} edgeset;

typedef struct node {
    edgeset              adj;
    int                  magiclabel;
    int                  pseudonumber;
    struct node         *next;
    int                  mark;
    int                  num;
    int                  number;
} node;

typedef struct nodeptr {
    struct node         *this;
    struct nodeptr      *next;
} nodeptr;

typedef struct nodeset {
    struct nodeptr      *head;
    struct nodeptr      *tail;
} nodeset;

typedef struct cuttree_node {
    struct cuttree_node *parent;
    struct cuttree_node *sibling;
    struct cuttree_node *child;
    double               cutval;
    int                  ndescendants;
    node                *special;
    nodeset              nlist;
    node                *pseudonode;
    int                  num;
    struct cuttree_node *next;
} cuttree_node;

typedef struct graph {
    int                  ncount;
    node                *nodelist;
    int                  ecount;
    edge                *edgelist;
    int                  magicnum;
    CCptrworld           edge_world;
    CCptrworld           edgeptr_world;
    CCptrworld           nodeptr_world;
} graph;


CC_PTRWORLD_ROUTINES (edge, edgealloc, edge_bulk_alloc, edgefree)
CC_PTRWORLD_LEAKS_ROUTINE (edge, edge_check_leaks, x, double)

CC_PTRWORLD_ROUTINES (edgeptr, edgeptralloc, edgeptr_bulk_alloc, edgeptrfree)
CC_PTRWORLD_LISTFREE_ROUTINE (edgeptr, edgeptr_listfree, edgeptrfree)
CC_PTRWORLD_LEAKS_ROUTINE (edgeptr, edgeptr_check_leaks, this, edge *)

CC_PTRWORLD_ROUTINES (nodeptr, nodeptralloc, nodeptr_bulk_alloc, nodeptrfree)
CC_PTRWORLD_LISTFREE_ROUTINE (nodeptr, nodeptr_listfree, nodeptrfree)
CC_PTRWORLD_LEAKS_ROUTINE (nodeptr, nodeptr_check_leaks, this, node *)


static void
    ghlink_free_world (graph *G),
    cuttree_free_work (cuttree_node *n, CCptrworld *nodeptr_world),
    delfromnodeset (nodeset *s, node *n, CCptrworld *nodeptr_world),
    delfromedgeset (edgeset *s, edge *e, CCptrworld *edgeptr_world),
    splitset (nodeset *s, nodeset *a, nodeset *b, int n),
    mergeset (nodeset *s, nodeset *a),
    unshrink (graph *G, node *pseudo, edgeset *esave, int num),
    initnodeset (nodeset *n),
    freenodeset (nodeset *s, CCptrworld *nodeptr_world),
    initgraph (graph *G),
    freegraph (graph *G),
    treefree_work (CC_GHnode *p),
    set_parents (cuttree_node *n),
    copy_dfs (cuttree_node *n, int *k, CC_GHnode *supply, int **listspace),
    cut_dfs (cuttree_node *n, int *k),
    print_tree_work (CC_GHnode *n);

static int
    addtonodeset (nodeset *s, node *n, CCptrworld *nodeptr_world),
    shrinkdown (graph *G, nodeset *a, node *pseudo, edgeset *esave, int num),
    gh_work (graph *G, cuttree_node *n, nodeset *nlist, nodeset *special,
        cuttree_node *supply, int *supplyhead, CCrandstate *rstate),
    countdescendants (cuttree_node *n),
    myrandnum (int n, CCrandstate *rstate),
    copy_cuttree (cuttree_node *root, int ncount, int markcount, CC_GHtree *T),
    buildgraph (graph *G, int ncount, int ecount, int *elist,
        double *ecap, int markcount, int *marks);


/*
 gomory_hu returns a pointer to the root of a rooted cut tree.  The tree is
   described by the parent, sibling, and child pointers.
 cutval is the weight on the cut (or, edge) between a node and its parent.
 ndescendants is the total number of nodes in the subtree rooted at that
   node, including itself.  This is the size of the cut between that node and
   its parent.
 special is the special node for this node of the tree.
 nlist is a nodeset containing all of the nodes of the graph grouped with the
   special node.  nlist contains special.
 pseudonode and next are work area.
 The calling routine should call cuttree_free(root) after it is done with the
   cuttree.

 gomory_hu constructs the min-cut tree for the nodes with mark = 1.
*/


int CCcut_gomory_hu (CC_GHtree *T, int ncount, int ecount, int *elist,
        double *ecap, int markcount, int *marks, CCrandstate *rstate)
{
    nodeset special;
    nodeset nlist;
    cuttree_node *root   = (cuttree_node *) NULL;
    cuttree_node *supply = (cuttree_node *) NULL;
    int supplyhead = 0;
    int i, rval = 0;
    graph G;

    initgraph (&G);
    initnodeset (&nlist);
    initnodeset (&special);

    rval = buildgraph (&G, ncount, ecount, elist, ecap, markcount, marks);
    if (rval) {
        fprintf (stderr, "buildgraph failed\n"); goto CLEANUP;
    }
    if (markcount == 0) markcount = ncount;   /* all nodes are terminals */

    supply = CC_SAFE_MALLOC (markcount + 1, cuttree_node);
    if (!supply) {
        fprintf (stderr, "out of memory in CCcut_gomory_hu\n");
        rval = 1; goto CLEANUP;
    }

    /* initialize root */
    root = &(supply[supplyhead++]);
    root->parent = (cuttree_node *) NULL;
    root->sibling = (cuttree_node *) NULL;
    root->child = (cuttree_node *) NULL;
    root->cutval = 0.0;
    root->pseudonode = (node *) NULL;
    root->nlist.head = root->nlist.tail = (nodeptr *) NULL;

    for (i = 0; i < G.ncount; i++) {
        rval = addtonodeset (&nlist, &(G.nodelist[i]), &G.nodeptr_world);
        if (rval) {
            fprintf (stderr, "addtonodeset failed\n");
            goto CLEANUP;
        }
        if (G.nodelist[i].mark == 1) {
            rval = addtonodeset (&special, &(G.nodelist[i]), &G.nodeptr_world);
            if (rval) {
                fprintf (stderr, "addtonodeset failed\n");
                goto CLEANUP;
            }
        }
    }

    rval = gh_work (&G, root, &nlist, &special, supply, &supplyhead, rstate);
    if (rval) {
        fprintf (stderr, "gh_work failed\n");
        goto CLEANUP;
    }

    freenodeset (&nlist, &G.nodeptr_world);
    freenodeset (&special, &G.nodeptr_world);

    countdescendants (root);
    root->cutval = GH_MAXDOUBLE;
    set_parents (root);

    rval = copy_cuttree (root, ncount, markcount, T);
    if (rval) {
        fprintf (stderr, "copy_cuttree failed\n"); fflush (stdout);
    }

    /* CCcut_GHtreeprint (T); */


CLEANUP:

    freegraph (&G);
    freenodeset (&nlist, &G.nodeptr_world);
    freenodeset (&special, &G.nodeptr_world);
    cuttree_free_work (root, &G.nodeptr_world);
    CC_IFFREE (supply, cuttree_node);
    ghlink_free_world (&G);

    return rval;
}

static void ghlink_free_world (graph *G)
{
    int total, onlist;

    if (edge_check_leaks (&G->edge_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding GH-edges\n",
                 total - onlist);
    }
    CCptrworld_delete (&G->edge_world);
    
    if (edgeptr_check_leaks (&G->edgeptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding GH-edgeptrs\n",
                 total - onlist);
    }
    CCptrworld_delete (&G->edgeptr_world);

    if (nodeptr_check_leaks (&G->nodeptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding GH-nodeptrs\n",
                 total - onlist);
    }
    CCptrworld_delete (&G->nodeptr_world);
}

static int addtonodeset (nodeset *s, node *n, CCptrworld *nodeptr_world)
{
    nodeptr *nnew;

    nnew = nodeptralloc (nodeptr_world);
    if (nnew == (nodeptr *) NULL) {
        fprintf (stderr, "nodeptralloc failed\n");
        return 1;
    }
    nnew->this = n;
    nnew->next = s->head;
    s->head = nnew;
    if (!s->tail) s->tail = nnew;
    return 0;
}

static void delfromnodeset (nodeset *s, node *n, CCptrworld *nodeptr_world)
{
    nodeptr *oldn, *n1;

    oldn = (nodeptr *) NULL;
    n1 = s->head;
    while (n1 && n1->this != n) {
        oldn = n1;
        n1 = n1->next;
    }
    if (!oldn) {                /* deleting head */
        s->head = n1->next;
        if (!s->head) s->tail = (nodeptr *) NULL;
    } else {
        oldn->next = n1->next;
        if (!oldn->next) s->tail = oldn;
    }
    nodeptrfree (nodeptr_world, n1);
}

static void delfromedgeset (edgeset *s, edge *e, CCptrworld *edgeptr_world)
{
    edgeptr *olde, *e1;

    olde = (edgeptr *) NULL;
    e1 = s->head;
    while (e1 && e1->this != e) {
        olde = e1;
        e1 = e1->next;
    }
    if (!olde) {                /* deleting head */
        s->head = e1->next;
        if (!s->head)
            s->tail = (edgeptr *) NULL;
    } else {
        olde->next = e1->next;
        if (!olde->next)
            s->tail = olde;
    }
    edgeptrfree (edgeptr_world, e1);
}

void CCcut_GHtreeinit (CC_GHtree *T)
{
    if (T) {
        T->root   = (CC_GHnode *) NULL;
        T->supply = (CC_GHnode *) NULL;
        T->listspace = (int *) NULL;
    }
}

void CCcut_GHtreefree (CC_GHtree *T)
{
    if (T) {
        treefree_work (T->root);
        CC_IFFREE (T->supply, CC_GHnode);
        CC_IFFREE (T->listspace, int);
    }
}

static void treefree_work (CC_GHnode *n)
{
    CC_GHnode *p, *pnext;

    if (n) {
        for (p = n->child; p; p = pnext) {
            pnext = p->sibling;
            treefree_work (p);
        }
    }
}

static void cuttree_free_work (cuttree_node *n, CCptrworld *nodeptr_world)
{
    cuttree_node *p, *pnext;

    if (n) {
        for (p = n->child; p; p = pnext) {
            pnext = p->sibling;
            cuttree_free_work (p, nodeptr_world);
        }
        freenodeset (&n->nlist, nodeptr_world);
    }
}

static int countdescendants (cuttree_node *n)
{
    cuttree_node *p;
    int i;

    i = 1;
    for (p = n->child; p; p = p->sibling)
        i += countdescendants (p);
    n->ndescendants = i;
    return i;
}

/* gh_work(n) expands the cuttree_node n. */

static int gh_work (graph *G, cuttree_node *n, nodeset *nlist,
        nodeset *special, cuttree_node *supply, int *supplyhead,
        CCrandstate *rstate)
{
    nodeset a_nlist, b_nlist;
    nodeset a_special, b_special;
    cuttree_node *a_cut, *b_cut;
    cuttree_node *cp, *nextcp;
    nodeptr *p;
    double cutvalue;
    edgeset esave;
    cuttree_node *newcut;
    node newnode;
    node  *anode = (node *) NULL;
    node  *bnode = (node *) NULL;
    node **names = (node **) NULL;
    int    *elist = (int *) NULL;
    double *ecap  = (double *) NULL;
    edgeptr *ep;
    int *cut = (int *) NULL;
    int i, k, ncount, ecount, cutcount;
    int rval = 0;

    /* termination cases for the recursion */
    if (!special->head->next) {
        ++(G->magicnum);
        if (n->parent) {
            n->parent->pseudonode->magiclabel = G->magicnum;
        }
        for (cp = n->child; cp; cp = cp->sibling) {
            cp->pseudonode->magiclabel = G->magicnum;
        }
        n->special = special->head->this;
        n->nlist.head = n->nlist.tail = (nodeptr *) NULL;
        for (p = nlist->head; p; p = p->next) {
            if (p->this->magiclabel != G->magicnum) {
                rval = addtonodeset (&n->nlist, p->this, &G->nodeptr_world);
                if (rval) {
                    fprintf (stderr, "addtonodeset failed\n");
                    goto CLEANUP;
                }
            }
        }
        return 0;
    }
    /* find the split */

    for (ncount = 0, p = nlist->head; p; p = p->next, ncount++) {
        p->this->num = ncount;
    }

    names = CC_SAFE_MALLOC (ncount, node *);
    if (!names) {
        fprintf (stderr, "out of memory in gh_work\n");
        rval = 1;  goto CLEANUP;
    }
    for (i = 0, p = nlist->head; p; p = p->next, i++) {
        names[i] = p->this;
    }
    for (ecount = 0, p = nlist->head; p; p = p->next) {
        for (ep = p->this->adj.head; ep; ep = ep->next) {
            ep->this->magiclabel = 0;
            ecount++;
        }
    }
    ecount /= 2;

    elist = CC_SAFE_MALLOC (2*ecount, int);
    ecap  = CC_SAFE_MALLOC (ecount, double);
    if (!elist || !ecap) {
        fprintf (stderr, "out of memory in gh_work\n");
        rval = 1; goto CLEANUP;
    }

    for (k = 0, p = nlist->head; p; p = p->next) {
        for (ep = p->this->adj.head; ep; ep = ep->next) {
            if (ep->this->magiclabel == 0) {
                elist[2*k]     = ep->this->ends[0]->num;
                elist[2*k + 1] = ep->this->ends[1]->num;
                ecap[k]        = ep->this->x;
                k++;
                ep->this->magiclabel = 1;
            }
        }
    }


    for (p = special->head, i = 0; p; p = p->next) {
        i++;
        if (myrandnum (i, rstate) == 0)
            anode = p->this;
    }
    for (p = special->head, i = 0; p; p = p->next) {
        if (p->this != anode) {
            i++;
            if (myrandnum (i, rstate) == 0)
                bnode = p->this;
        }
    }

    rval = CCcut_mincut_st (ncount, ecount, elist, ecap, anode->num,
                            bnode->num, &cutvalue, &cut, &cutcount);
    if (rval) {
        fprintf (stderr, "CCcut_mincut_st failed\n"); goto CLEANUP;
    }
    CC_IFFREE (elist, int);
    CC_IFFREE (ecap, double);

    /* make the two new cuttree_nodes */

    /* mark the cut */
    G->magicnum++;
    for (i = 0; i < cutcount; i++) {
        names[cut[i]]->magiclabel = G->magicnum;
    }
    CC_IFFREE (cut, int);
    CC_IFFREE (names, node *);

    /* divide them up */
    splitset (nlist, &a_nlist, &b_nlist, G->magicnum);
    splitset (special, &a_special, &b_special, G->magicnum);

    if (!a_special.head) {
        fprintf (stderr, "Yipes! a_special is null\n");
        if (!b_special.head)
            fprintf (stderr, "And so is b_special\n");
        rval = 1; goto CLEANUP;
    } else if (!b_special.head) {
        fprintf (stderr, "Yipes! b_special is null\n");
        rval = 1; goto CLEANUP;
    }
    newcut = &(supply[(*supplyhead)++]);
    newcut->parent = n;
    newcut->sibling = (cuttree_node *) NULL;
    newcut->child = (cuttree_node *) NULL;
    newcut->cutval = cutvalue;
    newcut->nlist.head = newcut->nlist.tail = (nodeptr *) NULL;

    if (n->parent && n->parent->pseudonode->magiclabel == G->magicnum) {
        a_cut = n;
        b_cut = newcut;
    } else {
        a_cut = newcut;
        b_cut = n;
    }

    /* divide the children up between the two sides */
    for (cp = n->child, n->child = (cuttree_node *) NULL,
                   newcut->child = (cuttree_node *) NULL; cp; cp = nextcp) {
        nextcp = cp->sibling;
        if (cp->pseudonode->magiclabel == G->magicnum) {
            cp->sibling = a_cut->child;
            a_cut->child = cp;
        } else {
            cp->sibling = b_cut->child;
            b_cut->child = cp;
        }
    }
    newcut->sibling = n->child;
    n->child = newcut;

    newnode.magiclabel = 0;
    rval = shrinkdown (G, &a_nlist, &newnode, &esave, G->magicnum);
    if (rval) {
        fprintf (stderr, "shrinkdown failed\n");
        goto CLEANUP;
    }

    a_cut->pseudonode = (node *) NULL;
    b_cut->pseudonode = &newnode;

    rval = addtonodeset (&a_nlist, &newnode, &G->nodeptr_world);
    if (rval) {
        fprintf (stderr, "addtonodeset failed\n");
        goto CLEANUP;
    }
    rval = gh_work (G, a_cut, &a_nlist, &a_special, supply, supplyhead,
                    rstate);
    if (rval) goto CLEANUP;
    delfromnodeset (&a_nlist, &newnode, &G->nodeptr_world);

    ++(G->magicnum);
    for (p = b_nlist.head; p; p = p->next)
        p->this->magiclabel = G->magicnum;

    unshrink (G, &newnode, &esave, G->magicnum);

    rval = shrinkdown (G, &b_nlist, &newnode, &esave, G->magicnum);
    if (rval) {
        fprintf (stderr, "shrinkdown failed\n");
        goto CLEANUP;
    }

    a_cut->pseudonode = &newnode;
    b_cut->pseudonode = (node *) NULL;

    rval = addtonodeset (&b_nlist, &newnode, &G->nodeptr_world);
    if (rval) {
        fprintf (stderr, "addtonodeset failed\n");
        goto CLEANUP;
    }
    rval = gh_work (G, b_cut, &b_nlist, &b_special, supply, supplyhead,
                    rstate);
    if (rval) goto CLEANUP;
    delfromnodeset (&b_nlist, &newnode, &G->nodeptr_world);

    ++(G->magicnum);
    for (p = a_nlist.head; p; p = p->next)
        p->this->magiclabel = G->magicnum;

    unshrink (G, &newnode, &esave, G->magicnum);

    nlist->head = a_nlist.head;
    nlist->tail = a_nlist.tail;
    mergeset (nlist, &b_nlist);

    special->head = a_special.head;
    special->tail = a_special.tail;
    mergeset (special, &b_special);

    rval = 0;

CLEANUP:

    CC_IFFREE (elist, int);
    CC_IFFREE (ecap, double);
    CC_IFFREE (cut, int);
    CC_IFFREE (names, node *);

    return rval;
}

static void splitset (nodeset *s, nodeset *a, nodeset *b, int n)
{
    nodeptr *p;
    nodeptr *pnext;

    a->head = a->tail = (nodeptr *) NULL;
    b->head = b->tail = (nodeptr *) NULL;
    for (p = s->head; p; p = pnext) {
        pnext = p->next;
        if (p->this->magiclabel == n) {
            p->next = a->head;
            a->head = p;
            if (!a->tail)
                a->tail = p;
        } else {
            p->next = b->head;
            b->head = p;
            if (!b->tail)
                b->tail = p;
        }
    }
}

static void mergeset (nodeset *s, nodeset *a)
{
    nodeptr *p;
    nodeptr *pnext;

    for (p = a->head; p; p = pnext) {
        pnext = p->next;
        p->next = s->head;
        s->head = p;
        if (!s->tail)
            s->tail = p;
    }
}

/* shrink everything with G->magicnum != num to pseudo. a is a list of
 * everything with G->magicnum == num. */

static int shrinkdown (graph *G, nodeset *a, node *pseudo, edgeset *esave,
        int num)
{
    edgeset enew;
    nodeptr *p;
    edgeptr *e, *enext;
    node *x;
    double cumx;

    pseudo->adj.head = pseudo->adj.tail = (edgeptr *) NULL;
    esave->head = esave->tail = (edgeptr *) NULL;
    for (p = a->head; p; p = p->next) {
        cumx = 0.0;
        enew.head = enew.tail = (edgeptr *) NULL;
        for (e = p->this->adj.head; e; e = enext) {
            enext = e->next;
            x = OTHEREND (e->this, p->this);
            if (x->magiclabel != num) {
                cumx += e->this->x;
                e->next = esave->head;
                esave->head = e;
                if (!esave->tail)
                    esave->tail = e;
            } else {
                e->next = enew.head;
                enew.head = e;
                if (!enew.tail)
                    enew.tail = e;
            }
        }
        if (cumx > 0.0) {
            e = edgeptralloc (&G->edgeptr_world);
            if (e == (edgeptr *) NULL) {
                fprintf (stderr, "edgeptralloc failed\n");
                return 1;
            }
            e->this = edgealloc (&G->edge_world);
            if (e->this == (edge *) NULL) {
                fprintf (stderr, "edgealloc failed\n");
                edgeptrfree (&G->edgeptr_world, e);
                return 1;
            }
            e->next = enew.head;
            enew.head = e;
            if (!enew.tail)
                enew.tail = e;
            e->this->ends[0] = p->this;
            e->this->ends[1] = pseudo;
            e->this->x = cumx;
            e = edgeptralloc (&G->edgeptr_world);
            if (e == (edgeptr *) NULL) {
                fprintf (stderr, "edgealloc failed\n");
                return 1;
            }
            e->this = enew.head->this;
            e->next = pseudo->adj.head;
            pseudo->adj.head = e;
            if (!pseudo->adj.tail)
                pseudo->adj.tail = e;
        }
        p->this->adj.head = enew.head;
        p->this->adj.tail = enew.tail;
    }
    return 0;
}

/* unshrink pseudo.  everything inside pseudo has been marked with num */

static void unshrink (graph *G, node *pseudo, edgeset *esave, int num)
{
    edgeptr *e, *enext;
    node *x;

    /* take out the new edges we added */
    for (e = pseudo->adj.head; e; e = enext) {
        enext = e->next;
        x = OTHEREND (e->this, pseudo);
        delfromedgeset (&x->adj, e->this, &G->edgeptr_world);
        edgefree (&G->edge_world, e->this);
        edgeptrfree (&G->edgeptr_world, e);
    }

    /* put back the edges we deleted */
    for (e = esave->head; e; e = enext) {
        enext = e->next;
        x = e->this->ends[0];
        if (x->magiclabel == num)
            x = e->this->ends[1];
        e->next = x->adj.head;
        x->adj.head = e;
        if (!x->adj.tail)
            x->adj.tail = e;
    }
}

static int myrandnum (int n, CCrandstate *rstate)
{
    return (CCutil_lprand (rstate) % n);
}

static void set_parents (cuttree_node *n)
{
    cuttree_node *c; 

    for (c = n->child; c; c = c->sibling) {
        c->parent = n;
        set_parents (c);
    }
}

static int copy_cuttree (cuttree_node *root, int ncount, int markcount,
        CC_GHtree *T)
{
    int rval = 0;
    int k;
    int *listspace;

    T->supply =    CC_SAFE_MALLOC (markcount + 1, CC_GHnode);
    T->listspace = CC_SAFE_MALLOC (ncount + 1, int);
    if (!T->supply || !T->listspace) {
        fprintf (stderr, "out of memory in copy_cuttree\n");
        rval = 1; goto CLEANUP;
    }
    k = 0;
    cut_dfs (root, &k);
    k = 0;
    listspace = T->listspace;
    copy_dfs (root, &k, T->supply, &listspace);
    T->root = &T->supply[0];

CLEANUP:

    if (rval) {
        CC_IFFREE (T->supply, CC_GHnode);
        CC_IFFREE (T->listspace, int);
    }
    return rval;
}

static void cut_dfs (cuttree_node *n, int *k)
{
    n->num = (*k)++;
    for (n = n->child; n; n = n->sibling) {
        cut_dfs (n, k);
    }
}

static void copy_dfs (cuttree_node *n, int *k, CC_GHnode *supply,
        int **listspace)
{
    CC_GHnode *c = &supply[*k];
    nodeptr *p;
    int count = 0;

    if (n->parent)  c->parent  = &supply[n->parent->num];
    else            c->parent  = (CC_GHnode *) NULL;
    if (n->sibling) c->sibling = &supply[n->sibling->num];
    else            c->sibling = (CC_GHnode *) NULL;
    if (n->child)   c->child   = &supply[n->child->num];
    else            c->child   = (CC_GHnode *) NULL;
    c->cutval  = n->cutval;
    c->ndescendants = n->ndescendants;
    c->special = n->special->number;
    c->num = *k;

    c->listcount = 0;
    for (p = n->nlist.head; p; p = p->next) {
        c->listcount++;
    }
    c->nlist = *listspace;
    for (p = n->nlist.head; p; p = p->next) {
        c->nlist[count++] = p->this->number;
    }
    *listspace = (*listspace) + c->listcount;

    (*k)++;
    for (n = n->child; n; n = n->sibling) {
        copy_dfs (n, k, supply, listspace);
    }
}

static int buildgraph (graph *G, int ncount, int ecount, int *elist,
        double *ecap, int markcount, int *marks)
{
    edge *e;
    edgeptr *e1;
    int i, k, n1, n2;
    int rval = 0;

    G->edgelist = (edge *) NULL;
    G->nodelist = (node *) NULL;
    G->magicnum = 0;
    G->ncount = 0;
    G->ecount = 0;

    G->nodelist = CC_SAFE_MALLOC (ncount, node);
    G->edgelist = CC_SAFE_MALLOC (ecount, edge);

    if (!G->nodelist || !G->edgelist) {
        fprintf (stderr, "out of memory in buildgraph\n");
        rval = 1; goto CLEANUP;
    }
    G->ncount = ncount;
    G->ecount = ecount;

    for (i = 0; i < ncount; i++) {
        G->nodelist[i].adj.head = G->nodelist[i].adj.tail = (edgeptr *) NULL;
        G->nodelist[i].magiclabel = 0;
        G->nodelist[i].mark = 0;
        G->nodelist[i].number = i;
    }
    for (i = 0, e = G->edgelist, k = 0; i < ecount; i++, e++) {
        n1 = elist[k++];
        n2 = elist[k++];
        e->ends[0] = G->nodelist + n1;
        e->ends[1] = G->nodelist + n2;
        e->x = ecap[i];
        e->magiclabel = 0;
    }

    for (i = ecount, e = G->edgelist; i; i--, e++) {
        e1 = edgeptralloc (&G->edgeptr_world);
        if (!e1) {
            fprintf (stderr, "out of memory in buildgraph\n");
            rval = 1; goto CLEANUP;
        }
        e1->next = e->ends[0]->adj.head;
        e1->this = e;
        e->ends[0]->adj.head = e1;
        if (e->ends[0]->adj.tail == (edgeptr *) NULL) {
            e->ends[0]->adj.tail = e1;
        }
        e1 = edgeptralloc (&G->edgeptr_world);
        if (!e1) {
            fprintf (stderr, "out of memory in buildgraph\n");
            rval = 1; goto CLEANUP;
        }
        e1->next = e->ends[1]->adj.head;
        e1->this = e;
        e->ends[1]->adj.head = e1;
        if (e->ends[1]->adj.tail == (edgeptr *) NULL) {
            e->ends[1]->adj.tail = e1;
        }
    }
    if (markcount) {
        for (i = 0; i < markcount; i++) {
            G->nodelist[marks[i]].mark = 1;
        }
    } else {
        for (i = 0; i < ncount; i++) {
            G->nodelist[i].mark = 1;
        }
    }

CLEANUP:

    if (rval) freegraph (G);
    return rval;
}

static void initnodeset (nodeset *n)
{
    if (n) {
        n->head = (nodeptr *) NULL;
        n->tail = (nodeptr *) NULL;
    }
}

static void freenodeset (nodeset *s, CCptrworld *nodeptr_world)
{
    if (s) {
        nodeptr_listfree (nodeptr_world, s->head);
        s->head = (nodeptr *) NULL;
        s->tail = (nodeptr *) NULL;
    }
}

static void initgraph (graph *G)
{
    if (G) {
        G->nodelist = (node *) NULL;
        G->edgelist = (edge *) NULL;
        G->ncount   = 0;
        G->ecount   = 0;
        G->magicnum = 0;
        CCptrworld_init (&G->edge_world);
        CCptrworld_init (&G->edgeptr_world);
        CCptrworld_init (&G->nodeptr_world);
    }
}


static void freegraph (graph *G)
{
    int i;
    node *n;

    if (G) {
        if (G->nodelist) {
            for (i = G->ncount, n = G->nodelist; i; i--, n++) {
                edgeptr_listfree (&G->edgeptr_world, n->adj.head);
                n->adj.head = n->adj.tail = (edgeptr *) NULL;
            }
            CC_FREE (G->nodelist, node);
        }
        CC_IFFREE (G->edgelist, edge);
    }
}

void CCcut_GHtreeprint (CC_GHtree *T)
{
    if (T) {
        printf ("GOMORY-HU TREE\n"); fflush (stdout);
        print_tree_work (T->root);
    }
}

static void print_tree_work (CC_GHnode *n)
{
    int i;

    printf ("T%d: ", n->num);
    printf ("Set (");
    for (i = 0; i < n->listcount; i++) {
        printf ("%d,", n->nlist[i]);
    }
    printf ("[%d]), ", n->special);
    if (n->parent) printf ("Parent %d, ", n->parent->num);
    else           printf ("Parent NULL, ");
    if (n->child)  printf ("Child %d, ", n->child->num);
    else           printf ("Child NULL, ");
    if (n->sibling) printf ("Sibling %d, ", n->sibling->num);
    else            printf ("Sibling NULL, ");
    if (n->parent) {
        printf ("Cnt %d, Val %.2f\n", n->ndescendants, n->cutval);
        fflush (stdout);
    } else {
        printf ("Cnt %d, ROOT\n", n->ndescendants); 
        fflush (stdout);
    }

    for (n = n->child; n; n = n->sibling) {
        print_tree_work (n);
    }
}

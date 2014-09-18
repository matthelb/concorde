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
/*  int CCverify_cut (CCtsp_lpcut_in *cut, int check_types, int *type)      */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCverify_classify (CCtsp_lpcut_in *cut, int check_types,            */
/*      CCverify_cutclass *class);                                          */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCverify_initcutclass (CCverify_cutclass *class)                   */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCverify_freecutclass (CCverify_cutclass *class)                   */
/*    MISSING                                                               */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "tsp.h"
#include "util.h"
#include "verify.h"
#include "macrorus.h"
#include "heldkarp.h"
/*#include "tinytsp.h"*/

#undef  DEBUG
#undef  DUMPSMALL
#define DUMPFAIL

#define HELDKARP_LIMIT 20000000

typedef struct atom {
    int mark;
    int representative;
    int cliquecount;
    int *cliquelist;
} atom;

typedef struct vclique {
    int mark;
    int atomcount;
    int *atomlist;
    int crosscount;
    int *crosslist;
    int color;
    int flipped;
    struct vclique *parent;
    struct vclique *child;
    struct vclique *sibling;
    struct vclique **prevsib;
} vclique;

#define CHILD_ADD(n,c) {                                       \
    (c)->sibling = (n)->child;                                 \
    if ((c)->sibling) (c)->sibling->prevsib = &((c)->sibling); \
    (c)->prevsib = &((n)->child);                              \
    (n)->child = (c);                                          \
    (c)->parent = (n);                                         \
}

#define CHILD_DEL(c) {                                         \
    *((c)->prevsib) = (c)->sibling;                            \
    if ((c)->sibling) (c)->sibling->prevsib = (c)->prevsib;    \
}

typedef struct node {
    int atom;
} node;

typedef struct graph {
    int nodecount;
    int edgecount;
    int *elist;
    int *elen;
} graph;

typedef struct atom_info {
    int atomcount;
    int cliquecount;
    int nodespace;
    char sense;
    int rhs;
    atom *atomlist;
    node *nodelist;
    vclique *cliquelist;
    vclique family[2]; /* dummy roots for the nested families of cliques */
} atom_info;

static int
    build_atom_info (CCtsp_lpcut_in *cut, atom_info *atoms),
    build_atoms (CCtsp_lpcut_in *cut, atom_info *atoms),
    build_clique_atoms (CCtsp_lpclique *clique, atom_info *atoms, int cnum),
    build_atom_cliques (atom_info *atoms),
    verify_atom_info (CCtsp_lpcut_in *cut, atom_info *atoms),
    verify_nodes (CCtsp_lpcut_in *cut, atom_info *atoms),
    verify_atom_clique (CCtsp_lpclique *clique, atom_info *atoms, int cnum),
    verify_atom_lists (atom_info *atoms),
    find_clique_atom (atom_info *atoms, int cnum, int anum),
    find_atom_clique (atom_info *atoms, int anum, int cnum),
    verify_subtour (atom_info *atoms, CCverify_cutclass *class),
    atom_singleton (atom *a, vclique *cliquelist, int c, int nflipped),
    atom_doubleton (atom *a, vclique *cliquelist, int c1, int c2,
        int nflipped),
    verify_chvatal_comb (atom_info *atoms, CCverify_cutclass *class),
    verify_tooth (atom_info *atoms, int tooth, int handle, int nflipped),
    build_binesting (atom_info *atoms),
    build_cross_cliques (atom_info *atoms, int cliq),
    binesting_dfs (atom_info *atoms, int cliq, int color),
    nest_atom_color (atom_info *atoms, int color),
    nest_clique (atom_info *atoms, int cliqnum, vclique *family),
    clique_complement (atom_info *atoms, vclique *cliq, int *catomcount,
        int **catomlist),
    verify_star (atom_info *atoms, int handnum, CCverify_cutclass *class),
    verify_bipartition (atom_info *atoms, int handnum),
    flip_disjoint (vclique *family),
    flip_disjoint_one (vclique *family),
    flip_concentric (vclique *family),
    child_count (vclique *root),
    family_count (vclique *root, int *nflipped),
    check_mark_concentric (vclique *root, int mark),
    collect_atom_mark_counts (atom_info *atoms, vclique *c, int *alist,
        int *p_acnt, int nflipped),
    verify_other (atom_info *atoms),
    build_complete_graph (int nodecount, graph **p_g),
    verify_rhs (atom_info *atoms, graph *g),
    tooth_has_cavity (atom_info *atoms, vclique *c, int nflipped),
    build_cutclass (CCverify_cutclass *class, int ncliques, int nfamilies);

static void
    free_atom_info (atom_info *atoms),
    renest_family (vclique *family),
    flip_marked_cliques (vclique *family),
    flip_clique (vclique *family, vclique *cliq),
    free_graph (graph **p_g),
    compute_lhs (atom_info *atoms, graph *g),
    cleanup_atom_info (CCtsp_lpcut_in *cut, atom_info *atoms),
    clear_atom_marks (atom_info *atoms),
    clear_clique_marks (atom_info *atoms);

static vclique
   *concentric_innermost (vclique *root);

int CCverify_cut (CCtsp_lpcut_in *cut, int check_types, int *type)
{
    atom_info atoms;
    int rval;

#ifdef DEBUG
    double sz = CCutil_zeit();
    printf ("Verifying cut, %d cliques %d atoms",cut->cliquecount,
            cut->skel.atomcount);
    fflush (stdout);
#endif /* DEBUG */

    if (type != (int *) NULL) *type = -1;
    
    if (cut->skel.atomcount == 0 || cut->skel.atoms == (int *) NULL) {
        fprintf (stderr, "Cut to verify has no skeleton\n");
        rval = 1; goto CLEANUP;
    }
    
    rval = build_atom_info (cut, &atoms);
    if (rval) {
        fprintf (stderr, "build_atom_info failed\n");
        goto CLEANUP;
    }

    rval = verify_atom_info (cut, &atoms);
    if (rval) {
        fprintf (stderr, "atom_info failed verification\n");
        goto CLEANUP;
    }

    if (check_types & CC_TYPE_SUBTOUR) {
        rval = verify_subtour (&atoms, (CCverify_cutclass *) NULL);
        if (rval == 0) {
            /* It's a valid subtour - we're done */
#ifdef DEBUG
            printf (" (subtour)");
            fflush (stdout);
#endif
            if (type != (int *) NULL) *type = CC_TYPE_SUBTOUR;
            goto CLEANUP;
        }
    }

    if (check_types & CC_TYPE_COMB) {
        rval = verify_chvatal_comb (&atoms, (CCverify_cutclass *) NULL);
        if (rval == 0) {
            /* It's a valid comb - we're done */
#ifdef DEBUG
            printf (" (comb 1,%d)", atoms.cliquecount-1);
            fflush (stdout);
#endif
            if (type != (int *) NULL) *type = CC_TYPE_COMB;
            goto CLEANUP;
        }
    }

    if (check_types & (CC_TYPE_STAR | CC_TYPE_BIPARTITION)) {
        rval = build_binesting (&atoms);
        if (rval == 0) {
            if (check_types & CC_TYPE_STAR) {
                rval = verify_star (&atoms, 0, (CCverify_cutclass *) NULL);
                if (rval == 0) {
                    /* It's a valid star, with handles marked 0 - we're done */
#ifdef DEBUG
                    int hcnt = family_count (&atoms.family[0],
                                             (int *) NULL) - 1;
                    printf (" (star %d,%d)", hcnt, atoms.cliquecount - hcnt);
                    fflush (stdout);
#endif
                    if (type != (int *) NULL) *type = CC_TYPE_STAR;
                    goto CLEANUP;
                }
                
                rval = verify_star (&atoms, 1, (CCverify_cutclass *) NULL);
                if (rval == 0) {
                    /* It's a valid star, with handles marked 1 - we're done */
#ifdef DEBUG
                    int hcnt = family_count (&atoms.family[1],
                                             (int *) NULL) - 1;
                    printf (" (star %d,%d)", hcnt, atoms.cliquecount - hcnt);
                    fflush (stdout);
#endif
                    if (type != (int *) NULL) *type = CC_TYPE_STAR;
                    goto CLEANUP;
                }
            }
            if (check_types & CC_TYPE_BIPARTITION) {
                rval = verify_bipartition (&atoms, 0);
                if (rval == 0) {
                    /* It's a valid bipartition, handles marked 0 - done */
#ifdef DEBUG
                    int hcnt = family_count (&atoms.family[0],
                                             (int *) NULL) - 1;
                    printf (" (bipartition %d,%d)", hcnt,
                            atoms.cliquecount - hcnt);
                    fflush (stdout);
#endif
                    if (type != (int *) NULL) *type = CC_TYPE_BIPARTITION;
                    goto CLEANUP;
                }
                
                rval = verify_bipartition (&atoms, 1);
                if (rval == 0) {
                    /* It's a valid bipartition, handles marked 1 - done */
#ifdef DEBUG
                    int hcnt = family_count (&atoms.family[1],
                                             (int *) NULL) - 1;
                    printf (" (bipartition %d,%d)", hcnt,
                            atoms.cliquecount - hcnt);
                    fflush (stdout);
#endif
                    if (type != (int *) NULL) *type = CC_TYPE_BIPARTITION;
                    goto CLEANUP;
                }
            }
            
#if 0
#ifdef DEBUG
            printf ("\nCurious, binested, but not star or bipartition:\n");
            CCtsp_print_lpcut_in (cut);
            fflush (stdout);
#endif
#endif
        }
    }

    if (check_types & CC_TYPE_OTHER) {
        rval = verify_other (&atoms);
        if (rval == 0) {
#ifdef DEBUG
            printf (" (other)");
            fflush (stdout);
#endif
            if (type != (int *) NULL) *type = CC_TYPE_OTHER;
            goto CLEANUP;
        }
        fprintf (stderr, "Unable to verify cut\n");
    }

    rval = -1;

  CLEANUP:
#ifdef DEBUG
    if (rval == 0) {
        printf (" in %.2f seconds\n", CCutil_zeit() - sz);
        fflush (stdout);
    } else {
        printf (" FAILED in %.2f seconds\n", CCutil_zeit() - sz);
        printf ("FAILED CUT:\n");
        CCtsp_print_lpcut_in (cut);
        fflush (stdout);
    }
#endif
    free_atom_info (&atoms);
    return rval;
}

static int build_atom_info (CCtsp_lpcut_in *cut, atom_info *atoms)
{
    int rval;
    int i;
    CCtsp_skeleton *skel = &cut->skel;
    int cliquecount = cut->cliquecount;
    int maxn;

    atoms->atomlist = (atom *) NULL;
    atoms->atomcount = 0;
    atoms->nodelist = (node *) NULL;
    atoms->cliquelist = (vclique *) NULL;

    maxn = CCtsp_max_node (cut);
    for (i=0; i<skel->atomcount; i++) {
        if (skel->atoms[i] > maxn) maxn = skel->atoms[i];
    }
    atoms->nodespace = maxn + 1;

    atoms->nodelist = CC_SAFE_MALLOC (atoms->nodespace, node);
    if (atoms->nodelist == (node *) NULL) {
        fprintf (stderr, "Out of memory in build_atom_info\n");
        goto FAILURE;
    }

    cleanup_atom_info (cut, atoms);

    rval = build_atoms (cut, atoms);
    if (rval) {
        fprintf (stderr, "build_atoms failed\n");
        goto FAILURE;
    }
    
    atoms->cliquecount = cliquecount;
    atoms->sense = cut->sense;
    atoms->rhs = cut->rhs;

    atoms->cliquelist = CC_SAFE_MALLOC (cliquecount, vclique);
    if (atoms->cliquelist == (vclique *) NULL) {
        fprintf (stderr, "Out of memory in build_atom_info\n");
        rval = -1;
        goto FAILURE;
    }
    for (i=0; i<cliquecount; i++) {
        atoms->cliquelist[i].atomlist = (int *) NULL;
        atoms->cliquelist[i].atomcount = 0;
        atoms->cliquelist[i].crosslist = (int *) NULL;
        atoms->cliquelist[i].crosscount = 0;
    }

    for (i=0; i<cliquecount; i++) {
        rval = build_clique_atoms (&cut->cliques[i], atoms, i);
        if (rval) {
            fprintf (stderr, "build_clique_atoms failed\n");
            goto FAILURE;
        }
    }

    rval = build_atom_cliques (atoms);
    if (rval) {
        fprintf (stderr, "build_atom_cliques failed\n");
        goto FAILURE;
    }

    return 0;

  FAILURE:
    free_atom_info (atoms);
    return -1;
}

static int build_atoms (CCtsp_lpcut_in *cut, atom_info *atoms)
{
    node *nodelist = atoms->nodelist;
    CCtsp_skeleton *skel = &cut->skel;
    atom *atomlist;
    int atomcount;
    int i;

    atomcount = skel->atomcount;
    atomlist = CC_SAFE_MALLOC (atomcount, atom);
    if (atomlist == (atom *) NULL) {
        goto FAILURE;
    }

    atoms->atomcount = atomcount;
    atoms->atomlist = atomlist;

    for (i=0; i<atomcount; i++) {
        atomlist[i].mark = 0;
        atomlist[i].representative = skel->atoms[i];
        atomlist[i].cliquecount = 0;
        atomlist[i].cliquelist = (int *) NULL;
        nodelist[atomlist[i].representative].atom = i;
    }

    return 0;

  FAILURE:
    CC_IFFREE (atoms->atomlist, atom);
    atoms->atomcount = 0;
    return -1;
}

static int build_clique_atoms (CCtsp_lpclique *clique, atom_info *atoms,
        int cnum)
{
    node *nodelist = atoms->nodelist;
    atom *atomlist = atoms->atomlist;
    vclique *newclique = &(atoms->cliquelist[cnum]);
    int i;
    int j;
    int cnt;

    clear_atom_marks (atoms);

    cnt = 0;
    CC_FOREACH_NODE_IN_CLIQUE (i, *clique, j) {
        if (nodelist[i].atom >= 0) {
            cnt++;
            atomlist[nodelist[i].atom].mark = 1;
        }
    }
    
    newclique->atomcount = cnt;
    newclique->atomlist = CC_SAFE_MALLOC (cnt, int);
    if (newclique->atomlist == (int *) NULL) {
        fprintf (stderr, "Out of memory in build_clique_atoms\n");
        return -1;
    }
    for (i=0, cnt=0; i<atoms->atomcount; i++) {
        if (atomlist[i].mark) {
            newclique->atomlist[cnt++] = i;
        }
    }
    if (cnt != newclique->atomcount) {
        fprintf (stderr, "SURPRISE - lost some atoms\n");
        CC_IFFREE (newclique->atomlist, int);
        return -1;
    }

    return 0;
}

static int build_atom_cliques (atom_info *atoms)
{
    atom *atomlist = atoms->atomlist;
    int atomcount = atoms->atomcount;
    vclique *cliquelist = atoms->cliquelist;
    int cliquecount = atoms->cliquecount;
    int i;
    int j;
    int a;

    for (i=0; i<atomcount; i++) {
        atomlist[i].cliquecount = 0;
    }

    for (i=0; i<cliquecount; i++) {
        for (j=0; j<cliquelist[i].atomcount; j++) {
            atomlist[cliquelist[i].atomlist[j]].cliquecount++;
        }
    }

    for (i=0; i<atomcount; i++) {
        if (atomlist[i].cliquecount > 0) {
            atomlist[i].cliquelist = CC_SAFE_MALLOC (atomlist[i].cliquecount,
                                                     int);
            if (atomlist[i].cliquelist == (int *) NULL) {
                fprintf (stderr, "Out of memory in build_atom_cliques\n");
                goto FAILURE;
            }
        } else {
            atomlist[i].cliquelist = (int *) NULL;
        }
        atomlist[i].cliquecount = 0;
    }

    for (i=0; i<cliquecount; i++) {
        for (j=0; j<cliquelist[i].atomcount; j++) {
            a = cliquelist[i].atomlist[j];
            atomlist[a].cliquelist[atomlist[a].cliquecount++] = i;
        }
    }

    return 0;

  FAILURE:
    for (i=0; i<atomcount; i++) {
        CC_IFFREE (atomlist[i].cliquelist, int);
    }
    return -1;
}

static void free_atom_info (atom_info *atoms)
{
    int i;

    if (atoms->atomlist) {
        for (i=0; i<atoms->atomcount; i++) {
            CC_IFFREE (atoms->atomlist[i].cliquelist, int);
        }
        CC_FREE (atoms->atomlist, atom);
    }
    if (atoms->cliquelist) {
        for (i=0; i<atoms->cliquecount; i++) {
            CC_IFFREE (atoms->cliquelist[i].atomlist, int);
            CC_IFFREE (atoms->cliquelist[i].crosslist, int);
        }
        CC_FREE (atoms->cliquelist, vclique);
    }

    CC_IFFREE (atoms->nodelist, node);

    atoms->nodespace = 0;
    atoms->atomcount = 0;
    atoms->cliquecount = 0;
}

/*    verify_atom_info verifies that the atom information in atoms
 *    corresponds to a subset of the atoms of cut.  In particular, the atom
 *    structure has the following properties:
 * (0) Every node of every clique is < nodespace (the size of nodelist)
 * (1) Every atom has a representative node (trivial by listing the
 *     representative for each atom)
 * (2) Every node that is a representative (ie, atom >= 0) is the
 *     representative for that atom.
 * (3) The representative node for every atom listed in a clique's
 *     atomlist is inside the clique
 * (4) The representative node for every atom not listed in a clique's
 *     atomlist is outside the clique
 * (5) Every atom listed in a clique's atomlist lists the clique in its
 *     cliquelist
 * (6) Every clique listed in an atom's cliquelist lists the atom in its
 *     atomlist
 * (7) The rhs, sense, and cliquecount are the same in the cut and the atom
 *
 */
static int verify_atom_info (CCtsp_lpcut_in *cut, atom_info *atoms)
{
    int i;
    int rval = 0;

/* (0) Every node of every clique is < nodespace (the size of nodelist)
 * (2) Every node that is a representative (ie, atom >= 0) is the
 *     representative for that atom.
 */
    rval = verify_nodes (cut, atoms);
    if (rval) goto CLEANUP;

/* (3) The representative node for every atom listed in a clique's
 *     atomlist is inside the clique
 * (4) The representative node for every atom not listed in a clique's
 *     atomlist is outside the clique
 */
    for (i=0; i<cut->cliquecount; i++) {
        rval = verify_atom_clique (&(cut->cliques[i]), atoms, i);
        if (rval) goto CLEANUP;
    }

/* (5) Every atom listed in a clique's atomlist lists the clique in its
 *     cliquelist
 * (6) Every clique listed in an atom's cliquelist lists the atom in its
 *     atomlist
 */
    rval = verify_atom_lists (atoms);
    if (rval) goto CLEANUP;

/* (7) The rhs, sense, and cliquecount are the same in the cut and the atom
 */
    if (cut->rhs != atoms->rhs) {
        fprintf (stderr, "cut rhs %d != atom rhs %d\n", cut->rhs,
                 atoms->rhs);
        rval = -1;
        goto CLEANUP;
    }
    if (cut->sense != atoms->sense) {
        fprintf (stderr, "cut sense %c != atom sense %c\n", cut->sense,
                 atoms->sense);
        rval = -1;
        goto CLEANUP;
    }
    if (cut->cliquecount != atoms->cliquecount) {
        fprintf (stderr, "cut cliquecount %d != atom cliquecount %d\n",
                 cut->cliquecount, atoms->cliquecount);
        rval = -1;
        goto CLEANUP;
    }

    rval = 0;

  CLEANUP:
    return rval;
}

/* (0) Every node of every clique is < nodespace (the size of nodelist)
 * (2) Every node that is a representative (ie, atom >= 0) is the
 *     representative for that atom.
 */
static int verify_nodes (CCtsp_lpcut_in *cut, atom_info *atoms)
{
    node *nodelist = atoms->nodelist;
    int nodespace = atoms->nodespace;
    CCtsp_lpclique *cliques = cut->cliques;
    int cliquecount = cut->cliquecount;
    int atomcount = atoms->atomcount;
    atom *atomlist = atoms->atomlist;
    int i;
    int k;
    int j;

    for (i=0; i<cliquecount; i++) {
        CC_FOREACH_NODE_IN_CLIQUE (j, cliques[i], k) {
            if (j < 0 || j >= nodespace) {
                fprintf (stderr, "clique %d contains bogus node %d\n",
                         i, j);
                return -1;
            }
            if (nodelist[j].atom < -1 || nodelist[j].atom >= atomcount) {
                fprintf (stderr, "node %d is in bogus atom %d\n", j,
                         nodelist[j].atom);
                return -1;
            }
            if (nodelist[j].atom >= 0 &&
                atomlist[nodelist[j].atom].representative != j) {
                fprintf (stderr, "node %d is not the representative of its atom %d\n",
                         j, nodelist[j].atom);
            }
        }
    }
    return 0;
}

/* (3) The representative node for every atom listed in a clique's
 *     atomlist is inside the clique
 * (4) The representative node for every atom not listed in a clique's
 *     atomlist is outside the clique
 */
static int verify_atom_clique (CCtsp_lpclique *clique, atom_info *atoms,
        int cnum)
{
    atom *atomlist = atoms->atomlist;
    int atomcount = atoms->atomcount;
    node *nodelist = atoms->nodelist;
    vclique *cliqueinfo = &(atoms->cliquelist[cnum]);
    int i;
    int j;
    int a;
    int rval;

    clear_atom_marks (atoms);

    CC_FOREACH_NODE_IN_CLIQUE (i, *clique, j) {
        a = nodelist[i].atom;
        if (a >= 0) {
            if (atomlist[a].mark) {
                fprintf (stderr, "Duplicate atom %d in clique\n", a);
                rval = -1;
                goto CLEANUP;
            }
            atomlist[a].mark = 1;
        }
    }

    for (i=0; i<cliqueinfo->atomcount; i++) {
        j = cliqueinfo->atomlist[i];
        if (atomlist[j].mark != 1) {
            fprintf (stderr, "Atom %d not in clique %d\n", j, cnum);
            rval = -1;
            goto CLEANUP;
        }
        atomlist[j].mark = 0;
    }

    for (i=0; i<atomcount; i++) {
        if (atomlist[i].mark != 0) {
            fprintf (stderr, "Atom %d not in cliquelist %d\n", j, cnum);
            rval = -1;
            goto CLEANUP;
        }
    }
    rval = 0;

  CLEANUP:
    return rval;
}

/* (5) Every atom listed in a clique's atomlist lists the clique in its
 *     cliquelist
 * (6) Every clique listed in an atom's cliquelist lists the atom in its
 *     atomlist
 */
static int verify_atom_lists (atom_info *atoms)
{
    int atomcount = atoms->atomcount;
    atom *atomlist = atoms->atomlist;
    int cliquecount = atoms->cliquecount;
    vclique *cliquelist = atoms->cliquelist;
    int i;
    int j;
    int rval;

    for (i=0; i<atomcount; i++) {
        for (j=0; j<atomlist[i].cliquecount; j++) {
            rval = find_clique_atom (atoms, atomlist[i].cliquelist[j], i);
            if (rval) {
                fprintf (stderr, "clique %d missing atom %d\n",
                         atomlist[i].cliquelist[j], i);
                return rval;
            }
        }
    }

    for (i=0; i<cliquecount; i++) {
        for (j=0; j<cliquelist[i].atomcount; j++) {
            rval = find_atom_clique (atoms, cliquelist[i].atomlist[j], i);
            if (rval) {
                fprintf (stderr, "atom %d missing clique %d\n",
                         cliquelist[i].atomlist[j], i);
                return rval;
            }
        }
    }
    return 0;
}

static int find_clique_atom (atom_info *atoms, int cnum, int anum)
{
    int i;
    vclique *cliquelist = atoms->cliquelist;

    for (i=0; i<cliquelist[cnum].atomcount; i++) {
        if (cliquelist[cnum].atomlist[i] == anum) return 0;
    }
    return 1;
}

static int find_atom_clique (atom_info *atoms, int anum, int cnum)
{
    int i;
    atom *atomlist = atoms->atomlist;

    for (i=0; i<atomlist[anum].cliquecount; i++) {
        if (atomlist[anum].cliquelist[i] == cnum) return 0;
    }
    return 1;
}

static int verify_subtour (atom_info *atoms, CCverify_cutclass *class)
{
    int rval;
    if (atoms->cliquecount == 1 && atoms->atomcount == 2 &&
        atoms->rhs == 2 && atoms->sense == 'G') {
        if (class != (CCverify_cutclass *) NULL) {
            rval = build_cutclass (class, 1, 1);
            if (rval) {
                fprintf (stderr, "build_cutclass failed\n");
                return rval;
            }
            class->type = CC_TYPE_SUBTOUR;
            class->nhandles = 1;
            class->nfamilies = 1;
            class->cliques[0] = 0;
            class->inverted[0] = 0;
            class->family_start[0] = 0;
        }
        return 0;
    } else {
        return -1;
    }
}

static int atom_singleton (atom *a, vclique *cliquelist, int c, int nflipped)
{
    int i;
    int j;
    int ccnt = 0;
    int fcnt = 0;

    if (cliquelist[c].flipped) {
        fcnt++;
        ccnt++;
    }
    
    for (i=0; i<a->cliquecount; i++) {
        j = a->cliquelist[i];
        if (cliquelist[j].flipped == 0) {
            if (j != c) return 0;
            ccnt++;
        } else {
            if (j == c) return 0;
            fcnt++;
        }
    }
    if (fcnt == nflipped && ccnt == 1) return 1;
    else return 0;
}

static int atom_doubleton (atom *a, vclique *cliquelist, int c1, int c2,
        int nflipped)
{
    int i;
    int j;
    int c1cnt = 0;
    int c2cnt = 0;
    int fcnt = 0;

    if (cliquelist[c1].flipped) {
        fcnt++;
        c1cnt++;
    }
    if (cliquelist[c2].flipped) {
        fcnt++;
        c2cnt++;
    }
    for (i=0; i<a->cliquecount; i++) {
        j = a->cliquelist[i];
        if (cliquelist[j].flipped == 0) {
            if      (j == c1) c1cnt++;
            else if (j == c2) c2cnt++;
            else              return 0;
        } else {
            if (j == c1 || j == c2) return 0;
            fcnt++;
        }
    }
    if (c1cnt == 1 && c2cnt == 1 && fcnt == nflipped) return 1;
    else return 0;
}

/* combs come in 4 varieties, depending on whether there is a handle cavity
 * (hc), and whether there is an outside cavity (oc).  Every comb has
 * (1) an even number of cliques (1 handle + odd teeth)
 * (2) at least 4 cliques
 * (3) rhs 3*(ncliques-1) + 1, sense 'G'
 * (4) 2*(ncliques-1) + hc + oc atoms
 * (5) a handle, which contains ncliques-1+hc atoms
 * (6) ncliques-1 teeth, which contain 2 atoms or atomcount-2 atoms
 * If we "invert" the teeth so that each tooth contains 2 atoms, then
 * (7) each tooth contains exactly 2 atoms, one of just the tooth, and one
 *     of just the tooth and the handle
 */
static int verify_chvatal_comb (atom_info *atoms, CCverify_cutclass *class)
{
    int cliquecount = atoms->cliquecount;
    vclique *cliquelist = atoms->cliquelist;
    int atomcount = atoms->atomcount;
    int i;
    int handlelo;
    int handlehi;
    int toothlo;
    int toothhi;
    int handlecount;
    int handle = -1;
    int nflipped;
    int rval;
    int j;

    if (cliquecount % 2 != 0) return 1;
    if (cliquecount < 4) return 1;
    if (atoms->rhs != 3 * cliquecount - 2) return 1;
    if (atoms->sense != 'G') return 1;

    if (atomcount == 2*cliquecount-2) {
        handlelo = cliquecount-1;
        handlehi = cliquecount-1;
    } else if (atomcount == 2*cliquecount-1) {
        handlelo = cliquecount-1;
        handlehi = cliquecount;
    } else if (atomcount == 2*cliquecount) {
        handlelo = cliquecount;
        handlehi = cliquecount;
    } else {
        return 1;
    }
    toothlo = 2;
    toothhi = atomcount-2;

    handlecount = 0;
    nflipped = 0;
    for (i=0; i<cliquecount; i++) {
        if (cliquelist[i].atomcount >= handlelo &&
            cliquelist[i].atomcount <= handlehi) {
            cliquelist[i].flipped = 0;
            handlecount++;
            handle = i;
        } else if (cliquelist[i].atomcount == toothlo) {
            cliquelist[i].flipped = 0;
        } else if (cliquelist[i].atomcount == toothhi) {
            cliquelist[i].flipped = 1;
            nflipped++;
        } else {
            return 1;
        }
    }
    if (handlecount != 1) return 1;

    for (i=0; i<cliquecount; i++) {
        if (i != handle) {
            if (verify_tooth (atoms, i, handle, nflipped)) {
                return 1;
            }
        }
    }

    if (class != (CCverify_cutclass *) NULL) {
        rval = build_cutclass (class, cliquecount, cliquecount);
        if (rval) {
            fprintf (stderr, "build_cutclass failed\n");
            return 1;
        }
        class->type = CC_TYPE_COMB;
        class->nhandles = 1;
        class->nfamilies = cliquecount;
        j=0;
        class->cliques[j] = handle;
        class->inverted[j] = 0;
        class->family_start[j] = 0;
        j++;
        for (i=0; i<cliquecount; i++) {
            if (i != handle) {
                class->cliques[j] = i;
                class->inverted[j] = cliquelist[i].flipped;
                class->family_start[j] = j;
                j++;
            }
        }
    }

    return 0;
}

static int verify_tooth (atom_info *atoms, int tooth, int handle, int nflipped)
{
    int i;
    int atomcount = atoms->atomcount;
    atom *atomlist = atoms->atomlist;
    vclique *cliquelist = atoms->cliquelist;
    int cliqueatomcount = cliquelist[tooth].atomcount;
    int *cliqueatomlist = cliquelist[tooth].atomlist;
    int inside = 0;
    int outside = 0;
    int inmark;

    clear_atom_marks (atoms);

    for (i=0; i<cliqueatomcount; i++) {
        atomlist[cliqueatomlist[i]].mark = 1;
    }

    if (cliquelist[tooth].flipped == 0) {
        inmark = 1;
    } else {
        inmark = 0;
    }

    for (i=0; i<atomcount; i++) {
        if (atomlist[i].mark == inmark) {
            if (atom_singleton (&atomlist[i], cliquelist, tooth, nflipped)) {
                if (outside == 1) return 1;
                outside = 1;
            } else if (atom_doubleton (&atomlist[i], cliquelist, handle,
                                       tooth, nflipped)) {
                if (inside == 1) return 1;
                inside = 1;
            } else {
                return 1;
            }
        }
    }

    if (inside == 1 && outside == 1) return 0;
    else return 1;
}

static int build_binesting (atom_info *atoms)
{
    int cliquecount = atoms->cliquecount;
    vclique *cliquelist = atoms->cliquelist;
    atom *a;
    int i;
    int rval;

    clear_clique_marks (atoms);

    for (i=0; i<cliquecount; i++) {
        rval = build_cross_cliques (atoms, i);
        if (rval) return rval;
    }
    
    for (i=0; i<cliquecount; i++) {
        cliquelist[i].color = -1;
    }

    rval = binesting_dfs (atoms, 0, 0);
    if (rval) return rval;

    for (i=0; i<cliquecount; i++) {
        if (cliquelist[i].color != 0 && cliquelist[i].color != 1) {
            return 1;
        }
    }

    /* every clique is flipped so 0 is outside, guaranteeing that
     * the family nests.
     */
    for (i=0; i<cliquecount; i++) {
        cliquelist[i].flipped = 0;
    }
    a = &atoms->atomlist[0];
    for (i=0; i<a->cliquecount; i++) {
        cliquelist[a->cliquelist[i]].flipped = 1;
    }
    
    rval = nest_atom_color (atoms, 0);
    if (rval) {
        fprintf (stderr, "nest_atom_color failed\n");
        return rval;
    }
    rval = nest_atom_color (atoms, 1);
    if (rval) {
        fprintf (stderr, "nest_atom_color failed\n");
        return rval;
    }
    return 0;
}

static int build_cross_cliques (atom_info *atoms, int cliq)
{
    vclique *cliquelist = atoms->cliquelist;
    vclique *c = &cliquelist[cliq];
    atom *atomlist = atoms->atomlist;
    int atomcount = atoms->atomcount;
    int catomcount = c->atomcount;
    int *catomlist = c->atomlist;
    int crosscount = 0;
    int *crosslist = (int *) NULL;
    atom *a;
    int i;
    int j;
    int k;
    int rval;

    crosslist = CC_SAFE_MALLOC (atoms->cliquecount, int);
    if (crosslist == (int *) NULL) {
        fprintf (stderr, "Out of memory in build_cross_cliques\n");
        rval = 1; goto CLEANUP;
    }
    
    for (i=0; i<catomcount; i++) {
        a = &atomlist[catomlist[i]];
        for (j=0; j<a->cliquecount; j++) {
            k = a->cliquelist[j];
            cliquelist[k].mark++;
        }
    }
    cliquelist[cliq].mark = 0;

    crosscount = 0;
    for (i=0; i<catomcount; i++) {
        a = &atomlist[catomlist[i]];
        for (j=0; j<a->cliquecount; j++) {
            k = a->cliquelist[j];
            if (cliquelist[k].mark > 0 &&
                cliquelist[k].mark < catomcount &&
                cliquelist[k].mark < cliquelist[k].atomcount &&
                catomcount + cliquelist[k].atomcount <
                cliquelist[k].mark + atomcount) {
                crosslist[crosscount++] = k;
            }
            cliquelist[k].mark = 0;
        }
    }

    if (crosscount == 0) {
        printf ("Warning: clique %d crosses no other cliques\n", cliq);
        c->crosslist = (int *) NULL;
    } else {
        c->crosslist = CC_SAFE_MALLOC (crosscount, int);
        if (c->crosslist == (int *) NULL) {
            fprintf (stderr, "Out of memory in build_cross_cliques\n");
            rval = 1; goto CLEANUP;
        }
        for (i=0; i<crosscount; i++) {
            c->crosslist[i] = crosslist[i];
        }
    }
    c->crosscount = crosscount;
    rval = 0;

 CLEANUP:
    CC_IFFREE (crosslist, int);
    return rval;
}

static int binesting_dfs (atom_info *atoms, int cliq, int color)
{
    vclique *c = &atoms->cliquelist[cliq];
    int crosscount = c->crosscount;
    int *crosslist = c->crosslist;
    int i;
    int rval;

    if (c->color == color) return 0;
    else if (c->color != -1) return 1;
    
    c->color = color;

    for (i=0; i<crosscount; i++) {
        rval = binesting_dfs (atoms, crosslist[i], 1-color);
        if (rval) return rval;
    }
    return 0;
}

static int nest_atom_color (atom_info *atoms, int color)
{
    int i;
    int cliquecount = atoms->cliquecount;
    vclique *cliquelist = atoms->cliquelist;
    int atomcount = atoms->atomcount;
    vclique *family = &atoms->family[color];
    vclique *c;
    int colorcount;
    int *perm = (int *) NULL;
    int *keys = (int *) NULL;
    int rval;

    perm = CC_SAFE_MALLOC (cliquecount, int);
    keys = CC_SAFE_MALLOC (cliquecount, int);
    if (perm == (int *) NULL ||
        keys == (int *) NULL) {
        fprintf (stderr, "Out of memory in nest_atom_color\n");
        rval = 1; goto CLEANUP;
    }

    family->parent = (vclique *) NULL;
    family->sibling = (vclique *) NULL;
    family->child = (vclique *) NULL;
    family->flipped = 0;
    family->mark = 0;
    family->atomcount = 0;
    family->atomlist = (int *) NULL;
    family->crosscount = 0;
    family->crosslist = (int *) NULL;
    family->color = color;
    
    colorcount = 0;
    for (i=0; i<cliquecount; i++) {
        c = &cliquelist[i];
        if (c->color == color) {
            if (c->flipped) keys[i] = atomcount - c->atomcount;
            else            keys[i] = c->atomcount;
            perm[colorcount++] = i;
        }
    }

    CCutil_int_perm_quicksort (perm, keys, colorcount);

    clear_atom_marks (atoms);
    clear_clique_marks (atoms);
    
    for (i=0; i<colorcount; i++) {
        rval = nest_clique (atoms, perm[i], family);
        if (rval) goto CLEANUP;
    }

    rval = 0;

  CLEANUP:
    CC_IFFREE (perm, int);
    CC_IFFREE (keys, int);
    return rval;
}

static int nest_clique (atom_info *atoms, int cliqnum, vclique *family)
{
    vclique *cliquelist = atoms->cliquelist;
    vclique *cliq = &cliquelist[cliqnum];
    atom *atomlist = atoms->atomlist;
    int catomcount;
    int *catomlist = (int *) NULL;
    vclique *c;
    atom *a;
    int i;
    int rval;

    if (cliq->flipped) {
        rval = clique_complement (atoms, cliq, &catomcount, &catomlist);
        if (rval) {
            fprintf (stderr, "clique_complement failed\n");
            goto CLEANUP;
        }
    } else {
        catomcount = cliq->atomcount;
        catomlist = cliq->atomlist;
    }
    
    cliq->child = (vclique *) NULL;
    
    for (i=0; i<catomcount; i++) {
        a = &atomlist[catomlist[i]];
        if (a->mark != 0) {
            c = &cliquelist[a->mark-1];
            if (c->mark == 1) {
                CHILD_DEL(c);
                c->mark = 0;
                CHILD_ADD (cliq, c);
            }
        }
        a->mark = cliqnum+1;
    }
    cliq->mark = 1;
    CHILD_ADD (family,cliq);

    rval = 0;
    
  CLEANUP:
    if (cliq->flipped) {
        CC_IFFREE (catomlist, int);
    }
    return rval;
}

static int clique_complement (atom_info *atoms, vclique *cliq, int *catomcount,
        int **catomlist)
{
    int atomcount = atoms->atomcount;
    int *cmark = (int *) NULL;
    int cnt;
    int i;
    int rval;

    *catomcount = atomcount - cliq->atomcount;
    *catomlist = CC_SAFE_MALLOC (*catomcount, int);
    cmark = CC_SAFE_MALLOC (atomcount, int);

    if (*catomlist == (int *) NULL ||
        cmark == (int *) NULL) {
        fprintf (stderr, "Out of memory in clique_complement\n");
        rval = 1; goto CLEANUP;
    }
    for (i=0; i<atomcount; i++) cmark[i] = 0;
    for (i=0; i<cliq->atomcount; i++) cmark[cliq->atomlist[i]] = 1;
    cnt = 0;
    for (i=0; i<atomcount; i++) {
        if (cmark[i] == 0) {
            (*catomlist)[cnt++] = i;
        }
    }
    if (cnt != *catomcount) {
        fprintf (stderr, "Lost some atoms\n");
        rval = 1; goto CLEANUP;
    }
    rval = 0;

 CLEANUP:
    if (rval) CC_IFFREE (*catomlist, int);
    CC_IFFREE (cmark, int);
    return rval;
}

/* try to flip to a more "natural" nesting */
static void renest_family (vclique *family)
{
    int cnt;
    vclique *cflip;
    vclique *c;

    for (;;) {
        cnt = 0;
        for (c = family->child; c; c = c->sibling) {
            if (c->flipped) {
                cnt++;
                cflip = c;
            }
        }
        if (cnt == 1) {
            flip_clique (family, cflip);
        } else {
            return;
        }
    }
}

static void flip_marked_cliques (vclique *family)
{
    int cnt;
    vclique *cflip = (vclique *) NULL;
    vclique *c;

    for (;;) {
        cnt = 0;
        for (c = family->child; c; c = c->sibling) {
            if (c->mark != 0) {
                cnt++;
                cflip = c;
            }
        }
        if (cnt == 1) {
            cflip->mark = 0;
            flip_clique (family, cflip);
        } else {
            return;
        }
    }
}

static void flip_clique (vclique *family, vclique *cliq)
{
    vclique *c;

    CHILD_DEL(cliq);
    c = cliq->child;

    cliq->child = family->child;
    if (cliq->child) cliq->child->prevsib = &(cliq->child);

    cliq->sibling = c;
    if (c) c->prevsib = &(cliq->sibling);

    family->child = cliq;
    cliq->prevsib = &(family->child);

    for (c = family->child; c; c = c->sibling) {
        c->parent = family;
    }

    for (c = cliq->child; c; c = c->sibling) {
        c->parent = cliq;
    }

    cliq->flipped = ! cliq->flipped;
}

static int verify_star (atom_info *atoms, int handnum,
        CCverify_cutclass *class)
{
    vclique *handles = &atoms->family[handnum];
    vclique *teeth = &atoms->family[1-handnum];
    int tcnt;
    int rval;
    int hcnt;
    int tsum;
    int tmult;
    vclique *c;
    vclique *t;
    int *alist = (int *) NULL;
    int acnt;
    int i;
    int nflipped;
    int cliquecount = atoms->cliquecount;

    clear_clique_marks (atoms);
    
    rval = flip_disjoint (teeth);
    if (rval) goto CLEANUP;

    rval = flip_concentric (handles);
    if (rval) goto CLEANUP;
    
    tcnt = child_count (teeth);
    if (tcnt < 3 || (tcnt % 2) != 1) {
        rval = 1; goto CLEANUP;
    }
    if (atoms->sense != 'G') {
        rval = 1; goto CLEANUP;
    }

    nflipped = 0;
    hcnt = family_count (handles, &nflipped) - 1;
    tsum = family_count (teeth, (int *) NULL) - 1;

    if (atoms->rhs != (tcnt+1) * hcnt + 2*tsum) {
        rval = 1; goto CLEANUP;
    }

    alist = CC_SAFE_MALLOC (atoms->atomcount, int);
    if (alist == (int *) NULL) {
        fprintf (stderr, "Out of memory in verify_star\n");
        rval = 1; goto CLEANUP;
    }
    
    clear_clique_marks (atoms);
    check_mark_concentric (handles->child, 1);

    for (c = teeth->child; c; c = c->sibling) {
        t = concentric_innermost (c);
        rval = collect_atom_mark_counts (atoms, t, alist, &acnt, nflipped);
        if (rval) goto CLEANUP;
        CCutil_int_array_quicksort (alist, acnt);
        if (alist[0] != 0 || alist[acnt-1] != hcnt) {
            rval = 1; goto CLEANUP;
        }
        tmult = family_count (c, (int *) NULL);
        for (i=1; i<acnt; i++) {
            if (alist[i] - alist[i-1] > tmult) {
                rval = 1; goto CLEANUP;
            }
        }
    }

    if (class != (CCverify_cutclass *) NULL) {
        rval = build_cutclass (class, cliquecount, tcnt+1);
        if (rval) {
            fprintf (stderr, "build_cutclass failed\n");
            goto CLEANUP;
        }
#if 0
        class->type = CC_TYPE_STAR;
        class->nhandles = 1;
        class->nfamilies = tcnt + 1;
        cnum = 0;
        fnum = 0;
        class->family_start[fnum] = cnum;
        fnum++;
        class->cliques[cnum] = handle;
        class->inverted[cnum] = 0;
        cnum++;
        for (i=0; i<cliquecount; i++) {
            if (i != handle) {
                class->family_start[fnum] = cnum;
                fnum++;
                class->cliques[cnum] = i;
                class->inverted[cnum] = cliquelist[i].flipped;
                cnum++;
            }
        }
#endif
    }

    rval = 0;

 CLEANUP:
    CC_IFFREE (alist, int);
    return rval;
}

/* verify_bipartition rounds the required coefficient on degenerate teeth
 * touching at least three handles up to the next integer.
 */
static int verify_bipartition (atom_info *atoms, int handnum)
{
    int teethnum = 1 - handnum;
    vclique *handles = &atoms->family[handnum];
    vclique *teeth = &atoms->family[teethnum];
    vclique *cliquelist = atoms->cliquelist;
    int rval;
    int tsum;
    int hsum;
    int tmult;
    vclique *c;
    vclique *t;
    int *alist = (int *) NULL;
    int i;
    int nflipped;

    clear_clique_marks (atoms);
    
    rval = flip_disjoint (teeth);
    if (rval) goto CLEANUP;

    rval = flip_disjoint_one (handles);
    if (rval) goto CLEANUP;

    if (atoms->sense != 'G') {
        rval = 1; goto CLEANUP;
    }

    clear_clique_marks (atoms);

    for (c = teeth->child; c; c = c->sibling) {
        t = concentric_innermost (c);
        for (i=0; i<t->crosscount; i++) {
            cliquelist[t->crosslist[i]].mark++;
        }
    }

    hsum = 0;
    for (c = handles->child; c; c = c->sibling) {
        tmult = family_count (c, (int *) NULL);
        if (tmult != 1 || (c->mark % 2) != 1) {
            rval = 1; goto CLEANUP;
        }
        hsum += c->mark + 1;
    }
    tsum = family_count (teeth, (int *) NULL) - 1;

    if (atoms->rhs != 2 * tsum + hsum) {
        rval = 1; goto CLEANUP;
    }

    clear_clique_marks (atoms);

    nflipped = 0;
    for (c = handles->child; c; c = c->sibling) {
        c->mark = 1;
        if (c->flipped) nflipped++;
    }

    for (c = teeth->child; c; c = c->sibling) {
        t = concentric_innermost (c);
        if (!tooth_has_cavity (atoms, t, nflipped)) {
            tmult = family_count (c, (int *) NULL);
            if (tmult < 2) {
                rval = 1; goto CLEANUP;
            }
        }
    }
    rval = 0;

 CLEANUP:
    CC_IFFREE (alist, int);
    return rval;
}

static int flip_disjoint (vclique *family)
{
    int ccnt;
    int rval;
    vclique *c;

    if (family == (vclique *) NULL || family->child == (vclique *) NULL) {
        return 1;
    }

    ccnt = child_count (family);
    if (ccnt == 1) {
        rval = check_mark_concentric (family->child, 1);
        if (rval) {
            flip_marked_cliques (family);
        }
    } else if (ccnt == 2) {
        rval = check_mark_concentric (family->child, 1);
        if (rval == 0) {
            check_mark_concentric (family->child, 0);
            rval = check_mark_concentric (family->child->sibling, 1);
        }
        if (rval) {
            flip_marked_cliques (family);
        }
    }

    for (c = family->child; c; c = c->sibling) {
        rval = check_mark_concentric (c, 0);
        if (rval) return 1;
    }

    return 0;
}

static int flip_disjoint_one (vclique *family)
{
    vclique *c;

    if (family == (vclique *) NULL || family->child == (vclique *) NULL) {
        return 1;
    }

    c = family->child;
    if (c->sibling == (vclique *) NULL) {
        if (c->child == (vclique *) NULL) {
            return 0;
        } else {
            flip_clique (family, c);
        }
    }

    for (c = family->child; c; c = c->sibling) {
        if (c->child != (vclique *) NULL) return 1;
    }

    return 0;
}

static int flip_concentric (vclique *family)
{
    int rval;

    if (family == (vclique *) NULL || family->child == (vclique *) NULL) {
        return 1;
    }

    if (child_count (family) == 2) {
        rval = check_mark_concentric (family->child, 1);
        if (rval == 0) {
            flip_marked_cliques (family);
        }
    }

    rval = check_mark_concentric (family, 0);
    return rval;
}

static int child_count (vclique *root)
{
    int n = 0;
    vclique *c;

    if (root == (vclique *) NULL) return -1;
    
    for (c = root->child; c; c = c->sibling) {
        n++;
    }
    return n;
}

/* family_count includes the root of the family.  Thus, when counting an
   entire family, subtract 1 for the dummy root */
static int family_count (vclique *root, int *nflipped)
{
    int n;
    vclique *c;

    if (root == (vclique *) NULL) return -1;

    n = 1;
    if (root->flipped && nflipped != (int *) NULL) (*nflipped) += 1;
    for (c = root->child; c; c = c->sibling) {
        n += family_count (c, nflipped);
    }
    return n;
}

static int check_mark_concentric (vclique *root, int mark)
{
    vclique *c;

    if (root == (vclique *) NULL) {
        return 0;
    }

    root->mark = mark;
    
    for (c = root->child; c; c = c->child) {
        if (c->sibling != (vclique *) NULL) {
            return 1;
        }
        c->mark = mark;
    }
    return 0;
}

static vclique *concentric_innermost (vclique *root)
{
    if (root == (vclique *) NULL) {
        return (vclique *) NULL;
    }
    while (root->child != (vclique *) NULL) {
        root = root->child;
    }
    return root;
}

static int collect_atom_mark_counts (atom_info *atoms, vclique *cliq,
        int *alist, int *p_acnt, int nflipped)
{
    int acnt = 0;
    atom *atomlist = atoms->atomlist;
    vclique *cliquelist = atoms->cliquelist;
    int catomcount;
    int *catomlist = (int *) NULL;
    int acliquecount;
    int *acliquelist;
    int i;
    int j;
    int a;
    int cnt;
    int fcnt;
    int rval;

    if (cliq->flipped) {
        rval = clique_complement (atoms, cliq, &catomcount, &catomlist);
        if (rval) {
            fprintf (stderr, "clique_complement failed\n");
            goto CLEANUP;
        }
    } else {
        catomcount = cliq->atomcount;
        catomlist = cliq->atomlist;
    }
        
    for (i=0; i<catomcount; i++) {
        a = catomlist[i];
        acliquecount = atomlist[a].cliquecount;
        acliquelist = atomlist[a].cliquelist;
        cnt = 0;
        fcnt = 0;
        for (j=0; j<acliquecount; j++) {
            if (cliquelist[acliquelist[j]].mark) {
                if (cliquelist[acliquelist[j]].flipped) fcnt++;
                else cnt++;
            }
        }
        alist[acnt++] = cnt + (nflipped - fcnt);
    }
    *p_acnt = acnt;

    rval = 0;

 CLEANUP:
    if (cliq->flipped) {
        CC_IFFREE (catomlist, int);
    }
    return rval;
}

static int tooth_has_cavity (atom_info *atoms, vclique *cliq, int nflipped)
{
    atom *atomlist = atoms->atomlist;
    vclique *cliquelist = atoms->cliquelist;
    int catomcount;
    int *catomlist = (int *) NULL;
    int acliquecount;
    int *acliquelist;
    int i;
    int j;
    int a;
    int cnt;
    int fcnt;
    int rval;

    if (cliq->flipped) {
        rval = clique_complement (atoms, cliq, &catomcount, &catomlist);
        if (rval) {
            fprintf (stderr, "clique_complement failed\n");
            goto CLEANUP;
        }
    } else {
        catomcount = cliq->atomcount;
        catomlist = cliq->atomlist;
    }
        
    for (i=0; i<catomcount; i++) {
        a = catomlist[i];
        acliquecount = atomlist[a].cliquecount;
        acliquelist = atomlist[a].cliquelist;
        cnt = 0;
        fcnt = 0;
        for (j=0; j<acliquecount; j++) {
            if (cliquelist[acliquelist[j]].mark) {
                if (cliquelist[acliquelist[j]].flipped) fcnt++;
                else cnt++;
            }
        }
        if (cnt + (nflipped - fcnt) == 0) {
            rval = 1; /* cavity found */
            goto CLEANUP;
        }
    }

    rval = 0;

 CLEANUP:
    if (cliq->flipped) {
        CC_IFFREE (catomlist, int);
    }
    return rval;
}

static int verify_other (atom_info *atoms)
{
    int rval;
    graph *g = (graph *) NULL;

    rval = build_complete_graph (atoms->atomcount, &g);
    if (rval) {
        fprintf (stderr, "build_complete_graph failed\n");
        goto CLEANUP;
    }

    compute_lhs (atoms, g);
    rval = verify_rhs (atoms, g);
    if (rval) {
        goto CLEANUP;
    }

    rval = 0;
  CLEANUP:
    if (g) free_graph (&g);
    return rval;
}

static int build_complete_graph (int nodecount, graph **p_g)
{
    graph *g = (graph *) NULL;
    int edgecount = (nodecount * (nodecount - 1)) / 2;
    int edgecount2;
    int *elist = (int *) NULL;
    int *elen = (int *) NULL;
    int i;
    int j;

    g = CC_SAFE_MALLOC (1, graph);
    if (!g) {
        fprintf (stderr, "Out of memory in build_graph\n");
        goto FAILURE;
    }
    g->nodecount = nodecount;
    g->edgecount = edgecount;
    g->elist = (int *) NULL;
    g->elen = (int *) NULL;

    g->elist = (int *) CC_SAFE_MALLOC (edgecount*2, int);
    if (!g->elist) {
        fprintf (stderr, "Out of memory in build_graph\n");
        goto FAILURE;
    }

    g->elen = (int *) CC_SAFE_MALLOC (edgecount, int);
    if (!g->elen) {
        fprintf (stderr, "Out of memory in build_graph\n");
        goto FAILURE;
    }

    edgecount2 = 0;
    elist = g->elist;
    elen = g->elen;
    for (i=0; i<nodecount; i++) {
        for (j=i+1; j<nodecount; j++) {
            elist[2*edgecount2] = i;
            elist[2*edgecount2+1] = j;
            elen[edgecount2] = 0;
            edgecount2++;
        }
    }
    if (edgecount2 != edgecount) {
        fprintf (stderr, "ERROR - found %d edges, expected %d\n",
                 edgecount2, edgecount);
        goto FAILURE;
    }

    *p_g = g;
    return 0;

  FAILURE:
    if (g) {
        CC_IFFREE (g->elen, int);
        CC_IFFREE (g->elist, int);
        CC_FREE (g, graph);
    }
    *p_g = (graph *) NULL;
    return -1;
}

static void free_graph (graph **p_g)
{
    if (*p_g) {
        CC_IFFREE ((*p_g)->elist, int);
        CC_IFFREE ((*p_g)->elen, int);
        CC_FREE (*p_g, graph)
    }
}

static void compute_lhs (atom_info *atoms, graph *g)
{
    int edgecount = g->edgecount;
    int *elen = g->elen;
    int *elist = g->elist;
    int cliquecount = atoms->cliquecount;
    vclique *cliquelist = atoms->cliquelist;
    atom *atomlist = atoms->atomlist;
    vclique *clique;
    int i;
    int j;

    for (i=0; i<edgecount; i++) {
        elen[i] = 0;
    }

    clear_atom_marks (atoms);
    for (i=0; i<cliquecount; i++) {
        clique = &cliquelist[i];
        for (j=0; j<clique->atomcount; j++) {
            atomlist[clique->atomlist[j]].mark = 1;
        }
        for (j=0; j<edgecount; j++) {
            if (atomlist[elist[2*j]].mark != atomlist[elist[2*j+1]].mark) {
                elen[j]++;
            }
        }
        for (j=0; j<clique->atomcount; j++) {
            atomlist[clique->atomlist[j]].mark = 0;
        }
    }
#ifdef DUMPSMALL
    if (atomcount < DUMPSMALL) {
        int k;
        printf ("\n%d\n", atomcount);
        for (i=0; i<atomcount; i++) {
            for (j=0; j<i; j++) {
                k = atomcount * (atomcount-1) / 2 -
                    (atomcount - j) * (atomcount - j - 1) / 2 + i - j - 1;
                if (elist[2*k] != j || elist[2*k+1] != i) {
                    fprintf (stderr, "Indexing problem, i %d j %d k %d elist[2k] %d elist[2k+1] %d\n",
                             i, j, k, elist[2*k], elist[2*k+1]);
                }
                printf ("%d ", elen[k]);
            }
            printf ("0");
            for (j=i+1; j<atomcount; j++) {
                k = atomcount * (atomcount-1) / 2 -
                    (atomcount - i) * (atomcount - i - 1) / 2 + j - i - 1;
                if (elist[2*k] != i || elist[2*k+1] != j) {
                    fprintf (stderr, "Indexing problem, i %d j %d k %d elist[2k] %d elist[2k+1] %d\n",
                             i, j, k, elist[2*k], elist[2*k+1]);
                }
                printf (" %d", elen[k]);
            }
            printf ("\n");
        }
        printf ("\n");
        fflush (stdout);
    }
#endif /* DUMPSMALL */
}

static int verify_rhs (atom_info *atoms, graph *g)
{
    int rval;
    double upbound;
    double optval;
    int foundtour;

    if (atoms->sense != 'G') {
        fprintf (stderr, "Constraint is not a >=\n");
        rval = 1;
        goto CLEANUP;
    }

    upbound = (double) atoms->rhs;

#ifdef DEBUG
    rval = CCheldkarp_small_elist (g->nodecount, g->edgecount, g->elist,
                                   g->elen, &upbound, &optval, &foundtour,
                                   1, (int *) NULL, HELDKARP_LIMIT, 1);
#else
    rval = CCheldkarp_small_elist (g->nodecount, g->edgecount, g->elist,
                                   g->elen, &upbound, &optval, &foundtour,
                                   1, (int *) NULL, HELDKARP_LIMIT, 2);
#endif
    if (rval == HELDKARP_SEARCHLIMITEXCEEDED) {
        fprintf (stderr, "CCheldkarp_small_elist search limit exceeded\n");
#ifdef DUMPFAIL
        {
            int i;
            printf ("Failed tsp\n");
            printf ("%d %d\n", g->nodecount, g->edgecount);
            for (i=0; i<g->edgecount; i++) {
                printf ("%d %d %d\n", g->elist[2*i], g->elist[2*i+1],
                        g->elen[i]);
            }
            printf ("%d\n", atoms->rhs);
            printf ("\n");
        }
#endif /* DUMPFAIL */
        rval = 1;
    } else if (rval) {
        fprintf (stderr, "CCheldkarp_small_elist failed\n");
        rval = 1;
    } else if (foundtour) {
        fprintf (stderr, "CCheldkarp_small_elist found better tour\n");
        rval = 1;
    } else {
        rval = 0;
    }

#if 0
    rval = CCtiny_bnc_tsp (g->nodecount, &dat, &upbound, (double *) NULL, 500);
    if (rval == CC_TINYTSP_INFEASIBLE) {
        rval = 0; goto CLEANUP;
    } else if (rval == CC_TINYTSP_ERROR) {
        fprintf (stderr, "Tinytsp failed\n");
        rval = 1; goto CLEANUP;
    } else if (rval == CC_TINYTSP_SEARCHLIMITEXCEEDED) {
        fprintf (stderr, "Tinytsp search limit exceeded\n");
        rval = 1; goto CLEANUP;
    } else if (rval) {
        fprintf (stderr, "Tinytsp invalid return code %d\n", rval);
        rval = 1; goto CLEANUP;
    } else {
        fprintf (stderr, "Tinytsp found better tour\n");
        rval = 1; goto CLEANUP;
    }
#endif

CLEANUP:
    return rval;
}
      
static void cleanup_atom_info (CCtsp_lpcut_in *cut, atom_info *atoms)
{
    node *nodelist = atoms->nodelist;
    CCtsp_lpclique *cliques = cut->cliques;
    int cliquecount = cut->cliquecount;
    int i;
    int k;
    int j;

    for (i=0; i<cliquecount; i++) {
        CC_FOREACH_NODE_IN_CLIQUE (j, cliques[i], k) {
            nodelist[j].atom = -1;
        }
    }
}

static void clear_atom_marks (atom_info *atoms)
{
    atom *atomlist = atoms->atomlist;
    int atomcount = atoms->atomcount;
    int i;

    for (i=0; i<atomcount; i++) {
        atomlist[i].mark = 0;
    }
}

static void clear_clique_marks (atom_info *atoms)
{
    vclique *cliquelist = atoms->cliquelist;
    int cliquecount = atoms->cliquecount;
    int i;

    for (i=0; i<cliquecount; i++) {
        cliquelist[i].mark = 0;
    }
}

#if 0
void CCverify_initcutclass (CCverify_cutclass *class)
{
    class->type         = 0;
    class->nhandles     = 0;
    class->nfamilies    = 0;
    class->cliques      = (int *) NULL;
    class->inverted     = (int *) NULL;
    class->family_start = (int *) NULL;
}

void CCverify_freecutclass (CCverify_cutclass *class)
{
    CC_IFFREE (class->cliques, int);
    CC_IFFREE (class->inverted, int);
    CC_IFFREE (class->family_start, int);
    CCverify_initcutclass (class);
}

int CCverify_classify (CCtsp_lpcut_in *cut, CCverify_cutclass *class)
{
    atom_info atoms;
    int rval;
    int handle;

#ifdef DEBUG
    double sz = CCutil_zeit();
    printf ("Classifying cut, %d cliques %d atoms",cut->cliquecount,
            cut->skel.atomcount);
    fflush (stdout);
#endif /* DEBUG */

    class->type = -1;
    
    if (cut->skel.atomcount == 0 || cut->skel.atoms == (int *) NULL) {
        fprintf (stderr, "Cut to classify has no skeleton\n");
        rval = 1; goto CLEANUP;
    }
    
    rval = build_atom_info (cut, &atoms);
    if (rval) {
        fprintf (stderr, "build_atom_info failed\n");
        goto CLEANUP;
    }

    rval = verify_atom_info (cut, &atoms);
    if (rval) {
        fprintf (stderr, "atom_info failed verification\n");
        goto CLEANUP;
    }

    rval = verify_subtour (&atoms, class);
    if (rval == 0) {
        /* It's a valid subtour - we're done */
#ifdef DEBUG
        printf (" (subtour)");
        fflush (stdout);
#endif
        goto CLEANUP;
    }

    rval = verify_chvatal_comb (&atoms, class);
    if (rval == 0) {
        /* It's a valid comb - we're done */
#ifdef DEBUG
        printf (" (comb 1,%d)", atoms.cliquecount-1);
        fflush (stdout);
#endif
        goto CLEANUP;
    }

    rval = build_binesting (&atoms);
    if (rval == 0) {
        rval = verify_star (&atoms, 0, class);
        if (rval == 0) {
            /* It's a valid star, with handles marked 0 - we're done */
#ifdef DEBUG
            int hcnt = family_count (&atoms.family[0], (int *) NULL) - 1;
            printf (" (star %d,%d)", hcnt, atoms.cliquecount - hcnt);
            fflush (stdout);
#endif
            if (type != (int *) NULL) *type = CC_TYPE_STAR;
            goto CLEANUP;
        }

        rval = verify_star (&atoms, 1, class);
        if (rval == 0) {
            /* It's a valid star, with handles marked 1 - we're done */
#ifdef DEBUG
            int hcnt = family_count (&atoms.family[1], (int *) NULL) - 1;
            printf (" (star %d,%d)", hcnt, atoms.cliquecount - hcnt);
            fflush (stdout);
#endif
            if (type != (int *) NULL) *type = CC_TYPE_STAR;
            goto CLEANUP;
        }

        rval = verify_bipartition (&atoms, 0);
        if (rval == 0) {
            /* It's a valid bipartition, with handles marked 0 - we're done */
#ifdef DEBUG
            int hcnt = family_count (&atoms.family[0], (int *) NULL) - 1;
            printf (" (bipartition %d,%d)", hcnt, atoms.cliquecount - hcnt);
            fflush (stdout);
#endif
            if (type != (int *) NULL) *type = CC_TYPE_BIPARTITION;
            goto CLEANUP;
        }

        rval = verify_bipartition (&atoms, 1);
        if (rval == 0) {
            /* It's a valid bipartition, with handles marked 1 - we're done */
#ifdef DEBUG
            int hcnt = family_count (&atoms.family[1], (int *) NULL) - 1;
            printf (" (bipartition %d,%d)", hcnt, atoms.cliquecount - hcnt);
            fflush (stdout);
#endif
            if (type != (int *) NULL) *type = CC_TYPE_BIPARTITION;
            goto CLEANUP;
        }

#if 0
#ifdef DEBUG
        printf ("\nCurious, binested, but not star or bipartition:\n");
        CCtsp_print_lpcut_in (cut);
        fflush (stdout);
#endif
#endif
    }

    

    rval = verify_other (&atoms);
    if (rval == 0) {
#ifdef DEBUG
        printf (" (other)");
        fflush (stdout);
#endif
        if (type != (int *) NULL) *type = CC_TYPE_OTHER;
        goto CLEANUP;
    }

    fprintf (stderr, "Unable to verify cut\n");
    rval = -1;

  CLEANUP:
#ifdef DEBUG
    if (rval == 0) {
        printf (" in %.2f seconds\n", CCutil_zeit() - sz);
        fflush (stdout);
    } else {
        printf (" FAILED in %.2f seconds\n", CCutil_zeit() - sz);
        printf ("FAILED CUT:\n");
        CCtsp_print_lpcut_in (cut);
        fflush (stdout);
    }
#endif
    free_atom_info (&atoms);
    return rval;
    
}
#endif

static int build_cutclass (CCverify_cutclass *class, int ncliques,
        int nfamilies)
{
    class->cliques = CC_SAFE_MALLOC (ncliques, int);
    class->inverted = CC_SAFE_MALLOC (ncliques, int);
    class->family_start = CC_SAFE_MALLOC (nfamilies, int);

    if (class->cliques == (int *) NULL ||
        class->inverted == (int *) NULL ||
        class->family_start == (int *) NULL) {
        fprintf (stderr, "Out of memory in build_cutclass\n");
        CC_IFFREE (class->cliques, int);
        CC_IFFREE (class->inverted, int);
        CC_IFFREE (class->family_start, int);
        return 1;
    }
    return 0;
}

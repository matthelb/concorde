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
/*               A PROGRAM TO DISPLAY RESTART FILES                         */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: January 6, 1998                                                   */
/*                                                                          */
/*  SEE short decsribtion in usage ().                                      */
/*                                                                          */
/*  This file currently has a copy of the restart reading code from         */
/*  TSP/bcontrol.c.  In the long run, this code should be shared instead    */
/*  of copied.                                                              */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"

typedef struct tsp_bbnode {
    int id;
    int status;
    int workstatus;
    int numtentative;
    struct tsp_bbnode *prev;
    struct tsp_bbnode *next;
    struct tsp_bbnode *parent;
    struct tsp_bbnode *child0;
    struct tsp_bbnode *child1;
    struct tsp_tnode *tentative_nodes;
    struct tsp_tnode *tparent;
    double lowerbound;
    double cputime;
    int number;
} tsp_bbnode;

typedef struct tsp_tnode {
    tsp_bbnode *parent;
    tsp_bbnode *child0;
    tsp_bbnode *child1;
} tsp_tnode;

#define BB_NEEDS_CUTTING             (1)
#define BB_NEEDS_TENTATIVE_CUTTING   (2)
#define BB_NEEDS_BRANCHING           (3)
#define BB_NEEDS_TENTATIVE_BRANCHING (4)
#define BB_DONE                      (5)
#define BB_IDLE              (1)
#define BB_WORKING           (2)
#define BB_PRUNED            (3)

#define TASK_WAIT_SECONDS  (60)

#define DRAWING_WIDTH (6.5*72)
#define DRAWING_HEIGHT (9*72)
#define DRAWING_BOTTOM (72)
#define DRAWING_LEFT (36+72)
#define DRAWING_TIC (7)
#define DRAWING_TICGAP (5)
#define DRAWING_TICMARGIN (5)
#define DRAWING_TITLEGAP (20)

CC_PTRWORLD_ROUTINES (tsp_bbnode, tsp_bbnode_alloc, tsp_bbnode_bulkalloc,
        tsp_bbnode_free)
CC_PTRWORLD_LEAKS_ROUTINE (tsp_bbnode, tsp_bbnode_check_leaks, id, int)

static char *resfname = (char *) NULL;
static int arg_maxdepth = -1;
static int leafsummary = 0;
static int nodelist = 0;
static int dig_before = -1;
static int dig_after = 2;
static int use_graphics = 0;

static int
    draw_tree (char *probname, double restart_upbound,
        CC_UNUSED int bbcount, CC_UNUSED double branchzeit, int maxdepth,
        tsp_bbnode *rootbbnode),
    report_leaves (char *probname, double restart_upbound,
        int bbcount, double branchzeit, double mod, char *format,
        tsp_bbnode *rootbbnode, int shownodes),
    report_tree (char *probname, double restart_upbound,
        int bbcount, double branchzeit, int maxdepth, double mod, char *format,
        tsp_bbnode *rootbbnode),
    read_restart (char *restart_name, char **p_probname,
        tsp_bbnode **p_rootbbnode, double *p_upbound, int *p_ncount,
        int *p_bbcount, double *p_branchzeit, CCptrworld *bbnode_world),
    read_bbtree (FILE *f, tsp_bbnode **p_b, CCptrworld *bbnode_world),
    read_tentative_nodes (FILE *f, int count, tsp_tnode **list,
        tsp_bbnode *parent, CCptrworld *bbnode_world),
    prob_name (char *buf, size_t buflen, const char *f, int n),
    parseargs (int ac, char **av);

static void
    number_tree (tsp_bbnode *bbnode, int *cnt, int depth, int maxdepth,
        double upbound),
    draw_prelude (int drawcnt, double lobound, double upbound, char *probname),
    draw_postlude (void),
    draw_node (tsp_bbnode *bbnode, int depth),
    draw_edge (tsp_bbnode *bfrom, tsp_bbnode *bto),
    draw_subtree_edges (tsp_bbnode *bbnode, int depth, int maxdepth,
        double upbound),
    draw_subtree_nodes (tsp_bbnode *bbnode, int depth, int maxdepth,
        double upbound),
    print_node (tsp_bbnode *bbnode, double mod, char *format),
    show_node (tsp_bbnode *bbnode, double mod, char *format),
    output_tree (tsp_bbnode *bbnode, char *buf, int buflen,
        int depth, int maxdepth, double mod, char *format),
    collect_leaves (tsp_bbnode *bbnode, double *lower, int *leafcount,
        double *leafvals, tsp_bbnode **leafbbs),
    init_bbnode (tsp_bbnode *bbnode),
    free_tree (tsp_bbnode **bbnode, CCptrworld *bbnode_world),
    usage (char *fname);

static double
    stepsize (double lo, double hi, int maxsteps);


int main (int ac, char **av)
{
    int rval = 0;
    tsp_bbnode *rootbbnode  = (tsp_bbnode *) NULL;
    char *probname = (char *) NULL;
    double restart_upbound = 0.0;
    int restart_ncount = 0;
    int bbcount = 0;
    double branchzeit = 0.0;
    char format[1024];
    double mod = -1.0;
    CCptrworld bbnode_world;

    CCptrworld_init (&bbnode_world);
    
    rval = parseargs (ac, av);
    if (rval) return 1;

    if (dig_before > 0) {
        int i;
        mod = 1.0;
        for (i=0; i<dig_before; i++) mod = mod * 10.0;
        sprintf (format, "%%0%d.%df", dig_before + dig_after + 1,
                 dig_after);
    } else {
        mod = -1.0;
        sprintf (format, "%%.%df", dig_after);
    }
        
    rval = read_restart (resfname, &probname, &rootbbnode, &restart_upbound,
                        &restart_ncount, &bbcount, &branchzeit, &bbnode_world);
    if (rval) {
        fprintf (stderr, "read_restart failed\n");
        goto CLEANUP;
    }

    if (use_graphics) {
        rval = draw_tree (probname, restart_upbound, bbcount, branchzeit,
                          arg_maxdepth, rootbbnode);
        if (rval) {
            fprintf (stderr, "draw_tree failed\n");
            goto CLEANUP;
        }
    } else if (leafsummary || nodelist) {
        rval = report_leaves (probname, restart_upbound, bbcount, branchzeit,
                              mod, format, rootbbnode, nodelist);
        if (rval) {
            fprintf (stderr, "report_leaves failed\n");
            goto CLEANUP;
        }
    } else {
        rval = report_tree (probname, restart_upbound, bbcount, branchzeit,
                            arg_maxdepth, mod, format, rootbbnode);
        if (rval) {
            fprintf (stderr, "report_tree failed\n");
            goto CLEANUP;
        }
    }

    rval = 0;

  CLEANUP:
    CC_IFFREE (probname, char);
    free_tree (&rootbbnode, &bbnode_world);
    CCptrworld_delete (&bbnode_world);
    return rval;
}

static void number_tree (tsp_bbnode *bbnode, int *cnt, int depth, int maxdepth,
        double upbound)
{
    if (bbnode->child0 != (tsp_bbnode *) NULL &&
        (maxdepth <= 0 || depth < maxdepth) &&
        bbnode->child0->lowerbound <= upbound - 1.0) {
        number_tree (bbnode->child0, cnt, depth + 1, maxdepth, upbound);
    }
    bbnode->number = *cnt;
    (*cnt)++;
    if (bbnode->child1 != (tsp_bbnode *) NULL &&
        (maxdepth <= 0 || depth < maxdepth) &&
        bbnode->child1->lowerbound <= upbound - 1.0) {
        number_tree (bbnode->child1, cnt, depth + 1, maxdepth, upbound);
    }
}
    
static int draw_tree (char *probname, double restart_upbound,
        CC_UNUSED int bbcount, CC_UNUSED double branchzeit, int maxdepth,
        tsp_bbnode *rootbbnode)
{
    int drawcnt = 0;

    number_tree (rootbbnode, &drawcnt, 0, maxdepth, restart_upbound);
    draw_prelude (drawcnt, rootbbnode->lowerbound, restart_upbound, probname);
    draw_subtree_edges (rootbbnode, 0, maxdepth, restart_upbound);
    draw_subtree_nodes (rootbbnode, 0, maxdepth, restart_upbound);
    draw_postlude ();
    return 0;
}

static double stepsize (double lo, double hi, int maxsteps)
{
    double step = 1.0;

    if (hi <= lo) hi = lo + 1.0;
    
    while ((hi - lo) / step <= maxsteps) {
        step = step / 10.0;
    }
    while ((hi - lo) / step > maxsteps) {
        if ((hi - lo) / (2.0 * step) <= maxsteps) return 2.0*step;
        if ((hi - lo) / (5.0 * step) <= maxsteps) return 5.0*step;
        step = step * 10.0;
    }
    return step;
}

static void draw_prelude (int drawcnt, double lobound, double upbound,
        char *probname)
{
    double labelstep = stepsize (lobound, upbound, 60);
    double labello;
    double m;
    
    printf ("%%!PS-Adobe-2.0\n");
    printf ("/yloc {%.2f exch sub %.2f mul %.2f div %.2f add} def\n", upbound,
            (double) DRAWING_HEIGHT, upbound - lobound,
            (double) DRAWING_BOTTOM);
    printf ("/xloc {%.2f mul %d div %.2f add} def\n", (double) DRAWING_WIDTH,
            drawcnt, (double) DRAWING_LEFT);
    printf ("/loc {yloc exch xloc exch} def\n");
    printf ("/n {pop loc 1 0 360 arc fill} def\n");
    printf ("/e {newpath 0.75 0.75 0.75 setrgbcolor loc moveto loc lineto stroke} def\n");
    printf ("/dn {0 1 0 setrgbcolor n} def\n");
    printf ("/cn {1 0 0 setrgbcolor n} def\n");
    printf ("/bn {1 0 1 setrgbcolor n} def\n");
    printf ("/un {1 1 0 setrgbcolor n} def\n");
    printf ("/tbn {un} def\n");
    printf ("/tcn {un} def\n");
    printf ("/dw {gsave newpath 0 0 moveto true charpath pathbbox\n");
    printf ("     exch 4 3 roll sub 3 1 roll sub grestore} def\n");
    printf ("/lbl {gsave loc translate newpath 0 0 moveto %.2f 0 lineto stroke\n",
            (double) -DRAWING_TIC);
    printf ("      dup dw %.2f 0 moveto 2 div exch neg exch rmoveto show grestore} def\n",
            (double) -(DRAWING_TIC + DRAWING_TICGAP));
    printf ("gsave\n");
    printf ("/Helvetica findfont 24 scalefont setfont\n");
    printf ("0 0 0 setrgbcolor\n");
    printf ("newpath\n");
    printf ("%.2f 2 div %.2f add %.2f %.2f add moveto\n",
            (double) DRAWING_WIDTH, (double) DRAWING_LEFT,
            (double) DRAWING_HEIGHT, (double) DRAWING_BOTTOM);
    printf ("(%s Branching Tree) dup stringwidth pop 2 div neg %.2f rmoveto show\n",
            probname, (double) DRAWING_TITLEGAP);
    printf ("grestore\n");
    
    printf ("gsave\n");
    printf ("0 setlinewidth\n");
    printf ("0 0 0 setrgbcolor\n");
    printf ("/Helvetica findfont 10 scalefont setfont\n");
    printf ("%.2f 0 translate\n", (double) -DRAWING_TICMARGIN);

    m = fmod (lobound, labelstep);
    if (m == 0.0) labello = lobound;
    else          labello = lobound + labelstep - m;

    printf ("newpath 0 %.2f loc moveto 0 %.2f loc lineto stroke\n",
            lobound, upbound);
    
    for (m = labello; m <= upbound; m += labelstep) {
        printf ("(%.0f) 0 %.2f lbl\n", m, m);
    }

    printf ("grestore\n");

    printf ("gsave\n");
    printf ("0.5 setlinewidth\n");
    printf ("0 0 0 setrgbcolor\n");
    printf ("\n");
}

static void draw_postlude (void)
{
    printf ("grestore\n");
    printf ("showpage\n");
    printf ("%%EOF\n");
}
    
static void draw_node (tsp_bbnode *bbnode, int depth)
{
    printf ("%d %.2f %d ", bbnode->number, bbnode->lowerbound, depth);
    if (bbnode->status == BB_NEEDS_CUTTING) {
        printf ("cn\n");
    } else if (bbnode->status == BB_NEEDS_BRANCHING) {
        printf ("bn\n");
    } else if (bbnode->status == BB_NEEDS_TENTATIVE_CUTTING) {
        printf ("tcn\n");
    } else if (bbnode->status == BB_NEEDS_TENTATIVE_BRANCHING) {
        printf ("tbn\n");
    } else if (bbnode->status == BB_DONE) {
        printf ("dn\n");
    } else {
        printf ("un\n");
    }
}

static void draw_edge (tsp_bbnode *bfrom, tsp_bbnode *bto)
{
    printf ("%d %.2f %d %.2f e\n", bfrom->number, bfrom->lowerbound,
            bto->number, bto->lowerbound);
}

static void draw_subtree_edges (tsp_bbnode *bbnode, int depth, int maxdepth,
        double upbound)
{
    if (bbnode->child0 != (tsp_bbnode *) NULL &&
        (maxdepth <= 0 || depth < maxdepth) &&
        bbnode->child0->lowerbound <= upbound - 1.0) {
        draw_edge (bbnode, bbnode->child0);
        draw_subtree_edges (bbnode->child0, depth+1, maxdepth, upbound);
    }
    if (bbnode->child1 != (tsp_bbnode *) NULL &&
        (maxdepth <= 0 || depth < maxdepth) &&
        bbnode->child1->lowerbound <= upbound - 1.0) {
        draw_edge (bbnode, bbnode->child1);
        draw_subtree_edges (bbnode->child1, depth+1, maxdepth, upbound);
    }
}

static void draw_subtree_nodes (tsp_bbnode *bbnode, int depth, int maxdepth,
        double upbound)
{
    draw_node (bbnode, depth);
    if (bbnode->child0 != (tsp_bbnode *) NULL &&
        (maxdepth <= 0 || depth < maxdepth) &&
        bbnode->child0->lowerbound <= upbound - 1.0) {
        draw_subtree_nodes (bbnode->child0, depth+1, maxdepth, upbound);
    }
    if (bbnode->child1 != (tsp_bbnode *) NULL &&
        (maxdepth <= 0 || depth < maxdepth) &&
        bbnode->child1->lowerbound <= upbound - 1.0) {
        draw_subtree_nodes (bbnode->child1, depth+1, maxdepth, upbound);
    }
}

static int report_leaves (char *probname, double restart_upbound,
        int bbcount, double branchzeit, double mod, char *format,
        tsp_bbnode *rootbbnode, int shownodes)
{
    double lower = 0.0;
    int leafcount = 0;
    double *leafvals = (double *) NULL;
    tsp_bbnode **leafbbs = (tsp_bbnode **) NULL;
    int *leafperm = (int *) NULL;
    int i;
    int rval = 0;
    char buf[80];
    size_t outcnt;
    double v;

    leafvals = CC_SAFE_MALLOC (bbcount, double);
    if (leafvals == (double *) NULL) {
        fprintf (stderr, "Out of memory in report_leaves\n");
        rval = 1; goto CLEANUP;
    }

    if (shownodes) {
        leafbbs = CC_SAFE_MALLOC (bbcount, tsp_bbnode *);
        if (leafbbs == (tsp_bbnode **) NULL) {
            fprintf (stderr, "Out of memory in report_leaves\n");
            rval = 1; goto CLEANUP;
        }
    }

    lower = restart_upbound;
    collect_leaves (rootbbnode, &lower, &leafcount, leafvals, leafbbs);

    printf ("%s: >= %.2f <= %.2f bb %d active %d time %.2f\n", probname, lower,
            restart_upbound, bbcount, leafcount, branchzeit);

    if (leafcount == 0) {
        rval = 0;
        goto CLEANUP;
    }

    leafperm = CC_SAFE_MALLOC (leafcount, int);
    if (leafperm == (int *) NULL) {
        fprintf (stderr, "Out of memory in report_leaves\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<leafcount; i++) leafperm[i] = i;

    CCutil_double_perm_quicksort (leafperm, leafvals, leafcount);

    outcnt = 0;
    for (i=0; i<leafcount; i++) {
        if (shownodes) {
            show_node (leafbbs[leafperm[i]], mod, format);
        } else {
            v = leafvals[leafperm[i]];
            if (mod > 0.0) v = fmod(v, mod);
            sprintf (buf, format, v);
            outcnt += strlen(buf) + 1;
            if (outcnt >= 75) {
                printf ("\n");
                outcnt = strlen(buf) + 1;
            }
            printf ("%s ",buf);
        }
    }
    if (!shownodes) {
        printf ("\n");
    }
    rval = 0;
    
 CLEANUP:
    CC_IFFREE (leafvals, double);
    CC_IFFREE (leafbbs, tsp_bbnode *);
    CC_IFFREE (leafperm, int);
    return rval;
}

static void print_node (tsp_bbnode *bbnode, double mod, char *format)
{
    double v = bbnode->lowerbound;
    char buf[80];
    
    switch (bbnode->status) {
    case BB_NEEDS_CUTTING:
        printf ("C"); break;
    case BB_NEEDS_BRANCHING:
        printf ("B"); break;
    case BB_NEEDS_TENTATIVE_CUTTING:
        printf ("T"); break;
    case BB_NEEDS_TENTATIVE_BRANCHING:
        printf ("t"); break;
    default:
        printf ("?%d", bbnode->status); break;
    }
    printf (" %5d ", bbnode->id);
    if (mod > 0.0) v = fmod (v, mod);
    printf (format, v);

    (void) prob_name (buf, sizeof (buf), "", bbnode->id);
    printf (" %s\n", buf);
}

static void show_node (tsp_bbnode *bbnode, double mod, char *format)
{
    tsp_bbnode *b;
    int j;

    print_node (bbnode, mod, format);

    if (bbnode->status == BB_NEEDS_BRANCHING &&
        bbnode->numtentative > 0) {
        for (j=0; j<bbnode->numtentative; j++) {
            b = bbnode->tentative_nodes[j].child0;
            if (b && b->status != BB_DONE) {
                printf ("   ");
                print_node (b, mod, format);
            }
            b = bbnode->tentative_nodes[j].child1;
            if (b && b->status != BB_DONE) {
                printf ("   ");
                print_node (b, mod, format);
            }
        }
    }
}

static int report_tree (char *probname, double restart_upbound,
        int bbcount, double branchzeit, int maxdepth, double mod, char *format,
        tsp_bbnode *rootbbnode)
{
    double lower = 0.0;
    int leafcount = 0;
    char buf[16384];

    lower = restart_upbound;
    collect_leaves (rootbbnode, &lower, &leafcount, (double *) NULL, (tsp_bbnode **) NULL);

    printf ("%s: >= %.2f <= %.2f bb %d active %d time %.2f\n", probname, lower,
            restart_upbound, bbcount, leafcount, branchzeit);

    output_tree (rootbbnode, buf, 0, 0, maxdepth, mod, format);

    return 0;
}

static void output_tree (tsp_bbnode *bbnode, char *buf, int buflen,
        int depth, int maxdepth, double mod, char *format)
{
    char mybuf[1024];
    double v = bbnode->lowerbound;
    size_t outcnt;
    size_t i;

    if (mod > 0.0) v = fmod(v, mod);
    sprintf (mybuf, "%d ", bbnode->id);
    sprintf (mybuf + strlen(mybuf), format, v);
    printf ("%s", mybuf);
    if (bbnode->child0 == (tsp_bbnode *) NULL &&
        bbnode->child1 == (tsp_bbnode *) NULL) {
        if (bbnode->status == BB_DONE) printf ("X");
        printf ("\n");
    } else if (maxdepth > 0 && depth >= maxdepth) {
        if (bbnode->status == BB_DONE) printf ("+");
        printf ("\n");
    } else {
        printf (" ");
        outcnt = strlen(mybuf) + 1;
        if (bbnode->parent == (tsp_bbnode *) NULL ||
            bbnode->parent->child1 == bbnode) {
            buf[buflen] = ' ';
        } else {
            buf[buflen] = '|';
        }
        for (i=1; i<outcnt; i++) buf[buflen+i] = ' ';
        buflen += outcnt;
        output_tree (bbnode->child0, buf, buflen, depth+1, maxdepth,
                     mod, format);
        buf[buflen] = '\0';
        printf ("%s", buf);
        output_tree (bbnode->child1, buf, buflen, depth+1, maxdepth,
                     mod, format);
    }
}

static int read_restart (char *restart_name, char **p_probname,
        tsp_bbnode **p_rootbbnode, double *p_upbound, int *p_ncount,
        int *p_bbcount, double *p_branchzeit, CCptrworld *bbnode_world)
{
    FILE *f = (FILE *) NULL;
    char *probname = (char *) NULL;
    int rval = 0;

    f = fopen (restart_name, "r");
    if (f == (FILE*) NULL) {
        perror (restart_name);
        fprintf (stderr, "Unable to open %s for input in read_restart\n",
                 restart_name);
        rval = 1; goto CLEANUP;
    }

    probname = CC_SAFE_MALLOC (CCtsp_PROB_FILE_NAME_LEN, char);
    if (probname == (char *) NULL) {
        fprintf (stderr, "Out of memory in read_restart\n");
        rval = 1; goto CLEANUP;
    }

    CCutil_readstr (f, probname, CCtsp_PROB_FILE_NAME_LEN);
    
    rval = fscanf (f, "%d%lf%d%lf\n", p_ncount, p_upbound, p_bbcount,
            p_branchzeit);
    if (rval <= 0) {
        perror (restart_name);
        fprintf (stderr, "fscanf from %s failed\n", restart_name);
        rval = 1; goto CLEANUP;
    }
    rval = read_bbtree (f, p_rootbbnode, bbnode_world);
    if (rval) {
        fprintf (stderr, "read_bbtree failed\n"); goto CLEANUP;
    }
    
    rval = fclose (f);
    if (rval) {
        perror (restart_name);
        fprintf (stderr, "fclose %s failed\n", restart_name);
        rval = 1; goto CLEANUP;
    }
    f = (FILE *) NULL;

    *p_probname = probname;
    
    rval = 0;

  CLEANUP:
    if (f != (FILE *) NULL) {
        fclose (f);
    }
    if (rval) {
        CC_IFFREE (probname, char);
        free_tree (p_rootbbnode, bbnode_world);
    }
    return rval;
}

static int read_bbtree (FILE *f, tsp_bbnode **p_b, CCptrworld *bbnode_world)
{
    int rval = 0;
    int child0, child1;
    tsp_bbnode *b = (tsp_bbnode *) NULL;

    b = tsp_bbnode_alloc (bbnode_world);
    if (b == (tsp_bbnode *) NULL) {
        fprintf (stderr, "tsp_bbnode_alloc failed\n");
        rval = 1; goto CLEANUP;
    }
    init_bbnode (b);
    
    rval = fscanf (f, "%d%d%d%d%d%lf%lf\n", &b->status, &b->id,
                    &child0, &child1, &b->numtentative, &b->lowerbound,
                    &b->cputime);
    if (rval <= 0) {
        perror ("restart_file");
        fprintf (stderr, "fscanf failed reading restart file\n");
        rval = 1; goto CLEANUP;
    }
    b->workstatus = BB_IDLE;

    if (b->numtentative > 0) {
        rval = read_tentative_nodes (f, b->numtentative, &b->tentative_nodes,
                                     b, bbnode_world);
        if (rval) {
            fprintf (stderr, "read_tentative_nodes failed\n");
            goto CLEANUP;
        }
    }

    if (child0 != -1) {
        rval = read_bbtree (f, &(b->child0), bbnode_world);
        if (rval) goto CLEANUP;
        if (b->child0->id != child0) {
            fprintf (stderr, "syntax error in restart file\n");
            rval = 1; goto CLEANUP;
        }
        b->child0->parent = b;
    }
    if (child1 != -1) {
        rval = read_bbtree (f, &(b->child1), bbnode_world);
        if (rval) goto CLEANUP;
        if (b->child1->id != child1) {
            fprintf (stderr, "syntax error in restart file\n");
            rval = 1; goto CLEANUP;
        }
        b->child1->parent = b;
    }

    *p_b = b;
    rval = 0;

  CLEANUP:
    if (rval) {
        tsp_bbnode_free (bbnode_world, b);
    }
    return rval;
}

static int read_tentative_nodes (FILE *f, int count, tsp_tnode **list,
        tsp_bbnode *parent, CCptrworld *bbnode_world)
{
    int rval = 0;
    int i;
    int obtained = 0;
    tsp_tnode *s;
    tsp_bbnode *child0 = (tsp_bbnode *) NULL;
    tsp_bbnode *child1 = (tsp_bbnode *) NULL;

    *list = CC_SAFE_MALLOC (count, tsp_tnode);
    if (*list == (tsp_tnode *) NULL) {
        fprintf (stderr, "out of memory in read_tentative_nodes\n");
        rval = 1; goto CLEANUP;
    }

    for (obtained = 0; obtained < count; obtained++) {
        s = &((*list)[obtained]);
        child0 = tsp_bbnode_alloc (bbnode_world);
        if (child0 == (tsp_bbnode *) NULL) {
            fprintf (stderr, "tsp_bbnode_alloc failed\n");
            rval = 1; goto CLEANUP;
        }
        init_bbnode (child0);
        child1 = tsp_bbnode_alloc (bbnode_world);
        if (child1 == (tsp_bbnode *) NULL) {
            fprintf (stderr, "tsp_bbnode_alloc failed\n");
            tsp_bbnode_free (bbnode_world, child0);
            rval = 1; goto CLEANUP;
        }
        init_bbnode (child1);

        rval = fscanf (f, "%d%d%lf%lf %d%d%lf%lf\n",
                    &child0->status, &child0->id,
                    &child0->lowerbound, &child0->cputime,
                    &child1->status, &child1->id,
                    &child1->lowerbound, &child1->cputime);
        if (rval <= 0) {
            perror ("restart_file");
            fprintf (stderr, "fscanf failed reading tentative line\n");
            rval = 1; goto CLEANUP;
        }
        child0->tparent = s;
        child1->tparent = s;
        s->child0 = child0;
        s->child1 = child1;
        s->parent = parent;
    }
    rval = 0;

CLEANUP:

    if (rval) {
        for (i = 0; i < obtained; i++) {
            tsp_bbnode_free (bbnode_world, (*list)[i].child0);
            tsp_bbnode_free (bbnode_world, (*list)[i].child1);
        }
        CC_IFFREE (*list, tsp_tnode);
    }
    return rval;
}

static void collect_leaves (tsp_bbnode *bbnode, double *lower, int *leafcount,
        double *leafvals, tsp_bbnode **leafbbs)
{
    if ((bbnode->child0 == (tsp_bbnode *) NULL ||
         bbnode->child1 == (tsp_bbnode *) NULL) &&
        bbnode->lowerbound < *lower) {
        *lower = bbnode->lowerbound;
    }
    if (bbnode->status != BB_DONE) {
        if (leafvals != (double *) NULL) {
            leafvals[*leafcount] = bbnode->lowerbound;
        }
        if (leafbbs != (tsp_bbnode **) NULL) {
            leafbbs[*leafcount] = bbnode;
        }
        (*leafcount)++;
    }
    if (bbnode->child0) {
        collect_leaves (bbnode->child0, lower, leafcount, leafvals, leafbbs);
    }
    if (bbnode->child1) {
        collect_leaves (bbnode->child1, lower, leafcount, leafvals, leafbbs);
    }
}

static int prob_name (char *buf, size_t buflen, const char *f, int n)
{
    int l = (int) strlen(f);
    int lastslash;
    int i;
    int d;

    if (l + 5 > (int) buflen || n < 0) {
        fprintf (stderr, "Cannot generate filename for %s node %d\n",
                 f, n);
        return -1;
    }

    for (i = 0, lastslash = -1; i < l; i++) {
        if (f[i] == '/') lastslash = i;
        buf[i] = f[i];
    }
    if (l > lastslash + 9) l = lastslash + 9;
    for (i = lastslash+1; i < l; i++) {
        if (buf[i] == '.') buf[i] = '_';
    }
    if (n < 1000) {
        buf[l++] = '.';
        d = n/100;
        buf[l++] = '0' + ((unsigned int) d);
        n -= d*100;
        d = n/10;
        buf[l++] = '0' + d;
        n -= d*10;
        d = n;
        buf[l++] = '0' + ((unsigned int) d);
    } else if (n < 1000 + (26*36*36 - 5)) {
        buf[l++] = '.';
#define NAMESTRNUM(xa,xb,xc) (((xa)-'a') * 1296 + ((xb)-'a'+10) * 36 + \
                              ((xc)-'a'+10))
        n -= 1000;
        if (n >= NAMESTRNUM('m','a','s')) n++;
        if (n >= NAMESTRNUM('p','u','l')) n++;
        if (n >= NAMESTRNUM('r','e','s')) n++;
        if (n >= NAMESTRNUM('s','a','v')) n++;
        if (n >= NAMESTRNUM('s','o','l')) n++;
        d = n/1296;
        buf[l++] = 'a' + ((unsigned int) d);
        n -= d*1296;
        d = n/36;
        buf[l++] = (d < 10) ? '0' + ((unsigned int) d)
                            : 'a' + ((unsigned int) (d-10));
        n -= d*36;
        d = n;
        buf[l++] = (d < 10) ? '0' + ((unsigned int) d)
                            : 'a' + ((unsigned int) (d-10));
    } else if (n < 1000 + (26*36*36 - 5) + 36*36*36*36) {
        n -= 1000;
        n -= 26*36*36 - 5;
        d = n/(36*36*36);
        buf[l++] = (d < 10) ? '0' + ((unsigned int) d)
                            : 'a' + ((unsigned int) (d-10));
        n -= d*36*36*36;
        d = n/(36*36);
        buf[l++] = (d < 10) ? '0' + ((unsigned int) d)
                            : 'a' + ((unsigned int) (d-10));
        n -= d*(36*36);
        d = n/36;
        buf[l++] = (d < 10) ? '0' + ((unsigned int) d)
                            : 'a' + ((unsigned int) (d-10));
        n -= d*36;
        d = n;
        buf[l++] = (d < 10) ? '0' + ((unsigned int) d)
                            : 'a' + ((unsigned int) (d-10));
    } else {
        fprintf (stderr, "Node number %d too large\n", n);
        return -1;
    }
        
    buf[l] = '\0';
    return 0;
}

static void init_bbnode (tsp_bbnode *bbnode)
{
    bbnode->id         = 0;
    bbnode->lowerbound = 0.0;
    bbnode->status     = BB_NEEDS_CUTTING;
    bbnode->workstatus = BB_IDLE;
    bbnode->prev       = (tsp_bbnode *) NULL;
    bbnode->next       = (tsp_bbnode *) NULL;
    bbnode->parent     = (tsp_bbnode *) NULL;
    bbnode->child0     = (tsp_bbnode *) NULL;
    bbnode->child1     = (tsp_bbnode *) NULL;
    bbnode->cputime    = 0.0;
}

static void free_tree (tsp_bbnode **bbnode, CCptrworld *bbnode_world)
{
    if (!(*bbnode))  return;
    free_tree (&((*bbnode)->child0), bbnode_world);
    free_tree (&((*bbnode)->child1), bbnode_world);
    tsp_bbnode_free (bbnode_world, *bbnode);
    *bbnode = (tsp_bbnode *) NULL;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "d:glnp:P:?", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'd':
            arg_maxdepth = atoi(boptarg);
            break;
        case 'l':
            leafsummary = 1;
            break;
        case 'n':
            nodelist = 1;
            break;
        case 'p':
            dig_before = atoi (boptarg);
            break;
        case 'P':
            dig_after = atoi (boptarg);
            break;
        case 'g':
            use_graphics = 1;
            break;
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind >= ac) {
        usage (av[0]);
        return 1;
    }

    resfname = av[boptind++];

    if (boptind < ac) {
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] restart_fname\n", fname);
    fprintf (stderr, "   -d n  dump only to depth n\n");
    fprintf (stderr, "   -l    output leaf summary\n");
    fprintf (stderr, "   -n    output leaf node summary\n");
    fprintf (stderr, "   -p n  only show last n digits before decimal\n");
    fprintf (stderr, "   -P N  only show first n digits after decimal\n");
    fprintf (stderr, "   -g    create graphical (postscript) picture\n");
}


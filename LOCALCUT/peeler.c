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

typedef struct stripnode {
    int label;
    int done;
    int nadj;
    int *adj;
    int *coef;
} stripnode;


static int
    strip_outside (int nodecount, int edgecount, stripnode *nodelist,
                   int *elist, int *coef, int rhs, int out),
    strip_triangular (int nodecount, int edgecount, stripnode *nodelist,
                      int *elist, int *coef, int rhs),
    strip_graph (int nnodes, stripnode *nodes, int rhs),
    strip_maxclique (int nnodes, stripnode *nodes, int *maxclique,
                     int *work1),
    strip_graph_triangular (int nnodes, stripnode *nodes, int rhs),
    strip_clique_triangular (int nnodes, stripnode *nodes,
                             int *maxclique);

static void
    subtract_clique (stripnode *nodes, int maxsize, int *maxclique,
                     int label),
    subtract_clique_outside (int nnodes, stripnode *nodes, int maxsize,
                             int *maxclique);


int main (CC_UNUSED int ac, CC_UNUSED char **av)
{
    int nodecount;
    int edgecount;
    int *elist = (int *) NULL;
    int *coef = (int *) NULL;
    int *adjspace = (int *) NULL;
    int *coefspace = (int *) NULL;
    stripnode *nodelist = (stripnode *) NULL;
    int rhs;
    int i;
    int rval;

    scanf ("%d%d", &nodecount, &edgecount);

    elist = CC_SAFE_MALLOC (edgecount*2, int);
    coef = CC_SAFE_MALLOC (edgecount, int);
    nodelist = CC_SAFE_MALLOC (nodecount, stripnode);
    adjspace = CC_SAFE_MALLOC (nodecount * nodecount, int);
    coefspace = CC_SAFE_MALLOC (nodecount * nodecount, int);
    if (elist == (int *) NULL ||
        coef == (int *) NULL ||
        nodelist == (stripnode *) NULL ||
        adjspace == (int *) NULL ||
        coefspace == (int *) NULL) {
        fprintf (stderr, "Out of memory\n");
        rval = -1;
        goto CLEANUP;
    }

    for (i=0; i<nodecount; i++) {
        nodelist[i].adj = adjspace + i*nodecount;
        nodelist[i].coef = coefspace + i*nodecount;
    }

    for (i=0; i<edgecount; i++) {
        scanf ("%d%d%d",&elist[2*i], &elist[2*i+1], &coef[i]);
        coef[i] = -coef[i];
    }

    scanf ("%d", &rhs);
    rhs = -rhs;

    for (i=0; i<nodecount; i++) {
        printf ("With %d outside:\n", i);
        strip_outside (nodecount, edgecount, nodelist, elist, coef, rhs, i);
        printf ("\n");
    }

    printf ("With tight triangular:\n");
    strip_triangular (nodecount, edgecount, nodelist, elist, coef, rhs);

    rval = 0;

  CLEANUP:
    CC_IFFREE (elist, int);
    CC_IFFREE (coef, int);
    CC_IFFREE (nodelist, stripnode);
    CC_IFFREE (adjspace, int);
    CC_IFFREE (coefspace, int);
    return rval;
}

static int strip_outside (int nodecount, int edgecount, stripnode *nodelist,
                          int *elist, int *coef, int rhs, int out)
{
    int i;
    int j;
    int minc;
    int xdiv;
    int c;

    for (i=0; i<nodecount; i++) {
        nodelist[i].label = 0;
        for (j=0; j<nodecount; j++) {
            nodelist[i].adj[j] = j;
            nodelist[i].coef[j] = 0;
        }
    }

    for (i=0; i<edgecount; i++) {
        nodelist[elist[2*i]].coef[elist[2*i+1]] = coef[i];
        nodelist[elist[2*i+1]].coef[elist[2*i]] = coef[i];
    }

    for (i=0; i<nodecount; i++) {
        if (i != out) {
            c = -nodelist[i].coef[out];
            for (j=0; j<nodecount; j++) {
                nodelist[i].coef[j] += c;
                nodelist[j].coef[i] += c;
            }
            rhs += 2*c;
        }
    }

    minc = 1000000;
    for (i=0; i<nodecount; i++) {
        if (i != out) {
            for (j=0; j<nodecount; j++) {
                if (j != out && i != j && nodelist[i].coef[j] < minc) {
                    minc = nodelist[i].coef[j];
                }
            }
        }
    }

#if 0
    for (i=0; i<nodecount; i++) {
        for (j=i+1; j<nodecount; j++) {
            if (nodelist[i].coef[j]) {
                printf ("%d %d %d\n", i, j, nodelist[i].coef[j]);
            }
            if (nodelist[i].coef[j] != nodelist[j].coef[i]) {
                printf ("ERROR: n[%d].c[%d]=%d n[%d].c[%d]=%d\n",
                        i,j,nodelist[i].coef[j],
                        j,i,nodelist[j].coef[i]);
            }
        }
    }

    printf ("minc %d\n", minc);
#endif /* 0 */

    for (i=0; i<nodecount; i++) {
        if (i != out) {
            for (j=0; j<nodecount; j++) {
                if (j != out) {
                    nodelist[i].coef[j] -= minc;
                }
            }
        }
    }
    rhs -= (nodecount-2) * minc;

    for (i=0; i<nodecount; i++) {
        nodelist[i].coef[i] = 0;
    }

#if 0
    for (i=0; i<nodecount; i++) {
        for (j=i+1; j<nodecount; j++) {
            if (nodelist[i].coef[j]) {
                printf ("%d %d %d\n", i, j, nodelist[i].coef[j]);
            }
            if (nodelist[i].coef[j] != nodelist[j].coef[i]) {
                printf ("ERROR: n[%d].c[%d]=%d n[%d].c[%d]=%d\n",
                        i,j,nodelist[i].coef[j],
                        j,i,nodelist[j].coef[i]);
            }
        }
    }
#endif /* 0 */

    for (i=0; i<nodecount; i++) {
        nodelist[i].nadj = 0;
        for (j=0; j<nodecount; j++) {
            if (nodelist[i].coef[j]) {
                nodelist[i].adj[nodelist[i].nadj] = j;
                nodelist[i].coef[nodelist[i].nadj] = nodelist[i].coef[j];
                nodelist[i].nadj++;
            }
        }
    }

    xdiv = rhs;
    for (i=0; i<nodecount; i++) {
        for (j=0; j<nodelist[i].nadj; j++) {
            xdiv = CCutil_our_gcd (xdiv, nodelist[i].coef[j]);
        }
    }

    rhs /= xdiv;
    for (i=0; i<nodecount; i++) {
        for (j=0; j<nodelist[i].nadj; j++) {
            nodelist[i].coef[j] /= xdiv;
        }
    }

    return strip_graph (nodecount, nodelist, rhs);
}

static int strip_graph (int nnodes, stripnode *nodes, int rhs)
{
    int deg;
    int *maxclique = (int *) NULL;
    int *work1 = (int *) NULL;
    int i;
    int rval = 0;
    int nrhs = -2 * rhs;

    maxclique = CC_SAFE_MALLOC (nnodes, int);
    work1 = CC_SAFE_MALLOC (nnodes, int);
    if (maxclique == (int *) NULL ||
        work1 == (int *) NULL) {
        rval = -1;
        goto CLEANUP;
    }

    while (1) {
        deg = strip_maxclique (nnodes, nodes, maxclique, work1);
        if (deg <= 1) {
            CC_FREE (work1, int);
            CC_FREE (maxclique, int);
            printf (" >= %d\n", nrhs);
            return 0;
        }
        for (i=0; i<deg; i++) {
            printf ("%d ", maxclique[i]);
        }
        printf ("\n");
        nrhs += 2*deg;
    }

  CLEANUP:
    CC_IFFREE (work1, int);
    CC_IFFREE (maxclique, int);
    return rval;
}

static int strip_maxclique (int nnodes, stripnode *nodes, int *maxclique,
                            int *work1)
{
    int simpsize = 0;
    int nonsize = 0;
    int i, ia, j, ja, k;
    int label = 0;
    int mindeg;
    int size;

    for (i=0; i<nnodes; i++) {
        nodes[i].label = 0;
        nodes[i].done = 0;
    }
    for (i=0; i<nnodes; i++) {
        if (!nodes[i].done) {
            mindeg = nodes[i].nadj;
            size = 1;
            label++;
            nodes[i].label = label;
            for (ia=0; ia<nodes[i].nadj; ia++) {
                j = nodes[i].adj[ia];
                nodes[j].label = label;
            }
            for (ia=0; ia<nodes[i].nadj; ia++) {
                j = nodes[i].adj[ia];
                if (nodes[j].label == label) {
                    size++;
                    if (nodes[j].nadj < mindeg) mindeg = nodes[j].nadj;
                    label++;
                    nodes[j].label = label;
                    for (ja=0; ja<nodes[j].nadj; ja++) {
                        k = nodes[j].adj[ja];
                        if (nodes[k].label == label-1) {
                            nodes[k].label = label;
                        }
                    }
                }
            }
            if (mindeg == size-1 && size > simpsize) {
                simpsize = size;
                maxclique[0] = i;
                for (ia=0, k=1; ia<nodes[i].nadj; ia++) {
                    j = nodes[i].adj[ia];
                    if (nodes[j].label == label) {
                        maxclique[k++] = j;
                    }
                }
                if (size != k) {
                    fprintf (stderr, "Size mismatch 1\n");
                    simpsize = k;
                }
            } else if (mindeg > size-1 && size > nonsize) {
                nonsize = size;
                work1[0] = i;
                for (ia=0, k=1; ia<nodes[i].nadj; ia++) {
                    j = nodes[i].adj[ia];
                    if (nodes[j].label == label) {
                        work1[k++] = j;
                    }
                }
                if (size != k) {
                    fprintf (stderr, "Size mismatch 2\n");
                    nonsize = k;
                }
                for (ia=0; ia<nodes[i].nadj; ia++) {
                    j = nodes[i].adj[ia];
                    if (nodes[j].label == label) {
                        nodes[j].done = 1;
                    }
                }
            }
        }
    }
    if (simpsize <= 1 && nonsize > 1) {
        for (i=0; i<nonsize; i++) {
            maxclique[i] = work1[i];
        }
        simpsize = nonsize;
    }

    label++;
    subtract_clique (nodes, simpsize, maxclique, label);

    return simpsize;
}

static void subtract_clique (stripnode *nodes, int maxsize, int *maxclique,
                             int label)
{
    int i, j;
    int k;

    for (i=0; i<maxsize; i++) {
        nodes[maxclique[i]].label = label;
    }

    for (i=0; i<maxsize; i++) {
        j = maxclique[i];
        k = 0;
        while (k < nodes[j].nadj) {
            if (nodes[nodes[j].adj[k]].label == label) {
                nodes[j].coef[k]--;
            }
            if (nodes[j].coef[k]) {
                k++;
            } else {
                nodes[j].nadj--;
                nodes[j].adj[k] = nodes[j].adj[nodes[j].nadj];
                nodes[j].coef[k] = nodes[j].coef[nodes[j].nadj];
            }
        }
    }
}

static int strip_triangular (int nodecount, int edgecount, stripnode *nodelist,
                             int *elist, int *coef, int rhs)
{
    int i;
    int j;
    int k;
    int minc;
    int xdiv;
    int c;

    for (i=0; i<nodecount; i++) {
        nodelist[i].label = 0;
        for (j=0; j<nodecount; j++) {
            nodelist[i].adj[j] = j;
            nodelist[i].coef[j] = 0;
        }
    }

    for (i=0; i<edgecount; i++) {
        nodelist[elist[2*i]].coef[elist[2*i+1]] = -coef[i];
        nodelist[elist[2*i+1]].coef[elist[2*i]] = -coef[i];
    }
    rhs = -rhs;

    for (i=0; i<nodecount; i++) {
        for (j=0; j<nodecount; j++) {
            nodelist[i].coef[j] *= 2;
        }
    }
    rhs *= 2;

    for (i=0; i<nodecount; i++) {
        minc = 1000000000;
        for (j=0; j<nodecount; j++) {
            for (k=0; k<nodecount; k++) {
                if (i != j && i != k && j != k) {
                    c = nodelist[i].coef[k] + nodelist[i].coef[j] -
                        nodelist[j].coef[k];
                    if (c < minc) minc = c;
                }
            }
        }
        if (minc % 2 != 0) {
            fprintf (stderr, "Whoa, parity problem in strip_triangular\n");
            return -1;
        }
        minc /= 2;
        for (j=0; j<nodecount; j++) {
            nodelist[i].coef[j] -= minc;
            nodelist[j].coef[i] -= minc;
        }
        rhs -= 2*minc;
    }

#if 1
    for (i=0; i<nodecount; i++) {
        for (j=i+1; j<nodecount; j++) {
            if (nodelist[i].coef[j]) {
                printf ("%d %d %d\n", i, j, nodelist[i].coef[j]);
            }
            if (nodelist[i].coef[j] != nodelist[j].coef[i]) {
                printf ("ERROR: n[%d].c[%d]=%d n[%d].c[%d]=%d\n",
                        i,j,nodelist[i].coef[j],
                        j,i,nodelist[j].coef[i]);
            }
        }
    }
    printf (">= %d\n", rhs);
#endif /* 0 */

    for (i=0; i<nodecount; i++) {
        nodelist[i].coef[i] = 0;
    }

    xdiv = rhs;
    for (i=0; i<nodecount; i++) {
        for (j=0; j<nodecount; j++) {
            xdiv = CCutil_our_gcd (xdiv, nodelist[i].coef[j]);
        }
    }

    rhs /= xdiv;
    for (i=0; i<nodecount; i++) {
        for (j=0; j<nodecount; j++) {
            nodelist[i].coef[j] /= xdiv;
        }
    }

    return strip_graph_triangular (nodecount, nodelist, rhs);
}

static int strip_graph_triangular (int nnodes, stripnode *nodes, int rhs)
{
    int deg;
    int *maxclique = (int *) NULL;
    int i;
    int rval = 0;

    maxclique = CC_SAFE_MALLOC (nnodes, int);
    if (maxclique == (int *) NULL) {
        rval = -1;
        goto CLEANUP;
    }

    while (1) {
        deg = strip_clique_triangular (nnodes, nodes, maxclique);
        if (deg == 0) {
            printf (" >= %d\n", rhs);
            rval = 0;
            goto CLEANUP;
        } else if (deg < 0) {
            fprintf (stderr, "strip_clique_triangular failed\n");
            rval = deg;
            goto CLEANUP;
        }
        for (i=0; i<deg; i++) {
            printf ("%d ", maxclique[i]);
        }
        printf ("\n");
    }

  CLEANUP:
    CC_IFFREE (maxclique, int);
    return rval;
}

static int strip_clique_triangular (int nnodes, stripnode *nodes,
                                    int *maxclique)
{
    int i;
    int j;
    int min_coef;
    int min_i;
    int min_j;
    int xdiv;
    int size;
    int a;
    int b;

    min_coef = 0;
    xdiv = 0;
    min_i = -1;
    min_j = -1;

    for (i=0; i<nnodes; i++) {
        for (j=0; j<nnodes; j++) {
            if (nodes[i].coef[j] != 0) {
                xdiv = CCutil_our_gcd (xdiv, nodes[i].coef[j]);
                if (min_i == -1 || nodes[i].coef[j] < min_coef) {
                    min_coef = nodes[i].coef[j];
                    min_i = i;
                    min_j = j;
                }
            }
        }
    }

    if (min_coef == 0) {
        return 0;
    }

    if (xdiv > 1) {
        printf ("Scaling down by %d\n", xdiv);
        for (i=0; i<nnodes; i++) {
            for (j=0; j<nnodes; j++) {
                nodes[i].coef[j] /= xdiv;
            }
        }
        min_coef /= xdiv;
    }

    if (min_coef > 1) {
        printf ("lost, minimum coefficient %d\n", min_coef);
        return -1;
    }

    size = 0;
    for (i=0; i<nnodes; i++) {
        a = nodes[i].coef[min_i];
        b = nodes[i].coef[min_j];
        if (a == b+1) {
            maxclique[size++] = i;
        } else if (a != b-1) {
            printf ("lost, triangle %d, %d, 1\n", a, b);
            return -2;
        }
    }

    subtract_clique_outside (nnodes, nodes, size, maxclique);

    return size;
}

static void subtract_clique_outside (int nnodes, stripnode *nodes, int maxsize,
        int *maxclique)
{
    int i, j;

    for (i=0; i<nnodes; i++) {
        nodes[i].label = 0;
    }
    for (i=0; i<maxsize; i++) {
        nodes[maxclique[i]].label = 1;
    }

    for (i=0; i<nnodes; i++) {
        for (j=0; j<nnodes; j++) {
            if (nodes[i].label != nodes[j].label) {
                nodes[i].coef[j]--;
            }
        }
    }
}

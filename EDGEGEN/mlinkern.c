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
/*           CHAINED LIN-KERNIGHAN FOR PERFECT MATCHINGS                    */
/*                                                                          */
/*                                                                          */
/*  Written by:  Andre Rohe (Bonn)                                          */
/*  Date: September, 1996                                                   */
/*  Changes: 11.22.1996 (bico) static functions, error handling             */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCedgegen_mlinkern (int ncount, CCdatagroup *dat, int wantlist,     */
/*      int *ecount, int **elist, CCkdtree *kt, int iterations,             */
/*      CCrandstate *rstate)                                                */
/*    COMPUTES the edge set of union of Chained Lin-Kernighan matchings.    */
/*       -ncount is the number of nodes in the graph                        */
/*       -dat gives the data to generate edge lengths                       */
/*       -if wantlist is nonzero, then the edgeset will be returned         */
/*       -ecount will return the number of edges found (if wantlist)        */
/*       -elist will return the matching in end end format (if wantlist)    */
/*       -kt is a pointer to a kd-tree                                      */
/*       -iterations is the number of Chained LK kicks matchings            */
/*                                                                          */
/*    NOTES:                                                                */
/*      Code is only set up for norms of CC_KD_NORM_TYPE.                   */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "kdtree.h"
#include "mlinkern.h"
#include "util.h"

#define MY_DOUBLE_MIN -999999999.0
#define MY_DOUBLE_MAX  999999999.0
#define MY_INT_MIN    -999999999
#define MY_INT_MAX     999999999
#define CutOff 4
#define CC_SWAP(x,y,t) ((t)=(x),(x)=(y),(y)=(t))
#define MAX_EXTRA_EDGES 100000

#define W_MAX 25        /* Max. Anzahl von w0 & w1, die moeglich ist */
#define Opt_Edg_Nod 10  /* Kanten in Graph, die gecheckt werden fuer Lin-Ker */

typedef struct node {
    int neigh;
    int weight;
} node;

typedef struct optedge {
    int weight;             /* weight of edge */
    int numbout;            /* edge to numbout from the vertex we came from */
    struct optedge *next;
} optedge;

typedef struct flip_coor {
    int a;   /* the 4 flip points */
    int b;
    int c;
    int d;
    int d1;
    int d2;
} flip_coor;

typedef struct w_T{
    int w;
    double weight;
} w_T;

typedef struct edge {
    int numbout;          /* edge to numbout from the vertex we came from */
    struct edge *next;
} edge;

typedef struct data {
    int nnodes;          /* The number of nodes in the Graph */
    int iteration;
    w_T w0[W_MAX], w1[W_MAX];
    optedge *edge_mem;
    optedge **adj;
    optedge *end;
    int *u;
    int *w;
    int v;
    int depth;
    int new_best;
    int print_level;
    CCdatagroup *dat;
} data;

typedef struct graph {
    edge *edge_graph,**adj_graph,*end_graph;
    int edge_count;
} graph;

typedef struct mqueue {
    int   in_Magic;
    int  *in;
    int  *ok_q;             /* Queue with all nodes which must be checked */
    int   ok_head, ok_tail; /* end and start for ok_q */
    char *ok;               /* Array, whether node i 2-opt optimal */
} mqueue;

typedef struct matchings{
    double p_weight;
    double match_weight;
    double best_match_weight;
    double very_best_match_weight;
    node  *match, *very_best_match;
} matchings;

typedef struct tflip {
    int max_deep;
    flip_coor *flip_stack;
    flip_coor *flip_best_stack;
    flip_coor *flip_all_stack;
    int flip_stack_count;
    int flip_best_stack_count;
    int flip_all_stack_count;
} tflip;



static void
    frei (data *D, tflip *fli, mqueue *qu),
    put_q (data *D, mqueue *qu, int v),
    init_q (mqueue *qu),
    flip (data *D, node *m, int a, int b, int c, int d),
    flipp (data *D, matchings *M, tflip *fli, int a, int b, int c, int d),
    unflipp (matchings *M, tflip *fli, int a, int b, int c, int d, int d1,
             int d2),
    new_opt_edge (data *D, int x, int y, int w, int z),
    init_queue (data *D, mqueue *qu),
    find_w0 (data *D, matchings *M, mqueue *qu, int u_0, int v),
    find_w1 (data *D, matchings *M, mqueue *qu, int u_1),
    edge_scan2 (data *D, matchings *M, tflip *fli, mqueue *qu),
    step345 (data *D, matchings *M, tflip *fli, mqueue *qu),
    lin_ker_do (data *D, matchings *M, tflip *fli, mqueue *qu);

static int
    in_graph (graph *G, int i, int j),
    init_graph (data *D, graph *G),
    get_q (mqueue *qu),
    non_empty_q (mqueue *qu),
    init (data *D, tflip *fli, mqueue *qu),
    compare_func (const void *e1, const void *e2),
    find_w (data *D, matchings *M, mqueue *qu, int u_depth, int w_depth_1),
    edge_scan1 (data *D, matchings *M, tflip *fli, mqueue *qu),
    restore_match (data *D, matchings *M, graph *G, tflip *fli),
    new_edge (data *D, graph *G, int x, int y),
    opt_makegraph (CCkdtree *kt, CCdatagroup *dat, data *D,
        CCrandstate *rstate),
    m_lin_ker (CCkdtree *kt, CCdatagroup *dat, data *D, matchings *M,
        graph *G, tflip *fli,mqueue *qu,int max_iteration,
        CCrandstate *rstate);

static node
    *nn_match (CCkdtree *kt, CCdatagroup *dat, data *D, matchings *M,
        CCrandstate *rstate);




/* ************************************************************/
/*                            ende                            */
/* ************************************************************/

static void frei (data *D, tflip *fli, mqueue *qu)
{
    CC_IFFREE (qu->in, int);
    CC_IFFREE (qu->ok, char);
    CC_IFFREE (qu->ok_q, int);

    CC_IFFREE (D->edge_mem, optedge);
    CC_IFFREE (D->adj, optedge*);
    CC_IFFREE (D->end, optedge);
    CC_IFFREE (D->u, int);
    CC_IFFREE (D->w, int);
    CC_IFFREE (fli->flip_stack, flip_coor);
    CC_IFFREE (fli->flip_all_stack, flip_coor);
    CC_IFFREE (fli->flip_best_stack, flip_coor);
}

/* ************************************************************/
/*               Nearest Neighbour Matching                   */
/* ************************************************************/

static node *nn_match (CCkdtree *kt, CCdatagroup *dat, data *D, matchings *M,
        CCrandstate *rstate)
{
    int i,nod;
    int sum;
    int nei;

    if (D->print_level)
        printf("\n  Start Nearest Neightbour Matching ...");

    M->very_best_match = CC_SAFE_MALLOC (D->nnodes, node);
    M->match           = CC_SAFE_MALLOC (D->nnodes, node);
    if (!M->very_best_match || !M->match) {
        fprintf (stderr, "out of memory in nn_match\n");
        return (node *) NULL;
    }

    CCkdtree_undelete_all (kt, D->nnodes);

    for (i = 0; i < D->nnodes; i++)
        M->match[i].neigh=-1;

    for (i = 0; i < D->nnodes/2; i++) {
        for (nod = CCutil_lprand(rstate) %D->nnodes; M->match[nod].neigh != -1;
             nod = (nod+1) % D->nnodes);
        nei=CCkdtree_node_nearest (kt, nod, dat, NULL);
        CCkdtree_delete (kt, nod);
        CCkdtree_delete (kt, nei);
        M->match[nod].neigh = nei;
        M->match[nod].weight = CCutil_dat_edgelen(nod, nei, dat);
        M->match[nei].neigh = nod;
        M->match[nei].weight = CCutil_dat_edgelen(nod, nei, dat);
    }

    if (D->print_level) printf("\n  ... Ready !! ");
    sum=0.0;
    for (i = 0; i < D->nnodes; i++)
        sum += M->match[i].weight;
    sum /= 2;
    if (D->print_level)
        printf("\n  The Nearest Neighbour Matching has weight %i",sum);

    CCkdtree_undelete_all (kt, D->nnodes);

    return (M->match);
}

static int in_graph (graph *G, int i, int j)
{
    edge *e;

    for (e = G->adj_graph[i]; e != G->end_graph; e = e->next) {
       if (e->numbout == j)
           return 1;
    }

    return 0;
}

static int new_edge (data *D, graph *G, int x, int y)
{
    edge *e;

    if (G->edge_count >= D->nnodes + MAX_EXTRA_EDGES) {
        fprintf (stderr, "\nToo many edges (%i)\n", G->edge_count);
        return 1;
    }
    e = G->edge_graph+G->edge_count;

    e->numbout = y;

    e->next = G->adj_graph[x];

    G->adj_graph[x] = e;
    G->edge_count++;
    return 0;
}

static int init_graph (data *D, graph *G)
{
    int nedges = D->nnodes + MAX_EXTRA_EDGES;
    int i;

    G->edge_graph = CC_SAFE_MALLOC (nedges, edge);
    G->end_graph  = CC_SAFE_MALLOC (1, edge);
    G->adj_graph  = CC_SAFE_MALLOC (D->nnodes, edge*);
    if (!G->edge_graph || !G->end_graph || !G->adj_graph) {
        fprintf (stderr, "out of memory in init_graph\n");
        return 1;
    }
    for (i = 0; i < D->nnodes; i++)
        G->adj_graph[i] = G->end_graph;
    G->edge_count=0;
    return 0;
}

static void put_q (data *D, mqueue *qu, int v)
{
    qu->ok_q[qu->ok_tail++] = v;
    qu->ok_tail %= (D->nnodes + 1);
}

static int get_q (mqueue *qu)
{
    return qu->ok_q[--qu->ok_tail];
}

static void init_q (mqueue *qu)
{
    qu->ok_tail = 0;
    qu->ok_head = 0;
}

static int non_empty_q (mqueue *qu)
{
    return qu->ok_head != qu->ok_tail;
}

static void flip (data *D, node *m, int a, int b, int c, int d)
{
    m[a].neigh = c;
    m[c].weight = m[a].weight = CCutil_dat_edgelen(a, c, D->dat);
    m[c].neigh = a;
    m[b].neigh = d;
    m[d].weight = m[b].weight = CCutil_dat_edgelen(b, d, D->dat);
    m[d].neigh = b;
}

/* ************************************************************/
/*                            flip                           */
/*             a,b,c,d sind die Punkte in the matching        */
/* ************************************************************/

static void flipp (data *D, matchings *M, tflip *fli, int a, int b, int c,
                   int d)
{
    int i;
    int d1 = CCutil_dat_edgelen (a, c, D->dat);
    int d2 = CCutil_dat_edgelen (b, d, D->dat);

    M->match_weight += d1 + d2-M->match[a].weight - M->match[c].weight;

    fli->flip_stack[fli->flip_stack_count].a = a;
    fli->flip_stack[fli->flip_stack_count].b = b;
    fli->flip_stack[fli->flip_stack_count].c = c;
    fli->flip_stack[fli->flip_stack_count].d = d;
    fli->flip_stack[fli->flip_stack_count].d1 = M->match[a].weight;
    fli->flip_stack[fli->flip_stack_count].d2 = M->match[c].weight;
    fli->flip_stack_count++;

    M->match[a].neigh = c;
    M->match[c].weight = (M->match[a].weight = d1);
    M->match[c].neigh = a;
    M->match[b].neigh = d;
    M->match[d].weight = (M->match[b].weight = d2);
    M->match[d].neigh = b;

    if (M->match_weight < M->best_match_weight) {
        D->new_best = 1;
        for (i = 0; i < fli->flip_stack_count; i++) {
            fli->flip_best_stack[i] = fli->flip_stack[i];
        }
        fli->flip_best_stack_count = fli->flip_stack_count;
        M->best_match_weight = M->match_weight;
    }
}

/* ************************************************************/
/*                            unflipp                         */
/*  this routine needs wb=wa+1, wd=wc+1 (%D->nnodes)          */
/*             a,b,c,d sind die Punkte in match               */
/* ************************************************************/

static void unflipp (matchings *M, tflip *fli, int a, int b, int c, int d,
                     int d1, int d2)
{
    fli->flip_stack_count--;

    M->match_weight += d1 + d2 - M->match[a].weight - M->match[c].weight;

    M->match[a].neigh = c;
    M->match[c].weight = (M->match[a].weight = d1);
    M->match[c].neigh = a;
    M->match[b].neigh = d;
    M->match[d].weight = (M->match[b].weight = d2);
    M->match[d].neigh = b;
}

static int restore_match (data *D, matchings *M, graph *G, tflip *fli)
{
    int i;

    if (fli->flip_all_stack_count > fli->max_deep) {
        fprintf (stderr, "\nfli->flip_all_stack_count = %i, over limit\n",
            fli->flip_all_stack_count);
        return 1;
    }

    for (i = 0; i < D->nnodes; i++) {
        if (i < M->match[i].neigh && in_graph (G, i, M->match[i].neigh) == 0) {
            if (new_edge (D, G, i, M->match[i].neigh)) {
                fprintf (stderr, "new_edge failed\n");
                return 1;
            }
        }
    }

    /* New Very Best match in LinKer ? */

    if (M->match_weight < M->very_best_match_weight) {
        M->very_best_match_weight = M->match_weight;
        memcpy((char *) M->very_best_match, (char *) M->match,
               D->nnodes*sizeof(node));
    } else {
        memcpy((char *) M->match, (char *) M->very_best_match,
               D->nnodes*sizeof(node));
        M->match_weight = M->very_best_match_weight;
    }
    return 0;
}

static void new_opt_edge (data *D, int x, int y, int w, int z)
{
    optedge *e, *e1, *e2 = (optedge *) NULL;

    e = D->edge_mem + z;

    e->numbout = y;
    e->weight = w;

    if (D->adj[x] == D->end) {   /* Die erste Kante */
        e->next = D->adj[x];
        D->adj[x] = e;
        return;
    }

    /* put in new edge in correct order; D->adj[x]->weight should be minimal */

    for (e1 = D->adj[x]; e1 != D->end; e1 = e1->next) {
        if (e1->weight >= w) {
            e->next = e1;
            if (e1 == D->adj[x])
                D->adj[x] = e;
            else
                e2->next = e;
            return;
        }
        e2 = e1;
    }
    e2->next = e;
    e->next = D->end;
}

static int opt_makegraph (CCkdtree *kt, CCdatagroup *dat, data *D,
        CCrandstate *rstate)
{
    int i,j;
    int ary[2*Opt_Edg_Nod];
    int nr = 0;

    D->end = CC_SAFE_MALLOC (1, optedge);
    D->adj = CC_SAFE_MALLOC (D->nnodes, optedge *);
    if (!D->end || !D->adj) {
        fprintf (stderr, "out of memory in opt_makegraph\n");
        CC_IFFREE (D->end, optedge);
        CC_IFFREE (D->adj, optedge *);
        return 1;
    }

    for (i = 0; i < D->nnodes; i++) {
        D->adj[i] = D->end;
    }

    CCkdtree_undelete_all (kt, D->nnodes);

    /* Edg_Opt_Nod Nearest Neighbours ! */

    if (D->print_level) {
        printf("\n  Ich build the %i Nearest Neighbour Graph ..", Opt_Edg_Nod);
        fflush(stdout);
    }

    D->edge_mem = CC_SAFE_MALLOC (Opt_Edg_Nod*D->nnodes, optedge);
    if (!D->edge_mem) {
        fprintf (stderr, "out of memory in opt_makegraph\n");
        CC_FREE (D->end, optedge);
        CC_FREE (D->adj, optedge *);
        return 1;
    }

    for (i = 0; i < D->nnodes; i++) {

        /* ************************************************** */
        /* Graph in Lin-Kernighan                             */
        /*  Einfach die Naechsten Nachbarn                    */
        /* ************************************************** */

        /* Graph mit Naechsten Nachbarn Irgendwo */
        CCkdtree_node_k_nearest (kt, D->nnodes, i, Opt_Edg_Nod, dat, NULL,
                                 ary, rstate);
        for (j = 0; j < Opt_Edg_Nod; j++) {
            new_opt_edge(D, i, ary[j], CCutil_dat_edgelen (i, ary[j],
                                                           D->dat), nr);
            nr++;
        }
    }

    if (D->print_level) {
        printf(" ... ready !! "); fflush(stdout);
    }
    return 0;
}

static void init_queue (data *D, mqueue *qu)
{
    int i, j;

    init_q (qu);

    j = 0;
    for (i = 0; i < D->nnodes; i++)
        qu->ok[i] = 1;
    for (i = 0;i < D->nnodes; i++) {
        qu->in[i] = -1;
        qu->ok[i] = 0;
        put_q (D, qu, i);
        j++;
    }
    qu->in_Magic=1;
}

static int init (data *D, tflip *fli, mqueue *qu)
{
    if (D->nnodes < 1000)
        fli->max_deep = 10*D->nnodes;
    else
        fli->max_deep = D->nnodes+20000;
    fli->flip_all_stack  = CC_SAFE_MALLOC (fli->max_deep, flip_coor);
    fli->flip_best_stack = CC_SAFE_MALLOC (fli->max_deep, flip_coor);
    fli->flip_stack      = CC_SAFE_MALLOC (fli->max_deep, flip_coor);

    D->u = CC_SAFE_MALLOC (fli->max_deep, int);
    D->w = CC_SAFE_MALLOC (fli->max_deep, int);

    qu->ok   = CC_SAFE_MALLOC (D->nnodes, char);
    qu->ok_q = CC_SAFE_MALLOC (D->nnodes+1, int);
    qu->in   = CC_SAFE_MALLOC (D->nnodes, int);

    if (!fli->flip_all_stack || !fli->flip_best_stack || !fli->flip_stack ||
        !D->u || !D->w || !qu->ok || !qu->ok_q || !qu->in) {
        fprintf (stderr, "out of memory in init\n");
        CC_IFFREE (fli->flip_all_stack, flip_coor);
        CC_IFFREE (fli->flip_best_stack, flip_coor);
        CC_IFFREE (fli->flip_stack, flip_coor);
        CC_IFFREE (D->u, int);
        CC_IFFREE (D->w, int);
        CC_IFFREE (qu->ok, char);
        CC_IFFREE (qu->ok_q, int);
        CC_IFFREE (qu->in, int);
        return 1;
    }
    return 0;
}

static int compare_func (const void *e1, const void *e2)
{
    int v1 = ((const w_T *)e1)->weight;
    int v2 = ((const w_T *)e2)->weight;
    return (v1 < v2) ? -1 : (v1 > v2) ? 1 : 0;
}

/* ************************************************************/
/*                          find_w0                           */
/*                  normal mit find_w0(u[0],v)                */
/* ************************************************************/

static void find_w0 (data *D, matchings *M, mqueue *qu, int u_0, int v)
{
    int u_1, count = 0;
    double d;
    optedge *e;

    d = M->match_weight - M->match[u_0].weight;
    e = D->adj[u_0];

    while (e != D->end && e->weight + d < M->best_match_weight) {
        /* is new p_weight <= M->best_match_weight ? */
        if (qu->in[u_1 = M->match[e->numbout].neigh] != qu->in_Magic)  {
            /* Man sucht w[0] und minimiert daher                      */
            /*  dist(D,u[0],w[0])-dist(D,u[1],w[0])                    */
            /* Wobei:                                                  */
            /*         e-> numbout                = w[0]               */
            /*         M->match[e->numbout].neigh = u[1]               */
            /*         u_0                        = u[0]               */
            if (e->numbout != v) {
                D->w0[count].weight = e->weight-M->match[u_1].weight;
                D->w0[count].w = e->numbout;
                count++;
            }
        }
        e=e->next;
    }
    D->w0[count].weight = MY_INT_MAX;
    D->w0[count].w = -1;
    if (count>1)
        qsort((void*)D->w0,(size_t)count,sizeof(w_T),compare_func);
}

/* ************************************************************/
/*                          find_w1                           */
/*           normalerweise mit find_w1(u[1])                  */
/* ************************************************************/

static void find_w1 (data *D, matchings *M, mqueue *qu, int u_1)
{
    int u_2, count = 0;
    double d;
    optedge *e;

    d = M->match_weight-M->match[u_1].weight;
    e = D->adj[u_1];

    while (e != D->end && (e->weight + d < M->best_match_weight)) {
        /* is new p_weight <= M->best_match_weight ? */
        if (qu->in[u_2 = M->match[e->numbout].neigh] != qu->in_Magic
            /* && e->numbout!=qu->in_Magic */) {
            /* Man sucht w[1] und minimiert daher                */
            /* dist(D,u[1],w[1])-dist(D,u[2],w[1])               */
            /* Wobei:                                            */
            /*      e-> numbout                = w[1]            */
            /*      M->match[e->numbout].neigh = u[2]            */
            /*      u_1                        = u[1]            */
            if (e->numbout != D->v) {
                D->w1[count].weight = e->weight - M->match[u_2].weight;
                D->w1[count].w = e->numbout;
                count++;
            }
        }
        e = e->next;
    }
    D->w1[count].weight = MY_INT_MAX;
    D->w1[count].w = -1;
    if (count > 1)
        qsort ((void*)D->w1, (size_t) count, sizeof(w_T), compare_func);
}

/* *****************************************************************/
/*                          find_w                                 */
/*  normaler Aufruf: w[D->depth]=find_w(u[D->depth],w[D->depth-1]) */
/*******************************************************************/

static int find_w (data *D, matchings *M, mqueue *qu, int u_depth,
        int w_depth_1)
{
    int minw, u_depth_1;
    double help, min, d;
    optedge *e;

    d = M->match_weight - M->match[u_depth].weight;
    min = MY_INT_MAX;
    minw = -1;

    /* Man sucht w[D->depth] und minimiert daher                           */
    /* dist(D,u[D->depth],w[D->depth])-dist(D,u[D->depth+1],w[D->depth])   */
    /* Wobei:                                                              */
    /*         e-> numbout                   = w[D->depth]                 */
    /*         M->match[e->numbout].neigh    = u[D->depth+1]               */
    /*         u_depth                       = u[D->depth]                 */
    /*         w_depth_1                     = w[D->depth-1]               */

    if (D->depth > 8) { /* Einstellung */
        int test = 0;
        for (e = D->adj[u_depth]; e != D->end &&
             e->weight + d < M->best_match_weight; e=e->next) {
            /* is new p_weight <= M->best_match_weight ? */
            u_depth_1 = M->match[e->numbout].neigh;
            if (e->numbout != w_depth_1 && u_depth_1 != w_depth_1) {
                if ((help = (e->weight - M->match[u_depth_1].weight +
                        (CCutil_dat_edgelen (D->v, u_depth_1, D->dat))))
                        < min) {
                    /* Einstellung */
                    test++;
                    if (qu->in[u_depth_1] != qu->in_Magic) {
                        test = 0;
                        min = help;
                        minw = e->numbout;
                    }
                }
            }
        }
        if (test > 0) return -1;
    } else if (D->depth > 4) { /* Einstellung */
        for (e = D->adj[u_depth]; e != D->end &&
             e->weight + d < M->best_match_weight; e=e->next) {
             /* is new p_weight <= M->best_match_weight ? */
            if (qu->in[u_depth_1 = M->match[e->numbout].neigh] !=
                qu->in_Magic ) {
                if ((help = (e->weight - M->match[u_depth_1].weight +
                    (CCutil_dat_edgelen(D->v,u_depth_1,D->dat) >> 1))) < min) {
                    /* Einstellung */
                    min = help;
                    minw = e->numbout;
                }
            }
        }
    } else {
        for (e = D->adj[u_depth]; e != D->end &&
            (e->weight + d < M->best_match_weight); e=e->next) {
            /* is new p_weight <= M->best_match_weight ? */
            if (qu->in[u_depth_1 = M->match[e->numbout].neigh] !=
                qu->in_Magic) {
                if ((help = (e->weight - M->match[u_depth_1].weight)) < min)  {
                    /* Einstellung */
                    min = help;
                    minw = e->numbout;
                }
            }
        }
    }
    return minw;
}

static void edge_scan2 (data *D, matchings *M, tflip *fli, mqueue *qu)
{
    double p_wei;

    /* step 3 : test match nicht noetig, wird in flipp gemacht ! */
    /* step 4 */

    D->w[D->depth] = find_w (D, M, qu, D->u[D->depth], D->w[D->depth-1]);
    p_wei = M->p_weight;
    if (D->w[D->depth] != -1) {
        D->u[D->depth + 1] = M->match[D->w[D->depth]].neigh;
        M->p_weight -= CCutil_dat_edgelen(D->w[D->depth-1], D->u[D->depth],
                                          D->dat);
        flipp (D,M, fli, D->u[D->depth], D->v,D->w[D->depth],
               D->u[D->depth+1]);
        M->p_weight += M->match[D->w[D->depth]].weight;
        qu->in[D->w[D->depth]]=qu->in_Magic;
        D->depth++;
        qu->in[D->u[D->depth]]=qu->in_Magic;

        edge_scan2(D,M,fli,qu);

        D->depth--;
        unflipp (M, fli, D->u[D->depth], D->w[D->depth], D->v,
                 D->u[D->depth+1],
                 fli->flip_stack[fli->flip_stack_count-1].d1,
                 fli->flip_stack[fli->flip_stack_count-1].d2);
        M->p_weight=p_wei;
    }
}

/* ************************************************************/
/*                         Step345                            */
/*  Vorher muessen u0 bestimmt sein !                         */
/* ************************************************************/

static void step345 (data *D, matchings *M, tflip *fli, mqueue *qu)
{
    double d0,d1;
    int z0,z1;
    double p_wei;


    fli->flip_stack_count = 0;

    qu->in_Magic++;
    d0 = M->p_weight - CCutil_dat_edgelen (D->u[0], D->v, D->dat);
    find_w0 (D, M, qu, D->u[0], D->v);
    for (z0 = 0; z0 < 4 && ((D->w[0] = D->w0[z0].w)!=-1) &&
        ((D->w0[z0].weight - D->w0[0].weight)*D->nnodes < M->match_weight);
        z0++) {
      if (d0 + CCutil_dat_edgelen (D->u[0], D->w[0], D->dat) <
          M->best_match_weight) {
            D->u[1] = M->match[D->w[0]].neigh;
            flipp (D, M, fli, D->u[0], D->v, D->w[0], D->u[1]);
            M->p_weight = d0 + M->match[D->u[0]].weight;
            p_wei = M->p_weight;
            qu->in_Magic++;
            qu->in[D->u[0]] = qu->in_Magic;
            qu->in[D->u[1]] = qu->in_Magic;
            qu->in[D->w[0]] = qu->in_Magic;
            qu->in[D->v] = qu->in_Magic;
            d1 = M->p_weight - CCutil_dat_edgelen (D->u[1], D->w[0], D->dat);
            find_w1 (D, M, qu, D->u[1]);
            for (z1 = 0; z1 < 2 && ((D->w[1] = D->w1[z1].w) != -1) &&
                ((D->w1[z1].weight - D->w1[0].weight) * 2 * D->nnodes <
                  M->match_weight); z1++) {
                if (d1 + CCutil_dat_edgelen (D->u[1], D->w[1], D->dat) <
                    M->best_match_weight) {
                    D->u[2]=M->match[D->w[1]].neigh;
                    qu->in_Magic++;
                    qu->in[D->u[0]]=qu->in_Magic;
                    qu->in[D->u[1]]=qu->in_Magic;
                    qu->in[D->w[0]]=qu->in_Magic;
                    qu->in[D->v]=qu->in_Magic;
                    qu->in[D->w[1]]=qu->in_Magic;
                    qu->in[D->u[2]]=qu->in_Magic;
                    flipp(D,M,fli,D->u[1],D->v,D->w[1],D->u[2]);
                    M->p_weight=d1+M->match[D->u[1]].weight;
                    D->depth=2;
                    /* w2 to w ..... */
                    edge_scan2(D,M,fli,qu);D->depth--;
                    /* step 5:  undo all flipps and M->p_weights exept first */
                    for (D->depth=D->depth;D->depth>0;D->depth--) {
                        unflipp (M, fli, D->u[D->depth], D->w[D->depth], D->v,
                                 D->u[D->depth+1],
                                 fli->flip_stack[fli->flip_stack_count-1].d1,
                                 fli->flip_stack[fli->flip_stack_count-1].d2);
                    }
                    M->p_weight=p_wei;
                }
            }
            /* step 5:  undo first flip  */
            unflipp (M, fli, D->u[0], D->w[0], D->v, D->u[1],
                    fli->flip_stack[fli->flip_stack_count-1].d1,
                    fli->flip_stack[fli->flip_stack_count-1].d2);
            M->p_weight = M->best_match_weight;
        }
    }
}

static int edge_scan1 (data *D, matchings *M, tflip *fli, mqueue *qu)
{
    M->p_weight = M->best_match_weight;
    M->match_weight = M->best_match_weight;
    D->new_best = 0;

    /* step 2, the first neighbour */

    D->u[0] = M->match[D->v].neigh;
    step345 (D, M, fli, qu);

    return(D->new_best);
}

static void lin_ker_do (data *D, matchings *M, tflip *fli, mqueue *qu)
{
    int i;

    fli->flip_all_stack_count = 0;

    /* step 1 */

    while (non_empty_q (qu)) {
        D->v = get_q (qu);
        qu->ok[D->v] = 1;
        if (edge_scan1 (D, M, fli, qu) == 1) {  /* best_match changed ? */
            for (i = 0; i < fli->flip_best_stack_count; i++) {
                if (D->iteration > 0)
                     fli->flip_all_stack[fli->flip_all_stack_count+i] =
                          fli->flip_best_stack[i];
                M->match[fli->flip_best_stack[i].a].neigh =
                     fli->flip_best_stack[i].c;
                M->match[fli->flip_best_stack[i].b].neigh =
                     fli->flip_best_stack[i].d;
                M->match[fli->flip_best_stack[i].c].neigh =
                     fli->flip_best_stack[i].a;
                M->match[fli->flip_best_stack[i].d].neigh =
                     fli->flip_best_stack[i].b;
                M->match[fli->flip_best_stack[i].a].weight =
                     CCutil_dat_edgelen (fli->flip_best_stack[i].c,
                                  fli->flip_best_stack[i].a, D->dat);
                M->match[fli->flip_best_stack[i].b].weight =
                    CCutil_dat_edgelen(fli->flip_best_stack[i].d,
                                fli->flip_best_stack[i].b, D->dat);
                M->match[fli->flip_best_stack[i].c].weight =
                    M->match[fli->flip_best_stack[i].a].weight;
                M->match[fli->flip_best_stack[i].d].weight =
                    M->match[fli->flip_best_stack[i].b].weight;
                if (qu->ok[fli->flip_best_stack[i].a] == 1)  {
                    put_q (D, qu, fli->flip_best_stack[i].a);
                    qu->ok[fli->flip_best_stack[i].a] = 0;
                }
                if (qu->ok[fli->flip_best_stack[i].b] == 1)  {
                    put_q (D, qu, fli->flip_best_stack[i].b);
                    qu->ok[fli->flip_best_stack[i].b] = 0;
                }
                if (qu->ok[fli->flip_best_stack[i].c]==1)  {
                    put_q (D, qu, fli->flip_best_stack[i].c);
                    qu->ok[fli->flip_best_stack[i].c] = 0;
                }
                if (qu->ok[fli->flip_best_stack[i].d]==1)  {
                    put_q (D, qu, fli->flip_best_stack[i].d);
                    qu->ok[fli->flip_best_stack[i].d] = 0;
                }
            }
            if (D->iteration > 0)
                fli->flip_all_stack_count += fli->flip_best_stack_count;
        }
        M->match_weight = M->best_match_weight;
    }
}

static int m_lin_ker (CCkdtree *kt, CCdatagroup *dat, data *D, matchings *M,
        graph *G, tflip *fli,mqueue *qu, int max_iteration,
        CCrandstate *rstate)
{
    int i, ende;

    D->print_level = 0;
    if (D->print_level) {
        printf("\n Initialization  ..");
        fflush (stdout);
    }
    if (init_graph (D, G)) {
        fprintf (stderr, "init_graph failed\n");
        return 1;
    }
    if (init (D, fli, qu)) {
        fprintf (stderr, "init failed\n");
        return 1;
    }
    if (D->print_level) {
        printf(".. ready !!! ");
        fflush(stdout);
    }

    M->match = nn_match (kt, dat, D, M, rstate);
    if (!M->match) {
        fprintf (stderr, "nn_match failed\n");
        frei (D, fli, qu);
        return 1;
    }
    memcpy ((char *) M->very_best_match, (char *) M->match,
            D->nnodes*sizeof(node));

    M->match_weight = 0;
    for (i = 0; i < D->nnodes; i++)
        M->match_weight += CCutil_dat_edgelen (i, M->match[i].neigh, D->dat);
    M->match_weight/=2;
    M->best_match_weight = M->match_weight;
    M->very_best_match_weight = M->match_weight;

    if (opt_makegraph (kt, dat, D, rstate)) {
        fprintf (stderr, "opt_makegraph failed\n");
        frei (D, fli, qu);
        return 1;
    }
    init_queue (D, qu);

    if (D->print_level) {
        printf("\n  Start LinKer with %i iterations ... ", max_iteration);
        fflush(stdout);
    }

    ende = 0;
    for (D->iteration = 0; ende == 0; D->iteration++) {
        if (D->print_level)
            printf("\n    %i. It., ..",D->iteration);
        else
            printf(".");
        fflush(stdout);
        lin_ker_do (D, M, fli, qu);
        if (D->print_level)
            printf(" | Length:%f |", M->best_match_weight);
        if (restore_match (D, M, G, fli)) {
            fprintf (stderr, "restore_match failed\n");
            frei (D, fli, qu);
            return 1;
        }
        if (D->iteration < max_iteration) { /* This is a big kick */
            int n1, n2, m1 = 0, m2, c1, c2;
            optedge *e;
            for (i = 0; i <= D->nnodes/5; i++) {
                n1 = CCutil_lprand(rstate) % D->nnodes;
                n2 = M->match[n1].neigh;
                c2 = CCutil_lprand (rstate) % 8;
                e = D->adj[n1];
                for (c1 = 0; c1 <= c2;) {
                    m1 = e->numbout;
                    e = e->next;
                    if (e == D->end) {
                        if (c1 == 0) {
                            n1 = CCutil_lprand (rstate) % D->nnodes;
                            n2 = M->match[n1].neigh;
                        }
                        e = D->adj[n1];
                    }
                    if (m1 != n1 && m1 != n2) c1++;
                }
                m2 = M->match[m1].neigh;
                flip (D, M->match, n1, n2, m1, m2);
            }
            for (i = 0; i < D->nnodes; i++) {
                qu->in[i] =- 1;
                qu->ok[i] = 0;
                put_q (D, qu, i);
            }
            M->match_weight = 0;
            for (i = 0; i < D->nnodes; i++)
                M->match_weight += CCutil_dat_edgelen (i, M->match[i].neigh,
                                                       D->dat);
            M->match_weight/=2;
            M->best_match_weight = M->match_weight;
        } else {
            ende=1;
        }
        if (D->print_level) {
            printf(" %i edg",G->edge_count);
            fflush (stdout);
        }
    }

    frei (D, fli, qu);
    return 0;
}

int CCedgegen_mlinkern (int ncount, CCdatagroup *dat, int wantlist,
        int *ecount, int **elist, CCkdtree *kt, int iterations,
        CCrandstate *rstate)
{
    data D;
    matchings M;
    graph G;
    mqueue qu;
    tflip fli;
    int i,j;
    edge *e;
    int rval = 0;
    int norm;

    D.nnodes = ncount;
    D.dat = dat;

    if (wantlist) {
        *ecount = 0;
        *elist = (int *) NULL;
    }

    CCutil_dat_getnorm (dat, &norm);
    if ((norm & CC_NORM_BITS) != CC_KD_NORM_TYPE) {
        printf ("Cannot run matching Lin-Kernighan with this norm\n");
        fflush (stdout);
        return 0;
    }

    /* Start linker */

    if (m_lin_ker (kt, dat, &D, &M, &G, &fli, &qu, iterations, rstate)) {
        fprintf (stderr, "m_lin_ker failed\n");
        rval = 1;
        goto CLEANUP;
    }

    /* Put edges in the correct format */

    *ecount = G.edge_count;
    *elist = CC_SAFE_MALLOC (2 *(*ecount), int);
    if (!(*elist)) {
        fprintf (stderr, "out of memory in mlinkern\n");
        rval = 1;
        goto CLEANUP;
    }

    for (j = 0, i = 0; i < D.nnodes; i++) {
        for (e = G.adj_graph[i]; e != G.end_graph; e = e->next) {
             (*elist)[j++] = i;
             (*elist)[j++] = e->numbout;
        }
    }

CLEANUP:

    CC_IFFREE (G.edge_graph, edge);
    CC_IFFREE (G.adj_graph, edge*);
    CC_IFFREE (G.end_graph, edge);
    CC_IFFREE (M.match, node);
    CC_IFFREE (M.very_best_match, node);

    return rval;
}

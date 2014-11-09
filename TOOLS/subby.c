#include <stdio.h>
#include "machdefs.h"
#include "macrorus.h"
#include "util.h"

typedef struct subby {
    double pi;
    int count;
    int *nodes;
} subby;

static int subby_cross (subby *s, subby *t, int *hit);
static int subby_equal (subby *s, subby *t, int *hit);
static int subby_uncross (int is, int it, int *scount, subby **slist, int *hit,
       int ncount, double *pi);
static int subby_complement (subby *s, int ncount);
static int print_the_subs (int ncount, double *pi, int *p_scount,
        subby **slist, subby *ssupply, CCdatagroup *dat);


int main (int ac, char **av)
{
    CCdatagroup dat;
    FILE *in = (FILE *) NULL;
    int i, norm, ncount, smax, scount = 0, k;
    subby **slist = (subby **) NULL;
    subby *ssupply = (subby *) NULL;
    int seed = 99;
    CCrandstate rstate;
    int rval = 0;
    double *pi = (double *) NULL;

    CCutil_init_datagroup (&dat);
    seed = (int) CCutil_real_zeit ();

    if (ac != 4) {
        fprintf (stderr, "Usage: %s tsp pi cuts\n", av[0]);
        rval = 1; goto CLEANUP;
    }

    CCutil_printlabel ();
    CCutil_sprand (seed, &rstate);

    rval = CCutil_gettsplib (av[1], &ncount, &dat);
    CCcheck_rval (rval, "CCutil_gettsplib failed");                     
    CCutil_dat_getnorm (&dat, &norm);

    if ((norm & CC_NORM_SIZE_BITS) != CC_D2_NORM_SIZE) {                    
        fprintf (stderr, "Only set up for 2D norms\n");
        rval = 1;  goto CLEANUP;                                            
    } 

    in = fopen (av[2], "r");
    if (!in) {
        fprintf (stderr, "could not open %s for reading\n", av[2]);
        rval = 1;  goto CLEANUP;
    }

    printf ("ncount = %d\n", ncount); fflush (stdout);


    fscanf (in, "%d", &k);
    if (k != ncount) {
        fprintf (stderr, "pi file does not match dat file\n");
        rval = 1;  goto CLEANUP;
    }
    pi = CC_SAFE_MALLOC (ncount, double);
    CCcheck_NULL (pi, "out of memory for pi");

    for (i = 0; i < ncount; i++) {
        fscanf (in, "%lf", &pi[i]);
    }
    fclose (in);

    in = fopen (av[3], "r");
    if (!in) {
        fprintf (stderr, "could not open %s for reading\n", av[2]);
        rval = 1;  goto CLEANUP;
    }

    fscanf (in, "%d %d", &k, &scount);
    if (k != ncount) {
        fprintf (stderr, "cut file does not match pi file\n");
        rval = 1;  goto CLEANUP;
    }
    smax = 100000 + scount + 2*ncount;

    printf ("%d Cuts\n", scount); fflush (stdout);

    slist = CC_SAFE_MALLOC (smax, subby *);
    CCcheck_NULL (slist, "out of memory for slist");
    ssupply = CC_SAFE_MALLOC (smax, subby);
    CCcheck_NULL (ssupply, "out of memory for ssupply");
    for (i = 0; i < smax; i++) slist[i] = &ssupply[i];
    
    for (k = 0; k < scount; k++) {
        fscanf (in, "%lf %d", &(slist[k]->pi), &(slist[k]->count));
        slist[k]->nodes = CC_SAFE_MALLOC (slist[k]->count, int);
        for (i = 0; i < slist[k]->count; i++) {
            fscanf (in, "%d", &(slist[k]->nodes[i]));
        }
    }
    fclose (in);

    rval = print_the_subs (ncount, pi,  &scount, slist, ssupply, &dat);
    CCcheck_rval (rval, "print_the_subs failed");


CLEANUP:

    CCutil_freedatagroup (&dat);
    if (in) fclose (in);

    for (i = 0; i < scount; i++) {
        CC_IFFREE (slist[i]->nodes, int);
    }
    CC_IFFREE (slist, subby *);
    CC_IFFREE (ssupply, subby);

    return rval;
}

static int print_the_subs (int ncount, double *pi, int *p_scount,
        subby **slist, subby *ssupply, CCdatagroup *dat)
{
    int rval = 0;
    int i, j, k, icnt;
    int *hit = (int *) NULL;
    int *sperm = (int *) NULL;
    int *sval = (int *) NULL;
    FILE *out = (FILE *) NULL;
    char buf[1024];
    double len, delta;
    int n0, n1, ecnt = 0;
    int scount = *p_scount;
    int *ar = (int *) NULL;

    printf ("print_the_subs ...\n"); fflush (stdout);


    hit = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (hit, "out of memory in print_the_subs");

    for (i = 0; i < ncount; i++)  {
        hit[i] = 0;
    }

    for (n0 = 0; n0 < ncount; n0++) {
        for (n1 = n0+1; n1 < ncount; n1++) {
            len = (double) CCutil_dat_edgelen (n0, n1, dat);
            if (pi[n0] + pi[n1] > len + 1.0) {
                delta = (pi[n0] + pi[n1] - len) / 2.0;
                if (pi[n0] < delta) delta = pi[n0];
                if (pi[n1] < delta) delta = pi[n1];
               
                ar = CC_SAFE_MALLOC (2, int);
                CCcheck_NULL (ar, "out of memory in print_the_subs");
                ar[0] = n0;
                ar[1] = n1;
                slist[scount]->pi = delta;
                slist[scount]->count = 2;
                slist[scount++]->nodes = ar;
                pi[n0] -= delta;
                pi[n1] -= delta;
                ecnt++;
            }
        }
    }
    printf ("%d 2-node subtours\n", ecnt); fflush (stdout);

    icnt = 0;

DOGGY:

    printf ("Check for crossings\n"); fflush (stdout);

    for (i = 0; i < scount; i++) {
        for (j = i+1; j < scount; j++) {
            if (subby_cross (slist[i], slist[j], hit)) {
                printf ("Cross X[%d,%d]\n", i, j);
                rval = subby_uncross (i, j, &scount, slist, hit, ncount, pi);
                CCcheck_rval (rval, "subby_uncross failed");
                icnt++;
                goto DOGGY;
            }
        }
    }
    if (icnt) printf ("\n");
    printf ("Number of intersections: %d\n", icnt); fflush (stdout);

    sperm = CC_SAFE_MALLOC (scount, int);
    CCcheck_NULL (sperm, "out of memory in print_the_subs");
    sval = CC_SAFE_MALLOC (scount, int);
    CCcheck_NULL (sval, "out of memory in print_the_subs");
    for (i = 0; i < scount; i++) {
        sperm[i] = i;
        sval[i] = slist[i]->count;
    }
    CCutil_int_perm_quicksort (sperm, sval, scount);

    sprintf (buf, "full.bnd");
    out = fopen (buf, "w");
    if (!out) {
        fprintf (stderr, "could not open %s for writing\n", buf);
        rval = 1;  goto CLEANUP;
    }
    fprintf (out, "%d\n", ncount);
    for (i = 0; i < ncount; i++) {
        fprintf (out, "%f\n", pi[i]);
    }
    printf ("scount = %d\n", scount); fflush (stdout);
    fprintf (out, "%d\n", scount);
    for (i = 0; i < scount; i++) {
        k = sperm[i];
        fprintf (out, "%f %d ", slist[k]->pi, slist[k]->count);
        for (j = 0; j < slist[k]->count; j++) {
            fprintf (out, "%d ", slist[k]->nodes[j]);
            if (j % 15 == 14) fprintf (out, "\n");
        }
        fprintf (out, "\n"); 
    }

    fclose (out);
    
    *p_scount = scount;

CLEANUP:

    return rval;
}

static int subby_complement (subby *s, int ncount)
{
    int i, k = 0,  rval = 0;
    char *hit = (char *) NULL;

    hit = CC_SAFE_MALLOC (ncount, char);
    CCcheck_NULL (hit, "out of memory in subby_complement");

    for (i = 0; i < ncount; i++) hit[i] = 0;
    for (i = 0; i < s->count; i++) hit[s->nodes[i]] = 1;
  
    CC_IFFREE (s->nodes, int);
    s->nodes = CC_SAFE_MALLOC (ncount - s->count, int);
    CCcheck_NULL (s->nodes, "out of memory in subby_complement");

    for (i = 0; i < ncount; i++) {
        if (!hit[i]) {
            s->nodes[k++] = i;
        }
    }
    if (k != ncount - s->count) {
        fprintf (stderr, "lost a node\n");
        rval = 1;  goto CLEANUP;
    }
    s->count = k;

CLEANUP:

    CC_IFFREE (hit, char);
    return rval;
}

static int subby_cross (subby *s, subby *t, int *hit)
{
    int i, in = 0, out = 0;
    subby *tmp;

    if (s->count < t->count) {
        CC_SWAP (s, t, tmp);
    }

    for (i = 0; i < s->count; i++) hit[s->nodes[i]] = 1;
    for (i = 0; i < t->count; i++) {
        if (hit[t->nodes[i]]) {
            in = 1;
        } else {
            out = 1;
        }
    }
    for (i = 0; i < s->count; i++) hit[s->nodes[i]] = 0;

    if (in && out) return 1;
    else           return 0;
}

static int subby_equal (subby *s, subby *t, int *hit)
{
    int i, val = 1;

    if (s->count != t->count) return 0;

    for (i = 0; i < s->count; i++) hit[s->nodes[i]] = 1;
    for (i = 0; i < t->count; i++) {
        if (!hit[t->nodes[i]]) {
            val = 0;  break;
        }
    }
    for (i = 0; i < s->count; i++) hit[s->nodes[i]] = 0;

    return val;
}


static int subby_uncross (int is, int it, int *scount, subby **slist, int *hit,
       int ncount, double *pi)
{
    int rval = 0, i, k, iu, iv, no_u = 0, no_v = 0;
    subby *s = slist[is];
    subby *t = slist[it];
    subby *tmp, *u, *v;
    int icount = 0, ucount = 0;

    printf ("uncross subtour %d (%d count) and subtour %d (%d count)\n",
              is, s->count, it, t->count);
    fflush (stdout);

    if (s->pi < t->pi) {
        CC_SWAP (s, t, tmp);
        CC_SWAP (is, it, i);
    }

    s->pi -= t->pi;
    
    for (i = 0; i < s->count; i++) hit[s->nodes[i]] = 1;
    ucount = s->count;
    for (i = 0; i < t->count; i++) {
        if (hit[t->nodes[i]]) icount++;
        else                  ucount++;
    }

    if (ucount == ncount) {
       fprintf (stderr, "union of subtours is entire node set\n");
       rval = 1;  goto CLEANUP;
    }

    iu = *scount;
    iv = *scount + 1;
    u = slist[iu];
    v = slist[iv];

    u->pi = t->pi;
    u->count = ucount;
    u->nodes = CC_SAFE_MALLOC (u->count, int);
    CCcheck_NULL (u->nodes, "out of memory in subby_uncross");

    v->pi = t->pi;
    v->count = icount;
    v->nodes = CC_SAFE_MALLOC (v->count, int);
    CCcheck_NULL (v->nodes, "out of memory in subby_uncross");


    icount = ucount = 0;
    for (i = 0; i < s->count; i++) {
        u->nodes[ucount++] = s->nodes[i];
    }

    for (i = 0; i < t->count; i++) {
        if (hit[t->nodes[i]]) {
            v->nodes[icount++] = t->nodes[i];
        } else {
            u->nodes[ucount++] = t->nodes[i];
        }
    }
    if (ucount != u->count || icount != v->count) {
        fprintf (stderr, "bad counts in uncross\n");
        rval = 1;  goto CLEANUP;
    }

    for (i = 0; i < s->count; i++) hit[s->nodes[i]] = 0;

    for (i = 0; i < *scount; i++) {
        if (subby_equal (slist[i], u, hit)) {
            slist[i]->pi += u->pi;
            no_u = 1;
            break;
        }
    }

    if (v->count == 1) {
        k = v->nodes[0];
        pi[k] += slist[i]->pi;
        no_v = 1;
    } else {
        for (i = 0; i < *scount; i++) {
            if (subby_equal (slist[i], v, hit)) {
                slist[i]->pi += v->pi;
                no_v = 1;
                break;
            }
        }
    }

    CC_IFFREE (t->nodes, int);
    if (no_u) {
       CC_IFFREE (u->nodes, int);
    }
    if (no_v) {
       CC_IFFREE (v->nodes, int);
    }

    if (no_u && no_v) {
        CC_SWAP (slist[it], slist[*scount-1], tmp);
        (*scount)--;
    } else if (no_u) {
        CC_SWAP (slist[it], slist[iv], tmp);
    } else if (no_v) {
        CC_SWAP (slist[it], slist[iu], tmp);
    } else {
        CC_SWAP (slist[it], slist[iv], tmp);
        (*scount)++;
    }

CLEANUP:

    for (i = 0; i < ncount; i++) hit[i] = 0;

    return rval;
}


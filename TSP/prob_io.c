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
/* Date: April 19, 1996  (dave)                                             */
/*       August 28, 1996 (bico)                                             */
/*       September 30, 1997 (dave)                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  CCtsp_PROB_FILE *CCtsp_prob_read (char *f, int n)                       */
/*    NONE                                                                  */
/*                                                                          */
/*  CCtsp_PROB_FILE *CCtsp_prob_read_name (char *f)                         */
/*    NONE                                                                  */
/*                                                                          */
/*  CCtsp_PROB_FILE *CCtsp_prob_write (char *f, int n)                      */
/*    NONE                                                                  */
/*                                                                          */
/*  CCtsp_PROB_FILE *CCtsp_prob_write_name (char *fname)                    */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_file_delete (char *f, int n)                             */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_getname (CCtsp_PROB_FILE *p, char *name)                 */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_getid (CCtsp_PROB_FILE *p, int *id)                      */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_getparent (CCtsp_PROB_FILE *p, int *parent)              */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_getub (CCtsp_PROB_FILE *p, double *ub)                   */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_getlb (CCtsp_PROB_FILE *p, double *lb)                   */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_getexactlb (CCtsp_PROB_FILE *p, CCbigguy *lb)            */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_getnnodes (CCtsp_PROB_FILE *p, int *nnodes)              */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_getchildren (CCtsp_PROB_FILE *p, int *child0,            */
/*      int *child1)                                                        */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_getreal (CCtsp_PROB_FILE *p, int *real)                  */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_getprocessed (CCtsp_PROB_FILE *p, int *processed)        */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_getinfeasible (CCtsp_PROB_FILE *p, int *infeasible)      */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_gettour (CCtsp_PROB_FILE *p, int ncount, int **tour,     */
/*      int silent)                                                         */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_getedges (CCtsp_PROB_FILE *p, int ncount, int *nedges,   */
/*      int **elist, int **elen, int silent)                                */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_getcuts (CCtsp_PROB_FILE *p, int *ncount,                */
/*      CCtsp_lpcuts *cuts, int silent)                                     */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_getwarmstart (CCtsp_PROB_FILE *p, CClp_warmstart **w,    */
/*      int silent)                                                         */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_getfulladj (CCtsp_PROB_FILE *p, int ncount,              */
/*      int *fullcount, CCtsp_genadj **adj,                                 */
/*      CCtsp_genadjobj **adjspace, int silent)                             */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_getfixed (CCtsp_PROB_FILE *p, int ncount, int *ecount,   */
/*      int **elist, int silent)                                            */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_getexactdual (CCtsp_PROB_FILE *p, ncount,                */
/*      CCtsp_bigdual **d, int silent)                                      */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_gethistory (CCtsp_PROB_FILE *p, int *depth,              */
/*      CCtsp_branchobj **history, int silent)                              */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_rclose (CCtsp_PROB_FILE *p)                              */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_putname (CCtsp_PROB_FILE *p, char *name)                 */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_putid (CCtsp_PROB_FILE *p, int id)                       */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_putparent (CCtsp_PROB_FILE *p, int parent)               */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_putub (CCtsp_PROB_FILE *p, double ub)                    */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_putlb (CCtsp_PROB_FILE *p, double lb)                    */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_putexactlb (CCtsp_PROB_FILE *p, CCbigguy lb)             */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_putnnodes (CCtsp_PROB_FILE *p, int nnodes)               */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_putchildren (CCtsp_PROB_FILE *p, int child0,             */
/*      int child1)                                                         */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_putreal (CCtsp_PROB_FILE *p, int real)                   */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_putprocessed (CCtsp_PROB_FILE *p, int processed)         */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_putinfeasible (CCtsp_PROB_FILE *p, int infeasible)       */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_puttour (CCtsp_PROB_FILE *p, int ncount, int *tour)      */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_putedges (CCtsp_PROB_FILE *p, int ncount, int nedges,    */
/*      int *elist, int *elen)                                              */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_putcuts (CCtsp_PROB_FILE *p, int ncount,                 */
/*      CCtsp_lpcuts *cuts)                                                 */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_putwarmstart (CCtsp_PROB_FILE *p, CClp_warmstart *w)     */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_putfulladj (CCtsp_PROB_FILE *p, int ncount,              */
/*      int fullcount, CCtsp_genadj *adj)                                   */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_prob_putfixed (CCtsp_PROB_FILE *p, int ncount,                */
/*      int ecount, int *elist)                                             */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_prob_putexact_dual (CCtsp_PROB_FILE *p,                       */
/*      CCtsp_bigdual *exact_dual, int ncount)                              */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_puthistory (CCtsp_PROB_FILE *p, int depth,               */
/*      CCtsp_branchobj *history)                                           */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_wclose (CCtsp_PROB_FILE *p)                              */
/*    NONE                                                                  */
/*                                                                          */
/*  int CCtsp_prob_copy_section (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t,    */
/*      char section, int silent)                                           */
/*    NONE                                                                  */
/*                                                                          */
/*  char *CCtsp_problabel (const char *probloc)                             */
/*    -RETURNS a copy of the probname portion of probfile, or NULL if       */
/*     unable to allocate space for the copy.                               */
/*                                                                          */
/*  CCtsp_PROB_FILE *CCtsp_prob_read_remote (char *hname, char *pname,      */
/*      int n)                                                              */
/*    Only exists if CC_NETREADY is defined                                 */
/*                                                                          */
/*  CCtsp_PROB_FILE *CCtsp_prob_write_remote (char *hname, char *pname,     */
/*      int n)                                                              */
/*    Only exists if CC_NETREADY is defined                                 */
/*                                                                          */
/*  CCtsp_PROB_FILE *CCtsp_prob_server (CC_SFILE *s)                        */
/*    Only exists if CC_NETREADY is defined                                 */
/*                                                                          */
/*  int CCtsp_prob_delete_remote (char *hname, char *pname, int n)          */
/*    Only exists if CC_NETREADY is defined                                 */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"
#include "bigguy.h"

#define PROB_LOCAL  1
#define PROB_REMOTE 2
#define PROB_SERVER 3

#define PROB_HEADER_VERSION    (1)
#define PROB_TOUR_VERSION      (1)
#define PROB_EDGES_VERSION     (1)
#define PROB_CUTS_VERSION      (1)
#define PROB_FULLADJ_VERSION   (1)
#define PROB_FIXED_VERSION     (1)
#define PROB_EXACTDUAL_VERSION (1)
#define PROB_HISTORY_VERSION   (1)


static int
    prob_getheader (CCtsp_PROB_FILE *p, CCtsp_PROB_FILE *h),
    prob_putheader (CCtsp_PROB_FILE *p, CCtsp_PROB_FILE *h),
    prob_name (char *buf, size_t buflen, char *f, int n),
    prob_copyheader (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t),
    begin_put (CCtsp_PROB_FILE *p, int *offset, char section),
    begin_get (CCtsp_PROB_FILE *p, int offset, char section, int silent),
    begin_copy (CCtsp_PROB_FILE *f, int foffset, CCtsp_PROB_FILE *t,
        int *toffset, char section, int silent),
    prob_copytour (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t, int silent),
    prob_copyedges (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t, int silent),
    prob_copycuts (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t, int silent),
    prob_copywarmstart (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t, int silent),
    prob_copyfulladj (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t, int silent),
    prob_copyfixed (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t, int silent),
    prob_copyexactdual (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t, int silent),
    prob_copyhistory (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t, int silent),
    remote_name (const char *f),
    split_name (const char *f, char *hostname, size_t hlen, char *probname, 
        size_t plen),
    copy_ints (CC_SFILE *f, CC_SFILE *t, int n),
    copy_int_val (CC_SFILE *f, CC_SFILE *t, int *n),
    copy_bits (CC_SFILE *f, CC_SFILE *t, int n, int nbits),
    copy_bits_val (CC_SFILE *f, CC_SFILE *t, int *n, int nbits),
    copy_chars (CC_SFILE *f, CC_SFILE *t, int n),
    copy_char_val (CC_SFILE *f, CC_SFILE *t, char *c),
    copy_bigguys (CC_SFILE *f, CC_SFILE *t, int n);

static void
    prob_init (CCtsp_PROB_FILE *p);



int CCtsp_prob_file_delete (char *f, int n)
{
    char nambuf[1024];
    int sval;

    if (remote_name (f)) {
#ifdef CC_NETREADY
        char hostbuf[1024];
        if (split_name (f, hostbuf, sizeof (hostbuf),
                        nambuf, sizeof (nambuf))) {
            fprintf (stderr, "Cannot split remote name\n");
            return -1;
        }
        return CCtsp_prob_delete_remote (hostbuf, nambuf, n);
#else /* CC_NETREADY */
        fprintf (stderr, "Remote problem deleting not enabled\n");
        return -1;
#endif
    }
    
    if (prob_name (nambuf, sizeof (nambuf), f, n)) return 1;
    /*
    printf ("Delete File %s at time %.0f\n", nambuf, CCutil_real_zeit());
    fflush (stdout);
    */
    sval = CCutil_sdelete_file (nambuf);
    if (sval) {
        fprintf (stderr, "Prob file %s could not be deleted\n", nambuf);
    }
    sval = CCutil_sdelete_file_backup (nambuf);
    /*
    if (!sval) {
        printf ("Deleted backup to file: %s\n", nambuf);
        fflush (stdout);
    }
    */

    return 0;
}

CCtsp_PROB_FILE *CCtsp_prob_read (char *f, int n)
{
    char nambuf[1024];

    if (remote_name (f)) {
#ifdef CC_NETREADY
        char hostbuf[1024];
        if (split_name (f, hostbuf, sizeof (hostbuf),
                        nambuf, sizeof (nambuf))) {
            fprintf (stderr, "Cannot split remote name\n");
            return (CCtsp_PROB_FILE *) NULL;
        }
        return CCtsp_prob_read_remote (hostbuf, nambuf, n);
#else /* CC_NETREADY */
        fprintf (stderr, "Remote problem reading not enabled\n");
        return (CCtsp_PROB_FILE *) NULL;
#endif
    }

    if (prob_name (nambuf, sizeof (nambuf), f, n))
        return (CCtsp_PROB_FILE *) NULL;

    return CCtsp_prob_read_name (nambuf);
}

CCtsp_PROB_FILE *CCtsp_prob_read_name (char *f)
{
    CCtsp_PROB_FILE *p = (CCtsp_PROB_FILE *) NULL;
    int rval;

#if 0
    if (remote_name (f)) {
        fprintf (stderr, "Cannot read remote problems by name\n");
        return (CCtsp_PROB_FILE *) NULL;
    }
#endif
    
    /*
    printf ("Read File %s at time %.0f\n", f, CCutil_real_zeit());
    fflush (stdout);
    */

    p = CC_SAFE_MALLOC (1, CCtsp_PROB_FILE);
    if (p == (CCtsp_PROB_FILE *) NULL) goto FAILURE;
    prob_init (p);

    p->f = CCutil_sopen (f, "r");
    if (!p->f) goto FAILURE;

    p->type = PROB_LOCAL;

    rval = prob_getheader (p, p);
    if (rval) {
        fprintf (stderr, "prob_getheader failed\n");
        goto FAILURE;
    }

    return p;
    
  FAILURE:
    if (p) {
        if (p->f) {
            CCutil_sclose (p->f);
        }
        CC_FREE (p, CCtsp_PROB_FILE);
    }
    return (CCtsp_PROB_FILE *) NULL;
}

#ifdef CC_NETREADY

CCtsp_PROB_FILE *CCtsp_prob_read_remote (char *hname, char *pname, int n)
{
    CCtsp_PROB_FILE *p = (CCtsp_PROB_FILE *) NULL;
    int rval;

    /* Hack for run at Princeton without .caam.rice.edu */

/*
    START HACK
    char tname[4048];
    sprintf (tname, "%s.caam.rice.edu", hname);
    END HACK - also change hname to tname
*/

    /*
    printf ("Read Remote Host %s name %s id %d\n", hname, pname, n);
    fflush (stdout);
    */

    p = CC_SAFE_MALLOC (1, CCtsp_PROB_FILE);
    if (p == (CCtsp_PROB_FILE *) NULL) goto FAILURE;
    prob_init (p);

    p->f = CCutil_snet_open (hname, CCtsp_PROB_PORT);
    if (p->f == (CC_SFILE *) NULL) {
        fprintf (stderr, "Unable to contact server\n");
        goto FAILURE;
    }
    p->type = PROB_REMOTE;

    rval = CCutil_swrite_char (p->f, CCtsp_Pread);
    if (rval) {
        fprintf (stderr, "CCutil_swrite_char failed\n"); goto FAILURE;
    }
    rval = CCutil_swrite_string (p->f, pname);
    if (rval) {
        fprintf (stderr, "CCutil_swrite_string failed\n"); goto FAILURE;
    }
    rval = CCutil_swrite_int (p->f, n);
    if (rval) {
        fprintf (stderr, "CCutil_swrite_int failed\n"); goto FAILURE;
    }

    rval = prob_getheader (p, p);
    if (rval) {
        fprintf (stderr, "prob_getheader failed\n");
        goto FAILURE;
    }

    return p;
    
  FAILURE:
    if (p) {
        if (p->f) {
            CCutil_sclose (p->f);
        }
        CC_FREE (p, CCtsp_PROB_FILE);
    }
    return (CCtsp_PROB_FILE *) NULL;
}

#endif /* CC_NETREADY */

/* reads header data from p into the struct h (p == h except from copyheader */
static int prob_getheader (CCtsp_PROB_FILE *p, CCtsp_PROB_FILE *h)
{
    char version;
    int i;

    if (p->type == PROB_REMOTE) {
        if (CCutil_swrite_char (p->f, CCtsp_Pheader)) return -1;
    }
    
    if (CCutil_sread_char (p->f, &version)) return 1;

    switch (version) {
    case 1:
        for (i = 0; i < CCtsp_PROB_FILE_NAME_LEN; i++) {
            if (CCutil_sread_char (p->f, &h->name[i]))           return 1;
        }
        if (CCutil_sread_int      (p->f, &h->parent))            return 1;
        if (CCutil_sread_int      (p->f, &h->id))                return 1;
        if (CCutil_sread_double   (p->f, &h->ub))                return 1;
        if (CCutil_sread_double   (p->f, &h->lb))                return 1;
        if (CCbigguy_sread        (p->f, &h->exactlb))           return 1;
        if (CCutil_sread_int      (p->f, &h->nnodes))            return 1;
        if (CCutil_sread_int      (p->f, &h->child0))            return 1;
        if (CCutil_sread_int      (p->f, &h->child1))            return 1;
        if (CCutil_sread_int      (p->f, &h->real))              return 1;
        if (CCutil_sread_int      (p->f, &h->processed))         return 1;
        if (CCutil_sread_int      (p->f, &h->infeasible))        return 1;
        if (p->type == PROB_LOCAL) {
            if (CCutil_sread_int  (p->f, &h->offsets.dat))       return 1;
            if (CCutil_sread_int  (p->f, &h->offsets.edge))      return 1;
            if (CCutil_sread_int  (p->f, &h->offsets.fulladj))   return 1;
            if (CCutil_sread_int  (p->f, &h->offsets.cut))       return 1;
            if (CCutil_sread_int  (p->f, &h->offsets.tour))      return 1;
            if (CCutil_sread_int  (p->f, &h->offsets.fix))       return 1;
            if (CCutil_sread_int  (p->f, &h->offsets.exactdual)) return 1;
            if (CCutil_sread_int  (p->f, &h->offsets.history))   return 1;
            if (CCutil_sread_int  (p->f, &h->offsets.warmstart)) return 1;
        }
        break;
    default:
        fprintf (stderr, "Unknown problem version %ud\n", (unsigned) version);
        return 1;
    }
    return 0;
}

CCtsp_PROB_FILE *CCtsp_prob_write (char *f, int n)
{
    char nambuf[1024];

    if (remote_name (f)) {
#ifdef CC_NETREADY
        char hostbuf[1024];
        if (split_name (f, hostbuf, sizeof (hostbuf),
                        nambuf, sizeof (nambuf))) {
            fprintf (stderr, "Cannot split remote name\n");
            return (CCtsp_PROB_FILE *) NULL;
        }
        return CCtsp_prob_write_remote (hostbuf, nambuf, n);
#else /* CC_NETREADY */
        fprintf (stderr, "Remote problem writing not enabled\n");
        return (CCtsp_PROB_FILE *) NULL;
#endif
    }

    if (prob_name (nambuf, sizeof (nambuf), f, n))
        return (CCtsp_PROB_FILE *) NULL;

    return CCtsp_prob_write_name (nambuf);
}

CCtsp_PROB_FILE *CCtsp_prob_write_name (char *fname)
{
    CCtsp_PROB_FILE *p = (CCtsp_PROB_FILE *) NULL;

#if 0
    if (remote_name (fname)) {
        fprintf (stderr, "Cannot write remote problems by name\n");
        return (CCtsp_PROB_FILE *) NULL;
    }
#endif

    /*
    printf ("Write File %s at time %.0f\n", fname, CCutil_real_zeit());
    fflush (stdout);
    */

    p = CC_SAFE_MALLOC (1, CCtsp_PROB_FILE);
    if (p == (CCtsp_PROB_FILE *) NULL) goto FAILURE;
    prob_init (p);

    p->f = CCutil_sopen (fname, "w");
    if (!p->f) goto FAILURE;

    p->type = PROB_LOCAL;
    
    if (prob_putheader (p, p)) {
        fprintf (stderr, "prob_putheader failed\n");
        goto FAILURE;
    }

    return p;

FAILURE:
    if (p) {
        if (p->f) {
            CCutil_sclose (p->f);
        }
        CC_FREE (p, CCtsp_PROB_FILE);
    }
    return (CCtsp_PROB_FILE *) NULL;
}

#ifdef CC_NETREADY

CCtsp_PROB_FILE *CCtsp_prob_write_remote (char *hname, char *pname, int n)
{
    CCtsp_PROB_FILE *p = (CCtsp_PROB_FILE *) NULL;
    int i;
    int rval;

    /* Hack for run at Princeton without .caam.rice.edu */

/*
    START HACK
    char tname[4048];
    sprintf (tname, "%s.caam.rice.edu", hname);
    END HACK - also change hname to tname
*/

    /*
    printf ("Write Remote Host %s name %s id %d\n", hname, pname, n);
    fflush (stdout);
    */

    p = CC_SAFE_MALLOC (1, CCtsp_PROB_FILE);
    if (p == (CCtsp_PROB_FILE *) NULL) goto FAILURE;
    prob_init (p);

    for (i = 0; pname[i] && i < CCtsp_PROB_FILE_NAME_LEN - 1; i++)
        p->name[i] = pname[i];
    p->name[i] = '\0';

    p->f = CCutil_snet_open (hname, CCtsp_PROB_PORT);
    if (p->f == (CC_SFILE *) NULL) {
        fprintf (stderr, "Unable to contact server\n");
        goto FAILURE;
    }
    p->type = PROB_REMOTE;

    rval = CCutil_swrite_char (p->f, CCtsp_Pwrite);
    if (rval) {
        fprintf (stderr, "CCutil_swrite_char failed\n"); goto FAILURE;
    }
    rval = CCutil_swrite_string (p->f, pname);
    if (rval) {
        fprintf (stderr, "CCutil_swrite_string failed\n"); goto FAILURE;
    }
    rval = CCutil_swrite_int (p->f, n);
    if (rval) {
        fprintf (stderr, "CCutil_swrite_int failed\n"); goto FAILURE;
    }
    
    return p;

FAILURE:
    if (p) {
        if (p->f) {
            CCutil_sclose (p->f);
        }
        CC_FREE (p, CCtsp_PROB_FILE);
    }
    return (CCtsp_PROB_FILE *) NULL;
}

int CCtsp_prob_delete_remote (char *hname, char *pname, int n)
{
    CC_SFILE *f = (CC_SFILE *) NULL;
    int rval;

    printf ("Delete Remote Host %s name %s id %d\n", hname, pname, n);
    fflush (stdout);

    f = CCutil_snet_open (hname, CCtsp_PROB_PORT);
    if (f == (CC_SFILE *) NULL) {
        fprintf (stderr, "Unable to contact server\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_swrite_char (f, CCtsp_Pdelete);
    if (rval) {
        fprintf (stderr, "CCutil_swrite_char failed\n"); goto CLEANUP;
    }
    rval = CCutil_swrite_string (f, pname);
    if (rval) {
        fprintf (stderr, "CCutil_swrite_string failed\n"); goto CLEANUP;
    }
    rval = CCutil_swrite_int (f, n);
    if (rval) {
        fprintf (stderr, "CCutil_swrite_int failed\n"); goto CLEANUP;
    }

    rval = CCutil_sclose (f);
    if (rval) {
        fprintf (stderr, "CCutil_sclose failed\n");
    }
    f = (CC_SFILE *) NULL;

    rval = 0;
  CLEANUP:
    if (f != (CC_SFILE *) NULL) {
        CCutil_sclose (f);
    }
    return rval;
}

CCtsp_PROB_FILE *CCtsp_prob_server (CC_SFILE *s)
{
    CCtsp_PROB_FILE *p = (CCtsp_PROB_FILE *) NULL;

    p = CC_SAFE_MALLOC (1, CCtsp_PROB_FILE);
    if (p == (CCtsp_PROB_FILE *) NULL) {
        return (CCtsp_PROB_FILE *) NULL;
    }
    prob_init (p);

    p->f = s;
    p->type = PROB_SERVER;

    return p;
}

#endif /* CC_NETREADY */

/* writes header data from struct h to p (p == h except from copyheader */
static int prob_putheader (CCtsp_PROB_FILE *p, CCtsp_PROB_FILE *h)
{
    int  i;

    if (p->type == PROB_LOCAL) {
        if (CCutil_srewind (p->f))                           return 1;
    } else if (p->type == PROB_REMOTE) {
        if (CCutil_swrite_char (p->f, CCtsp_Pheader))        return 1;
    }

    if (CCutil_swrite_char (p->f, PROB_HEADER_VERSION))
        return 1;
    for (i = 0; i < CCtsp_PROB_FILE_NAME_LEN; i++) {
        if (CCutil_swrite_char (p->f, h->name[i]))           return 1;
    }
    if (CCutil_swrite_int      (p->f, h->parent))            return 1;
    if (CCutil_swrite_int      (p->f, h->id))                return 1;
    if (CCutil_swrite_double   (p->f, h->ub))                return 1;
    if (CCutil_swrite_double   (p->f, h->lb))                return 1;
    if (CCbigguy_swrite        (p->f, h->exactlb))           return 1;
    if (CCutil_swrite_int      (p->f, h->nnodes))            return 1;
    if (CCutil_swrite_int      (p->f, h->child0))            return 1;
    if (CCutil_swrite_int      (p->f, h->child1))            return 1;
    if (CCutil_swrite_int      (p->f, h->real))              return 1;
    if (CCutil_swrite_int      (p->f, h->processed))         return 1;
    if (CCutil_swrite_int      (p->f, h->infeasible))        return 1;
    if (p->type == PROB_LOCAL) {
        if (CCutil_swrite_int  (p->f, h->offsets.dat))       return 1;
        if (CCutil_swrite_int  (p->f, h->offsets.edge))      return 1;
        if (CCutil_swrite_int  (p->f, h->offsets.fulladj))   return 1;
        if (CCutil_swrite_int  (p->f, h->offsets.cut))       return 1;
        if (CCutil_swrite_int  (p->f, h->offsets.tour))      return 1;
        if (CCutil_swrite_int  (p->f, h->offsets.fix))       return 1;
        if (CCutil_swrite_int  (p->f, h->offsets.exactdual)) return 1;
        if (CCutil_swrite_int  (p->f, h->offsets.history))   return 1;
        if (CCutil_swrite_int  (p->f, h->offsets.warmstart)) return 1;
    }

    return 0;
}

static int prob_name (char *buf, size_t buflen, char *f, int n)
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
        if (n >= NAMESTRNUM('x','x','x')) n++;
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

static void prob_init (CCtsp_PROB_FILE *p)
{
    int i;

    p->f = (CC_SFILE *) NULL;
    for (i = 0; i < CCtsp_PROB_FILE_NAME_LEN; i++)
        p->name[i] = '\0';
    p->type = 0;
    p->id = -1;
    p->parent = -1;
    p->lb = -1.0;
    p->ub = -1.0;
    p->exactlb = CCbigguy_ZERO;
    p->nnodes = -1;
    p->child0 = -1;
    p->child1 = -1;
    p->real = -1;
    p->processed = -1;
    p->infeasible = -1;
    p->offsets.dat = -1;
    p->offsets.edge = -1;
    p->offsets.fulladj = -1;
    p->offsets.cut = -1;
    p->offsets.tour = -1;
    p->offsets.basis = -1;
    p->offsets.norms = -1;
    p->offsets.fix = -1;
    p->offsets.exactdual = -1;
    p->offsets.history = -1;
    p->offsets.warmstart = -1;
}

int CCtsp_prob_rclose (CCtsp_PROB_FILE *p)
{
    int rval = 0;

    if (p->type == PROB_REMOTE) {
        rval |= CCutil_swrite_char (p->f, CCtsp_Pexit);
    }

    if (p->type != PROB_SERVER) {
        rval |= CCutil_sclose (p->f);
    }

    CC_FREE (p, CCtsp_PROB_FILE);

    if (rval) return 1;
    else return 0;
}

int CCtsp_prob_wclose (CCtsp_PROB_FILE *p)
{
    int rval = 0;

    rval |= prob_putheader (p, p);
    if (p->type == PROB_REMOTE) {
        rval |= CCutil_swrite_char (p->f, CCtsp_Pexit);
    }

    if (p->type != PROB_SERVER) {
        rval |= CCutil_sclose (p->f);
    }

    CC_FREE (p, CCtsp_PROB_FILE);

    return rval;
}

int CCtsp_prob_getname (CCtsp_PROB_FILE *p, char *name)
{
    int i;

    if (!p) return -1;

    for (i = 0; i < CCtsp_PROB_FILE_NAME_LEN; i++)
        name[i] = p->name[i];

    return 0;
}

int CCtsp_prob_putname (CCtsp_PROB_FILE *p, char *name)
{
    int i;

    if (!p) return 1;

    for (i = 0; name[i] && i < CCtsp_PROB_FILE_NAME_LEN - 1; i++)
        p->name[i] = name[i];
    for (; i < CCtsp_PROB_FILE_NAME_LEN; i++)
        p->name[i] = '\0';

    return 0;
}


int CCtsp_prob_getid (CCtsp_PROB_FILE *p, int *id)
{
    if (!p) return -1;

    *id = p->id;
    if (*id == -1) {
        printf ("Setting -1 ID to 0\n"); fflush (stdout);
        *id = 0;
    }
    return 0;
}

int CCtsp_prob_putid (CCtsp_PROB_FILE *p, int id)
{
    if (!p) return 1;

    p->id = id;
    return 0;
}

int CCtsp_prob_getparent (CCtsp_PROB_FILE *p, int *parent)
{
    if (!p) return -1;

    *parent = p->parent;
    return 0;
}

int CCtsp_prob_putparent (CCtsp_PROB_FILE *p, int parent)
{
    if (!p) return 1;

    p->parent = parent;
    return 0;
}

int CCtsp_prob_getub (CCtsp_PROB_FILE *p, double *ub)
{
    if (!p) return -1;

    *ub = p->ub;
    return 0;
}

int CCtsp_prob_putub (CCtsp_PROB_FILE *p, double ub)
{
    if (!p) return 1;

    p->ub = ub;
    return 0;
}

int CCtsp_prob_getlb (CCtsp_PROB_FILE *p, double *lb)
{
    if (!p) return -1;

    *lb = p->lb;
    return 0;
}

int CCtsp_prob_putlb (CCtsp_PROB_FILE *p, double lb)
{
    if (!p) return 1;

    p->lb = lb;
    return 0;
}

int CCtsp_prob_getexactlb (CCtsp_PROB_FILE *p, CCbigguy *lb)
{
    if (!p) return -1;

    *lb = p->exactlb;
    return 0;
}

int CCtsp_prob_putexactlb (CCtsp_PROB_FILE *p, CCbigguy lb)
{
    if (!p) return 1;

    p->exactlb = lb;
    return 0;
}

int CCtsp_prob_getnnodes (CCtsp_PROB_FILE *p, int *nnodes)
{
    if (!p) return -1;

    *nnodes = p->nnodes;
    return 0;
}

int CCtsp_prob_putnnodes (CCtsp_PROB_FILE *p, int nnodes)
{
    if (!p) return 1;

    p->nnodes = nnodes;
    return 0;
}

int CCtsp_prob_getchildren (CCtsp_PROB_FILE *p, int *child0, int *child1)
{
    if (!p) return -1;

    *child0 = p->child0;
    *child1 = p->child1;
    return 0;
}

int CCtsp_prob_putchildren (CCtsp_PROB_FILE *p, int child0, int child1)
{
    if (!p) return 1;

    p->child0 = child0;
    p->child1 = child1;
    return 0;
}

int CCtsp_prob_getreal (CCtsp_PROB_FILE *p, int *real)
{
    if (!p) return -1;

    *real = p->real;
    return 0;
}

int CCtsp_prob_putreal (CCtsp_PROB_FILE *p, int real)
{
    if (!p) return 1;

    p->real = real;
    return 0;
}

int CCtsp_prob_getprocessed (CCtsp_PROB_FILE *p, int *processed)
{
    if (!p) return -1;

    *processed = p->processed;
    return 0;
}

int CCtsp_prob_putprocessed (CCtsp_PROB_FILE *p, int processed)
{
    if (!p) return 1;

    p->processed = processed;
    return 0;
}

int CCtsp_prob_getinfeasible (CCtsp_PROB_FILE *p, int *infeasible)
{
    if (!p) return -1;

    *infeasible = p->infeasible;
    return 0;
}

int CCtsp_prob_putinfeasible (CCtsp_PROB_FILE *p, int infeasible)
{
    if (!p) return 1;

    p->infeasible = infeasible;
    return 0;
}

static int prob_copyheader (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t)
{
    int rval;

    if (f->type == PROB_SERVER && t->type == PROB_SERVER) {
        fprintf (stderr, "Cannot copy header from server to server\n");
        return -1;
    } else if (f->type == PROB_SERVER) {
        rval = prob_getheader (f, t);
        if (rval) fprintf (stderr, "prob_getheader failed\n");
        return rval;
    } else if (t->type == PROB_SERVER) {
        rval = prob_putheader (t, f);
        if (rval) fprintf (stderr, "prob_putheader failed\n");
        return rval;
    } else {
        strcpy (t->name, f->name);
        t->id = f->id;
        t->parent = f->parent;
        t->ub = f->ub;
        t->lb = f->lb;
        t->exactlb = f->exactlb;
        t->nnodes = f->nnodes;
        t->child0 = f->child0;
        t->child1 = f->child1;
        t->real = f->real;
        t->processed = f->processed;
        t->infeasible = f->infeasible;
        return 1;
    }
}

static int begin_put (CCtsp_PROB_FILE *p, int *offset, char section)
{
    if (p->type == PROB_LOCAL) {
        *offset = CCutil_stell (p->f);
    } else if (p->type == PROB_REMOTE) {
        if (CCutil_swrite_char (p->f, section)) return 1;
    }

    return 0;
}

static int begin_get (CCtsp_PROB_FILE *p, int offset, char section, int silent)
{
    if (p->type == PROB_LOCAL) {
        if (offset == -1) {
            if (!silent) {
                printf ("No section %c in file.\n", section);
                fflush (stdout);
            }
            return 1;
        }
        if (CCutil_sseek (p->f, offset)) {
            fprintf (stderr, "CCutil_sseek failed in begin_get\n");
            fflush (stdout);
            return -1;
        }
    } else if (p->type == PROB_REMOTE) {
        char exists;
        
        if (CCutil_swrite_char (p->f, section)) return -1;
        if (CCutil_sread_char (p->f, &exists)) return -1;
        if (exists == 0) {
            if (!silent) {
                printf ("No section %c in remote file\n", section);
                fflush (stdout);
            }
            return 1;
        }
    }
    return 0;
}

static int begin_copy (CCtsp_PROB_FILE *f, int foffset, CCtsp_PROB_FILE *t,
        int *toffset, char section, int silent)
{
    char exists;

    if (f->type == PROB_LOCAL) {
        if (foffset == -1) {
            if (!silent) {
                printf ("No section %c in file.\n", section);
                fflush (stdout);
            }
            exists = 0;
        } else {
            if (CCutil_sseek (f->f, foffset)) {
                fprintf (stderr, "CCutil_sseek failed in begin_copy\n");
                return -1;
            }
            exists = 1;
        }
    } else if (f->type == PROB_REMOTE) {
        if (CCutil_swrite_char (f->f, section)) return -1;
        if (CCutil_sread_char (f->f, &exists)) return -1;
        if (exists == 0) {
            if (!silent) {
                printf ("No section %c in remote file\n", section);
                fflush (stdout);
            }
        }
    } else if (f->type == PROB_SERVER) {
        exists = 1;
    }

    if (t->type == PROB_SERVER) {
        if (CCutil_swrite_char (t->f, exists)) return -1;
    }
    if (exists) {
        if (t->type == PROB_LOCAL) {
            *toffset = CCutil_stell (t->f);
        } else if (t->type == PROB_REMOTE) {
            if (CCutil_swrite_char (t->f, section)) return -1;
        }
        return 0;
    } else {
        return 1;
    }
}

int CCtsp_prob_puttour (CCtsp_PROB_FILE *p, int ncount, int *tour)
{
    int i;
    int rval;
    int nbits;

    if (!p) return 1;

    rval = begin_put (p, &p->offsets.tour, CCtsp_Ptour);
    if (rval) return rval;

    if (CCutil_swrite_char     (p->f, PROB_TOUR_VERSION)) return 1;
    if (CCutil_swrite_int      (p->f, ncount))            return 1;

    nbits = CCutil_sbits (ncount);
    
    for (i = 0; i < ncount; i++) {
        if (CCutil_swrite_bits (p->f, tour[i], nbits))    return 1;
    }

    return 0;
}

int CCtsp_prob_gettour (CCtsp_PROB_FILE *p, int ncount, int **tour, int silent)
{
    int i;
    int rval;
    int nbits;
    char version;
    int ncount2;

    if (!p) return -1;

    rval = begin_get (p, p->offsets.tour, CCtsp_Ptour, silent);
    if (rval) return rval;

    *tour = (int *) NULL;
    
    if (CCutil_sread_char (p->f, &version)) goto FAILURE;
    switch (version) {
    case 1:
        if (CCutil_sread_int (p->f, &ncount2)) goto FAILURE;
        if (ncount != ncount) {
            fprintf (stderr, "Wrong ncount in tour\n");
            goto FAILURE;
        }
        nbits = CCutil_sbits (ncount2);
        *tour = CC_SAFE_MALLOC (ncount2, int);
        if (!(*tour)) {
            fprintf (stderr, "out of memory in CCtsp_prob_gettour\n");
            goto FAILURE;
        }

        for (i = 0; i < ncount2; i++) {
            if (CCutil_sread_bits (p->f, &((*tour)[i]), nbits)) {
                goto FAILURE;
            }
        }
        break;
    default:
        fprintf (stderr, "Unknown tour version %ud\n", (unsigned) version);
        goto FAILURE;
    }
    return 0;

 FAILURE:
    CC_IFFREE (*tour, int);
    return -1;
}

static int prob_copytour (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t, int silent)
{
    int ncount;
    int rval;
    char version;
    int nbits;

    rval = begin_copy (f, f->offsets.tour, t, &t->offsets.tour, CCtsp_Ptour,
                       silent);
    if (rval == 1) return 0;
    else if (rval) return rval;

    rval = copy_char_val (f->f, t->f, (char *) &version);
    if (rval) return -1;
    
    switch (version) {
    case 1:
        rval = copy_int_val (f->f, t->f, &ncount);
        if (rval) return -1;
        nbits = CCutil_sbits (ncount);
        rval = copy_bits (f->f, t->f, ncount, nbits);
        if (rval) return rval;
        break;
    default:
        fprintf (stderr, "Unknown tour version %ud\n", (unsigned) version);
        return 1;
    }
    return 0;
}

int CCtsp_prob_putedges (CCtsp_PROB_FILE *p, int ncount, int ecount,
        int *elist, int *elen)
{
    int i;
    int rval;
    int nbits;

    if (!p) return 1;

    rval = begin_put (p, &p->offsets.edge, CCtsp_Pedges);
    if (rval) return rval;

    if (CCutil_swrite_char (p->f, PROB_EDGES_VERSION)) return 1;

    if (CCutil_swrite_int (p->f, ncount)) return 1;
    if (CCutil_swrite_int (p->f, ecount)) return 1;
    nbits = CCutil_sbits (ncount);

    for (i = 0; i < ecount; i++) {
        if (CCutil_swrite_bits (p->f, elist[2*i], nbits)) return 1;
        if (CCutil_swrite_bits (p->f, elist[2*i+1], nbits)) return 1;
        if (CCutil_swrite_int (p->f, elen[i])) return 1;
    }

    return 0;
}

int CCtsp_prob_getedges (CCtsp_PROB_FILE *p, int ncount, int *ecount,
        int **elist, int **elen, int silent)
{
    int i;
    int rval;
    char version;
    int nbits;
    int ncount2;

    if (!p) return -1;

    *elist = (int *) NULL;
    *elen = (int *) NULL;
    
    rval = begin_get (p, p->offsets.edge, CCtsp_Pedges, silent);
    if (rval) return rval;

    if (CCutil_sread_char (p->f, &version)) goto FAILURE;
    switch (version) {
    case 1:
        if (CCutil_sread_int (p->f, &ncount2)) goto FAILURE;
        if (ncount2 != ncount) {
            fprintf (stderr, "Wrong ncount in edges\n");
            goto FAILURE;
        }
        if (CCutil_sread_int (p->f, ecount)) goto FAILURE;
        nbits = CCutil_sbits (ncount2);
        *elist = CC_SAFE_MALLOC (2 * (*ecount), int);
        *elen = CC_SAFE_MALLOC (*ecount, int);
        if (!(*elist) || !(*elen)) {
            fprintf (stderr, "out of memory in CCtsp_prob_getedges\n");
            goto FAILURE;
        }

        for (i = 0; i < (*ecount); i++) {
            if (CCutil_sread_bits (p->f, &((*elist)[2*i]), nbits)) {
                goto FAILURE;
            }
            if (CCutil_sread_bits (p->f, &((*elist)[2*i+1]), nbits)) {
                goto FAILURE;
            }
            if (CCutil_sread_int (p->f, &((*elen)[i]))) {
                goto FAILURE;
            }
        }
        break;
    default:
        fprintf (stderr, "Unknown edges version %ud\n", (unsigned) version);
        goto FAILURE;
    }
    return 0;

 FAILURE:
    CC_IFFREE (*elist, int);
    CC_IFFREE (*elen, int);
    return -1;
    
}

static int prob_copyedges (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t, int silent)
{
    int ecount;
    int ncount;
    int rval;
    char version;
    int nbits;
    int i;

    rval = begin_copy (f, f->offsets.edge, t, &t->offsets.edge, CCtsp_Pedges,
                       silent);
    if (rval == 1) return 0;
    else if (rval) return rval;

    rval = copy_char_val (f->f, t->f, (char *) &version);
    if (rval) return rval;
    switch (version) {
    case 1:
        rval = copy_int_val (f->f, t->f, &ncount);
        if (rval) return rval;
        rval = copy_int_val (f->f, t->f, &ecount);
        if (rval) return rval;
        nbits = CCutil_sbits (ncount);
        for (i=0; i<ecount; i++) {
            rval = copy_bits (f->f, t->f, 2, nbits);
            if (rval) return rval;
            rval = copy_ints (f->f, t->f, 1);
            if (rval) return rval;
        }
        break;
    default:
        fprintf (stderr, "Unknown edges version %ud\n", (unsigned) version);
        return 1;
    }
    return 0;
}

int CCtsp_prob_putcuts (CCtsp_PROB_FILE *p, int ncount, CCtsp_lpcuts *cuts)
{
    int rval;

    if (!p) return 1;

    rval = begin_put (p, &p->offsets.cut, CCtsp_Pcuts);
    if (rval) return rval;

    rval = CCtsp_write_cuts (p->f, ncount, cuts, 1);
    if (rval) return rval;

    return 0;
}

int CCtsp_prob_getcuts (CCtsp_PROB_FILE *p, int *ncount, CCtsp_lpcuts *cuts,
        int silent)
{
    int rval;

    if (!p) return -1;

    rval = begin_get (p, p->offsets.cut, CCtsp_Pcuts, silent);
    if (rval) return rval;

    rval = CCtsp_read_cuts (p->f, ncount, cuts, 1, 0);
    if (rval) return rval;

    return 0;
}

static int prob_copycuts (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t, int silent)
{
    int rval;
    
    rval = begin_copy (f, f->offsets.cut, t, &t->offsets.cut, CCtsp_Pcuts,
                       silent);
    if (rval == 1) return 0;
    else if (rval) return rval;

    rval = CCtsp_copy_cuts (f->f, t->f, 1);
    if (rval) return rval;

    return 0;
}

int CCtsp_prob_putwarmstart (CCtsp_PROB_FILE *p, CClp_warmstart *w)
{
    int rval;
    
    if (!p) return 1;

    rval = begin_put (p, &p->offsets.warmstart, CCtsp_Pwarmstart);
    if (rval) return rval;

    if (CClp_swrite_warmstart (p->f, w)) {
        fprintf (stderr, "CClp_swrite_warmstart failed\n");
        return 1;
    }

    return 0;
}

int CCtsp_prob_getwarmstart (CCtsp_PROB_FILE *p, CClp_warmstart **w, int silent)
{
    int rval;

    if (!p) return -1;

    rval = begin_get (p, p->offsets.warmstart, CCtsp_Pwarmstart, silent);
    if (rval) return rval;


    if (CClp_sread_warmstart (p->f, w)) {
        fprintf (stderr, "CClp_sread_warmstart failed\n");;
        return -1;
    }
    return 0;
}

static int prob_copywarmstart (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t,
        int silent)
{
    CClp_warmstart *w = (CClp_warmstart *) NULL;
    int rval;
    
    rval = begin_copy (f, f->offsets.warmstart, t, &t->offsets.warmstart,
                       CCtsp_Pwarmstart, silent);
    if (rval == 1) return 0;
    else if (rval) return rval;

    rval = CClp_sread_warmstart (f->f, &w);
    if (rval) return rval;
    rval = CClp_swrite_warmstart (t->f, w);
    if (rval) return rval;

    CClp_free_warmstart (&w);

    return 0;
}

int CCtsp_prob_putfulladj (CCtsp_PROB_FILE *p, int ncount, int fullcount,
            CCtsp_genadj *adj)
{
    int i, j;
    int rval;
    int nbits;

    if (!p) return 1;

    rval = begin_put (p, &p->offsets.fulladj, CCtsp_Pfull);
    if (rval) return rval;

    if (CCutil_swrite_char (p->f, PROB_FULLADJ_VERSION)) return 1;

    if (CCutil_swrite_int (p->f, ncount)) return 1;
    if (CCutil_swrite_int (p->f, fullcount)) return 1;
    nbits = CCutil_sbits (ncount);
    for (i = 0; i < ncount; i++) {
        if (CCutil_swrite_bits (p->f, adj[i].deg, nbits)) return 1;
        for (j = 0; j < adj[i].deg; j++) {
            if (CCutil_swrite_bits (p->f, adj[i].list[j].end, nbits)) return 1;
            if (CCutil_swrite_int (p->f, adj[i].list[j].len)) return 1;
        }
    }
    return 0;
}

int CCtsp_prob_getfulladj (CCtsp_PROB_FILE *p, int ncount, int *fullcount,
            CCtsp_genadj **adj, CCtsp_genadjobj **adjspace, int silent)
{
    int i, j;
    CCtsp_genadj *a;
    CCtsp_genadjobj *s;
    int rval;
    int ncount2;
    char version;
    int nbits;

    *fullcount = 0;
    *adj = (CCtsp_genadj *) NULL;
    *adjspace = (CCtsp_genadjobj *) NULL;

    if (!p) return -1;

    rval = begin_get (p, p->offsets.fulladj, CCtsp_Pfull, silent);
    if (rval) return rval;

    if (CCutil_sread_char (p->f, &version)) goto FAILURE;
    switch (version) {
    case 1:
        if (CCutil_sread_int (p->f, &ncount2)) goto FAILURE;
        if (ncount != ncount2) {
            fprintf (stderr, "ncount incorrect in fulladj\n"); goto FAILURE;
        }
        nbits = CCutil_sbits (ncount2);
        if (CCutil_sread_int (p->f, fullcount)) goto FAILURE;
        *adjspace = CC_SAFE_MALLOC (*fullcount, CCtsp_genadjobj);
        *adj = CC_SAFE_MALLOC (ncount, CCtsp_genadj);
        if (!adjspace || !adj) {
            fprintf (stderr, "out of memory in CCtsp_prob_getfulladj\n");
            goto FAILURE;
        }

        s = *adjspace;
        a = *adj;
        for (i = 0; i < ncount; i++) {
            a[i].list = s;
            if (CCutil_sread_bits (p->f, &(a[i].deg), nbits)) goto FAILURE;
            for (j = 0; j < a[i].deg; j++) {
                if (CCutil_sread_bits (p->f, &(a[i].list[j].end), nbits))
                    goto FAILURE;
                if (CCutil_sread_int (p->f, &(a[i].list[j].len)))
                    goto FAILURE;
            }
            s += a[i].deg;
        }
        break;
    default:
        fprintf (stderr, "Unknown fulladj version %ud\n", (unsigned) version);
        goto FAILURE;
    }
    return 0;

FAILURE:
    CC_IFFREE (*adjspace, CCtsp_genadjobj);
    CC_IFFREE (*adj, CCtsp_genadj);
    *fullcount = 0;
    return -1;
}

static int prob_copyfulladj (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t, int silent)
{
    int rval;
    int deg;
    int ncount;
    char version;
    int nbits;
    int i;
    int j;

    rval = begin_copy (f, f->offsets.fulladj, t, &t->offsets.fulladj,
                       CCtsp_Pfull, silent);
    if (rval == 1) return 0;
    else if (rval) return rval;

    rval = copy_char_val (f->f, t->f, (char *) &version);
    if (rval) return rval;
    switch (version) {
    case 1:
        rval = copy_int_val (f->f, t->f, &ncount);
        if (rval) return rval;
        rval = copy_ints (f->f, t->f, 1);
        if (rval) return rval;
        nbits = CCutil_sbits (ncount);
        for (i = 0; i < ncount; i++) {
            rval = copy_bits_val (f->f, t->f, &deg, nbits);
            if (rval) return rval;
            for (j=0; j<deg; j++) {
                rval = copy_bits (f->f, t->f, 1, nbits);
                if (rval) return rval;
                rval = copy_ints (f->f, t->f, 1);
                if (rval) return rval;
            }
        }
        break;
    default:
        fprintf (stderr, "Unknown fulladj version %ud\n", (unsigned) version);
        return -1;
    }
    return 0;
}

int CCtsp_prob_putfixed (CCtsp_PROB_FILE *p, int ncount, int ecount,
        int *elist)
{
    int i;
    int rval;
    int nbits;

    if (!p) return 1;

    rval = begin_put (p, &p->offsets.fix, CCtsp_Pfixed);
    if (rval) return rval;

    if (CCutil_swrite_char (p->f, PROB_FIXED_VERSION)) return 1;

    if (CCutil_swrite_int (p->f, ncount)) return 1;
    if (CCutil_swrite_int (p->f, ecount)) return 1;

    nbits = CCutil_sbits (ncount);
    
    for (i = 0; i < 2*ecount; i++) {
        if (CCutil_swrite_bits (p->f, elist[i], nbits)) return 1;
    }

    return 0;
}

int CCtsp_prob_getfixed (CCtsp_PROB_FILE *p, int ncount, int *ecount,
        int **elist, int silent)
{
    int i;
    int rval;
    int ncount2;
    char version;
    int nbits;

    *ecount = 0;
    *elist = (int *) NULL;
    if (!p) return -1;

    rval = begin_get (p, p->offsets.fix, CCtsp_Pfixed, silent);
    if (rval) return rval;

    if (CCutil_sread_char (p->f, &version)) goto FAILURE;
    switch (version) {
    case 1:
        if (CCutil_sread_int (p->f, &ncount2)) goto FAILURE;
        if (ncount != ncount2) {
            fprintf (stderr, "wrong ncount in fixed edges\n");
            goto FAILURE;
        }
        if (CCutil_sread_int (p->f, ecount)) goto FAILURE;

        if (*ecount) {
            *elist = CC_SAFE_MALLOC (2*(*ecount), int);
            if (!(*elist)) {
                fprintf (stderr, "out of memory in CCtsp_prob_getfixed\n");
                goto FAILURE;
            }
            nbits = CCutil_sbits (ncount);
            for (i = 0; i < 2*(*ecount); i++) {
                if (CCutil_sread_bits (p->f, (&((*elist)[i])), nbits))
                    goto FAILURE;
            }
        }
        break;
    default:
        fprintf (stderr, "Unknown fixed version %ud\n", (unsigned) version);
        goto FAILURE;
    }
    return 0;

FAILURE:
    CC_IFFREE (*elist, int);
    *ecount = 0;
    return -1;
}

static int prob_copyfixed (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t, int silent)
{
    int ecount;
    int ncount;
    int nbits;
    char version;
    int rval;

    rval = begin_copy (f, f->offsets.fix, t, &t->offsets.fix, CCtsp_Pfixed,
                       silent);
    if (rval == 1) return 0;
    else if (rval) return rval;

    rval = copy_char_val (f->f, t->f, (char *) &version);
    if (rval) return rval;
    switch (version) {
    case 1:
        rval = copy_int_val (f->f, t->f, &ncount);
        if (rval) return rval;
        rval = copy_int_val (f->f, t->f, &ecount);
        if (rval) return rval;
        nbits = CCutil_sbits (ncount);
        rval = copy_bits (f->f, t->f, ecount*2, nbits);
        if (rval) return rval;
        break;
    default:
        fprintf (stderr, "Unknown fixed version %ud\n", (unsigned) version);
        return -1;
    }
    return 0;
}

int CCtsp_prob_putexactdual (CCtsp_PROB_FILE *p, CCtsp_bigdual *d, int ncount)
{
    int i;
    int rval;

    if (!p) return 1;

    rval = begin_put (p, &p->offsets.exactdual, CCtsp_Pdual);
    if (rval) return rval;

    if (CCutil_swrite_char (p->f, PROB_EXACTDUAL_VERSION)) return 1;
    if (CCutil_swrite_int (p->f, ncount)) return 1;
    if (CCutil_swrite_int (p->f, d->cutcount)) return 1;

    for (i = 0; i < ncount; i++) {
        if (CCbigguy_swrite (p->f, d->node_pi[i])) return 1;
    }
    for (i = 0; i < d->cutcount; i++) {
        if (CCbigguy_swrite (p->f, d->cut_pi[i])) return 1;
    }

    return 0;
}


int CCtsp_prob_getexactdual (CCtsp_PROB_FILE *p, int ncount, CCtsp_bigdual **d,
        int silent)
{
    CCtsp_bigdual *rd;
    int i;
    int rval;
    int ncount2;
    char version;

    *d = (CCtsp_bigdual *) NULL;

    if (!p) return -1;

    rval = begin_get (p, p->offsets.exactdual, CCtsp_Pdual, silent);
    if (rval) return rval;

    if (CCutil_sread_char (p->f, &version)) return -1;
    switch (version) {
    case 1:
        if (CCutil_sread_int (p->f, &ncount2)) goto FAILURE;
        if (ncount != ncount2) {
            fprintf (stderr, "wrong ncount in exact dual\n"); goto FAILURE;
        }
    
        *d = CC_SAFE_MALLOC (1, CCtsp_bigdual);
        if (!(*d)) {
            fprintf (stderr, "out of memory in CCtsp_prob_getexactdual\n");
            goto FAILURE;
        }
        rd = *d;
        rd->cutcount = 0;
        rd->node_pi = (CCbigguy *) NULL;
        rd->cut_pi = (CCbigguy *) NULL;
        
        if (CCutil_sread_int (p->f, &(rd->cutcount))) goto FAILURE;
        
        rd->node_pi = CC_SAFE_MALLOC (ncount, CCbigguy);
        if (!rd->node_pi) {
            fprintf (stderr, "out of memory in CCtsp_prob_getexactdual\n");
            goto FAILURE;
        }
        for (i = 0; i < ncount; i++) {
            if (CCbigguy_sread (p->f, &(rd->node_pi[i]))) goto FAILURE;
        }
        if (rd->cutcount) {
            rd->cut_pi = CC_SAFE_MALLOC (rd->cutcount, CCbigguy);
            if (!rd->cut_pi) {
                fprintf (stderr, "out of memory in CCtsp_prob_getexactdual\n");
                goto FAILURE;
            }
            for (i = 0; i < rd->cutcount; i++) {
                if (CCbigguy_sread (p->f, &(rd->cut_pi[i]))) goto FAILURE;
            }
        }
        break;
    default:
        fprintf (stderr, "Unknown exact dual version %ud\n", (unsigned) version);
        goto FAILURE;
    }
    return 0;

FAILURE:
    if (*d) {
        CC_IFFREE ((*d)->node_pi, CCbigguy);
        CC_IFFREE ((*d)->cut_pi, CCbigguy);
        CC_FREE (*d, CCtsp_bigdual);
    }
    return -1;
}

static int prob_copyexactdual (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t,
       int silent)
{
    int rval;
    int ncount;
    int cutcount;
    char version;

    rval = begin_copy (f, f->offsets.exactdual, t, &t->offsets.exactdual,
                       CCtsp_Pdual, silent);
    if (rval == 1) return 0;
    else if (rval) return rval;

    rval = copy_char_val (f->f, t->f, (char *) &version);
    if (rval) return rval;
    switch (version) {
    case 1:
        rval = copy_int_val (f->f, t->f, &ncount);
        if (rval) return rval;
        rval = copy_int_val (f->f, t->f, &cutcount);
        if (rval) return rval;
        rval = copy_bigguys (f->f, t->f, ncount);
        if (rval) return rval;
        rval = copy_bigguys (f->f, t->f, cutcount);
        if (rval) return rval;
        return 0;
    default:
        fprintf (stderr, "Unknown exact dual version %ud\n", (unsigned) version);
        return -1;
    }
}


int CCtsp_prob_puthistory (CCtsp_PROB_FILE *p, int depth,
        CCtsp_branchobj *history)
{
    int i, j;
    CCtsp_lpclique *c;
    int rval;

    if (!p) return 1;

    rval = begin_put (p, &p->offsets.history, CCtsp_Phistory);
    if (rval) return rval;

    if (CCutil_swrite_char (p->f, PROB_HISTORY_VERSION)) return 1;
    if (CCutil_swrite_int (p->f, depth)) return 1;

    for (i = 0; i < depth; i++) {
        if (CCutil_swrite_int (p->f, history[i].depth)) return 1;
        if (CCutil_swrite_int (p->f, history[i].rhs)) return 1;
        if (CCutil_swrite_int (p->f, history[i].ends[0])) return 1;
        if (CCutil_swrite_int (p->f, history[i].ends[1])) return 1;
        if (history[i].clique) {
            c = history[i].clique;
            if (CCutil_swrite_int (p->f, c->segcount)) return 1;
            for (j = 0; j < c->segcount; j++) {
                if (CCutil_swrite_int (p->f, c->nodes[j].lo)) return 1;
                if (CCutil_swrite_int (p->f, c->nodes[j].hi)) return 1;
            }
        } else {
            if (CCutil_swrite_int (p->f, 0)) return 1;
        }
        if (CCutil_swrite_char (p->f, history[i].sense)) return 1;
    }

    return 0;
}

int CCtsp_prob_gethistory (CCtsp_PROB_FILE *p, int *depth,
        CCtsp_branchobj **history, int silent)
{
    int i, j, nseg, rval;
    CCtsp_lpclique *c = (CCtsp_lpclique *) NULL;
    int *slist  = (int *) NULL;
    char version;

    *depth = 0;
    *history = (CCtsp_branchobj *) NULL;
    if (!p) return -1;

    rval = begin_get (p, p->offsets.history, CCtsp_Phistory, silent);
    if (rval) return rval;

    if (CCutil_sread_char (p->f, &version)) goto FAILURE;
    switch (version) {
    case 1:
        if (CCutil_sread_int (p->f, depth)) goto FAILURE;

        if (*depth) {
            *history = CC_SAFE_MALLOC (*depth, CCtsp_branchobj);
            if (!(*history)) {
                fprintf (stderr, "out of memory in CCtsp_prob_gethistory\n");
                goto FAILURE;
            }
            for (i = 0; i < (*depth); i++) {
                if (CCutil_sread_int (p->f, &((*history)[i].depth)))
                    goto FAILURE;
                if (CCutil_sread_int (p->f, &((*history)[i].rhs)))
                    goto FAILURE;
                if (CCutil_sread_int (p->f, &((*history)[i].ends[0])))
                    goto FAILURE;
                if (CCutil_sread_int (p->f, &((*history)[i].ends[1])))
                    goto FAILURE;
                if (CCutil_sread_int (p->f, &nseg))
                    goto FAILURE;
                if (nseg) {
                    slist = CC_SAFE_MALLOC (2*nseg, int);
                    if (!slist) goto FAILURE;
                    for (j = 0; j < nseg; j++) {
                        if (CCutil_sread_int (p->f, &slist[2*j]))
                            goto FAILURE;
                        if (CCutil_sread_int (p->f, &slist[2*j+1]))
                            goto FAILURE;
                    }
                    c = CC_SAFE_MALLOC (1, CCtsp_lpclique);
                    if (!c) goto FAILURE;
                    rval = CCtsp_seglist_to_lpclique (nseg, slist, c);
                    if (rval) {
                        fprintf (stderr, "CCtsp_seglist_to_lpclique failed\n");
                        CC_FREE (c, CCtsp_lpclique);
                        goto FAILURE;
                    }
                    (*history)[i].clique = c;
                } else {
                    (*history)[i].clique = (CCtsp_lpclique *) NULL;
                }
                if (CCutil_sread_char (p->f, &((*history)[i].sense)))
                    return 1;
            }
        }
        break;
    default:
        fprintf (stderr, "Unknown history version %ud\n", (unsigned) version);
        goto FAILURE;
    }
    return 0;

FAILURE:
    CC_IFFREE (slist, int);
    CC_IFFREE (*history, CCtsp_branchobj);
    *depth = 0;
    return -1;
}

static int prob_copyhistory (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t, int silent)
{
    int rval;
    int depth;
    int segcount;
    int i;
    char version;

    rval = begin_copy (f, f->offsets.history, t, &t->offsets.history,
                       CCtsp_Phistory, silent);
    if (rval == 1) return 0;
    else if (rval) return rval;

    rval = copy_char_val (f->f, t->f, (char *) &version);
    if (rval) return rval;
    switch (version) {
    case 1:
        rval = copy_int_val (f->f, t->f, &depth);
        if (rval) return rval;

        for (i=0; i<depth; i++) {
            rval = copy_ints (f->f, t->f, 4);
            if (rval) return rval;
            rval = copy_int_val (f->f, t->f, &segcount);
            if (rval) return rval;
            if (segcount != 0) {
                rval = copy_ints (f->f, t->f, 2*segcount);
                if (rval) return rval;
            }
            rval = copy_chars (f->f, t->f, 1);
            if (rval) return rval;
        }
        break;
    default:
        fprintf (stderr, "Unknown history version %ud\n", (unsigned) version);
        return -1;
    }
    return 0;
}

int CCtsp_prob_copy_section (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t,
        char section, int silent)
{
    switch (section) {
      case CCtsp_Pheader:
        return prob_copyheader (f, t);
      case CCtsp_Ptour:
        return prob_copytour (f, t, silent);
      case CCtsp_Pedges:
        return prob_copyedges (f, t, silent);
      case CCtsp_Pcuts:
        return prob_copycuts (f, t, silent);
      case CCtsp_Pwarmstart:
        return prob_copywarmstart (f, t, silent);
      case CCtsp_Pfull:
        return prob_copyfulladj (f, t, silent);
      case CCtsp_Pfixed:
        return prob_copyfixed (f, t, silent);
      case CCtsp_Pdual:
        return prob_copyexactdual (f, t, silent);
      case CCtsp_Phistory:
        return prob_copyhistory (f, t, silent);
      default:
        fprintf (stderr, "Invalid section %c in CCtsp_prob_copy_section\n",
                 section);
        return 1;
    }
}

static int remote_name (const char *f)
{
    return (CCutil_strchr_c (f, ':') != (const char *) NULL);
}

static int split_name (const char *f, char *hostname, size_t hlen,
        char *probname, size_t plen)
{
    size_t len;
    const char *p;
    
    p = CCutil_strchr_c (f, ':');
    if (p == (const char *) NULL) {
        fprintf (stderr, "non-net name in split_name\n");
        return 1;
    }

    len = p - f;
    if (len+1 > hlen) {
        fprintf (stderr, "hostname too long in split_name\n");
        return 1;
    }

    strncpy (hostname, f, len);
    hostname[len] = '\0';

    len = strlen (p+1);
    if (len+1 > plen) {
        fprintf (stderr, "filename too long in split_name\n");
        return 1;
    }

    strncpy (probname, p+1, len);
    probname[len] = '\0';
    return 0;
}

static int copy_ints (CC_SFILE *f, CC_SFILE *t, int n)
{
    int rval;
    int i;
    int x;

    for (i=0; i<n; i++) {
        rval = CCutil_sread_int (f, &x);
        if (rval) return rval;
        rval = CCutil_swrite_int (t, x);
        if (rval) return rval;
    }
    return 0;
}

static int copy_int_val (CC_SFILE *f, CC_SFILE *t, int *n)
{
    int rval;

    rval = CCutil_sread_int (f, n);
    if (rval) return rval;
    rval = CCutil_swrite_int (t, *n);
    if (rval) return rval;
    return 0;
}

static int copy_bits (CC_SFILE *f, CC_SFILE *t, int n, int nbits)
{
    int rval;
    int i;
    int x;

    for (i=0; i<n; i++) {
        rval = CCutil_sread_bits (f, &x, nbits);
        if (rval) return rval;
        rval = CCutil_swrite_bits (t, x, nbits);
        if (rval) return rval;
    }
    return 0;
}

static int copy_bits_val (CC_SFILE *f, CC_SFILE *t, int *n, int nbits)
{
    int rval;

    rval = CCutil_sread_bits (f, n, nbits);
    if (rval) return rval;
    rval = CCutil_swrite_bits (t, *n, nbits);
    if (rval) return rval;
    return 0;
}

static int copy_chars (CC_SFILE *f, CC_SFILE *t, int n)
{
    int rval;
    int i;
    char x;

    for (i=0; i<n; i++) {
        rval = CCutil_sread_char (f, &x);
        if (rval) return rval;
        rval = CCutil_swrite_char (t, x);
        if (rval) return rval;
    }
    return 0;
}

static int copy_char_val (CC_SFILE *f, CC_SFILE *t, char *n)
{
    int rval;

    rval = CCutil_sread_char (f, n);
    if (rval) return rval;
    rval = CCutil_swrite_char (t, *n);
    if (rval) return rval;
    return 0;
}

static int copy_bigguys (CC_SFILE *f, CC_SFILE *t, int n)
{
    int rval;
    int i;
    CCbigguy x;

    for (i=0; i<n; i++) {
        rval = CCbigguy_sread (f, &x);
        if (rval) return rval;
        rval = CCbigguy_swrite (t, x);
        if (rval) return rval;
    }
    return 0;
}

char *CCtsp_problabel (const char *probloc)
{
    const char *p;
    const char *problabel = probloc;
    char *probcopy = (char *) NULL;
    char *p2;

    p = CCutil_strrchr_c (problabel, ':');
    if (p != (const char *) NULL) problabel = p+1;
    p = CCutil_strrchr_c (problabel, '/');
    if (p != (const char *) NULL) problabel = p+1;
    probcopy = CCutil_strdup (problabel);
    if (probcopy == (char *) NULL) return (char *) NULL;
    p2 = CCutil_strchr (probcopy, '.');
    if (p2 != (char *) NULL) *p2 = '\0';
    return probcopy;
}

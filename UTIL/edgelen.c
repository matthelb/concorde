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
/*       FUNCTIONS FOR COMPUTING EDGE LENGTHS FOR GEOMETRIC PROBLEMS        */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: Summer 1994                                                       */
/*        Modified - March 2, 1995                                          */
/*                 - October 5, 1995 (Bico)                                 */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_dat_setnorm (CCdatagroup *dat, int norm)                     */
/*     NOTES:                                                               */
/*         Supported norms (with defs in edgelen.h) are:                    */
/*             CC_MAXNORM  -  the L-infinity norm                           */
/*             CC_EUCLIDEAN_CEIL - the norm for the plaXXXX problems        */
/*             CC_EUCLIDEAN - rounded L-2 norm                              */
/*             CC_EUCLIDEAN_3D - rounded L-2 norm in 3 space                */
/*             CC_USER - a norm specified by the user                       */
/*             CC_GEOGRAPHIC - distances on a sphere (Groetshel and         */
/*                             Holland)                                     */
/*             CC_GEOM - sphere in meters, coords in decimal degrees,       */
/*                             slight modification of CC_GEOGRAPHIC         */
/*             CC_ATT - pseudo-Euclidean norm for att532                    */
/*             CC_MATRIXNORM - complete graph (lower + diagonal matrix)     */
/*             CC_DSJRANDNORM - random edgelengths                          */
/*             CC_CRYSTAL - Bland-Shallcross xray norm                      */
/*                     - The coordinates generated for CC_CRYSTAL problems  */
/*                (in CCutil_getdata.c) have been diveded by the motor      */
/*                speeds (this makes the edgelen function faster) and       */
/*                scaled by CRYSTAL_SCALE (currently 10000) and rounded to  */
/*                the nearest integer (this lets the edgelen function       */
/*                produce integer lengths without further rounding). The    */
/*                result is a closer approximation to the Bland -           */
/*                Shallcross floating point length function than that       */
/*                given in TSPLIB_1.2.                                      */
/*             CC_SPARSE - a sparse graph                                   */
/*             CC_RHMAPx - where x = 1, 2, 3, 4, 5 one of 5 RH mapping      */
/*                norms.                                                    */
/*                                                                          */
/*         If CCUTIL_EDGELEN_FUNCTIONPTR has been defined in util.h,        */
/*         then CCutil_dat_edgelen is a pointer to a function instead of    */
/*         a function.  This saves a function call and results in           */
/*         improved performance on some machines for edgelen-intensive      */
/*         routines like linkern.  The function pointer is set by           */
/*         CCutil_dat_setnorm.                                              */
/*                                                                          */
/*         IMPORTANT: This means that if CCUTIL_EDGELEN_FUNCTIONPTR is set  */
/*         and you have more than one CCdatagroup, you must call            */
/*         CCutil_dat_setnorm whenever you switch from using one            */
/*         CCdatagroup to the other.  IF YOU DON'T DO THIS, EDGELEN WILL    */
/*         RETURN INCORRECT RESULTS.  For this reason,                      */
/*         CCUTIL_EDGELEN_FUNCTIONPTR should only be set with extreme       */
/*         caution.                                                         */
/*                                                                          */
/*         NOTE: CCUTIL_EDGELEN_FUNCTIONPTR does not work with the          */
/*         subdivision code for parallel TSP.                               */
/*                                                                          */
/*    To define a user norm, you must perform the following steps:          */
/*    1.  In util.h, define the struct CCdata_user to contain the data      */
/*        necessary for the computation of edge lengths.                    */
/*    2.  In edgelen.c, write the init_userdat and free_userdat functions   */
/*        which initialize and free a CCdata_user structure.                */
/*    3.  In edgelen.c, write the user_edgelen function which               */
/*        computes the length of the edge for node i to node j, using the   */
/*        userdat field of the CCdatagroup argument (userdat is of type     */
/*        CCdata_user).                                                     */
/*    4.  In getdata.c, write the build_user, read_user_text,               */
/*        read_user_binary, readmaster_user, and writemaster_user           */
/*        routines.  read_user_text reads the data file which provides      */
/*        the data for computing the edge lengths.  build_user and          */
/*        read_user_binary are optional routines which build random         */
/*        datasets and read binary datafiles.  writemaster_user writes a    */
/*        binary version of that data to the master file, and               */
/*        readmaster_user reads that same data.  See the comments before    */
/*        those routines in getdata for more details on what they should    */
/*        do.                                                               */
/*    5.  In getdata.c, write permute_user, which permutes the data to      */
/*        reflect a permutation of the nodes.                               */
/*                                                                          */
/*  int CCutil_dat_edgelen (int i, int j, CCdatagroup *dat)                 */
/*     compute the length of an edge                                        */
/*                                                                          */
/*  void CCutil_dsjrand_init (CCdatagroup *dat, int maxdist, int seed)      */
/*     initialize the dsjrand norm                                          */
/*                                                                          */
/*  void CCutil_dat_getnorm (CCdatagroup *dat, int *norm)                   */
/*     get the norm of a CCdatagroup                                        */
/*                                                                          */
/*  void CCutil_init_datagroup (CCdatagroup *dat)                           */
/*     initialize a CCdatagroup                                             */
/*                                                                          */
/*  void CCutil_freedatagroup (CCdatagroup *dat)                            */
/*     free a CCdatagroup                                                   */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "macrorus.h"


static double
    dtrunc (double);

static int
    edgelen_nonorm (int i, int j, CCdatagroup *dat),
    max_edgelen (int i, int j, CCdatagroup *dat),
    man_edgelen (int i, int j, CCdatagroup *dat),
    euclid_edgelen (int i, int j, CCdatagroup *dat),
    euclid_ceiling_edgelen (int i, int j, CCdatagroup *dat),
    euclid3d_edgelen (int i, int j, CCdatagroup *dat),
    geographic_edgelen (int i, int j, CCdatagroup *dat),
    geom_edgelen (int i, int j, CCdatagroup *dat),
    att_edgelen (int i, int j, CCdatagroup *dat),
    dsjrand_edgelen (int i, int j, CCdatagroup *dat),
    crystal_edgelen (int i, int j, CCdatagroup *dat),
    matrix_edgelen (int i, int j, CCdatagroup *dat),
    sparse_edgelen (int i, int j, CCdatagroup *dat),
    user_edgelen (int i, int j, CCdatagroup *dat),
    rhmap1_edgelen (int i, int j, CCdatagroup *dat),
    rhmap2_edgelen (int i, int j, CCdatagroup *dat),
    rhmap3_edgelen (int i, int j, CCdatagroup *dat),
    rhmap4_edgelen (int i, int j, CCdatagroup *dat),
    rhmap5_edgelen (int i, int j, CCdatagroup *dat),
    toroidal_edgelen (int i, int j, CCdatagroup *dat);

static void
    init_userdat (CCdata_user *userdat),
    free_userdat (CCdata_user *userdat),
    init_rhdata (CCdata_rhvector *rhdat),
    free_rhdata (CCdata_rhvector *rhdat);


#ifndef M_PI
#define M_PI 3.14159265358979323846264
#endif

#ifdef  CCUTIL_EDGELEN_FUNCTIONPTR

int (*CCutil_dat_edgelen) (int i, int j, CCdatagroup *dat) = edgelen_nonorm;

#else /* CCUTIL_EDGELEN_FUNCTIONPTR */

int CCutil_dat_edgelen (int i, int j, CCdatagroup *dat)
{
    if (dat->ndepot) {
        if (i >= dat->orig_ncount) {
            return dat->depotcost[j];
        } else if (j >= dat->orig_ncount) {
            return dat->depotcost[i];
        }
    }
    return (dat->edgelen)(i, j, dat);
}

#endif /* CCUTIL_EDGELEN_FUNCTIONPTR */


int CCutil_dat_setnorm (CCdatagroup *dat, int norm)
{
    switch (norm) {
    case CC_EUCLIDEAN_CEIL:
        dat->edgelen = euclid_ceiling_edgelen;
        break;
    case CC_EUCLIDEAN:
        dat->edgelen = euclid_edgelen;
        break;
    case CC_MAXNORM:
        dat->edgelen = max_edgelen;
        break;
    case CC_MANNORM:
        dat->edgelen = man_edgelen;
        break;
    case CC_EUCLIDEAN_3D:
        dat->edgelen = euclid3d_edgelen;
        break;
    case CC_USER:
        dat->edgelen = user_edgelen;
        break;
    case CC_GEOGRAPHIC:
        dat->edgelen = geographic_edgelen;
        break;
    case CC_GEOM:
        dat->edgelen = geom_edgelen;
        break;
    case CC_ATT:
        dat->edgelen = att_edgelen;
        break;
    case CC_MATRIXNORM:
        dat->edgelen = matrix_edgelen;
        break;
    case CC_DSJRANDNORM:
        dat->edgelen = dsjrand_edgelen;
        break;
    case CC_CRYSTAL:
        dat->edgelen = crystal_edgelen;
        break;
    case CC_SPARSE:
        dat->edgelen = sparse_edgelen;
        break;
    case CC_RHMAP1:
        dat->edgelen = rhmap1_edgelen;
        break;
    case CC_RHMAP2:
        dat->edgelen = rhmap2_edgelen;
        break;
    case CC_RHMAP3:
        dat->edgelen = rhmap3_edgelen;
        break;
    case CC_RHMAP4:
        dat->edgelen = rhmap4_edgelen;
        break;
    case CC_RHMAP5:
        dat->edgelen = rhmap5_edgelen;
        break;
    case CC_EUCTOROIDAL:
        dat->edgelen = toroidal_edgelen;
        break;
    default:
        fprintf (stderr, "ERROR:  Unknown NORM %d.\n", norm);
        return 1;
    }
    dat->norm = norm;

#ifdef CCUTIL_EDGELEN_FUNCTIONPTR
    CCutil_dat_edgelen = dat->edgelen;
#endif /* CCUTIL_EDGELEN_FUNCTIONPTR */

    return 0;
}

void CCutil_dat_getnorm (CCdatagroup *dat, int *norm)
{
    (*norm) = dat->norm;
}

static int edgelen_nonorm (int i, int j, CCdatagroup *dat)
{
    fprintf (stderr, "CCutil_dat_edgelen has been called with no norm set\n");
    fprintf (stderr, "This is a FATAL ERROR\n");
    if (i != 0 || j != 0 || dat != (CCdatagroup *) NULL) {
        /* so the compiler won't complain about unused variables */
        fprintf (stderr, "This is a FATAL ERROR\n");
        exit (1);
    }
    return -1;
}

/* Several variables that would normally be called y1 and y2 are called
   yy1 and yy2 to avoid conflict with the bessel functions */

static int max_edgelen (int i, int j, CCdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j], t2 = dat->y[i] - dat->y[j];

    if (t1 < 0)
        t1 *= -1;
    if (t2 < 0)
        t2 *= -1;
    t1 += 0.5;
    t2 += 0.5;

    return (int) (t1 < t2 ? t2 : t1);
}

static int man_edgelen (int i, int j, CCdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j], t2 = dat->y[i] - dat->y[j];

    if (t1 < 0)
        t1 *= -1;
    if (t2 < 0)
        t2 *= -1;

    return (int) (t1 + t2 + 0.5);
}


static int euclid_edgelen (int i, int j, CCdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j], t2 = dat->y[i] - dat->y[j];
    int temp;

    temp = (int) (sqrt (t1 * t1 + t2 * t2) + 0.5);
    return temp;
}

static int toroidal_edgelen (int i, int j, CCdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j];
    double t2 = dat->y[i] - dat->y[j];
    int temp;

    if (t1 < 0) t1 = -t1;
    if (t2 < 0) t2 = -t2;
    if (dat->gridsize - t1 < t1) t1 = dat->gridsize - t1;
    if (dat->gridsize - t2 < t2) t2 = dat->gridsize - t2;
    temp = (int) (sqrt (t1 * t1 + t2 * t2) + 0.5);
    return temp;
}

static int euclid3d_edgelen (int i, int j, CCdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j], t2 = dat->y[i] - dat->y[j];
    double t3 = dat->z[i] - dat->z[j];
    int temp;

    temp = (int) (sqrt (t1 * t1 + t2 * t2 + t3 * t3) + 0.5);
    return temp;
}

static int euclid_ceiling_edgelen (int i, int j, CCdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j], t2 = dat->y[i] - dat->y[j];
/*
    int rd;
    double max;

    max = sqrt (t1 * t1 + t2 * t2);
    rd = (int) max;
    return (((max - rd) > .000000001) ? rd + 1 : rd);
*/
    return (int) (ceil (sqrt (t1 * t1 + t2 * t2)));
}

#define GH_PI (3.141592)

static int geographic_edgelen (int i, int j, CCdatagroup *dat)
{
    double deg, min;
    double lati, latj, longi, longj;
    double q1, q2, q3;
    int dd;
    double x1 = dat->x[i], x2 = dat->x[j], yy1 = dat->y[i], yy2 = dat->y[j];

    deg = dtrunc (x1);
    min = x1 - deg;
    lati = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc (x2);
    min = x2 - deg;
    latj = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;

    deg = dtrunc (yy1);
    min = yy1 - deg;
    longi = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc (yy2);
    min = yy2 - deg;
    longj = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;

    q1 = cos (longi - longj);
    q2 = cos (lati - latj);
    q3 = cos (lati + latj);
    dd = (int) (6378.388 * acos (0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3))
                + 1.0);
    return dd;
}

static int geom_edgelen (int i, int j, CCdatagroup *dat)
{
    double lati, latj, longi, longj;
    double q1, q2, q3, q4, q5;

    lati = M_PI * dat->x[i] / 180.0;
    latj = M_PI * dat->x[j] / 180.0;

    longi = M_PI * dat->y[i] / 180.0;
    longj = M_PI * dat->y[j] / 180.0;

    q1 = cos (latj) * sin(longi - longj);
    q3 = sin((longi - longj)/2.0);
    q4 = cos((longi - longj)/2.0);
    q2 = sin(lati + latj) * q3 * q3 - sin(lati - latj) * q4 * q4;
    q5 = cos(lati - latj) * q4 * q4 - cos(lati + latj) * q3 * q3;
    return (int) (6378388.0 * atan2(sqrt(q1*q1 + q2*q2), q5) + 1.0);
}

#if 0
static int geom_edgelen (int i, int j, CCdatagroup *dat)
{
    double lati, latj, longi, longj;
    double q1, q2, q3;
    int dd;

    lati = M_PI * (dat->x[i] / 180.0);
    latj = M_PI * (dat->x[j] / 180.0);

    longi = M_PI * (dat->y[i] / 180.0);
    longj = M_PI * (dat->y[j] / 180.0);

    q1 = cos (longi - longj);
    q2 = cos (lati - latj);
    q3 = cos (lati + latj);
    dd = (int) (6378388.0 * acos (0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3))
                + 1.0);
    return dd;
}
#endif

static int att_edgelen (int i, int j, CCdatagroup *dat)
{
    double xd = dat->x[i] - dat->x[j];
    double yd = dat->y[i] - dat->y[j];
    double rij = sqrt ((xd * xd + yd * yd) / 10.0);
    double tij = dtrunc (rij);
    int dij;

    if (tij < rij)
        dij = (int) tij + 1;
    else
        dij = (int) tij;
    return dij;
}

static double dtrunc (double x)
{
    int k;

    k = (int) x;
    x = (double) k;
    return x;
}

void CCutil_dsjrand_init (CCdatagroup *dat, int maxdist, int seed)
{
    dat->dsjrand_factor = maxdist/2147483648.0;
    dat->dsjrand_param = 104*seed+1;
}

static int dsjrand_edgelen (int i, int j, CCdatagroup *dat)
{
    int di = (int) dat->x[i];
    int dj = (int) dat->x[j];
    int x, y, z;

    x = di&dj;
    y = di|dj;
    z = dat->dsjrand_param;

    x *= z;
    y *= x;
    z *= y;

    z ^= dat->dsjrand_param;

    x *= z;
    y *= x;
    z *= y;

    x = ((di+dj)^z)&0x7fffffff;
    return (int)(x * dat->dsjrand_factor);
}

#define CRYSTAL_SCALE 10000

#define CRYSTAL_FLIP_TOL ((180 * CRYSTAL_SCALE * 4) / 5)
#define CRYSTAL_NEEDS_FLIP(x) ((x) > (CRYSTAL_FLIP_TOL))
#define CRYSTAL_FLIP(x) ((2 * (CRYSTAL_FLIP_TOL)) - (x))

static int crystal_edgelen (int i, int j, CCdatagroup *dat)
{
    double w, w1;

    w = dat->x[i] - dat->x[j];
    if (w < 0)
        w = -w;
    w1 = dat->y[i] - dat->y[j];
    if (w1 < 0)
        w1 = -w1;
    if (CRYSTAL_NEEDS_FLIP (w1))
        w1 = CRYSTAL_FLIP (w1);
    if (w < w1)
        w = w1;
    w1 = dat->z[i] - dat->z[j];
    if (w1 < 0)
        w1 = -w1;
    if (w < w1)
        w = w1;

    return (int) w;
}

static int matrix_edgelen (int i, int j, CCdatagroup *dat)
{
    if (i > j)
        return (dat->adj[i])[j];
    else
        return (dat->adj[j])[i];
}

static int sparse_edgelen (int i, int j, CCdatagroup *dat)
{
    int *adj;
    int k, deg;

    if (i > j) {
        CC_SWAP (i, j, k);
    }
    adj = dat->adj[i];
    deg = dat->degree[i];

    for (k = 0; k < deg; k++) {
        if (adj[k] == j) {
            return dat->len[i][k];
        }
    }
    return dat->default_len;
}

void CCutil_init_datagroup (CCdatagroup *dat)
{
    dat->x = (double *) NULL;
    dat->y = (double *) NULL;
    dat->z = (double *) NULL;
    dat->adj = (int **) NULL;
    dat->adjspace = (int *) NULL;
    dat->len      = (int **) NULL;
    dat->lenspace = (int *) NULL;
    dat->degree   = (int *) NULL;
    dat->norm = 0;
    dat->dsjrand_param = 1;
    dat->dsjrand_factor = 1.0;
    dat->default_len = 100000;
    dat->sparse_ecount = 0;
    dat->edgelen = edgelen_nonorm;
    init_userdat (&dat->userdat);
    init_rhdata (&dat->rhdat);
    dat->ndepot = 0;
    dat->orig_ncount = 0;
    dat->depotcost = (int *) NULL;
    dat->orig_names = (int *) NULL;
}

void CCutil_freedatagroup (CCdatagroup *dat)
{
    CC_IFFREE (dat->x, double);
    CC_IFFREE (dat->y, double);
    CC_IFFREE (dat->z, double);
    CC_IFFREE (dat->adj, int *);
    CC_IFFREE (dat->adjspace, int);
    CC_IFFREE (dat->len, int *);
    CC_IFFREE (dat->lenspace, int);
    CC_IFFREE (dat->degree, int);
    free_userdat (&dat->userdat);
    free_rhdata (&dat->rhdat);
    CC_IFFREE (dat->depotcost, int);
    CC_IFFREE (dat->orig_names, int);
}

static void init_userdat (CCdata_user *userdat)
{
    userdat->x = (double *) NULL;
    userdat->y = (double *) NULL;
}

static void free_userdat (CCdata_user *userdat)
{
    CC_IFFREE (userdat->x, double);
    CC_IFFREE (userdat->y, double);
}

static int user_edgelen (int i, int j, CCdatagroup *dat)
{
    double dw = dat->userdat.x[i] - dat->userdat.x[j];
    double dw1 = dat->userdat.y[i] - dat->userdat.y[j];
    static const double ibm_xmult[7] = {1062.5,
        300.0,
        300.0,
        250.0,
        300.0,
        1000.0,
        154.6};
    static const double ibm_xadd[7] = {155.0 - 0.01 * 1062.5,
        197.5 - 0.05 * 300.0,
        212.5 - 0.10 * 300.0,
        227.5 - 0.15 * 250.0,
        240.5 - 0.20 * 300.0,
        255.0 - 0.25 * 1000.0,
        305.0 - 0.30 * 154.6};
    static const double ibm_ymult[7] = {1062.5,
        450.0,
        350.0,
        250.0,
        300.0,
        900.0,
        157.7};
    static const double ibm_yadd[7] = {150.0 - 0.01 * 1062.5,
        192.5 - 0.05 * 450.0,
        215.0 - 0.10 * 350.0,
        232.5 - 0.15 * 250.0,
        245.5 - 0.20 * 300.0,
        250.0 - 0.25 * 900.0,
        295.0 - 0.30 * 157.7};

    if (dw < 0.0)
        dw = -dw;
    dw /= 25400.0;
    if (dw <= 0.01) {
        dw *= 15500.0;
    } else if (dw >= 0.30) {
        dw = dw * 154.6 + (305.0 - 0.3 * 154.6);
    } else {
        dw = dw * ibm_xmult[(int) (dw / 0.05)] +
            ibm_xadd[(int) (dw / 0.05)];
    }
    if (dw1 < 0.0)
        dw1 = -dw1;
    dw1 /= 25400.0;
    if (dw1 <= 0.01) {
        dw1 *= 15000.0;
    } else if (dw1 >= 0.30) {
        dw1 = dw1 * 157.7 + (295.0 - 0.3 * 157.7);
    } else {
        dw1 = dw1 * ibm_ymult[(int) (dw1 / 0.05)] +
            ibm_yadd[(int) (dw1 / 0.05)];
    }
    if (dw < dw1)
        dw = dw1;
    return (int) dw;
}

static void init_rhdata (CCdata_rhvector *rhdat)
{
    rhdat->space = (char *) NULL;
    rhdat->vectors = (char **) NULL;
    rhdat->rhlength = 0;
    rhdat->dist_00 = 0;
    rhdat->dist_01 = 0;
    rhdat->dist_02 = 0;
    rhdat->dist_22 = 0;
    rhdat->p       = 0.0;
}

static void free_rhdata (CCdata_rhvector *rhdat)
{
    CC_IFFREE (rhdat->space, char);
    CC_IFFREE (rhdat->vectors, char *);
    rhdat->rhlength = 0;
}

static int rhmap1_edgelen (int i, int j, CCdatagroup *dat)
{
    char **vectors = dat->rhdat.vectors;
    int rhlength = dat->rhdat.rhlength;
    char *v1 = vectors[i];
    char *v2 = vectors[j];
    int n;
    int sum = 0;
    int dist_00 = dat->rhdat.dist_00;
    int dist_01 = dat->rhdat.dist_01;
    int dist_02 = dat->rhdat.dist_02;
    int dist_12 = dat->rhdat.dist_12;
    int dist_22 = dat->rhdat.dist_22;

    if (v1 == (char *) NULL || v2 == (char *) NULL) return 0;
    
    for (n=0; n<rhlength; n++) {
        if (v1[n] == 2) {
            if (v2[n] == 0) sum += dist_02;
            else if (v2[n] == 1) sum += dist_12;
            else sum += dist_22;
        } else {
            if (v1[n] == v2[n]) sum += dist_00;
            else if (v2[n] != 2) sum += dist_01;
            else if (v1[n] == 0) sum += dist_02;
            else sum += dist_12;
        }
    }
    return sum;
}

static int rhmap2_edgelen (int i, int j, CCdatagroup *dat)
{
    char **vectors = dat->rhdat.vectors;
    int rhlength = dat->rhdat.rhlength;
    char *v1 = vectors[i];
    char *v2 = vectors[j];
    int n;
    double sum = 0;
    double p = dat->rhdat.p;

    if (v1 == (char *) NULL || v2 == (char *) NULL) return 0;
    
    for (n=0; n<rhlength; n++) {
        if (v1[n] == 0) {
            if (v2[n] == 1) sum += 1;
            else if (v2[n] == 2) sum += p;
        } else if (v1[n] == 1) {
            if (v2[n] == 0) sum += 1;
            else if (v2[n] == 2) sum += (1-p);
        } else {
            if (v2[n] == 0) sum += p;
            else if (v2[n] == 1) sum += (1-p);
            else sum += 2*p*(1-p);
        }
    }
    return (int) (sum * 100.0 + 0.5);
}

#define MAX_DIST 1000

static int rhmap3_edgelen (int i, int j, CCdatagroup *dat)
{
    char **vectors = dat->rhdat.vectors;
    int len = dat->rhdat.rhlength;
    char *first = vectors[i];
    char *second = vectors[j];
    int xindex;
    int a = 0;
    int b = 0;
    int c = 0;
    int d = 0;
    int n = 0;
    double P = dat->rhdat.p;
    double Q = 1.0 - P;
    double trans;
    double temp;
    double theta;
    double term;

    if (first == (char *) NULL) {
        if (second == (char *) NULL) return 0;
        first = second;
        second = (char *) NULL;
    }
    
    if (second == (char *) NULL) {
        for (xindex = 0; xindex < len; xindex++) {
            if (first[xindex] == 1) a++;
            else if (first[xindex] == 0) b++;
        }
        trans = pow(sqrt(P), (double) a) * pow(sqrt(Q), (double) b);
    } else {
        for (xindex = 0; xindex < len; xindex++) {
            if ((first[xindex] != 2) || (second[xindex] != 2)) {
                n++;
                if ((first[xindex] == 1) && (second[xindex] == 1)) a++;
                if ((first[xindex] == 1) && (second[xindex] == 0)) b++;
                if ((first[xindex] == 0) && (second[xindex] == 1)) c++;
                if ((first[xindex] == 0) && (second[xindex] == 0)) d++;
            }
        }
        if (n == 0) return MAX_DIST;
        temp = (n - (a*P) - (d*Q));
        term = (4.0 * n * P * Q * (b+c));

        if (term >= (temp * temp)) return MAX_DIST;

        theta = (temp - sqrt((temp * temp) - term)) / (2.0 * n * P * Q);

        if (theta >= 1.0) return MAX_DIST;

        trans = pow((1.0 - (theta*P)),  (double) d) *
                pow((1.0 - (theta*Q)),  (double) a) *
                pow((theta * sqrt(P*Q)), (double)  (b+c));

    }
    return ((int) (-10.0 * log10(trans)));
}

static int rhmap4_edgelen (int i, int j, CCdatagroup *dat)
{
    char **vectors = dat->rhdat.vectors;
    int len = dat->rhdat.rhlength;
    char *first = vectors[i];
    char *second = vectors[j];
    int xindex;
    int a = 0;
    int b = 0;
    int c = 0;
    int d = 0;
    int n = 0;
    double P = dat->rhdat.p;
    double Q = 1.0 - P;
    double trans;
    double temp;
    double theta;
    double term;

    if (first == (char *) NULL) {
        if (second == (char *) NULL) return 0;
        first = second;
        second = (char *) NULL;
    }
    
    if (second == (char *) NULL) {
        for (xindex = 0; xindex < len; xindex++) {
            if (first[xindex] == 1) a++;
            else if (first[xindex] == 0) b++;
        }
        trans = pow(sqrt(P), (double) a) * pow(sqrt(Q), (double) b);
    } else {
        for (xindex = 0; xindex < len; xindex++) {
            if ((first[xindex] != 2) && (second[xindex] != 2)) {
                n++;
                if ((first[xindex] == 1) && (second[xindex] == 1)) a++;
                if ((first[xindex] == 1) && (second[xindex] == 0)) b++;
                if ((first[xindex] == 0) && (second[xindex] == 1)) c++;
                if ((first[xindex] == 0) && (second[xindex] == 0)) d++;
            }
        }
        if (n == 0) return MAX_DIST;
        temp = (n - (a*P) - (d*Q));
        term = (4.0 * n * P * Q * (b+c));

        if (term >= (temp * temp)) return MAX_DIST;

        theta = (temp - sqrt((temp * temp) - term)) / (2.0 * n * P * Q);

        if (theta >= 1.0) return MAX_DIST;

        trans = pow((1.0 - (theta*P)),  (double) d) *
                pow((1.0 - (theta*Q)),  (double) a) *
                pow((theta * sqrt(P*Q)),  (double) (b+c));

    }
    return ((int) (-10.0 * log10(trans)));
}

static int rhmap5_edgelen (int i, int j, CCdatagroup *dat)
{
    char **vectors = dat->rhdat.vectors;
    int rhlength = dat->rhdat.rhlength;
    char *v1 = vectors[i];
    char *v2 = vectors[j];
    int n;
    int mis = 0;
    int cnt = 0;

    if (v1 == (char *) NULL || v2 == (char *) NULL) return 0;
    
    for (n=0; n<rhlength; n++) {
        if (v1[n] != 2 && v2[n] != 2) {
            cnt++;
            if (v1[n] != v2[n]) mis++;
        }
    }
    if (cnt == 0) return 0;
    else return (int) (10.0 * rhlength * mis) / cnt;
}


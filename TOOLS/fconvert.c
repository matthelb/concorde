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
/*    UTIL PROGRAM FOR CONVERTING TO AND FROM BINARY DAT AND EDGE FILES     */
/*                                                                          */
/*                              TSP CODE                                    */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 7, 1995 (cofeb16)                                        */
/*                                                                          */
/*  SEE short decsription in usage ().                                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

static char *edgefilename = (char *) NULL;
static char *datfilename = (char *) NULL;
static int integerdata_in = 0;
static int binary_in = 0;
static int binary_out = 1;
static int reversed_int = 0;
static int double_len = 0;


int
    main (int ac, char **av);
static void
    usage (char *f);
static int
    parseargs (int ac, char **av),
    getgraph (char *edgename, char *datname, int *ncount, int *ecount,
              int **elist, int **elen, double **xcoord, double **ycoord),
    dumpedges (int ncount, int ecount, int *elist, int *elen),
    dumpdat (int ncount, double *xcoord, double *ycoord);



int main (int ac, char **av)
{
    double szeit;
    int ncount = 0, ecount = 0;
    int *elist = (int *) NULL, *elen = (int *) NULL;
    double *dlen = (double *) NULL;
    double *xcoord = (double *) NULL, *ycoord = (double *) NULL;

    if (parseargs (ac, av))
        return 1;
    if (edgefilename == (char *) NULL && datfilename == (char *) NULL) {
        usage (av[0]);
        return 1;
    }

    szeit = CCutil_zeit ();
    printf ("Reading files ... "); fflush (stdout);
    if (getgraph (edgefilename, datfilename, &ncount, &ecount, &elist, &elen,
                  &xcoord, &ycoord))
        goto CLEANUP;

    if (edgefilename && double_len) {
        if (CCutil_getedges_double (&ncount, edgefilename, &ecount, &elist,
                             &dlen, binary_in))
            goto CLEANUP;
    }
    printf ("DONE\n"); fflush (stdout);

    printf ("Writing files ... "); fflush (stdout);
    if (edgefilename) {
        if (!double_len) {
            if (dumpedges (ncount, ecount, elist, elen))
                goto CLEANUP;
        } else {
            char buf[256];
            sprintf (buf, "out.%d.%d", ncount, ecount);
            if (CCutil_writeedges_double (ncount, buf, ecount, elist,
                                   dlen, binary_out))

/*
            if (CCutil_writeedges_double (ncount, "edges.out", ecount, elist,
                                   dlen, binary_out))
*/
                goto CLEANUP;
        }
    }
    if (datfilename) {
        if (dumpdat (ncount, xcoord, ycoord))
            goto CLEANUP;
    }
    printf ("DONE\n"); fflush (stdout);

    printf ("Total Running Time: %.2f (seconds)\n", CCutil_zeit () - szeit);
    fflush (stdout);

CLEANUP:
    if (xcoord) {
        CC_FREE (xcoord, double);
    }
    if (ycoord) {
        CC_FREE (ycoord, double);
    }
    if (elist) {
        CC_FREE (elist, int);
    }
    if (elen) {
        CC_FREE (elen, int);
    }

    if (dlen) {
        CC_FREE (dlen, double);
    }

    return 0;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "abdie:n:r", &boptind, &boptarg)) != EOF)
        switch (c) {
        case 'a':
            binary_out = 0;
            break;
        case 'b':
            binary_in = 1;
            break;
        case 'd':
            double_len = 1;
            break;
        case 'i':
            integerdata_in = 1;
            break;
        case 'e':
            edgefilename = boptarg;
            break;
        case 'n':
            datfilename = boptarg;
            break;
        case 'r':
            reversed_int = 1;
            break;
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    if (boptind != ac) {
        usage (av[0]);
        return 1;
    }
    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s [-abir] [-e edgefile] [-n: datfile]\n", f);
    fprintf (stderr, "     a: produce ascii output files\n");
    fprintf (stderr, "     b: input files are in binary format\n");
    fprintf (stderr, "     i: input datfile has integer data\n");
    fprintf (stderr, "     e file: edgefile input\n");
    fprintf (stderr, "     d: the edgelengths are doubles\n");
    fprintf (stderr, "     n file: datafile input\n");
    fprintf (stderr, "     r: input files in our old (reversed) format\n");
}

static int getgraph (char *edgefile, char *datfile, int *ncount, int *ecount,
                     int **elist, int **elen, double **xcoord, double **ycoord)
{
    int xi, yi;

    *elist = (int *) NULL;
    *elen = (int *) NULL;
    *xcoord = (double *) NULL;
    *ycoord = (double *) NULL;

    if (edgefile != (char *) NULL && !double_len) {
        if (binary_in) {
            CC_SFILE *f = CCutil_sopen (edgefile, "r");
            int i, k;

            if (f == (CC_SFILE *) NULL)
                return 1;

            if (CCutil_sread_int (f, ncount)) {
                CCutil_sclose (f);
                return 1;
            }
            if (CCutil_sread_int (f, ecount)) {
                CCutil_sclose (f);
                return 1;
            }

            *elist = CC_SAFE_MALLOC(2 * (*ecount), int);
            if (!(*elist)) {
                CCutil_sclose (f);
                return 1;
            }
            *elen = CC_SAFE_MALLOC(*ecount, int);
            if (!(*elen)) {
                CCutil_sclose (f);
                CC_FREE (*elist, int);
                return 1;
            }
            for (i = 0, k = 0; i < *ecount; i++) {
                if (CCutil_sread_int (f, &((*elist)[k++]))) {
                    CCutil_sclose (f);
                    return 1;
                }
                if (CCutil_sread_int (f, &((*elist)[k++]))) {
                    CCutil_sclose (f);
                    return 1;
                }
                if (CCutil_sread_int (f, &((*elen)[i]))) {
                    CCutil_sclose (f);
                    return 1;
                }
            }
            if (CCutil_sclose (f))
                return 1;
        } else {
            FILE *in = fopen (edgefile, "r");
            int i, k;
            if (in == (FILE *) NULL) {
                perror (edgefile);
                fprintf (stderr, "Unable to open %s for input\n", edgefile);
                return 1;
            }
            *ncount = CCutil_readint (in);
            *ecount = CCutil_readint (in);

            *elist = CC_SAFE_MALLOC(2 * (*ecount), int);
            if (!(*elist)) {
                fclose (in);
                return 1;
            }
            *elen = CC_SAFE_MALLOC(*ecount, int);
            if (!(*elen)) {
                fclose (in);
                CC_FREE (*elist, int);
                return 1;
            }

            for (i = 0, k = 0; i < *ecount; i++) {
                (*elist)[k++] = CCutil_readint (in);
                (*elist)[k++] = CCutil_readint (in);
                (*elen)[i] = CCutil_readint (in);
            }
            fclose (in);
        }
    }

    if (datfile != (char *) NULL) {
        if (binary_in) {
            CC_SFILE *f = CCutil_sopen (datfile, "r");
            int i;

            if (f == (CC_SFILE *) NULL)
                return 1;
            if (reversed_int) {
                if (CCutil_sread_int_r (f, &i)) {
                    CCutil_sclose (f);
                    return 1;
                }
            } else {
                if (CCutil_sread_int (f, &i)) {
                    CCutil_sclose (f);
                    return 1;
                }
            }
            if (edgefile && i != *ncount) {
                fprintf (stderr, "dat file does not match edge file\n");
                if (CCutil_sclose (f))
                    fprintf (stderr, "could not close file\n");
                if (edgefile != (char *) NULL) {
                    CC_FREE (*elist, int);
                    CC_FREE (*elen, int);
                }
                return 1;
            } else {
                *ncount = i;
            }
            printf ("%d nodes ... ", *ncount); fflush (stdout);
            *xcoord = CC_SAFE_MALLOC (*ncount, double);
            if (!(*xcoord)) {
                if (CCutil_sclose (f))
                    fprintf (stderr, "could not close file\n");
                if (edgefile != (char *) NULL) {
                    CC_FREE (*elist, int);
                    CC_FREE (*elen, int);
                }
                return 1;
            }
            *ycoord = CC_SAFE_MALLOC (*ncount, double);
            if (!(*ycoord)) {
                if (CCutil_sclose (f))
                    fprintf (stderr, "could not close file\n");
                CC_FREE(*xcoord, double);
                if (edgefile != (char *) NULL) {
                    CC_FREE (*elist, int);
                    CC_FREE (*elen, int);
                }
                return 1;
            }
            if (reversed_int) {
                for (i = 0; i < *ncount; i++) {
                    if (CCutil_sread_int_r (f, &xi)) {
                        CCutil_sclose (f);
                        return 1;
                    }
                    if (CCutil_sread_int_r (f, &yi)) {
                        CCutil_sclose (f);
                        return 1;
                    }
                    (*xcoord)[i] = (double) xi;
                    (*ycoord)[i] = (double) yi;
                }
            } else {
                for (i = 0; i < *ncount; i++) {
                    if (CCutil_sread_int (f, &xi)) {
                        CCutil_sclose (f);
                        return 1;
                    }
                    if (CCutil_sread_int (f, &yi)) {
                        CCutil_sclose (f);
                        return 1;
                    }
                    (*xcoord)[i] = (double) xi;
                    (*ycoord)[i] = (double) yi;
                }
            }
            if (CCutil_sclose (f))
                return 1;
        } else {
            FILE *indat = fopen (datfile, "r");
            int i;
            if (indat == (FILE *) NULL) {
                perror (datfile);
                fprintf (stderr, "Unable to open %s\n", datfile);
                return 1;
            }
            fscanf (indat, " %d", &i);
            if (edgefile && i != *ncount) {
                fprintf (stderr, "dat file does not match edge file\n");
                fclose (indat);
                if (edgefile != (char *) NULL) {
                    CC_FREE (*elist, int);
                    CC_FREE (*elen, int);
                }
                return 1;
            } else {
                *ncount = i;
            }
            *xcoord = CC_SAFE_MALLOC (*ncount, double);
            if (!(*xcoord)) {
                fclose (indat);
                if (edgefile != (char *) NULL) {
                    CC_FREE (*elist, int);
                    CC_FREE (*elen, int);
                }
                return 1;
            }
            *ycoord = CC_SAFE_MALLOC (*ncount, double);
            if (!(*ycoord)) {
                fclose (indat);
                CC_FREE(*xcoord, double);
                if (edgefile != (char *) NULL) {
                    CC_FREE (*elist, int);
                    CC_FREE (*elen, int);
                }
                return 1;
            }
            if (integerdata_in) {
                for (i = 0; i < *ncount; i++) {
                    (*xcoord)[i] = (double) CCutil_readint (indat);
                    (*ycoord)[i] = (double) CCutil_readint (indat);
                }
            } else {
                printf ("WARNING: Doubles will be converted to ints\n");
                fflush (stdout);
                for (i = 0; i < *ncount; i++)
                    fscanf (indat, "%lf %lf", &((*xcoord)[i]),
                            &((*ycoord)[i]));
            }
            fclose (indat);
        }
    }
    return 0;
}

static int dumpedges (int ncount, int ecount, int *elist, int *elen)
{
    int i;

    if (binary_out) {
        CC_SFILE *f = CCutil_sopen ("edge.out", "w");

        if (f == (CC_SFILE *) NULL)
            return 1;
        if (CCutil_swrite_int (f, ncount)) {
            CCutil_sclose (f);
            return 1;
        }
        if (CCutil_swrite_int (f, ecount)) {
            CCutil_sclose (f);
            return 1;
        }

        for (i = 0; i < ecount; i++) {
            if (CCutil_swrite_int (f, elist[2 * i])) {
                CCutil_sclose (f);
                return 1;
            }
            if (CCutil_swrite_int (f, elist[(2 * i) + 1])) {
                CCutil_sclose (f);
                return 1;
            }
            if (CCutil_swrite_int (f, elen[i])) {
                CCutil_sclose (f);
                return 1;
            }
        }
        if (CCutil_sclose (f))
            return 1;
    } else {
        FILE *out = fopen ("edge.out", "w");
        if (out == (FILE *) NULL) {
            perror ("edge.out");
            fprintf (stderr, "Unable to open edge.out for output\n");
            return 1;
        }
        fprintf (out, "%d %d\n", ncount, ecount);
        for (i = 0; i < ecount; i++)
            fprintf (out, "%d %d %d\n", elist[2 * i], elist[(2 * i) + 1],
                     elen[i]);
        fclose (out);
    }

    return 0;
}

static int dumpdat (int ncount, double *xcoord, double *ycoord)
{
    int i;

    if (binary_out) {
        CC_SFILE *f = CCutil_sopen ("dat.out", "w");

        if (f == (CC_SFILE *) NULL)
            return 1;
        if (CCutil_swrite_int (f, ncount)) {
            CCutil_sclose (f);
            return 1;
        }

        for (i = 0; i < ncount; i++) {
            if (CCutil_swrite_int (f, (int) xcoord[i])) {
                CCutil_sclose (f);
                return 1;
            }
            if (CCutil_swrite_int (f, (int) ycoord[i])) {
                CCutil_sclose (f);
                return 1;
            }
        }
        if (CCutil_sclose (f))
            return 1;
    } else {
        FILE *out = fopen ("dat.out", "w");
        if (out == (FILE *) NULL) {
            perror ("dat.out");
            fprintf (stderr, "Unable to open dat.out for output\n");
            return 1;
        }
        fprintf (out, "%d\n", ncount);
        for (i = 0; i < ncount; i++)
            fprintf (out, "%f %f\n", xcoord[i], ycoord[i]);
        fclose (out);
    }
    return 0;
}


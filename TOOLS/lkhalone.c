#include <stdio.h>
#include "machdefs.h"
#include "util.h"

int main (int ac, char **av)
{
    int id, rval = 0;
    double newcost = -1;
    char *name, *tb;
    char buf[1024], buf2[1024], key[1024];
    FILE *p = (FILE *) NULL;
    FILE *out = (FILE *) NULL;

    if (ac < 3) {
        printf ("usage %s: name id\n", *av);
        rval = 1;  goto CLEANUP;
    }

    name = av[1];
    id = atoi(av[2]);

    printf ("LKH Subproblem: %s %d\n", name, id);
    fflush (stdout);

    sprintf (buf, "%s_par.%d", name, id);
    out = fopen (buf, "w");
    if (!out) {
        fprintf (stderr, "could not open %s for output\n", buf);
        rval = 1; goto CLEANUP;
    }

    fprintf (out, "PROBLEM_FILE = %s_lkh_tsp.%d\n", name, id);
    fprintf (out, "MAX_TRIALS = %d\n", 100);
    fprintf (out, "TOUR_FILE = %s_tour.%d\n", name, id);
    fprintf (out, "RUNS = 1\n");
    fprintf (out, "TRACE_LEVEL = 0\n");
    fclose (out);  out = (FILE *) NULL;
    

    sprintf (buf, "lkh %s_par.%d", name, id);
    p = popen (buf, "r");
    if (!p) {
        perror (buf);
        fprintf (stderr, "popen failed\n");
        rval = 1; goto CLEANUP;
    }

    while ((fgets (buf2, sizeof (buf2), p)) != NULL) {
        buf2[sizeof (buf2) - 1] = '\0';
        fputs (buf2, stdout);
        fflush (stdout);
    
        if (sscanf (buf2, "%s", key) != EOF) {
            if (!strcmp (key, "Z")) {
                tb = buf2;
                tb += strlen (key);
                while (*tb == ' ') tb++;
                if (sscanf (tb, "%lf", &newcost) == EOF) {
                    fprintf (stderr, "Could not read tourlen\n");
                    rval = 1;  goto CLEANUP;
                }
            }
        }
    }

    pclose (p);

    if (newcost == -1) {
       fprintf (stderr, "failed to produce a tourlen\n");
       rval = 1; goto CLEANUP;
    }

    printf ("New Tour Cost: %0.0f\n", newcost); fflush (stdout);

CLEANUP:

    if (out) fclose (out);
    return rval;
}

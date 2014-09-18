
/* global variables: xx, yy, nodes, tour_big and drin
    drin gives the thread/jog this node is currently assigned
    to, 0 if it is not assigned to any job)
    small_nr gives the number of a node in the small_tour it is in.
    this routine finds a rectangle which contains the 'right amount of
    nodes. The nodes in this rectangle are returned in xx_small, y
    and where_small.
    nodes, */


/* patch the small tour (tour_small with nodes_small nodes) into the
   big one (tour_big with nodes nodes). */
int doit( int scount, int *slist, dat, int *tour, int ncount )
{
  char *in_out;
  int i, j;
  int *stour, *fixed, *sfixed;
  int *big_nr;   /* array of length scount */
  int *big_tournr;   /* array of length scount */
  int fixed_count;
  double *x, *y;

  /* set up the tour for the subproblem */
  small_nr = malloc(sizeof(int)*ncount);
  in_out = malloc(ncount);
  for (i=0;i<ncount;i++) 
    in_out[i] = 0;
  for (i=0;i<scount;i++) 
    in_out[slist[i]] = 1;

  /** x, y and stour will be sorted like they appear in tour */
  x = malloc(sizeof(double) * scount);
  y = malloc(sizeof(double) * scount);
  stour = malloc(sizeof(int) * scount);
  
  /* calc x and y, don't use the isolated points */
  /* also calc stour, small_nr, & big_nr */
  for (scount=0,i=0;i<ncount;i++) {
    if (in_out[tour[i]]==1) {
      if (in_out[tour[(i-1+ncount)%ncount]]==1 ||
	  in_out[tour[(i+1)%ncount]]==1) {
	y[scount]=yy[tour[i]];
	x[scount]=xx[tour[i]];
	stour[scount] = scount;
	small_nr[tour[i]]=scount;
	big_nr[scount]=i;
	big_tournr[scount]=tour[i];
	scount++;
      }
      else {
	in_out[tour[i]]=0;
	small_nr[tour[i]]=-1;
      }
    }
  }
  
  /* Found out number of fixed edges. */
  for (fixed_count=0,i=0;i<scount;i++) {
    if ( ((big_nr[i]+1)%ncount) != big_nr[(i+1)%scount] ) {
      fixed_count++;
    }
  }
  /* malloc fixed array. */
  fixed = malloc(sizeof(int) * 2 * fixed_count);
  /* malloc fixed array, this one just indexes with scount . */
  sfixed = malloc(sizeof(int) * scount);
  for (i=0;i<scount;i++) 
    sfixed[i] = 0;
  /* set all fixed edges in array. */
  for (fixed_count=0,i=0;i<scount;i++) {
    if ( ((big_nr[i]+1)%ncount) != big_nr[(i+1)%scount] ) {
      fixed[2*fixed_count]=i;
      fixed[2*fixed_count+1]=(i+1)%scount;      
      fixed_count++;
      sfixed[i]=fixed_count;
      sfixed[(i+1)%scount]=fixed_count;
    }
  }

  /****************** CALL LINKERN **********/

  /* output of linkern is a new tour in an array new_stour */


  /* make the new_tour */  
  new_tour[0] = big_nr[new_stour[i]];
  j = 1; /* this is the counter for the big tour */
  for (i=0;i<scount;i++) {
    int snode = new_stour[i];
    int next_snode = new_stour[(i+1)%scount];
    if (sfixed[snode]==0 || sfixed[snode]!=sfixed[next_snode]) {	
      /* Take single edge from new_stour */
      new_tour[j++] = big_nr[next_snode];
    }
    else {
      int i2;
      /* Take set of edges from the big tour */
      int big_tour_snode = big_tournr[snode];
      if (small_nr[tour[big_tour_snode]] != snode) 
	printf("Something is wrong !!\n");

      i2 = big_tour_snode;
      do { /** loop until we are at a node inside subset again */
	i2 = (i2+1)%ncount;
      } while (in_out[tour[i2]]==0);
      if (small_nr[tour[i2]] == next_snode) {
	/* go forward */
	i2 = big_tour_snode;
	do { /** loop until we are at a node inside subset again */
	  i2 = (i2+1)%ncount;
	  new_tour[j++] = tour[i2];
	} while (in_out[tour[i2]]==0);
	if (small_nr[tour[i2]] != next_snode) {
	  printf("Something is wrong II !!\n");	  
	}
      }
      else {
	/* go backward */
	i2 = big_tour_snode;
	do { /** loop until we are at a node inside subset again */
	  i2 = (i2-1+ncount)%ncount;
	  new_tour[j++] = tour[i2];
	} while (in_out[tour[i2]]==0);
	if (small_nr[tour[i2]] != next_snode) {
	  printf("Something is wrong III !!\n");	  
	}
      }
    }
  }
}



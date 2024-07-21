/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "ppr_ffield.h"

extern int tokenize_string (char *, char ***);
extern void *scalloc (int, int, char *);

#ifdef DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x) 
#endif

/****************************************************************************
*  Function     : Write_Force_Field_PPR
*  Description  : Writes contents of the ffid_pointer into a ffield file 
                  with format of type associated with user-selected 
                  force field (-F ppr)
*  Parameters   : Name of force field file, ffield parameters structure, 
                  forcefield array, and MPI process rank
*  Effects      : Writes a new force field file per evaluation pass, per MPI
                  rank
****************************************************************************/

void Write_Force_Field_PPR (char *ffield_file, ppr_interaction * ffdata, double **ffid_pointer, int rank)
{
  int i, token;
  FILE *fo;

  /* open force field for writing */
  if ((fo = fopen (ffield_file, "w")) == NULL) {
    fprintf (stderr, "ERROR: [%d] Cannot write to %s...\n", rank,ffield_file);
    exit (1);
  }

  token = 0;
  fprintf (fo, "PPR\n\n");
  fprintf (fo, "%3d       ! Number of PPR atom types\n", ffdata->num_atom_types);

  /* print atom type section */
  for (i = 0; i < ffdata->num_atom_types; i++) 
    fprintf (fo, "%s %6.2f\n", ffdata->ppratypes[i].label, ffdata->ppratypes[i].mass);
  /* print pair type section */
  fprintf (fo, "\n%3d       ! Number of PPR pairwise parameters\n", ffdata->num_pair_types);
  for (i = 0; i < ffdata->num_pair_types; i++) {
    fprintf (fo, "%s\t%s\t", ffdata->pprptypes[i].type1, ffdata->pprptypes[i].type2);
    fprintf (fo, "%12.6f", *ffid_pointer[token++]);
    fprintf (fo, "%12.6f", *ffid_pointer[token++]);
    fprintf (fo, "%12.6f", *ffid_pointer[token++]);
    fprintf (fo, "%12.6f", *ffid_pointer[token++]);
    fprintf (fo, "%12.6f", *ffid_pointer[token++]);
    fprintf (fo, "%12.6f", *ffid_pointer[token++]);
    fprintf (fo, "%12.6f\n", *ffid_pointer[token++]);
  }
  fclose (fo);
}

/****************************************************************************
*  Function     : Read_Force_Field_PPR
*  Description  : Reads and parses contents of the user-provided ffield file
                  into the appropriate force field data structure format,
                  associated with user-selected force field (-F ppr)
*  Parameters   : Name of force field file, ffield parameters structure,
                  forcefield array, and MPI process rank
*  Effects      : Reads force field parameters per GA evaluation pass, per MPI 
                  rank
****************************************************************************/

int Read_Force_Field_PPR (char *ffield_file, ppr_interaction * ffdata, double ***ffid_pointer,int rank)
{
  FILE *fi;
  char *s;
  char **tmp;
  int a, token, i, j;
  double val;

  /* open force field file */
  if ((fi = fopen (ffield_file, "r")) == NULL) {
    fprintf (stderr, "ERROR: [%d] opening the force field file! terminating...\n",rank);
    exit (1);
  }

  s = (char *) scalloc (MAX_LINE, sizeof (char), "s3");
  tmp = (char **) scalloc (MAX_TOKENS, sizeof (char *), "tmp3");
  for (i = 0; i < MAX_TOKENS; i++)
    tmp[i] = (char *) scalloc (MAX_TOKEN_LEN, sizeof (char), "tmp3i");

  /* reading first header comment */
  fgets (s, MAX_LINE, fi);
  /* reading empty line */
  fgets (s, MAX_LINE, fi);
  /* line 2 is number of ppr atom types */
  fgets (s, MAX_LINE, fi);
  token = tokenize_string (s, &tmp);
  /* reading the number of ppr atom types*/
  ffdata->num_atom_types = atoi (tmp[0]);
  if (ffdata->num_atom_types < 1) {
    fprintf (stderr, "WARNING: [%d] number of ppr atoms in ffield file is 0!\n",rank);
    exit (1);
  }
  for (i = 0; i < ffdata->num_atom_types; i++)
    fgets (s, MAX_LINE, fi);

  /* reading empty line */
  fgets (s, MAX_LINE, fi);
  /* Number of ppr pair types */
  fgets (s, MAX_LINE, fi);
  token = tokenize_string (s, &tmp);
  /* reading the number of ppr pair types*/
  ffdata->num_pair_types = atoi (tmp[0]);
  if (ffdata->num_pair_types < 1) {
    fprintf (stderr, "WARNING: [%d] number of ppr pair types in ffield file is 0!\n",rank);
    exit (1);
  }

  ffdata->ppratypes = (ppr_atom_types *) scalloc (ffdata->num_atom_types, sizeof (ppr_atom_types), "ppr_atom_types");
  ffdata->pprptypes = (ppr_pair_types *) scalloc (ffdata->num_pair_types, sizeof (ppr_pair_types), "ppr_pair_types");
  *ffid_pointer = (double **) scalloc (ffdata->num_pair_types * 7, sizeof (double *), "ffid_pointer");

  // return to start
  rewind (fi);
  for (i = 0; i < 3; i++)
    fgets (s, MAX_LINE, fi);

  for (i = 0; i < ffdata->num_atom_types; i++) {
    fgets (s, MAX_LINE, fi);
    token = tokenize_string (s, &tmp);
    for (j = 0; j < (int) strlen (tmp[0]); ++j) {
      if (tmp[0][j]!=' ')  
        ffdata->ppratypes[i].label[j] = tmp[0][j];
      else
         break;
    }
    ffdata->ppratypes[i].mass = atof(tmp[1]);
  }

#if (DEBUG)
  int idx;
  if (rank == 0 && debug_level == 1) {
    printf ("%d atom types\n\n", ffdata->num_atom_types);
    for (idx = 0; idx < ffdata->num_atom_types; idx++) 
      printf ("%s %7.3f \n", ffdata->ppratypes[idx].label, ffdata->ppratypes[idx].mass);
    printf("\n");
  }
#endif

  fgets (s, MAX_LINE, fi);
  fgets (s, MAX_LINE, fi);
  a = 0;
  for (i = 0; i < ffdata->num_pair_types; i++) {
    fgets (s, MAX_LINE, fi);
    token = tokenize_string (s, &tmp);
    for (j = 0; j < (int) strlen (tmp[0]); ++j)
      ffdata->pprptypes[i].type1[j] = tmp[0][j];
    for (j = 0; j < (int) strlen (tmp[1]); ++j)
      ffdata->pprptypes[i].type2[j] = tmp[1][j];

    val = (double) atof (tmp[2]);
    ffdata->pprptypes[i].p1 = val;
    (*ffid_pointer)[a++] = &ffdata->pprptypes[i].p1;
    val = (double) atof (tmp[3]);
    ffdata->pprptypes[i].p2 = val;
    (*ffid_pointer)[a++] = &ffdata->pprptypes[i].p2;
    val = (double) atof (tmp[4]);
    ffdata->pprptypes[i].p3 = val;
    (*ffid_pointer)[a++] = &ffdata->pprptypes[i].p3;
    val = (double) atof (tmp[5]);
    ffdata->pprptypes[i].p4 = val;
    (*ffid_pointer)[a++] = &ffdata->pprptypes[i].p4;
    val = (double) atof (tmp[6]);
    ffdata->pprptypes[i].p5 = val;
    (*ffid_pointer)[a++] = &ffdata->pprptypes[i].p5;
    val = (double) atof (tmp[7]);
    ffdata->pprptypes[i].p6 = val;
    (*ffid_pointer)[a++] = &ffdata->pprptypes[i].p6;
    val = (double) atof (tmp[8]);
    ffdata->pprptypes[i].p7 = val;
    (*ffid_pointer)[a++] = &ffdata->pprptypes[i].p7;
  }

  ffdata->total_ppr_parameters=a;

#if (DEBUG)
  if (rank == 0 && debug_level == 1) {
    printf ("%d pair types\n\n", ffdata->num_pair_types);
    for (idx = 0; idx < ffdata->num_pair_types; idx++) {
      printf ("%s %s ", ffdata->pprptypes[idx].type1,ffdata->pprptypes[idx].type2);
      printf ("%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n", ffdata->pprptypes[idx].p1, ffdata->pprptypes[idx].p2, ffdata->pprptypes[idx].p3,ffdata->pprptypes[idx].p4, ffdata->pprptypes[idx].p5, ffdata->pprptypes[idx].p6, ffdata->pprptypes[idx].p7);
    }
  printf("\n");
  }
#endif
  
  if (rank == 0 && debug_level == 1) DEBUG_PRINT (("FFIELD: [%d] Read %d PPR parameters\n", rank,ffdata->total_ppr_parameters));

  // close file
  fclose (fi);

  /* deallocate helper storage */
  for (i = 0; i < MAX_TOKENS; i++)
    free (tmp[i]);
  free (tmp);
  free (s);

  if (rank == 0 && debug_level == 1) DEBUG_PRINT ((">> [%d] Finished reading force field\n",rank));

  return a;
}

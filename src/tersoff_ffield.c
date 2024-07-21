/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Code derived and contributed by Anurag Chaudhry, 2015
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "tersoff_ffield.h"
#include "mpi.h"

#ifdef DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x)
#endif

/****************************************************************************
*  Function     : Write_Force_Field_TERSOFF
*  Description  : Writes contents of the ffid_pointer into a ffield file 
                  with format of type TERSOFF
*  Parameters   : Name of force field file, ffield parameters structure, 
                  forcefield array, and MPI process rank
*  Effects      : Writes a new force field file per evaluation pass, per MPI
                  rank
****************************************************************************/

void Write_Force_Field_TERSOFF (char *ffield_file, tersoff_interaction * ffdata, double **ffid_pointer,
			     int rank)
{
  int i, token;
  FILE *fo;

  /* open force field for writing */
  if ((fo = fopen (ffield_file, "w")) == NULL) {
    fprintf (stderr, "ERROR: Cannot write to %s...\n", ffield_file);
    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  token = 0;
  fprintf (fo, "# TERSOFF force field\n");
  fprintf (fo, "# %3d (%3d): entries (parameters)\n", ffdata->num_entries,ffdata->num_entries*14);

  //TERSOFF entries
  for (i = 0; i < ffdata->num_entries; i++) {
    fprintf (fo, "%-3s ", ffdata->tersoff[i].element1);
    fprintf (fo, "%-3s ", ffdata->tersoff[i].element2);
    fprintf (fo, "%-3s ", ffdata->tersoff[i].element3);
    fprintf (fo, "%9.4f ", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f ", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f ", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f ", *ffid_pointer[token++]);     
    fprintf (fo, "%9.4f ", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f ", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f ", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f ", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f ", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f ", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f ", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f ", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f ", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f\n", *ffid_pointer[token++]);
   }
  fclose (fo);
}

/****************************************************************************
*  Function     : Read_Force_Field_TERSOFF
*  Description  : Reads and parses contents of the user-provided ffield file
                  into the appropriate TERSOFF force field data structure format,
*  Parameters   : Name of force field file, ffield parameters structure,
                  forcefield array, and MPI process rank
*  Effects      : Reads force field parameters per GA evaluation pass, per MPI 
                  rank
****************************************************************************/

int Read_Force_Field_TERSOFF (char *ffield_file, tersoff_interaction * tersoff, double ***ffid_pointer,
			   int rank)
{
  FILE *fi;
  char *s;
  char **strtmp;
  int token, i, j, a;
  double val;
  char file[MAX_LINE];

  /* open force field file */
  sprintf (file, "%s", ffield_file);

  if ((fi = fopen (file, "r")) == NULL) {
    fprintf (stderr, "ERROR: opening the force field file! terminating...\n");
    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  s = (char *) scalloc (MAX_LINE, sizeof (char), "s3");
  strtmp = (char **) scalloc (MAX_TOKENS, sizeof (char *), "tmp3");
  for (i = 0; i < MAX_TOKENS; i++)
    strtmp[i] = (char *) scalloc (MAX_TOKEN_LEN, sizeof (char), "tmp3i");
  char *tmpchar;

  /* Counting number of TERSOFF entries, each with 14 parameters plus interacting atoms */
  while (1) {
    if (!fgets(s, MAX_LINE, fi)) break;
    tmpchar=strchr(s,'#');
    if (tmpchar==NULL) 
      tersoff->num_entries+=1;
  }
  rewind (fi);
 
  if (initial_write) { 
    tersoff->tersoff = (tersoff_parameters *) scalloc (tersoff->num_entries, sizeof (tersoff_parameters), "tersoffparams");
    *ffid_pointer = (double **) scalloc (14*tersoff->num_entries, sizeof (double *), "ffid_pointer");
    tersoff->num_atypes=0;
    tersoff->masses = (double *) scalloc (20, sizeof (double), "tersoff_masses");
    tersoff->atypes = (char **) scalloc (20, sizeof (char *), "tersoff_types");
    for (i = 0; i < 101; i++)
      tersoff->atypes[i] = (char *) scalloc (5, sizeof (char), "tersoff_types_i");
  }

  a=0;
  i=0;
  /* Read in each entry with its 14 parameters plus interacting atoms */
  while (1) {
    if (!fgets(s, MAX_LINE, fi)) break;
    tmpchar=strchr(s,'#');
    if (tmpchar==NULL) { 
      token = tokenize_string(s, &strtmp);
      if (token<17) MPI_Abort (MPI_COMM_WORLD,1);
      for (j = 0; j < (int) strlen (strtmp[0]); ++j)
        tersoff->tersoff[i].element1[j] = strtmp[0][j];
      for (j = 0; j < (int) strlen (strtmp[1]); ++j)
        tersoff->tersoff[i].element2[j] = strtmp[1][j];
      for (j = 0; j < (int) strlen (strtmp[2]); ++j)
        tersoff->tersoff[i].element3[j] = strtmp[2][j];
      if ((strcmp(&tersoff->tersoff[i].element1[j],&tersoff->tersoff[i].element2[j])==0) &&
          (strcmp(&tersoff->tersoff[i].element1[j],&tersoff->tersoff[i].element3[j])==0) && initial_write) {
        strcpy(tersoff->atypes[tersoff->num_atypes],tersoff->tersoff[i].element1);
        if (strcmp(tersoff->atypes[tersoff->num_atypes],"Si")==0)
          tersoff->masses[tersoff->num_atypes]=28.0855;
        else if (strcmp(tersoff->atypes[tersoff->num_atypes],"C")==0)
          tersoff->masses[tersoff->num_atypes]=12.0107;
        else if (strcmp(tersoff->atypes[tersoff->num_atypes],"Ge")==0)
          tersoff->masses[tersoff->num_atypes]=72.64;
        else if (strcmp(tersoff->atypes[tersoff->num_atypes],"O")==0) 
          tersoff->masses[tersoff->num_atypes]=15.9994;
        else if (strcmp(tersoff->atypes[tersoff->num_atypes],"N")==0) 
          tersoff->masses[tersoff->num_atypes]=14.0067;
        else if (strcmp(tersoff->atypes[tersoff->num_atypes],"In")==0) 
          tersoff->masses[tersoff->num_atypes]=114.818;
        else if (strcmp(tersoff->atypes[tersoff->num_atypes],"Ga")==0) 
          tersoff->masses[tersoff->num_atypes]=69.723;
        else if (strcmp(tersoff->atypes[tersoff->num_atypes],"As")==0) 
          tersoff->masses[tersoff->num_atypes]=74.9216;
        else if (strcmp(tersoff->atypes[tersoff->num_atypes],"Ta")==0) 
          tersoff->masses[tersoff->num_atypes]=180.94788;
        else if (strcmp(tersoff->atypes[tersoff->num_atypes],"Cu")==0)
          tersoff->masses[tersoff->num_atypes]=63.546;
        else if (strcmp(tersoff->atypes[tersoff->num_atypes],"Ti")==0)
          tersoff->masses[tersoff->num_atypes]=47.867;
        else if (strcmp(tersoff->atypes[tersoff->num_atypes],"Zr")==0)
          tersoff->masses[tersoff->num_atypes]=91.224; 
        else if (strcmp(tersoff->atypes[tersoff->num_atypes],"Hf")==0)
          tersoff->masses[tersoff->num_atypes]=178.49;
        tersoff->num_atypes+=1;
      }
      val = atof(strtmp[3]);
      tersoff->tersoff[i].m=val;
      (*ffid_pointer)[a++] = &tersoff->tersoff[i].m;
      val = atof(strtmp[4]);
      tersoff->tersoff[i].gamma=val;
      (*ffid_pointer)[a++] = &tersoff->tersoff[i].gamma;
      val = atof(strtmp[5]);
      tersoff->tersoff[i].lambda3=val;
      (*ffid_pointer)[a++] = &tersoff->tersoff[i].lambda3;
      val = atof(strtmp[6]);
      tersoff->tersoff[i].c=val;
      (*ffid_pointer)[a++] = &tersoff->tersoff[i].c;
      val = atof(strtmp[7]);
      tersoff->tersoff[i].d=val;
      (*ffid_pointer)[a++] = &tersoff->tersoff[i].d;
      val = atof(strtmp[8]);
      tersoff->tersoff[i].costheta0=val;
      (*ffid_pointer)[a++] = &tersoff->tersoff[i].costheta0;
      val = atof(strtmp[9]);
      tersoff->tersoff[i].n=val;
      (*ffid_pointer)[a++] = &tersoff->tersoff[i].n;
      val = atof(strtmp[10]);
      tersoff->tersoff[i].beta=val;
      (*ffid_pointer)[a++] = &tersoff->tersoff[i].beta;
      val = atof(strtmp[11]);
      tersoff->tersoff[i].lambda2=val;
      (*ffid_pointer)[a++] = &tersoff->tersoff[i].lambda2;
      val = atof(strtmp[12]);
      tersoff->tersoff[i].B=val;
      (*ffid_pointer)[a++] = &tersoff->tersoff[i].B;
      val = atof(strtmp[13]);
      tersoff->tersoff[i].R=val;
      (*ffid_pointer)[a++] = &tersoff->tersoff[i].R;
      val = atof(strtmp[14]);
      tersoff->tersoff[i].D=val;
      (*ffid_pointer)[a++] = &tersoff->tersoff[i].D;
      val = atof(strtmp[15]);
      tersoff->tersoff[i].lambda1=val;
      (*ffid_pointer)[a++] = &tersoff->tersoff[i].lambda1;
      val = atof(strtmp[16]);
      tersoff->tersoff[i].A=val;
      (*ffid_pointer)[a++] = &tersoff->tersoff[i].A;

#if (DEBUG)
      if (rank == 0 && debug_level == 1) {
        if (tersoff->num_entries > 0)
          printf ("[%d]: %s %s %s %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n", 
            i, tersoff->tersoff[i].element1, tersoff->tersoff[i].element2,
            tersoff->tersoff[i].element3,tersoff->tersoff[i].m,tersoff->tersoff[i].gamma,tersoff->tersoff[i].lambda3,tersoff->tersoff[i].c,
            tersoff->tersoff[i].d,tersoff->tersoff[i].costheta0,tersoff->tersoff[i].n,tersoff->tersoff[i].beta,
            tersoff->tersoff[i].lambda2,tersoff->tersoff[i].B,tersoff->tersoff[i].R,tersoff->tersoff[i].D,tersoff->tersoff[i].lambda1,
            tersoff->tersoff[i].A);
      }
#endif
    i+=1;
    }
  }

  // close file
  fclose (fi);

  /* deallocate helper storage */
  for (i = 0; i < MAX_TOKENS; i++)
    free (strtmp[i]);
  free (strtmp);
  free (s);

  if (rank == 0 && debug_level == 1)
    DEBUG_PRINT ((">> Finished reading force field\n"));

  return a;
}

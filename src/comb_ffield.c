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

#include "comb_ffield.h"
#include "mpi.h"

#ifdef DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x)
#endif

/****************************************************************************
*  Function     : Write_Force_Field_COMB
*  Description  : Writes contents of the ffid_pointer into a ffield file 
                  with format of type COMB
*  Parameters   : Name of force field file, ffield parameters structure, 
                  forcefield array, and MPI process rank
*  Effects      : Writes a new force field file per evaluation pass, per MPI
                  rank
****************************************************************************/

void Write_Force_Field_COMB (char *ffield_file, comb_interaction * ffdata, double **ffid_pointer,
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
  fprintf (fo, "# COMB force field\n");
  fprintf (fo, "# %3d (%3d): entries (parameters)\n", ffdata->num_entries,ffdata->num_entries*46);

  //COMB entries
  for (i = 0; i < ffdata->num_entries; i++) {
    fprintf (fo, "%-3s ", ffdata->comb[i].element1);
    fprintf (fo, "%-3s ", ffdata->comb[i].element2);
    fprintf (fo, "%-3s ", ffdata->comb[i].element3);
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
*  Function     : Read_Force_Field_COMB
*  Description  : Reads and parses contents of the user-provided ffield file
                  into the appropriate COMB force field data structure format,
*  Parameters   : Name of force field file, ffield parameters structure,
                  forcefield array, and MPI process rank
*  Effects      : Reads force field parameters per GA evaluation pass, per MPI 
                  rank
****************************************************************************/

int Read_Force_Field_COMB (char *ffield_file, comb_interaction * comb, double ***ffid_pointer,
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

  /* Counting number of COMB entries, each with 46 parameters plus interacting atoms */
  while (1) {
    if (!fgets(s, MAX_LINE, fi)) break;
    tmpchar=strchr(s,'#');
    if (tmpchar==NULL) 
      comb->num_entries+=1;
  }
  rewind (fi);
 
  if (initial_write) { 
    comb->comb = (comb_parameters *) scalloc (comb->num_entries, sizeof (comb_parameters), "combparams");
    *ffid_pointer = (double **) scalloc (46*comb->num_entries, sizeof (double *), "ffid_pointer");
    comb->num_atypes=0;
    comb->masses = (double *) scalloc (20, sizeof (double), "comb_masses");
    comb->atypes = (char **) scalloc (20, sizeof (char *), "comb_types");
    for (i = 0; i < 101; i++)
      comb->atypes[i] = (char *) scalloc (5, sizeof (char), "comb_types_i");
  }

  a=0;
  i=0;
  /* Read in each entry with its 46 parameters plus interacting atoms */
  while (1) {
    if (!fgets(s, MAX_LINE, fi)) break;
    tmpchar=strchr(s,'#');
    if (tmpchar==NULL) { 
      token = tokenize_string(s, &strtmp);
      if (token<49) MPI_Abort (MPI_COMM_WORLD,1);
      for (j = 0; j < (int) strlen (strtmp[0]); ++j)
        comb->comb[i].element1[j] = strtmp[0][j];
      for (j = 0; j < (int) strlen (strtmp[1]); ++j)
        comb->comb[i].element2[j] = strtmp[1][j];
      for (j = 0; j < (int) strlen (strtmp[2]); ++j)
        comb->comb[i].element3[j] = strtmp[2][j];
      if ((strcmp(&comb->comb[i].element1[j],&comb->comb[i].element2[j])==0) &&
          (strcmp(&comb->comb[i].element1[j],&comb->comb[i].element3[j])==0) && initial_write) {
        strcpy(comb->atypes[comb->num_atypes],comb->comb[i].element1);
        if (strcmp(comb->atypes[comb->num_atypes],"Si")==0)
          comb->masses[comb->num_atypes]=28.0855;
        else if (strcmp(comb->atypes[comb->num_atypes],"O")==0) 
          comb->masses[comb->num_atypes]=15.9994;
        else if (strcmp(comb->atypes[comb->num_atypes],"Cu")==0)
          comb->masses[comb->num_atypes]=63.546;
        else if (strcmp(comb->atypes[comb->num_atypes],"Ti")==0)
          comb->masses[comb->num_atypes]=47.867;
        else if (strcmp(comb->atypes[comb->num_atypes],"Zr")==0)
          comb->masses[comb->num_atypes]=91.224; 
        else if (strcmp(comb->atypes[comb->num_atypes],"Hf")==0)
          comb->masses[comb->num_atypes]=178.49;
        comb->num_atypes+=1;
      }
      val = atof(strtmp[3]);
      comb->comb[i].m=val;
      (*ffid_pointer)[a++] = &comb->comb[i].m;
      val = atof(strtmp[4]);
      comb->comb[i].c=val;
      (*ffid_pointer)[a++] = &comb->comb[i].c;
      val = atof(strtmp[5]);
      comb->comb[i].d=val;
      (*ffid_pointer)[a++] = &comb->comb[i].d;
      val = atof(strtmp[6]);
      comb->comb[i].h=val;
      (*ffid_pointer)[a++] = &comb->comb[i].h;
      val = atof(strtmp[7]);
      comb->comb[i].n=val;
      (*ffid_pointer)[a++] = &comb->comb[i].n;
      val = atof(strtmp[8]);
      comb->comb[i].beta=val;
      (*ffid_pointer)[a++] = &comb->comb[i].beta;
      val = atof(strtmp[9]);
      comb->comb[i].lambda21=val;
      (*ffid_pointer)[a++] = &comb->comb[i].lambda21;
      val = atof(strtmp[10]);
      comb->comb[i].lambda22=val;
      (*ffid_pointer)[a++] = &comb->comb[i].lambda22;
      val = atof(strtmp[11]);
      comb->comb[i].B1=val;
      (*ffid_pointer)[a++] = &comb->comb[i].B1;
      val = atof(strtmp[12]);
      comb->comb[i].B2=val;
      (*ffid_pointer)[a++] = &comb->comb[i].B2;
      val = atof(strtmp[13]);
      comb->comb[i].R=val;
      (*ffid_pointer)[a++] = &comb->comb[i].R;
      val = atof(strtmp[14]);
      comb->comb[i].D=val;
      (*ffid_pointer)[a++] = &comb->comb[i].D;
      val = atof(strtmp[15]);
      comb->comb[i].lambda11=val;
      (*ffid_pointer)[a++] = &comb->comb[i].lambda11;
      val = atof(strtmp[16]);
      comb->comb[i].lambda12=val;
      (*ffid_pointer)[a++] = &comb->comb[i].lambda12;
      val = atof(strtmp[17]);
      comb->comb[i].A1=val;
      (*ffid_pointer)[a++] = &comb->comb[i].A1;
      val = atof(strtmp[18]);
      comb->comb[i].A2=val;
      (*ffid_pointer)[a++] = &comb->comb[i].A2;
      val = atof(strtmp[19]);
      comb->comb[i].K_LP_1=val;
      (*ffid_pointer)[a++] = &comb->comb[i].K_LP_1;
      val = atof(strtmp[20]);
      comb->comb[i].K_LP_3=val;
      (*ffid_pointer)[a++] = &comb->comb[i].K_LP_3;
      val = atof(strtmp[21]);
      comb->comb[i].K_LP_6=val;
      (*ffid_pointer)[a++] = &comb->comb[i].K_LP_6;
      val = atof(strtmp[22]);
      comb->comb[i].A123=val;
      (*ffid_pointer)[a++] = &comb->comb[i].A123;
      val = atof(strtmp[23]);
      comb->comb[i].Aconf=val;
      (*ffid_pointer)[a++] = &comb->comb[i].Aconf;
      val = atof(strtmp[24]);
      comb->comb[i].addrep=val;
      (*ffid_pointer)[a++] = &comb->comb[i].addrep;
      val = atof(strtmp[25]);
      comb->comb[i].R_omiga_b=val;
      (*ffid_pointer)[a++] = &comb->comb[i].R_omiga_b;
      val = atof(strtmp[26]);
      comb->comb[i].R_omiga_c=val;
      (*ffid_pointer)[a++] = &comb->comb[i].R_omiga_c;
      val = atof(strtmp[27]);
      comb->comb[i].R_omiga_d=val;
      (*ffid_pointer)[a++] = &comb->comb[i].R_omiga_d;
      val = atof(strtmp[28]);
      comb->comb[i].R_omiga_a=val;
      (*ffid_pointer)[a++] = &comb->comb[i].R_omiga_a;
      val = atof(strtmp[29]);
      comb->comb[i].QL1=val;
      (*ffid_pointer)[a++] = &comb->comb[i].QL1;
      val = atof(strtmp[30]);
      comb->comb[i].QU1=val;
      (*ffid_pointer)[a++] = &comb->comb[i].QU1;
      val = atof(strtmp[31]);
      comb->comb[i].DL1=val;
      (*ffid_pointer)[a++] = &comb->comb[i].DL1;
      val = atof(strtmp[32]);
      comb->comb[i].DU1=val;
      (*ffid_pointer)[a++] = &comb->comb[i].DU1;
      val = atof(strtmp[33]);
      comb->comb[i].QL2=val;
      (*ffid_pointer)[a++] = &comb->comb[i].QL2;
      val = atof(strtmp[34]);
      comb->comb[i].QU2=val;
      (*ffid_pointer)[a++] = &comb->comb[i].QU2;
      val = atof(strtmp[35]);
      comb->comb[i].DL2=val;
      (*ffid_pointer)[a++] = &comb->comb[i].DL2;
      val = atof(strtmp[36]);
      comb->comb[i].DU2=val;
      (*ffid_pointer)[a++] = &comb->comb[i].DU2;
      val = atof(strtmp[37]);
      comb->comb[i].chi=val;
      (*ffid_pointer)[a++] = &comb->comb[i].chi;
      val = atof(strtmp[38]);
      comb->comb[i].dJ=val;
      (*ffid_pointer)[a++] = &comb->comb[i].dJ;
      val = atof(strtmp[39]);
      comb->comb[i].dK=val;
      (*ffid_pointer)[a++] = &comb->comb[i].dK;
      val = atof(strtmp[40]);
      comb->comb[i].dL=val;
      (*ffid_pointer)[a++] = &comb->comb[i].dL;
      val = atof(strtmp[41]);
      comb->comb[i].dM=val;
      (*ffid_pointer)[a++] = &comb->comb[i].dM;
      val = atof(strtmp[42]);
      comb->comb[i].esm=val;
      (*ffid_pointer)[a++] = &comb->comb[i].esm;
      val = atof(strtmp[43]);
      comb->comb[i].cmn1=val;
      (*ffid_pointer)[a++] = &comb->comb[i].cmn1;
      val = atof(strtmp[44]);
      comb->comb[i].cml1=val;
      (*ffid_pointer)[a++] = &comb->comb[i].cml1;
      val = atof(strtmp[45]);
      comb->comb[i].cmn2=val;
      (*ffid_pointer)[a++] = &comb->comb[i].cmn2;
      val = atof(strtmp[46]);
      comb->comb[i].cml2=val;
      (*ffid_pointer)[a++] = &comb->comb[i].cml2;
      val = atof(strtmp[47]);
      comb->comb[i].coulcut=val;
      (*ffid_pointer)[a++] = &comb->comb[i].coulcut;
      val = atof(strtmp[48]);
      comb->comb[i].hfocor=val;
      (*ffid_pointer)[a++] = &comb->comb[i].hfocor;

#if (DEBUG)
      if (rank == 0 && debug_level == 1) {
        if (comb->num_entries > 0)
          printf ("[%d]: %s %s %s %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n", 
            i, comb->comb[i].element1, comb->comb[i].element2,
            comb->comb[i].element3,comb->comb[i].m,comb->comb[i].c,comb->comb[i].d,comb->comb[i].h,
            comb->comb[i].n,comb->comb[i].beta,comb->comb[i].lambda21,comb->comb[i].lambda22,
            comb->comb[i].B1,comb->comb[i].B2,comb->comb[i].R,comb->comb[i].D,comb->comb[i].lambda11,
            comb->comb[i].lambda12,comb->comb[i].A1,comb->comb[i].A2,comb->comb[i].K_LP_1,comb->comb[i].K_LP_3,
            comb->comb[i].K_LP_6,comb->comb[i].A123,comb->comb[i].Aconf,comb->comb[i].addrep,comb->comb[i].R_omiga_b,
            comb->comb[i].R_omiga_c,comb->comb[i].R_omiga_d,comb->comb[i].R_omiga_a,comb->comb[i].QL1,comb->comb[i].QU1,
            comb->comb[i].DL1,comb->comb[i].DU1,comb->comb[i].QL2,comb->comb[i].QU2,comb->comb[i].DL2,comb->comb[i].DU2,
            comb->comb[i].chi,comb->comb[i].dJ,comb->comb[i].dK,comb->comb[i].dL,comb->comb[i].dM,comb->comb[i].esm,
            comb->comb[i].cmn1,comb->comb[i].cml1,comb->comb[i].cmn2,comb->comb[i].cml2,comb->comb[i].coulcut,comb->comb[i].hfocor);
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

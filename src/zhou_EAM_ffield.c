/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
      Code derived and contributed by Anurag Chaudhry, 2015
 ------------------------------------------------------------------------ */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "zhou_EAM_ffield.h"
#include "mpi.h"

#ifdef DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x)
#endif

/*****************************
 * V : pair potential function
 * rho : density function
 * F : embedding function   
 */
double V (double r, double A, double B, double alpha, double beta, double re, double kappa, double lambda) {
   double pair;     
   if (r < 0.5 ) 
     pair = ( A*exp (-alpha * (0.5/re-1.) ) )  /  (1. + pow (0.5/re-kappa, 20.))
       - ( B*exp (-beta * (0.5/re-1.) ) )  /  (1. + pow (0.5/re-lambda, 20.));
   else 
     pair =  ( A*exp (-alpha * (r/re-1.) ) )  /  (1. + pow (r/re-kappa, 20.))
       - ( B*exp (-beta * (r/re-1.) ) )  /  (1. + pow (r/re-lambda, 20.));
   
   return pair;
}

double rho (double r, double fe, double beta, double lambda, double re) {
   double density;     
   if (r < 0.5)
           density = (fe * exp (-beta * (0.5/re-1.)))  /  (1. + pow (0.5/re-lambda, 20.));
   else
           density = (fe * exp (-beta * (r/re-1.)))  /  (1. + pow (r/re-lambda, 20.));

   return density;
}

double F (double rho_, double rhoe, double rhos, double Fn0, double Fn1, double Fn2, double Fn3, double F0, double F1, double F2, double F3, double Fe, double eta) {
   double embedding_func = 0.0;
   double rhon = .85*rhoe,
          rho0 = 1.15*rhoe;
   if (rho_ < rhon) {
      embedding_func =
            Fn0 * pow (rho_/rhon-1., 0.) +
            Fn1 * pow (rho_/rhon-1., 1.) +
            Fn2 * pow (rho_/rhon-1., 2.) +
            Fn3 * pow (rho_/rhon-1., 3.) ;
   }
   else if (rhon <= rho_ && rho_ < rho0)  {
      embedding_func =
            F0 * pow (rho_/rhoe-1., 0.) +
            F1 * pow (rho_/rhoe-1., 1.) +
            F2 * pow (rho_/rhoe-1., 2.) +
            F3 * pow (rho_/rhoe-1., 3.) ;
   }
   else if (rho0 <= rho_) {
      embedding_func =
            Fe*(1. - eta*log ( rho_/rhos )) * pow (rho_/rhos, eta) ;
   }
   return embedding_func;
}

/****************************************************************************
*  Function     : Write_Force_Field_ZHOU_EAM
*  Description  : Writes contents of the ffid_pointer into a ffield file 
                  with format of type ZHOU_EAM
*  Parameters   : Name of force field file, ffield parameters structure, 
                  forcefield array, and MPI process rank
*  Effects      : Writes a new force field file per evaluation pass, per MPI
                  rank
****************************************************************************/

void Write_Force_Field_ZHOU_EAM (char *ffield_file, zhou_EAM_interaction * ffdata, double **ffid_pointer,
			     int rank)
{
  int i, token;
  FILE *fo, *LAMMPSFile;
  char LAMMPSFilename[MAX_LINE], cmd[MAX_LINE]; 
  sprintf (LAMMPSFilename,"zhou.eam.%d", rank);
  int Nr = 2000;  // # of points at which potential and embedding function are evaluated
  int Nrho = 2000; // # of points at which electron density is evaluated
  int atomic_number = 29 ; // arbitrary value; LAMMPS will not read it.
  double mass;  // atomic mass
  double lattice_constant = 3.615; // arbitrary value; LAMMPS will not read it.
  char lattice_type[] = "FCC"; //arbitrary value; LAMMPS will not read it.
  double rcut, dr ;
  double rhomax ;
  double drho;  //distance between points where the electron density is evaluated.
  double A, B, alpha, beta, re, kappa, lambda, rhos, fe, rhoe, Fn0, Fn1, Fn2,Fn3, F0, F1, F2, F3, Fe, eta;
  int ntypes=1;   // number of element types in the potential. 
  /* open force field 'parameter' for writing */
  if ((fo = fopen (ffield_file, "w")) == NULL) {
    fprintf (stderr, "ERROR: Cannot write to %s...\n", ffield_file);
    MPI_Abort (MPI_COMM_WORLD, 1);
  }
  /* open EAM file for writing */
  if ((LAMMPSFile = fopen (LAMMPSFilename, "w")) == NULL) {
    fprintf (stderr, "ERROR: Cannot write to %s...\n", LAMMPSFilename);
    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  token = 0;
  fprintf (fo, "# ZHOU_EAM force field\n");
  fprintf (fo, "# %3d (%3d): entries (parameters)\n", ffdata->num_entries,ffdata->num_entries*20);

  // ZHOU_EAM parameters
  for (i = 0; i < ffdata->num_entries; i++) {
    fprintf (fo, "%-3s ", ffdata->zhou_EAM[i].element1);
    fprintf (fo, "%-3s ", ffdata->zhou_EAM[i].element2);
    fprintf (fo, "%-3s ", ffdata->zhou_EAM[i].element3);
    re = *ffid_pointer[token];
    rcut = sqrt(5.0)*re;  // ffield cutoff in Angstrom (corresponds to appox. 5th neighbor shell in FCC)
    dr = rcut/(double)(Nr-1); //distance between points where the interatomic potential and embedding function is evaluated.
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);
    fe = *ffid_pointer[token];
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);
    rhoe = *ffid_pointer[token];
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);
    rhos = *ffid_pointer[token];
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);     
    alpha = *ffid_pointer[token];
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);     
    beta = *ffid_pointer[token];
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);
    A = *ffid_pointer[token];
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);
    B = *ffid_pointer[token];
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);
    kappa = *ffid_pointer[token];
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);
    lambda = *ffid_pointer[token];
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);
    Fn0 = *ffid_pointer[token];
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);
    Fn1 = *ffid_pointer[token];
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);
    Fn2 = *ffid_pointer[token];
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);
    Fn3 = *ffid_pointer[token];
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);
    F0 = *ffid_pointer[token];
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);
    F1 = *ffid_pointer[token];
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);
    F2 = *ffid_pointer[token];
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);
    F3 = *ffid_pointer[token];
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);
    eta = *ffid_pointer[token];
    fprintf (fo, "%12.6f ", *ffid_pointer[token++]);
    Fe = *ffid_pointer[token];
    fprintf (fo, "%12.6f\n", *ffid_pointer[token++]);
    //printf ("\n\n\nfe =%12.6f rhoe=%12.6f F3=%12.6f Fe=%12.6f kappa=%12.6f lambda=%12.6f\n\n\n", fe,rhoe,F3,Fe,kappa,lambda);
   }
  fclose (fo);
  //prepare EAM data file
  // determine atomic mass
  if (strcmp(ffdata->zhou_EAM[0].element1,"Ni")==0) {
          mass=58.6934;atomic_number=28;  }
  else if (strcmp(ffdata->zhou_EAM[0].element1,"Mg")==0) {
          mass=24.3050; atomic_number=12;}
  else if (strcmp(ffdata->zhou_EAM[0].element1,"Al")==0) {
          mass=26.9815386; atomic_number=13;}
  else if (strcmp(ffdata->zhou_EAM[0].element1,"Co")==0) {
          mass=58.933195; atomic_number=27;}
  else if (strcmp(ffdata->zhou_EAM[0].element1,"Fe")==0) {
          mass=55.845; atomic_number=26;}
  else if (strcmp(ffdata->zhou_EAM[0].element1,"Ag")==0) {
          mass=107.8682; atomic_number=47;}
  else if (strcmp(ffdata->zhou_EAM[0].element1,"Au")==0) {
          mass=196.966569; atomic_number=79;}
  else if (strcmp(ffdata->zhou_EAM[0].element1,"W")==0) {
          mass=183.84; atomic_number=74;}
  else if (strcmp(ffdata->zhou_EAM[0].element1,"Pt")==0) {
          mass=195.084; atomic_number=78;}
  else if (strcmp(ffdata->zhou_EAM[0].element1,"Pb")==0) {
          mass=207.2; atomic_number=82;}
  else if (strcmp(ffdata->zhou_EAM[0].element1,"Mo")==0) {
          mass=95.96; atomic_number=42;}
  else if (strcmp(ffdata->zhou_EAM[0].element1,"Ru")==0) {
          mass=101.07; atomic_number=44;}
  else if (strcmp(ffdata->zhou_EAM[0].element1,"Pd")==0) {
          mass=106.42; atomic_number=46;}
  else if (strcmp(ffdata->zhou_EAM[0].element1,"Ta")==0) {
          mass=180.94788; atomic_number=73;}
  else if (strcmp(ffdata->zhou_EAM[0].element1,"Cu")==0) {
          mass=63.546; atomic_number=29;}
  else if (strcmp(ffdata->zhou_EAM[0].element1,"Ti")==0) {
          mass=47.867; atomic_number=22;}
  else if (strcmp(ffdata->zhou_EAM[0].element1,"Zr")==0) {
          mass=91.224; atomic_number=40;}
  else if (strcmp(ffdata->zhou_EAM[0].element1,"Cr")==0) {
          mass=51.9961; atomic_number=24;}
  else if (strcmp(ffdata->zhou_EAM[0].element1,"V")==0) {
          mass=50.9415; atomic_number=23;}
  else if (strcmp(ffdata->zhou_EAM[0].element1,"Hf")==0) {
          mass=178.49; atomic_number=72;}
  else {mass=0; printf("Mass undefined in EAM potential\n");}

  rhomax = rho (0.0,fe,beta,lambda,re);
  //rhomax = (fe * exp (-beta * (0.0/re-1.)))  /  (1. + pow (0./re-lambda, 20.));
  if (rhomax < 100.0)
          rhomax = 100.0;
  drho = rhomax/(double)(Nrho-1);
  // header for setfl format
  // assuming writing setfl data for single element
  fprintf (LAMMPSFile, " LAMMPS EAM Potential File in setfl Format \n \\
reference: X. W. Zhou, R. A. Johnson, H. N. G. Wadley PRB 69, 144113 (2006)\n\n \\
%5d %s  \n \\
%5d %-24.16E %d %-24.16E %24.16E\n \\
%5d %-15.5g %-15.5g %s\n",
         ntypes, ffdata->zhou_EAM[0].element1,
         Nrho, drho, Nr, dr, rcut,
         atomic_number, mass, lattice_constant, lattice_type);
  // Embedding function
   for (i = 0; i < Nrho; i++) {
      fprintf (LAMMPSFile, "%24.16E", F ((double)i*drho,rhoe,rhos,Fn0,Fn1,Fn2,Fn3,F0,F1,F2,F3,Fe,eta));
      if ((i+1) % 5 == 0)
              fprintf (LAMMPSFile, "\n");
   }
   // Density function
   for (i = 0; i < Nr; i++) {
      fprintf (LAMMPSFile, "%24.16E", rho ((double)i*dr,fe,beta,lambda,re));
      if ((i+1) % 5 == 0)
              fprintf (LAMMPSFile, "\n");
      }
   // Pair potential
   for (i = 0; i < Nr; i++)  {
      fprintf (LAMMPSFile, "%24.16E", V ((double)i*dr,A,B,alpha,beta,re,kappa,lambda) * (double)i*dr);
      if ((i+1) % 5 == 0)
              fprintf (LAMMPSFile, "\n");
   }

  fclose (LAMMPSFile);
  if (final_write && rank==0) {
    sprintf (cmd, "cp zhou.eam.0 zhou.eam.best");
    system (cmd);
  }
}

/****************************************************************************
*  Function     : Read_Force_Field_ZHOU_EAM
*  Description  : Reads and parses contents of the user-provided ffield file
                  into the appropriate ZHOU_EAM force field data structure format,
*  Parameters   : Name of force field file, ffield parameters structure,
                  forcefield array, and MPI process rank
*  Effects      : Reads force field parameters per GA evaluation pass, per MPI 
                  rank
****************************************************************************/

int Read_Force_Field_ZHOU_EAM (char *ffield_file, zhou_EAM_interaction * zhou_EAM, double ***ffid_pointer,
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

  /* Counting number of ZHOU_EAM entries, each with 20 parameters plus interacting atoms */
  while (1) {
    if (!fgets(s, MAX_LINE, fi)) break;
    tmpchar=strchr(s,'#');
    if (tmpchar==NULL) 
      zhou_EAM->num_entries+=1;
  }
  rewind (fi);
 
  if (initial_write) { 
    zhou_EAM->zhou_EAM = (zhou_EAM_parameters *) scalloc (zhou_EAM->num_entries, sizeof (zhou_EAM_parameters), "zhou_EAMparams");
    *ffid_pointer = (double **) scalloc (20*zhou_EAM->num_entries, sizeof (double *), "ffid_pointer");
    zhou_EAM->num_atypes=0;
    zhou_EAM->masses = (double *) scalloc (20, sizeof (double), "zhou_EAM_masses");
    zhou_EAM->atypes = (char **) scalloc (20, sizeof (char *), "zhou_EAM_types");
    for (i = 0; i < 20; i++)
      zhou_EAM->atypes[i] = (char *) scalloc (5, sizeof (char), "zhou_EAM_types_i");
  }

  a=0;
  i=0;
  /* Read in each entry with its 20 parameters plus interacting atoms */
  while (1) {
    if (!fgets(s, MAX_LINE, fi)) break;
    tmpchar=strchr(s,'#');
    if (tmpchar==NULL) { 
      token = tokenize_string(s, &strtmp);
      if (token<23) MPI_Abort (MPI_COMM_WORLD,1);
      for (j = 0; j < (int) strlen (strtmp[0]); ++j)
        zhou_EAM->zhou_EAM[i].element1[j] = strtmp[0][j];
      for (j = 0; j < (int) strlen (strtmp[1]); ++j)
        zhou_EAM->zhou_EAM[i].element2[j] = strtmp[1][j];
      for (j = 0; j < (int) strlen (strtmp[2]); ++j)
        zhou_EAM->zhou_EAM[i].element3[j] = strtmp[2][j];
      if ((strcmp(&zhou_EAM->zhou_EAM[i].element1[j],&zhou_EAM->zhou_EAM[i].element2[j])==0) &&
          (strcmp(&zhou_EAM->zhou_EAM[i].element1[j],&zhou_EAM->zhou_EAM[i].element3[j])==0) && initial_write) {
        strcpy(zhou_EAM->atypes[zhou_EAM->num_atypes],zhou_EAM->zhou_EAM[i].element1);
        if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"Ni")==0)
          zhou_EAM->masses[zhou_EAM->num_atypes]=58.6934;
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"Mg")==0)
          zhou_EAM->masses[zhou_EAM->num_atypes]=24.3050;
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"Al")==0)
          zhou_EAM->masses[zhou_EAM->num_atypes]=26.9815386;
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"Co")==0)
          zhou_EAM->masses[zhou_EAM->num_atypes]=58.933195;
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"Fe")==0)
          zhou_EAM->masses[zhou_EAM->num_atypes]=55.845;
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"Ag")==0) 
          zhou_EAM->masses[zhou_EAM->num_atypes]=107.8682;
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"Au")==0) 
          zhou_EAM->masses[zhou_EAM->num_atypes]=196.966569;
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"W")==0) 
          zhou_EAM->masses[zhou_EAM->num_atypes]=183.84;
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"Pt")==0) 
          zhou_EAM->masses[zhou_EAM->num_atypes]=195.084;
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"Pb")==0) 
          zhou_EAM->masses[zhou_EAM->num_atypes]=207.2;
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"Mo")==0) 
          zhou_EAM->masses[zhou_EAM->num_atypes]=95.96;
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"Ru")==0) 
          zhou_EAM->masses[zhou_EAM->num_atypes]=101.07;
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"Pd")==0) 
          zhou_EAM->masses[zhou_EAM->num_atypes]=106.42;
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"Ta")==0) 
          zhou_EAM->masses[zhou_EAM->num_atypes]=180.94788;
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"Cu")==0)
          zhou_EAM->masses[zhou_EAM->num_atypes]=63.546;
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"Ti")==0)
          zhou_EAM->masses[zhou_EAM->num_atypes]=47.867;
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"V")==0)
          zhou_EAM->masses[zhou_EAM->num_atypes]=50.9415;
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"Zr")==0)
          zhou_EAM->masses[zhou_EAM->num_atypes]=91.224; 
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"Cr")==0)
          zhou_EAM->masses[zhou_EAM->num_atypes]=51.9961; 
        else if (strcmp(zhou_EAM->atypes[zhou_EAM->num_atypes],"Hf")==0)
          zhou_EAM->masses[zhou_EAM->num_atypes]=178.49;
        zhou_EAM->num_atypes+=1;
      }
      val = atof(strtmp[3]);
      zhou_EAM->zhou_EAM[i].re=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].re;
      val = atof(strtmp[4]);
      zhou_EAM->zhou_EAM[i].fe=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].fe;
      val = atof(strtmp[5]);
      zhou_EAM->zhou_EAM[i].rhoe=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].rhoe;
      val = atof(strtmp[6]);
      zhou_EAM->zhou_EAM[i].rhos=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].rhos;
      val = atof(strtmp[7]);
      zhou_EAM->zhou_EAM[i].alpha=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].alpha;
      val = atof(strtmp[8]);
      zhou_EAM->zhou_EAM[i].beta=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].beta;
      val = atof(strtmp[9]);
      zhou_EAM->zhou_EAM[i].A=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].A;
      val = atof(strtmp[10]);
      zhou_EAM->zhou_EAM[i].B=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].B;
      val = atof(strtmp[11]);
      zhou_EAM->zhou_EAM[i].kappa=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].kappa;
      val = atof(strtmp[12]);
      zhou_EAM->zhou_EAM[i].lambda=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].lambda;
      val = atof(strtmp[13]);
      zhou_EAM->zhou_EAM[i].Fn0=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].Fn0;
      val = atof(strtmp[14]);
      zhou_EAM->zhou_EAM[i].Fn1=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].Fn1;
      val = atof(strtmp[15]);
      zhou_EAM->zhou_EAM[i].Fn2=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].Fn2;
      val = atof(strtmp[16]);
      zhou_EAM->zhou_EAM[i].Fn3=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].Fn3;
      val = atof(strtmp[17]);
      zhou_EAM->zhou_EAM[i].F0=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].F0;
      val = atof(strtmp[18]);
      zhou_EAM->zhou_EAM[i].F1=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].F1;
      val = atof(strtmp[19]);
      zhou_EAM->zhou_EAM[i].F2=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].F2;
      val = atof(strtmp[20]);
      zhou_EAM->zhou_EAM[i].F3=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].F3;
      val = atof(strtmp[21]);
      zhou_EAM->zhou_EAM[i].eta=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].eta;
      val = atof(strtmp[22]);
      zhou_EAM->zhou_EAM[i].Fe=val;
      (*ffid_pointer)[a++] = &zhou_EAM->zhou_EAM[i].Fe;

#if (DEBUG)
      if (rank == 0 && debug_level == 1) {
        if (zhou_EAM->num_entries > 0)
          printf ("[%d]: %s %s %s %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n", 
            i, zhou_EAM->zhou_EAM[i].element1, zhou_EAM->zhou_EAM[i].element2,
            zhou_EAM->zhou_EAM[i].element3,zhou_EAM->zhou_EAM[i].re,zhou_EAM->zhou_EAM[i].fe,zhou_EAM->zhou_EAM[i].rhoe,zhou_EAM->zhou_EAM[i].rhos,
            zhou_EAM->zhou_EAM[i].alpha,zhou_EAM->zhou_EAM[i].beta,zhou_EAM->zhou_EAM[i].A,zhou_EAM->zhou_EAM[i].B,
            zhou_EAM->zhou_EAM[i].kappa,zhou_EAM->zhou_EAM[i].lambda,zhou_EAM->zhou_EAM[i].Fn0,zhou_EAM->zhou_EAM[i].Fn1,zhou_EAM->zhou_EAM[i].Fn2,zhou_EAM->zhou_EAM[i].Fn3,
            zhou_EAM->zhou_EAM[i].F0,zhou_EAM->zhou_EAM[i].F1,zhou_EAM->zhou_EAM[i].F2,zhou_EAM->zhou_EAM[i].F3,zhou_EAM->zhou_EAM[i].eta,
            zhou_EAM->zhou_EAM[i].Fe);
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

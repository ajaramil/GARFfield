/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#ifndef TERSOFF_MOD_TYPES_H
#define TERSOFF_MOD_TYPES_H

#define TERSOFF_MOD_MAX_STR            512
#define TERSOFF_MOD_MAX_ATOM_TYPES     25
#define TERSOFF_MOD_MAX_MOLECULE_SIZE  20
#define MAX_LINE            512
#define MAX_TOKENS          512
#define MAX_TOKEN_LEN       512
#define SUCCESS  1

typedef struct
{
  char element1[5];             // center atom in 3body interaction
  char element2[5];             // the atom bonded to the center atom
  char element3[5];             // the atom influencing the 1-2 bond in a bond-order sense
  double beta;
  double alpha;         // beta, alpha, c3, c2, c1, c5, c4, and h are used only for 3-body interactions
  double h;
  double eta;
  double beta_ters;     //  must be 1.0
  double lambda2; 		//  (1/distance units)
  double B;			    //  (energy units)
  double R;			    // distance units, 
  double D;			    // distance units, 
  double lambda1; 		//  (1/distance units)
  double A;  			//  (energy units)
  double n;             // n, eta, lambda2, B, lambda1 and A are used only for 2-body interactions
  double c1;     		
  double c2;     		
  double c3;     		
  double c4;     		
  double c5;     		
} tersoff_mod_parameters;

typedef struct
{
  int num_entries;
  int num_atypes;
  char **atypes;
  double *masses;
  tersoff_mod_parameters *tersoff_mod;
} tersoff_mod_interaction;

#endif


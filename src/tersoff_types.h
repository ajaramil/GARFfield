/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
      Code derived and contributed by Anurag Chaudhry, 2015
 ------------------------------------------------------------------------ */

#ifndef TERSOFF_TYPES_H
#define TERSOFF_TYPES_H

#define TERSOFF_MAX_STR            512
#define TERSOFF_MAX_ATOM_TYPES     25
#define TERSOFF_MAX_MOLECULE_SIZE  20
#define MAX_LINE            512
#define MAX_TOKENS          512
#define MAX_TOKEN_LEN       512
#define SUCCESS  1

typedef struct
{
  char element1[5];             // center atom in 3body interaction
  char element2[5];             // the atom bonded to the center atom
  char element3[5];             // the atom influencing the 1-2 bond in a bond-order sense
  double m;             //  m must be either 3 or 1
  double gamma;         // m, gamma, lambda3, c, d and costheta0 are used only for 3-body interactions
  double d;
  double c;
  double n;             // n, beta, lambda2, B, lambda1 and A are used only for 2-body interactions
  double beta;
  double lambda1; 		//  (1/distance units)
  double lambda2; 		//  (1/distance units)
  double lambda3; 		// (1/distance units)
  double A;  			//  (energy units)
  double B;			    //  (energy units)
  double R;			    // distance units, 
  double D;			    // distance units, 
  double costheta0; 	// can be a value < -1 or > 1
} tersoff_parameters;

typedef struct
{
  int num_entries;
  int num_atypes;
  char **atypes;
  double *masses;
  tersoff_parameters *tersoff;
} tersoff_interaction;

#endif


/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
     Code derived and contributed by Anurag Chaudhry, 2015
 ------------------------------------------------------------------------ */

#ifndef ZHOU_EAM_TYPES_H
#define ZHOU_EAM_TYPES_H

#define ZHOU_EAM_MAX_STR            512
#define ZHOU_EAM_MAX_ATOM_TYPES     25
#define ZHOU_EAM_MAX_MOLECULE_SIZE  20
#define MAX_LINE            512
#define MAX_TOKENS          512
#define MAX_TOKEN_LEN       512
#define SUCCESS  1

typedef struct
{
  char element1[5];             // center atom in 3body interaction
  char element2[5];             // the atom bonded to the center atom
  char element3[5];             // the atom influencing the 1-2 bond in a bond-order sense
  double re;             //  NN distance
  double fe;         
  double rhoe;
  double rhos;
  double alpha;             
  double beta;
  double A; 		
  double B; 		
  double kappa; 		
  double lambda;  		
  double Fn0;			    
  double Fn1;			     
  double Fn2;			     
  double Fn3;			     
  double F0;			    
  double F1;			     
  double F2;			     
  double F3;			     
  double eta; 	
  double Fe; 	
} zhou_EAM_parameters;

typedef struct
{
  int num_entries;
  int num_atypes;
  char **atypes;
  double *masses;
  zhou_EAM_parameters *zhou_EAM;
} zhou_EAM_interaction;

#endif


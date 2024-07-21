/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012)-2014 California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#ifndef PQEQ_TYPES_H
#define PQEQ_TYPES_H

#define PQEQ_MAX_STR            512
#define PQEQ_MAX_ATOM_TYPES     25
#define PQEQ_MAX_MOLECULE_SIZE  20
#define MAX_LINE            512
#define MAX_TOKENS          512
#define MAX_TOKEN_LEN       512
#define SUCCESS  1

typedef struct
{
//  int type;
  char label[15];
  double mass;
} pqeq_atom_types;

typedef struct
{
  char type1[15];		// E
  char type2[15];		// P
  double p1;			// Xo	
  double p2;			// Jo
  double p3;			// Qc
  double p4;                    // Rc
  double p5;                    // Rs
  double p6;                    // K2
  double p7;                    // K4
} pqeq_pair_types;

typedef struct
{
  int num_atom_types;
  int num_pair_types;
  int total_pqeq_parameters;
  pqeq_atom_types *pqeqatypes;
  pqeq_pair_types *pqeqptypes;
} pqeq_interaction;

#endif

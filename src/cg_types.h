/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#ifndef CG_TYPES_H
#define CG_TYPES_H

#define CG_MAX_STR            512
#define CG_MAX_ATOM_TYPES     25
#define CG_MAX_MOLECULE_SIZE  20
#define MAX_LINE            512
#define MAX_TOKENS          512
#define MAX_TOKEN_LEN       512
#define SUCCESS  1

typedef struct
{
  int type;
  char label[15];
} cg_atom_types;

typedef struct
{
  int type1;
  int type2;
  double p1;			// D	
  double p2;			// alpha
  double p3;			// r0
} cg_pair_types;

typedef struct
{
  int num_atom_types;
  int num_pair_types;
  int total_cg_parameters;
  cg_atom_types *cgatypes;
  cg_pair_types *cgptypes;
} cg_interaction;

#endif

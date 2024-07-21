/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012)-2014 California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#ifndef PPR_TYPES_H
#define PPR_TYPES_H

#define PPR_MAX_STR            512
#define PPR_MAX_ATOM_TYPES     25
#define PPR_MAX_MOLECULE_SIZE  20
#define MAX_LINE            512
#define MAX_TOKENS          512
#define MAX_TOKEN_LEN       512
#define SUCCESS  1

typedef struct
{
//  int type;
  char label[15];
  double mass;
} ppr_atom_types;

typedef struct
{
  char type1[15];
  char type2[15];
  double p1;			// D0
  double p2;			// alpha
  double p3;			// R0
  double p4;                    // beta
  double p5;                    // gamma
  double p6;                    // eta
  double p7;                    // delta
} ppr_pair_types;

typedef struct
{
  int num_atom_types;
  int num_pair_types;
  int total_ppr_parameters;
  ppr_atom_types *ppratypes;
  ppr_pair_types *pprptypes;
} ppr_interaction;

#endif

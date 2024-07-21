/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#ifndef EFF_TYPES_H
#define EFF_TYPES_H

#define EFF_MAX_STR            512
#define EFF_MAX_ATOM_TYPES     25
#define EFF_MAX_MOLECULE_SIZE  20
#define MAX_LINE            512
#define MAX_TOKENS          512
#define MAX_TOKEN_LEN       512
#define SUCCESS  1

typedef struct
{
  int n;
} global_eff_parameters;

typedef struct
{
//  int atype;
  char name[15];		// Two character atom name
  char ecptype[15];		// One character ecp type (s, p, h, x or f)
  double ecpradius;
  double p1;
  double p2;
  double p3;
  double p4;
  double p5;
  double p6;
  double p7;
  double p8;
} ecp_parameters;

typedef struct
{
  int num_atom_types;
  global_eff_parameters gp;
  ecp_parameters *ecp;
} eff_interaction;

#endif

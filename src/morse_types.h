/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012)-2014 California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#ifndef MORSE_TYPES_H
#define MORSE_TYPES_H

#define MORSE_MAX_STR            512
#define MORSE_MAX_ATOM_TYPES     25
#define MORSE_MAX_MOLECULE_SIZE  20
#define MAX_LINE            512
#define MAX_TOKENS          512
#define MAX_TOKEN_LEN       512
#define SUCCESS  1

typedef struct
{
//  int type;
  char label[15];
  double mass;
} morse_atom_types;

typedef struct
{
  char type1[15];
  char type2[15];
  double p1;			// D	
  double p2;			// alpha
  double p3;			// r0
} morse_pair_types;

typedef struct
{
  int num_atom_types;
  int num_pair_types;
  int total_morse_parameters;
  morse_atom_types *morseatypes;
  morse_pair_types *morseptypes;
} morse_interaction;

#endif

/* ----------------------------------------------------------------------
 *    GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.
 *
 *       Copyright (2012-2014) California Institute of Technology
 *          Andres Jaramillo-Botero (ajaramil@caltech.edu)
 *------------------------------------------------------------------------- */

/* Derived from Carl Edward Rasmussen's fmin MATLAB routines
 
  Copyright (C) 2001 and 2002 by Carl Edward Rasmussen. Date 2002-02-13


  (C) Copyright 1999, 2000 & 2001, Carl Edward Rasmussen

  Permission is granted for anyone to copy, use, or modify these
  programs and accompanying documents for purposes of research or
  education, provided this copyright notice is retained, and note is
  made of any changes that have been made.

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "mpi.h"

#define COST_FUNC_DATATYPE_MIN (FLT_MIN*100)

#define RHO 0.01f
#define SIG 0.5f
#define INT 0.1f
#define EXT 3.0f
#define MAX 20
#define RATIO 100.0f

extern void *scalloc (int, int, char *);

int fmincg (void (*costFunc) (double * inputVector, double * cost, double * gradVector, int nDim, MPI_Comm lammps_comm), double * xVector, int nDim, int maxCostFuncCalls, MPI_Comm lammps_comm);

/* ----------------------------------------------------------------------
*    GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.
*
*       Copyright (2012-2014) California Institute of Technology
*          Andres Jaramillo-Botero (ajaramil@caltech.edu)
*------------------------------------------------------------------------- */

/* Derived from Carl Edward Rasmussen's Matlab fmincg code
   Copyright (C) 2001 and 2002 by Carl Edward Rasmussen. Date 2002-02-13
 

  (C) Copyright 1999, 2000 & 2001, Carl Edward Rasmussen
  
  Permission is granted for anyone to copy, use, or modify these
  programs and accompanying documents for purposes of research or
  education, provided this copyright notice is retained, and note is
  made of any changes that have been made.
*/

#include "fmincg.h"

/****************************************************************************
* Function     : fmincg
* Description  : Computes the conjugate gradient of the training set evaluation.
*                The Polack-Ribiere flavour of conjugate gradients is used to 
*                compute search directions, and a line search using quadratic 
*                and cubic polynomial approximations and the Wolfe-Powell stopping 
*                criteria is used together with the slope ratio method
*                for guessing initial step sizes.
* Parameters   : cost function costFunc with parameters vector (Vector), 
*                pointer to cost evaluation value (cost), gradient vector (gradVector)
*                of Vector, number of force field parameters optimized (numParams),
*		 and MPI communicator (lammps_comm), and initial guess parameters
*		 vector (initVector), number of force field paramters (numParams),
*		 maximum number of calls to the training set (maxcostFuncCalls), and MPI
*		 communicator (lammps_comm)
* Effects      : Returns 1 if the maximum number of cost function calls exceeds
*                maxcostFuncCalls and 2 if the line search fails
 *****************************************************************************/

int fmincg (void (*costFunc) (double *Vector,double* cost,double* gradVector,int numParams,MPI_Comm lammps_comm), 
        double *initVector,int numParams,int maxcostFuncCalls,MPI_Comm lammps_comm)
{
  int i, success = 0, costFuncCount = 0, lineSearchFuncCount = 0;

  double lineSearchFailed, f1, d1, z1, f0, f2, d2, f3, d3, z3, limit, z2, A, B, C;
  double *df0, *df1, *df2, *s, *x0, *tmp;

  df0 = (double *) scalloc (numParams, sizeof (double), "df0");
  df1 = (double *) scalloc (numParams, sizeof (double), "df1");
  df2 = (double *) scalloc (numParams, sizeof (double), "df2");
  s = (double *) scalloc (numParams, sizeof (double), "s");
  x0 = (double *) scalloc (numParams, sizeof (double), "x0");
  tmp = (double *) scalloc (numParams, sizeof (double), "tmp");
  double *x = initVector;

  lineSearchFailed = 0;

  if (costFuncCount >= maxcostFuncCalls)
    return 1;
  else
    costFuncCount++;
  (*costFunc) (initVector, &f1, df1, numParams, lammps_comm);

  for (i = 0; i < numParams; i++) {
    s[i] = -df1[i];
  }

  d1 = 0;
  for (i = 0; i < numParams; i++) {
    d1 += -s[i] * s[i];
  }
  z1 = 1.0f / (1 - d1);

  while (1) {
    for (i = 0; i < numParams; i++) {
      x0[i] = x[i];
      df0[i] = df1[i];
    }
    f0 = f1;

    for (i = 0; i < numParams; i++) {
      x[i] = x[i] + (z1) * s[i];
    }
    if (costFuncCount >= maxcostFuncCalls)
      return 1;
    else
      costFuncCount++;
    (*costFunc) (x, &f2, df2, numParams, lammps_comm);

    d2 = 0;
    for (i = 0; i < numParams; i++) {
      d2 += df2[i] * s[i];
    }

    f3 = f1;
    d3 = d1;
    z3 = -z1;

    success = 0;
    limit = -1;
    lineSearchFuncCount = 0;
    // begin line search
    while (1) {
      while ((((f2) > ((f1) + RHO * (z1) * (d1))) || ((d2) > -SIG * (d1)))
	     && lineSearchFuncCount < MAX) {
	limit = z1;
	if ((f2) > (f1)) {
	  z2 = z3 - (0.5f * (d3) * (z3) * (z3)) / ((d3) * (z3) + (f2) - (f3));
	}
	else {
	  A = 6 * ((f2) - (f3)) / (z3) + 3 * ((d2) + (d3));
	  B = 3 * ((f3) - (f2)) - (z3) * ((d3) + 2 * (d2));
	  z2 = (sqrt (B * B - A * (d2) * (z3) * (z3)) - B) / A;
	}
	if (isnan (z2) || isinf (z2)) {
	  z2 = (z3) * 0.5f;
	}
	A = ((z2 < INT * (z3)) ? z2 : INT * (z3));
	B = (1 - INT) * (z3);
	z2 = A > B ? A : B;
	z1 = z1 + z2;

	for (i = 0; i < numParams; i++) {
	  x[i] = x[i] + (z2) * s[i];
	}
	if (costFuncCount >= maxcostFuncCalls)
	  return 1;
	else
	  costFuncCount++;
	lineSearchFuncCount++;
	(*costFunc) (x, &f2, df2, numParams, lammps_comm);

	d2 = 0;
	for (i = 0; i < numParams; i++) {
	  d2 += df2[i] * s[i];
	}
	z3 = z3 - z2;
      }
      if ((f2 > f1 + (z1) * RHO * (d1)) || ((d2) > -SIG * (d1))) {
	break;			//failure
      }
      else if (d2 > SIG * (d1)) {
	success = 1;
	break;
      }
      else if (lineSearchFuncCount >= MAX) {
	break;
      }
      A = 6 * (f2 - f3) / z3 + 3 * (d2 + d3);
      B = 3 * (f3 - f2) - z3 * (d3 + 2 * d2);
      z2 = -d2 * z3 * z3 / (B + sqrt (B * B - A * d2 * z3 * z3));
      if (!(B * B - A * d2 * z3 * z3 >= 0) || isnan (z2) || isinf (z2) || z2 < 0) {
	if (limit < -0.5f) {
	  z2 = z1 * (EXT - 1);
	}
	else {
	  z2 = (limit - z1) / 2;
	}
      }
      else if ((limit > -0.5) && (z2 + z1 > limit)) {
	z2 = (limit - z1) / 2;
      }
      else if ((limit < -0.5) && (z2 + z1 > z1 * EXT)) {
	z2 = z1 * (EXT - 1.0);
      }
      else if (z2 < -z3 * INT) {
	z2 = -z3 * INT;
      }
      else if ((limit > -0.5) && (z2 < (limit - z1) * (1.0 - INT))) {
	z2 = (limit - z1) * (1.0 - INT);
      }
      f3 = f2;
      d3 = d2;
      z3 = -z2;
      z1 = z1 + z2;
      for (i = 0; i < numParams; i++) {
	x[i] = x[i] + z2 * s[i];
      }
      if (costFuncCount >= maxcostFuncCalls)
	return 1;
      else
	costFuncCount++;
      lineSearchFuncCount++;
      (*costFunc) (x, &f2, df2, numParams, lammps_comm);
      d2 = 0;
      for (i = 0; i < numParams; i++) {
	d2 += df2[i] * s[i];
      }
    }
    // line search ended
    if (success) {
      f1 = f2;
      //printf("Cost: %e\n", f1);

      A = 0;
      B = 0;
      C = 0;
      for (i = 0; i < numParams; i++) {
	A += df1[i] * df1[i];
	B += df2[i] * df2[i];
	C += df1[i] * df2[i];
      }
      for (i = 0; i < numParams; i++) {
	s[i] = ((B - C) / A) * s[i] - df2[i];
      }
      for (i = 0; i < numParams; i++) {
	tmp[i] = df1[i];
	df1[i] = df2[i];
	df2[i] = tmp[i];
      }
      d2 = 0;
      for (i = 0; i < numParams; i++) {
	d2 += df1[i] * s[i];
      }
      if (d2 > 0) {
	for (i = 0; i < numParams; i++) {
	  s[i] = -df1[i];
	}
	d2 = 0;
	for (i = 0; i < numParams; i++) {
	  d2 += -s[i] * s[i];
	}
      }
      A = d1 / (d2 - COST_FUNC_DATATYPE_MIN);
      z1 = z1 * ((RATIO < A) ? RATIO : A);
      d1 = d2;
      lineSearchFailed = 0;
    }
    else {
      f1 = f0;
      for (i = 0; i < numParams; i++) {
	x[i] = x0[i];
	df1[i] = df0[i];
      }
      if (lineSearchFailed) {
	break;
      }
      for (i = 0; i < numParams; i++) {
	tmp[i] = df1[i];
	df1[i] = df2[i];
	df2[i] = tmp[i];
      }
      for (i = 0; i < numParams; i++) {
	s[i] = -df1[i];
      }
      d1 = 0;
      for (i = 0; i < numParams; i++) {
	d1 += -s[i] * s[i];
      }
      z1 = 1 / (1 - d1);
      lineSearchFailed = 1;
    }
  }
  printf ("Best error: %4.5f\n", f1);

  free (df0);
  free (df1); 
  free (df2);
  free (s);
  free (x0);
  free (tmp);
  return 2;
}

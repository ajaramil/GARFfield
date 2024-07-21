/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#ifndef COMB_TYPES_H
#define COMB_TYPES_H

#define COMB_MAX_STR            512
#define COMB_MAX_ATOM_TYPES     25
#define COMB_MAX_MOLECULE_SIZE  20
#define MAX_LINE            512
#define MAX_TOKENS          512
#define MAX_TOKEN_LEN       512
#define SUCCESS  1

typedef struct
{
  char element1[5];             // center atom in 3body interaction
  char element2[5];             // the atom bonded to the center atom
  char element3[5];             // the atom influencing the 1-2 bond in a bond-order sense
  double m;
  double c;
  double d;
  double h;
  double n;
  double beta;
  double lambda21; 		// lambda2 of element 1 (1/distance units)
  double lambda22; 		// lambda2 of element 2 (1/distance units)
  double B1;			// B of element 1 (energy units)
  double B2; 			// B of element 2 (energy units)
  double R;			// cutoff, distance units, 0.5*(r_outer + r_inner)
  double D;			// cutoff, distance units, R - r_inner
  double lambda11; 		// lambda1 of element 1 (1/distance units)
  double lambda12; 		// lambda1 of element 2 (1/distance units)
  double A1;  			// A of element 1 (energy units)
  double A2; 			// A of element 2 (energy units)
  double K_LP_1;		// energy units, 1st order Legendre polynomial coefficient
  double K_LP_3; 		// energy units, 3rd order Legendre polynomial coefficient
  double K_LP_6;		// energy units, 6th order Legendre polynomial coefficient
  double A123;			// cos_theta, theta = equilibrium MOM or OMO bond angles
  double Aconf;			// cos_theta, theta = equilibrium MOM or OMO bond-bending coefficient
  double addrep;		// energy units, additional repulsion
  double R_omiga_b; 		// unit-less scaler for B
  double R_omiga_c;		// unit-less scaler for 0.5*(lambda21+lambda22)
  double R_omiga_d; 		// unit-less scaler for 0.5*(lambda11+lambda12)
  double R_omiga_a; 		// unit-less scaler for A
  double QL1;			// charge units, lower charge limit for element 1
  double QU1; 			// charge units, upper charge limit for element 1
  double DL1; 			// distance units, ion radius of element 1 with charge QL1
  double DU1; 			// distance units, ion radius of element 1 with charge QU1
  double QL2; 			// charge units, lower charge limit for element 2
  double QU2; 			// charge units, upper charge limit for element 2
  double DL2; 			// distance units, ion radius of element 2 with charge QL2
  double DU2; 			// distance units, ion radius of element 2 with charge QU2
  double chi; 			// energy units, self energy 1st power term
  double dJ;			// energy units, self energy 2nd power term
  double dK; 			// energy units, self energy 3rd power term
  double dL; 			// energy units, self energy 4th power term
  double dM; 			// energy units, self energy 6th power term
  double esm; 			// distance units, orbital exponent
  double cmn1; 			// self energy penalty, rho 1 of element 1
  double cml1; 			// self energy penalty, rho 1 of element 2
  double cmn2; 			// self energy penalty, rho 2 of element 1
  double cml2; 			// self energy penalty, rho 2 of element 2
  double coulcut; 		// long range Coulombic cutoff, distance units
  double hfocor;		// coordination term
} comb_parameters;

typedef struct
{
  int num_entries;
  int num_atypes;
  char **atypes;
  double *masses;
  comb_parameters *comb;
} comb_interaction;

#endif


/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#ifndef REAXC_TYPES_H
#define REAXC_TYPES_H

#define REAX_MAX_STR            512
#define REAX_MAX_NBRS           6
#define REAX_MAX_3BODY_PARAM    5
#define REAX_MAX_4BODY_PARAM    5
#define REAX_MAX_ATOM_TYPES     25
#define REAX_MAX_MOLECULE_SIZE  20

typedef struct
{
  int n_global;
  double *l;
  int vdw_type;
} global_parameters;


/* Single body parameters (32) */
typedef struct
{
  /* Line one in field file */
  char name[15];		// Two character atom name

  double r_s;			// 1
  double valency;		// 2 Valency of the atom
  double mass;			// 3 Mass of atom
  double r_vdw;			// 4
  double epsilon;		// 5
  double gamma;			// 6
  double r_pi;			// 7
  double valency_e;		// 8

  double nlp_opt;

  /* Line two in field file */
  double alpha;			// 9
  double gamma_w;		// 10
  double valency_boc;		// 11
  double p_ovun5;		// 12
  double p_ovun6;		// 13
  double chi;			// 14
  double eta;			// 15
  double p_hbond;			// 16: 1 for H, 2 for hbonding atoms (O,S,P,N), 0 for others
//  double sbp_16;

  /* Line three in field file */
  double r_pi_pi;		// 17
  double p_lp2;			// 18
  double heat;			// 19
  double p_boc3;		// 20
  double p_boc4;		// 21
  double p_boc5;		// 22
//  double b_o_131;
//  double b_o_132;
//  double b_o_133;
  double sbp_23;		// 23
  double sbp_24;		// 24

  /* Line four in the field file */
  double p_ovun2;		// 25
  double p_val3;		// 26
  double sbp_27;		// 27
  double valency_val;		// 28
  double p_val5;		// 29
  double rcore2;		// 30
  double ecore2;		// 31
  double acore2;		// 32

  /* Line five in force field file - only for lg */
  double lg1;
  double lg2;

} single_body_parameters;



/* Two Body Parameters (16 Bond; 6 off-diag) */
typedef struct
{
  /* type indice */
  int i;
  int j;
  int used;
  int used_offdiag;
  int ffid_ps;
  int od_ffid_ps;		// index of start vale of ffid_pointer

  /* Bond Energy parameters */
  double De_s;			// 1
  double De_p;			// 2
  double De_pp;			// 3
  double p_be1;			// 4
  double p_be2;			// 9

  /* Bond Order parameters */
  double p_bo1;			// 13
  double p_bo2;			// 14
  double p_bo3;			// 10
  double p_bo4;			// 11
  double p_bo5;			// 5
  double p_bo6;			// 7
  double r_s;			// off-diag 4
  double r_p;			// off-diag 5
  double r_pp;			// off-diag 6 r_o distances in BO formula
  double lg_cij;			// lg code

//  double p_boc3;      
//  double p_boc4;
//  double p_boc5;

  /* Over/Under coordination parameters */
  double p_ovun1;		// 8

  /* Van der Waal interaction parameters */
  double D;			// off-diag 1
  double r_vdW;			// off-diag 2
  double alpha;			// off-diag 3
//  double gamma_w;
//  double rcore;
//  double ecore;
//  double acore;

  /* electrostatic parameters */
  double gamma;			// note: this parameter is gamma^-3 and not gamma.

  double v13cor;		// 6
//  double ovc;

  /* unassigned */
  double tbp_12;		// 12
  double tbp_15;		// 15
  double tbp_16;		// 16
} two_body_parameters;



/* 3-body parameters (7)*/
typedef struct
{
  /* valence angle */
  double theta_00;		// 1
  double p_val1;		// 2
  double p_val2;		// 3
  double p_val4;		// 7
  double p_val7;		// 5

  /* penalty */
  double p_pen1;		// 6

  /* 3-body conjugation */
  double p_coa1;		// 4
} three_body_parameters;


typedef struct
{
  /* type index */
  int i;
  int j;
  int k;
  int used;
  int ffid_ps;

  int cnt;
  three_body_parameters prm[REAX_MAX_3BODY_PARAM];
} three_body_header;

/* hydrogen-bond parameters */
typedef struct
{
  /* type indice */
  int i;
  int j;
  int k;
  int l;
  int used;
  int ffid_ps;

  double r0_hb;			// 1
  double p_hb1;			// 2
  double p_hb2;			// 3
  double p_hb3;			// 4
} hbond_parameters;


/* 4-body parameters (5)*/
typedef struct
{
  double V1;			// 1
  double V2;			// 2
  double V3;			// 3

  /* torsion angle */
  double p_tor1;		// 4

  /* 4-body conjugation */
  double p_cot1;		// 5

  /* unassigned */
  double fbp_6;			// 6
  double fbp_7;			// 7

} four_body_parameters;


typedef struct
{
  /* type indice */
  int i;
  int j;
  int k;
  int l;
  int used;
  int ffid_ps;

  int cnt;
  four_body_parameters prm[REAX_MAX_4BODY_PARAM];
} four_body_header;

typedef struct
{
  int num_atom_types;
  int num_bond_types;
  int num_angle_types;
  int num_torsion_types;
  int num_off_diag_types;
  int num_hbond_types;
  global_parameters gp;
  single_body_parameters *sbp;
  two_body_parameters **tbp;
  three_body_header ***thbp;
  hbond_parameters ***hbp;
  four_body_header ****fbp;
} reax_interaction;

#endif

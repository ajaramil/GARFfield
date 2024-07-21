/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#ifndef TRAINSET_H
#define TRAINSET_H

extern int debug_level;
extern int normalizederrors;

#define num_objectives 13

enum {TOTAL,CHARGE,CELL,FREQUENCY,HEATFORM,GEOMETRY,STRUCTURE,ENERGY,FORCE,ATOM_FORCE,STRESS,PRESSURE};

typedef struct
{
  double weight;
  double nweight;
  int n;
  char **sname, line[512];
  char **op;
  double *factor;
  double eng, ff_val;
} eng_tset_item;

typedef struct
{
  int num_atom;
  int atom[4];
  char sname[512], line[512];
  double lit, ff_val;		//literature value
  double weight;
  double nweight;
} geom_tset_item;

typedef struct
{
  char sname[512], line[512];
  double dis, lit, ff_val;         //literature value
  double weight;
  double nweight;
} struc_tset_item;

typedef struct
{
  char sname[512], line[512];
  double weight;
  double nweight;
  int atomid;
  double charge, ff_val;
} charge_tset_item;

typedef struct
{
  char sname[512], line[512];
  double weight;
  double nweight;
  char type[4];
  double lit, ff_val;
} cell_tset_item;

typedef struct
{
  char sname[512], line[512];
  double weight;
  double nweight;
  double ff_val;
} force_tset_item;

typedef struct
{
  char sname[512], line[512];
  int atomid;
  double weight;
  double nweight;
  double fx,fy,fz,ff_val_x, ff_val_y, ff_val_z;
} atom_force_tset_item;

typedef struct
{
  char sname[512], line[512], val[6];
  int id;
  double weight;
  double nweight;
  double lit, ff_val;
} stress_tset_item;

typedef struct
{
  char sname[512], line[512], val[6];
  int id;
  double weight;
  double nweight;
  double lit, ff_val;
} press_tset_item;

struct tset_data
{
  double sumsecweights;	// 1st-level weights (objective functions)
  double *sumentryweights;	// 2nd-level weights (entry in objective function)
  eng_tset_item *eng;
  geom_tset_item *geom;
  struc_tset_item *struc;
  charge_tset_item *charge;
  cell_tset_item *cell;
  force_tset_item *force;
  atom_force_tset_item *atom_force;
  stress_tset_item *stress;
  press_tset_item *press;
};

struct tset_data *tset;
extern int *nfit;
extern double *secweight;
extern double *entryweight;
extern int random_wflag;
void parse_tset (char *, struct tset_data *, int **, int);
int get_num_entries (char *, char *, char *, int, int);
int get_stress_component_id (char **, int);
int get_press_component_id (char *, int);
int calc_initial_geom;
int check_int_or_float (char *, int);

#endif

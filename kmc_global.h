#ifndef KMC_GLOBAL_INCLUDED
#define KMC_GLOBAL_INCLUDED
#include <cstring>
#include "kmc_par.h"

using namespace std;

// Time parameters
//long long int timestep;
//double totaltime;

// System variables
const int nx=  par_nx;
const int ny=  par_ny;
const int nz=  par_nz;
//int nA, nB, nV, nI;

// Migration parameters
const double temp= par_temp; 
const double beta= par_beta;

const double muA= par_muA; 
const double muB= par_muB; 

const double emA= par_emA;
const double emB= par_emB;

// Ising model energy constants
extern double h0;
extern double c1_44, c1_43, c1_42, c1_41, c1_33, c1_32, c1_31, c1_22, c1_21, c1_11;
extern double c2_44, c2_43, c2_42, c2_41, c2_33, c2_32, c2_31, c2_22, c2_21, c2_11;
extern bool is_e2nbr;

// Global functions
extern void error(int nexit, string errinfo, int nnum=0, double num1=0, double num2=0);
extern void error(int nexit, string errinfo, char c[]);
extern double ran_generator();
extern int pbc(int x_, int nx_);

#endif // KMC_GLOBAL_INCLUDED

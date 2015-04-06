#ifndef KMC_GLOBAL_INCLUDED
#define KMC_GLOBAL_INCLUDED
#include <cstring>
#include <vector>
#include "kmc_par.h"

using namespace std;

// Time parameters
extern long long int timestep;
extern double totaltime;

// ltc variables
#define MAX_NNBR 20
extern int n1nbr, n2nbr;	// number of neighbors
extern int v1nbr[MAX_NNBR][3];	// indexes vectors of 1st neighbors
extern int v2nbr[MAX_NNBR][3];	// indexes vectors of 2nd neighbors
extern double vbra[3][3];	// coordinate vectors of bravice lattice

// System variables
const int nx=  par_nx;
const int ny=  par_ny;
const int nz=  par_nz;
extern int nA, nB, nV, nI;

// Files
extern FILE * his_sol;		// history file of solute atoms
extern FILE * his_vcc;		// history file of vacancy and time: record every several steps

// Defect lists
extern vector <int> list_vcc;	// A list contains indexs of all vacancies
extern vector <int> list_int;	// A list contains indexs of all interstitials
extern vector <int> ix, iy, iz; // image box of the vacancy ONE V

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
void write_conf(int *ptr_states);
void write_hissol(const vector<int> (&actions_sol)[2]);
void write_hisvcc();

#endif // KMC_GLOBAL_INCLUDED

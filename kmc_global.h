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

// System
const int nx=  par_nx;
const int ny=  par_ny;
const int nz=  par_nz;
extern int nA, nB, nV, nAA, nBB, nAB;
extern int sum_mag; // sum of magnitization; should be conserved
extern int states[nx][ny][nz];
extern bool itlAB[nx][ny][nz];

// Files
extern FILE * his_sol;		// history file of solute atoms
extern FILE * his_def;		// history file of vacancies and interstitials: record every several steps
extern FILE * out_engy;		// out file for energy calculations

// Solute information for his_sol
extern vector <int> actions_sol[2]; // A list contains solute atom moves from [0] to [1]

// Parameters for mechanisms
const double dis_rec= par_dis_rec;	// recombination distance
const double time_genr= 1.0/(par_dpasm1*nx*ny*nz);

// Defect lists
struct vcc{ // information of an vacancy
	int ltcp;
	int ix, iy, iz;
};
struct itl{ // information of an interstitial; can declare a vector to store all interstitials
	int ltcp;
	int dir; // direction
	int head; // the atom of the itl that in the front along the dir; useful for AB itl
	int ix, iy, iz;
};
extern vector <vcc> list_vcc;	// A list contains information of all vacancies
extern vector <itl> list_itl;  	// A list contains information of all interstitials

// Migration parameters
const double temp= par_temp; 
const double beta= par_beta;

const double muA= par_muA; 
const double muB= par_muB; 

const double emvA= par_emvA;
const double emvB= par_emvB;
const double emiA= par_emiA;
const double emiB= par_emiB;

const double erAA= par_erAA; // rotation energy barrier
const double erAB= par_erAB;
const double erBB= par_erBB;

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
void write_conf();
void write_hissol();
void write_hisdef();

#endif // KMC_GLOBAL_INCLUDED

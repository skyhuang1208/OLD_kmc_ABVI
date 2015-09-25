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
const int x_sink= (int) (nx/2);
extern int nA, nB, nV, nAA, nBB, nAB;
extern int sum_mag; // sum of magnitization; should be conserved
extern int states[nx][ny][nz];
extern bool itlAB[nx][ny][nz];

// Files
extern FILE * his_sol;		// history file of solute atoms
extern FILE * his_def;		// history file of vacancies and interstitials: record every several steps
extern FILE * out_engy;		// out file for energy calculations

// Parameters for mechanisms
const double dis_rec= par_dis_rec;	// recombination distance
const double rate_genr= par_dpasm1*nx*ny*nz;

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

extern int N_genr; // number of F-pair generation 
extern long long int Vja[2], Ija[2];

extern vector <int> list_sink;  // atoms in the sink
extern int nonconsv; // non-conseve magnitization
extern int n_noncsv; // number of vacancies fall into sink before list_sink has value
extern bool is_ncsv; // a flag indicating whether an non-conserved V fallen into sink

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
// Corrected H
extern double c1_44, c1_43, c1_42, c1_41, c1_33, c1_32, c1_31, c1_22, c1_21, c1_11, c1_40, c1_30, c1_20, c1_10;
extern double c1_00_ABA, c1_00_ABB, c1_00_A, c1_00_B, c1_00_V, c1_00_AA, c1_00_AB, c1_00_BB, c1_00_0;
extern double c2_44, c2_43, c2_42, c2_41, c2_33, c2_32, c2_31, c2_22, c2_21, c2_11, c2_40, c2_30, c2_20, c2_10;
extern double c2_00_ABA, c2_00_ABB, c2_00_A, c2_00_B, c2_00_V, c2_00_AA, c2_00_AB, c2_00_BB, c2_00_0;
// Uncorrected H
extern double unc1_44, unc1_43, unc1_42, unc1_41, unc1_33, unc1_32, unc1_31, unc1_22, unc1_21, unc1_11, unc1_40, unc1_30, unc1_20, unc1_10, unc1_00;
extern double unc2_44, unc2_43, unc2_42, unc2_41, unc2_33, unc2_32, unc2_31, unc2_22, unc2_21, unc2_11, unc2_40, unc2_30, unc2_20, unc2_10, unc2_00;

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

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <vector>
#include "kmc_global.h"

using namespace std;

////////// GLOBAL VARIABLES //////////
long long int timestep;
double totaltime;

#define MAX_NNBR 20
int n1nbr, n2nbr;	// number of neighbors
int v1nbr[MAX_NNBR][3];	// indexes vectors of 1st neighbors
int v2nbr[MAX_NNBR][3];	// indexes vectors of 2nd neighbors
double vbra[3][3];	// coordinate vectors of bravice lattice

int nA, nB, nV, nAA, nBB, nAB;
int sum_mag; // sum of magnitization; should be conserved
int states[nx][ny][nz];
bool itlAB[nx][ny][nz]= {false};

FILE * his_sol;		// history file of solute atoms
FILE * his_def;		// history file of defects
FILE * out_engy;	// out file of energy calculations

vector <vcc> list_vcc;	 // A list containing information of all vacancies
vector <itl> list_itl;   // A list containing information of all interstitials

long long int Vja[2]= {0};
long long int Ija[2]= {0};

vector <int> list_sink;  // atoms in the sink
int nonconsv= 0;
int n_noncsv= 0;
bool is_ncsv= false;

double h0;
double c1_44, c1_43, c1_42, c1_41, c1_33, c1_32, c1_31, c1_22, c1_21, c1_11, c1_40, c1_30, c1_20, c1_10, c1_00;
double c1_00_ABA, c1_00_ABB, c1_00_A, c1_00_B, c1_00_V, c1_00_AA, c1_00_AB, c1_00_BB, c1_00_0;
double c2_44, c2_43, c2_42, c2_41, c2_33, c2_32, c2_31, c2_22, c2_21, c2_11, c2_40, c2_30, c2_20, c2_10, c2_00;
double c2_00_ABA, c2_00_ABB, c2_00_A, c2_00_B, c2_00_V, c2_00_AA, c2_00_AB, c2_00_BB, c2_00_0;
double unc1_44, unc1_43, unc1_42, unc1_41, unc1_33, unc1_32, unc1_31, unc1_22, unc1_21, unc1_11, unc1_40, unc1_30, unc1_20, unc1_10, unc1_00;
double unc2_44, unc2_43, unc2_42, unc2_41, unc2_33, unc2_32, unc2_31, unc2_22, unc2_21, unc2_11, unc2_40, unc2_30, unc2_20, unc2_10, unc2_00;
bool is_e2nbr;

////////// GLOBAL FUNCTIONS //////////
void error(int nexit, string errinfo, int nnum, double num1, double num2){
	// exit number represents:
	// 0: main; 1: class_system; 2: class_events
	
	cout << "\nError: ";
	switch(nexit){
		case 0:  cout << "In main function "; break;
		case 1:  cout << "In class_system ";  break;
		case 2:  cout << "In class_events ";  break;
		default: cout << "UNKNOWN NEXIT" << endl;
	}

	cout << errinfo;
	switch(nnum){
		case 0:  cout << endl;				      break;
		case 1:  cout << ": " << num1 << endl;		      break;
		case 2:  cout << ": " << num1 << " " << num2 << endl; break;
		default: cout << "!!!ERROR FUNCTION MALFUNCTION!!! WRONG NNUM!!!" << endl;
	}
	cout << endl;
	exit(1); 
}

void error(int nexit, string errinfo, char c[]){
	// exit number represents:
	// 0: main; 1: class_system; 2: class_events
	
	cout << "\nError: ";
	switch(nexit){
		case 0:  cout << "In main function "; break;
		case 1:  cout << "In class_system ";  break;
		case 2:  cout << "In class_events ";  break;
		default: cout << "UNKNOWN NEXIT" << endl;
	}

	cout << errinfo << " " << c << endl;
	exit(1); 
}

double ran_generator(){
	static bool first= true;
	if(first){
		srand(time(NULL));
		first= false;
	}
	
	return rand()/((double) RAND_MAX+1.0);
}

int pbc(int x_, int nx_){ // Periodic Boundary Condition
	if	(x_<-nx_ || x_>=2*nx_)
		error(1, "(pbc) input x is out of bound", 2, x_, nx_);

	if	(x_<0) 		return (x_ + nx_);
	else if	(x_<nx_)	return  x_;
	else			return (x_ - nx_);
}

void write_conf(){
	ofstream of_xyz;
	ofstream of_ltcp;
	
	// determine the names of conf files
	if(0==timestep){
		of_xyz.open("t0.xyz");
		of_ltcp.open("t0.ltcp");
	}
	else{
		char name_xyz[40], name_ltcp[40];
		sprintf(name_xyz, "%lld", timestep);  strcat(name_xyz, ".xyz");
		sprintf(name_ltcp, "%lld", timestep); strcat(name_ltcp, ".ltcp");
		
		of_xyz.open(name_xyz);
		of_ltcp.open(name_ltcp);
	}

	if(!of_xyz.is_open()) error(1, "(write_conf) xyz file is not opened!");		// check
	if(!of_ltcp.is_open()) error(1, "(write_conf) ltcp file is not opened!");	// check
	
	// write out data
	of_xyz << nx*ny*nz << "\n" << "xyz " << timestep << " ";
	of_xyz << setprecision(15) << totaltime << "\n";
	of_ltcp << nx*ny*nz << "\n" << "ltcp " << timestep << " ";
	of_ltcp << setprecision(15) << totaltime << "\n";
	for(int i=0; i<nx; i ++){
		for(int j=0; j<ny; j ++){
			for(int k=0; k<nz; k ++){
				double x= i*vbra[0][0] + j*vbra[1][0] + k*vbra[2][0];
				double y= i*vbra[0][1] + j*vbra[1][1] + k*vbra[2][1];
				double z= i*vbra[0][2] + j*vbra[1][2] + k*vbra[2][2];
			
				if(-1==states[i][j][k] || 1==states[i][j][k]){
					of_xyz  << states[i][j][k] << " " << x << " " << y << " " << z << endl;
					of_ltcp << states[i][j][k] << " " << i << " " << j << " " << k << endl;
				}
				else if (0==states[i][j][k] && (! itlAB[i][j][k])){
					int id; for(id=0; list_vcc[id].ltcp != i*ny*nz+j*nz+k; id ++);
					
					of_xyz  << states[i][j][k] << " " << x << " " << y << " " << z << " "
						<< list_vcc[id].ix << " " << list_vcc[id].iy << " " << list_vcc[id].iz << endl;
					of_ltcp << states[i][j][k] << " " << i << " " << j << " " << k << " " 
						<< list_vcc[id].ix << " " << list_vcc[id].iy << " " << list_vcc[id].iz << endl;
				}
				else{
					int id; for(id=0; list_itl[id].ltcp != i*ny*nz+j*nz+k; id ++);

					int type= states[i][j][k];
					if(0==type) type= 3;
					of_xyz  << type << " " << x << " " << y << " " << z << " " 
						<< list_itl[id].ix << " " << list_itl[id].iy << " " << list_itl[id].iz << " "
						<< list_itl[id].dir << " " << list_itl[id].head << endl;
					of_ltcp << type << " " << i << " " << j << " " << k << " "
						<< list_itl[id].ix << " " << list_itl[id].iy << " " << list_itl[id].iz << " "
						<< list_itl[id].dir << " " << list_itl[id].head << endl;
				}
	}}}
	
	of_xyz.close();
	of_ltcp.close();
}

void write_hissol(){
	int ncheck= 0;

	fprintf(his_sol, "%d\n", nB);
	fprintf(his_sol, "T: %lld %e\n", timestep, totaltime);
	for(int i=0; i<nx*ny*nz; i++){
		if( -1== *(&states[0][0][0]+i) ){
			ncheck ++;
			fprintf(his_sol, "%d\n", i);
		}
	}

	if(ncheck != nB) error(0, "(write_hissol) nB inconsistent", 2, ncheck, nB); // delete it
}


void write_hisdef(){
	fprintf(his_def, "%lu\n", list_vcc.size()+list_itl.size());
	fprintf(his_def, "T: %lld %e\n", timestep, totaltime);
	for(int i=0; i<list_vcc.size(); i++){
		fprintf(his_def, "0 %d %d %d %d\n", list_vcc[i].ltcp, list_vcc[i].ix, list_vcc[i].iy, list_vcc[i].iz);
	}
	for(int i=0; i<list_itl.size(); i++){
		int type= *(&states[0][0][0]+list_itl[i].ltcp);
		if(0==type) type= 3;
		fprintf(his_def, "%d %d %d %d %d\n", type, list_itl[i].ltcp, list_itl[i].ix, list_itl[i].iy, list_itl[i].iz);
	}
}


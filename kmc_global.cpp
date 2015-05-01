#include <cstdlib>
#include <cstring>
#include <iostream>
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
vector <int> actions_sol[2]; // A list contains solute atom moves from [0] to [1]

vector <vcc> list_vcc;	 // A list contains information of all vacancies
vector <itl> list_itl;   // A list contains information of all interstitials

double h0;
double c1_44, c1_43, c1_42, c1_41, c1_33, c1_32, c1_31, c1_22, c1_21, c1_11;
double c2_44, c2_43, c2_42, c2_41, c2_33, c2_32, c2_31, c2_22, c2_21, c2_11;
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
	of_xyz << nx*ny*nz << "\n" << "xyz " << timestep << " " << totaltime << "\n";
	of_ltcp << nx*ny*nz << "\n" << "ltcp " << timestep << " " << totaltime << "\n";
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
	vector <int> sol_from;
	vector <int> sol_to;

	for(int i=0; i<actions_sol[0].size(); i ++){	// scan all actions of solutes and put them into from and to
		for(int j=0; j<sol_from.size(); j ++){	// see if the same solute atom moves
			if(actions_sol[0].at(i)==sol_to.at(j)){
				sol_to.at(j)=actions_sol[1].at(i);

				if(sol_from.at(j)==sol_to.at(j)){ // delete it if from and to are the same
					sol_from.erase(sol_from.begin()+j);
					sol_to.erase(sol_to.begin()+j);
				}
				
				goto skip_push_back;
			} 
		}
		sol_from.push_back(actions_sol[0].at(i));
		sol_to.push_back  (actions_sol[1].at(i));
skip_push_back:;
	}

	fprintf(his_sol, "%lu\n", sol_from.size());
	fprintf(his_sol, "T: %lld %e\n", timestep, totaltime);
	for(int j=0; j<sol_from.size(); j ++)
		fprintf(his_sol, "%d %d\n", sol_from.at(j), sol_to.at(j));
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


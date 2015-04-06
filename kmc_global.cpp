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

int nA, nB, nV, nI;

FILE * his_sol;		// history file of solute atoms
FILE * his_vcc;		// history file of vacancy and time: record every several steps

vector <int> list_vcc;	 // A list contains indexs of all vacancies
vector <int> list_int;	 // A list contains indexs of all interstitials
vector <int> ix, iy, iz; // image box of the vacancy ONE V

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

void write_conf(int *ptr_states){
	ofstream of_xyz;
	ofstream of_ltcp;
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

	if(!of_xyz.is_open()) error(1, "(write_conf) xyz file is not opened!");
	if(!of_ltcp.is_open()) error(1, "(write_conf) ltcp file is not opened!");
	
	of_xyz << nx*ny*nz << "\n" << "xyz " << timestep << " " << totaltime << "\n";
	of_ltcp << nx*ny*nz << "\n" << "ltcp " << timestep << " " << totaltime << "\n";
	for(int i=0; i<nx; i ++){
		for(int j=0; j<ny; j ++){
			for(int k=0; k<nz; k ++){
				double x= i*vbra[0][0] + j*vbra[1][0] + k*vbra[2][0];
				double y= i*vbra[0][1] + j*vbra[1][1] + k*vbra[2][1];
				double z= i*vbra[0][2] + j*vbra[1][2] + k*vbra[2][2];
				
				of_xyz  << *ptr_states << " " << x << " " << y << " " << z << "\n";
				of_ltcp << *ptr_states << " " << i << " " << j << " " << k << "\n";

				ptr_states ++;
	}}}
	
	of_xyz.close();
	of_ltcp.close();
}

void write_hissol(const vector<int> (&actions_sol)[2]){
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

void write_hisvcc(){
	fprintf(his_vcc, "%lu\n", list_vcc.size());
	fprintf(his_vcc, "T: %lld %e\n", timestep, totaltime);
	for(int i=0; i<list_vcc.size(); i++){
		fprintf(his_vcc, "%d %d %d %d\n", list_vcc.at(i), ix.at(i), iy.at(i), iz.at(i));
	}
}


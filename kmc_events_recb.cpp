#include <iostream>
#include <cmath>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"

using namespace std;

void class_events::rules_recb(int ii, int iv){ // execute the recombination
	int iltcp= list_itl[ii].ltcp;
	int vltcp= list_vcc[iv].ltcp;

	if(*(&states[0][0][0]+vltcp) != 0) error(2, "(rules_recb) the vacancy isnt 0", 1, *(&states[0][0][0]+vltcp));
	int ja; // the jumping atom
	double ran;
	switch(*(&states[0][0][0]+iltcp)){
		case  2:
			nAA --; nV --; nA +=2;
			ja= 1;
			break;
		case  0:
			nAB --; nV --; nA ++; nB ++;
			*(&itlAB[0][0][0]+iltcp)= false;
			ran= ran_generator();
			if(ran<0.5) ja= 1; 
			else	    ja=-1;
			break;
		case -2:
			nBB --; nV --; nB +=2;
			ja=-1;
			break;
		default: error(2, "(rules_recb) an unknown itl type", 1, *(&states[0][0][0]+iltcp));
	}
			
	*(&states[0][0][0]+iltcp) -=ja;
	*(&states[0][0][0]+vltcp) +=ja;
	
	if(-1==*(&states[0][0][0]+iltcp)){
		actions_sol[0].push_back(-1); // -1 meaning appear out of void
		actions_sol[1].push_back(iltcp);
	}
	if(-1==*(&states[0][0][0]+vltcp)){
		actions_sol[0].push_back(-1); // -1 meaning appear out of void
		actions_sol[1].push_back(vltcp);
	}
	
	list_itl.erase(list_itl.begin()+ii);
	list_vcc.erase(list_vcc.begin()+iv);
}

void class_events::recb_dir(int index){
	vector <int> list_rec; // recombination candidates

	int a1= (int) (list_itl[index].ltcp/nz)/ny; // The itl position
	int a2= (int) (list_itl[index].ltcp/nz)%ny;
	int a3= (int)  list_itl[index].ltcp%nz;
	
	int idir= list_itl[index].dir;
	for(int i=-rrecb_nnd; i<=rrecb_nnd; i ++){
		int b1= pbc(a1+i*v1nbr[idir][0], nx);
		int b2= pbc(a2+i*v1nbr[idir][1], ny);
		int b3= pbc(a3+i*v1nbr[idir][2], nz);
		int ltcp= b1*ny*nz + b2*nz + b3;

		if(0==states[b1][b2][b3] && (! itlAB[b1][b2][b3])) list_rec.push_back(ltcp);
	}
	
	if(list_rec.size() != 0){
		double ran= ran_generator();
		int vltcp= list_rec[(int) (ran*list_rec.size())];
	
		for(int i=0; i<list_vcc.size(); i ++){
			if(vltcp==list_vcc[i].ltcp){
				rules_recb(index, i);
				goto find_vcc;
			}
		}
		error(2, "(recb_dir) no vcc found");
find_vcc:;
	}
}

bool class_events::cal_dis(int d1, int d2, int d3){
	double x= d1*vbra[0][0] + d2*vbra[1][0] + d3*vbra[2][0];
	double y= d1*vbra[0][1] + d2*vbra[1][1] + d3*vbra[2][1];
	double z= d1*vbra[0][2] + d2*vbra[1][2] + d3*vbra[2][2];

	double dis2= x*x + y*y + z*z;

	if(dis2 <= rrecb2) return true;
	else		  return false;
}

void class_events::recb_randomV(int index){
	vector <int> list_rec; // recombination candidates

	int a1= (int) (list_vcc[index].ltcp/nz)/ny; // The vcc position
	int a2= (int) (list_vcc[index].ltcp/nz)%ny;
	int a3= (int)  list_vcc[index].ltcp%nz;

	for(int i= 0; i<list_itl.size(); i ++){ // searching in list_itl
		int b1= (int) (list_itl[i].ltcp/nz)/ny;
		int b2= (int) (list_itl[i].ltcp/nz)%ny;
		int b3= (int)  list_itl[i].ltcp%nz;
		
		bool is_inside= cal_dis(a1-b1, a2-b2, a3-b3);

		if(is_inside) list_rec.push_back(i);
	}

	if(list_rec.size() != 0){
		double ran= ran_generator();
		int index2= list_rec[(int) (ran*list_rec.size())];
	
		rules_recb(index2, index);
	}
}

bool class_events::recb_randomI(int index){
	vector <int> list_rec; // recombination candidates

	int a1= (int) (list_itl[index].ltcp/nz)/ny; // The itl position
	int a2= (int) (list_itl[index].ltcp/nz)%ny;
	int a3= (int)  list_itl[index].ltcp%nz;

	for(int i= 0; i<list_vcc.size(); i ++){ // searching in list_vcc
		int b1= (int) (list_vcc[i].ltcp/nz)/ny;
		int b2= (int) (list_vcc[i].ltcp/nz)%ny;
		int b3= (int)  list_vcc[i].ltcp%nz;
		
		bool is_inside= cal_dis(a1-b1, a2-b2, a3-b3);

		if(is_inside) list_rec.push_back(i);
	}

	if(list_rec.size() != 0){
		double ran= ran_generator();
		int index2= list_rec[(int) (ran*list_rec.size())];
	
		rules_recb(index, index2);

		return true;
	}
	else return false;
}



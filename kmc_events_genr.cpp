#include <cstdio>
#include <iostream>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"

using namespace std;

void class_events::genr(){
	double ran;
	int ltcp[2]; // the chosen 2 ltcp in generating frenkel pair; [0]->itl; [1]->vcc

	for(int i=0; i<2; i++){
		do{
			ran= ran_generator();
			ltcp[i]= (int) (ran*nx*ny*nz);
		}while(2==*(&states[0][0][0]+ltcp[i]) || 0==*(&states[0][0][0]+ltcp[i]) || -2==*(&states[0][0][0]+ltcp[i]));
	}
	int iid= list_itl.size();
	int vid= list_vcc.size();

	// initialize the itl in the list_itl
	itl temp_itl;
	list_itl.push_back(temp_itl);
	
	list_itl[iid].ltcp= ltcp[0];
	list_itl[iid].type= *(&states[0][0][0]+ltcp[0]) + *(&states[0][0][0]+ltcp[1]);
	ran= ran_generator();
	list_itl[iid].dir= (int) (ran*n1nbr);
	ran= ran_generator();
	list_itl[iid].head= *(&states[0][0][0]+ltcp[(int) (ran*2)]);
	list_itl[iid].ix= 0; 
	list_itl[iid].iy= 0; 
	list_itl[iid].iz= 0; 

	// initialize the vcc in the list_vcc
	vcc temp_vcc;
	list_vcc.push_back(temp_vcc);
	
	list_vcc[vid].ltcp= ltcp[1];
	list_vcc[vid].ix= 0;
	list_vcc[vid].iy= 0;
	list_vcc[vid].iz= 0;

	// edit numbers
	switch(*(&states[0][0][0]+ltcp[0])){
		case  1: nA --; break;
		case -1: nB --; break;
	}
	switch(*(&states[0][0][0]+ltcp[1])){
		case  1: nA --; break;
		case -1: nB --; break;
	}
	switch(list_itl[iid].type){
		case  2: nAA ++; break;
		case  0: nAB ++; break;
		case -2: nBB ++; break;
	}
	nV ++;
	
	// edit states
	*(&states[0][0][0]+ltcp[0])= *(&states[0][0][0]+ltcp[0]) + *(&states[0][0][0]+ltcp[1]);
	*(&states[0][0][0]+ltcp[1])= 0;

	// recombination check
	recb_randomI(iid);
	if(0==*(&states[0][0][0]+ltcp[1])) recb_randomV(vid); // check in case vid has recombined with iid
}

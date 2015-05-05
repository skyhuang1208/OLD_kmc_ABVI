#include <cstdio>
#include <iostream>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"

using namespace std;

void class_events::genr(){
	double ran;
	int ltcp[2]; // the chosen 2 ltcps in generating frenkel pair; [0]->itl; [1]->vcc

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
	ran= ran_generator();
	list_itl[iid].dir= (int) (ran*n1nbr);
	list_itl[iid].head= *(&states[0][0][0]+ltcp[0]); // choose the ltcp[0] because dir is randomly selected
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

	// Update numbers (before)
	switch(*(&states[0][0][0]+ltcp[0])){
		case  1: nA --; break;
		case -1: nB --; break;
			 // sol hash?
	}
	switch(*(&states[0][0][0]+ltcp[1])){
		case  1: nA --; break;
		case -1: nB --; break;
			 // sol hash?
	}
	
	// Update states
	*(&states[0][0][0]+ltcp[0])= *(&states[0][0][0]+ltcp[0]) + *(&states[0][0][0]+ltcp[1]);
	*(&states[0][0][0]+ltcp[1])= 0;
	if(0==*(&states[0][0][0]+ltcp[0])) *(&itlAB[0][0][0]+ltcp[0])= true;

	// update numbers (after)
	switch(*(&states[0][0][0]+ltcp[0])){
		case  2: nAA ++; break;
		case  0: nAB ++; break;
		case -2: nBB ++; break;
	}
	nV ++;

	// sink check
	int xi= (int) (ltcp[0]/nz)/ny;
	int xv= (int) (ltcp[1]/nz)/ny;
	bool i_in_sink= false;
	if(xi==x_sink){
		sink(false, iid);
		i_in_sink= true;
	}
	if(xv==x_sink)
		sink(true, vid);

	// recombination check
	bool has_recb= false;
	if(! i_in_sink) has_recb= recb_randomI(iid);

	if	(*(&states[0][0][0]+ltcp[1]) != 0) ;	// iid recombines with vid, DO NOT CHECK RECB
	else if	(has_recb) recb_randomV(vid-1);		// iid recombines with other vcc, erase 1 element in list_vcc, so -1
	else 		   recb_randomV(vid);		// iid recombines with nothing, normal check for vid
}

#include <iostream>
#include <cmath>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"

using namespace std;

double class_events::cal_ratesI(vector <bool> &isvcc, vector <double> &rates, vector <int> &ilist, vector <int> &inbrs, vector <int> &jatom){
	double sum_rate= 0;
	if(nAA + nAB + nBB != list_itl.size()) error(2, "(cal_ratesI) itl number inconsistent");
	
	for(int ii=0; ii < list_itl.size(); ii ++){ // ii: index of interstitial
		int ltcp= list_itl[ii].ltcp;
		int i= (int) (ltcp/nz)/ny;
		int j= (int) (ltcp/nz)%ny;
		int k= (int)  ltcp%nz;
		const int stateI= states[i][j][k]; // state of the itl

		if(1==abs(states[i][j][k])) error(2, "(itl_jump) there's an non-interstitial in the itl list");
		
		int dir, opp_dir; // the direction and opposite direction
		dir= list_itl[ii].dir;
		if(dir>=n1nbr/2) opp_dir= dir - n1nbr/2;
		else		 opp_dir= dir + n1nbr/2;
		for(int a=0; a<n1nbr; a ++){
			int x= pbc(i+v1nbr[a][0], nx);
			int y= pbc(j+v1nbr[a][1], ny);
			int z= pbc(k+v1nbr[a][2], nz);
			const int stateA= states[x][y][z]; // state of the atom
			
			if(1==stateA || -1==stateA){
				// determine the type of jumping atom
				int ja; //type of the jumping atom
				if     ( 2==stateI) ja= 1;
				else if(-2==stateI) ja=-1;
				else{
					if(a==opp_dir)  ja= -list_itl[ii].head;
					else if(a==dir) ja=  list_itl[ii].head;
					else{
						double ran= ran_generator();
						if(ran<0.5) ja= 1;
						else	    ja=-1;
					}
				}

				// calculate energy diff
				double e0= cal_energy(i, j, k, x, y, z);

				itl_rules(states[i][j][k], states[x][y][z], ja); // perform an imaginary jump
				if(stateI==states[i][j][k]) error(2, "imaginary jump error"); // delete it 
				//(check if direct using states +/- jatom is faster than the itl_rules function)

				double ediff= cal_energy(i, j, k, x, y, z) - e0;
				
				states[i][j][k]= stateI; //transit back
				states[x][y][z]= stateA;
				
				// give em, mu and calculate rate
				double em, mu;
				if(1==ja) { em= emiA; mu= muA;} 
				else	  { em= emiB; mu= muB;}
				
				if(a==dir || a==opp_dir)
					rates.push_back(mu * exp(-beta*(em+0.5*ediff)));
				else{
					double er; // rotation energy
					if     ( 2==stateI) er= erAA;
					else if(-2==stateI) er= erBB;
					else		    er= erAB;

					rates.push_back(mu * exp(-beta*(em+er+0.5*ediff)));
				}
							
				isvcc.push_back(false);
				ilist.push_back(ii);
				inbrs.push_back(a);
				jatom.push_back(ja);
				
				sum_rate += rates.back();
			}
		}
	}

	return sum_rate;
}

void class_events::itl_rules(int &itl, int &atom, int jatom){ // types of interstitial, atom, jumping atom
	// itl changes to atom, atom changes to interstitial
	if(1 == abs(itl)) error(2, "(itl_rules) input itl isn't an itl", 1, itl);
	if(1 != abs(atom)) error(2, "(itl_rules) input atom isn't an atom", 1, atom);
	if(1 != abs(jatom)) error(2, "(itl_rules) input jatom isn't an atom", 1, jatom);
	int sum_mag= itl + atom;

	if     ( 2==itl){
		switch(atom){
			case  1: itl= 1; atom= 2; break;
			case -1: itl= 1; atom= 0; break;
		}
	}
	else if(-2==itl){
		switch(atom){
			case  1: itl=-1; atom= 0; break;
			case -1: itl=-1; atom=-2; break;
		}
	}
	else if( 0==itl){
		switch(atom){
			case  1:
				if(1==jatom) {itl=-1; atom= 2;}
				else	     {itl= 1; atom= 0;}
				break;
			case -1:
				if(1==jatom) {itl=-1; atom= 0;}
				else	     {itl= 1; atom=-2;}
				break;
		}
	}

	if(sum_mag != itl+atom) error(2, "(itl_rules) magnitization isn't conserved"); // delete it
}

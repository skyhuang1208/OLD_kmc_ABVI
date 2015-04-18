#include <cstdio>
#include <iostream>
#include "kmc_global.h"
#include "kmc_events.h"

using namespace std;

double class_events::cal_energy(int x1, int y1, int z1, int x2, int y2, int z2){ // calculate the energy of the input 2 ltcps 
	int pos[2][3]= {{x1, y1, z1}, {x2, y2, z2}};

	// sums of 1st nn
	int s1um44=0, s1um42=0, s1um22=0;		// class 1
	int s1um43=0, s1um41=0, s1um32=0, s1um21=0;	// class 2
	int s1um33=0, s1um31=0, s1um11=0;		// class 3
	// sums of 2nd nn
	int s2um44=0, s2um42=0, s2um22=0;		// class 1
	int s2um43=0, s2um41=0, s2um32=0, s2um21=0;	// class 2
	int s2um33=0, s2um31=0, s2um11=0;		// class 3
	for(int i=0; i<2; i ++){
		int xi= pbc(pos[i][0], nx);
		int yi= pbc(pos[i][1], ny);
		int zi= pbc(pos[i][2], nz);

		int state_i= states[xi][yi][zi];
		if(0==state_i) continue;

		for(int a=0; a<n1nbr; a ++){ // 1st neighbors
			int xj= pbc(xi+(*(v1nbr+a))[0], nx);
			int yj= pbc(yi+(*(v1nbr+a))[1], ny);
			int zj= pbc(zi+(*(v1nbr+a))[2], nz);

			if((xj==x1) && (yj==y1) && (zj==z1)) continue; // avoid double count the (x1, y1, z1)-(x2, y2, z2) bond

			int state_j1= states[xj][yj][zj];
			if(0==state_j1) continue;

			// calculating the energy for ABVI system. 1st nn
			s1um44 += powc(state_i, 4) * powc(state_j1, 4);
			s1um43 += powc(state_i, 4) * powc(state_j1, 3) + powc(state_i, 3) * powc(state_j1, 4);
			s1um42 += powc(state_i, 4) * powc(state_j1, 2) + powc(state_i, 2) * powc(state_j1, 4);
			s1um41 += powc(state_i, 4) * powc(state_j1, 1) + powc(state_i, 1) * powc(state_j1, 4);
			s1um33 += powc(state_i, 3) * powc(state_j1, 3);
			s1um32 += powc(state_i, 3) * powc(state_j1, 2) + powc(state_i, 2) * powc(state_j1, 3);
			s1um31 += powc(state_i, 3) * powc(state_j1, 1) + powc(state_i, 1) * powc(state_j1, 3);
			s1um22 += powc(state_i, 2) * powc(state_j1, 2);
			s1um21 += powc(state_i, 2) * powc(state_j1, 1) + powc(state_i, 1) * powc(state_j1, 2);
			s1um11 += powc(state_i, 1) * powc(state_j1, 1);
		}
				
		if(! is_e2nbr) continue;

		for(int b=0; b<n2nbr; b ++){ // 2nd neighbors
			int xj= pbc(xi+(*(v2nbr+b))[0], nx);
			int yj= pbc(yi+(*(v2nbr+b))[1], ny);
			int zj= pbc(zi+(*(v2nbr+b))[2], nz);

			int state_j2= states[xj][yj][zj];
			if(0==state_j2) continue;

			// calculating the energy: 2nd nn
			s2um44 += powc(state_i, 4) * powc(state_j2, 4);
			s2um43 += powc(state_i, 4) * powc(state_j2, 3) + powc(state_i, 3) * powc(state_j2, 4);
			s2um42 += powc(state_i, 4) * powc(state_j2, 2) + powc(state_i, 2) * powc(state_j2, 4);
			s2um41 += powc(state_i, 4) * powc(state_j2, 1) + powc(state_i, 1) * powc(state_j2, 4);
			s2um33 += powc(state_i, 3) * powc(state_j2, 3);
			s2um32 += powc(state_i, 3) * powc(state_j2, 2) + powc(state_i, 2) * powc(state_j2, 3);
			s2um31 += powc(state_i, 3) * powc(state_j2, 1) + powc(state_i, 1) * powc(state_j2, 3);
			s2um22 += powc(state_i, 2) * powc(state_j2, 2);
			s2um21 += powc(state_i, 2) * powc(state_j2, 1) + powc(state_i, 1) * powc(state_j2, 2);
			s2um11 += powc(state_i, 1) * powc(state_j2, 1);
		}
	}

	double sum_energy= 
		c1_44 * s1um44 + c1_43 * s1um43 + c1_42 * s1um42 + c1_41 * s1um41 +
		c1_33 * s1um33 + c1_32 * s1um32 + c1_31 * s1um31 +
		c1_22 * s1um22 + c1_21 * s1um21 +
		c1_11 * s1um11 +
		c2_44 * s2um44 + c2_43 * s2um43 + c2_42 * s2um42 + c2_41 * s2um41 +
		c2_33 * s2um33 + c2_32 * s2um32 + c2_31 * s2um31 +
		c2_22 * s2um22 + c2_21 * s2um21 +
		c2_11 * s2um11;
	return sum_energy;
}

double class_events::cal_energy(){
	// sums of 1st nn
	int s1um44=0, s1um42=0, s1um22=0;		// class 1
	int s1um43=0, s1um41=0, s1um32=0, s1um21=0;	// class 2
	int s1um33=0, s1um31=0, s1um11=0;		// class 3
	// sums of 2nd nn
	int s2um44=0, s2um42=0, s2um22=0;		// class 1
	int s2um43=0, s2um41=0, s2um32=0, s2um21=0;	// class 2
	int s2um33=0, s2um31=0, s2um11=0;		// class 3
	for(int i=0; i<nx; i ++){
		for(int j=0; j<ny; j ++){
			for(int k=0; k<nz; k ++){
				int xi= pbc(i, nx);
				int yi= pbc(j, ny);
				int zi= pbc(k, nz);

				int state_i= states[xi][yi][zi];
				if(0==state_i) continue;

				for(int a=0; a<n1nbr; a ++){ // 1st neighbors
					int xj= pbc(xi+(*(v1nbr+a))[0], nx);
					int yj= pbc(yi+(*(v1nbr+a))[1], ny);
					int zj= pbc(zi+(*(v1nbr+a))[2], nz);

					int state_j1= states[xj][yj][zj];
					if(0==state_j1) continue;

					// calculating the energy for ABVI system. 1st nn
					s1um44 += powc(state_i, 4) * powc(state_j1, 4);
					s1um43 += powc(state_i, 4) * powc(state_j1, 3) + powc(state_i, 3) * powc(state_j1, 4);
					s1um42 += powc(state_i, 4) * powc(state_j1, 2) + powc(state_i, 2) * powc(state_j1, 4);
					s1um41 += powc(state_i, 4) * powc(state_j1, 1) + powc(state_i, 1) * powc(state_j1, 4);
					s1um33 += powc(state_i, 3) * powc(state_j1, 3);
					s1um32 += powc(state_i, 3) * powc(state_j1, 2) + powc(state_i, 2) * powc(state_j1, 3);
					s1um31 += powc(state_i, 3) * powc(state_j1, 1) + powc(state_i, 1) * powc(state_j1, 3);
					s1um22 += powc(state_i, 2) * powc(state_j1, 2);
					s1um21 += powc(state_i, 2) * powc(state_j1, 1) + powc(state_i, 1) * powc(state_j1, 2);
					s1um11 += powc(state_i, 1) * powc(state_j1, 1);
				}
				
				if(! is_e2nbr) continue;

				for(int b=0; b<n2nbr; b ++){ // 2nd neighbors
					int xj= pbc(xi+(*(v2nbr+b))[0], nx);
					int yj= pbc(yi+(*(v2nbr+b))[1], ny);
					int zj= pbc(zi+(*(v2nbr+b))[2], nz);

					int state_j2= states[xj][yj][zj];
					if(0==state_j2) continue;

					// calculating the energy: 2nd nn
					s2um44 += powc(state_i, 4) * powc(state_j2, 4);
					s2um43 += powc(state_i, 4) * powc(state_j2, 3) + powc(state_i, 3) * powc(state_j2, 4);
					s2um42 += powc(state_i, 4) * powc(state_j2, 2) + powc(state_i, 2) * powc(state_j2, 4);
					s2um41 += powc(state_i, 4) * powc(state_j2, 1) + powc(state_i, 1) * powc(state_j2, 4);
					s2um33 += powc(state_i, 3) * powc(state_j2, 3);
					s2um32 += powc(state_i, 3) * powc(state_j2, 2) + powc(state_i, 2) * powc(state_j2, 3);
					s2um31 += powc(state_i, 3) * powc(state_j2, 1) + powc(state_i, 1) * powc(state_j2, 3);
					s2um22 += powc(state_i, 2) * powc(state_j2, 2);
					s2um21 += powc(state_i, 2) * powc(state_j2, 1) + powc(state_i, 1) * powc(state_j2, 2);
					s2um11 += powc(state_i, 1) * powc(state_j2, 1);
				}
	}}}

	double sum_energy= 0.5*( // still need add constants
		c1_44 * s1um44 + c1_43 * s1um43 + c1_42 * s1um42 + c1_41 * s1um41 +
		c1_33 * s1um33 + c1_32 * s1um32 + c1_31 * s1um31 +
		c1_22 * s1um22 + c1_21 * s1um21 +
		c1_11 * s1um11 +
		c2_44 * s2um44 + c2_43 * s2um43 + c2_42 * s2um42 + c2_41 * s2um41 +
		c2_33 * s2um33 + c2_32 * s2um32 + c2_31 * s2um31 +
		c2_22 * s2um22 + c2_21 * s2um21 +
		c2_11 * s2um11	);
	return sum_energy;
}

int class_events::powc(int base, int index){
	if	(2==base){
		switch(index){
			case 1: return  2; break;
			case 2: return  4; break;
			case 3: return  8; break;
			case 4: return 16; break;
			default: error(2, "(pow_cal) an unknown index", 1, index);
		}
	}
	else if (1==base){return 1;}
	else if (0==base){return 0;}
	else if(-1==base){
		switch(index){
			case 1: return -1; break;
			case 2: return  1; break;
			case 3: return -1; break;
			case 4: return  1; break;
			default: error(2, "(pow_cal) an unknown index", 1, index);
		}
	}
	else if(-2==base){
		switch(index){
			case 1: return -2; break;
			case 2: return  4; break;
			case 3: return -8; break;
			case 4: return 16; break;
			default: error(2, "(pow_cal) an unknown index", 1, index);
		}
	}
	else error(2, "(pow_cal) an unknown base", 1, base);
}

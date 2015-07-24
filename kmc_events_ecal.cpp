#include <cstdio>
#include <iostream>
#include "kmc_global.h"
#include "kmc_events.h"

using namespace std;

double class_events::cal_energy(bool is_itl, int x1, int y1, int z1, int x2, int y2, int z2) const{ // calculate the energy of the input 2 ltcps 
	int pos[2][3]= {{x1, y1, z1}, {x2, y2, z2}};

	// AB-A, AB-B bonds, types of the 2 atoms (for C00)
	int ABA1= 0, ABB1= 0;
	int ABA2= 0, ABB2= 0;
	int type[2];
	
	// sums of 1st nn
	int s1um44=0, s1um42=0, s1um22=0;		// class 1
	int s1um43=0, s1um41=0, s1um32=0, s1um21=0;	// class 2
	int s1um33=0, s1um31=0, s1um11=0;		// class 3
	int s1um40=0, s1um30=0, s1um20=0, s1um10=0;
	// sums of 2nd nn
	int s2um44=0, s2um42=0, s2um22=0;		// class 1
	int s2um43=0, s2um41=0, s2um32=0, s2um21=0;	// class 2
	int s2um33=0, s2um31=0, s2um11=0;		// class 3
	int s2um40=0, s2um30=0, s2um20=0, s2um10=0;
	for(int i=0; i<2; i ++){
		int xi= pbc(pos[i][0], nx);
		int yi= pbc(pos[i][1], ny);
		int zi= pbc(pos[i][2], nz);

		int state_i= states[xi][yi][zi];

		if(! is_itl && 0==state_i) continue;

		// Calculations for 1st nn
		if(is_itl){ // because only in AB itl these terms appear in diff
			type[i]= state_i;
			s1um40 += n1nbr * state_i * state_i * state_i * state_i; // altho it's i**4 + j**4 inside n1nbr loop
			s1um30 += n1nbr * state_i * state_i * state_i;		 // actual diff is 8 times state i
			s1um20 += n1nbr * state_i * state_i;
			s1um10 += n1nbr * state_i;
		}

		for(int a=0; a<n1nbr; a ++){ // 1st neighbors
			int xj= pbc(xi+(*(v1nbr+a))[0], nx);
			int yj= pbc(yi+(*(v1nbr+a))[1], ny);
			int zj= pbc(zi+(*(v1nbr+a))[2], nz);

			if((xj==x1) && (yj==y1) && (zj==z1)) continue; // avoid double count the (x1, y1, z1)-(x2, y2, z2) bond

			int state_j1= states[xj][yj][zj];
			if(is_itl){
				if((0==state_i &&  1==state_j1) || (0==state_j1 &&  1==state_i)) ABA1 ++;
				if((0==state_i && -1==state_j1) || (0==state_j1 && -1==state_i)) ABB1 ++;
			}
			
			if(0==state_j1) continue;
			s1um44 += state_i * state_i * state_i * state_i * state_j1 * state_j1 * state_j1 * state_j1;
			s1um43 += state_i * state_i * state_i * state_i * state_j1 * state_j1 * state_j1
			        + state_i * state_i * state_i		* state_j1 * state_j1 * state_j1 * state_j1;
			s1um42 += state_i * state_i * state_i * state_i * state_j1 * state_j1
			        + state_i * state_i			* state_j1 * state_j1 * state_j1 * state_j1;
			s1um41 += state_i * state_i * state_i * state_i * state_j1
			        + state_i				* state_j1 * state_j1 * state_j1 * state_j1;
			s1um33 += state_i * state_i * state_i		* state_j1 * state_j1 * state_j1;
			s1um32 += state_i * state_i * state_i		* state_j1 * state_j1
			        + state_i * state_i			* state_j1 * state_j1 * state_j1;
			s1um31 += state_i * state_i * state_i		* state_j1
			        + state_i				* state_j1 * state_j1 * state_j1;
			s1um22 += state_i * state_i			* state_j1 * state_j1;
			s1um21 += state_i * state_i			* state_j1 
				+ state_i				* state_j1 * state_j1;
			s1um11 += state_i				* state_j1;
		}
				
if(! is_e2nbr) continue;

		// Calculations for 2nd nn
		if(is_itl){
			s2um40 += n2nbr * state_i * state_i * state_i * state_i; // altho it's i**4 + j**4 inside n1nbr loop
			s2um30 += n2nbr * state_i * state_i * state_i;		 // actual diff is 8 times state i
			s2um20 += n2nbr * state_i * state_i;
			s2um10 += n2nbr * state_i;
		}
		
		for(int b=0; b<n2nbr; b ++){ // 2nd neighbors
			int xj= pbc(xi+(*(v2nbr+b))[0], nx);
			int yj= pbc(yi+(*(v2nbr+b))[1], ny);
			int zj= pbc(zi+(*(v2nbr+b))[2], nz);

			int state_j2= states[xj][yj][zj];
			if(is_itl){
				if((0==state_i &&  1==state_j2) || (0==state_j2 &&  1==state_i)) ABA2 ++;
				if((0==state_i && -1==state_j2) || (0==state_j2 && -1==state_i)) ABB2 ++;
			}
			
			if(0==state_j2) continue;
			s2um44 += state_i * state_i * state_i * state_i * state_j2 * state_j2 * state_j2 * state_j2;
			s2um43 += state_i * state_i * state_i * state_i * state_j2 * state_j2 * state_j2
			        + state_i * state_i * state_i		* state_j2 * state_j2 * state_j2 * state_j2;
			s2um42 += state_i * state_i * state_i * state_i * state_j2 * state_j2
			        + state_i * state_i			* state_j2 * state_j2 * state_j2 * state_j2;
			s2um41 += state_i * state_i * state_i * state_i * state_j2
			        + state_i				* state_j2 * state_j2 * state_j2 * state_j2;
			s2um33 += state_i * state_i * state_i		* state_j2 * state_j2 * state_j2;
			s2um32 += state_i * state_i * state_i		* state_j2 * state_j2
			        + state_i * state_i			* state_j2 * state_j2 * state_j2;
			s2um31 += state_i * state_i * state_i		* state_j2
			        + state_i				* state_j2 * state_j2 * state_j2;
			s2um22 += state_i * state_i			* state_j2 * state_j2;
			s2um21 += state_i * state_i			* state_j2 
				+ state_i				* state_j2 * state_j2;
			s2um11 += state_i				* state_j2;
		}
	}

	// Summing up the terms
	double sum_energy= // pair terms
	  c1_44 * s1um44 + c1_43 * s1um43 + c1_42 * s1um42 + c1_41 * s1um41 +
	  c1_33 * s1um33 + c1_32 * s1um32 + c1_31 * s1um31 +
	  c1_22 * s1um22 + c1_21 * s1um21 + c1_11 * s1um11 +
	  c2_44 * s2um44 + c2_43 * s2um43 + c2_42 * s2um42 + c2_41 * s2um41 +
	  c2_33 * s2um33 + c2_32 * s2um32 + c2_31 * s2um31 +
	  c2_22 * s2um22 + c2_21 * s2um21 + c2_11 * s2um11;
	if(is_itl){ // C40, C30, C20, C10, C00
		double c00 = c1_00_ABA * ABA1 + c1_00_ABB * ABB1;
		if(is_e2nbr) c00+= c2_00_ABA * ABA2 + c2_00_ABB * ABB2;
		for(int i=0; i<2; i ++){
			switch(type[i]){
				case  2: c00 += c1_00_AA; break;
				case  1: c00 += c1_00_A ; break;
				case  0: c00 += c1_00_AB; break;
				case -1: c00 += c1_00_B ; break;
				case -2: c00 += c1_00_BB; break;
			}
			if(is_e2nbr){
				switch(type[i]){
					case  2: c00 += c2_00_AA; break;
					case  1: c00 += c2_00_A ; break;
					case  0: c00 += c2_00_AB; break;
					case -1: c00 += c2_00_B ; break;
					case -2: c00 += c2_00_BB; break;
				}
			}
		}
		
		sum_energy += c1_40 * s1um40 + c1_30 * s1um30 + c1_20 * s1um20 + c1_10 * s1um10;
		if(is_e2nbr) sum_energy += c2_40 * s2um40 + c2_30 * s2um30 + c2_20 * s2um20 + c2_10 * s2um10;
		sum_energy += c00;
	}

	return sum_energy;
}

int class_events::powc(int base, int index) const{
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

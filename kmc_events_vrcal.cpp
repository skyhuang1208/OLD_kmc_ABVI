#include <iostream>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"

using namespace std;

double class_events::cal_ratesV(vector <bool> &isvcc, vector <double> &rates, vector <int> &ilist, vector <int> &inbrs, vector <int> &jatom){
	double sum_rate= 0;
	if(nV != list_vcc.size()) error(2, "(cal_ratesV) vcc number inconsistent");
	
	for(int ivcc=0; ivcc < nV; ivcc ++){
		int i= (int) (list_vcc[ivcc].ltcp/nz)/ny;
		int j= (int) (list_vcc[ivcc].ltcp/nz)%ny;
		int k= (int)  list_vcc[ivcc].ltcp%nz;

		if(states[i][j][k] != 0) error(2, "(cal_ratesV) there's an non-vacancy in the vacancy list");
		
		for(int a=0; a<n1nbr; a ++){
			int x= pbc(i+v1nbr[a][0], nx);
			int y= pbc(j+v1nbr[a][1], ny);
			int z= pbc(k+v1nbr[a][2], nz);
			
			if(1==states[x][y][z] || -1==states[x][y][z]){
				double em, mu;

				if     (states[x][y][z]== 1) { em= emvA; mu= muA;} 
				else if(states[x][y][z]==-1) { em= emvB; mu= muB;} 
				else error(2, "(cal_ratesV) wrong type for V to jump into", 1, states[x][y][z]);

				double e0= cal_energy(i, j, k, x, y, z);

				states[i][j][k]= states[x][y][z];
				states[x][y][z]= 0;

				double ediff= cal_energy(i, j, k, x, y, z) - e0;
				
				states[x][y][z]= states[i][j][k];
				states[i][j][k]= 0;
				
				rates.push_back(mu * exp(-beta*(em+0.5*ediff)));
							
				isvcc.push_back(true);
				ilist.push_back(ivcc);
				inbrs.push_back(a);
				jatom.push_back(0);
				
				sum_rate += rates.back();
			}
		}
	}

	return sum_rate;
}

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "kmc_initial.h"
#include "kmc_par.h"

#define MAX_NNBR 20

using namespace std;

void class_initial::ltc_constructor(){
	double (*ptr_vbra)[3]; 
	int   (*ptr_v1nbr)[3];
	int   (*ptr_v2nbr)[3];
			
	// coordinate vectors of bravais lattices
	double vbra_bcc[3][3]= {{-0.5,  0.5,  0.5}, { 0.5, -0.5,  0.5}, { 0.5,  0.5, -0.5}};
		
	// BCC: index vectors of negibhor atoms: 1st nearest
	int v1nbr_bcc[8][3]= {{ 1,  0,  0}, { 0,  1,  0}, { 0,  0,  1}, { 1,  1,  1},
 		              {-1,  0,  0}, { 0, -1,  0}, { 0,  0, -1}, {-1, -1, -1}};
	int v2nbr_bcc[6][3]= {{ 0,  1,  1}, { 1,  0,  1}, { 1,  1,  0},
 		              { 0, -1, -1}, {-1,  0, -1}, {-1, -1,  0}};

	// Choose ltc structure
	if     (strcmp(type_ltc, "SC ")==0){}
	else if(strcmp(type_ltc, "BCC")==0){ptr_vbra= vbra_bcc;	n1nbr=8; ptr_v1nbr= v1nbr_bcc; n2nbr=6; ptr_v2nbr= v2nbr_bcc;}
	else if(strcmp(type_ltc, "FCC")==0){}
	else if(strcmp(type_ltc, "HCP")==0){}
	else	error(1, "(ltc_constructor) coldn't match the lattice type", type_ltc);
			
	// assign array values
	for(int i=0; i<3; i ++){     // Bravice vectors
		for(int j=0; j<3; j ++){
			vbra[i][j]= (*(ptr_vbra+i))[j];
	}}
	for(int i=0; i<n1nbr; i ++){ // v1nbr
		for(int j=0; j<3; j ++){
			v1nbr[i][j]= (*(ptr_v1nbr+i))[j];
			if(i>=n1nbr/2)
				if(v1nbr[i][j] != -v1nbr[i-n1nbr/2][j]) error(1, "(ltc_constructor) v1nbr isn't symmetry");
		}
	}
	for(int i=0; i<n2nbr; i ++){ // v2nbr
		for(int j=0; j<3; j ++){
			v2nbr[i][j]= (*(ptr_v2nbr+i))[j];
			if(i>=n2nbr/2)
				if(v2nbr[i][j] != -v2nbr[i-n2nbr/2][j]) error(1, "(ltc_constructor) v1nbr isn't symmetry");
		}
	}
}

void class_initial::init_states_array(int nVset, double compA){
	// STATE 0: vacancy, 1: A atom, -1: B atom

	int putV= (nx*ny*nz)/nVset;
	int Vcount= 0;

	for(int i=0; i<nx*ny*nz; i++){	
		if((i%putV==0) && (Vcount<nVset)){
			Vcount ++;
			*(&states[0][0][0] + i)= 0;
		}
		else{
			double ran= ran_generator();

			if(ran < compA) 	*(&states[0][0][0] + i)= +1;
			else			*(&states[0][0][0] + i)= -1;
		}
	}

	nV= 0; nA= 0; nB= 0; nAA= 0; nBB= 0; nAB= 0;
	////////// CHECK //////////
	for(int i=0; i<nx*ny*nz; i++){ 
		if(*(&states[0][0][0]+i)==  0){
			vcc temp_vcc;
			list_vcc.push_back(temp_vcc);
			list_vcc[nV].ltcp= i;
			list_vcc[nV].ix= 0;
			list_vcc[nV].iy= 0;
			list_vcc[nV].iz= 0;
							nV ++;
		}
		else if(*(&states[0][0][0]+i)== +1)	nA ++;
		else if(*(&states[0][0][0]+i)== -1)	nB ++;
		else				error(1, "(init_states_array) a state type is unrecognizable", 1, *(&states[0][0][0]+i));
	}
	if(nV != nVset) error(1, "(init_states_array) The number of vacancies is not nVset", 2, nV, nVset);
#define TOL 0.01
	int nAtotal= nA + nB;
	if(abs((double) nA/nAtotal-compA) > TOL) error(1, "(init_states_array) the composition of generated conf is inconsistent of compA", 1, nA);
	////////// CHECK //////////
	
	cout << "The random solution configuration has been generated!" << endl;
	cout << "Vacancy: " << nV << endl;
	cout << "Atype A: " << nA << ", pct: " << 100* (double) nA/(nAtotal) << "%" << endl;
	cout << "Atype B: " << nB << ", pct: " << 100* (double) nB/(nAtotal) << "%" << endl;
}

void class_initial::read_restart(char name_restart[], long long int &ts_initial, double &time_initial){
	ifstream if_re(name_restart, ios::in);
	if(!if_re.is_open()) error(1, "(read_restart) the file is not opened!");
	
	int timestep;
	double time;

	int ntotal;
	if_re >> ntotal;
	if(ntotal != nx*ny*nz) error(1, "(read_restart) the input total ltc number isnt consistent", 2, ntotal, nx*ny*nz);
	
	char c_ltcp[5];
	if_re >> c_ltcp >> timestep >> time;
	if(strcmp(c_ltcp, "ltcp") !=0) error(1, "(read_restart) please input a ltcp file (ltcp at the second line)"); // check
	ts_initial= timestep;
	time_initial= time;

	nV= 0; nA= 0; nB= 0; nAA= 0; nBB= 0; nAB= 0;
	for(int index=0; index<nx*ny*nz; index ++){	
		int type, i, j, k, ix, iy, iz, dir, head;
		if_re >> type >> i >> j >> k;
		if(index != i*ny*nz+j*nz+k) error(1, "(read_restart) the input index inconsistent");
	
		if( 0==type){
			if_re >> ix >> iy >> iz;
			
			vcc temp_vcc;
			list_vcc.push_back(temp_vcc);
			list_vcc[nV].ltcp= index;
			list_vcc[nV].ix= ix;
			list_vcc[nV].iy= iy;
			list_vcc[nV].iz= iz;

			nV ++;
		}
		else if	( 1==type) nA ++;
		else if	(-1==type) nB ++;
		else{
			if_re >> ix >> iy >> iz >> dir >> head;

			itl temp_itl;
			list_itl.push_back(temp_itl);
			list_itl[nAA+nAB+nBB].ltcp= index;
			list_itl[nAA+nAB+nBB].ix= ix;
			list_itl[nAA+nAB+nBB].iy= iy;
			list_itl[nAA+nAB+nBB].iz= iz;
			list_itl[nAA+nAB+nBB].dir= dir;
			list_itl[nAA+nAB+nBB].head= head;
			
			if	(type== 2) nAA ++;
			else if (type==-2) nBB ++;
			else if (type== 3){nAB ++; type= 0; *(&itlAB[0][0][0]+index)= true;}
		}

		*(&states[0][0][0]+index)= type;
	}
	if(nV+nA+nB+nAA+nBB+nAB != nx*ny*nz) error(1, "(read_restart) the number inconsistent", 2, nV+nA+nB+nAA+nBB+nAB, nx*ny*nz);

	cout << "The configuration has been generated from the restart file!" << endl;
	cout << "Vacancy: " << nV << endl;
	cout << "Atype A: " << nA << ", pct: " << 100* (double)nA / ntotal << "%" << endl;
	cout << "Atype B: " << nB << ", pct: " << 100* (double)nB / ntotal << "%" << endl;
	cout << " Itl AA: " << nAA << endl;
	cout << " Itl AB: " << nAB << endl;
	cout << " Itl BB: " << nBB << endl;
	
	if_re.close();
}


void class_initial::init_par(){
	// 1st-nn class 1
	c1_44= ( (eAA1AA -8*eAA1A -8*eAA1B +12*eAA1AB +2*eAA1BB +12*eAB1BB -8*eA1BB -8*eB1BB +eBB1BB) +
	       (48*eA1AB -12*eAB1AB +48*eAB1B) + (-48*eA1V +48*eV1V -48*eV1B) + (16*eA1A +32*eA1B +16*eB1B) )/576; 
	 
	c1_42= ( (-eAA1AA +20*eAA1A +20*eAA1B -36*eAA1AB -2*eAA1BB -36*eAB1BB +20*eA1BB +20*eB1BB -eBB1BB) +
	       (-216*eA1AB +36*eAB1AB -216*eAB1B) + (216*eA1V -216*eV1V +216*eV1B) + (-64*eA1A -128*eA1B -64*eB1B) )/576; 

	c1_22= ( (eAA1AA -32*eAA1A -32*eAA1B +60*eAA1AB +2*eAA1BB +60*eAB1BB -32*eA1BB -32*eB1BB +eBB1BB) +
	       (960*eA1AB -60*eAB1AB +960*eAB1B) + (-960*eA1V +960*eV1V -960*eV1B) + (256*eA1A +512*eA1B +256*eB1B) )/576; 
	
	// 1st-nn class 2
	c1_43= ( (eAA1AA -6*eAA1A -2*eAA1B +6*eAA1AB -6*eAB1BB +2*eA1BB +6*eB1BB -eBB1BB) +
	       (12*eA1AB -12*eAB1B) + (-12*eA1V +12*eV1B) + (8*eA1A -8*eB1B) )/288; 
	 
	c1_41= ( (-eAA1AA +12*eAA1A -4*eAA1B -6*eAA1AB +6*eAB1BB +4*eA1BB -12*eB1BB +eBB1BB) +
	       (-48*eA1AB +48*eAB1B) + (48*eA1V -48*eV1B) + (-32*eA1A +32*eB1B) )/288; 
	
	c1_32= ( (-eAA1AA +18*eAA1A +14*eAA1B -30*eAA1AB +30*eAB1BB -14*eA1BB -18*eB1BB +eBB1BB) +
	       (-60*eA1AB +60*eAB1B) + (60*eA1V -60*eV1B) + (-32*eA1A +32*eB1B) )/288; 
	
	c1_21= ( (eAA1AA -24*eAA1A -8*eAA1B +30*eAA1AB -30*eAB1BB +8*eA1BB +24*eB1BB -eBB1BB) +
	       (240*eA1AB -240*eAB1B) + (-240*eA1V +240*eV1B) + (128*eA1A -128*eB1B) )/288; 
	
	// 1st-nn class 3
	c1_33= ( (eAA1AA -4*eAA1A +4*eAA1B -2*eAA1BB +4*eA1BB -4*eB1BB +eBB1BB) +
	       (4*eA1A -8*eA1B +4*eB1B) )/144; 
	
	c1_31= ( (-eAA1AA +10*eAA1A -10*eAA1B +2*eAA1BB -10*eA1BB +10*eB1BB -eBB1BB) +
	       (-16*eA1A +32*eA1B -16*eB1B) )/144; 
	
	c1_11= ( (eAA1AA -16*eAA1A +16*eAA1B -2*eAA1BB +16*eA1BB -16*eB1BB +eBB1BB) +
	       (64*eA1A -128*eA1B +64*eB1B) )/144; 
	
	// 2nd-nn class 1
	c2_44= ( (eAA2AA -8*eAA2A -8*eAA2B +12*eAA2AB +2*eAA2BB +12*eAB2BB -8*eA2BB -8*eB2BB +eBB2BB) +
	       (48*eA2AB -12*eAB2AB +48*eAB2B) + (-48*eA2V +48*eV2V -48*eV2B) + (16*eA2A +32*eA2B +16*eB2B) )/576; 
	 
	c2_42= ( (-eAA2AA +20*eAA2A +20*eAA2B -36*eAA2AB -2*eAA2BB -36*eAB2BB +20*eA2BB +20*eB2BB -eBB2BB) +
	       (-216*eA2AB +36*eAB2AB -216*eAB2B) + (216*eA2V -216*eV2V +216*eV2B) + (-64*eA2A -128*eA2B -64*eB2B) )/576; 

	c2_22= ( (eAA2AA -32*eAA2A -32*eAA2B +60*eAA2AB +2*eAA2BB +60*eAB2BB -32*eA2BB -32*eB2BB +eBB2BB) +
	       (960*eA2AB -60*eAB2AB +960*eAB2B) + (-960*eA2V +960*eV2V -960*eV2B) + (256*eA2A +512*eA2B +256*eB2B) )/576; 
	
	// 2nd-nn class 2
	c2_43= ( (eAA2AA -6*eAA2A -2*eAA2B +6*eAA2AB -6*eAB2BB +2*eA2BB +6*eB2BB -eBB2BB) +
	       (12*eA2AB -12*eAB2B) + (-12*eA2V +12*eV2B) + (8*eA2A -8*eB2B) )/288; 
	 
	c2_41= ( (-eAA2AA +12*eAA2A -4*eAA2B -6*eAA2AB +6*eAB2BB +4*eA2BB -12*eB2BB +eBB2BB) +
	       (-48*eA2AB +48*eAB2B) + (48*eA2V -48*eV2B) + (-32*eA2A +32*eB2B) )/288; 
	
	c2_32= ( (-eAA2AA +18*eAA2A +14*eAA2B -30*eAA2AB +30*eAB2BB -14*eA2BB -18*eB2BB +eBB2BB) +
	       (-60*eA2AB +60*eAB2B) + (60*eA2V -60*eV2B) + (-32*eA2A +32*eB2B) )/288; 
	
	c2_21= ( (eAA2AA -24*eAA2A -8*eAA2B +30*eAA2AB -30*eAB2BB +8*eA2BB +24*eB2BB -eBB2BB) +
	       (240*eA2AB -240*eAB2B) + (-240*eA2V +240*eV2B) + (128*eA2A -128*eB2B) )/288; 
	
	// 2nd-nn class 3
	c2_33= ( (eAA2AA -4*eAA2A +4*eAA2B -2*eAA2BB +4*eA2BB -4*eB2BB +eBB2BB) +
	       (4*eA2A -8*eA2B +4*eB2B) )/144; 
	
	c2_31= ( (-eAA2AA +10*eAA2A -10*eAA2B +2*eAA2BB -10*eA2BB +10*eB2BB -eBB2BB) +
	       (-16*eA2A +32*eA2B -16*eB2B) )/144; 
	
	c2_11= ( (eAA2AA -16*eAA2A +16*eAA2B -2*eAA2BB +16*eA2BB -16*eB2BB +eBB2BB) +
	       (64*eA2A -128*eA2B +64*eB2B) )/144; 

	if(0==c2_44 && 0==c2_43 && 0==c2_42 && 0==c2_41 && 0==c2_33 && 0==c2_32 && 0==c2_31 && 0==c2_22 && 0==c2_21 && 0==c2_11) is_e2nbr= false;
	else is_e2nbr= true;

	// print out the parameters to log file
	cout << "\n##### Energy calculation parameters #####" << endl; 
	
	cout << "beta= " << beta << endl;
	printf("mu= %f %f\n", muA, muB);
	printf("Vacancy Em= %f %f\n", emvA, emvB);
	printf("Interstitial Em= %f %f\n", emiA, emiB);
	printf("Rotation Er(AA, AB, BB)= %f %f %f\n", erAA, erAB, erBB);
	
	cout << "Input epsilons:" << endl;
	cout << "(1st neigbor)" << endl;
	printf("AA-AA: %f, AA-A: %f, AA-V: %f, AA-AB: %f, AA-B: %f, AA-BB: %f\n", eAA1AA, eAA1A, eAA1V, eAA1AB, eAA1B, eAA1BB);
	printf("A-A:   %f, A-V:  %f, A-AB: %f, A-B:   %f, A-BB: %f\n", eA1A, eA1V, eA1AB, eA1B, eA1BB);
	printf("V-V:   %f, V-AB: %f, V-B:  %f, V-BB:  %f\n", eV1V, eV1AB, eV1B, eV1BB);
	printf("AB-AB: %f, AB-B: %f, AB-BB: %f\n", eAB1AB, eAB1B, eAB1BB);
	printf("B-B:   %f, B-BB: %f\n", eB1B, eB1BB);
	printf("BB-BB: %f\n", eBB1BB);
	cout << "(2nd neigbor)" << endl;
	printf("AA-AA: %f, AA-A: %f, AA-V: %f, AA-AB: %f, AA-B: %f, AA-BB: %f\n", eAA2AA, eAA2A, eAA2V, eAA2AB, eAA2B, eAA2BB);
	printf("A-A:   %f, A-V:  %f, A-AB: %f, A-B:   %f, A-BB: %f\n", eA2A, eA2V, eA2AB, eA2B, eA2BB);
	printf("V-V:   %f, V-AB: %f, V-B:  %f, V-BB:  %f\n", eV2V, eV2AB, eV2B, eV2BB);
	printf("AB-AB: %f, AB-B: %f, AB-BB: %f\n", eAB2AB, eAB2B, eAB2BB);
	printf("B-B:   %f, B-BB: %f\n", eB2B, eB2BB);
	printf("BB-BB: %f\n", eBB2BB);
	
	cout << "Ising formulation constants:" << endl;
	cout << "(1st neighbor)" << endl;
	printf("Class 1\nC44: %f, C42: %f, C22: %f\n", c1_44, c1_42, c1_22);
	printf("Class 2\nC43: %f, C41: %f, C32: %f, C21: %f\n", c1_43, c1_41, c1_32, c1_21);
	printf("Class 3\nC33: %f, C31: %f, C11: %f\n", c1_33, c1_31, c1_11);
	cout << "(2nd neighbor)" << endl;
	printf("Class 1\nC44: %f, C42: %f, C22: %f\n", c2_44, c2_42, c2_22);
	printf("Class 2\nC43: %f, C41: %f, C32: %f, C21: %f\n", c2_43, c2_41, c2_32, c2_21);
	printf("Class 3\nC33: %f, C31: %f, C11: %f\n", c2_33, c2_31, c2_11);

	if(is_e2nbr) 	cout << "\n2nd nn parameters are non-zero" << endl;
	else		cout << "\n2nd nn are 0, skip 2nd-nn calculations" << endl;
}


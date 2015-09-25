#ifndef KMC_EVENTS_INCLUDED
#define KMC_EVENTS_INCLUDED
#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

class class_events{
	public:
		class_events(){
			rrecb_int= (int) dis_rec;
			rrecb2= dis_rec * dis_rec;
			double nnd= sqrt(vbra[0][0]*vbra[0][0]+vbra[0][1]*vbra[0][1]+vbra[0][2]*vbra[0][2]);
			rrecb_nnd= (int) (dis_rec/nnd);

			cout << "##Generation parameters: (rate_genr)  " << rate_genr << " (damgae/s)" << endl;
			cout << "##Recombination parameters: (distance) " << dis_rec << ", (int) " << rrecb_int << ", (in nearest-neighbor distance) " << rrecb_nnd << endl;
		}
		
		// functions
		void genr();
		double jump();
		double ecal_whole() const; 
	
	private:
		// Variables for recombination check
		int rrecb_int;	// integer of recombination distance
		double rrecb2;	// square of recombination distance
		int rrecb_nnd;	// recombination distance in nearest neighbor distance

		////// functions of energy calculation //////
		double cal_energy(bool is_itl, int x1, int y1, int z1, int x2, int y2, int z2) const; 
		int powc(int base, int index) const;
		double cal_c00(int type[], int ABA1, int ABB1, int ABA2, int ABB2) const;
		double ecal_bond(bool is_itl, int x1, int y1, int z1, int x2, int y2, int z2) const; 
		double ecal_otf (bool is_itl, int x1, int y1, int z1, int x2, int y2, int z2) const; // corrected H on the fly 

		////// functions for jumps(jumping rate calculations) //////
		double cal_ratesV(vector <bool> &isvcc, vector <double> &rates, vector <int> &ilist, vector <int> &inbrs, vector <int> &jatom);
		double cal_ratesI(vector <bool> &isvcc, vector <double> &rates, vector <int> &ilist, vector <int> &inbrs, vector <int> &jatom);
		void itl_rules(int &itl, int &atom, int jatom);
		void actual_jumpV(int vid, int nid);
		void actual_jumpI(int iid, int nid, int jatom);
		
		///// functions for recombination /////
		void rules_recb(int ii, int iv);
		void recb_dir(int index);
		bool cal_dis(int d1, int d2, int d3);
		void recb_randomV(int index);
		bool recb_randomI(int index);
		void sink(bool isvcc, int index);
};

#endif // KMC_EVENTS_INCLUDED

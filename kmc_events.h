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

			cout << "recombination parameters: " << endl;
			cout << "recombination distance: " << dis_rec << ", (int) " << rrecb_int << ", (in nearest-neighbor distance) " << rrecb_nnd << endl;
		}
		
		// functions
		void genr();
		double jump();
		double cal_energy();
	
	private:
		// Variables for recombination check
		int rrecb_int;	// integer of recombination distance
		double rrecb2;	// square of recombination distance
		int rrecb_nnd;	// recombination distance in nearest neighbor distance

		
		////// functions of energy calculation //////
		double cal_energy(int x1, int y1, int z1, int x2, int y2, int z2); 
		int powc(int base, int index);

		
		////// functions for jumps(jumping rate calculations) //////
		double cal_ratesV(vector <bool> &isvcc, vector <double> &rates, vector <int> &ilist, vector <int> &inbrs, vector <int> &jatom);
		double cal_ratesI(vector <bool> &isvcc, vector <double> &rates, vector <int> &ilist, vector <int> &inbrs, vector <int> &jatom);
		void itl_rules(int &itl, int &atom, int jatom);
		void actual_jumpV(int dx, int dy, int dz, int x, int y, int z, int did);
		void actual_jumpI(int dx, int dy, int dz, int x, int y, int z, int did, int jatom, int inbr);
		
		///// functions for recombination /////
		void recb_rules(int ii, int iv);
		void recb_dir(int index);
		bool cal_dis(int d1, int d2, int d3);
		void recb_randomV(int index);
		void recb_randomI(int index);
};

#endif // KMC_EVENTS_INCLUDED

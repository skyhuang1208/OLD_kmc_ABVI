#ifndef KMC_EVENTS_INCLUDED
#define KMC_EVENTS_INCLUDED
#include <vector>
#include <iostream>

using namespace std;

class class_events{
	public:
		class_events(int step_write_his_, int *states_): step_write_his(step_write_his_), states(states_) 
		{
			actions_sol[0].reserve((int) (0.1*step_write_his)); // should be much larger than need
			actions_sol[1].reserve((int) (0.1*step_write_his));
		}
		
		// functions
		void events_main();
		double cal_energy(int* const states_ce);
	
	private:
		int step_write_his;	// timesteps for output to history files

		int* const states;	// Don't change the address of this!

		vector <int> actions_sol[2]; // A list contains solute atom moves from [0] to [1]
		
		//////functions for all events//////
		double cal_energy(int x1, int y1, int z1, int x2, int y2, int z2); 

		//////functions for vacancies //////
		double vac_jump(vector <double> &v_rate, vector <int> &v_ivcc, vector <int> &v_inbr);
		void vac_recb(int vpos[3]);

		//////functions for interstitials////
		void int_motions();
		int  int_eval(int x_int, int y_int, int z_int);
		void int_jump(int x_begin, int y_begin, int z_begin);
		
		/////sorts of rules /////
		void rules_int_jump(int &typei, int &typef);	// initial to final
		void     rules_recb(int &typei, int &typev);	// interstitial and vacancy
};


#endif // KMC_EVENTS_INCLUDED

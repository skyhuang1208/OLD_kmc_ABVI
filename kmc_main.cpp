#include <iostream>
#include <fstream>
#include <ctime>
#include "kmc_par.h"
#include "kmc_global.h"
#include "kmc_initial.h"
#include "kmc_events.h"
using namespace std;

int main(int nArg, char *Arg[]){
	int t0cpu= time(0);

	long long int ts_bg;
	double time_bg;
	int states[nx][ny][nz];
	

	cout << "########## Initializing System... ##########" << endl;
	class_initial init(&states[0][0][0], ts_bg, time_bg, Arg[1], par_name_sol, par_name_vcc);
	
	cout << "\n########## Initializing Events ... ##########" << endl;
	class_events events(par_step_write_his, &states[0][0][0]); 

	cout << "\n########## The Simulation Begins !! ##########" << endl;
	timestep=  ts_bg; totaltime= time_bg;
	int STEP_LOG= par_tstep/2000; cout << "TIMESTEP() TIME(s) ETOTAL(eV)";
	while((totaltime<= time_bg+par_tend) && (timestep != ts_bg+par_tstep)){
		// CALCULATION
		events.events_main();

		// OUTPUT DATA
		if(0==timestep%STEP_LOG)
			cout << endl << timestep << " " << totaltime << " " << events.cal_energy(&states[0][0][0]) << " ";
		
		if(0==timestep%par_confts){
			write_conf(&states[0][0][0]);
			cout << "<Output conf files at: " << timestep << ">";
		}
	}

	// finalizing
	if(timestep%STEP_LOG != 0)   
		cout << endl << timestep << " " << totaltime << " " << events.cal_energy(&states[0][0][0]) << endl;
	if(timestep%par_confts != 0){ 
		write_conf(&states[0][0][0]); 
		cout << "Output conf files at: " << timestep;
	}

	int tfcpu= time(0);
	cout << "\n**** The simulation is done! Total CPU time: " << tfcpu - t0cpu << " secs ****" << endl;

	return 0;
}

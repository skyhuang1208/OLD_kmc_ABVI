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

	cout << "########## Initializing System... ##########" << endl;
	class_initial init(ts_bg, time_bg, Arg[1]);
	
	cout << "\n########## Initializing Events ... ##########" << endl;
	class_events events; 

	cout << "\n########## The Simulation Begins !! ##########" << endl;
	timestep=  ts_bg; totaltime= time_bg;
	int STEP_LOG= par_step/2000; 
	int N_genr;
	cout << "TIMESTEP() TIME(s) ETOTAL(eV) #GENR() NAA() NA() NV() NAB() NB() NBB()";
	while((totaltime<= time_bg+par_time) && (timestep != ts_bg+par_step)){
		// CALCULATIONS
		int int_Tprev = (int) (totaltime/par_time_genr);
		totaltime += events.jump();
		int int_Tnow  = (int) (totaltime/par_time_genr);
		timestep ++;
		if(int_Tprev != int_Tnow || timestep==par_step_genr){
			events.genr();
			N_genr ++;
		}

		// OUTPUT DATA
		if(0==timestep%STEP_LOG)
			cout << endl << timestep << " " << totaltime << " " << events.cal_energy() << " " 
			     << N_genr << " " << nAA << " " << nA << " " << nV << " " << nAB << " " << nB << " " << nBB;
		if(0==timestep%par_confts){
			write_conf();
			cout << "<Output conf files at: " << timestep << ">";
		}
		if(0==timestep%par_his){
			write_hisdef();
 			write_hissol();
	
			actions_sol[0].clear(); actions_sol[1].clear();
		}
	}

	// finalizing
	if(timestep%STEP_LOG != 0)   
		cout << endl << timestep << " " << totaltime << " " << events.cal_energy() << endl;
	if(timestep%par_confts != 0){ 
		write_conf(); 
		cout << "Output conf files at: " << timestep;
	}

	int tfcpu= time(0);
	cout << "\n**** The simulation is done! Total CPU time: " << tfcpu - t0cpu << " secs ****" << endl;

	return 0;
}

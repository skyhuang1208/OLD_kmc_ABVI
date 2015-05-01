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
	class_initial init(ts_bg, time_bg, nArg, Arg[1]);
	
	cout << "\n########## Initializing Events ... ##########" << endl;
	class_events events; 

	cout << "\n########## The Simulation Begins !! ##########" << endl;
	timestep=  ts_bg; totaltime= time_bg;
	int N_genr= 0;
	cout << "TIMESTEP() TIME(s) GENR()	NA() NB()	NV() NAA() NAB() NBB()";
	while((totaltime<= time_bg+par_time) && (timestep != ts_bg+par_step)){
		// CALCULATIONS
//		int int_t0, int_t1;
//		int_t1  = (int) (totaltime/par_time_genr);
//		if(int_t0 < int_t1){ // Generations
//			for(int i=1; i<=(int_t1-int_t0); i ++) events.genr();
//			N_genr += (int_t1-int_t0);
//		}
//		int_t0 = int_t1;

		if(0==timestep%par_step_genr){
			events.genr();
			N_genr ++;
		} 
		
		totaltime += events.jump(); // Defect jumps
		timestep ++;

		// OUTPUT DATA
		if(0==timestep%step_log)
			printf("\n%lld %e %d	%d %d	%d %d %d %d ", timestep, totaltime, N_genr, nA, nB, nV, nAA, nAB, nBB);
		if(0==timestep%step_confts){
			write_conf();
			cout << "<Output conf files at: " << timestep << ">";
		}
		if(0==timestep%step_out)
			fprintf(out_engy, "%lld %e %f\n", timestep, totaltime, events.cal_energy());
		if(0==timestep%step_his){
			write_hisdef();
 			write_hissol();
	
			actions_sol[0].clear(); actions_sol[1].clear();
		}
		
	}

	// finalizing
	if(timestep%step_log != 0)   
		printf("\n%lld %f %d	%d %d	%d %d %d %d ", timestep, totaltime, N_genr, nA, nB, nV, nAA, nAB, nBB);
	if(timestep%step_confts != 0){ 
		write_conf(); 
		cout << "<Output conf files at: " << timestep << ">";
	}
	if(timestep%step_out != 0)
		fprintf(out_engy, "%lld %e %f\n", timestep, totaltime, events.cal_energy());

	int tfcpu= time(0);
	cout << "\n**** The simulation is done! Total CPU time: " << tfcpu - t0cpu << " secs ****" << endl;

	return 0;
}

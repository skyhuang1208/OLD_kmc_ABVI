#ifndef KMC_SYSTEM_INCLUDED
#define KMC_SYSTEM_INCLUDED
#include <iostream>
#include <cstring>
#include "kmc_global.h"
#include "kmc_par.h"

#define MAX_NSPE 20

using namespace std;

class class_initial{
	public:
		class_initial(int *states_, long long int &ts_bg, double &time_bg, char name_restart[], const char name_sol[], const char name_vcc[]): states000(states_)
		{
			stpcpy(type_ltc, par_ltc);

			ltc_constructor();	// execute constructor for vbra, n*nbr, and v*nbr
			
			// Printf the parameters
			cout << "Setting: " << endl << "nx= " << nx << ", ny= " << ny << ", nz= " << nz << endl;
			cout << "The crystal structure is: " << type_ltc << endl;
			for(int i=0; i<3; i++)
				cout << "v" << i << ": " << vbra[i][0] << " " << vbra[i][1] << " " << vbra[i][2] << endl;
			cout << "And the number of 1st neighbors: " << n1nbr << endl;
			cout << "    the number of 2nd neighbors: " << n2nbr << endl;
		
			cout << "\n####### Generating configuration... #######" << endl;
			if(par_isrestart){
				cout << "RESTART FROM restart file..." << endl;
				if(name_restart == NULL) error(0, "when restart flag is true, nArg must be 2");
				read_restart(name_restart, ts_bg, time_bg);
				his_sol= fopen(name_sol, "a");
				his_vcc= fopen(name_vcc, "a");
				cout << "Open " << name_sol << " & " << name_vcc << " with append mode" << endl;
			}
			else{
				cout << "START FROM a random configuration..." << endl;
				ts_bg= 0; time_bg= 0;
				init_states_array(par_nV, par_compA);
				write_conf(states000);
				cout << "Output t0 conf files" << endl;
				his_sol= fopen(name_sol, "w");
				his_vcc= fopen(name_vcc, "w");
				cout << "Open " << name_sol << " & " << name_vcc << "with write mode" << endl;
			}
			if(NULL==his_sol) error(2, "(class_events) the solute  history file was not opened!");
			if(NULL==his_vcc) error(2, "(class_events) the vacancy history file was not opened!");
			
			init_par();
			init_list_vcc();
		}
		
	private:
		char type_ltc[4];	// crystal structure type
		int* const states000;	// Don't change the address of this!

		// functions
		void ltc_constructor();
		void init_states_array(int nVset, double compA);
		void read_restart(char name_restart[], long long int &ts_initial, double &time_initial);
		void init_par();
		void init_list_vcc();
};
#endif // KMC_SYSTEM_INCLUDED

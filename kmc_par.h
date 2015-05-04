#ifndef KMC_PAR_INCLUDED
#define KMC_PAR_INCLUDED

const char   par_ltc[4]=		"BCC";

const int    par_nx=                      256;
const int    par_ny=                       64;
const int    par_nz=                       64;

const double par_compA=                  0.95; // composition of A atoms
const int    par_nV=                        0;

const double	    par_time=            1e04; 	// toal time (s)
const long long int par_step=            2e10; 	// toal timestep (give a minus step to ignore this quiterior to end the simulation)

const long long int step_log=	          5e4; 
const long long int step_confts= par_step/5e0;	// timestep that output a conf file for restart later
const long long int step_out=	 par_step/5e2; 
const long long int step_his=		  1e4;


const bool par_isrestart=		    false;

const char   par_name_sol[20]=      "history.sol";
const char   par_name_def[20]=      "history.def";
const char   par_name_engy[20]=      "out.energy";

// Parameters for events
const double par_dis_rec=0.866*3; // recombination distance
const double par_dpasm1=    1e-6; 
const long long int par_step_genr=     1e1; 

// Ising model energy calculation parameters
const double par_temp=                       500.0;
const double par_beta= 1.0/par_temp/8.617332478e-5; // 1/kbT, units: eV, K

const double par_muA=			   5.0e+15; // units: s^-1 
const double par_muB=			   5.0e+15; // units: s^-1 

const double par_emvA=			       1.0;
const double par_emvB=			       1.0;
const double par_emiA=			       0.5;
const double par_emiB=			       0.5;

const double par_erAA=				 0; // rotation energy barrier
const double par_erAB=				 0;
const double par_erBB=				 0;

// bonding energy parameters

// 1st nn
const double eAA1AA=                         0;
const double eAA1A=                          0;
const double eAA1V=                          0;
const double eAA1AB=                         0;
const double eAA1B=                          0;
const double eAA1BB=                         0;
// ---
const double eA1A=                      -1.070;
const double eA1V=                      -0.300;
const double eA1AB=                          0;
const double eA1B=                      -0.985;
const double eA1BB=                          0;
// ---
const double eV1V=                           0;
const double eV1AB=                          0;
const double eV1B=                      -0.300;
const double eV1BB=                          0;
// ---
const double eAB1AB=                         0;
const double eAB1B=                          0;
const double eAB1BB=                         0;
// ---
const double eB1B=                      -1.070;
const double eB1BB=                          0;
// ---
const double eBB1BB=                         0;

// 2nd nn
const double eAA2AA=                         0;
const double eAA2A=                          0;
const double eAA2V=                          0;
const double eAA2AB=                         0;
const double eAA2B=                          0;
const double eAA2BB=                         0;
// ---
const double eA2A=                           0;
const double eA2V=                           0;
const double eA2AB=                          0;
const double eA2B=                           0;
const double eA2BB=                          0;
// ---
const double eV2V=                           0;
const double eV2AB=                          0;
const double eV2B=                           0;
const double eV2BB=                          0;
// ---
const double eAB2AB=                         0;
const double eAB2B=                          0;
const double eAB2BB=                         0;
// ---
const double eB2B=                           0;
const double eB2BB=                          0;
// ---
const double eBB2BB=                         0;

// trapping number: solute atom trapping and intersitial trapping							 
// 	const int par_trNsol= 3; // trapping number by solute atoms
// 	const int par_trNint= 3; // trapping number by interstitials
#endif // KMC_PAR_INCLUDED

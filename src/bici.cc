// This is start code for BICI


// BICI (Bayesian individual-based compartmental inference)
// © C. M. Pooley and G. Marion

// Compilation: make
// Run: ./bici 

// ssh gaia.bioss.ac.uk  

//valgrind --tool=massif --exit-on-first-error=yes --error-exitcode=1 -s ./bici-para Execute/init.bici sim

// tar -xzf foo.tgz

// Load mpi: module load mpi/openmpi-x86_64
// ./bici-core.exe ../big2.bici sim
// mpirun -n 1 ./bici-para file.bici sim
// mpirun -n 1 ./bici-para Execute/init.bici sim
// mpirun --output :raw -n 1 ./bici-para Execute/init.bici sim
// mpirun -n 3 ./bici-para Execute/init.bici inf
// mpirun -n 10 ./bici-para Execute/init.bici inf
// mpirun -n 1 ./bici-para Examples/EX_A1.bici sim
// mpirun -n 2 ./bici-para Execute/init.bici inf
// mpirun -n 16 ./bici-para Jamie/scen-3-1.bici inf
// mpirun -n 16 ./bici-para Execute/init.bici inf
// mpirun -n 3 ./bici-para Grant/RealData_1000_DA_5000 inf
// mpirun -n 3 ./bici-para Grant/RealData_2000_DA_5000 inf
// mpirun -n 16 ./bici-para Jamie/scen-5-1.bici inf
// nohup mpirun -n 16 ./bici-para Jamie/scen-1-1b.bici inf > big2.txt&
// mpirun -n 10 ./bici-para Execute/init.bici inf
// ./bici-para Geno/SI_catbad_g1.0.bici inf > q1&
// ./bici-para Geno/farm_f1.0.bici inf > q1&
// ./bici-para Geno/SI_catbad_m1.bici inf > y1&

// ./bici-core Execute/init.bici inf

//valgrind --exit-on-first-error=yes --error-exitcode=1 --leak-check=yes -s ./bici-para Execute/init.bici sim

//valgrind --exit-on-first-error=yes --error-exitcode=1 --leak-check=yes -s ./bici-para Execute/init.bici inf

// nohup mpirun -n 16 ./bici-para Jamie/scen-1-1.bici inf > ddddd.txt&
// 06d_New_Inf_500.bici
// nohup mpirun -n 16 ./bici-para Jamie/simp_model.bici inf > eeee.txt&
// 06d_New_Inf_500.bici
// ./bici-para Jamie/test.bici post-sim
// ./bici-para Jamie/testped2.bici post-sim
// valgrind --exit-on-first-error=yes --error-exitcode=1 --leak-check=yes -s  ./bici-para Execute/init.bici inf
// valgrind --exit-on-first-error=yes --error-exitcode=1 --leak-check=yes -s  ./bici-para Execute/init.bici sim
// ./bici-para Jamie/temp4/scen-4-11.bici inf
// mpirun -n 10 ./bici-para Jamie/temp4/scen-4-11.bici inf
//  valgrind --exit-on-first-error=yes --error-exitcode=1 --leak-check=yes -s ./bici-para Jamie/temp4/scen-4-11.bici inf

// nohup valgrind --exit-on-first-error=yes --error-exitcode=1 --leak-check=yes -s ./bici-para Execute/mod13.bici inf > moddd13.txt

// nohup gdb -ex=r ./bici-para  --args Execute/mod14.bici inf > moddd14.txt

//nohup  mpirun -n 10  ./bici-para Execute/mod10.bici inf >  modd10.txt&
// mpirun -n 4  ./bici-para Execute/init.bici inf

// git clone https://github.com/theITEAM/BICI.git

// Different commands
// -co=0 / -core=0                            Runs the zeroth core
// -samp=0                                    Sets how sampling is done (for diagnostic purposes)
// -seed=100                                  This overides seed set in bici file 
// -sim=1                                     Selects a particular simulation when creating data

// When running in parallel then load mpi using: module load mpi/openmpi-x86_64

// Test: valgrind --exit-on-first-error=yes --error-exitcode=1 --leak-check=yes -s ./bici-para Execute/init.bici sim
// Test: valgrind --exit-on-first-error=yes --error-exitcode=1 --leak-check=yes -s ./bici-para Execute/init.bici inf
// Test:  valgrind --exit-on-first-error=yes --error-exitcode=1 --leak-check=yes -s  ./bici-para simple_model.bici inf

// Compiling for windows: .\make.bat from the directory D:\C_complier_BICI2.0\bin

// 45541 lines of code (15/11/24)
// 51000 lines of code (18/12/24)
// 56000 lines of code (16/05/25)
// 65634 lines of code (17/09/25)
// 81002 lines of code (25/02/26)

#include <iostream>
#include <sstream>
#include <math.h>
#include <map>
#include <algorithm>
#include <vector>
#include <iterator>
#include <signal.h>

using namespace std;

bool com_op = false;                                 // Set to true for command line output

#ifdef USE_MPI
#include "mpi.hh"
#endif

#include "model.hh"
#include "output.hh"
#include "input.hh"
#include "const.hh"
#include "utils.hh"
#include "simulate.hh"
#include "post_sim.hh"
#include "data_sim.hh"
#include "mcmc.hh"
#include "abc.hh"
#include "abc_smc.hh"
#include "mfa.hh"
#include "pas.hh"
#include "extend.hh"
#include "mpi.hh"
#include "lzw.hh"
#include "validation.hh"

vector <BICITag> get_tags(vector <string> &sec, Operation &mode, ExtFactor &ext_factor, string &file, vector <string> &data_sim_lines, bool &no_question, bool &test, string &scan_info);

int main(int argc, char** argv)
{	
#ifdef USE_MPI                            // This is for the parallel version of the code 
  //MPI_Init(&argc,&argv);       	
	 MPI_Init(NULL, NULL);
#endif
	
	vector <string> sec;
	for(auto i = 1u; i < (unsigned int)argc; i++) sec.push_back(argv[i]);
	
	//string te = "α this"; encode(te); return 0;

	//test(); return 0;

	auto total_time = clock();
	
	print_diag("start bici");
	
	print_diag("start sum");
	
	init_log_sum();
	
	//mvn_prior_check(); return 0;
	//mvn_jeffreys_check(); return 0; 
	//test_jeffreys(); return 0;
	//solve_cubic_spline_test(); return 0;
	//test_distribution(); return 0; 
	//test_incomplete_distribution(); return 0;
	//generate_data();  return 0; 
	//if(argc > 2){ simulate_trans_exp(); return 0;}
	
	auto core_spec = UNSET;
	auto seed = UNSET;
	auto sim_sel = 1u;
	auto mode = MODE_UNSET;
	
	ExtFactor ext_factor;
	string file;
	vector <string> data_sim_lines;
	auto no_question = false;
	auto test = false;
	string scan_info;
	
	auto tags = get_tags(sec,mode,ext_factor,file,data_sim_lines,no_question,test,scan_info);
	
	print_diag("start model");
	
	switch(mode){
	case TORNADO_SETUP:
		{
			Validation valid;
			valid.tornado_setup(mode,ext_factor,file,data_sim_lines,test);
		}
		return 0;
		
	case TORNADO_RESULT:
		{
			Validation valid;
			valid.tornado_result(mode,ext_factor,no_question,file);
		}
		return 0;
		
	case SCAN_SETUP:
		{
			Validation valid;
			valid.scan_setup(scan_info,mode,ext_factor,file,data_sim_lines,test);
		}
		return 0;
		
	case SCAN_RESULT:
		{
			Validation valid;
			valid.scan_result(scan_info,mode,ext_factor,no_question,file);
		}
		return 0;
	
	default: break;
	}
	
	Model model(mode,ext_factor,no_question);     // Used to store the model structure 
	
	model.samp_type = ALL_SAMP;
	
	for(const auto &tag : tags){
		switch(tag.type){
		case CORE: core_spec = tag.value; break;  
		case SEED: seed = tag.value; break;  
		case SIM_SEL: 
			if(model.mode != DATA_SIM) alert_input("The 'sim' tag can only be used with 'data-sim'.");
			sim_sel = tag.value; 
			break;  
		case OP: file = tag.te; break;
		case SAMP_TYPE:
			switch(tag.value){
			case 0: model.samp_type = LOCAL_SAMP; break;
			case 1: model.samp_type = SAMP_SAMP; break;
			case 2: model.samp_type = SIM_SAMP; break;
			case 3: model.samp_type = SIM_CL_SAMP; break;
			}
			break;
		case TAG_UNSET: alert_input("Tag must be set"); break; 
		}
	}
	
	if(true){
		switch(model.samp_type){
		case LOCAL_SAMP: cout << "Local samp" << endl; break;
		case SAMP_SAMP: cout << "Samp samp" << endl; break;
		case SIM_SAMP: cout << "Sim samp" << endl; break;
		case SIM_CL_SAMP: cout << "Sim cl samp" << endl; break;
		case ALL_SAMP: break;
		}
	}
	
	Mpi mpi(core_spec,model);               // Sets up mpi

	print_diag("start input");

	Output output(model,mpi);               // Sets up the class for model outputs

	{
		Input input(model,file,seed,mpi);     // Imports information from input into model
		output.init(input);
	}
	
	if(false && op()){
		cout << " turn off data summary" << endl;
		output.summary(model);                // Outputs a text file sumarising the model

		output.data_summary(model);           // Outputs a text file sumarising the data
	}

	switch(model.mode){
	case SIM:                               // Simulates from the model
		{
			Simulate simu(model,output,mpi);
			simu.run();
		}
		break;
		
	case INF:                               // Performs inference on the model
		{
			switch(model.details.algorithm){
			case PAS_MCMC:	
				{
					PAS pas(model,output,mpi);
					pas.run();
				}
				break;
				
			case MFA_ALG:
				{
					MFA mfa(model,output);
					mfa.run();
				}
				break;
				
			case DA_MCMC:
				{
					MCMC mcmc(model,output,mpi);
					mcmc.run();
				}
				break;
				
			case ABC_ALG:
				{
					ABC abc(model,output,mpi);
					abc.run();
				}
				break;
				
			case ABC_SMC_ALG:
				{
					ABC_SMC abcsmc(model,output,mpi);
					abcsmc.run();
				}
				break;
				
			default: 
				alert_input("This algorithm has not been implemented yet");
				return 0;
			}
		}
		break;
	
	case PPC:                               // Simulates from the posterior
		{
			PostSim post_simu(model,output,mpi);
			post_simu.run();
		}
		break;
		
	case EXT:                               // Simulates from the posterior
		{
			Extend extend(model,output,mpi);
			extend.run();
		}
		break;
		
	case DATA_SIM:                          // Simulates some data and adds to 
		{
			DataSim data_sim(model,output);
			for(const auto &va : data_sim_lines){
				data_sim.run(va,sim_sel);
			}
		}
		break;
		
	case DATA_SHOW:                         // Shows all the data sources
		{
			DataSim data_sim(model,output);
			data_sim.show();
		}
		break;
		
	case DATA_DEL:                          // Deletes data sources
		{
			DataSim data_sim(model,output);
			for(const auto &ds : data_sim_lines){
				data_sim.del(ds);
			}
		}
		break;
		
	case DATA_CLEAR:                          // Clear sim/inf/post-sim results
		{
			DataSim data_sim(model,output);
			for(const auto &ds : data_sim_lines){
				data_sim.clear(ds);
			}
		}
		break;
		
	case COMPRESS:
		{
			DataSim data_sim(model,output);
			data_sim.compress();
		}
		break;
		
	case DECOMPRESS:
		{
			DataSim data_sim(model,output);
			data_sim.decompress();
		}
		break;
		
	case TORNADO_RESULT: case TORNADO_SETUP: 
	case SCAN_RESULT: case SCAN_SETUP: 
		emsg("Should not be here");
		break;

	case MODE_UNSET:
		alert_input("The mode is unset");
		break;
	}
	
	if(model.mode != DATA_SHOW){	
		auto total_cpu = clock()-total_time;
	
		auto op_time = clock();
		output.end(file,total_cpu);
		auto op_cpu = clock()-op_time;
		
		if(!model.data_mode()){
			if(op() && !com_op) output.final_time(total_cpu+op_cpu,op_cpu);
		
			if(!com_op) output.final_memory_usage();
		}
	}
	
#ifdef USE_MPI
	MPI_Finalize();
#endif
}


/// Gets values for tags when BICI is run
vector <BICITag> get_tags(vector <string> &sec, Operation &mode, ExtFactor &ext_factor, string &file, vector <string> &data_sim_lines, bool &no_question, bool &test, string &scan_info)
{
	vector <BICITag> tags;
	
	// Extracts tags first
	
	auto i = 0u;
	while(i < sec.size()){         // Removes tags first
		auto te = sec[i];
		
		if(begin_str(te,"-")){
			BICITag tag; tag.type = TAG_UNSET; tag.processed = false;
		
			if(te == "-no-question"){
				no_question = true;
			}
			else{
				if(te == "-test"){
					test = true;
				}
				else{
					auto spl = split(te.substr(1),'=');
			
					if(spl.size() != 2) alert_input("Tag '"+te+"' is not recognised");
					
					if(spl[0] == "co" || spl[0] == "core"){
						tag.type = CORE;
						tag.value = number(spl[1]);
					}
					
					if(spl[0] == "seed"){
						tag.type = SEED;
						tag.value = number(spl[1])+1;
					}
					
					if(spl[0] == "op" || spl[0] == "output"){
						tag.type = OP;
						tag.te = spl[1];
					}
				
					if(spl[0] == "sim" || spl[0] == "simulation"){
						tag.type = SIM_SEL;
						tag.value = number(spl[1]);
					}
					
					
					if(spl[0] == "samp"){
						tag.type = SAMP_TYPE;
						tag.value = number(spl[1]);
					}
					
					if(tag.type == TAG_UNSET){
						alert_input("Tag '"+te+"' is not recognised");
					}
				
					tags.push_back(tag);
				}
			}
			
			sec.erase(sec.begin()+i);
		}
		else{
			i++;
		}
	}
	
	auto N = sec.size();
	for(auto i = 0u; i < N; i++){  // Processes non-tags
		auto te = sec[i];
		
		auto mode_new = MODE_UNSET;
		string file_new = "";
		if(te == "sim" || te == "simulate") mode_new = SIM;
	
		if(te == "inf" || te == "inference") mode_new = INF;

		if(te == "post-sim" || te == "posterior-simulation") mode_new = PPC;
			
		if(te == "ext" || te == "extend"){
			mode_new = EXT;
			
			if(i+1 == N) alert_input("There must be a value after 'ext'.");
			else{
				auto spl = split(sec[i+1],'%');
				auto val = number(spl[0]);
			
				if(val == UNSET){
					alert_input("'"+spl[0]+"' must be a number.");
				}
				else{
					ext_factor.value = val;
				
					if(spl.size() == 1){
						ext_factor.percent = false;
					}
					else{
						ext_factor.percent = true;
						
						if(val <= 100){
							alert_input("'"+te+"' must have a percentage greater than 100%.");	
						}
					}								
				}
			}
			i++;
		}
		
		if(te == "data-sim"){
			mode_new = DATA_SIM;
			if(i+1 == N) alert_input("There must be an expression after 'data-sim'.");
			else{
				i++; 
				while(i < N){ data_sim_lines.push_back(sec[i]); i++;}
			}
		} 
		
		if(te == "data-show") mode_new = DATA_SHOW;
		
		if(te == "data-del"){
			mode_new = DATA_DEL;
			if(i+1 == N) alert_input("There must be an expression after 'data-del'.");
			else{ 
				i++;
				while(i < N){ data_sim_lines.push_back(sec[i]); i++;}
			}
		}
		
		if(te == "clear"){
			mode_new = DATA_CLEAR;
			if(i+1 == N) alert_input("There must be an expression after 'clear'.");
			else{
				i++;
				while(i < N){ data_sim_lines.push_back(sec[i]); i++;}
			}
		}
		
		if(te == "tornado"){
			mode_new = TORNADO_SETUP;
			if(i+1 == N) alert_input("There must be an expression (or expressions) after 'torndo' that specify the data.");
			else{
				i++;
				while(i < N){ data_sim_lines.push_back(sec[i]); i++;}
			}
		}
			
		if(te == "tornado-res") mode_new = TORNADO_RESULT;
	
		if(te == "compress" || te == "comp"){
			mode_new = COMPRESS;
		}
		
		if(te == "decompress" || te == "decomp"){
			mode_new = DECOMPRESS;
		}
		
		if(begin_str(te,"scan:")){
			mode_new = SCAN_SETUP;
			if(i+1 == N) alert_input("There must be an expression (or expressions) after 'scan' that specify the data.");
			else{
				i++;
				while(i < N){ data_sim_lines.push_back(sec[i]); i++;}
			}
			scan_info = te.substr(5);
		}
		
		if(begin_str(te,"scan-res:")){
			mode_new = SCAN_RESULT;
			scan_info = te.substr(9);
		}
		
		if(end_str(te,".bici")){
			file_new = te;
		}
		
		if(mode_new == MODE_UNSET && file_new == ""){
			alert_input("Tag '"+te+"' is not recognised");
		}
		else{	
			if(mode_new != MODE_UNSET){
				if(mode != MODE_UNSET){
					alert_input("Cannot set multiple operations (sim/inf/post-sim/ext/data-sim/data-show/data-del/clear)");
				}
				else{
					mode = mode_new;
				}
			}
			
			if(file_new != ""){
				if(file != ""){
					alert_input("Cannot define multiple '.bici' files");
				}
				else{
					file = file_new;
				}
			}
		}
	}
	
	if(mode == MODE_UNSET){
		alert_input("A BICI operation must be specified to tell BICI what to do (sim/inf/post-sim/ext/data-sim/data-show/data-del)");
	}
	
	if(file == ""){
		alert_input("A '.bici' file must be specified.");
	}
	
	if(file == "default.bici"){
		file = default_file;
		com_op = true;
	}
	
	return tags;
}

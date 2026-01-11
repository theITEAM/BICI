// This is start code for BICI

// BICI (Bayesian individual-based compartmental inference)
// © C. M. Pooley and G. Marion

// Compilation: make
// Run: ./bici 

// ssh gaia.bioss.ac.uk  

// tar -xzf foo.tgz

// Load mpi: module load mpi/openmpi-x86_64
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

// When running in parallel then load mpi using: module load mpi/openmpi-x86_64

// Test: valgrind --exit-on-first-error=yes --error-exitcode=1 --leak-check=yes -s ./bici-para Execute/init.bici sim
// Test: valgrind --exit-on-first-error=yes --error-exitcode=1 --leak-check=yes -s ./bici-para Execute/init.bici inf
// Test:  valgrind --exit-on-first-error=yes --error-exitcode=1 --leak-check=yes -s  ./bici-para simple_model.bici inf

// Compiling for windows: .\make.bat from the directory D:\C_complier_BICI2.0\bin

// 45541 lines of code (15/11/24)
// 51000 lines of code (18/12/24)
// 56000 lines of code (16/05/25)
// 65634 lines of code (17/09/25)

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
#include "mcmc.hh"
#include "abc.hh"
#include "abc_smc.hh"
#include "mfa.hh"
#include "pas.hh"
#include "extend.hh"
#include "mpi.hh"

vector <BICITag> get_tags(int argc, char** argv, Operation &mode, ExtFactor &ext_factor,string &file);

int main(int argc, char** argv)
{	
	auto total_time = clock();
	
	//mvn_jeffreys_check(); return 0; 
	
#ifdef USE_MPI                            // This is for the parallel version of the code 
  MPI_Init(&argc,&argv);                 
#endif

	init_log_sum();
	
	//test_jeffreys(); return 0;
	//solve_cubic_spline_test(); return 0;
	//test_distribution(); return 0; 
	//test_incomplete_distribution(); return 0;
	//generate_data();  return 0; 
	//if(argc > 2){ simulate_trans_exp(); return 0;}
	
	auto core_spec = UNSET;
	auto seed = UNSET;
	auto mode = MODE_UNSET;
	
	ExtFactor ext_factor;
	string file="";
	
	auto tags = get_tags(argc,argv,mode,ext_factor,file);
	
	Model model(mode,ext_factor);       // Used to store the model structure 
	
	model.samp_type = ALL_SAMP;
	
	for(const auto &tag : tags){
		switch(tag.type){
		case CORE: core_spec = tag.value; break;  
		case SEED: seed = tag.value; break;  
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

	Input input(model,file,seed,mpi);// Imports information from input file into model

	Output output(model,input,mpi);         // Sets up the class for model outputs

	if(false && op()){
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
		
	case MODE_UNSET:
		alert_input("The mode is unset");
		break;
	}
	
	auto total_cpu = (clock()-total_time)/CLOCKS_PER_SEC;
	
	auto op_time = clock();
	output.end(file,total_cpu);
	auto op_cpu = (clock()-op_time)/CLOCKS_PER_SEC;
	
	if(op() && !com_op) output.final_time(total_cpu+op_cpu,op_cpu);
	
	if(!com_op) output.final_memory_usage();
	
#ifdef USE_MPI
	MPI_Finalize();
#endif
}


/// Gets values for tags when BICI is run
vector <BICITag> get_tags(int argc, char** argv, Operation &mode, ExtFactor &ext_factor, string &file)
{
	vector <BICITag> tags;
	
	auto N = (unsigned int)argc;
	for(auto i = 1u; i < N; i++){
		string te = argv[i];
		if(te.substr(0,1) != "-"){
			auto mode_new = MODE_UNSET;
			string file_new = "";
			if(te == "sim" || te == "simulate") mode_new = SIM;
			else{
				if(te == "inf" || te == "inference") mode_new = INF;
				else{
					if(te == "post-sim" || te == "posterior-simulation") mode_new = PPC;
					else{
						if(te == "ext" || te == "extend"){
							mode_new = EXT;
							
							if(i+1 == N) alert_input("There must be a value after 'ext'.");
							else{
								auto spl = split(argv[i+1],'%');
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
						else{
							if(end_str(te,".bici")){
								file_new = te;
							}
						}
					}
				}
			}
			
			if(mode_new == MODE_UNSET && file_new == ""){
				alert_input("Tag '"+te+"' is not recognised");
			}
			else{	
				if(mode_new != MODE_UNSET){
					if(mode != MODE_UNSET){
						alert_input("Cannot set multiple 'sim', 'inf' or 'post-sim'");
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
		else{
			BICITag tag; tag.type = TAG_UNSET; tag.processed = false;
		
			auto spl = split(te.substr(1),'=');
	
			if(spl.size() != 2) alert_input("Tag '"+te+"' is not recognised");
			
			if(spl[0] == "co" || spl[0] == "core"){
				tag.type = CORE;
				tag.value = number(spl[1]);
			}
			
			if(spl[0] == "seed"){
				tag.type = SEED;
				tag.value = number(spl[1]);
			}
			
			if(spl[0] == "op" || spl[0] == "output"){
				tag.type = OP;
				tag.te = spl[1];
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
	
	if(mode == MODE_UNSET){
		alert_input("Either 'sim', 'inf', 'post-sim' or 'ext' must be specified to tell BICI what to do.");
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

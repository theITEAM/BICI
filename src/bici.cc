// This is start code for BICI

// BICI (Bayesian individual-based compartmental inference)
// © C. M. Pooley and G. Marion

// Compilation: make
// Run: ./bici 

// Load mpi: module load mpi/openmpi-x86_64
// mpirun -n 3 ./bici 

// ssh cpooley@azog.bioss.ac.uk

// Different commands 
// -ch=0 / -chain=0                           Runs the zeroth chain
// -nch=3 / -nchain=3                         Runs multiple chains 
// -op=init.bici / -output=init.bici          Sets the output file
// -samp=0                                    Sets how sampling is done (for diagnostic purposes)

// When running in parallel then load mpi using: module load mpi/openmpi-x86_64

// Test: valgrind --leak-check=yes -s ./bici

// Compiling for windows: .\make.bat from the directory D:\C_complier_BICI2.0\bin

// 45541 lines of code (15/11/24)
// 51000 lines of code (18/12/25)

#include <iostream>
#include <sstream>
#include <math.h>
#include <map>
#include <algorithm>
#include <vector>
#include <iterator>
#include <signal.h>

using namespace std;

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
#include "mpi.hh"

vector <BICITag> get_tags(int argc, char** argv);

int main(int argc, char** argv)
{	
#ifdef USE_MPI                            // This is for the parallel version of the code 
  MPI_Init(&argc,&argv);                 
#endif
	
	init_log_sum();
	
	//test_distribution(); return 0; 
	//test_incomplete_distribution(); return 0;
	//generate_data();  return 0; 
	//if(argc > 2){ simulate_trans_exp(); return 0;}

	Model model;                            // Used to store the model structure 
	
	model.samp_type = ALL_SAMP;
	
	string file = "Execute/init.bici";
	
	auto chain = 0;
	auto nchain = 1;
	
	auto tags = get_tags(argc,argv);
	
	for(const auto &tag : tags){
		switch(tag.type){
		case CHAIN: chain = tag.value; break;  
		case NCHAIN: nchain = tag.value; break; 
		case OP: file = tag.te; break;
		case SAMP_TYPE:
			switch(tag.value){
			case 0: model.samp_type = LOCAL_SAMP; break;
			case 1: model.samp_type = SAMP_SAMP; break;
			case 2: model.samp_type = SIM_SAMP; break;
			case 3: model.samp_type = SIM_CL_SAMP; break;
			}
			break;
		case TAG_UNSET: emsg("Tag must be set"); break; 
		}
	}
	
	cout << "chain = " << chain << endl;
	
	switch(model.samp_type){
	case LOCAL_SAMP: cout << "Local samp" << endl; break;
	case SAMP_SAMP: cout << "Samp samp" << endl; break;
	case SIM_SAMP: cout << "Sim samp" << endl; break;
	case SIM_CL_SAMP: cout << "Sim cl samp" << endl; break;
	case ALL_SAMP: cout << "All samp" << endl; break;
	}
	
	Mpi mpi(model);                         // Sets up mpi

#ifdef USE_MPI
	chain = mpi.core;
#endif
	
	Input input(model,file,chain);           // Imports information from input file into model

	Output output(chain,model,input,mpi);   // Sets up the class for model outputs

	if(op()){
		output.summary(model);                // Outputs a text file sumarising the model

		output.data_summary(model);           // Outputs a text file sumarising the data
	}
	
	switch(model.mode){
	case SIM:                               // Simulates from the model
		{
			Simulate simu(model,output);
			simu.run();
		}
		break;
		
	case INF:                               // Performs inference on the model
		{
			//if(model.details.algorithm == DA_MCMC){
			//model.details.algorithm = MFA_ALG;
			//}

			switch(model.details.algorithm){
			case PAS_MCMC:	
				{
					PAS pas(model,mpi,output);
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
					MCMC mcmc(model,output);
					mcmc.run();
				}
				break;
				
			case ABC_ALG:
				{
					ABC abc(model,output);
					abc.run();
				}
				break;
				
			case ABC_SMC_ALG:
				{
					ABC_SMC abcsmc(model,output);
					abcsmc.run();
				}
				break;
				
			default: 
				emsg("This algorithm has not been implemented yet");
				return 0;
			}
		}
		break;
	
	case PPC:                               // Simulates from the posterior
		{
			PostSim post_simu(model,output);
			post_simu.run();
		}
		break;
		
	case MODE_UNSET:
		emsg("The mode is unset");
		break;
	}
	
	output.generate_files();
	
	output.updated_file(file);

#ifdef USE_MPI
	MPI_Finalize();
#endif
}


/// Gets values for tags when BICI is run
vector <BICITag> get_tags(int argc, char** argv)
{
	vector <BICITag> tags;
	
	for(auto i = 1u; i < (unsigned int)argc; i++){
		string te = argv[i];
		if(te.substr(0,1) != "-") emsg("Tag '"+te+"' is not recognised");
		auto spl = split(te.substr(1),'=');
		if(spl.size() != 2) emsg("Tag '"+te+"' is not recognised");
		
		BICITag tag; tag.type = TAG_UNSET; tag.processed = false;
		
		if(spl[0] == "nch" || spl[0] == "nchain"){
			tag.type = NCHAIN;
			tag.value = number(spl[1]);
		}
		
		if(spl[0] == "ch" || spl[0] == "chain"){
			tag.type = CHAIN;
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
			emsg("Tag '"+te+"' is not recognised");
		}
	
		tags.push_back(tag);
	}
	
	return tags;
}

// Implements an MCMC algorithm

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm> 
 
using namespace std;

#include "mcmc.hh"
#include "state.hh"
#include "utils.hh"
#include "matrix.hh"
#include "synchronise.hh"

MCMC::MCMC(const Model &model, Output &output, Mpi &mpi) : model(model), output(output), mpi(mpi)
{	
	const auto &de = model.details;
	
	if(de.anneal_type == ANNEAL_SCAN){
		nburnin = LARGE;
		nsample = LARGE;
	}
	else{
		auto f = de.burnin_frac/100.0;
		if(true){
			nburnin = f*de.sample;
			nsample = de.sample;
		}
		else{
			nburnin = (f/(1-f))*de.sample;
			nsample = de.sample + nburnin;
		}
	}
	
	num_per_core = de.num_per_core;
	output_param = de.output_param;
	output_state = de.output_state;

	nparam_op = 0;
	nstate_op = 0;
	
	for(auto ch = 0u; ch < num_per_core; ch++){
		chain.push_back(Chain(nburnin,nsample,model,output));
	}

	for(auto &ch : chain) ch.state.init();
}
 
 
/// Runs MCMC
void MCMC::run()
{
	auto time_init_start = clock();
	
	percentage_start(INIT_PER);
	
	for(auto ch = 0u; ch < chain.size(); ch++){
		chain[ch].init(ch,chain.size());
	}

	percentage_end();
	
	auto time_start = clock();
	
	percentage_start(RUN_PER);
	
	//chain[0].state.scan_param();
 	for(auto s = 0u; s < nsample; s++){
		percentage(s,nsample);
	
		if(model.sync_on) synchronise_proposal(s,chain,mpi);

		for(auto &ch : chain){
			ch.burn_update(s);

			ch.update(s);
		}		

		sample_op(s);	

#ifdef USE_MPI
		mpi.sample_barrier(s,nsample);
#endif
	}
	
	//for(auto &ch : chain) ch.cor_matrix.check();
	
	if(false){
		auto tot = chain[0].state.timer[TEMP1];
		cout << cpu_percent(chain[0].state.timer[TEMP2],tot) << "mbp" << endl;
		cout << cpu_percent(chain[0].state.timer[TEMP3],tot) << "update param" << endl;
		cout << cpu_percent(chain[0].state.species[0].timer[SSP_TEMP1],tot) << "markov" << endl;
		cout << cpu_percent(chain[0].state.species[0].timer[SSP_TEMP2],tot) << "mar calc" << endl;
	}
	
#ifdef USE_MPI
	mpi.barrier();
#endif

	output.set_inference_prop(double(100.0*nburnin)/nsample,"burnin-percent",BURNIN_FRAC_DEFAULT);
	
	double time_total = (clock()-time_start)/num_per_core;
	double init_time = (time_start-time_init_start)/num_per_core;
	for(auto ch = 0u; ch < num_per_core; ch++){
		auto ch_tot = mpi.core*num_per_core+ch;
		
		const auto &cha = chain[ch];
		auto diag = cha.diagnostics(time_total,init_time);
		output.set_diagnostics(ch_tot,diag);
		output.terminal_info.push_back(cha.get_terminal_info(ch_tot));
	}
	
	percentage_end();
	
	if(profiling) chain[0].state.profile_memory();
}


/// Determines if a state output is made
void MCMC::sample_op(unsigned int s)
{
	auto base = mpi.core*num_per_core;
	
	if(nparam_op <= double(s*output_param)/nsample){
		for(auto ch = 0u; ch < num_per_core; ch++){
			output.param_sample(s,base+ch,chain[ch].state);
		}	
		nparam_op++;
	}
	
	if(s >= nburnin){
		if(nstate_op <= double((s-nburnin)*output_state)/(nsample-nburnin)){
			for(auto ch = 0u; ch < num_per_core; ch++){
				output.state_sample(s,base+ch,chain[ch].state);
			}
			nstate_op++;
		}
	}
}

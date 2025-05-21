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
	percentage_start(INIT_PER);

	for(auto ch = 0u; ch < chain.size(); ch++){
		chain[ch].init(ch,chain.size());
	}
	percentage_end();
	
	percentage_start(RUN_PER);
	
	long time_start = clock();

	for(auto s = 0u; s < nsample; s++){
		//if(mpi.core == 0) cout << s << " samp" << endl;
		percentage(s,nsample);
		
		for(auto &ch : chain){
			ch.burn_update(s);
			ch.update(s);
		}		
		
		sample_op(s);
	}
	
	output.set_output_burnin(double(100.0*nburnin)/nsample);
	
	auto time_total = (clock()-time_start)/num_per_core;
	for(auto ch = 0u; ch < num_per_core; ch++){
		auto diag = chain[ch].diagnostics(time_total);
		output.set_diagnostics(mpi.core*num_per_core+ch,diag);
	}
	
	percentage_end();
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

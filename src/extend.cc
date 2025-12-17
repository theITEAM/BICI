// Simulates using the posterior distribution

#include <iostream>
#include <fstream>
#include <algorithm>
#include <assert.h>
#include <math.h>

using namespace std;

#include "extend.hh"
#include "state.hh"
#include "utils.hh"
#include "mpi.hh"

/// Initilaises the simulation
Extend::Extend(const Model &model, Output &output, Mpi &mpi) : model(model), output(output), mpi(mpi)
{	
	const auto &de = model.details;
	
	double f;

	const auto &ef = model.ext_factor;
	if(ef.percent) f = (ef.value-100)/100.0;
	else f = (ef.value-de.sample)/de.sample;
		
	nburnin = 0u;
	
	nsample = (unsigned int)(de.sample*f); 
	
	num_per_core = de.num_per_core;
	output_param = (unsigned int)(de.output_param*f);
	output_state = (unsigned int)(de.output_state*f);

	nparam_op = 0;
	nstate_op = 0;
	
	for(auto ch = 0u; ch < num_per_core; ch++){
		chain.push_back(Chain(nburnin,nsample,model,output));
	}

	for(auto &ch : chain) ch.state.init();
}


/// Performs a simulation
void Extend::run()
{	
	percentage_start(INIT_PER);
	
	// Copies loaded states
	auto first_samp = UNSET;
	if(op()){
		const auto &de = model.details;
		
		auto burn = double((nsample+de.sample)*de.burnin_frac)/100;
		
		State state(model);
		
		state.init();
		
		auto smax = model.sample.size();
		for(auto s = 0u; s < smax; s++){
			percentage(s,smax);
		
			const auto &samp = model.sample[s];	
			if(samp.num >= burn){
				auto param_val = model.post_param(samp);
				
				state.load_samp(param_val,samp);

				if(first_samp == UNSET) first_samp = samp.num;
				output.state_sample(samp.num,samp.ch,state);
			}
		}		
	}
	
	// Intialises chains
	for(auto ch = 0u; ch < chain.size(); ch++){
		auto ch_tot = mpi.core*num_per_core+ch;
		chain[ch].init(ch_tot,chain.size());
	}

	// Loads up proposal infomation
	for(auto ch = 0u; ch < chain.size(); ch++){
		auto &cha = chain[ch];
		auto ch_tot = mpi.core*num_per_core+ch;
		
		auto &cor_m = cha.cor_matrix;
		const auto &termi = model.terminal_info;
		auto i = 0u; while(i < termi.size() && termi[i].ch != ch_tot) i++;
		if(i == termi.size()) emsg("Cannot fint ti");
		
		const auto &ti = termi[i];
		cor_m.n = ti.n;
		cor_m.n_start = ti.n_start;
		cor_m.av = ti.av;
		cor_m.av2 = ti.av2;
		
		cha.join_proposal_update();
		cha.set_prop_info(ti.prop_info_store);
		
		for(auto &pro : cha.proposal) pro.update_sampler(cha.cor_matrix);
	}
	
	percentage_end();
	
	percentage_start(RUN_PER);
	
	long time_start = clock();
	
	for(auto s = 0u; s < nsample; s++){
		percentage(s,nsample);
		
		for(auto &ch : chain){
			ch.burn_update(s);
			ch.update(s);
		}		
		
		sample_op(s);

#ifdef USE_MPI
		mpi.sample_barrier(s,nsample);
#endif
	}
	
#ifdef USE_MPI
	mpi.barrier();
#endif

	// Updates quantities in inference	
	const auto &de = model.details;
	//auto frac_new = double(100.0*first_samp/(nsample+de.sample));
	//output.set_inference_prop(frac_new,"burnin-percent",BURNIN_FRAC_DEFAULT);
	output.set_inference_prop(nsample+de.sample,"update",MCMC_SAMPLE_DEFAULT);
	output.set_inference_prop(output_param+de.output_param,"param-output",MCMC_OP_PARAM_DEFAULT);
	output.set_inference_prop(output_state+de.output_state,"state-output",MCMC_OP_STATE_DEFAULT);
	
	//output.terminal_info = chain[0].get_terminal_info();
	
	double time_total = (clock()-time_start)/num_per_core;
	for(auto ch = 0u; ch < num_per_core; ch++){
		const auto &cha = chain[ch];
		auto diag = cha.diagnostics(time_total);
		output.set_diagnostics(mpi.core*num_per_core+ch,diag);
	}
	
	percentage_end();
}


/// Determines if a state output is made
void Extend::sample_op(unsigned int s)
{
	auto base = mpi.core*num_per_core;
	auto last = model.details.sample;
	
	if(nparam_op <= double(s*output_param)/nsample){
		for(auto ch = 0u; ch < num_per_core; ch++){
			output.param_sample(last+s,base+ch,chain[ch].state);
		}	
		nparam_op++;
	}
	
	if(s >= nburnin){
		if(nstate_op <= double((s-nburnin)*output_state)/(nsample-nburnin)){
			for(auto ch = 0u; ch < num_per_core; ch++){
				output.state_sample(last+s,base+ch,chain[ch].state);
			}
			nstate_op++;
		}
	}
}

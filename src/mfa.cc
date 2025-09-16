// Implements an MFA algorithm

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm> 
 
using namespace std;

#include "mfa.hh"
#include "state.hh"
#include "utils.hh"

MFA::MFA(const Model &model, Output &output) : model(model), output(output), state(model)
{	
	const auto &de = model.details;
	
	auto f = de.burnin_frac/100.0;
	nburnin = f*de.sample;
	nsample = de.sample;
}
 
 
/// Runs MFA
void MFA::run()
{
	print("Run...");
	
	auto param_val = model.param_sample();
	
	model.print_param(param_val);
	
	auto Li = state.like_integrate(param_val);
	
	cout << Li << endl;
	emsg("don");
	
	/*
	param_prop_init();
			
	if(false) output.print_param(state.param_val);
	
	Particle part;
	
	auto like_obs_max = -LARGE;
	
	print("Finding initial state");

	auto loop_max = 10u; 
	for(auto loop = 0u; loop < loop_max; loop++){ 
		print(to_string(loop)+" Initial state");
		
		auto param_val = model.param_sample();
		auto initc_val = model.initc_sample(param_val);
		
		if(false){
			output.print_initc(initc_val);
			output.print_param(param_val);
		}
		
		print("simulate");
		
		state.simulate(param_val,initc_val);
		
		print("simulate done");
	
		state.check("First check");
	
		print("Checked");
		
		state.resample_ind();   // Resamples individual such that fixed events become correct
		
		print("resampled");
		
		burn_info.add_L(state.like); 
		
		if(state.like.obs > like_obs_max){
			like_obs_max = state.like.obs;
			part = state.generate_particle();
		}
		
		state.check("before");
	}

	state.set_particle(part);
	
	state.check("before");
	
	print("Starting sampling");

	long time_start = clock();

	auto op_step = get_op_step(nsample);
	
	for(auto s = 0u; s < nsample; s++){
		if(s%100 == 0 && com_op) progress(s,nsample);		
		else{
			if(s%op_step == 0) cout << s << " " << nsample << " sample" << endl;
		}
	
		burn_info.add_L(state.like); 
		burn_info.setup(s,nburnin,nsample,model.details);
		
		if(param_op(s)) output.param_sample(s,state);
		
		if(state_op(s)) output.state_sample(s,state);
		
		if(burn_info.on){
			if(s%10 == 0) state.update_individual_sampler();
			
			for(auto &pro : proposal){
				if(pro.param_list.size() > 0){
					pro.add_sample(state.get_param_val_prop();
					if(s%10 == 0) pro.update_sampler();
				}
			}
	
			//if(adapt_prop_prob && s%UPDATE_PROP_PROB == UPDATE_PROP_PROB-1) set_proposal_prob();
		}
		
		update(s);
	
		if(s%100 == 0) state.check("periodic");
	}
	
	output.set_output_burnin(double(100.0*nburnin)/nsample);
	
	diagnostics(clock()-time_start);
	
	if(false){
		state.output_gen_dif();
		state.output_waifw("waifw2.txt");
	}
	*/
}


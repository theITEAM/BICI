// Simulates using the posterior distribution

#include <iostream>
#include <fstream>
#include <algorithm>
#include <assert.h>
#include <math.h>

using namespace std;

#include "post_sim.hh"
#include "state.hh"
#include "utils.hh"
#include "mpi.hh"

/// Initilaises the simulation
PostSim::PostSim(const Model &model, Output &output, Mpi &mpi) : model(model), output(output), mpi(mpi), state(model)
{	
	state.init();
}


/// Performs a simulation
void PostSim::run()
{
	percentage_start(RUN_PER);
	
	const auto &details = model.details;

	if(model.sample.size()==0) emsg("Sample number should not be zero");
	
	auto smax = details.num_per_core;

	// Precalculate parameter values
	auto param_av = calculate_param_av();
	
	for(auto s = 0u; s < smax; s++){
		percentage(s,smax);
		
		const auto &samp = model.sample[(unsigned int)(model.sample.size()*ran())];
		//const auto &samp = model.sample[0];

		PV param_val;
		if(model.details.ps_type == PS_SAMPLE) param_val = model.post_param(samp);
		else param_val = param_av;
	
		//model.print_param(param_val);
		
		if(false) model.print_param(param_val);

		if(model.details.param_only){
			auto ie_store = model.post_ie_store(samp);
			auto initc_val = model.initc_sample(param_val);
			state.simulate(param_val,initc_val,ie_store);
		}
		else{
			state.post_sim(param_val,samp);
		}
		
		state.check("Check state");
	
		output.param_sample(s,0,state);
		output.state_sample(s,0,state);
		
#ifdef USE_MPI
		mpi.sample_barrier(s,smax);
#endif
	}
	
#ifdef USE_MPI
	mpi.barrier();
#endif
	
	percentage_end();
}


/// Calculates an average value for 
PV PostSim::calculate_param_av() const
{
	PV pv;
	if(model.details.ps_type == PS_SAMPLE) return pv;
	
	auto N = model.param_vec.size();
	vector < vector <double> > value_arr;
	value_arr.resize(N);

	for(auto s = 0u; s < model.sample.size(); s++){			
		const auto &samp = model.sample[s];		
		auto param_val = model.post_param(samp);

		for(auto th = 0u; th < N; th++){
			const auto &pv = model.param_vec[th];	
			const auto &par = model.param[pv.th];
			if(par.variety == DIST_PARAM || par.variety == PRIOR_PARAM){
				value_arr[th].push_back(param_val.value[th]);
			}
		}
	}			
	
	vector <double> value_av(N,UNSET); 
	for(auto th = 0u; th < N; th++){
		auto M = value_arr[th].size();
		if(M > 0){
			switch(model.details.ps_type){
			case PS_MEAN:
				value_av[th] = sum(value_arr[th])/M;
				break;
				
			case PS_MEDIAN:
				value_av[th] = median(value_arr[th]);
				break;
			
			case PS_SAMPLE: emsg("Should not be"); break;
			}
		}
	}

	return model.post_param_av(value_av);
}
	
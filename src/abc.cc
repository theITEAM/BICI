/// Implements the basic ABC rejection algorithm

#include <iostream>
#include <fstream>
#include <algorithm> 
#include <sstream> 

using namespace std;

#include "abc.hh"
#include "utils.hh"

ABC::ABC(const Model &model, Output &output, Mpi &mpi) : model(model), output(output), mpi(mpi), state(model)
{	
	const auto &de = model.details;
	nsample = de.sample;
	acc_frac = de.accfrac;
	state.init();
}
 

/// Implements a version of simple ABC rejection algorithm
void ABC::run()
{
	vector <Particle> particle;
	
	auto smax = model.details.num_per_core;
	for(auto s = 0u; s < smax; s++){
		output.percentage(s,smax);
		progress(s,nsample/acc_frac);
		
		auto param_val = model.param_sample();
		auto initc_val = model.initc_sample(param_val);
	
		state.simulate(param_val,initc_val);
	
		particle.push_back(state.generate_particle(UNSET,0,true));
	}
	output.percentage(smax,smax);
	
	vector <double> obs_prob;                                // Calculates the cut-off
	for(auto &part : particle){
		obs_prob.push_back(part.like.obs);
	}
	
#ifdef USE_MPI 
	obs_prob = mpi.gather(obs_prob);
#endif
	
	sort(obs_prob.begin(),obs_prob.end());
	
	auto cutoff = obs_prob[obs_prob.size()-nsample];
	
#ifdef USE_MPI 
	mpi.bcast(cutoff);
#endif

	for(auto &part : particle){
		if(part.like.obs >= cutoff){
			output.param_sample(part);
			output.state_sample(part);
		}
	}
}
	
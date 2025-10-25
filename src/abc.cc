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
	
	percentage_start(RUN_PER);
	
	auto smax = model.details.num_per_core;
	for(auto s = 0u; s < smax; s++){
		percentage(s,smax);
		
		auto param_val = model.param_sample();
		auto initc_val = model.initc_sample(param_val);
	
		state.simulate(param_val,initc_val);
	
		particle.push_back(state.generate_particle(UNSET,0,true));
	}
	
	if(false){	
		auto tot = state.timer[SIM_CALC_POPNUM]+state.timer[SIM_PRECALC]+state.timer[SIM_POPIND]+state.timer[SIM_UPDATE]+state.timer[SIM_CHECK];
		cout << cpu_percent(state.timer[SIM_CALC_POPNUM],tot) << "sim calc" << endl;
		cout << cpu_percent(state.timer[SIM_PRECALC],tot) << "sim precalc" << endl;
		cout << cpu_percent(state.timer[SIM_POPIND],tot) << "sim popind" << endl;
		cout << cpu_percent(state.timer[SIM_UPDATE],tot) << "sim update" << endl;
		cout << cpu_percent(state.timer[SIM_CHECK],tot) << "sim check" << endl;
		
		const auto &ssp = state.species[0];
		auto up = state.timer[SIM_UPDATE];
		cout <<  cpu_percent(ssp.timer[UP_MARKOV],up) << "up markov" << endl;
		cout <<  cpu_percent(ssp.timer[SORT],up) << "sort" << endl;
		cout <<  cpu_percent(ssp.timer[ITER],up) << "iter" << endl;
		cout <<  cpu_percent(ssp.timer[CHECK],up) << "check" << endl;	
	}
	
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

	for(const auto &part : particle){
		if(part.like.obs >= cutoff){
			output.param_sample(part);
			output.state_sample(part);
		}
	}
	
	percentage_end();
}
	
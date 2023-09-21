/// Implements the basic ABC rejection algorithm

#include <iostream>
#include <fstream>
#include <algorithm> 
#include <sstream> 

using namespace std;

#include "abc.hh"
#include "utils.hh"

ABC::ABC(const Model &model, Output &output) : model(model), output(output), state(model)
{	
	const auto &de = model.details;
	nsample = de.sample;
	acc_frac = de.accfrac;
}
 

/// Implements a version of simple ABC rejection algorithm
void ABC::run()
{
	auto initc_val = model.initc_sample();
	
	vector <Particle> particle;
	
	for(auto s = 0u; s < nsample/acc_frac; s++){
		if(com_op) progress(s,nsample/acc_frac);
		else if(s%100 == 0) cout << s << " s\n";
		
		auto param_val = model.param_sample();
		
		state.simulate(param_val,initc_val);
	
		particle.push_back(state.generate_particle());
	}
	
	vector <double> obs_prob;                                           // Calculates the cut-off
	for(auto &part : particle){
		obs_prob.push_back(part.like.obs);
	}
	
	sort(obs_prob.begin(),obs_prob.end());
	
	auto cutoff = obs_prob[obs_prob.size()-nsample];
	
	auto num = 0u;
	for(auto &part : particle){
		if(part.like.obs >= cutoff){
			output.param_sample(num,part);
			output.state_sample(num,part);
			num++;
		}
	}
}
	
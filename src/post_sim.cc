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

/// Initilaises the simulation
PostSim::PostSim(const Model &model, Output &output) : model(model), output(output), state(model)
{	
}


/// Performs a simulation
void PostSim::run()
{
	const auto &details = model.details;

	auto smax = details.number;

	if(model.sample.size()==0) emsg("Sample number should not be zero");
	
	for(auto s = 0u; s < smax; s++){
		if(smax > 1) cout << "Simulation " << s+1 << endl;
		if(com_op) progress(s,smax);
	
		const auto &samp = model.sample[(unsigned int)(model.sample.size()*ran())];

		auto param_val = model.post_param(samp);
		
		if(false) output.print_param(param_val);

		state.post_sim(param_val,samp);
		
		state.check("Check state");
	
		output.param_sample(s,state);
		
		output.state_sample(s,state);
	}
}

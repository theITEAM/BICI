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

	for(auto s = 0u; s < smax; s++){
		percentage(s,smax);
		
		const auto &samp = model.sample[(unsigned int)(model.sample.size()*ran())];

		auto param_val = model.post_param(samp);
		
		if(false) model.print_param(param_val);

		state.post_sim(param_val,samp);
		
		state.check("Check state");
	
		output.param_sample(s,0,state);
		output.state_sample(s,0,state);
	}
	
#ifdef USE_MPI
	mpi.barrier();
#endif
	
	percentage_end();
}

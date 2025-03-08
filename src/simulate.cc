// Implements a modified Gillespie algorithm to simulate from the model
// A timestep is used over which populations are assumed constant

#include <iostream>
#include <fstream>
#include <algorithm>
#include <assert.h>
#include <math.h>

using namespace std;

#include "simulate.hh"
#include "state.hh"
#include "utils.hh"

/// Initilaises the simulation
Simulate::Simulate(const Model &model, Output &output, Mpi &mpi) : model(model), output(output), mpi(mpi), state(model)
{	
	state.init();
}


/// Performs simulation
void Simulate::run()
{
	auto smax = model.details.num_per_core;
	for(auto s = 0u; s < smax; s++){
		output.percentage(s,smax);
		progress(s,smax);
		
		auto param_val = model.param_sample();
		auto initc_val = model.initc_sample(param_val);

		state.simulate(param_val,initc_val);

		state.check("Check state");
	
		output.param_sample(UNSET,0,state);
		output.state_sample(UNSET,0,state);
		
		if(false) state.output_waifw("waifw.txt");
	}
	
#ifdef USE_MPI
	mpi.barrier();
#endif
	
	output.percentage(smax,smax);
}

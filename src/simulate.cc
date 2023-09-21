// Implements a modified Gillespie algorithm to simulate from the model

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
Simulate::Simulate(const Model &model, Output &output) : model(model), output(output), state(model)
{	
}


/// Performs a simulation
void Simulate::run()
{
	auto smax = model.details.number;
	for(auto s = 0u; s < smax; s++){
		if(com_op) progress(s,smax);
		
		auto initc_val = model.initc_sample();
		auto param_val = model.param_sample();
		
		output.print_param(param_val);
		
		output.print_initc(initc_val);
		
		state.simulate(param_val,initc_val);
	
		output.param_sample(s,state);
		
		output.state_sample(s,state);
	}
}

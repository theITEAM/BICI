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
Simulate::Simulate(const Model &model, Output &output, Mpi &mpi, bool sup_) : model(model), output(output), mpi(mpi), state(model)
{	
	sup = sup_;
	state.init();
}


/// Performs simulation
void Simulate::run()
{
	percentage_start(RUN_PER,sup);
	
	auto smax = model.details.num_per_core;
	for(auto s = 0u; s < smax; s++){
		percentage(s,smax,sup);
	
		auto param_val = model.param_sample();
		auto initc_val = model.initc_sample(param_val);

		state.simulate(param_val,initc_val,IEstore(),s,smax,sup);

		if(slow_check) state.check("Check state");

		output.param_sample(UNSET,0,state);
		output.state_sample(UNSET,0,state);
		
		/*
		for(auto i = 0u; i< 10; i++){
			output.param_sample(UNSET,0,state);
			output.state_sample(UNSET,0,state);
		}
		*/
		
#ifdef USE_MPI
		mpi.sample_barrier(s,smax);
#endif

		if(false) state.output_waifw("waifw.txt");
	}

	if(true){
		auto tot = state.timer[SIM_ITERATE];
		const auto &ssp = state.species[0];
		
		cout << cpu_percent(state.timer[SIM_PRECALC],tot) << "sim precalc" << endl;
		cout << cpu_percent(state.timer[SIM_POPIND],tot) << "sim popind" << endl;
		cout << cpu_percent(state.timer[SIM_MARKOV],tot) << "sim markov" << endl;
		cout << cpu_percent(state.timer[SIM_UPDATE],tot) << "sim update" << endl;
		cout << cpu_percent(state.timer[SIM_NEXTPOP],tot) << "sim nextpop" << endl;
		cout << "  -> " << cpu_percent(ssp.timer[UP_MARKOV],tot) << "up markov" << endl;
		cout << cpu_percent(state.timer[SIM_LIKE],tot) << "sim like" << endl;
		cout << cpu_percent(state.timer[SIM_CHECK],tot) << "sim check" << endl;
		cout << cpu_percent(state.timer[SIM_TEMP1],tot) << "sim check end" << endl;
		
		if(model.species[0].type == INDIVIDUAL){
			auto up = state.timer[SIM_UPDATE];
			cout << "Dividing up sim update:" << endl;
			cout <<  cpu_percent(ssp.timer[SORT],up) << "sort" << endl;
			cout <<  cpu_percent(ssp.timer[ITER],up) << "iter" << endl;
			cout <<  cpu_percent(ssp.timer[CHECK],up) << "check" << endl;	
		}
	}
	
#ifdef USE_MPI
	mpi.barrier();
#endif

	percentage_end(sup);
}

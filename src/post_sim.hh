#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "model.hh"
#include "state.hh"
#include "output.hh"

class PostSim                              // Simulates from the model
{
	public:
		PostSim(const Model &model, Output &output, Mpi &mpi);
		void run();
		
	private:
		PV calculate_param_av() const;
	
		const Model &model;
		Output &output;
		Mpi &mpi;
		
		State state;
};


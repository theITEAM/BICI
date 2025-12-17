#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "chain.hh"
#include "model.hh"
#include "state.hh"

#include "output.hh"

class Extend                           // Extends inference
{
	public:
		unsigned int num_per_core;             // Number of chains per core
		
		unsigned int nsample;                  // The total number of samples
		unsigned int nburnin;                  // The total burn-in time
	
		unsigned int nparam_op;                // Store the total nunber of parameters estimated
		unsigned int nstate_op;                // Store the total nunber of states estimated
		
		unsigned int output_param;             // The number of output parameter samples
		unsigned int output_state;             // The number of output state samples
		
		Extend(const Model &model, Output &output, Mpi &mpi);
		void run();
		
	private:
		void sample_op(unsigned int s);
	
		const Model &model;
		Output &output;
		Mpi &mpi;
		
		vector <Chain> chain;
};


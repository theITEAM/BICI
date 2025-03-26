#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "model.hh"
#include "proposal.hh"
#include "state.hh"
#include "output.hh"
#include "chain.hh"

class PAS                                 // Implements PAS algorithm
{
	public:
		unsigned int num_per_core;             // Number of chains per core
		
		unsigned int nsample;                  // The total number of samples
		unsigned int nburnin;                  // The total burn-in time
		unsigned int gen_update;               // Number of updates per generation
		
		unsigned int nparam_op;                // Store the total nunber of parameters estimated
		unsigned int nstate_op;                // Store the total nunber of states estimated
		
		unsigned int output_param;             // The number of output parameter samples
		unsigned int output_state;             // The number of output state samples
		
		PAS(const Model &model, Output &output, Mpi &mpi);
		
		void run();
		
	private:
		void sample_op(unsigned int s);
		void bootstrap();
		
		double phi;
		double phi_final;
		
		unsigned int step_gen_param;             // The step used for outputing generation param
		
		double quench_factor;
		
		vector <Chain> chain;
		
		const Model &model;
		Output &output;
		Mpi &mpi;
};


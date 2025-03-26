#pragma once

#include "struct.hh"
#include "model.hh"
#include "proposal.hh"
#include "state.hh"
#include "output.hh"

class ABC
{
	public:
		ABC(const Model &model, Output &output, Mpi &mpi);
		void run();
	
	private:
		unsigned int nsample;                  // The total number of samples
		
		double acc_frac;                       // The acceptance fraction
		
		const Model &model;
		Output &output;
		Mpi &mpi;
		
		State state;
};


#ifndef BICI__ABC_HH
#define BICI__ABC_HH

#include "struct.hh"
#include "model.hh"
#include "proposal.hh"
#include "state.hh"
#include "output.hh"

class ABC
{
	public:
		ABC(const Model &model, Output &output);
		void run();
	
	private:
		unsigned int nsample;               // The total number of samples
		
		double acc_frac;                    // The acceptance fraction
		
		const Model &model;
		Output &output;
		
		State state;
};
#endif

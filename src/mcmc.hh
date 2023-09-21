#ifndef BICI__MCMC_HH
#define BICI__MCMC_HH

#include <vector>

using namespace std;

#include "struct.hh"
#include "model.hh"
#include "proposal.hh"
#include "state.hh"
#include "output.hh"

class MCMC                              // Implements an MCMC chain
{
	public:
		unsigned int nsample;               // The total number of samples
		unsigned int nburnin;
		
		unsigned int thin_param;            // The thinning factor for parameters
		unsigned int thin_state;            // The thinning factor for parameters
		
		unsigned int burnin;
		
		MCMC(const Model &model, Output &output);
		void run();
		void update();
	
	private:
		void param_prop_init();
		void diagnostics(long total_time) const;
		
		double burnfac;                 // Allows for proposals to change quickly in initial burnin 
		
		vector <Proposal> proposal;
		
		const Model &model;
		Output &output;
		
		State state;
};
#endif
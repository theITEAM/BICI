#ifndef BICI__MCMC_HH
#define BICI__MCMC_HH

#include <vector>

using namespace std;

#include "struct.hh"
#include "model.hh"
#include "proposal.hh"
#include "state.hh"
#include "output.hh"
#include "update.hh"

class MCMC                                 // Implements an MCMC chain
{
	public:
		unsigned int nsample;                  // The total number of samples
		unsigned int nburnin;                  // The total burn-in time
		
		unsigned int nparam_op;                // Store the total nunber of parameters estimated
		unsigned int nstate_op;                // Store the total nunber of states estimated
		
		unsigned int output_param;             // The number of output parameter samples
		unsigned int output_state;             // The number of output state samples
		
		BurnInfo burn_info;                    // Information about burnin
		
		//CorMatrix cor_matrix;                  // Stores information about the correlation matrix

		bool all_events_correct;               // Flag which determines if all observed events are correct
		
		MCMC(const Model &model, Output &output);
		void run();
		void add_like_obs_affect(unsigned int p, Proposal &pp) const;
		void update_state(unsigned int s);
		bool param_op(unsigned int s);
		bool state_op(unsigned int s);
	
	private:
		
		//vector <Proposal> proposal;            // Stores all the proposals
		
		const Model &model;
		Output &output;
		
		State state;
		
		Update update;
};
#endif

#ifndef BICI__PAS_HH
#define BICI__PAS_HH

#include <vector>

using namespace std;

#include "struct.hh"
#include "model.hh"
#include "proposal.hh"
#include "state.hh"
#include "output.hh"

class PAS                                 // Implements PAS algorithm
{
	public:
		unsigned int nsample;                  // The total number of samples
		unsigned int nburnin;                  // The total burn-in time
		
		unsigned int nparam_op;                // Store the total nunber of parameters estimated
		unsigned int nstate_op;                // Store the total nunber of states estimated
		
		unsigned int output_param;             // The number of output parameter samples
		unsigned int output_state;             // The number of output state samples
		
		BurnInfo burn_info;                    // Information about burnin
		
		CorMatrix cor_matrix;                  // Stores information about the correlation matrix

		bool all_events_correct;               // Flag which determines
		
		PAS(const Model &model, Mpi &mpi, Output &output);
		void run();
		void add_like_obs_affect(unsigned int p, Proposal &pp) const;
		void update(unsigned int s);
		void set_proposal_prob();
		bool param_op(unsigned int s);
		bool state_op(unsigned int s);
		void check_join_proposal();
		void join_proposal(unsigned int th1, unsigned int th2);
	
	private:
		void param_prop_init();
		void add_parameter_prop(const vector <unsigned int> &vec);
		void diagnostics(long total_time) const;
			
		//vector <Proposal> proposal;            // Stores all the proposals
		
		const Model &model;
		Mpi &mpi;
		Output &output;
		
		State state;
};
#endif

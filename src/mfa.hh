#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "model.hh"
#include "proposal.hh"
#include "state.hh"
#include "output.hh"

class MFA                                  // Implements a MFA chain
{
	public:
		unsigned int nsample;                  // The total number of samples
		unsigned int nburnin;                  // The total burn-in time
		
		//unsigned int nparam_op;                // Store the total nunber of parameters estimated
		//unsigned int nstate_op;                // Store the total nunber of states estimated
		
		//unsigned int output_param;             // The number of output parameter samples
		//unsigned int output_state;             // The number of output state samples
		
		MFA(const Model &model, Output &output);
		
		void run();
		
		/*
		BurnInfo burn_info;                    // Information about burnin
		
		
		void run();
		void add_like_obs_affect(unsigned int p, Proposal &pp) const;
		void update(unsigned int s);
		void set_proposal_prob();
		bool param_op(unsigned int s);
		bool state_op(unsigned int s);
	
	private:
		unsigned int get_op_step(unsigned int nsample) const;
		void param_prop_init();
		void diagnostics(long total_time) const;
			
		vector <Proposal> proposal;            // Stores all the proposals
		*/
		
		const Model &model;
		Output &output;
		
		State state;
};


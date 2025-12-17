#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "model.hh"
#include "proposal.hh"
#include "state.hh"
#include "output.hh"
//#include "update.hh"

class Chain {
public:
	Chain(unsigned int nburnin_, unsigned int nsample_, const Model &model, Output &output);
	void init(unsigned int ch, unsigned int ch_max);
	void burn_update(unsigned int s);
	void update(unsigned int s);
	string diagnostics(double time, double anneal_time=UNSET) const;
	void pas_burn_update(unsigned int s, unsigned int gen_update, double phi);
	void join_proposal_update();
	unsigned int get_nburnin();
	void pas_burn_update_run(unsigned int s);
	double like_total_obs() const;
	void check_cor_matrix() const;
	vector <PropInfo> get_prop_info() const;
	void set_prop_info(const vector <PropInfo> &prop_info);
	TerminalInfo get_terminal_info(unsigned int ch_tot) const;

	State state;
	double Lobs_av;                        // Calculates variance in Lobs (for PAS)
	double Lobs_av2;                        
	double nLobs_av;                       

	CorMatrix cor_matrix;                  // Stores information about the correlation matrix

	vector <Proposal> proposal;            // Stores all the proposals 

private:
	unsigned int nsample;                  // The total number of samples
	unsigned int nburnin;                  // The total burn-in time
	
	BurnInfo burn_info;                    // Information about burnin
		
	void update_init();
	bool is_markov_pop(const AffectLike &al) const;
	void add_parameter_prop(const vector <unsigned int> &vec);
	void set_proposal_prob();
	void check_join_proposal();
	unsigned int find_which_list(unsigned int j, const vector < vector <unsigned int> > &par_list) const;
	string banner(string te) const;
	
	const Model &model;
	Output &output;
};


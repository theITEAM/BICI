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
	void pas_burn_update(unsigned int s, unsigned int g, unsigned int gen_update, double phi);
	void pas_burn_update_run(unsigned int s);
	double like_total_obs() const;
	void check_cor_matrix() const;
	
	State state;
	double Lobs_av;                        // Calculates Lobs (for PAS)
	double nLobs_av;                        // Calculates Lobs (for PAS)
	
private:
	unsigned int nsample;                  // The total number of samples
	unsigned int nburnin;                  // The total burn-in time
	
	BurnInfo burn_info;                    // Information about burnin
		
	CorMatrix cor_matrix;                  // Stores information about the correlation matrix

	vector <Proposal> proposal;            // Stores all the proposals

	void update_init();
	bool is_markov_pop(const AffectLike &al) const;
	void add_parameter_prop(const vector <unsigned int> &vec);
	void set_proposal_prob();
	void check_join_proposal();
	unsigned int find_which_list(unsigned int j, const vector < vector <unsigned int> > &par_list) const;
	string banner(string te) const;
	//void join_proposal(unsigned int th1, unsigned int th2);

	const Model &model;
	Output &output;
};


#ifndef BICI__CHAIN_HH
#define BICI__CHAIN_HH

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
	void init();
	void burn_update(unsigned int s);
	void update(unsigned int s);
	string diagnostics(unsigned int time, unsigned int anneal_time=UNSET) const;
	void pas_burn_update(unsigned int s, unsigned int g, unsigned int gen_update, double phi);
	double like_total_obs() const;
	
	State state;
	
private:
	unsigned int nsample;                  // The total number of samples
	unsigned int nburnin;                  // The total burn-in time
		
	BurnInfo burn_info;                    // Information about burnin
		
	CorMatrix cor_matrix;                  // Stores information about the correlation matrix

	vector <Proposal> proposal;            // Stores all the proposals

	//Update update;	
	void update_init();
	void add_parameter_prop(const vector <unsigned int> &vec);
	void set_proposal_prob();
	void check_join_proposal();
	void join_proposal(unsigned int th1, unsigned int th2) ;
	
	const Model &model;
	Output &output;
};
#endif

#ifndef BICI__UPDATE_HH
#define BICI__UPDATE_HH

#include <vector>

using namespace std;

#include "struct.hh"
#include "model.hh"
#include "proposal.hh"
#include "state.hh"
#include "output.hh"

class Update                                 // Implements update to a state
{
public:
	CorMatrix cor_matrix;                      // Stores information about the correlation matrix

	vector <Proposal> proposal;                // Stores all the proposals

	Update(const Model &model, const BurnInfo &burn_info, Output &output);
	void run(unsigned int s, State &state);
	void init();
	void diagnostics(long total_time, const State &state) const;
	void set_proposal_prob();
	void check_join_proposal();
		
private:		
	void add_parameter_prop(const vector <unsigned int> &vec);
	void join_proposal(unsigned int th1, unsigned int th2);
	
	const Model &model;
	const BurnInfo &burn_info;                    // Information about burnin
	Output &output;	

};
#endif

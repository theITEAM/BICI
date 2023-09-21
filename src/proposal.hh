#ifndef BICI__PARAMPROP_HH
#define BICI__PARAMPROP_HH

#include <vector>

using namespace std;

#include "struct.hh"
#include "model.hh"
#include "state.hh"
#include "output.hh"

class Proposal                              // Implements a proposal
{
	public:
		Proposal(PropType type_, vector <unsigned int> vec, const Model &model, const Output &output, double prob, const double &burnfac);
		void initialise_variables();
		void get_dependency();
		void get_affect_like();
		void set_mvn(vector <Particle> &particle, double si_);
		void setup_sum();
		void add_sample(const vector <double> &param_val);
		double mvn_probability(const vector <double> &param1, const vector <double> &param2) const;
		void update_sampler();
		vector <double> sample(vector <double> param_val) const;
		void print_info() const;
		void MH(State &state, bool burnin);
		
		PropType type;                   // The type of proposal
		
		string name;                      // Name of proposal
		
		// PARAM_PROP
		vector <unsigned int> param_list; // Lists parameter used in proposals
		
		vector <unsigned int> dependent;  // Lists parameters which are dependent (based on reparam)
		
		unsigned int N;                   // The number of parameters
		unsigned int ntr;                 // Number of times proposal has been performed
		unsigned int nac;                 // Number of times proposal has been accepted
		double si;                        // Jumping size of proposal
		vector <double> mean_sum;         // Sums up parameter samples
		vector < vector <double> > M_sum; // Sums up parameter samples
		unsigned int nsum;
		vector < vector <double> > M;     // The covariance matrix
		
		vector < vector <double> > Z;     // Cholesky matrix
		
		vector < vector <double> > M_inv;  // The covariance matrix
		
		double prop_prob;                 // The probability of performing proposal update
		
		long timer;                       // Times how long the proposal takes
		
		vector <AffectLike> affect_like;
		
		// IND_PROP
		unsigned int p_prop;                // The species for the proposals
		
		vector <IndSampler> ind_sampler;    // Information about the individual samplers
		
		void initialise_ind_variable();
		void MH_ind(State &state, bool burnin);
		void sample_ind(State &state);
		
		// IE_PROP
		unsigned int ie_prop;               // The individual effect being changed
		
		void MH_ie(State &state, bool burnin);
		
		void MH_ie_var(State &state, bool burnin);
		
		// MBP_PROP
		void mbp(MBPType type, State &state, bool burnin);
		
		double calculate_al(const Like &like_ch, double prop_prob) const;
		void mbp_population_affect();
		
		string diagnostics(long total_time) const;
		void update(State &state, bool burnin);
		
	private:
		const Model &model;
		const Output &output;
		const double &burnfac;
};
#endif

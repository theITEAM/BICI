#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "model.hh"
#include "state.hh"
#include "output.hh"

class Proposal                             // Implements a proposal
{
	public:
		Proposal(PropType type_, vector <unsigned int> vec, const Model &model, const Output &output, double w, const BurnInfo &burn_info);
		void initialise_variables();
		
		string print_info() const;
		void MH(State &state);
		
		unsigned int plot_index;
		
		PropType type;                         // The type of proposal
		
		bool on;                               // Determines if proposals is switched on
		
		string name;                           // Name of proposal
		
		// PARAM_PROP
		vector <unsigned int> param_list;      // Lists parameter used in proposals
		
		vector <unsigned int> dependent;       // Lists parameters which are dependent (based on reparam)
		
		unsigned int N;                        // The number of parameters
		unsigned int ntr;                      // Number of times proposal has been performed
		unsigned int nac;                      // Number of times proposal has been accepted
		unsigned int nfa;                      // Number of times proposal fails
		double si;                             // Jumping size of proposal
		vector < vector <double> > M;          // The covariance matrix
		
		vector < vector <double> > Z;          // Cholesky matrix
		
		vector < vector <double> > M_inv;      // The covariance matrix
		
		bool omega_check;                      // Set if omega needs to be checked
		double prop_weight;                    // The weight given to certain proposals
		
		double prop_prob;                      // The probability of performing proposal update
		
		unsigned int number;                   // The number of times the proposal is done
		
		vector <long> timer;                   // Times how long the proposal takes
		
		LocalSampler loc_samp;                 // Sampler used for local event changes
		
		vector <GeneticProp> gen_mut_info;     // Stores information about genetic mutation proposals
		
		vector <SwapProp> swap_info;           // Stores information about local swaps [#swap_rep]
		
		vector <AffectLike> affect_like;       // Stores how the likelihood is affected by proposal
		vector <AffectLike> affect_spline;     // Stores how splines are affected by proposal
			
		// IND_PROP
		unsigned int p_prop;                   // The species for the proposals
		
		unsigned int cl_prop;                  // The classifciation for the proposal (MBP_IC_POP_PROP)
		
		unsigned int c_prop;                   // The compartment for the proposal (MBP_IC_POP_PROP)
		
		vector <IndSampler> ind_sampler;       // Information about the individual samplers
		
		vector <bool> tr_change;               // Store transition change under par/event joint prop
			
		bool all_events_correct;               // Flag determining if all obs events correct
	
		void MH_event(State &state);
		void MH_multi_event(State &state);
		void MH_event_all(State &state);
		void MH_ind_local(State &state);
		void sample_ind_obs(State &state);
		void trans_tree(State &state);
		void trans_tree_swap_inf(State &state);
		void trans_tree_mut(State &state);
		void trans_tree_mut_local(State &state);
		void resimulate_ind_obs(State &state);
		void resimulate_single_ind_obs(State &state);
		void resimulate_ind_unobs(State &state);
		
		// IE_PROP
		unsigned int ie_prop;                 // The individual effect being changed
		IndEffGroupRef ind_eff_group_ref;     // Reference an element in individual effect (for covar proposals)
		
		void add_rem_ind(State &state); 
		void add_rem_tt_ind(State &state);	
		void MH_ie(State &state);
		void MH_ie_var(State &state);
		void MH_ie_covar(State &state);
		void MH_ie_var_cv(State &state);
		
		// MBP_PROP
		void mbp(State &state);
		void param_event_joint(Direction dir, State &state);
		void init_cond_frac(State &state);
		double calculate_al(const Like &like_ch, double prop_prob) const;
		string diagnostics(long total_time) const;
		void update(State &state);
		void switch_enter_source(State &state);
		void switch_leave_sink(State &state);
		void correct_obs_trans_events(State &state);
		void enforce_add_source(unsigned int i, double t, State &state);
		void enforce_add_sink(unsigned int i, double t, State &state);
		void enforce_add_trans(unsigned int i, const ObsData &ob, State &state);
		bool try_move_obs_trans_event(unsigned int i, unsigned int e, double t, const ObsData &ob, const vector <Event> &ev, State &state);
		void basic_ind_update(unsigned int i, vector <Event> &ev_new, State &state);
	
	private:
		const Model &model;
		const Output &output;
		const BurnInfo &burn_info;
	
	// Used in proposal_utils.cc
	public:
		void update_sampler(const CorMatrix &cor_matrix);
		void calculate_affect_spline();
		void mbp_population_affect();
		void set_mvn(double si_, const CorMatrix &cor_matrix);
		vector <double> sample(vector <double> param_val);
		double mvn_probability(const vector <double> &param1, const vector <double> &param2) const;
		
	private:
		void get_dependency();
		void get_affect_like();
		ICResult propose_init_cond(InitCondValue &icv, const State &state);
		vector <unsigned int> multinomial_resample(const vector <unsigned int> &x, const vector <double> &frac, double prob) const;
		void multinomial_reduce(unsigned int dN, vector <unsigned int> &x, const vector <double> &frac) const;
		void add_li_cha(unsigned int i, const vector <LocalIndChange> &add, vector <LocalIndChange> &licha, vector < vector <unsigned int> > &lc_ref) const;
		void remove_li_cha(unsigned int i, vector <LocalIndChange> &licha, vector < vector <unsigned int> > &lc_ref) const;
		void initialise_ind_variable();
		void initialise_ind_time_sampler();
		void initialise_swap_variable();
		void update_si(double fac);
		void update_ind_samp_si(unsigned int tr_gl, double fac);
		bool skip_proposal(double val) const;
		void set_omega_check();
	
	// Used in proposal_local.cc
	private:
		unsigned get_win() const;
		void update_win(unsigned int num, unsigned int num_opt);
		void pop_add_rem_local(State &state);
		void pop_move_local(State &state);
		void pop_ic_local(State &state);
		void pop_end_local(State &state);
		void pop_single_local(State &state);	
		void pop_ic(State &state);
		void pop_ic_swap(State &state);
};


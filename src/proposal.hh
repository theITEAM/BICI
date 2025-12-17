#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "model.hh"
#include "state.hh"
#include "output.hh"
#include "cor_matrix.hh"

class Proposal                             // Implements a proposal
{
	public:
		Proposal(PropType type_, vector <unsigned int> vec, const Model &model, const Output &output, double w, const BurnInfo &burn_info);
		void initialise_variables();
		void calculate_spec_precalc();
		vector <unsigned int> remove_precalc_done(const vector <unsigned int> &list_precalc, vector <bool> &mapl) const;
		
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
		
		vector <IEGref> ieg_check_pri;         // Stores which ind effect group priors need to be checked
		Hash hash_ieg_check_pri;
		
		double prop_weight;                    // The weight given to certain proposals
		
		double prop_prob;                      // The probability of performing proposal update
		
		unsigned int number;                   // The number of times the proposal is done
		
		vector <double> timer;                 // Times how long the proposal takes
		
		LocalSampler loc_samp;                 // Sampler used for local event changes
		
		vector <GeneticProp> gen_mut_info;     // Stores information about genetic mutation proposals
		
		vector <SwapProp> swap_info;           // Stores information about local swaps [#swap_rep]
		
		vector <AffectLike> affect_like;       // Stores how the likelihood is affected by proposal
		
		vector <SpecPrecalc> dependent_spec_precalc; // Works out update precalc for parameter sampling
		
		SpecPrecalc spec_precalc_after;        // Works out how to update precalculation
		
		// IND_PROP
		unsigned int p_prop;                   // The species for the proposals
		
		unsigned int cl_prop;                  // The classifciation for the proposal (MBP_IC_POP_PROP)
		
		unsigned int c_prop;                   // The compartment for the proposal (MBP_IC_POP_PROP)
		
		vector <IndSampler> ind_sampler;       // Information about the individual samplers
		
		vector <IndSimProb> ind_sim_prob;      // Probability of proposing ind sim prop [i]
		
		vector < vector <IndSimProb> > ind_sim_prob_cl;// Probability of proposing ind sim prop [i][cl]
		
		vector <bool> tr_change;               // Store transition change under par/event joint prop
			
		bool all_events_correct;               // Flag determining if all obs events correct
	
		MBPfast mbp_fast;                      // This is used to speed up calculation for MBPs
		
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
		string diagnostics(double total_time) const;
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
		void mbp_population_affect();
		void set_mvn(double si_, const CorMatrix &cor_matrix);
		double param_resample(PV &param_val, const vector < vector <double> > &popnum_t);
		double mvn_probability(const vector <double> &param_prop1, const vector <double> &param_prop2) const;
		bool prop_info_on() const;
		PropInfo get_prop_info() const;
		void set_prop_info(const PropInfo &pi);
		vector <unsigned int> get_prop_id() const;
		
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
		void update_si(PropResult res);
		void update_ind_samp_si(unsigned int tr_gl, PropResult res, unsigned int ntr);
		//bool skip_proposal(double val) const;
	
		void set_ieg_check_pri();
		//void set_omega_check();
		bool event_dif(const vector <Event> &ev1, const vector <Event> &ev2) const;
		void ind_obs_prob_update(IndSimProb &isp) const;
		void set_mbp_fast();
		double set_prop_prob();
		void add_sampler_info(vector <double> &vec, const Sampler &sa) const;
		void set_sampler_info(const vector <double> &vec, Sampler &sa) const;
		void add_sampler2D_info(vector <double> &vec, const Sampler2D &sa) const;
		void set_sampler2D_info(const vector <double> &vec, Sampler2D &sa) const;
		string print_ac(unsigned int nac, unsigned int ntr) const;
		string print_fa(unsigned int nfa, unsigned int ntr) const;
		string print_range(string te, const vector <double> &list) const;

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


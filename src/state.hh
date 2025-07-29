#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "model.hh"

class State                                // Stores information about the state
{
	public:
		unsigned int sample;                   // Sample number (for diagnostics)
		
		unsigned int T;                        // Total number of divisions
		
		vector <double> dtimepoint;            // The difference between timepoints
		
		vector <unsigned int> empty;           // An empty vector
		
		vector <double> param_val;             // The parameter value
		
		vector <SplineValue> spline_val;       // Values for the spline
		
		Like like;                             // Stores the overall likelihood and prior
		
		vector <double> prior_prob;            // The prior probability for each parameter
		
		vector <double> dist_prob;             // The distribution probability for each parameter
		
		vector <double> spline_prior;          // The prior probability a spline

		vector < vector <double> > popnum_t;   // The time variation in the population [ti][#]
		
		vector < vector < vector <PopIndRef> > > popnum_ind;// Individuals associated with population (for trans tree)
		
		vector <StateSpecies> species;         // State information about species
		unsigned int nspecies;

		GeneticDataValue genetic_value;        // Stores values for genetic data
		
		vector <Change> change;                // Stores changes in the state

		vector <double> update_timer;          // Stores timers for parameter updates
		vector <double> restore_timer;         // Stores timers for parameter restore
		vector <double> timer;                 // General purpose timers
		
		double dif_thresh;                     // Threshold for difference
		
		vector <BackPop> back_pop;             // Stores population info
		
		vector <AlgWarn> alg_warn;             // Stores any algorithm warnings
		
		State(const Model &model);
		void init();
		void simulate(const vector <double> &param_value, const vector <InitCondValue> &initc_val);
		void post_sim(const vector <double> &param_value, const Sample &samp);
		void simulate_iterate(unsigned int ti_start, unsigned int ti_end);
		void ensure_all_ind_event();
		vector <DeriveOutput> derive_calculate() const;
		void calculate_likelihood();
		void calculate_like();
		void add_like(Like like_ch);
		void change_add(const vector <double> &vec);
		void change_add(double num);
		void change_add();
		void update_spline(const vector <AffectLike> &affect_spline);
		void restore_spline(const vector <AffectLike> &affect_spline);
		void remove_store_spline(const vector <AffectLike> &affect_spline);
		Like update_param(const vector <AffectLike> &affect_like, const vector <double> &param_store);
		void restore(const vector <AffectLike> &affect_like);
		void likelihood_from_scratch();
		void resample_ind(bool do_pl = false);
		Particle generate_particle(unsigned int s, unsigned int chain, bool store_state) const;
		void set_particle(const Particle &part, bool calc_like=true);
		vector <double> prior_init_cond(double &like_ch);
		void update_individual_sampler();
		void update_popnum_t_init(unsigned int p, vector <unsigned int> cnum_i, vector <unsigned int> cnum_f);
		void back_init();
		void restore_back();
		unsigned int get_ti(double t) const;
		vector <double> get_param_val_prop() const;
		double get_trans_rate_est(unsigned int p, unsigned int tr, unsigned int ti) const;
		vector < vector <double> > get_population_rates(unsigned int p) const;
		void init_cor_matrix(); 
	
	private:
		vector <double> calculate_df(const DerFunc &df) const;
		string compact_vector(const vector <double> &value) const;
		void spline_init();
		void print_cpop(unsigned int ti) const;	
		
	// In state_update_ind.cc
	public:
		void update_markov_eqn_value(MarkovEqnDiv &div, double value, double &Li_markov, double &like_ch);
		Like update_ind(unsigned int p, unsigned int i, vector <Event> &ev_new, UpdateType type);
		void update_pop_change(unsigned int ti, unsigned int ti_next, const vector <PopChange> &pop_change, double &like_ch);
		void update_ie_population(unsigned int p, unsigned int i, unsigned int ie, double factor, double &like_ch);
		void update_ie_trans_tree(unsigned int p, unsigned int i, unsigned int ie, double ratio, double &like_ch);
		double get_w_from_indinffrom(const IndInfFrom &iif) const;
		
	private:
		void initialise_update_ind_maps();
		vector <unsigned int> pop_list;        // Used in state_update_ind
		vector <int> pop_map;
		vector < vector <bool> > markov_eqn_map; 
		vector <PopMarkovEqnRef> markov_eqn_list;
		vector < vector <bool> > trans_map; 
		vector <PopTransRef> trans_list;
		
	// In state_genetic.cc
	public:	
		GenChange update_tree(unsigned int p, unsigned int i, vector <Event> &ev_new);
		double sample_infection_source(Event &ev, unsigned int p) const;
		bool trg_contains_outside(unsigned int p, unsigned int tr_gl) const;
		double prob_infection_source(const Event &ev, unsigned int p) const;
		NodeRef pass_down(const NodeRef &nr, double t, double probif) const;
		NodeRef pass_up(const NodeRef &nr, double t, double probfi) const;
		InfNodeAlter get_inf_node_alteration(unsigned int p, unsigned int i,unsigned int n, const vector <Event> &ev_new, unsigned int e) const;
		void setup_transtree();
		void sample_genetic_value();
		double find_t_rec(const InfNode &in) const;
		vector <double> set_genetic_param();
		vector < vector <unsigned int> > calculate_gen_dif();
		void calculate_dif(unsigned int n, unsigned int e_start, unsigned int num,  vector <unsigned int> &dif) const;
		void calculate_dif_rev(unsigned int n, int e_start, unsigned int num, vector <unsigned int> &dif, unsigned int &root) const;
		void find_all_obs(unsigned int node, vector <unsigned int> &obs_list) const;
		vector < vector <Poss> > calculate_pop_ind() const;
		vector < vector < vector <Poss> > > calculate_pop_ind_total() const;
		double likelihood_genetic_process(); 
		double likelihood_genetic_obs(const vector < vector <unsigned int> > &gen_val);
		void trans_tree_proposal(const BurnInfo &burn_info, unsigned int &nac, unsigned int &ntr);
		bool set_ind_inf_from(double t, unsigned int p, vector <Event> &event, unsigned int p_from, unsigned int i_from, const vector <Event> &event_from);
		void trans_tree_swap_inf_proposal(const BurnInfo &burn_info, unsigned int &nfa, unsigned int &nac, unsigned int &ntr);
		void update_inf_node_ref(unsigned int k);
		InfNodeAlter no_alter(unsigned int n) const;
		void trans_tree_mut_proposal(const BurnInfo &burn_info, unsigned int &nac, unsigned int &ntr, double &si);
		void trans_tree_mut_local_proposal(const BurnInfo &burn_info, vector <GeneticProp> &gen_mut_info);
		void print_node(string te, unsigned int n) const;
		void print_inf_origin() const;
		double round_down(double t) const;
		double round_up(double t) const;
		void disconnect_inf_node(unsigned int n);
		void connect_inf_node(unsigned int n, const IndInfFrom &iif, unsigned int mut_num_store);
		
		void check_pop_ind(vector < vector <unsigned int> > pop_ind, const vector < vector <IndPop> > &ind_pop) const;
		NodeRef get_ind_noderef(unsigned int p, unsigned int i, double t) const;
		GenChange gen_change(GenChaType type, unsigned int n, const InfNodeAlter &node_alter, const IndInfFrom &iif_add) const; 
		void gen_change_update(const GenChange &gc);
		void remove_node(unsigned int n);
		void update_inf_ev(unsigned int n);
		void calculate_popnum_ind();
		void add_popnum_ind(unsigned int p, unsigned int i);
		void update_popnum_ind(unsigned int p, unsigned int i);
		void popnum_ind_recalc_w(unsigned int p, unsigned int i);
		void change_add_node(double t, const NodeRef &nr_add, GenChange &gc) const;
		void change_remove_node(double t, const NodeRef &nr_from, GenChange &gc) const;
		vector <unsigned int> get_gen_obs_list(unsigned int n, unsigned int j, unsigned int n_not_allow=UNSET) const;
		void get_gen_obs_list2(unsigned int n, unsigned int j, vector <unsigned int> &list, unsigned int n_not_allow) const;
		void output_gen_dif() const;
		void output_waifw(string file) const;
		vector <double> recalculate_iif_w();
		void recalculate_iif_w_undo(const vector <double> &vec);
		vector <double> recalculate_popnum_ind_w(unsigned int po);
		void recalculate_popnum_ind_w_undo(unsigned int po, const vector <double> &vec);
		 
	// In state_local.cc
	public:
		Like calculate_local_change(unsigned int p, vector <LocalChange> &local_change, int dir, bool &ill);
	
	private:
		void initialise_local_map(); 
		
		vector < vector <int> > comp_local_map;    // A map for local changes in compartment [p][c]
		vector < vector <double> > pop_trans_map;  // A map for changes in pop_trans_data [p][ref]
		vector < vector <double> > pop_data_map;   // A map for changes in pop_data [p][ref]
	
	// In like_integrate
	public:
		double like_integrate(const vector <double> &param_value);
	
	// In state_check.cc
	public:
		void check_trans_num(string ref);
		void check_simp(string ref);
		void check(string ref);
		void scan_variable(string name, double min, double max);
		void print_ev_data(string te, const vector <EventData> &event_data, unsigned int p) const;
		void check_popnum_t(string ref);
		void check_popnum_t2(string ref);
		void check_neg_rate(string name);
		void add_alg_warn(string te);
	
	private:
		void check_dependent_param(string ref);
		void check_ref(unsigned int p, string ref);
		void check_markov_trans(unsigned int p, string ref);
		void check_nm_trans(unsigned int p, string ref);
		void check_prior(string ref);
		void check_like(string ref);
		void check_spline(string ref);
		void check_ie(unsigned int p, string ref);
		void check_pop_like(unsigned int p, string ref);
		void check_cpop_st(string ref);
		void check_maps(string ref);
		void check_event_observed(unsigned int p, string ref);
		void check_obs_like(unsigned int p, string ref);
		void check_init_cond_like(unsigned int p, string ref);
		void check_init_cond_prior(string ref);
		void check_linearise();
		void check_spline_store(string ref);
		void check_genetic_value(string ref);
		void check_popnum_ind(string ref);
		void check_add_move_rem(string ref);
		void check_event_boundary(string ref);
		void check_effect_out_of_range();
		void output_dump() const;
		void get_closest_to_boundary();
		
	private:
		const Model &model;
};


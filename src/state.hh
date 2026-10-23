#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "model.hh"
#include "state_species.hh"

class State                                // Stores information about the state
{
	public:
		unsigned int sample;                   // Sample number (for diagnostics)
		
		unsigned int T;                        // Total number of divisions
		
		vector <unsigned int> empty;           // An empty vector
		
		PV param_val;                          // The parameter value
	
		Like like;                             // Stores the overall likelihood and prior
		
		vector <double> prior_prob;            // The prior probability for each parameter
	
		vector <double> prior_ieg;             // The prior for ind effect groups
		
		vector <double> dist_prob;             // The distribution probability for each parameter
		
		vector <double> spline_prior;          // The prior probability a spline

		vector < vector <double> > popnum_t;   // The time variation in the population [ti][#]
		
		vector <double> popcombw_value;        // The weights for different population combinations
		
		vector < vector <double> > popcomb_t;  // Time variation in population combinations
		
		vector <double> popcomb_store;         // Stores results so can be undone
	
		vector < vector < vector <PopIndRef> > > popnum_ind;// Individuals associated with population (for trans tree)
		
		vector <StateSpecies> species;         // State information about species
		unsigned int nspecies;

		GeneticDataValue genetic_value;        // Stores values for genetic data
		
		vector <Change> change;                // Stores changes in the state

		vector <double> update_timer;          // Stores timers for parameter updates
		vector <double> restore_timer;         // Stores timers for parameter restore
		vector <double> timer;                 // General purpose timers
		vector <double> check_timer;           // General purpose timers
		
		double dif_thresh;                     // Threshold for difference
		
		vector <BackPop> back_pop;             // Stores population info
		
		vector <AlgWarn> alg_warn;             // Stores any algorithm warnings
		
		bool trans_tree_on;                    // Determines if transmission tree is turned on
	
		// Used in simualtion
		vector <double> dpop;                  // Stores change in population
		vector <unsigned int> dpop_list;       // Lists changes
		
		vector <bool> map_pop;                 // Maps used in function popcomb_calc
		vector <bool> map_popcomb;
		
		RegraftObsType regraft_obs_type;       // This determines how regrafting is done using genetic observations (with transmission tree)
		
		//vector <double> dpopcomb;              // Stores change in popcomb
		//vector <unsigned int> dpopcomb_list;   // Lists changes in popcomb
		
		State(const Model &model);
		void init();
		void simulate(const PV &param_value, const vector <InitCondValue> &initc_val, IEstore ie_store= IEstore(), double val=UNSET, double val2=UNSET, bool sup=false);
		void determinisitic_resimulate();
		void individual_initial_state(double t);
		void post_sim(const PV &param_value, const Sample &samp);
		void load_samp(const PV &param_value, const Sample &samp);
		void simulate_iterate(unsigned int ti_start, unsigned int ti_end, double val=UNSET, double val2=UNSET, bool sup=false);
		void dpop_clear();
		void update_quantities(unsigned int ti, unsigned int ti_start, bool calc_markov, bool calc_precalc, bool include_derive);
		void ensure_all_ind_event();
		vector <DeriveOutput> derive_calculate(bool store_state);
		void calculate_likelihood();
		void calculate_like();
		void accept_update_ind(const UpdateIndInfo &uii);
		void accept(const Like &like_ch);
		void reject();
		void change_add(const vector <double> &vec);
		void change_add(double num);
		void change_add();
		Like update_param(const vector <AffectLike> &affect_like);
		void restore(const vector <AffectLike> &affect_like);
		void likelihood_from_scratch();
		void set_inconsistent_ind();
		Particle generate_particle(unsigned int s, unsigned int chain, bool store_state, bool dir_fl = true);
		void remove_trans_tree();
		void set_particle(const Particle &part);
		void regenerate(unsigned int ti_end, bool calc_markov, bool calc_precalc, bool include_derive=false);
		vector <double> prior_init_cond(double &like_ch);
		void update_individual_sampler();
		void update_pop_t_init(unsigned int p, vector <unsigned int> cnum_i, vector <unsigned int> cnum_f);
		void back_init();
		void restore_back();
		vector <double> get_param_val_prop() const;
		vector <double> get_trans_rate_est_para(const vector <unsigned int> &list, unsigned int p, unsigned int tr) const;
		void init_cor_matrix(); 
		void set_ie_from_samp(const IEstore &ie_store);
		void set_ie_from_data();
		void ie_finalise();
		double calculate(const EquationInfo &ei, unsigned int ti) const;
		void calculate_popcombw_value(bool include_derive);
		double popcombw_calc(unsigned int i);
		void mbp_update_initial_state(const vector <unsigned int> &mbp_popcombw_affect, vector <PopcombWStore> &popcombw_store);
		void pop_recalc(const PopChangeInfo &pop_change_info);
		void pop_restore();
		double frac_outside_CI();
		Like get_like_ch(const Like &like_st) const;
		bool extinction_criterion(const ExtinctInfo &ei) const;
		
	private:
		vector <double> calculate_df(const DerFunc &df) const;
		string compact_vector(const vector <double> &value) const;
		void spline_init();
		void print_cpop(unsigned int ti) const;	
		vector <double> calculate_popcomb(unsigned int ti, bool include_derive = false) const;
		vector <double> calculate_popnum(bool include_derive = false) const;
		void add_popnum(unsigned int c, vector <double> &popnum, const Individual &ind, const PopRef &pr) const;
		vector < vector <double> > calculate_popnum_t(unsigned int ti_end = UNSET);
		vector < vector <double> > calculate_popcomb_t(unsigned int ti_end = UNSET);
		//void calculate_expand_pop_derive();
		//void calculate_popcomb_derive_end();
		void recalculate_population(const vector <unsigned int> &list);
		void recalculate_popcomb(const vector <unsigned int> &popcomb_affect);
		//void recalculate_population_restore(const vector <unsigned int> &list, const vector <double> &vec);
		//void check_w() const;
		IndInfFrom get_ind_inf_from(double t, unsigned int p, unsigned int tr_gl, unsigned int p_from, unsigned int i_from, const vector <Event> &event_from) const;
		void set_ind_inf_from(const vector <string> &ind_key);
		
		
	// In state_update_ind.cc
	public:
		UpdateIndInfo update_ind(unsigned int p, unsigned int i, vector <Event> &ev_ne, UpdateType type, bool update_trans_tree = true);
		//Like update_ind(unsigned int p, unsigned int i, vector <Event> &ev_new, UpdateType type);
		void remove_trans_tree_ind_inf_likelihood(const vector <InfPoCha> &ind_po_cha, double &dmarkov);
		void add_trans_tree_ind_inf_likelihood(const vector <InfPoCha> &ind_po_cha, double &dmarkov);
		void change_pop_t(unsigned int ti, unsigned int ti_next, unsigned int k, double num);
		void change_population(unsigned int ti, unsigned int ti_next, unsigned int k, double num);
		void update_pop_change(unsigned int ti, unsigned int ti_next, double &like_ch);
		void update_ie_population(unsigned int p, unsigned int i, unsigned int ie, double factor, double &like_ch);
		void update_ie_trans_tree(unsigned int p, unsigned int i, unsigned int ie, double ratio, double &like_ch);
		double get_w_from_indinffrom(const IndInfFrom &iif) const;
		
	private:
		void initialise_update_ind_maps();
		vector <unsigned int> pop_list;        // Used in state_update_ind
		vector <int> pop_map;
		vector < vector <bool> > markov_eqn_map; 
		vector < vector <unsigned int> > markov_eqn_list;
		vector < vector <bool> > trans_map; 
		vector < vector <unsigned int> > trans_list;
		//vector <bool> spline_up_map; 
		//vector <unsigned int> spline_up_list;
		vector <bool> reparamth_up_map; 
		vector <unsigned int> reparamth_up_list;
		
	// In state_genetic.cc
	public:	
		TransTreeChange update_tree(unsigned int p, unsigned int i, vector <Event> &ev_new);
		double sample_infection_source(Event &ev, unsigned int p) const;
		bool set_iif(unsigned int p, unsigned int i, unsigned int e, IndInfFrom &iif, double t, const Event &ev, const Species &sp) const;
		bool trg_contains_outside(unsigned int p, unsigned int tr_gl) const;
		double prob_infection_source(const Event &ev, unsigned int p) const;
		NodeRef pass_down(const NodeRef &nr, double t, double &probif) const;
		NodeRef pass_up(const NodeRef &nr, double t, double &probfi) const;
		InfNodeAlter get_inf_node_alteration(unsigned int p, unsigned int i,unsigned int n, const vector <Event> &ev_new, unsigned int e) const;
		void setup_inf_node();
		void setup_transtree(bool set_iif);
		unsigned int get_nearest_gobs(double tdiv, const vector <unsigned int> &list) const;
		void setup_iif_with_genetic_data();
		void minimise_mut_num();
		void sample_genetic_value();
		double find_t_rec(const InfNode &in) const;
		vector <double> set_genetic_param();
		vector < vector <unsigned int> > calculate_gen_dif() const;
		void calculate_dif(unsigned int n, unsigned int e_start, unsigned int num,  vector <unsigned int> &dif) const;
		void calculate_dif_rev(unsigned int n, int e_start, unsigned int num, vector <unsigned int> &dif, unsigned int &root) const;
		void find_all_obs(unsigned int node, vector <unsigned int> &obs_list) const;
		vector < vector <Poss> > calculate_pop_ind() const;
		vector < vector < vector <Poss> > > calculate_pop_ind_total() const;
		double likelihood_genetic_process(); 
		double likelihood_genetic_obs(const vector < vector <unsigned int> > &gen_val) const;
		void update_like_markov(unsigned int m, unsigned int ti, double dlike_markov, unsigned int p);
		void trans_tree_proposal(const BurnInfo &burn_info, unsigned int &nac, unsigned int &ntr, bool gen_obs_max);
		NodeMove node_move(unsigned n, int e, double t, unsigned int num_max) const;
		void local_trans_tree_proposal(const BurnInfo &burn_info, unsigned int &nac, unsigned int &ntr, bool gen_obs_max);
		bool set_ind_inf_from(double t, unsigned int p, vector <Event> &event, unsigned int p_from, unsigned int i_from, const vector <Event> &event_from);
		void set_ind_inf_from_outside(double t, vector <Event> &event);
		void add_ind_inf_change(unsigned int p_ev, unsigned int i_ev, unsigned int e_ev, unsigned int c_bef, unsigned int c_aft, unsigned int p_aft, unsigned int i_aft, vector <IIFChange> &iif_change, bool &ill_fl) const;
		double get_branch_ev(vector <BranchEv> &branch, unsigned int n, int e, unsigned num, double t) const;
		double branch_sample_prob(vector <BranchEv> &branch, bool resample, string &warn) const;
		void print_new_inf_ev(string te, const vector <NewInfEv> &list) const;
		void check_branch(const vector <BranchEv> &branch) const;
		void accept_new_inf_ev(const vector <NewInfEv> &list);
		vector <NewInfEv> construct_new_inf_ev(const vector <BranchEv> &branch, vector <IIFChange> &iif_change, bool &ill_fl) const;
		void print_branch(string te, const vector <BranchEv> &branch) const;
		BranchMutInfo get_branch_mut_info(const vector <BranchEv> &branch, int sign) const;
		void set_branch_mut_info_after(const vector <BranchEv> &branch, vector <ObsMutSum> &obs_list, vector <AddMutSum> &add_list, int sign) const;
		void print_branch_mut_info(string te, const BranchMutInfo &mut_info) const;
		double branch_genetic_process(const vector <BranchEv> &branch) const;
		double branch_change_genetic_obs(const vector <BranchEv> &branch1, const vector <BranchEv> &branch2, const vector <BranchEv> &branch1_new, const vector <BranchEv> &branch2_new, vector <GenDifCha> &gen_dif_cha) const;
		double mean_branch_number() const;
		void trans_tree_swap_inf_ev_proposal(const BurnInfo &burn_info, vector <PropAcInfo> &prop_acc, bool use_time);
		void accept_iif_change(const vector <IIFChange> &iif_change);
		double iif_change_dmarkov(const vector <IIFChange> &iif_change) const;
		void trans_tree_swap_inf_proposal(const BurnInfo &burn_info, unsigned int &nfa, unsigned int &nac, unsigned int &ntr);
		unsigned int get_e_from_time(double t, const vector <Event> &ev) const;
		void update_inf_node_ref(unsigned int k);
		InfNodeAlter no_alter(unsigned int n) const;
		vector <unsigned int> get_not_obs(const vector <unsigned int> &obs_node1, const vector <unsigned int> &obs_node2, unsigned int N) const;
		void trans_tree_mut_proposal(const BurnInfo &burn_info, unsigned int &nac, unsigned int &ntr, double &si);
		void trans_tree_mut_origin_proposal(const BurnInfo &burn_info, unsigned int &nac, unsigned int &ntr, double &si);
		void trans_tree_mut_sample_proposal(const BurnInfo &burn_info, unsigned int &nac, unsigned int &ntr);
		void trans_tree_mut_local_proposal(const BurnInfo &burn_info, vector <GeneticProp> &gen_mut_info);
		void print_all_node() const;
		void print_node(string te, unsigned int n) const;
		void print_inf_origin() const;
		double round_down(double t) const;
		unsigned int inf_node_ti_from_tdiv(double tdiv) const;
		void disconnect_inf_node(unsigned int n);
		void connect_inf_node(unsigned int n, const IndInfFrom &iif, unsigned int mut_num_store);
		
		void check_pop_ind(const vector < vector <PopInd> > &pop_ind, const vector < vector <IndPop> > &ind_pop);
		NodeRef get_ind_noderef(unsigned int p, unsigned int i, double t) const;
		TransTreeChange no_trans_tree_change()  const;
		TransTreeChange trans_tree_fail() const;
		GenObsFunc get_genetic_obs_func_from(const vector <unsigned int> &obs_node, unsigned int n) const;
		GenObsFunc get_genetic_obs_func_to(const vector <unsigned int> &obs_node, unsigned int n, const NodeRef &nr_add, double t_add, TransTreeChange &gc, unsigned int &mutBC_num) const;
		double calc_genetic_obs_fast(int mutB_num, int mutA_num, const GenObsFunc &gof) const;
		TransTreeChange trans_tree_change(TransTreeChaType type, unsigned int n, const InfNodeAlter &node_alter, const IndInfFrom &iif_add, RegraftStore rs = REGRAFT_NO_STORE) const; 
		void accept_gen_dif_cha(const vector <GenDifCha> &gen_dif_cha);
		void gen_change_update(const TransTreeChange &gc);
		void remove_node(unsigned int n);
		void update_inf_ev(unsigned int n);
		void calculate_popnum_ind();
		void add_popnum_ind(unsigned int p, unsigned int i);
		void update_popnum_ind(unsigned int p, unsigned int i);
		double steadystate_genetic_obs();
		void popnum_ind_recalc_w(unsigned int p, unsigned int i);
		void change_add_node(double t, const NodeRef &nr_add, TransTreeChange &gc, RegraftStore rs = REGRAFT_NO_STORE) const;
		void change_remove_node(double t, const NodeRef &nr_from, TransTreeChange &gc, RegraftStore rs = REGRAFT_NO_STORE) const;
		void get_regraft_minimum(int &i, int &j, const GenObsFunc &gof) const;
		RGsamp get_regraft_sampler(const GenObsFunc &gof) const;
		void output_sampler_approximation(const RGsamp &rgs, const GenObsFunc &gof) const;
		bool change_add_node_with_gobs(unsigned int n, double t, const NodeRef &nr_add, TransTreeChange &gc, const GenObsFunc &gof, RegraftStore rs = REGRAFT_NO_STORE) const;
		bool change_remove_node_with_gobs(unsigned int n, double t, const NodeRef &nr_from, TransTreeChange &gc, const GenObsFunc &gof, RegraftStore rs = REGRAFT_NO_STORE) const;
		void binomial_split(int &B0, int &A0, double r1, double r2, string &warn) const;
		void change_add_node_gen_obs_max(unsigned int, double t, const NodeRef &nr_add, TransTreeChange &gc, const GenObsFunc &gof) const;
		vector <unsigned int> get_gen_obs_list(unsigned int n, unsigned int j, unsigned int n_not_allow=UNSET) const;
		void get_gen_obs_list2(unsigned int n, unsigned int j, vector <unsigned int> &list, unsigned int n_not_allow) const;
		void output_gen_dif() const;
		void output_waifw(string file) const;
		vector <double> recalculate_iif_w();
		void recalculate_iif_w_undo(const vector <double> &vec);
		vector <double> recalculate_popnum_ind_w(unsigned int po);
		void recalculate_popnum_ind_w_undo(unsigned int po, const vector <double> &vec);
		double harmomic_contribution(double var) const;
		double calc_genetic_obs_expected();
		void test_get_genetic_obs_func_from();
		void check_maximise_gen_obs_off(bool on);
		void check_gen_values(string te) const;
		void check_regraft();
		double get_infector_prediction_accuracy() const;
		 
	// In state_local.cc
	public:
		Like calculate_local_change(unsigned int p, vector <LocalChange> &local_change, bool &ill);
	
	private:
		void initialise_local_map(); 
		
		vector < vector <int> > comp_local_map;    // A map for local changes in compartment [p][c]
		vector < vector <double> > pop_trans_map;  // A map for changes in pop_trans_data [p][ref]
		vector < vector <double> > pop_data_map;   // A map for changes in pop_data [p][ref]
	
	// In like_integrate
	public:
		double like_integrate(const PV &param_value);
	
	// In state_check.cc
	public:
		void check_trans_num(string ref);
		void check_simp(string ref);
		void check(string ref);
		//void scan_variable(string name, double min, double max);
		void print_ev_data(string te, const vector <EventData> &event_data, unsigned int p) const;
		void check_pop_t(string ref);
		void check_precalc_eqn(string ref);
		void check_neg_rate(string name);
		void check_final_li_wrong();
		void check_markov_value_dif();
		void add_alg_warn(string te);
		void check_precalc_dif(string ref);
		void scan_param();
		void check_dpop();
		void check_regenerate(string ref);
		void check_inf_node_times(string ref) const;
		void check_markov_inf(string ref) const;
		void check_genetic_process();
		void check_iif(string ref) const;
		void check_back(string ref) const;
		void check_event_in_range(string ref) const;
	
	private:
		void check_dependent_param(string ref);
		void check_ref(unsigned int p, string ref);
		void check_markov_div_value(unsigned int p, string ref);
		void check_markov_trans(unsigned int p, string ref);
		void check_nm_trans(unsigned int p, string ref);
		void check_prior(string ref);
		void check_like(string ref);
		//void check_spline(string ref);
		void check_ie(unsigned int p, string ref);
		void check_pop_like(unsigned int p, string ref);
		void check_cpop_st(string ref);
		void check_maps(string ref);
		void check_event_observed(unsigned int p, string ref);
		void check_obs_like(unsigned int p, string ref);
		void check_init_cond_like(unsigned int p, string ref);
		void check_init_cond_prior(string ref);
		void check_linearise();
		void check_para_speedup(unsigned int p, string ref);
		//void check_spline_store(string ref);
		void check_genetic_value(string ref);
		void check_popnum_ind(string ref);
		void check_add_move_rem(string ref);
		void check_effect_out_of_range();
		void output_dump() const;
		void print_iif() const;
	
	// In state_profiler.cc
	public:
		void profile_memory() const;
	
	// In consistent.cc	
	public:
		vector <IndRef> force_consistent_solution();
		void resample_ind(bool if_wrong);
		void no_valid_state(const vector < vector <IndRef> > &not_possible) const;
		
	private:
		double calculate_basic_rate(unsigned int p, unsigned int tr, unsigned int ti) const;
		CompProb comp_prob_sample(const vector <CompProb> &cp_pos) const;
		ConsistTrans contr_sample(const vector <ConsistTrans> &list, const vector <ConsistNode> &node) const;
		ConsistComp conco_sample(const vector <ConsistComp> &list, const vector <ConsistNode> &node) const;

	private:
		const Model &model;
};


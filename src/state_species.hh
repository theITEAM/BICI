#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "species.hh"
#include "model.hh"
#include "source_sampler.hh"

class StateSpecies                         // Stores information about the state
{
	public:
		SpeciesType type;                      // Type of model (individual / population-based)

		double dt;                             // The time step
		
		unsigned int T;                        // Total number of divisions
		
		vector < vector <double> > Li_markov;  // The Markovian likelihood [markov_eq][ti]
		vector <double> Li_ie;                 // The likelihood of the individual effects [ie_group]
		vector <double> Li_obs_ind;            // Individual observation likelihood [#ind]
		vector <double> Li_obs_pop;            // Observation likelihood for population data
		vector <double> Li_obs_pop_trans;      // Observation likelihood for pop trans data
		double Li_init_cond;                   // Likelihood from initial conditions
	
		// USED IN POPULATION-BASED MODELS
		
		vector < vector <double> > trans_num;  // The number of transitions occuring [tr][time]
		vector < vector <double> > tnum_mean_st;// The expected number of transitions [tr][time]
		vector < vector <double> > Li_markov_pop;// The Markovian likelihood [tr][time]
		vector < vector <double> > cpop_st;    // The compartmental population [time][c]
		
		InitCondValue init_cond_val;           // The population values for the initial state
		
		// USED FOR MBPS
		vector < vector <double> > trans_num_f;// The number of transitions occuring [tr][time]
		vector < vector <double> > tnum_mean_st_f;// The expected number of transitions [tr][time]
		vector < vector <double> > cpop_st_f;  // The compartmental population [time][c]
		
		vector <Individual> individual;        // Provides information about individuals
		//Hash hash_individual;
		
		// USED IN OBSERVATION MODELS
		vector <double> pop_data_num;          // Values at pop data measurements ([pop meas])
		vector <double> pop_trans_data_num;    // Values at pop trans data meas. ([pop trans meas])
		vector < vector <double> > pop_data_cgl;      // Populations at data measurements ([pd meas][clgl])
		vector < vector <double> > pop_trans_data_tgl;// Transitions at data measurements ([td meas][trgl])
		
		vector <double> obs_eqn_value;         // Observation equation values (e.g. Se)
		
		vector <string> unused_name;           // Keeps a list of unused names (for addition)
		
		vector < vector < vector <NMEventRef> > > nm_trans_ev_ref; // Stores which events undergo NM transitions // [m][ti][#]

		vector < vector < vector <unsigned int> > > nm_trans_incomp_ref; // Stores incomplete NM [n][ti][#]
		
		vector < vector <double> > obs_trans_eqn_value; // The log of probability of not observing trans [#][ti] 
		vector < vector <unsigned int> > obs_trans_eqn_num;// The number of transition [m][ti]
		
		double prob_trans_tree;                // Used to check trans tree likelihood correct
		
		// These values are dynanamically updates as a simulation is performed
		vector <double> cpop;                  // The current population in the compartments (population-based models)
		
		vector < vector <double> > X_order;    // An ordered version of Xvector value
		
		vector <IEGroupSampler> ind_eff_group_sampler;// Stores samplers for individual effects
		
		unsigned int N;                        // The number of Markovian equations
	
		vector <MarkovEqnVariation> markov_eqn_vari; // Store variation in markov equations
		
		SourceSampler source_sampler;          // Used to sampling add/rem ind 
		
		RatePosteriorMean rate_mean;// Store transition rate means (for local proposals)
		
		vector <BackInd> back_ind;             // Store information to restore state
		vector <Back> back;
		
		// USED TO SIMULATE INDIVIDUAL-BASED MODELS
		vector <unsigned int> ind_sim_c;       // Under simulation stores compartment state
	
		vector <SimTrigEventDiv> trig_div;     // Stores future non-Markovian events

		vector <double> timer;                 // General purpose timers
		
		vector <AlgWarn> alg_warn;             // Stores any algorithm warnings
		
		StateSpecies(PV &param_val, const vector <Equation> &eqn, const Species &sp, const Model &model, const vector <unsigned int> &pop_affect_, Operation mode, const double &dif_thresh);
	
		void simulate_init();
		unsigned int get_cinit_to_use(vector <unsigned int> &cinit_to_use) const;
		void simulate_individual_init();
		void add_ind_source(unsigned int i, const vector <SourceSamp> &source_samp);
		TRange source_time_range(const IndData &ind) const;
		void error_load_sample(unsigned int num) const;
		void simulate_sample_init(unsigned int ti, const SampleSpecies &samp_sp);
		void set_tnum_mean(unsigned int ti_end, const vector < vector <double> > &popnum_t);
		void update_population_based(unsigned int ti, bool stoc, const vector <double> &popnum, const vector <double> &val_fast);
		vector <double> set_exp_fe(unsigned int f);
		void set_exp_fe_restore(unsigned int f, const vector <double> &store);
		vector <double> set_exp_ie(Individual &ind) const;
		vector <double> recalculate_exp_ie(unsigned int ie);
		void recalculate_exp_ie_restore(unsigned int ie, const vector <double> &store);
		void generate_A();
		void check_precalc_num(unsigned int n);
		void mbp(double sim_prob, vector < vector <double> > &popnum_t, const MBPfast &mbp_value_linear);
		void val_fast_update(unsigned int ti, vector <double> &val_fast, const vector < vector <double> > &popnum_t, const vector < vector <double> > &pop_grad, const LinearForm &lin_form) const;
		void factor_nopop_change(const vector <unsigned int> &list, vector < vector <double> > &val_store, const LinearForm &lin_form, bool &set) const;
		void linear_form_calculate(vector < vector <double> > &val_store, const vector <unsigned int> &list, const LinearForm &lin_form, const vector < vector <double> > &popnum_t) const;
		vector < vector <double> > pop_grad_calc(const LinearForm &lin_form) const;
		void mbp_accept(double &like_ch);
		vector <double> calc_val_fast_init(const LinearForm &lin_form, const vector < vector <double> > &pop_grad, const vector <double> &popnum) const;
		void set_tnum_mean_st_f(vector < vector <double> > &tnum_mean_st_f, unsigned int ti, const vector <double> &popnum, const vector <double> &cpop, double dt, const vector <double> &val_fast, const MBPfast &mbp_fast) const;
		unsigned int move_event(vector <Event> &ev, unsigned int index, double t_new) const;
		void make_consistent(vector <Event> &event) const;
		void set_cpop_st();
		vector < vector <double> > ibm_cpop_st() const;
		double get_trans_obs_prob(unsigned int trg, const ObsData &ob) const;
		void compare_covar(string te, const vector < vector <double> > &omega, unsigned int g) const;
		ListMove get_ev_link_listmove(unsigned int e, const vector <Event> &event, const Individual &ind, const vector < vector <double> > &popnum_t) const;
		TRange ev_link_trange(const ListMove &lm, const vector <Event> &event) const;
		bool ev_link(const Event &ev, const Individual &ind, const vector < vector <double> > &popnum_t) const;
		bool source_ind(unsigned int i) const;
		unsigned int source_num(unsigned int min, unsigned int max) const;
		vector < vector <InfPeriod> > get_inf_period(const vector <unsigned int> &ref) const;
		void tnum_mean_st_back(unsigned int tr, unsigned int ti1, unsigned int ti2);
		void cpop_st_update(unsigned int ti, unsigned int ti_next, unsigned int c, int ma);
		
		double get_indfac(const Individual &ind, const  MarkovEqn &mar_eqn) const;
		void update_ind_basic(const vector < vector <Event> > &ev_new);
		void calculate_tnum_mean_para(vector <double> &tnum_mean, vector <unsigned int> list, unsigned int i, const vector < vector <double> > &popnum_t, const vector < vector <double> > &cpop_t, double dt) const;
		double calculate_tnum_mean(unsigned int ti, unsigned int tr, const vector <double> &popnum, const vector <double> &cpop, double dt) const;
		
	private:	
		unsigned int ti_sort;                  // Trigger events which have been sorted
		
		unsigned int nnode;                    // The number of nodes
		vector <double> markov_tree_rate;      // Rates going up Markov tree
	
		void initialise_arrays();
		void reset_arrays();
		vector <double> sample_ie() const;
		void sample_ie_Amatrix();
		void ie_Amatrix_sampler_init();
		vector <double> calculate_tnum_mean_all(unsigned int ti, const vector <double> &popnum, const vector <double> &cpop, double dt) const;
		vector <double> calculate_tnum_mean_fast(unsigned int ti, const vector <double> &popnum, const vector <double> &cpop, double dt, const vector <double> &val_fast) const;
		vector <double> sample_trans_num(const vector <double> &tnum_mean, bool stochastic) const;
		void update_cpop(unsigned int ti, vector <double> &cpop, const vector <double> &tnum) const;
		bool enter_flat_dist(const IndData &ind) const;
		unsigned int ind_sample_init_c(const IndData &ind) const;
		
	// In 'state_species_ind.cc'
	public:
		void activate_initial_state(double t,const vector < vector <double> > &popnum_t);
		void update_individual_based(unsigned int ti, const vector < vector <Poss> > &pop_ind, const vector < vector <double> > &popnum_t, const vector <double> &val_fast);
		vector <double> calculate_omega(unsigned int g);
		void calculate_omega_restore(unsigned int g, const vector <double> &store);
		void remove_me_ev(const Event &ev);
		vector <double> likelihood_indfac_int();
		void likelihood_indfac_int_restore(const vector <double> &store);
		void likelihood_ind_trans();
		unsigned int add_individual(IndType ind_type, string name="");
		void remove_individual(unsigned int i, vector <InfNode> &inf_node);
		//string get_new_ind_name(string pre, unsigned int num) const;
		void calculate_N_unobs();
		void calculate_indfac_sum();
		Event get_event(EventType type, unsigned int i, unsigned int tr_gl, unsigned int move_c, unsigned int cl, unsigned int c_after, double t, const IndInfFrom &inf_from);
	
	private:
		void add_data_event(unsigned int i, double t, const vector <SimTrigEvent> &trig_vec, const vector < vector <double> > &popnum_t, const vector < vector <Poss> > &pop_ind);
		bool allow_event(double t, const IndTransRef &itr) const;
		void ie_sampler_init();
		void update_ind_c(unsigned int i, double t, unsigned int cl_trans, const vector < vector <double> > &popnum_t);
		void insert_trigger_event(const SimTrigEvent &trig);
		SimTrigEvent get_nm_trig_event(double t, unsigned int i,  unsigned int c, unsigned int cl, const vector < vector <double> > &popnum_t);
		bool try_insert_data_trans_event(double t, unsigned int cl, unsigned int i);
		void add_markov_transition(unsigned int i, unsigned int tgl);
		void markov_eqn_recalc(unsigned int e, unsigned int ti, const vector <double> &popnum);
		void markov_eqn_recalc_fast(unsigned int ti, const vector <double> &popnum, const vector <double> &val_fast);
		void update_markov_eqn_value(unsigned int ti, const vector < vector <double> > &popnum_t, const vector <double> &val_fast);
		void sort_trig_event(unsigned int ti);
		void sample_infecting_ind(unsigned int i, double t, unsigned int tr_gl, IndInfFrom &inf_from, const vector < vector <double> > &popnum_t, const vector < vector <Poss> > &pop_ind);
		IndInfFrom get_waifw(unsigned int i, double t) const;
		void update_ind_trans(double t, const IndTransRef &itr, const vector < vector <double> > &popnum_t, const vector < vector <Poss> > &pop_ind);
		void remove_markov_trans(unsigned int i);
		void update_ind_remove(unsigned int i, double t);
		void update_ind_move(unsigned int i, double t, unsigned int c_comp, unsigned int cl, const vector < vector <double> > &popnum_t);
		void update_markov_tree_rate(unsigned int e, double dif);
		void set_markov_tree_rate();
		void add_event(EventType type, unsigned int i, unsigned int tr_gl, unsigned int move_c, unsigned int cl, unsigned int c_after, double t, const IndInfFrom &inf_from);
		
	// In 'state_species_check.cc'
	public:
		void check(unsigned int ti, const vector < vector <double> > &popnum_t);
		void print_markov_eqn(unsigned int i) const;
		void print_ind(unsigned int i) const;
		void ensure_consistent(vector <Event> &ev) const;
		string print_event(string name, const Individual &ind) const;
		string print_event(const vector <Event> &event, bool str=false) const;
		void check_cl_ev(Individual ind) const;
		void set_ind_sim_c(unsigned int ti);
		bool all_events_correct(unsigned int i, unsigned int cl);
		
	private:
		void check_markov_eqn_ref() const;
		void check_cpop() const;
		void check_markov_indfac_sum();
		void check_markov_eqn(unsigned int ti, const vector <double> &popnum) const;
		void check_markov_tree_rate();
		void check_event() const;
		void check_erlang() const;
		void check_trans_num_neg() const;
		void add_alg_warn(string te);
		string check_prior(const Equation &eq) const;
		
		PV &param_val;            // Quantities reference from state
		const vector <Equation> &eqn;
		const Species &sp;                     
		const Model &model;
		const double &dif_thresh;
	
		vector <unsigned int> pop_affect;
		Operation mode;                        // SIM for simulation and INF for inference	
		
	// In 'state_species_like.cc'
	public:
		vector <double> likelihood_markov(unsigned int e, const vector <unsigned int> &list, double &like_ch);
		vector <double> likelihood_nm_trans(unsigned int m, const vector <unsigned int> &list, const vector < vector <double> > &popnum_t, double &like_ch);
		void likelihood_nm_trans_restore(unsigned int m, const vector <unsigned int> &list, const vector <double> store);
		vector <double> likelihood_nm_trans_incomp(unsigned int n, const vector <unsigned int> &list, const vector < vector <double> > &popnum_t, double &like_ch);
		void likelihood_nm_trans_incomp_restore(unsigned int m, const vector <unsigned int> &list, const vector <double> store);
		vector <double> likelihood_nm_trans_bp(unsigned int m, const vector <unsigned int> &list, const vector < vector <double> > &popnum_t, double &like_ch);
		void likelihood_nm_trans_bp_restore(unsigned int m, const vector <unsigned int> &list, const vector <double> store);
		double nm_trans_incomp_full_like(const vector <unsigned int> &nmtrans_ref, double dtdiv, double dt, const vector< vector <double> > &ref_val, const vector <double> &bp_val) const;
		double nm_trans_incomp_like(TransType type, double dtdiv, double dt, const vector <double> &ref_val) const;
		double nm_trans_incomp_like_no_log(TransType type, double dtdiv, double dt, const vector <double> &ref_val) const;
		vector <double> likelihood_markov_value(unsigned int e, const vector <unsigned int> &list, const vector < vector <double> > &popnum_t);
		vector <double> likelihood_markov_value_nopop(const vector <unsigned int> &me_list, const vector <unsigned int> &list, const vector < vector <double> > &popnum_t, const PV &param_val);
		void likelihood_markov_value_nopop_restore(const vector <unsigned int> &me_list, const vector <unsigned int> &list, const vector<double> &store);
		vector <double> likelihood_markov_value_linear(const vector <unsigned int> &list, const LinearForm &lin_form, const vector < vector <double> > &popnum_t);
		void likelihood_markov_value_linear_restore(const vector <unsigned int> &list, const LinearForm &lin_form, const vector<double> &store);
		MeanSD get_mean_sd(TransType type, const vector <double> &ref_val) const;
		vector <NMupdate> likelihood_ie_nm_trans_change(unsigned int i, unsigned int ie, double factor, const vector < vector <double> > &popnum_t, double &like_ch);
		void likelihood_ie_nm_trans_update(const vector <NMupdate> &nm_st);
		double likelihood_indeff_group(unsigned int g, double &like_ch);
		double likelihood_indeff_change(unsigned int i, unsigned int ie, double value_new, double value_old);
		void update_indfac_int(unsigned int c, double ti, double tf, const Individual &ind, double sign, double &like_ch);
		vector <double> likelihood_ie_change(unsigned int i, unsigned int ie, double factor, Like &like_ch);
		void likelihood_ie_change_restore(unsigned int i, unsigned int ie, double factor, const vector <double> &store);
		void likelihood_pop(const vector < vector <double> > &popnum_t);	
		void likelihood_pop_section(unsigned int tr, unsigned int ti1, unsigned int ti2, const vector < vector <double> > &popnum_t, const vector <PopChange> &pop_change, double &like_ch);
		void likelihood_pop_spline_section(unsigned int tr, unsigned int ti1, unsigned int ti2, const vector < vector <double> > &popnum_t, double &like_ch);
		void likelihood_ib_spline_section(unsigned int ee, unsigned int ti, unsigned int ti_next, const vector < vector <double> > &popnum_t, double &like_ch);
		double Li_update_tn(unsigned int tr, unsigned int ti, int sign);
		double Li_update_c(unsigned int c, int ma, unsigned int ti, unsigned int ti_next, const vector < vector <double> > &popnum_t);
		vector <double> likelihood_pop_change(unsigned int tr, const vector <unsigned int> &list, const vector < vector <double> > &popnum_t, double &like_ch);		
		void calculate_markov_Li(unsigned int tr, const vector <unsigned int> &list, double &like_ch, vector <double> &store);
		void likelihood_pop_change_restore(unsigned int tr, const vector <unsigned int> &list, const vector <double> &vec);
		vector <double> likelihood_pop_change_nopop(const vector <unsigned int> &list, const EqnNoPop &eq_nopop, double &like_ch);
		void likelihood_pop_change_nopop_restore(const vector <unsigned int> &list, const EqnNoPop &eq_nopop, const vector<double> &store);
		vector <double> likelihood_pop_change_linear(const vector <unsigned int> &list, const LinearForm &lin_form, const vector < vector <double> > &popnum_t, double &like_ch);
		bool likelihood_pop_change_linear_factor_nopop(const vector <unsigned int> &list, const LinearForm &lin_form, double &like_ch, vector <double> &store);
		void likelihood_pop_change_linear_restore(const vector <unsigned int> &list, const LinearForm &lin_form, const vector<double> &store);
		vector <double> likelihood_init_cond(double &like_ch);
		NMIncompVal get_nm_incomp_val(const NMTransIncomp &nmti, unsigned int ti, const Individual &ind, const vector < vector <double> > &popnum_t) const;
		void print_incomp_val(string te, NMIncompVal value) const;
		NMIncompVal get_nm_incomp_val_ie_factor(unsigned int ie, double factor, const NMTransIncomp &nmti, unsigned int ti, const Individual &ind, const vector < vector <double> > &popnum_t) const;
		double nm_obs_dprob(const Individual &ind) const;
		double nm_single_obs_dprob(unsigned int cl, const Individual &ind) const;
		
		double sum_markov_prob(double t1, double t2, unsigned int c, unsigned int tr_gl, unsigned int i, vector < vector <double> > &en) const;
		void calc_trans_diag(ParticleSpecies &ps, const vector < vector <double> > &popnum_t);
		void add_nm_prob(unsigned tr_gl_or, double t1, ParticleSpecies &ps, const vector < vector <double> > &popnum_t, const Individual &ind) const;
		void add_nm_cpd(unsigned tr_gl_or, double t1, double t2, ParticleSpecies &ps, const vector < vector <double> > &popnum_t, const Individual &ind) const;
	
	// In 'state_species_like_obs.cc'
	public:
		vector <double> calculate_obs_eqn(const vector <unsigned int> &list);
		vector <double> likelihood_obs_ind(const vector <unsigned int> &list, double &like_ch);
		double obs_mod_probability(double num, ObsModelVariety type, double value, double obs_mod_val) const;
		vector <double> likelihood_obs_pop(const vector <unsigned int> &list, double &like_ch);
		double calculate_pop_num(unsigned int i, double t, const vector <unsigned int> &c_nonzero, const vector <unsigned int> &comp_obs_mod_ref) const;
		vector <double> likelihood_obs_pop_trans(const vector <unsigned int> &list, double &like_ch);
		vector <double> calculate_pop_data_cgl_trgl();
		void restore_pop_data_cgl_trgl(const vector <double> &store);
		void likelihood_observation_restore(const vector <double> &store);
		double like_diag_test(unsigned int c, const ObsData &ob) const;
		double like_trans_prob_obs(unsigned int tr, const ObsData &ob) const;
		double like_comp_obs(unsigned int c, const ObsData &ob) const;
		double like_ind_init_obs(unsigned int c, unsigned int entref) const;
		double update_pop_trans_data(double t, unsigned int tr_gl, int sign);
		double update_pop_trans_data2(unsigned int ti, unsigned int tr_gl, int sign);
		double update_single_pop_trans_data(double dnum, unsigned int ref);
		double update_single_pop_data(double dnum, unsigned int ref);
		double update_pop_data(unsigned int i, unsigned int c_old, unsigned int c_new);
		vector <double> likelihood_unobs_trans(unsigned int trg, const vector <unsigned int> &list, double &like_ch);
		void set_all_obs_trans_value();
		void alter_obs_trans_eqn(unsigned int tr_gl, double t, int sign, bool observed, const vector <ObsData> &obs_data, double &like_ch);
		//void alter_obs_trans_eqn(unsigned int tr_gl, double t, int sign, double &like_ch);
		void calculate_obs_trans_eqn_num();
		double trans_logprob_notobs(unsigned int trg, unsigned int ti) const;
		void set_event_observed();
		
	// In 'species_update_ind.cc'
		vector <unsigned int> update_ind(unsigned int i, vector <Event> &event_new, const vector < vector <double> > &popnum_t, vector <PopUpdate> &pop_update, Like &like_ch);	
		void add_event_ref(unsigned int i, unsigned int ee,  const vector < vector <double> > &popnum_t, Like &like_ch);
		void set_m_ti_origin(vector <Event> &ev_new) const;
		void remove_event(Event &ev, const Individual &ind, Like &like_ch, unsigned int i, vector <Event> &event_old);
		void remove_all_event_ref(const vector <unsigned int> &ind_list);
		void add_all_event_ref(const vector <unsigned int> &ind_list, const vector < vector <Event> > &ev_new);
		void remove_all_nmincomp_ref(const vector <unsigned int> &ind_list);
		void add_all_nmincomp_ref(const vector <unsigned int> &ind_list, const vector < vector <IncompNMTransRef> > &incomp_ref_new);
		void adjust_incomp_ref(vector <IncompNMTransRef> &incomp_ref, const vector <Event> &ev);
		void set_e_origin(Event &ev,  unsigned int e, const vector <Event> &event);
		void set_incomp_ref(unsigned int i, const vector < vector <double> > &popnum_t, double &like_ch);
		void setup_nm_trans(const vector < vector <double> > &popnum_t);
		void incomp_turn_on(unsigned int i, unsigned int n, unsigned int ti, double dtdiv, unsigned int e_begin, double t_end, IncompNMTransRef &inm);
		void incomp_turn_off(IncompNMTransRef &inm);
		
		void likelihood_init_cond_change(unsigned int c_enter_old, unsigned int c_enter_new, Like &like_ch);
		void init_cond_change_back(unsigned int c_enter_old, unsigned int c_enter_new);
		void recalc_markov_value(unsigned int ee, unsigned int ti, unsigned int ti_next, const vector < vector <double> > &popnum_t, const vector <PopChange> &pop_change, double &like_ch);
		void restore_back();
		
	// In state_species_local.cc
		vector <LocalIndChange> local_ind_change(unsigned int i, unsigned int cl) const;
		void print_local_ind_change(const LocalIndChange &lich) const;
		double create_local_change(double &timefac, const LocalIndChange &lich, vector <Event> &ev_new, Event &enew, unsigned int nind_obs, const vector <Event> &ev_rev, LocalDir dir) const;
		void rate_store_init();
		void update_rate_mean(const vector < vector <double> > &popnum_t);
		void swap_check_rev_prob_try(double prob_try, double timefac, const LocalIndChange &lich, const vector <LocalIndChange> &add, const vector <Event> &ev_store) const; 
		bool tr_swap_same(const vector <TrSwap> &tr_swap, const vector <TrSwap> &tr_swap2) const;
};


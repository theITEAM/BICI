#ifndef BICI__STATE_SPECIES_HH
#define BICI__STATE_SPECIES_HH

#include <vector>

using namespace std;

#include "struct.hh"
#include "species.hh"
//#include "model.hh"

class StateSpecies                                             // Stores information about the state
{
	public:
		SpeciesType type;                      // Type of model (individual / population-based)

		vector < vector <double> > Li_markov;  // The Markovian likelihood [markov_eq][ti]
		vector <double> Li_nm_trans;           // The non-Markovian transition likelihood	[nm_trans]
		vector <double> Li_ie;                 // The likelihood of the individual effects [ie_group]
		double Li_obs;                         // The observation likelihood
		
		// USED IN POPULATION-BASED MODELS
		vector <double> cpop_init;             // The initial population in the compartments 
		
		vector < vector <double> > trans_num;  // The number of transitions occuring [tr][time]
		vector < vector <double> > tnum_mean_st;// The expected number of transitions [tr][time]
		vector < vector <double> > Li_markov_pop;// The Markovian likelihood [tr][time]
		vector < vector <double> > cpop_st;    // The compartmental population [time][c]
		
		// USED FOR MBPS
		vector < vector <double> > trans_num_f;  // The number of transitions occuring [tr][time]
		vector < vector <double> > tnum_mean_st_f;// The expected number of transitions [tr][time]
		vector < vector <double> > cpop_st_f;    // The compartmental population [time][c]
		
		// USED IN INDIVIDUAL-BASED MODELS
		vector <Individual> individual;        // Provides information about individuals
		
		vector < vector <NMTransInd> > nm_trans_ind; // Stores which individuals undergo NM transitions
		
		// USED IN OBSERVATION MODELS
		vector <double> pop_data_num;              // These correspond to the population data measurements
		vector <double> pop_trans_data_num;        // These correspond to the pop trans data measurements
		
		vector <double> obs_eqn_value;             // Observation equation values (e.g. Se)
		
		// These values are dynanamically updates as a simulation is performed
		vector <double> cpop;                        // The current population in the compartments (population-based models)
		
		vector < vector <double> > X_order;          // An ordered version of Xvector value
		
		vector <IEGroupSampler> ind_eff_group_sampler;// Stores samplers for individual effects
		
		unsigned int N;                               // The number of Markovian equations
	
		vector <MarkovEqnVariation> markov_eqn_vari;  // Store variation in markov equations
		
		// USED TO SIMULATE INDIVIDUAL-BASED MODELS
		vector <TrigEventDiv> trig_div;               // Stores future non-Markovian events
		vector <TrigEvent> end_nm_event;              // Stores non-Markovian events after end
		
	
		StateSpecies(const vector <double> &param_val, const vector <SplineValue> &spline_val, const vector <Equation> &eqn, const vector <Param> &param, const vector <ParamVecEle> &param_vec, const Species &sp, const Details &details, const	vector <double> &timepoint, const	vector <double> &dtimepoint,  const vector <unsigned int> &pop_affect, Operation mode);
		
		void simulate_init(const InitCondValue &initc_val);
		void update_population_based(unsigned int ti, bool stoc, const vector <double> &popnum);
		vector <double> set_exp_fe(unsigned int f);
		void set_exp_fe_restore(unsigned int f, const vector <double> &store);
		vector <double> set_exp_ie(Individual &ind) const;
		vector <double> recalculate_exp_ie(unsigned int ie);
		void recalculate_exp_ie_restore(unsigned int ie, const vector <double> &store);
		void generate_A();
		void mbp(double sim_prob, vector < vector <double> > &popnum_t);
		void mbp_accept(double &like_ch, const vector < vector <double> > &popnum_t);
	
	private:	
		unsigned int ti_sim;                                // Indexes time (used in simulation)
		
		unsigned int T;                                     // The maximum time
		
		unsigned int nnode;                                 // The number of nodes
		vector <double> markov_tree_rate;                   // Rates going up Markov tree
	
		void initialise_arrays();
		vector <double> sample_ie() const;
		void sample_ie_Amatrix();
		void ie_Amatrix_sampler_init();
		vector <double> calculate_tnum_mean(unsigned int ti, const vector <double> &popnum, const vector <double> &cpop, const vector <double> &param_val, double dt) const;
		double calculate_tnum_mean(unsigned int ti, unsigned int tr, const vector <double> &popnum, const vector <double> &cpop, const vector <double> &param_val, double dt) const;
		vector <double> sample_trans_num(const vector <double> &tnum_mean, bool stochastic) const;
		void update_cpop(vector <double> &cpop, const vector <double> &tnum) const;
		void add_individual(string name, unsigned int c, bool from_source);
		double get_indfac(const Individual &ind, const  MarkovEqn &mar_eqn) const;
		string get_new_ind_name(string pre, unsigned int num) const;
		
	// In 'state_species_ind.cc'
	public:
		void activate_initial_state(const vector <double> &popnum);
		void update_individual_based(unsigned int ti_, const vector <double> &popnum);
		vector <double> calculate_omega(unsigned int g);
		void calculate_omega_restore(unsigned int g, const vector <double> &store);
		void remove_me_ev(const Event &ev);
		vector <double> likelihood_indfac_int();
		void likelihood_indfac_int_restore(const vector <double> &store);
		void likelihood_ind_trans();
		
	private:
		void ie_sampler_init();
		void update_ind_c(unsigned int c, unsigned int i, double t, const vector <double> &popnum);
		void insert_trigger_event(const TrigEvent &trig);
		void add_markov_transition(unsigned int i, unsigned int tgl);
		void update_markov_eqn_value(const vector <double> &popnum);
		void sort_trig_event();
		void update_ind_trans(double t, const IndTransRef &itr, const vector <double> &popnum);
		void remove_markov_trans(unsigned int i);
		void update_ind_remove(unsigned int i, double t);
		void update_ind_move(unsigned int i, double t, unsigned int c_comp, unsigned int cl, const vector <double> &popnum);
		void update_markov_tree_rate(unsigned int e, double dif);
		void add_event(EventType type, unsigned int i, unsigned int tr_gl, unsigned int c_gl, unsigned int c_after, double t);
		
	// In 'state_species_check.cc'
	public:
		void check(const vector <double> &popnum) const;
		void print_markov_eqn(unsigned int i) const;
		void print_ind(unsigned int i) const;
		void ensure_consistent(unsigned int c, vector <Event> &ev) const;
		void recalculate_cl_ev(Individual &ind) const;
		void print_event(string name, const Individual &ind) const;
		unsigned int move_event(unsigned int cinit, vector <Event> &ev, unsigned int index, double t_new) const;
		void check_cl_ev(Individual ind) const;
		void check_transnum();
		
	private:
		void check_markov_eqn_ref() const;
		void check_cpop() const;
		void check_markov_eqn(const vector <double> &popnum) const;
		void check_markov_tree_rate() const;
		void check_event() const;
		void check_ind_ref() const;
		string check_prior(const Equation &eq) const;
		
		const vector <double> &param_val;                    // Quantities reference from state
		const vector <SplineValue> &spline_val;            // Values for the spline
		const vector <Equation> &eqn;
		const vector <Param> &param;
		const vector <ParamVecEle> &param_vec;
		const Species &sp;                           
		const Details &details;
		const vector <double> &timepoint;
		const vector <double> &dtimepoint;
		vector <unsigned int> pop_affect;
		Operation mode;                                    // SIM for simulation and INF for inference	
		
	// In 'state_species_like.cc'
	public:
		vector <double> likelihood_markov(unsigned int e, vector <unsigned int> list, AllTimeStatus as, double &like_ch);
		double likelihood_nm_trans(unsigned int i, const vector < vector <double> > &popnum_t, double &like_ch);
		vector <double> likelihood_markov_value(unsigned int e, vector <unsigned int> list, AllTimeStatus as, const vector < vector <double> > &popnum_t);
		double likelihood_indeff_group(unsigned int g, double &like_ch);
		double likelihood_indeff_change(unsigned int i, unsigned int ie, double value_new, double value_old);
		void update_indfac_int(unsigned int c, double ti, double tf, const Individual &ind, double sign, double &like_ch);
		vector <double> likelihood_ie_change(unsigned int i, unsigned int ie, double factor, Like &like_ch);
		void likelihood_ie_change_restore(unsigned int i, unsigned int ie, double factor, const vector <double> &store);
		void remove_ind_events_markov_eqn(unsigned int i, double &like_ch);
		void add_ind_events_markov_eqn(unsigned int i, double &like_ch);
		void likelihood_pop(const vector < vector <double> > &popnum_t);	
		void likelihood_pop_section(unsigned int tr, unsigned int ti1, unsigned int ti2, const vector < vector <double> > &popnum_t, double &like_ch);
		vector <double> likelihood_pop_change(unsigned int tr, const vector <unsigned int> &list, const vector < vector <double> > &popnum_t, double &like_ch);		
		void likelihood_pop_change_restore(unsigned int tr, const vector <unsigned int> &list, const vector <double> &vec);
		vector <double> likelihood_observation(double &like_ch);
		void likelihood_observation_restore(const vector <double> &store);
		double like_diag_test(unsigned int c, const ObsData &ob) const;
		double update_pop_trans(double t, unsigned int tr_gl, int sign);
};

#endif

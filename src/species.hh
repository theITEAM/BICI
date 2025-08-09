#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "equation.hh"
#include "hash.hh"

class Species   	                         // Stores details of a species
{ 	     
	public:
		Species(unsigned int p, const vector <double> &timepoint, const Details &details);
	
		unsigned int p_species;                // The number of the species
	 	
		string name;                           // The name of the species
		SpeciesType type;                      // Type of model (individual / population-based)

		vector <Classification> cla;           // Stores classification within the species
		unsigned int ncla;                     // Number of classifications

		unsigned int N;                        // The total number of compartments in species

		vector <IndEffect> ind_effect;         // Any individual effects within species
		
		vector <IEgroup> ind_eff_group;        // Groups of individual effects
		
		vector <FixedEffect> fix_effect;       // Any fixed effects within species
		
		bool ind_variation;                    // Set if ind_effect or fix_effect 
		
		vector <unsigned int> comp_mult;       // Used to convert from classification definition to global

		vector <CompGlobal> comp_gl;           // Global compartments

		vector <TransGlobal> tra_gl;           // Global transitions  
		
		bool tform_set;                        // Determines if tform is set  
		
		vector < vector <int> >	tr_shift;      // If tform not set this is used  
		
		vector < vector < vector <unsigned int> > >  trg_from_tr; // Converts from tr to trg [cgl][cl][tr]
		
		// Convert from a global compartment to another based on classification cl being in c
		vector < vector < vector <unsigned int> > > comp_global_convert; // [cgl][cl][c]
	
		bool comp_global_convert_set;
	
		vector < vector <bool> > tra_ie_affect;// Set if an ie affects a transition
		
		vector < vector <bool> > tra_fe_affect;// Set if an fe affects a transition
		
		vector <MarkovEqn> markov_eqn;         // List Markovian equations
			
		MarkovTree markov_tree;                // Constructs a tree of nodes to sample markov equations
			
		vector <bool> cgl_begin_nm;            // Determines if glob comp starts nm 
		
		bool trans_tree;                       // Determines if the transmission tree is turned on
		
		bool nm_flag;                          // Determines if non-Markovian transitions
		
		bool multi_trans_cl;                   // Set if transitions in more than one classification
		
		vector < vector <bool> > cgl_incomp_nmtrans_cl;// List if nm [cgl][#]
		
		unsigned int infection_cl;             // Determines which classification represents infection status
	
		bool contains_source;                  // Set to true if a source transition exists
		
		bool contains_sink;                    // Set to true if a sink transition exists
		
		vector <unsigned int> cgl_tr_source;   // Gives the transition entering  
		
		vector <unsigned int> cgl_tr_sink;     // Gives the transition leaving 
		
		unsigned int nindividual_in;           // The number of individuals in system
		unsigned int nindividual_obs;          // The number of observed individuals (including in A matrix)
		vector <IndData> individual;           // Stores data on individuals
			
		vector <Enter> enter;                  // Stores probability distribution for entering
		
		vector <DataSource> source;            // Sources of data
		InitCond init_cond;                    // Stores infomation about the initial conditions
		
		vector <NMTrans> nm_trans;             // A list of non-Markovian transitions
		
		vector <NMTransIncomp> nm_trans_incomp;// A list of incomplete non-Markovian transitions
		
		vector <unsigned int> obs_eqn;         // Lists observation equations
		
		// References observation transitions [trgl][ti][#]
		vector < vector < vector <unsigned int> > > obs_trans_eqn_ref;
		
		unsigned int last_obs_trans_ti;        // Last ti point for obs_trans
		
		vector <unsigned int> obs_trans_eqn;   // Lists all equations used in obs_trans 

		vector <ObsTrans> obs_trans;           // Any observations on transitions
		
		vector < vector <bool> > comp_period;  // Determines if comp begins or ends transition period [cgl][cl]
		
		vector < vector <bool> > comp_terminal_cl;  // Determines if compartment is terminal [cl][cgl]
		
		vector <bool> sink_exist;              // Determines if a sink exists [cl]
		
		vector <bool> comp_terminal;           // Determines if compartment is terminal [cgl]
		
		bool obs_trans_exist;                  // Set if there are observed transitons
		
		bool period_exist;                     // Set if a period transition exits
		
		bool add_rem_pop_on;                   // Set if there is any addition or removal
		vector < vector <double> > add_rem_pop;// Add / remove population [ti][cgl]
		vector < vector <unsigned int> > add_rem_pop_change; // Store when add_rem_pop is non-zero
		
		vector <PopFilter> pop_filter;         // Filters used for population data
		
		vector <PopData> pop_data;             // Stores population data
		
		vector < vector < vector <unsigned int> > > pop_data_ref; // References population transition data [ti][cgl]
		
		vector <PopTransFilter> pop_trans_filter; // Filters used for pop trans data  
		
		vector <PopTransData> pop_trans_data;  // Stores population transition data
		bool pop_trans_data_exist;
		
		vector < vector < vector <unsigned int> > > pop_trans_ref; // References population transition [ti][tr]
		
		vector <ParEventJointProp> par_event_joint;// Information about joint parameter event proposal
		
		vector <WarnData> warn;                // Stores any errors used to generate data

		unsigned int T;                        // The maximum time
		
		unsigned int line_num;                 // Stores the import line (for diagnostic error messages)
		
		// Determines if global transitions are connected [tr_from][tr_to]
		//vector < vector <bool> > tr_connected; 
		
		// List all potential transitions after tr [tr][#]
		vector < vector <unsigned int> > tr_after;
		
		// List all potential transitions before tr [tr][#]
		vector < vector <unsigned int> > tr_before;
		
		Hash hash_ind;                         // Stores individuals in a hash table
		
		Hash hash_enter;                       // Stores enter in a hash table

		void create_markov_tree();
		unsigned int find_c(unsigned int cl, string name) const;
		unsigned int update_c_comp(unsigned int c, unsigned int cl, unsigned int c_comp) const;
		void ind_ev_sampler_not_needed();
		bool correct_ev(unsigned int c, Event &enew) const;
		void create_ind_noobs();
		vector < vector <double> > get_focal_init_cond(const InitCondValue &icv) const;
		void set_ind_samp_needed(const vector <Equation> &eqn);
		double calc_enter_prob(unsigned int c, unsigned int entref, const vector <double> &obs_eqn_value) const;
		vector <unsigned int> get_cv_list(unsigned int ie, const vector <Param> &param, const vector <Equation> &eqn) const;
		vector < vector <double> > calculate_omega_basic(unsigned int g, const vector <double> &param_val, const vector <Param> &param) const;
		unsigned int tr_trans(unsigned int trg, unsigned int c) const;
		unsigned int get_comp_global_convert(unsigned int cgl, unsigned int cl, unsigned int c) const;
		void sampling_error(unsigned int trg, string warn) const;
		
	// In 'species_data.cc':
	public: 
		void initialise_data();
		void order_data_events();
		Filter set_comp_filt(string te, unsigned int cl_not_allow, BoundType bound, const DataSource &so);
		void init_pop_trans_ref();
		void init_pop_data_ref();
		vector <unsigned int> obs_eqn_add_vec(vector <EquationInfo> &eqn_info);
		EquationInfo he(EquationInfo eqn_inf, const DataSource &so);
		vector <AddRemLocal> find_add_rem_list() const;
		void add_tr_list(const vector <unsigned int> &tr_list, const vector <bool> &cp, vector <AddRemLocal> &add_rem_local) const;
		void add_unobs_Amatrix_ind();
		void jiggle_data();
		bool is_pop_num_fixed() const;
		void nm_trans_incomp_info();
		unsigned int get_tra_m(unsigned int tr, const Event &ev_orig) const;
		string tr_swap_name(unsigned int cl, const vector <TrSwap> &tswa) const;
		vector <unsigned int> get_vec_tr_swap(unsigned int c, const vector <TrSwap> &start) const;
		vector <unsigned int> get_vec_tr_swap_mid(unsigned int st, unsigned int num, const vector <EventCl> &timeline) const;
		vector < vector <double> > calc_nm_rate(bool calc_bp, const vector <double> &param_val, const vector <SplineValue> &spline_val, const vector < vector <double> > &popnum_t, const vector <Equation> &eqn, vector < vector <double> > &bp_store) const;
		void set_ob_trans_ev(const vector <Equation> &eqn);

 	private:
		void set_ind_tmin_tmax();
		void init_pop_data(const DataSource &so);
		vector < vector <double> > pop_get_percentage(const DataSource &so);
		void add_pop_data(const DataSource &so, int sign);
		void set_add_rem_pop_change();
		void add_ind_data(const DataSource &so);
		void remove_ind_data(const DataSource &so);
		void move_ind_data(const DataSource &so);
		void genetic_data(const DataSource &so);
		void comp_data(const DataSource &so);
		void test_data(const DataSource &so);
		void population_data(const DataSource &so);
		ObsModelVariety set_obs_mod_type(const ObsModel &om) const;
		double set_obs_mod_val(double value, unsigned int j, unsigned int col, const ObsModel &om, const Table &tab, const DataSource &so);
		void trans_data(const DataSource &so);
		void popu_trans_data(const DataSource &so);
		vector <string> global_convert(const Filter &filt) const;
		vector <string> trans_global_convert(unsigned int cl, const vector <string> &trans_name_str, const Filter &comp_filt);
		void set_default_enter();
		unsigned int find_individual(string name, bool create=true);
		void X_vector_order();
		void generate_data() const;
		void alert_source(string st, const DataSource &so, unsigned int c = UNSET, unsigned int r = UNSET);
		vector <string> find_comp_prob_str(unsigned int cl, string te, BoundType bound, string &emsg) const;
		vector <string> set_trans_filt(unsigned int cl, string te, BoundType bound, string &emsg) const;
		string add_bound(string te, BoundType bound) const;
		vector <EquationInfo> create_eqn_vector(const vector <string> &vec, EqnType type, const DataSource &so);
		unsigned int find_cl(string name) const; 
		unsigned int find_tr(unsigned int cl, string name) const;
		
	public:
	const vector <double> &timepoint;
	const Details &details;
};


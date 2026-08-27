#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "equation.hh"
#include "precalc.hh"
#include "species.hh"
#include "hash.hh"

class Model                                // Stores information about the model
{
	public:
		SampType samp_type;                    // Used to change individual sampling
	
		bool no_question;                      // If questions are asked (data-sim) 
	
		bool sync_on;                          // Determines if sychronisation of proposals is done
		
		vector <Species> species;              // Compartmental model /data for different species
		unsigned int nspecies;                 // The number of species
		
		vector <SpeciesSimp> species_simp;     // A simplifies version passed to equation
		
		bool deterministic;                    // Set if the model is deterministic
		
		vector <CompPos> comp_pos;             // Stores compartment possibilities (to speed up 
		Hash hash_comp_pos;
		
		vector <Prior> prior;                  // Vector of priors
		
		vector <Param> param;                  // Stores all the parameters in the model
			
		vector <ParamVecEle> param_vec;        // A vector of parameter references 
		// (this vector is used to store values and is ordered by dependency) 
		
		unsigned int nparam_vec;               // The total number of parameter
		
		vector <unsigned int> param_vec_latin; // Stores which parameter are latin hypercude
		
		vector <unsigned int> param_vec_prop;  // Vector of param_vec elements which undergo proposals
		unsigned int nparam_vec_prop;          // Total number of parameters under proposal
			
		vector <unsigned int> param_vec_tvreparam;// Vector of tvreparam (needed for op of param w/o state)
		unsigned int nparam_vec_tvreparam;     // Total number of tvreparam
		
		bool contains_tvreparam;               // Set if model contains some tvreparam
		
		vector <unsigned int> dynamic_param;   // Stores list of dynamic parameters (used for output)
		
		vector < vector < vector <unsigned int> > > pop_reparam_th; // Sets reference from pop to th [po][ti][#]
		bool pop_reparam_th_on;                // Determines if on
		
		vector <Population> pop;               // Stores information about populations of interest
		unsigned int npop;                     // Number of pop (excluding those used for derive)
		
		vector <PopComb> popcomb;              // Population combination
		unsigned int npopcomb;                 // Number of popcomb (excluding those used for derive)
		
		vector <PopCombWeight> popcombw;       // Stores weight used in popcomb
		unsigned int npopcombw;                // Number of popcomb weights (excluding those used for derive)
	
		vector <unsigned int> popcombw_not_const; // Store all non constant popcombw
	
		//vector <PopcombUpdate> popcomb_update_sample; // Determines how popcomb should be updated
	
		Hash hash_pop;                         // Hash table for populations
	
		GeneticData genetic_data;              // Stores genetic data
		 
		// This is a quick reference to see the population which causes infection
		vector < vector < vector < vector <InfCause> > > > inf_cause; // [p][tr_gl][p_from][c_from]
		
		vector <IEGref> ieg_ref;               // References all the ind effect groups in the model
		
		vector <Spline> spline;                // Stores information about splines
	 	Hash hash_spline;
		
		vector <double> timepoint;             // Discrete timepoints	
		
		vector <Derive> derive;                // Derived quantities

		vector <Define> define;                // Used to to definitions
		
		vector <Density> density;              // Used for the DEN and RDEN functions 
	
		vector <Equation> eqn;                 // Stores all the equations for the species
	
		vector <EqnCheck> eqn_check;           // Used to check that equations are correct
		
		vector < vector <string> > region;     // Stores any regions (used for dynamic-sim)
		
		bool trans_tree;                       // Set if trans_tree exists in the model
		bool trans_tree_output;                // Set if trans_tree is only applied at output (because no genetic data)
		
		vector <Table> param_samp_store;       // Stores parameter samples (used for EXT);
		
		Table inf_param_stats;                 // Stores parameter stats (used in tornado)
		
		Table inf_pred_acc;                    // Stores prediction accuracy (used in tornado)
		
		vector < vector <string> > ind_key_store; // Converts from number to name
		vector <Sample> sample;                // Stores inferred state (used for PPC/EXT/DATA_SIM)
		
		vector <TerminalInfo> terminal_info;   // Loaded proposal information [chain]
		
		Operation mode;                        // SIM for simulation, INF for inference, PPC for post sim
		
		ExtFactor ext_factor;                  // Extension factor                  
		
		Details details;                       // Details for simulation/inference
		
		Details inf_details;                   // Details for inference (used in tornado plot 
		
		Constant constant;                     // Stores all the contants in the model
		
		Precalc precalc_eqn;                   // The equation which precalculates quantities
		SpecPrecalc spec_precalc;              // Stores precalcultion
		SpecPrecalc spec_precalc_derive;       // Stores list for precalcultion of derived (for integrals)
		SpecPrecalc spec_precalc_sample;       // Collects all the precalc to calculate after sample        
		vector <unsigned int> spec_precalc_time_ref;  // References spec_precalc_list
		vector <SpecPrecalcTime> spec_precalc_list;   // Stores a time-dependent spec_precalc
		vector <double> precalc_init;          // Initial value for precalc
		
		vector <unsigned int> param_vec_refq;   // Stores where param are on precalc (using q)
		vector <unsigned int> spline_refq;      // Stores where spline is on precalc (using q)
	
		vector <unsigned long long> param_vec_refi;   // Stores where param are on precalc (using i)
		vector <unsigned long long> spline_refi;      // Stores where spline is on precalc (using i)
	
		Model(Operation mode_, ExtFactor ext_factor_, bool no_question_);
		void add_eq_ref(EquationInfo &eqi, Hash &hash_eqn, double tdiv = UNSET, bool keep_te = false);
		void param_val_init(PV &param_val) const;
		PV param_sample(bool no_precalc = false) const;
		vector < vector <double> > param_latin_hypercube_prior_sample(unsigned int N) const;
		PV param_latin_hypercube_sample(const vector <double> &prior_val) const;
		void param_spec_precalc_time(unsigned int ti, const vector < vector <double> > &popcomb_t, PV &param_val, bool store) const;
		//void dynamic_param_calc(unsigned int ti, const vector < vector <double> > &popcomb_t, PV &param_val, double &time, double &time2) const;
		bool sample_bounded() const;
		void sample_ieg_cv(PV &param_val) const;
		void param_spec_precalc_time_all(const vector < vector <double> > &popcomb_t, PV &param_val, bool store) const;
		//void param_spec_precalc_before(unsigned int th, PV &param_val, bool store) const;
		//void param_spec_precalc_after(unsigned int th, PV &param_val, bool store) const;
		PV post_param(const Sample &samp) const;
		PV post_param_av(const vector <double> &value_av) const;
		IEstore post_ie_store(const Sample &samp) const;
		vector <double> prior_prob(const PV &param_val) const;
		vector <double> dist_prob(const PV &param_val) const;
		double recalculate_ieg_prior(unsigned int j, vector <double> &prior_ieg, const PV &param_val, double &like_ch) const;
		double recalculate_prior(unsigned int th, vector <double> &prior_prob, const PV &param_val, double &prior_ch, double &prior_bounded_ch) const;
		double recalculate_dist(unsigned int th, vector <double> &prior_prob, const PV &param_val, double &like_ch) const;
		vector <double> prior_ieg_all(const PV &param_val) const;
		double prior_ieg_calculate(const IEGref &iegr, const PV &param_val) const;
		void set_ieg_ref();
		void set_omega_pv();
		string exchange_omega(string name) const;
		vector <double> spline_prior(const PV &param_val) const;
		double recalculate_spline_prior(unsigned int s, vector <double> &spline_prior, PV &param_val, double &like_ch) const;
		void create_species_simp();
		void setup_distgrid();
		double grid_geo_dist(double xi, double yi, double xi2, double yi2, double xmin, double xmax, double ymin, double ymax);
		void order_affect(vector <AffectLike> &vec) const;
		//AffectMap get_affect_map(vector <AffectLike> &vec, const vector <unsigned int> &param_list, const vector <unsigned int> &dependent, const SpecPrecalc &spec_precalc_after) const;
		void add_iif_w_affect(vector <AffectLike> &vec) const;
		void add_popnum_ind_w_affect(vector <AffectLike> &vec) const;
		void joint_affect_like(PropType type, const vector <bool> &tr_change, unsigned int p, vector <AffectLike> &vec) const;
		void add_ie_affect(unsigned int p, unsigned int ie, vector <AffectLike> &vec) const;
		void add_pop_affect(unsigned int po, vector <AffectLike> &vec) const;
		bool inbounds(const PV &param_val) const;
		double prior_total(const PV &param_val) const;
		double dist_total(const PV &param_val) const;
		vector <InitCondValue> initc_sample(const PV &param_val) const;
		vector <double> frac_combine(unsigned int p, const vector < vector <double> > &frac) const;
		void combine_cnum_reduce(unsigned p, InitCondValue &inc) const;
		void add_like_obs_affect(unsigned int p, vector <AffectLike> &affect_like) const;
		AllInd find_all_ind(string name) const;
		void set_hash_all_ind();
		//bool ie_cholesky_error(const PV &param_val) const;
		bool ieg_check_prior_error(const IEGref &iegr, const PV &param_val) const;
		void print_param(const PV &param_val) const;
		vector <double> get_param_val_prop(const PV &param_val) const;
		vector <double> get_param_val_tvreparam(const PV &param_val) const;
		vector < vector <double> > get_param_val_dynamic(const PV &param_val) const;
		void add_tvreparam(PV &param_val, const vector <double> &param_val_tvreparam) const;
		vector <double> compress_vec(const vector <double> &vec) const;
		vector <double> decompress_vec(const vector <double> &vc) const;
		PV get_param_val(const Particle &pa) const;
		PV set_param_val(const vector <double> &value_set) const;
		double calc_tdiv(double t) const; 
		double calc_t(double tdiv) const;
		void create_precalc_equation();
		vector <SpecPrecalc> get_dynamic_spec_list(const SpecPrecalcTime &spt, vector <bool> &mapQ) const;
		SpecPrecalcTime& get_spec_precalc_time(unsigned int ti);
		void set_precalc_init();
		void set_popcombw_not_const();
		void set_spec_precalc_sample();
		//void set_spec_precalc_all();
		void set_spec_precalc_time();
		void create_precalc_derive();
		unsigned get_ti_ref(const vector <unsigned int> &list, vector < vector <unsigned int> > &list_time, Hash &hash_list_time) const;
		void precalc_affectQ();
		//vector <PopcombUpdate> get_popcomb_update(const vector <unsigned int> &list) const;
		void add_dynamic_spline(const vector < vector <AffectQ> > &affectQ);
		void add_affect_likeQ(unsigned int q, vector < vector <unsigned int> > &ti_ref_list, vector <bool> &map_pop, vector <bool> &map_popcombw, vector <bool> &map_popcomb, ParamVecEle &pvec, const vector <unsigned int> &affect_popcombw, const vector < vector < vector <unsigned int> > > &me_precalc, const vector < vector <AffectQ> > &affectQ, const vector < vector <unsigned int> > &list_time);
		
		void print_precalc() const;
		vector <unsigned int> get_last_spline() const;
		
		bool in_bounds(double x, unsigned int j, const vector <double> &precalc) const;
		bool is_prior_bounded(unsigned int th) const;
		bool is_matrix(const Param &par) const;
		bool is_symmetric(const Param &par) const;
		void convert_fix_pr_const();
		void print_spec_precalc(string st, const SpecPrecalc &spre) const;
		string str_time_range(const vector <unsigned int> &lt) const;
		string str_spec_precalc(string st, const SpecPrecalc &spre) const;
		void set_param_spec_precalc();
		void set_pop_reparam_th();
		string load_prop_info(unsigned int ch, const vector <string> &lines);
		unsigned int get_end_bracket(string &te, unsigned int i);
		bool in_integral(unsigned int i, const vector <SumRange> &int_range) const;
		vector <string> equation_dep(string te, string &warn);
		bool time_specified_pop(unsigned int i, const string &te) const;
		unsigned int find_p(string name) const;
		unsigned int find_cl(unsigned int p, string name) const;
		bool data_mode() const;
		bool no_process() const;
		bool question(string st) const;
		unsigned int get_cl_from_trans(string name, unsigned int p) const;
		unsigned int get_cl_from_comp(string name, unsigned int p) const;
		double calculate_equation_zero_one(string te, double tdiv, string &err);
		double calculate_equation(string te, double tdiv, string &err);
		DiagTestSens get_diag_test_sens(string comp, unsigned int p, string &warn) const;
		void extract_popcomb(Hash &hashw, Hash &hashpc);
		void extract_popcomb_derive(Hash &hashw, Hash &hashpc);
		string print_it(const EqItem &it) const;
		void print_popcomb() const;
		void check_all_linear() const;
		void set_pop_list();
		void set_pop_grad_ref();
		void precalc_calculate_all(PV &param_val, const vector < vector <double> > &popcomb_t) const;
		void create_equation_check();
		void check_equation_values(const vector < vector <double> > &popnum_t, const vector < vector <double> > &popcomb_t, const PV &param_val, const vector < vector < vector <double> > > &derive_val) const;
		string param_vec_name(unsigned int th) const;
		string spline_name(unsigned int i) const;
		void set_eqn_precalcnum();
		bool par_ele_set(unsigned int th, unsigned int j) const;
		bool is_bernoulli(unsigned int th) const;
		bool is_strictly_positive(const Prior &pri) const;
	
	private:
		Hash hash_all_ind;                     // Stores individuals in a hash table
		
		vector <AllInd> all_ind;               // Stores reference for all individuals
	
		double prior_sample(const Prior &pri, const vector <double> &precalc) const;
		void prior_error(const Prior &pri, string warn) const;
		vector <double> latin_prior_sample(const Prior &pri, unsigned int N) const;
		long get_hash_number(const EquationInfo &eqi);	
};


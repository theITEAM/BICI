#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "equation.hh"
#include "precalc.hh"
#include "species.hh"
//#include "state_species.hh"
#include "hash.hh"

class Model                                // Stores information about the model
{
	public:
		SampType samp_type;                    // Used to change individual sampling
	
		vector <Species> species;              // Compartmental model /data for different species
		unsigned int nspecies;                 // The number of species
		
		vector <SpeciesSimp> species_simp;     // A simplifies version passed to equation
		
		vector <Prior> prior;                  // Vector of priors
		
		vector <Param> param;                  // Stores all the parameters in the model
			
		vector <ParamVecEle> param_vec;        // A vector of parameter references 
		// (this vector is used to store values and is ordered by dependency) 
		
		unsigned int nparam_vec;               // The total number of parameter
		
		vector <unsigned int> param_vec_prop;  // Vector of param_vec elements which undergo proposals
		unsigned int nparam_vec_prop;          // Total number of parameters under proposal
			
		vector <unsigned int> param_vec_tvreparam;// Vector of tvreparam (needed for op of param w/o state)
		unsigned int nparam_vec_tvreparam;     // Total number of tvreparam
		
		bool contains_tvreparam;               // Set if model contains some tvreparam
		
		vector < vector < vector <unsigned int> > > pop_reparam_th; // Sets reference from pop to th [po][ti][#]
		
		vector <Population> pop;               // Stores information about populations of interest
		
		Hash hash_pop;                         // Hash table for populations
	
		GeneticData genetic_data;              // Stores genetic data
		 
		// This is a quick reference to see the population which causes infection
		vector < vector < vector < vector <InfCause> > > > inf_cause; // [p][tr_gl][p_from][c_from]
		
		vector <IEGref> ieg_ref;               // References all the ind effect groups in the model
		
		vector <Spline> spline;                // Stores information about splines
	 	Hash hash_spline;
		
		vector <double> timepoint;             // Discrete timepoints	
		
		vector <Derive> derive;                // Derived quantities
		  
		vector <Density> density;              // Used for the DEN and RDEN functions 
	
		vector <Equation> eqn;                 // Stores all the equations for the species
	
		bool trans_tree;                       // Set if trans_tree exists in the model
		
		vector <Sample> sample;                // Stores inferred state (used for PPC)
		
		Operation mode;                        // SIM for simulation, INF for inference, PPC for post sim
		
		Details details;                       // Details for simulation/inference
		
		Constant constant;                     // Stores all the contants in the model
		
		Precalc precalc_eqn;                   // The equation which precalculates quantities
		SpecPrecalc spec_precalc;              // Stores precalcultion
		SpecPrecalc spec_precalc_derive;       // Stores list for precalcultion of derived (for integrals)
		SpecPrecalc spec_precalc_sample;       // Collects all the precalc to calculate after sample        
		SpecPrecalc spec_precalc_all;          // Collects all the precalc         
		vector <SpecPrecalcTime> spec_precalc_time; // Precalculation at different times
		vector <double> precalc_init;          // Initial value for precalc
		
		vector <unsigned int> param_vec_ref;   // Stores where param are on precalc
		vector <unsigned int> spline_ref;      // Stores where spline is on precalc
	
		Model(Operation mode_);
		void add_eq_ref(EquationInfo &eqi, Hash &hash_eqn, double tdiv = UNSET);
		void param_val_init(PV &param_val) const;
		PV param_sample() const;
		void param_spec_precalc_time(unsigned int ti, const vector < vector <double> > &popnum_t, PV &param_val, bool store) const;
		bool sample_bounded() const;
		void sample_ieg_cv(PV &param_val) const;
		void param_spec_precalc_time_all(const vector < vector <double> > &popnum_t, PV &param_val, bool store) const;
		//void param_spec_precalc_before(unsigned int th, PV &param_val, bool store) const;
		//void param_spec_precalc_after(unsigned int th, PV &param_val, bool store) const;
		PV post_param(const Sample &samp) const;
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
		void order_affect(vector <AffectLike> &vec) const;
		void affect_linearise_speedup(vector <AffectLike> &vec) const;
		AffectMap get_affect_map(vector <AffectLike> &vec, const vector <unsigned int> &param_list, const vector <unsigned int> &dependent, const SpecPrecalc &spec_precalc_after) const;
		void affect_nopop_speedup(vector <AffectLike> &vec, const vector <unsigned int> &param_list, const vector <unsigned int> &dependent, const SpecPrecalc &spec_precalc_after) const;
		void set_factor_nopop_only(vector <AffectLike> &vec, const vector <unsigned int> &param_list, const vector <unsigned int> &dependent, const SpecPrecalc &spec_precalc_after) const;
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
		void add_tvreparam(PV &param_val, const vector <double> &param_val_tvreparam) const;
		PV get_param_val(const Particle &pa) const;
		double calc_tdiv(double t) const; 
		double calc_t(double tdiv) const;
		void create_precalc_equation();
		void set_precalc_init();
		void set_spec_precalc_sample();
		void set_spec_precalc_all();
		void set_spec_precalc_time();
		void create_precalc_derive();
		void create_precalc_pop_grad();
		void precalc_affect();
		void print_precalc() const;
		vector <unsigned int> get_last_spline() const;
		void add_affect_like(unsigned int i, unsigned int i2, const vector <bool> &map_time, vector <bool> &map_PC, vector < vector <bool> > &map_me, ParamVecEle &pvec, const vector < vector <unsigned int> > &affect, const vector < vector <AffectME> > &affect_me);
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

	private:
		Hash hash_all_ind;                     // Stores individuals in a hash table
		
		vector <AllInd> all_ind;               // Stores reference for all individuals
	
		double prior_sample(const Prior &pri, const vector <double> &precalc) const;
		long get_hash_number(const EquationInfo &eqi);	
};


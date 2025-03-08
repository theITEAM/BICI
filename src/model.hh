#ifndef BICI__MODEL_HH
#define BICI__MODEL_HH

#include <vector>

using namespace std;

#include "struct.hh"
#include "equation.hh"
#include "state_species.hh"
#include "hash.hh"

class Model                                // Stores information about the model
{
	public:
		SampType samp_type;                    // Used to change individual sampling
	
		vector <Species> species;              // Compartmental model /data for different species
		unsigned int nspecies;                 // The number of species
		
		vector <SpeciesSimp> species_simp;     // A simplifies version passed to equation
			
		vector <Param> param;                  // Stores all the parameters in the model
			
		vector <ParamVecEle> param_vec;        // A vector of parameter references 
		// (this vector is used to store values and is ordered by dependency) 
		
		unsigned int nparam_vec;               // The total number of parameter
		
		vector <Population> pop;               // Stores information about populations of interest
		
		Hash hash_pop;                         // Hash table for populations
	
		GeneticData genetic_data;              // Stores genetic data
		 
		// This is a quick reference to see the population which causes infection
		vector < vector < vector < vector <InfCause> > > > inf_cause; // [p][tr_gl][p_from][c_from]
		
		vector <Spline> spline;                // Stores information about splines
	 	
		vector <double> timepoint;             // Discrete timepoints	
		unsigned int ntimepoint;               // The number of time points
		
		vector <Derive> derive;                // Derived quantities
		  
		vector <Equation> eqn;                 // Stores all the equations for the species
	
		bool trans_tree;                       // Set if trans_tree exists in the model
		
		vector <Sample> sample;                // Stores inferred state (used for PPC)
		
		Operation mode;                        // SIM for simulation, INF for inference, PPC for post sim
		
		Details details;                       // Details for simulation/inference
		
		Model(Operation mode_);
		void add_eq_ref(EquationInfo &eqi, double t = UNSET);
		string find_eq_dif(const EquationInfo &eqi) const;
		vector <double> param_sample() const;
		vector <double> post_param(const Sample &samp) const;
		vector <double> prior_prob(const vector <double> &param_val) const;
		vector <double> dist_prob(const vector <double> &param_val) const;
		double recalculate_prior(unsigned int th, vector <double> &prior_prob, const vector <double> &param_val, double &like_ch) const;
		double recalculate_dist(unsigned int th, vector <double> &dist_prob, const vector <double> &param_val, double &like_ch) const;
		vector <double> spline_prior(const vector <double> &param_val) const;
		double recalculate_spline_prior(unsigned int s, vector <double> &spline_prior, const vector <double> &param_val, double &like_ch) const;
		vector <double> calculate_popnum(vector <StateSpecies> &state_species) const;
		vector < vector <double> > calculate_popnum_t(vector <StateSpecies> &state_species, unsigned int ti_end = UNSET) const;
		vector <double> recalculate_population(vector < vector <double> > &popnum_t, const vector <unsigned int> &list, const vector <StateSpecies> &state_species) const;
		void recalculate_population_restore(vector < vector <double> > &popnum_t, const vector <unsigned int> &list, const vector <double> &vec) const;
		void create_species_simp();
		void order_affect(vector <AffectLike> &vec) const;
		bool div_value_fast_possible(const AffectLike &al) const;
		void affect_linearise_speedup(vector <AffectLike> &vec) const;
		void affect_linearise_speedup2(vector <AffectLike> &vec, const vector <unsigned int> &param_list, const vector <unsigned int> &dependent) const;
		void add_iif_w_affect(vector <AffectLike> &vec) const;
		void add_popnum_ind_w_affect(vector <AffectLike> &vec) const;
		void joint_affect_like(PropType type, const vector <bool> &tr_change, unsigned int p, vector <AffectLike> &vec) const;
		void add_ie_affect(unsigned int p, unsigned int ie, vector <AffectLike> &vec) const;
		void add_pop_affect(unsigned int po, vector <AffectLike> &vec) const;
		bool inbounds(const vector <double> &param_val) const;
		double prior_total(const vector <double> &param_val) const;
		double dist_total(const vector <double> &param_val) const;
		vector <InitCondValue> initc_sample(const vector <double> &param_val) const;
		vector <double> frac_combine(unsigned int p, const vector < vector <double> > &frac) const;
		void combine_cnum_reduce(unsigned p, InitCondValue &inc) const;
		void add_like_obs_affect(unsigned int p, vector <AffectLike> &affect_like) const;
		AllInd find_all_ind(string name) const;
		void set_hash_all_ind();
		
	private:
		Hash hash_eqn;                         // Stores a hash take for equations
		
		Hash hash_all_ind;                     // Stores individuals in a hash table
		
		vector <AllInd> all_ind;               // Stores reference for all individuals
	
		double prior_sample(const Prior &pri, const vector <double> &param_val) const;
		long get_hash_number(const EquationInfo &eqi);
};
#endif

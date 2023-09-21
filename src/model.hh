#ifndef BICI__MODEL_HH
#define BICI__MODEL_HH

#include <vector>

using namespace std;

#include "struct.hh"
#include "equation.hh"
#include "state_species.hh"

//#include "utils.hh"

class Model                                             // Stores information about the model
{
	public:
		vector <Species> species;                           // Compartmental model /data for different species
		unsigned int nspecies;                              // The number of species
		
		vector <SpeciesSimp> species_simp;                  // A simplifies version passed to equation
			
		vector <Param> param;                               // Stores all the parameters in the model
			
		vector <BranchParam> branch_param;                  // Groups together branching parameters
			
		vector <ParamVecEle> param_vec;                     // A vector of parameter references (this vector is used to store values and is ordered by dependency) 
		unsigned int nparam_vec;                            // The total number of parameter
		
		vector <Population> pop;                            // Stores information about populations of interest
	
		vector <Spline> spline;                             // Stores information about splines
	 	
		vector <double> timepoint;                          // Discrete timepoints	
		unsigned int ntimepoint;                            // The number of time points
		
		vector <Derive> derive;                             // Derived quantities
		  
		vector <Equation> eqn;                              // Stores all the equations for the species
	
		Operation mode;                                     // SIM for simulation and INF for inference	
		Details details;                                    // Details for simulation/inference
		
		Model();
		void add_eq_ref(EquationInfo &eqi);
		vector <double> param_sample() const;
		vector <double> prior_prob(const vector <double> &param_val) const;
		vector <double> dist_prob(const vector <double> &param_val) const;
		double recalculate_prior(unsigned int th, vector <double> &prior_prob, const vector <double> &param_val, double &like_ch) const;
		double recalculate_dist(unsigned int th, vector <double> &dist_prob, const vector <double> &param_val, double &like_ch) const;
		vector <double> spline_prior(const vector <double> &param_val) const;
		double recalculate_spline_prior(unsigned int s, vector <double> &spline_prior, const vector <double> &param_val, double &like_ch) const;
		
		vector <double> calculate_popnum(vector <StateSpecies> &state_species) const;
		vector < vector <double> > calculate_popnum_t(vector <StateSpecies> &state_species) const;
		vector <double> recalculate_population(vector < vector <double> > &popnum_t, const vector <unsigned int> &list, const vector <StateSpecies> &state_species) const;
		void recalculate_population_restore(vector < vector <double> > &popnum_t, const vector <unsigned int> &list, const vector <double> &vec) const;

		void create_species_simp();
		void set_branch_auto(vector <double> &param_val) const;
		void order_affect(vector <AffectLike> &vec) const;
		bool inbounds(const vector <double> &param_val) const;
		double prior_total(const vector <double> &param_val) const;
	
	private:
		double prior_sample(const Prior &pri, const vector <double> &param_val) const;
		double prior_probability(double x, const Prior &pri, const vector <double> &param_val) const;

	// In 'model_data.cc'
	public:
		vector <InitCondValue> initc_sample() const;
		
	
};
#endif

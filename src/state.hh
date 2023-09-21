#ifndef BICI__STATE_HH
#define BICI__STATE_HH

#include <vector>


using namespace std;

#include "struct.hh"
#include "model.hh"

class State                                             // Stores information about the state
{
	public:
		double t_start, t_end;                              // The time range
		double dt;                                          // The time step
		
		unsigned int T;                                     // Total number of divisions
		
		vector <double> dtimepoint;                         // The difference between timepoints
		
		vector <unsigned int> empty;
		
		vector <double> param_val;                          // The parameter value
		
		vector <SplineValue> spline_val;                    // Values for the spline
		
		Like like;                                          // Stores the overall likelihood and prior
		
		vector <double> prior_prob;                         // The prior probability for each parameter
		
		vector <double> dist_prob;                          // The distribution probability for each parameter
		
		vector <double> spline_prior;                       // The prior probability a spline
		
		vector < vector <double> > popnum_t;                // The time variation in the population
		
		vector <StateSpecies> species;                      // State information about species
		unsigned int nspecies;
		
		vector <Change> change;                             // Stores changes in the state
		
		long check_timer;                                   // CPU time for performing checks
		
		State(const Model &model);
		void simulate(const vector <double> &param_value, const vector <InitCondValue> &initc_val);
		vector <DeriveOutput> derive_calculate() const;
		void calculate_likelihood();
		void calculate_like();
		void add_like(Like like_ch);
		void change_add(const vector <double> &vec);
		void change_add(double num);
		Like update_param(const vector <AffectLike> &affect_like);
		void restore(const vector <AffectLike> &affect_like);
		bool check_invalid();
		void likelihood_from_scratch();
		void resample_ind();
		Particle generate_particle() const;
		void set_particle(const Particle &part);
		
	private:
		void spline_init();
		vector <double> recalculate_spline(unsigned int j, const vector <unsigned int> &list);
		void print_cpop(unsigned int ti) const;	
	
	// In state_update_ind.cc
	public:
		//void ind_proposal(unsigned int p,  unsigned int i, unsigned int cl);
		void update_markov_eqn_value(MarkovEqnDiv &div, double value, double &Li_markov, double &like_ch);
		Like update_ind(unsigned int p, unsigned int i, unsigned int cinit_new, vector <Event> &ev_new);
		void update_pop_change(unsigned int ti, unsigned int ti_next, double &like_ch);
		void update_ie_population(unsigned int p, unsigned int i, unsigned int ie, double factor, double &like_ch);
		double update_pop_trans(double t, unsigned int tr_gl, int sign);
		
	private:
		void initialise_update_ind_maps();
		vector <unsigned int> pop_list;                     // Used in state_update_ind
		vector <int> pop_map;
		vector < vector <bool> > markov_eqn_map; 
		vector <PopMarkovEqnRef> markov_eqn_list;
		vector < vector <bool> > trans_map; 
		vector <PopTransRef> trans_list;
		
		
	// In state_check.cc
	public:
		void check(string ref);

	private:
		void check_ref(unsigned int p, string ref);
		void check_markov_trans(unsigned int p, string ref);
		void check_nm_trans(unsigned int p, string ref);
		void check_prior(string ref);
		void check_popnum_t(string ref);
		void check_like(string ref);
		void check_spline(string ref);
		void check_ie(unsigned int p, string ref);
		void check_pop_like(unsigned int p, string ref);
		void check_cpop_st(string ref);
		void check_maps(string ref);
		void check_obs_like(unsigned int p, string ref);
			
	private:
		const Model &model;
};
#endif

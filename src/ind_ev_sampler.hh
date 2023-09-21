#ifndef BICI__IND_EVENT_SAMPLER_HH
#define BICI__IND_EVENT_SAMPLER_HH

#include <vector>

using namespace std;

#include "struct.hh"
#include "species.hh"
//#include "model.hh"

class IndEvSampler                       // Stores information about an individual event sampler
{
	public:
		IndEvSampler(vector <MarkovEqnVariation> &markov_eqn_vari, const vector <Individual> &individual, const Details &details, const Species &sp, vector <double> &obs_eqn_value, const vector <Equation> &eqn);
	
		bool illegal;                             // Determines if there is a problem with the sampler
		
		unsigned int i_store;                     // The individual undergoing the sampling
		unsigned int cl_store;                    // The classification
	
		unsigned int ti_start;                    // The time division the sampler starts
		unsigned int ti_end;                      // The time division the sampler ends
		
		unsigned int k_start;                     // The event being started at
	
		double t_start;                           // The start time 
		
		unsigned int c_start;                     // The starting compartment
		
		unsigned int nobs;                        // The number of observations
		
		unsigned int isl;                         // The island the ind class is on
	 
		unsigned int C;                           // Number of compartments in island
		
		unsigned int T;                           // The maximum time
		vector < vector <double> > ind_obs_prop;  // The observation probability for an individual
		vector <unsigned int> c_timeline;         // The timeline for g        
	
		bool calc_indfac;                         // Set to true if indfac needs to be calculated     
	
		vector <unsigned int> fixed_event_obs;    // Stores fixed observed event
	
		vector <Event> ev_new;
	
		bool needed(unsigned int i, unsigned int cl);
		void generate_ind_obs_timeline();
		unsigned int get_ti_lower(double t) const;
		unsigned int get_ti_upper(double t) const;
		vector <Event> sample_events(double &probif);
		double sample_events_prob(const vector <Event> &ev) const;
		double get_indfac(const Individual &ind, const MarkovEqn &mar_eqn) const;
		
		vector <MarkovEqnVariation> &markov_eqn_vari;
		const vector <Individual> &individual;
		const Details &details;
		const Species &sp;
		vector <double> &obs_eqn_value;
		const vector <Equation> &eqn;
};
#endif

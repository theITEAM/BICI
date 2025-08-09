#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "species.hh"

struct InitStateProb {                     // Used for sampling the initial state
	double prob;                             // The probability of compartment
	double prob_sum;                         // Cumulative probability
};


class IndEvSampler                         // Stores information about an individual event sampler
{
	public:
		IndEvSampler(vector <MarkovEqnVariation> &markov_eqn_vari, const vector <Individual> &individual, const Details &details, const Species &sp, vector <double> &obs_eqn_value, vector < vector <double> > &obs_trans_eqn_value, const vector <Equation> &eqn, const vector <InfNode> &inf_node, const vector <double> &param_val, const vector <SplineValue> &spline_val, vector < vector <double> > &popnum_t);
	
		bool illegal;                          // Determines if there is a problem with the sampler
		
		unsigned int i_store;                  // The individual undergoing the sampling
		unsigned int cl_store;                 // The classification
	
		unsigned int ti_start;                 // The time division the sampler starts
		unsigned int ti_end;                   // The time division the sampler ends
		
		double t_start;                        // The start time 
		
		double dt;                             // The timestep
		
		unsigned int c_start_raw;              // The start with the compartment in the selected classification set to zero
		
		unsigned int nobs;                     // The number of observations
		
		unsigned int isl;                      // The island the ind class is on
	 
		unsigned int C;                        // Number of compartments in island
		
		unsigned int T;                        // The maximum time
		vector < vector <double> > ind_obs_prob;// The observation probability for an individual
		vector <unsigned int> c_timeline;      // The timeline for global compartment number      
	
		vector < vector < vector <double> > > rate_store;// Stores expected number of events
		
		bool ind_variation;                    // Set if ind_effect or fix_effect 
		
		vector <unsigned int> fixed_event_obs; // Stores fixed observed event
	
		vector < vector <double> > nm_rate;    // Store the non-markovian rate [m][ti]
		
		vector <Event> ev_new;                 // New event sequence
	
		bool needed(unsigned int i, unsigned int cl);
		void rate_store_initialise();
		double get_trans_obs_prob(unsigned int trg, const ObsData &ob) const;
		void generate_ind_obs_timeline();
		vector <Event> sample_events(double &probif);
		double add_event(const Event &e, vector <Event> &ev_new, unsigned int nc_index, unsigned int fe_index, const vector <IslandComp> &island_comp, const vector <ObsData> &obs, const vector <TrigEventRef> &fixed_trans_ev) const;
		bool contain_period(const vector <Event> &event, unsigned int cl) const;
		
		vector <InitStateProb> init_state_sampler_obs(unsigned int i, unsigned int ti, const Event &ev) const;
		double sample_events_prob(const vector <Event> &ev) const;
		double get_indfac(const Individual &ind, const MarkovEqn &mar_eqn) const;
		void setup_nm();
		double calculate_rate(const Individual &ind, const IslandTrans &le, unsigned int ctime, int tii) const;
		double calculate_nm_rate(const Individual &ind, unsigned int m, int ti) const;
		
		vector <MarkovEqnVariation> &markov_eqn_vari;
		const vector <Individual> &individual;
		const Details &details;
		const Species &sp;
		vector <double> &obs_eqn_value;
		vector < vector <double> > &obs_trans_eqn_value;
		const vector <Equation> &eqn;
		const vector <InfNode> &inf_node;
		const vector <double> &param_val;
		const vector <SplineValue> &spline_val;
		vector < vector <double> > &popnum_t;
		const vector <NMTrans> &nm_trans;
		
	// In ind_ev_sampler_noobs.cc
	Event resample_init_event(unsigned int i, double &probif) const;
	double resample_init_event_prob(unsigned int i, const Event &e_init) const;
	vector <InitStateProb> source_sampler(unsigned int i, const Event &ev) const;
	unsigned int init_comp_sampler(unsigned int i, double &probif) const;
	double init_comp_sampler_prob(unsigned int c, unsigned int i) const;
	//vector <InitStateProb> init_state_sampler(unsigned int i, const Event &ev) const;
	
	// In ind_ev_simulate.cc
	vector <Event> simulate_events(unsigned int i, const Event &e_init, double &probif, const vector <TrigEventRef> &ste);
	double find_tnext(unsigned int m, const vector <TrigEventRef> &trig_event, unsigned int f, 	 const vector <FutureNMEvent> &future_nm) const;
	void add_future_nm_event(unsigned int c, unsigned int cl, double t, unsigned int i, vector <FutureNMEvent> &future_nm, unsigned int m, const vector <TrigEventRef> &trig_event);
	void sim_add_event(EventType type, unsigned int tr_gl, double t, unsigned int cl, unsigned int move_c, unsigned int c, bool observed, vector <Event> &ev_new) const; 
	unsigned int sim_add_data_event(const TrigEventRef &te, unsigned int &c, unsigned int i, vector <Event> &ev_new, double &probif);
	void get_trans_rate(unsigned int c, unsigned int ti, const Individual &ind, vector <double> &trans_rate, double &R) const;
	double simulate_events_prob(unsigned int i, const vector <Event> &ev, const vector <TrigEventRef> &trig_event) const;
	void print_trigger_event(string te, const vector <TrigEventRef> &trig_event) const;
	
	// In ind_ev_simulate_single.cc
	vector <Event> simulate_single_events(unsigned int i, unsigned int cl, double &probif, const vector <TrigEventRef> &trig_event);
	double simulate_single_events_prob(unsigned int i, unsigned int cl, const vector <Event> &ev, const vector <TrigEventRef> &trig_event) const;
	double find_tnext_single(unsigned int m, const vector <TrigEventRef> &trig_event, unsigned int f, const vector <FutureNMEvent> &future_nm, unsigned int e, const vector <Event> &ev) const;
	void get_trans_rate_single(unsigned int c, unsigned int cl, unsigned int ti, const Individual &ind, vector <double> &trans_rate, double &R) const;
};


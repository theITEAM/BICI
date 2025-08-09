#pragma once

#include <vector>

using namespace std;

#include "struct.hh"

class SourceSampler {
	public:
	double enter_frac;                       // Estimates fraction of individuals entering (as opposed to a source)	
	
	double update_num;                       // The number of updates to be done
	 
	private:
	double memory;                           // Stores how quick sampler "forgets"
	
	vector <double> enter_prob;              // Probability of entering in different compartments
	vector <double> enter_prob_sum;          // Sums probability of entering in different compartments
	
	vector < vector <double> > rate;         // Proportional to rate of entering [ti][#tr_gl_list]
	double rate_total;
	vector < vector <double> > rate_sum;     // Sums across index [ti][#tr_gl_list]
	
	vector <double> margin_sum;              // Sum of the marginal of the rate
	vector <unsigned int> me_list;           // Lists all source markov equations
	vector <unsigned int> me_map;            // Maps from tr_gl to list

	public:
		SourceSampler(const vector <MarkovEqn> &markov_eqn, const vector <TransGlobal> &tra_gl, const vector <CompGlobal> &comp_gl, const Details &details, const InitCond &init_cond);
		void update(unsigned int nind_obs, vector <Individual> &individual, const vector <MarkovEqnVariation> &markov_eqn_vari, bool contains_source);
		Event sample(double &probif) const;
		double sample_prob(const Event &ev) const;
		string print() const;
	
	private:
		const vector <MarkovEqn> &markov_eqn;
		const vector <TransGlobal> &tra_gl;
		const vector <CompGlobal> &comp_gl;
		const Details &details;
		const InitCond &init_cond;
};


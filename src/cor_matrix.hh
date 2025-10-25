#pragma once

#include <vector>
#include "model.hh"

using namespace std;

class CorMatrix 
{    
public:
	CorMatrix(const Model &model);
   
	unsigned int n;                  // Number of samples added to av and av2

	void init();
	void add_sample(vector <double> param_val_prop, unsigned int range);
	void check() const;
	vector < vector <double> > calculate_cor_matrix() const;
	void set_mvn_from_particle(vector <Particle> &particle);
	vector < vector <double> > find_covar(const vector <unsigned int> &param_list) const;
	unsigned int get_n();
	
private:
	unsigned int n_start;            // The start of samp

	unsigned int N;                  // The number of parameters
	
	vector <double> av;              // The sum of parameter values 
	vector < vector <double> > av2;  // The sum of param*param 
	vector < vector <double> > samp; // Stores historic parameter values
	
	const Model &model;
};
// Deals with correlation matrices

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
using namespace std;

#include "cor_matrix.hh"
#include "utils.hh"

/// Initialise the correlation matrix 			
CorMatrix::CorMatrix(const Model &model) : model(model)
{
	N = model.nparam_vec_prop;
}

			
/// Initialises the correlation matrix
void CorMatrix::init()
{
	av.clear(); av2.clear();
	n = 0;
	n_start = 0;

	av.resize(N,0);
	av2.resize(N); for(auto i = 0u; i < N; i++) av2[i].resize(N,0);
	samp.clear();
}


/// Adds a sample
void CorMatrix::add_sample(const vector <double> &param_value, unsigned int range)
{
	if(range < RANGE_MIN) range = RANGE_MIN;
	
	vector <double> paramv(N);
	for(auto i = 0u; i < N; i++){
		paramv[i] = param_value[model.param_vec_prop[i]];
	}
	
	// Shifting ensure that values are not constant
	auto sh = (ran()-0.5)*TINY;
	for(auto &val : paramv) val += sh;
	
	samp.push_back(paramv);
	n++;

	for(auto i = 0u; i < N; i++){
		auto vali = paramv[i];
		av[i] += vali;
		for(auto j = i; j < N; j++){
			auto valj = paramv[j];
			av2[i][j] += vali*valj;
		}
	}

	while(n-n_start > range){
		const auto &samp_val = samp[n_start];
		for(auto i = 0u; i < N; i++){
			auto vali = samp_val[i];
			av[i] -= vali;
			for(auto j = i; j < N; j++){
				auto valj = samp_val[j];
				av2[i][j] -= vali*valj;
			}
		}
		n_start++;
	}
}


/// Checks that correlation matrix is up-to-date
void CorMatrix::check() const
{
	vector <double> ave(N,0);
	vector < vector <double> > ave2;
	ave2.resize(N);
	for(auto j = 0u; j < N; j++) ave2[j].resize(N,0);
	
	for(auto i = n_start; i < n; i++){
		for(auto j = 0u; j < N; j++){
			ave[j] += samp[i][j];
			for(auto k = j; k < N; k++){
				ave2[j][k] += samp[i][j]*samp[i][k];
			}
		}
	}
	
	for(auto j = 0u; j < N; j++){
		if(dif(ave[j],av[j],DIF_THRESH)) emsg("Different in av");
		for(auto k = j; k < N; k++){
			if(dif(ave2[j][k],av2[j][k],DIF_THRESH)){
				emsg("Different in av2");
			}
		}	
	}
}


/// Calculates the correlation matrix
vector < vector <double> > CorMatrix::calculate_cor_matrix() const
{
	auto nn = n-n_start;
	
	vector <double> mu(N);
	for(auto i = 0u; i < N; i++) mu[i] = av[i]/nn;
	
	vector <double> var(N);
	for(auto i = 0u; i < N; i++) var[i] = av2[i][i]/nn - mu[i]*mu[i];
	
	vector < vector <double> > M;
	M.resize(N); for(auto i = 0u; i < N; i++) M[i].resize(N);
	
	for(auto i = 0u; i < N; i++){
		for(auto j = i; j < N; j++){
			M[i][j] = (av2[i][j]/nn - mu[i]*mu[j])/sqrt(var[i]*var[j]);
			M[j][i] = M[i][j]; 
		}
	}
	
	return M;
}


/// Sets up av and av2
void CorMatrix::set_mvn_from_particle(vector <Particle> &particle)
{
	init();
	for(const auto &part : particle){
		add_sample(part.param_val,LARGE);
	}
}


/// Finds covariance matrix for a subset of parameters
vector < vector <double> > CorMatrix::find_covar(const vector <unsigned int> &param_list) const
{
	vector <unsigned int> param_list_prop;
	for(auto th : param_list){
		auto ref = model.param_vec[th].ref;
		if(ref == UNSET) emsg("Param vec ref unset");
		param_list_prop.push_back(ref);
	}
	
	auto nsum = n-n_start;
	auto np = param_list_prop.size();
	
	vector < vector <double> > M;
	M.resize(np); for(auto j = 0u; j < np; j++) M[j].resize(np);

	for(auto j = 0u; j < np; j++){
		auto th_j = param_list_prop[j];
		double val;
		for(auto i = j; i < np; i++){
			auto th_i = param_list_prop[i];
			if(th_j < th_i) val = av2[th_j][th_i]/nsum - (av[th_j]/nsum)*(av[th_i]/nsum);
			else val = av2[th_i][th_j]/nsum - (av[th_j]/nsum)*(av[th_i]/nsum);
	
			M[j][i] = val;
			M[i][j] = val;
		}
	}
	
	return M;
}


/// Returns the number of samples gathered
unsigned int CorMatrix::get_n()
{
	return n-n_start;
}

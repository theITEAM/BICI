// Functions which (approximately) integrate the likelihood 
// This done by allowing a probabilistic allocation of individuals to compartments

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
using namespace std;

#include "state.hh"
#include "utils.hh"

/// Intergates the approximate likelihood 
double State::like_integrate(const vector <double> &param_value)
{
	param_val = param_value;
	spline_init();
	
	for(auto p = 0u; p < model.species.size(); p++){
		const auto &sp = model.species[p];
		
		auto N = sp.individual.size();
		cout << N << endl;
		for(auto i = 0u; i < N; i++){
			
		}
		
	}
	emsg("do");
	//for(auto &indd
	
	//enter_ref
	
	return 0;
}

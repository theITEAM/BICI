// This file defines the data

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
 
using namespace std;

#include "model.hh"
#include "utils.hh"

/// Samples the initial conditions for the state
vector <InitCondValue> Model::initc_sample() const
{
	vector <InitCondValue> initc(nspecies);
	
	for(auto p = 0u; p < nspecies; p++){
		const auto &ic = species[p].init_cond;
		
		initc[p].type = ic.type;
		
		switch(ic.type){
		case POP_INIT: case ZEROPOP_INIT:
			initc[p].cpop = ic.cpop;
			break;
		
		
		case UNSET_INIT: emsg("Should not be here"); break;
		}
	}
	
	return initc;
}

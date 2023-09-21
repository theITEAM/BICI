// This provides information about one species within the model

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
using namespace std;

#include "species.hh"
#include "utils.hh"
#include "matrix.hh"


/// Initialises the model model.param,model.param_vec
Species::Species(const vector <double> &timepoint, const Details &details) : timepoint(timepoint), details(details){};


/// Creates a node tree used to sample Markovian events
void Species::create_markov_tree()
{
	for(auto i = 0u; i < markov_eqn.size(); i++){
		MarkovNode mn; mn.parent = UNSET;
		markov_tree.node.push_back(mn);
	}
	
	do{
		vector <unsigned int> list;
		for(auto i = 0u; i < markov_tree.node.size(); i++){
			if(markov_tree.node[i].parent == UNSET) list.push_back(i);
		}
		
		if(list.size() <= 1) break;
		for(auto j = 0u; j < list.size(); j += 2){	
			MarkovNode mn; 
			mn.parent = UNSET;
			
			auto i = list[j];
			markov_tree.node[i].parent = markov_tree.node.size();
			mn.child.push_back(i);
			
			if(j+1 < list.size()){
				i = list[j+1];
				markov_tree.node[i].parent = markov_tree.node.size();
				mn.child.push_back(i);
			}
			
			markov_tree.node.push_back(mn);
		}
	}while(true);
	
	if(false){
		for(auto i = 0u; i < markov_tree.node.size(); i++){
			const auto &no = markov_tree.node[i];
			cout << "Node " << i << ": parent " << no.parent << "  ";
			cout << "child "; for(auto j : no.child) cout << j << ", ";
			cout << endl;		
		}
	}
}

/// Updates the global compartment based on a classification compartment being set to a given value
unsigned int Species::update_c_comp(unsigned int c, unsigned int cl, unsigned int c_comp) const
{
	int cco = comp_gl[c].cla_comp[cl];
	int cnew = c +(c_comp - cco)*comp_mult[cl];
	if(cnew < 0 || cnew >= (int)N) emsg("Out of range4");
	return (unsigned int) cnew;
}

	
/// Prints out observation events 
void Species::print_obs_data(string name, const vector <ObsData> &obs) const 
{
	cout << name << endl;
	cout << "Observations:" << endl;
	for(auto &ob : obs){
		cout << ob.t << " ";
		switch(ob.type){
		case OBS_TRANS_EV: 
			cout << "OBS_TRANS_EV " << cla[ob.cl].tra[ob.tr].name;
			break;
			
		case OBS_COMP_EV:
			cout << "OBS_COMP_EV " << cla[ob.cl].comp[ob.c].name;
			break;
			
		case OBS_TEST_EV:
			{
				cout << "OBS_TEST_EV res=";
				if(ob.test_res == true) cout << "+ve"; else cout << "-ve";
				const auto &ds = source[ob.so];
				cout << "  Se=" << ds.obs_model.Se.te << "  Sp=" << ds.obs_model.Sp.te;
			}
			break;
			
		case OBS_POP:
			{
				cout << "OBS_POP " << ob.ref;
			}
			break;
			
		default: emsg("event data not right:"+ob.type); break;
		}
		cout << "   ";
	}
	cout << endl << endl;
}


/// Gets ti div time from an actual time
unsigned int Species::get_ti(double t) const
{
	return (unsigned int)(OVER_ONE*(t- details.t_start)/details.dt);
}


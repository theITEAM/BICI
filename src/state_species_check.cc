// Performs checks used in State

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
 
using namespace std;

#include "state_species.hh"
#include "utils.hh"

void StateSpecies::check(const vector <double> &popnum) const
{
	if(type == INDIVIDUAL){
		check_ind_ref();
		check_markov_eqn_ref();
		check_cpop();
		check_markov_eqn(popnum);
		check_markov_tree_rate();
		//check_event();
	}
}


/// Checks that markov_eqn_ref correctly references individuals in markov_eqn
void StateSpecies::check_markov_eqn_ref() const
{  
	vector <vector <bool> > done;
	done.resize(N);
	for(auto e = 0u; e < N; e++){
		done[e].resize(markov_eqn_vari[e].ind_tra.size(),false);
	}
	
	for(auto i = 0u; i < individual.size(); i++){
		const auto &ind = individual[i];
		for(auto j = 0u; j < ind.markov_eqn_ref.size(); j++){
			const auto &mer = ind.markov_eqn_ref[j];
			auto e = mer.e; auto index = mer.index;
			if(index >= markov_eqn_vari[e].ind_tra.size()) emsg("Out of range6");
			const auto &mit = markov_eqn_vari[e].ind_tra[index];
			if(mit.i != i) emsg("individual does not agree");	
			if(mit.index != j) emsg("individual index does not agree");
			done[e][index] = true;
		}
	}
	 
	for(auto e = 0u; e < N; e++){
		const auto &me = sp.markov_eqn[e];
		const auto &me_vari = markov_eqn_vari[e];
		
		const auto &ind_tra = me_vari.ind_tra;
		
		for(auto i = 0u; i < ind_tra.size(); i++){
			if(done[e][i] != true) emsg("Markov equation problem");
		}
		
		auto sum = 0.0; 
		if(me.source == true) sum = 1;
		else{
			for(auto i = 0u; i < ind_tra.size(); i++){
				const auto &it = ind_tra[i];
				
				const auto &ind = individual[it.i];
				auto val = get_indfac(ind,me);
				sum += val;
			}
		}
		
		if(dif(markov_eqn_vari[e].indfac_sum,sum)){
			emsg("Problem with ie_sum");
		}
	}
}

/// Checks that cpop is correctly set
void StateSpecies::check_cpop() const
{
	auto N = cpop.size();
	if(N != sp.comp_gl.size()) emsg("Number of compartments is wrong");
	
	vector <double> cpop_check(N,0);
	for(const auto &ind : individual){
		auto c = ind.c;
		if(c != UNSET && c != SINK){
			if(c >= N) emsg("c out of range");
			cpop_check[c]++;
		}
	}
	
	for(auto c = 0u; c < N; c++){
		if(dif(cpop[c],cpop_check[c])){
			emsg("cpop wrong");
		}
	}				
}


/// Checks that markov_tree_rate is correctly set
void StateSpecies::check_markov_eqn(const vector <double> &popnum) const
{
	for(auto e = 0u; e < N; e++){
		const auto &me = sp.markov_eqn[e];
		const auto &me_vari = markov_eqn_vari[e];
		
		auto value = eqn[me.eqn_ref].calculate(ti_sim,popnum,param_val,spline_val);	
		if(dif(value,me_vari.value)) emsg("value problem");
		
		auto sum = 0.0; 	
		if(me.source == true) sum = 1;
		else{	
			for(auto i = 0u; i < me_vari.ind_tra.size(); i++){
				const auto &it = me_vari.ind_tra[i];
				
				const auto &ind = individual[it.i];
				auto val = get_indfac(ind,me);
				sum += val;
			}
		}
		
		if(dif(me_vari.indfac_sum,sum)) emsg("indfac_sum problem");
	}
};


/// Checks that markov_tree_rate is correctly set
void StateSpecies::check_markov_tree_rate() const
{
	const auto &node = sp.markov_tree.node;
	
	vector <double> markov_tree_rate_check(node.size());
	
	for(auto i = 0u; i < N; i++){
		auto &me_vari = markov_eqn_vari[i];
		markov_tree_rate_check[i] = me_vari.value*me_vari.indfac_sum;
	}
	
	for(auto i = N; i < node.size(); i++){
		const auto &no = node[i];
		
		auto sum = 0.0; for(auto j : no.child) sum += markov_tree_rate_check[j];
		markov_tree_rate_check[i] = sum;
	}
	
	for(auto i = 0u; i < markov_tree_rate.size(); i++){
		if(dif(markov_tree_rate_check[i],markov_tree_rate[i])){
			emsg("markov tree problem");
		}
	}		
}

/*
/// Checks events on markov_eqn and individual agree
void StateSpecies::check_event() const
{
	vector < vector < vector <bool> > > done;
	done.resize(N);	
	for(auto e = 0u; e < N; e++){
		const auto &me_vari = markov_eqn_vari[e];
		done[e].resize(T);
		for(auto ti = 0u; ti < T; ti++){
			done[e][ti].resize(me_vari.div[ti].ev.size(),false);
		}
	}
	
	for(auto i = 0u; i < individual.size(); i++){
		const auto &ind = individual[i];
		for(auto j = 0u; j < ind.ev.size(); j++){
			const auto &ev = ind.ev[j];
			
			for(const auto &meir : ev.markov_eq_ind_ref){
				auto e = meir.e;
				auto ti = meir.ti;
				auto index = meir.index;
			
				const auto &meev = markov_eqn_vari[e].div[ti].ev[index];
		
				done[e][ti][index] = true;
			
				if(meev.t != ev.t) emsg("time does not agree");
			}
		}
	}

	for(auto e = 0u; e < N; e++){
		const auto &me_vari = markov_eqn_vari[e];
		for(auto ti = 0u; ti < T; ti++){
			const auto &ev = me_vari.div[ti].ev;
			for(auto j = 0u; j < ev.size(); j++){
				if(done[e][ti][j] != true){
					emsg("event not done");
				}
			}
		}
	}
}
*/

/// Prints all the events in a given markov equation
void StateSpecies::print_markov_eqn(unsigned int e) const
{
	const auto &me = sp.markov_eqn[e];
	const auto &me_vari = markov_eqn_vari[e];
	
	cout << "For Markov Equation: " << eqn[me.eqn_ref].te_raw << endl;
	for(auto ti = 0u; ti < T; ti++){
		cout << timepoint[ti] << " time" << endl;
		const auto &div = me_vari.div[ti];
		cout << "Value: " << div.value << "  indfac: "<< div.indfac_int << endl; 
		
		cout << endl;
	}
}


/// Prints information about an individual
void StateSpecies::print_ind(unsigned int i) const
{
	const auto &ind = individual[i];
	
	cout << ind.name << ":" << endl;
	for(const auto &ev : ind.ev){
		cout << ev.t << " ";
		switch(ev.type){
		case ENTER_EV: cout << "Enter"; break;
		case LEAVE_EV: cout << "Leave"; break;
		case MOVE_EV: cout << "Move " << sp.comp_gl[ev.c_gl].name; break;
		case NM_TRANS_EV: cout << "Mar " << sp.tra_gl[ev.tr_gl].name; break;
		case M_TRANS_EV: cout << "NMar " << sp.tra_gl[ev.tr_gl].name; break;
		}
		cout << endl;
	}
	
	for(auto cl = 0u;  cl < sp.ncla; cl++){
		cout << sp.cla[cl].name << " - ";
		for(const auto &ev : ind.cl_ev[cl]){
			cout << ind.ev[ev.index_ev].t << " ";
			if(ev.index_nmtrans != UNSET){
				cout << "NMT ";
			}
		}
		cout << endl;
	}
	cout << endl;
}


/// Checks that individual references are correct
void StateSpecies::check_ind_ref() const
{
	// TO DO
}



/// Ensures that the event sequence is cosistent
void StateSpecies::ensure_consistent(unsigned int c, vector <Event> &ev) const
{
	for(auto k = 0u; k < ev.size(); k++){
		auto &tra = sp.tra_gl[ev[k].tr_gl];
		if(tra.i != c){
			if(c > CODE) emsg("HH");
			auto tr_new = tra.transform[c];
			if(tr_new == UNSET) emsg("Prob");
			ev[k].tr_gl = tr_new;
			auto c_after = sp.tra_gl[tr_new].f;
			ev[k].c_after =c_after;
			c = c_after;
		}
		else c = tra.f;
	}
}


/// Given the event sequence, this recalculates cl_ev (events in different classifications
void StateSpecies::recalculate_cl_ev(Individual &ind) const 
{
	auto &cl_ev = ind.cl_ev;
	
	for(auto cl = 0u; cl < sp.ncla; cl++) cl_ev[cl].clear();
	
	auto &ev = ind.ev;
	for(auto k = 0u; k < ev.size(); k++){
		ClEvent cle; cle.index_ev = k; cle.index_nmtrans = UNSET; cle.nm_ref = UNSET;

		auto tr_gl = ev[k].tr_gl;

		switch(ev[k].type){
		case ENTER_EV: case LEAVE_EV:
			for(auto cl = 0u; cl < sp.ncla; cl++){
				cl_ev[cl].push_back(cle);
			}
			break;

		case MOVE_EV:
			cl_ev[tr_gl].push_back(cle);
			break;
		
		case NM_TRANS_EV: 
			emsg("TO DO nm trans");
			break;
		
		case M_TRANS_EV:	
			cl_ev[sp.tra_gl[tr_gl].cl].push_back(cle);
			break;
		}
	}
}


/// Checks that cl_ev on individual is correctly specified
void StateSpecies::check_cl_ev(Individual ind) const 
{
	auto cl_ev_store = ind.cl_ev;
	
	recalculate_cl_ev(ind);
	
	auto &cl_ev = ind.cl_ev; 
	for(auto cl = 0u; cl < sp.ncla; cl++){
		if(cl_ev[cl].size() != cl_ev_store[cl].size()) emsg("si pr");
		for(auto k = 0u; k < cl_ev[cl].size(); k++){
			if(cl_ev[cl][k].index_ev != cl_ev_store[cl][k].index_ev) emsg("prob 1");
			if(cl_ev[cl][k].index_nmtrans != cl_ev_store[cl][k].index_nmtrans) emsg("prob 2");
			if(cl_ev[cl][k].nm_ref != cl_ev_store[cl][k].nm_ref) emsg("prob 3");
		}
	}
}


/// Prints a sequence of individual events
void StateSpecies::print_event(string name, const Individual &ind) const
{
	auto c = ind.cinit;
	cout << name << "  " << ind.name << ":"<< endl;
	if(c == SOURCE) cout << "Source   ";
	else{
		if(c == UNSET) cout << "Unset  ";
		else cout << "[" << sp.comp_gl[c].name << "]  ";
	}
	
	for(auto &ev : ind.ev){
		switch(ev.type){
		case ENTER_EV: 
			cout << " ENTER " << ev.t << ", ";
			c = ev.c_gl;
			break;
		
		case LEAVE_EV:
			cout << " LEAVE " << ev.t << ", ";
			c = UNSET;
			break;
		
		case MOVE_EV:
			cout << " MOVE " << ev.t << ", ";
			c = sp.update_c_comp(c,ev.tr_gl,ev.c_gl);
			break;
			
		case NM_TRANS_EV:
			emsg("Not done");
			break;
			
		case M_TRANS_EV:	
			cout << " >> " << ev.t << ", " << sp.tra_gl[ev.tr_gl].name << " >> ";
			c = sp.tra_gl[ev.tr_gl].f;
			break;
		}
		
		if(c == SINK) cout << " SINK"; 
		else{
			if(c == UNSET) cout << " UNSET";
			else cout << "[" << sp.comp_gl[c].name << "]  ";
		}
	}
	cout << endl;
}


/// Moves a single event to a new time and generates the new event sequence
unsigned int StateSpecies::move_event(unsigned int cinit, vector <Event> &ev, unsigned int index, double t_new) const 
{
	auto flag = false;
	auto ev_store = ev[index];
	auto t = ev_store.t;
	
	if(t_new > t){
		while(index+1 < ev.size() && ev[index+1].t < t_new){
			ev[index] = ev[index+1];
			index++;
			flag = true;
		}							
	}
	else{
		while(index > 0 && ev[index-1].t > t_new){
			ev[index] = ev[index-1];
			index--;
			flag = true;
		}			
	}
	
	auto tr_gl_new = UNSET;
	if(flag == true){   // If the move event crosses an event in another classification
		ev[index] = ev_store;
		ev[index].t = t_new;
		
		ensure_consistent(cinit,ev);
		tr_gl_new = ev[index].tr_gl;
	}
	else{
		ev[index].t = t_new;
	}
	
	return tr_gl_new;
}


/// Checks that trans num is correct
void StateSpecies::check_transnum()
{
	for(auto i = 0u; i < pop_data_num.size(); i++){
		const auto &popd = sp.pop_data[i];
		cout << popd.t << " " << pop_data_num[i] << " " << popd.value << " " << popd.sd << " pop\n";
	}
	
	for(auto i = 0u; i < pop_trans_data_num.size(); i++){
		const auto &poptd = sp.pop_trans_data[i];
		cout << poptd.tmin << " " << poptd.tmax <<"  " << pop_trans_data_num[i] << " " << poptd.value << " " << poptd.sd<< " trans\n";
	}
}


/// Because a transition rate/mean has become negative this 
string StateSpecies::check_prior(const Equation &eq) const
{
	string type = "prior"; if(mode == SIM) type = "parameter value";
	
	if(eq.param_ref.size() == 1){
		const auto &pr = eq.param_ref[0];
		return " Please check the "+type+" for "+param[pr.th].full_name+" is correctly specified.";
	}		
	return " Please check the "+type+" is correctly specified.";
}


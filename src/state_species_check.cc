// Performs checks for state species (used under simulation) 

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
 
using namespace std;

#include "state_species.hh"
#include "utils.hh"

/// Check under simulation
void StateSpecies::check(unsigned int ti, const vector < vector <double> > &popnum_t)
{
	const auto &popnum = popnum_t[ti];
	
	if(type == INDIVIDUAL){
		check_markov_eqn_ref();
		check_cpop();
		check_markov_eqn(ti,popnum);
		check_markov_tree_rate();
	}
	else{
		check_trans_num_neg(); 
	}
	check_erlang();
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
		auto &me_vari = markov_eqn_vari[e];
		
		const auto &ind_tra = me_vari.ind_tra;
		
		for(auto i = 0u; i < ind_tra.size(); i++){
			if(done[e][i] != true) emsg("Markov equation problem");
		}
		
		auto sum = 0.0; 
		if(me.source == true) sum = me.source_tr_gl.size();
		else{
			for(auto i = 0u; i < ind_tra.size(); i++){
				const auto &it = ind_tra[i];
				
				const auto &ind = individual[it.i];
				auto val = get_indfac(ind,me);
				sum += val;
			}
		}
		
		if(dif(me_vari.indfac_sum,sum,dif_thresh)){
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
	for(auto i = 0u; i < individual.size(); i++){
		auto c = ind_sim_c[i];
		if(c != UNSET){
			if(c >= N) emsg("c out of range");
			cpop_check[c]++;
		}
	}
	
	for(auto c = 0u; c < N; c++){
		if(dif(cpop[c],cpop_check[c],dif_thresh)){
			//cout << c << " " << cpop[c] << " "<< cpop_check[c] <<" dif" << endl;
			emsg("cpop wrong");
		}
	}				
}


/// Checks that markov_eqn is correctly set
void StateSpecies::check_markov_eqn(unsigned int ti, const vector <double> &popnum) const
{
	const auto &precalc = param_val.precalc;
	
	for(auto e = 0u; e < N; e++){
		const auto &me = sp.markov_eqn[e];
		
		const auto &me_vari = markov_eqn_vari[e];
		
		double value;
		if(me.rate){
			value = eqn[me.eqn_ref].calculate(ti,popnum,precalc);	
		}
		else{
			auto mean = eqn[me.eqn_ref].calculate(ti,popnum,precalc);
			value = 1.0/mean;
		}
		
		if(dif(value,me_vari.value,dif_thresh)){
			emsg("value problem");
		}
		
		auto sum = 0.0; 	
		if(me.source == true) sum = me.source_tr_gl.size();
		else{	
			for(auto i = 0u; i < me_vari.ind_tra.size(); i++){
				const auto &it = me_vari.ind_tra[i];
				
				const auto &ind = individual[it.i];
				auto val = get_indfac(ind,me);
				sum += val;
			}
		}
		
		if(dif(me_vari.indfac_sum,sum,dif_thresh)) emsg("indfac_sum problem");
	}
};


/// Checks that markov_tree_rate is correctly set
void StateSpecies::check_markov_tree_rate()
{
	auto markov_tree_rate_store = markov_tree_rate;
	
	set_markov_tree_rate();

	for(auto i = 0u; i < markov_tree_rate.size(); i++){
		if(dif(markov_tree_rate_store[i],markov_tree_rate[i],dif_thresh)){
			emsg("markov tree problem"); 
		}
	}		
}


/// Prints all the events in a given markov equation
void StateSpecies::print_markov_eqn(unsigned int e) const
{
	const auto &me = sp.markov_eqn[e];
	const auto &me_vari = markov_eqn_vari[e];
	
	cout << "For Markov Equation: " << eqn[me.eqn_ref].te_raw << endl;
	for(auto ti = 0u; ti < T; ti++){
		cout << calc_t(ti,model.details) << " time" << endl;
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
		cout << ev.tdiv << " ";
		switch(ev.type){
		case ENTER_EV: cout << "Enter"; break;
		case LEAVE_EV: cout << "Leave"; break;
		case MOVE_EV: cout << "Move " << sp.comp_gl[ev.cl].name; break;
		case NM_TRANS_EV: cout << "NM " << sp.tra_gl[ev.tr_gl].name; break;
		case M_TRANS_EV: cout << "Markovian " << sp.tra_gl[ev.tr_gl].name; break;
		}
		cout << endl;
	}
}


/// Ensures that the event sequence is cosistent
void StateSpecies::ensure_consistent(vector <Event> &ev) const
{
	auto c = UNSET;
	for(auto k = 0u; k < ev.size(); k++){	
		if(sp.correct_ev(c,ev[k]) == true) emsg("Inconsistent");
		c = ev[k].c_after; 
	}
}


/// Prints a sequence of individual events
string StateSpecies::print_event(string name, const Individual &ind) const
{
	stringstream ss;
	ss << name << "  " << ind.name << ":"<< endl;

	ss << print_event(ind.ev,true);
	
	return ss.str();
}


/// Prints a vector of events
string StateSpecies::print_event(const vector <Event> &event, bool str) const
{
	stringstream ss;
	
	auto c = UNSET;
	for(auto &ev : event){
		auto c_after = ev.c_after;
	
		auto t = calc_t(ev.tdiv,model.details);
		
		switch(ev.type){
		case ENTER_EV: 
			ss << " ENTER " << t << ", ";
			c = ev.c_after;
			break;
		
		case LEAVE_EV:
			ss << " LEAVE " << t << ", ";
			c = UNSET;
			break;
		
		case MOVE_EV:
			ss << " MOVE " << t << ", ";
			c = sp.update_c_comp(c,ev.cl,ev.move_c);
			break;
			
		case M_TRANS_EV:	
			ss << " >> " << t << ", " << sp.tra_gl[ev.tr_gl].name << " >> ";
			c = sp.tra_gl[ev.tr_gl].f;
			break;
			
		case NM_TRANS_EV:
			ss << " >> " << t << ", " << sp.tra_gl[ev.tr_gl].name;
			ss << " (origin ";
			if(ev.e_origin == UNSET) cout << "UNSET";
			else{
				ss << "[" << ev.e_origin << "]";
			}
			ss << ") >> ";
			c = sp.tra_gl[ev.tr_gl].f;
			break;
		}
		
		if(c != c_after){
			emsg("after problem ere");
		}
			
		if(c == UNSET) ss << " UNSET";
		else ss << "[" << sp.comp_gl[c].name << "]";

		const auto &iif = ev.ind_inf_from;
		ss << "<" << iif.p << " " << iif.i << ">";

		if(ev.inf_node_ref != UNSET){
			const auto &iif = ev.ind_inf_from;
			if(iif.p == ENTER_INF) ss << "(n=" << ev.inf_node_ref << ",ENT)";
			else{
				if(iif.p == OUTSIDE_INF) ss << "(n=" << ev.inf_node_ref << ",OUT)";
				else ss << "(n=" << ev.inf_node_ref << ",i=" << iif.i << ")";
			}
		}
		
		if(ev.observed) ss << "(obs)"; else ss << "(unobs)";
		ss << "  ";
	}
	ss << endl;
	
	if(str == false) cout << ss.str();
	
	return ss.str();
}


/// Checks the prior
string StateSpecies::check_prior(const Equation &eq) const
{
	string type = "prior is"; if(mode == SIM) type = "parameter values are";
	
	if(eq.param_ref.size() == 1){
		const auto &pr = eq.param_ref[0];
		return " Please check the "+type+" correctly set for "+model.param[pr.th].full_name;
	}		
	return " Please check the "+type+" correctly set.";
}


/// Sets ind_sim_c at time t based on individual timelines 
void StateSpecies::set_ind_sim_c(unsigned int ti)
{
	double t = ti;

	auto imax = individual.size();
	if(ind_sim_c.size() != imax) ind_sim_c.resize(imax);
	for(auto i = 0u; i < imax; i++){
		const auto &ind = individual[i];
		
		auto c = UNSET;
		auto e = 0u; while(e < ind.ev.size() && ind.ev[e].tdiv <= t){ c = ind.ev[e].c_after; e++;}
		
		ind_sim_c[i] = c;
	}
}


/// Determines if all events have been set correctly
bool StateSpecies::all_events_correct(unsigned int i, unsigned int cl)
{
	const auto &ev = individual[i].ev;
	
	for(const auto &eve : ev){
		if(eve.cl == cl && !eve.observed) return false;
	}
	
	return true;
}


/// Checks initial condtion has zero for internal state for IC
void StateSpecies::check_erlang() const 
{
	for(auto c = 0u; c < sp.comp_gl.size(); c++){
		const auto &co = sp.comp_gl[c];
		if(co.erlang_hidden){
			if(init_cond_val.cnum[c] != 0) emsg("Erlang must have zero initial state for internal compartments");
		}
	}
}


/// Checks if the transition number is negative
void StateSpecies::check_trans_num_neg() const
{
	for(auto tr = 0u; tr < trans_num.size(); tr++){
		for(auto ti = 0u; ti < trans_num[tr].size(); ti++){
			if(trans_num[tr][ti] < 0) emsg("Negative trans num");
		}
	}
}


/// Adds an algorithm warning
void StateSpecies::add_alg_warn(string te)
{
	add_alg_warning(te,UNSET,alg_warn);
}


// This deals with individual-based updates

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
// NOTES
// markov_eqn list all transitions which are Markovian. Attached to this list are the individuals 

using namespace std;

#include "state_species.hh"
#include "utils.hh"
#include "matrix.hh"


/// Updates the individuals in the system (using a modified Gillespie algorithm)
void StateSpecies::update_individual_based(unsigned int ti_, const vector <double> &popnum)
{	
	ti_sim = ti_;
	
	update_markov_eqn_value(popnum);
	
	sort_trig_event();
	
	auto i_trig = 0u;  // This indexes potential trigger events
	
	const auto &node = sp.markov_tree.node;
	
	auto &tri_ev = trig_div[ti_sim].ev;
	
	auto t = timepoint[ti_sim];
	auto tnext = timepoint[ti_sim+1];
	
	do{
		double R = 0, dt;
		
		if(nnode > 0) R = markov_tree_rate[nnode-1];

		if(R > -TINY && R < TINY) dt = LARGE;
		else{
			if(R < 0) emsg("Negative rate");
			auto ra = ran(); if(ra == 0) ra = TINY;
			dt = -log(ra)/R;
		}
		
		auto tf = t + dt;
		
		if(i_trig < tri_ev.size() && tf > tri_ev[i_trig].t){ // Trigger event
			const auto &trig = tri_ev[i_trig];
			tf = trig.t;
			auto i = trig.i;
			
			switch(trig.type){
			case ENTER_EV:
				add_event(ENTER_EV,i,UNSET,trig.c,trig.c,t);	
				update_ind_c(trig.c,i,tf,popnum);
				break;
				
			case LEAVE_EV:
				update_ind_remove(i,tf);
				break;
			
			case MOVE_EV:
				update_ind_move(i,tf,trig.c,trig.trg,popnum); 
				break;
				
			case M_TRANS_EV:
				emsg("Should not have event");
				break;
			
			case NM_TRANS_EV:
				IndTransRef itr; itr.i = trig.i; itr.index = UNSET; itr.tr_gl = trig.trg;	
				
				update_ind_trans(tf,itr,popnum);
				break;
			}
			i_trig++;
		}
		else{ 																													 // Markovian event
			if(tf > tnext) return;
			
			// Randomly selects an equation
			auto num = ran()*R;
			auto e = nnode-1;
			while(e >= N){
				auto &no = node[e];
				
				if(no.child.size() == 1) e = no.child[0];
				else{
					if(no.child.size() != 2) emsg("prob");
					auto val0 = markov_tree_rate[no.child[0]];
					auto val1 = markov_tree_rate[no.child[1]]; 
					if(num > val0+val1 || num < -TINY) emsg("Not in range");
					
					if(num < val0) e = no.child[0];
					else{ e = no.child[1]; num -= val0;}
				}
			}
			
			// Randomly selects an individual
			auto &me = sp.markov_eqn[e];
			auto &me_vari = markov_eqn_vari[e];
			
			if(me.source == true){
				//print("vec",me.source_tr_gl);
				auto si = me.source_tr_gl.size();
				if(si == 0) emsg("Source size zero");
				auto tr_gl = me.source_tr_gl[(unsigned int)(ran()*si)];
				
				auto name = get_new_ind_name("Source ",individual.size()+1);
				
				auto i = individual.size();
				
				add_individual(name,SOURCE,true);
	
				IndTransRef itr; itr.i = i; itr.index = UNSET, itr.tr_gl = tr_gl;
				update_ind_trans(tf,itr,popnum);
			}
			else{
				if(me_vari.ind_tra.size() == 0) emsg("should not be zero");
				
				unsigned int j;
				if(me.ind_variation == false){ // Randomly selects an individual
					j = (unsigned int)(ran()*me_vari.ind_tra.size());
				}	
				else{			                     // Selects weighted by individual factor
					auto &it = me_vari.ind_tra;
					auto sum = 0.0;
					vector <double> sum_st;
					for(auto k = 0u; k < it.size(); k++){
						const auto &itr = it[k];
						sum += get_indfac(individual[itr.i],me);
						sum_st.push_back(sum);
					}
					if(sum == 0) emsg("zero problem");
					auto z = ran()*sum;
					
					j = 0u; while(j < it.size() && z > sum_st[j]) j++;
					if(j == it.size()) emsg("erro select");
					
					if(dif(sum,me_vari.indfac_sum)) emsg("SHould be the same");
				}
			
				update_ind_trans(tf,me_vari.ind_tra[j],popnum);
			}
		}
		
		check(popnum);
		
		t = tf;
	}while(true);
	
}


/// Adds all the initial individuals to system
void StateSpecies::activate_initial_state(const vector <double> &popnum)
{
	const auto &node = sp.markov_tree.node;

	markov_tree_rate.clear();
	markov_tree_rate.resize(node.size(),0);
	for(auto i = 0u; i < N; i++){
		auto &me = sp.markov_eqn[i];
		auto &me_vari = markov_eqn_vari[i];
		
		const auto &eq = eqn[me.eqn_ref]; 
		auto value = eq.calculate(0,popnum,param_val,spline_val);
		
		me_vari.value = value;
		me_vari.indfac_sum = 0; if(me.source == true) me_vari.indfac_sum = 1;
		update_markov_tree_rate(i,value*me_vari.indfac_sum);
	}
	
	for(auto i = 0u; i < individual.size(); i++){
		const auto &ind = individual[i];
		if(ind.cinit != UNSET) update_ind_c(ind.cinit,i,details.t_start,popnum);
	}
	
	//for(auto i = 0u; i < N; i++){
		//auto &me_vari = markov_eqn_vari[i];
		//me_vari.ind_tra_init = me_vari.ind_tra;
	//}
}


/// Updates individual based on it entering compartment c
void StateSpecies::update_ind_c(unsigned int c, unsigned int i, double t, const vector <double> &popnum)
{
	auto &ind = individual[i];
	
	ind.c = c;
	if(c == SINK) return;

	if(c == UNSET) emsg("Should not be unset");

	cpop[c]++;
	
	const auto &co = sp.comp_gl[c];
	
	for(auto cl = 0u; cl < sp.ncla; cl++){
		const auto &tinfo = co.trainfo[cl];
		
		const auto &tref = tinfo.tra_ref;
		
		if(tref.size() > 0){
			if(tinfo.branch == false && tref.size() > 1){     // If markov rates leaving a compartment
				for(const auto &tre : tref){
					auto tgl = tre.tr;
					const auto &tr = sp.tra_gl[tgl];
					if(tr.type != EXP_RATE) emsg("Should be exponential"); 
					if(tr.i != c) emsg("Problem with selected");
			
					add_markov_transition(i,tgl);
				}
			}
			else{
				unsigned int tgl;
		
				if(tinfo.branch == true){                       // If braching then select branch
					if(tref.size() <= 1) emsg("Problem1");
					
					vector <double> sum_store;
					
					auto sum = 0.0;
					for(auto k = 0u; k < tref.size(); k++){
						auto th = tref[k].bp_th;			
						sum += param_val[th];
						sum_store.push_back(sum);
					}
					
					if(dif(sum,1)) emsg("Sum is not 1");
					
					auto z = ran();
					auto k = 0u; while(k < tref.size() && z > sum_store[k]) k++;
					if(k == tref.size()) emsg("Sampling problem");
					
					tgl = tref[k].tr;
				}
				else{                                          // If just a sinle branch
					if(tref.size() > 1) emsg("Problem2");
					tgl = tref[0].tr;
				}
				
				const auto &tr = sp.tra_gl[tgl];
				const auto &dp = tr.dist_param;
				
				if(tr.i != c) emsg("Problem with selected");
				
				switch(tr.type){
				case EXP_MEAN: emsg("SHould not be mean"); break;
				
				case EXP_RATE: 
					{
						add_markov_transition(i,tgl);
					}
					break;
			
				case GAMMA:
					{
						auto mean = eqn[dp[0].eq_ref].calculate_indfac(ind,ti_sim,popnum,param_val,spline_val);
						auto cv = eqn[dp[1].eq_ref].calculate_indfac(ind,ti_sim,popnum,param_val,spline_val);
						auto ts = t+gamma_sample(mean,cv);
						
						TrigEvent trig; 
						trig.type = NM_TRANS_EV; trig.i = i; trig.c = UNSET; trig.trg = tgl; trig.t = ts;
						insert_trigger_event(trig);
					}
					break;
				
				case ERLANG:
					{
						auto mean = eqn[dp[0].eq_ref].calculate_indfac(ind,ti_sim,popnum,param_val,spline_val);
						auto shape = eqn[dp[1].eq_ref].calculate_indfac(ind,ti_sim,popnum,param_val,spline_val);
						auto ts = t+gamma_sample(mean,sqrt(1.0/shape));
						TrigEvent trig; 
						trig.type = NM_TRANS_EV; trig.i = i; trig.c = UNSET; trig.trg = tgl; trig.t = ts;
						insert_trigger_event(trig);
					}
					break;

				case LOG_NORMAL:
					{
						auto mean = eqn[dp[0].eq_ref].calculate_indfac(ind,ti_sim,popnum,param_val,spline_val);
						auto cv = eqn[dp[1].eq_ref].calculate_indfac(ind,ti_sim,popnum,param_val,spline_val);
						auto ts = t+lognormal_sample(mean,cv);
					
						TrigEvent trig; 
						trig.type = NM_TRANS_EV; trig.i = i; trig.c = UNSET; trig.trg = tgl; trig.t = ts;
						insert_trigger_event(trig);
					}
					break;
					
				case WEIBULL:
					{
						auto scale = eqn[dp[0].eq_ref].calculate_indfac(ind,ti_sim,popnum,param_val,spline_val);
						auto shape = eqn[dp[1].eq_ref].calculate_indfac(ind,ti_sim,popnum,param_val,spline_val);
						auto ts = t+weibull_sample(shape,scale);
				
						TrigEvent trig; 
						trig.type = NM_TRANS_EV; trig.i = i; trig.c = UNSET; trig.trg = tgl; trig.t = ts;
						insert_trigger_event(trig);
					}
					break;
				}	
			}
		}
	}
}


/// Adds an individual the markovian transition corresponding to a global transition
void StateSpecies::add_markov_transition(unsigned int i, unsigned int tgl)
{
	auto &ind = individual[i];
	
	auto mef = sp.tra_gl[tgl].markov_eqn_ref;
	
	auto &me = sp.markov_eqn[mef];
	auto &me_vari = markov_eqn_vari[mef];
	
	MarkovEqnRef mer; mer.e = mef; mer.index = me_vari.ind_tra.size();
	IndTransRef itr; itr.i = i; itr.index = ind.markov_eqn_ref.size(); itr.tr_gl = tgl;

	ind.markov_eqn_ref.push_back(mer);
	me_vari.ind_tra.push_back(itr);
	
	auto val = get_indfac(ind,me);
	me_vari.indfac_sum += val;
	
	//add_markov_eqn_ev(i,t,mef,ADD_INDFAC);
	
	update_markov_tree_rate(mef,me_vari.value*val);
}


/// Updates markov_tree_rate
void StateSpecies::update_markov_tree_rate(unsigned int e, double dif)
{
	const auto &node = sp.markov_tree.node;

	do{
		markov_tree_rate[e] += dif;
		e = node[e].parent;
	}while(e != UNSET);
}


/// Inserts a future non-Markovian event
void StateSpecies::insert_trigger_event(const TrigEvent &trig)
{
	auto t = trig.t;
	if(t >= details.t_end){
		end_nm_event.push_back(trig);
	}
	else{
		auto j = (unsigned int)((trig.t-details.t_start)/details.dt);
		if(j >= trig_div.size()) emsg("nm out of range8");
		
		auto &ev = trig_div[j].ev;
		
		if(ti_sim == j){  // If adding to present time segment then must do it in an ordered way
			auto k = 0u; while(k < ev.size() && t > ev[k].t) k++;
			ev.insert(ev.begin()+k,trig);
		}
		else{         // If added to a later time segment ordering is done later
			ev.push_back(trig);
		}
	}
}

/// Updates the value of the markov equations
void StateSpecies::update_markov_eqn_value(const vector <double> &popnum)
{
	for(auto e = 0u; e < sp.markov_eqn.size(); e++){
		auto &me = sp.markov_eqn[e];
		auto &me_vari = markov_eqn_vari[e];
		
		auto value = eqn[me.eqn_ref].calculate(ti_sim,popnum,param_val,spline_val);
	
		if(me_vari.value != value){
			update_markov_tree_rate(e,(value-me_vari.value)*me_vari.indfac_sum);
			me_vari.value = value;
		}
	}
}


// Used to order trigger events
bool Trig_ord(TrigEvent ev1, TrigEvent ev2)                      
{ return (ev1.t < ev2.t); };  


/// Sorts the triger events in the time division
void StateSpecies::sort_trig_event()
{
	auto &tri_ev = trig_div[ti_sim].ev;
	sort(tri_ev.begin(),tri_ev.end(),Trig_ord);
}


/// Updates an individual with a transition
void StateSpecies::update_ind_trans(double t, const IndTransRef &itr, const vector <double> &popnum)
{
	auto i = itr.i;
	
	auto &ind = individual[i];

	/// Updates the individual compartment based on transition
	auto tr_gl = itr.tr_gl;
	
	if(ind.c != sp.tra_gl[tr_gl].i){
		tr_gl = sp.tra_gl[tr_gl].transform[ind.c];
		//if(!(ind.c == UNSET && tr.i == SOURCE)){
		cout << "initial comp trans not agree\n"; 
		//emsg("initial comp trans not agree");
		//}
	}
	
	const auto &tr = sp.tra_gl[tr_gl];	
	
	auto type = NM_TRANS_EV; if(tr.markov_eqn_ref != UNSET) type = M_TRANS_EV;
	
	if(ind.c != UNSET && ind.c != SOURCE) cpop[ind.c]--;
	
	auto c_new = tr.f;

	//if(type == NM_TRANS_EV) emsg("errpr2");
	add_event(type,i,tr_gl,UNSET,c_new,t);
	
	remove_markov_trans(i);
	
	/// Updates individual based on new compartment
	update_ind_c(c_new,i,t,popnum);
}


// Removes all current markov transitions for an individual
void StateSpecies::remove_markov_trans(unsigned int i)
{
	Individual &ind = individual[i];
	for(const auto &mef : ind.markov_eqn_ref){
		auto &me = sp.markov_eqn[mef.e];
		auto &me_vari = markov_eqn_vari[mef.e];
		
		auto indfac = get_indfac(ind,me);
		
		me_vari.indfac_sum -= indfac;
		
		update_markov_tree_rate(mef.e,-me_vari.value*indfac);
		
		auto &ind_tra = me_vari.ind_tra;
		auto k = mef.index;
	
		if(k+1 < ind_tra.size()){
			ind_tra[k] = ind_tra[ind_tra.size()-1];
			
			const auto &it = ind_tra[k];
			individual[it.i].markov_eqn_ref[it.index].index = k;
		}
		ind_tra.pop_back();
	}
	ind.markov_eqn_ref.clear();
}


/// Removes an individual from the system
void StateSpecies::update_ind_remove(unsigned int i, double t)
{
	auto &ind = individual[i];
	if(ind.c < CODE) cpop[ind.c]--;
	
	ind.c = UNSET;
	
	add_event(LEAVE_EV,i,UNSET,UNSET,UNSET,t);
	
	remove_markov_trans(i);
}


/// Moves an individual to a particular compartment in a given classification
void StateSpecies::update_ind_move(unsigned int i, double t, unsigned int c_comp, unsigned int cl, const vector <double> &popnum)
{
	auto &ind = individual[i];
	
	int c = ind.c;
	if(c == UNSET) return;
	
	cpop[c]--;
	
	auto cnew = sp.update_c_comp(c,cl,c_comp);
	
	add_event(MOVE_EV,i,cl,c_comp,cnew,t);	

	remove_markov_trans(i);
	
	update_ind_c(cnew,i,t,popnum);
}
		
				
/// Adds event individual timeline
void StateSpecies::add_event(EventType type, unsigned int i, unsigned int tr_gl, unsigned int c_gl, unsigned int c_after, double t)
{
	if(mode == INF && type == NM_TRANS_EV) emsg("errpr");
	Event ev; 
	ev.type = type;
	ev.tr_gl = tr_gl;
	ev.c_gl = c_gl;
	ev.t = t;
	ev.c_after = c_after;
	ev.fixed = false;
	
	auto &ind = individual[i];
	
	/*
	if(type == M_TRANS_EV){
		auto &tra = sp.tra_gl[tr_gl];
		if(tra.i == SOURCE) ind.birth = t;
		if(tra.f == SINK) ind.death = t;
	}
	*/
	
	ClEvent cl_ev; cl_ev.index_ev = ind.ev.size(); cl_ev.nm_ref = UNSET; cl_ev.index_nmtrans = UNSET;
	
	ind.ev.push_back(ev);
			
	switch(type){
	case ENTER_EV: case LEAVE_EV:
		for(auto cl = 0u; cl < sp.ncla; cl++){
			ind.cl_ev[cl].push_back(cl_ev);
		}
		break;
	
	case MOVE_EV:
		ind.cl_ev[tr_gl].push_back(cl_ev);
		break;
		
	case NM_TRANS_EV: 
		{	
			const auto &trang = sp.tra_gl[tr_gl];
			auto cl = trang.cl;
			
			NMTransInd nmti; 
			nmti.t = t; nmti.i = i; nmti.cl = cl; nmti.index = ind.cl_ev[cl].size();
			
			auto nm_ref = trang.nm_ref;
			if(nm_ref == UNSET) emsg("Problem nm ref");
			
			cl_ev.nm_ref = nm_ref;
			cl_ev.index_nmtrans = nm_trans_ind[nm_ref].size();
			
			nm_trans_ind[nm_ref].push_back(nmti);
			ind.cl_ev[cl].push_back(cl_ev);
		}
		break;
		
	case M_TRANS_EV:	
		ind.cl_ev[sp.tra_gl[tr_gl].cl].push_back(cl_ev);
		break;
	}
}

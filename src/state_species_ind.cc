// This deals with individual-based updates when performing simulation

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
using namespace std;

#include "state_species.hh"
#include "utils.hh"
#include "matrix.hh"

/// Updates the individuals in the system (using a modified Gillespie algorithm)
void StateSpecies::update_individual_based(unsigned int ti, const vector < vector <Poss> > &pop_ind, const vector < vector <double> > &popnum_t)
{
	auto dt = details.dt;

	update_markov_eqn_value(ti,popnum_t);

	sort_trig_event(ti);
	
	auto i_trig = 0u;                            // This indexes potential trigger events
	
	const auto &node = sp.markov_tree.node;
	
	auto &tri_ev = trig_div[ti].ev;

	double t = ti;
	double tnext = ti+1;
	do{
		double dtdiv;
		auto R = 0.0; if(nnode > 0) R = markov_tree_rate[nnode-1];
	
		if(R > -SMALL && R < TINY) dtdiv = LARGE;
		else{
			if(R < 0){		
				emsg("Negative rate2");
			}
			auto ra = ran(); if(ra == 0) ra = TINY;
			dtdiv = -log(ra)/(R*dt);
		}
		
		auto tf = t + dtdiv;
		
		if(tf < tnext){
			
			if(false){
				if(t == ti && ran() < 1){
					auto st = tf;
					auto tfi = (unsigned int)(tf);
					tf = tfi; 		
					cout << st << tf << endl;					
				}
			}
		}
		
		if(i_trig < tri_ev.size() && tf > tri_ev[i_trig].tdiv){ // Trigger event
			const auto &trig = tri_ev[i_trig];
			tf = trig.tdiv;
			auto i = trig.i;
			
			switch(trig.type){
			case ENTER_SIM_EV:
				{
					IndInfFrom iif;
					if(sp.comp_gl[trig.c].infected == true) iif.p = ENTER_INF;
					add_event(ENTER_EV,i,UNSET,UNSET,UNSET,trig.c,t,iif);	
					update_ind_c(i,tf,UNSET,popnum_t);
				}
				break;
				
			case LEAVE_SIM_EV:
				update_ind_remove(i,tf);
				break;
			
			case MOVE_EV:
				update_ind_move(i,tf,trig.c,trig.trg,popnum_t); 
				break;
				
			case M_TRANS_SIM_EV:
				emsg("Should not have event");
				break;
			
			case NM_TRANS_SIM_EV:
				{
					IndTransRef itr; itr.i = trig.i; itr.index = UNSET; itr.tr_gl = trig.trg;	
					update_ind_trans(tf,itr,popnum_t,pop_ind);
				}
				break;
				
			case DATA_TRANS_SIM_EV:
				{
					auto ty = trig.type;
					vector <SimTrigEvent> trig_vec;
					trig_vec.push_back(trig);
					
					while(i_trig+1 < tri_ev.size() && tri_ev[i_trig+1].tdiv == tf && 
					      tri_ev[i_trig+1].type == ty && tri_ev[i_trig+1].i == i){
						i_trig++;
						trig_vec.push_back(tri_ev[i_trig]);
					}
				
					add_data_event(i,tf,trig_vec,popnum_t,pop_ind);
				}
				break;
				
			case SOURCE_SIM_EV:
				{
					IndTransRef itr; itr.i = trig.i; itr.index = UNSET; itr.tr_gl = trig.trg;	
					update_ind_trans(tf,itr,popnum_t,pop_ind);
				}
				break;
			}
			i_trig++;
		}
		else{ 														   				    // Markovian event
			if(tf > tnext) return;
			
			auto num = ran()*R;
		
			auto e = nnode-1;                             // Randomly selects an equation
			
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
				auto si = me.source_tr_gl.size();
				if(si == 0) emsg("Source size zero");
				auto tr_gl = me.source_tr_gl[(unsigned int)(ran()*si)];
				
				auto i = individual.size();
				
				IndTransRef itr; itr.i = i; itr.index = UNSET, itr.tr_gl = tr_gl;
				
				if(allow_event(tf,itr)){
					add_individual(UNOBSERVED_IND);
					update_ind_trans(tf,itr,popnum_t,pop_ind);
				}
			}
			else{
				if(me_vari.ind_tra.size() == 0) emsg("should not be zero");
				
				unsigned int j;
				if(me.ind_variation == false){              // Randomly selects an individual
					j = (unsigned int)(ran()*me_vari.ind_tra.size());
				}	
				else{			                                  // Selects weighted by individual factor
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
					
					if(dif(sum,me_vari.indfac_sum,dif_thresh)){
						emsg("Should be the same");
					}
				}
			
				const auto &itr = me_vari.ind_tra[j];
				
				if(allow_event(tf,itr)){
					update_ind_trans(tf,itr,popnum_t,pop_ind);
				}
			}
		}
		
		if(slow_check) check(ti,popnum_t);
		//check(ti,popnum_t);
		
		t = tf;
	}while(true);

	if(testing) check(T,popnum_t);
}

						
/// Allows events or not
// Under interence this allows initial simulate state to lie close to data
bool StateSpecies::allow_event(double t, const IndTransRef &itr) const
{
	auto trg = itr.tr_gl;
	
	auto ti = get_ti(t);
	
	if(false){
		for(auto e : sp.obs_trans_eqn_ref[trg][ti]){
			cout << obs_trans_eqn_value[e][ti] << " " ;
		}
		cout << " hh" << endl;
	}

	auto prob_ac = exp(trans_logprob_notobs(trg,ti));

	// Don't allow event if sink and still trigger events
	const auto &tr = sp.tra_gl[trg];
	
	if(mode == INF && tr.variety == SINK_TRANS && itr.i < sp.individual.size() && t < sp.individual[itr.i].tdivmax){
		return false;
	}
	
	if(ran() < prob_ac) return true;
	
	return false;
}


/// Adds data event to the individual time line
void StateSpecies::add_data_event(unsigned int i, double t, const vector <SimTrigEvent> &trig_vec, const vector < vector <double> > &popnum_t, const vector < vector <Poss> > &pop_ind)
{
	auto c = ind_sim_c[i];

	vector <DataTransPos> tra_pos;
	auto sum = 0.0;
	for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
		const auto &trg = sp.tra_gl[tr];
	
		if(trg.i == c && trg.ev_type == M_TRANS_EV){ 
			auto prob = 0.0;
		
			auto flag = false;
			double pr;
			for(auto k = 0u; k < trig_vec.size(); k++){			
				pr = like_trans_prob_obs(tr,sp.individual[i].obs[trig_vec[k].c]);				 
				if(pr == LI_WRONG) flag = true;
				else prob += pr;
			}
			
			if(flag == false){						
				auto val = exp(prob);
				sum += val;
				
				DataTransPos tp; tp.trg = tr; tp.prob = val; tp.probsum = sum; 
				tra_pos.push_back(tp);
			}
		}
	}

	if(sum > 0){
		auto z = ran()*sum;
		auto k = 0u; while(k < tra_pos.size() && z > tra_pos[k].probsum) k++;
		if(k == tra_pos.size()) emsg("tra_pos proble");
		
		IndTransRef itr; itr.i = i; itr.index = UNSET; itr.tr_gl = tra_pos[k].trg;
		update_ind_trans(t,itr,popnum_t,pop_ind);
	}
}


/// Adds all the initial individuals to system
void StateSpecies::activate_initial_state(double t, const vector < vector <double> > &popnum_t)
{
	auto ti = get_ti(t);
	const auto &popnum = popnum_t[ti];
	
	const auto &precalc = param_val.precalc;
	
	const auto &node = sp.markov_tree.node;

	markov_tree_rate.clear();
	markov_tree_rate.resize(node.size(),0);
	for(auto i = 0u; i < N; i++){
		auto &me = sp.markov_eqn[i];
		auto &me_vari = markov_eqn_vari[i];
		
		const auto &eq = eqn[me.eqn_ref]; 
		
		double value;
		if(me.rate){
			value = eq.calculate(0,popnum,precalc);
		}
		else{
			auto mean = eq.calculate(0,popnum,precalc);
			value = 1/mean;
		}
		
		auto indfac_sum = 0;
		if(me.source == true) indfac_sum = me.source_tr_gl.size();
		
		me_vari.value = value;
		me_vari.indfac_sum = indfac_sum; 

		update_markov_tree_rate(i,value*indfac_sum);
	}
	
	for(auto i = 0u; i < individual.size(); i++){
		if(ind_sim_c[i] != UNSET){
			update_ind_c(i,t,UNSET,popnum_t);
		}
	}
}


/// Updates individual based on it entering compartment c
void StateSpecies::update_ind_c(unsigned int i, double t, unsigned int cl_trans, const vector < vector <double> > &popnum_t)
{
	auto &ind = individual[i];
	if(ind.ev.size() == 0) emsg("Should have an event");
	
	auto c = ind_sim_c[i];
	
	if(c == UNSET) emsg("Should not be unset");
	cpop[c]++;
	
	const auto &co = sp.comp_gl[c];
		
	for(auto cl = 0u; cl < sp.ncla; cl++){
		const auto &tlg = co.tra_leave_group[cl];
		
		const auto &tref = tlg.tr_list;
		
		if(tref.size() > 0){
			if(tlg.markov){                               // If markov rates leaving a compartment
				for(auto tgl : tref){
					const auto &tr = sp.tra_gl[tgl];
					if(tr.ev_type != M_TRANS_EV) emsg("Should be exponential"); 
					if(tr.i != c) emsg("Problem with selected");
			
					add_markov_transition(i,tgl);
				}
			}
			else{
				if(cl == cl_trans || cl_trans == UNSET){
					auto flag = false;

					// Tries to insert a data event
					if(mode == INF) flag = try_insert_data_trans_event(t,cl,i);
				
					if(flag == false){
						auto trig = get_nm_trig_event(t,i,c,cl,popnum_t);
						
						auto t = trig.tdiv;
						if(t < details.T){
							IndTransRef itr; itr.i = trig.i; itr.index = UNSET; itr.tr_gl = trig.trg;	
							if(allow_event(t,itr)){
								insert_trigger_event(trig);
							}
						}
					}
				}
			}
		}
	}
}

				
/// Gets a trigger event for a NM transition
SimTrigEvent StateSpecies::get_nm_trig_event(double t, unsigned int i, unsigned int c, unsigned int cl, const vector < vector <double> > &popnum_t)
{
	SimTrigEvent trig;
	
	auto &ind = individual[i];
	
	const auto &ev = ind.ev;
	
	if(ev.size() == 0) emsg("Should not be zero");
	
	auto e = ev.size()-1;
	
	auto t_begin = t;
	
	auto dt = details.dt;
	
	// If mid-way through NM trans then find origin (used for PPC)

	if(t_begin != ev[e].tdiv){
		while(e > 0){
			auto &eve = ev[e];
			if(eve.cl == cl && (eve.type == M_TRANS_EV || eve.type == NM_TRANS_EV)) break;
			e--;
		}
		t_begin = ev[e].tdiv;

		c = ev[e].c_after;
	}
	
	auto ti = get_ti(t_begin);
	const auto &popnum = popnum_t[ti];
	
	const auto &tlg = sp.comp_gl[c].tra_leave_group[cl];
		
	unsigned int tgl;
	
	const auto &tref = tlg.tr_list;
	
	const auto &precalc = param_val.precalc;
	
	auto loop = 0u;
	do{
		if(tlg.branch == true){                         // If branching then select branch	
			if(tref.size() <= 1) emsg("Problem1");
			
			vector <double> sum_store;
			
			auto sum = 0.0;
			for(auto k = 0u; k < tref.size(); k++){
				const auto &trg = sp.tra_gl[tref[k]];
				
				auto eq_ref = trg.bp.eq_ref;
				double prob;
				if(eq_ref == UNSET){
					prob = 1.0;
					for(auto e : trg.bp_other_eq){
						prob -= eqn[e].calculate_indfac(ind,ti,popnum,precalc);
					}
				}
				else{
					prob = eqn[eq_ref].calculate_indfac(ind,ti,popnum,precalc);
				}
			
				if(prob < 0) prob = 0;
				if(prob > 1) prob = 1;
									
				sum += prob;
				sum_store.push_back(sum);
			}
			
			if(!tlg.all_branches && dif(sum,1,dif_thresh)){
				const auto &trg = sp.tra_gl[tref[0]];
				const auto &tr = sp.cla[trg.cl].tra[trg.tr];
				const auto &co = sp.cla[trg.cl].comp[tr.i];
				run_error("Branching probabilities leaving compartment '"+co.name+"' sum to greater than one");
			}
			
			auto z = sum*ran();
		
			auto k = 0u; while(k < tref.size() && z > sum_store[k]) k++;
			if(k == tref.size()) emsg("Sampling problem");
			
			tgl = tref[k];		
		}
		else{                                           // If just a single branch
			if(tref.size() > 1) emsg("Problem22");
			tgl = tref[0];
		}
		
		const auto &tr = sp.tra_gl[tgl];

		const auto &dp = tr.dist_param;		
			
		string warn;
		
		switch(tr.type){	
		case EXP_RATE: case EXP_MEAN: emsg("Should not be rate2"); break;
		
		case EXP_RATE_NM: 	
			{
				auto rate = eqn[dp[0].eq_ref].calculate_indfac(ind,ti,popnum,precalc);
				auto ts = t_begin+exp_rate_sample(rate*dt,warn);
				trig.type = NM_TRANS_SIM_EV; trig.i = i; trig.c = UNSET; trig.trg = tgl; trig.tdiv = ts;
			}
			break;
			
		case EXP_MEAN_NM: 	
			{
				auto mean = eqn[dp[0].eq_ref].calculate_indfac(ind,ti,popnum,precalc);
				auto ts = t_begin+exp_mean_sample(mean/dt,warn);
				trig.type = NM_TRANS_SIM_EV; trig.i = i; trig.c = UNSET; trig.trg = tgl; trig.tdiv = ts;
			}
			break;

		case GAMMA:
			{
				auto mean = eqn[dp[0].eq_ref].calculate_indfac(ind,ti,popnum,precalc);
				auto cv = eqn[dp[1].eq_ref].calculate_indfac(ind,ti,popnum,precalc);
			
				auto ts = t_begin+gamma_sample(mean/dt,cv,warn);
				
				trig.type = NM_TRANS_SIM_EV; trig.i = i; trig.c = UNSET; trig.trg = tgl; trig.tdiv = ts;
			}
			break;
		
		case ERLANG:
			{
				auto mean = eqn[dp[0].eq_ref].calculate_indfac(ind,ti,popnum,precalc);
				auto shape = eqn[dp[1].eq_ref].calculate_indfac(ind,ti,popnum,precalc);
				auto ts = t_begin+gamma_sample(mean/dt,sqrt(1.0/shape),warn);
					
				trig.type = NM_TRANS_SIM_EV; trig.i = i; trig.c = UNSET; trig.trg = tgl; trig.tdiv = ts;
			}
			break;

		case LOG_NORMAL:
			{
				auto mean = eqn[dp[0].eq_ref].calculate_indfac(ind,ti,popnum,precalc);
				auto cv = eqn[dp[1].eq_ref].calculate_indfac(ind,ti,popnum,precalc);
				auto ts = t_begin+lognormal_sample(mean/dt,cv,warn);

				trig.type = NM_TRANS_SIM_EV; trig.i = i; trig.c = UNSET; trig.trg = tgl; trig.tdiv = ts;
			}
			break;
			
		case WEIBULL:
			{
				auto scale = eqn[dp[0].eq_ref].calculate_indfac(ind,ti,popnum,precalc);
				auto shape = eqn[dp[1].eq_ref].calculate_indfac(ind,ti,popnum,precalc);
				auto ts = t_begin+weibull_sample(scale/dt,shape,warn);
		
				trig.type = NM_TRANS_SIM_EV; trig.i = i; trig.c = UNSET; trig.trg = tgl; trig.tdiv = ts;
			}
			break;
			
		case PERIOD:
			{
				auto time = eqn[dp[0].eq_ref].calculate_indfac(ind,ti,popnum,precalc);
				auto ts = t_begin+period_sample(time/dt,warn);
	
				trig.type = NM_TRANS_SIM_EV; trig.i = i; trig.c = UNSET; trig.trg = tgl; trig.tdiv = ts;
			}
			break;
		}	
	
		if(warn != "") sp.sampling_error(tgl,warn);
		
		loop++; 
		if(loop > SIM_TRY_MAX){
			stringstream ss; 
			ss << "Could not simulate from the system (transition '" << tr.name << "' on individual '" << ind.name << "' from time " << t_begin << " did not exceed " << t << " after " << SIM_TRY_MAX << " tries).";
			emsg(ss.str());
		}
	}while(trig.tdiv < t);
	
	return trig;
}


/// Tries to insert an event from the data into the time line 
bool StateSpecies::try_insert_data_trans_event(double t, unsigned int cl, unsigned int i)
{
	const auto &ind = individual[i];
	const auto &co = sp.comp_gl[ind_sim_c[i]];
	const auto &tlg = co.tra_leave_group[cl];
	const auto &tref = tlg.tr_list;
	
	if(i >= sp.nindividual_obs) return false;
	
	auto &obs = sp.individual[i].obs;

	auto m = 0u; while(m < obs.size() && obs[m].tdiv <= t) m++;
	
	while(m < obs.size() && !((obs[m].type == OBS_TRANS_EV || obs[m].type == OBS_SINK_EV) && obs[m].cl == cl)) m++;

	if(m == obs.size()) return false;

	if(ind.ev.size() == 0) emsg("e origin prob");

	auto t_event = obs[m].tdiv;
	
	vector <SimTrigEvent> trig_pos;
	
	vector <double> val_sum_store;
	
	auto val_sum = 0.0;
	for(auto k = 0u; k < tref.size(); k++){
		auto tgl = tref[k];
		auto val = get_trans_obs_prob(tgl,obs[m]);
		if(val > 0){
			val_sum += val;
			val_sum_store.push_back(val_sum);
						
			SimTrigEvent trig; 
			trig.type = NM_TRANS_SIM_EV; trig.i = i; trig.c = UNSET; trig.trg = tgl; trig.tdiv = t_event;
			trig_pos.push_back(trig);
		}
	}
	
	if(trig_pos.size() == 0) return false;
	
	if(trig_pos.size() != val_sum_store.size()) emsg("Wrong size");
		
	auto z = val_sum*ran();
	auto k = 0u; while(k < val_sum_store.size() && z > val_sum_store[k]) k++;

	const auto &trig = trig_pos[k];

	insert_trigger_event(trig);

	return true;
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
void StateSpecies::insert_trigger_event(const SimTrigEvent &trig)
{
	auto t = trig.tdiv;
	
	if(t < details.T){
		auto ti = get_ti(trig.tdiv);
		if(ti >= trig_div.size()) emsg("nm out of range8");
		
		auto &ev = trig_div[ti].ev;
		
		if(ti_sort == ti){  // If adding to present time segment then must do it in an ordered way
			auto k = 0u; while(k < ev.size() && t > ev[k].tdiv) k++;
			ev.insert(ev.begin()+k,trig);
		}
		else{               // If added to a later time segment ordering is done later
			ev.push_back(trig);
		}
	}
}


/// Calculates indfac_sum on each Markov equation
void StateSpecies::calculate_indfac_sum()
{
	for(auto e = 0u; e < sp.markov_eqn.size(); e++){
		const auto &me = sp.markov_eqn[e];
		if(!me.source){
			markov_eqn_vari[e].indfac_sum = 0;
		}
	}
	
	for(auto i = 0u; i < individual.size(); i++){
		auto c = ind_sim_c[i];
		if(c != UNSET){
			
			for(const auto &tlg : sp.comp_gl[c].tra_leave_group){
				if(tlg.markov){
					for(auto tr : tlg.tr_list){
						auto e = sp.tra_gl[tr].markov_eqn_ref;
						const auto &me = sp.markov_eqn[e];
						markov_eqn_vari[e].indfac_sum += get_indfac(individual[i],me);
					}
				}
			}
		}
	}
}


/// Updates the value of the markov equations
void StateSpecies::update_markov_eqn_value(unsigned int ti, const vector < vector <double> > &popnum_t)
{
	const auto &popnum = popnum_t[ti];
	const auto &precalc = param_val.precalc;
	
	for(auto e = 0u; e < sp.markov_eqn.size(); e++){
		auto &me = sp.markov_eqn[e];
		auto &me_vari = markov_eqn_vari[e];
		
		const auto &eq = eqn[me.eqn_ref]; 
		double value;
		if(me.rate){		
			value = eq.calculate(ti,popnum,precalc);
			if(value < -TINY){
				run_error("The transition rate determined by equation '"+eq.te_raw+"' has become negative");
			}
		}
		else{
			auto mean = eq.calculate(ti,popnum,precalc);
			if(mean <= 0){
				if(mean < 0) run_error("The transition mean determined by equation '"+eq.te_raw+"' has become negative");
				if(mean == 0) run_error("The transition mean determined by equation '"+eq.te_raw+"' has become zero");
			}
			value = 1.0/mean;
		}

		if(me_vari.value != value){
			update_markov_tree_rate(e,(value-me_vari.value)*me_vari.indfac_sum);
			me_vari.value = value;
		}
	}
}


/// Used to order trigger events
bool Sim_Trig_ord(const SimTrigEvent &ev1, const SimTrigEvent &ev2)                      
{ 
	if(ev1.tdiv == ev2.tdiv){
		if(ev1.i == ev2.i) return ev1.i < ev2.i;
		return ev1.type < ev2.type;
	}
	return ev1.tdiv < ev2.tdiv; 
}


/// Sorts the triger events in the time division
void StateSpecies::sort_trig_event(unsigned int ti)
{
	ti_sort = ti;
	
	auto &tri_ev = trig_div[ti].ev;
	sort(tri_ev.begin(),tri_ev.end(),Sim_Trig_ord);
}


/// Samples the infecting individual from those individuals which are infected
void StateSpecies::sample_infecting_ind(unsigned int i, double t, unsigned int tr_gl, IndInfFrom &inf_from, const vector < vector <double> > &popnum_t, const vector < vector <Poss> > &pop_ind) 
{
	const auto &tra = sp.tra_gl[tr_gl];
	if(tra.infection.type == TRANS_INFECTION){ 
		if(mode == INF && false){ inf_from = get_waifw(i,t); return;}
		
		const auto &me = sp.markov_eqn[tra.markov_eqn_ref];
		const auto &eq = eqn[me.eqn_ref];
		const auto &lin = eq.linearise;
		
		if(!lin.on) emsg("Linearisation should be on");
		
		// First selects if infection comes from outside or from which population
		
		auto Npop = eq.pop_ref.size();
		
		auto j = UNSET;
		
		if(lin.multi_source){ // Samples from available sources (either populations of fro outside)		
			auto ti = get_ti(t);
			auto ss = eq.setup_source_sampler(ti,popnum_t[ti],param_val);
	
			j = ss.sample_inf_source();
		
			if(slow_check) prob_trans_tree += ss.prob_inf_source(j);
		}
		else{
			if(testing){
				if(Npop > 1) emsg("Npop wrong");
				if(Npop == 1 && lin.no_pop_precalc.type != ZERO) emsg("Npop wrong");
			}
			j = 0;
		}
		
		if(j == Npop){
			if(false){
				cout << "OUTSIDE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
			}
			inf_from.p = OUTSIDE_INF; 
			inf_from.i = UNSET;
			inf_from.pref = UNSET;
			inf_from.po = UNSET;
			inf_from.w = 1;
		}
		else{
			auto pr = eq.pop_ref[j];
			inf_from.p = pop[pr].sp_p;
			
			auto pos = sample_possibility(pop_ind[pr]);		
			inf_from.i = pos.i;
			inf_from.pref = j;
			inf_from.po = pr;
			inf_from.w = pos.weight;
			
			if(slow_check) prob_trans_tree += sample_probability(inf_from.i,pop_ind[pr]);
		}
	}
}


/// Gets waifw from a file created during simulation (for diagnostic purposes)
IndInfFrom StateSpecies::get_waifw(unsigned int i, double t) const 
{	
	ifstream fin("waifw.txt");

	vector < vector <double> > row;
	do{
		string line;
		getline(fin,line);
		
		if(fin.eof()) break;
		
		auto spl = split(line,',');
		vector <double> num;
		for(auto i = 0u; i < spl.size(); i++){
			num.push_back(number(spl[i]));
		}
		row.push_back(num);
	}while(true);
	
	auto j = 0u;
	while(j < row.size() && row[j][0] != i) j++;
	if(j == row.size()) emsg("could not find");
	
	IndInfFrom iif;
	iif.p = row[j][1];
	iif.i = row[j][2];
	iif.pref = row[j][3];
	
	auto dd = row[j][4] - t;
	if(dd*dd > 0.0001){
		cout << row[j][4] << " " << t << "compare time" << endl;
		emsg("error with time");
	}
	
	return iif;
}


/// Updates an individual with a transition
void StateSpecies::update_ind_trans(double t, const IndTransRef &itr, const vector < vector <double> > &popnum_t, const vector < vector <Poss> > &pop_ind)
{
	auto i = itr.i;
	
	if(ind_sim_c[i] == UNSET && individual[i].ev.size() != 0) return;

	auto tr_gl = itr.tr_gl;
	
	// If transmission tree then samples infecting individual 
	IndInfFrom inf_from;
	
	if(sp.trans_tree == true) sample_infecting_ind(i,t,tr_gl,inf_from,popnum_t,pop_ind);
	
	// Updates the individual compartment based on transition
	if(ind_sim_c[i] != sp.tra_gl[tr_gl].i){ // Case when intermediate event has occured
		tr_gl = sp.tr_trans(tr_gl,ind_sim_c[i]);
	}
	
	const auto &tr = sp.tra_gl[tr_gl];	
	
	auto type = NM_TRANS_EV; if(tr.markov_eqn_ref != UNSET) type = M_TRANS_EV;
	
	if(ind_sim_c[i] != UNSET) cpop[ind_sim_c[i]]--;
	
	add_event(type,i,tr_gl,UNSET,UNSET,tr.f,t,inf_from);
	
	remove_markov_trans(i);
	
	// Updates individual based on new compartment
	if(ind_sim_c[i] != UNSET){
		auto cl = tr.cl; if(tr.variety == SOURCE_TRANS) cl = UNSET;
		update_ind_c(i,t,cl,popnum_t);
	}
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
	if(ind_sim_c[i] != UNSET){
		cpop[ind_sim_c[i]]--;
	
		add_event(LEAVE_EV,i,UNSET,UNSET,UNSET,UNSET,t,IndInfFrom());
	
		remove_markov_trans(i);
	}
}


/// Moves an individual to a particular compartment in a given classification
void StateSpecies::update_ind_move(unsigned int i, double t, unsigned int c_comp, unsigned int cl, const vector < vector <double> > &popnum_t)
{
	int c = ind_sim_c[i];
	if(c == UNSET) return;
	
	cpop[c]--;
	
	auto cnew = sp.update_c_comp(c,cl,c_comp);
	
	IndInfFrom iif;
	if(sp.comp_gl[c].infected == false && sp.comp_gl[cnew].infected == true){
		iif.p = ENTER_INF;
	}
	
	add_event(MOVE_EV,i,UNSET,c_comp,cl,cnew,t,iif);	

	remove_markov_trans(i);
	
	update_ind_c(i,t,cl,popnum_t);
}
		
				
Event StateSpecies::get_event(EventType type, unsigned int i, unsigned int tr_gl, unsigned int move_c, unsigned int cl, unsigned int c_after, double t, const IndInfFrom &inf_from)
{
	if(tr_gl != UNSET) cl = sp.tra_gl[tr_gl].cl;
	
	Event e; 
	e.type = type;
	e.tr_gl = tr_gl;
	e.move_c = move_c;
	e.cl = cl;
	e.c_after = c_after;
	e.tdiv = t;
	e.inf_node_ref = UNSET;
	e.ind_inf_from = inf_from;
	
	e.observed = false;
	if(mode == INF && i < sp.nindividual_in){
		for(const auto &ob : sp.individual[i].obs){
			switch(ob.type){
			case OBS_TRANS_EV: case OBS_SOURCE_EV: case OBS_SINK_EV: 
				if(cl == ob.cl && ob.tdiv == t) e.observed = true;
				break;
			default: break;
			}
		}			
	}
	
	return e;
}
				
/// Adds event individual timeline
void StateSpecies::add_event(EventType type, unsigned int i, unsigned int tr_gl, unsigned int move_c, unsigned int cl, unsigned int c_after, double t, const IndInfFrom &inf_from)
{
	auto &ind = individual[i];
	auto &ev = ind.ev;
	
	auto e = get_event(type,i,tr_gl,move_c,cl,c_after,t,inf_from);
	
	ev.push_back(e);

	if(mode == INF){
		if(type == M_TRANS_EV || type == NM_TRANS_EV){
			auto iobs = i;
			if(iobs >= sp.nindividual_obs) iobs = sp.nindividual_obs;
			const auto &obs_data = sp.individual[iobs].obs;
			
			double temp;
			alter_obs_trans_eqn(tr_gl,t,1,e.observed,obs_data,temp);
		}
	}
	
	ind_sim_c[i] = c_after;
}


/// Adds an individual to the state
unsigned int StateSpecies::add_individual(IndType ind_type, string name)
{
	if(ind_type == UNOBSERVED_IND){
		if(name != "") emsg("Name should not be set");
		
		auto i = individual.size();
		
		name = sp.name+"-"+unobserved_str;
		if(mode == SIM) name = "Sim-"+sp.name+"-Ind-"+to_string(i+1);
		else{
			if(sp.trans_tree){
				if(unused_name.size() > 0){
					name = unused_name[unused_name.size()-1];
					unused_name.pop_back();
				}
				else{
					name = "Unobs-"+sp.name+"-Ind-"+to_string(i+1);
				}
			}
		}
	}
	
	Individual ind;
	ind.incomp_ref.resize(sp.ncla);
	ind.type = ind_type;
	ind.name = name;
	ind.ie = sample_ie();
	
	set_exp_ie(ind);
	
	auto i = individual.size();
	
	if(i < sp.nindividual_in) ind.init_c_set = sp.individual[i].init_c_set;
	else ind.init_c_set = false;

	switch(type){
	case UNOBSERVED_IND:
		if(sp.fix_effect.size() > 0) run_error("Cannot have fixed affects with unobserved individuals");
		break;
		
	case OBSERVED_IND:
		for(auto f = 0u; f < sp.fix_effect.size(); f++){
			const auto &fe = sp.fix_effect[f];
			const auto &X = fe.X_vector;
			
			if(i >= X.ind_list.size()) run_error("Individual '"+name+"' does not have a design matrix X specified for fixed effect '"+fe.name+"'");
			
			if(X.ind_list[i] != name) emsg("Names do not agree");
			ind.X.push_back(X.value[i]);
			ind.exp_fe.push_back(UNSET);
		}
		break;
	}
	
	individual.push_back(ind);
	
	ind_sim_c.push_back(UNSET);

	if(individual.size() > details.individual_max){
		run_error("Maximum number of "+tstr(details.individual_max)+" individuals reached");
	}
	
	return i;
}


/// Removes an unused individual
void StateSpecies::remove_individual(unsigned int i, vector <InfNode> &inf_node)
{
	auto &ind = individual[i];
	
	unused_name.push_back(ind.name);
	
	if(ind.ev.size() != 0){
		emsg("removed ind has events");
	}
	
	auto imax = individual.size()-1;
	if(i < imax){
		ind = individual[imax];
		for(const auto &ev : ind.ev){
			switch(ev.type){
			case M_TRANS_EV:
				markov_eqn_vari[ev.m].div[ev.ti].ind_trans[ev.index].i = i;
				
				if(sp.trans_tree){
					auto n = ev.inf_node_ref;
					if(n != UNSET){
						inf_node[n].i = i;
					}						
				}
				break;
				
			case NM_TRANS_EV:
				nm_trans_ev_ref[ev.m][ev.ti][ev.index].i = i;
				break;
				
			default: break;
			}
		}
		
		if(sp.nm_flag){
			for(const auto &ir : ind.incomp_ref){
				if(ir.on){
					nm_trans_incomp_ref[ir.n][ir.ti][ir.index] = i; 
				}
			}
		}
		
		individual[imax].ev.resize(0);
	}
	
	individual.pop_back();
}


/// Calculates the number of unobserved entering
void StateSpecies::calculate_N_unobs()
{
	const auto &ic = sp.init_cond;
	if(ic.type != INIT_POP_DIST) return;
	auto &icv = init_cond_val;

	auto foc_cl = ic.focal_cl;	
	if(foc_cl == UNSET){
		icv.N_total_unobs = 0;
	}
	else{
		auto C = sp.cla[foc_cl].ncomp;
		icv.N_focal_unobs.resize(C);
		for(auto c = 0u; c < C; c++) icv.N_focal_unobs[c] = 0;
	}
	
	for(auto i = sp.nindividual_obs; i < individual.size(); i++){
		const auto &ind = individual[i];
		if(ind.ev.size() > 0){
			if(ind.ev[0].type == ENTER_EV && ind.ev[0].tdiv == 0){
				if(foc_cl == UNSET) icv.N_total_unobs++;
				else{
					auto c = ind.ev[0].c_after;
					icv.N_focal_unobs[sp.comp_gl[c].cla_comp[foc_cl]]++;
				}
			}
		}
	}
}

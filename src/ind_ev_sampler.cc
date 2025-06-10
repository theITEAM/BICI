// This deals with individual-based updates where there are individual-based observations

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
using namespace std;

#include "ind_ev_sampler.hh"
#include "utils.hh"
//#include "matrix.hh"

IndEvSampler::IndEvSampler(vector <MarkovEqnVariation> &markov_eqn_vari, const vector <Individual> &individual, const Details &details, const Species &sp, vector <double> &obs_eqn_value, vector < vector <double> > &obs_trans_eqn_value, const vector <Equation> &eqn, const vector <InfNode> &inf_node, const vector <double> &param_val, const vector <SplineValue> &spline_val, vector < vector <double> > &popnum_t) : markov_eqn_vari(markov_eqn_vari), individual(individual), details(details), sp(sp), obs_eqn_value(obs_eqn_value), obs_trans_eqn_value(obs_trans_eqn_value), eqn(eqn), inf_node(inf_node), param_val(param_val),spline_val(spline_val),popnum_t(popnum_t),nm_trans(sp.nm_trans)
{
	T = sp.T;
	
	auto compmax = 0u;
	for(auto cl = 0u; cl < sp.ncla; cl++){
		if(sp.cla[cl].comp.size() > compmax) compmax = sp.cla[cl].comp.size();
	}
		
	ind_obs_prob.resize(T+1);
	for(auto ti = 0u; ti < T+1; ti++){
		ind_obs_prob[ti].resize(compmax);
	}
	
	c_timeline.resize(T);
	
	ind_variation = sp.ind_variation;
}

	
/// Determines if the sampler is needed
bool IndEvSampler::needed(unsigned int i, unsigned int cl)
{	
	const auto &ind = individual[i];
	const	auto &ev = ind.ev;
	
	if(sp.period_exist && contain_period(ev,cl)) return false;
	
	const auto &claa = sp.cla[cl];
	
	illegal = false;
	
	const auto &e = ind.ev[0];
	t_start = e.t;

	ti_start = get_ti_upper(t_start);
	auto c_start = e.c_after;
	
	auto ci = sp.comp_gl[c_start].cla_comp[cl];
	c_start_raw = c_start - ci*sp.comp_mult[cl];

	isl = claa.comp[ci].island_ref.i;
	if(isl == UNSET){
		cout << "island not set" << endl;
		return false;
	}
	
	const auto &comp = claa.island[isl].comp;
	C = comp.size();
	if(C == 1) return false;
	
	i_store = i; cl_store = cl;
	
	auto &obs = sp.individual[i].obs;
	nobs = obs.size();
	if(nobs > 0) ti_end	= 1+get_ti_lower(obs[nobs-1].t); 
	else ti_end = ti_start;
	
	if(sp.obs_trans_exist){ // If transitions are observed, extend to cover region
		auto ti_tobs = sp.last_obs_trans_ti;
		if(ti_tobs != UNSET && ti_tobs > ti_end){
			ti_end = ti_tobs;
			const auto &e_last = ind.ev[ind.ev.size()-1];
			if(e_last.c_after == UNSET){
				auto ti_leave = 1+get_ti_lower(e_last.t); 
				if(ti_end > ti_leave) ti_end = ti_leave;
			}
		}
	}
	
	// Works out how c changes over time within the sampler
	
	auto c = c_start;
	auto ti = ti_start;
	for(auto k = 1u; k < ev.size(); k++){
		const auto &e = ev[k];
		auto ti_next = get_ti_lower(e.t);
		
		while(ti <= ti_next){ c_timeline[ti] = c; ti++;}
		
		auto c_next = e.c_after;
		if(c_next != UNSET) c = c_next;
	}
	
	while(ti < ti_end){ c_timeline[ti] = c; ti++;}
	
	return true;
}
			
			
/// Sets up rate_store (which stores the rate of events) 
void IndEvSampler::rate_store_initialise()
{
	auto lmax = 0u, Cmax = 0u;
	for(const auto &claa : sp.cla){ 
		for(const auto &island : claa.island){
			if(island.comp.size() > Cmax) Cmax = island.comp.size();
			for(const auto &comp : island.comp){				
				if(comp.leave.size() > lmax) lmax = comp.leave.size();
			}
		}	
	}
	
	rate_store.resize(T);
	for(auto ti = 0u; ti < T; ti++){
		rate_store[ti].resize(Cmax);
		for(auto c = 0u; c < Cmax; c++){
			rate_store[ti][c].resize(lmax);
		}
	}		
}


/// Gets the observation probability on an observation for a transition
double IndEvSampler::get_trans_obs_prob(unsigned int trg, const ObsData &ob) const
{
	if(ob.time_vari) return obs_eqn_value[ob.obs_eqn_ref[trg]];
	else return obs_eqn_value[sp.obs_trans[ob.ref].obs_eqn_ref[trg]];
}


/// Generates a timeline giving the observation probability of an individual being in a state 
void IndEvSampler::generate_ind_obs_timeline()
{
	if(nm_trans.size() != nm_rate.size()) emsg("nm_rate not set");
	
	auto dt = details.dt;
	
	const auto &ind = individual[i_store];
	
	const auto &claa = sp.cla[cl_store];
	
	const auto &island = claa.island[isl];
	
	const auto &comp = island.comp;
	
	auto &tp = sp.timepoint;
	
	auto &obs = sp.individual[i_store].obs;
	int oi = nobs-1;
	
	vector <double> inf_other;
	if(sp.trans_tree && cl_store == sp.infection_cl){   // Incorporates transmission tree
		for(auto &ev : ind.ev){
			auto n = ev.inf_node_ref;
			if(n != UNSET && n != inf_node.size()){
				for(const auto &iev : inf_node[n].inf_ev){
					inf_other.push_back(iev.t);
				}
			}
		}
		
		sort(inf_other.begin(),inf_other.end());
	}
	
	int gi = inf_other.size()-1;
	
	// Works out the start and end of the sampler

	// Backpropogates the observation model
	
	fixed_event_obs.clear();
	
	for(auto c = 0u; c < C; c++) ind_obs_prob[ti_end][c] = 1;
	
	for(int tii = int(ti_end)-1; tii >= int(ti_start); tii--){
		auto ctime = c_timeline[tii];
		
		auto &iop = ind_obs_prob[tii];
		const auto &iop_next = ind_obs_prob[tii+1];
		
		for(auto c = 0u; c < C; c++){
			const auto &co = comp[c];
			
			auto sum = 0.0;

			auto lmax = co.leave.size();
		
			if(lmax == 0){
				iop[c] = iop_next[c];
			}
			else{
				vector <double> num(lmax);
				auto &rs = rate_store[tii][c];
				double rate;
				
				for(auto l = 0u; l < lmax; l++){
					const auto &le = co.leave[l];
				
					rate = calculate_rate(ind,le,ctime,tii);
			
					rs[l]	= rate;
					rate *= dt;
					num[l] = rate;
					sum += rate;
				}
			
				if(sum > IND_SAMP_THRESH){
					auto f = IND_SAMP_THRESH/sum;
					sum *= f;
					for(auto l = 0u; l < lmax; l++) num[l] *= f;
				}
				
				iop[c] = (1-sum)*iop_next[c];
				
				for(auto l = 0u; l < lmax; l++){
					const auto &le = co.leave[l];
				
					auto val = num[l]*iop_next[le.cf];
					
					if(sp.obs_trans_exist){  // Accounts for observation probability
						auto tr_gl = sp.trg_from_tr[c][cl_store][le.tr];		
						for(auto m : sp.obs_trans_eqn_ref[tr_gl][tii]){
							auto va = obs_trans_eqn_value[m][tii];
							if(va == LI_WRONG){ val = 0; break;}
							
							val *= exp(va);
						}
					}
				
					iop[c] += val;
				}
			}
		}
			
		auto normalise = false;
		
		// Accounts for observation model	
		while(oi >= 0 && obs[oi].t > tp[tii]){	
			const auto &ob = obs[oi];
			switch(ob.type){
			case OBS_SOURCE_EV: 
				break;
				
			case OBS_TRANS_EV: case OBS_SINK_EV: 
				if(ob.cl == cl_store){	
					switch(ob.type){
					case OBS_TRANS_EV:
						{
							vector <double> iop_new(C,0);
							for(auto c = 0u; c < C; c++){
								auto val = 0.0;
								for(const auto &le : comp[c].leave){
									auto trg = sp.trg_from_tr[ctime][cl_store][le.tr];
									val += get_trans_obs_prob(trg,ob)*iop[le.cf];
								}
								
								if(val < OBS_COMP_MIN) val = OBS_COMP_MIN;
								iop_new[c] = val; iop[c] *= val;
							}
							iop = iop_new;
						}
						break;
						
					case OBS_SINK_EV:
						{
							double val;
							for(auto c = 0u; c < C; c++){
								auto cgl = sp.get_comp_global_convert(ctime,cl_store,comp[c].c);
								auto trg = sp.cgl_tr_sink[cgl];
								if(trg == UNSET) val = OBS_COMP_MIN;
								else{
									val = get_trans_obs_prob(trg,ob);
									if(val < OBS_COMP_MIN) val = OBS_COMP_MIN;
								}
								
								iop[c] *= val;
							}
						}
						break;
						
					default: break;
					}
					
					// This check to see if there are multiple observations of the same event
					auto fl = false;
					if(fixed_event_obs.size() > 0){
						const auto &obs_last = obs[fixed_event_obs[fixed_event_obs.size()-1]];
						if(obs_last.t == obs[oi].t) fl = true;
					}
			
					if(fl == false) fixed_event_obs.push_back(oi);
				}
				break;
			
			case OBS_COMP_EV:
				if(ob.cl == cl_store){
					auto c_exact = ob.c_exact;
					if(c_exact != UNSET){
						for(auto j = 0u; j < C; j++){
							if(comp[j].c != c_exact) iop[j] *= OBS_COMP_MIN;
						}
					}
					else{
						for(auto j = 0u; j < C; j++){
							auto val = obs_eqn_value[ob.obs_eqn_ref[comp[j].c]];
							if(val < OBS_COMP_MIN) val = OBS_COMP_MIN;
							iop[j] *= val;
						}
					}
				}
				break;
				
			case OBS_TEST_EV:
				if(ob.cl == cl_store){
					const auto &om = sp.source[ob.so].obs_model;
					auto Se = obs_eqn_value[ob.Se_obs_eqn_ref];
					auto Sp = obs_eqn_value[ob.Sp_obs_eqn_ref];
				
					for(auto c = 0u; c < C; c++){
						const auto &co = comp[c];
						
						double val;
						if(om.diag_test_sens.comp[co.c] == true){ // Truely infected
							if(ob.test_res == true) val = Se; 
							else val = 1-Se;
						}
						else{
							if(ob.test_res == false) val = Sp; 
							else val = 1-Sp;
						}
						if(val < OBS_COMP_MIN) val = OBS_COMP_MIN;
						iop[c] *= val; 
					}					
				}
				break;
			}
			oi--;
			
			normalise = true;
		}
		
		while(gi >= 0 && inf_other[gi] > tp[tii]){
			for(auto j = 0u; j < C; j++){
				if(claa.comp[comp[j].c].infected == COMP_UNINFECTED){
					iop[j] = 0;
				}
			}
			gi--;
			
			normalise = true;
		}
		
		if(normalise == true){                                   // Normalises to avoid going to zero
			auto sum = 0.0; 
			for(auto c = 0u; c < C; c++) sum += iop[c];
			if(sum == 0) emsg("No observation probability");
			for(auto c = 0u; c < C; c++) iop[c] /= sum;
		}
	}

	if(false){
		cout << ti_start << " " << ti_end << " time range" << endl;
		for(auto ti = ti_start; ti <= ti_end; ti++){
			cout << ti << " " << c_timeline[ti] << "  obs prob " << endl;
			
			for(auto c = 0u; c < C; c++){ 
				const auto &co = comp[c];
			
				cout << claa.comp[co.c].name << ": prob=";
				cout << ind_obs_prob[ti][c];
				
				if(false){
					cout << ",  Leave: ";
					auto lmax = co.leave.size();
					for(auto l = 0u; l < lmax; l++){
						const auto &le = co.leave[l];
						
						auto ctime = c_timeline[ti];
		
						const auto &mer = le.markov_eqn_ref[ti][ctime];
						
						cout << claa.tra[le.tr].name << " " << eqn[sp.markov_eqn[mer.e].eqn_ref].te_raw << "," << mer.ti << " "<< markov_eqn_vari[mer.e].div[mer.ti].value << "  |   ";
					}
					cout << endl;
				}
				cout << "  ";
			}
			cout << "   time line" << endl;
		}
	
		if(isl == 1) emsg("Pp");
	}
}


/// Samples a new set of events based on observation probability
vector <Event> IndEvSampler::sample_events(double &probif)
{
	auto pif = 1.0;
	
	const auto &ind = individual[i_store];
	
	const auto &ev = ind.ev;
	
	const auto &claa = sp.cla[cl_store];
	
	const auto &island = claa.island[isl];
	
	const auto &comp = island.comp;
	
	auto &tp = sp.timepoint;
	
	auto t = t_start;
	auto ti = ti_start;
	
	auto e = ev[0];
	e.inf_node_ref = UNSET;
	e.ind_inf_from = IndInfFrom();
	
	auto c = e.c_after;
	if(!ind.init_c_set){
		auto sampler = init_state_sampler_obs(i_store,ti,e);
	
		if(sampler.size() > 0){
			auto z = ran();
			auto j = 0u; while(j < C && z > sampler[j].prob_sum) j++;
			if(j == C) emsg("should not be C");
			c = c_start_raw + claa.island[isl].comp[j].c*sp.comp_mult[cl_store];

			pif *= sampler[j].prob;
		}
	}
	
	ev_new.clear();
	switch(e.type){
	case ENTER_EV:
		e.c_after = c;
		break;
	
	case M_TRANS_EV:
		e.c_after = c;
		e.tr_gl = sp.cgl_tr_source[c];
		break;
		
	default:
		emsg("SOURCE EVENT");
		break;
	}
	ev_new.push_back(e);

	auto ci = sp.comp_gl[c].cla_comp[cl_store]; 
	auto cisland = claa.comp[ci].island_ref.c;
	
	auto &obs = sp.individual[i_store].obs;
	
	double time_fixed;

	auto feo_max = fixed_event_obs.size();
	int feo = int(feo_max)-1;
	if(feo >= 0) time_fixed = obs[fixed_event_obs[feo]].t; else time_fixed = LARGE;
	
	// Goes through events on the existing event sequence
	for(auto k = 1u; k <= ev.size(); k++){
		double t_next;
		
		// Events in the existing sequence for the selected classification are not transfered
		auto ignore = false; 
		
		if(k < ev.size()){
			const auto &e = ev[k];
			
			t_next = e.t;
			if(e.type == M_TRANS_EV || e.type == NM_TRANS_EV){
				const auto &tra = sp.tra_gl[e.tr_gl];
				if(tra.cl == cl_store && tra.variety != SINK_TRANS) ignore = true;
			}
		}
		else t_next = tp[T];
		
		if(ignore == false){
			// Goes down the time period up to the next event adding in transition events
			while(t < t_next){
				double tt;
				if(tp[ti+1] < t_next) tt = tp[ti+1];
				else tt = t_next;
				if(tt > time_fixed) tt = time_fixed;
				
				auto ddt = tt-t;
			
				const auto &co = comp[cisland];
				auto ctime = c_timeline[ti];
		
				auto lmax = co.leave.size();
				if(lmax > 0){	
					auto tii = ti; if(tii > ti_end) tii = ti_end;
					
					auto &iop = ind_obs_prob[tii];
					
					auto sum = 0.0, sum_op = 0.0;
					
					// Works out transition prob x obs prob for each possible transition
					vector <double> num(lmax), num_op(lmax);
					for(auto l = 0u; l < lmax; l++){
						const auto &le = co.leave[l];
						auto cf = le.cf;
					
						if(iop[cf] == 0) num[l] = 0;
						else{
							if(ti < ti_end) num[l] = ddt*rate_store[tii][cisland][l];
							else num[l] = ddt*calculate_rate(ind,le,ctime,ti);
							
							/*
							if(sp.obs_trans_exist){  // Accounts for observation probability
								auto tr_gl = sp.trg_from_tr[c][cl_store][le.tr];		
								for(auto m : sp.obs_trans_eqn_ref[tr_gl][ti]){
									num[l] *= exp(obs_trans_eqn_value[m][ti]);
								}
							}
							*/
		
							sum += num[l];
						
							num_op[l] = num[l]*iop[cf];
							sum_op += num_op[l];
						}
					}
				
					if(sum_op > 0){     // If any non-zero transitions leaving comprtment
						// Due to finite discretisation time, this ensure algorithm is stable
						if(sum > IND_SAMP_THRESH){
							auto f = IND_SAMP_THRESH/sum;
							sum *= f; sum_op *= f;
							for(auto l = 0u; l < lmax; l++){ num[l] *= f; num_op[l] *= f;}
						}
					
						auto stay_op = (1-sum)*iop[cisland];
						
						auto total = stay_op+sum_op;
						
						if(ran()*total < sum_op){    // An event happens and is added to timeline
							pif *= sum_op/total;
							
							// Selects which event 
							auto sum_op_sum = 0.0;
							vector <double> sum_store(lmax);
							for(auto l = 0u; l < lmax; l++){ 
								sum_op_sum += num_op[l];
								sum_store[l] = sum_op_sum;
							}
							
							auto z = ran()*sum_op_sum;
							auto l = 0u; while(l < lmax && z > sum_store[l]) l++;
							if(l == lmax){ 
								cout << "problem with sampler" << endl; illegal = true; 
								return ev_new;
							}
						
							pif *= num_op[l]/sum_op_sum;
						
							auto tr = co.leave[l].tr;
						
							auto tr_gl = sp.trg_from_tr[c][cl_store][tr];
							const auto &trg = sp.tra_gl[tr_gl];
							if(trg.i != c){ 
								cout << "Not consistent1" << endl; illegal = true; 
								return ev_new;
							}
							if(trg.type == PERIOD){ illegal = true; return ev_new;}
							
							auto t_new = t + ran()*(tt-t);
							pif *= 1.0/(tt-t);
					
							Event enew;
							if(trg.nm_trans_ref != UNSET) enew.type = NM_TRANS_EV;
							else enew.type = M_TRANS_EV;
							enew.move_c = UNSET;
							enew.cl = trg.cl;
							enew.tr_gl = tr_gl;
							enew.t = t_new;
							enew.observed = false;
							enew.inf_node_ref = UNSET;
							enew.ind_inf_from = IndInfFrom();
						
							c = trg.f; 
							ci = claa.tra[tr].f;
							cisland = claa.comp[ci].island_ref.c;
	
							enew.c_after = c;
							
							ev_new.push_back(enew);
						
							tt = t_new;
						}
						else{
							pif *= stay_op/total;
						}
					}
				}
				
				t = tt;
				
				if(t == tp[ti+1]) ti++;
				
				if(t == time_fixed){           // Adds a fixed event to the time line
					auto &ob = obs[fixed_event_obs[feo]];
					
					if(ob.type != OBS_SINK_EV){
						// Samples the fixed event
						const auto &co = comp[cisland];
						
						vector <Poss> poss;
						
						auto val_sum = 0.0;
						for(auto k = 0u; k < co.leave.size(); k++){
							const auto &le = co.leave[k];
						
							auto trg = sp.trg_from_tr[ctime][cl_store][le.tr];
						
							auto va = get_trans_obs_prob(trg,ob);
							if(va > 0){
								Poss po; po.i = k; po.weight = va;
								poss.push_back(po);
								val_sum += va;
							}
						}
						
						unsigned int tr;
						switch(poss.size()){
						case 0: illegal = true; return ev_new;
						case 1: tr = co.leave[poss[0].i].tr; break;
						default:
							{
								auto po = sample_possibility(poss);
								tr = co.leave[po.i].tr;
								pif *= po.weight/val_sum;
							}
							break;
						}
						
						auto tr_gl = sp.trg_from_tr[ctime][cl_store][tr];
						
						ci = claa.tra[tr].f;
						cisland = claa.comp[ci].island_ref.c;
				
						const auto &trg = sp.tra_gl[tr_gl];
						
						if(trg.i != c){
							illegal = true; 
							return ev_new;
						}
						
						Event enew;
						if(trg.nm_trans_ref != UNSET) enew.type = NM_TRANS_EV;
						else enew.type = M_TRANS_EV;
						enew.move_c = UNSET;
						enew.cl = cl_store;
						enew.tr_gl = tr_gl;
						enew.t = ob.t;
						enew.observed = true;
						enew.inf_node_ref = UNSET;
						enew.ind_inf_from = IndInfFrom();						
					
						c = trg.f; 
						
						enew.c_after = c;
						ev_new.push_back(enew);
					}
					feo--;
					if(feo >= 0) time_fixed = obs[fixed_event_obs[feo]].t; else time_fixed = LARGE;
				}
			}
			t = t_next;
			
			if(k < ev.size()){                // Looks to copy events from the existing sequence 
				
				auto enew = ev[k];
				
				if(sp.correct_ev(c,enew) == true){ illegal = true; return ev_new;}
				ev_new.push_back(enew);
			
				c = enew.c_after;
			}		
		}
		if(c == UNSET) break;
	}		
	
	probif = log(pif);
	
	if(sp.period_exist && contain_period(ev_new,cl_store)){
		illegal = true;
		return ev_new;
	}
	
	if(events_near_div(ev_new,details)) illegal = true;
	
	return ev_new;
}


/// Detemines if event sequence either begins or ends a period 
bool IndEvSampler::contain_period(const vector <Event> &event, unsigned int cl) const
{
	for(const auto &ev : event){
		auto c = ev.c_after;
		if(c != UNSET && sp.comp_period[c][cl] == true) return true;
	}
	return false;
}
		
		
/// Creates a sampler for the initial state of the individual
vector <InitStateProb> IndEvSampler::init_state_sampler_obs(unsigned int i, unsigned int ti, const Event &ev) const 
{
	vector <InitStateProb> init_state;
		
	const auto &indi = sp.individual[i];
	const auto &iop_start = ind_obs_prob[ti];
	const auto &claa = sp.cla[cl_store];
	
	auto prob_sum = 0.0;
	if(ev.type == ENTER_EV){  // Samples the original state
		const auto &ent = sp.enter[indi.enter_ref];
		if(ent.c_set != UNSET) return init_state;
		
		const auto &ent_cl = ent.cla[cl_store];
		if(ent_cl.c_set != UNSET) return init_state;
			
		const auto &oer = ent_cl.obs_eqn_ref;
		for(auto j = 0u; j < C; j++){
			auto c = claa.island[isl].comp[j].c;
			//auto cgl = c_start_raw + c*sp.comp_mult[cl_store];
			
			auto prob = obs_eqn_value[oer[c]]*iop_start[j];
			
			if(prob < 0) emsg("Cannot be negative");
			prob_sum += prob;
			
			InitStateProb isp; isp.prob = prob; isp.prob_sum = prob_sum;
			init_state.push_back(isp);
		}
	}
	else{
		for(auto j = 0u; j < C; j++){
			auto c = claa.island[isl].comp[j].c;
			auto cgl = c_start_raw + c*sp.comp_mult[cl_store];
			auto tr = sp.cgl_tr_source[cgl];
			auto prob = 0.0;
			if(tr != UNSET){
				auto &mev = markov_eqn_vari[sp.tra_gl[tr].markov_eqn_ref];
				if(mev.time_vari) prob = mev.div[ti].value;
				else prob = mev.div[0].value;
				if(prob < 0) emsg("Cannot be negative");
				prob_sum += prob;
			}
			InitStateProb isp; isp.prob = prob; isp.prob_sum = prob_sum;
			init_state.push_back(isp);
		}
	}
	
	if(prob_sum == 0.0) emsg("cannot be zero init event samp"); 
	
	for(auto &isp : init_state){
		isp.prob /= prob_sum; isp.prob_sum /= prob_sum;
	}
		
	return init_state;
}


/// Calculates the observation prob of sampling a set of events 
double IndEvSampler::sample_events_prob(const vector <Event> &ev) const
{
	auto pif = 1.0;
	
	const auto &ind = individual[i_store];
	
	const auto &claa = sp.cla[cl_store];
	
	const auto &island = claa.island[isl];
	
	const auto &comp = island.comp;
	
	auto &tp = sp.timepoint;
	
	auto t = t_start;
	auto ti = ti_start;

	const auto &e = ev[0];
	auto c = e.c_after;
	auto ci = sp.comp_gl[c].cla_comp[cl_store]; 

	if(!ind.init_c_set){	
		auto sampler = init_state_sampler_obs(i_store,ti,e);

		if(sampler.size() > 0){
			auto j = 0u; while(j < C && claa.island[isl].comp[j].c != ci) j++;
			if(j == C) emsg("Could not find");
	
			pif *= sampler[j].prob;
		}
	}
	
	auto cisland = claa.comp[ci].island_ref.c;
	
	auto &obs = sp.individual[i_store].obs;
	
	double time_fixed;
	auto feo_max = fixed_event_obs.size();
	int feo = int(feo_max)-1;
	if(feo >= 0) time_fixed = obs[fixed_event_obs[feo]].t; else time_fixed = LARGE;
	
	vector <unsigned int> list;

	// Goes through events on the existing event sequence
	for(auto k = 1u; k <= ev.size(); k++){
		double t_next;
		
		 // Events in the existing sequence for the selected classification are not transfered
		auto ignore = false; 
		
		if(k < ev.size()){
			const auto &e = ev[k];
			
			t_next = e.t;
			if(e.type == M_TRANS_EV || e.type == NM_TRANS_EV){
				const auto &tra = sp.tra_gl[e.tr_gl];
				if(tra.cl == cl_store && tra.variety != SINK_TRANS) ignore = true;
			}
		}
		else t_next = tp[T];
		
		if(ignore == true){
			list.push_back(k);
		}
		else{
			auto t_event = LARGE; if(list.size() > 0) t_event = ev[list[0]].t;
		
			while(t < t_next){
				double tt;
				if(tp[ti+1] < t_next) tt = tp[ti+1];
				else tt = t_next;
		
				if(tt > time_fixed) tt = time_fixed;
			
				auto ddt = tt-t;
				const auto &co = comp[cisland];
				auto ctime = c_timeline[ti];
		
				auto lmax = co.leave.size();
				if(lmax > 0){	
					auto tii = ti; if(tii > ti_end) tii = ti_end;
					
					auto &iop = ind_obs_prob[tii];
					
					auto sum = 0.0, sum_op = 0.0;

					// Works out transition prob x obs prob for each possible transition
					vector <double> num(lmax), num_op(lmax);
					for(auto l = 0u; l < lmax; l++){
						const auto &le = co.leave[l];
						auto cf = le.cf;
					
						if(iop[cf] == 0) num[l] = 0;
						else{
							if(ti < ti_end) num[l] = ddt*rate_store[tii][cisland][l];
							else{
								num[l] = ddt*calculate_rate(ind,le,ctime,ti);
							}
							
							/*
							if(sp.obs_trans_exist){  // Accounts for observation probability
								auto tr_gl = sp.trg_from_tr[c][cl_store][le.tr];		
								for(auto m : sp.obs_trans_eqn_ref[tr_gl][ti]){
									num[l] *= exp(obs_trans_eqn_value[m][ti]);
								}
							}
							*/
					
							sum += num[l];
					
							num_op[l] = num[l]*iop[cf];
							sum_op += num_op[l];
						}
					}
				
					if(sum_op > 0){     // If any non-zero transitions leaving comprtment
						// Due to finite discretisation time, this ensure algorithm is stable
						if(sum > IND_SAMP_THRESH){
							auto f = IND_SAMP_THRESH/sum;
							sum *= f; sum_op *= f;
							for(auto l = 0u; l < lmax; l++){ num[l] *= f; num_op[l] *= f;}
						}
					
						auto stay_op = (1-sum)*iop[cisland];
						
						auto total = stay_op+sum_op;
						
						if(t_event < tt){
							pif *= sum_op/total;
							
							const auto &e = ev[list[0]];
							const auto &trg = sp.tra_gl[e.tr_gl];
							auto tr = trg.tr;
							
							auto l = 0u; while(l < lmax && co.leave[l].tr != tr) l++;
							if(l == lmax) emsg("Problem with sampler");
							
							pif *= num_op[l]/sum_op;
							
							pif *= 1.0/(tt-t);
							c = trg.f; 
							ci = claa.tra[tr].f;
							cisland = claa.comp[ci].island_ref.c;
	
							list.erase(list.begin());
							t_event = LARGE; if(list.size() > 0) t_event = ev[list[0]].t;
							
							tt = e.t;
						}
						else{
							pif *= stay_op/total;
						}
					}
				}
				
				t = tt;
				
				if(t == tp[ti+1]) ti++;
				
				if(t == time_fixed){ 
					auto &ob = obs[fixed_event_obs[feo]];
						
					if(t == t_event){		
						const auto tr_gl = ev[list[0]].tr_gl;
						if(ob.type != OBS_SINK_EV){
							// Gets the probability of selecting this event
							const auto &co = comp[cisland];
							auto val = 0.0, val_sum = 0.0;
							for(auto k = 0u; k < co.leave.size(); k++){
								const auto &le = co.leave[k];
							
								auto trg = sp.trg_from_tr[ctime][cl_store][le.tr];
								
								auto va = get_trans_obs_prob(trg,ob);
								if(trg == tr_gl) val = va;
								val_sum += va;							
							}
						
							if(val == 0) return -LARGE;
							pif *= val/val_sum;
						
							const auto &trg = sp.tra_gl[tr_gl];
							
							auto tr = trg.tr;		
							c = trg.f; 
							ci = claa.tra[tr].f;
							cisland = claa.comp[ci].island_ref.c;
		
							list.erase(list.begin());
							t_event = LARGE; if(list.size() > 0) t_event = ev[list[0]].t;	
						}
						
						feo--;
						if(feo >= 0) time_fixed = obs[fixed_event_obs[feo]].t;
						else time_fixed = LARGE;
					}
					else{
						if(ob.type == OBS_SINK_EV) return log(pif);
						return -LARGE;
					}
				}
			}
			t = t_next;
			
			if(k < ev.size()){
				c = ev[k].c_after;
			}
		}
		if(c == UNSET) break;
	}		
	
	return log(pif);
}


/// Gets ti div time from an actual time (rounding down if equal)
unsigned int IndEvSampler::get_ti_lower(double t) const
{
	return (unsigned int)(ALMOST_ONE*(t-details.t_start)/details.dt);
}


/// Gets ti div time from an actual time (rounding up if equal)
unsigned int IndEvSampler::get_ti_upper(double t) const
{
	return (unsigned int)(OVER_ONE*(t-details.t_start)/details.dt);
}


/// Returns the individual-based factor which acts on Markov equation
double IndEvSampler::get_indfac(const Individual &ind, const MarkovEqn &mar_eqn) const
{
	auto val = 1.0;
	for(auto ie : mar_eqn.ind_eff_mult) val *= ind.exp_ie[ie];
	for(auto fe : mar_eqn.fix_eff_mult) val *= ind.exp_fe[fe];
	if(!mar_eqn.rate) val = 1/val;
	
	return val;
}


/// Sets up variables for dealing with nm transitions
void IndEvSampler::setup_nm()
{
	vector<vector<double> > bp_store;// Note, this is not calculated
	nm_rate = sp.calc_nm_rate(false,param_val,spline_val,popnum_t,eqn,bp_store);
	
	rate_store_initialise();
}


/// Calculates the rate for a transition (for markovian and nm)
double IndEvSampler::calculate_rate(const Individual &ind, const IslandTrans &le, unsigned int ctime, int tii) const
{
	if(le.nm_trans){ // Non markovian transition
		auto m = le.nm_trans_ref[ctime];
	
		if(!nm_trans[m].precalc_nm_rate){
			return calculate_nm_rate(ind,m,tii);
		}
		else{
			auto rate = nm_rate[m][tii];
			if(ind_variation == true){
				//for(const auto &ifr : sp.nm_trans[m].ind_fac_rate){
				for(const auto &ifr : nm_trans[m].ind_fac_rate){
					switch(ifr.type){
					case IND_EFF_MULT: rate *= ind.exp_ie[ifr.e]; break;
					case IND_EFF_DIV: rate /= ind.exp_ie[ifr.e]; break;
					case FIX_EFF_MULT: rate *= ind.exp_fe[ifr.e]; break;
					case FIX_EFF_DIV: rate /= ind.exp_fe[ifr.e]; break;
					}
				}
			}
			return rate;
		}
	}
	else{                   // Markovian transition
		const auto &mer = le.markov_eqn_ref[tii][ctime];

		auto rate = markov_eqn_vari[mer.e].div[mer.ti].value;
		if(ind_variation == true) rate *= get_indfac(ind,sp.markov_eqn[mer.e]);
		return rate;
	}
}

	
/// Calculates the rate of going down a given nm_trans at a given time
double IndEvSampler::calculate_nm_rate(const Individual &ind, unsigned int m, int ti) const
{
	const auto &nmt = nm_trans[m];
	
	auto bp = 1.0;
	
	auto bp_eq = nmt.bp_eq;
	if(bp_eq != UNSET){
		if(bp_eq == BP_FROM_OTHERS){
			for(auto e : nmt.bp_other_eq){
				bp -= eqn[e].calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
			}
		}
		else{
			bp = eqn[bp_eq].calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
			if(nmt.all_branches){
				auto div = 0.0;
				for(auto e : nmt.bp_all_eq){
					div += eqn[e].calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
				}
				bp /= div;
			}
		}
		check_bp(bp);
	}
	
	switch(nmt.type){
	case EXP_RATE: case EXP_MEAN: emsg("Should not be here"); break;

	case GAMMA: case ERLANG: case LOG_NORMAL: case PERIOD: 
		{
			const auto &eq = eqn[nmt.dist_param_eq_ref[0]];
			return bp/eq.calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
		}
		break;

	case EXP_RATE_NM:
		{
			const auto &eq = eqn[nmt.dist_param_eq_ref[0]];
			return bp*eq.calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
		}
		break;
		
	case EXP_MEAN_NM:
		{
			const auto &eq = eqn[nmt.dist_param_eq_ref[0]];
			return bp/eq.calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
		}
		break;
	
	case WEIBULL:
		{
			const auto &eq = eqn[nmt.dist_param_eq_ref[0]];	
			return bp/eq.calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
		}
		break;
		
	default: emsg("trans not recognised"); break;
	}
	
	return UNSET;
}

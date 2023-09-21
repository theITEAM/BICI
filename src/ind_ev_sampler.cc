// This deals with individual-based updates

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

IndEvSampler::IndEvSampler(vector <MarkovEqnVariation> &markov_eqn_vari, const vector <Individual> &individual, const Details &details, const Species &sp, vector <double> &obs_eqn_value, const vector <Equation> &eqn) : markov_eqn_vari(markov_eqn_vari), individual(individual), details(details), sp(sp), obs_eqn_value(obs_eqn_value), eqn(eqn)
{
	T = sp.T;
	
	auto compmax = 0u;
	for(auto cl = 0u; cl < sp.ncla; cl++){
		if(sp.cla[cl].comp.size() > compmax) compmax = sp.cla[cl].comp.size();
	}
		
	ind_obs_prop.resize(T+1);
	for(auto ti = 0u; ti < T+1; ti++){
		ind_obs_prop[ti].resize(compmax);
	}
	
	c_timeline.resize(T);
	
	calc_indfac = false; if(sp.ind_effect.size() > 0) calc_indfac = true;
}
	
	
/// Determines if the sampler is needed
bool IndEvSampler::needed(unsigned int i, unsigned int cl)
{	
	const auto &ind = individual[i];
	const	auto &ev = ind.ev;
	
	const auto &claa = sp.cla[cl];
	
	illegal = false;
	
	ev_new.clear();
	
	ti_start = 0u;
	k_start = 0u;
	c_start = ind.cinit;
	t_start = sp.timepoint[0];
	
	if(c_start > CODE){
		const auto &ev = ind.ev[0];
		t_start = ev.t;
		ev_new.push_back(ev);
		ti_start = get_ti_upper(ev.t);
		c_start = ev.c_after;
		k_start++;
	}
	
	auto ci = sp.comp_gl[c_start].cla_comp[cl];
	isl = claa.comp[ci].island_ref.i;
	if(isl == UNSET){
		cout << "island not set" << endl;
		return false;
	}
	
	C = claa.island[isl].comp.size();
	if(C == 1) return false;
	
	i_store = i; cl_store = cl;
	
	auto &obs = sp.individual[i].obs;
	nobs = obs.size();
	if(nobs > 0) ti_end	= 1+get_ti_lower(obs[nobs-1].t); 
	else ti_end = ti_start;
	
	// Works out how c changes over time within the sampler
	
	auto c = c_start;
	auto ti = ti_start;
	for(auto k = k_start; k < ev.size(); k++){
		const auto &e = ev[k];
		auto ti_next = get_ti_lower(e.t);
		
		while(ti <= ti_next){ c_timeline[ti] = c; ti++;}
		
		c = e.c_after;
	}
	
	while(ti < ti_end){ c_timeline[ti] = c; ti++;}
	
	return true;
}
			
				
/// Generates a timeline giving the observation probability of an individual being in a state 
void IndEvSampler::generate_ind_obs_timeline()
{
	auto dt = details.dt;
	
	const auto &ind = individual[i_store];
	
	const auto &claa = sp.cla[cl_store];
	
	const auto &island = claa.island[isl];
	
	const auto &comp = island.comp;
	
	auto &tp = sp.timepoint;
	
	auto &obs = sp.individual[i_store].obs;
	int oi = nobs-1;
	
	// Works out the start and end of the sampler

	// Backpropogates the observation model
	
	fixed_event_obs.clear();
	
	for(auto c = 0u; c < C; c++) ind_obs_prop[ti_end][c] = 1;
	
	for(int tii = int(ti_end)-1; tii >= int(ti_start); tii--){
		auto ctime = c_timeline[tii];
		
		auto &iop = ind_obs_prop[tii];
		const auto &iop_next = ind_obs_prop[tii+1];
		
		for(auto c = 0u; c < C; c++){
			const auto &co = comp[c];
			
			auto sum = 0.0;

			auto lmax = co.leave.size();
		
			if(lmax == 0){
				iop[c] = iop_next[c];
			}
			else{
				vector <double> num(lmax);
				for(auto l = 0u; l < lmax; l++){
					const auto &le = co.leave[l];
					
					const auto &mer = le.markov_eqn_ref[tii][ctime];
					
					num[l] = dt*markov_eqn_vari[mer.e].div[mer.ti].value;
					if(calc_indfac == true) num[l] *= get_indfac(ind,sp.markov_eqn[mer.e]);
			
					sum += num[l];
				}
			
				if(sum > IND_SAMP_THRESH){
					auto f = IND_SAMP_THRESH/sum;
					sum *= f;
					for(auto l = 0u; l < lmax; l++) num[l] *= f;
				}
				
				iop[c] = (1-sum)*iop_next[c];
				for(auto l = 0u; l < lmax; l++) iop[c] += num[l]*iop_next[co.leave[l].cf];
			}
		}
			
		// Accounts for observation model	
		while(oi >= 0 && obs[oi].t > tp[tii]){
			const auto &ob = obs[oi];
	
			switch(ob.type){
			case OBS_TRANS_EV:
				if(ob.cl == cl_store){
					const auto &tra = claa.tra[ob.tr];
					auto ci = tra.i;
					if(ci != SOURCE){
						auto c_sel = claa.comp[ci].island_ref.c;
						for(auto c = 0u; c < C; c++){
							if(c != c_sel) iop[c] = 0;
						}
	
						fixed_event_obs.push_back(oi);
					}
				}
				break;
			
			case OBS_COMP_EV:
				if(ob.cl == cl_store){
					for(auto c = 0u; c < C; c++){
						if(comp[c].c != ob.c) iop[c] = 0;
					}
				}
				break;
				
			case OBS_TEST_EV:
				if(ob.cl == cl_store){
					const auto &om = sp.source[ob.so].obs_model;
					auto Se = obs_eqn_value[ob.Se_obs_eqn];
					auto Sp = obs_eqn_value[ob.Sp_obs_eqn];
				
					for(auto c = 0u; c < C; c++){
						const auto &co = comp[c];
						
						if(om.diag_test_sens.comp[co.c] == true){ // Truely infected
							if(ob.test_res == true) iop[c] *= Se; 
							else iop[c] *= 1-Se;
						}
						else{
							if(ob.test_res == false) iop[c] *= Sp; 
							else iop[c] *= 1-Sp;
						}
					}					
				}
				break;
			
			case OBS_POP:
				break;
			}
			oi--;
		}
	}

	if(false){
		for(auto ti = ti_start; ti < ti_end; ti++){
			cout << ti << " " << c_timeline[ti] << "  obs prob " << endl;
			auto ctime = c_timeline[ti];
			
			for(auto c = 0u; c < C; c++){ 
				const auto &co = comp[c];
			
				cout << claa.comp[co.c].name << ": prob=";
				cout << ind_obs_prop[ti][c]<< ",  Leave: ";
				
				auto lmax = co.leave.size();
				for(auto l = 0u; l < lmax; l++){
					const auto &le = co.leave[l];
					
					const auto &mer = le.markov_eqn_ref[ti][ctime];
					cout << claa.tra[le.tr].name << " " << eqn[sp.markov_eqn[mer.e].eqn_ref].te_raw << "," << mer.ti << " "<< markov_eqn_vari[mer.e].div[mer.ti].value << "  |   ";
				}
				cout << endl;
			}
			cout << "   time line\n";
		}
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
	
	auto c = c_start;
	auto t = t_start;
	
	auto ci = sp.comp_gl[c].cla_comp[cl_store]; 
	auto cisland = claa.comp[ci].island_ref.c;
	
	auto ti = ti_start;
	
	auto &obs = sp.individual[i_store].obs;
	
	double time_fixed;
	
	auto feo_max = fixed_event_obs.size();
	int feo = int(feo_max)-1;
	if(feo >= 0) time_fixed = obs[fixed_event_obs[feo]].t; else time_fixed = LARGE;
	
	// Goes through events on the existing event sequence
	for(auto k = k_start; k <= ev.size(); k++){
		double t_next;
		
		 // Events in the existing sequence for the selected classification are not transfered
		auto ignore = false; 
		
		if(k < ev.size()){
			const auto &e = ev[k];
			
			t_next = e.t;
			if(e.type == M_TRANS_EV){
				if(sp.tra_gl[e.tr_gl].cl == cl_store) ignore = true;
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
					
					auto &iop = ind_obs_prop[tii];
					
					auto sum = 0.0, sum_op = 0.0;

					// Works out transition prob x obs prob for each possible transition
					vector <double> num(lmax), num_op(lmax);
					for(auto l = 0u; l < lmax; l++){
						const auto &le = co.leave[l];
						auto cf = le.cf;
					
						if(iop[cf] == 0) num[l] = 0;
						else{
							const auto &mer = le.markov_eqn_ref[ti][ctime];
							num[l] = ddt*markov_eqn_vari[mer.e].div[mer.ti].value;
							if(calc_indfac == true) num[l] *= get_indfac(ind,sp.markov_eqn[mer.e]);
				
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
							
							auto t_new = t + ran()*(tt-t);
							pif *= 1.0/(tt-t);
						
							Event enew;
							enew.type = M_TRANS_EV;
							enew.c_gl = UNSET;
							enew.tr_gl = tr_gl;
							enew.t = t_new;
							enew.fixed = false;
							
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
					auto tr = ob.tr;
					
					auto tr_gl = sp.trg_from_tr[c][cl_store][tr];
					const auto &trg = sp.tra_gl[tr_gl];
					
					if(trg.i != c){
						cout << "Not consistent2" << endl; illegal = true; 
						return ev_new;
					}
					
					Event enew;
					enew.type = M_TRANS_EV;
					enew.c_gl = UNSET;
					enew.tr_gl = tr_gl;
					enew.t = ob.t;
					enew.fixed = true;
					c = trg.f; 
					ci = claa.tra[tr].f;
					cisland = claa.comp[ci].island_ref.c;

					enew.c_after = c;
					ev_new.push_back(enew);
					feo--;
					if(feo >= 0) time_fixed = obs[fixed_event_obs[feo]].t; else time_fixed = LARGE;
				}
			}
			t = t_next;
			
			if(k < ev.size()){                // Looks to copy events from the existing sequence 
				auto enew = ev[k];
				const auto &trg = sp.tra_gl[enew.tr_gl];
				if(trg.i != c){   // Ensures that sequence is consistent
					if(c > CODE){
						cout << "Not consistent4" << endl; illegal = true; 
						return ev_new;
					}
					
					auto tr_new = trg.transform[c];
					if(tr_new == UNSET){
						cout << "Not consistent5" << endl; illegal = true; 
						return ev_new;
					}
					
					enew.tr_gl = tr_new;
					auto c_after = sp.tra_gl[tr_new].f;
					enew.c_after = c_after;
					c = c_after;
				}
				else c = trg.f;
				
				ev_new.push_back(enew);
			}
		}
	}		
	
	probif = log(pif);
	
	return ev_new;
}


/// Samples a new set of events based on observation probability
double IndEvSampler::sample_events_prob(const vector <Event> &ev) const
{
	auto pif = 1.0;
	
	const auto &ind = individual[i_store];
	
	const auto &claa = sp.cla[cl_store];
	
	const auto &island = claa.island[isl];
	
	const auto &comp = island.comp;
	
	auto &tp = sp.timepoint;
	
	auto c = c_start;
	auto t = t_start;
	
	auto ci = sp.comp_gl[c].cla_comp[cl_store]; 
	auto cisland = claa.comp[ci].island_ref.c;
	
	auto ti = ti_start;
	
	auto &obs = sp.individual[i_store].obs;
	
	double time_fixed;
	auto feo_max = fixed_event_obs.size();
	int feo = int(feo_max)-1;
	if(feo >= 0) time_fixed = obs[fixed_event_obs[feo]].t; else time_fixed = LARGE;
	
	vector <unsigned int> list;

	// Goes through events on the existing event sequence
	for(auto k = k_start; k <= ev.size(); k++){
		double t_next;
		
		 // Events in the existing sequence for the selected classification are not transfered
		auto ignore = false; 
		
		if(k < ev.size()){
			const auto &e = ev[k];
			
			t_next = e.t;
			if(e.type == M_TRANS_EV){
				if(sp.tra_gl[e.tr_gl].cl == cl_store) ignore = true;	
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
					
					auto &iop = ind_obs_prop[tii];
					
					auto sum = 0.0, sum_op = 0.0;

					// Works out transition prob x obs prob for each possible transition
					vector <double> num(lmax), num_op(lmax);
					for(auto l = 0u; l < lmax; l++){
						const auto &le = co.leave[l];
						auto cf = le.cf;
					
						if(iop[cf] == 0) num[l] = 0;
						else{
							const auto &mer = le.markov_eqn_ref[ti][ctime];
							num[l] = ddt*markov_eqn_vari[mer.e].div[mer.ti].value;
							if(calc_indfac == true) num[l] *= get_indfac(ind,sp.markov_eqn[mer.e]);
				
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
					if(t == t_event){
						const auto &e = ev[list[0]];
						const auto &trg = sp.tra_gl[e.tr_gl];
						auto tr = trg.tr;		
						c = trg.f; 
						ci = claa.tra[tr].f;
						cisland = claa.comp[ci].island_ref.c;
	
						list.erase(list.begin());
						t_event = LARGE; if(list.size() > 0) t_event = ev[list[0]].t;	
						
						feo--;
						if(feo >= 0) time_fixed = obs[fixed_event_obs[feo]].t;
						else time_fixed = LARGE;
					}
					else{
						return -LARGE;
					}
				}
			}
			t = t_next;
			
			if(k < ev.size()){
				c = ev[k].c_after;
			}
		}
	}		
	
	return log(pif);
}


/// Gets ti div time from an actual time (rounding down if equal)
unsigned int IndEvSampler::get_ti_lower(double t) const
{
	return (unsigned int)(ALMOST_ONE*(t- details.t_start)/details.dt);
}


/// Gets ti div time from an actual time (rounding up if equal)
unsigned int IndEvSampler::get_ti_upper(double t) const
{
	return (unsigned int)(OVER_ONE*(t- details.t_start)/details.dt);
}


/// Returns the individual-based factor which acts on Markov equation
double IndEvSampler::get_indfac(const Individual &ind, const MarkovEqn &mar_eqn) const
{
	auto val = 1.0;
	for(auto ie : mar_eqn.ind_eff_mult) val *= ind.exp_ie[ie];
	for(auto fe : mar_eqn.fix_eff_mult) val *= ind.exp_fe[fe];
	
	return val;
}
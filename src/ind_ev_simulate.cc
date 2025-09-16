// This deals with proposals which consist of simulating individuals

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

/// Simulates a new set of events (with certain triger events)
vector <Event> IndEvSampler::simulate_events(unsigned int i, const Event &e_init, double &probif, const vector <TrigEventRef> &trig_event)
{
	illegal = false;
	
	const auto &comp_terminal = sp.comp_terminal;
	
	vector <Event> ev_new;
	
	const auto &ind = individual[i];
	
	if(false) print_trigger_event("TRIG "+ind.name,trig_event);
	
	auto m = 0u;
	
	ev_new.push_back(e_init);
	auto c = e_init.c_after;

	if(c == UNSET) emsg("Cannot start unset");

	auto t = e_init.tdiv;

	vector <FutureNMEvent> future_nm;

	auto f = 0u;
	for(auto cl = 0u; cl < sp.ncla; cl++){
		add_future_nm_event(c,cl,t,i,future_nm,m,trig_event);
	}
	
	if(illegal) return ev_new;
	
	auto ti_start = get_ti(t);

	vector <double> trans_rate;
	double R;
	
	auto tnext_ev = find_tnext(m,trig_event,f,future_nm);

	for(auto ti = ti_start; ti < T; ti++){	
		double tend = ti+1;
		do{
			auto tnext = tend; if(tnext_ev < tnext) tnext = tnext_ev;
			
			do{
				get_trans_rate(c,ti,ind,trans_rate,R);
				
				if(R == 0){ t = tnext; break;}
			
				auto tnew = t-log(ran())/R;
				if(tnew > tnext) tnew = tnext;
				
				probif -= R*(tnew-t);
				t = tnew;
				
				if(tnew == tnext) break;
				
				vector <double> sum_store;
				auto sum = 0.0;
				for(auto va : trans_rate){
					sum += va;
					sum_store.push_back(sum);
				}
				
				auto z = ran()*R;
			
				auto k = 0u; while(k < sum_store.size() && z > sum_store[k]) k++;
				if(k == sum_store.size()) emsg("Problem simulate sampling");
				
				probif += log(trans_rate[k]);
				
				auto trg = sp.comp_gl[c].tr_leave_markov[k];
		
				const auto &tra = sp.tra_gl[trg];
				c = tra.f; 
				sim_add_event(M_TRANS_EV,trg,tnew,tra.cl,UNSET,c,false,ev_new);
				if(c == UNSET){
					return ev_new;
				}
				
				add_future_nm_event(c,tra.cl,t,i,future_nm,m,trig_event);
				tnext_ev = find_tnext(m,trig_event,f,future_nm);
				if(tnext_ev < tnext) tnext = tnext_ev;
			}while(true);
			
			if(tnext == tend) break;
			
			if(f < future_nm.size() && tnext == future_nm[f].tdiv){ // Does a future nm event
				const auto &fnm = future_nm[f];
				auto trg = fnm.trg;
				const auto &tra = sp.tra_gl[trg];
			
				if(tra.i != c) trg = sp.tr_trans(trg,c);
				
				const auto &tra2 = sp.tra_gl[trg];
				if(tra2.i != c) emsg("Shoudl be c");
			
				c = tra2.f;
				sim_add_event(NM_TRANS_EV,trg,t,tra.cl,UNSET,c,false,ev_new);
				if(c == UNSET){
					return ev_new;
				}
				
				add_future_nm_event(c,tra.cl,t,i,future_nm,m,trig_event);
				f++;
				tnext_ev = find_tnext(m,trig_event,f,future_nm);
			}
			else{                     
				if(m < trig_event.size()){                   		   // Does a trigger event
					auto &te = trig_event[m];
					switch(te.type){
					case TRIG_OBS_TRANS_EV: 
						{
							auto cl = sim_add_data_event(te,c,i,ev_new,probif);
							
							if(c == UNSET || illegal){
								return ev_new;	
							}
							
							add_future_nm_event(c,cl,t,i,future_nm,m,trig_event);
						}
						break;
						
					case TRIG_MOVE_EV: 
						{
							const auto &ev = sp.individual[i].ev[te.ref];
							c = sp.update_c_comp(c,ev.cl,ev.move_c);
							sim_add_event(MOVE_EV,UNSET,te.tdiv,ev.cl,ev.move_c,c,true,ev_new);
							add_future_nm_event(c,ev.cl,t,i,future_nm,m,trig_event);	
							if(c == UNSET) emsg("err");
						}
						break;
						
					case TRIG_LEAVE_EV:
						sim_add_event(LEAVE_EV,UNSET,te.tdiv,UNSET,UNSET,UNSET,true,ev_new);
						return ev_new;
					
					default:
						emsg("op should not be here");
						break;
					}
					
					m++; 
					
					tnext_ev = find_tnext(m,trig_event,f,future_nm);
				}
			}

			if(illegal) return ev_new;
		}while(true);
		
		if(comp_terminal[c]) break;
	}

	if(m != trig_event.size() || f != future_nm.size()){
		if(m != trig_event.size()){ illegal = true; } 
		else{
			emsg("f not right2");
		}
	}

	return ev_new;
}


/// Finds the next time (trigger or nm event) 
double IndEvSampler::find_tnext(unsigned int m, const vector <TrigEventRef> &trig_event, unsigned int f, const vector <FutureNMEvent> &future_nm) const
{
	double te_next = LARGE; 
	if(m < trig_event.size()) te_next = trig_event[m].tdiv; 
	if(f < future_nm.size()){
		if(future_nm[f].tdiv < te_next) te_next = future_nm[f].tdiv;
	}
	return te_next;
}


/// Adds a future non-Markovian event
void IndEvSampler::add_future_nm_event(unsigned int c, unsigned int cl, double t, unsigned int i, vector <FutureNMEvent> &future_nm, unsigned int m, const vector <TrigEventRef> &trig_event)
{
	const CompGlTransGroup &tlg = sp.comp_gl[c].tra_leave_group[cl];
	
	if(tlg.markov != false || tlg.tr_list.size() == 0) return;
	
	double tmax = details.T;
	
	// Looks to see if future event happens after trig event
	
	auto sim_nm_event = true; 
	auto end = false;
	for(auto j = m; j < trig_event.size(); j++){
		const auto &te = trig_event[j];
		if(te.tdiv > t){
			switch(te.type){
			case TRIG_OBS_TRANS_EV: 
				{
					const auto &ob = sp.individual[i].obs[te.ref];
					if(cl == ob.cl){
						for(auto tr : tlg.tr_list){
							if(get_trans_obs_prob(tr,ob) != 0){ sim_nm_event = false; break;}
						}
						end = true;
					}
				}
				break;
			
			case TRIG_MOVE_EV:
				emsg("Should not have a trig move with nm");
				break;
								
			case TRIG_LEAVE_EV: 
				end = true; 
				break;
		
			default: break;
			}
			
			if(end == true){ tmax = te.tdiv; break;}
		}
	}
	
	if(!sim_nm_event) return;

	const auto &ind = individual[i];
		
	auto ti = get_ti_lower(t);
			
	unsigned int k;
	if(tlg.branch){
		const auto &tr_list = tlg.tr_list;
		auto B = tr_list.size();
	
		// This selects the branch
		// Note tlg.all_branches is accounted for because the sum does not need to be 1
		
		auto sum = 0.0;
		vector <double> sum_store;
		for(auto b = 0u; b < B; b++){
			const auto &tra = sp.tra_gl[tr_list[b]];
			sum += eqn[tra.bp.eq_ref].calculate_indfac(ind,ti,popnum_t[ti],precalc);
			sum_store.push_back(sum);
		}
			
		auto z = sum*ran();
		k = 0u; while(k < B && z > sum_store[k]) k++;
		if(k == B) emsg("Should not be B");		
	}
	else k = 0;
	
	auto tr = tlg.tr_list[k];
	const auto &tra = sp.tra_gl[tr];
	
	const auto &dp = tra.dist_param;
	
	double ts;
	string warn;
	
	switch(tra.type){	
	case EXP_RATE: case EXP_MEAN: ts = UNSET; emsg("Should not be rate2"); break;
	
	case EXP_RATE_NM: 	
		{
			auto rate = eqn[dp[0].eq_ref].calculate_indfac(ind,ti,popnum_t[ti],precalc);
			ts = t+exp_rate_sample(rate*dt,warn);
		}
		break;
		
	case EXP_MEAN_NM: 	
		{
			auto mean = eqn[dp[0].eq_ref].calculate_indfac(ind,ti,popnum_t[ti],precalc);
			ts = t+exp_mean_sample(mean/dt,warn);
		}
		break;

	case GAMMA:
		{
			auto mean = eqn[dp[0].eq_ref].calculate_indfac(ind,ti,popnum_t[ti],precalc);
			auto cv = eqn[dp[1].eq_ref].calculate_indfac(ind,ti,popnum_t[ti],precalc);
			ts = t+gamma_sample(mean/dt,cv,warn);	
		}
		break;
	
	case ERLANG:
		{
			auto mean = eqn[dp[0].eq_ref].calculate_indfac(ind,ti,popnum_t[ti],precalc);
			auto shape = eqn[dp[1].eq_ref].calculate_indfac(ind,ti,popnum_t[ti],precalc);
			auto cv = sqrt(1.0/shape);
			ts = t+gamma_sample(mean/dt,cv,warn);
		}
		break;

	case LOG_NORMAL:
		{
			auto mean = eqn[dp[0].eq_ref].calculate_indfac(ind,ti,popnum_t[ti],precalc);
			auto cv = eqn[dp[1].eq_ref].calculate_indfac(ind,ti,popnum_t[ti],precalc);
			ts = t + lognormal_sample(mean/dt,cv,warn);
		}
		break;
		
	case WEIBULL:
		{
			auto scale = eqn[dp[0].eq_ref].calculate_indfac(ind,ti,popnum_t[ti],precalc);
			auto shape = eqn[dp[1].eq_ref].calculate_indfac(ind,ti,popnum_t[ti],precalc);
			ts = t+weibull_sample(scale/dt,shape,warn);
		}
		break;
		
	case PERIOD:
		{
			auto time = eqn[dp[0].eq_ref].calculate_indfac(ind,ti,popnum_t[ti],precalc);
			ts = t+period_sample(time,warn)/dt;
		}
		break;
		
	default: ts = UNSET; emsg("Option should not be here"); break;
	}	

	if(warn != "") sp.sampling_error(tr,warn);
	
	if(ts < tmax){
		FutureNMEvent fnme; fnme.tdiv = ts; fnme.trg = tr;
		
		auto k = future_nm.size();
		if(k == 0) future_nm.push_back(fnme);
		else{
			while(k > 0 && future_nm[k-1].tdiv > ts) k--;
			future_nm.insert(future_nm.begin()+k,fnme);
		}
	}
	else{
		if(end == true){
			illegal = true;
		}
	}		
}


/// Adds a new event to the event list
void IndEvSampler::sim_add_event(EventType type, unsigned int tr_gl, double t, unsigned int cl, unsigned int move_c, unsigned int c, bool observed, vector <Event> &ev_new) const 
{
	Event enew;
	enew.type = type;
	enew.move_c = move_c;
	enew.cl = cl;
	enew.tr_gl = tr_gl;
	enew.tdiv = t;
	enew.c_after = c;
	enew.observed = observed;
	enew.inf_node_ref = UNSET;
	enew.ind_inf_from = IndInfFrom();							
							
	ev_new.push_back(enew);
}
	

/// Adds data event to the individual time line
unsigned int IndEvSampler::sim_add_data_event(const TrigEventRef &te, unsigned int &c, unsigned int i, vector <Event> &ev_new, double &probif)
{
	const auto &ob = sp.individual[i].obs[te.ref];
	auto cl = ob.cl; 
	
	auto sum = 0.0;
	vector <double> sum_store, val_store; 
	for(auto tr : sp.comp_gl[c].tra_leave_group[cl].tr_list){
		auto val = get_trans_obs_prob(tr,ob);
		sum += val;
		val_store.push_back(val);
		sum_store.push_back(sum);
	}

	if(sum == 0){ illegal = true; return cl;}
	else{
		auto z = ran()*sum;
		auto k = 0u; while(k < sum_store.size() && z > sum_store[k]) k++;
		if(k == sum_store.size()){
			emsg("Cannot select");
		}
		
		auto tr = sp.comp_gl[c].tra_leave_group[ob.cl].tr_list[k];

		if(sum_store.size() > 0){
			probif += log(val_store[k]/sum);
		}
		
		const auto &tra = sp.tra_gl[tr];
		c = tra.f;
		auto type = M_TRANS_EV;	if(tra.nm_trans_ref != UNSET) type = NM_TRANS_EV;
		sim_add_event(type,tr,te.tdiv,cl,UNSET,c,true,ev_new);
	}
	
	return cl;
}


/// Gets all the possible Markovian transitions from c
void IndEvSampler::get_trans_rate(unsigned int c, unsigned int ti, const Individual &ind, vector <double> &trans_rate, double &R) const
{	
	R = 0.0;
	if(c == UNSET) return;
	
	auto &trl = sp.comp_gl[c].tr_leave_markov;
	auto N = trl.size();
	if(N == 0) return;
	
	trans_rate.resize(N);
	
	double r;
	for(auto k = 0u; k < N; k++){	
		auto trg = trl[k];
		auto m = sp.tra_gl[trg].markov_eqn_ref;
		
		const auto &mev = markov_eqn_vari[m];
		const auto &me = sp.markov_eqn[m];
			
		if(me.time_vari == true) r = mev.div[ti].value;
		else r = mev.div[0].value;
		if(me.ind_variation) r *= get_indfac(ind,me);
		if(r < TINY) r = 0;
		else{
			if(sp.obs_trans_exist){  // Accounts for observation probability
				for(auto m : sp.obs_trans_eqn_ref[trg][ti]){
					auto val = obs_trans_eqn_value[m][ti];
					if(val != 0){
						if(val == LI_WRONG) r = 0;
						else r *= exp(val);
					}
				}
			}
		}
		
		trans_rate[k] = r;
		R += r;
	}
}


/// Samples a new set of events based on observation probability
double IndEvSampler::simulate_events_prob(unsigned int i, const vector <Event> &ev, const vector <TrigEventRef> &trig_event) const
{
	double probif = 0.0;
	
	const auto &comp_terminal = sp.comp_terminal;
	
	const auto &ind = individual[i];

	const auto &e_init = ev[0];
	auto c = e_init.c_after;

	auto t = e_init.tdiv;
	
	auto ti_start = get_ti(t);
	
	vector <double> trans_rate;
	double R;
	
	auto E = ev.size();
	auto e = 1u;
	double t_ev_next;
	if(e < E) t_ev_next = ev[e].tdiv; else t_ev_next = LARGE;
	
	for(auto ti = ti_start; ti < T; ti++){
		double tend = ti+1;
	
		do{
			get_trans_rate(c,ti,ind,trans_rate,R);
			
			if(tend < t_ev_next){
				probif -= R*(tend-t);
				t = tend;
				break;
			}
			else{
				probif -= R*(t_ev_next-t);
				t = t_ev_next;
		
				if(e == E) emsg("E problem");
				const auto &eve = ev[e];
				
				if(eve.observed){  // An observed event
					if(eve.type == NM_TRANS_EV){
						auto k = 0u; while(k < trig_event.size() && trig_event[k].tdiv != eve.tdiv) k++;
						if(k == trig_event.size()) emsg("Trig event not allowed");
						
						const auto &ob = sp.individual[i].obs[trig_event[k].ref];
						if(ob.type != OBS_TRANS_EV && ob.type != OBS_SINK_EV){
							emsg("Should be an observed transition");
						}

						auto cl = ob.cl; 
						
						auto &tr_list =  sp.comp_gl[c].tra_leave_group[cl].tr_list;
						if(tr_list.size() > 1){
							auto sum = 0.0;
							auto val = LARGE;
							for(auto tr : tr_list){
								auto prob = get_trans_obs_prob(tr,ob);
								if(tr == eve.tr_gl) val = prob;
								sum += prob;		
							}
							if(val == LARGE) emsg("cannot find transition");
							if(val != sum) probif += log(val/sum);
						}
					}		
				}  
				else{
					auto trg = eve.tr_gl;
					
					switch(eve.type){
 					case M_TRANS_EV:  // A Markovian event
						{
							const auto &tlm = sp.comp_gl[c].tr_leave_markov;
							auto k = 0u; while(k < tlm.size() && tlm[k] != trg) k++;
							if(k == tlm.size()) emsg("Could not find tlm1");
							
							probif += log(trans_rate[k]);
						}
						break;
						
					case NM_TRANS_EV:  // A non-markovian event
						break;
						
					default: emsg("Should not be here"); break;
					}
				}
				
				c = eve.c_after;
				e++;
				if(e < E) t_ev_next = ev[e].tdiv; else t_ev_next = LARGE;
			}
		}while(true);
		if(c == UNSET || comp_terminal[c]) break;
	}
	
	return probif;
}


/// Prints trigger events
void IndEvSampler::print_trigger_event(string te, const vector <TrigEventRef> &trig_event) const
{
	cout << te << ":" << endl;
	
	for(const auto &ev : trig_event){
		switch(ev.type){
		case TRIG_OBS_TRANS_EV: cout << "OBS TRANS"; break;
		case TRIG_MOVE_EV: cout << "MOVE"; break;
		case TRIG_LEAVE_EV: cout << "LEAVE"; break;
		}
		cout << ", " << ev.tdiv << "   ";
	}
	cout << " Trig events" << endl;
}

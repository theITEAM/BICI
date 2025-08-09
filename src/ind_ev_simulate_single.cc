// This file deals with simulating individuals in a single classification

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
using namespace std;

#include "ind_ev_sampler.hh"
#include "utils.hh"

/// Simulates a new set of events for a given classification (with certain trigger events)
vector <Event> IndEvSampler::simulate_single_events(unsigned int i, unsigned int cl, double &probif, const vector <TrigEventRef> &trig_event)
{
	illegal = false;
	
	vector <Event> ev_new;
	
	const auto &ind = individual[i];
	
	const auto &event = ind.ev;
	
	if(false) print_trigger_event("TRIG "+ind.name,trig_event);
	
	auto m = 0u;
	
	const auto &e_init = event[0];
	ev_new.push_back(e_init);
	auto c = e_init.c_after;
	auto t = e_init.tdiv;
	
	auto ti_start = get_ti(t);
	
	vector <FutureNMEvent> future_nm;

	auto f = 0u;
	add_future_nm_event(c,cl,t,i,future_nm,m,trig_event);

	if(illegal) return ev_new;
	
	vector <double> trans_rate;
	double R;
	
	auto e = 1u;
	
	auto tnext_ev = find_tnext_single(m,trig_event,f,future_nm,e,event);		
	
	const auto &term = sp.comp_terminal_cl[cl];
	auto terminal = false;
	
	for(auto ti = ti_start; ti < T; ti++){		
		double tend = ti+1;

		do{
			auto tnext = tend; if(tnext_ev < tnext) tnext = tnext_ev;
			
			if(terminal) t = tnext;
			else{
				do{
					get_trans_rate_single(c,cl,ti,ind,trans_rate,R);
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
					if(k == sum_store.size()) emsg("Problem simulate sampling2");
					
					probif += log(trans_rate[k]);
					
					auto trg = sp.comp_gl[c].tra_leave_group[cl].tr_list[k];
				
					const auto &tra = sp.tra_gl[trg];
					c = tra.f; 
					sim_add_event(M_TRANS_EV,trg,tnew,tra.cl,UNSET,c,false,ev_new);
					if(c == UNSET) return ev_new;
					
					add_future_nm_event(c,tra.cl,t,i,future_nm,m,trig_event);
				
					tnext_ev = find_tnext_single(m,trig_event,f,future_nm,e,event);		
		
					if(tnext_ev < tnext) tnext = tnext_ev;
				}while(true);
				if(term[c]) terminal = true;
			}
			
			if(tnext == tend) break;
			
			if(e < event.size() && tnext == event[e].tdiv){   // Inserts event from existing sequence
				auto eve = event[e];
				auto cll = eve.cl;
				if(cll != cl || eve.c_after == UNSET){			

					switch(eve.type){
					case M_TRANS_EV: case NM_TRANS_EV:
						{
							auto trg = eve.tr_gl;
							const auto &tra = sp.tra_gl[trg];	
					
							if(tra.i != c){
								trg = sp.tr_trans(trg,c);
								if(trg == UNSET){ illegal = true; return ev_new;}	
								eve.tr_gl = trg;								
							}							
							
							const auto &tra2 = sp.tra_gl[trg];		
							
							if(tra2.i != c) emsg("Should be c other");
					
							c = tra2.f;
							eve.c_after = c;
						}
						break;
					
					case LEAVE_EV: 
						c = UNSET;
						break;
						
					case ENTER_EV:
						emsg("Should not enter");
						break;

					case MOVE_EV:
						{
							c = sp.update_c_comp(c,cll,eve.move_c);
							eve.c_after = c;
						}
						break;
					}
					
					ev_new.push_back(eve);
				}
				e++;
				tnext_ev = find_tnext_single(m,trig_event,f,future_nm,e,event);	
			}
			else{
				if(f < future_nm.size() && tnext == future_nm[f].tdiv){ // Does a future nm event
					const auto &fnm = future_nm[f];
				
					auto trg = fnm.trg;
					const auto &tra = sp.tra_gl[trg];
				
					if(tra.i != c){
						trg = sp.tr_trans(trg,c);
					}			
					const auto &tra2 = sp.tra_gl[trg];
				
					if(tra2.i != c) emsg("Shoudl be c future");
					c = tra2.f;
					sim_add_event(NM_TRANS_EV,trg,t,tra.cl,UNSET,c,false,ev_new);
					if(c == UNSET) return ev_new;
					
					add_future_nm_event(c,tra.cl,t,i,future_nm,m,trig_event);
					f++;
					tnext_ev = find_tnext_single(m,trig_event,f,future_nm,e,event);		
				}
				else{                     
					if(m < trig_event.size()){                   		   // Does a trigger event
						auto &te = trig_event[m];
						switch(te.type){
						case TRIG_OBS_TRANS_EV:
							{
								auto cl = sim_add_data_event(te,c,i,ev_new,probif);
								if(c == UNSET || illegal) return ev_new;
							
								add_future_nm_event(c,cl,t,i,future_nm,m,trig_event);
							}
							break;
							
						case TRIG_MOVE_EV: 
							{
								const auto &ev = sp.individual[i].ev[te.ref];
								c = sp.update_c_comp(c,ev.cl,ev.move_c);
								sim_add_event(MOVE_EV,UNSET,te.tdiv,ev.cl,ev.move_c,c,true,ev_new);
								add_future_nm_event(c,ev.cl,t,i,future_nm,m,trig_event);		
								tnext_ev = find_tnext(m,trig_event,f,future_nm);				
								if(c == UNSET) emsg("err");								
							}
							break;
							
						case TRIG_LEAVE_EV:
							sim_add_event(LEAVE_EV,UNSET,te.tdiv,UNSET,UNSET,UNSET,true,ev_new);
							return ev_new;
					
						default:
							emsg("op should not be here2");
							break;
						}
						
						m++; 
						tnext_ev = find_tnext_single(m,trig_event,f,future_nm,e,event);		
					}
				}
			}
			if(illegal) return ev_new;
			
			if(c == UNSET) break;
		}while(true);
		
		if(c == UNSET) break;
	}

	if(m != trig_event.size() || f != future_nm.size()){
		if(m != trig_event.size()) illegal = true;
		else{
			if(c != UNSET) emsg("f not right");
		}
	}

	return ev_new;
}


/// Simulates a new set of events based on observation probability
double IndEvSampler::simulate_single_events_prob(unsigned int i, unsigned int cl, const vector <Event> &ev, const vector <TrigEventRef> &trig_event) const
{
	double probif = 0.0;
	
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
	
	const auto &term = sp.comp_terminal_cl[cl];
	auto terminal = false;
	
	for(auto ti = ti_start; ti < T; ti++){
		double tend = ti+1;
		
		do{
			auto do_ev = false;
			
			if(terminal == false){
				get_trans_rate_single(c,cl,ti,ind,trans_rate,R);
			
				if(tend < t_ev_next){
					probif -= R*(tend-t);
					t = tend;
					break;
				}
				else{
					probif -= R*(t_ev_next-t);
					t = t_ev_next;
					do_ev = true;
				}
				
				if(term[c]) terminal = true;
			}
			else{
				if(tend < t_ev_next){ t = tend; break;}
				else{ t = t_ev_next; do_ev = true;}
			}
			
			if(do_ev){ 
				if(e == E) emsg("E problem");
				const auto &eve = ev[e];
				
				if(eve.cl != cl || eve.c_after == UNSET){		// Event from another classification		
				}
				else{
					if(eve.observed){  // An observed event
						if(eve.type == NM_TRANS_EV){
							auto k = 0u; while(k < trig_event.size() && trig_event[k].tdiv != eve.tdiv) k++;
							if(k == trig_event.size()) emsg("Trig event not allowed");
							
							const auto &ob = sp.individual[i].obs[trig_event[k].ref];
							if(ob.type != OBS_TRANS_EV) emsg("Should be an observed transition");
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
								const auto &tlm = sp.comp_gl[c].tra_leave_group[cl].tr_list;
								
								auto k = 0u; while(k < tlm.size() && tlm[k] != trg) k++;
								if(k == tlm.size()) emsg("Could not find tlm");
								
								probif += log(trans_rate[k]);
							}
							break;
							
						case NM_TRANS_EV:  // A non-markovian event
							break;
							
						default: emsg("Should not be here"); break;
						}
					}
				}
				
				c = eve.c_after;
				e++;
				if(e < E) t_ev_next = ev[e].tdiv; else t_ev_next = LARGE;
			}
		}while(true);
		if(c == UNSET) break;
	}

	return probif;
}


/// Finds the next time (trigger or nm event) 
double IndEvSampler::find_tnext_single(unsigned int m, const vector <TrigEventRef> &trig_event, unsigned int f, const vector <FutureNMEvent> &future_nm, unsigned int e, const vector <Event> &ev) const
{
	double te_next = LARGE; 
	if(m < trig_event.size()) te_next = trig_event[m].tdiv; 
	if(f < future_nm.size()){
		if(future_nm[f].tdiv < te_next) te_next = future_nm[f].tdiv;
	}
	if(e < ev.size()){
		if(ev[e].tdiv < te_next) te_next = ev[e].tdiv;
	}
	return te_next;
}


/// Gets all the possible Markovian transitions from c
void IndEvSampler::get_trans_rate_single(unsigned int c, unsigned int cl, unsigned int ti, const Individual &ind, vector <double> &trans_rate, double &R) const
{	
	R = 0.0;
	if(c == UNSET) return;
 
	const auto &tlg = sp.comp_gl[c].tra_leave_group[cl];
	if(!tlg.markov) return;

	const auto &trl = tlg.tr_list;
	auto N = trl.size();
	if(N == 0) return;
	
	trans_rate.resize(N);
	
	double r;
	for(auto k = 0u; k < N; k++){	
		auto trg = trl[k];
		
		const auto &tra = sp.tra_gl[trg];
		if(tra.f == UNSET) r = 0; // Avoid transitions which lead to sink
		else{
			auto m = tra.markov_eqn_ref;
			
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
		}
		
		trans_rate[k] = r;
		R += r;
	}
	
	auto sum = 0.0; for(auto va : trans_rate) sum += va;
}

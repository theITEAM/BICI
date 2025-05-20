// Functions related to local individual-based proposals

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


/// Finds all the possible local changes which can be made on an individual
vector <LocalIndChange> StateSpecies::local_ind_change(unsigned int i, unsigned int cl) const
{
	vector <LocalIndChange> licha;
	
	const auto &indd = sp.individual[i];
	
	if(!(i >= sp.nindividual_obs || indd.simulation_needed)) return licha;
		
	const auto &ev = individual[i].ev;
	
	const auto &claa = sp.cla[cl];
	
	// First constructs timeline in the given classification
	
	vector <EventCl> timeline;
	auto C = claa.comp.size();

	auto c = C;
	for(auto e = 0u; e < ev.size(); e++){
		const auto &eve = ev[e];
		
		auto c_aft = C;
		auto c_after = eve.c_after;
		if(c_after != UNSET) c_aft = sp.comp_gl[c_after].cla_comp[cl];
		
		EventCl evc;
		evc.c_bef = c;
		evc.c_aft = c_aft;
		evc.t = eve.t;
		evc.e = e;
		
		switch(eve.type){
		case ENTER_EV: 
			{
				evc.trs.type = ENTER_SW;
				evc.trs.trc = c_aft;
				timeline.push_back(evc);
			}
			break;
			
		case LEAVE_EV:
			{
				evc.trs.type = LEAVE_SW;
				evc.trs.trc = c;
				timeline.push_back(evc);
			}				
			break;
				
		case MOVE_EV:
			if(eve.cl == cl){
				evc.trs.type = MOVE_SW;
				evc.trs.trc = UNSET;
				timeline.push_back(evc);
			}
			break;
	
		case NM_TRANS_EV: case M_TRANS_EV:
			{
				const auto &tra = sp.tra_gl[eve.tr_gl];
				if(eve.cl == cl || tra.variety != NORMAL){
					evc.trs.trc = tra.tr;
					if(tra.i == UNSET) evc.trs.type = SOURCE_SW;
					else{
						if(tra.f == UNSET) evc.trs.type = SINK_SW;
						else evc.trs.type = TRANS_SW;
					}
					timeline.push_back(evc);
				}
			}
			break;
			
		default: emsg("Should not be event"); break;
		}
		
		c = c_aft;
	}
	
	if(c != C){                             // Adds in a final event time (if neccesary)
		EventCl evc;
		evc.c_bef = c;
		evc.c_aft = C;
		evc.t = details.t_end;
		evc.e = ev.size();
		
		evc.trs.type = LEAVE_SW;
		evc.trs.trc = c;
		timeline.push_back(evc);
	}
	
	if(false){
		for(auto evc : timeline){
			vector <TrSwap> trs; trs.push_back(evc.trs);
			cout << "t=" << evc.t << "  e=" << evc.e << "  c_bef=" << evc.c_bef 
						<< "  c_aft=" << evc.c_aft << "  trs=" << sp.tr_swap_name(cl,trs) << endl;  
		}
		emsg("WW");
	}
	
	auto no_enter = indd.init_c_set; 
	switch(sp.init_cond.type){
	case INIT_POP_FIXED: case INIT_POP_NONE: no_enter = true; break;
	default: break;
	}

	// Next works out which parts of timeline could be swapped for an alternative
	// These possibilities are added to licha and lc_ref stores a reference
	for(auto num = 0u; num <= IND_LOCAL_TRANS_MAX; num++){
		if(num <= timeline.size()){
			auto jmax = timeline.size()-num+1;
			if(num == 0) jmax = timeline.size();
			for(auto j = 0u; j < jmax; j++){
				auto c_beg = timeline[j].c_bef;
				unsigned int c_end;
				if(num == 0) c_end = c_beg;
				else c_end = timeline[j+num-1].c_aft;
				//if(claa.tr_swap[c_beg][c_end].size() > 1){		
				if(claa.tr_swap_possible[c_beg][c_end]){
					auto vec = sp.get_vec_tr_swap_mid(j,num,timeline);
					
					vector <TrSwap> trs; 
					for(auto jj = j; jj < j+num; jj++) trs.push_back(timeline[jj].trs);
					
					auto si = claa.hash_swap.existing(vec);
					
					if(si != UNSET){
						const auto &sw = claa.swap[si];
						if(!(no_enter == true && sw.start.size() > 0 && sw.start[0].type == ENTER_SW)){
							LocalIndChange lic;
							lic.i = i;
							lic.cl = cl;
							if(j > 0) lic.tmin = timeline[j-1].t; 
							else lic.tmin = timeline[j].t;
							if(j+num < timeline.size()) lic.tmax = timeline[j+num].t; 
							else lic.tmax = timeline[j+num-1].t;
							
							for(auto k = 0u; k < num; k++){
								lic.remove_e.push_back(timeline[j+k].e);
							}
							
							for(auto m = 0u; m < sw.swap.size(); m++){
								lic.tr_swap = sw.swap[m];
								lic.swap_rep_ref = sw.swap_rep_ref[m];
								licha.push_back(lic);					
							}
						}		
					}
				}
			}
		}
	}
	
	if(false){
		cout << "event" << endl;
		print_event(ev);
	
		cout << endl << endl << "LOCAL CHANGE" << endl;
		for(auto &lich : licha) print_local_ind_change(lich);
		emsg("done");
	}

	return licha;
}


/// Prints a local individual change
void StateSpecies::print_local_ind_change(const LocalIndChange &lich) const
{
	cout << individual[lich.i].name;
	cout << "  t=" << lich.tmin << " - " << lich.tmax;
	cout << "  " << sp.cla[lich.cl].swap_rep[lich.swap_rep_ref].name;
	if(false){
		cout << "  remove="; for(auto e : lich.remove_e) cout << e << ",";
		cout << "  SWAP=" << sp.tr_swap_name(lich.cl,lich.tr_swap);
	}
	cout << endl;
}	


/// Creates a new event sequence based on a given local change
double StateSpecies::create_local_change(double &timefac, const LocalIndChange &lich, vector <Event> &ev_new, Event &enew, unsigned int nind_obs, const vector <Event> &ev_rev, LocalDir dir) const
{
	auto i = lich.i;
	auto cl = lich.cl;
	const auto &claa = sp.cla[cl];
	
	const auto &ind = individual[i];
	const auto &event = ind.ev;

	auto iobs = i; if(iobs >= nind_obs) iobs = nind_obs;
	
	const auto &indd = sp.individual[iobs];
	
	const auto &obs = indd.obs;
	
	const auto &rem = lich.remove_e;
	auto tmin = lich.tmin;
	auto tmax = lich.tmax;
	auto dt = tmax-tmin;
	
	auto prob_try = 1.0;
	
	vector <double> time_list;

	for(const auto &trs : lich.tr_swap){
		switch(trs.type){
		case ENTER_SW: case SOURCE_SW:
			time_list.push_back(event[rem[0]].t);
			break;
			
		case LEAVE_SW:
			{
				auto e_end = rem[rem.size()-1];
				if(e_end == event.size()){
					time_list.push_back(details.t_end);
				}
				else{
					time_list.push_back(event[e_end].t);
				}
			}
			break;
		
		case SINK_SW:
			time_list.push_back(event[rem[rem.size()-1]].t);
			break;
		
		case TRANS_SW:
			switch(dir){
			case LOCAL_FORWARD: 
				time_list.push_back(tmin+ran()*dt); 
				break;
		
			case LOCAL_REVERSE:
				{
					auto ee_sel = UNSET;
					for(auto ee = 0u; ee < ev_rev.size(); ee++){
						const auto &evr = ev_rev[ee];
						if(evr.t > tmin && evr.t < tmax && (evr.type == M_TRANS_EV || evr.type == NM_TRANS_EV)){
							const auto &trg = sp.tra_gl[evr.tr_gl];
							if(trg.cl == cl && trg.tr == trs.trc && find_in(time_list,evr.t) == UNSET){
								if(ee_sel == UNSET) ee_sel = ee;
								else emsg("Two fit");
							}
						}							
					}
					if(ee_sel == UNSET) emsg("Could not find event");
						
					time_list.push_back(ev_rev[ee_sel].t); 
				}
				break;
			}
			break;
			
		case MOVE_SW: emsg("Should not be move"); break;
		}
	}
	
	auto K = time_list.size();
	
	if(K > 1) sort(time_list.begin(),time_list.end());
	
	int dE = int(lich.tr_swap.size()) - int(rem.size());

	// Constructs new event sequence
	
	auto E = event.size();
	auto R = rem.size();

	auto m = 0u;
	auto M = obs.size();
	
	auto k = 0u;                            // Indexes new events being added
	auto r = 0u;                            // Indexes removal of events
	const auto &tr_swap = lich.tr_swap;
	
	if(tr_swap.size() != K) emsg("SHould be TL");
	
	// Adds events before time period
	auto e = 0u;
	while(e < E && event[e].t < tmin){
		ev_new.push_back(event[e]); 
		e++;
	}
	
	auto c_new = UNSET;
	if(e > 0) c_new = ev_new[e-1].c_after;
	auto c_old = c_new;
	
	while(m < M && obs[m].t < tmin) m++;
	
	auto t_old = tmin, t_new = tmin;
	
	double t_event_old, t_event_new, t_obs;

	if(e < E) t_event_old = event[e].t; else t_event_old = LARGE;
	if(k < K) t_event_new = time_list[k]; else t_event_new = LARGE;
	if(m < M) t_obs = obs[m].t; else t_obs = LARGE;
	
	auto next = DO_NOTHING_LOC;
	
	auto Rint = 0.0;
	
	do{
		auto t_next = LARGE;
		if(t_event_new < t_next){ t_next = t_event_new; next = NEW_EVENT_LOC;}
		if(t_event_old < t_next){ t_next = t_event_old; next = OLD_EVENT_LOC;}
		if(t_obs < t_next){ t_next = t_obs; next = OBS_LOC;}
		
		if(t_next > tmax) break;
		
		switch(next){
		case NEW_EVENT_LOC:                   // Adds in new events
			{
				unsigned int c_after, tr_gl;
				const auto &trs = tr_swap[k];
				if(!(trs.type == LEAVE_SW && e == E)){
					switch(trs.type){
					case ENTER_SW:
						{
							enew.type = ENTER_EV;
							tr_gl = UNSET;
						
							if(event[e].type != ENTER_EV) emsg("Should be enter here");
							auto c_old = event[e].c_after;
							c_after = sp.comp_global_convert[c_old][cl][trs.trc];
							
							auto entref = indd.enter_ref;
							if(entref != UNSET){
								auto prob_new = sp.calc_enter_prob(c_after,entref,obs_eqn_value);
								if(prob_new == 0) return 0;	
								else{
									auto prob_old = sp.calc_enter_prob(c_old,entref,obs_eqn_value);
									prob_try *= (prob_new+PROB_MOD)/(prob_old+PROB_MOD);
								}
							}
						}
						break;
						
					case SOURCE_SW:
						{
							if(event[e].type != M_TRANS_EV) emsg("Should be M_TRANS_EV");
							enew.type = M_TRANS_EV;
							auto c_old = event[e].c_after;
							auto tr = trs.trc;
							auto cnew = claa.tra[tr].f;
							c_after = sp.comp_global_convert[c_old][cl][cnew];
							tr_gl = sp.trg_from_tr[c_after][cl][tr];
						}
						break;
					
					case LEAVE_SW:
						{
							if(event[e].type != LEAVE_EV) emsg("Should be M_TRANS_EV");
							enew.type = LEAVE_EV;
							tr_gl = UNSET;
							c_after = UNSET;
						}
						break;
						
					case SINK_SW:
						{		
							auto tr = trs.trc;
							tr_gl = sp.trg_from_tr[c_new][cl][tr];
							auto &tra = sp.tra_gl[tr_gl];
							enew.type = M_TRANS_EV;
							if(tra.nm_trans_ref != UNSET) enew.type = NM_TRANS_EV;
							c_after = UNSET;
						}
						break;
						
					 case TRANS_SW:
						{		
							auto tr = trs.trc;
							tr_gl = sp.trg_from_tr[c_new][cl][tr];
							auto &tra = sp.tra_gl[tr_gl];
							enew.type = M_TRANS_EV;
							if(tra.nm_trans_ref != UNSET) enew.type = NM_TRANS_EV;
							c_after = tra.f;
						}
						break;
					
					default: c_after = UNSET; tr_gl = UNSET; emsg("SHould not be move"); break;
					}
					
					auto t = time_list[k];
					enew.tr_gl = tr_gl;
					enew.c_after = c_after;
					enew.t = t;
					
					rate_mean.update_prob_try_Rint(t_new,t,tr_gl,cl,c_new,prob_try,Rint,ind,MULT_MOD);
					t_new = t;
					
					ev_new.push_back(enew);
					c_new = c_after;
				}
				else{
					rate_mean.update_prob_try_Rint(t_new,tmax,UNSET,cl,c_new,prob_try,Rint,ind,MULT_MOD);
					t_new = tmax;
				}
				k++; if(k < K) t_event_new = time_list[k]; else t_event_new = LARGE;
			}
			break;
			
		case OLD_EVENT_LOC:                   // Removes old event
			{
				const auto &eve = event[e];
				auto t = eve.t;
				auto tr_gl = UNSET;
				
				if(r < R && rem[r] == e){
					if(eve.observed) return 0;
					tr_gl = eve.tr_gl;
					r++;
				}
				else{
					auto ev_add = eve;
					
					if(sp.correct_ev(c_new,ev_add) == true) emsg("Inconsistent");

					if(c_new != UNSET && t > t_new){
						rate_mean.update_prob_try_Rint(t_new,t,UNSET,cl,c_new,prob_try,Rint,ind,MULT_MOD);
						t_new = t;
					}
					
					ev_new.push_back(ev_add);
					c_new = ev_add.c_after; 
				}
				
				if(c_old != UNSET && t > t_old){
					rate_mean.update_prob_try_Rint(t_old,t,tr_gl,cl,c_old,prob_try,Rint,ind,DIV_MOD);
					t_old = t;
				}
				
				c_old = eve.c_after;
				e++; if(e < E) t_event_old = event[e].t; else t_event_old = LARGE;
			}
			break;
		
		case OBS_LOC:                         // Processes an observation
			{
				const auto &ob = obs[m];
				
				if(c_old != c_new && ob.cl == cl){
					switch(ob.type){
					case OBS_COMP_EV:
						{
							auto c_exact = ob.c_exact;
							if(c_exact != UNSET){
								if(sp.comp_gl[c_new].cla_comp[cl] != c_exact) return 0;
								if(sp.comp_gl[c_old].cla_comp[cl] != c_exact){
									prob_try *= (1+PROB_MOD)/PROB_MOD;
								}
							}
							else{
								auto prob_new = obs_eqn_value[ob.obs_eqn_ref[sp.comp_gl[c_new].cla_comp[cl]]];
								if(prob_new == 0) return 0;
								auto prob_old = obs_eqn_value[ob.obs_eqn_ref[sp.comp_gl[c_old].cla_comp[cl]]];
								prob_try *= (prob_new+PROB_MOD)/(prob_old+PROB_MOD);
							}
						}
						break;
					
					case OBS_TRANS_EV: case OBS_SOURCE_EV: case OBS_SINK_EV:
						break;
					
					case OBS_TEST_EV:
						{
							const auto &om = sp.source[ob.so].obs_model;
							auto tsens_old = om.diag_test_sens.comp[sp.comp_gl[c_old].cla_comp[cl]];
							auto tsens_new = om.diag_test_sens.comp[sp.comp_gl[c_new].cla_comp[cl]];
							
							if(tsens_old != tsens_new){
								auto Se = obs_eqn_value[ob.Se_obs_eqn_ref];
								auto Sp = obs_eqn_value[ob.Sp_obs_eqn_ref];
								
								if(ob.test_res == true){
									if(tsens_new == true) prob_try *= (Se+LOG_THRESH)/(1-Sp+LOG_THRESH);
									else prob_try *= (1-Sp+LOG_THRESH)/(Se+LOG_THRESH);
								}
								else{
									if(tsens_new == true) prob_try *= (1-Se+LOG_THRESH)/(Sp+LOG_THRESH);
									else prob_try *= (Sp+LOG_THRESH)/(1-Se+LOG_THRESH);
								}
								
								if(prob_try == 0) return 0;
							}
						}
						break;
					}
				}
				
				m++; if(m < M) t_obs = obs[m].t; else t_obs = LARGE;
			}		
			break;
			
		case DO_NOTHING_LOC: emsg("SHould not be here"); break;
		}
	}while(true);
	
	// Accounts for the final part of the integral
	if(t_old != tmax && c_old != UNSET){ 
		rate_mean.update_prob_try_Rint(t_old,tmax,UNSET,cl,c_old,prob_try,Rint,ind,DIV_MOD);
	}
	
	if(t_new != tmax && c_new != UNSET){
		rate_mean.update_prob_try_Rint(t_new,tmax,UNSET,cl,c_new,prob_try,Rint,ind,MULT_MOD);		
	}

	prob_try *= exp(Rint);
	
	// Adds remaining events
	while(e < E){ ev_new.push_back(event[e]); e++;} 

	if(dE > 0){
		for(auto j = 0; j < dE; j++) timefac *= dt;
		timefac /= factori[dE];
	}
	else{
		if(dE < 0){
			for(auto j = 0; j < -dE; j++) timefac /= dt;
			timefac *= factori[-dE];
		}
	}
	
	prob_try *= timefac;
	
	if(prob_try < 0.05) prob_try = 0.05;
	if(prob_try > 20) prob_try = 20;			

	return prob_try;
}
	
	
/// Initialises the rate store
void StateSpecies::rate_store_init()
{
	rate_mean.first = true;
	
	auto M = sp.markov_eqn.size();
	auto N = sp.nm_trans.size();
	
	// Sets up Markovian transitions
	auto &tra_rate = rate_mean.tra_rate;
	for(auto i = 0u; i < M; i++){
		TraRate trate;
		trate.markov = true;
		trate.bp_on = false;
		const auto &mar_eqn = sp.markov_eqn[i];
		for(auto ie : mar_eqn.ind_eff_mult){
			IndFacRate ifr;
			if(mar_eqn.rate) ifr.type = IND_EFF_MULT; 
			else ifr.type = IND_EFF_DIV; 
			ifr.e = ie;
			trate.ind_fac_rate.push_back(ifr);
		}
		
		for(auto ie : mar_eqn.fix_eff_mult){
			IndFacRate ifr; 
			if(mar_eqn.rate) ifr.type = FIX_EFF_MULT; 
			else ifr.type = FIX_EFF_DIV; 
			ifr.e = ie;
			trate.ind_fac_rate.push_back(ifr);
		}
		trate.value.resize(T);
		trate.value_raw.resize(T,0);
		trate.integral.resize(T+1);
		
		tra_rate.push_back(trate);
	}
	
	// Sets up non-Markovian transitions
	for(auto i = 0u; i < N; i++){
		const auto &nmt = sp.nm_trans[i];
		TraRate trate;
		trate.markov = false;
		trate.bp_on = false; if(nmt.bp_eq != UNSET) trate.bp_on = true;
		trate.ind_fac_rate = nmt.ind_fac_rate;
		trate.value.resize(T,0);
		trate.value_raw.resize(T,0);
		trate.integral.resize(T+1);
		
		if(trate.bp_on){
			trate.bp.resize(T,0);
			trate.bp_raw.resize(T,0);
		}		
		tra_rate.push_back(trate);
	}
	
	// Sets if there is any individual variation
	for(auto &trate : rate_mean.tra_rate){
		trate.ind_variation = false;
		if(trate.ind_fac_rate.size() > 0) trate.ind_variation = true;
	}
			
	// Sets up tramean_ref
	for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
		const auto &tra = sp.tra_gl[tr];
		if(tra.markov_eqn_ref != UNSET){
			rate_mean.tramean_ref.push_back(tra.markov_eqn_ref);
		}
		else{
			rate_mean.tramean_ref.push_back(M+tra.nm_trans_ref);
		}
	}
	
	// Sets up compR
	auto &compR = rate_mean.compR;
	auto &compR_ref = rate_mean.compR_ref;
	compR_ref.resize(sp.ncla);
	for(auto cl = 0u; cl < sp.ncla; cl++){
		compR_ref[cl].resize(sp.comp_gl.size());
		for(auto c = 0u; c < sp.comp_gl.size(); c++){
			const auto &tlg = sp.comp_gl[c].tra_leave_group[cl];
	
			vector <unsigned int> list;
			for(auto tr : tlg.tr_list){
				const auto &tra = sp.tra_gl[tr];
				if(tra.markov_eqn_ref != UNSET) list.push_back(tra.markov_eqn_ref);
				else list.push_back(M+tra.nm_trans_ref); 
			}
			
			auto m = 0u; while(m < compR.size() && !equal_vec(compR[m].list,list)) m++;
			
			compR_ref[cl][c] = m;
			
			if(m == compR.size()){
				CompR compr;
				compr.list = list;
				
				switch(list.size()){ 
				case 0: 
					compr.update = NO_UP;
					break;
					
				case 1:
					compr.update = SINGLE_BRANCH_UP;
					for(auto m : list){
						if(rate_mean.tra_rate[m].ind_variation == true){
							compr.update = SINGLE_CALCALL_UP;
						}
					}
					break;
					
				default:
					if(tra_rate[list[0]].markov){
						compr.update = MARKOV_COMPR_UP;
						for(auto m : list){
							if(rate_mean.tra_rate[m].ind_variation == true){
								compr.update = MARKOV_CALCALL_UP;
							}
						}
					}
					else{
						compr.update = NM_CALCALL_UP;
					}
				}
				compr.integral.resize(T+1);
				compR.push_back(compr);
			}
		}
	}
	
	if(false){
		cout << "TR RATE" << endl;
		for(auto e = 0u; e < rate_mean.tra_rate.size(); e++){
			cout << e << " ";
			if(e < M) cout << eqn[sp.markov_eqn[e].eqn_ref].te;
			else cout << eqn[sp.nm_trans[e-M].dist_param_eq_ref[0]].te;
			
			cout << "ind vari: " << rate_mean.tra_rate[e].ind_variation;
				
			cout << " rate mean" << endl;
		}
		
		cout << "COMPR" << endl;
		for(auto e = 0u; e < rate_mean.compR.size(); e++){
			const auto &compr = rate_mean.compR[e];
			cout << e << " ";
			for(auto j : compr.list) cout << j << ",";
			
			switch(compr.update){
			case SINGLE_BRANCH_UP: cout << "SINGLE_BRANCH_UP"; break;
			case SINGLE_CALCALL_UP: cout << "SINGLE_CALCALL_UP"; break;
			case MARKOV_COMPR_UP: cout << "MARKOV_COMPR_UP"; break;
			case MARKOV_CALCALL_UP: cout << "MARKOV_CALCALL_UP"; break;
			case NM_CALCALL_UP: cout << "NM_CALCALL_UP"; break;
			case NO_UP: cout << "NO_UP"; break;
			}
			cout << endl;
		}
		emsg("P");
	}
}


/// Updates information which gives the rates of transition (for local proposals)
void StateSpecies::update_rate_mean(const vector < vector <double> > &popnum_t)
{
	auto fac = UPDATE_RATE_FAC;
	if(rate_mean.first == true){ fac = 0; rate_mean.first = false;}
	auto omf = 1-fac;
	
	for(auto e = 0u; e < sp.markov_eqn.size(); e++){
		const auto &mev = markov_eqn_vari[e];
		const auto &div = mev.div;
		auto &value = rate_mean.tra_rate[e].value_raw;
		if(mev.time_vari){
			for(auto ti = 0u; ti < T; ti++){
				value[ti] = fac*value[ti] + omf*div[ti].value;
			}
		}
		else{
			auto val = fac*value[0] + omf*div[0].value;
			for(auto ti = 0u; ti < T; ti++) value[ti] = val;
		}
	}
	
	auto M = sp.markov_eqn.size();
	
	vector < vector <double> > bp_store;
	auto nm_rate = sp.calc_nm_rate(true,param_val,spline_val,popnum_t,eqn,bp_store);
	
	for(auto m = 0u; m < sp.nm_trans.size(); m++){
		auto &trate = rate_mean.tra_rate[M+m];
		auto &value = trate.value_raw;
		for(auto ti = 0u; ti < T; ti++){
			value[ti] = fac*value[ti] + omf*nm_rate[m][ti];
		}  
		if(trate.bp_on){
			auto &bp = trate.bp_raw;
			for(auto ti = 0u; ti < T; ti++){
				bp[ti] = fac*bp[ti] + omf*bp_store[m][ti];
			}  
		}
	}
	
	auto dt = details.dt;
	
	// Generates a smoothed version of rates (this helps the proposals)
	for(auto tr = 0u; tr < rate_mean.tra_rate.size(); tr++){
		auto &trate = rate_mean.tra_rate[tr];
			
		const auto &value_old = trate.value_raw;
		auto mean = 0.0;
		for(auto ti = 0u; ti < T; ti++) mean += value_old[ti];
		mean /= T;
		
		auto &value_new = trate.value;
		for(auto ti = 0u; ti < T; ti++){
			value_new[ti] = LOCAL_SMOOTH*mean + (1-LOCAL_SMOOTH)*value_old[ti];
		}
		
		auto &integral = trate.integral;
		
		auto sum = 0.0;
		for(auto ti = 0u; ti < T; ti++){ 
			integral[ti] = sum;
			sum += value_new[ti]*dt;
		}
		integral[T] = sum;
		
		if(false){
			for(auto ti = 0u; ti < T; ti++){
				cout << ti << " " << value_old[ti] <<" " << value_new[ti] << " " << mean << " jj" << endl;
			}
		}
		
		const auto &bp_old = trate.bp_raw;
		if(bp_old.size() > 0){
			auto mean = 0.0;
			for(auto ti = 0u; ti < T; ti++) mean += bp_old[ti];
			mean /= T;
			
			auto &bp_new = trate.bp;
			for(auto ti = 0u; ti < T; ti++){
				bp_new[ti] = LOCAL_SMOOTH*mean + (1-LOCAL_SMOOTH)*bp_old[ti];
			}
		}
	}
	
	vector <double> value(T);
	
	auto &compR = rate_mean.compR;
	for(auto e = 0u; e < compR.size(); e++){
		auto &compr = compR[e]; 
		auto &integral = compr.integral;
			
		switch(compr.update){
		case SINGLE_BRANCH_UP: case MARKOV_COMPR_UP: 
			{	
				for(auto ti = 0u; ti < T; ti++) integral[ti] = 0;
				
				for(auto m : compr.list){
					const auto &integral_add = rate_mean.tra_rate[m].integral;
					for(auto ti = 0u; ti <= T; ti++) integral[ti] += integral_add[ti];
				}
			}
			break;
		
		default: break;
		}
		
		if(false){
			for(auto ti = 0u; ti <= T; ti++){
				cout << ti << " " << integral[ti] << "Integra" << endl;
			}
			emsg("do");
		}
	}
}


// Cnstructor
RatePosteriorMean::RatePosteriorMean(const Details &details)
{
	t_start = details.t_start; 
	dt = details.dt;
}


/// Gets the value from the vector 
double RatePosteriorMean::get_value(const vector <double> &vec, double v) const 
{
	auto i = (unsigned int) v;
	if(i+1 >= vec.size()) emsg("Vec out of range");
	auto f= v-i;
	return vec[i]*(1-f) + vec[i+1]*f;
}


/// Updates prob_try and Rint based on a transition
void RatePosteriorMean::update_prob_try_Rint(double t1, double t2, unsigned int tr_gl, unsigned int cl, unsigned int c, double &prob_try, double &Rint, const Individual &ind, Modify mod) const
{	
	const auto &compr = compR[compR_ref[cl][c]];
	if(compr.update == NM_CALCALL_UP) return;
	
	auto v1 = ALMOST_ONE*(t1-t_start)/dt;
	auto v2 = ALMOST_ONE*(t2-t_start)/dt;
		
	switch(compr.update){
	case SINGLE_BRANCH_UP: case MARKOV_COMPR_UP: 
		{
			const auto &integral = compr.integral;
			auto val = get_value(integral,v2) - get_value(integral,v1);
			if(mod == MULT_MOD) Rint -= val;
			else Rint += val;
		}
		break;
		
	case MARKOV_CALCALL_UP: case SINGLE_CALCALL_UP: 
		{
			auto val_total = 0.0;
			for(auto m : compr.list){
				const auto &trate = tra_rate[m];
				const auto &integral = trate.integral;
				auto val = get_value(integral,v2) - get_value(integral,v1);
			
				if(trate.ind_variation){
					for(const auto &ifr : trate.ind_fac_rate){
						switch(ifr.type){
						case IND_EFF_MULT: val *= ind.exp_ie[ifr.e]; break;
						case IND_EFF_DIV: val /= ind.exp_ie[ifr.e]; break;
						case FIX_EFF_MULT: val *= ind.exp_fe[ifr.e]; break;
						case FIX_EFF_DIV: val /= ind.exp_fe[ifr.e]; break;
						}
					}
				}	
				val_total += val;
			}
			if(mod == MULT_MOD) Rint -= val_total;
			else Rint += val_total;
		}
		break;
	
	case NM_CALCALL_UP: 
		emsg("Should not be here");
		break;
		
	case NO_UP: break;
	}
	
	if(tr_gl != UNSET){
		const auto &trate = tra_rate[tramean_ref[tr_gl]];
		
		auto ti = (unsigned int)(ALMOST_ONE*(t2-t_start)/dt);
		auto rate = trate.value[ti];

		if(trate.ind_variation){
			for(const auto &ifr : trate.ind_fac_rate){
				switch(ifr.type){
				case IND_EFF_MULT: rate *= ind.exp_ie[ifr.e]; break;
				case IND_EFF_DIV: rate /= ind.exp_ie[ifr.e]; break;
				case FIX_EFF_MULT: rate *= ind.exp_fe[ifr.e]; break;
				case FIX_EFF_DIV: rate /= ind.exp_fe[ifr.e]; break;
				}
			}
		}
		
		if(mod == MULT_MOD) prob_try *= (rate+TINY);
		else prob_try /= (rate+TINY);
	}
}


/// Works out if the reverse probability is correct
void StateSpecies::swap_check_rev_prob_try(double prob_try, double timefac, const LocalIndChange &lich, const vector <LocalIndChange> &add, const vector <Event> &ev_store) const 
{
	const auto &claa = sp.cla[lich.cl];
	const auto &swr = claa.swap_rep[lich.swap_rep_ref];
	const auto &sw = claa.swap[swr.swap_ref];

	// Looks for the reverse transition
	auto k_sel = UNSET;
	for(auto k = 0u; k < add.size(); k++){
		const auto &lich_new = add[k];
		//print_local_ind_change(lich_new);
	
		if(lich_new.cl == lich.cl && lich_new.tmin == lich.tmin && 
			lich_new.tmax == lich.tmax){
			if(tr_swap_same(lich_new.tr_swap,sw.start)){
				if(k_sel == UNSET) k_sel = k;
				else emsg("More than one");
			}
		}
	}
	if(k_sel == UNSET){
		print_local_ind_change(lich);
		for(auto k = 0u; k < add.size(); k++){
			print_local_ind_change(add[k]);
		}
		emsg("Cannot find reverese");
	}
	
	Event enew;
	vector <Event> ev_new;
	auto nind_obs = sp.nindividual_obs;
						
	auto timefac_rev = 1.0;

	auto prob_try_rev = create_local_change(timefac_rev,add[k_sel],ev_new,enew,nind_obs,ev_store,LOCAL_REVERSE);

	if(prob_try_rev != 0){
		if(dif(prob_try,1.0/prob_try_rev,dif_thresh)){
			emsg("prob_try dif");
		}
		
		if(dif(timefac,1.0/timefac_rev,dif_thresh)){
			cout << "Change:" << endl;
			print_local_ind_change(lich);
		
			cout << "Reverse:" << endl;
			print_local_ind_change(add[k_sel]);
			emsg("time_fac dif");
		}
	}
}


/// Determines if tr_swap is the same
bool StateSpecies::tr_swap_same(const vector <TrSwap> &tr_swap, const vector <TrSwap> &tr_swap2) const
{
	if(tr_swap.size() != tr_swap2.size()) return false;
	for(auto i = 0u; i < tr_swap.size(); i++){
		if(tr_swap[i].type != tr_swap2[i].type) return false;
		if(tr_swap[i].trc != tr_swap2[i].trc) return false;
	}
	return true;
}

// Updates an individual within the state species

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
 
using namespace std;

#include "state_species.hh"
#include "utils.hh"
#include "matrix.hh"

/// Works out how likelihood for nm transitions changes based on sequence alteration
// In the case of non-Markovian transitions m, ti, index and e_origin are automatically set  
// In the case of Markovian transitions m, index automatically set  
vector <unsigned int>  StateSpecies::update_ind(unsigned int i, vector <Event> &event_new, const vector < vector <double> > &popnum_t, vector <PopUpdate> &pop_update, Like &like_ch)
{
	auto &ind = individual[i];
	const auto &indi = sp.individual[i];
	
	BackInd bi; 		                                  // Stores individual (so can be recovered)
	bi.i = i; 
	bi.ev = ind.ev; 
	bi.incomp_ref = ind.incomp_ref;
	
	auto &identical_old = bi.identical_old;
	auto &identical_new = bi.identical_new;
	auto &remove_old = bi.remove_old;
	auto &add_new = bi.add_new;
	
	ind.ev = event_new;                               // Copies new event sequence into individual
	
	auto &ev_old = bi.ev; 
	auto &ev_new = ind.ev;     
	
	set_m_ti_origin(ev_new);
	
	if(false){
		cout << endl << endl;
		cout << "Event before:" << endl;
		print_event(bi.ev);
		cout << "Event after:" << endl;
		print_event(ev_new);
	}
	
	vector <unsigned int> add_ev;                     // Makes list of events to be added later
	
	if(false){
		cout << "event:" << endl; print_event(ind.ev);
		cout << "event new:" << endl; print_event(ev_new);
	}
	
	// If individual unobserved then sets obs to fake individual at end
	auto iobs = i; 
	auto observed = true;
	if(iobs >= sp.nindividual_obs){
		iobs = sp.nindividual_obs;
		observed = false;
	}
	
	const auto &obs = sp.individual[iobs].obs;
	
	const auto &popd = sp.pop_data;
	
	auto e_old = 0u;
	auto e_old_max = ev_old.size();
	auto e_new = 0u;
	auto e_new_max = ev_new.size();
	
	auto m = 0u; // Indexes observations
	auto k = 0u; // Indexes pop_data
	
	auto dLi_obs = 0.0, dLi_obs_ind = 0.0;
	
	auto c_old = UNSET; 
	auto c_new = UNSET;
	
	auto c_enter_old = UNSET;
	auto c_enter_new = UNSET;
	
	double t_old, t_new;
	
	auto t = 0.0;
	auto ti = 0u; 
	do{
		if(e_old < e_old_max) t_old = ev_old[e_old].tdiv; else t_old = T;
		if(e_new < e_new_max) t_new = ev_new[e_new].tdiv; else t_new = T;
		
		auto t_next = t_old; if(t_new < t_next) t_next = t_new;
		
		auto ti_next = get_ti_next(t_next,details);
	
		while(k < popd.size() && popd[k].tdiv <= t_next){   // Accounts for population data
			if(c_new != c_old) dLi_obs += update_pop_data(k,c_old,c_new);
			k++;
		}

		while(m < obs.size() && obs[m].tdiv <= t_next){
			const auto &ob = obs[m];
			
			switch(ob.type){
			case OBS_TRANS_EV:                              // Checks to see if agrees with observation
			case OBS_SOURCE_EV:
			case OBS_SINK_EV:  
				if(obs[m].tdiv == t_next){
					if(t_old != t_new){
						if(t_next == t_new) dLi_obs_ind -= LI_WRONG;
						else dLi_obs_ind += LI_WRONG;
					}
					
					// Updates probability asociated with observation of event
					auto tr_old = UNSET, tr_new = UNSET;
					if(t_next == t_old){
						const auto &ev = ev_old[e_old];
						if(ev.observed == true && (ev.type == M_TRANS_EV || ev.type == NM_TRANS_EV)){
							tr_old = ev.tr_gl;
						}
					}
					
					if(t_next == t_new){
						const auto &ev = ev_new[e_new];
						if(ev.observed == true && (ev.type == M_TRANS_EV || ev.type == NM_TRANS_EV)){
							tr_new = ev.tr_gl;
						}
					}
					
					if(tr_old != tr_new){
						if(tr_old != UNSET) dLi_obs_ind -= like_trans_prob_obs(tr_old,obs[m]);	
						if(tr_new != UNSET) dLi_obs_ind += like_trans_prob_obs(tr_new,obs[m]);
					}
				}
				break;
				
			case OBS_COMP_EV:
				{
					if(c_new != c_old){
						auto dLi = like_comp_obs(c_new,ob) - like_comp_obs(c_old,ob);
						dLi_obs_ind += dLi;
					}
				}
				break;
				
			case OBS_TEST_EV:
				{
					if(c_new != c_old){
						dLi_obs_ind += like_diag_test(c_new,ob) - like_diag_test(c_old,ob);
					}
				}
				break;
			}
			
			m++;
		}
		
		if(c_new != c_old && ti_next != ti){
			PopUpdate pup; 
			pup.ti = ti;
			pup.ti_next = ti_next;
			pup.c_old = c_old;
			pup.c_new = c_new;
			pop_update.push_back(pup);
		}
				
		// Updates the global rate part of Markov equations
		if(c_new != c_old && t != t_next){
			if(c_old != UNSET) update_indfac_int(c_old,t,t_next,ind,-1,like_ch.markov);
			if(c_new != UNSET) update_indfac_int(c_new,t,t_next,ind,1,like_ch.markov);
		}
		
		if(t_next == T) break;
		
		t = t_next; ti = ti_next;
		
		// Works out if events are identical
		auto identical = false;
		
		if(t_old == t_new){ 
			auto &eve_old = ev_old[e_old];
			auto &eve_new = ev_new[e_new];
			
			switch(eve_old.type){
			case M_TRANS_EV:
				if(eve_new.type == M_TRANS_EV && eve_old.m == eve_new.m){
					identical = true;
					eve_new.index = eve_old.index;
					if(e_old != e_new){
						markov_eqn_vari[eve_old.m].div[eve_old.ti].ind_trans[eve_old.index].index = e_new;
					}
				}
				break;
				
			case NM_TRANS_EV:	
				if(eve_new.type == NM_TRANS_EV && eve_old.m == eve_new.m){
					if(ev_old[eve_old.e_origin].tdiv == ev_new[eve_new.e_origin].tdiv){
						identical = true;
						
						eve_new.index = eve_old.index;
						eve_new.Li = eve_old.Li;
						eve_new.Li_bp = eve_old.Li_bp;
						
						// Updates nm_trans_ev_ref if neccesary
						if(e_old != e_new || eve_old.e_origin != eve_new.e_origin){ 
							auto &nmter = nm_trans_ev_ref[eve_new.m][eve_new.ti][eve_new.index];		
							if(nmter.i != i) emsg("Wrong nmter");
							nmter.e_end = e_new;
						}		
					}
				}
				break;
			
			default:
				break;
			}
			
			if(identical){
				identical_old.push_back(e_old);
				identical_new.push_back(e_new);
			
				c_old = eve_old.c_after; e_old++;
				c_new = eve_new.c_after; e_new++;
				
				auto tr_old = eve_old.tr_gl;
				auto tr_new = eve_new.tr_gl;
				if(tr_old != tr_new){
					if(sp.obs_trans_exist){
						alter_obs_trans_eqn(tr_old,t,-1,eve_old.observed,obs,dLi_obs);
						alter_obs_trans_eqn(tr_new,t,1,eve_new.observed,obs,dLi_obs);
					}
					
					if(sp.pop_trans_data_exist){
						dLi_obs += update_pop_trans_data(t,tr_old,-1);
						dLi_obs += update_pop_trans_data(t,tr_new,1);
					}			
				}
			}
		}
		
		if(!identical){
			if(t == t_old){                               // Removes old event
				auto &eve_old = ev_old[e_old];
				
				remove_old.push_back(e_old);
				remove_event(eve_old,ind,like_ch,i,ev_old); 
	
				c_old = eve_old.c_after; e_old++;

				switch(eve_old.type){
				case M_TRANS_EV: case NM_TRANS_EV:
					if(sp.obs_trans_exist){
						alter_obs_trans_eqn(eve_old.tr_gl,t,-1,eve_old.observed,obs,dLi_obs);
					}
				
					if(sp.pop_trans_data_exist){
						dLi_obs += update_pop_trans_data(t_next,eve_old.tr_gl,-1);
					}
					break;
				
				case ENTER_EV:
					if(observed == true && !indi.init_c_set){
						dLi_obs_ind -= like_ind_init_obs(c_old,indi.enter_ref);
					}
					c_enter_old = c_old;
					break;
					
				default: break;
				}
			}
			else{	                                        // Adds new event
				const auto &eve_new = ev_new[e_new];
				
				add_new.push_back(e_new);
				add_ev.push_back(e_new);
				
				c_new = eve_new.c_after; e_new++;
				
				switch(eve_new.type){
				case M_TRANS_EV: case NM_TRANS_EV:
					if(sp.obs_trans_exist){
						alter_obs_trans_eqn(eve_new.tr_gl,t,1,eve_new.observed,obs,dLi_obs);
					}
					
					if(sp.pop_trans_data_exist){
						dLi_obs += update_pop_trans_data(t_next,eve_new.tr_gl,1);
					}
					break;
				
				case ENTER_EV:
					if(observed == true && !indi.init_c_set){
						dLi_obs_ind += like_ind_init_obs(c_new,indi.enter_ref);
					}
					c_enter_new = c_new;
					break;
					
				default: break;
				}
			}
		}
	}while(true);

	// Updates the initial condition (if neccesary)
	if(c_enter_old != c_enter_new){
		switch(sp.init_cond.type){
		case INIT_POP_DIST: likelihood_init_cond_change(c_enter_old,c_enter_new,like_ch); break;
		case INIT_POP_NONE: break;
		case INIT_POP_FIXED: emsg("Cannot change pop fixed"); break;
		}
	}
	
	if(observed){
		back.push_back(Back(LI_OBS_IND,i,dLi_obs_ind));
		Li_obs_ind[i] += dLi_obs_ind;
	}
	dLi_obs += dLi_obs_ind;
	
	like_ch.obs += dLi_obs;
	
	// Updates incomp_ref
	set_incomp_ref(i,popnum_t,like_ch.nm_trans);
	
	back_ind.push_back(bi);
	
	return add_ev;
}


/// Sets ev.m, ev.ti, ev.e_origin for new sequenece
void StateSpecies::set_m_ti_origin(vector <Event> &ev_new) const
{
	if(sp.nm_flag){
		vector <unsigned int> e_vec(sp.ncla,0);
		for(auto e = 0u; e < ev_new.size(); e++){
			auto &ev = ev_new[e];
			switch(ev.type){
			case M_TRANS_EV:
				{
					const auto &tra = sp.tra_gl[ev.tr_gl];
					auto m = tra.markov_eqn_ref;
					ev.m = m;
					if(sp.markov_eqn[m].time_vari) ev.ti = get_ti(ev.tdiv);
					else ev.ti = 0;
					e_vec[tra.cl] = e;
				}
				break;
			
			case NM_TRANS_EV:
				{
					auto tr = ev.tr_gl;
					const auto &tra = sp.tra_gl[tr];
					
					auto cl = tra.cl;
					ev.e_origin = e_vec[cl];
					const auto &ev_orig = ev_new[ev.e_origin];
					
					ev.ti = get_ti(ev_orig.tdiv);
					e_vec[cl] = e;
					
					if(tra.i != ev_orig.c_after){             // Accounts for intermediate transition
						tr = sp.tr_trans(tr,ev_orig.c_after);
						if(tr == UNSET) emsg("trg unset2");
						ev.m = sp.tra_gl[tr].nm_trans_ref;
					}
					else ev.m = tra.nm_trans_ref;
					
					ev.index = UNSET;
				}
				break;
				
			default: break;
			}
		}
	}
	else{
		for(auto e = 0u; e < ev_new.size(); e++){
			auto &ev = ev_new[e];
			if(ev.type == M_TRANS_EV){
				auto m = sp.tra_gl[ev.tr_gl].markov_eqn_ref;
				ev.m = m;
				if(sp.markov_eqn[m].time_vari) ev.ti = get_ti(ev.tdiv);
				else ev.ti = 0;
				ev.index = UNSET;
			}
		}
	}
}
	

/// Removes event from old event sequence 
void StateSpecies::remove_event(Event &ev, const Individual &ind, Like &like_ch, unsigned int i, vector <Event> &event_old)
{
	switch(ev.type){
	case NM_TRANS_EV:
		{
			auto &nmter = nm_trans_ev_ref[ev.m][ev.ti];

			like_ch.nm_trans -= ev.Li + ev.Li_bp;
		
			auto nmt_ind = ev.index;
			
			if(nmt_ind+1 < nmter.size()){
				const auto &eref = nmter[nmter.size()-1];
				nmter[nmt_ind] = eref;	
				
				auto ii = eref.i;
				if(ii != i){
					individual[ii].ev[eref.e_end].index = nmt_ind;
				}
				else{
					auto index = eref.e_end;
					if(event_old[index].tdiv > ev.tdiv){	
						event_old[index].index = nmt_ind;
					}
					else{
						individual[ii].ev[index].index = nmt_ind;
					}
				}
			}
			nmter.pop_back();	
		}
		break;
		
	case M_TRANS_EV:
		{
			double Li;
	
			auto e = ev.m;
			auto ti = ev.ti;
			auto &it = markov_eqn_vari[e].div[ti].ind_trans;
			auto k = ev.index;
		
			if(k >= it.size()) emsg(" out of range a");
			
			if(k != it.size()-1){
				auto &it_shift = it[it.size()-1];
				it[k] = it_shift;
				auto ii = it_shift.i;
	
				if(ii != i){
					individual[ii].ev[it_shift.index].index = k;
				}
				else{
					auto index = it_shift.index;
					if(event_old[index].tdiv > ev.tdiv){	
						event_old[index].index = k;
					}
					else{
						individual[ii].ev[index].index = k;
					}
				}
			}
			it.pop_back();
		
			auto &me = sp.markov_eqn[e];
			
			// If an infection transition then works out specific probability on individual
			if(me.infection_trans){
				const auto &eq = eqn[me.eqn_ref];
				const auto &lin = eq.linearise; if(!lin.on) emsg("Linearisation should be on");
		
				const auto &iif = ev.ind_inf_from;
				auto p = iif.p;
				double va;
				switch(p){
				case OUTSIDE_INF: va = eq.calculate_calculation(lin.no_pop_calc,ti,param_val,spline_val); break;
				case ENTER_INF: va = UNSET; emsg("Should not be ENTER_INF"); break;
				default: va = eq.calculate_calculation(lin.pop_grad_calc[iif.pref],ti,param_val,spline_val); break;
				}
				
				if(me.ind_variation) va *= get_indfac(ind,me);
				Li = log(va*iif.w);
			}
			else{			
				auto val = markov_eqn_vari[e].div[ti].value;
				
				if(me.ind_variation) Li = log(get_indfac(ind,me)*(val+LOG_THRESH));
				else Li = log(val+LOG_THRESH);
			}
			
			back.push_back(Back(REMOVE_LI_MARKOV,e,ti,Li));
			
			Li_markov[e][ti] -= Li;
			
			like_ch.markov -= Li;
		}
		break;
		
	default: break;
	}
}


/// Adds event from individual to reference (either 
void StateSpecies::add_event_ref(unsigned int i, unsigned int ee, const vector < vector <double> > &popnum_t, Like &like_ch)
{
	auto &ind = individual[i];
	auto &event = ind.ev;
	auto &ev = event[ee];
	auto dt = details.dt;
	
	switch(ev.type){
	case NM_TRANS_EV:
		{
			const auto &ev_orig = event[ev.e_origin];
			auto t = ev_orig.tdiv;		
			auto ti = ev.ti;
			auto m = ev.m;
			auto dtdiv = ev.tdiv - t; 
			
			if(dtdiv <= 0){
				emsg("zero time1");
			}
			const auto &tra = sp.tra_gl[ev.tr_gl];	
			const auto &nmt = sp.nm_trans[m];
			const auto &ref = nmt.dist_param_eq_ref;
					
			vector <double> ref_val(ref.size());

			for(auto i = 0u; i < ref.size(); i++){
				const auto &eq = eqn[ref[i]];
				if(eq.ind_eff == true) ref_val[i] = eq.calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
				else{
					if(eq.time_vari) ref_val[i] = eq.calculate(ti,popnum_t[ti],param_val,spline_val);
					else ref_val[i] = eq.calculate_param_only(param_val);
				}			
			}
					 
			auto dLi = nm_trans_like(nmt.type,dtdiv,dt,ref_val);
			
			auto dLi_bp = 0.0;
			if(tra.branch){                               // Accounts for branching probability 
				auto bq_eq = nmt.bp_eq;
				double bp;
				if(bq_eq == BP_FROM_OTHERS){
					bp = 1;
					for(auto e : nmt.bp_other_eq){
						bp -= eqn[e].calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
					}
				}
				else{
					const auto &eq = eqn[nmt.bp_eq];
					bp = eq.calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
					
					if(nmt.all_branches){
						auto div = 0.0;
						for(auto e: nmt.bp_all_eq){
							div += eqn[e].calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
						}
						bp /= div;
					}			
				}
				
				dLi_bp = get_log_zero_one(bp);
			}	
			
			NMEventRef eref;
			eref.i = i;
			eref.e_end = ee;
			
			auto &nmter = nm_trans_ev_ref[m][ti];
			
			ev.Li = dLi; ev.Li_bp = dLi_bp;
			ev.index = nmter.size();
			
			nmter.push_back(eref);
			
			like_ch.nm_trans += dLi+dLi_bp;
		}
		break;

	case M_TRANS_EV:
		{
			double Li;
	
			// Add individual on to ind_trans and sets .markov_eq_ref for event
			auto e = ev.m;
			auto ti = ev.ti;
			
			auto &mev = markov_eqn_vari[e];
			auto &it = mev.div[ti].ind_trans;	
			ev.index = it.size();
				
			DivIndRef dir; dir.i = i; dir.index = ee;
			it.push_back(dir);
		
			auto &me = sp.markov_eqn[e];
			
			// If an infection transition then works out specific probability on individual
			if(me.infection_trans){ 
				const auto &eq = eqn[me.eqn_ref];
				const auto &lin = eq.linearise; if(!lin.on) emsg("Linearisation should be on");
		
				const auto &iif = ev.ind_inf_from;
				auto p = iif.p;
				double va;
				switch(p){
				case OUTSIDE_INF:
					va = eq.calculate_calculation(lin.no_pop_calc,ti,param_val,spline_val);
					break;
				case ENTER_INF:
					va = UNSET; emsg("Should not be ENTER_INFa"); 
					break;
				default: 
					va = eq.calculate_calculation(lin.pop_grad_calc[iif.pref],ti,param_val,spline_val); 
					break;
				}
				
				if(me.ind_variation) va *= get_indfac(ind,me);
				va *= dt;
				Li = log(va*iif.w);
			}
			else{
				auto val = markov_eqn_vari[e].div[ti].value;	
				if(me.ind_variation) Li = log(get_indfac(ind,me)*(val+LOG_THRESH));
				else Li = log(val+LOG_THRESH);
			}
			
			back.push_back(Back(ADD_LI_MARKOV,e,ti,Li));
			
			Li_markov[e][ti] += Li;
			like_ch.markov += Li;
		}
		break;
		
	default: break;
	}
}


/// Removes event references from individuals
void StateSpecies::remove_all_event_ref(const vector <unsigned int> &ind_list)
{
	for(auto j = 0u; j < ind_list.size(); j++){
		auto &event = individual[ind_list[j]].ev;
	
		for(auto e = 0u; e < event.size(); e++){
			const auto &ev = event[e];
			switch(ev.type){
			case NM_TRANS_EV:
				{
					auto &nmter = nm_trans_ev_ref[ev.m][ev.ti];

					auto nmt_ind = ev.index;
			
					if(nmt_ind+1 < nmter.size()){
						const auto &eref = nmter[nmter.size()-1];
						nmter[nmt_ind] = eref;
					
						individual[eref.i].ev[eref.e_end].index = nmt_ind;
					}
					nmter.pop_back();	
				}
				break;
				
			case M_TRANS_EV:
				{
					auto e = ev.m;
					auto ti = ev.ti;
					auto &it = markov_eqn_vari[e].div[ti].ind_trans;
					auto k = ev.index;
				
					if(k >= it.size()) emsg(" out of range b");
					if(k != it.size()-1){
						auto &it_shift = it[it.size()-1];
						it[k] = it_shift;
						individual[it_shift.i].ev[it_shift.index].index = k;
					}
					it.pop_back();
				}
				break;
				
			default: break;
			}
		}
	}
}


///  Adds event references to individuals
void StateSpecies::add_all_event_ref(const vector <unsigned int> &ind_list, const vector < vector <Event> > &ev_new)
{
	for(auto j = 0u; j < ind_list.size(); j++){
		auto i = ind_list[j];
		auto &event = individual[i].ev;
	
		event = ev_new[j];
		for(auto e = 0u; e < event.size(); e++){
			auto &ev = event[e];
			switch(ev.type){
			case NM_TRANS_EV:
				{
					NMEventRef eref;
					eref.i = i;
					eref.e_end = e;
					auto ti = ev.ti;
					auto m = ev.m;
				
					auto &nmter = nm_trans_ev_ref[m][ti];
					ev.index = nmter.size();
					nmter.push_back(eref);
				}
				break;
				
			case M_TRANS_EV:
				{
					auto e = ev.m;
					auto ti = ev.ti;
					
					auto &mev = markov_eqn_vari[e];
						
					auto &it = mev.div[ti].ind_trans;	
					ev.index = it.size();
						
					DivIndRef dir; dir.i = i; dir.index = e;
					it.push_back(dir);
				}
				break;
				
			default: break;
			}
		}
	}
}


/// Removes event references from individuals
void StateSpecies::remove_all_nmincomp_ref(const vector <unsigned int> &ind_list)
{
	for(auto j = 0u; j < ind_list.size(); j++){
		auto &incomp = individual[ind_list[j]].incomp_ref;
	
		for(auto cl = 0u; cl < sp.ncla; cl++){
			auto &inm = incomp[cl];
			if(inm.on){
				auto &nmter = nm_trans_incomp_ref[inm.n][inm.ti];
				auto index = inm.index;
				if(index+1 < nmter.size()){
					auto i = nmter[nmter.size()-1];
					nmter[index] = i;
					const auto &nmt = sp.nm_trans_incomp[inm.n];
					individual[i].incomp_ref[nmt.cl].index = index;
				}
				nmter.pop_back();
			}
		}
	}
}


///  Adds event references to individuals
void StateSpecies::add_all_nmincomp_ref(const vector <unsigned int> &ind_list, const vector < vector <IncompNMTransRef> > &incomp_ref_new)
{
	for(auto j = 0u; j < ind_list.size(); j++){
		auto i = ind_list[j];
		auto &incomp = individual[i].incomp_ref;
	
		incomp = incomp_ref_new[j];
		
		for(auto cl = 0u; cl < sp.ncla; cl++){
			auto &inm = incomp[cl];
			if(inm.on){
				auto &nmter = nm_trans_incomp_ref[inm.n][inm.ti];
				inm.index = nmter.size();
				nmter.push_back(i);
			}
		}
	}
}


/// Sets e_origin (and m and ti) for an event
void StateSpecies::set_e_origin(Event &ev, unsigned int e, const vector <Event> &event)
{
	auto e_bef = get_event_before(e,event);
	ev.e_origin = e_bef;
	
	const auto &ev_orig = event[e_bef];
	ev.ti = get_ti(ev_orig.tdiv);
	
	ev.m = sp.get_tra_m(ev.tr_gl,ev_orig);
}


/// Sets incomp_ref
void StateSpecies::set_incomp_ref(unsigned int i, const vector < vector <double> > &popnum_t, double &like_ch)
{ 
	if(std::isnan(like_ch)){
		emsg("Incomplete not a num before");
	}
	
	auto &ind = individual[i];
	auto &inmtr = ind.incomp_ref;
	
	const auto &ev = ind.ev;
	
	auto si = ev.size();
	
	if(si == 0){
		for(auto cl = 0u; cl < sp.ncla; cl++){
			auto &inm = inmtr[cl];
			if(inm.on){ 
				like_ch -= inm.Li;
				incomp_turn_off(inm);
			}
		}
		return;
	}
	
	auto el = ev.size()-1;
	
	double t_end = details.T;
	auto dt = details.dt;
	
	// Finds a list of classification which need incomplete likelihoods
	auto c_last = ev[el].c_after;
	
	if(c_last == UNSET){
		t_end = ev[el].tdiv;
		el--; 
		c_last = ev[el].c_after;
	}
	if(c_last == UNSET) emsg("c_last unset");
	
	const auto &cl_on = sp.cgl_incomp_nmtrans_cl[c_last];
	
	for(auto cl = 0u; cl < sp.ncla; cl++){
		auto &inm = inmtr[cl];
		if(cl_on[cl]){
			auto e = el; 
			while(e > 0 && ev[e].type != ENTER_EV &&
				!((ev[e].type == M_TRANS_EV || ev[e].type == NM_TRANS_EV) && ev[e].cl == cl)){
				e--;
			}

			auto &eve = ev[e];
			
			auto c = eve.c_after;
			
			auto n = sp.comp_gl[c].nmtransincomp_ref[cl];
			if(n == UNSET) emsg("Should not be unset");
		
			auto t = eve.tdiv;
			
			auto ti = get_ti(t);
	
			auto dtdiv = t_end - t; 
			if(dtdiv <= 0){
				emsg("zero time2");
			}
			
			if(!(inm.on && inm.n == n && inm.ti == ti && inm.dtdiv == dtdiv && inm.e_begin == e && inm.tdiv_end == t_end)){
				if(inm.on){ 
					like_ch -= inm.Li;
					incomp_turn_off(inm);
				}
				incomp_turn_on(i,n,ti,dtdiv,e,t_end,inm); 
				
				const auto &nmt = sp.nm_trans_incomp[n];
				auto val = get_nm_incomp_val(nmt,ti,ind,popnum_t);
				inm.Li = nm_trans_incomp_full_like(nmt.nmtrans_ref,dtdiv,dt,val.ref_val,val.bp_val);
				like_ch += inm.Li;
			}
		}
		else{	
			if(inm.on){ 
				like_ch -= inm.Li;
				incomp_turn_off(inm);
			}
		}
	}		
	if(std::isnan(like_ch)){
		emsg("Incomplete not a num");
	}
}


/// Adjusts incomp_ref for joint proposals
void StateSpecies::adjust_incomp_ref(vector <IncompNMTransRef> &inmtr, const vector <Event> &ev)
{
	auto el = ev.size()-1;
	
	double t_end = T;
	
	// Finds a list of classification which need incomplete likelihoods
	auto c_last = ev[el].c_after;
	if(c_last == UNSET){
		t_end = ev[el].tdiv;
		el--; 
		c_last = ev[el].c_after;
	}
	//const auto &cl_on = sp.cgl_incomp_nmtrans_cl[c_last];
	
	for(auto cl = 0u; cl < sp.ncla; cl++){
		auto &inm = inmtr[cl];
		if(inm.on){
			auto e = el; 
			while(e > 0 && ev[e].type != ENTER_EV &&
				!((ev[e].type == M_TRANS_EV || ev[e].type == NM_TRANS_EV) && ev[e].cl == cl)){
				e--;
			}

			auto &eve = ev[e];
			
			auto c = eve.c_after;
			
			auto n = sp.comp_gl[c].nmtransincomp_ref[cl];
			if(n == UNSET) emsg("Should not be unset");
		
			auto t = eve.tdiv;
			
			auto dtdiv = t_end - t; if(dtdiv <= 0) emsg("zero time3");
		
			inm.dtdiv = dtdiv;
			inm.e_begin = e;
			inm.n = n;
			inm.ti = get_ti(t);
			inm.tdiv_end = t_end;
		}
	}
}


/// Sets up non-Markovian transitions
void StateSpecies::setup_nm_trans(const vector < vector <double> > &popnum_t)
{
	if(sp.nm_flag == false) return;
	
	// Clears existing arrays
	for(auto m = 0u; m < sp.nm_trans.size(); m++){
		for(auto ti = 0u; ti < T; ti++) nm_trans_ev_ref[m][ti].clear();
	}
	
	nm_trans_incomp_ref.resize(sp.nm_trans_incomp.size());
	for(auto n = 0u; n < sp.nm_trans_incomp.size(); n++){
		for(auto ti = 0u; ti < T; ti++) nm_trans_incomp_ref[n][ti].clear();
	}
	
	for(auto &ind : individual){
		for(auto &ir : ind.incomp_ref) ir.on = false;
	}
		
	Like temp;
	for(auto i = 0u; i < individual.size(); i++){
		auto &ind = individual[i];
		auto &event = ind.ev;
		for(auto e = 0u; e < event.size(); e++){
			auto &ev = event[e];
			if(ev.type == NM_TRANS_EV){
				set_e_origin(ev,e,event);
				add_event_ref(i,e,popnum_t,temp);
			}
		}			
		
		{
			auto temp = 0.0;
			set_incomp_ref(i,popnum_t,temp);
		}
	}
}


/// Turns on a reference in nm_trans_incomp_ref 
void StateSpecies::incomp_turn_on(unsigned int i, unsigned int n, unsigned int ti, double dtdiv, unsigned int e_begin, double t_end, IncompNMTransRef &inm)
{   
	inm.n = n;
	inm.ti = ti;
	inm.dtdiv = dtdiv;
	inm.e_begin = e_begin;
	inm.tdiv_end = t_end;		
			
	auto &nmtr = nm_trans_incomp_ref[n][ti];
	inm.index = nmtr.size();
	nmtr.push_back(i);
	inm.on = true;
}


/// Turns off a reference in nm_trans_incomp_ref 
void StateSpecies::incomp_turn_off(IncompNMTransRef &inm)
{
	auto &nmter = nm_trans_incomp_ref[inm.n][inm.ti];
	auto index = inm.index;
	if(index+1 < nmter.size()){
		auto i = nmter[nmter.size()-1];
		nmter[index] = i;
		const auto &nmt = sp.nm_trans_incomp[inm.n];
		individual[i].incomp_ref[nmt.cl].index = index;
	}
	nmter.pop_back();
	inm.on = false;
}


/// Updates init_cond_val based on an individuals enter state changing
void StateSpecies::likelihood_init_cond_change(unsigned int c_enter_old, unsigned int c_enter_new, Like &like_ch)
{	
	const auto &ic = sp.init_cond;
	
	auto &icv = init_cond_val;
	
	auto foc_cl = ic.focal_cl;
	
	if(foc_cl == UNSET){
		int dN = 0;
		auto dLi_init_cond = 0.0;
		
		if(c_enter_old != UNSET){ 	
			icv.cnum[c_enter_old]--; 	
			if(sp.type != INDIVIDUAL){
				dLi_init_cond += log_int(icv.cnum[c_enter_old]+1);	
			}
			auto fr = icv.frac[c_enter_old];
			if(fr < TINY) dLi_init_cond += LARGISH;
			else dLi_init_cond -= log(fr);
			dN--;
		}
		
		if(c_enter_new != UNSET){ 
			icv.cnum[c_enter_new]++; 
			if(sp.type != INDIVIDUAL){
				dLi_init_cond -= log_int(icv.cnum[c_enter_new]);	
			}
			auto fr = icv.frac[c_enter_new];
			if(fr < TINY) dLi_init_cond -= LARGISH;
			else dLi_init_cond += log(fr);
			dN++;
		}
		
		if(dN != 0){
			if(sp.type == INDIVIDUAL){
				if(dN == 1) dLi_init_cond += log_int(icv.N_total_unobs+1);
				else dLi_init_cond -= log_int(icv.N_total_unobs);
				icv.N_total_unobs += dN;			
			}
			else{
				if(dN == 1) dLi_init_cond += log_int(icv.N_total+1);
				else dLi_init_cond -= log_int(icv.N_total);
			}
			
			auto dprior = -prior_probability(icv.N_total,ic.pop_prior,param_val,eqn);
			
			icv.N_total += dN;
			dprior += prior_probability(icv.N_total,ic.pop_prior,param_val,eqn);
		
			like_ch.init_cond_prior += dprior;
		}
		
		back.push_back(Back(IND_COND_VAL,c_enter_old,c_enter_new,dLi_init_cond));
	
		Li_init_cond += dLi_init_cond;
		like_ch.init_cond += dLi_init_cond;
	}
	else{
		auto dLi_init_cond = 0.0;
		int dN = 0;
		
		if(c_enter_old != UNSET){ 
			const auto &crr = ic.comp_reduce_ref[c_enter_old];
			auto c = crr.c;
			auto cred = crr.cred;
			
			like_ch.init_cond_prior -= prior_probability(icv.N_focal[c],ic.comp_prior[c],param_val,eqn);
			
			icv.N_focal[c]--;
			icv.cnum_reduce[c][cred]--;
			icv.cnum[c_enter_old]--;
			
			if(sp.type == INDIVIDUAL){
				dLi_init_cond -= log_int(icv.N_focal_unobs[c]);
				icv.N_focal_unobs[c]--;
			}
			else{
				dLi_init_cond -= log_int(icv.N_focal[c]+1) - log_int(icv.cnum_reduce[c][cred]+1);	
			}
			like_ch.init_cond_prior += prior_probability(icv.N_focal[c],ic.comp_prior[c],param_val,eqn);
			
			auto fr = icv.frac_comb[cred];
			if(fr < TINY) dLi_init_cond += LARGISH;
			else dLi_init_cond -= log(fr);
			dN--;
		}
		
		if(c_enter_new != UNSET){ 
			const auto &crr = ic.comp_reduce_ref[c_enter_new];
			auto c = crr.c;
			auto cred = crr.cred;
			
			like_ch.init_cond_prior -= prior_probability(icv.N_focal[c],ic.comp_prior[c],param_val,eqn);
			
			icv.N_focal[c]++;
			icv.cnum_reduce[c][cred]++;
			icv.cnum[c_enter_new]++;
				
			if(sp.type == INDIVIDUAL){
				icv.N_focal_unobs[c]++;
				dLi_init_cond += log_int(icv.N_focal_unobs[c]);
			}
			else{
				dLi_init_cond += log_int(icv.N_focal[c]) - log_int(icv.cnum_reduce[c][cred]);	
			}
			
			like_ch.init_cond_prior += prior_probability(icv.N_focal[c],ic.comp_prior[c],param_val,eqn);
			
			auto fr = icv.frac_comb[cred];
			if(fr < TINY) dLi_init_cond -= LARGISH;
			else dLi_init_cond += log(fr);
			dN++;
		}
	
		back.push_back(Back(IND_COND_VAL,c_enter_old,c_enter_new,dLi_init_cond));
	
		Li_init_cond += dLi_init_cond;
		like_ch.init_cond += dLi_init_cond;
	}	
}


/// Used for going back to the previous state
void StateSpecies::init_cond_change_back(unsigned int c_enter_old, unsigned int c_enter_new)
{		
	const auto &ic = sp.init_cond;
	auto &icv = init_cond_val;
	auto foc_cl = ic.focal_cl;
	
	if(foc_cl == UNSET){
		int dN = 0;	
		if(c_enter_old != UNSET){ 
			icv.cnum[c_enter_old]--; 	
			dN--;
		}
		
		if(c_enter_new != UNSET){ 
			icv.cnum[c_enter_new]++; 
			dN++;
		}
		
		if(dN != 0){
			if(sp.type == INDIVIDUAL){
				icv.N_total_unobs += dN;
			}
			icv.N_total += dN;
		}
	}
	else{
		int dN = 0;
		
		if(c_enter_old != UNSET){ 
			const auto &crr = ic.comp_reduce_ref[c_enter_old];
			auto c = crr.c;
			auto cred = crr.cred;
			
			icv.N_focal[c]--;
			icv.cnum_reduce[c][cred]--;
			icv.cnum[c_enter_old]--;
			
			if(sp.type == INDIVIDUAL){
				icv.N_focal_unobs[c]--;
			}
	
			dN--;
		}
		
		if(c_enter_new != UNSET){ 
			const auto &crr = ic.comp_reduce_ref[c_enter_new];
			auto c = crr.c;
			auto cred = crr.cred;
				
			icv.N_focal[c]++;
			icv.cnum_reduce[c][cred]++;
			icv.cnum[c_enter_new]++;
			if(sp.type == INDIVIDUAL) icv.N_focal_unobs[c]++;
			
			dN++;
		}
	}
}


/// Recalulates the value of Markov equation (because population change)
void StateSpecies::recalc_markov_value(unsigned int ee, unsigned int ti, unsigned int ti_next, const vector < vector <double> > &popnum_t, const vector <PopChange> &pop_change, double &like_ch)
{
	if(type != INDIVIDUAL) emsg("must be ind");
	
	const auto dt = details.dt;
	
	const auto &me = sp.markov_eqn[ee];
	auto &me_vari = markov_eqn_vari[ee];

	const auto &eq = eqn[me.eqn_ref];
	const auto &lin = eq.linearise;
	
	auto &Li_m = Li_markov[ee];
	
	if(me_vari.time_vari == false) emsg("must have time variation");
	
	// Different type of update:
	// RECALC This recalculates each time
	// USE_POP_DIF This uses just the change in populations (when equation can be linearised)
	// USE_POP_DIF_TIME This uses just the change in populations (when gradient time dependent)
	
	auto type = RECALC;
	if(linearise_speedup3 && lin.on){
		type = USE_POP_DIF;	
		if(lin.pop_grad_calc_time_dep) type = USE_POP_DIF_TIME;
	}
	
	vector <unsigned int> pr_store;
	auto diff = 0.0;
	switch(type){
	case RECALC:              // If recalculating
		break;
		
	case USE_POP_DIF:         // If population gradient not time dependant
		for(const auto &po_ch : pop_change){
			auto pr = lin.get_pop_ref(po_ch.po);
			if(pr != UNSET){
				auto va = dt*eq.calculate_calculation_notime(lin.pop_grad_calc[pr],param_val);
				diff += va*po_ch.num;
			}
		}
		break;
	
	case USE_POP_DIF_TIME:         // If population gradient time dependant
		// Makes a list of populations
		for(auto k = 0u; k < pop_change.size(); k++){
			auto pr = lin.get_pop_ref(pop_change[k].po);
			pr_store.push_back(pr);
		}
		break;
	}

	vector <double> dLi_store, value_store;
	double value;
	for(auto tii = ti; tii < ti_next; tii++){
		auto &div = me_vari.div[tii];
		
		auto val_old = div.value;
	
		switch(type){
		case RECALC:
			value = dt*eq.calculate(tii,popnum_t[tii],param_val,spline_val);	
			break;
			
		case USE_POP_DIF:
			value = val_old+diff;
			break;
			
		case USE_POP_DIF_TIME:
			{
				value = val_old;
				auto sum = 0.0;
				for(auto k = 0u; k < pop_change.size(); k++){
					auto pr = pr_store[k];
					if(pr != UNSET){
						const auto &po_ch = pop_change[k];
						auto va = eq.calculate_calculation(lin.pop_grad_calc[pr],tii,param_val,spline_val);
						sum += va*po_ch.num;
					}
				}
				value += dt*sum;
			
				if(false){
					auto val_ch = eq.calculate(tii,popnum_t[tii],param_val,spline_val);
					if(dif(value,val_ch,dif_thresh)){
						emsg("problem with value");
					}
				}
			}
			break;
		}
	
		auto dLi = -(value-val_old)*div.indfac_int;
		
		auto n = div.ind_trans.size();	
		if(n > 0 && !me.infection_trans){
			dLi += n*log((value+LOG_THRESH)/(val_old+LOG_THRESH));
		}
		dLi_store.push_back(dLi); value_store.push_back(div.value);
		
		Li_m[tii] += dLi;
		like_ch += dLi;
		div.value = value;
	}
	
	back.push_back(Back(LI_MARKOV,ee,ti,dLi_store));
	back.push_back(Back(VALUE_MARKOV,ti_next,value_store));
}

	
/// Restores the state to the orignal
void StateSpecies::restore_back()
{
	// Updates references
	for(auto &ba_ind : back_ind){
		// Removes references for new sequence
		auto &ind = individual[ba_ind.i];
		auto &ev_new = ind.ev;
		
		for(auto ee : ba_ind.add_new){
			const auto &ev = ev_new[ee];
			switch(ev.type){
			case M_TRANS_EV:
				{
					auto &it = markov_eqn_vari[ev.m].div[ev.ti].ind_trans;
					auto k = ev.index;
					if(k != it.size()-1){
						auto &it_shift = it[it.size()-1];
						it[k] = it_shift;
						individual[it_shift.i].ev[it_shift.index].index = k;
					}
					it.pop_back();
				}
				break;
				
			case NM_TRANS_EV:
				{
					auto &nmter = nm_trans_ev_ref[ev.m][ev.ti];
					auto nmt_ind = ev.index;
	
					if(nmt_ind+1 < nmter.size()){
						const auto &eref = nmter[nmter.size()-1];
						nmter[nmt_ind] = eref;
	
						individual[eref.i].ev[eref.e_end].index = nmt_ind;
					}
					nmter.pop_back();
				}
				break;
			
			default: break;
			}
		}
		
		auto &ev_old = ba_ind.ev;
		for(auto j = 0u; j < ba_ind.identical_old.size(); j++){
			auto e = ba_ind.identical_old[j];
			auto &eve_old = ev_old[e];
			auto ee = ba_ind.identical_new[j];
			eve_old.index = ev_new[ee].index; 
			
			if(e != ee){
				switch(eve_old.type){
				case M_TRANS_EV:
					markov_eqn_vari[eve_old.m].div[eve_old.ti].ind_trans[eve_old.index].index = e;
					break;
					
				case NM_TRANS_EV:
					nm_trans_ev_ref[eve_old.m][eve_old.ti][eve_old.index].e_end = e;
					break;
				
				default: emsg("SHoudl nbh"); break;
				}
			}
		}
		
		ev_new = ev_old;
			
		for(auto ee : ba_ind.remove_old){
			auto &ev = ev_new[ee];
		
			switch(ev.type){
			case M_TRANS_EV:
				{
					auto &it = markov_eqn_vari[ev.m].div[ev.ti].ind_trans;	
					ev.index = it.size();
			
					DivIndRef dir; dir.i = ba_ind.i; dir.index = ee;
					it.push_back(dir);
				}
				break;
				
			case NM_TRANS_EV:
				{
					auto &nmter = nm_trans_ev_ref[ev.m][ev.ti];
					ev.index = nmter.size();
					
					NMEventRef eref; eref.i = ba_ind.i; eref.e_end = ee;
					nmter.push_back(eref);	
				}
				break;
			
			default: break;
			}
		}
		
		auto &inm_list = ind.incomp_ref;
		const auto &inm_st_list = ba_ind.incomp_ref;
		for(auto cl = 0u; cl < sp.ncla; cl++){
			auto &inm = inm_list[cl];
			const auto &inm_st = inm_st_list[cl];
			
			if(inm_st.on){
				if(!(inm.on && inm.n == inm_st.n && inm.ti == inm_st.ti && inm.dtdiv == inm_st.dtdiv && inm.e_begin == inm_st.e_begin && inm.tdiv_end == inm_st.tdiv_end)){
					if(inm.on) incomp_turn_off(inm);
					incomp_turn_on(ba_ind.i,inm_st.n,inm_st.ti,inm_st.dtdiv,inm_st.e_begin,inm_st.tdiv_end,inm); 
					inm.Li = inm_st.Li;
				}
			}
			else{
				if(inm.on) incomp_turn_off(inm);
			}
		}
	}
	
	if(back.size() > 0){
		for(int b = back.size()-1; b >= 0; b--){
			const auto &ba = back[b];
			switch(ba.type){
			case REMOVE_LI_MARKOV:
				Li_markov[ba.i][ba.index] += ba.value;
				break;
				
			case ADD_LI_MARKOV:
				Li_markov[ba.i][ba.index] -= ba.value;
				break;
				
			case POP_DATA_NUM:
				pop_data_num[ba.i] = ba.value;
				break;
				
			case LI_OBS_POP:
				Li_obs_pop[ba.i] = ba.value;
				break;
			
			case POP_DATA_CGL:
				pop_data_cgl[ba.i][ba.index] -= ba.value;
				break;
				
			case POP_TRANS_DATA_NUM:
				{
					const auto &vec = ba.vec;
					auto sign = ba.value;
					auto tr_gl = ba.index;
					auto k = 0u;
					for(auto ref : sp.pop_trans_ref[ba.i][tr_gl]){
						pop_trans_data_tgl[ref][tr_gl] -= sign;	
						pop_trans_data_num[ref] -= vec[k]; k++;
					}
				}
				break;
				
			case LI_OBS_POP_TRANS:
				{
					const auto &vec = ba.vec;
					auto k = 0u;
					for(auto ref : sp.pop_trans_ref[ba.i][ba.index]){
						Li_obs_pop_trans[ref] -= vec[k]; k++;
					}
				}
				break;
				
			case LI_OBS_IND:
				Li_obs_ind[ba.i] -= ba.value;
				break;
				
			case TI_INDFAC:
				{
					auto e = ba.i;
					auto &mev = markov_eqn_vari[e];
					auto &Li_m = Li_markov[e];
					const auto &ti = ba.vec;
					const auto &dLi = back[b+1].vec;
					const auto &dif = back[b+2].vec;
					for(auto k = 0u; k < ti.size(); k++){
						unsigned int tii = ti[k]; 
						Li_m[tii] -= dLi[k];
						mev.div[tii].indfac_int -= dif[k];
					}
				}
				break;
				
			case DLI_INDFAC: case DIF_INDFAC: break;
			
			case LI_MARKOV: 
				{
					auto e = ba.i, ti = ba.index, ti_next = back[b+1].i;
					const auto &Li = ba.vec;
					const auto &value = back[b+1].vec;
					auto &Li_m = Li_markov[e];
					auto &me_vari = markov_eqn_vari[e];

					auto k = 0u;
					for(auto tii = ti; tii < ti_next; tii++){
						Li_m[tii] -= Li[k];
						me_vari.div[tii].value = value[k];
						k++;
					}
				}
				break;
				
			case OBS_TRANS_EQN_NUM:
				{
					auto tr_gl = ba.i;
					auto ti = ba.index;
					//auto sign = -ba.value;
					auto k = 0u;
					for(auto m : sp.obs_trans_eqn_ref[tr_gl][ti]){
						auto sign = -ba.vec[k]; k++;
						obs_trans_eqn_num[m][ti] += sign;
					}
				}
				break;
			
			case VALUE_MARKOV: break;
		
			case IND_COND_VAL:
				Li_init_cond -= ba.value;
				init_cond_change_back(ba.index,ba.i); 
				break;
			}
		}
	}
}

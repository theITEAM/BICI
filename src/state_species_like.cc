// Calculates different types of likelihood

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
 
using namespace std;

#include "state_species.hh"
#include "utils.hh"
#include "matrix.hh"


/// Markov likelihood for a specific equation
vector <double> StateSpecies::likelihood_markov(unsigned int e, const vector <unsigned int> &list, double &like_ch)
{
	auto &me = sp.markov_eqn[e];
	auto &me_vari = markov_eqn_vari[e];
	auto &Li_mark = Li_markov[e];
	
	vector <double> store;
	
	double va;
	
	for(auto ti : list){
		auto &div = me_vari.div[ti];
		
		auto val = div.value;
		if(val < 0){
			if(val > -SMALL){
				div.value = 0;
				val = 0;
			}
			else{
				emsg("Rate has become negative");
			}
		}
		
		auto Li = -val*div.indfac_int;
	
		auto si = div.ind_trans.size();
			
		if(me.infection_trans){ // If an infection transition then works out specific probability of specific individual
			const auto &eq = eqn[me.eqn_ref];
			const auto &lin = eq.linearise;
			if(!lin.on) emsg("Linearisation should be on");
			for(const auto &it : div.ind_trans){
				const auto &ind = individual[it.i];
				const auto &iif = ind.ev[it.index].ind_inf_from;
				auto p = iif.p;
				switch(p){
				case OUTSIDE_INF:
					va = eq.calculate_calculation(lin.no_pop_calc,ti,param_val,spline_val);
					break;
					
				case ENTER_INF: 
					va = UNSET;
					emsg("SHould not be ENTER_INFb");
					break;
					
				default:
					va = eq.calculate_calculation(lin.pop_grad_calc[iif.pref],ti,param_val,spline_val);
					break;
				}
				
				if(me.ind_variation) va *= get_indfac(ind,me);
				Li += log(va*iif.w);
			}
		}
		else{
			switch(si){
			case 0: break;
			
			case 1:
				if(me.ind_variation){
					Li += log(get_indfac(individual[div.ind_trans[0].i],me)*(val+LOG_THRESH));	
				}
				else Li += log(val+LOG_THRESH);	
				break;
				
			default:
				if(me.ind_variation){
					for(const auto &it : div.ind_trans) Li += log(get_indfac(individual[it.i],me));
				}
				Li += si*log(val+LOG_THRESH);
				break;
			}
		}
		
		like_ch += Li-Li_mark[ti];
		store.push_back(Li_mark[ti]);
		
		Li_mark[ti] = Li;
	}
	
	return store;
}


/// Calculate the value of a markov equation
vector <double> StateSpecies::likelihood_markov_value(unsigned int e, const vector <unsigned int> &list, const vector < vector <double> > &popnum_t)
{
	auto &me = sp.markov_eqn[e];
	auto &me_vari = markov_eqn_vari[e];

	vector <double> store;

	if(me.time_vari == false){
		auto value = eqn[me.eqn_ref].calculate_param_only(param_val);
		store.push_back(me_vari.div[0].value);
		me_vari.div[0].value = value;
	}		
	else{
		const auto &eq = eqn[me.eqn_ref];
		auto &div = me_vari.div;
		
		for(auto ti : list){	
			auto value = eq.calculate(ti,popnum_t[ti],param_val,spline_val);
			store.push_back(div[ti].value);
			div[ti].value = value;
		}
	}
	
	return store;
}


/// Calculate the value of a markov equation 
// This is the case when only the non-population part of the equation changes
vector <double> StateSpecies::likelihood_markov_value_fast(const vector <unsigned int> &me_list, const vector <unsigned int> &list, const vector < vector <double> > &popnum_t, const vector <double> &param_store, const vector <SplineValue> &spline_val)
{
	auto e_first = me_list[0];
	auto &me = sp.markov_eqn[e_first];
	
	vector <double> store;

	if(me.time_vari == false) emsg("time_vary should be true");
	
	const auto &eq = eqn[me.eqn_ref];
	const auto &calc = eq.linearise.no_pop_calc;
	
	if(eq.linearise.no_pop_calc_time_dep){ // Time dependent	
		for(auto ti : list){	
			auto d = eq.calculate_calculation(calc,ti,param_val,spline_val) - 
								 eq.calculate_calculation_spline_store(calc,ti,param_store,spline_val);
								
			for(auto e : me_list){		
				auto &div = markov_eqn_vari[e].div;			
				auto val = div[ti].value;
			
				if(slow_check){
					auto value = eq.calculate(ti,popnum_t[ti],param_val,spline_val);
					if(dif(value,val+d,dif_thresh)) emsg("markov fast problem");
				}
			
				store.push_back(val);
				div[ti].value = val+d;
			}
		}
	}
	else{                                            // Not time dependent
		auto d = eq.calculate_calculation(calc,UNSET,param_val,spline_val) - 
								 eq.calculate_calculation_spline_store(calc,UNSET,param_store,spline_val);
			
		for(auto ti : list){							
			for(auto e : me_list){		
				auto &div = markov_eqn_vari[e].div;			
				auto val = div[ti].value;
			
				if(slow_check){
					auto value = eq.calculate(ti,popnum_t[ti],param_val,spline_val);
					if(dif(value,val+d,dif_thresh)) emsg("markov fast problem");
				}
			
				store.push_back(val);
				div[ti].value = val+d;
			}
		}
	}
	
	return store;
}


/// Calculate the value of a series of markov equation
vector <double> StateSpecies::likelihood_markov_value_linear(const vector <unsigned int> &list, const LinearProp &linear_prop, const vector < vector <double> > &popnum_t)
{
	timer[LINEAR_TIMER] -= clock();
	
	if(list.size() == 0) emsg("list should not be zero");
	
	vector <double> store;
	
	auto N = linear_prop.me.size();
	vector <double> val(N);
	
	vector <double> no_pop;
	no_pop.resize(N);
	
	vector < vector <double> > pop_grad;
	pop_grad.resize(N);
	
	auto ti_last = UNSET;
	
	for(auto k = 0u; k < list.size(); k++){
		auto ti = list[k];
		
		if(k == 0){                            // Calculates the first value (along with gradients)
			timer[LINEARZERO_TIMER] -= clock();
	
			for(auto i = 0u; i < N; i++){
				auto e = linear_prop.me[i]; 			
				auto &me = sp.markov_eqn[e];

				const auto &eq = eqn[me.eqn_ref];
	
				val[i] = eq.calculate(ti,popnum_t[ti],param_val,spline_val);
				pop_grad[i] = eq.calculate_popnum_gradient(param_val);
				
				if(eq.linearise.no_pop_calc_time_dep){
					no_pop[i] = eq.calculate_no_pop(ti,param_val,spline_val);
				}
			}
			timer[LINEARZERO_TIMER] += clock();
		}
		else{                                    // Subsequent values uses differences in population
			const auto &popnum = popnum_t[ti];
			const auto &popnum_last = popnum_t[ti_last];
			
			for(const auto &pa : linear_prop.pop_affect){
				auto po = pa.po;
				if(popnum[po] != popnum_last[po]){
					auto dif = rectify(popnum[po]) - rectify(popnum_last[po]);
					for(auto &gr : pa.pop_grad_ref){
						val[gr.ref] += dif*pop_grad[gr.ref][gr.index];
					}
				}
			}
		}
		
		for(auto i = 0u; i < N; i++){
			auto e = linear_prop.me[i]; 
			auto &me = sp.markov_eqn[e];
			const auto &eq = eqn[me.eqn_ref];
			
			// Accounts for time variation in the non-population part
			if(eq.linearise.no_pop_calc_time_dep){    
				val[i] -= no_pop[i];
				no_pop[i] = eq.calculate_no_pop(ti,param_val,spline_val);
				val[i] += no_pop[i];
			}
			
			auto &me_vari = markov_eqn_vari[e];
			store.push_back(me_vari.div[ti].value);
			me_vari.div[ti].value = val[i];
			
			if(slow_check){
				if(dif(val[i],eq.calculate(ti,popnum_t[ti],param_val,spline_val),dif_thresh)){
					emsg("Not agree"); 
				}
			}
		}
		
		ti_last = ti;
	}
	
	timer[LINEAR_TIMER] += clock();
	
	return store;
}


/// Restores values 
void StateSpecies::likelihood_markov_value_linear_restore(const vector <unsigned int> &list, const LinearProp &linear_prop, const vector<double> &store)
{
	auto i = 0u;
	for(auto k = 0u; k < list.size(); k++){
		auto ti = list[k];
		
		for(auto e : linear_prop.me){
			auto &me_vari = markov_eqn_vari[e];
			me_vari.div[ti].value = store[i]; i++;
		}
	}
	if(i != store.size()) emsg("Store problem");
}


/// Sets indfac_int in markov_eqn
// This can potentially be sped up by only updating a certain selection of markov_eqns
vector <double> StateSpecies::likelihood_indfac_int()
{	
	vector <double> store;

	auto t_start = timepoint[0];
	auto t_end = timepoint[T];
	auto dt = details.dt;
	
	for(auto e = 0u; e < N; e++){
		auto &div = markov_eqn_vari[e].div;
		const auto &me = sp.markov_eqn[e];
		auto val = 0.0; 
		if(me.source == true) val = me.source_tr_gl.size()*dt*T/div.size();
		for(auto ti = 0u; ti < div.size(); ti++){
			store.push_back(div[ti].indfac_int);
			div[ti].indfac_int = val;
		}
	}
	
	for(auto i = 0u; i < individual.size(); i++){
		const auto &ind = individual[i];
		
		auto t = t_start;
		auto ti = 0u;
		
		auto c = UNSET;
		for(auto k = 0u; k <= ind.ev.size(); k++){
			auto t_ev = t_end; if(k < ind.ev.size()) t_ev = ind.ev[k].t;
	
			auto ti_new = (unsigned int)(ALMOST_ONE*(t_ev-t_start)/dt);
								
			if(c != UNSET){
				for(auto e : sp.comp_gl[c].me_ref){
					auto &mev = markov_eqn_vari[e];
					const auto &me = sp.markov_eqn[e];
					auto indfac = get_indfac(ind,me);
				
					if(mev.time_vari == false){
						mev.div[0].indfac_int += (t_ev-t)*indfac;
					}
					else{
						if(ti_new == ti){
							mev.div[ti].indfac_int += (t_ev-t)*indfac;
						}
						else{
							mev.div[ti].indfac_int += (timepoint[ti+1]-t)*indfac;
							for(auto j = ti+1; j < ti_new; j++){
								mev.div[j].indfac_int += dt*indfac;
							}
							mev.div[ti_new].indfac_int += (t_ev-timepoint[ti_new])*indfac;
						}
					}
				}
			}
			
			if(k < ind.ev.size()) c = ind.ev[k].c_after;
			
			t = t_ev;
			ti = ti_new;
		}
	}
	
	return store;
}


/// Looks at updating indfac_int based on a change to an ie
vector <double> StateSpecies::likelihood_ie_change(unsigned int i, unsigned int ie, double factor, Like &like_ch)
{
	vector <double> store;

	auto t_start = timepoint[0];
	auto t_end = timepoint[T];
	auto dt = details.dt;
	
	auto t = t_start;
	auto ti = 0u;
	
	const auto &ind = individual[i];
	const auto &ev = ind.ev;
	
	auto c = UNSET;
	auto dLi_sum = 0.0;
	
	double va, dLi;
			
	for(auto k = 0u; k <= ev.size(); k++){
		auto t_ev = t_end; if(k < ind.ev.size()) t_ev = ev[k].t;

		auto ti_new = (unsigned int)(ALMOST_ONE*(t_ev-t_start)/dt);
							
		if(c != UNSET){
			for(auto e : sp.comp_gl[c].me_ref){
				const auto &me = sp.markov_eqn[e];
				
				for(auto ie_me : me.ind_eff_mult){
					if(ie_me == ie){				
						auto &mev = markov_eqn_vari[e];
						auto &div = mev.div;
					
						auto indfac = (factor-1)*get_indfac(ind,me);
					
						auto &Li_mark = Li_markov[e];
					
						if(mev.time_vari == false){
							va = (t_ev-t)*indfac;
							div[0].indfac_int += va;
							dLi = -va*div[0].value;
							Li_mark[0] += dLi;
							dLi_sum += dLi;
							store.push_back(va);
							store.push_back(dLi);
						}
						else{
							if(ti_new == ti){
								va = (t_ev-t)*indfac;
								div[ti].indfac_int += va;
								dLi = -va*div[ti].value;
								Li_mark[ti] += dLi;
								dLi_sum += dLi;
								store.push_back(va);
								store.push_back(dLi);
							}
							else{
								va = (timepoint[ti+1]-t)*indfac;
								div[ti].indfac_int += va;
								dLi = -va*div[ti].value;
								Li_mark[ti] += dLi; 
								dLi_sum += dLi;
								store.push_back(va);
								store.push_back(dLi);
								
								for(auto j = ti+1; j < ti_new; j++){
									va = dt*indfac;
									div[j].indfac_int += va;
									dLi = -va*div[j].value;
									Li_mark[j] += dLi;
									dLi_sum += dLi;
									store.push_back(va);
									store.push_back(dLi);
								}
								
								va = (t_ev-timepoint[ti_new])*indfac;
								div[ti_new].indfac_int += va;
								dLi = -va*div[ti_new].value;
								Li_mark[ti_new] += dLi;
								dLi_sum += dLi;
								store.push_back(va);
								store.push_back(dLi);
							}
						}
					}
				}
			}
		}
		
		if(k < ev.size()){
			const auto &e = ev[k];
			
			if(e.type == M_TRANS_EV){
				const auto &me = sp.markov_eqn[e.m];
				
				for(auto ie_me : me.ind_eff_mult){
					if(ie_me == ie){
						dLi = log(factor);
						Li_markov[e.m][e.ti] += dLi;
						like_ch.markov += dLi;
					}
				}
			}	
			c = e.c_after;
		}
		
		t = t_ev;
		ti = ti_new;
	}
	
	like_ch.markov += dLi_sum;
	
	return store;
}


/// Looks at updating indfac_int based on a change to an ie
void StateSpecies::likelihood_ie_change_restore(unsigned int i, unsigned int ie, double factor, const vector <double> &store)
{
	auto t_start = timepoint[0];
	auto t_end = timepoint[T];
	auto dt = details.dt;
	
	auto ti = 0u;
	
	const auto &ind = individual[i];
	const auto &ev = ind.ev;
	
	auto c = UNSET;
	
	auto j = 0u;
	for(auto k = 0u; k <= ev.size(); k++){
		auto t_ev = t_end; if(k < ind.ev.size()) t_ev = ev[k].t;
		auto ti_new = (unsigned int)(ALMOST_ONE*(t_ev-t_start)/dt);
						
		if(c != UNSET){
			for(auto e : sp.comp_gl[c].me_ref){
				const auto &me = sp.markov_eqn[e];
				
				for(auto ie_me : me.ind_eff_mult){
					if(ie_me == ie){
						auto &mev = markov_eqn_vari[e];
						auto &div = mev.div;
					
						auto &Li_mark = Li_markov[e];
					
						if(mev.time_vari == false){
							div[0].indfac_int -= store[j]; j++;
							Li_mark[0] -= store[j]; j++;
						}
						else{
							if(ti_new == ti){
								div[ti].indfac_int -= store[j]; j++;
								Li_mark[ti] -= store[j]; j++;
							}
							else{
								div[ti].indfac_int -= store[j]; j++;
								Li_mark[ti] -= store[j]; j++;	
								for(auto k = ti+1; k < ti_new; k++){
									div[k].indfac_int -= store[j]; j++;
									Li_mark[k] -= store[j]; j++;	
								}
								div[ti_new].indfac_int -= store[j]; j++;
								Li_mark[ti_new] -= store[j]; j++;	
							}
						}
					}
				}
			}
		}
		
		if(k < ev.size()){
			const auto &e = ev[k];
			
			if(e.type == M_TRANS_EV){
				const auto &me = sp.markov_eqn[e.m];
				for(auto ie_me : me.ind_eff_mult){
					if(ie_me == ie){	
						Li_markov[e.m][e.ti] -= log(factor);
					}
				}
			}	
			c = e.c_after;
		}
		
		ti = ti_new;
	}
}


/// Updates indfac_int based on either adding individual in c (or removing) from the markov_eqn
void StateSpecies::update_indfac_int(unsigned int c, double t1, double t2, const Individual &ind, double sign, double &like_ch)
{
	auto me_ref = sp.comp_gl[c].me_ref;
	if(me_ref.size() == 0) return;
	
	auto t_start = timepoint[0];
	auto dt = details.dt;
	
	double dif, dLi;
	
	for(auto e : me_ref){
		auto &mev = markov_eqn_vari[e];
		const auto &me = sp.markov_eqn[e];
		auto indfac = sign*get_indfac(ind,me);
		auto &Li_m = Li_markov[e];
		
		vector <double> dLi_store;
		vector <double> dif_store;
		vector <double> ti_store;
		
		if(mev.time_vari == false){
			dif = (t2-t1)*indfac;
			auto &div = mev.div[0];
			div.indfac_int += dif;
			dLi = -div.value*dif;
			ti_store.push_back(0); dif_store.push_back(dif); dLi_store.push_back(dLi);
			Li_m[0] += dLi;
			like_ch += dLi;
		}
		else{
			auto ti = (unsigned int)(ALMOST_ONE*(t1-t_start)/dt);
			auto ti_new = (unsigned int)(ALMOST_ONE*(t2-t_start)/dt);
			if(ti_new == ti){
				dif = (t2-t1)*indfac;
				auto &div = mev.div[ti];
				div.indfac_int += dif;
				dLi = -div.value*dif;
				ti_store.push_back(ti); dif_store.push_back(dif); dLi_store.push_back(dLi);
				Li_m[ti] += dLi;
				like_ch += dLi;
			}
			else{
				dif = (timepoint[ti+1]-t1)*indfac;
				auto &div = mev.div[ti];
				div.indfac_int += dif;
				dLi = -div.value*dif;
				ti_store.push_back(ti); dif_store.push_back(dif); dLi_store.push_back(dLi);
				Li_m[ti] += dLi;
				like_ch += dLi;
			
				for(auto j = ti+1; j < ti_new; j++){
					dif = dt*indfac;
					auto &div = mev.div[j];
					div.indfac_int += dif;
					dLi = -div.value*dif;
					ti_store.push_back(j); dif_store.push_back(dif); dLi_store.push_back(dLi);
					Li_m[j] += dLi;
					like_ch += dLi;
				}
				
				dif = (t2-timepoint[ti_new])*indfac;
				auto &div2 = mev.div[ti_new];
				div2.indfac_int += dif;
				dLi = -div2.value*dif;
				ti_store.push_back(ti_new); dif_store.push_back(dif); dLi_store.push_back(dLi);
				Li_m[ti_new] += dLi;
				like_ch += dLi;
			}
		}
		
		back.push_back(Back(TI_INDFAC,e,ti_store));
		back.push_back(Back(DLI_INDFAC,dLi_store));
		back.push_back(Back(DIF_INDFAC,dif_store));
	}
}
		

/// Restores indfac_int
void StateSpecies::likelihood_indfac_int_restore(const vector <double> &store)
{	
	auto j = 0u;
	for(auto e = 0u; e < N; e++){
		auto &div = markov_eqn_vari[e].div;
		for(auto ti = 0u; ti < div.size(); ti++){
			div[ti].indfac_int = store[j]; j++;
		}
	}
}


/// Sets ind_trans in markov_eqn
void StateSpecies::likelihood_ind_trans()
{	
	for(auto e = 0u; e < N; e++){
		auto &div = markov_eqn_vari[e].div;
		for(auto ti = 0u; ti < div.size(); ti++){
			div[ti].ind_trans.clear();
		}
	}
	
	auto t_start = timepoint[0];
	
	for(auto i = 0u; i < individual.size(); i++){
		auto &ind = individual[i];
		
		for(auto k = 0u; k < ind.ev.size(); k++){
			auto &ev = ind.ev[k];
			if(ev.type == M_TRANS_EV){
				auto e = sp.tra_gl[ev.tr_gl].markov_eqn_ref;
				if(e == UNSET) emsg("SHould not be unset");
						
				auto &mev = markov_eqn_vari[e];
				auto dt = mev.dt;
				
				unsigned int ti;
				if(dt == UNSET) ti = 0;
				else ti = (unsigned int)((ev.t-t_start)/dt);
					
				DivIndRef dir; dir.i = i; dir.index = k;
				ev.m = e; ev.ti = ti; ev.index = mev.div[ti].ind_trans.size();
					
				mev.div[ti].ind_trans.push_back(dir);
			}
		}
	}
}


/// Calculate the likelihood of non-Markovian transitions
vector <double> StateSpecies::likelihood_nm_trans(unsigned int m, const vector <unsigned int> &list, const vector < vector <double> > &popnum_t, double &like_ch)
{
	vector <double> store;

	const auto &nmt = sp.nm_trans[m];
	
	const auto &ref = nmt.dist_param_eq_ref;
			
	vector <double> ref_val(ref.size());
	
	NmUpdate update = UPDATE_NO_TIME;
	for(auto i = 0u; i < ref.size(); i++){
		const auto &eq = eqn[ref[i]];
		if(eq.time_vari == true && update == UPDATE_NO_TIME) update = UPDATE_TIME;
		if(eq.ind_eff == true) update = UPDATE_IND;
	}
	
	if(update == UPDATE_NO_TIME){
		for(auto i = 0u; i < ref.size(); i++){
			ref_val[i] = eqn[ref[i]].calculate_param_only(param_val);
		}
	}
	
	auto &nm_ref = nm_trans_ev_ref[m];
	
	for(auto ti : list){
		auto &nmt_er = nm_ref[ti];
	
		if(nmt_er.size() > 0){	
			const auto &popnum = popnum_t[ti];
				
			if(update == UPDATE_TIME){
				for(auto i = 0u; i < ref.size(); i++){
					ref_val[i] = eqn[ref[i]].calculate(ti,popnum,param_val,spline_val);
				}
			}
		
			for(auto &nmti : nmt_er){
				auto &ind = individual[nmti.i];

				if(update == UPDATE_IND){
					for(auto i = 0u; i < ref.size(); i++){
						ref_val[i] = eqn[ref[i]].calculate_indfac(ind,ti,popnum,param_val,spline_val);
					}
				}
			
				auto &ev = ind.ev;
			
				auto e_end = nmti.e_end;
			
				auto &ev_end = ev[e_end];
				auto e_begin = ev_end.e_origin;
				if(e_begin == UNSET) emsg("SHould not be unset");
			
				auto t_trans = ev_end.t;
				auto t_from = ev[e_begin].t;
			
				auto dt = t_trans - t_from;
				if(dt <= 0) emsg("zero time4");

				auto Li_new = nm_trans_like(nmt.type,dt,ref_val);
				auto Li_old = ev_end.Li;
				store.push_back(Li_old);
				like_ch += Li_new-Li_old;
				ev_end.Li = Li_new; 
			}
		}
	}
	
	return store;
}


/// Restores nm_trans likelihood values
void StateSpecies::likelihood_nm_trans_restore(unsigned int m, const vector <unsigned int> &list, const vector <double> store) 
{
	auto k = 0u;
	auto &nm_ref = nm_trans_ev_ref[m];
	for(auto ti : list){
		auto &nmt_er = nm_ref[ti];
		for(auto &nmti : nmt_er){
			individual[nmti.i].ev[nmti.e_end].Li = store[k]; k++; 
		}
	}		
	
	if(k != store.size()) emsg("Wrong size");
}


/// Calculate the likelihood of non-Markovian transitions
vector <double> StateSpecies::likelihood_nm_trans_bp(unsigned int m, const vector <unsigned int> &list, const vector < vector <double> > &popnum_t, double &like_ch)
{
	vector <double> store;
	
	const auto &nmt = sp.nm_trans[m];

	auto bp_eq = nmt.bp_eq;
	
	switch(bp_eq){
	case UNSET:
		return store;
	
	case BP_FROM_OTHERS:
		{	
			auto &nm_ref = nm_trans_ev_ref[m];
			
			for(auto ti : list){
				if(nm_ref[ti].size() > 0){
					const auto &popnum = popnum_t[ti];
									
					for(auto &nmti : nm_ref[ti]){
						auto &ind = individual[nmti.i];
						auto prob = 1.0;
						for(auto e : nmt.bp_other_eq){
							prob -= eqn[e].calculate_indfac(ind,ti,popnum,param_val,spline_val);
						}
						auto &ev = ind.ev[nmti.e_end];
					
						auto Li_new = get_log_zero_one(prob);
						auto Li_old = ev.Li_bp;
						store.push_back(Li_old);
						like_ch += Li_new-Li_old;
						ev.Li_bp = Li_new;
					}
				}
			}
		}
		break;
		
	default:
		{	
			const auto &eq = eqn[bp_eq];
			
			NmUpdate update = UPDATE_NO_TIME;
			if(nmt.all_branches){
				for(auto e : nmt.bp_all_eq){
					const auto &eq = eqn[e];
					if(eq.time_vari == true && update == UPDATE_NO_TIME) update = UPDATE_TIME;
					if(eq.ind_eff == true) update = UPDATE_IND;
				}
			}
			else{
				if(eq.time_vari == true && update == UPDATE_NO_TIME) update = UPDATE_TIME;
				if(eq.ind_eff == true) update = UPDATE_IND;
			}
			
			auto Li_new = 0.0;
			if(update == UPDATE_NO_TIME){
				auto val = eq.calculate_param_only(param_val);
				if(nmt.all_branches){
					auto div = 0.0;
					for(auto e : nmt.bp_all_eq){
						div += eqn[e].calculate_param_only(param_val);
					}
				
					val /= div;
				}
				
				Li_new = get_log_zero_one(val);
			}
			
			auto &nmt_ref = nm_trans_ev_ref[m];
			
			for(auto ti : list){
				if(nmt_ref[ti].size() > 0){
					const auto &popnum = popnum_t[ti];
						
					if(update == UPDATE_IND){
						for(auto &nmti : nmt_ref[ti]){
							auto &ind = individual[nmti.i];
						
							auto val = eq.calculate_indfac(ind,ti,popnum,param_val,spline_val);
							
							if(nmt.all_branches){
								auto div = 0.0;
								for(auto e : nmt.bp_all_eq){
									div += eqn[e].calculate_indfac(ind,ti,popnum,param_val,spline_val);
								}
								val /= div;
							}
				
							auto &ev = ind.ev[nmti.e_end];
							
							Li_new = get_log_zero_one(val);
							auto Li_old = ev.Li_bp;
							store.push_back(Li_old);
							like_ch += Li_new-Li_old;
							ev.Li_bp = Li_new;
						}
					}
					else{
						if(update == UPDATE_TIME){
							auto val = eq.calculate(ti,popnum,param_val,spline_val);
							
							if(nmt.all_branches){
								auto div = 0.0;
								for(auto e : nmt.bp_all_eq){
									div += eqn[e].calculate(ti,popnum,param_val,spline_val);
								}
								val /= div;
							}
							
							Li_new = get_log_zero_one(val);
						}

						const auto nmti_old = nmt_ref[ti][0];
						auto Li_old =	individual[nmti_old.i].ev[nmti_old.e_end].Li_bp; 
						for(const auto &nmti : nmt_ref[ti]){
							store.push_back(Li_old);
							
							individual[nmti.i].ev[nmti.e_end].Li_bp = Li_new;
						}
						like_ch += (Li_new-Li_old)*nmt_ref[ti].size();
					}
				}
			}
		}
		break;
	}
	
	return store;
}


/// Restores Li_bp
void StateSpecies::likelihood_nm_trans_bp_restore(unsigned int m, const vector <unsigned int> &list, const vector <double> store)
{
	auto k = 0u;
	auto &nm_ref = nm_trans_ev_ref[m];
	for(auto ti : list){
		auto &nmt_er = nm_ref[ti];
		for(const auto &nmti : nmt_er){
			individual[nmti.i].ev[nmti.e_end].Li_bp = store[k]; k++; 
		}
	}		
	
	if(k != store.size()) emsg("Wrong size");
}


/// Calculate the likelihood of non-Markovian incomplete transitions
vector <double> StateSpecies::likelihood_nm_trans_incomp(unsigned int n, const vector <unsigned int> &list, const vector < vector <double> > &popnum_t, double &like_ch)
{
	vector <double> store;
	
	const auto &nmti = sp.nm_trans_incomp[n];
	auto B = nmti.nmtrans_ref.size();                 // The number of branches
	
	auto b_sel = UNSET;
	
	vector < vector <double> > ref_val;
	ref_val.resize(B);
	vector <double> bp_val(B);
	
	vector <NmUpdate> update(B,UPDATE_NO_TIME);
	auto update_bp = NO_UPDATE;
	for(auto b = 0u; b < B; b++){
		auto &up = update[b];
		
		auto m = nmti.nmtrans_ref[b];
		const auto &nmt = sp.nm_trans[m];
		
		const auto &ref = nmt.dist_param_eq_ref;
		ref_val[b].resize(ref.size());
		for(auto i = 0u; i < ref.size(); i++){
			const auto &eq = eqn[ref[i]];
			if(eq.time_vari == true && up == UPDATE_NO_TIME) up = UPDATE_TIME;
			if(eq.ind_eff == true) up = UPDATE_IND;
		}
			
		auto bp_eq = nmt.bp_eq;
		if(bp_eq == BP_FROM_OTHERS) b_sel = b;
		else{
			if(bp_eq != UNSET){
				if(update_bp == NO_UPDATE) update_bp = UPDATE_NO_TIME;
			
				if(nmti.all_branches){
					for(auto e : nmt.bp_all_eq){
						const auto &eq = eqn[e];
						if(eq.time_vari == true && update_bp == UPDATE_NO_TIME) update_bp = UPDATE_TIME;
						if(eq.ind_eff == true) update_bp = UPDATE_IND;
					}
				}
				else{
					const auto &eq = eqn[bp_eq];
					if(eq.time_vari == true && update_bp == UPDATE_NO_TIME) update_bp = UPDATE_TIME;
					if(eq.ind_eff == true) update_bp = UPDATE_IND;
				}
			}
		}
	}
	
	for(auto b = 0u; b < B; b++){
		if(update[b] == UPDATE_NO_TIME){
			auto m = nmti.nmtrans_ref[b];
			const auto &ref = sp.nm_trans[m].dist_param_eq_ref;
			for(auto i = 0u; i < ref.size(); i++){
				ref_val[b][i] = eqn[ref[i]].calculate_param_only(param_val);
			}
		}
	}
	
	if(update_bp == UPDATE_NO_TIME){		
		for(auto b = 0u; b < B; b++){
			if(b != b_sel){
				auto m = nmti.nmtrans_ref[b];
				auto bp_eq = sp.nm_trans[m].bp_eq;
				bp_val[b] = eqn[bp_eq].calculate_param_only(param_val);
			}
		}
		
		if(b_sel != UNSET){
			auto bp = 1.0; for(auto b = 0u; b < B; b++) if(b != b_sel) bp -= bp_val[b];
			bp_val[b_sel] = bp;
		}
		
		if(nmti.all_branches){
			auto div = 0.0;
			for(auto b = 0u; b < B; b++) div += bp_val[b];
			for(auto b = 0u; b < B; b++) bp_val[b] /= div;
		}
	}
	
	auto &nm_incomp_ref = nm_trans_incomp_ref[n];
	
	for(auto ti : list){
		auto &nmit_er = nm_incomp_ref[ti];
	
		const auto &popnum = popnum_t[ti];
		
		for(auto b = 0u; b < B; b++){
			if(update[b] == UPDATE_TIME){
				auto m = nmti.nmtrans_ref[b];
				const auto &ref = sp.nm_trans[m].dist_param_eq_ref;
				for(auto i = 0u; i < ref.size(); i++){
					ref_val[b][i] = eqn[ref[i]].calculate(ti,popnum,param_val,spline_val);
				}
			}
		}
		
		if(update_bp == UPDATE_TIME){	
			for(auto b = 0u; b < B; b++){
				if(b != b_sel){
					auto m = nmti.nmtrans_ref[b];			
					auto bp_eq = sp.nm_trans[m].bp_eq;
					bp_val[b] = eqn[bp_eq].calculate(ti,popnum,param_val,spline_val);
				}
			}
			
			if(b_sel != UNSET){
				auto bp = 1.0; for(auto b = 0u; b < B; b++) if(b != b_sel) bp -= bp_val[b];
				bp_val[b_sel] = bp;
			}
			
			if(nmti.all_branches){
				auto div = 0.0;
				for(auto b = 0u; b < B; b++) div += bp_val[b];
				for(auto b = 0u; b < B; b++) bp_val[b] /= div;
			}
		}
	
		for(auto i : nmit_er){
			auto &ind = individual[i];
			
			for(auto b = 0u; b < B; b++){
				if(update[b] == UPDATE_IND){	
					auto m = nmti.nmtrans_ref[b];
					const auto &ref = sp.nm_trans[m].dist_param_eq_ref;
					for(auto i = 0u; i < ref.size(); i++){
						ref_val[b][i] = eqn[ref[i]].calculate_indfac(ind,ti,popnum,param_val,spline_val);
					}
				}
			}
			
			if(update_bp == UPDATE_IND){	
				for(auto b = 0u; b < B; b++){
					if(b != b_sel){
						auto m = nmti.nmtrans_ref[b];			
						auto bp_eq = sp.nm_trans[m].bp_eq;			
						bp_val[b] = eqn[bp_eq].calculate_indfac(ind,ti,popnum,param_val,spline_val);
					}
				}
				
				if(b_sel != UNSET){
					auto bp = 1.0; for(auto b = 0u; b < B; b++) if(b != b_sel) bp -= bp_val[b];
					bp_val[b_sel] = bp;
				}
				
				if(nmti.all_branches){
					auto div = 0.0;
					for(auto b = 0u; b < B; b++) div += bp_val[b];
					for(auto b = 0u; b < B; b++) bp_val[b] /= div;
				}
			}
	
			const auto &ev = ind.ev;
			const auto &inmtr = ind.incomp_ref[nmti.cl];
			
			auto e_begin = inmtr.e_begin;	
			auto t_from = ev[e_begin].t;
			auto t_trans = inmtr.t_end;
			auto dt = t_trans - t_from;
			if(dt <= 0) emsg("zero time5");
			
			auto Li_new = nm_trans_incomp_full_like(nmti.nmtrans_ref,dt,ref_val,bp_val);
		
			auto &inmtref = ind.incomp_ref[nmti.cl];
			if(!inmtref.on) emsg("SHould be on");
			
			auto Li_old = inmtref.Li;
			store.push_back(Li_old);
			inmtref.Li = Li_new; 			
			like_ch += Li_new-Li_old;
		}
	}
	
	return store;
}


/// Restores nm_trans incomplete likelihood values
void StateSpecies::likelihood_nm_trans_incomp_restore(unsigned int n, const vector <unsigned int> &list, const vector <double> store) 
{
	const auto &nmti = sp.nm_trans_incomp[n];
	
	auto k = 0u;
	auto &nmi_ref = nm_trans_incomp_ref[n];
	for(auto ti : list){
		auto &nmit_er = nmi_ref[ti];
		for(auto i : nmit_er){
			individual[i].incomp_ref[nmti.cl].Li = store[k]; k++; 
		}
	}		
	
	if(k != store.size()) emsg("Wrong size");
}


/// Gets the likelihood for an incomplete transition
double StateSpecies::nm_trans_incomp_full_like(const vector <unsigned int> &nmtrans_ref, double dt, const vector< vector <double> > &ref_val, const vector <double> &bp_val) const
{
	auto B = nmtrans_ref.size();
	if(B == 1){
		const auto &nmt = sp.nm_trans[nmtrans_ref[0]];
		return nm_trans_incomp_like(nmt.type,dt,ref_val[0]);
	}
	else{
		auto p = 0.0;
		for(auto b = 0u; b < B; b++){
			auto bp = bp_val[b]; 
			check_bp(bp);
			
			const auto &nmt = sp.nm_trans[nmtrans_ref[b]];
			
			p += bp*nm_trans_incomp_like_no_log(nmt.type,dt,ref_val[b]);
		}
		
		return log(p);
	}
}


/// Gets the likelihood for a transition
double StateSpecies::nm_trans_incomp_like(TransType type, double dt, const vector <double> &ref_val) const
{		
	switch(type){
	case EXP_RATE: emsg("Should not be in NM"); return 0;
	case EXP_RATE_NM:	return exp_rate_upper_probability(dt,ref_val[0]); 
	case GAMMA: return gamma_upper_probability(dt,ref_val[0],ref_val[1]);
	case ERLANG: return gamma_upper_probability(dt,ref_val[0],sqrt(1.0/ref_val[1])); 
	case LOG_NORMAL: return lognormal_upper_probability(dt,ref_val[0],ref_val[1]);
	case WEIBULL:	return weibull_upper_probability(dt,ref_val[0],ref_val[1]);
	case PERIOD: if(dt < ref_val[0]) return 0; else return -LARGE;
	}
	return UNSET;
}


/// Gets the likelihood for a transition (with no log)
double StateSpecies::nm_trans_incomp_like_no_log(TransType type, double dt, const vector <double> &ref_val) const 
{		
	switch(type){
	case EXP_RATE: emsg("Should not be in NM"); return 0;
	case EXP_RATE_NM:	return exp_rate_upper_probability_no_log(dt,ref_val[0]); 
	case GAMMA: return gamma_upper_probability_no_log(dt,ref_val[0],ref_val[1]);
	case ERLANG: return gamma_upper_probability_no_log(dt,ref_val[0],sqrt(1.0/ref_val[1])); 
	case LOG_NORMAL: return lognormal_upper_probability_no_log(dt,ref_val[0],ref_val[1]);
	case WEIBULL:	return weibull_upper_probability_no_log(dt,ref_val[0],ref_val[1]);
	case PERIOD: return 1; 
	}
	return UNSET;
}


/// Gets the likelihood for a transition (with no log)
MeanSD StateSpecies::get_mean_sd(TransType type, const vector <double> &ref_val) const 
{		
	MeanSD msd;
	
	switch(type){
	case EXP_RATE: emsg("Should not be in NM"); break;
	case EXP_RATE_NM:	msd.mean = 1.0/ref_val[0]; msd.sd = sqrt(msd.mean); break; 
	case LOG_NORMAL: case GAMMA: msd.mean = ref_val[0]; msd.sd = ref_val[1]*msd.mean; break;
	case ERLANG: msd.mean = ref_val[0]; msd.sd = sqrt(1.0/ref_val[1])*msd.mean; break;
	case PERIOD: msd.mean = ref_val[0]; msd.sd = TINY; break;
	case WEIBULL:
		{
			auto lam = ref_val[0];
			auto k = ref_val[1];
			auto gam1 = tgamma(1+(1/k));
			auto gam2 = tgamma(1+(2/k));
			msd.mean = lam*gam1; 
			msd.sd = lam*sqrt(gam2-gam1*gam1);
		}
		break;	
	}
	
	return msd;
}


/// Works out how nm_trans likelihood changes as a result of ie changing
vector <NMupdate> StateSpecies::likelihood_ie_nm_trans_change(unsigned int i, unsigned int ie, double factor, const vector < vector <double> > &popnum_t, double &like_ch)
{
	vector <NMupdate> nm_update;
	
	const auto &ind = individual[i];
	const auto &event = ind.ev;

	for(auto e = 0u; e < event.size(); e++){
		const auto &ev = event[e];
	
		if(ev.type == NM_TRANS_EV){
			auto trg = ev.tr_gl;
			if(sp.tra_ie_affect[ev.tr_gl][ie]){
				auto t = ind.ev[ev.e_origin].t;
				auto m = ev.m;
				auto ti = ev.ti;
				
				const auto &nmt = sp.nm_trans[m];
				const auto &ref = nmt.dist_param_eq_ref;
			
				vector <double> ref_val(ref.size());
				for(auto i = 0u; i < ref.size(); i++){
					const auto &eq = eqn[ref[i]];
					ref_val[i] = eq.calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);	
					if(eq.ind_eff_exist[ie]) ref_val[i] *= factor;
				}
				
				auto Li = nm_trans_like(nmt.type,ev.t-t,ref_val);
				
				like_ch += Li-ev.Li;
				
				auto Li_bp = 0.0;
				
				if(sp.tra_gl[trg].branch){                 // Accounts for branching probability 
					auto bp_eq = nmt.bp_eq;
					double bp;
				
					if(bp_eq == BP_FROM_OTHERS){
						bp = 1.0;
						for(auto e : nmt.bp_other_eq){
							const auto &eq = eqn[e];
							auto value = eq.calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
							if(eq.ind_eff_exist[ie]) value *= factor;
							bp -= value;
						}
					}
					else{	
						const auto &eq = eqn[bp_eq];
						bp = eq.calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
						if(eq.ind_eff_exist[ie]) bp *= factor;
							
						if(nmt.all_branches){
							auto div = 0.0;
							for(auto e: nmt.bp_all_eq){
								const auto &eq = eqn[e];
								auto val = eq.calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
								if(eq.ind_eff_exist[ie]) val *= factor;
								div += val; 
							}
							bp /= div;
						}
					}
					
					Li_bp = get_log_zero_one(bp);
		
					like_ch += Li_bp-ev.Li_bp; 
				}	
				
				NMupdate nm; 
				nm.type = NM_EVENT_STORE; nm.i = i; nm.e = e; //nm.index = index;
				nm.Li = Li; nm.Li_bp = Li_bp; 
				nm_update.push_back(nm);
			}
		}
	}
	
	for(auto j = 0u; j < ind.incomp_ref.size(); j++){
		const auto &inmt = ind.incomp_ref[j];
		if(inmt.on){
			auto n = inmt.n;
			auto ti = inmt.ti;
			
			const auto &nmt = sp.nm_trans_incomp[n];
			
			auto dt = inmt.t_end - event[inmt.e_begin].t; if(dt <= 0) emsg("zero time");

			auto val = get_nm_incomp_val_ie_factor(ie,factor,nmt,ti,ind,popnum_t);
			
			auto Li = nm_trans_incomp_full_like(nmt.nmtrans_ref,dt,val.ref_val,val.bp_val);
		
			like_ch += Li-inmt.Li; 
					
			NMupdate nm; 
			nm.type = NM_INCOMP_STORE; nm.i = i; nm.e = j; //nm.index = index;
			nm.Li = Li;
			nm_update.push_back(nm);	
		}
	}
	
	return nm_update;
}


/// Updates nm_trans likelihoods if ie proposal is accepted
void StateSpecies::likelihood_ie_nm_trans_update(const vector <NMupdate> &nm_st)
{
	for(const auto &nm : nm_st){
		switch(nm.type){
		case NM_EVENT_STORE:
			{
				auto &ev = individual[nm.i].ev[nm.e];
				ev.Li = nm.Li;
				ev.Li_bp = nm.Li_bp;
			}
			break;
			
		case NM_INCOMP_STORE:
			{
				individual[nm.i].incomp_ref[nm.e].Li = nm.Li;
			}
			break;
		}
	}
}


/// Calculates individual effect group likelihood
double StateSpecies::likelihood_indeff_group(unsigned int g, double &like_ch)
{
	const auto &ieg = sp.ind_eff_group[g];
	
	auto Li = 0.0;
	
	const auto &omega_inv = ind_eff_group_sampler[g].omega_inv;
	
	auto log_det = determinant_fast(omega_inv);
	
	auto N = ieg.list.size();
		
	if(ieg.A_matrix.set == false){
		vector <double> samp(N);
		
		for(auto &ind : individual){
			for(auto j = 0u; j < N; j++){
				auto k = ieg.list[j].index;
				if(k >= ind.ie.size()) emsg("Out of range10");
				samp[j] = ind.ie[k];
			}
			
			auto sum = 0.0;
			for(auto j = 0u; j < N; j++){
				for(auto i = 0u; i < N; i++){
					sum += samp[j]*omega_inv[j][i]*samp[i];
				}
			}
			
			Li += 0.5*log_det - 0.5*sum;
		}
	}
	else{
		const auto &A_inv_nonzero = ind_eff_group_sampler[g].A_inv_nonzero;
	
		auto I = individual.size();
		
		vector < vector <double> > samp;
	
		samp.resize(I);
		for(auto i = 0u; i < I; i++){
			const auto &ind = individual[i];
			samp[i].resize(N);
			for(auto j = 0u; j < N; j++){
				auto k = ieg.list[j].index;
				if(k >= ind.ie.size()) emsg("Out of range10");
				samp[i][j] = ind.ie[k];
			}
		}
		
		if(false){
			for(auto j = 0u; j < I; j++){
				cout << individual[j].name << ": ";
				for(const auto &nz : A_inv_nonzero[j]){
					auto i = nz.i; 
					auto A_inv = nz.value;
					cout << individual[i].name << " " << A_inv << ", ";
				}
				cout << "j" << endl;
			}
			emsg("P");
		}
		
		auto sum = 0.0;
		for(auto j = 0u; j < I; j++){
			for(const auto &nz : A_inv_nonzero[j]){
				auto i = nz.i; 
				auto A_inv = nz.value;
			 
				for(auto jj = 0u; jj < N; jj++){
					for(auto ii = 0u; ii < N; ii++){
						sum += samp[j][jj]*omega_inv[jj][ii]*A_inv*samp[i][ii];
					}
				}
			}
		}
		
		Li += 0.5*I*log_det - 0.5*sum;
	}
	
	auto store = Li_ie[g];
	like_ch += Li-store;

	Li_ie[g] = Li;
	
	return store;
}


/// Calculates the change in Li_ie from a single individual effect changing
double StateSpecies::likelihood_indeff_change(unsigned int i, unsigned int ie, double value_new, double value_old)
{
	const auto &inde = sp.ind_effect[ie];
	auto g = inde.index;
	auto num = inde.num;
	
	const auto &ieg = sp.ind_eff_group[g];
	const auto &omega_inv = ind_eff_group_sampler[g].omega_inv;
		
	auto N = ieg.list.size();
	
	vector <unsigned int> list;	
	for(auto k = 0u; k < N; k++) list.push_back(ieg.list[k].index);
		
	auto sum = 0.0;
	
	if(ieg.A_matrix.set == false){                   // Works out the change in Li_ie
		const auto &ind = individual[i];
		for(auto k = 0u; k < N; k++){
			if(k == num) sum += omega_inv[num][num]*(value_new*value_new-value_old*value_old);
			else sum += 2*omega_inv[num][k]*(value_new-value_old)*ind.ie[list[k]];
		}
	}
	else{
		const auto &A_inv_nonzero = ind_eff_group_sampler[g].A_inv_nonzero;
	
		for(const auto &nz : A_inv_nonzero[i]){
			auto j = nz.i; 
			auto A_inv = nz.value;
			 
			for(auto jj = 0u; jj < N; jj++){
				if(j == i && jj == num){
					sum += omega_inv[num][num]*A_inv*(value_new*value_new-value_old*value_old);
				}
				else{
					sum += 2*omega_inv[jj][num]*A_inv*(value_new-value_old)*individual[j].ie[list[jj]];
				}
			}
		}
	}
	
	return -0.5*sum;
}


/// Calculates the likelihood for the population transition numbers
void StateSpecies::likelihood_pop(const vector < vector <double> > &popnum_t)
{
	for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
		auto &Li = Li_markov_pop[tr];
		auto &tn = trans_num[tr];
		auto &tnm = tnum_mean_st[tr];
		for(auto ti = 0u; ti < T; ti++){
			tnm[ti] = calculate_tnum_mean(ti,tr,popnum_t[ti],cpop_st[ti],param_val,dtimepoint[ti]);
			
			Li[ti] = poisson_probability(tn[ti],tnm[ti]);
		}
	}
}


/// Updates a section of likelihood
void StateSpecies::likelihood_pop_section(unsigned int tr, unsigned int ti1, unsigned int ti2, const vector < vector <double> > &popnum_t, double &like_ch)
{
	auto &Li = Li_markov_pop[tr];
	auto &tn = trans_num[tr];
	auto &tnm = tnum_mean_st[tr];
	
	for(auto ti = ti1; ti < ti2; ti++){
		tnm[ti] = calculate_tnum_mean(ti,tr,popnum_t[ti],cpop_st[ti],param_val,dtimepoint[ti]);
		like_ch -= Li[ti];
		Li[ti] = poisson_probability(tn[ti],tnm[ti]);
		like_ch += Li[ti];
	}
}


/// Updates the transition number and corresponding likelihood
double StateSpecies::Li_update_tn(unsigned int tr, unsigned int ti, int sign)
{
	auto &Li = Li_markov_pop[tr][ti];
	auto &tn = trans_num[tr][ti];
	auto &tnm = tnum_mean_st[tr][ti];
	
	tn += sign;
	
	auto Li_new = poisson_probability(tn,tnm);	
	auto d = Li_new-Li;	
	Li = Li_new;
	
	return d;
}


/// Updates the likelihood from change in c
double StateSpecies::Li_update_c(unsigned int c, int ma, unsigned int ti, unsigned int ti_next, const vector < vector <double> > &popnum_t)
{
	auto d = 0.0;
	for(auto tr : sp.comp_gl[c].tr_leave){
		auto &Li = Li_markov_pop[tr];
		auto &tn = trans_num[tr];
		auto &tnm = tnum_mean_st[tr];
		for(auto tii = ti; tii < ti_next; tii++){
			auto c_aft = cpop_st[tii][c];
			auto c_bef = c_aft - ma;
			if(c_bef != 0){                              // This is a faster way to calculate
				tnm[tii] *= cpop_st[tii][c]/c_bef;
			}
			else{
				tnm[tii] = calculate_tnum_mean(tii,tr,popnum_t[tii],cpop_st[tii],param_val,dtimepoint[tii]);
			}
			
			auto Li_new = poisson_probability(tn[tii],tnm[tii]);	
			d += Li_new-Li[tii];	
			Li[tii] = Li_new;
		}
	}
	
	return d;
}


/// Works out the change in the likelihood due to population change
vector <double> StateSpecies::likelihood_pop_change(unsigned int tr, const vector <unsigned int> &list, const vector < vector <double> > &popnum_t, double &like_ch)
{
	vector <double> store;

	auto &Li = Li_markov_pop[tr];
	auto &tn = trans_num[tr];
	auto &tnm = tnum_mean_st[tr];
	for(auto ti : list){
		store.push_back(tnm[ti]);
		tnm[ti] = calculate_tnum_mean(ti,tr,popnum_t[ti],cpop_st[ti],param_val,dtimepoint[ti]);
		
		auto old_val = Li[ti];
		auto new_val = poisson_probability(tn[ti],tnm[ti]);
		store.push_back(old_val);
		
		Li[ti] = new_val;
		
		like_ch += new_val - old_val;
	}
	
	return store;
}


/// Restores likelihood change from population
void StateSpecies::likelihood_pop_change_restore(unsigned int tr, const vector <unsigned int> &list, const vector <double> &vec)
{
	auto j = 0.0;
	auto &Li = Li_markov_pop[tr];
	auto &tnm = tnum_mean_st[tr];
	for(auto ti : list){
		tnm[ti] = vec[j]; j++;
		Li[ti] = vec[j]; j++;
	}
}


/// Calculate the likelihood of the initial condition
vector <double> StateSpecies::likelihood_init_cond(double &like_ch)
{
	vector <double> store;

	store.push_back(Li_init_cond);
		
	if(mode == SIM) return store;
	
	if(sp.init_cond.type != INIT_POP_DIST) return store;
	
	const auto &ic = sp.init_cond;
	const auto &icv = init_cond_val;
	
	auto foc_cl = ic.focal_cl;
	
	auto Li = 0.0;
	
	if(foc_cl == UNSET){
		if(sp.type == INDIVIDUAL){
			Li += factorial(icv.N_total_unobs);
			Li += multinomial_samp_probability(icv.cnum,icv.frac);
		}
		else{
			Li += multinomial_probability(icv.N_total,icv.cnum,icv.frac);
		}
	}
	else{                                            // A focal classification is set
		const auto &claa = sp.cla[foc_cl];
		if(sp.type == INDIVIDUAL){
			for(auto c = 0u; c < claa.ncomp; c++){
				Li += factorial(icv.N_focal_unobs[c]);
				Li += multinomial_samp_probability(icv.cnum_reduce[c],icv.frac_comb);
			}
		}
		else{
			for(auto c = 0u; c < claa.ncomp; c++){
				Li += multinomial_probability(icv.N_focal[c],icv.cnum_reduce[c],icv.frac_comb);
			}
		}
	}

	like_ch += Li-Li_init_cond;
	Li_init_cond = Li;

	return store;
}


/// Gets the values for the incomplete nm trans
NMIncompVal StateSpecies::get_nm_incomp_val(const NMTransIncomp &nmti, unsigned int ti, const Individual &ind, const vector < vector <double> > &popnum_t) const
{
	const auto &nmtrans_ref = nmti.nmtrans_ref;
	
	NMIncompVal value;
	auto &ref_val = value.ref_val;
	auto &bp_val = value.bp_val;
	
	auto B = nmtrans_ref.size();
	
	ref_val.resize(B);
	bp_val.resize(B,UNSET);
	
	auto b_sel = UNSET;
	for(auto b = 0u; b < B; b++){
		auto m = nmtrans_ref[b];
		
		const auto &ref = sp.nm_trans[m].dist_param_eq_ref;
		ref_val[b].resize(ref.size());
		for(auto i = 0u; i < ref.size(); i++){
			const auto &eq = eqn[ref[i]];
			if(eq.ind_eff == true) ref_val[b][i] = eq.calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
			else{
				if(eq.time_vari) ref_val[b][i] = eq.calculate(ti,popnum_t[ti],param_val,spline_val);
				else{
					ref_val[b][i] = eq.calculate_param_only(param_val);
				}
			}
		}
	
		auto bp_eq = sp.nm_trans[m].bp_eq;
		if(bp_eq != UNSET){
			if(bp_eq == BP_FROM_OTHERS) b_sel = b;
			else{
				const auto &eq = eqn[bp_eq];
				if(eq.ind_eff == true) bp_val[b] = eq.calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
				else{
					if(eq.time_vari) bp_val[b] = eq.calculate(ti,popnum_t[ti],param_val,spline_val);
					else bp_val[b] = eq.calculate_param_only(param_val);
				}
			}
		}
	}
	
	if(b_sel != UNSET){
		auto bp = 1.0;
		for(auto b = 0u; b < B; b++) if(b != b_sel) bp -= bp_val[b];
		bp_val[b_sel] = bp;
	}
	
	if(nmti.all_branches){
		auto div = 0.0;
		for(auto b = 0u; b < B; b++) div += bp_val[b];
		for(auto b = 0u; b < B; b++) bp_val[b] /= div;
	}
	
	return value;
}


/// Prints information about incomplete transition
void StateSpecies::print_incomp_val(string te, NMIncompVal value) const
{
	cout << te << endl;
	for(auto rval : value.ref_val){
		for(auto va : rval) cout << va << ", ";
	}
	cout << "trans val" << endl;
	for(auto va : value.bp_val) cout << va << ", " ;
	cout << "bp val" << endl;
	cout << endl;
}
		

/// Gets the values for the incomplete nm trans with individual effect changed by a factor
NMIncompVal StateSpecies::get_nm_incomp_val_ie_factor(unsigned int ie, double factor, const NMTransIncomp &nmti, unsigned int ti, const Individual &ind, const vector < vector <double> > &popnum_t) const
{
	const auto &nmtrans_ref = nmti.nmtrans_ref;
	 
	NMIncompVal value;
	auto &ref_val = value.ref_val;
	auto &bp_val = value.bp_val;
	
	auto B = nmtrans_ref.size();
	ref_val.resize(B);
	bp_val.resize(B,UNSET);
	
	auto b_sel = UNSET;
	for(auto b = 0u; b < B; b++){
		auto m = nmtrans_ref[b];
		
		const auto &ref = sp.nm_trans[m].dist_param_eq_ref;
		ref_val[b].resize(ref.size());
		for(auto i = 0u; i < ref.size(); i++){
			const auto &eq = eqn[ref[i]];
			ref_val[b][i] = eq.calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
			if(eq.ind_eff_exist[ie]) ref_val[b][i] *= factor;
		}
	
		auto bp_eq = sp.nm_trans[m].bp_eq;
		if(bp_eq != UNSET){
			if(bp_eq == BP_FROM_OTHERS) b_sel = b;
			else{
				const auto &eq = eqn[bp_eq];
				bp_val[b] = eq.calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
				if(eq.ind_eff_exist[ie]) bp_val[b] *= factor;
			}
		}
	}
			
	if(b_sel != UNSET){
		auto bp = 1.0;
		for(auto b = 0u; b < B; b++) if(b != b_sel) bp -= bp_val[b];
		bp_val[b_sel] = bp;
	}
	
	if(nmti.all_branches){
		auto div = 0.0;
		for(auto b = 0u; b < B; b++) div += bp_val[b];
		for(auto b = 0u; b < B; b++) bp_val[b] /= div;
	}
	
	return value;
}


/// Calculates the likelihood for observed nm events (used for simulation proposals) 
double StateSpecies::nm_obs_dprob(const Individual &ind) const
{
	const auto &event = ind.ev;
	
	auto dprob = 0.0;
	for(const auto &ev : event){	
		if(!ev.observed && ev.type == NM_TRANS_EV){
			dprob += ev.Li + ev.Li_bp;
		}
	}
	
	for(const auto &inmtr : ind.incomp_ref){
		if(inmtr.on) dprob += inmtr.Li;
	}

	return dprob;
}


/// Calculates the likelihood for observed nm events (used for simulation proposals) 
double StateSpecies::nm_single_obs_dprob(unsigned int cl, const Individual &ind) const
{
	const auto &event = ind.ev;
	
	auto dprob = 0.0;
	for(const auto &ev : event){
		if(!ev.observed && ev.type == NM_TRANS_EV && ev.cl == cl){
			dprob += ev.Li + ev.Li_bp;
		}
	}
	
	const auto &inmtr = ind.incomp_ref[cl]; 
	if(inmtr.on) dprob += inmtr.Li;//nmtir.Li;
	
	return dprob;
}


/// Used to sum up probability of an event
double StateSpecies::sum_markov_prob(double t1, double t2, unsigned int c, unsigned int tr_gl, unsigned int i, vector < vector <double> > &en) const
{
	auto tr_gl_cor = tr_gl;
	if(c != UNSET){
		tr_gl_cor = sp.tra_gl[tr_gl].transform[c];
	}
	
	const auto &tra_cor = sp.tra_gl[tr_gl_cor];
	
	auto &entr = en[tr_gl_cor];
	
	if(tra_cor.type != EXP_RATE) emsg("Should be exp_rate");
	
	auto e = tra_cor.markov_eqn_ref;
	
	const auto &me = sp.markov_eqn[e];
	const auto &mev = markov_eqn_vari[e];
	
	auto sum = 0.0;

	auto ti1 = get_ti(t1);
	auto ti2 = get_ti(t2);
	if(ti2 == entr.size()) ti2--;
	
	auto dt = details.dt;
	auto t_start = details.t_start;
	
	auto fac = 1.0;
	if(me.ind_variation){
		const auto &ind = individual[i];
		for(auto ie : me.ind_eff_mult) fac *= ind.exp_ie[ie];
		for(auto fe : me.fix_eff_mult) fac *= ind.exp_fe[fe];
	}
	
	for(auto ti = ti1; ti <= ti2; ti++){
		auto dtt = dt;
		if(ti == ti1) dtt = t_start+(ti+1)*dt - t1;
		else{
			if(ti == ti2) dtt = t2 - (t_start+ti*dt);
		}

		auto tii = ti; if(!mev.time_vari) tii = 0;
		auto val = fac*mev.div[tii].value*dtt;
		if(ti >= entr.size()) emsg("Prob");

		entr[ti] += val;
		sum -= val;
	}
 
	return sum;
}


/// Calculates the estimated and actual number of transitions in each division for each global transition
void StateSpecies::calc_trans_diag(ParticleSpecies &ps, const vector < vector <double> > &popnum_t) const
{
	auto T = timepoint.size()-1;
		
	auto &cpd = ps.cum_prob_dist;
	auto &en = ps.exp_num;
	auto t_end = details.t_end;
	auto dt = details.dt;
	
	cpd.resize(sp.tra_gl.size());
	for(auto tr = 0u; tr < sp.tra_gl.size(); tr++) cpd[tr].resize(H_BIN,0);
	
	switch(type){
	case POPULATION:
		en = tnum_mean_st;
		
		for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
			auto &tn = trans_num[tr];
			auto &tnm = tnum_mean_st[tr];
			for(auto ti = 0u; ti < T; ti++){
				auto n = tn[ti];
				auto lam = tnm[ti];
				auto pr = 1-(ran()*exp(poisson_probability(n,lam)) + poisson_upper_probability_no_log(n,lam));
				
				auto b = (unsigned int)(ALMOST_ONE*pr*H_BIN);
				if(b >= H_BIN) emsg("hbin range4");
							
				cpd[tr][b]++;
			}
		}
		break;
		
	case INDIVIDUAL:
		{
			en.resize(sp.tra_gl.size());
			for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
				en[tr].resize(T,0);
			}
			
			auto ncla = sp.ncla;

			// Accounts for sources
			for(auto m = 0u; m < sp.markov_eqn.size(); m++){
				const auto &me = sp.markov_eqn[m];
				if(me.source){
					for(auto ti = 0u; ti < T; ti++){
						const auto &mev = markov_eqn_vari[m];
						auto tii = ti; if(mev.time_vari == false) tii = 0;
						auto val = mev.div[tii].value*dt;
						for(auto tr : me.source_tr_gl){
							en[tr][ti] += val;
						}
					}
				}
			}

			for(auto i = 0u; i < individual.size(); i++){
				vector <unsigned int> origin(ncla,0);
				
				vector <double> cum_pr;
				
				const auto &ind = individual[i];
				
				const auto &eve = ind.ev;
				for(auto e = 0u; e < eve.size(); e++){
					const auto &ev = eve[e];
					
					switch(ev.type){
					case NM_TRANS_EV:
						{
							auto cl = ev.cl;
							auto tr_gl = ev.tr_gl;
							if(ev.e_origin != origin[cl]) emsg("e_origin wrong");
						
							const auto &ev_or = eve[ev.e_origin];
							auto tr_gl_cor = sp.tra_gl[tr_gl].transform[ev_or.c_after];
							
							add_nm_prob(tr_gl_cor,ev_or.t,ps,popnum_t,ind);
							add_nm_cpd(tr_gl_cor,ev_or.t,ev.t,ps,popnum_t,ind);
							
							origin[cl] = e;
						}
						break;
						
					case M_TRANS_EV:
						{
							auto cl = ev.cl;
							auto tr_gl = ev.tr_gl;
							
							auto sum = 0.0;
							for(auto ee = origin[cl]; ee < e; ee++){
								sum += sum_markov_prob(eve[ee].t,eve[ee+1].t,eve[ee].c_after,tr_gl,i,en);
							}
							
							auto pr = exp(sum);
							auto b = (unsigned int)(ALMOST_ONE*pr*H_BIN);
							if(b >= H_BIN) emsg("hbin range2");
							
							cpd[tr_gl][b]++;
							
							origin[cl] = e;
						}
						break;
			
					default: break;
					}
				}
				
				if(eve.size() > 0){
					// Works out when there is no event
					for(auto cl = 0u; cl < ncla; cl++){
						auto e_st = origin[cl];
						const auto &ev_or = eve[e_st];
						auto c = ev_or.c_after;
						if(c != UNSET){
							const auto &cgl = sp.comp_gl[c];
							
							const CompGlTransGroup &tlg = cgl.tra_leave_group[cl];
							if(tlg.markov){
								for(auto tr_gl : tlg.tr_list){
									auto sum = 0.0;
									for(auto ee = e_st; ee < eve.size(); ee++){
										double tnext;
										if(ee < eve.size()-1) tnext = eve[ee+1].t;
										else tnext = t_end;
										
										sum += sum_markov_prob(eve[ee].t,tnext,eve[ee].c_after,tr_gl,i,en);
									}
									
									auto pr = ran()*exp(sum);
									auto b = (unsigned int)(ALMOST_ONE*pr*H_BIN);
									if(b >= H_BIN) emsg("hbin range3");
									
									cpd[tr_gl][b]++;
								}
							}	
							else{
								const auto &tr_list = tlg.tr_list;
								
								auto B = tr_list.size();
								auto t1 = ev_or.t;
								
								if(B == 1){
									add_nm_prob(tlg.tr_list[0],t1,ps,popnum_t,ind);
									add_nm_cpd(tlg.tr_list[0],t1,UNSET,ps,popnum_t,ind);
								}
								else{
									auto ti = get_ti(t1);
			
									auto sum = 0.0;
									vector <double> sum_store;
									
									auto sum2 = 0.0;
									vector <double> sum_store2;
									
									vector <double> pr_st(B);
									
									for(auto b = 0u; b < B; b++){
										const auto &tra = sp.tra_gl[tr_list[b]];
			
									
										auto m = tra.nm_trans_ref;
										const auto &nmt = sp.nm_trans[m];
										const auto &ref = nmt.dist_param_eq_ref;

										vector <double> ref_val(ref.size());
										for(auto i = 0u; i < ref.size(); i++){
											ref_val[i] = eqn[ref[i]].calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
										}
			
										auto bp = eqn[tra.bp.eq_ref].calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
										pr_st[b] = nm_trans_incomp_like_no_log(nmt.type,t_end-t1,ref_val);
											
										sum += bp;
										sum_store.push_back(sum);
										
										sum2 += bp*pr_st[b];
										sum_store2.push_back(sum2);
									}
				
									// For future predition of transions rate the branch is randomly selected 
									{
										auto z = sum*ran();
										auto k = 0u; while(k < B && z > sum_store[k]) k++;
										if(k == B) emsg("Should not be B");	
										add_nm_prob(tlg.tr_list[k],t1,ps,popnum_t,ind);								
									}
									
									// For updating cumulative probability distribution, the branch must be 
									// selected according to bp multiplied probability of being branch
									{
										auto z = sum2*ran();
										auto k = 0u; while(k < B && z > sum_store2[k]) k++;
										if(k == B) emsg("Should not be B");		
										
										auto pr = ran()*pr_st[k];
										auto b = (unsigned int)(ALMOST_ONE*pr*H_BIN);
										if(b >= H_BIN) emsg("hbin range1");
			
										cpd[tlg.tr_list[k]][b]++;
									}
								}
							}
						}
					}
				}
			}
		}
		break;
	}
}


/// Adds probabilities for a non-Markovian transition
void StateSpecies::add_nm_prob(unsigned tr_gl_or, double t1, ParticleSpecies &ps, const vector < vector <double> > &popnum_t, const Individual &ind) const
{						
	auto &en = ps.exp_num;
	auto t_end = details.t_end;
	
	auto m = sp.tra_gl[tr_gl_or].nm_trans_ref;
	const auto &nmt = sp.nm_trans[m];
	const auto &ref = nmt.dist_param_eq_ref;

	auto ti = get_ti(t1);
	
	vector <double> ref_val(ref.size());
	for(auto i = 0u; i < ref.size(); i++){
		ref_val[i] = eqn[ref[i]].calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
	}

	if(nmt.type == PERIOD){
		auto t = t1 + ref_val[0];
		if(t < t_end){
			auto tii = get_ti(t);
			en[tr_gl_or][tii] += 1;
		}
		return;
	}

	auto msd = get_mean_sd(nmt.type,ref_val);
	auto tmin = t1 + msd.mean - 3*msd.sd; 
	if(tmin < t1) tmin = t1;
	if(tmin > t_end) tmin = t_end;
	auto tmax = t1 + msd.mean + 4*msd.sd;
	if(tmax > t_end) tmax = t_end;
	//tmin = t1; tmax = t_end;
	auto ti_min = get_ti(tmin);
	auto ti_max = get_ti(tmax);
	if(ti_min < ti+1) ti_min = ti+1;
	if(ti_max > T) ti_max = T;

	auto step = (unsigned int)(ti_max-ti_min)/20;
	if(step == 0) step = 1;
	//step = 1;
	
	// Does the first fraction of a division
	auto prob = 1.0;
	auto tot = 0.0;
	bool first = true;
	while(ti < ti_max){
		unsigned int tii;
		
		if(first){ tii = ti+1; first = false;}
		else{
			if(ti < ti_min) tii = ti_min;
			else{
				tii = ti + step; if(tii > ti_max) tii = ti_max;
			}
		}
		
		auto pr = nm_trans_incomp_like_no_log(nmt.type,timepoint[tii]-t1,ref_val);
		auto num = (prob-pr)/(tii-ti);
	
		while(ti < tii){
			tot += num;
			en[tr_gl_or][ti] += num;
			ti++;
		}
		
		prob = pr;
	}
}


/// Adds to cumulative probability distribution
void StateSpecies::add_nm_cpd(unsigned tr_gl_or, double t1, double t2, ParticleSpecies &ps, const vector < vector <double> > &popnum_t, const Individual &ind) const
{	
	auto &cpd = ps.cum_prob_dist;
	auto t_end = details.t_end;
	
	auto m = sp.tra_gl[tr_gl_or].nm_trans_ref;
	const auto &nmt = sp.nm_trans[m];
	const auto &ref = nmt.dist_param_eq_ref;

	auto ti = get_ti(t1);

	vector <double> ref_val(ref.size());
	for(auto i = 0u; i < ref.size(); i++){
		ref_val[i] = eqn[ref[i]].calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
	}
	
	double pr;
	if(t2 != UNSET){
		pr = nm_trans_incomp_like_no_log(nmt.type,t2-t1,ref_val);
	}
	else{
		pr = ran()*nm_trans_incomp_like_no_log(nmt.type,t_end-t1,ref_val);
	}
	
	auto b = (unsigned int)(ALMOST_ONE*pr*H_BIN);
	if(b >= H_BIN) emsg("hbin range1");
	
	cpd[tr_gl_or][b]++;
}


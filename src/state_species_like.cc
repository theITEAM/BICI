// Calculates the likelihood

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
vector <double> StateSpecies::likelihood_markov(unsigned int e, vector <unsigned int> list, AllTimeStatus as, double &like_ch)
{
	auto &me = sp.markov_eqn[e];
	auto &me_vari = markov_eqn_vari[e];
	auto &Li_mark = Li_markov[e];

	if(as == ALL_DIV){
		if(list.size() != 0) emsg("list pr");
		for(auto ti = 0u; ti < me_vari.div.size(); ti++) list.push_back(ti);
	}
	
	vector <double> store;
	for(auto ti : list){
		const auto &div = me_vari.div[ti];
		
		auto val = div.value;
		auto Li = -val*div.indfac_int;
	
		for(const auto &it : div.ind_trans){
			Li += log(get_indfac(individual[it.i],me)*(val+LOG_THRESH));			
		}
	
		like_ch += Li-Li_mark[ti];
		store.push_back(Li_mark[ti]);
		
		Li_mark[ti] = Li;
	}
	
	return store;
}


/// Calculate the value of a markov equation
vector <double> StateSpecies::likelihood_markov_value(unsigned int e, vector <unsigned int> list, AllTimeStatus as, const vector < vector <double> > &popnum_t)
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
		if(as == ALL_DIV){
			if(list.size() != 0) emsg("list pr");
			for(auto ti = 0u; ti < me_vari.div.size(); ti++) list.push_back(ti);
		}
		
		for(auto ti : list){	
			auto value = eqn[me.eqn_ref].calculate(ti,popnum_t[ti],param_val,spline_val);
			store.push_back(me_vari.div[ti].value);
			me_vari.div[ti].value = value;
		}
	}
	
	return store;
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
		
		auto val = 0.0; if(sp.markov_eqn[e].source == true) val = T/div.size();
		for(auto ti = 0u; ti < div.size(); ti++){
			store.push_back(div[ti].indfac_int);
			div[ti].indfac_int = val;
		}
	}
	
	for(auto i = 0u; i < individual.size(); i++){
		const auto &ind = individual[i];
		
		auto t = t_start;
		auto ti = 0u;
		
		auto c = ind.cinit;
		for(auto k = 0u; k <= ind.ev.size(); k++){
			auto t_ev = t_end; if(k < ind.ev.size()) t_ev = ind.ev[k].t;
	
			auto ti_new = (unsigned int)(ALMOST_ONE*(t_ev-t_start)/dt);
								
			if(c < CODE){
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
	
	auto c = ind.cinit;
	auto dLi_sum = 0.0;
	
	double va, dLi;
			
	for(auto k = 0u; k <= ev.size(); k++){
		auto t_ev = t_end; if(k < ind.ev.size()) t_ev = ev[k].t;

		auto ti_new = (unsigned int)(ALMOST_ONE*(t_ev-t_start)/dt);
							
		if(c < CODE){
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
				const auto &mef = e.markov_eq_ref;
				const auto &me = sp.markov_eqn[mef.e];
				
				for(auto ie_me : me.ind_eff_mult){
					if(ie_me == ie){
						dLi = log(factor);
						Li_markov[mef.e][mef.ti] += dLi;
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
	
	auto c = ind.cinit;
	
	auto j = 0u;
	for(auto k = 0u; k <= ev.size(); k++){
		auto t_ev = t_end; if(k < ind.ev.size()) t_ev = ev[k].t;
		auto ti_new = (unsigned int)(ALMOST_ONE*(t_ev-t_start)/dt);
						
		if(c < CODE){
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
				const auto &mef = e.markov_eq_ref;
				const auto &me = sp.markov_eqn[mef.e];
				for(auto ie_me : me.ind_eff_mult){
					if(ie_me == ie){	
						Li_markov[mef.e][mef.ti] -= log(factor);
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
			
		if(mev.time_vari == false){
			dif = (t2-t1)*indfac;
			auto &div = mev.div[0];
			div.indfac_int += dif;
			dLi = -div.value*dif;
			Li_markov[e][0] += dLi;
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
				Li_markov[e][ti] += dLi;
				like_ch += dLi;
			}
			else{
				dif = (timepoint[ti+1]-t1)*indfac;
				auto &div = mev.div[ti];
				div.indfac_int += dif;
				dLi = -div.value*dif;
				Li_markov[e][ti] += dLi;
				like_ch += dLi;
			
				for(auto j = ti+1; j < ti_new; j++){
					dif = dt*indfac;
					auto &div = mev.div[j];
					div.indfac_int += dif;
					dLi = -div.value*dif;
					Li_markov[e][j] += dLi;
					like_ch += dLi;
				}
				
				dif = (t2-timepoint[ti_new])*indfac;
				auto &div2 = mev.div[ti_new];
				div2.indfac_int += dif;
				dLi = -div2.value*dif;
				Li_markov[e][ti_new] += dLi;
				like_ch += dLi;
			}
		}
	}
}
		
		
// Removes events from markov eqn for old event sequence 
void StateSpecies::remove_ind_events_markov_eqn(unsigned int i, double &like_ch)
{
	const auto &ind = individual[i];
	
	for(auto ee = 0u; ee < ind.ev.size(); ee++){
		auto &ev = ind.ev[ee];
		
		if(ev.type == M_TRANS_EV){
			auto &mer = ev.markov_eq_ref;
			auto e = mer.e;
			auto ti = mer.ti;
			auto &it = markov_eqn_vari[e].div[ti].ind_trans;
			auto k = mer.index;
			if(k != it.size()-1){
				auto &it_shift = it[it.size()-1];
				it[k] = it_shift;
				individual[it_shift.i].ev[it_shift.index].markov_eq_ref.index = k;
			}
			it.pop_back();
		
			auto &me = sp.markov_eqn[e];
			auto val = markov_eqn_vari[e].div[ti].value;
			
			auto Li = log(get_indfac(ind,me)*(val+LOG_THRESH));
		
			Li_markov[e][ti] -= Li;
			like_ch -= Li;
		}
	}
}


// Adds events from markov eqn for old event sequence 
void StateSpecies::add_ind_events_markov_eqn(unsigned int i, double &like_ch)
{
	auto &ind = individual[i];
	
	auto t_start = timepoint[0];
	auto dt = details.dt;
	
	for(auto ee = 0u; ee < ind.ev.size(); ee++){
		auto &ev = ind.ev[ee];
				
		if(ev.type == M_TRANS_EV){
			// Add individual on to ind_trans and sets .markov_eq_ref for event
			auto e = sp.tra_gl[ev.tr_gl].markov_eqn_ref;
			auto &mev = markov_eqn_vari[e];
			
			auto ti = 0; if(mev.time_vari == true) ti = (unsigned int)((ev.t-t_start)/dt);
					
			auto &mer = ev.markov_eq_ref;
			auto &it = mev.div[ti].ind_trans;	
			mer.e = e; mer.ti = ti; mer.index = it.size();
				
			DivIndRef dir; dir.i = i; dir.index = ee;
			it.push_back(dir);
		
			auto &me = sp.markov_eqn[e];
			auto val = markov_eqn_vari[e].div[ti].value;
			
			auto Li = log(get_indfac(ind,me)*(val+LOG_THRESH));
			Li_markov[e][ti] += Li;
			like_ch += Li;
		}
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

/// Sets  ind_trans in markov_eqn
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
				IndDivRef idr; idr.e = e; idr.ti = ti; idr.index = mev.div[ti].ind_trans.size();
					
				mev.div[ti].ind_trans.push_back(dir);
				ev.markov_eq_ref = idr;
			}
			else ev.markov_eq_ref.e = UNSET;
		}
	}
}


/// Calculate the likelihood of non-Markovian transitions
double StateSpecies::likelihood_nm_trans(unsigned int i, const vector < vector <double> > &popnum_t, double &like_ch)
{
	auto t_start = details.t_start;
	auto timestep = details.dt;
	
	const auto &nmt = sp.nm_trans[i];
	const auto &nmt_ind = nm_trans_ind[i];
	
	// TO DO braching probability
	
	auto Li = 0.0;
	for(const auto &nmti : nmt_ind){
		const auto &ind = individual[nmti.i];
		auto &cl_ev = ind.cl_ev[nmti.cl];
		auto e = nmti.index;
		auto t_trans = ind.ev[cl_ev[e].index_ev].t;
		if(t_trans != nmti.t) emsg("Times do not agree");
		
		double t_from;
		if(e > 0) t_from = ind.ev[cl_ev[e-1].index_ev].t;
		else t_from = t_start;
		
		auto dt = t_trans - t_from;
		if(dt <= 0) emsg("zero time");
		
		// TO DO make faster by doing a non-time varying version
		
		auto ti = (unsigned int)((t_from-t_start)/timestep);
		const auto &popnum = popnum_t[ti];
		
		const auto &ref = nmt.dist_param_eq_ref;
		
		switch(nmt.type){
		case EXP_RATE: case EXP_MEAN: emsg("Should not be in NM"); break;
		case GAMMA:
			{
				auto mean = eqn[ref[0]].calculate_indfac(ind,ti,popnum,param_val,spline_val);
				auto cv = eqn[ref[1]].calculate_indfac(ind,ti,popnum,param_val,spline_val);
				Li += gamma_probability(dt,mean,cv);
			}
			break;
			
		case ERLANG:
			{
				auto mean = eqn[ref[0]].calculate_indfac(ind,ti,popnum,param_val,spline_val);
				auto shape = eqn[ref[1]].calculate_indfac(ind,ti,popnum,param_val,spline_val);
				Li += gamma_probability(dt,mean,sqrt(1.0/shape));
			}
			break;

		case LOG_NORMAL:
			{
				auto mean = eqn[ref[0]].calculate_indfac(ind,ti,popnum,param_val,spline_val);
				auto cv = eqn[ref[1]].calculate_indfac(ind,ti,popnum,param_val,spline_val);
				Li += lognormal_probability(dt,mean,cv);
			}
			break;
			
		case WEIBULL:
			{
				auto scale = eqn[ref[0]].calculate_indfac(ind,ti,popnum,param_val,spline_val);
				auto shape = eqn[ref[1]].calculate_indfac(ind,ti,popnum,param_val,spline_val);
				Li += weibull_probability(dt,shape,scale);
			}
			break;
		}
	}

	auto store = Li_nm_trans[i];
	like_ch = Li - store;
	Li_nm_trans[i] = Li;

	return store;
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
	
	if(ieg.A_matrix.set == false){   // Works out the change in Li_ie
		for(auto k = 0u; k < N; k++){
			if(k == num) sum += omega_inv[num][num]*(value_new*value_new-value_old*value_old);
			else sum += 2*omega_inv[num][k]*(value_new-value_old)*individual[i].ie[list[k]];
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


/// Works out the change in the likelihood
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


/// Works out the change in the likelihood
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


/// Works out the observation probability
vector <double> StateSpecies::likelihood_observation(double &like_ch)
{
	vector <double> store;
	
	if(mode == SIM){ Li_obs = 0; return store;}
	
	const auto &popd = sp.pop_data;
	const auto &ptransd = sp.pop_trans_data;
	
	auto Li = 0.0;
	
	switch(type){
	case POPULATION:
		{
			for(auto i = 0u; i < popd.size(); i++){
				store.push_back(pop_data_num[i]);
			}
			for(auto i = 0u; i < ptransd.size(); i++){
				store.push_back(pop_trans_data_num[i]);
			}
			
			for(auto i = 0u; i < popd.size(); i++){
				const auto &pd = popd[i];
				const auto &filt = pd.filt;
				
				auto ti = sp.get_ti(pd.t);
				if(ti == T) ti--;
				
				auto num = 0.0;
				for(auto c = 0u; c < sp.comp_gl.size(); c++){
					if(filt[c] == true){
						num += cpop_st[ti][c];
					}
				}
			
				pop_data_num[i] = num;
			}
			
			for(auto i = 0u; i < ptransd.size(); i++){
				pop_trans_data_num[i] = 0;
			}
			
			for(auto ti = 0u; ti < T; ti++){
				for(auto trg = 0u; trg < sp.tra_gl.size(); trg++){
					auto num = trans_num[trg][ti];
					if(num > 0){
						for(auto ref : sp.pop_trans_ref[ti][trg]){
							pop_trans_data_num[ref] += num;
						}
					}
				}
			}
		}
		break;
	
	case INDIVIDUAL:
		{	
			// Calculates the values for any equations associated with observation (e.g. Se, Sp)
			
			for(auto i = 0u; i < obs_eqn_value.size(); i++){
				obs_eqn_value[i] = eqn[sp.obs_eqn[i]].calculate_param_only(param_val);
			}
			
			for(auto i = 0u; i < popd.size(); i++) pop_data_num[i] = 0;
			for(auto i = 0u; i < ptransd.size(); i++) pop_trans_data_num[i] = 0;
			
			auto t_end = details.t_end;
			
			for(auto i = 0u; i < individual.size(); i++){
				const auto &obs = sp.individual[i].obs;
				const auto &ind = individual[i];
				
				auto c = ind.cinit;
				auto m = 0u;
				auto ev = ind.ev;
				
				for(auto k = 0u; k <= ev.size(); k++){
					double t;	
					if(k < ev.size()) t = ev[k].t; else t = t_end;
					
					while(m < obs.size() && obs[m].t < t){
						const auto &ob = obs[m];
						
						switch(ob.type){
						case OBS_TRANS_EV:  // Checks to see if event agrees with observation
							{
							}
							break;
							
						case OBS_COMP_EV:
							{
								if(sp.comp_gl[c].cla_comp[ob.cl] != ob.c) Li += LI_WRONG;
							}
							break;
							
						case OBS_TEST_EV:
							{
								Li += like_diag_test(c,ob);
							}
							break;
							
						case OBS_POP:
							{
								if(c < CODE){
									if(sp.pop_data[ob.ref].filt[c] == true){
										pop_data_num[ob.ref]++;
									}
								}
							}
							break;
						}
						
						m++;
					}
				
					if(k < ev.size()){
						const auto &e = ev[k];
						auto tr_gl = e.tr_gl;
						
						if(e.type == M_TRANS_EV){
							auto ti = sp.get_ti(e.t);
						
							if(sp.trans_not_allow[ti][tr_gl] == true){
								if(m < obs.size() && obs[m].t == e.t && sp.source[obs[m].so].trans_filt[tr_gl]){}
								else Li += LI_WRONG;
							}
							
							for(auto ref : sp.pop_trans_ref[ti][tr_gl]){
								pop_trans_data_num[ref]++;
							}
						}
						c = e.c_after;
					}	
				}
			}
		}
		break;
	}
	
	for(auto i = 0u; i < pop_data_num.size(); i++){
		const auto &popd = sp.pop_data[i];
		Li += normal_probability(pop_data_num[i],popd.value,popd.sd);
	}
	
	for(auto i = 0u; i < pop_trans_data_num.size(); i++){
		const auto &poptd = sp.pop_trans_data[i];
		Li += normal_probability(pop_trans_data_num[i],poptd.value,poptd.sd);
	}
		
	like_ch += Li-Li_obs;
	
	store.push_back(Li_obs);
	
	Li_obs = Li;

	return store;
}


/// Resores values
void StateSpecies::likelihood_observation_restore(const vector <double> &store)
{
	auto j = 0u;
	
	if(type == POPULATION){
		const auto &popd = sp.pop_data;
		const auto &ptransd = sp.pop_trans_data;
	
		for(auto i = 0u; i < popd.size(); i++){
			pop_data_num[i] = store[j]; j++;
		}
		for(auto i = 0u; i < ptransd.size(); i++){
			pop_trans_data_num[i] = store[j]; j++;
		}
	}
		
	Li_obs = store[j]; j++;
}


/// Gets the probability of a diagnostic test
double StateSpecies::like_diag_test(unsigned int c, const ObsData &ob) const
{
	if(c > CODE) return LI_WRONG;
	else{
		const auto &om = sp.source[ob.so].obs_model;
		const auto &cgl = sp.comp_gl[c];
		
		auto cl = om.diag_test_sens.cl;
		if(om.diag_test_sens.comp[cgl.cla_comp[cl]] == true){ // Truely infected
			auto Se = obs_eqn_value[ob.Se_obs_eqn];
			if(ob.test_res == true) return log(Se+LOG_THRESH); 
			else return log(1-Se+LOG_THRESH);
		}
		else{
			auto Sp = obs_eqn_value[ob.Sp_obs_eqn];
			if(ob.test_res == false) return log(Sp+LOG_THRESH); 
			else return log(1-Sp+LOG_THRESH);
		}
	}
}


/// Updates pop_trans and works out the change in observation likelihood
double StateSpecies::update_pop_trans(double t, unsigned int tr_gl, int sign)
{
	auto ti = sp.get_ti(t);
		
	auto dLi = 0.0; 
	for(auto ref : sp.pop_trans_ref[ti][tr_gl]){
		const auto &poptd = sp.pop_trans_data[ref];
	
		dLi -= normal_probability(pop_trans_data_num[ref],poptd.value,poptd.sd);
		pop_trans_data_num[ref] += sign;
		dLi += normal_probability(pop_trans_data_num[ref],poptd.value,poptd.sd);
	}
	
	return dLi;
}

// This code relates to updating the state based on a change in a given individual

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
// NOTES
// markov_eqn list all transitions which are Markovian. Attached to this list are the individuals 

using namespace std;

#include "state.hh"
#include "utils.hh"

/// Updates an individual with a new event sequence
Like State::update_ind(unsigned int p, unsigned int i, unsigned int cinit_new, vector <Event> &ev_new)
{
	Like like_ch;
	
	auto &sp = model.species[p];
	auto &ssp = species[p];
	
	auto &ind = ssp.individual[i];
	const auto &obs = sp.individual[i].obs;
	
	//ssp.print_event("start",ind);
	
	ssp.remove_ind_events_markov_eqn(i,like_ch.markov);
	
	auto c_old = ind.cinit;
	auto &ev_old = ind.ev;
	auto c_new = cinit_new;
	
	auto e_old = 0u;
	auto e_old_max = ev_old.size();
	auto e_new = 0u;
	auto e_new_max = ev_new.size();
	
	double t_old, t_new;
	
	const auto &tp = model.timepoint;
	auto t_start = tp[0];
	auto t_end = tp[T];
	auto dt = model.details.dt;

	auto m = 0u;
	
	auto dLi_obs = 0.0;
	
	auto t = tp[0];
	auto ti = 0u; 
	do{
		if(e_old < e_old_max) t_old = ev_old[e_old].t; else t_old = t_end;
		if(e_new < e_new_max) t_new = ev_new[e_new].t; else t_new = t_end;
		
		auto t_next = t_old; if(t_new < t_next) t_next = t_new;
	
		auto ti_next = 1+(unsigned int)(ALMOST_ONE*(t_next-t_start)/dt);
	
		
		while(m < obs.size() && obs[m].t < t_next){
			const auto &ob = obs[m];
			
			switch(ob.type){
			case OBS_TRANS_EV:  // Checks to see if event agrees with observation
				{
					//if(k < ev.size() && ev[k].t == t && sp.source[ob.so].trans_filt[ev[k].tr_gl] == true){}
					//else Li += LI_WRONG;
				}
				break;
				
			case OBS_COMP_EV:
				{
					if(c_new != c_old){
						if(sp.comp_gl[c_old].cla_comp[ob.cl] != ob.c) dLi_obs -= LI_WRONG;
						if(sp.comp_gl[c_new].cla_comp[ob.cl] != ob.c) dLi_obs += LI_WRONG;
					}
				}
				break;
				
			case OBS_TEST_EV:
				{
					if(c_new != c_old){
						dLi_obs += ssp.like_diag_test(c_new,ob) - ssp.like_diag_test(c_old,ob);
					}
				}
				break;
				
			case OBS_POP:
				{
					if(c_new != c_old){
						auto ref = ob.ref;
						const auto &popd = sp.pop_data[ref];
						const auto &filt = popd.filt;
						
						auto dnum = 0;
						if(c_old < CODE && filt[c_old] == true) dnum--;
						if(c_new < CODE && filt[c_new] == true) dnum++;
					
						if(dnum != 0){
							dLi_obs -= normal_probability(ssp.pop_data_num[ref],popd.value,popd.sd);
							ssp.pop_data_num[ob.ref] += dnum;
							dLi_obs += normal_probability(ssp.pop_data_num[ref],popd.value,popd.sd);
						}
					}
				}
				break;
			}
			
			m++;
		}
	
	/*
			if(k < ev.size()){
				const auto &e = ev[k];
				auto tr_gl = e.tr_gl;
				
				for(auto na : not_allow){
					if(sp.source[na].trans_filt[tr_gl] == true) Li += LI_WRONG;
				}

				for(auto tpo : trans_pop_on){
					if(sp.pop_trans_data[tpo].filt[tr_gl] == true) pop_trans_data_num[tpo]++;
				}
				c = e.c_after;
			}	
		}
	*/
	
	
		if(c_new != c_old && ti_next != ti){
			// Updates populations
			if(c_old < CODE){
				for(auto k : sp.comp_gl[c_old].pop_ref_simp){
					if(pop_map[k] == UNSET){ pop_list.push_back(k); pop_map[k] = -1;}
					else pop_map[k]--;
				}
			}
		
			if(c_new < CODE){
				for(auto k : sp.comp_gl[c_new].pop_ref_simp){
					if(pop_map[k] == UNSET){ pop_list.push_back(k); pop_map[k] = 1;}
					else pop_map[k]++;
				}
			}
		
			for(auto k : pop_list){
				auto ma = pop_map[k];
				const auto &po = model.pop[k];
				
				if(ma != 0){
					auto num = 1.0; 
					for(auto ie : po.ind_eff_mult) num *= ind.exp_ie[ie];
					for(auto fe : po.fix_eff_mult) num *= ind.exp_fe[fe];
					num *= ma;
		
					for(auto tii = ti; tii < ti_next; tii++) popnum_t[tii][k] += num;
				}
				
				for(const auto &mef : po.markov_eqn_ref){
					if(markov_eqn_map[mef.p][mef.e] == false){
						markov_eqn_map[mef.p][mef.e] = true;
						markov_eqn_list.push_back(mef);
					}
				}
				
				for(const auto &mtf : po.trans_ref){
					if(trans_map[mtf.p][mtf.tr] == false){
						trans_map[mtf.p][mtf.tr] = true;
						trans_list.push_back(mtf);
					}
				}
				
				pop_map[k] = UNSET;
			}		
			pop_list.clear();
			
			update_pop_change(ti,ti_next,like_ch.markov);
		}
				
		if(c_new != c_old && t != t_next){
			if(c_old < CODE) ssp.update_indfac_int(c_old,t,t_next,ind,-1,like_ch.markov);
			if(c_new < CODE) ssp.update_indfac_int(c_new,t,t_next,ind,1,like_ch.markov);
		}
		
		if(t_next == t_end) break;
		
		t = t_next; ti = ti_next;
		
		auto ti = sp.get_ti(t);
		
		if(t_next == t_old){   // Removes old event
			const auto &e = ev_old[e_old];
			
			if(e.type == M_TRANS_EV){
				if(sp.trans_not_allow[ti][e.tr_gl] == true){
					if(m < obs.size() && obs[m].t == t && sp.source[obs[m].so].trans_filt[e.tr_gl]){}
					else dLi_obs -= LI_WRONG;
				}
			
				dLi_obs += ssp.update_pop_trans(t_next,e.tr_gl,-1);
			}
			
			c_old = e.c_after; e_old++;
		}
		else{                  // Adds new event
			const auto &e = ev_new[e_new];

			if(e.type == M_TRANS_EV){
				if(sp.trans_not_allow[ti][e.tr_gl] == true){
					if(m < obs.size() && obs[m].t == t && sp.source[obs[m].so].trans_filt[e.tr_gl]){}
					else dLi_obs += LI_WRONG;
				}
	
				dLi_obs += ssp.update_pop_trans(t_next,e.tr_gl,1);
			}
			
			c_new = e.c_after;
			e_new++;
		}
	}while(true);
	
	like_ch.obs += dLi_obs;
	ssp.Li_obs += dLi_obs;
						
	ev_old = ev_new;  // Copies event sequence to new
	
	ssp.add_ind_events_markov_eqn(i,like_ch.markov);
	
	//ssp.print_event("end",ind);
	
	return like_ch;
}





/// Updates the likelihoods based on a set of population changes 
void State::update_pop_change(unsigned int ti, unsigned int ti_next, double &like_ch)
{
	/// Updates markov_eqn as a result of populations change in size
	for(const auto &mef : markov_eqn_list){
		auto pp = mef.p, ee = mef.e;
		auto &spp = model.species[pp];
		if(spp.type != INDIVIDUAL) emsg("must be ind");
	
		const auto &me = spp.markov_eqn[ee];
		auto &sspp = species[pp];
		auto &me_vari = sspp.markov_eqn_vari[ee];

		const auto &eq = model.eqn[me.eqn_ref];
		
		if(me_vari.time_vari == false) emsg("must have time variation");
		
		//for(auto tii = ti+1; tii <= ti_next; tii++){
		for(auto tii = ti; tii < ti_next; tii++){
			auto value = eq.calculate(tii,popnum_t[tii],param_val,spline_val);	
			
			auto &div = me_vari.div[tii];
			
			auto val_old = div.value;
			auto dLi = -(value-val_old)*div.indfac_int;
			auto n = div.ind_trans.size();
			if(n > 0){
				dLi += n*log((value+LOG_THRESH)/(val_old+LOG_THRESH));
			}
			sspp.Li_markov[ee][tii] += dLi;
			like_ch += dLi;
			div.value = value;
		}
		
		markov_eqn_map[pp][ee] = false;
	}
	markov_eqn_list.clear();
	
	/// Updates transitions in population species as a result of population change in size
	for(const auto &mtf : trans_list){
		auto pp = mtf.p, tr = mtf.tr;
		if(model.species[pp].type != POPULATION) emsg("must be ind");
		
		species[pp].likelihood_pop_section(tr,ti,ti_next,popnum_t,like_ch);
		
		trans_map[pp][tr] = false;
	}
	trans_list.clear();
}



/// Initilaises all the maps used to update quantities fast
void State::initialise_update_ind_maps()
{
	pop_map.resize(model.pop.size(),UNSET);
		
	markov_eqn_map.resize(model.nspecies);
	trans_map.resize(model.nspecies);
	for(auto p = 0u; p < model.nspecies; p++){
		switch(model.species[p].type){
		case INDIVIDUAL: markov_eqn_map[p].resize(model.species[p].markov_eqn.size(),false); break;
		case POPULATION: trans_map[p].resize(model.species[p].tra_gl.size(),false);
		}
	}
}


/// Updates an individual with a new event sequence
void State::update_ie_population(unsigned int p, unsigned int i, unsigned int ie, double dif, double &like_ch)
{
	const auto &tp = model.timepoint;
	auto t_start = tp[0];
	auto t_end = tp[T];
	auto dt = model.details.dt;
	
	auto &sp = model.species[p];
	auto &ssp = species[p];
	
	auto &ind = ssp.individual[i];
	const auto &ev = ind.ev;
	
	auto c = ind.cinit;
	auto ti = 0u;
	
	for(auto k = 0u; k <= ev.size(); k++){
		auto t_ev = t_end; if(k < ind.ev.size()) t_ev = ev[k].t;

		auto ti_next = 1+(unsigned int)(ALMOST_ONE*(t_ev-t_start)/dt);
							
		if(c < CODE){
			for(auto j : sp.comp_gl[c].pop_ref_simp){
				const auto &po = model.pop[j];
				for(auto iepop : po.ind_eff_mult){
					if(iepop == ie){
						auto num = dif; 
						for(auto fe : po.fix_eff_mult) num *= ind.exp_fe[fe];
			
						for(auto tii = ti; tii < ti_next; tii++) popnum_t[tii][j] += num;
					
						for(const auto &mef : po.markov_eqn_ref){
							if(markov_eqn_map[mef.p][mef.e] == false){
								markov_eqn_map[mef.p][mef.e] = true;
								markov_eqn_list.push_back(mef);
							}
						}
				
						for(const auto &mtf : po.trans_ref){
							if(trans_map[mtf.p][mtf.tr] == false){
								trans_map[mtf.p][mtf.tr] = true;
								trans_list.push_back(mtf);
							}
						}
					}
				}
			}
			
			update_pop_change(ti,ti_next,like_ch);
		}	
		
		if(k < ev.size()){
			c = ev[k].c_after;
			ti = ti_next;
		}
	}
}



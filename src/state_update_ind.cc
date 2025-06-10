// This code relates to updating the state based on a change in a given individual

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
using namespace std;

#include "state.hh"
#include "utils.hh"

/// Updates an individual with a new event sequence
// In the case of non-Markovian transitions m, ti, index and e_origin are automatically set  
// In the case of Markovian transitions m, index automatically set  
Like State::update_ind(unsigned int p, unsigned int i, vector <Event> &ev_ne, UpdateType type)
{
	timer[IND_TIMER] -= clock();
	
	if(type == UP_SINGLE) back_init();
	
	auto &sp = model.species[p];
	auto &ssp = species[p];
	
	auto &ind = ssp.individual[i];
	
	if(false){
		cout << " change" << endl;	
		ssp.print_event(ind.ev);
		ssp.print_event(ev_ne);
	}
	
	// Deals with event changes
	Like like_ch;
	
	vector <PopUpdate> pop_update;
	
	auto add_ev = ssp.update_ind(i,ev_ne,popnum_t,pop_update,like_ch);

	timer[IND_POP_UPDATE_TIMER] -= clock();
	
	// Deals with population changes
	for(const auto &pup : pop_update){
		auto ti = pup.ti, ti_next = pup.ti_next;
		auto c_old = pup.c_old, c_new = pup.c_new;
		
		// Updates populations
		if(c_old != UNSET){
			for(auto k : sp.comp_gl[c_old].pop_ref_simp){
				if(pop_map[k] == UNSET){ pop_list.push_back(k); pop_map[k] = -1;}
				else pop_map[k]--;
			}
		}
	
		if(c_new != UNSET){
			for(auto k : sp.comp_gl[c_new].pop_ref_simp){
				if(pop_map[k] == UNSET){ pop_list.push_back(k); pop_map[k] = 1;}
				else pop_map[k]++;
			}
		}
	
		vector <PopChange> pop_change;
		for(auto k : pop_list){
			auto ma = pop_map[k];
			const auto &po = model.pop[k];
			
			if(ma != 0){
				double num = ma; 
				for(auto ie : po.ind_eff_mult) num *= ind.exp_ie[ie];
				for(auto fe : po.fix_eff_mult) num *= ind.exp_fe[fe];
				
				back_pop.push_back(BackPop(POP_NUM_T,ti,ti_next,k,num));

				PopChange po_ch; po_ch.po = k; po_ch.num = num; 
				pop_change.push_back(po_ch);
				
				for(auto tii = ti; tii < ti_next; tii++) popnum_t[tii][k] += num;
				
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
				
			pop_map[k] = UNSET;
		}		
		pop_list.clear();

		if(false){
			cout << ti << " " << ti_next << " pop chang" << endl;
			for(auto pc : pop_change) cout << pc.po << " " << pc.num << "pc" << endl; 
		}
		
		update_pop_change(ti,ti_next,pop_change,like_ch.markov);
	}
	
	timer[IND_POP_UPDATE_TIMER] += clock();
	
	// Adds in new events
	for(auto e : add_ev) ssp.add_event_ref(i,e,popnum_t,like_ch);
	
	if(false) ssp.print_event("end",ind);
	
	timer[IND_TIMER] += clock();
	
	return like_ch;
}


/// Updates the likelihoods based on a set of population changes 
void State::update_pop_change(unsigned int ti, unsigned int ti_next, const vector <PopChange> &pop_change, double &like_ch)
{
	// This is the part which is slow for spatial models
	// Updates markov_eqn as a result of populations change in size
	for(const auto &mef : markov_eqn_list){
		auto pp = mef.p, ee = mef.e;
		species[pp].recalc_markov_value(ee,ti,ti_next,popnum_t,pop_change,like_ch);
		markov_eqn_map[pp][ee] = false;
	}
	markov_eqn_list.clear();
				
	// Updates transitions in population species as a result of population change in size
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


/// Updates an individual with a new individual effect value
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
	
	const auto &exp_fe = ind.exp_fe;
	const auto &exp_ie = ind.exp_ie;
	
	auto c = UNSET;
	auto ti = 0u;
	
	for(auto k = 0u; k <= ev.size(); k++){
		auto t_ev = t_end; if(k < ind.ev.size()) t_ev = ev[k].t;

		auto ti_next = (unsigned int)(ALMOST_ONE+(t_ev-t_start)/dt);
							
		if(c != UNSET){
			vector <PopChange> pop_change;
			for(auto j : sp.comp_gl[c].pop_ref_simp){
				const auto &po = model.pop[j];
				const auto &ind_eff_mult = po.ind_eff_mult;
				for(auto iepop : ind_eff_mult){
					if(iepop == ie){
						auto num = dif; 
						for(auto fe : po.fix_eff_mult) num *= exp_fe[fe];
						if(ind_eff_mult.size() > 1){
							for(auto iee : ind_eff_mult){ if(iee != ie) num *= exp_ie[iee];}
						}
						
						PopChange po_ch; po_ch.po = j; po_ch.num = num; 
						pop_change.push_back(po_ch);
				
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
		
			update_pop_change(ti,ti_next,pop_change,like_ch);
		}	
		
		if(k < ev.size()){
			c = ev[k].c_after;
			ti = ti_next;
		}
	}
}


/// Looks to update iif.w and corresponding change in markov likelihood
void State::update_ie_trans_tree(unsigned int p, unsigned int i, unsigned int ie, double ratio, double &like_ch)
{
	auto &ssp = species[p];
	
	auto &ind = ssp.individual[i];
	const auto &ev = ind.ev;
	
	const auto &inf_node = genetic_value.inf_node;
	for(auto eve : ev){
		auto node = eve.inf_node_ref;
		
		if(node != UNSET){
			const auto &in = inf_node[node];
			for(const auto &iev : in.inf_ev){
				if(iev.type == INFECT_OTHER){
					const auto &in2 = inf_node[iev.index];
					auto &ssp2 = species[in2.p];
					auto &ev = ssp2.individual[in2.i].ev[in2.e];
					auto &iif = ev.ind_inf_from;
					
					const auto &po = model.pop[iif.po];
					if(find_in(po.ind_eff_mult,ie) != UNSET){
						iif.w *= ratio;
						
						auto dLi = log(ratio);
						like_ch += dLi;
						const auto &sp2 = model.species[in2.p];

						auto ti = sp2.get_ti(ev.t);
						const auto &tr = sp2.tra_gl[ev.tr_gl];
						
						ssp2.Li_markov[tr.markov_eqn_ref][ti] += dLi;
					}
				}
			}
		}
	}
}


/// Based in IndInfFrom works out what the population weight is for infecting individual
double State::get_w_from_indinffrom(const IndInfFrom &iif) const 
{
	const auto &po = model.pop[iif.po];
	const auto &ind = species[iif.p].individual[iif.i];

	auto num = 1.0; 
	if(po.ind_variation){
		for(auto ie : po.ind_eff_mult) num *= ind.exp_ie[ie];
		for(auto fe : po.fix_eff_mult) num *= ind.exp_fe[fe];
	}
	
	return num;
}

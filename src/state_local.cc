// This generates local population-based proposals to the state

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
using namespace std;

#include "state.hh"
#include "utils.hh"

/// Initialises local_map (used for local proposals)
void State::initialise_local_map()
{
	comp_local_map.resize(model.species.size());
	pop_trans_map.resize(model.species.size());
	pop_data_map.resize(model.species.size());
	for(auto p = 0u; p < model.species.size(); p++){
		const auto &sp = model.species[p];
	
		comp_local_map[p].resize(sp.comp_gl.size(),UNSET);
		pop_trans_map[p].resize(sp.pop_trans_data.size(),UNSET);
		pop_data_map[p].resize(sp.pop_data.size(),UNSET);
	}
}


/// Used to local changes
bool LocalChange_ord (LocalChange lc1, LocalChange lc2)                      
{ return (lc1.ti < lc2.ti); };  


/// Calculates the likelihood change associated with a local change 
Like State::calculate_local_change(unsigned int p, vector <LocalChange> &local_change, bool &ill)
{
	Like like_ch;
	
	back_init();
	
	const auto &sp = model.species[p];
	auto &ssp = species[p];
	
	// Makes sure change is legal
	ill = false;
	for(const auto &lc : local_change){
		switch(lc.type){
		case REMOVE_EVENT:
			if(ssp.trans_num[lc.tr][lc.ti] == 0){ ill = true; return like_ch;}
			break;

		case ADD_EVENT:
			break;
		
		default: break;
		}
	}
	
	if(false){
		for(auto &lc : local_change){
			switch(lc.type){
			case ADD_EVENT:
				cout << "Add event " << sp.tra_gl[lc.tr].name << " " << lc.ti << endl; 
				break;
				
			case REMOVE_EVENT:
				cout << "Remove event " << sp.tra_gl[lc.tr].name << " " << lc.ti << endl; 
				break;
				
			case ADD_IC:
				cout << "Add IC from transition" << sp.tra_gl[lc.tr].name << " " << lc.ti << endl; 
				break;
				
			case REMOVE_IC:
				cout << "Remove IC from trantision " << sp.tra_gl[lc.tr].name << " " << lc.ti << endl; 
				break;
				
			case ADD_C_IC:
				cout << "Add IC from compartment " << sp.comp_gl[lc.tr].name << " " << lc.ti << endl; 
				break;
				
			case REMOVE_C_IC:
				cout << "Remove IC from compartment" << sp.comp_gl[lc.tr].name << " " << lc.ti << endl; 
				break;
			}
		}
	}
	
	if(local_change.size() == 0) emsg("No change");
	
	sort(local_change.begin(),local_change.end(),LocalChange_ord);    
	
	vector <unsigned int> cmap_list;
	
	auto cmap = comp_local_map[p];

	auto &pdmap = pop_data_map[p];
	vector <unsigned int> pd_ref_list;
	
	const auto &obs_eqn_value = ssp.obs_eqn_value;
	const auto &pop_filter = sp.pop_filter;
	const auto &pop_data = sp.pop_data;
	const auto &pop_data_ref = sp.pop_data_ref;

	// Updates compartments and populations
	for(auto j = 0u; j < local_change.size(); j++){
		const auto &lc = local_change[j];
		
		int sign = 0;
		
		switch(lc.type){
		case ADD_EVENT: sign = 1; break;
		case REMOVE_EVENT: sign = -1; break;
		case ADD_IC: sign = -1; break;
		case REMOVE_IC: sign = 1; break;
		case ADD_C_IC: sign = 1; break;
		case REMOVE_C_IC: sign = -1; break;
		}
		
		auto c_old = UNSET, c_new = UNSET;
		
		switch(lc.type){
		case ADD_EVENT: case REMOVE_EVENT: case ADD_IC: case REMOVE_IC: 
			{
				const auto &tr = sp.tra_gl[lc.tr];
				c_old = tr.i;	c_new = tr.f;
			}
			break;
			
		case ADD_C_IC: case REMOVE_C_IC: c_new = lc.tr; break;
		}
	
		if(c_old != UNSET){
			if(cmap[c_old] == UNSET){ cmap_list.push_back(c_old); cmap[c_old] = -sign;}
			else cmap[c_old] -= sign;
			
			for(auto k : sp.comp_gl[c_old].pop_ref_simp){
				if(pop_map[k] == UNSET){ pop_list.push_back(k); pop_map[k] = -sign;}
				else pop_map[k] -= sign;
			}
		}
	
		if(c_new != UNSET){
			if(cmap[c_new] == UNSET){ cmap_list.push_back(c_new); cmap[c_new] = sign;}
			else cmap[c_new] += sign; 
	
			for(auto k : sp.comp_gl[c_new].pop_ref_simp){
				if(pop_map[k] == UNSET){ pop_list.push_back(k); pop_map[k] = sign;}
				else pop_map[k] += sign;
			}
		}
		
		switch(lc.type){  // Makes change to initial condition
		case ADD_IC: case REMOVE_IC: case ADD_C_IC: case REMOVE_C_IC:
			if(sign	== 1) ssp.likelihood_init_cond_change(c_old,c_new,like_ch);
			else ssp.likelihood_init_cond_change(c_new,c_old,like_ch);
			break;
			
		default: break;
		}
		
		auto ti = (unsigned int)(lc.ti+1);
		auto ti_next = T; if(j+1 < local_change.size()) ti_next = local_change[j+1].ti+1;
		
		// Updates cpop_st
		for(auto c : cmap_list){
			auto ma = cmap[c];
			if(ma != 0){
				ssp.cpop_st_update(ti,ti_next,c,ma);
				for(auto tii = ti; tii < ti_next; tii++){
					double val;
					for(auto ref : pop_data_ref[tii][c]){
						const auto &pd = pop_data[ref];			
						
						if(pd.time_vari == true) val = ma*obs_eqn_value[pd.comp_obs_mod_ref[c]];
						else val = ma*obs_eqn_value[pop_filter[pd.ref].comp_obs_mod_ref[c]];
					
						if(pdmap[ref] == UNSET){ pd_ref_list.push_back(ref); pdmap[ref] = val;}
						else pdmap[ref] += val;
					}
				}	
				
				like_ch.markov += ssp.Li_update_c(c,ma,ti,ti_next,popnum_t);
			}
		}
		
		// Updates populations
		auto fl = false;
		
		vector <PopChange> pop_change;
		for(auto k : pop_list){
			auto ma = pop_map[k];
			//const auto &po = model.pop[k];
			if(ma != 0){
				back_pop.push_back(BackPop(POP_NUM_T,ti,ti_next,k,ma));

				PopChange po_ch; po_ch.po = k; po_ch.num = ma; 
				pop_change.push_back(po_ch);

				change_population(ti,ti_next,k,ma);
				fl = true;
			}
		}		
		
		if(fl == true) update_pop_change(ti,ti_next,pop_change,like_ch.markov);
	}
	
	for(auto k : pop_list) pop_map[k] = UNSET;
	pop_list.clear();
	
	for(auto c : cmap_list) cmap[c] = UNSET;

	auto &pmap = pop_trans_map[p];
	vector <unsigned int> ref_list;
	
	const auto &pop_trans_ref = sp.pop_trans_ref;
	const auto &pop_trans_filter = sp.pop_trans_filter;
	const auto &pop_trans_data = sp.pop_trans_data;
	auto &pop_trans_data_tgl = ssp.pop_trans_data_tgl;
	
	// Updates transition number along with poptrans observations
	for(auto j = 0u; j < local_change.size(); j++){
		const auto &lc = local_change[j];
		auto ti = lc.ti; 
		if(ti >= 0){
			auto tr = lc.tr;
			
			int sign = 0;
			
			switch(lc.type){
			case ADD_EVENT: sign = 1; break;
			case REMOVE_EVENT: sign = -1; break;
			case ADD_IC: case ADD_C_IC: case REMOVE_IC: case REMOVE_C_IC: emsg("Should not be here"); break;
			}
			like_ch.markov += ssp.Li_update_tn(tr,ti,sign);
			
			double val;
			ssp.back.push_back(Back(POP_TRANS_DATA_TGL,ti,tr,sign));
			
			for(auto ref : pop_trans_ref[ti][tr]){
				const auto &ptd = pop_trans_data[ref];
			
				pop_trans_data_tgl[ref][tr] += sign;
					
				if(ptd.time_vari == true){
					val = sign*obs_eqn_value[ptd.trans_obs_mod_ref[tr]];
				}
				else{
					val = sign*obs_eqn_value[pop_trans_filter[ptd.ref].trans_obs_mod_ref[tr]];
				}
				
				if(pmap[ref] == UNSET){ ref_list.push_back(ref); pmap[ref] = val;}
				else pmap[ref] += val;
			}
		}
	}
	
	// Updates event times
	
	// Updates pop_data measurements
	for(auto ref : pd_ref_list){
		auto dnum = pdmap[ref];
		if(dnum != 0) like_ch.obs += ssp.update_single_pop_data(dnum,ref);
		pdmap[ref] = UNSET;
	}
	
	// Updates pop_trans_measurements
	for(auto ref : ref_list){
		auto dnum = pmap[ref];
		if(dnum != 0) like_ch.obs += ssp.update_single_pop_trans_data(dnum,ref);
		pmap[ref] = UNSET;
	}
	
	return like_ch;
}

// This deals with general state species functions and population-based updates

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

/// Initialises the state species
StateSpecies::StateSpecies(const vector <double> &param_val, const vector <SplineValue> &spline_val, const vector <Equation> &eqn, const vector <Param> &param, const vector <ParamVecEle> &param_vec, const vector <Population> &pop, const Species &sp, const GeneticData &genetic_data, const Details &details, const	vector <double> &timepoint, const	vector <double> &dtimepoint, const vector <unsigned int> &pop_affect_, Operation mode_, const double &dif_thresh) : source_sampler(sp.markov_eqn,sp.tra_gl,sp.comp_gl,timepoint,dtimepoint,details,sp.init_cond), rate_mean(details), param_val(param_val), spline_val(spline_val), eqn(eqn), param(param), param_vec(param_vec), pop(pop), sp(sp), genetic_data(genetic_data), details(details), timepoint(timepoint), dtimepoint(dtimepoint), dif_thresh(dif_thresh)
{
	timer.resize(STSP_TIMER_MAX,0);
	
	pop_affect = pop_affect_;
	
	N = sp.markov_eqn.size();
	T = timepoint.size()-1;
	
	mode = mode_;
	
	nnode = sp.markov_tree.node.size();
	
	type = sp.type;
	
	initialise_arrays();
	reset_arrays();
}


/// Initialises all arrays
void StateSpecies::initialise_arrays()
{	
	auto G = sp.tra_gl.size();

	if(sp.type == POPULATION){
		trans_num.resize(G);
		tnum_mean_st.resize(G);
		trans_num_f.resize(G);
		tnum_mean_st_f.resize(G);
		for(auto tr = 0u; tr < G; tr++){
			trans_num[tr].resize(T);
			tnum_mean_st[tr].resize(T);
			trans_num_f[tr].resize(T);
			tnum_mean_st_f[tr].resize(T);
		}
	}		
	
	if(sp.type == INDIVIDUAL){
		nm_trans_ev_ref.resize(sp.nm_trans.size());
		for(auto m = 0u; m < sp.nm_trans.size(); m++){
			nm_trans_ev_ref[m].resize(T);
		}
		
		nm_trans_incomp_ref.resize(sp.nm_trans_incomp.size());
		for(auto n = 0u; n < sp.nm_trans_incomp.size(); n++){
			nm_trans_incomp_ref[n].resize(T);
		}
	}
	
	Li_ie.resize(sp.ind_eff_group.size());
	
	Li_markov_pop.resize(G);
	
	pop_data_num.resize(sp.pop_data.size());
	pop_trans_data_num.resize(sp.pop_trans_data.size());

	if(type == INDIVIDUAL){
		pop_data_cgl.resize(sp.pop_data.size());
		for(auto i = 0u; i < sp.pop_data.size(); i++){
			pop_data_cgl[i].resize(sp.comp_gl.size());
		}
	}
	
	pop_trans_data_tgl.resize(sp.pop_trans_data.size());
	for(auto i = 0u; i < sp.pop_trans_data.size(); i++){
		pop_trans_data_tgl[i].resize(sp.tra_gl.size());
	}
	
	obs_eqn_value.resize(sp.obs_eqn.size());
	
	auto E = sp.obs_trans_eqn.size();
	obs_trans_eqn_value.resize(E);
	obs_trans_eqn_num.resize(E);
	auto all_time = seq_vec(T);
	for(auto i = 0u; i < E; i++){
		obs_trans_eqn_value[i].resize(T,0);
		obs_trans_eqn_num[i].resize(T,0);
	}
 
	rate_store_init();
}


/// Resets arrays
void StateSpecies::reset_arrays()
{
	auto G = sp.tra_gl.size();
	trans_num.clear(); trans_num.resize(G);
	tnum_mean_st.clear(); tnum_mean_st.resize(G);
	ind_sim_c.clear();
	cpop_st.clear();
	
	for(auto m = 0u; m < sp.nm_trans.size(); m++){
		for(auto ti = 0u; ti < T; ti++) nm_trans_ev_ref[m][ti].clear();
	}
	
	for(auto n = 0u; n < sp.nm_trans_incomp.size(); n++){
		for(auto ti = 0u; ti < T; ti++) nm_trans_incomp_ref[n][ti].clear();
	}
	
	markov_eqn_vari.clear();
	markov_eqn_vari.resize(N);
	for(auto e = 0u; e < N; e++){
		auto &me_vari = markov_eqn_vari[e];
	
		me_vari.time_vari = sp.markov_eqn[e].time_vari;
		
		if(me_vari.time_vari == true){
			me_vari.dt = details.dt; 
			me_vari.div.resize(T);
		}		
		else{
			me_vari.dt = UNSET;
			me_vari.div.resize(1);
		}			
		
		me_vari.indfac_sum = UNSET;
		me_vari.value = UNSET;
	}
	
	trig_div.clear();
	trig_div.resize(T);
		
	type = sp.type;
	
	individual.clear();
	
	auto C = sp.comp_gl.size();
	
	cpop.clear();
	cpop.resize(C,0);	
	
	ti_sort = UNSET;
}


/// Initialises species ready for simulation be be performed
void StateSpecies::simulate_init()
{
	reset_arrays();
	
	ie_sampler_init();                            // Samples individual effects
}


/// Initialises individuals to be simulated
void StateSpecies::simulate_individual_init()
{
	// Loads observation equations (using in sampling individual initial state
	calculate_obs_eqn(seq_vec(obs_eqn_value.size()));
	
	auto C = sp.comp_gl.size();
	
	vector <unsigned int> pop_added(C,0);
	
	vector <unsigned int> unallocated;
	
	/// Adds any individuals from the data
	for(auto i = 0u; i < sp.nindividual_obs; i++){
		const auto &ind = sp.individual[i];
		
		auto added_flag = false;
		
		if(i < sp.nindividual_in){		
			for(const auto &ev : ind.ev){
				auto t = ev.t;
				
				if(t > details.t_start){
					SimTrigEvent trig; trig.i = i; trig.t = t;
					
					switch(ev.type){
					case ENTER_EV:
						trig.type = ENTER_SIM_EV;
						added_flag = true;
						trig.c = ind_sample_init_c(ind); trig.trg = UNSET;
						break;
						
					case LEAVE_EV:
						trig.type = LEAVE_SIM_EV;
						trig.c = UNSET; trig.trg = UNSET;
						break;
					
					case MOVE_EV:
						trig.type = MOVE_SIM_EV;
						trig.c = ev.move_c; trig.trg = ev.cl;
						break;
						
					default: emsg("event type not recog"); break;
					}	
					
					insert_trigger_event(trig);
				}
			}
			
			// Adds in event data (used during inference to set up initial state)
			for(auto m = 0u; m < ind.obs.size(); m++){
				const auto &ob = ind.obs[m];

				switch(ob.type){
				case OBS_TRANS_EV: case OBS_SOURCE_EV: case OBS_SINK_EV:
					{
						SimTrigEvent trig; 
						trig.type = DATA_TRANS_SIM_EV;
						trig.i = i; trig.t = ob.t; trig.c = m; trig.trg = UNSET;
						insert_trigger_event(trig);
						
						// Determines if individual is most likely to enter through an observed source
						if(added_flag == false){
							for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
								if(sp.tra_gl[tr].variety == SOURCE_TRANS){
									auto pr = like_trans_prob_obs(tr,ob);
									if(pr != LI_WRONG){ added_flag = true; break;}
								}
							}
						}						
					}
					break;

				default: break;
				}
			}
		}
		else{
			added_flag = true;
		}
		
		auto ii = add_individual(OBSERVED_IND,ind.name);
		
		// Directly adds ENTER_EV at start of simulation
		if(added_flag == false){                   // If individual not added then add at beginning
			unsigned int c;
			if(ind.ev.size() > 0 && ind.ev[0].type == ENTER_EV){
				if(ind.ev[0].t != details.t_start) emsg("Should be at start");
				c = ind_sample_init_c(ind);
				
				pop_added[c]++;
			
				IndInfFrom iif; if(sp.comp_gl[c].infected == true) iif.p = ENTER_INF;
			
				add_event(ENTER_EV,ii,UNSET,UNSET,UNSET,c,details.t_start,iif);
			}
			else{ 
				unallocated.push_back(ii);
			}
		}
	}

	
	ie_Amatrix_sampler_init();
	sample_ie_Amatrix();	
	
	if(sp.init_cond.type == INIT_POP_NONE){
		init_cond_val.cnum = pop_added;					
	}
		
	switch(type){
	case POPULATION:
		{
			for(auto c = 0u; c < C; c++){
				cpop[c] = double(init_cond_val.cnum[c]);
			}
		}			
		break;

	case INDIVIDUAL:
		{
			if(sp.init_cond.type != INIT_POP_NONE){
				string pre = ""; if(sp.nindividual_obs > 0) pre = "UO ";
				auto n = 0;
				for(auto c = 0u; c < sp.comp_gl.size(); c++){
					if(init_cond_val.cnum[c] > 0 && sp.fix_effect.size() > 0){
						string st = "";
						for(const auto &fe : sp.fix_effect){
							if(st != "") st += ","; 
							st += fe.name;
						}
						emsg("Cannot add unspecified individuals for a model with fixed effect(s) '"+st+"'"); 
					}
					
					auto imax = init_cond_val.cnum[c];	
					for(auto i = 0u; i < imax; i++){	
						unsigned int ii;
						if(unallocated.size() > 0){
							ii = unallocated[0];
							unallocated.erase(unallocated.begin());
						}
						else ii = add_individual(UNOBSERVED_IND);
						
						IndInfFrom iif; if(sp.comp_gl[c].infected == true) iif.p = ENTER_INF;
				
						add_event(ENTER_EV,ii,UNSET,UNSET,UNSET,c,details.t_start,iif);
						
						n++;
					}
				}
			}
		}
		break;
	}
	
	// Works out situation with individuals which do not have a specified starting position
	if(unallocated.size() > 0){
		vector <SourceSamp> source_samp;           // Generates a sampler for the source
		auto sum = 0.0;
		for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
			const auto &tra = sp.tra_gl[tr];
		
			if(tra.i == UNSET){
				const auto &eq = eqn[tra.dist_param[0].eq_ref];
				auto value = 1.0;
				if(eq.pop_ref.size() == 0) value = eq.calculate_no_popnum(0,param_val,spline_val);
				
				if(value != 0){
					SourceSamp ss; ss.tr_gl = tr; ss.prob_sum = sum;
					source_samp.push_back(ss);
				}
			}
		}
		
		switch(sp.init_cond.type){
		case INIT_POP_NONE: case INIT_POP_DIST:
			for(auto i : unallocated){
				const auto &ind = sp.individual[i];
				
				auto add_so = false;
				
				auto trange = source_time_range(ind);
				if(trange.tmin != details.t_start) add_so = true;
				else{
					if(source_samp.size() > 0){
						if(ran() < 0.5) add_so = true;
					}
				}
				
				if(add_so){ // Adds individual as a source
					add_ind_source(i,source_samp);
				}
				else{	      // Adds individual at start
			
					auto c = ind_sample_init_c(ind);
				
					pop_added[c]++;
			
					IndInfFrom iif; if(sp.comp_gl[c].infected == true) iif.p = ENTER_INF;
					add_event(ENTER_EV,i,UNSET,UNSET,UNSET,c,details.t_start,iif);
				}
			}
			break;
			
		case INIT_POP_FIXED:
			if(!sp.contains_source) emsg("Either an initial population distribution is added or the model must contain a source");
		
			for(auto i : unallocated){
				add_ind_source(i,source_samp);
			}
			break;
		}
	}
	
	for(auto f = 0u; f < sp.fix_effect.size(); f++) set_exp_fe(f);
	
	prob_trans_tree = 0;
	
	set_all_obs_trans_value();
}


/// Adds an individual through a source transition
void StateSpecies::add_ind_source(unsigned int i, const vector <SourceSamp> &source_samp)
{
	const auto &ind = sp.individual[i];
	auto trange = source_time_range(ind);
	
	if(source_samp.size() == 0) emsg("Cannot sample from source sampler");
	
	auto sum = source_samp[source_samp.size()-1].prob_sum;
	
	auto z = ran()*sum;
	auto j = 0u; 
	while(j < source_samp.size() && z > source_samp[j].prob_sum) j++;
	if(j == source_samp.size()) emsg("Could not simulate source");
	
	SimTrigEvent trig; 
	trig.type = SOURCE_SIM_EV;
	trig.i = i; trig.t = trange.tmin+ran()*(trange.tmax-trange.tmin); trig.c = UNSET; 
	trig.trg = source_samp[j].tr_gl;
	insert_trigger_event(trig);
}
	

/// Gets the potential time range to add an individual through a source
TRange StateSpecies::source_time_range(const IndData &ind) const
{
	const auto &obs = ind.obs;
	auto tmin = details.t_start;
	auto tmax = details.t_end;
	
	auto j = 0u;
	while(j < obs.size() && obs[j].not_alive == true) j++;
	
	if(j < obs.size()){
		auto t = obs[j].t;
		if(t < tmax) tmax = t;
		if(j > 0){
			auto tmi = obs[j-1].t;
			if(tmi > tmin) tmin = tmi;
			if(t < tmax) tmax = t;
		}
	}
	
	TRange tr; tr.tmin = tmin; tr.tmax = tmax;
	
	return tr;
}


/// Warning message when loading posterior sample
void StateSpecies::error_load_sample(unsigned int num) const
{
	emsg("Problem loading posterior sample "+to_string(num));
}


/// Initialises based on a posterior sample
void StateSpecies::simulate_sample_init(unsigned int ti_end, const SampleSpecies &samp_sp)
{
	auto t_end = timepoint[ti_end];
	
	string warn = "Problem loading posterior sample";
			
	switch(sp.type){
	case INDIVIDUAL:
		{
			Hash hash_compgl;
			for(auto c = 0u; c < sp.comp_gl.size(); c++){
				//auto hash_vec = hash_compgl.get_vec_string(sp.comp_gl[c].name);
				//hash_compgl.add(c,hash_vec);
				hash_compgl.add(c,sp.comp_gl[c].name);
			}
			
			Hash hash_tra;
			vector <TransRef> trans_ref;
			for(auto cl = 0u; cl < sp.ncla; cl++){
				const auto &claa = sp.cla[cl];
				for(auto tr = 0u; tr < claa.ntra; tr++){
					TransRef tref; tref.cl = cl; tref.tr = tr;
					//auto hash_vec = hash_tra.get_vec_string(claa.tra[tr].name);
					//hash_tra.add(trans_ref.size(),hash_vec);
					hash_tra.add(trans_ref.size(),claa.tra[tr].name);
					trans_ref.push_back(tref);
				}
			}
		
			Hash hash_ind; 
			for(auto i = 0u; i < sp.nindividual_obs; i++){
				const auto &ind = sp.individual[i];
				
				auto ii = add_individual(OBSERVED_IND,ind.name);
				//auto hash_vec = hash_ind.get_vec_string(ind.name);
				//hash_ind.add(ii,hash_vec);
				hash_ind.add(ii,ind.name);
			}
			
			const auto &ind_tab = samp_sp.ind_tab;
			
			// Checks individual effects are correct
			if(ind_tab.ncol-3 != sp.ind_effect.size()){
				error_load_sample(1);
			}
			
			if(ind_tab.ncol-3 != sp.ind_effect.size()){
				error_load_sample(20);
			}
			
			for(auto e = 0u; e < sp.ind_effect.size(); e++){
				if(ind_tab.heading[e+2] != sp.ind_effect[e].name){
					error_load_sample(2);
				}
			}
			
			for(auto r = 0u; r < ind_tab.nrow; r++){
				const auto &row = ind_tab.ele[r];
				
				auto name = row[0];
				auto hash_vec = hash_compgl.get_vec_string(name);
				auto i = hash_ind.existing(hash_vec);
				if(i == UNSET){
					i = add_individual(OBSERVED_IND,name);
					hash_ind.add(i,hash_vec);
				}					
			
				auto &ind = individual[i];
				
				auto source = false;
				if(row[1] == "yes") source = true;
				
				if(ind_tab.ncol != ind.ie.size()+3) error_load_sample(21);
				else{
					for(const auto &ieg : sp.ind_eff_group){
						if(ieg.ppc_resample == false){
							for(const auto &li : ieg.list){
								auto e = li.index;
								ind.ie[e] = log(number(row[2+e]));
							}
						}
					}
				}
				
				set_exp_ie(ind);
				
				auto ev_text = row[row.size()-1];
				if(ev_text != "unobserved"){
					auto col = split(ev_text,' ');		
				
					auto c = UNSET;
					for(auto e = 0u; e < col.size(); e++){
						auto tspl = split(col[e],':');
						
						if(tspl.size() != 2) error_load_sample(3);
						auto t = number(tspl[1]); if(t == UNSET) error_load_sample(4);
						
						if(t <= t_end){
							if(c == UNSET && source == false){
								//auto hash_vec = hash_compgl.get_vec_string(tspl[0]);
								//c = hash_compgl.existing(hash_vec);
								c = hash_compgl.find(tspl[0]);
								if(c == UNSET) error_load_sample(5);
							
								IndInfFrom inf_from;
								add_event(ENTER_EV,i,UNSET,UNSET,UNSET,c,t,inf_from);
							}
							else{
								auto tr = tspl[0];
						
								tr = replace(tr,"->","→");
							
								//auto hash_vec = hash_tra.get_vec_string(tr);
								//auto trr = hash_tra.existing(hash_vec);
								auto trr = hash_tra.find(tr);
								if(trr == UNSET) error_load_sample(8);
								
								const auto &tref = trans_ref[trr];
								
								auto trg = sp.trg_from_tr[c][tref.cl][tref.tr];
								if(trg == UNSET) error_load_sample(9);
							
								const auto &trgl = sp.tra_gl[trg];
								if(trgl.i != c) error_load_sample(10);
								
								auto ty = M_TRANS_EV;
								if(trgl.nm_trans_ref != UNSET) ty = NM_TRANS_EV;
								c = trgl.f;
								
								IndInfFrom inf_from;

								add_event(ty,i,trg,UNSET,UNSET,c,t,inf_from);				
							}
						}
					}
				}		
			}	
		
			// Looks to add any move leave events 
			// Loads observation equations (using in sampling individual initial state
			calculate_obs_eqn(seq_vec(obs_eqn_value.size()));
			
			// TO DO added individuals under PPC

			/// Adds any individuals from the data
			for(auto i = 0u; i < sp.nindividual_obs; i++){
				const auto &ind = sp.individual[i];
				
				//auto added_flag = false;
				
				if(1 == 1 && i < sp.nindividual_in){		
					for(const auto &ev : ind.ev){
						auto t = ev.t;
						
						if(t > t_end){
							SimTrigEvent trig; trig.i = i; trig.t = t;
							
							switch(ev.type){
							case ENTER_EV:
								trig.type = ENTER_SIM_EV;
								//added_flag = true;
								trig.c = ind_sample_init_c(ind); trig.trg = UNSET;
								break;
								
							case LEAVE_EV:
								trig.type = LEAVE_SIM_EV;
								trig.c = UNSET; trig.trg = UNSET;
								break;
							
							case MOVE_EV:
								trig.type = MOVE_SIM_EV;
								trig.c = ev.move_c; trig.trg = ev.cl;
								break;
								
							default: emsg("event type not recog"); break;
							}	
							
							insert_trigger_event(trig);
						}
					}
					
					/*
					/// Adds in event data (used during inference to set up initial state)
					for(auto m = 0u; m < ind.obs.size(); m++){
						const auto &ob = ind.obs[m];
						if(ob.type == OBS_TRANS_EV){
							SimTrigEvent trig; 
							trig.type = DATA_TRANS_SIM_EV;
							trig.i = i; trig.t = ob.t; trig.c = m; trig.trg = UNSET;
							insert_trigger_event(trig);
							
							/// Determines if individual is most likely to enter through an observed source
							if(added_flag == false){
								for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
									if(sp.tra_gl[tr].variety == SOURCE_TRANS){
										auto pr = like_trans_prob_obs(ob.obs_eqn_ref[tr]);
										if(pr != LI_WRONG){ added_flag = true; break;}
									}
								}
							}						
						}
					}
					*/
				}
				else{
					//added_flag = true;
				}
				
				
					
				//auto ii = add_individual(OBSERVED_IND,ind.name);
				//auto hash_vec = hash_ind.get_vec_string(row[0]);
				//hash_tra.add(ii,hash_vec);
				//auto ii = i;
				
				// TO DO
				
				/*
				auto ii = add_individual(OBSERVED_IND,ind.name);
			
				*/
		
				
				/*
				// Directly adds ENTER_EV at start of simulation
				// If individual not added then add at beginning
				if(added_flag == false && individual[ii].ev.size() == 0){
					unsigned int c;
					if(ind.ev.size() > 0 && ind.ev[0].type == ENTER_EV){
						if(ind.ev[0].t < t_end){
							stringstream ss;
							ss << "Cannot add individual '" << ind.name << "' at time " << ind.ev[0].t << " prior to posterior simulation start time " << to_string(t_end) << ".";
							emsg(ss.str());
						}
						
						if(ind.ev[0].t != t_end) emsg("Should be at start");
						c = ind_sample_init_c(ind);
						
						IndInfFrom iif; if(sp.comp_gl[c].infected == true) iif.p = ENTER_INF;
					
						add_event(ENTER_EV,ii,UNSET,UNSET,UNSET,c,details.t_start,iif);
					}
				}
				*/
			}
		
			
			for(auto f = 0u; f < sp.fix_effect.size(); f++) set_exp_fe(f);
		}
		
		ie_Amatrix_sampler_init();
		sample_ie_Amatrix();	
		//emsg("do");
		break;
	
	case POPULATION:
		{
			const auto &tab = samp_sp.cpop_init_tab;
			auto C = sp.comp_gl.size();
			if(tab.nrow != C) error_load_sample(11);

			for(auto c = 0u; c < C; c++){
				if(tab.ele[c][0] != sp.comp_gl[c].name) error_load_sample(12);
				cpop[c] = number(tab.ele[c][1]); 
			}
			
			init_cond_val.cnum.resize(C);
			for(auto c = 0u; c < C; c++) init_cond_val.cnum[c] = (unsigned int)(cpop[c]);
		}
		
		{
			const auto &tab = samp_sp.trans_num_tab;
		
			for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
				if(tab.ele[tr][0] != sp.tra_gl[tr].name) error_load_sample(15);
				
				auto &tnum = trans_num[tr];
				
				auto spl = split(tab.ele[tr][1],',');
				auto k = 0u;
				while(k < spl.size() && tnum.size() < ti_end){
					auto val = spl[k];
					if(val.substr(0,1) == "Z"){
						auto num = number(val.substr(1));
						if(num == UNSET || num != (unsigned int)(num)) error_load_sample(15);
						for(auto j = 0u; j < num; j++){
							tnum.push_back(0);
							if(tnum.size() == ti_end) break;
						}
					}
					else{
						auto num = number(val);
						tnum.push_back(num);
					}
					k++;
				}
			}
			
			for(auto ti = 0u; ti < ti_end; ti++){
				vector <double> tnum(sp.tra_gl.size());
				for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
					tnum[tr] = trans_num[tr][ti];
				}
			
				cpop_st.push_back(cpop);
	
				update_cpop(ti,cpop,tnum);
			}
			
			cpop_st.push_back(cpop);
		}
		break;
	}
}


/// Sets values for tnum_mean
void StateSpecies::set_tnum_mean(unsigned int ti_end, const vector < vector <double> > &popnum_t)
{
	for(auto ti = 0u; ti < ti_end; ti++){
		auto ddt = timepoint[ti+1] - timepoint[ti];
	
		auto tnum_mean = calculate_tnum_mean(ti,popnum_t[ti],cpop_st[ti],param_val,ddt);
		
		for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
			tnum_mean_st[tr].push_back(tnum_mean[tr]);
		}
	}
}


/// Samples the initial state for the individual
unsigned int StateSpecies::ind_sample_init_c(const IndData &ind) const
{
	const auto &ent = sp.enter[ind.enter_ref];
	
	if(ent.c_set != UNSET) return ent.c_set;
	
	auto c_after = 0u;
	for(auto cl = 0u; cl < sp.ncla; cl++){
		const auto &ent_cl = ent.cla[cl];
		if(ent_cl.c_set != UNSET) c_after += sp.comp_mult[cl]*ent_cl.c_set;
		else{
			auto sum = 0u; 
			vector <double> prob, prob_sum;
			const auto &claa = sp.cla[cl];
			const auto &oer = ent_cl.obs_eqn_ref;
			auto C = claa.ncomp;
			for(auto c = 0u; c < C; c++){
				auto pr = obs_eqn_value[oer[c]];
				sum += pr;
				prob.push_back(pr); prob_sum.push_back(sum);
			}
			
			auto z = ran()*sum;
			auto c = 0u; while(c < C && z > prob_sum[c]) c++;
			if(c == C) emsg("Zero problem");
			
			c_after += sp.comp_mult[cl]*c;
		}
	}
	
	return c_after;
}


/// Gets a new unique name for an individual
string StateSpecies::get_new_ind_name(string pre, unsigned int num) const
{
	auto loop = 0u;
	do{
		string name;
		if(loop == 0) name = pre+"Ind."+tstr(num); 
		else name = pre+"Ind."+tstr(num)+"("+alphabet[loop-1]+")";
		
		auto i = 0u; while(i < individual.size() && individual[i].name != name) i++;
		if(i == individual.size()) return name;
		loop++;
	}while(loop < 27);
	
	emsg("Could not get individual name");
	return "";
}


/// Sets up the exponential of the fixed effects
vector <double> StateSpecies::set_exp_fe(unsigned int f)
{	
	const auto &fe = sp.fix_effect[f];
		
	const auto &par = param[fe.th];

	double pval;
	if(par.variety == CONST_PARAM){
		pval = par.get_value(0);
	}
	else{
		auto th2 = par.get_param_vec(0);
		if(th2 == UNSET) emsg("Should not be unset1");		
		pval = param_val[th2];
	}
	
	vector <double> store;

	auto sum = 0.0;
	auto num = 0u;
	for(auto &ind : individual){
		store.push_back(ind.exp_fe[f]);
		
		auto X = ind.X[f];
		if(X != UNSET){
			auto val = exp(X*pval);
			ind.exp_fe[f] = val;
			sum += val;
			num++;
		}
	}
	
	if(sum == 0) emsg("Fixed effect problem");
	
	auto fac = num/sum;

	for(auto &ind : individual){
		auto X = ind.X[f];
		if(X == UNSET) ind.exp_fe[f] = 1;
		else ind.exp_fe[f] *= fac;
	}
	
	return store;
}


/// Restores the values for exp_fe
void StateSpecies::set_exp_fe_restore(unsigned int f, const vector <double> &store)
{
	auto j = 0u;
	for(auto &ind : individual){
		ind.exp_fe[f] = store[j]; j++;
	}
}


/// Sets the values for exponentiated individual effects
vector <double> StateSpecies::set_exp_ie(Individual &ind) const
{
	const auto &indeff = sp.ind_effect;
	const auto &iegs = ind_eff_group_sampler;
	
	auto N = indeff.size();
	
	if(ind.ie.size() != N) emsg("Sizes do not agree");
	
	vector <double> store;
	
	ind.exp_ie.resize(N);
	for(auto i = 0u; i < indeff.size(); i++){
		const auto &ie = indeff[i];
		auto var = iegs[ie.index].omega[ie.num][ie.num];
		store.push_back(ind.exp_ie[i]);
		
		ind.exp_ie[i] = exp(ind.ie[i]-0.5*var);
		//check_ie_out_of_range(i,ind);
	}
	
	return store;
}


/// Recalcualted exp_ie for a given individual effect
vector <double> StateSpecies::recalculate_exp_ie(unsigned int ie)
{
	auto &ind_eff = sp.ind_effect[ie];
	auto g = ind_eff.index;
	auto j = ind_eff.num;
	const auto &iegs = ind_eff_group_sampler[g];
	auto var = iegs.omega[j][j];
		
	vector <double> store;
	for(auto &ind : individual){	
		store.push_back(ind.exp_ie[ie]);
		ind.exp_ie[ie] = exp(ind.ie[ie]-0.5*var);
	}
	
	return store;
}


/// Restores values for exp_ie
void StateSpecies::recalculate_exp_ie_restore(unsigned int ie, const vector <double> &store)
{
	auto j = 0u;
	for(auto &ind : individual){
		ind.exp_ie[ie] = store[j]; j++;
	}
}


/// Initialises sampler used to sample individual effects
void StateSpecies::ie_sampler_init()
{
	ind_eff_group_sampler.clear();
	ind_eff_group_sampler.resize(sp.ind_eff_group.size());
	
	for(auto i = 0u; i < sp.ind_eff_group.size(); i++){
		const auto &ieg = sp.ind_eff_group[i];
		
		auto &iegs = ind_eff_group_sampler[i];
		auto N = ieg.list.size();
		iegs.omega.resize(N); iegs.omega_Z.resize(N); iegs.omega_inv.resize(N);
		for(auto j = 0u; j < N; j++){
			iegs.omega[j].resize(N); iegs.omega_Z[j].resize(N); iegs.omega_inv[j].resize(N);
		}
		
		calculate_omega(i);
	}
}



/// Initialises sampler for individual effects which use A martrix
void StateSpecies::ie_Amatrix_sampler_init()
{
	for(auto ig = 0u; ig < sp.ind_eff_group.size(); ig++){
		const auto &ieg = sp.ind_eff_group[ig];
		auto &iegs = ind_eff_group_sampler[ig];
		
		const auto &Aieg = ieg.A_matrix;
		if(Aieg.set == true){
			auto I = individual.size();
			vector <unsigned int> map(I);
			for(auto i = 0u; i < I; i++){
				auto name = individual[i].name;
				auto ii = Aieg.hash_ind_list.find(name);
				if(ii == UNSET) emsg("Cannot find individual '"+name+"' in the A matrix");
				map[i] = ii;
			}
			
			auto &A = iegs.A;
			A.resize(I);
			for(auto j = 0u; j < I; j++){
				A[j].resize(I);
				for(auto i = 0u; i < I; i++){
					A[j][i] = Aieg.value[map[j]][map[i]];
				}
			}
			
			auto illegal = false;
			iegs.A_Z = calculate_cholesky(A,illegal);
			if(illegal) emsg("Cholesky decomposition of A matrix not possible");
			
			auto A_inv = invert_matrix(A);
			tidy(A_inv);
			//print_matrix("Ainv final",A_inv);
					
			iegs.A_inv = A_inv;
			iegs.A_inv_diag.resize(I);
			iegs.A_inv_nonzero.resize(I);
			for(auto j = 0u; j < I; j++){
				iegs.A_inv_diag[j] = A_inv[j][j];
				for(auto i = 0u; i < I; i++){
					if(A_inv[j][i] != 0){
						NZElement nze; nze.i = i; nze.value = A_inv[j][i];
						iegs.A_inv_nonzero[j].push_back(nze);
					}
				}
			}				
		}
	}
}
			

/// Calculates the omega matrix for individual effect
vector <double> StateSpecies::calculate_omega(unsigned int g)
{
	auto &iegs = ind_eff_group_sampler[g];
	auto N = iegs.omega.size();
	
	vector <double> store;
	for(auto j = 0u; j < N; j++){
		for(auto i = 0u; i < N; i++){	
			store.push_back(iegs.omega[j][i]);
			store.push_back(iegs.omega_Z[j][i]);
			store.push_back(iegs.omega_inv[j][i]);
		}
	}
	
	iegs.omega = sp.calculate_omega_basic(g,param_val,param);
	
	auto illegal = false;
	iegs.omega_Z = calculate_cholesky(iegs.omega,illegal);
	if(illegal) emsg("Cholesky should not be illegal");
	
	iegs.omega_inv = invert_matrix(iegs.omega);
	
	return store;
}


/// Restore the values in omega
void StateSpecies::calculate_omega_restore(unsigned int g, const vector <double> &store) 
{
	auto &iegs = ind_eff_group_sampler[g];
	const auto &ieg = sp.ind_eff_group[g];
	auto N = ieg.list.size();
	
	auto k = 0u;
	for(auto j = 0u; j < N; j++){
		for(auto i = 0u; i < N; i++){
			iegs.omega[j][i] = store[k]; k++;
			iegs.omega_Z[j][i] = store[k]; k++;
			iegs.omega_inv[j][i] = store[k]; k++;
		}
	}
}


/// Samples the individual effects 
vector <double> StateSpecies::sample_ie() const
{
	auto N = sp.ind_effect.size();               // The number of individual effects
	
	vector <double> vec(N,UNSET);
	
	auto flag = false;
	for(auto i = 0u; i < sp.ind_eff_group.size(); i++){
		const auto &ieg = sp.ind_eff_group[i];
		
		const auto &iegs = ind_eff_group_sampler[i];
		
		if(ieg.A_matrix.set == false){
			auto samp = sample_mvn(iegs.omega_Z);
		
			if(samp.size() != ieg.list.size()) emsg("Wrong size2");
			
			for(auto j = 0u; j < ieg.list.size(); j++){
				auto k = ieg.list[j].index;
				if(k >= N) emsg("Out of range5");
				if(vec[k] != UNSET) emsg("Should not be set");
				vec[k] = samp[j];
			}
		}
		else flag = true;
	}
	
	if(flag == false){
		for(auto i = 0u; i < N; i++){
			if(vec[i] == UNSET) emsg("Should not be unset3");
		}
	}
	
	return vec;
}


/// Samples the individual effects with A matrix
void StateSpecies::sample_ie_Amatrix()
{
	auto flag = false;
	string warn;
	for(auto i = 0u; i < sp.ind_eff_group.size(); i++){
		const auto &ieg = sp.ind_eff_group[i];
		
		const auto &iegs = ind_eff_group_sampler[i];
		
		if(ieg.A_matrix.set == true && !(mode == PPC && ieg.ppc_resample == false)){
			flag = true;
			
			auto N = ieg.list.size();
			auto I = individual.size();
			auto M = N*I;
			
			vector <double> vec(M);
			for(auto i = 0u; i < M; i++) vec[i] = normal_sample(0,1,warn);
			
			for(auto j = 0u; j < M; j++){
				auto j_A = (unsigned int)(j/N);
				auto j_O = j%N;
				
				auto sum = 0.0;
				for(auto i = 0u; i < M; i++){
					auto i_A = (unsigned int)(i/N);
					auto i_O = i%N;
				
					sum += iegs.A_Z[j_A][i_A]*iegs.omega_Z[j_O][i_O]*vec[i];
				}
				
				individual[j_A].ie[ieg.list[j_O].index] = sum;
			}				
		}
	}
	
	for(auto &ind : individual){		
		set_exp_ie(ind);
		for(auto i = 0u; i < sp.ind_effect.size(); i++){
			if(ind.ie[i] == UNSET) emsg("Should be set1");
		}

		if(flag == true) set_exp_ie(ind);
	}
}


/// Returns the individual-based factor which acts on Markov equation
double StateSpecies::get_indfac(const Individual &ind, const MarkovEqn &mar_eqn) const
{
	auto val = 1.0;
	for(auto ie : mar_eqn.ind_eff_mult) val *= ind.exp_ie[ie];
	for(auto fe : mar_eqn.fix_eff_mult) val *= ind.exp_fe[fe];
	
	if(!mar_eqn.rate) val = 1.0/val;
		
	if(false){
		cout << ind.name << "   ";
		cout << eqn[mar_eqn.eqn_ref].te << " ";
		for(auto ie : mar_eqn.ind_eff_mult) cout << "mult IE " << ie << endl;
		for(auto fe : mar_eqn.fix_eff_mult) cout << "mult FE " << fe << endl;
		cout << val << "val" << endl;
	}
	
	return val;
}


/// Updates population-based species
void StateSpecies::update_population_based(unsigned int ti, bool stoc, const vector <double> &popnum)
{
	auto t = timepoint[ti];
	auto ddt = timepoint[ti+1] - t;
		
	auto tnum_mean = calculate_tnum_mean(ti,popnum,cpop,param_val,ddt);
	
	if(false){
		for(auto val : cpop) cout << val << ","; 
		cout << "cpop" << endl;
		
		for(auto val : param_val) cout << val << ","; 
		cout << "paramval" << endl;
		
		for(auto val : popnum) cout << val << ","; 
		cout << "popnum" << endl;
		
		for(auto val : tnum_mean) cout << val << ","; 
		cout << " val" << endl;
	}
	
	auto tnum = sample_trans_num(tnum_mean,stoc);
	
	for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
		tnum_mean_st[tr].push_back(tnum_mean[tr]);
		trans_num[tr].push_back(tnum[tr]);
	}
	
	cpop_st.push_back(cpop);
	
	update_cpop(ti,cpop,tnum);
}


/// Updates the compartmental populations based on transition numbers
void StateSpecies::update_cpop(unsigned int ti, vector <double> &cpop, const vector <double> &tnum) const
{
	for(auto i = 0u; i < sp.tra_gl.size(); i++){
		const auto &tra = sp.tra_gl[i];
		auto num = tnum[i];
	
		if(tra.i != UNSET) cpop[tra.i] -= num;
		if(tra.f != UNSET) cpop[tra.f] += num;
	}
	
	if(sp.add_rem_pop_on){
		for(auto c : sp.add_rem_pop_change[ti]){
			cpop[c] += sp.add_rem_pop[ti][c];
		}
	}
}


/// Samples transition number given transition rates
vector <double> StateSpecies::sample_trans_num(const vector <double> &tnum_mean, bool stochastic) const
{
	auto N = tnum_mean.size();
	vector <double> tnum(N);
	
	string warn;
	for(auto tr = 0u; tr < N; tr++){
		if(stochastic == true){
			auto val = poisson_sample(tnum_mean[tr],warn);	
			if(val == UNSET) sp.sampling_error(tr,warn);
			tnum[tr] = val;
		}
		else tnum[tr] = tnum_mean[tr];
	}
	
	return tnum;
}


/// Calculates the rate for different transitions
vector <double> StateSpecies::calculate_tnum_mean(unsigned int ti, const vector <double> &popnum, const vector <double> &cpop, const vector <double> &param_val, double dt) const
{
	auto N = sp.tra_gl.size();
	vector <double> tnum_mean(N);
	
	for(auto i = 0u; i < N; i++){
		tnum_mean[i] = calculate_tnum_mean(ti,i,popnum,cpop,param_val,dt);
	}
	
	return tnum_mean;
}

/// Calculates the rate for different transitions
double StateSpecies::calculate_tnum_mean(unsigned int ti, unsigned int i, const vector <double> &popnum, const vector <double> &cpop, const vector <double> &param_val, double dt) const
{
	const auto &tr = sp.tra_gl[i];
	if(tr.ev_type != M_TRANS_EV) emsg("Should be exponential");
	
	const auto &eq = eqn[tr.dist_param[0].eq_ref];
	double rate;

	switch(tr.type){
	case EXP_RATE: 
		{
			rate = eq.calculate(ti,popnum,param_val,spline_val);
			if(rate < 0){	
				emsg("The transition rate for '"+tr.name+"' through equation '"+eq.te_raw+"' has become negative."+check_prior(eq));
			}
		}
		break;
		
	case EXP_MEAN: 
		{
			auto mean = eq.calculate(ti,popnum,param_val,spline_val);
			if(mean <= 0){	
				if(mean < 0){	
					emsg("The transition mean for '"+tr.name+"' through equation '"+eq.te_raw+"' has become negative."+check_prior(eq));
				}
				else{
					emsg("The transition mean for '"+tr.name+"' through equation '"+eq.te_raw+"' has become zero."+check_prior(eq));
				}
			}
			rate = 1/mean;
		}
		break;
		
	default: rate = 0; emsg("SHould not be here"); break;
	}
	
	auto ci = tr.i;
	double pop;
	
	if(ci == UNSET) pop = 1.0;
	else{ pop = cpop[ci]; if(pop < 0) pop = 0;}
	
	return pop*rate*dt;
}


/// Generates an A matrix (for testing)
void StateSpecies::generate_A()
{
	ofstream fout("Testing/M.csv");
	for(auto i = 0u; i < individual.size(); i++){
		if(i != 0) fout << ",";
		fout << "\"" << individual[i].name << "\""; 
	}		
	fout << endl;
	
	for(auto j = 0u; j < individual.size(); j++){
		for(auto i = 0u; i < individual.size(); i++){
			if(i != 0) fout << ",";
			if(i == j) fout << 1;
			else{
				if((unsigned int)(i/10) == (unsigned int)(j/10)) fout << 0.5;
				else fout << 0;
			}
		}
		fout << endl;
	}
}


/// Performs a MBP
void StateSpecies::mbp(double sim_prob, vector < vector <double> > &popnum_t)
{
	cpop_st_f.clear();
	
	auto C = sp.comp_gl.size();
	vector <double> cpop(C);
	for(auto c = 0u; c < C; c++) cpop[c] = double(init_cond_val.cnum[c]);	
		
	auto N = sp.tra_gl.size();
	
	string warn;
	
	vector <int> num(N);
	for(auto ti = 0u; ti < T; ti++){
		const auto &popnum = popnum_t[ti];
	
		cpop_st_f.push_back(cpop);
	
		auto ddt = dtimepoint[ti];
	
		for(auto tr = 0u; tr < N; tr++){
			auto tnum_mean_i = tnum_mean_st[tr][ti];
			auto tnum_mean_f = calculate_tnum_mean(ti,tr,popnum,cpop,param_val,ddt);

			auto num_i = trans_num[tr][ti];

			unsigned int num_f;	
			if(sim_prob == 0 || ran() > sim_prob){  // Does a MBP
				if(tnum_mean_f > tnum_mean_i){
					auto val = poisson_sample(tnum_mean_f-tnum_mean_i,warn);
					if(val == UNSET) sp.sampling_error(tr,warn);
					num_f = num_i + val;
				}
				else{
					if(tnum_mean_i == tnum_mean_f) num_f = num_i;
					else{
						num_f = binomial_sample(tnum_mean_f/tnum_mean_i,num_i,warn);
						if(num_f == UNSET) sp.sampling_error(tr,warn);
					}
				}
			}
			else{                                   // Does a simulation
				num_f = poisson_sample(tnum_mean_f,warn);	
				if(num_f == UNSET) sp.sampling_error(tr,warn);
			}
		
			num[tr] = num_f;
			tnum_mean_st_f[tr][ti] = tnum_mean_f; 
			trans_num_f[tr][ti] = num_f;
		}
		
		/// Updates the populations
		if(ti < T-1){
			auto &popnum_new = popnum_t[ti+1];
			
			for(auto po : pop_affect) popnum_new[po] = popnum[po];
		
			for(auto tr = 0u; tr < N; tr++){
				auto n = num[tr];
				if(n != 0){
					const auto &tra = sp.tra_gl[tr];
					auto ci = tra.i;
					if(ci != UNSET){
						cpop[ci] -= n;
						for(auto po : sp.comp_gl[ci].pop_ref_simp) popnum_new[po] -= n;
					}
					
					auto cf = tra.f;
					if(cf != UNSET){
						cpop[cf] += n;
						for(auto po : sp.comp_gl[cf].pop_ref_simp) popnum_new[po] += n;
					}
				}
			}
			
			if(sp.add_rem_pop_on){
				for(auto c : sp.add_rem_pop_change[ti]){
					auto val = sp.add_rem_pop[ti][c];
					cpop[c] += val;
					for(auto po : sp.comp_gl[c].pop_ref_simp) popnum_new[po] += val;
				}
			}
		}
	}
	
	trans_num = trans_num_f;
	tnum_mean_st = tnum_mean_st_f;
	cpop_st = cpop_st_f;
}


/// Accepts a MBP
void StateSpecies::mbp_accept(double &like_ch, const vector < vector <double> > &popnum_t)
{
	like_ch -= sum(Li_markov_pop);
	for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
		auto &Li = Li_markov_pop[tr];
		auto &tn = trans_num[tr];
		auto &tnm = tnum_mean_st[tr];
		for(auto ti = 0u; ti < T; ti++){
			tnm[ti] = calculate_tnum_mean(ti,tr,popnum_t[ti],cpop_st[ti],param_val,dtimepoint[ti]);
			Li[ti] = poisson_probability(tn[ti],tnm[ti]);
		}
	}
	like_ch += sum(Li_markov_pop);
}


/// Moves a single event to a new time and generates the new event sequence
unsigned int StateSpecies::move_event(vector <Event> &ev, unsigned int index, double t_new) const 
{
	auto flag = false;
	auto ev_store = ev[index];
	auto E = ev.size();
	auto t = ev_store.t;
	
	if(t_new > t){
		while(index+1 < E && ev[index+1].t < t_new){
			ev[index] = ev[index+1];
			index++;
			flag = true;
		}							
	}
	else{
		while(index > 0 && ev[index-1].t > t_new){
			ev[index] = ev[index-1];
			index--;
			flag = true;
		}			
	}
	
	auto tr_gl_new = UNSET;

	// If the move event crosses an event in another classification
	if(flag == true){  
		ev[index] = ev_store;
		ev[index].t = t_new;
		
		ensure_consistent(ev);
		tr_gl_new = ev[index].tr_gl;
	}
	else{
		ev[index].t = t_new;
	}

	return tr_gl_new;
}


/// Gets ti div time from an actual time
unsigned int StateSpecies::get_ti(double t) const
{
	return (unsigned int)(OVER_ONE*(t-details.t_start)/details.dt);
}


/// Used to order events
bool Event_ord (const Event &ev1, const Event &ev2)                      
{ return (ev1.t < ev2.t); };  


/// Makes sure event sequence
void StateSpecies::make_consistent(vector <Event> &event) const
{
	if(sp.multi_trans_cl){	
		// Sorts events by time 
		sort(event.begin(),event.end(),Event_ord);    
		
		auto c = UNSET;
		for(auto k = 0u; k < event.size(); k++){	
			auto &ev = event[k];
			if(sp.correct_ev(c,ev) == true) emsg("Inconsistent");
			c = ev.c_after; 
		}
	}
	
	set_m_ti_origin(event);
}


/// Sets cpop_st for population-based models (based on initial condition and trans_num)
void StateSpecies::set_cpop_st()
{	
	if(sp.type != POPULATION) emsg("cpop_st only defined for population-based model");
		
	auto C = sp.comp_gl.size();
	vector <double> cpop(C);
	
	switch(sp.init_cond.type){
	case INIT_POP_FIXED: case INIT_POP_DIST:	
		for(auto c = 0u; c < C; c++) cpop[c] = double(init_cond_val.cnum[c]);
		break;
		
	case INIT_POP_NONE:
		for(auto c = 0u; c < C; c++) cpop[c] = 0;
		break;
	}
	
	auto N = sp.tra_gl.size();
			
	cpop_st.resize(T);
	for(auto ti = 0u; ti < T; ti++){
		auto C = sp.comp_gl.size();
			
		for(auto c = 0u; c < C; c++) cpop_st[ti] = cpop;
		
		for(auto tr = 0u; tr < N; tr++){
			auto n = trans_num[tr][ti];
		
			if(n != 0){
				const auto &tra = sp.tra_gl[tr];
				auto ci = tra.i;
				if(ci != UNSET) cpop[ci] -= n;
				
				auto cf = tra.f;
				if(cf != UNSET) cpop[cf] += n;
			}
		}
		
		if(sp.add_rem_pop_on){
			for(auto c : sp.add_rem_pop_change[ti]){
				cpop[c] += sp.add_rem_pop[ti][c];
			}
		}
	}
}


/// Gets the observation probability on an observation for a transition
double StateSpecies::get_trans_obs_prob(unsigned int trg, const ObsData &ob) const
{
	if(ob.time_vari) return obs_eqn_value[ob.obs_eqn_ref[trg]];
	else return obs_eqn_value[sp.obs_trans[ob.ref].obs_eqn_ref[trg]];
}


/// Compares a covariance matrix with one from individual effects
void StateSpecies::compare_covar(string te, const vector < vector <double> > &omega, unsigned int g) const
{
	const auto &ieg = sp.ind_eff_group[g];
	
	auto N = ieg.list.size();
	vector <double> av;
	vector < vector <double> > av2;
	av.resize(N); av2.resize(N);
	for(auto j = 0u; j < N; j++){	
		av2[j].resize(N,0);
	}
	
	for(const auto &ind : individual){
		for(auto i = 0u; i < N; i++){
			av[i] += ind.ie[ieg.list[i].index];
			for(auto j = 0u; j < N; j++){
				av2[i][j] += ind.ie[ieg.list[i].index]*ind.ie[ieg.list[j].index];
			}
		}
	}
	
	auto T = individual.size();
	
	vector < vector <double> > M;
	M.resize(N);
	for(auto j = 0u; j < N; j++){
		M[j].resize(N);
		for(auto i = 0u; i < N; i++){
			M[j][i] = av2[j][i]/T - (av[j]/T)*(av[i]/T);
		}
	}
	
	cout << te << endl;
	print_matrix("omega",omega);
	print_matrix("omega ind",M);
}


/// Gets a list of linked events
ListMove StateSpecies::get_ev_link_listmove(unsigned int e, const vector <Event> &event,  const Individual &ind, const vector < vector <double> > &popnum_t) const
{
	ListMove lm;

	auto E = event.size();
	
	// Keeps a map which stores all event moves considered so far
	auto &move = lm.move;
	move.resize(E,false);

	const auto &ev_ro = event[e];
	if(ev_ro.observed == true) return lm;

	if(ev_ro.type != NM_TRANS_EV && ev_ro.type != M_TRANS_EV) return lm;
		
	// Makes sure that first event is not linked to something earlier	
	if(ev_ro.type == NM_TRANS_EV){
		auto ee = ev_ro.e_origin;
		if(ev_link(ev_ro,ind,popnum_t)){
			if(!event[ee].observed) return lm;
			else{
				if(sp.tra_gl[ev_ro.tr_gl].type == PERIOD) return lm;
			}
		}
	}
	
	// Works out which transitions come after
	vector < vector <unsigned int> > after;
	after.resize(E);
	
	for(auto ee = 0u; ee < E; ee++){
		const auto &ev = event[ee];
		if(ev.observed == false && ev.type == NM_TRANS_EV){
			after[ev.e_origin].push_back(ee);
		}
	}
	
	// Makes a list of all events which are being moved
	auto &list = lm.list;

	list.push_back(e); move[e] = true;

	// Expands the number of events which are simultaneously moved
	auto i = 0u; 
	while(i < list.size()){
		auto ee = list[i];
		for(auto eee : after[ee]){
			if(!move[eee] && !event[eee].observed && ev_link(event[eee],ind,popnum_t)){
				list.push_back(eee);
				move[eee] = true;
			}
		}
		i++;
	}
	
	return lm;
}


/// Works out the time range over which the multi-event proposals can be moved 
TRange StateSpecies::ev_link_trange(const ListMove &lm, const vector <Event> &event) const 
{
	auto ncla = sp.ncla;
	auto E = event.size();
	auto t_start = details.t_start;
	auto t_end = details.t_end;

	const auto &move = lm.move;
	
	// Works out the size of any potential change which is valid
	vector <double> time(ncla,t_start);
	vector <unsigned int> e_st(ncla,UNSET);
	
	double dtmin = -LARGE, dtmax = LARGE;
	for(auto e = 0u; e < E; e++){
		const auto &ev = event[e];
		auto t = ev.t;
		
		if(ev.c_after == UNSET){
			for(auto cl = 0u; cl < ncla; cl++){
				auto ee = e_st[cl];
				if(ee != UNSET && move[ee]){
					auto dtt = t-time[cl];
					if(dtt < dtmax) dtmax = dtt;
				}
			}
		}
		else{
			switch(ev.type){ 
			case ENTER_EV:
				for(auto cl = 0u; cl < ncla; cl++){
					time[cl] = t;
					e_st[cl] = e;
				}
				break;
				
			case LEAVE_EV:
				emsg("Should be dealt with");
				break;
				
			case MOVE_EV: case NM_TRANS_EV: case M_TRANS_EV: 
				{
					auto cl = ev.cl;
					auto ee = e_st[cl];
					
					if(move[e] && !(ee != UNSET && move[ee])){
						auto dtt = time[cl]-t;
						if(dtt > dtmin) dtmin = dtt;
					}
					
					if(!move[e] && (ee != UNSET && move[ee])){
						auto dtt = t-time[cl];
						if(dtt < dtmax) dtmax = dtt;
					}
				
					if(move[e]){
						auto dtt = t_end-t;
						if(dtt < dtmax) dtmax = dtt;
					}
					
					time[cl] = ev.t;
					e_st[cl] = e;
				}
				break;
			}	
		}
	}
	
	TRange trange; trange.tmin = dtmin; trange.tmax = dtmax;
	return trange;	
}


/// Determines if two events are linked, either because of period or a small CV
bool StateSpecies::ev_link(const Event &ev, const Individual &ind, const vector < vector <double> > &popnum_t) const 
{
	const auto &tra = sp.tra_gl[ev.tr_gl];
	auto m = tra.nm_trans_ref;
	const auto &nmt = sp.nm_trans[m];
	const auto &ref = nmt.dist_param_eq_ref;
	auto ti = ev.ti;

	switch(tra.type){
	case EXP_RATE: case EXP_MEAN: emsg("SHould not be"); return false;
	case EXP_RATE_NM: case EXP_MEAN_NM: return false;
	case GAMMA: case LOG_NORMAL:
		{
			auto cv = eqn[ref[1]].calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
			if(cv < LINK_CV_THRESH) return true;
		}
		return false;
	
	case ERLANG:
		{
			auto shape = eqn[ref[1]].calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
			if(1.0/sqrt(shape) < LINK_CV_THRESH) return true;
			return true;
		}
		return false;
		
	case WEIBULL: 
		{
			auto shape = eqn[ref[1]].calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
			if(shape > LINK_WEIBULL_SHAPE) return true;
		}
		return false;
		
	case PERIOD: return true;
	}
	
	return true;
}



/// Determines if individual is started by a source
bool StateSpecies::source_ind(unsigned int i) const 
{
	const auto &ind = individual[i];
	if(ind.ev.size() > 0){
		if(ind.ev[0].type == M_TRANS_EV) return true;
	}
	return false;
}


/// Works out how many source individual 
unsigned int StateSpecies::source_num(unsigned int min, unsigned int max) const
{
	auto num = 0u;
	for(auto i = min; i < max; i++){
		if(source_ind(i)) num++;
	}	
	return num;
}



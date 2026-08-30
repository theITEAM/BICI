// Routines to make enforce an legal set of events to be consitent with observations

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm> 
 
using namespace std;

#include "state.hh"
#include "utils.hh"
#include "ind_ev_sampler.hh"

/// Used to order total_obs
bool total_obs_ord(const TotalObs &to1, const TotalObs &to2) 
{ return (to1.tdiv < to2.tdiv); };  


/// This forces a solution for an individual
vector <IndRef> State::force_consistent_solution()
{
	auto pl = false;

	vector <IndRef> not_possible_list;

	string warn;
	
	const auto prob_factor = 0.5; // This is used to promote event sequences with fewer events
	
	for(auto p = 0u; p < model.species.size(); p++){
		const auto &sp = model.species[p];
	
		if(sp.type == INDIVIDUAL){
			auto &ssp = species[p];
			
			vector <unsigned int> ind_list;
			
			for(auto i = 0u; i < sp.nindividual_in; i++){
				if(ssp.inconsistent(i)) ind_list.push_back(i);	
				//ind_list.push_back(i); // change
			}
			//cout << ind_list.size() << " force" << endl;
			
			//ind_list.push_back(0);
			//ind_list.clear();
			
			if(ind_list.size() == 0){
				ssp.inconsistent_ind = false;
			}
			else{
				// Works out regions where certain transitions are surpressed
				
				vector <unsigned int> apply_illegal_trans;
				vector < vector <bool> > illegal_trans;
				
				if(sp.obs_trans_exist){ 
					const auto &obs_trans_eqn_value = ssp.obs_trans_eqn_value;
					
					auto M = obs_trans_eqn_value.size();
					
					illegal_trans.resize(T);
					for(auto ti = 0u; ti < T; ti++){
						illegal_trans[ti].resize(M);
						for(auto m = 0u; m < M; m++){
							if(obs_trans_eqn_value[m][ti] == LI_WRONG) illegal_trans[ti][m] = true;
							else illegal_trans[ti][m] = false;
						}
					}
					
					if(false){
						for(auto m = 0u; m < M; m++){
							for(auto ti = 0u; ti < T; ti++){
								cout << illegal_trans[ti][m] << ",";
							}
							cout << " ill" << endl;
						}
						emsg("ill");
					}
					
					vector <unsigned int> list;
					for(auto i = 0u; i < sp.obs_trans.size(); i++){
						auto &ot = sp.obs_trans[i];
						list.push_back(ot.ti_min);
						list.push_back(ot.ti_max);
					}
					sort(list.begin(),list.end());
					
					auto j = 0u;
					while(j < list.size()){
						auto ti = list[j];
						if(ti < T) apply_illegal_trans.push_back(ti);
						while(j+1 < list.size() && list[j+1] == list[j]) j++;
						j++;
					}
				}	
				
				const auto &comp_gl = sp.comp_gl;
				const auto &tra_gl = sp.tra_gl;
				
				vector <unsigned int> map_co(comp_gl.size(),UNSET);
				
				for(auto i : ind_list){
					auto &ind = ssp.individual[i];
					const auto &indi = sp.individual[i];
					
					if(pl) cout << endl << ind.name << " individual not consistent <<<<<<<<<<<<<<<<" << endl;
					
					// Combines all the observations together
					
					vector <TotalObs> total_obs;
					
					for(auto j = 0u; j < indi.ev.size(); j++){   // Leave and move
						const auto &ev = indi.ev[j];
						if(ev.type != ENTER_EV){
							TotalObs ot;
							ot.type = EVENT_TOTAL;
							ot.index = j;
							ot.tdiv = ev.tdiv;
							ot.pos_tr = false;
							total_obs.push_back(ot);
						}
					}
						
					for(auto j = 0u; j < indi.obs.size(); j++){   // Individual observations
						TotalObs ot;
						ot.type = OBS_TOTAL;
						ot.index = j;
						ot.tdiv = indi.obs[j].tdiv;
						ot.pos_tr = false;
						total_obs.push_back(ot);
					}
				
					auto tobs_max = 0.0;
					for(auto &ot : total_obs){
						if(ot.tdiv > tobs_max) tobs_max = ot.tdiv;
					}
					
					auto ill_tra_ti = UNSET;
					vector <unsigned int> apply_illegal_gillespie;
					
					for(auto ti : apply_illegal_trans){       // Illegal transitions
						if(ti == 0) ill_tra_ti = ti;
						else{
							if(ti < tobs_max){
								TotalObs ot;
								ot.type = ILLEGAL_TRANS_TOTAL;
								ot.index = UNSET;
								ot.tdiv = ti;
								ot.pos_tr = false;
								total_obs.push_back(ot);
							}
							else{
								apply_illegal_gillespie.push_back(ti);
							}
						}
					}
					
					sort(total_obs.begin(),total_obs.end(),total_obs_ord);
					
					vector <CompProb> c_enter_pos;
					vector <CompProb> c_source_pos;
					
					if(ind.init_c_set){
						const auto &enter = sp.enter[indi.enter_ref];
						CompProb cp; cp.c = enter.c_set; cp.prob = 1;
						c_enter_pos.push_back(cp);
					}
					else{
						const auto &enter = sp.enter[indi.enter_ref];
						if(enter.set){
							vector < vector <double> > prob;
							prob.resize(sp.ncla);
							
							for(auto cl = 0u; cl < sp.ncla; cl++){
								const auto &cla = sp.cla[cl];
								prob[cl].resize(cla.comp.size(),0);
								const auto &ecl = enter.cla[cl];
								if(ecl.c_set != UNSET) prob[cl][ecl.c_set] = 1;
								else{
									for(auto cc = 0u; cc < cla.comp.size(); cc++){
										prob[cl][cc] = ssp.obs_eqn_value[ecl.obs_eqn_ref[cc]];
									}
								}
							}
							
							for(auto c = 0u; c < sp.comp_gl.size(); c++){
								const auto &cgl = sp.comp_gl[c];
								auto prob_all = 1.0;
								for(auto cl = 0u; cl < sp.ncla; cl++) prob_all *= prob[cl][cgl.cla_comp[cl]];
								if(prob_all > 0){
									CompProb cp; cp.c = c; cp.prob = prob_all;
									c_enter_pos.push_back(cp);
								}
							}
						}
						else{ // If entry not set then all possible values are possible
							for(auto c = 0u; c < sp.comp_gl.size(); c++){
								CompProb cp; cp.c = c; cp.prob = 1;
								c_enter_pos.push_back(cp);
							}
						}
					}
					
					if(c_enter_pos.size() == 0) emsg("center pos");
					
					auto kmax = total_obs.size();
					
					for(auto k = 0u; k < kmax; k++){
						auto &to = total_obs[k];
						
						auto &node = to.node;
						auto &pos_tr_list = to.pos_tr_list;
						auto &pos_co_list = to.pos_co_list;
						
						if(k == 0){   // Works out possiblities for compartment individual could enter
							for(const auto &cp : c_enter_pos){
								auto c = cp.c;
							
								ConsistNode no;
								no.c = c;
								no.prob = cp.prob;
								node.push_back(no);
							}
						}
						else{
							const auto &to_last = total_obs[k-1];
							if(to_last.pos_tr){
								for(const auto &ptl : to_last.pos_tr_list){
									auto c = tra_gl[ptl.trg].f;
									
									auto prob_new = to_last.node[ptl.node].prob*prob_factor;
									
									auto nn = map_co[c];
									if(nn == UNSET){
										map_co[c] = node.size();
										
										ConsistNode no;
										no.c = c;
										no.prob = prob_new;
										node.push_back(no);
									}
									else{
										if(prob_new > node[nn].prob){
											node[nn].prob = prob_new;
										}
									}
								}
							}
							else{
								for(const auto &pc : to_last.pos_co_list){
									auto c = pc.c;
				
									auto prob_new = pc.prob;
									
									auto nn = map_co[c];
									if(nn == UNSET){
										map_co[c] = node.size();
										
										ConsistNode no;
										no.c = c;
										no.prob = prob_new;
										node.push_back(no);
									}
									else{
										if(prob_new > node[nn].prob) node[nn].prob = prob_new;
									}
								}
							}
							
							for(auto &no : node) map_co[no.c] = UNSET;
						}
				
						auto n = 0u;
						while(n < node.size()){
							auto c = node[n].c;
							for(auto tr : comp_gl[c].tr_leave){
								if(!sp.trans_illegal(tr,ill_tra_ti,illegal_trans)){ 
									auto cf = tra_gl[tr].f;
									if(cf != UNSET){
										ConsistTrans cfr; cfr.node = n; cfr.trg = tr;
										
										auto &no = node[n];
										
										auto nn = map_co[cf];
										auto prob_new = no.prob*prob_factor;
										if(nn == UNSET){
											ConsistTrans ctr; ctr.node = node.size(); ctr.trg = tr;
											no.leave.push_back(ctr);
										
											map_co[cf] = node.size();
										
											ConsistNode no_new;
											no_new.c = cf;
											no_new.prob = prob_new;
											no_new.from.push_back(cfr);
										
											node.push_back(no_new);				
										}
										else{
											ConsistTrans ctr; ctr.node = nn; ctr.trg = tr;
											no.leave.push_back(ctr);
											
											node[nn].from.push_back(cfr);
											if(prob_new > node[nn].prob){
												node[nn].prob = prob_new;
											}
										}
									}
								}
							}
							n++;
						}
							
						switch(to.type){
						case EVENT_TOTAL:
							{
								const auto &ev = indi.ev[to.index];
								switch(ev.type){
								case LEAVE_EV:
									for(auto n = 0u; n < node.size(); n++){
										const auto &no = node[n];
										auto c = no.c;
									
										ConsistComp pc; pc.node = n; pc.c = c; pc.prob = 1;
										pos_co_list.push_back(pc);
									}
									break;
									
								case MOVE_EV:
									{
										const auto &ev = indi.ev[to.index];
										auto cl = ev.cl;
										auto move_c = ev.move_c;
										for(auto n = 0u; n < node.size(); n++){
											const auto &no = node[n];
											auto c = no.c;
										
											c = sp.get_comp_global_convert(c,cl,move_c);
									
											ConsistComp pc; pc.node = n; pc.c = c; pc.prob = 1;
											pos_co_list.push_back(pc);
										}
									}
									break;
							
								default: emsg("option not supported"); break;
								}
							}
							break;
							
						case OBS_TOTAL:
							{
								const auto &ob = indi.obs[to.index];
								switch(ob.type){
								case OBS_SOURCE_EV:
									emsg("SHould not be source");
									break;
									
								case OBS_TRANS_EV: case OBS_SINK_EV:
									{
										to.pos_tr = true;
										for(auto n = 0u; n < node.size(); n++){
											const auto &no = node[n];
											auto c = no.c;
											for(auto tr : comp_gl[c].tr_leave){
												if(ssp.get_trans_obs_prob(tr,ob) != 0){
													ConsistTrans ptr; ptr.node = n; ptr.trg = tr; 
													pos_tr_list.push_back(ptr);
												}
											}
										}			
									}
									break;
									
								case OBS_COMP_EV:
									{
										for(auto n = 0u; n < node.size(); n++){
											const auto &no = node[n];
											auto c = no.c;
											auto lprob = ssp.like_comp_obs(c,ob);
											if(lprob != LI_WRONG){
												ConsistComp pc; pc.node = n; pc.c = c; pc.prob = exp(lprob)*no.prob;
												pos_co_list.push_back(pc);
											}
										}
									}
									break;
								}
							}
							break;
						
						case ILLEGAL_TRANS_TOTAL:
							{
								for(auto n = 0u; n < node.size(); n++){
									const auto &no = node[n];
									ConsistComp pc; pc.node = n; pc.c = no.c; pc.prob = 1;
									pos_co_list.push_back(pc);
								}
								ill_tra_ti = to.tdiv;
							}
							break;
						}
						
						for(const auto &no : node) map_co[no.c] = UNSET;
					}
					
					if(pl){
						cout << endl << "Total node" << endl; 	
						for(auto k = 0u; k < kmax; k++){
							cout << "SECTION " << k << " time " << total_obs[k].tdiv << ":" << endl;
							const auto &to = total_obs[k];
							if(to.pos_tr) cout << "pos_tr" << endl;
							else cout << "comp_tr" << endl;
							
							for(const auto &no : to.node){
								cout << add_escape_char(comp_gl[no.c].name) << ": ";
								cout << "prob: " << no.prob << "  ";
								
								cout << "from ";
								for(const auto &fr : no.from) cout << fr.node << " " << add_escape_char(tra_gl[fr.trg].name) << "  ";
								
								cout << "leave ";
								for(const auto &le : no.leave) cout << le.node << " " << add_escape_char(tra_gl[le.trg].name) << "  ";
								cout << endl;
							}							
					
							if(to.pos_tr_list.size() > 0){
								for(const auto &ptl : to.pos_tr_list){
									cout << add_escape_char(tra_gl[ptl.trg].name) << ", ";
								}
								cout << "possible trans" << endl;
							}
								
							if(to.pos_co_list.size() > 0){
								for(const auto &pcl : to.pos_co_list){
									cout << add_escape_char(comp_gl[pcl.c].name) << " " << pcl.prob << ", ";
								}
								cout << "possible comp" << endl;
							}
						}
					}
							
					// Constructs event sequence
					auto not_possible = false;
					
					vector <unsigned int> tr_sel(kmax,UNSET);
					vector < vector <unsigned int> > path;
					path.resize(kmax);
					
					auto cinit = UNSET;
					
					if(kmax == 0){ // If no observations then just samples in accordance with enter probability
						auto cp = comp_prob_sample(c_enter_pos);
						cinit = cp.c;
					}
					else{
						for(int k = kmax-1; k >= 0; k--){
							auto &to = total_obs[k];
							auto &node = to.node;
							
							unsigned int n;
							if(to.pos_tr){ // Selects an observed transition
								vector <ConsistTrans> list;
								for(const auto &prl : to.pos_tr_list){
									if(cinit == UNSET || tra_gl[prl.trg].f == cinit){
										list.push_back(prl);
									}
								}
											
								if(list.size() == 0){ not_possible = true; break;}
							
								// Sample accornding to probability
								auto ct = contr_sample(list,node);
								tr_sel[k] = ct.trg;
								n = ct.node;
							}
							else{
								vector <ConsistComp> list;
								for(const auto &pc : to.pos_co_list){
									if(cinit == UNSET || pc.c == cinit){
										list.push_back(pc);
									}
								}
								
								if(list.size() == 0){ not_possible = true; break;}
							
								auto cc = conco_sample(list,node);
								n = cc.node;
							}
							
							while(node[n].from.size() != 0){
								auto fr = contr_sample(node[n].from,node);
								path[k].push_back(fr.trg);
								n = fr.node;
							}
							
							cinit = node[n].c;
						}
					}
						
					if(not_possible){
						IndRef ir; ir.p = p; ir.i = i;
						not_possible_list.push_back(ir);
					}
					else{
						if(cinit == UNSET) emsg("cinit should be set");
				
						vector <Event> ev_new;
						
						{  // Adds start event
							auto ev = ssp.get_event(ENTER_EV,i,UNSET,UNSET,UNSET,cinit,ind.ev[0].tdiv,IndInfFrom());
							ev_new.push_back(ev);
						}
						
						auto t = ev_new[ev_new.size()-1].tdiv;
					
						for(auto k = 0u; k < kmax; k++){
							const auto &to = total_obs[k];
								
							auto tnext = to.tdiv;
							if(tnext <= t) emsg("Equal time");
							
							auto rate_min = 0.2/(tnext-t);
							
							auto jmax = path[k].size();  // Adds events between observations
							if(jmax > 0){ 	
								vector <double> t_list;
								vector <unsigned int> tr_list;
								for(auto j = 0u; j < jmax; j++){
									t_list.push_back(t + ((j+0.5)/jmax)*(tnext-t));
									tr_list.push_back(path[k][jmax-1-j]);
								}
								
								if(to.pos_tr){  // The gap ends in a transtion
									t_list.push_back(tnext);
									tr_list.push_back(tr_sel[k]);
								}

								vector <double> dt_list;
								auto dt_tot = 0.0;
								for(auto j = 0u; j < tr_list.size(); j++){
									auto ti = get_ti(t_list[j]);
									auto rate = calculate_basic_rate(p,tr_list[j],ti);
									if(rate < rate_min) rate = rate_min;
									
									auto dt = exp_rate_sample(rate,warn);
									if(warn != "") emsg("Problem with rate");
									
									dt_list.push_back(dt);
									dt_tot += dt;
								}
								
								if(!to.pos_tr){ // When no final transition exists then makes update_ind
									auto dt_add = (tnext-t)-dt_tot;
									if(dt_add < dt_tot/10) dt_add = dt_tot/10;
									
									dt_list.push_back(dt_add);
									dt_tot += dt_add;
								}
								
								if(dt_list.size() != jmax+1) emsg("dt list size wrong");
							
								auto sum = 0.0;
								for(auto j = 0u; j < jmax; j++){
									sum += dt_list[j];
									auto tt = t + (sum/dt_tot)*(tnext-t);
									auto tr = tr_list[j];
									
									const auto &tra = tra_gl[tr];
									auto ty = M_TRANS_EV; if(tra.nm_trans_ref != UNSET) ty = NM_TRANS_EV;
									
									auto ev = ssp.get_event(ty,i,tr,UNSET,UNSET,tra.f,tt,IndInfFrom());
									ev_new.push_back(ev);
								}
							}
							
							switch(to.type){
							case EVENT_TOTAL:
								{
									const auto &ev = indi.ev[to.index];
									switch(ev.type){
									case LEAVE_EV:
										{
											auto ev = ssp.get_event(LEAVE_EV,i,UNSET,UNSET,UNSET,UNSET,tnext,IndInfFrom());
											ev_new.push_back(ev);
										}
										break;
										
									case MOVE_EV:
										{
											auto c = ev_new[ev_new.size()-1].c_after;
											
											const auto &iev = indi.ev[to.index];
											auto cl = ev.cl;
											auto move_c = ev.move_c;
											c = sp.get_comp_global_convert(c,cl,move_c);
											auto ev = ssp.get_event(MOVE_EV,i,UNSET,iev.move_c,cl,c,tnext,IndInfFrom());
											ev_new.push_back(ev);
										}
										break;
										
									default:
										emsg("SHould not be here2");
										break;
									}
								}
								break;
								
							case OBS_TOTAL:
								{
									auto trsel = tr_sel[k];
									if(trsel != UNSET){
										const auto &tra = tra_gl[trsel];
										auto ty = M_TRANS_EV; if(tra.nm_trans_ref != UNSET) ty = NM_TRANS_EV;
										auto ev = ssp.get_event(ty,i,trsel,UNSET,UNSET,tra.f,tnext,IndInfFrom());
										ev_new.push_back(ev);
									}
								}
								break;
							
							case ILLEGAL_TRANS_TOTAL:						
								break;
							}
							
							t = tnext;
						}
					
						// Performs Gillespie algorithm to fill in any subsequent events
						if(t < T){
							if(pl) cout << t << " START GILLESPIE" << endl;
							
							auto a = 0u;
							auto amax = apply_illegal_gillespie.size();
							
							auto c = ev_new[ev_new.size()-1].c_after;
							while(c != UNSET && t < T){
								const auto &tr_leave = comp_gl[c].tr_leave;
								
								auto ti = get_ti(t);
								
								vector <double> rate_sum;
								auto R = 0.0;
								auto jmax = tr_leave.size();
								for(auto j = 0u; j < jmax; j++){
									auto tr = tr_leave[j];
									if(!sp.trans_illegal(tr,ill_tra_ti,illegal_trans)){ 
										auto rate = calculate_basic_rate(p,tr,ti);
										R += rate;
									}
									rate_sum.push_back(R);
								}
								
								if(R <= RATE_MIN) break;
								
								auto t_next = t + exp_rate_sample(R,warn);
								if(warn != "") emsg("Rate wearning:"+warn);
								
								if(a < amax && t_next >= apply_illegal_gillespie[a]){
									t = apply_illegal_gillespie[a];
									a++;
									ill_tra_ti = t;
								}
								else{
									if(t_next >= T) break;
									
									auto z = ran()*R;
									auto j = 0u;
									while(j < jmax && z > rate_sum[j]) j++;
									if(j == jmax) emsg("Problem sampling");
									
									auto trsel = tr_leave[j];
									const auto &tra = tra_gl[trsel];
									auto ty = M_TRANS_EV; if(tra.nm_trans_ref != UNSET) ty = NM_TRANS_EV;
									auto ev = ssp.get_event(ty,i,trsel,UNSET,UNSET,tra.f,t_next,IndInfFrom());
									ev_new.push_back(ev);
									
									t = t_next;
									c = tra.f;
								}
							}
						}
						
						if(pl){
							cout << "OLD event" << endl;
							ssp.print_event(ind.ev);
					
							//cout << "NEW event" << endl;
							//ssp.print_event(ev_new);
						}
						
						auto gc = update_tree(p,i,ev_new);
						if(gc.type != GENCHA_FAIL){
							auto like_ch = update_ind(p,i,ev_new,UP_SINGLE);

							if(pl){
								cout << endl;
								cout << "After:" << endl; ssp.print_event(ssp.individual[i].ev);
								cout << "Liobs after: " << ssp.Li_obs_ind[i] << endl;
							}

							if(ssp.inconsistent(i)) emsg("Still not consistent");
	
							accept(like_ch);

							//if(pl) check("make consistent");

							gen_change_update(gc); 	
							if(sp.trans_tree) update_popnum_ind(p,i);
						}
						
						if(pl){
							cout << endl;
							cout << "PATH" << endl;
							cout << "cinit: " << add_escape_char(comp_gl[cinit].name) << endl;
							for(auto k = 0u; k < kmax; k++){
								cout << "SECTION " << k << " time " << total_obs[k].tdiv << ":" << endl;
								if(path[k].size() > 0){
									cout << "path - ";
									for(auto tr : path[k]) cout << add_escape_char(tra_gl[tr].name) << "  ";
								}
								
								if(tr_sel[k] != UNSET){
									cout << "  sel trans - " << add_escape_char(tra_gl[tr_sel[k]].name) << "  ";
								}
								cout << endl;
							}
							//emsg("zpath");
						}
					}
				}
			
				if(testing){
					for(auto va : map_co){
						if(va != UNSET) emsg("Problem with map_co");
					}
				}
			}
		}
	}
	
	return not_possible_list;
}


/// Calculates the characteristic transition time for a given transition
double State::calculate_basic_rate(unsigned int p, unsigned int tr, unsigned int ti) const
{
	const auto &sp = model.species[p];
	const auto &ssp = species[p];
	auto dt = model.details.dt;
	
	const auto &tra = sp.tra_gl[tr];
	
	switch(tra.type){
	case EXP_RATE:
		{
			const auto &val_t = ssp.markov_eqn_vari[tra.markov_eqn_ref].value_t;
			if(val_t.size() == 1) return val_t[0];
			else return val_t[ti];
		}
		break;
	
	case EXP_MEAN: 
		{
			const auto &val_t = ssp.markov_eqn_vari[tra.markov_eqn_ref].value_t;
			if(val_t.size() == 1) ti = 0;
			auto val = val_t[ti];
			if(ti < TINY) ti = TINY;
			return 1.0/val;
		}
		break;

	case EXP_RATE_NM:
		{
			const auto &eq = model.eqn[tra.dist_param[0].eq_ref];
			auto rate = eq.calculate(ti,popcomb_t[ti],param_val.precalc);
			return dt*rate;
		}
		break;
		
	case GAMMA: case ERLANG: case LOG_NORMAL: case WEIBULL: case PERIOD: case EXP_MEAN_NM: 
		{
			const auto &eq = model.eqn[tra.dist_param[0].eq_ref];
			auto mean = eq.calculate(ti,popcomb_t[ti],param_val.precalc);
			if(mean < TINY) mean = TINY;
			return dt/mean;
		}	
		break;
	}
	
	emsg("Rate not set");
	return UNSET;
}


/// Resamples individual using the observation sampler (this gets fixed events correct)
void State::resample_ind(bool if_wrong)
{
	const auto &precalc = param_val.precalc;
	
	auto pl = false; 
	//if(if_wrong) pl = true;
	
	if(pl) check(" before resample");
			
	for(auto p = 0u; p < species.size(); p++){
		auto &sp = model.species[p];
		if(sp.type == INDIVIDUAL){
			auto &ssp = species[p];
		
			vector <unsigned int> ind_list;
		
			for(auto i = 0u; i < sp.nindividual_in; i++){
				if(if_wrong == false || ssp.inconsistent(i)) ind_list.push_back(i);
			}
			//cout << ind_list.size() << " resample" << endl;
			//ind_list.push_back(500);
			
			if(ind_list.size() > 0){
				IndEvSampler ind_ev_samp(ssp.markov_eqn_vari,ssp.individual,model.details,sp,ssp.obs_eqn_value,ssp.obs_trans_eqn_value,model.eqn,genetic_value.inf_node,precalc,popcomb_t,if_wrong);

				ind_ev_samp.setup_nm();
			
				for(auto i : ind_list){
					//cout << i << " " << 
					//cout << i << " " << ssp.individual[i].name << " resample" << endl;
					
					if(pl) cout << i << " " << ssp.individual[i].name << " resample" << endl;
					for(auto cl = 0u; cl < sp.ncla; cl++){
						if(ind_ev_samp.needed(i,cl) == true){
							if(false && pl){
								cout << endl << endl << endl << endl;
								cout << ssp.individual[i].name << endl;
								cout << "Before:" << endl; ssp.print_event(ssp.individual[i].ev);
								cout << ssp.Li_obs_ind[i] << " Liobs befor" << endl;
							}
							
							ind_ev_samp.generate_ind_obs_timeline();
							//ind_ev_samp.print_ind_obs_timeline();
						
							auto probif = 0.0;
							auto ev_new = ind_ev_samp.sample_events(probif);

							if(pl){
								cout << endl << endl << endl << endl;
								cout << ssp.individual[i].name << endl;
								cout << "Before:" << endl; ssp.print_event(ssp.individual[i].ev);
								cout << endl;
								cout << "Propose:" << endl; ssp.print_event(ev_new);
								cout << "Liobs befor: " << ssp.Li_obs_ind[i] << endl;
					
								cout << "Illegal:" << ind_ev_samp.illegal << endl;
							}
					
							if(ind_ev_samp.illegal == false){
								auto gc = update_tree(p,i,ev_new);
								if(gc.type != GENCHA_FAIL){
									auto like_ch = update_ind(p,i,ev_new,UP_SINGLE);

									if(pl){
										cout << endl;
										cout << "After:" << endl; ssp.print_event(ssp.individual[i].ev);
										cout << "Liobs after: " << ssp.Li_obs_ind[i] << endl;
										
										/*
										if(ssp.inconsistent(i)){
											ind_ev_samp.print_ind_obs_timeline();
											ind_ev_samp.pr_generate_ind_obs_timeline();
											
											model.print_param(param_val);
											emsg("done");
										}
										*/
									}
		
									accept(like_ch);

									gen_change_update(gc); 	
									if(sp.trans_tree) update_popnum_ind(p,i);
								
									//if(pl) check("during resample");
									
									if(std::isnan(like.markov)){
										cout << endl << endl << endl << endl;
										cout << ssp.individual[i].name << endl;
										cout << "Before:" << endl; ssp.print_event(ssp.individual[i].ev);
										cout << endl;
										cout << "Propose:" << endl; ssp.print_event(ev_new);
								
										emsg("markov is nan");
									}
									//emsg("done");
									//if(ssp.individual[i].name == "GLOC1") emsg("done");
									break;
								}
							}	
						}
					}
				}
			}	
		}
	}
}


/// Samples from a list of comp_prob
CompProb State::comp_prob_sample(const vector <CompProb> &cp_pos) const 
{
	vector <double> sum_st;
	auto sum = 0.0; 
	auto jmax = cp_pos.size();
	for(auto j = 0u; j < jmax; j++){
		sum += cp_pos[j].prob;
		sum_st.push_back(sum);
	}
	if(sum == 0) emsg("Zero probability of entering");
	auto z = ran()*sum;
	auto j = 0u; while(j < jmax && z > sum_st[j]) j++;
	if(j == jmax) emsg("Problem sample");
	
	return cp_pos[j];
}


/// Samples from a list of possible
ConsistTrans State::contr_sample(const vector <ConsistTrans> &list, const vector <ConsistNode> &node) const
{
	vector <double> sum_st;								
	auto jmax = list.size();
	auto sum = 0.0;
	for(auto j = 0u; j < jmax; j++){
		sum += node[list[j].node].prob;
		sum_st.push_back(sum);
	}
	if(sum == 0) emsg("Sample problem");

	auto z = ran()*sum;
	auto j = 0u; while(j < jmax && z > sum_st[j]) j++;
	if(j == jmax) emsg("problem jmax");
								
	return list[j];
}


/// Samples from a list of possible
ConsistComp State::conco_sample(const vector <ConsistComp> &list, const vector <ConsistNode> &node) const
{
	vector <double> sum_st;								
	auto jmax = list.size();
	auto sum = 0.0;
	for(auto j = 0u; j < jmax; j++){
		sum += node[list[j].node].prob*list[j].prob;
		sum_st.push_back(sum);
	}
	if(sum == 0) emsg("Sample problem");

	auto z = ran()*sum;
	auto j = 0u; while(j < jmax && z > sum_st[j]) j++;
	if(j == jmax) emsg("problem jmax");
								
	return list[j];
}

struct IndNum {
	string name;
	unsigned int num;
};


/// Used to order genetic data
bool IndNum_ord (const IndNum &in1, const IndNum &in2)
{ return (in1.num > in2.num); };


/// This gets called if there is not valid starting condition for the chain
void State::no_valid_state(const vector < vector <IndRef> > &not_possible) const
{
	if(not_possible.size() > 0){ // Looks for individuals that consistently fail
		vector <IndNum> ind_num;
		
		Hash hash;
		for(auto i = 0u; i < not_possible.size(); i++){
			for(const auto &np : not_possible[i]){
				auto name = model.species[np.p].individual[np.i].name;
				auto j = hash.find(name);
				if(j == UNSET){
					hash.add(ind_num.size(),name);
					IndNum in;
					in.name = name;
					in.num = 1;
					ind_num.push_back(in);
				}
				else{
					ind_num[j].num++;
				}
			}
		}
		
		sort(ind_num.begin(),ind_num.end(),IndNum_ord);
		
		auto max_len = 100u;
		string prob_ind;
		for(auto &in : ind_num){
			if(prob_ind != "") prob_ind += ", ";
			prob_ind += in.name;
			if(prob_ind.length() > max_len) break;
		}
		
		run_error("Could not find consistent initial state for the following individuals (please check individual-based data): "+trunc(prob_ind,max_len));
	}
	
	run_error("Could not find consistent initial state");	
}

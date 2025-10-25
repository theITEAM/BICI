// This provides information about a single species within the model
// All properties of Species are fixed during simulation or inference

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
using namespace std;

#include "species.hh"
#include "utils.hh"
#include "matrix.hh"


/// Initialises the model
Species::Species(unsigned int p, const vector <double> &timepoint, const Details &details) : timepoint(timepoint), details(details)
{
	p_species = p;
	obs_trans_exist = false;
	add_rem_pop_on = false;
	pop_trans_data_exist = false;
};


/// Creates a node tree used to sample Markovian events (used during simulation)
void Species::create_markov_tree()
{
	for(auto i = 0u; i < markov_eqn.size(); i++){
		MarkovNode mn; mn.parent = UNSET;
		markov_tree.node.push_back(mn);
	}
	
	do{
		vector <unsigned int> list;
		for(auto i = 0u; i < markov_tree.node.size(); i++){
			if(markov_tree.node[i].parent == UNSET) list.push_back(i);
		}
		
		if(list.size() <= 1) break;
		for(auto j = 0u; j < list.size(); j += 2){	
			MarkovNode mn; 
			mn.parent = UNSET;
			
			auto i = list[j];
			markov_tree.node[i].parent = markov_tree.node.size();
			mn.child.push_back(i);
			
			if(j+1 < list.size()){
				i = list[j+1];
				markov_tree.node[i].parent = markov_tree.node.size();
				mn.child.push_back(i);
			}
			
			markov_tree.node.push_back(mn);
		}
	}while(true);
	
	if(false){
		for(auto i = 0u; i < markov_tree.node.size(); i++){
			const auto &no = markov_tree.node[i];
			cout << "Node " << i << ": parent " << no.parent << "  ";
			cout << "child "; for(auto j : no.child) cout << j << ", ";
			cout << endl;		
		}
	}
}


/// Updates the global compartment based on a classification and comp being set to a given value
unsigned int Species::update_c_comp(unsigned int c, unsigned int cl, unsigned int c_comp) const
{
	int cco = comp_gl[c].cla_comp[cl];
	int cnew = c +(c_comp - cco)*comp_mult[cl];
	if(cnew < 0 || cnew >= (int)N) emsg("Out of range4");
	return (unsigned int) cnew;
}


/// Potentially corrects enew based on current compartment
bool Species::correct_ev(unsigned int c, Event &enew) const
{
	switch(enew.type){
	case ENTER_EV: 
		c = enew.c_after;
		break;
	
	case LEAVE_EV:
		c = UNSET;
		enew.c_after = c;
		break;

	case MOVE_EV: 
		{
			auto cl = enew.cl, c_comp = enew.move_c;
			c = update_c_comp(c,cl,c_comp);
			enew.c_after = c;
		}
		break;
		
		
	case NM_TRANS_EV:
	case M_TRANS_EV:
		{
			auto tr = enew.tr_gl;
			const auto &trg = tra_gl[tr];
			if(trg.i != c){                               // Ensures that sequence is consistent
				if(c == UNSET){
					emsg("Not consistent4");
					return true;
				}
		
				tr = tr_trans(tr,c); if(tr == UNSET) return true;
				enew.tr_gl = tr;
				auto c_after = tra_gl[tr].f;
				enew.c_after = c_after;
				c = c_after;
			}
			else c = trg.f;
		}
		break;
	}
	
	return false;
}


/// Creates a "fake" individual a the end of the list 
// (used when dealing with unobserved individuals)
void Species::create_ind_noobs()
{
	nindividual_obs = individual.size();
	if(type == INDIVIDUAL){
		IndData ind;
		ind.name = "No observation";
		ind.enter_ref = UNSET;
		ind.tdivmin = -LARGE;
		ind.tdivmax = -LARGE;
		individual.push_back(ind);
	}
}


/// Determines if the population is fixed
bool Species::is_pop_num_fixed() const
{
	auto fl = true;                                   // Sets if the total pop number is fixed
	if(init_cond.type == INIT_POP_DIST){	
		if(init_cond.focal_cl == UNSET){
			if(init_cond.pop_prior.type != FIX_PR) fl = false;
		}
		else{
			const auto &claa = cla[init_cond.focal_cl];
			for(auto c = 0u; c < claa.comp.size(); c++){
				if(init_cond.comp_prior[c].type != FIX_PR) fl = false;
			}
		}
	}
	
	return fl;
}


/// Sets up nm_flag and cgl_incomp_nmtrans
void Species::nm_trans_incomp_info()
{
	nm_flag = false;
	for(auto &tr : tra_gl){
		switch(tr.type){
			case EXP_RATE_NM: case EXP_MEAN_NM: case GAMMA: case ERLANG: case LOG_NORMAL: 
			case WEIBULL: case PERIOD:
				nm_flag = true;
				break;
			default: break;
		}
	}
	
	cgl_incomp_nmtrans_cl.resize(comp_gl.size());
	
	for(auto c = 0u; c < comp_gl.size(); c++){
		cgl_incomp_nmtrans_cl[c].resize(ncla);
		
		const auto &cgl = comp_gl[c];
		for(auto cl = 0u; cl < ncla; cl++){
			if(cgl.tra_leave_group[cl].markov == false){
				cgl_incomp_nmtrans_cl[c][cl] = true;
			}
		}			
	}
}


/// Gets nm transition number from event 
unsigned int Species::get_tra_m(unsigned int tr, const Event &ev_orig) const
{
	const auto &tra = tra_gl[tr];
	if(tra.i != ev_orig.c_after){ // accounts for intermediate transition
		auto trg = tr_trans(tr,ev_orig.c_after);
		if(trg == UNSET){
			emsg("trg unset1");
		}
		return tra_gl[trg].nm_trans_ref;
	}
	return tra.nm_trans_ref;
}
	
	
/// Finds a name for a sequence of events
string Species::tr_swap_name(unsigned int cl, const vector <TrSwap> &tswa) const 
{
	const auto &claa = cla[cl];
	
	stringstream ss;
	for(auto k = 0u; k < tswa.size(); k++){
		const auto &tsw = tswa[k];
		switch(tsw.type){
		case ENTER_SW: ss << "Enter->" << claa.comp[tsw.trc].name; break;
		case LEAVE_SW: ss << claa.comp[tsw.trc].name << "->Leave"; break;
		case SOURCE_SW:	ss << "+->" << claa.tra[tsw.trc].name; break;
		case SINK_SW: ss << claa.tra[tsw.trc].name << "->-"; break;
		case TRANS_SW: ss << claa.tra[tsw.trc].name; break;
		case MOVE_SW: ss << "MOVE"; break;
		}
		ss << " ... ";
	}
	
	return ss.str();
}


/// Converts information about a potential local individual swap of events to a hash vec
vector <unsigned int> Species::get_vec_tr_swap(unsigned int c, const vector <TrSwap> &start) const
{
	vector <unsigned int> vec;
	vec.push_back(c);
	for(const auto &trs : start){
		vec.push_back(trs.type);
		vec.push_back(trs.trc);
	}
	
	return vec;
}


/// Converts information about a potential local individual swap of events to a hash vec
vector <unsigned int> Species::get_vec_tr_swap_mid(unsigned int st, unsigned int num, const vector <EventCl> &timeline) const
{
	vector <unsigned int> vec;
	vec.push_back(timeline[st].c_bef);
	
	for(auto j = st; j < st+num; j++){
		const auto &trs =	timeline[j].trs;
		vec.push_back(trs.type);
		vec.push_back(trs.trc);
	}
	
	return vec;
}


/// Calcucalculate_no_popnum(lates non-Markovian rate
vector < vector <double> > Species::calc_nm_rate(bool calc_bp, const vector <double> &precalc, const vector < vector <double> > &popnum_t, const vector <Equation> &eqn, vector < vector <double> > &bp_store) const
{
	vector < vector <double> > nm_rate;
	
	auto dt = details.dt;
	
	auto M = nm_trans.size();
	nm_rate.resize(M);
	for(auto m = 0u; m < M; m++){
		nm_rate[m].resize(T);
		
		const auto &nmt = nm_trans[m];
	
		if(nmt.precalc_nm_rate || calc_bp){		
			// Calculates branching probabilities
			vector <double> bp(T,1);
			
			auto bp_eq = nmt.bp_eq;
			if(bp_eq != UNSET){
				if(bp_eq == BP_FROM_OTHERS){
					auto time_vari = false;
					for(auto e : nmt.bp_other_eq){
						if(eqn[e].time_vari) time_vari = true;
					}
							
					if(time_vari){
						for(auto ti = 0u; ti < T; ti++){
							const auto &popnum = popnum_t[ti];
							auto bp_f = 1.0;
							for(auto e : nmt.bp_other_eq){
								bp_f -= eqn[e].calculate(ti,popnum,precalc);
							}
							check_bp(bp_f);
							bp[ti] = bp_f;
						}
					}
					else{
						auto bp_f = 1.0;
						for(auto e : nmt.bp_other_eq){
							bp_f -= eqn[e].calculate_param(precalc);
						}
						check_bp(bp_f);
						for(auto ti = 0u; ti < T; ti++) bp[ti] = bp_f;
					}
				}
				else{
					auto update = UPDATE_NO_TIME;
					if(nmt.all_branches){
						for(auto e : nmt.bp_all_eq) if(eqn[e].time_vari) update = UPDATE_TIME;
					}
					else{
						if(eqn[bp_eq].time_vari) update = UPDATE_TIME;
					}
					
					const auto &eq = eqn[bp_eq];
					if(update == UPDATE_TIME){
						for(auto ti = 0u; ti < T; ti++){
							const auto &popnum = popnum_t[ti];
							auto bp_f = eq.calculate(ti,popnum,precalc);
							if(nmt.all_branches){
								auto div = 0.0;
								for(auto e : nmt.bp_all_eq){
									div += eqn[e].calculate(ti,popnum,precalc);
								}
								bp_f /= div;
							}
							check_bp(bp_f);
							bp[ti] = bp_f;
						}
					}
					else{
						auto bp_f = eq.calculate_param(precalc);
						if(nmt.all_branches){
							auto div = 0.0;
							for(auto e : nmt.bp_all_eq){
								div += eqn[e].calculate_param(precalc);
							}
							bp_f /= div;
						}
						check_bp(bp_f);
						for(auto ti = 0u; ti < T; ti++) bp[ti] = bp_f;
					}
				}
			}
			if(calc_bp) bp_store.push_back(bp);
			
			switch(nmt.type){
			case EXP_RATE: case EXP_MEAN: emsg("Should not be here"); break;

			case GAMMA: case ERLANG: case LOG_NORMAL: case PERIOD: case EXP_MEAN_NM:
				{
					const auto &eq = eqn[nmt.dist_param_eq_ref[0]];
					if(eq.time_vari){
						for(auto ti = 0u; ti < T; ti++){	
							nm_rate[m][ti] = dt*bp[ti]/eq.calculate(ti,popnum_t[ti],precalc);
						}
					}
					else{
						auto val = 1.0/eq.calculate_param(precalc); 
						for(auto ti = 0u; ti < T; ti++){
							nm_rate[m][ti] = dt*bp[ti]*val;
						}
					}
				}
				break;

			case EXP_RATE_NM:
				{
					const auto &eq = eqn[nmt.dist_param_eq_ref[0]];
					if(eq.time_vari){
						for(auto ti = 0u; ti < T; ti++){	
							nm_rate[m][ti] = dt*bp[ti]*eq.calculate(ti,popnum_t[ti],precalc);
						}		
					}
					else{
						auto val = eq.calculate_param(precalc); 
						for(auto ti = 0u; ti < T; ti++){
							nm_rate[m][ti] = dt*bp[ti]*val;
						}
					}
				}
				break;
			
			case WEIBULL:
				{
					// The gamma function correction to mean is relatively small so here it is ignored
					const auto &eq = eqn[nmt.dist_param_eq_ref[0]];	
					if(eq.time_vari){
						for(auto ti = 0u; ti < T; ti++){	
							nm_rate[m][ti] = dt*bp[ti]/eq.calculate(ti,popnum_t[ti],precalc);
						}
					}
					else{
						auto val = 1.0/eq.calculate_param(precalc); 
						for(auto ti = 0u; ti < T; ti++){
							nm_rate[m][ti] = dt*bp[ti]*val;
						}
					}
					
					if(false){ // This more complicate version didn't perform much better
						const auto &eq = eqn[nmt.dist_param_eq_ref[0]];
						const auto &eq_sh = eqn[nmt.dist_param_eq_ref[1]];
					
						if(eq.time_vari || eq_sh.time_vari){
							for(auto ti = 0u; ti < T; ti++){
								nm_rate[m][ti] = dt/eq.calculate_no_popnum(ti,precalc);
							}
						}
						else{
							auto val = dt/(eq.calculate_param(precalc)); 
							
							for(auto ti = 0u; ti < T; ti++) nm_rate[m][ti] = val;
						}
					}
				}
				break;
				
			default: emsg("op not here2"); break;
			}
		}
		else{
			for(auto ti = 0u; ti < T; ti++) nm_rate[m][ti] = UNSET;
		}
	}
	
	return nm_rate;
}


/// Determines if individual samplers are needed (i.e. if all transitions are known)
void Species::set_ind_samp_needed(const vector <Equation> &eqn)
{
	// Determines if transition is not observed at any point over the entire time period
	vector <bool> tra_not_obs(tra_gl.size(),false);
	for(auto trg = 0u; trg < tra_gl.size(); trg++){
		auto not_obs = false;
		for(auto ti = 0u; ti < T; ti++){
			auto fl = false;
			for(auto ter : obs_trans_eqn_ref[trg][ti]){
				const auto &ot = obs_trans[ter];
				if(ot.is_one[trg] && ot.single_trans){ fl = true; break;}
			}
			
			if(fl == false){ not_obs = true; break;}
		}
		tra_not_obs[trg] = not_obs;
	}
	
	if(false){
		for(auto trg = 0u; trg < tra_gl.size(); trg++){
			cout << tra_gl[trg].name << " " << tra_not_obs[trg] << " not obs" << endl;
		}
		emsg("done");
	}

	// This stores if there is a link between a classification compartment and trans
	vector < vector < vector <bool> > > comp_tra_link; // [cl][c][tr]
	comp_tra_link.resize(ncla);
	for(auto cl = 0u; cl < ncla; cl++){
		const auto &claa = cla[cl];
		comp_tra_link[cl].resize(claa.ncomp);
		for(auto c = 0u; c < claa.ncomp; c++){
			comp_tra_link[cl][c].resize(claa.ntra,false);
		}
	}

	for(auto trg = 0u; trg < tra_gl.size(); trg++){
		const auto &tra = tra_gl[trg];
		auto cl = tra.cl;
		
		auto c = cla[cl].tra[tra.tr].i;
		if(c != UNSET){
			comp_tra_link[cl][c][tra.tr] = true;
			
			for(auto k = 0u; k < tr_after[trg].size(); k++){
				auto trg_aft = tr_after[trg][k];
				const auto &tra_aft = tra_gl[trg_aft];
				comp_tra_link[cl][c][tra_aft.tr] = true;
			}
		}
	}		
	
	if(false){
		for(auto cl = 0u; cl < ncla; cl++){
			const auto &claa = cla[cl];
			for(auto c = 0u; c < claa.ncomp; c++){
				for(auto tr = 0u; tr < claa.ntra; tr++){
					cout << claa.comp[c].name << " " << claa.tra[tr].name << " " 
							<< comp_tra_link[cl][c][tr] << endl;
				}
			}
		}
		emsg("done");
	}
	
	// Makes a list of global not obseved transitions stratified by cl
	vector < vector <unsigned int> > trg_not_obs_cl; 
	trg_not_obs_cl.resize(ncla);
	for(auto trg = 0u; trg < tra_gl.size(); trg++){
		const auto &tra = tra_gl[trg];
		if(tra_not_obs[trg]) trg_not_obs_cl[tra.cl].push_back(trg);
	}

	if(false){
		for(auto cl = 0u; cl < ncla; cl++){
			for(auto trg : trg_not_obs_cl[cl]) cout << cl << " "<< tra_gl[trg].name << endl;
		}
		emsg("done");
	}
	
	for(auto i = 0u; i < individual.size(); i++){
		auto &ind = individual[i];
		if(i < nindividual_in){	
			ind.sample_needed.resize(ncla,true);
				
			auto e = ind.enter_ref;

			auto num_pos = 0u;

			if(e != UNSET){
				const auto &ent = enter[e];
				
				vector <bool> cgl_pos(comp_gl.size(),false);
				if(ent.c_set != UNSET){ cgl_pos[ent.c_set] = true; num_pos++;}
				else{
					vector < vector <bool> > c_pos;
					c_pos.resize(ncla);
					for(auto cl = 0u; cl < ncla; cl++){
						const auto &claa = cla[cl];
						c_pos[cl].resize(claa.ncomp,false);
						const auto &ent_cl = ent.cla[cl];
						if(ent_cl.c_set != UNSET) c_pos[cl][ent_cl.c_set] = true;
						else{
							for(auto c = 0u; c < claa.ncomp; c++){
								if(!eqn[ent_cl.eqn[c].eq_ref].is_zero()) c_pos[cl][c] = true;
							}
						}
					}
					
					for(auto c = 0u; c < comp_gl.size(); c++){
						const auto &cla_comp = comp_gl[c].cla_comp;
						auto cl = 0u; while(cl < ncla && c_pos[cl][cla_comp[cl]] == true) cl++;
						
						if(cl == ncla){
							cgl_pos[c] = true;
							num_pos++;
						}
					}
				}

				for(auto cl = 0u; cl < ncla; cl++){
					// First determines if initial state is set within classification
	
					const auto &claa = cla[cl];
					vector <bool> non_zero(claa.ncomp,false);
					
					for(auto c = 0u; c < comp_gl.size(); c++){
						const auto &cgl = comp_gl[c];
						if(cgl_pos[c]) non_zero[cgl.cla_comp[cl]] = true;
					}
					
					auto unobs = false;
					for(auto c = 0u; c < claa.ncomp; c++){
						if(non_zero[c]){
							for(auto trg : trg_not_obs_cl[cl]){
								const auto &tra = tra_gl[trg];
								if(comp_tra_link[cl][c][tra.tr]){ unobs = true; break;}
							}
						}
						if(unobs == true) break;
					}
								
					if(unobs == false) ind.sample_needed[cl] = false;
				}
			}
			
			ind.simulation_needed = false;
			for(auto need : ind.sample_needed){
				if(need) ind.simulation_needed = true;
			}
		
			// If uncertainty in initial condition
			if(num_pos > 1) ind.simulation_needed = true; 
			
			ind.move_needed = ind.simulation_needed;
			if(contains_source) ind.move_needed = true;
		}
		else{
			ind.sample_needed.resize(ncla,false);
			ind.simulation_needed = false;
			ind.move_needed = false;
		}
	}
}


/// Calcualtes the probability of an individual entering a given state
double Species::calc_enter_prob(unsigned int c, unsigned int entref, const vector <double> &obs_eqn_value) const
{
	if(c == UNSET) return 0;
	const auto &ent = enter[entref];
	if(ent.c_set != UNSET){
		if(c == ent.c_set) return 1;
		return 0;
	}
	else{
		const auto &cla_comp = comp_gl[c].cla_comp;
		
		auto prob = 1.0;
		for(auto cl = 0u; cl < ncla; cl++){
			const auto &ent_cl = ent.cla[cl];
			auto cc = cla_comp[cl];
			if(ent_cl.c_set){
				if(cc != ent_cl.c_set) return 0;
			}
			else{
				auto val = obs_eqn_value[ent_cl.obs_eqn_ref[cc]];
				if(val < 0) emsg("Initial compartment observational likelihood has become negative"); 
				prob *= val;
			}
		}
		return prob;
	}
}


/// Gets parameters cv's which correspond individual effect in mean (used for IE_VAR_CV_PROP)
vector <unsigned int> Species::get_cv_list(unsigned int ie, const vector <Param> &param, const vector <Equation> &eqn) const
{
	auto fl = true;
	
	const auto &ind_eff = ind_effect[ie]; 
	
	vector <unsigned int> list;
	HashSimp hash_list;
	
	if(ind_eff.markov_eqn_ref.size() > 0) fl = false;
	else{
		for(auto ref : ind_eff.nm_trans_ref){
			const auto &nmt = nm_trans[ref];
			
			switch(nmt.type){
			case GAMMA: case LOG_NORMAL:
				{
					const auto &mean_eqn = eqn[nmt.dist_param_eq_ref[0]];
					const auto &cv_eqn = eqn[nmt.dist_param_eq_ref[1]];
					if(find_in(mean_eqn.ind_eff_mult,ie) == UNSET) fl = false;
					if(find_in(cv_eqn.ind_eff_mult,ie) != UNSET) fl = false;
					if(cv_eqn.param_ref.size() != 1) fl = false;
					else{
						const auto &pr = cv_eqn.param_ref[0];
						const auto &par = param[pr.th];
						if(par.variety == CONST_PARAM) fl = false;
						else{
							auto th = par.get_param_vec(pr.index);	
							if(th == UNSET) fl = false;
							else add_to_vec(list,th,hash_list);
						}
					}
				}
				break;
				
			default: fl = false; break;
			}
		}			
	}
	
	if(fl == false) list.resize(0);
	return list;
}


/// Calculates constant omega
vector < vector <double> > Species::calculate_omega_const(unsigned int g, const vector <Param> &param) const
{
	const auto &ieg = ind_eff_group[g];
	const auto &par = param[ieg.th];
	auto N = ieg.list.size();
	
	vector < vector <double> > omega;
	omega.resize(N);
	for(auto j = 0u; j < N; j++) omega[j].resize(N);
	
	for(auto j = 0u; j < N; j++){
		for(auto i = 0u; i < N; i++){
			omega[j][i] = par.get_value(j*N+i);
		}
	}
	
	return omega;
}


/// Calculate the actual omega matrix
vector < vector <double> > Species::calculate_omega_basic(unsigned int g, const PV &param_val, const vector <Param> &param) const
{
	const auto &ieg = ind_eff_group[g];		
	const auto &par = param[ieg.th];
	vector < vector <double> > omega;
	
	if(par.variety == CONST_PARAM){
		omega = calculate_omega_const(g,param);
	}
	else{
		auto N = ieg.list.size();
	
		const auto &value = param_val.value;

		omega.resize(N);
		for(auto j = 0u; j < N; j++) omega[j].resize(N);

		for(auto j = 0u; j < N; j++){
			for(auto i = 0u; i < N; i++){
				omega[j][i] = value[ieg.omega_pv[j][i]];
			}
		}
	}
	
	convert_cor_var(omega);
	
	//print_matrix("om",omega);
	
	return omega;
}	


/// Converts correlations to variances
void Species::convert_cor_var(vector < vector <double> > &omega) const
{
	auto N = omega.size();
	for(auto j = 0u; j < N; j++){
		for(auto i = 0u; i < N; i++){
			if(j != i){
				omega[j][i] *= sqrt(omega[i][i]*omega[j][j]);
			}
		}
	}
}

	
/// Shifts a global transition such that fits with a given c
unsigned int Species::tr_trans(unsigned int trg, unsigned int c) const
{
	if(tform_set) return tra_gl[trg].tform[c];
	
	const auto &tra = tra_gl[trg];
	auto cl_tra = tra.cl;
	const auto &co = comp_gl[c];
	const auto &coi = comp_gl[tra.i];
	
	for(auto cl = 0u; cl < ncla; cl++){
		int dif = co.cla_comp[cl]-coi.cla_comp[cl];
		if(dif != 0){
			if(cl == cl_tra) emsg("Should not differ in cl");
			trg += tr_shift[cl_tra][cl]*dif;
		}
	}
	if(tra_gl[trg].i != c) emsg("could not do tform");
	
	return trg;
}


/// Convert from a global compartment to another based on classification cl being in c
unsigned int Species::get_comp_global_convert(unsigned int cgl, unsigned int cl, unsigned int c) const
{
	if(comp_global_convert_set) return comp_global_convert[cgl][cl][c];
	
	const auto &co = comp_gl[cgl]; 
	cgl += ((int)c-co.cla_comp[cl])*comp_mult[cl];
	if(comp_gl[cgl].cla_comp[cl] != c) emsg("Wrong c");
	
	return cgl;
}


/// Outputs a sampling error message
void Species::sampling_error(unsigned int trg, string warn) const
{
	const auto &tra = tra_gl[trg];
	run_error("Sampling error for transition '"+cla[tra.cl].tra[tra.tr].name+"': "+warn+". Transition distribution quantities (means, sds etc...) have threshold limits to ensure numerical accuracy. Consider changing the values or restricting the priors on model parameter which determine this transition."); 
}


/// Works out how to update equations using changes in populations
void Species::sim_linear_speedup_init(const vector <Equation> &eqn)
{
	if(!sim_linearise_speedup){
		switch(type){
		case INDIVIDUAL:
			for(auto m = 0u; m < markov_eqn.size(); m++){
				sim_linear_speedup.calc.push_back(m);
			}
			break;
			
		case POPULATION:
			for(auto tr = 0u; tr < tra_gl.size(); tr++){
				sim_linear_speedup.calc.push_back(tr);
			}
			break;
		}
		return;
	}
	
	vector <LinearFormInit> lfinit;
	
	switch(type){
	case INDIVIDUAL:
		for(auto m = 0u; m < markov_eqn.size(); m++){
			const auto &me = markov_eqn[m];
			auto e = me.eqn_ref;
			const auto &eq = eqn[e]; 
				if(eq.linearise.on && eq.linearise.pop_grad_time_dep == false){
				LinearFormInit lfi; lfi.m = m; lfi.e = e;
				lfinit.push_back(lfi);
			}
			else{
				sim_linear_speedup.calc.push_back(m);
			}
		}
		break;
		
	case POPULATION:
		for(auto tr = 0u; tr < tra_gl.size(); tr++){
			const auto &tra = tra_gl[tr];
			auto e = tra.dist_param[0].eq_ref;
			const auto &eq = eqn[e];
			
			if(eq.linearise.on && eq.linearise.pop_grad_time_dep == false){
				LinearFormInit lfi; lfi.m = tr; lfi.e = e;
				lfinit.push_back(lfi);
			}
			else{
				sim_linear_speedup.calc.push_back(tr);
			}
		}
		break;
	}
	
	set_linear_form(sim_linear_speedup.lin_form,lfinit,eqn);
}


/// Sets a linear form such that proposals can be made more quickly
void Species::set_linear_form(LinearForm &lin_form, const vector <LinearFormInit> &lfinit, const vector <Equation> &eqn) const
{
	lin_form.factor_nopop_only = false;
	
	Hash hash;
	for(auto k = 0u; k < lfinit.size(); k++){
		const auto &lfi = lfinit[k];
		
		auto e = lfi.e;
		const auto &eq = eqn[e];
		
		auto &lin = eq.linearise;
		
		// Generates a hash vector
		vector <unsigned int> vec;
		if(lin.pop_grad_precalc.size() != eq.pop_ref.size()) emsg("pop_ref size not right");
		
		for(auto i = 0u; i < eq.pop_ref.size(); i++){
			const auto &it = lin.pop_grad_precalc[i];
			vec.push_back(it.type);
			vec.push_back(it.num);
			vec.push_back(eq.pop_ref[i]);
		}
		
		auto n = hash.existing(vec);
		if(n == UNSET){
			n = lin_form.sum_e.size();
			hash.add(n,vec);
			lin_form.sum_e.push_back(e);
			
			for(auto i = 0u; i < eq.pop_ref.size(); i++){
				auto po = eq.pop_ref[i];
							
				{ // turn off
					auto k = 0u; while(k < lin_form.pop_affect.size() && lin_form.pop_affect[k].po != po) k++;
					if(k == lin_form.pop_affect.size()) k = UNSET;
					if(k != lin_form.hash_po.find(po)) emsg("lf prob");
				}
				
				auto k = lin_form.hash_po.find(po);
				if(k == UNSET){
					k = lin_form.pop_affect.size();
					PopAffect pa; pa.po = po;
					lin_form.hash_po.add(k,po);
					lin_form.pop_affect.push_back(pa);
				}
			
				GradRef gr; gr.ref = n; gr.index = i;
				lin_form.pop_affect[k].pop_grad_ref.push_back(gr);
			}
		}
		
		LinearFormItem lf;
		lf.m = lfi.m; lf.e = e; lf.sum_e_ref = n;
		lf.factor_precalc = eq.linearise.factor_precalc;
		lf.no_pop_precalc = eq.linearise.no_pop_precalc;

		lin_form.list.push_back(lf);
	}
	
	lin_form.nopop_same = true;
	lin_form.factor_same = true;
	
	const auto &list = lin_form.list;
	
	if(list.size() > 0){
		const auto &factor_precalc_fir = list[0].factor_precalc;
		const auto &no_pop_precalc_fir = list[0].no_pop_precalc;

		for(auto k = 1u; k < list.size(); k++){
			if(!item_equal(factor_precalc_fir,list[k].factor_precalc)){
				lin_form.factor_same = false;
			}				
			
			if(!item_equal(no_pop_precalc_fir,list[k].no_pop_precalc)){
				lin_form.nopop_same = false;
			}				
		}			
	}
	
	if(false){
		cout << lin_form.nopop_same << " nopop_same" << endl;
		cout << lin_form.factor_same << " factor same" << endl;

		for(const auto &lf : lin_form.list){
			cout << tra_gl[lf.m].name << ", ";
		}
		cout << "calculate fast" << endl;
		
		cout << lin_form.list.size() << " " << lin_form.sum_e.size() << "num" << endl;
		/*
		for(auto k = 0u; k < lin_form.pop_affect.size(); k++){
			const auto &pa = lin_form.pop_affect[k];
			cout << "  " << pa.po << ": ";
			for(const auto &pgr : pa.pop_grad_ref) cout << pgr.ref << " " << pgr.index << ", ";
			cout << "  grad ref" << endl;
		}
		*/
		emsg("mbp speed up");
	}
}


/// Determines if two equation items are equal
bool Species::item_equal(const EqItem &it1, const EqItem &it2) const
{
	if(it1.type != it2.type) return false;
	if(it1.num != it2.num) return false;
	return true;
}

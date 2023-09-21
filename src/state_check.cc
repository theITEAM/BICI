// This file file stores a state of the model

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

// Checks everything used in inference
void State::check(string ref)
{
	if(testing == false) return;
	
	check_timer -= clock();
	
	//cout << "CHECK\n";
	for(auto p = 0u; p < nspecies; p++){
	check_obs_like(p,ref);
	}
	
	//check_cpop_st(ref);
	check_popnum_t(ref);
		
	for(auto p = 0u; p < nspecies; p++){
		const auto &sp = model.species[p];
		switch(sp.type){
		case INDIVIDUAL:
			check_ref(p,ref);
			check_ie(p,ref);
			check_markov_trans(p,ref);
			check_nm_trans(p,ref);
			break;
			
		case POPULATION:
			check_pop_like(p,ref);
			break;
		}
		
		check_obs_like(p,ref);
	}

	check_like(ref);
	check_spline(ref);
	check_prior(ref);
	
	check_maps(ref);

	check_timer += clock();
}


// Checks references between markov_eqn_vari,   ind.ev  and   nm_trans
void State::check_ref(unsigned int p, string ref)
{
	const auto &sp = model.species[p];
	auto &ssp = species[p];
	for(auto i = 0u; i < ssp.individual.size(); i++){
		auto &ind = ssp.individual[i];
		ssp.check_cl_ev(ind);
	}
	
	// Checks individual ev consistent
	for(auto i = 0u; i < ssp.individual.size(); i++){
		auto &ind = ssp.individual[i];
		auto c = ind.cinit;
		for(auto k = 0u; k < ind.ev.size(); k++){
			auto &ev = ind.ev[k];
			
			switch(ev.type){
			case ENTER_EV: 
				if(c != UNSET) emsg("SHould start unset");
				c = ev.c_gl;
				break;
			
			case LEAVE_EV:
				c = UNSET;
				break;
				
			case MOVE_EV:
				c = sp.update_c_comp(c,ev.tr_gl,ev.c_gl);
				break;
				
			case NM_TRANS_EV: 
				// TO DO
				break;
		
			case M_TRANS_EV:	
				const auto &tra = sp.tra_gl[ev.tr_gl];
				if(tra.i != c) emsg("hh");
				if(tra.f != ind.ev[k].c_after) emsg("JJ");
				c = tra.f;
				break;
			}
		}
	}
	

	// Checks individual ev correctly time ordered
	for(auto i = 0u; i < ssp.individual.size(); i++){
		auto &ind = ssp.individual[i];

		if(ind.ev.size() > 0){
			for(auto k = 0u; k < ind.ev.size()-1; k++){
				if(ind.ev[k].t > ind.ev[k+1].t) emsg("time ordering problem"+ref);
			}
		}
	}
		
	// Checks that ind cl_ev is correctly specified
	for(auto i = 0u; i < ssp.individual.size(); i++){
		auto &ind = ssp.individual[i];
		for(auto cl = 0u; cl < sp.ncla; cl++){
			for(auto k = 1u; k < ind.cl_ev[cl].size(); k++){
				if(ind.ev[ind.cl_ev[cl][k-1].index_ev].t >= ind.ev[ind.cl_ev[cl][k].index_ev].t){
					emsg("cl order wrong");
				}
			}
			
			for(auto k = 0u; k < ind.cl_ev[cl].size(); k++){
				const auto &cl_ev =ind.cl_ev[cl][k];
			
				auto &ev = ind.ev[cl_ev.index_ev];
				
				if(cl_ev.nm_ref != UNSET){
					const auto &trang = sp.tra_gl[ev.tr_gl];			
					if(cl_ev.nm_ref != trang.nm_ref) emsg("nm_ref problem");
	
					auto &nti = ssp.nm_trans_ind[cl_ev.nm_ref][cl_ev.index_nmtrans];
					if(nti.t != ev.t) emsg("DOes not match");		
					if(nti.i != i) emsg("Does not match");
					if(nti.cl != cl) emsg("Does not match");
					if(nti.index != k) emsg("Does not match");
				}
			}
		}
	}
	
	for(auto j = 0u; j < sp.nm_trans.size(); j++){
		const auto &nti = ssp.nm_trans_ind[j];
		for(auto k = 0u; k < nti.size(); k++){
			auto &nt = nti[k];
			auto &cl_ev = ssp.individual[nt.i].cl_ev[nt.cl][nt.index];
			
			if(cl_ev.index_nmtrans != k) emsg("Match");
		}
	}
}


// Checks that the likelihood is correctly specified
void State::check_markov_trans(unsigned int p, string ref)
{
	const auto &sp = model.species[p]; 
	auto &ssp = species[p];

	const auto individual_store = ssp.individual;

	// Checks that exp_fe is set correctly
	for(auto f = 0u; f < sp.fix_effect.size(); f++) ssp.set_exp_fe(f);
	
	for(auto i = 0u; i < ssp.individual.size(); i++){
		const auto &ind = ssp.individual[i];
		const auto &ind_st = individual_store[i];
		
		for(auto f = 0u; f < sp.fix_effect.size(); f++){
			if(dif(ind.exp_fe[f],ind_st.exp_fe[f])) emsg("Error in exp_fe");
		}
	}
	
	// Checks that exp_ie is calculated correctly
	for(auto i = 0u; i < ssp.individual.size(); i++){
		auto &ind = ssp.individual[i];
		const auto ind_store = ind;
		ssp.set_exp_ie(ind);
		
		for(auto j = 0u; j < sp.ind_effect.size(); j++){
			if(dif(ind.exp_ie[j],ind_store.exp_ie[j])) emsg("Error in exp_ie: "+ref);
		}
	}
	
	auto markov_eqn_vari_store = ssp.markov_eqn_vari;
	ssp.likelihood_indfac_int();
	ssp.likelihood_ind_trans();
	
	// Checks that markov equations are correctly specified
	for(auto e = 0u; e < ssp.N; e++){
		auto Li_markov_store = ssp.Li_markov[e];
		ssp.likelihood_markov_value(e,empty,ALL_DIV,popnum_t);
		auto temp = 0.0;
		ssp.likelihood_markov(e,empty,ALL_DIV,temp);
			
		const auto &divi = ssp.markov_eqn_vari[e].div;
		for(auto ti = 0u; ti < divi.size(); ti++){
			auto &div = divi[ti];
			const auto &divc = markov_eqn_vari_store[e].div[ti];
				
			if(dif(div.value,divc.value)){
				cout << ti <<  "should be: " << div.value << "  but is: " << divc.value << "val\n";
				emsg("value error");
			}
			
			if(dif(div.indfac_int,divc.indfac_int)){
				cout << ti <<  "should be: " << div.indfac_int << "  but is: " << divc.indfac_int << "val\n";
				emsg("indfac_int");
			}
			
			if(div.ind_trans.size() != divc.ind_trans.size()) emsg("Wrong si");
			
			for(auto j = 0u; j < div.ind_trans.size(); j++){
				auto k = 0u; while(k < divc.ind_trans.size() && !(div.ind_trans[j].i == divc.ind_trans[k].i && div.ind_trans[j].index == divc.ind_trans[k].index)) k++;
				
				if(k == divc.ind_trans.size()) emsg("Cannot find");
			}
			
			if(dif(ssp.Li_markov[e][ti],Li_markov_store[ti])){
				emsg("Li_markov error");
			}
		}
	}
	
	// Checks references between markov eqn div   and individual events
	for(auto e = 0u; e < ssp.N; e++){
		auto &mev = ssp.markov_eqn_vari[e];
		for(auto ti = 0u; ti < mev.div.size(); ti++){
			auto &div = mev.div[ti];
			for(auto j = 0u; j < div.ind_trans.size(); j++){
				auto &dir = div.ind_trans[j];
				auto &ind = ssp.individual[dir.i];
				auto &ev = ind.ev[dir.index];
				
				auto &mef = ev.markov_eq_ref;
				if(mef.e != e) emsg("ref prob");
				if(mef.ti != ti) emsg("ref prob");
				if(mef.index != j) emsg("ref prob");
			}
		}
	}
	
	for(auto i = 0u; i < ssp.individual.size(); i++){
		auto &ind = ssp.individual[i];
		for(auto k = 0u; k < ind.ev.size(); k++){
			auto &ev = ind.ev[k];
			if(ev.type == M_TRANS_EV){
				auto &mer = ev.markov_eq_ref;
				if(mer.e == UNSET) emsg("prob");
				auto &it = ssp.markov_eqn_vari[mer.e].div[mer.ti].ind_trans[mer.index];
				if(it.i != i) emsg("prob");
				if(it.index != k) emsg("prob");
			}
		}
	}
}


/// Checks likelihoods for non-Markovian transitions
void State::check_nm_trans(unsigned int p, string ref)
{
	auto &ssp = species[p];
	for(auto i = 0u; i < ssp.nm_trans_ind.size(); i++){
		auto Li_nm_trans_store = ssp.Li_nm_trans[i];
		
		double temp;
		ssp.likelihood_nm_trans(i,popnum_t,temp);
		
		if(dif(ssp.Li_nm_trans[i],Li_nm_trans_store)) emsg("nm trans error"+ref);
	}
}


// Checks that the quantities in like are correctly specified
void State::check_like(string ref)
{
	auto like_st = like;
	
	calculate_like();
	
	if(dif(like.obs,like_st.obs)) emsg("like obs");
		
	if(dif(like.markov,like_st.markov)){
		cout << like.markov << " " << like_st.markov << "Pp\n";
		emsg("like markov "+ref);
	}
	
	if(dif(like.nm_trans,like_st.nm_trans)) emsg("like nm_trans");
	
	if(dif(like.prior,like_st.prior)) emsg("like prior");
		
	if(dif(like.ie,like_st.ie)) emsg("like ie");
	
	if(dif(like.dist,like_st.dist)) emsg("like dist");
	
	if(dif(like.spline_prior,like_st.spline_prior)) emsg("like spline prior "+ref);
}


// Checks that the splines are correctly specified
void State::check_spline(string ref)
{	
	auto spline_val_st = spline_val;
	if(spline_val.size() != model.spline.size()) emsg("Wrong size"+ref);
			
	spline_init();

	for(auto j = 0u; j < model.spline.size(); j++){
		const auto &val = spline_val[j].val;
		const auto &val_st = spline_val_st[j].val;
		for(auto ti = 0u; ti < T; ti++){
			if(val[ti] != val_st[ti]){
				cout << val[ti] << " " << val_st[ti] << " jj\n";
				emsg("Spline error");
			}
		}
	}
}


/// Checks that priors are correctly specified
void State::check_prior(string ref)
{
	// Checks prior_prob is correct
	auto prior_prob_st = prior_prob;
	prior_prob = model.prior_prob(param_val);
	for(auto th = 0u; th < prior_prob.size(); th++){
		if(dif(prior_prob[th],prior_prob_st[th])){
			cout << prior_prob[th] << " " << prior_prob_st[th] << " jj\n";
			emsg("Prior probability wrong"+ref);
		}
	}
	
	// Checks dist_prob is correct
	auto dist_prob_st = dist_prob;
	dist_prob = model.dist_prob(param_val);
	for(auto th = 0u; th < dist_prob.size(); th++){
		if(dif(dist_prob[th],dist_prob_st[th])) emsg("Dist prob probability wrong");
	}
	
	// Checks spline_prior is correct
	auto spline_prior_st = spline_prior;
	spline_prior = model.spline_prior(param_val);
	for(auto s = 0u; s < spline_prior.size(); s++){
		if(dif(spline_prior[s],spline_prior_st[s])) emsg("Spline Prior probability wrong");
	}
}


/// Checks that popnum_t is correctly specified
void State::check_popnum_t(string ref)
{
	auto popnum_t_store = popnum_t;
	
	popnum_t = model.calculate_popnum_t(species);
	
	if(popnum_t.size() != T) emsg("Wrong size"+ref);
	
	for(auto ti = 0u; ti < T; ti++){
		if(popnum_t[ti].size() != model.pop.size()) emsg("Wrong size");
		for(auto k = 0u; k < model.pop.size(); k++){
			if(dif(popnum_t[ti][k],popnum_t_store[ti][k])){
				cout << popnum_t[ti][k] << " " << popnum_t_store[ti][k] << "\n";
				cout << "ti: " << ti << " " << k << "\n";
				emsg("popnum_t problem"+ref);
			}
		}
	}
}
	

/// Checks that cpop_st is correctly set	
void State::check_cpop_st(string ref)
{
	for(auto p = 0u; p < model.nspecies; p++){
		auto &ssp = species[p];
		ssp.cpop = ssp.cpop_init;
	}
	
	auto popnum = popnum_t[0];
	
	for(auto ti = 0u; ti < T; ti++){
		for(auto po = 0u; po < popnum.size(); po++){
			if(popnum[po] != popnum_t[ti][po]){
				emsg("Proplem with popnum"+ref);
			}
		}
		
		for(auto p = 0u; p < model.nspecies; p++){
			auto &sp = model.species[p];
			auto &ssp = species[p];
			
			auto C = sp.comp_gl.size();
			auto N = sp.tra_gl.size();
	
			switch(ssp.type){
			case INDIVIDUAL:
				{	auto cpop_st_store = ssp.cpop_st;
		
					emsg("to do");
				}
				break;
			
			case POPULATION:
				{
					for(auto c = 0u; c < C; c++){
						if(ssp.cpop_st[ti][c] != ssp.cpop[c]) emsg("error");
					}
					
					for(auto tr = 0u; tr < N; tr++){
						auto n = ssp.trans_num[tr][ti];
					
						if(n != 0){
							const auto &tra = sp.tra_gl[tr];
							auto ci = tra.i;
							if(ci != SOURCE){
								ssp.cpop[ci] -= n;
								for(auto po : sp.comp_gl[ci].pop_ref_simp) popnum[po] -= n;
							}
							
							auto cf = tra.f;
							if(cf != SINK){
								ssp.cpop[cf] += n;
								for(auto po : sp.comp_gl[cf].pop_ref_simp) popnum[po] += n;
							}
						}
					}
				}
				break;
			}
		}
	}
}


/// Checks likelihoods for individual effect
void State::check_ie(unsigned int p, string ref)
{
	auto &sp = model.species[p];
	auto &ssp = species[p];

	for(auto g = 0u; g < sp.ind_eff_group.size(); g++){
		const auto &iesg = ssp.ind_eff_group_sampler[g];
		auto iesg_store = iesg;
		
		ssp.calculate_omega(g);
		
		if(dif(iesg_store.omega, iesg.omega)){
			print("om store",iesg_store.omega); print("om2",iesg.omega);
			emsg("Omega problem");
		}
		if(dif(iesg_store.omega_Z, iesg.omega_Z)) emsg("Omega_Z problem");
		if(model.mode == INF){
			if(dif(iesg_store.omega_inv,iesg.omega_inv)) emsg("Omega_inv problem");
		}	
		
		auto Li_ie_store = ssp.Li_ie[g];
		auto temp = 0.0;
		ssp.likelihood_indeff_group(g,temp);
	
		if(dif(Li_ie_store,ssp.Li_ie[g])){
			cout << Li_ie_store - ssp.Li_ie[g] << " dif\n";
			emsg("L_ie problem"+ref);
		}
	}
}


// Checks the population likelihood is correct
void State::check_pop_like(unsigned int p, string ref)
{
	auto &sp = model.species[p];
	auto &ssp = species[p];
	
	auto Li_markov_pop_store = ssp.Li_markov_pop;
	auto trans_num_store = ssp.trans_num;
	auto tnum_mean_st_st = ssp.tnum_mean_st;
	
	ssp.likelihood_pop(popnum_t);
	
	for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
		for(auto ti = 0u; ti < T; ti++){
			if(dif(ssp.tnum_mean_st[tr][ti],tnum_mean_st_st[tr][ti])){
				cout << p << " " <<tr << " " << ti << "  true:"<<  ssp.tnum_mean_st[tr][ti] << " " << tnum_mean_st_st[tr][ti] << " jj\n";
				emsg("trans num mean prob"+ref);
			}			
			
			if(dif(ssp.trans_num[tr][ti],trans_num_store[tr][ti])) emsg("trans num prob"+ref);
			
			if(dif(ssp.Li_markov_pop[tr][ti],Li_markov_pop_store[tr][ti])){
				cout << p << " " <<tr << " " << ti << "  true:"<<  ssp.Li_markov_pop[tr][ti] << " " << Li_markov_pop_store[tr][ti] << " jj\n";
				
				emsg("Li mark prob");
			}
		}
	}
}


/// Checks all the maps are correcly set
void State::check_maps(string ref)
{
	for(auto p = 0u; p < species.size(); p++){
		if(pop_list.size() != 0) emsg("pop list"+ref);
		if(markov_eqn_list.size() != 0) emsg("pop list");
		if(trans_list.size() != 0) emsg("pop list");
		
		for(auto i = 0u; i < model.pop.size(); i++){
			if(pop_map[i] != UNSET) emsg("prob");
		}
	
		switch(model.species[p].type){
		case INDIVIDUAL: 
			for(auto i = 0u; i < model.species[p].markov_eqn.size(); i++){
				if(markov_eqn_map[p][i] != false) emsg("prob");
			}
			break;
		
		case POPULATION: 
			for(auto i = 0u; i < model.species[p].tra_gl.size(); i++){
				if(trans_map[p][i] != false) emsg("prob");
			}
			break;
		}
	}
}


/// Checks that the observation likelihood is correctly specified
void State::check_obs_like(unsigned int p, string ref)
{
	auto &sp = model.species[p];
	auto &ssp = species[p];

	auto pop_data_num_store = ssp.pop_data_num;
	auto pop_trans_data_num_store = ssp.pop_trans_data_num;
	
	auto Li_obs_store = ssp.Li_obs;
	
	auto temp = 0.0;
	ssp.likelihood_observation(temp);
	
	if(dif(Li_obs_store,ssp.Li_obs)){
		cout << Li_obs_store << "   true: " << ssp.Li_obs << "j\n";
		emsg("Li_obs error"+ref);
	}
	
	for(auto i = 0u; i < sp.pop_data.size(); i++){
		if(dif(pop_data_num_store[i],ssp.pop_data_num[i])){
			cout << i << " i   " << pop_data_num_store[i] << " " << ssp.pop_data_num[i] << "uu\n";
			emsg("pop_data_num error");
		}
	}
	
	for(auto i = 0u; i < sp.pop_trans_data.size(); i++){
		if(dif(pop_trans_data_num_store[i],ssp.pop_trans_data_num[i])) emsg("pop_trans_data_num error");
	}
}

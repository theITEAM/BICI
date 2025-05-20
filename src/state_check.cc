// This file checks the state is correctly specified (it looks for any inconsistencies)

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
using namespace std;

#include "state.hh"
#include "utils.hh"

/// Checks for any transition numbers becoming negative
void State::check_trans_num(string ref)
{
	for(auto p = 0u; p < model.species.size(); p++){	
		const auto &sp = model.species[p]; 
		auto &ssp = species[p];
		if(sp.type == POPULATION){
			for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
				for(auto ti = 0u; ti < T; ti++){
					if(ssp.trans_num[tr][ti] < 0) emsg("trans_num has become negative"+ref);
					//if(trans_num_st[tr][ti] < 0) emsg("trans_num has become negative"+ref);
				}
			}
		}
	}
}


/// Performs a simple check
void State::check_simp(string ref)
{
	for(auto p = 0u; p < model.species.size(); p++){	
		const auto &sp = model.species[p]; 
		auto &ssp = species[p];

		const auto individual_store = ssp.individual;

		// Checks that ev.m is correctly specified
		for(auto i = 0u; i < ssp.individual.size(); i++){
			const auto &ind = ssp.individual[i];
			for(auto e = 0u; e < ind.ev.size(); e++){
				const auto &ev = ind.ev[e];
				switch(ev.type){
				case M_TRANS_EV:
					{
						const auto &tra = sp.tra_gl[ev.tr_gl];
						if(ev.m != tra.markov_eqn_ref) emsg("Markov equation not correct"+ref);
					
						const auto &mev = ssp.markov_eqn_vari[tra.markov_eqn_ref];
						if(mev.time_vari){
							if(ev.ti != ssp.get_ti(ev.t)) emsg("ti wrong1");
						}
						else{
							if(ev.ti != 0) emsg("ti wrong2");
						}
						
						if(mev.div[ev.ti].ind_trans[ev.index].i != i){
							emsg("i is wrong");
						}
					}
					break;
				
				case NM_TRANS_EV:
					{
						auto e_bef = get_event_before(e,ind.ev);
						if(ev.e_origin != e_bef){
							emsg("e_bef not right");
						}
						
						const auto &tra = sp.tra_gl[ev.tr_gl];
						
						const auto &ev_orig = ind.ev[e_bef];
	
						if(ev.m != sp.get_tra_m(tra,ev_orig)){
							emsg("NMarkov equation not correct2"+ref);
						}
						
						if(ev.ti != ssp.get_ti(ev_orig.t)){
							emsg("ti wrong3"+ref);
						}
					
						if(ssp.nm_trans_ev_ref[ev.m][ev.ti][ev.index].i != i){
							emsg("i is wrong nm");
						}
						
						if(ssp.nm_trans_ev_ref[ev.m][ev.ti][ev.index].e_end != e){
							emsg("e_end is wrong");
						}
					}
					break;
				default: break;
				}
			}
		}

		for(auto e = 0u; e < ssp.markov_eqn_vari.size(); e++){
			const auto &mev = ssp.markov_eqn_vari[e];
			for(auto ti = 0u; ti < mev.div.size(); ti++){
				const auto &dit = mev.div[ti].ind_trans;
				for(auto k = 0u; k < dit.size(); k++){ 
					const auto &it = dit[k];
					const auto &ev = ssp.individual[it.i].ev[it.index];	
					if(ev.index != k){
						//cout << it.i << " "<< it.index " " << " " 
						emsg("k wrong");
					}
				}
			}
		}
	}
}


// Checks everything used in inference
void State::check(string ref)
{
	if(testing == false) return;
	
	timer[CHECK_TIMER] -= clock();

	check_trans_num(ref);

	check_event_boundary(ref);

	check_effect_out_of_range();

	check_dependent_param(ref);
	
	check_simp(ref);

	check_popnum_ind(ref);
	
	check_genetic_value(ref);
	
	check_spline_store(ref);
	
	check_cpop_st(ref);
	
	check_popnum_t(ref);
	
	for(auto p = 0u; p < nspecies; p++){
		const auto &sp = model.species[p];
		
		switch(sp.type){
		case INDIVIDUAL:
			check_event_observed(p,ref);
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
		check_init_cond_like(p,ref);
	}
	check_init_cond_prior(ref);
	
	check_like(ref);
	
	check_spline(ref);
	
	check_prior(ref);
	
	check_maps(ref);

	check_add_move_rem(ref);
	
	if(check_lin) check_linearise();
	
	timer[CHECK_TIMER] += clock();
}


/// Checks that dependent parameters are correctly set
void State::check_dependent_param(string ref)
{
	for(auto th = 0u; th < model.param_vec.size(); th++){	
		const auto &pv = model.param_vec[th];
		if(pv.variety == REPARAM_PARAM){
			auto re = model.param[pv.th].value[pv.index].eq_ref;
			if(re == UNSET) emsg("Reparam is not set");	
			auto value = model.eqn[re].calculate_param_only(param_val);
		
			if(dif(value,param_val[th],dif_thresh)){
				cout << value << " " << param_val[th] << "compare" << endl;
				emsg("Parameter value different "+ref);
			}
		}
	}
}


/// Checks references between markov_eqn_vari, ind.ev and nm_trans
void State::check_ref(unsigned int p, string refst)
{
	const auto &sp = model.species[p];
	auto &ssp = species[p];
	
	// Checks that event references are all correct
	for(auto i = 0u; i < ssp.individual.size(); i++){
		const auto &ind = ssp.individual[i];
		for(auto e = 0u; e < ind.ev.size(); e++){
			const auto &ev = ind.ev[e]; 
			
			if(ev.type != NM_TRANS_EV){
			}
			else{	
				if(ev.e_origin != get_event_before(e,ind.ev)){
					emsg("e_origin is wrong");
				}
				
				const auto &ev_orig = ind.ev[ev.e_origin];
				const auto &tra = sp.tra_gl[ev.tr_gl];
				auto m = sp.get_tra_m(tra,ev_orig);
	
				if(m != ev.m){
					emsg("m is not right");
				}
				
				auto t = ind.ev[ev.e_origin].t;		
				auto ti = sp.get_ti(t);
				if(ti != ev.ti){
					emsg("ti is not right");
				}
			}
		
			const auto nmt_ind = ev.index;
			
			if(ev.type == NM_TRANS_EV){
				const auto &ev_orig = ind.ev[ev.e_origin];
				
				auto t = ev_orig.t;
			
				auto ti = sp.get_ti(t);
		
				const auto &tra = sp.tra_gl[ev.tr_gl];
				auto m = tra.nm_trans_ref;
				
				if(tra.i != ev_orig.c_after){               // Accounts for intermdiate transition
					m = sp.tra_gl[tra.transform[ev_orig.c_after]].nm_trans_ref;
				}

				if(nmt_ind == UNSET) emsg("index should be set");
			
				if(nmt_ind >= ssp.nm_trans_ev_ref[m][ti].size()){
					emsg("index out of range");
				}
				
				const auto &ev_ref = ssp.nm_trans_ev_ref[m][ti][nmt_ind];
				
				if(ev_ref.i != i){
					emsg("Ref i prob");
				}
			
				if(ev_ref.e_end != e){
					emsg("Ref e_end prob   "+refst);
				}
				
				const auto &nmt = sp.nm_trans[m];
				const auto &ref = nmt.dist_param_eq_ref;
			
				vector <double> ref_val(ref.size());
				for(auto i = 0u; i < ref.size(); i++){
					ref_val[i] = model.eqn[ref[i]].calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
				}
				
				auto dt = ev.t - t;				
				
				if(dif(ev.Li,nm_trans_like(nmt.type,dt,ref_val),dif_thresh)){
					emsg("nm like not agree"+refst);
				}
				
				auto bp_eq = nmt.bp_eq;
				if(bp_eq != UNSET){
					double bp_val;
					if(bp_eq == BP_FROM_OTHERS){
						bp_val = 1.0;
						for(auto e : nmt.bp_other_eq){
							bp_val -=  model.eqn[e].calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
						}
					}
					else{
						bp_val = model.eqn[bp_eq].calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
						
						if(nmt.all_branches){
							auto div = 0.0;
							for(auto e : nmt.bp_all_eq){
								div += model.eqn[e].calculate_indfac(ind,ti,popnum_t[ti],param_val,spline_val);
							}
							bp_val /= div;
						}
					}
					
					if(dif(ev.Li_bp,get_log_zero_one(bp_val),dif_thresh)){
						emsg("Libp like not agree");
					}
				}
				else{
					if(ev.Li_bp != 0) emsg("Li_bp should be zero"+refst);
				}
			}
		}
	}
	
	// Checks incomp_ref
	{
		for(auto i = 0u; i < ssp.individual.size(); i++){
			const auto &ind = ssp.individual[i];
			for(auto cl = 0u; cl < sp.ncla; cl++){
				const auto &inmtr = ind.incomp_ref[cl];
				if(inmtr.on){
					const auto &nmt = sp.nm_trans_incomp[inmtr.n];
					auto ii = ssp.nm_trans_incomp_ref[inmtr.n][inmtr.ti][inmtr.index];
	
					if(ii != i){
						emsg("i does not agree "+refst);
					}
				
					auto dt = inmtr.t_end - ind.ev[inmtr.e_begin].t; if(dt <= 0) emsg("zero time3");
				
					auto val = ssp.get_nm_incomp_val(nmt,inmtr.ti,ind,popnum_t);

					auto Li_new = ssp.nm_trans_incomp_full_like(nmt.nmtrans_ref,dt,val.ref_val,val.bp_val);
					if(dif(Li_new,inmtr.Li,dif_thresh)){
						emsg("Li nm_trans_incomp does not agree");
					}			
				}
			}
		}
		
		for(auto i = 0u; i < ssp.individual.size(); i++){
			auto &ind = ssp.individual[i];
			auto inmtr_store = ind.incomp_ref;
			double temp;
			ssp.set_incomp_ref(i,popnum_t,temp);
			
			for(auto cl = 0u; cl < sp.ncla; cl++){
				auto it1 = inmtr_store[cl];
				auto it2 = ind.incomp_ref[cl];
				if(it1.on != it2.on) emsg("On not right");
				if(it1.on){
					if(it1.n != it2.n) emsg("n not right");
					if(it1.ti != it2.ti){
						emsg("ti not right");
					}
					if(dif(it1.dt,it2.dt,dif_thresh)){
						emsg("dt not right");
					}
					
					if(it1.index != it2.index){
						emsg("index not right");
					}
					
					if(dif(it1.Li,it2.Li,dif_thresh)){
						emsg("Li not right");
					}
				}
			}
		}
	}
	
	// Checks individual ev consistent
	for(auto i = 0u; i < ssp.individual.size(); i++){
		auto &ind = ssp.individual[i];
			
		if(i >= sp.nindividual_in && i < sp.nindividual_obs){  // Individuals only in A matrix
			if(ind.ev.size() != 0) emsg("Individual should have no events");
		}
		else{
			if(ind.ev.size() == 0){
				cout << ind.name << endl;
				emsg("Individual has no events "+refst);
			}
			
			auto c = UNSET;
			for(auto k = 0u; k < ind.ev.size(); k++){
				auto &ev = ind.ev[k];
				
				switch(ev.type){
				case ENTER_EV: 
					if(c != UNSET) emsg("SHould start unset");
					c = ev.c_after;
					if(ev.cl != UNSET) emsg("ev cl should be unset");
					break;
				
				case LEAVE_EV:
					c = UNSET;
					if(ev.cl != UNSET) emsg("ev cl should be unset");
					break;
					
				case MOVE_EV:
					c = sp.update_c_comp(c,ev.cl,ev.move_c);
					break;
			
				case M_TRANS_EV: case NM_TRANS_EV: 
					{
						const auto &tra = sp.tra_gl[ev.tr_gl];
						if(ev.cl != tra.cl){
							emsg("ev cl should be set1");
						}				
						if(tra.i != c){
							emsg("Individual event sequence problem");
						}
						if(tra.f != ind.ev[k].c_after) emsg("JJ");
						c = tra.f;
					}
					break;
				}
				
				if(c != ev.c_after) emsg("SHould be c_after");
			}
		}
	}
	
	// Checks individual ev correctly time ordered
	for(auto i = 0u; i < ssp.individual.size(); i++){
		auto &ind = ssp.individual[i];

		if(ind.ev.size() > 0){
			for(auto k = 0u; k < ind.ev.size()-1; k++){
				if(ind.ev[k].t > ind.ev[k+1].t) emsg("time ordering problem"+refst);
			}
		}
	}
}


// Checks the likelihood is correctly specified
void State::check_markov_trans(unsigned int p, string ref)
{
	const auto &sp = model.species[p]; 
	auto &ssp = species[p];

	const auto individual_store = ssp.individual;

	// Checks that ev.m is correctly specified
	for(auto i = 0u; i < ssp.individual.size(); i++){
		const auto &ind = ssp.individual[i];
		for(auto e = 0u; e < ind.ev.size(); e++){
			const auto &ev = ind.ev[e];
			switch(ev.type){
			case M_TRANS_EV:
				{
					const auto &tra = sp.tra_gl[ev.tr_gl];
					if(ev.m != tra.markov_eqn_ref) emsg("Markov equation not correct");
				}
				break;
			
			case NM_TRANS_EV:
				{
					const auto &tra = sp.tra_gl[ev.tr_gl];
					auto e_bef = get_event_before(e,ind.ev);
					if(ev.e_origin != e_bef){
						emsg("e_bef not right");
					}
					
					const auto &ev_orig = ind.ev[e_bef];
	
					if(ev.m != sp.get_tra_m(tra,ev_orig)){
						emsg("NMarkov equation not correct"+ref);
					}
				}
				break;
				
			default: break;
			}
		}
	}
	
	// Checks that exp_fe is set correctly
	for(auto f = 0u; f < sp.fix_effect.size(); f++) ssp.set_exp_fe(f);
	
	for(auto i = 0u; i < ssp.individual.size(); i++){
		const auto &ind = ssp.individual[i];
		const auto &ind_st = individual_store[i];
		
		for(auto f = 0u; f < sp.fix_effect.size(); f++){
			if(dif(ind.exp_fe[f],ind_st.exp_fe[f],dif_thresh)){
				emsg("Error in exp_fe");
			}
		}
	}
	
	// Checks that exp_ie is calculated correctly
	for(auto i = 0u; i < ssp.individual.size(); i++){
		auto &ind = ssp.individual[i];
		const auto ind_store = ind;
		ssp.set_exp_ie(ind);
		
		for(auto j = 0u; j < sp.ind_effect.size(); j++){
			if(dif(ind.exp_ie[j],ind_store.exp_ie[j],dif_thresh)){
				emsg("Error in exp_ie: "+ref);
			}
		}
	}
	
	auto markov_eqn_vari_store = ssp.markov_eqn_vari;
	ssp.likelihood_indfac_int();
	ssp.likelihood_ind_trans();
	
	// Checks that markov equations are correctly specified
	for(auto e = 0u; e < ssp.N; e++){
		auto Li_markov_store = ssp.Li_markov[e];
		
		auto list = seq_vec(ssp.markov_eqn_vari[e].div.size());
		ssp.likelihood_markov_value(e,list,popnum_t);
		auto temp = 0.0;
		ssp.likelihood_markov(e,list,temp);
			
		const auto &divi = ssp.markov_eqn_vari[e].div;
		for(auto ti = 0u; ti < divi.size(); ti++){
			auto &div = divi[ti];
			const auto &divc = markov_eqn_vari_store[e].div[ti];
				
			if(dif(div.value,divc.value,dif_thresh)){
				emsg("value error");
			}
			
			if(dif(div.indfac_int,divc.indfac_int,dif_thresh)){
				emsg("indfac_int");
			}
			
			if(div.ind_trans.size() != divc.ind_trans.size()){
				emsg("Wrong siz");
			}
			
			for(auto j = 0u; j < div.ind_trans.size(); j++){
				auto k = 0u; while(k < divc.ind_trans.size() && !(div.ind_trans[j].i == divc.ind_trans[k].i && div.ind_trans[j].index == divc.ind_trans[k].index)) k++;
				
				if(k == divc.ind_trans.size()){
					emsg("Cannot find here");
				}
			}
			
			if(dif(ssp.Li_markov[e][ti],Li_markov_store[ti],dif_thresh)){
				emsg("Li_markov error"+ref);
			}
		}
	}
	
	// Checks references between markov eqn div and individual events
	for(auto e = 0u; e < ssp.N; e++){
		auto &mev = ssp.markov_eqn_vari[e];
		for(auto ti = 0u; ti < mev.div.size(); ti++){
			auto &div = mev.div[ti];
			for(auto j = 0u; j < div.ind_trans.size(); j++){
				auto &dir = div.ind_trans[j];
				auto &ind = ssp.individual[dir.i];
				auto &ev = ind.ev[dir.index];
				
				if(ev.m != e) emsg("ref prob");
				if(ev.ti != ti) emsg("ref prob");
				if(ev.index != j) emsg("ref prob");
			}
		}
	}
	
	for(auto i = 0u; i < ssp.individual.size(); i++){
		auto &ind = ssp.individual[i];
		for(auto k = 0u; k < ind.ev.size(); k++){
			auto &ev = ind.ev[k];
			if(ev.type == M_TRANS_EV){
				if(ev.m == UNSET) emsg("prob");
				auto &it = ssp.markov_eqn_vari[ev.m].div[ev.ti].ind_trans[ev.index];
				if(it.i != i) emsg("prob");
				if(it.index != k) emsg("prob");
			}
		}
	}
}


/// Checks likelihoods for non-Markovian transitions
void State::check_nm_trans(unsigned int p, string ref)
{
	const auto &sp = model.species[p];
	
	auto all_vec = seq_vec(T);
	
	auto &ssp = species[p];
	
	// Checks likelihoods
	
	auto ind_store = ssp.individual;
	for(auto i = 0u; i < sp.nm_trans.size(); i++){
		double temp;
		ssp.likelihood_nm_trans(i,all_vec,popnum_t,temp);
	}
	
	for(auto i = 0u; i < ind_store.size(); i++){
		for(auto e = 0u; e < ind_store[i].ev.size(); e++){
			const auto &ev_store = ind_store[i].ev[e];	
			const auto &ev = ssp.individual[i].ev[e];
			if(ev.type == NM_TRANS_EV){
				if(dif(ev_store.Li,ev.Li,dif_thresh)){
					emsg("ev Li problem");
				}
				if(dif(ev_store.Li_bp,ev.Li_bp,dif_thresh)){
					emsg("ev Li_bp problem");
				}
			}
		}
	}
		
	for(auto i = 0u; i < sp.nm_trans_incomp.size(); i++){
		double temp;
		ssp.likelihood_nm_trans_incomp(i,all_vec,popnum_t,temp);
	}
	
	for(auto i = 0u; i < ind_store.size(); i++){
		for(auto e = 0u; e < ind_store[i].incomp_ref.size(); e++){
			const auto &inmtr_store = ind_store[i].incomp_ref[e];	
			const auto &inmtr = ssp.individual[i].incomp_ref[e];
			if(inmtr.on){
				if(dif(inmtr_store.Li,inmtr.Li,dif_thresh)){
					emsg("nm_ref Li problem"+ref);
				}
			}
		}
	}
}


// Checks the quantities in "like" are correctly specified
void State::check_like(string ref)
{
	if(std::isnan(like.init_cond)) emsg("init_cond nan");
	if(std::isnan(like.init_cond_prior)) emsg("init_cond_prior nan");
	if(std::isnan(like.obs)) emsg("obs nan");
	if(std::isnan(like.prior)) emsg("prior nan");
	if(std::isnan(like.spline_prior)) emsg("spline_prior nan");
	if(std::isnan(like.dist)) emsg("dist nan");
	if(std::isnan(like.markov)) emsg("markov nan");
	if(std::isnan(like.nm_trans)) emsg("nm_trans nan");
	if(std::isnan(like.genetic_process)) emsg("genetic_process nan");
	if(std::isnan(like.genetic_obs)) emsg("genetic_obs nan");
	if(std::isnan(like.ie)) emsg("ie nan");
	
	auto like_st = like;
	
	calculate_like();

	if(dif(like.init_cond,like_st.init_cond,dif_thresh)){
		emsg("like init_cond;"+ref);
	}
	
	if(dif(like.init_cond_prior,like_st.init_cond_prior,dif_thresh)){
		emsg("like init_cond_prior");
	}
	
	if(dif(like.obs,like_st.obs,dif_thresh)){
		cout << like.obs << " " << like_st.obs << " like" << endl;
		emsg("like obs "+ref);
	}
	
	if(dif(like.markov,like_st.markov,dif_thresh)){
		emsg("like markov "+ref);
	}
	
	if(dif(like.nm_trans,like_st.nm_trans,dif_thresh)){ 
		//cout << like.nm_trans << " " << like_st.nm_trans << " " << like.nm_trans-like_st.nm_trans <<  " dif" << endl;
		emsg("like nm_trans"+ref);
	}
	
	if(dif(like.prior,like_st.prior,dif_thresh)){
		emsg("like prior");
	}
		
	if(dif(like.ie,like_st.ie,dif_thresh)){
		emsg("like ie");
	}
	
	if(dif(like.dist,like_st.dist,dif_thresh)){
		emsg("like dist");
	}
	
	if(dif(like.spline_prior,like_st.spline_prior,dif_thresh)){
		emsg("like spline prior "+ref);
	}
	
	const auto &gen_data = model.genetic_data;
	
	if(gen_data.on){
		if(dif(like.genetic_process,like_st.genetic_process,dif_thresh)){
			emsg("like genetic_process: "+ref);
		}
		auto gen_dif = calculate_gen_dif();
		
		if(difi(genetic_value.gen_dif,gen_dif)){
			emsg("gendif prob");
		}
		
		if(dif(like.genetic_obs,like_st.genetic_obs,dif_thresh)){
			emsg("like genetic_obs: "+ref);
		}
	}
	else{
		if(dif(0,like_st.genetic_process,dif_thresh)){
			emsg("genetic_process Should be zero");
		}
		
		if(dif(0,like_st.genetic_obs,dif_thresh)){
			emsg("genetic_obs Should be zero");
		}
	}
	
	if(std::isnan(like.init_cond)) emsg("init_cond nan");
	if(std::isnan(like.init_cond_prior)) emsg("init_cond_prior nan");
	if(std::isnan(like.obs)) emsg("obs nan");
	if(std::isnan(like.prior)) emsg("prior nan");
	if(std::isnan(like.spline_prior)) emsg("spline_prior nan");
	if(std::isnan(like.dist)) emsg("dist nan");
	if(std::isnan(like.markov)) emsg("markov nan");
	if(std::isnan(like.nm_trans)) emsg("nm_trans nan");
	if(std::isnan(like.genetic_process)) emsg("genetic_process nan");
	if(std::isnan(like.genetic_obs)) emsg("genetic_obs nan");
	if(std::isnan(like.ie)) emsg("ie nan");
}


/// Checks the splines are correctly specified
void State::check_spline(string ref)
{	
	auto spline_val_st = spline_val;
	if(spline_val.size() != model.spline.size()) emsg("Wrong size1 "+ref);
			
	spline_init();

	for(auto j = 0u; j < model.spline.size(); j++){
		const auto &val = spline_val[j].val;
		const auto &val_st = spline_val_st[j].val;
		for(auto ti = 0u; ti < T; ti++){
			if(val[ti] != val_st[ti]){
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
		if(dif(prior_prob[th],prior_prob_st[th],dif_thresh)){
			emsg("Prior probability wrong"+ref);
		}
	}
	
	// Checks dist_prob is correct
	auto dist_prob_st = dist_prob;
	dist_prob = model.dist_prob(param_val);
	for(auto th = 0u; th < dist_prob.size(); th++){
		if(dif(dist_prob[th],dist_prob_st[th],dif_thresh)){
			emsg("Dist prob probability wrong");
		}
	}
	
	// Checks spline_prior is correct
	auto spline_prior_st = spline_prior;
	spline_prior = model.spline_prior(param_val);
	for(auto s = 0u; s < spline_prior.size(); s++){
		if(dif(spline_prior[s],spline_prior_st[s],dif_thresh)){
			emsg("Spline Prior probability wrong");
		}
	}
}


/// Checks that popnum_t is correctly specified
void State::check_popnum_t(string ref)
{
	auto popnum_t_store = popnum_t;
	
	popnum_t = model.calculate_popnum_t(species);
	
	if(popnum_t.size() != T){
		emsg("Wrong size2 "+ref);
	}
	
	for(auto ti = 0u; ti < T; ti++){
		if(popnum_t[ti].size() != model.pop.size()) emsg("Wrong size3");
		for(auto k = 0u; k < model.pop.size(); k++){
			if(dif(popnum_t[ti][k],popnum_t_store[ti][k],dif_thresh)){
				//output_dump();
				cout << core() << " " << popnum_t[ti][k] << " " << popnum_t_store[ti][k] << " popn" << endl;
				get_closest_to_boundary();
				emsg("popnum_t problem"+ref);
			}
		}
	}
}


/// Checks that popnum_t is correctly specified
void State::check_popnum_t2(string ref)
{
	auto popnum_t_store = popnum_t;
	
	popnum_t = model.calculate_popnum_t(species);
	
	if(popnum_t.size() != T){
		emsg("Wrong size2 "+ref);
	}
	
	for(auto ti = 0u; ti < T; ti++){
		if(popnum_t[ti].size() != model.pop.size()) emsg("Wrong size3");
		for(auto k = 0u; k < model.pop.size(); k++){
			if(dif(popnum_t[ti][k],popnum_t_store[ti][k],dif_thresh)){
				//output_dump();
				cout << core() << " " << popnum_t[ti][k] << " " << popnum_t_store[ti][k] << " popn" << endl;
				emsg("popnum_t problem2"+ref);
			}
		}
	}
	popnum_t = popnum_t_store; 
}


/// Checks that cpop_st is correctly set	
void State::check_cpop_st(string ref)
{
	for(auto p = 0u; p < model.nspecies; p++){
		auto &sp = model.species[p];
		auto &ssp = species[p];
			
		switch(sp.type){
		case INDIVIDUAL:
			{
				
				auto C = sp.comp_gl.size();
				
				vector <double> cpop(C,0); 
			
				switch(sp.init_cond.type){
				case INIT_POP_FIXED: case INIT_POP_DIST:
					{		
						for(auto c = 0u; c < C; c++){
							cpop[c] = double(ssp.init_cond_val.cnum[c]);
						}
						
						ssp.set_ind_sim_c(0);
						vector <double> cpo(C,0);
						for(auto i = 0u; i < ssp.individual.size(); i++){
							auto c = ssp.ind_sim_c[i];
							if(c != UNSET) cpo[c]++;
						}
								
						for(auto c = 0u; c < C; c++){
							if(cpop[c] != cpo[c]){
								emsg("error init_cond_valh "+ref);
							}
						}			
					}
					break;
				case INIT_POP_NONE: break;
				}
			}
			break;
			
		case POPULATION:
			{
				auto cpop_st_st = ssp.cpop_st;
				ssp.set_cpop_st();
				if(dif(cpop_st_st,ssp.cpop_st,dif_thresh)){
					emsg("error cpop_st"+ref);
				}
			}
			break;
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
		
		if(dif(iesg_store.omega,iesg.omega,dif_thresh)){
			print("om store",iesg_store.omega); print("om2",iesg.omega);
			emsg("Omega problem");
		}
		
		if(dif(iesg_store.omega_Z, iesg.omega_Z,dif_thresh)){
			emsg("Omega_Z problem");
		}
		
		if(model.mode == INF){
			if(dif(iesg_store.omega_inv,iesg.omega_inv,dif_thresh)){
				emsg("Omega_inv problem");
			}
		}	
		
		auto Li_ie_store = ssp.Li_ie[g];
		auto temp = 0.0;
		ssp.likelihood_indeff_group(g,temp);
	
		if(dif(Li_ie_store,ssp.Li_ie[g],dif_thresh)){
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
			if(dif(ssp.tnum_mean_st[tr][ti],tnum_mean_st_st[tr][ti],dif_thresh)){
				emsg("trans num mean prob"+ref);
			}			
			
			if(dif(ssp.trans_num[tr][ti],trans_num_store[tr][ti],dif_thresh)){
				emsg("trans num prob"+ref);
			}
			
			if(dif(ssp.Li_markov_pop[tr][ti],Li_markov_pop_store[tr][ti],dif_thresh)){	
				emsg("Li mark prob");
			}
		}
	}
}


/// Checks all the maps are correcly set
void State::check_maps(string ref)
{
	if(pop_list.size() != 0) emsg("pop list"+ref);
	if(markov_eqn_list.size() != 0) emsg("markov_eqn_list");
	if(trans_list.size() != 0) emsg("trans_list");
		
	for(auto i = 0u; i < model.pop.size(); i++){
		if(pop_map[i] != UNSET) emsg("prob");
	}
		
	for(auto p = 0u; p < species.size(); p++){
		for(auto c = 0u; c < comp_local_map[p].size(); c++){
			if(comp_local_map[p][c] != UNSET) emsg("comp_local_map");
		}
	
		for(auto r = 0u; r <  pop_trans_map[p].size(); r++){ 
			if(pop_trans_map[p][r] != UNSET) emsg("pop_trans_map");
		}
		
		for(auto r = 0u; r <  pop_data_map[p].size(); r++){ 
			if(pop_data_map[p][r] != UNSET) emsg("pop_data_map");
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


/// Checks that the "Observed" property of individual events is correctly set
void State::check_event_observed(unsigned int p, string ref)
{
	auto &ssp = species[p];
	auto individual_store = ssp.individual;
	
	ssp.set_event_observed();
	
	for(auto i = 0u; i < ssp.individual.size(); i++){
		const auto &ind = ssp.individual[i];
		for(auto e = 0u; e < ind.ev.size(); e++){	
			if(ind.ev[e].observed != individual_store[i].ev[e].observed){
				emsg("Problem with observed"+ref);
			}
		}
	}
}


/// Checks that the observation likelihood is correctly specified
void State::check_obs_like(unsigned int p, string ref)
{
	auto &sp = model.species[p];
	auto &ssp = species[p];

	// Checks obs_eqn correctly set
	auto obs_eqn_value_store = ssp.obs_eqn_value;
	ssp.calculate_obs_eqn(seq_vec(sp.obs_eqn.size()));
	
	if(dif(obs_eqn_value_store,ssp.obs_eqn_value,dif_thresh)){
		emsg("obs eqn value difference"+ref);
	}

	// Checks pop_data_cpop and pop_trans_data_tr are correctly set
	auto pop_data_cgl_store = ssp.pop_data_cgl;
	auto pop_trans_data_tgl_store = ssp.pop_trans_data_tgl;
	ssp.calculate_pop_data_cgl_trgl();
	
	if(dif(pop_data_cgl_store,ssp.pop_data_cgl,dif_thresh)){
		emsg("pop_data_cgl error"+ref);
	}
	
	if(dif(pop_trans_data_tgl_store,ssp.pop_trans_data_tgl,dif_thresh)){
		//auto vec = which_dif(pop_trans_data_tgl_store,ssp.pop_trans_data_tgl);
		//auto j = vec[0], i = vec[1];
		//cout << j << " " << i << endl;
		emsg(" pop_trans_data_tgl error"+ref);
	}
	
	auto temp = 0.0;
	
	// Individual observation model
	auto Li_obs_ind_store = ssp.Li_obs_ind;
	ssp.likelihood_obs_ind(seq_vec(sp.nindividual_in),temp);
	
	for(auto i = 0u; i < ssp.Li_obs_ind.size(); i++){
		if(dif(Li_obs_ind_store[i],ssp.Li_obs_ind[i],dif_thresh)){
			emsg("Li_obs_ind error"+ref);
		}
	}
	
	// Population observation model
	auto Li_obs_pop_store = ssp.Li_obs_pop;
	auto pop_data_num_store = ssp. pop_data_num;
	ssp.likelihood_obs_pop(seq_vec(sp.pop_data.size()),temp);
	
	for(auto i = 0u; i < ssp.Li_obs_pop.size(); i++){		
		if(dif(Li_obs_pop_store[i],ssp.Li_obs_pop[i],dif_thresh)){
			emsg("Li_obs_pop error"+ref);
		}
	}
	
	for(auto i = 0u; i < sp.pop_data.size(); i++){
		if(dif(pop_data_num_store[i],ssp.pop_data_num[i],dif_thresh)){
			emsg("pop_data_num error");
		}
	}
	
	// Population transition observation model
	auto Li_obs_pop_trans_store = ssp.Li_obs_pop_trans;
	auto pop_trans_data_num_store = ssp.pop_trans_data_num;

	ssp.likelihood_obs_pop_trans(seq_vec(sp.pop_trans_data.size()),temp);
	
	for(auto i = 0u; i < sp.pop_trans_data.size(); i++){
		if(dif(pop_trans_data_num_store[i],ssp.pop_trans_data_num[i],dif_thresh)){
			emsg("pop_trans_data_num error");
		}
	}
	
	for(auto i = 0u; i < sp.pop_trans_data.size(); i++){
		if(dif(Li_obs_pop_trans_store[i],ssp.Li_obs_pop_trans[i],dif_thresh)){
			emsg("Li_obs_pop_trans error"+ref);
		}
	}
	
	// Parts realted to probability of unobserved transitions
	auto obs_trans_eqn_value_store = ssp.obs_trans_eqn_value;
	auto vec_all = seq_vec(T);
	for(auto e = 0u; e < ssp.obs_trans_eqn_value.size(); e++){
		double temp;
		ssp.likelihood_unobs_trans(e,vec_all,temp);
		
		for(auto ti = 0u; ti < T; ti++){
			if(dif(obs_trans_eqn_value_store[e][ti],ssp.obs_trans_eqn_value[e][ti],dif_thresh)){
				emsg("obs_trans_eqn_value error");		
			}
		}
	}
	
	auto obs_trans_eqn_num_store = ssp.obs_trans_eqn_num;
	ssp.calculate_obs_trans_eqn_num();
	
	if(difi(obs_trans_eqn_num_store,ssp.obs_trans_eqn_num)){
		for(auto m = 0u; m < sp.obs_trans_eqn.size(); m++){
			for(auto ti = 0u; ti < T; ti++){
				if(obs_trans_eqn_num_store[m][ti] != ssp.obs_trans_eqn_num[m][ti]){
						cout << model.eqn[sp.obs_trans_eqn[m]].te << "eqn" << endl;
		
					cout << m << " ti=" << ti << " " << obs_trans_eqn_num_store[m][ti] << "  true="<< ssp.obs_trans_eqn_num[m][ti] << " different" << endl;
				}
			}
		}	
		emsg("obs_trans_eqn_num error "+ref);		
	}
}


/// Checks that the initial condition likelihood is correctly specified
void State::check_init_cond_like(unsigned int p, string ref)
{
	auto &sp = model.species[p];
	auto &ssp = species[p];
	
	const auto &ic = sp.init_cond;
	auto &icv = ssp.init_cond_val;

	if(ic.type != INIT_POP_DIST) return;

	auto foc_cl = ic.focal_cl;

	if(sp.type == INDIVIDUAL){                       // Check cnum agrees with individuals
		if(foc_cl == UNSET){
			auto N_total_unobs_store = icv.N_total_unobs;
			ssp.calculate_N_unobs();
		
			if(dif(N_total_unobs_store,icv.N_total_unobs,dif_thresh)){
				emsg("N_total_unobs problem");
			}
		}
		else{
			auto N_focal_unobs_store = icv.N_focal_unobs;
			ssp.calculate_N_unobs();
		
			if(difi(N_focal_unobs_store,icv.N_focal_unobs)){
				emsg("N_focal_unobs problem");
			}
			
			auto cnum_reduce_store = icv.cnum_reduce;
			
			const auto &claa = sp.cla[foc_cl];
			for(auto c = 0u; c < claa.ncomp; c++){
				for(auto cred = 0u; cred < ic.N_reduce; cred++){
					icv.cnum_reduce[c][cred] = 0;
				}
			}
			
			for(const auto &ind : ssp.individual){
				if(ind.ev.size() > 0){
					if(ind.ev[0].type == ENTER_EV && ind.ev[0].t == model.details.t_start){
						auto c = ind.ev[0].c_after;
						const auto &crr = ic.comp_reduce_ref[c];
						icv.cnum_reduce[crr.c][crr.cred]++; 
					}
				}					
			}
			
			if(difi(cnum_reduce_store,icv.cnum_reduce)){
				for(auto c = 0u; c < cnum_reduce_store.size(); c++){
					cout << c << "c" << endl;
					for(auto val : cnum_reduce_store[c]) cout << val <<","; 
					cout << "store" << endl;
					for(auto val : icv.cnum_reduce[c]) cout << val <<","; 
					cout << "true" << endl;
				}
				
				emsg("cnum_reduce error");
			}
		}
		
		auto cnum_store = icv.cnum;
	
		for(auto cgl = 0u; cgl < sp.comp_gl.size(); cgl++){
			icv.cnum[cgl] = 0;
		}	
			
		for(const auto &ind : ssp.individual){
			if(ind.ev.size() > 0){
				if(ind.ev[0].type == ENTER_EV && ind.ev[0].t == model.details.t_start){
					auto c = ind.ev[0].c_after;
					icv.cnum[c]++; 
				}
			}					
		}
			
		if(difi(cnum_store,icv.cnum)){
			emsg("cnum error");
		}
	}		
	
	if(foc_cl != UNSET){
		auto frac_comb_store = icv.frac_comb;
		icv.frac_comb = model.frac_combine(p,icv.frac_focal);
			
		if(dif(frac_comb_store,icv.frac_comb,dif_thresh)){
			emsg("frac_comp problem");
		}
		
		auto cnum_store = icv.cnum;
		model.combine_cnum_reduce(p,icv);
		if(difi(cnum_store,icv.cnum)){
			emsg("combine_cnum_reduce problem");
		}
	}
	
	auto Li_init_cond_store = ssp.Li_init_cond;
	double temp;
	ssp.likelihood_init_cond(temp);
	
	if(dif(Li_init_cond_store,ssp.Li_init_cond,dif_thresh)){
		emsg("Li_init_cond error"+ref);
	}
}

/// Checks that the initial condition likelihood is correctly specified
void State::check_init_cond_prior(string ref)
{
	auto init_cond_prior_store = like.init_cond_prior;
	prior_init_cond(like.init_cond_prior);
	
	if(dif(init_cond_prior_store,like.init_cond_prior,dif_thresh)){
		emsg("init_cond_prior prob"+ref);	
	}
}


/// Checks that the linearised version of the equations give the same as non-linearised
void State::check_linearise() const 
{
	auto ti = (unsigned int)(ran()*T);

	const auto &popnum = popnum_t[ti];
	
	for(auto e = 0u; e < model.eqn.size(); e++){
		const auto &eq = model.eqn[e];
		if(eq.linearise.on){
			auto val = eq.calculate(ti,popnum,param_val,spline_val);
			auto val2 = eq.calculate_linearise_check(ti,popnum,param_val,spline_val);
			if(dif(val,val2,dif_thresh)){
				cout << ti << "ti" << endl;
				cout << eq.te << endl;
				cout << val << " " << val2 << "val" << endl;
				emsg("Equation different");
			}
		}
	}
}


/// Checks the spline store value are all UNSET
void State::check_spline_store(string ref)
{
	for(auto &spl : spline_val){
		for(auto &val : spl.store){
			if(val != UNSET){
				emsg("spline store error"+ref);
			}
		}
	}
}



/// Checks that quantities in genetic_value are correctly specified
void State::check_genetic_value(string ref)
{
	const auto &gen_data = model.genetic_data;

	if(model.trans_tree == false) return;

	if(gen_data.on){                                 // Checks mutation rate 
		auto mut_rate = model.eqn[gen_data.mut_rate.eq_ref].calculate_param_only(param_val);
		auto seq_var = model.eqn[gen_data.seq_var.eq_ref].calculate_param_only(param_val);
		
		if(dif(mut_rate,genetic_value.mut_rate,dif_thresh)){
			emsg("mut rate is wrong");
		}
		
		if(dif(seq_var,genetic_value.seq_var,dif_thresh)){
			emsg("seq var is wrong");
		}
	}
	
	/// Checks that ev.ind_inf_from is correctly specified
	for(auto p = 0u; p < model.species.size(); p++){
		const auto &sp = model.species[p];
		if(sp.trans_tree){
			const auto &ssp = species[p];

			for(const auto &ind : ssp.individual){
				for(auto ev : ind.ev){
					const auto &iif = ev.ind_inf_from;
				
					switch(ev.type){
					case ENTER_EV:
						if(sp.comp_gl[ev.c_after].infected == true){
							if(iif.p != ENTER_INF) emsg("SHould be ENTER_INF1a");
							if(ev.inf_node_ref > CODE) emsg("Should not be code1");
						}
						else{
							if(iif.p != UNSET){
								emsg("SHould not be ENTER_INF1b");
							}
							if(ev.inf_node_ref != UNSET) emsg("Should not be code2");
						}
						break;
						
					case M_TRANS_EV:
						{
							const auto &tra = sp.tra_gl[ev.tr_gl];
							if(tra.infection.type == TRANS_INFECTION){
								if(iif.p != OUTSIDE_INF && iif.p > CODE){
									species[p].print_event(ind.ev);
									emsg("SHould be ENTER_INF2");
								}

								if(ev.inf_node_ref > CODE){
									emsg("Should not be code3");
								}
								
								if(iif.p == OUTSIDE_INF){
									if(iif.po != UNSET) emsg("SHould be UNSET1");
									if(iif.pref != UNSET) emsg("SHould be UNSET2");
									if(iif.w != 1){
										emsg("SHould be UNSET3");
									}
								}
								else{                              // Checks pref and po consistent
									auto m = tra.markov_eqn_ref;
									const auto &me = sp.markov_eqn[m];
									const auto &eq = model.eqn[me.eqn_ref];
									const auto &lin = eq.linearise;
									if(!lin.on) emsg("Linearisation should be on");
			
									if(iif.po != eq.pop_ref[iif.pref]){
										emsg("pref and po not consistent");
									}
									
									if(iif.p != model.pop[iif.po].sp_p){
										emsg("iif species is not right");
									}
									
									// Checks that the weight is correct
									auto w = get_w_from_indinffrom(iif);
									if(dif(w,iif.w,dif_thresh)){ 
										emsg("iif.w not right");
									}
								}
							}
							else{
								if(iif.p != UNSET) emsg("SHould not be ENTER_INF2");
								if(ev.inf_node_ref != UNSET){
									emsg("Should not be code4");
								}
							}
						}
						break;
					
					default:
						if(iif.p != UNSET) emsg("SHould be ENTER_INF 3");
						if(ev.inf_node_ref != UNSET) emsg("Should not be code5");
						break;
					}
				}
			}
		}
		
		// Check iif from the point of view of nodes
		for(auto n = 0u; n < genetic_value.inf_node.size(); n++){
			const auto &in = genetic_value.inf_node[n];
			
			const auto &event = species[in.p].individual[in.i].ev;
			
			const auto &ev = event[in.e];
			
			if(in.t_start != ev.t) emsg("t_start not agree");
		
			auto e = 0u;
			auto c = UNSET;
			for(auto k = 0u; k < in.inf_ev.size(); k++){
				const auto &ie = in.inf_ev[k];
				
				if(ie.type == INFECT_OTHER){
					auto t = ie.t;
					
					const auto &in_inf = genetic_value.inf_node[ie.index];
					
					auto ev = species[in_inf.p].individual[in_inf.i].ev[in_inf.e];
					if(ev.type != M_TRANS_EV) emsg("SHould be markovian");
					if(ev.t != t) emsg("Porblem with time");
					
					const auto &iif = ev.ind_inf_from;
					
					const auto &tra = model.species[in_inf.p].tra_gl[ev.tr_gl]; 
					if(tra.infection.type != TRANS_INFECTION) emsg("SHould be and infection");
					
					while(e < event.size() && event[e].t < t){
						c = event[e].c_after;
						e++;
					}
				
					const auto &inf_c = model.inf_cause[in_inf.p][ev.tr_gl][in.p][c];
					if(inf_c.eq_ref == UNSET) emsg("Cannot fint equation for iif");
						
					if(iif.p != in.p) emsg("check p not agree");
					if(iif.i != in.i) emsg("check i not agree");
					if(iif.pref != inf_c.pref){
						const auto &eq = model.eqn[inf_c.eq_ref];
						for(auto j = 0u; j <eq.pop_ref.size(); j++){
							cout << j << " " <<model.pop[eq.pop_ref[j]].name << "pop ref" << endl;
						}
						cout << iif.pref << " " << inf_c.pref <<" jj" << endl;
						emsg("check pref not agree");
					}
					if(iif.po != inf_c.po) emsg("check po not agree");
					
					auto w = get_w_from_indinffrom(iif);
					if(dif(w,iif.w,dif_thresh)){
						emsg("check w not agree");				
					}
				}
			}
		}
	}
	
	// Checks obs_node_ref
	if(genetic_value.obs_node_ref.size() != gen_data.obs.size()){
		emsg("obs_node_ref wrong size4");
	}
	
	auto nobs_not_infected = 0u;
	for(auto m = 0u; m < genetic_value.obs_node_ref.size(); m++){
		const auto &onr = genetic_value.obs_node_ref[m];
		if(onr.node == UNSET) nobs_not_infected++;
		else{		
			const auto &node_ev = genetic_value.inf_node[onr.node].inf_ev[onr.index];
			if(node_ev.type != GENETIC_OBS) emsg("Should be genetic obs");
			if(node_ev.index != m){
				emsg("node ev index agree");
			}
			if(node_ev.t != gen_data.obs[m].t){
				emsg("gen data obs t not agree");
			}
		}
	}
	
	// Checks nobs_not_infected 
	if(genetic_value.nobs_not_infected != nobs_not_infected){
		emsg("nobs_not_infected not agree"+ref);
	}
	
	// Checks iif is specified properly for for ENTER_INF and OUTSIDE_INF
	for(auto p = 0u; p < model.species.size(); p++){
		const auto &sp = model.species[p];
		const auto &ssp = species[p];
	
		if(sp.trans_tree){
			for(const auto &ind : ssp.individual){ 	
				for(const auto &ev : ind.ev){
					const auto &iif = ev.ind_inf_from;
					if(iif.p == UNSET){
						if(iif.i != UNSET) emsg("i unset1");
						if(iif.pref != UNSET) emsg("pref unset");
						if(iif.po != UNSET) emsg("po unset");
						if(iif.w != UNSET) emsg("w unset1");
					}
					else{
						if(iif.p == ENTER_INF){
							if(iif.i != UNSET) emsg("i unset2");
							if(iif.pref != UNSET) emsg("pref unset");
							if(iif.po != UNSET) emsg("po unset");
							if(iif.w != UNSET) emsg("w unset1a");
						}
						else{
							if(iif.p == OUTSIDE_INF){
								if(iif.i != UNSET) emsg("i unset3");
								if(iif.pref != UNSET) emsg("pref unset");
								if(iif.po != UNSET) emsg("po unset");
								if(iif.w != 1) emsg("w should be 1a");
							}
							else{
								if(iif.i == UNSET) emsg("i unset4");
								if(iif.pref == UNSET) emsg("pref unset");
								if(iif.po == UNSET) emsg("po unset");
								if(iif.w == UNSET) emsg("w unset2");
							}
						}
					}
				}
			}
		}
		else{
			for(const auto &ind : ssp.individual){ 	
				for(const auto &ev : ind.ev){
					const auto &iif = ev.ind_inf_from;
	
					if(iif.p != UNSET) emsg("p unset");
					if(iif.i != UNSET) emsg("i unset5");
					if(iif.pref != UNSET) emsg("pref unset");
					if(iif.po != UNSET) emsg("po unset");
					if(iif.w != UNSET) emsg("w unset3");
				}
			}
		}
	}
						
	// Checks inf_node
	const auto &inf_node = genetic_value.inf_node;
	
	if(!gen_data.on){ // Checks that mut_num is set to UNSET if no genetic data
		for(auto n = 0u; n < inf_node.size(); n++){
			const auto &in = inf_node[n];
			for(const auto &iev : in.inf_ev){
				if(iev.mut_num != UNSET) emsg("Mut_num should be UNSET");
			}
		}
		
		for(const auto &o : genetic_value.inf_origin){
			if(o.mut_num != UNSET) emsg("Mut_num should be UNSET");
		}
	}
			
	for(auto n = 0u; n < inf_node.size(); n++){
		const auto &in = inf_node[n];
		const auto &eve = species[in.p].individual[in.i].ev[in.e];
		const auto &sp = model.species[in.p];
		
		const auto &from = in.from;
		if(from.node != ENTER_INF && from.node != OUTSIDE_INF){
			if(inf_node[from.node].inf_ev[from.index].index != n) emsg("from problem");
		}
		
		if(in.t_start != eve.t){
			cout << in.p << " " << in.i << " "<< in.e  << "  inpie" << endl;
			cout << n << " for node n" << endl;
			cout << in.t_start << " " << eve.t << "  time" << endl;
			cout << inf_node.size() <<" nn" << endl;
			emsg("inf_node Times do not agree");
		}
		
		if(eve.inf_node_ref != n) emsg("inf_node_ref wrong");
		
		if(eve.ind_inf_from.p == UNSET) emsg("should not be unset");
		
		if(eve.ind_inf_from.p == ENTER_INF){
			if(eve.type != ENTER_EV) emsg("Should be infected entry");	
			if(sp.comp_gl[eve.c_after].infected != true) emsg("Shoudl be infection");
			if(in.from.node != ENTER_INF) emsg("SHould not have a from node1");
			
			if(in.from.index == ENTER_INF || in.from.index == OUTSIDE_INF) emsg("index is wrong1");
			const auto &io = genetic_value.inf_origin[in.from.index];
			if(io.node != n) emsg("inforigin ref wrong1");
		}		
		else{
			if(eve.ind_inf_from.p == OUTSIDE_INF){
				if(eve.type != M_TRANS_EV) emsg("Should be outside transition");	
				if(sp.comp_gl[eve.c_after].infected != true) emsg("Shoudl be infection");
				if(in.from.node != OUTSIDE_INF) emsg("SHould not have a from node2");
				
				if(in.from.index == ENTER_INF || in.from.index == OUTSIDE_INF) emsg("index is wrong2");
			
				const auto &io = genetic_value.inf_origin[in.from.index];
				if(io.node != n){
					emsg("inforigin ref wrong2");
				}
			}
			else{
				const auto &in_from = genetic_value.inf_node[in.from.node];
			
				if(eve.ind_inf_from.p == UNSET) emsg("Species should be set");
				
				if(in_from.p != eve.ind_inf_from.p) emsg("in from p not agree");
				if(in_from.i != eve.ind_inf_from.i) emsg("in from i not agree");
			}
		}
		
		for(auto e = 1u; e < in.inf_ev.size(); e++){
			if(in.inf_ev[e-1].t > in.inf_ev[e].t){
				emsg("inf_ev time order problem");
			}
		}
		
		for(auto e = 0u; e < in.inf_ev.size(); e++){
			const auto &iev = in.inf_ev[e];
			switch(iev.type){
			case INFECT_OTHER:
				{
					const auto &in_to = genetic_value.inf_node[iev.index];
					
					if(in_to.from.node != n) emsg("node not match");
					if(in_to.from.index != e) emsg("node not match index");
					
					const auto &eve_to = species[in_to.p].individual[in_to.i].ev[in_to.e];
					if(eve_to.ind_inf_from.p != in.p) emsg("ind_inf_from p prob");
					if(eve_to.ind_inf_from.i != in.i){
						emsg("ind_inf_from i prob");
					}
				}
				break;
				
			case GENETIC_OBS:
				{
					const auto &onr = genetic_value.obs_node_ref[iev.index];
					if(onr.node != n) emsg("obs_node_ref node is wrong");
					if(onr.index != e) emsg("obs_node_ref event is wrong");
					if(iev.t != gen_data.obs[iev.index].t) emsg("obs_node_ref t wrong");
				}
				break;
			}
		}
 	}
	
	// Checks inf_origin
	for(auto &io : genetic_value.inf_origin){
		const auto &in = genetic_value.inf_node[io.node];
		const auto &eve = species[in.p].individual[in.i].ev[in.e];
		const auto &sp = model.species[in.p];
		
		if(sp.comp_gl[eve.c_after].infected != true) emsg("Shoudl be infection");
	}

	// Checks individual time lines
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
				
		if(sp.trans_tree){		
			for(auto i = 0u; i < species[p].individual.size(); i++){
				const auto &ind = species[p].individual[i];
				for(auto e = 0u; e < ind.ev.size(); e++){
					const auto &eve = ind.ev[e];
					const auto &fr = eve.ind_inf_from;
					
					switch(eve.type){
					case ENTER_EV:
						if(sp.comp_gl[eve.c_after].infected == true){
							if(eve.inf_node_ref == UNSET) emsg("inf_node_ref should be set2");	
							const auto &in = genetic_value.inf_node[eve.inf_node_ref];
							if(in.p != p) emsg("inf_node_ref p wrong");
							if(in.i != i) emsg("inf_node_ref i wrong");
							if(in.e != e) emsg("inf_node_ref e wrong");
							if(fr.p != ENTER_INF) emsg("ind_inf_from.p should not be set1");
						}
						else{
							if(fr.p != UNSET) emsg("ind_inf_from.p should not be set2");
						}
						
						if(fr.i != UNSET) emsg("ind_inf_from.i should not be set");
						if(fr.pref != UNSET) emsg("ind_inf_from.i should not be set");
						break;
					
					case M_TRANS_EV:
						{
							const auto &tra = sp.tra_gl[eve.tr_gl];
							if(tra.infection.type == TRANS_INFECTION){	
								auto m = tra.markov_eqn_ref;
								const auto &me = sp.markov_eqn[m];
								const auto &eq = model.eqn[me.eqn_ref];
								const auto &lin = eq.linearise;
								if(!lin.on) emsg("Linearisation should be on");
						
								if(fr.p == UNSET) emsg("fr.p should not be unset");
								
								if(fr.p == OUTSIDE_INF){
									if(eve.inf_node_ref == UNSET) emsg("inf_node_ref should be set1");	
									if(fr.pref != UNSET){
										emsg("pref prob");
									}
								}
								else{
									if(eve.inf_node_ref == UNSET) emsg("inf_node_ref should be set1");	
									const auto &in = genetic_value.inf_node[eve.inf_node_ref];
									if(in.p != p) emsg("inf_node_ref p wrong");
									if(in.i != i) emsg("inf_node_ref i wrong");
									if(in.e != e) emsg("inf_node_ref e wrong");
									
									if(fr.pref >= eq.pop_ref.size()) emsg("pref prob2");		
								}
								break;
							}
						}
						break;
						
					default:
						if(eve.inf_node_ref != UNSET) emsg("inf_node_ref should be unset");
						if(fr.p != UNSET) emsg("fr.p should be unset");
						if(fr.i != UNSET) emsg("fr.i should be unset");
						if(fr.pref != UNSET) emsg("fr.pref should be unset");
						break;
					}
				}
			}
		}
	}
	
	// Checks gen_dif
	auto gen_dif = calculate_gen_dif();
	auto nob = gen_dif.size();
	auto fl = false;
	for(auto j = 0u; j < nob; j++){
		for(auto i = 0u; i < nob; i++){
			if(dif(genetic_value.gen_dif[j][i],gen_dif[j][i],dif_thresh)){
				fl = true;
				emsg("gen_dif problem1a");
			}
		}
	}
	if(fl == true) emsg("gen_dif problem1b");
	
	
	// Checks that inf_node are time ordered
	for(auto n = 0u; n < genetic_value.inf_node.size(); n++){
		const auto &in = genetic_value.inf_node[n];
		auto t = in.t_start;
		for(auto e = 0u; e < in.inf_ev.size(); e++){
			auto tt = in.inf_ev[e].t;
			if(tt < t) emsg("inf_node time oprder prob");
			tt = t;
		}
	}
	
	// Checks t_rec correctly specified
	for(auto n = 0u; n < genetic_value.inf_node.size(); n++){
		const auto &in = genetic_value.inf_node[n];
		if(in.t_rec != find_t_rec(in)) emsg("t_rec is wrong");
	}
	
	// Checks that infections are within the infected periods of individuals
	if(model.genetic_data.on){   
		for(auto n = 0u; n < genetic_value.inf_node.size(); n++){
			const auto &in = genetic_value.inf_node[n];
			if(in.inf_ev.size() > 0){
				auto t_rec = round_up(find_t_rec(in));

				for(auto k = 0u; k < in.inf_ev.size(); k++){
					const auto &iev = in.inf_ev[k];
					if(iev.type == INFECT_OTHER && iev.t > t_rec){
						const auto &ssp = species[in.p];
						const auto &ev = ssp.individual[in.i].ev;
						ssp.print_event(ev);
						emsg("Infection not within infectious period "+ref);
					}
				}
			}
		}
		
		// Checks that all observations are correctly within infection period
		for(auto n = 0u; n < genetic_value.inf_node.size(); n++){
			const auto &in = genetic_value.inf_node[n];
			if(in.i >= model.species[in.p].nindividual_obs){ // Unobserved individual
				for(const auto &iev : in.inf_ev){
					if(iev.type == GENETIC_OBS) emsg("Should not be genetic observation");
				}		
			}				
			else{	
				auto t_start = in.t_start;
				auto t_rec = in.t_rec;
				
				vector <unsigned int> list_in;
				for(const auto &iev : in.inf_ev){
					if(iev.type == GENETIC_OBS) list_in.push_back(iev.index);
				}			
				
				vector <unsigned int> list_obs;
				for(const auto &igo : gen_data.ind_gen_obs[in.p][in.i]){
					if(igo.t >= t_start && igo.t <= t_rec) list_obs.push_back(igo.m);
				}
				
				if(!equal_vec(list_in,list_obs)){
					cout << t_start << " " << t_rec << " time range" << endl;
					cout << "Individual " << in.i << ":" << endl;
					
					cout << "List from data: ";
					for(auto m : list_obs) cout <<  gen_data.obs[m].t << ", ";
					cout << endl;
					
					cout << "List from inf_node: ";
					for(auto m : list_in) cout <<  gen_data.obs[m].t << ", ";
					cout << endl;
					
					emsg("Observations not correctly implemented");
				}			
			}
		}
	}
}


/// CHecks that popnum_ind is correctly specified
void State::check_popnum_ind(string ref)
{
	auto popnum_ind_store = popnum_ind;
	auto species_store = species;
	
	calculate_popnum_ind();
	
	// Checks popnum_ind is the same
	for(auto ti = 0u; ti < T; ti++){
		for(auto po = 0u; po < model.pop.size(); po++){
			if(popnum_ind[ti][po].size() != popnum_ind_store[ti][po].size()){
				emsg("popnum_ind size wrong:"+ref);
			}
			
			for(auto j = 0u; j < popnum_ind[ti][po].size(); j++){
				const auto &pir = popnum_ind[ti][po][j];
				auto i = pir.i;
				
				auto k = 0u; while(k < popnum_ind_store[ti][po].size() && popnum_ind_store[ti][po][k].i != i) k++;
				if(k == popnum_ind_store[ti][po].size()){
					emsg("popnum_ind_store could not find");
				}
				
				const auto &pir2 = popnum_ind_store[ti][po][k];
				if(pir.i != pir2.i) emsg("popnum_ind i wrong");
				if(dif(pir.w,pir2.w,dif_thresh)){
					emsg("popnum_ind w wrong "+ref);
				}
			}
		}
	}
	
	// Checks popnum_ind reference is correct
	for(auto ti = 0u; ti < T; ti++){
		for(auto po = 0u; po < model.pop.size(); po++){
			auto p = model.pop[po].sp_p;
			for(auto j = 0u; j < popnum_ind[ti][po].size(); j++){
				const auto &pir = popnum_ind[ti][po][j];
				const auto &pnir = species[p].individual[pir.i].popnum_ind_ref[pir.index];
				if(pnir.po != po) emsg("popnum_ind_ref po wrong"+ref);
				if(pnir.ti != ti) emsg("popnum_ind_ref ti wrong");
				if(pnir.index != j) emsg("popnum_ind_ref j wrong");
			}
		}
	}
}


/// Checks that add rem and move events are strictly followed
void State::check_add_move_rem(string ref)
{
	for(auto p = 0u; p < species.size(); p++){
		const auto &sp = model.species[p];
		const auto &ssp = species[p];
		
		if(sp.type == INDIVIDUAL){
			for(auto i = 0u; i < sp.nindividual_in; i++){
				const auto &indd = sp.individual[i];
				const auto &ind = ssp.individual[i];
				vector <EventData> ed;
				for(auto &ev : ind.ev){
					switch(ev.type){
					case ENTER_EV: case LEAVE_EV: case MOVE_EV:
						{
							auto fl_add = true;
							if(ev.type == ENTER_EV){
								fl_add = false;
								for(const auto &evd : indd.ev){
									if(evd.type == ENTER_EV) fl_add = true;
								}
							}
							
							if(fl_add == true){
								EventData evd;
								evd.type = ev.type;
								evd.move_c = ev.move_c;
								evd.cl = ev.cl;
								evd.t = ev.t;
								evd.tr = ev.tr_gl;  //???Not sure about this
								ed.push_back(evd);
							}
						}
						break;
						
					default: break;
					}
				}
				
				auto fl = false;
				if(indd.ev.size() != ed.size()){ 
					// TO DO SORT OUT ERROR
					
					//cout << indd.ev.size() << " "<<ed.size() << "indd ev size wrong" << endl;
					//fl = true;
				}
				else{
					for(auto e = 0u; e < indd.ev.size(); e++){
						const auto &ev1 = indd.ev[e];
						const auto &ev2 = ed[e];
						
						if(ev1.type != ev2.type){ cout << e << ": yype wrong" << endl; fl = true;}
						if(ev1.move_c != ev2.move_c){
							cout << e << ": Mmve_c wrong" << endl; fl = true;
						}
						
						if(ev1.cl != ev2.cl){
							cout << e << ": cl wrong" << endl; fl = true;
						}
						
						if(dif(ev1.t,ev2.t,dif_thresh)){ 
							cout << e << ": t wrong" << " " << ev1.t << " " << ev2.t << " " << ev1.t - ev2.t  <<endl; 
							fl = true;
						}
					}
				}
				
				if(fl == true){
					cout << "Individual " << ind.name << ":" << endl;
					print_ev_data("Data",indd.ev,p);
					print_ev_data("Ind",ed,p);
					
					ssp.print_event("Event",ind);		
					emsg("PROBLEM WITH DATA"+ref);
				}
			}
		}
	}
}


/// Prints a set of event data
void State::print_ev_data(string te, const vector <EventData> &event_data, unsigned int p) const 
{
	const auto &sp = model.species[p];
	
	cout << te << ": ";
	for(const auto &ev : event_data){
		switch(ev.type){
		case ENTER_EV: cout << "enter"; break;
		case LEAVE_EV: cout << "leave"; break;
		case MOVE_EV: cout << "move to " << sp.cla[ev.cl].comp[ev.move_c].name; break;
		case M_TRANS_EV: cout << "markov trans " << sp.cla[ev.cl].tra[ev.tr].name; break;	
		case NM_TRANS_EV: cout << "non-markov trans " << sp.cla[ev.cl].tra[ev.tr].name; break;
		}
		cout << "," << ev.t << "  ";
	}
	cout << endl;
}


/// Scans a variable and plots the variation in likelihood
void State::scan_variable(string name, double min, double max) 
{
	const auto imax = 100u;
	const auto &vec = model.param_vec;
	auto th = 0u; while(th < vec.size() && vec[th].name != name) th++;
	if(th == vec.size()) emsg("Could not find");

	ofstream scan("Output/scan.txt");
	for(auto i = 0u; i < imax; i++){
		cout << i << " scan" << endl;
		auto val = min+(max-min)*double(i)/imax;
	
		auto param_store = param_val;
		
		param_val[th] = val;
		
		likelihood_from_scratch();
		
		auto Li = like.init_cond+
		like.init_cond_prior+
		like.obs+
		like.prior+ 
		like.spline_prior+ 
		like.dist+
		like.markov+ 
		like.nm_trans+
		like.genetic_process+
		like.genetic_obs+
		like.ie;
		
		scan << val << " " << Li << endl;
	}
}


/// Checks if any events are near to a division boundary
void State::check_event_boundary(string ref)
{
	for(auto p = 0u; p < model.species.size(); p++){	
		const auto &sp = model.species[p];
		
		if(sp.type == INDIVIDUAL){	
			const auto &ssp = species[p];
			for(const auto &ind : ssp.individual){
				if(events_near_div(ind.ev,model.details)){
					cout << ind.name <<":" << endl;
					for(auto ev  : ind.ev) cout << ev.t << " t" << endl; 
					emsg("Event near a boundary"+ref);
				}
			}
		}
	}
}



/// Checks if individual or fixed effect are too large or small
void State::check_effect_out_of_range()
{
	string ty = "value"; if(model.mode == INF) ty = "prior";
		
	for(auto p = 0u; p < model.species.size(); p++){	
		const auto &sp = model.species[p];
		
		if(sp.type == INDIVIDUAL){	
			const auto &ssp = species[p];
			for(const auto &ind : ssp.individual){
				for(auto e = 0u; e < sp.ind_effect.size(); e++){
					auto val = ind.exp_ie[e];
					if(val > EFFECT_MAX || val < EFFECT_MIN){
						auto ie_name = sp.ind_effect[e].name;
				
						stringstream ss;
						ss << "Individual effect '[" << ie_name << "]' has a value '" << val << "' which is out of the maximum allowable range " << EFFECT_MIN << "-" << EFFECT_MAX << ". ";
						ss << "Please restrict the " << ty << " for the variance parameter Ω^" << ie_name << "," << ie_name;
								
						emsg(ss.str()); 
					}
				}
				
				for(auto e = 0u; e < sp.fix_effect.size(); e++){
					auto val = ind.exp_fe[e];
					if(val > EFFECT_MAX || val < EFFECT_MIN){
						auto fe_name = sp.fix_effect[e].name;
				
						stringstream ss;
						ss << "Fixed effect '〈" << fe_name << "〉' has a value '" << val << "' which is out of the maximum allowable range " << EFFECT_MIN << "-" << EFFECT_MAX << ". ";
						ss << "Please restrict the " << ty << " for the fixed effect parameter " << fe_char << "^" << fe_name;
								
						emsg(ss.str()); 
					}
				}
			}
		}
	}
}


/// Dumps state information (for diagnostic purposes)
void State::output_dump() const
{
	model.print_param(param_val);
	for(const auto &ssp : species){
		for(const auto &ind : ssp.individual){
			cout << ind.name << ":  exp_ie:";
			for(auto val : ind.exp_ie) cout << val << ",";
			cout << "  exp_fe:";
			for(auto val : ind.exp_fe) cout << val << ",";
			cout << endl;
		}
	}	
}


/// Gets the closest distance from a transition event and a division boundary
void State::get_closest_to_boundary()
{
	auto t_start = model.details.t_start, dt = model.details.dt;
		
	for(const auto &ssp : species){
		if(ssp.type == INDIVIDUAL){
			auto dmin = LARGE;
			for(const auto &ind : ssp.individual){
				for(const auto &ev : ind.ev){
					switch(ev.type){
					case NM_TRANS_EV: case M_TRANS_EV:
						{
							auto f = (ev.t-t_start)/dt;
							auto d = f-(unsigned int)(f+0.5);
							if(d < 0) d = -d;
							if(d < dmin) dmin = d;
						}
						break;
						
					default: break;
					}
				}
			}
			
			cout << dmin << " smallest distance" << endl;
		}
	}
}



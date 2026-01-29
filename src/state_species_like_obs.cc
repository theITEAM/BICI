// Functions associated with observation likelihoods

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
 
using namespace std;

#include "state_species.hh"
#include "utils.hh"
#include "matrix.hh"

/// Calculates values for the observation equations
vector <double> StateSpecies::calculate_obs_eqn(const vector <unsigned int> &list)
{
	const auto &precalc = param_val.precalc;
	
	vector <double> store;
	for(auto i : list){
		store.push_back(obs_eqn_value[i]);
		obs_eqn_value[i] = eqn[sp.obs_eqn[i]].calculate_param_ti_fix(precalc);
	}

	return store;
}


/// Calculates the observation likelihood for a list of individuals
vector <double> StateSpecies::likelihood_obs_ind(const vector <unsigned int> &list, double &like_ch)
{	
	vector <double> store;
	
	if(mode == SIM) return store;
		
	double t_end = T;
	for(auto i : list){
		auto Li = 0.0;
		
		const auto &indi = sp.individual[i];
		const auto &obs = indi.obs;
		const auto &ind = individual[i];
				
		auto c = UNSET;
		auto m = 0u;
		auto ev = ind.ev;
		
		for(auto k = 0u; k <= ev.size(); k++){
			double t;	
			if(k < ev.size()) t = ev[k].tdiv; else t = t_end;
			
			while(m < obs.size() && obs[m].tdiv <= t){
				const auto &ob = obs[m];
				
				switch(ob.type){
				case OBS_TRANS_EV: case OBS_SOURCE_EV: case OBS_SINK_EV: 
					if(obs[m].tdiv != t) Li += LI_WRONG; 
					else{
						if(k < ev.size()){
							const auto &e = ev[k];
							if(e.observed == true){
								Li += like_trans_prob_obs(e.tr_gl,obs[m]);
							}
						}
					}
					break;	
				case OBS_COMP_EV: Li += like_comp_obs(c,ob); break;
				case OBS_TEST_EV: Li += like_diag_test(c,ob); break;
				}
				
				m++;
			}
		
			if(k < ev.size()){
				const auto &e = ev[k];
				c = e.c_after;
				
				if(e.type == ENTER_EV) Li += like_ind_init_obs(c,indi.enter_ref);
			}	
		}
		
		like_ch += Li-Li_obs_ind[i];
		store.push_back(Li_obs_ind[i]);
		Li_obs_ind[i] = Li;
	}
	
	return store;
}


/// Returns the observation model probability
double StateSpecies::obs_mod_probability(double num, ObsModelVariety type, double value, double obs_mod_val) const
{
	switch(type){
	case NORMAL_OBS:
		return normal_probability(value,num,obs_mod_val);
		
	case POISSON_OBS:
		{
			auto lam = num; if(lam < 0.5) lam = 0.5;
			auto i = round_int(value); if(i < 0) emsg("Data cannot be negative");
			return poisson_probability(i,lam);
		}
		 
	case NEGBIN_OBS:
		{
			auto i = round_int(value); if(i < 0) emsg("Data cannot be negative");
			return neg_binomial_probability(i,num,obs_mod_val);		
		}
	}
	
	return UNSET;
}


/// The observation probability for population data
vector <double> StateSpecies::likelihood_obs_pop(const vector <unsigned int> &list, double &like_ch)
{	
	vector <double> store;
	
	if(mode == SIM) return store;
	
	const auto &popd = sp.pop_data;
	const auto &pop_filt = sp.pop_filter;

	for(auto i : list){
		const auto &pd = popd[i];
		
		const auto &pf = pop_filt[pd.ref];
		
		double num;
		if(pd.time_vari) num = calculate_pop_num(i,pd.tdiv,pf.c_nonzero,pd.comp_obs_mod_ref);
		else num = calculate_pop_num(i,pd.tdiv,pf.c_nonzero,pf.comp_obs_mod_ref);
		
		auto val = obs_mod_probability(num,pd.type,pd.value,pd.obs_mod_val);
		
		like_ch += val - Li_obs_pop[i];
		
		store.push_back(Li_obs_pop[i]);
		store.push_back(pop_data_num[i]);
		Li_obs_pop[i] = val;
		pop_data_num[i] = num;
	}
	
	return store;
}


/// Calculates the population number based on  
double StateSpecies::calculate_pop_num(unsigned int i, double t, const vector <unsigned int> &c_nonzero, const vector <unsigned int> &comp_obs_mod_ref) const
{
	auto num = 0.0;	
	switch(type){
	case POPULATION: case DETERMINISTIC:
		{
			auto ti = get_ti(t);
			if(ti == T) ti--;

			for(auto c : c_nonzero){
				auto val = obs_eqn_value[comp_obs_mod_ref[c]];
				if(val != 0) num += val*cpop_st[ti][c];
			}
		}
		break;
		
	case INDIVIDUAL:
		{
			for(auto c : c_nonzero){
				auto val = obs_eqn_value[comp_obs_mod_ref[c]];
				if(val != 0) num += val*pop_data_cgl[i][c];
			}
		}
		break;
	}
	
	return num;
}


/// The observation probability for population transition data
vector <double> StateSpecies::likelihood_obs_pop_trans(const vector <unsigned int> &list, double &like_ch)
{	
	vector <double> store;
	
	if(mode == SIM) return store;
	
	const auto &poptd = sp.pop_trans_data;
	const auto &poptf = sp.pop_trans_filter;

	for(auto i : list){
		const auto &ptd = poptd[i];
		const auto &ptf = poptf[ptd.ref];
		
		auto num = 0.0;	
		double va;
		for(auto tr : ptf.tr_nonzero){
			if(ptd.time_vari) va = obs_eqn_value[ptd.trans_obs_mod_ref[tr]];
			else va = obs_eqn_value[ptf.trans_obs_mod_ref[tr]];
			if(va != 0) num += va*pop_trans_data_tgl[i][tr];
		}  
	
		auto val = obs_mod_probability(num,ptd.type,ptd.value,ptd.obs_mod_val);
		
		like_ch += val - Li_obs_pop_trans[i];
		
		store.push_back(Li_obs_pop_trans[i]);
		store.push_back(pop_trans_data_num[i]);
		Li_obs_pop_trans[i] = val;
		pop_trans_data_num[i] = num;
	}
	
	return store;
}


/// Calculates the populations in each cgl at each popdata measurement
// and calculates the number of transitions at each poptransdata meaasurement
vector <double> StateSpecies::calculate_pop_data_cgl_trgl()
{
	vector <double> store;
	
	if(mode == SIM) return store;
	
	const auto &popd = sp.pop_data;
	if(type == INDIVIDUAL){
		for(auto i = 0u; i < popd.size(); i++){
			for(auto c = 0u; c < sp.comp_gl.size(); c++){
				store.push_back(pop_data_cgl[i][c]);
				pop_data_cgl[i][c] = 0;
			}
		}
	}
	
	const auto &poptd = sp.pop_trans_data;
	for(auto i = 0u; i < poptd.size(); i++){
		for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){		
			store.push_back(pop_trans_data_tgl[i][tr]);
			pop_trans_data_tgl[i][tr] = 0;
		}
	}
	
	switch(type){
	case POPULATION: case DETERMINISTIC:
		for(auto ti = 0u; ti < T; ti++){
			for(auto trg = 0u; trg < sp.tra_gl.size(); trg++){
				auto num = trans_num[trg][ti];
				if(num > 0){
					for(auto ref : sp.pop_trans_ref[ti][trg]){
						pop_trans_data_tgl[ref][trg] += num;
					}
				}
			}
		}
		break;
		
	case INDIVIDUAL:
		{
			auto N = individual.size();
			vector <unsigned int> indc(N), index(N,0);
			
			for(auto i = 0u; i < N; i++) indc[i] = UNSET;
			
			double t;
			for(auto j = 0u; j <= popd.size(); j++){
				if(j < popd.size()) t = popd[j].tdiv; else t = LARGE;
				
				for(auto i = 0u; i < N; i++){
					//const auto &ind = individual[i];
					auto c = indc[i];
					
					const auto &ev = individual[i].ev;
					while(index[i] < ev.size() && ev[index[i]].tdiv < t){
						const auto &e = ev[index[i]];
					
						if(e.type == M_TRANS_EV || e.type == NM_TRANS_EV){
							auto ti = get_ti(e.tdiv);
							auto tr_gl = e.tr_gl;
					
							for(auto ref : sp.pop_trans_ref[ti][tr_gl]){
								pop_trans_data_tgl[ref][tr_gl]++;
							}
						}
									
						c = e.c_after;
						index[i]++;
					}
					
					if(c != UNSET && j < popd.size()){
						pop_data_cgl[j][c]++;
					}
					
					indc[i] = c;
				}
			}
		}
		break;
	}

	if(false){
		for(auto i = 0u; i < poptd.size(); i++){
			for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
				cout << pop_trans_data_tgl[i][tr] <<" ,";
			}
			cout << "pp" << endl;
		}
	}
	
	return store;
}


/// Restore values of pop_data_cgl and pop_trans_data_tgl
void StateSpecies::restore_pop_data_cgl_trgl(const vector <double> &store)
{
	auto k = 0u;
	
	const auto &popd = sp.pop_data;
	if(type == INDIVIDUAL){
		for(auto i = 0u; i < popd.size(); i++){
			for(auto c = 0u; c < sp.comp_gl.size(); c++){
				pop_data_cgl[i][c] = store[k]; k++;
			}
		}
	}
	
	const auto &poptd = sp.pop_trans_data;
	for(auto i = 0u; i < poptd.size(); i++){
		for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
			pop_trans_data_tgl[i][tr] = store[k]; k++;
		}
	}
}


/// Gets the probability of a diagnostic test
double StateSpecies::like_diag_test(unsigned int c, const ObsData &ob) const
{
	if(c == UNSET) return LI_WRONG;
	else{
		const auto &om = sp.source[ob.so].obs_model;
		const auto &cgl = sp.comp_gl[c];
		
		auto cl = om.diag_test_sens.cl;
		if(om.diag_test_sens.comp[cgl.cla_comp[cl]] == true){ // Truely infected
			auto Se = obs_eqn_value[ob.Se_obs_eqn_ref];
			if(ob.test_res == true) return log(Se+LOG_THRESH); 
			else return log(1-Se+LOG_THRESH);
		}
		else{
			auto Sp = obs_eqn_value[ob.Sp_obs_eqn_ref];
			if(ob.test_res == false) return log(Sp+LOG_THRESH); 
			else return log(1-Sp+LOG_THRESH);
		}
	}
}

/// Gets the observation likelihood of transition probability observation
double StateSpecies::like_trans_prob_obs(unsigned int tr, const ObsData &ob) const
{
	auto val = get_trans_obs_prob(tr,ob);

	if(val == 1) return 0;
	
	if(val == 0) return LI_WRONG;
	else{
		if(val < 0) emsg("Compartment observational likelihood has become negative"); 
		return log(val);
	}						
}


/// Gets the observation likelihood of a compartmental observation
double StateSpecies::like_comp_obs(unsigned int c, const ObsData &ob) const
{
	if(ob.not_alive){ 
		if(c != UNSET) return LI_WRONG;
		return 0;
	}
	else{
		if(c == UNSET) return LI_WRONG;
		else{
			auto c_exact = ob.c_exact;
			if(c_exact != UNSET){
				if(sp.comp_gl[c].cla_comp[ob.cl] != c_exact) return LI_WRONG;
			}
			else{
				auto val = obs_eqn_value[ob.obs_eqn_ref[sp.comp_gl[c].cla_comp[ob.cl]]];
				if(val != 1){
					if(val == 0) return LI_WRONG;
					else{
						if(val < 0) emsg("Compartment observational likelihood has become negative"); 
						return log(val);
					}						
				}
			}
			return 0;
		}
	}
}


/// Calcualtes the probability of an individual entering a given state
double StateSpecies::like_ind_init_obs(unsigned int c, unsigned int entref) const
{
	if(c == UNSET) return LI_WRONG;
	
	const auto &ent = sp.enter[entref];

	if(ent.c_set != UNSET){
		if(c == ent.c_set) return 0;
		return LI_WRONG;
	}
	else{
		const auto &cla_comp = sp.comp_gl[c].cla_comp;
		
		auto prob = 1.0;
		for(auto cl = 0u; cl < sp.ncla; cl++){
			const auto &ent_cl = ent.cla[cl];
			auto cc = cla_comp[cl];
			if(ent_cl.c_set != UNSET){
				if(cc != ent_cl.c_set) return LI_WRONG;
			}
			else{
				auto val = obs_eqn_value[ent_cl.obs_eqn_ref[cc]];
				if(val < 0) emsg("Initial compartment observational likelihood has become negative"); 
				prob *= val;
			}
		}

		if(prob == 0) return LI_WRONG;
		return log(prob);
	}
}


/// Updates pop_trans and works out the change in observation likelihood
double StateSpecies::update_pop_trans_data(double t, unsigned int tr_gl, int sign)
{
	auto ti = get_ti(t);
	return update_pop_trans_data2(ti,tr_gl,sign);
}


/// Updates pop_trans and works out the change in observation likelihood
double StateSpecies::update_pop_trans_data2(unsigned int ti, unsigned int tr_gl, int sign)
{
	auto dLi = 0.0; 
	
	vector <double> back_vec, Li_vec;
	
	double val;
	for(auto ref : sp.pop_trans_ref[ti][tr_gl]){
		const auto &ptd = sp.pop_trans_data[ref];
		
		//if(ptd.eqn_zero[tr_gl] == true) emsg("Should not be zero");			
	
		auto dL = -obs_mod_probability(pop_trans_data_num[ref],ptd.type,ptd.value,ptd.obs_mod_val);
		//auto dL = -normal_probability(pop_trans_data_num[ref],ptd.value,ptd.sd);
	
		if(ptd.time_vari){
			val = sign*obs_eqn_value[ptd.trans_obs_mod_ref[tr_gl]];
		}
		else{
			const auto &ptf = sp.pop_trans_filter[ptd.ref];
			val = sign*obs_eqn_value[ptf.trans_obs_mod_ref[tr_gl]];
		}
		
		back_vec.push_back(val);
	
		pop_trans_data_tgl[ref][tr_gl] += sign;	
		pop_trans_data_num[ref] += val;
		
		dL += obs_mod_probability(pop_trans_data_num[ref],ptd.type,ptd.value,ptd.obs_mod_val);
	
		//dL += normal_probability(pop_trans_data_num[ref],ptd.value,ptd.sd);
		
		Li_obs_pop_trans[ref] += dL;
		Li_vec.push_back(dL);
		
		dLi += dL;
	}

	back.push_back(Back(POP_TRANS_DATA_NUM,ti,tr_gl,sign,back_vec));
	back.push_back(Back(LI_OBS_POP_TRANS,ti,tr_gl,Li_vec));

	return dLi;
}


/// Changes a population transition measurement
double StateSpecies::update_single_pop_trans_data(double dnum, unsigned int ref)
{
	const auto &ptd = sp.pop_trans_data[ref];
	
	//auto dL = -normal_probability(pop_trans_data_num[ref],ptd.value,ptd.sd);
	auto dL = -obs_mod_probability(pop_trans_data_num[ref],ptd.type,ptd.value,ptd.obs_mod_val);
	
	back.push_back(Back(POP_TRANS_DATA_NUM_SINGLE,ref,pop_trans_data_num[ref]));
	
	pop_trans_data_num[ref] += dnum;
	//dL += normal_probability(pop_trans_data_num[ref],ptd.value,ptd.sd);
	dL += obs_mod_probability(pop_trans_data_num[ref],ptd.type,ptd.value,ptd.obs_mod_val);
	
	back.push_back(Back(LI_OBS_POP_TRANS_SINGLE,ref,Li_obs_pop_trans[ref]));
	
	Li_obs_pop_trans[ref] += dL;
	
	return dL;
}


/// Changes a single population measurement
double StateSpecies::update_single_pop_data(double dnum, unsigned int ref)
{
	const auto &ptd = sp.pop_data[ref];
	
	//auto dL = -normal_probability(pop_data_num[ref],ptd.value,ptd.sd);
	auto dL = -obs_mod_probability(pop_data_num[ref],ptd.type,ptd.value,ptd.obs_mod_val);
	
	back.push_back(Back(POP_DATA_NUM,ref,pop_data_num[ref]));
	
	pop_data_num[ref] += dnum;
	//dL += normal_probability(pop_data_num[ref],ptd.value,ptd.sd);
	dL += obs_mod_probability(pop_data_num[ref],ptd.type,ptd.value,ptd.obs_mod_val);
	
	back.push_back(Back(LI_OBS_POP,ref,Li_obs_pop[ref]));
	
	Li_obs_pop[ref] += dL;
	
	return dL;
}


/// Updates pop_data based on a change in compartment for an individual
double StateSpecies::update_pop_data(unsigned int i, unsigned int c_old, unsigned int c_new)
{
	const auto &pd = sp.pop_data[i];
	const auto &pf = sp.pop_filter[pd.ref];
	
	auto num = pop_data_num[i];
	
	//auto dL = -normal_probability(num,pd.value,pd.sd);
	auto dL = -obs_mod_probability(num,pd.type,pd.value,pd.obs_mod_val);
	
	if(c_old != UNSET){
		back.push_back(Back(POP_DATA_CGL,i,c_old,-1));
		
		pop_data_cgl[i][c_old]--;
		if(pf.eqn_zero[c_old] == false){
			if(pd.time_vari == true) num -= obs_eqn_value[pd.comp_obs_mod_ref[c_old]];
			else num -= obs_eqn_value[pf.comp_obs_mod_ref[c_old]];
		}	
	}
	
	if(c_new != UNSET){
		back.push_back(Back(POP_DATA_CGL,i,c_new,1));
	
		pop_data_cgl[i][c_new]++;
		if(pf.eqn_zero[c_new] == false){
			if(pd.time_vari == true) num += obs_eqn_value[pd.comp_obs_mod_ref[c_new]];
			else num += obs_eqn_value[pf.comp_obs_mod_ref[c_new]];
		}
	}
	
	//dL += normal_probability(num,pd.value,pd.sd);
	dL += obs_mod_probability(num,pd.type,pd.value,pd.obs_mod_val);
	
	back.push_back(Back(POP_DATA_NUM,i,pop_data_num[i]));
	
	pop_data_num[i] = num;
	
	back.push_back(Back(LI_OBS_POP,i,Li_obs_pop[i]));
	
	Li_obs_pop[i] += dL;
	
	return dL;
}


/// Sets the observation probability for individual transition data
vector <double> StateSpecies::likelihood_unobs_trans(unsigned int e, const vector <unsigned int> &list, double &like_ch)
{
	vector <double> store;
	
	const auto &eq = eqn[sp.obs_trans_eqn[e]];
	auto &otev = obs_trans_eqn_value[e];
	const auto &oten = obs_trans_eqn_num[e];
	const auto &precalc = param_val.precalc;

	double log_val;
	if(eq.time_vari == false){
		auto val = 1 - eq.calculate_param(precalc);
		if(val < TINY) log_val = LI_WRONG;
		else{
			if(val > 1) emsg("Trans prob out of range");
			log_val = log(val);
		}
		
		for(auto ti : list){
			auto log_val_old = otev[ti];
			store.push_back(log_val_old);
			like_ch += (log_val-log_val_old)*oten[ti];
			otev[ti] = log_val;
		}
	}
	else{
		auto num=0u;
		for(auto ti : list){
			auto val = 1 - eq.calculate_no_popnum(ti,precalc);	
		
			if(val < TINY) log_val = LI_WRONG;
			else{
				if(val > 1) emsg("Trans prob out of range");
				log_val = log(val);
			}
			
			auto log_val_old = otev[ti];
			store.push_back(log_val_old);
			like_ch += (log_val-log_val_old)*oten[ti];
			num+=oten[ti];
			otev[ti] = log_val;
		}
	}

	return store;
}


/// Sets all observed transition values
void StateSpecies::set_all_obs_trans_value()
{
	auto all_time = seq_vec(T);
	for(auto m = 0u; m < sp.obs_trans_eqn.size(); m++){
		auto temp = 0.0;
		likelihood_unobs_trans(m,all_time,temp);
	}
}
				

/// Alters obs_trans_eqn_num based on a change in event
void StateSpecies::alter_obs_trans_eqn(unsigned int tr_gl, double t, int sign, bool observed, const vector <ObsData> &obs_data, double &like_ch)
{							
	auto ti = get_ti(t);
	const auto &oter = sp.obs_trans_eqn_ref[tr_gl][ti];
		
	if(observed == false){
		vector <double> vec;
		for(auto m : oter){
			like_ch += sign*obs_trans_eqn_value[m][ti];
			vec.push_back(sign);
			obs_trans_eqn_num[m][ti] += sign;
		}
		back.push_back(Back(OBS_TRANS_EQN_NUM,tr_gl,ti,vec));
	}
	else{
		// Checks for the case when one transition data is observed and the other not
		if(oter.size() > 1){
			vector <double> vec;
			for(auto m : oter){
				auto fl = false;
				auto eq = sp.obs_trans_eqn[m];
	
				for(const auto &od : obs_data){
					auto tt = od.tdiv;
					if(tt == t){
						switch(od.type){
						case OBS_SOURCE_EV: case OBS_TRANS_EV: case OBS_SINK_EV:
							{
								const auto &ot = sp.obs_trans[od.ref];
								if(eq == ot.tra_prob_eqn[tr_gl].eq_ref) fl = true;
							}
							break;
						default: break;
						}
						
						if(fl == true) break;
					}
					else{
						if(tt > t) break;
					}
				}
			
				if(fl == false){
					like_ch += sign*obs_trans_eqn_value[m][ti];
					vec.push_back(sign);
					obs_trans_eqn_num[m][ti] += sign;
				}
				else{
					vec.push_back(0);
				}
			}
			back.push_back(Back(OBS_TRANS_EQN_NUM,tr_gl,ti,vec));
		}
	}
}
		

/// Calculates number of not observed transition in obs_trans_eqn (for diagnostic purposes)	
void StateSpecies::calculate_obs_trans_eqn_num()
{
	
	for(auto m = 0u; m < sp.obs_trans_eqn.size(); m++){
		for(auto ti = 0u; ti < T; ti++) obs_trans_eqn_num[m][ti] = 0;
	}
	
	for(auto i = 0u; i < individual.size(); i++){
		const auto &ind = individual[i];
		for(const auto &ev : ind.ev){
			if(ev.type == M_TRANS_EV || ev.type == NM_TRANS_EV){
				auto temp = 0.0;
				alter_obs_trans_eqn(ev.tr_gl,ev.tdiv,1,ev.observed,sp.individual[i].obs,temp);
			}
		}
	}
}


/// Works out the log of the probability of not observing a transition
double StateSpecies::trans_logprob_notobs(unsigned int trg, unsigned int ti) const 
{
	auto sum = 0.0;
	for(auto e : sp.obs_trans_eqn_ref[trg][ti]){
		sum += obs_trans_eqn_value[e][ti];
	}
	
	return sum;
}


/// Sets which events are observed (used when calculating transition probability)
void StateSpecies::set_event_observed()
{
	for(auto &ind : individual){
		for(auto &ev : ind.ev) ev.observed = false;
	}
		
	auto num = 0u;
	for(auto i = 0u; i < sp.nindividual_in; i++){
		auto &ind = individual[i];
		const auto &indi = sp.individual[i];
		for(auto &ev : ind.ev){
			switch(ev.type){
			case ENTER_EV: case LEAVE_EV: case MOVE_EV: ev.observed = true; break;
			default: break;
			}
		}		
		
		for(auto &ob : indi.obs){
			switch(ob.type){
			case OBS_TRANS_EV: case OBS_SOURCE_EV: case OBS_SINK_EV: 
				{
					for(auto &ev : ind.ev){
						if(ev.tdiv == ob.tdiv && (ev.type == M_TRANS_EV || ev.type == NM_TRANS_EV)){	
							ev.observed = true;
							num++;
						}						
					}
				}
				break;

			default:
				break;
			}
		}
	}
	
	if(testing){
		auto num_obs = 0u;
		auto i_sel = UNSET;
		for(auto i = 0u; i < sp.nindividual_in; i++){
			const auto &indi = sp.individual[i];
			const auto &ind = individual[i];
			
			auto numi = 0u;
			for(auto &ob : indi.obs){
				switch(ob.type){
				case OBS_TRANS_EV: case OBS_SOURCE_EV: case OBS_SINK_EV: 
					for(auto &ev : ind.ev){
						if(ev.tdiv == ob.tdiv && (ev.type == M_TRANS_EV || ev.type == NM_TRANS_EV)){	
							numi++;
						}
					}
					break;
				default: break;
				}
			}
			
			auto numob = 0u;
			for(auto &ob : indi.obs){
				switch(ob.type){
				case OBS_TRANS_EV: case OBS_SOURCE_EV: case OBS_SINK_EV:
					num_obs++;
					numob++;
					break;
				default: break;
				}
			}
			
			if(numi != numob) i_sel = i;
		}
		
		if(num != num_obs && false){
			cout << num << " / " << num_obs << "NUM OBS" << " Ind:" << sp.individual[i_sel].name << endl;
		}
	}
}

// This deals with individual-based updates where there are no observations

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
using namespace std;

#include "ind_ev_sampler.hh"
#include "utils.hh"

Event IndEvSampler::resample_init_event(unsigned int i, double &probif) const
{
	const auto &ind = individual[i];
	
	auto e_init = ind.ev[0];
	e_init.inf_node_ref = UNSET;
	e_init.ind_inf_from = IndInfFrom();
	
	if(!ind.init_c_set){
		switch(e_init.type){
		case ENTER_EV: 
			e_init.c_after = init_comp_sampler(i,probif); 
			break;
			
		case M_TRANS_EV: 
			{
				auto sampler = source_sampler(i,e_init);
				auto z = ran();
				auto k = 0u; while(k < sampler.size() && z > sampler[k].prob_sum) k++;
				
				probif += log(sampler[k].prob);
				
				e_init.tr_gl = k;
				e_init.c_after = sp.tra_gl[k].f; 
			}
			break;
			
		default: emsg("OP prob"); break;
		}
	}
	
	return e_init;
}


/// Finds the probability of selecting the initial state
double IndEvSampler::resample_init_event_prob(unsigned int i, const Event &e_init) const 
{	
	const auto &ind = individual[i];
	
	if(ind.init_c_set) return 0;
	
	switch(e_init.type){
	case ENTER_EV:
		return init_comp_sampler_prob(e_init.c_after,i);
	
	case M_TRANS_EV: 
		{
			auto sampler = source_sampler(i,e_init);
			return log(sampler[e_init.tr_gl].prob);
		}
		
	default: emsg("OP prob"); return 0; 
	}
}


/// Creates a sampler for the initial state of the individual
unsigned int IndEvSampler::init_comp_sampler(unsigned int i, double &probif) const 
{
	const auto &indi = sp.individual[i];
		
	const auto &ent = sp.enter[indi.enter_ref];
	
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
			
			probif += log(prob[c]/sum);
			c_after += sp.comp_mult[cl]*c;
		}
	}
	
	return c_after;
}	
	
	
/// Gets the probability of an initial compartmental sampler
double IndEvSampler::init_comp_sampler_prob(unsigned int c_after, unsigned int i) const 
{
	const auto &indi = sp.individual[i];
		
	const auto &ent = sp.enter[indi.enter_ref];
	
	if(ent.c_set != UNSET) return 0;
	
	const auto &cla_comp = sp.comp_gl[c_after].cla_comp;
	
	auto probif = 0.0;
	for(auto cl = 0u; cl < sp.ncla; cl++){
		const auto &ent_cl = ent.cla[cl];
		if(ent_cl.c_set == UNSET){
			auto sum = 0u; 
			const auto &claa = sp.cla[cl];
			const auto &oer = ent_cl.obs_eqn_ref;
			auto C = claa.ncomp;
			
			auto cc = cla_comp[cl];
			
			double pr_st = UNSET;
			for(auto c = 0u; c < C; c++){
				auto pr = obs_eqn_value[oer[c]];
				sum += pr;
				if(c == cc) pr_st = pr;
			}
				
			probif += log(pr_st/sum);
		}
	}
	
	return probif;
}
	
	
/// Creates a sampler for the initial state of the individual
vector <InitStateProb> IndEvSampler::source_sampler(unsigned int i, const Event &ev) const 
{
	vector <InitStateProb> init_state;
	
	auto t = ev.t;
	auto ti = sp.get_ti(t); 
	auto N = sp.tra_gl.size();

	auto prob_sum = 0.0;
	if(ev.observed){ // If the event is observed
		const auto &ob = sp.individual[i].obs[0];
		if(ob.t != t) emsg("Should have the same time");
		
		for(auto tr = 0u; tr < N; tr++){	
			auto prob = 0.0;

			const auto &tra = sp.tra_gl[tr];
			if(tra.variety == SOURCE_TRANS){
				const auto &mev = markov_eqn_vari[tra.markov_eqn_ref];
				auto tii = ti; if(mev.time_vari == false) tii = 0;
			
				prob = mev.div[tii].value;
				
				prob *= get_trans_obs_prob(tr,ob);
				if(prob < 0) emsg("rate is negative");
				prob_sum += prob;
			}
		
			InitStateProb isp; isp.prob = prob; isp.prob_sum = prob_sum;
			init_state.push_back(isp);
		}	
	}
	else{
		for(auto tr = 0u; tr < N; tr++){	
			auto prob = 0.0;

			const auto &tra = sp.tra_gl[tr];
			if(tra.variety == SOURCE_TRANS){
				const auto &mev = markov_eqn_vari[tra.markov_eqn_ref];
				auto tii = ti; if(mev.time_vari == false) tii = 0;
				prob = mev.div[tii].value;
				if(prob < 0) emsg("rate is negative");
				prob_sum += prob;
			}
			
			InitStateProb isp; isp.prob = prob; isp.prob_sum = prob_sum;
			init_state.push_back(isp);
		}		
	}

	for(auto &isp : init_state){
		isp.prob /= prob_sum; isp.prob_sum /= prob_sum;
	}
		
	return init_state;
}

// Allows for samples to be sampled from the source

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
using namespace std;

#include "source_sampler.hh"
#include "utils.hh"
#include "matrix.hh"

SourceSampler::SourceSampler(const vector <MarkovEqn> &markov_eqn, const vector <TransGlobal> &tra_gl, const vector <CompGlobal> &comp_gl, const Details &details, const InitCond &init_cond) : markov_eqn(markov_eqn), tra_gl(tra_gl), comp_gl(comp_gl), details(details), init_cond(init_cond)
{};
	
	
/// Used to generate a sampler for adding new individuals
void SourceSampler::update(unsigned int nind_obs, vector <Individual> &individual, const vector <MarkovEqnVariation> &markov_eqn_vari, bool contains_source)
{
	auto T = details.T;
	auto I = individual.size();
	auto C = comp_gl.size();
		
	if(rate.size() == 0){
		memory = 1.0 - 500.0/(details.sample*details.burnin_frac);
		if(memory < 0.5) memory = 0.5;
		
		enter_frac = 0.5;
		update_num = 10;
		
		me_map.resize(markov_eqn.size(),UNSET);
	
		for(auto e = 0u; e < markov_eqn.size(); e++){
			const auto &me = markov_eqn[e];
			if(me.source == true){
				me_map[e] = me_list.size();
				me_list.push_back(e);
			}
		}
		
		enter_prob.resize(C,1.0/C);
		enter_prob_sum.resize(C);
		
		rate.resize(T); rate_sum.resize(T);
		for(auto ti = 0u; ti < T; ti++){
			rate[ti].resize(me_list.size(),0);
			rate_sum[ti].resize(me_list.size());
		}
		
		margin_sum.resize(T);
	}
	
	auto N = me_list.size();
	
	// Updates the fraction of individuals which enter
	{
		auto val = 0.0; 
		for(auto i = nind_obs; i < I; i++){
			const auto &e_init = individual[i].ev[0];
			if(e_init.type == ENTER_EV) val++;
		}
		
		enter_frac = memory*enter_frac + (1-memory)*(val/(I-nind_obs+0.5));
		if(enter_frac < 0.1) enter_frac = 0.1;
		
		if(init_cond.type != INIT_POP_DIST) enter_frac = 0;
		else{
			if(!contains_source) enter_frac = 1;
		}
	}
	
	// Updates the number of updates per proposal
	{
		update_num = memory*update_num + (1-memory)*(I-nind_obs);
		if(update_num < 1) update_num = 1;
	}
	
	// Updates sampler for compartment being entered
	if(init_cond.type == INIT_POP_DIST){
		vector <double> frac_unobs(C,0);
		auto ntot = 0u;
		for(auto i = nind_obs; i < individual.size(); i++){
			const auto &ind = individual[i];
			if(ind.ev.size() > 0){
				if(ind.ev[0].type == ENTER_EV && ind.ev[0].tdiv == 0){
					frac_unobs[ind.ev[0].c_after]++;
					ntot++;
				}
			}
		}
		
		if(ntot > 0){
			auto sum = 0.0;
			for(auto c = 0u; c < C; c++){
				enter_prob[c] = memory*enter_prob[c] + (1-memory)*(frac_unobs[c]/ntot);
				sum += enter_prob[c];
				enter_prob_sum[c] = sum;
			}
			
			for(auto c = 0u; c < C; c++){
				enter_prob[c] /= sum;
				enter_prob_sum[c] /= sum;
			}
		}
	}
		
	// Updates sampler for source
	if(contains_source){
		for(auto i = 0u; i < N; i++){
			auto e = me_list[i];
			const auto &me = markov_eqn[me_list[i]];
			const auto &mev = markov_eqn_vari[e];
			
			auto fac = me.source_tr_gl.size();
			for(auto ti = 0u; ti < T; ti++){
				auto tii = ti; if(mev.time_vari == false) tii = 0;
				auto val = mev.div[tii].value+TINY;
				rate[ti][i] = rate[ti][i]*memory + fac*val;
			}
		}
		
		auto sum = 0.0;
		for(auto ti = 0u; ti < T; ti++){
			auto sum2 = 0.0;
			for(auto i = 0u; i < N; i++){
				sum2 += rate[ti][i];
				rate_sum[ti][i] = sum2;
			}
			if(sum2 > 0){
				for(auto i = 0u; i < N; i++) rate_sum[ti][i] /= sum2;
			}
			
			sum += sum2;
			margin_sum[ti] = sum;
		}
		rate_total = sum;
		if(sum == 0) emsg("zero rate"); 
		
		for(auto ti = 0u; ti < T; ti++){
			margin_sum[ti] /= sum;
		}
	}
}


/// Samples an initial source 
Event SourceSampler::sample(double &probif) const
{
	Event ev;
	ev.move_c = UNSET;
	ev.cl = UNSET;
	ev.observed = false;
	
	auto T = details.T;
	auto N = me_list.size();
	
	if(ran() < enter_frac){   // An enter event
		auto C = comp_gl.size();
	
		ev.type = ENTER_EV;
		auto z = ran();
		auto c = 0u; while(c < C && z > enter_prob_sum[c]) c++;
		if(c == C) emsg("samp prob");
		 
		ev.c_after = c;
		ev.tr_gl = UNSET;
		ev.tdiv = 0;
	
		auto val = enter_frac*enter_prob[c];
		if(val == 0) probif = -LARGE;
		else probif += log(val);
	}
	else{
		auto z = ran();
		auto ti = 0u; while(ti < T && z > margin_sum[ti]) ti++;
		if(ti == T) emsg("Out of range");

		const auto &rs = rate_sum[ti];
		z = ran();
		auto k = 0u; while(k < N && z > rs[k]) k++;
		if(k == N) emsg("select prob");
		
		auto &me = markov_eqn[me_list[k]];
		auto S = me.source_tr_gl.size();
		if(S == 0) emsg("zero S");
		
		auto s = 0u; if(S > 1) s = (unsigned int)(ran()*S);
		
		auto tr_gl = me.source_tr_gl[s];
		const auto &tr = tra_gl[tr_gl];
		
		ev.type = M_TRANS_EV;
		ev.c_after = tr.f;
		ev.cl = tr.cl;
		ev.tr_gl = tr_gl;
		ev.tdiv = ti+ran();
		
		auto val = (1-enter_frac)*rate[ti][k]/(S*rate_total);
		if(val == 0) probif = -LARGE;
		else probif += log(val);
	}
	ev.inf_node_ref = UNSET;
	ev.ind_inf_from = IndInfFrom();
											
	return ev;
}


/// Calculate the sample probability
double SourceSampler::sample_prob(const Event &ev) const
{
	switch(ev.type){
	case ENTER_EV: 
		if(enter_frac == 0) emsg("KK");
		return log(enter_frac*enter_prob[ev.c_after]);
	
	case M_TRANS_EV:
		{
			auto ti = get_ti(ev.tdiv);
			auto e = tra_gl[ev.tr_gl].markov_eqn_ref;
			auto &me = markov_eqn[e];
			auto S = me.source_tr_gl.size();
			return log((1-enter_frac)*rate[ti][me_map[e]]/(S*rate_total));
		}
	
	default:
		emsg("Op prob");
		return 0;
	}
}


/// Prints the sampler
string SourceSampler::print() const 
{
	auto T = details.T;
	auto C = comp_gl.size();
	auto N = me_list.size();
	
	stringstream ss;
	
	ss << "Source Sampler:" << endl;
	if(rate.size() == 0){ ss << "Not used" << endl; return ss.str();}
	
	ss << "# update: " << update_num << endl;
	ss << "Enter fraction: " << enter_frac << endl;
	
	ss << "Enter prob:" << endl;
	for(auto c = 0u; c < C; c++){
		ss << comp_gl[c].name << " " << enter_prob[c] << endl;
	}
	
	ss << "Source prob:" << endl;
	for(auto ti = 0u; ti < T; ti++){
		ss << ti << ": ";
		for(auto i = 0u; i < N; i++){
			auto &me = markov_eqn[me_list[i]];
			for(auto tr_gl : me.source_tr_gl) ss << tra_gl[tr_gl].name << " ";
			ss << " " << rate[ti][i]/rate_total << ", ";
		}
		ss << endl;
	}
	
	return ss.str();
}

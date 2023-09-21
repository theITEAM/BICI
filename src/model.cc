// This file defines the compartmental model 

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
//#include "math.h"

using namespace std;

#include "model.hh"
#include "utils.hh"

/// Initialises the model 
Model::Model()
{
	mode = MODE_UNSET;  
	nspecies = 0;
};


/// Adds an equation reference to an EquationInfo
void Model::add_eq_ref(EquationInfo &eqi)
{
	//auto val = number(eqi.te);

	//if(val == UNSET){
	auto e = 0u;
	auto emax = eqn.size();
	while(e < emax){
		const auto &eq = eqn[e]; 
		if(eq.te_init == eqi.te && eq.type == eqi.type && eq.sp_p == eqi.p && eq.sp_cl == eqi.cl) break;
		e++;
	}
	
	if(e < emax) eqi.eq_ref = e;
	else{
		eqi.eq_ref = emax;		
		
		Equation eq(eqi.te,eqi.type,eqi.p,eqi.cl,eqi.line_num,species_simp,param,spline,param_vec,pop,timepoint);
		eqn.push_back(eq);
	}
}


/// Samples a set of parameter values from the model
vector <double> Model::param_sample() const
{
	vector <double> param_val(nparam_vec);
	
	for(auto th = 0u; th < nparam_vec; th++){
		const auto &pv = param_vec[th];
		
		const auto &par = param[pv.th];
		switch(par.variety){
		case CONST_PARAM:
			param_val[th] = par.value[pv.index].value;
			break;
			
		case REPARAM_PARAM:
			{
				auto eq_ref = par.value[pv.index].eq_ref;
				if(eq_ref == UNSET) emsg("eq_ref should be set");
				param_val[th] = eqn[eq_ref].calculate_param_only(param_val);
			}
			break;
			
		case DIST_PARAM: case PRIOR_PARAM:
			param_val[th] = prior_sample(pv.prior,param_val);
			break;
			
		case UNSET_PARAM: emsg("error param"); break;
		}
	}
	
	if(mode == SIM) set_branch_auto(param_val);
	
	return param_val;
}


/// Determines if a set of parameters is within the bounds of the prior
bool Model::inbounds(const vector <double> &param_val) const
{
	auto total = prior_total(param_val);
	if(total < -LARGE/2) return false;
	return true;
}


/// The total prior for all parameters
double Model::prior_total(const vector <double> &param_val) const 
{
	auto prior = prior_prob(param_val);
	return sum(prior);
}


/// Calculate the prior for the parameters
vector <double> Model::prior_prob(const vector <double> &param_val) const 
{
	vector <double> prior(nparam_vec,0);
	
	for(auto th = 0u; th < nparam_vec; th++){
		const auto &pv = param_vec[th];
		if(pv.variety == PRIOR_PARAM){
			prior[th] = prior_probability(param_val[th],pv.prior,param_val);
		}
	}
	
	return prior;
}


/// Calculate the prior for the parameters
vector <double> Model::dist_prob(const vector <double> &param_val) const
{
	vector <double> dist(nparam_vec,0);
	
	for(auto th = 0u; th < nparam_vec; th++){
		const auto &pv = param_vec[th];
		if(pv.variety == DIST_PARAM){
			dist[th] = prior_probability(param_val[th],pv.prior,param_val);
		}
	}
	
	return dist;
}


/// Recalculates the prior for a given parameter
double Model::recalculate_prior(unsigned int th, vector <double> &prior_prob, const vector <double> &param_val, double &like_ch) const
{
	auto store = prior_prob[th];
	const auto &pv = param_vec[th];
	prior_prob[th] = prior_probability(param_val[th],pv.prior,param_val);
	like_ch += prior_prob[th]-store;
	
	return store;
}


/// Calculates the prior for the spline
vector <double> Model::spline_prior(const vector <double> &param_val) const
{
	vector <double> spline_prior;
	
	for(auto s = 0u; s < spline.size(); s++){
		const auto &spl = spline[s];
		
		auto Li = 0.0;
		if(spl.constant == false && spl.info.smooth == true){
			if(spl.info.on != true) emsg("spline should be on");
			
			auto type = spl.info.smooth_type;
			auto sd = spl.info.smooth_value;
			auto cv = sqrt(exp(sd*sd)-1);
			
			for(auto t = 1u; t < spl.param_ref.size(); t++){
				auto val = param_val[spl.param_ref[t]];
				auto last = param_val[spl.param_ref[t-1]];
				
				switch(type){
				case LOG_NORMAL_SMOOTH: Li += lognormal_probability(val,last,cv); break;
				case NORMAL_SMOOTH: Li += normal_probability(val,last,sd); break;
				}	
			}
		}
		spline_prior.push_back(Li);
	}
	
	return spline_prior;
}


/// Recalculates the spline prior for a given parameter
double Model::recalculate_spline_prior(unsigned int s, vector <double> &spline_prior, const vector <double> &param_val, double &like_ch) const
{
	const auto &spl = spline[s];
	
	auto type = spl.info.smooth_type;
	auto sd = spl.info.smooth_value;
	auto cv = sqrt(exp(sd*sd)-1);
	
	auto Li = 0.0;
	for(auto t = 1u; t < spl.param_ref.size(); t++){
		auto val = param_val[spl.param_ref[t]];
		auto last = param_val[spl.param_ref[t-1]];
		
		switch(type){
		case LOG_NORMAL_SMOOTH: Li += lognormal_probability(val,last,cv); break;
		case NORMAL_SMOOTH: Li += normal_probability(val,last,sd); break;
		}	
	}
	
	auto store = spline_prior[s];
	like_ch += Li-store;
	spline_prior[s] = Li;
	
	return store;
}
	

/// Calculates pop from cpop
vector <double> Model::calculate_popnum(vector <StateSpecies> &state_species) const
{
	vector <double> popnum(pop.size(),0);
	
	for(auto i = 0u; i < pop.size(); i++){     // Calculates population based on cpop
		const auto &po = pop[i];
		const auto &ssp = state_species[po.sp_p];
		if(ssp.type == POPULATION){
			auto sum = 0.0;
			for(auto &te : po.term) sum += ssp.cpop[te.c]*te.w;
			popnum[i] = sum;
		}
	}

	for(auto p = 0u; p < nspecies; p++){
		const auto &ssp = state_species[p];
		if(ssp.type == INDIVIDUAL){
			const auto &sp = species[p];
			for(const auto &ind : ssp.individual){
				auto c = ind.c;
				if(c < CODE){
					for(const auto &pr : sp.comp_gl[c].pop_ref){
						const auto &po = pop[pr.po];
				
						auto num = 1.0; 
						for(auto ie : po.ind_eff_mult) num *= ind.exp_ie[ie];
						for(auto fe : po.fix_eff_mult) num *= ind.exp_fe[fe];
							
						if(po.term[pr.index].c != c) emsg("Problem");
						popnum[pr.po] += po.term[pr.index].w*num;
					}
				}
			}
		}
	}
		
	return popnum;
}

/// Recalculates popnum_t for populations identified by list and individuals in ssp 
// Used for diagnostics
vector < vector <double> > Model::calculate_popnum_t(vector <StateSpecies> &state_species) const
{
	vector < vector <double> > popnum_t;

	auto T = ntimepoint-1;
	
	for(auto ti = 0u; ti < T; ti++){
		auto t = timepoint[ti];
		for(auto p = 0u; p < nspecies; p++){
			auto &ssp = state_species[p];
			switch(ssp.type){
			case INDIVIDUAL:
				{
					for(auto &ind : ssp.individual){
						auto c = ind.cinit;
						auto e = 0u; while(e < ind.ev.size() && ind.ev[e].t < t){ c = ind.ev[e].c_after; e++;}
						ind.c = c;
					}
				}
				break;
				
			case POPULATION:	
				ssp.cpop = ssp.cpop_st[ti];
				break;
			}
		}

		popnum_t.push_back(calculate_popnum(state_species));	
	}
	
	return popnum_t;
}


/// Recalculates popnum_t for populations identified by list and individuals in ssp 
vector <double> Model::recalculate_population(vector < vector <double> > &popnum_t, const vector <unsigned int> &list, const vector <StateSpecies> &state_species) const
{
	vector <double> store;
	
	if(pop.size() == 0) return store;
	
	auto T = ntimepoint-1;
	
	for(auto k : list){
		for(auto t = 0u; t < T; t++){
			store.push_back(popnum_t[t][k]);
			popnum_t[t][k] = 0;
		}
	}
	
	// This map shows which populations need to be recalulated
	vector <bool> map(pop.size(),false); for(auto k : list) map[k] = true;
	
	// This stores individual factors for different populations
	vector <double> inffac(pop.size(),UNSET);	

	for(auto p = 0u; p < nspecies; p++){
		const auto &ssp = state_species[p];
		if(ssp.type == INDIVIDUAL){
			const auto &sp = species[p];
		
			// Goes through each individual and works out how pop is updated				
			for(const auto &ind : ssp.individual){
				auto ti = 0u;
		
				vector <unsigned int> indfac_calc;
		
				auto c = ind.cinit;	
				for(auto e = 0u; e <= ind.ev.size(); e++){
					double t;
					if(e < ind.ev.size()) t = ind.ev[e].t;
					else t = timepoint[T];
						
					while(timepoint[ti] < t){
						if(c < CODE){						
							for(const auto &pr : sp.comp_gl[c].pop_ref){
								auto k = pr.po;
								if(map[k] == true){
									auto num = inffac[k];
									if(inffac[k] == UNSET){
										const auto &po = pop[k];
										
										num = 1.0; 
										for(auto ie : po.ind_eff_mult) num *=  ind.exp_ie[ie];
										for(auto fe : po.fix_eff_mult) num *= ind.exp_fe[fe];
							
										indfac_calc.push_back(k);
										inffac[k] = num;
									}
								
									popnum_t[ti][k] += num; 
								}
							}
						}
						ti++;			
					}
					
					if(e < ind.ev.size()) c = ind.ev[e].c_after;
				}
				
				if(ti != T) emsg("Problem with ti");
				
				for(auto k : indfac_calc) inffac[k] = UNSET;
			}
		}
	}
	
	return store;
}


/// Recalculates popnum_t for populations identified by list and individuals in ssp 
void Model::recalculate_population_restore(vector < vector <double> > &popnum_t, const vector <unsigned int> &list, const vector <double> &vec) const
{
	auto T = ntimepoint-1;
	
	auto j = 0u;
	for(auto k : list){
		for(auto t = 0u; t < T; t++){
			popnum_t[t][k] = vec[j]; j++; 
		}
	}
}



/// Creates species_simp
void Model::create_species_simp()
{
	for(auto &sp :species){
		SpeciesSimp ss(sp.name,sp.cla,sp.ind_effect,sp.fix_effect,sp.comp_mult,sp.comp_gl);
		species_simp.push_back(ss);
	}
}


/// Looks to set branching probabiltit
void Model::set_branch_auto(vector <double> &param_val) const
{
	for(auto &bpar : branch_param){
		vector <unsigned int> list;
		
		auto sum = 0.0;
		for(auto i = 0u; i < bpar.group.size(); i++){
			auto th = bpar.group[i];
			if(param_val[th] == UNSET) list.push_back(th);
			else sum += param_val[th];
		}
		
		if(list.size() != 1) emsg("Should have one unset");
		auto val = 1-sum;
	
		if(val < 0 || val > 1) emsg("BP out of range3");
	
		param_val[list[0]] = val;
	}
}


/// Samples a value from a prior distribution
double Model::prior_sample(const Prior &pri, const vector <double> &param_val) const
{				
	switch(pri.type){
	case UNIFORM_PR:
		{	
			auto min = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			auto max = eqn[pri.dist_param[1].eq_ref].calculate_param_only(param_val);
			if(min > max) emsg("Ordering of uniform prior is wrong");
			return min+ran()*(max-min);
		}
		break;
		
	case EXP_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			return exp_mean_sample(mean);
		}
		break;
	
	case NORMAL_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			auto sd = eqn[pri.dist_param[1].eq_ref].calculate_param_only(param_val);
			return normal_sample(mean,sd);
		}
		break;
		
	case GAMMA_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			auto cv = eqn[pri.dist_param[1].eq_ref].calculate_param_only(param_val);
			return gamma_sample(mean,cv);
		}
		break;
	
	case LOG_NORMAL_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			auto cv = eqn[pri.dist_param[1].eq_ref].calculate_param_only(param_val);
			return lognormal_sample(mean,cv);
		}
		break;
		
	case BETA_PR:
		{	
			auto alpha = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			auto beta = eqn[pri.dist_param[1].eq_ref].calculate_param_only(param_val);
			return beta_sample(alpha,beta);
		}
		break;
	
	case BERNOULLI_PR:
		{	
			auto z = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			return bernoulli_sample(z);
		}
		break;
		
	case FIX_PR:
		{
			auto val = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			return val;
		}
		break;
		
	case FLAT_PR:
		{
			return gamma_alpha_sample(1);
		}
		break;
		
	case DIRICHLET_PR:
		{
			auto alpha = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			return gamma_alpha_sample(alpha);
		}
		break;
	}
	
	return UNSET;
}


/// The log probability of sampling from a distribution
double Model::prior_probability(double x, const Prior &pri, const vector <double> &param_val) const
{				
	switch(pri.type){
	case UNIFORM_PR:
		{	
			auto min = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			auto max = eqn[pri.dist_param[1].eq_ref].calculate_param_only(param_val);
			if(x < min || x > max) return -LARGE;
			return 1.0/(max-min);
		}
		break;
		
	case EXP_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			if(x < 0 || mean <= 0) return -LARGE;
			return exp_mean_probability(x,mean);
		}
		break;
	
	case NORMAL_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			auto sd = eqn[pri.dist_param[1].eq_ref].calculate_param_only(param_val);
			if(sd <= 0) return -LARGE;
			return normal_probability(x,mean,sd);
		}
		break;
		
	case GAMMA_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			auto cv = eqn[pri.dist_param[1].eq_ref].calculate_param_only(param_val);
			if(x <= 0 || mean <= 0 || cv <= 0) return -LARGE;
			return gamma_probability(x,mean,cv);
		}
		break;
	
	case LOG_NORMAL_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			auto cv = eqn[pri.dist_param[1].eq_ref].calculate_param_only(param_val);
			if(x <= 0 || mean <= 0 || cv <= 0) return -LARGE;
			return lognormal_probability(x,mean,cv);
		}
		break;
		
	case BETA_PR:
		{	
			auto alpha = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			auto beta = eqn[pri.dist_param[1].eq_ref].calculate_param_only(param_val);
			if(x <= 0 || x >= 1 || alpha <= 0 || beta <= 0) return -LARGE;
			return beta_probability(x,alpha,beta);
		}
		break;
	
	case BERNOULLI_PR:
		{	
			auto z = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			if(x != 0 && x != 1) return -LARGE;
			if(z < 0 || z > 1) return -LARGE;
			return bernoulli_probability(x,z);
		}
		break;
		
	case FIX_PR:
		{
			auto val = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			if(x != val) return -LARGE;
			return 0;
		}
		break;
		
	case FLAT_PR:
		{
			if(x <= 0) return -LARGE;
			return gamma_alpha_probability(x,1);
		}
		break;
		
	case DIRICHLET_PR:
		{
			auto alpha = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			if(x <= 0 || alpha <= 0) return -LARGE; 
			return gamma_alpha_probability(x,alpha);
		}
		break;
	}
	
	return UNSET;
}


// Used to order affects
bool AL_ord(AffectLike al1, AffectLike al2)                      
{ return (al1.order_num < al2.order_num); };  


/// Orders a list of AffectLike based on the priority in the calculation
void Model::order_affect(vector <AffectLike> &vec) const
{
	for(auto &al : vec){
		switch(al.type){
			case SPLINE_PRIOR_AFFECT: case PRIOR_AFFECT: case DIST_AFFECT: case OMEGA_AFFECT:
			al.order_num = 0; 
			break;
			
		case SPLINE_AFFECT: case EXP_FE_AFFECT: case POP_AFFECT: case EXP_IE_AFFECT:
			al.order_num = 1; 
			break;
		
		case INDFAC_INT_AFFECT: 
			al.order_num = 2; 
			break;
		
		case DIV_VALUE_AFFECT:
			al.order_num = 3; 
			break;
			
		case MARKOV_LIKE_AFFECT: case NM_TRANS_AFFECT: case LIKE_IE_AFFECT: 
		case MARKOV_POP_AFFECT: case LIKE_OBS_AFFECT:
			al.order_num = 4; 
			break;
		}
	}
	
	sort(vec.begin(),vec.end(),AL_ord);
}

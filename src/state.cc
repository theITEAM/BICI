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
#include "matrix.hh"
#include "ind_ev_sampler.hh"

/// Initialises the model 
State::State(const Model &model) : model(model)
{
	T = model.ntimepoint-1;
	
	check_timer = 0;
	
	dtimepoint.resize(T);
	for(auto ti = 0u; ti < T; ti++) dtimepoint[ti] = model.timepoint[ti+1] - model.timepoint[ti];
	
	spline_val.clear();
	spline_val.resize(model.spline.size());
	for(auto &spl : spline_val) spl.val.resize(T,UNSET);
	
	popnum_t.resize(T);
	
	nspecies = model.nspecies;
	for(auto p = 0u; p < nspecies; p++){
		vector <unsigned int> pop_affect;                        // Lists populations affect by species
		for(auto po = 0u; po < model.pop.size(); po++){
			if(model.pop[po].sp_p == p) pop_affect.push_back(po);
		}
		
		StateSpecies ss(param_val,spline_val,model.eqn,model.param,model.param_vec,model.species[p],model.details,model.timepoint,dtimepoint,pop_affect,model.mode);
		species.push_back(ss);
	}
	
	initialise_update_ind_maps();
}


/// Simulate the state from parameter values, individual effects and initial conditions
void State::simulate(const vector <double> &param_value, const vector <InitCondValue> &initc_val) 
{
	param_val = param_value;
	
	spline_init();
	
	for(auto p = 0u; p < nspecies; p++){
		species[p].simulate_init(initc_val[p]);
	}
	
	auto popnum = model.calculate_popnum(species);  
	
	for(auto p = 0u; p < nspecies; p++){
		auto &ssp = species[p];
		if(ssp.type == INDIVIDUAL){
			ssp.activate_initial_state(popnum);
		}
	}
	
	for(auto p = 0u; p < nspecies; p++){
		auto &ssp = species[p];
		ssp.check(popnum);
	}
	
	for(auto ti = 0u; ti < T; ti++){	
		if(model.mode == SIM) print_cpop(ti);
	
		auto &pop = popnum_t[ti];
		
		pop = model.calculate_popnum(species);         // Calculates the population numbers
	
		for(auto p = 0u; p < nspecies; p++){
			auto &sp = species[p];
			switch(sp.type){
			case INDIVIDUAL: sp.update_individual_based(ti,pop); break;
			case POPULATION: sp.update_population_based(ti,model.details.stochastic,pop);
			}
			sp.check(pop);
		}
	}
	
	if(false){
		cout << "PRINT\n";	
		for(auto p = 0u; p < nspecies; p++){
			auto &ssp = species[p];
			for(auto i = 0u; i < ssp.individual.size(); i++){
				ssp.print_ind(i);
			}	
		}	
	}
	
	likelihood_from_scratch();
}


/// Calcualtes derived quantities
vector <DeriveOutput> State::derive_calculate() const
{
	vector <DeriveOutput> output;
	
	for(const auto &der : model.derive){
		DeriveOutput op;
		for(auto i = 0u; i < der.eq.size(); i++){
			auto &eqn = model.eqn[der.eq[i].eq_ref];
			
			string val_str;
			if(der.time_dep == false){
				auto val = eqn.calculate_param_only(param_val);
				val_str = tstr(val);
			}
			else{
				for(auto ti = 0u; ti < T; ti++){	
					if(ti != 0) val_str += "|";
					auto val = eqn.calculate(ti,popnum_t[ti],param_val,spline_val);
					val_str += tstr(val);
				}
			}
			op.value_str.push_back(val_str);
		}
		
		output.push_back(op);
	}
	
	return output;
}


/// Works out values for the spline based on the parameter vector
void State::spline_init()
{
	for(auto j = 0u; j < model.spline.size(); j++){
		const auto &spline = model.spline[j];
		const auto &div = spline.div;
		auto &val = spline_val[j].val;
		
		if(spline.constant == false){
			for(auto ti = 0u; ti < T; ti++){
				auto f = div[ti].f;
				val[ti] = param_val[div[ti].th1]*f + param_val[div[ti].th2]*(1-f);
			}
		}
		else val = spline.const_val;
	}
}


/// Recalculates the values for a portion of a spline 
vector <double> State::recalculate_spline(unsigned int j, const vector <unsigned int> &list)
{
	const auto &spline = model.spline[j];
	const auto &div = spline.div;
	auto &val = spline_val[j].val;
		
	if(spline.constant == true) emsg("SHould not recalculate");
	
	vector <double> store;
	for(auto ti : list){
		auto f = div[ti].f;
		
		store.push_back(val[ti]);
		val[ti] = param_val[div[ti].th1]*f + param_val[div[ti].th2]*(1-f);
	}
	
	return store;
}


/// Prints cpop to the console
void State::print_cpop(unsigned int ti) const 
{
	if(com_op == true) return;
	
	cout << model.timepoint[ti] << ": ";
	
	for(auto p = 0u; p < nspecies; p++){
		auto &ssp = species[p];
		cout << model.species[p].name << " - ";
		for(auto i = 0u; i < ssp.cpop.size(); i++) cout << ssp.cpop[i] << " ";
		cout << "   ";
	}
	cout << "\n";
}


/// Calculates the likelihood for all parts of the model
void State::calculate_likelihood()
{
	prior_prob = model.prior_prob(param_val);
	
	dist_prob = model.dist_prob(param_val);
	
	spline_prior = model.spline_prior(param_val);
	
	for(auto p = 0u; p < nspecies; p++){
		auto &ssp = species[p];

		auto temp = 0.0;
		ssp.likelihood_observation(temp);		
		
		switch(ssp.type){
		case INDIVIDUAL:
			{
				for(auto g = 0u; g < ssp.ind_eff_group_sampler.size(); g++){
					double temp = 0.0;
					ssp.likelihood_indeff_group(g,temp);
				}
			
				ssp.Li_markov.resize(ssp.N);
				for(auto e = 0u; e < ssp.N; e++){
					auto &me_vari = ssp.markov_eqn_vari[e];
					ssp.Li_markov[e].resize(me_vari.div.size());
					auto temp = 0.0;
					ssp.likelihood_markov(e,empty,ALL_DIV,temp);
				}
				
				ssp.nm_trans_ind.resize(model.species[p].nm_trans.size());
				for(auto i = 0u; i < ssp.nm_trans_ind.size(); i++){
					double temp;
					ssp.likelihood_nm_trans(i,popnum_t,temp);
				}
			}
			break;
			
		case POPULATION:
			{	
				if(model.details.stochastic == true) ssp.likelihood_pop(popnum_t);
			}
			break;
		}
	}
	
	calculate_like();
}


/// Calculates the overall likelihood 
void State::calculate_like()
{
	like.prior = sum(prior_prob);
	like.dist = sum(dist_prob);
	like.spline_prior = sum(spline_prior);
	
	like.markov = 0;
	like.nm_trans = 0;
	like.ie = 0;
	like.obs = 0;
	for(auto p = 0u; p < nspecies; p++){
		auto &ssp = species[p];
		
		like.obs += ssp.Li_obs;
		
		switch(ssp.type){
		case INDIVIDUAL:
			like.markov += sum(ssp.Li_markov);
			like.nm_trans += sum(ssp.Li_nm_trans);
			like.ie += sum(ssp.Li_ie);
			break;
			
		case POPULATION:
			like.markov += sum(ssp.Li_markov_pop);
			break;
		}
	}
}


/// Adds a change in likelihood to the state
void State::add_like(Like like_ch)
{
	like.obs += like_ch.obs;
	like.prior += like_ch.prior;
	like.dist += like_ch.dist;
	like.ie += like_ch.ie;
	like.spline_prior += like_ch.spline_prior;
	like.markov += like_ch.markov;
	like.nm_trans += like_ch.nm_trans;
}


/// Updates that state based on a change in parameters
Like State::update_param(const vector <AffectLike> &affect_like)
{	
	Like like_ch;
	like_ch.prior = 0.0; like_ch.dist = 0.0; 
	like_ch.spline_prior = 0.0;
	like_ch.ie = 0.0;
	like_ch.markov = 0.0; like_ch.nm_trans = 0.0;
	
	change.clear();
	
	auto K = affect_like.size();
	for(auto i = 0u; i < K; i++){
		const auto &alike = affect_like[i];
		
		switch(alike.type){
		case LIKE_OBS_AFFECT:
			{		
				change_add(species[alike.num].likelihood_observation(like_ch.obs));
			}
			break;
			
		case MARKOV_POP_AFFECT:
			{
				change_add(species[alike.num].likelihood_pop_change(alike.num2,alike.list,popnum_t,like_ch.markov));
			}
			break;
			
		case LIKE_IE_AFFECT:
			{
				change_add(species[alike.num].likelihood_indeff_group(alike.num2,like_ch.ie));
			}
			break;
			
		case EXP_IE_AFFECT:
			{
				change_add(species[alike.num].recalculate_exp_ie(alike.num2));
			}
			break;
			
		case OMEGA_AFFECT:
			{
				change_add(species[alike.num].calculate_omega(alike.num2));
			}
			break;
			
		case NM_TRANS_AFFECT:
			{
				change_add(species[alike.num].likelihood_nm_trans(alike.num2,popnum_t,like_ch.nm_trans));
			}
			break;
			
		case POP_AFFECT:
			{	
				change_add(model.recalculate_population(popnum_t,alike.list,species));
			}
			break;
			
		case SPLINE_PRIOR_AFFECT:
			{
				change_add(model.recalculate_spline_prior(alike.num,spline_prior,param_val,like_ch.spline_prior));
			}
			break;
			
		case PRIOR_AFFECT:
			{
				change_add(model.recalculate_prior(alike.num,prior_prob,param_val,like_ch.prior));
			}
			break;
			
		case DIST_AFFECT:
			{
				change_add(model.recalculate_prior(alike.num,dist_prob,param_val,like_ch.dist));
			}
			break;
			
		case SPLINE_AFFECT:
			{
				change_add(recalculate_spline(alike.num,alike.list));
			}
			break;
			
		case DIV_VALUE_AFFECT:   // Updates div.value on Markov transition
			{
				auto p = alike.num, e = alike.num2;
				change_add(species[p].likelihood_markov_value(e,alike.list,SEL_DIV,popnum_t));
			}
			break;
			
		case INDFAC_INT_AFFECT:      // Updates div.indfac_int on Markov transition
			{
				auto p = alike.num;
				change_add(species[p].likelihood_indfac_int());
			}
			break;
			
		case EXP_FE_AFFECT:      // Updates exp_fe on individual
			{
				auto p = alike.num, f = alike.num2;
				change_add(species[p].set_exp_fe(f));
			}
			break;
			
		case MARKOV_LIKE_AFFECT: // Updates the Markov liklihood
			{
				auto p = alike.num, e = alike.num2;	
				change_add(species[p].likelihood_markov(e,alike.list,SEL_DIV,like_ch.markov));
			}
			break;
		}
	}
	
	return like_ch;
}


/// Restores state to as before proposal
void State::restore(const vector <AffectLike> &affect_like)
{	
	for(auto i = 0u; i < change.size(); i++){
		const auto &alike = affect_like[i];
		const auto &vec = change[i].vec;
		const auto &li = alike.list;
		
		switch(alike.type){
		case LIKE_OBS_AFFECT:
			{
				species[alike.num].likelihood_observation_restore(vec);
			}
			break;
			
		case MARKOV_POP_AFFECT:
			{
				species[alike.num].likelihood_pop_change_restore(alike.num2,alike.list,vec);
			}
			break;
			
		case LIKE_IE_AFFECT:
			{
				species[alike.num].Li_ie[alike.num2] = change[i].num;
			}
			break;
			
		case EXP_IE_AFFECT:
			{
				species[alike.num].recalculate_exp_ie_restore(alike.num2,vec);
			}
			break;
			
		case OMEGA_AFFECT:
			{
				species[alike.num].calculate_omega_restore(alike.num2,vec);
			}
			break;
			
		case NM_TRANS_AFFECT:
			{
				species[alike.num].Li_nm_trans[alike.num2] = change[i].num;
			}
			break;
			
		case POP_AFFECT:
			{
				model.recalculate_population_restore(popnum_t,alike.list,vec);
			}
			break;
			
		case SPLINE_PRIOR_AFFECT:
			{
				spline_prior[alike.num] = change[i].num;
			}
			break;
			
		case PRIOR_AFFECT:
			{
				prior_prob[alike.num] = change[i].num;
			}
			break;
			
		case DIST_AFFECT:
			{
				dist_prob[alike.num] = change[i].num;
			}
			break;
			
		case SPLINE_AFFECT:
			{
				auto &val = spline_val[alike.num].val;
				for(auto k = 0u; k < li.size(); k++) val[li[k]] = vec[k];
			}
			break;
			
		case DIV_VALUE_AFFECT:   // Restore div.value on Markov transition
			{
				auto p = alike.num, e = alike.num2;
				auto &div = species[p].markov_eqn_vari[e].div;
				for(auto k = 0u; k < li.size(); k++) div[li[k]].value = vec[k];
			}
			break;
			
		case INDFAC_INT_AFFECT:      // Restore div.indfac_int on Markov transition
			{
				auto p = alike.num;
				species[p].likelihood_indfac_int_restore(vec);
			}
			break;
			
		case EXP_FE_AFFECT:  
			{
				auto p = alike.num, f = alike.num2;
				
				species[p].set_exp_fe_restore(f,vec);
			}
			break;
			
		case MARKOV_LIKE_AFFECT: // Updates the Markov liklihood
			{
				auto p = alike.num, e = alike.num2;
				for(auto k = 0u; k < li.size(); k++) species[p].Li_markov[e][li[k]] = vec[k];	
			}
			break;
		}
	}
}


/// Adds a change to the state
void State::change_add(const vector <double> &vec)
{
	Change ch; ch.vec = vec; ch.num = UNSET;
	change.push_back(ch);
}


/// Adds a change to the state
void State::change_add(double num)
{
	Change ch; ch.num = num;
	change.push_back(ch);
}


/// CHecks if the parameter set is invalid (either outside of prior or Cholesky problem)
bool State::check_invalid()
{
	// TO DO
	return true;
}


/// Sets up everything needed to represent the Markov likelihood
void State::likelihood_from_scratch()
{
	// Finds the populations
	popnum_t = model.calculate_popnum_t(species);

	// Initialises markov_eqn with value calculated and inffac_int set to zero
	for(auto p = 0u; p < species.size(); p++){
		auto &sp = model.species[p];
		auto &ssp = species[p];
		
		switch(sp.type){
		case INDIVIDUAL:
			for(auto g = 0u; g < sp.ind_eff_group.size(); g++){
				ssp.calculate_omega(g);
			}
			
			ssp.likelihood_indfac_int();
			ssp.likelihood_ind_trans();
			for(auto e = 0u; e < ssp.N; e++){
				ssp.likelihood_markov_value(e,empty,ALL_DIV,popnum_t);
			}
			break;
			
		case POPULATION:
			ssp.Li_markov_pop.clear();
			ssp.Li_markov_pop.resize(sp.tra_gl.size());
			for(auto tr = 0u; tr < sp.tra_gl.size(); tr++) ssp.Li_markov_pop[tr].resize(T,0);
			break;
		}
	}

	calculate_likelihood();
}


/// Resamples individual using the observation sampler (this gets fixed event correct
void State::resample_ind()
{
	for(auto p = 0u; p < species.size(); p++){
		auto &sp = model.species[p];
		if(sp.type == INDIVIDUAL){
			auto &ssp = species[p];
	
			IndEvSampler ind_ev_samp(ssp.markov_eqn_vari,ssp.individual,model.details,sp,ssp.obs_eqn_value,model.eqn);
			
			for(auto i = 0u; i < ssp.individual.size(); i++){
				auto &ind = ssp.individual[i];
	
				//cout << ssp.Li_obs << " Li obs\n";
				//ssp.print_event("before",ind);
				//sp.print_obs_data("obs",sp.individual[i].obs);
				//auto cinit = ind.cinit;
	
				for(auto cl = 0u; cl < sp.ncla; cl++){
					if(ind_ev_samp.needed(i,cl) == true){
						//ssp.print_event("start",ind);
						//sp.print_obs_data(ind.name,sp.individual[i].obs);
						
						ind_ev_samp.generate_ind_obs_timeline();
					
						auto probif = 0.0;
						auto ev_new = ind_ev_samp.sample_events(probif);
						
						if(ind_ev_samp.illegal == false){
							auto like_ch = update_ind(p,i,ind.cinit,ev_new);
							//cout << like_ch.obs <<" dd\n";
							//ssp.print_event("mid",ind);
				
							add_like(like_ch);
					
							ssp.recalculate_cl_ev(ind); 
						}
					}
					//ssp.print_event("after",ind);
					
					//check("  resample");
				}
			}
		}
	}		
	
	check("resample");
}


/// Generates a particle from the state
Particle State::generate_particle() const
{
	Particle part;
	part.param_val = param_val;
	part.like = like;
	part.dir_out = derive_calculate();
	
	for(auto p = 0u; p < species.size(); p++){
		const auto &ssp = species[p];
		ParticleSpecies part_sp;
		part_sp.cpop_init = ssp.cpop_init;
		part_sp.trans_num = ssp.trans_num;
		part_sp.individual = ssp.individual;
	
		part.species.push_back(part_sp);
	}
		
	return part;
}


/// Generates a particle from the state
void State::set_particle(const Particle &part)
{
	param_val = part.param_val;
	like = part.like;

	for(auto p = 0u; p < species.size(); p++){
		auto &ssp = species[p];
		const auto &part_sp = part.species[p];
		ssp.cpop_init = part_sp.cpop_init;
		ssp.trans_num = part_sp.trans_num;
		ssp.individual = part_sp.individual;
	}
	
	likelihood_from_scratch();
}


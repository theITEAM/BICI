// Checks that inputs are correctly specified

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

#include "state.hh"
#include "utils.hh"


/// Works out where memory is being used in the state
void State::profile_memory() const
{
	cout << "STATE" << endl;
	auto f = 100000.0/(1024*total_memory());
	auto sum_min = 5.0/f;
	
	auto sum_tot = 0.0;
	
	{	
		auto sum = sizeof(State);
		sum_tot += sum;
		if(sum > sum_min) cout << "State: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = mem(empty);
		sum_tot += sum;
		if(sum > sum_min) cout << "empty: " << (unsigned int)(sum*f) << endl;
	}
		vector <unsigned int> empty;           // An empty vector
		
	{	
		auto sum = 0.0;
		sum += mem(param_val.value);
		sum += mem(param_val.precalc);
		sum += mem(param_val.value_old);
		sum += mem(param_val.precalc_old);
		sum += mem(param_val.value_ch);
		sum += mem(param_val.precalc_ch);
		sum_tot += sum;
		if(sum > sum_min) cout << "param_val: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = mem(prior_prob)+mem(prior_ieg)+mem(dist_prob)+mem(spline_prior);
		sum += mem(update_timer)+mem(restore_timer)+mem(timer)+mem(check_timer);
		sum_tot += sum;
		if(sum > sum_min) cout << "various: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = mem(popnum_t);
		sum_tot += sum;
		if(sum > sum_min) cout << "popnum_t: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = mem(popcomb_t);
		sum_tot += sum;
		if(sum > sum_min) cout << "popcomb_t: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = 0.0;
		for(auto i = 0u; i < popnum_ind.size(); i++){
			for(auto j = 0u; j < popnum_ind[i].size(); j++){
				sum += popnum_ind[i][j].size()*sizeof(PopIndRef);
			}
		}
		sum_tot += sum;
		if(sum > sum_min) cout << "popnum_ind: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = sizeof(GeneticDataValue);
		sum += genetic_value.obs_node_ref.size()*sizeof(NodeRef);
		sum += genetic_value.inf_origin.size()*sizeof(InfOrigin);
		sum += mem(genetic_value.gen_dif);
		for(const auto &in : genetic_value.inf_node){
			sum += sizeof(InfNode) + in.inf_ev.size()*sizeof(InfEvent);			
		}
		sum_tot += sum;
		if(sum > sum_min) cout << "genetic_value: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = 0.0;
		for(auto &val : change) sum += sizeof(Change) + mem(val.vec);
		sum_tot += sum;
		if(sum > sum_min) cout << "change: " << (unsigned int)(sum*f) << endl;
	}

	{
		auto sum = 0.0;
		for(const auto &val : back_pop){
			sum += sizeof(BackPop)+mem(val.vec);
		}
		sum_tot += sum;
		
		if(sum > sum_min) cout << "back_pop: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = alg_warn.size()*sizeof(AlgWarn);
		sum_tot += sum;
		if(sum > sum_min) cout << "alg_warn: " << (unsigned int)(sum*f) << endl;
	}
	
	for(const auto &ssp : species){
		cout << "FOR SPECIES " << endl;
		{
			auto sum = sizeof(StateSpecies);
			sum_tot += sum;
			if(sum > sum_min) cout << "StateSpecies: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.Li_markov)+mem(ssp.Li_ie)+mem(ssp.Li_obs_ind)
			          +mem(ssp.Li_obs_pop)+mem(ssp.Li_obs_pop_trans);
			sum_tot += sum;
			if(sum > sum_min) cout << "Li: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.trans_num);
			sum_tot += sum;
			if(sum > sum_min) cout << "trans_num: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.tnum_mean_st);
			sum_tot += sum;
			if(sum > sum_min) cout << "tnum_mean_st: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.Li_markov_pop);
			sum_tot += sum;
			if(sum > sum_min) cout << "Li_markov_pop: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.cpop_st);
			sum_tot += sum;
			if(sum > sum_min) cout << "cpop_st: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			const auto &ic = ssp.init_cond_val;
			auto sum = sizeof(InitCondValue)+mem(ic.N_focal)+mem(ic.N_focal_unobs)+mem(ic.cnum_reduce)
			          +mem(ic.frac_focal)+mem(ic.frac_comb)+mem(ic.cnum)+mem(ic.frac);
			sum_tot += sum;
			if(sum > sum_min) cout << "cpop_st: " << (unsigned int)(sum*f) << endl;
		}
	
		{
			auto sum = mem(ssp.trans_num_f);
			sum_tot += sum;
			if(sum > sum_min) cout << "trans_num_f: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.tnum_mean_st_f);
			sum_tot += sum;
			if(sum > sum_min) cout << "tnum_mean_st_f: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.cpop_st_f);
			sum_tot += sum;
			if(sum > sum_min) cout << "cpop_st_f: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &ind : ssp.individual){
				sum += sizeof(Individual)+mem(ind.name)+mem(ind.ie)+mem(ind.exp_ie)+mem(ind.X)+mem(ind.exp_fe);
				sum += ind.ev.size()*sizeof(Event);
				sum += ind.popnum_ind_ref.size()*sizeof(PopnumIndRef);
				sum += ind.incomp_ref.size()*sizeof(IncompNMTransRef);
				sum += ind.markov_eqn_ref.size()*sizeof(MarkovEqnRef);
			}				
			sum_tot += sum;
			if(sum > sum_min) cout << "individual: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.pop_data_num);
			sum_tot += sum;
			if(sum > sum_min) cout << "pop_data_num: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.pop_trans_data_num);
			sum_tot += sum;
			if(sum > sum_min) cout << "pop_trans_data_num: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.pop_data_cgl);
			sum_tot += sum;
			if(sum > sum_min) cout << "pop_data_cgl: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.pop_trans_data_tgl);
			sum_tot += sum;
			if(sum > sum_min) cout << "pop_trans_data_tgl: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.obs_eqn_value);
			sum_tot += sum;
			if(sum > sum_min) cout << "obs_eqn_value: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.unused_name);
			sum_tot += sum;
			if(sum > sum_min) cout << "unused_name: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : ssp.nm_trans_ev_ref){
				for(const auto &val2 : val){
					sum += val2.size()*sizeof(NMEventRef);
				}
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "nm_trans_ev_ref: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.nm_trans_incomp_ref);
			sum_tot += sum;
			if(sum > sum_min) cout << "nm_trans_incomp_ref: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.obs_trans_eqn_value);
			sum_tot += sum;
			if(sum > sum_min) cout << "obs_trans_eqn_value: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.obs_trans_eqn_num);
			sum_tot += sum;
			if(sum > sum_min) cout << "obs_trans_eqn_num: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.cpop);
			sum_tot += sum;
			if(sum > sum_min) cout << "cpop: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.X_order);
			sum_tot += sum;
			if(sum > sum_min) cout << "X_order: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : ssp.ind_eff_group_sampler){
				sum += mem(val.omega)+mem(val.omega_Z)+mem(val.omega_inv)+mem(val.A)
				      +mem(val.A_Z)+mem(val.A_inv)+mem(val.A_inv_diag);
				for(const auto &val2 : val.A_inv_nonzero){
					sum += val2.size()*sizeof(NZElement);
				}
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "ind_eff_group_sampler: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : ssp.markov_eqn_vari){
				sum += sizeof(MarkovEqnVariation);
				sum += val.ind_tra.size()*sizeof(IndTransRef);
				for(const auto &val2 : val.div){
					sum += sizeof(MarkovEqnDiv)+val2.ind_trans.size()*sizeof(DivIndRef);
				}
			}
			sum_tot += sum;
			if(sum > sum_min) cout << " markov_eqn_vari: " << (unsigned int)(sum*f) << endl;
		}
	
		{
			auto sum = 0.0;
			const auto &ss = ssp.source_sampler;
			sum += sizeof(SourceSampler)+ss.get_mem();
			sum_tot += sum;
			if(sum > sum_min) cout << "source_sampler: " << (unsigned int)(sum*f) << endl;
		}
	
		{
			auto sum = 0.0;
			const auto &rm = ssp.rate_mean;
			sum += sizeof(RatePosteriorMean);
			sum += mem(rm.tramean_ref)+mem(rm.compR_ref);
			
			for(const auto &val : rm.tra_rate){
				sum += sizeof(TraRate)+mem(val.value_raw)+mem(val.bp_raw)+mem(val.value)
				                      +mem(val.integral)+mem(val.bp);
				sum += val.ind_fac_rate.size()*sizeof(IndFacRate);
			}				
			
			for(const auto &val : rm.compR){
				sum += sizeof(CompR)+mem(val.list)+mem(val.integral);
			} 
			sum_tot += sum;
			if(sum > sum_min) cout << "rate_mean: " << (unsigned int)(sum*f) << endl;
		}
	
		{
			auto sum = 0.0;
			for(const auto &val : ssp.back_ind){
				sum += mem(val.identical_old)+mem(val.identical_new)
				      +mem(val.remove_old)+mem(val.add_new);
				sum += val.ev.size()*sizeof(Event);
				sum += val.incomp_ref.size()*sizeof(IncompNMTransRef);
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "back_ind: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : ssp.back){
				sum += sizeof(Back)+mem(val.vec);
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "back: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(ssp.ind_sim_c);
			sum_tot += sum;
			if(sum > sum_min) cout << "ind_sim_c: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : ssp.trig_div){
				sum += val.ev.size()*sizeof(SimTrigEvent);
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "trig_div: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(timer);
			sum_tot += sum;
			if(sum > sum_min) cout << "timer: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = ssp.alg_warn.size()*sizeof(AlgWarn);
			sum_tot += sum;
			if(sum > sum_min) cout << "alg_warn: " << (unsigned int)(sum*f) << endl;
		}
	}
}


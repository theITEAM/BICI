// Implements an MCMC algorithm

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm> 
 
using namespace std;

#include "mcmc.hh"
#include "state.hh"
#include "utils.hh"
#include "matrix.hh"

MCMC::MCMC(const Model &model, Output &output) : model(model), output(output), state(model), update(model,burn_info,output)
{	
	const auto &de = model.details;
	
	if(de.anneal_type == ANNEAL_SCAN){
		nburnin = LARGE;
		nsample = LARGE;
		burn_info.prop_join_step = 20;
	}
	else{
		auto f = de.burnin_frac/100.0;
		if(true){
			nburnin = f*de.sample;
			nsample = de.sample;
		}
		else{
			nburnin = (f/(1-f))*de.sample;
			nsample = de.sample + nburnin;
		}
		burn_info.prop_join_step = nburnin/4;
	}
	
	output_param = de.output_param;
	output_state = de.output_state;

	nparam_op = 0;
	nstate_op = 0;
	
	burn_info.on = true;
	burn_info.range = LARGE;
}
 
 
/// Runs MCMC
void MCMC::run()
{
	print("Run...");
	
	update.init();
	
	if(false) output.print_param(state.param_val);
	
	Particle part;
	
	auto like_obs_max = -LARGE;
	
	print("Finding initial state");

	auto loop_max = 10u;
	for(auto loop = 0u; loop < loop_max; loop++){ 
		print(to_string(loop)+" Initial state");
		
		auto param_val = model.param_sample();
		auto initc_val = model.initc_sample(param_val);
		
		if(false){
			output.print_initc(initc_val);
			output.print_param(param_val);
		}
		
		print("simulate");
		
		state.simulate(param_val,initc_val);
	
		print("simulate done");
	
		state.check("First check");
	
		print("Checked");
		
		state.resample_ind();    // Resamples individual such that fixed events become correct

		state.check("Resample check");
	
		print("resampled");
		
		burn_info.add_L(state.like); 
		
		auto like_obs = state.like.obs + state.like.genetic_obs;
		if(like_obs > like_obs_max){
			like_obs_max = like_obs;
			part = state.generate_particle();
		}
	
		state.check("after resample");
	}

	state.set_particle(part);
	
	state.check("before");
	
	print("Starting sampling");

	long time_start = clock();

	auto op_step = get_op_step(nsample);
	
	for(auto s = 0u; s < nsample; s++){
		if(s%100 == 0 && com_op) progress(s,nsample);		
		else{
			if(s%op_step == 0 && op()) cout << s << " " << nsample << " sample" << endl;
		}
	
		burn_info.add_L(state.like); 
		burn_info.setup(s,nburnin,nsample,model.details);
		
		if(param_op(s)) output.param_sample(s,state);
		
		if(state_op(s)) output.state_sample(s,state);
		
		if(burn_info.on){
			update.cor_matrix.add_sample(state.param_val);
		
			if(s%10 == 0){
				state.update_individual_sampler();
				for(auto &pro : update.proposal) pro.update_sampler(update.cor_matrix);
			}
			
			//if(s != 0 && s%burn_info.prop_join_step == 0) update.check_join_proposal();
	
			if(adapt_prop_prob && s%UPDATE_PROP_PROB == UPDATE_PROP_PROB-1) update.set_proposal_prob();
		}
	
		update.run(s,state);
	}
	
	output.set_output_burnin(double(100.0*nburnin)/nsample);
	
	update.diagnostics(clock()-time_start,state);
	
	if(false){
		state.output_gen_dif();
		state.output_waifw("waifw2.txt");
	}
}


/// Sets up quantities during the burn-in phase
void BurnInfo::setup(unsigned int s, unsigned int &nburnin,  unsigned int &nsample, const Details &details)
{
	dprob_suppress = false;
	
	if(s < nburnin){
		on = true; 	
		if(s < 100 && s < nburnin/8) fac = 20;
		else{ if(s < 200 && s < nburnin/2) fac = 10; else fac = 1;}
		
		if(s < nburnin/8) dprob_suppress = true;
	}
	else on = false;
	
	auto type = details.anneal_type;
	
	// Estimates the current standard deviation in the log likelihood
	auto dphi = 0.0;
	if(s == 0 || type == ANNEAL_SCAN){
		auto n = L_samp.size()-start;
		auto var = L_samp_sum2/n - (L_samp_sum/n)*(L_samp_sum/n);
		if(var < TINY) var = TINY;
		dphi = 1.0/sqrt(var);
		if(n == 0) emsg("Should not be zero");
	}
	
	switch(type){
	case ANNEAL_NONE:
		phi = 1.0;
		range = 0;
		break;
		
	case ANNEAL_SCAN:
		if(false){
			cout << s << " " << phi << " " << details.anneal_rate << "  rate" << endl;
		}
		
		if(s == 0){
			phi = 0.0;
			range = 50;
		}
		else{
			if(on == true && phi < PHI_LIM){
				phi += details.anneal_rate*dphi;
				if(phi > PHI_LIM){
					phi = PHI_LIM;
					nburnin = (unsigned int) (1.5*s);
					nsample = details.sample + nburnin;
				}
			}
		}
		break;
		
	case ANNEAL_POWERAUTO:
		{
			if(s == 0){
				auto val = 0.01*dphi;
			
				auto fr_cut = 0.05;
				auto power_thresh = -1.0/log(fr_cut);
				auto peak = power_thresh*pow(fr_cut,power_thresh-1);
				
				if(val > peak) power = power_thresh;
				else{
					power = 2*power_thresh;
				
					// Uses Newtons method to solve for power 
					auto loop = 0u, loopmax = 100u;
					for(loop = 0; loop < loopmax; loop++){
						auto power_st = power;
						auto f = log(power/val)+(power-1)*log(fr_cut);
						auto df_dpower = 1.0/power +log(fr_cut);
						power -= f/df_dpower;
						if(power < 1.1*power_thresh) power = 1.1*power_thresh;
						if((power-power_st)*(power-power_st) < TINY) break;
					}
					if(loop == loopmax) emsg("Convergence problem");
				}
			
				range = 0;
				phi = 0;
			}
			else{
				if(on == true){
					auto fr = double(s)/nburnin;
					phi = pow(fr,power);
				}
				else phi = 1;
			}
		}
		break;
		
	case ANNEAL_LOGAUTO:
		{
			if(s == 0){
				auto val = 0.01*dphi;
	
				if(val > exp(-1)) alpha = 1;
				else{
					alpha = 2;
				
					// Uses Newtons method to solve for alpha 
					auto loop = 0u, loopmax = 100u;
					for(loop = 0; loop < loopmax; loop++){
						auto alpha_st = alpha;
						auto f = log(alpha/val)-alpha;
						auto df_dalpha = 1.0/alpha - 1;
						alpha -= f/df_dalpha;
						if(alpha < 1.1) alpha = 1.1;
						if((alpha-alpha_st)*(alpha-alpha_st) < TINY) break;
					}
					if(loop == loopmax) emsg("Convergence problem");
				}
			
				range = 0;
				epsilon = 1/(exp(alpha)-1);
				phi = 0;
			}
			else{
				if(on == true){
					auto fr = double(s)/nburnin;
					phi = (1+epsilon)*exp(alpha*(fr-1)) - epsilon;
				}
				else phi = 1;
			}
		}
		break;
		
	case ANNEAL_POWER:
		{
			if(s == 0) power = details.anneal_power;
			
			if(on == true){
				auto fr = double(s)/nburnin;
				phi = pow(fr,power);
			}
			else phi = 1;
		}
		break;
	}
	
	if(false){ // This code was used when testing different annealing schemes
		double origin, log_ratio_av;
		if(s == 0){
			phi = 0.0;
			range = nburnin/20; if(range < 2) range = 2;
			power = 8;
			origin = 0;
			log_ratio_av = UNSET;
		}
		else{
			if(on == true){
				auto fr = double(s-origin)/(nburnin-origin);
			
				// The power is updated to make try to make changes based on standard deviation in likelihood 
				
				// This is the change based on SD
				auto n = L_samp.size()-start;
				auto var = L_samp_sum2/n - (L_samp_sum/n)*(L_samp_sum/n);
				if(var < TINY) var = TINY;
				auto dphi = 1.0/sqrt(var);
				
				// This is the change based on power law model
				auto dphi2 = power*pow(fr,power-1);
				
				if(true){  // This scans up the inverse temperature
					phi += 0.003*dphi;
					if(phi > 1) phi = 1;
				}
				else{
					phi = pow(fr,power);
				
					// Calculate ratio between power law update and update based on sd of L
					auto log_ratio = log(dphi/dphi2);
					
					if(log_ratio_av == UNSET){
						log_ratio_av = log_ratio;
					}
					else{
						auto f = 1-1.0/range;
						log_ratio_av = f*log_ratio_av + (1-f)*log_ratio;
					}
					
					//auto dpower = log_ratio-log_ratio_av;
					auto dpower = 0.0;
					auto dif = -(log_ratio-log_ratio_av)/50;
					if(dif > 1) dif = 1;
					if(dif < -1) dif = -1;
					
					dpower += dif;
						
					if(power < 1) power = 1; 
					if(power > 100) power = 100;
					
					auto fr_new = pow(phi,1.0/power);
					
					origin = (s-fr_new*nburnin)/(1-fr_new);
				}
			}
		}
	}
	
	init_cond = phi;
	init_cond_prior = 1;
	obs = phi;
	prior = 1;
	spline_prior = 1;
	dist = 1;
	markov = phi;
	nm_trans = phi;
	genetic_process = phi;
	genetic_obs = phi;
	ie = 1;
	
	if(s < nburnin){
		auto fr = 1.2*double(s)/nburnin;
		if(fr > 1) fr = 1;
		genetic_obs = pow(fr,4);
	}
	else genetic_obs = 1;
}


/// Adds a likelihood sample to burnin 
void BurnInfo::add_L(const Like &like)
{
	if(!on) return;
	
	auto L = like.obs + like.markov + like.nm_trans;
		
	L_samp.push_back(L);
	L_samp_sum += L;
	L_samp_sum2 += L*L;
	while(L_samp.size()-start > range){
		L_samp_sum -= L_samp[start];
		L_samp_sum2 -= L_samp[start]*L_samp[start];
		start++;
	}
}



/// Determines if a parameter output is made
bool MCMC::param_op(unsigned int s)
{
	if(nparam_op <= double(s*output_param)/nsample){
		nparam_op++;
		return true;
	}
	return false;
}


/// Determines if a state output is made
bool MCMC::state_op(unsigned int s)
{
	if(!burn_info.on){
		if(nstate_op <= double((s-nburnin)*output_state)/(nsample-nburnin)){
			nstate_op++;
			return true;
		}
	}
	return false;
}



// Implements an MCMC chain

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm> 
 
using namespace std;

#include "chain.hh"
#include "state.hh"
#include "utils.hh"
#include "matrix.hh"
#include "synchronise.hh"

Chain::Chain(unsigned int nburnin_, unsigned int nsample_, const Model &model, Output &output) : state(model), cor_matrix(model), model(model), output(output)
{
	nburnin = nburnin_;
	nsample = nsample_;
	
	const auto &de = model.details;
	
	if(de.anneal_type == ANNEAL_SCAN){
		burn_info.prop_join_step = 20;
	}
	else{
		auto val = nburnin/4; if(val == 0) val = 1;
		burn_info.prop_join_step = val;
	}
	
	burn_info.on = true;
	burn_info.range = LARGE;
};


/// Initialises chain
void Chain::init(unsigned int ch, unsigned int ch_max)
{
	print_diag("update_init");
		
	update_init();
	
	print_diag("start init_chain");
	
	vector < vector <IndRef> > not_possible;
	
	if(model.mode == EXT){    // Loads initial state
		if(model.sample.size() == 0) emsg("Must be a sample");
	
		auto smax = 0u, imax = 0u;
		for(auto i = 0u; i < model.sample.size(); i++){
			const auto &sa = model.sample[i];
			
			if(sa.ch == ch && sa.num > smax){ smax = sa.num; imax = i;}
		}
		
		const auto &samp = model.sample[imax];
		
		auto param_val = model.post_param(samp);
		state.load_samp(param_val,samp);
	}
	else{                     // Simulates initial state
		auto fit_max = -LARGE;
	
		print_diag("Finding initial state...");

		// Determines which metric to use (likelihood of observation or resuidual method 
		auto use_li_obs = true;
		for(const auto &sp : model.species){ 
			if(sp.type == INDIVIDUAL) use_li_obs = false;
		}
		
		Particle part;

		auto loop_max = model.details.chain_nsiminit;

		auto param_latin = model.param_latin_hypercube_prior_sample(loop_max);
		
		for(auto loop = 0u; loop < loop_max; loop++){ 
			//cout << endl << loop << " loop" << endl;
			print_diag(to_string(loop)+" Initial state");
			
			percentage(loop+ch*loop_max,loop_max*ch_max);
			
			auto param_val = model.param_latin_hypercube_sample(param_latin[loop]);
			
			auto initc_val = model.initc_sample(param_val);
			
			if(false){
				output.print_initc(initc_val);
				model.print_param(param_val);
			}
		
			print_diag("simulate");
			
			//cout << "switch off" << endl; model.precalc_eqn.print_calc(); 
	
			state.simulate(param_val,initc_val);

			print_diag("simulate done");
		
			//state.check("First check");
		
			print_diag("Checked");
			
			print_diag("resample");
			 	
			state.resample_ind(false);    // Resamples individual such that fixed events become correct

			auto np = state.force_consistent_solution();
		
			if(np.size() > 0){  // Some individuals are not possible
				not_possible.push_back(np);
			}
			else{
				print_diag("resampled");

				//state.check("Resample check");
			
				print_diag("resampled");
				
				auto L = like_total_obs();
				burn_info.add_L(L); 

				double fit;
				if(use_li_obs){  // This fit metric uses the log-likelihood of observations (good for population-based)
					fit = L;
				}
				else{            // This fit metric uses the fit between the transitions and the model (good for individual-based data)
					fit = -(state.frac_outside_CI()-0.1);
					if(fit > 0) fit *= -1;
				}
			
				if(fit > fit_max){
					fit_max = fit;
					part = state.generate_particle(UNSET,UNSET,true,false);
				}
			
				if(false){
					for(auto th : model.param_vec_latin){
						cout << model.param_vec_name(th) << " " << param_val.value[th] << ",";
					}
					cout << fit << " fit   best:" << fit_max << endl;
				}
				
				print_diag("resampled after");
				//state.check("after resample");
			}
		}
		
		if(fit_max == -LARGE) state.no_valid_state(not_possible);

		state.set_particle(part);
	}
	
	for(const auto &ssp: state.species) ssp.check_obs_inconsitent();
	
	state.check("init");
	
	print_diag("done");
}


/// Iterates chain
void Chain::burn_update(unsigned int s)
{
	burn_info.add_L(like_total_obs()); 
	burn_info.setup(s,nburnin,nsample,model.details);
	
	if(burn_info.on){
		state.dif_thresh = DIF_THRESH_BURNIN;
		
		cor_matrix.add_sample(state.get_param_val_prop(),cor_matrix.n*FRAC_COR_UNUSED);
	
		if(s%10 == 0){
			state.update_individual_sampler();
			for(auto &pro : proposal) pro.update_sampler(cor_matrix);
		}
		
		if(s != 0 && s%burn_info.prop_join_step == 0) check_join_proposal();
		
		//if(adapt_prop_prob && s%UPDATE_PROP_PROB == UPDATE_PROP_PROB-1) set_proposal_prob();
	}
	else{
		state.dif_thresh = DIF_THRESH;
	}
}


/// Gets the value for nburnin
unsigned int Chain::get_nburnin()
{
		return nburnin;
}


/// Iterates chain
void Chain::pas_burn_update(unsigned int s, unsigned int gen_update, double phi)
{
	burn_info.add_L(like_total_obs()); 
	
	burn_info.pas_setup(gen_update,phi);
	
	cor_matrix.add_sample(state.get_param_val_prop(),gen_update);

	if(s%10 == 0){
		state.update_individual_sampler();
		for(auto &pro : proposal) pro.update_sampler(cor_matrix);
	}
	
	if(s != 0 && s%burn_info.prop_join_step == 0) check_join_proposal();

	//if(adapt_prop_prob && s%UPDATE_PROP_PROB == UPDATE_PROP_PROB-1) set_proposal_prob();
}


/// Iterates chain
void Chain::pas_burn_update_run(unsigned int s)
{
	burn_info.add_L(like_total_obs()); 
	
	burn_info.pas_setup_run(s,nburnin);
	
	if(burn_info.on){
		cor_matrix.add_sample(state.get_param_val_prop(),cor_matrix.n*FRAC_COR_UNUSED);
	
		if(s%10 == 0){
			state.update_individual_sampler();
			for(auto &pro : proposal) pro.update_sampler(cor_matrix);
		}
		
		if(s != 0 && s%burn_info.prop_join_step == 0) check_join_proposal();

		//if(adapt_prop_prob && s%UPDATE_PROP_PROB == UPDATE_PROP_PROB-1) set_proposal_prob();
	}
	
	state.dif_thresh = DIF_THRESH;
}


/// Updates the size of the proposal (if in burnin phase)
void BurnInfo::update_si(double &si, PropResult res, unsigned int ntr) const 
{
	if(on){
		auto fac = 0.01;
		switch(res){
		case ACCEPT: fac = 0.01; break;
		case ACCEPT50: fac = 0.005; break;
		case REJECT: fac = -0.005; break;
		case ACCEPT_SMALL: fac = 0.001; break;
		case ACCEPT50_SMALL: fac = 0.0005; break;
		case REJECT_SMALL: fac = -0.0005; break;
		}
		
		switch(res){
		case ACCEPT: case ACCEPT50: case REJECT:
			if(ntr < 20) fac *= 50;
			else{
				if(ntr < 200) fac *= 10;
			}
			break;
			
		case ACCEPT_SMALL: case ACCEPT50_SMALL: case REJECT_SMALL:
			if(ntr < 20) fac *= 500;
			else{
				if(ntr < 100) fac *= 100;
				else{
					if(ntr < 200) fac *= 10;
				}
			}
			break;
		}
	
		si *= 1+fac;
		if(si < 0.001) si = 0.001;
		if(si > 1000) si = 1000;
	}
}


/// Sets up quantities during the burn-in phase
void BurnInfo::setup(unsigned int s, unsigned int &nburnin,  unsigned int &nsample, const Details &details)
{
	dprob_suppress = false;
	
	if(s < nburnin){
		on = true; 	
		//if(s < 100 && s < nburnin/8) fac = 20;
		//else{ if(s < 200 && s < nburnin/2) fac = 10; else fac = 1;}
		
		//if(s < nburnin/8) dprob_suppress = true;
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
					//auto dpower = 0.0;
					auto dif = -(log_ratio-log_ratio_av)/50;
					if(dif > 1) dif = 1;
					if(dif < -1) dif = -1;
					
					//dpower += dif;
						
					if(power < 1) power = 1; 
					if(power > 100) power = 100;
					
					auto fr_new = pow(phi,1.0/power);
					
					origin = (s-fr_new*nburnin)/(1-fr_new);
				}
			}
		}
	}
	
	set_phi();
}


/// Sets up quantities during the burn-in phase
void BurnInfo::pas_setup(unsigned int gen_update, double phi_)
{
	phi = phi_;
	dprob_suppress = false;
	
	on = true;
	//if(g == 0 && s < gen_update/2) fac = 20;
	//else fac = 1;
	
	set_phi();
	
	range = gen_update/2;
}


/// Sets up quantities during the burn-in phase
void BurnInfo::pas_setup_run(unsigned int s, unsigned int &nburnin)
{
	phi = 1;
	dprob_suppress = false;
	
	//fac = 1;
	if(s < nburnin) on = true; 
	else on = false;
	
	set_phi();
}


/// Sets inverse temperatures based on global value phi
void BurnInfo::set_phi()
{
	init_cond = phi;
	init_cond_prior = 1;
	obs = phi;
	prior = 1;
	prior_bounded = phi;
	spline_prior = 1;
	dist = 1;
	markov = 1;
	nm_trans = 1;
	genetic_process = 1;
	genetic_obs = phi;
	ie = 1;
	
	if(false){
		init_cond = 0;
		init_cond_prior = 0;
		obs = 0;
		prior = 1;
		prior_bounded = 1;
		spline_prior = 1;
		dist = 0;
		markov = 0;
		nm_trans = 0;
		genetic_process = 0;
		genetic_obs = 0;
		ie = 1;
	}
}


/// Adds a likelihood sample to burnin 
void BurnInfo::add_L(double L)
{
	if(!on) return;
	
	L_samp.push_back(L);
	L_samp_sum += L;
	L_samp_sum2 += L*L;
	while(L_samp.size()-start > range){
		L_samp_sum -= L_samp[start];
		L_samp_sum2 -= L_samp[start]*L_samp[start];
		start++;
	}
}


/// Calculates the total observation likelihood
double Chain::like_total_obs() const
{
	const auto &like = state.like;
	
	return like.init_cond + like.prior_bounded + like.obs + like.genetic_obs;		 
	//return like.init_cond + like.prior_bounded + like.obs + like.markov + like.nm_trans + 
	//     like.genetic_process + like.genetic_obs;
}	


/// Updates all the proposals
void Chain::update(unsigned int s)
{
	state.sample = s;

	auto pl = false;
	
	for(auto &pro : proposal){ 
		//state.check_markov_value_dif();
		//if(core() == 4) cout << pro.name  << endl;
	
		if(pro.on){
			//if(op()) state.check_markov_value_dif();
				
			//if(op()) cout << s << " " << core() << " " << add_escape_char(pro.name) << " " <<  state.sample << " proposal" << endl;
			
			auto ch_on = false;
		//cout << s << " " << core() << " " << pro.name << " " <<  state.sample << " proposal" << endl; ch_on = true;
			
			if(true){
				if(pl) cout << s << " " << core() << " " << pro.name << " " <<  state.sample << " proposal" << endl;
		
				pro.update(state);
			}
			else{
				switch(pro.type){
				case IND_EVENT_TIME_PROP:
				//case DET_IC_POP_PROP:
				//case INIT_COND_FRAC_PROP:
					if(pl) cout << s << " " << core() << " " << pro.name << " " <<  state.sample << " proposal" << endl;
			
					pro.update(state);
					break;
			
				default:
					break;
				}
			}

			//if(core() == 4) 
			if(ch_on) state.check(" After prop check"); 
		 
			if(pl) state.check(" After prop check");
			if(pl) state.check_pop_t("hhh");
		
			//state.check_precalc_eqn("up");
			//state.check_neg_rate(pro.name);
			//state.check_markov_div_value(0," after update");
		}
	}
	
	// This resamples to remove any inconsistent indiviudals 
	if(burn_info.on){
		auto fl = false;
		for(const auto &ssp : state.species){
			if(ssp.inconsistent_ind) fl = true;
		}
		if(fl) state.resample_ind(true);   
	}
	
	//if(testing == true && s%10 == 0) state.check("check");
	if(testing == true && s%CHECK_THIN == CHECK_THIN-1) state.check("end");
}


/// Initialises proposal distributions
void Chain::update_init()
{	
	print_diag("update_init 1");

	cor_matrix.init();

	for(auto i = 0u; i < model.nparam_vec; i++){     // Univariate distributions
		const auto &pv = model.param_vec[i];
		if(pv.prop_pos){
			vector <unsigned int> vec;
			vec.push_back(i);
			
			add_parameter_prop(vec,false);
			
			if(pv.strictly_positive) add_parameter_prop(vec,true);
		}
	}
	
	print_diag("update_init 2");
	
	// Proposals on individual effects
	for(auto p = 0u; p < model.nspecies; p++){
		auto &sp = model.species[p];
		if(sp.type == INDIVIDUAL){
			for(auto i = 0u; i < sp.ind_effect.size(); i++){
				vector <unsigned int> vec;
				vec.push_back(p); vec.push_back(i);
				Proposal pp(IE_PROP,vec,model,output,0.3,burn_info);
				proposal.push_back(pp);
			}
		}
	}
	
	// Does proposals which simultaneously change variances and ies
	for(auto p = 0u; p < model.nspecies; p++){
		auto &sp = model.species[p];
		
		for(auto i = 0u; i < sp.ind_eff_group.size(); i++){
			auto &ieg = sp.ind_eff_group[i];
			const auto &par = model.param[ieg.th];
			if(par.variety != CONST_PARAM && par.variety != DYNAMIC_PARAM){
				for(auto j = 0u; j < ieg.list.size(); j++){	
					auto ie = ieg.list[j].index;
					vector <unsigned int> vec;
					vec.push_back(ieg.omega_pv[j][j]);
					vec.push_back(p);
					vec.push_back(ie);
					
					Proposal pp(IE_VAR_PROP,vec,model,output,1,burn_info);
					proposal.push_back(pp);

					// Adds proposals for joint ie, variance and cv
					auto list = sp.get_cv_list(ie,model.param,model.eqn);
					if(list.size() > 0){
						for(auto th : list) vec.push_back(th);
										
						Proposal pp(IE_VAR_CV_PROP,vec,model,output,1,burn_info);
						proposal.push_back(pp);
					}
				}
			
				if(ieg.list.size() > 0){
					for(auto j = 0u; j < ieg.list.size(); j++){
						for(auto k = j+1; k < ieg.list.size(); k++){
							vector <unsigned int> vec;
							vec.push_back(ieg.omega_pv[j][k]);
							vec.push_back(p);
							vec.push_back(i);
							vec.push_back(j);
							vec.push_back(k);		
							Proposal pp(IE_COVAR_PROP,vec,model,output,1,burn_info);
							proposal.push_back(pp);
						}
					}
				}
			}
		}
	}
	
	print_diag("update_init 3");
	
	// Adds MBPs which don't change parameters
	for(auto p = 0u; p < model.nspecies; p++){
		if(model.species[p].type == POPULATION){
			vector <unsigned int> vec;
			vec.push_back(p);
			Proposal pp(MBPII_PROP,vec,model,output,PROP_MBPII_W,burn_info);
			model.add_like_obs_affect(p,pp.affect_like);
			proposal.push_back(pp);
		}
	}
	
	print_diag("update_init 4");
	
	// Adds MBPs which change initial conditions
	for(auto loop = 0u; loop < 1; loop++){
		for(auto p = 0u; p < model.nspecies; p++){
			const auto &sp = model.species[p];
			const auto &ic = sp.init_cond;

			if(ic.type == INIT_POP_DIST){
				switch(sp.type){
				case POPULATION:
					{
						if(sp.ncla > 1){
							vector <unsigned int> vec;
							vec.push_back(p);
							
							Proposal pp(MBP_IC_RESAMP_PROP,vec,model,output,PROP_MBPII_W,burn_info);		
							model.add_like_obs_affect(p,pp.affect_like);
							proposal.push_back(pp);
						}					
							
						auto foc_cl = ic.focal_cl;
						if(foc_cl == UNSET){
							vector <unsigned int> vec;
							vec.push_back(p);
							
							Proposal pp(MBP_IC_POPTOTAL_PROP,vec,model,output,PROP_MBPIC_W,burn_info);
							AffectLike al; al.type = PRIOR_INIT_COND_AFFECT; al.num = UNSET; al.num2 = UNSET;
							param_vec_add_affect(pp.affect_like,al);		
							model.add_like_obs_affect(p,pp.affect_like);
							proposal.push_back(pp);
						}
						else{
							const auto &claa = sp.cla[foc_cl];
							for(auto c = 0u; c < claa.comp.size(); c++){
								if(!claa.comp[c].erlang_hidden){
									vector <unsigned int> vec;
									vec.push_back(p); vec.push_back(foc_cl); vec.push_back(c);
									
									Proposal pp(MBP_IC_POP_PROP,vec,model,output,PROP_MBPIC_W,burn_info);
									AffectLike al; al.type = PRIOR_INIT_COND_AFFECT; al.num = UNSET; al.num2 = UNSET;
									param_vec_add_affect(pp.affect_like,al);		
									model.add_like_obs_affect(p,pp.affect_like);
									proposal.push_back(pp);
								}
							}
						}
					}
					break;
					
				case DETERMINISTIC:
					{
						if(sp.ncla > 1){
							vector <unsigned int> vec;
							vec.push_back(p);
							
							Proposal pp(DET_IC_RESAMP_PROP,vec,model,output,PROP_MBPII_W,burn_info);		
							proposal.push_back(pp);
						}					
							
						auto foc_cl = ic.focal_cl;
						if(foc_cl == UNSET){
							vector <unsigned int> vec;
							vec.push_back(p);
							
							Proposal pp(DET_IC_POPTOTAL_PROP,vec,model,output,PROP_MBPIC_W,burn_info);
							proposal.push_back(pp);
						}
						else{
							const auto &claa = sp.cla[foc_cl];
							for(auto c = 0u; c < claa.comp.size(); c++){
								if(!claa.comp[c].erlang_hidden){
									vector <unsigned int> vec;
									vec.push_back(p); vec.push_back(foc_cl); vec.push_back(c);
									
									Proposal pp(DET_IC_POP_PROP,vec,model,output,PROP_MBPIC_W,burn_info);
									proposal.push_back(pp);
								}
							}
						}
					}
					break;
					
				case INDIVIDUAL:
					break;
				}
			}
		}
	}
	
	print_diag("update_init 5");
	
	// Makes proposals which make local changes to population-based model
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		if(sp.type == POPULATION){
			// Finds potential local transition combinations
			auto popN_fixed = sp.is_pop_num_fixed();
		
			if(!popN_fixed){
				vector <unsigned int> vec; vec.push_back(p);
				Proposal pp(POP_IC_PROP,vec,model,output,PROP_LOCAL_W,burn_info);
				proposal.push_back(pp);
			}
			
			if(sp.tra_gl.size() > 0){
				{
					vector <unsigned int> vec; vec.push_back(p); 
					Proposal pp(POP_SINGLE_LOCAL_PROP,vec,model,output,PROP_LOCAL_W,burn_info);
					proposal.push_back(pp); 
				}
				
				{
					auto list = sp.find_add_rem_list();
					for(auto &li : list){
						auto vec = li.tr_list;
						vec.push_back(p);
						Proposal pp(POP_ADD_REM_LOCAL_PROP,vec,model,output,PROP_LOCAL_W,burn_info);
						proposal.push_back(pp);
					}
				}

				if(sp.init_cond.type == INIT_POP_DIST){	
					vector <unsigned int> vec; vec.push_back(p);
					Proposal pp(POP_IC_SWAP_PROP,vec,model,output,PROP_LOCAL_W,burn_info);
					proposal.push_back(pp);
				}
				
				{
					vector <unsigned int> vec; vec.push_back(p);
					Proposal pp(POP_MOVE_LOCAL_PROP,vec,model,output,PROP_LOCAL_W,burn_info);
					proposal.push_back(pp);
				}
				
				{
					vector <unsigned int> vec; vec.push_back(p);
					Proposal pp(POP_END_LOCAL_PROP,vec,model,output,PROP_LOCAL_W,burn_info);
					proposal.push_back(pp);
				}
			}
		}
	}
	
	print_diag("update_init 6");
	
	// Does proposals which change parameter and resimulate 
	
	
	// Makes poposals in init cont which change initial condition fractions
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		const auto &ic = sp.init_cond;
		if(ic.type == INIT_POP_DIST){
			vector <unsigned int> vec; vec.push_back(p);
			Proposal pp(INIT_COND_FRAC_PROP,vec,model,output,PROP_LOCAL_W,burn_info);
			proposal.push_back(pp);
		}
	}
	
	print_diag("update_init 7");
	
	if(true){ // Event changes 
		for(auto p = 0u; p < model.nspecies; p++){
			auto &sp = model.species[p];
			auto st = model.samp_type;
			if(sp.type == INDIVIDUAL){	
				vector <unsigned int> vec; vec.push_back(p);
				{				
					Proposal pp(IND_EVENT_TIME_PROP,vec,model,output,1,burn_info);
					if(st == LOCAL_SAMP || st == ALL_SAMP) proposal.push_back(pp); 
					//proposal.push_back(pp); 
				}
				
				{				
					Proposal pp(IND_MULTI_EVENT_PROP,vec,model,output,1,burn_info);
					if(st == LOCAL_SAMP || st == ALL_SAMP) proposal.push_back(pp); 
					//proposal.push_back(pp); 
				}
				
				{				
					Proposal pp(IND_EVENT_ALL_PROP,vec,model,output,1,burn_info);
					if(st == LOCAL_SAMP || st == ALL_SAMP) proposal.push_back(pp); 
					//proposal.push_back(pp); 
				}
				
				for(auto cl = 0u; cl < sp.ncla; cl++){
					if(sp.cla[cl].swap_rep.size() > 0){
						vector <unsigned int> vec2; vec2.push_back(p); vec2.push_back(cl);
						Proposal pp(IND_LOCAL_PROP,vec2,model,output,1,burn_info);
						if(st == LOCAL_SAMP || st == ALL_SAMP) proposal.push_back(pp); 
					}
				}	
				
				if(sp.nindividual_in > 0){	 
					Proposal pp2(IND_OBS_SAMP_PROP,vec,model,output,0.1,burn_info);
					if(st == SAMP_SAMP || st == ALL_SAMP) proposal.push_back(pp2);
					
					Proposal pp3(IND_OBS_RESIM_PROP,vec,model,output,0.1,burn_info);
					if(st == SIM_SAMP || st == ALL_SAMP) proposal.push_back(pp3);
					
					Proposal pp4(IND_OBS_RESIM_SINGLE_PROP,vec,model,output,0.1,burn_info);
					if(st == SIM_CL_SAMP || st == ALL_SAMP) proposal.push_back(pp4);
				}
				
				{
					Proposal pp(IND_UNOBS_RESIM_PROP,vec,model,output,1,burn_info);
					if(st == SIM_SAMP || st == ALL_SAMP) proposal.push_back(pp);
				}
				
				if(sp.contains_source || sp.init_cond.type == INIT_POP_DIST){
					if(sp.trans_tree){ // For trans-tree only add/rem a single individual
						Proposal pp(IND_ADD_REM_TT_PROP,vec,model,output,1,burn_info);
						proposal.push_back(pp);
					}
					else{              // Otherwise add/rem multiple individuals
						Proposal pp(IND_ADD_REM_PROP,vec,model,output,1,burn_info);
						proposal.push_back(pp);
					} 
				}
				
				if(sp.contains_source && sp.init_cond.type != INIT_POP_FIXED){
					Proposal pp(IND_OBS_SWITCH_ENTER_SOURCE_PROP,vec,model,output,1,burn_info);
					proposal.push_back(pp);
				}
				
				if(sp.contains_sink){
					Proposal pp(IND_OBS_SWITCH_LEAVE_SINK_PROP,vec,model,output,1,burn_info);
					proposal.push_back(pp);
				}
			}
		}
	}

	print_diag("update_init 8");	

	if(false){
		// Joint proposals on events and parameters
		for(auto p = 0u; p < model.nspecies; p++){
			auto &sp = model.species[p];
			if(sp.type == INDIVIDUAL){
				for(const auto &pej : sp.par_event_joint){
					vector <unsigned int> vec;
					vec.push_back(p);
					vec.push_back(pej.th);
					for(auto tr : pej.tr_list) vec.push_back(tr);
					switch(pej.dir){
					case FORWARD:
						{
							Proposal pp(PAR_EVENT_FORWARD_PROP,vec,model,output,0.3,burn_info);
							//proposal.push_back(pp);
						}
						break;
						
					case FORWARD_SQ:
						{
							Proposal pp(PAR_EVENT_FORWARD_SQ_PROP,vec,model,output,0.3,burn_info);
							//proposal.push_back(pp);
						}
						break;
						
					case BACKWARD_SQ:
						{
							Proposal pp(PAR_EVENT_BACKWARD_SQ_PROP,vec,model,output,0.3,burn_info);
							//proposal.push_back(pp);
						}
						break;
					}
				}
			}
		}		
	}
	
	
	// Adds a transmission tree proposal
	if(model.trans_tree){
		vector <unsigned int> vec;
		{
			Proposal pp(TRANS_TREE_PROP,vec,model,output,1,burn_info);
			proposal.push_back(pp);
		}
		
		{
			Proposal pp(TRANS_TREE_SWAP_INF_PROP,vec,model,output,1,burn_info);
			proposal.push_back(pp);
		}
		
		if(model.genetic_data.on){
			{
				Proposal pp(TRANS_TREE_MUT_PROP,vec,model,output,1,burn_info);
				proposal.push_back(pp);
			}
			
			{
				Proposal pp(TRANS_TREE_MUT_LOCAL_PROP,vec,model,output,1,burn_info);
				proposal.push_back(pp);
			}
		}
	}
		
	/*
	// Adds proposals to correct observed transitions
	for(auto p = 0u; p < model.nspecies; p++){
		auto &sp = model.species[p];
		if(sp.type == INDIVIDUAL){
			auto fl = false;
			for(auto i = 0u; i < sp.nindividual_obs; i++){
				for(const auto &ob : sp.individual[i].obs){
					switch(ob.type){
					case OBS_SINK_EV: case OBS_TRANS_EV: fl = true; break;
					default: break;
					}
				}
			}
			
			if(fl){			
				vector <unsigned int> vec; vec.push_back(p);
				Proposal pp(CORRECT_OBS_TRANS_PROP,vec,model,output,UNSET,burn_info);
				proposal.push_back(pp);
			}
		}
	}	
	*/	
	 
	print_diag("update_init 9");

	if(model.mode != EXT){
		for(auto i = 0u; i < INF_INIT_NPARAM_SAMPLE; i++){
			auto param_val = model.param_sample(true);
			cor_matrix.add_sample(model.get_param_val_prop(param_val),LARGE);
		}
	}
	
	for(auto &pro : proposal) pro.update_sampler(cor_matrix);
	
	//print_diag("update_init 10");
	
	if(1 == 0) cout << "print_info_turned off in code" << endl;
	else{
		string te = "";
		for(auto &pro : proposal) te += pro.print_info();

		output.prop_summary(te);
	}
	
	print_diag("update_init end");
}


/// Determines if affect type if Markov pop
bool Chain::is_markov_pop(const AffectLike &al) const
{
	switch(al.type){
	case MARKOV_POP_AFFECT:
		return true;
	default: break;
	}
	
	return false;
}	


/// Adds a parameter proposal
void Chain::add_parameter_prop(const vector <unsigned int> &vec, bool log_trans)
{
	//cout << "LOG TRANS" << endl; if(log_trans) return; 
	
	auto bern_fl = false;
	for(auto th : vec){
		if(model.is_bernoulli(th)) bern_fl = true;
	}
	
	if(model.species[0].type == DETERMINISTIC){
		auto type = DET_PARAM_PROP;
		if(log_trans) type = DET_LOG_PARAM_PROP;
		if(bern_fl){
			type = DET_BERNOULLI_PROP;
			if(vec.size() != 1) emsg("Cannot update more than one bernoulli");
		}
	
		Proposal pp(type,vec,model,output,1,burn_info);
		//for(auto p = 0u; p < model.species.size(); p++){
			//model.add_like_obs_affect(p,pp.affect_like);
		//}
		proposal.push_back(pp);
		return;
	}

	auto type = PARAM_PROP;
	if(log_trans) type = LOG_PARAM_PROP;
	if(bern_fl){
		type = BERNOULLI_PROP;
		if(vec.size() != 1) emsg("Cannot update more than one bernoulli");
	}
	
	Proposal pp(type,vec,model,output,1,burn_info);
	proposal.push_back(pp);
	
	/// Looks at adding MBPs
	vector <unsigned int> list;
	HashSimp hash_list;
	for(const auto &al : pp.affect_like){
		if(is_markov_pop(al)) add_to_vec(list,al.num,hash_list);
	}

	for(auto k = 0u; k < list.size(); k++){    // Adds one MBP per species affected
		auto pp2 = pp;
		pp2.p_prop = list[k];
		
		auto j = 0u;                             // Removes MARKOV_POP_AFFECT etc...
		while(j < pp2.affect_like.size()){
			auto al = pp2.affect_like[j];
			
			auto fl = false;
		
			if(al.num == list[k]) fl = is_markov_pop(al);
			
			switch(al.type){
			case MARKOV_VALUE_AFFECT: fl = true; break;
			default: break;
			}
			
			if(fl) pp2.affect_like.erase(pp2.affect_like.begin()+j);
			else j++;
		}

		model.add_like_obs_affect(pp2.p_prop,pp2.affect_like);
	
		pp2.mbp_population_affect();
	
		if(bern_fl){
			if(vec.size() != 1) emsg("Cannot update more than one bernoulli");
			pp2.type = MBP_BERNOULLI_PROP;
			pp2.name = "MBP BERNOULLI "+pp2.name;
			
			proposal.push_back(pp2);
		}
		else{
			pp2.type = MBP_PROP; if(log_trans) pp2.type = LOG_MBP_PROP;
			pp2.name = "MBP "+pp2.name;
			proposal.push_back(pp2);
		}
	}
}


/// Used to order affects
bool prop_time_ord(const PropTime &pt1, const PropTime &pt2)                      
{ return (pt1.time > pt2.time); };  


/// Generates a banner
string Chain::banner(string te) const 
{
	stringstream ss;
	
	auto len = 70;
	auto num = (len-te.length())/2-1;
	
	for(auto i = 0u; i < num; i++) ss << "*";
	ss << " " << te << " ";
	for(auto i = 0u; i < len-num-te.length(); i++) ss << "*";
	
	return ss.str();
}


/// Outputs diagnostics about the proposals
string Chain::diagnostics(double total_time, double init_time, double anneal_time) const
{                    
	stringstream ss;
	
	auto num = state.alg_warn.size();
	for(const auto &ssp : state.species) num += ssp.alg_warn.size();
	
	if(num > 0){
		ss << ALG_WARN << ":" << endl;
		for(const auto &aw : state.alg_warn){
			ss << aw.te << " (core: " << aw.core << ", sample: ";
			if(aw.sample == UNSET) ss << "unset"; else ss << aw.sample;
			ss << ", number: " << aw.num << ")" << endl;
		}
		
		for(auto p = 0u; p < model.species.size(); p++){
			const auto &ssp = state.species[p];
			if(ssp.alg_warn.size() > 0){
				ss << "Warning in species " << model.species[p].name << endl;
				for(const auto &aw : ssp.alg_warn){
					ss << aw.te << " (core: " << aw.core << ", number: " << aw.num << ")" << endl;
				}
			}
		}		
		ss << endl;
	}
	
	{	
		auto sum_param_prop = 0.0, sum_event_prop = 0.0;
		auto sum_joint_prop = 0.0, sum_tree_prop = 0.0;
		auto sum_param_prop_precalc = 0.0, sum_event_prop_precalc = 0.0;
		auto sum_joint_prop_precalc = 0.0, sum_tree_prop_precalc = 0.0;
		auto sum_prop = 0.0;
		
		auto sum_resample = 0.0, sum_ie = 0.0, sum_pop = 0.0, sum_timedep = 0.0, sum_calc = 0.0, sum_t_precalc = 0.0;
		for(auto &pro : proposal){
			auto t = pro.timer[PROP_TIMER];
			auto t_precalc = pro.timer[PARAM_RESAMPLE_TIMER]+pro.timer[PARAM_IE_TIMER]+pro.timer[PARAM_POP_TIMER]+pro.timer[PARAM_TIMEDEP_TIMER]+pro.timer[PARAM_CALC_TIMER];
			
			sum_resample += pro.timer[PARAM_RESAMPLE_TIMER];
			sum_ie += pro.timer[PARAM_IE_TIMER];
			sum_pop += pro.timer[PARAM_POP_TIMER];
			sum_timedep += pro.timer[PARAM_TIMEDEP_TIMER];
			sum_calc += pro.timer[PARAM_CALC_TIMER];
			sum_t_precalc += t_precalc;
		
			sum_prop += t;
			
			switch(pro.type){
			case PARAM_PROP: case LOG_PARAM_PROP: case BERNOULLI_PROP: 
			case DET_PARAM_PROP: case DET_LOG_PARAM_PROP: case DET_BERNOULLI_PROP: 
			case DET_IC_POP_PROP: case DET_IC_POPTOTAL_PROP: case DET_IC_RESAMP_PROP:
			case INIT_COND_FRAC_PROP: case IE_PROP: case IE_VAR_PROP:
			case IE_COVAR_PROP: case IE_VAR_CV_PROP:
				sum_param_prop += t;
				sum_param_prop_precalc += t_precalc;
				break;
			
			case IND_EVENT_TIME_PROP: case IND_MULTI_EVENT_PROP: case IND_EVENT_ALL_PROP:
			case IND_OBS_SAMP_PROP: case IND_OBS_RESIM_PROP: case IND_OBS_RESIM_SINGLE_PROP:
			case IND_UNOBS_RESIM_PROP: case IND_ADD_REM_PROP: case IND_ADD_REM_TT_PROP:
			case POP_ADD_REM_LOCAL_PROP: case POP_MOVE_LOCAL_PROP: case POP_IC_LOCAL_PROP:
			case POP_END_LOCAL_PROP: case POP_SINGLE_LOCAL_PROP: case POP_IC_PROP:
			case POP_IC_SWAP_PROP:
			case IND_LOCAL_PROP: case CORRECT_OBS_TRANS_PROP:
			case IND_OBS_SWITCH_ENTER_SOURCE_PROP: case IND_OBS_SWITCH_LEAVE_SINK_PROP:
			case MBPII_PROP: case MBP_IC_POP_PROP:
			case MBP_IC_POPTOTAL_PROP: case MBP_IC_RESAMP_PROP:
				sum_event_prop += t;
				sum_event_prop_precalc += t_precalc;
				break;
				
			case MBP_PROP: case LOG_MBP_PROP: case MBP_BERNOULLI_PROP: 
			case PAR_EVENT_FORWARD_PROP: case PAR_EVENT_FORWARD_SQ_PROP:
			case PAR_EVENT_BACKWARD_SQ_PROP:
				sum_joint_prop += t;
				sum_joint_prop_precalc += t_precalc;
				break;
			
			case TRANS_TREE_PROP: case TRANS_TREE_SWAP_INF_PROP: case TRANS_TREE_MUT_PROP:
			case TRANS_TREE_MUT_LOCAL_PROP:
				sum_tree_prop += t;
				sum_tree_prop_precalc += t_precalc;
				break;
			}
		}
		
		
		ss << banner("PROPOSALS TYPES") << endl;
		ss << "Shows how run-time CPU time is distributed into different proposal types:" << endl << endl;
			
		if(sum_param_prop != 0){
			ss << "Parameters - " << cpu_percent(sum_param_prop,total_time) << " (" << cpu_percent(sum_param_prop_precalc,total_time) << " from sampling)" << endl;
		}
		
		if(sum_event_prop != 0){
			ss << "Events - " << cpu_percent(sum_event_prop,total_time);

			ss << " (ind update:" << cpu_percent(state.timer[IND_TIMER],total_time) << " of which " << cpu_percent(state.timer[IND_POP_UPDATE_TIMER],total_time) << " is from population)" << " (" << cpu_percent(sum_event_prop_precalc,total_time) << " from sampling)" << endl;

		}
		
		if(sum_joint_prop != 0){
			ss << "Joint event/parameter - " << cpu_percent(sum_joint_prop,total_time) << " (" << cpu_percent(sum_joint_prop_precalc,total_time) << " from sampling)" << endl;
		}

		if(sum_tree_prop != 0){
			ss << "Transmission tree - " << cpu_percent(sum_tree_prop,total_time) << " (" << cpu_percent(sum_tree_prop_precalc,total_time) << " from sampling)" << endl;
		}
		
		auto other = double(total_time)-double(sum_prop);
		ss << "Non-proposal - " << cpu_percent(other,total_time);
		ss << " (samplers: " << cpu_percent(state.timer[UPDATE_SAMPLER_TIMER],total_time);
    ss << ", o/p params: " << cpu_percent(output.timer[PARAM_OUTPUT],total_time);
		ss << ", o/p state: " << cpu_percent(output.timer[STATE_OUTPUT],total_time);
		ss << ", der-pre: " << cpu_percent(state.timer[DERIVE_PRECALC_TIMER],total_time);
		ss << ", derive: " << cpu_percent(state.timer[DERIVE_TIMER],total_time);
		ss << ", checking: " << cpu_percent(state.timer[CHECK_TIMER],total_time);
	
		ss << ")" << endl;
		
		ss << "Initialisation percentage - " 
				<< cpu_percent(init_time,init_time+total_time) << endl;
				
		if(anneal_time != UNSET){
			ss << "Annealing percentage - " 
				<< cpu_percent(anneal_time,total_time) << endl << endl;
		}

		ss << endl;
		ss << "Sampling split - ";
		ss << "resample " << cpu_percent(sum_resample,sum_t_precalc) << ", ";
		ss << "ie " << cpu_percent(sum_ie,sum_t_precalc) << ", ";
		ss << "pop " << cpu_percent(sum_pop,sum_t_precalc) << ", ";
		ss << "timedep " << cpu_percent(sum_timedep,sum_t_precalc) << ", ";
		ss << "calc " << cpu_percent(sum_calc,sum_t_precalc);
		ss << endl;
	
		ss << endl;
		ss <<  "Total CPU running time: " << get_cpu_time(total_time) << endl;
	
	  ss << endl;
		
		if(!model.deterministic){
			ss << banner("CHECKING") << endl;
			ss << "Shows how run-time CPU checking time is distributed:" << endl << endl;
			auto ttot = state.timer[CHECK_TIMER];
			auto fl = false;
			for(auto j = 0u; j < CHECK_MAX-1; j++){
				auto t = state.check_timer[j];
				if(t >= 0.01*ttot){
					if(fl == true) ss << ", ";
					ss << check_name[j] << " "  << cpu_percent(t,ttot);
					fl = true;
				}
			}

			ss << endl << endl;
		}
	}
	
	{
		ss << banner("PARAMETER UPDATE") << endl;
		ss << "Determines which parts of parameter update use most CPU" << endl << endl;
		
		//auto sum_restore = 0.0;
		for(auto i = 0u; i < AFFECT_MAX; i++){	
			auto update_num = cpu_percent(state.update_timer[i],total_time);
			auto restore_num = cpu_percent(state.restore_timer[i],total_time);
			
			if(!(update_num == "0%" && restore_num == "0%")){
				switch(i){
				case SPLINE_PRIOR_AFFECT: ss << "Spline prior"; break; 
				case PRIOR_AFFECT: ss << "Prior"; break; 
				case IEG_PRIOR_AFFECT: ss << "IEG Prior"; break; 
				case DIST_AFFECT: ss << "Distribution"; break; 
				case EXP_FE_AFFECT: ss << "Exp FE"; break; 
				case MARKOV_VALUE_AFFECT: ss << "Markov value"; break; 
				case MARKOV_LIKE_AFFECT: ss << "Markov like"; break; 
				//case POP_AFFECT: ss << "Population"; break; 
				case NM_TRANS_AFFECT: ss << "NM trans"; break; 
				case NM_TRANS_BP_AFFECT: ss << "NM trans bp"; break; 
				case NM_TRANS_INCOMP_AFFECT: ss << "NM trans incomp"; break; 
				case OMEGA_AFFECT: ss << "Omega"; break; 
				case EXP_IE_AFFECT: ss << "Exp IE"; break; 
				case LIKE_IE_AFFECT: ss << "Like IE"; break; 
				case INDFAC_INT_AFFECT: ss << "Individual factor"; break; 
				case MARKOV_POP_AFFECT: ss << "Markov pop"; break; 
				case OBS_EQN_AFFECT: ss << "Obs Eqn"; break; 
				case LIKE_OBS_IND_AFFECT: ss << "Like Obs Ind"; break; 
				case LIKE_OBS_POP_AFFECT: ss << "Like Obs Pop"; break; 
				case LIKE_OBS_POP_TRANS_AFFECT: ss << "Like Obs Pop Trans"; break; 
				case LIKE_UNOBS_TRANS_AFFECT: ss << "Like unobs trans"; break;
				case POP_DATA_CGL_TGL_AFFECT: ss << "Pop data cgl tgl"; break;
				case LIKE_INIT_COND_AFFECT: ss << "Like init cond"; break;
				case PRIOR_INIT_COND_AFFECT: ss << "prior init cond"; break;
				case LIKE_GENETIC_PROCESS_AFFECT: ss << "Like genetic process"; break; 
				case GENETIC_VALUE_AFFECT: ss << "Genetic value"; break; 
				case LIKE_GENETIC_OBS_AFFECT: ss << "Like genetic obs"; break; 
				}
				ss << "  ";
				ss << "update - " << update_num << "     ";
		
				ss << "restore - " << restore_num << endl;
			}
			
			//sum_update += state.update_timer[i];
			//sum_restore += state.restore_timer[i];
		}
		
		ss << endl << endl;
	}
	
	{
		vector <PropTime> prop_time;
		for(auto i = 0u; i < proposal.size(); i++){
			PropTime pt; pt.i = i; pt.time = proposal[i].timer[PROP_TIMER];
			prop_time.push_back(pt); 
		}
		
		sort(prop_time.begin(),prop_time.end(),prop_time_ord);
		
		ss << banner("MCMC PROPOSAL DIAGNOSTICS") << endl;
		ss << "Diagnostics for different proposals (ordered by CPU time):";
		ss << endl << endl;
		
		for(const auto &pt : prop_time){
			const auto &pro = proposal[pt.i];
			ss << pro.diagnostics(total_time);
			ss << endl;
		}
	}
	
	ss << endl;

	if(false){
		for(const auto &ssp : state.species){
			ss << ssp.source_sampler.print();
		}
	}
	
	return ss.str();
}


/// Checks if new proposals are added which join parameters
void Chain::check_join_proposal()
{
	auto pl = false;
	
	auto n = cor_matrix.get_n();
	
	auto nn_max = 400.0;
	auto nn = n-INF_INIT_NPARAM_SAMPLE;
	if(nn > nn_max) nn = nn_max;
	
	auto num = 1.0;
	auto ndecay = nn_max/num;
	auto sh = exp(-num);
	auto f = (exp(-(nn/ndecay))-sh)/(1-sh);
	
	auto thresh = f + (1-f)*PROP_JOIN_COR_MIN;
	if(n < 50) thresh = 1;
	
	// Switches off all multinomial proposals
	vector <bool> on_st;
	for(auto &prop : proposal){	
		on_st.push_back(prop.on);
		if(prop.is_param_prop() && prop.param_list.size() > 1){
			prop.on = false;
		}
	}
	
	auto M_norm = cor_matrix.calculate_cor_matrix(false);
	auto M_log_trans = cor_matrix.calculate_cor_matrix(true);
	
	{  // This looks for whatever has the most extreame correlation (log or normal) to construct multivariable proposals
		const auto &list = cor_matrix.str_pos_list;
		for(auto j : list){
			for(auto i : list){
				if(j != i){
					auto cor1 = M_norm[j][i]; if(cor1 < 0) cor1 *= -1;
					auto cor2 = M_log_trans[j][i]; if(cor2 < 0) cor2 *= -1;
					if(cor1 > cor2) M_log_trans[j][i] = 0;
					else M_norm[j][i] = 0;
				}
			}
		}
	}
	
	//print_matrix("M",M_norm);
	//print_matrix("M_log",M_log_trans);

	for(auto loop = 0u; loop < 2; loop++){ // Goes through normal then log transformmed
		if(pl) cout << loop << " loop" << endl;	
		auto log_trans = false; 
		auto M = M_norm;
		
		if(loop == 1){ log_trans = true; M = M_log_trans;}
		
		// Works out all multivariate proposals
		vector <unsigned int> pvp_list; // List of parameter to be considered
		
		vector < vector <unsigned int> > par_list;
		for(auto i = 0u; i < model.nparam_vec_prop; i++){	
			auto th = model.param_vec_prop[i];
			if(!model.is_bernoulli(th)){
				if(!(loop == 1 && !model.param_vec[th].strictly_positive)){
					pvp_list.push_back(i);
					
					vector <unsigned int> list;
					list.push_back(th);
					par_list.push_back(list);
				}
			}
		}
		
		for(auto j : pvp_list){
			for(auto i: pvp_list){
				if(i > j){
					if(M[j][i] > thresh || M[j][i] < -thresh){
						auto th_j = model.param_vec_prop[j];
						auto th_i = model.param_vec_prop[i];
						
						auto lj = find_which_list(th_j,par_list);
						auto li = find_which_list(th_i,par_list);
						if(lj != li){
							for(auto va : par_list[lj]){
								par_list[li].push_back(va);
							}
									
							if(lj+1 != par_list.size()){
								par_list[lj] = par_list[par_list.size()-1];
							}
							par_list.pop_back();
						}
					}
				}
			}
		}
		
		// Switch on those which are in list
		vector <unsigned int> add;
		for(auto j = 0u; j < par_list.size(); j++){
			auto fl = false;
			for(auto i = 0u; i < proposal.size(); i++){ // Tries to re-turn on an existing proposal
				auto &prop = proposal[i];
				if(prop.is_param_prop()){
					if(equal_vec(prop.param_list,par_list[j]) && prop.log_trans == log_trans){
						fl = true;
						prop.on = true;
						if(on_st[i] == false){
							if(pl) cout << "RE TURN ON: " << prop.name << endl;
							prop.update_sampler(cor_matrix);
						}
					}
				}
			}
			
			if(fl == false){
				add.push_back(j);
			}
		}
		
		for(auto i = 0u; i < proposal.size(); i++){
			auto &prop = proposal[i];
			if(on_st[i] && prop.on == false){
				if(pl) cout << "REMOVE: " << prop.name << endl;
			}
		}
		
		// Adds any new proposals
		auto pr_st = proposal.size();
		for(auto j : add) add_parameter_prop(par_list[j],log_trans);
		
		for(auto j = pr_st; j < proposal.size(); j++){
			auto &prop = proposal[j];
			if(pl) cout << "JOIN ADD: " << prop.name << endl;
			prop.update_sampler(cor_matrix);
		}
		
		if(pl){
			cout << "after" << endl;
			for(auto &prop : proposal){	
				cout << prop.name << " " << prop.on << endl;
			}
		}
	}
}


/// Finds which list a parameter is on
unsigned int Chain::find_which_list(unsigned int j, const vector < vector <unsigned int> > &par_list) const 
{
	for(auto k = 0u; k < par_list.size(); k++){
		for(auto i = 0u; i < par_list[k].size(); i++){
			if(par_list[k][i] == j) return k;
		}
	}
	emsg("Cannot find on list");
	
	return UNSET;
}


/// Check correlation matrix is specified correctly
void Chain::check_cor_matrix() const
{
	cor_matrix.check();
}


/// Gathers all the infomation used in proposals
vector <PropInfo> Chain::get_prop_info() const
{
	vector <PropInfo> prop_info;
	
	for(auto i = 0u; i < proposal.size(); i++){
		const auto &pro = proposal[i];
		if(pro.on && pro.prop_info_on()){
			prop_info.push_back(pro.get_prop_info());
		}
	}
	
	return prop_info;
}


/// Sets the propsal info
void Chain::set_prop_info(const vector <PropInfo> &prop_info)
{
	for(auto i = 0u; i < proposal.size(); i++){
		auto &pro = proposal[i];
		if(pro.on && pro.prop_info_on()){
			auto id = pro.get_prop_id();
			
			auto j = 0u; 
			while(j < prop_info.size() && !(prop_info[j].type == pro.type && equal_vec(prop_info[j].id,id))) j++;
			if(j == prop_info.size()){
				emsg("Coudl not find prop_sync"+pro.name+" "+print_prop_info(prop_info,model));
			}
			pro.set_prop_info(prop_info[j]);
		}
	}
}


/// Gets terminal information from chain
TerminalInfo Chain::get_terminal_info(unsigned int ch_tot) const 
{
	TerminalInfo ti;
	ti.ch = ch_tot;
	ti.prop_info_store = get_prop_info();
	ti.n = cor_matrix.n; 
	ti.n_start = cor_matrix.n_start; 
	ti.av = cor_matrix.av;
	ti.av2 = cor_matrix.av2;
	ti.log_av = cor_matrix.log_av;
	ti.log_av2 = cor_matrix.log_av2;
	
	return ti;
}

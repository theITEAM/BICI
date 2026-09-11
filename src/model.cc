// This file defines the compartmental model 
 
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
using namespace std;

#include "model.hh"
#include "utils.hh"
#include "matrix.hh"

/// Initialises the model 
Model::Model(Operation mode_, ExtFactor ext_factor_, bool no_question_) : precalc_eqn(species_simp,spline,param,param_vec,pop,constant,timepoint,details)
{
	mode = mode_; 
	ext_factor = ext_factor_;	
	no_question = no_question_;
	sync_on = true;
	deterministic = false;
	nspecies = 0;
	genetic_data.on = false;
};


/// Adds an equation reference to EquationInfo
void Model::add_eq_ref(EquationInfo &eqi, Hash &hash_eqn, double tdiv, bool keep_te)
{
	auto ti = UNSET;
	if(tdiv != UNSET) ti = get_ti(tdiv);
	
	auto te_st = eqi.te;
	
	if(te_st == "") emsg_input("Equation does not have any text");
	
	auto vec = hash_eqn.get_vec_eqn(te_st,(unsigned int)eqi.type,eqi.p,eqi.cl,eqi.infection_trans,ti);

	auto e = hash_eqn.existing(vec);
	if(e != UNSET){
		eqi.eq_ref = e;
	}
	else{
		eqi.eq_ref = eqn.size();	
		
		hash_eqn.add(eqi.eq_ref,vec);
	
		Equation eq(eqi,ti,species_simp,comp_pos,hash_comp_pos,param,prior,region,derive,spline,param_vec,density,pop,hash_pop,popcomb,popcombw,constant,timepoint,details,define);
		
		if(false && eq.warn != ""){ 
			cout << "TURN OFF" << endl;
			cout << eq.warn << endl; 
			emsg_input("warning");
		}
		
		eqn.push_back(eq);
	}
	
	if(keep_te) eqi.te = te_st;
}


/// Initialises a parameter value
void Model::param_val_init(PV &param_val) const
{
	auto N = nparam_vec;
	auto M = precalc_eqn.pcsize;
	param_val.value.resize(N,UNSET); 
	param_val.precalc = precalc_init;
	param_val.value_old.resize(N,UNSET);
	param_val.precalc_old.resize(M,UNSET);
}


/// Samples a set of parameter values from the model
PV Model::param_sample(bool no_precalc) const
{
	PV param_val;
	param_val_init(param_val);
	
	auto &value = param_val.value;
	auto &precalc = param_val.precalc;
	
	for(auto th = 0u; th < nparam_vec; th++){
		const auto &pv = param_vec[th];
		
		precalc_eqn.calculate(pv.spec_precalc_before,param_val,false);
	
		const auto &par = param[pv.th];
	
		switch(par.variety){
		case CONST_PARAM:
			value[th] = par.get_value(pv.index);
			break;
		
		case REPARAM_PARAM: 
			{
				if(pv.reparam_time_dep == false){
					auto eq_ref = par.get_eq_ref(pv.index);
					if(eq_ref == UNSET) emsg("eq_ref should be set");
					value[th] = eqn[eq_ref].calculate_param(precalc);
				}
			}
			break;
			
		case DIST_PARAM: case PRIOR_PARAM:	
			value[th] = prior_sample(prior[pv.prior_ref],precalc);
			break;
		
		case DYNAMIC_PARAM: emsg("Should not sample dynamic param"); break;
		
		case UNSET_PARAM: emsg("error param"); break;
		}
		
		if(!pv.reparam_time_dep){
			precalc_eqn.calculate(pv.set_param_spec_precalc,param_val,false);
		}		
	}

	sample_ieg_cv(param_val); // Samples any individual effect covariances 
	
	if(no_precalc) return param_val;
	
	precalc_eqn.calculate(spec_precalc_sample,param_val,false);
	
	//print_param(param_val);
	
	if(false){
		for(auto val: param_val.value) cout << val << ","; 	
		cout << " value" << endl;
		
		//for(auto val: param_val.precalc) cout << val << ","; 	
		auto imax = param_val.precalc.size();
		if(imax > 1000) imax = 1000;
		for(auto i = 0u; i < imax; i++){
			cout << i << " " << param_val.precalc[i] << " prec" << endl;
		}
		if(imax == 1000) cout << "..." << endl;
		cout << " precalc" << endl;
		emsg("sample");
	}

	return param_val;
}


/// Generates Latin Hypercube sampling for prior distribution with constant hyper paerameters
vector < vector <double> > Model::param_latin_hypercube_prior_sample(unsigned int N) const
{
	vector < vector <double> > param_latin;
	
	auto M = param_vec_latin.size();
	
	vector < vector <double> > par_samp;
	for(auto i = 0u; i < M; i++){
		auto th = param_vec_latin[i];
		const auto &pv = param_vec[th];
		par_samp.push_back(latin_prior_sample(prior[pv.prior_ref],N));
	}
	
	auto store = par_samp;
	
	for(auto j = 0u; j < N; j++){
		vector <double> vec;
		for(auto i = 0u; i < M; i++){
			auto ns = par_samp[i].size();
			auto k = (unsigned int)(ran()*ns);
			vec.push_back(par_samp[i][k]);
			
			if(k+1 < ns) par_samp[i][k] = par_samp[i][ns-1];
			par_samp[i].pop_back();
		}
		param_latin.push_back(vec);
	}
	
	if(false){
		for(auto j = 0u; j < N; j++){
			cout << j << ": ";
			for(auto i = 0u; i < M; i++){
				auto th = param_vec_latin[i];
				cout << param_vec_name(th) << " " << store[i][j] << ",";
			}
			cout << " before" << endl;
		}
		
		for(auto j = 0u; j < N; j++){
			cout << j << ": ";
			for(auto i = 0u; i < M; i++){
				auto th = param_vec_latin[i];
				cout << param_vec_name(th) << " " << param_latin[j][i] << ",";
			}
			cout << " latin" << endl;
		}
	}
	
	return param_latin;
}


/// Samples a set of parameter values from the model using some priors set by latin hypercube
PV Model::param_latin_hypercube_sample(const vector <double> &prior_val) const
{
	PV param_val;
	param_val_init(param_val);
	
	auto &value = param_val.value;
	auto &precalc = param_val.precalc;
	
	auto M = param_vec_latin.size();
	if(prior_val.size() != M) emsg("Problemn witrh prior_val");
	
	for(auto i = 0u; i < M; i++){
		value[param_vec_latin[i]] = prior_val[i];
	}
	
	for(auto th = 0u; th < nparam_vec; th++){
		const auto &pv = param_vec[th];
		
		precalc_eqn.calculate(pv.spec_precalc_before,param_val,false);
	
		const auto &par = param[pv.th];
	
		switch(par.variety){
		case CONST_PARAM:
			value[th] = par.get_value(pv.index);
			break;
		
		case REPARAM_PARAM: 
			{
				if(pv.reparam_time_dep == false){
					auto eq_ref = par.get_eq_ref(pv.index);
					if(eq_ref == UNSET) emsg("eq_ref should be set");
					value[th] = eqn[eq_ref].calculate_param(precalc);
				}
			}
			break;
			
		case DIST_PARAM: case PRIOR_PARAM:	
			if(pv.latin_sample){
				if(value[th] == UNSET) emsg("Latin should be set");
			}
			else value[th] = prior_sample(prior[pv.prior_ref],precalc);
			break;
		
		case DYNAMIC_PARAM: emsg("Should not sample dynamic param"); break;
		
		case UNSET_PARAM: emsg("error param"); break;
		}
		
		if(!pv.reparam_time_dep){
			precalc_eqn.calculate(pv.set_param_spec_precalc,param_val,false);
		}		
	}

	sample_ieg_cv(param_val); // Samples any individual effect covariances 
	
	precalc_eqn.calculate(spec_precalc_sample,param_val,false);
	
	//print_param(param_val);
	
	if(false){
		for(auto val: param_val.value) cout << val << ","; 	
		cout << " value" << endl;
		
		//for(auto val: param_val.precalc) cout << val << ","; 	
		auto imax = param_val.precalc.size();
		if(imax > 1000) imax = 1000;
		for(auto i = 0u; i < imax; i++){
			cout << i << " " << param_val.precalc[i] << " prec" << endl;
		}
		if(imax == 1000) cout << "..." << endl;
		cout << " precalc" << endl;
		emsg("latin sample");
	}

	return param_val;
}


/// Determines if bounded priors are sampled or they use a uniform distribution
// For PAS and DA it makes sense to use a uniform, such that there is diversity in initial state
bool Model::sample_bounded() const
{
	switch(details.algorithm){
	case DA_MCMC: case PAS_MCMC: return false;
	default: return true;	
	}
}


/// Samples covariance matrices for indidivudal effects using Wishart distribution
void Model::sample_ieg_cv(PV &param_val) const
{
	auto sbound = sample_bounded();
	
	for(auto p = 0u; p < species.size(); p++){
		const auto &sp = species[p];
		for(auto g = 0u; g < sp.ind_eff_group.size(); g++){
			const auto &ieg = sp.ind_eff_group[g];
			
			const auto &par = param[ieg.th];
			if(par.variety == PRIOR_PARAM){
				const auto &pri = prior[ieg.prior_ref];
				switch(pri.type){
				case MVN_DEFAULT_PR: case MVN_NORM_LKJ_PR: case MVN_UNIFORM_LKJ_PR:
				case MVN_INV_WISH_PR: case MVN_JEF_PR: case MVN_UNIFORM_PR:
					break;
					
				default:
					run_error("Should be a covariance prior");
					break;
				}
				
				auto N = ieg.list.size();
				
				switch(pri.type){
				case MVN_DEFAULT_PR: case MVN_NORM_LKJ_PR: case MVN_UNIFORM_LKJ_PR:
					{
						string warn;
						for(auto j = 0u; j < N; j++){
							double val;
							if(pri.type == MVN_UNIFORM_LKJ_PR){
								val = ieg.var_min + ran()*(ieg.var_max-ieg.var_min);
							}
							else{
								val = normal_sample(0,ieg.lkj_sd,warn);
								if(warn != "") emsg("Problem with sampling covariance matrix");
								if(val < 0) val *= -1;
							}
							
							auto th = ieg.omega_pv[j][j];
							param_val.value[th] = val;
							precalc_eqn.calculate(param_vec[th].set_param_spec_precalc,param_val,false);
						}

						auto R = LKJ_sample(ieg.lkj_eta,N);
						
						for(auto j = 0u; j < N; j++){
							for(auto i = 0u; i < j; i++){
								auto th = ieg.omega_pv[j][i];
								param_val.value[th] = R[j][i];
								precalc_eqn.calculate(param_vec[th].set_param_spec_precalc,param_val,false);
							}
						}
					
						auto omega = sp.calculate_omega_basic(g,param_val,param);
						auto det = determinant_fast(omega);
						if(det == UNSET) emsg("Sampling determinant problem");
					}
					break;
					
				case MVN_INV_WISH_PR:
					{
						auto M = inv_wishart_sample(ieg.inv_wish_S,ieg.inv_wish_nu,N);
						for(auto j = 0u; j < N; j++){
							for(auto i = 0u; i <= j; i++){
								auto th = ieg.omega_pv[j][i];
								auto val = M[j][i];
								if(j != i) val /= sqrt(M[j][j]*M[i][i]);
								param_val.value[th] = val;
								precalc_eqn.calculate(param_vec[th].set_param_spec_precalc,param_val,false);
							}
						}
					
						auto omega = sp.calculate_omega_basic(g,param_val,param);
						auto det = determinant_fast(omega);
						if(det == UNSET) emsg("Sampling determinant problem");
					}
					break;
				
				case MVN_JEF_PR: case MVN_UNIFORM_PR:
					{
						auto loop = 0u, loopmax = 1000u;
						for(loop = 0; loop < loopmax; loop++){
							if(sbound){
								emsg("Jeffreys problem");
								auto det = 0.0;
								auto det_min = exp(ieg.log_det_min);
								//auto det_max = exp(  ieg.log_det_max);
								auto det_max = UNSET;//exp(  ieg.log_det_max);
						
								switch(pri.type){
								case MVN_UNIFORM_PR:
									det = det_min+ran()*(det_max-det_min);
									break;
									
								case MVN_JEF_PR:
									if(N == 1) det = det_min*exp(ran()*log(det_max/det_min));
									else{
										auto power = -0.5*N-0.5;				
										det = pow(pow(det_min,power+1)+ran()*(pow(det_max,power+1)-pow(det_min,power+1)),1.0/(power+1));
									}
									break;
									
								default: emsg("op er"); break;
								}
								
								auto var = pow(det,1.0/N);
								
								for(auto j = 0u; j < N; j++){
									for(auto i = 0u; i < N; i++){
										auto valu = 0.0; if(i == j) valu = var;
										
										auto th = ieg.omega_pv[j][i];
										param_val.value[th] = valu;
										precalc_eqn.calculate(param_vec[th].set_param_spec_precalc,param_val,false);
									}
								}
							}
							else{
								for(auto j = 0u; j < N; j++){
									for(auto i = 0u; i < N; i++){
										auto valu = 0.0; 
										if(i == j) valu = ran()*(ieg.var_max);
										else valu = COR_MAX*(2*ran()-1);
											
										auto th = ieg.omega_pv[j][i];
										param_val.value[th] = valu;
										
										precalc_eqn.calculate(param_vec[th].set_param_spec_precalc,param_val,false);
									}
								}
							}
							
							auto omega = sp.calculate_omega_basic(g,param_val,param);
							
							auto det = determinant_fast(omega);
							if(det != UNSET && det > ieg.log_det_min) break;
						}
						
						if(loop == loopmax){
							run_error("Could not sample a covariance matrix.");
						}
					}
					break;
					
				default: emsg("wrong distribution"); break;
				}
			}			
		}
	}
}

/// Updates precalc for all times
void Model::param_spec_precalc_time_all(const vector < vector <double> > &popcomb_t, PV &param_val, bool store) const
{
	for(auto ti = 0u; ti < details.T; ti++){
		param_spec_precalc_time(ti,popcomb_t,param_val,store);
	}
}

	
/// Updates parameter precalculation at given time point
void Model::param_spec_precalc_time(unsigned int ti, const vector < vector <double> > &popcomb_t, PV &param_val, bool store) const
{
	auto ref = spec_precalc_time_ref[ti];
	if(ref == UNSET) return;

	const auto &upt = spec_precalc_list[ref];
	
	auto &value = param_val.value;
	auto &precalc = param_val.precalc;

	for(auto th : upt.pv){
		const auto &pv = param_vec[th];
	
		const auto &par = param[pv.th];
		auto eq_ref = par.get_eq_ref(pv.index);
		if(eq_ref == UNSET) emsg("eq_ref should be set");
			
		precalc_eqn.calculate(pv.spec_precalc_before,param_val,false);
	
		if(store) param_val.value_change(th);
		
		value[th] = eqn[eq_ref].calculate_all_time(ti,popcomb_t,precalc);
	
		//precalc_eqn.print_spec_precalc("spec",pv.set_param_spec_precalc);
		precalc_eqn.calculate(pv.set_param_spec_precalc,param_val,store);
	}
	
	for(const auto &dpu : upt.dynamic_param_update){
		auto th = dpu.th;
		const auto &par = param[th];
		if(par.variety == DYNAMIC_PARAM){
			auto s_start = par.spline_ref;
			
			auto ntimes = par.spline_info.knot_tdiv.size();
			auto N = par.N/ntimes;
		
			auto ti = dpu.ti;
			auto ti_end = dpu.ti_end;
		
			for(auto j = 0u; j < N; j++){
				auto s = s_start+j;				
			
				auto i = spline_refi[s];
				
				auto eq_ref = par.get_eq_ref(j*ntimes);
			
				if(eq_ref == UNSET) emsg("eq_ref should be set");
				
				//time -= clock();
				auto val = eqn[eq_ref].calculate_all_time(ti,popcomb_t,precalc);
				//time += clock();
				
				const auto &di = par.dynamic_info;
				switch(di.type){
				case BIN_THRESH: case BIN_THRESH_DIST: case BIN_THRESH_REGION: case BIN_THRESH_EQN:
					{
						auto thresh = eqn[di.thresh.eq_ref].calculate_all_time(ti,popcomb_t,precalc);
						if(thresh < 0){
							run_error("For dynamic parameter '"+par.full_name+"' the threshold has become negative.");
						}
						
						if(val > thresh) val = 1;
						else val = 0;
					}
					break;
					
				case BIN_MIN_MAX: case BIN_MIN_MAX_DIST: case BIN_MIN_MAX_REGION: case BIN_MIN_MAX_EQN:
					{
						auto threshmin = eqn[di.threshmin.eq_ref].calculate_all_time(ti,popcomb_t,precalc);
						auto threshmax = eqn[di.threshmax.eq_ref].calculate_all_time(ti,popcomb_t,precalc);
					
						if(threshmin < 0){
							run_error("For dynamic parameter '"+par.full_name+"' the minimum threshold has become negative.");
						}
						
						if(threshmax < 0){
							run_error("For dynamic parameter '"+par.full_name+"' the maximum threshold has become negative.");
						}
						
						if(threshmin > threshmax){
							run_error("For dynamic parameter '"+par.full_name+"' the minimum threshold has become larger than the maximum threshold.");
						}
					
						if(ti == 0){ // At the start use a simple threshold model
							if(val > threshmax) val = 1;
							else val = 0;
						}
						else{ // Otherwise use previous value 
							auto val_last = precalc[i+ti-1];
							if(val_last == 0){
								if(val > threshmax) val = 1;
								else val = 0;
							}
							else{
								if(val_last == 1){
									if(val <= threshmin) val = 0;
									else val = 1;
								}
								else{
									emsg("val_last problem");
								}	
							}
						}
					}
					break;
				}
				
				for(auto tii = ti; tii < ti_end; tii++){
					precalc[i+tii] = val;
				}
			}
		}
	}
	
	//precalc_eqn.print_spec_precalc("sp.txt",upt.spec_precalc);
	//print_spec_precalc("timevary end ",upt.spec_precalc);
	
	precalc_eqn.calculate(upt.spec_precalc,param_val,store);
}


/// Gather a parameter set from a posterior sample
PV Model::post_param(const Sample &samp) const
{
	PV param_val;
	param_val_init(param_val);
	auto &value = param_val.value;
	auto &precalc = param_val.precalc;
	
	/*
	for(auto th = 0u; th < nparam_vec; th++){
		const auto &pv = param_vec[th];	
		const auto &par = param[pv.th];
		cout << par.name << " " << samp.param_value[pv.th][pv.index] << endl;
	}
	*/
	
	for(auto th = 0u; th < nparam_vec; th++){
		const auto &pv = param_vec[th];	
		const auto &par = param[pv.th];
		
		precalc_eqn.calculate(pv.spec_precalc_before,param_val,false);
			
		switch(par.variety){
		case CONST_PARAM:
			value[th] = par.get_value(pv.index);
			break;
			
		case REPARAM_PARAM:
			if(pv.reparam_time_dep == false){
				auto eq_ref = par.get_eq_ref(pv.index);
				if(eq_ref == UNSET) emsg("eq_ref should be set");
			
				value[th] = eqn[eq_ref].calculate_param(precalc);
			}
			break;
			
		case DIST_PARAM: case PRIOR_PARAM:	
			if(pv.ppc_resample){
				value[th] = prior_sample(prior[pv.prior_ref],precalc);
			}
			else{
				value[th] = samp.param_value[pv.th][pv.index];
				if(value[th] == UNSET) emsg("Parameter is unspecified");
			}
			break;
			
		case DYNAMIC_PARAM: emsg("Should not sample dynamic param"); break;
		
		case UNSET_PARAM: emsg("error param"); break;
		}
		
		if(!pv.reparam_time_dep){
			precalc_eqn.calculate(pv.set_param_spec_precalc,param_val,false);
		}		
	}
	
	precalc_eqn.calculate(spec_precalc,param_val,false);
	
	return param_val;
}


/// Gather average over posterior samples
PV Model::post_param_av(const vector <double> &value_av) const
{
	PV param_val;
	param_val_init(param_val);
	auto &value = param_val.value;
	auto &precalc = param_val.precalc;
	
	/*
	for(auto th = 0u; th < nparam_vec; th++){
		const auto &pv = param_vec[th];	
		const auto &par = param[pv.th];
		cout << par.name << " " << samp.param_value[pv.th][pv.index] << endl;
	}
	*/
	
	for(auto th = 0u; th < nparam_vec; th++){
		const auto &pv = param_vec[th];	
		const auto &par = param[pv.th];
		
		precalc_eqn.calculate(pv.spec_precalc_before,param_val,false);
			
		switch(par.variety){
		case CONST_PARAM:
			value[th] = par.get_value(pv.index);
			break;
			
		case REPARAM_PARAM:
			if(pv.reparam_time_dep == false){
				auto eq_ref = par.get_eq_ref(pv.index);
				if(eq_ref == UNSET) emsg("eq_ref should be set");
			
				value[th] = eqn[eq_ref].calculate_param(precalc);
			}
			break;
			
		case DIST_PARAM: case PRIOR_PARAM:	
			value[th] = value_av[th];
			break;
		
		case DYNAMIC_PARAM: emsg("Should not sample dynamic param"); break;
		
		case UNSET_PARAM: emsg("error param"); break;
		}
		
		if(!pv.reparam_time_dep){
			precalc_eqn.calculate(pv.set_param_spec_precalc,param_val,false);
		}
	}
	
	precalc_eqn.calculate(spec_precalc,param_val,false);
	
	return param_val;
}


/// Create a store for ie effects
IEstore Model::post_ie_store(const Sample &samp) const
{
	IEstore ie_store;
	
	const auto &ind_key = ind_key_store[samp.ind_key_ref];
	
	for(auto p = 0u; p < nspecies; p++){
		const auto &sp = species[p];
		
		IEstoreSpecies iess;
		iess.on = true;
		const auto &tab = samp.species[p].ind_tab;
		auto N = tab.nrow;
		auto nie = sp.ind_effect.size();
		
		iess.ie_value.resize(N);
		for(auto i = 0u; i < N; i++){
			iess.ie_value[i].resize(nie,UNSET);
		}
	
		vector <unsigned int> ie_col;
		for(auto i = 0u; i < sp.ind_effect.size(); i++){
			ie_col.push_back(find_in(tab.heading,sp.ind_effect[i].name));
		}
		
		for(auto r = 0u; r < N; r++){
			const auto &row = tab.ele[r];
			
			auto n = number(row[0]);
			if(n == UNSET) emsg("Cannot find ind_key");
				
			auto name = ind_key[n];
			iess.hash.add(r,name);
			
			for(auto i = 0u; i < nie; i++){
				auto c = ie_col[i];
				if(c != UNSET){
					auto n = number(row[c]);
					if(n == UNSET) emsg("IE is not a number");
					iess.ie_value[r][i] = n;
				}
			}
		}
		
		ie_store.species.push_back(iess);
	}
	
	return ie_store;
}


/// Determines if a set of parameters is within the bounds of the prior
bool Model::inbounds(const PV &param_val) const
{
	auto total = prior_total(param_val) + dist_total(param_val);
	
	if(total < -LARGE/2) return false;
	return true;
}


/// The total prior for all parameters
double Model::prior_total(const PV &param_val) const 
{
	auto prior = prior_prob(param_val);
	return sum(prior);
}


/// Calculate the prior for the parameters
vector <double> Model::prior_prob(const PV &param_val) const 
{
	vector <double> pri(nparam_vec,0);
	
	const auto &precalc = param_val.precalc;
	const auto &value = param_val.value;
	
	for(auto th = 0u; th < nparam_vec; th++){
		const auto &pv = param_vec[th];
		if(pv.variety == PRIOR_PARAM){
			pri[th] = prior_probability(value[th],prior[pv.prior_ref],precalc,eqn);
		}
	}
	
	return pri;
}


/// The total prior for all parameters
double Model::dist_total(const PV &param_val) const 
{
	auto dist = dist_prob(param_val);
	return sum(dist);
}


/// Calculate the prior for the parameters
vector <double> Model::dist_prob(const PV &param_val) const
{
	vector <double> dist(nparam_vec,0);
	
	const auto &precalc = param_val.precalc;
	const auto &value = param_val.value;
	
	for(auto th = 0u; th < nparam_vec; th++){
		const auto &pv = param_vec[th];
		if(pv.variety == DIST_PARAM){
			dist[th] = prior_probability(value[th],prior[pv.prior_ref],precalc,eqn);
		}
	}
	
	return dist;
}


/// Recalculates the prior for a given ind effect group
double Model::recalculate_ieg_prior(unsigned int j, vector <double> &prior_ieg, const PV &param_val, double &like_ch) const
{
	auto store = prior_ieg[j];
	prior_ieg[j] = prior_ieg_calculate(ieg_ref[j],param_val);
	like_ch += prior_ieg[j]-store;
	
	return store;
}


/// Recalculates the prior for a given parameter
double Model::recalculate_prior(unsigned int th, vector <double> &prior_prob, const PV &param_val, double &prior_ch, double &prior_bounded_ch) const
{
	const auto &precalc = param_val.precalc;
	const auto &value = param_val.value;
	
	auto store = prior_prob[th];

	const auto &pv = param_vec[th];
	prior_prob[th] = prior_probability(value[th],prior[pv.prior_ref],precalc,eqn);
	
	if(is_prior_bounded(th)) prior_bounded_ch += prior_prob[th]-store;
	else prior_ch += prior_prob[th]-store;
	
	return store;
}


/// Recalculates the distribution for a given parameter
double Model::recalculate_dist(unsigned int th, vector <double> &prior_prob, const PV &param_val, double &like_ch) const
{
	const auto &precalc = param_val.precalc;
	const auto &value = param_val.value;
	
	auto store = prior_prob[th];
	const auto &pv = param_vec[th];
	prior_prob[th] = prior_probability(value[th],prior[pv.prior_ref],precalc,eqn);
	
	like_ch += prior_prob[th]-store;
	
	return store;
}


/// The prior for all the ind effect groups
vector <double> Model::prior_ieg_all(const PV &param_val) const 
{
	vector <double> pri;
	for(auto i = 0u; i < ieg_ref.size(); i++){
		pri.push_back(prior_ieg_calculate(ieg_ref[i],param_val));
	}
	
	return pri;
}


/// The prior for an ind effect group
double Model::prior_ieg_calculate(const IEGref &iegr, const PV &param_val) const
{
	const auto &sp = species[iegr.p];
	
	const auto &ieg = sp.ind_eff_group[iegr.i];
	
	const auto &par = param[ieg.th];
	
	if(par.variety == CONST_PARAM || par.variety == DYNAMIC_PARAM) return UNSET;
	
	auto omega = sp.calculate_omega_basic(iegr.i,param_val,param);

	auto N = ieg.list.size();

	for(auto j = 0u; j < N; j++){
		auto val = omega[j][j];
		if(val <= VAR_MIN || val > VAR_MAX) return -LARGE;
	}	
		
	auto log_det = determinant_fast(omega);
	if(log_det == UNSET) return -LARGE;
	
	const auto &pri = prior[ieg.prior_ref];
	
	switch(pri.type){
	case MVN_DEFAULT_PR:
	case MVN_NORM_LKJ_PR:
	case MVN_UNIFORM_LKJ_PR:
		{
			auto sum = 0.0;
			for(auto j = 0u; j < N; j++){
				auto val = omega[j][j];	
				if(pri.type == MVN_UNIFORM_LKJ_PR){
					if(val < ieg.var_min || omega[j][j] > ieg.var_max) return -LARGE;
				}
				else{
					sum += normal_probability(val,0,ieg.lkj_sd);
				}	
			}

			auto R = sp.calculate_R(iegr.i,param_val,param);
			auto log_det = determinant_fast(R);
			if(log_det == UNSET) return -LARGE;
			sum += log_det*(ieg.lkj_eta-1);
			return sum;
		}
	
	case MVN_INV_WISH_PR:
		{
			auto omega_inv = invert_matrix(omega);
			auto tr = diag_sum(omega_inv);
			return (0.5*N-0.5)*diag_log_sum(omega) - 0.5*(ieg.inv_wish_nu+N+1)*log_det - 0.5*ieg.inv_wish_S*tr;
		}
	
	case MVN_JEF_PR: 
		{
			for(auto i = 0u; i < omega.size(); i++){
				auto val = omega[i][i];
				if(val > ieg.var_max) return -LARGE;
			}
			if(log_det < ieg.log_det_min) return -LARGE;
			return (0.5*N-0.5)*diag_log_sum(omega)-(0.5*N+0.5)*log_det;
		}
		
	case MVN_UNIFORM_PR: 
		for(auto i = 0u; i < omega.size(); i++){
			auto val = omega[i][i];
			if(val < ieg.var_min || val > ieg.var_max) return -LARGE;
		}
		return 0;
		
	default: emsg("opp prob"); break;
	}

	return UNSET;
}


/// Sets ieg_ref which references ind effect groups
void Model::set_ieg_ref()
{
	for(auto p = 0u; p < species.size(); p++){
		const auto &sp = species[p];
		for(auto i = 0u; i < sp.ind_eff_group.size(); i++){
			IEGref iegr; iegr.p = p; iegr.i = i;
			ieg_ref.push_back(iegr);
		}
	}
}


/// Sets param vec which determines covariance matrices
void Model::set_omega_pv()
{
	for(auto &iegr : ieg_ref){
		auto &sp = species[iegr.p];
		auto &ieg = sp.ind_eff_group[iegr.i];
		
		const auto &par = param[ieg.th];
		if(par.variety == CONST_PARAM){ // Checks that set matrix is valid
			auto omega = sp.calculate_omega_const(iegr.i,param);
			auto N = omega.size();
			for(auto j = 0u; j < N; j++){
				for(auto i = 0u; i < N; i++){
					auto val = omega[j][i];
					auto name = exchange_omega(par.name);
					if(ieg.list.size() > 1) name += "_"+ieg.list[j].name+","+ieg.list[i].name; 
					if(j == i){
						if(val <= 0) alert_input("The variance "+name+" must have a positive value");
					}
					else{
						if(val >= 1) alert_input("The correlation "+name+" must be less than 1");
						if(val <= -1) alert_input("The correlation "+name+" must be greater than -1");
					}
				}
			}
			
			sp.convert_cor_var(omega);
			auto log_det = determinant_fast(omega);
			if(log_det == UNSET){
				alert_input("The matrix "+par.name+" must have a positive determinant");
			}						
		}
		else{
			auto N = ieg.list.size();
			ieg.prior_ref = par.element[0].prior_ref;
			
			const auto &pri = prior[ieg.prior_ref];
	
			ieg.log_det_min = UNSET;
			ieg.lkj_eta = ETA_DEFAULT;
			ieg.lkj_sd = SD_DEFAULT;
			ieg.var_min = UNSET;
			ieg.var_max = UNSET;
			ieg.inv_wish_S = UNSET;
			ieg.inv_wish_nu = UNSET;
			
			switch(pri.type){
			case MVN_DEFAULT_PR:
				break;
			
			case MVN_NORM_LKJ_PR:
				ieg.lkj_sd = eqn[pri.dist_param[0].eq_ref].calculate_constant();
				ieg.lkj_eta = eqn[pri.dist_param[1].eq_ref].calculate_constant();
				break;
			
			case MVN_UNIFORM_LKJ_PR:
				ieg.var_min = eqn[pri.dist_param[0].eq_ref].calculate_constant();
				ieg.var_max = eqn[pri.dist_param[1].eq_ref].calculate_constant();
				ieg.lkj_eta = eqn[pri.dist_param[2].eq_ref].calculate_constant();
				break;
	
			case MVN_INV_WISH_PR:
				ieg.inv_wish_S = eqn[pri.dist_param[0].eq_ref].calculate_constant();
				ieg.inv_wish_nu = eqn[pri.dist_param[1].eq_ref].calculate_constant();
				break;

			case MVN_JEF_PR:
				ieg.log_det_min = eqn[pri.dist_param[0].eq_ref].calculate_constant();
				ieg.var_max = eqn[pri.dist_param[1].eq_ref].calculate_constant();
				break;
				
			case MVN_UNIFORM_PR:
				ieg.var_min = eqn[pri.dist_param[0].eq_ref].calculate_constant();
				ieg.var_max = eqn[pri.dist_param[1].eq_ref].calculate_constant();
				break;
				
			default:
				emsg("Prior not regonised");
				break;
			}
				 
			if(ieg.var_min != UNSET){
				if(ieg.var_min <= 0){
					alert_input("For the prior '"+pri.name+"' the minumum value '"+tstr(ieg.var_min)+"' must be positive.");
				}
			}
			
			if(ieg.var_max != UNSET){
				if(ieg.var_max <= 0){
					alert_input("For the prior '"+pri.name+"' the maximum value '"+tstr(ieg.var_max)+"' must be positive.");
				}
			}
			
			if(ieg.var_min != UNSET){
				if(ieg.var_min >= ieg.var_max){
					alert_input("For the prior '"+pri.name+"' the minumum value '"+tstr(ieg.var_min)+"' is larger than the maximum value '"+tstr(ieg.var_max)+"'.");
				}
			}
			
			if(ieg.lkj_sd != UNSET){
				if(ieg.lkj_sd <= 0) alert_input("For the prior '"+pri.name+"' the standard deviation '"+tstr(ieg.lkj_sd)+"' must be positive.");
			}
			
			if(ieg.lkj_eta != UNSET){
				if(ieg.lkj_eta < 1) alert_input("For the prior '"+pri.name+"' eta '"+tstr(ieg.lkj_eta)+"' must not be less than one.");
			}
			
			if(ieg.inv_wish_S != UNSET){
				if(ieg.inv_wish_S <= 0) alert_input("For the prior '"+pri.name+"' the scale S '"+tstr(ieg.inv_wish_S)+"' must be positive.");
			}
			
			if(ieg.inv_wish_nu != UNSET){
				if(ieg.inv_wish_nu <= N+1) alert_input("For the prior '"+pri.name+"' the degrees of freedom '"+tstr(ieg.inv_wish_nu)+"' must be greater than the number of dimentions of the matrix plus one.");
			}
			
			ieg.omega_pv.resize(N);
			for(auto j = 0u; j < N; j++){
				ieg.omega_pv[j].resize(N);
				for(auto i = 0u; i < N; i++){
					auto val = par.get_param_vec(j*N+i);
					if(j > i) val = par.get_param_vec(i*N+j);
					if(val == UNSET) emsg("Problem with omega matrix");
					ieg.omega_pv[j][i] = val;
				}
			}
		}
	}
}


/// Exchanges \Omega with \omega
string Model::exchange_omega(string name) const 
{
	if(begin_str(name,"\\Omega")) return "\\omega"+name.substr(6);
	if(begin_str(name,"Ω")) return "ω"+name.substr(1);
	
	emsg("Omega exhange problem");
	return name;
}


/// Calculates the prior for the spline
vector <double> Model::spline_prior(const PV &param_val) const
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
			
			const auto &tdiv = spl.info.knot_tdiv;
	
			const auto &value = param_val.value;
	
			for(auto t = 1u; t < spl.param_ref.size(); t++){
				const auto &pr1 = spl.param_ref[t]; if(pr1.cons) emsg("SHould not be cons");
				const auto &pr2 = spl.param_ref[t-1]; if(pr2.cons) emsg("SHould not be cons2");
				
				auto val = value[pr1.index];
				auto last = value[pr2.index];
				
				auto dt = tdiv[t]-tdiv[t-1];
				if(dt > 1){  // Only applies spline prior if time difference greater than dt	
					switch(type){
					case LOG_NORMAL_SMOOTH: Li += lognormal_probability(val,last,cv); break;
					case NORMAL_SMOOTH: Li += normal_probability(val,last,sd); break;
					}
				}					
			}
		}
		spline_prior.push_back(Li);
	}
	
	return spline_prior;
}


/// Recalculates the spline prior for a given parameter
double Model::recalculate_spline_prior(unsigned int s, vector <double> &spline_prior, PV &param_val, double &like_ch) const
{
	const auto &spl = spline[s];
	const auto &info = spl.info;
	const auto &tdiv = info.knot_tdiv;
	const auto &value = param_val.value;
	
	auto type = info.smooth_type;
	auto sd = info.smooth_value;
	auto cv = sqrt(exp(sd*sd)-1);
	
	auto Li = 0.0;
	for(auto t = 1u; t < spl.param_ref.size(); t++){
		const auto &pr1 = spl.param_ref[t]; if(pr1.cons) emsg("SHould not be cons");
		const auto &pr2 = spl.param_ref[t-1]; if(pr2.cons) emsg("SHould not be cons2");
				
		auto val = value[pr1.index];
		auto last = value[pr2.index];
		
		auto dt = tdiv[t]-tdiv[t-1];
		if(dt > 1){  // Only applies spline prior if time difference greater than dt	
			switch(type){
			case LOG_NORMAL_SMOOTH: Li += lognormal_probability(val,last,cv); break;
			case NORMAL_SMOOTH: Li += normal_probability(val,last,sd); break;
			}
		}			
	}
	
	auto store = spline_prior[s];
	like_ch += Li-store;
	spline_prior[s] = Li;
	
	return store;
	
}
	

/// Creates species_simp
void Model::create_species_simp()
{
	for(auto &sp :species){
		SpeciesSimp ss(sp.name,sp.cla,sp.ind_effect,sp.fix_effect,sp.comp_mult,sp.comp_gl,sp.tra_gl,sp.trans_tree);
		species_simp.push_back(ss);
	}
}


/// Sets up a distgrid, which groups compartments and provides a minimum distance matrix
void Model::setup_distgrid()
{
	auto N = N_DISTGRID*N_DISTGRID;
	
	for(auto &sp : species){
		for(auto &cla : sp.cla){
			if(cla.comp.size() > 10*N){
				vector < vector <double> > M;
				M.resize(N); for(auto c = 0u; c < N; c++) M[c].resize(N,UNSET);
				
				auto xmin = LARGE, xmax = -LARGE;
				auto ymin = LARGE, ymax = -LARGE;
			
				switch(cla.coord){
				case CARTESIAN: 
					for(auto &co : cla.comp){
						auto x = co.x, y = co.y;
						if(x > xmax) xmax = x; 
						if(x < xmin) xmin = x;
						if(y > ymax) ymax = y; 
						if(y < ymin) ymin = y;
					}
					
					if(xmin == xmax) xmax += TINY;
					if(ymin == ymax) ymax += TINY;
					
					for(auto &co : cla.comp){
						auto x = co.x, y = co.y;
						auto xi = (unsigned int)(ALMOST_ONE*N_DISTGRID*(x-xmin)/(xmax-xmin));
						auto yi = (unsigned int)(ALMOST_ONE*N_DISTGRID*(y-ymin)/(ymax-ymin));
						co.dist_grid_ref = yi*N_DISTGRID + xi;
					}
					
					for(auto yi = 0u; yi < N_DISTGRID; yi++){
						for(auto xi = 0u; xi < N_DISTGRID; xi++){
							auto c = yi*N_DISTGRID + xi;
							for(auto yi2 = yi; yi2 < N_DISTGRID; yi2++){
								for(auto xi2 = 0u; xi2 < N_DISTGRID; xi2++){
									auto cc = yi2*N_DISTGRID + xi2;
									if(cc >= c){
										double xii, yii, xii2, yii2;
										if(xi == xi2){ xii = xi+0.5; xii2 = xii;}
										else{
											if(xi < xi2){ xii = xi+1; xii2 = xi2;}
											else{ xii = xi; xii2 = xi2+1;}
										}
										
										if(yi == yi2){ yii = yi+0.5; yii2 = yii;}
										else{
											if(yi < yi2){ yii = yi+1; yii2 = yi2;}
											else{ yii = yi; yii2 = yi2+1;}
										}
										
										auto dx = (xii2-xii)*(xmax-xmin)/N_DISTGRID;
										auto dy = (yii2-yii)*(ymax-ymin)/N_DISTGRID;
										auto d = sqrt(dx*dx + dy*dy);
										M[c][cc] = d;
										M[cc][c] = d;
									}
								}
							}
						}
					}
					
					if(false){ // Checks that matrix is correct 
						vector < vector <double> > M_ch;
						M_ch.resize(N); for(auto c = 0u; c < N; c++) M_ch[c].resize(N,LARGE);
			
						for(auto c = 0u; c < cla.comp.size(); c++){
							const auto &co = cla.comp[c];
							for(auto cc = 0u; cc < cla.comp.size(); cc++){
								const auto &co2 = cla.comp[cc];
								
								auto x = co.x, y = co.y;
								auto xx = co2.x, yy = co2.y;
								auto dx = x-xx, dy = y-yy;
								auto d = sqrt(dx*dx + dy*dy);
		
								if(d < M_ch[co.dist_grid_ref][co2.dist_grid_ref]) M_ch[co.dist_grid_ref][co2.dist_grid_ref] = d;
							}
						}
						
						for(auto c = 0u; c < N; c++){
							for(auto cc = 0u; cc < N; cc++){
								if(M[c][cc] == UNSET) emsg("Problem unset");
								if(M_ch[c][cc] != LARGE){
									if(M[c][cc] > M_ch[c][cc]) emsg("prob"); 
								}
							}
						}
						
						print_matrix("Mch",M_ch);
					}
					break;
					
				case LATLNG:
					for(auto &co : cla.comp){
						auto x = co.lng, y = co.lat;
						if(x > xmax) xmax = x; 
						if(x < xmin) xmin = x;
						if(y > ymax) ymax = y; 
						if(y < ymin) ymin = y;
					}
					
					if(xmin == xmax) xmax += TINY;
					if(ymin == ymax) ymax += TINY;
					
					for(auto &co : cla.comp){
						auto x = co.lng, y = co.lat;
						auto xi = (unsigned int)(ALMOST_ONE*N_DISTGRID*(x-xmin)/(xmax-xmin));
						auto yi = (unsigned int)(ALMOST_ONE*N_DISTGRID*(y-ymin)/(ymax-ymin));
						auto k = yi*N_DISTGRID + xi;
						if(k >= N) emsg("Out of range");
						co.dist_grid_ref = k;
					}
					
					for(auto yi = 0u; yi < N_DISTGRID; yi++){
						for(auto xi = 0u; xi < N_DISTGRID; xi++){			
							auto c = yi*N_DISTGRID + xi;
							for(auto yi2 = yi; yi2 < N_DISTGRID; yi2++){
								for(auto xi2 = 0u; xi2 < N_DISTGRID; xi2++){	
									auto cc = yi2*N_DISTGRID + xi2;
									if(cc >= c){
										auto d = 0.0, d2 = 0.0;
										if(c != cc){
											if(xi == xi2){
												auto yii = yi+1, yii2 = yi2;
												if(yi2 < yi){ yii = yi, yii2 = yi2+1;}
												d = grid_geo_dist(xi+0.5,yii,xi2+0.5,yii2,xmin,xmax,ymin,ymax);
											}
											else{
												d = grid_geo_dist(xi+1,yi,xi2,yi2,xmin,xmax,ymin,ymax);
												d2 = grid_geo_dist(xi+1,yi+1,xi2,yi2,xmin,xmax,ymin,ymax); if(d2 < d) d = d2; 
												d2 = grid_geo_dist(xi+1,yi,xi2,yi2+1,xmin,xmax,ymin,ymax); if(d2 < d) d = d2; 
												d2 = grid_geo_dist(xi+1,yi+1,xi2,yi2+1,xmin,xmax,ymin,ymax); if(d2 < d) d = d2; 
												d2 = grid_geo_dist(xi,yi,xi2+1,yi2,xmin,xmax,ymin,ymax); if(d2 < d) d = d2; 
												d2 = grid_geo_dist(xi,yi+1,xi2+1,yi2,xmin,xmax,ymin,ymax); if(d2 < d) d = d2; 
												d2 = grid_geo_dist(xi,yi,xi2+1,yi2+1,xmin,xmax,ymin,ymax); if(d2 < d) d = d2; 
												d2 = grid_geo_dist(xi,yi+1,xi2+1,yi2+1,xmin,xmax,ymin,ymax); if(d2 < d) d = d2; 
											}
										}
					
										if(false){
											auto dd = grid_geo_dist(xi,yi,xi2,yi2,xmin,xmax,ymin,ymax);
											auto ty = 0;
											d2 = grid_geo_dist(xi,yi,xi2+1,yi2,xmin,xmax,ymin,ymax); if(d2 < dd){ dd = d2; ty = 1;} 
											d2 = grid_geo_dist(xi,yi,xi2,yi2+1,xmin,xmax,ymin,ymax); if(d2 < dd){ dd = d2; ty = 2;} 
											d2 = grid_geo_dist(xi,yi,xi2+1,yi2+1,xmin,xmax,ymin,ymax); if(d2 < dd){ dd = d2; ty = 3;} 
											d2 = grid_geo_dist(xi+1,yi,xi2,yi2,xmin,xmax,ymin,ymax); if(d2 < dd){ dd = d2; ty = 4;} 
											d2 = grid_geo_dist(xi+1,yi,xi2+1,yi2,xmin,xmax,ymin,ymax); if(d2 < dd){ dd = d2; ty = 5;} 
											d2 = grid_geo_dist(xi+1,yi,xi2,yi2+1,xmin,xmax,ymin,ymax); if(d2 < dd){ dd = d2; ty = 6;} 
											d2 = grid_geo_dist(xi+1,yi,xi2+1,yi2+1,xmin,xmax,ymin,ymax); if(d2 < dd){ dd = d2; ty = 7;} 
											d2 = grid_geo_dist(xi,yi+1,xi2,yi2,xmin,xmax,ymin,ymax); if(d2 < dd){ dd = d2; ty = 8;} 
											d2 = grid_geo_dist(xi,yi+1,xi2+1,yi2,xmin,xmax,ymin,ymax); if(d2 < dd){ dd = d2; ty = 9;} 
											d2 = grid_geo_dist(xi,yi+1,xi2,yi2+1,xmin,xmax,ymin,ymax); if(d2 < dd){ dd = d2; ty = 10;} 
											d2 = grid_geo_dist(xi,yi+1,xi2+1,yi2+1,xmin,xmax,ymin,ymax); if(d2 < dd){ dd = d2; ty = 11;} 
											d2 = grid_geo_dist(xi+1,yi+1,xi2,yi2,xmin,xmax,ymin,ymax); if(d2 < dd){ dd = d2; ty = 12;} 
											d2 = grid_geo_dist(xi+1,yi+1,xi2+1,yi2,xmin,xmax,ymin,ymax); if(d2 < dd){ dd = d2; ty = 13;} 
											d2 = grid_geo_dist(xi+1,yi+1,xi2,yi2+1,xmin,xmax,ymin,ymax); if(d2 < dd){ dd = d2; ty = 14;} 
											d2 = grid_geo_dist(xi+1,yi+1,xi2+1,yi2+1,xmin,xmax,ymin,ymax); if(d2 < dd){ dd = d2; ty = 15;} 
											if(d != dd)	emsg("problem"+tstr(ty));
										}										
										
										M[c][cc] = d;
										M[cc][c] = d;
									}
								}
							}
						}
					}
					
					if(false){ // Checks that matrix is correct 
						vector < vector <double> > M_ch;
						M_ch.resize(N); for(auto c = 0u; c < N; c++) M_ch[c].resize(N,LARGE);
						auto fac = MM_PI/180.0;
			
						for(auto c = 0u; c < cla.comp.size(); c++){
							const auto &co = cla.comp[c];
							for(auto cc = 0u; cc < cla.comp.size(); cc++){
								const auto &co2 = cla.comp[cc];
								auto d = geo_dist(fac*co.lat,fac*co.lng,fac*co2.lat,fac*co2.lng);
								if(d < M_ch[co.dist_grid_ref][co2.dist_grid_ref]) M_ch[co.dist_grid_ref][co2.dist_grid_ref] = d;
							}
						}
						
						for(auto c = 0u; c < N; c++){
							for(auto cc = 0u; cc < N; cc++){
								if(M[c][cc] == UNSET) emsg("Problem unset");
								if(M_ch[c][cc] != LARGE){
									if(M[c][cc] > M_ch[c][cc]) emsg("prob"); 
								}
							}
						}
						
						print_matrix("Mch",M_ch);
					}
					break;
				}
			
				//print_matrix("M",M);
				cla.dist_grid.on = true;
				cla.dist_grid.M = M;
			}
		}			
	}
}


/// Finds the distance between two grid points 
double Model::grid_geo_dist(double xi, double yi, double xi2, double yi2, double xmin, double xmax, double ymin, double ymax)
{		
	auto fr_lat = (ymax-ymin)/N_DISTGRID;
	auto fr_lng = (xmax-xmin)/N_DISTGRID;
	auto fac = MM_PI/180.0;
			
	return geo_dist(fac*(ymin+yi*fr_lat),fac*(xmin+xi*fr_lng),fac*(ymin+yi2*fr_lat),fac*(xmin+xi2*fr_lng));
}


/// Used to order affects
bool AL_ord(const AffectLike &al1, const AffectLike &al2)                      
{ return (al1.order_num < al2.order_num); };  


/// Orders a list of AffectLike based on the priority in the calculation
void Model::order_affect(vector <AffectLike> &vec) const
{
	for(auto &al : vec){
		switch(al.type){
		case SPLINE_PRIOR_AFFECT: case IEG_PRIOR_AFFECT:  case PRIOR_AFFECT: 
			al.order_num = 0; 
			break;
		
		case DIST_AFFECT: case OMEGA_AFFECT: 
		case LIKE_INIT_COND_AFFECT: case PRIOR_INIT_COND_AFFECT:
		case GENETIC_VALUE_AFFECT:
			al.order_num = 1; 
			break;
			
		case EXP_FE_AFFECT: case EXP_IE_AFFECT:
			al.order_num = 2; 
			break;
			
		//case POP_AFFECT:
			//al.order_num = 3; 
			//break;
			
		case INDFAC_INT_AFFECT: 
			al.order_num = 6; 
			break;
		
		case MARKOV_VALUE_AFFECT:
		case POP_DATA_CGL_TGL_AFFECT:
			al.order_num = 7; 
			break;
			
		case OBS_EQN_AFFECT: //case OBS_TRANS_EQN_AFFECT:
		case IIF_W_AFFECT:
		case POPNUM_IND_W_AFFECT:
			al.order_num = 8; 
			break;
			
		case MARKOV_LIKE_AFFECT: case NM_TRANS_AFFECT:  case NM_TRANS_BP_AFFECT: 
		case NM_TRANS_INCOMP_AFFECT: case LIKE_IE_AFFECT: 
		case MARKOV_POP_AFFECT:
		case LIKE_OBS_IND_AFFECT: case LIKE_OBS_POP_AFFECT: 
		case LIKE_OBS_POP_TRANS_AFFECT: case LIKE_UNOBS_TRANS_AFFECT:
		case LIKE_GENETIC_PROCESS_AFFECT: case LIKE_GENETIC_OBS_AFFECT:
			al.order_num = 9; 
			break;
			
		case AFFECT_MAX: break;
		}
	}
	
	sort(vec.begin(),vec.end(),AL_ord);
}

	
/// Works out when to add IIF_W_AFFECT
void Model::add_iif_w_affect(vector <AffectLike> &vec) const
{	
	auto fl = false;

	for(auto i = 0u; i < vec.size(); i++){
		if(vec[i].type == EXP_IE_AFFECT){
			auto p = vec[i].num;
			auto ie = vec[i].num2;
			
			for(auto po : species[p].ind_effect[ie].pop_ref){
				for(const auto &pme : pop[po].markov_eqn_ref){
					const auto &markov_eqn = species[pme.p].markov_eqn;
					for(auto e : pme.list){
						if(markov_eqn[e].infection_trans == true) fl = true;
					}
				}
			}
		}
		
		if(vec[i].type == EXP_FE_AFFECT){
			auto p = vec[i].num;
			auto fe = vec[i].num2;
			
			for(auto po : species[p].fix_effect[fe].pop_ref){
				for(const auto &pme : pop[po].markov_eqn_ref){
					const auto &markov_eqn = species[pme.p].markov_eqn;
					for(auto e : pme.list){
						if(markov_eqn[e].infection_trans == true) fl = true;
					}
				}
			}
		}
	}
	
	if(fl == true){
		AffectLike al; 	
		al.type = IIF_W_AFFECT; al.num = UNSET; al.num2 = UNSET;
		param_vec_add_affect(vec,al);	
	}
}


/// Works out when to add POPNUM_IND_W_AFFECT
void Model::add_popnum_ind_w_affect(vector <AffectLike> &vec) const
{
	for(auto i = 0u; i < vec.size(); i++){
		if(vec[i].type == EXP_IE_AFFECT){
			auto p = vec[i].num;
			const auto &sp = species[p];
			if(sp.trans_tree){
				auto ie = vec[i].num2;
				
				for(auto po : sp.ind_effect[ie].pop_ref){
					AffectLike al; 	
					al.type = POPNUM_IND_W_AFFECT; al.num = po; al.num2 = UNSET;
					param_vec_add_affect(vec,al);	
				}
			}
		}
		
		if(vec[i].type == EXP_FE_AFFECT){
			auto p = vec[i].num;
			const auto &sp = species[p];
			if(sp.trans_tree){
				auto fe = vec[i].num2;
				
				for(auto po : species[p].fix_effect[fe].pop_ref){
					AffectLike al; 	
					al.type = POPNUM_IND_W_AFFECT; al.num = po; al.num2 = UNSET;
					param_vec_add_affect(vec,al);	
				}
			}
		}
	}
}


/// For joint proposals works out the things which need updating
void Model::joint_affect_like(PropType type, const vector <bool> &tr_change, unsigned int p, vector <AffectLike> &vec) const
{
	add_like_obs_affect(p,vec);
	
	const auto &sp = species[p];
	
	auto N = sp.tra_gl.size();
	
	// Gets all the transitions affected by change
	vector <bool> tr_affect(N,false);
	
	switch(type){
	case PAR_EVENT_FORWARD_PROP: case PAR_EVENT_FORWARD_SQ_PROP:
		for(auto tr = 0u; tr < N; tr++){
			if(tr_change[tr]){
				for(auto trf : sp.tr_after[tr]) tr_affect[trf] = true;
			}
		}
		break;
	
	case PAR_EVENT_BACKWARD_SQ_PROP:
		for(auto tr = 0u; tr < N; tr++){
			if(tr_change[tr]){
				for(auto trf : sp.tr_before[tr]) tr_affect[trf] = true;
			}
		}
		break;
	
	default: emsg("Should not be op"); break;
	}
	
	// Works out which and markov equation nmtrans are affected
	auto M = sp.nm_trans.size();
	vector <bool> nm_trans_affect(M,false);
	auto E = sp.markov_eqn.size();
	vector <bool> markov_affect(E,false);
	
	for(auto tr = 0u; tr < N; tr++){
		if(tr_affect[tr]){
			const auto &tra = sp.tra_gl[tr]; 
			if(tra.nm_trans_ref != UNSET) nm_trans_affect[tra.nm_trans_ref] = true;
			
			if(tra.markov_eqn_ref != UNSET) markov_affect[tra.markov_eqn_ref] = true;
		}
	}
	
	auto T = details.T;
	
	// Adds changes to nm_trans
	for(auto m = 0u; m < M; m++){
		if(nm_trans_affect[m]){
			{
				AffectLike al; 
				al.type = NM_TRANS_AFFECT; al.num = p; al.num2 = m;
				al.map.resize(T,true);
				param_vec_add_affect(vec,al);	
			}
			
			auto &nm = sp.nm_trans[m];
			
			if(nm.bp_eq != UNSET){
				AffectLike al; 
				al.type = NM_TRANS_BP_AFFECT; al.num = p; al.num2 = m;
				al.map.resize(T,true);
				param_vec_add_affect(vec,al);
			}	
		
			if(nm.trans_incomp_ref != UNSET){
				AffectLike al; 
				al.type = NM_TRANS_INCOMP_AFFECT; al.num = p; al.num2 = nm.trans_incomp_ref;
				al.map.resize(T,true);
				param_vec_add_affect(vec,al);	
			}
		}
	}
	
	// Adds changes to markov_eqn
	for(auto e = 0u; e < E; e++){
		if(markov_affect[e]){
			const auto &me = sp.markov_eqn[e];
			//const auto &eq = eqn[me.eqn_ref];
			
			AffectLike al; 			
			if(me.time_vari == true) al.map.resize(T,true);
			else al.map.resize(1,true);
				
			al.type = MARKOV_VALUE_AFFECT; 
			al.num = p; al.num2 = e;
			param_vec_add_affect(vec,al);
		
			if(sp.type == INDIVIDUAL){
				al.type = MARKOV_LIKE_AFFECT;
				param_vec_add_affect(vec,al);
			}
		}
	}
	
	// Update individual observations
	{
		auto N = sp.individual.size();
			
		auto fl = false;
		
		vector <bool> map(N,false);
		for(auto i = 0u; i < N; i++){
			const auto &ind = sp.individual[i];
			
			for(auto m = 0u; m < ind.obs.size(); m++){
				const auto &ob = ind.obs[m];
				switch(ob.type){
				case OBS_TRANS_EV: case OBS_SOURCE_EV: case OBS_SINK_EV:
				case OBS_COMP_EV:
					fl = true;
					map[i] = true;	
					break;
				default: break;
				}
			}
		}
		
		if(fl){
			AffectLike al; al.type = LIKE_OBS_IND_AFFECT; al.num = p; al.num2 = UNSET; al.map = map;
			param_vec_add_affect(vec,al);
		}
	}
	
	if(false){
		for(auto tr = 0u; tr < N; tr++){
			cout << sp.tra_gl[tr].name << " " << tr_change[tr] <<" " << tr_affect[tr] << " af " << endl;
		}
	}
//	POP_DATA_CGL_TGL_AFFECT
	
	//AffectLike al; al.type = POP_DATA_CGL_TGL_AFFECT; al.num = p; al.num2 = UNSET;
	//param_vec_add_affect(pp.affect_like,al);
		
}
	

/// Adds all the downstream effects of making a change to a given individual effect
void Model::add_ie_affect(unsigned int p, unsigned int ie, vector <AffectLike> &vec) const
{
	auto T = details.T;
	
	const auto &sp = species[p];
	
	auto &ind_eff = sp.ind_effect[ie];
		
	{
		AffectLike al; al.type = EXP_IE_AFFECT; al.num = p; al.num2 = ie;
		param_vec_add_affect(vec,al);	
	}
	
	if(ind_eff.markov_eqn_ref.size() > 0){ 
		{
			AffectLike al; al.type = INDFAC_INT_AFFECT; al.num = p; al.num2 = UNSET;
			param_vec_add_affect(vec,al);		
		}
		
		for(auto e : ind_eff.markov_eqn_ref){	
			if(sp.type == INDIVIDUAL){
				AffectLike al; al.type = MARKOV_LIKE_AFFECT; al.num = p; al.num2 = e;
				auto &me = sp.markov_eqn[e];
				if(me.time_vari == true) al.map.resize(T,true);
				else al.map.resize(1,true);

				param_vec_add_affect(vec,al);	
			}
		}
	}
	
	for(auto e : ind_eff.nm_trans_ref){	
		{
			AffectLike al; 
			al.type = NM_TRANS_AFFECT; al.num = p; al.num2 = e;
			al.map.resize(T,true);
			param_vec_add_affect(vec,al);	
		}
		
		auto &nm = sp.nm_trans[e];
		
		if(nm.bp_eq != UNSET){
			AffectLike al; 
			al.type = NM_TRANS_BP_AFFECT; al.num = p; al.num2 = e;
			al.map.resize(T,true);
			param_vec_add_affect(vec,al);
		}	
	
		if(nm.trans_incomp_ref != UNSET){
			AffectLike al; 
			al.type = NM_TRANS_INCOMP_AFFECT; al.num = p; al.num2 = nm.trans_incomp_ref;
			al.map.resize(T,true);
			param_vec_add_affect(vec,al);	
		}
	}
	
	for(auto po : ind_eff.pop_ref){	
		add_pop_affect(po,vec);
	}	
}


/// Incorporate the effect of altering population po for parameter k
void Model::add_pop_affect(unsigned int po, vector <AffectLike> &vec) const
{
	const auto &popu = pop[po];
	auto T = details.T;

	/*
	AffectLike al; 
	al.type = POP_AFFECT; al.num = UNSET; al.num2 = UNSET; 
	//al.map.resize(pop.size(),false); al.map[po] = true;
	al.map.resize(npop,false); al.map[po] = true;
	param_vec_add_affect(vec,al);		
	*/
	
	for(const auto &mer : popu.markov_eqn_ref){
		auto p = mer.p;
		const auto &sp = species[p];
		if(sp.type != DETERMINISTIC){
			for(auto e : mer.list){
				AffectLike al; al.type = MARKOV_VALUE_AFFECT;
				al.num = p; al.num2 = e; 
		
				auto &me = sp.markov_eqn[e];
				if(me.time_vari == true) al.map.resize(T,true);
				else al.map.resize(1,true);
				
				param_vec_add_affect(vec,al);		
				
				if(sp.type == INDIVIDUAL){
					al.type = MARKOV_LIKE_AFFECT; 
					param_vec_add_affect(vec,al);
				}		
			}			
		}
	}
	
	for(const auto &tref : popu.trans_ref){
		auto p = tref.p;
		for(auto tr : tref.tr_list){
			AffectLike al; al.map.resize(T,true);
			al.type = MARKOV_POP_AFFECT; al.num = p; al.num2 = tr;
			param_vec_add_affect(vec,al);
		}
	}
}

	
/// Samples the initial conditions for the state
vector <InitCondValue> Model::initc_sample(const PV &param_val) const
{
	vector <InitCondValue> initc(nspecies);
	
	const auto &precalc = param_val.precalc;
	
	for(auto p = 0u; p < nspecies; p++){
		const auto &sp = species[p];
		const auto &ic = sp.init_cond;
		auto &inc = initc[p];
		
		inc.N_total_unobs = UNSET;
		
		auto N = sp.comp_gl.size();
		
		switch(ic.type){
		case INIT_POP_FIXED:
			inc.cnum = ic.cnum;
			break;
		
		case INIT_POP_DIST:
			{
				auto foc_cl = ic.focal_cl;
				if(foc_cl == UNSET){
					inc.N_total = round_int(prior_sample(ic.pop_prior,precalc));
					inc.frac = dirichlet_sample(ic.alpha);
					inc.cnum = multinomial_sample(inc.N_total,inc.frac);
				}
				else{
					const auto &claa = sp.cla[foc_cl];
					inc.N_focal.resize(claa.ncomp);
					for(auto c = 0u; c < claa.ncomp; c++){
						inc.N_focal[c] = round_int(prior_sample(ic.comp_prior[c],precalc));
					}
					
					inc.frac_focal.resize(sp.ncla);
					for(auto cl = 0u; cl < sp.ncla; cl++){
						if(cl != foc_cl){
							inc.frac_focal[cl] = dirichlet_sample(ic.alpha_focal[cl]);
						}
					}
					
					inc.frac_comb = frac_combine(p,inc.frac_focal);
			
					inc.cnum_reduce.resize(claa.ncomp);
					for(auto c = 0u; c < claa.ncomp; c++){
						inc.cnum_reduce[c] = multinomial_sample(inc.N_focal[c],inc.frac_comb);
					}
					
					combine_cnum_reduce(p,inc);
				}
			}
			break;
	
		case INIT_POP_NONE: break;
	
		default: emsg("Should not be here init pop"); break;
		}
	
		if(false){
			for(auto cgl = 0u; cgl < N; cgl++){
				cout << sp.comp_gl[cgl].name << " " << inc.cnum[cgl] << endl;
			}	
			emsg("Initial population");
		}
	}
	
	return initc;
}


/// Combines fractions together
vector <double> Model::frac_combine(unsigned int p, const vector < vector <double> > &frac) const
{
	const auto &sp = species[p];
	const auto &ic = sp.init_cond;
	
	auto foc_cl = ic.focal_cl;
	
	vector <double> frac_comb;
	for(auto cr = 0u; cr < ic.N_reduce; cr++){
		auto &co = sp.comp_gl[ic.comp_reduce[cr]];
		
		auto fac = 1.0;
		for(auto cl = 0u; cl < sp.ncla; cl++){
			if(cl != foc_cl){
				fac *= frac[cl][co.cla_comp[cl]];
			}
		}
		
		frac_comb.push_back(fac);
	}
	
	return frac_comb;
}


/// Combines together the reduced initial populations to get cnum
void Model::combine_cnum_reduce(unsigned p, InitCondValue &inc) const
{
	const auto &sp = species[p];
	const auto &ic = sp.init_cond;
	auto foc_cl = ic.focal_cl;
	
	auto &cnum = inc.cnum;
	cnum.resize(sp.comp_gl.size());
	for(auto c = 0u; c < sp.cla[foc_cl].ncomp; c++){
		auto shift = c*sp.comp_mult[foc_cl];
		for(auto j = 0u; j < ic.N_reduce; j++){
			cnum[ic.comp_reduce[j]+shift] = inc.cnum_reduce[c][j];
		}
	}
}


/// Samples from the prior
double Model::prior_sample(const Prior &pri, const vector <double> &precalc) const
{				
	double val=UNSET;
	string warn;
	
	switch(pri.type){
	case MVN_DEFAULT_PR: case MVN_NORM_LKJ_PR: case MVN_UNIFORM_LKJ_PR: case MVN_INV_WISH_PR:
	case MVN_JEF_PR: case MVN_UNIFORM_PR: case MVN_COR_PR:  // Sampling for this is done collectively
		return UNSET;
	
	case INVERSE_PR:
		{	
			auto min = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			auto max = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
			if(min <= 0) run_error("For the prior '"+pri.name+"' the minumum value '"+tstr(min)+"' has become non-positive.");
			
			if(max <= 0) run_error("For the prior '"+pri.name+"' the maximum value '"+tstr(max)+"' has become non-positive.");
			
			if(min >= max) run_error("For the prior '"+pri.name+"' the minumum value '"+tstr(min)+"' is larger than the maximum value '"+tstr(max)+"'.");
			
			val = inverse_sample(min,max,warn);
		}
		break;
		
	case POWER_PR:
		{	
			auto min = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			auto max = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
			if(min <= 0) run_error("For the prior '"+pri.name+"' the minumum value '"+tstr(min)+"' has become non-positive.");
			
			if(max <= 0) run_error("For the prior '"+pri.name+"' the maximum value '"+tstr(max)+"' has become non-positive.");
			
			if(min >= max) run_error("For the prior '"+pri.name+"' the minumum value '"+tstr(min)+"' is larger than the maximum value '"+tstr(max)+"'.");
			
			auto power = eqn[pri.dist_param[2].eq_ref].calculate_param(precalc);
			
			val = power_sample(min,max,power,warn);
		}
		break;
		
	case UNIFORM_PR:
		{	
			auto min = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			auto max = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
		
			if(min >= max) run_error("For the prior '"+pri.name+"' the minumum value '"+tstr(min)+"' is larger than the maximum value '"+tstr(max)+"'.");
			
			val = uniform_sample(min,max,warn);
		}
		break;
		
	case EXP_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			val = exp_mean_sample(mean,warn);
		}
		break;
	
	case NORMAL_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			auto sd = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
			val = normal_sample(mean,sd,warn);
		}
		break;
		
	case GAMMA_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			auto cv = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
			val = gamma_sample(mean,cv,warn);
		}
		break;
	
	case LOG_NORMAL_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			auto cv = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
			val = lognormal_sample(mean,cv,warn);
		}
		break;
		
	case BETA_PR:
		{	
			auto alpha = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			auto beta = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
			val = beta_sample(alpha,beta,warn);
		}
		break;
	
	case BERNOULLI_PR:
		{	
			auto z = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			val = bernoulli_sample(z,warn);
		}
		break;
		
	case FIX_PR:
		{
			val = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
		}
		break;
		
	case DIRICHLET_PR:
		{
			auto alpha = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			val = gamma_alpha_sample(alpha,warn);
		}
		break;
		
	case MDIR_PR:	
		val = UNSET;
		emsg("Should not sample from MDIR");
		break;
	
	case UNSET_PR:
		val = UNSET;
		emsg("Prior should be set");
		break;
	}
	
	if(val == UNSET) prior_error(pri,warn);
	
	return val;
}


/// Outputs an error message for the prior
void Model::prior_error(const Prior &pri, string warn) const 
{
	run_error("For prior '"+pri.name+"' the following error occurred: "+warn+". This prior is for '"+pri.in+"'. Prior distribution quantities (means, sds etc...) have threshold limits to ensure numerical accuracy. Consider changing the values or restricting the priors on model parameter which determine this distribution"); 
}


/// Samples from N divisions of the prior using Latin hypercube sampling 
vector <double> Model::latin_prior_sample(const Prior &pri, unsigned int N) const
{				

	vector <double> val_vec;
	
	string warn;
		
	switch(pri.type){
	case MVN_DEFAULT_PR: case MVN_NORM_LKJ_PR: case MVN_UNIFORM_LKJ_PR: case MVN_INV_WISH_PR:
	case MVN_JEF_PR: case MVN_UNIFORM_PR: case MVN_COR_PR: 
		emsg("Should not be here9");
		break;
	
	case INVERSE_PR:
		{	
			auto min = eqn[pri.dist_param[0].eq_ref].calculate_value();
			auto max = eqn[pri.dist_param[1].eq_ref].calculate_value();
			if(min <= 0) run_error("For the prior '"+pri.name+"' the minumum value '"+tstr(min)+"' has become non-positive.");
			
			if(max <= 0) run_error("For the prior '"+pri.name+"' the maximum value '"+tstr(max)+"' has become non-positive.");
			
			if(min >= max) run_error("For the prior '"+pri.name+"' the minumum value '"+tstr(min)+"' is larger than the maximum value '"+tstr(max)+"'.");
			
			for(auto i = 0u; i < N; i++){
				auto val = inverse_from_cdf((i+ran())/N,min,max,warn);
				if(val == UNSET) prior_error(pri,warn);
				val_vec.push_back(val);
			}
		}
		break;
		
	case POWER_PR:
		{	
			auto min = eqn[pri.dist_param[0].eq_ref].calculate_value();
			auto max = eqn[pri.dist_param[1].eq_ref].calculate_value();
			if(min <= 0) run_error("For the prior '"+pri.name+"' the minumum value '"+tstr(min)+"' has become non-positive.");
			
			if(max <= 0) run_error("For the prior '"+pri.name+"' the maximum value '"+tstr(max)+"' has become non-positive.");
			
			if(min >= max) run_error("For the prior '"+pri.name+"' the minumum value '"+tstr(min)+"' is larger than the maximum value '"+tstr(max)+"'.");
			
			auto power = eqn[pri.dist_param[2].eq_ref].calculate_value();
			
			for(auto i = 0u; i < N; i++){
				auto val = power_from_cdf((i+ran())/N,min,max,power,warn);
				if(val == UNSET) prior_error(pri,warn);
				val_vec.push_back(val);
			}
		}
		break;
		
	case UNIFORM_PR:
		{	
			auto min = eqn[pri.dist_param[0].eq_ref].calculate_value();
			auto max = eqn[pri.dist_param[1].eq_ref].calculate_value();
			if(min >= max) run_error("For the prior '"+pri.name+"' the minumum value '"+tstr(min)+"' is larger than the maximum value '"+tstr(max)+"'.");
			
			for(auto i = 0u; i < N; i++){
				auto val = uniform_from_cdf((i+ran())/N,min,max,warn);
				if(val == UNSET) prior_error(pri,warn);
				val_vec.push_back(val);
			}
		}
		break;
		
	case EXP_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_value();
			
			for(auto i = 0u; i < N; i++){
				auto val = exp_mean_from_cdf((i+ran())/N,mean,warn);
				if(val == UNSET) prior_error(pri,warn);
				val_vec.push_back(val);
			}
		}
		break;
	
	case NORMAL_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_value();
			auto sd = eqn[pri.dist_param[1].eq_ref].calculate_value();
			
			for(auto i = 0u; i < N; i++){
				auto val = normal_from_cdf((i+ran())/N,mean,sd,warn);
				if(val == UNSET) prior_error(pri,warn);
				val_vec.push_back(val);
			}
		}
		break;
		
	case GAMMA_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_value();
			auto cv = eqn[pri.dist_param[1].eq_ref].calculate_value();
		
			for(auto i = 0u; i < N; i++){
				auto val = gamma_from_cdf((i+ran())/N,mean,cv,warn);
				if(val == UNSET) prior_error(pri,warn);
				val_vec.push_back(val);
			}
		}
		break;
	
	case LOG_NORMAL_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_value();
			auto cv = eqn[pri.dist_param[1].eq_ref].calculate_value();
			
			for(auto i = 0u; i < N; i++){
				auto val = lognormal_from_cdf((i+ran())/N,mean,cv,warn);
				if(val == UNSET) prior_error(pri,warn);
				val_vec.push_back(val);
			}
		}
		break;
		
	case BETA_PR:
		{	
			auto alpha = eqn[pri.dist_param[0].eq_ref].calculate_value();
			auto beta = eqn[pri.dist_param[1].eq_ref].calculate_value();
		
			for(auto i = 0u; i < N; i++){
				auto val = beta_from_cdf((i+ran())/N,alpha,beta,warn);
				if(val == UNSET) prior_error(pri,warn);
				val_vec.push_back(val);
			}
		}
		break;
	
	case BERNOULLI_PR:
		{	
			auto z = eqn[pri.dist_param[0].eq_ref].calculate_value();
			
			for(auto i = 0u; i < N; i++){
				auto val = bernoulli_from_cdf((i+ran())/N,z,warn);
				if(val == UNSET) prior_error(pri,warn);
				val_vec.push_back(val);
			}
		}
		break;
		
	case FIX_PR:
		{
			auto val = eqn[pri.dist_param[0].eq_ref].calculate_value();
			for(auto i = 0u; i < N; i++){
				val_vec.push_back(val);
			}
		}
		break;
		
	case DIRICHLET_PR:
		{
			auto alpha = eqn[pri.dist_param[0].eq_ref].calculate_value();
			
			for(auto i = 0u; i < N; i++){
				auto val = gamma_alpha_from_cdf((i+ran())/N,alpha,warn);
				if(val == UNSET) prior_error(pri,warn);
				val_vec.push_back(val);
			}
		}
		break;
		
	case MDIR_PR:	
		emsg("Should not sample from MDIR");
		break;
	
	case UNSET_PR:
		emsg("Prior should be set");
		break;
	}
	
	return val_vec;
}


/// In pop_data or trans_pop_data then this need to be recalculated under MBP
void Model::add_like_obs_affect(unsigned int p, vector <AffectLike> &affect_like) const
{
	const auto &sp = species[p];
	auto npop_data = sp.pop_data.size();
	auto npop_trans_data = sp.pop_trans_data.size();

	if(npop_data > 0 || npop_trans_data > 0){
		AffectLike al; al.type = POP_DATA_CGL_TGL_AFFECT; al.num = p; al.num2 = UNSET;
		param_vec_add_affect(affect_like,al);
	}
	
	if(npop_data > 0){
		AffectLike al; al.type = LIKE_OBS_POP_AFFECT; al.num = p; al.num2 = UNSET; al.map = true_vec(npop_data);
		param_vec_add_affect(affect_like,al);
	}
	
	if(npop_trans_data > 0){		
		AffectLike al; al.type = LIKE_OBS_POP_TRANS_AFFECT; al.num = p; al.num2 = UNSET; al.map = true_vec(npop_trans_data);
		param_vec_add_affect(affect_like,al);
	}
}


/// Finds an individual reference based on its name
AllInd Model::find_all_ind(string name) const
{
	auto j = hash_all_ind.find(name);
	if(j == UNSET) emsg("Could not find individual '"+name+"'");
	return all_ind[j];
}


/// Sets a hash table for all individuals
void Model::set_hash_all_ind()
{
	for(auto p = 0u; p < species.size(); p++){
		const auto &sp = species[p];
		if(sp.type == INDIVIDUAL){
			for(auto i = 0u; i < sp.individual.size(); i++){
				hash_all_ind.add(all_ind.size(),sp.individual[i].name);
				AllInd ai; ai.i = i; ai.p = p;
				all_ind.push_back(ai);
			}
		}
	}
}


/// Checks that the prior for the ind effect group is within bounds
bool Model::ieg_check_prior_error(const IEGref &iegr, const PV &param_val) const
{
	const auto &sp = species[iegr.p];

	const auto &ieg = sp.ind_eff_group[iegr.i];
	
	const auto &par = param[ieg.th];
	
	if(par.variety == CONST_PARAM || par.variety == DYNAMIC_PARAM) return false;
	
	auto omega = sp.calculate_omega_basic(iegr.i,param_val,param);
	
	auto N = omega.size();
	
	for(auto j = 0u; j < N; j++){
		auto val = omega[j][j];
		if(val <= VAR_MIN || val > VAR_MAX) return true;
	}	
		
	auto log_det = determinant_fast(omega);
	if(log_det == UNSET) return true;
	
	const auto &pri = prior[ieg.prior_ref];
	
	switch(pri.type){
	case MVN_DEFAULT_PR:
	case MVN_NORM_LKJ_PR:
	case MVN_UNIFORM_LKJ_PR:
		{
			for(auto j = 0u; j < N; j++){
				auto val = omega[j][j];
				if(pri.type == MVN_UNIFORM_LKJ_PR){
					if(val < ieg.var_min || omega[j][j] > ieg.var_max) return true;
				}
			}

			auto R = sp.calculate_R(iegr.i,param_val,param);
			auto log_det = determinant_fast(R);
			if(log_det == UNSET) return true;
		}
		break;
	
	case MVN_INV_WISH_PR:
		break;
		
	case MVN_JEF_PR: 
		{
			for(auto i = 0u; i < omega.size(); i++){
				auto val = omega[i][i];
				if(val > ieg.var_max) return true;
			}
			if(log_det < ieg.log_det_min) return true;
		}
		break;
		
	case MVN_UNIFORM_PR: 
		for(auto i = 0u; i < omega.size(); i++){
			auto val = omega[i][i];
			if(val < ieg.var_min || val > ieg.var_max) return true;
		}
		break;
		
	default: emsg("opp prob"); break;
	}
	
	return false;
}


/// Prints a set of parameters
void Model::print_param(const PV &param_val) const
{
	if(com_op == true) return;
	
	const auto &value = param_val.value;
	
	auto max = 100u;
	
	cout << "Parameters:" << endl;
	auto thmax = param_vec.size();
	if(thmax > max) thmax = max;
	for(auto th = 0u; th < thmax; th++){
		cout << param_vec_name(th) << " "<< value[th] << endl;
	}
	if(thmax == max) cout << "..." << endl;
	cout << endl;
	
	const auto &precalc = param_val.precalc;
	
	cout << "Precalc:" << endl;
	auto premax = precalc.size();
	if(premax > max) premax = max;
	for(auto i = 0u; i < premax; i++){
		cout << precalc[i] << ",";
	}
	if(premax == max) cout << "..." << endl;
	cout << endl;
}


/// Gets param_val_prop from a full parameter vector (i.e. removes reparam)
vector <double> Model::get_param_val_prop(const PV &param_val) const
{
	const auto &value = param_val.value;
	
	vector <double> vec(nparam_vec_prop);
	for(auto i = 0u; i < nparam_vec_prop; i++){
		vec[i] = value[param_vec_prop[i]];
	}
	
	return vec;
}


/// Gets tvreparam from a full parameter vector (i.e. removes other parameters)
vector <double> Model::get_param_val_tvreparam(const PV &param_val) const
{
	const auto &value = param_val.value;
	
	vector <double> vec(nparam_vec_tvreparam);
	for(auto i = 0u; i < nparam_vec_tvreparam; i++){
		vec[i] = value[param_vec_tvreparam[i]];
	}
	
	return vec;
}


/// Stores the values of dynamic parameters
vector < vector <double> > Model::get_param_val_dynamic(const PV &param_val) const
{
	vector < vector <double> > value;
	for(auto th : dynamic_param){
		vector <double> val;
		
		const auto &par = param[th];
		if(par.trace_output || par.state_output){
			auto s_start = par.spline_ref;
			
			vector <unsigned int> ti_list;
			for(const auto &seg : par.spline_segment) ti_list.push_back(seg.ti_min);
			
			auto ntimes = par.spline_info.knot_tdiv.size();
			if(ntimes != ti_list.size()) emsg("Problem with list size");
		
			auto N = par.N/ntimes;
			
			auto &precalc = param_val.precalc;
			for(auto j = 0u; j < N; j++){
				auto i = spline_refi[s_start+j];
				for(auto k = 0u; k < ntimes; k++){
					val.push_back(precalc[i+ti_list[k]]);
				}
			}
		}
		value.push_back(compress_vec(val));
	}

	return value;
}


/// Compresses a vector by combining together the same value
vector <double> Model::compress_vec(const vector <double> &vec) const
{
	vector <double> vc;
	
	auto si = vec.size();
	
	auto i = 0u; 
	while(i < si){
		auto val = vec[i];
		auto n = 0u;
		while(i+1 < si && vec[i+1] == val){ i++; n++;}
		vc.push_back(LARGE*n+val);
		i++;
	}

	return vc;
}


/// Decompresses a vector
vector <double> Model::decompress_vec(const vector <double> &vc) const
{ 
	vector <double> vec;
	
	for(auto &va : vc){
		unsigned int n = va/LARGE;
		auto val = va - n*LARGE;
		for(auto i = 0u; i <= n; i++) vec.push_back(val);
	}

	return vec;
}


/// For the output of parameters tvreparam are placed into param_val
void Model::add_tvreparam(PV &param_val, const vector <double> &param_val_tvreparam) const 
{
	if(param_val_tvreparam.size() == 0) return;
	
	auto &value = param_val.value;
	for(auto i = 0u; i < nparam_vec_tvreparam; i++){
		value[param_vec_tvreparam[i]] = param_val_tvreparam[i];
	}
}


/// Reconstructs param_val from param_val_prop
PV Model::get_param_val(const Particle &pa) const
{
	PV param_val;
	param_val_init(param_val);
	auto &value = param_val.value;
	auto &precalc = param_val.precalc;

	for(auto &va : value) va = UNSET;
	
	for(auto i = 0u; i < nparam_vec_prop; i++){
		auto th = param_vec_prop[i];
	
		value[th] = pa.param_val_prop[i];
	
		precalc_eqn.calculate(param_vec[th].set_param_spec_precalc,param_val,false);
	}

	for(auto th = 0u; th < nparam_vec; th++){
		const auto &pv = param_vec[th];
		
		precalc_eqn.calculate(pv.spec_precalc_before,param_val,false);
	
		if(value[th] == UNSET){
			const auto &par = param[pv.th];
				
			switch(par.variety){
			case PRIOR_PARAM: case DIST_PARAM:
				{
					const auto &pri = prior[par.get_prior_ref(pv.index)];
					if(pri.type != FIX_PR) emsg("Prior should be fixed");
					
					value[th] = prior_sample(pri,precalc);
					precalc_eqn.calculate(pv.set_param_spec_precalc,param_val,false);
				}
				break;
				
			case REPARAM_PARAM:
				if(pv.reparam_time_dep == false){
					auto eq_ref = par.get_eq_ref(pv.index);
					if(eq_ref == UNSET) emsg("eq_ref should be set");
					
					value[th] = eqn[eq_ref].calculate_param(precalc);
					precalc_eqn.calculate(pv.set_param_spec_precalc,param_val,false);
				}
				break;
			
			default: emsg("Option prob"); break;
			}
		}
	}
	
	precalc_eqn.calculate(spec_precalc,param_val,false);

	if(false){
		for(auto val: param_val.value) cout << val << ","; 	
		cout << " value" << endl;
		for(auto val: param_val.precalc) cout << val << ","; 	
		cout << " precalc" << endl;
		emsg("get param val");
	}
	
	return param_val;
}


/// Sets PV from a set of parameter valuesparam_val
PV Model::set_param_val(const vector <double> &value_set) const
{
	PV param_val;
	param_val_init(param_val);
	auto &value = param_val.value;
	auto &precalc = param_val.precalc;

	for(auto &va : value) va = UNSET;
	
	for(auto i = 0u; i < nparam_vec_prop; i++){
		auto th = param_vec_prop[i];
	
		value[th] = value_set[th];
	
		precalc_eqn.calculate(param_vec[th].set_param_spec_precalc,param_val,false);
	}

	for(auto th = 0u; th < nparam_vec; th++){
		const auto &pv = param_vec[th];
		
		precalc_eqn.calculate(pv.spec_precalc_before,param_val,false);
	
		if(value[th] == UNSET){
			const auto &par = param[pv.th];
				
			switch(par.variety){
			case PRIOR_PARAM: case DIST_PARAM:
				{
					const auto &pri = prior[par.get_prior_ref(pv.index)];
					if(pri.type != FIX_PR) emsg("Prior should be fixed");
					
					value[th] = prior_sample(pri,precalc);
					precalc_eqn.calculate(pv.set_param_spec_precalc,param_val,false);
				}
				break;
				
			case REPARAM_PARAM:
				if(pv.reparam_time_dep == false){
					auto eq_ref = par.get_eq_ref(pv.index);
					if(eq_ref == UNSET) emsg("eq_ref should be set");
					
					value[th] = eqn[eq_ref].calculate_param(precalc);
					precalc_eqn.calculate(pv.set_param_spec_precalc,param_val,false);
				}
				break;
			
			default: emsg("Option prob"); break;
			}
		}
	}
	
	precalc_eqn.calculate(spec_precalc,param_val,false);

	if(false){
		for(auto val: param_val.value) cout << val << ","; 	
		cout << " value" << endl;
		for(auto val: param_val.precalc) cout << val << ","; 	
		cout << " precalc" << endl;
		emsg("get param val");
	}
	
	return param_val;
}


Param::Param(Constant &constant) : constant(constant) 
{
}

/// Gets the value for param vec
unsigned int Param::get_param_vec(unsigned int i) const
{
	const auto &er = element_ref[i];
	auto ind = er.index; 
	if(ind == UNSET) emsg("Cannot get element ref2");
	if(er.cons != false) emsg("Should not get constant1: "+name);
	return element[ind].param_vec_ref;
}


/// Gets a vector of children from an element
const vector <ParamRef>& Param::get_child(unsigned int i) const
{
	const auto &er = element_ref[i];
	auto ind = er.index; 
	if(ind == UNSET) emsg("Cannot get element ref2");
	if(er.cons != false) emsg("Should not get constant2");
	return element[ind].child;
}


/// Gets a vector of parents from an element
const vector <ParamRef>& Param::get_parent(unsigned int i) const
{
	const auto &er = element_ref[i];
	auto ind = er.index; 
	if(ind == UNSET) emsg("Cannot get element ref");
	if(er.cons != false) emsg("Should not get constant3");
	return element[ind].parent;
}


/// Gets equation reference
unsigned int Param::get_eq_ref(unsigned int i) const
{
	const auto &er = element_ref[i];
	auto ind = er.index; 
	if(ind == UNSET) emsg_input("Cannot get element ref3");
	if(er.cons != false) emsg("Should not get constant4");
	return element[ind].value.eq_ref;
}


/// Adds all elements
void Param::all_elements()
{
	for(auto i = 0u; i < N; i++) add_element(i);
}


/// Adds an element to the parameter
void Param::add_element(unsigned int i, bool all_spline) 
{
	auto &er = element_ref[i];
	if(er.index == UNSET){
		er.index = element.size();
		er.cons = false;
		
		ParamElement ele;
		ele.prior_ref = default_prior_ref;
		ele.used = false;
		ele.value.eq_ref = UNSET;
		element.push_back(ele);
	}
	
	// If a reparameterised spline then also add spline points
	if(all_spline && variety == REPARAM_PARAM && time_dep){
		const auto &dep_la = dep[dep.size()-1];
		auto N = dep_la.list.size();
	
		if(i%N == 0){
			for(auto j = 1u; j < N; j++){
				auto ii = i+j;
				auto &er2 = element_ref[ii];
				if(er2.index == UNSET){
					er2.index = element.size();
					er2.cons = false;
					
					ParamElement ele;
					ele.prior_ref = default_prior_ref;
					ele.used = false;
					element.push_back(ele);
				}
			}
		}
	}
}


/// Sets the prior
void Param::set_prior(unsigned int i, unsigned int prior_ref)
{
	add_element(i);
	element[element_ref[i].index].prior_ref = prior_ref;
}


/// Adds a parent onto the list
void Param::add_parent(unsigned int i, const ParamRef &pr)
{
	const auto &er = element_ref[i];
	auto ind = er.index; if(ind == UNSET) emsg_input("Cannot get element ref4");
	if(er.cons != false) emsg("Should not get constant5");
	return element[ind].parent.push_back(pr);
}


/// Adds a child onto the list
void Param::add_child(unsigned int i, const ParamRef &pr)
{
	const auto &er = element_ref[i];
	auto ind = er.index; if(ind == UNSET) emsg_input("Cannot get element ref");
	if(er.cons != false) emsg("Should not get constant6");
	return element[ind].child.push_back(pr);
}


/// Gets the value of an element
double Param::get_value(unsigned int i) const 
{
	const auto &er = element_ref[i];
	auto ind = er.index; if(ind == UNSET) emsg("Cannot get element ref5");
	if(er.cons) return constant.value[ind];
	return element[ind].value.value;
}


/// Gets the value of an element
unsigned int Param::get_prior_ref(unsigned int i) const 
{
	const auto &er = element_ref[i];
	auto ind = element_ref[i].index; if(ind == UNSET) emsg_input("Cannot get element ref6");
	if(er.cons != false) emsg("Should not get constant7");
	return element[ind].prior_ref;
}



/// Sets an element as used
void Param::set_used(unsigned int i) 
{
	const auto &er = element_ref[i];
	if(!er.cons && variety != DYNAMIC_PARAM){
		auto ind = er.index; if(ind == UNSET) emsg_input("Cannot get element ref7");
		element[ind].used = true;
	}
}


/// Gets the value of an element
string Param::get_value_te(unsigned int i) const
{
	const auto &er = element_ref[i];
	auto ind = er.index; if(ind == UNSET) emsg_input("Cannot get element ref8");
	if(er.cons) emsg_input("Should not be cons");
	return element[ind].value.te;
}


/// Sets the value of an element
void Param::set_value_eqn(unsigned int i, const EquationInfo &val) 
{
	add_element(i);
	element[element_ref[i].index].value = val;
}


/// Gets the value of an element
bool Param::exist(unsigned int i) const
{
	if(element_ref[i].index == UNSET) return false;
	return true;
}


/// Adds a constant value to the list
unsigned int Param::add_cons(double val)
{					
	return constant.add(val); 
	/* turn off
	auto num = cons.size();
	if(num < 10){
		auto i = 0u; while(i < num && cons[i] != val) i++;
		if(i < num) return i;
	}
	cons.push_back(val);
	return num;
	*/
}


/// Adds a constant value to the list
void Param::set_cons(unsigned int i, double val)
{					
	element_ref[i].index = add_cons(val);
	element_ref[i].cons = true;
}


/// Calculates tdiv from t
double Model::calc_tdiv(double t) const 
{
	return (t - details.t_start)/details.dt;
}




/// Calculates t from tdiv
double Model::calc_t(double tdiv) const 
{
	return details.t_start + tdiv*details.dt;
}


/// Create precalculation
void Model::create_precalc_equation()
{		
	precalc_eqn.set_all_time();

	// Extracts any precalculation with in the parameter definitions
	param_vec_refq.resize(param_vec.size(),UNSET);
	spline_refq.resize(spline.size(),UNSET);
	
	auto last_spline = get_last_spline();
	
	for(auto th = 0u; th < param_vec.size(); th++){
		auto &pv = param_vec[th];
		const auto &par = param[pv.th];

		switch(par.variety){
		case DIST_PARAM: case PRIOR_PARAM:
			{
				const auto &pri = prior[pv.prior_ref];
				for(const auto &dp : pri.dist_param){
					auto &eq = eqn[dp.eq_ref];
					precalc_eqn.add_eqn(eq.calcu,param_vec_refq,spline_refq,pv.spec_precalc_before);
					eq.precalc_done = true;
					//for(auto i : pv.list_precalc_before) list_precalc.push_back(i);
				}
			}
			break;
				
		case REPARAM_PARAM:
			{
				auto ind = pv.index;
				auto eq_ref = par.get_eq_ref(ind);
			
				if(eq_ref == UNSET) emsg("eq_ref should be set");
				auto &eq = eqn[eq_ref];
				
				precalc_eqn.add_eqn(eq.calcu,param_vec_refq,spline_refq,pv.spec_precalc_before);
		
				if(pv.reparam_spl_ti != UNSET){
					auto d = ind%par.spline_info.knot_tdiv.size();  // Distance down spline
					if(d > 0){  // If not at the begining of the spline then gets information from beginning 
						const auto &er = par.element_ref[ind-d];
						if(er.cons == false){			
							pv.spec_precalc_before = param_vec[par.get_param_vec(ind-d)].spec_precalc_before;
						}
					}
					pv.spec_precalc_before = precalc_eqn.shrink_sprec(pv.reparam_spl_ti,pv.spec_precalc_before);
					//print_spec_precalc("after",pv.spec_precalc_before);
				}
				//for(auto i : pv.list_precalc_before) list_precalc.push_back(i);
			}
			break;
			
		default: break;
		}
		
		param_vec_refq[th] = precalc_eqn.add_param(th);
		
		auto s = last_spline[th];
		if(s != UNSET){
			spline_refq[s] = precalc_eqn.add_spline(s,spec_precalc);
		}
	}

	// Adds any splines not already added (these have no parameters)
	for(auto s = 0u; s < spline_refq.size(); s++){
		if(spline_refq[s] == UNSET){
			spline_refq[s] = precalc_eqn.add_spline(s,spec_precalc);
		}
	}
	
	auto cl = clock();
	
	precalc_eqn.clear_timer();
		
	// Extracts precalculation from all other equations apart from reparameterisation and distribution
	// Multipilier of populations are not combine to allow for pop gradient factor to be extracted

	for(auto &eq : eqn){
		if(eq.precalc_done == false && eq.type != DERIVE_EQN && eq.type != DEFINE_EQN){
			precalc_eqn.add_eqn(eq.calcu,param_vec_refq,spline_refq,spec_precalc);
			eq.precalc_done = true;
		}
	}

	if(false && profiling){
		for(auto i = 0u; i < 20; i++){
			cout << i << " " << get_cpu_time(precalc_eqn.timer[i])  << " time" << endl;
		}
		cout << get_cpu_time(clock()-cl) << " time for precalc" << endl;
	}
}


/// Used to order genetic data
bool Info_ord (const PrecalcInfo &in1, const PrecalcInfo &in2)
{ return (in1.q < in2.q); };  



/// Sets update in precalc at different times (this allows for populations in reparameterised eqns)
void Model::set_spec_precalc_time()
{
	auto T = details.T;
	
	print_diag("h16c1");
	spec_precalc_time_ref.resize(T,UNSET);
	
	for(auto th = 0u; th < param_vec.size(); th++){
		auto &pv = param_vec[th];
		if(pv.reparam_time_dep){
			auto ti = pv.reparam_spl_ti;
			
			if(ti != UNSET){
				auto &spt = get_spec_precalc_time(ti);	
				spt.pv.push_back(th);
			}
		}
	}
	
	vector <unsigned int> par_reparam_time_precalc_affect;
	
	for(auto th = 0u; th < param.size(); th++){
		const auto &par = param[th];
		
		if(par.variety == DYNAMIC_PARAM || par.reparam_time_dep){
			for(auto j = 0u; j < par.spline_segment.size(); j++){
				const auto &seg = par.spline_segment[j];
			
				DynamicParamUpdate dpu;
				dpu.th = th;
				dpu.i = j;
				dpu.ti = seg.ti_min;
				dpu.ti_end = seg.ti_max;
				
				auto &spt = get_spec_precalc_time(dpu.ti);	
				spt.dynamic_param_update.push_back(dpu);
			}
		}
	}

	vector <bool> mapQ(precalc_eqn.pcalcu.size(),false);

	print_diag("h16c2");
	for(auto &spl : spec_precalc_list){   // Update precalc for dynamic/reparam parameters
		auto spec_list = get_dynamic_spec_list(spl,mapQ);
		spl.spec_precalc = precalc_eqn.combine_spec_precalc(spec_list);
		//spl.spec_precalc = precalc_eqn.combine_spec_precalc(vector <unsigned int> (),spec_list);
		//print_spec_precalc("NEW ",spl.spec_precalc);
	}
	
	// Removes those in spec_precalc_time from list_precalc
	auto M = precalc_eqn.pcalcu.size();
	vector <bool> map(M,false);
	
	for(const auto &spt : spec_precalc_list){
		for(const auto &in : spt.spec_precalc.info) map[in.q] = true;
	}
	
	vector <PrecalcInfo> info_new;
	for(const auto &in : spec_precalc.info){
		if(map[in.q] == false) info_new.push_back(in);
	}
	
	print_diag("h16c4");
	
	spec_precalc.info = info_new;
	
	if(false){ 
		for(auto ti = 0u; ti < T; ti++){ 
			auto k = spec_precalc_time_ref[ti];
			if(k != UNSET){
				auto &spt = spec_precalc_list[k];
				
				cout << ti << ": ";		
				for(auto th :  spt.pv) cout << param_vec_name(th) << ",";
				cout << "vector update ";

				for(const auto &dpu : spt.dynamic_param_update){
					cout << param[dpu.th].name << " (" << dpu.ti << " - " << dpu.ti_end << ")" << ",";
				}	
				cout << "dynamic param" << endl;
			
				const auto &spre = spt.spec_precalc;
				if(spre.info.size() > 0){
					print_spec_precalc("precalc",spre);
				}
			}
		}
		emsg("UU");
	}
}
 
 
/// Gets specprecalc for dynamic parameter
vector <SpecPrecalc> Model::get_dynamic_spec_list(const SpecPrecalcTime &spt, vector <bool> &mapQ) const
{		
	const auto &pcalcu = precalc_eqn.pcalcu;
	auto Q = pcalcu.size();

	vector <SpecPrecalc> spec_list;
	
	// Splits them up into the relavant time steps
	vector < vector <unsigned int> > group;
	
	const auto &dpu = spt.dynamic_param_update;
	
	for(auto k = 0u; k < dpu.size(); k++){
		auto i = 0u;
		while(i < group.size() && !(dpu[group[i][0]].ti == dpu[k].ti && dpu[group[i][0]].ti_end == dpu[k].ti_end)) i++;
	
		if(i < group.size()) group[i].push_back(k);
		else{
			vector <unsigned int> gr; gr.push_back(k);
			group.push_back(gr);
		}
	}
	
	if(false){
		for(const auto &gr : group){
			cout << "group: ";
			for(auto i : gr) cout << i << " " << dpu[i].ti << "-" << dpu[i].ti_end << ",";
			cout << endl;
		}
		
		emsg("groups");
	}

	for(const auto &gr : group){
		const auto &dp = dpu[gr[0]]; 
		auto ltime = seq_vec(dp.ti,dp.ti_end);
		SpecPrecalc spec;
		spec.list_time.push_back(ltime);
		PrecalcInfo pi; pi.tlist = 0;
		auto &info = spec.info;
		
		if(gr.size() == 1){
			const auto &par = param[dp.th];
			auto s_start = par.spline_ref;
	
			auto ntimes = par.spline_info.knot_tdiv.size();
			auto N = par.N/ntimes;
		
			for(auto j = 0u; j < N; j++){
				const auto &spl = spline[s_start+j];
				for(auto q : spl.dynamic_precalc){
					pi.q = q;
					info.push_back(pi);
				}
			}
		}
		else{
			for(auto g = 0u; g < gr.size(); g++){
				const auto &dpg = dpu[gr[g]]; 
				
				const auto &par = param[dpg.th];
				auto s_start = par.spline_ref;
		
				auto ntimes = par.spline_info.knot_tdiv.size();
				auto N = par.N/ntimes;
			
				for(auto j = 0u; j < N; j++){
					const auto &spl = spline[s_start+j];
					for(auto q : spl.dynamic_precalc){
						mapQ[q] = true;
					}
				}
			}
			
			for(auto q = 0u; q < Q; q++){
				if(mapQ[q]){
					mapQ[q] = false;
					pi.q = q;
					info.push_back(pi);
				}
			}
		}
		
		spec_list.push_back(spec);

		//print_spec_precalc("here",spec);
	}
	
	return spec_list;
}

	
/// Gets spec_precalc_time (or creates if it doesn't exist)
SpecPrecalcTime& Model::get_spec_precalc_time(unsigned int ti)
{
	auto k = spec_precalc_time_ref[ti];
	if(k == UNSET){
		k = spec_precalc_list.size();
		spec_precalc_time_ref[ti] = k;
		SpecPrecalcTime spre;
		spre.spec_precalc.hash.off();
		spec_precalc_list.push_back(spre);
	}
	
	return spec_precalc_list[k];
}

				
/// Sets the initial value for precalc
void Model::set_precalc_init()
{        
	auto T = details.T;

	precalc_init = precalc_eqn.calculate_precalc_init(spec_precalc);
	
	if(false){
		for(auto i = 0u; i < precalc_eqn.pcalcu.size(); i++){
			cout << i << ": ";
			const auto &ca = precalc_eqn.pcalcu[i];
			auto ii = ca.iref;
			
			if(ca.time_dep){
				auto num = 0u;
				for(auto ti = 0u; ti < T; ti++){
					if(precalc_init[ii+ti] == UNSET) num++;
				}
				if(num == T) cout << "all unset";
				else{
					if(num == 0) cout << "set";
					else cout << num << " unset";
				}
			}
			else{
				if(precalc_init[ii] == UNSET) cout << "unset";
				else cout << "set";
			}
			cout << endl;
		}
		
		cout << "precalc init" << endl;
		emsg("precalc init");
	}
}
	
	
/// Sets popcombw_not_const which give all non-constant popcombw
void Model::set_popcombw_not_const()
{
	for(auto i = 0u; i < npopcombw; i++){
		if(popcombw[i].it.type != NUMERIC) popcombw_not_const.push_back(i);
	}
}

	
/// Sets precalculation to be done after sampling 
void Model::set_spec_precalc_sample()
{
	spec_precalc_sample = precalc_eqn.calculate_spec_precalc_sample(spec_precalc);
}

	
/// Creates any precalulation for derived quantities
void Model::create_precalc_derive()
{	
	for(auto &eq : eqn){
		if(eq.precalc_done == false && eq.type == DERIVE_EQN){		
			precalc_eqn.add_eqn(eq.calcu,param_vec_refq,spline_refq,spec_precalc_derive);
			
			for(auto &inte : eq.integral){		
				precalc_eqn.add_eqn(inte.calc,param_vec_refq,spline_refq,spec_precalc_derive);
			}
			
			eq.precalc_done = true;
		}
	}
}


/// Adds list ot list_time and returns reference
unsigned int Model::get_ti_ref(const vector <unsigned int> &list, vector < vector <unsigned int> > &list_time, Hash &hash_list_time) const 
{
	if(list.size() == details.T) return ALL_TIME_STEP;
		
	auto j = hash_list_time.existing(list);
	if(j == UNSET){
		j = list_time.size();
		list_time.push_back(list);
		hash_list_time.add(j,list);
	}
	
	return j;
}


/// Works out how precalculation is affected by changes in parameters
void Model::precalc_affectQ()
{
	print_diag("h15d1");

	const auto &pcalcu = precalc_eqn.pcalcu;
	auto Q = pcalcu.size();
	auto T = details.T;
	
	// Works out which precalc affect popcombw
	vector <unsigned int> affect_popcombw(Q,UNSET);
	for(auto i = 0u; i < popcombw.size(); i++){
		const auto &it = popcombw[i].it;
		switch(it.type){
		case REG_PRECALC: affect_popcombw[it.num] = i; break;
		case REG_PRECALC_TIME: emsg("should not be time dep"); break;
		default: break;
		}
	}

	// Works out which precalc affect markov equation
	vector < vector < vector <unsigned int> > > me_precalc; // Stores where me gets information from
	me_precalc.resize(species.size());
	for(auto p = 0u; p < species.size(); p++){
		const auto &sp = species[p];
		
		me_precalc[p].resize(sp.markov_eqn.size());
		for(auto i = 0u; i < sp.markov_eqn.size(); i++){
			const auto &me = sp.markov_eqn[i];
						
			const auto &eq = eqn[me.eqn_ref];
			for(const auto &ca : eq.calcu){
				for(const auto &it : ca.item){
					switch(it.type){
						case REG_PRECALC: case REG_PRECALC_TIME:
							me_precalc[p][i].push_back(it.num);
							break;
						
						case POPCOMB:
							{
								const auto &pc = popcomb[it.num];
								for(const auto &el : pc.ele){
									const auto &pcw = popcombw[el.wref];
									const auto &it2 = pcw.it;
									
									switch(it2.type){
									case REG_PRECALC: case REG_PRECALC_TIME:
										me_precalc[p][i].push_back(it2.num);
										break;
										
									default:
										break;
									}
								}
							}
							break;
							
						default:
							break;
					}
				}					
			}
		}
	}
	
	if(false){
		for(auto i = 0u; i < param_vec.size(); i++){
			cout << param_vec_name(i) << " " << param_vec_refq[i] <<  " param_vec_ref" << endl;
		}
		
		for(auto i = 0u; i < spline.size(); i++){
			cout << spline_name(i) << " " << spline_refq[i] <<  " spline_ref" << endl;
		}
	}
	
	// Creates a list of which elements are subsequently affected
		
	Hash hash_list_time;
	vector < vector <unsigned int> > list_time;
	
	vector < vector <AffectQ> > affectQ;
	affectQ.resize(Q);

	for(auto s = 0u; s < spline.size(); s++){ // Makes connection between spline and parameter vec
		auto q = spline_refq[s];
		const auto &spl = spline[s];
		
		if(!spl.constant){
			for(auto ind = 0u; ind < spl.param_ref.size(); ind++){
				const auto &pr = spl.param_ref[ind];
				if(!pr.cons){
					auto th = pr.index;
					vector <unsigned int> list;
					
					switch(spl.type){
					case LINEAR_SPL: 	
						for(auto ti = 0u; ti < T; ti++){ 
							auto ind2 = spl.div[ti].index;
							if(ind2 == ind || ind2+1 == ind) list.push_back(ti);
						}
						break;
						
					case SQUARE_SPL:
						for(auto ti = 0u; ti < T; ti++){ 
							auto ind2 = spl.div[ti].index;
							if(ind2 == ind) list.push_back(ti);
						}
						break;
						
					case CUBICPOS_SPL: case CUBIC_SPL:	
						for(auto ti = 0u; ti < T; ti++){ 
							list.push_back(ti);
						}
						break;
						
					default: emsg("not here"); break;
					}
					
					auto qp = param_vec_refq[th];
					
					auto ti_ref = get_ti_ref(list,list_time,hash_list_time);
					
					AffectQ afQ;
					afQ.q = q;
					afQ.ti_set = UNSET;
					afQ.ti_ref = ti_ref;
		
					affectQ[qp].push_back(afQ);			
				}
			}
		}
	}

	for(auto q = 0u; q < Q; q++){
		const auto &ca = pcalcu[q];
	
		for(const auto &it : ca.pre_item){
			switch(it.type){	
			case PRE_SPLINEREF:
				break;
			
			case PRE_CONSTSPLINEREF:
				{
					auto s = it.num;
					if(spline[s].dynamic){
						auto qp = spline_refq[s];
						AffectQ afQ;
						afQ.q = q;
						afQ.ti_set = UNSET;
						afQ.ti_ref = UNSET;
						affectQ[qp].push_back(afQ);	
					}
				}
				break;
				
			case PRE_PARAMVEC:
				{
					const auto &pv = param_vec[it.num];
					const auto &par = param[pv.th];
		
					if(par.variety == REPARAM_PARAM){
						auto ti_set = pv.reparam_spl_ti;
						
						auto eq_ref = par.get_eq_ref(pv.index);
						if(eq_ref == UNSET) emsg("eq_ref should be set");
						for(const auto &ca2 : eqn[eq_ref].calcu){
							for(const auto &it2 : ca2.item){
								switch(it2.type){
								case REG_PRECALC: 
									{
										auto qp = it2.num;
										AffectQ afQ;
										afQ.q = q;
										afQ.ti_set = UNSET;
										afQ.ti_ref = UNSET;
										affectQ[qp].push_back(afQ);	
									}
									break;
									
								case REG_PRECALC_TIME: 
									if(ti_set != UNSET){
										auto qp = it2.num;
										AffectQ afQ;
										afQ.q = q;
										afQ.ti_set = ti_set;
										afQ.ti_ref = ALL_TIME_STEP;// ti_ref_set;
										affectQ[qp].push_back(afQ);	
									}
									else emsg("Should not be time dependent");
									break;
									
								case POPCOMB:
									{
										const auto &pc = popcomb[it2.num];
										for(const auto &el : pc.ele){
											const auto &pcw = popcombw[el.wref];
											const auto &it3 = pcw.it;
								
											switch(it3.type){
											case REG_PRECALC: 
												{
													auto qp = it3.num;
													AffectQ afQ;
													afQ.q = q;
													afQ.ti_set = UNSET;
													afQ.ti_ref = UNSET;
													affectQ[qp].push_back(afQ);	
												}	
												break;
								
											case REG_PRECALC_TIME:
												if(ti_set != UNSET){
													auto qp = it3.num;
													AffectQ afQ;
													afQ.q = q;
													afQ.ti_set = ti_set;
													afQ.ti_ref = ALL_TIME_STEP;//ti_ref_set;
													affectQ[qp].push_back(afQ);	
												}
												else emsg("Should not be time dependent");
												break;
											
											default:
												break;
											}
										}
									}
									break;
									
								default: 
									break;
								}
							}
						}
					}
				}
				break;
				
			case PRE_REG: 
				{
					auto qp = it.num;
					if(ca.time_dep){
						AffectQ afQ;
						afQ.q = q;
						afQ.ti_set = UNSET;
						afQ.ti_ref = ALL_TIME_STEP;
						affectQ[qp].push_back(afQ);	
					}
					else{
						AffectQ afQ;
						afQ.q = q;
						afQ.ti_set = UNSET;
						afQ.ti_ref = UNSET;
						affectQ[qp].push_back(afQ);	
					}
				}
				break;
				
			case PRE_REGTIME:
				{
					auto qp = it.num;
					if(ca.time_dep){
						AffectQ afQ;
						afQ.q = q;
						afQ.ti_set = UNSET;
						afQ.ti_ref = ALL_TIME_STEP;
						affectQ[qp].push_back(afQ);	
					}
					else{
						emsg("Must be time dep");
					}
				}
				break;
			
			default: break;
			}	
		}
	}
	
	if(false){
		for(auto q = 0u; q < Q; q++){
			cout << q << ": ";
			for(const auto &aQ : affectQ[q]){
				cout << aQ.q << " " << aQ.ti_set << " " << aQ.ti_ref << ",  ";
			}
			cout << endl;
		}
		emsg("affectQ");
	}
	
	if(false){
		for(auto i = 0u; i < param_vec_refq.size(); i++){
			cout << param_vec_name(i) << ":";
			auto qp = param_vec_refq[i];
			for(auto afQ : affectQ[qp]) cout << afQ.q << " " << afQ.ti_set << " " << afQ.ti_ref << "   ";
			cout << endl;
		}
		
		for(auto i = 0u; i < list_time.size(); i++){
			cout << i << ": ";
			for(auto ti : list_time[i]) cout << ti << ",";
			cout << "time" << endl;
		}
	}
	
	vector < vector <unsigned int> > ti_ref_list;
	ti_ref_list.resize(Q);
	vector <bool> map_pop(pop.size(),false);
	vector <bool> map_popcombw(popcombw.size(),false);
	vector <bool> map_popcomb(popcomb.size(),false);

	print_diag("h15d2");
	
	// Works out precalc after a param_vec element
	for(auto j = 0u; j < param_vec.size(); j++){
		add_affect_likeQ(param_vec_refq[j],ti_ref_list,map_pop,map_popcombw,map_popcomb,param_vec[j],affect_popcombw,me_precalc,affectQ,list_time);
		//print_spec_precalc(pv.name,pv.spec_precalc_after);
	}

	/*
	// For constants in reparm works out precalc
	for(auto th = 0u; th < param.size(); th++){
		auto &par = param[th];
		if(par.reparam_time_dep){
			auto cons = false;
			for(const auto &ele : par.element_ref){
				if(ele.cons) cons = true;
			}
			
			if(cons){
				vector <bool> map(Q,false);
				vector <unsigned int> list;
				
				auto s_start = par.spline_ref;
		
				auto ntimes = par.spline_info.knot_tdiv.size();
				auto N = par.N/ntimes;
	
				for(auto j = 0u; j < N; j++){
					auto q = spline_refq[s_start+j];
					list.push_back(q); map[q] = true;
				}
				
				auto k = 0u;
				while(k < list.size()){
					auto q = list[k];
					for(const auto &afQ : affectQ[q]){
						auto qq = afQ.q;
						if(map[qq] == false){
							par.reparam_time_precalc_affect.push_back(qq);
							list.push_back(qq); map[qq] = true;
						}
					}
					k++;
				}
				
				
				cout << par.name << " ADD CON: ";
				for(auto q : par.reparam_time_precalc_affect) cout << q << ",";
				cout << endl;
			}
		}
	}
	*/
	
	print_diag("h15d3");
	
	add_dynamic_spline(affectQ);
	
	/*
	{
		auto vec = seq_vec(param_vec.size());
		popcomb_update_sample = get_popcomb_update(vec);
	}
	*/
	
	print_diag("h15d4");

	// Makes sure list_precalc does not have repeated value
	vector <bool> mapl(Q,false);
	
	auto info_st = spec_precalc.info;
	spec_precalc.info.clear();
	for(const auto &in : info_st){
		if(mapl[in.q] == false){
			mapl[in.q] = true;
			spec_precalc.info.push_back(in);
		}
	}		
	
	if(testing){  // Checks that there is no overlap between different precalc 
		vector <bool> mapl(Q,false);
		for(const auto &pv : param_vec){
			for(const auto &in : pv.spec_precalc_before.info){
				if(mapl[in.q] == true) emsg("double1");
				mapl[in.q] = true;
			}
			
			for(const auto &in : pv.set_param_spec_precalc.info){
				if(mapl[in.q] == true) emsg("double2");
				mapl[in.q] = true;
			}
			
			for(const auto &in : pv.spec_precalc_after.info){
				if(mapl[in.q] == true) emsg("double3");
				mapl[in.q] = true;
			}
			
			for(const auto &in : pv.spec_precalc_before.info) mapl[in.q] = false;
			for(const auto &in : pv.set_param_spec_precalc.info) mapl[in.q] = false;
			for(const auto &in : pv.spec_precalc_after.info) mapl[in.q] = false;
		}
	}
		
	//if(debugging && false){
	if(false){
		cout << "EXTRACT calculations" << endl;
		for(auto e = 0u; e < eqn.size(); e++){
			if(e < 2){
				const auto &eq = eqn[e];
				if(eq.calcu.size() > 0) eq.print_calculation();
			}
		}
	}

	if(slow_check){ // Checks maps are correctly defined
		cout << "test maps" << endl;
		for(const auto &va : ti_ref_list){ if(va.size() != false) emsg("ti_ref_list prob");}
		for(auto va : map_pop){ if(va != false) emsg("map_pop prob");}
		for(auto va : map_popcombw){ if(va != false) emsg("map_popcombw prob");}
		for(auto va : map_popcomb){ if(va != false) emsg("map_popcomb prob");}
	}
	
	if(false) print_precalc();
	
	print_diag("h15d5");
}


/// Adds information for updating dynamic parameters (as well as constants in reparam
void Model::add_dynamic_spline(const vector < vector <AffectQ> > &affectQ)
{
	const auto &pcalcu = precalc_eqn.pcalcu;
	
	vector <bool> map(pcalcu.size(),false);
	
	for(auto s = 0u; s < spline.size(); s++){ // Makes connection between spline and parameter vec
		auto &spl = spline[s];
		const auto &par = param[spl.th];
		if(spl.dynamic || par.reparam_time_dep){
			auto q = spline_refq[s];
			
			vector <unsigned int> list;
			list.push_back(q); map[q] = true;
			
			auto k = 0u;
			while(k < list.size()){
				auto q = list[k];
				const auto &ca = pcalcu[q];
				if(!ca.time_dep) emsg("Should be time dependant");
				for(const auto &afQ : affectQ[q]){
					auto qq = afQ.q;
					if(map[qq] == false){
						spl.dynamic_precalc.push_back(qq);
						list.push_back(qq); map[qq] = true;
					}
				}
				k++;
			}
			
			for(auto q : list) map[q] = false;
		}
	}
	
	if(false){
		for(auto s = 0u; s < spline.size(); s++){ 
			cout << spline_name(s) << " dynamic_precalc: ";
			for(auto q : spline[s].dynamic_precalc) cout << q << ",";
			cout << endl;
		}		
		emsg("dynam");
	}
	
	if(testing){
		for(auto va : map) if(va != false) emsg("map problem");
	}
}
	
	

/// Adds how parameter change results in change in likelihood
void Model::add_affect_likeQ(unsigned int q, vector < vector <unsigned int> > &ti_ref_list, vector <bool> &map_pop, vector <bool> &map_popcombw, vector <bool> &map_popcomb, ParamVecEle &pvec, const vector <unsigned int> &affect_popcombw, const vector < vector < vector <unsigned int> > > &me_precalc, const vector < vector <AffectQ> > &affectQ, const vector < vector <unsigned int> > &list_time)
{		
	auto T = details.T;
	const auto &pcalcu = precalc_eqn.pcalcu;
	
	// Constructs a tree map_PC from the original parameter to see all those precalc affected
	
	vector <unsigned int> list;
	list.push_back(q); ti_ref_list[q].push_back(UNSET);
	auto k = 0u;
	while(k < list.size()){
		auto q = list[k];
		const auto &ca = pcalcu[q];
		if(ca.time_dep){    // Propagates time
			for(const auto &afQ : affectQ[q]){
				auto qq = afQ.q;
				
				auto fl = true;
				auto ti_set = afQ.ti_set;
				if(ti_set != UNSET){ // Checks to see if ti_set is within current times
					fl= false;
					for(auto ti_ref : ti_ref_list[q]){
						if(ti_ref == ALL_TIME_STEP || ti_ref == UNSET){ fl = true; break;}
						else{
							for(auto ti : list_time[ti_ref]){
								if(ti == ti_set){ fl = true; break;}
							}
							if(fl == true) break;
						}
					}
				}
				
				if(fl){
					if(ti_ref_list[qq].size() == 0) list.push_back(qq);
					
					auto ti_ref = afQ.ti_ref;
					if(ti_ref != ALL_TIME_STEP){
						emsg("should be ALL_TIME_STEP");
					}
					
					for(auto ti_ref : ti_ref_list[q]){
						add_to_vec(ti_ref_list[qq],ti_ref);
					}
				}
			}
		}
		else{            // Goes from one parameter to another or spline
			for(const auto &afQ : affectQ[q]){
				auto qq = afQ.q;
				
				if(afQ.ti_set != UNSET) emsg("should not be ti_set");
			
				if(ti_ref_list[qq].size() == 0) list.push_back(qq);
				auto ti_ref = afQ.ti_ref;
				if(ti_ref == ALL_TIME_STEP) ti_ref_list[qq].clear();
				add_to_vec(ti_ref_list[qq],ti_ref);
			}
		}
		k++;
	}
	
	pvec.spec_precalc_after = precalc_eqn.spec_from_ti_ref_list(list,ti_ref_list,list_time);
	
	// Deals with affect on pop, popcombw and popcomb
	vector <unsigned> list_pop;
	vector <unsigned> list_popcombw;
	vector <unsigned> list_popcomb;
	
	// Deals individual effects
	const auto &par = param[pvec.th];
	for(const auto &iegr : par.ieg_ref){
		const auto &sp = species[iegr.p];
		
		const auto &ieg = sp.ind_eff_group[iegr.i];
		auto N = ieg.list.size();
		auto i = pvec.index%N;
		auto j = pvec.index/N;
		
		if(i == j){ // The populations are only affected if parameter comes from diagonal
			const auto &ie = sp.ind_effect[ieg.list[i].index];
			for(auto po : ie.pop_ref){
				if(map_pop[po] == false){ map_pop[po] = true; list_pop.push_back(po);}
			}
		}
	}
	
	for(const auto &sp : species){
		for(const auto &fe : sp.fix_effect){
			if(fe.th == pvec.th){
				for(auto po : fe.pop_ref){
					if(map_pop[po] == false){ map_pop[po] = true; list_pop.push_back(po);}
				}
			}
		}
	}

	for(auto po : list_pop){
		for(const auto &pcr : pop[po].popcomb_ref){
			auto pc = pcr.pcref;
			if(map_popcomb[pc] == false){
				map_popcomb[pc] = true;
				list_popcomb.push_back(pc);
			}
		}
	}
	
	for(auto q : list){
		auto j = affect_popcombw[q];
		if(j != UNSET){
			if(map_popcombw[j] == false){
				map_popcombw[j] = true;
				list_popcombw.push_back(j);
				
				for(auto &pcref : popcombw[j].pcref){
					auto pc = pcref.pc;
					if(map_popcomb[pc] == false){
						map_popcomb[pc] = true;
						list_popcomb.push_back(pc);
					}
				}
			}
		}
	}
	
	pvec.pop_affect = list_pop;
	pvec.popcombw_affect = list_popcombw;
	pvec.popcomb_affect = list_popcomb;
	
	for(auto j : list_pop) map_pop[j] = false;
	for(auto j : list_popcombw) map_popcombw[j] = false;
	for(auto j : list_popcomb) map_popcomb[j] = false;
	
	if(mode == INF || mode == EXT){
		for(auto p = 0u; p < species.size(); p++){
			const auto &sp = species[p];
			for(auto e = 0u; e < sp.markov_eqn.size(); e++){
				const auto &me = sp.markov_eqn[e];
				
				auto fl = false;
				
				vector <bool> map_time(T,false);
				for(auto q : me_precalc[p][e]){
					for(auto ti_ref : ti_ref_list[q]){
						fl = true;
						if(ti_ref == ALL_TIME_STEP || ti_ref == UNSET){
							for(auto ti = 0u; ti < T; ti++) map_time[ti] = true;
						}
						else{
							for(auto ti : list_time[ti_ref]) map_time[ti] = true;
						}
					}
				}

				AffectLike al; 
				if(sp.type != DETERMINISTIC){
					if(me.time_vari == true) al.map = map_time;
					else al.map.resize(1,true);
		
					al.type = MARKOV_VALUE_AFFECT; al.num = p; al.num2 = e;
					param_vec_add_affect(pvec.affect_like,al);
				}
				
				if(fl){				
					switch(sp.type){
					case INDIVIDUAL:
						{
							al.type = MARKOV_LIKE_AFFECT;
							param_vec_add_affect(pvec.affect_like,al);
						}
						break;
						
					case POPULATION:
						{
							AffectLike al; al.map = map_time;
							al.type = MARKOV_POP_AFFECT; al.num = p; al.num2 = e;
							param_vec_add_affect(pvec.affect_like,al);
						}
						break;
						
					case DETERMINISTIC:
						break;
					}
				}
			}
		}
	}
	
	for(auto q : list) ti_ref_list[q].clear();
}





/* This is a note about precalculation quantities
In model
spec_precalc = Stores precalcultion from equations and gradients (not reparam_time_dep) 
spec_precalc_derive = Stores list for precalcultion of derived (for integrals)
spec_precalc_sample = Collects all the precalc to calculate after sample (not reparam_time_dep)     
     
Within param_vec
spec_precalc_before = Precalculation which need to be done before parameter evaluated (when reparam)
set_param_spec_precalc = Sets the parameter/spline values within precalc
spec_precalc_after = All precalculation which need to be done after parameter

spec_precalc_list stores infomation about precalc at different ti
 -> pv = which param_vec are updated
 -> spec_precalc = combines together spec_precalc_after from all pv
 */
			
/// Prints information about precalculation
void Model::print_precalc() const
{
	auto T = details.T;
	
	precalc_eqn.print_calc();
	
	{
		ofstream fout("soec_precalc.txt");
		
		fout << "MODEL" << endl;
		fout << str_spec_precalc("spec_precalc",spec_precalc);
		fout << str_spec_precalc("spec_precalc_derive",spec_precalc_derive);
		fout << str_spec_precalc("spec_precalc_sample",spec_precalc_sample);
		//str_spec_precalc("spec_precalc_all",spec_precalc_all);
		
		//return;
		
		//precalc_eqn.print_calc(); 
		fout << "PARAMVEC" << endl;
		auto thmax = param_vec.size();
		if(thmax > 100) thmax = 100;
		
		for(auto th = 0u; th < thmax; th++){
			const auto &pv = param_vec[th];
		
			fout << "<<" << param_vec_name(th) << ">>" << endl;
			fout << str_spec_precalc("spec_precalc_before",pv.spec_precalc_before);
			fout << str_spec_precalc("set_param_spec_precalc",pv.set_param_spec_precalc);
			fout << str_spec_precalc("spec_precalc_after",pv.spec_precalc_after);
			fout << endl;
		} 
		if(thmax == 100) fout << "..." << endl;
		
		
		//print_spec_precalc("spec_precalc",spec_precalc);
		if(spec_precalc_time_ref.size() == 0) return;
		
		for(auto ti = 0u; ti < T; ti++){
			auto r = spec_precalc_time_ref[ti];
			if(r != UNSET){
				const auto &spl = spec_precalc_list[r];
				fout << "UPDATE time " << ti << endl;
				for(auto th : spl.pv){
					const auto &pv = param_vec[th];
					fout << param_vec_name(th) << " update" << endl;
					fout << str_spec_precalc("spec_precalc_before",pv.spec_precalc_before);
					fout << str_spec_precalc("set_param_spec_precalc",pv.set_param_spec_precalc);
					fout << endl;
				}
				fout << str_spec_precalc("timevary end ",spl.spec_precalc);
				fout << endl;
			}
		}
	}
	
	emsg("PRECALC");
}


/// Detemines if a parameter is the last on along a spline
vector <unsigned int> Model::get_last_spline() const
{
	vector <unsigned int> last_spline(param_vec.size(),UNSET);
	
	for(auto s = 0u; s < spline.size(); s++){
		const auto &spl = spline[s];
		
		if(!spl.constant){
			for(int i = spl.param_ref.size()-1; i >= 0; i--){
				const auto &pr = spl.param_ref[i]; 
				if(!pr.cons){
					last_spline[pr.index] = s;
					break;
				}
			}
		}
	}
	
	return last_spline;
}


/// Adds a constant to the model
unsigned int Constant::add(double val)
{
	auto vec = hash.get_vec_double(val);
	auto j = hash.existing(vec);
	if(j == UNSET){
		j = value.size();
		hash.add(j,vec);
		value.push_back(val);
	}
	
	return j;
}


/// The log probability of sampling from a distribution
bool Model::in_bounds(double x, unsigned int j, const vector <double> &precalc) const
{				
	const auto &pv = param_vec[j];
	const auto &par = param[pv.th];
	
	switch(par.variety){
	case CONST_PARAM: case REPARAM_PARAM: 
		return true;
	
	case DIST_PARAM: case PRIOR_PARAM:
		{
			const auto &pri = prior[pv.prior_ref];
		
			switch(pri.type){
			case MVN_DEFAULT_PR: case MVN_NORM_LKJ_PR: case MVN_UNIFORM_LKJ_PR: case MVN_INV_WISH_PR:
				if(x <= VAR_MIN || x > VAR_MAX) return false;
				return true;
				
			case MVN_JEF_PR: case MVN_UNIFORM_PR: 
				if(x < TINY) return false;
				break;
				
			case MVN_COR_PR: 
				if(x < -COR_MAX || x > COR_MAX) return false;
				break;
				
			case INVERSE_PR: case POWER_PR:
				{	
					auto min = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
					auto max = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
					if(min <= 0 || max <= 0) return false;
					if(x < min || x > max) return false;
				}
				break;
				
			case UNIFORM_PR:
				{	
					auto min = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
					auto max = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
					if(x < min || x > max) return false;
				}
				break;
				
			case EXP_PR:
				{	
					auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
					if(x < 0 || mean <= 0) return false;
				}
				break;
			
			case NORMAL_PR:
				{	
					auto sd = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
					if(sd <= 0) return false;
				}
				break;
				
			case GAMMA_PR:
				{	
					auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
					auto cv = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
					if(x <= 0 || mean <= 0 || cv <= 0) return false;
				}
				break;
			
			case LOG_NORMAL_PR:
				{	
					auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
					auto cv = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
					if(x <= 0 || mean <= 0 || cv <= 0) return false;
				}
				break;
				
			case BETA_PR:
				{	
					auto alpha = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
					auto beta = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
					if(x <= 0 || x >= 1 || alpha <= 0 || beta <= 0) return false;
				}
				break;
			
			case BERNOULLI_PR:
				{	
					auto z = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
					if(x != 0 && x != 1) return false;
					if(z < 0 || z > 1) return false;
				}
				break;
				
			case FIX_PR:
				{
					auto val = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
					if(x != val) return false;
				}
				break;
				
			case DIRICHLET_PR:
				{
					auto alpha = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
					if(x <= 0 || alpha <= 0) return false; 
				}
				break;
				
			case MDIR_PR:
				emsg("Should not be mdir");
				break;
				
			case UNSET_PR:
				emsg("Should not be unset"); 
				break;
			}
		}
		break;
	
	case DYNAMIC_PARAM: emsg("error dynamic param"); break;
	case UNSET_PARAM: emsg("error param"); break;
	}
	
	return true;
}


/// Determines if the prior is bounded or not
bool Model::is_prior_bounded(unsigned int th) const
{
	auto pref = param_vec[th].prior_ref;
	if(pref == UNSET) return false;
	
	switch(prior[param_vec[th].prior_ref].type){
		// SORT
	case INVERSE_PR: case UNIFORM_PR: case POWER_PR: case MVN_JEF_PR: case MVN_UNIFORM_PR:
		return true;	
		
	default:
		return false;	
	}
}


/// Determines if a parameter is a matrix
bool Model::is_matrix(const Param &par) const 
{
	if(par.dep.size() == 2 && par.dep[0].index == par.dep[1].index) return true;
	return false;
}


/// Determines if parameter is symmetric
bool Model::is_symmetric(const Param &par) const
{
	if(is_matrix(par) && param_omega(par.name)) return true;
	return false;
}


/// Converts any fixed priors to constants
void Model::convert_fix_pr_const()
{
	for(auto &par : param){
		for(auto i = 0u; i < par.N; i++){
			auto &er = par.element_ref[i];
			auto ind = er.index;
			if(ind != UNSET && !er.cons){
				const auto &el = par.element[ind];
				auto prr = el.prior_ref;
				if(prr != UNSET){
					const auto &pri = prior[prr];	
					if(pri.type == FIX_PR){
						const auto &dp = pri.dist_param[0];
						auto num = number(dp.te_raw);
						if(num != UNSET){
							er.cons = true;
							er.index = par.add_cons(num);
						}
					}
				}
			}
		}
	}
}


/// Prints a specified precalculation
void Model::print_spec_precalc(string st, const SpecPrecalc &spre) const 
{
	cout << str_spec_precalc(st,spre);
}	


/// Converts a time range into a string
string Model::str_time_range(const vector <unsigned int> &lt) const
{
	stringstream ss;
	if(lt.size() == 0) ss << "no time";
	else{	
		auto start = lt[0];
		ss << start;
		for(auto j = 1u; j < lt.size(); j++){
			if(lt[j] != lt[j-1]+1) ss << "-" << lt[j-1] << "|" << lt[j];
		}
		ss << "-" << lt[lt.size()-1];
	}
	
	return ss.str();
}


/// Prints a specified precalculation
string Model::str_spec_precalc(string st, const SpecPrecalc &spre) const 
{
	stringstream ss;
	ss << st << ": ";
	
	auto imax = spre.info.size();
	if(imax > 100) imax = 100;
	
	for(auto i = 0u; i < imax; i++){
		const auto &in = spre.info[i];
		ss << in.q;
		if(in.tlist != UNSET){
			if(in.tlist == ALL_TIME_STEP) ss << "(all time)";
			else{
				const auto &lt = spre.list_time[in.tlist];	
				ss << "(" << str_time_range(lt) << ")";
			}
		}
		ss << ",";
	}
	ss << endl;
	if(imax == 100) ss << "..." << endl;
	
	return ss.str();
}


/// Sets precalculation for parameter value
void Model::set_param_spec_precalc()
{
	for(auto th = 0u; th < param_vec.size(); th++){
		auto &pv = param_vec[th];
		auto spl_fl = false; if(pv.spline_ref != UNSET) spl_fl = true;
		precalc_eqn.set_param(pv.set_param_spec_precalc,pv.spec_precalc_after,spl_fl);
	}
}


/// If model contains tvreparam then create a reference from pop to th
void Model::set_pop_reparam_th()
{
	pop_reparam_th_on = false;
	if(!contains_tvreparam) return;
	
	//auto P = pop.size();
	auto P = npop;
	auto T = details.T;
	
	pop_reparam_th.resize(P);
	for(auto po = 0u; po < P; po++){
		pop_reparam_th[po].resize(T);
	}
	
	for(auto th = 0u; th < param_vec.size(); th++){
		const auto &pv = param_vec[th];
		if(pv.reparam_time_dep){
			const auto &par = param[pv.th];
			auto eq_ref = par.get_eq_ref(pv.index);
			if(eq_ref == UNSET) emsg("eq_ref should be set");

			const auto &eq = eqn[eq_ref];
			for(const auto &ca : eq.calcu){
				for(const auto &it : ca.item){
					switch(it.type){
					case POPNUM:
						add_to_vec(pop_reparam_th[it.num][pv.reparam_spl_ti],th);
						break;
						
					case POPNUMTIME:
						{
							const auto &ptr = eq.pop_time_ref[it.num];
							if(ptr.ti > pv.reparam_spl_ti){
								auto t1 = calc_t(ptr.ti);
								auto t2 = calc_t(pv.reparam_spl_ti);
								alert_input("Equation '"+eq.te_raw+"' cannot be evaluated because time "+tstr(t1)+" is after the time "+tstr(t2)+" when this expressions is calculated");
							}	
							add_to_vec(pop_reparam_th[ptr.po][ptr.ti],th);		
						}
						break;
						
					default:
						break;
					}
				}
			}
		}						
	}

	for(auto po = 0u; po < P; po++){
		for(auto ti = 0u; ti < T; ti++){
			if(pop_reparam_th[po][ti].size() > 0) pop_reparam_th_on = true;
		}
	}		
	
	if(!pop_reparam_th_on){
		pop_reparam_th.clear();
	}
			
	if(false){
		for(auto po = 0u; po < P; po++){
			for(auto ti = 0u; ti < T; ti++){
				for(auto th : pop_reparam_th[po][ti]){
					cout << ti << " " << pop[po].name << " " << param_vec_name(th) << endl;
				}
			}
		}
		emsg("pop reparam th");
	}
}


/// Loads up proposal information from the file
string Model::load_prop_info(unsigned int ch, const vector <string> &lines)
{
	TerminalInfo ti;
	ti.ch = ch;
	for(auto i = 0u; i < lines.size(); i++){
		auto spl = split(lines[i],'|');
			
		auto type = spl[0];
		
		if(type == "covar_matrix"){	
			if(spl.size() != 7) return "not 5";
		
			ti.n = number(spl[1]);
			if(ti.n == UNSET) return "n unset";

			ti.n_start = number(spl[2]);
			if(ti.n_start == UNSET) return "n_start unset";

			auto num = ti.n-ti.n_start;

			{
				auto vec = get_number_vec(spl[3]);
				for(auto va : vec) if(va == UNSET) return "av unset";
				auto N = vec.size();
				ti.av = vec;
				for(auto &va : ti.av) va *= num;
				
				auto vec2 = get_number_vec(spl[4]);
				for(auto va : vec2) if(va == UNSET) return "av2 unset";
				
				if(vec2.size() != N*N) return "Not square";
				
				ti.av2.resize(N);
				for(auto j = 0u; j < N; j++){
					for(auto k = 0u; k < N; k++) ti.av2[j].push_back((vec2[j*N+k]+vec[j]*vec[k])*num);
				}
			}
			
			{
				auto vec = get_number_vec(spl[5]);
				for(auto va : vec) if(va == UNSET) return "log_av unset";
				auto N = vec.size();
				ti.log_av = vec;
				for(auto &va : ti.log_av) va *= num;
				
				auto vec2 = get_number_vec(spl[6]);
				for(auto va : vec2) if(va == UNSET) return "log_av2 unset";
				
				if(vec2.size() != N*N) return "Not square";
				
				ti.log_av2.resize(N);
				for(auto j = 0u; j < N; j++){
					for(auto k = 0u; k < N; k++) ti.log_av2[j].push_back((vec2[j*N+k]+vec[j]*vec[k])*num);
				}
			}
			
			/*
			print_vector("in av",ti.av);
			print_matrix("in av2",ti.av2);
			print_vector("in logav",ti.log_av);
			print_matrix("in logav2",ti.log_av2);
			*/
		}
		else{
			if(spl.size() != 4) return "not 4";
		
			auto j = 0u;
			while(j < prop_info_str.size() && prop_info_str[j] != type) j++;
			if(j == prop_info_str.size()) emsg("could not find");
		
			PropInfo pi;
			pi.type = prop_info_list[j];
			auto vec = get_number_vec(spl[1]);
			for(auto va : vec){
				if(va == UNSET) return "va uns";
				pi.id.push_back(va);
			}
			pi.value = number(spl[2]);
			pi.vec = get_number_vec(spl[3]);
			for(auto va : pi.vec){
				if(va == UNSET) return "va unse2";
			}
			
			ti.prop_info_store.push_back(pi);
		}
	}
	
	terminal_info.push_back(ti);
	
	return "";
}


/// Works out where the end bracket is
unsigned int Model::get_end_bracket(string &te, unsigned int i)
{
	while(i < te.length() && te.substr(i,1) != "(") i++;
	if(i == te.length()) return i;
	i++;
	
	auto num = 1u;
	while(i < te.length()){
		while(i < te.length() && te.substr(i,1) != "(" && te.substr(i,1) != ")") i++;
		if(i == te.length()) return i;
	
		if(te.substr(i,1) == "(") num++;
		else{
			num--;
			if(num == 0) return i;
		}
		i++;
	}
	
	return i;
}


/// Determines if in integral
bool Model::in_integral(unsigned int i, const vector <SumRange> &int_range) const
{
	for(const auto &ir : int_range){
		if(i > ir.i_start && i < ir.i_end) return true;
	}
	return false;
}
		
		
/// Works out dependency for an equation
vector <string> Model::equation_dep(string te, string &warn)
{
	vector <string> dep;
	
	// Sum range
	vector <SumRange> sum_range;
	
	{
		auto i = 0u;
		while(i < te.length()){
			while(i < te.length() && !in_text(te,i,"Σ")) i++;
			if(i < te.length()){
				while(i < te.length() && te.substr(i,1) != "_") i++;
				if(i < te.length()){
					i++;
					auto di = get_dep_info(te,i,"(");

					if(di.warn != ""){
						warn = "An error occured: "+di.warn;
						return dep;
					}

					for(auto de : di.spl){				
						SumRange ra;
						ra.i_start = i;
						ra.i_end = get_end_bracket(te,i);
						ra.index = de;
						sum_range.push_back(ra);
					}
				}
			}
		}
	}
	
	// Integral range
	vector <SumRange> int_range;
	{
		auto i = 0u;
		while(i < te.length()){
			while(i < te.length() && !in_text(te,i,"∫")) i++;
			if(i < te.length()){
				SumRange ra;
				ra.i_start = i;
				ra.i_end = get_end_bracket(te,i);
				ra.index = "t";
				int_range.push_back(ra);
				i++;
			}
		}
	}
	
	if(false){
		for(auto sr : sum_range){
			cout << sr.i_start << " " << sr.i_end << " "<< sr.index << "  sum range" << endl;
		}
		
		for(auto sr : int_range){
			cout << sr.i_start << " " << sr.i_end << " "<< sr.index << "  int range" << endl;
		}
	}
	
	auto time_dep = false;
	
	auto i = 0u;
	while(i < te.length()){
		while(i < te.length() && te.substr(i,1) != "%" && te.substr(i,1) != "{" && !in_text(te,i,"Σ")) i++;
		if(i < te.length()){
			auto istore = i; 

			string type = "";
			if(in_text(te,i,"Σ")){  // Changes in a sum max
				while(i < te.length() && te.substr(i,1) != "[" && te.substr(i,1) != "(") i++;
				if(i < te.length() && te.substr(i,1) == "["){
					type = "sum";
					i++;
				}
			}
			else{
				if(te.substr(i,1) == "{"){  // Changes in a population
					if(!in_integral(i,int_range) && !time_specified_pop(i,te)) time_dep = true;

					type = "pop";
					i++;
					auto ist = i;
					while(i < te.length() && te.substr(i,1) != ":" && te.substr(i,1) != "}") i++;
					if(i < te.length() && te.substr(i,1) == ":") i++;
					else i = ist;
				}
				else{                               // Changes in a parameter
					auto j = i; 	
					while(j <= te.length()-3 && te.substr(j,3) != "(t)" && te.substr(j,1) != "$") j++;
					if(j <= te.length()-3 && te.substr(j,3) == "(t)"){
						if(!in_integral(i,int_range)) time_dep = true;		
					}
				
					auto ist = i;
					while(i < te.length() && te.substr(i,1) != "_" && te.substr(i,1) != "$") i++;
					if(i < te.length()){
						if(te.substr(i,1) == "_"){
							type = "param";
							i++;
						}
						else{
							if(te.substr(ist+1,i-ist-1) == "t"){
								if(!in_integral(i,int_range)) time_dep = true;		
							}								
						}
					}
				}
			}
			
			if(type != ""){
				auto end = "$}[<(;"; if(type == "sum") end = ",";
				
				auto ist = i;
				auto di = get_dep_info(te,i,end);

				if(di.warn != ""){
					warn = "An error occured: ";
					if(type == "pop"){
						warn = "Population '{"+te.substr(ist,di.iend-ist)+"}' has an error: ";
					}
					if(type == "param"){
						warn = "Parameter '"+te.substr(istore+1,i-istore-2)+"' has a misspecified dependency: ";
					}
					if(type == "sum"){
						warn = "Sum '"+te.substr(istore+1,i-istore-2)+"' has a misspecified dependency: ";
					}
					warn += di.warn;
					return dep;
				}
		
				for(auto k = 0u; k < di.spl.size(); k++){
					auto de = di.spl[k];
					
					if(type == "pop"){ // Checks to see if compartment
						auto spl = split(de,'|');
						if(spl.size() > 1) de = "";
						else{
							for(const auto &sp : species){
								for(const auto &claa : sp.cla){ 	
									if(claa.hash_comp.find(de) != UNSET) de = "";
								}
							}
						}
					}
					
					if(de != ""){
						auto k = 0u; 
						while(k < sum_range.size() && !(i > sum_range[k].i_start && i < sum_range[k].i_end && sum_range[k].index == de)) k++;
						
						if(k == sum_range.size()) add_to_vec(dep,de);
					}
				}
				i = di.iend;
			}
		}
	}

	if(time_dep) dep.push_back("t");

	if(false){
		cout << endl;
		cout << te << " Equation" << endl;
		for(auto de : dep) cout << "," << de;
		cout << endl;
		emsg("eqdo");
	}
	
	return dep;
}


/// Determines it time-specified pop, e.g. {I ; t=50} 
bool Model::time_specified_pop(unsigned int i, const string &te) const
{
	auto ist = i;
	while(i < te.length() && te.substr(i,1) != "}") i++;
	if(i == te.length()) return false;
	
	auto spl = split(te.substr(ist,i-ist),';');
	
	if(spl.size() > 1){
		const auto &last = spl[spl.size()-1];
		auto sple = split(last,'=');
			
		if(sple.size() == 2){
			if(sple[0] == "t") return true;
		}
	}
	return false;
}
		
		
/// Returns the species number from its name
unsigned int Model::find_p(string name) const 
{
	for(auto p = 0u; p < nspecies; p++){
		if(toLower(species[p].name) == toLower(name)) return p;
	}
	return UNSET;
}


/// Returns the classification number from its name
unsigned int Model::find_cl(unsigned int p, string name) const 
{
	const auto &sp = species[p];
	for(auto cl = 0u; cl < sp.ncla; cl++){
		if(toLower(sp.cla[cl].name) == toLower(name)) return cl;
	}
	return UNSET;
}


/// Determines if in a data mode
bool Model::data_mode() const 
{
	switch(mode){
	case DATA_SIM: case DATA_SHOW: case DATA_DEL: case DATA_CLEAR: 
	case COMPRESS: case DECOMPRESS:
		return true;
	default: break;
	}
	
	return false;
}


/// Determines if model doesn't need to be processed
bool Model::no_process() const 
{
	switch(mode){
	case DATA_SHOW: case DATA_DEL: case DATA_CLEAR: 
	case COMPRESS: case DECOMPRESS:
		return true;
	default: break;
	}
	
	return false;
}


/// Asks a question
bool Model::question(string st) const 
{
	if(no_question) return true;
	
	char type;
	do{
		cout << st << " [y/n]" << endl;
		cin >> type;
		if(type == 'y' || type == 'Y') return true;
		if(type == 'n' || type == 'N') return false;
	}while( !cin.fail());
	return false;
}

/// Gets the classification from the trans name (or if formated using S->E|E->I or S->E:0.5|E->I:1) 
unsigned int Model::get_cl_from_trans(string name, unsigned int p) const
{
	const auto &sp = species[p];
	
	name = replace(name,"->","→");

	auto spl = split_with_bracket(name,'|');
	
	if(spl.size() == 1){	
		auto spl2 = split(name,':');
		for(auto cl = 0u; cl < sp.ncla; cl++){
			const auto &claa = sp.cla[cl];
			for(auto c = 0u; c < claa.ntra; c++){
				if(claa.tra[c].name == spl2[0]) return cl;
			}
		}
	}
	else{
		auto cl_st = UNSET;
		for(auto i = 0u; i < spl.size(); i++){
			auto cl = get_cl_from_trans(spl[i],p);
			if(cl_st == UNSET) cl_st = cl;
			else{
				if(cl_st != cl) return UNSET;
			}
		}
		return cl_st;
	}
	
	return UNSET;
}


/// Gets the classification from the compartment name (or if formated using S|E or S:0.5|E:1) 
unsigned int Model::get_cl_from_comp(string name, unsigned int p) const
{
	const auto &sp = species[p];

	auto spl = split_with_bracket(name,'|');
	
	if(spl.size() == 1){	
		auto spl2 = split(name,':');
		for(auto cl = 0u; cl < sp.ncla; cl++){
			const auto &claa = sp.cla[cl];
			for(auto c = 0u; c < claa.ncomp; c++){
				if(claa.comp[c].name == spl2[0]) return cl;
			}
		}
	}
	else{
		auto cl_st = UNSET;
		for(auto i = 0u; i < spl.size(); i++){
			auto cl = get_cl_from_comp(spl[i],p);
			if(cl_st == UNSET) cl_st = cl;
			else{
				if(cl_st != cl) return UNSET;
			}
		}
		return cl_st;
	}
	
	return UNSET;
}


/// Calculate the equation for a string
double Model::calculate_equation_zero_one(string te, double tdiv, string &err) 
{
	auto val = calculate_equation(te,tdiv,err);
	if(err == ""){
		if(val < 0 || val > 1) err = "Fraction must be between zero and one.";
	}
	
	return val;
}


/// Calculate the equation for a string
double Model::calculate_equation(string te, double tdiv, string &err) 
{
	auto ti = get_ti(tdiv);
	
	auto p = UNSET, cl = UNSET;
	
	auto eqi = add_equation_info(te,MODEL_CALC,p,cl);
	if(eqi.error){
		err = eqi.emsg;
		return UNSET;
	}
	
	Equation eq(eqi,ti,species_simp,comp_pos,hash_comp_pos,param,prior,region,derive,spline,param_vec,density,pop,hash_pop,popcomb,popcombw,constant,timepoint,details,define);
	
	if(eq.warn != ""){
		err = eq.warn;
		return UNSET;
	}
	else{
		vector <double> precalc;
		return eq.calculate_no_popcomb(ti,precalc);
	}
}


/// Uses a string to get information about diag test data (e.g. I[Se:0.5])
DiagTestSens Model::get_diag_test_sens(string comp, unsigned int p, string &warn) const
{
	vector <string> co_list;         // The test-sensitive compartment list
	vector <string> Se_list;         // Sensitivity list
	
	auto cl_sel = UNSET;

	auto spl = split(comp,',');
	for(auto i = 0u; i < spl.size(); i++){
		auto spl2 = split(spl[i],'[');				
		string co, Se;
		if(spl2.size() == 2){
			co = spl2[0];
			auto bra = trim(spl2[1]);
			auto len = bra.length();
			if(len > 0){
				if(bra.substr(len-1,1) == "]"){
					bra = bra.substr(0,len-1);
					auto spl3 = split(bra,':');
					if(spl3.size() == 2){
						if(trim(spl3[0]) == "Se"){
							Se = trim(spl3[1]);
						}
					}
				}
			}						
		}
			
		if(co == "" || Se == ""){	
			warn = "In 'comp' the value '"+spl[i]+"' should be in the format 'compartment(Se=sensitivity)'";
		}
		else{
			auto cl = get_cl_from_comp(co,p);
			if(cl == UNSET){
				warn = "Value '"+co+"' is not a compartment";
			}
			
			if(cl_sel == UNSET) cl_sel = cl;
			else{
				if(cl != cl_sel){
					warn = "In 'comp' cannot have compartments from different classifications.";
				}
			}
	
			co_list.push_back(co);
			Se_list.push_back(Se);
		}
	}
	
	if(cl_sel == UNSET){
		warn = "In 'comp' there is a syntax error.";
	}
	
	DiagTestSens dts;
	dts.cl = cl_sel;
			
	const auto &claa = species[p].cla[cl_sel];
			
	for(auto c = 0u; c < claa.ncomp; c++){
		auto name = claa.comp[c].name;
		
		TestComp tc;
		auto k = find_in(co_list,name);
		
		if(k != UNSET){ tc.on = true; tc.Se_str = Se_list[k];}
		else tc.on = false;
		
		dts.comp.push_back(tc);
	}
		
	return dts;
}


/// Extracts linear combinations of populations from equations
void Model::extract_popcomb(Hash &hashw, Hash &hashpc)
{
	for(auto &eq : eqn){
		if(eq.type != DERIVE_EQN) eq.extract_popcomb(eq.calcu,hashw,hashpc,false);
	}
	
	npopcomb = popcomb.size();
	npopcombw = popcombw.size();
}


/// Prints an item
string Model::print_it(const EqItem &it) const
{	
	stringstream ss; 
	switch(it.type){
	case NUMERIC: ss << constant.value[it.num]; break;
	case REG_PRECALC: ss << "Rpre" << it.num; break;
	case ONE: ss << "1"; break;
	case ZERO: ss << "0"; break;
	default: emsg("Cannot print: "+tstr(it.type)); break;	
	}
	
	return ss.str();
}


/// Prints information about popcombs
void Model::print_popcomb() const
{
	cout << "POPCOMB" << endl;
	for(auto i = 0u; i < popcomb.size(); i++){
		if(i >= npopcomb) cout << "DERIVE";
		cout << i << ": ";
		const auto &popc = popcomb[i].ele;
		for(const auto &pc : popc){
			const auto &pcw = popcombw[pc.wref];
			cout << print_it(pcw.it) << "*" << pop[pc.po].name << " w=" << pc.wref << ", ";
		}				
		
		for(const auto &mer : popcomb[i].markov_eqn_ref){
			cout << mer.p << " species  ";
			for(auto e : mer.list){
				cout << eqn[species[mer.p].markov_eqn[e].eqn_ref].te_raw << ",";
			}
		}
		
		cout << endl;
	}
	
	cout << "popcombw" << endl;
	for(auto i = 0u; i < popcombw.size(); i++){
		if(i >= npopcombw) cout << "DERIVE";
		const auto &pcw = popcombw[i];
		const auto &it = pcw.it;
		cout << "w=" << i << ": " <<  print_it(it) << " ";
		for(const auto &va : pcw.pcref) cout << va.pc << " " << va.index << ", ";
		cout << endl;
	}
	//emsg("G");
	
	for(auto th = 0u; th < param_vec.size(); th++){
		const auto &pv = param_vec[th];
		cout << add_escape_char(param_vec_name(th)) << "  ";
	
		cout << "pop_affect: ";
		for(auto val : pv.pop_affect) cout << val << ",";
		
		cout << "  popcombw_affect: ";
		for(auto val : pv.popcombw_affect) cout << val << ",";
		
		cout << "  popcomb_affect: ";
		for(auto val : pv.popcomb_affect) cout << val << ",";
	
		cout << endl;
	}
	
	/*
	for(auto th = 0u; th < param_vec.size(); th++){
		const auto &pus = popcomb_update_sample[th];
		
		cout << "after " << param_vec[th].name << ": ";
		cout << "popcombw_list: ";
		for(auto val : pus.popcombw_list) cout << val << ", ";
		cout << "  ";
		cout << "popcomb_list: ";
		for(auto val : pus.popcomb_list) cout << val << ", ";
		cout << endl;
	}
	*/
	
	emsg("popcombw");
}


/// Extracts linear combinations of populations from equations for derived quantities
void Model::extract_popcomb_derive(Hash &hashw, Hash &hashpc)
{	
	for(auto &eq : eqn){
		if(eq.type == DERIVE_EQN){	
			eq.extract_popcomb(eq.calcu,hashw,hashpc,true);
			
			for(auto &inte : eq.integral){		
				eq.extract_popcomb(inte.calc,hashw,hashpc,true);
			}
		}
	}
}


/// Check that everything is linear and with only one popcomb
void Model::check_all_linear() const
{
	return;
	
	for(const auto &eq : eqn){
		auto num = 0u; 
		for(const auto &ca : eq.calcu){
			for(const auto &it : ca.item){
				if(it.type == POPCOMB) num++;
			}				
		}
	
		if(eq.markov_eqn_ref != UNSET && eq.lin.on != true){
			emsg("Equation not linear!");
		}
		if(num > 1) emsg("More than one popcomb");
	}
}


/// Sets pop_list and popcomb_list
void Model::set_pop_list()
{
	//for(auto po = 0u; po < pop.size(); po++){
	for(auto po = 0u; po < npop; po++){
		species[pop[po].p].pop_list.push_back(po);
	}
	
	//for(auto po = 0u; po < popcomb.size(); po++){
	for(auto po = 0u; po < npopcomb; po++){
		species[popcomb[po].p].popcomb_list.push_back(po);
	}
}


/// Sets population gradient references
void Model::set_pop_grad_ref()
{
	//vector <unsigned int> map(pop.size(),UNSET);
	vector <unsigned int> map(npop,UNSET);
	vector <unsigned int> list;

	for(auto &eq : eqn){
		auto &lin = eq.lin;
		if(lin.on){
			lin.pop_grad_ref.resize(eq.pop_ref.size());
			
			for(auto i = 0u; i < eq.pop_ref.size(); i++){
				auto po = eq.pop_ref[i];
				map[po] = i;
				list.push_back(po);
			}

			for(auto j = 0u; j < eq.popcomb_ref.size(); j++){
				const auto &pc = popcomb[eq.popcomb_ref[j]];
				for(const auto &el : pc.ele){
					auto i = map[el.po];
					if(i == UNSET) emsg("problem with pop_grad_ref");
					
					PopGradRef pgr;
					pgr.wref = el.wref;
					pgr.popcomb_grad = lin.popcomb_grad_precalc[j];

					lin.pop_grad_ref[i].push_back(pgr);
				}
			}
			
			for(auto va : list) map[va] = UNSET;
			list.clear();
		}
	}
	
	for(auto va : map) if(va != UNSET) emsg("Should be unset");
}



/// Calculates all the precalc quantities (used in diangostics)
void Model::precalc_calculate_all(PV &param_val, const vector < vector <double> > &popcomb_t) const 
{
	for(auto th = 0u; th < param_vec.size(); th++){
		const auto &pv = param_vec[th];
		precalc_eqn.calculate(pv.spec_precalc_before,param_val,false);
		if(!pv.reparam_time_dep){
			precalc_eqn.calculate(pv.set_param_spec_precalc,param_val,false);
		}		
	}
	precalc_eqn.calculate(spec_precalc,param_val,false);
	
	param_spec_precalc_time_all(popcomb_t,param_val,false);
}


// Creates eqn_check
void Model::create_equation_check()
{
	for(const auto &eq : eqn){
		EqnCheck eq_ch;
		eq_ch.calcu = eq.calcu;
		eq_ch.integral = eq.integral;
		eq_ch.param_ref = eq.param_ref;
		eq_ch.derive_ref = eq.derive_ref;
		eqn_check.push_back(eq_ch);
	}
}


/// Checks all equation values against stored simplified version
void Model::check_equation_values(const vector < vector <double> > &popnum_t, const vector < vector <double> > &popcomb_t, const PV &param_val, const vector < vector < vector <double> > > &derive_val) const
{
	cout << "CHECK EQUATIONS" << endl;
	auto T = details.T;

	for(auto e = 0u; e < eqn.size(); e++){
		const auto &eq = eqn[e];	
		cout << e << " / " <<  eqn.size() << " "<< eq.te_raw << endl;
		if(eq.type != DEFINE_EQN){
			const auto &eq_ch = eqn_check[e];
			
			if(eq.time_vari){
				for(auto ti = 0u; ti < T; ti++){
					auto val = eq.calculate_derive(ti,popcomb_t,param_val.precalc,derive_val);
					auto val2 = eq.calculate_check(ti,param_val,popnum_t,derive_val,eq_ch.calcu,eq_ch.integral,eq_ch.param_ref,eq_ch.derive_ref,spline_refi);
		
					if(dif(val,val2,DIF_THRESH)){
						cout << eq.te_raw << " " << ti << " " << val << " " << val2 << " comp ti" << endl;
						emsg("Equation value problem");
					}
				}
			}
			else{
				auto val = eq.calculate_derive(0,popcomb_t,param_val.precalc,derive_val);
				auto val2 = eq.calculate_check(0,param_val,popnum_t,derive_val,eq_ch.calcu,eq_ch.integral,eq_ch.param_ref,eq_ch.derive_ref,spline_refi);
				if(dif(val,val2,DIF_THRESH)){
					for(auto ti = 0u; ti < T; ti++){
						cout << ti << " ";
						for(auto va : popnum_t[ti]) cout << va << ",";
						cout << "    ";
						for(auto va : popcomb_t[ti]) cout << va << ",";
						cout << " pop" << endl;
					}
					cout << eq.te_raw << " " << val << " " << val2 << " comp" << endl;
				
					emsg("Equation value problem");
				}
			}
		}
	}
}


/// Gets the parameter vector name
string Model::param_vec_name(unsigned int th) const
{
	const auto &pv = param_vec[th];
	const auto &par = param[pv.th]; 
	return add_escape_char(get_param_name_with_dep(par,par.dep,pv.index));
}


/// Gets the spline name
string Model::spline_name(unsigned int i) const
{
	const auto &spl = spline[i];
	const auto &par = param[spl.th]; 
	
	auto dep_reduce = par.dep;
	dep_reduce.pop_back();
	
	return get_param_name_with_dep(par,dep_reduce,spl.index);
}


/// Sets the correct numbers for precalculation from q to i
void Model::set_eqn_precalcnum()
{
	for(auto &eq : eqn){
		if(eq.type != DEFINE_EQN){		
			precalc_eqn.adjust_num(eq.calcu);
			for(auto &inte : eq.integral) precalc_eqn.adjust_num(inte.calc);

			auto &lin = eq.lin;
			if(lin.on){			
				precalc_eqn.adjust_it(lin.no_pop_precalc);
				for(auto &it : lin.popcomb_grad_precalc) precalc_eqn.adjust_it(it);
			}
		}
	}
	
	for(auto &pcw : popcombw){
		precalc_eqn.adjust_it(pcw.it);
	}
	
	for(auto &pvr : param_vec_refq){
		param_vec_refi.push_back(precalc_eqn.pcalcu[pvr].iref);
	}
	
	for(auto &sr : spline_refq){
		spline_refi.push_back(precalc_eqn.pcalcu[sr].iref);
	}
	
	param_vec_refq.clear();
	spline_refi.clear();
}


/// Determines if an element is set or not (used for outputting)
bool Model::par_ele_set(unsigned int th, unsigned int j) const
{
	const auto &par = param[th];
	if(par.time_dep){ // If spline then shifts to beginning of spline
		auto ntimes = par.spline_info.knot_tdiv.size();
		j -= j%ntimes;
	}
	
	const auto &ele = par.element_ref[j];
	if(ele.index == UNSET) return false;
	return true;
}


/// Determines if a variable is Bernoulli
bool Model::is_bernoulli(unsigned int th) const 
{
	if(prior[param_vec[th].prior_ref].type == BERNOULLI_PR) return true;
	return false;
}


/// Detemines if variable is strictly positive
bool Model::is_strictly_positive(const Prior &pri) const 
{
	switch(pri.type){
	case INVERSE_PR: case POWER_PR: case EXP_PR: case GAMMA_PR: 
	case LOG_NORMAL_PR: case BETA_PR: case DIRICHLET_PR:
		return true;
		
	case UNIFORM_PR:
		{
			auto val_min = eqn[pri.dist_param[0].eq_ref].is_num();
			if(val_min != UNSET){
				if(val_min > 0) return true;
				else{
					if(val_min == 0){
						auto val_max = eqn[pri.dist_param[1].eq_ref].is_num();
						if(val_max != 0) return true;
					}
				}
			}
		}
		break;
			
	default:
		break;
	}
	
	return false;
}



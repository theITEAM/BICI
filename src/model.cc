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
Model::Model(Operation mode_) : precalc_eqn(species_simp,spline,param_vec,pop,constant,timepoint,details)
{
	mode = mode_;  
	nspecies = 0;
};


/// Adds an equation reference to EquationInfo
void Model::add_eq_ref(EquationInfo &eqi, Hash &hash_eqn, double tdiv)
{
	auto ti = UNSET;
	if(tdiv != UNSET) ti = get_ti(tdiv);
	
	if(eqi.te == "") emsg_input("Equation does not have any text");
	
	auto vec = hash_eqn.get_vec_eqn(eqi.te,(unsigned int)eqi.type,eqi.p,eqi.cl,eqi.infection_trans,ti);

	auto e = hash_eqn.existing(vec);
	if(e != UNSET){
		eqi.eq_ref = e;
	}
	else{
		eqi.eq_ref = eqn.size();	
		
		hash_eqn.add(eqi.eq_ref,vec);
	
		Equation eq(eqi.te,eqi.type,eqi.p,eqi.cl,eqi.infection_trans,ti,eqi.line_num,species_simp,param,prior,derive,spline,param_vec,pop,hash_pop,constant,timepoint,details);
		
		if(false && eq.warn != ""){ cout << eq.warn << endl; emsg_input("warning");}
		
		eqn.push_back(eq);
	}
}


/// Initialises a parameter value
void Model::param_val_init(PV &param_val) const
{
	auto N = nparam_vec;
	auto M = precalc_eqn.calcu.size();
	param_val.value.resize(N,UNSET); 
	//param_val.precalc.resize(M,UNSET);
	param_val.precalc = precalc_init;
	param_val.value_old.resize(N,UNSET);
	param_val.precalc_old.resize(M,UNSET);
	
	//precalc_eqn.calc_spline_const(param_val,spline_ref);
}


/// Samples a set of parameter values from the model
PV Model::param_sample() const
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
			precalc_eqn.calculate(pv.set_param_spec_precalc,param_val,false);
			break;
		
		case REPARAM_PARAM: 
			{
				if(pv.reparam_time_dep == false){
					auto eq_ref = par.get_eq_ref(pv.index);
					if(eq_ref == UNSET) emsg("eq_ref should be set");

					value[th] = eqn[eq_ref].calculate_param(precalc);
					precalc_eqn.calculate(pv.set_param_spec_precalc,param_val,false);
				}
			}
			break;
			
		case DIST_PARAM: case PRIOR_PARAM:	
			value[th] = prior_sample(prior[pv.prior_ref],precalc);
			precalc_eqn.calculate(pv.set_param_spec_precalc,param_val,false);
			break;
			
		case UNSET_PARAM: emsg("error param"); break;
		}
	}
	
	sample_ieg_cv(param_val); // Samples any individual effect covariances 
	
	//precalc_eqn.calculate(spec_precalc,param_val,false);
	precalc_eqn.calculate(spec_precalc_sample,param_val,false);
	
	
	//print_param(param_val);
	
	if(false){
		for(auto val: param_val.value) cout << val << ","; 	
		cout << " value" << endl;
		
		//for(auto val: param_val.precalc) cout << val << ","; 	
		auto imax =  param_val.precalc.size();
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
				if(pri.type != MVN_JEF_PR && pri.type != MVN_UNIFORM_PR){
					emsg("Should be mvn-jeffreys or mvn-uniform");
				}
				
				auto N = ieg.list.size();
				
				auto loop = 0u, loopmax = 1000u;
				for(loop = 0; loop < loopmax; loop++){
					if(sbound){
						emsg("To do");
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
		}
	}
}

/// Updatews precalc for all times
void Model::param_spec_precalc_time_all(const vector < vector <double> > &popnum_t, PV &param_val, bool store) const
{
	for(auto ti = 0u; ti < details.T; ti++){
		param_spec_precalc_time(ti,popnum_t,param_val,store);
	}
}

	
/// Updates parameter precalculation at given time point
void Model::param_spec_precalc_time(unsigned int ti, const vector < vector <double> > &popnum_t, PV &param_val, bool store) const
{
	const auto &upt = spec_precalc_time[ti];
	if(upt.pv.size() > 0){
		auto &value = param_val.value;
		auto &precalc = param_val.precalc;

		for(auto th : upt.pv){
			const auto &pv = param_vec[th];
			const auto &par = param[pv.th];
			auto eq_ref = par.get_eq_ref(pv.index);
			if(eq_ref == UNSET) emsg("eq_ref should be set");
				
			precalc_eqn.calculate(pv.spec_precalc_before,param_val,false);
		
			if(store) param_val.value_change(th);
			value[th] = eqn[eq_ref].calculate_all_time(ti,popnum_t,precalc);
			
			precalc_eqn.calculate(pv.set_param_spec_precalc,param_val,store);
		}
	}
	
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
			precalc_eqn.calculate(pv.set_param_spec_precalc,param_val,false);
			break;
			
		case REPARAM_PARAM:
			if(pv.reparam_time_dep == false){
				auto eq_ref = par.get_eq_ref(pv.index);
				if(eq_ref == UNSET) emsg("eq_ref should be set");
			
				value[th] = eqn[eq_ref].calculate_param(precalc);
				precalc_eqn.calculate(pv.set_param_spec_precalc,param_val,false);
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
			precalc_eqn.calculate(pv.set_param_spec_precalc,param_val,false);
			break;
			
		case UNSET_PARAM: emsg("error param"); break;
		}
	}
	
	precalc_eqn.calculate(spec_precalc,param_val,false);
	
	return param_val;
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
	
	if(par.variety == CONST_PARAM) return UNSET;
	
	auto omega = sp.calculate_omega_basic(iegr.i,param_val,param);
	
	for(auto i = 0u; i < omega.size(); i++){
		if(omega[i][i] > ieg.var_max) return -LARGE;
	}
	
	auto N = ieg.list.size();

	auto log_det = determinant_fast(omega);
	
	if(log_det == UNSET) return -LARGE;
	if(log_det < ieg.log_det_min) return -LARGE;
		
	const auto &pri = prior[ieg.prior_ref];
	
	switch(pri.type){
	case MVN_JEF_PR: 
		{
			return (0.5*N-0.5)*diag_log_sum(omega)-(0.5*N+0.5)*log_det;
		}
		
	case MVN_UNIFORM_PR: return 0;
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
					if(j != i){
						auto val = omega[j][i];
						auto name = exchange_omega(par.name)+"_"+ieg.list[j].name+","+ieg.list[i].name; 
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
	
			auto min = eqn[pri.dist_param[0].eq_ref].calculate_constant();
			auto max = eqn[pri.dist_param[1].eq_ref].calculate_constant();

			if(min <= 0) alert_input("For the prior '"+pri.name+"' the minumum value '"+tstr(min)+"' must be positive.");
			
			if(max <= 0) alert_input("For the prior '"+pri.name+"' the maximum value '"+tstr(max)+"' must be positive.");
			
			if(min >= max){
				alert_input("For the prior '"+pri.name+"' the minumum value '"+tstr(min)+"' is larger than the maximum value '"+tstr(max)+"'.");
			}
			
			ieg.log_det_min = N*log(min);
			ieg.var_max = max;
			
			ieg.omega_pv.resize(N);
			for(auto j = 0u; j < N; j++){
				ieg.omega_pv[j].resize(N);
				for(auto i = 0u; i < N; i++){
					auto val =  par.get_param_vec(j*N+i);
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


/// Used to order affects
bool AL_ord(const AffectLike &al1, const AffectLike &al2)                      
{ return (al1.order_num < al2.order_num); };  


/// Orders a list of AffectLike based on the priority in the calculation
void Model::order_affect(vector <AffectLike> &vec) const
{
	for(auto &al : vec){
		switch(al.type){
		case SPLINE_PRIOR_AFFECT: case IEG_PRIOR_AFFECT:  case PRIOR_AFFECT: case DIST_AFFECT: case OMEGA_AFFECT: 
		case LIKE_INIT_COND_AFFECT: case PRIOR_INIT_COND_AFFECT:
		case GENETIC_VALUE_AFFECT:
			al.order_num = 1; 
			break;
			
		case EXP_FE_AFFECT: case EXP_IE_AFFECT:
			al.order_num = 2; 
			break;
			
		case POP_AFFECT:
			al.order_num = 3; 
			break;
			
		case INDFAC_INT_AFFECT: 
			al.order_num = 4; 
			break;
		
		case DIV_VALUE_AFFECT: case DIV_VALUE_NOPOP_AFFECT: case DIV_VALUE_LINEAR_AFFECT: 
		case POP_DATA_CGL_TGL_AFFECT:
			al.order_num = 5; 
			break;
			
		case OBS_EQN_AFFECT: //case OBS_TRANS_EQN_AFFECT:
		case IIF_W_AFFECT:
		case POPNUM_IND_W_AFFECT:
			al.order_num = 6; 
			break;
			
		case MARKOV_LIKE_AFFECT: case NM_TRANS_AFFECT:  case NM_TRANS_BP_AFFECT: case NM_TRANS_INCOMP_AFFECT: case LIKE_IE_AFFECT: 
		case MARKOV_POP_AFFECT: case MARKOV_POP_NOPOP_AFFECT: case MARKOV_POP_LINEAR_AFFECT: 
		case LIKE_OBS_IND_AFFECT: case LIKE_OBS_POP_AFFECT: 
		case LIKE_OBS_POP_TRANS_AFFECT: case LIKE_UNOBS_TRANS_AFFECT:
		case LIKE_GENETIC_PROCESS_AFFECT: case LIKE_GENETIC_OBS_AFFECT:
			al.order_num = 7; 
			break;
			
		case AFFECT_MAX: break;
		}
	}
	
	sort(vec.begin(),vec.end(),AL_ord);
}


/// Gets a map which stores how the proposal changes the system
AffectMap Model::get_affect_map(vector <AffectLike> &vec, const vector <unsigned int> &param_list, const vector <unsigned int> &dependent, const SpecPrecalc &spec_precalc_after) const
{
	AffectMap amap;
	
	auto &map = amap.precalc_map;
	map.resize(precalc_eqn.calcu.size(),false);
	
	// Creates a map which shows all precalc which are affected by proposal
	
	for(auto th : param_list){
		const auto &pv = param_vec[th];
		for(const auto &in : pv.set_param_spec_precalc.info) map[in.i] = true;
	}
	
	for(auto th : dependent){
		const auto &pv = param_vec[th];
		for(const auto &in : pv.spec_precalc_before.info) map[in.i] = true;
		for(const auto &in : pv.set_param_spec_precalc.info) map[in.i] = true;
	}
	
	for(const auto &in : spec_precalc_after.info) map[in.i] = true;
	
	// Creates maps for fixed and individual effects which change	
	auto &ie_map = amap.ie_map;
	ie_map.resize(nspecies);
	for(auto p = 0u; p < nspecies; p++){
		const auto &sp = species[p];
		ie_map[p].resize(sp.ind_effect.size(),false);
	}
	
	auto &fe_map = amap.fe_map;
	fe_map.resize(nspecies);
	for(auto p = 0u; p < nspecies; p++){
		const auto &sp = species[p];
		fe_map[p].resize(sp.fix_effect.size(),false);
	}
	
	for(auto i = 0u; i < vec.size(); i++){
		switch(vec[i].type){
		case EXP_IE_AFFECT: ie_map[vec[i].num][vec[i].num2] = true; break;
		case EXP_FE_AFFECT: fe_map[vec[i].num][vec[i].num2] = true; break;
		default: break;
		}
	}
	
	return amap;
}
	
	
/// Looks to speed up calculation of DIV_VALUE_AFFECT for parameter proposals
// Calculates changes in non-population term as a way to speed up div value calculation 
void Model::affect_nopop_speedup(vector <AffectLike> &vec, const vector <unsigned int> &param_list, const vector <unsigned int> &dependent, const SpecPrecalc &spec_precalc_after) const
{			
	auto affect_map = get_affect_map(vec,param_list,dependent,spec_precalc_after);
	
	// Generates a list of elements in vec which can be converted	
	vector <Listie>	vec_list;
	
	for(auto i = 0u; i < vec.size(); i++){
		auto &al = vec[i];
		
		auto e = UNSET;
	
		switch(al.type){
		case DIV_VALUE_AFFECT:
			{
				const auto &meq = species[al.num].markov_eqn[al.num2];
				if(meq.rate) e = meq.eqn_ref;
			}
			break;
			
		case MARKOV_POP_AFFECT:
			{
				const auto &tra = species[al.num].tra_gl[al.num2];						
				if(tra.type == EXP_RATE) e = tra.dist_param[0].eq_ref;
			}
			break;
			
		default: break;
		}
		
		if(e != UNSET){
			const auto &eq = eqn[e];
			
			if(eq.linearise.on && eq.time_vari == true && eq.pop_ref.size() > 0){
				auto fl = false;
			
				for(const auto &pgc : eq.linearise.pop_grad_precalc){
					switch(pgc.type){
					case REG_PRECALC: case REG_PRECALC_TIME: if(affect_map.precalc_map[pgc.num]) fl = true; break;
					default: break;
					}
				}
				
				{
					const auto &fac = eq.linearise.factor_precalc;
					switch(fac.type){
					case REG_PRECALC: case REG_PRECALC_TIME: if(affect_map.precalc_map[fac.num]) fl = true; break;
					default: break;
					}
				}
				
				// Checks to see if update includes change to populations
				// (in which case div fast is not possible)
				if(fl == false){
					for(auto po : eq.pop_ref){
						auto p = pop[po].sp_p;
						
						for(auto ie : pop[po].ind_eff_mult){
							if(affect_map.ie_map[p][ie]) fl = true;
						}
			
						for(auto fe : pop[po].fix_eff_mult){
							if(affect_map.fe_map[p][fe]) fl = true;
						}
					}
				}
				
				if(fl == false){ 
					Listie ie; ie.i = i; ie.e = e;
					vec_list.push_back(ie);
				}
			}
		}
	}
	
	if(vec_list.size() == 0) return;
	
	vector <bool> remove(vec.size(),false);
	
	for(auto k = 0u; k < vec_list.size(); k++){
		auto i = vec_list[k].i;
		if(remove[i] == false){
			auto &al = vec[i];
			switch(al.type){
			case DIV_VALUE_AFFECT: al.type = DIV_VALUE_NOPOP_AFFECT; break;
			case MARKOV_POP_AFFECT: al.type = MARKOV_POP_NOPOP_AFFECT; break;
			default: emsg("no op"); break;
			}
			al.eq_nopop.list.push_back(al.num2);
			
			// Looks to combine multiple changes on to one update	
			const auto &eq = eqn[vec_list[k].e];
			const auto &it = eq.linearise.no_pop_precalc;
			
			for(auto kk = k+1; kk < vec_list.size(); kk++){
				auto ii = vec_list[kk].i;
				if(remove[ii] == false){
					const auto &al2 = vec[ii];
					auto eq2 = eqn[vec_list[kk].e];
					const auto &it2 = eq2.linearise.no_pop_precalc;
			
					if(it.type == it2.type && it.num == it2.num){
						al.eq_nopop.list.push_back(al2.num2);
						remove[ii] = true;
					}
				}
			}
		}
	}
	
	vector <AffectLike> vec_new;
	for(auto i = 0u; i < vec.size(); i++){
		if(remove[i] == false) vec_new.push_back(vec[i]);
	}
	vec = vec_new;
}


/// Looks to speed up calculation of DIV_VALUE_LINEAR_AFFECT/MARKOV_POP_LINEAR_AFFECT
// Makes changes to factor or no-pop value 
void Model::set_factor_nopop_only(vector <AffectLike> &vec, const vector <unsigned int> &param_list, const vector <unsigned int> &dependent, const SpecPrecalc &spec_precalc_after) const
{				
	auto affect_map = get_affect_map(vec,param_list,dependent,spec_precalc_after);
	
	for(auto i = 0u; i < vec.size(); i++){
		auto &ve = vec[i];
		if(ve.type == DIV_VALUE_LINEAR_AFFECT || ve.type == MARKOV_POP_LINEAR_AFFECT){	
			auto fl = false;
			for(const auto &lf : ve.lin_form.list){
				const auto &eq = eqn[lf.e];
				
				for(const auto &pgc : eq.linearise.pop_grad_precalc){
					switch(pgc.type){
					case REG_PRECALC: case REG_PRECALC_TIME: if(affect_map.precalc_map[pgc.num]) fl = true; break;
					default: break;
					}
				}
			
				// Checks to see if update includes change to populations
				// (in which case div fast is not possible)
				if(fl == false){
					for(auto po : eq.pop_ref){
						auto p = pop[po].sp_p;
						
						for(auto ie : pop[po].ind_eff_mult){
							if(affect_map.ie_map[p][ie]) fl = true;
						}
			
						for(auto fe : pop[po].fix_eff_mult){
							if(affect_map.fe_map[p][fe]) fl = true;
						}
					}
				}
				
				if(fl == true) break;
			}
			if(fl == false){
				ve.lin_form.factor_nopop_only = true;
			}
		}
	}
}

		
/// Looks to speed up calculation of DIV_VALUE_AFFECT 
// Useful when: (1) Time varying, (2) linearly seperatble into population terms
// Specifically it calculates gradients in population and uses
// changes in populations each time step to update the entire Markov div values 
void Model::affect_linearise_speedup(vector <AffectLike> &vec) const
{
	auto i = 0u; 
	while(i < vec.size()){
		while(i < vec.size()){
			if(vec[i].type == DIV_VALUE_AFFECT){
				if(species[vec[i].num].markov_eqn[vec[i].num2].time_vari == true) break;
			}
			if(vec[i].type == MARKOV_POP_AFFECT) break;
			i++;
		}		
					
		if(i < vec.size()){
			// Combines together all eqns linearisable and with the same time span and species
			auto p = vec[i].num;
			
			auto pop_ref_max = 0u;
			
			auto type = vec[i].type;
			
			auto imin = i;
			do{
				// Same species
				if(vec[i].num != p) break;
				
				// Must be of the same type
				if(vec[i].type != type) break;
			
			 	// Must be a rate
				unsigned int e;
				if(type == DIV_VALUE_AFFECT){
					const auto &meq = species[p].markov_eqn[vec[i].num2];	
					if(!meq.rate) break;  
					e = meq.eqn_ref;
				}
				else{
					const auto &tra = species[p].tra_gl[vec[i].num2];	
					if(tra.type != EXP_RATE) break;
					e = tra.dist_param[0].eq_ref;
				}
				
				auto &eq = eqn[e];
				const auto &lin = eq.linearise;
				
				// Must be linearisable
				if(lin.on == false) break;
				
				// Gradients without factor must be time independent 
				if(lin.pop_grad_time_dep == true) break;
				
				if(!equal_vec(vec[i].list,vec[imin].list)) break;
				
				auto pr = eq.pop_ref.size();
				if(pr > pop_ref_max) pop_ref_max = pr;
				
				i++;
			}while(i < vec.size());
			
			auto imax = i; 
			
			if(pop_ref_max > 0 && imax > imin){
				AffectLike af; 
				if(type == DIV_VALUE_AFFECT) af.type = DIV_VALUE_LINEAR_AFFECT;
				else af.type = MARKOV_POP_LINEAR_AFFECT;
				af.num = p; af.num2 = UNSET; af.list = vec[imin].list;
				af.map = vec[imin].map;
				
				const auto &sp = species[p];
				
				Hash hash_no_pop;
				
				vector <LinearFormInit> lfinit;
				for(auto j = imin; j < imax; j++){
					auto me = vec[j].num2;
					
					LinearFormInit lfi; 
					lfi.m = me;
					if(type == DIV_VALUE_AFFECT) lfi.e = sp.markov_eqn[me].eqn_ref;
					else lfi.e = sp.tra_gl[me].dist_param[0].eq_ref;
					lfinit.push_back(lfi);
				}
						
				species[p].set_linear_form(af.lin_form,lfinit,eqn);
				
				vec[imin] = af;
				
				if(imin+1 != imax){
					vec.erase(vec.begin()+imin+1,vec.begin()+imax);
					i = imin;
				}
			}
			i++;
		}
	}
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
				for(auto pme : pop[po].markov_eqn_ref){
					if(species[pme.p].markov_eqn[pme.e].infection_trans == true) fl = true;
				}
			}
		}
		
		if(vec[i].type == EXP_FE_AFFECT){
			auto p = vec[i].num;
			auto fe = vec[i].num2;
			
			for(auto po : species[p].fix_effect[fe].pop_ref){
				for(auto pme : pop[po].markov_eqn_ref){
					if(species[pme.p].markov_eqn[pme.e].infection_trans == true) fl = true;
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
			auto ie = vec[i].num2;
			
			for(auto po : species[p].ind_effect[ie].pop_ref){
				AffectLike al; 	
				al.type = POPNUM_IND_W_AFFECT; al.num = po; al.num2 = UNSET;
				param_vec_add_affect(vec,al);	
			}
		}
		
		if(vec[i].type == EXP_FE_AFFECT){
			auto p = vec[i].num;
			auto fe = vec[i].num2;
			
			for(auto po : species[p].fix_effect[fe].pop_ref){
				AffectLike al; 	
				al.type = POPNUM_IND_W_AFFECT; al.num = po; al.num2 = UNSET;
				param_vec_add_affect(vec,al);	
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
			
			AffectLike al; 			
			if(me.time_vari == true) al.map.resize(T,true);
			else al.map.resize(1,true);
					
			al.type = DIV_VALUE_AFFECT; al.num = p; al.num2 = e;
			param_vec_add_affect(vec,al);
		
			al.type = MARKOV_LIKE_AFFECT;
			param_vec_add_affect(vec,al);
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
				case OBS_COMP_EV: case OBS_TEST_EV: 
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
			AffectLike al; al.type = MARKOV_LIKE_AFFECT; al.num = p; al.num2 = e;
			auto &me = sp.markov_eqn[e];
			if(me.time_vari == true) al.map.resize(T,true);
			else al.map.resize(1,true);

			param_vec_add_affect(vec,al);	
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

	AffectLike al; 	
	al.type = POP_AFFECT; al.num = UNSET; al.num2 = UNSET; 
	al.map.resize(pop.size(),false); al.map[po] = true;
	param_vec_add_affect(vec,al);		
	
	for(const auto &mer : popu.markov_eqn_ref){
		AffectLike al; al.type = DIV_VALUE_AFFECT;
		al.num = mer.p; al.num2 = mer.e; 
		
		auto &me = species[mer.p].markov_eqn[mer.e];
		if(me.time_vari == true) al.map.resize(T,true);
		else al.map.resize(1,true);
		
		param_vec_add_affect(vec,al);		
			
		al.type = MARKOV_LIKE_AFFECT; 
		param_vec_add_affect(vec,al);	
	}
	
	for(const auto &tref : popu.trans_ref){
		AffectLike al; al.map.resize(T,true);
		al.type = MARKOV_POP_AFFECT; al.num = tref.p; al.num2 = tref.tr;
		param_vec_add_affect(vec,al);
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
	case MVN_JEF_PR: case MVN_UNIFORM_PR: case MVN_COR_PR:  // Sampling for this is done collectively
		return UNSET;
	
	case INVERSE_PR:
		{	
			auto min = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			auto max = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
			if(min <= 0) run_error("For the prior '"+pri.name+"' the minumum value '"+tstr(min)+"' has become non-positive.");
			
			if(max <= 0) run_error("For the prior '"+pri.name+"' the maximum value '"+tstr(max)+"' has become non-positive.");
			
			if(min >= max) run_error("For the prior '"+pri.name+"' the minumum value '"+tstr(min)+"' is larger than the maximum value '"+tstr(max)+"'.");
			//val = min*exp(ran()*log(max/min));
			
			val = min+ran()*(max-min);
		}
		break;
		
	case POWER_PR:
		{	
			auto min = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			auto max = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
			if(min <= 0) run_error("For the prior '"+pri.name+"' the minumum value '"+tstr(min)+"' has become non-positive.");
			
			if(max <= 0) run_error("For the prior '"+pri.name+"' the maximum value '"+tstr(max)+"' has become non-positive.");
			
			if(min >= max) run_error("For the prior '"+pri.name+"' the minumum value '"+tstr(min)+"' is larger than the maximum value '"+tstr(max)+"'.");
			
			//auto power = eqn[pri.dist_param[2].eq_ref].calculate_param(precalc);
			//if(power == -1) val = min*exp(ran()*log(max/min));
			//else val = pow(pow(min,power+1)+ ran()*(pow(max,power+1)-pow(min,power+1)),1.0/(power+1));
			val = min+ran()*(max-min);
		}
		break;
		
	case UNIFORM_PR:
		{	
			auto min = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			auto max = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
		
			if(min >= max) run_error("For the prior '"+pri.name+"' the minumum value '"+tstr(min)+"' is larger than the maximum value '"+tstr(max)+"'.");
			val = min+ran()*(max-min);
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
	
	if(val == UNSET){
		run_error("For prior '"+pri.name+"' the following error occurred: "+warn+". This prior is for '"+pri.in+"'. Prior distribution quantities (means, sds etc...) have threshold limits to ensure numerical accuracy. Consider changing the values or restricting the priors on model parameter which determine this distribution"); 
	}
	
	return val;
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
	
	if(par.variety == CONST_PARAM) return false;
	
	auto omega = sp.calculate_omega_basic(iegr.i,param_val,param);
	
	auto L = omega.size();
	for(auto j = 0u; j < L; j++){
		auto val = omega[j][j]; 
		if(val > ieg.var_max) return true;
	}
	
	auto det = determinant_fast(omega);

	if(det == UNSET) return true;
	if(det < ieg.log_det_min) return true;
	
	return false;
}

/*
/// Checks that the cholesky matrices can all be specified
bool Model::ie_cholesky_error(const PV &param_val) const
{
	for(auto p = 0u; p < species.size(); p++){
		const auto &sp = species[p];
		for(auto g = 0u; g < sp.ind_eff_group.size(); g++){
			auto omega = sp.calculate_omega_basic(g,param_val,param);
			auto illegal = false;
			calculate_cholesky(omega,illegal);
			if(illegal) return true;
		}
	}

	return false;
}
*/


/// Prints a set of parameters
void Model::print_param(const PV &param_val) const
{
	if(com_op == true) return;
	
	const auto &value = param_val.value;
	
	cout << "Parameters:" << endl;
	for(auto th = 0u; th < param_vec.size(); th++){
		cout << add_escape_char(param_vec[th].name) << " "<< value[th] << endl;
	}
	cout << endl;
	
	const auto &precalc = param_val.precalc;
	
	cout << "Precalc:" << endl;
	for(auto i = 0u; i < precalc.size(); i++){
		cout << precalc[i] << ",";
	}
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


/// Gets param_val_prop from a full parameter vector (i.e. removes reparam)
vector <double> Model::get_param_val_tvreparam(const PV &param_val) const
{
	const auto &value = param_val.value;
	
	vector <double> vec(nparam_vec_tvreparam);
	for(auto i = 0u; i < nparam_vec_tvreparam; i++){
		vec[i] = value[param_vec_tvreparam[i]];
	}
	
	return vec;
}


/// For the output of parametes tvreparam are placed into param_val
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

Param::Param(Constant &constant) : constant(constant) 
{
}

/// Gets the value for param vec
unsigned int Param::get_param_vec(unsigned int i) const
{
	const auto &er = element_ref[i];
	auto ind = er.index; 
	if(ind == UNSET) emsg("Cannot get element ref2");
	if(er.cons != false) emsg("Should not get constant1");
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
	if(!er.cons){
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

	// Extracts any precalculation within the parameter definitions
	param_vec_ref.resize(param_vec.size(),UNSET);
	spline_ref.resize(spline.size(),UNSET);
	
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
					precalc_eqn.add_eqn(eq.calcu,param_vec_ref,spline_ref,pv.spec_precalc_before);
					eq.precalc_done = true;
					//for(auto i : pv.list_precalc_before) list_precalc.push_back(i);
				}
			}
			break;
				
		case REPARAM_PARAM:
			{
				auto eq_ref = par.get_eq_ref(pv.index);
				if(eq_ref == UNSET) emsg("eq_ref should be set");
				auto &eq = eqn[eq_ref];
				
				precalc_eqn.add_eqn(eq.calcu,param_vec_ref,spline_ref,pv.spec_precalc_before);
				//for(auto i : pv.list_precalc_before) list_precalc.push_back(i);
			}
			break;
			
		default: break;
		}
		
		param_vec_ref[th] = precalc_eqn.add_param(th);
		
		auto s = last_spline[th];
		if(s != UNSET){
			spline_ref[s] = precalc_eqn.add_spline(s,spec_precalc);
		}
	}

	// Adds any splines not already added (these have no parameters)
	for(auto s = 0u; s < spline_ref.size(); s++){
		if(spline_ref[s] == UNSET){
			spline_ref[s] = precalc_eqn.add_spline(s,spec_precalc);
		}
	}
	
	// Extracts precalculation from all other equations apart from reparameterisation and distribution
	// Mulitpilier of populations are not combine to allow for pop gradient factor to be extracted
	for(auto &eq : eqn){
		if(eq.precalc_done == false && eq.type != DERIVE_EQN){
			eq.stop_combine_fl = precalc_eqn.add_eqn(eq.calcu,param_vec_ref,spline_ref,spec_precalc,PRECALC_STOP_COMBINE_MULT);
			eq.precalc_done = true;
		}
	}
}


/// Sets update in precalc at different times (this allows for populations in reparameterised eqns)
void Model::set_spec_precalc_time()
{
	auto T = details.T;
	
	SpecPrecalc spre_before;
	auto s_before = UNSET;
	
	spec_precalc_time.resize(T);
	for(auto th = 0u; th < param_vec.size(); th++){
		auto &pv = param_vec[th];
	
		if(pv.reparam_time_dep){
			auto s = pv.spline_ref;
			const auto &spl = spline[s];
			if(spl.type == SQUARE_SPL){
				if(s != s_before){
					s_before = s; 
					spre_before = pv.spec_precalc_before;
				}
				
				const auto &list_time = pv.spec_precalc_after.list_time;
				if(list_time.size() > 0){
					const auto &lt = list_time[0];
					if(lt.size() > 0){
						auto ti = lt[0];
						
						if(pv.reparam_spl_ti != ti) emsg("reparam_spl_ti problem");
						spec_precalc_time[ti].pv.push_back(th);
						
						pv.spec_precalc_before = precalc_eqn.shrink_sprec(lt,spre_before);
					}
				}
			}
		}
	}
	
	const auto S = spline.size();
	vector < vector <unsigned int> > spline_precalc;
	
	vector < vector <unsigned int> > spline_dep;
	auto C = precalc_eqn.calcu.size();
	spline_dep.resize(C);
	
	spline_precalc.resize(S);
	for(auto s = 0u; s < S; s++){
		const auto &spl = spline[s];
		if(spl.type == SQUARE_SPL){
			auto i = spline_ref[s];
			spline_precalc[s].push_back(i);
			spline_dep[i].push_back(s);
		}
	}
	
	for(const auto &sp : spec_precalc.info){
		auto i = sp.i;
		
		const auto &ca = precalc_eqn.calcu[i];
		
		if(ca.op != SINGLE){
			const auto &item = ca.item;
			const auto N = item.size();
			
			if(ca.time_dep){
				for(auto j = 0u; j < N; j++){
					const auto &it = item[j];
					
					if(it.type == REG_PRECALC_TIME){
						auto ii = it.num;
						for(auto s : spline_dep[ii]){
							spline_precalc[s].push_back(i);
							spline_dep[i].push_back(s);
						}
					}
				}
			}
		}
  }
	
	//vector <Hash> hash_ti(T);
	
	for(auto s = 0u; s < S; s++){
		const auto &list = spline_precalc[s];
		if(list.size() > 0){
			const auto &spl = spline[s];
			
			auto R = spl.param_ref.size();
			
			vector < vector <unsigned int> > ti_list;
			ti_list.resize(R);
			for(auto ti = 0u; ti < T; ti++){
				ti_list[spl.div[ti].index].push_back(ti);
			}
			
			if(false){
				for(auto r = 0u; r < R; r++){
					cout << r << ": ";
					for(auto ti : ti_list[r]) cout << ti << ",";
					cout << "r" << endl;
				}
			}
			
			for(auto r = 0u; r < R; r++){
				auto ti_first = ti_list[r][0];
				
				for(auto i : list){
					vector <unsigned int> ti_not_const;
					for(auto ti : ti_list[r]){
						if(precalc_init[i+ti] == UNSET) ti_not_const.push_back(ti);
					}
					
					if(ti_not_const.size() > 0){	
						precalc_eqn.sp_add(spec_precalc_time[ti_first].spec_precalc,i,ti_not_const);
					}
					/*
						auto &sp = spec_precalc_time[ti_first].spec_precalc;
						
						//auto k = hash_ti[ti_first].existing(ti_not_const);
						//if(k == UNSET){
							k = spec_precalc_time[ti_first].update.size();
							hash_ti[ti_first].add(k, ti_not_const);
							
							SpecPrecalc up;
							up.list_precalc_time = ti_not_const;
							
							spec_precalc_time[ti_first].update.push_back(up);
						}
						
						spec_precalc_time[ti_first].update[k].list_precalc.push_back(i);
					}
					*/
				}
			}
		}
	}

	// Removes those in spec_precalc_time from list_precalc
	auto M = precalc_eqn.calcu.size();
	vector <bool> map(M,false);
	
	for(const auto &spt : spec_precalc_time){
		for(const auto &in : spt.spec_precalc.info) map[in.i] = true;
	}
	
	vector <PrecalcInfo> info_new;
	for(const auto &in : spec_precalc.info){
		if(map[in.i] == false) info_new.push_back(in);
	}
	
	spec_precalc.info = info_new;
	
	if(false){ 
		for(auto ti = 0u; ti < T; ti++){ 
			auto &spt = spec_precalc_time[ti];
			if(spt.pv.size() > 0){
				cout << ti << ": ";
				for(auto th :  spt.pv) cout << param_vec[th].name << ",";
				cout << "vector update" << endl;
				
				//for(auto i : upt.pv_list_precalc) cout << i << ",";
				//cout << "pv_list_precalc" << endl;
				cout << endl;
			}
			
			const auto &spre = spt.spec_precalc;
			if(spre.info.size() > 0){
				print_spec_precalc("precalc",spre);
			}
		}
		//emsg("UU");
	}
}


/// Sets the initial value for precalc
void Model::set_precalc_init()
{        
	precalc_init = precalc_eqn.calculate_precalc_init(spec_precalc);
	
	if(false){
		for(auto i = 0u; i < precalc_init.size(); i++){
			cout << i << ": " << precalc_init[i] << endl;
		}
		cout << "precalc init" << endl;
		emsg("precalc init");
	}
}
	
	
/// Sets precalculation to be done after sampling 
void Model::set_spec_precalc_sample()
{
	spec_precalc_sample = precalc_eqn.calculate_spec_precalc_sample(spec_precalc);
}

	
/// Sets precalculation to be done after sampling 
void Model::set_spec_precalc_all()
{
	spec_precalc_all = precalc_eqn.calculate_spec_precalc_all(spec_precalc,spec_precalc_time);
}

	
/// Creates any precalulation for derived quantities
void Model::create_precalc_derive()
{	
	for(auto &eq : eqn){
		if(eq.precalc_done == false && eq.type == DERIVE_EQN){	
			precalc_eqn.add_eqn(eq.calcu,param_vec_ref,spline_ref,spec_precalc_derive);
		
			for(auto &inte : eq.integral){		
				precalc_eqn.add_eqn(inte.calc,param_vec_ref,spline_ref,spec_precalc_derive);
			}
			
			eq.precalc_done = true;
		}
	}
}


/// Shifts calculations in pop_grad into precalc	
void Model::create_precalc_pop_grad()
{
	for(auto &eq : eqn){
		if(eq.stop_combine_fl){
			precalc_eqn.add_eqn(eq.calcu,param_vec_ref,spline_ref,spec_precalc);
		}
	}
	
	for(auto &eq : eqn){
		if(eq.linearise.on){
			auto &lin = eq.linearise;
			
			for(auto j = 0u; j < lin.pop_grad_calc_store.size(); j++){
				precalc_eqn.add_eqn(lin.pop_grad_calc_store[j],param_vec_ref,spline_ref,spec_precalc);
			}
			precalc_eqn.add_eqn(lin.factor_calc,param_vec_ref,spline_ref,spec_precalc);
			precalc_eqn.add_eqn(lin.no_pop_calc_store,param_vec_ref,spline_ref,spec_precalc);
		}
	}
}

	

/// Works out how precalculation is affected by changes in parameters
void Model::precalc_affect()
{
	const auto &calcu = precalc_eqn.calcu;
	auto C = calcu.size();

	// Works out which affect markov equation
	vector < vector <AffectME> > affect_me;
	affect_me.resize(C);
	
	vector < vector <bool> > map_me;
	map_me.resize(species.size());
	for(auto p = 0u; p < species.size(); p++){
		const auto &sp = species[p];
		switch(sp.type){
		case INDIVIDUAL:
			{
				map_me[p].resize(sp.markov_eqn.size(),false);
				for(auto i = 0u; i < sp.markov_eqn.size(); i++){
					const auto &me = sp.markov_eqn[i];
						
					const auto &eq = eqn[me.eqn_ref];
					for(const auto &ca : eq.calcu){
						for(const auto &it : ca.item){
							switch(it.type){
								case REG_PRECALC: case REG_PRECALC_TIME:
									AffectME ame; ame.p = p; ame.e = i;
									affect_me[it.num].push_back(ame);
									break;
								
								default:
									break;
							}
						}					
					}
				}
			}
			break;
			
		case POPULATION:
			{
				map_me[p].resize(sp.tra_gl.size(),false);
				for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
					const auto &tra = sp.tra_gl[tr];
					const auto &eq = eqn[tra.dist_param[0].eq_ref];
				
					for(const auto &ca : eq.calcu){
						for(const auto &it : ca.item){
							switch(it.type){
								case REG_PRECALC: case REG_PRECALC_TIME:
									AffectME ame; ame.p = p; ame.e = tr;
									affect_me[it.num].push_back(ame);
									break;
								
								default:
									break;
							}
						}					
					}
				}
			}
			break;
		}		
	}
	
	// Creates a list of which elements are subsequently affected
	vector < vector <unsigned int> > affect;
	affect.resize(C);
	
	for(auto i = 0u; i < C; i++){
		const auto &ca = calcu[i];
		for(const auto &it : ca.item){
			switch(it.type){
			case PARAMVEC:
				{
					const auto &pv = param_vec[it.num];
					const auto &par = param[pv.th];
					
					if(par.variety == REPARAM_PARAM){
						auto eq_ref = par.get_eq_ref(pv.index);
						if(eq_ref == UNSET) emsg("eq_ref should be set");
						for(const auto &ca2 : eqn[eq_ref].calcu){
							for(const auto &it2 : ca2.item){
								switch(it2.type){
								case REG_PRECALC: case REG_PRECALC_TIME:
									affect[it2.num].push_back(i);	
									break;
								default: break;
								}
							}
						}
					}
				}
				break;
				
			case REG_PRECALC: case REG_PRECALC_TIME:
				affect[it.num].push_back(i);	
				break;
			
			default: break;
			}
		}			
	}
	
	auto T = details.T;
	
	vector <bool> map_PC(C,false);
		
	vector <bool> map_time(T,true);
	for(auto j = 0u; j < param_vec.size(); j++){
		// Spline dependency
		const auto &pv = param_vec[j];
		
		auto th = pv.th, ind = pv.index;
	
		const auto &par = param[th];
		
		if(par.spline_info.on == true){	
			auto nknot = par.spline_info.knot_tdiv.size();
		
			vector <unsigned int> vec; vec.push_back(th); vec.push_back(ind/nknot); 
			auto s = hash_spline.existing(vec);
			if(s == UNSET) emsg_input("Cannot find spline");
		
			{ // turn off
				auto ss = 0u;
				while(ss < spline.size() && !(th == spline[ss].th && ind/nknot == spline[ss].index)) ss++;
				if(ss != s) emsg("not agree");
			}
			
			const auto &spl = spline[s];
			
			vector <bool> map(T,false);
			switch(spl.type){
			case LINEAR_SPL: 
				for(auto ti = 0u; ti < T; ti++){ 
					auto ind = spl.div[ti].index;
					const auto &pr = spl.param_ref[ind];
					if(!pr.cons && pr.index == j) map[ti] = true;
					
					if(ind+1 >= spl.param_ref.size()) emsg("Out of range");
					const auto &pr2 = spl.param_ref[ind+1];
					if(!pr2.cons && pr2.index == j) map[ti] = true;
				}
				break;
				
			case SQUARE_SPL:
				for(auto ti = 0u; ti < T; ti++){ 
					auto ind = spl.div[ti].index;
					const auto &pr = spl.param_ref[ind];
					if(!pr.cons && pr.index == j) map[ti] = true;
				}
				break;
				
			case CUBICPOS_SPL: case CUBIC_SPL:
				for(auto ti = 0u; ti < T; ti++) map[ti] = true;
				break;
			}
	
			add_affect_like(spline_ref[s],param_vec_ref[j],map,map_PC,map_me,param_vec[j],affect,affect_me);
		}		
		else{
			// Basic parameter dependency
			add_affect_like(param_vec_ref[j],UNSET,map_time,map_PC,map_me,param_vec[j],affect,affect_me);
		}
	}

	// Makes sure list_precalc does not have repeated value
	vector <bool> mapl(C,false);
	
	auto info_st = spec_precalc.info;
	spec_precalc.info.clear();
	for(const auto &in : info_st){
		if(mapl[in.i] == false){
			mapl[in.i] = true;
			spec_precalc.info.push_back(in);
		}
	}		
	
	if(testing){  // Checks that there is no overlap between  
		vector <bool> mapl(C,false);
		for(const auto &pv : param_vec){
			for(const auto &in : pv.spec_precalc_before.info){
				if(mapl[in.i] == true) emsg("double1");
				mapl[in.i] = true;
			}
			
			/*
			for(const auto &in : pv.set_param_spec_precalc.info){
				if(mapl[in.i] == true) emsg("double2");
				mapl[in.i] = true;
			}
			*/
			
			for(const auto &in : pv.spec_precalc_after.info){
				if(mapl[in.i] == true) emsg("double3");
				mapl[in.i] = true;
			}
			
			for(const auto &in : pv.spec_precalc_before.info) mapl[in.i] = false;
			for(const auto &in : pv.set_param_spec_precalc.info) mapl[in.i] = false;
			for(const auto &in : pv.spec_precalc_after.info) mapl[in.i] = false;
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
	
	if(false) print_precalc();
}


/// Prints information about precalculation
void Model::print_precalc() const
{		
	precalc_eqn.print_calc(); 
	
	cout << "PARAM PRECALC" << endl;
	auto thmax = param_vec.size();
	if(thmax > 100) thmax = 100;
	for(auto th = 0u; th < thmax; th++){
		const auto &pv = param_vec[th];
	
		cout << "<<" << pv.name << ">>" << endl;
		print_spec_precalc("before",pv.spec_precalc_before);
		print_spec_precalc("set param",pv.set_param_spec_precalc);
		print_spec_precalc("after",pv.spec_precalc_after);
		cout << endl;
	} 
	if(thmax == 100) cout << "..." << endl;
	
	//print_spec_precalc("spec_precalc",spec_precalc);
		
	print_spec_precalc("HERE spec_precalc_sample",spec_precalc_sample);
		
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


/// Adds how parameter change results in change in likelihood
void Model::add_affect_like(unsigned int i, unsigned int i2, const vector <bool> &map_time, vector <bool> &map_PC, vector < vector <bool> > &map_me, ParamVecEle &pvec, const vector < vector <unsigned int> > &affect, const vector < vector <AffectME> > &affect_me)
{
	vector <unsigned int> list;
	if(i2 != UNSET) list.push_back(i2);
	list.push_back(i);
	auto k = 0u;
	while(k < list.size()){
		for(auto ii : affect[list[k]]){
			if(map_PC[ii] == false){
				map_PC[ii] = true;
				list.push_back(ii);
			}
		}
		k++;
	}
	
	sort(list.begin(),list.end());

	vector <unsigned int> list_time;
	for(auto ti = 0u; ti < map_time.size(); ti++){
		if(map_time[ti] == true) list_time.push_back(ti);
	}

	precalc_eqn.sp_add(pvec.spec_precalc_after,list,list_time);
	
	for(auto i : list) map_PC[i] = false;
	
	vector <AffectME> list_me;
	for(auto i : list){
		for(const auto &me : affect_me[i]){
			if(map_me[me.p][me.e] == false){
				map_me[me.p][me.e] = true;
				list_me.push_back(me);
			}
		}
	}
	
	for(const auto &me_ref : list_me){
		auto p = me_ref.p, e = me_ref.e;
	
		const auto &sp = species[p];
		switch(sp.type){
		case INDIVIDUAL:
			{
				const auto &me = sp.markov_eqn[e];
					
				AffectLike al; 
				if(me.time_vari == true) al.map = map_time;
				else al.map.resize(1,true);

				al.type = DIV_VALUE_AFFECT; al.num = p; al.num2 = e;
				param_vec_add_affect(pvec.affect_like,al);
			
				al.type = MARKOV_LIKE_AFFECT;
				param_vec_add_affect(pvec.affect_like,al);
				
				map_me[p][e] = false;
			}
			break;
			
		case POPULATION:
			{
				AffectLike al; al.map = map_time;
				al.type = MARKOV_POP_AFFECT; al.num = p; al.num2 = e;
				param_vec_add_affect(pvec.affect_like,al);
				map_me[p][e] = false;
			}
			break;
		}
	}
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
	if(is_matrix(par) && begin_str(par.name,"Ω")) return true;
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
		ss << in.i;
		if(in.tlist != UNSET){
			const auto &lt = spre.list_time[in.tlist];
			ss << "(" << str_time_range(lt) << ")";
		
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
	vector <unsigned int> no_time;
	
	for(auto th = 0u; th < param_vec.size(); th++){
		auto &pv = param_vec[th];
		auto spl_fl = false; if(pv.spline_ref != UNSET) spl_fl = true;
		precalc_eqn.set_param(pv.set_param_spec_precalc,pv.spec_precalc_after,spl_fl);
		//precalc_eqn.sp_add(pv.set_param_spec_precalc,param_vec_ref[th],no_time);
	}
}


/// If model contains tvreparam then create a reference from pop to th
void Model::set_pop_reparam_th()
{
	if(!contains_tvreparam) return;
	
	auto P = pop.size();
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
						
					case POPTIMENUM:
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

	if(false){
		for(auto po = 0u; po < P; po++){
			for(auto ti = 0u; ti < T; ti++){
				for(auto th : pop_reparam_th[po][ti]){
					cout << ti << " " << pop[po].name << " " << param_vec[th].name << endl;
				}
			}
		}
		emsg("pop reparam th");
	}
}

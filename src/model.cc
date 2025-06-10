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
Model::Model(Operation mode_)
{
	mode = mode_;  
	nspecies = 0;
};


/// Adds an equation reference to EquationInfo
void Model::add_eq_ref(EquationInfo &eqi, Hash &hash_eqn, double t)
{
	auto ti = UNSET;
	if(t != UNSET){
		ti = 0; while(ti < timepoint.size()-2 && timepoint[ti] < t) ti++;
	}	
	
	if(eqi.te == "") emsg("Equation does not have any text");
	
	// If "max" is used in sums then this means the equation depends the state of eqi.c
	string c_te = "";
	if(eqi.c != UNSET && find_in(eqi.te,"max") != UNSET) c_te = find_eq_dif(eqi);
	
	auto vec = hash_eqn.get_vec_eqn(eqi.te,(unsigned int)eqi.type,eqi.p,eqi.cl,eqi.infection_trans,ti,c_te);
	
	auto e = hash_eqn.existing(vec);
	if(e != UNSET){
		eqi.eq_ref = e;
		/*
		if(eqi.te != eqn[e].te_init){
			cout << eqi.te << "| " << eqn[e].te_init <<" init" << endl;
			emsg("wrong");
		}
		*/
	}
	else{
		eqi.eq_ref = eqn.size();	
		
		hash_eqn.add(eqi.eq_ref,vec);
		Equation eq(eqi.te,eqi.type,eqi.p,eqi.cl,eqi.c,eqi.infection_trans,ti,eqi.line_num,species_simp,param,spline,param_vec,pop,hash_pop,timepoint);
		
		eqn.push_back(eq);
	}
}


/// When an equation contains "max" this works out an equation differentiator
string Model::find_eq_dif(const EquationInfo &eqi) const 
{
	const auto &sp = species[eqi.p];
	
	vector <bool> cl_fl(sp.ncla,false);
	
	auto te = eqi.te;
	auto i = 0u;
	auto len = te.length();
	while(i < len){
		while(i < len-4 && te.substr(i,4) != "max:") i++;
		if(i < len-4){
			while(i < len && te.substr(i,1) != "_") i++;
			if(i < len){
				i++;
				auto di = get_dep_info(te,i,notparam_list);
				for(auto ve : di.spl){
					auto ind = remove_prime(ve);
					for(auto cl = 0u; cl < sp.ncla; cl++){
						if(sp.cla[cl].index == ind) cl_fl[cl] = true;
					}
				}
				i = di.iend;
			}				
		}
		i++;
	}
	
	const auto &cgl = sp.comp_gl[eqi.c];
	
	string c_te = "";
	for(auto cl = 0u; cl < sp.ncla; cl++){
		if(cl_fl[cl] == true){
			if(c_te != "") c_te += "_";
			c_te += sp.cla[cl].comp[cgl.cla_comp[cl]].name;
		}
	}
	
	return c_te;
}


/// Samples a set of parameter values from the model
vector <double> Model::param_sample() const
{
	vector <double> param_val(nparam_vec,0);
	
	auto loop = 0u, loopmax = 100u;
	do{
		for(auto th = 0u; th < nparam_vec; th++){
			const auto &pv = param_vec[th];
			
			const auto &par = param[pv.th];
			switch(par.variety){
			case CONST_PARAM:
				param_val[th] = par.get_value(pv.index);
				break;
				
			case REPARAM_PARAM: 
				{
					auto eq_ref = par.get_eq_ref(pv.index);
					if(eq_ref == UNSET) emsg("eq_ref should be set");
					param_val[th] = eqn[eq_ref].calculate_param_only(param_val);
				}
				break;
				
			case DIST_PARAM: case PRIOR_PARAM:
				param_val[th] = prior_sample(prior[pv.prior_ref],param_val);
				break;
				
			case UNSET_PARAM: emsg("error param"); break;
			}
		}
	
		if(!ie_cholesky_error(param_val)) break;
		loop++;
	}while(loop < loopmax);
	if(loop == loopmax) emsg("Could not sample parameters due to Cholesky error");
	
	return param_val;
}


/// Gather a parameter set from a posterior sample
vector <double> Model::post_param(const Sample &samp) const
{
	vector <double> param_val(nparam_vec);
	
	for(auto th = 0u; th < nparam_vec; th++){
		const auto &pv = param_vec[th];	
		const auto &par = param[pv.th];
		switch(par.variety){
		case CONST_PARAM:
			param_val[th] = par.get_value(pv.index);
			break;
			
		case REPARAM_PARAM:
			{
				auto eq_ref = par.get_eq_ref(pv.index);
				if(eq_ref == UNSET) emsg("eq_ref should be set");
				param_val[th] = eqn[eq_ref].calculate_param_only(param_val);
			}
			break;
			
		case DIST_PARAM: case PRIOR_PARAM:	
			if(pv.ppc_resample){
				param_val[th] = prior_sample(prior[pv.prior_ref],param_val);
			}
			else{
				param_val[th] = samp.param_value[pv.th][pv.index];
				if(param_val[th] == UNSET) emsg("Parameter is unspecified");
			}
			break;
			
		case UNSET_PARAM: emsg("error param"); break;
		}
	}
	
	return param_val;
}


/// Determines if a set of parameters is within the bounds of the prior
bool Model::inbounds(const vector <double> &param_val) const
{
	auto total = prior_total(param_val) + dist_total(param_val);
	
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
	vector <double> pri(nparam_vec,0);
	
	for(auto th = 0u; th < nparam_vec; th++){
		const auto &pv = param_vec[th];
		if(pv.variety == PRIOR_PARAM){
			pri[th] = prior_probability(param_val[th],prior[pv.prior_ref],param_val,eqn);
		}
	}
	
	return pri;
}


/// The total prior for all parameters
double Model::dist_total(const vector <double> &param_val) const 
{
	auto dist = dist_prob(param_val);
	return sum(dist);
}


/// Calculate the prior for the parameters
vector <double> Model::dist_prob(const vector <double> &param_val) const
{
	vector <double> dist(nparam_vec,0);
	
	for(auto th = 0u; th < nparam_vec; th++){
		const auto &pv = param_vec[th];
		if(pv.variety == DIST_PARAM){
			dist[th] = prior_probability(param_val[th],prior[pv.prior_ref],param_val,eqn);
		}
	}
	
	return dist;
}


/// Recalculates the prior for a given parameter
double Model::recalculate_prior(unsigned int th, vector <double> &prior_prob, const vector <double> &param_val, double &like_ch) const
{
	auto store = prior_prob[th];
	const auto &pv = param_vec[th];
	prior_prob[th] = prior_probability(param_val[th],prior[pv.prior_ref],param_val,eqn);
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
			
			const auto &knot = spl.info.knot_times;
	
			for(auto t = 1u; t < spl.param_ref.size(); t++){
				auto val = param_val[spl.param_ref[t]];
				auto last = param_val[spl.param_ref[t-1]];
				
				auto dt = knot[t]-knot[t-1];
				if(dt > details.dt){  // Only applies spline prior if time difference greater than dt	
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
double Model::recalculate_spline_prior(unsigned int s, vector <double> &spline_prior, const vector <double> &param_val, double &like_ch) const
{
	const auto &spl = spline[s];
	const auto &info = spl.info;
	const auto &knot = info.knot_times;
		
	auto type = info.smooth_type;
	auto sd = info.smooth_value;
	auto cv = sqrt(exp(sd*sd)-1);
	
	auto Li = 0.0;
	for(auto t = 1u; t < spl.param_ref.size(); t++){
		auto val = param_val[spl.param_ref[t]];
		auto last = param_val[spl.param_ref[t-1]];
		
		auto dt = knot[t]-knot[t-1];
		if(dt > details.dt){  // Only applies spline prior if time difference greater than dt	
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
	

/// Calculates pop from cpop
vector <double> Model::calculate_popnum(vector <StateSpecies> &state_species) const
{
	vector <double> popnum(pop.size(),0);
	
	for(auto i = 0u; i < pop.size(); i++){            // Calculates population based on cpop
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
			for(auto i = 0u; i < ssp.individual.size(); i++){
				const auto &ind = ssp.individual[i];
				
				auto c = ssp.ind_sim_c[i];
				if(c != UNSET){
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
// Used only for diagnostics
vector < vector <double> > Model::calculate_popnum_t(vector <StateSpecies> &state_species, unsigned int ti_end) const
{
	vector < vector <double> > popnum_t;

	if(ti_end == UNSET) ti_end = ntimepoint-1;
	for(auto ti = 0u; ti < ti_end; ti++){
		for(auto p = 0u; p < nspecies; p++){
			auto &ssp = state_species[p];
			switch(ssp.type){
			case INDIVIDUAL: ssp.set_ind_sim_c(ti); break;
			case POPULATION: ssp.cpop = ssp.cpop_st[ti]; break;
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
		
				auto c = UNSET;
				for(auto e = 0u; e <= ind.ev.size(); e++){
					double t;
					if(e < ind.ev.size()) t = ind.ev[e].t;
					else t = timepoint[T];
						
					while(timepoint[ti] < t){
						if(c != UNSET){						
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
		SpeciesSimp ss(sp.name,sp.cla,sp.ind_effect,sp.fix_effect,sp.comp_mult,sp.comp_gl,sp.trans_tree);
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
		case SPLINE_PRIOR_AFFECT: case PRIOR_AFFECT: case DIST_AFFECT: case OMEGA_AFFECT: 
		case LIKE_INIT_COND_AFFECT: case PRIOR_INIT_COND_AFFECT:
		case GENETIC_VALUE_AFFECT:
			al.order_num = 1; 
			break;
			
		case SPLINE_AFFECT: case EXP_FE_AFFECT: case EXP_IE_AFFECT:
			al.order_num = 2; 
			break;
			
		case POP_AFFECT: 
			al.order_num = 3; 
			break;
			
		case INDFAC_INT_AFFECT: 
			al.order_num = 4; 
			break;
		
		case DIV_VALUE_AFFECT: case DIV_VALUE_FAST_AFFECT: case DIV_VALUE_LINEAR_AFFECT: 
		case POP_DATA_CGL_TGL_AFFECT:
			al.order_num = 5; 
			break;
			
		case OBS_EQN_AFFECT: //case OBS_TRANS_EQN_AFFECT:
		case IIF_W_AFFECT:
		case POPNUM_IND_W_AFFECT:
			al.order_num = 6; 
			break;
			
		case MARKOV_LIKE_AFFECT: case NM_TRANS_AFFECT:  case NM_TRANS_BP_AFFECT: case NM_TRANS_INCOMP_AFFECT: case LIKE_IE_AFFECT: 
		case MARKOV_POP_AFFECT: case LIKE_OBS_IND_AFFECT: case LIKE_OBS_POP_AFFECT: 
		case LIKE_OBS_POP_TRANS_AFFECT: case LIKE_UNOBS_TRANS_AFFECT:
		case LIKE_GENETIC_PROCESS_AFFECT: case LIKE_GENETIC_OBS_AFFECT:
			al.order_num = 7; 
			break;
			
		case AFFECT_MAX: break;
		}
	}
	
	sort(vec.begin(),vec.end(),AL_ord);
}


/// Determines if possible to convert DIV_VALUE_AFFECT -> DIV_VALUE_FAST_AFFECT
bool Model::div_value_fast_possible(const AffectLike &al) const
{
	if(al.type == DIV_VALUE_AFFECT){
		const auto &meq = species[al.num].markov_eqn[al.num2];
		auto eq = eqn[meq.eqn_ref];
		if(eq.linearise.on && eq.time_vari == true && eq.pop_ref.size() > 0 && meq.rate) return true;
	}
	return false;
}


/// Looks to speed up calculation of DIV_VALUE_AFFECT for parameter proposals
// Useful when then thing changing can be seperated from population terms
// Specifically it calculates changes in non-population term as a way to speed up div value calculation 
void Model::affect_linearise_speedup2(vector <AffectLike> &vec, const vector <unsigned int> &param_list, const vector <unsigned int> &dependent) const
{			
	auto N = param_vec.size();
	
	// Creates a map which shows which parameters change under proposal
	vector <bool> param_map(N,false);
	for(auto k : param_list) param_map[k] = true;
	for(auto k : dependent) param_map[k] = true;
		
	// Creates maps for fixed and individual effects which change	
	vector < vector <bool> > ie_map;
	ie_map.resize(nspecies);
	for(auto p = 0u; p < nspecies; p++){
		const auto &sp = species[p];
		ie_map[p].resize(sp.ind_effect.size(),false);
	}
	
	vector < vector <bool> > fe_map;
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
	
	// Generates a list of elements in vec which can be converted	
	vector <unsigned int>	vec_list;
	
	for(auto i = 0u; i < vec.size(); i++){
		auto &ve = vec[i];
		if(div_value_fast_possible(ve)){	
			const auto &meq = species[ve.num].markov_eqn[ve.num2];
			auto eq = eqn[meq.eqn_ref];
	
			auto fl = false;
		
			const auto &pop_grad_calc = eq.linearise.pop_grad_calc;
			for(auto j = 0u; j < pop_grad_calc.size(); j++){
				for(const auto &ca : pop_grad_calc[j]){
					for(const auto &it : ca.item){ 
						switch(it.type){
						case PARAMVEC:
							if(param_map[it.num] == true){ fl = true; break;}
							break;
							
						case SPLINEREF:
							{
								const auto &spl = spline[it.num];
								for(auto th : spl.param_ref){
									if(param_map[th] == true){ fl = true; break;}
								}
							}
							break;
						
						case PARAMETER: case SPLINE: emsg("SHould not be"); break;
						
						default: break;
						}
					}
					if(fl == true) break;
				}
				if(fl == true) break;
			}
			
			// Checks to see if update includes change to populations
			// (in which case div fast is not possible)
			if(fl == false){
				for(auto po : eq.pop_ref){
					auto p = pop[po].sp_p;
					
					for(auto ie : pop[po].ind_eff_mult){
						if(ie_map[p][ie]) fl = true;
					}
		
					for(auto fe : pop[po].fix_eff_mult){
						if(fe_map[p][fe]) fl = true;
					}
				}
			}
			
			if(fl == false) vec_list.push_back(i);
		}
	}
	
	if(vec_list.size() == 0) return;
	
	vector <bool> remove(vec.size(),false);
	
	for(auto k = 0u; k < vec_list.size(); k++){
		auto i = vec_list[k];
		if(remove[i] == false){
			auto &ve = vec[i];
			ve.type = DIV_VALUE_FAST_AFFECT;
			ve.me_list.push_back(ve.num2);
			
			// Looks to combine multiple changes on to one update	
			const auto &meq = species[ve.num].markov_eqn[ve.num2];
			auto eq = eqn[meq.eqn_ref];
			for(auto kk = k+1; kk < vec_list.size(); kk++){
				auto ii = vec_list[kk];
				if(remove[ii] == false){
					auto &ve2 = vec[ii];
					
					const auto &meq2 = species[ve2.num].markov_eqn[ve2.num2];
					auto eq2 = eqn[meq2.eqn_ref];
					if(eq.equal_calc(eq.linearise.no_pop_calc,eq.cons,eq2.linearise.no_pop_calc,eq2.cons)){
						ve.me_list.push_back(ve2.num2);
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

					
/// Looks to speed up calculation of DIV_VALUE_AFFECT 
// Useful when: (1) Time varying, (2) linearly seperatble into population terms
// Specifically it calculates gradients in population and uses
// changes in populations each time step to update the entire Markov div values 
void Model::affect_linearise_speedup(vector <AffectLike> &vec) const
{
	auto i = 0u; 
	while(i < vec.size()){
		while(i < vec.size() && !(vec[i].type == DIV_VALUE_AFFECT && species[vec[i].num].markov_eqn[vec[i].num2].time_vari == true)) i++;
		if(i < vec.size()){
			// Combines together all eqns linearisable and with the same time span and species
			auto p = vec[i].num;
			
			auto imin = i;
			do{
				if(vec[i].num != p) break;
				
				const auto &meq = species[p].markov_eqn[vec[i].num2];
				if(vec[i].type != DIV_VALUE_AFFECT) break;
				if(!meq.rate) break;
				
				auto eq = eqn[meq.eqn_ref];
				if(eq.linearise.on == false) break;
				
				if(i == imin && eq.pop_ref.size() == 0) break;
				
				if(eq.linearise.pop_grad_calc_time_dep == true) break;
				if(!equal_vec(vec[i].list,vec[imin].list)) break;
				i++;
			}while(i < vec.size());
			auto imax = i; 
			
			if(imax > imin){
				AffectLike af; 
				af.type = DIV_VALUE_LINEAR_AFFECT;
				af.num = p; af.num2 = UNSET; af.list = vec[imin].list;
				af.map = vec[imin].map;
				auto &lin = af.linear_prop;
				
				const auto &sp = species[p];
				
				Hash hash_no_pop;
				
				for(auto j = imin; j < imax; j++){
					auto me = vec[j].num2;
					const auto &eq = eqn[sp.markov_eqn[me].eqn_ref];
				
					for(auto i = 0u; i < eq.pop_ref.size(); i++){
						auto po = eq.pop_ref[i];
						
						auto k = 0u; while(k < lin.pop_affect.size() && lin.pop_affect[k].po != po) k++;
						if(k == lin.pop_affect.size()){
							PopAffect pa; pa.po = po;
							lin.pop_affect.push_back(pa);
						}
						
						GradRef gr; gr.ref = lin.me.size(); gr.index = i;
						lin.pop_affect[k].pop_grad_ref.push_back(gr);
					}
					
					// Adds a no_pop calculation
					if(eq.linearise.no_pop_calc_time_dep){ 
						auto hnum = eq.get_calc_hash_num(eq.linearise.no_pop_calc);
						auto vec = hash_no_pop.get_vec_double(hnum);
						
						auto fl = false;
						
						auto n = hash_no_pop.existing(vec);
						if(n != UNSET){
							auto me_sa = lin.me_no_pop[n];
							const auto &eq_sa = eqn[sp.markov_eqn[me_sa].eqn_ref];
				
							if(eq.equal_calc(eq.linearise.no_pop_calc,eq.cons,eq_sa.linearise.no_pop_calc,eq_sa.cons)){
								lin.no_pop_ref.push_back(n);
								fl = true;
							}
						}
						
						if(fl == false){ // Adds a new equation if not already present
							auto k = lin.me_no_pop.size();
							hash_no_pop.add(k,vec);
							lin.me_no_pop.push_back(me);
						
							lin.no_pop_ref.push_back(k);
						}
					}
					else lin.no_pop_ref.push_back(UNSET);
					
					lin.me.push_back(me);
				}
				
				if(false){
					for(auto i = 0u; i < lin.me.size(); i++){
						const auto &eq = eqn[sp.markov_eqn[lin.me[i]].eqn_ref];
						cout << i << " " << lin.me[i] << " " << eq.te_raw <<"  ME" << endl;
					}
					
					for(auto &pa : lin.pop_affect){
						cout << pa.po << " " << pop[pa.po].name << " ";
						for(auto &gr : pa.pop_grad_ref) cout << gr.ref << "," << gr.index <<"   ";
						cout << endl;
					}
					emsg("affect");
				}
					
				vec[imin] = af;
				
				if(imin+1 != imax){
					vec.erase(vec.begin()+imin+1,vec.begin()+imax);
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


/// For joint proposals works out the tihings which need updating
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
	
	default: emsg("SHould not be op"); break;
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
	
	auto T = ntimepoint-1;
	
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
	auto T = ntimepoint-1;
	
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
	auto T = ntimepoint-1;

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
vector <InitCondValue> Model::initc_sample(const vector <double> &param_val) const
{
	vector <InitCondValue> initc(nspecies);
	
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
					inc.N_total = round_int(prior_sample(ic.pop_prior,param_val));

					inc.frac = dirichlet_sample(ic.alpha);
					inc.cnum = multinomial_sample(inc.N_total,inc.frac);
				}
				else{
					const auto &claa = sp.cla[foc_cl];
					inc.N_focal.resize(claa.ncomp);
					for(auto c = 0u; c < claa.ncomp; c++){
						inc.N_focal[c] = round_int(prior_sample(ic.comp_prior[c],param_val));
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
		
	case DIRICHLET_PR:
		{
			auto alpha = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			return gamma_alpha_sample(alpha);
		}
		break;
		
	case MDIR_PR:	
		emsg("SHould not sample from MDIR");
		break;
	
	case UNSET_PR:
		emsg("Prior should be set");
		break;
	}
	
	return UNSET;
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
	//auto vec = hash_all_ind.get_vec_string(name);
	//auto j = hash_all_ind.existing(vec);	
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
				//auto vec = hash_all_ind.get_vec_string(sp.individual[i].name);
				//hash_all_ind.add(all_ind.size(),vec);
				hash_all_ind.add(all_ind.size(),sp.individual[i].name);
				AllInd ai; ai.i = i; ai.p = p;
				all_ind.push_back(ai);
			}
		}
	}
}


/// Checks that the cholesky matrices can all be specified
bool Model::ie_cholesky_error(const vector <double> &param_val) const
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


/// Prints a set of parameters
void Model::print_param(const vector <double> &vec) const
{
	if(com_op == true) return;
	
	cout << "Parameters:" << endl;
	for(auto th = 0u; th < param_vec.size(); th++){
		cout << param_vec[th].name << " "<< vec[th] << endl;
	}
	cout << endl;
}


/// Gets param_val_prop from a full parameter vector (i.e. removes reparam)
vector <double> Model::get_param_val_prop(const vector <double> &param_val) const
{
	vector <double> vec(nparam_vec_prop);
	for(auto i = 0u; i < nparam_vec_prop; i++){
		vec[i] = param_val[param_vec_prop[i]];
	}
	
	return vec;
}


/// Reconstructs param_val from param_val_prop
vector <double> Model::get_param_val(const vector <double> &param_val_prop) const
{
	vector <double> param_val(nparam_vec,UNSET);
	for(auto i = 0u; i < nparam_vec_prop; i++){
		param_val[param_vec_prop[i]] = param_val_prop[i];
	}
	
	for(auto th = 0u; th < nparam_vec; th++){
		const auto &pv = param_vec[th];
		
		if(param_val[th] == UNSET){
			const auto &par = param[pv.th];
			switch(par.variety){
			case PRIOR_PARAM:
				{
					const auto &pri = prior[par.get_prior_ref(pv.index)];
					if(pri.type != FIX_PR) emsg("Prior should be fixed");
					param_val[th] = prior_sample(pri,param_val);
				}
				break;
				
			case REPARAM_PARAM:
				{
					auto eq_ref = par.get_eq_ref(pv.index);
					if(eq_ref == UNSET) emsg("eq_ref should be set");
					param_val[th] = eqn[eq_ref].calculate_param_only(param_val);
				}
				break;
			
			default: emsg("Option prob"); break;
			}
		}
	}
	
	return param_val;
}


/// Gets the value for param vec
unsigned int Param::get_param_vec(unsigned int i) const
{
	auto ref = element_ref[i]; if(ref == UNSET) emsg("Cannot get element ref1");
	return element[ref].param_vec_ref;
}


/// Gets a vector of children from an element
const vector <ParamRef>& Param::get_child(unsigned int i) const
{
	auto ref = element_ref[i]; if(ref == UNSET) emsg("Cannot get element ref2");
	return element[ref].child;
}


/// Gets a vector of parents from an element
const vector <ParamRef>& Param::get_parent(unsigned int i) const
{
	auto ref = element_ref[i]; if(ref == UNSET) emsg("Cannot get element ref");
	return element[ref].parent;
}


/// Gets equation reference
unsigned int Param::get_eq_ref(unsigned int i) const
{
	auto ref = element_ref[i]; if(ref == UNSET) emsg("Cannot get element ref3");
	return element[ref].value.eq_ref;
}


/// Adds all elements
void Param::all_elements()
{
	for(auto i = 0u; i < N; i++) add_element(i);
}


/// Adds an element to the parameter
void Param::add_element(unsigned int i) 
{
	if(element_ref[i] == UNSET){
		element_ref[i] = element.size();
		ParamElement ele;
		//if(variety == REPARAM_PARAM) ele.value = add_equation_info(default_text,REPARAM);
		//else ele.value.te = default_text;
		ele.prior_ref = default_prior_ref;
		
		element.push_back(ele);
	}
}


/// Sets the prior
void Param::set_prior(unsigned int i, unsigned int prior_ref)
{
	add_element(i);
	element[element_ref[i]].prior_ref = prior_ref;
}


/// Adds a parent onto the list
void Param::add_parent(unsigned int i, const ParamRef &pr)
{
	auto ref = element_ref[i]; if(ref == UNSET) emsg("Cannot get element ref4");
	return element[ref].parent.push_back(pr);
}


/// Adds a child onto the list
void Param::add_child(unsigned int i, const ParamRef &pr)
{
	auto ref = element_ref[i]; if(ref == UNSET) emsg("Cannot get element ref");
	return element[ref].child.push_back(pr);
}


/// Gets the value of an element
double Param::get_value(unsigned int i) const 
{
	auto ref = element_ref[i]; if(ref == UNSET) emsg("Cannot get element ref5");
	if(variety == CONST_PARAM) return cons[ref];
	return element[ref].value.value;
}


/// Gets the value of an element
unsigned int Param::get_prior_ref(unsigned int i) const 
{
	auto ref = element_ref[i]; if(ref == UNSET) emsg("Cannot get element ref6");
	return element[ref].prior_ref;
}



/// Sets an element as used
void Param::set_used(unsigned int i) 
{
	auto ref = element_ref[i]; if(ref == UNSET) emsg("Cannot get element ref7");
	if(variety != CONST_PARAM) element[ref].used = true;
}


/// Sets the value of an element
void Param::set_value_te(unsigned int i, string te) 
{
	add_element(i);
	element[element_ref[i]].value.te = te;
}


/// Gets the value of an element
string Param::get_value_te(unsigned int i) const
{
	auto ref = element_ref[i]; if(ref == UNSET) emsg("Cannot get element ref8");
	return element[ref].value.te;
}


/// Sets the value of an element
void Param::set_value_eqn(unsigned int i, const EquationInfo &val) 
{
	add_element(i);
	element[element_ref[i]].value = val;
}

/*
/// Sets the value of a weight
void Param::set_weight(unsigned int i, double w) 
{
	auto ref = element_ref[i]; if(ref == UNSET) emsg("Cannot get element ref");
	element[element_ref[i]].weight = w;
}
*/

/*
/// Gets the value of a weight
double Param::get_weight(unsigned int i) 
{
	auto ref = element_ref[i]; if(ref == UNSET) emsg("Cannot get element ref");
	return element[element_ref[i]].weight;
}
*/


/// Gets the value of an element
bool Param::exist(unsigned int i) const
{
	if(element_ref[i] == UNSET) return false;
	return true;
}


/// Adds a constant value to the list
unsigned int Param::add_cons(double val)
{					
	auto num = cons.size();
	if(num < 10){
		auto i = 0u; while(i < num && cons[i] != val) i++;
		if(i < num) return i;
	}
	cons.push_back(val);
	return num;
}


/// Adds a constant value to the list
void Param::set_cons(unsigned int i, double val)
{					
	element_ref[i] = add_cons(val);
}

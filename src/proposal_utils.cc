/// Utility functions which are used in proposals

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <algorithm>

using namespace std;

#include "proposal.hh"
#include "utils.hh"
#include "matrix.hh"

/// Looks at population changes under and mbp and how they might affect other species
void Proposal::mbp_population_affect()
{
	for(auto po = 0u; po < model.pop.size(); po++){
		const auto &pop = model.pop[po];
		
		for(auto &pmr : pop.markov_eqn_ref){
			auto p = pmr.p, e = pmr.e;

			AffectLike al; al.map.resize(model.details.T,true);
			
			al.type = DIV_VALUE_AFFECT; al.num = p; al.num2 = e;
			param_vec_add_affect(affect_like,al);		
	
			al.type = MARKOV_LIKE_AFFECT;
			param_vec_add_affect(affect_like,al);
		}
	
		for(auto &tre : pop.trans_ref){
			auto p = tre.p, tr = tre.tr;
			if(p != p_prop){
				AffectLike al; al.map.resize(model.details.T,true);
				al.type = MARKOV_POP_AFFECT; al.num = p; al.num2 = tr;
				param_vec_add_affect(affect_like,al);
			}
		}
	}
	
	set_mbp_fast();
}


/// Gets how other parameter and priors are affect by change
void Proposal::get_dependency()
{
	vector <unsigned int> map(model.param_vec.size(),0);
	
	vector <unsigned int> list;
	for(auto i = 0u; i < N; i++){
		auto j = param_list[i];
		list.push_back(j);
		map[j] = 1;
	}
	
	for(auto i = 0u; i < list.size(); i++){
		auto j = list[i]; 
		
		const auto &pv = model.param_vec[j];
		const auto &par = model.param[pv.th];
	
		for(const auto &ch : par.get_child(pv.index)){
			const auto &par_child = model.param[ch.th];
			
			if(par_child.variety == REPARAM_PARAM){
				auto k = par_child.get_param_vec(ch.index);
				if(k != UNSET){
					if(map[k] == 0){ map[k] = 2; list.push_back(k);}
				}
			}
		}
	}
	
	// Adds all parameters which are dependent
	for(auto i = 0u; i < model.param_vec.size(); i++){
		if(map[i] == 2) dependent.push_back(i);
	}
	
	if(false){
		cout << name << ":";
		for(auto th : dependent) cout << model.param_vec[th].name <<","; 
		cout << endl;
		emsg("done");
	}
}


// Gets how the likelihood is altered under the change
void Proposal::get_affect_like()
{
	auto vec = param_list;
	for(auto th : dependent) vec.push_back(th);
	
	for(auto i = 0u; i < vec.size(); i++){
		auto th = vec[i];
		for(const auto &al : model.param_vec[th].affect_like){
			param_vec_add_affect(affect_like,al);
		}
	}
	
	switch(type){
	case IE_COVAR_PROP:
		{
			const auto &ieg = model.species[p_prop].ind_eff_group[ind_eff_group_ref.ieg];
			for(const auto &li : ieg.list){
				model.add_ie_affect(p_prop,li.index,affect_like);
			}
		}
		break;
		
	case PAR_EVENT_FORWARD_PROP: case PAR_EVENT_FORWARD_SQ_PROP: 
	case PAR_EVENT_BACKWARD_SQ_PROP: 
		model.joint_affect_like(type,tr_change,p_prop,affect_like);
		break;
	default: break;
	}

	model.add_iif_w_affect(affect_like);
		
	model.add_popnum_ind_w_affect(affect_like);

	if(nopop_speedup && type == PARAM_PROP){   
		model.affect_nopop_speedup(affect_like,dependent_update_precalc,update_precalc);
	}
	
	model.order_affect(affect_like);
	
	if(linearise_speedup){
		model.affect_linearise_speedup(affect_like);
	}
	
	if(linearise_factor_nopop_speedup && type == PARAM_PROP){ 
		model.set_factor_nopop_only(affect_like,dependent_update_precalc,update_precalc);
	}
	
	model.order_affect(affect_like);
}


/// Sets up the MVN sampling distribution from particle samples
void Proposal::set_mvn(double si_, const CorMatrix &cor_matrix)
{
	update_sampler(cor_matrix);	

	M_inv = invert_matrix(M);
	si = si_;
}


/// Initialise sampler
void Proposal::update_sampler(const CorMatrix &cor_matrix)
{
	if(!on || param_list.size() == 0) return;
	
	M = cor_matrix.find_covar(param_list);
	auto Mst = M;

	auto loop = 0u, loopmax = 1000u;
	do{
		auto illegal = false;
		Z = calculate_cholesky(M,illegal);
		if(!illegal) break;
		
		for(auto i = 0u; i < M.size(); i++){
			M[i][i] *= 1.1;
		}
		loop++;
	}while(loop < loopmax);
	if(loop == loopmax){
		/*
		//cout << name << "  proposal" << endl;
		{
			cout << "check" << endl;
			cor_matrix.check();
			
			ofstream fout("pout.txt");
			auto th = param_list[0];
			
			auto av = 0.0, av2 = 0.0, nav = 0.0;
			for(auto i = 0u; i < cor_matrix.n; i++){
				auto val = cor_matrix.samp[i][th];
				fout << i << " " << cor_matrix.n_start << " " << val << endl;
				if(i >= cor_matrix.n_start){
					av += val; av2 += val; nav++;
				}				
			}
			cout << av/nav << " " << (av2/nav) - (av/nav)*(av/nav) << " av" << endl;
		
		}
		
		
		*/
		//print_matrix("A",M);
		//print_matrix("B",Mst);
		
		//for(auto th : param_list){
			//auto pv = model.param_vec[th];
			//cout << pv.name << "pv" << endl;
		//}
		print_matrix("B",Mst);
		emsg("Cholesky problem");
	}
}

/// Samples from the covariance matrix
bool Proposal::param_resample(PV &param_val, const vector < vector <double> > &popnum_t)
{
	timer[PARAM_RESAMPLE_TIMER] -= clock();
	
	auto vec = sample_mvn(Z);

	auto &value = param_val.value;
	auto &precalc = param_val.precalc;
	
	for(auto i = 0u; i < N; i++){
		auto j = param_list[i];
		param_val.value_change(j);
		value[j] += si*vec[i];
		
		if(model.in_bounds(value[j],j,precalc) == false){
			param_val.restore(); 
			timer[PARAM_RESAMPLE_TIMER] += clock();
			return false;
		}
		
		model.param_update_precalc_after(j,param_val,true);
	}	

	for(auto k = 0u; k < dependent.size(); k++){
		auto j = dependent[k];

		const auto &pv = model.param_vec[j];
		const auto &par = model.param[pv.th]; 
		
		auto ref = par.get_eq_ref(pv.index);
		if(ref == UNSET) emsg("Reparam is not set");	
		
		const auto &dup = dependent_update_precalc[k];
		if(dup.list_precalc.size() > 0){
			model.precalc_eqn.calculate(dup.list_precalc,dup.list_precalc_time,param_val,true);
		}
		
		//model.param_update_precalc(j,param_val,true);
		
		
		param_val.value_change(j);
		if(pv.reparam_time_dep == false) value[j] = model.eqn[ref].calculate_param(precalc);
		else{
			auto ti = pv.reparam_spl_ti;
			value[j] = model.eqn[ref].calculate(ti,popnum_t[ti],precalc);
		}
		
		if(model.in_bounds(value[j],j,precalc) == false){
			param_val.restore(); 
			timer[PARAM_RESAMPLE_TIMER] += clock();
			return false;
		}
	
		model.param_update_precalc_after(j,param_val,true);
	}

	if(omega_check && model.ie_cholesky_error(param_val)){
		param_val.restore(); 
		timer[PARAM_RESAMPLE_TIMER] += clock();
		return false;
	}
	
	for(const auto &pr_up : update_precalc){
		model.precalc_eqn.calculate(pr_up.list_precalc,pr_up.list_precalc_time,param_val,true);
	}
	
	timer[PARAM_RESAMPLE_TIMER] += clock();
	
	return true;
}


/// Returns the MVN probability (or something proportional to it
double Proposal::mvn_probability(const vector <double> &param_prop1, const vector <double> &param_prop2) const
{
	vector <unsigned int> list;
	for(auto th : param_list) list.push_back(model.param_vec[th].ref);

	//for(auto i = 0u; i < param_list.size(); i++) cout << param_list[i] << " " << list[i] << " comp" << endl;
	//emsg("J");

	double sum = 0;
	for(auto v1 = 0u; v1 < N; v1++){
		for(auto v2 = 0u; v2 < N; v2++){
			double val1 = param_prop1[list[v1]] - param_prop2[list[v1]];
			double val2 = param_prop1[list[v2]] - param_prop2[list[v2]];
			sum += (val1/si)*M_inv[v1][v2]*(val2/si);
		}
	}	
	return -0.5*sum;
}


/// Proposes a new initial condition
ICResult Proposal::propose_init_cond(InitCondValue &icv, const State &state)
{
	vector < vector <double> > num;
	
	const	auto &sp = model.species[p_prop];

	const auto &ic = sp.init_cond;
	
	auto foc_cl = ic.focal_cl;
	
	if(foc_cl == UNSET){		
		switch(type){
		case MBP_IC_POPTOTAL_PROP:
			{	
				string warn;
				auto dN_total = normal_int_sample(si,warn);
				if(dN_total == UNSET) emsg("MCMC proposal problem: "+warn);
		
				if(-dN_total > int(icv.N_total)) return IC_FAIL;
				
				icv.N_total += dN_total;
		
				if(prior_probability(icv.N_total,ic.pop_prior,state.param_val.precalc,model.eqn) == -LARGE){
					return IC_FAIL;
				}
		
				if(dN_total > 0){
					auto dcnum = multinomial_sample(dN_total,icv.frac);
					for(auto c = 0u; c < sp.comp_gl.size(); c++) icv.cnum[c] += dcnum[c];
				}
				else{
					if(dN_total < 0){
						multinomial_reduce(-dN_total,icv.cnum,icv.frac);
					}
				}
			}
			break;
			
		case MBP_IC_RESAMP_PROP: // Resamples initial 
			icv.cnum = multinomial_resample(icv.cnum,icv.frac,si);
			break;
			
		default: emsg("prop prob"); break;
		}
	}
	else{
		switch(type){
		case MBP_IC_POP_PROP:
			{	
				auto cpr = c_prop;
				
				string warn;
				auto dN_focal = normal_int_sample(si,warn);
				if(dN_focal == UNSET) emsg("MCMC proposal problem: "+warn);
		
				if(-dN_focal > int(icv.N_focal[cpr])) return IC_FAIL;
				
				icv.N_focal[cpr] += dN_focal;
		
				if(prior_probability(icv.N_focal[cpr],ic.comp_prior[cpr],state.param_val.precalc,model.eqn) == -LARGE){
					return IC_FAIL;
				}
			
				if(dN_focal == 0) return IC_ZERO;
				
				if(dN_focal > 0){
					auto dcnum_reduce = multinomial_sample(dN_focal,icv.frac_comb);
					for(auto c = 0u; c < ic.N_reduce; c++) icv.cnum_reduce[cpr][c] += dcnum_reduce[c];
				}
				else{
					if(dN_focal < 0){
						multinomial_reduce(-dN_focal,icv.cnum_reduce[cpr],icv.frac_comb);
					}
				}
				
				model.combine_cnum_reduce(p_prop,icv);
			}
			break;
			
			
		case MBP_IC_RESAMP_PROP: // Resamples initial 
			{
				auto foc_cl = ic.focal_cl;
				const auto &claa = sp.cla[foc_cl]; 
				
				for(auto c = 0u; c < claa.ncomp; c++){
					icv.cnum_reduce[c] = multinomial_resample(icv.cnum_reduce[c],icv.frac_comb,si);
				}
				model.combine_cnum_reduce(p_prop,icv);
			}
			break;
			
		default: emsg("prop prob"); break;
		}
	}
	
	return IC_SUCCESS;
}

/// Resamples from a multinomial distribution
vector <unsigned int> Proposal::multinomial_resample(const vector <unsigned int> &x, const vector <double> &frac, double prob) const
{
	auto N = frac.size();
	auto x_new = x;
	
	if(false){ // Direct way
		vector <unsigned int> list;
		for(auto i = 0u; i < N; i++){
			for(auto j = 0u; j < x[i]; j++) list.push_back(i);
		}
		
		auto sum = 0.0;
		vector <double> frac_sum(N);
		for(auto i = 0u; i < N; i++){
			sum += frac[i];	
			frac_sum[i] = sum;
		}
		if(dif(sum,1,DIF_THRESH)) emsg("sum prob");
		
		for(auto k = 0u; k < list.size(); k++){
			if(ran() < prob){
				auto z = ran();
				auto i = 0u; while(i < N && z > frac_sum[i]) i++;
				if(i == N) emsg("Samp prob");
				
				x_new[list[k]]--;
				x_new[i]++;
			}
		}
	}
	else{  // Faster way
		string warn;
		
		auto total = 0u;
		for(auto i = 0u; i < N; i++){
			auto n = binomial_sample(prob,x[i],warn);
			if(n == UNSET) emsg("Problem sampling from bionomial distribution: "+warn);
			x_new[i] -= n;
			total += n;
		}
		
		auto x_resamp = multinomial_sample(total,frac);
		
		for(auto i = 0u; i < N; i++){
			x_new[i] += x_resamp[i];
		}	
		
		auto totaft = 0u;	for(auto i = 0u; i < N; i++) totaft += x_new[i];
	}
	
	return x_new;
}


/// Reduces the multinomial distribution by a certain number of individuals
void Proposal::multinomial_reduce(unsigned int dN, vector <unsigned int> &x, const vector <double> &frac) const
{
	auto N = frac.size();
	if(true){
		vector <unsigned int> list;
		for(auto i = 0u; i < N; i++){
			for(auto j = 0u; j < x[i]; j++) list.push_back(i);
		}
	
		auto nlist = list.size();
		for(auto j = 0u; j < dN; j++){
			if(nlist == 0) emsg("problem reduce");
			auto k = (unsigned int)(ran()*nlist);
			x[list[k]]--;
			if(k < nlist-1) list[k] = list[nlist-1];
			list.pop_back();
			nlist--;
		}
	}
	else{ // TO DO Need a faster way 	
	}
}


/// Adds possibilities for local changes from add to the total list
void Proposal::add_li_cha(unsigned int i, const vector <LocalIndChange> &add, vector <LocalIndChange> &licha, vector < vector <unsigned int> > &lc_ref) const
{
	for(auto ad : add){
		ad.ref = lc_ref[i].size();
		lc_ref[i].push_back(licha.size());
		licha.push_back(ad);
	}
}


/// Removes possibilities for local changes for individual i
void Proposal::remove_li_cha(unsigned int i, vector <LocalIndChange> &licha, vector < vector <unsigned int> > &lc_ref) const
{
	for(auto j : lc_ref[i]){
		if(j+1 < licha.size()){
			auto &lc = licha[j];
			lc = licha[licha.size()-1];
			lc_ref[lc.i][lc.ref] = j;
		}
		licha.pop_back();
	}
	lc_ref[i].clear();
}


/// Sets up a list and map
ListMap::ListMap(unsigned int min_, unsigned int max_)
{
	min = min_;
	map.resize(min,UNSET);
	for(auto i = min_; i < max_; i++){
		map.push_back(list.size());
		list.push_back(i);
	}
}


/// Adds an element
void ListMap::add(unsigned int i)
{
	map.push_back(list.size());
	list.push_back(i);
}


/// Removes an element (but leaves a hole)
void ListMap::remove(unsigned int i)
{
	auto m = map[i];
	map[i] = UNSET;
	if(m < list.size()-1){
		list[m] = list[list.size()-1];
		map[list[m]] = m;
	}
	list.pop_back();
}


/// Puts an element back in hole
void ListMap::add_hole(unsigned int i)
{
	map[i] = list.size();
	list.push_back(i);
}


/// Removes a hole
void ListMap::remove_hole(unsigned int i)
{
	if(map[i] != UNSET) emsg("not hole");
	if(i < map.size()-1){
		map[i] = map[map.size()-1];
		list[map[i]] = i;
	}
	map.pop_back();
}


/// Checks that every thing is correct
void ListMap::check(unsigned int num)
{
	if(list.size() != num) emsg("unobs prob");

	for(auto j = 0u; j < list.size(); j++){
		if(map[list[j]] != j){
			emsg("inconstent1");
		}
	}			
		
	for(auto i = 0u; i < list.size(); i++){
		if(i < min){
			if(map[i] != UNSET) emsg("should be unset");
		}
		else{
			if(list[map[i]] != i) emsg("inconstent2");
		}
	}
}


/// Initialises sampling variables for individual proposals
void Proposal::initialise_ind_variable()
{
	const auto &sp = model.species[p_prop];

	ind_sampler.resize(sp.tra_gl.size());
	for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
		auto &is = ind_sampler[tr];
		is.si = 1;
		is.ntr = 0;
		is.nac = 0;
		is.nfa = 0;
	}
}


/// Initialises sampling variables for individual proposals
void Proposal::initialise_ind_time_sampler()
{
	const auto &sp = model.species[p_prop];

	ind_sampler.resize(sp.nindividual_obs);
	for(auto i = 0u; i < sp.nindividual_obs; i++){
		auto &is = ind_sampler[i];
		is.si = 1;
		is.ntr = 0;
		is.nac = 0;
		is.nfa = 0;
	}
}


/// Initialises sampling variables for swap proposals
void Proposal::initialise_swap_variable()
{
	const auto &claa = model.species[p_prop].cla[cl_prop];

	swap_info.resize(claa.swap_rep.size());
	for(auto j = 0u; j < claa.swap_rep.size(); j++){
		auto &si = swap_info[j];
		si.ntr = 0; si.ntr2 = 0;
		si.nac = 0;
		si.nfa = 0;
		si.nzero = 0;
		si.nfilt = 0;
	}
}



/// Updates the size of the proposal (if in burnin phase)
void Proposal::update_si(double fac)
{
	if(burn_info.on){
		auto f = 1+fac*burn_info.fac;
		if(f < 0.5) f = 0.5; 
		if(f > 10) f = 10;
		si *= f;
		if(si < 0.001) si = 0.001;
		if(si > 1000) si = 1000;
	}
}


/// Updates the size of individual sampler (if in burnin phase)
void Proposal::update_ind_samp_si(unsigned int tr_gl, double fac)
{
	if(burn_info.on){
		auto &si = ind_sampler[tr_gl].si;
		si *= 1+fac*burn_info.fac;
		if(si < SD_MIN) si = SD_MIN;
	}
}


/*
/// Deterines probability proposal won't be performed (assuming no prob adaptation)
bool Proposal::skip_proposal(double val) const 
{
	if(!adapt_prop_prob && ran() < val) return true;
	return false;
}
*/


/// Determines if omega needs to be checked over proposal
void Proposal::set_omega_check()
{
	omega_check = false;
	
	for(auto th : param_list){
		if(model.param_vec[th].omega_fl) omega_check = true;
	}
}


/// Determines if two sampled sequences are different
bool Proposal::event_dif(const vector <Event> &ev1, const vector <Event> &ev2) const
{
	if(ev1.size() != ev2.size()) return true;
	
	for(auto i = 0u; i < ev1.size(); i++){
		if(ev1[i].c_after != ev2[i].c_after) return true;
		if(ev1[i].tdiv != ev2[i].tdiv) return true;
	}
	
	return false;
}


/// Update probability of doing proposal
void Proposal::ind_obs_prob_update(IndSimProb &isp) const
{
	isp.ntr *= PROP_SIM_PROB_FADE;
	isp.nac *= PROP_SIM_PROB_FADE;
	auto pr = 2*isp.nac/isp.ntr;	
	if(pr > 1) pr = 1; 
	if(pr < 0.05) pr = 0.05;
	isp.prob = pr;
	isp.done = true;
}


/// Sets quantities used to speed up MBPs
void Proposal::set_mbp_fast()
{
	const auto &sp = model.species[p_prop];
	
	vector <LinearFormInit> lfinit;
	
	for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
		const auto &tra = sp.tra_gl[tr];
		auto e = tra.dist_param[0].eq_ref;
		const auto &eq = model.eqn[e];
		
		if(eq.linearise.on && eq.linearise.pop_grad_time_dep == false){
			LinearFormInit lfi; lfi.m = tr; lfi.e = e;
			lfinit.push_back(lfi);
		}
		else{
			mbp_fast.calc_tr.push_back(tr);
		}
	}
	
	model.set_linear_form(p_prop,mbp_fast.lin_form,lfinit);
}


/// Sets the proposal probability 
double Proposal::set_prop_prob()
{
	switch(type){
	case IND_OBS_SAMP_PROP:
		{
			// This works out a metric for the time taken to do generate_ind_obs_timeline()
			const auto &sp = model.species[p_prop];
			
			auto effort = 0.0;
			for(auto cl = 0u; cl < sp.ncla; cl++){
				const auto &claa = sp.cla[cl];
		
				auto av = 0.0, nav = 0.0;
				for(auto i = 0u; i < claa.island.size(); i++){
					auto sum = 0u;
					const auto &comp = claa.island[i].comp;
					if(comp.size() > 1){
						for(const auto &co : comp) sum += (1+co.leave.size());
					}
					av += sum; nav++;
				}				
				av /= nav;
				effort += av; 
			}
			
			auto pr = 0.5/effort;
			if(pr > 0.4) pr = 0.4;
			if(pr < 0.01) pr = 0.01;
			return pr;
		}
		
	case IE_PROP: case IE_VAR_PROP: case IE_COVAR_PROP: case IE_VAR_CV_PROP:
		{
			const auto &sp = model.species[p_prop];
			
			auto g = UNSET;
			switch(type){
			case IE_PROP: case IE_VAR_PROP: case IE_VAR_CV_PROP: g = sp.ind_effect[ie_prop].index; break;
			case IE_COVAR_PROP: g = ind_eff_group_ref.ieg; break;
			default: break;
			}			
			
			const auto &ieg = sp.ind_eff_group[g];
	
			auto N = ieg.list.size();
			auto effort = N;
			const auto &A = ieg.A_matrix;
			if(A.set == true) effort *= A.sparsity*A.value.size();
			
			auto pr = 4000.0/(2000+effort);
			
			switch(type){
			case IE_PROP: 
				pr *= 0.2;
				if(pr > 0.2) pr = 0.2;
				if(pr < 0.05) pr = 0.05;
				break;
				
			default:
				if(pr > 1) pr = 1;
				if(pr < 0.1) pr = 0.1;
				break;
			}
			
			return pr;
		}
		
	case PARAM_PROP: return 0.5;

	case IND_OBS_RESIM_SINGLE_PROP: return 0.5;
	case IND_OBS_RESIM_PROP: return 0.5;
	case IND_MULTI_EVENT_PROP: return 0.5;
	case IND_EVENT_ALL_PROP: return 0.5;
	case TRANS_TREE_PROP: return 0.1;
	case TRANS_TREE_SWAP_INF_PROP: return 0.1;
	case TRANS_TREE_MUT_PROP: return 0.1;
	case TRANS_TREE_MUT_LOCAL_PROP: return 0.1;
	case POP_SINGLE_LOCAL_PROP: return 0.5;
	case POP_ADD_REM_LOCAL_PROP: return 0.5;
	default: break;
	}
	return 1;
}

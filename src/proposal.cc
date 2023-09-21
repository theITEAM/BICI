/// Implements MCMC proposals

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include "proposal.hh"
#include "ind_ev_sampler.hh"
#include "utils.hh"
#include "matrix.hh"

Proposal::Proposal(PropType type_, vector <unsigned int> vec, const Model &model, const Output &output, double prob, const double &burnfac) : model(model), output(output), burnfac(burnfac)
{	
	type = type_;
	prop_prob = prob;
	
	timer = 0;
	name = "PROPOSAL ";
	
	switch(type){
	case PARAM_PROP: case MBP_PROP: case IE_VAR_PROP:   // Parameter proposals
		if(type == PARAM_PROP) name += "PARAM_PROP ";
		if(type == MBP_PROP) name += "MBP_PROP ";
		if(type == IE_VAR_PROP){
			name += "IE_VAR_PROP ";
			if(vec.size() != 3) emsg("ert");
			p_prop = vec[1]; ie_prop = vec[2];
			vec.resize(1);
		}

		param_list = vec;
	
		initialise_variables();                 // Initialises the basic sampling variables
	
		get_dependency();                       // Gets how other parameter and priors are affect by change
	
		get_affect_like();                      // Gets how the likelihood is altered under the change
		break;
		
	case IND_PROP:                            // Individual proposals
		name += "IND_PROP";
		if(vec.size() != 1) emsg("error");
		p_prop = vec[0];
		
		initialise_ind_variable();              // Initialises sampling variables for individual proposals
		break;
		
	case INDSAMP_PROP:                            // Individual proposals
		name += "INDSAMP_PROP";
		if(vec.size() != 1) emsg("error");
		p_prop = vec[0];
		ntr = 0; nac = 0; si = 1; N = 0;
		break;
		
	case MBPII_PROP:
		name += "MBPII_PROP";
		if(vec.size() != 1) emsg("error");
		p_prop = vec[0];
		ntr = 0; nac = 0; si = 1; N = 0;
		mbp_population_affect();
		break;
		
	case IE_PROP:                            // Individual proposals
		name += "IE_PROP";
		if(vec.size() != 2) emsg("error");
		p_prop = vec[0];
		ie_prop = vec[1];
		ntr = 0; nac = 0; si = 1; N = 0;
		break;
	}
}


/// Looks at population changes under and mbp and how they might affect other species
void Proposal::mbp_population_affect()
{
	for(auto po = 0u; po < model.pop.size(); po++){
		const auto &pop = model.pop[po];
		
		for(auto &pmr : pop.markov_eqn_ref){
			auto p = pmr.p, e = pmr.e;

			AffectLike al; al.map.resize(model.ntimepoint-1,true);
			
			al.type = DIV_VALUE_AFFECT; al.num = p; al.num2 = e;
			param_vec_add_affect(affect_like,al);		
	
			al.type = MARKOV_LIKE_AFFECT;
			param_vec_add_affect(affect_like,al);
		}
	
		for(auto &tre : pop.trans_ref){
			auto p = tre.p, tr = tre.tr;
			if(p != p_prop){
				AffectLike al; al.map.resize(model.ntimepoint-1,true);
				al.type = MARKOV_POP_AFFECT; al.num = p; al.num2 = tr;
				param_vec_add_affect(affect_like,al);
			}
		}
	}
}


/// Initialises the basic sampling variables
void Proposal::initialise_variables()
{
	N = param_list.size();
	
	name += " for: ";
	for(auto j = 0u; j < N; j++){
		if(j != 0) name += ",";
		name += model.param_vec[param_list[j]].name;
	}
	
	ntr = 0; nac = 0; si = 1;
	
	setup_sum();

	M.resize(N);
	for(auto j = 0u; j < N; j++) M[j].resize(N,UNSET);
}

	
// Gets how other parameter and priors are affect by change
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
		
		for(const auto &ch : par.child[pv.index]){
			const auto &par_child = model.param[ch.th];
			
			if(par_child.variety == REPARAM_PARAM){
				auto k = par_child.param_vec_ref[ch.index];
				if(k == UNSET) emsg("Should be set2");
			
				if(map[k] == 0){ map[k] = 2; list.push_back(k);}
			}
		}
	}
	
	// Adds all parameters which are dependent
	for(auto i = 0u; i < model.param_vec.size(); i++){
		if(map[i] == 2) dependent.push_back(i);
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
	
	model.order_affect(affect_like);
}


/// Sets up the MVN sampling distribution from particle samples
void Proposal::set_mvn(vector <Particle> &particle, double si_)
{
	setup_sum();
	for(const auto &part : particle){
		add_sample(part.param_val);
	}
	update_sampler();	
	M_inv = invert_matrix(M);
	
	si = si_;
}


/// Zeros the sums used to calculate the covariance matrix
void Proposal::setup_sum()
{
	M_sum.clear();
	M_sum.resize(N);
	for(auto j = 0u; j < N; j++) M_sum[j].resize(N,0);
	
	mean_sum.clear();
	mean_sum.resize(N,0);
	nsum = 0;
}
	
	
/// Adds sample to proposal distribution (to allow for adaptive proposals
void Proposal::add_sample(const vector <double> &param_val)
{
	for(auto j = 0u; j < N; j++){
		mean_sum[j] += param_val[param_list[j]];
		for(auto i = 0u; i < N; i++){
			M_sum[j][i] += param_val[param_list[j]]*param_val[param_list[i]];
		}
	}
	nsum++;
}


/// Initialise sampler
void Proposal::update_sampler()
{
	for(auto j = 0u; j < N; j++){
		for(auto i = 0u; i < N; i++){
			auto val = M_sum[j][i]/nsum - (mean_sum[j]/nsum)*(mean_sum[i]/nsum);
			if(val < 0) val = 0;
			M[j][i] = val;
		}
	}
	
	Z = calculate_cholesky(M);
}

/// Samples from the covariance matrix
vector <double> Proposal::sample(vector <double> param_val) const
{
	auto vec = sample_mvn(Z);
	for(auto i = 0u; i < N; i++){
		param_val[param_list[i]] += si*vec[i];
	}	
	
	for(auto j : dependent){
		const auto &pv = model.param_vec[j];
		auto ref = model.param[pv.th].value[pv.index].eq_ref;
		if(ref == UNSET) emsg("Reparam is not set");	
		param_val[j] = model.eqn[ref].calculate_param_only(param_val);
	}
	
	return param_val;
}


/// Returns the MVN probability (or something proportional to it
double Proposal::mvn_probability(const vector <double> &param1, const vector <double> &param2) const
{
	double sum = 0;
	for(auto v1 = 0u; v1 < N; v1++){
		for(auto v2 = 0u; v2 < N; v2++){
			double val1 = param1[param_list[v1]] - param2[param_list[v1]];
			double val2 = param1[param_list[v2]] - param2[param_list[v2]];
			sum += (val1/si)*M_inv[v1][v2]*(val2/si);
		}
	}	
	return -0.5*sum;
}


/// Prints information about the sampler
void Proposal::print_info() const
{
	if(com_op == true) return;
	
	switch(type){
	case PARAM_PROP: case MBP_PROP: case IE_VAR_PROP:
		{
			if(type == PARAM_PROP) cout << "PARAMETER SAMPLER for ";
			if(type == MBP_PROP) cout << "MBP SAMPLER for ";
			if(type == IE_VAR_PROP) cout << "IE VAR SAMPLER for ";
			
			for(auto th : param_list) cout << model.param_vec[th].name << ", ";
			cout << endl;
			
			for(auto th : dependent) cout << model.param_vec[th].name << ", ";
			cout << "dependent" << endl;
		}
		break;
		
	case IND_PROP:
		{
			cout << "INDIVIDUAL event time sampler for species " << model.species[p_prop].name << endl;
		}
		break;
		
	case INDSAMP_PROP:
		{
			cout << "INDIVIDUAL sampler for species " << model.species[p_prop].name << endl;
		}
		break;
		
	case MBPII_PROP:
		{
			cout << "MBP TYPE II for species "  << model.species[p_prop].name << endl;
		}
		break;
		
	case IE_PROP:
		{
			cout << "IE sampler for " << model.species[p_prop].ind_effect[ie_prop].name << endl;
		}
		break;
	}
			
	for(const auto &al : affect_like){
		cout << output.print_affect_like(al) << endl;
	}
	
	cout << endl;
}


/// Performs a Metropolis-Hastings proposal update
void Proposal::MH(State &state, bool burnin)
{
	auto param_store = state.param_val;
	
	state.param_val = sample(state.param_val);
	
	ntr++; //cout << "try\n";
	
	if(model.inbounds(state.param_val) == false){ state.param_val = param_store; return;}
	
	auto like_ch = state.update_param(affect_like);

	auto al = calculate_al(like_ch,0);
	
	if(ran() < al){ //cout << "ac\n";
		nac++;
		state.add_like(like_ch);
		if(burnin) si *= 1+0.01*burnfac;
	}
	else{
		if(burnin) si *= 1-0.005*burnfac;
		
		state.param_val = param_store;
		state.restore(affect_like);
	}
}


/// Performs a MBP proposal update
void Proposal::mbp(MBPType type, State &state, bool burnin)
{
	auto param_store = state.param_val;
	
	auto sim_prob = 0.0;
	switch(type){
	case TYPE_I: state.param_val = sample(state.param_val); break;
	case TYPE_II: sim_prob = si; break;
	}
	
	ntr++;//cout << "try\n";
	
	if(model.inbounds(state.param_val) == false){ state.param_val = param_store; return;}
	
	auto &ssp = state.species[p_prop];
	if(ssp.type != POPULATION) emsg("must be population");
	
	auto popnum_t_st = state.popnum_t;
	auto trans_num_st = ssp.trans_num;
	auto tnum_mean_st_st = ssp.tnum_mean_st;
	auto cpop_st_st = ssp.cpop_st;
		
	ssp.mbp(sim_prob,state.popnum_t);
	
	auto like_ch = state.update_param(affect_like);
	
	auto al = calculate_al(like_ch,0);
	
	if(ran() < al){ //cout << "ac\n";
		nac++;
		state.add_like(like_ch);
		ssp.mbp_accept(state.like.markov,state.popnum_t);
		
		if(burnin){ si *= 1+0.01*burnfac; if(si > 1 && type == TYPE_II) si = 1;}
	}
	else{
		if(burnin) si *= 1-0.005*burnfac;
		state.param_val = param_store;
		state.popnum_t = popnum_t_st;
		ssp.trans_num = trans_num_st;
		ssp.tnum_mean_st = tnum_mean_st_st;
		ssp.cpop_st = cpop_st_st;

		state.restore(affect_like);		
	}
}


/// Performs a Metropolis-Hastings proposal update for an individual event time
void Proposal::MH_ind(State &state, bool burnin)
{
	auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	
	auto t_start = model.details.t_start;
	auto t_end = model.details.t_end;
	
	for(auto i = 0u; i < ssp.individual.size(); i++){
		auto &ind = ssp.individual[i];
	
		auto cinit = ind.cinit;
	
		for(auto cl = 0u; cl < sp.ncla; cl++){
			auto jmax = ind.cl_ev[cl].size();
			if(jmax > 0){
				auto j = (unsigned int) (ran()*jmax);
				
				const auto &cl_ev = ind.cl_ev[cl][j];
				auto index = cl_ev.index_ev;
			
				auto &ev = ind.ev[index];
			
				if(ev.type == M_TRANS_EV && ev.fixed == false){
					auto t = ev.t;
			
					auto tmin = t_start; 
					if(j > 0){
						const auto &cl_ev_min = ind.cl_ev[cl][j-1];
						tmin = ind.ev[cl_ev_min.index_ev].t; 
					}
					
					auto tmax = t_end; 
					if(j < jmax-1){
						const auto &cl_ev_max = ind.cl_ev[cl][j+1];
						tmax = ind.ev[cl_ev_max.index_ev].t; 
					}
				
					auto tr_gl = ev.tr_gl;
					auto &samp = ind_sampler[tr_gl];
					
					auto t_new = t + normal_sample(0,samp.si);
		
					// Stops events bassing through birth 
					if(index > 0 && t_new < ind.ev[0].t && cinit == SOURCE) tmin = LARGE;
					
					auto emax = ind.ev.size()-1;
					const auto &evmax = ind.ev[emax];
					if(index < emax && t_new > evmax.t && evmax.type == M_TRANS_EV){
						if(sp.tra_gl[evmax.tr_gl].variety == SINK_TRANS) tmax = -LARGE;
					}

					samp.ntr++;
					if(t_new <= tmin || t_new >= tmax){
						if(burnin) samp.si *= 1-0.01*burnfac;
					}
					else{
						auto ev_store = ind.ev;
						auto ev_new = ev_store;
						
						auto tr_gl_new = ssp.move_event(cinit,ev_new,index,t_new);
						
						auto like_ch = state.update_ind(p_prop,i,ind.cinit,ev_new);
					
						auto dprob = 0.0;
						if(tr_gl_new != UNSET){  // Accounts for differeces in time sampler
							const auto &samp_new = ind_sampler[tr_gl_new];
							auto dt = t_new-t;
							dprob = normal_probability(dt,0,samp_new.si) - normal_probability(dt,0,samp.si); 
						}
						
						auto al = calculate_al(like_ch,dprob);
				
						//cout << "try ind\n";
						if(ran() < al){ //cout << "ac ind\n";
							samp.nac++;
							state.add_like(like_ch);
							
							if(tr_gl_new != UNSET) ssp.recalculate_cl_ev(ind); 
							
							if(burnin) samp.si *= 1+0.005*burnfac;
						}
						else{
							if(burnin) samp.si *= 1-0.01*burnfac;
							
							state.update_ind(p_prop,i,ind.cinit,ev_store);
						}
						
						//state.check("  ind prop");
					}
				}
			}
		}
	}
}


/// Performs a sampling proposal update for an individual event time
void Proposal::sample_ind(State &state)
{
	auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	
	IndEvSampler ind_ev_samp(ssp.markov_eqn_vari,ssp.individual,model.details,sp,ssp.obs_eqn_value,model.eqn);
		
	for(auto i = 0u; i < ssp.individual.size(); i++){
		auto &ind = ssp.individual[i];
	
		for(auto cl = 0u; cl < sp.ncla; cl++){
			//ssp.print_event("start",ind);
			//sp.print_obs_data(ind.name,sp.individual[i].obs);
						
			if(ind_ev_samp.needed(i,cl) == true){
				ntr++;
				
				ind_ev_samp.generate_ind_obs_timeline();
				auto probif = 0.0;
				auto ev_new = ind_ev_samp.sample_events(probif);
			
				if(ind_ev_samp.illegal == false){
					if(testing == true){ // Diagnostic
						auto prob = ind_ev_samp.sample_events_prob(ev_new);
						if(dif(probif,prob)) emsg("sampler different");
					}
					
					auto ev_store = ind.ev;
					
					auto like_ch = state.update_ind(p_prop,i,ind.cinit,ev_new);
					
					//ssp.print_event("mid",ind);
					
					ind_ev_samp.generate_ind_obs_timeline();
					auto probfi = ind_ev_samp.sample_events_prob(ev_new);
					
					auto al = calculate_al(like_ch,probfi-probif);
				
					//cout << "try ind\n";
					if(ran() < al){ //cout << "ac ind\n";
						nac++;
						state.add_like(like_ch);
						
						ssp.recalculate_cl_ev(ind); 
					}
					else{
						state.update_ind(p_prop,i,ind.cinit,ev_store);
					}
					
					state.check("  ind prop");
				}
			}
		}
	}
}


/// Performs a Metropolis-Hastings proposal update for individual effects
void Proposal::MH_ie(State &state, bool burnin)
{
	auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	
	const auto &inde = sp.ind_effect[ie_prop];
	
	auto in_me = false; 
	if(inde.markov_eqn_ref.size() > 0) in_me = true;
	
	auto in_pop = false;
	for(const auto &po : model.pop){
		for(auto ie_pop : po.ind_eff_mult){
			if(ie_pop == ie_prop) in_pop = true;
		}
	}
	
	auto g = inde.index;
	auto num = inde.num;
	
	auto &iegs = ssp.ind_eff_group_sampler[g];

	auto var = iegs.omega[num][num];
	auto sd = sqrt(var);
	
	vector <double> store;
	
	for(auto i = 0u; i < ssp.individual.size(); i++){
		auto &ind = ssp.individual[i];
	
		auto value_st = ind.ie[ie_prop];
		auto value_prop = value_st + normal_sample(0,sd*si);
		
		auto exp_old = ind.exp_ie[ie_prop];
		auto exp_prop = exp(value_prop-0.5*var);
		
		auto factor = exp_prop/exp_old;
	
		Like like_ch;
	
		if(in_me == true){                       // If in a markov_eqn then updates likelihood
			store = ssp.likelihood_ie_change(i,ie_prop,factor,like_ch);
		}
		
		if(in_pop == true){                      // If in a population then updates likelihoods from that
 			state.update_ie_population(p_prop,i,ie_prop,exp_prop-exp_old,like_ch.markov);
		}
		
		auto dLi_ie = ssp.likelihood_indeff_change(i,ie_prop,value_prop,value_st);

		like_ch.ie += dLi_ie;
		
		auto al = calculate_al(like_ch,0);
	
		//cout << "try ind\n";
		ntr++;
		if(ran() < al){ //cout << "ac ind\n";
			nac++;
			state.add_like(like_ch);
			
			ind.ie[ie_prop] = value_prop;
			ind.exp_ie[ie_prop] = exp_prop;
			ssp.Li_ie[g] += dLi_ie;
			
			if(burnin) si *= 1+0.005*burnfac;
		}
		else{ 
			if(burnin) si *= 1-0.01*burnfac;
			
			if(in_me == true){
				ssp.likelihood_ie_change_restore(i,ie_prop,factor,store);
			}
			
			if(in_pop == true){                      // If in a population then updates likelihoods from that
				state.update_ie_population(p_prop,i,ie_prop,exp_old-exp_prop,like_ch.markov);
			}
		}
	
		//state.check("  ind prop");
	}
}


/// Performs a Metropolis-Hastings joint proposal update for individual effects and a variance
void Proposal::MH_ie_var(State &state, bool burnin)
{
	auto param_store = state.param_val;
	
	state.param_val = sample(state.param_val);
	
	auto th = param_list[0];
	
	if(state.param_val[th] <= 0){ state.param_val = param_store; return;}
	
	auto factor = sqrt(state.param_val[th]/param_store[th]);

	auto &ssp = state.species[p_prop];
	for(auto &ind : ssp.individual) ind.ie[ie_prop] *= factor;
	
	auto like_ch = state.update_param(affect_like);

	auto al = calculate_al(like_ch,-like_ch.ie);
	
	ntr++; //cout << "try\n";
	if(ran() < al){ //cout << "ac\n";
		nac++;
		state.add_like(like_ch);
		if(burnin) si *= 1+0.01*burnfac;
	}
	else{ 
		if(burnin) si *= 1-0.005*burnfac;
		
		state.param_val = param_store;
		state.restore(affect_like);
		for(auto &ind : ssp.individual) ind.ie[ie_prop] /= factor;
	}
}

/// Calculates acceptance probability base don changes in likelihoods
double Proposal::calculate_al(const Like &like_ch, double dprob) const
{					
	return exp(like_ch.obs + like_ch.prior + like_ch.spline_prior +
						like_ch.dist + like_ch.markov + like_ch.nm_trans + like_ch.ie + dprob);
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
	}
}


/// Provides diagnostic information about the proposals
string Proposal::diagnostics(long total_time) const
{
	stringstream ss;
	
	switch(type){
	case PARAM_PROP:
		{
			ss << "Parameter proposal: ";
			for(auto i = 0u; i < param_list.size(); i++){
				auto th =  param_list[i];
				if(i != 0) ss << ", ";
				ss << model.param_vec[th].name; 
			}
			ss << endl;
			ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "%   Size:" << si;
			ss << endl;
		}
		break;
		
	case MBP_PROP:
		{
			ss << "Type I MBP parameter proposal: ";
			for(auto i = 0u; i < param_list.size(); i++){
				auto th =  param_list[i];
				if(i != 0) ss << ", ";
				ss << model.param_vec[th].name; 
			}
			ss << endl;
			ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "%   Size:" << si;
			ss << endl;
		}
		break;
		
	case MBPII_PROP:
		{
			ss << "Type II MBP proposal for species " << model.species[p_prop].name << endl;
			ss << endl;
			ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "%   Probability:" << si;
			ss << endl;
		}
		break;
	
	case IND_PROP:
		{
			ss << "Individual event proposals: ";
			ss << endl;
			const auto &sp = model.species[p_prop];
			for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
				const auto &tra = sp.tra_gl[tr];
				const auto &samp = ind_sampler[tr];
				ss << tra.name << "   ";
				ss << " Acceptance: " << (unsigned int) (100.0*samp.nac/samp.ntr) << "%   Size:" << samp.si;
				ss << endl;
			}
		}
		break;
	
	case INDSAMP_PROP:
		{
			ss << "Individual sampler proposals: ";
			ss << endl;
			ss << endl;
			ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "% ";
			ss << endl;
		}
		break;
		
	case IE_PROP:
		{
			auto ie_name =  model.species[p_prop].ind_effect[ie_prop].name;
			ss << "Individual effect proposal for " << ie_name;
			ss << endl;
			ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "%   Size:" << si;
			ss << endl;
		}
		break;
		
	case IE_VAR_PROP:
		{
			auto j = param_list[0];
			auto th = model.param_vec[j].th;
			auto name = model.param[th].name;
			auto ie_name = model.species[p_prop].ind_effect[ie_prop].name;

			ss << "Joint proposal on " << name << " and " << ie_name;
			ss << endl;
			ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "%   Size:" << si;
			ss << endl;
		}
		break;
	}
	
	ss << "CPU time: " << cpu_percent(timer,total_time) << endl;
	
	return ss.str();
}


/// Updates the system based on the proposal type
void Proposal::update(State &state, bool burnin)
{
	if(ran() > prop_prob) return;
	
	//if(testing == true) cout << "proposal " << name << endl;
	
	timer -= clock();
	switch(type){
	case PARAM_PROP: MH(state,burnin); break;
	case IND_PROP: MH_ind(state ,burnin); break;
	case INDSAMP_PROP: sample_ind(state); break;
	case MBP_PROP: mbp(TYPE_I,state,burnin); break;
	case MBPII_PROP: mbp(TYPE_II,state,burnin); break;
	case IE_PROP: MH_ie(state,burnin); break;
	case IE_VAR_PROP: MH_ie_var(state,burnin); break;
	}	
		
	if(testing == true) state.check("end");
	
	timer += clock();
}

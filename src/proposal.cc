/// Implements MCMC proposals

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <algorithm>

using namespace std;

#include "proposal.hh"
#include "ind_ev_sampler.hh"
#include "utils.hh"
#include "matrix.hh"

Proposal::Proposal(PropType type_, vector <unsigned int> vec, const Model &model, const Output &output, double w,  const BurnInfo &burn_info) : model(model), output(output), burn_info(burn_info)
{	
	timer.resize(PROPTIMER_MAX,0);

	type = type_;
	prop_prob = 1;
	prop_weight = w;
	number = 0;
	
	name = "PROPOSAL ";
	
	auto win = (unsigned int)(model.ntimepoint/10); if(win == 0) win = 1;
	loc_samp.win = win;
	
	p_prop = UNSET;
	cl_prop = UNSET; 
	c_prop = UNSET;     
	
	on = true;
	
	N = 0;
	
	switch(type){
	case IND_OBS_SWITCH_ENTER_SOURCE_PROP:
		p_prop = vec[0];
		name += "IND_OBS_SWITCH_ENTER_SOURCE_PROP";	
		break;
		
	case IND_OBS_SWITCH_LEAVE_SINK_PROP:
		p_prop = vec[0];
		name += "IND_OBS_SWITCH_LEAVE_SINK_PROP";	
		break;
		
	case CORRECT_OBS_TRANS_PROP: 
		p_prop = vec[0];
		all_events_correct = false;
		name += "CORRECT_OBS_TRANS_PROP";	
		break;
		
	case PAR_EVENT_FORWARD_PROP: case PAR_EVENT_FORWARD_SQ_PROP:
	case PAR_EVENT_BACKWARD_SQ_PROP:
		switch(type){
		case PAR_EVENT_FORWARD_PROP: name += "PAR_EVENT_FORWARD_PROP "; break;
		case PAR_EVENT_FORWARD_SQ_PROP: name += "PAR_EVENT_FORWARD_SQ_PROP "; break;
		case PAR_EVENT_BACKWARD_SQ_PROP: name += "PAR_EVENT_BACKWARD_SQ_PROP "; break;
		default: emsg("not here"); break;
		}
	
		p_prop = vec[0];
		param_list.push_back(vec[1]);
		tr_change.resize(model.species[p_prop].tra_gl.size(),false);
		cl_prop = UNSET;
		for(auto k = 2u; k < vec.size(); k++){
			auto tr = vec[k];
			const auto &tra = model.species[p_prop].tra_gl[tr];
			if(cl_prop == UNSET) cl_prop = tra.cl;
			else{
				if(cl_prop != tra.cl) emsg("Different classifications");
			}
			tr_change[tr] = true;
		}		
			
		initialise_variables();                 // Initialises the basic sampling variables
	
		get_dependency();                       // Gets how other parameter and priors are affect by change
	
		get_affect_like();                      // Gets how the likelihood is altered under the change
	
		si = 0.1;
		break;
		
	case TRANS_TREE_PROP:
		name += "TRANS_TREE_PROP ";
		si = 0.1;
		break;
		
	case TRANS_TREE_SWAP_INF_PROP:
		name += "TRANS_TREE_SWAP_INF_PROP ";
		si = 0.1;
		break;
		
	case TRANS_TREE_MUT_PROP:
		name += "TRANS_TREE_MUT_PROP ";
		si = 1;
		break;
		
	case TRANS_TREE_MUT_LOCAL_PROP:
		name += "TRANS_TREE_MUT_LOCAL_PROP ";
		gen_mut_info.resize(GEN_UPDATE_MAX);
		break;
		
	case PARAM_PROP: case MBP_PROP:    // Parameter proposals
	case IE_VAR_PROP: case IE_COVAR_PROP: case IE_VAR_CV_PROP: 
		if(type == PARAM_PROP) name += "PARAM_PROP ";
		if(type == MBP_PROP) name += "MBP_PROP ";
		if(type == IE_VAR_PROP){
			name += "IE_VAR_PROP ";
			if(vec.size() != 3) emsg("ert");
			p_prop = vec[1]; ie_prop = vec[2];
			vec.resize(1);
		}
		
		if(type == IE_COVAR_PROP){
			name += "IE_COVAR_PROP ";
			if(vec.size() != 5) emsg("ert");
			p_prop = vec[1]; 
			ind_eff_group_ref.ieg = vec[2];
			ind_eff_group_ref.i = vec[3];
			ind_eff_group_ref.j = vec[4];
			vec.resize(1);
		}
		
		if(type == IE_VAR_CV_PROP){
			name += "IE_VAR_CV_PROP ";
			p_prop = vec[1]; ie_prop = vec[2];
			vec.erase(vec.begin()+1,vec.begin()+3);
		}

		param_list = vec;
	
		initialise_variables();                 // Initialises the basic sampling variables
	
		get_dependency();                       // Gets how other parameter and priors are affect by change
	
		get_affect_like();                      // Gets how the likelihood is altered under the change
		break;
		
	case IND_EVENT_TIME_PROP:                 // Individual event time proposals
		name += "IND_EVENT_TIME_PROP";	
		if(vec.size() != 1) emsg("error vec num");
		p_prop = vec[0];
		
		initialise_ind_variable();              // Initialises sampling variables for individual proposals
		break;
		
	case IND_MULTI_EVENT_PROP:                 // Individual multi event time proposals
		name += "IND_MULTI_EVENT_PROP";	
		if(vec.size() != 1) emsg("error vec num");
		p_prop = vec[0];
		
		initialise_ind_variable();              // Initialises sampling variables for individual proposals
		break;
	
	case IND_EVENT_ALL_PROP:                 // Individual event time proposals for all events
		name += "IND_EVENT_ALL_PROP";	
		if(vec.size() != 1) emsg("error vec num");
		p_prop = vec[0];
		
		initialise_ind_time_sampler();         // Initialises sampling variables for individuals
		break;
		
	case IND_LOCAL_PROP:                      // Individual event time proposals
		name += "IND_LOCAL_PROP";	
		if(vec.size() != 2) emsg("error vec num");
		p_prop = vec[0]; cl_prop = vec[1];
		
		initialise_swap_variable();             // Initialises sampling variables for swap proposals
		break;
	
	case IND_OBS_RESIM_PROP:  
		{
			name += "IND_OBS_RESIM_PROP"; 
			if(vec.size() != 1) emsg("error vec num");
			p_prop = vec[0];
			ntr = 0; nac = 0; nfa = 0;
			
			auto N = model.species[p_prop].nindividual_in;
			ind_sim_prob.resize(N);
			for(auto &va : ind_sim_prob){
				va.ntr = PROP_SIM_NAC_INIT; va.nac = PROP_SIM_NAC_INIT; va.prob = 1;
				va.done = false;
			}
		}
		break;
		
	case IND_OBS_SAMP_PROP:                   // Individual proposals
	case IND_OBS_RESIM_SINGLE_PROP:  
		{
			switch(type){
			case IND_OBS_SAMP_PROP: name += "IND_OBS_SAMP_PROP"; break;
			case IND_OBS_RESIM_SINGLE_PROP: name += "IND_OBS_RESIM_SINGLE_PROP"; break;
			default: emsg("op prob"); break;
			}
			
			if(vec.size() != 1) emsg("error vec num");
			p_prop = vec[0];
			ntr = 0; nac = 0; nfa = 0;
			
			const auto &sp = model.species[p_prop];
			auto N = sp.nindividual_in;
			auto K = sp.ncla;
			ind_sim_prob_cl.resize(N);
			for(auto i = 0u; i < N; i++){
				ind_sim_prob_cl[i].resize(K);
				for(auto cl = 0u; cl < K; cl++){
					auto &va = ind_sim_prob_cl[i][cl];
					va.ntr = PROP_SIM_NAC_INIT; va.nac = PROP_SIM_NAC_INIT; va.prob = 1;
					va.done = false;
				}
			}
		}
		break;
		
	case IND_UNOBS_RESIM_PROP:         
	case IND_ADD_REM_PROP:    
	case IND_ADD_REM_TT_PROP:  	
		switch(type){
		case IND_UNOBS_RESIM_PROP: name += "IND_UNOBS_RESIM_PROP"; break;
		case IND_ADD_REM_PROP: name += "IND_ADD_REM_PROP"; break;
		case IND_ADD_REM_TT_PROP: name += "IND_ADD_REM_TT_PROP"; break;
		default: emsg("op prob"); break;
		}
			
		if(vec.size() != 1) emsg("error vec num");
		p_prop = vec[0];
		ntr = 0; nac = 0; nfa = 0; si = 5; 
		break;
		
	case MBPII_PROP:
		name += "MBPII_PROP";
		if(vec.size() != 1) emsg("error vec num");
		p_prop = vec[0];
		ntr = 0; nac = 0; si = 1;
		mbp_population_affect();
		break;
		
	case MBP_IC_POP_PROP:
		name += "MBP_IC_POP_PROP";
		if(vec.size() != 3) emsg("error vec num");
		p_prop = vec[0];
		cl_prop = vec[1];
		c_prop = vec[2];
		ntr = 0; nac = 0; si = 10;
		mbp_population_affect();
		break;
		
	case MBP_IC_POPTOTAL_PROP:
		name += "MBP_IC_POPTOTAL_PROP";
		if(vec.size() != 1) emsg("error vec num");
		p_prop = vec[0];
		cl_prop = UNSET;
		c_prop = UNSET;
		ntr = 0; nac = 0; si = 10;
		mbp_population_affect();
		break;
		
	case MBP_IC_RESAMP_PROP:
		name += "MBP_IC_RESAMP_PROP";
		if(vec.size() != 1) emsg("error vec num");
		p_prop = vec[0];
		ntr = 0; nac = 0; si = 1;
		mbp_population_affect();
		break;
		
	case INIT_COND_FRAC_PROP:
		name += "INIT_COND_FRAC_PROP";
		if(vec.size() != 1) emsg("error vec num");
		p_prop = vec[0];
		ntr = 0; nac = 0; 
		break;
		
	case IE_PROP:                            // Individual proposals
		name += "IE_PROP";
		if(vec.size() != 2) emsg("error vec num");
		p_prop = vec[0];
		ie_prop = vec[1];
		ntr = 0; nac = 0; si = 1;
		break;
		
	case POP_ADD_REM_LOCAL_PROP:
		name += "POP_ADD_REM_LOCAL_PROP";
		loc_samp.tr_list = vec; loc_samp.tr_list.pop_back();
		p_prop = vec[vec.size()-1];
		ntr = 0; nac = 0;
		break;
		
	case POP_MOVE_LOCAL_PROP:
		name += "POP_MOVE_LOCAL_PROP";
		if(vec.size() != 1) emsg("error vec num");
		p_prop = vec[0];
		ntr = 0; nac = 0;
		break;
	
	case POP_IC_LOCAL_PROP:
		{
			name += "POP_IC_LOCAL_PROP";
			if(vec.size() != 1) emsg("error vec num");
			p_prop = vec[0];
			ntr = 0; nac = 0;
			auto N = model.species[p_prop].tra_gl.size();
			loc_samp.tr_ic.init(N);
		}
		break;
		
	case POP_END_LOCAL_PROP:
		{
			name += "POP_END_LOCAL_PROP";
			if(vec.size() != 1) emsg("error vec num");
			p_prop = vec[0];
			ntr = 0; nac = 0;
			auto N = model.species[p_prop].tra_gl.size();
			loc_samp.tr_ic.init(N);
		}
		break;
		
	case POP_SINGLE_LOCAL_PROP:
		{
			name += "POP_SINGLE_LOCAL_PROP";
			if(vec.size() != 1) emsg("error vec num");
			p_prop = vec[0];
			ntr = 0; nac = 0;
			auto T = model.ntimepoint-1;
			auto ntr = model.species[p_prop].tra_gl.size();
			loc_samp.tr_samp.init(ntr,T);
		}
		break;
		
	case POP_IC_PROP:
		{
			name += "POP_IC_PROP";
			if(vec.size() != 1) emsg("error vec num");
			p_prop = vec[0];
			ntr = 0; nac = 0;
			auto C = model.species[p_prop].comp_gl.size();
			loc_samp.comp_ic.init(C);
		}
		break;
		
	case POP_IC_SWAP_PROP:
		{
			name += "POP_IC_SWAP_PROP";
			if(vec.size() != 1) emsg("error vec num");
			p_prop = vec[0];
			ntr = 0; nac = 0;
			auto C = model.species[p_prop].comp_gl.size();
			loc_samp.comp_ic.init(C);
		}
		break;
	}
	
	plot_index = 0;
	
	set_omega_check();
}


/// Initialises basic sampling variables
void Proposal::initialise_variables()
{
	N = param_list.size();
	
	name += " for: ";
	for(auto j = 0u; j < N; j++){
		if(j != 0) name += ",";
		name += model.param_vec[param_list[j]].name;
	}
	
	ntr = 0; nac = 0; si = 1;
	M.resize(N);
	for(auto j = 0u; j < N; j++) M[j].resize(N,UNSET);
}


/// Prints information about the sampler
string Proposal::print_info() const
{
	stringstream ss;

	switch(type){
	case TRANS_TREE_PROP:
		ss << "TRANSMISSION TREE SAMPLER";	
		break;
		
	case TRANS_TREE_SWAP_INF_PROP:
		ss << "TRANSMISSION TREE SWAP INF SAMPLER";	
		break;
		
	case TRANS_TREE_MUT_PROP:
		ss << "TRANSMISSION TREE MUTATION NUMBER SAMPLER";	
		break;
		
	case TRANS_TREE_MUT_LOCAL_PROP:
		ss << "TRANSMISSION TREE LOCAL MUTATION NUMBER SAMPLER";	
		break;
		
	case PARAM_PROP: case MBP_PROP: case IE_VAR_PROP: case IE_COVAR_PROP: case IE_VAR_CV_PROP: 
		{
			if(type == PARAM_PROP) ss << "PARAMETER SAMPLER for ";
			if(type == MBP_PROP) ss << "MBP SAMPLER for ";
			if(type == IE_VAR_PROP) ss << "IE VAR SAMPLER for ";
			if(type == IE_COVAR_PROP) ss << "IE COVAR SAMPLER for ";
			if(type == IE_VAR_CV_PROP) ss << "IE VAR CV SAMPLER for ";
			
			{
				string str = ""; 
				for(auto th : param_list) str +=  model.param_vec[th].name+", ";
				ss << "param: " << trunc(str) << endl;
			}
			
			{
				string str = ""; 
				for(auto th : dependent) str += model.param_vec[th].name+", ";
				ss << "dependent" << trunc(str) << endl;
			}
		}
		break;
		
	case PAR_EVENT_FORWARD_PROP:
		{
			ss << "INDIVIDUAL forward joint parameter event for species " << model.species[p_prop].name << endl;
		}
		break;
		
	case PAR_EVENT_FORWARD_SQ_PROP:
		{
			ss << "INDIVIDUAL forward squash joint parameter event for species " << model.species[p_prop].name << endl;
		}
		break;
		
	case PAR_EVENT_BACKWARD_SQ_PROP:
		{
			ss << "INDIVIDUAL backward squash joint parameter event for species " << model.species[p_prop].name << endl;
		}
		break;
		
	case IND_EVENT_TIME_PROP:
		{
			ss << "INDIVIDUAL event time sampler for species " << model.species[p_prop].name << endl;
		}
		break;
		
	case IND_MULTI_EVENT_PROP:
		{
			ss << "INDIVIDUAL multi event time sampler for species " << model.species[p_prop].name << endl;
		}
		break;
		
	case IND_EVENT_ALL_PROP:
		{
			ss << "INDIVIDUAL all event time sampler for species " << model.species[p_prop].name << endl;
		}
		break;
		
	case IND_OBS_SWITCH_ENTER_SOURCE_PROP:
		ss << "INDIVIDUAL correct switch enter/source propsal for species " << model.species[p_prop].name << endl;
		break;
		
	case IND_OBS_SWITCH_LEAVE_SINK_PROP:
		ss << "INDIVIDUAL correct switch leave/sink propsal for species " << model.species[p_prop].name << endl;
		break;
	
	case CORRECT_OBS_TRANS_PROP:
		ss << "INDIVIDUAL correct obs trans propsal for species " << model.species[p_prop].name << endl;
		break;
		
	case IND_LOCAL_PROP:
		{
			ss << "INDIVIDUAL local propsal for species " << model.species[p_prop].name 
					<< " classification " << model.species[p_prop].cla[cl_prop].name << endl;
		}
		break;
		
	case IND_OBS_SAMP_PROP:
		{
			ss << "INDIVIDUAL sampler for species (observed) " << model.species[p_prop].name << endl;
		}
		break;
		
	case IND_OBS_RESIM_PROP:
		{
			ss << "INDIVIDUAL resimulate for species (observed) " << model.species[p_prop].name << endl;
		}
		break;
		
	case IND_OBS_RESIM_SINGLE_PROP:
		{
			ss << "INDIVIDUAL resimulate single cl for species (observed) " << model.species[p_prop].name << endl;
		}
		break;
		
	case IND_UNOBS_RESIM_PROP:
		{
			ss << "INDIVIDUAL resimulate for species (unobserved) " << model.species[p_prop].name << endl;
		}
		break;
	
	case IND_ADD_REM_PROP:
		{
			ss << "INDIVIDUAL add/remove (unobserved) " << model.species[p_prop].name << endl;
		}
		break;
		
	case IND_ADD_REM_TT_PROP:
		{
			ss << "INDIVIDUAL add/remove with trans-tree (unobserved) " << model.species[p_prop].name << endl;
		}
		break;
		
	case MBPII_PROP:
		{
			ss << "MBP TYPE II for species "  << model.species[p_prop].name << endl;
		}
		break;
	
	case MBP_IC_POP_PROP:
		{
			const auto &sp = model.species[p_prop];
			ss << "MBP for initial population for species "  << sp.name;
			ss << " for " << sp.cla[cl_prop].name << " compartment " << sp.cla[cl_prop].comp[c_prop].name << endl;
		}
		break;
		
	case MBP_IC_POPTOTAL_PROP:
		{
			const auto &sp = model.species[p_prop];
			ss << "MBP for total initial population for species "  << sp.name << endl;
		}
		break;
		
	case MBP_IC_RESAMP_PROP:
		{
			const auto &sp = model.species[p_prop];
			ss << "MBP resampling init_cond for "  << sp.name << endl;
		}
		break;
		
		case INIT_COND_FRAC_PROP:
		{
			const auto &sp = model.species[p_prop];
			ss << "Gibbs proposal for init_cond fractions "  << sp.name << endl;
		}
		break;
	
	case IE_PROP:
		{
			ss << "IE sampler for " << model.species[p_prop].ind_effect[ie_prop].name << endl;
		}
		break;
		
	case POP_ADD_REM_LOCAL_PROP:
		{
			ss << "LOCAL sampler add remove evnets for " << model.species[p_prop].name << endl;
		}
		break;
		
	case POP_MOVE_LOCAL_PROP:
		{
			ss << "LOCAL sampler move events for " << model.species[p_prop].name << endl;
		}
		break;
		
	case POP_IC_LOCAL_PROP:
		{
			ss << "LOCAL sampler init cond for " << model.species[p_prop].name << endl;
		}
		break;
		
	case POP_END_LOCAL_PROP:
		{
			ss << "LOCAL sampler end for " << model.species[p_prop].name << endl;
		}
		break;
		
	case POP_SINGLE_LOCAL_PROP:
		{
			ss << "LOCAL sampler single events for " << model.species[p_prop].name << endl;
		}
		break;
		
	case POP_IC_PROP:
		{
			ss << "GLOBAL sampler init cond for " << model.species[p_prop].name << endl;
		}
		break;
		
	case POP_IC_SWAP_PROP:
		{
			ss << "GLOBAL sampler init cond swap for " << model.species[p_prop].name << endl;
		}
		break;
	}

	if(1 == 0) cout << "Affect like output is turned off in code" << endl;		
	else{
		for(const auto &al : affect_like){
			ss << output.print_affect_like(al) << endl;
		}
	}
	
	if(!on) ss << "Proposal switched off" << endl;
	
	ss << endl;
	ss << endl;
	
	return ss.str();
}


/// Performs a Metropolis-Hastings proposal update
void Proposal::MH(State &state)
{
	auto pl = false;
	
	auto param_store = state.param_val;
	
	state.param_val = sample(state.param_val);

	if(false){
		for(auto j = 0u; j < model.param_vec.size(); j++){
			cout << name << endl;
			
			const auto &pv = model.param_vec[j];
			if(param_store[j] != state.param_val[j]){
				cout << pv.name << " " << param_store[j] << " -> " << state.param_val[j] << endl;
			}
		}
	}
	
	ntr++; 
	
	if(model.inbounds(state.param_val) == false || 
		(omega_check && model.ie_cholesky_error(state.param_val))){ 
		state.param_val = param_store; 
		update_si(-0.005);
		return;
	}

	state.update_spline(affect_spline);
	
	auto like_ch = state.update_param(affect_like,param_store);
	
	auto al = calculate_al(like_ch,0);
	
	if(pl) cout << al << "al" << endl;
	if(ran() < al){ 
		if(pl) cout << "accept" << endl;
		nac++;
			
		state.add_like(like_ch);
		update_si(0.01);
		state.remove_store_spline(affect_spline);
	}
	else{ 
		if(pl) cout << "reject" << endl;
		state.param_val = param_store;
		state.restore(affect_like);
		state.restore_spline(affect_spline);
		
		update_si(-0.005); if(si > 3) si *= 0.75;
	}

	if(pl) state.check("ev");
}


/// Performs a MBP proposal update
void Proposal::mbp(State &state)
{
	auto pl = false;
	
	auto param_store = state.param_val;
	InitCondValue init_cond_store;
	
	auto &ssp = state.species[p_prop];
	
	ntr++;
	
	auto popnum_t_st = state.popnum_t;

	auto sim_prob = 0.0;
	switch(type){
	case MBP_PROP: 
		state.param_val = sample(state.param_val); 

		if(model.inbounds(state.param_val) == false){
			state.param_val = param_store; 
			update_si(-0.005);
			return;
		}
		break;
		
	case MBPII_PROP:
		sim_prob = si; 
		break;
		
	case MBP_IC_POP_PROP: case MBP_IC_POPTOTAL_PROP: case MBP_IC_RESAMP_PROP:
		{
			init_cond_store = ssp.init_cond_val;

			auto res = propose_init_cond(ssp.init_cond_val,state);
			
			if(res != IC_SUCCESS){
				ssp.init_cond_val = init_cond_store;
				if(res != IC_ZERO){
					update_si(-0.005);
					if(type == MBP_IC_POP_PROP || type == MBP_IC_POPTOTAL_PROP){
						if(si < POP_SI_LIM) si = POP_SI_LIM;
					}
				}
				return;
			}
			
			state.update_popnum_t_init(p_prop,init_cond_store.cnum,ssp.init_cond_val.cnum);
		}
		break;
		
	default: emsg("op err"); break;
	}
	
	if(ssp.type != POPULATION) emsg("must be population");
	
	auto trans_num_st = ssp.trans_num;
	auto tnum_mean_st_st = ssp.tnum_mean_st;
	auto cpop_st_st = ssp.cpop_st;
	
	state.update_spline(affect_spline);

	ssp.mbp(sim_prob,state.popnum_t);

	auto like_ch = state.update_param(affect_like,param_store);
	
	auto al = calculate_al(like_ch,0);
	
	if(pl) cout << al << " al" << endl; 
	
	if(ran() < al){
		if(pl) cout << "accept" << endl;
		nac++;
		
		switch(type){
		case MBP_IC_POP_PROP: case MBP_IC_POPTOTAL_PROP: case MBP_IC_RESAMP_PROP:
			//state.likelihood_init_cond(p_prop,like_ch.init_cond); CHA
			ssp.likelihood_init_cond(like_ch.init_cond);
			break;
		default: break;
		}
		
		state.add_like(like_ch);
		ssp.mbp_accept(state.like.markov,state.popnum_t);
	
		update_si(0.005); 
		if(si > 1 && (type == MBPII_PROP || type == MBP_IC_RESAMP_PROP)) si = 1;
		state.remove_store_spline(affect_spline);
	}
	else{ 
		if(pl) cout << "reject" << endl;
		state.param_val = param_store;
		state.popnum_t = popnum_t_st;
		ssp.trans_num = trans_num_st;
		ssp.tnum_mean_st = tnum_mean_st_st;
		ssp.cpop_st = cpop_st_st;
		switch(type){
		case MBP_IC_POP_PROP: case MBP_IC_POPTOTAL_PROP: case MBP_IC_RESAMP_PROP:
			ssp.init_cond_val = init_cond_store;
			break;
		default: break;
		}
		
		state.restore(affect_like);		
		state.restore_spline(affect_spline);
		update_si(-0.005);
		if(type == MBP_IC_POP_PROP || type == MBP_IC_POPTOTAL_PROP){
			if(si < POP_SI_LIM) si = POP_SI_LIM;
		}
	}
}


/// Performs a Metropolis-Hastings proposal update for an individual event time
void Proposal::MH_event(State &state)
{
	if(skip_proposal(0.5)) return;

	auto pl = false;

	auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	
	auto t_start = model.details.t_start;
	auto t_end = model.details.t_end;

	for(auto i = 0u; i < ssp.individual.size(); i++){	
		if(i >= sp.nindividual_obs || sp.individual[i].move_needed){
			auto &ind = ssp.individual[i];
		
			auto loop_max = ind.ev.size();
		
			for(auto loop = 0u; loop < loop_max; loop++){
				auto &event = ind.ev;
			
				auto e = (unsigned int)(ran()*event.size());
				if(e == event.size()) emsg("event not in range");
				auto &ev = event[e];
				
				if(ev.observed == false){
					switch(ev.type){
					case NM_TRANS_EV: case M_TRANS_EV:
						if(!(sp.period_exist && ev.c_after != UNSET && sp.comp_period[ev.c_after][ev.cl])){							
							auto t = ev.t;
					
							auto tmin = t_start; 
							if(ev.c_after == UNSET) tmin = event[e-1].t;
							else{
								auto e_before = get_event_before(e,event);
								if(e_before != UNSET) tmin = event[e_before].t;
							}
							
							auto tmax = t_end;
							if(e == 0){
								if(event.size() > 1) tmax = event[1].t;
								if(sp.tra_gl[ ev.tr_gl].variety != SOURCE_TRANS) emsg("Should be source trans");
							}
							else{	
								auto e_after = get_event_after(e,event);
								if(e_after != UNSET) tmax = event[e_after].t;
							}
							
							auto tr_gl = ev.tr_gl;
							auto &samp = ind_sampler[tr_gl];
							
							string warn;
							auto dt = normal_sample(0,samp.si,warn);
							if(dt == UNSET) emsg("Problem with MCMC proposal: "+warn);
							
							auto t_new = t + dt;
							
							samp.ntr++;
							if(t_new <= tmin || t_new >= tmax || event_near_div(t_new,model.details)){
								samp.nfa++;
								update_ind_samp_si(tr_gl,-0.005);
							}
							else{
								auto ev_new = ind.ev;
									
								auto tr_gl_new = ssp.move_event(ev_new,e,t_new);
		
								if(pl){
									cout << endl << endl;
									cout << t << " -> " << t_new << " time change" << endl;
									cout << ind.name << endl;
									cout << "event old:" << endl; ssp.print_event(ind.ev);
									cout << "event new:" << endl; ssp.print_event(ev_new);
								}
								
								auto gc = state.update_tree(p_prop,i,ev_new);
								if(gc.type == GENCHA_FAIL){
									samp.nfa++;
									update_ind_samp_si(tr_gl,-0.005);
								}
								else{
									auto like_ch = state.update_ind(p_prop,i,ev_new,UP_SINGLE);
								
									auto dprob = 0.0;
									if(tr_gl_new != UNSET){  // Accounts for differeces in time sampler
										const auto &samp_new = ind_sampler[tr_gl_new];
										auto dt = t_new-t;
										dprob = normal_probability(dt,0,samp_new.si) - normal_probability(dt,0,samp.si); 
									}
									
									gc.update_like_ch(like_ch,dprob);
									
									auto al = calculate_al(like_ch,dprob);
							
									if(pl) cout << i << " " <<  al << " " << e << " al" << endl;
								
									if(ran() < al){
										if(pl) cout << " accept" << endl;
										samp.nac++;
										state.add_like(like_ch);
										state.gen_change_update(gc); 
										if(sp.trans_tree) state.update_popnum_ind(p_prop,i);
										
										update_ind_samp_si(tr_gl,0.005);
									}
									else{ 
										if(pl) cout << "reject" << endl;
										state.restore_back();
										update_ind_samp_si(tr_gl,-0.005);
									}

									if(pl) state.check("ind prop");
	
									if(pl){
										state.check_popnum_t2("Pnum");
									}
								}
							}
						}
						break;
						
					default: break;
					}
				}
			}
		}
	}
	if(pl) cout << " done" << endl;
}


/// Performs a Metropolis-Hastings proposal update for multiple individual event times
void Proposal::MH_multi_event(State &state)
{
	if(skip_proposal(0.5)) return;

	auto pl = false;

	auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	
	for(auto i = 0u; i < ssp.individual.size(); i++){	
		if(i >= sp.nindividual_obs || sp.individual[i].move_needed){
			auto &ind = ssp.individual[i];
		
			auto loop_max = ind.ev.size();
			for(auto loop = 0u; loop < loop_max; loop++){
				auto &event = ind.ev;
				auto E = event.size();
			
				// Samples to root event
				auto e = (unsigned int)(ran()*E);
				if(e == E) emsg("event not in range");
				
				auto lm = ssp.get_ev_link_listmove(e,event,ind,state.popnum_t);
				if(lm.list.size() > 1){
					auto trange = ssp.ev_link_trange(lm,event);
					
					auto tr_gl = event[lm.list[0]].tr_gl;
				
					auto &samp = ind_sampler[tr_gl];
					
					vector <EventMove> evmo;
					
					auto ill = false;
					
					string warn;
					auto dt = normal_sample(0,samp.si,warn);
					if(dt == UNSET) emsg("Problem with MCMC proposal: "+warn);
							
					if(dt >= trange.tmax || dt < trange.tmin) ill = true;
					else{
						// Reorders events 
						for(auto e = 0u; e < E; e++){
							auto tnew = event[e].t;
							if(lm.move[e]){
								tnew += dt;
								if(event_near_div(tnew,model.details)) ill = true;
							}
							
							EventMove em;
							em.t = tnew;
							em.e = e;
							em.move = lm.move[e];
							evmo.push_back(em);
						}
					}
					
					vector <Event> ev_new;
					unsigned int tr_gl_new;
					if(!ill){
						vector <bool> move_new;
						for(auto e = 0u; e < E; e++){
							const auto &evm = evmo[e];
							auto ev = event[evm.e];
							ev.t = evm.t;
							ev_new.push_back(ev);
							move_new.push_back(evm.move);
						}
					
						ssp.ensure_consistent(ev_new);
					
						if(pl){
							cout << endl << endl;
							cout << ind.name << endl;
							cout << "event old:" << endl; ssp.print_event(ind.ev);
							cout << "event new:" << endl; ssp.print_event(ev_new);
						}
						
						auto e_new = 0u; while(e_new < E && move_new[e_new] == false) e_new++;
						if(e_new == E) emsg("Should not be E");
						
						auto lm_new = ssp.get_ev_link_listmove(e_new,ev_new,ind,state.popnum_t);
						
						if(lm.list.size() != lm_new.list.size())  ill = true;
						else{
							if(!equal_vec(lm_new.move,move_new)) ill = true;
							else{
								tr_gl_new = event[lm_new.list[0]].tr_gl;
							}
						}
					}
					
					samp.ntr++;
					if(ill){
						samp.nfa++;
						update_ind_samp_si(tr_gl,-0.005);
					}
					else{		
						auto gc = state.update_tree(p_prop,i,ev_new);
						if(gc.type == GENCHA_FAIL){
							samp.nfa++;
							update_ind_samp_si(tr_gl,-0.005);
						}
						else{
							auto like_ch = state.update_ind(p_prop,i,ev_new,UP_SINGLE);
						
							auto dprob = 0.0;
							if(tr_gl_new != tr_gl){  // Accounts for differeces in time sampler
								const auto &samp_new = ind_sampler[tr_gl_new];
								dprob = normal_probability(dt,0,samp_new.si) - normal_probability(dt,0,samp.si); 
							}
							
							gc.update_like_ch(like_ch,dprob);
							
							auto al = calculate_al(like_ch,dprob);
					
							if(pl) cout << i << " " <<  al << " " << e << " al" << endl;
						
							if(ran() < al){
								if(pl) cout << " accept" << endl;
								samp.nac++;
								state.add_like(like_ch);
								state.gen_change_update(gc); 
								if(sp.trans_tree) state.update_popnum_ind(p_prop,i);
								
								update_ind_samp_si(tr_gl,0.005);
							}
							else{ 
								if(pl) cout << "reject" << endl;
								state.restore_back();
								update_ind_samp_si(tr_gl,-0.005);
							}

							if(pl) state.check("ind multi event");
						}
					}
				}
			}
		}
	}
	if(pl) cout << " done" << endl;
}
	

/// Performs a Metropolis-Hastings proposal update for all individual event times
void Proposal::MH_event_all(State &state)
{
	if(skip_proposal(0.5)) return;

	auto pl = false;

	auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	
	auto t_start = model.details.t_start;
	auto t_end = model.details.t_end;

	for(auto i = 0u; i < sp.nindividual_obs; i++){	
		if(sp.individual[i].move_needed){
			auto &ind = ssp.individual[i];
		
			auto &event = ind.ev;
			auto E = event.size();
			
			auto emin = 0u; auto emax = E;
			
			auto tmin = t_start;
			auto tmax = t_end;
			
			if(event[0].type == ENTER_EV){ emin = 1; tmin = event[0].t;}
			if(event[E-1].type == LEAVE_EV){ emax = E-1; tmax = event[E-1].t;}
			
			auto ill = false;
			for(auto e = emin; e < emax; e++){
				if(event[e].observed) ill = true;
			}
			
			if(ill == false){
				auto &samp = ind_sampler[i];
				
				string warn;
				auto dt = normal_sample(0,samp.si,warn);
				if(dt == UNSET) emsg("Problem with MCMC proposal: "+warn);
				
				if(dt > 0){
					if(event[emax-1].t+dt > tmax) ill = true;
				}
				else{
					if(event[emin].t+dt < tmin) ill = true;
				}
			
				samp.ntr++;
				if(ill){
					samp.nfa++;
					update_ind_samp_si(i,-0.005);
				}
				else{
					auto ev_new = event;
					for(auto e = emin; e < emax; e++){
						ev_new[e].t += dt;
					}
					
					if(!events_near_div(ev_new,model.details)){
						if(pl){
							cout << endl << endl;
							cout << ind.name << endl;
							cout << "event old:" << endl; ssp.print_event(ind.ev);
							cout << "event new:" << endl; ssp.print_event(ev_new);
						}
						
						auto gc = state.update_tree(p_prop,i,ev_new);
						if(gc.type == GENCHA_FAIL){
							samp.nfa++;
							update_ind_samp_si(i,-0.005);
						}
						else{
							auto like_ch = state.update_ind(p_prop,i,ev_new,UP_SINGLE);
						
							auto dprob = 0.0;
							gc.update_like_ch(like_ch,dprob);
							
							auto al = calculate_al(like_ch,dprob);
					
							if(pl) cout << i << " " <<  al << " al" << endl;
						
							if(ran() < al){
								if(pl) cout << " accept" << endl;
								samp.nac++;
								state.add_like(like_ch);
								state.gen_change_update(gc); 
								if(sp.trans_tree) state.update_popnum_ind(p_prop,i);
								
								update_ind_samp_si(i,0.005);
							}
							else{ 
								if(pl) cout << "reject" << endl;
								state.restore_back();
								update_ind_samp_si(i,-0.005);
							}

							if(pl) state.check("ind event all prop");

							if(pl){
								state.check_popnum_t2("Pnum");
							}
						}
					}
				}
			}
		}
	}
	if(pl) cout << " done" << endl;
}


/// Performs Metropolis-Hastings local updates for individuals 
void Proposal::MH_ind_local(State &state)
{
	auto pl = false;     
			
	auto &sp = model.species[p_prop];
	auto cl = cl_prop;
		
	auto &ssp = state.species[p_prop];
		
	const auto &individual = ssp.individual;
		
	vector <LocalIndChange> licha;
	
	vector < vector <unsigned int> > lc_ref;
	lc_ref.resize(individual.size());

	// Adds all possible local changes
	for(auto i = 0u; i < individual.size(); i++){
		auto add = ssp.local_ind_change(i,cl);
		add_li_cha(i,add,licha,lc_ref);
	}
	
	if(false){
		for(auto &lich : licha) ssp.print_local_ind_change(lich);
		emsg("do");	
	}
		
	const vector <Event> empty_event;           

	if(licha.size() > 0){
		auto loopmax = (unsigned int)(licha.size());
		if(loopmax > individual.size()) loopmax= individual.size();
		if(loopmax == 0) loopmax = 1;

		Event enew;
		enew.move_c = UNSET;
		enew.cl = cl;
		enew.observed = false;
		enew.inf_node_ref = UNSET;
		enew.ind_inf_from = IndInfFrom();	

		auto nind_obs = sp.nindividual_obs;
	
		for(auto loop = 0u; loop < loopmax; loop++){
			auto index = (unsigned int)(ran()*licha.size());

			if(pl) cout << endl << endl << "START" << endl;
			
			const auto &lich = licha[index];
			
			vector <Event> ev_new;
			
			if(pl){
				ssp.print_local_ind_change(lich);
				cout << "event older:" << endl; ssp.print_event(ssp.individual[lich.i].ev);
			}
			//cout << "event older:" << endl; ssp.print_event(ssp.individual[lich.i].ev);
			//ssp.print_local_ind_change(lich);
		
			auto timefac = 1.0;
		
			auto prob_try = ssp.create_local_change(timefac,lich,ev_new,enew,nind_obs,empty_event,LOCAL_FORWARD);
		
			//prob_try = 1.0/timefac;
	
			auto &si = swap_info[lich.swap_rep_ref];
			si.ntr++;
		
			if(prob_try == 0){
				si.nzero++;
			}
			else{
				if(ran() > prob_try){
					si.nfilt++;
				}
				else{
					si.ntr2++;
					
					auto ev_store = ssp.individual[lich.i].ev;
					if(pl){
						cout << endl << endl << "Local change:" << endl;
						ssp.print_local_ind_change(lich);
				
						cout << "event old:" << endl; ssp.print_event(ssp.individual[lich.i].ev);
						cout << "event new:" << endl; ssp.print_event(ev_new);
					}
				
					if(!events_near_div(ev_new,model.details)){
						auto i = lich.i;
					
						auto gc = state.update_tree(p_prop,i,ev_new);
						if(gc.type == GENCHA_FAIL){
							si.nfa++;
						}
						else{
							auto like_ch = state.update_ind(p_prop,i,ev_new,UP_SINGLE);
						
							auto add = ssp.local_ind_change(i,cl);
						
							auto dprob = 0.0;
							gc.update_like_ch(like_ch,dprob);
							
							auto al = calculate_al(like_ch,dprob);
							
							// Accounts for change in number of potential proposals
							auto Ni = licha.size();
							auto Nf = Ni + add.size() - lc_ref[i].size();
							if(Nf != Ni) al *= double(Ni)/Nf; 
											
							al *= timefac/prob_try;
							
							if(pl) cout << al << " al" << endl;
							
							if(ran() < al){ 
								if(pl) cout << "accept" << endl;
							
								if(testing) ssp.swap_check_rev_prob_try(prob_try,timefac,lich,add,ev_store);
								si.nac++;
								remove_li_cha(i,licha,lc_ref);
								add_li_cha(i,add,licha,lc_ref);
								state.add_like(like_ch);
								state.gen_change_update(gc); 
								if(sp.trans_tree) state.update_popnum_ind(p_prop,i);
								if(testing){ if(Nf != licha.size()) emsg("Wrong number");}
							}
							else{
								if(pl) cout << "reject" << endl;
								state.restore_back();
							}
							
							if(pl) state.check("ind prop");
						}
					}		
				}
			}
		}		
	}

	if(false){  // Used for testing
		for(auto j = 0u; j < licha.size(); j++){
			const auto &lic = licha[j];
			if(lc_ref[lic.i][lic.ref] != j) emsg("Problem with reference"); 
		}
		
		for(auto i = 0u; i < individual.size(); i++){
			for(auto j : lc_ref[i]){
				if(licha[j].i != i) emsg("Problem with reference i");
			}
		}
		
		auto lc_ref_st = lc_ref;
		auto licha_st = licha;
	
		for(auto i = 0u; i < individual.size(); i++) lc_ref[i].clear();
		licha.clear();
		
		for(auto i = 0u; i < individual.size(); i++){	
			auto add = ssp.local_ind_change(i,cl);
			add_li_cha(i,add,licha,lc_ref);
		}
		
		if(licha.size() != licha_st.size()) emsg("licha prob");
		for(auto i = 0u; i < individual.size(); i++){
			if(lc_ref[i].size() != lc_ref_st[i].size()) emsg("lc_ref problem");
		}
		
		state.check("  ind prop");
	}
}

	
/// Performs a transmission tree proposal
void Proposal::trans_tree(State &state)
{
	if(skip_proposal(0.9)) return;
	
	state.trans_tree_proposal(burn_info,nac,ntr);
}


/// Performs a transmission tree proposal swapping infectee and infector
void Proposal::trans_tree_swap_inf(State &state)
{
	if(skip_proposal(0.9)) return;
	
	state.trans_tree_swap_inf_proposal(burn_info,nfa,nac,ntr);
}


/// Performs a transmission tree mutation proposal (simple change in mutation number on a link)
void Proposal::trans_tree_mut(State &state)
{
	if(skip_proposal(0.9)) return;
	
	state.trans_tree_mut_proposal(burn_info,nac,ntr,si);
}


/// Performs a local tree proposal (shifts events around a branching point
void Proposal::trans_tree_mut_local(State &state)
{
	if(skip_proposal(0.9)) return;
	
	state.trans_tree_mut_local_proposal(burn_info,gen_mut_info);
}


/// Performs a sampling proposal update for an observed individual 
void Proposal::sample_ind_obs(State &state)
{
	if(skip_proposal(0.8)) return;
	
	auto pl = false;
	
	auto burn = burn_info.on;

	if(pl) state.check("before");
	
	auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	
	IndEvSampler ind_ev_samp(ssp.markov_eqn_vari,ssp.individual,model.details,sp,ssp.obs_eqn_value,ssp.obs_trans_eqn_value,model.eqn,state.genetic_value.inf_node,state.param_val,state.spline_val,state.popnum_t);
	
	ind_ev_samp.setup_nm();
	
	for(auto i = 0u; i < sp.nindividual_in; i++){
		auto &ind = ssp.individual[i];
		
		if(sp.individual[i].simulation_needed){
			for(auto cl = 0u; cl < sp.ncla; cl++){
				if(sp.individual[i].sample_needed[cl] || ssp.all_events_correct(i,cl) == false){
					if(pl) cout << endl << endl << i << "update    cl=" << cl << endl;
			
					if(ind_ev_samp.needed(i,cl) == true){
						auto &isp = ind_sim_prob_cl[i][cl];
					
						if(ran() < isp.prob){
							ntr++; if(burn) isp.ntr++;
							ind_ev_samp.generate_ind_obs_timeline();
							
							auto probif = 0.0;
						
							auto ev_new = ind_ev_samp.sample_events(probif);
							
							if(pl){	
								cout << i << " " << ind.name << "  individual" << endl;
								cout << "old: "; ssp.print_event(ind.ev);
								cout << "new: "; ssp.print_event(ev_new);
								cout << ind_ev_samp.illegal << "illegal" << endl;
								if(ind_ev_samp.illegal) emsg("stop");
							}
								
							if(ind_ev_samp.illegal == true) nfa++;
							else{
								if(testing == true){ // Diagnostic
									auto prob = ind_ev_samp.sample_events_prob(ev_new);
									if(dif(probif,prob,DIF_THRESH)) emsg("sampler different here");
								}
								
								auto ev_store = ind.ev;
								
								auto gc = state.update_tree(p_prop,i,ev_new);
								
								if(gc.type == GENCHA_FAIL) nfa++;
								else{
									auto like_ch = state.update_ind(p_prop,i,ev_new,UP_SINGLE);
								
									if(pl){
										cout << "modi: "; ssp.print_event(ev_new);
									}
									
									ind_ev_samp.generate_ind_obs_timeline();
									auto probfi = ind_ev_samp.sample_events_prob(ev_store);
									
									auto dprob = probfi-probif;
								
									gc.update_like_ch(like_ch,dprob);
									
									auto al = calculate_al(like_ch,dprob);
								
									if(pl){
										cout << al << " " << dprob <<  "al" << endl;	
										print_like(like_ch);
									}
								
									if(ran() < al){
										if(pl) cout << "ac ind" << endl;
										nac++; if(burn){ if(event_dif(ev_store,ind.ev)) isp.nac++;}
										state.add_like(like_ch);
										state.gen_change_update(gc); 	
										if(sp.trans_tree) state.update_popnum_ind(p_prop,i);
									}
									else{ 
										if(pl) cout << "reject" << endl;
										state.restore_back();
									}
								}
								
								if(pl) state.check(" indupdate");
							}
							if(burn) ind_obs_prob_update(isp);
						}
					}
				}
			}
		}
	}
}


/// Performs a resimulation proposal update for an observed individual 
void Proposal::resimulate_ind_obs(State &state)
{
	//if(skip_proposal(0.5)) return;
	 
	auto pl = false;
	auto burn = burn_info.on;

	auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	
	if(sp.nindividual_in == 0) return;

	IndEvSampler ind_ev_samp(ssp.markov_eqn_vari,ssp.individual,model.details,sp,ssp.obs_eqn_value,ssp.obs_trans_eqn_value,model.eqn,state.genetic_value.inf_node,state.param_val,state.spline_val,state.popnum_t);
		
	for(auto i = 0u; i < sp.nindividual_in; i++){	
		const auto &indd = sp.individual[i];
		
		if(indd.simulation_needed){
			auto &isp = ind_sim_prob[i];
			if(ran() < isp.prob){
				auto &ind = ssp.individual[i];
				
				auto probif = 0.0, probfi = 0.0;
				
				auto e_init = ind_ev_samp.resample_init_event(i,probif);
		
				if(pl){
					cout << endl << endl << ind.name << "ind" << endl; 
					cout << "older: "; ssp.print_event(ind.ev); 
				}
						
				auto ev_new = ind_ev_samp.simulate_events(i,e_init,probif,indd.trig_ev_ref);
		
				if(pl){
					cout << "new: "; ssp.print_event(ind.ev); 	
					if(ind_ev_samp.illegal) emsg("do");
				}
				ntr++; if(burn) isp.ntr++;
				if(ind_ev_samp.illegal){
					nfa++;
				}
				else{
					if(testing == true){ // Diagnostic
						auto prob = ind_ev_samp.resample_init_event_prob(i,e_init) + ind_ev_samp.simulate_events_prob(i,ev_new,indd.trig_ev_ref);
						if(dif(probif,prob,DIF_THRESH)){
							emsg("sampler different3");
						}
					}
					
					auto ev_store = ind.ev;	
					
					auto gc = state.update_tree(p_prop,i,ev_new);
					if(gc.type == GENCHA_FAIL){
						nfa++;
					}
					else{
						if(pl){
							cout << endl << endl << ind.name << "ind " << endl; 
							cout << "old: "; ssp.print_event(ind.ev);
							for(auto i = 0u; i < ind.ev.size(); i++){
								cout << ind.ev[i].t - int(ind.ev[i].t) << ",";
							}
							cout << " dif old " << endl;
						}
						
						probfi += ssp.nm_obs_dprob(ind);
						auto like_ch = state.update_ind(p_prop,i,ev_new,UP_SINGLE);
						probif += ssp.nm_obs_dprob(ind);
						
						if(pl){
							cout << "new: "; ssp.print_event(ind.ev);
						}
						
						probfi += ind_ev_samp.resample_init_event_prob(i,ev_store[0]) + 
											ind_ev_samp.simulate_events_prob(i,ev_store,indd.trig_ev_ref);
				
						auto dprob = probfi-probif;
						
						gc.update_like_ch(like_ch,dprob);
						if(pl) print_like(like_ch);
						
						auto al = calculate_al(like_ch,dprob);

						if(pl) cout << al << " al try ind" << endl;
						if(ran() < al){ 
							if(pl) cout << "ac ind" << endl;
							nac++; if(burn){ if(event_dif(ev_store,ind.ev)) isp.nac++;}
							state.add_like(like_ch);
							state.gen_change_update(gc); 	
							if(sp.trans_tree) state.update_popnum_ind(p_prop,i);
						}
						else{ 
							if(pl) cout << "rej ind" << endl;
							state.restore_back();
						}
					}
					
					if(pl) state.check(" After ind prop");
				}
				
				if(burn) ind_obs_prob_update(isp);
			}
		}
	}
}


/// Performs resimulation proposal update for an observed individual in different classifications
void Proposal::resimulate_single_ind_obs(State &state)
{
	if(skip_proposal(0.5)) return;
	
	auto pl = false;
	auto burn = burn_info.on;

	auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	
	if(sp.nindividual_in == 0) return;

	IndEvSampler ind_ev_samp(ssp.markov_eqn_vari,ssp.individual,model.details,sp,ssp.obs_eqn_value,ssp.obs_trans_eqn_value,model.eqn,state.genetic_value.inf_node,state.param_val,state.spline_val,state.popnum_t);
	
	for(auto i = 0u; i < sp.nindividual_in; i++){	
		const auto &indd = sp.individual[i];
		for(auto cl = 0u; cl < sp.ncla; cl++){
			if((indd.sample_needed[cl] || ssp.all_events_correct(i,cl) == false) && !sp.sink_exist[cl]){		
				auto &isp = ind_sim_prob_cl[i][cl];
				if(ran() < isp.prob){
					auto &ind = ssp.individual[i];
					
					const auto &trig = indd.cl_trig_ev_ref[cl];
					
					auto probif = 0.0, probfi = 0.0;
		
					if(pl){
						cout << endl << " " << endl << ind.name << "ind" << endl; 
						cout << "older: "; ssp.print_event(ind.ev);
						cout << cl << " cl" << endl;
					}
							
					auto ev_new = ind_ev_samp.simulate_single_events(i,cl,probif,trig);
			
					if(pl){
						cout << "newer: "; ssp.print_event(ev_new);
					}
					
					ntr++; if(burn) isp.ntr++;
					if(ind_ev_samp.illegal){
						nfa++;
					}
					else{
						if(testing == true){ // Diagnostic
							auto prob = ind_ev_samp.simulate_single_events_prob(i,cl,ev_new,trig);
							if(dif(probif,prob,DIF_THRESH)){
								cout <<  probif << " " << prob << " prob" << endl;
								emsg("sampler different4");
							}
						}
						
						auto ev_store = ind.ev;	
						
						auto gc = state.update_tree(p_prop,i,ev_new);
						if(gc.type == GENCHA_FAIL){
							nfa++;
						}
						else{
							if(pl){
								cout << endl << " " << endl << ind.name << "ind" << endl; 
								cout << "old: "; ssp.print_event(ind.ev);
								cout << "new: "; ssp.print_event(ev_new);
							}
							
							probfi += ssp.nm_single_obs_dprob(cl,ind);
							auto like_ch = state.update_ind(p_prop,i,ev_new,UP_SINGLE);
							probif += ssp.nm_single_obs_dprob(cl,ind);
							
							probfi += ind_ev_samp.simulate_single_events_prob(i,cl,ev_store,trig);
					
							auto dprob = probfi-probif;
							
							gc.update_like_ch(like_ch,dprob);
										
							auto al = calculate_al(like_ch,dprob);
						
							if(pl) cout << al << " al try ind" << endl;
							if(ran() < al){
								if(pl) cout << "ac ind" << endl;
								nac++; if(burn){ if(event_dif(ev_store,ind.ev)) isp.nac++;}
								state.add_like(like_ch);
								state.gen_change_update(gc); 	
								if(sp.trans_tree) state.update_popnum_ind(p_prop,i);
							}
							else{
								state.restore_back();
							}
						}
						
						if(pl) state.check(" After ind prop");
					}
					
					if(burn) ind_obs_prob_update(isp);
				}
			}
		}
	}
}


/// Performs a resimulation proposal update for an unobserved individual 
void Proposal::resimulate_ind_unobs(State &state)
{
	auto pl = false;
	
	auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	
	if(sp.nindividual_obs == ssp.individual.size()) return;

	IndEvSampler ind_ev_samp(ssp.markov_eqn_vari,ssp.individual,model.details,sp,ssp.obs_eqn_value,ssp.obs_trans_eqn_value,model.eqn,state.genetic_value.inf_node,state.param_val,state.spline_val,state.popnum_t);
	
	const auto &so_samp = ssp.source_sampler;
	
	auto copy_init = false;
	if(sp.init_cond.type != INIT_POP_DIST && sp.contains_source == false) copy_init = true;
	
	vector <TrigEventRef> ste;
	
	for(auto i = sp.nindividual_obs; i < ssp.individual.size(); i++){	
		auto &ind = ssp.individual[i];
		
		if(pl){
			cout << endl << endl;
			cout << "older:"; ssp.print_event(ind.ev);
		}
		
		auto probif = 0.0, probfi = 0.0;
	
		Event e_init;
		if(copy_init == true) e_init = ind.ev[0]; 
		else e_init = so_samp.sample(probif);
		auto ev_new = ind_ev_samp.simulate_events(i,e_init,probif,ste);
		
		ntr++;
		if(ind_ev_samp.illegal == true) nfa++;
		else{
			if(testing == true){                         // Diagnostic
				auto prob = ind_ev_samp.simulate_events_prob(i,ev_new,ste);
				if(copy_init == false) prob += so_samp.sample_prob(e_init);
				if(dif(probif,prob,DIF_THRESH)){
					emsg("sampler different3");
				}
			}
			
			auto ev_store = ind.ev;			
			
			if(pl){
				cout << endl << endl << "Start";
				cout << "begin:"; ssp.print_event(ev_store);
				cout << "end:"; ssp.print_event(ev_new);
			}
			
			auto gc = state.update_tree(p_prop,i,ev_new);
				
			if(gc.type == GENCHA_FAIL){
				nfa++;
			}
			else{
				probfi += ssp.nm_obs_dprob(ind);
				auto like_ch = state.update_ind(p_prop,i,ev_new,UP_SINGLE);
				probif += ssp.nm_obs_dprob(ind);
			
				probfi += ind_ev_samp.simulate_events_prob(i,ev_store,ste);
				if(copy_init == false) probfi += so_samp.sample_prob(ev_store[0]);
			
				auto dprob = probfi-probif;
				gc.update_like_ch(like_ch,dprob);
				
				auto al = calculate_al(like_ch,dprob);
			
				if(pl) cout << al << " try ind" << endl;
				
				if(ran() < al){ 
					if(pl) cout << "ac ind" << endl;
					nac++;
					state.add_like(like_ch);
					state.gen_change_update(gc); 	
					if(sp.trans_tree) state.update_popnum_ind(p_prop,i);
				}
				else{
					state.restore_back();
				}	
			}
			if(pl) state.check("unobs ev");
		}
	}
}


/// Add and removes indiviual from the system
void Proposal::add_rem_ind(State &state)
{
	auto pl = false;
	
	const auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	const auto &so_samp = ssp.source_sampler;
	auto &inf_node = state.genetic_value.inf_node;
	
	IndEvSampler ind_ev_samp(ssp.markov_eqn_vari,ssp.individual,model.details,sp,ssp.obs_eqn_value,ssp.obs_trans_eqn_value,model.eqn,state.genetic_value.inf_node,state.param_val,state.spline_val,state.popnum_t);
	
	vector <Event> ev_empty;
	
	auto nobs = sp.nindividual_obs;
	
	// We need to keep track of the total number of source transition and observed number of source transitions
	
	auto Nso_obs = ssp.source_num(0,nobs);
	auto Nso_unobs = ssp.source_num(nobs,ssp.individual.size());
	
	ListMap unobs(nobs,ssp.individual.size());
	if(testing) unobs.check(ssp.individual.size()-nobs);

	vector <TrigEventRef> ste;

	for(auto loop = 0u; loop < 6; loop++){
		auto nunobs = unobs.list.size();
		state.back_init();
		
		string warn;
		int dN = normal_int_sample(si,warn);
		if(dN == UNSET) emsg("MCMC proposal problem: "+warn);
			
		auto Nso_unobs_new = Nso_unobs;
	
		if(dN != 0){
			ntr++;
			
			if(pl) cout << endl << endl << "Start" << endl;
				
			if(dN > 0){                                  // Adds individuals
				vector <unsigned int> i_store;
				
				Like like_ch;
				
				auto probif = 0.0;
				for(auto k = 0; k < dN; k++){
					auto i = ssp.add_individual(UNOBSERVED_IND);
					
					auto probif_st = probif;
					
					auto e_init = so_samp.sample(probif);
		
					auto ev_new = ind_ev_samp.simulate_events(i,e_init,probif,ste);

					if(pl){
						cout << endl << "New event add individual" << endl;					
						ssp.print_event(ev_new);
					}
					
					if(ind_ev_samp.illegal) emsg("SHould not be illegal");
					
					if(testing == true){                     // Diagnostic
						auto prob = so_samp.sample_prob(e_init) +
            						ind_ev_samp.simulate_events_prob(i,ev_new,ste);
						if(dif(probif-probif_st,prob,DIF_THRESH)){
							emsg("sampler different4a");
						}
					}
			
					i_store.push_back(i);
					
					auto li_ch = state.update_ind(p_prop,i,ev_new,UP_MULTI);
					probif += ssp.nm_obs_dprob(ssp.individual[i]);
					
					if(ssp.source_ind(i)) Nso_unobs_new++;
					
					add_on_like(li_ch,like_ch);
				}
				
				auto probfi = -factorial(nunobs+dN) + factorial(nunobs);
					
				// ADDITION
				// This assumes likelihood has a factor Nso_total!/Nso_unobs! implicitly included 
				// This was found to be neccesary, e.g. for a simple model with a source and some observed individuals
				if(Nso_unobs_new != Nso_unobs){
					probfi -= factorial(Nso_obs + Nso_unobs_new) - factorial(Nso_obs + Nso_unobs);
					probfi += factorial(Nso_unobs_new) - factorial(Nso_unobs); 	
				}
				
				auto al = calculate_al(like_ch,probfi-probif);
				
				if(pl) cout << al << "al add" << endl;
				if(ran() < al){ 
					if(pl) cout << "ac ind" << endl;
					nac++;
					state.add_like(like_ch);
		
					for(auto i : i_store){
						unobs.add(i);
					}
					
					Nso_unobs = Nso_unobs_new;
					
					update_si(0.005); if(si > ADD_REM_IND_MAX) si = ADD_REM_IND_MAX;
				}
				else{
					state.restore_back();
					for(int k = dN-1; k >= 0; k--){
						auto i = i_store[k];	
						ssp.remove_individual(i,inf_node);
					}
					update_si(-0.005); if(si < 1) si = 1;
				}
			}
			else{                                        // Removes individuals
				dN *= -1;
			
				if(int(nunobs) < dN){
					update_si(-0.005); if(si < 1) si = 1;
				}
				else{
					vector < vector <Event> > event_store;
					vector <unsigned int> i_store;

					Like like_ch;
				
					auto probif = -factorial(nunobs) + factorial(nunobs-dN);
				
					auto probfi = 0.0;
					for(auto k = 0; k < dN; k++){
						auto m = (unsigned int)(ran()*nunobs);
						auto i = unobs.list[m];
						unobs.remove(i);
						nunobs--;
						
						if(ssp.source_ind(i)) Nso_unobs_new--;
						
						auto ev_store = ssp.individual[i].ev;
						i_store.push_back(i);
						
						if(pl){
							cout << endl << "Remove event individual" << endl;		
							ssp.print_event(ev_store);
						}
					
						const auto &e_init = ev_store[0];
						if(so_samp.enter_frac == 0 && e_init.type == ENTER_EV) emsg("prob samp");
						if(so_samp.enter_frac == 1 && e_init.type != ENTER_EV) emsg("prob samp2");
							
						probfi += ssp.nm_obs_dprob(ssp.individual[i]);
						auto li_ch = state.update_ind(p_prop,i,ev_empty,UP_MULTI);
						add_on_like(li_ch,like_ch);
					
						probfi += so_samp.sample_prob(ev_store[0]) + 
						       ind_ev_samp.simulate_events_prob(i,ev_store,ste);
						
						event_store.push_back(ev_store);			
					}		
					
					// ADDITION
					// This assumes likelihood has a factor Nso_total!/Nso_unobs! implicitly included 
					// This was found to be neccesary, e.g. for a simple model with a source and some observed individuals
					if(Nso_obs != 0 && Nso_unobs_new != Nso_unobs){
						probif += factorial(Nso_obs + Nso_unobs_new) - factorial(Nso_obs + Nso_unobs);
						probif -= factorial(Nso_unobs_new) - factorial(Nso_unobs); 
					}
				
					auto al = calculate_al(like_ch,probfi-probif);
					
					if(pl) cout << al << "al remove" << endl;
					if(ran() < al){ 
						if(pl) cout << "ac ind" << endl;
						nac++;
						state.add_like(like_ch);
						
						sort(i_store.begin(),i_store.end());   
						for(int k = dN-1; k >= 0; k--){
							auto i = i_store[k];
							unobs.remove_hole(i);
							ssp.remove_individual(i,inf_node);
						}
						
						Nso_unobs = Nso_unobs_new;
					
						update_si(0.005); if(si > ADD_REM_IND_MAX) si = ADD_REM_IND_MAX;
					}
					else{
						for(int k = dN-1; k >= 0; k--){
							auto i = i_store[k];
							unobs.add_hole(i);
						}
						state.restore_back();
						update_si(-0.005); if(si < 1) si = 1;
					}
				}				
			}
		}
		if(pl) state.check(" Add rem prop");
	}
	
	if(Nso_unobs != ssp.source_num(nobs,ssp.individual.size())) emsg("Source number problem");
	
	if(testing) unobs.check(ssp.individual.size()-nobs);
}


/// Add and removes single indiviual from the system (accounts for transmission tree)
void Proposal::add_rem_tt_ind(State &state)
{
	auto pl = false;
	
	const auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	const auto &so_samp = ssp.source_sampler;
	
	IndEvSampler ind_ev_samp(ssp.markov_eqn_vari,ssp.individual,model.details,sp,ssp.obs_eqn_value,ssp.obs_trans_eqn_value,model.eqn,state.genetic_value.inf_node,state.param_val,state.spline_val,state.popnum_t);
	
	auto &inf_node = state.genetic_value.inf_node;
	
	vector <Event> ev_empty;
	vector <TrigEventRef> ste;
	
	auto nobs = sp.nindividual_obs;
	
	auto Nso_obs = ssp.source_num(0,nobs);
	auto Nso_unobs = ssp.source_num(nobs,ssp.individual.size());
	
	for(auto loop = 0u; loop < 6; loop++){
		ntr++;
	
		auto Nso_unobs_new = Nso_unobs;
	
		if(ran() < 0.5){                               // Add individual	
			auto probif = 0.0;
			auto i = ssp.add_individual(UNOBSERVED_IND);
					
			auto e_init = so_samp.sample(probif);
			auto ev_new = ind_ev_samp.simulate_events(i,e_init,probif,ste);	
			
			if(ind_ev_samp.illegal) emsg("SHould not be illegal");
		
			if(testing == true){                         // Diagnostic
				auto prob = so_samp.sample_prob(e_init) +
    				ind_ev_samp.simulate_events_prob(i,ev_new,ste);
				if(dif(probif,prob,DIF_THRESH)) emsg("sampler different4b");
			}
		
			auto gc = state.update_tree(p_prop,i,ev_new);
					
			if(gc.type == GENCHA_FAIL) nfa++;
			else{	
				auto li_ch = state.update_ind(p_prop,i,ev_new,UP_SINGLE);
				probif += ssp.nm_obs_dprob(ssp.individual[i]);
				if(ssp.source_ind(i)) Nso_unobs_new++;
			
				auto nunobs = ssp.individual.size()-nobs;
				auto probfi = -log(nunobs);
					
				// ADDITION
				// This assumes likelihood has a factor Nso_total!/Nso_unobs! implicitly included 
				// This was found to be neccesary, e.g. for a simple model with a source and some observed individuals
				if(Nso_obs != 0 && Nso_unobs_new != Nso_unobs){
					probfi -= factorial(Nso_obs + Nso_unobs_new) - factorial(Nso_obs + Nso_unobs);
					probfi += factorial(Nso_unobs_new) - factorial(Nso_unobs); 	
				}
				
				auto dprob = probfi-probif;
				gc.update_like_ch(li_ch,dprob);
	
				auto al = calculate_al(li_ch,dprob);
				
				if(pl) cout << al << " al" << endl; 
				if(ran() < al){ 
					if(pl) cout << "ac ind" << endl;
					nac++;
					state.add_like(li_ch);
					state.gen_change_update(gc); 	
					if(sp.trans_tree) state.update_popnum_ind(p_prop,i);
					
					Nso_unobs = Nso_unobs_new;
					update_si(0.005);
				}
				else{
					state.restore_back();
					ssp.remove_individual(i,inf_node);
					update_si(-0.005); if(si < 1) si = 1;
				}
			}
		}
		else{                                          // Removes individual
			auto nunobs = ssp.individual.size()-sp.nindividual_obs;
			if(nunobs == 0) nfa++;
			else{
				auto i = sp.nindividual_obs + (unsigned int)(ran()*nunobs);
				if(ssp.source_ind(i)) Nso_unobs_new--;
						
				auto fl = false;                           // Checks to see for onward transmission
				for(const auto &ev : ssp.individual[i].ev){
					if(ev.inf_node_ref != UNSET){
						if(inf_node[ev.inf_node_ref].inf_ev.size() > 0){ fl = true; break;}
					}
				}
				
				if(fl == false){
					auto probif = log(1.0/nunobs);
					
					// ADDITION
					// This assumes likelihood has a factor Nso_total!/Nso_unobs! implicitly included 
					// This was found to be neccesary, e.g. for a simple model with a source and some observed individuals
					if(Nso_unobs_new != Nso_unobs){
						probif += factorial(Nso_obs + Nso_unobs_new) - factorial(Nso_obs + Nso_unobs);
						probif -= factorial(Nso_unobs_new) - factorial(Nso_unobs); 
					}
					
					auto probfi = 0.0;
						
					auto ev_store = ssp.individual[i].ev;
						
					auto gc = state.update_tree(p_prop,i,ev_empty);
						
					if(gc.type == GENCHA_FAIL) nfa++;
					else{	
						const auto &e_init = ev_store[0];
						if(so_samp.enter_frac == 0 && e_init.type == ENTER_EV) emsg("prob samp");
						if(so_samp.enter_frac == 1 && e_init.type != ENTER_EV) emsg("prob samp2");
							
						probfi += ssp.nm_obs_dprob(ssp.individual[i]);
						auto li_ch = state.update_ind(p_prop,i,ev_empty,UP_SINGLE);
						
						probfi += so_samp.sample_prob(ev_store[0]) + 
						       ind_ev_samp.simulate_events_prob(i,ev_store,ste);
							
						auto dprob = probfi-probif;
						gc.update_like_ch(li_ch,dprob);
					
						auto al = calculate_al(li_ch,dprob);
						
						if(pl) cout << al << " al" << endl;
						if(ran() < al){ 
							if(pl) cout << "ac ind" << endl;
							nac++;
							state.add_like(li_ch);
							state.gen_change_update(gc); 	
							if(sp.trans_tree) state.update_popnum_ind(p_prop,i);
												
							ssp.remove_individual(i,inf_node);
							Nso_unobs = Nso_unobs_new;
							update_si(0.005);
						}
						else{
							state.restore_back();
							update_si(-0.005); if(si < 1) si = 1;
						}
					}				
				}
			}
		}
		if(pl) state.check(" Add rem prop");
	}
	
	if(Nso_unobs != ssp.source_num(nobs,ssp.individual.size())) emsg("Source number problem");
}


/// Performs a Metropolis-Hastings proposal update for individual effects
void Proposal::MH_ie(State &state)
{	
	auto pl = false;
	
	auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	
	const auto &inde = sp.ind_effect[ie_prop];
	
	auto in_me = false; 
	if(inde.markov_eqn_ref.size() > 0) in_me = true;
	
	auto in_nm = false; 
	if(inde.nm_trans_ref.size() > 0 || inde.nm_trans_incomp_ref.size() > 0) in_nm = true;
	
	auto in_pop = false;
	for(const auto &po : model.pop){
		for(auto ie_pop : po.ind_eff_mult){
			if(ie_pop == ie_prop) in_pop = true;
		}
	}
	
	auto in_transtree = false;
	if(sp.trans_tree){
		for(auto po : inde.pop_ref){
			for(const auto &pme : model.pop[po].markov_eqn_ref){
				if(model.species[pme.p].markov_eqn[pme.e].infection_trans) in_transtree = true;
			}
		}
	}

	auto g = inde.index;
	auto num = inde.num;
	
	auto &iegs = ssp.ind_eff_group_sampler[g];

	auto var = iegs.omega[num][num];
	auto sd = sqrt(var);
	
	vector <double> store;
	vector <NMupdate> update_nm;
	
	for(auto i = 0u; i < ssp.individual.size(); i++){
		auto &ind = ssp.individual[i];
		if(pl) cout << "Individual: " << ind.name << endl;
		
		auto value_st = ind.ie[ie_prop];
		
		string warn;
		auto d = normal_sample(0,sd*si,warn);
		if(d == UNSET) emsg("Problem with MCMC proposal: "+warn);
		
		auto value_prop = value_st + d;
		
		auto exp_old = ind.exp_ie[ie_prop];
		auto exp_prop = exp(value_prop-0.5*var);
		
		auto factor = exp_prop/exp_old;
		
		Like like_ch;
		
		if(in_me == true){                       // If in a markov_eqn then updates likelihood
			store = ssp.likelihood_ie_change(i,ie_prop,factor,like_ch);	
		}
		
		if(in_pop == true){                // If in a population then updates likelihoods from that
			state.update_ie_population(p_prop,i,ie_prop,exp_prop-exp_old,like_ch.markov);
		}
		
		if(in_nm == true){                 // If in a nm_trans then updates likelihood
			update_nm = ssp.likelihood_ie_nm_trans_change(i,ie_prop,factor,state.popnum_t,like_ch.nm_trans);	
		}
		
		if(in_transtree == true){
			state.update_ie_trans_tree(p_prop,i,ie_prop,exp_prop/exp_old,like_ch.markov);
		}
	
		auto dLi_ie = ssp.likelihood_indeff_change(i,ie_prop,value_prop,value_st);
	
		like_ch.ie += dLi_ie;
		auto al = calculate_al(like_ch,0);
	
		if(pl) cout << al << " al" << endl;
		ntr++;
		if(ran() < al){ 
			if(pl) cout << "ac ind" << endl;
			nac++;
			state.add_like(like_ch);
			
			ind.ie[ie_prop] = value_prop;
			ind.exp_ie[ie_prop] = exp_prop;
			ssp.Li_ie[g] += dLi_ie;
			
			ssp.likelihood_ie_nm_trans_update(update_nm);
		
			if(inde.pop_ref.size() > 0) state.popnum_ind_recalc_w(p_prop,i);
			
			update_si(0.005);
		}
		else{ 
			if(pl) cout << "rej ind" << endl;
			
			if(in_me == true){
				ssp.likelihood_ie_change_restore(i,ie_prop,factor,store);
			}
			
			if(in_pop == true){             // If in a population then updates likelihoods from that
				state.update_ie_population(p_prop,i,ie_prop,exp_old-exp_prop,like_ch.markov);
			}
			
			if(in_transtree == true){
				state.update_ie_trans_tree(p_prop,i,ie_prop,exp_old/exp_prop,like_ch.markov);
			}

			update_si(-0.01);
		}
		if(pl) state.check("ind prop");
	}
}


/// Performs a Metropolis-Hastings joint proposal update for individual effects and a variance
void Proposal::MH_ie_var(State &state)
{
	auto pl = false;
	
	auto param_store = state.param_val;
	
	state.param_val = sample(state.param_val);
	
	auto th = param_list[0];
	
	if(state.param_val[th] <= 0 || model.inbounds(state.param_val) == false || 
		model.ie_cholesky_error(state.param_val)){ 
		update_si(-0.005);
		state.param_val = param_store; 
		return;
	}
	
	auto factor = sqrt(state.param_val[th]/param_store[th]);

	auto &ssp = state.species[p_prop];
	for(auto &ind : ssp.individual) ind.ie[ie_prop] *= factor;
	
	state.update_spline(affect_spline);
	auto like_ch = state.update_param(affect_like,param_store);

	auto al = calculate_al(like_ch,-like_ch.ie);
	
	if(pl) cout << al << " al" << endl;
	ntr++;
	if(ran() < al){
		if(pl) cout << "accept" << endl;
		nac++;
		state.add_like(like_ch);
		update_si(0.01);
		state.remove_store_spline(affect_spline);
	}
	else{ 
		state.param_val = param_store;
		state.restore(affect_like);
		state.restore_spline(affect_spline);
		for(auto &ind : ssp.individual) ind.ie[ie_prop] /= factor;
		update_si(-0.005);
	}
	
	if(pl) state.check("covar ie");
}


/// Performs a Metropolis-Hastings joint proposal update for individual effects and a variance
void Proposal::MH_ie_var_cv(State &state)
{
	auto pl = false;//true;
	
	auto param_store = state.param_val;
	
	state.param_val = sample(state.param_val);
	
	if(model.inbounds(state.param_val) == false || 
		 model.ie_cholesky_error(state.param_val)){ 
		update_si(-0.005);
		state.param_val = param_store; 
		return;
	}
	
	auto th = param_list[0];
	
	if(state.param_val[th] < TINY){ state.param_val = param_store; return;}
	
	auto factor = sqrt(state.param_val[th]/param_store[th]);

	auto &ssp = state.species[p_prop];
	for(auto &ind : ssp.individual) ind.ie[ie_prop] *= factor;
	
	state.update_spline(affect_spline);
	auto like_ch = state.update_param(affect_like,param_store);

	auto al = calculate_al(like_ch,-like_ch.ie);
	
	if(pl) cout << al << " al" << endl;
	ntr++;
	if(ran() < al){
		if(pl) cout << "accept" << endl;
		nac++;
		state.add_like(like_ch);
		update_si(0.01);
		state.remove_store_spline(affect_spline);
	}
	else{ 
		state.param_val = param_store;
		state.restore(affect_like);
		state.restore_spline(affect_spline);
		for(auto &ind : ssp.individual) ind.ie[ie_prop] /= factor;
		update_si(-0.005);
	}
	
	if(pl) state.check("covar ie");
}


/// Performs a Metropolis-Hastings joint proposal update for individual effects and a variance
void Proposal::MH_ie_covar(State &state)
{
	auto pl = false;
	
	auto param_store = state.param_val;

	state.param_val = sample(state.param_val);
	
	if(model.inbounds(state.param_val) == false){ 
		state.param_val = param_store; 
		update_si(-0.005);
		return;
	}
	
	auto th = param_list[0];
	auto g = ind_eff_group_ref.ieg;
	
	auto &ssp = state.species[p_prop];
	const auto &iegs = ssp.ind_eff_group_sampler[g];
	
	const auto &sp = model.species[p_prop];
	const auto &ieg = sp.ind_eff_group[ind_eff_group_ref.ieg];

	auto E = ieg.list.size();
	
	auto fac = state.param_val[th]/param_store[th]; 
		
	auto omega_bef = iegs.omega;
	auto omega_aft = omega_bef;

	{		
		auto i = ind_eff_group_ref.i;
		auto j = ind_eff_group_ref.j;
		omega_aft[i][j] *= fac;
		omega_aft[j][i] *= fac;
	}
	
	auto illegal_aft = false;
	auto Z_omega_aft = calculate_cholesky(omega_aft,illegal_aft);
	if(illegal_aft){
		state.param_val = param_store; 
		update_si(-0.005);
		return;
	}
	
	auto illegal_bef = false;
	auto Z_omega_bef = calculate_cholesky(omega_bef,illegal_bef);
	if(illegal_bef) emsg("SHould not have cholesky error");
		
	auto inv_Z_omega_bef = invert_matrix(Z_omega_bef);
			
	auto M = matrix_mult(Z_omega_aft,inv_Z_omega_bef);
	
	vector <unsigned int> index;
	for(auto e = 0u; e < E; e++){
		index.push_back(ieg.list[e].index);
	}	
	
	//ssp.compare_covar("BEFORE",omega_bef,g);
	
	vector <double> vec(E), vec_new(E);
	vector < vector <double> > store;
	for(auto i = 0u; i < ssp.individual.size(); i++){
		auto &ie = ssp.individual[i].ie;
		
		for(auto e = 0u; e < E; e++) vec[e] = ie[index[e]];
		store.push_back(vec);
		
		vec_new = matrix_mult(M,vec);
		for(auto e = 0u; e < E; e++) ie[index[e]] = vec_new[e];
	}

	//ssp.compare_covar("AFTER",omega_aft,g);

	state.update_spline(affect_spline);
	auto like_ch = state.update_param(affect_like,param_store);

	auto al = calculate_al(like_ch,-like_ch.ie);
	
	if(pl) cout << al << " al" << endl;
	ntr++;
	if(ran() < al){ 
		if(pl) cout << "ac" << endl;
		nac++;
		state.add_like(like_ch);
		update_si(0.01);
		state.remove_store_spline(affect_spline);
	}
	else{ 
		state.param_val = param_store;
		state.restore(affect_like);
		state.restore_spline(affect_spline);
		
		for(auto i = 0u; i < ssp.individual.size(); i++){
			auto &ie = ssp.individual[i].ie;
			for(auto e = 0u; e < E; e++) ie[index[e]] = store[i][e];
		}
		update_si(-0.005);
	}
	if(pl) state.check("uuu");
}


/// Proposals which jointly changes parameter and events
void Proposal::param_event_joint(Direction dir, State &state)
{
	auto pl = false;
	
	auto th = param_list[0];
	
	auto param_new = sample(state.param_val);
	
	ntr++; 	
	if(model.inbounds(param_new) == false){ 
		nfa++;
		update_si(-0.005);
		return;
	}

	auto illegal = false;

	auto fac = param_new[th]/state.param_val[th];
	if(fac < 0) illegal = true;

	const auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];

	// Finds new events 	
	auto tmin = model.details.t_start;
	auto tmax = model.details.t_end;
	
	vector < vector <Event> > ev_new, ev_old;
	vector < vector <IncompNMTransRef> > incomp_ref_old, incomp_ref_new;
	
	vector <unsigned int> ind_list;
	auto factor = 1.0;

	auto incomp_fl = false; if(sp.nm_flag && dir == FORWARD) incomp_fl = true;
		
	for(auto i = 0u; i < ssp.individual.size(); i++){
		const auto &ind = ssp.individual[i];
	
		auto event = ind.ev;
		
		auto cl_sel = UNSET;
		
		switch(dir){
		case FORWARD:  // Scales selected trans and shift all events after
			{
				double sh = UNSET;
				for(auto e = 1u; e < event.size(); e++){
					auto &ev = event[e];
						
					if((ev.type == M_TRANS_EV || ev.type == NM_TRANS_EV) && tr_change[ev.tr_gl]){
						if(sh != UNSET || ev.observed) illegal = true; 
						else{
							auto e_origin = ev.e_origin;
							if(ev.type == M_TRANS_EV) e_origin = get_event_before(e,event);
							
							sh = (ev.t-ind.ev[e_origin].t)*(fac-1);
							
							ev.t += sh;
							if(ev.t > tmax) illegal = true;
							cl_sel = ev.cl;
						}
					}
					else{
						if(ev.cl == cl_sel){
							if(ev.observed) illegal = true;
							ev.t += sh;
						}
					}
				}
			}
			break;
	
		case FORWARD_SQ:
			{
				auto tend = tmax;
				
				double fac_squash = UNSET; 
				
				for(auto e = 1u; e < event.size(); e++){
					auto &ev = event[e];
						
					if((ev.type == M_TRANS_EV || ev.type == NM_TRANS_EV) && tr_change[ev.tr_gl]){
						if(fac_squash != UNSET || ev.observed) illegal = true; 
						else{
							auto e_origin = ev.e_origin;
							if(ev.type == M_TRANS_EV) e_origin = get_event_before(e,event);
							
							auto sh = (ev.t-ind.ev[e_origin].t)*(fac-1);
							fac_squash = (tend-(ev.t+sh))/(tend-ev.t);
							
							if(fac_squash <= 0) illegal = true;
				
							ev.t += sh;
							cl_sel = ev.cl;
						}
					}
					else{
						if(ev.cl == cl_sel){
							if(ev.observed) illegal = true;
							ev.t = tend-(tend-ev.t)*fac_squash;
							factor *= fac_squash;
						}
					}
				}
			}
			break;
			
		case BACKWARD_SQ:
			{	
				double fac_squash = UNSET; 
				
				for(int e = event.size()-1; e >= 0; e--){
					auto &ev = event[e];
						
					if((ev.type == M_TRANS_EV || ev.type == NM_TRANS_EV) && tr_change[ev.tr_gl]){
						if(fac_squash != UNSET) illegal = true; 
						else{
							auto e_origin = ev.e_origin;
							if(ev.type == M_TRANS_EV) e_origin = get_event_before(e,event);
							
							auto t_ori = ind.ev[e_origin].t;
							auto sh = (ev.t-t_ori)*(fac-1);
							fac_squash = (t_ori-tmin-sh)/(t_ori-tmin);
							
							if(fac_squash <= 0) illegal = true; 
							
							cl_sel = ev.cl;
						}
					}
					else{
						if(ev.cl == cl_sel){
							ev.t = tmin+(ev.t-tmin)*fac_squash;
							factor *= fac_squash;
						}
					}
				}
			}
			break;
	
		default: emsg("SHould not be here");
		}
		
		if(illegal) break;
		
		if(cl_sel != UNSET){
			ind_list.push_back(i);
			ev_new.push_back(event);
		}
	}
	
	if(illegal){ nfa++; update_si(-0.005); return;}

	ssp.remove_all_event_ref(ind_list);
	if(incomp_fl) ssp.remove_all_nmincomp_ref(ind_list);

	for(auto j = 0u; j < ind_list.size(); j++){
		ssp.make_consistent(ev_new[j]);
		
		const auto &ind = ssp.individual[ind_list[j]];
		ev_old.push_back(ind.ev);
		if(incomp_fl){
			auto inrefnew = ind.incomp_ref;
			incomp_ref_old.push_back(inrefnew);
			ssp.adjust_incomp_ref(inrefnew,ev_new[j]);
			incomp_ref_new.push_back(inrefnew);
		}
	}
	
	ssp.add_all_event_ref(ind_list,ev_new);
	if(incomp_fl) ssp.add_all_nmincomp_ref(ind_list,incomp_ref_new);

	auto param_store = state.param_val;
	state.param_val = param_new;
	
	// Update param
	state.update_spline(affect_spline);
	auto like_ch = state.update_param(affect_like,param_store);
	
	double n = ind_list.size();

	auto al = calculate_al(like_ch,n*log(fac)+log(factor));
	
	if(pl) cout << al << " al" << endl;
		
	if(ran() < al){ 
		if(pl) cout << "accept" << endl;
		nac++;
			
		state.add_like(like_ch);
		update_si(0.01);
		state.remove_store_spline(affect_spline);
	}
	else{ 
		if(pl) cout << "reject" << endl;
		
		state.param_val = param_store;
		state.restore(affect_like);
		state.restore_spline(affect_spline);
		
		ssp.remove_all_event_ref(ind_list);
		ssp.add_all_event_ref(ind_list,ev_old);
	
		if(incomp_fl){
			ssp.remove_all_nmincomp_ref(ind_list);
			ssp.add_all_nmincomp_ref(ind_list,incomp_ref_old);
		}
		
		update_si(-0.005);
	}

	if(pl) state.check("ev");
}


/// Proposal which changes init_cond fractions
void Proposal::init_cond_frac(State &state)
{
	const auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	const auto &ic = sp.init_cond;
	auto &icv = ssp.init_cond_val;
	auto C = sp.comp_gl.size();
	
	auto foc_cl = ic.focal_cl;
	
	if(foc_cl == UNSET){
		auto alp = ic.alpha;
		for(auto c = 0u; c < C; c++){
			if(alp[c] != ALPHA_ZERO) alp[c] += icv.cnum[c];
		}
		
		icv.frac = dirichlet_sample(alp);
	}
	else{
		auto alp_focal = ic.alpha_focal;
		
		for(auto cgl = 0u; cgl < C; cgl++){
			const auto &co = sp.comp_gl[cgl];
			auto n = icv.cnum[cgl];
			if(n != 0){
				for(auto cl = 0u; cl < sp.ncla; cl++){
					if(cl != foc_cl){
						auto c = co.cla_comp[cl];
						if(alp_focal[cl][c] != ALPHA_ZERO) alp_focal[cl][c] += n;
					}
				}
			}
		}
			
		for(auto cl = 0u; cl < sp.ncla; cl++){
			if(cl != foc_cl){
				icv.frac_focal[cl] = dirichlet_sample(alp_focal[cl]);
			}
		}
		
		icv.frac_comb = model.frac_combine(p_prop,icv.frac_focal);
	}
	
	ssp.likelihood_init_cond(state.like.init_cond);
	state.prior_init_cond(state.like.init_cond_prior);
}


/// Calculates acceptance probability base don changes in likelihoods
double Proposal::calculate_al(const Like &like_ch, double dprob) const
{						
	return calc_al(like_ch,dprob,burn_info);
}


/// Provides diagnostic information about the proposals
string Proposal::diagnostics(long total_time) const
{
	stringstream ss;
	
	if(!on) ss << "<Switched off>" << endl;
	switch(type){
	case TRANS_TREE_PROP:
		{
			ss << "Transmission tree proposal: " << endl;
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "% ";
			ss << endl;  
		}
		break;
		
	case TRANS_TREE_SWAP_INF_PROP:
		{
			ss << "Transmission tree swap inf proposal: " << endl;
			if(ntr == 0) ss << "No proposals";
			else{
				ss << " Acceptance: " << (unsigned int) (100.0*nac/(ntr-nfa+TINY)) << "% ";
				ss << " Fail: " << (unsigned int) (100.0*nfa/ntr) << "% ";
			}
			ss << endl;  
		}
		break;
		
	case TRANS_TREE_MUT_PROP:
		{
			ss << "Transmission tree mutation proposal: " << endl;
			if(ntr == 0) ss << "No proposals";
			else	ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "%  Size:" << si;	
			ss << endl;  
		}
		break;
		
	case TRANS_TREE_MUT_LOCAL_PROP:
		{
			ss << "Transmission tree local mutation proposal: " << endl;
			for(auto loop = 0u; loop < GEN_UPDATE_MAX; loop++){
				auto type = GenBranchUpdate(loop);
				switch(type){
				case OBS_GEN_UPDATE: ss << "Obs"; break;
				case SINGLE_BRANCH: ss << "Single branch"; break;
				case TWO_BRANCH1: ss << "Two branch1"; break;
				case TWO_BRANCH2: ss << "Two branch2"; break;
				case TWO_BRANCH3: ss << "Two branch3"; break;
				case ORIGIN_GEN_UPDATE: ss << "Origin"; break;
				default: emsg("not opt"); break;
				}
				ss << " - ";
				const auto &gmi = gen_mut_info[type];
				if(gmi.ntr == 0) ss << "No proposals";
				else ss << " Acceptance: " << (unsigned int) (100.0*gmi.nac/gmi.ntr) << "%  Size:" << gmi.si;	
			}
			ss << endl;  
		}
		break;
		
		
	case PARAM_PROP:
		{
			ss << "Parameter proposal: ";
			for(auto i = 0u; i < param_list.size(); i++){
				auto th = param_list[i];
				if(i != 0) ss << ", ";
				ss << model.param_vec[th].name; 
			}
			ss << endl;
			if(ntr == 0) ss << "No proposals";
			else{
				ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "%   Size:" << si;
			}
			ss << endl;
		}
		break;
	
	case PAR_EVENT_FORWARD_PROP: case PAR_EVENT_FORWARD_SQ_PROP:
	case PAR_EVENT_BACKWARD_SQ_PROP:
		{
			switch(type){
			case PAR_EVENT_FORWARD_PROP: ss << "Parameter event forward proposal: "; break;
			case PAR_EVENT_FORWARD_SQ_PROP: ss << "Parameter event forward squash proposal: "; break;
			case PAR_EVENT_BACKWARD_SQ_PROP: ss << "Parameter event backward squash proposal: "; break;
			default: emsg("not here"); break;
			}
				
			for(auto i = 0u; i < param_list.size(); i++){
				auto th = param_list[i];
				if(i != 0) ss << ", ";
				ss << model.param_vec[th].name; 
			}
			ss << endl;
			if(ntr == 0) ss << "No proposals";
			else{
				ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "%   Size:" << si;
			}
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
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "%   Size:" << si;
			ss << endl;
		}
		break;
		
	case MBPII_PROP:
		{
			ss << "Type II MBP proposal for species " << model.species[p_prop].name << endl;
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "%   Probability:" << si;
			ss << endl;
		}
		break;
		
	case MBP_IC_POP_PROP:
		{
			const auto &sp = model.species[p_prop];
			ss << "MBP proposal initial population for species "  << sp.name;
			ss << " for " << sp.cla[cl_prop].name << " compartment " << sp.cla[cl_prop].comp[c_prop].name << endl;
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "%   Size:" << si;
			ss << endl;
		}
		break;
		
	case MBP_IC_POPTOTAL_PROP:
		{
			const auto &sp = model.species[p_prop];
			ss << "MBP proposal total initial population for species "  << sp.name << endl;
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "%   Size:" << si;
			ss << endl;
		}
		break;
		
	case MBP_IC_RESAMP_PROP:
		{
			const auto &sp = model.species[p_prop];
			ss << "MBP proposal for resampling init_cond for species " << sp.name << endl;
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "%   Size:" << si;
			ss << endl;
		}
		break;
		
	case INIT_COND_FRAC_PROP:
		{
			const auto &sp = model.species[p_prop];
			ss << "Proposal for init_cond fractions: "  << sp.name << endl;
		}
		break;
	
	case IND_EVENT_TIME_PROP: case IND_MULTI_EVENT_PROP:
		{
			const auto &sp = model.species[p_prop];
			
			if(type == IND_EVENT_TIME_PROP) ss << "Individual event proposals ";
			else ss << "Individual multi-event proposals ";
			ss << endl;
			
			vector < vector < vector <double> > > ac_list,  fa_list, si_list;
			ac_list.resize(sp.ncla); fa_list.resize(sp.ncla); si_list.resize(sp.ncla);
			for(auto cl = 0u; cl < sp.ncla; cl++){
				const auto &claa = sp.cla[cl];
				ac_list[cl].resize(claa.ntra); fa_list[cl].resize(claa.ntra); si_list[cl].resize(claa.ntra);
			}
		
			auto ntr_tot = 0u;
			for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
				const auto &tra = sp.tra_gl[tr];
				const auto &samp = ind_sampler[tr];
				if(samp.ntr > 0){
					ac_list[tra.cl][tra.tr].push_back((unsigned int) (100.0*samp.nac/(samp.ntr+TINY)));
					fa_list[tra.cl][tra.tr].push_back((unsigned int) (100.0*samp.nfa/(samp.ntr+TINY)));
					si_list[tra.cl][tra.tr].push_back(samp.si);
				}
				
				ntr_tot += samp.ntr;
			}
		
			if(ntr_tot == 0) ss << "No proposals" << endl;
			else{
				for(auto cl = 0u; cl < sp.ncla; cl++){
					const auto &claa = sp.cla[cl];
					for(auto tr = 0u; tr < claa.ntra; tr++){
						const auto &tra = claa.tra[tr];
						
						if(ac_list[cl][tr].size() > 0){
							ss << replace_arrow(tra.name) << "   ";
							ss << " Acceptance: " << min(ac_list[cl][tr]) << " - " << max(ac_list[cl][tr]) << "%   ";
							ss << " Fail: " << min(fa_list[cl][tr]) << " - " << max(fa_list[cl][tr]) << "%   ";
							ss << "Size:" << min(si_list[cl][tr]) << " - " << max(si_list[cl][tr]) ;
							ss << endl;
						}
					}
				}
			}
		}
		break;
		
	case IND_EVENT_ALL_PROP: 
		{
			const auto &sp = model.species[p_prop];
			
			ss << "Individual all event proposals" << endl;
			
			vector <double> ac_list, fa_list, si_list;
		
			auto ntot = 0u;
			for(auto i = 0u; i < sp.nindividual_obs; i++){
				const auto &samp = ind_sampler[i];
				if(samp.ntr > 0){
					ac_list.push_back((unsigned int) (100.0*samp.nac/(samp.ntr+TINY)));
					fa_list.push_back((unsigned int) (100.0*samp.nfa/(samp.ntr+TINY)));
					si_list.push_back(samp.si);
				}
				ntot += samp.ntr;
			}
		
			if(ntot == 0) ss << "No proposals" << endl;
			else{
				ss << " Acceptance: " << int(mean(ac_list)) << " (" << min(ac_list) << " - " << max(ac_list) << ")    ";
				ss << " Fail: " << int(mean(fa_list)) << " (" << min(fa_list) << " - " << max(fa_list) << ")   ";
				ss << "Size:" << mean(si_list) << " (" << min(si_list) << " - " << max(si_list) << ")";
				ss << endl;
			}
		}
		break;
	
	case IND_OBS_RESIM_PROP:
		{
			ss << "Individual resimulation proposals (for observed)";
		
			ss << ":" << endl;	
					
			if(ntr == 0) ss << "No proposals";
			else{
				ss << " Acceptance: " << (unsigned int) (100.0*nac/(ntr+TINY)) << "%    ";
				ss << " Fail: " << (unsigned int) (100.0*nfa/(ntr+TINY)) << "%    ";
				
				vector <double> prob_list;
				for(const auto &isp : ind_sim_prob){
					if(isp.done) prob_list.push_back(isp.prob);
					else prob_list.push_back(0);
				}
				
				ss << " Prob: " << tstr(mean(prob_list),3) << " (" << tstr(min(prob_list),3) << " - " << tstr(max(prob_list),3) << ")  ";
			}
			ss << endl;
		}
		break;
		
	case IND_OBS_SAMP_PROP:
	case IND_OBS_RESIM_SINGLE_PROP:
		{
			ss << "Individual ";
			switch(type){
			case IND_OBS_SAMP_PROP: ss << "sampler proposals (for observed)"; break;
			case IND_OBS_RESIM_SINGLE_PROP: ss << "resimulation single cl proposals (for observed)"; break;
			default: emsg("def op"); break;
			}
		
			ss << ":" << endl;	
			if(ntr == 0) ss << "No proposals";
			else{
				ss << " Acceptance: " << (unsigned int) (100.0*nac/(ntr+TINY)) << "%    ";
				ss << " Fail: " << (unsigned int) (100.0*nfa/(ntr+TINY)) << "%    ";
				
				const auto &sp = model.species[p_prop];
				for(auto cl = 0u; cl < sp.ncla; cl++){
					vector <double> prob_list;
					for(const auto &isp : ind_sim_prob_cl){
						if(isp[cl].done) prob_list.push_back(isp[cl].prob);
						else prob_list.push_back(0);
					}			
					ss << " " << sp.cla[cl].name << " prob: " << tstr(mean(prob_list),3) << " (" << tstr(min(prob_list),3) << " - " << tstr(max(prob_list),3) << ")  ";
				}
			}
			ss << endl;
		}			
		break;
		
	case IND_UNOBS_RESIM_PROP:
		{
			ss << "Individual resimulation proposals (for unobserved)"; 
	
			ss << ":" << endl;	
			if(ntr == 0) ss << "No proposals";
			else{
				ss << " Acceptance: " << (unsigned int) (100.0*nac/(ntr+TINY)) << "%    ";
				ss << " Fail: " << (unsigned int) (100.0*nfa/(ntr+TINY)) << "%    ";
			}
			ss << endl;
		}
		break;
		
	case IND_OBS_SWITCH_ENTER_SOURCE_PROP:
		ss << "Switch enter source prop" << endl;
		if(ntr == 0) ss << "No proposals";
		else ss << " Acceptance: " << (unsigned int) (100.0*nac/(ntr+TINY)) << "%";
		ss << endl;
		break;

	case IND_OBS_SWITCH_LEAVE_SINK_PROP:
		ss << "Switch leave sink prop" << endl;
		if(ntr == 0) ss << "No proposals";
		else ss << " Acceptance: " << (unsigned int) (100.0*nac/(ntr+TINY)) << "%";
		ss << endl;
		break;
		
	case CORRECT_OBS_TRANS_PROP:
		ss << "Correct obs trans prop" << endl;
		break;
		
	case IND_LOCAL_PROP: 
		{
			auto &claa = model.species[p_prop].cla[cl_prop];
			ss << "Individual local proposals for " << claa.name << ":" << endl;
			for(auto j = 0u; j < claa.swap_rep.size(); j++){
				const auto &si = swap_info[j];
				ss << replace_arrow(claa.swap_rep[j].name) << " ";
				if(si.ntr == 0) ss << "No proposals";
				else{
					ss << " Acceptance: " << (unsigned int) (100.0*si.nac/(si.ntr2+TINY)) << "%  ";
					ss << " Fail: " << (unsigned int) (100.0*si.nfa/(si.ntr+TINY)) << "%  ";
					ss << " Zero: " << (unsigned int) (100.0*si.nzero/(si.ntr+TINY)) << "%  ";
					ss << " Filt: " << (unsigned int) (100.0*si.nfilt/(si.ntr+TINY)) << "%  ";
				}
				ss << endl;
			}
		}
		break;
		
	case IND_ADD_REM_PROP: 
		{
			ss << "Individual add/remove proposals (for unobserved):" << endl;
			
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/(ntr+TINY)) << "%  Size:" << si;;
			ss << endl;
		}
		break;
			
	case IND_ADD_REM_TT_PROP: 
		{
			ss << "Individual add/remove proposals with trans-tree (for unobserved):" << endl;
			
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/(ntr+TINY)) << "%";
			ss << endl;
		}
		break;
		
	case IE_PROP:
		{
			auto ie_name =  model.species[p_prop].ind_effect[ie_prop].name;
			ss << "Individual effect proposal for " << ie_name;
			ss << endl;
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "%   Size:" << si;
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
			
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "%   Size:" << si;
			ss << endl;
		}
		break;
		
	case IE_COVAR_PROP:
		{
			auto j = param_list[0];
			auto th = model.param_vec[j].th;
			auto name = model.param[th].name;
			
			ss << "Joint proposal on " << name << " and individual effects";
			ss << endl;
			
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "%   Size:" << si;
			ss << endl;
		}
		break;
		
	case IE_VAR_CV_PROP:
		{
			auto j = param_list[0];
			auto th = model.param_vec[j].th;
			auto name = model.param[th].name;
			auto ie_name = model.species[p_prop].ind_effect[ie_prop].name;

			ss << "Joint proposal on " << name << " and " << ie_name;
			
			ss << " chainging ";
			for(auto k = 1u; k < param_list.size(); k++){
				auto th = model.param_vec[param_list[k]].th;
				if(k != 1) ss << ",";
				ss << model.param[th].name;
			}
			ss << endl;
			
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/ntr) << "%   Size:" << si;
			ss << endl;
		}
		break;
		
	case POP_ADD_REM_LOCAL_PROP: 
		{
			ss << "Population add/remove local proposals: ";
			for(auto tr : loc_samp.tr_list) ss << model.species[p_prop].tra_gl[tr].name << ",";
			ss << endl;
			
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/(ntr+TINY)) << "%";
			ss << "  Window: " << loc_samp.win << endl;
		}
		break;
		
	case POP_MOVE_LOCAL_PROP: 
		{
			ss << "Population move local proposals:" << endl;
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/(ntr+TINY)) << "%";
			ss << endl;
		}
		break;
		
	case POP_IC_LOCAL_PROP: 
		{
			ss << "Population init cond local proposals:" << endl;
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/(ntr+TINY)) << "%";
			ss << "  Window: " << loc_samp.win << endl;
		}
		break;
		
	case POP_END_LOCAL_PROP: 
		{
			ss << "Population end local proposals:" << endl;
			
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/(ntr+TINY)) << "%";
			ss << "  Window: " << loc_samp.win << endl;
		}
		break;
		
	case POP_SINGLE_LOCAL_PROP: 
		{
			ss << "Population single add/remove local proposals:" << endl;
			
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/(ntr+TINY)) << "%";
			ss << endl;
		}
		break;
	
	case POP_IC_PROP: 
		{
			ss << "Population IC proposals:" << endl;
			
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/(ntr+TINY)) << "%";
			ss << endl;
		}
		break;
		
	case POP_IC_SWAP_PROP: 
		{
			ss << "Population IC swap proposals:" << endl;
			
			if(ntr == 0) ss << "No proposals";
			else ss << " Acceptance: " << (unsigned int) (100.0*nac/(ntr+TINY)) << "%";
			ss << endl;
		}
		break;
	}
	
	ss << "CPU time: " << cpu_percent(timer[PROP_TIMER],total_time) << endl;
	//	 << "  Prob: " << prop_prob << endl;
	
	return ss.str();
}


/// Switchs betweeen an enter and source event
void Proposal::switch_enter_source(State &state)
{
	auto pl = false;

	const auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];

	const auto &det = sp.details;
	auto t_start = det.t_start;
	auto t_end = det.t_end;
	const auto tlim = t_start+0.05*(t_end-t_start);
	
	for(auto i = 0u; i < sp.nindividual_obs; i++){
		const auto &indd = sp.individual[i];
		const auto &ind = ssp.individual[i];
		const auto &ev = ind.ev;
		const auto e_init = ev[0];
	
		if(indd.enter_ref == UNSET) emsg("enter_ref should be set");
		if(!sp.enter[indd.enter_ref].set){
			// Get the maximum time over which a swap could be done
			auto tmax = tlim;
			if(ev.size() > 1){
				if(ev[1].t < tmax) tmax = ev[1].t; 
			}
			
			if(indd.obs.size() > 0){
				const auto &ob = indd.obs[0];
				auto tma = ob.t;
				if(indd.obs[0].not_alive){ // If there are not alive observation account for
					auto trange = ssp.source_time_range(indd);
					tma = trange.tmax;
				}
				
				if(tma < tmax) tmax = tma;
			}
		
			Event e_new;
			
			double dprob = UNSET;
			switch(e_init.type){
			case ENTER_EV:  // Switches from enter to source
				{
					auto c = e_init.c_after;
					auto tr = sp.cgl_tr_source[c];
					if(tr != UNSET){
						const auto &tra = sp.tra_gl[tr];
						
						auto t = t_start + ran()*(tmax-t_start);
						
						dprob = log(tmax-t_start);

						IndInfFrom inf_from;
						e_new = ssp.get_event(M_TRANS_EV,i,tr,UNSET,UNSET,tra.f,t,inf_from);
						
						if(event_near_div(t,model.details)) dprob = UNSET;
					}
				}
				break;
				
			case M_TRANS_EV:  // Swiches from source to enter
				if(e_init.observed == false && e_init.t < tmax){
					const auto &tra = sp.tra_gl[e_init.tr_gl];
						
					dprob = -log(tmax-t_start);
					
					IndInfFrom inf_from;
					e_new = ssp.get_event(ENTER_EV,i,UNSET,UNSET,UNSET,tra.f,t_start,inf_from);
					e_new.observed = true;
				}
				break;
				
			default: emsg("Wrong start"); break;
			}
				
			if(dprob != UNSET){
				ntr++;
				
				auto ev_new = ev;
				ev_new[0] = e_new;
				
				if(pl){
					cout << endl << endl;
					cout << ind.name << endl;
					cout << "event old:" << endl; ssp.print_event(ind.ev);
					cout << "event new:" << endl; ssp.print_event(ev_new);
				}
			
				if(!events_near_div(ev_new,model.details)){
					auto gc = state.update_tree(p_prop,i,ev_new);
					if(gc.type == GENCHA_FAIL){
						nfa++;
					}
					else{
						auto like_ch = state.update_ind(p_prop,i,ev_new,UP_SINGLE);
					
						gc.update_like_ch(like_ch,dprob);
						
						auto al = calculate_al(like_ch,dprob);
				
						if(pl) print_like(like_ch);
				
						if(pl) cout << i << " " <<  al << " al" << endl;
					
						if(ran() < al){
							if(pl) cout << " accept" << endl;
							nac++;
							state.add_like(like_ch);
							state.gen_change_update(gc); 
							if(sp.trans_tree) state.update_popnum_ind(p_prop,i);
						}
						else{ 
							if(pl) cout << "reject" << endl;
							state.restore_back();
						}

						if(pl) state.check("ind prop");
					}
				}
			}
		}
	}
}


/// Switchs betweeen an leave and sink event
void Proposal::switch_leave_sink(State &state)
{
	auto pl = false;
	
	const auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];

	const auto &det = sp.details;
	auto t_start = det.t_start;
	auto t_end = det.t_end;
	const auto tlim = t_end-0.05*(t_end-t_start);
	
	for(auto i = 0u; i < sp.nindividual_obs; i++){
		const auto &indd = sp.individual[i];
		const auto &ind = ssp.individual[i];
		const auto &ev = ind.ev;
	
		// Get the maximum time over which a swap could be done
		auto tmin = tlim;
		if(indd.tmax > tmin) tmin = indd.tmax;
		
		vector <Event> ev_new;
		
		double dprob = UNSET;
		
		const auto &e_last = ev[ev.size()-1];
		auto c = e_last.c_after;
		if(c == UNSET){                 // Switches from sink to leave
			if(!e_last.observed && e_last.t > tlim){		
				auto tt = ev[ev.size()-2].t;
				if(tt > tmin) tmin = tt;
		
				dprob = -log(t_end-tmin);
			
				ev_new = ev;
				ev_new.pop_back();
			}
		}
		else{                            // Switches from leave to sink
			ev_new = ev;
			
			IndInfFrom inf_from;
			
			auto tr_gl = sp.cgl_tr_sink[c];
			if(tr_gl != UNSET){
				auto tt = e_last.t;
				if(tt > tmin) tmin = tt;
			
				if(tmin < t_end){
					const auto &tra = sp.tra_gl[tr_gl];

					auto type = tra.ev_type; 
					auto t = tmin + ran()*(t_end-tmin);
					
					dprob = log(t_end-tmin);

					auto e_new = ssp.get_event(type,i,tr_gl,UNSET,UNSET,tra.f,t,inf_from);
					ev_new.push_back(e_new);
				}
			}
		}
		
		if(dprob != UNSET && !events_near_div(ev_new,model.details)){
			if(pl){
				cout << endl << endl;
				cout << ind.name << endl;
				cout << "event old:" << endl; ssp.print_event(ind.ev);
				cout << "event new:" << endl; ssp.print_event(ev_new);
			}
								
			ntr++;
			
			auto gc = state.update_tree(p_prop,i,ev_new);
			if(gc.type == GENCHA_FAIL){
				nfa++;
			}
			else{
				auto like_ch = state.update_ind(p_prop,i,ev_new,UP_SINGLE);
			
				gc.update_like_ch(like_ch,dprob);
				
				auto al = calculate_al(like_ch,dprob);
		
				if(pl) cout << i << " " <<  al << " al" << endl;
			
				if(ran() < al){
					if(pl) cout << " accept" << endl;
					nac++;
					state.add_like(like_ch);
					state.gen_change_update(gc); 
					if(sp.trans_tree) state.update_popnum_ind(p_prop,i);
				}
				else{ 
					if(pl) cout << "reject" << endl;
					state.restore_back();
				}

				if(pl) state.check("switch leave sink");
			}
		}
	}
}


/// Tries to correct any observed transition events
void Proposal::correct_obs_trans_events(State &state)
{
	if(all_events_correct) return;
	
	const auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	
	all_events_correct = true;
	for(auto i = 0u; i < sp.nindividual_obs; i++){
		const auto &indd = sp.individual[i];
		const auto &ev = ssp.individual[i].ev;
		
		auto e = 0u;
		for(auto m = 0u; m < indd.obs.size(); m++){
			const auto &ob = indd.obs[m];
			switch(ob.type){
			case OBS_SINK_EV: case OBS_TRANS_EV:
				{
					auto t = ob.t;
					while(e < ev.size() && ev[e].t < t) e++;
					if(!(e < ev.size() && ev[e].t == t)){
						switch(ob.type){
						case OBS_SINK_EV: enforce_add_sink(i,t,state); break;
						case OBS_TRANS_EV: enforce_add_trans(i,ob,state); break;
						default: break;
						}
						all_events_correct = false;
					}
				}
				break;	
			default: break;
			}
		}
	}
}


/// Tries to put a source at a given time
void Proposal::enforce_add_sink(unsigned int i, double t, State &state)
{
	const auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	
	const auto &ev = ssp.individual[i].ev;
	
	vector <Event> ev_new;
	auto e = 0u;
	auto c = UNSET;
	while(e < ev.size() && ev[e].t < t){
		ev_new.push_back(ev[e]);
		c = ev[e].c_after;
		e++;
	}
	
	if(c != UNSET){
		auto tr = sp.cgl_tr_sink[c];
		if(tr != UNSET){
			const auto &tra = sp.tra_gl[tr];
			auto type = tra.ev_type;
		
			IndInfFrom inf_from;
			auto e = ssp.get_event(type,i,tr,UNSET,UNSET,tra.f,t,inf_from);
			
			ev_new.push_back(e);
			basic_ind_update(i,ev_new,state);
		}
	}
}


/// Tries to put a source at a given time
void Proposal::enforce_add_trans(unsigned int i, const ObsData &ob, State &state)
{
	auto &ssp = state.species[p_prop];	
	const auto &ev = ssp.individual[i].ev;
	
	auto t = ob.t;
	
	auto e_mid = 0u; while(e_mid < ev.size() && ev[e_mid].t) e_mid++;

	for(auto e = e_mid; e < ev.size(); e++){
		if(try_move_obs_trans_event(i,e,t,ob,ev,state)) return;
	}
	
	for(int e = int(e_mid)-1; e >= 0; e--){
		if(try_move_obs_trans_event(i,e,t,ob,ev,state)) return;
	}
}


/// Tries to move an event to an observed transition
bool Proposal::try_move_obs_trans_event(unsigned int i, unsigned int e, double t, const ObsData &ob, const vector <Event> &ev, State &state)
{
	const auto &eve = ev[e];
	if(eve.type != NM_TRANS_EV && eve.type != M_TRANS_EV) return false;
	
	auto trg = eve.tr_gl;
	
	auto &ssp = state.species[p_prop];
	
	if(ssp.get_trans_obs_prob(trg,ob) == 0) return false; 
	
	if(eve.t > t){
		auto e_before = get_event_before(e,ev);
		if(e_before != UNSET) 
		if(ev[e_before].t > t) return false;
	}
	else{
		auto e_after = get_event_after(e,ev);
		if(e_after != UNSET) 
		if(ev[e_after].t < t) return false;
	}
	
	auto ev_new = ev;
						
	ssp.move_event(ev_new,e,t);
		
	auto ee = 0u; while(ee < ev_new.size() && ev_new[ee].t < t) ee++;
	if(ee == ev_new.size()) emsg("Proble with ev new");
	ev_new[ee].observed = true;	
	
	basic_ind_update(i,ev_new,state);

	return true;
}


/// This provides a basic update of an individual to a new 
void Proposal::basic_ind_update(unsigned int i, vector <Event> &ev_new, State &state)
{
	auto pl = false;
	
	auto gc = state.update_tree(p_prop,i,ev_new);
	if(gc.type == GENCHA_FAIL){
	}
	else{
		auto like_ch = state.update_ind(p_prop,i,ev_new,UP_SINGLE);
	
		double dprob = 0;
		gc.update_like_ch(like_ch,dprob);
		
		auto al = calculate_al(like_ch,dprob);

		if(pl) cout << i << " " <<  al << " al" << endl;

		if(ran() < al){
			if(pl) cout << " accept" << endl;
			state.add_like(like_ch);
			state.gen_change_update(gc); 
			if(model.species[p_prop].trans_tree) state.update_popnum_ind(p_prop,i);
		}
		else{ 
			if(pl) cout << "reject" << endl;
			state.restore_back();
		}

		if(pl) state.check("basic ind");
	}
}


/// Updates the system based on the proposal type
void Proposal::update(State &state)
{
	if(ran() > prop_prob) return;
	
	number++;
	
	timer[PROP_TIMER] -= clock();
	switch(type){
	case TRANS_TREE_PROP: trans_tree(state); break;
	case TRANS_TREE_SWAP_INF_PROP: trans_tree_swap_inf(state); break;
	case TRANS_TREE_MUT_PROP: trans_tree_mut(state); break;
	case TRANS_TREE_MUT_LOCAL_PROP: trans_tree_mut_local(state); break;
	case PARAM_PROP: MH(state); break;
	case PAR_EVENT_FORWARD_PROP: param_event_joint(FORWARD,state); break;
	case PAR_EVENT_FORWARD_SQ_PROP: param_event_joint(FORWARD_SQ,state); break;
	case PAR_EVENT_BACKWARD_SQ_PROP: param_event_joint(BACKWARD_SQ,state); break;
	case IND_EVENT_TIME_PROP: MH_event(state); break;
	case IND_MULTI_EVENT_PROP: MH_multi_event(state); break;
	case IND_EVENT_ALL_PROP: MH_event_all(state); break;
	case IND_OBS_SWITCH_ENTER_SOURCE_PROP: switch_enter_source(state); break;
	case IND_OBS_SWITCH_LEAVE_SINK_PROP: switch_leave_sink(state); break;
	case CORRECT_OBS_TRANS_PROP: correct_obs_trans_events(state); break;
	case IND_LOCAL_PROP: MH_ind_local(state); break;
	case IND_OBS_SAMP_PROP: sample_ind_obs(state); break;
	case IND_OBS_RESIM_PROP: resimulate_ind_obs(state); break;
	case IND_OBS_RESIM_SINGLE_PROP: resimulate_single_ind_obs(state); break;
	case IND_UNOBS_RESIM_PROP: resimulate_ind_unobs(state); break;
	case IND_ADD_REM_PROP: add_rem_ind(state); break;
	case IND_ADD_REM_TT_PROP: add_rem_tt_ind(state); break;
	case MBP_PROP: mbp(state); break;
	case MBPII_PROP: mbp(state); break;
	case MBP_IC_POP_PROP: mbp(state); break;
	case MBP_IC_POPTOTAL_PROP: mbp(state); break;
	case MBP_IC_RESAMP_PROP: mbp(state); break;
	case INIT_COND_FRAC_PROP: init_cond_frac(state); break;
	case IE_PROP: MH_ie(state); break;
	case IE_VAR_PROP: MH_ie_var(state); break;
	case IE_COVAR_PROP: MH_ie_covar(state); break;
	case IE_VAR_CV_PROP: MH_ie_var_cv(state); break;
	case POP_ADD_REM_LOCAL_PROP: pop_add_rem_local(state); break;
	case POP_MOVE_LOCAL_PROP: pop_move_local(state); break;
	case POP_IC_LOCAL_PROP: pop_ic_local(state); break;
	case POP_END_LOCAL_PROP: pop_end_local(state); break;
	case POP_SINGLE_LOCAL_PROP: pop_single_local(state); break;
	case POP_IC_PROP: pop_ic(state); break;	
	case POP_IC_SWAP_PROP: pop_ic_swap(state); break;	
	}	
		 
	timer[PROP_TIMER] += clock();
}

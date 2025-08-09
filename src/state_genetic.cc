// Functions related to genetic observations

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
using namespace std;

#include "state.hh"
#include "utils.hh"

/// Used to time order infection nodes
bool InfNode_ord (const InfNode &in1, const InfNode &in2)                      
{ return (in1.tdiv_start < in2.tdiv_start); };  


/// Used to time order infection events
bool InfEv_ord (const InfEvent &ev1, const InfEvent &ev2)                      
{ return (ev1.tdiv < ev2.tdiv); };  


/// Given an initial event sequence works out how to adapt the transmission tree
GenChange State::update_tree(unsigned int p, unsigned int i, vector <Event> &ev_new)
{
	auto pl = false;
	
	const auto &sp = model.species[p];
	
	if(sp.trans_tree == false) return GenChange(NO_GENETIC_CHANGE);
	
	auto &inf_node = genetic_value.inf_node;
	
	const auto &ev_old = species[p].individual[i].ev;
	
	vector <unsigned int> node_old;
	
	const auto &ssp =species[p];
	
	if(pl){
		cout << endl << endl << "BEGIN" << endl;
		ssp.print_event(ev_old);
		ssp.print_event(ev_new);
	}
	
	for(const auto &ev : ev_old){
		auto n = ev.inf_node_ref;
		if(n != UNSET) node_old.push_back(n);
	}
	
	if(node_old.size() > 1) return GenChange(GENCHA_FAIL);
	
	vector <unsigned int> nodee_new;
	
	for(auto e = 0u; e < ev_new.size(); e++){
		const auto &ev = ev_new[e];
		switch(ev.type){
		case ENTER_EV:
			if(sp.comp_gl[ev.c_after].infected == true){
				nodee_new.push_back(e);
			}
			break;
			
		case M_TRANS_EV:
			if(sp.tra_gl[ev.tr_gl].infection.type == TRANS_INFECTION){
				nodee_new.push_back(e);
			}
			break;
			
		default: break;
		}
	}

	if(nodee_new.size() > 1) return GenChange(GENCHA_FAIL);
	
	switch(node_old.size()){
	case 0: 
		switch(nodee_new.size()){
		case 0:                                        // No infection -> no infection
			return GenChange(NO_GENETIC_CHANGE);
	
		case 1:                                        // No infection -> infection
			{		
				auto e_add = nodee_new[0];
				
				auto &ev = ev_new[e_add];
				
				auto probif = sample_infection_source(ev,p);
				
				if(probif == UNSET) return GenChange(GENCHA_FAIL); // Fails to find infection
					
				if(testing){
					auto prob = prob_infection_source(ev,p);
					if(dif(prob,probif,dif_thresh)){
						emsg("Source sampler problem");
					}
				}
			
				auto ina = get_inf_node_alteration(p,i,UNSET,ev_new,e_add);
				if(ina.possible == true){ 
					auto gc = gen_change(ADD_NODE,UNSET,ina,ev.ind_inf_from);
				
					ev.inf_node_ref = inf_node.size();
					
					gc.node_add.p = p;
					gc.node_add.i = i;
					gc.node_add.e = e_add;
					
					gc.probif += probif;
					return gc;
				}			
			}
			break;
		}
		break;
		
	case 1:
		switch(nodee_new.size()){
		case 0:                                         // Infection -> no infection
			{		
				auto n = node_old[0];
				const auto &in = inf_node[n];
				
				if(in.inf_ev.size() == 0){         // Only remove if no observations or onward infections
					const auto &ev = species[in.p].individual[in.i].ev[in.e];
					
					auto probfi = prob_infection_source(ev,in.p);
					
					auto gc = gen_change(REM_NODE,n,InfNodeAlter(),ev.ind_inf_from);
					
					gc.probfi += probfi;
					return gc;
				}
			}
			break;
	
		case 1:                                         // Infection -> Infection
			{		
				auto n = node_old[0];
				const auto &in_old = inf_node[n];
				
				auto e_add = nodee_new[0];
				auto &ev = ev_new[e_add];
				
				auto t_old = in_old.tdiv_start;
				auto t_new = ev.tdiv;
		
				auto ina = get_inf_node_alteration(p,i,n,ev_new,e_add);
				if(ina.possible == true){ 
					if(ina.unchanged){
						ev = ssp.individual[in_old.i].ev[in_old.e];
						return GenChange(NO_GENETIC_CHANGE);
					}
					
					auto probif = 0.0, probfi = 0.0;
				
					NodeRef nr;
					
					const auto &from = in_old.from;
					if(from.node == ENTER_INF || from.node == OUTSIDE_INF){
						nr = from;
					}
					else{
						if(t_old == t_new){
							nr = from;
						}
						else{
							if(t_old < t_new){
								auto nr_start = in_old.from; nr_start.index++;
								nr = pass_down(nr_start,t_new,probif);
							}
							else{
								nr = pass_up(in_old.from,t_new,probfi);
							}
						}
					}
					
					auto nn = nr.node;
					if(nn == UNSET) return GenChange(GENCHA_FAIL);
					
					IndInfFrom iif_prop;
					
					if(nn == ENTER_INF || nn == OUTSIDE_INF){
						auto nn_new = ENTER_INF;
						if(t_new != 0) nn_new = OUTSIDE_INF;
						if(nn_new != nn && nn_new == OUTSIDE_INF){
							if(trg_contains_outside(in_old.p,ev.tr_gl) == false) return GenChange(GENCHA_FAIL);	
						}
					
						iif_prop.p = nn_new;
						iif_prop.i = UNSET;
						iif_prop.pref = UNSET;
						iif_prop.po = UNSET;
						if(nn_new == OUTSIDE_INF) iif_prop.w = 1;
						else iif_prop.w = UNSET;		
					}
					else{
						const auto &in_new = inf_node[nn];
					
						// Checks if new connection point in within infectious range
						if(t_new < in_new.tdiv_start || t_new > round_up(in_new.tdiv_rec)){
							return GenChange(GENCHA_FAIL);
						}
						
						// Works out IndInfFrom for new position on tree
						
						auto pp = in_new.p;
						
						iif_prop.p = pp;
						iif_prop.i = in_new.i;
					
						// Get the c of new infecting individual
						auto t_inf = round_down(t_new);
						const auto &ev_inf = species[in_new.p].individual[in_new.i].ev;
						
						auto e = in_new.e;
						
						auto c = UNSET;
						while(e < ev_inf.size() && ev_inf[e].tdiv <= t_inf){ c = ev_inf[e].c_after; e++;}
						if(c == UNSET) return GenChange(GENCHA_FAIL);
						
						// Works out pref from c
						
						const auto &tra = sp.tra_gl[ev.tr_gl];
						auto m = tra.markov_eqn_ref;
						const auto &me = sp.markov_eqn[m];
						const auto &eq = model.eqn[me.eqn_ref];
						auto pref = eq.comp_pref_convert[pp][c];
						if(pref == UNSET) return GenChange(GENCHA_FAIL);
						//if(pref == UNSET) emsg("Could not get pref"); 
						
						iif_prop.pref = pref;
						iif_prop.po = eq.pop_ref[pref];
						iif_prop.w = get_w_from_indinffrom(iif_prop);
					}
					ev.ind_inf_from = iif_prop;
					ev.inf_node_ref = n;
					
					auto gc = gen_change(REGRAFT_NODE,node_old[0],ina,iif_prop);
				
					gc.probif += probif;
					gc.probfi += probfi;
					return gc;
				}			
			}
			break;
		}
		break;
	}
	
	return GenChange(GENCHA_FAIL);
}


/// Determines is a transition contains a term for outside infection
bool State::trg_contains_outside(unsigned int p, unsigned int tr_gl) const
{
	const auto &sp = model.species[p];
	
	const auto &tra = sp.tra_gl[tr_gl];
	if(tra.infection.type != TRANS_INFECTION) emsg("Should be transition");
	const auto &me = sp.markov_eqn[tra.markov_eqn_ref];
	const auto &eq = model.eqn[me.eqn_ref];
	const auto &lin = eq.linearise;
	if(lin.no_pop_calc.size() == 0) return false;
	return true;
}


/// Samples the infection source for a new infection 
double State::sample_infection_source(Event &ev, unsigned int p) const
{
	auto probif = 0.0;
	
	auto &iif = ev.ind_inf_from;
		
	if(ev.type == ENTER_EV){
		if(false){
			cout << "ENTER INF <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
		}
		iif.p = ENTER_INF; 
		iif.i = UNSET;
		iif.pref = UNSET;
		iif.po = UNSET;
		iif.w = UNSET;
	}
	else{
		const auto &sp = model.species[p];
		auto ti = get_ti(ev.tdiv);
		
		const auto &tra = sp.tra_gl[ev.tr_gl];
		
		if(tra.infection.type != TRANS_INFECTION) emsg("Should be transition");
		
		const auto &me = sp.markov_eqn[tra.markov_eqn_ref];

		const auto &eq = model.eqn[me.eqn_ref];
			
		const auto &lin = eq.linearise;
			
		if(!lin.on) emsg("Linearisation should be on");
			
		// First selects if infection comes from outside or from which population
			
		auto Npop = eq.pop_ref.size();
			
		auto j = UNSET;
			
		if(lin.multi_source){        // Samples from available sources (either populations of from outside)		
			auto ss = eq.setup_source_sampler(ti,popnum_t[ti],param_val,spline_val);	
			j = ss.sample_inf_source();
			if(j == UNSET) return UNSET;
			
			probif += ss.prob_inf_source(j);
		}
		else{
			if(testing){
				if(Npop > 1) emsg("Npop wrong");
				if(Npop == 1 && lin.no_pop_calc.size() > 0) emsg("Npop wrong");
			}
			j = 0;
		}
		
		if(j == Npop){
			if(false){
				cout << "OUTSIDE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
			}
			iif.p = OUTSIDE_INF; 
			iif.i = UNSET;
			iif.pref = UNSET;
			iif.po = UNSET;
			iif.w = 1;
		}
		else{
			auto pr = eq.pop_ref[j];
			const auto &pop = model.pop[pr];

			const auto &pi = popnum_ind[ti][pr];
			
			auto si = pi.size();
			if(si == 0) return UNSET;
			
			if(pop.ind_variation){ // Samples in proportion to weight
				auto sum = 0.0;
				vector <double> sum_st(si);
				for(auto k = 0u; k < si; k++){
					sum += pi[k].w;
					sum_st[k] = sum;
				}
				auto z = ran()*sum;
				auto k = 0u; while(k < si && z > sum_st[k]) k++;
				if(k == si) emsg("Problem sampling");
				
				const auto &pi_sel = pi[k];
				
				probif += log(pi_sel.w/sum);
				iif.i = pi_sel.i;
				iif.w = pi_sel.w;
			}
			else{
				probif += log(1.0/si);
				const auto &pi_sel = pi[(unsigned int)(ran()*si)];
				
				iif.i = pi_sel.i;
				iif.w = pi_sel.w;
			}
		
			iif.p = pop.sp_p;
			
			iif.pref = j;
			iif.po = pr;
		}
	}
	
	return probif;
}


/// The probability of an infection source
double State::prob_infection_source(const Event &ev, unsigned int p) const
{
	if(ev.type == ENTER_EV) return 0;
		
	const auto &sp = model.species[p];
	const auto &ssp = species[p];
	
	auto ti = get_ti(ev.tdiv);
	auto &iif = ev.ind_inf_from;

	const auto &tra = sp.tra_gl[ev.tr_gl];
	if(tra.infection.type != TRANS_INFECTION) emsg("Should be transition");
	
	auto m = tra.markov_eqn_ref;
	const auto &me = sp.markov_eqn[m];
	const auto &mev = ssp.markov_eqn_vari[m];
	
	const auto &eq = model.eqn[me.eqn_ref];
		
	const auto &lin = eq.linearise;
		
	if(!lin.on) emsg("Linearisation should be on");
		
	// First selects if infection comes from outside or from which population
		
	auto probfi = 0.0;
	
	auto pref = iif.pref;
		
	if(lin.multi_source){        // Prop from available sources (either populations or outside)		
		auto value = mev.div[ti].value;
		if(value < TINY) probfi -= LARGE;
		else{		
			double val;
			if(pref == UNSET) val = eq.calculate_calculation(lin.no_pop_calc,ti,param_val,spline_val);
			else{
				auto pr = eq.pop_ref[pref];
				val = popnum_t[ti][pr]*eq.calculate_calculation(lin.pop_grad_calc[pref],ti,param_val,spline_val);
			}
			probfi += log(val/value);
		}
	}
		
	if(pref != UNSET){
		auto pr = eq.pop_ref[pref];
		const auto &pop = model.pop[pr];

		const auto &pi = popnum_ind[ti][pr];
		
		auto si = pi.size();
		if(si == 0){
			return 0;
			//emsg("Cannot find infector2"); LOOK AT
		}
		
		if(pop.ind_variation){
			auto sum = 0.0; for(auto k = 0u; k < si; k++) sum += pi[k].w;
			probfi += log(iif.w/sum);
		}
		else{
			probfi += log(1.0/si);
		}
	}
	
	return probfi;
}
	
	
/// Passes down the tree until a certain time
NodeRef State::pass_down(const NodeRef &nref, double t, double &probif) const
{
	auto n = nref.node;
	const auto &inf_ev = genetic_value.inf_node[n].inf_ev;
	
	for(auto e = nref.index; e < inf_ev.size(); e++){
		const auto &ev = inf_ev[e];
		if(ev.tdiv > t){
			NodeRef nr; nr.node = n; nr.index = e; return nr;
		}
		if(ev.type == INFECT_OTHER){
			probif += LOG_HALF;
			if(ran() < 0.5){
				NodeRef nr_new; nr_new.node = ev.index; nr_new.index = 0;
				return pass_down(nr_new,t,probif);
			}
		}
	}
	
	NodeRef nr; nr.node = n; nr.index = inf_ev.size(); return nr;
}


/// Passes up the tree until a certain time is reached
NodeRef State::pass_up(const NodeRef &nref, double t, double &probfi) const
{
	const auto &in = genetic_value.inf_node[nref.node];
	const auto &inf_ev = in.inf_ev;
	
	for(int e = nref.index; e >= 0; e--){
		double t_last; if(e > 0) t_last = inf_ev[e-1].tdiv; else t_last = in.tdiv_start;
		
		if(t_last < t){
			NodeRef nr; nr.node = nref.node; nr.index = e; return nr;
		}
		if(inf_ev[e].type == INFECT_OTHER){
			probfi += LOG_HALF;
		}
	}
	
	auto n = in.from.node;
	if(n != ENTER_INF && n != OUTSIDE_INF){
		return pass_up(in.from,t,probfi);
	}
	
	NodeRef nr; nr.node = UNSET; nr.index = UNSET; return nr;
}


/// Works out how the node is altered for the new timeline
InfNodeAlter State::get_inf_node_alteration(unsigned int p, unsigned int i, unsigned int n, const vector <Event> &ev_new, unsigned int e_add) const
{			
	const auto &sp = model.species[p]; 
	const auto &gen_data = model.genetic_data;

	auto on = gen_data.on; if(i >= sp.nindividual_obs) on = false;
	
	auto t_start = ev_new[e_add].tdiv;
	double t_rec = model.details.T;;
	for(auto ee = e_add+1; ee < ev_new.size(); ee++){
		const auto &ev = ev_new[ee];
		if(ev.type == M_TRANS_EV){
			if(sp.tra_gl[ev.tr_gl].infection.type == TRANS_RECOVERY){
				t_rec = ev.tdiv;
				break;
			}
		}
	}

	InfNodeAlter ina;
	ina.e_add = e_add;
	
	ina.tdiv_start = t_start;
	ina.tdiv_rec = t_rec;

	if(n != UNSET){
		const auto &in_old = genetic_value.inf_node[n];
	
		if(t_start == in_old.tdiv_start && t_rec == in_old.tdiv_rec){
			ina.unchanged = true;
			return ina;
		}	
		
		// Makes sure that any onward infections are within the time span
		auto t_rec_ro = round_up(t_rec);
		auto t_start_ro = round_up(t_start);
		
		const auto &inf_ev = in_old.inf_ev;
		
		for(auto k = 0u; k < inf_ev.size(); k++){
			const auto &iev = inf_ev[k];
			if(iev.type == INFECT_OTHER){
				auto t = iev.tdiv;
				if(t < t_start_ro || t > t_rec_ro){
				//if(t < t_start || t > t_rec_ro){
					ina.possible = false;
				}
			}
		}
	
		if(on){
			const auto &obs = gen_data.ind_gen_obs[p][i];
			
			auto t_start_old = in_old.tdiv_start;
			if(t_start < t_start_old){                    // Looks to add observations at the beginning
				auto k = 0u; 
				while(k < obs.size() && obs[k].tdiv < t_start) k++;

				while(k < obs.size() && obs[k].tdiv <= t_start_old){
					ina.obs_add_begin.push_back(obs[k].m);
					k++;
				}
			}
			else{                                         // Impossible if removing observations
				auto k = 0u; 
				while(k < obs.size() && obs[k].tdiv < t_start_old) k++;
				if(k < obs.size() && obs[k].tdiv <= t_start){
					ina.possible = false;
				}
			}
			
			auto t_rec_old = in_old.tdiv_rec;
			if(t_rec > t_rec_old){                        // Looks to add observations at the end
				auto k = 0u; 
				while(k < obs.size() && obs[k].tdiv < t_rec_old) k++;

				while(k < obs.size() && obs[k].tdiv <= t_rec){
					ina.obs_add_end.push_back(obs[k].m);
					k++;
				}
			}
			else{                                         // Impossible if removing observations
				auto k = 0u; 
				while(k < obs.size() && obs[k].tdiv <= t_rec) k++;
				if(k < obs.size() && obs[k].tdiv < t_rec_old) ina.possible = false;
			}
		}
	}
	else{ // For newly introduced infections use complete time span
		if(on){
			const auto &obs = gen_data.ind_gen_obs[p][i];
		
			auto k = 0u; 
			while(k < obs.size() && obs[k].tdiv < t_start) k++;

			while(k < obs.size() && obs[k].tdiv <= t_rec){
				ina.obs_add_end.push_back(obs[k].m);
				k++;
			}
		}
	}

	return ina;
}


/// Sets up the transition tree
void State::setup_transtree()
{
	auto &inf_node = genetic_value.inf_node;
	inf_node.clear();
	
	const auto &gen_data = model.genetic_data;
	
	// Sets up the nodes
	
	for(auto p = 0u; p < model.nspecies; p++){ 
		const auto &sp = model.species[p];
		auto &ssp = species[p];

		for(auto i = 0u; i < ssp.individual.size(); i++){
			auto &ind = ssp.individual[i];
			
			for(auto e = 0u; e < ind.ev.size(); e++){
				auto &eve = ind.ev[e];
				eve.inf_node_ref = UNSET;
				
				if(eve.type == ENTER_EV){
					if(sp.comp_gl[eve.c_after].infected == true){
						InfNode in; in.tdiv_start = eve.tdiv; in.p = p; in.i = i; in.e = e; 
						in.tdiv_rec = find_t_rec(in);
						in.from.node = UNSET; in.from.index = UNSET;
						inf_node.push_back(in);
						
						if(eve.ind_inf_from.p != ENTER_INF){
							emsg("ind_inf_from.p should not be seta");
						}
					}
					else{
						if(eve.ind_inf_from.p != UNSET) emsg("ind_inf_from.p should not be setb");
					}
				}
				else{
					if(eve.type == MOVE_EV){
						const auto &eve_last = ind.ev[e-1];
						if(sp.comp_gl[eve_last.c_after].infected == false &&
						   sp.comp_gl[eve.c_after].infected == true){
							
							InfNode in; in.tdiv_start = eve.tdiv; in.p = p; in.i = i; in.e = e; 
							in.tdiv_rec = find_t_rec(in);
							in.from.node = UNSET; in.from.index = UNSET;
							inf_node.push_back(in);
						
							if(eve.ind_inf_from.p != ENTER_INF){
								emsg("ind_inf_from.p should not be seta");
							}
						}
						else{
							if(eve.ind_inf_from.p != UNSET) emsg("ind_inf_from.p should not be setb");
						}
					}
					else{
						if(eve.ind_inf_from.p != UNSET){
							InfNode in; in.tdiv_start = eve.tdiv; in.p = p; in.i = i; in.e = e; 
							in.tdiv_rec = find_t_rec(in);
							in.from.node = UNSET; in.from.index = UNSET;
							inf_node.push_back(in);
						}
					}
				}
			}
		}
	}
	
	sort(inf_node.begin(),inf_node.end(),InfNode_ord);    
	
	for(auto n = 0u; n < inf_node.size(); n++){
		const auto &in = inf_node[n];
		species[in.p].individual[in.i].ev[in.e].inf_node_ref = n;
	}
	
	// Adds which nodes infect which other ones and observations
	genetic_value.nobs_not_infected = 0;
	genetic_value.obs_node_ref.clear();
	genetic_value.inf_origin.clear();
	
	// Initialises ind_node_ref[p][ind]
	vector < vector <unsigned int> > ind_node_ref;
	ind_node_ref.resize(model.nspecies);
	for(auto p = 0u; p < model.nspecies; p++){ 
		auto &ssp = species[p];	
		ind_node_ref[p].resize(ssp.individual.size(),UNSET);
	}
	
	// Adds observations to inf_ev timelines on nodes
	// Adds infection to inf_ev timelines on nodes
	// Adds origins to inf_origin
	auto m = 0u;
	for(auto k = 0u; k <= inf_node.size(); k++){
		double t;
		if(k < inf_node.size()) t = inf_node[k].tdiv_start; else t = LARGE;
		
		while(m < gen_data.obs.size() && gen_data.obs[m].tdiv < t){  // Adds observations
			const auto &obs = gen_data.obs[m];
			
			NodeRef nr;
			
			auto kk = ind_node_ref[obs.p][obs.i];
			nr.node = kk;
			if(kk == UNSET){                               // Individual is not infected
				genetic_value.nobs_not_infected++;
				nr.index = UNSET;
			}
			else{
				auto &in2 = inf_node[kk];
				if(obs.tdiv < in2.tdiv_start || obs.tdiv > in2.tdiv_rec){// Individual not infectious during observation
					genetic_value.nobs_not_infected++;
					nr.node = UNSET;
					nr.index = UNSET;
				}
				else{
					InfEvent ie; ie.type = GENETIC_OBS; ie.tdiv = obs.tdiv; ie.index = m; ie.mut_num = UNSET;
					
					nr.index = in2.inf_ev.size();
					in2.inf_ev.push_back(ie);
				}
			}
			
			genetic_value.obs_node_ref.push_back(nr);
			m++;
		}
		
		if(k < inf_node.size()){                        // Creates origin
			auto &in = inf_node[k];
		
			ind_node_ref[in.p][in.i] = k;
			const auto &eve = species[in.p].individual[in.i].ev[in.e];
		
			auto pp = eve.ind_inf_from.p;
			if(pp == UNSET) emsg("Should not be unset");
		
			if(pp == ENTER_INF || pp == OUTSIDE_INF){	
				in.from.node = pp;
				in.from.index = genetic_value.inf_origin.size();
			
				InfOrigin io; io.node = k; io.mut_num = UNSET;
		
				genetic_value.inf_origin.push_back(io);
			}
			else{
				auto kk = ind_node_ref[eve.ind_inf_from.p][eve.ind_inf_from.i];
			
				if(kk == UNSET){
					emsg("Problem with inf_from");
				}
	
				auto &in2 = inf_node[kk];
				
				in.from.node = kk;
				in.from.index = in2.inf_ev.size();
				
				InfEvent ie; ie.type = INFECT_OTHER; ie.tdiv = eve.tdiv; ie.index = k; ie.mut_num = UNSET;
				in2.inf_ev.push_back(ie);
			}
		}
	}
		
	if(false){
		for(auto n = 0u; n < inf_node.size(); n++){
			cout << "node " << n << endl;
			const auto &no = inf_node[n];
			cout << no.tdiv_start << " " << no.tdiv_rec << " time range" << endl;
			
			for(const auto &iev : no.inf_ev){
				cout << iev.tdiv << ",";
			}
			cout << endl;
		}
		emsg("done");
	}
}


/// Finds the recovery time for an infection node
double State::find_t_rec(const InfNode &in) const
{
	const auto &sp = model.species[in.p];
	const auto &ev = species[in.p].individual[in.i].ev;
	 
	for(auto e = in.e; e < ev.size(); e++){
		if(ev[e].type == M_TRANS_EV && sp.tra_gl[ev[e].tr_gl].infection.type == TRANS_RECOVERY){
			return ev[e].tdiv;
		}
	}

	return model.details.T;
}


/// Sets the values used for the genetic parameters
vector <double> State::set_genetic_param()
{
	vector <double> store;
	store.push_back(genetic_value.mut_rate);
	store.push_back(genetic_value.seq_var);
	
	const auto &gen_data = model.genetic_data;
	genetic_value.mut_rate = model.eqn[gen_data.mut_rate.eq_ref].calculate_param_only(param_val);
	
	genetic_value.seq_var = model.eqn[gen_data.seq_var.eq_ref].calculate_param_only(param_val);

	return store;
}


/// Samples quantities in sample_genetic_value()
void State::sample_genetic_value()
{
	auto mut_rate = genetic_value.mut_rate*model.details.dt;
	
	auto seq_var = genetic_value.seq_var;
	
	auto &inf_node = genetic_value.inf_node;
	
	string warn;
	
	// Samples mutation events down the timelines for ind_nodes
	for(auto k = 0u; k < inf_node.size(); k++){
		auto &in = inf_node[k];
		
		auto t = in.tdiv_start;
		for(auto e = 0u; e < in.inf_ev.size(); e++){
			auto &ev = in.inf_ev[e];
			
			auto val = poisson_sample(mut_rate*(ev.tdiv-t),warn);
			if(val == UNSET) emsg("Problem sampling from a Poisson distribution: "+warn);
			ev.mut_num = val;
			t = ev.tdiv;
		}
	}
	
	for(auto &io : genetic_value.inf_origin){
		auto &in = inf_node[io.node];
		auto val = poisson_sample(seq_var+mut_rate*in.tdiv_start,warn);
		if(val == UNSET) emsg("Problem sampling from a Poisson distribution: "+warn);
		io.mut_num = val;
	}
	
	// Calculates the genetic difference between different observed sequences
	genetic_value.gen_dif = calculate_gen_dif();
	
	if(false){
		for(auto &io : genetic_value.inf_origin){
			cout << "ORIGIN" << " " << io.node << " " << io.mut_num << endl;
		}
			
		for(auto k = 0u; k < inf_node.size(); k++){
			const auto &in = inf_node[k];
			cout << k << " " << in.tdiv_start << " " << model.species[in.p].individual[in.i].name 
					 << "  from: " << in.from.node <<" " << in.from.index << endl;
			
			for(auto e = 0u; e < in.inf_ev.size(); e++){
				const auto &ev = in.inf_ev[e];
				cout << "  " << ev.tdiv << " " << ev.mut_num << "  ";
				switch(ev.type){
				case INFECT_OTHER: cout << "infect other "; break;
				case GENETIC_OBS: cout << "observation "; break;
				}
				cout << ev.index << ",   ";
			}
			cout << endl;
			cout << endl;
		}
	}
}


/// Calculates the difference in genetic sequences
vector < vector <unsigned int> > State::calculate_gen_dif()
{
	const auto &gen_data = model.genetic_data;

	vector < vector <unsigned int> > gen_dif;
	
	auto nobs = gen_data.obs.size();
	
	vector <unsigned int> root_num;
	
	for(auto m = 0u; m < nobs; m++){
		//for(auto m = 5u; m < 6; m++){
	
		const auto &nr = genetic_value.obs_node_ref[m];
		
		vector <unsigned int> dif(nobs,UNSET);
		dif[m] = 0;
		auto r_num = UNSET;
		
		auto n = nr.node;
		
		if(n != UNSET && n != ENTER_INF && n != OUTSIDE_INF){			
			calculate_dif(n,nr.index+1,0,dif);
			calculate_dif_rev(n,nr.index,0,dif,r_num);
		}
		
		gen_dif.push_back(dif);
		root_num.push_back(r_num);
	}

	// Calculates observations across different origins

	const auto &io = genetic_value.inf_origin;
	auto norigin = io.size();
	
	vector < vector <unsigned int> > origin_obs_list; 
	origin_obs_list.resize(norigin);
	for(auto i = 0u; i < norigin; i++){
		find_all_obs(io[i].node,origin_obs_list[i]);
	}
	
	for(auto i = 0u; i < norigin; i++){
		for(auto j = 0u; j < norigin; j++){
			if(j != i){
				auto num = io[i].mut_num + io[j].mut_num;
				for(auto ob1 : origin_obs_list[i]){
					for(auto ob2 : origin_obs_list[j]){
						if(gen_dif[ob1][ob2] != UNSET) emsg("Should not be unset");
						gen_dif[ob1][ob2] = num+root_num[ob1]+root_num[ob2];
					}
				}
			}
		}
	}
	
	return gen_dif;
}


/// Calculates genetic difference going down events in a node
void State::calculate_dif(unsigned int n, unsigned int e_start, unsigned int num, vector <unsigned int> &dif) const
{
	const auto &inf_node = genetic_value.inf_node;
	const auto &in = inf_node[n];

	for(auto e = e_start; e < in.inf_ev.size(); e++){
		const auto &ev = in.inf_ev[e]; 

		num += ev.mut_num;
		
		auto i = ev.index;
		switch(ev.type){
		case INFECT_OTHER:
			calculate_dif(i,0,num,dif);
			break;
		
		case GENETIC_OBS:
			if(dif[i] != UNSET) emsg("Should not be gen obs");
			dif[i] = num;
			break;
		}
	}			
}


/// Calculates genetic difference going up events in a node
void State::calculate_dif_rev(unsigned int n, int e_start, unsigned int num,  vector <unsigned int> &dif, unsigned int &r_num) const
{
	const auto &inf_node = genetic_value.inf_node;
	const auto &in = inf_node[n];

	for(int e = e_start; e >= 0; e--){
		const auto &ev = in.inf_ev[e]; 
		auto i = ev.index;
		
		if(e != e_start){
			switch(ev.type){
			case INFECT_OTHER:
				calculate_dif(i,0,num,dif);
				break;
			
			case GENETIC_OBS:
				if(dif[i] != UNSET) emsg("Should not be genobs2");
				dif[i] = num;
				break;
			}
		}
		num += ev.mut_num;
	}			
	
	const auto &fr = in.from;
	
	if(fr.node == ENTER_INF || fr.node == OUTSIDE_INF) r_num = num;
	else{		
		calculate_dif(fr.node,fr.index+1,num,dif); 
		calculate_dif_rev(fr.node,fr.index,num,dif,r_num); 
	}
}


/// Finds all the observations from the beginning of a given node
void State::find_all_obs(unsigned int node, vector <unsigned int> &obs_list) const
{
	const auto &in = genetic_value.inf_node[node];
	for(auto e = 0u; e < in.inf_ev.size(); e++){
		const auto &ie = in.inf_ev[e];
		switch(ie.type){
		case INFECT_OTHER:
			find_all_obs(ie.index,obs_list);
			break;
			
		case GENETIC_OBS:
			obs_list.push_back(ie.index);
			break;
		}
	}
}	


/// Calculate the genetic process likelihood
double State::likelihood_genetic_process()
{
	const auto &gen_data = model.genetic_data;
	
	if(model.mode == SIM  || !gen_data.on) return 0;
	
	auto Li = 0.0;

	auto mut_rate = genetic_value.mut_rate*model.details.dt;
	auto seq_var = genetic_value.seq_var;
	
	for(const auto &in : genetic_value.inf_node){
		auto t = in.tdiv_start;
		for(auto e = 0u; e < in.inf_ev.size(); e++){
			auto &ev = in.inf_ev[e];
			Li += poisson_probability(ev.mut_num,mut_rate*(ev.tdiv-t));
			t = ev.tdiv;
		}
	}
	
	for(auto &io : genetic_value.inf_origin){
		const auto &in = genetic_value.inf_node[io.node];
		Li += poisson_probability(io.mut_num,seq_var+mut_rate*in.tdiv_start);
	}
	
	return Li;
}


/// Calculate the genetic process likelihood
double State::likelihood_genetic_obs(const vector < vector <unsigned int> > &gen_val)
{
	const auto &gen_data = model.genetic_data;
	
	if(model.mode == SIM  || !gen_data.on) return 0;
	
	auto N = gen_data.obs.size();
	
	const auto &gen_dif = gen_data.gen_dif;
	
	auto Li = 0.0;
	
	for(auto j = 0u; j < N; j++){
		for(auto i = 0u; i < N; i++){	
			double val1 =	gen_val[j][i], val2 = gen_dif[j][i];
			if(val1 != UNSET) Li += -(val1-val2)*(val1-val2);
		}
	}
	
	Li *= gen_data.obs_scale;
	
	Li += genetic_value.nobs_not_infected*GEN_OBS_WRONG;
	
	return Li;
}


/// Used to order global events
bool GlobalEvent_ord(const GlobalEvent &ge1, const GlobalEvent &ge2)                      
{ return (ge1.tdiv < ge2.tdiv); };  


/// Performs a transmission tree proposal
void State::trans_tree_proposal(const BurnInfo &burn_info, unsigned int &nac, unsigned int &ntr)
{
	auto pl = false;
	
	// Makes a global list of all events
	if(pl) cout << "Start trans tree proposal" << endl;
	
	vector <GlobalEvent> glob_ev;
	
	for(auto p = 0u; p < model.species.size(); p++){
		const auto &sp = model.species[p];
		if(sp.trans_tree){
			const auto &ssp = species[p];
			for(auto i = 0u; i < ssp.individual.size(); i++){	
				const auto &ind = ssp.individual[i];
				for(auto e = 0u; e < ind.ev.size(); e++){
					GlobalEvent ge; ge.tdiv = ind.ev[e].tdiv; ge.p = p; ge.i = i; ge.e = e;
					glob_ev.push_back(ge);
				}
			}
		}
	}
	
	sort(glob_ev.begin(),glob_ev.end(),GlobalEvent_ord);    

	// Initialises ind_pop
	vector < vector <IndPop> > ind_pop;
	ind_pop.resize(model.species.size());
	for(auto p = 0u; p < model.species.size(); p++){
		const auto &sp = model.species[p];
		if(sp.type == INDIVIDUAL){
			const auto &ssp = species[p];
	
			ind_pop[p].resize(ssp.individual.size());
			for(auto i = 0u; i < ssp.individual.size(); i++){	
				ind_pop[p][i].c = UNSET;
			}
		}			
	}
	
	// Goes through time keeping track of which individuals are in which populations

	vector < vector <unsigned int> > pop_ind;
	pop_ind.resize(model.pop.size());
	
	auto egl = 0u;
	for(auto ti = 0u; ti < T; ti++){
		if(pl) cout << ti << "ti" << endl;
	
		// Updates populations
		//while(egl < glob_ev.size() && glob_ev[egl].tdiv <= tdiv){
		while(egl < glob_ev.size() && glob_ev[egl].tdiv <= ti){
			const auto &ge = glob_ev[egl];
			auto p = ge.p;
			auto i = ge.i;
			
			const auto &sp = model.species[p];
			const auto &ssp = species[p];
			
			const auto &ev = ssp.individual[i].ev[ge.e];
			
			auto &ip = ind_pop[p][i];
			
			auto c_before = ip.c;
			if(c_before != UNSET){
				for(const auto &pr : ip.pop_ref){          // Removes existing population references
					auto po = pr.po;
					auto ind = pr.index;
					if(pop_ind[po][ind] != i) emsg("Should be i");
					
					auto N = pop_ind[po].size();
					if(ind+1 < N){
						auto ii = pop_ind[po][N-1];
						pop_ind[po][ind] = ii;
						
						auto &iip = ind_pop[p][ii];
						auto k = 0u; while(k < iip.pop_ref.size() && !(iip.pop_ref[k].po == po && iip.pop_ref[k].index == N-1)) k++;
						if(k == iip.pop_ref.size()) emsg("cannot find");
						
						iip.pop_ref[k].index = ind;
					}
					pop_ind[po].pop_back();
				}
				
				ip.pop_ref.clear();
			}
			
			auto c_after = ev.c_after;
			if(c_after != UNSET){                        // Adds individual to population
				for(auto po : sp.comp_gl[c_after].pop_ref_simp){
					PopRef pr; pr.po = po; pr.index = pop_ind[po].size();
					ip.pop_ref.push_back(pr);
					pop_ind[po].push_back(i);
				}
			}
			ip.c = c_after;
			
			egl++;
		}
		
		// Tries to reassign infections in the next time period
		auto t_next = ti+1;
	
		auto egl2 = egl;
		while(egl2 < glob_ev.size() && glob_ev[egl2].tdiv < t_next){
			const auto &ge = glob_ev[egl2];
			auto p = ge.p;
			auto i = ge.i;
			
			const auto &sp = model.species[p];
			auto &ssp = species[p];
			
			auto &ev = ssp.individual[i].ev[ge.e];
			if(ev.type == M_TRANS_EV){
				const auto &tra = sp.tra_gl[ev.tr_gl];
				if(tra.infection.type == TRANS_INFECTION){
					const auto &popnum = popnum_t[ti];  
					
					if(pl){
						cout << endl << endl << "START" << endl;
						cout << "p:" << p << " individual:" << i << "   t:" << model.calc_t(ev.tdiv) <<" t" << endl;
						ssp.print_event(ssp.individual[i].ev);
					}
					
					// Randomly samples a new infection source
				
					auto m = tra.markov_eqn_ref;
					const auto &me = sp.markov_eqn[m];
					const auto &mev = ssp.markov_eqn_vari[m];

					const auto &eq = model.eqn[me.eqn_ref];
			
					const auto &lin = eq.linearise;
					if(!lin.on) emsg("Linearisation should be on");
			
					const auto &pop_ref = eq.pop_ref;
			
					auto Npop = pop_ref.size();
					
					// Works out the probability of the reverse transition
					auto probfi = 0.0;
					auto ii_from = ev.ind_inf_from;
				
					{
						auto k_from = ii_from.pref;
						
						if(lin.multi_source){ 	
							if(k_from == UNSET){
								auto val = eq.calculate_calculation(lin.no_pop_calc,ti,param_val,spline_val);	probfi += log(val/mev.div[ti].value);
							}
							else{
								auto pr = pop_ref[k_from];
								auto val = eq.calculate_calculation(lin.pop_grad_calc[k_from],ti,param_val,spline_val);
							
								probfi += log(val*popnum[pr]/mev.div[ti].value);
							}
						}
						
						if(k_from != UNSET){
							const auto &ind_from = species[ii_from.p].individual[ii_from.i];
				
							auto pr = pop_ref[k_from];
								
							auto num = 1.0;
							const auto &po = model.pop[pr];
						
							if(po.ind_variation){
								for(auto ie : po.ind_eff_mult) num *= ind_from.exp_ie[ie];
								for(auto fe : po.fix_eff_mult) num *= ind_from.exp_fe[fe];
							}
							
							probfi += log(num/popnum[pr]);
						}
					}
						
					// Samples a new event source
					auto probif = 0.0;
					IndInfFrom ii_prop;
					
					{
						auto k_prop = 0u;
					
						if(lin.multi_source){ 	
							auto ss = eq.setup_source_sampler(ti,popnum,param_val,spline_val);

							k_prop = ss.sample_inf_source();
							//if(k_prop == UNSET) emsg("Cannot sample source");
							probif += ss.prob_inf_source(k_prop);
						}
						
						if(k_prop != UNSET){
							if(k_prop == Npop){	
								ii_prop.p = OUTSIDE_INF;
								ii_prop.pref = UNSET;
								ii_prop.po = UNSET;
								ii_prop.i = UNSET;
								ii_prop.w = 1;
							}
							else{
								auto pr = pop_ref[k_prop];
									
								const auto &po = model.pop[pr];
							
								auto p = po.sp_p;
								
								ii_prop.p = p;
								ii_prop.pref = k_prop;
								ii_prop.po = pr;
								
								const auto &list = pop_ind[pr];
								if(po.ind_variation){
									const auto &ssp = species[p];
									
									vector <Poss> poss;
									for(auto i : list){
										const auto &ind  = ssp.individual[i];
										auto num = 1.0;
										for(auto ie : po.ind_eff_mult) num *= ind.exp_ie[ie];
										for(auto fe : po.fix_eff_mult) num *= ind.exp_fe[fe];
										Poss pos; pos.i = i; pos.weight = num;	
										poss.push_back(pos);									
									}
								
									auto pos_sel = sample_possibility(poss);
									ii_prop.i = pos_sel.i;
									ii_prop.w = pos_sel.weight;
									probif += sample_probability(pos_sel.i,poss);
								}
								else{
									if(list.size() == 0) emsg("zero list");
									auto kk = (unsigned int)(ran()*list.size());
									ii_prop.i = list[kk];
									ii_prop.w = 1;
								
									probif += log(1.0/list.size());
								}
							}
						}
					}
					
					auto dlike_markov = probif - probfi;
		
					if(probif != -LARGE && (ii_from.i != ii_prop.i || ii_from.p != ii_prop.p || ii_from.pref != ii_prop.pref)){
						auto n = ev.inf_node_ref;
					
						auto gc = gen_change(REGRAFT_NODE,n,no_alter(n),ii_prop);
				
						auto al = exp(burn_info.markov*dlike_markov 
									+ burn_info.genetic_process*gc.dlike_genetic_process
									+ burn_info.genetic_obs*gc.dlike_genetic_obs
									+ probfi-probif+gc.probfi-gc.probif);
						
						if(pl) cout << al << " al" << endl;
						
						ntr++;
						if(ran() < al){
							if(pl) cout << "accept" << endl;
							nac++;
							gen_change_update(gc);
							
							like.markov += dlike_markov;
							ssp.Li_markov[m][ti] += dlike_markov;
						}
						
						if(pl){
							cout << "AFTER" << endl;
							for(auto &ir : ssp.markov_eqn_vari[m].div[ti].ind_trans){
								cout << ir.i << " ";
								auto ev = ssp.individual[ir.i].ev[ir.index];
								auto iif = ev.ind_inf_from;
								cout << model.calc_t(ev.tdiv) << " p" << iif.p << " i"<< iif.i << " pref" << iif.pref << " po" << iif.po << "  yy" << endl;
							}
						}
						
						if(pl) check(" transtree update2");
					}
				}
			}
		
			egl2++;
		}
		
		if(pl) check_pop_ind(pop_ind,ind_pop);
	}
}


/// Sets inf_ev_from in event at time t based on infection coming from event_from
// If not possible then returns false
bool State::set_ind_inf_from(double t, unsigned int p, vector <Event> &event, unsigned int p_from, unsigned int i_from, const vector <Event> &event_from)
{
	auto e = 0u; while(e < event.size() && event[e].tdiv < t) e++;
	if(e == event.size() || t != event[e].tdiv){
		emsg("Problem");
	}
	
	auto c_from = UNSET;
	auto e2 = 0u; 
	while(e2 < event_from.size() && event_from[e2].tdiv < t){ c_from = event_from[e2].c_after; e2++;}
	
	if(c_from == UNSET) emsg("Should not be unset");
	
	auto &ev = event[e];
	const auto &inf_c = model.inf_cause[p][ev.tr_gl][p_from][c_from];
	auto pref = inf_c.pref;
	//if(pref == UNSET) emsg("Should be p ref");
	if(pref == UNSET) return false;
	
	auto &iif = event[e].ind_inf_from;
	iif.p = p_from;
	iif.i = i_from;
	iif.pref = pref;
	iif.po = inf_c.po;
	iif.w = get_w_from_indinffrom(iif);
	return true;
}


/// Sets inf_ev_from in event at time t based on infection coming from event_from
// If not possible then returns false
void State::set_ind_inf_from_outside(double t, vector <Event> &event)
{
	auto e = 0u; while(e < event.size() && event[e].tdiv < t) e++;
	if(e == event.size() || t != event[e].tdiv) emsg("Problem");
	
	auto &iif = event[e].ind_inf_from;
	iif.p = OUTSIDE_INF;
	iif.i = UNSET;
	iif.pref = UNSET;
	iif.po = UNSET;
	iif.w = 1;
}


/// Performs swapping of infected individual proposal
// Here we look to replace node k_A with its infector node k_B
void State::trans_tree_swap_inf_proposal(const BurnInfo &burn_info, unsigned int &nfa, unsigned int &nac, unsigned int &ntr)
{
	auto pl = false;
	
	auto &inf_node = genetic_value.inf_node;
	if(inf_node.size() == 0) return;
	
	auto &ind_gen_obs = model.genetic_data.ind_gen_obs;
	
	auto loop_max = (unsigned int)(0.3*inf_node.size());
	if(loop_max < 1) loop_max = 1;
	
	for(auto loop = 0u; loop < loop_max; loop++){
		auto k_A = (unsigned int)(inf_node.size()*ran());
		auto &in_A = inf_node[k_A];
		auto &fr_A = in_A.from;
		auto e_index = in_A.from.index;
		auto t_A = in_A.tdiv_start;
		
		if(fr_A.node != ENTER_INF && fr_A.node != OUTSIDE_INF){
			auto k_B = fr_A.node;
			auto &in_B = inf_node[k_B];
			const auto &inf_evB = in_B.inf_ev;
			auto t_B = in_B.tdiv_start;
			
			if(pl) cout << endl << endl << "START " << t_A << " "<< t_B << endl;
			
			// Checks to see if there is an observation
			auto fl = false;
			for(auto e = 0u; e < e_index; e++){
				if(inf_evB[e].type == GENETIC_OBS) fl = true;
			}	
			
			// Checks to see if there is an onward infection which
			// cannot be accounted for deu to finite step size

			auto t_inf_min = round_up(t_A);
			for(auto e = e_index+1; e < inf_evB.size(); e++){
				if(inf_evB[e].tdiv < t_inf_min){
					if(inf_evB[e].type == INFECT_OTHER){
						fl = true;
					}
				}
				else break;
			}
			
			// Checks to see if external infection
			//if(in_B.from.node == ENTER_INF || in_B.from.node == OUTSIDE_INF) fl = true;
			if(in_B.from.node == ENTER_INF) fl = true;
		
			ntr++;
			if(fl == true) nfa++;
			else{				
				// Looks at constructing new timelines
				auto p_A = in_A.p;
				auto &sp_A = species[p_A];
				auto i_A = in_A.i;
				auto ev_A_new = sp_A.individual[i_A].ev;
				
				auto p_B = in_B.p;
				auto &sp_B = species[p_B];
				auto i_B = in_B.i;
				auto ev_B_new = sp_B.individual[i_B].ev;
				
				// Checks if in the time range
				auto fl = false;
				auto e_B_after = get_event_after(in_B.e,ev_B_new);
				if(e_B_after != UNSET && ev_B_new[e_B_after].tdiv < t_A) fl = true;
				 
				auto e_A_before = get_event_before(in_A.e,ev_A_new);
				if(e_A_before != UNSET && ev_A_new[e_A_before].tdiv > t_B) fl = true;
				
				// Check for unaccounted observatins on A
				if(i_A < ind_gen_obs[p_A].size()){
					const auto &igo = ind_gen_obs[p_A][i_A];
					if(igo.size() > 0 && igo[0].tdiv < t_A) fl = true;
				}
			
				if(fl == true) nfa++;
				else{
					if(pl){
						cout << k_A << " " << k_B << " " << species[p_A].individual[i_A].name << " " << species[p_B].individual[i_B].name << "  kkA B" << endl;
						cout << endl << endl << "BEFORE" << endl;
						cout << "Ind" << endl;
						cout << "Event" << endl;
						
						species[p_A].print_event(species[p_A].individual[i_A].ev);
						print_node("node",k_A);
					
						species[p_A].print_event(species[p_B].individual[i_B].ev);
						print_node("node",k_B);
						cout << endl << endl;
					}
					
					auto &ev_B = ev_B_new[in_B.e]; 
				 
					auto fl = false;
					if(ev_B.type == ENTER_EV){
						fl = true;
					}
					else{
						sp_A.move_event(ev_A_new,in_A.e,t_B);
						sp_B.move_event(ev_B_new,in_B.e,t_A);
						
						if(in_B.from.node == OUTSIDE_INF){
							set_ind_inf_from_outside(t_B,ev_A_new);
						}
						else{
							const auto &in_from = inf_node[in_B.from.node];
							auto p_fr = in_from.p, i_fr = in_from.i;
						
							if(set_ind_inf_from(t_B,p_A,ev_A_new,p_fr,i_fr,species[p_fr].individual[i_fr].ev) == false){
								fl = true;
							}
						}
					}
					
					if(fl == true || set_ind_inf_from(t_A,p_B,ev_B_new,p_A,i_A,ev_A_new) == false){
						nfa++;
					}
					else{
						vector <IIFChange> iif_change;
						
						auto ill_fl = false;
						
						// Looks at switch the source for other infections on B's time line
						if(e_index > 0){
							const auto &ev_bef = species[p_B].individual[i_B].ev;
							const auto &ev_aft = ev_A_new;
							auto c_bef = UNSET;
							auto c_aft = UNSET;
							auto e_bef = 0u, e_aft = 0u;;
							
							for(auto e = 0u; e < e_index; e++){
								const auto &ie = in_B.inf_ev[e];
								auto tt = ie.tdiv;
								
								while(e_bef < ev_bef.size() && ev_bef[e_bef].tdiv < tt){
									c_bef = ev_bef[e_bef].c_after;
									e_bef++;
								}
								
								while(e_aft < ev_aft.size() && ev_aft[e_aft].tdiv < tt){
									c_aft = ev_aft[e_aft].c_after;
									e_aft++;
								}
								
								if(p_B != p_A || c_bef != c_aft){  // The source has switched
									const auto &in_inf = inf_node[ie.index];
								
									auto p_ev = in_inf.p, i_ev = in_inf.i, e_ev = in_inf.e; 
									const auto &ev = species[p_ev].individual[i_ev].ev[e_ev];
									
									auto ti = get_ti(tt);
									
									IIFChange iif_ch;                // Used to store a change in source
									iif_ch.p_ev = p_ev; 
									iif_ch.i_ev = i_ev; 
									iif_ch.e_ev = e_ev; 
									iif_ch.ti = ti;
									
									// Looks at before
									const auto &iif_bef = ev.ind_inf_from;
									const auto &inf_c_bef = model.inf_cause[p_ev][ev.tr_gl][p_B][c_bef];
									if(p_B != iif_bef.p) emsg("p not agree");
									if(i_B != iif_bef.i) emsg("p not agree");
									if(inf_c_bef.pref != iif_bef.pref) emsg("pref not agree");
									if(inf_c_bef.po != iif_bef.po) emsg("po not agree");
									
									const auto &eq_bef = model.eqn[inf_c_bef.eq_ref];
									const auto &lin_bef = eq_bef.linearise;
									auto va_bef = eq_bef.calculate_calculation(lin_bef.pop_grad_calc[iif_bef.pref],ti,param_val,spline_val);
									
									iif_ch.Li_markov_bef = log(va_bef*iif_bef.w);
									
									// Looks at after
									const auto &inf_c_aft = model.inf_cause[p_ev][ev.tr_gl][p_A][c_aft];
									auto pref = inf_c_aft.pref;
									if(pref == UNSET) ill_fl = true;
									else{
										//if(pref == UNSET) emsg("Should be p ref");
								
										const auto &eq_aft = model.eqn[inf_c_aft.eq_ref];
										const auto &lin_aft = eq_aft.linearise;
										auto va_aft = eq_aft.calculate_calculation(lin_aft.pop_grad_calc[pref],ti,param_val,spline_val);
		
										auto &iif_aft = iif_ch.ind_inf_from;
										iif_aft.p = p_A;
										iif_aft.i = i_A;
										iif_aft.pref = pref;
										iif_aft.po = inf_c_aft.po;
										iif_aft.w = get_w_from_indinffrom(iif_aft);
										
										iif_ch.Li_markov_aft = log(va_aft*iif_aft.w);
								
										iif_change.push_back(iif_ch);							
									}
								}
							}
						}
						
						if(ill_fl == false){
							back_init();
							auto like_ch = update_ind(p_A,i_A,ev_A_new,UP_MULTI);
							
							auto li_ch = update_ind(p_B,i_B,ev_B_new,UP_MULTI);
							
							add_on_like(li_ch,like_ch);
								
							for(const auto &iif_ch : iif_change){
								like_ch.markov += iif_ch.Li_markov_aft - iif_ch.Li_markov_bef;
							}
							
							auto al = calc_al(like_ch,0,burn_info);
							
							if(pl) print_like(li_ch);
							
							if(pl) cout << al << " al" << endl;
							if(ran() < al){
								if(pl) cout << "accept" << endl;
								add_like(like_ch);
							
								// Updates quantities for inf_node
								InfNode in_A_new;
								in_A_new.tdiv_start = t_B;
								in_A_new.tdiv_rec = in_A.tdiv_rec;
								in_A_new.p = p_A;
								in_A_new.i = i_A;
								in_A_new.e = in_A.e;
								in_A_new.from = in_B.from;
								for(auto e = 0u; e <= e_index; e++){
									in_A_new.inf_ev.push_back(in_B.inf_ev[e]);
								}
								
								in_A_new.inf_ev[e_index].index = fr_A.node;

								for(auto e = 0u; e < in_A.inf_ev.size(); e++){
									in_A_new.inf_ev.push_back(in_A.inf_ev[e]);
								}
								
								// Updates quantities for inf_node from
								InfNode in_B_new;
								in_B_new.tdiv_start = t_A;
								in_B_new.tdiv_rec = in_B.tdiv_rec;
								in_B_new.p = p_B;
								in_B_new.i = i_B;
								in_B_new.e = in_B.e;
								in_B_new.from.node = k_B;
								in_B_new.from.index = e_index;
								for(auto e = e_index+1; e < in_B.inf_ev.size(); e++){
									in_B_new.inf_ev.push_back(in_B.inf_ev[e]);
								}
								
								in_A = in_A_new;
								in_B = in_B_new;
		
								for(const auto &iif_ch : iif_change){
									auto p_ev = iif_ch.p_ev;
									auto &ssp = species[p_ev];
									auto &ev = ssp.individual[iif_ch.i_ev].ev[iif_ch.e_ev];
									ev.ind_inf_from = iif_ch.ind_inf_from;
									
									auto m = model.species[p_ev].tra_gl[ev.tr_gl].markov_eqn_ref;
									ssp.Li_markov[m][iif_ch.ti] += iif_ch.Li_markov_aft - iif_ch.Li_markov_bef;
								}
								
								update_inf_node_ref(k_A);
								update_inf_node_ref(k_B);
							
								if(pl){
									cout << "AFTER" << endl;
									cout << "Ind" << endl;
									cout << "Event" << endl;
									species[0].print_event(species[p_A].individual[i_A].ev);
									print_node("node",k_A);
								
									species[0].print_event(species[p_B].individual[i_B].ev);
									print_node("node",k_B);
									cout << endl << endl;
								}
								
								update_popnum_ind(p_A,i_A);
								update_popnum_ind(p_B,i_B);
								
								nac++;
							}
							else{
								if(pl) cout << "reject" << endl;
								restore_back();
							}
						}
					}
				}
			}
			if(pl) check(" transtree update2");
		}
	}
}


/// Updates reference from inf_node
void State::update_inf_node_ref(unsigned int k)
{
	auto &inf_node = genetic_value.inf_node;
	auto &obs_node_ref = genetic_value.obs_node_ref;
	
	const auto &in = inf_node[k];
	auto &from = in.from;
	
	if(from.node == ENTER_INF || from.node == OUTSIDE_INF){
		genetic_value.inf_origin[from.index].node = k;
	}
	else{
		auto &in_from = inf_node[from.node];
		in_from.inf_ev[from.index].index = k;
		
		auto &ev = species[in.p].individual[in.i].ev[in.e];
		ev.inf_node_ref = k;
	}
	
	for(auto e = 0u; e < in.inf_ev.size(); e++){
		const auto &ev = in.inf_ev[e];
	
		switch(ev.type){	
		case INFECT_OTHER:
			{
				auto &in2 = inf_node[ev.index];
				in2.from.node = k;
				in2.from.index = e;
				auto &ev = species[in2.p].individual[in2.i].ev[in2.e];
				auto &iif = ev.ind_inf_from;
				if(!(iif.p == in.p && iif.i == in.i)){ 
					iif.p = in.p;
					iif.i = in.i;
				}
			}
			break;
			
		case GENETIC_OBS:
			obs_node_ref[ev.index].index = e;
			break;
		}
	}
}


/// Gets an inf_node n without alteration
InfNodeAlter State::no_alter(unsigned int n) const
{
	auto na = InfNodeAlter();
	const auto &in = genetic_value.inf_node[n];
	na.tdiv_start = in.tdiv_start;
	na.tdiv_rec = in.tdiv_rec;
	return na;
}
	
	
/// Performs a transmission tree mutation proposal
void State::trans_tree_mut_proposal(const BurnInfo &burn_info, unsigned int &nac, unsigned int &ntr, double &si)
{
	auto pl = false;
	
	auto &inf_node = genetic_value.inf_node;
	const auto &gen_data = model.genetic_data;
	
	auto mut_rate = genetic_value.mut_rate*model.details.dt;
	auto seq_var = genetic_value.seq_var;

	auto &gen_dif = genetic_value.gen_dif;
	const auto &gen_val = gen_data.gen_dif;
	
	// Goes through the infection nodes
	for(auto n = 0u; n < inf_node.size(); n++){
		auto &in = inf_node[n];
		
		for(auto e = 0u; e < in.inf_ev.size(); e++){
			auto &iev = in.inf_ev[e];
	
			string warn;
			int dn = normal_int_sample(si,warn); 
			if(dn == UNSET) emsg("MCMC proposal problem: "+warn);
		
			if(dn != 0){
				int num = iev.mut_num;
				auto num_new = num + dn;
				
				if(num_new >= 0){
					double dtdiv;
					if(e == 0) dtdiv = iev.tdiv-in.tdiv_start;
					else dtdiv = iev.tdiv-in.inf_ev[e-1].tdiv;
					auto lam = mut_rate*dtdiv;
					
					auto dlike_genetic_process = poisson_probability(num_new,lam) - poisson_probability(num,lam);
			
					auto obs_node = get_gen_obs_list(n,e);
					
					auto nobs = genetic_value.obs_node_ref.size();
					vector <bool> obs_done(nobs,false);	
					for(auto m : obs_node) obs_done[m] = true;

					vector <unsigned int> obs_other;
					for(auto m = 0u; m < nobs; m++){
						if(obs_done[m] == false && genetic_value.obs_node_ref[m].node != UNSET) obs_other.push_back(m);
					}
					
					auto sum = 0.0;
					double d, d_update;
					for(auto m1 : obs_node){
						const auto &gen_dif2 = gen_dif[m1];
						const auto &gen_val2 = gen_val[m1];
						for(auto m2 : obs_other){
							d = int(gen_dif2[m2]) - int(gen_val2[m2]);
							d_update = d + dn;
							sum += d_update*d_update - d*d;
						}
					}
					auto dlike_genetic_obs = -2*sum*model.genetic_data.obs_scale;
					
					auto al = exp(burn_info.genetic_process*dlike_genetic_process
									+ burn_info.genetic_obs*dlike_genetic_obs);
			
					auto fac = -0.005;
			
					ntr++;
					if(ran() < al){ 
						fac = 0.01;	
						nac++;
					
						like.genetic_process += dlike_genetic_process;
						like.genetic_obs += dlike_genetic_obs;
					
						for(auto m1 : obs_node){
							for(auto m2 : obs_other){
								gen_dif[m1][m2] += dn;
								gen_dif[m2][m1] += dn;
							}
						}
						
						iev.mut_num = num_new;	
					}
						
					if(burn_info.on){ si *= 1+fac*burn_info.fac; if(si < 1) si = 1;}
			
					if(pl) check("update");
				}
			}
		}
	}
	
	// Goes through the infection origins
	for(auto &io : genetic_value.inf_origin){ 
		const auto &in = inf_node[io.node];
		
		string warn;
		int dn = normal_int_sample(si,warn);
		if(dn == UNSET) emsg("MCMC proposal problem: "+warn);
		
		if(dn != 0){
			int num = io.mut_num;
			auto num_new = num + dn;
			
			if(num_new >= 0){
				auto lam = seq_var+mut_rate*in.tdiv_start;
			
				auto dlike_genetic_process = poisson_probability(num_new,lam) - poisson_probability(num,lam);
		
				auto obs_node = get_gen_obs_list(io.node,0);
				
				auto nobs = genetic_value.obs_node_ref.size();
				vector <bool> obs_done(nobs,false);	
				for(auto m : obs_node) obs_done[m] = true;

				vector <unsigned int> obs_other;
				for(auto m = 0u; m < nobs; m++){
					if(obs_done[m] == false && genetic_value.obs_node_ref[m].node != UNSET) obs_other.push_back(m);
				}
				
				auto sum = 0.0;
				double d, d_update;
				for(auto m1 : obs_node){
					const auto &gen_dif2 = gen_dif[m1];
					const auto &gen_val2 = gen_val[m1];
					for(auto m2 : obs_other){
						d = int(gen_dif2[m2]) - int(gen_val2[m2]);
						d_update = d + dn;
						sum += d_update*d_update - d*d;
					}
				}
				auto dlike_genetic_obs = -2*sum*model.genetic_data.obs_scale;
				
				auto al = exp(burn_info.genetic_process*dlike_genetic_process
								+ burn_info.genetic_obs*dlike_genetic_obs);
	
				auto fac = -0.005;
				
				ntr++;
				if(ran() < al){ 
					fac = 0.01;	
					nac++;
				
					like.genetic_process += dlike_genetic_process;
					like.genetic_obs += dlike_genetic_obs;
				
					for(auto m1 : obs_node){
						for(auto m2 : obs_other){
							gen_dif[m1][m2] += dn;
							gen_dif[m2][m1] += dn;
						}
					}
					
					io.mut_num = num_new;	
				}
					
				if(burn_info.on){ si *= 1+fac*burn_info.fac; if(si < 1) si = 1;}
		
				if(pl) check("update");
			}
		}
	}
}


/// Performs local transmission tree mutation proposals (the interchange adjacent mutation events)
void State::trans_tree_mut_local_proposal(const BurnInfo &burn_info, vector <GeneticProp> &gen_mut_info)
{
	auto pl = false;// true;
	
	auto &inf_node = genetic_value.inf_node;
	const auto &gen_data = model.genetic_data;
	
	auto mut_rate = genetic_value.mut_rate*model.details.dt;
	auto seq_var = genetic_value.seq_var;
	
	auto nobs = genetic_value.obs_node_ref.size();
	
	for(auto loop = 0u; loop < GEN_UPDATE_MAX; loop++){
		if(pl) cout << loop << " loop" << endl;
		auto type = OBS_GEN_UPDATE;
		switch(loop){
		case 0: type = OBS_GEN_UPDATE; break;
		case 1: type = SINGLE_BRANCH; break;
		case 2: type = TWO_BRANCH1; break;
		case 3: type = TWO_BRANCH2; break;
		case 4: type = TWO_BRANCH3; break;
		case 5: type = ORIGIN_GEN_UPDATE; break;
		default: emsg("Not op"); break;
		}
	
		auto &gmi = gen_mut_info[type];

		if(type != ORIGIN_GEN_UPDATE){
			for(auto n = 0u; n < inf_node.size(); n++){
				auto &in = inf_node[n];

				for(auto e = 0u; e < in.inf_ev.size(); e++){
					string warn;
					int dn = normal_int_sample(gmi.si,warn);
					if(dn == UNSET) emsg("MCMC proposal problem: "+warn);
		
					if(dn != 0){
						auto &ev = in.inf_ev[e];
					
						if(type == OBS_GEN_UPDATE){
							if(ev.type == GENETIC_OBS){
								if(e+1 < in.inf_ev.size()){
									auto &ev2 = in.inf_ev[e+1];
				
									int num1_old = ev.mut_num;
									int num2_old = ev2.mut_num;
									int num1_new = num1_old+dn;
									int num2_new = num2_old-dn;
										
									if(num1_new >= 0 && num2_new  >= 0){
										auto msel = ev.index;
										
										// Finds change in genetic process
										double t_last; 
										if(e > 0) t_last = in.inf_ev[e-1].tdiv; 
										else t_last = in.tdiv_start;
									
										auto lam1 = mut_rate*(ev.tdiv-t_last);
										auto lam2 = mut_rate*(ev2.tdiv-ev.tdiv);
							
										auto dlike_genetic_process = poisson_probability(num1_new,lam1) - poisson_probability(num1_old,lam1)
																								+poisson_probability(num2_new,lam2) - poisson_probability(num2_old,lam2);
									
										// Finds change in gen_dif and genetic_obs
										auto obs_node = get_gen_obs_list(n,e+1);

										vector <bool> obs_done(nobs,false);	
										for(auto m : obs_node) obs_done[m] = true;
										obs_done[msel] = true; 

										vector <unsigned int> obs_other;
										for(auto m = 0u; m < nobs; m++){
											if(obs_done[m] == false && genetic_value.obs_node_ref[m].node != UNSET) obs_other.push_back(m);
										}

										auto &gen_dif = genetic_value.gen_dif;
										
										const auto &gen_val2 = gen_dif[msel];
										const auto &gen_dif2 = gen_data.gen_dif[msel];
										
										auto sum = 0.0;
										double d, d_update;
										for(auto m : obs_node){
											d = int(gen_dif2[m]) - int(gen_val2[m]);
											d_update = d + dn;
											sum += d_update*d_update - d*d;
										}
										
										for(auto m : obs_other){
											d = int(gen_dif2[m]) - int(gen_val2[m]);
											d_update = d - dn;
											sum += d_update*d_update - d*d;
										}
								
										auto dlike_genetic_obs = -2*sum*model.genetic_data.obs_scale;
							
										auto al = exp(burn_info.genetic_process*dlike_genetic_process
															+ burn_info.genetic_obs*dlike_genetic_obs);
				
										if(pl) cout << al << " al" << endl;
										auto fac = -0.005;
									
										gmi.ntr++;
										if(ran() < al){ 
											fac = 0.01;	
											gmi.nac++;
										
											like.genetic_process += dlike_genetic_process;
											like.genetic_obs += dlike_genetic_obs;
										
											for(auto m : obs_node){
												gen_dif[msel][m] -= dn;
												gen_dif[m][msel] -= dn;
											}
											
											for(auto m : obs_other){
												gen_dif[msel][m] += dn;
												gen_dif[m][msel] += dn;
											}
											
											ev.mut_num = num1_new;	
											ev2.mut_num = num2_new;	
										}
											
										if(burn_info.on){ gmi.si *= 1+fac*burn_info.fac; if(gmi.si < 1) gmi.si = 1;}
										
										if(pl) check("update");
									}
								}
							}
						}
						else{
							if(ev.type == INFECT_OTHER){
								auto &in2 = inf_node[ev.index];
							
								auto nbranch = 0u;
								if(e+1 < in.inf_ev.size()) nbranch++;
								if(in2.inf_ev.size() > 0) nbranch++;
								
								if((nbranch == 1 && type == SINGLE_BRANCH) ||
						 	 		 (nbranch == 2 && type == TWO_BRANCH1) ||
								 	 (nbranch == 2 && type == TWO_BRANCH2) ||
									 (nbranch == 2 && type == TWO_BRANCH3)){
								
									auto dn2 = -dn, dn3 = -dn;
									switch(type){
									case TWO_BRANCH1: dn2 = -dn, dn3 = -dn; break;
									case TWO_BRANCH2: dn2 = dn, dn3 = -dn; break;
									case TWO_BRANCH3: dn2 = -dn, dn3 = dn; break;
									default: break;
									}
								
									auto ill = false;
									
									// Works out change in likelihood for genetic process
									auto dlike_genetic_process = 0.0;
									
									{
										auto num_old = ev.mut_num;
										auto num_new = int(num_old)+dn;
										if(num_new < 0) ill = true;
										else{
											double t_last; 
											if(e > 0) t_last = in.inf_ev[e-1].tdiv;
											else t_last = in.tdiv_start;
										
											auto lam2 = mut_rate*(ev.tdiv-t_last);
											dlike_genetic_process += poisson_probability(num_new,lam2) - poisson_probability(num_old,lam2);
										}
									}
									
									if(ill == false && e+1 < in.inf_ev.size()){
										auto &ev2 = in.inf_ev[e+1];
										auto num_old = ev2.mut_num;
										auto num_new = int(num_old)+dn2;
										if(num_new < 0) ill = true;
										else{
											auto lam2 = mut_rate*(ev2.tdiv-ev.tdiv);
											dlike_genetic_process += poisson_probability(num_new,lam2) - poisson_probability(num_old,lam2);
										}
									}
									
									if(ill == false && in2.inf_ev.size() > 0){
										auto &ev2 = in2.inf_ev[0];
										auto num_old = ev2.mut_num;
										auto num_new = int(num_old)+dn3;
										if(num_new < 0) ill = true;
										else{
											auto lam2 = mut_rate*(ev2.tdiv-ev.tdiv);
											dlike_genetic_process += poisson_probability(num_new,lam2) - poisson_probability(num_old,lam2);
										}
									}
									
									if(ill == false){	
										// Works out change is genetic_obs (if necessary)
										auto dlike_genetic_obs = 0.0;
										
										vector <unsigned int> obs_node1, obs_node2;
										auto &gen_dif = genetic_value.gen_dif;
										
										auto dif = 2*dn;
																				
										if(nbranch == 2){
											// Finds change in gen_dif and genetic_obs
											auto sum = 0.0;
											double d, d_update;
											
											switch(type){
											case TWO_BRANCH1:
												obs_node1 = get_gen_obs_list(n,e+1);
												obs_node2 = get_gen_obs_list(ev.index,0);
												break;
												
											case TWO_BRANCH2:
												{
													obs_node1 = get_gen_obs_list(n,e+1);
													auto obs_node_not = get_gen_obs_list(ev.index,0);
													
													auto N = gen_data.obs.size();
													vector <bool> map(N,false);
													for(auto m : obs_node_not) map[m] = true;
													for(auto m : obs_node1) map[m] = true;
													for(auto m = 0u; m < N; m++){
														if(map[m] == false && genetic_value.obs_node_ref[m].node != UNSET){
															obs_node2.push_back(m);
														}
													}
													dif *= -1;
												}
												break;
												
											case TWO_BRANCH3:
												{
													auto obs_node_not = get_gen_obs_list(n,e+1);
													obs_node2 = get_gen_obs_list(ev.index,0);
													auto N = gen_data.obs.size();
													vector <bool> map(N,false);
													for(auto m : obs_node_not) map[m] = true;
													for(auto m : obs_node2) map[m] = true;
													for(auto m = 0u; m < N; m++){
														if(map[m] == false && genetic_value.obs_node_ref[m].node != UNSET){
															obs_node1.push_back(m);
														}
													}
													dif *= -1;
												}
												break;
											default: emsg("not def"); break;
											}
											
											for(auto m1 : obs_node1){
												const auto &gen_val2 = gen_dif[m1];
												const auto &gen_dif2 = gen_data.gen_dif[m1];
											
												for(auto m2 : obs_node2){	
													d = int(gen_dif2[m2]) - int(gen_val2[m2]);
													d_update = d + dif;
													sum += d_update*d_update - d*d;
												}
											}
											dlike_genetic_obs += -2*sum*model.genetic_data.obs_scale;
										}
										
										auto al = exp(burn_info.genetic_process*dlike_genetic_process
																+ burn_info.genetic_obs*dlike_genetic_obs);
										
										if(pl) cout << al << " " << dlike_genetic_process << " " << dlike_genetic_obs << " al" << endl;
										auto fac = -0.005;
									
										gmi.ntr++;
										if(ran() < al){
											fac = 0.01;	
											gmi.nac++;
										
											like.genetic_process += dlike_genetic_process;
											like.genetic_obs += dlike_genetic_obs;
										
											ev.mut_num += dn;	
											if(e+1 < in.inf_ev.size()) in.inf_ev[e+1].mut_num += dn2;
											if(in2.inf_ev.size() > 0) in2.inf_ev[0].mut_num += dn3;
											
											if(nbranch == 2){
												for(auto m1 : obs_node1){
													for(auto m2 : obs_node2){	
														gen_dif[m1][m2] -= dif;
														gen_dif[m2][m1] -= dif;
													}
												}
											}
										}
											
										if(burn_info.on){
											gmi.si *= 1+fac*burn_info.fac; if(gmi.si < 1) gmi.si = 1;
										}
										
										if(pl) check("update21");
									}
								}
							}
							break;
						}
					}
				}
			}
		}
		else{
			// Interchanges the value at the origin with the mutations in the first event
			for(auto &io : genetic_value.inf_origin){
				auto &in = inf_node[io.node];
				
				string warn;
				int dn = normal_int_sample(gmi.si,warn);
				if(dn == UNSET) emsg("MCMC tree proposal problem: "+warn);
		
				if(dn != 0 && in.inf_ev.size() > 0){
					auto &ev = in.inf_ev[0];
					
					int num_old = io.mut_num;
					int num_new = num_old + dn;
					int num2_old = ev.mut_num;
					int num2_new = num2_old - dn;
					
					if(num_new >= 0 && num2_new >= 0){
						auto lam = seq_var+mut_rate*in.tdiv_start;
					
						auto dlike_genetic_process = poisson_probability(num_new,lam) - poisson_probability(num_old,lam);
				
						auto lam2 = mut_rate*(ev.tdiv-in.tdiv_start);
						dlike_genetic_process += poisson_probability(num2_new,lam2) - poisson_probability(num2_old,lam2);
					
						auto al = exp(burn_info.genetic_process*dlike_genetic_process);
						auto fac = -0.005;
						
						gmi.ntr++;
						if(ran() < al){ 
							fac = 0.01;	
							gmi.nac++;
						
							like.genetic_process += dlike_genetic_process;
						
							io.mut_num = num_new;	
							ev.mut_num = num2_new;
						}
							
						if(burn_info.on){ gmi.si *= 1+fac*burn_info.fac; if(gmi.si < 1) gmi.si = 1;}
				
						if(pl) check("update31");
					}
				}
			}
		}
	}
}


/// Prints information about an infection node
void State::print_node(string te, unsigned int n) const 
{
	if(n == OUTSIDE_INF){ cout << te << ": OUTSIDE_INF" << endl; return;}
	if(n == ENTER_INF){ cout << te << ": ENTER_INF" << endl; return;}
	
	const auto &in = genetic_value.inf_node[n];
		
	cout << te << ": Infection node " << n << endl;
	cout << "t_start=" << model.calc_t(in.tdiv_start);
	cout << "  t_rec=" << model.calc_t(in.tdiv_rec);
	cout << "  p=" << in.p << "  i=" << in.i << "  e="<< in.e << endl;
	
	if(in.from.node == OUTSIDE_INF){
		cout << "from outside " <<  in.from.index << " ";
		cout << "mut_num = " << genetic_value.inf_origin[in.from.index].mut_num;
	}
	else{
		if(in.from.node == ENTER_INF){
			cout << "from enter " <<  in.from.index << " ";
			cout << "mut_num = " << genetic_value.inf_origin[in.from.index].mut_num;
		}
		else cout << "from: node=" << in.from.node << "   index=" << in.from.index;
	}
	cout << endl;
	
	cout << "inf_ev" << endl;
	for(const auto &ie : in.inf_ev){
		cout << "t=" << model.calc_t(ie.tdiv) << "  ";
		switch(ie.type){
		case INFECT_OTHER: cout << "INFECT_OTHER (" << ie.index << ")"; break;
		case GENETIC_OBS: cout << "GENETIC OBS (" << ie.index << ")"; break;
		}
		
		cout << "  mut_num=" << ie.mut_num << endl;
	}
	cout << endl;
}


/// Prints ionformation about an infection node
void State::print_inf_origin() const 
{
	for(const auto &io : genetic_value.inf_origin){
		cout << io.node << " " << io.mut_num << " Inf origin" << endl;
	}
}


/// Rounds down to the begining of the division
double State::round_down(double t) const 
{
	return (unsigned int)(t);
}


/// Rounds up to the end of the division
double State::round_up(double t) const 
{
	return 1+(unsigned int)(t);
}


/// In the infection tree removes line from node to node which infected it
void State::disconnect_inf_node(unsigned int n)
{
	auto &inf_node = genetic_value.inf_node;
	
	auto &in = inf_node[n];
	
	auto &fr = in.from;
	if(fr.node == ENTER_INF || fr.node == OUTSIDE_INF){
		auto &inf_origin = genetic_value.inf_origin;
		
		auto ind = fr.index;
		auto si = inf_origin.size();
		if(ind >= si) emsg("inf_origin out of range");
		if(ind+1 < si){
			inf_origin[ind] = inf_origin[si-1];
			if(inf_node[inf_origin[ind].node].from.index != si-1) emsg("prob6");
			inf_node[inf_origin[ind].node].from.index = ind;
		}			
		inf_origin.pop_back();
	}
	else{	
		auto ind = fr.index;
		auto &in2 = inf_node[fr.node];
		auto &inf_ev2 = in2.inf_ev;
	
		double t_last; if(ind > 0) t_last = inf_ev2[ind-1].tdiv; else t_last = in2.tdiv_start;
		auto t1 = inf_ev2[ind].tdiv;
		
		if(ind+1 < inf_ev2.size()){                    // Combines together mutations
			auto t2 = inf_ev2[ind+1].tdiv;
			
			auto frac = (t1-t_last)/(t2-t_last);
			if(frac < 0 || frac > 1) emsg("frac out of range");
		
			if(model.genetic_data.on) inf_ev2[ind+1].mut_num += inf_ev2[ind].mut_num;
			inf_ev2.erase(inf_ev2.begin()+ind);
		}
		else{                                          // Removes last entry      
			inf_ev2.pop_back();
		}
		
		for(auto i = ind; i < inf_ev2.size(); i++){
			const auto &ie2 = inf_ev2[i];
			switch(ie2.type){
			case INFECT_OTHER:
				{
					auto &in3 = inf_node[ie2.index];
					if(in3.from.index != i+1) emsg("prob7");
					in3.from.index--;
				}
				break;
				
			case GENETIC_OBS:
				{
					auto &onr = genetic_value.obs_node_ref[ie2.index];
					if(onr.index != i+1) emsg("prob8");
					onr.index--;
				}
				break;
			}
		}
	}
	
	fr.node = UNSET;
	fr.index = UNSET;
	
	auto &ev_inf = species[in.p].individual[in.i].ev[in.e];
	auto &iif = ev_inf.ind_inf_from;
	iif.p = UNSET;
	iif.i = UNSET;
	iif.pref = UNSET;
}


/// In the infection tree adds line from node to node which infected it
void State::connect_inf_node(unsigned int n, const IndInfFrom &iif, unsigned int mut_num_store)
{
	if(model.genetic_data.on){
		if(mut_num_store == UNSET) emsg("Should not be unset");
	}
	else{
		if(mut_num_store != UNSET) emsg("Should be unset "+to_string(mut_num_store));
	}
	
	auto &inf_node = genetic_value.inf_node;
	
	auto &in = inf_node[n];
	auto &fr = in.from;

	if(iif.p == ENTER_INF || iif.p == OUTSIDE_INF){
		auto &inf_origin = genetic_value.inf_origin;
	
		fr.node = iif.p;
		fr.index = inf_origin.size();
		
		InfOrigin io;
		io.node = n;		
		io.mut_num = mut_num_store;
		inf_origin.push_back(io);
	}
	else{
		auto t = in.tdiv_start;
		
		// Works out the node which instigates the infection
		
		auto nr = get_ind_noderef(iif.p,iif.i,t);
		auto &in2 = genetic_value.inf_node[nr.node];
		
		auto &inf_ev2 = in2.inf_ev;
		auto j = nr.index;
		
		double t_last;
		
		if(j > 0) t_last = inf_ev2[j-1].tdiv; else t_last = in2.tdiv_start;
		auto t1 = t;

		fr = nr;
			
		InfEvent ie_add;
		ie_add.tdiv = t;
		ie_add.type = INFECT_OTHER;
		ie_add.index = n;
		
		if(j < inf_ev2.size()){                        // Inserting between two events
			auto t2 = inf_ev2[j].tdiv;
			
			auto frac = (t1-t_last)/(t2-t_last);
			if(frac < 0 || frac > 1) emsg("frac out of range");
			
			if(model.genetic_data.on){
				auto mut_num = inf_ev2[j].mut_num;
				inf_ev2[j].mut_num = mut_num-mut_num_store;
				ie_add.mut_num = mut_num_store;
			}
			else{
				ie_add.mut_num = UNSET;
			}
			
			inf_ev2.insert(inf_ev2.begin()+j,ie_add);
		}
		else{                                          // Adds event to end
			ie_add.mut_num = mut_num_store;
			inf_ev2.push_back(ie_add);
		}
		

		for(auto i = j+1; i < inf_ev2.size(); i++){
			const auto &ie2 = inf_ev2[i];
			switch(ie2.type){
			case INFECT_OTHER:
				{
					auto &in3 = inf_node[ie2.index];
					if(in3.from.index != i-1) emsg("prob9");
					in3.from.index++;
				}
				break;
				
			case GENETIC_OBS:
				{
					auto &onr = genetic_value.obs_node_ref[ie2.index];
					if(onr.index != i-1) emsg("prob10");
					onr.index++;
				}
				break;
			}
		}
	}

	species[in.p].individual[in.i].ev[in.e].ind_inf_from = iif;
}


/// Finds the node reference for the individual and time
NodeRef State::get_ind_noderef(unsigned int p, unsigned int i, double t) const 
{
	if(p == OUTSIDE_INF || p == ENTER_INF){
		NodeRef nr; nr.node = p; nr.index = UNSET;
		return nr;
	}
	
	auto &ev = species[p].individual[i].ev;
	if(ev.size() == 0) emsg("problem ev");
		
	auto t_inf = round_down(t);

	int e = ev.size()-1;
	while(e >= 0 && ev[e].tdiv > t_inf) e--;
	if(e < 0) emsg("Problem e");
		
	while(e >= 0 && ev[e].inf_node_ref == UNSET) e--;
	if(e < 0){
		emsg("Problem e2");
	}

	auto n2 = ev[e].inf_node_ref;
		
	auto &inf_ev2 = genetic_value.inf_node[n2].inf_ev;
		
	auto j = 0u; while(j < inf_ev2.size() && inf_ev2[j].tdiv < t) j++;
	
	NodeRef nr; nr.node = n2; nr.index = j;
	
	return nr;
}


/// Checks that pop_ind is correctly specified
void State::check_pop_ind(vector < vector <unsigned int> > pop_ind, const vector < vector <IndPop> > &ind_pop) const 
{
	vector < vector <unsigned int> > pop_ind_compare;
	pop_ind_compare.resize(model.pop.size());
	
	for(auto p = 0u; p < model.species.size(); p++){
		const auto &sp = model.species[p];
		const auto &ssp = species[p];
		for(auto i = 0u; i < ssp.individual.size(); i++){	
			auto c = ind_pop[p][i].c;
			
			if(c != UNSET){
				for(auto po : sp.comp_gl[c].pop_ref_simp){
					pop_ind_compare[po].push_back(i);
				}
			}
		}
	}
	
	for(auto po = 0u; po < model.pop.size(); po++){
		sort(pop_ind[po].begin(),pop_ind[po].end());
		sort(pop_ind_compare[po].begin(),pop_ind_compare[po].end());
		
		if(!equal_vec(pop_ind[po],pop_ind_compare[po])) emsg("pop_ind problem");
	}
}


/// Works out how genetics change as a result of a node being regraphed / added removed
GenChange State::gen_change(GenChaType type, unsigned int n, const InfNodeAlter &node_alter, const IndInfFrom &iif_add) const 
{			
	auto data_on = model.genetic_data.on;
		
	auto mut_rate = genetic_value.mut_rate*model.details.dt;
	auto seq_var = genetic_value.seq_var;
	
	GenChange gc(type);
	gc.n = n;
	gc.iif_add = iif_add;
	gc.node_alter = node_alter;
	gc.regraft_type = NORMAL_REGRAFT;
	
	auto t_add = node_alter.tdiv_start;
	
	switch(type){
	case GENCHA_FAIL: emsg("Should not fail"); break;
	
	case REGRAFT_NODE:
		{
			const auto &inf_node = genetic_value.inf_node;
			auto nr_add = get_ind_noderef(iif_add.p,iif_add.i,t_add); 
			gc.nr_add = nr_add;
		
			const auto &in = inf_node[n];
	
			auto t_now = in.tdiv_start;
			
			const auto &nr_from = in.from;

			if((nr_from.node == ENTER_INF || nr_from.node == OUTSIDE_INF)
				&& nr_from.node > CODE && nr_add.node > CODE){
				nr_add.index = nr_from.index;
				gc.regraft_type = ORIGIN_REGRAFT;
		
				if(t_add != t_now && data_on){
					auto mut_num = genetic_value.inf_origin[nr_from.index].mut_num;
					
					gc.dlike_genetic_process += poisson_probability(mut_num,seq_var+mut_rate*t_add)
								-poisson_probability(mut_num,seq_var+mut_rate*t_now);
				}
			}					
			else{	
				if(nr_from.node == nr_add.node && (nr_from.index == nr_add.index || nr_from.index+1 == nr_add.index)){ 
					if(data_on){
						const auto &in_from = inf_node[nr_from.node];
						const auto &inf_ev_from = in_from.inf_ev;
							
						// Keeps all the mutation numbers the same
						
						auto j_from = nr_from.index;
						
						double t1; 
						if(j_from > 0) t1 = inf_ev_from[j_from-1].tdiv; 
						else t1 = in_from.tdiv_start;
				
						auto numA = inf_ev_from[j_from].mut_num;
				
						gc.dlike_genetic_process += poisson_probability(numA,mut_rate*(t_add-t1)) 
																			- poisson_probability(numA,mut_rate*(t_now-t1)); 
						
						if(j_from+1 < inf_ev_from.size()){
							auto numB = inf_ev_from[j_from+1].mut_num;
							auto t2 = inf_ev_from[j_from+1].tdiv;
							gc.dlike_genetic_process += poisson_probability(numB,mut_rate*(t2-t_add)) 
																			- poisson_probability(numB,mut_rate*(t2-t_now)); 
						}
					}
					
					gc.regraft_type = UNCHANGED_REGRAFT;
				}
				else{		
					if(data_on){  // Works out changes to gen_dif
						change_remove_node(t_now,nr_from,gc);  // Removes node
					
						change_add_node(t_add,nr_add,gc);      // Adds new node
				
						auto obs_node = get_gen_obs_list(n,0); 
					
						if(obs_node.size() > 0){
							vector <NodeRef> path_from, path_add;
				
							path_add.push_back(nr_add);
						
							auto n_from = n, n_add = nr_add.node;
							auto term_from = false, term_add = false;
							
							vector <ObsMut> obs_mut;
							auto mut_sum = 0u, mut_sum_root = UNSET;
								
							auto from_con = true; if(nr_from.node == OUTSIDE_INF || nr_from.node == ENTER_INF) from_con = false;
							auto add_con = true; if(n_add == OUTSIDE_INF || n_add == ENTER_INF) add_con = false;
							
							ConType con_type;
							if(from_con && add_con) con_type = FROM_CON_ADD_CON;
							if(from_con && !add_con) con_type = FROM_CON_ADD_UNCON;
							if(!from_con && add_con) con_type = FROM_UNCON_ADD_CON;
							if(!from_con && !add_con) con_type = FROM_UNCON_ADD_UNCON;
			
							switch(con_type){
							case FROM_CON_ADD_CON:               // Both from and add connected
								{	
									while(n_from != n_add){          // Works out the path back to the root
										const auto &fr_from = inf_node[n_from].from;
										if(fr_from.node == ENTER_INF || fr_from.node == OUTSIDE_INF) term_from = true;
										
										const auto &fr_add = inf_node[n_add].from;
										if(fr_add.node == ENTER_INF || fr_add.node == OUTSIDE_INF) term_add = true;
										
										if(term_from == true){
											if(term_add == true) break;
											else{ path_add.push_back(fr_add); n_add = fr_add.node;}
										}
										else{
											if(term_add == true){ path_from.push_back(fr_from); n_from = fr_from.node;}
											else{
												if(inf_node[n_from].tdiv_start < inf_node[n_add].tdiv_start){ 
													path_add.push_back(fr_add); n_add = fr_add.node;
												}
												else{
													path_from.push_back(fr_from); n_from = fr_from.node;
												}
											}
										}
									}
										
									auto plug_in = false;
									
									if(n_from == n_add && path_add.size() == 1 && t_add < t_now && path_add[0].index <= path_from[path_from.size()-1].index) plug_in = true;
										
									for(auto k = 0u; k < path_from.size(); k++){ // Goes up path_from
										const auto &nr = path_from[k];
									
										int jmax = nr.index;
										
										auto flag = false;
										
										auto jmin = 0;
										if(k == path_from.size()-1 && n_from == n_add){
											jmin = path_add[path_add.size()-1].index;
											if(plug_in == false) jmin++;
											flag = true;
										}
											
										if(jmin <= jmax){
											auto list = get_gen_obs_list(nr.node,nr.index+1,n); 
											for(auto m : list){
												ObsMut om; om.mut_num = mut_sum; om.obs = m;
												obs_mut.push_back(om);
											}
											
											const auto &inf_ev = inf_node[nr.node].inf_ev;
												
											for(int j = jmax; j >= jmin; j--){
												if(j < jmax){
													switch(inf_ev[j].type){
													case GENETIC_OBS:
														ObsMut om; om.mut_num = mut_sum; om.obs = inf_ev[j].index;
														obs_mut.push_back(om);
														break;
														
													case INFECT_OTHER:
														{
															auto list = get_gen_obs_list(inf_ev[j].index,0,n); 
															for(auto m : list){
																ObsMut om; om.mut_num = mut_sum; om.obs = m;
																obs_mut.push_back(om);
															}
														}
														break;
													}
												}
									
												if(flag == true && plug_in == true && j == jmin){
													mut_sum += inf_ev[j].mut_num - gc.mut_num;
												}												
												else{
													mut_sum += inf_ev[j].mut_num;
												}
											}
										}
									}
									
									// This is the special case of going up and plugging directly into new inf_node
									if(plug_in){
										mut_sum_root = mut_sum;
									}
									else{
										if(n_from != n_add){  // Neet to add jumping from one origin to another
											const auto &root1_from = inf_node[path_from[path_from.size()-1].node].from;
											if(root1_from.node != ENTER_INF && root1_from.node != OUTSIDE_INF){
												emsg("not root node3");		
											}
											
											mut_sum += genetic_value.inf_origin[root1_from.index].mut_num;
										
											mut_sum_root = mut_sum;
											
											const auto &root2_from = inf_node[path_add[path_add.size()-1].node].from;
											if(root2_from.node != ENTER_INF && root2_from.node != OUTSIDE_INF){
												emsg("not root node4");	
											}
											
											mut_sum += genetic_value.inf_origin[root2_from.index].mut_num;
										}
										else{
											mut_sum_root = mut_sum;							
										}
								
										for(int k = int(path_add.size())-1; k >= 0; k--){    // Goes down path_add
											const auto &nr = path_add[k];
										
											int jmin = 0;
											if(k == int(path_add.size())-1 && n_from == n_add){
												jmin = path_from[path_from.size()-1].index+1;
											}			
											
											int jmax = nr.index;
											if(jmin <= jmax){
												const auto &inf_ev = inf_node[nr.node].inf_ev;
												
												for(auto j = jmin; j <= jmax; j++){
													if(j == jmax && k == 0){
														mut_sum += gc.mut_num;
													}
													else mut_sum += inf_ev[j].mut_num;
													
													if(j < jmax){
														switch(inf_ev[j].type){
														case GENETIC_OBS:
															ObsMut om; om.mut_num = mut_sum; om.obs = inf_ev[j].index;
															obs_mut.push_back(om);
															break;
															
														case INFECT_OTHER:
															{
																auto list = get_gen_obs_list(inf_ev[j].index,0,n); 			
																for(auto m : list){
																	ObsMut om; om.mut_num = mut_sum; om.obs = m;
																	obs_mut.push_back(om);
																}
															}
															break;
														}
													}
												}
											
												if(k != 0){
													auto list = get_gen_obs_list(nr.node,nr.index+1,n); 
													for(auto m : list){
														ObsMut om; om.mut_num = mut_sum; om.obs = m;
														obs_mut.push_back(om);
													}
												}
											}
										}
								
										{ // Adds bit on add node
											auto list = get_gen_obs_list(nr_add.node,nr_add.index,n); 
											for(auto m : list){
												ObsMut om; om.mut_num = mut_sum; om.obs = m;
												obs_mut.push_back(om);
											}
										}
									}
								}
								break;
							
							case FROM_CON_ADD_UNCON:             // From connected and add unconnected
								{
									do{
										const auto &fr_from = inf_node[n_from].from;
										if(fr_from.node == ENTER_INF || fr_from.node == OUTSIDE_INF) term_from = true;
								
										if(term_from == false){ path_from.push_back(fr_from); n_from = fr_from.node;}
									}while(term_from == false);
									
									for(auto k = 0u; k < path_from.size(); k++){    // Goes up path_from
										const auto &nr = path_from[k];
									
										auto list = get_gen_obs_list(nr.node,nr.index+1); 
										for(auto m : list){
											ObsMut om; om.mut_num = mut_sum; om.obs = m;
											obs_mut.push_back(om);
										}
										
										const auto &inf_ev = inf_node[nr.node].inf_ev;
											
										int jmax = nr.index;
										for(int j = jmax; j >= 0; j--){
											if(j < jmax){
												switch(inf_ev[j].type){
												case GENETIC_OBS:
													ObsMut om; om.mut_num = mut_sum; om.obs = inf_ev[j].index;
													obs_mut.push_back(om);
													break;
													
												case INFECT_OTHER:
													{
														auto list = get_gen_obs_list(inf_ev[j].index,0); 
														for(auto m : list){
															ObsMut om; om.mut_num = mut_sum; om.obs = m;
															obs_mut.push_back(om);
														}
													}
													break;
												}
											}
											mut_sum += inf_ev[j].mut_num;
										}
									}
								
									const auto &root_from = inf_node[path_from[path_from.size()-1].node].from;
									
									if(root_from.node != ENTER_INF && root_from.node != OUTSIDE_INF){
										emsg("not root node1");
									}
									
									mut_sum += genetic_value.inf_origin[root_from.index].mut_num;
									
									mut_sum_root = mut_sum;
									
									mut_sum += gc.mut_num;
								}
								break;
							
							case FROM_UNCON_ADD_CON:             // From unconnected and add connected
								{
									do{
										const auto &fr_add = inf_node[n_add].from;
										if(fr_add.node == ENTER_INF || fr_add.node == OUTSIDE_INF) term_add = true;
								
										if(term_add == false){ path_add.push_back(fr_add); n_add = fr_add.node;}
									}while(term_add == false);
								
									if(nr_from.node != ENTER_INF && nr_from.node != OUTSIDE_INF) emsg("not root node2");
								
									mut_sum += genetic_value.inf_origin[nr_from.index].mut_num;
									
									const auto &root_from = inf_node[path_add[path_add.size()-1].node].from;
								
									if(root_from.node != ENTER_INF && root_from.node != OUTSIDE_INF){
										emsg("not root node2");
									}
									
									mut_sum_root = mut_sum;
								
									mut_sum += genetic_value.inf_origin[root_from.index].mut_num;
									
									for(int k = path_add.size()-1; k >= 0; k--){    // Goes down path_add
										const auto &nr = path_add[k];
													
										int jmax = nr.index;
										
										const auto &inf_ev = inf_node[nr.node].inf_ev;
										
										for(auto j = 0; j <= jmax; j++){
											if(j == jmax && k == 0) mut_sum += gc.mut_num;
											else mut_sum += inf_ev[j].mut_num;
											
											if(j < jmax){
												switch(inf_ev[j].type){
												case GENETIC_OBS:
													ObsMut om; om.mut_num = mut_sum; om.obs = inf_ev[j].index;
													obs_mut.push_back(om);
													break;
													
												case INFECT_OTHER:
													{
														auto list = get_gen_obs_list(inf_ev[j].index,0); 
														for(auto m : list){
															ObsMut om; om.mut_num = mut_sum; om.obs = m;
															obs_mut.push_back(om);
														}
													}
													break;
												}
											}
										}
									
										if(k != 0){
											auto list = get_gen_obs_list(nr.node,nr.index+1); 
											for(auto m : list){
												ObsMut om; om.mut_num = mut_sum; om.obs = m;
												obs_mut.push_back(om);
											}
										}
									}
									
									{ // Adds bit on add node			
										auto list = get_gen_obs_list(nr_add.node,nr_add.index); 
										for(auto m : list){
											ObsMut om; om.mut_num = mut_sum; om.obs = m;
											obs_mut.push_back(om);
										}
									}
								}
								break;
							
							case FROM_UNCON_ADD_UNCON:
								{
								}
								break;
							}
								
							if(mut_sum > 0){		
								{ // Works out observations branching from root
									auto nobs = genetic_value.obs_node_ref.size();
									vector <bool> obs_done(nobs,false);
									for(auto om : obs_mut) obs_done[om.obs] = true;
									for(auto m : obs_node) obs_done[m] = true;
									
									for(auto m = 0u; m < nobs; m++){
										if(obs_done[m] == false && genetic_value.obs_node_ref[m].node != UNSET){
											ObsMut om; om.mut_num = mut_sum_root; om.obs = m;
											obs_mut.push_back(om);
										}
									}
								}

								const auto &gen_val = model.genetic_data.gen_dif;
								const auto &gen_dif = genetic_value.gen_dif;
								
								auto &gen_dif_cha = gc.gen_dif_cha;
								
								auto sum = 0.0;
								for(const auto &om : obs_mut){     // Works out change to gen_dif
									int ch = mut_sum-2*om.mut_num;
									if(ch != 0){						
										GenDifCha gdc; 
										gdc.ob1 = om.obs;
										gdc.ch = ch;
										const auto &gen_val2 = gen_val[om.obs];
										const auto &gen_dif2 = gen_dif[om.obs];
										
										double d, d_update;
										for(auto m : obs_node){
											d = int(gen_dif2[m]) - int(gen_val2[m]);
											d_update = d + ch;
											sum += d_update*d_update - d*d;
											gdc.ob2 = m;
											
											if(gdc.ob1 == gdc.ob2) emsg("Should not be same");
								
											gen_dif_cha.push_back(gdc);
										}			
									}								
								}
							
								gc.dlike_genetic_obs -= 2*sum*model.genetic_data.obs_scale;
							}
						}
					}
				}
			}
			
			{  // Works out change as a result of alteration	
				if(data_on){
					if(t_now != t_add && in.inf_ev.size() > 0){
						const auto &iev = in.inf_ev[0];
						gc.dlike_genetic_process += poisson_probability(iev.mut_num,mut_rate*(iev.tdiv-t_add)) 
																			- poisson_probability(iev.mut_num,mut_rate*(iev.tdiv-t_now));
					}
					
					gc.dlike_genetic_obs -= GEN_OBS_WRONG*(node_alter.obs_add_end.size() + 
																					 node_alter.obs_add_begin.size());
				}
			}
		}
		break;
	
	case ADD_NODE:
		{
			auto nr_add = get_ind_noderef(iif_add.p,iif_add.i,t_add); 
			gc.nr_add = nr_add;
		
			if(data_on) change_add_node(t_add,nr_add,gc);    
		}
		break;

	case REM_NODE:
		{
			if(data_on){
				const auto &in = genetic_value.inf_node[n]; 
				change_remove_node(in.tdiv_start,in.from,gc);
			}
		}
		break;
		
	case NO_GENETIC_CHANGE:
		break;
	}
	
	return gc;
}
					

/// Adding node to new branch
void State::change_add_node(double t, const NodeRef &nr_add, GenChange &gc) const 
{
	const auto &inf_node = genetic_value.inf_node;
	auto mut_rate = genetic_value.mut_rate*model.details.dt;
	
	string warn;
	
	if(nr_add.node != OUTSIDE_INF && nr_add.node != ENTER_INF){
		const auto &in_add = inf_node[nr_add.node];
		const auto &inf_ev_add = in_add.inf_ev;
		auto j_add = nr_add.index;
		
		double t1; 
		if(j_add > 0) t1 = inf_ev_add[j_add-1].tdiv;
		else t1 = in_add.tdiv_start;
		auto t2 = t;
		
		if(j_add == inf_ev_add.size()){                // Adds at the end
			auto val = poisson_sample(mut_rate*(t2-t1),warn);
			if(val == UNSET) emsg("Problem sampling from a Poisson distribition: "+warn);
			gc.mut_num = val;
			auto prob = poisson_probability(gc.mut_num,mut_rate*(t2-t1));
			gc.dlike_genetic_process += prob;
			gc.probif += prob;
		}
		else{
			auto t3 = inf_ev_add[j_add].tdiv;
			auto frac = (t2-t1)/(t3-t1);
			auto num = inf_ev_add[j_add].mut_num;
			
			auto val = binomial_sample(frac,num,warn);
			if(val == UNSET) emsg("Problem sampling from a Poisson distribition: "+warn);
			gc.mut_num = val;
		
			gc.dlike_genetic_process += poisson_probability(gc.mut_num,mut_rate*(t2-t1)) 
																+ poisson_probability(num-gc.mut_num,mut_rate*(t3-t2))
																- poisson_probability(num,mut_rate*(t3-t1));                   
			gc.probif += binomial_probability(gc.mut_num,frac,num);
		}
	}
	else{
		auto seq_var = genetic_value.seq_var;
		
		auto lam = seq_var+mut_rate*t;
		auto val = poisson_sample(lam,warn);
		if(val == UNSET) emsg("Problem sampling from a Poisson distribition: "+warn);
		gc.mut_num = val;
		auto prob = poisson_probability(gc.mut_num,lam);
		gc.dlike_genetic_process += prob;
		gc.probif += prob;
	}
}
				
			
/// Removing node from branch
void State::change_remove_node(double t, const NodeRef &nr_from, GenChange &gc) const 		
{	
	const auto &inf_node = genetic_value.inf_node;
	auto mut_rate = genetic_value.mut_rate*model.details.dt;

	if(nr_from.node != OUTSIDE_INF && nr_from.node != ENTER_INF){ // Removing current node position
		const auto &in_from = inf_node[nr_from.node];
		const auto &inf_ev_from = in_from.inf_ev;
		auto j_from = nr_from.index;
		
		double t1; 
		if(j_from > 0) t1 = inf_ev_from[j_from-1].tdiv; 
		else t1 = in_from.tdiv_start;
		auto t2 = inf_ev_from[j_from].tdiv;
		
		if(j_from+1 == inf_ev_from.size()){            // End event
			auto num = inf_ev_from[j_from].mut_num;
			auto prob = poisson_probability(num,mut_rate*(t2-t1));
			gc.dlike_genetic_process -= prob; 
			gc.probfi += prob;
		}
		else{
			auto t3 = inf_ev_from[j_from+1].tdiv;
			auto numA = inf_ev_from[j_from].mut_num;
			auto numB = inf_ev_from[j_from+1].mut_num;

			gc.dlike_genetic_process += poisson_probability(numA+numB,mut_rate*(t3-t1)) 
																- poisson_probability(numA,mut_rate*(t2-t1)) 
																- poisson_probability(numB,mut_rate*(t3-t2)); 
			auto frac = (t2-t1)/(t3-t1);
			gc.probfi += binomial_probability(numA,frac,numA+numB);
		}
	}
	else{
		auto seq_var = genetic_value.seq_var;
	
		auto lam = seq_var+mut_rate*t;
		auto num = genetic_value.inf_origin[nr_from.index].mut_num;
	
		auto prob = poisson_probability(num,lam);
		gc.dlike_genetic_process -= prob;
		gc.probfi += prob;
	}
}


/// Gets a list of all genetic observations from the begining a given point on a node
vector <unsigned int> State::get_gen_obs_list(unsigned int n, unsigned int j, unsigned int n_not_allow) const 
{
	vector <unsigned int> list;
	get_gen_obs_list2(n,j,list,n_not_allow);
	return list;
}

/// Gets a list of all genetic observations from the begining of a node
void State::get_gen_obs_list2(unsigned int n, unsigned int j, vector <unsigned int> &list, unsigned int n_not_allow) const 
{
	const auto &inf_ev = genetic_value.inf_node[n].inf_ev;
	for(auto jj = j; jj < inf_ev.size(); jj++){
		const auto &ev = inf_ev[jj];
		switch(ev.type){
		case GENETIC_OBS: list.push_back(ev.index); break;
		case INFECT_OTHER: if(ev.index != n_not_allow) get_gen_obs_list2(ev.index,0,list,n_not_allow); break; 
		}
	}
}


/// Outputs the genetic diference (for diagnostic purposes)
void State::output_gen_dif() const
{
	cout << "OUTPUT GEN_DIF" << endl; 
	
	const auto &gen_dif = genetic_value.gen_dif;
	const auto &gen_val = model.genetic_data.gen_dif;

	auto nobs = gen_dif.size();
		
	ofstream fout("gen_dif.txt");
	fout << like.genetic_obs << " gen_obs" << endl;
	for(auto j = 0u; j < nobs; j++){
		for(auto i = 0u; i < nobs; i++){
			fout << gen_dif[j][i] << "(" << gen_val[j][i] <<"  ";
		}
		fout << endl;
	}
	
	const auto &inf_node = genetic_value.inf_node;
	
	ofstream fout2("node.txt");
	for(auto p = 0u; p < model.species.size(); p++){
		auto &ssp = species[p];
		for(auto i = 0u; i < ssp.individual.size(); i++){
			for(auto n = 0u; n < inf_node.size(); n++){
				auto &in = inf_node[n];
				if(in.p == p && in.i == i){
					fout2 << i << " " << in.tdiv_start << "  ";
					for(auto &ev : in.inf_ev){
						fout2 << ev.tdiv << " ";
						switch(ev.type){
						case INFECT_OTHER: fout2 << "infect " << ev.index << " "; break;
						case GENETIC_OBS: fout2 << "obs "; break;
						}
						fout2 << ev.mut_num << "  ";
					}
					fout2 << endl;
				}
			}
		}
	}
}


/// Output swho acquired infection from whom (for diagnostic purposes)
void State::output_waifw(string file) const 
{	
	ofstream fout(file);
	fout << "i,p,i,pref,t" << endl;
	for(auto p = 0u; p < model.species.size(); p++){
		auto &ssp = species[p];
		for(auto i = 0u; i < ssp.individual.size(); i++){
			for(auto &ev : ssp.individual[i].ev){
				auto &iif  = ev.ind_inf_from;
				if(iif.p != UNSET){
					fout << i << "," << iif.p << "," << iif.i << "," << iif.pref << "," << ev.tdiv << endl;
				}
			}
		}
	}
}


/// Recalculates iif.w for weights on trans_tree infection events
vector <double> State::recalculate_iif_w()
{
	vector <double> store;
	
	const auto &inf_node = genetic_value.inf_node;
	
	for(auto &in : inf_node){
		auto node = in.from.node;
		if(node != ENTER_INF && node != OUTSIDE_INF){
			auto &iif = species[in.p].individual[in.i].ev[in.e].ind_inf_from;
		
			const auto &po = model.pop[iif.po];
			const auto &ind = species[iif.p].individual[iif.i];

			auto num = 1.0; 
			if(po.ind_variation){
				for(auto ie : po.ind_eff_mult) num *= ind.exp_ie[ie];
				for(auto fe : po.fix_eff_mult) num *= ind.exp_fe[fe];
			}
			
			store.push_back(iif.w);
			
			iif.w = num; 
		}
	}
	
	return store;
}


/// Undos changes to iif_w
void State::recalculate_iif_w_undo(const vector <double> &vec)
{
	const auto &inf_node = genetic_value.inf_node;
	
	auto k = 0u;
	for(auto &in : inf_node){
		auto node = in.from.node;
		if(node != ENTER_INF && node != OUTSIDE_INF){
			auto &iif = species[in.p].individual[in.i].ev[in.e].ind_inf_from;
		
			iif.w = vec[k]; k++;
		}
	}
	if(k != vec.size()) emsg("Wrong size");
}


/// Recalculates popnum_ind.w for a given population
vector <double> State::recalculate_popnum_ind_w(unsigned int po)
{
	vector <double> store;

	const auto &pop = model.pop[po];
	
	const auto &individual = species[pop.sp_p].individual;
	
	for(auto ti = 0u; ti < T; ti++){
		for(auto &pir : popnum_ind[ti][po]){
			const auto &ind = individual[pir.i];
			auto num = 1.0; 
			for(auto ie : pop.ind_eff_mult) num *= ind.exp_ie[ie];
			for(auto fe : pop.fix_eff_mult) num *= ind.exp_fe[fe];
			store.push_back(pir.w);
			pir.w = num;
		}
	}

	return store;
}


/// Undos changes to popnum_ind.w
void State::recalculate_popnum_ind_w_undo(unsigned int po, const vector <double> &vec)
{
	auto k = 0u;
	for(auto ti = 0u; ti < T; ti++){
		for(auto &pir : popnum_ind[ti][po]){
			pir.w = vec[k]; k++;
		}
	}
}


/// Updates the likelihood change and dprob based on changes to transmission tree
void GenChange::update_like_ch(Like &like_ch, double &dprob)
{
	if(type != NO_GENETIC_CHANGE){
		dprob += probfi-probif;
									
		like_ch.genetic_process += dlike_genetic_process;
		like_ch.genetic_obs += dlike_genetic_obs; 
		
		dlike_genetic_process = 0;
		dlike_genetic_obs = 0;
	}
}


/// If the update is accepted then this updates state
void State::gen_change_update(const GenChange &gc) 
{
	if(gc.type == NO_GENETIC_CHANGE) return;
	
	auto &inf_node = genetic_value.inf_node;
	const auto &na = gc.node_alter;	
	
	auto n = gc.n;
	
	// Need to create a new node
	if(gc.type == ADD_NODE){
		n = inf_node.size();
		inf_node.push_back(gc.node_add);
	}
		
	if(na.e_add != UNSET) inf_node[n].e = na.e_add;
		
	// Makes alterations to the connecting infnode
	auto &in = inf_node[n];
														
	in.tdiv_start = na.tdiv_start;
	in.tdiv_rec = na.tdiv_rec;
	
	string warn;
	
	// Adds observations
	if(na.obs_add_begin.size() > 0 || na.obs_add_end.size() > 0){
		auto &inf_ev = in.inf_ev;
			
		for(auto k = 0u; k < na.obs_add_begin.size(); k++){ // Adds observations at beginning
			auto m = na.obs_add_begin[k];
			const auto &obs = model.genetic_data.obs[m];
			
			InfEvent iev; iev.tdiv = obs.tdiv; iev.type = GENETIC_OBS; iev.index = m; iev.mut_num = 0;
			inf_ev.insert(inf_ev.begin()+k,iev);
			genetic_value.nobs_not_infected--;
		}
			
		auto mut_rate = genetic_value.mut_rate*model.details.dt;
		
		for(auto m : na.obs_add_end){                  // Adds observations at the end
			const auto &obs = model.genetic_data.obs[m];
			
			double t_last;
			if(inf_ev.size() == 0) t_last = in.tdiv_start; 
			else t_last = inf_ev[inf_ev.size()-1].tdiv;

			auto num = poisson_sample(mut_rate*(obs.tdiv-t_last),warn);
			if(num == UNSET) emsg("Problem sampling from a Poisson distribition: "+warn);
				
			InfEvent iev; iev.tdiv = obs.tdiv; iev.type = GENETIC_OBS; iev.index = m; iev.mut_num = num;
			inf_ev.push_back(iev);
			
			genetic_value.nobs_not_infected--;
		}

		update_inf_ev(n);
	}
	
	switch(gc.type){
	case REGRAFT_NODE:
		switch(gc.regraft_type){ // The reconnection position is the same place on the tree
		case NORMAL_REGRAFT:
			disconnect_inf_node(n);
			connect_inf_node(n,gc.iif_add,gc.mut_num);
			break;
			
		case ORIGIN_REGRAFT:
			{
				auto &nr_from = in.from;
				if(na.tdiv_start == 0) nr_from.node = ENTER_INF;
				else nr_from.node = OUTSIDE_INF;
			}
			break;
			
		case UNCHANGED_REGRAFT:
			{
				const auto &nr_from = in.from;
				if(nr_from.node != ENTER_INF){
					inf_node[nr_from.node].inf_ev[nr_from.index].tdiv = na.tdiv_start;	
				}
			}
			break;
		}
		break;
		
	case ADD_NODE:
		connect_inf_node(n,gc.iif_add,gc.mut_num);
		break;
		
	case REM_NODE:
		disconnect_inf_node(n);
		remove_node(n);
		break;
		
	case NO_GENETIC_CHANGE: break;
		
	case GENCHA_FAIL: emsg("Should not be here"); break;
	}
	
	auto &gen_dif = genetic_value.gen_dif;
		
	// If new observations have appeared then recalculate genetic difference and likelihood
	if(na.obs_add_end.size() > 0 || na.obs_add_begin.size() > 0){
		gen_dif = calculate_gen_dif();
		like.genetic_obs = likelihood_genetic_obs(gen_dif);
		like.genetic_process = likelihood_genetic_process();
	}
	else{
		like.genetic_process += gc.dlike_genetic_process;
		like.genetic_obs += gc.dlike_genetic_obs;
	
		for(const auto &cha : gc.gen_dif_cha){
			gen_dif[cha.ob1][cha.ob2] += cha.ch;
			gen_dif[cha.ob2][cha.ob1] += cha.ch;
		}
	}
}


/// Removes node from the system 
void State::remove_node(unsigned int n)
{
	auto &inf_node = genetic_value.inf_node;
	auto &inf_origin = genetic_value.inf_origin;
	
	if(n+1 < inf_node.size()){
		auto &in = inf_node[n];
		in = inf_node[inf_node.size()-1];
	
		species[in.p].individual[in.i].ev[in.e].inf_node_ref = n;
		
		const auto &from = in.from;
		if(from.node == ENTER_INF || from.node == OUTSIDE_INF){
			inf_origin[from.index].node = n;
		}
		else{
			inf_node[from.node].inf_ev[from.index].index = n;
		}
		
		for(auto &iev : in.inf_ev){
			switch(iev.type){
			case INFECT_OTHER:
				inf_node[iev.index].from.node = n;
				break;
			
			case GENETIC_OBS:
				genetic_value.obs_node_ref[iev.index].node = n;
				break;
			}
		}
	}
	inf_node.pop_back();
}
	

/// Updates all the references along a node
void State::update_inf_ev(unsigned int n)
{
	auto &inf_ev = genetic_value.inf_node[n].inf_ev;
	
	for(auto k = 0u; k < inf_ev.size(); k++){
		NodeRef nr; nr.node = n; nr.index = k;
			
		const auto &iev = inf_ev[k];
		switch(iev.type){
		case GENETIC_OBS: genetic_value.obs_node_ref[iev.index] = nr; break;
		case INFECT_OTHER: genetic_value.inf_node[iev.index].from = nr; break;
		}
	}
}


/// Calculates popnum_ind
void State::calculate_popnum_ind()
{
	popnum_ind.clear();
	popnum_ind.resize(T);
	for(auto ti = 0u; ti < T; ti++){
		popnum_ind[ti].resize(model.pop.size());
	}
	
	for(auto p = 0u; p < model.species.size(); p++){
		const auto &sp = model.species[p];
		if(sp.trans_tree){
			auto &ssp = species[p];
			for(auto i = 0u; i < ssp.individual.size(); i++){			
				add_popnum_ind(p,i);		
			}
		}
	}
}
				

/// Adds an individual to popnum_ind
void State::add_popnum_ind(unsigned int p, unsigned int i)
{
	auto &ind = species[p].individual[i];

	const auto &comp_gl = model.species[p].comp_gl;

	ind.popnum_ind_ref.clear();
	
	auto e = 0u;

	auto c = UNSET;
	for(auto ti = 0u; ti < T; ti++){
		while(e < ind.ev.size() && ind.ev[e].tdiv <= ti){
			c = ind.ev[e].c_after; e++;
		}
		
		if(c != UNSET){
			for(const auto &po : comp_gl[c].pop_ref_simp){	
				const auto &pop = model.pop[po];
			
				auto num = 1.0; 
				for(auto ie : pop.ind_eff_mult) num *= ind.exp_ie[ie];
				for(auto fe : pop.fix_eff_mult) num *= ind.exp_fe[fe];

				PopIndRef pir;
				pir.i = i;
				pir.w = num;
				pir.index = ind.popnum_ind_ref.size();
				
				PopnumIndRef pnir;
				pnir.po = po;
				pnir.ti = ti;
				pnir.index = popnum_ind[ti][po].size();
				
				popnum_ind[ti][po].push_back(pir);
				ind.popnum_ind_ref.push_back(pnir);
			}
		}
	}
}
			
			
/// Updates popnum_ind after a proposal has been accepted
void State::update_popnum_ind(unsigned int p, unsigned int i)
{
	auto &ind = species[p].individual[i];
	
	// Removes current references
	for(const auto &pir : ind.popnum_ind_ref){
		auto &pi = popnum_ind[pir.ti][pir.po];
		auto k = pir.index;
		if(k+1 < pi.size()){
			auto &pi_end = pi[pi.size()-1];
			pi[k] = pi_end;
			auto pp = model.pop[pir.po].sp_p;
			species[pp].individual[pi_end.i].popnum_ind_ref[pi_end.index].index = k;
		}		 
		pi.pop_back();	
	}
	
	// Adds new references
	add_popnum_ind(p,i);
}


/// Updates popnum_ind after a proposal has been accepted
void State::popnum_ind_recalc_w(unsigned int p, unsigned int i)
{
	auto &ind = species[p].individual[i];
	
	// Removes current references
	for(const auto &pir : ind.popnum_ind_ref){
		const auto &pop = model.pop[pir.po];
			
		auto num = 1.0; 
		for(auto ie : pop.ind_eff_mult) num *= ind.exp_ie[ie];
		for(auto fe : pop.fix_eff_mult) num *= ind.exp_fe[fe];

		popnum_ind[pir.ti][pir.po][pir.index].w = num;
	}
}

// Functions related to genetic observations

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
using namespace std;

#include "state.hh"
#include "matrix.hh"
#include "utils.hh"

/// Used to time order infection nodes
bool InfNode_ord (const InfNode &in1, const InfNode &in2)                      
{ return (in1.tdiv_start < in2.tdiv_start); };  


/// Used to time order infection events
bool InfEv_ord (const InfEvent &ev1, const InfEvent &ev2)                      
{ return (ev1.tdiv < ev2.tdiv); };  


/// Given an initial event sequence works out how to adapt the transmission tree
TransTreeChange State::update_tree(unsigned int p, unsigned int i, vector <Event> &ev_new)
{
	auto pl = false;

	const auto &sp = model.species[p];
	
	if(sp.trans_tree == false) return no_trans_tree_change();
	
	auto &inf_node = genetic_value.inf_node;
	
	const auto &ev_old = species[p].individual[i].ev;
	
	// node_old stores all the inf_node currently attached to the individual
	
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
	
	if(node_old.size() > 1) run_error("BICI can only currently perform inference with pathogen genetic data when individuals cannot become reinfected.");
		
	// node_new stores events where new nodes need to be attached
	
	vector <unsigned int> node_new;
	
	for(auto e = 0u; e < ev_new.size(); e++){
		const auto &ev = ev_new[e];
		switch(ev.type){
		case ENTER_EV:
			if(sp.comp_gl[ev.c_after].infected == true){
				node_new.push_back(e);
			}
			break;
			
		case M_TRANS_EV:
			if(sp.tra_gl[ev.tr_gl].infection.type == TRANS_INFECTION){
				node_new.push_back(e);
			}
			break;
			
		default: break;
		}
	}
	
	if(node_new.size() > 1) run_error("BICI can only currently perform inference with pathogen genetic data when individuals cannot become reinfected.");
	
	// Works out how to connect nodes
	
	switch(node_old.size()){
	case 0: 
		switch(node_new.size()){
		case 0:                                        // No infection -> no infection
			return no_trans_tree_change();
	
		case 1:                                        // No infection -> infection
			{		
				// Find where infection comes from
				auto e_add = node_new[0];
				
				auto &ev = ev_new[e_add];
				
				auto probif = sample_infection_source(ev,p);
				
				if(probif == UNSET) return trans_tree_fail(); // Fails to find infection
					
				if(testing){
					auto prob = prob_infection_source(ev,p);
					if(dif(prob,probif,dif_thresh)){
						emsg("Source sampler problem here");
					}
				}
			
				auto ina = get_inf_node_alteration(p,i,UNSET,ev_new,e_add);
				if(ina.possible == true){ 
					auto gc = trans_tree_change(ADD_NODE,UNSET,ina,ev.ind_inf_from);
				
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
		switch(node_new.size()){
		case 0:                                         // Infection -> no infection
			{		
				auto n = node_old[0];
				const auto &in = inf_node[n];
				
				if(in.inf_ev.size() == 0){         // Only remove if no observations or onward infections
					const auto &ev = species[in.p].individual[in.i].ev[in.e];
					
					auto probfi = prob_infection_source(ev,in.p);
					
					auto gc = trans_tree_change(REM_NODE,n,InfNodeAlter(),ev.ind_inf_from);
					
					gc.probfi += probfi;
					return gc;
				}
			}
			break;
	
		case 1:                                         // Infection -> Infection
			{		
				auto n = node_old[0];
				const auto &in_old = inf_node[n];
				
				auto e_add = node_new[0];
				auto &ev = ev_new[e_add];
				
				auto t_old = in_old.tdiv_start;
				auto t_new = ev.tdiv;
		
				auto ina = get_inf_node_alteration(p,i,n,ev_new,e_add);
				if(ina.possible == true){ 
					if(ina.unchanged){
						ev = ssp.individual[in_old.i].ev[in_old.e];
						return no_trans_tree_change();
					}
					
					auto probif = 0.0, probfi = 0.0;
				
					// Finds where the attachment is going to occur
					// Accounts for the branching in probif and probfi
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
					if(nn == UNSET) return trans_tree_fail();
					
					// Works out ind_inf_from for the connection
					
					IndInfFrom iif_prop;
					
					if(nn == ENTER_INF || nn == OUTSIDE_INF){
						auto nn_new = ENTER_INF;
						if(ev.type != ENTER_EV) nn_new = OUTSIDE_INF;
							
						if(nn_new != nn && nn_new == OUTSIDE_INF){
							if(trg_contains_outside(in_old.p,ev.tr_gl) == false) return trans_tree_fail();	
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
						if(t_new < in_new.ti_start || t_new > in_new.ti_rec){
							return trans_tree_fail();
						}
						
						// Works out IndInfFrom for new position on tree
						
						auto pp = in_new.p;
						
						if(!set_iif(pp,in_new.i,in_new.e,iif_prop,t_new,ev,sp)){
							return trans_tree_fail();
						}
					}
					ev.ind_inf_from = iif_prop;
					ev.inf_node_ref = n;
					
					auto gc = trans_tree_change(REGRAFT_NODE,node_old[0],ina,iif_prop);
				
					gc.probif += probif;
					gc.probfi += probfi;
					return gc;
				}			
			}
			break;
		}
		break;
	}
	
	return trans_tree_fail();
}


/// Calculcate extra properties of ind_inf_from
bool State::set_iif(unsigned int p, unsigned int i, unsigned int e, IndInfFrom &iif, double t, const Event &ev, const Species &sp) const 
{
	iif.p = p;
	iif.i = i;
					
	// Get the c of new infecting individual
	auto t_inf = round_down(t);
	const auto &ev_inf = species[p].individual[i].ev;
	
	auto c = UNSET;
	while(e < ev_inf.size() && ev_inf[e].tdiv <= t_inf){ c = ev_inf[e].c_after; e++;}
	if(c == UNSET) return false;
	
	// Works out pref from c
	const auto &tra = sp.tra_gl[ev.tr_gl];
	auto m = tra.markov_eqn_ref;
	const auto &me = sp.markov_eqn[m];
	const auto &eq = model.eqn[me.eqn_ref];
	
	auto pref = eq.comp_pref_convert[p][c];
	if(pref == UNSET) return false;
	
	iif.pref = pref;
	iif.po = eq.pop_ref[pref];
	iif.w = get_w_from_indinffrom(iif);
	
	return true;
}


/// Determines is a transition contains a term for outside infection
bool State::trg_contains_outside(unsigned int p, unsigned int tr_gl) const
{
	const auto &sp = model.species[p];
	
	const auto &tra = sp.tra_gl[tr_gl];
	if(tra.infection.type != TRANS_INFECTION) emsg("Should be transition");
	const auto &me = sp.markov_eqn[tra.markov_eqn_ref];
	const auto &eq = model.eqn[me.eqn_ref];
	const auto &lin = eq.lin;
	
	if(lin.no_pop_precalc.type == ZERO) return false;
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
			
		const auto &lin = eq.lin;
			
		if(!lin.on) emsg("Linearisation should be on");
			
		// First selects if infection comes from outside or from which population
			
		auto Npop = eq.pop_ref.size();
			
		auto j = UNSET;
			
		if(lin.multi_source){        // Samples from available sources (either populations of from outside)		
			auto ss = eq.setup_source_sampler(ti,popnum_t[ti],param_val,popcombw_value);	
			j = ss.sample_inf_source();
			if(j == UNSET) return UNSET;
			
			probif += ss.prob_inf_source(j);
		}
		else{
			if(testing){
				if(Npop > 1) emsg("Npop wrong");
				if(Npop == 1 && lin.no_pop_precalc.type != ZERO) emsg("Npop wrong");
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
		
			iif.p = pop.p;
			
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
		
	const auto &precalc = param_val.precalc;
	
	const auto &lin = eq.lin;
	
	auto dt = model.details.dt;
	
	if(!lin.on) emsg("Linearisation should be on");
		
	// First selects if infection comes from outside or from which population
		
	auto probfi = 0.0;
	
	auto pref = iif.pref;
	
	//const auto &precalc = param_val.precalc;
	
	if(lin.multi_source){        // Prop from available sources (either populations or outside)		
		auto value = mev.value_t[ti];
		if(value < TINY) probfi -= LARGE;
		else{		
			double val;
			if(pref == UNSET) val = dt*eq.calculate_no_pop(ti,precalc);
			else{
				auto pr = eq.pop_ref[pref];
				val = popnum_t[ti][pr]*dt*eq.calculate_pop_grad(pref,ti,popcombw_value,precalc);
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
// E.g. Time range changem, onward transmission change, adding in genetic observations
InfNodeAlter State::get_inf_node_alteration(unsigned int p, unsigned int i, unsigned int n, const vector <Event> &ev_new, unsigned int e_add) const
{			
	const auto &sp = model.species[p]; 
	const auto &gen_data = model.genetic_data;
	const auto &inf_node = genetic_value.inf_node;

	auto on = gen_data.on; if(i >= sp.nindividual_obs) on = false;
	
	// Works out start and end times for infection
	
	const auto &evn = ev_new[e_add];
	auto nev = ev_new.size();
	
	auto t_start = evn.tdiv;
	double t_rec = model.details.T;
	for(auto ee = e_add+1; ee < nev; ee++){
		const auto &ev = ev_new[ee];
		if(ev.type == M_TRANS_EV || ev.type == NM_TRANS_EV){
			if(sp.tra_gl[ev.tr_gl].infection.type == TRANS_RECOVERY){
				t_rec = ev.tdiv;
				break;
			}
		}
	}

	auto ti_start = inf_node_ti_from_tdiv(t_start);
	auto ti_rec = inf_node_ti_from_tdiv(t_rec);
	
	InfNodeAlter ina;
	ina.e_add = e_add;
	
	ina.tdiv_start = t_start;
	ina.tdiv_rec = t_rec;
	ina.ti_start = ti_start;
	ina.ti_rec = ti_rec;
	
	if(n != UNSET){
		const auto &in_old = inf_node[n];
	
		// Makes sure that any onward infections are within the time span
		
		const auto &inf_ev = in_old.inf_ev;
		
		auto c = evn.c_after;
	
		auto ei = e_add+1;
		for(auto k = 0u; k < inf_ev.size(); k++){
			const auto &iev = inf_ev[k];
			if(iev.type == INFECT_OTHER){
				auto t = iev.tdiv;
				if(t < ti_start || t > ti_rec){
					ina.possible = false;
				}
				else{ // Looks for any changes in which population is causing infection
					auto t_do = round_down(t);
					while(ei < nev && ev_new[ei].tdiv <= t_do){ c = ev_new[ei].c_after; ei++;}
					
					const auto &no_inf = inf_node[iev.index];
					
					auto pp = no_inf.p;
					auto ii = no_inf.i;
					auto ee = no_inf.e;
					const auto &ev = species[pp].individual[ii].ev[ee];
					
					const auto &inf_c = model.inf_cause[pp][ev.tr_gl][p][c];
					
					if(inf_c.eq_ref == UNSET) ina.possible = false;
					else{
						const auto &iif = ev.ind_inf_from;
						
						if(iif.pref != inf_c.pref){
							InfPoCha ipc;
							ipc.p = pp;
							ipc.i = ii;
							ipc.e = ee;
							ipc.pref_new = inf_c.pref;
							ipc.po_new = inf_c.po;
							ina.ind_po_cha.push_back(ipc);
						}
					}
				}
			}
		}
	
		// Looks to add in any observation at begining or end
		if(on && !(t_start == in_old.tdiv_start && t_rec == in_old.tdiv_rec)){
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


/// Sets up inf_node
void State::setup_inf_node()
{
	auto &inf_node = genetic_value.inf_node;
	inf_node.clear();

	for(auto p = 0u; p < model.nspecies; p++){ 
		const auto &sp = model.species[p];
		auto &ssp = species[p];

		for(auto i = 0u; i < ssp.individual.size(); i++){
			auto &ind = ssp.individual[i];

			for(auto e = 0u; e < ind.ev.size(); e++){
				auto &ev = ind.ev[e];
				ev.inf_node_ref = UNSET;
				auto &iif = ev.ind_inf_from;
				
				switch(ev.type){
				case ENTER_EV:
					if(sp.comp_gl[ev.c_after].infected == true){
						InfNode in; 
						in.tdiv_start = ev.tdiv; in.p = p; in.i = i; in.e = e; 
						in.tdiv_rec = find_t_rec(in);
						in.ti_start = inf_node_ti_from_tdiv(in.tdiv_start);
						in.ti_rec = inf_node_ti_from_tdiv(in.tdiv_rec);
						in.from.node = UNSET; in.from.index = UNSET;
						inf_node.push_back(in);
						
						iif.p = ENTER_INF;
					}
					else{
						if(iif.p != UNSET) emsg("ind_inf_from.p should not be setb");
					}
					break;
				
				case MOVE_EV:
					{
						const auto &ev_last = ind.ev[e-1];
						if(sp.comp_gl[ev_last.c_after].infected == false &&
						   sp.comp_gl[ev.c_after].infected == true){
							
							InfNode in;
							in.p = p; in.i = i; in.e = e;
							in.tdiv_start = ev.tdiv; 							
							in.tdiv_rec = find_t_rec(in);
							in.ti_start = inf_node_ti_from_tdiv(in.tdiv_start);
							in.ti_rec = inf_node_ti_from_tdiv(in.tdiv_rec);
							in.from.node = UNSET; in.from.index = UNSET;
							inf_node.push_back(in);
						
							iif.p = ENTER_INF;
						}
						else{
							if(ev.ind_inf_from.p != UNSET) emsg("ind_inf_from.p should not be setb");
						}
					}
					break;
					
				case M_TRANS_EV: 
					{	
						const auto &trg = sp.tra_gl[ev.tr_gl];
						if(trg.infection.type == TRANS_INFECTION){
							InfNode in; 
							in.p = p; in.i = i; in.e = e; 
							in.tdiv_start = ev.tdiv; 
							in.tdiv_rec = find_t_rec(in);
							in.ti_start = inf_node_ti_from_tdiv(in.tdiv_start);
							in.ti_rec = inf_node_ti_from_tdiv(in.tdiv_rec);
							in.from.node = UNSET; in.from.index = UNSET;
							inf_node.push_back(in);
						}
					}
					break;
					
				case NM_TRANS_EV:
					{
						const auto &trg = sp.tra_gl[ev.tr_gl];
						if(trg.infection.type == TRANS_INFECTION) emsg("Non-Markovian transition cannot be used for infection");
					}
					break;
					
				default: break;
				}
			}
		}
	}
	
	sort(inf_node.begin(),inf_node.end(),InfNode_ord);    
	
	for(auto n = 0u; n < inf_node.size(); n++){
		const auto &in = inf_node[n];
		auto p = in.p;
		if(p != ENTER_INF && p != OUTSIDE_INF){
			species[in.p].individual[in.i].ev[in.e].inf_node_ref = n;
		}
	}
	
	
	if(false){	
		for(auto n = 0u; n < inf_node.size(); n++){
			cout << "node " << n << endl;
			const auto &no = inf_node[n];
			cout << no.tdiv_start << " " << no.tdiv_rec << " time range" << endl;
		}
		emsg("inf node");
	}
}

	
/// Sets up the transition tree
void State::setup_transtree(bool set_iif)
{
	if(!model.trans_tree && !model.trans_tree_output) return;
	
	if(trans_tree_on != false) emsg("Trans tree should not be set");
	
	trans_tree_on = true;
	
	// Intialises inf_node (without doing connections between nodes)
	setup_inf_node();

	// Creates ind_inf_from on individuals 
	if(set_iif){
		if(model.genetic_data.on){
			setup_iif_with_genetic_data();
		}
		else{
			auto pop_ind = calculate_pop_ind_total(); 
			for(auto p = 0u; p < species.size(); p++){
				const auto &sp = model.species[p];
				if(sp.trans_tree || sp.trans_tree_output) species[p].setup_iif(popnum_t,pop_ind);
			}
		}		
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
	
	auto &inf_node = genetic_value.inf_node;
	const auto &gen_data = model.genetic_data;
	
	if(false){
		cout << "OUTPUT" << endl;
		for(auto p = 0u; p < model.nspecies; p++){ 
			auto &ssp = species[p];	
			for(auto i = 0u; i < ssp.individual.size(); i++){
				const auto &ind = ssp.individual[i];
				cout << ind.name << ": " << endl;
				ssp.print_event(ind.ev);
			}
		}
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
			
			if(pp == UNSET) emsg("Should not be unset3");
		
			if(pp == ENTER_INF || pp == OUTSIDE_INF){	
				in.from.node = pp;
				in.from.index = genetic_value.inf_origin.size();
			
				InfOrigin io; io.node = k; io.mut_num = UNSET;
		
				genetic_value.inf_origin.push_back(io);
			}
			else{
				auto kk = ind_node_ref[pp][eve.ind_inf_from.i];
			
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

	// If genetic data then adds and samples mutation numbers 
	if(model.genetic_data.on){
		set_genetic_param();
		sample_genetic_value();
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


/// Gets the nearest genetic observation to a given 
unsigned int State::get_nearest_gobs(double tdiv, const vector <unsigned int> &list) const
{
	if(list.size() == 0) return UNSET;
	
	const auto &gen_data = model.genetic_data;
	auto dmin = LARGE;
	auto msel = UNSET;
	for(auto m : list){
		auto dt = gen_data.obs[m].tdiv-tdiv;
		if(dt*dt < dmin){
			dmin = dt*dt;
			msel = m;
		}
	}
	
	return msel;	
}

		
/// Sets up ind_inf_from using genetic observations (this is used to 
void State::setup_iif_with_genetic_data()
{
	auto &inf_node = genetic_value.inf_node;
		
	auto nnode = inf_node.size();

	// Find all potential links between individuals
	vector < vector <PossLink> > poss_link; // [node][#]
	poss_link.resize(nnode);
	
	// Creates ind_inf_from on individuals 	
	auto pop_ind = calculate_pop_ind_total(); 
	for(auto p = 0u; p < species.size(); p++){
		const auto &sp = model.species[p];
		if(sp.trans_tree || sp.trans_tree_output){
			species[p].calc_poss_link(poss_link,popnum_t,pop_ind);
		}
	}		
	
	// Stores connections between individuals and observations
	const auto &gen_data = model.genetic_data;

	vector < vector < vector <unsigned> > > connected_gobs;  
	connected_gobs.resize(model.species.size());
	for(auto p = 0u; p < model.species.size(); p++){
		const auto &individual = species[p].individual;
		connected_gobs[p].resize(individual.size());
		
		for(auto i = 0u; i < individual.size(); i++){
			for(const auto &ob : gen_data.ind_gen_obs[p][i]){
				connected_gobs[p][i].push_back(ob.m);
			}
		}
	}

	for(int n = nnode-1; n >= 0; n--){
		const auto &no = inf_node[n];
		auto p = no.p, i = no.i;
		
		auto &iif = species[p].individual[i].ev[no.e].ind_inf_from; 
		
		if(iif.p != OUTSIDE_INF && iif.p != ENTER_INF){
			auto distmin = LARGE;
			auto plsel = UNSET;
			
			const auto &possl = poss_link[n];
			
			auto nposs = possl.size();
			if(nposs == 0) emsg("nposs prob");
			
			auto gobs = get_nearest_gobs(no.tdiv_start,connected_gobs[p][i]);
			if(gobs != UNSET){
				for(auto j = 0u; j < nposs; j++){
					const auto &pl = possl[j];
					const auto &iif = pl.iif;
					if(iif.p != OUTSIDE_INF && pl.prob > 0){
						auto gobs_fr = get_nearest_gobs(no.tdiv_start,connected_gobs[iif.p][iif.i]);
						if(gobs_fr != UNSET){
							auto d = gen_data.gen_dif[gobs][gobs_fr];
							if(d < distmin){ distmin = d; plsel = j;}
						}
					}
				}
			}
			
			if(plsel == UNSET){   // If not selected then randomly select in prortion to prob
				vector <double> prob_sum;
				auto sum = 0.0;
				for(auto j = 0u; j < nposs; j++){
					sum += possl[j].prob;
					prob_sum.push_back(sum);
				}		
				if(sum == 0) emsg("zero prob");
			
				auto z = ran()*sum;
				plsel = 0u; while(plsel < nposs && z > prob_sum[plsel]) plsel++;
				if(plsel == nposs) emsg("nposs prob");
			}
			
			iif = possl[plsel].iif;
		
			if(iif.p != OUTSIDE_INF){
				for(auto m : connected_gobs[p][i]){
					connected_gobs[iif.p][iif.i].push_back(m);
				}
			}
		}
	}
	
	if(false){
		for(auto p = 0u; p < model.species.size(); p++){
			const auto &individual = species[p].individual;
			for(auto i = 0u; i < individual.size(); i++){
				cout << individual[i].name << ":";
				for(auto m : connected_gobs[p][i]) cout << m << ",";
				cout << endl;
			}	
		}
	}
}


/// Finds the recovery time for an infection node
double State::find_t_rec(const InfNode &in) const
{
	const auto &sp = model.species[in.p];
	const auto &ev = species[in.p].individual[in.i].ev;
	 
	for(auto e = in.e; e < ev.size(); e++){
		if((ev[e].type == M_TRANS_EV || ev[e].type == NM_TRANS_EV) && sp.tra_gl[ev[e].tr_gl].infection.type == TRANS_RECOVERY){
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
	
	const auto &precalc = param_val.precalc;
	
	const auto &gen_data = model.genetic_data;
	genetic_value.mut_rate = model.eqn[gen_data.mut_rate.eq_ref].calculate_param(precalc);
	
	genetic_value.seq_var = model.eqn[gen_data.seq_var.eq_ref].calculate_param(precalc);

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
	
	if(model.genetic_data.on){
		minimise_mut_num();
	}
	
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


/// Minimises mut_num for the initial state
void State::minimise_mut_num()
{
	const auto &gen_data = model.genetic_data;
	
	auto &gen_dif = genetic_value.gen_dif;
	const auto &gen_val = gen_data.gen_dif;
	
	auto glike = likelihood_genetic_obs(gen_dif);

	auto &inf_node = genetic_value.inf_node;
	
	auto loopmax = 1000u;
	for(auto loop = 0u; loop < loopmax; loop++){
		auto fl = false;
		for(auto n = 0u; n < inf_node.size(); n++){
			auto &no = inf_node[n];
			for(auto e = 0u; e < no.inf_ev.size(); e++){
				auto obs_node = get_gen_obs_list(n,e);	
				auto nobs = genetic_value.obs_node_ref.size();
				auto obs_other = get_not_obs(obs_node,vector <unsigned int> (),nobs);
				
				auto lin = 0.0, quad = 0.0;
				for(auto m1 : obs_node){
					const auto &gen_dif2 = gen_dif[m1];
					const auto &gen_val2 = gen_val[m1];
					for(auto m2 : obs_other){
						auto d = int(gen_dif2[m2]) - int(gen_val2[m2]);
						lin += 2*d;
						quad++;
					}
				}
			
				if(quad > 0){
					lin *= -gen_data.obs_scale*2;
					quad *= -gen_data.obs_scale*2;
						
					auto &mut_num = no.inf_ev[e].mut_num;
					auto mut_num_new = round_int(mut_num-lin/(2*quad));
					
					if(mut_num_new < 0) mut_num_new = 0;
					
					if(mut_num_new != (int)mut_num){ 
						fl = true;
						int dn = mut_num_new - mut_num;
						auto dL = dn*lin + quad*dn*dn;
					
						mut_num = mut_num_new;
						
						for(auto m1 : obs_node){
							for(auto m2 : obs_other){
								gen_dif[m1][m2] += dn;
								gen_dif[m2][m1] += dn;
							}
						}
						
						glike += dL;
					}
				}
			}
		}
		
		if(fl == false) break;
	}
	
	if(testing){
		auto gen_dif_check = calculate_gen_dif();
		if(difi(gen_dif_check,gen_dif)) emsg("gen_dif_prob");
		
		auto glike_check = likelihood_genetic_obs(gen_dif);
		if(dif(glike_check,glike,0.1)) emsg("glike prob");
	}
}

	
/// Calculates the difference in genetic sequences
vector < vector <unsigned int> > State::calculate_gen_dif() const
{
	const auto &gen_data = model.genetic_data;

	vector < vector <unsigned int> > gen_dif;
	
	auto nobs = gen_data.obs.size();
	
	vector <unsigned int> root_num;
	
	for(auto m = 0u; m < nobs; m++){
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
						if(gen_dif[ob1][ob2] != UNSET) emsg("Should not be unset3");
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
			if(dif[i] != UNSET) emsg("Should not be gen obs2");
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
	
	if(model.mode == SIM || !gen_data.on) return 0;
	
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
double State::likelihood_genetic_obs(const vector < vector <unsigned int> > &gen_val) const 
{
	const auto &gen_data = model.genetic_data;

	if(model.mode == SIM || !gen_data.on || !trans_tree_on) return 0;
	
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


/// Updates the markov likelihood
void State::update_like_markov(unsigned int m, unsigned int ti, double dlike_markov, unsigned int p)
{
	if(isinf(dlike_markov)){  // This can be infinite if OUTSIDE_INF probability is zero	
		vector <unsigned int> list; list.push_back(ti);
		species[p].likelihood_markov(m,list,like.markov);
	}
	else{
		like.markov += dlike_markov;
					
		species[p].Li_markov[m][ti] += dlike_markov;	
	}
}
							
							
/// Performs a global transmission tree proposal
void State::trans_tree_proposal(const BurnInfo &burn_info, unsigned int &nac, unsigned int &ntr, bool gen_obs_max)
{
	if(gen_obs_max && !genetic_value.maximise_gen_obs) return;
	
	auto pl = false;
	
	// Makes a global list of all events
	if(pl) cout << "Start trans tree proposal" << endl;
	
	const auto &precalc = param_val.precalc;

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

	vector < vector <PopInd> > pop_ind;
	pop_ind.resize(model.npop);
	
	auto dt = model.details.dt;
	
	auto egl = 0u;
	for(auto ti = 0u; ti < T; ti++){
		if(pl) cout << ti << "ti" << endl;
	
		// Updates populations
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
					if(pop_ind[po][ind].i != i) emsg("Should be i");
					
					auto N = pop_ind[po].size();
					if(ind+1 < N){
						auto &pi = pop_ind[po][ind];
						pi = pop_ind[po][N-1];
						ind_pop[p][pi.i].pop_ref[pi.index].index = ind;
					}
					pop_ind[po].pop_back();
				}
				
				ip.pop_ref.clear();
			}
			
			auto c_after = ev.c_after;
			if(c_after != UNSET){                        // Adds individual to population
				for(auto po : sp.comp_gl[c_after].pop_ref_simp){
					
					PopRef pr; pr.po = po; pr.index = pop_ind[po].size();
					PopInd pi; pi.i = i; pi.index = ip.pop_ref.size();
					
					ip.pop_ref.push_back(pr);
					pop_ind[po].push_back(pi);
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
			
			//check_w();
			
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
			
					const auto &lin = eq.lin;
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
								auto val = eq.calculate_no_pop(ti,precalc);	
								probfi += log(val*dt/mev.value_t[ti]);
							}
							else{
								auto pr = pop_ref[k_from];
								auto val = eq.calculate_pop_grad(k_from,ti,popcombw_value,precalc);
								probfi += log(val*popnum[pr]*dt/mev.value_t[ti]);
							}
						}
						
						if(k_from != UNSET){
							const auto &ind_from = species[ii_from.p].individual[ii_from.i];
			
							auto pr = ii_from.po;
							
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
							auto ss = eq.setup_source_sampler(ti,popnum,param_val,popcombw_value);

							k_prop = ss.sample_inf_source();
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
							
								auto p = po.p;
								
								ii_prop.p = p;
								ii_prop.pref = k_prop;
								ii_prop.po = pr;
								
								const auto &list = pop_ind[pr];
								if(po.ind_variation){
									const auto &ssp = species[p];
									
									vector <Poss> poss;
									for(const auto &pi : list){
										auto i = pi.i;
										const auto &ind = ssp.individual[i];
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
									ii_prop.i = list[kk].i;
									ii_prop.w = 1;
								
									probif += log(1.0/list.size());
								}
							}
						}
					}
					
					auto dlike_markov = probif - probfi;
					
					if(probif != -LARGE && (ii_from.i != ii_prop.i || ii_from.p != ii_prop.p || ii_from.pref != ii_prop.pref)){
						auto n = ev.inf_node_ref;
						
						auto type = REGRAFT_NODE; if(gen_obs_max) type = REGRAFT_GEN_OBS_MAX;
					
						auto gc = trans_tree_change(type,n,no_alter(n),ii_prop);
		
						auto al = 0.0;
						if(gc.type != TRANS_TREE_FAIL){
							if(gen_obs_max){
								if(gc.regraft_type != REGRAFT_NOT_NEEDED && gc.dlike_genetic_obs > 0 && dlike_markov > LI_WRONG) al = 1; 
							}
							else{
								al = exp(burn_info.markov*dlike_markov 
										+ burn_info.genetic_process*gc.dlike_genetic_process
										+ burn_info.genetic_obs*gc.dlike_genetic_obs
										+ probfi-probif+gc.probfi-gc.probif);
							}
						}
				
						if(pl) cout << al << " al" << endl;
						
						ntr++;
						if(ran() < al){
							if(pl) cout << "accept" << endl;
							nac++;
							gen_change_update(gc);
							update_like_markov(m,ti,dlike_markov,p);
						}
						
						if(pl){
							cout << "AFTER" << endl;
							ssp.print_event(ssp.individual[i].ev);
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
		
		if(pl || slow_check) check_pop_ind(pop_ind,ind_pop);
	}
}


/// Moves from one node location up the tree and back down again
NodeMove State::node_move(unsigned n, int e, double t, unsigned int num_max) const 
{
	NodeMove nm;
	nm.success = false;
	 
	auto &inf_node = genetic_value.inf_node;
	
	auto num = 0u;
	do{ // Goes up the tree
		if(n > CODE) break;
			
		e--;
		if(e < 0){ // Goes up to the next node
			const auto &fr = inf_node[n].from;
			n = fr.node;
			e = fr.index;
			
			num++;
		
			if(ran() < 0.5) break;  // Determines if goes down the next branch
		}						
		else{
			const auto &ie = inf_node[n].inf_ev[e];
			if(ie.type == INFECT_OTHER){
				num++;
				if(ran() < 0.5){
					n = ie.index; e = -1;
					break;
				}
			}
		}
	}while(num < num_max);
	
	auto fl = false;
	
	if(n > CODE){  // Picks a new origin to go down
		const auto &inf_origin = genetic_value.inf_origin;
	
		auto k = (unsigned int)(ran()*inf_origin.size());
		if((int)k == e) fl = true; // Stops going down the same origin
		else{		
			n = inf_origin[k].node;
			e = -1;
		}
	}
	
	if(t < inf_node[n].ti_start) fl = true;
		
	if(num < num_max && !fl){
		do{ // Goes down the tree
			e++;
			const auto &no2 = inf_node[n];
			if((unsigned int)e == no2.inf_ev.size()){
				if(t > no2.ti_rec) fl = true;
				break;
			}
			else{
				const auto &ie = inf_node[n].inf_ev[e];
			
				if(ie.tdiv >= t) break;
			
				if(ie.type == INFECT_OTHER){
					num++;
					if(ran() < 0.5){ n = ie.index; e = -1;}	
				}
			}
		}while(num < num_max);
		
		if(num < num_max && !fl){
			nm.success = true;
			nm.n = n;
			nm.e = e;
			nm.num = num;
		}
	}	

	return nm;
}


/// Performs a local transmission tree proposal
void State::local_trans_tree_proposal(const BurnInfo &burn_info, unsigned int &nac, unsigned int &ntr, bool gen_obs_max)
{
	auto pl = false;
	
	if(pl) check("local transtree update start");
	
	if(gen_obs_max && !genetic_value.maximise_gen_obs) return;
		
	auto &inf_node = genetic_value.inf_node;

	auto nnode = inf_node.size();
	for(auto k = 0u; k < nnode; k++){
		auto nod = (unsigned int)(ran()*nnode);
		
		const auto &no = inf_node[nod]; 
		
		auto n = no.from.node;
		if(n != ENTER_INF && n != OUTSIDE_INF){
			auto t = no.tdiv_start;
			auto e = no.from.index;
			auto nm = node_move(n,e,t,NODE_MAX_DIST);
			
			if(nm.success && testing){
				if(slow_check){  // Used for testing reverse transitions 
					auto loop = 0u, loopmax = 100000u;
					while(loop < loopmax){
						auto nm2 = node_move(nm.n,nm.e,t,NODE_MAX_DIST);
				
						if(nm2.n == n && nm2.e == e){
							if(nm2.num != nm.num) emsg("num rev wrong");
							break;
						}
					
						loop++;
					}
					if(loop == loopmax){
						cout << nod << " " << n << " " << e << " -> "<< nm.n << " " << t << " node move" << endl;					
						emsg("cannot get back");
					}
				}
							
				const auto &ind = species[no.p].individual[no.i];
				const auto &ev = ind.ev[no.e];
				auto tr_gl = ev.tr_gl;
				
				const auto &no_new = inf_node[nm.n];		
				const auto &ev_from = species[no_new.p].individual[no_new.i].ev;
				
				auto ii_prop = get_ind_inf_from(t,no.p,tr_gl,no_new.p,no_new.i,ev_from);

				if(ii_prop.pref != UNSET){
					// Works out change in like_markov as a result of event source changing
					auto &ssp = species[no.p];
					ssp.back.clear();
					
					auto dlike_markov = ssp.change_trans_tree_Li_markov_iif(tr_gl,ev.ind_inf_from,t,ind,-1);
					dlike_markov += ssp.change_trans_tree_Li_markov_iif(tr_gl,ii_prop,t,ind,1);
			
					auto type = REGRAFT_NODE; if(gen_obs_max) type = REGRAFT_GEN_OBS_MAX;
					
					auto gc = trans_tree_change(type,nod,no_alter(nod),ii_prop);

					auto al = 0.0;
					if(gc.type != TRANS_TREE_FAIL){
						if(gen_obs_max){
							if(gc.regraft_type != REGRAFT_NOT_NEEDED && gc.dlike_genetic_obs > 0) al = 1; 
						}
						else{
							al = exp(burn_info.markov*dlike_markov 
								+ burn_info.genetic_process*gc.dlike_genetic_process
								+ burn_info.genetic_obs*gc.dlike_genetic_obs
								+ gc.probfi-gc.probif);
						}
					}
				
					if(pl) cout << al << " al" << endl;
					
					ntr++;
					if(ran() < al){
						if(pl) cout << "accept" << endl;
						nac++;
						gen_change_update(gc);
						like.markov += dlike_markov;
					}
					else{
						if(pl) cout << "reject" << endl;
						ssp.restore_back2();
					}
					ssp.back.clear();
						
					if(pl) check("local transtree update");
				}
			}
		}
	}
}


/// Sets inf_ev_from in event at time t based on infection coming from event_from
// If not possible then returns false
bool State::set_ind_inf_from(double t, unsigned int p, vector <Event> &event, unsigned int p_from, unsigned int i_from, const vector <Event> &event_from)
{
	auto e = 0u; while(e < event.size() && event[e].tdiv < t) e++;
	if(e == event.size() || t != event[e].tdiv){
		emsg("Problem iif");
	}
	
	auto t_do = round_down(t);
	
	auto c_from = UNSET;
	auto e2 = 0u;
	while(e2 < event_from.size() && event_from[e2].tdiv <= t_do){ c_from = event_from[e2].c_after; e2++;}
	
	if(c_from == UNSET) emsg("Should not be unset4");
	
	auto &ev = event[e];
	const auto &inf_c = model.inf_cause[p][ev.tr_gl][p_from][c_from];
	auto pref = inf_c.pref;
	if(pref == UNSET) return false;
	
	auto &iif = event[e].ind_inf_from;
	iif.p = p_from;
	iif.i = i_from;
	iif.pref = pref;
	iif.po = inf_c.po;
	iif.w = get_w_from_indinffrom(iif);
	return true;
}


/// Gets inf_inf_from for a particular transition
// If not possible then returns false
IndInfFrom State::get_ind_inf_from(double t, unsigned int p, unsigned int tr_gl, unsigned int p_from, unsigned int i_from, const vector <Event> &event_from) const 
{
	auto t_do = round_down(t);
	
	auto c_from = UNSET;
	auto e2 = 0u;
	while(e2 < event_from.size() && event_from[e2].tdiv <= t_do){ c_from = event_from[e2].c_after; e2++;}
	if(c_from == UNSET) emsg("Should not be unset4");

	const auto &inf_c = model.inf_cause[p][tr_gl][p_from][c_from];
	
	IndInfFrom iif;
	iif.p = p_from;
	iif.i = i_from;
	iif.pref = inf_c.pref;
	if(iif.pref != UNSET){
		iif.po = inf_c.po;
		iif.w = get_w_from_indinffrom(iif);
	}
	
	return iif;
}


/// Sets inf_ev_from in event at time t based on infection coming from event_from
// If not possible then returns false
void State::set_ind_inf_from_outside(double t, vector <Event> &event)
{
	auto e = 0u; while(e < event.size() && event[e].tdiv < t) e++;
	if(e == event.size() || t != event[e].tdiv) emsg("Problem iif outside");
	
	auto &iif = event[e].ind_inf_from;
	iif.p = OUTSIDE_INF;
	iif.i = UNSET;
	iif.pref = UNSET;
	iif.po = UNSET;
	iif.w = 1;
}


/// Adds a change to an iif
void State::add_ind_inf_change(unsigned int p_ev, unsigned int i_ev, unsigned int e_ev, unsigned int c_bef, unsigned int c_aft, unsigned int p_aft, unsigned int i_aft, vector <IIFChange> &iif_change, bool &ill_fl) const
{	
	const auto &ssp = species[p_ev];
	const auto &indd = ssp.individual[i_ev];
	const auto &ev = indd.ev[e_ev];
	
	const auto &iif_bef = ev.ind_inf_from;
	auto ti = get_ti(ev.tdiv);
		
	IIFChange iif_ch;                // Used to store a change in source
	iif_ch.p_ev = p_ev; 
	iif_ch.i_ev = i_ev; 
	iif_ch.e_ev = e_ev; 
	iif_ch.ti = ti;
	
	auto &iif_aft = iif_ch.ind_inf_from;
		
	if(c_bef == UNSET){
		c_bef = species[iif_bef.p].get_ind_c(iif_bef.i,ti);
		c_aft = species[p_aft].get_ind_c(i_aft,ti);

		if(c_bef == UNSET) emsg(" cannot be");
		if(c_aft == UNSET){
			ill_fl = true;
			return;
		}
		
		if(c_bef == c_aft && iif_bef.p == p_aft && !model.ind_pop_variation){
			if(i_aft != iif_bef.i){
				iif_aft = iif_bef;
				iif_aft.i = i_aft;
				iif_ch.dLi_markov = 0;
				
				iif_change.push_back(iif_ch);		
			}
			return;
		}
	}
	
	auto tr_gl = ev.tr_gl;
	const auto &sp = model.species[p_ev];
	const auto &me = sp.markov_eqn[sp.tra_gl[tr_gl].markov_eqn_ref];

	auto inffac = 1.0; if(me.ind_variation) inffac = ssp.get_indfac(indd,me);
		
	if(true){ // TURN OFF
		const auto &inf_c_bef = model.inf_cause[p_ev][tr_gl][iif_bef.p][c_bef];
		if(inf_c_bef.pref != iif_bef.pref) emsg("pref not agree");
		if(inf_c_bef.po != iif_bef.po) emsg("po not agree");
		if(inf_c_bef.eq_ref != me.eqn_ref) emsg("eqn not right");
	}
	
	const auto &eq = model.eqn[me.eqn_ref];
	
	const auto &precalc = param_val.precalc;
	auto dt = model.details.dt;
	
	auto va_bef = eq.calculate_pop_grad(iif_bef.pref,ti,popcombw_value,precalc);
	iif_ch.dLi_markov = -log(iif_bef.w*dt*log_thresh(va_bef)*inffac);
	
	// Looks at after
	const auto &inf_c_aft = model.inf_cause[p_ev][tr_gl][p_aft][c_aft];
	auto pref = inf_c_aft.pref;
	if(pref == UNSET) ill_fl = true;
	else{
		iif_aft.p = p_aft;
		iif_aft.i = i_aft;
		iif_aft.pref = pref;
		iif_aft.po = inf_c_aft.po;
		iif_aft.w = get_w_from_indinffrom(iif_aft);

		auto va_aft = eq.calculate_pop_grad(pref,ti,popcombw_value,precalc);
		iif_ch.dLi_markov += log(iif_aft.w*dt*log_thresh(va_aft)*inffac);

		iif_change.push_back(iif_ch);							
	}
}
					
					
/// Accept changes to iif
void State::accept_iif_change(const vector <IIFChange> &iif_change)
{
	for(const auto &iif_ch : iif_change){
		auto p_ev = iif_ch.p_ev;
		auto &ssp = species[p_ev];
		auto &ev = ssp.individual[iif_ch.i_ev].ev[iif_ch.e_ev];
		ev.ind_inf_from = iif_ch.ind_inf_from;
		
		auto m = model.species[p_ev].tra_gl[ev.tr_gl].markov_eqn_ref;
		ssp.Li_markov[m][iif_ch.ti] += iif_ch.dLi_markov;
	}
}
					

/// The overall change in likelihood as a result of 
double State::iif_change_dmarkov(const vector <IIFChange> &iif_change) const
{
	auto sum = 0.0;
	for(const auto &iif_ch : iif_change) sum += iif_ch.dLi_markov;

	return sum;
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
			
			// Checks to see if there is an observation on B before branch to A
			auto fl = false;
			for(auto e = 0u; e < e_index; e++){
				if(inf_evB[e].type == GENETIC_OBS) fl = true;
			}	
			
			// Checks to see if there is an onward infection which cannot be accounted for due to finite step size
			auto t_inf_min = in_A.ti_start;
			for(auto e = e_index+1; e < inf_evB.size(); e++){
				if(inf_evB[e].tdiv < t_inf_min){
					if(inf_evB[e].type == INFECT_OTHER){
						fl = true;
					}
				}
				else break;
			}
			
			// Checks to see if B is entering the system
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
				
				// Makes sure that neither of the events are directly observed
				if(!ev_A_new[in_A.e].observed && !ev_B_new[in_B.e].observed){
			
					auto fl = false;
					
					// Checks if any events in B in critical classification before t_A
					auto e_B_after = get_event_after(in_B.e,ev_B_new);
					if(e_B_after != UNSET && ev_B_new[e_B_after].tdiv < t_A) fl = true;
					 
					// Checks if any events in A in critical classification between t_B and t_A
					auto e_A_before = get_event_before(in_A.e,ev_A_new);
					if(e_A_before != UNSET && ev_A_new[e_A_before].tdiv > t_B) fl = true;
					
					// Check for unaccounted observations on A
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
							
								// Sets new iif for A
								if(set_ind_inf_from(t_B,p_A,ev_A_new,p_fr,i_fr,species[p_fr].individual[i_fr].ev) == false){
									fl = true;
								}
							}
						}
						
						// Sets new iif for B
						if(fl == false){
							if(set_ind_inf_from(t_A,p_B,ev_B_new,p_A,i_A,ev_A_new) == false) fl = true;
						}
						
						if(fl == true){
							nfa++;
						}
						else{
							vector <IIFChange> iif_change;
							
							auto ill_fl = false;
								
							// Looks at switching the source for other infections on B's time line
							if(e_index > 0){
								const auto &ev_bef = species[p_B].individual[i_B].ev;
								const auto &ev_aft = ev_A_new;
								auto c_bef = UNSET;
								auto c_aft = UNSET;
								auto e_bef = 0u, e_aft = 0u;
								
								for(auto e = 0u; e < e_index; e++){
									const auto &ie = in_B.inf_ev[e];
									auto ti = get_ti(ie.tdiv);
									
									while(e_bef < ev_bef.size() && ev_bef[e_bef].tdiv < ti){
										c_bef = ev_bef[e_bef].c_after;
										e_bef++;
									}
									
									while(e_aft < ev_aft.size() && ev_aft[e_aft].tdiv < ti){
										c_aft = ev_aft[e_aft].c_after;
										e_aft++;
									}
									
									if(p_B != p_A || c_bef != c_aft || model.ind_pop_variation){  // The source has switched	
										const auto &in_inf = inf_node[ie.index];
									
										add_ind_inf_change(in_inf.p,in_inf.i,in_inf.e,c_bef,c_aft,p_A,i_A,iif_change,ill_fl);
									}
								}
							}
							
							if(ill_fl == false){
								back_init();
								auto uii = update_ind(p_A,i_A,ev_A_new,UP_MULTI,false);
								auto &like_ch = uii.like_ch;  
								
								auto uii2 = update_ind(p_B,i_B,ev_B_new,UP_MULTI,false);
								auto &li_ch = uii2.like_ch;
								
								add_on_like(li_ch,like_ch);
									
								like_ch.markov += iif_change_dmarkov(iif_change);
								
								auto al = calc_al(like_ch,0,burn_info);
								
								if(pl) print_like(li_ch);
								
								if(pl) cout << al << " al" << endl;
								if(ran() < al){
									if(pl) cout << "accept" << endl;
									accept(like_ch);
							
									// Updates quantities for inf_node
									InfNode in_A_new;
									in_A_new.tdiv_start = t_B;
									in_A_new.tdiv_rec = in_A.tdiv_rec;
									in_A_new.ti_start = inf_node_ti_from_tdiv(t_B);
									in_A_new.ti_rec = inf_node_ti_from_tdiv(in_A_new.tdiv_rec);
									in_A_new.p = p_A;
									in_A_new.i = i_A;
									in_A_new.e = get_e_from_time(t_B,ev_A_new);
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
									in_B_new.ti_start = inf_node_ti_from_tdiv(t_A);
									in_B_new.ti_rec = inf_node_ti_from_tdiv(in_B_new.tdiv_rec);
									in_B_new.p = p_B;
									in_B_new.i = i_B;
									in_B_new.e = get_e_from_time(t_A,ev_B_new);
									in_B_new.from.node = k_B;
									in_B_new.from.index = e_index;
									for(auto e = e_index+1; e < in_B.inf_ev.size(); e++){
										in_B_new.inf_ev.push_back(in_B.inf_ev[e]);
									}
									
									in_A = in_A_new;
									in_B = in_B_new;
			
									accept_iif_change(iif_change);
								
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
			}
			if(pl) check(" transtree update2");
		}
	}
}


/// Gets a branch for trans_tree_swap_inf_ev_proposal
double State::get_branch_ev(vector <BranchEv> &branch, unsigned int n, int e, unsigned num, double t) const
{
	auto &inf_node = genetic_value.inf_node;
	
	// Initial state 
	BranchEv be; be.type = START_SI; be.n = n; be.e = e; be.mut_sum = 0u;
	if(e == -1) be.tdiv = inf_node[n].tdiv_start;
	else be.tdiv = inf_node[n].inf_ev[e].tdiv;
	branch.push_back(be);
	
	auto num_branch = 0u;
	auto num_obs = 0u;
	auto mut_sum = 0u;
	do{     // Goes down branch 
		const auto &br = branch[branch.size()-1];
		auto n = UNSET;
		int e = UNSET;
		switch(br.type){
		case START_SI: 
			n = br.n; e = br.e+1;
			break;
			
		case SHIFT_IND_SI:
			n = inf_node[br.n].inf_ev[br.e].index;
			e = 0;
			break;
		
		case BRANCH_SI: case GENETIC_OBS_SI:
			emsg("SI problem");
			break;
		}
		
		auto num_branch_st = num_branch;
		
		const auto &inf_ev = inf_node[n].inf_ev;
		while(e < (int)inf_ev.size()){
			const auto &ie = inf_node[n].inf_ev[e];
			mut_sum += ie.mut_num;
			BranchEv be; be.n = n; be.e = e; be.tdiv = ie.tdiv; be.mut_sum = mut_sum;
			switch(inf_ev[e].type){
			case INFECT_OTHER: be.type = BRANCH_SI; num_branch++;	break;
			case GENETIC_OBS: be.type = GENETIC_OBS_SI; num_obs++; break;
			}
			branch.push_back(be);
			e++;
		}
		
		if(num != UNSET){
			if(num_branch > num) return UNSET;
		
			if(num_branch == num) break;
		}
		else{
			if(inf_node[n].ti_rec >= t) return t;
		}
		
		if(num_branch == num_branch_st) return UNSET;
		
		auto &br_last = branch[branch.size()-1];
		if(br_last.type != BRANCH_SI) return UNSET;
		br_last.type = SHIFT_IND_SI;
		num_branch--;
	}while(true);
	
	return branch[branch.size()-1].tdiv;
}


/// Gets probability of generating mutation numbers
double State::branch_sample_prob(vector <BranchEv> &branch, bool resample, string &warn) const
{
	auto prob = 0.0;
	
	auto kmax = branch.size();
	
	auto k = 1u;
	while(k < kmax){
		auto kk = k; while(kk < kmax && branch[kk].type != BRANCH_SI) kk++;

		if(kk < kmax){          // Samples between branchs
			vector <unsigned int> mnum; 
						
			auto N = branch[kk].mut_sum-branch[k-1].mut_sum;	
			if(kk > k){           // Multinomial sampling	
				if(N > 0){
					auto t_tot = branch[kk].tdiv-branch[k-1].tdiv;
				
					vector <double> frac;
					for(auto kkk = k; kkk <= kk; kkk++){
						const auto &br = branch[kkk];
						const auto &br_last = branch[kkk-1];
						frac.push_back((br.tdiv-br_last.tdiv)/t_tot);
						if(!resample) mnum.push_back(br.mut_sum-br_last.mut_sum);
					}
					
					if(resample){
						mnum = multinomial_sample(N,frac);
					}
					
					prob += multinomial_probability(N,mnum,frac);
					
					/*
					for(auto va : frac) cout << va << ","; cout << " frac";
					for(auto va : mnum) cout << va << ","; cout << " mnum";
					cout <<probfi << " " << N << " N";
					emsg("multi");
					*/
				}
				else{
					if(resample){
						for(auto kkk = k; kkk < kk; kkk++) mnum.push_back(0);
					}
				}
			}
			else{
				if(resample) mnum.push_back(N);
			}
			
			if(resample){
				auto j = 0u;
				for(auto kkk = k; kkk < kk; kkk++){
					branch[kkk].mut_sum = branch[kkk-1].mut_sum + mnum[j];
					j++;
				}
			}
		}
		else{    // Samples mutation number for any observations at the end
			auto mut_rate = genetic_value.mut_rate*model.details.dt;
			for(auto kkk = k; kkk < kmax; kkk++){
				auto &be = branch[kkk];
				if(be.type == BRANCH_SI) emsg("should not be branch");
			
				const auto &be_last = branch[kkk-1];
			
				auto lam = mut_rate*(be.tdiv-be_last.tdiv);
				if(resample){
					be.mut_sum = be_last.mut_sum + poisson_sample(lam,warn);
				}
			
				prob += poisson_probability(be.mut_sum-be_last.mut_sum,lam);
			}
		}
		k = kk+1;
	}
	
	return prob;
}


/// Prints information about a branch
void State::print_branch(string te, const vector <BranchEv> &branch) const
{
	cout << endl << te << ":   " << endl;
	for(const auto &ev : branch){
		switch(ev.type){
		case START_SI: cout << "START_SI"; break;
		case BRANCH_SI: cout << "BRANCH_SI"; break;
		case SHIFT_IND_SI: cout << "SHIFT_IND_SI"; break;
		case GENETIC_OBS_SI: cout << "GENETIC_OBS_SI"; break;
		}
		
		cout << " n=" << ev.n << " e=" << ev.e << " tdiv=" << ev.tdiv << " mut_sum=" << ev.mut_sum << endl;
	}		
}


/// Constructs a sew inf_ev based on a branch 
vector <NewInfEv> State::construct_new_inf_ev(const vector <BranchEv> &branch, vector <IIFChange> &iif_change, bool &ill_fl) const 
{
	vector <NewInfEv> list;
	
	const auto &inf_node = genetic_value.inf_node;
	
	auto kmax = branch.size(); 
	
	auto k = 1u;
	while(k <= kmax){
		auto kk = k; while(kk < kmax && branch[kk].type != SHIFT_IND_SI) kk++;

		NewInfEv nie;
		auto &inf_ev = nie.inf_ev;
		
		const auto &br = branch[k-1];
		auto &n = nie.n;
		
		auto estart = 0u;
		
		switch(br.type){
		case START_SI:
			{		
				n = br.n;
				auto e = br.e;
				if(e != -1){ // Copies events before infection
					const auto &inf_ev2 = inf_node[n].inf_ev;
					for(auto ee = 0u; ee <= (unsigned int)e; ee++){
						inf_ev.push_back(inf_ev2[ee]);
					}
					estart = inf_ev.size();
				}
			}
			break;
			
		case SHIFT_IND_SI:
			n = inf_node[br.n].inf_ev[br.e].index;
			break;
			
		default:	
			emsg("Option error2");
			break;
		}
			
		for(auto kkk = k; kkk <= kk; kkk++){
			if(kkk < kmax){
				const auto &br = branch[kkk];
				const auto &br_last = branch[kkk-1];
				auto ie = inf_node[br.n].inf_ev[br.e];
				ie.mut_num = br.mut_sum - br_last.mut_sum;
			
				inf_ev.push_back(ie);
			}
		}		

		// Works out changes in iif
		const auto &no = inf_node[n];
		for(auto e = estart; e < inf_ev.size(); e++){
			const auto &iev = inf_ev[e];
			if(iev.type == INFECT_OTHER){
				const auto &in_inf = inf_node[iev.index];
				add_ind_inf_change(in_inf.p,in_inf.i,in_inf.e,UNSET,UNSET,no.p,no.i,iif_change,ill_fl);
			}
		}

		list.push_back(nie);
		
		k = kk+1;		
	}
	
	return list;
}


/// Accepts new inf_ev 
void State::accept_new_inf_ev(const vector <NewInfEv> &list)
{
	auto &inf_node = genetic_value.inf_node;
	auto &obs_node_ref = genetic_value.obs_node_ref;
	
	for(const auto &nie : list){
		auto n = nie.n;
		auto &inf_ev = inf_node[n].inf_ev;
		inf_ev = nie.inf_ev;
		
		for(auto e = 0u; e < inf_ev.size(); e++){		
			const auto &ev = inf_ev[e];
			
			switch(ev.type){
			case INFECT_OTHER:
				{
					auto &no = inf_node[ev.index];
					auto &from = no.from;
					from.node = n;
					from.index = e;
				}
				break;
				
			case GENETIC_OBS:
				{
					auto m = ev.index;	
					auto &onr = obs_node_ref[m];
					onr.index = e;
				}
				break;
			}
		}
		
	}
}


/// Prints any create infection events
void State::print_new_inf_ev(string te, const vector <NewInfEv> &list) const
{
	cout << endl << endl << "New infection events: " << te << endl;
	
	const auto &inf_node = genetic_value.inf_node;
	for(const auto &nie : list){
		const auto &no = inf_node[nie.n];
		cout << "Node " << nie.n << " " << species[no.p].individual[no.i].name << ":" << endl;
		for(const auto &ev : nie.inf_ev){
			cout << "t=" << ev.tdiv << " index=" << ev.index << " mut=" << ev.mut_num << "  ";
			switch(ev.type){
			case INFECT_OTHER: cout << "INFECT_OTHER"; break;
			case GENETIC_OBS: cout << "GENETIC_OBS"; break;
			}
			cout << endl;
		}
		cout << endl;
	}
}


/// Checks branch has been correctly specified 
void State::check_branch(const vector <BranchEv> &branch) const
{
	for(auto k = 0u; k < branch.size(); k++){
		const auto &br = branch[k];
		if(br.mut_sum == UNSET) emsg("mut_sum unset");
		if(k > 0){
			if(br.mut_sum < branch[k-1].mut_sum) emsg("mut_sum order");
		}
	}
}


/// Converts from a brahc to branch mutation info
BranchMutInfo State::get_branch_mut_info(const vector <BranchEv> &branch, int sign) const
{
	BranchMutInfo bmi;
	
	auto &inf_node = genetic_value.inf_node;
	
	for(auto e = 1u; e < branch.size(); e++){
		const auto &br = branch[e];
		switch(br.type){
		case BRANCH_SI:
			{
				AddMutSum add;
				auto n = inf_node[br.n].inf_ev[br.e].index;
				add.ob = get_gen_obs_list(n,0);
				add.mut_sum = sign*br.mut_sum;
				add.mut_sum_after = UNSET;
				bmi.add_list.push_back(add);	
			}
			break;
		
		case GENETIC_OBS_SI:
			{
				ObsMutSum obs;
				obs.ob = inf_node[br.n].inf_ev[br.e].index;
				obs.mut_sum = sign*br.mut_sum;
				obs.mut_sum_after = UNSET;
				bmi.obs_list.push_back(obs);	
			}
			break;
			
		default: break;
		}
	}
	
	return bmi;
}


/// Converts from a brahc to branch mutation info
void State::set_branch_mut_info_after(const vector <BranchEv> &branch, vector <ObsMutSum> &obs_list, vector <AddMutSum> &add_list, int sign) const
{
	auto nbr = 0u, nob = 0u;
	for(auto e = 1u; e < branch.size(); e++){
		const auto &br = branch[e];
		switch(br.type){
		case BRANCH_SI:
			add_list[nbr].mut_sum_after = sign*br.mut_sum;
			nbr++;
			break;
		
		case GENETIC_OBS_SI:
			obs_list[nob].mut_sum_after = sign*br.mut_sum;
			nob++;
			break;
			
		default: break;
		}
	}
		
	if(nbr != add_list.size()) emsg("nbr prob");
	if(nob != obs_list.size()) emsg("nob prob");
}


/// Prints branch mutation information
void State::print_branch_mut_info(string te, const BranchMutInfo &bmi) const
{
	cout << endl << endl << "Branch mut info " << te << ":" << endl;
	for(const auto &obs : bmi.obs_list){
		cout << "OBS mut_sum=" << obs.mut_sum << " " << "mut_sum_after=" << obs.mut_sum_after << " ";
		cout << "  m=" << obs.ob << endl;
	}
	
	for(const auto &add : bmi.add_list){
		cout << "ADD mut_sum=" << add.mut_sum << " " << "mut_sum_after=" << add.mut_sum_after << " ";
		cout << "  m=";
		for(auto m : add.ob) cout << m << ",";
		cout << endl;
	}
}


/// Works out the change in genetic observation likelihood as well as 
double State::branch_genetic_process(const vector <BranchEv> &branch) const
{
	auto mut_rate = genetic_value.mut_rate*model.details.dt;
	
	auto sum = 0.0;
	for(auto i = 1u; i < branch.size(); i++){
		const auto &br = branch[i];
		const auto &br_last = branch[i-1];
		auto lam = mut_rate*(br.tdiv-br_last.tdiv);
		
		sum += poisson_probability(br.mut_sum-br_last.mut_sum,lam);
	}
	
	return sum;
}

		
/// Works out the change in genetic observation likelihood as well as gen_dif
double State::branch_change_genetic_obs(const vector <BranchEv> &branch1, const vector <BranchEv> &branch2, const vector <BranchEv> &branch1_new, const vector <BranchEv> &branch2_new, vector <GenDifCha> &gen_dif_cha) const
{
	/*
	print_branch("branch1",branch1);
	print_branch("branch2",branch2);
	print_branch("branch1_new",branch1_new);
	print_branch("branch2_new",branch2_new);
	*/
	
	auto mut_info1 = get_branch_mut_info(branch1,1);
	auto mut_info2 = get_branch_mut_info(branch2,-1);
	
	set_branch_mut_info_after(branch1_new,mut_info1.obs_list,mut_info2.add_list,1);
	set_branch_mut_info_after(branch2_new,mut_info2.obs_list,mut_info1.add_list,-1);
	
	/*
	print_branch_mut_info("mut_info1",mut_info1);
	print_branch_mut_info("mut_info2",mut_info2);
	*/
	
	// Combines information from both branchs into a single list
	auto &obs_list = mut_info1.obs_list;
	for(const auto &ob : mut_info2.obs_list) obs_list.push_back(ob);

	auto &add_list = mut_info1.add_list;
	for(const auto &ad : mut_info2.add_list) add_list.push_back(ad);
	
	// Accounts for all the observation not included above
	{
		const auto &obs_node_ref = genetic_value.obs_node_ref;
		auto N = obs_node_ref.size();
		
		vector <bool> map(N,false);
		for(const auto &obs : obs_list) map[obs.ob] = true;
		
		for(const auto &add : add_list){
			for(auto m : add.ob) map[m] = true;
		}
		
		AddMutSum add;
		auto &ob = add.ob;
		
		for(auto m = 0u; m < N; m++){
			if(map[m] == false && obs_node_ref[m].node != UNSET) ob.push_back(m);
		}
		
		if(ob.size() > 0){
			add.mut_sum = 0;
			add.mut_sum_after = 0;
			
			add_list.push_back(add);
		}
	}

	const auto &gen_dif = genetic_value.gen_dif;
	const auto &gen_val = model.genetic_data.gen_dif;

	auto sum = 0.0;
	int d, d_update;
	// Accounts for differences between observations and added branchs
	for(const auto &add : add_list){
		if(add.ob.size() > 0){
			auto mut_sum = add.mut_sum;
			auto mut_sum_after = add.mut_sum_after;
			
			GenDifCha gdc;
			gdc.ob2 = add.ob;
			gdc.sign = 1;
			for(const auto &obs : obs_list){
				auto ch = abs(obs.mut_sum_after-mut_sum_after) - abs(obs.mut_sum-mut_sum);
				if(ch != 0){
					auto m = obs.ob;
					
					const auto &gen_dif2 = gen_dif[m];
					const auto &gen_val2 = gen_val[m];
					for(auto m2 : add.ob){
						d = int(gen_dif2[m2]) - int(gen_val2[m2]);
						d_update = d + ch;
						sum += d_update*d_update - d*d;
					}
					
					gdc.ob1 = m;
					gdc.ch = ch; 
					gen_dif_cha.push_back(gdc);
				}
			}
		}
	}
	
	// Accounts for difference between observations
	auto N = obs_list.size();
	if(N > 1){
		for(auto j = 0u; j < N; j++){
			const auto &ob1 = obs_list[j];
			auto mut_sum = ob1.mut_sum;
			auto mut_sum_after = ob1.mut_sum_after;
			auto m = ob1.ob;
			const auto &gen_dif2 = gen_dif[m];
			const auto &gen_val2 = gen_val[m];
			
			GenDifCha gdc;
			gdc.ob2.push_back(m);
			gdc.sign = 1;
			for(auto i = j+1; i < N; i++){
				const auto &ob2 = obs_list[i];
				auto ch = abs(ob2.mut_sum_after-mut_sum_after) - abs(ob2.mut_sum-mut_sum);
				if(ch != 0){
					auto m2 = ob2.ob;
					d = int(gen_dif2[m2]) - int(gen_val2[m2]);
					d_update = d + ch;
					sum += d_update*d_update - d*d;
				
					gdc.ob1 = m2;
					gdc.ch = ch; 
					gen_dif_cha.push_back(gdc);
				}
			}
		}
	}
	
	return -2*sum*model.genetic_data.obs_scale;
}


/// Used to order branch events
bool BranchEv_ord (const BranchEv &be1, const BranchEv &be2)  
 { return (be1.tdiv < be2.tdiv); };  
 
 
/// Returns the mean branch number
double State::mean_branch_number() const 
{
	const auto &inf_node = genetic_value.inf_node;
	auto sum = 0.0;
	for(const auto &in : inf_node){
		auto nbr = 0u;
		for(const auto &ev : in.inf_ev){
			if(ev.type == INFECT_OTHER) nbr++;	
		}	
		sum += nbr;
	}
	sum /= inf_node.size();
	if(sum < 1) sum = 1;
	
	return sum;
}

 
/// Performs swapping of infected individuals prior to an indiv dual becoming infected
// Here we look to swap infection events from infected to infector
void State::trans_tree_swap_inf_ev_proposal(const BurnInfo &burn_info, vector <PropAcInfo> &prop_acc, bool use_time)
{
	auto pl = false;

	auto &inf_node = genetic_value.inf_node;
	if(inf_node.size() == 0) return;
	
	auto loop_max = (unsigned int)(0.3*inf_node.size());
	if(loop_max < 1) loop_max = 1;
	
	auto num_mean = 0.0;
	if(!use_time) num_mean = mean_branch_number(); 
	
	string warn;
	for(auto loop = 0u; loop < loop_max; loop++){
		auto num = UNSET;
		if(!use_time){	
			num = 1+floor_int(exp_mean_sample(num_mean,warn)); // The number of branches to select
		}        
		
		auto nsel = (unsigned int)(inf_node.size()*ran());

		const	auto &no = inf_node[nsel];
		
		const auto &no_fr = no.from;

		auto nfr = no_fr.node;
		if(nfr != OUTSIDE_INF && nfr != ENTER_INF){
			auto jfr = no_fr.index;

			vector <BranchEv> branch1, branch2;
			
			double t;
			if(use_time){       // This uses the maximum recovery time of the two individuals
				t = no.ti_rec;
				auto tt = inf_node[nfr].ti_rec;
				if(tt > t) t = tt;
				
				t = get_branch_ev(branch1,nsel,-1,UNSET,t);
				if(t != UNSET){
					t = get_branch_ev(branch2,nfr,jfr,UNSET,t);
					if(t != UNSET){
						num = 0;
						for(const auto &br : branch1){
							if(br.type == BRANCH_SI) num++;
						}
						for(const auto &br : branch2){
							if(br.type == BRANCH_SI) num++;
						}
					
						if(num == 0) t = UNSET;
					}
				}
			}
			else{
				if(ran() < 0.5){
					t = get_branch_ev(branch1,nsel,-1,num,UNSET);
					if(t != UNSET) get_branch_ev(branch2,nfr,jfr,UNSET,t);
				}
				else{
					t = get_branch_ev(branch1,nfr,jfr,num,UNSET);
					if(t != UNSET) get_branch_ev(branch2,nsel,-1,UNSET,t);
				}
			}
			
			if(t != UNSET){
				if(pl){
					cout << endl << endl << "usetime=" << use_time << " num=" << num << "START" << endl;
					print_node("node sel",nsel);
					print_node("node from",nfr);
				}
				
				if(pl){
					print_branch("branch1",branch1);
					print_branch("branch2",branch2);
				}	
				
				auto probfi = 0.0;
				probfi += branch_sample_prob(branch1,false,warn); 
				probfi += branch_sample_prob(branch2,false,warn); 
				
				// Switches BRANCH_SI from one branch to the other
				vector <BranchEv> branch1_new, branch2_new;
				for(const auto &el : branch1){
					if(el.type == BRANCH_SI) branch2_new.push_back(el);
					else branch1_new.push_back(el);
				}
				
				for(const auto &el : branch2){
					if(el.type == BRANCH_SI) branch1_new.push_back(el);
					else branch2_new.push_back(el);
				}
			
				sort(branch1_new.begin(), branch1_new.end(), BranchEv_ord);
				sort(branch2_new.begin(), branch2_new.end(), BranchEv_ord);
				
				if(pl){
					print_branch("branch1",branch1);
					print_branch("branch2",branch2);
					print_branch("branch1_new",branch1_new);
					print_branch("branch2_new",branch2_new);
				}
				
				auto probif = 0.0;
				probif += branch_sample_prob(branch1_new,true,warn); 
				probif += branch_sample_prob(branch2_new,true,warn);
				if(warn != ""){
					run_error("Problem with sampling: "+warn);
				}

				if(testing){
					check_branch(branch1_new);
					check_branch(branch2_new);
				}
	
				if(pl){
					print_branch("branch1 after",branch1_new);
					print_branch("branch2 after",branch2_new);
				}
				
				if(prop_acc.size() <= num) prop_acc.resize(num+1);
				prop_acc[num].ntr++;
				
				vector <IIFChange> iif_change;
				auto ill_fl = false;
			
				auto new_inf_ev1 = construct_new_inf_ev(branch1_new,iif_change,ill_fl);
				auto new_inf_ev2 = construct_new_inf_ev(branch2_new,iif_change,ill_fl);
				
				if(ill_fl) prop_acc[num].nfa = true;
				else{
					if(pl){
						print_new_inf_ev("ev1",new_inf_ev1);
						print_new_inf_ev("ev2",new_inf_ev2);
					}
					
					vector <GenDifCha> gen_dif_cha;
					
					Like like_ch;
					like_ch.markov = iif_change_dmarkov(iif_change);
					like_ch.genetic_obs = branch_change_genetic_obs(branch1,branch2,branch1_new,branch2_new,gen_dif_cha);
					like_ch.genetic_process = probif-probfi;
					
					if(false){
						auto dgen = branch_genetic_process(branch1_new) + branch_genetic_process(branch2_new) - 
												branch_genetic_process(branch1) - branch_genetic_process(branch2);
						cout << dgen << " " << probif-probfi <<  "dgen" << endl;
					}
					
					auto al = calc_al(like_ch,probfi-probif,burn_info);
					if(pl) cout << al << " al" << endl;
					//if(pl) print_like(li_ch);
			
					if(ran() < al){
						prop_acc[num].nac++;
						accept(like_ch);
						
						accept_gen_dif_cha(gen_dif_cha);
						accept_new_inf_ev(new_inf_ev1);
						accept_new_inf_ev(new_inf_ev2);
						accept_iif_change(iif_change);
					}					
				
					if(pl) check("swap inf_ev");
				}
			}
		}
	}
}


/// Gets the event number from a time
unsigned int State::get_e_from_time(double t, const vector <Event> &ev) const
{
	auto e = 0u; while(e <ev.size() && ev[e].tdiv != t) e++;
	if(e == ev.size()) emsg("Cannot find time");
	
	return e;
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
	na.ti_start = in.ti_start;
	na.ti_rec = in.ti_rec;
	return na;
}
	
	
/// Performs a transmission tree mutation proposal
void State::trans_tree_mut_proposal(const BurnInfo &burn_info, unsigned int &nac, unsigned int &ntr, double &si)
{
	auto pl = false;
	
	auto &inf_node = genetic_value.inf_node;
	const auto &gen_data = model.genetic_data;
	
	auto mut_rate = genetic_value.mut_rate*model.details.dt;
	
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
					auto obs_other = get_not_obs(obs_node,vector <unsigned int> (),nobs);
					
					auto sum = 0.0;
					int d, d_update;
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
			
					ntr++;
					if(ran() < al){ 
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
						
						burn_info.update_si(si,ACCEPT,ntr);
					}
					else{
						burn_info.update_si(si,REJECT,ntr); if(si < 1) si = 1;
					}
						
					if(pl) check("update");
				}
			}
		}
	}
}


/// Performs a transmission tree origin mutation proposal
void State::trans_tree_mut_origin_proposal(const BurnInfo &burn_info, unsigned int &nac, unsigned int &ntr, double &si)
{
	auto pl = false;
	
	auto seq_var = genetic_value.seq_var;
	
	auto &inf_node = genetic_value.inf_node;
	const auto &gen_data = model.genetic_data;
	
	auto mut_rate = genetic_value.mut_rate*model.details.dt;
	
	auto &gen_dif = genetic_value.gen_dif;
	const auto &gen_val = gen_data.gen_dif;
	
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
				int d, d_update;
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
	
				ntr++;
				if(ran() < al){ 
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
					
					burn_info.update_si(si,ACCEPT,ntr);
				}
				else{
					burn_info.update_si(si,REJECT,ntr); if(si < 1) si = 1;
				}
							
				if(pl) check("update");
			}
		}
	}
}


/// Initialsies sampler for normal integer distribution
void NormalIntSampler::init(double lin, double quad, double lam)
{
	auto var1 = -1/(2*(quad-VTINY));
	auto mean1 = lin*var1;
	
	auto var2 = lam;
	auto mean2 = lam;
	
	auto var = var1*var2/(var1+var2);
	auto mean = (mean1*var2 + mean2*var1)/(var1+var2);
	
	auto sigma = sqrt(var);
					
	auto f = 4u;
	int imin = floor_int(mean-f*sigma);
	int imax = 1+floor_int(mean+f*sigma);
	if(imin < 0) imin = 0;
	if(imax < 0) imax = 0;
	
	ele.clear();
	
	auto sum = 0.0;
	for(int i = imin; i <= imax; i++){
		auto pr = exp(-(i-mean)*(i-mean)/(2*var));
		sum += pr;
		
		NormalIntEle el;
		el.num = i;
		el.prob = pr;
		el.prob_sum = sum;
		
		ele.push_back(el);
	}
}


/// Samplers from distibution
unsigned int NormalIntSampler::sample(double &probif) const
{
	auto kmax = ele.size();
	
	auto z = ran()*ele[kmax-1].prob_sum;
	
	auto k = 0u; while(k < kmax && z > ele[k].prob_sum) k++;
	if(k == kmax) emsg("kmax problem");

	const auto &el = ele[k];
	
	probif += log(el.prob);
	
	return el.num;
}


/// Samplers from distibution
double NormalIntSampler::prob(unsigned int i) const
{
	int k = (int)i-ele[0].num;
	if(k < 0) k = 0;
	if(k >= (int)ele.size()) k = ele.size()-1;

	return log(ele[k].prob);
}


/// Performs a transmission tree mutation proposal
void State::trans_tree_mut_sample_proposal(const BurnInfo &burn_info, unsigned int &nac, unsigned int &ntr)
{
	auto pl = false;
	
	NormalIntSampler sampler;
	
	auto &inf_node = genetic_value.inf_node;
	const auto &gen_data = model.genetic_data;
	
	auto mut_rate = genetic_value.mut_rate*model.details.dt;
	auto seq_var = genetic_value.seq_var;

	auto &gen_dif = genetic_value.gen_dif;
	const auto &gen_val = gen_data.gen_dif;
	
	auto sca = model.genetic_data.obs_scale;
	
	// Goes through the infection nodes
	for(auto n = 0u; n < inf_node.size(); n++){
		auto &in = inf_node[n];
		
		for(auto e = 0u; e < in.inf_ev.size(); e++){
			auto &iev = in.inf_ev[e];
				
			double dtdiv;
			if(e == 0) dtdiv = iev.tdiv-in.tdiv_start;
			else dtdiv = iev.tdiv-in.inf_ev[e-1].tdiv;
			auto lam = mut_rate*dtdiv;
			
			auto obs_node = get_gen_obs_list(n,e);
			auto nobs = genetic_value.obs_node_ref.size();
			auto obs_other = get_not_obs(obs_node,vector <unsigned int> (),nobs);
			
			auto lin = 0.0, quad = 0.0;
			for(auto m1 : obs_node){
				const auto &gen_dif2 = gen_dif[m1];
				const auto &gen_val2 = gen_val[m1];
				for(auto m2 : obs_other){
					auto d = int(gen_dif2[m2]) - int(gen_val2[m2]);
					lin += 2*d; 
					quad++;
				}
			}
			lin *= -2*sca; quad *= -2*sca;
		
			sampler.init(lin,quad,lam);
			
			auto num = iev.mut_num;
				
			auto probif = 0.0;
			auto num_new = sampler.sample(probif);
			
			if(num_new != num){
				auto probfi = sampler.prob(num);
				
				auto dlike_genetic_process = poisson_probability(num_new,lam) - poisson_probability(num,lam);

				int dn = (int)num_new-(int)num;

				auto dlike_genetic_obs = lin*dn+quad*dn*dn;
				
				auto al = exp(burn_info.genetic_process*dlike_genetic_process
								+ burn_info.genetic_obs*dlike_genetic_obs + probfi - probif);

				ntr++;
				if(ran() < al){ 
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
				
				if(pl) check("update");
			}
		}
	}
	
	
	// Goes through the infection origins
	for(auto &io : genetic_value.inf_origin){ 
		const auto &in = inf_node[io.node];
	
		auto lam = seq_var+mut_rate*in.tdiv_start;
	
		auto obs_node = get_gen_obs_list(io.node,0);
		auto nobs = genetic_value.obs_node_ref.size();
		auto obs_other = get_not_obs(obs_node,vector <unsigned int> (),nobs);
	
		auto lin = 0.0, quad = 0.0;
		for(auto m1 : obs_node){
			const auto &gen_dif2 = gen_dif[m1];
			const auto &gen_val2 = gen_val[m1];
			for(auto m2 : obs_other){
				auto d = int(gen_dif2[m2]) - int(gen_val2[m2]);
				lin += 2*d; 
				quad++;
			}
		}
		lin *= -2*sca; quad *= -2*sca;
	
		sampler.init(lin,quad,lam);
		
		auto num = io.mut_num;
			
		auto probif = 0.0;
		auto num_new = sampler.sample(probif);
		
		if(num_new != num){
			auto probfi = sampler.prob(num);

			auto dlike_genetic_process = poisson_probability(num_new,lam) - poisson_probability(num,lam);

			int dn = (int)num_new-(int)num;

			auto dlike_genetic_obs = lin*dn+quad*dn*dn;
				
			auto al = exp(burn_info.genetic_process*dlike_genetic_process
								+ burn_info.genetic_obs*dlike_genetic_obs + probfi - probif);

			ntr++;
			if(ran() < al){ 
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
			
			if(pl) check("update");
		}
	}
	
			/*	
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
				int d, d_update;
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
	
				ntr++;
				if(ran() < al){ 
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
					
					burn_info.update_si(si,ACCEPT,ntr);
				}
				else{
					burn_info.update_si(si,REJECT,ntr); if(si < 1) si = 1;
				}
							
				if(pl) check("update");
			}
		}
	}
	*/
}


/// Gets list of observations not on two lists
vector <unsigned int> State::get_not_obs(const vector <unsigned int> &obs_node1, const vector <unsigned int> &obs_node2, unsigned int N) const
{
	vector <unsigned int> obs_node;
	
	vector <bool> map(N,false);
	for(auto m : obs_node1) map[m] = true;
	for(auto m : obs_node2) map[m] = true;
	for(auto m = 0u; m < N; m++){
		if(map[m] == false && genetic_value.obs_node_ref[m].node != UNSET){
			obs_node.push_back(m);
		}
	}
	
	return obs_node;
}

													
/// Performs local transmission tree mutation proposals (the interchange adjacent mutation events)
void State::trans_tree_mut_local_proposal(const BurnInfo &burn_info, vector <GeneticProp> &gen_mut_info)
{
	auto pl = false;
	
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
		case 5: type = TWO_BRANCH4; break;
		case 6: type = TWO_BRANCH5; break;
		case 7: type = TWO_BRANCH6; break;
		case 8: type = ORIGIN_GEN_UPDATE; break;
		default: emsg("Not op"); break;
		}
		
		auto N = gen_data.obs.size();
	
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
										int d, d_update;
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
										
										gmi.ntr++;
										if(ran() < al){ 
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
													
											burn_info.update_si(gmi.si,ACCEPT,gmi.ntr);
										}
										else{
											burn_info.update_si(gmi.si,REJECT,gmi.ntr); if(gmi.si < 1) gmi.si = 1;
										}
										
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
						 	 		 (nbranch == 2 && (type == TWO_BRANCH1 || type == TWO_BRANCH2 || type == TWO_BRANCH3 || 
															 	 		 type == TWO_BRANCH4 || type == TWO_BRANCH5 || type == TWO_BRANCH6))){
								
									int dn1 = dn, dn2 = -dn, dn3 = -dn;
									switch(type){
									case TWO_BRANCH1: dn1 = dn; dn2 = -dn; dn3 = -dn; break;
									case TWO_BRANCH2: dn1 = dn; dn2 = dn; dn3 = -dn; break;
									case TWO_BRANCH3: dn1 = dn; dn2 = -dn; dn3 = dn; break;
									case TWO_BRANCH4: dn1 = 0; dn2 = dn; dn3 = -dn; break;
									case TWO_BRANCH5: dn1 = -dn; dn2 = dn; dn3 = 0; break;
									case TWO_BRANCH6: dn1 = -dn; dn2 = 0; dn3 = dn; break;
									default: break;
									}
								
									auto ill = false;
									
									// Works out change in likelihood for genetic process
									auto dlike_genetic_process = 0.0;
									
									{
										auto num_old = ev.mut_num;
										auto num_new = int(num_old)+dn1;
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
										
										vector <unsigned int> obs_node1, obs_node2, obs_node3;
										auto &gen_dif = genetic_value.gen_dif;
										
										int dif = 0;
																			
										if(nbranch == 2){
											// Finds change in gen_dif and genetic_obs
										
											switch(type){
											case TWO_BRANCH1:
												obs_node1 = get_gen_obs_list(n,e+1);
												obs_node2 = get_gen_obs_list(ev.index,0);
												dif = 2*dn;
												break;
												
											case TWO_BRANCH2:
												{
													obs_node1 = get_gen_obs_list(n,e+1);
													auto obs_node_not = get_gen_obs_list(ev.index,0);
													obs_node2 = get_not_obs(obs_node1,obs_node_not,N);
													dif = -2*dn;
												}
												break;
												
											case TWO_BRANCH3:
												{
													auto obs_node_not = get_gen_obs_list(n,e+1);
													obs_node2 = get_gen_obs_list(ev.index,0);
													obs_node1 = get_not_obs(obs_node2,obs_node_not,N);
													dif = -2*dn;
												}
												break;
												
											case TWO_BRANCH4:
												{
													obs_node1 = get_gen_obs_list(n,e+1);
													obs_node2 = get_gen_obs_list(ev.index,0);
													obs_node3 = get_not_obs(obs_node1,obs_node2,N);
													dif = dn;
												}
												break;
											
											case TWO_BRANCH5:
												{
													obs_node1 = get_gen_obs_list(n,e+1);
													obs_node3 = get_gen_obs_list(ev.index,0);
													obs_node2 = get_not_obs(obs_node1,obs_node3,N);
													dif = dn;
												}
												break;
											
											case TWO_BRANCH6:
												{
													obs_node1 = get_gen_obs_list(ev.index,0);
													obs_node3 = get_gen_obs_list(n,e+1);
													obs_node2 = get_not_obs(obs_node1,obs_node3,N);
													dif = dn;
												}
												break;
												
											default: emsg("not def"); break;
											}
											
											auto sum = 0.0;
											int d, d_update;
											
											switch(type){
											case TWO_BRANCH1: case TWO_BRANCH2: case TWO_BRANCH3: 
												for(auto m1 : obs_node1){
													const auto &gen_val2 = gen_dif[m1];
													const auto &gen_dif2 = gen_data.gen_dif[m1];
												
													for(auto m2 : obs_node2){	
														d = int(gen_dif2[m2]) - int(gen_val2[m2]);
														d_update = d + dif;
														sum += d_update*d_update - d*d;
													}
												}
												break;

											case TWO_BRANCH4: case TWO_BRANCH5: case TWO_BRANCH6:
												for(auto m1 : obs_node1){
													const auto &gen_val2 = gen_dif[m1];
													const auto &gen_dif2 = gen_data.gen_dif[m1];
												
													for(auto m2 : obs_node3){	
														d = int(gen_dif2[m2]) - int(gen_val2[m2]);
														d_update = d - dif;
														sum += d_update*d_update - d*d;
													}
												}
												
												for(auto m1 : obs_node2){
													const auto &gen_val2 = gen_dif[m1];
													const auto &gen_dif2 = gen_data.gen_dif[m1];
												
													for(auto m2 : obs_node3){	
														d = int(gen_dif2[m2]) - int(gen_val2[m2]);
														d_update = d + dif;
														sum += d_update*d_update - d*d;
													}
												}
												break;
												
											default: emsg("not def"); break;
											}
						
											dlike_genetic_obs += -2*sum*model.genetic_data.obs_scale;
										}
										
										auto al = exp(burn_info.genetic_process*dlike_genetic_process
																+ burn_info.genetic_obs*dlike_genetic_obs);
										
										if(pl) cout << al << " " << dlike_genetic_process << " " << dlike_genetic_obs << " al" << endl;
									
										gmi.ntr++;
										if(ran() < al){
											gmi.nac++;
										
											like.genetic_process += dlike_genetic_process;
											like.genetic_obs += dlike_genetic_obs;
										
											ev.mut_num += dn1;	
											if(e+1 < in.inf_ev.size()) in.inf_ev[e+1].mut_num += dn2;
											if(in2.inf_ev.size() > 0) in2.inf_ev[0].mut_num += dn3;
											
											if(nbranch == 2){
												switch(type){
												case TWO_BRANCH1: case TWO_BRANCH2: case TWO_BRANCH3: 
													for(auto m1 : obs_node1){
														for(auto m2 : obs_node2){	
															gen_dif[m1][m2] -= dif;
															gen_dif[m2][m1] -= dif;
														}
													}
													break;
													
												case TWO_BRANCH4: case TWO_BRANCH5: case TWO_BRANCH6:
													for(auto m1 : obs_node1){
														for(auto m2 : obs_node3){	
															gen_dif[m1][m2] += dif;
															gen_dif[m2][m1] += dif;
														}
													}
													
													for(auto m1 : obs_node2){
														for(auto m2 : obs_node3){	
															gen_dif[m1][m2] -= dif;
															gen_dif[m2][m1] -= dif;
														}
													}
													break;
												
												default: emsg("not def"); break;
												}
											}
											
											burn_info.update_si(gmi.si,ACCEPT,gmi.ntr);
										}
										else{
											burn_info.update_si(gmi.si,REJECT,gmi.ntr); if(gmi.si < 1) gmi.si = 1;
										}
										
										if(pl) check("update21");
									}
								}
							}
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
						
						gmi.ntr++;
						if(ran() < al){ 
							gmi.nac++;
						
							like.genetic_process += dlike_genetic_process;
						
							io.mut_num = num_new;	
							ev.mut_num = num2_new;
							
							burn_info.update_si(gmi.si,ACCEPT,gmi.ntr);
						}
						else{
							burn_info.update_si(gmi.si,REJECT,gmi.ntr); if(gmi.si < 1) gmi.si = 1;
						}
				
						if(pl) check("update31");
					}
				}
			}
		}
	}
}


/// Prints information about all the nodes
void State::print_all_node() const
{
	for(auto n = 0u; n < genetic_value.inf_node.size(); n++){
		print_node("Node"+tstr(n),n);
	}
}
				

/// Prints information about an infection node
void State::print_node(string te, unsigned int n) const 
{
	if(n == OUTSIDE_INF){ cout << te << ": OUTSIDE_INF" << endl; return;}
	if(n == ENTER_INF){ cout << te << ": ENTER_INF" << endl; return;}
	
	const auto &in = genetic_value.inf_node[n];
		
	cout << te << ": Infection node " << n << endl;
	/*
	cout << "t_start=" << model.calc_t(in.tdiv_start);
	cout << "  t_rec=" << model.calc_t(in.tdiv_rec);
	cout << "  ti_start=" << model.calc_t(in.ti_start);
	cout << "  ti_rec=" << model.calc_t(in.ti_rec);
	*/
	cout << "t_start=" << in.tdiv_start;
	cout << "  t_rec=" << in.tdiv_rec;
	cout << "  ti_start=" << in.ti_start;
	cout << "  ti_rec=" << in.ti_rec;
	
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
		//cout << "t=" << model.calc_t(ie.tdiv) << " (tdiv=" << ie.tdiv << ")  ";
		cout << "t=" << ie.tdiv << " (tdiv=" << ie.tdiv << ")  ";
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
unsigned int State::inf_node_ti_from_tdiv(double tdiv) const 
{
	if(tdiv == 0) return 0;
	return 1+(unsigned int)(tdiv);
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
		if(mut_num_store == UNSET) emsg("Should not be unset1");
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
	
	auto jmax = inf_ev2.size();
	auto j = 0u; while(j < jmax && inf_ev2[j].tdiv < t) j++;
	
	NodeRef nr; nr.node = n2; nr.index = j;
	
	return nr;
}


/// Checks that pop_ind is correctly specified
void State::check_pop_ind(const vector < vector <PopInd> > &pop_ind, const vector < vector <IndPop> > &ind_pop) 
{
	cout << "chec pop_ind" << endl;
	check_timer[CHECK_POP_IND_GENTIC] -= clock();
	
	vector < vector <unsigned int> > pop_ind_compare;
	pop_ind_compare.resize(model.npop);
	
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
	
	for(auto po = 0u; po < model.npop; po++){
		vector <unsigned int> pop_ind_list;
		for(const auto &pi : pop_ind[po]) pop_ind_list.push_back(pi.i);
		
		sort(pop_ind_list.begin(),pop_ind_list.end());
		sort(pop_ind_compare[po].begin(),pop_ind_compare[po].end());
		
		if(!equal_vec(pop_ind_list,pop_ind_compare[po])) emsg("pop_ind problem");
	}
	
	for(auto po = 0u; po < model.npop; po++){
		const auto &pop = model.pop[po];
		auto p = pop.p;
		
		for(auto k = 0u; k < pop_ind[po].size(); k++){
			const auto &pi = pop_ind[po][k];
			const auto &pr = ind_pop[p][pi.i].pop_ref[pi.index];
			if(pr.po != po) emsg("Population not right");
			if(pr.index != k) emsg("index wring");
		}
	}
	
	for(auto p = 0u; p < model.species.size(); p++){
		const auto &ssp = species[p];
		for(auto i = 0u; i < ssp.individual.size(); i++){	
			const auto &ip = ind_pop[p][i];
			for(auto k = 0u; k < ip.pop_ref.size(); k++){
				const auto &pr = ip.pop_ref[k];
				const auto &pi = pop_ind[pr.po][pr.index];
				if(pi.i != i) emsg("ind wring");
				if(pi.index != k) emsg("index wring");
			}
		}
	}
	
	check_timer[CHECK_POP_IND_GENTIC] += clock();
}


/// Indicates no tranmssion tree change
TransTreeChange State::no_trans_tree_change() const
{	
	TransTreeChange gc;
	gc.probif = 0; 
	gc.probfi = 0; 
	gc.mut_num = UNSET; 
	gc.mutA_num = UNSET; 
	gc.dlike_genetic_process = 0; 
	gc.dlike_genetic_obs = 0;
	gc.type = NO_TRANS_TREE_CHANGE;
	
	return gc;
}


/// Indicates no tranmssion tree change
TransTreeChange State::trans_tree_fail() const
{	
	TransTreeChange gc;
	gc.type = TRANS_TREE_FAIL;
	
	return gc;
}


/// Returns a function which can be used to calculate the genetic observation likelihood for when the reconnection is
GenObsFunc State::get_genetic_obs_func_from(const vector <unsigned int> &obs_node, unsigned int n) const
{
	GenObsFunc gof;
	
	auto &inf_node = genetic_value.inf_node;
	const auto &in = inf_node[n];
	
	auto nobs = genetic_value.obs_node_ref.size();
	
	const auto &gen_val = model.genetic_data.gen_dif;
	const auto &gen_dif = genetic_value.gen_dif;
			
	const auto &nr_from = in.from;
	
	auto fac = -2*model.genetic_data.obs_scale;
	
	auto con_mut = 0.0, lin_mutA = 0.0, lin_mutB = 0.0, quad_mutAB = 0.0, quad_mut = 0.0;
	
	vector <unsigned int> obs_nodeC;
	unsigned int mutB;
	
	auto mutA = in.inf_ev[0].mut_num;
	
	auto n_from = nr_from.node;
	auto ind = nr_from.index;
	
	if(n_from == OUTSIDE_INF || n_from == ENTER_INF){
		mutB = genetic_value.inf_origin[ind].mut_num;
	}
	else{		
		obs_nodeC = get_gen_obs_list(n_from,ind+1); 
		mutB = inf_node[n_from].inf_ev[ind].mut_num;
	}
	
	auto obs_nodeB = get_not_obs(obs_node,obs_nodeC,nobs);
		
	int d, d_update;
	if(obs_nodeB.size() > 0){
		int ch = -mutB-mutA;
		for(auto mA : obs_node){
			const auto &gen_val2 = gen_val[mA];
			const auto &gen_dif2 = gen_dif[mA];
			
			for(auto m : obs_nodeB){
				auto d = int(gen_dif2[m]) - int(gen_val2[m]);
				d_update = d + ch;
				
				con_mut += d_update*d_update-d*d;
				lin_mutA += 2*d_update;
				lin_mutB += 2*d_update;
				quad_mutAB++;
				quad_mut++;
			}
		}
	}
	
	if(obs_nodeC.size() > 0){
		//auto mutC = inf_node[n_from].inf_ev[ind+1].mut_num;
	
		//int ch = mutB+mutC-mutA;
		int ch = mutB-mutA;
		for(auto mA : obs_node){
			const auto &gen_val2 = gen_val[mA];
			const auto &gen_dif2 = gen_dif[mA];
			
			for(auto m : obs_nodeC){
				d = int(gen_dif2[m]) - int(gen_val2[m]);
				d_update = d + ch;
				
				con_mut += d_update*d_update-d*d;
				lin_mutA += 2*d_update;
				lin_mutB -= 2*d_update;
				quad_mutAB--;
				quad_mut++;
			}
		}
	}
	
	gof.con_mut = con_mut*fac;
	gof.lin_mutA = lin_mutA*fac;
	gof.lin_mutB = lin_mutB*fac;
	gof.quad_mutAB = quad_mutAB*fac;
	gof.quad_mut = quad_mut*fac;

	return gof;
}


/// Returns a function which can be used to calculate the genetic observation likelihood for when the reconnection is
GenObsFunc State::get_genetic_obs_func_to(const vector <unsigned int> &obs_node, unsigned int n, const NodeRef &nr_add, double t_add, TransTreeChange &gc, unsigned int &mutBC_num) const 
{
	GenObsFunc gof;
	
	auto &inf_node = genetic_value.inf_node;
	const auto &in = inf_node[n];
	auto t_now = in.tdiv_start;
	auto nobs = genetic_value.obs_node_ref.size();
	
	const auto &gen_val = model.genetic_data.gen_dif;
	const auto &gen_dif = genetic_value.gen_dif;
			
	const auto &nr_from = in.from;
	
	auto fac = -2*model.genetic_data.obs_scale;

	// Works out information for the end point of the node
	
	vector <NodeRef> path_from, path_add;

	path_add.push_back(nr_add);

	auto n_from = n, n_add = nr_add.node;
	auto term_from = false, term_add = false;
	
	vector <ObsMut> obs_mut, obs_mutC;
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
			
			if(n_from == n_add && path_add.size() == 1 && t_add < t_now && path_add[0].index <= path_from[path_from.size()-1].index){
				plug_in = true;
			}
			
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
							//end_type = C_GEN_END;
							mutBC_num = inf_ev[j].mut_num;
							//mut_sum += inf_ev[j].mut_num - gc.mut_num;
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
				if(n_from != n_add){  // Needs to add jumping from one origin to another
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
								//end_type = B_GEN_END;
								//mut_sum += gc.mut_num;
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
						obs_mutC.push_back(om);
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
			
			//mut_sum += gc.mut_num;
			//end_type = B_GEN_END;
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
					if(j == jmax && k == 0){
						//end_type = B_GEN_END;
						//mut_sum += gc.mut_num;
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
					obs_mutC.push_back(om);
				}
			}
		}
		break;
	
	case FROM_UNCON_ADD_UNCON:
		emsg("Should not be uncon");
		break;
	}
	
	{ // Works out observations branching from root
		vector <bool> obs_done(nobs,false);
		for(auto om : obs_mut) obs_done[om.obs] = true;
		for(auto om : obs_mutC) obs_done[om.obs] = true;
			
		for(auto m : obs_node) obs_done[m] = true;
		
		for(auto m = 0u; m < nobs; m++){
			if(obs_done[m] == false && genetic_value.obs_node_ref[m].node != UNSET){
				ObsMut om; om.mut_num = mut_sum_root; om.obs = m;
				if(mutBC_num != UNSET) obs_mutC.push_back(om);
				else obs_mut.push_back(om);
			}
		}
	}
	
	gof.mut_sum = mut_sum;
	
	auto &gen_dif_cha = gc.gen_dif_cha;
	
	auto mutA = in.inf_ev[0].mut_num;

	auto con_mut = 0.0, lin_mutA = 0.0, lin_mutB = 0.0, quad_mutAB = 0.0, quad_mut = 0.0;
	
	for(const auto &om : obs_mut){     // Works out change to gen_dif
		int ch = mut_sum-mutA-2*om.mut_num;
	
		GenDifCha gdc; 
		gdc.ob1 = om.obs;
		gdc.ch = ch;
		gdc.sign = 1;
		
		const auto &gen_val2 = gen_val[om.obs];
		const auto &gen_dif2 = gen_dif[om.obs];
	
		int d, d_update;
		for(auto m : obs_node){
			d = int(gen_dif2[m]) - int(gen_val2[m]);
			d_update = d + ch;

			if(mutBC_num == UNSET){
				con_mut += d_update*d_update - d*d;	
				lin_mutA += 2*d_update;
				lin_mutB += 2*d_update;
				quad_mutAB++;
				quad_mut++;
			}
			else{
				d_update += mutBC_num;
				con_mut += d_update*d_update - d*d;
				lin_mutA += 2*d_update;
				lin_mutB -= 2*d_update;
				quad_mutAB--;
				quad_mut++;
			}
			
			gdc.ob2.push_back(m);
			if(gdc.ob1 == m) emsg("Should not be same");
		}		
		gen_dif_cha.push_back(gdc);
	}				

	for(const auto &om : obs_mutC){     // Works out change to gen_dif
		int ch = -mut_sum-mutA;
		
		GenDifCha gdc; 
		gdc.ob1 = om.obs;
		gdc.ch = ch;
		gdc.sign = -1;
		const auto &gen_val2 = gen_val[om.obs];
		const auto &gen_dif2 = gen_dif[om.obs];
		
		int d, d_update;
		for(auto m : obs_node){
			d = int(gen_dif2[m]) - int(gen_val2[m]);
			d_update = d + ch;
			
			if(mutBC_num == UNSET){
				con_mut += d_update*d_update - d*d;
				lin_mutA += 2*d_update;
				lin_mutB -= 2*d_update;
				quad_mutAB--;
				quad_mut++;
			}
			else{
				d_update -= mutBC_num;
				con_mut += d_update*d_update - d*d;
				lin_mutA += 2*d_update;
				lin_mutB += 2*d_update;
				quad_mutAB++;
				quad_mut++;
			}
			
			gdc.ob2.push_back(m);
			if(gdc.ob1 == m) emsg("Should not be same");	
		}			
		gen_dif_cha.push_back(gdc);				
	}
	
	gof.con_mut = con_mut*fac;
	gof.lin_mutA = lin_mutA*fac;
	gof.lin_mutB = lin_mutB*fac;
	gof.quad_mutAB = quad_mutAB*fac;
	gof.quad_mut = quad_mut*fac;
	
	return gof;
}


/// Calculate the change in genetic_obs for a given mutation number
double State::calc_genetic_obs_fast(int mutB_num, int mutA_num, const GenObsFunc &gof) const
{
	return gof.con_mut + gof.lin_mutA*mutA_num + gof.lin_mutB*mutB_num + 2*gof.quad_mutAB*mutA_num*mutB_num + gof.quad_mut*(mutA_num*mutA_num + mutB_num*mutB_num);
}


/// Works out how genetics change as a result of a node being regraphed / added / removed
TransTreeChange State::trans_tree_change(TransTreeChaType type, unsigned int n, const InfNodeAlter &node_alter, const IndInfFrom &iif_add, RegraftStore rs) const 
{		
	auto data_on = model.genetic_data.on;
		
	auto mut_rate = genetic_value.mut_rate*model.details.dt;
	auto seq_var = genetic_value.seq_var;
	
	TransTreeChange gc;
	gc.probif = 0; 
	gc.probfi = 0; 
	gc.mut_num = UNSET;
	gc.mutA_num = UNSET; 
	gc.dlike_genetic_process = 0; 
	gc.dlike_genetic_obs = 0;
	gc.type = type; 

	gc.n = n;
	gc.iif_add = iif_add;
	gc.node_alter = node_alter;
	
	auto &regraft_type = gc.regraft_type;				
	regraft_type = NORMAL_REGRAFT;
	if(rs == REGRAFT_RETRIEVE) gc.rev_store = node_alter.rev_store;
	
	auto t_add = node_alter.tdiv_start;

	const auto &inf_node = genetic_value.inf_node;

	auto &nr_add = gc.nr_add;

	switch(type){
	case TRANS_TREE_FAIL: emsg("Should not fail"); break;
	
	case REGRAFT_NODE: 
		{
			nr_add = get_ind_noderef(iif_add.p,iif_add.i,t_add); 
		
			const auto &in = inf_node[n];
			const auto &inf_ev = in.inf_ev;
	
			auto t_now = in.tdiv_start;
			
			const auto &nr_from = in.from;

			auto mutA_start = UNSET; if(inf_ev.size() > 0) mutA_start = inf_ev[0].mut_num;
			gc.mutA_num = mutA_start;
			
			if(nr_from.node > CODE && nr_add.node > CODE){ // Starts and ends independent from tree
				nr_add.index = nr_from.index;
				regraft_type = ORIGIN_REGRAFT;
			
				if(t_add != t_now && data_on){
					auto mut_num = genetic_value.inf_origin[nr_from.index].mut_num;
					
					gc.dlike_genetic_process += poisson_probability(mut_num,seq_var+mut_rate*t_add)
								-poisson_probability(mut_num,seq_var+mut_rate*t_now);
				}
			}					
			else{	
				if(nr_from.node == nr_add.node && (nr_from.index == nr_add.index || nr_from.index+1 == nr_add.index)){ 
					// Attaches to the same point on the tree (here mutation numbers are unchanged)
					
					regraft_type = UNCHANGED_REGRAFT;
					
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
				}
				else{		
					// Either from or to are attached to the tree
				
					if(data_on){  // Works out changes to gen_dif 
						auto obs_node = get_gen_obs_list(n,0); 
						
						if(rs == REGRAFT_RETRIEVE){
							regraft_type = RegraftType(gc.rev_store[0]);
						}
						else{
							if(obs_node.size() > 0){ 
								switch(regraft_obs_type){
								case RG_OBS: regraft_type = REGRAFT_WITH_OBS; break;
								case RG_OBS2: regraft_type = REGRAFT_WITH_OBS2; break;
								case RG_OBS3: regraft_type = REGRAFT_WITH_OBS3; break;
								case RG_NOT_SET: emsg("RG obs not set"); break;
								}
							}
							
							if(rs == REGRAFT_STORE) gc.rev_store.push_back(regraft_type);
						}
							
						if(regraft_type == NORMAL_REGRAFT){
							change_add_node(t_add,nr_add,gc,rs);      // Adds new node
							change_remove_node(t_now,nr_from,gc,rs);  // Removes node				
						}
						else{
							auto mutBC_num = UNSET;
							
							// Gets genetic_obs as function of the (potentially) two mutation rates
							auto gof_to = get_genetic_obs_func_to(obs_node,n,nr_add,t_add,gc,mutBC_num); 
																
							switch(regraft_type){
							case REGRAFT_WITH_OBS:    // This keeps mutA_num fixed
								{
									change_add_node(t_add,nr_add,gc,rs);      // Adds new node
									change_remove_node(t_now,nr_from,gc,rs);  // Removes node
								}
								break;
									
							case REGRAFT_WITH_OBS2:   // This sample mutA_num from a range detemined by mut_sum
								{
									change_add_node(t_add,nr_add,gc,rs);      // Adds new node
									change_remove_node(t_now,nr_from,gc,rs);  // Removes node
									
									if(rs == REGRAFT_RETRIEVE){               // Samples gc.mutA_num
										gc.mutA_num = gc.rev_store[2];
									}
									else{
										if(rs == REGRAFT_STORE) gc.rev_store.push_back(gc.mutA_num);
										
										auto ms = gof_to.mut_sum;
										if(mutBC_num == UNSET) ms += gc.mut_num;
										else ms += mutBC_num-gc.mut_num;
		
										int sh = int(ran()*(2*ms+1))-(int)ms;
										
										if(-sh > (int)gc.mutA_num) return trans_tree_fail();
										gc.mutA_num += sh;
									}
								}
								break;
								
							case REGRAFT_WITH_OBS3:		
								{
									if(!change_add_node_with_gobs(n,t_add,nr_add,gc,gof_to,rs)){ // Adds new node
										return trans_tree_fail();
									}
								
									auto gof_from = get_genetic_obs_func_from(obs_node,n); 
									
									if(!change_remove_node_with_gobs(n,t_now,nr_from,gc,gof_from,rs)){     // Removes node	
										return trans_tree_fail();
									}
								}
								break;
								
							default:
								emsg("SHould not be");
								break;
							}

							auto mut = gc.mut_num;                 // Corrects gen_dif changes based on sampled mutation
							auto mutA = gc.mutA_num;
							for(auto &cha : gc.gen_dif_cha){ 
								auto &ch = cha.ch;
								ch += mutA;
								if(mutBC_num == UNSET) ch += mut*cha.sign;
								else ch += (mutBC_num - mut)*cha.sign;
							}
							
							gc.dlike_genetic_obs += calc_genetic_obs_fast(mut,mutA,gof_to);
						}
					}
				}
			}
			
			{  // Works out change as a result of alteration	
				if(data_on){
					if((t_now != t_add || gc.mutA_num != mutA_start) && gc.mutA_num != UNSET){
						const auto &iev = inf_ev[0];
						gc.dlike_genetic_process += poisson_probability(gc.mutA_num,mut_rate*(iev.tdiv-t_add)) 
																			- poisson_probability(mutA_start,mut_rate*(iev.tdiv-t_now));
					}
					
					gc.dlike_genetic_obs -= GEN_OBS_WRONG*(node_alter.obs_add_end.size() + 
																					 node_alter.obs_add_begin.size());
				}
			}
		}
		break;
		
	case REGRAFT_GEN_OBS_MAX:
		{
			regraft_type = REGRAFT_NOT_NEEDED;
		
			nr_add = get_ind_noderef(iif_add.p,iif_add.i,t_add); 
		
			const auto &in = inf_node[n];
			const auto &inf_ev = in.inf_ev;
	
			auto t_now = in.tdiv_start;
			
			const auto &nr_from = in.from;

			if(!data_on) emsg("Gen data should be on");

			auto mutA_start = UNSET; if(inf_ev.size() > 0) mutA_start = inf_ev[0].mut_num;
			gc.mutA_num = mutA_start;
			
			if(!(nr_from.node > CODE && nr_add.node > CODE) &&
				 !(nr_from.node == nr_add.node && (nr_from.index == nr_add.index || nr_from.index+1 == nr_add.index))){			
				auto obs_node = get_gen_obs_list(n,0); 
				
				if(obs_node.size() > 0){
					regraft_type = REGRAFT_WITH_OBS;
						
					auto mutBC_num = UNSET;
					
					// Gets genetic_obs as function of the (potentially) two mutation rates
					auto gof_to = get_genetic_obs_func_to(obs_node,n,nr_add,t_add,gc,mutBC_num); 

					change_remove_node(t_now,nr_from,gc);  	

					change_add_node_gen_obs_max(n,t_add,nr_add,gc,gof_to);        // Adds new node and maximises observation prob

					auto mut = gc.mut_num;                 // Corrects gen_dif changes based on sampled mutation
					auto mutA = gc.mutA_num;
					for(auto &cha : gc.gen_dif_cha){ 
						auto &ch = cha.ch;
						ch += mutA;
						if(mutBC_num == UNSET) ch += mut*cha.sign;
						else ch += (mutBC_num - mut)*cha.sign;
					}
					
					gc.dlike_genetic_obs += calc_genetic_obs_fast(mut,mutA,gof_to);

					if((t_now != t_add || gc.mutA_num != mutA_start) && gc.mutA_num != UNSET){
						const auto &iev = inf_ev[0];
						gc.dlike_genetic_process += poisson_probability(gc.mutA_num,mut_rate*(iev.tdiv-t_add)) 
																			- poisson_probability(mutA_start,mut_rate*(iev.tdiv-t_now));
					}
					
					gc.dlike_genetic_obs -= GEN_OBS_WRONG*(node_alter.obs_add_end.size() + 
																					 node_alter.obs_add_begin.size());
				}
			}
		}
		break;
	
	case ADD_NODE:
		{
			nr_add = get_ind_noderef(iif_add.p,iif_add.i,t_add); 
		
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
		
	case NO_TRANS_TREE_CHANGE:
		break;
	}
	
	// Makes sure that inf_ev times are not equal
	switch(type){ 
	case REGRAFT_NODE: case REGRAFT_GEN_OBS_MAX: case ADD_NODE:
		switch(regraft_type){
		case NORMAL_REGRAFT: case REGRAFT_WITH_OBS: case REGRAFT_WITH_OBS2: case REGRAFT_WITH_OBS3:
			{
				auto no = nr_add.node;
				if(no != OUTSIDE_INF && no != ENTER_INF){
					auto ind = nr_add.index;
					if(ind < inf_node[no].inf_ev.size()){
						if(t_add == inf_node[no].inf_ev[ind].tdiv) return trans_tree_fail();
					}
				}
			}
			break;
		default: break;
		}
		break;
	default: break;
	}
		
	return gc;
}
					

/// Adding node to new branch
void State::change_add_node(double t, const NodeRef &nr_add, TransTreeChange &gc, RegraftStore rs) const 
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
			if(rs == REGRAFT_RETRIEVE) val = gc.rev_store[1];
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
		
			if(rs == REGRAFT_RETRIEVE) val = gc.rev_store[1];
			if(val == UNSET) emsg("Problem sampling from a Poisson distribition: "+warn);
			gc.mut_num = val;
					
			gc.dlike_genetic_process += poisson_probability(gc.mut_num,mut_rate*(t2-t1)) 
																+ poisson_probability(num-gc.mut_num,mut_rate*(t3-t2))
																- poisson_probability(num,mut_rate*(t3-t1));            
			if(gc.mut_num > num) emsg("h1");																
			gc.probif += binomial_probability(gc.mut_num,frac,num);
		}
	}
	else{
		auto seq_var = genetic_value.seq_var;
		
		auto lam = seq_var+mut_rate*t;
		auto val = poisson_sample(lam,warn);
		if(rs == REGRAFT_RETRIEVE) val = gc.rev_store[1];
		if(val == UNSET) emsg("Problem sampling from a Poisson distribition: "+warn);
		gc.mut_num = val;
		auto prob = poisson_probability(gc.mut_num,lam);
		gc.dlike_genetic_process += prob;
		gc.probif += prob;
	}
}
				
			
/// Removing node from branch
void State::change_remove_node(double t, const NodeRef &nr_from, TransTreeChange &gc, RegraftStore rs) const 		
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
			if(rs == REGRAFT_STORE) gc.rev_store.push_back(num);
			auto prob = poisson_probability(num,mut_rate*(t2-t1));
			gc.dlike_genetic_process -= prob; 
			gc.probfi += prob;
		}
		else{
			auto t3 = inf_ev_from[j_from+1].tdiv;
			auto numA = inf_ev_from[j_from].mut_num;
			auto numB = inf_ev_from[j_from+1].mut_num;

			if(rs == REGRAFT_STORE) gc.rev_store.push_back(numA);

			gc.dlike_genetic_process += poisson_probability(numA+numB,mut_rate*(t3-t1)) 
																- poisson_probability(numA,mut_rate*(t2-t1)) 
																- poisson_probability(numB,mut_rate*(t3-t2)); 
															
			auto frac = (t2-t1)/(t3-t1);
			if(numA> numA+numB) emsg("h2");			
			gc.probfi += binomial_probability(numA,frac,numA+numB);
		}
	}
	else{
		auto seq_var = genetic_value.seq_var;
	
		auto lam = seq_var+mut_rate*t;
		auto num = genetic_value.inf_origin[nr_from.index].mut_num;
		if(rs == REGRAFT_STORE) gc.rev_store.push_back(num);
	
		auto prob = poisson_probability(num,lam);
		gc.dlike_genetic_process -= prob;
		gc.probfi += prob;
	}
}


/// Based on a gof, this returns the approximate sampler
void State::get_regraft_minimum(int &i, int &j, const GenObsFunc &gof) const 
{
	auto b = gof.lin_mutA, c = gof.lin_mutB, d = gof.quad_mutAB, e = gof.quad_mut;
	
	auto denom = e*e-d*d;
	if(denom == 0){
		auto x = -0.5*b/e;
		i = round_int(x);
		j = UNSET;
	}
	else{		
		auto x = 0.5*(b*d-c*e)/denom;
		auto y = 0.5*(c*d-b*e)/denom;
		
		i = round_int(x);
		j = round_int(y);
	}
}


/// Based on a gof, this returns the approximate sampler
RGsamp State::get_regraft_sampler(const GenObsFunc &gof) const 
{
	auto b = gof.lin_mutA, c = gof.lin_mutB, d = gof.quad_mutAB, e = gof.quad_mut;
	
	RGsamp rgs;
	auto &ele = rgs.ele;
	vector <RGsamp> rg_samp;
	
	auto denom = e*e-d*d;
	if(denom == 0){
		auto x = -0.5*b/e;
		auto i = round_int(x);
		
		auto xi = (unsigned int)((x-i+0.5)*REGRAPFT_SAMPLER_DIV);
		
		auto var = -0.5/e;	

		auto vi = floor_int(((log(var)-REGRAPFT_SAMPLER_LOGVAR_MIN)/(REGRAPFT_SAMPLER_LOGVAR_MAX-REGRAPFT_SAMPLER_LOGVAR_MIN))*REGRAPFT_SAMPLER_NVAR);
		if(vi < 0) vi = 0;
		if(vi >= (int)REGRAPFT_SAMPLER_NVAR) return rgs;

		const auto &rsamp = model.genetic_data.regraft_sampler[model.genetic_data.regraft_sampler_1d_ref[xi][vi]];

		auto sum = 0.0;
		for(const auto &el : rsamp.ele){ 
			auto B = i + el.di;
			if(B >= 0){
				sum += el.prob;
				
				RGsampEle rgel;
				rgel.B = B;
				rgel.prob = el.prob;
				rgel.sum = sum;
				ele.push_back(rgel);
			}
		}
		rgs.sum_tot = sum;
		
		rgs.type = rsamp.type;
	}
	else{		
		auto x = 0.5*(b*d-c*e)/denom;
		auto y = 0.5*(c*d-b*e)/denom;
		
		auto i = round_int(x);
		auto j = round_int(y);
		
		auto xi = (unsigned int)((x-i+0.5)*REGRAPFT_SAMPLER_DIV);
		auto yi = (unsigned int)((y-j+0.5)*REGRAPFT_SAMPLER_DIV);
		
		auto var = -0.5*e/denom;	
		
		auto vi = floor_int(((log(var)-REGRAPFT_SAMPLER_LOGVAR_MIN)/(REGRAPFT_SAMPLER_LOGVAR_MAX-REGRAPFT_SAMPLER_LOGVAR_MIN))*REGRAPFT_SAMPLER_NVAR);
		
		if(vi < 0) vi = 0;
		if(vi >= (int)REGRAPFT_SAMPLER_NVAR) return rgs;
		
		auto cor = -d/e;
		auto ci = (unsigned int)(((cor+1)/2)*REGRAPFT_SAMPLER_NCOR);
	
		const auto &rsamp = model.genetic_data.regraft_sampler[model.genetic_data.regraft_sampler_2d_ref[xi][yi][vi][ci]];

		auto sum = 0.0;
		for(const auto &el : rsamp.ele){ 
			auto B = i + el.di;
			if(B >= 0){
				auto A = j + el.dj;
				if(A >= 0){
					sum += el.prob;
					
					RGsampEle rgel;
					rgel.B = B;
					rgel.A = A;
					rgel.prob = el.prob;
					rgel.sum = sum;
					ele.push_back(rgel);
				}
			}
		}
		rgs.sum_tot = sum;
		
		rgs.type = rsamp.type;
	}

	return rgs;
}


/// Samples from a 1D regraft sampler
double RGsamp::sample1D(unsigned int &B, unsigned int &A, double r1, double r2, string &warn, const TransTreeChange &gc, RegraftStore rs) const
{	
	unsigned int k;
	
	if(rs == REGRAFT_RETRIEVE) k = gc.rev_store[1];
	else{
		auto z = ran()*sum_tot;
		auto kmax = ele.size();
		k = 0u; 
		while(k < kmax && z > ele[k].sum) k++;
		if(k == kmax) emsg("k prob1");
	}


	const auto &el = ele[k];
	
	auto N = el.B;

	if(N == 0){ 
		B = 0; A = 0;
		return log(el.prob/sum_tot);
	}
	
	if(r1 < 0 || r2 < 0){
		cout << r1 << " " << r2 << " r12" << endl;
		emsg("neg num");
	}
	auto pr = r1/(r1+r2);
	
	if(rs == REGRAFT_RETRIEVE) B = gc.rev_store[2];
	else B = binomial_sample(pr,N,warn);
	
	A = N-B;
	
	if(B > N) emsg("h3");			
	
	return log(el.prob/sum_tot)+binomial_probability(B,pr,N);
}


/// Samples from a 2D regraft sampler
double RGsamp::sample2D(unsigned int &B, unsigned int &A, const TransTreeChange &gc, RegraftStore rs) const
{
	unsigned int k;
	
	if(rs == REGRAFT_RETRIEVE) k = gc.rev_store[1];
	else{
		auto z = ran()*sum_tot;
		auto kmax = ele.size();
		k = 0u; 
		while(k < kmax && z > ele[k].sum) k++;
		if(k == kmax) emsg("k prob2");
	}
	
	const auto &el = ele[k];
	
	B = el.B; 
	A = el.A;
		
 	return log(el.prob/sum_tot);
}


/// Samples from a 1D regraft sampler
double RGsamp::prob1D(unsigned int B, unsigned int A, double r1, double r2, bool &err, TransTreeChange &gc, RegraftStore rs) const
{
	auto N = B+A;
	
	auto kmax = ele.size(); 
	auto k = 0u;
	while(k < kmax && ele[k].B != N) k++;

	if(k == kmax){ err = true; return UNSET;}

	if(rs == REGRAFT_STORE) gc.rev_store.push_back(k);
	
	const auto &el = ele[k];
	
	if(N == 0) return log(el.prob/sum_tot);
	
	if(r1 <= 0 || r2 <= 0) emsg("neg num");
	auto pr = r1/(r1+r2);

	if(rs == REGRAFT_STORE) gc.rev_store.push_back(B);

	return log(el.prob/sum_tot)+binomial_probability(B,pr,N);
}


/// Samples from a 2D regraft sampler
double RGsamp::prob2D(unsigned int B, unsigned int A, bool &err, TransTreeChange &gc, RegraftStore rs) const
{
	auto kmax = ele.size(); 
	auto k = 0u;
	while(k < kmax && !(ele[k].B == B && ele[k].A == A)) k++;
	if(k == kmax){ err = true; return UNSET;}
	
	if(rs == REGRAFT_STORE) gc.rev_store.push_back(k);
	
 	return log(ele[k].prob/sum_tot);
}


/// Outputs a comparison between approximate sampler and true function
void State::output_sampler_approximation(const RGsamp &rgs, const GenObsFunc &gof) const
{
	cout << endl << "sampler approximation:" << endl;
	if(rgs.type == REGRAFT_2D_SAMPLER){
		auto sum = 0.0;
		for(auto &el : rgs.ele){
			sum += exp(calc_genetic_obs_fast(el.B,el.A,gof));
		}

		auto sh = log(sum);

		for(auto &el : rgs.ele){
			cout << el.B << " " << el.A << " " << calc_genetic_obs_fast(el.B,el.A,gof)-sh << " " << log(el.prob/rgs.sum_tot) << endl;
		}
	}
	else{
		auto sum = 0.0;
		for(auto &el : rgs.ele){
			sum += exp(calc_genetic_obs_fast(el.B,0,gof));
		}

		auto sh = log(sum);
		
		for(auto &el : rgs.ele){
			cout << el.B << " " << calc_genetic_obs_fast(el.B,0,gof)-sh << " " << log(el.prob/rgs.sum_tot) << endl;
		}
	}
}


/// Adding node to new branch (accounting for genetic observation)
bool State::change_add_node_with_gobs(unsigned int n, double t, const NodeRef &nr_add, TransTreeChange &gc, const GenObsFunc &gof, RegraftStore rs) const 
{
	const auto &inf_node = genetic_value.inf_node;
	auto mut_rate = genetic_value.mut_rate*model.details.dt;
	
	auto rgs = get_regraft_sampler(gof);
	if(rgs.ele.size() == 0) return false;
	
	if(false) output_sampler_approximation(rgs,gof);
		
	auto type = rgs.type;
	
	string warn;
	unsigned int B, A;
	
	if(nr_add.node != OUTSIDE_INF && nr_add.node != ENTER_INF){
		const auto &in_add = inf_node[nr_add.node];
		const auto &inf_ev_add = in_add.inf_ev;
		auto j_add = nr_add.index;
		
		double t1; 
		if(j_add > 0) t1 = inf_ev_add[j_add-1].tdiv;
		else t1 = in_add.tdiv_start;
		auto t2 = t;
		
		if(j_add == inf_ev_add.size()){                // Adds at the end
			//cout << "add1" << endl;
			if(type != REGRAFT_1D_SAMPLER) emsg("Should be 1d");
			if(inf_node[n].inf_ev[0].tdiv-t2 == 0) emsg("zero0");
			gc.probif += rgs.sample1D(B,A,t2-t1,inf_node[n].inf_ev[0].tdiv-t2,warn,gc,rs);
			
			auto prob = poisson_probability(B,mut_rate*(t2-t1));
			gc.dlike_genetic_process += prob;
		}
		else{	
			//cout << "add2" << endl;
			auto num = inf_ev_add[j_add].mut_num;	
			auto t3 = inf_ev_add[j_add].tdiv;
			
			if(type == REGRAFT_1D_SAMPLER){
				gc.probif += rgs.sample1D(B,A,t2-t1,t3-t2,warn,gc,rs);
			}
			else{
				gc.probif += rgs.sample2D(B,A,gc,rs);
			}
			if(B > num) return false;
			
			gc.dlike_genetic_process += poisson_probability(B,mut_rate*(t2-t1)) 
																+ poisson_probability(num-B,mut_rate*(t3-t2))
																- poisson_probability(num,mut_rate*(t3-t1));                   
		}
	}
	else{
		//cout << "add3" << endl;
		auto seq_var = genetic_value.seq_var;
		auto lam = seq_var+mut_rate*t;
	
		if(type != REGRAFT_1D_SAMPLER) emsg("Should be 1d");
		if(mut_rate*(inf_node[n].inf_ev[0].tdiv-t) == 0) emsg("zero2");
		gc.probif += rgs.sample1D(B,A,lam,mut_rate*(inf_node[n].inf_ev[0].tdiv-t),warn,gc,rs);
		
		auto prob = poisson_probability(B,lam);
		gc.dlike_genetic_process += prob;
	}
	
	gc.mut_num = B;
	gc.mutA_num =	A;	
	
	return true;
}


/// Removing node from branch
bool State::change_remove_node_with_gobs(unsigned int n, double t, const NodeRef &nr_from, TransTreeChange &gc, const GenObsFunc &gof, RegraftStore rs) const 		
{	
	const auto &inf_node = genetic_value.inf_node;
	auto mut_rate = genetic_value.mut_rate*model.details.dt;

	auto rgs = get_regraft_sampler(gof);
	if(rgs.ele.size() == 0) return false;
	
	auto A = inf_node[n].inf_ev[0].mut_num;
	
	auto err = false;
	
	auto type = rgs.type;
	
	if(nr_from.node != OUTSIDE_INF && nr_from.node != ENTER_INF){ // Removing current node position
		const auto &in_from = inf_node[nr_from.node];
		const auto &inf_ev_from = in_from.inf_ev;
		auto j_from = nr_from.index;
		
		double t1; 
		if(j_from > 0) t1 = inf_ev_from[j_from-1].tdiv; 
		else t1 = in_from.tdiv_start;
		auto t2 = inf_ev_from[j_from].tdiv;
		
		auto B = inf_ev_from[j_from].mut_num;
		
		if(j_from+1 == inf_ev_from.size()){            // End event
			//cout << "rem1" << endl;
			gc.probfi += rgs.prob1D(B,A,t2-t1,inf_node[n].inf_ev[0].tdiv-t2,err,gc,rs);
			if(err) return false;
			
			auto prob = poisson_probability(B,mut_rate*(t2-t1));
			gc.dlike_genetic_process -= prob; 
		}
		else{
			//cout << "rem2" << endl; 
			auto t3 = inf_ev_from[j_from+1].tdiv;
			
			auto C = inf_ev_from[j_from+1].mut_num;

			if(type == REGRAFT_1D_SAMPLER){
				gc.probfi += rgs.prob1D(B,A,t2-t1,t3-t2,err,gc,rs);
			}
			else{
				gc.probfi += rgs.prob2D(B,A,err,gc,rs);
			}
			if(err) return false;
		
			gc.dlike_genetic_process += poisson_probability(B+C,mut_rate*(t3-t1)) 
																- poisson_probability(B,mut_rate*(t2-t1)) 
																- poisson_probability(C,mut_rate*(t3-t2)); 						
		}
	}
	else{
		//cout << "rem3" << endl;
		auto B = genetic_value.inf_origin[nr_from.index].mut_num;
		
		auto seq_var = genetic_value.seq_var;
	
		auto lam = seq_var+mut_rate*t;
		
		if(type != REGRAFT_1D_SAMPLER) emsg("Should be 1d");
		
		gc.probfi += rgs.prob1D(B,A,lam,mut_rate*(inf_node[n].inf_ev[0].tdiv-t),err,gc,rs);
		if(err) return false;
		
		auto prob = poisson_probability(B,lam);
	
		gc.dlike_genetic_process -= prob;
	}
	
	return true;
}


/// Splits up using the binomial distribution
void State::binomial_split(int &B0, int &A0, double r1, double r2, string &warn) const
{
	auto N = B0;
	
	if(N == 0){ 
		B0 = 0; A0 = 0;
		return;
	}
	
	if(r1 <= 0 || r2 <= 0) emsg("neg num");
	auto pr = r1/(r1+r2);
	B0 = binomial_sample(pr,N,warn);
	A0 = N-B0;
}


/// Adding node to new branch and maximising observation probability
void State::change_add_node_gen_obs_max(unsigned int n, double t, const NodeRef &nr_add, TransTreeChange &gc, const GenObsFunc &gof) const 
{
	const auto &inf_node = genetic_value.inf_node;
	auto mut_rate = genetic_value.mut_rate*model.details.dt;
	
	int B0, A0;
	get_regraft_minimum(B0,A0,gof);
	
	auto fl_2d = true;  // Determines if two dimensional
	if(A0 == UNSET) fl_2d = false;
	
	if(B0 < 0) B0 = 0;
	
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
			if(fl_2d) emsg("should not be 2d");
			binomial_split(B0,A0,t2-t1,inf_node[n].inf_ev[0].tdiv-t2,warn);
			
			auto prob = poisson_probability(B0,mut_rate*(t2-t1));
			gc.dlike_genetic_process += prob;
		}
		else{	
			auto num = inf_ev_add[j_add].mut_num;	
			auto t3 = inf_ev_add[j_add].tdiv;
			
			if(!fl_2d){
				binomial_split(B0,A0,t2-t1,t3-t2,warn);		
			}
			
			if(B0 > (int)num) B0 = num;
			if(A0 < 0) A0 = 0;

			gc.dlike_genetic_process += poisson_probability(B0,mut_rate*(t2-t1)) 
																+ poisson_probability(num-B0,mut_rate*(t3-t2))
																- poisson_probability(num,mut_rate*(t3-t1));                   
		}
	}
	else{
		auto seq_var = genetic_value.seq_var;
		auto lam = seq_var+mut_rate*t;
		
		if(fl_2d) emsg("should not be 2d");
		binomial_split(B0,A0,lam,mut_rate*(inf_node[n].inf_ev[0].tdiv-t),warn);
		
		auto prob = poisson_probability(B0,lam);
		gc.dlike_genetic_process += prob;
	}
	
	gc.mut_num = B0;
	gc.mutA_num =	A0;	
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


/// Outputs who acquired infection from whom (for diagnostic purposes)
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
	
	const auto &individual = species[pop.p].individual;
	
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
void TransTreeChange::update_like_ch(Like &like_ch, double &dprob)
{
	if(type != NO_TRANS_TREE_CHANGE){
		dprob += probfi-probif;
									
		like_ch.genetic_process += dlike_genetic_process;
		like_ch.genetic_obs += dlike_genetic_obs; 
		
		dlike_genetic_process = 0;
		dlike_genetic_obs = 0;
	}
}


/// Accepts a change to a gentic difference
void State::accept_gen_dif_cha(const vector <GenDifCha> &gen_dif_cha)
{
	auto &gen_dif = genetic_value.gen_dif;
		
	for(const auto &cha : gen_dif_cha){
		auto ob1 = cha.ob1;
		auto ch = cha.ch;
		if(ch != 0){
			for(auto ob2 : cha.ob2){
				gen_dif[ob1][ob2] += ch;
				gen_dif[ob2][ob1] += ch;
			}
		}
	}
}

				
/// If the update is accepted then this updates state
void State::gen_change_update(const TransTreeChange &gc)
{
	if(gc.type == NO_TRANS_TREE_CHANGE) return;
	
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
	in.ti_start = na.ti_start;
	in.ti_rec = na.ti_rec;
	
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
	
	auto &gen_dif = genetic_value.gen_dif;

	switch(gc.type){
	case REGRAFT_NODE: case REGRAFT_GEN_OBS_MAX:
		switch(gc.regraft_type){ // The reconnection position is the same place on the tree
		case NORMAL_REGRAFT: case REGRAFT_WITH_OBS: case REGRAFT_WITH_OBS2: case REGRAFT_WITH_OBS3:
			disconnect_inf_node(n);
			connect_inf_node(n,gc.iif_add,gc.mut_num);
			
			if(gc.regraft_type != NORMAL_REGRAFT){ 
				like.genetic_obs += gc.dlike_genetic_obs;
				 
				accept_gen_dif_cha(gc.gen_dif_cha);
			}
			
			if(gc.mutA_num != UNSET){
				inf_node[n].inf_ev[0].mut_num = gc.mutA_num;
			}			
			break;
			
			
		case ORIGIN_REGRAFT:
			{
				auto &nr_from = in.from;
				if(species[in.p].individual[in.i].ev[in.e].type == ENTER_EV) nr_from.node = ENTER_INF;
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
			
		case REGRAFT_NOT_NEEDED:
			emsg("SHouold not be REGRAFT_NOT_NEEDED");
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
		
	case NO_TRANS_TREE_CHANGE: break;
		
	case TRANS_TREE_FAIL: emsg("Should not be here20"); break;
	}
			
	// If new observations have appeared then recalculate genetic difference and likelihood
	if(na.obs_add_end.size() > 0 || na.obs_add_begin.size() > 0){
		gen_dif = calculate_gen_dif();
		like.genetic_obs = likelihood_genetic_obs(gen_dif);
		like.genetic_process = likelihood_genetic_process();
	}
	else{	
		like.genetic_process += gc.dlike_genetic_process;
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
		//popnum_ind[ti].resize(model.pop.size());
		popnum_ind[ti].resize(model.npop);
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
			auto pp = model.pop[pir.po].p;
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


/// Estimates expected steady state value for genetic obs 
double State::steadystate_genetic_obs()
{
	return 0.5*model.genetic_data.obs.size();
}


/// The contribution coming from a harmonic oscillator with va certain variance
double State::harmomic_contribution(double var) const 
{
	auto val = (1.0/var)*exp(-1.0/(2*var));
	if(val > 0.5 || var > 1) val = 0.5;
	return val;
}


/// Calculates the expected value for genetic_obs
double State::calc_genetic_obs_expected()
{
	if(!model.genetic_data.on) return 0;
		
	const auto &inf_node = genetic_value.inf_node;
	const auto &gen_data = model.genetic_data;
	
	if(true){ // This is just a simple upper limit
		auto sum = 0.0;
		for(auto n = 0u; n < inf_node.size(); n++){
			sum += inf_node[n].inf_ev.size();
		}
		return -sum/2;
	}
	else{
		auto nobs = genetic_value.obs_node_ref.size();
		
		auto c = gen_data.obs_scale;
		
		vector < vector <double> > cont_max;
		cont_max.resize(inf_node.size());
		for(auto n = 0u; n < inf_node.size(); n++){
			cont_max[n].resize(inf_node[n].inf_ev.size(),0);
		}
		
		auto nums = 0u;
		for(auto n = 0u; n < inf_node.size(); n++){
			auto &no = inf_node[n];
			for(auto e = 0u; e < no.inf_ev.size(); e++){
				// Mutation update
				{
					auto obs_node = get_gen_obs_list(n,e);	
			
					auto obs_other = get_not_obs(obs_node,vector <unsigned int> (),nobs);

					auto var = 1.0/(4*c*obs_node.size()*obs_other.size());
					auto val = harmomic_contribution(var);
					
					if(val > cont_max[n][e]) cont_max[n][e] = val;
				}
				
				// OBS_GEN_UPDATE
				if(no.inf_ev[e].type == GENETIC_OBS && e+1 < no.inf_ev.size()){
					auto num = no.inf_ev[e].mut_num + no.inf_ev[e+1].mut_num;
					
					if(num > 0){
						auto var = 1.0/(c*nobs);
						auto val = 0.5*harmomic_contribution(var);
						if(val > cont_max[n][e]) cont_max[n][e] = val;
						if(val > cont_max[n][e+1]) cont_max[n][e+1] = val;
						nums++;
					}
				}
				
				// Node possibilities
				{
					const auto &ev = no.inf_ev[e];
					if(ev.type == INFECT_OTHER){
						auto n2 = ev.index;
						const	auto &in2 = inf_node[n2];
								
						auto nbranch = 0u;
						if(e+1 < no.inf_ev.size()) nbranch++;
						if(in2.inf_ev.size() > 0) nbranch++;
						
						if(nbranch == 2){					
							auto obs_node1 = get_gen_obs_list(n,e+1);
							auto obs_node2 = get_gen_obs_list(ev.index,0);
							auto NB = obs_node1.size();
							auto NC = obs_node2.size();
							auto NA = nobs-NB-NC;
							
							auto num = no.inf_ev[e].mut_num + no.inf_ev[e+1].mut_num + in2.inf_ev[0].mut_num;
							if(num > 0){
								{
									auto var = 1.0/(8*c*NB*NC);
									auto val = (1.0/3)*harmomic_contribution(var);
									if(val > cont_max[n][e]) cont_max[n][e] = val;
									if(val > cont_max[n][e+1]) cont_max[n][e+1] = val;
									if(val > cont_max[n2][0]) cont_max[n2][0] = val;
								}
								
								{
									auto var = 1.0/(8*c*NC*NA);
									auto val = (1.0/3)*harmomic_contribution(var);
									if(val > cont_max[n][e]) cont_max[n][e] = val;
									if(val > cont_max[n][e+1]) cont_max[n][e+1] = val;
									if(val > cont_max[n2][0]) cont_max[n2][0] = val;
								}
								
								{
									auto var = 1.0/(8*c*NA*NB);
									auto val = (1.0/3)*harmomic_contribution(var);
									if(val > cont_max[n][e]) cont_max[n][e] = val;
									if(val > cont_max[n][e+1]) cont_max[n][e+1] = val;
									if(val > cont_max[n2][0]) cont_max[n2][0] = val;
								}
								
								{
									auto var = 1.0/(2*c*NA*(nobs-NA));
									auto val = (1.0/2)*harmomic_contribution(var);
									if(val > cont_max[n][e+1]) cont_max[n][e+1] = val;
									if(val > cont_max[n2][0]) cont_max[n2][0] = val;
								}
								
								{
									auto var = 1.0/(2*c*NB*(nobs-NB));
									auto val = (1.0/2)*harmomic_contribution(var);
									if(val > cont_max[n][e]) cont_max[n][e] = val;
									if(val > cont_max[n2][0]) cont_max[n2][0] = val;
								}
								
								{
									auto var = 1.0/(2*c*NC*(nobs-NC));
									auto val = (1.0/2)*harmomic_contribution(var);
									if(val > cont_max[n][e]) cont_max[n][e] = val;
									if(val > cont_max[n][e+1]) cont_max[n][e+1] = val;
								}
							}
						}					
					}
				}		
			}
		}

		auto sum = 0.0;
		auto ma = 0.0;
		for(auto n = 0u; n < cont_max.size(); n++){
			for(auto e = 0u; e < cont_max[n].size(); e++){
				sum += cont_max[n][e];
				ma -= 0.5;
			}
		}

		return -sum;
	}
}


/// Tests that the function get_genetic_obs_func_from is working correctly
void State::test_get_genetic_obs_func_from()
{
	cout << "test_get_genetic_obs_func_from" << endl;
	auto &inf_node = genetic_value.inf_node;
	
	int sh = 30;
	
	for(auto n = 0u; n < inf_node.size(); n++){
		cout << n << " n" << endl;
		auto obs_node = get_gen_obs_list(n,0); 
						
		if(obs_node.size() > 0){ 
			auto gof_from = get_genetic_obs_func_from(obs_node,n); 
		
			const auto &in = inf_node[n];
			const auto &nr_from = in.from;
			auto mutA_num = in.inf_ev[0].mut_num;
			
			auto mutC = UNSET;
			auto indC = nr_from.index+1;
			
			unsigned int num;
			if(nr_from.node > CODE) num = genetic_value.inf_origin[nr_from.index].mut_num;
			else{
				num = inf_node[nr_from.node].inf_ev[nr_from.index].mut_num;
				
				if(indC < inf_node[nr_from.node].inf_ev.size()){
					mutC = inf_node[nr_from.node].inf_ev[indC].mut_num;
				}		
			}
	
			auto gobs_st = like.genetic_obs;
			
			for(int mutA = (int)mutA_num-sh; mutA < (int)mutA_num+sh; mutA++){
				for(int mutB = num-sh; mutB < (int)num+sh; mutB++){
					int mutCC = UNSET; if(mutC != UNSET) mutCC = (int)num+(int)mutC-(int)mutB;
					
					if(mutA >= 0 && mutB >= 0 && (mutCC == UNSET || mutCC >= 0)){
						if(nr_from.node > CODE) genetic_value.inf_origin[nr_from.index].mut_num = mutB;
						else inf_node[nr_from.node].inf_ev[nr_from.index].mut_num = mutB;
				
						inf_node[n].inf_ev[0].mut_num = mutA;
				
						if(mutC != UNSET) inf_node[nr_from.node].inf_ev[indC].mut_num = num+mutC -mutB;
					
						genetic_value.gen_dif = calculate_gen_dif();	
						like.genetic_obs = likelihood_genetic_obs(genetic_value.gen_dif);

						auto est = calc_genetic_obs_fast(mutB,mutA,gof_from);
						auto act = like.genetic_obs - gobs_st;
						
						if(dif(est,act,0.000000001)){
							cout <<  mutA << " " << mutA_num <<"   " << mutB << " " << num << " " << est << " " << act << "   mut" << endl;
							emsg("Different");
						}
					}
				}
			}
			
			if(nr_from.node > CODE) genetic_value.inf_origin[nr_from.index].mut_num = num;
			else inf_node[nr_from.node].inf_ev[nr_from.index].mut_num = num;
			inf_node[n].inf_ev[0].mut_num = mutA_num;
				
			if(mutC != UNSET) inf_node[nr_from.node].inf_ev[indC].mut_num = mutC;
			
			genetic_value.gen_dif = calculate_gen_dif();	
			like.genetic_obs = likelihood_genetic_obs(genetic_value.gen_dif);
		}
	}
}


/// Determines if maximise_gen_obs_off is turned off (this is used to help maximise the genetic observation likelihood)
void State::check_maximise_gen_obs_off(bool on) 
{
	auto &mgi = genetic_value.maximise_gen_obs;
	if(!mgi) return;
	
	if(!on || !model.genetic_data.on){ mgi = false; return;}

	auto exp = calc_genetic_obs_expected();
	if(like.genetic_obs > exp) mgi = false;
}


/// Checks genetic values have been correctly updated
void State::check_gen_values(string te) const
{
	auto gen_dif = calculate_gen_dif();
	
	if(difi(genetic_value.gen_dif,gen_dif)) emsg("gen dif prob:"+te);
	
	auto obs = likelihood_genetic_obs(genetic_value.gen_dif); 
	if(dif(like.genetic_obs,obs,SMALL)){
		cout << like.genetic_obs << " " << obs << "obs" << endl;
		emsg("gobs wrong:"+te);	
	}
}


/// Checks that regraft is correctly working
void State::check_regraft()
{		
	auto pl = false;
	
	if(like.genetic_obs > -0.1 && false){
			print_matrix("M",model.genetic_data.gen_dif);
		print_matrix("GD",genetic_value.gen_dif);

	const auto &ssp = species[0];		
		for(auto i = 0u; i < ssp.individual.size(); i++){
			const auto &ind = ssp.individual[i];
			cout << ind.name << " ";
			for(auto e =0u; e < ind.ev.size(); e++){
				auto &iif = ind.ev[e].ind_inf_from;
				if(iif.i != UNSET){
					cout << " from " << ssp.individual[iif.i].name;
				}
			}
			cout << "  infect" << endl;
		}
		
		for(auto ob = 0u; ob < model.genetic_data.obs.size(); ob++){
			auto obs = model.genetic_data.obs[ob];
	
			cout << ob << " " << obs.name << " " << ssp.individual[obs.i].name << endl; 
		}
		print_all_node();
		
		auto gen_dif = calculate_gen_dif();
		print_matrix("gd",gen_dif);
	
		auto nobs = model.genetic_data.obs.size();
		for(auto m = 0u; m < nobs; m++){
			const auto &nr = genetic_value.obs_node_ref[m];
			cout << nr.node << " " << nr.index << " node ref" << endl;
		}
		
		emsg("don");
	}
	
	if(ran() > 0.1) return;
	
	if(pl || ran() < 0.1) cout << "check regraft" << endl;

	auto ch_fl = false;

	auto change_time = false;
	if(ran() < 0.5) change_time = true;

	const auto &inf_node = genetic_value.inf_node;

	for(auto loop = 0u; loop < 3; loop++){
		switch(loop){
		case 0: regraft_obs_type = RG_OBS; break;
		case 1: regraft_obs_type = RG_OBS2; break;
		case 2: regraft_obs_type = RG_OBS3; break;
		}
		
		for(auto n = 0u; n < inf_node.size(); n++){
			const auto &in = inf_node[n];
			
			auto nn = (unsigned int)(ran()*inf_node.size());
			
			if(nn != n){
				auto t = in.tdiv_start;
				auto tstart = t;
				
				const auto &in_new = inf_node[nn];
				
				auto &ev = species[in.p].individual[in.i].ev;
				auto &eve = ev[in.e]; 
				
				if(eve.type != ENTER_EV){
					if(eve.tr_gl == UNSET) emsg("prob");
					const auto &sp = model.species[in.p];

					if(pl) cout << change_time << " change_time" << endl;
					
					if(change_time){
						auto tmi = 0.0;
						auto tma = in.tdiv_rec;
						if(in.inf_ev.size() > 0) tma = in.inf_ev[0].tdiv;
						if(in.e+1 < ev.size()){
							if(ev[in.e+1].tdiv < tma) tma = ev[in.e+1].tdiv;
						}
						if(in.e > 0){
							if(ev[in.e-1].tdiv > tmi) tmi = ev[in.e-1].tdiv;
						}
						
						t = tmi + ran()*(tma-tmi);
						
						if(pl) cout << tstart << " -> " << t << " tt" << endl;
					}
		
					IndInfFrom iif_prop;
					if(set_iif(in_new.p,in_new.i,in_new.e,iif_prop,t,eve,sp)){
						auto alt = no_alter(n);
						auto alt_back = alt;
						
						if(change_time){
							alt.tdiv_start = t;
							alt.ti_start = inf_node_ti_from_tdiv(alt.tdiv_start);
						}
						
						if(pl) species[0].print_event(ev);
						//print_node("node",n);
						
						if(pl) cout << endl << endl << "start" << endl;
						
						auto gc = trans_tree_change(REGRAFT_NODE,n,alt,iif_prop,REGRAFT_STORE);
			
						if(gc.type != TRANS_TREE_FAIL){					
							auto iif = eve.ind_inf_from;
						
							gen_change_update(gc);
					
							eve.tdiv = t;
					
							if(ch_fl) check_gen_values("mid");
							
							if(pl) cout << endl << "rev" << endl;
							
							alt_back.rev_store = gc.rev_store;
				
							auto gc2 = trans_tree_change(REGRAFT_NODE,n,alt_back,iif,REGRAFT_RETRIEVE);
							if(gc2.type != gc.type){
								emsg("Type does not agree");
							}
							if(gc2.regraft_type != gc.regraft_type) emsg("regraft_type does not agree");
						
							if(dif(gc.probif,gc2.probfi,TINY)){
								cout << gc.probif << " " << gc2.probfi << " hh" << endl;
								emsg("Difference probif1");
							}
							
							if(dif(gc.probfi,gc2.probif,TINY)){
								cout << gc.probfi << " " << gc2.probif << "chfi" << endl;
								emsg("Difference probif2");
							}
							
							if(dif(gc.dlike_genetic_process,-gc2.dlike_genetic_process,TINY)){
								cout << gc.dlike_genetic_process << " " << -gc2.dlike_genetic_process <<" val" << endl;
								emsg("dlike_genetic_process");
							}				
						
							if(dif(gc.dlike_genetic_obs,-gc2.dlike_genetic_obs,TINY)) emsg("dlike_genetic_obs");
							
							gen_change_update(gc2);
							
							eve.tdiv = tstart;
							
							if(ch_fl) check_gen_values("end");
							
							//check("OK");
						}
					}
				}
			}
		}
	}
	
	if(pl) cout << "done" << endl;
	
	regraft_obs_type = RG_NOT_SET;
}


/// Gets the prediction accuracy of identifying the correct infector
double State::get_infector_prediction_accuracy() const
{
	auto num = 0.0;
	for(const auto &da : model.genetic_data.waifw_data){
		const auto &ind = species[da.to_p].individual[da.to_i];
		for(const auto &ev : ind.ev){
			const auto &iif = ev.ind_inf_from;
			if(iif.p == da.from_p && iif.i == da.from_i){ num++; break;}
		}
	}
	
	return num/model.genetic_data.waifw_data.size();
}

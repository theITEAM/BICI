// Functions for updating the model based on command lines in the input

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <algorithm>

using namespace std;

#include "input.hh"
#include "utils.hh"

/// Adds a new species to the model
void Input::add_species(string name, SpeciesType type, bool trans_tree)
{
	check_name_input(name,"Species name",true);
	
	auto p = find_p(name);
	if(p != UNSET){ alert_import("There is already a species with the name '"+name+"'"); return;}
	
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		
		auto cl = find_cl(p,name);
		if(cl != UNSET){ alert_import("Cannot have the same name as classification '"+sp.cla[cl].name+"' "); return;}
			
		cl = find_cl_index(p,name);
		if(cl != UNSET){ alert_import("Cannot have the same name as index for classification '"+sp.cla[cl].name+"' "); return;}
	
		for(auto cl = 0u; cl < sp.ncla; cl++){
			auto c = find_c(p,cl,name);
			if(c != UNSET){
				alert_import("Cannot have the same name as compartment '"+sp.cla[cl].comp[c].name+"' in classification '"+sp.cla[cl].name+"'");
				return;
			}
		}
	}
	
	Species sp(model.species.size(),model.timepoint,model.details); 
	sp.name = name; sp.ncla = 0; sp.type = type; sp.trans_tree = trans_tree; sp.infection_cl = UNSET;
	sp.line_num = line_num;
	
	model.species.push_back(sp);	
	model.nspecies = model.species.size();  
}


/// Adds a new classification to the model
void Input::add_classification(unsigned int p, string name, string index, Coord coord)
{
	check_name_input(name,"Classification name",true);

	if(p == UNSET){ 
		alert_import("Cannot add a classification when the species is not set"); 
		return;
	}
	
	auto cl = find_cl(p,name);
	if(cl != UNSET){ alert_import("There is already a classification with the name '"+name+"'"); return;}
	
	if(name == index){ alert_import("Name and index must be different"); return;}
	
	Classification clas; clas.name = name; clas.index = index; clas.ntra = 0; clas.ncomp = 0; clas.coord = coord;
	
	auto &sp = model.species[p];
	
	sp.cla.push_back(clas);
	sp.ncla++;
}


/// Adds a new compartment to the model
void Input::add_compartment(string name, unsigned int p, unsigned int cl, double x, double y, double lat, double lng, bool markov_branch, CompInfected infected, string erlang_source)
{
	if(check_comp_exist(name,p) == true) return;
	
	auto &sp = model.species[p];
	auto &claa = sp.cla[cl];
	
	Compartment co; 
	co.name = name; 
	co.x = x; co.y = y; co.lat = lat; co.lng = lng; 
	co.markov_branch = markov_branch;
	co.all_branches = false;
	co.infected = infected;
	co.erlang_hidden = false;
	co.erlang_c_start = UNSET;
	if(erlang_source != ""){
		co.erlang_hidden = true;
		auto cc = 0u; while(cc < claa.comp.size() && claa.comp[cc].name != erlang_source) cc++;
		if(cc == claa.comp.size()) emsg_input("Erlang source problem");
		co.erlang_c_start = cc;
	}	
	co.erlang_source = erlang_source;
	
	if(infected == COMP_INFECTED){
		if(sp.infection_cl != cl && sp.infection_cl != UNSET) alert_import("Cannot sect 'infected' for more than one classification"); 
		sp.infection_cl = cl;
	}
	
	//auto vec = claa.hash_comp.get_vec_string(co.name);
	claa.hash_comp.add(claa.ncomp,co.name);

	claa.comp.push_back(co);
	claa.ncomp++;
}
	

/// Adds a transition to the model
void Input::add_transition(unsigned int p, unsigned int cl, unsigned int i, unsigned int f, TransType type)
{
	auto &sp = model.species[p];
	auto &claa = sp.cla[cl];
	auto &tra = claa.tra;
	
	for(auto j = 0u; j < tra.size(); j++){
		if(tra[j].i == i && tra[j].f == f){
			if(i == UNSET){ alert_import("Cannot have two sources into the same compartment"); return;}
			else{
				if(f == UNSET){ alert_import("Cannot have two sinks leaving the same compartment"); return;}
				else{
					alert_import("Cannot have two transitions with the same source and destination compartments");
					return;
				}
			}
		}
	}
	
	if(i == f){ alert_import("The 'from' and 'to' compartments must be different"); return;}
		
	Transition tr;
	tr.i = i; tr.f = f;
	tr.type = type; 
	tr.bp_set = BP_UNSET;
	tr.branch = false;
	tr.line_num = line_num;
	tr.erlang_hidden = false;

	if(i == UNSET){ tr.variety = SOURCE_TRANS; tr.name = "+→"+claa.comp[f].name;}
	else{
		if(f == UNSET){ tr.variety = SINK_TRANS; tr.name = claa.comp[i].name+"→-";}
		else{ tr.variety = NORMAL; tr.name = claa.comp[i].name+"→"+claa.comp[f].name;}
	}
	
	if(sp.type == POPULATION){
		string type = "";
		switch(tr.type){
		case GAMMA: type = "gamma"; break;
		case LOG_NORMAL: type = "log-normal"; break;
		case WEIBULL: type = "Weibull"; break;
		case PERIOD: type = "period"; break;
		default: break;
		}
		
		if(type != ""){
			alert_import("Transition '"+tr.name+"' cannot have a '"+type+"' distribution because '"+sp.name+"' uses a population-based model."); 
		}
	}
	
	tra.push_back(tr);	
	claa.ntra++;
}
	

/// Determines if compartments/transition are branching
void Input::determine_branching() const 
{
	for(auto p = 0u; p < model.nspecies; p++){
		auto &sp = model.species[p];
		for(auto cl = 0u; cl < sp.ncla; cl++){
			auto &claa = sp.cla[cl];
			for(auto tr = 0u; tr < claa.ntra; tr++) claa.tra[tr].branch = false;

			for(auto c = 0u; c < claa.ncomp; c++){
				auto &co = claa.comp[c];
				
				vector <unsigned int> list; 
				auto flag = false;
				for(auto tr = 0u; tr < claa.ntra; tr++){
					const auto &trr = claa.tra[tr];
					if(trr.i == c){
						list.push_back(tr);
						if(trr.type != EXP_RATE && trr.type != EXP_MEAN) flag = true;
					}
				}
				
				co.branching = false; if(list.size() > 1) co.branching = true;
				
				if(list.size() > 1 && (flag == true || co.markov_branch == true)) co.branch = true;
				else co.branch = false;
				
				for(auto i = 0u; i < list.size(); i++){
					claa.tra[list[i]].branch = co.branch;
				}
			}					
		}
	}
}
	
	
/// Clones a classification from another species
void Input::clone_class(unsigned int p_to, unsigned int p_from, unsigned int cl_from)
{
	auto cla_from = model.species[p_from].cla[cl_from];
	
	add_classification(p_to,cla_from.name,cla_from.index,cla_from.coord);

	auto cl = model.species[p_to].ncla-1;

	auto &claa = model.species[p_to].cla[cl];

	claa.ncomp = cla_from.ncomp;
	claa.comp = cla_from.comp;
	claa.hash_comp = cla_from.hash_comp;
}


/// Creates equations based on EquationInfo specifications
void Input::create_equations(unsigned int per_start, unsigned int per_end)
{
	model.create_species_simp();

	Hash hash_eqn;
	
	for(auto &der : model.derive){                       // Derived quantities
		for(auto &eq : der.eq) model.add_eq_ref(eq,hash_eqn);
	}

	auto dper = double(per_end-per_start)/model.species.size();
	for(auto p = 0u; p < model.species.size(); p++){    // Transitions
		auto &sp = model.species[p];  
		
		auto per_st = per_start+dper*p;
		auto fr = dper/sp.tra_gl.size();
		for(auto i = 0u; i < sp.tra_gl.size(); i++){
			auto &tr_gl = sp.tra_gl[i];
				
			if(i%100 == 0) percentage(per_st+fr*i,100);
			if(tr_gl.branch == true && tr_gl.bp_set == BP_SET) model.add_eq_ref(tr_gl.bp,hash_eqn);
			for(auto &dp : tr_gl.dist_param){
				model.add_eq_ref(dp,hash_eqn);
			}
		}
	}
	
	for(auto &sp : model.species){                       // Initial conditions
		auto &ic = sp.init_cond;
		if(ic.type == INIT_POP_DIST){
			if(ic.focal_cl == UNSET){
				for(auto &eqi : ic.pop_prior.dist_param){
					model.add_eq_ref(eqi,hash_eqn);
				}	
			}
			else{
				for(auto c = 0u; c < sp.cla[ic.focal_cl].ncomp; c++){
					for(auto &eqi : ic.comp_prior[c].dist_param){
						model.add_eq_ref(eqi,hash_eqn);
					}	
				}
			}  
		}
	}
	
	for(auto &sp : model.species){                       // Observation model
		for(auto &pf : sp.pop_filter){                     // Population filter   
			pf.time_vari = false;
			for(auto &cpe : pf.comp_prob_eqn){		
				model.add_eq_ref(cpe,hash_eqn);
				const auto &eq = model.eqn[cpe.eq_ref];
				if(eq.time_vari == true) pf.time_vari = true;
			}
			if(pf.time_vari == false){
				pf.comp_obs_mod_ref = sp.obs_eqn_add_vec(pf.comp_prob_eqn);
			}
			pf.eqn_zero = set_eqn_zero(pf.comp_prob_eqn);
			pf.c_nonzero = false_list(pf.eqn_zero);
		}
		
		for(auto &pd : sp.pop_data){                       // Population data (time vari obs mod)
			auto &pf = sp.pop_filter[pd.ref];
			pd.time_vari = pf.time_vari;
			if(pd.time_vari){
				for(auto &cpe : pf.comp_prob_eqn) model.add_eq_ref(cpe,hash_eqn,pd.t);
				pd.comp_obs_mod_ref = sp.obs_eqn_add_vec(pf.comp_prob_eqn);
			}
		}
		
		for(auto &ptf : sp.pop_trans_filter){              // Population transition filter   
			ptf.time_vari = false;
			for(auto &tpe : ptf.trans_prob_eqn){		
				model.add_eq_ref(tpe,hash_eqn);
				const auto &eq = model.eqn[tpe.eq_ref];
				if(eq.time_vari == true) ptf.time_vari = true;
			}
			if(ptf.time_vari == false){
				ptf.trans_obs_mod_ref = sp.obs_eqn_add_vec(ptf.trans_prob_eqn);
			}
			ptf.eqn_zero = set_eqn_zero(ptf.trans_prob_eqn);
			ptf.tr_nonzero = false_list(ptf.eqn_zero);
		}
		
		for(auto &ptd : sp.pop_trans_data){                // Population transition data
			auto &ptf = sp.pop_trans_filter[ptd.ref];
			ptd.time_vari = ptf.time_vari;
			if(ptd.time_vari){
				for(auto &tpe : ptf.trans_prob_eqn){		
					model.add_eq_ref(tpe,hash_eqn,0.5*(ptd.tmin+ptd.tmax));
				}
				ptd.trans_obs_mod_ref = sp.obs_eqn_add_vec(ptf.trans_prob_eqn);
			}
		}
		
		for(auto &ecp : sp.enter){                         // Equations for entering probability
			if(ecp.c_set == UNSET){
				for(auto cl = 0u; cl < sp.ncla; cl++){
					auto &ecp_cl = ecp.cla[cl];
					if(ecp_cl.c_set == UNSET){
						for(auto &cop : ecp_cl.eqn) model.add_eq_ref(cop,hash_eqn,ecp.time);					
						ecp_cl.obs_eqn_ref = sp.obs_eqn_add_vec(ecp_cl.eqn);
					}
				}
			}
		}
		
		// Deals with individual observation probabilities with obs model time independant
		for(auto &ot : sp.obs_trans){                      // Observed transition data 
			ot.time_vari = false;
			for(auto &tpe : ot.tra_prob_eqn){
				model.add_eq_ref(tpe,hash_eqn);
				
				const auto &eq = model.eqn[tpe.eq_ref];
				if(eq.time_vari == true) ot.time_vari = true;
			}
		
			if(ot.time_vari == false){
				ot.obs_eqn_ref = sp.obs_eqn_add_vec(ot.tra_prob_eqn);
			}
		}
	
		for(auto &ind : sp.individual){                    // Individual data 
			for(auto &ob : ind.obs){
				switch(ob.type){
				case OBS_TRANS_EV: case OBS_SOURCE_EV: case OBS_SINK_EV:  // Ind. transition
					{
						if(ob.ref < sp.obs_trans.size()){
							auto &ot = sp.obs_trans[ob.ref];
							ob.time_vari = ot.time_vari;
							
							if(ob.time_vari == true){
								auto tra_prob_eqn = ot.tra_prob_eqn;
								for(auto &tpe : tra_prob_eqn){
									model.add_eq_ref(tpe,hash_eqn,ob.t);
								}
								
								ob.obs_eqn_ref = sp.obs_eqn_add_vec(tra_prob_eqn);
							}
						}
					}
					break;
					
				case OBS_COMP_EV:                              // Ind. Compartment
					if(ob.c_exact == UNSET){
						for(auto &cop : ob.c_obs_prob_eqn){
							model.add_eq_ref(cop,hash_eqn,ob.t);
						}
				
						ob.obs_eqn_ref = sp.obs_eqn_add_vec(ob.c_obs_prob_eqn);
						ob.eqn_zero = set_eqn_zero(ob.c_obs_prob_eqn);
					}
					break;
					
				case OBS_TEST_EV:                              // Ind. test
					model.add_eq_ref(ob.Se_eqn,hash_eqn,ob.t);
					ob.Se_obs_eqn_ref = add_to_vec(sp.obs_eqn,ob.Se_eqn.eq_ref);
					
					model.add_eq_ref(ob.Sp_eqn,hash_eqn,ob.t);
					ob.Sp_obs_eqn_ref = add_to_vec(sp.obs_eqn,ob.Sp_eqn.eq_ref);
					break;
				}
			}
		}
	}

	if(model.genetic_data.on){                  // Genetic data
		model.add_eq_ref(model.genetic_data.mut_rate,hash_eqn);
		model.add_eq_ref(model.genetic_data.seq_var,hash_eqn);
	}

	for(auto p = 0u; p < model.species.size(); p++){
		const auto &sp = model.species[p];
		if(sp.type == POPULATION){
			if(sp.ind_effect.size() > 0){
				const auto &ie = sp.ind_effect[0];
				alert_line("A population-based model cannot have individual effects such as '["+ie.name+"]'",ie.line_num); 
			}
			
			if(sp.fix_effect.size() > 0){
				const auto &fe = sp.fix_effect[0];
				alert_line("A population-based model cannot have fixed effects such as '<"+fe.name+">'",fe.line_num); 
			}
		}
	}
	
	// Adds in any parameter reparameterisation 
	// This is done at the end so it is only aplied to any "used" parameter  
	bool flag;
	do{
		flag = false;
		auto num = 0u;
		for(auto th = 0u; th < model.param.size(); th++){
			auto &par = model.param[th];
		
			if(par.variety != CONST_PARAM){
				switch(par.variety){
				case DIST_PARAM: case PRIOR_PARAM:
					for(auto &ele : par.element){
						auto &pri = model.prior[ele.prior_ref];
						for(auto &eqi : pri.dist_param){
							if(eqi.eq_ref == UNSET){
								flag = true;
								num++;
								model.add_eq_ref(eqi,hash_eqn);
							}
						}
					}
					break;
				
				case REPARAM_PARAM:
					if(par.reparam_eqn != ""){   // Reparameterisation set by equation
						add_reparam_eqn(par,hash_eqn); 
						flag = true;
					}
					else{                        // Reparametersisation set individually
						for(auto &ele : par.element){
							auto &eqi = ele.value;	
							if(eqi.eq_ref == UNSET){
								flag = true;
								num++;
								model.add_eq_ref(eqi,hash_eqn);
							}
						}
					}
					break;
					
				case CONST_PARAM: emsg_input("Should not be const"); break;
					
				case UNSET_PARAM: emsg_input("Should not be unset7"); break;
				}
			}
		}
	}while(flag == true);
	
	for(auto eq : model.eqn){
		if(eq.warn != ""){
			stringstream ss; ss << endl << "For equation '" << eq.te_raw << "':" << endl << " " << eq.warn;
			alert_line(ss.str(),eq.line_num);
		}
	}
}


/// Tries to linearise equations in terms of populations
void Input::linearise_eqn(unsigned int per_start, unsigned int per_end)
{
	auto fl = false;
	if(model.mode != INF){
		fl = true;
		for(const auto &sp : model.species){
			if(sp.trans_tree) fl = false;
		}
	}

	if(fl == true) return;
	
	auto dper = double(per_end-per_start)/model.species.size();
	
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		
		auto per_st = per_start+dper*p;
		auto fr = dper/sp.tra_gl.size();
		for(auto e = 0u; e < sp.markov_eqn.size(); e++){
			const auto &me = sp.markov_eqn[e];
		
			if(e%100 == 0) percentage(per_st+fr*e,100);
			model.eqn[me.eqn_ref].calculate_linearise();    
		}			
	}
	
	//for(auto &eq : model.eqn){
	//eq.calculate_linearise();        // Tries to linearise equations in terms of populations
	//}
}


/// Sets a vector which determines if an equation is zero
vector <bool> Input::set_eqn_zero(const vector <EquationInfo> &eq_info)
{
	vector <bool> zero;
	for(auto i = 0u; i < eq_info.size(); i++){
		const auto &eq = model.eqn[eq_info[i].eq_ref];
		zero.push_back(eq.is_zero());
	}
	
	return zero;
}


/// Adds a parent child relationship if there is a reparameterisation or distribution
void Input::add_parent_child(const EquationInfo eqi, unsigned int i, unsigned int th, Hash &hash)
{
	ParamRef parref; parref.th = th; parref.index = i;
					
	auto &par = model.param[th];
	
	auto ref = eqi.eq_ref;
	if(ref != UNSET){
		const auto &eq = model.eqn[ref];

		for(auto pr : eq.param_ref){
			vector <unsigned int> vec;
			vec.push_back(th); vec.push_back(i);
			vec.push_back(pr.th); vec.push_back(pr.index);
		
			auto p = hash.existing(vec);
			if(p == UNSET){
				hash.add(0,vec);
		
				par.add_parent(i,pr);
				model.param[pr.th].add_child(pr.index,parref);
			}
			
			//add_to_list(par.parent[i],pr);
			//add_to_list(model.param[pr.th].child[pr.index],parref);
		}
	}
}

						
/// Works out which equations are Markovian
void Input::create_markov_eqn()
{
	for(auto e = 0u; e < model.eqn.size(); e++){
		auto &eqn = model.eqn[e];
	
		switch(eqn.type){
		case TRANS_RATE: case TRANS_MEAN: case SOURCE_RATE: case SOURCE_MEAN:
			if(eqn.sp_p == UNSET) alert_import("Equation problem"); 
			else{
				auto &sp = model.species[eqn.sp_p];
				
				eqn.markov_eqn_ref = sp.markov_eqn.size();
				
				MarkovEqn me; 
				me.eqn_ref = e; 
				me.rate = true;	if(eqn.type == TRANS_MEAN || eqn.type == SOURCE_MEAN) me.rate = false;
				me.time_vari = eqn.time_vari;
				me.infection_trans = eqn.infection_trans;
				
				const auto &eq = model.eqn[e];	
				me.source = false; 
				if(eq.type == SOURCE_RATE || eq.type == SOURCE_MEAN){
					me.source = true;	
					me.source_tr_gl = eq.source_tr_gl;
				}
			
				me.ind_eff_mult = eq.ind_eff_mult;
				me.fix_eff_mult = eq.fix_eff_mult;
				me.ind_variation = false;
				if(me.ind_eff_mult.size() > 0 || me.fix_eff_mult.size() > 0) me.ind_variation = true;
		
				sp.markov_eqn.push_back(me);
			}
			break;
			
		default: break;
		}
	}
	
	// Add a reference from global transitions to markovian equations
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
		
		for(auto &tr : sp.tra_gl){
			tr.markov_eqn_ref = UNSET;
			if(tr.type == EXP_RATE || tr.type == EXP_MEAN){
				auto e = tr.dist_param[0].eq_ref; if(e == UNSET) emsg_input("Prob1");
				const auto &eq = model.eqn[e];
				auto mef = eq.markov_eqn_ref; if(mef == UNSET) emsg_input("Prob2");
				tr.markov_eqn_ref = mef;
			}
		}
	}

	if(false){
		for(const auto &sp : model.species){ 
			for(const auto &me : sp.markov_eqn){
				cout << model.eqn[me.eqn_ref].te << "  infection trans:";
				if(me.infection_trans) cout << "true"; else cout << "false";
				cout << endl;
			}
		}
	}
}


/// Works out which populations affect which markov eqns (and vice versa)
void Input::create_markov_eqn_pop_ref()
{
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
				
		switch(sp.type){
		case INDIVIDUAL:
			{
				
				for(auto e = 0u; e < sp.markov_eqn.size(); e++){
					auto eq = sp.markov_eqn[e].eqn_ref;
					if(eq == UNSET) emsg_input("Problem with ME");
				
					for(auto k : model.eqn[eq].pop_ref){
						PopMarkovEqnRef mer; mer.p = p; mer.e = e;
						model.pop[k].markov_eqn_ref.push_back(mer);
					}
				}
				
			}
			break;
			
		case POPULATION:
			{
				for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
					const auto &tra = sp.tra_gl[tr];
					
					vector <unsigned int> list;
					
					for(auto &ei : tra.dist_param) list.push_back(ei.eq_ref);
					if(tra.bp.eq_ref != UNSET) list.push_back(tra.bp.eq_ref);
					
					for(auto eq : list){
						for(auto k : model.eqn[eq].pop_ref){
							PopTransRef mtr; mtr.p = p; mtr.tr = tr;
							model.pop[k].trans_ref.push_back(mtr);
						}
					}
				}
			}
			break;
		}
	}		
}


/// Intialises global compartments and transitions
void Input::global_comp_trans_init()
{
	for(auto &sp : model.species){
		auto N = 1u;                                   // Works out global conversion
		for(auto cl = 0u; cl < sp.ncla; cl++){
			sp.comp_mult.push_back(N);
			N *= sp.cla[cl].comp.size();
		}
		
		sp.N = N;
	
		for(auto c = 0u; c < N; c++){                  // Constructs global compartments
			CompGlobal co;
			co.erlang_hidden = false;
			co.erlang_c_start = UNSET;
			co.infected = false;
			
			string na = "";
			for(auto cl = 0u; cl < sp.ncla; cl++){
				auto &claa = sp.cla[cl];
				auto cc = (c/sp.comp_mult[cl])%claa.comp.size();
				co.cla_comp.push_back(cc);
				
				const auto &comp = claa.comp[cc];
				if(comp.erlang_hidden == true) co.erlang_hidden = true;
				
				if(comp.infected == COMP_INFECTED) co.infected = true;
					
				if(na != "") na += "|";
				na += comp.name;
				
				CompGlTransGroup ti; ti.branch = comp.branch; 
				ti.all_branches = comp.all_branches;
				ti.markov = UNSET;
				co.tra_leave_group.push_back(ti);
			}
			co.name = na;
			
			sp.comp_gl.push_back(co);
		}

		for(auto cl = 0u; cl < sp.ncla; cl++){         // Constructs global transitions
			const auto &cla = sp.cla[cl];
			
			vector <DepConv> dep_conv;
			for(auto cl2 = 0u; cl2 < sp.ncla; cl2++){
				if(cl2 != cl){
					DepConv dc; 
					dc.before = sp.cla[cl2].index;
					dep_conv.push_back(dc);
				}
			}
						
			for(auto j = 0u; j < cla.ntra; j++){
				const auto &tr = cla.tra[j];

				auto bp_swap_temp = swap_template(tr.bp.te,dep_conv);
				if(bp_swap_temp.warn != ""){ alert_equation(tr.bp,bp_swap_temp.warn); return;}
	
				vector <SwapResult> dist_swap_temp;
				for(auto k = 0u; k < tr.dist_param.size(); k++){
					auto swap_temp = swap_template(tr.dist_param[k].te,dep_conv);
					if(swap_temp.warn != ""){ alert_equation(tr.dist_param[k],swap_temp.warn); return;
					}
					dist_swap_temp.push_back(swap_temp);
				}
				
				for(auto c = 0u; c < N; c++){
					const auto &co = sp.comp_gl[c];		

					auto cc = co.cla_comp[cl];
					auto ii = UNSET, ff = UNSET;

					if(tr.i == UNSET){  
						if(tr.f == cc){ ii = UNSET; ff = c;}
					}
					else{
						if(tr.f == UNSET){
							if(tr.i == cc){ ii = c; ff = UNSET;}
						}
						else{
							if(tr.i == cc){
								int dif = (int)tr.f - (int)tr.i;
								ii = c; ff = c + dif*sp.comp_mult[cl];
							}
						}
					}

					if(!(ii == UNSET && ff == UNSET)){
						string name;
						for(auto cl2 = 0u; cl2 < sp.ncla; cl2++){
							const auto &claa2 = sp.cla[cl2];
							if(name != "") name += "|";
							if(cl2 == cl) name += tr.name;
							else name += claa2.comp[co.cla_comp[cl2]].name;
						}
			
						TransGlobal tr_gl;
						tr_gl.name = name;
						tr_gl.cl = cl;
						tr_gl.tr = j;
						tr_gl.i = ii;
						tr_gl.f = ff;
						tr_gl.erlang_hidden = tr.erlang_hidden;
						tr_gl.type = tr.type;
						tr_gl.ev_type = M_TRANS_EV;
						if(tr.type != EXP_RATE && tr.type != EXP_MEAN) tr_gl.ev_type = NM_TRANS_EV;
						tr_gl.variety = tr.variety;
						tr_gl.branch = tr.branch;
						tr_gl.all_branches = tr.all_branches;
						tr_gl.markov_eqn_ref = UNSET;
						tr_gl.nm_trans_ref = UNSET;
							
						tr_gl.infection.type = TRANS_INF_UNCHANGE;
						if(sp.trans_tree == true){
							if((ii == UNSET || sp.comp_gl[ii].infected == false) && (ff != UNSET && sp.comp_gl[ff].infected == true)){
								tr_gl.infection.type = TRANS_INFECTION;
							}
							
							if((ii != UNSET && sp.comp_gl[ii].infected == true) && (ff == UNSET || sp.comp_gl[ff].infected == false)){
								tr_gl.infection.type = TRANS_RECOVERY;
							}
						}
	
						tr_gl.bp = tr.bp;
						tr_gl.bp_set = tr.bp_set;
						tr_gl.dist_param = tr.dist_param;
						tr_gl.line_num = tr.line_num;
						
						auto kk = 0u;
						for(auto cl2 = 0u; cl2 < sp.ncla; cl2++){
							if(cl2 != cl){
								dep_conv[kk].after = sp.cla[cl2].comp[co.cla_comp[cl2]].name;
								kk++;
							}
						}
	
						if(tr_gl.infection.type == TRANS_INFECTION){
							for(auto &dp : tr_gl.dist_param) dp.infection_trans = true;
						}
						
						tr_gl.bp.c = c;
						
						tr_gl.bp.te = swap_index_temp(dep_conv,bp_swap_temp);
						if(check_swap){
							auto te_ch = tr.bp.te;
							auto res = swap_index(te_ch,dep_conv); 
							if(res.warn != ""){ alert_equation(tr_gl.bp,res.warn); return;}
							if(te_ch != tr_gl.bp.te) emsg_input("Swap not agree");
						}
						
						for(auto i = 0u; i < tr_gl.dist_param.size(); i++){
							tr_gl.dist_param[i].c = c;
							
							tr_gl.dist_param[i].te = swap_index_temp(dep_conv,dist_swap_temp[i]);
							if(check_swap){
								auto te_ch = tr.dist_param[i].te;
								auto res = swap_index(te_ch,dep_conv); 	
								if(res.warn != ""){ alert_equation(tr_gl.dist_param[i],res.warn); return;}
								if(te_ch != tr_gl.dist_param[i].te) emsg_input("Swap not agree");
							}
						}

						if(ii != UNSET){
							auto &ti = sp.comp_gl[ii].tra_leave_group[cl];
							if(ti.branch != tr.branch) emsg_input("Branch does not agree");
					
							ti.tr_list.push_back(sp.tra_gl.size());
						}
					
						sp.tra_gl.push_back(tr_gl);
					}
				}
			}
		}
		
		// Sets Markov property of global compartment
		for(auto &cgl : sp.comp_gl){
			for(auto &tlg : cgl.tra_leave_group){
				tlg.markov = true;
				if(tlg.branch == true) tlg.markov = false;
				for(auto tr : tlg.tr_list){
					if(sp.tra_gl[tr].type != EXP_RATE && sp.tra_gl[tr].type != EXP_MEAN) tlg.markov = false;
				}
			}
		}
		
		// Sets erlang_c_start
		for(auto c = 0u; c < sp.comp_gl.size(); c++){
			auto &co = sp.comp_gl[c];
			if(co.erlang_hidden == true){
				auto vec = co.cla_comp;
				for(auto cl = 0u; cl < sp.ncla; cl++){
					const auto &claa = sp.cla[cl];
					
					const auto &com = claa.comp[vec[cl]];
					if(com.erlang_hidden == true) vec[cl] = com.erlang_c_start;
				}
				
				auto c_new = 0u; for(auto cl = 0u; cl < sp.ncla; cl++) c_new += sp.comp_mult[cl]*vec[cl];
				co.erlang_c_start = c_new;
			}
		}
	}
	
	// Calculates "tform" 
	// This allows for transtions o be changed based on a transition in another cl
	// Either the initial of final compartment can be used for the transformation
	// If memory requirement is too high, this is not used
	for(auto &sp : model.species){
		sp.tform_set = false;
		if(sp.type == INDIVIDUAL){
			if(1 == 0 && sp.tra_gl.size()*sp.comp_gl.size() < TFORM_MAX){ // TURN OFF
				sp.tform_set = true;
				
				for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
					sp.tra_gl[tr].tform.resize(sp.comp_gl.size(),UNSET);
				}
				
				for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
					auto &t1 = sp.tra_gl[tr];
					for(auto tr2 = tr; tr2 < sp.tra_gl.size(); tr2++){
						auto &t2 = sp.tra_gl[tr2];
						if(t1.cl == t2.cl && t1.tr == t2.tr){
							if(t2.i != UNSET) t1.tform[t2.i] = tr2;
							if(t2.f != UNSET) t1.tform[t2.f] = tr2;
							if(t1.i != UNSET) t2.tform[t1.i] = tr;
							if(t1.f != UNSET) t2.tform[t1.f] = tr;
						}
					}
				}
				
				if(false){
					for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
						auto &t1 = sp.tra_gl[tr];
						cout << t1.name << ":" << endl;
						for(auto c = 0u; c < sp.comp_gl.size(); c++){
							if(t1.tform[c] != UNSET){
								auto &t2 = sp.tra_gl[t1.tform[c]];
								cout << t2.name << ",";
							}
						}
						cout << endl;
					}
				}
				
				if(false){
					for(auto c = 0u; c < sp.comp_gl.size(); c++){
						const auto &cgl = sp.comp_gl[c];
						cout << cgl.name << " name" << endl;
						for(auto tr : cgl.tr_leave){
							cout << sp.tra_gl[tr].name << ",";
						}
						cout << "trr" << endl;
					}
				}
			}
			else{
				sp.tr_shift.resize(sp.ncla);
				for(auto cl = 0u; cl < sp.ncla; cl++){
					sp.tr_shift[cl].resize(sp.ncla);
					auto mult = 1u;
					for(auto cl2 = 0u; cl2 < sp.ncla; cl2++){
						if(cl2 == cl) sp.tr_shift[cl][cl2] = UNSET;
						else{
							sp.tr_shift[cl][cl2] = mult;
							mult *= sp.cla[cl2].ncomp;
						}
					}
				}
			}
		}
	}
}
	
	
/// Sets tr_leave_markov from global compartments
void Input::set_tr_leave_markov()
{
	for(auto &sp : model.species){
		for(auto c = 0u; c < sp.comp_gl.size(); c++){
			auto &cgl = sp.comp_gl[c];
			for(const auto &tlg : cgl.tra_leave_group){
				for(auto tr : tlg.tr_list){
					cgl.tr_leave.push_back(tr);
					if(sp.tra_gl[tr].markov_eqn_ref != UNSET){
						cgl.tr_leave_markov.push_back(tr);
					}
				}
			}
		}
	}
}


/// Sets if transitions in multiple classification
void Input::set_multi_trans_cl()
{
	for(auto &sp : model.species){
		auto num = 0u;
		for(const auto &claa : sp.cla){
			if(claa.ntra > 0) num++;
		}
		if(num > 1) sp.multi_trans_cl = true;
		else sp.multi_trans_cl = false;
	}
}


/// Determines if global transitions are connected
void Input::set_tr_connected()
{
	for(auto &sp : model.species){
		auto N = sp.tra_gl.size();		
		
		// Determines if global transitions are connected [tr_from][tr_to]
		vector < vector <bool> > tr_con; 
		
		//auto &tr_con = sp.tr_connected;
		auto &tr_aft = sp.tr_after;
		
		// First gets all direct transitions
		tr_con.resize(N);
		tr_aft.resize(N);
		for(auto tr = 0u; tr < N; tr++){
			tr_con[tr].resize(N,false);
			
			const auto &tra = sp.tra_gl[tr];
			auto f = tra.f;
			if(f != UNSET){
				for(auto trf : sp.comp_gl[f].tr_leave){
					if(tr_con[tr][trf] == false){
						tr_con[tr][trf] = true;
						tr_aft[tr].push_back(trf);
					}
				}
			}
		}	
		
		auto tr_aft_dir = tr_aft;
		
		vector <unsigned int> tr_aft_index(N,0);
		
		// Expands to get all subsequent transitions
		bool fl;
		do{
			fl = false;
			for(auto tr = 0u; tr < N; tr++){
				while(tr_aft_index[tr] < tr_aft[tr].size()){
					auto k = tr_aft_index[tr];
				
					auto trf = tr_aft[tr][k];
			
					for(auto trf2 : tr_aft_dir[trf]){
						if(tr_con[tr][trf2] == false && trf2 != trf){
							tr_con[tr][trf2] = true;
							tr_aft[tr].push_back(trf2);
							fl = true;
						}
					}
					
					tr_aft_index[tr]++;
				}
			}
		}while(fl == true);
			
		// Gets transitions before
		auto &tr_before = sp.tr_before;
		tr_before.resize(N);
		for(auto tr = 0u; tr < N; tr++){
			for(auto trf : tr_aft[tr]){
				tr_before[trf].push_back(tr);
			}
		}
	}

	if(false){
		cout << "Connected" << endl;
		for(const auto &sp : model.species){
			const auto &tr_aft = sp.tr_after;
		
			auto N = sp.tra_gl.size();		
			for(auto tr = 0u; tr < N; tr++){
				cout << sp.tra_gl[tr].name << ": ";
				for(auto tr2 : tr_aft[tr]){
					cout << sp.tra_gl[tr2].name << ", ";
				}
				cout << "trans after    ";
				
				for(auto tr2 : sp.tr_before[tr]){
					cout << sp.tra_gl[tr2].name << ", ";
				}
				cout << "trans before" << endl;
			}				
		}
		emsg_input("done");
	}
}
	

/// Determines if nm_rate can be pre-calculated
void Input::set_precalc_nm_rate()
{
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
	
		for(auto m = 0u; m < sp.nm_trans.size(); m++){
			auto &nmt = sp.nm_trans[m];

			auto precalc = true;
			
			{		
				const auto &eq = model.eqn[nmt.dist_param_eq_ref[0]];
				for(auto pr : eq.pop_ref){
					if(model.pop[pr].sp_p == p){ 
						alert_import("Non-Markovian transition cannot depend on a population in the same species");
					}
				}
			}
			
			auto bp_eq = nmt.bp_eq;
			if(bp_eq != UNSET){
				if(bp_eq == BP_FROM_OTHERS){
					for(auto e : nmt.bp_other_eq){
						const auto &eq = model.eqn[e];
						
						for(auto pr : eq.pop_ref){
							if(model.pop[pr].sp_p == p){ 
								alert_import("Non-Markovian transition cannot depend on a population in the same species");
							}
						}
						
						if(eq.ind_eff_mult.size() > 0 || eq.fix_eff_mult.size() > 0){
							precalc = false;
						}
					}
				}
				else{
					const auto &eq = model.eqn[bp_eq];	
					if(nmt.all_branches){
						if(eq.pop_ref.size() > 0){
							precalc = false;
						}
					}
					else{
						if(eq.pop_ref.size() > 0 || eq.ind_eff_mult.size() > 0 
							|| eq.fix_eff_mult.size() > 0){
							precalc = false;
						}
					}
				}
			}
			
			nmt.precalc_nm_rate = precalc;
		}
	}
}


/// Checks that transition are either markovian or non-markovian
void Input::check_markov_or_nm() const
{
	for(const auto &sp : model.species){
		for(const auto &tr : sp.tra_gl){
			auto num = 0u;
			if(tr.nm_trans_ref != UNSET) num++;
			if(tr.markov_eqn_ref != UNSET) num++;

			if(num != 1) emsg_input("Transition should be markovian or nm");
		}
	}
}


/// Checks no population in species for nm trans (otherwise simulation proposals not correct)
void Input::check_nm_pop()
{
	for(auto p = 0u; p < model.species.size(); p++){
		const auto &sp = model.species[p];
	
		for(const auto &tr : sp.tra_gl){
			if(tr.type != EXP_RATE && tr.type != EXP_MEAN){
				for(const auto &dp : tr.dist_param){
					const auto &eq = model.eqn[dp.eq_ref];
					for(auto &pop : eq.pop_ref){
						const auto &po = model.pop[pop];
						if(po.sp_p == p){
							alert_import("Non-Markovian parameter cannot depend on population from the same species.");
						}
					}
				}
				
				if(tr.bp.eq_ref != UNSET){
					const auto &eq = model.eqn[tr.bp.eq_ref];
					for(auto &pop : eq.pop_ref){
						const auto &po = model.pop[pop];
						if(po.sp_p == p){
							alert_import("Branching probability cannot depend on population from the same species.");
						}
					}
				}
			}
		}
	}
}


/// Divides an Erlang distribution into multiple compartments
void Input::create_population_erlang()
{
	for(auto p = 0u; p < model.nspecies; p++){
		auto &sp = model.species[p];
		if(sp.type == POPULATION){
			for(auto cl = 0u; cl < sp.ncla; cl++){
				auto &claa = sp.cla[cl];
				auto tr = 0u;
				while(tr < claa.ntra){
					auto &tra = claa.tra[tr];
					if(tra.type == ERLANG){
						auto ci = tra.i;
						auto cf = tra.f;
						auto cname = claa.comp[ci].name;
						auto infected = claa.comp[ci].infected;
						const auto &mean_eqn = tra.dist_param[0];
						const auto &shape_eqn = tra.dist_param[1];
						auto shape = number(shape_eqn.te);
						auto shapei = (unsigned int)shape;
						auto tname = tra.name;
						auto lnum = tra.line_num;
						
						if(shape == UNSET || shape <= 0 || shapei != shape){
							alert_equation(shape_eqn,"Shape not a positive integer"); return;
						}
						
						auto rate_new = tstr(shapei)+"/("+mean_eqn.te_raw+")";
						for(auto j = 1u; j <= shapei; j++){
							unsigned int cnew;
							if(j == shapei) cnew = cf;
							else{
								cnew = claa.ncomp;
								auto name = cname+"_";
								name += tstr(j);
								name += "_erlang";
							
								add_compartment(name,p,cl,UNSET,UNSET,UNSET,UNSET,false,infected,cname);
							}
							
							add_transition(p,cl,ci,cnew,EXP_RATE);
							
							auto &tra2 = claa.tra[claa.ntra-1];
							tra2.dist_param.push_back(he(add_equation_info(rate_new,TRANS_RATE,p,cl)));
							tra2.line_num = lnum;
							if(j < shapei) tra2.erlang_hidden = true;
							else tra2.name = tname;
							
							string bp = "";
							if(bp != "") tra2.bp_set = BP_SET;
							if(bp == "") bp = "b";
							tra2.bp = he(add_equation_info(bp,BP,p,cl));
							
							ci = cnew;
						}
						
						claa.tra.erase(claa.tra.begin()+tr);
						claa.ntra--;
					}
					else tr++;
				}				
			}
		}	
	}
}


/// Attaches compartments to source equations
void Input::source_equation_comp()
{
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
			const auto &tra = sp.tra_gl[tr];
			if(tra.type == EXP_RATE || tra.type == EXP_MEAN){
				auto &eq = model.eqn[tra.dist_param[0].eq_ref];
				if(eq.type == SOURCE_RATE || eq.type == SOURCE_MEAN){
					eq.source_tr_gl.push_back(tr);
				}
			}
		}
	}	
}


/// Calculates the timepoints used to solve the equations
void Input::calculate_timepoint()
{
	const auto &de = model.details;
	
	auto t = de.t_start;
	model.timepoint.push_back(t);
	t += de.dt+TINY;

	auto t_end = de.t_end;
	if(model.mode == PPC) t_end = de.ppc_t_end;
	
	while(t < t_end-TINY){
		model.timepoint.push_back(t);
		t += de.dt;
	}
	model.timepoint.push_back(t_end);
	model.ntimepoint = model.timepoint.size();
	
	if(false){
		for(auto i = 0u; i < model.timepoint.size(); i++){
			cout << i << " "<<  model.timepoint[i] << " timepoint" << endl;
		}
	}
}


/// Used for converting splines to paramvec
unsigned int Input::get_spline_i(const ParamRef &pr)
{					
	const auto &spl = model.spline;		
	
	auto th = pr.th, index = pr.index;
	auto ntimes = model.param[th].spline_info.knot_times.size();
	if(index%ntimes != 0) emsg_input("Should be zero1");
	index /= ntimes;
	auto i = 0u; while(i < spl.size() && !(spl[i].th == th && spl[i].index == index)) i++;
	if(i == spl.size()) emsg_input("Spline could not be found");
	
	return i;
}

					
/// Creates splices for the relevant parameters
void Input::create_spline()
{
	const auto &tp = model.timepoint;
	auto ntp = tp.size();
	
	for(auto th = 0u; th < model.param.size(); th++){
		const auto &par = model.param[th];
		
		if(par.time_dep == true){
			auto &sinfo = par.spline_info;
			if(sinfo.on != true) emsg_input("Problem loading spline");
			
			const auto &times = sinfo.knot_times;
			auto ntimes = times.size();
			
			auto dep_reduce = par.dep;
			dep_reduce.pop_back();
			
			auto ndep = dep_reduce.size();
			for(auto &dep : dep_reduce) dep.mult /= ntimes;
				
			auto N = 1u;
			for(auto d = 0u; d < ndep; d++) N *= dep_reduce[d].list.size();

			for(auto j = 0u; j < N; j++){
				Spline spl; 
				spl.name = get_param_name_with_dep(par,dep_reduce,j);
				spl.th = th; 
				spl.index = j;
				spl.constant = false;
				spl.info = par.spline_info;
		
				vector <unsigned int> param_ref;
			
				if(par.variety == CONST_PARAM) spl.constant = true;
				else{	
					for(auto i = 0u; i < ntimes; i++){
						auto k = 0u;
						while(k < model.param_vec.size()){
							const auto &pv = model.param_vec[k];
							
							if(pv.th == th && pv.index == j*ntimes + i) break;
							k++;
						}
						if(k == model.param_vec.size()) emsg_input("Could not find");
						
						param_ref.push_back(k);
					}
					
					spl.param_ref = param_ref;
				}
				
				auto i = 0u;                             // Indexes times
				for(auto ti = 0u; ti < ntp-1; ti++){
					if(tp[ti] < times[0]){ // Extends parameter factors to start time
						spl.const_val.push_back(1);
						if(model.mode != PPC) emsg_input("Should be PPC");
						if(!begin_str(spl.name,"f~")) emsg_input("Should be a parameter factor");
					}
					else{
						auto tmid = (tp[ti]+tp[ti+1])/2;
						while(i+1 < times.size() && times[i+1] < tmid) i++;
						if(i+1 >= times.size()){ // This extends spline for PPC
							if(model.mode != PPC) emsg_input("Should be PPC");
							
							if(spl.constant == false){
								SplineDiv sd;	sd.th1 = param_ref[i]; sd.th2 = param_ref[i]; sd.f = 1;	
								spl.div.push_back(sd);
							}
							else{
								spl.const_val.push_back(par.get_value(j*ntimes+i));
							}
						}
						else{
							auto f = (times[i+1]-tmid)/(times[i+1]-times[i]);
							if(spl.constant == false){
								SplineDiv sd;	sd.th1 = param_ref[i]; sd.th2 = param_ref[i+1]; sd.f = f;	
								spl.div.push_back(sd);
							}
							else{
								auto val1 = par.get_value(j*ntimes+i);
								auto val2 = par.get_value(j*ntimes+i+1);
								spl.const_val.push_back(val1*f+val2*(1-f));
							}
						}
					}
				}
					
				model.spline.push_back(spl);
 			}
		}			
	}
	
	// Converts equations such that the spline is referenced	
	for(auto &eq : model.eqn){
		for(auto &ca : eq.calc){
			for(auto &it : ca.item){
				if(it.type == SPLINE){
					const auto &pr = eq.param_ref[it.num];
					it.type = SPLINEREF;
					it.num = get_spline_i(pr);
				}
				
				if(it.type == SPLINE_TIME){
					auto &tr = eq.time_ref[it.num];
					const auto &pr = eq.param_ref[tr.num];
					it.type = SPLINEREF_TIME;
					tr.num = get_spline_i(pr);
				}
			}
		}
	
		if(eq.ans.type == SPLINE){
			const auto &pr = eq.param_ref[eq.ans.num];
			eq.ans.type = SPLINEREF;
			eq.ans.num = get_spline_i(pr);
		}
		
		if(eq.ans.type == SPLINE_TIME){
			auto &tr = eq.time_ref[eq.ans.num];
			const auto &pr = eq.param_ref[tr.num];
			eq.ans.type = SPLINEREF_TIME;
			tr.num = get_spline_i(pr); 
		}
	}
	
	if(false){
		for(auto i = 0u; i < model.spline.size(); i++){
			auto spl = model.spline[i];
			
			cout << spl.name << "name" << endl;
			cout << spl.constant << "const" << endl;
			cout << model.param[spl.th].name << " param" << endl;
			cout << spl.index << " index" << endl;
			for(auto j : spl.param_ref) cout << j << ","; 
			cout << "param_ref" << endl;
			
			for(auto d = 0u; d < spl.div.size(); d++){
				const auto &di = spl.div[d];
				cout << d << " " << model.param_vec[di.th1].name << " " <<  model.param_vec[di.th2].name << " " << di.f << " div"<< endl;
			}
		}
		emsg_input("Done");
	}
}


/// Creates a list of non-Markovian transitions
void Input::create_nm_trans()
{			
	for(auto &sp : model.species){                 // Sets nm_trans
		sp.nm_trans_incomp_info();                   // Sets information about incomplete transitions
		
		auto &nm = sp.nm_trans;
	
		for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
			auto &trg = sp.tra_gl[tr];
		
			auto type = trg.type;
		
			if(type == EXP_RATE || type == EXP_MEAN){
				trg.nm_trans_ref = UNSET;
			}
			else{	
				NMTrans nmtra;
				nmtra.name = trg.name;
				nmtra.type = type;
				nmtra.bp_eq = trg.bp.eq_ref; 
				nmtra.all_branches = trg.all_branches;
				nmtra.trans_incomp_ref = UNSET;
				
				if(nmtra.bp_eq == UNSET){
					const auto &trr = sp.cla[trg.cl].tra[trg.tr];
				
					if(trr.bp_set == BP_DERIVED){
						nmtra.bp_eq = BP_FROM_OTHERS;
					}
				}
				
				// Gets other transitions
				auto ci = trg.i;
				auto cl = trg.cl;
				
				const auto &tr_list = sp.comp_gl[ci].tra_leave_group[cl].tr_list;
				
				for(auto tr : tr_list){
					auto e = sp.tra_gl[tr].bp.eq_ref;
					if(nmtra.all_branches){
						nmtra.bp_all_eq.push_back(e);
					}
					
					if(nmtra.bp_eq == BP_FROM_OTHERS && e != BP_FROM_OTHERS){
						nmtra.bp_all_eq.push_back(e);
					}
				}
				
				for(auto i = 0u; i < trg.dist_param.size(); i++){
					nmtra.dist_param_eq_ref.push_back(trg.dist_param[i].eq_ref);
				}
				
				// Finds how individual factor alter the rate (used in individual sampling)
				switch(type){
				case EXP_RATE: case EXP_MEAN: emsg_input("Should not be here"); break;
				
				case GAMMA: case ERLANG: case LOG_NORMAL: case PERIOD: case WEIBULL:
				case EXP_MEAN_NM:
					{
						const auto &eq = model.eqn[nmtra.dist_param_eq_ref[0]];
						
						for(auto e : eq.ind_eff_mult){
							IndFacRate ifr; ifr.type = IND_EFF_DIV; ifr.e = e;
							nmtra.ind_fac_rate.push_back(ifr);
						}
						
						for(auto e : eq.fix_eff_mult){
							IndFacRate ifr; ifr.type = FIX_EFF_DIV; ifr.e = e;
							nmtra.ind_fac_rate.push_back(ifr);
						}
					}
					break;
				
				case EXP_RATE_NM:
					{
						const auto &eq = model.eqn[nmtra.dist_param_eq_ref[0]];
						
						for(auto e : eq.ind_eff_mult){
							IndFacRate ifr; ifr.type = IND_EFF_MULT; ifr.e = e;
							nmtra.ind_fac_rate.push_back(ifr);
						}
						
						for(auto e : eq.fix_eff_mult){
							IndFacRate ifr; ifr.type = FIX_EFF_MULT; ifr.e = e;
							nmtra.ind_fac_rate.push_back(ifr);
						}
					}
					break;
				}
				
				auto k = 0u;
				while(k < nm.size() && 
				!(nm[k].type == nmtra.type && 
				nm[k].bp_eq == nmtra.bp_eq &&
				equal_vec(nm[k].bp_all_eq,nmtra.bp_all_eq) == true &&
				equal_vec(nm[k].bp_other_eq,nmtra.bp_other_eq) == true &&
				equal_vec(nm[k].dist_param_eq_ref,nmtra.dist_param_eq_ref) == true)) k++;
				
				if(k == nm.size()) nm.push_back(nmtra);
				
				trg.nm_trans_ref = k;
			}
		}
	}

	for(auto &sp : model.species){  // Sets nm_trans_incomp
		for(auto c = 0u; c < sp.comp_gl.size(); c++){
			auto &cgl = sp.comp_gl[c];
		
			cgl.nmtransincomp_ref.resize(sp.ncla);
			for(auto cl = 0u; cl < sp.ncla; cl++){
				const auto &tlg = cgl.tra_leave_group[cl];
				if(tlg.markov == false){
					NMTransIncomp nmtrin;
					nmtrin.name = sp.comp_gl[c].name;
					nmtrin.cl = cl;
					nmtrin.all_branches = false;
					for(auto tr : tlg.tr_list){
						const auto &tra = sp.tra_gl[tr];
						if(tra.all_branches) nmtrin.all_branches = true;
						auto m = tra.nm_trans_ref;	
						if(m == UNSET) emsg_input("nm_trans_ref should be set");
						nmtrin.nmtrans_ref.push_back(m);
					}
		
					auto n = 0u;
					while(n < sp.nm_trans_incomp.size() && 
					  !equal_vec(sp.nm_trans_incomp[n].nmtrans_ref,nmtrin.nmtrans_ref)) n++;
					
					cgl.nmtransincomp_ref[cl] = n;
					
					if(n == sp.nm_trans_incomp.size()){	
						for(auto m : nmtrin.nmtrans_ref){
							sp.nm_trans[m].trans_incomp_ref = n;
						}							
					
						sp.nm_trans_incomp.push_back(nmtrin);
					}
				}
			}	
		}
		
		for(const auto &nmt : sp.nm_trans){
			if(nmt.trans_incomp_ref == UNSET) emsg_input("trans_incomp_ref not set");
		}
	}
	
	if(false){
		for(auto &sp : model.species){              // Sets nm_trans
			cout << "nm trans" << endl;
			for(auto &nm : sp.nm_trans){
				cout << nm.name << "  trans" << endl;
			}
		}
		
		for(auto &sp : model.species){               // Sets nm_trans
			cout << "nm trans incomp" << endl;
			for(auto &nm : sp.nm_trans_incomp){
				cout << nm.name << "  trans incomp" << endl;
			}
		}
		
		emsg_input("done");
	}
}


/// Sets a matrix which determines which ie and fe affect a given transition
void Input::set_tra_ie_fe_affect()
{
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
		
		auto N = sp.tra_gl.size();
		auto E = sp.ind_effect.size();
		auto F = sp.fix_effect.size();
		
		sp.tra_ie_affect.resize(N);
		sp.tra_fe_affect.resize(N);
		for(auto tr = 0u; tr < N; tr++){
			sp.tra_ie_affect[tr].resize(E,false);
			sp.tra_fe_affect[tr].resize(F,false);
		}
		
		for(auto e = 0u; e < E; e++){
			vector <bool> me(sp.markov_eqn.size(),false);
			vector <bool> nm(sp.nm_trans.size(),false);
			
			const auto &ie = sp.ind_effect[e];
			for(auto i : ie.markov_eqn_ref) me[i] = true;
			for(auto i : ie.nm_trans_ref) nm[i] = true;
			
			for(auto tr = 0u; tr < N; tr++){
				const auto &tra = sp.tra_gl[tr];
				if(tra.markov_eqn_ref != UNSET && me[tra.markov_eqn_ref]) sp.tra_ie_affect[tr][e] = true;
				if(tra.nm_trans_ref != UNSET && nm[tra.nm_trans_ref]) sp.tra_ie_affect[tr][e] = true;
			}
		}
		
		for(auto f = 0u; f < F; f++){
			vector <bool> me(sp.markov_eqn.size(),false);
			vector <bool> nm(sp.nm_trans.size(),false);
			
			const auto &fe = sp.fix_effect[f];
			for(auto i : fe.markov_eqn_ref) me[i] = true;
			for(auto i : fe.nm_trans_ref) nm[i] = true;
			
			for(auto tr = 0u; tr < N; tr++){
				const auto &tra = sp.tra_gl[tr];
				if(tra.markov_eqn_ref != UNSET && me[tra.markov_eqn_ref]) sp.tra_fe_affect[tr][f] = true;
				if(tra.nm_trans_ref != UNSET && nm[tra.nm_trans_ref]) sp.tra_fe_affect[tr][f] = true;
			}
		}
		
		if(false){
			for(auto tr = 0u; tr < N; tr++){
				const auto &tra = sp.tra_gl[tr];
				cout << tra.name << "  IE: ";
				for(auto e = 0u; e < E; e++) cout << sp.tra_ie_affect[tr][e] << ", ";
				cout << "  FE: ";
				for(auto f = 0u; f < F; f++) cout << sp.tra_fe_affect[tr][f] << ", ";
				cout << " Affect" << endl;
			}
			emsg ("Affect");
		}
	}
}


/// References markov eqns and nm_trans ind_eff_group
void Input::ind_fix_eff_group_trans_ref()
{
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
	
		for(auto &ieg : sp.ind_eff_group){
			vector <unsigned int> ie_list;
			for(const auto &li : ieg.list) ie_list.push_back(li.index);
			
			for(auto e = 0u; e < sp.markov_eqn.size(); e++){
				for(auto ie : sp.markov_eqn[e].ind_eff_mult){
					if(find_in(ie_list,ie) != UNSET) ieg.markov_eqn_ref.push_back(e);
				}
			}
			
			for(auto e = 0u; e < sp.nm_trans.size(); e++){
				const auto &nmt = sp.nm_trans[e];
				vector <unsigned int> eq_list;
				for(auto eq : nmt.dist_param_eq_ref) eq_list.push_back(eq);
				if(nmt.bp_eq != UNSET){
					if(nmt.bp_eq == BP_FROM_OTHERS){
						for(auto e : nmt.bp_other_eq) eq_list.push_back(e);
					}
					else{
						if(nmt.all_branches){
							for(auto e : nmt.bp_all_eq) eq_list.push_back(e);
						}
						else eq_list.push_back(nmt.bp_eq);
					}
				}
				
				for(auto eq : eq_list){
					for(auto ie : model.eqn[eq].ind_eff_mult){
						if(find_in(ie_list,ie) != UNSET){
							ieg.nm_trans_ref.push_back(e);
							add_to_vec(sp.ind_effect[ie].nm_trans_ref,e);
							break;
						}
					}
					
					for(auto ie : model.eqn[eq].ind_eff_mult){
						add_to_vec(sp.ind_effect[ie].nm_trans_ref,e);
						add_to_vec(sp.ind_effect[ie].nm_trans_incomp_ref,nmt.trans_incomp_ref);
					}
					
					for(auto fe : model.eqn[eq].fix_eff_mult){
						add_to_vec(sp.fix_effect[fe].nm_trans_ref,e);
						
						add_to_vec(sp.fix_effect[fe].nm_trans_incomp_ref,nmt.trans_incomp_ref);
					}
				}
			}
		}
		
		for(auto e = 0u; e < sp.markov_eqn.size(); e++){
			for(auto ie : sp.markov_eqn[e].ind_eff_mult){
				add_to_vec(sp.ind_effect[ie].markov_eqn_ref,e);
			}
			
			for(auto fe : sp.markov_eqn[e].fix_eff_mult){
				add_to_vec(sp.fix_effect[fe].markov_eqn_ref,e);
			}
		}
	}
	
	if(false){
		for(const auto &sp : model.species){ 	
			for(auto ie =0u; ie < sp.ind_effect.size(); ie++){
				const auto &ind_eff = sp.ind_effect[ie];
				cout << ind_eff.name <<"  ind eff" << endl;
				for(auto k = 0u; k < ind_eff.nm_trans_ref.size(); k++){
					cout << ind_eff.nm_trans_ref[k] << " j" << endl;
				}					
			}
		}
		emsg_input("J");
	}
}


/// References populations in individual effect
void Input::ind_fix_eff_pop_ref()
{
	for(auto k = 0u; k < model.pop.size(); k++){
		const auto &pop = model.pop[k];
		for(auto ie : pop.ind_eff_mult){
			auto &sp = model.species[pop.sp_p];
			add_to_vec(sp.ind_effect[ie].pop_ref,k);
		}
		
		for(auto fe : pop.fix_eff_mult){
			auto &sp = model.species[pop.sp_p];
			add_to_vec(sp.fix_effect[fe].pop_ref,k);
		}
	}
}


/// Works out how markov eqns change on c->c'  
void Input::create_markov_comp_gl()
{
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];

		if(sp.type == INDIVIDUAL){
			auto N = sp.comp_gl.size();
			for(auto c = 0u; c < N; c++){
				auto &co = sp.comp_gl[c];
				
				for(auto cl = 0u; cl < sp.ncla; cl++){
					for(const auto &tgl :  co.tra_leave_group[cl].tr_list){
						const auto &tr = sp.tra_gl[tgl];
						if(tr.type == EXP_RATE || tr.type == EXP_MEAN){
							if(tr.markov_eqn_ref == UNSET) emsg_input("Should not be unset4");
							co.me_ref.push_back(tr.markov_eqn_ref);
						}
						else{
							if(tr.markov_eqn_ref != UNSET) emsg_input("Should not be set3");
						}
					}
				}
			}
		}
	}
}


// Works out how interconnected compartments are by making separate "islands"
void Input::create_island()
{
	for(auto p = 0u; p < model.nspecies; p++){
		auto &sp = model.species[p];
		
		if(sp.type == INDIVIDUAL){
			for(auto cl = 0u; cl < sp.ncla; cl++){
				auto &claa = sp.cla[cl];
				
				// Creates isolated islands for each compartment
				vector <Island> island;
				for(auto c = 0u; c < claa.ncomp; c++){
					Island isl; 
					IslandComp ico; ico.c = c; ico.nm_trans = false;
					isl.comp.push_back(ico);
					
					island.push_back(isl);
				}
					
				// Joins together islands based on transitions			
				for(auto tr = 0u; tr < claa.ntra; tr++){
					const auto &tra = claa.tra[tr];
					if(tra.variety != SOURCE_TRANS && tra.variety != SINK_TRANS){
						auto c1 = tra.i;
						auto c2 = tra.f;
						
						unsigned int isl1, j = UNSET; 
						for(isl1 = 0u; isl1 < island.size(); isl1++){
							const auto &co = island[isl1].comp;
							j = 0; while(j < co.size() && co[j].c != c1) j++;
							if(j < co.size()) break;
						}
						if(isl1 == island.size()) emsg_input("Cannot find island1");
							
						IslandTrans itr;
						itr.tr = tr;
						itr.nm_trans = false;
						island[isl1].comp[j].leave.push_back(itr);
					
						auto isl2 = UNSET;
						if(c2 != UNSET){
							for(isl2 = 0u; isl2 < island.size(); isl2++){
								const auto &co = island[isl2].comp;
								j = 0; while(j < co.size() && co[j].c != c2) j++;
								if(j < co.size()) break;
							}
							if(isl2 == island.size()) emsg_input("Cannot find island2");
						}
						
						if(isl1 != isl2){
							for(const auto &co : island[isl2].comp) island[isl1].comp.push_back(co);
							island.erase(island.begin()+isl2);
						}
					}
				}
				
				// Links from island transitions to the equations determined by markov_eqn
				const auto T = model.ntimepoint-1;
				for(auto &isl : island){
					for(auto &co : isl.comp){
						auto comp_nm_flag = false;
						for(auto &le : co.leave){
							auto tr = le.tr;
							const auto &tra = claa.tra[tr];
							
							// Sets cf which gives the final
							auto f = tra.f;
							auto j = 0u; while(j < isl.comp.size() && isl.comp[j].c != f) j++;
							if(j == isl.comp.size()) emsg_input("problem");
							le.cf = j;
							
							// Given a global compartment c  works put the me which should be used
							// Note, the contribution to in cl is ignored) hence this account for trans in other cl
							vector <unsigned int> me_ref(sp.comp_gl.size(),UNSET);
						
							auto nm_flag = false;
						
							for(auto k = 0u; k < sp.tra_gl.size(); k++){
								const auto &trg = sp.tra_gl[k];
								
								if(trg.cl == cl && trg.tr == tr){
									auto me = trg.markov_eqn_ref;
									if(me == UNSET){
										switch(trg.type){
										case EXP_RATE: case EXP_MEAN:
											emsg_input("Should not be unset5");
											break;
										default:
											nm_flag = true;
											comp_nm_flag = true;
											break;
										}
									}
									
									auto cc = trg.i - tra.i*sp.comp_mult[cl];
									if(trg.i == UNSET) cc = trg.f - tra.f*sp.comp_mult[cl];
									for(auto c = 0u; c < claa.comp.size(); c++){
										me_ref[cc+c*sp.comp_mult[cl]] = me;
									}	
								}
							}
							
							le.nm_trans = nm_flag;
							
							if(nm_flag == false){
								for(auto c = 0u; c < sp.comp_gl.size(); c++){
									if(me_ref[c] == UNSET) emsg_input("me ref prob");
								}							
								
								// Accounts for time variation  
								le.markov_eqn_ref.resize(T);
								for(auto ti = 0u; ti < T; ti++){
									auto &mer = le.markov_eqn_ref[ti];
									mer.resize(sp.comp_gl.size());
									
									for(auto c = 0u; c < sp.comp_gl.size(); c++){
										auto e = me_ref[c];
										mer[c].e = e;
										if(sp.markov_eqn[e].time_vari == false) mer[c].ti = 0;
										else mer[c].ti = ti;					
									}
								}
							}
							else{
								le.nm_trans_ref.resize(sp.comp_gl.size(),UNSET);
								
								for(auto k = 0u; k < sp.tra_gl.size(); k++){
									const auto &trg = sp.tra_gl[k];
								
									if(trg.cl == cl && trg.tr == tr){
										auto cc = trg.i - tra.i*sp.comp_mult[cl];
										if(trg.i == UNSET) cc = trg.f - tra.f*sp.comp_mult[cl];
										for(auto c = 0u; c < claa.comp.size(); c++){
											le.nm_trans_ref[cc+c*sp.comp_mult[cl]] = trg.nm_trans_ref;
										}	
									}
								}
							}
						}
						co.nm_trans	= comp_nm_flag;
					}
				}
				
				for(auto i = 0u; i < island.size(); i++){
					for(auto c = 0u; c < island[i].comp.size(); c++){
						auto &iref = claa.comp[island[i].comp[c].c].island_ref;
						iref.i = i; iref.c = c;
					}
				}
				
				if(false){		
					cout << "For classification " << claa.name << endl;
					for(const auto &isl : island){
						cout << "Island: " << endl;
						for(const auto &co : isl.comp){
							cout << claa.comp[co.c].name << ", ";
							cout << "Leave: ";
							for(const auto &le : co.leave){
								cout << claa.tra[le.tr].name << ","<< le.cf << "  ";
							}
							cout << endl;
						}
						cout << endl;
						cout << endl;
					}
					emsg_input("Islands");
				}
				
				claa.island = island;
			}
		}
	}
}


/// Works out how changes to parameters affect likelihoods
void Input::param_affect_likelihood()
{
	auto T = model.ntimepoint-1;

	// Stores spline time maps for each of the model variables
	vector < vector <bool> > spline_map;
	spline_map.resize(model.param_vec.size());
	for(auto k = 0u; k < model.param_vec.size(); k++){
		const auto &pv = model.param_vec[k];
		
		auto th = pv.th, ind = pv.index;
	
		const auto &par = model.param[th];
		
		if(par.spline_info.on == true){	
			auto nknot = par.spline_info.knot_times.size();
		
			auto s = 0u;
			while(s < model.spline.size() && !(th == model.spline[s].th && ind/nknot == model.spline[s].index)) s++;
			if(s == model.spline.size()) emsg_input("Cannot find spline");
		
			const auto &spl = model.spline[s];
			
			vector <bool> map(T,false);
			for(auto ti = 0u; ti < T; ti++){ 
				if(spl.div[ti].th1 == k || spl.div[ti].th2 == k) map[ti] = true;
			}
	
			spline_map[k] = map;
			AffectLike al; al.type = SPLINE_AFFECT; al.num = s; al.num2 = UNSET; al.map = map;
			param_vec_add_affect(model.param_vec[k].affect_like,al);		
		}		
	}

	/// The effect of parameters on Markov equations
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp =  model.species[p];
		switch(sp.type){
		case INDIVIDUAL:
			{
				for(auto i = 0u; i < sp.markov_eqn.size(); i++){
					const auto &me = sp.markov_eqn[i];
				
					const auto &eqn = model.eqn[me.eqn_ref];
				
					for(auto &pref : eqn.param_ref){
						auto th = pref.th, ind = pref.index;
						auto &par = model.param[th];
					
						if(par.variety != CONST_PARAM){
							auto k = par.get_param_vec(ind);	
							if(k != UNSET){
								if(par.spline_info.on == true){	
									auto nknot = par.spline_info.knot_times.size();
									
									for(auto t = 0u; t < nknot; t++){ // Goes down the spline
										k = par.get_param_vec(ind+t);
										AffectLike al; al.map = spline_map[k];
										al.type = DIV_VALUE_AFFECT; al.num = p; al.num2 = i;
										param_vec_add_affect(model.param_vec[k].affect_like,al);		
								
										al.type = MARKOV_LIKE_AFFECT;
										param_vec_add_affect(model.param_vec[k].affect_like,al);
									}
								}
								else{
									AffectLike al; 
									if(me.time_vari == true) al.map.resize(T,true);
									else al.map.resize(1,true);
						
									al.type = DIV_VALUE_AFFECT; al.num = p; al.num2 = i;
									param_vec_add_affect(model.param_vec[k].affect_like,al);
								
									al.type = MARKOV_LIKE_AFFECT;
									param_vec_add_affect(model.param_vec[k].affect_like,al);
								}
							}
						}
					}
				}
			}
			break;
			
		case POPULATION:
			{
				for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
					const auto &tra = sp.tra_gl[tr];
					const auto &eqn = model.eqn[tra.dist_param[0].eq_ref];
					for(auto &pref : eqn.param_ref){
						auto th = pref.th, ind = pref.index;
					
						auto &par = model.param[th];
						if(par.variety != CONST_PARAM){
							auto k = par.get_param_vec(ind);	
							if(k != UNSET){
								if(par.spline_info.on == true){
									auto nknot = par.spline_info.knot_times.size();
								
									for(auto t = 0u; t < nknot; t++){ // Goes down the spline
										k = par.get_param_vec(ind+t);
										AffectLike al; al.map = spline_map[k];
										al.type = MARKOV_POP_AFFECT; al.num = p; al.num2 = tr; 
										param_vec_add_affect(model.param_vec[k].affect_like,al);
									}
								}
								else{
									AffectLike al; al.map.resize(T,true);
									al.type = MARKOV_POP_AFFECT; al.num = p; al.num2 = tr;
									param_vec_add_affect(model.param_vec[k].affect_like,al);
								}
							}
						}
					}
				}
			}
			break;
		}
	}
	
	// Works out how the parameter affect the prior
	for(auto th = 0u; th < model.nparam_vec; th++){
		const auto &pv = model.param_vec[th];
		const auto &par = model.param[pv.th];
		if(par.variety == PRIOR_PARAM || par.variety == DIST_PARAM){
			AffectLike al;  al.num = th; al.num2 = UNSET; 
			if(par.variety == PRIOR_PARAM) al.type = PRIOR_AFFECT; else al.type = DIST_AFFECT;
			
			param_vec_add_affect(model.param_vec[th].affect_like,al);
			
			for(const auto &pref : par.get_parent(pv.index)){
				const auto th_par = model.param[pref.th].get_param_vec(pref.index);				
				if(th_par != UNSET){
					param_vec_add_affect(model.param_vec[th_par].affect_like,al);
				}
			}
		}
	}
	
	// Works out how parameter affects the spline priors
	for(auto s = 0u; s < model.spline.size(); s++){
		const auto &spl = model.spline[s];
		
		if(spl.constant == false && spl.info.smooth == true){
			AffectLike al; al.type = SPLINE_PRIOR_AFFECT; al.num = s; al.num2 = UNSET; 
			for(auto t = 0u; t < spl.param_ref.size(); t++){
				param_vec_add_affect(model.param_vec[spl.param_ref[t]].affect_like,al);
			}
		}
	}
	
	// Looks at individual fixed effects in Markov equations
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp =  model.species[p];
		if(sp.type == INDIVIDUAL){
			for(auto i = 0u; i < sp.markov_eqn.size(); i++){
				const auto &me = sp.markov_eqn[i];
				
				const auto &eqn = model.eqn[me.eqn_ref];
				for(auto j = 0u; j < eqn.fix_eff_mult.size(); j++){
					auto f = eqn.fix_eff_mult[j];
					
					const auto &fe = sp.fix_effect[f];
					
					const auto &par = model.param[fe.th];
					if(par.name != fe_char+"^"+fe.name) emsg_input("names do not match");
					if(par.N != 1) emsg_input("Should be univariate");
					
					if(par.variety != CONST_PARAM){
						auto k = par.get_param_vec(0);
						
						if(k != UNSET){				
							AffectLike al; 	
							al.type = EXP_FE_AFFECT; al.num = p; al.num2 = f;
							param_vec_add_affect(model.param_vec[k].affect_like,al);		
							
							al.type = INDFAC_INT_AFFECT; al.num = p; al.num2 = UNSET;
							param_vec_add_affect(model.param_vec[k].affect_like,al);		
							
							al.type = MARKOV_LIKE_AFFECT; al.num = p; al.num2 = i;
							if(me.time_vari == true) al.map.resize(T,true);
							else al.map.resize(1,true);
							param_vec_add_affect(model.param_vec[k].affect_like,al);		
						}
					}
				}
			}
		}
	}

	// Looks at the effect of fixed effects in populations
	for(auto po = 0u; po < model.pop.size(); po++){
		const auto &pop = model.pop[po];
		auto p = pop.sp_p;
		const auto &sp = model.species[p];
		
		if(sp.type == INDIVIDUAL){
			for(auto f : pop.fix_eff_mult){		
				const auto &fe = sp.fix_effect[f];
				const auto &par = model.param[fe.th];
				if(par.name != fe_char+"^"+fe.name) emsg_input("names do not match");
				if(par.N != 1) emsg_input("Should be univariate");
				
				if(par.variety != CONST_PARAM){
					auto k = par.get_param_vec(0);
					if(k == UNSET) emsg_input("Problem with fix_eff_mult");
					
					AffectLike al; 	
					al.type = EXP_FE_AFFECT; al.num = p; al.num2 = f;
					param_vec_add_affect(model.param_vec[k].affect_like,al);		
					
					model.add_pop_affect(po,model.param_vec[k].affect_like);
				}
			}
		}
	}
	
	// Looks at the affect of parameters which determine non-Marovian transitions
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
	
		if(sp.type == INDIVIDUAL){
			for(auto i = 0u; i < sp.nm_trans.size(); i++){
				const auto &nmt = sp.nm_trans[i];
				
				for(auto k = 0u; k < nmt.dist_param_eq_ref.size(); k++){
					add_nm_trans_affect(p,i,nmt.dist_param_eq_ref[k],NM_TRANS_AFFECT,spline_map);
				}
				
				if(nmt.bp_eq != UNSET){
					if(nmt.bp_eq == BP_FROM_OTHERS){
						for(auto e : nmt.bp_other_eq) add_nm_trans_affect(p,i,e,NM_TRANS_BP_AFFECT,spline_map);
					}
					else{
						if(nmt.all_branches){
							for(auto e : nmt.bp_all_eq){
								add_nm_trans_affect(p,i,e,NM_TRANS_BP_AFFECT,spline_map);
							}
						}
						else{
							add_nm_trans_affect(p,i,nmt.bp_eq,NM_TRANS_BP_AFFECT,spline_map);
						}
					}
				}
			}
			
			for(auto i = 0u; i < sp.nm_trans_incomp.size(); i++){
				const auto &nmti = sp.nm_trans_incomp[i];
				
				for(auto m : nmti.nmtrans_ref){
					const auto &nmt = sp.nm_trans[m];
					for(auto k = 0u; k < nmt.dist_param_eq_ref.size(); k++){
						add_nm_trans_affect(p,i,nmt.dist_param_eq_ref[k],NM_TRANS_INCOMP_AFFECT,spline_map);
					}
					
					if(nmt.bp_eq != UNSET){
						if(nmt.bp_eq == BP_FROM_OTHERS){
							for(auto e : nmt.bp_other_eq) add_nm_trans_affect(p,i,e,NM_TRANS_INCOMP_AFFECT,spline_map);
						}
						else{
							if(nmt.all_branches){
								for(auto e : nmt.bp_all_eq){
									add_nm_trans_affect(p,i,e,NM_TRANS_INCOMP_AFFECT,spline_map);
								}
							}
							else{
								add_nm_trans_affect(p,i,nmt.bp_eq,NM_TRANS_INCOMP_AFFECT,spline_map);
							}
						}
					}
				}
			}
		}
	}
	
	// Variance and correlation parameters which affect ie likelihood
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
	
		if(sp.type == INDIVIDUAL){
			for(auto g = 0u; g < sp.ind_eff_group.size(); g++){
				const auto &ieg = sp.ind_eff_group[g];
				auto N = ieg.list.size();
				
				for(auto j = 0u; j < N; j++){
					for(auto i = 0u; i < N; i++){
						auto &par = model.param[ieg.omega[j][i]];
						if(par.variety != CONST_PARAM){
							auto k = par.get_param_vec(0);
							if(k == UNSET) emsg_input("Should not be unset6");
							
							AffectLike al; 	
							al.type = OMEGA_AFFECT; al.num = p; al.num2 = g;
							param_vec_add_affect(model.param_vec[k].affect_like,al);	
							
							if(j == i){  // If on the di
								auto ie = ieg.list[i].index;
								
								model.add_ie_affect(p,ie,model.param_vec[k].affect_like);
							}
							
							{
								AffectLike al; 	
								al.type = LIKE_IE_AFFECT; al.num = p; al.num2 = g;
								param_vec_add_affect(model.param_vec[k].affect_like,al);	
							}
						}						
					}
				}
			}
		}
	}
	
	// Effect of parameters on obs_eqn
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
		auto E = sp.obs_eqn.size();
		vector < vector <bool> > map;
		map.resize(model.param_vec.size());
		for(auto k = 0u; k < model.param_vec.size(); k++) map[k].resize(E,false);
		
		for(auto i = 0u; i < E; i++) add_to_map(sp.obs_eqn[i],i,map,spline_map);
		
		add_map_list(map,p,OBS_EQN_AFFECT);
	}
	
	// Effect of parameters on individual observation probability
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
		if(sp.type == INDIVIDUAL){
			auto N = sp.individual.size();
			
			vector < vector <bool> > map;
			map.resize(model.param_vec.size());
			for(auto k = 0u; k < model.param_vec.size(); k++) map[k].resize(N,false);
			// Here the map has dimensions of param and individual
			for(auto i = 0u; i < N; i++){
				const auto &ind = sp.individual[i];
				
				auto ref = ind.enter_ref;
				if(ref != UNSET){    // Probability when adding individuals
					const auto &enter = sp.enter[ref];
					if(enter.c_set == UNSET){
						for(auto cl = 0u; cl < sp.ncla; cl++){
							const auto &ent_cl = enter.cla[cl];
							if(ent_cl.c_set == UNSET){
								const auto &claa = sp.cla[cl];
								for(auto c = 0u; c < claa.ncomp; c++){
									auto eq = ent_cl.eqn[c].eq_ref;
									add_to_map(eq,i,map,spline_map);
								}
							}
						}
					}
				}
				
				for(auto m = 0u; m < ind.obs.size(); m++){
					const auto &ob = ind.obs[m];
				
					switch(ob.type){
					case OBS_TRANS_EV: case OBS_SOURCE_EV: case OBS_SINK_EV:
						if(ob.time_vari){
							const auto &obs_eqn_ref = ob.obs_eqn_ref;
							for(auto tr = 0u; tr < obs_eqn_ref.size(); tr++){
								auto eq = sp.obs_eqn[obs_eqn_ref[tr]];
								add_to_map(eq,i,map,spline_map);
							}
						}
						else{
							const auto &obs_eqn_ref = sp.obs_trans[ob.ref].obs_eqn_ref;
							for(auto tr = 0u; tr < obs_eqn_ref.size(); tr++){
								auto eq = sp.obs_eqn[obs_eqn_ref[tr]];
								add_to_map(eq,i,map,spline_map);
							}
						}
						break;
					
					case OBS_COMP_EV:
						for(auto c = 0u; c < ob.c_obs_prob_eqn.size(); c++){
							auto eq = ob.c_obs_prob_eqn[c].eq_ref;
							if(eq == UNSET) emsg_input("Equation is not set");
							add_to_map(eq,i,map,spline_map);
						}
						break;
						
					case OBS_TEST_EV:
						{
							auto eq = ob.Se_eqn.eq_ref;
							add_to_map(eq,i,map,spline_map);
							eq = ob.Sp_eqn.eq_ref;
							add_to_map(eq,i,map,spline_map);
						}
						break;
					}
				}
			}
			
			if(false){			
				for(auto th = 0u; th < model.param_vec.size(); th++){
					cout <<  model.param_vec[th].name <<": ";
					for(auto i = 0u; i < N; i++){
						cout << map[th][i] <<",";
					}
					cout << endl;
				}
				emsg_input("done");
			}
			
			add_map_list(map,p,LIKE_OBS_IND_AFFECT);
		}
	}
	
	// Effect of parameters on obs_trans_eqn_value
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
		for(auto e = 0u; e < sp.obs_trans_eqn.size(); e++){
			add_obs_trans_eqn(p,e,spline_map);
		}
	}
	
	// Effect of parameters on pop_data observation probability
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
		
		auto N = sp.pop_data.size();
			
		vector < vector <bool> > map;
		map.resize(model.param_vec.size());
		for(auto k = 0u; k < model.param_vec.size(); k++) map[k].resize(N,false);
	
		for(auto i = 0u; i < N; i++){
			auto &popd = sp.pop_data[i];
			if(popd.time_vari == true) omr_add_to_map(popd.comp_obs_mod_ref,i,map,spline_map,sp);
			else omr_add_to_map(sp.pop_filter[popd.ref].comp_obs_mod_ref,i,map,spline_map,sp);
		}
		add_map_list(map,p,LIKE_OBS_POP_AFFECT);
	}
	
	// Effect of parameters on pop_trans_data observation probability
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
		
		auto N = sp.pop_trans_data.size();
			
		vector < vector <bool> > map;
		map.resize(model.param_vec.size());
		for(auto k = 0u; k < model.param_vec.size(); k++) map[k].resize(N,false);
	
		for(auto i = 0u; i < N; i++){
			auto &poptd = sp.pop_trans_data[i];
			
			if(poptd.time_vari == true){
				omr_add_to_map(poptd.trans_obs_mod_ref,i,map,spline_map,sp);
			}
			else{
				omr_add_to_map(sp.pop_trans_filter[poptd.ref].trans_obs_mod_ref,i,map,spline_map,sp);
			}	
		}
		add_map_list(map,p,LIKE_OBS_POP_TRANS_AFFECT);
	}
	
	// Effect of parameters on genetic data 
	const auto &gen_data = model.genetic_data;
	if(gen_data.on){
		{
			const auto &eqn = model.eqn[gen_data.mut_rate.eq_ref];
			for(auto &pref : eqn.param_ref){
				const auto &par = model.param[pref.th];
				if(par.variety != CONST_PARAM){
					auto k = par.get_param_vec(pref.index);
					if(k == UNSET) emsg_input("Parameter prob");
					
					{
						AffectLike al; 	
						al.type = GENETIC_VALUE_AFFECT; al.num = UNSET; al.num2 = UNSET;
						param_vec_add_affect(model.param_vec[k].affect_like,al);	
					}
					
					{
						AffectLike al; 	
						al.type = LIKE_GENETIC_PROCESS_AFFECT; al.num = UNSET; al.num2 = UNSET;
						param_vec_add_affect(model.param_vec[k].affect_like,al);	
					}
				}
			}
		}
		
		{
			const auto &eqn = model.eqn[gen_data.seq_var.eq_ref];
			for(auto &pref : eqn.param_ref){
				const auto &par = model.param[pref.th];
				if(par.variety != CONST_PARAM){
					auto k = par.get_param_vec(pref.index);
					if(k == UNSET) emsg_input("Parameter prob");
					
					{
						AffectLike al; 	
						al.type = GENETIC_VALUE_AFFECT; al.num = UNSET; al.num2 = UNSET;
						param_vec_add_affect(model.param_vec[k].affect_like,al);	
					}
					
					{
						AffectLike al; 	
						al.type = LIKE_GENETIC_PROCESS_AFFECT; al.num = UNSET; al.num2 = UNSET;
						param_vec_add_affect(model.param_vec[k].affect_like,al);	
					}
				}
			}
		}
	}
}


/// Adds a vector of observation probabilites to map
void Input::omr_add_to_map(vector <unsigned int> &omr, unsigned int i, vector < vector <bool> > &map, const vector < vector <bool> > &spline_map, const Species &sp) const
{
	for(auto k = 0u; k < omr.size(); k++){
		auto eq = sp.obs_eqn[omr[k]];
		add_to_map(eq,i,map,spline_map);
	}
}


/// Add affect terms for NM_TRANS_AFFECT, NM_TRANS_BP_AFFECT and NM_TRANS_INCOMP
void Input::add_nm_trans_affect(unsigned int p, unsigned int i, unsigned int eq, AffectType type, const vector < vector <bool> > &spline_map)
{
	auto T = model.ntimepoint-1;
	auto &sp = model.species[p];

	auto &eqn = model.eqn[eq];
	for(const auto &pref : eqn.param_ref){
		auto th = pref.th, ind = pref.index;

		auto &par = model.param[th];
		if(par.variety != CONST_PARAM){
			auto k = par.get_param_vec(ind);	
			if(k != UNSET){
				if(par.spline_info.on == true){	
					auto nknot = par.spline_info.knot_times.size();
					
					for(auto t = 0u; t < nknot; t++){ // Goes down the spline
						k = par.get_param_vec(ind+t);
						AffectLike al; al.type = type; al.num = p; al.num2 = i;
						al.map = spline_map[k];
						param_vec_add_affect(model.param_vec[k].affect_like,al);		
					}
				}
				else{
					AffectLike al; al.type = type; al.num = p; al.num2 = i;
					al.map.resize(T,true);
					param_vec_add_affect(model.param_vec[k].affect_like,al);		
				}
			}
		}
	}
	
	// Looks at the effects of fixed effects on transitions
	for(auto j = 0u; j < eqn.fix_eff_mult.size(); j++){
		auto f = eqn.fix_eff_mult[j];
	
		const auto &fe = sp.fix_effect[f];
	
		const auto &par = model.param[fe.th];
		
		if(par.name != fe_char+"^"+fe.name) emsg_input("names do not match");
		if(par.N != 1) emsg_input("Should be univariate");
		
		if(par.variety != CONST_PARAM){
			auto k = par.get_param_vec(0);
		
			if(k != UNSET){				
				AffectLike al; 	
				al.type = EXP_FE_AFFECT; al.num = p; al.num2 = f;
				param_vec_add_affect(model.param_vec[k].affect_like,al);		
			
				al.type = type; al.num = p; al.num2 = i;
				al.map.resize(T,true);
				param_vec_add_affect(model.param_vec[k].affect_like,al);		
			}
		}
	}
}

					
/// Adds the effect parameters have on obs_trans value
void Input::add_obs_trans_eqn(unsigned int p, unsigned int e, const vector < vector <bool> > &spline_map)
{
	auto &sp = model.species[p];
	auto eq = sp.obs_trans_eqn[e];
	
	const auto &eqn = model.eqn[eq];
	
	vector <bool> map_all(sp.T,true);

	for(auto &pref : eqn.param_ref){
		const auto &par = model.param[pref.th];
		if(par.variety != CONST_PARAM){	
			auto k = par.get_param_vec(pref.index);
			if(k == UNSET) emsg_input("Parameter prob");
			
			if(spline_map[k].size() > 0){                // There is a spline in the equation	
				if(!par.spline_info.on) emsg_input("Should be spline");
				auto nknot = par.spline_info.knot_times.size();
				
				for(auto j = 0u; j < nknot; j++){
					if(k+j >= spline_map.size()) emsg_input("Prob0");
					
					{
						AffectLike al; 
						al.type = LIKE_UNOBS_TRANS_AFFECT; al.num = p; al.num2 = e; al.map = spline_map[k+j];
						param_vec_add_affect(model.param_vec[k+j].affect_like,al);
					}
				}
			}
			else{
				{
					AffectLike al; al.type =  LIKE_UNOBS_TRANS_AFFECT; al.num = p; al.num2 = e; al.map = map_all;
					param_vec_add_affect(model.param_vec[k].affect_like,al);
				}
			}
		}		
	}
}
		
		
/// Based on an equation and index adds all parameters affected into map
void Input::add_to_map(unsigned int eq, unsigned int i, vector < vector <bool> > &map, const vector < vector <bool> > &spline_map) const
{
	const auto &eqn = model.eqn[eq];
	for(auto &pref : eqn.param_ref){
		const auto &par = model.param[pref.th];
		if(par.variety != CONST_PARAM){
			auto k = par.get_param_vec(pref.index);
			if(k == UNSET) emsg_input("Parameter prob");
			
			if(spline_map[k].size() > 0){                // There is a spline in the equation
				auto ti = eqn.ti_fix;
				if(ti == UNSET) emsg_input("Time should be fixed");
				
				if(!par.spline_info.on) emsg_input("SHould be spline");
				auto nknot = par.spline_info.knot_times.size();
				
				for(auto j = 0u; j < nknot; j++){
					if(k+j >= spline_map.size()) emsg_input("Prob0");
					if(ti >= spline_map[k+j].size()){
						emsg_input("Prob-1");
					}
					if(spline_map[k+j][ti] == true){
						if(k+j >= map.size()) emsg_input("Prob1");
						if(i >= map[k+j].size()) emsg_input("Prob2");
						map[k+j][i] = true;
					}
				}
			}
			else{
				if(k >= map.size()) emsg_input("Prob3");
				if(i >= map[k].size()) emsg_input("Prob4");
				map[k][i] = true;
			}
		}
	}
}


/// Based on a map of changes adds in how parameters are affected
void Input::add_map_list(const vector < vector <bool> > &map, unsigned int p, AffectType type) const
{
	if(false){
		cout << "add map list" << endl;
		for(auto k = 0u; k < map.size(); k++){
			for(auto i = 0u; i < map[k].size(); i++) cout << map[k][i] << ",";
			cout << "kk" << endl;
		}
	}
	
	for(auto k = 0u; k < model.param_vec.size(); k++){
		if(find_in(map[k],true) != UNSET){		
			AffectLike al; al.type = type; al.num = p; al.num2 = UNSET; al.map = map[k];
			param_vec_add_affect(model.param_vec[k].affect_like,al);
		}
	}
}


/// Works out a convertion from tr to trg
void Input::create_trg_from_tr()
{
	for(auto p = 0u; p < model.nspecies; p++){
		auto &sp = model.species[p];
		
		sp.trg_from_tr.resize(sp.comp_gl.size());
		
		for(auto c = 0u; c < sp.comp_gl.size(); c++){
			sp.trg_from_tr[c].resize(sp.ncla);
			for(auto cl = 0u; cl < sp.ncla; cl++){
				const auto &claa = sp.cla[cl];
				sp.trg_from_tr[c][cl].resize(claa.tra.size(),UNSET);
			}
		}
		
		for(auto trg = 0u; trg < sp.tra_gl.size(); trg++){
			const auto &tra = sp.tra_gl[trg];
			auto tr = tra.tr;
			auto cl = tra.cl;
			
			auto c = tra.i; if(c == UNSET) c = tra.f;
			
			auto cgl = sp.comp_gl[c];
			auto ci = cgl.cla_comp[cl];
			auto mult = sp.comp_mult[cl];
			c -= ci*mult;
			
			for(auto ci = 0u; ci < sp.cla[cl].comp.size(); ci++){
				sp.trg_from_tr[c+ci*mult][cl][tr] = trg;
			}
		}
	}
}	


/// Divides source rate (e.g. when both sexes then rate needs to be half) 
void Input::source_rate_divide() const
{
	for(auto p = 0u; p < model.nspecies; p++){
		auto &sp = model.species[p];
		
		for(auto cl = 0u; cl < sp.ncla; cl++){
			const auto &claa = sp.cla[cl];
			for(auto i = 0u; i < claa.ntra; i++){
				const auto &tr = claa.tra[i];
				if(tr.i == UNSET){
					vector <SourceTransRate> list;
					for(auto j = 0u; j < sp.tra_gl.size(); j++){
						const auto &trg  = sp.tra_gl[j];
						if(trg.cl == cl && trg.tr == i){
							auto rate = trg.dist_param[0].te;
							
							auto k = 0u; while(k < list.size() && list[k].rate != rate) k++;
							if(k == list.size()){
								SourceTransRate str; str.rate = rate;
								list.push_back(str);
							}
							list[k].tr_gl.push_back(j);
						}
					}
					for(auto val:list){
						if(val.tr_gl.size() > 1){
							auto fac = 1.0/val.tr_gl.size();
							for(auto j : val.tr_gl){
								auto &trg  = sp.tra_gl[j];
								switch(trg.type){
								case EXP_RATE: trg.dist_param[0].te = to_string(fac)+"*("+val.rate+")"; break;
								case EXP_MEAN: trg.dist_param[0].te = "("+val.rate+")/"+to_string(fac); break;
								default: emsg("Source not exponential"); break;
								}
							}
						}
					}
				}
			}
		}
	}		
}


/// Sets up the property is_one in obs_trans
void Input::setup_obs_trans_is_one()
{
	// Sets is_one in obs_trans
	for(auto &sp : model.species){
		for(auto &ot : sp.obs_trans){
			ot.is_one.resize(sp.tra_gl.size(),false);
			
			vector <unsigned> tr_list;
			for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
				const auto &eq = model.eqn[ot.tra_prob_eqn[tr].eq_ref];
				if(eq.is_one()) ot.is_one[tr] = true;
				if(!eq.is_zero()){
					const auto &tra = sp.tra_gl[tr];
					add_to_vec(tr_list,tra.tr);
				}
			}
			if(tr_list.size() <= 1) ot.single_trans = true;
			else ot.single_trans = false;
		}
	}
}


/// Sets up obs_trans (this references obs_trans objects from [ti][tr_gl][#])
void Input::setup_obs_trans()
{
	for(auto &sp : model.species){
		sp.last_obs_trans_ti = UNSET;
		
		sp.obs_trans_eqn_ref.resize(sp.tra_gl.size());
		for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
			sp.obs_trans_eqn_ref[tr].resize(sp.T);
		}
	
		for(auto i = 0u; i < sp.obs_trans.size(); i++){
			auto &ot = sp.obs_trans[i];
	
			for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
				auto eq = ot.tra_prob_eqn[tr].eq_ref;
				
				if(!model.eqn[eq].is_zero()){
					auto index = add_to_vec(sp.obs_trans_eqn,eq);
					for(auto ti = ot.ti_min; ti < ot.ti_max; ti++){
						sp.obs_trans_eqn_ref[tr][ti].push_back(index);
					}
					
					if(sp.last_obs_trans_ti == UNSET) sp.last_obs_trans_ti = ot.ti_max;
					else{
						if(ot.ti_max > sp.last_obs_trans_ti){
							sp.last_obs_trans_ti = ot.ti_max;
						}
					}
				}
			}
		}
	}
}


/// Create the equation for the unset branch probability
void Input::bp_create_unset()
{
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
		for(auto cl = 0u; cl < sp.ncla; cl++){
			auto &claa = sp.cla[cl];
			for(auto &tra : claa.tra){
				tra.all_branches = false;
			}
			
			for(auto c = 0u; c < claa.ncomp; c++){
				vector <unsigned int> set, unset;
				for(auto tr = 0u; tr < claa.ntra; tr++){
					const auto &tra = claa.tra[tr];
					if(tra.i == c && tra.branch == true){
						if(tra.bp_set == true) set.push_back(tr);
						else unset.push_back(tr);
					}
				}
				
				if(set.size() + unset.size() > 1){	
					switch(unset.size()){
					case 0:
						for(auto tr : set){
							auto &tra = claa.tra[tr];
							tra.all_branches = true; 
							if(tra.i != UNSET) claa.comp[tra.i].all_branches = true; 
						}
						break;
						
					case 1:
						{
							auto &tra = claa.tra[unset[0]];
							
							string st = "1 - (";
							auto flag = false;
							for(auto tr : set){
								if(flag == true) st += "+";
								flag = true;
								st += "("+claa.tra[tr].bp.te_raw+")";
							}
							st += ")";
							
							if(!includes(st,"<") && !includes(st,"[")){ 
								tra.bp = he(add_equation_info(st,BP,p,cl));
								tra.bp_set = BP_SET;
							}
							else{			
								tra.bp_set = BP_DERIVED;
							}
						}
						break;
					
					default:
						emsg_input("Branching probabilities should be set");
						break;
					}
				}
			}
		}
	}
}


/// If doesn't contain a population then EXP_RATE -> EXP_RATE_NM
void Input::exp_nm_convert()
{
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
		if(sp.type == INDIVIDUAL){  
			for(auto cl = 0u; cl < sp.ncla; cl++){
				auto &claa = sp.cla[cl];
				for(auto tr = 0u; tr < claa.ntra; tr++){
					auto &tra = claa.tra[tr];
					
					if(tra.type == EXP_RATE && tra.i != UNSET){
						const auto &co = claa.comp[tra.i];
				
						if(co.branching == false){
							// Looks to see if transition contains any time variation
							
							vector <unsigned int> list;
							
							auto time_vari = false;
							for(auto k = 0u; k < sp.tra_gl.size(); k++){
								auto &trg = sp.tra_gl[k];
								
								if(trg.cl == cl && trg.tr == tr){
									const auto &eq = model.eqn[trg.dist_param[0].eq_ref];
									if(eq.time_vari) time_vari = true;
									list.push_back(k);
								}
							}
								
							if(time_vari == false){
								tra.type = EXP_RATE_NM;
								tra.dist_param[0].type = TRANS_NM_RATE;
								for(auto k : list){
									auto &trg = sp.tra_gl[k];
									trg.type = EXP_RATE_NM;
									trg.dist_param[0].type = TRANS_NM_RATE;
									auto c = trg.i;
									if(c != UNSET){
										sp.comp_gl[c].tra_leave_group[trg.cl].markov = false;
										sp.comp_gl[c].tr_leave_markov.clear();
									}
								}
							}
						}
					}
				}
			}
		}
	}	
}


/// If all branches are rates then removes branching
void Input::markov_bp_convert()
{
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
		for(auto cl = 0u; cl < sp.ncla; cl++){
			auto &claa = sp.cla[cl];
			for(auto c = 0u; c < claa.ncomp; c++){
				vector <unsigned int> list;
				auto flag = false;
				for(auto tr = 0u; tr < claa.ntra; tr++){
					const auto &tra = claa.tra[tr];
					if(tra.i == c && tra.branch == true){
						if(tra.type != EXP_RATE && tra.type != EXP_MEAN) flag = true;	
						if(tra.bp_set != BP_SET) flag = true;
						else{
							auto st = tra.bp.te;
							if(includes(st,"<") || includes(st,"[")) flag = true;
						}
						
						list.push_back(tr);
					}
				}
				
				if(list.size() > 0){
					if(flag == false){
						string div = "";
						for(auto tr : list){
							if(div != "") div += "+";
							div += "("+claa.tra[tr].bp.te+")";
						}					
					
						for(auto tr : list){	
							auto &tra = claa.tra[tr];
							tra.branch = false;
							tra.bp_set = BP_UNSET;
							auto bp_te = tra.bp.te;
							if(tra.all_branches) bp_te = "("+tra.bp.te+")/("+div+")";
							
							if(tra.type == EXP_RATE){
								tra.dist_param[0].te = "("+bp_te+")*("+tra.dist_param[0].te+")";
							}
							else{
								tra.dist_param[0].te = "("+tra.dist_param[0].te+")/("+bp_te+")";
							}
							tra.bp.te = "";
						}
						auto &co = claa.comp[c];
						co.markov_branch = false;
						co.branch = false; 
					}
					else{
						for(auto tr : list){
							auto &tra = claa.tra[tr];
							if(tra.type == EXP_RATE){
								tra.type = EXP_RATE_NM;
								tra.dist_param[0].type = TRANS_NM_RATE;
							}
							
							if(tra.type == EXP_MEAN){
								tra.type = EXP_MEAN_NM;
								tra.dist_param[0].type = TRANS_NM_MEAN;
							}
						}					
					}
				}
			}
		}
	}
}


/// Sets up information about infection down transitions (for trans-tree)
void Input::setup_trans_infection()
{
	for(auto &sp : model.species){
		for(auto &tra : sp.tra_gl){
			auto &inf = tra.infection;
			if(inf.type == TRANS_INFECTION){
				if(tra.type != EXP_RATE && tra.type != EXP_MEAN){
					alert_line("This infection transition must be Markovian.",tra.line_num);
				}
				const auto &eqn = model.eqn[tra.dist_param[0].eq_ref];
				if(eqn.linearise.on != true){
					alert_line("For transmission trees the rate must be linearly expressed in terms of populations.",tra.line_num);
				}
			}
		}
	}
}


/// Sets if there is a transition tree somewhere in the model
void Input::set_trans_tree()
{
	model.trans_tree = false;
	for(const auto &sp : model.species){
		if(sp.trans_tree)	model.trans_tree = true;
	}
}


/// Used to order genetic data
bool ObsGeneticData_ord (const ObsGeneticData &in1, const ObsGeneticData &in2)
{ return (in1.t < in2.t); };  


/// Used to order individual genetic observations 
bool ObsGenRef_ord (const ObsGenRef &ogr1, const ObsGenRef &ogr2)                      
{ return (ogr1.t < ogr2.t); };  


/// Adds genetic data to the model (this is potentially multi-species)
void Input::add_genetic_data()
{
	auto &gd = model.genetic_data;
				
	gd.on = false;
	
	// Makes a list of individuals across all species
	
	for(auto &sp : model.species){	
		for(const auto &so : sp.source){ 
			if(so.cname == GENETIC_DATA){
				if(model.genetic_data.on == true){
					alert_import("Cannot have more than one genetic data."); return;
				}
			
				gd.on = true;
			
				const auto &tab = so.table;
				
				auto nobs = tab.nrow;
				
				gd.mut_rate = he(add_equation_info(so.mut_rate_str,MUT_RATE));
				
				gd.seq_var = he(add_equation_info(so.seq_var_str,SEQ_VAR));
				
				switch(so.gen_data_type){
				case SNP_DATA:
					{
						auto nSNP = tab.ncol-2;
					
						gd.nSNP = nSNP;
						
						for(auto r = 0u; r < tab.nrow; r++){
							auto name = tab.ele[r][0];
							
							auto ai = model.find_all_ind(name);
							
							ObsGeneticData ogd;
							ogd.name = "";
							ogd.p = ai.p;
							ogd.i = ai.i;
						
							auto t_str = tab.ele[r][1];
							ogd.t = number(t_str);
							if(ogd.t == UNSET){ alert_import("'"+t_str+"' is not a number"); return;}
						
							for(auto j = 0u; j < nSNP; j++){
								auto ch = tab.ele[r][j+2];
								if(ch == "A") ogd.snp.push_back(A_CH);
								else{
									if(ch == "C") ogd.snp.push_back(C_CH);
									else{	
										if(ch == "T") ogd.snp.push_back(T_CH);
										else{
											if(ch == "G") ogd.snp.push_back(G_CH);
											else{
												alert_import("SNP character '"+ch+"' not recognised");
											}
										}
									}
								}
							}
							
							gd.obs.push_back(ogd);
						}
						
						sort(gd.obs.begin(),gd.obs.end(),ObsGeneticData_ord);    
							
						gd.gen_dif.resize(nobs);
						for(auto j = 0u; j < nobs; j++){
							gd.gen_dif[j].resize(nobs);
							for(auto i = 0u; i < nobs; i++){
								auto num = 0u;
								for(auto c = 0u; c < nSNP; c++){
									if(gd.obs[j].snp[c] != gd.obs[i].snp[c]) num++;
								}
								gd.gen_dif[j][i] = num;
							}
						}
					}
					break;
				
				case MATRIX_DATA:
					{
						gd.nSNP = UNSET;
						
						for(auto r = 0u; r < tab.nrow; r++){
							auto name = tab.ele[r][0];
							
							auto ai = model.find_all_ind(name);
							
							ObsGeneticData ogd;
							ogd.name = tab.ele[r][2];
							ogd.p = ai.p;
							ogd.i = ai.i;
							ogd.r = r;
						
							auto t_str = tab.ele[r][1];
							ogd.t = number(t_str);
							if(ogd.t == UNSET){ alert_import("'"+t_str+"' is not a number"); return;}
							
							gd.obs.push_back(ogd);
						}
						
						sort(gd.obs.begin(),gd.obs.end(),ObsGeneticData_ord);    
						 
						gd.gen_dif.resize(nobs);
						for(auto j = 0u; j < nobs; j++){
							gd.gen_dif[j].resize(nobs);
							for(auto i = 0u; i < nobs; i++){
								auto rj = gd.obs[j].r;
								auto ri = gd.obs[i].r;
								gd.gen_dif[j][i] = number(tab.ele[rj][ri+3]);
							}
						}
					}
					break;
				}					

				if(false){
					for(auto j = 0u; j < 10; j++){
						for(auto i = 0u; i < 10; i++){
							cout << gd.gen_dif[j][i] << ",";
						}
						cout << endl;
					}
					emsg_input("Genetic difference");
				}
			}
		}
	}

	{                                              // Sets up ind_gen_obs
		gd.ind_gen_obs.resize(model.species.size());
		for(auto p = 0u; p < model.species.size(); p++){
			const auto &sp = model.species[p];
			if(sp.trans_tree){	
				auto nind = sp.individual.size();
				gd.ind_gen_obs[p].resize(nind);
			}
		}
		
		for(auto m = 0u; m < gd.obs.size(); m++){
			const auto &ob = gd.obs[m];
			ObsGenRef ogr; ogr.m = m; ogr.t = ob.t;
			gd.ind_gen_obs[ob.p][ob.i].push_back(ogr);
		}		
		
		for(auto p = 0u; p < model.species.size(); p++){
			for(auto &ind : gd.ind_gen_obs[p]){
				sort(ind.begin(),ind.end(),ObsGenRef_ord);    
			}
		}			
	}
	
	gd.obs_scale = UNSET;
	if(gd.on){
		auto nobs = gd.gen_dif.size();
		
		auto mean = 0.0;
		for(auto j = 0u; j < nobs; j++){
			for(auto i = 0u; i < nobs; i++) mean += gd.gen_dif[j][i];
		}
		mean /= (nobs*nobs);
		gd.obs_scale = 10.0/(mean*mean);
	}
}


/// Determines if source transition exists in model
void Input::set_contains_source_sink()
{
	for(auto &sp : model.species){
		sp.contains_source = false;
		sp.contains_sink = false;
		for(const auto &claa : sp.cla){
			for(const auto &tr : claa.tra){
				if(tr.variety == SOURCE_TRANS) sp.contains_source = true;
				if(tr.variety == SINK_TRANS) sp.contains_sink = true;
			}
		}
	}
}


/// Used to order trigger events
bool TrigEvRef_ord(const TrigEventRef &st1, const TrigEventRef &st2)                      
{ return (st1.t < st2.t); };  


/// Creates timeline of events used for individual simulation 
void Input::create_trig_ev_ref()
{
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
		
		if(sp.type == INDIVIDUAL){
			for(auto i = 0u; i < sp.nindividual_in; i++){
				auto &trig_event = sp.individual[i].trig_ev_ref;

				// Adds ovserved 
				const auto &obs = sp.individual[i].obs;
				for(auto m = 0u; m < obs.size(); m++){
					switch(obs[m].type){
					case OBS_TRANS_EV:case OBS_SINK_EV:
						{
							if(!trig_ev_exist(TRIG_OBS_TRANS_EV,obs[m].t,trig_event)){
								TrigEventRef ste; 
								ste.type = TRIG_OBS_TRANS_EV; 
								ste.ref = m; ste.t = obs[m].t;
								trig_event.push_back(ste);
							}
						}
						break;
		
					default:
						break;
					}
				}
		
				const auto &ev = sp.individual[i].ev;
				for(auto m = 0u; m < ev.size(); m++){
					const auto &eve = ev[m];
					switch(eve.type){
					case LEAVE_EV:
						{
							if(!trig_ev_exist(TRIG_LEAVE_EV,eve.t,trig_event)){
								TrigEventRef ste; 
								ste.type = TRIG_LEAVE_EV; ste.ref = UNSET; ste.t = eve.t;
								trig_event.push_back(ste);
							}
						}
						break;
						
					case MOVE_EV:
						{
							TrigEventRef ste; 
							ste.type = TRIG_MOVE_EV; ste.ref = m; ste.t = eve.t;
							trig_event.push_back(ste);
						}
						break;
					
					default:
						break;
					}
				}	

				sort(trig_event.begin(),trig_event.end(),TrigEvRef_ord);		
			}
		}
	}		
}


/// Determines if a trigger event already exists at a given time
bool Input::trig_ev_exist(TrigEventType type, double t, const vector <TrigEventRef> &trig_event) const 
{
	auto len = trig_event.size();
	if(len > 0){
		auto &te = trig_event[len-1];
		if(te.t == t && te.type == type) return true;
	}
	
	return false;
}


/// Creates timeline of events used for individual sampling within a classification
void Input::create_cl_trig_ev_ref()
{
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
		
		if(sp.type == INDIVIDUAL){
			for(auto i = 0u; i < sp.nindividual_in; i++){
				auto &cl_trig_event = sp.individual[i].cl_trig_ev_ref;
				cl_trig_event.resize(sp.ncla);
				
				auto &fixed_trans_ev = sp.individual[i].fixed_trans_ev;
				fixed_trans_ev.resize(sp.ncla);

				// Adds observed 
				const auto &obs = sp.individual[i].obs;
				for(auto m = 0u; m < obs.size(); m++){
					const auto &ob = obs[m];
					switch(ob.type){
					case OBS_TRANS_EV:
						{
							if(!trig_ev_exist(TRIG_OBS_TRANS_EV,obs[m].t,cl_trig_event[ob.cl])){
								TrigEventRef ste; 
								ste.type = TRIG_OBS_TRANS_EV; ste.ref = m; ste.t = obs[m].t;
								cl_trig_event[ob.cl].push_back(ste);
								fixed_trans_ev[ob.cl].push_back(ste);
							}
						}
						break;
				
					default:
						break;
					}
				}
		
				const auto &ev = sp.individual[i].ev;
				for(auto m = 0u; m < ev.size(); m++){
					const auto &eve = ev[m];
					switch(eve.type){			
					case MOVE_EV:
						{
							TrigEventRef ste; 
							ste.type = TRIG_MOVE_EV; ste.ref = m; ste.t = eve.t;
							cl_trig_event[eve.cl].push_back(ste);
							fixed_trans_ev[eve.cl].push_back(ste);
						}
						break;
					
					default:
						break;
					}
				}	

				for(auto cl = 0u; cl < sp.ncla; cl++){
					sort(cl_trig_event[cl].begin(),cl_trig_event[cl].end(),TrigEvRef_ord);	
					sort(fixed_trans_ev[cl].begin(),fixed_trans_ev[cl].end(),TrigEvRef_ord);	
				}					
			}
		}
	}		
}


/// Sets any joint event/parameter proposals
void Input::set_joint_param_event()
{
	vector <ParEventJointInfo> par_joint_info(model.param_vec.size());
	
	for(auto p = 0u; p < model.species.size(); p++){
		const auto &sp = model.species[p];
		for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
			const auto &tra = sp.tra_gl[tr];
			if(tra.type != EXP_RATE && tra.type != EXP_MEAN){
				const auto &eq = model.eqn[tra.dist_param[0].eq_ref];
				
				for(const auto &pref : eq.param_ref){
					const auto &par = model.param[pref.th];
					
					if(par.variety != CONST_PARAM){
						auto th = par.get_param_vec(pref.index);  
						
						if(eq.param_linear(th)){	
							auto &pji = par_joint_info[th];
							
							if(pji.p == UNSET) pji.p = p;
							else{
								if(pji.p != p) pji.possible = false;
							}
							pji.tr_list.push_back(tr);
						}
					}
				}
			}
		}
	}
	
	for(auto th = 0u; th < model.param_vec.size(); th++){
		const auto &pji = par_joint_info[th];
		if(pji.possible && pji.tr_list.size() > 0){
			auto &sp = model.species[pji.p];
			if(sp.type == INDIVIDUAL){
				ParEventJointProp jp; 
				jp.th = th;
				jp.tr_list = pji.tr_list;
				jp.dir = FORWARD;
				sp.par_event_joint.push_back(jp); 

				auto after_fl = false;
				for(auto tr : pji.tr_list){
					const auto &tra = sp.tra_gl[tr];
					if(tra.f != UNSET){
						const auto &co = sp.comp_gl[tra.f];
						if(co.tr_leave.size() > 0){ after_fl = true; break;}
					}
				}
				
				if(after_fl == true){ 
					jp.dir = FORWARD_SQ;
					sp.par_event_joint.push_back(jp); 
				}
				
				if(false){ // Backward proposals were not found to be very effective
					auto before_fl = false;
					for(auto tr : pji.tr_list){
						const auto &tra = sp.tra_gl[tr];
						if(tra.i != UNSET){
							const auto &co = sp.comp_gl[tra.i];
							if(co.tr_enter.size() > 0){ before_fl = true; break;}
						}
					}
					
					if(before_fl == true){ 
						jp.dir = BACKWARD_SQ;
						sp.par_event_joint.push_back(jp); 
					}
				}
			}
		}
	}
	
	if(false){
		cout << "Joint proposals" << endl;
		for(const auto &sp : model.species){ 
			for(const auto &pej : sp.par_event_joint){
				cout << model.param_vec[pej.th].name << ":" << endl; 
				for(const auto tr : pej.tr_list){
					cout << sp.tra_gl[tr].name << ",";
				}
				cout << " transitions" << endl;
			}
		}
		emsg_input("done");
	}
}


/// Determines if initial compartment is set
void Input::set_init_c_set()
{
	for(auto &sp : model.species){
		for(auto &indd : sp.individual){
			indd.init_c_set = true;
			if(indd.enter_ref != UNSET){
				const auto &enter = sp.enter[indd.enter_ref];
				if(enter.c_set == UNSET) indd.init_c_set = false;
			}
		}
	}
}


/// Initialises local individual changes
void Input::set_local_ind_init()
{
	for(auto &sp : model.species){
		if(sp.type == INDIVIDUAL){
			auto no_enter = false;
			switch(sp.init_cond.type){
			case INIT_POP_FIXED: no_enter = true; break;
			case INIT_POP_NONE:                          // If all individual have their initial 
				no_enter = true;
				for(const auto &indd : sp.individual){
					if(indd.init_c_set == false){ no_enter = false; break;}
				}						
				break;
			case INIT_POP_DIST: break;
			}
			
			for(auto cl = 0u; cl < sp.ncla; cl++){
				auto &claa = sp.cla[cl];
				for(auto tr = 0u; tr < claa.ntra; tr++){
					const auto &tra = claa.tra[tr];
					if(tra.i != UNSET) claa.comp[tra.i].tr_leave.push_back(tr);
					if(tra.f != UNSET) claa.comp[tra.f].tr_enter.push_back(tr);
				}
				auto C = claa.ncomp;
				vector < vector <unsigned int> > tr_swap_ref; 
				vector < vector < vector <TrSwap> > > tr_swap_list;
				
				tr_swap_ref.resize(C+1);
				for(auto c = 0u; c <= C; c++) tr_swap_ref[c].resize(C+1,UNSET);
				
				vector < vector <unsigned int> > tr_swap_leave;
				tr_swap_leave.resize(C+1);
				
				for(auto c = 0u; c < C; c++){              // No events needed to go between same comp
					vector <TrSwap> vec;
					add_tr_swap(c,c,vec,tr_swap_leave,tr_swap_ref,tr_swap_list);
				}

				if(no_enter != true){
					for(auto c = 0u; c < C; c++){             // Individual entering
						vector <TrSwap> vec;
						TrSwap tsw; tsw.type = ENTER_SW; tsw.trc = c; 
						vec.push_back(tsw);
						add_tr_swap(C,c,vec,tr_swap_leave,tr_swap_ref,tr_swap_list);
					}
				}
				
				for(auto c = 0u; c < C; c++){              // Individual leaving
					vector <TrSwap> vec;
					TrSwap tsw; tsw.type = LEAVE_SW; tsw.trc = c; 
					vec.push_back(tsw);
					add_tr_swap(c,C,vec,tr_swap_leave,tr_swap_ref,tr_swap_list);
				}
				
				for(auto c = 0u; c < C; c++){              // Individual leaving end
					vector <TrSwap> vec;
					TrSwap tsw; tsw.type = LEAVE_SW; tsw.trc = c; 
					vec.push_back(tsw);
					add_tr_swap(c,C,vec,tr_swap_leave,tr_swap_ref,tr_swap_list);
				}
				
				for(auto tr = 0u; tr < claa.ntra; tr++){
					const auto &tra = claa.tra[tr];
					if(tra.type != PERIOD){
						if(tra.i == UNSET){
							vector <TrSwap> vec;
							TrSwap tsw; tsw.type = SOURCE_SW; tsw.trc = tr;
							vec.push_back(tsw);
							add_tr_swap(C,tra.f,vec,tr_swap_leave,tr_swap_ref,tr_swap_list);
						}
						else{
							if(tra.f == UNSET){
								vector <TrSwap> vec;
								TrSwap tsw; tsw.type = SINK_SW; tsw.trc = tr;
								vec.push_back(tsw);
								add_tr_swap(tra.i,C,vec,tr_swap_leave,tr_swap_ref,tr_swap_list);
							}
							else{
								vector <TrSwap> vec;
								TrSwap tsw; tsw.type = TRANS_SW; tsw.trc = tr;
								vec.push_back(tsw);
								add_tr_swap(tra.i,tra.f,vec,tr_swap_leave,tr_swap_ref,tr_swap_list);
							}
						}
					}
				}

				if(false) print_tr_swap(cl,tr_swap_ref,tr_swap_list,sp);
		
				auto tr_swap_ref_dir = tr_swap_ref;
				auto tr_swap_list_dir = tr_swap_list;
				
				auto tr_swap_leave_dir = tr_swap_leave;
				
				// Iteratively works out all the possible movements between different compartments
				do{
					auto flag = false;
					for(auto ci = 0u; ci <= C; ci++){
						for(auto j = 0u; j < tr_swap_leave[ci].size(); j++){ 
							auto cf = tr_swap_leave[ci][j];
						
							if(cf != C && cf != ci){
								if(ci >= tr_swap_ref.size()) emsg_input("prob");
								if(cf >= tr_swap_ref[ci].size()){
									emsg_input("prob2here");
								}
								
								auto ref = tr_swap_ref[ci][cf];
								if(ref != UNSET){
									for(auto k = 0u; k < tr_swap_list[ref].size(); k++){
										auto num_ev = tr_swap_list[ref][k].size();		
										if(num_ev < IND_LOCAL_TRANS_MAX){	
											for(auto cg : tr_swap_leave_dir[cf]){
												auto ref_dir = tr_swap_ref_dir[cf][cg];
												
												if(ref_dir != UNSET){
													const auto &tsw_fg = tr_swap_list_dir[ref_dir];
													for(auto j = 0u; j < tsw_fg.size(); j++){
														const auto &tsw_fg_j = tsw_fg[j];
														
														if(num_ev + tsw_fg_j.size() <= IND_LOCAL_TRANS_MAX){
															auto vec = tr_swap_list[ref][k];//tsw_if_k;
															for(const auto &tsw : tsw_fg_j) vec.push_back(tsw);
															
															bool end = false;
															if(vec.size() > 2){
																if(vec[0].type == ENTER_SW) end = true;
																if(vec[0].type == SOURCE_SW) end = true;
																if(vec[vec.size()-1].type == SINK_SW) end = true;
																if(vec[vec.size()-1].type == LEAVE_SW) end = true;
															}
														
															if(end == false){
																if(add_tr_swap(ci,cg,vec,tr_swap_leave,tr_swap_ref,tr_swap_list) == true) flag = true;
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
					if(flag == false) break;
				}while(true);
				
				if(false) print_tr_swap(cl,tr_swap_ref,tr_swap_list,sp);
		
				// Creates all potential swaps
				for(auto ci = 0u; ci <= C; ci++){
					for(auto cf : tr_swap_leave[ci]){
						auto ref = tr_swap_ref[ci][cf];
						if(ref != UNSET){
							const auto &tsw_if = tr_swap_list[ref];
							if(tsw_if.size() > 1 && !(ci == C && cf == C)){
								for(auto k = 0u; k < tsw_if.size(); k++){
									const auto &tsw_if_k = tsw_if[k];
									
									vector < vector <TrSwap> > swap;
									for(auto kk = 0u; kk < tsw_if.size(); kk++){
										if(kk != k){
											const auto &tsw_if_kk = tsw_if[kk];
											
											// Checks sequences do not share the same event
											auto fl = false;
											for(auto j = 0u; j < tsw_if_k.size(); j++){
												for(auto jj = 0u; jj < tsw_if_kk.size(); jj++){
													if(tsw_if_k[j].type == tsw_if_kk[jj].type && 
														tsw_if_k[j].trc == tsw_if_kk[jj].trc){ fl = true; break;}
												}
												if(fl == true) break;
											}
								
											// Cannot interchange between enter/source and leave/sink
											if(tsw_if_k.size() > 0 && tsw_if_kk.size() > 0){
												auto ty_k = tsw_if_k[0].type;
												auto ty_kk = tsw_if_kk[0].type;
												if(ty_k == ENTER_SW && ty_kk == SOURCE_SW) fl = true;
												if(ty_k == SOURCE_SW && ty_kk == ENTER_SW) fl = true;
												
												ty_k = tsw_if_k[tsw_if_k.size()-1].type;
												ty_kk = tsw_if_kk[tsw_if_kk.size()-1].type;
												if(ty_k == LEAVE_SW && ty_kk == SINK_SW) fl = true;
												if(ty_k == SINK_SW && ty_kk == LEAVE_SW) fl = true;
											}
									
											if(fl == false){
												swap.push_back(tsw_if_kk);
											}
										}
									}
									
									if(swap.size() > 0){
										string name;
										if(ci == C) name = "UNSET"; else name = claa.comp[ci].name;
										name += ": ";
										name += sp.tr_swap_name(cl,tsw_if_k);
										name += " ====> ";
								
										for(const auto &sw :swap){
											name += sp.tr_swap_name(cl,sw) + " | ";
										}
										
										Swap sw;
										sw.name = name;
										sw.c = ci;
										sw.start = tsw_if_k;
								
										sw.swap = swap;
										
										for(auto kk = 0u; kk < swap.size(); kk++){
											sw.swap_rep_ref.push_back(claa.swap_rep.size());
											
											SwapRep swap_rep;
											swap_rep.name = sp.tr_swap_name(cl,sw.start)+" --> "+sp.tr_swap_name(cl,swap[kk]);
											swap_rep.swap_ref = claa.swap.size();
											swap_rep.index = kk;
											claa.swap_rep.push_back(swap_rep);
										}
										
										auto hash_vec = sp.get_vec_tr_swap(ci,tsw_if_k);
										claa.hash_swap.add(claa.swap.size(),hash_vec);
										
										claa.swap.push_back(sw);
									}
								}
							}
						}
					}
				}
				
				// Sets up array which dtermines if swap is possible
				auto &sw_pos = claa.tr_swap_possible;
				sw_pos.resize(C+1);
				for(auto c = 0u; c <= C; c++){
					sw_pos[c].resize(C+1);
					for(auto cc = 0u; cc <= C; cc++){
						auto ref = tr_swap_ref[c][cc];
						if(ref != UNSET){
							if(tr_swap_list[ref].size() > 1) sw_pos[c][cc] = true;
							else sw_pos[c][cc] = false;
						}
					}
				}
			}
		}
	}
	
	if(false){
		for(auto &sp : model.species){
			for(auto cl = 0u; cl < sp.ncla; cl++){
				auto &claa = sp.cla[cl];
				auto C = claa.ncomp;
				cout << claa.swap.size() << " " << claa.swap_rep.size() << " " << C << " swap" << endl;
			}
		}
		emsg_input("gg");
	}
	
	if(false){
		cout << endl << endl << "SWAP" << endl;
		for(auto &sp : model.species){
			for(auto cl = 0u; cl < sp.ncla; cl++){
				auto &claa = sp.cla[cl];
				for(const auto &sw : claa.swap){
					cout << sw.name << endl;
				}
			}
		}
		emsg_input("done");
	}
	
	if(false){
		for(auto &sp : model.species){
			for(auto cl = 0u; cl < sp.ncla; cl++){
				auto &claa = sp.cla[cl];
				cout << endl << endl << "CLASSIFICATION " << claa.name << endl;
				for(const auto &swe : claa.swap_rep){
					cout << swe.name << endl; 
				}					
			}
		}
		emsg_input("done");
	}
}


/// Prints tr_swap (for diagnostic purposes)
void Input::print_tr_swap(unsigned int cl, const vector < vector <unsigned int> > &tr_swap_ref, const vector < vector < vector <TrSwap> > > &tr_swap_list, const Species &sp) const
{
	auto &claa = sp.cla[cl];
	auto C = claa.ncomp;
		
	cout << endl << endl <<"CLASSIFICATION " << claa.name << endl;
	for(auto ci = 0u; ci <= C; ci++){
		for(auto cf = 0u; cf <= C; cf++){
			if(ci == C) cout << "UNSET"; else cout << claa.comp[ci].name;
			cout << "->";
			if(cf == C) cout << "UNSET"; else cout << claa.comp[cf].name;
			cout << ":" << endl;
			
			auto ref = tr_swap_ref[ci][cf];
			if(ref != UNSET){
				auto &tswa = tr_swap_list[ref];
				for(auto j = 0u; j < tswa.size(); j++){
					cout << "  ";
					cout << sp.tr_swap_name(cl,tswa[j]);
					cout << endl;
				}
			}
		}
	}
	//emsg_input("done print_tr_swap");
}


/// Adds a new possibility 
bool Input::add_tr_swap(unsigned int ci, unsigned int cf, const vector <TrSwap> &vec, vector < vector <unsigned int> > &tr_swap_leave, vector < vector <unsigned int> > &tr_swap_ref, vector < vector < vector <TrSwap> > > &tr_swap_list)
{
	auto ref = tr_swap_ref[ci][cf];
	if(ref == UNSET){
		ref = tr_swap_list.size();
		tr_swap_ref[ci][cf] = ref;
		
		vector < vector <TrSwap> > add;
		tr_swap_list.push_back(add);
		
		tr_swap_leave[ci].push_back(cf);
	}
	
	//auto &trsw = tr_swap[ci][cf];
	//if(trsw.size() == 0) tr_swap_leave[ci].push_back(cf);
	auto &trsw = tr_swap_list[ref];
	
	// Checks if already exists
	auto fl = false;
	for(auto k = 0u; k < trsw.size(); k++){
		const auto &trs = trsw[k];
		if(trs.size() == vec.size()){
			auto j = 0u; 
			while(j < vec.size() && trs[j].type == vec[j].type && trs[j].trc == vec[j].trc) j++;
			if(j == vec.size()){ fl = true; break;}
		}
	}
	
	if(fl == false){
		trsw.push_back(vec);
		return true;
	}
	return false;
}


/// Sets comp_global_convert (which allows for fast conversion of global compartments)
void Input::set_comp_global_convert()
{
	for(auto &sp : model.species){
		auto num = 0u;
		for(auto cl = 0u; cl < sp.ncla; cl++) num += sp.cla[cl].ncomp;
		
		sp.comp_global_convert_set = false;
		if(1 == 0 && sp.comp_gl.size()*num < CGC_MAX){ // TURN ON
			sp.comp_global_convert_set = true;
			
			auto &cgc = sp.comp_global_convert;
			cgc.resize(sp.comp_gl.size());
			for(auto cgl = 0u; cgl < sp.comp_gl.size(); cgl++){
				cgc[cgl].resize(sp.ncla);
				for(auto cl = 0u; cl < sp.ncla; cl++){
					auto cold = sp.comp_gl[cgl].cla_comp[cl];				
					const auto &claa = sp.cla[cl];
					cgc[cgl][cl].resize(claa.comp.size());
					for(auto c = 0u; c < claa.comp.size(); c++){
						cgc[cgl][cl][c] = cgl + (c-cold)*sp.comp_mult[cl];
					}
				}
			}			
		}
	}			
}


// Sets the time varying parameter for transitions
void Input::set_tra_global_time_vari()
{
	for(auto &sp : model.species){
		for(auto &tra : sp.tra_gl){
			tra.time_vari = false;
			for(const auto &eqi : tra.dist_param){
				if(model.eqn[eqi.eq_ref].time_vari) tra.time_vari = true;
			}
			
			if(tra.bp.eq_ref != UNSET){
				if(model.eqn[tra.bp.eq_ref].time_vari) tra.time_vari = true;
			}
		}
	}
}


/// Sets flag determining if individual variation
void Input::set_ind_variation()
{
	for(auto &sp : model.species){
		sp.ind_variation = false;
		if(sp.ind_effect.size() > 0 || sp.fix_effect.size() > 0) sp.ind_variation = true;
	}
}


/// Determines if glob comp starts nm 
void Input::set_cgl_begin_nm()
{
	for(auto &sp : model.species){
		sp.cgl_begin_nm.resize(sp.comp_gl.size(),false);
		for(auto c = 0u; c < sp.comp_gl.size(); c++){
			if(sp.comp_gl[c].nmtransincomp_ref.size() > 0){
				sp.cgl_begin_nm[c] = true;
			}
		}
	}
}	
	
	
/// Sets trg reference from compartments
void Input::set_cgl_tr_source_sink()
{
	for(auto &sp : model.species){
		sp.cgl_tr_source.resize(sp.comp_gl.size(),UNSET);
		sp.cgl_tr_sink.resize(sp.comp_gl.size(),UNSET);
		for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
			const auto &tra = sp.tra_gl[tr];
			if(tra.i == UNSET) sp.cgl_tr_source[tra.f] = tr;
			if(tra.f == UNSET) sp.cgl_tr_sink[tra.i] = tr;
		}
	}
}


/// Sets transitions entering compartment
void Input::set_tr_enter()
{
	for(auto &sp : model.species){
		for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
			const auto &tra = sp.tra_gl[tr];
			if(tra.f != UNSET){
				sp.comp_gl[tra.f].tr_enter.push_back(tr);
			}
		}
	}
}
 	
	
/// Reference to see the population which causes infection 		
void Input::set_inf_cause()
{
	model.inf_cause.resize(model.species.size());
	for(auto p = 0u; p < model.species.size(); p++){
		const auto &sp = model.species[p];
		
		auto &ic = model.inf_cause[p]; ic.resize(sp.tra_gl.size());
		for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
			auto &ic_tr = ic[tr]; ic_tr.resize(model.species.size());
			
			for(auto p2 = 0u; p2 < model.species.size(); p2++){
				const auto &sp2 = model.species[p2];
				auto &ic_tr_p = ic_tr[p2]; ic_tr_p.resize(sp2.comp_gl.size());
			}
			
			const auto &tra = sp.tra_gl[tr];
			auto m = tra.markov_eqn_ref;
			if(m != UNSET){
				const auto &me = sp.markov_eqn[m];
				const auto &eq = model.eqn[me.eqn_ref];
				if(eq.infection_trans){
					for(auto p2 = 0u; p2 < model.species.size(); p2++){
						const auto &sp2 = model.species[p2];
						for(auto c = 0u; c < sp2.comp_gl.size(); c++){
							auto pref = eq.comp_pref_convert[p2][c];
							if(pref != UNSET){
								auto &inf_c = ic_tr[p2][c];
								
								inf_c.eq_ref = me.eqn_ref;
								inf_c.pref = pref;
								inf_c.po = eq.pop_ref[pref];
							}
						}
					}
				}
			}		
		}
	}
}


/// Sets the flag ppc_resample (determines if quantities are sampled or taken from posterior)
void Input::set_ppc_resample()
{
	auto st = model.details.ppc_resample;
	if(st.length() == 0) return;
	
	auto spl = split(st,',');
	for(auto i = 0u; i < spl.size(); i++){
		auto name = spl[i];
		auto fl = false;
		
		// Looks at parameters
		for(auto &pv : model.param_vec){
			if(model.param[pv.th].name == name){
				pv.ppc_resample = true;
				fl = true;
			}
		}
		
		// Looks and individual effects
		auto spl2 = split(name,'-');
		for(auto &sp : model.species){
			for(auto &ieg : sp.ind_eff_group){
				if(ieg.list.size() == spl2.size()){
					auto k = 0u;
					while(k < ieg.list.size() && find_in(spl2,ieg.list[k].name) != UNSET) k++;
					if(k == ieg.list.size()){
						ieg.ppc_resample = true;
						fl = true;
					}
				}
			}
		}
		
		if(fl == false) emsg_input("For 'resample' the value '"+name+"' is not within the model"); 
	}

	if(false){
		for(auto &pv : model.param_vec){
			cout << pv.name << " " << pv.ppc_resample << "  resamp" << endl;
		}
		
		for(const auto &sp : model.species){
			for(const auto &ieg : sp.ind_eff_group){
				for(const auto &li : ieg.list) cout << li.name << ",";
				cout << " " << ieg.ppc_resample << "  resamp" << endl;
			}
		}
		
		emsg_input("Resample");
	}
}


/// Determines if compartemnt begins or ends a transition period
void Input::set_comp_period()
{
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
		sp.period_exist = false;
		
		auto C = sp.comp_gl.size();
		sp.comp_period.resize(C);
		for(auto c = 0u; c < C; c++) sp.comp_period[c].resize(sp.ncla,false);
		
		for(const auto &tra : sp.tra_gl){
			if(tra.type == PERIOD){
				sp.period_exist = true;
				
				if(tra.i != UNSET) sp.comp_period[tra.i][tra.cl] = true;
				if(tra.f != UNSET) sp.comp_period[tra.f][tra.cl] = true;
			}
		}
	}
}


/// Determines if a compartment is terminal
void Input::set_comp_terminal()
{ 
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
		auto C = sp.comp_gl.size();
		sp.comp_terminal.resize(C,true);
		sp.comp_terminal_cl.resize(sp.ncla);
		for(auto cl = 0u; cl < sp.ncla; cl++){
			sp.comp_terminal_cl[cl].resize(C,true);
		}
		
		for(const auto &tra : sp.tra_gl){
			auto c = tra.i;
			if(c != UNSET){
				sp.comp_terminal[c] = false;
				sp.comp_terminal_cl[tra.cl][c] = false;
			}
		}
			
		if(false){
			for(auto c = 0u; c < C; c++){
				cout << sp.comp_gl[c].name <<" " <<  sp.comp_terminal[c] << " terminal" << endl;
			}
				
			for(auto cl = 0u; cl < sp.ncla; cl++){
				for(auto c = 0u; c < C; c++){
					cout << cl << " " << sp.comp_gl[c].name <<" " << sp.comp_terminal_cl[cl][c] << " terminal_cl" << endl;
				}
			}
		}
	}
}


/// Determines if a sink exists (used to stop single simulation proposals)
void Input::set_sink_exist()
{
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
		sp.sink_exist.resize(sp.ncla,false);
		for(auto cl = 0u; cl < sp.ncla; cl++){
			for(auto &tra : sp.cla[cl].tra){
				if(tra.variety == SINK_TRANS) sp.sink_exist[cl] = true;
			}
		}
		
		if(false){
			for(auto cl = 0u; cl < sp.ncla; cl++){
				cout << cl << " " << sp.sink_exist[cl] << "sink exist" << endl;
			}
			emsg_input("done");
		}
	}
}


/// Sets if ind_eff exist in equations 
void Input::set_eqn_ind_eff_exist()
{
	for(auto &eqn : model.eqn){
		switch(eqn.type){
		case BP: case SOURCE_RATE: case SOURCE_MEAN: case TRANS_RATE: case TRANS_MEAN: 
		case TRANS_NM_RATE: case TRANS_NM_MEAN: case TRANS_SHAPE: case TRANS_SCALE: case TRANS_CV:
			{
				auto p = eqn.sp_p;
				if(p == UNSET) emsg_input("Species is unset"); 
				
				const auto &sp = model.species[p];
				
				auto nie = sp.ind_effect.size();
				eqn.ind_eff_exist.resize(nie,false);
				for(auto ie : eqn.ind_eff_mult) eqn.ind_eff_exist[ie] = true;
			}
			break;

		default: break;
		}
	}		
}


/// Sets parent child relationships for parameters
void Input::set_param_parent_child()
{
	Hash hash;
	
	for(auto th = 0u; th < model.param.size(); th++){
		const auto &par = model.param[th];
		switch(par.variety){
		case DIST_PARAM: case PRIOR_PARAM:
			for(auto i = 0u; i < par.N; i++){
				auto ref = par.element_ref[i];
				if(ref != UNSET){
					const auto &ele = par.element[ref];
					if(ele.used){
						const auto &pri = model.prior[ele.prior_ref];
						for(const auto &eqi : pri.dist_param){
							add_parent_child(eqi,i,th,hash);
						}
					}
				}
			}
			break;
			
		case REPARAM_PARAM:
			for(auto i = 0u; i < par.N; i++){
				auto ref = par.element_ref[i];
				if(ref != UNSET){
					const auto &ele = par.element[ref];
					if(ele.used){
						const auto &eqi = ele.value;	
						add_parent_child(eqi,i,th,hash); 
					}
				}
			}
			break;
			
		default: break;
		}
	}
}


/// Sets which parameters are used in equations
// This writes over previously used, because of potnetial equation simplification
void Input::set_param_use()
{
	auto N = model.param.size();
	
	// Finds which parameters are actually used in the equations
	//vector < vector <bool> > eqn_used;
	//eqn_used.resize(N);
	for(auto th = 0u; th < N; th++){
		auto &par = model.param[th];
		for(auto &ele : par.element) ele.used = false;
	}
	
	for(auto &eq : model.eqn){
		for(auto &ca : eq.calc){
			for(auto &it : ca.item){
				if(it.type == PARAMETER){
					const auto &pr = eq.param_ref[it.num]; 
					if(pr.th >= model.param.size()) emsg_input("Out of range1");
					if(pr.index >= model.param[pr.th].N) emsg_input("Out of range2");
					model.param[pr.th].set_used(pr.index);
				}
			}
		}
	
		if(eq.ans.type == PARAMETER){
			if(eq.ans.num == UNSET) emsg_input("done");
			const auto &pr = eq.param_ref[eq.ans.num]; 
			model.param[pr.th].set_used(pr.index);
		}
	}
	
	for(const auto &sp : model.species){
		for(const auto &ieg : sp.ind_eff_group){
			auto N = ieg.list.size(); 
			for(auto j = 0u; j < N; j++){
				for(auto i = 0u; i < N; i++){
					model.param[ieg.omega[j][i]].set_used(0);
				}
			}
		}
		
		for(const auto &fe : sp.fix_effect){
			model.param[fe.th].set_used(0);
		}
	}
	
	if(false){
		for(auto th = 0u; th < N; th++){
			const auto &par = model.param[th];
			cout << par.name << endl;
			auto imax = par.element.size(); if(imax > 100) imax = 100;
			for(auto i = 0u; i < imax; i++) cout << par.element[i].used << " ";
			cout << endl;
		}
		emsg_input("Shows parameters used");
	}
}


/// Sets a flag if param is in omega
void Input::set_omega_fl()
{
	for(auto &pv : model.param_vec) pv.omega_fl = false;
	
	for(auto p = 0u; p < model.species.size(); p++){
		const auto &sp = model.species[p];
		for(auto g = 0u; g < sp.ind_eff_group.size(); g++){
			const auto &ieg = sp.ind_eff_group[g];
			auto N = ieg.list.size();
			for(auto j = 0u; j < N; j++){
				for(auto i = 0u; i < N; i++){
					auto &par = model.param[ieg.omega[j][i]];
					if(par.variety != CONST_PARAM){
						auto th = par.get_param_vec(0); if(th == UNSET) emsg("Should not be unset2");
						model.param_vec[th].omega_fl = true;
					}
				}
			}
		}
	}
	
	if(false){
		for(auto &pv : model.param_vec){
			cout << pv.name << " " << pv.omega_fl << endl;
		}
		emsg("done");
	}
}

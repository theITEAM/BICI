/// Functions for updating the model based on command lines in the input

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>

using namespace std;

#include "input.hh"
#include "utils.hh"

/// Adds a new species to the model
void Input::add_species(string name, SpeciesType type)
{
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
	
	Species sp(model.timepoint,model.details); sp.name = name; sp.ncla = 0; sp.type = type; 
	
	model.species.push_back(sp);	
	model.nspecies = model.species.size();  
}


/// Adds a new classification to the model
void Input::add_classification(unsigned int p, string name, string index, Coord coord)
{
	if(p == UNSET){ alert_import("Cannot add a classification when the species is not set"); return;}
	
	auto cl = find_cl(p,name);
	if(cl != UNSET){ alert_import("There is already a classification with the name '"+name+"'"); return;}
	
	if(name == index){ alert_import("Name and index must be different"); return;}
	
	Classification clas; clas.name = name; clas.index = index; clas.ntra = 0; clas.ncomp = 0; clas.coord = coord;
	
	auto &sp = model.species[p];
	
	sp.cla.push_back(clas);
	sp.ncla++;
	
	/*
	for(auto pp = 0u; pp < model.nspecies; pp++){
		if(pp != p){
			auto cl = find_cl(p,name);
			if(cl != UNSET){
				clone_compartments(pp,cl,p,model.species[p].ncla-1);
			}
		}
	}
	*/
}


/// Adds a new compartment to the model
void Input::add_compartment(string name, unsigned int p, unsigned int cl, double x, double y, double lat, double lng, bool markov_branch, string erlang_source)
{
	if(check_comp_exist(name,p) == true) return;
	
	Compartment co; 
	co.name = name; 
	co.x = x; co.y = y; co.lat = lat; co.lng = lng; 
	co.markov_branch = markov_branch;
	co.erlang_hidden = false; if(erlang_source != "") co.erlang_hidden = true;
	co.erlang_source = erlang_source;
	
	auto &claa = model.species[p].cla[cl];
	claa.comp.push_back(co);
	claa.ncomp++;
}
	

/// Adds a transition to the model
void Input::add_transition(unsigned int p, unsigned int cl, unsigned int i, unsigned int f, TransType type)
{
	auto &claa = model.species[p].cla[cl];
	auto &tra = claa.tra;

	for(auto j = 0u; j < tra.size(); j++){
		if(tra[j].i == i && tra[j].f == f){
			if(i == SOURCE){ alert_import("Cannot have two sources into the same compartment"); return;}
			else{
				if(f == SINK){ alert_import("Cannot have two sinks leaving the same compartment"); return;}
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
	tr.bp_set = false;
	tr.line_num = line_num;
	tr.erlang_hidden = false;
	
	if(i == SOURCE){ tr.variety = SOURCE_TRANS; tr.name = "+→"+claa.comp[f].name;}
	else{
		if(f == SINK){ tr.variety = SINK_TRANS; tr.name = claa.comp[i].name+"→-";}
		else{ tr.variety = NORMAL; tr.name = claa.comp[i].name+"→"+claa.comp[f].name;}
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
						if(trr.type != EXP_RATE && trr.type != ERLANG) flag = true;
					}
				}
				
				if(list.size() > 1 && (flag == true || co.markov_branch == true)) co.branch = true;
				else co.branch = false;
				
				for(auto i = 0u; i < list.size(); i++){
					claa.tra[list[i]].branch = co.branch;
				}
			}					
		}
	}
}
	

/// Determines groups of parameters which 
void Input::branch_param_group()
{
	for(auto p = 0u; p < model.nspecies; p++){
		auto &sp = model.species[p];
		
		for(auto c = 0u; c < sp.comp_gl.size(); c++){
			vector <unsigned int> group;
			for(auto i = 0u; i < sp.tra_gl.size(); i++){
				const auto &tr = sp.tra_gl[i];
				if(tr.branch == true && tr.i == c){
					const auto &eqn = model.eqn[tr.bp.eq_ref];
					if(eqn.param_ref.size() != 1){
						alert_line("'bp' should only contain a single parameter.",tr.line_num);
					}
					
					const auto &pref = eqn.param_ref[0];
					auto th = model.param[pref.th].param_vec_ref[pref.index];
					auto cl = tr.cl;
					auto &tref = sp.comp_gl[c].trainfo[cl].tra_ref;
					
					auto k = 0u; while(k < tref.size() && tref[k].tr != i) k++;
					if(k == tref.size()) emsg(" problem finding");
					
					tref[k].bp_th = th;
				
					group.push_back(th);
				}
			}
			
			if(group.size() > 0){	
				if(group.size() == 1) emsg("Cannot have one branch");
				
				auto &bp = model.branch_param;
				auto j = 0u; while(j < bp.size() && equal_vec(group,bp[j].group) == false) j++;
				
				if(j == bp.size()){
					BranchParam bpar; bpar.group = group;
					bp.push_back(bpar);
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
}


/// Adds an equation to the model
EquationInfo Input::add_equation_info(string _te, EqnType _type, unsigned int _p, unsigned int _cl)                  
{
	EquationInfo info;
	info.te_raw = _te;
	info.te = basic_equation_check(_te);
	info.type = _type;
	info.p = _p;
	info.cl = _cl;
	info.value = UNSET;
	info.eq_ref = UNSET;
	info.line_num = line_num;

	return info;
}


/// Creates equations based on EquationInfo specifications
void Input::create_equations()
{
	model.create_species_simp();
	
	for(auto th = 0u; th < model.param.size(); th++){
		auto &par = model.param[th];
	
		if(par.auto_value == true){
			for(auto &eqi : par.value){
				eqi.value = UNSET;
			}
		}
		else{
			switch(par.variety){
			case DIST_PARAM: case PRIOR_PARAM:
				for(auto i = 0u; i < par.N; i++){
					auto &pri = par.prior[i];
					for(auto &eqi : pri.dist_param){
						model.add_eq_ref(eqi);
						
						add_parent_child(eqi,i,th);
					}
				}
				break;
				
			case REPARAM_PARAM:
				for(auto i = 0u; i < par.N; i++){
					auto &eqi = par.value[i];
					model.add_eq_ref(eqi);	
					add_parent_child(eqi,i,th);
				}
				break;
				
			case CONST_PARAM: 
				for(auto &eqi : par.value){
					eqi.value = number(eqi.te);
					if(eqi.value == UNSET){
						emsg("Value not set");
					}
				}
				break;
				
			case UNSET_PARAM: emsg("Should not be unset7"); break;
			}
		}
	}
	
	for(auto &der : model.derive){                       // Derived quantities
		for(auto &eq : der.eq) model.add_eq_ref(eq);
	}
	
	for(auto &sp : model.species){                       // Transitions
		for(auto &tr_gl : sp.tra_gl){ 
			if(tr_gl.branch == true) model.add_eq_ref(tr_gl.bp);
			for(auto &dp : tr_gl.dist_param) model.add_eq_ref(dp);
		}
	}
	
	for(auto &sp : model.species){                       // Observation model
		for(auto &ds : sp.source){                     
			auto &om = ds.obs_model;
			if(om.Se.te != "") model.add_eq_ref(om.Se);
			if(om.Sp.te != "") model.add_eq_ref(om.Sp);
			if(om.trapprob.te != "") model.add_eq_ref(om.trapprob);
		}
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
				alert_line("A population-based model cannot have fixed effects such as '["+fe.name+"]'",fe.line_num); 
			}
		}
		
	
	}
	
	for(auto eq : model.eqn){
		if(eq.warn != "") alert_line("\nFor equation '"+eq.te_raw+"':\n "+eq.warn,eq.line_num); 
	}
}


/// Adds a parent child relationship if there is a reparameterisation or distribution
void Input::add_parent_child(const EquationInfo eqi, unsigned int i, unsigned int th)
{
	ParamRef parref; parref.th = th; parref.index = i;
					
	auto &par = model.param[th];
	
	auto ref = eqi.eq_ref;
	if(ref != UNSET){
		const auto &eq = model.eqn[ref];

		if(eq.param_ref.size() > 0){								
			for(auto pr : eq.param_ref){
				add_to_list(par.parent[i],pr);
				add_to_list(model.param[pr.th].child[pr.index],parref);
			}
		}
	}
}

						
/// Works out which equations are Markovian
void Input::create_markov_eqn()
{
	for(auto e = 0u; e < model.eqn.size(); e++){
		auto &eqn = model.eqn[e];
		if(eqn.type == TRANS_RATE || eqn.type == SOURCE_RATE){
			if(eqn.sp_p == UNSET) alert_import("Equation problem"); 
			else{
				auto &sp = model.species[eqn.sp_p];
				
				eqn.markov_eqn_ref = sp.markov_eqn.size();
				
				MarkovEqn me; 
				me.eqn_ref = e; 
		
				me.time_vari = false;
				if(eqn.pop_ref.size() > 0) me.time_vari = true;
				for(const auto &pr : eqn.param_ref){
					if(model.param[pr.th].time_dep == true) me.time_vari = true;
				}
			
				const auto &eq = model.eqn[e];	
				me.source = false; 
				if(eq.type == SOURCE_RATE){
					me.source = true;	
					me.source_tr_gl = eq.source_tr_gl;
				}
		
				me.ind_eff_mult = eq.ind_eff_mult;
				me.fix_eff_mult = eq.fix_eff_mult;
				me.ind_variation = false;
				if(me.ind_eff_mult.size() > 0 || me.fix_eff_mult.size() > 0) me.ind_variation = true;
		
				sp.markov_eqn.push_back(me);
			}
		}
	}
	
	// Add a reference from global transitions to markovian equations
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];
		
		if(sp.type == INDIVIDUAL){
			for(auto &tr : sp.tra_gl){
				tr.markov_eqn_ref = UNSET;
				if(tr.type == EXP_RATE || tr.type == EXP_MEAN){
					auto e = tr.dist_param[0].eq_ref; if(e == UNSET) emsg("Prob");
					const auto &eq = model.eqn[e];
					auto mef = eq.markov_eqn_ref; if(mef == UNSET) emsg("Prob");
					tr.markov_eqn_ref = mef;
				}
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
					if(eq == UNSET) emsg("Problem with ME");
				
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
		auto N = 1u;                                                  // Works out global conversion
		for(auto cl = 0u; cl < sp.ncla; cl++){
			sp.comp_mult.push_back(N);
			N *= sp.cla[cl].comp.size();
		}
		
		sp.N = N;
	
		for(auto c = 0u; c < N; c++){                                 // Constructs global compartments
			CompGlobal co;
			co.erlang_hidden = false;
			
			string na = "";
			for(auto cl = 0u; cl < sp.ncla; cl++){
				auto &claa = sp.cla[cl];
				auto cc = (c/sp.comp_mult[cl])%claa.comp.size();
				co.cla_comp.push_back(cc);
				if(claa.comp[cc].erlang_hidden == true) co.erlang_hidden = true;
				
				if(na != "") na += "|";
				na += claa.comp[cc].name;
				
				CompGlTrans ti; ti.branch = claa.comp[cc].branch;
				co.trainfo.push_back(ti);
			}
			co.name = na;
			
			sp.comp_gl.push_back(co);
		}

		for(auto cl = 0u; cl < sp.ncla; cl++){                              // Constructs global transitions
			const auto &cla = sp.cla[cl];
			for(auto j = 0u; j < cla.ntra; j++){
				const auto &tr = cla.tra[j];

				for(auto c = 0u; c < N; c++){
					const auto &co = sp.comp_gl[c];		

					auto cc = co.cla_comp[cl];
					auto ii = UNSET, ff = UNSET;

					if(tr.i == SOURCE){
						if(tr.f == cc){ ii = SOURCE; ff = c;}
					}
					else{
						if(tr.f == SINK){
							if(tr.i == cc){ ii = c; ff = SINK;}
						}
						else{
							if(tr.i == cc){
								int dif = (int)tr.f - (int)tr.i;
								ii = c; ff = c + dif*sp.comp_mult[cl];
							}
						}
					}

					if(ii != UNSET){
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
						tr_gl.variety = tr.variety;
						tr_gl.bp_set = tr.bp_set;
						tr_gl.branch = tr.branch;

						tr_gl.bp = tr.bp;
						tr_gl.dist_param = tr.dist_param;
						tr_gl.transform.resize(sp.comp_gl.size(),UNSET);
						tr_gl.line_num = tr.line_num;
						
						vector <DepConv> dep_conv;
						for(auto cl2 = 0u; cl2 < sp.ncla; cl2++){
							if(cl2 != cl){
								DepConv dc; 
								dc.before = sp.cla[cl2].index;
								dc.after = sp.cla[cl2].comp[co.cla_comp[cl2]].name;
								dep_conv.push_back(dc);
							}
						}

						auto res = swap_index(tr_gl.bp.te,dep_conv); 
						if(res.warn != ""){ alert_equation(tr_gl.bp,res.warn); return;}
 
						for(auto i = 0u; i < tr_gl.dist_param.size(); i++){
							res = swap_index(tr_gl.dist_param[i].te,dep_conv); 	
							if(res.warn != ""){ alert_equation(tr_gl.dist_param[i],res.warn); return;}
						}
 
						if(ii != SOURCE){
							auto &ti = sp.comp_gl[ii].trainfo[cl];
							if(ti.branch != tr.branch) emsg("Branch does not agree");
							TransBPRef tbpref; tbpref.tr = sp.tra_gl.size(); tbpref.bp_th = UNSET;
							ti.tra_ref.push_back(tbpref);
						}
						
						sp.tra_gl.push_back(tr_gl);
					}
				}
			}
		}
	}
	
	/// Calculates "transform" allows for transtiont o be changed based on a transition in another cl
	// Either the initial of final compartment cna be used for the transformation
	for(auto &sp : model.species){
		for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
			auto &t1 = sp.tra_gl[tr];
			for(auto tr2 = tr; tr2 < sp.tra_gl.size(); tr2++){
				auto &t2 = sp.tra_gl[tr2];
				if(t1.cl == t2.cl && t1.tr == t2.tr){
					if(t2.i != SOURCE) t1.transform[t2.i] = tr2;
					if(t2.f != SINK) t1.transform[t2.f] = tr2;
					if(t1.i != SOURCE) t2.transform[t1.i] = tr;
					if(t1.f != SINK) t2.transform[t1.f] = tr;
				}
			}
		}
		
		if(false){
			for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
				auto &t1 = sp.tra_gl[tr];
				cout << t1.name << ":" << endl;
				for(auto c = 0u; c < sp.comp_gl.size(); c++){
					if(t1.transform[c] != UNSET){
						auto &t2 = sp.tra_gl[t1.transform[c]];
						cout << t2.name << ",";
					}
				}
				cout << endl;
			}
		}
	}
}


/// Divides an Erland distribution into multiple compartments
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
								auto name = get_erlang_name(cname,claa.comp);
								add_compartment(name,p,cl,UNSET,UNSET,UNSET,UNSET,false,cname);
							}
							
							add_transition(p,cl,ci,cnew,EXP_RATE);
							
							auto &tra2 = claa.tra[claa.ntra-1];
							tra2.dist_param.push_back(add_equation_info(rate_new,TRANS_RATE,p,cl));
							tra2.line_num = lnum;
							if(j < shapei) tra2.erlang_hidden = true;
							else tra2.name = tname;
							
							string bp = "";
							if(bp != "") tra2.bp_set = true;
							if(bp == "") bp = "b";
							tra2.bp = add_equation_info(bp,BP,p,cl);
							
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


/// Find a new compartmental name for an erlang compartment
string Input::get_erlang_name(string name, const vector <Compartment> &comp) const
{
	auto num = 1u;
	do{
		auto na = name+tstr(num);
		auto c = 0u; while(c < comp.size() && comp[c].name != na) c++;
		if(c == comp.size()) return na;
		num++;
	}while(true);
}


// Combines branching probabilities with rates in population model
void Input::population_bp_rate_combine()
{
	for(auto p = 0u; p < model.nspecies; p++){
		auto &sp = model.species[p];
		if(sp.type == POPULATION){
			for(auto cl = 0u; cl < sp.ncla; cl++){
				auto &claa = sp.cla[cl];
				for(auto tr = 0u; tr < claa.ntra; tr++){
					auto &tra = claa.tra[tr];
					
					if(tra.branch == true){
						if(tra.type != EXP_RATE) emsg("Should be rate");
						
						auto &dp = tra.dist_param[0];
						
						auto rate = dp.te_raw;
						auto bp = tra.bp.te_raw;
						string prod;
						if(rate.length() > 4 && rate.substr(0,3) == "1/(" && rate.substr(rate.length()-1,1) == ")"){
							prod = bp+rate.substr(1); 
						}
						else{
							prod = bp+"*("+rate+")";
						}
						dp = add_equation_info(prod,TRANS_RATE,dp.p,dp.cl);
					}
				}
			}
		}
	}
}


// Attaches compartments to source equations
void Input::source_equation_comp()
{
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
			const auto &tra = sp.tra_gl[tr];
			if(tra.type == EXP_RATE){
				auto &eq = model.eqn[tra.dist_param[0].eq_ref];
				if(eq.type == SOURCE_RATE){
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
		
	while(t < de.t_end-TINY){
		model.timepoint.push_back(t);
		t += de.dt;
	}
	model.timepoint.push_back(de.t_end);
	model.ntimepoint = model.timepoint.size();
	
	if(false){
		for(auto i = 0u; i < model.timepoint.size(); i++){
			cout << i << " "<<  model.timepoint[i] << " timepoint" << endl;
		}
	}
}


/// Creates splices for the relavent parameters
void Input::create_spline()
{
	const auto &tp = model.timepoint;
	auto ntp = tp.size();
	
	for(auto th = 0u; th < model.param.size(); th++){
		const auto &par = model.param[th];
		
		if(par.time_dep == true){
			auto &sinfo = par.spline_info;
			if(sinfo.on != true) emsg("Problem loading spline");
			
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
						if(k == model.param_vec.size()) emsg("Could not find");
						
						param_ref.push_back(k);
					}
					
					spl.param_ref = param_ref;
				}
							
				auto i = 0u; // Indexes times
				for(auto ti = 0u; ti < ntp-1; ti++){
					auto tmid = (tp[ti]+tp[ti+1])/2;
					while(i+1 < times.size() && times[i+1] < tmid) i++;
					if(i+1 >= times.size()) emsg("Range problem");
					
					auto f = (times[i+1]-tmid)/(times[i+1]-times[i]);
					if(spl.constant == false){
						SplineDiv sd;	sd.th1 = param_ref[i]; sd.th2 = param_ref[i+1]; sd.f = f;	
						spl.div.push_back(sd);
					}
					else{
						auto val1 = par.value[j*ntimes+i].value;
						auto val2 = par.value[j*ntimes+i+1].value;
						spl.const_val.push_back(val1*f+val2*(1-f));
					}
				}
				
				model.spline.push_back(spl);
 			}
		}			
	}
	
	/// Converts equations such that the spline is referenced
	const auto &spl = model.spline;			
	for(auto &eq : model.eqn){
		for(auto &ca : eq.calc){
			for(auto &it : ca.item){
				if(it.type == SPLINE){
					auto th = it.num, index = it.index;
					auto ntimes = model.param[th].spline_info.knot_times.size();
					if(index%ntimes != 0) emsg("Should be zero");
					index /= ntimes;
					auto i = 0u; while(i < spl.size() && !(spl[i].th == th && spl[i].index == index)) i++;
					if(i == spl.size()) emsg("Spline could not be found");
					it.num = i; it.index = UNSET;
				}
			}
		}
	
		if(eq.ans.type == SPLINE){
			auto th = eq.ans.num, index = eq.ans.index;
			auto ntimes = model.param[th].spline_info.knot_times.size();
			if(index%ntimes != 0) emsg("Should be zero");
			index /= ntimes;
			auto i = 0u; while(i < spl.size() && !(spl[i].th == th && spl[i].index == index)) i++;
			if(i == spl.size()) emsg("Spline could not be found");
			eq.ans.num = i; eq.ans.index = UNSET;
		}
	}
}


/// If any parameters are called "dist" then put values for the distance matrix into them
void Input::set_dist(Param &par)
{
	if(par.dep.size() != 2){	
		alert_import("To be a distance matrix the parameter must have two dependencies."); 
		return;
	}
	
	auto index = remove_prime(par.dep[0].index);
	
	if(index != remove_prime(par.dep[1].index)){
		alert_import("To be a distance matrix the dependencies have to both refer to the same classification."); 
		return;
	}

	unsigned int p, cl = UNSET;
	for(p = 0; p < model.species.size(); p++){
		const auto &sp = model.species[p];
		for(cl = 0; cl < sp.ncla; cl++){
			if(sp.cla[cl].index == index) break;
		}
		if(cl < sp.ncla) break;
	}
	
	if(p == model.species.size()) emsg("Classification is undefined");
	
	const auto &claa = model.species[p].cla[cl];
	
	auto N = claa.ncomp;

	vector <unsigned int> ind(2);

	switch(claa.coord){
	case CARTESIAN:
		{
			vector <double> px, py;
	
			for(auto i = 0u; i < N; i++){
				px.push_back(claa.comp[i].x); py.push_back(claa.comp[i].y); 
			}
			
			for(auto j = 0u; j < N; j++){
				for(auto i = j; i < N; i++){
					auto val = 0.0;
					if(i != j){
						auto dx = px[i]-px[j], dy = py[i]-py[j];
						val = sqrt(dx*dx + dy*dy);
					}
					
					ind[0] = j; ind[1] = i; set_element(par.value,par.dep,ind,tstr(val));
					ind[0] = i; ind[1] = j; set_element(par.value,par.dep,ind,tstr(val));
				}
			}
		}
		break;
		
	case LATLNG:
		{
			vector <double> lat, lng;
			
			for(auto i = 0u; i < N; i++){
				lat.push_back(claa.comp[i].lat*M_PI/180.0);
				lng.push_back(claa.comp[i].lng*M_PI/180.0);
			}
			
			auto r = 6371.0;

			for(auto j = 0u; j < N; j++){
				for(auto i = j; i < N; i++){
					auto val = 0.0;
					if(i != j){
						auto si = sin(0.5*(lat[j]-lat[i]));
						auto si2 = sin(0.5*(lng[j]-lng[i]));
						val = 2*r*asin(sqrt(si*si + cos(lat[i])*cos(lat[j])*si2*si2));
					}
					
					ind[0] = j; ind[1] = i; set_element(par.value,par.dep,ind,tstr(val));
					ind[0] = i; ind[1] = j; set_element(par.value,par.dep,ind,tstr(val));
				}
			}
		}
		break;
	}
	par.variety = CONST_PARAM;
}


/// Creates a list of non-Markovian transitions
void Input::create_nm_trans()
{			
	for(auto &sp : model.species){
		auto &nm = sp.nm_trans;
	
		for(auto &trg : sp.tra_gl){
			auto type = trg.type;
			if(type == EXP_RATE || type == EXP_MEAN){
				trg.nm_ref = UNSET;
			}
			else{
				NMTrans nmtra;
				nmtra.type = type;
				nmtra.bp_eq_ref = UNSET; 
				if(trg.branch == true) nmtra.bp_eq_ref = trg.bp.eq_ref;
				for(auto i = 0u; i < trg.dist_param.size(); i++){
					nmtra.dist_param_eq_ref.push_back(trg.dist_param[i].eq_ref);
				}
				
				auto k = 0u;
				while(k < nm.size() && !(nm[k].type == nmtra.type && nm[k].bp_eq_ref == nmtra.bp_eq_ref && equal_vec(nm[k].dist_param_eq_ref,nmtra.dist_param_eq_ref) == true)) k++;
				
				if(k == nm.size()) nm.push_back(nmtra);
				
				trg.nm_ref = k;
			}
		}
	}

}


/// References markov eqns and nm_trans ind_eff_group
void Input::ind_eff_group_trans_ref()
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
				if(nmt.bp_eq_ref != UNSET) eq_list.push_back(nmt.bp_eq_ref);
				for(auto eq : nmt.dist_param_eq_ref) eq_list.push_back(eq);
				
				for(auto eq : eq_list){
					for(auto ie : model.eqn[eq].ind_eff_mult){
						if(find_in(ie_list,ie) != UNSET){
							ieg.nm_trans_ref.push_back(e);
							add_to_vec(sp.ind_effect[ie].nm_trans_ref,e);
							break;
						}
					}
				}
			}
		}
		
		for(auto e = 0u; e < sp.markov_eqn.size(); e++){
			for(auto ie : sp.markov_eqn[e].ind_eff_mult){
				add_to_vec(sp.ind_effect[ie].markov_eqn_ref,e);
			}
		}
	}
}


/// References populations in individual effect
void Input::ind_eff_pop_ref()
{
	for(auto k = 0u; k < model.pop.size(); k++){
		const auto &pop = model.pop[k];
		for(auto ie : pop.ind_eff_mult){
			auto &sp = model.species[pop.sp_p];
			add_to_vec(sp.ind_effect[ie].pop_ref,k);
		}
	}
}




/// Works out how markov eqns change on c->c'  
void Input::create_markov_comp_gl()
{
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp = model.species[p];

		auto N = sp.comp_gl.size();
		for(auto c = 0u; c < N; c++){
			auto &co = sp.comp_gl[c];
			
			for(auto cl = 0u; cl < sp.ncla; cl++){
				const auto &tinfo = co.trainfo[cl];
		
				const auto &tref = tinfo.tra_ref;
				for(const auto &tre : tref){
					auto tgl = tre.tr;
					const auto &tr = sp.tra_gl[tgl];
					if(tr.type == EXP_RATE){
						if(tr.markov_eqn_ref == UNSET) emsg("Should not be unset4");
						co.me_ref.push_back(tr.markov_eqn_ref);
					}
					else{
						if(tr.markov_eqn_ref != UNSET) emsg("Should not be set3");
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
					IslandComp ico; ico.c = c;
					isl.comp.push_back(ico);
					
					island.push_back(isl);
				}
					
				// Joins together islands based on transitions			
				for(auto tr = 0u; tr < claa.ntra; tr++){
					const auto &tra = claa.tra[tr];
					if(tra.variety != SOURCE_TRANS){
						auto c1 = tra.i;
						auto c2 = tra.f;
						
						unsigned int isl1, j = UNSET; 
						for(isl1 = 0u; isl1 < island.size(); isl1++){
							const auto &co = island[isl1].comp;
							j = 0; while(j < co.size() && co[j].c != c1) j++;
							if(j < co.size()) break;
						}
						if(isl1 == island.size()) emsg("Cannot find island1");
							
						IslandTrans itr;
						itr.tr = tr;
						island[isl1].comp[j].leave.push_back(itr);
					
						auto isl2 = UNSET;
						if(c2 != SINK){
							for(isl2 = 0u; isl2 < island.size(); isl2++){
								const auto &co = island[isl2].comp;
								j = 0; while(j < co.size() && co[j].c != c2) j++;
								if(j < co.size()) break;
							}
							if(isl2 == island.size()) emsg("Cannot find island2");
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
						for(auto &le : co.leave){
							auto tr = le.tr;
							const auto &tra = claa.tra[tr];
							
							// Sets cf which gives the final
							auto f = tra.f;
							auto j = 0u; while(j < isl.comp.size() && isl.comp[j].c != f) j++;
							if(j == isl.comp.size()) emsg("problem");
							le.cf = j;
							
							// Given a global compartment c  works put the me which should be used
							// Note, the contribution to in cl is ignored )hence this account for trans in other cl
							vector <unsigned int> me_ref(sp.comp_gl.size(),UNSET);
						
							for(auto k = 0u; k < sp.tra_gl.size(); k++){
								const auto &trg = sp.tra_gl[k];
								
								if(trg.cl == cl && trg.tr == tr){
									auto me = trg.markov_eqn_ref;
									if(me == UNSET) emsg("Should not be unset5");
									
									auto cc = trg.i - tra.i*sp.comp_mult[cl];
									for(auto c = 0u; c < claa.comp.size(); c++){
										me_ref[cc+c*sp.comp_mult[cl]] = me;
									}	
								}
							}
							
							for(auto c = 0u; c < sp.comp_gl.size(); c++){
								if(me_ref[c] == UNSET) emsg("me ref prob");
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
					emsg("Islands");
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
						
						auto k = par.param_vec_ref[ind];	
						if(k != UNSET){
							if(par.spline_info.on == true){
								if(ind != 0) emsg("Index problem");
								
								const auto &dep = par.dep[par.dep.size()-1];
								if(dep.mult != 1) emsg("mult prob");
								
								auto s = 0u; while(s < model.spline.size() && th != model.spline[s].th) s++;
								if(s == model.spline.size()) emsg("Cannot find spline");
								
								const auto &spl = model.spline[s];
								if(spl.div.size() != T) emsg("Wrong size1");
							
								for(auto t = 0u; t < dep.list.size(); t++){ // Goes down the spline
									k = par.param_vec_ref[ind+t];
							
									AffectLike al; al.map.resize(T,false);
									for(auto ti = 0u; ti < T; ti++){ 
										if(spl.div[ti].th1 == k || spl.div[ti].th2 == k) al.map[ti] = true;
									}
									
									al.type = SPLINE_AFFECT; al.num = s; al.num2 = UNSET;
									param_vec_add_affect(model.param_vec[k].affect_like,al);		
									
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
			break;
			
		case POPULATION:
			{
				for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
					const auto &tra = sp.tra_gl[tr];
					const auto &eqn = model.eqn[tra.dist_param[0].eq_ref];

					for(auto &pref : eqn.param_ref){
						auto th = pref.th, ind = pref.index;
					
						auto &par = model.param[th];
						auto k = par.param_vec_ref[ind];	
						if(k != UNSET){
							if(par.spline_info.on == true){
								const auto &dep = par.dep[par.dep.size()-1];
								if(dep.mult != 1) emsg("mult prob");
								
								auto s = 0u; while(s < model.spline.size() && th != model.spline[s].th) s++;
								if(s == model.spline.size()) emsg("Cannot find spline");
								
								const auto &spl = model.spline[s];
								if(spl.div.size() != T) emsg("Wrong size1");
							
								for(auto t = 0u; t < dep.list.size(); t++){ // Goes down the spline
									k = par.param_vec_ref[ind+t];
							
									AffectLike al; al.map.resize(T,false);
									for(auto ti = 0u; ti < T; ti++){ 
										if(spl.div[ti].th1 == k || spl.div[ti].th2 == k) al.map[ti] = true;
									}
									
									al.type = SPLINE_AFFECT; al.num = s; al.num2 = UNSET;
									param_vec_add_affect(model.param_vec[k].affect_like,al);		
									
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
			break;
		}
	}
	
	
	/// Effect of parameters on the observation mode
	for(auto eq = 0u; eq < model.eqn.size(); eq++){
		const auto &eqn = model.eqn[eq];
		
		switch(eqn.type){			
		case SE: case SP: case TRAP_PROB:
			for(auto &pref : eqn.param_ref){
				auto th = pref.th, ind = pref.index;
			
				auto &par = model.param[th];
				auto k = par.param_vec_ref[ind];	
				if(k != UNSET){
					AffectLike al; al.type = LIKE_OBS_AFFECT; al.num = eqn.sp_p; al.num2 = UNSET;
					param_vec_add_affect(model.param_vec[k].affect_like,al);
				}
			}
			break;
			
		default: break;
		}
	}
	
	/// Works out how the parameter affect the prior
	for(auto th = 0u; th < model.nparam_vec; th++){
		const auto &pv = model.param_vec[th];
		const auto &par = model.param[pv.th];
		if(par.variety == PRIOR_PARAM || par.variety == DIST_PARAM){
			AffectLike al;  al.num = th; al.num2 = UNSET; 
			if(par.variety == PRIOR_PARAM) al.type = PRIOR_AFFECT; else al.type = DIST_AFFECT;
			
			param_vec_add_affect(model.param_vec[th].affect_like,al);
			
			for(const auto &pref : par.parent[pv.index]){
				const auto th_par = model.param[pref.th].param_vec_ref[pref.index];				
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
					if(par.name != fe.name) emsg("names do not match");
					if(par.N != 1) emsg("Should be univariate");
					auto k = par.param_vec_ref[0];
					
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
	

	/// Looks at the effect of fixed effects in populations
	for(auto po = 0u; po < model.pop.size(); po++){
		const auto &pop = model.pop[po];
		auto p = pop.sp_p;
		const auto &sp = model.species[p];
		
		if(sp.type == INDIVIDUAL){
			for(auto f : pop.fix_eff_mult){		
				const auto &fe = sp.fix_effect[f];
				const auto &par = model.param[fe.th];
				if(par.name != fe.name) emsg("names do not match");
				if(par.N != 1) emsg("Should be univariate");
				auto k = par.param_vec_ref[0];
					
				AffectLike al; 	
				al.type = EXP_FE_AFFECT; al.num = p; al.num2 = f;
				param_vec_add_affect(model.param_vec[k].affect_like,al);		
				
				add_pop_affect(k,po);
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
					auto &eqn = model.eqn[nmt.dist_param_eq_ref[k]];
					for(const auto &pref : eqn.param_ref){
						auto th = pref.th, ind = pref.index;
				
						auto &par = model.param[th];
						auto k = par.param_vec_ref[ind];	
						if(k != UNSET){
							AffectLike al; al.type = NM_TRANS_AFFECT; al.num = p; al.num2 = i;
							param_vec_add_affect(model.param_vec[k].affect_like,al);		
						}
					}
					
					/// Looks at the effects of fixed effects on transitions
					for(auto j = 0u; j < eqn.fix_eff_mult.size(); j++){
						auto f = eqn.fix_eff_mult[j];
					
						const auto &fe = sp.fix_effect[f];
					
						const auto &par = model.param[fe.th];
						if(par.name != fe.name) emsg("names do not match");
						if(par.N != 1) emsg("Should be univariate");
						auto k = par.param_vec_ref[0];
					
						if(k != UNSET){				
							AffectLike al; 	
							al.type = EXP_FE_AFFECT; al.num = p; al.num2 = f;
							param_vec_add_affect(model.param_vec[k].affect_like,al);		
						
							al.type = NM_TRANS_AFFECT; al.num = p; al.num2 = i;
							param_vec_add_affect(model.param_vec[k].affect_like,al);		
						}
					}
				}
			}
		}
	}
	
	/// Variance and correlation parameters which affect ie likelihood
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
							auto k = par.param_vec_ref[0];
							if(k == UNSET) emsg("Should not be unset6");
							
							AffectLike al; 	
							al.type = OMEGA_AFFECT; al.num = p; al.num2 = g;
							param_vec_add_affect(model.param_vec[k].affect_like,al);	
							
							if(j == i){  // If on the di
								auto ie = ieg.list[i].index;
								
								auto &ind_eff = sp.ind_effect[ie];
								
								al.type = EXP_IE_AFFECT; al.num = p; al.num2 = ie;
								param_vec_add_affect(model.param_vec[k].affect_like,al);	
								
								if(ind_eff.markov_eqn_ref.size() > 0){ 
									al.type = INDFAC_INT_AFFECT; al.num = p; al.num2 = UNSET;
									param_vec_add_affect(model.param_vec[k].affect_like,al);		
								
									for(auto e : ind_eff.markov_eqn_ref){	
										AffectLike al; al.type = MARKOV_LIKE_AFFECT; al.num = p; al.num2 = e;
										auto &me = model.species[p].markov_eqn[e];
										if(me.time_vari == true) al.map.resize(T,true);
										else al.map.resize(1,true);
							
										param_vec_add_affect(model.param_vec[k].affect_like,al);	
									}
								}
								
								for(auto e : ind_eff.nm_trans_ref){	
									AffectLike al; 
									al.type = NM_TRANS_AFFECT; al.num = p; al.num2 = e;
									param_vec_add_affect(model.param_vec[k].affect_like,al);	
								}
								
								for(auto po : ind_eff.pop_ref){	
									add_pop_affect(k,po);
								}
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
}


/// Incorporate the effect of altering population po for parameter k
void Input::add_pop_affect(unsigned int k, unsigned int po)
{
	const auto &pop = model.pop[po];
	auto T = model.ntimepoint-1;

	AffectLike al; 	
	al.type = POP_AFFECT; al.num = UNSET; al.num2 = UNSET; 
	al.map.resize(model.pop.size(),false); al.map[po] = true;
	param_vec_add_affect(model.param_vec[k].affect_like,al);		
	
	for(const auto &mer : pop.markov_eqn_ref){
		AffectLike al; al.type = DIV_VALUE_AFFECT;
		al.num = mer.p; al.num2 = mer.e; 
		
		auto &me = model.species[mer.p].markov_eqn[mer.e];
		if(me.time_vari == true) al.map.resize(T,true);
		else al.map.resize(1,true);
		
		param_vec_add_affect(model.param_vec[k].affect_like,al);		
			
		al.type = MARKOV_LIKE_AFFECT; 
		param_vec_add_affect(model.param_vec[k].affect_like,al);	
	}
	
	for(const auto &tref : pop.trans_ref){
		AffectLike al; al.map.resize(T,true);
		al.type = MARKOV_POP_AFFECT; al.num = tref.p; al.num2 = tref.tr;
		param_vec_add_affect(model.param_vec[k].affect_like,al);
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
			const auto &tra =  sp.tra_gl[trg];
			auto tr = tra.tr;
			auto cl = tra.cl;
			
			auto c = tra.i; if(c == SOURCE) c = tra.f;
			
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

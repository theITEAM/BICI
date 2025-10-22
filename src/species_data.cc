// Initialises species data structure from sources

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
using namespace std;

#include "utils.hh"
#include "species.hh"

/// Initialises the data sources
void Species::initialise_data(Operation mode)
{
	T = details.T;
	
	init_cond.type = INIT_POP_NONE;
	
	Hash hash_pop_filter; 
	Hash hash_pop_trans_filter; 
	
	print_diag("load data start");
	
	for(const auto &so : source){
		print_diag("load "+so.name);
		
		switch(so.cname){
		case INIT_POP: init_pop_data(so); break;
		case ADD_POP: add_pop_data(so,1); break;
		case REMOVE_POP: add_pop_data(so,-1); break;
		case ADD_IND: add_ind_data(so); break;
		case REMOVE_IND: remove_ind_data(so); break;
		case MOVE_IND: move_ind_data(so); break;
		case COMP_DATA: comp_data(so); break;
		case TEST_DATA: test_data(so); break;
		case POP_DATA: population_data(so,hash_pop_filter); break;
		case TRANS_DATA: trans_data(so); break;
		case POP_TRANS_DATA: popu_trans_data(so,hash_pop_trans_filter); break;
		case GENETIC_DATA: genetic_data(so); break;
		default:
			emsg_input("data type not added:"); break;
		}
	}

	print_diag("load data end");
	
	if(type == INDIVIDUAL) jiggle_data(mode);
	
	nindividual_in = individual.size();

	add_unobs_Amatrix_ind();

	X_vector_order();
	
	order_data_events();
	
	set_default_enter();
	
	set_add_rem_pop_change();
	
	set_ind_tmin_tmax();
}


/// Sets the timerange over which measurements are made on individuals
void Species::set_ind_tmin_tmax()
{
	for(auto &ind : individual){
		auto tmin = LARGE;
		auto tmax = -LARGE;
		for(const auto &ev : ind.ev){
			auto t = ev.tdiv;
			if(t < tmin) tmin = t;
			if(t > tmax) tmax = t;
		}
			
		for(const auto &ob : ind.obs){
			if(!(ob.type == OBS_COMP_EV && ob.not_alive)){
				auto t = ob.tdiv;
				if(t < tmin) tmin = t;
				if(t > tmax) tmax = t;
			}
		}
		 
		ind.tdivmin = tmin; ind.tdivmax = tmax;
	}
}



/// Used to order genetic data
bool Fraction_ord (const Fraction &frac1, const Fraction &frac2)                      
{ return (frac1.fr > frac2.fr); };  


/// init-pop command
void Species::init_pop_data(const DataSource &so)
{
	const auto &tab = so.table;
	
	if(init_cond.type != INIT_POP_NONE) alert_source("Cannot set initial conditions twice",so);
	
	auto foc_cl = so.focal_cl; 
	
	init_cond.type = so.init_pop_type;
	init_cond.focal_cl = foc_cl;
	
	const auto &comp = comp_gl;
	switch(init_cond.type){
	case INIT_POP_FIXED:
		if(foc_cl != UNSET){	
			auto per = pop_get_percentage(so);
		
			if(false){
				for(auto cl = 0u; cl < ncla; cl++){
					const auto &claa = cla[cl];	
					for(auto c = 0u; c < claa.ncomp; c++){
						cout << claa.comp[c].name << " " << per[cl][c] << " Percent" << endl;
					}
				}
			}
			
			init_cond.cnum.resize(comp.size(),0);
			
			if(false){ // This just rounds based on fraction
				for(auto c = 0u; c < comp.size(); c++){
					const auto &co = comp[c];
					
					auto fac = 1.0;
					for(auto cl = 0u; cl < ncla; cl++){
						if(cl == foc_cl){
							fac *= per[cl][co.cla_comp[cl]];
						}
						else{
							fac *= per[cl][co.cla_comp[cl]]/100;
						}
					}
					init_cond.cnum[c] = fac;
				}
			}
			else{  // This gets populations right in the focal classification
				const auto &claa = cla[foc_cl];
				
				for(auto j = 0u; j < claa.comp.size(); j++){
					auto num = per[foc_cl][j];
			
					auto n = 0u;
			
					vector <Fraction> frac_list;
					for(auto c = 0u; c < comp.size(); c++){
						const auto &co = comp[c];
						
						if(co.cla_comp[foc_cl] == j){
							auto fac = 1.0;
							for(auto cl = 0u; cl < ncla; cl++){
								if(cl != foc_cl) fac *= per[cl][co.cla_comp[cl]]/100;
							}
						
							fac *= num;
							auto faci = (unsigned int)(fac);
							
							init_cond.cnum[c] = faci;
							Fraction frac; frac.c = c; frac.fr = fac-faci;
							frac_list.push_back(frac);
							
							n += faci;
						}
					}
					
					if(n > num) emsg_input("Should not be larger");
					if(n < num){
						sort(frac_list.begin(),frac_list.end(),Fraction_ord);
						for(auto k = 0u; k < num-n; k++){
							init_cond.cnum[frac_list[k].c]++; 
						}
					}
				}
			}
		}
		else{
			if(tab.ncol != ncla+1){
				alert_source("Does not have the right number of columns",so);
				return;
			}
			
			init_cond.cnum.resize(comp.size(),UNSET);
			
			for(auto c = 0u; c < comp.size(); c++){
				if(comp[c].erlang_hidden == true) init_cond.cnum[c] = 0;
			}
			
			for(auto r = 0u; r < tab.nrow; r++){
				auto c = 0u;
				for(auto cl = 0u; cl < ncla; cl++){
					auto name = tab.ele[r][cl];
				
					const auto &claa = cla[cl];
				
					auto j = claa.hash_comp.find(name);
					if(j == UNSET){ alert_source("Could not find '"+name+"'",so); return;}
			
					c += comp_mult[cl]*j;
				}
				
				auto ele = tab.ele[r][ncla];
				auto num = number(ele);
				if(num == UNSET){ alert_source("Not a number '"+ele+"'",so); return;}
				init_cond.cnum[c] = num;
			}
			
			for(auto i = 0u; i < comp.size(); i++){
				if(init_cond.cnum[i] == UNSET){
					alert_source("Population for '"+comp[i].name+"' not set",so); return;
				}
			}
		}
		
		if(false){
			cout << "init pop" << endl;
			for(auto c = 0u; c < comp.size(); c++){
				const auto &co = comp[c];
				cout << co.name << " " << init_cond.cnum[c] << " pop" << endl;
			}
		}
		break;
	
	case INIT_POP_DIST:
		if(foc_cl != UNSET){
			const auto &claa = cla[foc_cl];
			
			init_cond.focal_cl = foc_cl;
			init_cond.mult = comp_mult[foc_cl];
			
			for(auto cgl = 0u; cgl < comp_gl.size(); cgl++){
				const auto &co = comp_gl[cgl];
				if(co.cla_comp[foc_cl] == 0){
					init_cond.comp_reduce.push_back(cgl);
				}
			}
			init_cond.N_reduce = init_cond.comp_reduce.size();
			
			init_cond.comp_reduce_ref.resize(comp_gl.size());
			for(auto c = 0u; c < claa.ncomp; c++){
				for(auto cred = 0u; cred < init_cond.N_reduce; cred++){
					auto cgl = init_cond.comp_reduce[cred] + c*init_cond.mult;
					init_cond.comp_reduce_ref[cgl].c = c;
					init_cond.comp_reduce_ref[cgl].cred = cred;
				}
			}
		
			init_cond.comp_prior.resize(claa.ncomp);
			
			vector <bool> prior_set(claa.ncomp,false);
			
			// Sets internal Erlang 
			for(auto c = 0u; c < claa.ncomp; c++){
				const auto &co = claa.comp[c];
				if(co.erlang_hidden == true){
					init_cond.comp_prior[c] = convert_text_to_prior("fix(0)",UNSET,"Erlang",false);
					prior_set[c] = true;
				}
			}
			
			auto &af = init_cond.alpha_focal;
			
			af.resize(ncla);
			for(auto cl = 0u; cl < ncla; cl++){
				if(cl != foc_cl){
					const auto C = cla[cl].ncomp;
					af[cl].resize(C,UNSET);
					for(auto c = 0u; c < C; c++){
						const auto &co = cla[cl].comp[c];
						if(co.erlang_hidden) af[cl][c] = ALPHA_ZERO;
					}
				}
			}
			
			if(tab.ncol != 2) alert_source("Should have two columns",so);
			for(auto r = 0u; r < tab.nrow; r++){
				auto name = tab.ele[r][0];
				auto val = tab.ele[r][1];
							
				auto flag = false;
				for(auto cl = 0u; cl < ncla; cl++){
					for(auto c = 0u; c < cla[cl].ncomp; c++){
						if(cla[cl].comp[c].name == name){
							if(cl == foc_cl){
								auto pri = convert_text_to_prior(val,so.line_num,"Population for "+name,false);
								if(pri.error != "") alert_source(pri.error,so,1,r);
								else{
									init_cond.comp_prior[c] = pri;
									prior_set[c] = true;
								}
							}
							else{
								if(af[cl][c] != UNSET){
									alert_source("The alpha value for compartment '"+name+"' is set more than once",so,0,r);
								}	
								
								af[cl][c] = number(val);
								if(af[cl][c] == UNSET){
									alert_source("The value '"+name+"' must be a number",so,1,r);
								}
							}
							
							flag = true;
							break;
						}
					}
					if(flag == true) break;
				}
				
				if(flag == false){
					alert_source("The compartment '"+name+"' is not recognised",so,0,r);
				}
			}		

			for(auto c = 0u; c < claa.ncomp; c++){
				if(prior_set[c] == false){
					alert_source("The prior for '"+claa.comp[c].name+"' is not set",so);
				}				
			}			
			
			// By default sets alpha values to 1
			for(auto cl = 0u; cl < ncla; cl++){
				if(cl != foc_cl){
					for(auto c = 0u; c < cla[cl].ncomp; c++){
						if(af[cl][c] == UNSET) af[cl][c] = 1;
					}
				}
			}
		}
		else{
			init_cond.pop_prior = so.pop_prior;
			
			if(tab.ncol != ncla+1){
				alert_source("Does not have the right number of columns",so);
				return;
			}
			
			auto C = comp.size();
			init_cond.alpha.resize(C,UNSET);
			
			for(auto c = 0u; c < C; c++){
				const auto &co = comp[c];
				if(co.erlang_hidden) init_cond.alpha[c] = ALPHA_ZERO;
			}
			
			for(auto r = 0u; r < tab.nrow; r++){
				auto c = 0u;
				for(auto cl = 0u; cl < ncla; cl++){
					auto name = tab.ele[r][cl];
				
					const auto &claa = cla[cl];
					auto j = claa.hash_comp.find(name);
					if(j == UNSET){ alert_source("Could not find '"+name+"'",so); return;}
			
					c += comp_mult[cl]*j;
				}
				
				auto ele = tab.ele[r][ncla];
				auto num = number(ele);
				if(num == UNSET){ alert_source("Not a number '"+ele+"'",so); return;}
				init_cond.alpha[c] = num;
			}
			
			for(auto i = 0u; i < comp.size(); i++){
				if(init_cond.alpha[i] == UNSET){
					alert_source("Alpha value for '"+comp[i].name+"' not set",so); return;
				}
			}
		
			if(false){
				cout << init_cond.pop_prior.type << "prior type" << endl;
				cout << "init pop prior" << endl;
				for(auto c = 0u; c < comp.size(); c++){
					const auto &co = comp[c];
					cout << co.name << " " << init_cond.alpha[c] << " alpha" << endl;
				}
			}
		}
		break;
		
	default: emsg_input("de f prob"); break;
	}
}
		

/// Fills in missing percentage when defining initial population
vector < vector <double> > Species::pop_get_percentage(const DataSource &so)
{
	const auto &tab = so.table;
	auto foc_cl = so.focal_cl; 
	
	vector < vector <double> > per;

	per.resize(ncla);
	for(auto cl = 0u; cl < ncla; cl++){
		for(auto c = 0u; c < cla[cl].ncomp; c++){
			double val = UNSET; if(cla[cl].comp[c].erlang_hidden == true) val = 0;
			per[cl].push_back(val);
		}
	}					
	
	if(tab.ncol != 2) alert_source("Should have two columns",so);
	for(auto r = 0u; r < tab.nrow; r++){
		auto name = tab.ele[r][0];
		
		auto flag = false;
		for(auto cl = 0u; cl < ncla; cl++){
			for(auto c = 0u; c < cla[cl].ncomp; c++){
				if(cla[cl].comp[c].name == name){
					if(per[cl][c] != UNSET){
						alert_source("The comartment '"+name+"' is set more than once",so,0,r);
					}					
					
					auto val = tab.ele[r][1];
					if(cl == foc_cl){
						per[cl][c] = number(val);
						if(per[cl][c] == UNSET) alert_source("The value '"+val+"' is not a number1",so,1,r);
					}
					else{
						if(is_percent(val) == false){
							alert_source("The value '"+val+"' is not a percentage",so,1,r);
						}
						else{
							per[cl][c] = number(val.substr(0,val.length()-1));
						}
					}
					
					flag = true;
					break;
				}
			}
			if(flag == true) break;
		}
		
		if(flag == false){
			alert_source("The compartment '"+name+"' is not recognised",so,0,r);
		}
	}		

	for(auto cl = 0u; cl < ncla; cl++){
		const auto &claa = cla[cl];
		auto sum = 0.0;
		vector <unsigned int> list;
		
		for(auto c = 0u; c < claa.ncomp; c++){
			if(per[cl][c] == UNSET) list.push_back(c);
			else sum += per[cl][c];
		}
		
		string st = ""; 
		for(auto j = 0u; j < list.size(); j++){
			if(j != 0) st += ", ";
			st += claa.comp[list[j]].name;
		}
		
		if(cl == foc_cl){
			if(list.size() != 0){
				alert_source("Populations for the compartment(s) '"+st+"' must be set",so);
			}
		}
		else{
			if(list.size() == 0){
				alert_source("Not all percentages for compartmental population percentages in '"+claa.name+"' should be set",so);
			}
			
			if(list.size() > 2){
				alert_source("All but one of the compartmental population percentages in '"+st+"' should be set",so);
			}
			
			if(list.size() == 1){
				if(sum > 100.0){
					alert_source("Population percentages in '"+claa.name+"' add to over 100%",so);
				}
				else{
					per[cl][list[0]] = 100-sum; 
				}
			}
		}
	}	
	
	return per;
}


/// add-pop / remove-pop command
void Species::add_pop_data(const DataSource &so, int sign) 
{	
	const auto &tab = so.table;
	
	if(add_rem_pop_on == false){
		add_rem_pop.resize(T);
		for(auto ti = 0u; ti < T; ti++){
			add_rem_pop[ti].resize(comp_gl.size(),0);
		}
	}
	
	add_rem_pop_on = true;
	
	auto warn = false;
	for(auto j = 0u; j < tab.nrow; j++){
		auto t_str = tab.ele[j][0];
		double t = number(t_str);
		if(t_str == "start") t = details.t_start;
		
		if(t < details.t_start || t >= details.t_end){
			warn = true;
		}
		else{
			auto ti = get_ti(calc_tdiv(t,details)); 
			
			if(ncla+2 != tab.ncol) emsg_input("Columns not right1");
			
			auto name = t_str;
			
			auto c = 0u;
			for(auto cl = 0u; cl < ncla; cl++){
				auto name = tab.ele[j][cl+1];
							
				const auto &claa = cla[cl];
				
				auto j = claa.hash_comp.find(name);
				if(j == UNSET){ alert_source("Could not find '"+name+"'",so); return;}
						
				c += comp_mult[cl]*j;
			}

			auto value_str = tab.ele[j][tab.ncol-1];
			double value = number(value_str);
		
			add_rem_pop[ti][c] += value*sign; 
		}
	}
	
	if(warn) data_ignored(so);
	
	if(false){
		for(auto ti = 0u; ti < T; ti++){
			cout << ti << ": ";
			for(auto va : add_rem_pop[ti]) cout << va << ",";
			cout << " add rem" << endl;
		}
		emsg_input("don");
	}
}


/// Creates warning mnessage if data is outside of time range
void Species::data_ignored(const DataSource &so)
{
	data_warning.push_back("Data outside of time range is ignored in source '"+so.name+"'."); 
}


/// Sets add_rem_pop_change
void Species::set_add_rem_pop_change()
{
	if(add_rem_pop_on == false) return;
	
	add_rem_pop_change.resize(T);
	for(auto ti = 0u; ti < T; ti++){
		for(auto c = 0u; c < comp_gl.size(); c++){
			if(add_rem_pop[ti][c] != 0){
				add_rem_pop_change[ti].push_back(c);
			}
		}
	}
}


/// add-ind command
void Species::add_ind_data(const DataSource &so)
{
	const auto &tab = so.table;

	for(auto j = 0u; j < tab.nrow; j++){
		auto i = find_individual(tab.ele[j][0]);
		auto &ind = individual[i];

		auto t_str = tab.ele[j][1];
		double t = number(t_str);
		if(t_str == "start") t = details.t_start;
		
		if(t < details.t_start || t >= details.t_end){
			alert_source("The time '"+t_str+"' must be between the start and end times",so,1,j); 
			return;
		}
	
		if(ncla+2 != tab.ncol) emsg_input("Columns not right2");
		
		auto name = t_str;
		vector <string> val;
		for(auto cl = 0u; cl < ncla; cl++){
			const auto &el = tab.ele[j][cl+2]; 
			name += "&"+el;
			val.push_back(el);
		}

		auto vec = hash_enter.get_vec_string(name);
		auto k = hash_enter.existing(vec);
		if(k == UNSET){
			k = enter.size();
			hash_enter.add(k,vec);
	
			// Determines if initial compartment is exactly specified
			auto c_set = 0u;
			for(auto cl = 0u; cl < ncla; cl++){
				const auto &claa = cla[cl];

				auto c = claa.hash_comp.find(val[cl]);
				if(c == UNSET){ c_set = UNSET; break;}
				c_set += comp_mult[cl]*c;
			}
		
			Enter ep;
			ep.name = name;
			ep.tdiv = calc_tdiv(t,details); 
			ep.c_set = c_set;
			
			if(c_set == UNSET){
				ep.cla.resize(ncla);
				for(auto cl = 0u; cl < ncla; cl++){
					const auto &claa = cla[cl];

					auto &ent_cl = ep.cla[cl];
					
					auto c = claa.hash_comp.find(val[cl]);
					
					ent_cl.c_set = c;
					if(c == UNSET){
						string emg;
						auto comp_prob_str = find_comp_prob_str(cl,val[cl],LOWER_BOUND,emg);
						if(emg != ""){ alert_source(emg,so,cl,j); return;}
						
						ent_cl.eqn = create_eqn_vector(comp_prob_str,COMP_PROB,so);
					}
				}
			}
			ep.set = true;
			
			enter.push_back(ep);
		}
	
		ind.enter_ref = k;
	
		EventData ev; 
		ev.so = so.index;
		ev.type = ENTER_EV;
		ev.move_c = UNSET;
		ev.cl = UNSET;
		ev.tr = UNSET;
		ev.tdiv = calc_tdiv(t,details);
		ind.ev.push_back(ev);
	}		
}


/// remove-ind command
void Species::remove_ind_data(const DataSource &so)
{
	const auto &tab = so.table;
	
	auto warn = false;
	for(auto j = 0u; j < tab.nrow; j++){
		auto i = find_individual(tab.ele[j][0]);
		auto &ind = individual[i];

		auto t_str = tab.ele[j][1];
		double t = number(t_str);
		if(t_str == "start") t = details.t_start;
		
		if(t < details.t_start || t >= details.t_end){
			warn = true;
		}
		else{
			if(tab.ncol != 2) emsg_input("Columns not right3");
		
			EventData ev; 
			ev.so = so.index;
			ev.type = LEAVE_EV;
			ev.cl = UNSET;
			ev.tr = UNSET;
			ev.tdiv = calc_tdiv(t,details);
			ind.ev.push_back(ev);
		}
	}		
	
	if(warn) data_ignored(so);
}


/// move-ind command
void Species::move_ind_data(const DataSource &so)
{
	auto cl = so.cl;
		
	const auto &claa = cla[cl];
	for(const auto &tr : claa.tra){
		switch(tr.type){
		case EXP_MEAN_NM: case GAMMA: case ERLANG: case LOG_NORMAL: case WEIBULL: case PERIOD:
			alert_source("Move data cannot be applied to classification '"+claa.name+"' because this contains non-Markovian transitions.",so);
			break;
		default: break;
		}
	}
	
	const auto &tab = so.table;
	auto warn = false;
	for(auto j = 0u; j < tab.nrow; j++){
		auto i = find_individual(tab.ele[j][0]);
		auto &ind = individual[i];

		auto t_str = tab.ele[j][1];
		double t = number(t_str);
		if(t_str == "start") t = details.t_start;
		
		if(t <= details.t_start || t >= details.t_end){
			warn = true;
		}
		else{
			if(tab.ncol != 3) emsg_input("Columns not right4");
			
			auto val = tab.ele[j][2];
			auto c = find_c(cl,val);
			if(c == UNSET) alert_source("Value '"+val+"' is not a compartment2",so,2,j);
				
			EventData ev; 
			ev.so = so.index;
			ev.type = MOVE_EV;
			ev.move_c = c;
			ev.cl = cl;
			ev.tr = UNSET;
			ev.tdiv = calc_tdiv(t,details);
			
			// Checks that 
			for(const auto &eve : ind.ev){
				if(eve.type == MOVE_EV && eve.tdiv == ev.tdiv){
					alert_source("Individual '"+ind.name+"' cannot move more than once at time '"+t_str+"'",so,1,j); 
				}
			}
			
			ind.ev.push_back(ev);			
		}
	}		
	
	if(warn) data_ignored(so);
}


/// comp-data command
void Species::comp_data(const DataSource &so)
{
	const auto &tab = so.table;
	
	auto warn = false;
	for(auto j = 0u; j < tab.nrow; j++){
		auto i = find_individual(tab.ele[j][0]);
		auto &ind = individual[i];

		auto val = tab.ele[j][1];
		auto t = number(val);
		if(val == "start") t = details.t_start;
		if(t == UNSET){
			alert_source("Value '"+val+"' is not a number2",so,1,j);
			return;
		}
		
		if(t < details.t_start || t > details.t_end){
			warn = true;
		}
		else{
			auto cl = so.cl;
			
			auto comp = tab.ele[j][2];
			if(comp != missing_str){
				auto c_exact = find_c(cl,comp);

				ObsData ob; 
				ob.so = so.index;
				ob.type = OBS_COMP_EV;
				ob.c_exact = c_exact;
				ob.not_alive = false;
				if(c_exact == UNSET){
					if(comp == not_alive_str) ob.not_alive = true;
					
					string emsg;
					auto prob_str = find_comp_prob_str(cl,comp,LOWER_BOUND,emsg);
					if(emsg != ""){ alert_source(emsg,so,2,j); return;}
					ob.c_obs_prob_eqn = create_eqn_vector(prob_str,COMP_PROB,so);
				}
				ob.cl = so.cl;
				ob.tdiv = calc_tdiv(t,details);
				ob.time_vari = false;
				ind.obs.push_back(ob);		
			}
		}
	}

	if(warn) data_ignored(so);
}


/// Genetic data command (this makes sure individual is infected when an observation is made
void Species::genetic_data(const DataSource &so)
{
	if(trans_tree == false) return;
	
	const auto &tab = so.table;
	
	const auto &claa = cla[infection_cl];
	string str = "";
	for(auto c = 0u; c < claa.ncomp; c++){
		if(claa.comp[c].infected == COMP_INFECTED){
			if(str != "") str += "|";
			str += claa.comp[c].name;
		}
	}
	
	string ems;
	auto prob_str = find_comp_prob_str(infection_cl,str,LOWER_BOUND,ems);
	if(ems != "") emsg_input("Should not be error");
		
	auto p_eqn = create_eqn_vector(prob_str,COMP_PROB,so);
	
	auto warn = false;
	for(auto j = 0u; j < tab.nrow; j++){
		auto i = find_individual(tab.ele[j][0],false);
		if(i != UNSET){
			auto &ind = individual[i];

			auto val = tab.ele[j][1];
			auto t = number(val);
			if(t == UNSET){
				alert_source("Value '"+val+"' is not a number2",so,1,j);
				return;
			}
			
			if(t < details.t_start || t > details.t_end){
				warn = true;
			}
			else{
				ObsData ob; 
				ob.so = so.index;
				ob.type = OBS_COMP_EV;
				ob.c_exact = UNSET;
				ob.not_alive = false;
				ob.c_obs_prob_eqn = p_eqn;
				ob.cl = infection_cl;
				ob.tdiv = calc_tdiv(t,details);
				ob.time_vari = false;
				ind.obs.push_back(ob);	
			}				
		}
	}
	
	if(warn) data_ignored(so);
}


/// test-data command
void Species::test_data(const DataSource &so)
{
	const auto &tab = so.table;
	
	auto warn = false;
	for(auto j = 0u; j < tab.nrow; j++){
		auto i = find_individual(tab.ele[j][0]);
		auto &ind = individual[i];

		auto val = tab.ele[j][1];
		auto t = number(val);
		if(t == UNSET){
			alert_source("Value '"+val+"' is not a number3",so,1,j);
			return;
		}
		
		if(t < details.t_start || t > details.t_end){
			warn = true;
		}
		else{
			ObsData ob; 
			ob.type = OBS_TEST_EV;
			ob.so = so.index;
			ob.cl = so.cl;
			ob.tdiv = calc_tdiv(t,details);
			ob.not_alive = false;
			ob.time_vari = false;
			ob.Se_eqn = he(add_equation_info(so.obs_model.Se_str,SE),so);
			ob.Sp_eqn = he(add_equation_info(so.obs_model.Sp_str,SP),so);
			
			auto res = tab.ele[j][2];
			if(res == so.obs_model.diag_pos) ob.test_res = true;
			else{
				if(res == so.obs_model.diag_neg) ob.test_res = false;
				else{
					alert_source("Value '"+res+"' is not a positive or negative test result",so,2,j);
					return;
				}
			}
		
			ind.obs.push_back(ob);
		}
	}		
	
	if(warn) data_ignored(so);
}


/// pop-data command
void Species::population_data(const DataSource &so, Hash &hash_pop_filter)
{
	const auto &tab = so.table;
	
	auto cf = set_comp_filt(so.filter_str,UNSET,LOWER_BOUND,so);

	auto warn = false;
	for(auto j = 0u; j < tab.nrow; j++){
		auto time_str = tab.ele[j][0];
		auto t = number(time_str);
		if(t == UNSET){
			alert_source("The time '"+time_str+"' is not a number4",so,0,j);
			return;
		}
		
		if(t < details.t_start || t > details.t_end){
			warn = true;
		}
		else{
			auto col = 1;
			auto name = so.filter_str;
			for(auto cl = 0u; cl < ncla; cl++){
				if(cf.cla[cl].type == FILE_FILT){
					name += "|file:"+tab.ele[j][col];
					col++;
				}
			}	
			
			auto pf = hash_pop_filter.find(name);
			if(pf == UNSET){
				pf = pop_filter.size();
				hash_pop_filter.add(pf,name);
				
				auto co = 1;
				auto cf2 = cf;
				for(auto cl = 0u; cl < ncla; cl++){
					if(cf.cla[cl].type == FILE_FILT){
						cf2.cla[cl].type = COMP_FILT;
						string emsg;
						cf2.cla[cl].comp_prob_str = find_comp_prob_str(cl,tab.ele[j][co],LOWER_BOUND,emsg);
						if(emsg != ""){ alert_source(emsg,so,co,j); return;}
						co++;
					}
				}	
				
				PopFilter pofi; 
				pofi.name = name;
				pofi.comp_prob_eqn = create_eqn_vector(global_convert(cf2),COMP_PROB,so);
				pop_filter.push_back(pofi);
			}
			
			auto val_str = tab.ele[j][col];
			auto value = number(val_str);
			if(value == UNSET){
				alert_source("The value '"+val_str+"' is not a number5",so,1,j);
				return;
			}
			col++;
			
			PopData pd;
			pd.so = so.index;
			pd.tdiv = calc_tdiv(t,details);
			pd.ref = pf;
			pd.type = set_obs_mod_type(so.obs_model);
			pd.value = value;
			pd.obs_mod_val = set_obs_mod_val(value,j,col,so.obs_model,tab,so);
			pd.time_vari = false;
			pop_data.push_back(pd);
		}
	}		
	
	if(warn) data_ignored(so);
}


/// Converts from obsmod type to variety 
ObsModelVariety Species::set_obs_mod_type(const ObsModel &om) const 
{
	switch(om.type){
	case NORMAL_PERCENT_OBSMOD: case NORMAL_SD_OBSMOD: case NORMAL_FILE_OBSMOD: return NORMAL_OBS;
	case POISSON_OBSMOD: return POISSON_OBS;
	case NEGBIN_OBSMOD: case NEGBIN_FILE_OBSMOD: return NEGBIN_OBS;
	}
	emsg_input("Option not recognised");
	return NORMAL_OBS;
}

 
/// Gets the standard deviation (either from observation model or data table)
double Species::set_obs_mod_val(double value, unsigned int j, unsigned int col, const ObsModel &om, const Table &tab, const DataSource &so)
{
	switch(om.type){
	case NORMAL_PERCENT_OBSMOD: 
		{
			auto sd = value*om.percent/100; 
			if(sd < 1) sd = 1;
			return sd;
		}
		break;
		
	case NORMAL_SD_OBSMOD: 
		{
			auto sd = om.sd;
			if(sd < 1) sd = 1;
			return sd;
		}
		break;
		
	case NORMAL_FILE_OBSMOD: 
		{
			auto sd_str = tab.ele[j][col];
			auto sd = number(sd_str);
			if(sd == UNSET){
				alert_source("The standard deviation '"+sd_str+"' is not a number6",so,2,j);
				return sd;	
			}
			if(sd < 1) sd = 1;
			return sd;
		}
		break;
		
	case POISSON_OBSMOD:
		return UNSET;
		
	case NEGBIN_OBSMOD:
		return om.p;
		
	case NEGBIN_FILE_OBSMOD:
		{
			auto p_str = tab.ele[j][col];
			auto p = number(p_str);
			if(p == UNSET){
				alert_source("The standard deviation '"+p_str+"' is not a number6",so,2,j);
				return p;	
			}
			return p;
		}
		break;
	}
	
	return UNSET;
}

 
/// trans_data  command
void Species::trans_data(const DataSource &so)
{
	const auto &tab = so.table;
	
	// Accounts for the time period measurements are made
	auto so_ts = details.t_start;
	auto so_te = details.t_end;
	if(so.time_range == SPEC_TIME){
		so_ts = so.time_start;
		so_te = so.time_end;
	}
	
	if(so_ts > details.t_end || so_te < details.t_start){
		data_warning.push_back("Data entirely outside of time range ignored in source '"+so.name+"'."); 
		return;
	}

	// Adds the times of events
	auto warn = false;
	for(auto j = 0u; j < tab.nrow; j++){
		auto i = find_individual(tab.ele[j][0]);
		auto &ind = individual[i];

		auto val = tab.ele[j][1];
		auto t = number(val);
		
		if(t == UNSET){
			alert_source("Value '"+val+"' is not a number",so,1,j);
			return;
		}
		
		if(so.time_range == SPEC_TIME && (t < so_ts || t > so_te)){
			alert_source("Value '"+val+"' is outside the time range over which transitions are observed.",so,1,j);
		}			
		
		if(t < details.t_start || t >= details.t_end){
			warn = true;
		}
		else{
			ObsData ob; 
			ob.so = so.index;
			ob.ref = obs_trans.size();
			ob.type = OBS_TRANS_EV;
			ob.cl = so.cl;
			ob.tdiv = calc_tdiv(t,details);
			ob.time_vari = false;
			ob.not_alive = false;
			
			ind.obs.push_back(ob);
		}
	}		

	auto cf = set_comp_filt(so.filter_str,UNSET,LOWER_UPPER_BOUND,so);

	string errmsg;
	auto trans_filt = set_trans_filt(so.cl,so.filter_trans_str,LOWER_UPPER_BOUND,errmsg);
	
	if(errmsg != ""){ alert_source(errmsg,so); return;}
	
	if(so_ts < details.t_start) so_ts = details.t_start;
	if(so_te > details.t_end) so_ts = details.t_end;
	
	auto ti_min = get_ti(calc_tdiv(so_ts,details));
	auto ti_max = get_ti(calc_tdiv(so_te,details));
		
	auto prob_str = trans_global_convert(so.cl,trans_filt,cf);
	ObsTrans ob_tr;
	ob_tr.name = so.filter_trans_str;
	ob_tr.type = OBS_TRANS_EV; 
	ob_tr.tra_prob_eqn = create_eqn_vector(prob_str,TRANS_PROB,so);
	ob_tr.ti_min = ti_min;
	ob_tr.ti_max = ti_max;

	obs_trans.push_back(ob_tr);
	
	obs_trans_exist = true;
	
	if(warn) data_ignored(so);
}


/// Sets if observation is OBS_SOURCE_EV, OBS_TRANS_EV or OBS_SINK_EV
void Species::set_ob_trans_ev(const vector <Equation> &eqn)
{
	for(auto &ot : obs_trans){
		auto source_fl = false, trans_fl = false, sink_fl = false;

		for(auto tr = 0u; tr < tra_gl.size(); tr++){
			const auto &eq = eqn[ot.tra_prob_eqn[tr].eq_ref];
			if(!eq.is_zero()){
				switch(tra_gl[tr].variety){
				case NORMAL: trans_fl = true; break;
				case SOURCE_TRANS: source_fl = true; break; 
				case SINK_TRANS: sink_fl = true; break;
				}
			}
		}

		if(source_fl && trans_fl){
			alert_input("Transition data '"+ot.name+"' cannot include source and non-source transitions");
		}

		if(sink_fl && trans_fl){
			alert_input("Transition data '"+ot.name+"' cannot include sink and non-sink transitions");
		}

		if(source_fl && sink_fl){
			alert_input("Transition data '"+ot.name+"' cannot include source and sink transitions");
		}

		ot.type = OBS_TRANS_EV;
		if(source_fl) ot.type = OBS_SOURCE_EV;
		if(sink_fl) ot.type = OBS_SINK_EV;
	}

	for(auto &ind : individual){
		for(auto &ob : ind.obs){
			if(ob.type == OBS_TRANS_EV) ob.type = obs_trans[ob.ref].type;
		}
	}
}


/// trans_data and pop-trans-data command
void Species::popu_trans_data(const DataSource &so, Hash &hash_pop_trans_filter)
{
	const auto &tab = so.table;
	
	auto cf = set_comp_filt(so.filter_str,UNSET,LOWER_BOUND,so);
	
	string emg;
	auto trans_filt = set_trans_filt(so.cl,so.filter_trans_str,LOWER_BOUND,emg);
	if(emg != ""){ alert_source(emg,so); return;}
	
	auto warn = false;
	for(auto j = 0u; j < tab.nrow; j++){
		auto tstart_str = tab.ele[j][0];
		auto tstart = number(tstart_str);
		if(tstart == UNSET){
			alert_source("The start time '"+tstart_str+"' is not a number7",so,0,j);
			return;
		}
		
		auto tend_str = tab.ele[j][1];
		auto tend = number(tend_str);
		if(tend == UNSET){
			alert_source("The end time '"+tend_str+"' is not a number7",so,0,j);
			return;
		}
		
		if(tstart < details.t_start || tend > details.t_end){
			warn = true;
		}
		else{
			auto col = 2;
			auto name = so.filter_trans_str+"|"+so.filter_str;
			for(auto cl = 0u; cl < ncla; cl++){
				if(cl != so.cl && cf.cla[cl].type == FILE_FILT){	
					name += "|file:"+tab.ele[j][col];
					col++;
				}
			}	
			
			auto pf = hash_pop_trans_filter.find(name);
			
			if(pf == UNSET){
				pf = pop_trans_filter.size();
				hash_pop_trans_filter.add(pf,name);
				
				auto co = 2;
				auto cf2 = cf;
				for(auto cl = 0u; cl < ncla; cl++){
					if(cl != so.cl && cf.cla[cl].type == FILE_FILT){
						cf2.cla[cl].type = COMP_FILT;
						string emsg;
						cf2.cla[cl].comp_prob_str = find_comp_prob_str(cl,tab.ele[j][co],LOWER_BOUND,emsg);
						if(emsg != ""){ alert_source(emsg,so,col,j); return;}
						co++;
					}
				}
			
				auto prob_str = trans_global_convert(so.cl,trans_filt,cf2);
		
				PopTransFilter pofi; 
				pofi.name = name;
				pofi.trans_prob_eqn = create_eqn_vector(prob_str,TRANS_PROB,so);
				pop_trans_filter.push_back(pofi);
			}
			
			auto val_str = tab.ele[j][col];
			auto value = number(val_str);

			if(value == UNSET){
				alert_source("The value '"+val_str+"' is not a number9",so,1,j);
				return;
			}
			
			PopTransData ptd;
			ptd.so = so.index;
			ptd.tdivmin = calc_tdiv(tstart,details);
			ptd.tdivmax = calc_tdiv(tend,details);
			ptd.ref = pf;		
			ptd.type = set_obs_mod_type(so.obs_model);
			ptd.value = value;
			ptd.obs_mod_val = set_obs_mod_val(value,j,col,so.obs_model,tab,so);
			ptd.time_vari = false;
			pop_trans_data.push_back(ptd);
			
			pop_trans_data_exist = true;
		}
	}		
	
	if(warn) data_ignored(so);
}


/// Finds the compartment from the name
unsigned int Species::find_c(unsigned int cl, string name) const
{
	//const auto &claa = cla[cl];
	//auto vec = claa.hash_comp.get_vec_string(name);
	//return claa.hash_comp.existing(vec);
	return cla[cl].hash_comp.find(name);
}


/// Finds an individual with a given name, overwise creates a new individual
unsigned int Species::find_individual(string name, bool create)
{
	auto vec = hash_ind.get_vec_string(name);
	auto i = hash_ind.existing(vec);
	
	if(i == UNSET && create){
		i = individual.size();
		if(i >= details.individual_max){
			alert_input("The number of individuals exceeds the limit of "+to_string(details.individual_max)+".");
		}
		
		hash_ind.add(i,vec);
		
		IndData ind;
		ind.name = name;
		ind.enter_ref = UNSET;
		individual.push_back(ind);
	}
	
	return i;
}


/// Alerts a problem with a data source
void Species::alert_source(string st, const DataSource &so, unsigned int c, unsigned int r)
{
	if(c != UNSET) st += " (col '"+so.table.heading[c]+"', row "+tstr(r+1)+")";
	WarnData wa; wa.te = st; wa.line_num = so.line_num;
	warn.push_back(wa);
}


/// Orders X_vector such that order agrees with individual
void Species::X_vector_order()
{
	for(auto &fe : fix_effect){
		const auto &xvec_old = fe.X_vector;
		
		Xvector xvec_new;
		for(auto i = 0u; i < individual.size(); i++){
			auto name = individual[i].name;
			xvec_new.ind_list.push_back(name);
			
			auto j = xvec_old.hash_ind_list.find(name);
			if(j == UNSET) xvec_new.value.push_back(UNSET);
			else xvec_new.value.push_back(xvec_old.value[j]);			
		}
		xvec_new.hash_ind_list.create(xvec_new.ind_list);

		fe.X_vector = xvec_new;		
	}
}


/// Converts from a filter in different classification to one for each global compartment
/// This consists of multiplying together a series of equations 
vector <string> Species::global_convert(const Filter &filt) const
{
	const auto N = comp_gl.size();
	vector <string> gfilt(N);
	
	for(auto c = 0u; c < N; c++){
		const auto &cgl = comp_gl[c];
		
		string te="";
		for(auto cl = 0u; cl < ncla; cl++){
			auto &fcl = filt.cla[cl];
			if(fcl.type != COMP_FILT) emsg_input("Should be COMP_FILT");
			if(te != "") te += "*";
		
			const auto cc = cgl.cla_comp[cl];
			te += "("+fcl.comp_prob_str[cc]+")";
		}
	
		gfilt[c] = te;
	}
	
	return gfilt;
}


/// Converts from transition filter in classification to global filter
// For each transition gives a list of equations which multiply to give fraction observed
vector <string> Species::trans_global_convert(unsigned int cl, const vector <string> &trans_str, const Filter &comp_filt)
{
	vector <string> gfilt;
	auto filtg = global_convert(comp_filt); 
	for(auto tr = 0u; tr < tra_gl.size(); tr++){
		const auto &trg = tra_gl[tr];
		if(trg.cl != cl) gfilt.push_back("0");
		else{
			auto c = trg.i; if(c == UNSET) c = trg.f;
			auto str = "("+trans_str[trg.tr]+")*("+filtg[c]+")";
			gfilt.push_back(str);
		}
	}	
	
	return gfilt;
}


/// Used to order events
bool EventData_ord (const EventData &ev1, const EventData &ev2)                      
{ return (ev1.tdiv < ev2.tdiv); };  


/// Used to order observations
bool ObsData_ord (const ObsData &ob1, const ObsData &ob2)                      
{ return (ob1.tdiv < ob2.tdiv); };  


/// Used to order population data times
bool PopData_ord (const PopData &pd1, const PopData &pd2)                      
{ return (pd1.tdiv < pd2.tdiv); }; 


/// Used to order poptrans data times
bool PopTransData_ord (const PopTransData &pd1, const PopTransData &pd2)                      
{ return (pd1.tdivmin < pd2.tdivmin); }; 


/// Orders individual events by time
void Species::order_data_events()
{
	sort(pop_data.begin(),pop_data.end(),PopData_ord);    
	sort(pop_trans_data.begin(),pop_trans_data.end(),PopTransData_ord);    
	
	// Sorts events by time
	for(auto i = 0u; i < individual.size(); i++){
		auto &ind = individual[i];
		sort(ind.ev.begin(),ind.ev.end(),EventData_ord);    
		sort(ind.obs.begin(),ind.obs.end(),ObsData_ord);    
	}
		
	if(false){
		for(const auto &pd : pop_data){
			cout << pd.tdiv << " " << pd.value << " " << pd.obs_mod_val << " Population data" << endl;
		}
		
		for(const auto &ptd : pop_trans_data){
			cout << ptd.tdivmin << " " << ptd.tdivmax << " " << ptd.value << " " << ptd.obs_mod_val << " Population trans data" << endl;
		}
	}
}


/// Given a string (e.g. "Sex=M,Loc=file") this constructs the filter
Filter Species::set_comp_filt(string te, unsigned int cl_not_allow, BoundType bound, const DataSource &so)
{
	Filter filt;
	
	filt.cla.resize(ncla);

	// Sets default values for filter
	for(auto cl = 0u; cl < ncla; cl++){
		const auto &claa = cla[cl];
				
		filt.cla[cl].type = COMP_FILT;
		filt.cla[cl].comp_prob_str.resize(claa.ncomp);
		for(auto c = 0u; c < claa.ncomp; c++){
			filt.cla[cl].comp_prob_str[c] = "1";
		}
	}
	
	if(te != ""){
		auto spl = split(te,',');
		for(auto j = 0u; j < spl.size(); j++){
			auto spl2 = split(spl[j],'=');
			if(spl2.size() != 2){ 
				alert_source("In 'filter' the expression '"+spl[j]+"' must contain an equals sign",so); 
				return filt;
			}
			
			if(spl2[0] == ""){
				alert_source("In 'filter' the value '"+spl[j]+"' does not specify a classification",so);
				return filt;
			}
			
			auto cl = find_cl(spl2[0]);
			if(cl == UNSET){ 
				alert_source("In 'filter' the value '"+spl2[0]+"' is not a classification",so); 
				return filt;
			}
			
			if(cl == cl_not_allow){ 
				alert_source("In 'filter' the value '"+spl2[0]+"' cannot be the same as the transition classification",so); 
				return filt;
			}
						
			if(toLower(spl2[1]) == "file"){
				filt.cla[cl].type = FILE_FILT;
			}
			else{
				filt.cla[cl].type = COMP_FILT;
				string emsg;
				filt.cla[cl].comp_prob_str = find_comp_prob_str(cl,spl2[1],bound,emsg);
				if(emsg != ""){
					alert_source(emsg,so); 
				}
			}
		}
	}
	
	if(false){
		cout << te << " Filter text" << endl;
		for(auto cl = 0u; cl < ncla; cl++){
			const auto &claa = cla[cl];
			for(auto c = 0u; c < claa.ncomp; c++){
				cout << claa.comp[c].name <<" " << filt.cla[cl].comp_prob_str[c] << " comp prob string" << endl;
			}
		}
	}
	
	return filt;
}


/// Given a string, e.g. S|E or M:0.5|F:0.7, works out equations for probability in compartments
vector <string> Species::find_comp_prob_str(unsigned int cl, string te, BoundType bound, string &emsg) const
{
	const auto &claa = cla[cl];
	
	emsg = "";
	
	vector <string> fil(claa.ncomp);
	
	for(auto c = 0u; c < fil.size(); c++){
		fil[c] = "0";
	}	
	
	if(te == not_alive_str){
		return fil;
	}
	
	if(te == "."){
		te = "";
		const auto &claa = cla[cl];
		for(auto c = 0u; c < claa.comp.size(); c++){
			if(te != "") te += "|";
			te += claa.comp[c].name;
		}
	}
	
	auto spl = split_with_bracket(te,'|');
	for(auto k = 0u; k < spl.size(); k++){
		auto spl2 = split(spl[k],':');
		auto c = find_c(cl,spl2[0]);
		if(c == UNSET){
			emsg = "In 'filter' the compartment '"+spl2[0]+"' is not in classification '"+claa.name+"'";
			return fil;
		}
				
		switch(spl2.size()){
		case 1:
			fil[c] = "1";
			break;
			
		case 2:
			fil[c] = add_bound(spl2[1],bound);
			break;
			
		default:
			emsg = "In 'filter' the expression '"+spl[k]+"' is not understood";
			break;
		}
	}
	
	for(auto c = 0u; c < claa.ncomp; c++){
		const auto &co = claa.comp[c];
		if(co.erlang_hidden){
			//auto vec = claa.hash_comp.get_vec_string(co.erlang_source);
			//auto cc = claa.hash_comp.existing(vec);
			auto cc = claa.hash_comp.find(co.erlang_source);
			if(cc == UNSET) { emsg = "Erlang problem"; return fil;}
			fil[c] = fil[cc];
		}
	}

	return fil;
}


/// Given a string, e.g. S->E or S->E:0.5|S->I:0.7, works out equations for prob in transitions
vector <string> Species::set_trans_filt(unsigned int cl, string te, BoundType bound, string &errmsg) const
{
	const auto &claa = cla[cl];
	
	vector <string> fil(claa.ntra);
	
	for(auto tr = 0u; tr < fil.size(); tr++){
		fil[tr] = "0";
	}	
	
	auto spl = split_with_bracket(te,'|');
	for(auto k = 0u; k < spl.size(); k++){
		auto spl2 = split(spl[k],':');
		auto tr_name = replace(spl2[0],"->","→");
		auto tr = find_tr(cl,tr_name);
		if(tr == UNSET){
			errmsg = "The expression '"+spl2[0]+"' is not a transition in '"+claa.name+"'";
			return fil;
		}
				
		switch(spl2.size()){
		case 1:
			fil[tr] = "1";
			break;
			
		case 2:
			fil[tr] = add_bound(spl2[1],bound);
			break;
			
		default:
			errmsg = "The expression '"+spl[k]+"' is not understood";
			break;
		}
	}
		
	return fil;
}

			
/// Adds a boundary to a numeric value 
string Species::add_bound(string te, BoundType bound) const 
{
	stringstream ss;
	switch(bound){
	case LOWER_BOUND: ss << "max(" << te << "|" << LOW_BOUND << ")"; break;
	case LOWER_UPPER_BOUND: ss << "min(max(" << te << "|" << LOW_BOUND << ")|" << UP_BOUND<< ")"; break;
	default: emsg_input("Bound not reconginsed"); break;
	}
	return ss.str();
}

			
/// Creates a vector of euqations from a vector of strings
vector <EquationInfo> Species::create_eqn_vector(const vector <string> &vec, EqnType type, const DataSource &so)
{	
	auto N = vec.size();
	vector <EquationInfo> eqn(N);
	for(auto c = 0u; c < N; c++){
		eqn[c] = he(add_equation_info(vec[c],type),so);
	}
	
	return eqn;
}


/// For each time div and global transition references any transition data associated with it
void Species::init_pop_trans_ref()
{
	T = details.T;
	
	pop_trans_ref.resize(T);
	for(auto ti = 0u; ti < T; ti++){
		pop_trans_ref[ti].resize(tra_gl.size());
	}
	
	for(auto i = 0u; i < pop_trans_data.size(); i++){
		const auto &ptd = pop_trans_data[i];
		const auto &ptf = pop_trans_filter[ptd.ref];
		
		auto ti_min = get_ti(ptd.tdivmin);
		auto ti_max = get_ti(ptd.tdivmax);
		
		for(auto tr : ptf.tr_nonzero){
			for(auto ti = ti_min; ti < ti_max; ti++){
				pop_trans_ref[ti][tr].push_back(i);
			}
		}
	}
}


/// For each time div and global transition references any population data associated with it
void Species::init_pop_data_ref()
{
	T = details.T;
	
	pop_data_ref.resize(T);
	for(auto ti = 0u; ti < T; ti++){
		pop_data_ref[ti].resize(comp_gl.size());
	}
	
	for(auto i = 0u; i < pop_data.size(); i++){
		const auto &ptd = pop_data[i];
		const auto &pf = pop_filter[ptd.ref];

		auto ti = get_ti(ptd.tdiv); if(ti == T) ti--;
		
		for(auto c : pf.c_nonzero){
			pop_data_ref[ti][c].push_back(i);
		}
	}
}


/// Returns the classification number from its name
unsigned int Species::find_cl(string name) const 
{
	for(auto cl = 0u; cl < ncla; cl++){
		if(toLower(cla[cl].name) == toLower(name)) return cl;
	}
	return UNSET;
}


/// Returns the transition number from its name
unsigned int Species::find_tr(unsigned int cl, string name) const 
{
	const auto &tra = cla[cl].tra;
	for(auto tr = 0u; tr < tra.size(); tr++){
		if(tra[tr].name == name) return tr;
	}
	return UNSET;
}


/// Adds a vector of eqn into onto obs_eqn
vector <unsigned int> Species::obs_eqn_add_vec(vector <EquationInfo> &eqn_info) 
{
	vector <unsigned int> vec;
	for(auto i = 0u; i < eqn_info.size(); i++){
		vec.push_back(add_to_vec(obs_eqn,eqn_info[i].eq_ref,hash_obs_eqn));
	}	
	return vec;
}


/// For individuals not added to the system set an uninformative prior on initial state
void Species::set_default_enter()
{
	auto fl = false;
	for(auto &ind : individual){
		if(ind.enter_ref == UNSET){ ind.enter_ref = enter.size(); fl = true;}
	}
	
	if(fl == true){
		Enter ep;
		ep.name = "Uninformative"; 
		ep.tdiv = 0; 
		ep.c_set = UNSET;
		ep.cla.resize(ncla);
		for(auto cl = 0u; cl < ncla; cl++){
			auto &ep_cl = ep.cla[cl];
			ep_cl.c_set = UNSET;
			
			auto N = cla[cl].ncomp;
			auto val = 1.0/N;
			auto ei = add_equation_info(to_string(val),CONST_EQN,p_species,UNSET);

			for(auto c = 0u; c < N; c++) ep_cl.eqn.push_back(ei);
		}
		ep.set = false;
		
		enter.push_back(ep);
	}
}


/// Handles any error message
EquationInfo Species::he(EquationInfo eqn_inf, const DataSource &so)
{
	eqn_inf.line_num = so.line_num;
	
	if(eqn_inf.error) alert_source(eqn_inf.emsg,so);	
	return eqn_inf;
}	


/// Generates a potential list of events which can be added and removed from the system
vector <AddRemLocal> Species::find_add_rem_list() const
{
	vector <AddRemLocal> add_rem_local;
	
	// Works out which global compartments contribute to a population 
	vector <bool> cp(comp_gl.size(),false);
	for(auto c = 0u; c < comp_gl.size(); c++){
		if(comp_gl[c].pop_ref.size() > 0) cp[c] = true;
	}
	
	for(auto c = 0u; c < comp_gl.size(); c++){
		if(cp[c] == false){
			for(auto tr : comp_gl[c].tr_leave){
				vector <unsigned int> tr_list;
				tr_list.push_back(tr);
				
				add_tr_list(tr_list,cp,add_rem_local);
			}
		}
	}
	
	if(false){
		for(const auto &arl : add_rem_local){
			for(auto tr :arl.tr_list) cout << tra_gl[tr].name << ","; 
			cout << "ADD REM" << endl; 
		}
	}
	
	return add_rem_local;
}


/// Add new items to add remove local
void Species::add_tr_list(const vector <unsigned int> &tr_list, const vector <bool> &cp, vector <AddRemLocal> &add_rem_local) const 
{
	auto tr = tr_list[tr_list.size()-1];
	auto cf = tra_gl[tr].f;
	
	if(cf == UNSET || cp[cf] == false){
		AddRemLocal arl; arl.tr_list = tr_list;
		add_rem_local.push_back(arl);
	}
	else{
		for(auto tr : comp_gl[cf].tr_leave){
			auto tr_list_new = tr_list;
			tr_list_new.push_back(tr);		
			add_tr_list(tr_list_new,cp,add_rem_local);			
		}
	}
}


/// Adds any individuals which are unobserved 
void Species::add_unobs_Amatrix_ind()
{
	for(const auto &ieg : ind_eff_group){
		const auto &A = ieg.A_matrix;
		
		if(A.set){
			for(auto na : A.ind_list){
				find_individual(na);
			}
		}
	}
}


/// Makes tiny shifts in time ensure that transition times are not equal
void Species::jiggle_data(Operation mode)
{
	auto range = TINY;
	
	// Shifts events away from division boundaries
	for(auto &ind : individual){	
		for(auto &ev : ind.ev){
			if(ev.tdiv == (unsigned int)(ev.tdiv)){
				switch(ev.type){
				case ENTER_EV: 
					if(ev.tdiv > TINY) ev.tdiv -= TINY;
					break;
				case LEAVE_EV: case MOVE_EV: default:
					if(ev.tdiv < T-TINY) ev.tdiv += TINY;
					break;
				}
			}
		}

		for(auto &ob : ind.obs){
			if(ob.tdiv == (unsigned int)(ob.tdiv)){
				switch(ob.type){
				case OBS_SOURCE_EV: 
					if(ob.tdiv > TINY) ob.tdiv -= TINY;
					break;
					
				case OBS_SINK_EV: case OBS_TRANS_EV: 
				default:
					if(ob.tdiv < T-TINY) ob.tdiv += TINY;
					break;
				}
			}
		}
	}
	
	Hash hash_t;
	
	// Adds times to hash time when NOT a transition event
	for(auto &ind : individual){	
		for(const auto &ev : ind.ev){
			auto vec = hash_t.get_vec_double(ev.tdiv);
			hash_t.add(0,vec);
		}
		
		for(const auto &ob : ind.obs){
			switch(ob.type){
			case OBS_SOURCE_EV: case OBS_SINK_EV: case OBS_TRANS_EV: break;
			default:
				{
					auto vec = hash_t.get_vec_double(ob.tdiv);
					hash_t.add(0,vec);
				}
				break;
			}
		}
	}
	
	for(auto &ind : individual){	
		double t_enter = UNSET, t_leave = UNSET;
		auto enter_so = UNSET, leave_so = UNSET;
		for(auto &ev : ind.ev){
			switch(ev.type){
			case ENTER_EV:
				if(t_enter != UNSET){
					auto extra = post_sim_define_twice("Enter",mode);
					alert_input("Cannot set entry time twice for individual '"+ind.name+"'."+extra);
				}				
				t_enter = ev.tdiv; enter_so = ev.so;
				break;
				
			case LEAVE_EV:
				if(t_leave != UNSET){
					auto extra = post_sim_define_twice("Leave",mode);
					alert_input("Cannot set leave time twice for individual '"+ind.name+"'."+extra);
				}				
		
				t_leave = ev.tdiv; leave_so = ev.so;
				break;
			case MOVE_EV: break;
			default: emsg_input("Should not have diferent type"); break;
			}
		}
			
		double tmin = 0.0, tmax = details.T;
	
		// Makes sure move events are in time range
		for(const auto &ev : ind.ev){
			switch(ev.type){
			case MOVE_EV:
				if(t_enter != UNSET && ev.tdiv <= t_enter){
					alert_input("The move event at time "+tstr(calc_t(ev.tdiv,details))+" for '"+ind.name+"' "+fr_ds(ev.so)+" must be after the individual enters the system "+fr_ds(enter_so)+".");
				}
				
				if(t_leave != UNSET && ev.tdiv >= t_leave){
					alert_input("The move event at time "+tstr(calc_t(ev.tdiv,details))+" for '"+ind.name+"' "+fr_ds(ev.so)+" must be before the individual leaves the system "+fr_ds(leave_so)+".");
				}
			
				if(ev.tdiv <= tmin || ev.tdiv >= tmax){	
					alert_input("The move event at time "+tstr(calc_t(ev.tdiv,details))+" for '"+ind.name+"' is not within the time range.");
				}
				break;
				
			default: break;
			}
		}
		
		for(auto &ob : ind.obs){
			if(t_enter != UNSET && ob.tdiv < t_enter){
				alert_input("The observation at time "+tstr(calc_t(ob.tdiv,details))+" on individual '"+ind.name+"' "+fr_ds(ob.so)+" must be after the individual enters the system "+fr_ds(enter_so)+".");
			}
			
			if(t_leave != UNSET && ob.tdiv > t_leave){
				alert_input("The observation at time "+tstr(calc_t(ob.tdiv,details))+" on individual '"+ind.name+"' "+fr_ds(ob.so)+" must be before the individual leaves the system "+fr_ds(leave_so)+".");
			}
			
			if(ob.tdiv < tmin || ob.tdiv > tmax){
				alert_input("Observation at time "+tstr(calc_t(ob.tdiv,details))+" on individual '"+ind.name+"' is out of range");
			}
			
			switch(ob.type){
			case OBS_SOURCE_EV: case OBS_SINK_EV: case OBS_TRANS_EV: 
				{
					double t;
					const auto loopmax = 100;
					auto loop = 0u;
					
					do{
						auto d = (double(loop)*range)/loopmax;
						
						t = ob.tdiv + d;
							
						if(t >= tmin && t <= tmax &&
   						hash_t.existing(hash_t.get_vec_double(t)) == UNSET){
							break;
						}
						
						t = ob.tdiv - d;
						if(t >= tmin && t <= tmax &&
   						hash_t.existing(hash_t.get_vec_double(t)) == UNSET){
							break;
						}
						loop++;
					}while(loop < loopmax);
					
					if(loop == loopmax){
						alert_input("Could not deal with transition event time '"+tstr(t)+"' on individual '"+ind.name+"'");
					}
					
					ob.tdiv = t;
					auto vec = hash_t.get_vec_double(t);
					hash_t.add(0,vec);
				}
				break;
				
			default: break;
			}
		}
	}
	
	if(false){
		for(auto &ind : individual){
			cout << ind.name << ": ";		
			for(auto &ev : ind.ev){
				cout << ev.tdiv << "(" << int(ev.tdiv) << "),";
			}
			cout << "time" << endl;
		}
	}
}


/// Gives a string which describes where the data comes from
string Species::fr_ds(unsigned int so) const
{
	return "(from data '"+source[so].name+"')";
}


/// Adds a waring message about adding/removing individuals twice
string Species::post_sim_define_twice(string type, Operation mode) const
{
	if(mode != PPC) return "";
	
	auto fl = false;
	for(const auto &so : source){
		switch(so.cname){
		case ADD_IND: if(type == "Enter") fl = true; break;
		case REMOVE_IND: if(type == "Leave") fl = true; break;
		default: break;
		}
	}
	
	if(fl == false) return "";
	
	if(type == "Enter"){
		return " Note, under posterior simulation individuals are added based on the posterior. The commands 'add-ind-post-sim' is only used to add additional individuals.";
	}
	
	if(type == "Leave"){
		return " Note, under posterior simulation individuals are removed based on the posterior. The commands 'remove-ind-post-sim' is only used to remove additional individuals.";
	}
	
	return "";
}

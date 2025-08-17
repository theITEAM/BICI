// Checks that inputs are correctly specified

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

#include "input.hh"
#include "utils.hh"

/// Check for initial pop (more checks get done at the end)
void Input::check_initial_pop_error(bool end)
{
	for(const auto &sp : model.species){
		auto name = sp.name;
		const auto &source = sp.source;

		auto ninitpop = 0u, naddind = 0u, naddpop = 0u;
		for(auto i = 0u; i < source.size(); i++){
			const auto &so = source[i];
			
			switch(so.cname){
			case INIT_POP: ninitpop++; break;
			case ADD_IND: naddind++; break;
			case ADD_POP: case ADD_POP_SIM: naddpop++; break;
			default: break;
			}
		}
		
		switch(model.mode){
		case SIM:
			if(ninitpop > 1 && !sp.contains_source){
				alert_import("Only one initial population should be set for species '"+name+"'");
				return;
			}
			break;
		
		case INF:
			if(ninitpop > 1){
				alert_import("Only one initial population should be set for species '"+name+"'");
				return;
			}
			
			if(end == true){
				switch(sp.type){
				case POPULATION:
					if(ninitpop == 0 && naddpop == 0 && !sp.contains_source){
						alert_import("Either 'Init. Pop.' or 'Add Pop.' must be set for species '"+name+"'");
						return;
					}
					break;
				
				case INDIVIDUAL:
					if(ninitpop == 0 && naddind == 0 && !sp.contains_source){
						alert_import("Either 'Init. Pop.' or 'Add Ind.' must be set for species '"+name+"'");
						return;
					}
					break;
				}
			}
			break;
		
		case PPC:
			break;
			
		default: emsg_input("Should not be default3a"); return;
		}
	}
}


/// Before adding a data source check that everything is correct
void Input::data_source_check_error(const DataSource &ds)
{
	const auto &tab = ds.table;
	
	if(tab.file != ""){
		// Checks that elements are correctly specified
		for(auto c = 0u; c < tab.ncol; c++){
			for(auto r = 0u; r < tab.nrow; r++){
				auto result = check_element(ds,r,c);
				if(result != ""){
					auto te = in_data_source(ds); 
					te += " the element '"+tab.ele[r][c]+"' (col '"+tab.heading[c]+"', row "+tstr(r+2)+") has the following error: "+result;
					alert_import(te);
					return;
				}
			}
		}
	}
	
	switch(ds.cname){
	case POP_TRANS_DATA:
		for(auto r = 0u; r < tab.nrow; r++){
			auto sta_str = tab.ele[r][0], end_str = tab.ele[r][1];
			auto sta = number(sta_str), end = number(end_str);
			if(sta >= end){
				alert_import(in_data_source(ds)+" on row "+tstr(r+2)+", the start time '"+sta_str+"' must be less than the end time '"+end_str+"'.");
				return;
			}
		}
		break;
		
	default: break;
	}
}


/// Checks the value of an element is correct
string Input::check_element(const DataSource &ds, unsigned int r, unsigned int c)
{
	auto col = ds.load_col[c];
	
	auto te = ds.table.ele[r][c];
	
	switch(col.type){
	case FLOAT_EL:
		{
			auto num = number(te);
			if(num == UNSET) return "The value '"+te+"' must be a number";
		}
		break;
		
	case TIME_FLOAT_EL:
		{
			if(te != "start" && te != "end" && te != "no"){
				auto num = number(te);
				if(num == UNSET) return "The value '"+te+"' must be a number";
				
				switch(ds.cname){
				case TRANS_DATA:	
					if(ds.time_start != UNSET){
						if(num <= ds.time_start || num >= ds.time_end){
							return "Time '"+te+"' must be between the observation start and end times";
						}
					}
					break;
				
				case ADD_IND:
					if(num < model.details.t_start || num > model.details.t_end){
						return "Time '"+te+"' must be between the start and end times";
					}
					break;
					
				default:
					break;
				}
			}
		}
		break;
	
	case POS_FLOAT_EL:
		{
			auto num = number(te);
			if(num == UNSET) return "The value '"+te+"' must be a number";
			if(num <= 0) return "The value '"+te+"' must be positive";
		}
		break;
		
	case ZERO_ONE_EL:
		{
			auto num = number(te);
			if(num == UNSET) return "The value '"+te+"' must be a number";
			if(num <= 0 || num >= 1) return "The value '"+te+"' must be between zero and one";
		}
		break;
		
	case POS_ZERO_FLOAT_EL:
		{
			auto num = number(te);
			if(num == UNSET) return "The value '"+te+"' must be a number";
			if(num < 0) return "The value '"+te+"' must be positive or zero";
		}
		break;
		
	case POS_INT_EL:
		{
			auto num = number(te);
			if(num == UNSET) return "The value '"+te+"' must be a number";
			if(num != int(num)) return "The value '"+te+"' must be an integer";	
			if(num < 0) return "The value '"+te+"' must be non-negative";
		}
		break;

	case POS_FLOAT_PERCENT_EL:
		{
			if(te.length() > 0 && te.substr(te.length()-1,1) == "%"){
				te = te.substr(0,te.length()-1);
			}
			
			auto num = number(te);
			if(num == UNSET) return "The value '"+te+"' must be a number";
			if(num < 0) return "The value '"+te+"' must be non-negative";
		}
		break;
	
	case TEXT_EL:
		if(te == "") return "Element empty";
		if(includes(trim(te)," ")) return "The value '"+te+"' should not include a space";
		break;
		
	case COMP_EL: 
		{  // This allows 'S', 'E' etc...
			if(find_c(ds.p,col.cl,te) == UNSET){
				auto name = model.species[ds.p].cla[col.cl].name;
				return "Compartment '"+te+"' not in classification '"+name+"'";
			}
		}
		break;
		
	case COMP_ALL_EL: 
		{
			auto result = find_comp_from_name(ds.p,te);
			if(result.error != "") return result.error;
		}
		break;
		
	case COMP_SOURCE_SINK_EL: 
		{  // Allows for compartment as well as potentially source or sink
			if(te != "+" && te != "-"){
				auto result = find_comp_from_name(ds.p,te);
				if(result.error != "") return result.error;
			}
		}
		break;
		
	case COMP_PROB_EL: case COMP_PROB_NA_EL: 
		{ // This allows 'S', 'E', 'S|E', 'S:0.4|E:0.6' etc...
			if(te != missing_str && !(col.type == COMP_PROB_NA_EL && te == not_alive_str)){
				const auto &claa = model.species[ds.p].cla[col.cl];
			
				auto spl = split_with_bracket(te,'|');
			
				auto syntax_error = false;
					
				auto colon_flag = false, normal_flag = false;
				for(auto i = 0u; i < spl.size(); i++){
					if(spl[i] == ""){ syntax_error = true; break;}
					else{
						auto spl2 = split(spl[i],':');
						
						if(spl2[0] == ""){ syntax_error = true; break;}
						
						if(find_c(ds.p,col.cl,spl2[0]) == UNSET){
							return "Compartment '"+spl2[0]+"' not in classification '"+claa.name+"'";
						}
						
						if(spl2.size() == 1){
							normal_flag = true;
						}							
						else{
							colon_flag = true; 
							if(spl2[1] == ""){ syntax_error = true; break;}
							
							auto valid = check_eqn_valid(spl2[1]);
							if(valid != SUCCESS) return "Error with equation '"+spl2[1]+"'";
							
						}
						
						if(spl2.size() > 2){ syntax_error = true; break;}
					}
				}
					
				if(normal_flag == true && colon_flag == true) syntax_error = true;

				if(syntax_error == true){
					return "There is a syntax error in table element with value '"+te+"'.";
				}
			}
		}
		break;
		
	case DIAG_TEST_EL:
		{
			auto pos = ds.obs_model.diag_pos;
			auto neg = ds.obs_model.diag_neg;
			if(te != pos && te != neg){
				return "The value '"+te+"' is neither '"+pos+"' or '"+neg+"'";
			}
		}
		break;
		
	case PRIOR_DIR_EL:
		{
			auto num = number(te);
			if(num == UNSET){	
				auto pri = convert_text_to_prior(te,line_num,"",false);
				if(pri.error != "") return pri.error;
			}	
		}
		break;
	}
	
	return "";
}


/// Determines if an equation is valid
Result Input::check_eqn_valid(string te) const 
{
	if(false) cout << te << "Check" << endl;
	return SUCCESS;
}


/// Used to order compartments
bool Comp_ord (Compartment co1, Compartment co2)                      
{ return (co1.name < co2.name); };  


/// Checks compartmental structure (looks for repeated names etc...) 
void Input::check_comp_structure()
{
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		for(auto cl = 0u; cl < sp.ncla; cl++){
			const auto &claa = sp.cla[cl];
			
			for(auto p2 = p; p2 < model.nspecies; p2++){
				const auto &sp2 = model.species[p2];
				auto cl2min = 0u; if(p == p2) cl2min = cl+1;
				for(auto cl2 =  cl2min; cl2 < sp2.ncla; cl2++){
					const auto &claa2 = sp2.cla[cl2];
					
					if(!(p == p2 && cl == cl2)){
						if(claa.name == claa2.name){
							if(p == p2){
								alert_import("The classification '"+claa.name+"' is repeated twice in '"+sp.name+"'"); return;
							}
							else{							
								if(claa.index != claa2.index){
									alert_import("The classification '"+claa.name+"' on species '"+sp.name+"' and '"+sp2.name+"' have different indices."); return;
								}
								else{
									vector <string> vec1, vec2;
									
									for(auto c = 0u; c < claa.ncomp; c++) vec1.push_back(claa.comp[c].name);
									for(auto c = 0u; c < claa2.ncomp; c++) vec2.push_back(claa2.comp[c].name);
									
									if(equal_vec(vec1,vec2) == false){
										alert_import("For classification '"+claa.name+"' the compartments do not agree between species '"+sp.name+"' and '"+sp2.name+"'"); return;
									}
								}
							}
						}
						else{
							if(claa.index == claa2.index){
								if(p == p2){
									alert_import("In species '"+sp.name+"' the classifications '"+claa.name+"' and '"+claa2.name+"' should not share the same index '"+claa.index+"'"); return;
								}
								else{
									alert_import("Species '"+sp.name+" classifications '"+claa.name+"' and species '"+sp2.name+"' classification '"+claa2.name+"' should not share the same index '"+claa.index+"'"); return;
								}
							}
						}
					}
				}
			}
		}
	}
}


/// Checks that the imported model is correct
void Input::check_import_correct()
{
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		
		if(sp.trans_tree == true && sp.infection_cl == UNSET){
			alert_line("Because 'trans-tree' is set so at least one compartment must be set as 'infected'.",sp.line_num);
			return;
		}
		
		// Checks that branching probability has not been specified in error 
		for(auto cl = 0u; cl < sp.ncla; cl++){	
			const auto &claa = sp.cla[cl];
			for(const auto &tr : claa.tra){
				if(tr.branch == false && tr.bp_set == true){
					alert_line("'bp' should not be set for this transition.",tr.line_num);
					return;
				}
			}
		}
		
		// Checks that non-Markovian transitions not used for population models
		if(sp.type == POPULATION){
			for(auto cl = 0u; cl < sp.ncla; cl++){
				const auto &claa = sp.cla[cl];
				for(const auto &tr : claa.tra){
					if(tr.type != EXP_RATE && tr.type != EXP_MEAN && tr.type != ERLANG){
						alert_line("Because species '"+sp.name+"' has a population-based model this transition must have a type 'exp(rate)', 'exp(mean)' or 'erlang'.",tr.line_num);
						return;
					}						
				}
			}
		}
	}
}


/// Fuction used for diagnostic work
void Input::temp_check(unsigned int num)
{	
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		for(auto cl = 0u; cl < sp.ncla; cl++){
			const auto &claa = sp.cla[cl];
			for(const auto &tr : claa.tra){
				if(tr.line_num > 1000){
					cout << num;
					exit (EXIT_FAILURE);
				}
			}
		}
	}
}


/// Checks if reparameterised equations involve time
void Input::check_reparam_time()
{
	for(auto th = 0u; th < model.param.size(); th++){
		auto &par = model.param[th];
		if(par.variety == REPARAM_PARAM){
			for(const auto &ele : par.element){
				if(ele.used){
					const auto &eqi = ele.value;	
					
					const auto &eq = model.eqn[eqi.eq_ref]; 
					if(eq.time_vari){
						if(eq.pop_ref.size() > 0){
							alert_import("The reparamerised parameter '"+par.full_name+"' element '"+eqi.te_raw+"' should not contain populations.");
						}
						else{
							alert_import("The reparamerised parameter '"+par.full_name+"' element '"+eqi.te_raw+"' should not contain time variation.");
						}
					}
				}
			}
		}
	}
}


/// Checks the name is valid
void Input::check_name_input(string name, string te, bool sup_not_allow) 
{
	if(name.length() > name_ch_max){
		alert_import(te+" '"+name+"' cannot be more than "+tstr(name_ch_max)+" characters",true);
		return;
	}
	
	vector <string> invalid_name = {"Compartment","Population","Alpha","Distribution","file"};
	
	for(auto i = 0u; i < invalid_name.size(); i++){
		if(toLower(name) == toLower(invalid_name[i])){
			alert_import(te+" cannot use reserved word '"+name+"'",true); 
			return;
		}
	}
	
	for(auto i = 0u; i < name.length(); i++){
		auto ch = name.substr(i,1);
		if(includes(name_notallow,ch)){
			if(ch == " ") alert_import(te+" '"+name+"' cannot contain any spaces",true);
			else alert_import(te+" '"+name+"' cannot use the character '"+ch+"'",true);
			return;
		}
		
		if(str_eq(name,i,sigma)){
			alert_import(te+" '"+name+"' cannot use the character '"+sigma+"'",true);
			return;
		}
		
		if(sup_not_allow){
			if(ch == "^"){
				alert_import(te+" '"+name+"' cannot contain the superscript character '"+ch+"'",true);
				return;
			}
		}
	}
}


/// Determines if a valid colour 
bool Input::is_Color(string color) const
{
	string allow = "0123456789abcdefABCDEF";
	color = trim(color);
	if(color.substr(0,1) == "#"){
		if(color.length() != 7) return false;
		for(auto i = 1u; i < 7; i++){
			auto ch = color.substr(i,1);
			if(!includes(allow,ch)) return false;
		}
	}
	else{
		if(color.substr(0,4) != "rgb(") return false;
		if(color.substr(color.length()-1,1) != ")") return false;
		
		auto cont = color.substr(4,color.length()-5);
		auto spl = split(cont,',');
		if(spl.size() != 3) return false;
		for(auto i = 0u; i < 3; i++){
			auto num = number(spl[i]);
			if(num == UNSET) return false;
			if(num != (int)(num)) return false;
			if(num < 0 || num > 255) return false;
		}
	}
	
	return true;
}


/// Calls when cannot find a compartment in 
void Input::cannot_find_trans_comp(string te, unsigned int p, unsigned int cl, string co)
{
	
	string st = "In transition '"+te+"' cannot find compartment '"+co+"'";
	
	if(co == ""){
		st = "In transition '"+te+"' the initial and final compartment aren't specified";
	}
	else{
		const auto &sp = model.species[p];
		if(sp.ncla > 1) st += " in classification '"+sp.cla[cl].name+"'";
		st += ".";
		
		for(auto cl2 = 0u; cl2 < sp.ncla; cl2++){
			if(cl2 != cl){
				auto c = find_c(p,cl2,co);
				if(c != UNSET){
					st += " Perhaps this transition should be placed in classification '"+sp.cla[cl2].name+"'?";
				}
			}
		}	
	}
	
	alert_import(st); 
}


/// Works out the memory used for each parameter
void Input::param_eqn_mem() const 
{
	cout << " mem use" << endl;
	for(const auto &par : model.param){
		if(par.variety == REPARAM_PARAM){
			vector <unsigned int> vec;
			for(auto &ele : par.element) vec.push_back(ele.value.eq_ref);
			eqn_mem_usage(par.name,vec);
		}
	}
	
	for(auto th = 0u; th < model.param.size(); th++){
		vector <unsigned int> vec;
		vec.push_back(th);
		param_mem_usage(model.param[th].name,vec);
	}
	
	{
		vector <unsigned int> vec;
		for(auto th = 0u; th < model.param.size(); th++) vec.push_back(th);
		param_mem_usage("Total",vec);
	}
	
	{
		vector <unsigned int> vec;
		for(auto i = 0u; i < model.eqn.size(); i++) vec.push_back(i);
		eqn_mem_usage("Total",vec);
	}
	
	param_vec_mem_usage();
	
	comp_global_mem_usage();
	
	trans_global_mem_usage();
	
	species_mem_usage();
}


/// The memory usage for the parameter vector
void Input::param_vec_mem_usage() const 
{
	vector <double> ve;
	auto vecsize = sizeof(ve)/4; 
	string ss;
	auto strsize = sizeof(ss)/4;
	
	auto norm = 0.0, affect_like = 0.0,	affect_like_me_list = 0.0, affect_like_list = 0.0, affect_like_map = 0.0, affect_like_linear_prop = 0.0;
		
	for(const auto &pele : model.param_vec){
		norm += strsize+pele.name.length();
		norm += 8;
	
		for(const auto &al : pele.affect_like){
			affect_like += 4;
			affect_like_me_list += vecsize + al.me_list.size();
			affect_like_list += vecsize + al.list.size();
			affect_like_map += vecsize + al.map.size();
			affect_like_linear_prop += vecsize + al.linear_prop.me.size();
			affect_like_linear_prop += vecsize + al.linear_prop.me_no_pop.size();
			affect_like_linear_prop += vecsize + al.linear_prop.no_pop_ref.size();
			
			affect_like_linear_prop += vecsize;
			for(const auto &pa : al.linear_prop.pop_affect){
				affect_like_linear_prop += 1 + vecsize + 2*pa.pop_grad_ref.size();
			}
		}
	}
	
	auto sum = norm + affect_like + affect_like_me_list + affect_like_list + affect_like_map +affect_like_linear_prop;
	
	cout << "Param vec	memory usage = " << 4*sum/(1000000) << "(";
	
	cout << "norm:" << int(100*norm/sum) << ",";
	cout << "affect_like:" << int(100*affect_like/sum) << ")";
	cout << "affect_like_me_list:" << int(100*affect_like_me_list/sum) << ")";
	cout << "affect_like_list:" << int(100*affect_like_list/sum) << ")";
	cout << "affect_like_map:" << int(100*affect_like_map/sum) << ")";
	cout << "affect_like_linear_prop:" << int(100*affect_like_linear_prop/sum) << ")";
	cout << endl << endl;
}


/// Estimates the memory usage for a set of equations
void Input::eqn_mem_usage(string name, const vector <unsigned int> &vec) const
{
	vector <double> ve;
	auto vecsize = sizeof(ve)/4; 
	string ss;
	auto strsize = sizeof(ss)/4;
	
	auto uni = 0.0, tex =0.0, timer = 0.0, calc = 0.0, cons = 0.0, param_ref = 0.0;
	auto pop_ref = 0.0, source_tr_gl = 0.0, ind_eff_mult = 0.0, ind_eff_exist = 0.0;
	auto fix_eff_mult = 0.0,  comp_pref = 0.0, time = 0.0, linear = 0.0;
	
	for(auto eq : vec){
		const auto &eqn = model.eqn[eq];
		
		uni += 17;
		tex += strsize*3+eqn.te.length() + eqn.te_raw.length() + eqn.warn.length();
		timer += vecsize + 2*eqn.timer.size();
		calc += vecsize; for(const auto &ca : eqn.calcu) calc += 2 + vecsize + 2*ca.item.size();
		cons += vecsize+2*eqn.cons.size();
		param_ref += vecsize+3*eqn.param_ref.size();
		pop_ref += vecsize+eqn.pop_ref.size();
		source_tr_gl += vecsize+eqn.source_tr_gl.size();
		ind_eff_mult += vecsize+eqn.ind_eff_mult.size();
		ind_eff_exist += vecsize+eqn.ind_eff_exist.size();
		fix_eff_mult += vecsize+eqn.fix_eff_mult.size();
		comp_pref += vecsize;
		for(const auto &val : eqn.comp_pref_convert) comp_pref += vecsize+val.size();
		time += vecsize+2*eqn.time.size();
		
		const auto &lin = eqn.linearise;
		
		linear += 4;
		for(auto &ca : lin.no_pop_calc) linear += 2 + vecsize + 2*ca.item.size();
		for(const auto &pgc : lin.pop_grad_calc){
			linear += vecsize;
			for(const auto &ca : pgc) linear += 2 + vecsize + 2*ca.item.size();
		}
		linear += vecsize; 
		for(const auto &prfp : lin.pop_ref_from_po){
			linear += vecsize + 2*prfp.size(); 
		}
	}
	
	auto sum = uni+tex+timer+calc+cons+param_ref+pop_ref+source_tr_gl+
							ind_eff_mult+ind_eff_exist+fix_eff_mult+comp_pref+time+linear;
	
	cout << "Memory usage for '" << name << "' = " << 4*sum/(1000000) << "(";
	
	cout << "uni:" << int(100*uni/sum) << ",";
	cout << "tex:" << int(100*tex/sum) << ",";
	cout << "timer:" << int(100*timer/sum) << ",";
	cout << "calc:" << int(100*calc/sum) << ",";
	cout << "cons:" << int(100*cons/sum) << ",";
	cout << "param_ref:" << int(100*param_ref/sum) << ",";
	cout << "pop_ref:" << int(100*pop_ref/sum) << ",";
	cout << "source_tr_gl:" << int(100*source_tr_gl/sum) << ",";
	cout << "ind_eff_mult:" << int(100*ind_eff_mult/sum) << ",";
	cout << "ind_eff_exist:" << int(100*ind_eff_exist/sum) << ",";
	cout << "fix_eff_mult:" << int(100*fix_eff_mult/sum) << ",";
	cout << "comp_pref:" << int(100*comp_pref/sum) << ",";
	cout << "time:" << int(100*time/sum) << ",";
	cout << "linear:" << int(100*linear/sum) << ")";
	cout << endl << endl;
}


/// Memory size of eqninfo
unsigned int Input::sizeof_eqninfo(const EquationInfo &ei) const
{
	string ss;
	auto strsize = sizeof(ss)/4;
	return 10 + strsize*3 + ei.te_raw.length() + ei.te.length() + ei.emsg.length();
}


/// Works out amount of memory used in parameter definition
void Input::param_mem_usage(string name, const vector <unsigned int> &vec) const
{
	vector <double> ve;
	auto vecsize = sizeof(ve)/4; 

	auto element_ref = 0.0, element = 0.0;
	for(auto th : vec){
		auto &par = model.param[th];
		
		element_ref += vecsize + par.element_ref.size();
		
		element += vecsize;
		for(auto &ele : par.element){
			// From value
			element += sizeof_eqninfo(ele.value);
			
			element++;
			
			// from parents and children
			element += vecsize + 3*ele.parent.size();
			element += vecsize + 3*ele.child.size();
			
			// others
			element += 4;
		}
	}
	
	auto sum = element + element_ref;
	
	cout << "Param usage for '" << name << "' = " << 4*sum/(1000000) << "(";
	
	cout << "element_ref:" << int(100*element_ref/sum) << ",";
	cout << "element:" << int(100*element/sum) << ")";
	cout << endl << endl;
}


/// Works out memory usage for comp model
void Input::comp_global_mem_usage() const
{
	vector <double> ve;
	auto vecsize = sizeof(ve)/4; 
	string ss;
	auto strsize = sizeof(ss)/4;
	
	auto uni = 0.0, te = 0.0, vec = 0.0, popref = 0.0, leave = 0.0;
	
	for(const auto &sp : model.species){
		for(const auto &co : sp.comp_gl){
			uni += 3;
			te += strsize + co.name.length();
			vec += 6*vecsize +  co.pop_ref_simp.size() + co.tr_enter.size() + co.tr_leave.size() 
			+ co.tr_leave_markov.size() + co.me_ref.size() + co.nmtransincomp_ref.size();
			popref += vecsize + 2*co.pop_ref.size();
			leave += vecsize;
			for(const auto &tlg : co.tra_leave_group){
				leave += 3+vecsize+tlg.tr_list.size();
			}	
		}
	}

	auto sum = uni+te+vec+popref+leave;
	
	cout << "Comp global usage = " << 4*sum/(1000000) << "(";
	cout << "uni:" << int(100*uni/sum) << ",";
	cout << "te:" << int(100*te/sum) << ",";
	cout << "vec:" << int(100*vec/sum) << ",";
	cout << "popref:" << int(100*popref/sum) << ",";
	cout << "leave:" << int(100*leave/sum) << ")";
	cout << endl;
}


/// Works out memory usage for comp model
void Input::trans_global_mem_usage() const
{
	vector <double> ve;
	auto vecsize = sizeof(ve)/4; 
	string ss;
	auto strsize = sizeof(ss)/4;
	
	auto uni = 0.0, te = 0.0, vec = 0.0, inf = 0.0, eq = 0.0;
	
	for(const auto &sp : model.species){
		for(const auto &tr : sp.tra_gl){
			uni += 16;
			te += strsize + tr.name.length();
			vec += 2*vecsize + tr.bp_other_eq.size() + tr.tform.size();
			inf += 1 + vecsize + tr.infection.pop_ref.size();
			eq += sizeof_eqninfo(tr.bp);
			for(const auto &dp : tr.dist_param) eq += sizeof_eqninfo(dp);
		}
	}

	auto sum = uni+te+vec+inf+eq;
	
	cout << "Trans global usage = " << 4*sum/(1000000) << "(";
	cout << "uni:" << int(100*uni/sum) << ",";
	cout << "te:" << int(100*te/sum) << ",";
	cout << "vec:" << int(100*vec/sum) << ",";
	cout << "eq:" << int(100*eq/sum) << ",";
	cout << endl;
}


/// Works out memory usage for comp model
void Input::species_mem_usage() const
{
	vector <double> ve;
	auto vecsize = sizeof(ve)/4; 

	auto trg_from_tr = 0.0, comp_global_convert = 0.0;
	auto obs_trans_eqn_ref = 0.0, pop_data_ref = 0.0, pop_trans_ref = 0.0;
		
	for(const auto &sp : model.species){
		for(auto ve1 : sp.trg_from_tr){
			for(auto ve2 : ve1) trg_from_tr += vecsize+ve2.size();
		}
		
		for(auto ve1 : sp.comp_global_convert){
			for(auto ve2 : ve1) comp_global_convert += vecsize+ve2.size();
		}
		
		for(auto ve1 : sp.obs_trans_eqn_ref){
			for(auto ve2 : ve1) obs_trans_eqn_ref += vecsize+ve2.size();
		}
		
		for(auto ve1 : sp.pop_data_ref){
			for(auto ve2 : ve1) pop_data_ref += vecsize+ve2.size();
		}
		
		for(auto ve1 : sp.pop_trans_ref){
			for(auto ve2 : ve1) pop_trans_ref += vecsize+ve2.size();
		}
	}


	auto sum = trg_from_tr+comp_global_convert+obs_trans_eqn_ref
			+pop_data_ref+pop_trans_ref;
	
	cout << "Species usage = " << 4*sum/(1000000) << "(";
	cout << "trg_from_tr:" << int(100*trg_from_tr/sum) << ",";
	cout << "comp_global_convert:" << int(100*comp_global_convert/sum) << ",";
	cout << "obs_trans_eqn_ref:" << int(100*obs_trans_eqn_ref/sum) << ",";
	cout << "pop_data_ref:" << int(100*pop_data_ref/sum) << ",";
	cout << "pop_trans_ref:" << int(100*pop_trans_ref/sum) << ",";
	cout << endl;
}


/// Checks if too much memory is being used
void Input::check_memory_too_large()
{
	auto mem = memory_usage();
	
#ifdef USE_MPI	
	mpi.sum(mem);
#endif

	auto f_mem = 0.0, f_free = 1.0;
	
	if(op()){
		f_mem = mem/total_memory();
		f_free = memory_available()/mem;
	}

#ifdef USE_MPI	
	mpi.bcast(f_mem);
	mpi.bcast(f_free);
#endif
	
	auto fl = false;
	/*
	if(f_mem > MEM_FRAC_MAX){
		fl = true;
		alert_import("Out of memory error. BICI currently requires '"+mem_print(mem)+"' and this computer only has '"+mem_print(mem*MEM_FRAC_MAX/f_mem)+"' available");
	}
	else{
		if(f_free < MEM_FREE_MIN){
			fl = true;
			alert_import("Out of memory error. This computer potentially has enough memory but other processes need to be terminated. BICI currently requires '"+mem_print(mem)+"' and this computer has '"+mem_print(mem*MEM_FRAC_MAX/f_mem)+"' available");
		}
	}
	*/
	
	if(fl) output_error_messages(err_mess,true);
}


/// Checks to see if derived quantities do not use derived quanties yet to be calculated  
void Input::check_derived_order()
{
	for(auto d = 0u; d < model.derive.size(); d++){
		const auto &der = model.derive[d];
		
		for(const auto &ei : der.eq){
			const auto &eqn = model.eqn[ei.eq_ref]; 
			for(const auto &ca : eqn.calcu){
				for(const auto &it : ca.item){
					if(it.type == DERIVE){
						const auto &dr = eqn.derive_ref[it.num];
						if(dr.i == d){
							alert_line("Derived quantity '"+der.full_name+"' cannot depend on itself",der.line_num);
							return;
						}
						if(dr.i > d){
							auto &der2 = model.derive[dr.i];
							alert_line("Derived quantity '"+der.full_name+"' cannot depend on derived quantity '"+der2.full_name+"' which is defined later.",der.line_num);
							return;
						}
					}
				}
			}
		}			
	} 
}

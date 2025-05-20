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
			
		default: alert_import("Should not be default3a"); return;
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
					string te = in_file_text(tab.file);
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
				alert_import(in_file_text(tab.file)+" on row "+tstr(r+2)+", the start time '"+sta_str+"' must be less than the end time '"+end_str+"'.");
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
							return "Time '"+te+"' must be between the start and end times";
						}
					}
					else{
						if(num <= model.details.t_start || num >= model.details.t_end){
							return "Time '"+te+"' must be between the start and end times";
						}
					}
					break;
				
				default:
					if(num < model.details.t_start || num > model.details.t_end){
						return "Time '"+te+"' must be between the start and end times";
					}
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
				auto pri = convert_text_to_prior(te,line_num,false);
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
			for(auto k = 0u; k < par.value.size(); k++){
				const auto &eqi = par.value[k];	
				const auto &eq = model.eqn[eqi.eq_ref]; 
				if(eq.time_vari){
					if(eq.pop_ref.size() > 0){
						alert_import("The reparamerised parameter '"+par.full_name+"' element '"+eqi.te_raw+"' should not contain populations.");
					}
					else{
						alert_import("The reparamerised parameter '"+par.full_name+"' element '"+eqi.te_raw+"' should not contain contain time variation.");
					}
				}
			}
		}
	}
}


/// Checks the name is valid
bool Input::check_name_input(string name, string te) 
{
	if(name.length() > name_ch_max){
		alert_import(te+" '"+name+"' cannot be more than "+tstr(name_ch_max)+" characters");
	}
	
	vector <string> invalid_name = {"Compartment","Population","Alpha","Distribution","file"};
	
	for(auto i = 0u; i < invalid_name.size(); i++){
		if(toLower(name) == toLower(invalid_name[i])){
			alert_import(te+" cannot use reserved word '"+name+"'"); 
		}
	}
	
	for(auto i = 0u; i < name.length(); i++){
		auto ch = name.substr(i,1);
		if(includes(name_notallow,ch)){
			alert_import(te+" '"+name+"' cannot use character '"+ch+"'");
			return true;
		}
		
		if(str_eq(name,i,sigma)){
			alert_import(te+" '"+name+"' cannot use character '"+sigma+"'");
			return true;
		}
	}
	
	return false;
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


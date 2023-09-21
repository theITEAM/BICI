///Checks that inputs are correctly specified

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

#include "input.hh"
#include "utils.hh"


/// Check for initial pop (more bech get done at end
void Input::check_initial_pop_error(bool end)
{
	for(const auto &sp : model.species){
		auto name = sp.name;
		const auto &source = sp.source;

		auto ninitpopprior = 0u, ninitpop = 0u, naddind = 0u;
		for(auto i = 0u; i < source.size(); i++){
			const auto &so = source[i];
			
			switch(so.cname){
			case INIT_POP_PRIOR: ninitpopprior++; break;
			case INIT_POP: ninitpop++; break;
			case ADD_IND: naddind++; break;
			default: break;
			}
		}
		
		switch(model.mode){
		case SIM:
			if(ninitpopprior > 0){
				alert_import("The initial population prior should not be set for species '"+name+"'");
				return;
			}
			
			/*
			if(end == true && ninitpop == 0){
				alert_import("The initial population should be set for species '"+name+"'");
				return;
			}
			*/
			
			if(ninitpop > 1){
				alert_import("Only one initial population should be set for species '"+name+"'");
				return;
			}
			break;
		
		case INF:
			if(ninitpopprior > 1){
				alert_import("Only one initial population prior should be set for species '"+name+"'");
				return;
			}
			
			if(ninitpop > 1){
				alert_import("Only one initial population should be set for species '"+name+"'");
				return;
			}
			
			if(end == true && ninitpopprior == 0 && ninitpop == 0 && naddind == 0){
				alert_import("Either 'Init. Pop.', 'Init. Pop. Prior' or 'Add Ind.' must be set for species '"+name+"'");
				return;
			}
			
			if(ninitpopprior == 1 && ninitpop == 1){
				alert_import("'Init. Pop.' and 'Init. Pop. Prior' cannot both be set for species '"+name+"'");
				return;
			}
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
					alert_import("In data file '"+tab.file+"' (col '"+tab.heading[c]+"', row "+tstr(r+2)+") the following error is found: "+result);
					return;
				}
			}
		}

// TO DO
/*
		// Checks not adding individual twice
		if(so.type == "Add Ind." || so.type == "Remove Ind."){ 
			let list=[];
			for(let r = 0u; r < tab.nrow; r++) list.push(tab.ele[r][0]);
			list.sort();
			for(let r = 0u; r < list.length-1; r++){
				if(list[r] == list[r+1]){
					let te; if(so.type == "Add Ind.") te = "added"; else te = "removed";
					
					inter.help = {title:"Error with data table", te:"Individual '"+list[r]+"' is "+te+" more than once."};
					return true;
				}
			}
		}
		
		if((so.type == "Transition" || so.type == "Source" || so.type == "Sink") && 
				so.spec.time_radio.value == "File"){
			for(let r = 0u; r < tab.nrow; r++){
				if(Number(tab.ele[r][2]) >= Number(tab.ele[r][3])){				
					inter.help = {title:"Error with data table", te:"For the observations the 'start' time must be before the 'end' time (row "+(r+2)+")."};
					return true;
				}
			}
		}
		
		if(so.type == "Diag. Test"){
			for(let r = 0u; r < tab.nrow; r++){
				let res = tab.ele[r][2];
				
				let neg_res = so.spec.neg_result;
				let pos_res = so.spec.pos_result;
				
				if(res != neg_res && res != pos_res){
					inter.help = {title:"Error with data table", te:"The test result '"+res+"' is neither negative '"+neg_res+"' or positive '"+pos_res+"' (col '"+tab.heading[2]+"', row "+(r+2)+")."};
					return true;
				}
			}
		}
		
		if(so.type == "Set Traps"){
			let list=[];
			for(let r = 0u; r < tab.nrow; r++) list.push(tab.ele[r][0]);
			list.sort();
			for(let r = 0u; r < list.length-1; r++){
				if(list[r] == list[r+1]){	
					inter.help = {title:"Error with data table", te:"In column '"+tab.heading[0]+"' the trapping event '"+list[r]+"' is set more than once."};
					return true;
				}
			}
		}
		*/
	}
	
	/*
	let desc = check_data_source_valid(so);
	if(so.error == true){
		if(desc == undefined) desc = so.error_mess;
		inter.help = {title:"Error with data table", te:desc};
		return true;
	}
	
	return false;
	*/
}


/// Checks the value of an element is correct
string Input::check_element(const DataSource &ds, unsigned int r, unsigned int c) const
{
	auto col = ds.load_col[c];
	
	auto te = ds.table.ele[r][c];
	
	switch(col.type){
	case FLOAT_EL:
		{
			auto num = number(te);
			if(num == UNSET) return "Must be a number";
		}
		break;
		
	case TIME_FLOAT_EL:
		{
			if(te != "start" && te != "end" && te != "no"){
				auto num = number(te);
				if(num == UNSET) return "Must be a number";
			}
		}
		break;
	
	case POS_FLOAT_EL:
		{
			auto num = number(te);
			if(num == UNSET) return "Must be a number";
			if(num <= 0) return "Must be positive";
		}
		break;
		
	case POS_ZERO_FLOAT_EL:
		{
			auto num = number(te);
			if(num == UNSET) return "Must be a number";
			if(num < 0) return "Must be positive or zero";
		}
		break;
		
	case POS_INT_EL:
		{
			auto num = number(te);
			if(num == UNSET) return "The Must be a number";
			if(num != int(num)) return "Must be an integer";	
			if(num < 0) return "Must be non-negative";
		}
		break;

	case POS_FLOAT_PERCENT_EL:
		{
			if(te.length() > 0 && te.substr(te.length()-1,1) == "%"){
				te = te.substr(0,te.length()-1);
			}
			
			auto num = number(te);
			if(num == UNSET) return "Must be a number";
			if(num < 0) return "Must be non-negative";
		}
		break;
	
	case TEXT_EL:
		if(te == "") return "Element empty";
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
		
	case COMP_PROB_EL: 
		{ // This allows 'S', 'E', 'S|E', 'S:0.4|E:0.6' etc...
			const auto &claa = model.species[ds.p].cla[col.cl];
		
			auto spl = split(te,'|');
		
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
						if(valid != SUCCESS) return "Error with equation'"+spl2[1]+"'";
					}
					
					if(spl2.size() > 2){ syntax_error = true; break;}
				}
			}
				
			if(normal_flag == true && colon_flag == true) syntax_error = true;

			if(syntax_error == true){
				return "There is a syntax error in table element.";
			}
		}
		break;
	
	case MULTI_COMP_EL:
		{
			const auto &claa = model.species[ds.p].cla[col.cl];
			auto spl = split(te,'|');
			
			for(auto j = 0u; j < spl.size(); j++){
				if(find_c(ds.p,col.cl,spl[j]) == UNSET){
					return "Value '"+spl[j]+"' not in classification '"+claa.name+"'";
				}
			}
		}
		break;
		
	case DIAG_TEST_EL:
		{
			auto pos = ds.obs_model.diag_pos;
			auto neg = ds.obs_model.diag_neg;
			if(te != pos && te != neg){
				return "Value '"+te+"' is neither '"+pos+"' or '"+neg+"'";
			}
		}
		break;
	}
	
	return "";
}


/// Determines if an equation is valid
Result Input::check_eqn_valid(string te) const
{
	if(false) cout << te;
	// TO DO
	return SUCCESS;
}


// Used to order compartments
bool Comp_ord (Compartment co1, Compartment co2)                      
{ return (co1.name < co2.name); };  


/// Checks compartmental structure (looks for repeated names etc...) 
void Input::check_comp_structure()
{
	// TO DO temporarily turned off
/*
	// Order compartments alphabetically
	for(auto p = 0u; p < model.nspecies; p++){
		auto &sp = model.species[p];
		for(auto cl = 0u; cl < sp.ncla; cl++){
			auto &claa = sp.cla[cl];
			
			sort(claa.comp.begin(),claa.comp.end(),Comp_ord);    
		}
	}
*/
	
	
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
		
		/// Checks that non-Markovian transitions not used for population models
		if(sp.type == POPULATION){
			for(auto cl = 0u; cl < sp.ncla; cl++){
				const auto &claa = sp.cla[cl];
				for(const auto &tr : claa.tra){
					if(tr.type != EXP_RATE && tr.type != ERLANG){
						alert_line("Because species '"+sp.name+"' has a population-based model this transition must have a type 'exp(rate)', 'exp(mean)' or 'erlang'.",tr.line_num);
						return;
					}						
				}
			}
		}
	}
}


// Fuction used for diagnostic work
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


/// Checks that the equation is formatted correctly and puts '%' and '$' around parameters
string Input::basic_equation_check(string te)
{
	if(check_brackets_match(te) == false) return te;

	if(check_chnotallowed(te) == true) return te;

	auto two_variable_flag = false;
	
	auto i = 0u; 
	while(i < te.length()){
		string type, right;
		while(i < te.length()){
			if(i < te.length()-4){
				if(te.substr(i,4) == "exp(" || te.substr(i,4) == "cos(" ||
					te.substr(i,4) == "sin("|| te.substr(i,4) == "log("){
					i += 4;
				}
			}
			
			if(i < te.length()-5 && te.substr(i,5) == "step(") i += 5;			
			
			if(i < te.length()-4 && te.substr(i,4) == "pow("){ two_variable_flag = true; i += 4;}			
			
			if(i < te.length()-7 && te.substr(i,7) == "thresh("){ two_variable_flag = true; i += 7;}
			
			if(i < te.length()-4 && te.substr(i,4) == "max("){ two_variable_flag = true; i += 4;}			
			
			if(i < te.length()-4 && te.substr(i,4) == "min("){ two_variable_flag = true; i += 4;}			
			
			if(i < te.length()-4 && te.substr(i,4) == "abs(") i += 4;			
			
			if(i < te.length()-5 && te.substr(i,5) == "sqrt(") i += 5;			
			
			if(te.substr(i,1) == "|" && two_variable_flag == true) i++;
			
			if(str_eq(te,i,"Σ")){
				while(i < te.length() && te.substr(i,1) != " " && te.substr(i,1) != "(") i++;
			}
			
			auto ch = te.substr(i,1);
			
			if(ch == "{"){ type = "pop"; right = "}"; break;}
			if(ch == "["){ type = "ie"; right = "]"; break;}
			if(ch == "<"){ type = "fe"; right = ">"; break;}
			if(ch == "}"){
				alert_import("There should not be a right bracket '}'"); //, cur:icur2+i, len:1});
				return te;
			}

			if(ch == "]"){
				alert_import("There should not be a right bracket ']'");//, cur:icur2+i, len:1});
				return te;
			}
		
			if(ch == ">"){
				alert_import("There should not be a right bracket '>'");//, cur:icur2+i, len:1});
				return te;
			}
			
			if(includes(notparam_list,ch) == false){
				type = "param"; break;
			}
			i++;
		}
			
		if(type != ""){
			auto ist = i;
			
			if(type == "param"){
				string warn;
				i = param_end(te,i,warn); 
				if(i == UNSET){
					alert_import(warn);
					return te;
				}
				//auto tex = te.substr(ist,i-ist);
			
				//auto icur3 = ist;
				//TO DO
				//check_parameter(tex,icur3,filter,eqn);

				
				te = te.substr(0,ist)+"%"+te.substr(ist,i-ist)+"$"+te.substr(i);
				
				i += 1;
			}
			else{
				while(i < te.length() && te.substr(i,1) != right) i++;
				if(i == te.length()){
					alert_import("There should be a right bracket '"+right+"'");
					return te;
				}
				else{
					//auto tex = te.substr(ist+1,i-ist-1);
					//auto icur3 = ist+1;

					if(type == "pop"){
//TO DO
					//	check_population(tex,icur3,filter,eqn);
					}

					if(type == "ie"){
//TO DO
					//	check_ie(tex,icur3,filter,eqn); 
					}
					
					if(type == "fe"){
//TO DO
					//	check_ie(tex,icur3,filter,eqn); 
					}
				}
			}
		}
		i++;
	}

	return te;
}


/// Determines when a parameter definition ends
unsigned int Input::param_end(const string &st, unsigned int i, string &warn) const 
{
	auto sub_brack = false;
	auto sup_brack = false;
	auto type = NORMAL_TEXT;
	
	do{
		auto ch = st.substr(i,1);
		if(ch == "^"){
			if(type != NORMAL_TEXT){
				warn = "Problem understanding parameter definition.";
				return UNSET;
			}
			type = SUP_TEXT;
		}
			
		if(ch == "_"){
			if(type != NORMAL_TEXT && type != SUP_TEXT){
				warn = "Problem understanding parameter definition.";
				return UNSET;
			}
			type = SUB_TEXT;
		}
	
		if(ch == "("){
			if(i > 0 && st.substr(i-1,1) == "^") sup_brack = true;
			else{
				if(i > 0 && st.substr(i-1,1) == "_") sub_brack = true;
				else break;
			}
		}
		
		if(ch == ")"){
			if(!(type == SUP_TEXT && sup_brack == true) && !(type == SUB_TEXT && sub_brack == true)) break;
			sub_brack = false;
			sup_brack = false;
		}
		
		if(includes(paramend_list,ch)) break;
		i++;
	}while(i < st.length());
	if(i <= st.length() - 3){
		if(st.substr(i,3) == "(t)" || st.substr(i,3) == "(a)") i += 3;
	}

	return i;
}			


/// Checks that all brackets match in the equation
bool Input::check_brackets_match(string te)
{
	vector <Bracket> brack_list;
	for(auto i = 0u; i < te.length(); i++){
		auto ch = te.substr(i,1);
		if(ch == "(" || ch == "[" || ch == "{"){
			Bracket br; br.i = i; br.ch = ch; brack_list.push_back(br);
		}
		
		if(ch == ")" || ch == "]" || ch == "}"){
			if(brack_list.size() == 0){
				alert_import("For equation '"+te+"': "+"The bracket '"+ch+"' does not match up.");
				return false;
			}
			
			auto chlast = brack_list[brack_list.size()-1].ch;
			if((chlast == "(" && ch == ")") || (chlast == "[" && ch == "]") ||(chlast == "{" && ch == "}")){
				brack_list.pop_back();
			}
			else{
				alert_import("For equation '"+te+"': "+"The bracket '"+chlast+"' does not match up.");
				return false;
			}
		}		
	}
	
	if(brack_list.size() != 0){
		auto last = brack_list[brack_list.size()-1];
		alert_import("For equation '"+te+"': "+"The bracket '"+last.ch+"' does not match up.");
		return false;
	}
	
	return true;
}


/// Checks that the epression does not include ceratin characters which are not allowed
bool Input::check_chnotallowed(string te)
{
	for(auto i = 0u; i < te.length(); i++){
		auto num = te.at(i);
		if(num >= 0){  // This ensures not a unicode character
			auto ch = te.substr(i,1);
			
			if(includes(chnotallowed,ch)){			
				alert_import("For equation '"+te+"': "+"The character '"+ch+"' is not allowed.");
				return true;
			}
		}
	}
	
	return false;
}

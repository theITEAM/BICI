// Utility function for equations

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <random>
#include <stdexcept>
#include "stdlib.h"
#include "math.h"
#include <sys/stat.h>
#include <cstring>
#include <signal.h>
#include <math.h>
#include <algorithm>

using namespace std;

#include "const.hh"
#include "utils.hh"

/// Adds an equation to the model
EquationInfo add_equation_info(string _te, EqnType _type, unsigned int _p, unsigned int _cl)
{
	EquationInfo info;
	info.te_raw = _te;
	info.te = _te;
	info.type = _type;
	info.p = _p;
	info.cl = _cl;
	info.value = UNSET;
	info.eq_ref = UNSET;
	info.line_num = UNSET;

	info.error = false;
	info.emsg = basic_equation_check(info.te,_type);
	if(info.emsg != "") info.error = true;
		
	return info;
}


/// Checks that a string is a valid name for a parameter / fe or ie
string check_valid_name(string name, string type)
{
	auto escape_char = get_escape_char();
	
	if(name == "" && type != "population") return "The "+type+" name is not set";	

	if(type == "parameter"){
		if(end_str(name,"(t)")) name = name.substr(0,name.length()-3);
		
		auto spl = split(name,'_');
		if(spl.size() > 1) name = spl[0];
	}
	
	if(type == "population"){
		auto spl = split(name,'<');
		if(spl.size() > 1) name = spl[0];
		
		auto spl2 = split(name,'[');
		if(spl2.size() > 1) name = spl2[0];
	}

	auto j = 0u;
	while(j < name.length()){
		auto done = false;
		
		// Checks to see if greek character
		if(j+1 < name.length()){
			for(auto k = 0u; k < escape_char.size(); k++){
				auto st = escape_char[k][1];
				auto len = st.length();
				if(name.substr(j,len) == st){
					if(st == "Σ"){
						return "In "+type+" name '"+name+"' the character '\\Sigma' is not expected";
					}
					done = true;
					j += len;
					break;
				}
			}
		}
		
		if(done == false){
			auto ch = name.substr(j,1);
			if(includes(chnotallowed,ch) || (type != "population" && includes(name_notallow,ch))){
				if(ch == "_") ch = "underscore";
				return "In "+type+" name '"+name+"' the character '"+ch+"' is not expected";
			}
			j++;
		}
	}
	return "";
}


/// Checks that the equation is formatted correctly and puts '%' and '$' around parameters
string basic_equation_check(string &te, EqnType eqn_type)
{
	auto emg = check_brackets_match(te); if(emg != "") return emg;

	emg = check_chnotallowed(te); if(emg != "") return emg;

	auto two_variable_flag = false;
	
	auto i = 0u; 
	while(i < te.length()){
		string type, right;
		while(i < te.length()){
			auto ist = i;
			
			if(i < te.length()-4){
				if(te.substr(i,4) == "exp(" || te.substr(i,4) == "cos(" ||
					te.substr(i,4) == "sin("|| te.substr(i,4) == "log("){
					i += 4;
				}
			}
			
			if(i < te.length()-5 && te.substr(i,5) == "step(") i += 5;			
			
			if(i < te.length()-4 && te.substr(i,4) == "pow("){ two_variable_flag = true; i += 4;}			
			
			if(i < te.length()-7 && te.substr(i,7) == "thresh("){ two_variable_flag = true; i += 7;}
			
			if(i < te.length()-7 && te.substr(i,7) == "ubound("){ two_variable_flag = true; i += 7;}
			
			if(i < te.length()-4 && te.substr(i,4) == "max("){ two_variable_flag = true; i += 4;}			
			
			if(i < te.length()-4 && te.substr(i,4) == "min("){ two_variable_flag = true; i += 4;}			
			
			if(i < te.length()-4 && te.substr(i,4) == "abs(") i += 4;			
			
			if(i < te.length()-5 && te.substr(i,5) == "sqrt(") i += 5;			
			
			if(i == ist){			
				if(te.substr(i,1) == "|" && two_variable_flag == true) i++;
				
				if(str_eq(te,i,"Σ")){
					while(i < te.length() && te.substr(i,1) != " " && te.substr(i,1) != "(") i++;
				}
				
				auto ch = te.substr(i,1);
				
				if(ch == "{"){ type = "pop"; right = "}"; break;}
				if(ch == "["){ type = "ie"; right = "]"; break;}
				if(ch == "<"){ type = "fe"; right = ">"; break;}
				if(ch == "}"){
					return "There should not be a right bracket '}'";
				}

				if(ch == "]"){
					return "There should not be a right bracket ']'";
				}
			
				if(ch == ">"){
					return "There should not be a right bracket '>'";
				}
				
				if(includes(notparam_list,ch) == false && !str_eq(te,i,sigma)){
					type = "param"; break;
				}
				i++;
			}
		}
	
		if(type != ""){
			auto ist = i;
			
			if(type == "param"){
				i = param_end(te,i,emg);

				if(emg != "") return emg;
				
				if(i != ist){
					auto name = te.substr(ist,i-ist);
					emg = check_valid_name(name,"parameter"); if(emg != "") return emg;
					
					te = te.substr(0,ist)+"%"+name+"$"+te.substr(i);
				}
				
				i += 1;
			}
			else{
				while(i < te.length() && te.substr(i,1) != right) i++;
				if(i == te.length()){
					return "There should be a right bracket '"+right+"'";
				}
				else{
					auto tex = te.substr(ist+1,i-ist-1);
					if(type == "pop"){
						emg = check_valid_name(tex,"population"); if(emg != "") return emg;
					}

					if(type == "ie"){
						emg = check_valid_name(tex,"individual effect"); if(emg != "") return emg;
						
						if(eqn_type == SOURCE_RATE){
							return "Source rate '"+te+"' should not contain individual effects";
						}
						
						if(eqn_type == SOURCE_MEAN){
							return "Source mean '"+te+"' should not contain individual effects";
						}
					}
					
					if(type == "fe"){
						emg = check_valid_name(tex,"fixed effect"); if(emg != "") return emg;
						
						if(eqn_type == SOURCE_RATE){
							return "Source rate '"+te+"' should not contain fixed effects";
						}
						
						if(eqn_type == SOURCE_MEAN){
							return "Source mean '"+te+"' should not contain fixed effects";
						}
					}
				}
			}
		}
		i++;
	}

	return "";
}


/// Checks that all brackets match in the equation
string check_brackets_match(string te)
{
	vector <Bracket> brack_list;
	for(auto i = 0u; i < te.length(); i++){
		auto ch = te.substr(i,1);
		if(ch == "(" || ch == "[" || ch == "{"){
			Bracket br; br.i = i; br.ch = ch; brack_list.push_back(br);
		}
		
		if(ch == ")" || ch == "]" || ch == "}"){
			if(brack_list.size() == 0){
				return "For equation '"+te+"': "+"The bracket '"+ch+"' does not match up.";
			}
			
			auto chlast = brack_list[brack_list.size()-1].ch;
			if((chlast == "(" && ch == ")") || (chlast == "[" && ch == "]") ||(chlast == "{" && ch == "}")){
				brack_list.pop_back();
			}
			else{
				return "For equation '"+te+"': "+"The bracket '"+chlast+"' does not match up.";
			}
		}		
	}
	
	if(brack_list.size() != 0){
		auto last = brack_list[brack_list.size()-1];
		return "For equation '"+te+"': "+"The bracket '"+last.ch+"' does not match up.";
	}
	
	return "";
}


/// Checks that the epression does not include certain characters which are not allowed
string check_chnotallowed(string te)
{
	for(auto i = 0u; i < te.length(); i++){
		auto num = te.at(i);
		if(num >= 0){                         // This ensures not a unicode character
			auto ch = te.substr(i,1);
			
			if(includes(chnotallowed,ch)){			
				return "For equation '"+te+"': "+"The character '"+ch+"' is not allowed.";
			}
		}
	}
	
	return "";
}


/// Determines when a parameter definition ends
unsigned int param_end(const string &st, unsigned int i, string &emg)
{
	auto sub_brack = false;
	auto sup_brack = false;
	auto type = NORMAL_TEXT;
	
	auto ist = i;
	auto ch_first = st.substr(i,1);
	
	do{
		auto ch = st.substr(i,1);
		if(ch == "^"){
			if(type != NORMAL_TEXT){
				emg = "Problem understanding parameter definition.";
				return i;
			}
			type = SUP_TEXT;
		}
			
		if(ch == "_"){
			if(type != NORMAL_TEXT && type != SUP_TEXT){
				emg = "Problem understanding parameter definition.";
				return i;
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
		
		if(includes(paramend_list,ch) || str_eq(st,i,sigma)) break;
		i++;
	}while(i < st.length());
	
	if(ch_first == "^"){
		emg = "Superscript '"+st.substr(ist,i-ist)+"' doesn't have a parameter associated with it.";
		return i;
	}
	
	if(ch_first == "_"){
		emg = "Subscript '"+st.substr(ist,i-ist)+"' doesn't have a parameter associated with it.";
		return i;
	}
	
	if(i <= st.length() - 3){
		if(st.substr(i,3) == "(t)") i += 3;
	}

	return i;
}			

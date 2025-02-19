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
	info.emsg = basic_equation_check(info.te);
	if(info.emsg != "") info.error = true;
		
	return info;
}


/// Checks that the equation is formatted correctly and puts '%' and '$' around parameters
string basic_equation_check(string &te)
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
					te = te.substr(0,ist)+"%"+te.substr(ist,i-ist)+"$"+te.substr(i);
				}
				
				i += 1;
			}
			else{
				while(i < te.length() && te.substr(i,1) != right) i++;
				if(i == te.length()){
					return "There should be a right bracket '"+right+"'";
				}
				else{
					if(type == "pop"){
						// TO DO
						// check_population(tex,icur3,filter,eqn);
					}

					if(type == "ie"){
						// TO DO
						// check_ie(tex,icur3,filter,eqn); 
					}
					
					if(type == "fe"){
						// TO DO
						// check_ie(tex,icur3,filter,eqn); 
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
				return "For equation '"+te+"': "+"The character '"+ch+"' is not allowed2.";
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
	if(i <= st.length() - 3){
		if(st.substr(i,3) == "(t)") i += 3;
	}

	return i;
}			

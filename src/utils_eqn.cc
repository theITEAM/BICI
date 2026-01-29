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
					if(st == "∫"){
						return "In "+type+" name '"+name+"' the character '\\int' is not expected";
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


/// Checks for derived functions
unsigned int der_func_check(unsigned int i, string te, string name, string &warn_der)
{
	auto namet = name+"("; 
	if(te.substr(i,namet.length()) == namet){
		auto before = trim(te.substr(0,i));
		if(before != ""){
			warn_der = "The function '"+name+"' must be on its own line. Unexpected '"+before+"' at the beginning of the line.";
		}
		
		while(i < te.length() && te.substr(i,1) != ")") i++;
		if(i == te.length()){
			warn_der = "There should be a right bracket ')' after '"+name+"('.";
		}
		
		auto after = trim(te.substr(i+1));
		if(after != ""){
			warn_der = "Unexpected '"+after+"' at end of line.";
		}
	}

	return i;
}

					
/// Checks that the equation is formatted correctly and puts '%' and '$' around parameters
string basic_equation_check(string &te, EqnType eqn_type)
{
	auto emg = check_brackets_match(te); if(emg != "") return emg;

	emg = check_chnotallowed(te); if(emg != "") return emg;

	vector <Sum> sum;

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
				
				if(str_eq(te,i,sigma)){
					auto ist = i;
					
					i += sigma.length();
					
					auto ibra = i;
					while(ibra < te.length() && te.substr(ibra,1) != "(") ibra++;
					if(ibra == te.length()){
						return "There should be a left bracket '(' after the sum";
					}
					
					auto tex = te.substr(ist,ibra-ist);
					
					if(te.substr(i,1) == "^"){
						return "The sum '"+tex+"' cannot contain a superscript '^'"; 
					}
			
					if(te.substr(i,1) != "_"){
						return "The character '_' must follow the sum 'Σ'"; 
					}
					i++;
				
					// Check distance specification
					auto ik = ist;
					while(ik < ibra && te.substr(ik,1) != "[") ik++;
					
					string te_dep;
					if(ik < ibra) te_dep = te.substr(0,ik);
					else te_dep = te;
					
					auto di = get_dep_info(te_dep,i,notparam_list);
					if(di.warn != ""){
						return "The sum 'Σ' has misspecified indices: "+di.warn;
					}
					
					if(ik < ibra){
						auto dist = trim(te.substr(ik,ibra-ik));					
						if(di.spl.size() > 1){
							return "Because the sum '"+tex+"' is limited by distance, it can only have one index";
						}
			
						if(dist.substr(0,1) != "[" || dist.substr(dist.length()-1,1) != "]"){
							return "The expression '"+dist+"' should have the format '[index,dist]'.";
						}

						auto spl = split(dist.substr(1,dist.length()-2),',');
						if(spl.size() != 2){
							return "The expression '"+dist+"' should have the format '[index,dist]'.";
						}
					
						if(di.spl[0] == spl[0]){
							return "In the sum '"+tex+"', the index and that inside the square brackets [...] cannot be the same.";
						}
				
						if(remove_prime(di.spl[0]) != remove_prime(spl[0])){
							return "In the sum '"+tex+"', the values '"+di.spl[0]+"' and '"+spl[0]+"' do not correspond to the same index.";
						}
				
						auto num = number(spl[1]);
						if(num == UNSET || num < 0){
							return "In the sum '"+tex+"', the distance '"+spl[1]+"' must be a positive number.";
						}
					}

					auto iend = ibra;
					auto num = 0;
					while(iend < te.length()){
						if(te.substr(iend,1) == "(") num++;
						if(te.substr(iend,1) == ")") num--;
						if(num == 0) break;
						iend++;
					}	
					if(iend == te.length()){
						return "The sum 'Σ' does not have a final right bracket ')'";
					}
					
					auto dep = di.spl;
					
					// Checks that index is not on surrounding sum
					for(auto k = 0u; k < sum.size(); k++){
						const auto &su = sum[k];
						if(su.i_start < ibra && su.i_end > ibra){
							for(auto d = 0u; d < dep.size(); d++){
								if(find_in(su.dep,dep[d]) != UNSET){
									return "The index '"+dep[d]+"' is already being used in a surrounding sum. Perhaps a prime is required?";
								}
							}
						}
					}
									
					Sum su; su.dep = dep; su.i_start = ibra; su.i_end = iend+1;
					sum.push_back(su);
					
					i = ibra;
				}
				
				if(str_eq(te,i,"∫")){
					while(i < te.length() && te.substr(i,1) != "(") i++;
					if(i == te.length()){
						return "There should not be a left bracket '(' after the integral";
					}
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
				string warn_der;
				i = der_func_check(i,te,RN_name,warn_der);
				i = der_func_check(i,te,RNE_name,warn_der);
				i = der_func_check(i,te,RNC_name,warn_der);
				i = der_func_check(i,te,GT_name,warn_der);
				i = der_func_check(i,te,GTE_name,warn_der);
				i = der_func_check(i,te,GTC_name,warn_der);
				if(warn_der != "") return warn_der;
				
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
						
						//if(eqn_type == REPARAM_EQN){
							//return "Reparameterisation equation cannot contain populations.";
						//}
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


/// Returns a syntax error for a line
CommandLine syntax_error(string emsg, bool fatal)
{
	CommandLine cline; 
	cline.command = EMPTY; 
	cline.emsg = emsg;
	cline.fatal = fatal;
	return cline;
}


/// Gets all the command tags from a given line
CommandLine get_command_tags(string trr, unsigned int line_num)
{
	auto spl = split(trr,' ');
	
	auto type = spl[0];
	
	Command com = EMPTY;
	if(type == "species") com = SPECIES;
	if(type == "classification" || type == "class") com = CLASS;
	if(type == "set") com = SET;
	if(type == "camera" || type == "view") com = CAMERA;
	if(type == "compartment" || type == "comp") com = COMP;
	if(type == "compartment-all" || type == "comp-all") com = COMP_ALL;
	if(type == "transition" || type == "trans") com = TRANS;
	if(type == "transition-all" || type == "trans-all") com = TRANS_ALL;
	if(type == "data-dir") com = DATA_DIR;
	if(type == "description" || type == "desc") com = DESC;
	if(type == "label") com = LABEL;
	if(type == "box") com = BOX;
	if(type == "parameter" || type == "param") com = PARAM;
	if(type == "define") com = DEFINE;
	if(type == "derived" || type == "der") com = DERIVED;
	if(type == "init-pop-inf") com = INIT_POP;
	if(type == "add-pop-inf") com = ADD_POP;
	if(type == "remove-pop-inf") com = REMOVE_POP;
	if(type == "add-ind-inf") com = ADD_IND;
	if(type == "remove-ind-inf") com = REMOVE_IND;
	if(type == "move-ind-inf") com = MOVE_IND;
	if(type == "init-pop-sim") com = INIT_POP_SIM;
	if(type == "add-pop-sim") com = ADD_POP_SIM;
	if(type == "remove-pop-sim") com = REMOVE_POP_SIM;
	if(type == "add-ind-sim") com = ADD_IND_SIM;
	if(type == "remove-ind-sim") com = REMOVE_IND_SIM;
	if(type == "move-ind-sim") com = MOVE_IND_SIM;
	if(type == "add-pop-post-sim") com = ADD_POP_POST_SIM;
	if(type == "remove-pop-post-sim") com = REMOVE_POP_POST_SIM;
	if(type == "add-ind-post-sim") com = ADD_IND_POST_SIM;
	if(type == "remove-ind-post-sim") com = REMOVE_IND_POST_SIM;
	if(type == "move-ind-post-sim") com = MOVE_IND_POST_SIM;
	if(type == "comp-data") com = COMP_DATA;
	if(type == "trans-data") com = TRANS_DATA;
	if(type == "test-data") com = TEST_DATA;
	if(type == "pop-data") com = POP_DATA;
	if(type == "pop-trans-data") com = POP_TRANS_DATA;
	if(type == "ind-effect-data") com = IND_EFFECT_DATA;
	if(type == "ind-group-data") com = IND_GROUP_DATA;
	if(type == "genetic-data") com = GENETIC_DATA;
	if(type == "simulation" || type == "sim") com = SIMULATION;
	if(type == "inference" || type == "inf") com = INFERENCE; 
	if(type == "posterior-simulation" || type == "post-sim") com = POST_SIM;
	if(type == "ind-effect") com = IND_EFFECT;
	if(type == "fixed-effect") com = FIXED_EFFECT;
	if(type == "param-sim") com = SIM_PARAM;
	if(type == "state-sim") com = SIM_STATE;
	if(type == "param-inf") com = INF_PARAM;
	if(type == "param-stats-inf") com = INF_PARAM_STATS;
	if(type == "generation-inf") com = INF_GEN;
	if(type == "state-inf") com = INF_STATE;
	if(type == "param-post-sim") com = POST_SIM_PARAM;
	if(type == "state-post-sim") com = POST_SIM_STATE;
	if(type == "diagnostics-inf") com = INF_DIAGNOSTICS;
	if(type == "map") com = MAP;
	if(type == "post-sim" || type == "post-simulation" ) com = POST_SIM;
	
	if(type == "param-mult") com = PARAM_MULT;
	
	if(type == "proposal-inf") com = PROPOSAL_INFO;
	
	if(type == "trans-diag-inf") com = TRANS_DIAG;
	
	if(type == "warning-sim") com = SIM_WARNING;
	
	if(type == "warning-inf") com = INF_WARNING;
	
	if(type == "warning-post-sim") com = PPC_WARNING;
	
	if(type == "ic") com = IC_DATA;
	
	if(com == EMPTY){
		return syntax_error("Command '"+type+"' not recognised.",true);
	}
	
	auto must_term = false; if(find_in(must_term_str,type) != UNSET) must_term = true;
		
	auto num_quote = 0u;
	vector <unsigned int> quote_pos;
			
	for(auto i = 0u; i < trr.length(); i++){
		if(trr.substr(i,1) == "\""){ num_quote++; quote_pos.push_back(i);}
	}
	
	if(num_quote%2 != 0) return syntax_error("Syntax error: Quotes do not match up.",true);
	
	for(auto i = 0u; i < num_quote; i += 2){
		if(quote_pos[i]+1 == quote_pos[i+1]){
			return syntax_error("Syntax error: No content within the quotation marks.",true);
		}
	}

	vector <Fragment> frag; 
	
	auto i = 0u; 
	auto quote = 0u;
	auto after_eq = false;
	
	do{
		auto ist = i; i++;
		
		while(i < trr.length() && !(trr.substr(i,1) == "=" && quote == 0) && 
				!(quote == 0 && trr.substr(i,1) == "\"")&& 
				!(quote == 1 && trr.substr(i,1) == "\"") && 
				!(quote == 0 && trr.substr(i,1) == " ") &&
				!(quote == 0 && after_eq == true)){
			i++;
		}

		after_eq = false;
		if(trr.substr(i,1) == "=") after_eq = true;
		
		auto te = trr.substr(ist,i-ist);
		while(te.substr(0,1) == endli) te = te.substr(1);
		while(te.substr(te.length()-1,1) == endli) te = te.substr(0,te.length()-1);

		te = trim(te); if(quote == 0) te = toLower(te);

		if(te != ""){
			auto pos = ist; while(pos < i && trr.substr(pos,1) == " ") pos++;
			Fragment fr; fr.text = trim(trr.substr(ist,i-ist)); fr.pos = pos; fr.pos_end = i; fr.quote = quote;
			frag.push_back(fr);
		}
		
		while(i < trr.length() && quote == 0 && trr.substr(i,1) == " ") i++;
		if(i < trr.length()){
			if(trr.substr(i,1) == "\""){ 
				quote = 1-quote; i++;
			}
		}
	}while(i < trr.length());
	
	if(frag.size() == 0) return syntax_error("Does not contain any content",true);
	
	if(frag[0].quote == 1) return syntax_error("Should not start with quotes",true);
	
	auto num = double(frag.size()-1)/3;
	auto numi = (unsigned int)(num);

	for(auto n = 0u; n < num; n++){
		auto ii = 1+n*3;
		if(frag[ii].text == "="){
			return syntax_error("An equals sign '=' is misplaced.",must_term);
		}
	}
	
	for(auto n = 0u; n < num; n++){
		auto ii = 1+n*3;

		if(ii+2 >= frag.size()){
			if(ii+1 < frag.size() && frag[ii+1].text == "="){
				return syntax_error("The tag '"+frag[ii].text+"' is unset",must_term);
			}
			else{
				return syntax_error("The text '"+frag[ii].text+"' cannot be understood",must_term);
			}
		}
		
		if(frag[ii+1].text != "="){
			return syntax_error("The tag '"+frag[ii].text+"' is missing an equals sign",must_term);
		}
		
		if(ii+2 < frag.size() && frag[ii+2].text == "="){
			return syntax_error("The tag '"+frag[ii].text+"' cannot be followed by '=='",must_term);
		}
		
		if(ii+3 < frag.size() && frag[ii+3].text == "="){
			return syntax_error("The tag '"+frag[ii].text+"' is unset",must_term);
		}
	}

	if(num != numi) return syntax_error("Syntax error",must_term);
	
	vector <Tag> tags;

	for(auto n = 0u; n < num; n++){
		if(frag[1+3*n+0].quote != 0) return syntax_error("Syntax error",must_term);
		if(frag[1+3*n+1].text != "=") return syntax_error("Syntax error",must_term); 
		if(frag[1+3*n+0].text == "") return syntax_error("Syntax error: Tag name not specified",must_term); 
		if(frag[1+3*n+2].text == "") return syntax_error("Tag "+frag[1+3*n+0].text+" must have content",must_term); 
		
		const auto &fr = frag[1+3*n+0];
		
		Tag tag; 
		tag.name = fr.text; tag.pos = fr.pos; tag.pos_end = fr.pos_end; 
		tag.value = frag[1+3*n+2].text; tag.done = 0;
		tags.push_back(tag);
	}
	
	if(tags.size() > 0){
		for(auto n = 0u; n < tags.size()-1; n++){
			for(auto nn = n+1; nn < tags.size(); nn++){
				if(tags[n].name == tags[nn].name){
					return syntax_error("The tag '"+tags[n].name+"' is set more than once",must_term);
				}
			}
		}
	}
		
	CommandLine cline; 
	cline.command = com; cline.command_name = type; cline.type_pos = frag[0].pos; cline.tags = tags;
	cline.line_num = line_num;
	
	return cline;
}


/// Gets a default name for a data source
void get_default_name(string &name, Command cname)
{
	switch(cname){
	case INIT_POP: name = "Initial population"; break;
	case ADD_IND: name = "Added individuals"; break;
	case REMOVE_IND: name = "Removed individuals"; break;
	case ADD_POP: name = "Added populations"; break;
	case REMOVE_POP: name = "Removed populations"; break;
	default: break;
	}
}

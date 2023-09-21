// This file deals with equations in BICI. 
// It takes a text expression and converts this to a mathematical algorithm.

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
 
#include "math.h"

using namespace std;

#include "equation.hh"
#include "utils.hh"

// Equations support the following functions:
// exp(val)
// sin(val)
// cos(val)
// log(val)
// step(val)             Goes to 1 if val > 0 otherwise zero
// pow(val|val2)         The power function
// thresh(val|val2)      if val < val2 then set to zero
// max(val|val2)         Takes the maximum of two numbers
// min(val|val2)         Takes the minimum of two numbers
// abs(val)              Takes the absolute value |val|
// sqrt(val)
    
		
/// Initialises the equation 
Equation::Equation(string tex, EqnType ty, unsigned int p, unsigned int cl, unsigned int li_num, vector <SpeciesSimp> &species, vector <Param> &param, vector <Spline> &spline, vector <ParamVecEle> &param_vec, vector <Population> &pop, const vector <double> &timepoint) : species(species), param(param), spline(spline), param_vec(param_vec), pop(pop), timepoint(timepoint)
{
	plfl = false;  // Set to true to print operations to terminal (used for diagnostics)
	
	if(plfl == true) cout << endl << tex << "  start equation\n";
	nspecies = species.size();
	
	te = tex;
	te_init = te;
	te_raw = replace(tex,"%",""); te_raw = replace(te_raw,"$",""); 

	type = ty; sp_p = p; sp_cl = cl; line_num = li_num; warn = "";
	markov_eqn_ref = UNSET;
	
	//ind_eff_mult = UNSET;
	
	te = trim(te); if(te == ""){ warn = "There is no equation"; return;}

	unravel_sum(); 
	if(warn != "") return;

	minus_sign_adjust();                             // Accounts for no number before minus sign
	
	te = replace(te,"×","*");                        // Converts × to *
	
	auto op = extract_operations();                  // Extracts the operations in the expression

  if(plfl == true) print_operations(op);

	if(warn != "") return;  

	//combine_populations(op);                         // Combines together populations (to speed up calculation)

	create_calculation(op);                          // Works out the sequence of calculation to generate result

	if(warn != "") return;  

	extract_ind_eff();                               // Extracts individual effects from equation

	if(warn != "") return;  
	
	extract_fix_eff();                               // Extracts fixed effects from equation

	if(warn != "") return;  

	check();

	if(plfl == true) print_calculation();
}


/// Checks the equation have been set up correctly
void Equation::check()
{
	for(const auto &ca : calc){
		if(ca.item.size() == 0){
			warn = "Calculation item with no elements";
			return;
		}
	}
}


/// Gets all compartments consistent with a string
// 'I'  all infectious individuls
// 'I,M' all infectious individuals
// 'E|I,M' alld infectous or exposed individuals
vector <unsigned int> Equation::get_all_comp(unsigned int p, string st)
{
	vector <unsigned int> state;

	const auto &sp = species[p]; 

	auto ncla = sp.cla.size();
	
	vector <bool> filt_set(ncla,false);
	
	vector < vector <bool> > filt;
	filt.resize(ncla);
	for(auto cl = 0u; cl < ncla; cl++){
		const auto &claa = sp.cla[cl];
		filt[cl].resize(claa.comp.size(),true);
	}

	auto vec = split(st,',');
	for(auto k = 0u; k < vec.size(); k++){
		auto st = trim(vec[k]);
		if(st != ""){
			auto vec2 = split(st,'|');
			
			vector <CompRef> comp_ref;
			auto cl = UNSET;
			for(auto j = 0u; j < vec2.size(); j++){
				auto name = vec2[j];
				
				auto cref = find_comp_from_name(p,name);
				if(cref.error != ""){ warn = cref.error; return state;}
				if(cl == UNSET) cl = cref.cl;
				else{
					if(cl != cref.cl){ 
						warn = "In population '"+st+"' cannt mix up compartments from different classifications"; 
						return state;
					}
				}
				
				comp_ref.push_back(cref);
			
			// Look for potential Erlang hidden compartments
				for(auto c = 0u; c < sp.cla[cl].ncomp; c++){
					if(sp.cla[cl].comp[c].erlang_source == name){
						cref.c = c; comp_ref.push_back(cref);
					}
				}
			}
	
			if(filt_set[cl] == true){
				warn = "Cannot put a filter on classification '"+sp.cla[cl].name+"' more than once"; 
				return state;
			}
			
			filt_set[cl] = true;
			for(auto c = 0u; c < sp.cla[cl].ncomp; c++) filt[cl][c] = false;
			for(const auto &cref : comp_ref) filt[cl][cref.c] = true;		
		}
	}
	
	for(auto c = 0u; c < sp.comp_gl.size(); c++){
		const auto &co = sp.comp_gl[c];
		
		auto cl = 0u; while(cl < ncla && filt[cl][co.cla_comp[cl]] == true) cl++;
		if(cl == ncla) state.push_back(c);
	}
/*
  auto vec = split(st,',');
  for(auto k = 0u; k < vec.size(); k++){
		auto cl = 0u;
		while(cl < ncla){
			const auto &claa = sp.cla[cl];
			auto c = 0u; while(c < claa.comp.size() && vec[k] != claa.comp[c].name) c++;
			
			if(c < claa.comp.size()){
				filt[cl] = c;
				break;
			}
			
			cl++;
		}
		
		if(cl == ncla){ warn = "Cannot find the compartment '"+vec[k]+"' in the population"; return state;}
	}
		
	vector <int> count(ncla);
  for(auto cl = 0u; cl < ncla; cl++){ if(filt[cl] == -1) count[cl] = 0u; else count[cl] = filt[cl];}

	auto j = 0u;
  do{
    auto c = 0u; for(auto cl = 0u; cl < ncla; cl++) c += sp.comp_mult[cl]*count[cl];
    state.push_back(c);

		bool flag;
    j = 0;
    do{
      flag = false;
      if(filt[j] != -1){ j++; flag = true;}
      else{
        count[j]++; if(count[j] >= (int)sp.cla[j].comp.size()){ count[j] = 0; j++; flag = true;}
      }
    }while(flag == true && j < ncla);
  }while(j < ncla);
*/
  return state;
	
}


/// Prints all the operation for a calculation
void Equation::print_operations(vector <EqItem> op) const
{
	cout << "List of operations:" << endl;
  for(auto i = 0u; i < op.size(); i++){
    switch(op[i].type){
      case LEFTBRACKET: cout << "("; break;
      case RIGHTBRACKET: cout << ")"; break;
			case FUNCDIVIDE: cout << "|"; break;
      case PARAMETER: cout << param[op[i].num].name << op[i].index; break;
			case SPLINE: cout << "Spline " << param[op[i].num].name << op[i].index; break;
			case IE: cout << species[sp_p].ind_effect[op[i].num].name; break;
			case ONE: cout << "1"; break;
			case FE: cout << species[sp_p].fix_effect[op[i].num].name; break;
      case POPNUM: cout << pop[op[i].num].name; break;
      case EXPFUNC: cout << "exp"; break;
      case SINFUNC: cout << "sin"; break;
      case COSFUNC: cout << "cos"; break;
      case LOGFUNC: cout << "log"; break;
      case STEPFUNC: cout << "step"; break;
			case POWERFUNC: cout << "pow"; break;
      case THRESHFUNC: cout << "thresh"; break; 
      case MAXFUNC: cout << "max"; break; 
			case MINFUNC: cout << "min"; break; 
			case ABSFUNC: cout << "abs"; break;
			case SQRTFUNC: cout << "sqrt"; break;
      case ADD: cout << "+"; break;
      case TAKE: cout << "-"; break;
      case MULTIPLY: cout << "*"; break;
			case DIVIDE: cout << "/"; break;
      case REG: cout << "R" << op[i].num; break;
      case NUMERIC: cout << "numeric" << op[i].constant; break;
			case TIME: cout << "time"; break;
		}
    cout << " ";
  }
  cout << endl;
}
 
 
/// Prints steps used for a calculation
void Equation::print_calculation() const   
{
	cout << "For equation '" << te_raw << "' calculation:" << endl;
	cout << "Unravelled: " << te << endl;

  for(auto i = 0u; i < calc.size(); i++){
		const auto &ca = calc[i];
		cout << ca.op << " operation\n";
		switch(ca.op){
		case EXPFUNC: cout <<  "exp("; break;
		case SINFUNC: cout <<  "sin("; break;
		case COSFUNC: cout <<  "cos("; break;
		case LOGFUNC: cout <<  "log("; break;
		case STEPFUNC: cout <<  "step("; break;
		case POWERFUNC: cout <<  "power("; break;
		case THRESHFUNC: cout <<  "thresh("; break;
		case MAXFUNC: cout << "max("; break; 
		case MINFUNC: cout << "min("; break; 
		case ABSFUNC: cout << "abs("; break;
		case SQRTFUNC: cout << "sqrt("; break;
		case ADD: break;
		case TAKE: break;
		case MULTIPLY: break;
		case DIVIDE: break;
		default: emsg("Eq problem1"); break;
		}
			
		for(auto j = 0u; j < ca.item.size(); j++){
			const auto &it = ca.item[j];
			switch(it.type){
			case PARAMETER: 
				if(it.index == UNSET) cout << param_vec[it.num].name; 
				else cout << param[it.num].name << it.index; 
				break;
			case SPLINE: 
				if(it.index == UNSET) cout << "Spline " << spline[it.num].name; 
				else{
					auto par = param[it.num]; 
					cout << "Spline " << get_param_name_with_dep(par,par.dep,it.index);
				}					
				break;
			case POPNUM: cout << "'" << pop[it.num].name << "'"; break;
			case IE: cout << species[sp_p].ind_effect[it.num].name; break;
			case ONE: cout << "1"; break;
			case FE: cout << species[sp_p].fix_effect[it.num].name; break;
			case REG: cout << "R" << it.num; break;
			case NUMERIC: cout << it.constant; break;
			case TIME: cout << "time"; break;
			default: emsg("Eq problem2"); break;
			}
			
			if(j != ca.item.size() -1){
				switch(ca.op){
				case ADD: cout <<  "+"; break;
				case TAKE: cout <<  "-"; break;
				case MULTIPLY: cout <<  "*"; break;
				case DIVIDE: cout <<  "/"; break;
				case POWERFUNC: cout << "|"; break;
				case THRESHFUNC: cout << "|"; break;
				default: break;
				}
			}
		}
	
		switch(ca.op){
		case EXPFUNC: cout << ")"; break;
		case SINFUNC: cout << ")"; break;
		case COSFUNC: cout << ")"; break;
		case LOGFUNC: cout << ")"; break;
		case STEPFUNC: cout << ")"; break;
		case POWERFUNC: cout << ")"; break;
		case THRESHFUNC: cout << ")"; break;
		default: break;
		}
			
    cout <<  " > ";

    cout <<  "R" << ca.reg_store << endl;
  }

  switch(ans.type){
    case PARAMETER:
			if(ans.index == UNSET) cout << param_vec[ans.num].name; 
			else cout << param[ans.num].name << ans.index; 
			break;
		case SPLINE:
			if(ans.index == UNSET) cout << "Spline " << spline[ans.num].name; 
			else{
				auto par = param[ans.num]; 
				cout << "Spline " << get_param_name_with_dep(par,par.dep,ans.index);
			}	
			break;
    case POPNUM: cout << pop[ans.num].name; break;
		case IE: cout << species[sp_p].ind_effect[ans.num].name; break;
    case ONE: cout << "1"; break;
		case FE: cout << species[sp_p].fix_effect[ans.num].name; break;
    case REG: cout <<  "R" << ans.num; break;
    case NUMERIC: cout << ans.constant; break;
		case TIME: cout << "time"; break;
		default: cout << " Problem here\n"; emsg("Eq problem3"); break;
  }
  cout <<  " Answer" << endl << endl;
	
	for(auto k = 0u; k < ind_eff_mult.size(); k++){
		cout << "multiplied by IE " << species[sp_p].ind_effect[ind_eff_mult[k]].name << endl;
	}
	
	for(auto k = 0u; k < fix_eff_mult.size(); k++){
		cout << "multiplied by FE " << species[sp_p].fix_effect[fix_eff_mult[k]].name << endl;
	}
}

 
/// Unravels all the sums in the equation
void Equation::unravel_sum() 
{
	const string sigma = "Σ";
	auto i = 0u;

	while(i < te.length()){
		while(i < te.length() && !str_eq(te,i,sigma)) i++;

		if(i < te.length()){
			auto istart = i;

			i += sigma.length();
			if(i == te.length()){ warn = "The character '_' must follow the sum 'Σ'"; return;}
			if(te.substr(i,1) != "_"){ warn = "The character '_' must follow the sum 'Σ'"; return;}
			i++;
			
			auto di = get_dep_info(te,i,notparam_list);
			if(di.result == FAIL){ warn = "The sum has 'Σ' has misspecified indices"; return;}

			i = di.iend;
			
			while(i < te.length() && te.substr(i,1) == " ") i++;
			
			if(i == te.length() || te.substr(i,1) != "("){
				warn = "Terms in a sum 'Σ' must be enclosed by left '(' and right ')' brackets."; 
				return;
			}
			i++;

			auto num = 1u;

			auto ist = i;
			while(i < te.length()){
				auto ch = te.substr(i,1);
				
				if(ch == ")"){ num--; if(num == 0) break; else i++;}
				else{
					if(ch == "{"){
						while(i < te.length() && te.substr(i,1) != "}") i++;
						if(i == te.length()){ warn = "Problem getting '{' bracket to match up."; return;}
					}
					else{
						if(ch == "["){
							while(i < te.length() && te.substr(i,1) != "]") i++;
							if(i == te.length()){ warn = "Problem getting '[' bracket to match up."; return;}
						}
						else{
							if(ch == "%"){
								while(i < te.length() && te.substr(i,1) != "$") i++;
								if(i == te.length()){ warn = "Problem getting '%' bracket to match up."; return;}
							}
							else{
								if(ch == "(") num++;
								i++;
							}
						}
					}
				}
			}

			if(i == te.length()){ warn = "Problem getting right ')' bracket which encloses sum."; return;}
			
			auto iend = i;

			auto content = te.substr(ist,i-ist);
	
			vector <CompPos> comp_pos;
 
			auto ndep = di.spl.size();

			for(auto j = 0u; j < ndep; j++){
				auto ind = di.spl[j];
				auto ind2 = remove_prime(ind);
				if(ind == ind2){ warn = "The index '"+ind+"' in the sum must have a prime"; return;}

				auto list = find_list_from_index(ind2);

				if(list.size() == 0){ warn = "The index '"+ind2+"' is not found within the model"; return;}
				
				CompPos cp; cp.list = list; cp.index = 0;
				comp_pos.push_back(cp);
			}

			auto te_new = te.substr(0,istart)+"( ";

			auto first = true;
			
			bool fl;			
			do{
				vector <DepConv> dep_conv;
				for(auto j = 0u; j < ndep; j++){
					DepConv dc; dc.before = di.spl[j]; dc.after = comp_pos[j].list[comp_pos[j].index];
					dep_conv.push_back(dc);
				}

				auto cont_new = content;
				
				auto res = swap_index(cont_new,dep_conv); 
				if(res.warn != ""){ warn = res.warn; return;}
				
				auto kk = find_in(res.done,false); 
				if(kk != UNSET){
					warn = "The sum index '"+dep_conv[kk].before+"' is not found in the equation";
					return;
				}

				auto k = 0u;
				do{
					fl = false;
					
					comp_pos[k].index++; 
					if(comp_pos[k].index >= comp_pos[k].list.size()){ 
						comp_pos[k].index = 0; k++; fl = true;
					}
				}while(fl == true && k < ndep);

				if(first == false) te_new += " + ";
				te_new += "("+cont_new+")";
				first = false;
			}while(fl == false);
			
			i = te_new.length();

			te_new += te.substr(iend);
			
			te = te_new;
		}	
	}
}


/// Finds the list of compartments from a given index
vector <string> Equation::find_list_from_index(string ind) const
{
	vector <string> vec;

	for(auto p = 0u; p < nspecies; p++){
		const auto &sp = species[p];
		for(auto cl = 0u; cl < sp.cla.size(); cl++){
			const auto &claa = sp.cla[cl];
			if(claa.index == ind){
				for(const auto &co : claa.comp){
					vec.push_back(co.name);
				} 
				return vec;
			}
		}
	}

	return vec;
}


/// Tries to get a number from a string
double Equation::get_float(unsigned int i, unsigned int &raend) const
{
	auto sign = 1;
	if(te.substr(i,1) == "-"){ sign = -1; i++;}

  int num = te.at(i)-48; if(!(num >= 0 && num <= 9)) return UNSET;

  auto num2 = 0.0;
  do{
    num = te.at(i)-48;
    if(num >= 0 && num <= 9) num2 = num2*10+num;
    else break;
    i++;
  }while(i < te.length());

  if(te.substr(i,1) == "."){   // looks for a decimal point
    i++;
    auto fac = 0.1;
    do{
			num = te.at(i)-48;
      if(num >= 0 && num <= 9){ num2 += fac*num; fac /= 10;}
      else break;
      i++;
    }while(i < te.length());
  }

  raend = i;
  return sign*num2;
}


/// Tries to get a parameter name from a string
ParamRef Equation::get_param_name(unsigned int i, unsigned int &raend)
{
	ParamRef pref; pref.found = false;
	
	if(te.substr(i,1) == "%"){
		auto ist = i;
		while(i < te.length() && te.substr(i,1) != "$") i++;
	
		if(i == te.length()){ warn = "Could not find right bracket '$'"; return pref;}
		
		auto content = trim(te.substr(ist+1,i-ist-1));
		
		if(content == "t"){
			pref.th = TIME_VAR;
		}
		else{
			auto pp = get_param_prop(content);
			
			/*
			auto age_dep = false, time_dep = false;
			if(content.length() > 3){
				auto end = content.substr(content.length()-3,3);
				auto fl = false;
				if(end == "(t)"){ time_dep = true; fl = true;}
				if(end == "(a)"){ age_dep = true; fl = true;}
				if(fl == true) content = content.substr(0,content.length()-3);
			}
			
			auto spl = split(content,'_');
			
			auto name = spl[0];
			*/
			auto name = pp.name;
			
			auto th = 0u; while(th < param.size() && param[th].name != name) th++;
			if(th == param.size()){ warn = "Could not find parameter '"+name+"'"; return pref;}
			
			const auto &par = param[th];
			
			if(pp.age_dep != par.age_dep){ 
				warn = "The age dependency in the equation does not agree with the definition for parameter '"+par.name+"'."; 
				return pref;
			}
			
			if(pp.time_dep != par.time_dep){ 
				warn = "The time dependency in the equation does not agree with the definition for parameter '"+par.name+"'."; 
				return pref;
			}
			
			if(pp.dep.size() != par.dep.size()){
				warn = "The dependency in the equation does not agree with the definition for parameter '"+par.name+"'."; 
				return pref;
			}
			
			if(pp.prime == true){
				warn = "\""+content+"\" should not have a prime (perhaps a sum is missing?)."; 
				return pref;
			}
			
			auto ind = 0u;
			for(auto de = 0u; de < par.dep.size(); de++){
				if(!((par.time_dep == true || par.age_dep == true) && de == par.dep.size()-1)){
					auto j = find_in(par.dep[de].list,remove_prime(pp.dep[de]));
					if(j == UNSET){
						warn = "'"+content+"' does not agree with the definition '"+par.full_name+"'."; 
						return pref;
					}
					ind += j*par.dep[de].mult;
				}
			}
			
			/*
		
			if(spl.size() > 2){ warn = "More than one '_'"; return pref;}
		
			
			if(spl.size() == 2){
				auto spl2 = split(spl[1],',');
				
				if(spl2.size() != par.dep.size()){
					warn = "The number of dependencies in the equation does not agree with the definition for parameter '"+par.name+"'."; 
					return pref;
				}

				for(auto de = 0u; de < par.dep.size(); de++){
					auto j = find_in(par.dep[de].list,spl2[de]);
					if(j == UNSET){
						warn = "'"+content+"' in the equation does not agree with the definition '"+par.full_name+"'."; 
						return pref;
					}
					ind += j*par.dep[de].mult;
				}
			}
			else{
				if(par.dep.size() != 0){
					warn = "Parameter '"+par.name+"' should have dependency defined by '"+par.full_name+"'";
					return pref;
				}
				//for(auto de = 0u; de < par.dep.size(); de++){
			}
			*/
			pref.th = th; pref.index = ind;
		}
		pref.found = true; 
		raend = i+1;
	}
	
  return pref;
}


/// Tries to get a population from the 
unsigned int Equation::get_pop(unsigned int i, unsigned int &raend)            
{
	auto p = UNSET;
	
	if(te.substr(i,1) == "{"){
		i++;
		auto ist = i;
		while(i < te.length() && te.substr(i,1) != "}") i++;
		if(i == te.length()){ warn = "The '{' bracket does not match up"; return p;}

		
		auto sp_p2 = sp_p;
		if(sp_p2 == UNSET && species.size() == 1) sp_p2 = 0;
			
		if(sp_p2 == UNSET){ warn = "A population can only be set if there is a species"; return p;}

		
		auto cont = te.substr(ist,i-ist);
		if(includes(cont,":")){
			auto spl = split(cont,':');
			if(spl.size() != 2){ warn = "A population contains more than one ':'"; return p;}
			
			auto name = spl[0];
			sp_p2 = 0; while(sp_p2 < nspecies && species[sp_p2].name != name) sp_p2++;
			if(sp_p2 == nspecies){ warn = "Species '"+name+"' does not exist"; return p;}

			cont = spl[1];
		}

		auto name = "{"+cont+"}";
	
		p = 0; while(p < pop.size() && !(name == pop[p].name && pop[p].sp_p == sp_p2)) p++;
		
		if(p == pop.size()){
			Population po;
			po.name = name;
			//po.ind_eff_mult = UNSET;
			
			auto k = 0u; while(k < cont.length() && cont.substr(k,1) != "[" && cont.substr(k,1) != "<") k++;
			
			if(k < cont.length()){
				auto extra = cont.substr(k);
			
				cont = cont.substr(0,k);	
				k = 0;
				while(k < extra.length()){
					while(k < extra.length() && extra.substr(k,1) != "[" && extra.substr(k,1) != "<") k++;
					if(k < extra.length()){
						auto ch = extra.substr(k,1);
						auto kst = k;
						
						while(k < extra.length() && extra.substr(k,1) != "]" && extra.substr(k,1) != ">") k++;
						if(k == extra.length()){ warn = "Brackets do not match up in '"+name+"'"; return p;}
						
						auto ch2 = extra.substr(k,1);
						if(ch == "["){
							if(ch2 != "]"){ warn = "Brackets do not match up in '"+name+"'"; return p;}
							
							auto ie_name = extra.substr(kst+1,k-(kst+1));
							auto &sp = species[sp_p2];
							auto ie = 0u; while(ie < sp.ind_effect.size() && sp.ind_effect[ie].name != ie_name) ie++;
							if(ie == sp.ind_effect.size()){
								IndEffect ind_eff; 
								ind_eff.name = ie_name; ind_eff.index = UNSET; ind_eff.num = UNSET; 
								ind_eff.line_num = line_num;
								sp.ind_effect.push_back(ind_eff);
							}
							po.ind_eff_mult.push_back(ie);
						}
						else{
							if(ch2 != ">"){ warn = "Brackets do not match up in '"+name+"'"; return p;}
							
							auto fe_name = extra.substr(kst+1,k-(kst+1));
	
							auto &sp = species[sp_p2];
							auto fe = 0u; while(fe < sp.fix_effect.size() && sp.fix_effect[fe].name != fe_name) fe++;
							if(fe == sp.fix_effect.size()){
								FixedEffect fix_eff; 
								fix_eff.name = fe_name;
								fix_eff.line_num = line_num;
								sp.fix_effect.push_back(fix_eff);
							}
							po.fix_eff_mult.push_back(fe);
						}
					}
				}
			}
			
			auto state = get_all_comp(sp_p2,cont);
	
			vector <unsigned int> vec;
			vector <double> vecf;

			for(auto st = 0u; st < state.size(); st++){
				vec.push_back(state[st]);
				vecf.push_back(1);
			}

			po.sp_p = sp_p2;
			for(auto i = 0u; i < vec.size(); i++){
				PopulationTerm pt; pt.c = vec[i]; pt.w = vecf[i];
				po.term.push_back(pt);
			}
			
			pop.push_back(po);
		}
		i++;
		raend = i;
	}
	
	return p;
}


/// Tries to get an individual effect from the equation
unsigned int Equation::get_ie(unsigned int i, unsigned int &raend)            
{
	auto ie = UNSET;

	if(te.substr(i,1) == "["){
		i++;
		auto ist = i;
		while(i < te.length() && te.substr(i,1) != "]") i++;
		if(i == te.length()){ warn = "The '[' bracket does not match up"; return ie;}

		if(sp_p == UNSET){ warn = "An individual effect can only be set if there is a species"; return ie;}

		auto name = te.substr(ist,i-ist);

		auto &sp = species[sp_p];
		ie = 0; while(ie < sp.ind_effect.size() && name != sp.ind_effect[ie].name) ie++;
		
		if(ie == sp.ind_effect.size()){
			IndEffect ind_eff; 
			ind_eff.name = name; ind_eff.index = UNSET; ind_eff.num = UNSET;
			ind_eff.line_num = line_num;
			sp.ind_effect.push_back(ind_eff);
		}
	
		i++;
		raend = i;
	}
	
	return ie;
}


/// Tries to get a fixed effect from the equation
unsigned int Equation::get_fe(unsigned int i, unsigned int &raend)            
{
	auto fe = UNSET;

	if(te.substr(i,1) == "<"){
		i++;
		auto ist = i;
		while(i < te.length() && te.substr(i,1) != ">") i++;
		if(i == te.length()){ warn = "The '<' bracket does not match up"; return fe;}

		if(sp_p == UNSET){ warn = "A fixed effect can only be set if there is a species"; return fe;}

		auto name = te.substr(ist,i-ist);

		const auto &fix_eff = species[sp_p].fix_effect;

		fe = 0; while(fe < fix_eff.size() && name != fix_eff[fe].name) fe++;
		if(fe == fix_eff.size()){ warn = "Fixed effect '"+name+"' is not specified"; return fe;}
		
		i++;
		raend = i;
	}
	
	return fe;
}


/*
// Works out if a particular population can be combined to others
bool Equation::pop_combine(const vector <EqItem> &op, unsigned int i,	unsigned int &popcomb, double &popcombnum) const   
{
  if(i > long(op.size())-3) return false;

	popcomb = UNSET;

  if(op[i].type == POPNUM && op[i+1].type == MULTIPLY && op[i+2].type == NUMERIC){
		popcomb = op[i].num; popcombnum = numeric[op[i+2].num];
	}

  if(op[i].type == NUMERIC && op[i+1].type == MULTIPLY && op[i+2].type == POPNUM){ 
		popcomb = op[i+2].num; popcombnum = numeric[op[i].num];
	}

  if(op[i].type == POPNUM && op[i+1].type == DIVIDE && op[i+2].type == NUMERIC){ 
		popcomb = op[i].num; popcombnum = 1.0/numeric[op[i+2].num];
	}

  if(popcomb == UNSET) return false;
  if(i > 0){
    switch(op[i-1].type){
      case MULTIPLY: case DIVIDE: return false;
      case TAKE: popcombnum *= -1; break;
			default: break;
    }
  }

  if(i+3 < op.size()){
    switch(op[i+3].type){
		case MULTIPLY: case DIVIDE: return false;
		default: break;
    }
  }
	
  return true;
}
*/


// Determines if a quantity exists
bool Equation::quant(const vector <EqItem> &op, int i) const
{
  if(i < 0 || i >= (int)op.size()) return false;
  switch(op[i].type){
	case PARAMETER: case SPLINE: case POPNUM: case IE: case FE: case REG: case NUMERIC: case TIME: return true;
  default: return false;
	}
}


// Determines if a function exists
bool Equation::is_func(const vector <EqItem> &op, int i) const
{
  if(i < 0 || i >= (int)op.size()) return false;
  switch(op[i].type){
	case EXPFUNC: case SINFUNC: case COSFUNC: case LOGFUNC: case POWERFUNC: 
	case THRESHFUNC: case STEPFUNC: case MAXFUNC: case MINFUNC: case ABSFUNC: case SQRTFUNC:
		return true;
	default: return false;
	}
}


// Determines if a operator is a particular type
bool Equation::optype(const vector <EqItem> &op, int i, EqItemType type) const
{
	if(i < 0 || i >= (int)op.size()) return false;
  if(op[i].type == type) return true;
	return false;
}


// Gives the relative priority of an operator
int Equation::prio(const vector <EqItem> &op, int i) const
{
  if(i < 0 || i >= (int)op.size()) return -1;
  switch(op[i].type){
    case ADD: return 0;
    case TAKE: return 0;
    case MULTIPLY: return 1;
    case DIVIDE: return 1;
    default: return -1;
  }
}


// Copes with minus signs which don't have a number before them
void Equation::minus_sign_adjust()
{
  auto i = 0u;                                             
  while(i < te.length()){
    if(te.substr(i,1) == "-"){
      if(i == 0 && te.substr(i,1) == "-") te = "0"+te;
      else{
        if(te.substr(i-1,1) == "(" || te.substr(i-1,1) == "|") te = te.substr(0,i)+"0"+te.substr(i);
        else i++;
      }
    }
    else i++;
  }
}


/// Takes the text and converts it into a series of operations and objects
vector <EqItem> Equation::extract_operations()
{
	vector <EqItem> op;
	
	auto i = 0u;
  while(i < te.length()){ 
		auto ch = te.substr(i,1);
		char cha = te.at(i);
		
		EqItem item;
	
		switch(cha){
      case '(': item.type = LEFTBRACKET; op.push_back(item); break;
      case ')': item.type = RIGHTBRACKET; op.push_back(item); break;
      case '|': item.type = FUNCDIVIDE; op.push_back(item); break;
      case '*': item.type = MULTIPLY; op.push_back(item); break;
      case '/': item.type = DIVIDE; op.push_back(item); break;
			case '+': item.type = ADD; op.push_back(item); break;
      case '-': item.type = TAKE; op.push_back(item); break;
			case ' ': break;

      default:
        auto doneflag = false;
				
        if(te.substr(i,4) == "exp(" && doneflag == false){
					item.type = EXPFUNC; op.push_back(item); 
          doneflag = true;
          i += 2;
        }

        if(te.substr(i,4) == "sin(" && doneflag == false){
					item.type = SINFUNC; op.push_back(item); 
          doneflag = true;
          i += 2;
        }

        if(te.substr(i,4) == "cos(" && doneflag == false){
					item.type = COSFUNC; op.push_back(item); 
          doneflag = true;
          i += 2;
        }

        if(te.substr(i,4) == "log(" && doneflag == false){
					item.type = LOGFUNC; op.push_back(item); 
          doneflag = true;
          i += 2;
        }

        if(te.substr(i,5) == "step(" && doneflag == false){
					item.type = STEPFUNC; op.push_back(item); 
					doneflag = true;
          i += 3;
        }

				if(te.substr(i,4) == "pow(" && doneflag == false){
					item.type = POWERFUNC; op.push_back(item); 
          doneflag = true;
          i += 2;
        }
				
				if(te.substr(i,7) == "thresh(" && doneflag == false){
					item.type = THRESHFUNC; op.push_back(item); 
          doneflag = true;
          i += 5;
        }
				
				if(te.substr(i,4) == "max(" && doneflag == false){
					item.type = MAXFUNC; op.push_back(item); 
          doneflag = true;
          i += 2;
        }
				
				if(te.substr(i,4) == "min(" && doneflag == false){
					item.type = MINFUNC; op.push_back(item); 
          doneflag = true;
          i += 2;
        }
				
				if(te.substr(i,4) == "abs(" && doneflag == false){
					item.type = ABSFUNC; op.push_back(item); 
          doneflag = true;
          i += 2;
        }
				
				if(te.substr(i,5) == "sqrt(" && doneflag == false){
					item.type = SQRTFUNC; op.push_back(item); 
          doneflag = true;
          i += 3;
        }
				
        if(doneflag == false){
					unsigned int raend;
          auto num = get_float(i,raend);
          if(num != UNSET){
						i = raend-1; 
						item.type = NUMERIC; 
						item.constant = num;
						op.push_back(item); 
            doneflag = true;
          }
        }

        if(doneflag == false){
					unsigned int raend;
					
          auto pref = get_param_name(i,raend); if(warn != "") return op;
          if(pref.found == true){
            i = raend-1;
						if(pref.th == TIME_VAR){
							item.type = TIME; op.push_back(item); 
						}
						else{
							auto &par = param[pref.th];
							
							if(par.time_dep == true){
								if(par.spline_info.on != true) emsg("Spline should be on");
								item.type = SPLINE;
							}							
							else item.type = PARAMETER; 
							
							if(item.type == PARAMETER && param[pref.th].variety == CONST_PARAM){
								item.type = NUMERIC; 
								item.constant = param[pref.th].value[pref.index].value;
								op.push_back(item); 
							}
							else{
								item.num = pref.th; 
								item.index = pref.index; 
								op.push_back(item); 
							}
							
							auto ii = 0u; 
							while(ii < param_ref.size() && !(param_ref[ii].th == pref.th && param_ref[ii].index == pref.index)) ii++;
			
							if(ii == param_ref.size()){
								ParamRef parref; parref.th = pref.th; parref.index = pref.index; 
								param_ref.push_back(parref);
							}
							
							par.used = true;
						}
            doneflag = true;
          }
        }

        if(doneflag == false){
					unsigned int raend;
          auto numi = get_pop(i,raend); if(warn != "") return op;
          if(numi != UNSET){
            i = raend-1;
						item.type = POPNUM; item.num = numi; op.push_back(item); 
						
						add_to_vec(pop_ref,numi);
						
						doneflag = true;
          }
        }
				
				if(doneflag == false){
					unsigned int raend;
          auto ie = get_ie(i,raend); if(warn != "") return op;
          if(ie != UNSET){
            i = raend-1;
						item.type = IE; item.num = ie; op.push_back(item); 
						doneflag = true;
          }
        }

				if(doneflag == false){
					unsigned int raend;
          auto fe = get_fe(i,raend); if(warn != "") return op;
          if(fe != UNSET){
						i = raend-1;
            item.type = FE; item.num = fe; op.push_back(item); 
						doneflag = true;
          }
        }
				
				if(doneflag == false){
					cout << "Problem with expression. The character '"+ch+"' was not expected.\n";
				}
        if(doneflag == false){ warn = "Problem with expression. The character '"+ch+"' was not expected."; return op;}
        break;
    }
    i++;
  }

  i = 0;                                            // Adds in unspecified multiply signs
  while(i+1 < op.size()){
		if((quant(op,i) || optype(op,i,RIGHTBRACKET)) && 
		    (quant(op,i+1) || is_func(op,i+1) || optype(op,i+1,LEFTBRACKET))){		 
			EqItem item; item.type = MULTIPLY;
      op.insert(op.begin()+i+1,item);
    }
    else i++;
  }
	
	return op;
}


/// Looks to combine together multiple populations into one (where possible) to help speed up calculation
// TO DO
/*
void Equation::combine_populations(vector <EqItem> &op) const
{

	if(sp_p != 
	
  i = 0;                                            expressions
  while(i < int(op.size())-3){
		unsigned int popcomb;
		double popcombnum;
    
    if(pop_combine(op,i,popcomb,popcombnum) == true){
      auto ist = i;
			vector <double> w(N,0);

      stringstream ss;
		  while(pop_combine(op,i,popcomb,popcombnum) == true){
        ss << pop[popcomb].name << "*" << popcombnum;

        for(j = 0; j < pop[popcomb].term.size(); j++){
          w[pop[popcomb].term[j].c] += popcombnum*pop[popcomb].term[j].w;
        }
        i += 3;

        if(i == op.size()) break;
        if(op[i].type != ADD && op[i].type != TAKE) break;
        i++;
        ss << "+";
      }
      
			auto name = ss.str();
    
			auto p = 0u; while(p < pop.size() && name != pop[p].name) p++;
			if(p == pop.size()){
				Population po;
				po.name = name;
				
				for(auto c = 0u; c < N; c++){
					if(w[c] != 0){
						PopulationTerm popt; popt.c = c; popt.w = w[c];
						po.term.push_back(popt); 
					}
        }
				
				pop.push_back(po);
      }
			
      op[ist].type = POPNUM; op[ist].num = p; op[ist].index = UNSET;
      op.erase(op.begin()+ist+1,op.begin()+i);
      i = ist;
    }
    i++;
  }
}
*/

/// Works out the sequence of calculation to generate result
void Equation::create_calculation(vector <EqItem> &op)
{
	nreg = 0;
	
	auto flag = false;
  do{                                               // Breaks down the equations into simple pieces
    flag = false;
		
    for(auto i = 0; i < int(op.size())-1; i++){
      switch(op[i].type){
			case LEFTBRACKET: case RIGHTBRACKET:
				break;
				
			case ADD: case TAKE: case MULTIPLY: case DIVIDE:
				if(i > 0 && quant(op,i-1) == true && quant(op,i+1) == true && prio(op,i) >= prio(op,i-2) && prio(op,i) >= prio(op,i+2)){
					auto add_flag = false;
					if(calc.size() > 0){
						auto &cal = calc[calc.size()-1];
						if(cal.op == op[i].type && (cal.op == ADD || cal.op == MULTIPLY)){
							if(op[i-1].type == REG && op[i-1].num == cal.reg_store){
								cal.item.push_back(op[i+1]);
								add_flag = true;
							} 
						}
					}

					if(add_flag == true){   // Combines with previous operation
						op.erase(op.begin()+i,op.begin()+i+2);
						i--;
					}
					else{                   // Inserts a new calculation
						Calculation cal;
						cal.item.push_back(op[i-1]);
						cal.item.push_back(op[i+1]);
					
						cal.op = op[i].type;
						cal.reg_store = nreg;
					
						calc.push_back(cal);
					
						i--;
						op[i].type = REG;
						op[i].num = nreg;
				
						op.erase(op.begin()+i+1,op.begin()+i+3);
						nreg++;
					}
					flag = true;
				}
				break;
				
			case PARAMETER: case SPLINE: case POPNUM:               // Removes brackets around a quantity
			case IE: case FE: case REG: case NUMERIC: case TIME: 
				if(optype(op,i-1,LEFTBRACKET) && optype(op,i+1,RIGHTBRACKET) && !is_func(op,i-2)){
					i--;
					op[i] = op[i+1];
					op.erase(op.begin()+i+1,op.begin()+i+3);
			 
					flag = true;
				}
				break;
				
			case ONE: break;
			
			case FUNCDIVIDE: break;
			
			case EXPFUNC: case LOGFUNC: case SINFUNC: case COSFUNC: case STEPFUNC: case ABSFUNC: case SQRTFUNC: // Univariate
				if(optype(op,i+1,LEFTBRACKET) && quant(op,i+2) == true && optype(op,i+3,RIGHTBRACKET)){		
					Calculation cal;
					cal.item.push_back(op[i+2]);
					
					cal.op = op[i].type;
					cal.reg_store = nreg;
					
					calc.push_back(cal);
					
					op[i].type = REG;
					op[i].num = nreg;
				
					op.erase(op.begin()+i+1,op.begin()+i+4);
					nreg++;
					
					flag = true;
				}
				break;
				
			case POWERFUNC: case THRESHFUNC: case MAXFUNC: case MINFUNC:  // Functions with two variables
				if(optype(op,i+1,LEFTBRACKET) && quant(op,i+2) == true && optype(op,i+3,FUNCDIVIDE) && quant(op,i+4) == true && optype(op,i+5,RIGHTBRACKET)){		
					Calculation cal;
					cal.item.push_back(op[i+2]);
					cal.item.push_back(op[i+4]);
					
					cal.op = op[i].type;
					cal.reg_store = nreg;
					
					calc.push_back(cal);
					
					op[i].type = REG;
					op[i].num = nreg;
					
					op.erase(op.begin()+i+1,op.begin()+i+6);
					nreg++;
					
					flag = true;
				}
      }
    }

    if(plfl == true){ cout << "OPERERATIONS AFTER \n"; print_operations(op);}
  }while(flag == true && op.size() > 1);

	if(op.size() != 1){ warn = "Problem with expression. Calculation cannot be performed."; return;}

	ans = op[0];
}
	

/// Calculates the value for an equation
double Equation::calculate_param_only(const vector <double> &param_val) const 
{
 	vector <double> regcalc(nreg);

  for(auto i = 0u; i < calc.size(); i++){
		const auto &ca = calc[i];
		
		const auto &item = ca.item;
		const auto N = item.size();
		
		vector <double> num(N);
		
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
			
			switch(it.type){
				case PARAMETER:
					num[j] = param_val[it.num]; 
					if(it.index != UNSET) emsg("SHould be unset1"); 
					if(num[j] == UNSET) emsg("Param must be set");
					break;
				case SPLINE: emsg("Should not be here"); break;
				case IE: emsg("Should not include ind effect"); break;
				case ONE: num[j] = 1; break;
				case FE: emsg("Should not include fixed effect"); break;
				case POPNUM: emsg("Should not include population"); break;
				case REG: num[j] = regcalc[it.num]; break;
				case NUMERIC: num[j] = it.constant; break;
				case TIME: emsg("Should not include time"); break;
				default: emsg("Equation error"); break;
			}
		}

    regcalc[ca.reg_store] = calculate_operation(ca.op,num);
  }

  switch(ans.type){
    case PARAMETER: 
			if(ans.index != UNSET) emsg("Should be unset2"); 
			if(param_val[ans.num] == UNSET) emsg("Param must be set");
			return param_val[ans.num];
    case SPLINE: emsg("SHould not be here"); return UNSET;
    case POPNUM: emsg("SHould not be here"); return UNSET;
    case REG: return regcalc[ans.num];
    case NUMERIC: return ans.constant;
		case IE: emsg("Should not include ind effect"); return UNSET;
		case ONE: return 1; 		
		case FE: emsg("Should not include fixed effect"); return UNSET;
		case TIME: emsg("Should not include time"); return UNSET;
		default: emsg("Equation error"); return UNSET;
  }
}


/// Calculates the value for an equation
double Equation::calculate(unsigned int ti, const vector <double> &popnum, const vector <double> &param_val, const vector <SplineValue> &spline_val) const 
{
 	vector <double> regcalc(nreg);

  for(auto i = 0u; i < calc.size(); i++){
		const auto &ca = calc[i];
		
		const auto &item = ca.item;
		const auto N = item.size();
		
		vector <double> num(N);
		
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
			
			switch(it.type){
				case PARAMETER:
					num[j] = param_val[it.num]; 
					if(it.index != UNSET) emsg("SHould be unset3"); 
					break;
				case SPLINE:
					num[j] = spline_val[it.num].val[ti];
					if(it.index != UNSET) emsg("SHould be unset4"); 
					break;
				case IE: emsg("Should not include ind effect"); break;
				case ONE: num[j] = 1; break;
				case FE: emsg("Should not include fixed effect"); break;
				case POPNUM: num[j] = rectify(popnum[it.num]); break;
				case REG: num[j] = regcalc[it.num]; break;
				case NUMERIC: num[j] = it.constant; break;
				case TIME: num[j] = timepoint[ti]; break;
				default: emsg("Equation error"); break;
			}
		}

    regcalc[ca.reg_store] = calculate_operation(ca.op,num);
  }

  switch(ans.type){
    case PARAMETER:
			if(ans.index != UNSET) emsg("Should be unset5"); 
			return param_val[ans.num];
    case SPLINE: 
			if(ans.index != UNSET) emsg("SHould be unset6"); 
			return spline_val[ans.num].val[ti];
    case POPNUM: return rectify(popnum[ans.num]);
    case REG: return regcalc[ans.num];
    case NUMERIC: return ans.constant;
		case IE: emsg("Should not include ind effect"); return UNSET;
		case ONE: return 1; 		
		case FE: emsg("Should not include fixed effect"); return UNSET;
		case TIME: return timepoint[ti];
		default: emsg("Equation error"); return UNSET;
  }
}

/// Calculates the value for an equation
double Equation::calculate_indfac(const Individual &ind, unsigned int ti, const vector <double> &popnum, const vector <double> &param_val, const vector <SplineValue> &spline_val) const 
{	
	return indfac(ind)*calculate(ti,popnum,param_val,spline_val);
}


/// Calculates the individual factor for the equation
double Equation::indfac(const Individual &ind) const
{
	auto num = 1.0; 
	for(auto ie : ind_eff_mult) num *= ind.exp_ie[ie];
	for(auto fe : fix_eff_mult) num *= ind.exp_fe[fe];
	
	return num;
}


/// Performs an operation of a set of numbers 
double Equation::calculate_operation(EqItemType op, vector <double> &num) const 
{
	auto N = num.size();
	
	switch(op){
	case ADD:
		{
			auto ans = 0.0; 
			for(auto val : num) ans += val;
			return ans;
		}
		
	case TAKE:
		if(N != 2) emsg("For TAKE should be 2");
		return num[0]-num[1]; 
	
	case MULTIPLY:
		{
			auto ans = 1.0; 
			for(auto val : num) ans *= val;
			return ans;
		}
		
	case DIVIDE:
		if(N != 2) emsg("For DIVIDE should be 2");
		if(num[1] == 0) emsg("Equation '"+te_raw+"' caused a division by zero."); 
		return num[0]/num[1]; 
		
	case EXPFUNC:
		if(N != 1) emsg("For EXPFUNC should be 1");
		return exp(num[0]); 
			
	case SINFUNC:
		if(N != 1) emsg("For SINFUNC should be 1");
		return sin(num[0]); 
	
	case COSFUNC: 
		if(N != 1) emsg("For COSFUNC should be 1");
		return cos(num[0]); 
		
	case LOGFUNC: 
		if(N != 1) emsg("For LOGFUNC should be 1");
		if(num[0] <= 0) emsg("For equation '"+te_raw+"' the quantity inside the log became negative."); 
		return log(num[0]);
		
	case STEPFUNC:
		if(N != 1) emsg("For STEPFUNC should be 1");
		if(num[0] > 0) return 1;
		return 0; 
		
	case POWERFUNC:
		if(N != 2) emsg("For POWERFUNC should be 2");
		return pow(num[0],num[1]);
	
	case THRESHFUNC:
		if(N != 2) emsg("For POWERFUNC should be 2");
		if(num[0] < num[1]) return 0;
		return num[0];
		
	case MAXFUNC:
		if(N != 2) emsg("For MAXFUNC should be 2");
		if(num[0] > num[1]) return num[0];
		return num[1]; 
		
	case MINFUNC:
		if(N != 2) emsg("For MINFUNC should be 2");
		if(num[0] < num[1]) return num[0];
		return num[1]; 
		
	case ABSFUNC:
		if(N != 1) emsg("For MINFUNC should be 1");
		if(num[0] > 0) return num[0]; 
		return -num[0];
		
	case SQRTFUNC:
		if(N != 1) emsg("For SQRTFUNC should be 1");
		if(num[0] < 0) emsg("For equation '"+te_raw+"' the quantity inside the square root became negative."); 
		return sqrt(num[0]);
		
	default: emsg("Equation error"); break;
	}
	
	return 0;
}


/// Extracts individual effects from the equation
void Equation::extract_ind_eff()
{
	do{
		auto ie_num = UNSET;
		auto reg_used  = UNSET;
		auto ist = UNSET, jst = UNSET;
		
		for(auto i = 0u; i < calc.size(); i++){
			const auto &ca = calc[i];
			
			const auto &item = ca.item;
			const auto N = item.size();
			
			for(auto j = 0u; j < N; j++){
				const auto &it = item[j];
				
				if(it.type == IE){
					ie_num = it.num;
					ist = i; jst = j;
					
					if(ca.op != MULTIPLY && !(ca.op == DIVIDE && j == 0)){
						warn = "Equation must by multiplied by individual effect";
						return;
					}
					reg_used = ca.reg_store;
				}
				
				if(it.type == REG && it.num == reg_used){
					if(ca.op != MULTIPLY && !(ca.op == DIVIDE && j == 0)){
						warn = "Equation must by multiplied by individual effect";
						return;
					}
					reg_used = ca.reg_store;
				}
			}
		}
		
		if(ie_num != UNSET){
			ind_eff_mult.push_back(ie_num);
			if(calc[ist].op == MULTIPLY && calc[ist].item.size() > 1){
				calc[ist].item.erase(calc[ist].item.begin()+jst);
			}
			else calc[ist].item[jst].type = ONE;
		}
		else break;
	}while(true);
}


/// Extracts fixed effects from the equation
void Equation::extract_fix_eff()
{
	do{
		auto fe_num = UNSET;
		auto reg_used  = UNSET;
		auto ist = UNSET, jst = UNSET;
	
		for(auto i = 0u; i < calc.size(); i++){
			const auto &ca = calc[i];
			
			const auto &item = ca.item;
			const auto N = item.size();
			
			for(auto j = 0u; j < N; j++){
				const auto &it = item[j];
				
				if(it.type == FE && fe_num == UNSET){
					fe_num = it.num;
					ist = i; jst = j;
					
					if(ca.op != MULTIPLY && !(ca.op == DIVIDE && j == 0)){
						warn = "Equation must by multiplied by fixed effect";
						return;
					}
					reg_used = ca.reg_store;
					break;
				}
				
				if(it.type == REG && it.num == reg_used){
					if(ca.op != MULTIPLY && !(ca.op == DIVIDE && j == 0)){
						warn = "Equation must by multiplied by fixed effect";
						return;
					}
					reg_used = ca.reg_store;
				}
			}
		}
		
		if(fe_num != UNSET){
			fix_eff_mult.push_back(fe_num);
			if(calc[ist].op == MULTIPLY && calc[ist].item.size() > 1){
				calc[ist].item.erase(calc[ist].item.begin()+jst);
			}
			else calc[ist].item[jst].type = ONE;
		}
		else break;
	}while(true);
}


/// Finds a compartment based on a name
CompRef Equation::find_comp_from_name(unsigned int p, string te) const 
{
	te = trim(te);
	if(te != ""){ 
		const auto &sp = species[p];
		for(auto cl = 0u; cl < sp.cla.size(); cl++){
			const auto &claa = sp.cla[cl];
			for(auto c = 0u; c < claa.ncomp; c++){
				if(claa.comp[c].name == te){
					CompRef cr; cr.p = p; cr.cl = cl; cr.c = c; cr.error = "";
					return cr;
				}
			}
		}
	}
	
	CompRef cr; cr.p = UNSET; cr.cl = UNSET; cr.c = UNSET; cr.error = "Compartment '"+te+"' not specfied";
	return cr;
}


// This file deals with equations in BICI. 
// It takes a text expression and converts this to a mathematical algorithm.

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm> 
 
#include "math.h"

using namespace std;

#include "equation.hh"
#include "utils.hh"

// Equations support the following functions:
// exp(val)              Exponential
// sin(val)              Sine function
// cos(val)              Cosine function
// log(val)              Logarithm function
// step(val)             Goes to 1 if val > 0 otherwise zero
// pow(val|val2)         The power function
// thresh(val|val2)      if val < val2 then set to zero
// ubound(val|val2)      if val > val2 then set to infinity
// max(val|val2)         Takes the maximum of two numbers
// min(val|val2)         Takes the minimum of two numbers
// abs(val)              Takes the absolute value |val|
// sqrt(val)             Takes the square root of a value
    
		
/// Initialises the equation 
Equation::Equation(string tex, EqnType ty, unsigned int p, unsigned int cl, unsigned int c, bool inf_trans, unsigned int tif, unsigned int li_num, vector <SpeciesSimp> &species, vector <Param> &param, vector <Spline> &spline, vector <ParamVecEle> &param_vec, vector <Population> &pop, Hash &hash_pop, const vector <double> &timepoint) : species(species), param(param), spline(spline), param_vec(param_vec), pop(pop), hash_pop(hash_pop), timepoint(timepoint)
{
	plfl = false;  // Set to true to print operations to terminal (used for diagnostics)
	
	if(plfl) cout << endl << tex << "  start equation" << endl;

	nspecies = species.size();
	
	type = ty; sp_p = p; sp_cl = cl; sp_c = c; ti_fix = tif, line_num = li_num; warn = "";
	infection_trans = inf_trans;
	markov_eqn_ref = UNSET;
	ans.type = NUMERIC; ans.num = UNSET;
	
	tex = trim(tex); if(tex == ""){ warn = "There is no equation"; return;}

	te = tex;
	te_init = te;
	te_raw = te;
	te_raw = replace(te_raw,"%","");
	te_raw = replace(te_raw,"$","");
	//if(ti_fix != UNSET) te_raw += "(t="+to_string(timepoint[ti_fix])+")";
	
	if(warn != "") return;

	unravel_sum(); 
	if(warn != "") return;

	te = replace(te,"×","*");                        // Converts × to *
	
	auto op = extract_operations();                  // Extracts the operations in the 	expression
	
	check_repeated_operator(op);                     // Checks for repeated operators e.g. "**"
	
	if(warn != "") return; 
	
	replace_minus(op);                               // Replaces minums sign with (-1)*
	
	if(plfl) print_operations(op);

	create_calculation(op);                          // Works out the sequence of calculation to generate result

	if(warn != "") return; 
	
	if(plfl == true) print_calculation();

	if(simplify_eqn == true) simplify();             // Simplifies by combining constants

	if(warn != "") return;  

	extract_ind_eff();                               // Extracts individual effects from equation

	if(warn != "") return;  
	
	extract_fix_eff();                               // Extracts fixed effects from equation

	if(warn != "") return;  

	check();

	if(plfl == true) print_calculation();

	linearise.on = false;

	ind_eff = false; 
	if(ind_eff_mult.size() > 0 || fix_eff_mult.size() > 0) ind_eff = true;
	
	if(plfl) print_calculation();
	
	set_time_vari();
	
	// Truncates strings to save memory
	te = trunc(te);
	te_raw = trunc(te_raw);
	
	//print_calculation();
}


/// Sets up a conversion from p,c to pref (used for creates proposal in the transition tree)
void Equation::setup_comp_pref_convert()
{
	comp_pref_convert.resize(species.size());
	for(auto p = 0u; p < species.size(); p++){
		const auto &sp = species[p];
		if(sp.trans_tree){
			comp_pref_convert[p].resize(sp.comp_gl.size(),UNSET);
		}
	}
	
	for(auto pref = 0u; pref < pop_ref.size(); pref++){
		const auto &po = pop[pop_ref[pref]];
		auto p = po.sp_p;
		
		for(const auto &te : po.term){
			auto c = te.c;
			
			if(comp_pref_convert[p][c] != UNSET){
				warn = "For transition a compartment cannot appear in more than one population";
			}
			
			comp_pref_convert[p][c] = pref; 
		}
	}
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
vector <unsigned int> Equation::get_all_comp(unsigned int p, string te)
{
	vector <unsigned int> state;
	
	const auto &sp = species[p]; 

	auto ncla = sp.cla.size();
	
	vector <bool> filt_set(ncla,false);

	vector < vector <unsigned int> > list_all;
	list_all.resize(ncla);

	auto vec = split(te,',');
	for(auto k = 0u; k < vec.size(); k++){
		auto st = trim(vec[k]);
		if(st != ""){
			auto vec2 = split(st,'|');
			
			vector <unsigned int> list;
			auto cl = UNSET;
			for(auto j = 0u; j < vec2.size(); j++){
				auto name = vec2[j];
				
				auto cref = find_comp_from_name(p,name);
				if(cref.error != ""){ 
					warn = "In population '{"+te+"}': "+cref.error; 
					return state;
				}
				
				if(cl == UNSET) cl = cref.cl;
				else{
					if(cl != cref.cl){ 
						warn = "In population '{"+te+"}' cannot mix up compartments from different classifications"; 
						return state;
					}
				}
				
				list.push_back(cref.c);
			
				// Look for potential Erlang hidden compartments
				for(auto c = 0u; c < sp.cla[cl].ncomp; c++){
					if(sp.cla[cl].comp[c].erlang_source == name){
						list.push_back(c);
					}
				}
			}
	
			if(filt_set[cl] == true){
				warn = "In population '{"+te+"}' cannot put a filter on classification '"+sp.cla[cl].name+"' more than once"; 
				return state;
			}
			
			filt_set[cl] = true;
			list_all[cl] = list;
		}
	}
	
	for(auto cl = 0u; cl < ncla; cl++){
		if(filt_set[cl] == false) list_all[cl] = seq_vec(sp.cla[cl].ncomp);
	}
	
	state = list_all[0];
	
	const auto &comp_mult = sp.comp_mult;
	for(auto cl = 1u; cl < ncla; cl++){
		auto li = list_all[cl];
		auto mu = comp_mult[cl];
		for(auto &val : li) val *= mu;

		vector <unsigned int> state_new;
		for(auto val : state){
			for(auto val2 : li){
				state_new.push_back(val+val2);
			}
		}
		state = state_new;
	}
	
  return state;
}


/// Prints all the operations for a calculation
void Equation::print_operations(const vector <EqItem> &op) const
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
      case UBOUNDFUNC: cout << "ubound"; break; 
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
			case NOOP: cout << "No operation"; break;
		}
    cout << " ";
  }
  cout << endl;
}
 
 
/// Prints steps used for a calculation
void Equation::print_calculation() const   
{
	cout << "For equation '" << te_raw << "' calculation:" << endl;

  for(auto i = 0u; i < calc.size(); i++){
		print_ca(calc[i]);
		cout << endl;
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
		case IE: cout << "[" << species[sp_p].ind_effect[ans.num].name << "]"; break;
    case ONE: cout << "1"; break;
		case FE: cout << "<" << species[sp_p].fix_effect[ans.num].name << ">"; break;
    case REG: cout <<  "R" << ans.num; break;
    case NUMERIC: cout << ans.constant; break;
		case TIME: cout << "time"; break;
		default: emsg("Eq problem3"); break;
  }
  cout <<  " Answer" << endl << endl;
	
	for(auto k = 0u; k < ind_eff_mult.size(); k++){
		cout << "multiplied by IE " << species[sp_p].ind_effect[ind_eff_mult[k]].name << endl;
	}
	
	for(auto k = 0u; k < fix_eff_mult.size(); k++){
		cout << "multiplied by FE " << species[sp_p].fix_effect[fix_eff_mult[k]].name << endl;
	}
}

void Equation::print_ca(const Calculation &ca) const
{  
	switch(ca.op){
	case EXPFUNC: cout <<  "exp("; break;
	case SINFUNC: cout <<  "sin("; break;
	case COSFUNC: cout <<  "cos("; break;
	case LOGFUNC: cout <<  "log("; break;
	case STEPFUNC: cout <<  "step("; break;
	case POWERFUNC: cout <<  "power("; break;
	case THRESHFUNC: cout <<  "thresh("; break;
	case UBOUNDFUNC: cout <<  "ubound("; break;
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
		case NUMERIC: 
			if(it.constant == INFY) cout << "INFY";
			else{
				if(it.constant == UNDEF) cout << "UNDEF";
				else cout << it.constant; 
			}
			break;
		case TIME: cout << "time"; break;
		default: emsg("Eq problem2"); break;
		}
		
		if(j != ca.item.size()-1){
			switch(ca.op){
			case ADD: cout <<  "+"; break;
			case TAKE: cout <<  "-"; break;
			case MULTIPLY: cout <<  "*"; break;
			case DIVIDE: cout <<  "/"; break;
			case POWERFUNC: cout << "|"; break;
			case THRESHFUNC: cout << "|"; break;
			case UBOUNDFUNC: cout << "|"; break;
			case MAXFUNC: cout << "|"; break;
			case MINFUNC: cout << "|"; break;
		
			default: break;
			}
		}
	}

	switch(ca.op){
	case EXPFUNC: case SINFUNC: case COSFUNC: case LOGFUNC: case STEPFUNC: case POWERFUNC: 
	case THRESHFUNC: case UBOUNDFUNC: case MAXFUNC: case MINFUNC: case ABSFUNC: case SQRTFUNC:
		cout << ")"; 
		break;
	default: break;
	}
		
	cout <<  " > ";
	cout <<  "R" << ca.reg_store;
}
 
 
/// Unravels all the sums in the equation
void Equation::unravel_sum() 
{
	auto i = 0u;

	while(i < te.length()){
		while(i < te.length() && !str_eq(te,i,sigma)) i++;

		if(i < te.length()){
			auto istart = i;

			i += sigma.length();
			if(i == te.length()){ warn = "The character '_' must follow the sum 'Σ'"; return;}

			double distmax = UNSET;
			if(te.substr(i,1) == "^"){
				i++;

				auto ist = i;
				while(i < te.length() && te.substr(i,1) != "_") i++;
				if(i == te.length()){ warn = "The character '_' must follow the sum 'Σ'"; return;}
				
				auto sup = te.substr(ist,i-ist);
				
				auto fl = false;
				if(sup.length() < 4) fl = true;
				else{
					if(sup.substr(0,4) != "max:") fl = true;
					auto num = sup.substr(4);
					distmax = number(num);
					if(distmax == UNSET){
						warn = "Sum superscript number '"+num+"' not recognised.";
						return;
					}
				}
				
				if(fl == true){
					warn = "Sum superscript '"+sup+"' not recognised.";
					return;
				}
			}
			
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

				auto list = find_list_from_index(ind2,distmax);
	
				if(list.size() == 0){ warn = "The index '"+ind2+"' is not found within the model"; return;}
				
				CompPos cp; cp.list = list; cp.index = 0;
				comp_pos.push_back(cp);
			}

			auto te_new = te.substr(0,istart)+"( ";

			auto first = true;
			
			vector <DepConv> dep_conv;
			for(auto j = 0u; j < ndep; j++){
				DepConv dc; dc.before = di.spl[j];
				dep_conv.push_back(dc);
			}
			
			auto swap_temp = swap_template(content,dep_conv);
			if(swap_temp.warn != ""){ warn = swap_temp.warn; return;}
	
			auto kk = find_in(swap_temp.done,false); 
			if(kk != UNSET){
				warn = "The sum index '"+dep_conv[kk].before+"' is not found in the equation";
				return;
			}

			bool fl;			
			do{
				for(auto j = 0u; j < ndep; j++){
					dep_conv[j].after = comp_pos[j].list[comp_pos[j].index];
				}

				auto cont_new = swap_index_temp(dep_conv,swap_temp);
		
				if(check_swap){
					auto cont_new_ch = content;
					auto res = swap_index(cont_new_ch,dep_conv); 
					if(res.warn != ""){ warn = res.warn; return;}
					if(cont_new != cont_new_ch){
						cout << cont_new << " " << cont_new_ch << " compare" << endl; 
						emsg("Swap index dif res");
					}
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
vector <string> Equation::find_list_from_index(string ind, double dist_max) const
{
	vector <string> vec;

	if(dist_max != UNSET){
		if(sp_p == UNSET || sp_c == UNSET) emsg("Error with 'max' function. The sum index '"+ind+"' does not appear in the equation.");

		const auto &sp = species[sp_p];
		for(auto cl = 0u; cl < sp.cla.size(); cl++){
			const auto &claa = sp.cla[cl];
			if(claa.index == ind){
				const auto &co = sp.comp_gl[sp_c];
				auto c = co.cla_comp[cl];
				const auto &comp = claa.comp;
				for(auto cc = 0u; cc < comp.size(); cc++){
					auto d = find_dist(c,cc,comp,claa.coord); 
					if(d < dist_max) vec.push_back(comp[cc].name);
				}
				return vec;
			}
		}	
	}
					
	for(auto p = 0u; p < nspecies; p++){
		const auto &sp = species[p];
		for(auto cl = 0u; cl < sp.cla.size(); cl++){
			const auto &claa = sp.cla[cl];
			if(claa.index == ind){
				for(const auto &co : claa.comp) vec.push_back(co.name);
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
ParamRef Equation::get_param_name(unsigned int i, double &dist, unsigned int &raend)
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
		
			auto name = pp.name;
			
			auto th = 0u; while(th < param.size() && param[th].name != name) th++;
			if(th == param.size()){
				if(name == dist_matrix_name){
					dist = get_distance(pp);
				}
				else{
					warn = "Could not find parameter '"+name+"'"; return pref;
				}
			}
			
			if(dist == UNSET){
				const auto &par = param[th];
				
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
					if(!(par.time_dep == true && de == par.dep.size()-1)){
						//const auto &hash_list = par.dep[de].hash_list;
						//auto vec = hash_list.get_vec_string(remove_prime(pp.dep[de]));
						//auto j = hash_list.existing(vec);
						auto j = par.dep[de].hash_list.find(remove_prime(pp.dep[de]));
						if(j == UNSET){
							warn = "'"+content+"' does not agree with the definition '"+par.full_name+"'."; 
							emsg("PPP ");
							return pref;
						}
						ind += j*par.dep[de].mult;
					}
				}
				
				pref.th = th; pref.index = ind;
			}
		}
		pref.found = true; 
		raend = i+1;
	}
	
  return pref;
}


/// Tries to get a population from the equation
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

		string name;
		if(species.size() > 1) name += species[sp_p2].name+":";
		name += "{"+cont+"}";

		auto vec = hash_pop.get_vec_string(name);
		p = hash_pop.existing(vec);
		if(p == UNSET){
			p = pop.size();
			hash_pop.add(pop.size(),vec);
			
			Population po;
			po.name = name;
			
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
								if(ie_name  == ""){ warn = "Individual effect in population '[]' must contain text"; return p;}
									
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
	
							if(fe_name  == ""){ warn = "Fixed effect in population '<>' must contain text"; return p;}
							
							auto &sp = species[sp_p2];
							auto fe = 0u; 
							while(fe < sp.fix_effect.size() && sp.fix_effect[fe].name != fe_name) fe++;
							if(fe == sp.fix_effect.size()){
								warn = "Fixed effect '<"+fe_name+">' unspecified. This is specified through the 'fixed-effect' command"; return p;
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
			
			po.ind_variation = false;
			if(po.ind_eff_mult.size() > 0 || po.fix_eff_mult.size() > 0) po.ind_variation = true;
			
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
		if(name == ""){
			warn = "Individual effect '[]' must contain text";
			return ie;
		}
		
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

		if(name == ""){
			warn = "Fixed effect '<"+name+">' must contain text";
			return fe;
		}
		
		const auto &fix_eff = species[sp_p].fix_effect;

		fe = 0; while(fe < fix_eff.size() && name != fix_eff[fe].name) fe++;
		if(fe == fix_eff.size()){ warn = "Fixed effect '<"+name+">' is not specified"; return fe;}
		
		i++;
		raend = i;
	}
	
	return fe;
}


/// Determines if a quantity exists
bool Equation::quant(const vector <EqItem> &op, int i) const
{
  if(i < 0 || i >= (int)op.size()) return false;
  switch(op[i].type){
	case PARAMETER: case SPLINE: case POPNUM: case IE: case FE: case REG: case NUMERIC: case TIME: return true;
  default: return false;
	}
}


/// Determines if a function exists
bool Equation::is_func(const vector <EqItem> &op, int i) const
{
  if(i < 0 || i >= (int)op.size()) return false;
  switch(op[i].type){
	case EXPFUNC: case SINFUNC: case COSFUNC: case LOGFUNC: case POWERFUNC: 
	case THRESHFUNC: case UBOUNDFUNC: case STEPFUNC: case MAXFUNC: case MINFUNC: case ABSFUNC: case SQRTFUNC:
		return true;
	default: return false;
	}
}


/// Determines if a operator is a particular type
bool Equation::optype(const vector <EqItem> &op, int i, EqItemType type) const
{
	if(i < 0 || i >= (int)op.size()) return false;
  if(op[i].type == type) return true;
	return false;
}


/// Gives the relative priority of an operator
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
      case '-': item.type = TAKE; op.push_back(item);	break;
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
				
				if(te.substr(i,7) == "ubound(" && doneflag == false){
					item.type = UBOUNDFUNC; op.push_back(item); 
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
					
					double dist = UNSET;
				
          auto pref = get_param_name(i,dist,raend); if(warn != "") return op;
				
          if(pref.found == true){
            i = raend-1;
						
						if(dist != UNSET){
							item.type = NUMERIC; item.constant = dist; op.push_back(item); 
						}
						else{
							if(pref.th == TIME_VAR){
								if(ti_fix != UNSET){
									item.type = NUMERIC; item.constant = timepoint[ti_fix]; op.push_back(item); 
								}
								else{
									item.type = TIME; op.push_back(item); 
								}
							}
							else{
								auto &par = param[pref.th];
								
								if(par.time_dep == true){
									if(par.spline_info.on != true) emsg("Spline should be on");
									item.type = SPLINE;
								}							
								else item.type = PARAMETER; 
								
								auto done = false;
								
								if(done == false){ 
									if(item.type == PARAMETER && param[pref.th].variety == CONST_PARAM){
										item.type = NUMERIC; 
										const auto &par = param[pref.th];
										item.constant = par.value[pref.index].value;
										op.push_back(item); 
										done = true;
									}
								}
									
								if(done == false){
									if(item.type == PARAMETER && param[pref.th].variety == REPARAM_PARAM){
										auto num = number(param[pref.th].value[pref.index].te);
										if(num != UNSET){
											item.type = NUMERIC; 
											item.constant = num;
											op.push_back(item); 
											done = true;
										}
									}
								}
								
								if(done == false){
									item.num = pref.th; 
									item.index = pref.index; 
									op.push_back(item); 
									
									// Adds in factor to equation
									if(par.param_mult != UNSET){
										{
											EqItem item2; item2.type = MULTIPLY; op.push_back(item2); 
										}
										item.type = SPLINE; 
										item.num = par.param_mult;
										if(par.time_dep){
											const auto &dep = par.dep;
											item.index /= dep[dep.size()-1].list.size();
										}
										
										const auto &par_mult = param[par.param_mult];
										const auto &dep = par_mult.dep;
										
										item.index *= dep[dep.size()-1].list.size();
										
										op.push_back(item); 
									}
								}
								
								par.used = true;
							}
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
				
				//if(doneflag == false){
				//cout << "Problem with expression. The character '"+ch+"' was not expected." << endl;
				//}
				
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


/// Works out the sequence of calculations to generate result
void Equation::create_calculation(vector <EqItem> &op)
{
	nreg = 0;
	
	auto flag = false;
  do{                                             // Breaks down the equations into simple pieces
    flag = false;
		
    for(auto i = 0; i < int(op.size())-1; i++){
      switch(op[i].type){
			case LEFTBRACKET: case RIGHTBRACKET:
				break;
				
			case ADD: case TAKE: case MULTIPLY: case DIVIDE:
				if(op[i].type == ADD && (i == 0 || op[i-1].type == LEFTBRACKET)){ // Removes (+...)
					op.erase(op.begin()+i);
					flag = true;
				}
				else{
					if(i > 0 && quant(op,i-1) == true && quant(op,i+1) == true && 
					    prio(op,i) >= prio(op,i-2) && prio(op,i) >= prio(op,i+2)){
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

						if(add_flag == true){                 // Combines with previous operation
							op.erase(op.begin()+i,op.begin()+i+2);
							i--;
						}
						else{                                 // Inserts a new calculation
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
				}
				break;
				
			case PARAMETER: case SPLINE: case POPNUM: // Removes brackets around a quantity
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
				
			case POWERFUNC: case THRESHFUNC:            // Functions with two variables
			case UBOUNDFUNC: case MAXFUNC: case MINFUNC:
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
				break;
				
			case NOOP: emsg("Hould not be no op"); break;
      }
    }

    if(plfl == true) print_operations(op);
  }while(flag == true && op.size() > 1);

 	string wa = "";
		
	if(op.size() != 1){
		vector <EqItemType> oper;
		for(const auto &val : op){
			switch(val.type){
			case ADD: case TAKE: case MULTIPLY: case DIVIDE:
				oper.push_back(val.type);
				break;
			default: break;
			}
		}
	
		if(oper.size() == 1){
			wa = "There is potentially a stray "+op_name(oper[0])+" in the expression. ";
		}
	}

	if(wa != ""){
		warn = "Problem with expression. ";
		if(wa != "error") warn += wa;	
		warn += "Calculation cannot be performed.";
		return;
	}

	ans = op[0];
}


/// Calculates the value for an equation (assuming equation only has parameters)
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
				case SPLINE: emsg("Should not be here1"); break;
				case IE: emsg("Should not include ind effect"); break;
				case ONE: num[j] = 1; break;
				case FE: emsg("Should not include fixed effect"); break;
				case POPNUM: emsg("Should not include population"); break;
				case REG: num[j] = regcalc[it.num]; break;
				case NUMERIC: num[j] = it.constant; break;
				case TIME: emsg("Should not include time1"); break;
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
    case SPLINE: emsg("SHould not be here2"); return UNSET;
    case POPNUM: emsg("SHould not be here3"); return UNSET;
    case REG: return regcalc[ans.num];
    case NUMERIC: return ans.constant;
		case IE: emsg("Should not include ind effect"); return UNSET;
		case ONE: return 1; 		
		case FE: emsg("Should not include fixed effect"); return UNSET;
		case TIME: emsg("Should not include time2"); return UNSET;
		default: emsg("Equation error"); return UNSET;
  }
}


/// Calculates the value for an equation using only parameters and splines (if ti_fix is set)
double Equation::calculate_param_only_ti_fix(const vector <double> &param_val, const vector <SplineValue> &spline_val) const 
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
				case SPLINE:
					if(ti_fix == UNSET) emsg("Should not be here4");
					num[j] = spline_val[it.num].val[ti_fix];
					if(it.index != UNSET) emsg("SHould be unset4"); 
					break;
				case IE: emsg("Should not include ind effect"); break;
				case ONE: num[j] = 1; break;
				case FE: emsg("Should not include fixed effect"); break;
				case POPNUM: emsg("Should not include population"); break;
				case REG: num[j] = regcalc[it.num]; break;
				case NUMERIC: num[j] = it.constant; break;
				case TIME: emsg("Should not include time3"); break;
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
    case SPLINE: 
			if(ti_fix == UNSET) emsg("Should not be here5");
			return spline_val[ans.num].val[ti_fix];
    case POPNUM: emsg("SHould not be here6"); return UNSET;
    case REG: return regcalc[ans.num];
    case NUMERIC: return ans.constant;
		case IE: emsg("Should not include ind effect"); return UNSET;
		case ONE: return 1; 		
		case FE: emsg("Should not include fixed effect"); return UNSET;
		case TIME: emsg("Should not include time4"); return UNSET;
		default: emsg("Equation error"); return UNSET;
  }
}



/// Calculates the value for an equation (but does not allow population number)
double Equation::calculate_no_popnum(unsigned int ti, const vector <double> &param_val, const vector <SplineValue> &spline_val) const 
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
				case POPNUM: emsg("Should not include popnum"); break;
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
    case POPNUM: emsg("Should not include popnum"); return UNSET;
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


/// Calculates the value for an equation (including individual-based factor)
double Equation::calculate_indfac(const Individual &ind, unsigned int ti, const vector <double> &popnum, const vector <double> &param_val, const vector <SplineValue> &spline_val) const 
{	
	return indfac(ind)*calculate(ti,popnum,param_val,spline_val);
}


/// Performs an operation on a set of numbers 
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
		if(N != 2){
			emsg("For TAKE should be 2");
		}
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
		if(N != 2) emsg("For THRESHFUNC should be 2");
		if(num[0] < num[1]) return 0;
		return num[0];
		
	case UBOUNDFUNC:
		if(N != 2) emsg("For UBOUNDFUNC should be 2");
		if(num[0] > num[1]) return INFINITY;
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
		
	case NOOP: 
		return 0;
	
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
						warn = "Entire equation must by multiplied by individual effect ["+species[sp_p].ind_effect[it.num].name+"]";
						return;
					}
					reg_used = ca.reg_store;
				}
				
				if(it.type == REG && it.num == reg_used){
					if(ca.op != MULTIPLY && !(ca.op == DIVIDE && j == 0)){
						warn = "Entire equation must by multiplied by individual effect ["+species[sp_p].ind_effect[it.num].name+"]";
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
						warn = "Entire equation must by multiplied by fixed effect <"+species[sp_p].fix_effect[it.num].name+">";
						return;
					}
					reg_used = ca.reg_store;
					break;
				}
				
				if(it.type == REG && it.num == reg_used){
					if(ca.op != MULTIPLY && !(ca.op == DIVIDE && j == 0)){
						warn = "Entire equation must by multiplied by fixed effect <"+species[sp_p].fix_effect[it.num].name+">";
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


/// Finds a compartment from a name
CompRef Equation::find_comp_from_name(unsigned int p, string te) const 
{
	te = trim(te);
	
	if(te != ""){ 
		const auto &sp = species[p];
		for(auto cl = 0u; cl < sp.cla.size(); cl++){
			const auto &claa = sp.cla[cl];

			//auto vec = claa.hash_comp.get_vec_string(te);
			//auto c = claa.hash_comp.existing(vec);
			auto c = claa.hash_comp.find(te);
			if(c != UNSET){
				CompRef cr; cr.p = p; cr.cl = cl; cr.c = c; cr.error = "";
				return cr;
			}
		}
	}
	
	CompRef cr; cr.p = UNSET; cr.cl = UNSET; cr.c = UNSET; cr.error = "Compartment '"+te+"' not specfied";
	return cr;
}


/// Replaces a register with a constant 
void Equation::replace_reg(unsigned int num, EqItem it_change)
{
	for(auto i = num+1; i < calc.size(); i++){
		auto &ca = calc[i];
		
		for(auto j = 0u; j < ca.item.size(); j++){
			auto &it = ca.item[j];
			if(it.type == REG){
				if(it.num == num){
					it = it_change;
				}
				else{
					if(it.num > num){
						it.num--;
					}
				}
			}
		}
		calc[i].reg_store--;
	}	

	if(ans.type == REG){
		if(ans.num == num){
			ans = it_change;
		}
		else{
			if(ans.num > num) ans.num--;
		}
	}
	nreg--;
}


/// Simplifies equations as much as possible
// (1) By removing zeros
// (2) Making used of infinity / undefined
// (3) Combining together addition and multiplication terms 
void Equation::simplify()
{
	if(false){
		cout << " START SIMPLIFY" << endl; print_calculation();
	}
	
	auto loop = 0u;
	bool flag;
	do{		
		flag = false;
	
		// Converts a negative TAKE constant to an ADD
		for(auto i = 0u; i < calc.size(); i++){
			auto &ca = calc[i];
			if(ca.op == TAKE && ca.item[1].type == NUMERIC){
				ca.op = ADD; ca.item[1].constant *= -1;
				flag = true;
			}
		}
		
		// Combines together ADD, MULT from two different calculations and removes
		if(calc.size() > 0){
			auto reg_used = get_reg_used(calc);
		
			for(int i = ((int) calc.size())-1; i >= 0; i--){
				auto i_un = (unsigned int) i;
				auto &ca = calc[i];
			
				auto reg_rem = false;
				
				switch(ca.op){
				case ADD: case MULTIPLY:
					if(!(ans.type == REG && ans.num == ca.reg_store)){		
						if(reg_used[ca.reg_store].size() == 1){	
							auto ii = reg_used[ca.reg_store][0];
							auto &ca2 = calc[ii];
							if(ca2.op == ca.op){
								reg_rem = true;
								for(auto j = 0u; j < ca.item.size(); j++){
									const auto &it = ca.item[j];
									if(it.type == REG){
										auto &ru = reg_used[it.num];
										auto k = 0u; 
										while(k < ru.size() && ru[k] != i_un) k++;
										if(k == ru.size()) emsg("cannot find");
										ru[k] = ii;
									}
									ca2.item.push_back(it);
								}
							}
						}
					}
					break;
				
				default: break;
				}
				
				if(reg_rem == true){
					flag = true;
					EqItem rep; rep.type = NUMERIC; 
					switch(ca.op){
					case ADD: rep.constant = 0; break;
					case MULTIPLY: rep.constant = 1; break;
					default: emsg("prob"); break;
					}
					rep.index = UNSET; rep.num = UNSET;
					replace_reg(ca.reg_store,rep);
					calc.erase(calc.begin()+i);
						
					reg_used.erase(reg_used.begin()+i);		
					for(auto ii = 0u; ii < reg_used.size(); ii++){
						for(auto &va : reg_used[ii]) if(va > i_un) va--;
					}
				}
			}
			
			if(false){
				for(auto ii = 0u; ii < reg_used.size(); ii++){
					cout << ii << ": ";
					for(auto k : reg_used[ii]) cout << k << ",";
					cout << "regused" << endl;
				}
			}
			
			if(testing) check_reg_used(calc,reg_used);
		}
		
		for(int i = ((int) calc.size())-1; i >= 0; i--){
			auto &ca = calc[i];
			
			EqItem rep; rep.type = NOOP;
			double rep_con = UNSET;

			switch(ca.op){
			case MULTIPLY:
				{
					// Combines together constant values
					{
						auto j = 0u; while(j < ca.item.size() && ca.item[j].type != NUMERIC) j++;
						if(j < ca.item.size()){
							auto jj = j+1;
							while(jj < ca.item.size()){
								if(ca.item[jj].type == NUMERIC){
									ca.item[j].constant = mult_const(ca.item[j],ca.item[jj]);
									ca.item.erase(ca.item.begin()+jj);
									flag = true;
								}
								else jj++;
							}
						}
					}
					
					{
						// Removes anything with one
						auto j = 0u;
						while(j < ca.item.size()){
							if(ca.item[j].type == NUMERIC && ca.item[j].constant == 1){
								ca.item.erase(ca.item.begin()+j);
								flag = true;
							}
							else j++;
						}
					}
					
					// If contains zero then entire sum is zero
					auto j = 0u;
					while(j < ca.item.size() && !(ca.item[j].type == NUMERIC && ca.item[j].constant == 0)) j++;
					if(j < ca.item.size()) rep_con = 0;
					else{
						if(ca.item.size() == 0) rep_con = 1;
						else{
							if(ca.item.size() == 1) rep = ca.item[0];
						}
					}
				}
				break;
				
			case ADD:
				{			
					// Combines together constant values
					auto j = 0u; while(j < ca.item.size() && ca.item[j].type != NUMERIC) j++;
					if(j < ca.item.size()){
						auto jj = j+1;
						while(jj < ca.item.size()){
							if(ca.item[jj].type == NUMERIC){
								ca.item[j].constant = add_const(ca.item[j],ca.item[jj]);
								ca.item.erase(ca.item.begin()+jj);
							}
							else jj++;
						}
					}
					
					// Removes anything with zero
					j = 0;
					while(j < ca.item.size()){
						if(ca.item[j].type == NUMERIC && ca.item[j].constant == 0){
							ca.item.erase(ca.item.begin()+j);
						}
						else j++;
					}
					
					if(ca.item.size() == 0) rep_con = 0;
					if(ca.item.size() == 1) rep = ca.item[0];
				}
				break;
				
			case TAKE:
				break;
				
			case DIVIDE:
				{
					double con1 = UNSET, con2 = UNSET;
					if(ca.item[0].type == NUMERIC) con1 = ca.item[0].constant;
					if(ca.item[1].type == NUMERIC) con2 = ca.item[1].constant;
				
					if(con2 == 0){
						warn = "Equation contains a division by zero";
					}
					else{
						if(con1 == UNDEF || con2 == UNDEF) rep_con = UNDEF;
						else{
							if(con1 != UNSET){
								if(con1 == 0) rep_con = 0;
								else{
									if(con1 == INFY) rep_con = INFY;
									else{
										if(con2 != UNSET){
											if(con2 == INFY) rep_con = 0;
											else rep_con = con1/con2;
										}
									}
								}								
							}
						}
					}
				}
				break;
				
			case EXPFUNC:
				if(ca.item[0].type == NUMERIC){
					auto con = ca.item[0].constant;
					if(con == INFY) rep_con = INFY;
					else{
						if(con == UNDEF) rep_con = UNDEF;
						else rep_con = exp(con);
					}
				}					
				break;
			
			case SINFUNC:
				if(ca.item[0].type == NUMERIC){
					auto con = ca.item[0].constant;
					if(con == INFY) rep_con = UNDEF;
					else{
						if(con == UNDEF) rep_con = UNDEF;
						else rep_con = sin(con);
					}
				}					
				break;
				
			case COSFUNC:
				if(ca.item[0].type == NUMERIC){
					auto con = ca.item[0].constant;
					if(con == INFY) rep_con = UNDEF;
					else{
						if(con == UNDEF) rep_con = UNDEF;
						else rep_con = cos(con);
					}	
				}					
				break;
			
			case LOGFUNC:
				if(ca.item[0].type == NUMERIC){
					auto con = ca.item[0].constant;
					if(con <= 0){
						warn = "Equation contains the logarithm of a non-positive";
					}
					else{
						if(con == INFY) rep_con = INFY;
						else{
							if(con == UNDEF) rep_con = UNDEF;
							else rep_con = log(con);	
						}
					}
				}
				break;
				
			case POWERFUNC:
				{
					double con1 = UNSET, con2 = UNSET;
					if(ca.item[0].type == NUMERIC) con1 = ca.item[0].constant;
					if(ca.item[1].type == NUMERIC) con2 = ca.item[1].constant;
				
					if(con1 == UNDEF || con2 == UNDEF) rep_con = UNDEF;
					else{
						if(con2 == 0) rep_con = 1;
						else{
							if(con1 != UNSET){
								if(con1 == INFY || con2 == INFY) rep_con = INFY;
								else{	
									if(con1 == 0) rep_con = 0;
									else{
										if(con2 != UNSET) rep_con = pow(con1,con2);
									}
								}
							}
						}
					}
				}				
				break;
		
			case THRESHFUNC:
				if(ca.item[0].type == NUMERIC && ca.item[1].type == NUMERIC){
					auto con1 = ca.item[0].constant, con2 = ca.item[1].constant;
					if(con1 == UNDEF || con2 == UNDEF) rep_con = UNDEF;
					else{
						if(con1 < con2) rep_con = 0;
						else rep_con = con1;
					}
				}
				break;
				
			case UBOUNDFUNC:
				if(ca.item[0].type == NUMERIC && ca.item[1].type == NUMERIC){
					auto con1 = ca.item[0].constant, con2 = ca.item[1].constant;
					if(con1 == UNDEF || con2 == UNDEF) rep_con = UNDEF;
					else{
						if(con1 > con2) rep_con = INFY;
						else rep_con = con1;
					}
				}
				break;
			
			case STEPFUNC:
				if(ca.item[0].type == NUMERIC){
					auto con = ca.item[0].constant;
					if(con == UNDEF) rep_con = UNDEF;
					else{
						if(con > 0) rep_con = 1; else rep_con = 0;
					}
				}
				break;
				
			case MAXFUNC:
				if(ca.item[0].type == NUMERIC && ca.item[1].type == NUMERIC){
					auto con1 = ca.item[0].constant, con2 = ca.item[1].constant;
					if(con1 == UNDEF || con2 == UNDEF) rep_con = UNDEF;
					else{
						if(con1 > con2) rep_con = con1;
						else rep_con = con2;
					}
				}
				break;
				
			case MINFUNC:
				if(ca.item[0].type == NUMERIC && ca.item[1].type == NUMERIC){
					auto con1 = ca.item[0].constant, con2 = ca.item[1].constant;
					if(con1 == UNDEF || con2 == UNDEF) rep_con = UNDEF;
					else{
						if(con1 < con2) rep_con = con1;
						else rep_con = con2;
					}
				}
				break;
				
			case ABSFUNC:
				if(ca.item[0].type == NUMERIC){
					auto con = ca.item[0].constant;
					if(con == UNDEF) rep_con = UNDEF;
					else{
						if(con < 0) rep_con = -con;
						else rep_con = con;
					}
				}
				break;
			
			case SQRTFUNC:
				if(ca.item[0].type == NUMERIC){
					auto con = ca.item[0].constant;
					if(con < 0){
						warn = "Equation contains the square root of a non-positive";
					}
					else{
						if(con == UNDEF) rep_con = UNDEF;
						else{
							if(con == INFY) rep_con = INFY;
							else rep_con = sqrt(con);
						}
					}
				}
				break;

			default: emsg("Should be simplification"); break;
			}
			
			if(rep.type != NOOP || rep_con != UNSET){
				flag = true;
				if(rep.type != NOOP && rep_con != UNSET) emsg("Cannot be both here");
				
				if(rep_con != UNSET){
					rep.type = NUMERIC; rep.constant = rep_con; rep.index = UNSET; rep.num = UNSET;
				}
				
				replace_reg(ca.reg_store,rep);
				calc.erase(calc.begin()+i);
			}
		}
		
		loop++; if(loop > 100) emsg("Simplify equation problem");
	}while(flag == true);
	
	remove_unused();
	
	// Checks that INFY and UNDEF do to exist in the final equation
	for(auto i = 0u; i < calc.size(); i++){
		const auto &ca = calc[i];
		for(auto j = 0u; j < ca.item.size(); j++){
			const auto &it = ca.item[j];
			if(it.type == NUMERIC){
				if(it.constant == UNDEF) warn = "The equation contains an undefined quantity";
				if(it.constant == INFY) warn = "The equation contains an infinite quantity";
			}
		}
	}
	
	if(ans.type == NUMERIC){
		if(ans.constant == UNDEF) warn = "The equation contains an undefined quantity";
		if(ans.constant == INFY) warn = "The equation contains an infinite quantity";
	}
	
	if(false) print_calculation();
}


/// Gets a list showing how many times a register is used in calculation
vector < vector <unsigned int> > Equation::get_reg_used(const vector <Calculation> &calc) const
{
	vector < vector <unsigned int> > reg_used(calc.size());
	for(auto ii = 0u; ii < calc.size(); ii++){
		const auto &item = calc[ii].item;
		for(auto j = 0u; j < item.size(); j++){
			if(item[j].type == REG){ reg_used[item[j].num].push_back(ii);}
		}
	}

	return reg_used;
}


/// Checks that reg_used is still correctly specified
void Equation::check_reg_used(const vector <Calculation> &calc, vector < vector <unsigned int> > &reg_used) const
{
	auto reg_used_check = get_reg_used(calc);
	if(reg_used.size() != reg_used_check.size()) emsg("reg used size wrong");

	for(auto i = 0u; i < reg_used.size(); i++){
		sort(reg_used[i].begin(),reg_used[i].end());
		sort(reg_used_check[i].begin(),reg_used_check[i].end());
		if(!equal_vec(reg_used[i],reg_used_check[i])){
			for(auto va : reg_used[i]) cout << va << " "; 
			cout << " va"<< endl;
			for(auto va : reg_used_check[i]) cout << va << " "; 
			cout << " va"<< endl;
			
			emsg("Problem with get_reg");
		}
	}
}


/// Removes any unused registers
void Equation::remove_unused()
{
	if(false) print_calculation();
	
	if(ans.type != REG){
		nreg = 0;
		calc.resize(0);
		return;
	}
	
	vector < vector <unsigned int> > reg_dep;
	reg_dep.resize(nreg);
	for(auto j = 0u; j < calc.size(); j++){
		auto r = calc[j].reg_store;
		reg_dep[r].push_back(j);
		for(const auto &it : calc[j].item){
			if(it.type == REG){
				for(auto rr : reg_dep[it.num]) reg_dep[r].push_back(rr);
			}
		}
	}
	
	if(false){
		for(auto r = 0u; r < nreg; r++){
			cout << r << ": ";
			for(auto k : reg_dep[r]) cout << k << ",";
			cout << "reg dep" << endl;
		}
	}
	
	vector <bool> map(nreg,false);
	for(auto r : reg_dep[ans.num]) map[r] = true;
	
	EqItem rep; rep.type = NOOP;
	for(int r = nreg-1; r >= 0; r--){
		if(map[r] == false){
			if(calc[r].reg_store != (unsigned int)r) emsg("Reg_store needs to be searched");
			replace_reg(r,rep);
			calc.erase(calc.begin()+r);	
		}
	}
}

 
/// Adds constant values
double Equation::add_const(EqItem item1, EqItem item2) const
{
	auto val1 = item1.constant, val2 = item2.constant;
	if(val1 == UNDEF || val2 == UNDEF) return UNDEF; 

	if(val1 == INFY || val2 == INFY) return INFY;
	return val1+val2;
}


/// Multiplies constant values
double Equation::mult_const(EqItem item1, EqItem item2) const
{
	auto val1 = item1.constant, val2 = item2.constant;

	if(val1 == UNDEF || val2 == UNDEF) return UNDEF; 
	if(val1 == INFY){
		if(val2 == 0) return UNDEF;
		return INFY;
	}
	else{
		if(val2 == INFY){
			if(val1 == 0) return UNDEF;
			return INFY;
		}
	}
	
	return val1*val2;
}


/// Calculates a list of all population used 
void Equation::calculate_pop_ref()
{				
	for(auto i = 0u; i < calc.size(); i++){
		const auto &ca = calc[i];
		for(const auto &it : ca.item){
			if(it.type == POPNUM) add_to_vec(pop_ref,it.num);
		}
	}
	
	if(ans.type == POPNUM) add_to_vec(pop_ref,ans.num); 
}		


/// Determines if the equation is zero
bool Equation::is_zero() const
{
	if(calc.size() == 0 && ans.type == NUMERIC && ans.constant == 0) return true;
	return false;
}


/// Determines if the equation is one
bool Equation::is_one() const
{
	if(calc.size() == 0 && ans.type == NUMERIC && ans.constant == 1) return true;
	return false;
}


/// Determines whether equation is time varying or not
void Equation::set_time_vari()
{
	time_vari = false;
	if(pop_ref.size() > 0) time_vari = true;
	for(const auto &pr : param_ref){
		if(param[pr.th].time_dep == true) time_vari = true;
	}
	
  for(const auto &ca : calc){
		for(const auto &it : ca.item){
			if(it.type == TIME) time_vari = true;
		}
		if(ans.type == TIME) time_vari = true;
	}
}


/// Calculates the individual factor for the equation
double Equation::indfac(const Individual &ind) const
{
	auto num = 1.0; 
	for(auto ie : ind_eff_mult) num *= ind.exp_ie[ie];
	for(auto fe : fix_eff_mult) num *= ind.exp_fe[fe];
	
	return num;
}


/// Determines if equation in linear in a given parameter
bool Equation::param_linear(unsigned int th) const
{
	switch(ans.type){
	case PARAMETER:
		if(ans.num == th) return true;
		break;
		
	case REG:
		{
			const auto &last = calc[calc.size()-1];
			if(last.reg_store != ans.num) emsg("Reg does not agree");
			
			switch(last.op){
			case MULTIPLY:
				for(const auto &it : last.item){
					if(it.type == PARAMETER && it.num == th) return true;
				}
				break;
				
			case ADD:
				if(last.item.size() == 1){
					const auto &it = last.item[0];
					if(it.type == PARAMETER && it.num == th) return true;
				}
				break;
			
			case DIVIDE:
				{
					const auto &it = last.item[0];
					if(it.type == PARAMETER && it.num == th) return true;
				}
				break;
			
			default: break;
			}
		}
		break;
		
	default: break;
	}
	
	return false; 
}
	

/// Sets any references to parameters used in equations
void Equation::setup_param_ref()
{
	for(auto i = 0u; i < calc.size(); i++){
		const auto &ca = calc[i];
		for(const auto &it : ca.item){
			switch(it.type){
			case PARAMETER: case SPLINE:
				add_to_vec(param_ref,it.num,it.index);
				break;
			default: break;
			}
		}
	}
	
	switch(ans.type){
	case PARAMETER: case SPLINE:
		add_to_vec(param_ref,ans.num,ans.index); 
		break;
		default: break;
	}
}


/// Sets up references (pop_ref, param_ref) in equations
void Equation::setup_references()
{
	setup_param_ref();
	calculate_pop_ref();       
	if(infection_trans) setup_comp_pref_convert();
	set_time_vari();
	
	if(false){
		print_calculation();
		cout << param_ref.size() << " " << pop_ref.size() << " param pop" << endl;
	}
}


/// Finds the distance between two compartments
double Equation::find_dist(unsigned int c, unsigned int cc, const vector <Compartment> &comp, Coord coord) const
{
	if(c == cc) return 0;

	switch(coord){
	case CARTESIAN:
		{
			auto x = comp[c].x, y = comp[c].y;
			auto xx = comp[cc].x, yy = comp[cc].y;

			auto dx = x-xx, dy = y-yy;
			return sqrt(dx*dx + dy*dy);
		}
		break;
		
	case LATLNG:
		{
			auto lat = comp[c].lat*MM_PI/180.0;
			auto lng = comp[c].lng*MM_PI/180.0;

			auto lat2 = comp[cc].lat*MM_PI/180.0;
			auto lng2 = comp[cc].lng*MM_PI/180.0;

			auto r = 6371.0;

			auto dlng = lng2-lng;
			if(dlng > MM_PI) dlng -= 2*MM_PI;
			if(dlng < -MM_PI) dlng += 2*MM_PI;
	
			auto si = sin(0.5*(lat2-lat));
			auto si2 = sin(0.5*dlng);
			return 2*r*asin(sqrt(si*si + cos(lat)*cos(lat2)*si2*si2));
		}
		break;
	}

	emsg("Cannot find distance");
	return UNSET;
}

				
/// Gets the distance between two compartments
double Equation::get_distance(const ParamProp &pp)
{
	if(pp.dep.size() != 2){	
		warn = "The distance matrix '"+pp.name+"' must have two indices."; 
		return UNSET;
	}
	
	auto dep1 = pp.dep[0];
	auto dep2 = pp.dep[1];

	for(auto p = 0u; p < species.size(); p++){
		const auto &sp = species[p];
		for(auto cl = 0u; cl < sp.cla.size(); cl++){
			const auto &claa = sp.cla[cl];
			
			//auto vec = claa.hash_comp.get_vec_string(dep1);
			//auto c = claa.hash_comp.existing(vec);
			auto c = claa.hash_comp.find(dep1);
			if(c != UNSET){
				//auto vec = claa.hash_comp.get_vec_string(dep2);
				//auto cc = claa.hash_comp.existing(vec);
				auto cc = claa.hash_comp.find(dep2);
				if(cc == UNSET){
					warn = "Compartment '"+dep2+"' is not recognised in distance matrix.";
					return UNSET;
				}
				else{
					return find_dist(c,cc,claa.comp,claa.coord); 
				}
			}
		}
	}
	
	warn = "Compartment '"+dep1+"' is not recognised in distance matrix.";
	return UNSET;
}


/// Provides a name for an operator
string Equation::op_name(EqItemType type) const 
{
	switch(type){
	case ADD: return "addition '+'";
	case TAKE: return "subtraction '-'";
	case MULTIPLY: return "multiplication '×'";
	case DIVIDE: return "division '/'"; 
	default: break;
	}
	return "";
}


/// Checks for any repeated operators in expression
void Equation::check_repeated_operator(const vector <EqItem> &op)
{
	if(op.size() == 1){
		switch(op[0].type){
		case ADD: case TAKE: case MULTIPLY: case DIVIDE:
			warn = "Expression cannot be simply "+op_name(op[0].type)+".";
			return;
		default: break;
		}
	}
	
	if(op.size() > 0){
		for(auto i = 0u; i < op.size()-1; i++){
			switch(op[i].type){
			case ADD: case TAKE: case MULTIPLY: case DIVIDE:
				switch(op[i+1].type){
				case ADD: case TAKE: case MULTIPLY: case DIVIDE:
					warn = "Problem with the expression. Cannot have operators "+op_name(op[i].type)+" and "+op_name(op[i+1].type)+" next to each other";
					
					break;
				default: break;
				}
			default: break;
			}
		}
	}
}


/// Replaces minus sign with plus one multiply
void Equation::replace_minus(vector <EqItem> &op)
{
	auto i = 0u;
	while(i < op.size()){
		if(op[i].type == TAKE){
			op[i].type = ADD;
			
			EqItem item2; item2.type = NUMERIC; item2.constant = -1;
			op.insert(op.begin()+i+1,item2);
	
			EqItem item3; item3.type = MULTIPLY;
			op.insert(op.begin()+i+2,item3);
		}
		else i++;
	}
}

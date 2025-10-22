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
Equation::Equation(string tex, EqnType ty, unsigned int p, unsigned int cl, bool inf_trans, unsigned int tif, unsigned int li_num, const vector <SpeciesSimp> &species, vector <Param> &param, vector <Prior> &prior, const vector <Derive> &derive, const vector <Spline> &spline, const vector <ParamVecEle> &param_vec, vector <Population> &pop, Hash &hash_pop, Constant &constant, const vector <double> &timepoint, const Details &details) : species(species), param(param), prior(prior), derive(derive), spline(spline), param_vec(param_vec), pop(pop), hash_pop(hash_pop), constant(constant), timepoint(timepoint), details(details)
{
	plfl = false;  // Set to true to print operations to terminal (used for diagnostics)

	precalc_done = false;
	stop_combine_fl = false;
	contain_population = false;

	if(plfl) cout << endl << tex << "  start equation" << endl;

	nspecies = species.size();
	type = ty; sp_p = p; sp_cl = cl; ti_fix = tif, line_num = li_num; warn = "";
	infection_trans = inf_trans;
	markov_eqn_ref = UNSET;
	//ans.type = NUMERIC; ans.num = UNSET;
	
	tex = trim(tex); if(tex == ""){ warn = "There is no equation"; return;}

	te = tex;
	te_raw = te;
	te_raw = replace(te_raw,"%","");
	te_raw = replace(te_raw,"$","");
	
	if(warn != "") return;

	unravel_sum(); 

	if(warn != "") return;

	te = replace(te,"×","*");                        // Converts × to *
	
	auto op = extract_operations();                  // Extracts the operations in the 	expression
	
	if(warn != "") return; 
	
	check_repeated_operator(op);                     // Checks for repeated operators e.g. "**"
	
	//if(integral.size() > 0){ print_operations(op); emsg("op");}

	if(warn != "") return; 
	
	replace_minus(op);                               // Replaces minums sign with (-1)*
		
	if(plfl) print_operations(op);

	time_integral(op);                            // Incorporates time integral

	if(warn != "") return; 
	
	calcu = create_calculation(op);               // Works out the sequence of calculation to generate result
	
	if(warn != "") return; 
	
	if(plfl == true) print_calculation();

	if(simplify_eqn == true){
		simplify(calcu);     // Simplifies by combining constants
	}
	
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
	te.clear();
	if(!debugging){
		//if(type == REPARAM || type == REPARAM_EQN) te_raw = trunc(te_raw,20);
		//else te_raw = trunc(te_raw,20);	
		te_raw = trunc(te_raw,20);	
	}
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
	if(calcu.size() == 0){
		warn = "No equation";
		return;
	}
	
	for(auto i = 0u; i < calcu.size(); i++){
		const auto &ca = calcu[i];
		if(i+1 != calcu.size() && ca.item.size() == 0){
			warn = "Calculation item with no elements";
			return;
		}
	}
}


/// Checks that operators are correctly specified
void Equation::check_opl(const vector <EqItemList> &opl) const 
{
	auto op = create_op_from_opl(opl);
	for(auto i = 0u; i < op.size(); i++){
		if(op[i].type == POWERFUNC){
			if(i+5 >= op.size()) emsg_input("op prob3");
			if(op[i+1].type != LEFTBRACKET) emsg_input("op prob1");
			if(op[i+3].type != FUNCDIVIDE) emsg_input("op prob2");
			if(op[i+5].type != RIGHTBRACKET) emsg_input("op prob2");
		}
	}
}


/// Gets all compartments consistent with a string
// 'I'  all infectious individuls
// 'I,M' all infectious individuals
// 'E|I,M' alld infectous or exposed individuals
vector <unsigned int> Equation::get_all_comp(unsigned int p, string te)
{
	if(toLower(trim(te)) == "all") te = "";
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
	auto imax = op.size();
	if(imax > 100) imax = 100;
  for(auto i = 0u; i < imax; i++){
    switch(op[i].type){
			case TINT: cout << "\\int dt"; break;
      case LEFTBRACKET: cout << "("; break;
      case RIGHTBRACKET: cout << ")"; break;
			case FUNCDIVIDE: cout << "|"; break;
			case INTEGRAL: 
				{
					const auto &inte = integral[op[i].num];
					cout << "\\int dt[" << inte.ti_min << "," << inte.ti_max << "]";
				}
				break;
				
      case PARAMETER: 
				{
					const auto &pr = param_ref[op[i].num];
					cout << add_escape_char(param[pr.th].name) << pr.index;
				}
				break;
				
			case PARAMVEC:
				emsg_input("Should not be param vec"); 
				break;
				
			case SPLINE:
				{
					const auto &pr = param_ref[op[i].num];
					cout << "Spline " << add_escape_char(param[pr.th].name) << pr.index;
				}
				break;
			
			case DERIVE: 
				{
					const auto &dr = derive_ref[op[i].num];
					cout << "Derive(" << derive[dr.i].name << dr.index << ")";
					if(dr.ti != UNSET) cout << "(ti=" << dr.ti << ")";
				}
				break;
				
			case SPLINEREF: case CONSTSPLINEREF: emsg_input("spline ref should not"); break;
			case IE: cout << species[sp_p].ind_effect[op[i].num].name; break;
			case ONE: cout << "1"; break;
			case ZERO: cout << "0"; break;
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
      case SIGFUNC: cout << "sigmoid"; break;
			case ADD: cout << "+"; break;
      case TAKE: cout << "-"; break;
      case MULTIPLY: cout << "*"; break;
			case SINGLE: break;
			case DIVIDE: cout << "/"; break;
      case REG: cout << "R" << op[i].num; break;
      case NUMERIC: cout << "numeric" << constant.value[op[i].num]; break;
			case TIME: cout << "time"; break;
			case REG_FAC: cout << "Rfac" << op[i].num; break;
			case REG_PRECALC: cout << "Rpre" << op[i].num; break;
			case REG_PRECALC_TIME: cout << "Rpretime" << op[i].num; break;
			case NOOP: cout << "No operation"; break;
		}
    cout << " ";
  }
	if(imax < op.size()) cout << "...";
  cout << endl;
}

	
/// Prints steps used for a calculation
void Equation::print_calculation() const   
{
	cout << "For equation '" << add_escape_char(te_raw) << "' calculation:" << endl;
	auto imax = calcu.size();
	if(imax > 10) imax = 10;
  for(auto i = 0u; i < imax; i++){
		const auto &ca = calcu[i];
		print_ca(i,ca);
		cout << endl;
  }
	if(imax < calcu.size()) cout << "..." << endl;

  cout << endl << endl;
	
	for(auto k = 0u; k < ind_eff_mult.size(); k++){
		cout << "multiplied by IE " << species[sp_p].ind_effect[ind_eff_mult[k]].name << endl;
	}
	
	for(auto k = 0u; k < fix_eff_mult.size(); k++){
		cout << "multiplied by FE " << species[sp_p].fix_effect[fix_eff_mult[k]].name << endl;
	}
}


/// Prints the calculation from a specific calculation
void Equation::print_ca(unsigned int i, const Calculation &ca) const
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
	case SIGFUNC: cout << "sigmoid("; break;
	case ADD: break;
	case TAKE: break;
	case MULTIPLY: break;
	case SINGLE: break;
	case DIVIDE: break;
	default: emsg_input("Eq problem1"); break;
	}
		
	const auto &item = ca.item;
		
	auto jmax = item.size();
	if(jmax > 10) jmax = 10;
	for(auto j = 0u; j < jmax; j++){
		print_item(item[j]);
		
		if(j != item.size()-1){
			switch(ca.op){
			case ADD: cout <<  "+"; break;
			case TAKE: cout <<  "-"; break;
			case MULTIPLY: cout << "*"; break;
			case SINGLE: cout <<  ""; break;
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
	if(jmax < item.size()) cout << "...";

	switch(ca.op){
	case EXPFUNC: case SINFUNC: case COSFUNC: case LOGFUNC: case STEPFUNC: case POWERFUNC: 
	case THRESHFUNC: case UBOUNDFUNC: case MAXFUNC: case MINFUNC: case ABSFUNC: case SQRTFUNC:
	case SIGFUNC:
		cout << ")"; 
		break;
	default: break;
	}
		
	cout <<  " > ";
	cout <<  "R" << i;
}
 
 
/// Prints a multiplication
void Equation::print_mult(string name, const vector <EqItem> &item) const
{  
	cout << name << ": ";
	for(auto j = 0u; j < item.size(); j++){
		if(j != 0) cout << "*";
		print_item(item[j]);
	}	
	cout << endl;
}

 
 /// Prints a simgle item
void Equation::print_item(const EqItem &it) const
{
	switch(it.type){
	case INTEGRAL:
		{
			const auto &inte = integral[it.num];
			cout << "\\int dt[" << inte.ti_min << "," << inte.ti_max << "]";
		}
		break;
		
	case PARAMETER: 
		{
			const auto &pr = param_ref[it.num];
			cout << add_escape_char(param[pr.th].name) << pr.index;
		}
		break;
		
	case PARAMVEC: 
		cout << add_escape_char(param_vec[it.num].name);
		break;
	
	case SPLINE: 
		{
			const auto &pr = param_ref[it.num];
			auto par = param[pr.th]; 
			cout << "Spline " << add_escape_char(get_param_name_with_dep(par,par.dep,pr.index));
		}
		break;
		
	case SPLINEREF:
		cout << "Spline " << add_escape_char(spline[it.num].name); 
		break;
		
	case CONSTSPLINEREF:
		cout << "Const Spline " << add_escape_char(spline[it.num].name); 
		break;
	
	case DERIVE: 
		{
			const auto &dr = derive_ref[it.num];
			cout << "Derive(" << derive[dr.i].name << dr.index << ")";
			if(dr.ti != UNSET) cout << "(ti=" << dr.ti << ")";
		}
		break;
			
	case POPNUM: 
		//cout << "P" << it.num;
		cout << "'" << pop[it.num].name << "'";
		break;

	case IE: cout << species[sp_p].ind_effect[it.num].name; break;
	case ONE: cout << "1"; break;
	case ZERO: cout << "0"; break;
	case FE: cout << species[sp_p].fix_effect[it.num].name; break;
	case REG: cout << "R" << it.num; break;
	case REG_FAC: cout << "Rfac" << it.num; break;
	case REG_PRECALC: cout << "Rpre" << it.num; break;
	case REG_PRECALC_TIME: cout << "Rpretime" << it.num; break;
	case NUMERIC: 
		{
			auto val = constant.value[it.num];
			if(val == INFY) cout << "INFY";
			else{
				if(val == UNDEF) cout << "UNDEF";
				else cout << val; 
			}
		}
		break;
	case TIME: cout << "time"; break;
	default: emsg_input("Eq problem2a"); break;
	}
}

		
/// Unravels all the sums in the equation
void Equation::unravel_sum() 
{
	auto fl = true;
	while(fl){ // Iterates to do nested sums
		fl = false;
		string te_new = "";
		
		auto i = 0u;

		while(i < te.length()){
			auto ist = i;
			while(i < te.length() && !str_eq(te,i,sigma)) i++;
			if(i < te.length()){
				fl = true;
				te_new += te.substr(ist,i-ist);
				
				auto istart = i;

				i += sigma.length();
				if(i == te.length()){ warn = "The character '_' must follow the sum 'Σ'"; return;}

				auto ibra = i;
				while(ibra < te.length() && te.substr(ibra,1) != "(") ibra++;
		
				if(ibra == te.length()){
					warn = "Terms in a sum 'Σ' must be enclosed by left '(' and right ')' brackets."; 
					return;
				}
				
				auto tex = te.substr(istart,ibra-istart);
				
				double distmax = UNSET;
				string comp_max;
				
				if(te.substr(i,1) == "^"){
					warn = "The sum '"+tex+"' cannot contain a superscript '^'"; return;
				}
				
				if(te.substr(i,1) != "_"){ warn = "The character '_' must follow the sum 'Σ'"; return;}
				i++;
				
				auto di = get_dep_info(te,i,notparam_list);
				if(di.warn != ""){
					warn = "The sum 'Σ' has misspecified indices: "+di.warn;
					return;
				}

				auto dist = trim(te.substr(di.iend,ibra-(di.iend)));
				if(dist != ""){
					auto spl = split(dist.substr(1,dist.length()-2),',');
					if(spl.size() != 2) emsg_input("Split prob");
					
					comp_max = trim(spl[0]);			
					distmax = number(spl[1]);
				}
				
				i = ibra+1;
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
			
				auto content = te.substr(ist,i-ist);
		
				vector <CompPos> comp_pos;
	 
				auto ndep = di.spl.size();

				for(auto j = 0u; j < ndep; j++){
					auto ind = di.spl[j];
					auto ind2 = remove_prime(ind);
					
					auto list = find_list_from_index(ind2,distmax,comp_max);
		
					if(list.size() == 0){ warn = "The index '"+ind2+"' is not found within the model"; return;}
					
					CompPos cp; cp.list = list; cp.index = 0;
					comp_pos.push_back(cp);
				}

				te_new += "(";

				auto first = true;
				
				vector <DepConv> dep_conv;
				for(auto j = 0u; j < ndep; j++){
					DepConv dc; dc.before = di.spl[j];
					dep_conv.push_back(dc);
				}
				
				auto swap_temp = swap_template(content,dep_conv);
				if(swap_temp.warn != ""){ warn = swap_temp.warn; return;}
		
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
							emsg_input("Swap index dif res");
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

					if(first == false) te_new += "+";
					te_new += "("+cont_new+")";
					first = false;
				}while(fl == false);
			}	
			else{
				te_new += te.substr(ist);
			}
		}
		
		if(fl){
			te = te_new;
		}
	}
}


/// Finds the list of compartments from a given index
vector <string> Equation::find_list_from_index(string ind, double dist_max, string comp_max) const
{
	vector <string> vec;

	if(dist_max != UNSET){
		if(comp_max == ""){
			emsg_input("Error with 'max' function. The compartment isn't specified");
		}
		
		for(auto p = 0u; p < species.size(); p++){
			const auto &sp = species[p];
			for(auto cl = 0u; cl < sp.cla.size(); cl++){
				const auto &claa = sp.cla[cl];
				if(claa.index == ind){
					const auto &comp = claa.comp;
					
					unsigned int c;
					for(c = 0u; c < comp.size(); c++){
						if(comp[c].name == comp_max) break;
					}
					if(c == comp.size()) emsg_input("Could not find compartment '"+comp_max+"'");
					
					for(auto cc = 0u; cc < comp.size(); cc++){
						auto d = find_dist(c,cc,comp,claa.coord); 
						if(d < dist_max) vec.push_back(comp[cc].name);
					}
					return vec;
				}
			}	
		}
		emsg_input("Could not find");
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
    while(i < te.length()){
			num = te.at(i)-48;
      if(num >= 0 && num <= 9){ num2 += fac*num; fac /= 10;}
      else break;
      i++;
    }
  }

  raend = i;
  return sign*num2;
}


/// Tries to get a parameter name from a string
ParamRef Equation::get_param_name(unsigned int i, double &dist, unsigned int &raend)
{
	ParamRef pref; 
	
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
					if(name == iden_matrix_name || name == iden_matrix_name2){
						dist = get_identity(pp);
					}
					else{
						warn = "Could not find parameter '"+name+"'"; return pref;
					}
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
							return pref;
						}
						ind += j*par.dep[de].mult;
					}
				}
				
				pref.th = th; pref.index = ind;
			}
		}
	
		raend = i+1;
	}
	
  return pref;
}


/// Tries to get a parameter name from a string
DeriveRef Equation::get_derive_name(unsigned int i, unsigned int &raend)
{
	DeriveRef dref; 
	
	if(te.substr(i,1) == "%"){
		auto ist = i;
		while(i < te.length() && te.substr(i,1) != "$") i++;
	
		if(i == te.length()){ warn = "Could not find right bracket '$'"; return dref;}
		
		auto content = trim(te.substr(ist+1,i-ist-1));
	
		auto pp = get_param_prop(content);
		
		auto name = pp.name;
	
		auto d = 0u; while(d < derive.size() && derive[d].name != name) d++;	
		if(d < derive.size()){
			const auto &der = derive[d];
		
			if(pp.time_dep != der.time_dep){ 
				warn = "The time dependency in the equation does not agree with the definition for parameter '"+der.name+"'."; 
				return dref;
			}
			
			auto dep = pp.dep;
			if(dep.size() > 0 && dep[dep.size()-1] == "t") dep.pop_back();
			
			if(dep.size() != der.dep.size()){
				warn = "The dependency in the equation does not agree with the definition for parameter '"+der.name+"'."; 
				return dref;
			}
			
			if(pp.prime == true){
				warn = "\""+content+"\" should not have a prime (perhaps a sum is missing?)."; 
				return dref;
			}
			
			auto ind = 0u;
			for(auto de = 0u; de < der.dep.size(); de++){
				auto j = der.dep[de].hash_list.find(remove_prime(pp.dep[de]));
				if(j == UNSET){
					warn = "'"+content+"' does not agree with the definition '"+der.name+"'."; 
					return dref;
				}
				ind += j*der.dep[de].mult;
			}
			dref.i = d; dref.index = ind;
		}
	
		raend = i+1;
	}
	
  return dref;
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
			hash_pop.add(p,vec);
			
			string start = "In population {"+cont+"}: ";
	
			Population po;
			po.name = name;
			
			auto spl = split(cont,';');
			
			if(spl.size() > 2){
				warn = start+"More than one ';' unexpected"; return p;
			}
			
			if(spl.size() == 2){
				auto extra = spl[1];
			
				cont = spl[0];	
				
				auto k = 0u;
				while(k < extra.length()){
					while(k < extra.length() && extra.substr(k,1) == " ") k++;
					if(k == extra.length()) break;
					
					auto ch = extra.substr(k,1);
					if(ch == "[" || ch == "<"){
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
								warn = "Fixed effect '<"+fe_name+">' unspecified. This must be specified through the 'fixed-effect' command"; return p;
							}
							po.fix_eff_mult.push_back(fe);
						}
					}
					else{
						if(ch != "×" && ch != "*"){
							warn = start+"Character '"+ch+"' is unexpected"; return p;
						}
					}
					k++;
				}
			}
	
			// Checks that population individual/fixed effects are divided correctly
			for(auto k = 0u; k < cont.length(); k++){
				auto ch = cont.substr(k,1);
				if(ch == "<"){
					warn = start+"Population fixed effects must be placed after ';' divider";
					return p;
				}
				
				if(ch == "["){
					warn = start+"Population individual effects must be placed after ';' divider";
					return p;
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
			
			contain_population = true;
			
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
		
		if(type == DERIVE_EQN){
			warn = "Derived equations cannot contain individual effects (except in populations)"; return ie;
		}
		
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

		if(type == DERIVE_EQN){
			warn = "Derived equations cannot contain fixed effects (except in populations)"; return fe;
		}
		
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
	case PARAMETER: case PARAMVEC:
	case SPLINE: case SPLINEREF: case CONSTSPLINEREF:
	case DERIVE:
	case INTEGRAL:
	case POPNUM: case IE: case FE: case REG: case NUMERIC: case TIME: 
		return true;
  default: return false;
	}
}


/// Determines if a quantity exists
bool Equation::quantl(const vector <EqItemList> &opl, unsigned int i) const
{
  if(i == UNSET_LIST) return false;
  switch(opl[i].type){
	case PARAMETER: case PARAMVEC:
	case SPLINE: case SPLINEREF: case CONSTSPLINEREF:
	case DERIVE:
	case INTEGRAL:
	case POPNUM: case IE: case FE: case REG: case NUMERIC: case TIME: 
		return true;
  default: return false;
	}
}


/// Determines if a function exists
bool Equation::is_func(const vector <EqItem> &op, int i) const
{
  if(i < 0 || i >= (int)op.size()) return false;
  switch(op[i].type){
	case EXPFUNC: case SINFUNC: case COSFUNC: case LOGFUNC: case POWERFUNC: 
	case THRESHFUNC: case UBOUNDFUNC: case STEPFUNC: case MAXFUNC: case MINFUNC:
	case ABSFUNC: case SQRTFUNC: case SIGFUNC:
		return true;
	default: return false;
	}
}


/// Determines if a function exists
bool Equation::is_funcl(const vector <EqItemList> &opl, unsigned int i) const
{
  if(i == UNSET_LIST) return false;
  switch(opl[i].type){
	case EXPFUNC: case SINFUNC: case COSFUNC: case LOGFUNC: case POWERFUNC: 
	case THRESHFUNC: case UBOUNDFUNC: case STEPFUNC: case MAXFUNC: case MINFUNC:
	case ABSFUNC: case SQRTFUNC: case SIGFUNC:
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


/// Determines if a operator is a particular type
bool Equation::optypel(const vector <EqItemList> &opl, unsigned int i, EqItemType type) const
{
	if(i == UNSET_LIST) return false;
  if(opl[i].type == type) return true;
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


/// Gives the relative priority of an operator list
int Equation::priol(const vector <EqItemList> &opl, unsigned int i) const
{
  if(i == UNSET_LIST) return -1;
  switch(opl[i].type){
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
		if(str_eq(te,i,tint)){
			i = extract_integral(te,i,op);
			
			if(warn != "") return op;
		}
		else{
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
					
					if(te.substr(i,4) == "sig(" && doneflag == false){
						item.type = SIGFUNC; op.push_back(item); 
						doneflag = true;
						i += 2;
					}
					
					if(doneflag == false){
						unsigned int raend;
						auto num = get_float(i,raend);
						if(num != UNSET){
							i = raend-1; 
							item.type = NUMERIC; 
							item.num = constant.add(num);
							op.push_back(item); 
							doneflag = true;
						}
					}

					// Looks for derived parameters
					if(doneflag == false){
						unsigned int raend;
						
						auto dref = get_derive_name(i,raend); if(warn != "") return op;
					
						if(dref.i != UNSET){
							if(type != DERIVE_EQN){
								warn = "Derived quantity '"+derive[dref.i].full_name+"' cannot be in a non-derived equation."; 
							}
							
							i = raend-1;
							
							item.type = DERIVE;
							item.num = derive_ref.size();
							op.push_back(item);
							
							derive_ref.push_back(dref);
							doneflag = true;
						}
					}
				
					// Looks for parameters
					if(doneflag == false){
						unsigned int raend;
						
						double dist = UNSET;
					
						auto pref = get_param_name(i,dist,raend); if(warn != "") return op;
					
						if(pref.th != UNSET || dist != UNSET){
							i = raend-1;
							
							if(dist != UNSET){
								item.type = NUMERIC; item.num = constant.add(dist); op.push_back(item); 
							}
							else{
								if(pref.th == TIME_VAR){
									if(ti_fix != UNSET){
										item.type = NUMERIC; 
										item.num = constant.add(timepoint[ti_fix]); 
										op.push_back(item); 
									}
									else{
										item.type = TIME; op.push_back(item); 
									}
								}
								else{
									auto &par = param[pref.th];
									
									if(par.time_dep == true){
										if(par.spline_info.on != true) emsg_input("Spline should be on");
										item.type = SPLINE;
									}							
									else item.type = PARAMETER; 
									
									auto done = false;
									
									if(done == false){ 
										//if(item.type == PARAMETER && param[pref.th].variety == CONST_PARAM){
										const auto &ele =	par.element_ref[pref.index];
										if(item.type == PARAMETER){
											if(ele.cons){
												item.type = NUMERIC; 
												item.num = ele.index; //constant.add(par.cons[ele.index]);
												op.push_back(item); 
												done = true;
											}
										}
									}
									
									if(done == false){
										if(item.type == PARAMETER && par.variety == REPARAM_PARAM && par.exist(pref.index)){
											auto num = number(par.get_value_te(pref.index));
											if(num != UNSET){
												item.type = NUMERIC; 
												item.num = constant.add(num);
												op.push_back(item); 
												done = true;
											}
										}
									}
									
									if(done == false){
										item.num = add_param_ref(pref);
										op.push_back(item); 
									
										par.add_element(pref.index,true);
									
										// Adds in factor to equation
										if(par.param_mult != UNSET){
											{
												EqItem item2; item2.type = MULTIPLY; op.push_back(item2); 
											}
												
											auto pref_fac = param_ref[item.num];
											
											if(par.time_dep){
												const auto &dep = par.dep;
												pref_fac.index /= dep[dep.size()-1].list.size();
											}
											
											const auto &par_mult = param[par.param_mult];
											const auto &dep = par_mult.dep;
											pref_fac.index *= dep[dep.size()-1].list.size();
											pref_fac.th = par.param_mult;
											
											EqItem item3;
											item3.type = SPLINE; 
											item3.num = add_param_ref(pref_fac);
											
											op.push_back(item3); 
										}
									}
									
									par.used = true;
								}
							}
							doneflag = true;
						}
					}

					// Looks for populations
					if(doneflag == false){
						unsigned int raend;
						auto numi = get_pop(i,raend); if(warn != "") return op;
						if(numi != UNSET){
							i = raend-1;
							item.type = POPNUM; item.num = numi; op.push_back(item); 
							doneflag = true;
						}
					}
				
					// Looks for individual effects
					if(doneflag == false){
						unsigned int raend;
						auto ie = get_ie(i,raend); if(warn != "") return op;
						if(ie != UNSET){
							i = raend-1;
							item.type = IE; item.num = ie; op.push_back(item); 
							doneflag = true;
						}
					}

					// Looks for fixed effects
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
						warn = "Problem with expression. The character '"+ch+"' is not expected."; 
						return op;
					}
					break;
			}
			i++;
		}
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
	
	if(op.size() > OP_MAX) warn = "The calculation is too long.";
	
	return op;
}


/// Gets an integral bound
unsigned int Equation::get_integral_bound(string st)
{
	auto t = number(st);
	if(t == UNSET){ 
		warn = "Problem with expression. Integral bound '"+st+"' is not a number.";
		return UNSET;
	}

	if(t < details.t_start || t > details.t_end){
		warn = "Problem with expression. Integral bound '"+st+"' is not within the system time period.";
		return UNSET;
	}		
	
	if(!on_timestep(t,details)){
		warn = "Problem with expression. Integral bound '"+st+"' is not at on a time-step.";
		return UNSET;
	}
	
	return get_ti_over(calc_tdiv(t,details));
}


/// Extracts information about an integral 
unsigned int Equation::extract_integral(const string &te, unsigned int i, vector <EqItem> &op)
{
	i += tint.length();
	auto ist = i;
	while(i < te.length() && te.substr(i,1) != "(") i++;
	if(i == te.length()){
		warn = "Problem with expression. There must be a right bracket '(' after the integral."; 
		return i;
	}
	
	auto cont = trim(te.substr(ist,i-ist));
	
	if(!end_str(cont,"dt")){
		if(end_str(cont,"]")) warn = "'dt' should be after the integral bounds '∫[min,max]dt'.";
		else warn = "Problem with expression. An integral must always be followed by 'dt'."; 
		return i;
	}
	
	cont = trim(cont.substr(0,cont.length()-2));
	
	auto ti_min = 0u;
	auto ti_max = details.T; 
	
	if(cont.length() != 0){
		auto fl = false;
		if(cont.substr(0,1) != "[" || cont.substr(cont.length()-1,1) != "]") fl = true;
		else{
			auto spl = split(cont.substr(1,cont.length()-2),',');
			if(spl.size() != 2) fl = true;
			else{
				auto sup = spl[1];
				auto sub = spl[0];
				ti_max = get_integral_bound(sup);
				ti_min = get_integral_bound(sub);
				
				if(ti_max == UNSET || ti_min == UNSET) return i;
				if(ti_min == ti_max){
					warn = "Problem with expression. The integral bounds cannot have the same value '"+sub+"'.";
					return i;					
				}
				else{
					if(ti_min > ti_max){
						warn = "Problem with expression. Integral bound '"+sub+"' must be less than '"+sup+"'.";
						return i;					
					}
				}
			}
		}
		
		if(fl == true){
			warn = "Problem with expression. Integral bounds should be set using the format '∫[min,max] dt (...)'."; 
			return i;
		}
	}
	
	Integral inte;
	inte.ti_min = ti_min;
	inte.ti_max = ti_max;

	EqItem item;
	item.type = TINT;
	item.num = integral.size();
	
	integral.push_back(inte);
	op.push_back(item); 

	return i;
}


/// Creates a list of operations (speeds up initialisation because doesn't need erase)
vector <EqItemList> Equation::create_opl(const vector <EqItem> &op) const
{
	vector <EqItemList> opl;
	for(auto i = 0u; i < op.size(); i++){
		const auto &ob = op[i];
		EqItemList add; 
		add.type = ob.type;
		add.num = ob.num;
		opl.push_back(add);
	}
	
	for(auto i = 0u; i < opl.size(); i++){
		if(i == 0) opl[i].prev = UNSET_LIST; else opl[i].prev = i-1;
		if(i+1 == opl.size()) opl[i].next = UNSET_LIST; else opl[i].next = i+1;
	}
	
	return opl;
}


/// Converts back to a list of operations
vector <EqItem> Equation::create_op_from_opl(const vector <EqItemList> &opl) const
{
	vector <EqItem> op;
	
	auto i = get_opl_start(opl);
	do{
		EqItem add;
		add.type = opl[i].type;
		add.num = opl[i].num;
		op.push_back(add);
		i = opl[i].next;
	}while(i != UNSET_LIST);
	
	return op;
}


/// Erases members of list
void Equation::erase_opl(unsigned int i, unsigned int n, vector <EqItemList> &opl) const
{
	auto iend = i;
	for(auto j = 0u; j < n; j++){
		if(i == UNSET_LIST) emsg_input("Problem with erase");
		iend = opl[iend].next;
	}
	
	auto iprev = opl[i].prev;
	if(iprev != UNSET_LIST) opl[iprev].next = iend;
	if(iend != UNSET_LIST) opl[iend].prev = iprev;
	
	iend = i;
	for(auto j = 0u; j < n; j++){
		auto ne = opl[iend].next;
		opl[iend].prev = CUT_LIST;
		opl[iend].next = CUT_LIST;
		iend = ne;
	}
}


/// Gets an operation from an operation list
EqItem Equation::get_op(const vector <EqItemList> &opl, unsigned int i) const
{
	EqItem op; 
	op.type = opl[i].type;
	op.num = opl[i].num;
	
	return op;
}


/// Gets the start of opl list
unsigned int Equation::get_opl_start(const vector <EqItemList> &opl) const
{
	auto i = 0u; 
	while(i < opl.size() && opl[i].prev != UNSET_LIST) i++;
	if(i == opl.size()) emsg_input("Cannot find start of equation");
	
	return i;
}		
	

/// Gets the number of elements in the opl
unsigned int Equation::get_opl_num(const vector <EqItemList> &opl) const
{
	auto i = get_opl_start(opl);
	auto num = 1u;
	while(opl[i].next != UNSET_LIST){
		i = opl[i].next;
		num++;
	}
	
	return num;
}

	
/// Works out the sequence of calculations to generate result
vector <Calculation> Equation::create_calculation(vector <EqItem> &op)
{
	vector <Calculation> calc;

	auto opl = create_opl(op);                       // Creates a list of operations
	
	auto loop = 0u;
	
	auto flag = false;
  do{                                             // Breaks down the equations into simple pieces
    flag = false;
		
		auto i = get_opl_start(opl);
		do{
			auto iprev = opl[i].prev;
			auto inext = opl[i].next;
			auto iprev2 = UNSET_LIST; if(iprev != UNSET_LIST) iprev2 = opl[iprev].prev;
			auto inext2 = UNSET_LIST; if(inext != UNSET_LIST) inext2 = opl[inext].next;
			auto inext3 = UNSET_LIST; if(inext2 != UNSET_LIST) inext3 = opl[inext2].next;
			auto inext4 = UNSET_LIST; if(inext3 != UNSET_LIST) inext4 = opl[inext3].next;
			auto inext5 = UNSET_LIST; if(inext4 != UNSET_LIST) inext5 = opl[inext4].next;
			
  		switch(opl[i].type){
			case TINT: emsg_input("Integral should not be here"); break;
			
			case LEFTBRACKET: case RIGHTBRACKET:
				break;
				
			case ADD: case TAKE: case MULTIPLY: case DIVIDE:
				if(opl[i].type == ADD && (iprev == UNSET_LIST || (opl[iprev].type == LEFTBRACKET))){ // Removes (+...)
					erase_opl(i,1,opl);
					i = inext;
					flag = true;
				}
				else{
					if(quantl(opl,iprev) && quantl(opl,inext) && 
							priol(opl,i) >= priol(opl,iprev2) && priol(opl,i) >= priol(opl,inext2)){
						auto add_flag = false;
						if(calc.size() > 0){
							auto ii = calc.size()-1;
							auto &cal = calc[ii];
							if(cal.op == opl[i].type && (cal.op == ADD || cal.op == MULTIPLY)){
								if(opl[iprev].type == REG && opl[iprev].num == ii){	
									cal.item.push_back(get_op(opl,inext));
									add_flag = true;
								} 
							}
						}

						if(add_flag == true){                 // Combines with previous operation
							erase_opl(i,2,opl);
							i = iprev;
						}
						else{                                 // Inserts a new calculation
							Calculation cal;
							cal.item.push_back(get_op(opl,iprev));
							cal.item.push_back(get_op(opl,inext));
						
							auto nreg = calc.size();
							
							cal.op = opl[i].type;
							
							calc.push_back(cal);
						
							erase_opl(i,2,opl);
							i = iprev;
							opl[i].type = REG;
							opl[i].num = nreg;
						}
						flag = true;
					}
				}
				break;
			
			// Removes brackets around a quantity
			case PARAMETER: case PARAMVEC:
			case SINGLE:
			case INTEGRAL:
			case SPLINE: case SPLINEREF: case CONSTSPLINEREF:
			case DERIVE:
			case POPNUM:
			case IE: case FE: case REG: case NUMERIC: case TIME: 
				if(optypel(opl,iprev,LEFTBRACKET) && optypel(opl,inext,RIGHTBRACKET) && !is_funcl(opl,iprev2)){
					opl[iprev].type = opl[i].type;
					opl[iprev].num = opl[i].num;
					erase_opl(i,2,opl);
					i = iprev;
					flag = true;
				}
				break;
				
			case ONE: case ZERO: break;
			
			case FUNCDIVIDE: break;
			
			case EXPFUNC: case LOGFUNC: case SINFUNC: case COSFUNC: case STEPFUNC: case ABSFUNC: case SQRTFUNC: case SIGFUNC: // Univariate
				if(optypel(opl,inext,LEFTBRACKET) && quantl(opl,inext2) == true && optypel(opl,inext3,RIGHTBRACKET)){		
					Calculation cal;
					cal.item.push_back(get_op(opl,inext2));
					
					auto nreg = calc.size();
					
					cal.op = opl[i].type;
					
					calc.push_back(cal);
					
					opl[i].type = REG;
					opl[i].num = nreg;
				
					erase_opl(inext,3,opl);
					
					flag = true;
				}
				break;
				
			case POWERFUNC: case THRESHFUNC:            // Functions with two variables
			case UBOUNDFUNC: case MAXFUNC: case MINFUNC:
				if(optypel(opl,inext,LEFTBRACKET) && quantl(opl,inext2) == true && optypel(opl,inext3,FUNCDIVIDE) && quantl(opl,inext4) == true && optypel(opl,inext5,RIGHTBRACKET)){		
					Calculation cal;
					cal.item.push_back(get_op(opl,inext2));
					cal.item.push_back(get_op(opl,inext4));
					
					auto nreg = calc.size();
					
					cal.op = opl[i].type;
					
					calc.push_back(cal);
					
					opl[i].type = REG;
					opl[i].num = nreg;
					
					erase_opl(inext,5,opl);
					
					flag = true;
				}
				break;
			
			case REG_FAC: emsg_input("Should not be reg fac"); break;
			case REG_PRECALC: emsg_input("Should not be reg precalc"); break;
			case REG_PRECALC_TIME: emsg_input("Should not be reg precalc time"); break;
			
			case NOOP: emsg_input("Hould not be no op"); break;
			}
			
			if(i == UNSET_LIST || i == CUT_LIST){ flag = true; break;}
			
			i = opl[i].next;
		}while(i != UNSET_LIST);

    if(plfl == true) print_operations(create_op_from_opl(opl));
		
		loop++;
  }while(flag == true);
	
 	string wa = "";

	if(get_opl_num(opl) != 1){
		vector <EqItemType> oper;
		auto i = 0u;
		do{
			switch(opl[i].type){
			case ADD: case TAKE: case MULTIPLY: case DIVIDE:
				oper.push_back(opl[i].type);
				break;
			default: break;
			}
			i = opl[i].next;
		}while(i != UNSET_LIST);
	
		if(oper.size() == 1){
			wa = "There is potentially a stray "+op_name(oper[0])+" in the expression. ";
		}
		else{
			wa = "The equation syntax is incorrect. ";
		}
		//emsg_input("done");
	}

	if(wa != ""){
		warn = "Problem with expression. ";
		if(wa != "error") warn += wa;	
		warn += "Calculation cannot be performed.";
		return calc;
	}

	int i = opl.size()-1;
	while(i >= 0 && opl[i].next != UNSET_LIST) i--;
	if(i < 0) emsg_input("Could not find last");
	
	auto last = get_op(opl,i);
	if(last.type != REG){
		Calculation ca;
		ca.op = ADD;
		ca.item.push_back(last);
		calc.push_back(ca);
	}
	
	return calc;
}


/// Extracts individual effects from the equation
void Equation::extract_ind_eff()
{
	auto C = calcu.size();
	
	do{
		auto ie_num = UNSET;
		auto reg_used  = UNSET;
		auto ist = UNSET, jst = UNSET;
		
		for(auto i = 0u; i < C; i++){
			const auto &ca = calcu[i];
			
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
					reg_used = i;
				}
				
				if(it.type == REG && it.num == reg_used){
					if(ca.op != MULTIPLY && !(ca.op == DIVIDE && j == 0)){
						warn = "Entire equation must by multiplied by individual effect ["+species[sp_p].ind_effect[it.num].name+"]";
						return;
					}
					reg_used = i;
				}
			}
		}
		
		if(ie_num != UNSET){
			ind_eff_mult.push_back(ie_num);
			auto &ca = calcu[ist];
			if(ca.op == MULTIPLY && ca.item.size() > 1){
				if(jst+1 < ca.item.size()){
					ca.item[jst] = ca.item[ca.item.size()-1];
				}
				ca.item.pop_back();
			}
			else calcu[ist].item[jst].type = ONE;
		}
		else break;
	}while(true);
}


/// Extracts fixed effects from the equation
void Equation::extract_fix_eff()
{
	auto C = calcu.size();
	
	do{
		auto fe_num = UNSET;
		auto reg_used  = UNSET;
		auto ist = UNSET, jst = UNSET;
	
		for(auto i = 0u; i < C; i++){
			const auto &ca = calcu[i];
			
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
					reg_used = i;
					break;
				}
				
				if(it.type == REG && it.num == reg_used){
					if(ca.op != MULTIPLY && !(ca.op == DIVIDE && j == 0)){
						warn = "Entire equation must by multiplied by fixed effect <"+species[sp_p].fix_effect[it.num].name+">";
						return;
					}
					reg_used = i;
				}
			}
		}
		
		if(fe_num != UNSET){
			fix_eff_mult.push_back(fe_num);
			auto &ca = calcu[ist];
			if(ca.op == MULTIPLY && ca.item.size() > 1){
				if(jst+1 < ca.item.size()){
					ca.item[jst] = ca.item[ca.item.size()-1];
				}
				ca.item.pop_back();
			}
			else calcu[ist].item[jst].type = ONE;
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

			auto c = claa.hash_comp.find(te);
			if(c != UNSET){
				CompRef cr; cr.p = p; cr.cl = cl; cr.c = c; cr.error = "";
				return cr;
			}
		}
	}
	
	
	CompRef cr; cr.p = UNSET; cr.cl = UNSET; cr.c = UNSET;
	if(te == "") cr.error = "No compartment specfied";
	else cr.error = "Compartment '"+te+"' not specfied";
	return cr;
}


/// Replaces a register with a constant 
void Equation::replace_reg(const vector <EqItem> &reg_replace, vector <Calculation> &calc, const vector <bool> &calc_on, bool pl)
{
	for(auto i = 0u; i < calc.size(); i++){
		if(calc_on[i]){
			auto &ca = calc[i];
	
			for(auto j = 0u; j < ca.item.size(); j++){
				auto &it = ca.item[j];
				if(it.type == REG){
					if(calc_on[it.num] == false){
						if(pl) cout << it.num << " reg replace" << endl;
						it = reg_replace[it.num];
					}
				}
			}
		}
	}	
	
	if(pl){
		cout << "after replace:" << endl;
		
		for(auto i = 0u; i < calcu.size(); i++){
			const auto &ca = calcu[i];
			print_ca(i,ca);
			cout << "    " << calc_on[i];
			cout << endl;
		}	
	}
}


/// For ADD, MULT inserts equation for one reg into a later calculation
void Equation::insert_reg(vector <Calculation> &calc, vector <bool> &calc_on)
{
	if(calc.size() == 0) return;
		
	// Works out how many times a register is used
	vector <unsigned int> reg_used(calc.size(),0);
	for(auto i = 0u; i < calc.size(); i++){
		if(calc_on[i]){
			const auto &item = calc[i].item;
			for(auto j = 0u; j < item.size(); j++){
				if(item[j].type == REG){
					auto ii = item[j].num;
					if(!calc_on[ii]) emsg_input("Calc should be on");
					reg_used[ii]++;
				}
			}
		}
	}
	
	for(auto i = 0u; i < calc.size(); i++){
		if(calc_on[i]){
			auto &ca = calc[i];
			switch(ca.op){
			case ADD: case MULTIPLY:
				{
					auto j = 0u;
					while(j < ca.item.size()){
						auto fl = false;
						
						const auto &it = ca.item[j];
						if(it.type == REG){
							auto ii = it.num;
							if(reg_used[ii] == 1){
								const auto &ca2 = calc[ii];
								if(ca2.op == ca.op){
									for(const auto &itadd : ca2.item){
										ca.item.push_back(itadd);
									}
									calc_on[ii] = false;
									fl = true;
								}
							}
						}
						
						if(fl == true){
							if(j+1 < ca.item.size()){
								ca.item[j] = ca.item[ca.item.size()-1];
							}
							ca.item.pop_back();
						}
						else{
							j++;
						}
					}
				}	
				break;
				
			default: break;
			}
		}
	}
}


/// Simplifies equations as much as possible
// (1) By removing zeros
// (2) Making used of infinity / undefined
// (3) Combining together addition and multiplication terms 
void Equation::simplify(vector <Calculation> &calc)
{
	auto pl = false;
	
	if(pl){ cout << " START SIMPLIFY" << endl; print_calculation();}
		
	// This stores if the line of the calculation is on or off
	// This is faster than using erase
	vector <bool> calc_on(calc.size(),true);
	
	vector <EqItem> reg_replace(calc.size());
	
	const auto &cval = constant.value;
	
	auto loop = 0u;
	bool flag;
	do{		
		if(pl) cout << loop << "Loop" << endl;
		
		flag = false;
	
		// Converts a negative TAKE constant to an ADD
		for(auto i = 0u; i < calc.size(); i++){
			if(calc_on[i]){
				auto &ca = calc[i];
				if(ca.op == TAKE){
					auto &it = ca.item[1];
					if(it.type == NUMERIC){
						auto val = cval[it.num];
						
						ca.op = ADD; 
						it.num = constant.add(val);
						flag = true;
					}
				}
			}
		}
		
		// Converts DIVIDE my constant to MULTIPLY
		for(auto i = 0u; i < calc.size(); i++){
			if(calc_on[i]){
				auto &ca = calc[i];
				if(ca.op == DIVIDE){
					if(ca.item[0].type != NUMERIC && ca.item[1].type == NUMERIC){
						auto val = cval[ca.item[1].num];
						if(val != 0){
							ca.op = MULTIPLY;
							ca.item[1].num = constant.add(1.0/val);
						}
					}
				}
			}
		}
		
		insert_reg(calc,calc_on);
	
		for(auto i = 0u; i < calc.size(); i++){
			if(calc_on[i]){
				auto &ca = calc[i];

				// Checks to see if there are any registers which have been turned off
				auto ill = false;
				for(auto j = 0u; j < ca.item.size(); j++){
					if(ca.item[j].type == REG && calc_on[ca.item[j].num] == false){ ill = true; break;} 
				}
				
				if(ill == false){
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
											ca.item[j].num = mult_const(ca.item[j],ca.item[jj]);
											if(jj+1 < ca.item.size()){
												ca.item[jj] = ca.item[ca.item.size()-1];
											}
											ca.item.pop_back();
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
									if(ca.item[j].type == NUMERIC && cval[ca.item[j].num] == 1){
										if(j+1 < ca.item.size()){
											ca.item[j] = ca.item[ca.item.size()-1];
										}
										ca.item.pop_back();
										flag = true;
									}
									else j++;
								}
							}
							
							// If contains zero then entire sum is zero
							auto j = 0u;
							while(j < ca.item.size() && !(ca.item[j].type == NUMERIC && cval[ca.item[j].num] == 0)) j++;
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
										ca.item[j].num = add_const(ca.item[j],ca.item[jj]);
										if(jj+1 < ca.item.size()){
											ca.item[jj] = ca.item[ca.item.size()-1];
										}
										ca.item.pop_back();
									}
									else jj++;
								}
							}
							
							// Removes anything with zero
							j = 0;
							while(j < ca.item.size()){
								if(ca.item[j].type == NUMERIC && cval[ca.item[j].num] == 0){
									if(j+1 < ca.item.size()){
										ca.item[j] = ca.item[ca.item.size()-1];
									}
									ca.item.pop_back();
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
							if(ca.item[0].type == NUMERIC) con1 = cval[ca.item[0].num];
							if(ca.item[1].type == NUMERIC) con2 = cval[ca.item[1].num];
						
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
							auto con = cval[ca.item[0].num];
							if(con == INFY) rep_con = INFY;
							else{
								if(con == UNDEF) rep_con = UNDEF;
								else rep_con = exp(con);
							}
						}					
						break;
					
					case SINFUNC:
						if(ca.item[0].type == NUMERIC){
							auto con = cval[ca.item[0].num];
							if(con == INFY) rep_con = UNDEF;
							else{
								if(con == UNDEF) rep_con = UNDEF;
								else rep_con = sin(con);
							}
						}					
						break;
						
					case COSFUNC:
						if(ca.item[0].type == NUMERIC){
							auto con = cval[ca.item[0].num];
							if(con == INFY) rep_con = UNDEF;
							else{
								if(con == UNDEF) rep_con = UNDEF;
								else rep_con = cos(con);
							}	
						}					
						break;
					
					case LOGFUNC:
						if(ca.item[0].type == NUMERIC){
							auto con = cval[ca.item[0].num];
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
							if(ca.item[0].type == NUMERIC) con1 = cval[ca.item[0].num];
							if(ca.item[1].type == NUMERIC) con2 = cval[ca.item[1].num];
						
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
							auto con1 = cval[ca.item[0].num], con2 = cval[ca.item[1].num];
							if(con1 == UNDEF || con2 == UNDEF) rep_con = UNDEF;
							else{
								if(con1 < con2) rep_con = 0;
								else rep_con = con1;
							}
						}
						break;
						
					case UBOUNDFUNC:
						if(ca.item[0].type == NUMERIC && ca.item[1].type == NUMERIC){
							auto con1 = cval[ca.item[0].num], con2 = cval[ca.item[1].num];
							if(con1 == UNDEF || con2 == UNDEF) rep_con = UNDEF;
							else{
								if(con1 > con2) rep_con = INFY;
								else rep_con = con1;
							}
						}
						break;
					
					case STEPFUNC:
						if(ca.item[0].type == NUMERIC){
							auto con = cval[ca.item[0].num];
							if(con == UNDEF) rep_con = UNDEF;
							else{
								if(con > 0) rep_con = 1; else rep_con = 0;
							}
						}
						break;
						
					case MAXFUNC:
						if(ca.item[0].type == NUMERIC && ca.item[1].type == NUMERIC){
							auto con1 = cval[ca.item[0].num], con2 = cval[ca.item[1].num];
							if(con1 == UNDEF || con2 == UNDEF) rep_con = UNDEF;
							else{
								if(con1 > con2) rep_con = con1;
								else rep_con = con2;
							}
						}
						break;
						
					case MINFUNC:
						if(ca.item[0].type == NUMERIC && ca.item[1].type == NUMERIC){
							auto con1 = cval[ca.item[0].num], con2 = cval[ca.item[1].num];
							if(con1 == UNDEF || con2 == UNDEF) rep_con = UNDEF;
							else{
								if(con1 < con2) rep_con = con1;
								else rep_con = con2;
							}
						}
						break;
						
					case ABSFUNC:
						if(ca.item[0].type == NUMERIC){
							auto con = cval[ca.item[0].num];
							if(con == UNDEF) rep_con = UNDEF;
							else{
								if(con < 0) rep_con = -con;
								else rep_con = con;
							}
						}
						break;
					
					case SQRTFUNC:
						if(ca.item[0].type == NUMERIC){
							auto con = cval[ca.item[0].num];
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
						
					case SIGFUNC:
						if(ca.item[0].type == NUMERIC){
							auto con = cval[ca.item[0].num];
							
							if(con == UNDEF) rep_con = UNDEF;
							else{
								if(con == INFY) rep_con = 1;
								else rep_con = 1/(1+exp(-con));
							}
						}
						break;

					default: emsg_input("Should be simplification"); break;
					}
					
					if(rep.type != NOOP || rep_con != UNSET){
						if(i+1 != calc.size()){					
							flag = true;
							if(rep.type != NOOP && rep_con != UNSET) emsg_input("Cannot be both here");
							
							if(rep_con != UNSET){
								rep.type = NUMERIC; rep.num = constant.add(rep_con);
							}
							
							reg_replace[i] = rep;
							calc_on[i] = false;
						}
					}
				}
			}
		}
		
		replace_reg(reg_replace,calc,calc_on,pl);
		
		// If last entry just refers to previous register then remove
		if(calcu.size() > 0){
			const auto &la = calcu[calcu.size()-1];
			if((la.op == ADD || la.op == MULTIPLY) && la.item.size() == 1){
				if(la.item[0].type == REG && la.item[0].num+2 == calcu.size()){
					calcu.pop_back();
					flag = true;
				}
			}
		}
		
		loop++; if(loop > 100) emsg_input("Simplify equation problem");
	}while(flag == true);

	remove_unused(calc,calc_on);
	
	// Checks that INFY and UNDEF do to exist in the final equation
	for(auto i = 0u; i < calc.size(); i++){
		const auto &ca = calc[i];
		for(auto j = 0u; j < ca.item.size(); j++){
			const auto &it = ca.item[j];
			if(it.type == NUMERIC){
				if(cval[it.num] == UNDEF) warn = "The equation contains an undefined quantity";
				if(cval[it.num] == INFY) warn = "The equation contains an infinite quantity";
			}
		}
	}
	
	if(false) print_calculation();
}


/// Removes any unused registers
void Equation::remove_unused(vector <Calculation> &calc, vector <bool> &calc_on)
{
	if(false) print_calculation();

	auto C = calc.size();
	if(C == 0) emsg_input("calc no size");
	
	vector <bool> used(C,false);
	
	int j = C-1;
	if(calc_on[j] != true) emsg_input("Last should not change");
	
	used[j] = true;
	while(j >= 0){
		if(used[j] == true){
			for(const auto &it : calc[j].item){
				if(it.type == REG){
					auto jj = it.num;
					if(calc_on[jj] != true) emsg_input("calc should be on");
					used[jj] = true;
				}
			}
		}
		j--;
	}
	
	// Turns off any unused
	for(auto j = 0u; j < calc.size(); j++){
		if(calc_on[j] && !used[j]) calc_on[j] = false;
	}
	
	auto calc_old = calc;
	calc.clear();
	
	vector <unsigned int> map;
	for(auto i = 0u; i < calc_old.size(); i++){
		if(calc_on[i] == false){
			map.push_back(UNSET_LIST);
		}
		else{
			map.push_back(calc.size());
			calc.push_back(calc_old[i]);
		}
	}
	
	for(auto i = 0u; i < calc.size(); i++){
		auto &ca = calc[i];
		for(auto j = 0u; j < ca.item.size(); j++){
			auto &it = ca.item[j];
			if(it.type == REG){
				auto num = map[it.num];
				if(num == UNSET_LIST) emsg_input("calc map problem");
				it.num = num;
			}
		}
	}
}

 
/// Adds constant values
unsigned int Equation::add_const(EqItem item1, EqItem item2)
{
	auto val1 = constant.value[item1.num], val2 = constant.value[item2.num];
	
	double num;
	if(val1 == UNDEF || val2 == UNDEF) num = UNDEF; 
	else{
		if(val1 == INFY || val2 == INFY) num = INFY;
		else num = val1+val2;
	}
	
	return constant.add(num);
}


/// Multiplies constant values
unsigned int Equation::mult_const(EqItem item1, EqItem item2)
{
	auto val1 = constant.value[item1.num], val2 = constant.value[item2.num];

	double num;
	if(val1 == UNDEF || val2 == UNDEF) num = UNDEF;
	else{
		if(val1 == INFY){
			if(val2 == 0) num = UNDEF;
			else num = INFY;
		}
		else{
			if(val2 == INFY){
				if(val1 == 0) num = UNDEF;
				else num = INFY;
			}
			else{
				num = val1*val2;
			}
		}
	}
	
	return constant.add(num);
}


/// Calculates a list of all population used 
void Equation::calculate_pop_ref()
{				
	pop_ref.clear();
	HashSimp hash_pr;
	for(auto i = 0u; i < calcu.size(); i++){
		const auto &ca = calcu[i];
		for(const auto &it : ca.item){
			if(it.type == POPNUM) add_to_vec(pop_ref,it.num,hash_pr);
		}
	}
	
	for(const auto &inte : integral){
		const auto &calc = inte.calc;
		for(auto i = 0u; i < calc.size(); i++){
			const auto &ca = calc[i];
			for(const auto &it : ca.item){
				if(it.type == POPNUM) add_to_vec(pop_ref,it.num,hash_pr);
			}
		}
	}
}		


/// Set if equation is a simple number
double Equation::is_num() const
{
	if(calcu.size() == 1){
		const auto &ca = calcu[0];
		if(ca.op == ADD || ca.op == MULTIPLY){ 	
			if(ca.item.size() == 1){
				const auto &it = ca.item[0];
				if(it.type == NUMERIC) return constant.value[it.num];
			}
		}
	}
	
	return UNSET;
}


/// Determines if the equation is zero
bool Equation::is_zero() const
{
	if(calcu.size() == 1){
		const auto &ca = calcu[0];
		if(ca.op == ADD || ca.op == MULTIPLY){ 
			if(ca.op == ADD && ca.item.size() == 0) return true;
			
			if(ca.item.size() == 1){
				const auto &it = ca.item[0];
				if(it.type == NUMERIC && constant.value[it.num] == 0) return true;
			}
		}
	}
	
	return false;
}


/// Determines if the equation is one
bool Equation::is_one() const
{
	if(calcu.size() == 1){
		const auto &ca = calcu[0];
		if(ca.op == ADD || ca.op == MULTIPLY){ 	
			if(ca.op == MULTIPLY && ca.item.size() == 0) return true;
			
			if(ca.item.size() == 1){
				const auto &it = ca.item[0];
				if(it.type == NUMERIC && constant.value[it.num] == 1) return true;
			}
		}
	}
	
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
	
  for(const auto &ca : calcu){
		for(const auto &it : ca.item){
			if(it.type == TIME) time_vari = true;
		}
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


/// Determines if equation is linear in a given parameter
bool Equation::param_linear(unsigned int th) const
{
	if(false) cout << th << endl;
	// This is only used for joint parameter/event proposals and these are turned off
	/*
	if(calcu.size() == 1){
		const auto &ca = calcu[0];
		if(ca.op == ADD || ca.op == MULTIPLY){ 	
			if(ca.item.size() == 1){
				const auto &it = ca.item[0];
				
	switch(ans.type){
	case PARAMETER: case PARAMVEC:
		if(ans.num == th) return true;
		break;
		
	case REG:
		{
			const auto &last = calc[calc.size()-1];
			if(last.reg_store != ans.num) emsg("Reg does not agree");
			
			switch(last.op){
			case MULTIPLY:
				for(const auto &it : last.item){
					if(it.type == PARAMVEC && it.num == th) return true;
				}
				break;
				
			case ADD:
				if(last.item.size() == 1){
					const auto &it = last.item[0];
					if(it.type == PARAMVEC && it.num == th) return true;
				}
				break;
			
			case DIVIDE:
				{
					const auto &it = last.item[0];
					if(it.type == PARAMVEC && it.num == th) return true;
				}
				break;
			
			default: break;
			}
		}
		break;
		
	default: break;
	}
	*/
	return false; 
}
	

/// Removes any unused parameter references (based on simplifications)
void Equation::remove_unused_param_ref()
{
	auto N = param_ref.size();
	vector <bool> on(N,false);
	
	for(auto i = 0u; i < calcu.size(); i++){
		const auto &ca = calcu[i];
		for(const auto &it : ca.item){
			switch(it.type){
			case PARAMETER: case PARAMVEC: case SPLINE: case SPLINEREF: on[it.num] = true; break;
			default: break;
			}
		}
	}
	
	for(const auto &inte : integral){
		const auto &calc = inte.calc;
		for(auto i = 0u; i < calc.size(); i++){
			const auto &ca = calc[i];
			for(const auto &it : ca.item){
				switch(it.type){
				case PARAMETER: case PARAMVEC: case SPLINE: case SPLINEREF: 
					on[it.num] = true; 
					break;
				default: break;
				}
			}
		}
	}
	
	auto j = 0u; while(j < N && on[j]) j++;
	if(j == N) return;
	
	// Removes parameters which are not used
	vector <unsigned int> map;
	auto param_ref_old = param_ref;
	param_ref.clear();
	for(auto i = 0u; i < N; i++){
		if(on[i]){
			map.push_back(param_ref.size());
			param_ref.push_back(param_ref_old[i]);
		}
		else{
			map.push_back(UNSET);
		}
	}
	
	for(auto i = 0u; i < calcu.size(); i++){
		auto &ca = calcu[i];
		for(auto &it : ca.item){
			switch(it.type){
			case PARAMETER: case PARAMVEC: case SPLINE: case SPLINEREF: 
				it.num = map[it.num];
				break;
			default: break;
			}
		}
	}
	
	for(auto &inte : integral){
		auto &calc = inte.calc;
		for(auto i = 0u; i < calc.size(); i++){
			auto &ca = calc[i];
			for(auto &it : ca.item){
				switch(it.type){
				case PARAMETER: case PARAMVEC: case SPLINE: case SPLINEREF: 
					it.num = map[it.num];
					break;
				default: break;
				}
			}
		}
	}
}


/// Sets up references (pop_ref, param_ref) in equations
void Equation::setup_references()
{
	remove_unused_param_ref();
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

	emsg_input("Cannot find distance");
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
			
			auto c = claa.hash_comp.find(dep1);
			if(c != UNSET){
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


/// Gets the distance between two compartments
double Equation::get_identity(const ParamProp &pp)
{
	if(pp.dep.size() != 2){	
		warn = "The identity matrix '"+pp.name+"' must have two indices."; 
		return UNSET;
	}
	
	auto dep1 = pp.dep[0];
	auto dep2 = pp.dep[1];

	for(auto p = 0u; p < species.size(); p++){
		const auto &sp = species[p];
		for(auto cl = 0u; cl < sp.cla.size(); cl++){
			const auto &claa = sp.cla[cl];
			
			auto c = claa.hash_comp.find(dep1);
			if(c != UNSET){
				auto cc = claa.hash_comp.find(dep2);
				if(cc == UNSET){
					warn = "Compartment '"+dep2+"' is not recognised in distance matrix.";
					return UNSET;
				}
				else{
					if(c == cc) return 1;
					else return 0;
				}
			}
		}
	}
	
	warn = "Compartment '"+dep1+"' is not recognised in identity matrix.";
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
	// Cannot be empty
	if(op.size() == 0){
		warn = "Expression is empty.";
		return;
	}
	
	// Cannot have just an operator
	if(op.size() == 1){
		switch(op[0].type){
		case ADD: case TAKE: case MULTIPLY: case DIVIDE:
			warn = "Expression cannot be simply "+op_name(op[0].type)+".";
			return;
		default: break;
		}
	}
	
	// Cannot have "(*..." or "(/..." or "*..." or "/..."
	// Cannot have "...+-*/)" or "...+-*/" 
	for(auto i = 0u; i < op.size(); i++){
		auto ty = op[i].type;
		switch(ty){
		case ADD: case TAKE: case MULTIPLY: case DIVIDE:
			if(i+1 == op.size()){
				warn = "Problem with the expression. Cannot have operator "+op_name(ty)+" at the end";
				return;
			}
			else{
				if(op[i+1].type == RIGHTBRACKET){
					warn = "Problem with the expression. Cannot have operator "+op_name(ty)+" next to right bracket ')'";
					break;
				}
			}
			
			if(i == 0){
				if(ty == MULTIPLY || ty == DIVIDE){
					warn = "Problem with the expression. Cannot have operator "+op_name(ty)+" at the start";
					return;
				}
			}
			else{
				if(op[i-1].type == LEFTBRACKET && (ty == MULTIPLY || ty == DIVIDE)){
					warn = "Problem with the expression. Cannot have operator "+op_name(ty)+" after left bracket '('";
					return;
				}
			}
			break;
		
		default: break;
		}
	}
	
	if(op.size() > 0){
		// Cannot have operators next to  each other
		for(auto i = 0u; i < op.size()-1; i++){
			switch(op[i].type){
			case ADD: case TAKE: case MULTIPLY: case DIVIDE:
				switch(op[i+1].type){
				case ADD: case TAKE: case MULTIPLY: case DIVIDE:
					warn = "Problem with the expression. Cannot have operators "+op_name(op[i].type)+" and "+op_name(op[i+1].type)+" next to each other";
					return;
				default: break;
				}
			default: break;
			}
		}
	
		// Cannot have empty brackets
		for(auto i = 0u; i < op.size()-1; i++){
			if(op[i].type == LEFTBRACKET && op[i+1].type == RIGHTBRACKET){
				warn = "Problem with the expression. Cannot have brackets '()' with no content.";
				return;
			}
		}			
	}
}


/// Replaces time integral with separate calculation
void Equation::time_integral(vector <EqItem> &op)
{
	auto i = 0u; 
	while(i < op.size()){
		while(i < op.size() && op[i].type != TINT) i++;
		
		if(i < op.size()){
			auto e = op[i].num;
			auto &inte = integral[e];
			
			// Works out content of integral
			auto ist = i+1;
			
			auto j = ist;
			if(!(j < op.size() && op[j].type == LEFTBRACKET)){
				warn = "Problem with the expression. A time integral must be followed by a left bracket '('.";
				return;
			}
			
			auto num = 0u;
			while(j < op.size()){
				auto ty = op[j].type;
				if(ty == LEFTBRACKET) num++;
				if(ty == RIGHTBRACKET){
					num--;
					if(num == 0) break;
				}
				if(ty == TINT){
					warn = "Problem with the expression. Cannot have nested integrals.";
					return;
				}
				j++;
			}
			
			if(j == op.size()){
				warn = "Problem with the expression. Content in time integral must be enclosed by round brackets '(...)'.";
				return;
			}
			
			if(ist+1 == j){
				warn = "Problem with the expression. The time integral does not contain any content.";
				return;
			}
			
			// Creates mini calculation within integral
			vector <EqItem> op_mini;
			for(auto k = ist+1; k < j; k++){
				op_mini.push_back(op[k]);
			}
			
			auto calc = create_calculation(op_mini);
			if(simplify_eqn == true) simplify(calc);
			
			inte.calc = calc;
			
			// Removes calculation
			op.erase(op.begin()+i,op.begin()+j+1);
			
			// Adds intergal reference term
			EqItem add;
			add.type = INTEGRAL;
			add.num = e;
			op.insert(op.begin()+i,add);
		
			if(type != DERIVE_EQN){
				warn = "Time integrals can only appear in derived equations.";
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
			
			EqItem item2; item2.type = NUMERIC; item2.num = constant.add(-1);
			op.insert(op.begin()+i+1,item2);
	
			EqItem item3; item3.type = MULTIPLY;
			op.insert(op.begin()+i+2,item3);
		}
		else i++;
	}
}


/// Adds a parameter reference
unsigned int Equation::add_param_ref(const ParamRef &pref)
{
	auto num = param_ref.size();
	param_ref.push_back(pref);
	return num;
}


/// Gets a number to represent an equation
double Equation::get_calc_hash_num(const vector <Calculation> &calc) const
{
	const auto &cval = constant.value;
	
	auto val = 0.0;
	for(auto i = 0u; i < calc.size(); i++){
		const auto &ca = calc[i];
		for(const auto &it : ca.item){
			val += MM_PI*it.type;
			if(it.type == NUMERIC) val += cval[it.num];
			else val += 0.3*it.num;
		}
		val += 0.1*ca.op + i;
		val *= MM_PI/3;
	}
	
	return val;
}

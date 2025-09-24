// This file deals with linearisation of equations (in terms of populations and splines).
// This was introduces as a way to speed up the farm-based models 

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
 
#include "math.h"

using namespace std;

#include "equation.hh"
#include "utils.hh"

/// Works out if possible to linearise the equation in terms of populations/spline/terms involving time 
/// This is to help speed up likelihood calculation.
void Equation::calculate_linearise()
{
	if(linearise.on == true) return;

	auto C = calcu.size();

	auto pl = false;

	if(pl){
		cout << endl << endl << "START" << endl;
		print_calculation();
	}
	
	vector <LinearCalculation> lin_calc;
	for(auto i = 0u; i < C; i++){
		const auto &ca = calcu[i];
		
		if(i != 0 && ca.item.size() == 0) emsg_input("zero item problem");
		
		LinearCalculation lin;
		switch(ca.op){
		case ADD: 
			for(auto j = 0u; j < ca.item.size(); j++){
				auto lin_it = convert_to_linear_calculation(ca.item[j],ADD,lin_calc);
				if(j == 0) lin = lin_it;
				else{
					calc_add(lin.no_pop_calc.calc,lin_it.no_pop_calc.calc);
					if(pl){
						print_calc("calcfin",lin.no_pop_calc.calc);
					}
					for(auto k = 0u; k < lin_it.pop_calc.size(); k++){
						const auto &pc = lin_it.pop_calc[k];
						auto po = pc.po;
						auto m = 0u; while(m < lin.pop_calc.size() && lin.pop_calc[m].po != po) m++;
						
						if(m == lin.pop_calc.size()){
							lin.pop_calc.push_back(pc);
						}
						else{
							calc_add(lin.pop_calc[m].calc,lin_it.pop_calc[k].calc);
						}
					}
				}
			}
			break;
		
		case MULTIPLY: 
			for(auto j = 0u; j < ca.item.size(); j++){
				auto lin_it = convert_to_linear_calculation(ca.item[j],MULTIPLY,lin_calc);
				if(j == 0) lin = lin_it;
				else{		
					if(lin.pop_calc.size() > 0 && lin_it.pop_calc.size() > 0){ 
						return; // Stops quadratic terms in population
					}
					
					if(false){
						print_linear_calc("Lin1 ",lin);
						print_linear_calc("Lin2 ",lin_it);
					}
					
					// Multiplies population terms in lin with non_pop in lin_it
					for(auto k = 0u; k < lin.pop_calc.size(); k++){
						auto &pc = lin.pop_calc[k];
						calc_mult(pc.calc,lin_it.no_pop_calc.calc);
					}
					
					// Multiplies population terms in lin_it with non_pop in lin
					for(auto k = 0u; k < lin_it.pop_calc.size(); k++){
						auto &pc = lin_it.pop_calc[k];
						calc_mult(pc.calc,lin.no_pop_calc.calc);
						lin.pop_calc.push_back(pc);
					}
					
					// Multiplies non_pop in both
					calc_mult(lin.no_pop_calc.calc,lin_it.no_pop_calc.calc);
					
					if(false){
						print_linear_calc("Lin after ",lin);
						emsg_input("fin");
					}
				}
			}
			break;
			
		case DIVIDE:
			{
				lin = convert_to_linear_calculation(ca.item[0],MULTIPLY,lin_calc);
				auto lin_div = convert_to_linear_calculation(ca.item[1],MULTIPLY,lin_calc);
				if(lin_div.pop_calc.size() > 0) return;
				if(lin_div.no_pop_calc.calc.size() == 0) emsg_input("Problem with div");
				
				calc_div(lin.no_pop_calc.calc,lin_div.no_pop_calc.calc);
				for(auto j = 0u; j < lin.pop_calc.size(); j++){
					calc_div(lin.pop_calc[j].calc,lin_div.no_pop_calc.calc);
				}
				
				if(false){
					print_linear_calc("Lin after ",lin);
					emsg_input("fin");
				}
			}
			break;
		
		// Single parameter functions
		case EXPFUNC: case SINFUNC: case COSFUNC: case LOGFUNC:  
		case STEPFUNC: case ABSFUNC: case SQRTFUNC: case SIGFUNC:
			if(single_param_func(ca,lin,lin_calc) == false) return;
			break;
			
		// Two parameter functions
		case POWERFUNC: case THRESHFUNC: case UBOUNDFUNC: case MAXFUNC: case MINFUNC:
			if(two_param_func(ca,lin,lin_calc) == false) return;
			break;
	
		case TAKE: emsg_input("Should not have take"); break;
		
		default: emsg_input("Eq problem10"); break;
		}
	
		if(pl && false){
			cout << endl;
			print_linear_calc("LIN FINAL "+to_string(i),lin);
			cout << endl;
		}
		
		lin_calc.push_back(lin);
	}

	if(pl){
		for(const auto &lin : lin_calc){
			print_linear_calc("LIN FINAL ",lin);
		}
	}

	//auto lc_final = convert_to_linear_calculation(ans,ADD,lin_calc);
	auto lc_final = lin_calc[lin_calc.size()-1];
	
	{
		//linearise.no_pop_calc.type = ZERO; 
		
		auto &calc = lc_final.no_pop_calc.calc;
		
		if(calc.size() == 0){ // Sets up a zero sum
			Calculation ca; ca.op = ADD;
			EqItem it; it. type = ZERO; it.num = UNSET;
			ca.item.push_back(it);
			calc.push_back(ca);
		}
		else simplify(calc); 
		
		linearise.no_pop_calc_store = calc;
	}
	
	if(false) print_linear_final();
	
	auto Npop = pop_ref.size();
	if(Npop != lc_final.pop_calc.size()){
		emsg_input("Pop_calc problem");
	}
	
	// Orders gradients in the same way as pop_ref
	for(auto i = 0u; i < Npop; i++){
		auto po = pop_ref[i];
		
		auto k = 0u; while(k < Npop && lc_final.pop_calc[k].po != po) k++;
		if(k == Npop) emsg_input("Could not find population");
		
		auto &calc = lc_final.pop_calc[k].calc; 
		linearise.pop_grad_calc_store.push_back(calc);
		
		/*
		if(calc.size() != 1) return;
		if(calc[0].op != MULTIPLY) return;
		linearise.pop_grad_calc.push_back(calc[0].item);
		*/
	}

	linearise.multi_source = false;
	
	auto num = linearise.pop_grad_calc_store.size();
	if(linearise.no_pop_calc_store.size() > 0) num++;
	if(num > 1) linearise.multi_source = true; 
	
	linearise.on = true;
	
	linearise.init_pop_ref_from_po(pop_ref);
		
	if(linearise.pop_grad_calc_store.size() != pop_ref.size()){
		emsg_input("Population number does not agree");
	}
	
	simplify(linearise.no_pop_calc_store);
	for(auto &ca : linearise.pop_grad_calc_store) simplify(ca);
		
	if(false){
		print_calculation();
		cout << endl << "AFTER LINEARISE" << endl;
		cout << endl << endl << te << ":" << endl;
		print_calc("No pop dep: ",linearise.no_pop_calc_store);

		for(auto k = 0u; k < linearise.pop_grad_calc_store.size(); k++){
			print_calc(pop[pop_ref[k]].name,linearise.pop_grad_calc_store[k]);	
		}
		emsg_input("Linear done");
	}
	
	get_pop_grad_calc_factorise();
	
	// Adds a one term if there is no term
	for(auto &calc : linearise.pop_grad_calc_store){
		if(calc.size() == 1){
			auto &ca = calc[0];
			if(ca.op == MULTIPLY && ca.item.size() == 0){
				EqItem it; it.type = ONE; it.num = UNSET;
				ca.item.push_back(it);
			}
		}
	}
}

	
/// Prints the calculation and the linear final version (for diagnostics)
void Equation::print_linear_final() const
{
	cout << te << " te" << endl;
	print_calculation();
	cout << endl << endl;

	cout << endl << endl << te << ":" << endl;
	print_calc("No pop dep: ",linearise.no_pop_calc_store);

	for(auto k = 0u; k < linearise.pop_grad_calc_store.size(); k++){
		print_calc(pop[pop_ref[k]].name,linearise.pop_grad_calc_store[k]);	
	}
}
	

/// Converts from a item to a linear calculation
LinearCalculation Equation::convert_to_linear_calculation(const EqItem &it, EqItemType op, const vector <LinearCalculation> &lin_calc)
{
	LinearCalculation lin;
	if(warn != "") return lin;
	
	lin.no_pop_calc.po = UNSET;
	switch(it.type)
	{
		case REG: 
			lin = lin_calc[it.num];
			break;
			
		case POPNUM: 
			{
				Calculation ca; ca.op = op;
				EqItem it2; it2.type = NUMERIC; it2.num = constant.add(1);
				ca.item.push_back(it2);
				
				PopCalculation pop_c;
				pop_c.po = it.num;
				pop_c.calc.push_back(ca);
				lin.pop_calc.push_back(pop_c);
			}
			break;
		
		case PARAMVEC: case SPLINEREF: case CONSTSPLINEREF: case NUMERIC: case TIME:
		case REG_PRECALC: case REG_PRECALC_TIME:
			{
				Calculation ca; ca.op = op; ca.item.push_back(it);
				lin.no_pop_calc.calc.push_back(ca);
			}
			break;
		
		case IE: case FE: case ONE: case ZERO: emsg_input("Should not be here7"); break;			
		default: emsg_input("Equation error11"); break;
	}
	
	return lin;
}


/// Multiplies two calculation
void Equation::calc_mult(vector <Calculation> &calc, const vector <Calculation> &calc2) 
{
	if(false){
		cout << "CALC MULT" << endl;
		print_calc("calc",calc);
		print_calc("calc2",calc2);
	}
	
	if(calc.size() == 0) return;
	if(calc2.size() == 0){ calc.resize(0); return;}
	
	const auto &cval = constant.value;
	
	if(true){ // Extra code to make more efficient
		if(calc.size() == 1){
			const auto &ca = calc[0];
			if(ca_is_one(ca)){ calc = calc2; return;}	
			else{			
				if(ca.op == MULTIPLY){
					if(calc2[calc2.size()-1].op == MULTIPLY){
						auto it_store = ca.item;
						calc = calc2;
						auto &ca_last = calc[calc.size()-1];
						for(auto j = 0u; j < it_store.size(); j++){
							const auto &it = it_store[j];
							if(it.type == NUMERIC){
								auto k = 0u; while(k < ca_last.item.size() && ca_last.item[k].type != NUMERIC) k++;
								if(k == ca_last.item.size()) ca_last.item.push_back(it);
								else ca_last.item[k].num = mult_const(ca_last.item[k],it);
							}
							else ca_last.item.push_back(it);
						}
						return;
					}
					else{
						if(ca.item.size() == 1 && ca.item[0].type == NUMERIC){
							auto con = cval[ca.item[0].num];
							calc = calc2;
							Calculation ca; ca.op = MULTIPLY; 
							EqItem it1; it1.type = REG; it1.num = calc.size()-1;
							EqItem it2; it2.type = NUMERIC; it2.num = constant.add(con);
							ca.item.push_back(it1);
							ca.item.push_back(it2);
							calc.push_back(ca);
							return;
						}
					}
				}
			}
		}
		
		if(calc2.size() == 1){
			const auto &ca2 = calc2[0];
		
			if(ca2.op == MULTIPLY){
				if(ca2.item.size() == 1 && ca2.item[0].type == NUMERIC && cval[ca2.item[0].num] == 1){
					return;
				}					
				else{
					auto &ca_last = calc[calc.size()-1];
					if(ca_last.op == MULTIPLY){
					for(auto j = 0u; j < ca2.item.size(); j++){
							const auto &it2 = ca2.item[j];
							if(it2.type == NUMERIC){
								auto k = 0u; while(k < ca_last.item.size() && ca_last.item[k].type != NUMERIC) k++;
								if(k == ca_last.item.size()) ca_last.item.push_back(it2);
								else ca_last.item[k].num = mult_const(ca_last.item[k],it2);
							}
							else ca_last.item.push_back(it2);
						}
						return;
					}
					else{
						if(ca2.item.size() == 1 && ca2.item[0].type == NUMERIC){
							auto con = cval[ca2.item[0].num];
								
							Calculation ca; ca.op = MULTIPLY;
							EqItem it1; it1.type = REG; it1.num = calc.size()-1;
							EqItem it2; it2.type = NUMERIC; it2.num = constant.add(con);
							ca.item.push_back(it1);
							ca.item.push_back(it2);
							calc.push_back(ca);
							return;
						}
					}
				}
			}
		}
	}
	
	auto num = calc2.size();

	// Shifts forward registers by num
	for(auto i = 0u; i < calc.size(); i++){
		for(auto &it : calc[i].item){
			if(it.type == REG) it.num += num;
		}
		//calc[i].reg_store += num;
	}

	for(auto i = 0u; i < num; i++){
		calc.insert(calc.begin()+i,calc2[i]);
	}
	
	auto &ca_last = calc[calc.size()-1];
	if(ca_last.op == MULTIPLY){
		EqItem it; it.type = REG; it.num = num-1;
		ca_last.item.push_back(it);
	}
	else{
		Calculation ca; ca.op = MULTIPLY; 
		EqItem it1; it1.type = REG; it1.num = num-1;
		EqItem it2; it2.type = REG; it2.num = calc.size()-1;
		ca.item.push_back(it1);
		ca.item.push_back(it2);
		calc.push_back(ca);
	}
	
	if(false) print_calc("combine",calc);
}


/// Multiplies two calculation
void Equation::calc_add(vector <Calculation> &calc, const vector <Calculation> &calc2) const 
{
	if(false){
		print_calc("calc",calc);
		print_calc("calc2",calc2);
		cout << "ADD" << endl;
		emsg_input("here add");
	}
	
	if(calc2.size() == 0) return;
	
	if(calc.size() == 0){
		calc = calc2; return;
	}
	
	if(true){ // Extra code to make more efficient
		if(calc2.size() == 1 && calc.size() >= 1){
			const auto &ca2 = calc2[0];
			auto &ca_last = calc[calc.size()-1];
			if(ca2.op == ADD && ca_last.op == ADD){
				for(const auto &it : ca2.item) ca_last.item.push_back(it);
				return;
			}
		}
		
		if(calc.size() == 1 && calc2.size() >= 1){
			if(calc[0].op == ADD && calc2[calc2.size()-1].op == ADD){
				auto it_store = calc[0].item;
				calc = calc2;
				auto &ca_last = calc[calc.size()-1];
				for(const auto &it : it_store) ca_last.item.push_back(it);
				return;
			}
		}
	}
	
	auto num = calc2.size();

	// Shifts forward registers by num
	for(auto i = 0u; i < calc.size(); i++){
		for(auto &it : calc[i].item){
			if(it.type == REG) it.num += num;
		}
		//calc[i].reg_store += num;
	}

	for(auto i = 0u; i < num; i++){
		calc.insert(calc.begin()+i,calc2[i]);
	}
	
	auto &ca_last = calc[calc.size()-1];
	if(ca_last.op == ADD){
		EqItem it; it.type = REG; it.num = num-1;
		ca_last.item.push_back(it);
	}
	else{
		Calculation ca; ca.op = ADD;
		EqItem it1; it1.type = REG; it1.num = num-1;
		EqItem it2; it2.type = REG; it2.num = calc.size()-1;
		ca.item.push_back(it1);
		ca.item.push_back(it2);
		calc.push_back(ca);
	}
	
	if(false) print_calc("combine add",calc);
}


/// Detemines if a calculation is one
bool Equation::ca_is_one(const Calculation &ca) const
{
	if(ca.item.size() == 1 && ca.item[0].type == NUMERIC && constant.value[ca.item[0].num] == 1){
		if(ca.op == MULTIPLY || ca.op == ADD) return true;
	}
	return false;
}
	

/// Divides two calculations
void Equation::calc_div(vector <Calculation> &calc, const vector <Calculation> &calc2) const 
{
	if(false){
		print_calc("calc",calc);
		print_calc("calc2",calc2);
	
		cout << "DIV" << endl;
		emsg_input("here div");
	}
	
	if(calc2.size() == 0) emsg_input("Problem");
	if(calc.size() == 0) return;
		
	auto num = calc.size();
	
	for(auto ca : calc2){
		for(auto &it : ca.item){
			if(it.type == REG) it.num += num;
		}
		//ca.reg_store += num;
		
		calc.push_back(ca);
	}
	
	Calculation ca; ca.op = DIVIDE;
	EqItem it1; it1.type = REG; it1.num = num-1;
	EqItem it2; it2.type = REG; it2.num = calc.size()-1;
	ca.item.push_back(it1);
	ca.item.push_back(it2);
	calc.push_back(ca);
	
	if(false) print_calc("combine div",calc);
}


/// Prints a linear calculation 	
void Equation::print_linear_calc(string te, const LinearCalculation &lin) const 
{
	cout << te << ":" << endl;
	print_calc("No pop dep",lin.no_pop_calc.calc);
	
	for(auto i = 0u; i < lin.pop_calc.size(); i++){
		const auto &pc = lin.pop_calc[i];
		print_calc("Pop "+pop[pc.po].name,pc.calc);
	}
}


/// Prints a calculation
void Equation::print_calc(string st, const vector <Calculation> &calc) const
{
	cout << st << ": "; 
	for(auto i = 0u; i < calc.size(); i++){
		const auto &ca = calc[i];
		print_ca(i,ca); 
		cout << ",  ";
	}
	cout << endl;
}


/// Finds the linear calculation for a single parameter function
bool Equation::single_param_func(Calculation ca, LinearCalculation &lin, const vector <LinearCalculation> &lin_calc) const
{
	const auto &it = ca.item[0];

	switch(it.type){
	case REG:
		{
			lin = lin_calc[it.num];
			if(lin.pop_calc.size() > 0) return false;
			auto &cal = lin.no_pop_calc.calc;
			ca.item[0].num = cal.size()-1;
			cal.push_back(ca);
		
			if(false) print_linear_calc("single func", lin);
		}
		break;
		
	case PARAMVEC: case SPLINEREF: case CONSTSPLINEREF: case NUMERIC: case TIME: 
		lin.no_pop_calc.calc.push_back(ca);
		break;
		
	case POPNUM: return false;
	case IE: case ONE: case ZERO: case FE: emsg_input("Eq Lin should not be"); break;
	default: emsg_input("Eq problem2"); break;
	}
	
	return true;
}


/// Finds the linear calculation for a two parameter function
bool Equation::two_param_func(Calculation ca, LinearCalculation &lin, const vector <LinearCalculation> &lin_calc) const
{
	auto &it1 = ca.item[0];
	auto &it2 = ca.item[1];
	
	switch(it1.type){
	case REG:
		switch(it2.type){
		case REG:
			{		
				lin = lin_calc[it1.num];
				if(lin.pop_calc.size() > 0) return false;
				
				const auto &lin2 = lin_calc[it2.num];
				if(lin2.pop_calc.size() > 0) return false;
				
				auto &calc = lin.no_pop_calc.calc;
				
				auto num = calc.size();
				
				const auto &calc2 = lin2.no_pop_calc.calc;
				for(auto ca_add : calc2){
					for(auto &it : ca_add.item){
						if(it.type == REG) it.num += num;
					}
					//ca_add.reg_store += num;
					calc.push_back(ca_add);
				}
				
				it1.num = num-1; 
				it2.num = calc.size()-1;
				calc.push_back(ca);
				if(false){
					print_linear_calc("calc",lin); 
					emsg_input("K");
				}
			}
			break;
			
		case PARAMVEC: case SPLINEREF: case CONSTSPLINEREF: case NUMERIC: case TIME:
			{		
				lin = lin_calc[it1.num];
				if(lin.pop_calc.size() > 0) return false;
				auto &cal = lin.no_pop_calc.calc;
				it1.num = cal.size()-1;
				cal.push_back(ca);
			}
			break;
		case POPNUM: return false;
		case IE: case ONE: case ZERO: case FE: emsg_input("Eq Lin should not be"); break;
		default: emsg_input("Eq problem2"); break;
		}
		break;
		
	case PARAMVEC: case SPLINEREF: case CONSTSPLINEREF: case NUMERIC: case TIME: 
		switch(it2.type){
		case REG:
			{
				lin = lin_calc[it2.num];
				if(lin.pop_calc.size() > 0) return false;
				auto &cal = lin.no_pop_calc.calc;
				it2.num = cal.size()-1;
				cal.push_back(ca);
			}
			break;
			
		case PARAMVEC: case SPLINEREF: case CONSTSPLINEREF: case NUMERIC: case TIME:
			{		
				lin.no_pop_calc.calc.push_back(ca);
			}
			break;
		case POPNUM: return false;
		case IE: case ONE: case ZERO: case FE: emsg_input("Eq Lin should not be"); break;
		default: emsg_input("Eq problem2"); break;
		}
		break;
	case POPNUM: return false;
	case IE: case ONE: case ZERO: case FE: emsg_input("Eq Lin should not be"); break;
	default: emsg_input("Eq problem2"); break;
	}
	
	return true;
}


/// Calculates the vector of gradients against population number
vector <double> Equation::calculate_popnum_gradient_without_factor(const vector <double> &precalc) const
{
	vector <double> pop_grad;

	const auto &pop_grad_precalc = linearise.pop_grad_precalc;

	for(auto i = 0u; i < pop_grad_precalc.size(); i++){
		pop_grad.push_back(calculate_item_no_time(pop_grad_precalc[i],precalc));
	}
	
	return pop_grad;
}


/// Determines if an equation is time dependent
bool Equation::calc_time_dep(const vector <Calculation> &calc) const
{
  for(auto i = 0u; i < calc.size(); i++){
		const auto &ca = calc[i];
		
		if(item_time_dep(ca.item)) return true;
  }
	
	return false;
}


/// Determines if a vector of equation items is time dependent
bool Equation::item_time_dep(const vector <EqItem> item) const 
{
	const auto N = item.size();
		
	for(auto j = 0u; j < N; j++){
		const auto &it = item[j];
		if(it_time_dep(it)) return true;
	}
	
	return false;
}


/// Determines if item is time dependent
bool Equation::it_time_dep(const EqItem &it) const
{
	switch(it.type){
	case PARAMETER: case PARAMVEC: case IE: case ONE: case ZERO: 
	case FE:	case REG: case NUMERIC: case REG_PRECALC:
		return false;
	case SPLINE: case SPLINEREF: case CONSTSPLINEREF:
	case TIME: case REG_PRECALC_TIME:
		return true;
	case POPNUM: emsg_input("Should not have a population"); break;
	default: emsg_input("Equation error8"); break;
	}
	
	return false;
}

		
/// Calculates the value for an equation using the linearised form (for diagnostic checking)
double Equation::calculate_linearise_check(unsigned int ti, const vector <double> &popnum, const vector <double> &precalc) const 
{
	auto val = calculate_item(linearise.no_pop_precalc,ti,precalc);
	auto factor = calculate_item(linearise.factor_precalc,ti,precalc);
	for(auto j = 0u; j < pop_ref.size(); j++){
		val += rectify(popnum[pop_ref[j]])*factor*calculate_item(linearise.pop_grad_precalc[j],ti,precalc);
	}
	
	return val;
}


/// Determines if two calculations are the same
bool Equation::equal_calc(const vector <Calculation> &calc1, const vector <Calculation> &calc2) const 
{
	if(calc1.size() != calc2.size()) return false;
		
	for(auto i = 0u; i < calc1.size(); i++){
		const auto &ca1 = calc1[i]; 
		const auto &ca2 = calc2[i];
		if(ca1.op != ca2.op) return false;
		if(ca1.item.size() != ca2.item.size()) return false;
		for(auto j = 0u; j < ca1.item.size(); j++){
			const auto &it1 = ca1.item[j];
			const auto &it2 = ca2.item[j];
			
			if(it1.type != it2.type) return false;
		}
	}

	return true;
}
	
				
/// Creates a source sampler
InfSourceSampler Equation::setup_source_sampler(unsigned int ti, const vector <double> &popnum, const PV &param_val) const 
{
	InfSourceSampler ss;
	auto &val_store = ss.val_store;
	auto &val_sum_store = ss.val_sum_store;
	
	auto Npop = pop_ref.size();

	const auto &precalc = param_val.precalc;

	auto factor = 1.0;  emsg("fa");// calculate_mult(linearise.factor,ti,precalc);
	
	auto val_sum = 0.0;
	for(auto j = 0u; j < Npop; j++){
		auto po = pop_ref[j];
		auto val = popnum[po]*factor; emsg("H");// *calculate_mult(linearise.pop_grad_calc[j],ti,precalc);
		
		val_sum += val;
		val_store.push_back(val);
		val_sum_store.push_back(val_sum);
	}

	{
		auto val = calculate_no_pop(ti,precalc);
		val_sum += val;
		
		val_store.push_back(val);
		val_sum_store.push_back(val_sum);
	}
	
	ss.sum = val_sum;
		
	return ss;
}


/// Samples from the source sampler
unsigned int InfSourceSampler::sample_inf_source() const
{
	if(sum == 0) return UNSET;
	
	auto N = val_sum_store.size();
	auto z = ran()*sum;
			
	auto j = 0u; while(j < N && z > val_sum_store[j]) j++;
	if(j == N) emsg("Select prob");
	
	return j;
}


/// Samples from the source sampler
double InfSourceSampler::prob_inf_source(unsigned int j) const
{
	if(j == UNSET) return -LARGE;
	return log(val_store[j]/sum);
}

/// Gets the population reference from the population
void Linearise::init_pop_ref_from_po(const vector <unsigned int> &pop_ref)
{
	pop_ref_from_po.resize(EQ_POP_REF);
	
	for(auto i = 0u; i < pop_ref.size(); i++){
		auto po = pop_ref[i];
		PopRefFromPo pr; pr.i = i; pr.po = po;
		pop_ref_from_po[po%EQ_POP_REF].push_back(pr);
	}
}


/// Gets a population ref from a population number
unsigned int Linearise::get_pop_ref(unsigned int po) const 
{
	for(const auto &pr : pop_ref_from_po[po%EQ_POP_REF]){
		if(pr.po == po) return pr.i;
	}

	return UNSET;
}


/// Tries to take any common factors from pop_grad_calc
void Equation::get_pop_grad_calc_factorise()
{	
	auto &pgc = linearise.pop_grad_calc_store;
	
	Calculation fac_ca;
	fac_ca.op = MULTIPLY;
	
	if(pgc.size() == 0){
		if(fac_ca.item.size() == 0){
			EqItem it; it.type = ONE; it.num = UNSET; fac_ca.item.push_back(it);
		}
		linearise.factor_calc.push_back(fac_ca);
		return;
	}
	
	auto j = 0u;
	while(j < pgc[0][pgc[0].size()-1].item.size()){
		const auto &it_fi = pgc[0][pgc[0].size()-1].item[j];
		
		switch(it_fi.type){
		case ONE: case ZERO: case REG_PRECALC: case REG_PRECALC_TIME: 
		case NUMERIC: case TIME: case CONSTSPLINEREF: 
			{
				auto fl = false;
				
				vector <unsigned int> list;
				for(auto i = 0u; i < pgc.size(); i++){
					const auto &item = pgc[i][pgc[i].size()-1].item;
					
					unsigned int k;
					for(k = 0; k < item.size(); k++){
						const auto &it = item[k];
						if(it_fi.type == it.type && it_fi.num == it.num){
							list.push_back(k);
							break;
						}
					}
					if(k == item.size()){ fl = true; break;}
				}
				
				if(fl == false){
					fac_ca.item.push_back(it_fi);
					for(auto i = 0u; i < pgc.size(); i++){
						auto &ite = pgc[i][pgc[i].size()-1].item;
						
						auto k = list[i];
						if(k >= ite.size()) emsg("wrong");
						if(k+1 < ite.size()){
							ite[k] = ite[ite.size()-1];
						}
						ite.pop_back();
					}
				}
				else j++;
			}
			break;

		default:
			j++;
			break;
		}
	}
	
	if(fac_ca.item.size() == 0){
		EqItem it; it.type = ONE; it.num = UNSET; fac_ca.item.push_back(it);
	}
	linearise.factor_calc.push_back(fac_ca);
	
	if(false){
		//print_calculation();
	
		cout << " AFTE FACTOR" << endl;
		for(auto k = 0u; k < pgc.size(); k++){
			print_calc(pop[pop_ref[k]].name,pgc[k]);	
		}
		print_calc("factor",linearise.factor_calc);
		print_calc("no pop",linearise.no_pop_calc_store);
	}
	//emsg("jj");
}


/// References precalculation in linearisation
void Equation::set_precalc()
{	
	linearise.no_pop_precalc = get_precalc(linearise.no_pop_calc_store);
	
	for(auto i = 0u; i < linearise.pop_grad_calc_store.size(); i++){
		linearise.pop_grad_precalc.push_back(get_precalc(linearise.pop_grad_calc_store[i]));
	}

	linearise.factor_precalc = get_precalc(linearise.factor_calc);

	linearise.pop_grad_time_dep = false;
	for(const auto &it : linearise.pop_grad_precalc){
		if(it_time_dep(it)){ linearise.pop_grad_time_dep = true; break;}
	}
	
	linearise.factor_time_dep = it_time_dep(linearise.factor_precalc);
	linearise.no_pop_calc_time_dep = it_time_dep(linearise.no_pop_precalc);
		
	if(false){
		cout << "LINEARISATION" << endl;
		print_calculation();
		cout << "no_pop: "; print_item(linearise.no_pop_precalc); cout << endl;
		cout << "factor: "; print_item(linearise.factor_precalc); cout << endl;
		for(const auto &it : linearise.pop_grad_precalc){
			cout << "grad: "; print_item(it); cout << endl;
		}
	}
	
	// Deletes calculations as they are no longer required
	linearise.no_pop_calc_store.clear();
	linearise.pop_grad_calc_store.clear();
	linearise.factor_calc.clear();
}


/// Converts from a calculation to a precalc
EqItem Equation::get_precalc(const vector <Calculation> &calc) const
{
	//print_calc("pc",calc);
	if(calc.size() != 1) emsg("precalc should be one");
	const auto &ca = calc[0];
	
	if(ca.op != ADD && ca.op != MULTIPLY) emsg("precalc op");
	if(ca.item.size() != 1) emsg("precalc item wrong");
	
	const auto &it = ca.item[0];
	switch(it.type){
	case REG_PRECALC: case REG_PRECALC_TIME: case ONE: case ZERO: 
	case NUMERIC: case TIME: case CONSTSPLINEREF:  
		break;
	default: 
		emsg("precalc should be reg"); 
		break;
	}
	
	return it;
}

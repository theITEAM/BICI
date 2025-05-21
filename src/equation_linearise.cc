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
	linearise.on = false;

	auto ncalc = calc.size();

	auto pl = false;
	
	vector <LinearCalculation> lin_calc;
	for(auto i = 0u; i < ncalc; i++){
		const auto &ca = calc[i];
		
		if(ca.item.size() == 0) emsg("zero item problem");
		
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
						emsg("fin");
					}
				}
			}
			break;
			
		case DIVIDE:
			{
				lin = convert_to_linear_calculation(ca.item[0],MULTIPLY,lin_calc);
				auto lin_div = convert_to_linear_calculation(ca.item[1],MULTIPLY,lin_calc);
				if(lin_div.pop_calc.size() > 0) return;
				if(lin_div.no_pop_calc.calc.size() == 0) emsg("Problem with div");
				
				calc_div(lin.no_pop_calc.calc,lin_div.no_pop_calc.calc);
				for(auto j = 0u; j < lin.pop_calc.size(); j++){
					calc_div(lin.pop_calc[j].calc,lin_div.no_pop_calc.calc);
				}
				
				if(false){
					print_linear_calc("Lin after ",lin);
					emsg("fin");
				}
			}
			break;
		
		// Single parameter functions
		case EXPFUNC: case SINFUNC: case COSFUNC: case LOGFUNC:  
		case STEPFUNC: case ABSFUNC: case SQRTFUNC:
			single_param_func(ca,lin,lin_calc);
			break;
			
		// Two parameter functions
		case POWERFUNC: case THRESHFUNC: case UBOUNDFUNC: case MAXFUNC: case MINFUNC:
			two_param_func(ca,lin,lin_calc);
			break;
	
		case TAKE: emsg("SHould not have take"); break;
		
		default: emsg("Eq problem10"); break;
		}
		
		if(pl){
			cout << endl;
			print_linear_calc("LIN FINAL "+to_string(i),lin);
			cout << endl;
		}
		
		lin_calc.push_back(lin);
	}

	//if(true){
	if(pl){
		for(const auto &lin : lin_calc){
			print_linear_calc("LIN FINAL ",lin);
		}
		emsg("do");
	}

	auto lc_final = convert_to_linear_calculation(ans,ADD,lin_calc);
	
	linearise.no_pop_calc = lc_final.no_pop_calc.calc;
	linearise.no_pop_calc_time_dep = calc_time_dep(linearise.no_pop_calc);
	
	if(false) print_linear_final();
	
	auto Npop = pop_ref.size();
	if(Npop != lc_final.pop_calc.size()){
		emsg("Pop_calc problem");
	}
	
	// Orders gradients in the same way as pop_ref
	linearise.pop_grad_calc_time_dep = false;
	for(auto i = 0u; i < Npop; i++){
		auto po = pop_ref[i];
		
		auto k = 0u; while(k < Npop && lc_final.pop_calc[k].po != po) k++;
		if(k == Npop) emsg("Could not find population");
		
		const auto &calc = lc_final.pop_calc[k].calc; 
		if(calc_time_dep(calc) == true) linearise.pop_grad_calc_time_dep = true;
		linearise.pop_grad_calc.push_back(calc);
	}
	
	linearise.multi_source = false;
	
	auto num = linearise.pop_grad_calc.size();
	if(linearise.no_pop_calc.size() > 0) num++;
	if(num > 1) linearise.multi_source = true; 
	
	linearise.on = true;

	if(linearise.pop_grad_calc.size() != pop_ref.size()) emsg("Population number does not agree");

	if(false){
		cout << endl << endl << te << ":" << endl;
		print_calc("No pop dep",linearise.no_pop_calc);

		for(auto k = 0u; k < linearise.pop_grad_calc.size(); k++){
			print_calc(pop[pop_ref[k]].name,linearise.pop_grad_calc[k]);	
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
	print_calc("No pop dep",linearise.no_pop_calc);

	for(auto k = 0u; k < linearise.pop_grad_calc.size(); k++){
		print_calc(pop[pop_ref[k]].name,linearise.pop_grad_calc[k]);	
	}
}
	

/// Converts from a item to a linear calculation
LinearCalculation Equation::convert_to_linear_calculation(const EqItem &it, EqItemType op, const vector <LinearCalculation> &lin_calc) const
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
				Calculation ca; ca.op = op; ca.reg_store = 0;
				EqItem it2; it2.type = NUMERIC; it2.constant = 1;
				ca.item.push_back(it2);
				
				PopCalculation pop_c;
				pop_c.po = it.num;
				pop_c.calc.push_back(ca);
				lin.pop_calc.push_back(pop_c);
			}
			break;
		
		case PARAMETER: case SPLINE: case NUMERIC: case TIME: 
			{
				Calculation ca; ca.op = op; ca.item.push_back(it); ca.reg_store = 0;
				lin.no_pop_calc.calc.push_back(ca);
			}
			break;
		
		case IE: case FE: case ONE: emsg("Should not be here7"); break;
				
		default: emsg("Equation error11"); break;
	}
	
	return lin;
}


/// Multiplies two calculation
void Equation::calc_mult(vector <Calculation> &calc, const vector <Calculation> &calc2) const 
{
	if(false){
		cout << "CALC MULT" << endl;
		print_calc("calc",calc);
		print_calc("calc2",calc2);
	}
	
	if(calc.size() == 0) return;
	if(calc2.size() == 0){ calc.resize(0); return;}
	
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
								else ca_last.item[k].constant *= it.constant;
							}
							else ca_last.item.push_back(it);
						}
						return;
					}
					else{
						if(ca.item.size() == 1 && ca.item[0].type == NUMERIC){
							auto con = ca.item[0].constant;
							calc = calc2;
							Calculation ca; ca.op = MULTIPLY; ca.reg_store = calc.size();
							EqItem it1; it1.type = REG; it1.num = calc.size()-1;
							EqItem it2; it2.type = NUMERIC; it2.constant = con ;
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
				if(ca2.item.size() == 1 && ca2.item[0].type == NUMERIC && ca2.item[0].constant == 1){
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
								else ca_last.item[k].constant *= it2.constant;
							}
							else ca_last.item.push_back(it2);
						}
						return;
					}
					else{
						if(ca2.item.size() == 1 && ca2.item[0].type == NUMERIC){
							auto con = ca2.item[0].constant;
								
							Calculation ca; ca.op = MULTIPLY; ca.reg_store = calc.size();
							EqItem it1; it1.type = REG; it1.num = calc.size()-1;
							EqItem it2; it2.type = NUMERIC; it2.constant = con ;
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
		calc[i].reg_store += num;
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
		Calculation ca; ca.op = MULTIPLY; ca.reg_store = calc.size();
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
		emsg("here add");
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
		calc[i].reg_store += num;
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
		Calculation ca; ca.op = ADD; ca.reg_store = calc.size();
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
	if(ca.item.size() == 1 && ca.item[0].type == NUMERIC && ca.item[0].constant == 1){
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
		emsg("here div");
	}
	
	if(calc2.size() == 0) emsg("Problem");
	if(calc.size() == 0) return;
		
	auto num = calc.size();
	
	for(auto ca : calc2){
		for(auto &it : ca.item){
			if(it.type == REG) it.num += num;
		}
		ca.reg_store += num;
		
		calc.push_back(ca);
	}
	
	Calculation ca; ca.op = DIVIDE; ca.reg_store = calc.size();
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
	for(const auto &ca : calc){ print_ca(ca); cout << ",  ";}
	cout << endl;
}


/// Finds the linear calculation for a single parameter function
void Equation::single_param_func(Calculation ca, LinearCalculation &lin, const vector <LinearCalculation> &lin_calc) const
{
	const auto &it = ca.item[0];
	
	switch(it.type){
	case REG:
		{
			lin = lin_calc[it.num];
			if(lin.pop_calc.size() > 0) return;
			auto &cal = lin.no_pop_calc.calc;
			ca.reg_store = cal.size();
			ca.item[0].num = cal.size()-1;
			cal.push_back(ca);
		
			if(false) print_linear_calc("single func", lin);
		}
		break;
		
	case PARAMETER: case SPLINE: case NUMERIC: case TIME: 
		ca.reg_store = 0;
		lin.no_pop_calc.calc.push_back(ca);
		break;
		
	case POPNUM: return;
	case IE: case ONE: case FE: emsg("Eq Lin should not be"); break;
	default: emsg("Eq problem2"); break;
	}
}


/// Finds the linear calculation for a two parameter function
void Equation::two_param_func(Calculation ca, LinearCalculation &lin, const vector <LinearCalculation> &lin_calc) const
{
	auto &it1 = ca.item[0];
	auto &it2 = ca.item[1];
	
	switch(it1.type){
	case REG:
		switch(it2.type){
		case REG:
			{		
				lin = lin_calc[it1.num];
				if(lin.pop_calc.size() > 0) return;
				
				const auto &lin2 = lin_calc[it2.num];
				if(lin2.pop_calc.size() > 0) return;
				
				auto &calc = lin.no_pop_calc.calc;
				
				auto num = calc.size();
				
				const auto &calc2 = lin2.no_pop_calc.calc;
				for(auto ca_add : calc2){
					for(auto &it : ca_add.item){
						if(it.type == REG) it.num += num;
					}
					ca_add.reg_store += num;
					calc.push_back(ca_add);
				}
				
				it1.num = num-1; 
				it2.num = calc.size()-1;
				ca.reg_store = calc.size();
				calc.push_back(ca);
				if(false){
					print_linear_calc("calc",lin); 
					emsg("K");
				}
			}
			break;
			
		case PARAMETER: case SPLINE: case NUMERIC: case TIME:
			{		
				lin = lin_calc[it1.num];
				if(lin.pop_calc.size() > 0) return;
				auto &cal = lin.no_pop_calc.calc;
				it1.num = cal.size()-1;
				ca.reg_store = cal.size();
				cal.push_back(ca);
			}
			break;
		case POPNUM: return;
		case IE: case ONE: case FE: emsg("Eq Lin should not be"); break;
		default: emsg("Eq problem2"); break;
		}
		break;
		
	case PARAMETER: case SPLINE: case NUMERIC: case TIME: 
		switch(it2.type){
		case REG:
			{
				lin = lin_calc[it2.num];
				if(lin.pop_calc.size() > 0) return;
				auto &cal = lin.no_pop_calc.calc;
				it2.num = cal.size()-1;
				ca.reg_store = cal.size();
				cal.push_back(ca);
			}
			break;
			
		case PARAMETER: case SPLINE: case NUMERIC: case TIME:
			{		
				ca.reg_store = 0;
				lin.no_pop_calc.calc.push_back(ca);
			}
			break;
		case POPNUM: return;
		case IE: case ONE: case FE: emsg("Eq Lin should not be"); break;
		default: emsg("Eq problem2"); break;
		}
		break;
	case POPNUM: return;
	case IE: case ONE: case FE: emsg("Eq Lin should not be"); break;
	default: emsg("Eq problem2"); break;
	}
}


/// Calculates the vector of gradients against population number
vector <double> Equation::calculate_popnum_gradient(const vector <double> &param_val) const
{
	const auto &pop_grad_calc = linearise.pop_grad_calc;
	
	vector <double> pop_grad;
	
	for(auto i = 0u; i < pop_grad_calc.size(); i++){
		pop_grad.push_back(calculate_calculation_notime(pop_grad_calc[i],param_val));
	}
	
	return pop_grad;
}


/// Calculates the term with no population
double Equation::calculate_no_pop(unsigned int ti, const vector <double> &param_val, const vector <SplineValue> &spline_val) const
{
	return calculate_calculation(linearise.no_pop_calc,ti,param_val,spline_val);
}


/// Determines if an equation is time dependent
bool Equation::calc_time_dep(const vector <Calculation> &calc) const
{
  for(auto i = 0u; i < calc.size(); i++){
		const auto &ca = calc[i];
		
		const auto &item = ca.item;
		const auto N = item.size();
		
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
			
			switch(it.type){
				case PARAMETER: case IE: case ONE: case FE:	case REG: case NUMERIC: break;
				case SPLINE: case TIME: return true;
				case POPNUM: emsg("Should not have a population"); break;
				default: emsg("Equation error"); break;
			}
		}
  }
	
	return false;
}


/// Calculates based on calc
double Equation::calculate_calculation(const vector <Calculation> &calc, unsigned int ti, const vector <double> &param_val, const vector <SplineValue> &spline_val) const 
{
	auto imax = calc.size();
	if(imax == 0) return 0;
	
	vector <double> regcalc(imax);

  for(auto i = 0u; i < imax; i++){
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
				case POPNUM: emsg("Should not have a population"); break;
				case REG: num[j] = regcalc[it.num]; break;
				case NUMERIC: num[j] = it.constant; break;
				case TIME: num[j] = timepoint[ti]; break;
				default: emsg("Equation error"); break;
			}
		}

    regcalc[ca.reg_store] = calculate_operation(ca.op,num);
  }
	
	return regcalc[imax-1];
}

/// Calculates based on calc (using spline store if that is set)
double Equation::calculate_calculation_spline_store(const vector <Calculation> &calc, unsigned int ti, const vector <double> &param_val, const vector <SplineValue> &spline_val) const 
{
	auto imax = calc.size();
	if(imax == 0) return 0;
	
	vector <double> regcalc(imax);

  for(auto i = 0u; i < imax; i++){
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
					{
						const auto &sv =  spline_val[it.num];
						num[j] = sv.store[ti];
						if(num[j] == UNSET) num[j] = sv.val[ti];
						if(it.index != UNSET) emsg("SHould be unset4"); 
					}
					break;
				case IE: emsg("Should not include ind effect"); break;
				case ONE: num[j] = 1; break;
				case FE: emsg("Should not include fixed effect"); break;
				case POPNUM: emsg("Should not have a population"); break;
				case REG: num[j] = regcalc[it.num]; break;
				case NUMERIC: num[j] = it.constant; break;
				case TIME: num[j] = timepoint[ti]; break;
				default: emsg("Equation error"); break;
			}
		}

    regcalc[ca.reg_store] = calculate_operation(ca.op,num);
  }
	
	return regcalc[imax-1];
}


/// Calculates based on calc (no time variation)
double Equation::calculate_calculation_notime(const vector <Calculation> &calc, const vector <double> &param_val) const 
{
	auto imax = calc.size();
	if(imax == 0) return 0;
	
	vector <double> regcalc(imax);

  for(auto i = 0u; i < imax; i++){
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
				case SPLINE: emsg("Should not include spline"); break;
				case IE: emsg("Should not include ind effect"); break;
				case ONE: num[j] = 1; break;
				case FE: emsg("Should not include fixed effect"); break;
				case POPNUM: emsg("Should not have a population"); break;
				case REG: num[j] = regcalc[it.num]; break;
				case NUMERIC: num[j] = it.constant; break;
				case TIME: emsg("Should not have time"); break;
				default: emsg("Equation error"); break;
			}
		}

    regcalc[ca.reg_store] = calculate_operation(ca.op,num);
  }
	
	return regcalc[imax-1];
}


/// Calculates the value for an equation using the linearised form (for diagnostic checking)
double Equation::calculate_linearise_check(unsigned int ti, const vector <double> &popnum, const vector <double> &param_val, const vector <SplineValue> &spline_val) const 
{
	auto val = 0.0;
	
	const auto &calc = linearise.no_pop_calc;
	if(calc.size() > 0) val += calculate_calculation(calc,ti,param_val,spline_val);
	
	for(auto j = 0u; j < pop_ref.size(); j++){
		const auto &pc = linearise.pop_grad_calc[j];
		val += rectify(popnum[pop_ref[j]])*calculate_calculation(pc,ti,param_val,spline_val);
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
		if(ca1.reg_store != ca2.reg_store) return false;
		if(ca1.item.size() != ca2.item.size()) return false;
		for(auto j = 0u; j < ca1.item.size(); j++){
			const auto &it1 = ca1.item[j];
			const auto &it2 = ca2.item[j];
			
			if(it1.type != it2.type) return false;
			if(it1.num != it2.num) return false;
			if(it1.index != it2.index) return false;
			if(it1.constant != it2.constant) return false;
		}
	}

	return true;
}
	


/// Creates a source sampler
InfSourceSampler Equation::setup_source_sampler(unsigned int ti, const vector <double> &popnum, const vector <double> &param_val, const vector <SplineValue> &spline_val) const 
{
	InfSourceSampler ss;
	auto &val_store = ss.val_store;
	auto &val_sum_store = ss.val_sum_store;
	
	auto Npop = pop_ref.size();

	auto val_sum = 0.0;
	for(auto j = 0u; j < Npop; j++){
		auto po = pop_ref[j];
		
		auto val = popnum[po]*(calculate_calculation(linearise.pop_grad_calc[j],ti,param_val,spline_val));
		
		val_sum += val;
		val_store.push_back(val);
		val_sum_store.push_back(val_sum);
	}

	{
		auto val = calculate_calculation(linearise.no_pop_calc,ti,param_val,spline_val);
		val_sum += val;
		
		val_store.push_back(val);
		val_sum_store.push_back(val_sum);
	}
	
	if(val_sum == 0) emsg("Problem selecting");
	
	return ss;
}


/// Samples from the source sampler
unsigned int InfSourceSampler::sample_inf_source() const
{
	auto N = val_sum_store.size();
	auto z = ran()*val_sum_store[N-1];
			
	auto j = 0u; while(j < N && z > val_sum_store[j]) j++;
	if(j == N) emsg("Select prob");
	
	return j;
}


/// Samples from the source sampler
double InfSourceSampler::prob_inf_source(unsigned int j) const
{
	return log(val_store[j]/val_sum_store[val_sum_store.size()-1]);
}

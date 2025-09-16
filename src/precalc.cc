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

#include "precalc.hh"
#include "utils.hh"

Precalc::Precalc(const vector <SpeciesSimp> &species, const vector <Spline> &spline, const vector <ParamVecEle> &param_vec, const vector <Population> &pop, Constant &constant, const vector <double> &timepoint, const Details &details) : species(species), spline(spline), param_vec(param_vec), pop(pop), constant(constant), timepoint(timepoint), details(details)
{
}


/// Adds an equation onto the precalculated equation calculation
bool Precalc::add_eqn(vector <Calculation> &calc, const vector <double> &param_vec_ref, const vector <double> &spline_ref, vector <unsigned int> &list_precalc, PrecalcAddType add_type)
{
	vector <unsigned int> become_Rrecalc(calc.size(),UNSET);

	// Removes any reference to parameter and spline 
	for(auto &ca : calc){
		for(auto &it : ca.item){
			switch(it.type){
			case PARAMVEC: 
				it.type = REG_PRECALC; 
				it.num = param_vec_ref[it.num];
				//list_precalc.push_back(it.num);
				break;
				
			case SPLINEREF: 
				{
					auto &spl = spline[it.num];
					if(spl.constant == true){
						it.type = CONSTSPLINEREF;
					}
					else{
						it.type = REG_PRECALC_TIME; 
						it.num = spline_ref[it.num];
						//list_precalc.push_back(it.num);
					}
				}
				break;
			
			default:
				break;
			}
		}
	}
	
	auto stop_combine_fl = false;
	
	if(add_type == PRECALC_PARAM_ONLY) return stop_combine_fl;
	
	bool fl;
	do{
		for(auto i = 0u; i < calc.size(); i++){
			if(become_Rrecalc[i] == UNSET){
				auto &ca = calc[i];
				
				auto &item = ca.item;
				
				auto nregderpop = 0u, npop = 0u;
				for(auto j = 0u; j < item.size(); j++){
					switch(item[j].type){
					case REG: case DERIVE: case INTEGRAL: nregderpop++; break;
					case POPNUM: nregderpop++; npop++; break;
					default: break;
					}
				}
			
				// This replaces calculation with that in non-pop 
				if(nregderpop == 0 && i+1 != calc.size()){  	
					PreCalc pc; pc.op = ca.op; pc.item = ca.item;
					auto it = add(pc,list_precalc);
					become_Rrecalc[i] = it.num;
				}
				else{ // Attempts to put part of calculation into non-pop
					if(ca.op == DIVIDE){ // If dividing by something not 
						auto &denom = item[1];
						if(denom.type != REG && denom.type != POPNUM && 
						   denom.type != DERIVE && denom.type != INTEGRAL){
							ca.op = MULTIPLY;
							
							PreCalc pcalc;
							pcalc.op = DIVIDE;
								
							EqItem one; one.type = ONE;
							pcalc.item.push_back(one);	
							pcalc.item.push_back(denom);
							
							denom = add(pcalc,list_precalc);
						}
					}
					
					auto flag = false;
					
					if(ca.op == MULTIPLY && i+1 < calc.size() && calc[i+1].op == DIVIDE){
						if(item.size() - nregderpop >= 1){
							auto &ca_next = calc[i+1];
							const auto &itdiv = ca_next.item;
							if(itdiv[0].type == REG && itdiv[0].num == i){
								auto &denom = itdiv[1];
								if(denom.type != REG && denom.type != POPNUM && 
								   denom.type != DERIVE && denom.type != INTEGRAL){
									flag = true;
								
									PreCalc pcalc;
									pcalc.op = DIVIDE;
								
									EqItem one; one.type = ONE;
									pcalc.item.push_back(one);
									pcalc.item.push_back(denom);
									
									auto inew = add(pcalc,list_precalc);
									
									PreCalc pcalc2;
									pcalc2.op = MULTIPLY;
							
									vector <EqItem> item_new;
							
									for(auto j = 0u; j < item.size(); j++){
										const auto &it = item[j];
										switch(it.type){
										case REG: case POPNUM: case DERIVE: case INTEGRAL: item_new.push_back(it); break;
										default: pcalc2.item.push_back(it); break;
										}
									}
									pcalc2.item.push_back(inew);
								
									auto inew2 = add(pcalc2,list_precalc);
									item_new.push_back(inew2);
									
									ca_next.op = MULTIPLY;
									ca_next.item = item_new;
									
									become_Rrecalc[i] = USINT_MAX;
								}
							}
						}
					}
					
					if(flag == false && (ca.op == MULTIPLY || ca.op == ADD)){
						if(item.size() - nregderpop >= 2){
							if(add_type == PRECALC_STOP_COMBINE_MULT && npop > 0){ // Stops being combined so we can get factor
								stop_combine_fl = true;
							}
							else{
								PreCalc pcalc;
								pcalc.op = ca.op;
								
								vector <EqItem> item_new;
								
								for(auto j = 0u; j < item.size(); j++){
									const auto &it = item[j];
									switch(it.type){
									case REG: case POPNUM: case DERIVE: case INTEGRAL: item_new.push_back(it); break;
									default: pcalc.item.push_back(it); break;
									}
								}
						
								auto inew = add(pcalc,list_precalc);
								item_new.push_back(inew);
								
								ca.item = item_new;
							}
						}
					}
				}
			}
		}
	
		fl = false;
		// Replaces register with non_popcalc
		
		for(auto i = 0u; i < calc.size(); i++){
			if(become_Rrecalc[i] == UNSET){
				for(auto &it : calc[i].item){
					if(it.type == REG && become_Rrecalc[it.num] != UNSET){
						it.type = REG_PRECALC; 
						it.num = become_Rrecalc[it.num];
						if(calcu[it.num].time_dep) it.type = REG_PRECALC_TIME; 
						if(it.num == USINT_MAX) emsg("Problem used");
						fl = true;
					}
				}
			}
		}
	}while(fl == true);
	
	// Removes become_Rfac

	auto calc_old = calc;

	vector <unsigned int> map(calc_old.size(),UNSET);
	
	calc.clear();
	for(auto i = 0u; i < calc_old.size(); i++){
		if(become_Rrecalc[i] == UNSET){
			map[i] = calc.size();
			calc.push_back(calc_old[i]);
		}
	}	
	
	for(auto &ca : calc){
		for(auto &it : ca.item){
			if(it.type == REG){
				it.num = map[it.num];
				if(it.num == UNSET) emsg("Problem with map");
			}
		}
	}
	
	return stop_combine_fl;
}


/// Adds a parameter to precalculation
unsigned int Precalc::add_param(unsigned int th)
{
	PreCalc ca;
	ca.op = SINGLE;
	ca.time_dep = false;
	
	EqItem it;
	it.type = PARAMVEC;
	it.num = th;
	ca.item.push_back(it);

	auto vec = get_vec(ca);
	
	auto j = hash_ca.existing(vec);
	if(j == UNSET){
		j = calcu.size();
		//list_precalc.push_back(j);
		
		calcu.push_back(ca);
		hash_ca.add(j,vec);
	}

	return j;
}


/// Adds a spline to precalculation
unsigned int Precalc::add_spline(unsigned int s, vector <unsigned int> &list_precalc)
{
	PreCalc ca;
	ca.op = SINGLE;
	ca.time_dep = true;
	
	EqItem it;
	it.type = SPLINEREF;
	it.num = s;
	ca.item.push_back(it);

	const auto &spl = spline[s];
	for(auto i = 0u; i < spl.param_ref.size(); i++){
		EqItem it;
		it.type = PARAMVEC;
		it.num = spl.param_ref[i];
		ca.item.push_back(it);
	}
	auto vec = get_vec(ca);
	
	auto j = hash_ca.existing(vec);
	if(j == UNSET){
		j = calcu.size();
		list_precalc.push_back(j);
		
		calcu.push_back(ca);
		hash_ca.add(j,vec);
		
		for(auto ti = 1u; ti < details.T; ti++){
			PreCalc ca;	
			ca.op = NOOP;
			ca.time_dep = true;
			calcu.push_back(ca);
		}
	}
	
	return j;
}


/// Adds a calculation 
EqItem Precalc::add(PreCalc &pcalc, vector <unsigned int> &list_precalc)
{
	auto vec = get_vec(pcalc);
	
	// Determines if time dependent
	pcalc.time_dep = false;
	
	for(auto &it : pcalc.item){
		switch(it.type){
		case SPLINEREF: case CONSTSPLINEREF: case POPNUM: case TIME: case REG_PRECALC_TIME:
			pcalc.time_dep = true; 
			break;
		default: break;
		}
	}
						
	auto j = hash_ca.existing(vec);
	if(j == UNSET){
		j = calcu.size();
		list_precalc.push_back(j);
		
		hash_ca.add(j,vec);
		
		calcu.push_back(pcalc);	
		
		if(pcalc.time_dep == true){
			PreCalc ca2; ca2.op = NOOP; ca2.time_dep = true;
			auto T = details.T;
			for(auto ti = 1u; ti < T; ti++) calcu.push_back(ca2);
		}
	}	
	
	EqItem inew;
	inew.type = REG_PRECALC; if(pcalc.time_dep == true) inew.type = REG_PRECALC_TIME;
	inew.num = j;
	
	return inew;
}


/// Gets a hash vector for a calulation item
vector <unsigned int> Precalc::get_vec(const PreCalc &ca) const 
{
	const auto &item = ca.item;
	
	vector <unsigned int> vec;
		
	vec.push_back(ca.op);
	for(auto j = 0u; j < item.size(); j++){
		vec.push_back(item[j].type);
		vec.push_back(item[j].num);
	}

	return vec;
}	


/// Prints steps used for a calculation
void Precalc::print_calc() const   
{
	cout << "PRE-CALCULATION:" << endl;
	
	auto T = details.T;
	auto imax = calcu.size(); 
	if(imax > 1000) imax = 1000;
  for(auto i = 0u; i < imax; i++){
		const auto &ca = calcu[i];
		print_ca(i,ca);
		cout << endl;
		
		if(ca.time_dep == true) i += T-1;
  }
	
	if(imax < calcu.size()) cout << "...";

  cout << endl << endl;
}


/// Calculates the entire precacluation
void Precalc::calculate_all(const vector <unsigned int> &list_recalc, PV &param_val) const 
{
	vector <unsigned int> list_time;
	for(auto ti = 0u; ti < details.T; ti++) list_time.push_back(ti);
	
	calculate(list_recalc,list_time,param_val,false);
}


/// Calculates the value for an equation
void Precalc::calculate(const vector <unsigned int> &list_recalc, const vector <unsigned int> &list_time, PV &param_val, bool store) const 
{
	const auto &value = param_val.value;
	const auto &cval = constant.value;
	auto &precalc = param_val.precalc;
	
	for(auto i : list_recalc){
		const auto &ca = calcu[i];
		
		if(store){
			if(ca.time_dep){ 
				for(auto ti : list_time) param_val.precalc_change(i+ti);
			}
			else{
				param_val.precalc_change(i);
			}
		}
		
		if(ca.op == SINGLE){
			const auto &it = ca.item[0];
			switch(it.type){
			case PARAMVEC: 
				precalc[i] = value[it.num];
				break;
			
			case SPLINEREF:
				{
					const auto &spl = spline[it.num];
					
					switch(spl.type){
					case LINEAR_SPL:
						for(auto ti : list_time){
							const auto &div = spl.div[ti];
					
							auto f = div.f;
							precalc[i+ti] = value[div.th1]*f+value[div.th2]*(1-f);
						}
						break;
						
					case SQUARE_SPL:
						for(auto ti : list_time){
							const auto &div = spl.div[ti];
							precalc[i+ti] = value[div.th1];
						}
						break;
						
					case CUBICPOS_SPL: case CUBIC_SPL: 
						{
							vector <double> val;
							for(auto pr : spl.param_ref) val.push_back(value[pr]);
						
							auto cspl = solve_cubic_spline(spl.info.knot_tdiv,val,spl.type);
						
							calculate_cubic_spline_precalc(precalc,i,list_time,spl.cubic_div,cspl);
						}
						break;
					}
				}
				break;
				
			default: emsg("Should not be single"); break;
			}
		}
		else{
			const auto &item = ca.item;
			const auto N = item.size();
			
			vector <double> num(N);
			
			if(ca.time_dep){
				for(auto ti : list_time){
					for(auto j = 0u; j < N; j++){
						const auto &it = item[j];
						
						switch(it.type){
							case REG_PRECALC: num[j] = precalc[it.num]; break;
							case REG_PRECALC_TIME: num[j] = precalc[it.num+ti]; break;							
							case ONE: num[j] = 1; break;
							case NUMERIC: num[j] = constant.value[it.num]; break;
							case CONSTSPLINEREF: num[j] = spline[it.num].const_val[ti]; break;
							case TIME: num[j] = timepoint[ti]; break;
							default: eqn_type_error(it.type,6); break;
						}
					}

					precalc[i+ti] = calculate_operation(ca.op,num);
				}
			}
			else{			
				for(auto j = 0u; j < N; j++){
					const auto &it = item[j];	
					switch(it.type){	
						case ONE: num[j] = 1; break;
						case REG_PRECALC: num[j] = precalc[it.num]; break;
						case NUMERIC: num[j] = cval[it.num]; break;
						default: eqn_type_error(it.type,8); break;
					}
				}

					
				precalc[i] = calculate_operation(ca.op,num);
			}
		}
  }
}


/// Prints the calculation from a specific calculation
void Precalc::print_ca(unsigned int i, const PreCalc &ca) const
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
	
	for(auto j = 0u; j < item.size(); j++){
		const auto &it = item[j];
		switch(it.type){
		case INTEGRAL:
			emsg("Should not have integral");
			break;
			
		case PARAMETER: 
			emsg("Should not be parameter");
			break;
			
		case PARAMVEC: 
			cout << add_escape_char(param_vec[it.num].name);
			break;
		
		case SPLINE: 
			emsg("Should not be spline");
			break;
			
		case SPLINEREF:
			cout << "Spline " << add_escape_char(spline[it.num].name); 
			break;
		
		case CONSTSPLINEREF:
			cout << "Const Spline " << add_escape_char(spline[it.num].name); 
			break;
			
		case DERIVE: 
			emsg("Should not be derive");
			break;
				
		case POPNUM: 
			cout << "'" << pop[it.num].name << "'";
			break;
	
		case IE: emsg("SHould not be ie"); break;
		case ONE: cout << "1"; break;
		case FE: emsg("SHould not be fe"); break;
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
		default: cout << it.type << endl; emsg_input("Precalc Eq "); break;
		}
		
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

	switch(ca.op){
	case EXPFUNC: case SINFUNC: case COSFUNC: case LOGFUNC: case STEPFUNC: case POWERFUNC: 
	case THRESHFUNC: case UBOUNDFUNC: case MAXFUNC: case MINFUNC: case ABSFUNC: case SQRTFUNC:
	case SIGFUNC:
		cout << ")"; 
		break;
	default: break;
	}
		
	cout <<  " > ";
	cout <<  "Rpre" << i;
	
	if(ca.time_dep) cout << " Time dep";
}
 

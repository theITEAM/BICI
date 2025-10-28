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
bool Precalc::add_eqn(vector <Calculation> &calc, const vector <unsigned int> &param_vec_ref, const vector <unsigned int> &spline_ref, SpecPrecalc &spec_precalc, PrecalcAddType add_type)
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
					auto it = add(pc,spec_precalc);
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
							
							denom = add(pcalc,spec_precalc);
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
									
									auto inew = add(pcalc,spec_precalc);
									
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
								
									auto inew2 = add(pcalc2,spec_precalc);
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
						
								auto inew = add(pcalc,spec_precalc);
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
	it.num = th; if(th == UNSET) emsg("Parem vec problem");
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
unsigned int Precalc::add_spline(unsigned int s, SpecPrecalc &spec_precalc)
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
		const auto &pr = spl.param_ref[i];
		if(pr.cons){
			EqItem it;
			it.type = NUMERIC;
			it.num = pr.index;
			ca.item.push_back(it);
		}
		else{
			EqItem it;
			it.type = PARAMVEC;
			it.num = pr.index;
			ca.item.push_back(it);
		}
	}
	auto vec = get_vec(ca);
	
	auto j = hash_ca.existing(vec);
	if(j == UNSET){
		j = calcu.size();
		
		calcu.push_back(ca);
		
		sp_add(spec_precalc,j,all_time);
		
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
EqItem Precalc::add(PreCalc &pcalc, SpecPrecalc &spec_precalc)
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
		
		hash_ca.add(j,vec);
		
		calcu.push_back(pcalc);	
		
		sp_add(spec_precalc,j,all_time);
		
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
	auto imax_lim = LARGE;
	//auto imax_lim = 1000;
	if(imax > imax_lim) imax = imax_lim;
  for(auto i = 0u; i < imax; i++){
		const auto &ca = calcu[i];
		print_ca(i,ca);
		cout << endl;
		
		if(ca.time_dep == true) i += T-1;
  }
	
	if(imax < calcu.size()) cout << "...";

  cout << endl << endl;
}


/// Gets a value for a spline knot
double Precalc::get_splineval(const ElementRef &er, const vector <double> &value, const vector <double> &cval) const
{
	if(er.cons) return cval[er.index];
	return value[er.index];
}


/*
/// Works out any constants within splines
void Precalc::calc_spline_const(PV &param_val, const vector <unsigned int> &spline_ref) const
{
	for(auto i = 0u; i < spline.size(); i++){
		const auto &spl = spline[i];
		
		vector <unsigned int> list_time;
		switch(spl.type){
		case SQUARE_SPL:
			for(auto ti = 0u; ti < details.T; ti++){
				const auto &div = spl.div[ti];
				if(spl.param_ref[div.index].cons){
					list_time.push_back(ti);
				}
			}
			break;
		
		default: break;
		}
		
		if(list_time.size() > 0){
			vector <unsigned int> list_recalc;
			list_recalc.push_back(spline_ref[i]);
			calculate(list_recalc,list_time,param_val,false);
		}
	}
}
*/


/// Calulcates initial value for precalc (using any constant values_)
/// Calculates the value for an equation
vector <double> Precalc::calculate_precalc_init(const SpecPrecalc &spec_precalc) const 
{
	vector <double> precalc(calcu.size(),UNSET);
	const auto &cval = constant.value;
	
	for(const auto &in : spec_precalc.info){
		auto i = in.i;
		
		const auto &ca = calcu[i];
		
		if(ca.op == SINGLE){
			const auto &it = ca.item[0];
			switch(it.type){
			case PARAMVEC: 
				break;
			
			case SPLINEREF:
				{
					const auto &spl = spline[it.num];
					if(spl.constant){
						for(auto ti : all_time){
							precalc[i+ti] = spl.const_val[ti];
						}
					}
					else{
						switch(spl.type){
						case LINEAR_SPL:
							for(auto ti : all_time){
								const auto &div = spl.div[ti];
						
								auto f = div.f;
								auto ind = div.index;
							
								const auto &pr = spl.param_ref[ind];
								const auto &pr2 = spl.param_ref[ind];
								
								if(pr.cons && pr2.cons){
									auto val1 = cval[pr.index];
									auto val2 = cval[pr2.index]; 
									precalc[i+ti] = val1*f+val2*(1-f);
								}
							}
							break;
							
						case SQUARE_SPL:
							for(auto ti : all_time){
								const auto &div = spl.div[ti];
								const auto &pr = spl.param_ref[div.index];
								if(pr.cons) precalc[i+ti] = cval[pr.index]; 
							}
							break;
							
						case CUBICPOS_SPL: case CUBIC_SPL: 
							{
								vector <double> val;
								auto fl = false;
								for(const auto &pr : spl.param_ref){
									if(pr.cons) val.push_back(cval[pr.index]);
									else fl = true;
								}
								
								if(fl == false){
									auto cspl = solve_cubic_spline(spl.info.knot_tdiv,val,spl.type);
									calculate_cubic_spline_precalc(precalc,i,all_time,spl.cubic_div,cspl);
								}
							}
							break;
						}
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
				for(auto ti : all_time){
					auto fl = false;
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
						if(num[j] == UNSET){ fl = true; break;}
					}

					if(fl == false) precalc[i+ti] = calculate_operation(ca.op,num);
				}
			}
			else{			
				auto fl = false;
				for(auto j = 0u; j < N; j++){
					const auto &it = item[j];	
					switch(it.type){	
						case ONE: num[j] = 1; break;
						case REG_PRECALC: num[j] = precalc[it.num]; break;
						case NUMERIC: num[j] = cval[it.num]; break;
						default: eqn_type_error(it.type,8); break;
					}
					if(num[j] == UNSET){ fl = true; break;}
				}

				if(fl == false) precalc[i] = calculate_operation(ca.op,num);
			}
		}
  }
	
	return precalc;
}


/// Calculates the value for an equation
void Precalc::calculate(const SpecPrecalc &spec_calc, PV &param_val, bool store) const 
{
	const auto &info = spec_calc.info;
	if(info.size() == 0) return;
	
	const auto &value = param_val.value;
	const auto &cval = constant.value;
	auto &precalc = param_val.precalc;
	
	for(const auto &ci : info){
		auto i = ci.i;

		const auto &ca = calcu[i];
		
		if(store){
			if(ca.time_dep){ 
				for(auto ti : spec_calc.list_time[ci.tlist]) param_val.precalc_change(i+ti);
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
					if(!spl.constant){
						const auto &list_time = spec_calc.list_time[ci.tlist];
						
						switch(spl.type){
						case LINEAR_SPL:
							for(auto ti : list_time){
								const auto &div = spl.div[ti];
						
								auto f = div.f;
								auto ind = div.index;
								
								auto val1 = get_splineval(spl.param_ref[ind],value,cval);
								auto val2 = get_splineval(spl.param_ref[ind+1],value,cval);
								precalc[i+ti] = val1*f+val2*(1-f);
							}
							break;
							
						case SQUARE_SPL:
							for(auto ti : list_time){
								const auto &div = spl.div[ti];
								precalc[i+ti] = get_splineval(spl.param_ref[div.index],value,cval);
							}
							break;
							
						case CUBICPOS_SPL: case CUBIC_SPL: 
							{
								vector <double> val;
								for(const auto &pr : spl.param_ref){
									val.push_back(get_splineval(pr,value,cval));
								}
								
								auto cspl = solve_cubic_spline(spl.info.knot_tdiv,val,spl.type);
							
								calculate_cubic_spline_precalc(precalc,i,list_time,spl.cubic_div,cspl);
							}
							break;
						}
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
				const auto &list_time = spec_calc.list_time[ci.tlist];
						
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
	
	if(item.size() > 0 && (item[0].type == SPLINEREF || item[0].type == CONSTSPLINEREF)){
		const auto &it = item[0];
		switch(it.type){
		case SPLINEREF:
			cout << "Spline " << spline[it.num].name; 
			break;
		
		case CONSTSPLINEREF:
			cout << "Const Spline " << spline[it.num].name; 
			break;
			
		default:
			break;
		}
	}
	else{
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
				cout << param_vec[it.num].name;
				break;
			
			case SPLINE: 
				emsg("Should not be spline");
				break;
				
			case SPLINEREF:
				cout << "Spline " << spline[it.num].name; 
				break;
			
			case CONSTSPLINEREF:
				cout << "Const Spline " << spline[it.num].name; 
				break;
				
			case DERIVE: 
				emsg("Should not be derive");
				break;
					
			case POPNUM: 
				cout << "'" << pop[it.num].name << "'";
				break;
				
			case POPTIMENUM: 
				emsg("Should not be poptime num");
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
	}
		
	cout <<  " > ";
	cout <<  "Rpre" << i;
	
	if(ca.time_dep) cout << " Time dep";
}
 

/// Adds an element to a precalculation
void Precalc::sp_add(SpecPrecalc &sprec, unsigned int i, const vector <unsigned int> &ltime) const 
{
	auto &list_time = sprec.list_time;
			
	auto k = sprec.hash.find(i);
	if(k != UNSET){ // Already exists
		emsg("already");
		auto &in = sprec.info[k];
		
		if(calcu[i].time_dep){
			if(!equal_vec(list_time[in.tlist],ltime)){
				emsg("create new");
			}
		}
	}
	else{
		PrecalcInfo pi; pi.i = i; 
		
		if(calcu[i].time_dep){
			auto j = 0u;
			while(j < list_time.size() && !equal_vec(ltime,list_time[j])) j++;
			if(j == list_time.size()) list_time.push_back(ltime);
		
			pi.tlist = j;
		}
		else{
			pi.tlist = UNSET;
		}
		
		sprec.hash.add(sprec.info.size(),i);
		sprec.info.push_back(pi);
	}
}


/// Adds a list of elements to a precalculation
void Precalc::sp_add(SpecPrecalc &sprec, const vector <unsigned int> &i_list, const vector <unsigned int> &ltime) const
{
	if(i_list.size() == 0) return;
	
	for(auto i : i_list) sp_add(sprec,i,ltime);
	
/*
	auto &list_time = sprec.list_time;
		
	auto j = 0u;
	while(j < list_time.size() && !equal_vec(ltime,sprec.list_time[j])) j++;
	if(j == list_time.size()) list_time.push_back(ltime);

	for(auto i : i_list){
		if(calcu[i].time_dep){
			PrecalcInfo pi; pi.i = i; pi.tlist = j;
			sprec.info.push_back(pi);
		}
		else{
			PrecalcInfo pi; pi.i = i; pi.tlist = UNSET;
			sprec.info.push_back(pi);
		}
	}
	*/
}


/// Sets up a spec precalc bu combining several together
SpecPrecalc Precalc::combine_spec_precalc(const vector <unsigned int> &param_list_tot) const 
{
	auto C = param_list_tot.size();
	
	if(C == 1) return param_vec[param_list_tot[0]].spec_precalc_after;
	
	auto T = details.T;
	
	SpecPrecalc spre;

	vector <unsigned int> index(C,0);
	
	do{
		auto imin = LARGE;
		for(auto j = 0u; j < C; j++){
			const auto &info = param_vec[param_list_tot[j]].spec_precalc_after.info;
			if(index[j] < info.size()){
				auto i = info[index[j]].i;
				if(i < imin) imin = i;
			}
		}
			
		if(imin == LARGE) break;
		
		if(calcu[imin].time_dep == false){	
			sp_add(spre,imin,all_time); 
		}
		else{
			vector <bool> map_ti(T,false);
			
			for(auto j = 0u; j < C; j++){
				const auto &spa = param_vec[param_list_tot[j]].spec_precalc_after;
				const auto &in = spa.info[index[j]];
				if(in.i == imin){
					for(auto ti : spa.list_time[in.tlist]) map_ti[ti] = true;
				}
			}
			
			vector <unsigned int> ti_vec;
			for(auto ti = 0u; ti < T; ti++){
				if(map_ti[ti]) ti_vec.push_back(ti);
			}
		
			sp_add(spre,imin,ti_vec); 
		}
		
		for(auto j = 0u; j < C; j++){
			const auto &spa = param_vec[param_list_tot[j]].spec_precalc_after;
			if(spa.info[index[j]].i == imin) index[j]++;
		}
	}while(true);
		
	return spre;
}


/// This shrinks a precalculation based on only being needed for a given time range
SpecPrecalc Precalc::shrink_sprec(const vector <unsigned int> &lt, SpecPrecalc spre) const
{
	auto T = details.T;
	
	vector <bool> map(T,false);
	for(auto ti : lt) map[ti] = true;
	for(auto &listt : spre.list_time){
		vector <unsigned int> listt_new;
		for(auto ti : listt){
			if(map[ti] == true) listt_new.push_back(ti);
		}			
		listt = listt_new;
	}
	
	vector <PrecalcInfo> info_new;
	for(const auto &in : spre.info){
		if(in.tlist == UNSET) info_new.push_back(in);
		else{
			if(spre.list_time[in.tlist].size() > 0) info_new.push_back(in);
		}
	}
	
	if(info_new.size() != spre.info.size()) spre.info = info_new;
	
	return spre;
}


/// 
void Precalc::set_param(SpecPrecalc &set_param_spec_precalc, SpecPrecalc &spec_precalc_after, bool spl_fl) const
{
	// Transfers over parameter
	auto &info = spec_precalc_after.info;

	if(info.size() == 0) emsg("zero size");
	const auto &list_time = spec_precalc_after.list_time;
	set_param_spec_precalc.info.push_back(info[0]);
	info.erase(info.begin());

	// Transfers the spline (if it exists)
	if(spl_fl){
		const auto &in = info[0];
		sp_add(set_param_spec_precalc,in.i,list_time[in.tlist]);
		info.erase(info.begin());
	}
	
	/*
	// Transfers over parameter
	auto &info = spec_precalc_after.info;

	if(info.size() == 0) emsg("zero size");
	const auto &list_time = spec_precalc_after.list_time;
	set_param_spec_precalc.info.push_back(info[0]);
	//info.erase(info.begin());

	// Transfers the spline (if it exists)
	if(spl_fl){
		const auto &in = info[1];
		sp_add(set_param_spec_precalc,in.i,list_time[in.tlist]);
		//info.erase(info.begin());
	}
	*/
}


// Sets precalculation to be done after sampling 
SpecPrecalc Precalc::calculate_spec_precalc_sample(const SpecPrecalc &spec_precalc) const
{
	auto C = calcu.size();
	vector <bool> map(C,false);
	
	for(auto &in : spec_precalc.info) map[in.i] = true;
	
	for(const auto &pv : param_vec){
		
		for(auto &in : pv.spec_precalc_before.info) map[in.i] = true;
		
		if(!pv.reparam_time_dep){
			for(auto &in : pv.spec_precalc_after.info) map[in.i] = true;
		}
	}
	
	SpecPrecalc spec;
	spec.list_time.push_back(all_time);
	
	for(auto i = 0u; i < C; i++){
		if(map[i]){
			PrecalcInfo pi; pi.i = i;
			if(calcu[i].time_dep) pi.tlist = 0;
			else pi.tlist = UNSET;	
			spec.info.push_back(pi);			
		}
	}
	
	return spec;
}


// Sets all precalculation to be done (apart from derived)
SpecPrecalc Precalc::calculate_spec_precalc_all(const SpecPrecalc &spec_precalc, const vector <SpecPrecalcTime> &spec_precalc_time) const
{
	auto C = calcu.size();
	vector <bool> map(C,false);
	
	for(auto &in : spec_precalc.info) map[in.i] = true;
	
	for(const auto &pv : param_vec){
		for(auto &in : pv.spec_precalc_before.info) map[in.i] = true;
		
		if(!pv.reparam_time_dep){
			for(auto &in : pv.set_param_spec_precalc.info) map[in.i] = true;
			for(auto &in : pv.spec_precalc_after.info) map[in.i] = true;
		}
	}
	
	SpecPrecalc spec;
	spec.list_time.push_back(all_time);
	
	for(auto i = 0u; i < C; i++){
		if(map[i]){
			PrecalcInfo pi; pi.i = i;
			if(calcu[i].time_dep) pi.tlist = 0;
			else pi.tlist = UNSET;	
			spec.info.push_back(pi);			
		}
	}
	
	return spec;
}


/// Sets varaible all time
void Precalc::set_all_time()
{
	all_time = seq_vec(details.T);
}


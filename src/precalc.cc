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

Precalc::Precalc(const vector <SpeciesSimp> &species, const vector <Spline> &spline, const vector <Param> &param, const vector <ParamVecEle> &param_vec, const vector <Population> &pop, Constant &constant, const vector <double> &timepoint, const Details &details) : species(species), spline(spline), param(param), param_vec(param_vec), pop(pop), constant(constant), timepoint(timepoint), details(details)
{
	//clear_timer();
	num = 0;
	pcsize = 0;
}


/// Clears the timer (used for diagnostic purposes)
void Precalc::clear_timer()
{
	timer.clear(); timer.resize(20,0);
}


/// Adds an equation onto the precalculated equation but only for param and spline
void Precalc::add_eqn_simp(vector <Calculation> &calc, const vector <unsigned int> &param_vec_refq, const vector <unsigned int> &spline_refq)
{
	// Removes any reference to parameter and spline 
	for(auto &ca : calc){
		for(auto &it : ca.item){
			switch(it.type){
			case PARAMVEC: 
				it.type = REG_PRECALC; 
				it.num = param_vec_refq[it.num];
				break;
				
			case SPLINEREF: 
				it.type = REG_PRECALC_TIME; 
				it.num = spline_refq[it.num];
				if(it.num == UNSET) emsg("spline_ref unset");
				break;
			
			default:
				break;
			}
		}
	}
}


/// Adds an equation onto the precalculated equation calculation
void Precalc::add_eqn(vector <Calculation> &calc, const vector <unsigned int> &param_vec_refq, const vector <unsigned int> &spline_refq, SpecPrecalc &spec_precalc)
{
	// Removes any reference to parameter and spline 
	for(auto &ca : calc){
		for(auto &it : ca.item){		
			switch(it.type){
			case PARAMVEC: 
				it.type = REG_PRECALC; 
				it.num = param_vec_refq[it.num];
				if(it.num == UNSET) emsg("param_vec_ref unset");
				break;
				
			case SPLINEREF: 
				{
					it.type = REG_PRECALC_TIME; 
					it.num = spline_refq[it.num];
					if(it.num == UNSET) emsg("spline_ref unset");
				}
				break;
			
			default:
				break;
			}
		}
	}
	
	//auto stop_combine_fl = false;
	
	//if(add_type == PRECALC_PARAM_ONLY) return stop_combine_fl;
	
	vector <unsigned int> become_Rrecalc(calc.size(),UNSET);

	bool fl;
	do{
		fl = combine_multiply_add(calc,become_Rrecalc); // Tries to simplify multiplications and additions
		//fl = false;
		
		for(auto i = 0u; i < calc.size(); i++){
			if(become_Rrecalc[i] == UNSET){
				auto &ca = calc[i];
				
				auto &item = ca.item;
				
				auto nregderpop = 0u;
				for(auto j = 0u; j < item.size(); j++){
					switch(item[j].type){
					case POPNUM: case REG: case DERIVE: case INTEGRAL: case POPNUMTIME: nregderpop++; break;
					default: break;
					}
				}
			
				// This replaces calculation with that in non-pop 
				if(nregderpop == 0 && i+1 != calc.size()){  	
					auto it = add(ca,spec_precalc);
					become_Rrecalc[i] = it.num;
				}
				else{ // Attempts to put part of calculation into non-pop
					if(ca.op == DIVIDE){ // If dividing by something not 
						auto &denom = item[1];
						switch(denom.type){
						case REG: case POPNUM: case DERIVE: case INTEGRAL: case POPNUMTIME: break;
						default:
							{
								ca.op = MULTIPLY;
								
								Calculation calc;
								calc.op = DIVIDE;
									
								EqItem one; one.type = ONE;
								calc.item.push_back(one);	
								calc.item.push_back(denom);
								
								denom = add(calc,spec_precalc);
							}
							break;
						}
					}
					
					auto flag = false;
					
					if(ca.op == MULTIPLY && i+1 < calc.size() && calc[i+1].op == DIVIDE){
						if(item.size() - nregderpop >= 1){
							auto &ca_next = calc[i+1];
							const auto &itdiv = ca_next.item;
							if(itdiv[0].type == REG && itdiv[0].num == i){
								auto &denom = itdiv[1];
								
								switch(denom.type){
								case REG: case POPNUM: case DERIVE: case INTEGRAL: case POPNUMTIME: break;
								default:
									{
										flag = true;
									
										Calculation calc;
										calc.op = DIVIDE;
									
										EqItem one; one.type = ONE;
										calc.item.push_back(one);
										calc.item.push_back(denom);
										
										auto inew = add(calc,spec_precalc);
										
										Calculation calc2;
										calc2.op = MULTIPLY;
								
										vector <EqItem> item_new;
								
										for(auto j = 0u; j < item.size(); j++){
											const auto &it = item[j];
											switch(it.type){
											case REG: case POPNUM: case DERIVE: case INTEGRAL: case POPNUMTIME: item_new.push_back(it); break;
											default: calc2.item.push_back(it); break;
											}
										}
										calc2.item.push_back(inew);
									
										auto inew2 = add(calc2,spec_precalc);
										item_new.push_back(inew2);
										
										ca_next.op = MULTIPLY;
										ca_next.item = item_new;
										
										become_Rrecalc[i] = USINT_MAX;
									}
									break;
								}
							}
						}
					}
					
					if(flag == false && (ca.op == MULTIPLY || ca.op == ADD)){
						if(item.size() - nregderpop >= 2){
							Calculation calc;
							calc.op = ca.op;
							
							vector <EqItem> item_new;
							
							for(auto j = 0u; j < item.size(); j++){
								const auto &it = item[j];
								switch(it.type){
								case REG: case POPNUM: case DERIVE: case INTEGRAL: case POPNUMTIME: item_new.push_back(it); break;
								default: calc.item.push_back(it); break;
								}
							}
					
							auto inew = add(calc,spec_precalc);
							item_new.push_back(inew);
							
							ca.item = item_new;
						}
					}
					
					if(flag == false && nregderpop == 0){
						switch(ca.op){
						case EXPFUNC: case SINFUNC: case COSFUNC: case LOGFUNC: 	// Single parameter functions
						case STEPFUNC: case ABSFUNC: case SQRTFUNC: case SIGFUNC:
						case POWERFUNC: case THRESHFUNC: case UBOUNDFUNC:         // Two parameter functions
						case MAXFUNC: case MINFUNC:
							{
								//Calculation calc;
								//calc.op = ca.op;
								//pcalc.item = ca.item;
							
								//auto inew = add(pcalc,spec_precalc);
								auto inew = add(ca,spec_precalc);
								
								vector <EqItem> item_new;
								item_new.push_back(inew);
								ca.op = ADD;
								ca.item = item_new;
							}
							break;
					
						default: break;
						}
					}
				}
			}
		}
	
		// Replaces register with non_popcalc	
		for(auto i = 0u; i < calc.size(); i++){
			if(become_Rrecalc[i] == UNSET){
				for(auto &it : calc[i].item){
					if(it.type == REG && become_Rrecalc[it.num] != UNSET){
						it.type = REG_PRECALC; 
						it.num = become_Rrecalc[it.num];
						if(it.num == CODE) emsg("Should not be code");
						if(pcalcu[it.num].time_dep) it.type = REG_PRECALC_TIME; 
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
	
	//print_calc_section(iii);
}


/// Combines together multiplication and addition operations to simplify
bool Precalc::combine_multiply_add(vector <Calculation> &calcu, vector <unsigned int> &become_Rrecalc) const 
{
	auto imax = calcu.size();
	vector <unsigned int> used(imax,0);

	auto fl = false;
	
	for(auto i = 0u; i < imax; i++){       // Determines how many times each resigter is used
		if(become_Rrecalc[i] == UNSET){	
			const auto &ca = calcu[i];
			for(const auto &it : ca.item){
				if(it.type == REG) used[it.num]++;
			}
		}
	}
	
	for(auto i = 0u; i < imax; i++){
		if(become_Rrecalc[i] == UNSET){	
			auto &ca = calcu[i];
			if(ca.op == MULTIPLY){
				auto cont = false;
				for(auto j = 0u; j < ca.item.size(); j++){
					auto &it = ca.item[j];
					if(it.type == REG && used[it.num] == 1 && calcu[it.num].op == MULTIPLY) cont = true;
				}
				
				if(cont){
					vector <EqItem> item;
					for(auto j = 0u; j < ca.item.size(); j++){
						auto &it = ca.item[j];
						auto num = it.num;
						if(it.type == REG && used[num] == 1 && calcu[num].op == MULTIPLY){
							for(const auto &it_fr : calcu[num].item) item.push_back(it_fr);
							become_Rrecalc[num] = CODE;
						}
						else{
							item.push_back(it);
						}
					}
					ca.item = item;
					fl = true;
				}
			}
			
			if(ca.op == ADD){
				auto cont = false;
				for(auto j = 0u; j < ca.item.size(); j++){
					auto &it = ca.item[j];
					if(it.type == REG && used[it.num] == 1 && calcu[it.num].op == ADD) cont = true;
				}
				
				if(cont){
					vector <EqItem> item;
					for(auto j = 0u; j < ca.item.size(); j++){
						auto &it = ca.item[j];
						auto num = it.num;
						if(it.type == REG && used[num] == 1 && calcu[num].op == ADD){
							for(const auto &it_fr : calcu[num].item) item.push_back(it_fr);
							become_Rrecalc[num] = CODE;
						}
						else{
							item.push_back(it);
						}
					}
					ca.item = item;
					fl = true;
				}
			}
		}
	}
	
	return fl;
}


/// Adds a parameter to precalculation
unsigned int Precalc::add_param(unsigned int th)
{
	PreCalc ca;
	ca.op = SINGLE;
	ca.time_dep = false;
	
	PreEqItem it;
	it.type = PRE_PARAMVEC;
	it.num = th; if(th == UNSET) emsg("Parem vec problem");
	ca.pre_item.push_back(it);

	auto vec = get_vec(ca);
	
	auto q = hash_ca.existing(vec);
	if(q == UNSET){
		q = pcalcu.size();
		calcu_add(ca);
		hash_ca.add(q,vec);
	}

	return q;
}


/// Adds a calculation to the list
void Precalc::calcu_add(PreCalc &ca)
{
	if(ca.op == NOOP){ 
		emsg("SHould not add noop");
	}
	else{
		ca.iref = pcsize;
		pcalcu.push_back(ca);
		if(pcalcu.size() >= USINT_MAX) emsg("BICI ran out of memory");
		if(ca.time_dep) pcsize += details.T;
		else pcsize++;
	}
}


/// Adds a spline to precalculation
unsigned int Precalc::add_spline(unsigned int s, SpecPrecalc &spec_precalc)
{
	PreCalc ca;
	ca.op = SINGLE;
	ca.time_dep = true;
	
	PreEqItem it;
	it.type = PRE_SPLINEREF;
	it.num = s;
	ca.pre_item.push_back(it);

	auto vec = get_vec(ca);
	
	auto q = hash_ca.existing(vec);
	if(q == UNSET){
		q = pcalcu.size();
		calcu_add(ca);
		
		sp_add(spec_precalc,q,all_time);
		
		hash_ca.add(q,vec);
	}
	
	return q;
}


/// Adds a calculation 
EqItem Precalc::add(Calculation &ca, SpecPrecalc &spec_precalc)
{
	PreCalc pcalc; pcalc.op = ca.op; 
	pcalc.time_dep = false;
	
	for(const auto &it : ca.item){
		PreEqItem pi; pcalc.op = ca.op; 
		switch(it.type){
		case PARAMVEC: pi.type = PRE_PARAMVEC; break;
		case SPLINEREF: pi.type = PRE_SPLINEREF; break;
		case CONSTSPLINEREF: pi.type = PRE_CONSTSPLINEREF; break;
		case REG_PRECALC: pi.type = PRE_REG; break;
		case REG_PRECALC_TIME: pi.type = PRE_REGTIME; break;
		case ONE: pi.type = PRE_ONE; break;
		case ZERO: pi.type = PRE_ZERO; break;
		case NUMERIC: pi.type = PRE_NUMERIC; break;
		case TIME: pi.type = PRE_TIME; break;
		default: emsg("Type not recognised:"+tstr(it.type)); break;
		}
		pi.num = it.num;
		pcalc.pre_item.push_back(pi);
		
		// Determines if time dependent
		switch(pi.type){
		case PRE_SPLINEREF: case PRE_CONSTSPLINEREF: case PRE_TIME: case PRE_REGTIME:
			pcalc.time_dep = true; 
			break;
		default: break;
		}
	}

	auto vec = get_vec(pcalc);
			
	auto q = hash_ca.existing(vec);
	if(q == UNSET){
		q = pcalcu.size();
		
		hash_ca.add(q,vec);
		calcu_add(pcalc);
	
		num++;
	
		sp_add(spec_precalc,q,all_time);
	}	
	
	EqItem inew;
	inew.type = REG_PRECALC; if(pcalc.time_dep == true) inew.type = REG_PRECALC_TIME;
	inew.num = q;
	
	return inew;
}


/// Gets a hash vector for a calulation item
vector <unsigned int> Precalc::get_vec(const PreCalc &ca) const 
{
	const auto &pre_item = ca.pre_item;
	
	vector <unsigned int> vec;
		
	vec.push_back(ca.op);
	for(auto j = 0u; j < pre_item.size(); j++){
		vec.push_back(pre_item[j].type);
		vec.push_back(pre_item[j].num);
	}

	return vec;
}	


/// Prints steps used for a calculation
void Precalc::print_calc() const   
{
	stringstream ss;
	ss << "PRE-CALCULATION:" << endl;
	
	auto li = 0u, limax = 10000000u;
  for(auto q = 0u; q < pcalcu.size(); q++){
		const auto &ca = pcalcu[q];
		ss << q << " " << "(i" << ca.iref << ")" << " ....";
		ss << print_ca(q,ca);
		ss << endl;
	
		li++;
		if(li > limax) break;
  }
	
	if(li > limax) ss << "...";

  ss << endl << endl;
	
	cout << "TURN OFF PRECALC" << endl;
	ofstream fout("precalc.txt");
	fout << ss.str();
}


/// Prints a SpecPrecalc to a file
void Precalc::print_spec_precalc(string name, const SpecPrecalc &spec) const 
{
	auto te = print_spec_precalc_str(name,spec);
	
	if(true){
		cout << te;
	}
	else{
		ofstream fout(name);
		fout << te;
	}
}


string Precalc::print_spec_precalc_str(string name, const SpecPrecalc &spec) const 
{
	stringstream ss;
	cout << "Output " << name << ": " << endl;
	for(const auto &va : spec.info){
		ss << va.q << " " << va.tlist << endl;
	}
	ss << "Time lines" << endl;
	for(const auto &tl : spec.list_time){
		for(auto ti : tl) ss << ti << ",";
		ss << endl;		
	}
	
	return ss.str();
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


/// Calculcates initial value for precalc (using any constant values_)
/// Calculates the value for an equation
vector <double> Precalc::calculate_precalc_init(const SpecPrecalc &spec_precalc) const 
{
	vector <double> precalc(pcsize,UNSET);
	const auto &cval = constant.value;
	
	auto calc_err = NO_ERROR;
	
	for(const auto &in : spec_precalc.info){
		auto q = in.q;
		
		const auto &ca = pcalcu[q];
		auto i = ca.iref;
		
		if(ca.op == SINGLE){
			const auto &it = ca.pre_item[0];
			switch(it.type){
			case PRE_PARAMVEC: 
				break;
			
			case PRE_SPLINEREF:
				{
					const auto &spl = spline[it.num];
					if(spl.constant){
						for(auto ti : all_time){
							if(spl.dynamic) precalc[i+ti] = UNSET;
							else precalc[i+ti] = spl.const_val[ti];
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
			const auto &pre_item = ca.pre_item;
			const auto N = pre_item.size();
			
			vector <double> num(N);
			
			if(ca.time_dep){
				for(auto ti : all_time){
					auto fl = false;
					for(auto j = 0u; j < N; j++){
						const auto &it = pre_item[j];
						
						switch(it.type){
							case PRE_REG: num[j] = precalc[pcalcu[it.num].iref]; break;
							case PRE_REGTIME: num[j] = precalc[pcalcu[it.num].iref+ti]; break;							
							case PRE_ONE: num[j] = 1; break;
							case PRE_NUMERIC: num[j] = constant.value[it.num]; break;
							case PRE_CONSTSPLINEREF: 
								{
									const auto &spl = spline[it.num];
									if(spl.dynamic) num[j] = UNSET;
									else num[j] = spl.const_val[ti]; 
								}
								break;
							case PRE_TIME: num[j] = timepoint[ti]; break;
							default: eqn_type_error(it.type,6); break;
						}
						if(num[j] == UNSET){ fl = true; break;}
					}

					if(fl == false) precalc[i+ti] = calculate_operation(ca.op,num,calc_err);
				}
			}
			else{			
				auto fl = false;
				for(auto j = 0u; j < N; j++){
					const auto &it = pre_item[j];	
					switch(it.type){	
						case PRE_ONE: num[j] = 1; break;
						case PRE_REG: num[j] = precalc[pcalcu[it.num].iref]; break;
						case PRE_NUMERIC: num[j] = cval[it.num]; break;
						default: eqn_type_error(it.type,8); break;
					}
					if(num[j] == UNSET){ fl = true; break;}
				}

				if(fl == false) precalc[i] = calculate_operation(ca.op,num,calc_err);
			}
		}
  }

	return precalc;
}


/// Gets a list_time from a spec precalculation
const vector <unsigned int>&	Precalc::get_list_time(unsigned int cit, const SpecPrecalc &spec) const
{
	if(cit == ALL_TIME_STEP) return all_time;
	return spec.list_time[cit];
}					
						

/// Calculates the value for an equation
void Precalc::calculate(const SpecPrecalc &spec_calc, PV &param_val, bool store) const 
{
	const auto &info = spec_calc.info;
	if(info.size() == 0) return;
	
	const auto &value = param_val.value;
	const auto &cval = constant.value;
	auto &precalc = param_val.precalc;
	
	auto calc_err = NO_ERROR;
	
	for(const auto &ci : info){
		auto q = ci.q;
	
		const auto &ca = pcalcu[q];
		auto i = ca.iref;
		
		if(store){
			if(ca.time_dep){ 
				const auto &list_time = get_list_time(ci.tlist,spec_calc);
				for(auto ti : list_time) param_val.precalc_change(i+ti);
			}
			else{
				param_val.precalc_change(i);
			}
		}
		
		if(ca.op == SINGLE){
			const auto &it = ca.pre_item[0];
		
			switch(it.type){
			case PRE_PARAMVEC: 
				precalc[i] = value[it.num];
				break;
			
			case PRE_SPLINEREF:
				{
					const auto &spl = spline[it.num];
					if(!spl.constant){
						const auto &list_time = get_list_time(ci.tlist,spec_calc);
				
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
			if(ca.time_dep){
				//calc_time_varying(i,ca,get_list_time(ci.tlist,spec_calc),precalc,calc_err);
				calc_time_varying_fast(i,ca,get_list_time(ci.tlist,spec_calc),precalc,calc_err);
			}
			else{			
				const auto &pre_item = ca.pre_item;
				const auto N = pre_item.size();
				
				vector <double> num(N);
			
				for(auto j = 0u; j < N; j++){
					const auto &it = pre_item[j];	
					switch(it.type){	
						case PRE_ONE: num[j] = 1; break;
						case PRE_REG: num[j] = precalc[pcalcu[it.num].iref]; break;
						case PRE_NUMERIC: num[j] = cval[it.num]; break;
						default: eqn_type_error(it.type,8); break;
					}
					
					if(num[j] == UNSET) emsg("Precalc unset2");
				}

					
				precalc[i] = calculate_operation(ca.op,num,calc_err);
			}
		}
  }
	
	if(calc_err != NO_ERROR) calc_error(calc_err,"Within precalculation");
}


/// Calcuates time varying precalculation
void Precalc::calc_time_varying(unsigned int i, const PreCalc &ca, const vector <unsigned int> &list_time, vector <double> &precalc, CalcError &calc_err) const
{				
	const auto &pre_item = ca.pre_item;
	const auto N = pre_item.size();
	
	vector <double> num(N);
	for(auto ti : list_time){
		for(auto j = 0u; j < N; j++){
			const auto &it = pre_item[j];
			
			switch(it.type){
				case PRE_REG: 
					num[j] = precalc[pcalcu[it.num].iref]; 
					if(num[j] == UNSET){ // Just needed for checking
						cout << it.num << endl; 
						emsg("Reg Precalc unset");
					}
					break;
				case PRE_REGTIME:
					num[j] = precalc[pcalcu[it.num].iref+ti]; 
					if(num[j] == UNSET){  // Just needed for checking
						cout << it.num << " " << ti << endl; 
						emsg("Regt Precalc unset");
					}
					break;							
				case PRE_ONE: num[j] = 1; break;
				case PRE_ZERO: num[j] = 0; break;
				case PRE_NUMERIC: num[j] = constant.value[it.num]; break;
				case PRE_CONSTSPLINEREF: num[j] = spline[it.num].const_val[ti]; break;
				case PRE_TIME: num[j] = timepoint[ti]; break;
				default: eqn_type_error(it.type,7); break;
			}
		}

		precalc[i+ti] = calculate_operation(ca.op,num,calc_err);
	}
}


/// Get a non-timevarying quantity
double Precalc::get_val(const PreEqItem &it, const vector <double> &precalc) const 
{
	switch(it.type){
	case PRE_REG: return precalc[pcalcu[it.num].iref]; 
	case PRE_NUMERIC: return constant.value[it.num];
	case PRE_ONE: return 1;
	case PRE_ZERO: return 0;
	default: return UNSET;
	}
}


/// Calcuates time varying precalculation
void Precalc::calc_time_varying_fast(unsigned int i, const PreCalc &ca, const vector <unsigned int> &list_time, vector <double> &precalc, CalcError &calc_err) const
{				
	const auto &item = ca.pre_item;
	const auto N = item.size();
	auto *prec = &precalc[i];
	
	switch(ca.op){
	case ADD:
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
		
			auto val = get_val(it,precalc);
			if(val != UNSET){
				if(j == 0){ for(auto ti : list_time) prec[ti] = val;}
				else{ for(auto ti : list_time) prec[ti] += val;}
			}
			else{
				auto num = it.num;
				switch(it.type){
				case PRE_REGTIME:
					{
						auto *prec2 = &precalc[pcalcu[num].iref];
						if(j == 0){ for(auto ti : list_time) prec[ti] = prec2[ti];}
						else{ for(auto ti : list_time) prec[ti] += prec2[ti];}
					}
					break;
					
				case PRE_CONSTSPLINEREF:
					{
						const auto &const_val = spline[num].const_val;
						if(j == 0){ for(auto ti : list_time) prec[ti] = const_val[ti];}
						else{ for(auto ti : list_time) prec[ti] += const_val[ti];}
					}
					break;
					
				case PRE_TIME: 
					{
						if(j == 0){ for(auto ti : list_time) prec[ti] = timepoint[ti];}
						else{ for(auto ti : list_time) prec[ti] += timepoint[ti];}
					}
					break;
					
				default: emsg("calc_time_varying_fast problem"); break;
				}
			}
		}
		break;
		
	case MULTIPLY:
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
		
			auto val = get_val(it,precalc);
			if(val != UNSET){
				if(j == 0){ for(auto ti : list_time) prec[ti] = val;}
				else{ for(auto ti : list_time) prec[ti] *= val;}
			}
			else{
				auto num = it.num;
				switch(it.type){
				case PRE_REGTIME:
					{
						auto *prec2 = &precalc[pcalcu[num].iref];
						if(j == 0){ for(auto ti : list_time) prec[ti] = prec2[ti];}
						else{ for(auto ti : list_time) prec[ti] *= prec2[ti];}
					}
					break;
					
				case PRE_CONSTSPLINEREF:
					{
						const auto &const_val = spline[num].const_val;
						if(j == 0){ for(auto ti : list_time) prec[ti] = const_val[ti];}
						else{ for(auto ti : list_time) prec[ti] *= const_val[ti];}
					}
					break;
					
				case PRE_TIME: 
					{
						if(j == 0){ for(auto ti : list_time) prec[ti] = timepoint[ti];}
						else{ for(auto ti : list_time) prec[ti] *= timepoint[ti];}
					}
					break;
					
				default: emsg("calc_time_varying_fast problem"); break;
				}
			}
		}
		break;
	
	case TAKE:
		emsg("Should not be take"); 
		break;
	
	case DIVIDE:
		for(auto j = 0u; j < 2; j++){
			const auto &it = item[j];
		
			auto val = get_val(it,precalc);
			if(val != UNSET){
				if(j == 0){ for(auto ti : list_time) prec[ti] = val;}
				else{ 
					if(val <= 0){ calc_err = DIV_BY_ZERO_ERROR; val = TINY;}
					for(auto ti : list_time) prec[ti] /= val;
				}
			}
			else{
				auto num = it.num;
				switch(it.type){
				case PRE_REGTIME:
					{
						auto *prec2 = &precalc[pcalcu[num].iref];
						if(j == 0){ for(auto ti : list_time) prec[ti] = prec2[ti];}
						else{ 
							for(auto ti : list_time){
								auto val = prec2[ti];
								if(val <= 0){ calc_err = DIV_BY_ZERO_ERROR; val = TINY;}
								prec[ti] /= val;
							}
						}
					}
					break;
					
				case PRE_CONSTSPLINEREF:
					{
						const auto &const_val = spline[num].const_val;
						if(j == 0){ for(auto ti : list_time) prec[ti] = const_val[ti];}
						else{ 
							for(auto ti : list_time){
								auto val = const_val[ti];
								if(val <= 0){ calc_err = DIV_BY_ZERO_ERROR; val = TINY;}
								prec[ti] /= val;
							}
						}
					}
					break;
					
				case PRE_TIME: 
					{
						if(j == 0){ for(auto ti : list_time) prec[ti] = timepoint[ti];}
						else{ 
							for(auto ti : list_time){
								auto val = timepoint[ti];
								if(val <= 0){ calc_err = DIV_BY_ZERO_ERROR; val = TINY;}
								prec[ti] /= val;
							}
						}
					}
					break;
					
				default: emsg("calc_time_varying_fast problem"); break;
				}
			}
		}
		break;
	
	case EXPFUNC:
		{
			const auto &it = item[0];
		
			auto num = it.num;
			switch(it.type){
			case PRE_REGTIME:
				{
					auto *prec2 = &precalc[pcalcu[num].iref];
					for(auto ti : list_time) prec[ti] = exp(prec2[ti]);
				}
				break;
				
			case PRE_CONSTSPLINEREF:
				{
					const auto &const_val = spline[num].const_val;
					for(auto ti : list_time) prec[ti] = exp(const_val[ti]);
				}
				break;
				
			case PRE_TIME: 
				{
					for(auto ti : list_time) prec[ti] = exp(timepoint[ti]);
				}
				break;
				
			default: emsg("calc_time_varying_fast problem"); break;
			}
		}
		break;
		
	case SINFUNC:
		{
			const auto &it = item[0];
		
			auto num = it.num;
			switch(it.type){
			case PRE_REGTIME:
				{
					auto *prec2 = &precalc[pcalcu[num].iref];
					for(auto ti : list_time) prec[ti] = sin(prec2[ti]);
				}
				break;
				
			case PRE_CONSTSPLINEREF:
				{
					const auto &const_val = spline[num].const_val;
					for(auto ti : list_time) prec[ti] = sin(const_val[ti]);
				}
				break;
				
			case PRE_TIME: 
				{
					for(auto ti : list_time) prec[ti] = sin(timepoint[ti]);
				}
				break;
				
			default: emsg("calc_time_varying_fast problem"); break;
			}
		}
		break;
		
	case COSFUNC:
		{
			const auto &it = item[0];
		
			auto num = it.num;
			switch(it.type){
			case PRE_REGTIME:
				{
					auto *prec2 = &precalc[pcalcu[num].iref];
					for(auto ti : list_time) prec[ti] = cos(prec2[ti]);
				}
				break;
				
			case PRE_CONSTSPLINEREF:
				{
					const auto &const_val = spline[num].const_val;
					for(auto ti : list_time) prec[ti] = cos(const_val[ti]);
				}
				break;
				
			case PRE_TIME: 
				{
					for(auto ti : list_time) prec[ti] = cos(timepoint[ti]);
				}
				break;
				
			default: emsg("calc_time_varying_fast problem"); break;
			}
		}
		break;
		
	case LOGFUNC:
		{
			const auto &it = item[0];
		
			auto num = it.num;
			switch(it.type){
			case PRE_REGTIME:
				{
					auto *prec2 = &precalc[pcalcu[num].iref];
					for(auto ti : list_time){
						auto val = prec2[ti]; if(val <= 0){ calc_err = LOG_NEG_ERROR; val = TINY;}
						prec[ti] = log(val);
					}
				}
				break;
				
			case PRE_CONSTSPLINEREF:
				{
					const auto &const_val = spline[num].const_val;
					for(auto ti : list_time){
						auto val = const_val[ti]; if(val <= 0){ calc_err = LOG_NEG_ERROR; val = TINY;}
						prec[ti] = log(val);
					}
				}
				break;
				
			case PRE_TIME: 
				{
					for(auto ti : list_time){
						auto val = timepoint[ti]; if(val <= 0){ calc_err = LOG_NEG_ERROR; val = TINY;}
						prec[ti] = log(val);
					}
				}
				break;
				
			default: emsg("calc_time_varying_fast problem"); break;
			}
		}
		break;
	
	case STEPFUNC:
		{
			const auto &it = item[0];
		
			auto num = it.num;
			switch(it.type){
			case PRE_REGTIME:
				{
					auto *prec2 = &precalc[pcalcu[num].iref];
					for(auto ti : list_time){
						if(prec2[ti] > 0) prec[ti] = 1; 
						else prec[ti] = 0;
					}
				}
				break;
				
			case PRE_CONSTSPLINEREF:
				{
					const auto &const_val = spline[num].const_val;
					for(auto ti : list_time){
						if(const_val[ti] > 0) prec[ti] = 1; 
						else prec[ti] = 0;
					}
				}
				break;
				
			case PRE_TIME: 
				{
					for(auto ti : list_time){
						if(timepoint[ti] > 0) prec[ti] = 1; 
						else prec[ti] = 0;
					}
				}
				break;
				
			default: emsg("calc_time_varying_fast problem"); break;
			}
		}
		break;
	
	case POWERFUNC:
		for(auto j = 0u; j < 2; j++){
			const auto &it = item[j];
		
			auto val = get_val(it,precalc);
			if(val != UNSET){
				if(j == 0){ for(auto ti : list_time) prec[ti] = val;}
				else{ for(auto ti : list_time) prec[ti] = pow(prec[ti],val);}
			}
			else{
				auto num = it.num;
				switch(it.type){
				case PRE_REGTIME:
					{
						auto *prec2 = &precalc[pcalcu[num].iref];
						if(j == 0){ for(auto ti : list_time) prec[ti] = prec2[ti];}
						else{ for(auto ti : list_time) prec[ti] = pow(prec[ti],prec2[ti]);}
					}
					break;
					
				case PRE_CONSTSPLINEREF:
					{
						const auto &const_val = spline[num].const_val;
						if(j == 0){ for(auto ti : list_time) prec[ti] = const_val[ti];}
						else{ for(auto ti : list_time) prec[ti] = pow(prec[ti],const_val[ti]);}
					}
					break;
					
				case PRE_TIME: 
					{
						if(j == 0){ for(auto ti : list_time) prec[ti] = timepoint[ti];}
						else{ for(auto ti : list_time) prec[ti] = pow(prec[ti],timepoint[ti]);}
					}
					break;
					
				default: emsg("calc_time_varying_fast problem"); break;
				}
			}
		}
		break;
		
	case THRESHFUNC:
		for(auto j = 0u; j < 2; j++){
			const auto &it = item[j];
		
			auto val = get_val(it,precalc);
			if(val != UNSET){
				if(j == 0){ for(auto ti : list_time) prec[ti] = val;}
				else{ 
					for(auto ti : list_time){ 
						if(prec[ti] < val) prec[ti] = 0;
					}
				}
			}
			else{
				auto num = it.num;
				switch(it.type){
				case PRE_REGTIME:
					{
						auto *prec2 = &precalc[pcalcu[num].iref];
						if(j == 0){ for(auto ti : list_time) prec[ti] = prec2[ti];}
						else{ 
							for(auto ti : list_time){
								if(prec[ti] < prec2[ti]) prec[ti] = 0;
							}
						}
					}
					break;
					
				case PRE_CONSTSPLINEREF:
					{
						const auto &const_val = spline[num].const_val;
						if(j == 0){ for(auto ti : list_time) prec[ti] = const_val[ti];}
						else{ 
							for(auto ti : list_time){
								if(prec[ti] < const_val[ti]) prec[ti] = 0; 
							}
						}
					}
					break;
					
				case PRE_TIME: 
					{
						if(j == 0){ for(auto ti : list_time) prec[ti] = timepoint[ti];}
						else{ 
							for(auto ti : list_time){
								if(prec[ti] < timepoint[ti]) prec[ti] = 0; 
							}
						}
					}
					break;
					
				default: emsg("calc_time_varying_fast problem"); break;
				}
			}
		}
		break;
		
	case UBOUNDFUNC:
		for(auto j = 0u; j < 2; j++){
			const auto &it = item[j];
		
			auto val = get_val(it,precalc);
			if(val != UNSET){
				if(j == 0){ for(auto ti : list_time) prec[ti] = val;}
				else{ 
					for(auto ti : list_time){ 
						if(prec[ti] > val) prec[ti] = INFINITY;
					}
				}
			}
			else{
				auto num = it.num;
				switch(it.type){
				case PRE_REGTIME:
					{
						auto *prec2 = &precalc[pcalcu[num].iref];
						if(j == 0){ for(auto ti : list_time) prec[ti] = prec2[ti];}
						else{ 
							for(auto ti : list_time){
								if(prec[ti] > prec2[ti]) prec[ti] = INFINITY;
							}
						}
					}
					break;
					
				case PRE_CONSTSPLINEREF:
					{
						const auto &const_val = spline[num].const_val;
						if(j == 0){ for(auto ti : list_time) prec[ti] = const_val[ti];}
						else{ 
							for(auto ti : list_time){
								if(prec[ti] > const_val[ti]) prec[ti] = INFINITY; 
							}
						}
					}
					break;
					
				case PRE_TIME: 
					{
						if(j == 0){ for(auto ti : list_time) prec[ti] = timepoint[ti];}
						else{ 
							for(auto ti : list_time){
								if(prec[ti] > timepoint[ti]) prec[ti] = INFINITY; 
							}
						}
					}
					break;
					
				default: emsg("calc_time_varying_fast problem"); break;
				}
			}
		}
		break;
		
	case MAXFUNC:
		for(auto j = 0u; j < 2; j++){
			const auto &it = item[j];
		
			auto val = get_val(it,precalc);
			if(val != UNSET){
				if(j == 0){ for(auto ti : list_time) prec[ti] = val;}
				else{ 
					for(auto ti : list_time){ 
						if(val > prec[ti]) prec[ti] = val;
					}
				}
			}
			else{
				auto num = it.num;
				switch(it.type){
				case PRE_REGTIME:
					{
						auto *prec2 = &precalc[pcalcu[num].iref];
						if(j == 0){ for(auto ti : list_time) prec[ti] = prec2[ti];}
						else{ 
							for(auto ti : list_time){
								if(prec2[ti] > prec[ti]) prec[ti] = prec2[ti];
							}
						}
					}
					break;
					
				case PRE_CONSTSPLINEREF:
					{
						const auto &const_val = spline[num].const_val;
						if(j == 0){ for(auto ti : list_time) prec[ti] = const_val[ti];}
						else{ 
							for(auto ti : list_time){
								if(const_val[ti] > prec[ti]) prec[ti] = const_val[ti]; 
							}
						}
					}
					break;
					
				case PRE_TIME: 
					{
						if(j == 0){ for(auto ti : list_time) prec[ti] = timepoint[ti];}
						else{ 
							for(auto ti : list_time){
								if(timepoint[ti] > prec[ti]) prec[ti] = timepoint[ti]; 
							}
						}
					}
					break;
					
				default: emsg("calc_time_varying_fast problem"); break;
				}
			}
		}
		break;
		
	case MINFUNC:
		for(auto j = 0u; j < 2; j++){
			const auto &it = item[j];
		
			auto val = get_val(it,precalc);
			if(val != UNSET){
				if(j == 0){ for(auto ti : list_time) prec[ti] = val;}
				else{ 
					for(auto ti : list_time){ 
						if(val < prec[ti]) prec[ti] = val;
					}
				}
			}
			else{
				auto num = it.num;
				switch(it.type){
				case PRE_REGTIME:
					{
						auto *prec2 = &precalc[pcalcu[num].iref];
						if(j == 0){ for(auto ti : list_time) prec[ti] = prec2[ti];}
						else{ 
							for(auto ti : list_time){
								if(prec2[ti] < prec[ti]) prec[ti] = prec2[ti];
							}
						}
					}
					break;
					
				case PRE_CONSTSPLINEREF:
					{
						const auto &const_val = spline[num].const_val;
						if(j == 0){ for(auto ti : list_time) prec[ti] = const_val[ti];}
						else{ 
							for(auto ti : list_time){
								if(const_val[ti] < prec[ti]) prec[ti] = const_val[ti]; 
							}
						}
					}
					break;
					
				case PRE_TIME: 
					{
						if(j == 0){ for(auto ti : list_time) prec[ti] = timepoint[ti];}
						else{ 
							for(auto ti : list_time){
								if(timepoint[ti] < prec[ti]) prec[ti] = timepoint[ti]; 
							}
						}
					}
					break;
					
				default: emsg("calc_time_varying_fast problem"); break;
				}
			}
		}
		break;
		
	case ABSFUNC:
		{
			const auto &it = item[0];
		
			auto num = it.num;
			switch(it.type){
			case PRE_REGTIME:
				{
					auto *prec2 = &precalc[pcalcu[num].iref];
					for(auto ti : list_time){
						prec[ti] = prec2[ti];
						if(prec[ti] < 0) prec[ti] = -prec[ti];
					}
				}
				break;
				
			case PRE_CONSTSPLINEREF:
				{
					const auto &const_val = spline[num].const_val;
					for(auto ti : list_time){
						prec[ti] = const_val[ti];
						if(prec[ti] < 0) prec[ti] = -prec[ti];
					}
				}
				break;
				
			case PRE_TIME: 
				{
					for(auto ti : list_time){
						prec[ti] = timepoint[ti];
						if(prec[ti] < 0) prec[ti] = -prec[ti];
					}
				}
				break;
				
			default: emsg("calc_time_varying_fast problem"); break;
			}
		}
		break;
		
	case SQRTFUNC:
		{
			const auto &it = item[0];
		
			auto num = it.num;
			switch(it.type){
			case PRE_REGTIME:
				{
					auto *prec2 = &precalc[pcalcu[num].iref];
					for(auto ti : list_time){
						auto val = prec2[ti]; if(val < 0){ calc_err = SQRT_NEG_ERROR; val = TINY;}
						prec[ti] = sqrt(val);
					}
				}
				break;
				
			case PRE_CONSTSPLINEREF:
				{
					const auto &const_val = spline[num].const_val;
					for(auto ti : list_time){
						auto val = const_val[ti]; if(val < 0){ calc_err = SQRT_NEG_ERROR; val = TINY;}
						prec[ti] = sqrt(val);
					}
				}
				break;
				
			case PRE_TIME: 
				{
					for(auto ti : list_time){
						auto val = timepoint[ti]; if(val < 0){ calc_err = SQRT_NEG_ERROR; val = TINY;}
						prec[ti] = sqrt(val);
					}
				}
				break;
				
			default: emsg("calc_time_varying_fast problem"); break;
			}
		}
		break;
		
	case SIGFUNC:
		{
			const auto &it = item[0];
		
			auto num = it.num;
			switch(it.type){
			case PRE_REGTIME:
				{
					auto *prec2 = &precalc[pcalcu[num].iref];
					for(auto ti : list_time) prec[ti] = 1/(1+exp(-prec2[ti]));
				}
				break;
				
			case PRE_CONSTSPLINEREF:
				{
					const auto &const_val = spline[num].const_val;
					for(auto ti : list_time) prec[ti] = 1/(1+exp(-const_val[ti]));
				}
				break;
				
			case PRE_TIME: 
				{
					for(auto ti : list_time){
						prec[ti] = 1/(1+exp(-timepoint[ti]));
					}
				}
				break;
				
			default: emsg("calc_time_varying_fast problem"); break;
			}
		}
		break;
		
	default:
		emsg("fast problem");
		break;
	}
	
	if(false){ // Testing
		cout << " test" << endl;
		auto T = details.T;
		vector <double> store;
		for(auto ti = 0u; ti < T; ti++) store.push_back(precalc[i+ti]);
		calc_time_varying(i,ca,list_time,precalc,calc_err);
		for(auto ti = 0u; ti < T; ti++){
			if(dif(precalc[i+ti],store[ti],TINY)) emsg("problem with precalcfast");
		}
	}
}


/// Prints the calculation from a specific calculation
string Precalc::print_ca(unsigned int q, const PreCalc &ca) const
{ 
	stringstream ss;
	
	ss << q << ": ";
	
	switch(ca.op){
	case EXPFUNC: ss <<  "exp("; break;
	case SINFUNC: ss <<  "sin("; break;
	case COSFUNC: ss <<  "cos("; break;
	case LOGFUNC: ss <<  "log("; break;
	case STEPFUNC: ss <<  "step("; break;
	case POWERFUNC: ss <<  "power("; break;
	case THRESHFUNC: ss <<  "thresh("; break;
	case UBOUNDFUNC: ss <<  "ubound("; break;
	case MAXFUNC: ss << "max("; break; 
	case MINFUNC: ss << "min("; break; 
	case ABSFUNC: ss << "abs("; break;
	case SQRTFUNC: ss << "sqrt("; break;
	case SIGFUNC: ss << "sigmoid("; break;
	case ADD: break;
	case TAKE: break;
	case MULTIPLY: break;
	case SINGLE: break;
	case DIVIDE: break;
	default: emsg_input("Eq problem1"); break;
	}
		
	const auto &item = ca.pre_item;
	
	if(item.size() > 0 && (item[0].type == PRE_SPLINEREF || item[0].type == PRE_CONSTSPLINEREF)){
		const auto &it = item[0];
		switch(it.type){
		case PRE_SPLINEREF:
			ss << "Pre Spline " << it.num;//spline[it.num].name; 
			break;
		
		case PRE_CONSTSPLINEREF:
			ss << "Pre Const Spline " << it.num;//spline[it.num].name; 
			break;
			
		default:
			break;
		}
	}
	else{
		for(auto j = 0u; j < item.size(); j++){
			const auto &it = item[j];
	
			switch(it.type){
			case PRE_PARAMVEC:
				{
					const auto &pv = param_vec[it.num];
					ss << add_escape_char(param[pv.th].name) << pv.index;
				}
				break;
			
			case PRE_SPLINEREF:
				{
					const auto &spl = spline[it.num];
					ss << "Pre Spline " << add_escape_char(param[spl.th].name) << it.num;
				}
				//ss << "Spline " << it.num; //spline[it.num].name; 
				break;
			
			case PRE_CONSTSPLINEREF:
				ss << "Pre Const Spline " << it.num;//spline[it.num].name; 
				break;
			
			case PRE_ONE: ss << "1"; break;
			case PRE_REG: ss << "Rpre" << it.num; break;
			case PRE_REGTIME: ss << "Rpretime" << it.num; break;
			case PRE_NUMERIC: 
				{
					auto val = constant.value[it.num];
					if(val == INFY) ss << "INFY";
					else{
						if(val == UNDEF) ss << "UNDEF";
						else ss << val; 
					}
				}
				break;
			case PRE_TIME: ss << "time"; break;
			default: ss << it.type << endl; emsg_input("Precalc Eq "); break;
			}
			
			if(j != item.size()-1){
				switch(ca.op){
				case ADD: ss <<  "+"; break;
				case TAKE: ss <<  "-"; break;
				case MULTIPLY: ss << "*"; break;
				case SINGLE: ss <<  ""; break;
				case DIVIDE: ss <<  "/"; break;
				case POWERFUNC: ss << "|"; break;
				case THRESHFUNC: ss << "|"; break;
				case UBOUNDFUNC: ss << "|"; break;
				case MAXFUNC: ss << "|"; break;
				case MINFUNC: ss << "|"; break;
			
				default: break;
				}
			}
		}

		switch(ca.op){
		case EXPFUNC: case SINFUNC: case COSFUNC: case LOGFUNC: case STEPFUNC: case POWERFUNC: 
		case THRESHFUNC: case UBOUNDFUNC: case MAXFUNC: case MINFUNC: case ABSFUNC: case SQRTFUNC:
		case SIGFUNC:
			ss << ")"; 
			break;
		default: break;
		}
	}
		
	ss <<  " > ";
	ss <<  "Rpre" << q;
	
	if(ca.time_dep) ss << " Time dep";
	
	return ss.str();
}
 

/// Adds a new list_list to a specprecalc
unsigned int Precalc::add_list_time(SpecPrecalc &sprec, const vector <unsigned int> &ltime) const 
{
	if(ltime.size() == details.T) return ALL_TIME_STEP;
	
	auto &list_time = sprec.list_time;

	auto val = sum(ltime);
	
	auto j = sprec.hash_time.find(val);
	if(j == UNSET){
		j = list_time.size();
		sprec.hash_time.add(j,val);
		list_time.push_back(ltime);
		return j;
	}
	else{
		if(equal_vec(ltime,list_time[j])) return j;
		
		j = 0;
		while(j < list_time.size() && !equal_vec(ltime,list_time[j])) j++;	
		if(j == list_time.size()) list_time.push_back(ltime);
	}
	
	return j;
}


/// Converts from a spec_precalc to a map
vector <bool> Precalc::map_empty() const
{
	vector <bool> map(pcsize,false);
	return map;
}

	
/// Converts from a spec_precalc to a map
vector <bool> Precalc::map_from_spec(const SpecPrecalc &spec) const
{
	auto T = details.T;
	
	vector <bool> map(pcsize,false);
	for(const auto &inf : spec.info){
		auto q = inf.q;
		auto i = pcalcu[q].iref;
		auto k = inf.tlist;
		if(k == UNSET) map[i] = true;
		else{
			if(k == ALL_TIME_STEP){
				for(auto ti = 0u; ti < T; ti++) map[i+ti] = true;
			}
			else{
				for(auto ti : spec.list_time[k]) map[i+ti] = true;
			}
		}
	}
	
	return map;
}


/// Adds spec_precalc to a map
void Precalc::map_add_spec(vector <bool> &map, const SpecPrecalc &spec) const
{
	auto T = details.T;
	
	for(const auto &inf : spec.info){
		auto q = inf.q;
		auto i = pcalcu[q].iref;
		auto k = inf.tlist;
		if(k == UNSET) map[i] = true;
		else{
			if(k == ALL_TIME_STEP){
				for(auto ti = 0u; ti < T; ti++) map[i+ti] = true;
			}
			else{
				for(auto ti : spec.list_time[k]) map[i+ti] = true;
			}
		}
	}
}


/// Converts from a spec_precalc to a map
void Precalc::map_remove_spec(vector <bool> &map, const SpecPrecalc &spec) const
{
	auto T = details.T;
	
	for(const auto &inf : spec.info){
		auto q = inf.q;
		auto i = pcalcu[q].iref;
		auto k = inf.tlist;
		if(k == UNSET) map[i] = false;
		else{
			if(k == ALL_TIME_STEP){
				for(auto ti = 0u; ti < T; ti++) map[i+ti] = false;
			}
			else{
				for(auto ti : spec.list_time[k]) map[i+ti] = false;
			}
		}
	}
}


/// Converts from a spec_precalc to a map
const SpecPrecalc Precalc::spec_from_map(const vector <bool> &map) const
{
	SpecPrecalc spec;
	spec.hash.off();
	
	auto T = details.T;
	for(auto q = 0u; q < pcalcu.size(); q++){
		const auto &ca = pcalcu[q];
		auto i = ca.iref;
		if(ca.time_dep){
			vector <unsigned int> ltime;
			for(auto ti = 0u; ti < T; ti++){
				if(map[i+ti]) ltime.push_back(ti);
			}
			if(ltime.size() > 0){
				PrecalcInfo pi; pi.q = q; pi.tlist = add_list_time(spec,ltime);
				spec.info.push_back(pi);
			}
		}
		else{
			if(map[i]){
				PrecalcInfo pi; pi.q = q; pi.tlist = UNSET;
				spec.info.push_back(pi);
			}
		}
	}
	
	return spec;
}


/// Converts from a spec_precalc to a map
const SpecPrecalc Precalc::spec_from_ti_ref_list(vector <unsigned int> list, const vector < vector <unsigned int> > &ti_ref_list, const vector < vector <unsigned int> > &list_time) const
{
	auto T = details.T;
	
	SpecPrecalc spec;
	
	for(auto k = 0u; k < list.size(); k++){
		auto q = list[k];
		const auto &ca = pcalcu[q];
		
		if(ca.time_dep){
			const auto &trlist = ti_ref_list[q];
			switch(trlist.size()){
			case 0: emsg("SHould not be zero"); break;
			case 1:
				{
					auto j = trlist[0];
					if(j == UNSET) emsg("isss");
					if(j != ALL_TIME_STEP) j = add_list_time(spec,list_time[j]);
					PrecalcInfo pi; pi.q = q; pi.tlist = j;
					spec.hash.add(spec.info.size(),q);
					spec.info.push_back(pi);
				}
				break;
			
			default: // Combines together multiple sections
				{
					vector <bool> map(T,false);
					for(auto tr_ref : trlist){
						if(tr_ref == ALL_TIME_STEP){
							for(auto ti = 0u; ti < T; ti++) map[ti] = true;
						}
						else{
							for(auto ti : list_time[tr_ref]) map[ti] = true;
						}
					}
					
					vector <unsigned int> ltime;
					for(auto ti = 0u; ti < T; ti++){
						if(map[ti]) ltime.push_back(ti);
					}
					
					PrecalcInfo pi; pi.q = q; pi.tlist = add_list_time(spec,ltime);
					spec.hash.add(spec.info.size(),q);
					spec.info.push_back(pi);
				}
				break;
			}
		}		
		else{
			PrecalcInfo pi; pi.q = q; pi.tlist = UNSET;
			spec.hash.add(spec.info.size(),q);
			spec.info.push_back(pi);
		}
	}
	
	return spec;
}

	
/// Adds an element to a precalculation
void Precalc::sp_add(SpecPrecalc &sprec, unsigned int q, const vector <unsigned int> &ltime) const 
{
	auto &list_time = sprec.list_time;

	auto k = sprec.hash.find(q);
	if(k != UNSET){ // Already exists
		auto &in = sprec.info[k];
		
		if(pcalcu[q].time_dep){
			if(in.tlist != ALL_TIME_STEP){  // Need to combine together two timelines
				const auto &ltime_now = list_time[in.tlist];
				if(!equal_vec(ltime_now,ltime)){
					auto ltime_new = combine_list_time(ltime_now,ltime);
					in.tlist = add_list_time(sprec,ltime_new);
				}
			}
		}
	}
	else{
		PrecalcInfo pi; pi.q = q; 
		
		if(pcalcu[q].time_dep){
			pi.tlist = add_list_time(sprec,ltime);
		}
		else{
			pi.tlist = UNSET;
		}
		
		sprec.hash.add(sprec.info.size(),q);
		sprec.info.push_back(pi);
	}
}


/// Combines together two lists;
vector <unsigned int> Precalc::combine_list_time(const vector <unsigned int> &lt1,const vector <unsigned int> &lt2) const 
{
	auto T = details.T;
	vector <bool> map(T,false);
	for(auto i : lt1) map[i] = true;
	for(auto i : lt2) map[i] = true;
	
	vector <unsigned int> list;
	for(auto i = 0u; i < T; i++){
		if(map[i]) list.push_back(i);
	}
	
	for(auto i : lt1) cout << i << ","; 
	cout << "lt1" << endl;
	for(auto i : lt2) cout << i << ","; 
	cout << "lt2" << endl;
	for(auto i : list) cout << i << ",";
	cout << "list" << endl;
	emsg(" This is just testing if this part is working correctly");
	return list;
}


/* 
/// Adds a list of elements to a precalculation
void Precalc::sp_add(SpecPrecalc &sprec, const vector <unsigned int> &q_list, const vector <unsigned int> &ltime) const
{
	if(q_list.size() == 0) return;
	
	for(auto q : q_list) sp_add(sprec,q,ltime);
}
*/



bool PrecalcInfo_ord (const PrecalcInfo &pi1, const PrecalcInfo &pi2)                      
{ return (pi1.q < pi2.q); };


/// Sets up a spec precalc by combining several param_vec together
SpecPrecalc Precalc::combine_pv_spec_precalc(const vector <unsigned int> &param_list_tot) const 
{
	auto C = param_list_tot.size();
	
	if(C == 0){
		SpecPrecalc spec; 
		return spec;
	}
	
	if(C == 1){
		return param_vec[param_list_tot[0]].spec_precalc_after;
	}
	
	vector <bool> time_map(details.T,false);  

	// Combines all into a big list
	vector <PrecalcInfo> info; 
	vector < vector <unsigned int > > list_time; 
	
	for(auto k = 0u; k < C; k++){
		const auto &spa = param_vec[param_list_tot[k]].spec_precalc_after;
		
		auto lt_start = list_time.size();
		for(auto j = 0u; j < spa.list_time.size(); j++){
			list_time.push_back(spa.list_time[j]);
		}
		
		for(auto inf : spa.info){
			if(lt_start > 0){
				auto &k = inf.tlist;
				if(k != UNSET && k != ALL_TIME_STEP) k += lt_start;
			}
			info.push_back(inf);
		}
	}
	
	sort(info.begin(),info.end(),PrecalcInfo_ord);
		
	SpecPrecalc spec;
	auto jmax = info.size();
	auto j = 0u;
	while(j < jmax){
		auto j_st = j;
		auto inf = info[j];
	
		auto k = inf.tlist;
		auto q = inf.q;
		j++;
		while(j < jmax && info[j].q == q) j++;
		
		if(k == UNSET){ // Not time dependent
			for(auto jj = j_st+1; jj < j; jj++){
				if(info[jj].tlist != UNSET) emsg("Should all be unset");
			}
		}
		else{
			auto fl = false;
			for(auto jj = j_st; jj < j; jj++){
				if(info[jj].tlist == ALL_TIME_STEP) fl = true;
			}
			
			if(fl == true){ // At least one element shows all time
				inf.tlist = ALL_TIME_STEP;
			}
			else{           
				if(j == j_st+1){ // Only one timeline 
					inf.tlist = add_list_time(spec,list_time[k]);
				}
				else{            // Combine multiple time lines
					vector <unsigned int> time_list;
		
					for(auto jj = j_st; jj < j; jj++){
						for(auto ti : list_time[info[jj].tlist]){
							if(time_map[ti] == false){
								time_map[ti] = true;
								time_list.push_back(ti);
							}
						}	
					}

					sort(time_list.begin(),time_list.end());
					inf.tlist = add_list_time(spec,time_list);
					
					for(auto ti : time_list) time_map[ti] = false;
				}
			}
		}
		
		spec.info.push_back(inf);
	}

	//print_spec_precalc("after",spec);
	
	if(false){ // Checks if correct
		cout << "check map spec" << endl;
		for(auto ti = 0u; ti < details.T; ti++){
			if(time_map[ti] != false) emsg("Should not be false");
		}				

		auto map = map_from_spec(spec);
		
		auto map_true = map_empty();
		for(auto i = 0u; i < param_list_tot.size(); i++){
			const auto &pv2 = param_vec[param_list_tot[i]];
			map_add_spec(map_true,pv2.spec_precalc_after);
		}
		
		if(map.size() != map_true.size()) emsg("map size wrong");
		for(auto i = 0u; i < map.size(); i++){
			if(map[i] != map_true[i]) emsg("map problem precalc");
		}
	}
	
	return spec;
}


/// Sets up a spec precalc by combining several together
SpecPrecalc Precalc::combine_spec_precalc(const vector <SpecPrecalc> &spec_list) const 
{
	auto S = spec_list.size();
	
	if(S == 0){
		SpecPrecalc spec; 
		return spec;
	}
	
	if(S == 1){
		return spec_list[0];
	}
	
	vector <bool> time_map(details.T,false);  

	// Combines all into a big list
	vector <PrecalcInfo> info; 
	vector < vector <unsigned int > > list_time; 

	for(auto k = 0u; k < S; k++){
		const auto &spa = spec_list[k];
		auto lt_start = list_time.size();
		for(auto j = 0u; j < spa.list_time.size(); j++){
			list_time.push_back(spa.list_time[j]);
		}
		
		for(auto inf : spa.info){
			if(lt_start > 0){
				auto &k = inf.tlist;
				if(k != UNSET && k != ALL_TIME_STEP) k += lt_start;
			}
			info.push_back(inf);
		}
	}

	sort(info.begin(),info.end(),PrecalcInfo_ord);
		
	SpecPrecalc spec;
	auto jmax = info.size();
	auto j = 0u;
	while(j < jmax){
		auto j_st = j;
		auto inf = info[j];
	
		auto k = inf.tlist;
		auto q = inf.q;
		j++;
		while(j < jmax && info[j].q == q) j++;
		
		if(k == UNSET){ // Not time dependent
			for(auto jj = j_st+1; jj < j; jj++){
				if(info[jj].tlist != UNSET) emsg("Should all be unset");
			}
		}
		else{
			auto fl = false;
			for(auto jj = j_st; jj < j; jj++){
				if(info[jj].tlist == ALL_TIME_STEP) fl = true;
			}
			
			if(fl == true){ // At least one element shows all time
				inf.tlist = ALL_TIME_STEP;
			}
			else{           
				if(j == j_st+1){ // Only one timeline 
					inf.tlist = add_list_time(spec,list_time[k]);
				}
				else{            // Combine multiple time lines
					vector <unsigned int> time_list;
		
					for(auto jj = j_st; jj < j; jj++){
						for(auto ti : list_time[info[jj].tlist]){
							if(time_map[ti] == false){
								time_map[ti] = true;
								time_list.push_back(ti);
							}
						}	
					}

					sort(time_list.begin(),time_list.end());
					inf.tlist = add_list_time(spec,time_list);
					
					for(auto ti : time_list) time_map[ti] = false;
				}
			}
		}
		
		spec.info.push_back(inf);
	}

	//print_spec_precalc("after",spec);
	
	if(false){ // Checks if correct
		cout << "spec check" << endl;
		for(auto ti = 0u; ti < details.T; ti++){
			if(time_map[ti] != false) emsg("Should not be false");
		}				

		auto map = map_from_spec(spec);
		
		auto map_true = map_empty();
		for(const auto &spec : spec_list){
			map_add_spec(map_true,spec);
		}
		
		if(map.size() != map_true.size()) emsg("map size wrong");
		for(auto i = 0u; i < map.size(); i++){
			if(map[i] != map_true[i]) emsg("map problem precalc");
		}
	}
	
	return spec;
}

	
/// Sets up a spec precalc by combining several together and removing one
SpecPrecalc Precalc::combine_pv_spec_precalc_remove(const vector <unsigned int> &param_list_tot, const SpecPrecalc &remove) const 
{
	auto C = param_list_tot.size();
	if(C == 0){
		SpecPrecalc spec; 
		return spec;
	}

	const auto &pv = param_vec[param_list_tot[0]];
	
	auto map = map_from_spec(pv.spec_precalc_after);
	for(auto i = 1u; i < param_list_tot.size(); i++){
		const auto &pv2 = param_vec[param_list_tot[i]];
		map_add_spec(map,pv2.spec_precalc_after);
	}
	
	map_remove_spec(map,remove);
	
	return spec_from_map(map);
}


/*
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
*/


/// Shrinks to just a single time value (used for precalculation before reparameterise square spline
SpecPrecalc Precalc::shrink_sprec(unsigned int ti, const SpecPrecalc &spre_old) const
{
	vector <unsigned int> ltime; ltime.push_back(ti);
	
	SpecPrecalc spec;
	auto j = add_list_time(spec,ltime);
	for(auto pi : spre_old.info){
		if(pi.tlist != UNSET) pi.tlist = j;
		
		spec.hash.add(spec.info.size(),pi.q);
		spec.info.push_back(pi);
	}
	
	return spec;
}


/// Transfers parameter definition in precalc from spec_precalc_after into set_param_spec_precalc
void Precalc::set_param(SpecPrecalc &set_param_spec_precalc, SpecPrecalc &spec_precalc_after, bool spl_fl) const
{
	// Transfers over parameter
	auto &info = spec_precalc_after.info;

	if(info.size() == 0) emsg("zero size1");
	
	set_param_spec_precalc.info.push_back(info[0]);
	info.erase(info.begin());

	// Transfers the spline (if it exists)
	if(spl_fl){
		if(info.size() == 0) emsg("zero size2");
		const auto &in = info[0];	
		sp_add(set_param_spec_precalc,in.q,get_list_time(in.tlist,spec_precalc_after));
		info.erase(info.begin());
	}
}


// Sets precalculation to be done after sampling 
SpecPrecalc Precalc::calculate_spec_precalc_sample(const SpecPrecalc &spec_precalc) const
{
	auto map = map_from_spec(spec_precalc);
	
	for(const auto &pv : param_vec){
		if(pv.reparam_time_dep){
			map_remove_spec(map,pv.set_param_spec_precalc);
			map_remove_spec(map,pv.spec_precalc_after);
		}
	}
	
	for(const auto &spl : spline){
		if(spl.dynamic){
			for(auto q : spl.dynamic_precalc){
				auto i = pcalcu[q].iref;
				for(auto ti : all_time) map[i+ti] = false;
			}
		}
	}
	
	return spec_from_map(map);
}


/// Sets varaible all time
void Precalc::set_all_time()
{
	all_time = seq_vec(details.T);
}


/// Returns memory of hash_ca
double Precalc::hash_ca_mem() const
{
	return hash_ca.mem();	
}


/// Switches off hash table
void Precalc::hash_off()
{
	hash_ca.off();
}


/// Gets a string output from an equation type
void Precalc::eqn_type_error(PreEqItemType type, unsigned int ref) const
{
	string st;
	switch(type){	
	case PRE_PARAMVEC: st = "PRE_PARAMVEC"; break;
	case PRE_SPLINEREF: st = "PRE_SPLINEREF"; break;
	case PRE_CONSTSPLINEREF: st = "PRE_CONSTSPLINEREF"; break;
	default: st = "DEFAULT"+tstr(type); break;
	}
	
	emsg("EQN TYPE ERROR: "+st+" "+tstr(ref));
}


/// Sets the correct numbers for precalculation from q to i
void Precalc::adjust_num(vector <Calculation> &calc) const 
{
	for(auto &ca : calc){
		for(auto &it : ca.item) adjust_it(it);
	}
}


/// Sets the correct numbers for precalculation from q to i
void Precalc::adjust_it(EqItem &it) const 
{
	switch(it.type){
	case REG_PRECALC: case REG_PRECALC_TIME:
		it.num = pcalcu[it.num].iref;
		break;
	default: break;
	}
}

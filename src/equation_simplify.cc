// This file deals with ways to simplify equations

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

/// Simplifies based on numerical 
void Equation::simplify_operations(vector <EqItem> &op)
{
	auto loop = 0u;

	bool change, change2, change3;
	do{
		change = false;
		// Combines numeric values
		if(true){
			vector <EqItem> op_new;
		
			auto i = 0u;
			auto imax = op.size(); 
			while(i < imax){
				double addnum = UNSET;
				
				auto ty = op[i].type;
				switch(ty){
				case ADD:
					if(i > 0 && op[i-1].type == LEFTBRACKET){
						i++;
						addnum = CODE;
					}
					break;
					
				case LEFTBRACKET:
					if(i+2 < imax){
						if(op[i+1].type == NUMERIC && op[i+2].type == RIGHTBRACKET){
							if(!(i > 0 && (op[i-1].type == SUM || op[i-1].type == TINT))){
								op_new.push_back(op[i+1]);	
								i += 3;
								addnum = CODE;
							}
						}
					}
					break;
			
				case EXPFUNC: case SINFUNC: case COSFUNC: case LOGFUNC:    // Single parameter function 
				case ABSFUNC: case SQRTFUNC: case STEPFUNC: case SIGFUNC:
					if(i+3 < imax){
						if(op[i+1].type == LEFTBRACKET && op[i+2].type == NUMERIC && op[i+3].type == RIGHTBRACKET){
							auto con = constant.value[op[i+2].num];
							addnum = numeric_one_function(con,ty);
							i += 4;
						}
					}
					break;
					
				case POWERFUNC: case THRESHFUNC:            // Functions with two variables
				case UBOUNDFUNC: case MAXFUNC: case MINFUNC:
					if(i+5 < imax){
						if(op[i+1].type == LEFTBRACKET && op[i+2].type == NUMERIC && 
							op[i+3].type == FUNCDIVIDE && op[i+4].type == NUMERIC && op[i+5].type == RIGHTBRACKET){
							auto con1 = constant.value[op[i+2].num];
							auto con2 = constant.value[op[i+4].num];
							addnum = numeric_two_function(con1,con2,ty);
							i += 6;
						}
					}
					break;
			 
				case NUMERIC:
					if(i+2 < imax){
						switch(op[i+1].type){
						case ADD:
							if(op[i+2].type == NUMERIC){
								auto sum = constant.value[op[i].num];
								auto ist = i;
								while(ist+2 < op.size() && op[ist+1].type == ADD && op[ist+2].type == NUMERIC && !(ist+3 < op.size() && (op[ist+3].type == MULTIPLY || op[ist+3].type == DIVIDE))){
									sum += constant.value[op[ist+2].num];
									ist += 2;
								}
								
								if(ist != i){
									addnum = sum;
									i = ist+1;
								}
							}
							break;
							
						case MULTIPLY:
							if(op[i+2].type == NUMERIC){
								addnum = constant.value[op[i].num];
								while(i+2 < op.size() && op[i+1].type == MULTIPLY && op[i+2].type == NUMERIC){
									addnum *= constant.value[op[i+2].num];
									i += 2;
								}
								i += 1;
							}
							break;
							
						case DIVIDE:
							if(op[i+2].type == NUMERIC){
								auto con1 = constant.value[op[i].num];
								auto con2 = constant.value[op[i+2].num];
								if(con2 == 0){
									warn = "Equation contains a division by zero";
								}
								else{
									addnum = con1/con2;	
									i += 3;
								}
							}
							break;
		 
						default: break;
						}
					}
					
					if(addnum == UNSET){
						auto val = constant.value[op[i].num];
						if(val == 0){
							if(i+2 < imax){
								switch(op[i+1].type){
								case MULTIPLY:  // Removes 0*a*b
									{
										auto ist = i;
										while(ist+2 < imax && op[ist+1].type == MULTIPLY && quant(op,ist+2)) ist += 2;
										if(ist != i){
											op_new.push_back(op[i]);
											i = ist+1;
											addnum = CODE;
										}
										
										if(addnum == UNSET){ // Removes 0*(a*b)
											if(op[ist+2].type == LEFTBRACKET){
												auto iend = get_other_bracket(ist+2,op);
												if(iend != UNSET){
													op_new.push_back(op[i]);
													i = iend+1;
													addnum = CODE;
												}
											}
										}
									}
									break;
								
								case DIVIDE: // Removes 0/a
									if(quant(op,i+2)){
										op_new.push_back(op[i]);
										i += 3;
										addnum = CODE;
									}
									else{
										if(op[i+2].type == LEFTBRACKET){   // Removes 0/(a+b)
											auto iend = get_other_bracket(i+2,op);
											if(iend != UNSET){
												op_new.push_back(op[i]);
												i = iend+1;
												addnum = CODE;
											}
										}
									}
									break;
									
								default:
									break;
								}
							}	
						
							if(addnum == UNSET){
								auto len = op_new.size();
								if(len >= 2){
									switch(op_new[len-1].type){
									case MULTIPLY: // Removes a*0
										{
											auto ist = len;
											while(ist >= 2 && op_new[ist-1].type == MULTIPLY && quant(op_new,ist-2)) ist -= 2;
											if(ist != len){
												for(auto k = 0u; k < len-ist; k++){
													op_new.pop_back();
												}
												op_new.push_back(op[i]);
												i++;
												addnum = CODE;
											}
										
											if(addnum == UNSET){ // Removes (a*b)*0
												if(op_new[ist-2].type == RIGHTBRACKET){
													auto iend = get_other_bracket(ist-2,op_new);
													if(iend != UNSET){
														for(auto k = 0u; k < len-iend; k++){
															op_new.pop_back();
														}
														op_new.push_back(op[i]);
														i++;
														addnum = CODE;
													}
												}
											}
										}
										break;	
			
									default:
										break;									
									}
								}
							}
							
							if(addnum == UNSET){ 
								if(i+2 < imax){
									auto len = op_new.size();
									switch(op[i+1].type){
									case ADD:       // Removes 0+a
										if(!(len > 0 && (op_new[len-1].type == MULTIPLY || op_new[len-1].type == DIVIDE))){
											i += 2;
											addnum = CODE;
										}
										break;
									
									default: break;
									}
								}
							}
						
							if(addnum == UNSET){
								auto len = op_new.size();
								if(len >= 2){
									switch(op_new[len-1].type){
									case ADD:    // Removes a+0
										if(!(i+1 < imax && (op[i+1].type == MULTIPLY || op[i+1].type == DIVIDE))){
											op_new.pop_back();
											i++;
											addnum = CODE;
										}
										break;
									
									default:
										break;
									}
								}
							}
						}
						
						if(addnum == UNSET){
							if(val == 1){
								if(i+2 < imax){
									if(op[i+1].type == MULTIPLY){   // Removes 1*a
										i+=2;
										addnum = CODE;
									}
								}
								
								if(addnum == UNSET){
									auto len = op_new.size();
									if(len >= 2){
										if(op_new[len-1].type == MULTIPLY){   // Removes a*1
											op_new.pop_back();
											i++;
											addnum = CODE;
										}
									}
								}
							}
						}
					}
					break;
				
				case DIVIDE:
					if(i+1 < imax){ // Converts /numeric to *numeric
						if(op[i+1].type == NUMERIC){  
							auto con = constant.value[op[i+1].num];
							if(con == 0) warn = "Equation contains a division by zero";
							{
								EqItem ei; 
								ei.type = MULTIPLY;
								op_new.push_back(ei); 
							}
							i+=2;
							addnum = 1.0/con;
						}
					}
					break;
					
				default: break;
				}
				
				if(addnum == UNSET){
					op_new.push_back(op[i]);
					i++;
				}
				else{
					if(addnum != CODE){
						auto inew = op_new.size();
						
						// Removes brackets around numeric quantity
						if(inew > 0 && op_new[inew-1].type == LEFTBRACKET && i < imax && op[i].type == RIGHTBRACKET){
							if(!(inew > 1 && (op_new[inew-2].type == SUM || op_new[inew-2].type == TINT))){	
								op_new.pop_back(); i++;
							}
						}
						
						EqItem ei; 
						ei.type = NUMERIC;
						ei.num = constant.add(addnum);
						op_new.push_back(ei); 
					}
					change = true;
				}
			}
			
			if(change) op = op_new;
		}
		
		change2 = false;
		if(true){    // Looks to remove brackets
			vector <EqItem> op_new;
		
			auto i = 0u;
			auto imax = op.size(); 
			while(i < imax){
				auto fl = false;
				auto iend = UNSET;
					
				const auto &it = op[i];
				if(it.type == LEFTBRACKET && !(i+2 < imax && op[i+1].type == ONE && op[i+2].type == DIVIDE)){
					auto mult_fl = false;
					
					auto fl_left = false;
					if(i == 0) fl_left = true;
					else{
						switch(op[i-1].type){
						case LEFTBRACKET: case ADD: fl_left = true; break;
						case MULTIPLY: fl_left = true; mult_fl = true; break;
						default: break;
						}
					}
					
					if(fl_left){
						iend = get_other_bracket(i,op);
						if(iend != UNSET){
							auto fl_right = false;
							if(iend+1 == imax) fl_right = true;
							else{
								switch(op[iend+1].type){
								case RIGHTBRACKET: case ADD: fl_right = true; break;
								case MULTIPLY: fl_right = true; mult_fl = true; break;
								default: break;
								}
							}
							
							if(fl_right){
								fl = true;
								if(mult_fl == true){
									for(auto ii = i+1; ii < iend; ii++){
										switch(op[ii].type){
										case ADD: case TAKE: fl = false; break;
										default: break;
										}
									}									
								}
							}
						}						
					}
				}
				
				if(fl){
					for(auto ii = i+1; ii < iend; ii++){
						op_new.push_back(op[ii]);
					}
					change2 = true;
					i = iend+1;
				}
				else{
					op_new.push_back(it);
					i++;
				}
			}
			
			if(change2) op = op_new;
		}
		
		change3 = false;
		if(true){    // Looks to multiply out brackets a*(b+c) = (a*b+a*c). For popcomb
			vector <EqItem> op_new;
		
			auto i = 0u;
			auto imax = op.size(); 
			while(i < imax){
				auto fl = false;
					
				const auto &it = op[i];
				
				// a*(b+c) = (a*b+a*c)
				if(true && it.type == LEFTBRACKET && i > 1 && op[i-1].type == MULTIPLY){ 
					auto iend = get_other_bracket(i,op);
					if(iend != UNSET){
						if(contain_linear_pop(i+1,iend,op)){
							auto sec = eqn_split_add(i+1,iend,op);
							if(sec.size() > 1){
								auto ele = get_element_from_end(i-2,op);
								if(ele.start != UNSET && !contain_time(ele,op) && !contain_ie_fe(ele,op)){
									fl = true;
									change3 = true;
									
									for(auto ii = ele.start; ii <= ele.end; ii++) op_new.pop_back();
									
									add_it(LEFTBRACKET,op_new);
									
									for(auto k = 0u; k < sec.size(); k++){
										for(auto ii = ele.start; ii < ele.end; ii++) op_new.push_back(op[ii]);
										
										add_it(MULTIPLY,op_new);
										
										for(auto ii = sec[k].start; ii < sec[k].end; ii++) op_new.push_back(op[ii]);
										
										if(k+1 < sec.size()) add_it(ADD,op_new);
									}
									
									add_it(RIGHTBRACKET,op_new);
									
									i = iend+1; 
								}
							}
						}
					}
				}
				
				// (b+c)*a = (a*b+a*c)
				if(true && !fl && it.type == RIGHTBRACKET && i+2 < imax && op[i+1].type == MULTIPLY){ 
					auto istart = get_other_bracket(i,op);
					if(istart != UNSET){
						if(contain_linear_pop(istart+1,i,op)){
							auto sec = eqn_split_add(istart+1,i,op);
					
							if(sec.size() > 1){
								auto ele = get_element_from_start(i+2,op);
								if(ele.end != UNSET && !contain_time(ele,op) && !contain_ie_fe(ele,op)){
									fl = true;
									change3 = true;
									
									for(auto ii = istart; ii < i; ii++) op_new.pop_back();
									
									add_it(LEFTBRACKET,op_new);
									
									for(auto k = 0u; k < sec.size(); k++){
										for(auto ii = ele.start; ii < ele.end; ii++) op_new.push_back(op[ii]);
										
										add_it(MULTIPLY,op_new);
										
										for(auto ii = sec[k].start; ii < sec[k].end; ii++) op_new.push_back(op[ii]);
										
										if(k+1 < sec.size()) add_it(ADD,op_new);
									}
									
									add_it(RIGHTBRACKET,op_new);
									
									i = ele.end; 
								}
							}
						}
					}
				}
				
				// (b+c)/a = (b+c)*(1/a)
				if(true && !fl && it.type == RIGHTBRACKET && i+2 < imax && op[i+1].type == DIVIDE){ 
					auto ele = get_element_from_start(i+2,op);
					if(ele.end != UNSET){
						fl = true;
						change3 = true;
						
						add_it(RIGHTBRACKET,op_new);
						add_it(MULTIPLY,op_new);
						add_it(LEFTBRACKET,op_new);
						add_it(ONE,op_new);
						add_it(DIVIDE,op_new);
						for(auto ii = ele.start; ii < ele.end; ii++) op_new.push_back(op[ii]);
						add_it(RIGHTBRACKET,op_new);
									
						i = ele.end; 
					}
				}
				
				// Looks to pull factors out of a sum \sum (a(t)*(b+c)) = a(t)* \sum ((b+c))
				if(true && !fl && it.type == SUM && i+1 < imax && op[i+1].type == LEFTBRACKET){ 
					auto si = sum_info[op[i].num];
					
					auto iend = get_other_bracket(i+1,op);
					if(iend != UNSET){
						auto sec = eqn_split_mult(i+2,iend,op);
						if(sec.size() > 0){
							for(auto se : sec){
								print_eqn_range(se,op);
								cout << endl;
							}
							
							vector <unsigned int> list, list_not;
							for(auto k = 0u; k < sec.size(); k++){
								const auto &se = sec[k];
								if(contain_time(se,op) && !contain_dep(se,si.dep,op)) list.push_back(k);
								else list_not.push_back(k);
							}
							
							if(list.size() > 0){
								fl = true;
								change3 = true;
						
								for(auto k = 0u; k < list.size(); k++){
									const auto &se = sec[list[k]];
									if(k != 0) add_it(MULTIPLY,op_new);
									for(auto ii = se.start; ii < se.end; ii++) op_new.push_back(op[ii]);
								}
								add_it(MULTIPLY,op_new);
							
								op_new.push_back(it);
								
								add_it(LEFTBRACKET,op_new);
								if(list_not.size() == 0){
									add_it(ONE,op_new);
								}
								else{
									for(auto k = 0u; k < list_not.size(); k++){
										const auto &se = sec[list_not[k]];
										if(k != 0) add_it(MULTIPLY,op_new);
										for(auto ii = se.start; ii < se.end; ii++) op_new.push_back(op[ii]);
									}
								}
								add_it(RIGHTBRACKET,op_new);
								
								print_operations(op_new);
								i = iend+1; 
							}
						}
					}
				}
				
				if(!fl){
					op_new.push_back(it);
					i++;
				}
			}
			
			if(change3){
				op = op_new;
					print_operations(op_new);
						//emsg("HH");
			}
		}
	
		loop++;
	}while(change || change2 || change3);
}


/// Adds an operator item
void Equation::add_it(EqItemType type, vector <EqItem> &op) const
{
	EqItem ei; ei.type = type;
	op.push_back(ei);
}


/// Gets an element fro m the end point
EqnRange Equation::get_element_from_end(unsigned int end, const vector <EqItem> &op) const 
{
	EqnRange er;
	er.end = end+1;
	
	if(op[end].type == RIGHTBRACKET){
		auto start = get_other_bracket(end,op);
		if(start > 0 && is_func(op,start-1)) start--;
		er.start = start;
	}
	else{	
		er.start = end;
	}
	return er;
}


/// Gets an element from the start point
EqnRange Equation::get_element_from_start(unsigned int start, const vector <EqItem> &op) const 
{
	EqnRange er;
	er.start = start;
	if(is_func(op,start)){
		if(start+1 >= op.size() || op[start+1].type != LEFTBRACKET) er.end = UNSET;
		else{
			auto iend = get_other_bracket(start+1,op);
			if(iend == UNSET) er.end = UNSET;
			else er.end = iend+1;
		}
	}
	else{
		if(op[start].type == LEFTBRACKET){
			auto iend = get_other_bracket(start,op);
			if(iend == UNSET) er.end = UNSET;
			else er.end = iend+1;
		}
		else{
			er.end = start+1;
		}
	}
	return er;
}


/// Prints a range within an equation
void Equation::print_eqn_range(const EqnRange &er, const vector <EqItem> &op) const
{
	vector <EqItem> op2;
	for(auto i = er.start; i < er.end; i++) op2.push_back(op[i]);
	
	cout << "Range:";
	print_operations(op2);
}


/// Determines if linear population exist within a section
bool Equation::contain_time(const EqnRange &er, const vector <EqItem> &op) const 
{
	for(auto i = er.start; i < er.end; i++){
		switch(op[i].type){
		case POP_INDEX: case POPNUM: case POPTIMENUM: case TIME:
		case SPLINE: case SPLINEREF: case CONSTSPLINEREF: case REG_PRECALC_TIME:
		case POPCOMB: case POPCOMBTIME:
			return true;
			
		case PARAM_INDEX: 
			{
				const auto &pind = param_index[op[i].num];
				if(param[pind.th].time_dep) return true;
			}
			break;
		default: break;
		}
	}
	return false;
}


/// Determines if individual or fixed effect exist within a section
bool Equation::contain_ie_fe(const EqnRange &er, const vector <EqItem> &op) const 
{
	for(auto i = er.start; i < er.end; i++){
		switch(op[i].type){
		case IE: case FE: return true;
		default: return false;
		}
	}
	return false;
}


/// Determines if linear population exist within a section
bool Equation::contain_dep(const EqnRange &er, const vector <string> &dep, const vector <EqItem> &op) const 
{
	for(auto i = er.start; i < er.end; i++){
		switch(op[i].type){
		case POP_INDEX:
			{
				const auto &pind = pop_index[op[i].num];
				for(const auto &ins : pind.dep){
					if(find_in(dep,ins.index) != UNSET) return true;
				}
			}
			break;
		
		case PARAM_INDEX: 
			{
				const auto &pind = param_index[op[i].num];
				for(const auto &ins : pind.index_not_set){
					if(find_in(dep,ins.index) != UNSET) return true;
				}
			}
			break;
		default: break;
		}
	}
	return false;
}


/// Determines if linear population exist within a section
bool Equation::contain_linear_pop(unsigned int start, unsigned int end, const vector <EqItem> &op) const 
{
	auto num = 0u;
	for(auto i = start; i < end; i++){
		switch(op[i].type){
		case LEFTBRACKET: num++; break;
		case RIGHTBRACKET: num--; break;
		case POP_INDEX: if(num == 0) return true; break;
		default: break;
		}
	}
	return false;
}


/// Splits an section into pieces that are added up 
vector <EqnRange> Equation::eqn_split_add(unsigned int start, unsigned int end, const vector <EqItem> &op) const 
{
	vector <EqnRange> sec;
	
	auto num = 0u;
	auto ii = start;
	for(auto i = start; i < end; i++){
		switch(op[i].type){
		case LEFTBRACKET: num++; break;
		case RIGHTBRACKET: num--; break;
		case ADD:
			if(num == 0){
				EqnRange er;
				er.start = ii;
				er.end = i;
				sec.push_back(er);
				ii = i+1;
			}
			break;
		default: break;
		}
	}
	
	EqnRange er;
	er.start = ii;
	er.end = end;
	sec.push_back(er);

	return sec;
}


/// Splits an section into pieces that are added up 
vector <EqnRange> Equation::eqn_split_mult(unsigned int start, unsigned int end, const vector <EqItem> &op) const 
{
	{
		EqnRange er;
		er.start = start; er.end = end;
		print_eqn_range(er,op);
	}
	
	vector <EqnRange> sec;
	
	auto num = 0u;
	auto ii = start;
	for(auto i = start; i < end; i++){
		switch(op[i].type){
		case LEFTBRACKET: num++; break;
		case RIGHTBRACKET: num--; break;
		case DIVIDE:
			if(num == 0){
				sec.clear();
				return sec;
			}
			break;
		case MULTIPLY:
			if(num == 0){
				EqnRange er;
				er.start = ii;
				er.end = i;
				sec.push_back(er);
				ii = i+1;
			}
			break;
		default: break;
		}
	}
	
	EqnRange er;
	er.start = ii;
	er.end = end;
	sec.push_back(er);

	return sec;
}


/// Simplifies equations as much as possible (OLD CODE)
// (1) By removing zeros
// (2) Making used of infinity / undefined
// (3) Combining together addition and multiplication terms 
void Equation::simplify(vector <Calculation> &calc)
{
	return;
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
		
		//print_calc("calc "+tstr(loop),calc); cout << endl;
		
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
		
		// Converts DIVIDE by constant to MULTIPLY
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

					auto ty = ca.op;
					switch(ty){
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
					
					case EXPFUNC: case SINFUNC: case COSFUNC: case LOGFUNC:
					case STEPFUNC: case ABSFUNC: case SQRTFUNC: case SIGFUNC:
						if(ca.item[0].type == NUMERIC){
							auto con = cval[ca.item[0].num];
							rep_con = numeric_one_function(con,ty);
						}					
						break;
					
	
					case POWERFUNC: case THRESHFUNC:	case UBOUNDFUNC:
					case MAXFUNC: case MINFUNC:
						{
							double con1 = UNSET, con2 = UNSET;
							if(ca.item[0].type == NUMERIC) con1 = cval[ca.item[0].num];
							if(ca.item[1].type == NUMERIC) con2 = cval[ca.item[1].num];
							rep_con = numeric_two_function(con1,con2,ty);
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


/// Evalulates a function with one variable
double Equation::numeric_one_function(double con, EqItemType ty)
{
	switch(ty){
	case EXPFUNC:
		if(con == INFY) return INFY;
		else{
			if(con == UNDEF) return UNDEF;
			else return exp(con);
		}
		break;
		
	case SINFUNC:
		if(con == INFY) return UNDEF;
		else{
			if(con == UNDEF) return	UNDEF;
			else return sin(con);
		}
		break;
		
	case COSFUNC:
		if(con == INFY) return UNDEF;
		else{
			if(con == UNDEF) return UNDEF;
			else return cos(con);
		}	
		break;

	case LOGFUNC:
		if(con <= 0){
			warn = "Equation contains the logarithm of a non-positive";
		}
		else{
			if(con == INFY) return INFY;
			else{
				if(con == UNDEF) return UNDEF;
				else return log(con);	
			}
		}
		break;
		
	case STEPFUNC:
		if(con == UNDEF) return UNDEF;
		else{
			if(con > 0) return 1; 
			else return 0;
		}
		break;
		
	case ABSFUNC:
		if(con == UNDEF) return UNDEF;
		else{
			if(con < 0) return -con;
			else return con;
		}
		break;
							
	case SQRTFUNC:
		if(con < 0){
			warn = "Equation contains the square root of a non-positive";
		}
		else{
			if(con == UNDEF) return UNDEF;
			else{
				if(con == INFY) return INFY;
				else return sqrt(con);
			}
		}
		break;
		
	case SIGFUNC:
		if(con == UNDEF) return UNDEF;
		else{
			if(con == INFY) return 1;
			else return 1/(1+exp(-con));
		}
		break;
		
	default: emsg("Function not here"); break;
	}
	
	return UNSET;
}


/// Evalulates a function with one variable
double Equation::numeric_two_function(double con1, double con2, EqItemType ty)
{
	switch(ty){
	case POWERFUNC:
		if(con1 == UNDEF || con2 == UNDEF) return UNDEF;
		else{
			if(con2 == 0) return 1;
			else{
				if(con1 != UNSET){
					if(con1 == INFY || con2 == INFY) return INFY;
					else{	
						if(con1 == 0) return 0;
						else{
							if(con2 != UNSET) return pow(con1,con2);
						}
					}
				}
			}
		}
		break;
			
	case THRESHFUNC:
		if(con1 == UNDEF || con2 == UNDEF) return UNDEF;
		else{
			if(con1 < con2) return 0;
			else return con1;
		}
		break;
		
	case UBOUNDFUNC:
		if(con1 == UNDEF || con2 == UNDEF) return UNDEF;
		else{
			if(con1 > con2) return INFY;
			else return con1;
		}
		break;
		
	case MAXFUNC:
		if(con1 == UNDEF || con2 == UNDEF) return UNDEF;
		else{
			if(con1 > con2) return con1;
			else return con2;
		}
		break;
		
	case MINFUNC:
		if(con1 == UNDEF || con2 == UNDEF) return UNDEF;
		else{
			if(con1 < con2) return con1;
			else return con2;
		}
		break;
	
	default: emsg("Function not here"); break;
	}
	
	return UNSET;
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

/// Extracts a linear combination from 
void Equation::extract_popcomb(vector <Calculation> &calcu, Hash &hashw, Hash &hashpc, bool derive)
{
	//cout << "EXTRACT\n"; print_calculation();
	
	bool lc_exist = false;
	
	auto imax = calcu.size();
	vector < vector <unsigned int> > used;
	used.resize(imax);
	
	vector <bool> fl(imax,false); 
	vector <bool> remove(imax,false); 
	
	for(auto i = 0u; i < imax; i++){       // Determines how many times each resigter is used
		const auto &ca = calcu[i];
		
		vector <unsigned int> list;
		for(const auto &it : ca.item){
			if(it.type == REG){
				auto num = it.num;
				if(fl[num] == false){
					used[num].push_back(i);
					fl[num] = true;
					list.push_back(num);
				}
			}
		}
		
		for(auto num : list) fl[num] = false;
	}
	
	for(auto i = 0u; i < imax; i++){ 
		auto &ca = calcu[i];
		
		if(ca.op == ADD && ca.item.size() > 0){
			// Determines which elements can be added
			auto jmax = ca.item.size();
		
			vector <PopCombTemp> pop_co_temp;
			for(auto j = 0u; j < jmax; j++){
				const auto &it = ca.item[j];
				switch(it.type){
				case REG:
					{
						auto r = it.num;
						if(used[r].size() == 1){
							const auto &ca2 = calcu[r];
							if(ca2.op == MULTIPLY && ca2.item.size() == 2){
								const auto &it1 = ca2.item[0];
								const auto &it2 = ca2.item[1];
							
								if(it2.type == POPNUM){
									if(it1.type == NUMERIC || it2.type == REG_PRECALC){
										PopCombTemp pct;
										pct.po = it2.num;
										pct.it = it1;
										pct.r = r;
										pct.j = j;
										pop_co_temp.push_back(pct);
									}
								}
								else{
									if(it1.type == POPNUM){
										if(it2.type == NUMERIC || it2.type == REG_PRECALC){
											PopCombTemp pct;
											pct.po = it1.num;
											pct.it = it2;
											pct.r = r;
											pct.j = j;
											pop_co_temp.push_back(pct);
										}
									}
								}
							}
						}
					}
					break;
					
				case POPNUM:
					{
						PopCombTemp pct;
						pct.po = it.num;
						pct.it.type = ONE;
						pct.r = UNSET;
						pct.j = j;
						pop_co_temp.push_back(pct);
					}
					break;
					
				case POP_INDEX: 
					emsg("pop index"); // CHECKON
					break;
					
				default: break;
				}
			}	

			if(pop_co_temp.size() > 1){
				add_popcomb(i,pop_co_temp,remove,used,calcu,hashw,hashpc,derive);
				lc_exist = true;
			}				
		}
	}

	if(lc_exist){
		//print_calculation();
		
		// Removes unused registers
		vector <unsigned int> reg_map(imax,UNSET);
		
		vector <Calculation> calcu_new;
		for(auto i = 0u; i < imax; i++){
			if(remove[i] == false){
				reg_map[i] = calcu_new.size();
				auto ca = calcu[i];
				for(auto &it : ca.item){
					if(it.type == REG){
						it.num = reg_map[it.num];
						if(it.num == UNSET) emsg("SHould not be unset");					
					}
				}
				calcu_new.push_back(ca);
			}
		}
		
		calcu = calcu_new;
	}

	// Converts any other population to a popcomb
	for(auto &ca : calcu){
		for(auto &it : ca.item){
			if(it.type == POPNUM){
				add_popcomb_single(it,hashw,hashpc,derive);
			}
		}
	}
	
	//cout << " ex\n"; print_calculation();
	//emsg("H");
}

/// Adds a popcomb to the model
void Equation::add_popcomb(unsigned int i, const vector <PopCombTemp> &pop_co_temp, vector <bool> &remove, const vector < vector <unsigned int> > &used, vector <Calculation> &calcu, Hash &hashw, Hash &hashpc, bool derive)
{
	auto &item = calcu[i].item;
	
	auto jmax = item.size();
	vector <bool> jmap(jmax,false);
	
	vector <PopComb> pop_co;
	vector <unsigned int> vecpc;
	for(auto i = 0u; i < pop_co_temp.size(); i++){
		const auto &pct = pop_co_temp[i];
		
		auto po = pct.po;
		
		if(pct.r != UNSET) remove[pct.r] = true;
		jmap[pct.j] = true;
		
		const auto &it = pct.it;
		
		vector <unsigned int> vec;
		vec.push_back(it.type); vec.push_back(it.num);

		auto k = hashw.existing(vec);
		if(k == UNSET){
			k = popcombw.size();
			PopCombWeight pcw; 
			pcw.it = it;
			popcombw.push_back(pcw);
			
			hashw.add(k,vec);
		}
			
		PopComb pc;
		pc.po = po;
		pc.wref = k;
		pop_co.push_back(pc);
		
		vecpc.push_back(po); vecpc.push_back(k); 
	}
	
	auto pcnum = hashpc.existing(vecpc);
	if(pcnum == UNSET){
		pcnum = popcomb.size();
		popcomb.push_back(pop_co);
		
		hashpc.add(pcnum,vecpc);
	
		if(!derive){
			for(auto i = 0u; i < pop_co.size(); i++){
				auto po = pop_co[i].po;
				auto k = pop_co[i].wref;
				
				PopCombIn pci;
				pci.pc = pcnum;
				pci.index = i;
				popcombw[k].pcref.push_back(pci);
			
				PopCombRef pcr; 
				pcr.pcref = pcnum;
				pcr.wref = k;
				pop[po].popcomb_ref.push_back(pcr);
			}
		}
	}
	
	vector <EqItem> item_new;
	for(auto j = 0u; j < jmax; j++){
		if(jmap[j] == false) item_new.push_back(item[j]);
	}
	
	{
		EqItem it; it.type = POPCOMB; it.num = pcnum;
		item_new.push_back(it);
	}
	
	// If only popcomb left then substitute directly into other equations
	if(item_new.size() == 1 && i+1 < calcu.size()){ 
		if(used[i].size() == 0) emsg("There are no used");
		for(auto r : used[i]){
			for(auto &it : calcu[r].item){
				if(it.type == REG && it.num == i){
					it = item_new[0];
				}
			}
		}
		remove[i] = true;
	}
	else{
		item = item_new;
	}
}


/// Adds a popcomb to the model
void Equation::add_popcomb_single(EqItem &item, Hash &hashw, Hash &hashpc, bool derive)
{
	EqItem it; it.type = ONE; //it.num = constant.add(1);

	vector <unsigned int> vecpc;
	
	auto po = item.num;
	
	vector <unsigned int> vec;
	vec.push_back(it.type); vec.push_back(it.num);

	auto k = hashw.existing(vec);
	if(k == UNSET){
		k = popcombw.size();
		PopCombWeight pcw; 
		pcw.it = it;
		popcombw.push_back(pcw);
			
		hashw.add(k,vec);
	}
	
	vector <PopComb> pop_co;
	
	PopComb pc;
	pc.po = po;
	pc.wref = k;
	pop_co.push_back(pc);
	vecpc.push_back(po); vecpc.push_back(k); 
	
	auto pcnum = hashpc.existing(vecpc);
	if(pcnum == UNSET){
		pcnum = popcomb.size();
		popcomb.push_back(pop_co);
		hashpc.add(pcnum,vecpc);
		
		if(!derive){
			PopCombIn pci;
			pci.pc = pcnum;
			pci.index = 0;
			popcombw[k].pcref.push_back(pci);

			PopCombRef pcr; 
			pcr.pcref = pcnum;
			pcr.wref = k;
			pop[po].popcomb_ref.push_back(pcr);
		}
	}
	
	item.type = POPCOMB; item.num = pcnum;
}


/// Outputs popcomb name
string Equation::popcomb_name(unsigned int i) const
{
	stringstream ss;
	ss << "POPCOMB["; 
	const auto &pc = popcomb[i];
	auto jmax = pc.size(); if(jmax > 3) jmax = 3;
	for(auto j = 0u; j < jmax; j++){
		if(j != 0) ss << " + ";
		//print_item(popcombw[pc[j].wref].it);
		ss << "W" << pc[j].wref << "*";
		ss << pop[pc[j].po].name;
	}
	if(jmax == 3) ss << "...";	
	ss << "]";
	
	return ss.str();
}



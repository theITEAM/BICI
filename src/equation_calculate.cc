// This file deals with calculation of equations

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm> 

using namespace std;

#include "equation.hh"
#include "utils.hh"

/// Calculates a constant 
double Equation::calculate_constant() const 
{
	if(calcu.size() != 1) emsg("Constant problem1");
	const auto &ca = calcu[0];
	if(ca.op != ADD && ca.op != MULTIPLY) emsg("Constant problem2");
	if(ca.item.size() != 1) emsg("Constant problem3");
	
	const auto &it = ca.item[0];
	if(it.type != NUMERIC) emsg("Constant problem4");
	return constant.value[it.num];
}


/// Calculates the value for an equation using only parameters and splines (if ti_fix is set)
double Equation::calculate_param_ti_fix(const vector <double> &precalc) const 
{
	auto C = calcu.size();
 	vector <double> regcalc(C);

	const auto &cval = constant.value;
	
	auto calc_err = NO_ERROR;
	
  for(auto i = 0u; i < C; i++){
		const auto &ca = calcu[i];
		
		const auto &item = ca.item;
		const auto N = item.size();
	
		vector <double> num(N);
		
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
			
			switch(it.type){
				case REG: num[j] = regcalc[it.num]; break;
				case REG_PRECALC: num[j] = precalc[it.num]; break;
				case REG_PRECALC_TIME: num[j] = precalc[it.num+ti_fix]; break;
				case NUMERIC: num[j] = cval[it.num]; break;
				case CONSTSPLINEREF: num[j] = spline[it.num].const_val[ti_fix]; break;
				case TIME: num[j] = timepoint[ti_fix]; break;
				case ONE: num[j] = 1; break;
				case ZERO: num[j] = 0; break;
				default: eqn_type_error(it.type,2); break;
			}
		}

    regcalc[i] = calculate_operation(ca.op,num,calc_err);
  }

	if(calc_err != NO_ERROR) calc_error(calc_err,te_raw);
	
	return regcalc[C-1];
}


/// Calculates the value for an equation (but does not allow population number)
double Equation::calculate_no_popcomb(unsigned int ti, const vector <double> &precalc) const 
{
	auto C = calcu.size();
 	vector <double> regcalc(C);

	auto calc_err = NO_ERROR;

  for(auto i = 0u; i < calcu.size(); i++){
		const auto &ca = calcu[i];
		
		const auto &item = ca.item;
		const auto N = item.size();
		
		vector <double> num(N);
		
		const auto &cval = constant.value;
		
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
			
			switch(it.type){
				case ONE: num[j] = 1; break;
				case ZERO: num[j] = 0; break;
				case REG: num[j] = regcalc[it.num]; break;
				case REG_PRECALC: num[j] = precalc[it.num]; break;
				case REG_PRECALC_TIME: num[j] = precalc[it.num+ti]; break;
				case NUMERIC: num[j] = cval[it.num]; break;
				case CONSTSPLINEREF: num[j] = spline[it.num].const_val[ti]; break;	
				case TIME: num[j] = timepoint[ti]; break;
				default: eqn_type_error(it.type,3); break;
			}
		}

    regcalc[i] = calculate_operation(ca.op,num,calc_err);
  }

	if(calc_err != NO_ERROR) calc_error(calc_err,te_raw);
	
	return regcalc[C-1];
}


/// Calculate the value of an integral
double Equation::calculate_integral(unsigned int i, const vector < vector <double> > &popcomb_t, const vector <double> &precalc, const vector < vector < vector <double> > > &derive_val) const
{
	auto dt = details.dt;
	const auto &inte = integral[i];
	
	auto su = 0.0;
	if(calc_para_speedup){
		vector <unsigned int> list;
		for(auto ti = inte.ti_min; ti < inte.ti_max; ti++) list.push_back(ti);

		auto vec = calculate_para(inte.calc,list,popcomb_t,precalc,derive_val);
		su = sum(vec);  
	}
	else{
		for(auto ti = inte.ti_min; ti < inte.ti_max; ti++){
			auto val = calculate_calc(inte.calc,ti,popcomb_t[ti],precalc,derive_val);
		
			//print_working(inte.calc,ti,popcomb_t[ti],precalc,derive_val);
			//emsg("work");
			su += val;
		}		
	}
	
	return dt*su;
}


/// Calculates a derived value
double Equation::calculate_derive_item(unsigned int ti, unsigned int num, const vector < vector < vector <double> > > &derive_val) const
{
	const auto &dr = derive_ref[num];
	const auto &dv = derive_val[dr.i][dr.index];
	if(dv.size() == 1) return dv[0];  // Not time dependent
	else{
		if(dr.ti != UNSET) return dv[dr.ti];
		else{
			if(ti != UNSET) return dv[ti];
			else emsg("Derive problem");
		}
	}
	return UNSET;
}

					
/// Calculates derived equations
double Equation::calculate_derive(unsigned int ti, const vector < vector <double> > &popcomb_t, const vector <double> &precalc, const vector < vector < vector <double> > > &derive_val) const 
{
	auto C = calcu.size();
 	vector <double> regcalc(C);

	const auto &cval = constant.value;

	auto calc_err = NO_ERROR;

  for(auto i = 0u; i < C; i++){
		const auto &ca = calcu[i];
		
		const auto &item = ca.item;
		const auto N = item.size();
		
		vector <double> num(N);
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
			
			switch(it.type){
				case INTEGRAL:
					num[j] = calculate_integral(it.num,popcomb_t,precalc,derive_val);
					break;
					
				case DERIVE:
					num[j] = calculate_derive_item(ti,it.num,derive_val);
					break;
					
				case ONE: num[j] = 1; break;	
				case ZERO: num[j] = 0; break;	
				case POPCOMB:
					{
						if(ti == UNSET) emsg("ti should be set");
						num[j] = popcomb_t[ti][it.num]; 
					}
					break;
				case POPCOMBTIME:
					{
						const auto &ptr = popcomb_time_ref[it.num];
						num[j] = popcomb_t[ptr.ti][ptr.pc]; 
					}
					break;
		
				case REG: num[j] = regcalc[it.num]; break;
				case REG_PRECALC: num[j] = precalc[it.num]; break;
				case REG_PRECALC_TIME: num[j] = precalc[it.num+ti]; break;
				case NUMERIC: num[j] = cval[it.num]; break;
				case CONSTSPLINEREF: num[j] = spline[it.num].const_val[ti]; break;
				case TIME: num[j] = timepoint[ti]; break;
				default: eqn_type_error(it.type,151); break;
			}
		}

    regcalc[i] = calculate_operation(ca.op,num,calc_err);
  }

	if(calc_err != NO_ERROR) calc_error(calc_err,te_raw);
	
	return regcalc[C-1];
}


/// Calculates the value for an equation (including individual-based factor)
double Equation::calculate_indfac(const Individual &ind, unsigned int ti, const vector <double> &popcomb, const vector <double> &precalc) const 
{	
	return indfac(ind)*calculate(ti,popcomb,precalc);
}


/// Calculates the value for an equation
double Equation::calculate(unsigned int ti, const vector <double> &popcomb, const vector <double> &precalc) const 
{
	auto C = calcu.size();
 	vector <double> regcalc(C);

	const auto &cval = constant.value;

	auto calc_err = NO_ERROR;

  for(auto i = 0u; i < C; i++){
		const auto &ca = calcu[i];
		
		const auto &item = ca.item;
		const auto N = item.size();
		
		vector <double> num(N);
		
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
			
			switch(it.type){
				case ONE: num[j] = 1; break;
				case ZERO: num[j] = 0; break;
				case POPNUM: emsg("pnum not possible"); break;
				case POPNUMTIME: emsg("poptime num not possible"); break;
				case POPCOMB: num[j] = popcomb[it.num]; break;				
				case POPCOMBTIME: emsg("popcombtime num not possible"); break;
				case REG: num[j] = regcalc[it.num]; break;
				case REG_PRECALC: num[j] = precalc[it.num]; break;
				case REG_PRECALC_TIME: num[j] = precalc[it.num+ti]; break;
				case NUMERIC: num[j] = cval[it.num]; break;
				case CONSTSPLINEREF: num[j] = spline[it.num].const_val[ti]; break;
				case TIME: num[j] = timepoint[ti]; break;
				default: eqn_type_error(it.type,22); break;
			}
		}

    regcalc[i] = calculate_operation(ca.op,num,calc_err);
  }
		
	if(calc_err != NO_ERROR) calc_error(calc_err,te_raw);
	
	return regcalc[C-1];
}


/// Calculates the value for an equation
double Equation::calculate_reg(unsigned int ti, const vector <double> &popcomb, const vector <double> &precalc, const vector <unsigned int> &list, vector <double> &regcalc) const 
{
	const auto &cval = constant.value;
	
	auto calc_err = NO_ERROR;
	
	for(auto i : list){
		const auto &ca = calcu[i];
		
		const auto &item = ca.item;
		const auto N = item.size();
		
		vector <double> num(N);
		
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
			
			switch(it.type){
				case ONE: num[j] = 1; break;
				case ZERO: num[j] = 0; break;
				case POPNUM: emsg("pnum not possible"); break;
				case POPNUMTIME: emsg("poptime num not possible"); break;
				case POPCOMB: num[j] = popcomb[it.num]; break;				
				case POPCOMBTIME: emsg("popcombtime num not possible"); break;
				case REG: num[j] = regcalc[it.num]; break;
				case REG_PRECALC: num[j] = precalc[it.num]; break;
				case REG_PRECALC_TIME: num[j] = precalc[it.num+ti]; break;
				case NUMERIC: num[j] = cval[it.num]; break;
				case CONSTSPLINEREF: num[j] = spline[it.num].const_val[ti]; break;
				case TIME: num[j] = timepoint[ti]; break;
				default: eqn_type_error(it.type,18); break;
			}
		}

    regcalc[i] = calculate_operation(ca.op,num,calc_err);
  }

	if(calc_err != NO_ERROR) calc_error(calc_err,te_raw);
	
	return regcalc[regcalc.size()-1];
}


/// Calculates the value for an equation (but allows for other times for population to be accessed)
double Equation::calculate_all_time(unsigned int ti, const vector < vector <double> > &popcomb_t, const vector <double> &precalc) const 
{
	auto C = calcu.size();
 	vector <double> regcalc(C);

	const auto &cval = constant.value;

	auto calc_err = NO_ERROR;

  for(auto i = 0u; i < C; i++){
		const auto &ca = calcu[i];
		
		const auto &item = ca.item;
		const auto N = item.size();
		
		vector <double> num(N);
		
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
			
			switch(it.type){
				case ONE: num[j] = 1; break;
				case ZERO: num[j] = 0; break;
				case POPCOMB: num[j] = popcomb_t[ti][it.num]; break;
				case POPCOMBTIME:
					{
						const auto &ptr = popcomb_time_ref[it.num];
						num[j] = popcomb_t[ptr.ti][ptr.pc]; 
					}
					break;
				case REG: num[j] = regcalc[it.num]; break;
				case REG_PRECALC: num[j] = precalc[it.num]; break;
				case REG_PRECALC_TIME: num[j] = precalc[it.num+ti]; break;
				case NUMERIC: num[j] = cval[it.num]; break;
				case CONSTSPLINEREF: num[j] = spline[it.num].const_val[ti]; break;
				case TIME: num[j] = timepoint[ti]; break;
				default: eqn_type_error(it.type,19); break;
			}
		}

    regcalc[i] = calculate_operation(ca.op,num,calc_err);
  }

	if(calc_err != NO_ERROR) calc_error(calc_err,te_raw);
	
	return regcalc[C-1];
}

/// Calculates the value for an equation
double Equation::calculate_calc(const vector <Calculation> &calc, unsigned int ti, const vector <double> &popcomb, const vector <double> &precalc, const vector < vector < vector <double> > > &derive_val) const
{
	auto C = calc.size();
 	vector <double> regcalc(C);

	const auto &cval = constant.value;

	auto calc_err = NO_ERROR;

  for(auto i = 0u; i < C; i++){
		const auto &ca = calc[i];
		
		const auto &item = ca.item;
		const auto N = item.size();
		
		vector <double> num(N);
		
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
			
			switch(it.type){
				case DERIVE:
					{
						const auto &dr = derive_ref[it.num];
						const auto &dv = derive_val[dr.i][dr.index];
						if(dv.size() == 1) num[j] = dv[0];  // Not time dependent
						else{
							if(dr.ti != UNSET) num[j] = dv[dr.ti];
							else{
								if(ti != UNSET) num[j] = dv[ti];
								else emsg("Derive problem");
							}
						}
					}
					break;
						
				case ONE: num[j] = 1; break;
				case ZERO: num[j] = 0; break;
				case POPNUM: emsg("Cannot calculate pnum");  break;
				case POPNUMTIME: emsg("Cannot calculate poptimenum"); break;
				case POPCOMB: num[j] = popcomb[it.num]; break;
				case POPCOMBTIME: emsg("Cannot calculate poptimenum"); break;
				case REG: num[j] = regcalc[it.num]; break;
				case REG_PRECALC: num[j] = precalc[it.num]; break;
				case REG_PRECALC_TIME: num[j] = precalc[it.num+ti]; break;
				case NUMERIC: num[j] = cval[it.num]; break;
				case CONSTSPLINEREF: num[j] = spline[it.num].const_val[ti]; break;
				case TIME: num[j] = timepoint[ti]; break;
				default: eqn_type_error(it.type,20); break;
			}
		}
		
    regcalc[i] = calculate_operation(ca.op,num,calc_err);
  }

	if(calc_err != NO_ERROR) calc_error(calc_err,te_raw);
	
	return regcalc[C-1];
}


/// Calculates the value for an equation
double Equation::print_working(const vector <Calculation> &calc, unsigned int ti, const vector <double> &popcomb, const vector <double> &precalc, const vector < vector < vector <double> > > &derive_val) const
{
	auto C = calc.size();
 	vector <double> regcalc(C);

	const auto &cval = constant.value;

	auto imax = C; if(imax > 100) imax = 100;
	
	auto calc_err = NO_ERROR;
	
  for(auto i = 0u; i < imax; i++){
		const auto &ca = calc[i];
		
		const auto &item = ca.item;
		const auto N = item.size();
		
		vector <double> num(N);
		
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
			
			switch(it.type){
				case DERIVE:
					{
						const auto &dr = derive_ref[it.num];
						const auto &dv = derive_val[dr.i][dr.index];
						if(dv.size() == 1) num[j] = dv[0];  // Not time dependent
						else{
							if(dr.ti != UNSET) num[j] = dv[dr.ti];
							else{
								if(ti != UNSET) num[j] = dv[ti];
								else emsg("Derive problem");
							}
						}
					}
					break;
						
				case ONE: num[j] = 1; break;
				case ZERO: num[j] = 0; break;
				case POPNUM: emsg("Cannot calculate pnum"); break;
				case POPNUMTIME: emsg("Cannot calculate poptimenum"); break;
				case POPCOMB: num[j] = popcomb[it.num]; break;				
				case POPCOMBTIME: emsg("popcombtime num not possible"); break;
				case REG: num[j] = regcalc[it.num]; break;
				case REG_PRECALC: num[j] = precalc[it.num]; break;
				case REG_PRECALC_TIME: num[j] = precalc[it.num+ti]; break;
				case NUMERIC: num[j] = cval[it.num]; break;
				case CONSTSPLINEREF: num[j] = spline[it.num].const_val[ti]; break;
				case TIME: num[j] = timepoint[ti]; break;
				default: eqn_type_error(it.type,21); break;
			}
		}
		
		
    regcalc[i] = calculate_operation(ca.op,num,calc_err);
		
		print_ca(i,ca); 
		cout << " = " << regcalc[i] << " ca" <<  endl;
		for(auto va :num) cout << va << ","; 
		cout << "Num" << endl;
  }

	if(calc_err != NO_ERROR) calc_error(calc_err,te_raw);
	
	return regcalc[C-1];
}


/// Tests that equation_para is working
void Equation::test_calculate_para(const vector <Calculation> &calc, const vector <unsigned int> &list, const vector < vector <double> > &popcomb_t, const vector <double> &precalc, const vector < vector < vector <double> > > &derive_val, string ref) const
{
	auto vec = calculate_para(calc,list,popcomb_t,precalc,derive_val);
	
	auto j = 0u;
	for(auto ti : list){
		auto val = calculate_calc(calc,ti,popcomb_t[ti],precalc,derive_val);
		if(dif(vec[j],val,DIF_THRESH)){
			emsg("Calculate para problem"+ref);
		}
		j++;
	}
}


/// Calculates the value for an equation in parallel across several times
vector <double> Equation::calculate_para(const vector <Calculation> &calc, const vector <unsigned int> &list, const vector < vector <double> > &popcomb_t, const vector <double> &precalc, const vector < vector < vector <double> > > &derive_val) const 
{
	/*
	if(lin.on){
		auto value_para = calculate_linear_list_derive(list,popcomb_t,precalc,derive_val);
		
		for(auto &va : value_para){
			if(va < 0 && va > -TINY) va = 0;
		}
		
		return value_para;
	}
	*/
	
	auto K = list.size();
	auto C = calc.size();
 	vector < vector <double> > regcalc;
	regcalc.resize(C);
	
	const auto &cval = constant.value;

  for(auto i = 0u; i < C; i++){
		regcalc[i].resize(K);
		auto &rc = regcalc[i];
		
		const auto &ca = calc[i];
		const auto &item = ca.item;
		const auto N = item.size();
		
		switch(ca.op){
		case MULTIPLY:
			{
				vector <unsigned int> tdj;
				
				auto val = 1.0;
				for(auto j = 0u; j < N; j++){
					const auto &it = item[j];
						
					switch(it.type){
						case DERIVE:	
							{
								if(derive_val.size() == 0) emsg("Derive val does not exist");
								const auto &dr = derive_ref[it.num];
								const auto &dv = derive_val[dr.i][dr.index];
								if(dv.size() == 1) val *= dv[0];  // Not time dependent
								else{
									if(dr.ti != UNSET) val *= dv[dr.ti];
									else tdj.push_back(j);
								}
							}
							break;
							
						case ONE: break;
						case ZERO: val = 0; break;
						
						case REG_PRECALC: val *= precalc[it.num]; break;
						
						case NUMERIC: val *= cval[it.num]; break;
						
						case POPCOMB: case POPCOMBTIME: 
						case REG: case REG_PRECALC_TIME: case CONSTSPLINEREF: case TIME: 
							tdj.push_back(j);
							break;
							
						default:
							eqn_type_error(it.type,180); 
							break;
					}
				}
				
				if(tdj.size() == 0){               // If no time dependency then set all to same
					for(auto k = 0u; k < K; k++) rc[k] = val;
				}
				else{
					auto fl = true;                // For the first it sets the value and later modify
					for(auto j : tdj){
						const auto &it = item[j];
					
						switch(it.type){
						case DERIVE:	
							{
								const auto &dr = derive_ref[it.num];
								const auto &dv = derive_val[dr.i][dr.index];
							
								if(fl){
									if(val == 1){ for(auto k = 0u; k < K; k++) rc[k] = dv[list[k]];}
									else{ for(auto k = 0u; k < K; k++) rc[k] = val*dv[list[k]];}
								}
								else{ for(auto k = 0u; k < K; k++) rc[k] *= dv[list[k]];}
							}								
							break;
							
						case POPNUM: case POPNUMTIME: emsg("should not be pnum"); break;
						
						case POPCOMB: 
							{
								auto itn = it.num;
								if(fl){
									if(val == 1){ for(auto k = 0u; k < K; k++) rc[k] = popcomb_t[list[k]][itn];}
									else{ for(auto k = 0u; k < K; k++) rc[k] = val*popcomb_t[list[k]][itn];}
								}
								else{ for(auto k = 0u; k < K; k++) rc[k] *= popcomb_t[list[k]][itn];}
							}
							break;
							
						case POPCOMBTIME: 
							{
								const auto &ptr = popcomb_time_ref[it.num];
								auto value = popcomb_t[ptr.ti][ptr.pc];
								if(fl){
									if(val == 1){ for(auto k = 0u; k < K; k++) rc[k] = value;}
									else{ for(auto k = 0u; k < K; k++) rc[k] = val*value;}
								}
								else{ for(auto k = 0u; k < K; k++) rc[k] *= value;}
							}
							break;
							
						case REG: 
							{
								const auto &regc = regcalc[it.num];  
								if(fl){
									if(val == 1){ for(auto k = 0u; k < K; k++) rc[k] = regc[k];}
									else{ for(auto k = 0u; k < K; k++) rc[k] = val*regc[k];}
								}
								else{ for(auto k = 0u; k < K; k++) rc[k] *= regc[k];}
							}
							break;
						
						case REG_PRECALC_TIME: 
							{
								auto itn = it.num;
								if(fl){ 
									if(val == 1){ for(auto k = 0u; k < K; k++) rc[k] = precalc[itn+list[k]];}
									else{ for(auto k = 0u; k < K; k++) rc[k] = val*precalc[itn+list[k]];}
								}
								else{ for(auto k = 0u; k < K; k++) rc[k] *= precalc[itn+list[k]];}
							}
							break;
							
						case CONSTSPLINEREF: 
							{			
								const auto &cv = spline[it.num].const_val;
								if(fl){ 
									if(val == 1){ for(auto k = 0u; k < K; k++) rc[k] = cv[list[k]];}
									else{ for(auto k = 0u; k < K; k++) rc[k] = val*cv[list[k]];}
								}
								else{ for(auto k = 0u; k < K; k++) rc[k] *= cv[list[k]];}
							}
							break;
							
						case TIME:
							{
								if(fl){ 
									if(val == 1){ for(auto k = 0u; k < K; k++) rc[k] = timepoint[list[k]];}
									else{ for(auto k = 0u; k < K; k++) rc[k] = val*timepoint[list[k]];}
								}
								else{ for(auto k = 0u; k < K; k++) rc[k] *= timepoint[list[k]];}
							}
							break;
							
						default: break;
						}
						fl = false;
					}
				}
			}	
			break;

		case ADD:
			{
				vector <unsigned int> tdj;
				
				auto val = 0.0;
				for(auto j = 0u; j < N; j++){
					const auto &it = item[j];
						
					switch(it.type){
						case DERIVE:	
							{
								if(derive_val.size() == 0) emsg("Derive val does not exist");
								const auto &dr = derive_ref[it.num];
								const auto &dv = derive_val[dr.i][dr.index];
								if(dv.size() == 1) val += dv[0];  // Not time dependent
								else{
									if(dr.ti != UNSET) val += dv[dr.ti];
									else tdj.push_back(j);
								}
							}
							break;
							
						case ONE: val += 1; break;
						case ZERO: break;
						
						case REG_PRECALC: val += precalc[it.num]; break;
						
						case NUMERIC: val += cval[it.num]; break;
						
						case POPNUM: case POPNUMTIME: emsg("not pn"); break;
						
						case POPCOMB: case POPCOMBTIME: 
						case REG: case REG_PRECALC_TIME: case CONSTSPLINEREF: case TIME: 
							tdj.push_back(j);
							break;
							
						default:
							eqn_type_error(it.type,181); 
							break;
					}
				}
				
				if(tdj.size() == 0){               // If no time dependency then set all to same
					for(auto k = 0u; k < K; k++) rc[k] = val;
				}
				else{
					auto fl = true;                // For the first it sets the value and later modify
					for(auto j : tdj){
						const auto &it = item[j];
					
						switch(it.type){
						case DERIVE:	
							{
								const auto &dr = derive_ref[it.num];
								const auto &dv = derive_val[dr.i][dr.index];
							
								if(fl){ 
									if(val == 0){ for(auto k = 0u; k < K; k++) rc[k] = dv[list[k]];}
									else{ for(auto k = 0u; k < K; k++) rc[k] = val+dv[list[k]];}
								}
								else{ for(auto k = 0u; k < K; k++) rc[k] += dv[list[k]];}
							}								
							break;
							
						case POPNUM: case POPNUMTIME: emsg("not pn"); break;
							
						case POPCOMB: 
							{
								auto itn = it.num;
								if(fl){
									if(val == 0){ for(auto k = 0u; k < K; k++) rc[k] = popcomb_t[list[k]][itn];}
									else{ for(auto k = 0u; k < K; k++) rc[k] = val+popcomb_t[list[k]][itn];}
								}
								else{ for(auto k = 0u; k < K; k++) rc[k] += popcomb_t[list[k]][itn];}
							}
							break;
							
						case POPCOMBTIME: 
							{
								const auto &ptr = popcomb_time_ref[it.num];
								auto value = popcomb_t[ptr.ti][ptr.pc];
								if(fl){
									if(val == 1){ for(auto k = 0u; k < K; k++) rc[k] = value;}
									else{ for(auto k = 0u; k < K; k++) rc[k] = val+value;}
								}
								else{ for(auto k = 0u; k < K; k++) rc[k] += value;}
							}
							break;
							
						case REG: 
							{
								const auto &regc = regcalc[it.num];  
								if(fl){ 
									if(val == 0){ for(auto k = 0u; k < K; k++) rc[k] = regc[k];}
									else{ for(auto k = 0u; k < K; k++) rc[k] = val+regc[k];}
								}
								else{ for(auto k = 0u; k < K; k++) rc[k] += regc[k];}
							}
							break;
						
						case REG_PRECALC_TIME: 
							{
								auto itn = it.num;
								if(fl){
									if(val == 0){ for(auto k = 0u; k < K; k++) rc[k] = precalc[itn+list[k]];}
									else{ for(auto k = 0u; k < K; k++) rc[k] = val+precalc[itn+list[k]];}
								}
								else{ for(auto k = 0u; k < K; k++) rc[k] += precalc[itn+list[k]];}
							}
							break;
							
						case CONSTSPLINEREF: 
							{			
								const auto &cv = spline[it.num].const_val;
								if(fl){
									if(val == 0){ for(auto k = 0u; k < K; k++) rc[k] = cv[list[k]];}
									else{ for(auto k = 0u; k < K; k++) rc[k] = val+cv[list[k]];}
								}
								else{ for(auto k = 0u; k < K; k++) rc[k] += cv[list[k]];}
							}
							break;
							
						case TIME:
							{
								if(fl){ 
									if(val == 0){ for(auto k = 0u; k < K; k++) rc[k] = timepoint[list[k]];}
									else{ for(auto k = 0u; k < K; k++) rc[k] = val+timepoint[list[k]];}
								}
								else{ for(auto k = 0u; k < K; k++) rc[k] += timepoint[list[k]];}
							}
							break;
							
						default: break;
						}
						
						fl = false;
					}
				}
			}	
			break;

		default:
			{
				auto kmax = N*K;
				vector <double> num(kmax);
			
				for(auto j = 0u; j < N; j++){
					const auto &it = item[j];
					
					double val = UNSET;
					switch(it.type){
						case DERIVE:	
							{
								if(derive_val.size() == 0) emsg("Derive val does not exist");
								const auto &dr = derive_ref[it.num];
								const auto &dv = derive_val[dr.i][dr.index];
								if(dv.size() == 1) val = dv[0];  // Not time dependent
								else{
									if(dr.ti != UNSET) val = dv[dr.ti];
									else{
										auto k = j;
										for(auto ti : list){ num[k] = dv[ti]; k += N;}
									}
								}
							}
							break;
						case ONE: val = 1; break;
						case ZERO: val = 0; break;
						case REG_PRECALC: val = precalc[it.num]; break;
						case NUMERIC: val = cval[it.num]; break;
						case POPNUM: case POPNUMTIME: emsg("not pnptn"); break;
						case POPCOMB: 
							{
								auto itn = it.num;
								auto k = j;
								for(auto ti : list){ num[k] = popcomb_t[ti][itn]; k += N;}
							}
							break;
							
						case POPCOMBTIME: 
							{
								const auto &ptr = popcomb_time_ref[it.num];
								auto value = popcomb_t[ptr.ti][ptr.pc];
								auto k = j;
								for(auto ii = 0u; ii < list.size(); ii++){ num[k] = value; k += N;}
								//for(auto ti : list){ num[k] = value; k += N;}
							}
							break;
							
						case REG: 
							{
								const auto &regc = regcalc[it.num];
								auto kk = j;
								for(auto k = 0u; k < K; k++){
									num[kk] = regc[k];
									kk += N;
								}
							}
							break;
						
						case REG_PRECALC_TIME: 
							{
								auto itn = it.num;
								auto k = j;
								for(auto ti : list){ num[k] = precalc[itn+ti]; k +=N;}
							}
							break;
							
						case CONSTSPLINEREF:
							{			
								const auto &cv = spline[it.num].const_val;
								auto k = j;
								for(auto ti : list){ num[k] = cv[ti]; k += N;}	
							}
							break;
							
						case TIME:
							{
								auto k = j;
								for(auto ti : list){ num[k] = timepoint[ti]; k += N;}	
							}
							break;
							
						default:
							eqn_type_error(it.type,182); 
							break;
					}
					
					if(val != UNSET){
						for(auto k = j; k < kmax; k += N) num[k] = val;
					}		
				}
					
				auto op = ca.op;
				
				for(auto k = 0u; k < K; k++){
					double *nu = &num[k*N]; 
					
					switch(op){
					case ADD:
						{
							auto ans = 0.0; 
							for(auto j = 0u; j < N; j++) ans += nu[j];
							rc[k] = ans;
						}
						break;
					
					case MULTIPLY:
						{
							auto ans = 1.0; 
							for(auto j = 0u; j < N; j++) ans *= nu[j];
							rc[k] = ans;
						}
						break;
						
					case TAKE: 
						rc[k] = nu[0]-nu[1]; 
						break;
						
					case DIVIDE:
						if(nu[1] == 0) run_error("Equation caused a division by zero1."); 
						rc[k] = nu[0]/nu[1]; 
						break;
						
					case EXPFUNC: 
						rc[k] = exp(nu[0]); 
						break; 			
					
					case SINFUNC: 
						rc[k] = sin(nu[0]); 
						break; 
					
					case COSFUNC: 
						rc[k] = cos(nu[0]); 
						break; 
					
					case LOGFUNC: 
						if(nu[0] <= 0) run_error("The quantity inside a log function became negative."); 
						rc[k] = log(nu[0]);
						break;
						
					case STEPFUNC:
						if(nu[0] > 0) rc[k] = 1;
						else rc[k] = 0; 
						break;
						
					case POWERFUNC:
						rc[k] = pow(nu[0],nu[1]);
						break;
						
					case THRESHFUNC:
						if(nu[0] < nu[1]) rc[k] = 0;
						else rc[k] = nu[0];
						break;
						
					case UBOUNDFUNC:
						if(nu[0] > nu[1]) rc[k] = INFINITY;
						else rc[k] = nu[0];
						break;
						
					case MAXFUNC:
						if(nu[0] > nu[1]) rc[k] = nu[0];
						else rc[k] = nu[1]; 
						break;
						
					case MINFUNC:
						if(nu[0] < nu[1]) rc[k] = nu[0];
						else rc[k] = nu[1];
						break;
						
					case ABSFUNC:
						if(nu[0] > 0) rc[k] = nu[0]; 
						else rc[k] = -nu[0];
						break;
						
					case SQRTFUNC:
						if(nu[0] < 0) run_error("A square root of a negative number was found in an equation."); 
						rc[k] = sqrt(nu[0]);
						break;
						
					case SIGFUNC:
						rc[k] = 1/(1+exp(-nu[0]));
						break;
				
					default: emsg("Equation error7"); break;
					}
				}
			}
			break;
		}
  }
	
	auto &value_para = regcalc[C-1];
	
	for(auto &va : value_para){
		if(va < 0 && va > -TINY) va = 0;
	}
	
	return value_para;
}


/// Calculates the value for an equation which only depends on parameters
double Equation::calculate_param(const vector <double> &precalc) const 
{
	auto C = calcu.size();
 	vector <double> regcalc(C);

	const auto &cval = constant.value;

	auto calc_err = NO_ERROR;

  for(auto i = 0u; i < C; i++){
		const auto &ca = calcu[i];
		
		const auto &item = ca.item;
		const auto N = item.size();
		
		vector <double> num(N);
		
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
			
			switch(it.type){
				case ONE: num[j] = 1; break;
				case ZERO: num[j] = 0; break;
				case REG: num[j] = regcalc[it.num]; break;
				case REG_PRECALC: num[j] = precalc[it.num]; break;
				case NUMERIC: num[j] = cval[it.num]; break;
				default: eqn_type_error(it.type,183); break;
			}
		}

    regcalc[i] = calculate_operation(ca.op,num,calc_err);
  }

	if(calc_err != NO_ERROR) calc_error(calc_err,te_raw);

	return regcalc[C-1];
}


/// Calculates the value for an equation which only depends on numbers (e.g. prior hyperparamters)
double Equation::calculate_value() const 
{
	auto C = calcu.size();
 	vector <double> regcalc(C);

	const auto &cval = constant.value;

	auto calc_err = NO_ERROR;

  for(auto i = 0u; i < C; i++){
		const auto &ca = calcu[i];
		
		const auto &item = ca.item;
		const auto N = item.size();
		
		vector <double> num(N);
		
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
			
			switch(it.type){
				case ONE: num[j] = 1; break;
				case ZERO: num[j] = 0; break;
				case REG: num[j] = regcalc[it.num]; break;
				case NUMERIC: num[j] = cval[it.num]; break;
				default: eqn_type_error(it.type,183); break;
			}
		}

    regcalc[i] = calculate_operation(ca.op,num,calc_err);
  }

	if(calc_err != NO_ERROR) calc_error(calc_err,te_raw);

	return regcalc[C-1];
}


/// Calculates a population gradient
double Equation::calculate_pop_grad(unsigned int pref, unsigned int ti, const vector <double> &popcombw_value, const vector <double> &precalc) const
{
	const auto &pgr = lin.pop_grad_ref[pref];
	if(pgr.size() == 1){
		return popcombw_value[pgr[0].wref]*calculate_item(pgr[0].popcomb_grad,ti,precalc);
	}
	else{
		auto sum = 0.0;
		for(const auto &pg : pgr){
			sum += popcombw_value[pg.wref]*calculate_item(pg.popcomb_grad,ti,precalc);
		}
		return sum;
	}
}


/// Calculates the term with no population
double Equation::calculate_no_pop(unsigned int ti, const vector <double> &precalc) const
{
	return calculate_item(lin.no_pop_precalc,ti,precalc);
}


/// Uses linearity to calculate equation value
double Equation::calculate_linear(unsigned int ti, const vector <double> &popcomb, const vector <double> &precalc) const
{
	auto val = calculate_item(lin.no_pop_precalc,ti,precalc);
	
	const auto &pgp = lin.popcomb_grad_precalc;
	for(auto i = 0u; i < pgp.size(); i++){
		val += calculate_item(pgp[i],ti,precalc)*popcomb[popcomb_ref[i]];
	}
	if(val < 0 && val > -TINY) val = 0;

	return val;
}


/// Uses linearity to calculate equation value at multiple time points
vector <double> Equation::calculate_linear_list(const vector <unsigned int> &list, const vector < vector <double> > &popcomb_t, const vector <double> &precalc) const
{
	auto N = list.size();
	
	vector <double> val(N);
	
	{
		const auto &it = lin.no_pop_precalc;
		if(item_timedep(it)){
			for(auto i = 0u; i < N; i++) val[i] = calculate_item(it,list[i],precalc);
		}
		else{
			auto va = calculate_item_no_time(it,precalc);
			for(auto i = 0u; i < N; i++) val[i] = va;
		}
	}
	
	const auto &pgp = lin.popcomb_grad_precalc;
	for(auto j = 0u; j < pgp.size(); j++){
		auto pc = popcomb_ref[j];
		
		const auto &it = pgp[j];
		if(item_timedep(it)){
			for(auto i = 0u; i < N; i++){
				auto ti = list[i];
				val[i] += calculate_item(it,ti,precalc)*popcomb_t[ti][pc];
			}
		}
		else{
			auto va = calculate_item_no_time(it,precalc);
			for(auto i = 0u; i < N; i++){
				val[i] += va*popcomb_t[list[i]][pc];
			}
		}
	}
	
	for(auto &va : val){ 
		if(va < 0 && va > -TINY) va = 0;
	}
	
	return val;
}


/// Calculates a single item
double Equation::calculate_item(const EqItem &it, unsigned int ti, const vector <double> &precalc) const
{
	switch(it.type){
	case ONE: return 1;
	case ZERO: return 0;
	case REG_PRECALC: return precalc[it.num];
	case REG_PRECALC_TIME: return precalc[it.num+ti];
	case NUMERIC: return constant.value[it.num]; 
	case CONSTSPLINEREF: return spline[it.num].const_val[ti];		
	case TIME: return timepoint[ti]; 
	default: eqn_type_error(it.type,12); break;
	}
	
	return UNSET;
}


/// Calculates a single item
bool Equation::item_timedep(const EqItem &it) const
{
	switch(it.type){
	case REG_PRECALC_TIME: case CONSTSPLINEREF: case TIME: return true;
	default: return false;
	}
}


/// Calculates a single item
double Equation::calculate_item_no_time(const EqItem &it, const vector <double> &precalc) const
{
	switch(it.type){
	case ONE: return 1;
	case ZERO: return 0;
	case REG_PRECALC: return precalc[it.num];
	case NUMERIC: return constant.value[it.num]; 
	default: eqn_type_error(it.type,13); break;
	}
	
	return UNSET;
}


/// Calculates a single item (but uses the old value if it is available
double Equation::calculate_item_old(const EqItem &it, unsigned int ti, const vector <double> &precalc, const vector <double> &precalc_old) const
{
	switch(it.type){
	case ONE: return 1;
	case ZERO: return 0;	
	case REG_PRECALC: 
		{
			auto val = precalc_old[it.num];
			if(val != UNSET) return val;
			return precalc[it.num];
		}
		break;
	case REG_PRECALC_TIME: 
		{
			auto val = precalc_old[it.num+ti];
			if(val != UNSET) return val;
			return precalc[it.num+ti];
		}
		break;
	case NUMERIC: return constant.value[it.num]; 
	case CONSTSPLINEREF: return spline[it.num].const_val[ti];		
	case TIME: return timepoint[ti]; 
	default: eqn_type_error(it.type,14); break;
	}
	
	return UNSET;
}


// Calculates a single item (but uses the old value if it is available
double Equation::calculate_item_old_no_time(const EqItem &it, const vector <double> &precalc, const vector <double> &precalc_old) const
{
	switch(it.type){
	case ONE: return 1;
	case ZERO: return 0;	
	case REG_PRECALC: 
		{
			auto val = precalc_old[it.num];
			if(val != UNSET) return val;
			return precalc[it.num];
		}
		break;
	case NUMERIC: return constant.value[it.num]; 
	default: eqn_type_error(it.type,15); break;
	}
	
	return UNSET;
}


/// Calculates equation for checking
double Equation::calculate_check(unsigned int ti, const PV &param_val, const vector < vector <double> > &popnum_t, const vector < vector < vector <double> > > &derive_val, const vector <Calculation> &calcu, const vector <Integral> &integral, const vector <ParamRef> &param_ref, const vector <DeriveRef> &derive_ref, const vector <unsigned long long> &spline_refi) const 
{
	auto C = calcu.size();
 	vector <double> regcalc(C);

	const auto &cval = constant.value;

	auto calc_err = NO_ERROR;

  for(auto i = 0u; i < C; i++){
		const auto &ca = calcu[i];
		
		const auto &item = ca.item;
		const auto N = item.size();
		
		vector <double> num(N);
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
			
			switch(it.type){
				case INTEGRAL:
					{
						auto dt = details.dt;
						const auto &inte = integral[it.num];
	
						auto su = 0.0;
						for(auto ti = inte.ti_min; ti < inte.ti_max; ti++){
							auto val = calculate_check(ti,param_val,popnum_t,derive_val,inte.calc,integral,param_ref,derive_ref,spline_refi);
							su += val;
						}		
						num[j] =dt*su;
					}			
					break;
					
				case DERIVE:
					{
						const auto &dr = derive_ref[it.num];
						const auto &dv = derive_val[dr.i][dr.index];
						if(dv.size() == 1) num[j] = dv[0];  // Not time dependent
						else{
							if(dr.ti != UNSET) num[j] = dv[dr.ti];
							else{
								if(ti != UNSET) num[j] = dv[ti];
								else emsg("Derive problem");
							}
						}
					}
					break;
					
				case POPNUM:
					{
						if(ti == UNSET) emsg("ti should be set");
						num[j] = popnum_t[ti][it.num]; 
						if(num[j] < 0) num[j] = 0;
					}
					break;
		
				case POPNUMTIME:
					{
						const auto &ptr = pop_time_ref[it.num];
						num[j] = popnum_t[ptr.ti][ptr.po]; 
						if(num[j] < 0) num[j] = 0;
					}
					break;
					
				case PARAMETER:	
					{
						const auto &pr = param_ref[it.num];
						const auto &par = param[pr.th];
						const auto &er = par.element_ref[pr.index];
						if(er.cons) num[j] = cval[er.index];
						else{
							const auto &ele = par.element[er.index];
							num[j] = param_val.value[ele.param_vec_ref];
						}
					}
					break;
					
				case SPLINE:
					{
						const auto &pr = param_ref[it.num];
						const auto &par = param[pr.th];
						const auto &er = par.element_ref[pr.index];
						if(er.cons){
							auto s = 0u; 
							for(s = 0u; s < spline.size(); s++){
								 const auto &spl = spline[s];
								 const auto &sinfo = spl.info;
								 if(spl.th == pr.th && spl.index*sinfo.knot_tdiv.size() == pr.index) break;
							}			
							if(s == spline.size()) emsg("cannot find spline");

							num[j] = spline[s].const_val[ti];
						}
						else{
							const auto &ele = par.element[er.index];
							const auto &pv = param_vec[ele.param_vec_ref];
							num[j] = param_val.precalc[spline_refi[pv.spline_ref]+ti];
						}
					}
					break;
					
				case ONE: num[j] = 1; break;	
				case ZERO: num[j] = 0; break;	
				case REG: num[j] = regcalc[it.num]; break;
				case NUMERIC: num[j] = cval[it.num]; break;
				case CONSTSPLINEREF: num[j] = spline[it.num].const_val[ti]; break;
				case TIME: num[j] = timepoint[ti]; break;
		
				case REG_PRECALC: emsg("no precalc"); break;
				case REG_PRECALC_TIME: emsg("no precalc"); break;
				case POPCOMB: emsg("Should not be popcomb"); break;
				case POPCOMBTIME: emsg("Should not be popcombtime"); break;		
				default: eqn_type_error(it.type,150); break;
			}
		}

    regcalc[i] = calculate_operation(ca.op,num,calc_err);
  }

	if(calc_err != NO_ERROR) calc_error(calc_err,te_raw);
	
	return regcalc[C-1];
}

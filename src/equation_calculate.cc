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
				default: eqn_type_error(it.type,2); break;
			}
		}

    regcalc[i] = calculate_operation(ca.op,num);
  }

	return regcalc[C-1];
}


/// Calculates the value for an equation (but does not allow population number)
double Equation::calculate_no_popnum(unsigned int ti, const vector <double> &precalc) const 
{
	auto C = calcu.size();
 	vector <double> regcalc(C);

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
				case REG: num[j] = regcalc[it.num]; break;
				case REG_PRECALC: num[j] = precalc[it.num]; break;
				case REG_PRECALC_TIME: num[j] = precalc[it.num+ti]; break;
				case NUMERIC: num[j] = cval[it.num]; break;
				case CONSTSPLINEREF: num[j] = spline[it.num].const_val[ti]; break;	
				case TIME: num[j] = timepoint[ti]; break;
				default: eqn_type_error(it.type,3); break;
			}
		}

    regcalc[i] = calculate_operation(ca.op,num);
  }

	return regcalc[C-1];
}


/// Calculate the value of an integral
double Equation::calculate_integral(unsigned int i, const vector < vector <double> > &popnum_t, const vector <double> &precalc, const vector < vector < vector <double> > > &derive_val) const
{
	auto dt = details.dt;
	const auto &inte = integral[i];
	
	auto su = 0.0;
	if(calc_para_speedup){
		vector <unsigned int> list;
		for(auto ti = inte.ti_min; ti < inte.ti_max; ti++) list.push_back(ti);

		auto vec = calculate_para(inte.calc,list,popnum_t,precalc,derive_val);
		su = sum(vec);  
	}
	else{
		for(auto ti = inte.ti_min; ti < inte.ti_max; ti++){
			auto val = calculate_calc(inte.calc,ti,popnum_t[ti],precalc,derive_val);
			su += val;
		}		
	}
	
	return dt*su;
}


/// Calculates derived equations
double Equation::calculate_derive(unsigned int ti, const vector < vector <double> > &popnum_t, const vector <double> &precalc, const vector < vector < vector <double> > > &derive_val) const 
{
	auto C = calcu.size();
 	vector <double> regcalc(C);

	const auto &cval = constant.value;

  for(auto i = 0u; i < C; i++){
		const auto &ca = calcu[i];
		
		const auto &item = ca.item;
		const auto N = item.size();
		
		vector <double> num(N);
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
			
			switch(it.type){
				case INTEGRAL:
					num[j] = calculate_integral(it.num,popnum_t,precalc,derive_val);
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
					
				case ONE: num[j] = 1; break;	
				case POPNUM:
					{
						if(ti == UNSET) emsg("ti should be set");
						num[j] = rectify(popnum_t[ti][it.num]); 
					}
					break;
				case POPTIMENUM:
					{
						const auto &ptr = pop_time_ref[it.num];
						num[j] = rectify(popnum_t[ptr.ti][ptr.po]); 
					}
					break;
				case REG: num[j] = regcalc[it.num]; break;
				case REG_PRECALC: num[j] = precalc[it.num]; break;
				case REG_PRECALC_TIME: num[j] = precalc[it.num+ti]; break;
				case NUMERIC: num[j] = cval[it.num]; break;
				case CONSTSPLINEREF: num[j] = spline[it.num].const_val[ti]; break;
				case TIME: num[j] = timepoint[ti]; break;
				default: eqn_type_error(it.type,5); break;
			}
		}

    regcalc[i] = calculate_operation(ca.op,num);
  }

	return regcalc[C-1];
}


/// Calculates the value for an equation (including individual-based factor)
double Equation::calculate_indfac(const Individual &ind, unsigned int ti, const vector <double> &popnum, const vector <double> &precalc) const 
{	
	return indfac(ind)*calculate(ti,popnum,precalc);
}


/// Calculates the value for an equation
double Equation::calculate(unsigned int ti, const vector <double> &popnum, const vector <double> &precalc) const 
{
	auto C = calcu.size();
 	vector <double> regcalc(C);

	const auto &cval = constant.value;

  for(auto i = 0u; i < C; i++){
		const auto &ca = calcu[i];
		
		const auto &item = ca.item;
		const auto N = item.size();
		
		vector <double> num(N);
		
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
			
			switch(it.type){
				case ONE: num[j] = 1; break;
				case POPNUM: num[j] = rectify(popnum[it.num]); break;
				case POPTIMENUM: emsg("poptime num not possible"); break;
				case REG: num[j] = regcalc[it.num]; break;
				case REG_PRECALC: num[j] = precalc[it.num]; break;
				case REG_PRECALC_TIME: num[j] = precalc[it.num+ti]; break;
				case NUMERIC: num[j] = cval[it.num]; break;
				case CONSTSPLINEREF: num[j] = spline[it.num].const_val[ti]; break;
				case TIME: num[j] = timepoint[ti]; break;
				default: eqn_type_error(it.type,17); break;
			}
		}

    regcalc[i] = calculate_operation(ca.op,num);
  }

	return regcalc[C-1];
}


/// Calculates the value for an equation (but allows for other times for population to be accessed)
double Equation::calculate_all_time(unsigned int ti, const vector < vector <double> > &popnum_t, const vector <double> &precalc) const 
{
	auto C = calcu.size();
 	vector <double> regcalc(C);

	const auto &cval = constant.value;

  for(auto i = 0u; i < C; i++){
		const auto &ca = calcu[i];
		
		const auto &item = ca.item;
		const auto N = item.size();
		
		vector <double> num(N);
		
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
			
			switch(it.type){
				case ONE: num[j] = 1; break;
				case POPNUM: num[j] = rectify(popnum_t[ti][it.num]); break;
				case POPTIMENUM: 
					{
						const auto &ptr = pop_time_ref[it.num];
						if(ti < ptr.ti) emsg("Cannot calculate population time before it is set");
						num[j] = rectify(popnum_t[ptr.ti][ptr.po]); 
					}
					break;
				case REG: num[j] = regcalc[it.num]; break;
				case REG_PRECALC: num[j] = precalc[it.num]; break;
				case REG_PRECALC_TIME: num[j] = precalc[it.num+ti]; break;
				case NUMERIC: num[j] = cval[it.num]; break;
				case CONSTSPLINEREF: num[j] = spline[it.num].const_val[ti]; break;
				case TIME: num[j] = timepoint[ti]; break;
				default: eqn_type_error(it.type,17); break;
			}
		}

    regcalc[i] = calculate_operation(ca.op,num);
  }

	return regcalc[C-1];
}

/// Calculates the value for an equation
double Equation::calculate_calc(const vector <Calculation> &calc, unsigned int ti, const vector <double> &popnum, const vector <double> &precalc, const vector < vector < vector <double> > > &derive_val) const
{
	auto C = calc.size();
 	vector <double> regcalc(C);

	const auto &cval = constant.value;

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
				case POPNUM: num[j] = rectify(popnum[it.num]); break;
				case POPTIMENUM: emsg("Cannot calculate poptimenum"); break;
				case REG: num[j] = regcalc[it.num]; break;
				case REG_PRECALC: num[j] = precalc[it.num]; break;
				case REG_PRECALC_TIME: num[j] = precalc[it.num+ti]; break;
				case NUMERIC: num[j] = cval[it.num]; break;
				case CONSTSPLINEREF: num[j] = spline[it.num].const_val[ti]; break;
				case TIME: num[j] = timepoint[ti]; break;
				default: eqn_type_error(it.type,17); break;
			}
		}
		
    regcalc[i] = calculate_operation(ca.op,num);
  }

	return regcalc[C-1];
}


/// Tests that equation_para is working
void Equation::test_calculate_para(const vector <Calculation> &calc, const vector <unsigned int> &list, const vector < vector <double> > &popnum_t, const vector <double> &precalc, const vector < vector < vector <double> > > &derive_val, string ref) const
{
	auto vec = calculate_para(calc,list,popnum_t,precalc,derive_val);
	
	auto j = 0u;
	for(auto ti : list){
		auto val = calculate_calc(calc,ti,popnum_t[ti],precalc,derive_val);
		if(dif(vec[j],val,DIF_THRESH)){
			emsg("Calculate para problem"+ref);
		}
		j++;
	}
}


/// Calculates the value for an equation in parallel across several times
vector <double> Equation::calculate_para(const vector <Calculation> &calc, const vector <unsigned int> &list, const vector < vector <double> > &popnum_t, const vector <double> &precalc, const vector < vector < vector <double> > > &derive_val) const 
{
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
						
						case REG_PRECALC: val *= precalc[it.num]; break;
						
						case NUMERIC: val *= cval[it.num]; break;
						
						case POPNUM: case POPTIMENUM: case REG: case REG_PRECALC_TIME: case CONSTSPLINEREF: case TIME: 
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
							
						case POPNUM: 
							{
								auto itn = it.num;
								if(fl){
									if(val == 1){ for(auto k = 0u; k < K; k++) rc[k] = rectify(popnum_t[list[k]][itn]);}
									else{ for(auto k = 0u; k < K; k++) rc[k] = val*rectify(popnum_t[list[k]][itn]);}
								}
								else{ for(auto k = 0u; k < K; k++) rc[k] *= rectify(popnum_t[list[k]][itn]);}
							}
							break;
							
						case POPTIMENUM: 
							{
								const auto &ptr = pop_time_ref[it.num];
								auto value = rectify(popnum_t[ptr.ti][ptr.po]);
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
						
						case REG_PRECALC: val += precalc[it.num]; break;
						
						case NUMERIC: val += cval[it.num]; break;
						
						case POPNUM: case POPTIMENUM: case REG: case REG_PRECALC_TIME: case CONSTSPLINEREF: case TIME: 
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
							
						case POPNUM: 
							{
								auto itn = it.num;
								if(fl){
									if(val == 0){ for(auto k = 0u; k < K; k++) rc[k] = rectify(popnum_t[list[k]][itn]);}
									else{ for(auto k = 0u; k < K; k++) rc[k] = val+rectify(popnum_t[list[k]][itn]);}
								}
								else{ for(auto k = 0u; k < K; k++) rc[k] += rectify(popnum_t[list[k]][itn]);}
							}
							break;
							
						case POPTIMENUM: 
							{
								const auto &ptr = pop_time_ref[it.num];
								auto value = rectify(popnum_t[ptr.ti][ptr.po]);
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
						case REG_PRECALC: val = precalc[it.num]; break;
						case NUMERIC: val = cval[it.num]; break;
						
						case POPNUM: 
							{
								auto itn = it.num;
								auto k = j;
								for(auto ti : list){ num[k] = rectify(popnum_t[ti][itn]); k += N;}
							}
							break;
							
						case POPTIMENUM: 
							{
								const auto &ptr = pop_time_ref[it.num];
								auto value = rectify(popnum_t[ptr.ti][ptr.po]);
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
	
	return regcalc[C-1];
}


/// Calculates the value for an equation which only depends on parameters
double Equation::calculate_param(const vector <double> &precalc) const 
{
	auto C = calcu.size();
 	vector <double> regcalc(C);

	const auto &cval = constant.value;

  for(auto i = 0u; i < C; i++){
		const auto &ca = calcu[i];
		
		const auto &item = ca.item;
		const auto N = item.size();
		
		vector <double> num(N);
		
		for(auto j = 0u; j < N; j++){
			const auto &it = item[j];
			
			switch(it.type){
				case ONE: num[j] = 1; break;
				case REG: num[j] = regcalc[it.num]; break;
				case REG_PRECALC: num[j] = precalc[it.num]; break;
				case NUMERIC: num[j] = cval[it.num]; break;
				default: eqn_type_error(it.type,183); break;
			}
		}

    regcalc[i] = calculate_operation(ca.op,num);
  }

	return regcalc[C-1];
}


/// Calculates a population gradient
double Equation::calculate_pop_grad(unsigned int pref, unsigned int ti, const vector <double> &precalc) const
{
	return calculate_item(linearise.factor_precalc,ti,precalc)*
	       calculate_item(linearise.pop_grad_precalc[pref],ti,precalc);	
}


/// Calculates a population gradient
double Equation::calculate_pop_grad_no_time(unsigned int pref, const vector <double> &precalc) const
{
	return calculate_item_no_time(linearise.factor_precalc,precalc)*
	       calculate_item_no_time(linearise.pop_grad_precalc[pref],precalc);	
}


/// Calculates a population gradient without factor term
double Equation::calculate_pop_grad_without_factor(unsigned int pref, const vector <double> &precalc) const
{
	return calculate_item_no_time(linearise.pop_grad_precalc[pref],precalc);	
}


/// Calculates factor part of population gradient
double Equation::calculate_factor(unsigned int ti, const vector <double> &precalc) const
{
	return calculate_item(linearise.factor_precalc,ti,precalc);
}
      

/// Calculates the term with no population
double Equation::calculate_no_pop(unsigned int ti, const vector <double> &precalc) const
{
	return calculate_item(linearise.no_pop_precalc,ti,precalc);
}


/// Calculates a single item
double Equation::calculate_item(const EqItem &it, unsigned int ti, const vector <double> &precalc) const
{
	switch(it.type){
	case ZERO: return 0;
	case ONE: return 1;
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
double Equation::calculate_item_no_time(const EqItem &it, const vector <double> &precalc) const
{
	switch(it.type){
	case ZERO: return 0;
	case ONE: return 1;
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
	case ZERO: return 0;
	case ONE: return 1;
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
	case ZERO: return 0;
	case ONE: return 1;
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



"use strict";
// Functions which evaluate equations

const LEFTBRACKET=0, RIGHTBRACKET=1, FUNCDIVIDE=2, ADD=3, TAKE=4, MULTIPLY=5, DIVIDE=6, REG=7, TIME=11, NUMERIC=15, EXPFUNC=16, SINFUNC=17, COSFUNC=18, LOGFUNC=19, POWERFUNC=20, THRESHFUNC=21, UBOUNDFUNC=22, STEPFUNC=23, MAXFUNC=24, MINFUNC=25, ABSFUNC=26, SQRTFUNC=27, SIGFUNC=28;

/// Calculates an equation based on a text string
function equation_calculate(te,t)
{
	te = String(te);
	let eqn = setup_eqn(te);
	
	if(eqn.err == true) return 0;
		
	return calculate(eqn,t);
}


/// Sets up the calculation for the equation
function setup_eqn(te)
{
	te = te.replace(/×/g,"*");                      // Converts × to *

	te = minus_sign_adjust(te)

	let eqn;

	let op = extract_operations(te);
	if(typeof op == "string") return err(op)
	else{
		if(false) print_operations(op);
		
		eqn = create_eqn(op);
		eqn.te = te;
	
		if(false) print_eqn(eqn);
	}
	
	return eqn;
}


/// Takes the text and converts it into a series of operations and objects
function extract_operations(te)
{
	let op = [];
	
	let i = 0;
  while(i < te.length){ 
		let ch = te.substr(i,1);
		
		switch(ch){
			case "\n": break;
      case "(": op.push({type:LEFTBRACKET}); break;
      case ")": op.push({type:RIGHTBRACKET}); break;
      case "|": op.push({type:FUNCDIVIDE}); break;
      case "*": op.push({type:MULTIPLY}); break;
      case "/": op.push({type:DIVIDE}); break;
			case "+": op.push({type:ADD}); break;
      case "-": op.push({type:TAKE}); break;
			case " ": break;

      default:
        let doneflag = false;
				
        if(te.substr(i,4) == "exp(" && doneflag == false){
					op.push({type:EXPFUNC}); 
          doneflag = true;
          i += 2;
        }

        if(te.substr(i,4) == "sin(" && doneflag == false){
					op.push({type:SINFUNC}); 
          doneflag = true;
          i += 2;
        }

        if(te.substr(i,4) == "cos(" && doneflag == false){
					op.push({type:COSFUNC}); 
          doneflag = true;
          i += 2;
        }

        if(te.substr(i,4) == "log(" && doneflag == false){
					op.push({type:LOGFUNC}); 
          doneflag = true;
          i += 2;
        }

        if(te.substr(i,5) == "step(" && doneflag == false){
					op.push({type:STEPFUNC}); 
					doneflag = true;
          i += 3;
        }

				if(te.substr(i,4) == "pow(" && doneflag == false){
					op.push({type:POWERFUNC}); 
          doneflag = true;
          i += 2;
        }
				
				if(te.substr(i,7) == "thresh(" && doneflag == false){
					op.push({type:THRESHFUNC}); 
          doneflag = true;
          i += 5;
        }
				
				if(te.substr(i,7) == "ubound(" && doneflag == false){
					op.push({type:UBOUNDFUNC}); 
          doneflag = true;
          i += 5;
        }
				
				if(te.substr(i,4) == "max(" && doneflag == false){
					op.push({type:MAXFUNC}); 
          doneflag = true;
          i += 2;
        }
				
				if(te.substr(i,4) == "min(" && doneflag == false){
					op.push({type:MINFUNC}); 
          doneflag = true;
          i += 2;
        }
				
				if(te.substr(i,4) == "abs(" && doneflag == false){
					op.push({type:ABSFUNC}); 
          doneflag = true;
          i += 2;
        }
				
				if(te.substr(i,5) == "sqrt(" && doneflag == false){
					op.push({type:SQRTFUNC}); 
          doneflag = true;
          i += 3;
        }
				
				if(te.substr(i,4) == "sig(" && doneflag == false){
					op.push({type:SIGFUNC}); 
          doneflag = true;
          i += 2;
        }
				
        if(doneflag == false){
			    let res = get_float(i,te);
          if(res != undefined){
						i = res.i-1; 
						op.push({type:NUMERIC, constant:res.num}); 
            doneflag = true;
          }
        }

				if(doneflag == false){
					if(te.substr(i,1) == "t" && (i == te.length-1 || notparam_list.includes(te.substr(i+1,1)))){
						op.push({type:TIME}); 
						doneflag = true;
					}
				}
				
				if(doneflag == false){
					return "Problem with equation '"+te+"'. The character '"+ch+"' was not expected.";
				}
        if(doneflag == false){
					return "Problem with equation '"+te+"'. The character '"+ch+"' was not expected.";
				}					
        break;
    }
    i++;
  }

  i = 0;                                           // Adds in unspecified multiply signs
  while(i+1 < op.length){
		if((quant(op,i) || optype(op,i,RIGHTBRACKET)) && 
		    (quant(op,i+1) || is_func(op,i+1) || optype(op,i+1,LEFTBRACKET))){		 
			let item; item.type = MULTIPLY;
      op.splice(i+1,0,item);
    }
    else i++;
  }
	
	return op;
}


/// Prints all the operation for a eqn
function print_operations(op)
{
	error("List of operations:");
	error(op);
	
	let st="";
  for(let i = 0; i < op.length; i++){
		
    switch(op[i].type){
      case LEFTBRACKET: st += "("; break;
      case RIGHTBRACKET: st += ")"; break;
			case FUNCDIVIDE: st += "|"; break;
      case EXPFUNC: st += "exp"; break;
      case SINFUNC: st += "sin"; break;
      case COSFUNC: st += "cos"; break;
      case LOGFUNC: st += "log"; break;
      case STEPFUNC: st += "step"; break;
			case POWERFUNC: st += "pow"; break;
      case THRESHFUNC: st += "thresh"; break; 
      case UBOUNDFUNC: st += "ubound"; break; 
      case MAXFUNC: st += "max"; break; 
			case MINFUNC: st += "min"; break; 
			case ABSFUNC: st += "abs"; break;
			case SQRTFUNC: st += "sqrt"; break;
			case SIGFUNC: st += "sig"; break;
      case ADD: st += "+"; break;
      case TAKE: st += "-"; break;
      case MULTIPLY: st += "*"; break;
			case DIVIDE: st += "/"; break;
      case REG: st += "R"+op[i].num; break;
      case NUMERIC: st += "numeric"+op[i].constant; break;
			case TIME: st += "time"; break;
		}
    st += " ";
  }
	error(st);
}


/// Works out the sequence of eqn to generate result
function create_eqn(op)
{
	let nreg = 0;
	
	let calc=[];
	
	let flag = false;
  do{                                              // Breaks down the equations into simple pieces
    flag = false;
		
    for(let i = 0; i < op.length-1; i++){
      switch(op[i].type){
			case LEFTBRACKET: case RIGHTBRACKET:
				break;
				
			case ADD: case TAKE: case MULTIPLY: case DIVIDE:
				if(i > 0 && quant(op,i-1) == true && quant(op,i+1) == true && prio(op,i) >= prio(op,i-2) && prio(op,i) >= prio(op,i+2)){
					let add_flag = false;
					if(calc.length > 0){
						let cal = calc[calc.length-1];
						if(cal.op == op[i].type && (cal.op == ADD || cal.op == MULTIPLY)){
							if(op[i-1].type == REG && op[i-1].num == cal.reg_store){
								cal.item.push(copy(op[i+1]));
								add_flag = true;
							} 
						}
					}

					if(add_flag == true){                    // Combines with previous operation
						op.splice(i,2);
						i--;
					}
					else{                                    // Inserts a new eqn
						let cal = {op:op[i].type, reg_store:nreg, item:[]};
						cal.item.push(copy(op[i-1]));
						cal.item.push(copy(op[i+1]));
					
						calc.push(cal);
					
						i--;
						op[i].type = REG;
						op[i].num = nreg;
						op[i].constant = undefined;
				
						op.splice(i+1,2);
						nreg++;
					}
					flag = true;
				}
				break;
				
			case REG: case NUMERIC: case TIME: 
				if(optype(op,i-1,LEFTBRACKET) && optype(op,i+1,RIGHTBRACKET) && !is_func(op,i-2)){
					i--;
					op[i] = op[i+1];
					op.splice(i+1,2);
			 
					flag = true;
				}
				break;
				
			case FUNCDIVIDE: break;
			
			case EXPFUNC: case LOGFUNC: case SINFUNC: case COSFUNC: case STEPFUNC: case ABSFUNC: 
			case SQRTFUNC: case SIGFUNC: // Univariate
				if(optype(op,i+1,LEFTBRACKET) && quant(op,i+2) == true && optype(op,i+3,RIGHTBRACKET)){		
					let cal = {op:op[i].type, reg_store:nreg, item:[]} ;
					cal.item.push(copy(op[i+2]));
					calc.push(cal);
					
					op[i].type = REG;
					op[i].num = nreg;
				
					op.splice(i+1,3);
					nreg++;
					
					flag = true;
				}
				break;
				
			// Functions with two variables
			case POWERFUNC: case THRESHFUNC: case UBOUNDFUNC: case MAXFUNC: case MINFUNC:  
				if(optype(op,i+1,LEFTBRACKET) && quant(op,i+2) == true && optype(op,i+3,FUNCDIVIDE) && quant(op,i+4) == true && optype(op,i+5,RIGHTBRACKET)){		
					let cal = {op:op[i].type, reg_store:nreg, item:[]};
					cal.item.push(copy(op[i+2]));
					cal.item.push(copy(op[i+4]));
					calc.push(cal);
					
					op[i].type = REG;
					op[i].num = nreg;
					
					op.splice(i+1,5);
					nreg++;
					
					flag = true;
				}
				break;
      }
    }
  }while(flag == true && op.length > 1);

	if(op.length != 1) return err("Problem with expression. Calculation cannot be performed.");
	
	let ans = op[0];
	
	return {nreg:nreg, ans:ans, calc:calc}
}


/// Prints steps used for a eqn
function print_eqn(eqn)  
{
	error("For equation '"+eqn.te+"' eqn:");
	if(eqn.err == true){ error("Could not evaluate equation: "+eqn.msg); return;}

	let calc = eqn.calc;
  for(let i = 0; i < calc.length; i++){
		let st = "";
		let ca = calc[i];
	
		switch(ca.op){
		case EXPFUNC: st += "exp("; break;
		case SINFUNC: st += "sin("; break;
		case COSFUNC: st += "cos("; break;
		case LOGFUNC: st += "log("; break;
		case STEPFUNC: st += "step("; break;
		case POWERFUNC: st += "power("; break;
		case THRESHFUNC: st += "thresh("; break;
		case UBOUNDFUNC: st += "ubound("; break;
		case MAXFUNC: st += "max("; break; 
		case MINFUNC: st += "min("; break; 
		case ABSFUNC: st += "abs("; break;
		case SQRTFUNC: st += "sqrt("; break;
		case SIGFUNC: st += "sig("; break;
		case ADD: break;
		case TAKE: break;
		case MULTIPLY: break;
		case DIVIDE: break;
		default: alert_help("Equation error","Code 1"); break;
		}
			
		for(let j = 0; j < ca.item.length; j++){
			let it = ca.item[j];
			switch(it.type){
			case REG: st += "R"+it.num; break;
			case NUMERIC: st += it.constant; break;
			case TIME: st += "time"; break;
			default: error("Eq problem2"); break;
			}
			
			if(j != ca.item.length -1){
				switch(ca.op){
				case ADD: st += "+"; break;
				case TAKE: st += "-"; break;
				case MULTIPLY: st += "*"; break;
				case DIVIDE: st += "/"; break;
				case POWERFUNC: st += "|"; break;
				case THRESHFUNC: st += "|"; break;
				case UBOUNDFUNC: st += "|"; break;
				default: break;
				}
			}
		}
	
		switch(ca.op){
		case EXPFUNC: st += ")"; break;
		case SINFUNC: st += ")"; break;
		case COSFUNC: st += ")"; break;
		case LOGFUNC: st += ")"; break;
		case STEPFUNC: st += ")"; break;
		case POWERFUNC: st += ")"; break;
		case THRESHFUNC: st += ")"; break;
		case UBOUNDFUNC: st += ")"; break;
		default: break;
		}
			
    st += " > ";

		error(st+"R"+ca.reg_store);
  }

	let an="";
  switch(ans.type){
    case REG: an += "R"+ans.num; break;
		case NUMERIC: an += ans.constant; break;
		case TIME: an += "time"; break;
		default: error("Eq problem3: "+ans.type); break;
  }
  error(an+" Answer");
	error("");
}


/// Tries to get a number from a string
function get_float(i,te)
{
	let sign = 1;
	if(te.substr(i,1) == "-"){ sign = -1; i++;}

  let num = te.charCodeAt(i)-48; if(!(num >= 0 && num <= 9)) return;

  let num2 = 0.0;
  do{
    num = te.charCodeAt(i)-48;
    if(num >= 0 && num <= 9) num2 = num2*10+num;
    else break;
    i++;
  }while(i < te.length);

  if(te.substr(i,1) == "."){                       // Looks for a decimal point
    i++;
    let fac = 0.1;
    do{
			num = te.charCodeAt(i)-48;
      if(num >= 0 && num <= 9){ num2 += fac*num; fac /= 10;}
      else break;
      i++;
    }while(i < te.length);
  }

	return {num:sign*num2, i:i};
}


/// Determines if a quantity exists
function quant(op,i)
{
  if(i < 0 || i >= op.length) return false;
  switch(op[i].type){
	case REG: case NUMERIC: case TIME: return true;
  default: return false;
	}
}


/// Determines if a operator is a particular type
function optype(op,i,type)
{
	if(i < 0 || i >= op.length) return false;
  if(op[i].type == type) return true;
	return false;
}


/// Determines if a function exists
function is_func(op,i)
{
  if(i < 0 || i >= op.length) return false;
  switch(op[i].type){
	case EXPFUNC: case SINFUNC: case COSFUNC: case LOGFUNC: case POWERFUNC: 
	case THRESHFUNC: case UBOUNDFUNC: 
	case STEPFUNC: case MAXFUNC: case MINFUNC: case ABSFUNC: case SQRTFUNC: case SIGFUNC:
		return true;
	default: return false;
	}
}


/// Gives the relative priority of an operator
function prio(op,i)
{
  if(i < 0 || i >= op.length) return -1;
  switch(op[i].type){
    case ADD: return 0;
    case TAKE: return 0;
    case MULTIPLY: return 1;
    case DIVIDE: return 1;
    default: return -1;
  }
}


/// Copies with minus signs which don't have a number before them
function minus_sign_adjust(te)
{
  let i = 0;                                             
  while(i < te.length){
    if(te.substr(i,1) == "-"){
      if(i == 0 && te.substr(i,1) == "-") te = "0"+te;
      else{
        if(te.substr(i-1,1) == "(" || te.substr(i-1,1) == "|") te = te.substr(0,i)+"0"+te.substr(i);
        else i++;
      }
    }
    else i++;
  }

	return te;
}


/// Performs an operation of a set of numbers 
function calculate_operation(op,num)
{
	let N = num.length;
	
	switch(op){
	case ADD:
		{
			let ans = 0.0; 
			for(let i = 0; i < N; i++) ans += num[i];
			return ans;
		}
		
	case TAKE:
		if(N != 2) error("For TAKE should be 2");
		return num[0]-num[1]; 
	
	case MULTIPLY:
		{
			let ans = 1.0; 
			for(let i = 0; i < N; i++) ans *= num[i];
			return ans;
		}
		
	case DIVIDE:
		if(N != 2) error("For DIVIDE should be 2");
		if(num[1] == 0) error("Equation '"+te_raw+"' caused a division by zero."); 
		return num[0]/num[1]; 
		
	case EXPFUNC:
		if(N != 1) error("For EXPFUNC should be 1");
		return Math.exp(num[0]); 
			
	case SINFUNC:
		if(N != 1) error("For SINFUNC should be 1");
		return Math.sin(num[0]); 
	
	case COSFUNC: 
		if(N != 1) error("For COSFUNC should be 1");
		return Math.cos(num[0]); 
		
	case LOGFUNC: 
		if(N != 1) error("For LOGFUNC should be 1");
		if(num[0] <= 0) error("For equation '"+te_raw+"' the quantity inside the log became negative."); 
		return Math.log(num[0]);
		
	case STEPFUNC:
		if(N != 1) error("For STEPFUNC should be 1");
		if(num[0] > 0) return 1;
		return 0; 
		
	case POWERFUNC:
		if(N != 2) error("For POWERFUNC should be 2");
		return Math.pow(num[0],num[1]);
	
	case THRESHFUNC:
		if(N != 2) error("For THRESHFUNC should be 2");
		if(num[0] < num[1]) return 0;
		return num[0];
		
	case UBOUNDFUNC:
		if(N != 2) error("For UBOUNDFUNC should be 2");
		if(num[0] > num[1]) return Infinity;
		return num[0];
		
	case MAXFUNC:
		if(N != 2) error("For MAXFUNC should be 2");
		if(num[0] > num[1]) return num[0];
		return num[1]; 
		
	case MINFUNC:
		if(N != 2) error("For MINFUNC should be 2");
		if(num[0] < num[1]) return num[0];
		return num[1]; 
		
	case ABSFUNC:
		if(N != 1) error("For MINFUNC should be 1");
		if(num[0] > 0) return num[0]; 
		return -num[0];
		
	case SQRTFUNC:
		if(N != 1) error("For SQRTFUNC should be 1");
		if(num[0] < 0) alert_help("Equation error","Cannot have a negative quantity in a sqaure root");
		return Math.sqrt(num[0]);
		
	case SIGFUNC:
		if(N != 1) error("For SIGFUNC should be 1");
		return 1/(1+Math.exp(-num[0]));
		
	default: alert_help("Equation error","Code 2"); break;
	}
	
	return 0;
}


/// Calculates the value for an equation
function calculate(eqn,t) 
{
 	let regcalc=[];

	let calc = eqn.calc;
  for(let i = 0; i < calc.length; i++){
		let ca = calc[i];
		
		let item = ca.item;
		let N = item.length;
		
		let num=[];
		
		for(let j = 0; j < N; j++){
			let it = item[j];
			
			switch(it.type){
				case REG: num[j] = regcalc[it.num]; break;
				case NUMERIC: num[j] = it.constant; break;
				case TIME: num[j] = t; break;
				default: alert_help("Equation error","Code3"); break;
			}
		}

	  regcalc[ca.reg_store] = calculate_operation(ca.op,num);
	}

	let ans	= eqn.ans;
  switch(ans.type){
    case REG: return regcalc[ans.num];
    case NUMERIC: return ans.constant;
		case TIME: return t;
		default: alert_help("Equation error","Code 4"); return UNSET;
  }
}

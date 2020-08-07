double ratecalc(long eq, double *popnum, double *param)   // Calculates the rate of an event
{
  double r;

  r = calculate(eq,popnum,param);
  if(r < ratesmall) r = ratesmall;
  return r;
}

double ratecalcdep(long d, double *popnum, double *param)    // uses ncompind and popnum from the reduced set for a link
{
  long i, eq;
  double num, num1, num2, r;

  eq = transdepeq[d];
  for(i = 0; i < neqncalc[eq]; i++){
    switch(eqncalc1[eq][i]){
    case PARAM: num1 = param[eqncalcnum1[eq][i]]; break;
    case POPNUM: num1 = popnum[eqncalcfastnum1[d][i]]; break;
    case REG: num1 = regcalc[eqncalcnum1[eq][i]]; break;
    case NUMERIC: num1 = numeric[eqncalcnum1[eq][i]]; break;
    }

    switch(eqncalc2[eq][i]){
    case PARAM: num2 = param[eqncalcnum2[eq][i]]; break;
    case POPNUM: num2 = popnum[eqncalcfastnum2[d][i]]; break;
    case REG: num2 = regcalc[eqncalcnum2[eq][i]]; break;
    case NUMERIC: num2 = numeric[eqncalcnum2[eq][i]]; break;
    }

    switch(eqncalcop[eq][i]){
    case ADD: num = num1+num2; break;
    case TAKE: num = num1-num2; break;
    case MULTIPLY: num = num1*num2; break;
    case DIVIDE: if(num2 == 0) emsg("A division by zero!"); num = num1/num2; break;
    case POWER: num = pow(num1,num2); break;
    case EXPFUNC: num = exp(num2); break;
    case SINFUNC: num = sin(num2); break;
    case COSFUNC: num = cos(num2); break;
    case LOGFUNC: if(num2 <= 0) emsg("Log cannot be calculated"); num = log(num2); break;
    case STEPFUNC: if(num2 > 0) num = 1; else num = 0; break;
    }

    regcalc[eqncalcstore[eq][i]] = num;
  }

  switch(eqncalcans[eq]){
  case PARAM: r = param[eqncalcansnum[eq]]; break;
  case POPNUM: r = popnum[eqncalcansfastnum[d]]; break;
  case REG: r = regcalc[eqncalcansnum[eq]]; break;
  case NUMERIC: r = numeric[eqncalcansnum[eq]]; break;
  }

  if(r < ratesmall) r = ratesmall;

  return r;
}

double ratecalcnotdep(long eq, double *paramval)    // uses ncompind and popnum from the reduced set for a link
{
  long i;
  double num, num1, num2, r;

  for(i = 0; i < neqncalc[eq]; i++){
    switch(eqncalc1[eq][i]){
    case PARAM: num1 = paramval[eqncalcnum1[eq][i]]; break;
    case POPNUM: emsg("not dep prob"); break;
    case REG: num1 = regcalc[eqncalcnum1[eq][i]]; break;
    case NUMERIC: num1 = numeric[eqncalcnum1[eq][i]]; break;
    }

    switch(eqncalc2[eq][i]){
    case PARAM: num2 = paramval[eqncalcnum2[eq][i]]; break;
    case POPNUM: emsg("not dep prob"); break;
    case REG: num2 = regcalc[eqncalcnum2[eq][i]]; break;
    case NUMERIC: num2 = numeric[eqncalcnum2[eq][i]]; break;
    }

    switch(eqncalcop[eq][i]){
    case ADD: num = num1+num2; break;
    case TAKE: num = num1-num2; break;
    case MULTIPLY: num = num1*num2; break;
    case DIVIDE: if(num2 == 0) emsg("A division by zero!"); num = num1/num2; break;
    case POWER: num = pow(num1,num2); break;
    case EXPFUNC: num = exp(num2); break;
    case SINFUNC: num = sin(num2); break;
    case COSFUNC: num = cos(num2); break;
    case LOGFUNC: if(num2 <= 0) emsg("Log cannot be calculated"); num = log(num2); break;
    case STEPFUNC: if(num2 > 0) num = 1; else num = 0; break;
    }

    regcalc[eqncalcstore[eq][i]] = num;
  }

  switch(eqncalcans[eq]){
  case PARAM: r = paramval[eqncalcansnum[eq]]; break;
  case POPNUM: emsg("not dep prob"); break;
  case REG: r = regcalc[eqncalcansnum[eq]]; break;
  case NUMERIC: r = numeric[eqncalcansnum[eq]]; break;
  }

  if(r < ratesmall) r = ratesmall;

  return r;
}

/*
double ratecalcdeptakeoff(long d, double *popnum, double *param, long c)    // uses ncompind and popnum from the reduced set for a link
{
  long i, eq, p;
  double num, num1, num2, r;

  eq = transdepeq[d];

  for(i = 0; i < neqncalc[eq]; i++){
    switch(eqncalc1[eq][i]){
    case PARAM: num1 = param[eqncalcnum1[eq][i]]; break;
    case POPNUM:
      p = eqncalcfastnum1[d][i];
      num1 = popnum[p]; if(c != -1) num1 += popnumtake[eq_popnum[eq][p]][c];
      break;
    case REG: num1 = regcalc[eqncalcnum1[eq][i]]; break;
    case NUMERIC: num1 = numeric[eqncalcnum1[eq][i]]; break;
    default:
      switch(eqncalcop[eq][i]){
      case ADD: case TAKE: case MULTIPLY: case DIVIDE: case POWER: emsg(eq+"probratecalca"); break;
      }
      break;
    }

    switch(eqncalc2[eq][i]){
    case PARAM: num2 = param[eqncalcnum2[eq][i]]; break;
    case POPNUM: p = eqncalcfastnum2[d][i]; num2 = popnum[p] + popnumtake[eq_popnum[eq][p]][c]; break;
    case REG: num2 = regcalc[eqncalcnum2[eq][i]]; break;
    case NUMERIC: num2 = numeric[eqncalcnum2[eq][i]]; break;
    default: emsg("probratecalcb"); break;
    }

    switch(eqncalcop[eq][i]){
    case ADD: num = num1+num2; break;
    case TAKE: num = num1-num2; break;
    case MULTIPLY: num = num1*num2; break;
    case DIVIDE: if(num2 == 0) emsg("A division by zero!"); num = num1/num2; break;
    case POWER: num = pow(num1,num2); break;
    case EXPFUNC: num = exp(num2); break;
    case SINFUNC: num = sin(num2); break;
    case COSFUNC: num = cos(num2); break;
    case LOGFUNC: if(num2 <= 0) emsg("Log cannot be calculated"); num = log(num2); break;
    case STEPFUNC: if(num2 > 0) num = 1; else num = 0; break;
    default: emsg("probratecalc1"); break;
    }

    regcalc[eqncalcstore[eq][i]] = num;
  }

  switch(eqncalcans[eq]){
  case PARAM: r = param[eqncalcansnum[eq]]; break;
  case POPNUM: p = eqncalcansfastnum[d]; r = popnum[p] + popnumtake[eq_popnum[eq][p]][c]; break;
  case REG: r = regcalc[eqncalcansnum[eq]]; break;
  case NUMERIC: r = numeric[eqncalcansnum[eq]]; break;
  default: emsg("probratecalc2"); break;
  }

  if(r < ratesmall) r = ratesmall;

  return r;
}
*/

double ratecalcdeptakeoff(long d, double *popnum, double *param, long c)    // Removes the individual at c when the quantity is calculated
{
  long i, eq, p;
  double num, num1, num2, r;

  eq = transdepeq[d];

  for(i = 0; i < neqncalc[eq]; i++){
    switch(eqncalc1[eq][i]){
      case PARAM: num1 = param[eqncalcnum1[eq][i]]; break;
      case POPNUM:
        p = eqncalcfastnum1[d][i];
        num1 = popnum[p]; if(c != -1) num1 += popnumtake[eq_popnum[eq][p]][c];
        break;
        case REG: num1 = regcalc[eqncalcnum1[eq][i]]; break;
        case NUMERIC: num1 = numeric[eqncalcnum1[eq][i]]; break;
      default:
        switch(eqncalcop[eq][i]){
          case ADD: case TAKE: case MULTIPLY: case DIVIDE: case POWER: emsg(eq+"probratecalca"); break;
        }
        break;
    }

    switch(eqncalc2[eq][i]){
      case PARAM: num2 = param[eqncalcnum2[eq][i]]; break;
      case POPNUM: p = eqncalcfastnum2[d][i]; num2 = popnum[p] + popnumtake[eq_popnum[eq][p]][c]; break;
      case REG: num2 = regcalc[eqncalcnum2[eq][i]]; break;
      case NUMERIC: num2 = numeric[eqncalcnum2[eq][i]]; break;
      default: emsg("probratecalcb"); break;
    }

    switch(eqncalcop[eq][i]){
      case ADD: num = num1+num2; break;
      case TAKE: num = num1-num2; break;
      case MULTIPLY: num = num1*num2; break;
      case DIVIDE: if(num2 == 0) emsg("A division by zero!"); num = num1/num2; break;
      case POWER: num = pow(num1,num2); break;
      case EXPFUNC: num = exp(num2); break;
      case SINFUNC: num = sin(num2); break;
      case COSFUNC: num = cos(num2); break;
      case LOGFUNC: if(num2 <= 0) emsg("Log cannot be calculated"); num = log(num2); break;
      case STEPFUNC: if(num2 > 0) num = 1; else num = 0; break;
      default: emsg("probratecalc1"); break;
    }

    regcalc[eqncalcstore[eq][i]] = num;
  }

  switch(eqncalcans[eq]){
    case PARAM: r = param[eqncalcansnum[eq]]; break;
    case POPNUM: p = eqncalcansfastnum[d]; r = popnum[p] + popnumtake[eq_popnum[eq][p]][c]; break;
    case REG: r = regcalc[eqncalcansnum[eq]]; break;
    case NUMERIC: r = numeric[eqncalcansnum[eq]]; break;
    default: emsg("probratecalc2"); break;
  }

  if(r < ratesmall) r = ratesmall;

  return r;
}

double Chain::transra(long eq, double t, long cnow)                      // Calculates the transition rate but removes the individual at cnow
{
  if(transdep[eq] == 1) return transradep(eq,t,cnow);
  else return transnotdepeq_val[transdepref[eq]];
}

double Chain::transradep(long eq, double t, long cnow)                      // Calculates the transition rate for dependent case but removes the individual at cnow
{
  long d;
  long e, ee;

  //t += tiny;

  if(eq == -1 || transdep[eq] != 1) emsg("should be dep");                          // calculates the value of the link rate

  d = transdepref[eq];
  t *= onefac;

  e = depeqdiv[d][long(t*depeqdivmax/tmax)];
  do{
    ee = depeq_evnext[e]; if(ee == -1) emsg("neg ee");
    if(ee == -1 || depeq_evt[ee] > t) break;
    e = ee;
  }while(1 == 1);

  /*
  if(checkon == 1){ // used for checking
    long c=NOTALIVE,k = 0; while(k < nindev[ist] && indev[ist][k].t < t){ c = tra[indev[ist][k].tr].cf; k++;}
    if(c != cnow) emsg("cnow prob2");
  }
  */

  if(transdepeqrecalc[d][cnow] == 0) return depeq_evval[e];
  else{ return ratecalcdeptakeoff(d,depeq_evpopnum[e],param,cnow);}
}

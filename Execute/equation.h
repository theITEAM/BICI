// Parses an equation string and 

const long LEFTBRACKET = 0, RIGHTBRACKET = 1, ADD = 4, TAKE = 5, MULTIPLY = 6;                 // Diferent elements to an equation
const long DIVIDE = 3, REG = 8, NUMERIC = 9, PARAM = 10, POPNUM = 11;
const long EXPFUNC = 20, SINFUNC = 21, COSFUNC = 22, LOGFUNC = 23, POWER = 24, STEPFUNC = 25;
const double NULLVAL = -1.0001;

vector<string> eqnstr;                              // Stores the string which describes an equation

vector<long> neqncalc;                              // The number of basic operations to calculate an equation
vector< vector<long> > eqncalc1;                    // Stores equation calculations
vector< vector<long> > eqncalcnum1;
vector< vector<long> > eqncalc2;
vector< vector<long> > eqncalcnum2;
vector< vector<long> > eqncalcstore;
vector< vector<long> > eqncalcop;
vector<long> eqncalcans;
vector<long> eqncalcansnum;

vector< vector<long> > eqncalcfastnum1;             // Look up table used so that the eq timelines can easily calculate population number
vector< vector<long> > eqncalcfastnum2;
vector<long> eqncalcansfastnum;

long raend;                                         // Temporary global variables
vector<long> op;
vector<long> opnum;
vector<double> regcalc;
vector<double> numeric;
long popcomb;
double popcombnum;

//void printmodel();
void printop();                                     // Outputs the operations for an equation calculation
double getfloat(string text, long i);               // Tries to get a number from a string
long getparamname(string text, long i);             // Tries to get a parameter from a string
long getpopnumname(string text, long i);            // Tries to get a population name from a string
long quant(long i);                                 // Determines if a quantitiy exists
long prio(long i);                                  // Gives the relative priority of an operator
void printcalc(long eq);                            // Outputs the calculation
long popcombine(long i);                            // Works out if a particular population can be combined to others

long addequation(string s)                          // Adds a new equatrion to the list
{
  long len, i, j, jst, jst2, doneflag, numi, flag, eq, st, p, ist, npopnum, plfl = 0, c;
  double num;
  string a, name;
  vector<long> state, vec;
  vector<float> vecf;
  vector<double> reg;
  vector<long> calc1;
  vector<long> calcnum1;
  vector<long> calc2;
  vector<long> calcnum2;
  vector<long> calcstore;
  vector<long> calcop;
  vector <double> w;

  if(s == "") return -1;

  i = 0;                                          // Cope with a negate minus sign
  while(i < s.length()){
    if(s.substr(i,1) == "-"){
      if(i == 0 && s.substr(i,1) == "-") s = "0"+s;
      else{
        if(s.substr(i-1,1) == "(") s = s.substr(0,i)+"0"+s.substr(i);
        else i++;
      }
    }
    else i++;
  }

  op.clear(); opnum.clear();

  for(eq = 0; eq < eqnstr.size(); eq++){
    if(eqnstr[eq] == s) return eq;
  }

  j = 0;                                             // Extracts popnum
  while(j < long(s.length())-1){
    if(s.substr(j,1) == "{"){
      j += 1;
      jst = j;
      while(j < s.length() && s.substr(j,1) != "}") j++;
      if(j == s.length()) emsg("popnumprob");

      state = getallcomp(s.substr(jst,j-jst));

      vec.clear(); vecf.clear();

      for(st = 0; st < state.size(); st++){
        vec.push_back(state[st]);
        vecf.push_back(1);
      }

      name = s.substr(jst-1,j-(jst-1)+1);
      npopnum = popnumterm.size();
      for(p = 0; p < npopnum; p++) if(popnumname[p] == name) break;
      if(p == npopnum){
        popnumname.push_back(name);
        popnumterm.push_back(vec);
        popnumtermweight.push_back(vecf);
      }
    }
    j++;
  }

  //for(i = 0; i < popnumname.size(); i++) cout << popnumname[i] << " popnum\n";

  len = s.length();

  i = 0;
  while(i < len){ 
    const char *cstr = s.substr(i,1).c_str();
    switch(cstr[0]){
      case '(': op.push_back(LEFTBRACKET); opnum.push_back(-1); break;
      case ')': op.push_back(RIGHTBRACKET); opnum.push_back(-1); break;
      case '*':  op.push_back(MULTIPLY); opnum.push_back(-1); break;
      case '/': op.push_back(DIVIDE); opnum.push_back(-1); break;
      case '^': op.push_back(POWER); opnum.push_back(-1); break;
      case '+': op.push_back(ADD); opnum.push_back(-1); break;
      case '-': op.push_back(TAKE); opnum.push_back(-1); break;
      default:
        doneflag = 0;
        if(s.substr(i,4) == "exp(" && doneflag == 0){
          op.push_back(EXPFUNC); opnum.push_back(-1);
          doneflag = 1;
          i += 2;
        }

        if(s.substr(i,4) == "sin(" && doneflag == 0){
          op.push_back(SINFUNC); opnum.push_back(-1);
          doneflag = 1;
          i += 2;
        }

        if(s.substr(i,4) == "cos(" && doneflag == 0){
          op.push_back(COSFUNC); opnum.push_back(-1);
          doneflag = 1;
          i += 2;
        }

        if(s.substr(i,4) == "log(" && doneflag == 0){
          op.push_back(LOGFUNC); opnum.push_back(-1);
          doneflag = 1;
          i += 2;
        }

        if(s.substr(i,5) == "step(" && doneflag == 0){
          op.push_back(STEPFUNC); opnum.push_back(-1);
          doneflag = 1;
          i += 3;
        }

        if(doneflag == 0){
          num = getfloat(s,i);
          if(num != NULLVAL){
            i = raend-1; 
            op.push_back(NUMERIC); opnum.push_back(numeric.size());
            numeric.push_back(num);
            doneflag = 1;
          }
        }

        if(doneflag == 0){
          numi = getparamname(s,i);
          if(numi >= 0){
            i = raend-1;
            op.push_back(PARAM); opnum.push_back(numi);
            doneflag = 1;
          }
        }

        if(doneflag == 0){
          numi = getpopnumname(s,i);
          if(numi >= 0){
            i = raend-1;
            op.push_back(POPNUM); opnum.push_back(numi);
            doneflag = 1;
          }
        }

        if(doneflag == 0){   // reads in a new parameter
          a = s.substr(i,1);
          if(a == "["){
            ist = i;
            do{
              a = s.substr(i,1);
              i++;
            }while(i < len && a != "]");
            if(a != "]") emsg("eqnprob");

            op.push_back(PARAM); opnum.push_back(paramname.size());

            paramname.push_back(s.substr(ist+1,i-ist-2));
            paramprior.push_back(-1);
            paramclass.push_back("");
            parampriordep.push_back(vector<long>());
            doneflag = 1;
            i--;
          }
        }

        if(doneflag == 0) emsg(s.substr(i)+"Problem with expression1");
        break;
    }
    i++;
  }

  i = 0;                                            // Adds in unspecified multiply signs
  while(i < long(op.size())-1){
    if((op[i] == PARAM || op[i] == POPNUM || op[i] == RIGHTBRACKET) &&
        (op[i+1] == PARAM || op[i+1] == POPNUM || op[i+1] == EXPFUNC || op[i+1] == SINFUNC || op[i+1] == COSFUNC || op[i+1] == LOGFUNC || 
         op[i+1] == STEPFUNC || op[i+1] == LEFTBRACKET)){
      op.insert(op.begin()+i+1,MULTIPLY);
      opnum.insert(opnum.begin()+i+1,-1);
    }
    else i++;
  }

  if(plfl == 1) printop();

  i = 0; w.resize(ncomp);                               // Looks to combine together multiple population expressions
  while(i < long(op.size())-3){
    if(popcombine(i) == 1){
      ist = i; for(c = 0; c < ncomp; c++) w[c] = 0;
      stringstream ss;
      while(popcombine(i) == 1){
        ss << popnumname[popcomb] << "*" << popcombnum;

        for(j = 0; j < popnumterm[popcomb].size(); j++){
          w[popnumterm[popcomb][j]] += popcombnum*popnumtermweight[popcomb][j];
        }
        i += 3;

        if(i == op.size()) break;
        if(op[i] != ADD && op[i] != TAKE) break;
        i++;
        ss << "+";
      }

      vec.clear(); vecf.clear();
      for(c = 0; c < ncomp; c++){
        if(w[c] != 0){
          vec.push_back(c);
          vecf.push_back(w[c]);
        }
      }

      name = ss.str();
      npopnum = popnumterm.size();
      for(p = 0; p < npopnum; p++){ if(popnumname[p] == name) break;}
      if(p == npopnum){
        popnumname.push_back(ss.str());
        popnumterm.push_back(vec);
        popnumtermweight.push_back(vecf);
      }

      op[ist] = POPNUM; opnum[ist] = p;
      op.erase(op.begin()+ist+1,op.begin()+i);
      opnum.erase(opnum.begin()+ist+1,opnum.begin()+i);
      i = ist;
    }
    i++;
  }

  do{                                               // Breaks down the equations into simple pieces
    flag = 0;
    for(i = 0; i < long(op.size())-1; i++){
      switch(op[i]){
        case ADD: case TAKE: case MULTIPLY: case DIVIDE: case POWER:
          if(i > 0 && quant(i-1) == 1 && quant(i+1) == 1 && prio(i) >= prio(i-2) && prio(i) >= prio(i+2)){
            calc1.push_back(op[i-1]); calcnum1.push_back(opnum[i-1]);
            calc2.push_back(op[i+1]); calcnum2.push_back(opnum[i+1]);
            calcstore.push_back(reg.size());
            calcop.push_back(op[i]);

            i--;
            op[i] = REG;
            opnum[i] = reg.size();

            op.erase(op.begin()+i+1,op.begin()+i+3);
            opnum.erase(opnum.begin()+i+1,opnum.begin()+i+3);
            reg.push_back(0);
            flag = 1;
          }
          break;
        case PARAM: case POPNUM: case REG: case NUMERIC:
          if(i > 0){
            if(op[i-1] == LEFTBRACKET && op[i+1] == RIGHTBRACKET){
              i--;
              op[i] = op[i+1];
              opnum[i] = opnum[i+1];

              op.erase(op.begin()+i+1,op.begin()+i+3);
              opnum.erase(opnum.begin()+i+1,opnum.begin()+i+3);

              i--;                                      // Looks for a fuction
              if(i >= 0){
                switch(op[i]){
                  case EXPFUNC: case SINFUNC: case COSFUNC: case LOGFUNC: case STEPFUNC:
                    if(op[i+1] == PARAM || op[i+1] == POPNUM || op[i+1] == REG || op[i+1] == NUMERIC){
                      calc1.push_back(-1); calcnum1.push_back(-1);
                      calc2.push_back(op[i+1]); calcnum2.push_back(opnum[i+1]);
                      calcstore.push_back(reg.size());
                      calcop.push_back(op[i]);

                      op[i] = REG;
                      opnum[i] = reg.size();

                      op.erase(op.begin()+i+1);
                      opnum.erase(opnum.begin()+i+1);
                      reg.push_back(0);
                    }
                    break;
                }
              }
              flag = 1;
            }
          }
          break;
      }
    }

    if(plfl == 1) printop();
  }while(flag == 1 && op.size() > 1);

  if(op.size() != 1){ cout << s << " s\n"; emsg("Problem with expression2");}

  neqncalc.push_back(calc1.size());
  eqncalc1.push_back(calc1);
  eqncalcnum1.push_back(calcnum1);
  eqncalc2.push_back(calc2);
  eqncalcnum2.push_back(calcnum2);
  eqncalcstore.push_back(calcstore);
  eqncalcop.push_back(calcop);
  eqncalcans.push_back(op[0]);
  eqncalcansnum.push_back(opnum[0]);

  eqnstr.push_back(s);

  if(reg.size() > regcalc.size()) regcalc.resize(reg.size());

  //printcalc(eqncalcans.size()-1);

  return long(eqncalcans.size())-1;
}

long popcombine(long i)   // works out if a particular population can be combined to others
{
  if(i > long(op.size())-3) return 0;

  popcomb = -1;

  if(op[i] == POPNUM && op[i+1] == MULTIPLY && op[i+2] == NUMERIC){ popcomb = opnum[i]; popcombnum = numeric[opnum[i+2]];}

  if(op[i] == NUMERIC && op[i+1] == MULTIPLY && op[i+2] == POPNUM){ popcomb = opnum[i+2]; popcombnum = numeric[opnum[i]];}

  if(op[i] == POPNUM && op[i+1] == DIVIDE && op[i+2] == NUMERIC){ popcomb = opnum[i]; popcombnum = 1.0/numeric[opnum[i+2]];}

  if(popcomb == -1) return 0;
  if(i > 0){
    switch(op[i-1]){
      case MULTIPLY: case DIVIDE: case POWER: return 0;
      case TAKE: popcombnum *= -1; break;
    }
  }

  if(i+3 < op.size()){
    switch(op[i+3]){
      case MULTIPLY: case DIVIDE: case POWER: return 0;
    }
  }

  return 1;
}

long quant(long i)                           // Determines if a quantitiy exists
{
  if(i < 0 || i >= long(op.size())) return 0;
  if(op[i] == PARAM || op[i] == POPNUM || op[i] == REG || op[i] == NUMERIC) return 1;
  else return 0;
}

long prio(long i)                            // Gives the relative priority of an operator
{
  if(i < 0 || i >= op.size()) return -1;
  switch(op[i]){
    case ADD: return 0;
    case TAKE: return 0;
    case MULTIPLY: return 1;
    case DIVIDE: return 1;
    case POWER: return 2;
    default: return -1;
  }
}

void printop()                                // Outputs the operations for an equation calculation
{
  long i;

  for(i = 0; i < op.size(); i++){
    switch(op[i]){
      case LEFTBRACKET: cout << "("; break;
      case RIGHTBRACKET: cout << ")"; break;
      case PARAM: cout << paramname[opnum[i]]; break;
      case POPNUM: cout << popnumname[opnum[i]]; break;
      case EXPFUNC: cout << "exp"; break;
      case SINFUNC: cout << "sin"; break;
      case COSFUNC: cout << "cos"; break;
      case LOGFUNC: cout << "log"; break;
      case STEPFUNC: cout << "step"; break;
      case ADD: cout << "+"; break;
      case TAKE: cout << "-"; break;
      case MULTIPLY: cout << "*"; break;
      case POWER: cout << "^"; break;
      case DIVIDE: cout << "/"; break;
      case REG: cout << "R" << opnum[i]; break;
      case NUMERIC: cout << numeric[opnum[i]]; break;
    }
    cout << " ";
  }
  cout << "\n";
}

long getparamname(string text, long i)       // Tries to get a parameter from a string
{
  long p, l, num, maxlen, pmax;

  maxlen = 0; pmax = -1;
  for(p = 0; p < paramname.size(); p++){
    l = paramname[p].length();
    if(text.substr(i,l+2) == "["+paramname[p]+"]"){
      pmax = p; raend = i+l+2; break;
    }
  }

  return pmax;
}

long getpopnumname(string text, long i)                    // Tries to get a population name from a string
{
  long j, l, num, maxlen, jmax;

  maxlen = 0; jmax = -1;
  for(j = 0; j < popnumname.size(); j++){
    l = popnumname[j].length();
    if(text.substr(i,l) == popnumname[j]){
      const char *cstr = text.substr(i+l,1).c_str();
      num = cstr[0];
      if(i+l == text.length() || (!(num >= 48 && num < 58) && !(num >= 65 && num < 91) && !(num >= 97 && num < 123))){
        if(l > maxlen){ maxlen = l; jmax = j; raend = i+l;}
      }
    }
  }

  return jmax;
}

double getfloat(string text, long i)                           // Tries to get a number from a string
{
  long num, sign = 1;
  double fac, num2;

  if(text.substr(i,1) == "-"){ sign = -1; i++;}

  const char *cstr = text.substr(i,1).c_str();
  num = cstr[0]-48; if(!(num >= 0 && num <= 9)) return NULLVAL;

  num2 = 0;
  do{
    const char *cstr = text.substr(i,1).c_str();
    num = cstr[0]-48;
    if(num >= 0 && num <= 9) num2 = num2*10+num;
    else break;
    i++;
  }while(i < text.length());

  if(text.substr(i,1) == "."){   // looks for a decimal point
    i++;
    fac = 0.1;
    do{
      const char *cstr = text.substr(i,1).c_str();
      num = cstr[0]-48;
      if(num >= 0 && num <= 9){ num2 += fac*num; fac /= 10;}
      else break;
      i++;
    }while(i < text.length());
  }

  raend = i;
  return sign*num2;
}

void printcalc(long eq)                             // Outputs the calculation
{
  long i;

  cout << "Calculation: " << eqnstr[eq] << "\n";

  for(i = 0; i < neqncalc[eq]; i++){
    if(eqncalcop[eq][i] == ADD || eqncalcop[eq][i] == TAKE || eqncalcop[eq][i] == MULTIPLY || eqncalcop[eq][i] == DIVIDE || eqncalcop[eq][i] == POWER){
      switch(eqncalc1[eq][i]){
        case PARAM: cout << paramname[eqncalcnum1[eq][i]]; break;
        case POPNUM: cout << "'" << popnumname[eqncalcnum1[eq][i]] << "'"; break;
        case REG: cout << "R" << eqncalcnum1[eq][i]; break;
        case NUMERIC: cout <<  numeric[eqncalcnum1[eq][i]]; break;
      }
    }

    cout <<  " ";

    switch(eqncalcop[eq][i]){
      case ADD: cout <<  "+"; break;
      case TAKE: cout <<  "-"; break;
      case MULTIPLY: cout <<  "*"; break;
      case DIVIDE: cout <<  "/"; break;
      case POWER: cout <<  "^"; break;
      case EXPFUNC: cout <<  "exp"; break;
      case SINFUNC: cout <<  "sin"; break;
      case COSFUNC: cout <<  "cos"; break;
      case LOGFUNC: cout <<  "log"; break;
      case STEPFUNC: cout <<  "step"; break;
    }

    cout <<  " ";

    switch(eqncalc2[eq][i]){ 
      case PARAM: cout << paramname[eqncalcnum2[eq][i]]; break;
      case POPNUM: cout << "'" << popnumname[eqncalcnum2[eq][i]] << "'"; break;
      case REG: cout <<  "R" << eqncalcnum2[eq][i]; break;
      case NUMERIC: cout <<  numeric[eqncalcnum2[eq][i]]; break;
    }

    cout <<  "  > ";

    cout <<  "R" << eqncalcstore[eq][i] << ", ";
  }

  switch(eqncalcans[eq]){
    case PARAM: cout << paramname[eqncalcansnum[eq]]; break;
    case POPNUM: cout << popnumname[eqncalcansnum[eq]]; break;
    case REG: cout <<  "R" << eqncalcansnum[eq]; break;
    case NUMERIC: cout <<  numeric[eqncalcansnum[eq]]; break;
  }
  cout <<  " Answer\n";
}

double calculate(long eq, double *popnum, double *paramval)
{
  long i;
  double num, num1, num2;

  for(i = 0; i < neqncalc[eq]; i++){
    switch(eqncalc1[eq][i]){
      case PARAM: num1 = paramval[eqncalcnum1[eq][i]]; break;
      case POPNUM: num1 = popnum[eqncalcnum1[eq][i]]; break;
      case REG: num1 = regcalc[eqncalcnum1[eq][i]]; break;
      case NUMERIC: num1 = numeric[eqncalcnum1[eq][i]]; break;
    }

    switch(eqncalc2[eq][i]){
      case PARAM: num2 = paramval[eqncalcnum2[eq][i]]; break;
      case POPNUM: num2 = popnum[eqncalcnum2[eq][i]]; break;
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
    case PARAM: return paramval[eqncalcansnum[eq]]; break;
    case POPNUM: return popnum[eqncalcansnum[eq]]; break;
    case REG: return regcalc[eqncalcansnum[eq]]; break;
    case NUMERIC: return numeric[eqncalcansnum[eq]]; break;
  }
}

double calculatenotdep(long eq, double *param)
{
  long i;
  double num, num1, num2;

  for(i = 0; i < neqncalc[eq]; i++){
    switch(eqncalc1[eq][i]){
      case PARAM: num1 = param[eqncalcnum1[eq][i]]; break;
      case POPNUM: emsg("not dep1"); break;
      case REG: num1 = regcalc[eqncalcnum1[eq][i]]; break;
      case NUMERIC: num1 = numeric[eqncalcnum1[eq][i]]; break;
    }

    switch(eqncalc2[eq][i]){
      case PARAM: num2 = param[eqncalcnum2[eq][i]]; break;
      case POPNUM: cout << eqnstr[eq] << " eq\n"; emsg("not depw"); break;
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
    case PARAM: return param[eqncalcansnum[eq]]; break;
    case POPNUM: cout << eqnstr[eq] << " eq\n"; emsg("not dep3"); break;
    case REG: return regcalc[eqncalcansnum[eq]]; break;
    case NUMERIC: return numeric[eqncalcansnum[eq]]; break;
  }
  emsg("eqn prob");
}

vector<long> getallcomp(string st)                                       // Gets all the comprtments from a string
{
  //long k, j, cl2, flag, c, filt[nclass], count[nclass];
  long k, j, flag, c, cl, filt[nclass], count[nclass];
  vector<string> vec;
  vector<long> state;

  for(cl = 0; cl < nclass; cl++) filt[cl] = -1;

  vec = getcommasep(st);
  for(k = 0; k < vec.size(); k++){
    for(cl = 0; cl < nclass; cl++){ j = getclassval(cl,vec[k]); if(j >= 0){ filt[cl] = j; break;}}
    if(cl == nclass) emsg(st+"cannot find filter");
  }

  for(cl = 0; cl < nclass; cl++){ if(filt[cl] == -1) count[cl] = 0; else count[cl] = filt[cl];}

  do{
    c = 0; for(cl = 0; cl < nclass; cl++) c += classmult[cl]*count[cl];
    state.push_back(c);

    j = 0;
    do{
      flag = 0;
      if(filt[j] != -1){ j++; flag = 1;}
      else{
        count[j]++; if(count[j] >= nclassval[j]){ count[j] = 0; j++; flag = 1;}
      }
    }while(flag == 1 && j < nclass);
  }while(j < nclass);

  return state;
}

vector<string> getcommasep(string a)                                                    // takes a string and splits it up by the commas
{
  long j, jst;
  vector<string> commasep;

  j = 0;
  while(j < a.length()){
    jst = j;
    while(j < a.length() && a.substr(j,1) != ",") j++;
    commasep.push_back(a.substr(jst,j-jst));
    j++;
  }

  return commasep;
}


long getclassval(long cl, string st)                                                // Returns the value of a classification from a string
{
  long j;

  for(j = 0; j < nclassval[cl]; j++) if(classval[cl][j] == st) return j;
  return -1;
}


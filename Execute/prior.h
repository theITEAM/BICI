double Chain::priorcalc()                                                           // The total prior probability
{
  long pr, i, c, s;
  double Lpri = 0;

  if(tbirthfl == 1){ for(i = 0; i < nindtot; i++) Lpri += agprior[compval[tra[indev[i][0].tr].cf][agecl]];}
  for(pr = 0; pr < nprior; pr++) Lpri += priorprob(pr);

  for(s = 0; s < nsmooth; s++){ for(i = 0; i < nsmoothparam[s]-1; i++) Lpri += priorsmooth(s,i);}
  
  return Lpri;
}

double Chain::priorsmooth(long s, long j)
{
  double d;

  if(j < 0 || j >= nsmoothparam[s]-1) return 0;

  switch(smoothtype[s]){
    case SMOOTH:
      d = param[smoothparam[s][j+1]] - param[smoothparam[s][j]];
      return -smoothval[s][j]*d*d;

    case LOGSMOOTH:
      if(param[smoothparam[s][j+1]] < 0 || param[smoothparam[s][j]] < 0){
        stringstream ss; ss << "When log smoothing '" << smoothname[s] << "' cannot have negative parameter value.";
      }
      d = log(param[smoothparam[s][j+1]]) - log(param[smoothparam[s][j]]);
      return -smoothval[s][j]*d*d;
  }
  emsg("Smooth type problem");
}

double boundprob(double val, string st, string s)
{
  if(simon == 1){
    if(st.substr(0,1) == "[") st = st.substr(1,st.length()-1);
    if(st.substr(st.length()-1,1) == "]") st = st.substr(0,st.length()-1);
    stringstream ss;
    if(val < 0) ss << "The value " << st << " = " << val << " in the " << s << " distribution is negative.";
    else ss << "The value for " << st << " in the " << s << " distribution is too small.";
    emsg(ss.str());
  }
  else sampill = 1;

  return 1;
}

double Chain::priorsamp(long pr)
{
  double mean, sd, alpha, beta, lam, k, val;

  switch(priortype[pr]){ 
    case FIX:
      return priorminval[pr];

    case FLAT:
      return priorminval[pr] + ran()*(priormaxval[pr]-priorminval[pr]);

    case GAMMA:
      mean = calculatenotdep(prioreq1[pr],param); sd = calculatenotdep(prioreq2[pr],param);
      if(sd < 0) return boundprob(sd,eqnstr[prioreq2[pr]],"gamma");
      alpha = mean*mean/(sd*sd); beta = mean/(sd*sd);
      val = gammasamp(alpha,beta);//
      if(alpha < gammalim || val == 0){
        if(simon == 1) emsg("The shape parameter for the gamma distribution is too small");
        else{ sampill = 1; return 1;}
      }
      return val;

    case NORMAL:
      mean = calculatenotdep(prioreq1[pr],param); sd = calculatenotdep(prioreq2[pr],param);
      if(sd < 0) return boundprob(sd,eqnstr[prioreq2[pr]],"normal");
      return normal(mean,sd);

    case LOGNORMAL:
      mean = calculatenotdep(prioreq1[pr],param); sd = calculatenotdep(prioreq2[pr],param);
      if(sd < 0) return boundprob(sd,eqnstr[prioreq2[pr]],"log-normal");
      val = exp(normal(mean,sd));
      if(val == 0 || isnan(val)){
        if(simon == 1) emsg("The value for the log-normal distribution is out of range.");
        else{ sampill = 1; return 1;}
      }
      return val;

    case EXPO:
      return -log(ran())/calculatenotdep(prioreq1[pr],param);

    case BETA:
      alpha = calculatenotdep(prioreq1[pr],param); beta = calculatenotdep(prioreq2[pr],param);

      if(alpha < gammalim) return boundprob(alpha,eqnstr[prioreq1[pr]],"beta");
      if(beta < gammalim) return boundprob(beta,eqnstr[prioreq2[pr]],"beta");
      val = betasamp(alpha,beta);
      if(val == 0 || val == 1 || isnan(val)){
        if(simon == 1) emsg("The value for the beta distribution is out of range.");
        else{ sampill = 1; return 1;}
      }
      return val;

    case WEIBULL:
      lam = calculatenotdep(prioreq1[pr],param); k = calculatenotdep(prioreq2[pr],param); 
      if(lam <= 0) return boundprob(lam,eqnstr[prioreq1[pr]],"Weibull");
      if(k < gammalim) return boundprob(k,eqnstr[prioreq2[pr]],"Weibull");
      val = weibullsamp(lam,k);
      if(val == 0 || isnan(val)){
        if(simon == 1) emsg("The value for the Weibull distribution is out of range.");
        else{ sampill = 1; return 1;}
      }
      return val;

    default: emsg("Error sampling"); break;
  }
}

double Chain::priorprob(long pr)                                                   // Calculates the prior probability
{
  long p;
  double mean, sd, rate, alpha, beta, lam, k;

  if(pr == -1) emsg("prior problem");
  p = priorparam[pr];

  switch(priortype[pr]){
    case FIX:
      return 0;
      break;

    case FLAT:
      if(param[p] < priorminval[pr] || param[p] > priormaxval[pr]) return -large;
      return 1.0/(priormaxval[pr]-priorminval[pr]);
      break;

    case GAMMA:
      mean = calculatenotdep(prioreq1[pr],param); sd = calculatenotdep(prioreq2[pr],param);
      alpha = mean*mean/(sd*sd); beta = mean/(sd*sd);
      if(mean <= 0 || sd <= 0 || param[p] <= 0 || alpha < 0.1) return -large;
      return gammaprob(param[p],alpha,beta);
      break;

    case NORMAL:
      mean = calculatenotdep(prioreq1[pr],param); sd = calculatenotdep(prioreq2[pr],param);
      if(sd <= 0) return -large;
      return normalprob(param[p],mean,sd*sd);
      break;

    case LOGNORMAL:
      mean = calculatenotdep(prioreq1[pr],param); sd = calculatenotdep(prioreq2[pr],param); 
      if(sd <= 0 || param[p] <= 0) return -large;
      return lognormalprob(param[p],mean,sd*sd);
      break;

    case EXPO:
      rate = calculatenotdep(prioreq1[pr],param);
      return expprob(param[p],rate);
      break;

    case BETA:
      alpha = calculatenotdep(prioreq1[pr],param); beta = calculatenotdep(prioreq2[pr],param);
      if(alpha < 0.1 || beta < 0.1 || param[p] <= 0 || param[p] >= 1) return -large;
      return betaprob(param[p],alpha,beta);
      break;

    case WEIBULL:
      lam = calculatenotdep(prioreq1[pr],param); k = calculatenotdep(prioreq2[pr],param);
      if(lam < 0 || k < 0 || param[p] <= 0) return -large;
      return weibullprob(param[p],lam,k);
      break;

    default: emsg("prior problem"); break;
  }
}

void addparampriordep(long p, long pr)
{
  long j;
  for(j = 0; j < parampriordep[p].size(); j++) if(parampriordep[p][j] == pr) return;
  parampriordep[p].push_back(pr);
}

void priorinit()
{
  long pr, eq, p, j;
  double dt;

  parampriordep.resize(nparam);
  for(pr = 0; pr < nprior; pr++){
    switch(priortype[pr]){
      case GAMMA: case NORMAL: case LOGNORMAL: case  BETA: case WEIBULL:
        eq = prioreq1[pr]; for(j = 0; j < neq_param[eq]; j++) addparampriordep(eq_param[eq][j],pr);
        eq = prioreq2[pr]; for(j = 0; j < neq_param[eq]; j++) addparampriordep(eq_param[eq][j],pr);
        break;

      case EXPO:
        eq = prioreq1[pr]; for(j = 0; j < neq_param[eq]; j++) addparampriordep(eq_param[eq][j],pr);
        break;
    }
  }

  for(p = 0; p < nparam; p++) nparampriordep.push_back(parampriordep[p].size());

  agprior.resize(nclassval[agecl]);             // Provides a prior such that the prior probability of being in any one of the age classifications is the same
  for(j = 0; j < nclassval[agecl]; j++){
    if(j == 0) dt = age[j]; else dt = age[j]-age[j-1];
    //agprior[j] = log(1.0/(nclassval[agecl]*dt));
    agprior[j] = log(1.0/(dt));
  }
}


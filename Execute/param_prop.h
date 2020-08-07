void Chain::param_prop()                                        // Makes random walk proposals on model parameters
{
  long p, pr, c, s, j;
  double param_st, param_new, al, Ltoti, Ltotf, sum, num[ncompswa], dparam;

  for(p = 0; p < nparam; p++){
    param_st = param[p];
    param_new = param_st + normal(0,jump_param[p]);

    if(notinsideprior(p,param_new) == 1) al = -1;
    else{
      Ltoti = L();
      changeparam(1,p,param_new);
      Ltotf = L();
      al = exp(Ltotf-Ltoti);
      //cout << samp << " " << jump_param[p] << " " << al << " " << Linm << " " << Linmst << " " << Ltoti << " " << Ltotf  <<"al\n"; 
    }

    ntr_param[p]++;
    if(ran() < al){
      nac_param[p]++;
      if(samp < burnin){ if(samp < 50) jump_param[p] *= 1.4; else jump_param[p] *= 1.1;}
    }
    else{
      if(al >= 0) changeparam(0,p,param_st); 
      if(samp < burnin){ if(samp < 50) jump_param[p] *= 0.8; else jump_param[p] *= 0.95;}
    }
  }

  for(s = 0; s < nsmooth; s++){             // makes changes to all smooth parameters at the same time
    dparam = normal(0,jump_smooth[s]);

    for(j = 0; j < nsmoothparam[s]; j++){
      p = smoothparam[s][j]; if(notinsideprior(p,param[p]+dparam) == 1) break;
    }

    if(j < nsmoothparam[s]) al = -1;
    else{
      Ltoti = L();
      for(j = 0; j < nsmoothparam[s]; j++){ p = smoothparam[s][j]; changeparam(1,p,param[p]+dparam);}
      Ltotf = L();
      al = exp(Ltotf-Ltoti);
    }

    ntr_smooth[s]++;
    if(ran() < al){
      nac_smooth[s]++;
      if(samp < burnin){ if(samp < 50) jump_smooth[s] *= 1.4; else jump_smooth[s] *= 1.1;}
    }
    else{
      if(al >= 0){ for(j = 0; j < nsmoothparam[s]; j++){ p = smoothparam[s][j]; changeparam(0,p,param[p]-dparam);}}
      if(samp < burnin){ if(samp < 50) jump_smooth[s] *= 0.8; else jump_smooth[s] *= 0.95;}
    }
  }

  sum = 0; for(c = 0; c < ncompswa; c++){ num[c] = gammasamp(ncompinit[c]+1+dirchal[c]-1,1); sum += num[c];}

  for(c = 0; c < ncompswa; c++) probinit[c] = num[c]/sum;
  //for(c = 0; c < ncompswa; c++) cout << c<< " " << probinit[c] << " init\n";

  Liinit = logsum[nindinit]; for(c = 0; c < ncompswa; c++) Liinit += ncompinit[c]*log(probinit[c]);
}

void Chain::paramnorm_prop()                                        // Makes a proposal which is normal fitted to posterior
{
  long p, pr;
  double param_st, param_new, al, Ltoti, Ltotf, probif, probfi;

  for(p = 0; p < nparam; p++){
    pr = paramprior[p];
    if(priortype[pr] != FIX){
      ntr_paramnorm[p]++;

      param_st = param[p];
      paramgrad(p);
      if(pvar > 0){
        param_new = normal(pmean,sqrt(pvar));

        if(notinsideprior(p,param_new) == 0){
          probif = normalprob(param_new,pmean,pvar);

          Ltoti = L();
          changeparam(1,p,param_new);
          Ltotf = L();

          if(Ltotf - Ltoti < -large/2) al = 0;
          else{
            paramgrad(p);
            if(pvar < 0) al = 0;
            else{
              probfi = normalprob(param_st,pmean,pvar);
              al = exp(Ltotf-Ltoti+probfi-probif);
            }
          }

          if(ran() < al) nac_paramnorm[p]++;
          else changeparam(0,p,param_st);
        }
      }
    }
  }
}

void Chain::param_sampprop()                                        // Makes proposals on model parameters by sampling from the prior
{
  long p, pr;
  double  param_st, param_new, al, Ltoti, Ltotf, probfi, probif;

  for(p = 0; p < nparam; p++){
    pr = paramprior[p];
    param_st = param[p]; probfi = priorprob(pr);
    param_new = priorsamp(pr);

    Ltoti = L();
    changeparam(1,p,param_new);
    Ltotf = L();
    probif = priorprob(pr);
    al = exp(Ltotf-Ltoti+probfi-probif);

    ntr_paramsamp[p]++;
    if(ran() < al) nac_paramsamp[p]++;
    else changeparam(0,p,param_st);
  }
}

void Chain::param_gradprop()                                        // Makes proposals to increase posterior probability
{
  long p, pr, j;
  double  param_st, param_new, al, Ltoti, Ltotf, probfi, probif, dLpri;

  for(p = 0; p < nparam; p++){
    pr = paramprior[p];
    if(priortype[pr] != FIX){
      param_st = param[p];
      paramgrad(p);
      param_new = param_st + pgrad*pvel[p];

      dLpri = 0;
      if(pr >= 0) dLpri -= priorprob(pr); for(j = 0; j < nparampriordep[p]; j++) dLpri -= priorprob(parampriordep[p][j]);
      param[p] = param_new;
      if(pr >= 0) dLpri += priorprob(pr); for(j = 0; j < nparampriordep[p]; j++) dLpri += priorprob(parampriordep[p][j]);
      param[p] = param_st;

      if(dLpri > -100){ 
        Ltoti = L();
        changeparam(1,p,param_new);
        Ltotf = L();
        if(Ltotf < Ltoti){ changeparam(0,p,param_st); pvel[p] *= 0.5;}
        else pvel[p] *= 2;
      }
    }
  }
}

void Chain::changeparam(long dir, long p, double param_new)              // Makes changes to parameters
{
  long j, d, eq, pr, m, s;
  double t, tt, valold, valnew, valold2, valnew2, pd, op;
  long e, ee, i;
  
  //timeparam[p] -=clock();
  pr = paramprior[p];

  if(pr >= 0) Lpri -= priorprob(pr);
  for(j = 0; j < nparampriordep[p]; j++) Lpri -= priorprob(parampriordep[p][j]);
  for(j = 0; j < paramsmooth[p].size(); j++){ s = paramsmooth[p][j]; i = paramsmoothi[p][j]; Lpri -= priorsmooth(s,i)+priorsmooth(s,i-1);}

  param[p] = param_new;

  if(pr >= 0) Lpri += priorprob(pr);
  for(j = 0; j < nparampriordep[p]; j++) Lpri += priorprob(parampriordep[p][j]);
  for(j = 0; j < paramsmooth[p].size(); j++){ s = paramsmooth[p][j]; i = paramsmoothi[p][j]; Lpri += priorsmooth(s,i)+priorsmooth(s,i-1);}

  for(j = 0; j < nparam_dep[p]; j++){
    d = param_dep[p][j];
    eq = transdepeq[d];

    e = depeqdiv[d][0];
    t = 0;
    do{
      ee = depeq_evnext[e];
      tt = depeq_evt[ee];

      valold = depeq_evval[e];
      valnew = ratecalcdep(d,depeq_evpopnum[e],param);
      depeq_evval[e] = valnew;

      Liexp += (valnew-valold)*depeq_evn[e]*(tt-t);
      Lir += (log(valnew)-log(valold))*ndepeq_ev_ev[e];

      e = ee; t = tt;
    }while(t < tmax);
  }

  for(j = 0; j < nparam_notdep[p]; j++){
    d = param_notdep[p][j]; 
    valold = transnotdepeq_val[d];
    valnew = ratecalcnotdep(transnotdepeq[d],param);
    transnotdepeq_val[d] = valnew;
    Lir += (log(valnew)-log(valold))*transnotdepeq_num[d];
    Liexp += (valnew-valold)*transnotdepeq_dt[d];
  }

  for(j = 0; j < nparam_cope[p]; j++){
    d = param_cope[p][j];

    valold = cope_val[d];
    op = calcobsprob(cope[d]);
    valnew = log(op);
    cope_val[d] = valnew;
    cope_op[d] = op;

    Lob_st += (valnew-valold)*cope_num[d];
  }

  for(j = 0; j < nparam_cpe[p]; j++){
    d = param_cpe[p][j];

    valold = cpe_val[d];
    valold2 = cpe_valoneminus[d];

    pd = calcobsprob(cpe[d]);

    valnew = log(pd);
    valnew2 = log(1-pd);

    cpe_pd[d] = pd;
    cpe_val[d] = valnew;
    cpe_valoneminus[d] = valnew2;

    Lob_pd += (valnew-valold)*cpe_num[d] + (valnew2-valold2)*cpe_oneminus[d];
  }

  for(j = 0; j < nparam_capev[p]; j++){
    d = param_capev[p][j];

    valold = capev_val[d];
    valold2 = capev_valoneminus[d];

    pd = calcobsprob(capevprobeqn[d]);

    valnew = log(pd);
    valnew2 = log(1-pd);

    capev_val[d] = valnew;
    capev_valoneminus[d] = valnew2;

    if(capev_num[d] > 0) Lev_pd += (valnew-valold)*capev_num[d];
    if(capev_oneminus[d] > 0) Lev_pd += (valnew2-valold2)*capev_oneminus[d];
  }

  for(j = 0; j < nparam_derm[p]; j++){
    m = param_derm[p][j];
    Lder -= derL[m]; derL[m] = dercalc(m,dermpopnum[m],param); Lder += derL[m];
  }

  if(param_nmfl[p] == 1){
    if(dir == 1){
      Linmst = Linm;
      setnmeq();
      Linm = 0; for(i = 0; i < nindtot; i++) Linm += likenm(i);
    }
    else{
      Linm = Linmst;
      setnmeq();
    }
  }
  //timeparam[p] += clock();
}

void Chain::paramgrad(long p)              // Calculates the gradient and curvature on a parameter
{
  long j, d, eq, pr,  e, ee, i;
  double t, tt, valold, valnewup, valnewdo, valold2, valnewup2, valnewdo2, pd, op, dLup, dLdo, va, val, pcurve, dd, Linmtemp;

  //timegradparam[p] -=clock();
  dLup = 0; dLdo = 0;

  val = param[p]; dd = val/1000;
  pr = paramprior[p];

  if(pr >= 0){ va = priorprob(pr); dLup -= va; dLdo -= va;}
  for(j = 0; j < nparampriordep[p]; j++){ va = priorprob(parampriordep[p][j]); dLup -= va; dLdo -= va;}

  param[p] = val + dd;
  if(pr >= 0) dLup += priorprob(pr); 
  for(j = 0; j < nparampriordep[p]; j++) dLup += priorprob(parampriordep[p][j]);

  param[p] = val - dd; 
  if(pr >= 0) dLdo += priorprob(pr); 
  for(j = 0; j < nparampriordep[p]; j++) dLdo += priorprob(parampriordep[p][j]);

  for(j = 0; j < nparam_dep[p]; j++){
    d = param_dep[p][j];
    eq = transdepeq[d];

    e = depeqdiv[d][0];
    t = 0;
    do{
      ee = depeq_evnext[e];
      tt = depeq_evt[ee];

      valold = depeq_evval[e];
      param[p] = val + dd; valnewup = ratecalcdep(d,depeq_evpopnum[e],param);
      param[p] = val - dd; valnewdo = ratecalcdep(d,depeq_evpopnum[e],param);

      dLup += -(valnewup-valold)*depeq_evn[e]*(tt-t) + (log(valnewup)-log(valold))*ndepeq_ev_ev[e];
      dLdo += -(valnewdo-valold)*depeq_evn[e]*(tt-t) + (log(valnewdo)-log(valold))*ndepeq_ev_ev[e];
      e = ee; t = tt;
    }while(t < tmax);
  }

  for(j = 0; j < nparam_notdep[p]; j++){
    d = param_notdep[p][j];
    valold = transnotdepeq_val[d];
    param[p] = val + dd; valnewup = ratecalcnotdep(transnotdepeq[d],param);
    param[p] = val - dd; valnewdo = ratecalcnotdep(transnotdepeq[d],param);

    dLup += (log(valnewup)-log(valold))*transnotdepeq_num[d] - (valnewup-valold)*transnotdepeq_dt[d];
    dLdo += (log(valnewdo)-log(valold))*transnotdepeq_num[d] - (valnewdo-valold)*transnotdepeq_dt[d];
  }

  for(j = 0; j < nparam_cope[p]; j++){
    d = param_cope[p][j];

    valold = cope_val[d];
    param[p] = val + dd; valnewup = log(calcobsprob(cope[d]));
    param[p] = val - dd; valnewdo = log(calcobsprob(cope[d]));
    dLup += (valnewup-valold)*cope_num[d];
    dLdo += (valnewdo-valold)*cope_num[d];
  }

  for(j = 0; j < nparam_cpe[p]; j++){
    d = param_cpe[p][j];

    valold = cpe_val[d]; valold2 = cpe_valoneminus[d];
    param[p] = val + dd; pd = calcobsprob(cpe[d]); valnewup = log(pd); valnewup2 = log(1-pd);
    param[p] = val - dd; pd = calcobsprob(cpe[d]); valnewdo = log(pd); valnewdo2 = log(1-pd);

    dLup += (valnewup-valold)*cpe_num[d] + (valnewup2-valold2)*cpe_oneminus[d];
    dLdo += (valnewdo-valold)*cpe_num[d] + (valnewdo2-valold2)*cpe_oneminus[d];
  }

  for(j = 0; j < nparam_capev[p]; j++){
    d = param_capev[p][j];

    valold = capev_val[d]; valold2 = capev_valoneminus[d];
    param[p] = val + dd; pd = calcobsprob(capevprobeqn[d]); valnewup = log(pd); valnewup2 = log(1-pd);
    param[p] = val - dd; pd = calcobsprob(capevprobeqn[d]); valnewdo = log(pd); valnewdo2 = log(1-pd);

    if(capev_num[d] > 0){ dLup += (valnewup-valold)*capev_num[d];  dLdo += (valnewdo-valold)*capev_num[d];}
    if(capev_oneminus[d] > 0){ dLup += (valnewup2-valold2)*capev_oneminus[d]; dLdo += (valnewdo2-valold2)*capev_oneminus[d];}
  }

  if(param_nmfl[p] == 1){
    param[p] = val + dd; setnmeq(); Linmtemp = 0; for(i = 0; i < nindtot; i++) Linmtemp += likenm(i); dLup += Linmtemp - Linm;
    param[p] = val - dd; setnmeq(); Linmtemp = 0; for(i = 0; i < nindtot; i++) Linmtemp += likenm(i); dLdo += Linmtemp - Linm;
  }

  pgrad = (dLup-dLdo)/(2*dd);
  pcurve = (dLup+dLdo)/(dd*dd);
 
  if(pcurve > 0){ pmean = 0; pvar = 0;}
  else{
    pmean = val - pgrad/pcurve;
    pvar = -1.0/pcurve;
    if(pmean*pmean > 100*val*val){ pmean = 0; pvar = 0;}
  }
  param[p] = val; if(param_nmfl[p] == 1) setnmeq();
  //timegradparam[p] +=clock();
}

void paraminit()                                                 // Initialisation for parameter proposals
{
  long p, j, d, eq, k, m;

  for(p = 0; p < nparam; p++){
    param_dep.push_back(vector<long>());
    for(d = 0; d < ntransdepeq; d++){
      eq = transdepeq[d];
      for(j = 0; j < neq_param[eq]; j++) if(eq_param[eq][j] == p) break;

      if(j < neq_param[eq]) param_dep[p].push_back(d);
    }
    nparam_dep.push_back(param_dep[p].size());

    param_notdep.push_back(vector<long>());
    for(d = 0; d < ntransnotdepeq; d++){
      eq = transnotdepeq[d];
      for(j = 0; j < neq_param[eq]; j++) if(eq_param[eq][j] == p) break;

      if(j < neq_param[eq]) param_notdep[p].push_back(d);
    }
    nparam_notdep.push_back(param_notdep[p].size());

    param_cope.push_back(vector<long>());
    for(d = 0; d < ncope; d++){
      eq = cope[d];
      for(j = 0; j < neq_param[eq]; j++) if(eq_param[eq][j] == p) break;

      if(j < neq_param[eq]) param_cope[p].push_back(d);
    }
    nparam_cope.push_back(param_cope[p].size());

    param_cpe.push_back(vector<long>());
    for(d = 0; d < ncpe; d++){
      eq = cpe[d];
      for(j = 0; j < neq_param[eq]; j++) if(eq_param[eq][j] == p) break;

      if(j < neq_param[eq]) param_cpe[p].push_back(d);
    }
    nparam_cpe.push_back(param_cpe[p].size());

    param_capev.push_back(vector<long>());
    for(d = 0; d < ncapev; d++){
      eq = capevprobeqn[d];
      for(j = 0; j < neq_param[eq]; j++) if(eq_param[eq][j] == p) break;

      if(j < neq_param[eq]) param_capev[p].push_back(d);
    }
    nparam_capev.push_back(param_capev[p].size());

    param_derm.push_back(vector<long>());
    for(d = 0; d < nderive; d++){
      eq = derive[d];
      for(j = 0; j < neq_param[eq]; j++) if(eq_param[eq][j] == p) break;

      if(j < neq_param[eq]){
        for(m = 0; m < nderm; m++){
          if(derm[m] == d) param_derm[p].push_back(m);
        }
      }
    }
    nparam_derm.push_back(param_derm[p].size());
  }

  for(p = 0; p < nparam; p++) param_nmfl.push_back(0);
  for(k = 0; k < nnmeq; k++){ eq = nmeq[k]; for(j = 0; j < neq_param[eq]; j++) param_nmfl[eq_param[eq][j]] = 1;}

  
  /*
  for(p = 0; p < -nparam; p++){
  logout << paramname[p] << ": ";
  for(j = 0; j < nparam_dep[p]; j++) logout << eqnstr[transdepeq[param_dep[p][j]]] << ", ";
  for(j = 0; j < nparam_notdep[p]; j++) logout << eqnstr[transnotdepeq[param_notdep[p][j]]] << ", ";
  for(j = 0; j < nparam_cope[p]; j++) logout << eqnstr[cope[param_cope[p][j]]] << ", ";
  for(j = 0; j < nparam_cpe[p]; j++) logout << eqnstr[cpe[param_cpe[p][j]]] << ", ";
  logout << "\n";
}
  */
}

long notinsideprior(long p, double val)       // Determines if a value is within the prior
{
  long pr;

  pr = paramprior[p];
  if(pr >= 0){
    switch(priortype[pr]){
      case FIX: return 1;
      case FLAT: if(val > priormaxval[pr] || val < priorminval[pr]) return 1; break;
      case GAMMA: case EXPO: case WEIBULL: if(val < 0) return 1; break;
      case BETA: if(val > 1|| val < 0) return 1; break;
      case NORMAL: break;
      case LOGNORMAL: break;
      default: emsg("prior prob param"); break;
    }
  }
  return 0;
}


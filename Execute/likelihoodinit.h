void likelihoodinit()                                          // Initialises the quantities need to make fast likelihood calculations
{
  long eq, p, i, tr, k, j, c, d, jmax, ci, cf, n, flag, kmax, nn, cl, f, ki, kf, cff;
  long popnum_flag[npopnum], param_flag[nparam];
  vector< vector<double> > dpopnum;
  vector<double> dpopvec;
  double dd;
  vector<long> list;
  vector< vector<long> > mapdep, mapnotdep;
  vector< vector<long> > mapdepi, mapnotdepi;

  long ref, refi, reff;
  vector< vector <long> > transdepeq_n, transnotdepeq_n;

  for(eq = 0; eq < eqnstr.size(); eq++){
    eq_popnum.push_back(vector<long>());
    eq_param.push_back(vector<long>());

    for(p = 0; p < npopnum; p++) popnum_flag[p] = 0;
    for(p = 0; p < nparam; p++) param_flag[p] = 0;

    for(i = 0; i < eqncalc1[eq].size(); i++){
      switch(eqncalcop[eq][i]){
        case ADD: case TAKE: case MULTIPLY: case DIVIDE: case POWER:
          switch(eqncalc1[eq][i]){
            case PARAM: param_flag[eqncalcnum1[eq][i]] = 1; break;
            case POPNUM: popnum_flag[eqncalcnum1[eq][i]] = 1; break;
          }
          break;
      }

      switch(eqncalc2[eq][i]){
        case PARAM: param_flag[eqncalcnum2[eq][i]] = 1; break;
        case POPNUM: popnum_flag[eqncalcnum2[eq][i]] = 1; break;
      }
    }

    switch(eqncalcans[eq]){
      case PARAM: param_flag[eqncalcansnum[eq]] = 1; break;
      case POPNUM: popnum_flag[eqncalcansnum[eq]] = 1; break;
    }

    for(p = 0; p < npopnum; p++){
      if(popnum_flag[p] == 1) eq_popnum[eq].push_back(p);
    }
    neq_popnum.push_back(eq_popnum[eq].size());

    for(p = 0; p < nparam; p++){
      if(param_flag[p] == 1) eq_param[eq].push_back(p);
    }
    neq_param.push_back(eq_param[eq].size());
  }

  /*(
  for(eq = 0; eq < eqnstr.size(); eq++){
    cout << eqnstr[eq] << ":  popnum ";
    for(i = 0; i < eq_popnum[eq].size(); i++) cout << popnumname[eq_popnum[eq][i]] << ", "; 

    cout << "    params: ";
    for(i = 0; i < eq_param[eq].size(); i++) cout << paramname[eq_param[eq][i]] << ",  ";
    cout << "\n";
  }
  */

  // Creates look up table
  for(eq = 0; eq < eqnstr.size(); eq++){ transdep.push_back(-1); transdepref.push_back(-1);}

  for(tr = 0; tr < ntra; tr++){
    if(tra[tr].type == EXP_TR && tra[tr].like == 1){
      eq = tra[tr].eq;
      if(neq_popnum[eq] > 0){
        k = transdep[eq];
        if(k == -1){
          k = transdepeq.size();
          transdepeq.push_back(eq);
          transdep[eq] = 1;
          transdepref[eq] = k; 
        }
      }
      else{
        k = transdep[eq];
        if(k == -1){
          k = transnotdepeq.size();
          transnotdepeq.push_back(eq);
          transdep[eq] = 0;
          transdepref[eq] = k;
        }
      }
    }
  }

  ntransdepeq = transdepeq.size();
  ntransnotdepeq = transnotdepeq.size();

  //for(j = 0; j < ntransdepeq; j++) cout << eqnstr[transdepeq[j]] << " dep\n";
  //for(j = 0; j < ntransnotdepeq; j++) cout << eqnstr[transnotdepeq[j]] << " notdep\n";

  dpopnum.resize(ncomp+1);     // the change in population number
  for(c = 0; c <= ncomp; c++){
    dpopnum[c].resize(npopnum);
    for(p = 0; p < npopnum; p++){ 
      dd = 0;
      for(j = 0; j < popnumterm[p].size(); j++){
        if(popnumterm[p][j] == c) dd += popnumtermweight[p][j];
      }
      dpopnum[c][p] = dd;
    }
  }

  transdepeq_n.resize(ntransdepeq); for(d = 0; d < ntransdepeq; d++){ transdepeq_n[d].resize(ncomp); for(c = 0; c < ncomp; c++) transdepeq_n[d][c] = 0;}
  transnotdepeq_n.resize(ntransnotdepeq); for(d = 0; d < ntransnotdepeq; d++){ transnotdepeq_n[d].resize(ncomp); for(c = 0; c < ncomp; c++) transnotdepeq_n[d][c] = 0;}

  //comptransdep.resize(ncomp);
  for(tr = 0; tr < ntra; tr++){
    if(tra[tr].type == EXP_TR && tra[tr].like == 1){
      eq = tra[tr].eq; ci = tra[tr].ci; d = transdepref[eq];
      //cout << eq << " " << traci[tr] << " "<<tracf[tr] << " eq\n";
      //cout << ci << " " << eqnstr[eq] <<" out\n";
      if(d >= 0 && ci != NOTALIVE){
        if(transdep[eq] == 1) transdepeq_n[d][ci]++;
        else transnotdepeq_n[d][ci]++;
      }
    }
  }

  transchref.resize(ncomp+1);
  for(ci = 0; ci <= ncomp; ci++) transchref[ci].resize(ncomp+1);

  ref = 0;
  ci = NOTALIVE;
  for(cf = 0; cf < ncomp; cf++){                                                       // first looks at the case when individuals are simply created
    transchref[ci][cf].push_back(ref);

    transdepeqch.push_back(vector<EQCH> ());
    transnotdepeqch.push_back(vector<NDEQCH> ());

    for(d = 0; d < ntransdepeq; d++){
      eq = transdepeq[d];

      n = transdepeq_n[d][cf];

      flag = 0;                         // works out change in population numbers
      if(n != 0) flag = 1;
      kmax = neq_popnum[eq];
      if(kmax > 0){
        dpopvec.resize(kmax);
        for(k = 0; k < kmax; k++){
          p = eq_popnum[eq][k];
          dd = dpopnum[cf][p];
          if(dd > -tiny && dd < tiny) dd = 0;
          if(dd != 0) flag = 2;
          dpopvec[k] = dd;
        }
      }

      if(flag != 0){
        EQCH eqch;
        eqch.d = d; eqch.n = n; if(flag == 2) eqch.valch = 1; else eqch.valch = 0; 
        if(kmax > 0){ eqch.popnum = new double[kmax]; for(k = 0; k < kmax; k++)  eqch.popnum[k] = dpopvec[k];}
        transdepeqch[ref].push_back(eqch);
      }
    }

    for(d = 0; d < ntransnotdepeq; d++){                                         // goes through not dependent equations
      eq = transnotdepeq[d];

      n = transnotdepeq_n[d][cf];
      if(n != 0){ NDEQCH eqch; eqch.d = d; eqch.n = n; transnotdepeqch[ref].push_back(eqch);}
    }
    ref++;
  }

  cf = NOTALIVE;
  for(ci = 0; ci < ncomp; ci++){                                                       // first looks at the case when individuals are simply created
    transchref[ci][cf].push_back(ref);

    transdepeqch.push_back(vector<EQCH> ());
    transnotdepeqch.push_back(vector<NDEQCH> ());

    for(d = 0; d < ntransdepeq; d++){
      eq = transdepeq[d];

      n = -transdepeq_n[d][ci];

      flag = 0;                         // works out change in population numbers
      if(n != 0) flag = 1;
      kmax = neq_popnum[eq];
      if(kmax > 0){
        dpopvec.resize(kmax);
        for(k = 0; k < kmax; k++){
          p = eq_popnum[eq][k];
          dd = -dpopnum[ci][p];
          if(dd > -tiny && dd < tiny) dd = 0;
          if(dd != 0) flag = 2;
          dpopvec[k] = dd;
        }
      }

      if(flag != 0){
        EQCH eqch;
        eqch.d = d; eqch.n = n; if(flag == 2) eqch.valch = 1; else eqch.valch = 0; 
        if(kmax > 0){ eqch.popnum = new double[kmax]; for(k = 0; k < kmax; k++) eqch.popnum[k] = dpopvec[k];}
        transdepeqch[ref].push_back(eqch);
      }
    }

    for(d = 0; d < ntransnotdepeq; d++){                                         // goes through not dependent equations
      eq = transnotdepeq[d];

      n = -transnotdepeq_n[d][ci];
      if(n != 0){ NDEQCH eqch; eqch.d = d; eqch.n = n; transnotdepeqch[ref].push_back(eqch);}
    }
    ref++;
  }

  mapdep.resize(ntransdepeq); mapdepi.resize(ntransdepeq);
  mapnotdep.resize(ntransnotdepeq); mapnotdepi.resize(ntransnotdepeq);

  for(ci = 0; ci < ncomp; ci++){                                                                        // considered transitions in the same classification
    for(cl = 0; cl < nclass; cl++){
      i = compval[ci][cl];

      for(f = 0; f < nclassval[cl]; f++){
        if(f != i){
          cf = ci + classmult[cl]*(f-i);

          if(nclassval[cl] > 200){  // If lots of compartments then does a simple death alive combination
            transchref[ci][cf].push_back(transchref[ci][NOTALIVE][0]);
            transchref[ci][cf].push_back(transchref[NOTALIVE][cf][0]);
          }
          else{                     // Otherwise optimises for a particular transition
            refi = transchref[NOTALIVE][ci][0];
            reff = transchref[NOTALIVE][cf][0];

            transchref[ci][cf].push_back(ref);

            transdepeqch.push_back(vector<EQCH> ());
            transnotdepeqch.push_back(vector<NDEQCH> ());

            list.clear();
            for(k = 0; k < transdepeqch[refi].size(); k++){
              d = transdepeqch[refi][k].d; if(mapdep[d].size() == 0) list.push_back(d); mapdep[d].push_back(refi); mapdepi[d].push_back(k);
            }
            for(k = 0; k < transdepeqch[reff].size(); k++){
              d = transdepeqch[reff][k].d; if(mapdep[d].size() == 0) list.push_back(d); mapdep[d].push_back(reff); mapdepi[d].push_back(k);
            }

            for(j = 0; j < list.size(); j++){
              d = list[j];
              ki = -1; kf = -1;
              for(k = 0; k < mapdep[d].size(); k++){
                if(mapdep[d][k] == refi) ki = mapdepi[d][k]; else{ if(mapdep[d][k] == reff) kf = mapdepi[d][k]; else emsg("pr");}
              }

              EQCH eqch;
              if(kf >= 0) n = transdepeqch[reff][kf].n; else n = 0;
              if(ki >= 0) n -= transdepeqch[refi][ki].n;
              eqch.n = n;

              flag = 0; if(eqch.n != 0) flag = 1;
              kmax = neq_popnum[transdepeq[d]];
              if(kmax > 0){
                eqch.popnum = new double[kmax];
                for(k = 0; k < kmax; k++){
                  if(kf >= 0) dd = transdepeqch[reff][kf].popnum[k]; else dd = 0;
                  if(ki >= 0) dd -= transdepeqch[refi][ki].popnum[k];
                  if(dd > -tiny && dd < tiny) dd = 0; if(dd != 0) flag = 2;
                  eqch.popnum[k] = dd;
                }
              }

              if(flag != 0){
                eqch.d = d; if(flag == 2) eqch.valch = 1; else eqch.valch = 0;
                transdepeqch[ref].push_back(eqch);
              }

              mapdep[d].clear(); mapdepi[d].clear();
            }

            list.clear();
            for(k = 0; k < transnotdepeqch[refi].size(); k++){
              d = transnotdepeqch[refi][k].d; if(mapnotdep[d].size() == 0) list.push_back(d); mapnotdep[d].push_back(refi); mapnotdepi[d].push_back(k);
            }
            for(k = 0; k < transnotdepeqch[reff].size(); k++){
              d = transnotdepeqch[reff][k].d; if(mapnotdep[d].size() == 0) list.push_back(d); mapnotdep[d].push_back(reff); mapnotdepi[d].push_back(k);
            }
  
            for(j = 0; j < list.size(); j++){
              d = list[j];
              ki = -1; kf = -1;
              for(k = 0; k < mapnotdep[d].size(); k++){ if(mapnotdep[d][k] == refi) ki = mapnotdepi[d][k];else{ if(mapnotdep[d][k] == reff) kf = mapnotdepi[d][k];else emsg("pr");}}
  
              if(kf >= 0) n = transnotdepeqch[reff][kf].n; else n = 0;
              if(ki >= 0) n -= transnotdepeqch[refi][ki].n;
              if(n != 0){ NDEQCH eqch; eqch.n = n; eqch.d = d; transnotdepeqch[ref].push_back(eqch);}
              mapnotdep[d].clear(); mapnotdepi[d].clear();
            }
            ref++;
          }
        }
      }
    }
  }

  for(ci = 0; ci < ncomp; ci++){       // goes through all the other transitions and makes them from a composite of those above
    for(cf = 0; cf < ncomp; cf++){
      if(ci != cf){
        if(transchref[ci][cf].size() == 0){
          c = ci;
          for(cl = 0; cl < nclass; cl++){
            i = compval[ci][cl]; f = compval[cf][cl];
            if(i != f){
              cff = c + (f-i)*classmult[cl];
              if(transchref[c][cff].size() != 1) break;
              transchref[ci][cf].push_back(transchref[c][cff][0]);
              c = cff;
            }
          }
          if(cl == nclass){
            if(c != cf) emsg("pro");
          }
          else{
            transchref[ci][cf].clear();
            transchref[ci][cf].push_back(transchref[ci][NOTALIVE][0]);
            transchref[ci][cf].push_back(transchref[NOTALIVE][cf][0]);
          }
        }
        if(transchref[ci][cf].size() == 0) emsg("no");
      }
    }
  }

  eqncalcfastnum1.resize(ntransdepeq);                       // Look up take used so that the eq timelines can easily calculate population number
  eqncalcfastnum2.resize(ntransdepeq);
  eqncalcansfastnum.resize(ntransdepeq);
  for(d = 0; d < ntransdepeq; d++){
    eq = transdepeq[d];

    eqncalcfastnum1[d].resize(neqncalc[eq]);
    eqncalcfastnum2[d].resize(neqncalc[eq]);

    for(i = 0; i < neqncalc[eq]; i++){
      switch(eqncalc1[eq][i]){
        case POPNUM: 
          p = eqncalcnum1[eq][i];
          j = 0; while(j < neq_popnum[eq] &&  p != eq_popnum[eq][j]) j++;
          if(j == neq_popnum[eq]) emsg("cannot find");
          eqncalcfastnum1[d][i] = j;
          break;
      }

      switch(eqncalc2[eq][i]){
        case POPNUM:
          p = eqncalcnum2[eq][i];
          j = 0; while(j < neq_popnum[eq] &&  p != eq_popnum[eq][j]) j++;
          if(j == neq_popnum[eq]) emsg("cannot find");
          eqncalcfastnum2[d][i] = j;
          break;
      }
    }

    switch(eqncalcans[eq]){
      case POPNUM:
        p = eqncalcansnum[eq];
        j = 0; while(j < neq_popnum[eq] && p != eq_popnum[eq][j]) j++;
        if(j == neq_popnum[eq]) emsg("cannot find");
        eqncalcansfastnum[d] = j;
        break;
    }
  }

 for(d = 0; d < nderive; d++){ if(neq_popnum[derive[d]] > 0) derivetemporal.push_back(d);}
 nderivetemporal = derivetemporal.size();
}

void Chain::eqtimelineinit()                           // Initialises the timelines use for the dependent and not dependent equations
{
  long d, loop, eq, nn, nnew, tr, jj, e, ti;
  long div, divnew, stim;
  double val, sum, dt, t, tnext;

  popnummax = 0;                                // works out the maximum number of popnum
  for(d = 0; d < ntransdepeq; d++){
    eq = transdepeq[d];
    if(neq_popnum[eq] > popnummax) popnummax = neq_popnum[eq];
  }

  ndepeq_evstack = 0;

  ndepeq_ev = 0;
  depeqdiv.resize(ntransdepeq);

  for(d = 0; d < ntransdepeq; d++){
    eq = transdepeq[d];

    //depeq_evt[ndepeq_ev] = t;
    //tdepeq_evnext[ndepeq_ev] = ndepeq_ev+1;
    //depeq_evback[ndepeq_ev] = -1;

    div = 0; t = 0; nn = -1;
    for(ti = 0; ti <= settime.size(); ti++){
      if(ti == settime.size()) tnext = tmax;
      else{ tnext = settime[ti]; if(tnext > tmax) tnext = tmax;}

      if(tnext > t){
        nnew = 0;
        for(tr = 0; tr < ntra; tr++){ 
          if(tra[tr].eq == eq && tra[tr].ci == NOTALIVE && compval[tra[tr].cf][settimecl] == ti && tra[tr].like == 1) nnew++;
        }
        if(nnew != nn || tnext == tmax){
          divnew = long(1.000000000000001*depeqdivmax*tnext/tmax);
          while(div < divnew){ depeqdiv[d].push_back(ndepeq_ev); div++;}

          nn = nnew;

          depeq_evn.push_back(nn);
          depeq_evpopnum.push_back(new double[popnummax]); for(jj = 0; jj < neq_popnum[eq]; jj++) depeq_evpopnum[ndepeq_ev][jj] = 0;
          depeq_evt.push_back(t);
          depeq_evnext.push_back(ndepeq_ev+1);
          if(t == 0) depeq_evback.push_back(-1); else depeq_evback.push_back(ndepeq_ev-1);
          val = ratecalcdep(d,depeq_evpopnum[ndepeq_ev],param); depeq_evval.push_back(val);
          ndepeq_ev++;

          Liexp += val*nn*(tnext-t);
          t = tnext;
        }
      }
      if(tnext == tmax) break;
    }

    depeq_evn.push_back(0);
    depeq_evpopnum.push_back(new double[popnummax]); for(jj = 0; jj < neq_popnum[eq]; jj++) depeq_evpopnum[ndepeq_ev][jj] = 0;
    depeq_evt.push_back(tmax);
    depeq_evnext.push_back(-1);
    depeq_evback.push_back(ndepeq_ev-1);
    depeq_evval.push_back(0);
    ndepeq_ev++;
  }

  for(e = 0; e < ndepeq_ev; e++){ ndepeq_ev_ev.push_back(0); depeq_ev_evt.push_back(vector<double> ());}

  if(depeq_evn.size() != ndepeq_ev) cout << "num wrong\n";
  if(depeq_evpopnum.size() != ndepeq_ev) cout << "num wrong\n";
  if(depeq_evt.size() != ndepeq_ev) cout << "num wrong\n";
  if(depeq_evnext.size() != ndepeq_ev) cout << "num wrong\n";
  if(depeq_evback.size() != ndepeq_ev) cout << "num wrong\n";
  if(ndepeq_ev_ev.size() != ndepeq_ev) cout << "num wrong\n";
  if(depeq_ev_evt.size() != ndepeq_ev) cout << "num wrong\n";

  for(d = 0; d < ntransnotdepeq; d++){
    eq = transnotdepeq[d];

    sum = 0;
    nn = 0;
    for(tr = 0; tr < ntra; tr++){ 
      if(tra[tr].eq == eq && tra[tr].ci == NOTALIVE && tra[tr].like == 1){
        nn++;
        if(nsettime == 0) dt = tmax;
        else{
          stim = compval[tra[tr].cf][settimecl];
          if(stim == 0) dt = settime[stim];
          else{ if(stim == nsettime) dt = tmax-settime[stim-1]; else dt = settime[stim]-settime[stim-1];}
        }
        sum += dt;
      }
    }
    val = ratecalcnotdep(eq,param);

    Liexp += val*sum;

    transnotdepeq_num.push_back(0);
    transnotdepeq_dt.push_back(sum);
    transnotdepeq_val.push_back(val);
  }
}

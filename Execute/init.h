void init()
{
  long ob, fi, c, tr, cl, eq, j, d;
  long i, ntrast;
  double sum, t;
  vector <long> eqfl;

  for(ob = 0; ob < nobs; ob++){                                                                                           // generates individual timelines
    i = obsi[ob];
    indobs[i].push_back(ob);
    indobst[i].push_back(obst[ob]);
  }

  indfixev.resize(nind);
  for(fi = 0; fi < nfixev; fi++){
    i = fixev[fi].i;
    indfixev[i].push_back(fi);
  }

  for(i = 0; i < nind; i++){
    nindobs.push_back(indobs[i].size());
    nindfixev.push_back(indfixev[i].size());
    if(nindobs[i] == 0 && nindfixev[i] == 0 && indfixtenter[i] == large && indfixtleave[i] == large) emsg("no observation");
  }

  sum = 0; logsum[0] = 0; for(i = 1; i < logsummax; i++){ sum += log(i); logsum[i] = sum;}

  compleave.resize(ncomp+1); for(tr = 0; tr < ntra; tr++) compleave[tra[tr].ci].push_back(tr);
  ncompleave.resize(ncomp+1); for(c = 0; c <= ncomp; c++) ncompleave[c] = compleave[c].size();

  compclleave.resize(ncomp+1); for(c = 0; c <= ncomp; c++) compclleave[c].resize(nclass); 
  for(cl = 0; cl < nclass; cl++){ for(tr = 0; tr < ntra; tr++) if(tra[tr].cl == cl && tra[tr].type != EXP_TR) compclleave[tra[tr].ci][cl].push_back(tr);}
  ncompclleave.resize(ncomp+1); for(c = 0; c <= ncomp; c++){ ncompclleave[c].resize(nclass); for(cl = 0; cl < nclass; cl++) ncompclleave[c][cl] = compclleave[c][cl].size();}

  trabeg = ntra;     // adds in some transition which introduce an individual at the beginning
  for(c = 0; c < ncomp; c++){
    TT tt; tt.type = BEG_TR; tt.ci = NOTALIVE; tt.cf = c; tt.cl = -1; tt.i = -1; tt.f = -1; tt.eq = -1; tt.capev = -1; tt.like = 1;
    tra.push_back(tt);
    ntra++;
  }

  traend = ntra;     // adds in some transition which introduce an individual at the beginning
  for(c = 0; c < ncomp; c++){
    TT tt; tt.type = END_TR; tt.ci = c; tt.cf = NOTALIVE; tt.cl = -1; tt.i = -1; tt.f = -1; tt.eq = -1; tt.capev = -1; tt.like = 1; 
    tra.push_back(tt);
    ntra++;
  }

  tranull = ntra;
  TT tt; tt.type = NULL_TR; tt.ci = NOTALIVE; tt.cf = NOTALIVE; tt.cl = -1; tt.i = -1; tt.f = -1; tt.eq = -1; tt.capev = -1; tt.like = 1;
  tra.push_back(tt);
  ntra++;

  for(tr = 0; tr < ntra; tr++){
    tra[tr].dc = tra[tr].cf-tra[tr].ci;
    switch(tra[tr].type){
      case FIXED_TR: case GAMMA_TR: case WEI_TR: tra[tr].nm = 1; break;
      default: tra[tr].nm = 0;
    }
    //if(tra[tr].nm == 1 && tra[tr].fix == 1) fixnmfl = 1;
    if(tra[tr].nm == 1){ cl = tra[tr].cl; for(j = 0; j < clnm.size(); j++) if(clnm[j] == cl) break; if(j == clnm.size()) clnm.push_back(cl);}
  }
  nclnm = clnm.size();
  cstnm.resize(nclass);
  tstnm.resize(nclass);

  eqfl.resize(eqnstr.size());                          // works out which equtions relate to non-Karkovian transitions
  for(eq = 0; eq < eqnstr.size(); eq++) eqfl[eq] = 0;

  for(tr = 0; tr < ntra; tr++){
    switch(tra[tr].type){
      case GAMMA_TR: case WEI_TR: eqfl[tra[tr].eq] = 1; eqfl[tra[tr].eqshape] = 1; break;
      case FIXED_TR: eqfl[tra[tr].eq] = 1; break;
    }
  }

  nnmeq = 0; for(eq = 0; eq < eqnstr.size(); eq++) if(eqfl[eq] == 1){ nmeq.push_back(eq); nnmeq++;}

  if(simon == 0){
    firstobst.resize(nind); lastobst.resize(nind);
    for(i = 0; i < nind; i++){
      t = tmax;
      if(nindobs[i] > 0){ if(indobst[i][0] < t) t = indobst[i][0];}
      if(nindfixev[i] > 0){ if(fixev[indfixev[i][0]].t < t) t = fixev[indfixev[i][0]].t;}
      if(indfixtenter[i] != large){ if(indfixtenter[i] < t) t = indfixtenter[i];}
      if(indfixtleave[i] != large){ if(indfixtleave[i] < t) t = indfixtleave[i];}
      if(t == tmax) emsg("prob z");                                      // Need to add in fixed event observations
      firstobst[i] = t;

      t = 0;
      if(nindobs[i] > 0){ if(indobst[i][nindobs[i]-1] > t) t = indobst[i][nindobs[i]-1];}
      if(nindfixev[i] > 0){ if(fixev[indfixev[i][nindfixev[i]-1]].t > t) t = fixev[indfixev[i][nindfixev[i]-1]].t;}
      if(indfixtenter[i] != large){ if(indfixtenter[i] > t) t = indfixtenter[i];}
      if(indfixtleave[i] != large){ if(indfixtleave[i] > t) t = indfixtleave[i];}
      if(t == 0) emsg("prob y");
      lastobst[i] = t;
    }
  }
}

void moveinit()
{
  long tr;

  if(movefl == 1){  // Copies all the transitions for "move" transitions which have a likelihood of zero
    moventra = ntra;
    for(tr = 0; tr < moventra; tr++){
      tra.push_back(tra[tr]);
      tra[tra.size()-1].capev = -1; tra[tra.size()-1].like = 0;
      ntra++;
    }
  }
}

void Chain::initparamsamp()
{
  long loop, loopmax = 1000, looptry, looptrymax = 1000;
  long paramset[nparam], nparamset, flag, eq, p, pr, j;

  param = new double[nparam];

  looptry = 0;
  do{
    sampill = 0;
    for(p = 0; p < nparam; p++) paramset[p] = 0;

    loop = 0; nparamset = 0;
    do{
      for(p = 0; p < nparam; p++){     // samples parameters from the prior
        if(paramset[p] == 0){
          pr = paramprior[p]; if(pr < 0){ cout << paramname[p] << " na\n"; emsg("No prior!");}
          flag = 0;                                // checks that the required other parameters are set

          switch(priortype[pr]){
            case EXPO: case GAMMA: case NORMAL: case LOGNORMAL: case BETA: case WEIBULL:
              eq = prioreq1[pr];
              for(j = 0; j < neq_param[eq]; j++) if(paramset[eq_param[eq][j]] == 0) flag = 1;
              break;
          }

          switch(priortype[pr]){
            case GAMMA: case NORMAL: case LOGNORMAL: case BETA: case WEIBULL:
              eq = prioreq2[pr];
              for(j = 0; j < neq_param[eq]; j++) if(paramset[eq_param[eq][j]] == 0) flag = 1;
              break;
          }

          if(flag == 0){
            param[p] = priorsamp(pr);
            //cout << paramname[p] << " " <<  param[p] <<" samp\n";
            paramset[p] = 1;
            nparamset++;
          }
        }
      }
      loop++; if(sampill == 1) break;
    }while(nparamset < nparam && loop < loopmax);
    if(loop == loopmax) emsg("Cound not find initial parameters");
    looptry++;
  }while(looptry < looptrymax && sampill == 1);
  if(looptry == looptrymax) emsg("Could not sample consistent parameters from prior.");

  //param[0] = 0.003; param[1] = 0.1;

  //for(p = 0; p < nparam; p++) cout << paramname[p] << " " << param[p] << "\n";
}

void Chain::initchain()
{
  long c, cl, d, k, i, p, j;
  double av;

  nindinit = 0;
  ncompinit.resize(ncompswa); for(c = 0; c < ncompswa; c++) ncompinit[c] = 0;
  probinit.resize(ncompswa); for(c = 0; c < ncompswa; c++) probinit[c] = 1.0/ncompswa;

  Lir = 0; Liexp = 0; Liinit = 0; Linm = 0; Lev_pd = 0; Lob_pd = 0; Lob_st = nobs*notobsdL; Lpop = 0; Lder = 0;

  pop.resize(npopm); popL.resize(npopm);
  for(p = 0; p < npopm; p++){ pop[p] = 0; popL[p] = popcalc(p,pop[p]); Lpop += popL[p];}

  dermpopnum.resize(nderm); derL.resize(nderm);
  for(j = 0; j < nderm; j++){
    dermpopnum[j] = new double[npopnum];
    for(k = 0; k < npopnum; k++) dermpopnum[j][k] = 0;
    derL[j] = dercalc(j,dermpopnum[j],param);
    Lder += derL[j];
  }

  chainobsinit();
  chainpartinit();
  addreminit();
  simsumcalc();

  nmeq_val.resize(eqnstr.size()); setnmeq();

  eqtimelineinit();

  pairac.resize(nclass-3); for(cl = 0; cl < nclass-3; cl++){ pairac[cl].resize(2); pairac[cl][0] = 1; pairac[cl][1] = 1;}

  ntr_param.resize(nparam); nac_param.resize(nparam); jump_param.resize(nparam);
  for(p = 0; p < nparam; p++){ ntr_param[p] = 0.1; nac_param[p] = 0; jump_param[p] = param[p]/10; if(jump_param[p] == 0) jump_param[p] = 0.1;}

  ntr_paramsamp.resize(nparam); nac_paramsamp.resize(nparam);
  for(p = 0; p < nparam; p++){ ntr_paramsamp[p] = 0.1; nac_paramsamp[p] = 0;}

  ntr_paramnorm.resize(nparam); nac_paramnorm.resize(nparam);
  for(p = 0; p < nparam; p++){ ntr_paramnorm[p] = 0.1; nac_paramnorm[p] = 0;}

  ntr_smooth.resize(nsmooth); nac_smooth.resize(nsmooth); jump_smooth.resize(nsmooth);
  for(p = 0; p < nsmooth; p++){
    ntr_smooth[p] = 0.1; nac_smooth[p] = 0;
    av = 0; for(j = 0; j < nsmoothparam[p]; j++) av += param[smoothparam[p][j]];
    jump_smooth[p] = av/(10*nsmoothparam[p]); if(jump_smooth[p] == 0) jump_smooth[p] = 0.1;
  }

  ntr_part.resize(nind); nac_part.resize(nind); nfa_part.resize(nind); nde_part.resize(nind);
  for(i = 0; i < nind; i++){ 
    ntr_part[i].resize(clall+1); nac_part[i].resize(clall+1); nfa_part[i].resize(clall+1); nde_part[i].resize(clall+1);
    for(cl = 0; cl <= clall; cl++){ ntr_part[i][cl] = 0.1; nac_part[i][cl] = 0; nfa_part[i][cl] = 0; nde_part[i][cl] = 0;}
  }

  ntr_indsim.resize(nind); nac_indsim.resize(nind); nfa_indsim.resize(nind);
  for(i = 0; i < nind; i++){ ntr_indsim[i] = 0.1; nac_indsim[i] = 0; nfa_indsim[i] = 0;}
  ntr_indsimuo = 0.1; nac_indsimuo = 0; nfa_indsimuo = 0;

  ntr_move.resize(nind); nac_move.resize(nind); jump_move.resize(nind);
  for(i = 0; i < nind; i++){ 
    for(cl = 0; cl < clall; cl++){ ntr_move[i].push_back(0.1); nac_move[i].push_back(0.1);}
    jump_move[i].resize(ntrans); for(j = 0; j < ntrans; j++) jump_move[i][j] = tmax/10;
  }

  ntr_sing.resize(nind); nac_sing.resize(nind); nfa_sing.resize(nind);
  for(i = 0; i < nind; i++){ 
    for(cl = 0; cl < clall; cl++){ ntr_sing[i].push_back(0.1); nac_sing[i].push_back(0.1); nfa_sing[i].push_back(0.1);}
  }

  ntr_twothree.resize(nind); nac_twothree.resize(nind); nfa_twothree.resize(nind);
  for(i = 0; i < nind; i++){ 
    for(cl = 0; cl < clall; cl++){ ntr_twothree[i].push_back(0.1); nac_twothree[i].push_back(0.1); nfa_twothree[i].push_back(0.1);}
  }

  ntr_pair.resize(nind); nac_pair.resize(nind); nfa_pair.resize(nind);
  for(i = 0; i < nind; i++){
    for(cl = 0; cl < clall; cl++){ntr_pair[i].push_back(0.1); nac_pair[i].push_back(0.1); nfa_pair[i].push_back(0.1);}
  }

  ntr_gap.resize(nind); nac_gap.resize(nind); nfa_gap.resize(nind);
  for(i = 0; i < nind; i++){ 
    for(cl = 0; cl < clall; cl++){ntr_gap[i].push_back(0.1); nac_gap[i].push_back(0.1); nfa_gap[i].push_back(0.1);}
  }

  ntr_tbirth.resize(nind); nac_tbirth.resize(nind); jump_tbirth.resize(nind);                      // Diagnostics for proposals on life of observed individuals
  for(i = 0; i < nind; i++){ ntr_tbirth[i] = 0.1; nac_tbirth[i] = 0; jump_tbirth[i] = tmax/10;}

  ntr_tbirthent.resize(nind); nac_tbirthent.resize(nind); jump_tbirthent.resize(nind);
  for(i = 0; i < nind; i++){ ntr_tbirthent[i] = 0.1; nac_tbirthent[i] = 0; jump_tbirthent[i] = tmax/10;}

  ntr_tent.resize(nind); nac_tent.resize(nind); jump_tent.resize(nind);
  for(i = 0; i < nind; i++){ ntr_tent[i] = 0.1; nac_tent[i] = 0; jump_tent[i] = tmax/10;}

  ntr_tlea.resize(nind); nac_tlea.resize(nind); jump_tlea.resize(nind);
  for(i = 0; i < nind; i++){ ntr_tlea[i] = 0.1; nac_tlea[i] = 0; jump_tlea[i] = tmax/10;}

  ntr_entswitch.resize(nind); nac_entswitch.resize(nind);
  for(i = 0; i < nind; i++){ ntr_entswitch[i] = 0.1; nac_entswitch[i] = 0;}

  ntr_birthentswitch.resize(nind); nac_birthentswitch.resize(nind);
  for(i = 0; i < nind; i++){ ntr_birthentswitch[i] = 0.1; nac_birthentswitch[i] = 0;}

  ntr_leaswitch.resize(nind); nac_leaswitch.resize(nind);
  for(i = 0; i < nind; i++){ ntr_leaswitch[i] = 0.1; nac_leaswitch[i] = 0;}

  ntr_tbirthuo = 0.1; nac_tbirthuo = 0; jump_tbirthuo = tmax/10;                      // Diagnostics for proposals on life of uoobserved individuals
  ntr_tbirthentuo = 0.1; nac_tbirthentuo = 0; jump_tbirthentuo = tmax/10;
  ntr_tentuo = 0.1; nac_tentuo = 0; jump_tentuo = tmax/10;
  ntr_tleauo = 0.1; nac_tleauo = 0; jump_tleauo = tmax/10;
  ntr_entswitchuo = 0.1; nac_entswitchuo = 0;
  ntr_birthentswitchuo = 0.1; nac_birthentswitchuo = 0;
  ntr_leaswitchuo = 0.1; nac_leaswitchuo = 0;

  ntr_add = 0.1; nac_add = 0; nfa_add = 0;
  ntr_rem = 0.1; nac_rem = 0;

  for(i = 0; i < 30; i++){ timeprop[i] = 0; ntimeprop[i] = 0;}
}

long Chain::start()                                                           // Works out the starting event sequence
{
  long c, cc, tr, j, k, ob, per, pernew, loop, i, ctop, ag, agpos[nage+1], ob1, ob2, ca, ca1, ca2, numobs, numnotobs, nlifet;
  double t, tt, prob, probsum, Ltoti, Ltotf, val, lifet, probmax;//z;
  vector <double> Lst;
  vector <short> agpostlist;
  vector < vector <EV> > posev;
  vector <double> posevprob;
  vector <double> posevprobsum;
  vector <double> posevtbirth;

  timeprop[INIT] -= clock();

  // Makes an estimate for capture probability and mortality
  numobs = 0; numnotobs = 0; lifet = 0; nlifet = 0;
  for(i = 0; i < nind; i++){
    if(nindobs[i] >= 2){
      lifet += obst[indobs[i][nindobs[i]-1]] - obst[indobs[i][0]];
      nlifet++;

      numobs += nindobs[i]-2;
      for(j = 0; j < nindobs[i]-1; j++){
        ob1 = indobs[i][j]; ob2 = indobs[i][j+1];
        ca1 = obscap[ob1]; ca2 = obscap[ob2];
        if(ca1 >= 0 && ca2 >= 0){
          for(ca = 0; ca < ncap; ca++){
            t = capt[ca];
            if(t > obst[ob1] && t < obst[ob2]){
              for(c = 0; c < ncomp; c++){
                if(capprobeqn[ca][c] >= 0 && 
                  ((obsprobeqn[ob1][c] >= 0 && capprobeqn[ca1][c] >= 0) ||
                   (obsprobeqn[ob2][c] >= 0 && capprobeqn[ca2][c] >= 0))) break;
              }
              if(c < ncomp){
                numnotobs++;
                while(ca < ncap-1 && capt[ca+1] == t) ca++;
              }
            }
          }
        }
      }
    }
  }
  if(numobs+numnotobs > 0) capprobest = double(numobs)/(numobs+numnotobs); else capprobest = 0.95;

  if(nlifet  > 0) muest = 1.0/(lifet/nlifet); else muest = 10.0/tmax;

  nindtot = 0;

  samp = -1; per = 60;
  for(i = 0; i < nind; i++){
    //cout << i << " " << indid[i] << " i\n";
    //if(i < 5) cout << i << " "<< Liexp << " Li\n";
    pernew = long(60+(30.0*i)/nind); if(pernew != per){ cout << "3|" << pernew << "|\n"; cout.flush(); per = pernew;}

    addemptyind(0);

    capwarn = "";

    probsum = 0; posev.clear(); posevprob.clear(); posevprobsum.clear(); posevtbirth.clear();
    for(loopstart = 0; loopstart < loopstartmax; loopstart++){
      if(loopstart%4 == 0){ uocapmin = -large; uocapmax = large;}

      if(sourcefl == 0) tent = 0;
      else{
        if(indfixtenter[i] != large) tent = indfixtenter[i];
        else{
          if(tbirthfl == 0){
            if(ran() < 0.1 || notimerange == 1) tent = 0;
            else{
              if(ran() < 0.1) tent = ran()*firstobst[i];
              else tent = firstobst[i] + log(ran())/muest;
            }
          }
          else{
            if(notimerange == 1) tent = 0;
            else{
              tent = firstobst[i] - ran()*age[nage];
              if(tent < 0) tent = 0;
            }
          }
        }
      }
      if(tent < 0) tent = 0;
      if(tent < uocapmin){ tent = uocapmin + ran()*(firstobst[i]-uocapmin);}

      if(tbirthfl == 0) indtbirth[i] = tent; 
      else{
        ctop = getctop(tent,tent);  // Works out the possible age classifications

        for(ag = 0; ag <= nage; ag++){ agpos[ag] = 0; if(tent - age[ag] < 0) agpos[ag] = 1;}

        for(c = 0; c < ncomps; c++){
          for(ag = 0; ag <= nage; ag++){
            tr = compiftra[NOTALIVE][c+ctop+ag*classmult[agecl]];
            if(tr >= 0) agpos[ag] = 1;
          }
        }
        agpostlist.clear();
        for(ag = 0; ag <= nage; ag++){ if(agpos[ag] == 1) agpostlist.push_back(ag);}

        if(agpostlist.size() == 0) indtbirth[i] =  -ran()*age[ag];
        else{
          ag = agpostlist[long(ran()*agpostlist.size())];
          if(ag == 0) indtbirth[i] = tent - ran()*age[ag];
          else indtbirth[i] = tent - age[ag-1] - ran()*(age[ag]-age[ag-1]);
        }
      }

      if(sinkfl == 0) tlea = tmax;
      else{
        if(indfixtleave[i] != large) tlea = indfixtleave[i];
        else{
          if(ran() < 0.1 || notimerange == 1) tlea = tmax;
          else{
            if(ran() < 0.1) tlea = lastobst[i] + ran()*(tmax-lastobst[i]);
            else{
              tlea = lastobst[i] - log(ran())/muest;
              if(tlea > tmax) tlea = tmax;
            }
          }
        }
      }

      if(tlea > uocapmax){ tlea = lastobst[i] + ran()*(uocapmax-lastobst[i]);}

      prob = startseq(i);
      if(prob > 0){
        posev.push_back(evnew);
        probsum += prob;
        posevprob.push_back(prob);
        posevprobsum.push_back(probsum);
        posevtbirth.push_back(indtbirth[i]);
      }
      if(posev.size() > 10 || probsum > 3) break;
    }

    if(loopstart == loopstartmax){
      if(capwarn != "" && capwarn != "not") emsg2(capwarn);
      else{
        stringstream ss; ss << "Could not get initial condition for individual '" << indid[i] << "'.";
        emsg2(ss.str());
      }
    }
    else{
      //cout << probsum << " " << posev.size()  << "U\n";
      //z = ran()*probsum; j = 0; while(j < posev.size() && z > posevprobsum[j]) j++; if(j == posev.size()) emsg("posev");
      probmax =0; for(k = 0; k < posev.size(); k++){ if(posevprob[k] > probmax){ probmax = posevprob[k]; j = k;}}
      evnew = posev[j];
      indtbirth[i] = posevtbirth[j];
      indcha(i);
    }
    //cout << loopstart  << " loopst\n";
    if(warning.size() > 10) break;
  }

  if(warning.size() > 0){
    cout << "e|";
    for(j = 0; j < warning.size(); j++) cout << warning[j] << "*";
    cout << "\n";
    exit (EXIT_FAILURE);
  }

  Lpri = priorcalc();

  //long p; for(p = 0; p < nparam; p++) cout << paramname[p] << " " << param[p] << ", "; cout << " param\n";

  //long pr; for(pr = 0; pr < nprior; pr++) cout << paramname[priorparam[pr]] << " " << priorprob(pr) << " pr\n";
  //emsg("P");

  //ch[0]->check(-2);
  //emsg("P");

  if(simon == 0){
    long p; pvel.resize(nparam); for(p = 0; p < nparam; p++) pvel[p] = 0.001;

    samp = 0;
    do{                                              // Finds a reasonable parameter set for the initial event sequence
      Lst.push_back(L());
      param_prop();
      param_gradprop(); 
      //for(p = 0; p < nparam; p++) cout << paramname[p] << " " << param[p] << ", "; cout << " param\n";
      //cout << L() << " " << Lst[Lst.size()-1] << " "<<   L() -Lst[Lst.size()-1]  << " " << Linm << " " << Lst.size() << "Pa\n";

      if(Lst.size() > 20){
        if(Lst.size() > 10000) emsg("Initial like problem");
        //Lout();
        
        //cout <<  Lpri << " " <<  (Liinit + Lir - Liexp + Linm)*invTLi << " " << Lob_pd*invT_pd + Lob_st*invT_st + Lpop*invT_pop + Lder*invT_der << " h\n";
        //cout << L()-tiny  << " " << Lst[Lst.size()-20] << " yy\n";
        if(L()-tiny <= Lst[Lst.size()-20]) break;
      }
    }while(1 == 1);
  }
  for(loop = 0; loop < 50; loop++) param_prop();
  timeprop[INIT] += clock();

  return 0;
}

void Chain::addemptyind(double tbirth)                    // Adds a new empty individual (i.e. it is in the notalive category)
{
  indtbirth.push_back(tbirth);
  nindev.push_back(2);
  evnew.clear();
  EV evbeg; evbeg.t = 0; evbeg.tr = tranull; evnew.push_back(evbeg);
  EV evend; evend.t = tmax; evend.tr = tranull; evnew.push_back(evend);
  indev.push_back(evnew);

  nindtot++;
  if(nindtot > indmax){
    stringstream ss; ss << "Maximum number of " << indmax << " individuals exceeded.";
    emsg(ss.str());
  }
}

void Chain::rememptyind(long i)                    // Removes an individual
{
  nindtot--;
  if(i < nindtot){
    nindev[i] = nindev[nindtot];
    indtbirth[i] = indtbirth[nindtot];
    indev[i] = indev[nindtot];
  }
  nindev.pop_back();
  indev.pop_back();
  indtbirth.pop_back();
}

double getlasttime(long cl, long e, vector<EV> &ev)   // Gets the last time a transition occured in a given classification
{
  long ee;
  long tr;

  for(ee = e; ee > 0; ee--){ 
    tr = ev[e].tr; if(tra[tr].cl == cl && tra[tr].type >= 0) break;
  }
  return ev[e].t;
}


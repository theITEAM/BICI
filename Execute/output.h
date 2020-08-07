vector <long> traoutref;
vector <string> traoutname;
long ntraout;
vector <long> ntraev;
long evtot, nindtotst;

void traceinit()                                                       // Initialises the trace plot
{
  long p, pp, c, cl, i, f, j, tr, d;

  ntraout = 0;
  traoutref.resize(ntra);

  for(tr = 0; tr < ntra; tr++){
    traoutref[tr] = -1;

    cl = tra[tr].cl;
    if(cl >= 0){
      i = tra[tr].i; f = tra[tr].f;
      if(i >= 0 && i < nclassval[cl] && f >=0 && f < nclassval[cl]){
        stringstream ss; ss << classval[cl][i] << " → " << classval[cl][f];
        j = 0; while(j < ntraout && traoutname[j] != ss.str()) j++;
        if(j == ntraout){ traoutname.push_back(ss.str()); ntraout++;}
        traoutref[tr] = j;
      }
    }
  }
  ntraev.resize(ntraout);

  if(noout == 1) return;

  if(simon == 1) cout << "c|"; else cout << "1|";

  for(pp = 0; pp < nparam; pp++){
    p = paramorder[pp];
    cout << paramclass[p] << "|" << paramname[p] << "|Transition rate parameter " << paramname[p] << "|";
  }

  if(simon == 0){
    for(c = 0; c < ncompswa; c++){
      cout << "Init. Prob.|";
      for(cl = 0; cl < nclass-2; cl++){ if(!(cl == agecl && nclassval[agecl] == 1)){ if(cl != 0)  cout << ","; cout << classval[cl][compval[c][cl]];}}
      cout << "|";
      cout << "Probability of initially being in ";
      for(cl = 0; cl < nclass-2; cl++){ if(!(cl == agecl && nclassval[agecl] == 1)){ if(cl != 0)  cout << ","; cout << classval[cl][compval[c][cl]];}}
      cout << "|";
    }
  }

  for(j = 0; j < ntraout; j++) cout << "Trans.|" << traoutname[j] << "|Number of " <<  traoutname[j] << " transitions|";

  for(d = 0; d < nderive; d++){
    if(neq_popnum[derive[d]] == 0) cout << "Der.|" << derivename[d] << "|Derived quantity: " << derivename[d] << "|";
  }

  if(simon == 0){
    if(checkon == 1){
      cout << "Misc.|Lev|Event likelihood|";
      cout << "Misc.|Lnm|Non Markovian event likelihood|";
      cout << "Misc.|Lob|Observation probability|";
      cout << "Misc.|Lcap|Capture probability|";
      cout << "Misc.|Levcap|Event capture probability|";
      cout << "Misc.|Linit|Initial condition probability|";
      cout << "Misc.|Lpop|Population probability|";
      cout << "Misc.|Lpop|Derived probability|";
      cout << "Misc.|Pr|Prior probability|";
      cout << "Misc.|Nerr|The number of observation|";
    }
    else{
      cout << "Misc.|log(π(θ,ξ¦D))|The log of the posterior probability|";
      cout << "Misc.|log(π(D¦ξ,θ))|The log of the observation probability|";
      cout << "Misc.|log(π(ξ¦θ))|The log of latent process likelihood|";
      cout << "Misc.|log(π(θ))|The log of the prior probability|";
    }
  }
  cout << "Misc.|# Ind.|Total number of individuals|";
  cout << "Misc.|# Transitions|Total number transitions|";
  cout << "\n";
  cout.flush();
}

void Chain::traceplot()                                                 // Plots traces for the parameter values
{
  long p, pp, cl, k, d, c, e;
  long i, j, tr, evtotst, nindtotst;
  double Ltot, Lobtot;

  if(noout == 1) return;

  timeprop[TRACEPLOT] -= clock();

  evtotst = 0; 
  if(simon == 1){ nindtotst = nindtot_sim; for(i = 0; i < nindtot_sim; i++) evtotst += nindev_sim[i];}
  else{ nindtotst = nindtot; for(i = 0; i < nindtot; i++) evtotst += nindev[i];}

  if(simon == 1) cout << "d|"; else cout << "0|";

  for(pp = 0; pp < nparam; pp++){ p = paramorder[pp]; cout << param[p] << "|";}

  if(simon == 0){
    for(c = 0; c < ncompswa; c++) cout << probinit[c] << "|";
  }

  for(j = 0; j < ntraout; j++) cout << ntraev[j] << "|";

  for(d = 0; d < nderive; d++){ if(neq_popnum[derive[d]] == 0) cout << calculatenotdep(derive[d],param) << "|";}

  if(simon == 0){
    Ltot = Lir-Liexp + Linm + Liinit;
    Lobtot = Lob_st+Lob_pd +Lev_pd + Lpop+Lder;
    if(checkon == 1) cout << Lir-Liexp << "|" << Linm << "|" << Lob_st << "|" << Lob_pd << "|" << Lev_pd << "|" << Liinit << "|"<<  Lpop << "|"<< Lder << "|" << Lpri << "|" << nwrong << "|";
    else cout << Ltot+Lobtot+Lpri << "|" << Lobtot << "|" << Ltot << "|" << Lpri << "|";
  }
  //cout << nindtotst << "|" << evtot << "|" ;
  cout << nindtotst << "|" << evtotst << "|" ;

  //if(nindtot == 0) cout << "0"; else cout << indtbirth[0];
  //cout << "|";
  cout << "\n";
  cout.flush();

  timeprop[TRACEPLOT] += clock();
}

void oe(string st, vector<EV> ev)
{
  long e, ci, cf;
  long tr;
  double t;

  timeprop[EVENTPLOT] -= clock();

  if(ev.size() == 0) cout << "NA ";
  else{
    if(1 == 1){
      for(e = 0; e < ev.size(); e++){
        tr = ev[e].tr; t = ev[e].t;
        //outcomp(tra[tr].ci);
        if(tra[tr].ci == NOTALIVE) cout << "NA";
        cout << " -> "; outcomp(tra[tr].cf); cout << " " << t << ",       "; 
      }
    }
    else{
      outcomp(tra[ev[0].tr].ci); cout << ":";
      for(e = 0; e < ev.size(); e++){
        tr = ev[e].tr; t = ev[e].t;
        if(tra[tr].type < 0){
          cout << "END";
        }
        else{
          ci = tra[tr].ci; cf = tra[tr].cf;
          cout << " -> ";
          if(ci == NOTALIVE || cf == NOTALIVE) outcomp(cf);
          else cout << classval[tra[tr].cl][tra[tr].f];
        }
        cout << " " << t << ",  ";
      }
    }
  }
  cout << st << "\n";

  timeprop[EVENTPLOT] += clock();
}

void oprop()
{
  cout << "\n";
  oe("BEF",evrev);
  oe("AFT",evnew);
}

void outputmodel()                                     // Outputs the model (used in debugging)
{
  long cl, tr, d, p, j, c, eq, ob, cap, f;

  ofstream model("GP/model.txt");

  for(p = 0; p < npopm; p++){
    model << p << " t:" <<  popmt[p] << "   c:";
    for(j = 0; j < popmcomp[p].size(); j++) model << popmcomp[p][j] << ",";
    model << "\n";
  }

  model << "Trasitions:\n";
  for(tr = 0; tr < ntra; tr++){
    cl = tra[tr].cl;
    if(cl < 0){
      c = tra[tr].ci; if(c == NOTALIVE) model << "NA"; else{ for(cl = 0; cl < nclass; cl++){ model << classval[cl][compval[c][cl]]; if(cl < nclass-1) model << ",";}}
      model << " -> ";
      c = tra[tr].cf; if(c == NOTALIVE) model << "NA"; else{ for(cl = 0; cl < nclass; cl++){ model << classval[cl][compval[c][cl]]; if(cl < nclass-1) model << ",";}}
    }
    else model << classname[cl] << ": " << classval[cl][tra[tr].i] << " -> " << classval[cl][tra[tr].f];
    model << "     ";
    eq = tra[tr].eq; if(eq >= 0) model << eqnstr[tra[tr].eq];
    model << "  Capev:" << tra[tr].capev << "   Like:" << tra[tr].like << "\n";
  }

  for(d = 0; d < ntransdepeq; d++) model << eqnstr[transdepeq[d]] << " depeqn\n";
  for(d = 0; d < ntransnotdepeq; d++) model << eqnstr[transnotdepeq[d]] << " notdepeqn\n";

  for(p = 0; p < paramname.size(); p++){
    model << paramname[p] << "\n";
  }

  for(p = 0; p < popnumname.size(); p++){
    model << popnumname[p] << ": ";
    //for(j = 0; j < popnumterm[p].size(); j++) model << popnumterm[p][j] << " " << popnumtermweight[p][j] << ",  ";
    model << "\n"; 
  }

  for(cl = 0; cl < nclass; cl++){
    model << classname[cl] << ":\n";
    for(j = 0; j < nclassval[cl]; j++){
      model << classval[cl][j] << "\n";
    }
  }

  for(cap = 0; cap < ncap; cap++){
    model << capt[cap] << " " << capname[cap] << "; ";
    model << " cap\n\n";
  }
  model << "\n";

  for(ob = 0; ob < nobs; ob++){
    model << obst[ob] << " " << indid[obsi[ob]] << " ";
    if(obscap[ob] >= 0)  model << capname[obscap[ob]] << "; ";
    //for(c = 0; c < ncomp; c++) if(obsprobeqn[ob][c] != -1) model << compname[c] << ",";
    model << " obs\n\n";
  }

  for(eq = 0; eq < eqnstr.size(); eq++){
    model << "Eqn:" << eqnstr[eq] << "\n";
  }

  model << "Prior:\n";
  for(p = 0; p < priortype.size(); p++){
    model << p << " ";
    switch(priortype[p]){
      case FLAT: model << "FLAT" << priorminval[p] << " " <<  priormaxval[p] << "\n"; break;
      case GAMMA: model << "GAMMA" << eqnstr[prioreq1[p]] << " " << eqnstr[prioreq2[p]] << "\n"; break;
      case NORMAL: model << "NORMAL" << eqnstr[prioreq1[p]] << " " << eqnstr[prioreq2[p]] << "\n"; break;
      case LOGNORMAL: model << "LOGNORMAL" << eqnstr[prioreq1[p]] << " " << eqnstr[prioreq2[p]] << "\n"; break;
      case EXPO: model << "EXPO" << eqnstr[prioreq1[p]] << "\n"; break;
      case BETA: model << "BETA" << eqnstr[prioreq1[p]] << " " << eqnstr[prioreq2[p]] << "\n"; break;
      case WEIBULL: model << "WEIBULL" << eqnstr[prioreq1[p]] << " " << eqnstr[prioreq2[p]] << "\n"; break;

      default: model << "error prior type\n"; break;
    }
  }

  model << "\n";
  for(f = 0; f <= ncapevtrange; f++) model << capevtrange[f] << ","; model << "capevtrange\n";
  model << "\n";

  model << "fixfl: " << fixfl << "\n";
  model << "sinkfl: " << sinkfl << "\n";
  model << "sourcefl: " << sourcefl << "\n";
  model << "gammafl: " << gammafl << "\n";
  model << "nmfl: " << nmfl << "\n";
  model << "capfl: " << capfl << "\n";
  model << "popfl: " << popfl << "\n";
  model << "tbirthfl: " << tbirthfl << "\n";
  model << "pairfl: " << pairfl << "\n";
  model << "singeventfl: " << singeventfl << "\n";
  model << "twothreefl: " << twothreefl << "\n";
  model << "gapfl: " << gapfl << "\n";
 
  model << "\n";
}

void outcomp(long c)                // outputs the state of a compartment
{
  long cl;
  if(c == NOTALIVE) cout << "NA";
  else{
    cout << compname[c];
   // for(cl = 0; cl < nclass; cl++){ cout << classval[cl][compval[c][cl]]; if(cl < nclass-1) cout << ",";}
  }
}

void Chain::eventplot()  // outputs an events sample
{
  long cl, k, loop, tr, cf, i, j, e;
  double fac;
  vector <long> list;
  vector <double> listt;

  for(loop = 0; loop < 2; loop++){
    if(simon == 1){ sim(0); cout << "b|";}
    else{
      switch(loop){
        case 0: sim(tmax); cout << "5|"; break;
        case 1: sim(0); cout << "a|"; break;
      }
    }

    fac = (tmaxactual-tminactual)/tmax2; 
    if(fac > 1.01 || fac < 0.99){
      cout << tminactual << " "<< tmaxactual << " " << tmax2 << " "<< fac << "\n"; emsg("facprob",fac,tmax,tmax2);
    }

    if(noout == 0){
      cout << nindtot_sim << "|";
      for(i = 0; i < nindtot_sim; i++){
        list.clear(); listt.clear();
        for(k = 0; k < nindev_sim[i]; k++){
          tr = indev_sim[i][k].tr;
          if(tra[tr].type >= 0){
            cf = tra[tr].cf;
            if(cf == ncomp) cf = classmult[capevcl]; else cf = cf%classmult[capevcl];
            list.push_back(cf); listt.push_back(tminactual + fac*indev_sim[i][k].t);
          }
        }
        cout << indtbirth_sim[i] << "|";
        cout << list.size() << "|"; for(k = 0; k < list.size(); k++) cout << list[k] << "|" << listt[k] << "|";
      }
      cout << "\n";
    }

    if(noout == 2){
      for(i = 0; i < nindtot_sim; i++){
        cout << "ind" << i << ": ";
        for(k = 0; k < nindev_sim[i]; k++){
          tr = indev_sim[i][k].tr;
          cf = tra[tr].cf;  cout.flush();
          if(cf == ncomp) cout << classmult[capevcl];
          else cout << tra[tr].cf%classmult[capevcl];
          cout << " " << tminactual + fac*indev_sim[i][k].t << ",  ";
        }
        cout << "\n";
      }
    }

    if(noout == 0 && loop == 0 && nderive > 0){ deriveplotcalc(); derivepl();}

    if(loop == 0){
      nindtotst = nindtot_sim;

      evtot = 0;
      for(j = 0; j < ntraout; j++) ntraev[j] = 0;

      for(i = 0; i < nindtot_sim; i++){
        for(e = 0; e < nindev_sim[i]; e++){
          tr = indev_sim[i][e].tr;
          if(traoutref[tr] >= 0) ntraev[traoutref[tr]]++;
          if(indev_sim[i][e].t != 0 && indev_sim[i][e].t != tmax && tra[tr].cl < nclass-3) evtot++;
        }
      } 
    }
    if(simon == 1) break;
  }
  cout.flush();
}

void Chain::diagnosticschain()        // Outputs the success of MCMC proposals
{
  long p, e, cl;
  long i;
  double f, fmin, fmax, fav;

  if(noout == 1) return;

  cout << "6|";
  cout << "Random walk parameter samples:|";
  for(p = 0; p < nparam; p++) cout << paramname[p] << " ac: " << nac_param[p]/ntr_param[p] << "  jump: " << jump_param[p] << "|";
  //for(p = 0; p < nparam; p++) cout << paramname[p] << " ac samp: " << nac_paramsamp[p]/ntr_paramsamp[p] << "|";

  cout << "|Normal fitted parameter samples :|";
  for(p = 0; p < nparam; p++) cout << paramname[p] << " ac norm: " << nac_paramnorm[p]/ntr_paramnorm[p] << "|";

  cout << "|Random walk global smoothing proposals :|";
  for(p = 0; p < nsmooth; p++) cout << smoothname[p] << " ac: " << nac_smooth[p]/ntr_smooth[p] << "  jump: " << jump_smooth[p] << "|";

  for(cl = 0; cl <= clall; cl++){
    cout << "|Particle proposals in ";
    if(cl == clall) cout << "all classifications";
    else cout << classname[cl];
    cout << ":|";

    fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nac_part[i][cl]/ntr_part[i][cl]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
    cout << "part ac:" << fav << ", " << fmin << "-" << fmax << "      |";

    fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nfa_part[i][cl]/ntr_part[i][cl]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
    cout << "part fa:" << fav << ", " << fmin << "-" << fmax << "      |";

    fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nde_part[i][cl]/ntr_part[i][cl]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
    cout << "part de:" << fav << ", " << fmin << "-" << fmax << "      |";

    fmin = npartmax; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nindpart[i][cl]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
    cout << "nindpart:" << fav << ", " << fmin << "-" << fmax << "      |";
  }

  cout << "|Simulate proposals (used on unobserved):|";
  fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nac_indsim[i]/ntr_indsim[i]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
  cout << "indsim ac:" << fav << ", " << fmin << "-" << fmax << "      |";

  fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nfa_indsim[i]/ntr_indsim[i]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
  cout << "indsim fa:" << fav << ", " << fmin << "-" << fmax << "      |";

  cout << "indsimuo ac:" <<  nac_indsimuo/ntr_indsimuo << "      |";
  cout << "indsimuo fa:" <<  nfa_indsimuo/ntr_indsimuo << "      |";

  for(cl = 0; cl < clall; cl++){
    cout << "|Local changes in " <<  classname[cl] << ":|";
    fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nac_move[i][cl]/ntr_move[i][cl]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
    cout << "move ac:" << fav << ", " << fmin << "-" << fmax << "      |";

    fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nac_sing[i][cl]/ntr_sing[i][cl]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
    cout << "sing ac:" << fav << ", " << fmin << "-" << fmax << "      |";

    fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nfa_sing[i][cl]/ntr_sing[i][cl]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
    cout << "sing fa:" << fav << ", " << fmin << "-" << fmax << "      |";

    fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nac_twothree[i][cl]/ntr_twothree[i][cl]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
    cout << "twothree ac:" << fav << ", " << fmin << "-" << fmax << "      |";

    fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nfa_twothree[i][cl]/ntr_twothree[i][cl]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
    cout << "twothree fa:" << fav << ", " << fmin << "-" << fmax << "      |";

    fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nac_pair[i][cl]/ntr_pair[i][cl]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
    cout << "pair ac:" << fav << ", " << fmin << "-" << fmax << "      |";

    fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nfa_pair[i][cl]/ntr_pair[i][cl]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
    cout << "pair fa:" << fav << ", " << fmin << "-" << fmax << "      |";

    fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nac_gap[i][cl]/ntr_gap[i][cl]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
    cout << "gap ac:" << fav << ", " << fmin << "-" << fmax << "      |";

    fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nfa_gap[i][cl]/ntr_gap[i][cl]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
    cout << "gap fa:" << fav << ", " << fmin << "-" << fmax << "      |";
  }

  // life obsevred
  cout << "|Life span proposals:|";
  fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nac_tbirth[i]/ntr_tbirth[i]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
  cout << "tbirth ac:" << fav << ", " << fmin << "-" << fmax << "      |";

  fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nac_tbirthent[i]/ntr_tbirthent[i]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
  cout << "tbirthent ac:" << fav << ", " << fmin << "-" << fmax << "      |";

  fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nac_tent[i]/ntr_tent[i]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
  cout << "tent ac:" << fav << ", " << fmin << "-" << fmax << "      |";

  fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nac_tlea[i]/ntr_tlea[i]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
  cout << "tlea ac:" << fav << ", " << fmin << "-" << fmax << "      |";

  fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nac_entswitch[i]/ntr_entswitch[i]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
  cout << "entswitch ac:" << fav << ", " << fmin << "-" << fmax << "      |";

  fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nac_birthentswitch[i]/ntr_birthentswitch[i]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
  cout << "birthentswitch ac:" << fav << ", " << fmin << "-" << fmax << "      |";

  fmin = 1; fmax = 0; fav = 0; for(i = 0; i < nind; i++){ f = nac_leaswitch[i]/ntr_leaswitch[i]; if (f > fmax) fmax = f; if(f < fmin) fmin = f; fav += f/nind;}
  cout << "leaswitch ac:" << fav << ", " << fmin << "-" << fmax << "      |";

  // life unobserved
  cout << "tbirthuo ac:" << nac_tbirthuo/ntr_tbirthuo << "      |";
  cout << "tbirthentuo ac:" << nac_tbirthentuo/ntr_tbirthentuo << "      |";
  cout << "tentuo ac:" << nac_tentuo/ntr_tentuo << "      |";
  cout << "tleauo ac:" <<  nac_tleauo/ntr_tleauo << "      |";
  cout << "entswitchuo ac:" << nac_entswitchuo/ntr_entswitchuo << "      |";
  cout << "birthentswitchuo ac:" <<  nac_birthentswitchuo/ntr_birthentswitchuo << "      |";
  cout << "leaswitchuo ac:" << nac_leaswitchuo/ntr_leaswitchuo << "      |"; 
  cout << "add ac:" << nac_add/ntr_add << "    |";
  cout << "rem ac:" << nac_rem/ntr_rem << "    |";
  cout << "add fa:" << nfa_add/ntr_add << "    |";

  cout << "|Probability of proposals:|";
  cout << "part:" << part_prob << ", induo:" << induo_prob << ", move:" << move_prob << ", sing:" << sing_prob;
  cout << ", twothree:" << twothree_prob << ", pair:" << pair_prob << ", gap:" << gap_prob << ", life:" << life_prob << ", addrem:" << addrem_prob <<"   |";

  cout << "|Time of proposals:|";
  //for(p = 0; p < nparam; p++) cout << paramname[p] << " " << long(100*timeparam[p]/totaltime) << "%|";
  //for(p = 0; p < nparam; p++) cout << paramname[p] << " " << long(100*timegradparam[p]/totaltime) << "% grad|";
  
  cout << "Paramameters:" << long(100*timeprop[PARAM_PROP]/totaltime) << "%|";
  cout << "Particle:" << long(100*timeprop[PART_PROP]/totaltime) << "%|";
  cout << "Simualtion of unobserved:" << long(100*timeprop[INDSIMUO_PROP]/totaltime) << "%|";
  cout << "Moving events:" << long(100*timeprop[MOVE_PROP]/totaltime) << "%|";
  cout << "Lifespan proposals:" << long(100*timeprop[LIFE_PROP]/totaltime) << "%|";
  cout << "Adding/removing indivduals:" << long(100*timeprop[ADDREM_PROP]/totaltime) << "%|";
  cout << "Adding/removing single event:" << long(100*timeprop[SINGEVENT_PROP]/totaltime) << "%|";
  cout << "Adding/removing pairs of events:" << long(100*timeprop[PAIR_PROP]/totaltime) << "%|";
  cout << "Local changes to events:" << long(100*timeprop[GAP_PROP]/totaltime) << "%|";
  cout << "Diagnostic checks:" << long(100*timeprop[CHECK]/totaltime) << "%|";
/*
  cout << "Time indchange:" << long(100*timeprop[INDCHANGE]/totaltime) << "%|";
  cout << "Time trace:" << long(100*timeprop[TRACEPLOT]/totaltime) << "%|";
  cout << "Time event plot:" << long(100*timeprop[EVENTPLOT]/totaltime) << "%|";
  cout << "Time part sim:" << long(100*timeprop[PARTSIM]/totaltime) << "%|";
  cout << "Time part init:" << long(100*timeprop[PARTINIT]/totaltime) << "%|";
  cout << "Time init:" << long(100*timeprop[INIT]/totaltime) << "%|";
  cout << "Time initback:" << long(100*timeprop[INITBA]/totaltime) << "%|";
*/
  cout << "\n";

  //for(i = 0; i < nind; i++){ cout << i << ":";  for(e = 0; e < nindobs[i]; e++) cout << obs_part[indobs[i][e]] << ","; cout << "\n";}
  cout.flush();
}

double Chain::Lout()                                                                                    // Outputs elements of the likelihood
{
   cout << "Lpri: " << Lpri << "   Liinit: " <<  Liinit << "   Lir: " << Lir << "   Liexp: " << Liexp << "    Linm: " << Linm
        << "   Lev_pd: " << Lev_pd << "   Lob_pd: " << Lob_pd << "    Lob_st: " << Lob_st << "   Lpop: " << Lpop << "   Lder: " << Lder << "\n";
}

void Chain::plotparamlike(long p, double min, double max)  // plots the marginal posterior distribution
{
  double v;

  ofstream likeplot("GP/likeplot");
  for(v = min; v < max; v += (max-min)/1000){
    changeparam(1,p,v);
    likeplot << v << " " << L() << "\n";
  }
}

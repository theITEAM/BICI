const short CHC = 0, BEG = 2, END = 3, OBS = 4, FIXEV = 5;
struct EVPART { double t; short ty; short ch;};

short nevp;                                                      // when doing simulation from particles this shows all the fixed times
vector <EVPART> evp;
vector<short> evpref;                                            // Gives a reference to cexist

long irev;                                                       // Used to keep track of event on the exisiting sequence
short erev, nevrev;
vector <EV> evrev;

short ppbegst, ppendst;                                          // Temporarily stores variables

short nexist;
vector<short> cexist;                                            // Keeps track of changes to c in the existing event sequence
vector<double> texist;

short nRdep, nRnotdep;                                           // Used when doing simulations
double Rdep,  Rnotdep;
vector <double>  Rdepst, Rnotdepst, rdep, rnotdep;

bool evpcomp(EVPART lhs, EVPART rhs) { return lhs.t < rhs.t; }

short nfuteref;
vector <short> futeref;

vector <short> evc;
vector <double> evt;

short nevplist;
vector<short> evplist, evplistpart;                              // Divides up into sections which are simulated

double part_sum;                                                 // Used for sampling the initial compartment for the particle
vector<double> part_sumst;

double wsum_part;

double partdt;                                                   // The Minimum time step used when doing simulation

double tent, tlea;                                               // stores the entry and leaving time for an individual

vector< vector<short> > cifcl;                                   // [ci][cf] give the classification from a change in conpartment;

short ctop;                                                       // This represents age time and fix part of the compartment on individual entry


void Chain::part_prop(short i)                                    // Makes a particle proposal on an individual
{
  short c;
  double probif, probfi, Ltoti, Ltotf, al;

  ntr_part[i]++;

  evrev = indev[i];

  ist = i;
  tent = evrev[0].t; tlea = evrev[nindev[i]-1].t;

  part_init(1,i,indtbirth[i]);
  probif = part_sim(i); if(probif == -large){ nfa_part[i]++; return;}
  probfi = part_rev(i);

  Ltoti = L();
  indchange(i);
  Ltotf = L();

  if(npartmax == 1){   // Used to check if the reverese transition is being calculated correctly
    double dd, probif2;
    vector<EV> evrevst;
    evrevst = evrev;
    evrev = indev[i];
    part_init(1,i,indtbirth[i]);
    probif2 = part_rev(i);
    dd = probif2 - probif;
    if(dd*dd > tiny){ cout << probif << " "<< probif2 << " rev\n"; emsg("rev");}
    evrev = evrevst;
  }

  al = exp(Ltotf - Ltoti + probfi - probif); //cout << al << " " << Ltoti << " " << Ltotf << " " << probif << " " << probfi << "al\n";

  if(ran() < al){
    nac_part[i]++;

    if(samp < burnin){ c = tra[indev[i][0].tr].cf; if(c == NOTALIVE) emsg("pr"); part_start[i][c%ncomps]++;}
  }
  else{
    evnew = evrev;
    indchange(i);
  }
}

double Chain::part_sim(short i)                                   // Simulates using particles
{
  short c, li, p, pp, ppbeg, ppend, evpbeg, evpend, k, nev, e, cf;
  double z, prob = 0, fac, pf;
  vector <short> plist;

  if(nmfl == 1) futev.clear();

  evpbeg = 0; ppbeg = 0;
  for(li = 0; li < nevplist; li++){
    evpend = evplist[li]; ppend = evplistpart[li]; if(ppbeg == ppend || evpbeg == evpend) emsg("npart");

    for(p = ppbeg; p < ppend; p++) prob += part_simsec(i,evpbeg,evpend,part[p]);

    partwcalc(ppbeg,ppend,evp[evpend].ty,evp[evpend].ch);
    ppbegst = ppbeg; ppendst = ppend;

    if(samp < burnin){                   // aims to adapt the number of particles so that not too many failures
      if(wsum_part == 0) fac = 1.1; else fac = 0.95;
      switch(evp[evpend].ty){
        case OBS: pf = fac*obs_part[evp[evpend].ch]; if(pf > npartmax) pf = npartmax; if(pf < 1) pf = 1; obs_part[evp[evpend].ch] = pf; break;
        //case FIXEV: pf = fac*fix_part[evp[evpend].ch]; if(pf > npartmax) pf = npartmax; if(pf < 1) pf = 1; fix_part[evp[evpend].ch] = pf; break;
      }
    }

    if(wsum_part == 0) return -large;

    evpbeg = evpend; ppbeg = ppend;
  }

  // generates the new event sequence
  z = ran()*wsum_part; pp = ppbegst; while(pp < ppendst && z > part[pp].wsum) pp++; if(pp == ppendst) emsg("pro");
  prob += log(part[pp].w/wsum_part);

  do{ plist.push_back(pp); pp = part[pp].pback;}while(pp != -1);

  evnew.clear();
  c = part[plist[plist.size()-1]].cstart;
  EV evbeg; evbeg.t = tent; if(tent == 0) evbeg.tr = trabeg+c; else{ evbeg.tr = compiftra[NOTALIVE][c]; if(evbeg.tr < 0) emsg("warn");} evnew.push_back(evbeg);

  for(k = plist.size()-1; k >= 0; k--){
    p = plist[k];

    nev = part[p].evcst.size();
    for(e = 0; e < nev; e++){
      cf = part[p].evcst[e];
      EV ev; ev.t = part[p].evtst[e]; ev.tr = compiftra[c][cf]; if(ev.tr < 0) emsg("prob"); evnew.push_back(ev);
      c = cf;
    }
  }

  EV evend; evend.t = tlea; if(tlea == tmax) evend.tr = traend+c; else{ evend.tr = compiftra[c][NOTALIVE]; if(evend.tr < 0) emsg("warn");} evnew.push_back(evend);

  return prob;
}

double Chain::part_rev(short i)                                   // Calculates the probability of create the exisiting event sequence
{
  short li, ppbeg, ppend, evpbeg, evpend, p;
  double prob = 0;

  irev = i; erev = 1; nevrev = evrev.size();

  evpbeg = 0; ppbeg = 0;
  for(li = 0; li < nevplist; li++){
    evpend = evplist[li]; ppend = evplistpart[li]; if(ppbeg == ppend || evpbeg == evpend) emsg("npart");

    prob += part_revsec(i,evpbeg,evpend,part[ppbeg]);
    for(p = ppbeg+1; p < ppend; p++) prob += part_simsec(i,evpbeg,evpend,part[p]);

    partwcalc(ppbeg,ppend,evp[evpend].ty,evp[evpend].ch);
    ppbegst = ppbeg; ppendst = ppend;
    if(wsum_part == 0) emsg("rev prob");

    evpbeg = evpend; ppbeg = ppend;
  }
  prob += log(part[ppbegst].w/wsum_part);
  if(erev != nevrev-1) emsg("erev prob");

  return prob;
}

double Chain::part_simsec(long i, short evpbeg, short evpend, PART &pa)     // simulates new sequence
{
  short e, k, co, cl, j, ch, c, pp, m, n, cc, cnew, fl, nev, ee;
  double prob = 0, t, tt, tst, tnext, tnextold, z, R, r, val, kshape, futnext, tbeg, tsta, tinside, tsto;
  long tr;

  tbeg = evp[evpbeg].t; t = tbeg;

  nfuteref = 0; futeref.clear(); evc.clear(); evt.clear();

  if(evpbeg == 0){                 // initialise
    z = ran()*part_sum; c = 0; while(c < ncomps && z > part_sumst[c]) c++; if(c == ncomps) emsg("ncomps");
    prob += log(part_start[i][c]/part_sum);
    c += ctop;
    pa.cstart = c;
    pa.pback = -1;

    if(nmfl == 1){ for(k = 0; k < nnonexpc[c]; k++) addfute(nonexpc[c][k],t);}
  }
  else{            // Does the particle filter
    z = ran()*wsum_part; pp = ppbegst; while(pp < ppendst && z > part[pp].wsum) pp++; if(pp == ppendst) emsg("pro");
    prob += log(part[pp].w/wsum_part);
    c = part[pp].cend;
    pa.pback = pp;
    if(nmfl == 1){ futeref = part[pp].futerefst; nfuteref = futeref.size();}
  }

  k = evpref[evpbeg]; if(k >= 0){ if(k == 0) co = NOTALIVE; else co = cexist[k-1]; tnextold = texist[k]-tiny;}else{ tnextold = tmax; co = NOTALIVE;}

  for(e = evpbeg; e < evpend; e++){
    tnext = evp[e+1].t;

    if(evp[e].ty == CHC){ c += evp[e].ch; evc.push_back(c); evt.push_back(t);}
    getnotdep(compleavenotdep[c]);

    do{
      while(tnextold < t){ co = cexist[k]; k++; if(k < nexist) tnextold = texist[k]-tiny; else tnextold = tmax;}
      getdep(compleavedep[c],co,t+tiny);
      R = Rdep + Rnotdep;

      if(nRdep > 0){ tt = t+partdt; if(tt > tnext) tt = tnext;}
      else tt = tnext;

      if(nmfl == 1){ if(nfuteref > 0){ futnext = futev[futeref[nfuteref-1]].t; if(tt > futnext) tt = futnext;} else futnext = large;} else futnext = large;

      if(R == 0) t = tt;
      else{
        tst = t;
        t += -log(ran())/R;
        if(t > tt) t = tt;
        else{ // actual event
          z = ran()*R;
          if(z < Rdep){
            j = 0; while(j < nRdep && z > Rdepst[j]) j++; if(j == nRdep) emsg("nRdep");
            r = rdep[j]; tr = compleavedep[c][j];
          }
          else{
            z -= Rdep;
            j = 0; while(j < nRnotdep && z > Rnotdepst[j]) j++; if(j == nRnotdep) emsg("nRnotdep");
            r = rnotdep[j]; tr = compleavenotdep[c][j];
          }
          prob += log(r);
          c = tra[tr].cf;

          evc.push_back(c); evt.push_back(t);

          if(nmfl == 1){ for(n = 0; n < nnonexptra[tr]; n++) addfute(nonexptra[tr][n],t);}    // adds on non-Markovian future events

          getnotdep(compleavenotdep[c]);
        }
        prob -= (t-tst)*R;
      }

      if(t == futnext){   // Future event
        n = futeref[nfuteref-1]; futeref.pop_back(); nfuteref--;
 
        tsta = futev[n].tnow; tinside = large;
        cl = futev[n].cl;
        cc = c;
        ee = evc.size()-1;                                                                         // Checks if there has been some intervening event in that classification
        while(ee >= 0 && evt[ee] > tsta){ cnew = evc[ee]; if(cifcl[cc][cnew] == cl) tinside = tsto; cc = cnew; tsto = evt[ee]; ee--;}

        if(tsta < tbeg){
          pp = pa.pback;
          do{
            nev = part[pp].evcst.size();
            if(nev > 0){
              ee = nev-1;
              while(ee >= 0 && part[pp].evtst[ee] > tsta){ cnew = part[pp].evcst[ee]; if(cifcl[cc][cnew] == cl) tinside = tsto; cc = cnew; tsto = part[pp].evtst[ee]; ee--;}
              if(part[pp].evtst[ee] <= tsta) break;
            }
            pp = part[pp].pback;
          }while(pp != -1);
        }

        if(tinside == large){   // Adds the future event
          tr = futev[n].tr;
          switch(tra[tr].type){
            case GAMMA_TR: prob += gammaprob(futev[n].t-futev[n].tnow,futev[n].val1,futev[n].val2); break;
            case WEI_TR: prob += weibullprob(futev[n].t-futev[n].tnow,futev[n].val1,futev[n].val2); break;
          }

          c += tra[tr].dc;
          evc.push_back(c); evt.push_back(t);
          for(n = 0; n < nnonexptra[tr]; n++) addfute(nonexptra[tr][n],t);
        }
        else{                   // There is an intervening event
          switch(tra[tr].type){
            case GAMMA_TR: prob += gammaup(tinside-futev[n].tnow,futev[n].val1,futev[n].val2); break;
            case WEI_TR: prob += weibullup(tinside-futev[n].tnow,futev[n].val1,futev[n].val2); break;
          }
        }
      }
    }while(t < tnext);
  }

  pa.cend = c; pa.evcst = evc; pa.evtst = evt;
  if(nmfl == 1){
    if(evpend < nevp-1) pa.futerefst = futeref;
    else{   // Gives probability of remaining events
      for(j = 0; j < nfuteref; j++){
        n = futeref[j];
        switch(tra[futev[n].tr].type){
          case GAMMA_TR: prob += gammaup(tmax-futev[n].tnow,futev[n].val1,futev[n].val2); break;
          case WEI_TR: prob += weibullup(tmax-futev[n].tnow,futev[n].val1,futev[n].val2); break;
        }
      }
    }
  }

  return prob;
}

double Chain::part_revsec(long i, short evpbeg, short evpend, PART &pa)     // probability of simulating the existing event sequence
{
  short e, k, co, cl, j, ch, c, pp, m, n, cc, cnew, fl, nev, ee;
  double prob = 0, t, tt, tst, tnext, tnextold, z, R, r, val, kshape, futnext, tbeg, tsta, tinside, tsto;
  long tr;

  tbeg = evp[evpbeg].t; t = tbeg;

  nfuteref = 0; futeref.clear(); evc.clear(); evt.cleart();
  if(evpbeg == 0){                 // initialise
    c = (tra[evrev[0].tr].cf)%ncomps;
    prob += log(part_start[i][c]/part_sum);
    c += ctop;
    pa.cstart = c;
    pa.pback = -1;

    if(nmfl == 1){ for(k = 0; k < nnonexpc[c]; k++) addfuterev(nonexpc[c][k],t);}
  }
  else{            // Does the particle filter
    pp = ppbegst;
    prob += log(part[pp].w/wsum_part);
    c = part[pp].cend;
    pa.pback = pp;
    if(nmfl == 1){ futeref = part[pp].futerefst; nfuteref = futeref.size();}
  }

  k = evpref[evpbeg]; if(k >= 0){ if(k == 0) co = NOTALIVE; else co = cexist[k-1]; tnextold = texist[k]-tiny;}else{ tnextold = tmax; co = NOTALIVE;}
  for(e = evpbeg; e < evpend; e++){
    tnext = evp[e+1].t;

    if(evp[e].ty == CHC){ c += evp[e].ch; evc.push_back(c); evt.push_back(t);}
    getnotdep(compleavenotdep[c]);

    do{
      while(tnextold < t){ co = cexist[k]; k++; if(k < nexist) tnextold = texist[k]-tiny; else tnextold = tmax;}

      getdep(compleavedep[c],co,t+tiny);

      R = Rdep + Rnotdep;

      if(nRdep > 0){ tt = t+partdt; if(tt > tnext) tt = tnext;}
      else tt = tnext;

      if(nmfl == 1){ if(nfuteref > 0){ futnext = futev[futeref[nfuteref-1]].t; if(tt > futnext) tt = futnext;} else futnext = large;} else futnext = large;

      if(R == 0){
        if(erev < nevrev && evrev[erev].t < tt) emsg("event missing");
        t = tt;
      }
      else{
        tst = t;
        if(erev < nevrev && evrev[erev].t < tt){ // Actual event occurs
          t = evrev[erev].t; tr = evrev[erev].tr; erev++;
          if(transdep[tra[tr].eq] == 1){
            j = 0; while(j < nRdep && compleavedep[c][j] != tr) j++; if(j == nRdep) emsg("rev prob");
            r = rdep[j];
          }
          else{
            j = 0; while(j < nRnotdep && compleavenotdep[c][j] != tr) j++; if(j == nRnotdep) emsg("rev prob");
            r = rnotdep[j];
          }
          prob += log(r);
          c = tra[tr].cf;

          evc.push_back(c); evt.push_back(t);

          if(nmfl == 1){ for(n = 0; n < nnonexptra[tr]; n++) addfuterev(nonexptra[tr][n],t);}    // adds on non-Markovian future events

          getnotdep(compleavenotdep[c]);
        }
        else t = tt;
        prob -= (t-tst)*R;
      }

      if(t == futnext){   // Future event
        n = futeref[nfuteref-1]; futeref.pop_back(); nfuteref--;

        tsta = futev[n].tnow; tinside = large;
        cl = futev[n].cl;
        cc = c;
        ee = evc.size()-1;                                                                         // Checks if there has been some intervening event in that classification
        while(ee >= 0 && evt[ee] > tsta){ cnew = evc[ee]; if(cifcl[cc][cnew] == cl) tinside = tsto; cc = cnew; tsto = evt[ee]; ee--;}

        if(tsta < tbeg){
          pp = pa.pback;
          do{
            nev = part[pp].evcst.size();
            if(nev > 0){
              ee = nev-1;
              while(ee >= 0 && part[pp].evtst[ee] > tsta){ cnew = part[pp].evcst[ee]; if(cifcl[cc][cnew] == cl) tinside = tsto; cc = cnew; tsto = part[pp].evtst[ee]; ee--;}
              if(part[pp].evtst[ee] <= tsta) break;
            }
            pp = part[pp].pback;
          }while(pp != -1);
        }

        if(tinside == large){   // Adds the future event
          tr = futev[n].tr;

          if(erev == nevrev || tra[tr].dc != tra[evrev[erev].tr].dc) emsg("fit prob");
          erev++;

          switch(tra[tr].type){
            case GAMMA_TR: prob += gammaprob(futev[n].t-futev[n].tnow,futev[n].val1,futev[n].val2); break;
            case WEI_TR: prob += weibullprob(futev[n].t-futev[n].tnow,futev[n].val1,futev[n].val2); break;
          }

          c += tra[tr].dc;
          evc.push_back(c); evt.push_back(t);
          for(n = 0; n < nnonexptra[tr]; n++) addfute(nonexptra[tr][n],t);
        }
        else{                   // There is an intervening event
          switch(tra[tr].type){
            case GAMMA_TR: prob += gammaup(tinside-futev[n].tnow,futev[n].val1,futev[n].val2); break;
            case WEI_TR: prob += weibullup(tinside-futev[n].tnow,futev[n].val1,futev[n].val2); break;
          }
        }
      }
    }while(t < tnext);
  }

  pa.cend = c; pa.evcst = evc; pa.evtst = evt;
  if(evpend < nevp-1) pa.futerefst = futeref;
  else{   // Gives probability of remaining events
    for(j = 0; j < nfuteref; j++){
      n = futeref[j];
      switch(tra[futev[n].tr].type){
        case GAMMA_TR: prob += gammaup(tmax-futev[n].tnow,futev[n].val1,futev[n].val2); break;
        case WEI_TR: prob += weibullup(tmax-futev[n].tnow,futev[n].val1,futev[n].val2); break;
      }
    }
  }

  return prob;
}

void Chain::addfute(long tr, double t)                                             // Adds a future event
{
  short m;
  double tt, val, kshape;
  FUTEV fev;

  m = 0; while(m < nfuteref && t < futev[futeref[m]].t) m++;
  futeref.insert(futeref.begin()+m,futev.size());
  nfuteref++;

  fev.tr = tr; fev.tnow = t;
  val = nmeq_val[tra[tr].eq];
  switch(tra[tr].type){
    case GAMMA_TR: kshape = nmeq_val[tra[tr].eqshape]; tt = t+gammasamp(kshape,kshape/val); fev.val1 = kshape; fev.val2 = kshape/val; break;
    case WEI_TR: kshape = nmeq_val[tra[tr].eqshape]; tt = t+weibullsamp(val,kshape); fev.val1 = val; fev.val2 = kshape; break;
    case FIXED_TR: tt = t+val; break;
  }
  fev.t = tt;
  futev.push_back(fev);
}

void Chain::addfuterev(long tr, double t)                                             // Adds the future event consistent with exisiting event sequence
{
  short m, e, cl;
  double tt, val, kshape, tmin;
  FUTEV fev;

  m = 0; while(m < nfuteref && tt < futev[futeref[m]].t) m++;
  futeref.insert(futeref.begin()+m,futev.size());
  nfuteref++;

  val = nmeq_val[tra[tr].eq];
  switch(tra[tr].type){
    case GAMMA_TR: kshape = nmeq_val[tra[tr].eqshape]; fev.val1 = kshape; fev.val2 = kshape/val; break;
    case WEI_TR: kshape = nmeq_val[tra[tr].eqshape]; fev.val1 = val; fev.val2 = kshape; break;
  }

  fev.tr = tr; fev.tnow = t;
  cl = tra[tr].cl;
  e = erev; while(e < nevrev-1 && tra[evrev[e].tr].cl != cl) e++;
  if(e == nevrev-1 || tra[evrev[e].tr].dc != tra[tr].dc){
    tmin = evrev[e].t;
    switch(tra[tr].type){
      case GAMMA_TR: do{ tt = t+gammasamp(kshape,kshape/val);}while(tt < tmin); break;
      case WEI_TR: do{ tt = t+weibullsamp(val,kshape);}while(tt < tmin);  break;
      case FIXED_TR: tt = t+val; if(tt < tmin) emsg("PROB"); break;
    }
  }
  else tt = evrev[e].t;

  fev.t = tt;
  futev.push_back(fev);
}

void Chain::partwcalc(short ppbeg, short ppend, short ty, short ch)   // Calculates the particle weightings at the observations
{
  short p, k, cl, tr, c, vi;
  double w;

  wsum_part = 0;
  for(p = ppbeg; p < ppend; p++){
    c = part[p].cend;

    w = 1;
    switch(ty){
      case OBS:
        k = obsprob_cope[ch][c];
        if(k < 0) w = 0; else w *= exp(cope_val[k]);
        break;

      case FIXEV:
        /*
        cl = fixevcl[ch]; tr = fixevtr[ch];
        vi = (c/classmult[cl])%nclassval[cl];
        if(vi != transi[cl][tr]) w[p] = 0;
        else{
          c += (transf[cl][tr] - vi)*classmult[cl];
          indev_part[p].push_back(c); indevt_part[p].push_back(fixevt[ch]);
          end[p] = c;
        }
        */
        break;

        case END: //makes sure final state can die
          //if(ch == 1 && sinktrans[c] == -1) w[p] = 0;
          break;
    }
    wsum_part += w;
    part[p].w = w; part[p].wsum = wsum_part;
  }
}

void Chain::getdep(vector <short> &vec, short co, double t)                   // Gets the total rate of dependent transitions
{
  short j, eq;
  Rdep = 0; nRdep = vec.size();
  for(j = 0; j < nRdep; j++){
    eq = tra[vec[j]].eq; if(eq < 0) emsg("new eq");
 
    rdep[j] = transra(eq,t,co);//+rshift; if(startmode == 1){ if(rdep[j] > 100/tmax) rdep[j] = 100/tmax;}
    Rdep += rdep[j];
    Rdepst[j] = Rdep;
  }
}

void Chain::getnotdep(vector <short> &vec)                                  // Gets the total rate of not dependent transition
{
  short j, eq;

  Rnotdep = 0; nRnotdep = vec.size();
  for(j = 0; j < nRnotdep; j++){
    eq = tra[vec[j]].eq; if(eq < 0) emsg("new eq");
    rnotdep[j] = transnotdepeq_val[transdepref[eq]];
    Rnotdep += rnotdep[j];
    Rnotdepst[j] = Rnotdep;
  }
}

void Chain::part_init(short exist, short i, double tbirth)          // Initialises a particle proposal
{
  short cl, e, tr, v, vi, vf, ag, ti, ob, f, nparttot, k, c;
  double t, pf, dt;

  cexist.clear(); texist.clear();
  if(exist == 1){ for(e = 0; e < nindev[i]; e++){ cexist.push_back(tra[indev[i][e].tr].cf); texist.push_back(indev[i][e].t);}}
  nexist = cexist.size();

  evp.clear(); nevp = 0;
  evp.push_back(EVPART ()); evp[nevp].ty = BEG; evp[nevp].t = tent; nevp++;

  ctop = 0;
  if(tbirthfl == 1){
    ag = 0; while(ag < nage && tbirth + age[ag] < tent) ag++;
    ctop += ag*classmult[agecl];
    while(ag < nage && tbirth + age[ag] < tlea){ evp.push_back(EVPART ()); evp[nevp].ty = CHC; evp[nevp].ch = classmult[agecl]; evp[nevp].t = tbirth + age[ag]; nevp++; ag++;}
  }

  ti = 0; while(ti < nsettime && settime[ti] < tent) ti++;
  ctop += ti*classmult[settimecl];
  while(ti < nsettime && settime[ti] < tlea){ evp.push_back(EVPART ()); evp[nevp].ty = CHC; evp[nevp].ch = classmult[settimecl]; evp[nevp].t = settime[ti]; nevp++; ti++;}

  if(i < nind){
    for(e = 0; e < nindfixev[i]; e++){ evp.push_back(EVPART ()); evp[nevp].ty = FIXEV; evp[nevp].ch = indfixev[i][e]; evp[nevp].t = fixevt[indfixev[i][e]]; nevp++;}
    for(ob = 0; ob < nindobs[i]; ob++){ evp.push_back(EVPART ()); evp[nevp].ty = OBS; evp[nevp].ch = indobs[i][ob]; evp[nevp].t = indobst[i][ob]; nevp++;}
  }

  f = 0; while(f < nfixtrange && fixtrange[f] < tent) f++;
  ctop += f*classmult[fixcl];
  while(f < nfixtrange && fixtrange[f] < tlea){ evp.push_back(EVPART ()); evp[nevp].ty = CHC; evp[nevp].ch = classmult[fixcl]; evp[nevp].t = fixtrange[f]; nevp++; f++;}

  evp.push_back(EVPART ()); evp[nevp].ty = END; evp[nevp].t = tlea; if(tlea < tmax) evp[nevp].ch = 1; else evp[nevp].ch = 0; nevp++;

  sort(evp.begin(),evp.end(),evpcomp);

  //for(e = 0; e < nevp; e++) cout << e << " " << evp[e].t << " "<<  evp[e].ty << " " << evp[e].ch << " evp\n";
  //short fi; for(e = 0; e < nindfixev[i]; e++){ fi = indfixev[i][e]; cout <<  fixevt[fi] << ","; cout << "fi\n"; }

  nevplist = 0; evplist.clear(); evplistpart.clear(); evpref.clear();
  k = 0;
  nparttot = 0;
  for(e = 0; e < nevp; e++){
    if(exist == 1){
      t = evp[e].t;
      while(k < cexist.size() && texist[k] < t) k++;
      evpref.push_back(k);
    }
    else evpref.push_back(-1);

    if(e == nevp-1 || evp[e].ty == OBS || evp[e].ty == FIXEV){
      evplist.push_back(e);
      switch(evp[e].ty){
        case OBS: pf = obs_part[evp[e].ch]; break;
        case FIXEV: pf = fix_part[evp[e].ch]; break;
        default: pf = 1; break;
      }
      nparttot += short(pf); evplistpart.push_back(nparttot);

     nevplist++;
    }
  }

  if(i < nind){ if(nevplist != nindobs[i]+nindfixev[i]+1) emsg("wrongh");}
  else{ if(nevplist != 1) emsg("wrongg");}

  partdt = tmax/50.0; dt = (tlea-tent)/10; if(dt < partdt) partdt = dt;

  part.resize(nparttot);
  part_sum = 0; for(c = 0; c < ncomps; c++){ part_sum += part_start[i][c]; part_sumst[c] = part_sum;}
}

void Chain::chainpartinit()             // Initialises quantities used when sampling from partivles
{
  short ob, fi, c;
  long i;

  obs_part.resize(nobs); for(ob = 0; ob < nobs; ob++) obs_part[ob] = npartstart;    // initialises the particle numbets

  fix_part.resize(nfixev); for(fi = 0; fi < nfixev; fi++) fix_part[fi] = 1;//npartmax;

  part_start.resize(nind); for(i = 0; i < nind; i++){ part_start[i].resize(ncomps); for(c = 0; c < ncomps; c++) part_start[i][c] = 50.0/ncomps;}
}

void partinit()                                 // Initialises quantities used when sampling from partivles
{
  short eq, c, cc, cl, i, f, ci, cf, p, j, d, depmax, ndmax;
  long tr, tr2;

  compleavedep.resize(ncomp+1); compleavenotdep.resize(ncomp+1);    // Used when simulating particles
  for(tr = 0; tr < ntra; tr++){
    if(tra[tr].type == EXP_TR && tra[tr].fix == 0){
      eq = tra[tr].eq;
      switch(transdep[eq]){
        case 0: compleavenotdep[tra[tr].ci].push_back(tr); break;
        case 1: compleavedep[tra[tr].ci].push_back(tr); break;
      }
    }
  }
  depmax = 0; ndmax = 0;
  ncompleavedep.resize(ncomp+1); for(c = 0; c <= ncomp; c++){ ncompleavedep[c] = compleavedep[c].size(); if(ncompleavedep[c] > depmax) depmax = ncompleavedep[c];}
  ncompleavenotdep.resize(ncomp+1); for(c = 0; c <= ncomp; c++){ ncompleavenotdep[c] = compleavenotdep[c].size(); if(ncompleavenotdep[c] > ndmax) ndmax = ncompleavenotdep[c];}

  Rdepst.resize(depmax); Rnotdepst.resize(ndmax);
  rdep.resize(depmax); rnotdep.resize(ndmax);

  part_sumst.resize(ncomps);

  for(p = 0; p < npopnum; p++){                                                               // works out which rates need to be recalculated
    popnumtake.push_back(vector<double>());
    for(c = 0; c <= ncomp; c++) popnumtake[p].push_back(0);
    for(j = 0; j < popnumterm[p].size(); j++) popnumtake[p][popnumterm[p][j]] -= popnumtermweight[p][j];
  }

  for(d = 0; d < ntransdepeq; d++){
    transdepeqrecalc.push_back(vector<short>());
    eq = transdepeq[d];
    for(c = 0; c < ncomp; c++){
      for(p = 0; p < neq_popnum[eq]; p++) if(popnumtake[eq_popnum[eq][p]][c] != 0) break;
      if(p < neq_popnum[eq]) transdepeqrecalc[d].push_back(1);
      else transdepeqrecalc[d].push_back(0);
    }
  }

  nonexpc.resize(ncomp+1); nnonexpc.resize(ncomp+1);
  for(tr = 0; tr < ntra; tr++){
    switch(tra[tr].type){
      case FIXED_TR: case GAMMA_TR: case WEI_TR:
        nonexpc[tra[tr].ci].push_back(tr);
        break;
    }
  }
  for(c = 0; c <= ncomp; c++) nnonexpc[c] = nonexpc[c].size();

  nonexptra.resize(ntra); nnonexptra.resize(ntra);
  for(tr = 0; tr < ntra; tr++){
    cf = tra[tr].cf;
    for(j = 0; j < nnonexpc[cf]; j++){
      tr2 = nonexpc[cf][j];
      if(tra[tr2].cl == tra[tr].cl) nonexptra[tr].push_back(tr2);
    }
  }
  for(tr = 0; tr < ntra; tr++) nnonexptra[tr] = nonexptra[tr].size();

  nonexptrarev.resize(ntra); nnonexptrarev.resize(ntra);
  for(tr = 0; tr < ntra; tr++){
    ci = tra[tr].ci;
    for(j = 0; j < nnonexpc[ci]; j++){
      tr2 = nonexpc[ci][j];
      if(tra[tr2].cl == tra[tr].cl) nonexptrarev[tr].push_back(tr2);
    }
  }
  for(tr = 0; tr < ntra; tr++) nnonexptrarev[tr] = nonexptrarev[tr].size();

  cifcl.resize(ncomp+1);
  for(c = 0; c <= ncomp; c++){ cifcl[c].resize(ncomp+1); for(cc = 0; cc <= ncomp; cc++) cifcl[c][cc] = -1;}
  for(c = 0; c < ncomp; c++){
    for(cl = 0; cl < nclass; cl++){
      i = compval[c][cl]; for(f = 0; f < nclassval[cl]; f++) if(i != f) cifcl[c][c+(f-i)*classmult[cl]] = cl;
    }
  }
}

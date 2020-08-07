

struct EVI { long i; double t; long tr;};
bool compEVI(EVI lhs, EVI rhs) { return lhs.t < rhs.t; }

void Chain::check(long num)                       // Used to chech the algorithm is working
{
  long e, nev, j, f, cl, c, p, ce;
  long i, fev, tr;
  double t;

  errorcode = num;
  //cout << "check:" << num << "\n";

  timeprop[CHECK] -= clock();

  if(isnan(Lpri)) emsg("Lpri is NaN",num);
  if(isnan(Liinit)) emsg("Liinit is NaN");
  if(isnan(Lir)) emsg("Lir is NaN");
  if(isnan(Liexp)) emsg("Liexp is NaN");
  if(isnan(Linm)) emsg("Linm is NaN");
  if(isnan(Lev_pd)) emsg("Lev_pd is NaN");
  if(isnan(Lob_pd)) emsg("Lob_pd is NaN");
  if(isnan(Lob_st)) emsg("Lob_st is NaN");
  if(isnan(Lpop)) emsg("Lpop is NaN");

  for(i = 0; i < nindtot; i++) checkevseqsimp(indev[i]);

  for(p = 0; p < npopm; p++) if(popmap[p] != 0) emsg("Popmap");
  for(p = 0; p < nderm; p++) if(dermap[p] != 0) emsg("Popmap");

  if(nindev.size() != nindtot) emsg("nindevp prb");
  if(indev.size() != nindtot) emsg("indevp prb");

  for(i = 0; i < nind; i++){  // checks fixed events
    for(j = 0; j < nindfixev[i]; j++){
      fev = indfixev[i][j];
      t = fixev[fev].t;
      //for(e = 0; e < nindev[i]; e++) cout << indev[i][e].t << " "<< t << "com\n";

      nev = nindev[i]; e = 0; while(e < nev && indev[i][e].t < t) e++;
      if(e == nev){ cout << t << " t\n"; emsg("fixev not there");}
      if(indev[i][e].t != t){ cout << indid[i] << " " << indev[i][e].t << " " << t <<" y\n"; emsg("fixev not there2");}
      if(tra[indev[i][e].tr].dc != fixev[fev].dc) emsg("wrong fixed event");
      if(tra[indev[i][e].tr].like != fixev[fev].like){
        cout << i << " " << e << " " << indid[i] << " " << e <<" " << t << " " << tra[indev[i][e].tr].like  << " " << fixev[fev].like << " " << nindfixev[i] << " li\n";
        emsg("fev.like does not match");
      }

      ce = fixev[fev].capev;
      if(ce >= 0){
        if(tra[indev[i][e].tr].capev != ce) emsg("ce prob");

        c = tra[indev[i][e].tr].ci;
        if(capevfilt[ce][c] != 1) emsg("capev prob");
      }
    }

    tr = indev[i][0].tr; if(tra[tr].like != indlikeenter[i]) emsg("likeenter problem");
    tr = indev[i][nindev[i]-1].tr; if(tra[tr].like != indlikeleave[i]) emsg("likeleave problem");

    for(e = 0; e < nindev[i]; e++){
      tr = indev[i][e].tr;
      if(tra[tr].capev >= 0 && tra[tr].type >= 0){
        if(capevall[tra[tr].capev] == 1){
          if(tra[tr].ci == NOTALIVE){
            if(indev[i][e].t != indfixtenter[i]){ cout << i << " "<< e << " " << indev[i][e].t << "entt\n"; emsg("fixed enter time wrong");}
            if(tra[tr].capev != indfixentercapev[i]) emsg("capev not match");
          }
          else{
            if(tra[tr].cf == NOTALIVE){
              if(indev[i][e].t != indfixtleave[i]){ cout << indid[i] << " " << indev[i][e].t << " " <<indfixtleave[i] << " lea\n";  emsg("fixed leave time wrong");}
              if(tra[tr].capev != indfixleavecapev[i]) emsg("capev not match");
            }
            else{
              for(j = 0; j < nindfixev[i]; j++){ 
                if(fixev[indfixev[i][j]].t == indev[i][e].t){
                  if(fixev[indfixev[i][j]].like != tra[tr].like) emsg("Like do not agree");
                  break;
                }
              }
              if(j == nindfixev[i]){ cout <<indid[i] << " "  <<  i << " " << tra[tr].ci << " " <<  tra[tr].cf << "i\n"; emsg("fixed event should not be there");}
            }
          }
        }
      }
    }
  }

  for(i = nind; i < nindtot; i++){
    for(e = 0; e < nindev[i]; e++){
      tr = indev[i][e].tr;
      if(tra[tr].capev >= 0 && tra[tr].type >= 0){
        if(capevall[tra[tr].capev] == 1){
          for(ce = 0; ce < ncapev; ce++) cout << capev_val[ce] << " " << capev_valoneminus[ce] << " ce\n";

          oe("TT",indev[i]);
          cout << i << " " << samp << " " << num << " samp\n"; 
          emsg("unobserved fixed");
        }
      }
    }
  }

  for(i = 0; i < nindtot; i++){
    if(nindev[i] != indev[i].size()) emsg("si2");

    getlifespan(i);

    if(i < nind){
      if(tent > firstobst[i]) emsg("firstobst");
      if(tlea < lastobst[i]) emsg("lastobst");
    }

    if(tbirthfl == 1){
      if(indtbirth[i] < tent-age[nage]){ cout << i << " " << indid[i] << "Id\n"; emsg("too early");}
      if(indtbirth[i] > tent) emsg("too late");
    }
    else{ if(indtbirth[i] != tent){ cout <<indtbirth[i] << "  "<< tent << " yy\n";  emsg("tbirth prob");}}

    checkevseq(indev[i]);
  }

  checklike(num);
  checkderive();
  //checkparamgrad();

  for(i = 0; i < nindtot; i++) if(nindev[i] < 2) emsg("prob 1");

  /*
  for(e = 0; e < fixnmlist.size(); e++){ f = fixnmlist[e]; fixnmmap[fixev[f].cl][fixev[f].trai][fixev[f].traf].pop_back();}
  fixnmlist.clear();
  for(cl = 0; cl < nclass; cl++){
    for(i = 0; i < nclassval[cl]; i++){
      for(j = 0; j < nclassval[cl]; j++){ if(fixnmmap[cl][i][j].size() != 0) emsg("pr");}
    }
  }
  */

  //cout << "checkendn\n";
  timeprop[CHECK] += clock();
  errorcode = -1;
}

void checkevseq(vector<EV> &vec)       // CHecks that an event sequnece is consistent and correct
{
  long e, ee, cl, inside, c, li, ctop, ti, f, k, n;
  long tr;
  double t;

  if(vec.size() == 2 && vec[0].tr ==  tranull && vec[1].tr == tranull) return;

  tr = vec[0].tr; if(tra[tr].ci != NOTALIVE) emsg("no begin1");
  c = tra[tr].cf;
  if(vec[0].t == 0){
    if(tr != trabeg+c) emsg("beg prob h");
  }
  else{
    if(tr == trabeg+c) emsg("beg prob 3");
    if(tr != compiftra[NOTALIVE][c] && tr != compiftra[NOTALIVE][c]+moventra) emsg("beg prob2");
  }
  tr = vec[long(vec.size())-1].tr; if(tra[tr].cf != NOTALIVE) emsg("end prob");

  /*
  e = 0; while(e < vec.size() && tra[vec[e].tr].cf != NOTALIVE) e++;
  if(e == vec.size()) emsg("END prob");

  tr = vec[e].tr; c = tra[tr].ci;
  if(vec[e].t == tmax){
    if(tr != traend+c) emsg("end prob");
  }
  else{
    if(tr != compiftra[c][NOTALIVE]) emsg("end prob");
  }
  */
  
  //e++; while(e <  vec.size()){ if(tra[vec[e].tr].type >= 0) emsg("prob after end"); e++;}

  for(e = 0; e < long(vec.size())-1; e++) if(vec[e].tr < 0 || vec[e].tr >= ntra) emsg("tra out of range");
  for(e = 0; e < long(vec.size())-1; e++){  // does a simple check
    if(vec[e].t > vec[e+1].t){ oe("ord",vec); emsg("order1");}
    if(vec[e].t > vec[e+1].t-evdtmin){ emsg("events within evdtmin range",simnum,samp);}
  }
  c = NOTALIVE;
  for(e = 0; e < vec.size(); e++){
    tr = vec[e].tr;
    if(e > 0 && tra[tr].ci != c) emsg("consi vec1");
    c = tra[tr].cf;

    if(c != NOTALIVE){
      t = vec[e].t;
      ti = 0; while(ti < nsettime && settime[ti] <= t) ti++;
      if(compval[c][settimecl] != ti) emsg("setime wrong h");

      f = 0; while(f < ncapevtrange && capevtrange[f] <= t) f++;
      if(compval[c][capevcl] != f) emsg("fix wrong");
    }
  }

 
  // checks all nm transitions have an end point
  /*
  c = tra[vec[0].tr].cf;
  //for(k = 0; k < nnonexpc[c]; k++) checkfute(vec[0].t,nonexpc[c][k],1,vec);
  for(e = 0; e < vec.size(); e++){
    tr = vec[e].tr;
    for(n = 0; n < nnonexptra[tr]; n++) checkfute(vec[e].t,nonexptra[tr][n],e+1,vec);
  }

  for(e = 0; e < vec.size(); e++){
    switch(tra[vec[e].tr].type){
      case FIXED_TR: case GAMMA_TR: case WEI_TR:
        t = vec[e].tbeg; cl = tra[vec[e].tr].cl; inside = 0;
        ee = e-1; while(ee > 0 && vec[ee].t > t){ if(tra[vec[ee].tr].cl == cl){ cout << t << " " << vec[ee].t << "tt\n";inside++;} ee--;}
        if(vec[ee].t != t) emsg("tbeg prob\n");
        if(inside != 0) emsg("inside h");
        break;

      case FIXEDEND_TR: case GAMMAEND_TR: case WEIEND_TR:
        break;
    }
  }
  */
}


long checkevdt()
{
  long e;
  for(e = 0; e < long(evnew.size())-1; e++){ if(evnew[e+1].t-evnew[e].t < evdtmin) return 1;}
  return 0;
}

void checkevseqsimp(vector<EV> &vec)       // CHecks that an event sequnece is consistent and correct
{
  long e, ee, cl, inside, c, li, ctop, ti, f, k, n;
  long tr;
  double t;

  if(vec.size() == 2 && vec[0].tr ==  tranull && vec[1].tr == tranull) return;

  for(e = 0; e < long(vec.size())-1; e++) if(vec[e].tr < 0 || vec[e].tr >= ntra){ oe("jj",vec); cout<<  vec.size() <<  "si\n"; emsg("tra out of range");}
  for(e = 0; e < long(vec.size())-1; e++){  // does a simple check
    if(vec[e].t >= vec[e+1].t) emsg("order ch");
  }

  c = NOTALIVE;
  for(e = 0; e < vec.size(); e++){
    tr = vec[e].tr;
    if(e > 0 && tra[tr].ci != c){
      oe("HH",vec);
      emsg("consi vec2");
    }
    c = tra[tr].cf;
  }
}

void Chain::checklike(long num)                           // Checks that the likelihood is correctly calculated
{
  long cap, tr, c, p, flag, j, ob, eq, k, d, cf, ci, stim, ce, fl;
  long nev, e, i, ee, ed, nindinitch, ncompinitch[ncompswa], numd, numnd, cev_num[ncapev], cev_oneminus[ncapev];                          // The number of events observed
  double popnum[npopnum], t, tt, pd, R, dd;
  double Lr=0, Lexp=0, Linit=0, Lpr=0, Lpopch=0, Lnm=0, Le_pd=0, Lo_pd=0, Lo_st=0;
  vector <EVI> ev;
  vector <long> ncompind;
  vector <long> stat;
  vector< vector<long> > depeq_ev_evflag;

  for(ce = 0; ce < ncapev; ce++){ cev_num[ce] = 0; cev_oneminus[ce] = 0;}

  for(i = 0; i < nindtot; i++){
    for(e = 0; e < nindev[i]; e++){
      EV eve = indev[i][e];
      EVI evi; evi.i = i; evi.t = eve.t; evi.tr = eve.tr;
      if(tra[eve.tr].type >= 0) ev.push_back(evi);
    }
  }
  sort(ev.begin(),ev.end(),compEVI);

  ncompind.resize(ncomp+1); for(c = 0; c < ncomp; c++) ncompind[c] = 0; ncompind[NOTALIVE] = nindtot;
  stat.resize(nindtot); for(i = 0; i < nindtot; i++) stat[i] = NOTALIVE;
 
  for(e = 0; e < ndepeq_ev; e++) depeq_ev_evflag.push_back(vector<long>());

  for(d = 0; d < ntransdepeq; d++){
    e = depeqdiv[d][0];
    eq = transdepeq[d];

    do{
      for(k = 0; k < ndepeq_ev_ev[e]; k++) depeq_ev_evflag[e].push_back(0);
      e = depeq_evnext[e];
    }while(e != -1);
  }

  cap = 0; t = 0; stim = 0;
  nev = ev.size();
  numd = 0; numnd = 0;
  for(e = 0; e <= nev; e++){
    if(e == nev) tt = tmax; else tt = ev[e].t;
    if(stim < nsettime && t >= settime[stim]) stim++;

    for(p = 0; p < npopnum; p++){ popnum[p] = 0; for(i = 0; i < popnumterm[p].size(); i++) popnum[p] += ncompind[popnumterm[p][i]]*popnumtermweight[p][i];}
    //for(p = 0; p < npopnum; p++) cout << t << " " << popnumname[p] << " " << popnum[p] << " j\n"; 

    while(cap < ncap && capt[cap] < tt){    // does the observations
      for(i = 0; i < nindtot; i++){
        flag = 0;
        if(i < nind){
          for(j = 0; j < nindobs[i]; j++){ ob = indobs[i][j]; if(obscap[ob] == cap){ flag = 1; break;}}
        }
        c = stat[i];

        if(stat[i] != NOTALIVE && flag == 0){
          eq = capprobeqn[cap][c];
          if(eq >= 0){
            pd = calcobsprob(eq);
            if(pd == 1) emsg("notcap");
          }
        }

        if(stat[i] == NOTALIVE){
          //if(flag == 1) emsg("individual not at observation");
        }
        else{
          eq = capprobeqn[cap][c];
          if(eq >= 0){
            pd = calcobsprob(eq);
            if(flag == 1){ numd++; Lo_pd += log(pd);}
            else{ numnd++; Lo_pd += log(1-pd);}
          }
          else{
            //if(flag == 1) emsg("cap prob zero");
          }
        }
      }
      cap++;
    }
 
    //cout << Lo_pd << " Lo\n";
    R = 0;
    for(i = 0; i < nindtot; i++){                    // transitions
      c = stat[i];
      if(c != NOTALIVE){
        for(j = 0; j < ncompleave[c]; j++){ tr = compleave[c][j]; if(tra[tr].type == EXP_TR) R += ratecalc(tra[tr].eq,popnum,param);}
      }
    }

    c = NOTALIVE;                                    // sources
    for(j = 0; j < ncompleave[c]; j++){
      tr = compleave[c][j];
      if(compval[tra[tr].cf][settimecl] == stim) R += ratecalc(tra[tr].eq,popnum,param);
    }

    Lexp += R*(tt-t);

    t = tt;
    if(e == nev) break;

    i = ev[e].i;
    tr = ev[e].tr;

    ce = tra[tr].capev;     // capture of events

    if(ce >= 0){
      if(i < nind){
        fl = 0;
        for(j = 0; j < nindfixev[i]; j++){ if(fixev[indfixev[i][j]].t == t && fixev[indfixev[i][j]].capev == ce){ fl = 1; break;}}
        if(t == indfixtenter[i] && indfixentercapev[i] == ce) fl = 1;
        if(t == indfixtleave[i] && indfixleavecapev[i] == ce) fl = 1;

        if(fl == 1) cev_num[ce]++; else cev_oneminus[ce]++;
      }
      else cev_oneminus[ce]++;
    }

    //cout << i << " " <<  tra[tr].ci << " " << tra[tr].cf << " " << t << " "<< tra[tr].like <<  "trans\n";

    if(tra[tr].type == EXP_TR){
      eq = tra[tr].eq;
      if(transdep[eq] == 1){
        if(tra[tr].like == 1) Lr += log(ratecalc(eq,popnum,param));

        d = transdepref[eq];
        ed = depeqdiv[d][long(onefac*t*depeqdivmax/tmax)];
        do{
          ee = depeq_evnext[ed]
              ; if(ee == -1) emsg("neg ee");
          if(ee == -1 || depeq_evt[ee] >= t) break;
          ed = ee;
        }while(1 == 1);

        for(k = 0; k < ndepeq_ev_ev[ed]; k++) if(depeq_ev_evt[ed][k] == t) break;
        if(k == ndepeq_ev_ev[ed]){
          if(tra[tr].like == 1) emsg("cannot find");
        }
        else{
          if(tra[tr].like == 0) emsg("this should not be here");

          dd = ratecalc(eq,popnum,param)-depeq_evval[ed];
          if(dd*dd > tiny){
            cout << " " << depeq_ev_evt[ed][k] << " " << t << " t\n";
            cout << ratecalc(eq,popnum,param) << " " << depeq_evval[ed] << " " << ed << " fz\n";
            cout << t << " t\n";
            emsg("problem with rate h");
          }
          if(depeq_ev_evflag[ed][k] == 1) emsg("duplicate");
          depeq_ev_evflag[ed][k] = 1;
        }
      }
      else{
        //cout << eq << " " <<  eqnstr[eq] << "eqn\n"; nnn++;
        if(tra[tr].like == 1) Lr += log(ratecalc(eq,popnum,param));
      }
    }

    ci = stat[i];
    if(tra[tr].ci != ci) emsg("not consistent",num,samp);
    cf = tra[tr].cf; stat[i] = cf;
    ncompind[ci]--; if(ncompind[ci] < 0) emsg("neg pop ch");
    ncompind[cf]++;
  }

  if(capevfl == 0){
    if(Lev_pd != 0) emsg("Lev_pd prob"); 
  }
  else{
    for(ce = 0; ce < ncapev; ce++){
      if(cev_num[ce] != capev_num[ce]) emsg("capevnum",cev_num[ce],capev_num[ce]);
      if(cev_oneminus[ce] != capev_oneminus[ce]) emsg("capevoneminus");
      pd = calcobsprob(capevprobeqn[ce]);
      if(capev_val[ce] != log(pd)){ cout <<  capev_val[ce] << " "<< log(pd) << "val\n"; emsg("capev val");}
      if(capev_valoneminus[ce] != log(1-pd)) emsg("capev valoneminus");

      if(capev_num[ce] > 0){
        if(pd == 0) emsg("Cannot have zero transition capture probability");
        Le_pd += capev_num[ce]*log(pd);
      }
      if(capev_oneminus[ce] > 0){
        if(pd == 1) emsg("Cannot have transition capture probability of one");
        Le_pd += capev_oneminus[ce]*log(1-pd);
      }
    }
  }

  //cout << nnn << " " << tott << " nnn\n";
  for(d = 0; d < ntransdepeq; d++){
    e = depeqdiv[d][0];
    eq = transdepeq[d];

    do{
      for(k = 0; k < ndepeq_ev_ev[e]; k++){
        if(depeq_ev_evflag[e][k] == 0) emsg("event not accounted for");
      }
      e = depeq_evnext[e];
    }while(e != -1);
  }

  Lpr = priorcalc();

  nwrong = 0;
  Lo_st = 0;
  for(i = 0; i < nind; i++){
    if(i < nindtot){
      c = NOTALIVE; e = 0;
      for(j = 0; j < nindobs[i]; j++){
        ob = indobs[i][j];

        t = obst[ob]; //cout << t << " t\n";
        while(e < nindev[i] && indev[i][e].t < t){ if(tra[indev[i][e].tr].type >= 0) c = tra[indev[i][e].tr].cf; e++;}
        if(c == NOTALIVE){ Lo_st += notobsdL; nwrong++;}
        else{
          eq = obsprobeqn[ob][c];
          if(eq == -1){ Lo_st += notobsdL; //cout << i << ","; 
            nwrong++;}
          else Lo_st += log(calcobsprob(eq));
        }
      }
    }
  }
  //cout << "  " << nwrong << "nwrong\n";

  Lpopch = 0;                                                                 // checks popm
  for(p = 0; p < npopm; p++){ 
    short numm = 0;
    t = popmt[p];
    for(i = 0; i < nindtot; i++){
      c = NOTALIVE; e = 0; while(e < nindev[i] && indev[i][e].t < t){ c = tra[indev[i][e].tr].cf; e++;}
      if(popmcomp[p][c] == 1) numm++;
    }
    if(numm != pop[p]) emsg("pop not right");

    dd = popL[p]-popcalc(p,numm); if(dd > tiny) emsg("popL prob");
    Lpopch += popcalc(p,numm);
  }

  for(k = 0; k < nnmeq; k++){
    eq = nmeq[k];
    if(nmeq_val[eq] != calculatenotdep(eq,param)){
     emsg("nmeq1",nmeq_val[eq],calculatenotdep(eq,param),k);
    }
  }

  Lnm = 0; for(i = 0; i < nindtot; i++) Lnm += likenm(i);             // Caluculates the likelihood for the non Markovian transitions

  nindinitch = 0;
  for(c = 0; c < ncompswa; c++) ncompinitch[c] = 0;
  for(i = 0; i < nindtot; i++){ if(indev[i][0].t == 0){ nindinitch++; ncompinitch[tra[indev[i][0].tr].cf%ncompswa]++;}}
  for(c = 0; c < ncompswa; c++){ if(ncompinitch[c] != ncompinit[c]){ cout << ncompinitch[c] << " " << ncompinit[c] << " ch\n"; emsg("proo");}}
  if(nindinitch != nindinit) emsg("indinit");

  Linit = logsum[nindinit]; for(c = 0; c < ncompswa; c++) Linit += ncompinit[c]*log(probinit[c]);

  dd = Liexp-Lexp; if(dd*dd > tiny){ cout << Liexp << "   true:" << Lexp << "\n"; emsg("Liexp prob");}
  dd = Lir-Lr; if(dd*dd > 0.00001) emsg("Lir propb");
  dd = Linm-Lnm; if(dd*dd > 0.0001){ cout << num << " " <<  Linm << " True:" << Lnm  << "\n"; emsg("Lnm prob");}
  Linm = Lnm;

  dd = Liinit-Linit; if(dd*dd > tiny){ cout << Liinit << " "<< Linit << "\n"; emsg("Linit prob");}
  dd = Lev_pd-Le_pd; if(dd*dd > tiny){ cout << Lev_pd << " " << Le_pd << "\n"; emsg("Lev_pd prob");}
  dd = Lob_pd-Lo_pd; if(dd*dd > tiny){ cout << Lob_pd << " " << Lo_pd << "\n"; emsg("Lob_pd prob");}
  dd = Lob_st-Lo_st; if(dd*dd > 0.0001){ cout << Lob_st << " " << Lo_st << " jj\n"; emsg("Lob_st prob");}
  dd = Lpop-Lpopch; if(dd*dd > 0.0001){ cout << Lpop << " " << Lpopch << " jj\n"; emsg("Lpop prob");}
  dd = Lpri-Lpr; if(dd*dd > tiny){ cout << num << " " << Lpri << " " << Lpr << "\n"; emsg("Lpri prob");}
}

void checkdist()
{
  const long nbin = 100;
  double bin[nbin];
  long bi;
  long j;
  double max = 20,mean = 2, kshape = 5, lam= 10;
  for(bi = 0; bi < nbin; bi++) bin[bi] = 0;
 
  for(j = 0; j < 100000; j++){
    // bi = long(double(gammasamp(kshape,kshape/mean)*nbin)/max);
    // bi = long(double(gammasamptail(kshape,kshape/mean,4)*nbin)/max);
    bi = long(double(weibullsamp(lam,kshape)*nbin)/max);
    // bi = long(double(weibullsamptail(lam,kshape,3.5)*nbin)/max);

    if(bi < nbin) bin[bi]++;
  }
  ofstream dis("GP/dis");
  cout << nbin << "h\n";
  for(bi = 0; bi < nbin; bi++) dis << bi*max/nbin << " " << bin[bi] << "\n";

  //gammasamptail(double a, double b, double tmin)   
}

void Chain::checkderive()                                             // CHecks that the derived likelihood is correct
{
  long d, dm, i, j, k, e, tr, p, c;
  double t, popnum[npopnum], sum, L, LL, dd;
  vector <EV> ev;
  vector <long> nc;

  for(i = 0; i < nindtot; i++){
    for(e = 0; e < nindev[i]; e++) ev.push_back(indev[i][e]);
  }
  sort(ev.begin(),ev.end(),compareev);

  nc.resize(ncomp+1); for(c = 0; c < ncomp; c++) nc[c] = 0; nc[NOTALIVE] = nindtot_sim;

  L = 0; j = 0;
  for(dm = 0; dm < nderm; dm++){
    t = dermt[dm];
    while(j < ev.size() && ev[j].t <= t){
      tr = ev[j].tr; nc[tra[tr].ci]--; nc[tra[tr].cf]++;
      j++;
    }

    for(p = 0; p < npopnum; p++){
      sum = 0; for(k = 0; k < popnumterm[p].size(); k++) sum += popnumtermweight[p][k]*nc[popnumterm[p][k]];
      popnum[p] = sum;
    }

    d = derm[dm];
    LL = dercalc(dm,popnum,param);
    dd = LL-derL[dm]; if(dd*dd > tiny){ cout << dm << " " << LL << " " << derL[dm]<< "\n"; emsg("derL prob");}
    L += LL;
  }
  dd = L-Lder; if(dd*dd > tiny){ cout << L << " " << Lder << " " <<  "\n"; emsg("Lder prob");}

  for(dm = 0; dm < nderm; dm++){ if(dermap[dm] != 0) emsg("dermapprob");}
}

void Chain::checkparamgrad()   // Checks that paramgrad is working correctly
{
  long p;
  double val, dd, Lst, Lup, Ldo, d, pmean2, pvar2, pgrad, pcurve;

  Lst = L();
  for(p = 0; p < nparam; p++){
    val = param[p]; dd = val/1000; if(dd == 0) dd = 0.0000001;
    paramgrad(p);
    changeparam(1,p,val+dd); Lup = L();
    changeparam(1,p,val-dd); Ldo = L();

    pgrad = (Lup-Ldo)/(2*dd);
    pcurve = (Lup+Ldo-2*Lst)/(dd*dd);
    cout << Lst << " " << Lup << " " << Ldo << " " <<  pgrad << " " << pcurve << " j\n";

    if(pcurve > 0){ pmean2 = 0; pvar2 = 0;}
    else{
      pmean2 = val - pgrad/pcurve;
      pvar2 = -1.0/pcurve;
      if(pmean2*pmean2 > 100*val*val){ pmean2 = 0; pvar2 = 0;}
    }
    d = (pmean2 - pmean)/val; 
    if(d*d > 0.0001){
      cout << "z|";
      cout << pmean2 << " " << pgrad << " " << pcurve <<" " << pgrad/pcurve << "ppp\n";
      cout << " " << val << " " << Lup << " " << Ldo << " " << Lst << " yy\n";
      emsg("paramgrad prob");
    }
    d = (pvar2 - pvar)/(val*val); 
    if(d*d > 0.0001){ 
      cout << Lup << " " << Ldo << " " << Lst << " " << Lup -Lst <<" " << Ldo-Lst << " L\n";
      cout << pvar2 << " " << pvar << " " << paramname[p] << " " <<  param[p] << " cop\n"; emsg("paramgrad prob2");
    }

    changeparam(1,p,val);
  }
}




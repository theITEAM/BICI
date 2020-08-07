vector< vector< vector<short> > > nthreeev;
vector< vector< vector< vector<short> > > > threeev1, threeev2;                       // [cl][i][f]  This gives all the combination of three events starting at i and ending f

vector< vector< vector<short> > > ntwoev;
vector< vector< vector< vector<short> > > > twoev;                          // [cl][i][f]  This gives all the combination of two events starting at i and ending f


vector <short> addfl;
short ntrcl;
vector <long> trcl;
vector <short> trcle;
vector <EV> evrem, evadd;


void Chain::movegroup_prop(long i)                                        // moves a group of connected events in individual i
{
  short e, nev, g, gg, cl, ee, k, kmax, tran;
  long tr, trr;
  double tmi, tma, dt, t, al, Ltoti, Ltotf;

  /*
  nev = nindev[i];
  if(nev <= 3) return;

  e = 1 + short(ran()*(nev-3));
  tr = indev[i][e].tr;
  cl = tra[tr].cl;

  addfl.resize(nev); for(ee = 1; ee < nev-1; ee++) addfl[ee] = 0; addfl[e] = 1;

  g = -1;   // Makes timeline just in cl
  ntrcl = 0; trcl.clear(); trcle.clear();
  for(ee = 1; ee < nev-1; ee++){ tr = indev[i][ee].tr; if(tra[tr].cl == cl){ trcl.push_back(tr); trcle.push_back(ee); if(ee == e) g = ntrcl; ntrcl++;}}
  if(g == -1) emsg("g");

  gg = g+1;
  while(gg < ntrcl){
    trr = trcl[gg-1]; kmax =nnonexptra[trr]; k = 0; while(k < kmax && tra[nonexptra[trr][k]].dc != tra[trcl[gg]].dc) k++;
    if(k < kmax){ addfl[trcle[gg]] = 1; gg++;}
    else break;
  }
  if(gg < ntrcl) tma = indev[i][trcle[gg]].t; else tma = indev[i][nev-1].t;

  gg = g-1;
  while(gg >= 0){
    trr = trcl[gg]; kmax =nnonexptra[trr]; k = 0; while(k < kmax && tra[nonexptra[trr][k]].dc != tra[trcl[gg+1]].dc) k++;
    if(k < kmax){ addfl[trcle[gg]] = 1; g--;}
    else break;
  }
  if(gg >= 0) tmi = indev[i][trcle[gg]].t; else tmi = indev[i][0].t;

  ee = 1; while(ee < nev && addfl[ee] != 1) ee++; tran = tra[indev[i][ee].tr].trans;
  dt = normal(0,jump_movegroup[i][tran]);

  evnew = indev[i];
  for(ee = 1; ee < nev-1; ee++) if(addfl[ee] == 1){ t = evnew[ee].t +dt; if(t < tmi || t > tma) break; evnew[ee].t = t;}

  if(ee < nev-1) al = -1;
  else{
    sort(evnew.begin(),evnew.end(),compareev);
    //zz TOD

    evrev = indev[i];
    Ltoti = Lpri + (Liinit + Lir - Liexp + Linm)*invTLi + Lob_pd*invT_pd + Lob_st*invT_st + Lpop*invT_pop;
    indchange(i);
    Ltotf = Lpri + (Liinit + Lir - Liexp + Linm)*invTLi + Lob_pd*invT_pd + Lob_st*invT_st + Lpop*invT_pop;

    al = exp(Ltotf - Ltoti); //cout << al << " " << Ltoti << " " << Ltotf << " " << probif << " " << probfi << "al\n";
  }

  ntr_movegroup[i]++;
  if(ran() < al){
    nac_movegroup[i]++;

    if(samp < burnin){ jump_movegroup[i][tran] *= 1.05;}
  }
  else{
    if(al >= 0){ evnew = evrev; indchange(i);}

    if(samp < burnin){ jump_movegroup[i][tran] *= 0.975;}
  }
  */
}

void Chain::singevent_prop(short i)    // looks at inserting or removing single events
{
  short e, ee, eee, nev, cl, ii, ff, mm, c, ntree, ci, cm, cmm, cf, dc1, dc2, ntwo, nthree, sel, num, ei, ef, cm2, mm2, mmsel;
  long tr, tr2, tr3, trnew;
  double al, ti, tf, Ltoti, Ltotf, tnew;

  nev = nindev[i];
  e = short(ran()*(nev-1));
  tr = indev[i][e].tr; cl = tra[tr].cl;
  if(cl >= nclass-3 || tra[tr].cf == NOTALIVE) return;

  if(ran() < 0.5){  // Tries to add a single event
    if(cl == -1){ cl = short(ran()*(nclass-3)); num = nclass-3;} else num = 1;

    for(ee = e+1; ee < nev; ee++){ tr2 = indev[i][ee].tr; if(tra[tr2].cl == cl || tra[tr2].cl < 0) break;}
    if(ee == nev) emsg("nevprob");

    ci = tra[tr].ci; if(ci == NOTALIVE) ii = nclassval[cl]; else ii = compval[ci][cl];
    cm = tra[tr].cf; mm = compval[cm][cl];
    cf = tra[tr2].cf; if(cf == NOTALIVE) ff = nclassval[cl]; else ff = compval[cf][cl];

    nthree = nthreeev[cl][ii][ff];
    if(nthree > 0){
      ntr_addsing[i]++;

      //evrev = indev[i];
      //evnew = evrev;

      evnew.clear();

      sel = short(ran()*nthree);
      dc1 = (threeev1[cl][ii][ff][sel]-mm)*classmult[cl];
      dc2 = (threeev2[cl][ii][ff][sel]-mm)*classmult[cl];

      ti = evnew[e].t; tf = evnew[ee].t; if(tf == ti) emsg("pr");
      tnew = ti + ran()*(tf-ti);

      cm += dc1;
      if(ti == 0) trnew = trabeg+cm; else trnew = compiftra[ci][cm];
      if(trnew < 0 || tra[tr].fix != tra[trnew].fix){ nfa_addsing[i]++; return;}
      EV ev; ev.tr = trnew; ev.t = ti; evnew.push_back(ev); if(tra[tr].nm != tra[trnew].nm) emsg("nm cha");
      //}

      eee = e+1;
      while(eee < ee && evnew[eee].t < tnew){
        cmm = tra[evnew[eee].tr].cf+dc1; trnew = compiftra[cm][cmm]; cm = cmm; if(trnew < 0){ nfa_addsing[i]++; return;};
        EV ev; ev.tr = trnew; ev.t = evnew[eee].t; evnew.push_back(ev);
        evnew[eee].tr = trnew;
        eee++;
      }

      cmm = cm + dc2 - dc1; trnew = compiftra[cm][cmm]; cm = cmm;
      if(trnew < 0 || tra[trnew].fix == 1){ nfa_addsing[i]++; return;}
      EV ev; ev.t = tnew; ev.tr = trnew; evnew.insert(evnew.begin()+eee,ev);
      eee++;

      if(dc2 != 0){
        while(eee <= ee && evnew[eee].t < tnew){ 
          cmm = tra[evnew[eee].tr].cf+dc2; trnew = compiftra[cm][cmm]; cm = cmm; if(trnew < 0){ nfa_addsing[i]++; return;};
          evnew[eee].tr = trnew;
          eee++;
        }

        if(tf == tmax) trnew = traend+cm; else{ trnew = compiftra[cm][cf]; if(trnew < 0){ nfa_addsing[i]++; return;}}
        if(fixfl == 1){ if(tra[tr2].fix != tra[trnew].fix == 0){ nfa_addsing[i]++; return;}}
        evnew[ee+1].tr = trnew; if(tra[trnew].nm == 1) emsg("prob");
      }

      Ltoti = L(); oe("bef",indev[i]); oe("add",evnew); cout << ti << " " << tf << " tif\n";
      indchange(i,e,ti,tf);
      Ltotf = L(); oe("aft",indev[i]);

      al = exp(Ltotf-Ltoti)*(nthree*(tf-ti)*num*(nev-1))/(nev*ntwoev[cl][ii][ff]);//cout << al << " " << Ltotf << " "<< Ltoti <<   "al add \n";

      if(ran() < al){
        nac_addsing[i]++;
      }
      else{ evnew = evrev; indchange(i,e,ti,tf);}
    }
  }
  else{    // Tries to remove a single event
    if(cl < 0 ||  tra[tr].nm == 1 || tra[tr].fix == 1 || nnonexptra[tr] > 0) return; if(e == 0 || e == nev-1) emsg("should not be");

    for(ef = e+1; ef < nev; ef++){ tr2 = indev[i][ef].tr; if(tra[tr2].cl == cl || tra[tr2].cl < 0) break;}
    if(ef == nev) emsg("nevprob");

    for(ei = e-1; ei >= 0; ei--){ tr3 = indev[i][ei].tr; if(tra[tr3].cl == cl || tra[tr3].cl < 0) break;}
    if(ei < 0) emsg("nevprob");
    if(ei == 0) num = nclass-3; else num = 1;

    ntr_remsing[i]++;

    ci = tra[tr3].ci; if(ci == NOTALIVE) ii = nclassval[cl]; else ii = compval[ci][cl];
    cm = tra[tr3].cf; mm = compval[cm][cl];
    cm2 = tra[tr].cf; mm2 = compval[cm2][cl];
    cf = tra[tr2].cf; if(cf == NOTALIVE) ff = nclassval[cl]; else ff = compval[cf][cl];

    ntwo = ntwoev[cl][ii][ff];
    if(ntwo > 0){
      //evrev = indev[i];
      evnew.clear();

      ti = evnew[ei].t; tf = evnew[ef].t;

      sel = short(ran()*ntwo); mmsel = twoev[cl][ii][ff][sel];
      dc1 = (mmsel-mm)*classmult[cl];
      dc2 = (mmsel-mm2)*classmult[cl];

      if(dc1 != 0){
        cm += dc1;

        if(nnonexptra[tr3] > 0){ nfa_remsing[i]++; return;};
        if(ti == 0) trnew = trabeg + cm; else trnew = compiftra[ci][cm]; if(trnew < 0 || nnonexptra[trnew] > 0){ nfa_remsing[i]++; return;}
        if(fixfl == 1){ if(tra[tr3].fix == 1 && tra[trnew].fix == 0){ nfa_remsing[i]++; return;}}
        evnew[ei].tr = trnew;

        for(eee = ei+1; eee < e; eee++){
          cmm = tra[evnew[eee].tr].cf+dc1; trnew = compiftra[cm][cmm]; cm = cmm; if(trnew < 0){ nfa_remsing[i]++; return;};
          evnew[eee].tr = trnew;
        }
      }

      if(dc2 != 0){
        for(eee = e+1; eee < ef; eee++){
          cmm = tra[evnew[eee].tr].cf+dc2; trnew = compiftra[cm][cmm]; cm = cmm; if(trnew < 0){ nfa_remsing[i]++; return;};
          evnew[eee].tr = trnew;
        }

        if(tf == tmax) trnew = traend + cm; else{ trnew = compiftra[cm][cf]; if(trnew < 0){ nfa_remsing[i]++; return;}}
        if(fixfl == 1){ if(tra[tr2].fix == 1 && tra[trnew].fix == 0){ nfa_remsing[i]++; return;}}
        evnew[ef].tr = trnew;
      }

      evnew.erase(evnew.begin()+e);

      //outputevents("rem",indev[i]); outputevents("rem2",evnew);
      //oprop();

      Ltoti = L();
      indchange(i,ei,ti,tf);
      Ltotf = L();

      al = exp(Ltotf-Ltoti)*(ntwo*(nev-1))/(nthreeev[cl][ii][ff]*(tf-ti)*num*(nev-2)); //cout << al <<  "al\n";

      if(ran() < al){
        nac_remsing[i]++;
      }
      else{ evnew = evrev; indchange(i,ei,ti,tf);}
    }
  }
  //cout << "don\n";
}

void eventpropinit()    // This initialise tra.trans such that jumping size for events can be adjusted for each trans
{
  short j, cl, nval, ii, iii, iiii, iiiii, jj, jjj, jjjj;
  long tr;
  vector< vector<short> > transleave;

  ntrans = 0;
  for(tr = 0; tr < ntra; tr++){
    cl = tra[tr].cl;
    j = 0; while(j < ntrans && (trans[j].cl != cl || trans[j].i != tra[tr].i || trans[j].f != tra[tr].f)) j++;
    tra[tr].trans = j;
    if(j == ntrans){ TRANS tran; tran.cl = cl; tran.i = tra[tr].i; tran.f = tra[tr].f; trans.push_back(tran); ntrans++;}
  }

  nthreeev.resize(nclass-3); threeev1.resize(nclass-3); threeev2.resize(nclass-3); ntwoev.resize(nclass-3); twoev.resize(nclass-3);
  for(cl = 0; cl < nclass-3; cl++){
    nval = nclassval[cl];

    transleave.clear();
    transleave.resize(nval+1);
    for(j = 0; j < ntrans; j++){ if(trans[j].cl == cl){ transleave[trans[j].i].push_back(trans[j].f);}}
    for(ii = 0; ii < nval; ii++){ transleave[nval].push_back(ii); transleave[ii].push_back(nval);}

    nthreeev[cl].resize(nval+1); threeev1[cl].resize(nval+1); threeev2[cl].resize(nval+1);
    for(ii = 0; ii <= nval; ii++){ 
      nthreeev[cl][ii].resize(nval+1); for(iii = 0; iii <= nval; iii++) nthreeev[cl][ii][iii] = 0;
      threeev1[cl][ii].resize(nval+1); threeev2[cl][ii].resize(nval+1);
    }

    for(ii = 0; ii <= nval; ii++){
      for(j = 0; j < transleave[ii].size(); j++){
        iii = transleave[ii][j];
        if(iii != nval){
          for(jj = 0; jj < transleave[iii].size(); jj++){
            iiii = transleave[iii][jj];
            if(iiii != nval){
              for(jjj = 0; jjj < transleave[iiii].size(); jjj++){
                iiiii = transleave[iiii][jjj];
                threeev1[cl][ii][iiiii].push_back(iii);
                threeev2[cl][ii][iiiii].push_back(iiii);
                nthreeev[cl][ii][iiiii]++;
              }
            }
          }
        }
      }
    }

    ntwoev[cl].resize(nval+1); twoev[cl].resize(nval+1);
    for(ii = 0; ii <= nval; ii++){ 
      ntwoev[cl][ii].resize(nval+1); for(iii = 0; iii <= nval; iii++) ntwoev[cl][ii][iii] = 0;
      twoev[cl][ii].resize(nval+1);
    }

    for(ii = 0; ii <= nval; ii++){
      for(j = 0; j < transleave[ii].size(); j++){
        iii = transleave[ii][j];
        if(iii != nval){
          for(jj = 0; jj < transleave[iii].size(); jj++){
            iiii = transleave[iii][jj];
            twoev[cl][ii][iiii].push_back(iii);
            ntwoev[cl][ii][iiii]++;
          }
        }
      }
    }

    if(1 == 0){
      for(ii = 0; ii <= nval; ii++){
        for(iii = 0; iii <= nval; iii++){
          cout << ii << " " << iii << ": ";
          for(j = 0; j < nthreeev[cl][ii][iii]; j++) cout << threeev1[cl][ii][iii][j] << "," << threeev2[cl][ii][iii][j] << "   ";
          cout << "  threeev\n";
        }
      }

      for(ii = 0; ii <= nval; ii++){
        for(iii = 0; iii <= nval; iii++){
          cout << ii << " " << iii << ": "; for(j = 0; j < ntwoev[cl][ii][iii]; j++) cout << twoev[cl][ii][iii][j] << "   ";
          cout << "  twoev\n";
        }
      }
    }
  }
}

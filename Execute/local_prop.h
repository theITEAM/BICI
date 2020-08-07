// makes local changes to the event sequence

double prswitch = 0.7;

vector< vector< vector< vector< vector <long> > > > > locpath;
vector< vector< vector< vector <long> > > > nlocpath;

vector< vector <long> > locpathst;
vector <long> nlocpathst;

void Chain::sing_prop(long i)
{
  long e, ee, nev, cl, res;

  nev = nindev[i];

  e = long(ran()*(nev-1));
  if(e == 0) cl = long(ran()*(nclass-3)); else{ cl = tra[indev[i][e].tr].cl; if(cl >= nclass-3) return; if(cl < 0) emsg("uu");}

  ee = e+1; while(ee < nev && tra[indev[i][ee].tr].cl != cl && tra[indev[i][ee].tr].cl >= 0) ee++;
  if(ran() < 0.5){
  }
  else{
    ee++; if(ee == nev) return;
    while(ee < nev && tra[indev[i][ee].tr].cl != cl && tra[indev[i][ee].tr].cl >= 0) ee++;
  }

  ntr_sing[i][cl]++;
  if(ran() < 0.5){  // Fixes right event
    res = localresamp(i,cl,e+1,indev[i][e].t+tiny,indev[i][ee].t,3,1);
    switch(res){
      case 1: nac_sing[i][cl]++; break;
      case -1: nfa_sing[i][cl]++; break;
    }
  }
  else{ // Fixes left event
    res = localresamp(i,cl,e,indev[i][e].t,indev[i][ee].t-tiny,2,1);
    switch(res){
      case 1: nac_sing[i][cl]++; break;
      case -1: nfa_sing[i][cl]++;  break;
    }
  }
}

void Chain::twothree_prop(long i)                                        // Interchanges between two and three events
{
  long e, ee, nev, cl;

  nev = nindev[i];
  e = long(ran()*(nev-1));
  if(e == 0) cl = long(ran()*(nclass-3)); else{ cl = tra[indev[i][e].tr].cl; if(cl >= nclass-3) return; if(cl < 0) emsg("uu");}

  ee = e+1; while(ee < nev && tra[indev[i][ee].tr].cl != cl && tra[indev[i][ee].tr].cl >= 0) ee++;
  if(ran() < 0.5){
  }
  else{
    ee++; if(ee == nev) return;
    while(ee < nev && tra[indev[i][ee].tr].cl != cl && tra[indev[i][ee].tr].cl >= 0) ee++;
  }

  ntr_twothree[i][cl]++;
  switch(localresamp(i,cl,e,indev[i][e].t,indev[i][ee].t,1,1)){
    case 1: nac_twothree[i][cl]++; break;
    case -1: nfa_twothree[i][cl]++; break;
  }
}

void Chain::pair_prop(long i)
{
  long e, ee, nev, cl, fl;

  nev = nindev[i];
  e = long(ran()*(nev-1));
  if(e == 0) cl = long(ran()*(nclass-3)); else{ cl = tra[indev[i][e].tr].cl; if(cl >= nclass-3) return;}

  ee = e+1; while(ee < nev && tra[indev[i][ee].tr].cl != cl && tra[indev[i][ee].tr].cl >= 0) ee++;
  if(ran() < 0.5){
    fl = 0;
  }
  else{
    ee++; if(ee == nev) return;
    while(ee < nev && tra[indev[i][ee].tr].cl != cl && tra[indev[i][ee].tr].cl >= 0) ee++;
    ee++; if(ee == nev) return;
    while(ee < nev && tra[indev[i][ee].tr].cl != cl && tra[indev[i][ee].tr].cl >= 0) ee++;
    if(compval[tra[indev[i][e].tr].cf][cl] != compval[tra[indev[i][ee].tr].ci][cl]) return;
    fl = 1;
  }

  ntr_pair[i][cl]++;
  switch(localresamp(i,cl,e+1,indev[i][e].t+tiny,indev[i][ee].t-tiny,0,1)){
    case 1: nac_pair[i][cl]++; if(samp < burnin) pairac[cl][fl]++; break;
    case -1: nfa_pair[i][cl]++; break;
  }
}

void Chain::gap_prop(long i)                  // Resample event sequences between the gaps in obsrvations
{
  long e, nev, cl, o, ty;
  double tmi, tma;

  getlifespan(i);

  e = 0;
  cl = long(ran()*(nclass-3));
  for(o = -1; o < nindobs[i]; o++){
    ty = 0;
    if(o == -1){ ty = 2; tmi = tent;} else tmi = indobst[i][o];
    if(o < nindobs[i]-1) tma = indobst[i][o+1]; else{ ty = 3; tma = tlea;}

    nev = nindev[i];
    while(e < nev && indev[i][e].t < tmi) e++;

    ntr_gap[i][cl]++;
    switch(localresamp(i,cl,e,tmi,tma,ty,0)){
      case 1: nac_gap[i][cl]++;  break;
      case -1: nfa_gap[i][cl]++; break;
    }
  }
}

long Chain::localresamp(long i, long cl, long estart, double t, double tend, long type, long evfac)  // estart is first event in time window
{
  long e, nval, c, eend, ii, ff, np, ci, cf, nev, nevadd, ji, jf, k, numi, n, fl, nevnew;
  long tr, tr2, p;
  double z, tt, Ltoti, Ltotf, al, fac = 1, t1, t2, addtime[3];
  vector <long> inter;

  //oe("seq",indev[i]); cout << estart << " " << t << " " << tend << " " << type << " " << cl << "uu\n";

  if(tend <= t){ if(t - tend > 2*tiny) emsg("t wrong2"); return 0;}

  nval = nclassval[cl];

  tr = indev[i][estart].tr;
  ci = tra[tr].ci;
  if(ci == NOTALIVE) ii = nval; else ii = compval[ci][cl];
  c = ci;

  nev = nindev[i];
  eend = estart;

  if(fixfl == 1 && i < nind){   // checks no fixed events in region
    e = 0; while(e < nindfixev[i] && fixev[indfixev[i][e]].t <= t) e++;
    while(e < nindfixev[i] && fixev[indfixev[i][e]].t < tend){
      if(fixev[indfixev[i][e]].cl == cl) return -1;
      e++;
    }
  }

  fl = 0;
  while(eend < nev && indev[i][eend].t <= tend){
    tr = indev[i][eend].tr; //if(fixfl == 1){ if(tra[tr].fix == 1 && indev[i][eend].t > t &&  indev[i][eend].t < tend) return -1;}
    if(tra[tr].cl == cl || tra[tr].cl < 0){ if(fl == 0) fl = 1; else inter.push_back(compval[tra[tr].ci][cl]);}
    eend++; c = tra[tr].cf;
  }
  if(c == NOTALIVE) ff = nval; else ff = compval[c][cl];

  numi = long(inter.size())+1; if(numi == 1 && ii == ff) numi = 0; if(numi > 3) return -1;

  /*
  if(checkon == 1){
    for(ji = 0; ji < nlocpath[cl][ii][ff][numi]; ji++){
      p = locpath[cl][ii][ff][numi][ji];
      for(k = 0; k < nlocpathst[p]; k++) if(locpathst[p][k] != inter[k]) break;
      if(k == nlocpathst[p]) break;
    }
    if(ji == nlocpath[cl][ii][ff][numi]){ for(k = 0; k < inter.size(); k++) cout << inter[k] <<","; cout << ii << "  "<< ff << "hh\n"; emsg("cannot find");}
  }
  */
  /*
  if(type == 2 && cl == 1){
    oe("ind",indev[i]);
    cout << indid[i] << " " << cl << " " << estart << " " << t << " " << tend << " yy\n";
  }
  */

  p = -1;
  switch(type){
    case 1:   // when there are events at the begining and end
      switch(numi){
        case 2:
          if(ran() < prswitch){  // 2 to 3 transition
           n = nlocpath[cl][ii][ff][3];
            if(n == 0) return -1; if(nlocpath[cl][ii][ff][2] == 0) emsg("p");
            fac *= double(n)/nlocpath[cl][ii][ff][2];

            p = locpath[cl][ii][ff][3][long(ran()*n)];

            nevadd = 3; addtime[0] = t; addtime[1] = t+ran()*(tend-t); addtime[2] = tend;
            fac *= (tend-t);
          }
          else{                // 2 to 2 transition
            n = nlocpath[cl][ii][ff][2]; if(n <= 1) return -1;
            p = locpath[cl][ii][ff][2][long(ran()*n)];

            nevadd = 2; addtime[0] = t; addtime[1] = tend;
          }
          break;

        case 3:
          if(ran() < prswitch){  // 3 to 2 transition
            n = nlocpath[cl][ii][ff][2]; if(n == 0) return -1; if(nlocpath[cl][ii][ff][3] == 0) emsg("p");
            fac *= double(n)/nlocpath[cl][ii][ff][3];
            p = locpath[cl][ii][ff][2][long(ran()*n)];

            nevadd = 2; addtime[0] = t; addtime[1] = tend;
            fac /= (tend-t);
          }
          else{                // 3 to 3 transition
            n = nlocpath[cl][ii][ff][3]; if(n <= 1) return -1;
            p = locpath[cl][ii][ff][3][long(ran()*n)];

            nevadd = 3; addtime[0] = t; addtime[1] = t+ran()*(tend-t); addtime[2] = tend;
          }
          break;

          default: emsg("pp1"); return -1;
      }
      break;

    case 0:
      switch(numi){
        case 0:
          if(ii != ff) emsg("proby");

          n = nlocpath[cl][ii][ff][2]; if(n == 0) return -1; if(nlocpath[cl][ii][ff][0] != 1) emsg("P");
          fac *= n;
          p = locpath[cl][ii][ff][2][long(ran()*n)];

          nevadd = 2;
          t1 = t+ran()*(tend-t); t2 = t+ran()*(tend-t);
          if(t1 < t2){ addtime[0] = t1; addtime[1] = t2;} else{ addtime[0] = t2; addtime[1] = t1;}
          fac *= 0.5*(tend-t)*(tend-t);
          break;

        case 1:
          n = nlocpath[cl][ii][ff][2]; if(n == 0) return -1; if(nlocpath[cl][ii][ff][1] != 1) emsg("P");
          fac *= n;
          p = locpath[cl][ii][ff][2][long(ran()*n)];

          nevadd = 2;
          t1 = t+ran()*(tend-t); t2 = t+ran()*(tend-t);
          if(t1 < t2){ addtime[0] = t1; addtime[1] = t2;} else{ addtime[0] = t2; addtime[1] = t1;}
          fac *= 0.5*(tend-t);
          break;

        case 2:
          if(ii == ff){
            n = nlocpath[cl][ii][ff][0]; if(n != 1) emsg("pp2");
            fac /= nlocpath[cl][ii][ff][2];
            p = locpath[cl][ii][ff][0][0];

            nevadd = 0;
            fac /= 0.5*(tend-t)*(tend-t);
          }
          else{
            n = nlocpath[cl][ii][ff][1]; if(n == 0) return -1; if(n != 1) emsg("PP");
            fac /= nlocpath[cl][ii][ff][2];
            p = locpath[cl][ii][ff][1][0];

            nevadd = 1; addtime[0] = t+ran()*(tend-t);
            fac /= 0.5*(tend-t);
          }
          break;

        default: return -1;
      }
      break;

    case 2:  // includes the begin event
      switch(numi){
        case 1:
          n = nlocpath[cl][ii][ff][2]; if(n == 0) return -1; if(nlocpath[cl][ii][ff][1] != 1) emsg("P");
          fac *= n;
          p = locpath[cl][ii][ff][2][long(ran()*n)];

          nevadd = 2; addtime[0] = t; addtime[1] = t+ran()*(tend-t);
          fac *= (tend-t);
          break;

        case 2:
          n = nlocpath[cl][ii][ff][1]; if(n == 0) return -1; if(n != 1) emsg("PP");
          fac /= nlocpath[cl][ii][ff][2];
          p = locpath[cl][ii][ff][1][0];

          nevadd = 1; addtime[0] = t;
          fac /= (tend-t);
          break;

          default:return -1;
      }
      break;

   case 3:  // includes the end event
      //if(ii != nval) emsg("val");
      switch(numi){
        case 1:
          n = nlocpath[cl][ii][ff][2]; if(n == 0) return -1; if(nlocpath[cl][ii][ff][1] != 1) emsg("P");
          fac *= n;
          p = locpath[cl][ii][ff][2][long(ran()*n)];

          nevadd = 2; addtime[0] = t+ran()*(tend-t); addtime[1] = tend;
          fac *= (tend-t);
          break;

        case 2:
          n = nlocpath[cl][ii][ff][1]; if(n == 0) return -1; if(n != 1) emsg("PP");
          fac /= nlocpath[cl][ii][ff][2];
          p = locpath[cl][ii][ff][1][0];
          nevadd = 1; addtime[0] = tend;
          fac /= (tend-t);
          break;

        default: return -1;
      }
      break;
  }

  if(p == -1 || nevadd < 0 || estart > nev){
    stringstream ss; emsg("UU");
  }

  if(nevadd != 0 && locpathst[p].size() != nevadd-1) emsg(" addpro");
  evnew.clear();
  e = estart;
  ji = ii;
  for(k = 0; k <= nevadd; k++){
    if(k < nevadd) tt = addtime[k]; else tt = tend;

    while(e < nev && indev[i][e].t <= tt){
      tr = indev[i][e].tr;
      if(tra[tr].cl != cl && tra[tr].cl >= 0){
        cf = ci + tra[tr].dc; tr2 = compiftra[ci][cf]; if(tr2 < 0) emsg("P"); 
        if(fixfl == 1){ if(tra[tr2].capev != tra[tr].capev) return -1;}
        EV ev; ev.t = indev[i][e].t; ev.tr = tr2; evnew.push_back(ev);
        ci = cf;
      }
      e++;
    }
    if(k == nevadd) break;

    if(k < nevadd-1) jf = locpathst[p][k]; else jf = ff;

    if(ji == nval){ 
      if(estart != 0) emsg("h"); cf = tra[indev[i][0].tr].cf; cf = cf + (jf-compval[cf][cl])*classmult[cl]; 
      if(tt == 0) tr = trabeg+cf; 
      else{ tr = compiftra[NOTALIVE][cf]; if(i < nind && indlikeenter[i] == 0) tr += moventra;}
    }
    else{
      if(jf == nval){ 
        if(tt == tmax) tr = traend+ci;
        else{ tr = compiftra[ci][NOTALIVE]; if(i < nind && indlikeleave[i] == 0) tr += moventra;}
      }
      else{ cf = ci + (jf-ji)*classmult[cl]; tr = compiftra[ci][cf];}
    }
    if(tr < 0){ if(tr != -1) emsg("err3"); return -1;}
    if(fixfl == 1){ if(tra[tr].capev >= 0 && capevall[tra[tr].capev] == 1 && tt > t && t < tend) return -1;}

    EV ev; ev.t = tt; ev.tr = tr; evnew.push_back(ev);
    ci = cf; ji = jf;
  }

  nevnew = evnew.size();       // Checks events are not too close
  if(nevnew > 0){
    for(e = 0; e < nevnew-1; e++){ if(evnew[e+1].t-evnew[e].t < evdtmin) return -1;}
    if(estart > 0){ if(evnew[0].t-indev[i][estart-1].t < evdtmin) return -1;}
    if(eend < nev-1){ if(indev[i][eend+1].t - evnew[nevnew-1].t < evdtmin) return -1;}
  }

  if(checkevdt() == 1) return -1;

  Ltoti = L();
  indchange(i,estart,t,tend);
  Ltotf = L();

  if(evfac == 1) fac *= double(nev-1)/(nindev[i]-1);

  al = exp(Ltotf-Ltoti)*fac; if(plfl == 1) cout << estart << " " <<  al << " " << fac << " " << Ltotf << "  "<< Ltoti<<  "al\n\n\n";
  //if(type == 2 && cl == 1) cout << al << "al\n";

  if(ran() < al){
    for(e = 0; e < nindev[i]-1; e++){
      if(indev[i][e+1].t-indev[i][e].t < evdtmin){ evnew = evrev; indchange(i,estart,t,tend); return 0;}
    }
    return 1;
  }
  else{ evnew = evrev; indchange(i,estart,t,tend); return 0;}
}

void localpropinit()
{
  long cl, nval, j, ii, jj, iii, jjj, iiii, ci, ff, k, kk;
  long tr, p;
  double sum;
  vector<long> vec;
  vector< vector<long> > transleave;

  ntrans = 0;
  for(tr = 0; tr < ntra; tr++){
    cl = tra[tr].cl;
    j = 0; while(j < ntrans && (trans[j].cl != cl || trans[j].i != tra[tr].i || trans[j].f != tra[tr].f)) j++;
    tra[tr].trans = j;
    if(j == ntrans){ TRANS tran; tran.cl = cl; tran.i = tra[tr].i; tran.f = tra[tr].f; trans.push_back(tran); ntrans++;}
  }

  locpath.resize(nclass-3); nlocpath.resize(nclass-3);
  for(cl = 0; cl < nclass-3; cl++){
    nval = nclassval[cl];
    locpath[cl].resize(nval+1); nlocpath[cl].resize(nval+1);
    for(j = 0; j <= nval; j++){
      locpath[cl][j].resize(nval+1); nlocpath[cl][j].resize(nval+1);
      for(jj = 0; jj <= nval; jj++){
        locpath[cl][j][jj].resize(4); nlocpath[cl][j][jj].resize(4);
        for(k = 0; k < 4; k++) nlocpath[cl][j][jj][k] = 0;
      }
    }

    transleave.clear();
    transleave.resize(nval+1);
    for(j = 0; j < ntrans; j++){ 
      if(trans[j].cl == cl){
        for(tr = 0; tr < ntra; tr++){   // Removes those trans which are all fixed
          if(tra[tr].cl == cl && tra[tr].i == trans[j].i && tra[tr].f == trans[j].f){
            if(tra[tr].capev < 0) break;
            else{
              if(capevall[tra[tr].capev] == 0) break;
            }
          }
        }
        if(tr < ntra) transleave[trans[j].i].push_back(trans[j].f);
      }
    }

    for(ii = 0; ii < nval; ii++){ transleave[nval].push_back(ii); transleave[ii].push_back(nval);}

    for(ii = 0; ii <= nval; ii++){     // no events 
      ff = ii;
      vec.clear();
      locpath[cl][ii][ff][0].push_back(locpathst.size()); nlocpath[cl][ii][ff][0]++;
      locpathst.push_back(vec); nlocpathst.push_back(vec.size());
    }

    for(ii = 0; ii <= nval; ii++){    // 1 event
      for(j = 0; j < transleave[ii].size(); j++){
        ff = transleave[ii][j];
        vec.clear();
        locpath[cl][ii][ff][1].push_back(locpathst.size()); nlocpath[cl][ii][ff][1]++;
        locpathst.push_back(vec); nlocpathst.push_back(vec.size());
      }
    }

    for(ii = 0; ii <= nval; ii++){    // 2 events
      for(j = 0; j < transleave[ii].size(); j++){
        iii = transleave[ii][j];
        if(iii != nval){
          for(jj = 0; jj < transleave[iii].size(); jj++){
            ff = transleave[iii][jj];
            vec.clear(); vec.push_back(iii);
            locpath[cl][ii][ff][2].push_back(locpathst.size()); nlocpath[cl][ii][ff][2]++;
            locpathst.push_back(vec); nlocpathst.push_back(vec.size());
          }
        }
      }
    }

    for(ii = 0; ii <= nval; ii++){    // 3 events
      for(j = 0; j < transleave[ii].size(); j++){
        iii = transleave[ii][j];
        if(iii != nval){
          for(jj = 0; jj < transleave[iii].size(); jj++){
            iiii = transleave[iii][jj];
            if(iiii != nval){
              for(jjj = 0; jjj < transleave[iiii].size(); jjj++){
                ff = transleave[iiii][jjj];
                vec.clear(); vec.push_back(iii); vec.push_back(iiii);
                locpath[cl][ii][ff][3].push_back(locpathst.size()); nlocpath[cl][ii][ff][3]++;
                locpathst.push_back(vec); nlocpathst.push_back(vec.size());
              }
            }
          }
        }
      }
    }

    for(ii = 0; ii <= nval; ii++){
      if(nlocpath[cl][ii][ii][2] > 0 && ii != nval) pairfl = 1;
      for(ff = 0; ff <= nval; ff++){
        if(nlocpath[cl][ii][ff][1] > 0 && nlocpath[cl][ii][ff][2] > 0) singeventfl = 1;
        if(nlocpath[cl][ii][ff][2] > 0 && nlocpath[cl][ii][ff][3] > 0) twothreefl = 1;
        j = 0; if(nlocpath[cl][ii][ff][0] > 0) j++; if(nlocpath[cl][ii][ff][1] > 0) j++; if(nlocpath[cl][ii][ff][2] > 0) j++;
        if(j >= 2 && nobs > 0) gapfl = 1;
      }
    }

    if(1 == 0){
      for(ii = 0; ii <= nval; ii++){
        for(iii = 0; iii <= nval; iii++){
          cout << ii << " " << iii << ":\n";
          for(kk = 0; kk < 4; kk++){
            cout << kk << " events: ";
            for(j = 0; j < nlocpath[cl][ii][iii][kk]; j++){
              p = locpath[cl][ii][iii][kk][j];
              for(k = 0; k < nlocpathst[p]; k++) cout << locpathst[p][k] << ",";
              cout << "*     ";
            }
            cout << "\n";
          }
        }
      }
    }
  }
}

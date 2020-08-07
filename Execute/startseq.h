// looks to make a sensible choice for the starting configuration
vector< vector< vector< vector<long> > > >  path;
vector< vector<double> >  cpos;
vector< vector<long> > cposto, cposfr;

struct PT { long dc; double t;};
bool comparePT(PT lhs, PT rhs) { return lhs.t < rhs.t; }
vector <PT> pt;

double Chain::startseq(long i)
{
  long e, c, cc, ob, j, p, cf, cs, csf, fev, cl, ti, tf, si, k, ctopp;
  double val, val2, sum, t, z;
  long tr, timeback, timefor;

  //if(indid[i] == "Ind. 4") cout << tent <<" " << tlea << " " << indtbirth[i] << " minmax:" << uocapmin << " " << uocapmax << "start\n";

  vector< vector <double> > prob;
  vector <long> pos;
  vector <double> possum;

  partcl = clall;
  indsim_init(i,indtbirth[i],1);

  timeprop[INITBA] -= clock();
  prob.resize(nevp);
  for(e = nevp-1; e >= 0; e--){
    prob[e].resize(ncomps);
    if(e == nevp-1){ for(c = 0; c < ncomps; c++) prob[e][c] = 1;}
    else{
      ctopp = getctop(evp[e].t,indtbirth[i]);

      for(c = 0; c < ncomps; c++) prob[e][c] = 0;

      for(cc = 0; cc < ncomps; cc++){
        val = prob[e+1][cc];
        if(val > 0){
          for(j = 0; j < cposto[cc].size(); j++){
            c = cposto[cc][j];
            val2 = val*cpos[c][cc];
            //if(e == 2) cout << val2 << " " << cc << " " << c << " " <<  checkforfixedinpath(c,cc,ctopp) << " hh\n";
            if(prob[e][c] < val2){ if(fixfl== 0 || checkforfixedinpath(c,cc,ctopp) == 0) prob[e][c] = val2;}
          }
        }
      }

      switch(evp[e].ty){
        case OBS:
          ob = evp[e].ch;
          for(c = 0; c < ncomps; c++){ if(obsprobeqn[ob][c+ctopp] < 0) prob[e][c] = 0;}
          break;

        case FIXEV:
          ti = fixev[evp[e].ch].trai; tf = fixev[evp[e].ch].traf; cl = fixev[evp[e].ch].cl;
          for(c = 0; c < ncomps; c++){ if(compval[c][cl] == ti) prob[e][c] = prob[e][c+(tf-ti)*classmult[cl]];}
          for(c = 0; c < ncomps; c++){ if(compval[c][cl] != ti) prob[e][c] = 0;}
          break;
      }
    }
  }
  timeprop[INITBA] += clock();

  //if(i == 31){
  /*
  if(indid[i] == "Ind. 4"){
    cout << indid[i] << " " << i <<  " id\n";
    for(e = 0; e < nevp; e++) cout << e << " " << evp[e].t << " "<<  evp[e].ty << " " << evp[e].ch << " evp\n";
    for(e = 0; e < nevp-1; e++){
      cout << "evp " << e << ":\n"; 
      for(c = 0; c < ncomps; c++){ if(prob[e][c] != 0) cout << compname[c] << ":" << prob[e][c] << ", ";} cout<< "\n";
    }
    emsg("P");
  }
  */
  evnew.clear();

  sum = 0; pos.clear(); possum.clear();
  t = evp[0].t;
  for(c = 0; c < ncomps; c++){
    val = prob[0][c];
    if(val > 0){
      if(t == 0 || compiftra[NOTALIVE][c+ctop] >= 0){
        if(t == 0) val *= dirchal[c%ncompswa]; // takes account of Diriclet prior on initial state
        sum += val;
        pos.push_back(c); possum.push_back(sum);
      }
    }
  }

  if(sum == 0) return 0;
  //if(indid[i] == "Ind. 4") cout << "Here\n";

  z = ran()*sum; p = 0; while(p < pos.size() && z > possum[p]) p++;
  c = pos[p]+ctop;

  if(t == 0) tr = trabeg+c;
  else{ 
    tr = compiftra[NOTALIVE][c]; if(tr < 0) return 0;
    if(i < nind && indlikeenter[i] == 0) tr += moventra;
    if(indfixtenter[i] != large && indfixtenter[i] != t) return 0;
    if(tra[tr].capev >= 0){
      if(capevall[tra[tr].capev] == 1 && indfixtenter[i] != t) return 0;
    }
  }
  EV evbeg; evbeg.t = t; evbeg.tr = tr; evnew.push_back(evbeg);

  for(e = 0; e < nevp-1; e++){
    t = evp[e].t; ctopp = getctop(t,indtbirth[i]);
    switch(evp[e].ty){
      case CHC:
        cf = c + evp[e].ch;
        tr = compiftra[c][cf]; if(tr < 0) return 0;
        EV ev; ev.tr = tr; ev.t = t; evnew.push_back(ev);
        c = cf;
        break;

      case FIXEV:
        fev = evp[e].ch; cl = fixev[fev].cl;
        if(compval[c][cl] != fixev[fev].trai) return 0;
        cf = c + fixev[fev].dc;
        tr = compiftra[c][cf]; if(tr < 0) return 0;
        if(fixev[fev].like == 0) tr += moventra;
        EV evf; evf.tr = tr; evf.t = t; evnew.push_back(evf);
        c = cf;
        break;
    }

    sum = 0; pos.clear(); possum.clear();
    cs = c%ncomps;// ctop = c-cs;

    if(evp[e+1].t == evp[e].t) csf = cs;   // if there is no time then must be the same (e.g. observation at same time as change as age/time/fixed time lines
    else{
      for(j = 0; j < cposfr[cs].size(); j++){
        csf = cposfr[cs][j];
        val = prob[e+1][csf];
        if(val > 0){
          if(fixfl == 0 || checkforfixedinpath(cs,csf,ctopp) == 0){ sum += val*cpos[cs][csf]; pos.push_back(csf); possum.push_back(sum);}
        }
      }
      if(sum == 0) return 0;
      z = ran()*sum; p = 0; while(p < pos.size() && z > possum[p]) p++;
      csf = pos[p];
    }

    pt.clear();
    for(cl = 0; cl < nclass-3; cl++){
      ti = compval[cs][cl]; tf = compval[csf][cl];
      if(ti != tf){
        if(evp[e+1].t == evp[e].t){ cout << e << " " << nevp << " " << i << "\n"; emsg("same time");}
        si = path[cl][ti][tf].size(); if(si == 0) emsg("tt");
        for(k = 0; k < si; k++){
          tr = path[cl][ti][tf][k];
          PT pta; pta.t =  evp[e].t + (k+0.1+0.9*ran())*(evp[e+1].t -evp[e].t)/si; pta.dc = (trans[tr].f-trans[tr].i)*classmult[cl];
          pt.push_back(pta);
        }
      }
    }
    sort(pt.begin(),pt.end(),comparePT);

    for(k = 0; k < pt.size(); k++){
      cf = c + pt[k].dc;
      tr = compiftra[c][cf]; if(tr < 0) return 0;
      EV evf; evf.tr = tr; evf.t = pt[k].t; evnew.push_back(evf);
      c = cf;
    }
  }
  t = evp[nevp-1].t;

  if(t == tmax) tr = traend+c; 
  else{
    tr = compiftra[c][NOTALIVE]; if(tr < 0) return 0;
    if(i < nind && indlikeleave[i] == 0) tr += moventra;
    if(indfixtleave[i] != large && indfixtleave[i] != t) return 0;
    if(tra[tr].capev >= 0){
      if(capevall[tra[tr].capev] == 1 && indfixtleave[i] != t) return 0;
    }
  }
  EV evend; evend.t = t; evend.tr = tr; evnew.push_back(evend);

  //cout << tent << " " << tlea << "tlea\n";
  if(numwrong(i,evnew) > 0){ cout << indid[i] << "  id\n"; oe("r",evnew); emsg("num wrong");}
  //if(indid[i] == "Ind. 4") cout << probuocap(i,evnew) << " here2\n";

  return probuocap(i,evnew);
}

long getctop(double t, double tbirth)                    // Returns the compartment part related to the top three classifications
{
  long ctop, ag, ti, f;

  ctop = 0;
  if(tbirthfl == 1){
    ag = 0; while(ag < nage && tbirth + age[ag] < t) ag++;
    ctop += ag*classmult[agecl];
  }

  ti = 0; while(ti < nsettime && settime[ti] < t) ti++;
  ctop += ti*classmult[settimecl];

  f = 0; while(f < ncapevtrange && capevtrange[f] <= t) f++;
  ctop += f*classmult[capevcl];

  return ctop;
}

long checkforfixedinpath(long c, long cc, long ctopp)           // Considers if fixed events stop a path from occuring
{
  long cl, ti, tf, si, k, tr, tran, ci, cf;

  for(cl = 0; cl < nclass-3; cl++){
    ti = compval[c][cl];
    tf = compval[cc][cl];
    if(ti != tf){
       si = path[cl][ti][tf].size();
       ci = c+ctopp;
       for(k = 0; k < si; k++){
         tran = path[cl][ti][tf][k];
         cf = ci + (trans[tran].f-trans[tran].i)*classmult[cl];
         tr = compiftra[ci][cf];
         if(tr < 0) return 1;
         if(tra[tr].capev >= 0){ if(capevall[tra[tr].capev] == 1) return 1;}
         ci = cf;
       }
    }
  }
  return 0;
}

void startseqinit()        // Initialises quantitites used in seting up the initial sequence
{
  long cl, nval, j, jj, jjj, k, si, sinew, fl, tr, c, cc, num, i, f;
  vector< vector<long> > transleave;

  path.resize(nclass-3);
  for(cl = 0; cl < nclass-3; cl++){
    nval = nclassval[cl];

    path[cl].resize(nval);
    for(j = 0; j < nclassval[cl]; j++) path[cl][j].resize(nval);
 
    for(tr = 0; tr < ntrans; tr++){
      if(trans[tr].cl == cl){ path[cl][trans[tr].i][trans[tr].f].push_back(tr);}
    }

    transleave.clear();
    transleave.resize(nval);
    for(j = 0; j < ntrans; j++){
      if(trans[j].cl == cl){
        for(tr = 0; tr < ntra; tr++){   // Removes those trans which are all fixed
          if(tra[tr].cl == cl && tra[tr].i == trans[j].i && tra[tr].f == trans[j].f && tra[tr].like == 1){
            if(tra[tr].capev < 0) break;
            else{
              if(capevall[tra[tr].capev] != 1) break;
            }
          }
        }
        if(tr < ntra) transleave[trans[j].i].push_back(j);
      }
    }
    //for(c = 0; c < nclassval[cl]; c++){ for(j = 0; j < transleave[c].size(); j++) cout <<  transleave[c][j] << " tr\n";}

    do{
      fl = 0;
      for(j = 0; j < nval; j++){
        for(jj = 0; jj < nval; jj++){
          si = path[cl][j][jj].size();
          if(j != jj && si > 0){
            for(k = 0; k < transleave[jj].size(); k++){
              tr = transleave[jj][k];
              jjj = trans[tr].f;
              if(jjj != j){
                sinew = path[cl][j][jjj].size();
                if(sinew == 0 || si+1 < sinew){ path[cl][j][jjj] = path[cl][j][jj]; path[cl][j][jjj].push_back(tr); fl = 1;}
              }
            }
          }
        }
      }
    }while(fl == 1);
    /*
    for(j = 0; j < nval; j++){
      for(jj = 0; jj < nval; jj++){
        cout << j << " " << jj << "  "; for(k = 0; k < path[cl][j][jj].size(); k++) cout << trans[path[cl][j][jj][k]].i << "," << trans[path[cl][j][jj][k]].f << "   ,"; 
        cout << "\n";
      }
    }
    */
  }

  cpos.resize(ncomps);
  for(c = 0; c < ncomps; c++){
    cpos[c].resize(ncomps);
    for(cc = 0; cc < ncomps; cc++){
      num = 0;
      for(cl = 0; cl < nclass-3; cl++){
        i = compval[c][cl];
        f = compval[cc][cl];
        if(i != f){ si = path[cl][i][f].size(); if(si == 0) break; num += si;}
      }
      if(cl < nclass-3) cpos[c][cc] = 0;
      else{
        if(notimerange == 1) cpos[c][cc] = exp(-100*num); // If no data range then make transistion very unlikely
        else cpos[c][cc] = exp(-1*num);
      }
    }
  }

  cposfr.resize(ncomps); cposto.resize(ncomps);
  for(cc = 0; cc < ncomps; cc++){
    for(c = 0; c < ncomps; c++){
      if(cpos[c][cc] > 0){ cposto[cc].push_back(c); cposfr[c].push_back(cc);}
    }
  }
}

double Chain::probuocap(long i, vector<EV> &vec)                                       // Gets the probability of captures which are not observed
{
  long nev, ca, c, e, j, eq, ob;
  double tbeg, tend, t, prob;

  nev = vec.size();
  tbeg = vec[0].t; tend = vec[nev-1].t;

  ca = 0; while(ca < ncap && capt[ca] < tbeg) ca++;

  prob = 1;
  c = NOTALIVE; e = 0;
  while(ca < ncap && capt[ca] < tend){
    t = capt[ca];
    while(e < nev && vec[e].t < t){ c = tra[vec[e].tr].cf; e++;}

    if(c != NOTALIVE){
      eq = capprobeqn[ca][c];
      if(eq != -1){
        for(j = 0; j < nindobs[i]; j++){ if(t == obst[indobs[i][j]]) break;}
        if(j == nindobs[i]){
          if(eqnstr[eq] == "1"){
            if(t < firstobst[i]){ if(t > uocapmin) uocapmin = t;}
            if(t > lastobst[i]){ if(t < uocapmax) uocapmax = t;}
            //cout <<  indid[i] << " " << uocapmin << " " << uocapmax << "  uoran\n";;

            stringstream ss; ss << "Individual '" << indid[i] << "' should be observed in capture '" << capname[ca] << "' but it is not.";
            if(capwarn != "not"){
              if(capwarn == "") capwarn = ss.str();
              else{
                if(capwarn != ss.str()){
                  stringstream ss2; ss2 << "Individual '" << indid[i] << "' should be observed in other captures but it is not.";
                  capwarn = ss.str();
                }
              }
            }
            return 0;
          }
          //if(prob > 0.2) 
          prob *= capprobest;
        }
      }
    }
    ca++;
  }

  c = NOTALIVE; e = 0;
  for(j = 0; j < nindobs[i]; j++){                       // Takes into accout observations
    ob = indobs[i][j];
    t = obst[ob];
    while(e < nev && vec[e].t < t){ c = tra[vec[e].tr].cf; e++;}

    eq = obsprobeqn[ob][c];
    if(eq < 0){ prob = 0; break;}
    prob *= calcobsprob(eq);
  }

  for(e = 0; e < nev-1; e++){ if(vec[e+1].t - vec[e].t < evdtmin) prob = 0;} // makes sure not too close

  capwarn = "not";

  //prob *= exp(-(tend-tbeg)*muest); 
  //if(loopstart > loopstartmax/2 && prob < 0.05) prob = 0.05;
  //if(loopstart > 0.9*loopstartmax) prob = 1;

  return prob*exp(-(tend-tbeg)*muest);
}

short checkiffixed(long i, double t)
{
  long e;

  if(i >= nind) return 0;
  e = 0; while(e < nindfixev[i] && fixev[indfixev[i][e]].t < t) e++;
  if(e < nindfixev[i] && fixev[indfixev[i][e]].t == t) return 1;
  if(t == indfixtenter[i] || t == indfixtleave[i]) return 1;
  return 0;
}

void Chain::move_prop(long i)                                        // moves the events in individual i
{
  long loop, e, nev, cl, fl, ee, dc;
  long tr, tr2, trr;
  double t, tst, Ltoti, Ltotf, al;

  getlifespan(i);
  nev = nindev[i];
  for(loop = 0; loop < nev-2; loop++){
    e = 1+long(ran()*(nev-2));
    tr = indev[i][e].tr;
    tst = indev[i][e].t;
    if(checkiffixed(i,tst) == 0){  // Makes sure the event is not fixed
      fl = 0;
      switch(tra[tr].type){                                             // checks for a valid move (not crossing with another in the same classification)
        case EXP_TR: case GAMMA_TR: case WEI_TR:
          cl = tra[tr].cl; if(cl < 0) emsg("pr1");

          t = normal(tst,jump_move[i][tra[tr].trans]);
          dc = tra[tr].dc;
          if(t > tlea || t < tent) fl = 1;
          else{
            if(t > tst){
              ee = e+1;
              while(ee < nev-1 && indev[i][ee].t < t){
                tr2 = indev[i][ee].tr; if(tra[tr2].cl == cl){ fl = 1; break;}
                trr = compiftra[tra[tr2].ci-dc][tra[tr2].cf-dc]; if(trr < 0){ fl = 1; break;}
                if(fixfl == 1){ if(tra[trr].capev != tra[tr2].capev){ fl = 1; break;}}
                ee++;
              }
              if(indev[i][ee].t-t < evdtmin || t-indev[i][ee-1].t < evdtmin) fl = 1;
            }
            else{
              ee = e-1;
              while(ee > 0 && indev[i][ee].t > t){
                tr2 = indev[i][ee].tr; if(tra[tr2].cl == cl){ fl = 1; break;}
                trr = compiftra[tra[tr2].ci+dc][tra[tr2].cf+dc]; if(trr < 0){ fl = 1; break;}
                if(fixfl == 1){ if(tra[trr].capev != tra[tr2].capev){ fl = 1; break;}}
                ee--;
              }
              if(t-indev[i][ee].t < evdtmin || indev[i][ee+1].t-t < evdtmin) fl = 1;
            }
            if(fixfl == 1 && fl == 0){
              while(ee > 0 && (indev[i][ee].t > t || ee == e)) ee--;
              tr2 = indev[i][ee].tr; trr = compiftra[tra[tr2].cf][tra[tr2].cf+dc]; 
              if(trr < 0 || (tra[trr].capev >= 0 && capevall[tra[trr].capev] == 1)) fl = 1;
            }
          }
          break;

        default: fl = 2; break;
      }

      if(fl != 2){
        if(fl == 1) al = -1;
        else{
          Ltoti = L(); //oe("T",indev[i]); cout << e  << " " << t << "  move\n";
          e = eventmove(i,e,t);
          Ltotf = L(); //oe("A",indev[i]);

          al = exp(Ltotf - Ltoti); //if(samp > 1000 && e == 3) cout << tst << " "  << t << " " << jump_move[i][tra[tr].trans] <<  " " << al  << "al\n";
        }

        ntr_move[i][cl]++;
        if(ran() < al){
          nac_move[i][cl]++;

          if(samp < burnin){ jump_move[i][tra[tr].trans] *= 1.05;}
        }
        else{
          if(al >= 0) eventmove(i,e,tst);

          if(samp < burnin){ jump_move[i][tra[tr].trans] *= 0.95;}
        }
      }
    }
  }
  if(nindev[i] != nev) emsg("pr");
}

long Chain::eventmove(long i, long e, double tnew)        // Does fast movement of events
{
  long dc, ob, nob, nev, c, o, coi, cof, n;
  double t, tev, tobs;
  long tr, tr2;

  t = indev[i][e].t; tr = indev[i][e].tr; dc = tra[tr].dc;

  if(capevfl == 1){
    long ce = tra[tr].capev;
    if(ce >= 0){ capev_oneminus[ce]--; Lev_pd -= capev_valoneminus[ce];}
  }

  if(nmfl == 1){                  // works out any changes to non-Markovian likelihood
    long ee, cl;
    long trr;
    double tbeg, tend;

    if(tra[tr].nm == 1 || tra[tr].ntraend > 0){
      cl = tra[tr].cl;
      ee = e-1; while(ee > 0 && tra[indev[i][ee].tr].cl != cl) ee--;
      tbeg = indev[i][ee].t;
      nmchange(tr,t-tbeg,tnew-tbeg); for(n = 0; n < tra[tr].ntraend; n++) nmupchange(tra[tr].traend[n],t-tbeg,tnew-tbeg);
    }

    if(tra[tr].nnonexptra > 0){
      cl = tra[tr].cl;
      ee = e+1; while(ee < nindev[i]-1 && tra[indev[i][ee].tr].cl != cl) ee++;
      tend = indev[i][ee].t; trr = indev[i][ee].tr;
      if(tra[trr].nm == 1) nmchange(trr,tend-t,tend-tnew);
      for(n = 0; n < tra[trr].ntraend; n++){ tr2 = tra[trr].traend[n]; if(tra[tr2].cl == cl) nmupchange(tr2,tend-t,tend-tnew);}
    }
  }

  indremevent(i,e);

  o = 0; nob = nindobs[i]; while(o < nob && indobst[i][o] < t) o++;
  nev = nindev[i];
  if(tnew > t){
    c = tra[tr].cf;
    if(o < nob) tobs = indobst[i][o]; else tobs = large;

    if(e < nev) tev = indev[i][e].t; else tev = large;
    while(tobs < tnew || tev < tnew){
      if(tev < tobs){
        secchange(i,c,c-dc,t,tev); t = tev;
        EV evcha; evcha = indev[i][e]; tr2 = evcha.tr;
        evcha.tr = compiftra[tra[tr2].ci-dc][tra[tr2].cf-dc];
        if(evcha.tr == -1) emsg("not there1");
        if(tra[tr2].like == 0) evcha.tr += moventra;
        indremevent(i,e); indaddevent(i,e,evcha);
        c = tra[tr2].cf;
        e++; if(e < nev) tev = indev[i][e].t; else tev = large;
      }
      else{
        ob = indobs[i][o]; coi = obsprob_cope[ob][c]; cof = obsprob_cope[ob][c-dc];
        if(coi != cof){
          if(coi < 0) Lob_st -= notobsdL; else{ Lob_st -= cope_val[coi]; cope_num[coi]--;}
          if(cof < 0) Lob_st += notobsdL; else{ Lob_st += cope_val[cof]; cope_num[cof]++;}
        }
        o++; if(o < nob) tobs = indobst[i][o]; else tobs = large;
      }
    }
    secchange(i,c,c-dc,t,tnew);
    tr = compiftra[c-dc][c];
    EV evnew; evnew.t = tnew; evnew.tr = tr;
    indaddevent(i,e,evnew);
  }
  else{
    c = tra[tr].ci; o--; e--;
    if(o >= 0) tobs = indobst[i][o]; else tobs = -large;

    if(e >= 0) tev = indev[i][e].t; else tev = -large;

    while(tobs > tnew || tev > tnew){
      if(tev > tobs){
        secchange(i,c,c+dc,tev,t); t = tev;
        EV evcha; evcha = indev[i][e]; tr2 = evcha.tr; evcha.tr = compiftra[tra[tr2].ci+dc][tra[tr2].cf+dc]; 
        if(evcha.tr == -1) emsg("not there2");
        if(tra[tr2].like == 0) evcha.tr += moventra;
        indremevent(i,e); indaddevent(i,e,evcha);
        c = tra[tr2].ci;
        e--; if(e >= 0) tev = indev[i][e].t; else tev = -large;
      }
      else{
        ob = indobs[i][o]; coi = obsprob_cope[ob][c]; cof = obsprob_cope[ob][c+dc];
        if(coi != cof){
          if(coi < 0) Lob_st -= notobsdL; else{ Lob_st -= cope_val[coi]; cope_num[coi]--;}
          if(cof < 0) Lob_st += notobsdL; else{ Lob_st += cope_val[cof]; cope_num[cof]++;}
        }
        o--; if(o >= 0) tobs = indobst[i][o]; else tobs = -large;
      }
    }
    secchange(i,c,c+dc,tnew,t);
    tr = compiftra[c][c+dc];
    EV evnew; evnew.t = tnew; evnew.tr =tr; e++;
    indaddevent(i,e,evnew);
  }

  if(capevfl == 1){
    long ce = tra[tr].capev;
    if(ce >= 0){ capev_oneminus[ce]++; Lev_pd += capev_valoneminus[ce];}
  }

  return e;
}

void Chain::nmupchange(long tr, double dt, double dt2)        // Calculates the chnage in Linm from moving an evvent
{
  double mean, kshape, lam;

  switch(tra[tr].type){
    case GAMMA_TR:
      mean = nmeq_val[tra[tr].eq]; kshape = nmeq_val[tra[tr].eqshape];
      Linm += gammaup(dt2,kshape,kshape/mean) - gammaup(dt,kshape,kshape/mean);
      break;

    case WEI_TR:
      lam = nmeq_val[tra[tr].eq]; kshape = nmeq_val[tra[tr].eqshape];
      Linm += weibullup(dt2,lam,kshape) - weibullup(dt,lam,kshape);
      break;
      case FIXED_TR: emsg("cannot nove"); break;
  }
}

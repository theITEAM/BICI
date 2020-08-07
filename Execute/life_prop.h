void Chain::life_prop(long i)                    // Makes changes to entry and leaving times as well as time of birth
{
  long type, fl, elea, nev, e;
  double p, tbirth, tnew, facup = 1.04, facdown = 0.98, facup2 = 1.004, facdown2 = 0.998, dt, tlimmin, tlimmax, fot, lot;

  p = ran(); type = 0; while(type <  7 && p > problifepropsum[type]) type++; if(type == 7) emsg("err");

  if(sourcefl == 0 && (type == 0 ||type == 2 || type == 4 || type == 5)) return;
  if(sinkfl == 0 && (type == 1 || type == 6)) return;
  if(tbirthfl == 0 && (type == 2 || type == 3 || type == 4)) return;

  elea = getlifespan(i);
  switch(type){
    case 0:    // change birth and enter time
      if(tent == 0) return;
      if(i < nind && indfixtenter[i] != large) return;

      if(i < nind) ntr_tbirthent[i]++; else ntr_tbirthentuo++; 
      if(i < nind) dt = normal(0,jump_tent[i]); else dt = normal(0,jump_tentuo);
      tnew = tent + dt;

      nev = nindev[i];
      e = 1; while(e < nev && tra[indev[i][e].tr].cl >= nclass-3) e++; if(e < nev) tlimmax = indev[i][e].t; else tlimmax = large;

      if(i < nind) fot = firstobst[i]; else fot = large;
      if(tnew < 0 || tnew > tmax || (tnew > tent && tnew > tlimmax) || tnew > fot) fl = 0;
      else fl = lifechange(i,tnew,tlea,indtbirth[i]+dt,0);

      if(fl == 1){
        if(i < nind) nac_tbirthent[i]++; else nac_tbirthentuo++;
        if(samp < burnin){ if(i < nind) jump_tbirthent[i] *= facup; else jump_tbirthentuo *= facup2;}
      }
      else{
        if(samp < burnin){ if(i < nind) jump_tbirthent[i] *= facdown; else jump_tbirthentuo *= facdown2;}
      }
      break;

    case 1:    // change leave time
      if(tlea == tmax) return;
      if(i < nind && indfixtleave[i] != large) return;

      if(i < nind) ntr_tlea[i]++; else ntr_tleauo++;
      if(i < nind) tnew = tlea + normal(0,jump_tlea[i]); else tnew = tlea + normal(0,jump_tleauo);

      nev = nindev[i];
      e = elea+1; while(e < nev && tra[indev[i][e].tr].cl >= nclass-3) e++; if(e < nev) tlimmax = indev[i][e].t; else tlimmax = large;
      e = elea-1; while(e >= 0 && tra[indev[i][e].tr].cl >= nclass-3) e--; if(e >= 0) tlimmin = indev[i][e].t; else tlimmin = -large;

      if(i < nind) lot = lastobst[i]; else lot = -large;
      if(tnew < 0 || tnew > tmax || tnew < lot || tnew < tlimmin || tnew > tlimmax) fl = 0;
      else fl = lifechange(i,tent,tnew,indtbirth[i],0);

      if(fl == 1){
        if(i < nind) nac_tlea[i]++; else nac_tleauo++;
        if(samp < burnin){ if(i < nind) jump_tlea[i] *= facup; else jump_tleauo *= facup2;}
      }
      else{
        if(samp < burnin){ if(i < nind) jump_tlea[i] *= facdown; else jump_tleauo *= facdown2;}
      }
      break;

    case 2:    // change enter time
      if(tent == 0) return;
      if(i < nind && indfixtenter[i] != large) return;

      if(i < nind) ntr_tent[i]++; else ntr_tentuo++;
      if(i < nind) tnew = tent + normal(0,jump_tent[i]); else tnew = tent + normal(0,jump_tentuo);

      nev = nindev[i];
      e = 1; while(e < nev && tra[indev[i][e].tr].cl >= nclass-3) e++; if(e < nev) tlimmax = indev[i][e].t; else tlimmax = large;

      if(i < nind) fot = firstobst[i]; else fot = large;
      if(tnew < 0 || tnew > tmax || (tnew > tent && tnew > tlimmax) || tnew < indtbirth[i] || tnew > fot || tnew > indtbirth[i]+age[nage]) fl = 0;
      else fl = lifechange(i,tnew,tlea,indtbirth[i],0);

      if(fl == 1){
        if(i < nind) nac_tent[i]++; else nac_tentuo++;
        if(samp < burnin){ if(i < nind) jump_tent[i] *= facup; else jump_tentuo *= facup2;}
      }
      else{
        if(samp < burnin){ if(i < nind) jump_tent[i] *= facdown; else jump_tentuo *= facdown2;}
      }
      break;

    case 3:    // change birth time
      if(i < nind) ntr_tbirth[i]++; else ntr_tbirthuo++;
      if(i < nind) tbirth = indtbirth[i] + normal(0,jump_tbirth[i]); else tbirth = indtbirth[i] + normal(0,jump_tbirthuo);

      if(tbirth >= tent || tbirth < tent - age[nage]) fl = 0;
      else fl = lifechange(i,tent,tlea,tbirth,0);

      if(fl == 1){
        if(i < nind) nac_tbirth[i]++; else nac_tbirthuo++;
        if(samp < burnin){ if(i < nind) jump_tbirth[i] *= facup; else jump_tbirthuo *= facup2;}
      }
      else{
        if(samp < burnin){ if(i < nind) jump_tbirth[i] *= facdown; else jump_tbirthuo *= facdown2;}
      }
      break;

    case 4:    // switch enter
      if(i < nind && indfixtenter[i] != large) return;
      if(indtbirth[i] <= 0){
       nev = nindev[i]; e = 1; while(e < nev && tra[indev[i][e].tr].cl >= nclass-3) e++; if(e < nev) tlimmax = indev[i][e].t; else tlimmax = large;

       if(i < nind) dt = firstobst[i]; else dt = large;
       if(dt > tmax/20) dt = tmax/20; if(tlimmax < dt) dt = tlimmax;
       if(dt > indtbirth[i] + age[nage]) dt = indtbirth[i] + age[nage];

        if(tent == 0){
          tnew = ran()*dt;
          if(i < nind) ntr_entswitch[i]++; else ntr_entswitchuo++;
          if(lifechange(i,tnew,tlea,indtbirth[i],log(dt)) == 1){
            if(i < nind) nac_entswitch[i]++; else nac_entswitchuo++;
          }
        }
        else{
          if(tent > dt) return;
          if(i < nind) ntr_entswitch[i]++; else ntr_entswitchuo++;
          if(lifechange(i,0,tlea,indtbirth[i],-log(dt)) == 1){
            if(i < nind) nac_entswitch[i]++; else nac_entswitchuo++;
          }
        }
      }
      break;

    case 5:    // switch birth enter
      if(i < nind && indfixtenter[i] != large) return;
      nev = nindev[i]; e = 1; while(e < nev && tra[indev[i][e].tr].cl >= nclass-3) e++; if(e < nev) tlimmax = indev[i][e].t; else tlimmax = large;
      if(i < nind) dt = firstobst[i]; else dt = large;
      if(dt > tmax/20) dt = tmax/20; if(tlimmax < dt) dt = tlimmax;
      if(tent == 0){
        tnew = ran()*dt;
        if(i < nind) ntr_birthentswitch[i]++; else ntr_birthentswitchuo++;
        if(lifechange(i,tnew,tlea,indtbirth[i]+tnew,log(dt)) == 1){
          if(i < nind) nac_birthentswitch[i]++; else nac_birthentswitchuo++;
        }
      }
      else{
        if(tent > dt) return;
        if(i < nind) ntr_birthentswitch[i]++; else ntr_birthentswitchuo++;
        if(lifechange(i,0,tlea,indtbirth[i]-tent,-log(dt)) == 1){
          if(i < nind) nac_birthentswitch[i]++; else nac_birthentswitchuo++;
        }
      }
      break;

    case 6:    // switch leave
      if(i < nind && indfixtleave[i] != large) return;

      nev = nindev[i];
      e = elea+1; while(e < nev && tra[indev[i][e].tr].cl >= nclass-3) e++; if(e < nev) tlimmax = indev[i][e].t; else tlimmax = large;
      e = elea-1; while(e >= 0 && tra[indev[i][e].tr].cl >= nclass-3) e--; if(e >= 0) tlimmin = indev[i][e].t; else tlimmin = -large;

      if(i < nind) dt = tmax-lastobst[i]; else dt = large;
      if(dt > tmax/20) dt = tmax/20; if(tlimmax < tmax) return; 
      if(tmax-tlimmin < dt) dt = tmax-tlimmin;

      if(tlea == tmax){
        tnew = tmax-ran()*dt;
        if(i < nind) ntr_leaswitch[i]++; else ntr_leaswitchuo++;
        if(lifechange(i,tent,tnew,indtbirth[i],log(dt)) == 1){
          if(i < nind) nac_leaswitch[i]++; else nac_leaswitchuo++;
        }
      }
      else{
        if(tlea < tmax-dt) return;
        if(i < nind) ntr_leaswitch[i]++; else ntr_leaswitchuo++;
        if(lifechange(i,tent,tmax,indtbirth[i],-log(dt)) == 1){
          if(i < nind) nac_leaswitch[i]++; else nac_leaswitchuo++;
        }
      }
      break;
  }
}

long Chain::lifechange(long i, double tentnew, double tleanew, double tbirthnew, double dprob)    // Makes a change to the life
{
  long e, ag, agold, ti, nev, c, cf, f;
  double al, Ltoti, Ltotf, dLpri;
  long tr, trr;

  //evrev = indev[i];

  ag = 0; while(ag < nage && tbirthnew + age[ag] < tentnew) ag++;
  ctop = ag*classmult[agecl];

  ti = 0; while(ti < nsettime && settime[ti] <= tentnew) ti++;
  ctop += ti*classmult[settimecl];

  f = 0; while(f < ncapevtrange && capevtrange[f] <= tentnew) f++;
  ctop += f*classmult[capevcl];

  tr = indev[i][0].tr;
  if(tbirthfl == 1){ agold = compval[tra[tr].cf][agecl]; dLpri = agprior[ag] - agprior[agold];}else dLpri = 0;

  c = tra[tr].cf%ncomps + ctop;

  evnew.clear();

  if(tentnew == 0) tr = trabeg+c; 
  else{
    tr = compiftra[NOTALIVE][c]; if(tr < 0) return 0;
    if(sourcefixfl == 1){
      if(i < nind && indlikeenter[i] == 0) tr += moventra;
      if(i < nind && indfixtenter[i] != large){
        if(indfixtenter[i] != tentnew) return 0;
      }
      if(tra[tr].capev >= 0 && capevall[tra[tr].capev] == 1){
        if(i >= nind) return 0; else{ if(indfixtenter[i] != tentnew) return 0;}
      }
    }
  }
  EV evbeg; evbeg.t = tentnew; evbeg.tr = tr; evnew.push_back(evbeg);

  for(e = 0; e < nindev[i]; e++){
    switch(tra[indev[i][e].tr].type){
      //case BEG_TR: evnew.push_back(evrev[e]); evnew[evnew.size()-1].t = tentnew; break;
      //case END_TR: evnew.push_back(evrev[e]); evnew[evnew.size()-1].t = tleanew; break;
      case BEG_TR: case END_TR: case  GROW_TR: case SETTIME_TR: break;
      default: if(tra[indev[i][e].tr].cl >= 0) evnew.push_back(indev[i][e]); break;
    }
  }

  while(ag < nage && tbirthnew + age[ag] < tleanew){
    EV ev; ev.t =  tbirthnew + age[ag]; ev.tr = compiftra[ag*classmult[agecl]][(ag+1)*classmult[agecl]]; evnew.push_back(ev);
    ag++;
  }

  while(ti < nsettime && settime[ti] < tleanew){
    EV ev; ev.t = settime[ti]; ev.tr = compiftra[ti*classmult[settimecl]][(ti+1)*classmult[settimecl]]; evnew.push_back(ev);
    ti++;
  }

  while(f < ncapevtrange && capevtrange[f] < tleanew){
    EV ev; ev.t = settime[ti]; ev.tr = compiftra[f*classmult[capevcl]][(f+1)*classmult[capevcl]]; evnew.push_back(ev);
    f++;
  }

  EV evend; evend.t = tleanew; evend.tr = tranull; evnew.push_back(evend);

  sort(evnew.begin(),evnew.end(),compareev);

  nev = evnew.size();
  for(e = 1; e < nev; e++){
    tr = evnew[e].tr;
    if(tra[tr].cf == NOTALIVE){
      if(evnew[e].t == tmax) trr = traend+c;
      else{ 
        trr = compiftra[c][NOTALIVE]; if(trr < 0) return 0;
        if(sinkfixfl == 1){
          if(i < nind && indlikeleave[i] == 0) tr += moventra;

          if(i < nind && indfixtleave[i] != large){
            if(indfixtleave[i] != evnew[e].t) return 0;
          }
          if(tra[trr].capev >= 0 && capevall[tra[trr].capev] == 1){
            if(i >= nind) return 0; else{ if(indfixtleave[i] != evnew[e].t) return 0;}
          }
        }
      }
      evnew[e].tr = trr;
      break;
    }
    else{
      if(tra[tr].type >= 0){ cf = c + tra[tr].dc; evnew[e].tr = compiftra[c][cf]; c = cf;}
    }
  }
  //outputevents("BEF",indev[i]); outputevents("AFT",evnew);

  if(e == nev) emsg("prob e");

  if(checkevdt() == 1) return 0;

  Ltoti = L(); //cout << "bef:"; Lout();
  indcha(i);
  Ltotf = L(); //cout << "aft:"; Lout();

  al = exp(Ltotf - Ltoti + dprob + dLpri); //cout << al << " " << dprob << " " << dLpri << "al\n";
  if(ran() < al){ indtbirth[i] = tbirthnew; Lpri += dLpri; return 1;}
  else{ indrev(i); return 0; }
}


void lifepropinit()     // Initialises the life proposals
{
  long j;
  double sum;

  if(sourcefl == 0){ problifeprop[0] = 0; problifeprop[2] = 0; problifeprop[4] = 0; problifeprop[5] = 0;}
  if(sinkfl == 0){ problifeprop[1] = 0; problifeprop[6] = 0;}
  if(tbirthfl == 0){ problifeprop[2] = 0; problifeprop[3] = 0; problifeprop[4] = 0;}

  sum = 0; for(j = 0; j < 7; j++) sum += problifeprop[j];
  for(j = 0; j < 7; j++) problifeprop[j] /=sum;
  sum = 0; for(j = 0; j < 7; j++){ sum += problifeprop[j]; problifepropsum[j] = sum;}
}


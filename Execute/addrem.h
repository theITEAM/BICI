vector <double> induoent;                                    // The distribution for adding individuals at  t > 0
vector <double> induoinit;                                    // The distribution for adding individuals at t = 0
double induolifetime;                                      // The expected lifetime for added individuals
double nuoinit;                                           // The number of indivviduals which start at t = 0
double nuoend;                                            // The number of indivviduals which end at t = tmax
double nuosource;                                         // The number of indivviduals which enter the populatiom
double nuosink;                                         // The number of indivviduals leave the population
 

void  Chain::addrem_prop()                                                   // Adds and remove indiduals from the population
{
  //short i, c, ti, ag, cl, tr, j;
  //double tenter, tleave, tbirth, sum, z, probif, probfi, sumst[ncomp], tmi, tma, bval, lval, fra, al, Ltoti, Ltotf, dLpri;;

  j = 0; for(i = 0; i < nindtot; i++) if(indev[i][0] == 0) j++;

  if(ran() < 0.5){     // adds an individual#
    //cout << "add\n";
    ntr_add++;
    i = nindtot;

    probfi = log(1.0/(nindtot-nind+1));
    probif = 0;

    fra = nuoinit/(nuoinit+nuosource);
    if(sourcefl == 0 || ran() < fra){
      if(sourcefl != 0) probif += log(fra);

      sum = 0; for(c = 0; c < ncompswa; c++){ sum += induoinit[c]; sumst[c] = sum;}
      z = ran()*sum; cnew = 0; while(cnew < ncompswa && z > sumst[cnew]) cnew++; if(cnew == ncompswa) emsg("ncomps");
      probif += log(induoinit[cnew]/sum);

      tent = 0;
    }
    else{
      probif += log(1-fra);
      sum = 0; for(c = 0; c < ncompswf; c++){ sum += induoent[c]; sumst[c] = sum;}
      z = ran()*sum; cnew = 0; while(cnew < ncompswf && z > sumst[cnew]) cnew++; if(cnew == ncompswf) emsg("ncompsswf");

      probif += log(induoent[cnew]/sum);

      ti = (cnew/classmult[settimecl])%nclassval[settimecl];

      if(ti == 0) tmi = 0; else tmi = settime[ti-1];
      if(ti == nsettime) tma = tmax; else tma = settime[ti];

      tenter = tmi + ran()*(tma-tmi);
      probif += log(1.0/(tma-tmi));
    }

    bval = 0;
    if(tbirthfl == 1){
      ag = (cnew/classmult[agecl])%nclassval[agecl];

      if(ag == nage-1){
        bval = -log(ran())*induotbirthtime;
        probif += -log(induotbirthtime) - bval/induotbirthtime;
        tbirth = tenter - age[ag] - bval;
      }
      else{
        tmi = tenter - age[ag+1]; tma = tenter - age[ag];
        tbirth = tmi + ran()*(tma-tmi);
        probif += log(1.0/(tma-tmi));
      }
    }
    else tbirth = tenter;

    lval = 0;
    fra = nuoend/(nuoend+nuosink);
    if(sinkfl == 0 || ran() < fra){
      tleave = tmax;
      if(sinkfl == 1) probif += log(fra);
    }
    else{
      probif += log(1-fra);
      lval = -log(ran())*induolifetime; tleave = tenter + lval; if(tleave >= tmax){ nfa2_add++; return;}
      probif += -log(induolifetime) - lval/induolifetime;
    }

    indsiminit(0,i,tbirth,tenter,tleave);

    probif += indsim(i); if(probif < -large/2){ nfa_add++; return;}

    for(cl = 0; cl < nclass; cl++) statinit[cl] = (evinit_fin/classmult[cl])%nclassval[cl];

    addemptyind(statinit,tbirth);

    if(tbirthfl == 1) dLpri = logtbirthrate - (tenter-tbirth)*tbirthrate; else dLpri = 0;

    if(resampleind(i,dLpri + probfi - probif,tenter,tleave)==1){
      nac_add++;
      Lpri += dLpri;
      if(samp < burnin){    //  tunes proposals
        if(tenter == 0){
          induoinit[cnew%ncompswa]++; nuoinit++;
        } 
        else{
          if(bval != 0) induotbirthtime = 0.999*induotbirthtime + 0.001*bval;
          induoent[cnew]++; nuosource++;
        }

        if(tleave == tmax) nuoend++;
        else{
          if(lval != 0) induolifetime = 0.999*induolifetime + 0.001*lval;
          nuosink++;
        }
      }
    }
    else rememptyind(i);
  }
  else{     // removes an event
    //cout << "rem\n";
    ntr_rem++;
    if(nindtot == nind){ nfa_rem++; return;}

    probif = log(1.0/(nindtot-nind));
    probfi = 0;

    i = nind+short(ran()*(nindtot-nind));

    if(indinit[i][agecl] == -1) tenter = indevt[i][agecl][0]; else tenter = 0;
    if(nindev[i][agecl] > 0 && transf[agecl][indev[i][agecl][nindev[i][agecl]-1]] == -1) tleave = indevt[i][agecl][nindev[i][agecl]-1]; else tleave = tmax;
    tbirth = indtbirth[i];

    cnew = 0; for(cl = 0; cl < nclass-2; cl++) cnew += indinit[i][cl]*classmult[cl];
    indsiminit(1,i,tbirth,tenter,tleave);

    fra = nuoinit/(nuoinit+nuosource);
    if(tenter == 0){
      if(sourcefl != 0) probfi += log(fra);
      sum = 0; for(c = 0; c < ncompswa; c++) sum += induoinit[c];
      probfi += log(induoinit[cnew%ncom]/sum);
    }
    else{
      probfi += log(1-fra);
      sum = 0; for(c = 0; c < ncompswf; c++) sum += induoent[c];
      probfi += log(induoent[cnew]/sum);

      ti = (cnew/classmult[settimecl])%nclassval[settimecl];
      if(ti == 0) tmi = 0; else tmi = settime[ti-1];
      if(ti == nsettime) tma = tmax; else tma = settime[ti];
      if(tenter < tmi || tenter > tma) emsg("pro");
      probfi += log(1.0/(tma-tmi));
    }

    bval = 0;
    if(tbirthfl == 1){
      ag = (cnew/classmult[agecl])%nclassval[agecl];
      if(ag == nage-1) probfi += -log(induotbirthtime) - (tenter - age[ag] - tbirth)/induotbirthtime;
      else{
        tmi = tenter - age[ag+1]; tma = tenter - age[ag];
        if(tbirth < tmi || tbirth > tma) emsg("pro");
        probfi += log(1.0/(tma-tmi));
      }
    }

    fra = nuoend/(nuoend+nuosink);
    if(tleave == tmax){
      if(sinkfl == 1) probfi += log(fra);
    }
    else{
      probfi += log(1-fra);
      probfi += -log(induolifetime) - (tleave-tenter)/induolifetime;
    }

    probfi += probexistingevseq(i); if(probfi < -large/2) emsg("cannot be");

    if(tbirthfl == 1) dLpri = -(logtbirthrate - (tenter-tbirth)*tbirthrate); else dLpri = 0;

    Ltoti = Lpri + (Liinit + Lir - Liexp + Ligam)*invTLi + Lob_pd*invT_pd + Lob_st*invT_st + Lpop*invT_pop + dLwrong*nwrong;

    indinit_st = indinit[i][agecl];
    nindev_st = nindev[i][agecl]; for(j = 0; j < nindev_st; j++){ indev_st[j] = indev[i][agecl][j]; indevt_st[j] = indevt[i][agecl][j];}
    nindev_new = 0; indinit_new = -1;

    indchange(i,agecl);

    Ltotf = dLpri + Lpri + (Liinit + Lir - Liexp + Ligam)*invTLi + Lob_pd*invT_pd + Lob_st*invT_st + Lpop*invT_pop + dLwrong*nwrong;

    al = exp(Ltotf - Ltoti + probfi - probif); //cout << al << " " <<  Ltotf - Ltoti <<" " <<  probfi << " " << probif<< " " << al  <<  "al\n";

    if(ran() < al){
      nac_rem++;
      Lpri += dLpri;
      rememptyind(i);
    }
    else{
      indinit_new = indinit_st;
      nindev_new = nindev_st; for(j = 0; j < nindev_st; j++){ indev_new[j] = indev_st[j]; indevt_new[j] = indevt_st[j];}

      indchange(i,agecl);
    }
  }
}

void Chain::addreminit()                                     // Initialising proposals for adding and removing individuals
{
  short c;
  double sum;

  induoinit.resize(ncompswa); for(c = 0; c < ncompswa; c++) induoinit[c] = 50.0/ncompswa;

  induoent.resize(ncompswf);
  sum = 0;
  for(c = 0; c < ncompswf; c++){
    if(compiftra[NOTALIVE][c] >= 0) induoent[c] = 1; else induoent[c] = 0;
    sum += induoent[c];
  }
  for(c = 0; c < ncompswf; c++) induoent[c] *= 50/sum;

  induolifetime = tmax/10;

  nuoinit = 2; nuosource = 10;
  nuoend = 2; nuosink = 10;
}

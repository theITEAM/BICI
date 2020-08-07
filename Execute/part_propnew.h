void Chain::part_prop(long i)                                    // Makes a particle proposal on an individual
{
  long c, cc, p, li, li2, e, evpbeg, evpend, ty, ch, k, nev, ee, pp, pbeg, pend, numzero, co, fl;
  double probif, probfi, Ltoti, Ltotf, al, tma, tt, w, z, t, wi, wf, fac, facup = 1.02, facdown = 0.95, initfac;
  vector <long> bootstrap;
  long tr;

  num++;

  partcl = long(ran()*(clall+0.99999));

  ntr_part[i][partcl]++;

  timeprop[PARTINIT] -= clock();
  ist = i;

  getlifespan(i);

  indsim_init(i,indtbirth[i],1);

  nevpsec = 0; evpsec.clear(); evpref.clear();   // divides up into sections defined by observations or fixed events
  k = 0; ee = 0;
  for(e = 0; e < nevp; e++){
    t = evp[e].t;
    while(k < cexist.size() && texist[k] < t) k++;
    evpref.push_back(k);

    fl = 0;
    switch(evp[e].ty){
      case OBS: case FIXEV: case OTHCL: fl = 1; break;
      default: if(e == nevp-1) fl = 1; else fl = 0;
    }

    if(fl == 1){
      EVPSEC evpl;
      evpl.evpbeg = ee; evpl.tbeg = evp[ee].t; evpl.evpend = e; evpl.tend = evp[e].t;
      evpsec.push_back(evpl); nevpsec++;
      ee = e;
    }
  }
  bootstrap.resize(nevpsec);
  npart = long(nindpart[i][partcl]);

  part.resize(npart*nevpsec);
  wsumst.resize(npart*nevpsec);

  if(nmfl == 1) futev.clear();

  insertexisting(indev[i]);
  timeprop[PARTINIT] += clock();

  timeprop[PARTSIM] -= clock();
  erev = 0;
  for(li = 0; li < nevpsec; li++){
    evpbeg = evpsec[li].evpbeg; evpend = evpsec[li].evpend; if(evpbeg == evpend) emsg("npart");
    pbeg = li*npart; pend = pbeg+npart;

    for(p = pbeg+1; p < pend; p++) part_simsec(i,evpbeg,evpend,part[p]);

    ty = evp[evpend].ty; ch = evp[evpend].ch;
    if((ty == END && ch == 1) || ty == FIXEV || ty == OTHCL){ t = onefac*evpsec[li].tend; k = evpref[evpend]; if(k == 0) co = NOTALIVE; else co = cexist[k-1];}

    numzero = 0;
    for(p = pbeg; p < pend; p++){
      c = part[p].cend;
      w = part[p].w;
      if(w > 0){
        switch(ty){
          case OBS: k = obsprob_cope[ch][c]; if(k < 0) w = 0; else w *= exp(cope_val[k]); break;

          case FIXEV:
            if(compval[c][fixev[ch].cl] != fixev[ch].trai) w = 0;
            else{
              if(fixev[ch].like == 1){
                cc = c+fixev[ch].dc; if(cc >= ncomp) emsg("ncomp");
                tr = compiftra[c][cc];
                if(tr < 0) w = 0; else{ if(tra[tr].type == EXP_TR) w *= transra(tra[tr].eq,t,co);}
              }
            }
            break;

          case OTHCL:
            if(c+ch >= ncomp) emsg("OTHCL prob");
            tr = compiftra[c][c+ch];
            if(tr < 0) w = 0; else{ if(tra[tr].type == EXP_TR) w *= transra(tra[tr].eq,t,co);}
            break;

          case END:
            if(ch == 1){
              tr = compiftra[c][NOTALIVE];
              if(tr < 0) w = 0; else{ if(tra[tr].type == EXP_TR) w *= transra(tra[tr].eq,t,co);}
            }
            break;
        }
      }
      if(w == 0) numzero++;
      part[p].w = w;
    }
    /*
    for(p = pbeg; p < pend; p++){
      cout << p <<" " << part[p].w << " " <<  evp[evpend].t << " " << ty <<  "wei\n";
    }
    */

    //if(part[pbeg].w == 0){ for(p = pbeg; p < pend; p++) part[p].w += 0.0001;}

    if(li < nevpsec-1){
      if(numzero >= 0.5*npart){  // does a bootstrap
        bootstrap[li] = 1;
        p = pbeg+npart; part[p].pback = pbeg; part[p].w = 1;

        wsum = 0; for(p = pbeg; p < pend; p++){ wsum += part[p].w; wsumst[p] = wsum;}
        if(wsum == 0){
          nfa_part[i][partcl]++;
          if(samp < burnin){ nindpart[i][partcl] *= facup; if(nindpart[i][partcl] >= npartmax) nindpart[i][partcl] = npartmax*onefac;}
          timeprop[PARTSIM] += clock();
          return;
        }

        for(p = pbeg+npart+1; p < pend+npart; p++){
          z = ran()*wsum; pp = pbeg; while(pp < pend && z > wsumst[pp]) pp++; if(pp == pend) emsg("pro");
          part[p].pback = pp; part[p].w = 1;
        }
      }
      else{
        bootstrap[li] = 0;
        for(p = pbeg+npart; p < pend+npart; p++){
          pp = p-npart;
          part[p].pback = pp; part[p].w = part[pp].w;
        }
      }
    }
    else bootstrap[li] = 0;
  }
  timeprop[PARTSIM] += clock();

  // generates the new event sequence

  wi = part[pbeg].w;
  wsum = 0; for(p = pbeg+1; p < pend; p++){ wsum += part[p].w; wsumst[p] = wsum;}

  if(wsum == 0){ 
    nfa_part[i][partcl]++;
    if(samp < burnin){ nindpart[i][partcl] *= facup; if(nindpart[i][partcl] >= npartmax) nindpart[i][partcl] = npartmax*onefac;}
    return;
  }

  z = ran()*wsum; p = pbeg+1; while(p < pend && z > wsumst[p]) p++; if(p == pend) emsg("pro");
  wf = part[p].w;
  fac = wsum*wi/((wsum+wi-wf)*wf);

  psel.resize(nevpsec); li = nevpsec-1;
  do{
    if(bootstrap[li] == 1) fac *= part[li*npart].w/part[p].w;
    psel[li] = p; li--; p = part[p].pback;
  }while(p != -1); 
  if(li != -1) emsg("ii");

  evnew.clear();
  //c = part[psel[0]].cstart;
  //EV evbeg; evbeg.t = tent; if(tent == 0) evbeg.tr = trabeg+c; else{ evbeg.tr = compiftra[NOTALIVE][c]; if(evbeg.tr < 0) emsg("warnbeg");} evnew.push_back(evbeg);

  for(li = 0; li < nevpsec; li++){ p = psel[li]; nev = part[p].ev.size(); for(e = 0; e < nev; e++) evnew.push_back(part[p].ev[e]);}

  c = part[psel[nevpsec-1]].cend;
  EV evend; evend.t = tlea; 
  if(tlea == tmax) tr = traend+c; 
  else{ tr = compiftra[c][NOTALIVE]; if(tr < 0) emsg("warnend"); if(i < nind && indlikeleave[i] == 0) tr += moventra;}
  evend.tr = tr;
  evnew.push_back(evend);

  if(checkevdt() == 1) return;

  probif = probstart(i,evnew[0].t,tra[evnew[0].tr].cf%ncomps) + probindsim(evnew);
  probfi = probstart(i,indev[i][0].t,tra[indev[i][0].tr].cf%ncomps) + probindsim(indev[i]);
 
  Ltoti = Lwnm();
  indcha(i);
  Ltotf = Lwnm();

  al = exp(Ltotf - Ltoti + probfi - probif)*fac; //cout << al << " " << psel[0] << " " << Ltoti << " " << Ltotf << " " << probif << " " << probfi << " " << fac <<  "al\n";

  if(ran() < al){
    nac_part[i][partcl]++;
    if(psel[0] == 0) nde_part[i][partcl]++;

    initfac = 1;
    if(samp < burnin){
      if(psel[0] != 0){ nindpart[i][partcl] *= facdown; if(nindpart[i][partcl] < 2) nindpart[i][partcl] = 2;}     // if particle does not go all the way throgu then cound a failurre
      else { nindpart[i][partcl] *= facup; if(nindpart[i][partcl] >= npartmax) nindpart[i][partcl] = npartmax*onefac;}
    }
  }
  else{
    indrev(i);
    initfac = 0.2;
    if(samp < burnin){ nindpart[i][partcl] *= facup; if(nindpart[i][partcl] >= npartmax) nindpart[i][partcl] = npartmax*onefac;}
  }
  
  if(samp < burnin){
    c = tra[indev[i][0].tr].cf; if(c == NOTALIVE) emsg("pr");
    if(partcl == clall){
      if(indev[i][0].t == 0) sim_initc[i][c%ncomps] += initfac; else sim_enterc[i][c%ncomps] += initfac;
    }
    else{
      c = compval[c][partcl];
      if(indev[i][0].t == 0) simcl_initc[i][partcl][c] += initfac; else simcl_enterc[i][partcl][c] += initfac;
    }
  }
}

void Chain::insertexisting(vector<EV> &ev)                                      // populates particles with the existing sequence
{
  long e, ee, nev, li, li2, p, c, cl, j;
  long tr, trr;
  double t, tt, ttt, tma, mean, kshape, lam;

  nev = ev.size();
  tr = ev[0].tr; tt = ev[0].t; c = tra[tr].cf; part[0].cstart = c;
  if(tt == 0) part[0].w = 1; else part[0].w = transra(tra[tr].eq,tt,NOTALIVE);

  part[0].pback = -1;
  e = 0;
  for(li = 0; li < nevpsec; li++){
    p = li*npart;
    tma = evpsec[li].tend;
    part[p].ev.clear();
    while(e < nev && ev[e].t < tma){
      tr = ev[e].tr; t = ev[e].t;

      if(nmfl == 1 && li > 0){   // Adds in the future
        if(tra[tr].nm == 1 && !(tra[tr].capev >= 0 && capevall[tra[tr].capev] == 1)){
          cl = tra[tr].cl;
          if(partcl == clall && cl == partcl){
            ee = e-1; while(ee > 0 && tra[ev[ee].tr].cl != cl) ee--;
            tt = ev[ee].t;
  
            if(tt < evpsec[li].tbeg){
              li2 = li-1; while(li2 >= 0 && tt < evpsec[li2].tend){ addexistingfute(part[li2*npart],tr,tt,t); li2--;};
            }
          }
        }

        if(tra[tr].ntraend > 0){
          for(j = 0; j < tra[tr].ntraend; j++){
            trr = tra[tr].traend[j];
            cl = tra[trr].cl;
            ee = e-1; while(ee > 0 && tra[ev[ee].tr].cl != cl) ee--;
            tt = ev[ee].t;
            if(tt < evpsec[li].tbeg){ // Need to sample from
              switch(tra[trr].type){
                case GAMMA_TR:
                  mean = nmeq_val[tra[trr].eq]; kshape = nmeq_val[tra[trr].eqshape];
                  ttt = tt + gammasamptail(kshape,kshape/mean,t-tt);
                  break;

                case WEI_TR:
                  lam = nmeq_val[tra[trr].eq]; kshape = nmeq_val[tra[trr].eqshape];
                  ttt = tt + weibullsamptail(lam,kshape,t-tt);
                  break;

                case FIXED_TR: emsg("prob fi"); break;
              }
              if(ttt < t) emsg("order pp");
              li2 = li-1; while(li2 >= 0 && tt < evpsec[li2].tend){ addexistingfute(part[li2*npart],tr,tt,ttt); li2--;};
            }
          }
        }
      }

      c = tra[tr].cf; if(c == NOTALIVE) emsg("pr");
      part[p].ev.push_back(ev[e]);
      e++;
    }

    /*
    for(li2 = evpsec[li].evpbeg; li2 < evpsec[li].evpend; li2++){   // looks at then case when there is an event on the boundary
      if(evp[li2].t == tma && evp[li2].ty == CHC){
        if(e >= nev) emsg("prob h");
        if(tra[ev[e].tr].cf != c) emsg("prob");
        c = tra[ev[e].tr].cf;
        part[p].ev.push_back(ev[e]); e++;
      }
    }
    */

    part[p].cend = c;
    if(nmfl == 1) part[p].futerefst.clear();
  }
  if(e != nev-1) emsg("over");

  if(checkon == 1){
    for(li = 0; li < nevpsec; li++){
      p = li*npart; e = evpsec[li].evpbeg; if(evp[e].ty == CHC){ if(tra[part[p].ev[0].tr].dc != evp[e].ch) emsg("dc");}
    }
  }
}

void Chain::addexistingfute(PART &pa, long tr, double tbeg, double tend)
{
  long m;

  nfuteref = pa.futerefst.size();
  m = 0; while(m < nfuteref && tend < futev[pa.futerefst[m]].t) m++;
  pa.futerefst.insert(pa.futerefst.begin()+m,futev.size());

  FUTEV fev; fev.tr = tr; fev.tnow = tbeg; fev.t = tend; futev.push_back(fev);
}

void Chain::part_simsec(long i, long evpbeg, long evpend, PART &pa)     // simulates new sequence
{
  long e, k, co, j, ch, c, cf, p, m, n, fl, a, fev, cl;
  double t, tt, tst, tnext, tnextold, z, R, r, val, kshape, futnext, tmid;
  long tr;

  t = evp[evpbeg].t;

  evnew.clear();

  if(evpbeg == 0){                 // Initialise
    c = selectstart(i,t);
    c += ctop;
    pa.cstart = c;
    pa.pback = -1;
    if(t == 0){ pa.w = 1; tr = trabeg+c;} 
    else{ 
      tr = compiftra[NOTALIVE][c]; if(tr < 0){ pa.w = 0; return;};
      if(i < nind && indlikeenter[i] == 0) tr += moventra;
      else pa.w = transra(tra[tr].eq,t,NOTALIVE);
    }
    if(nmfl == 1){ nfuteref = 0; futeref.clear();}
    addevnew(tr,t);
  }
  else{                            // Does the particle filter
    if(pa.w == 0) return;
    p = pa.pback;
    c = part[p].cend;
    getnotdep(compleavesimclnotdep[partcl][c]);
    if(nmfl == 1){ futeref = part[p].futerefst; nfuteref = futeref.size();}
  }

  k = evpref[evpbeg]; tnextold = texist[k]; if(k == 0) co = NOTALIVE; else co = cexist[k-1];

  for(e = evpbeg; e < evpend; e++){
    switch(evp[e].ty){
      case CHC: case OTHCL:
        tr = compiftra[c][c + evp[e].ch]; if(tr < 0){ pa.w = 0; return;};
        c = addevnew(tr,t);
        break;

      case FIXEV:
        fev = evp[e].ch; cl = fixev[fev].cl; if(compval[c][cl] != fixev[fev].trai){ pa.w = 0; return;}
        tr = compiftra[c][c + fixev[fev].dc]; if(tr < 0){ pa.w = 0; return;};
        if(fixev[fev].like == 0) tr += moventra;
        c = addevnew(tr,t);
        break;
    }

    tnext = evp[e+1].t;
    if(t < tnext){
      do{
        if(compleavesimcldep[partcl][c].size() > 0){
          tt = t+partdt; if(tt > tnext) tt = tnext;
          tmid = 0.5*(t+tt); while(tnextold < tmid){ co = cexist[k]; k++; if(k < nexist) tnextold = texist[k]; else tnextold = tmax;}
          getdep(compleavesimcldep[partcl][c],co,tmid);
        }
        else{ Rdep = 0; tt = tnext;}

        R = Rdep + Rnotdep;

        if(nmfl == 1){ if(nfuteref > 0){ futnext = futev[futeref[nfuteref-1]].t; if(tt > futnext) tt = futnext;} else futnext = large;} else futnext = large;

        if(R == 0) t = tt;
        else{
          tst = t;
          t += -log(ran())/R;
          if(t > tt) t = tt;
          else{ // actual event
            z = ran()*R;
            if(z <= Rdep){
              j = 0; while(j < nRdep && z > Rdepst[j]) j++; if(j == nRdep) emsg("nRdep");
              r = rdep[j]; tr = compleavesimcldep[partcl][c][j];
            }
            else{
              z -= Rdep;
              j = 0; while(j < nRnotdep && z > Rnotdepst[j]) j++; if(j == nRnotdep){ 
                stringstream ss;
                for(j = 0; j < nRdep; j++) ss <<  Rdepst[j] << "   "; ss << "dep  ";
                for(j = 0; j < nRnotdep; j++) ss << Rnotdepst[j] << "   "; ss << "notdep  ";

                ss  << z << " " << Rdep << " " << Rnotdep << " R ";
                cout << "z|" << ss.str() << "|\n"; emsg(ss.str(),simnum,num);}
              r = rnotdep[j]; tr = compleavesimclnotdep[partcl][c][j];
            }
            //prob += log(r);
            c = addevnew(tr,t);
          }
          //prob -= (t-tst)*R;
        }

        if(t == futnext){   // Future event
          n = futeref[nfuteref-1]; //futeref.pop_back(); nfuteref--;
          tr = compiftra[c][c + tra[futev[n].tr].dc]; if(tr < 0){ pa.w = 0; return;};
          c = addevnew(tr,t);
        }
      }while(t < tnext);
    }
  }

  pa.cend = c;
  pa.ev = evnew;
  if(nmfl == 1 && evpend < nevp-1) pa.futerefst = futeref;
}

void Chain::chainpartinit()             // Initialises quantities used when sampling from partivles
{
  long ob, fi, c, a, enterpos[ncomps], j, cl;
  long i;

  nindpart.resize(nind); 
  for(i = 0; i < nind; i++){ nindpart[i].resize(clall+1); for(cl = 0; cl <= clall; cl++) nindpart[i][cl] = npartstart;}

  for(c = 0; c < ncomps; c++) enterpos[c] = 0; //  works out if an individual can enter into a particular compartment
  for(c = 0; c < ncomp; c++) if(compiftra[NOTALIVE][c] >= 0) enterpos[c%ncomps] = 1;

  sim_initcsumtot.resize(nind); sim_initc.resize(nind); sim_initcsum.resize(nind); logsim_initc.resize(nind); 
  sim_entercsumtot.resize(nind); sim_enterc.resize(nind); sim_entercsum.resize(nind); logsim_enterc.resize(nind); 
  for(i = 0; i < nind; i++){
    sim_initc[i].resize(ncomps); for(c = 0; c < ncomps; c++) sim_initc[i][c] = 400.0/ncomps;
    sim_initcsum[i].resize(ncomps); logsim_initc[i].resize(ncomps);
    sim_enterc[i].resize(ncomps); for(c = 0; c < ncomps; c++){ if(enterpos[c] == 1) sim_enterc[i][c] = 400.0/ncomps; else sim_enterc[i][c] = 0;}
    sim_entercsum[i].resize(ncomps); logsim_enterc[i].resize(ncomps);
  }

  simuo_initcsumtot.resize(naddremdiv+1); simuo_initc.resize(naddremdiv+1); simuo_initcsum.resize(naddremdiv+1); logsimuo_initc.resize(naddremdiv+1);
  simuo_entercsumtot.resize(naddremdiv+1); simuo_enterc.resize(naddremdiv+1); simuo_entercsum.resize(naddremdiv+1); logsimuo_enterc.resize(naddremdiv+1); 
  for(a = 0; a <= naddremdiv; a++){
    simuo_initc[a].resize(ncomps); for(c = 0; c < ncomps; c++) simuo_initc[a][c] = 400.0/ncomps;
    simuo_initcsum[a].resize(ncomps); logsimuo_initc[a].resize(ncomps);
    simuo_enterc[a].resize(ncomps); for(c = 0; c < ncomps; c++){ if(enterpos[c] == 1) simuo_enterc[a][c] = 400.0/ncomps; else simuo_enterc[a][c] = 0;}
    simuo_entercsum[a].resize(ncomps); logsimuo_enterc[a].resize(ncomps);
  }

  simcl_initcsumtot.resize(nind); simcl_initc.resize(nind); simcl_initcsum.resize(nind); logsimcl_initc.resize(nind);
  simcl_entercsumtot.resize(nind); simcl_enterc.resize(nind); simcl_entercsum.resize(nind); logsimcl_enterc.resize(nind);
  for(i = 0; i < nind; i++){
    simcl_initcsumtot[i].resize(clall); simcl_initc[i].resize(clall); simcl_initcsum[i].resize(clall); logsimcl_initc[i].resize(clall);
    simcl_entercsumtot[i].resize(clall); simcl_enterc[i].resize(clall); simcl_entercsum[i].resize(clall); logsimcl_enterc[i].resize(clall);
    for(cl = 0; cl < clall; cl++){
      simcl_initc[i][cl].resize(nclassval[cl]); for(c = 0; c < nclassval[cl]; c++) simcl_initc[i][cl][c] = 400.0/nclassval[cl];
      simcl_initcsum[i][cl].resize(nclassval[cl]); logsimcl_initc[i][cl].resize(nclassval[cl]);
      simcl_enterc[i][cl].resize(nclassval[cl]); for(c = 0; c < nclassval[cl]; c++) simcl_enterc[i][cl][c] = 400.0/nclassval[cl];
      simcl_entercsum[i][cl].resize(nclassval[cl]); logsimcl_enterc[i][cl].resize(nclassval[cl]);
    }
  }
}

void Chain::simsumcalc()                         // Used to initialise sampling of initial compartment
{
  long c, a, ag, cl;
  long i;

  for(i = 0; i < nind; i++){
    sim_initcsumtot[i] = 0;
    for(c = 0; c < ncomps; c++){ sim_initcsumtot[i] += sim_initc[i][c]; sim_initcsum[i][c] = sim_initcsumtot[i];}
    if(sim_initcsumtot[i] == 0){ for(c = 0; c < ncomps; c++) logsim_initc[i][c] = -large;}
    else{ for(c = 0; c < ncomps; c++) logsim_initc[i][c] = log(sim_initc[i][c]/sim_initcsumtot[i]);}

    sim_entercsumtot[i] = 0;
    for(c = 0; c < ncomps; c++){ sim_entercsumtot[i] += sim_enterc[i][c]; sim_entercsum[i][c] = sim_entercsumtot[i];}
    if(sim_entercsumtot[i] == 0){ for(c = 0; c < ncomps; c++) logsim_enterc[i][c] = -large;}
    else{ for(c = 0; c < ncomps; c++) logsim_enterc[i][c] = log(sim_enterc[i][c]/sim_entercsumtot[i]);}

    for(cl = 0; cl < clall; cl++){
      simcl_initcsumtot[i][cl] = 0;
      for(c = 0; c < nclassval[cl]; c++){ simcl_initcsumtot[i][cl] += simcl_initc[i][cl][c]; simcl_initcsum[i][cl][c] = simcl_initcsumtot[i][cl];}
      if(simcl_initcsumtot[i][cl] == 0){ for(c = 0; c < nclassval[cl]; c++) logsimcl_initc[i][cl][c] = -large;}
      else{ for(c = 0; c < nclassval[cl]; c++) logsimcl_initc[i][cl][c] = log(simcl_initc[i][cl][c]/simcl_initcsumtot[i][cl]);}

      simcl_entercsumtot[i][cl] = 0;
      for(c = 0; c < nclassval[cl]; c++){ simcl_entercsumtot[i][cl] += simcl_enterc[i][cl][c]; simcl_entercsum[i][cl][c] = simcl_entercsumtot[i][cl];}
      if(simcl_entercsumtot[i][cl] == 0){ for(c = 0; c < nclassval[cl]; c++) logsimcl_enterc[i][cl][c] = -large;}
      else{ for(c = 0; c < nclassval[cl]; c++) logsimcl_enterc[i][cl][c] = log(simcl_enterc[i][cl][c]/simcl_entercsumtot[i][cl]);}
    }
  }

  for(a = 0; a <= naddremdiv; a++){
    simuo_initcsumtot[a] = 0;
    for(c = 0; c < ncomps; c++){ simuo_initcsumtot[a] += simuo_initc[a][c]; simuo_initcsum[a][c] = simuo_initcsumtot[a];}
    if(simuo_initcsumtot[a] == 0){ for(c = 0; c < ncomps; c++) logsimuo_initc[a][c] = -large;}
    else{ for(c = 0; c < ncomps; c++) logsimuo_initc[a][c] = log(simuo_initc[a][c]/simuo_initcsumtot[a]);}

    simuo_entercsumtot[a] = 0;
    for(c = 0; c < ncomps; c++){ simuo_entercsumtot[a] += simuo_enterc[a][c]; simuo_entercsum[a][c] = simuo_entercsumtot[a];}
    if(simuo_entercsumtot[a] == 0){ for(c = 0; c < ncomps; c++) logsimuo_enterc[a][c] = -large;}
    else{ for(c = 0; c < ncomps; c++) logsimuo_enterc[a][c] = log(simuo_enterc[a][c]/simuo_entercsumtot[a]);}
  }

  addremdivsumtot = 0; for(a = 0; a <= naddremdiv; a++){ addremdivsumtot += addremdiv[a]; addremdivsum[a] = addremdivsumtot;}
  if(addremdivsumtot == 0){ for(a = 0; a <= naddremdiv; a++) logaddremdiv[a] = -large;}
  else{ for(a = 0; a <= naddremdiv; a++) logaddremdiv[a] = log(addremdiv[a]/addremdivsumtot);}

  for(a = 0; a <= naddremdiv; a++){
    induoagesumtot[a] = 0;  for(ag = 0; ag < nclassval[agecl]; ag++){ induoagesumtot[a] += induoage[a][ag]; induoagesum[a][ag] = induoagesumtot[a];}
    if(induoagesumtot[a] == 0){ for(ag = 0; ag < nclassval[agecl]; ag++) loginduoage[a][ag] = -large;}
    else{ for(ag = 0; ag < nclassval[agecl]; ag++) loginduoage[a][ag] = log(induoage[a][ag]/induoagesumtot[a]);}
  }
}

void partinit()                                 // Initialises quantities used when sampling from particles
{
  long eq, c, cc, cl, i, f, ci, cf, p, j, d, depmax, ndmax;
  long tr, tr2;
  vector <long> nnonexpc;                             // [c][#] gives the non exponential transition activated when individual enters c
  vector< vector <long> > nonexpc;

  //compleavesimdep.resize(ncomp+1); compleavesimnotdep.resize(ncomp+1);    // Used when simulating particles
  compleavesimcldep.resize(clall+1); compleavesimclnotdep.resize(clall+1);
  for(cl = 0; cl <= clall; cl++){ compleavesimcldep[cl].resize(ncomp+1); compleavesimclnotdep[cl].resize(ncomp+1);}

  for(tr = 0; tr < ntra; tr++){
     if(tra[tr].type == EXP_TR && tra[tr].cf != NOTALIVE && !(tra[tr].capev >= 0 && capevall[tra[tr].capev] == 1) && tra[tr].like == 1){
      eq = tra[tr].eq;
      switch(transdep[eq]){
        case 0: compleavesimclnotdep[clall][tra[tr].ci].push_back(tr); break;
        case 1: compleavesimcldep[clall][tra[tr].ci].push_back(tr); break;
      }

      cl = tra[tr].cl;
      if(cl >= 0 && cl < clall){
        switch(transdep[eq]){
          case 0: compleavesimclnotdep[cl][tra[tr].ci].push_back(tr); break;
          case 1: compleavesimcldep[cl][tra[tr].ci].push_back(tr); break;
        }
      }
    }
  }
  //ncompleavesimdep.resize(ncomp+1); for(c = 0; c <= ncomp; c++) ncompleavesimdep[c] = compleavesimdep[c].size();
  //ncompleavesimnotdep.resize(ncomp+1); for(c = 0; c <= ncomp; c++) ncompleavesimnotdep[c] = compleavesimnotdep[c].size();

  compleavedep.resize(ncomp+1); compleavenotdep.resize(ncomp+1);    // Used when simulating in sim.h
  for(tr = 0; tr < ntra; tr++){
    if(tra[tr].type == EXP_TR && tra[tr].like == 1){
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

  for(p = 0; p < npopnum; p++){                                                               // works out which rates need to be recalculated
    popnumtake.push_back(vector<double>());
    for(c = 0; c <= ncomp; c++) popnumtake[p].push_back(0);
    for(j = 0; j < popnumterm[p].size(); j++) popnumtake[p][popnumterm[p][j]] -= popnumtermweight[p][j];
  }

  for(d = 0; d < ntransdepeq; d++){
    transdepeqrecalc.push_back(vector<long>());
    eq = transdepeq[d];
    for(c = 0; c < ncomp; c++){
      for(p = 0; p < neq_popnum[eq]; p++) if(popnumtake[eq_popnum[eq][p]][c] != 0) break;
      if(p < neq_popnum[eq]) transdepeqrecalc[d].push_back(1);
      else transdepeqrecalc[d].push_back(0);
    }
  }

  nonexpc.resize(ncomp+1); nnonexpc.resize(ncomp+1);
  for(tr = 0; tr < ntra; tr++){
    if(tra[tr].like == 1){
      switch(tra[tr].type){
        case FIXED_TR: case GAMMA_TR: case WEI_TR:
          nonexpc[tra[tr].ci].push_back(tr);
          break;
      }
    }
  }
  for(c = 0; c <= ncomp; c++) nnonexpc[c] = nonexpc[c].size();

  //nonexptra.resize(ntra); nnonexptra.resize(ntra);
  for(tr = 0; tr < ntra; tr++){
    cf = tra[tr].cf;
    if(cf != NOTALIVE){
      for(j = 0; j < nnonexpc[cf]; j++){
        tr2 = nonexpc[cf][j];
        //if(tra[tr2].cl == tra[tr].cl || tra[tr].cl < 0) nonexptra[tr].push_back(tr2);
        if(tra[tr2].cl == tra[tr].cl || tra[tr].cl < 0) tra[tr].nonexptra.push_back(tr2);
      }
    }
  }
  for(tr = 0; tr < ntra; tr++) tra[tr].nnonexptra = tra[tr].nonexptra.size();

  for(tr = 0; tr < ntra; tr++){
    tra[tr].ntraend = 0;

    ci = tra[tr].ci;
    if(ci != NOTALIVE){
      for(j = 0; j < nnonexpc[ci]; j++){
        tr2 = nonexpc[ci][j];
        if(tr2 != tr && (tra[tr2].cl == tra[tr].cl || tra[tr].cl < 0)){
          tra[tr].traend.push_back(tr2);
          tra[tr].ntraend++;
        }
      }
    }
  }

  /*
  nonexptrain.resize(ntra); nnonexptrain.resize(ntra);
  for(tr = 0; tr < ntra; tr++){
    ci = tra[tr].ci;
    if(ci != NOTALIVE){
      for(j = 0; j < nnonexpc[ci]; j++){
        tr2 = nonexpc[ci][j];
        if(tra[tr2].cl == tra[tr].cl || tra[tr].cl < 0) nonexptrain[tr].push_back(tr2);
      }
    }
  }
  for(tr = 0; tr < ntra; tr++) nnonexptrain[tr] = nonexptrain[tr].size();
        */
            
  /*
  fixnmlist.clear();
  fixnmmap.resize(nclass);
  for(cl = 0; cl < nclass; cl++){
    fixnmmap[cl].resize(nclassval[cl]);
    for(i = 0; i < nclassval[cl]; i++)  fixnmmap[cl][i].resize(nclassval[cl]);
  }
  */
}

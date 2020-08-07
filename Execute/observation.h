void observationinit()                                                                               // Initialises the observation model
{
  long i, j, ob, c, eq, jj, cap, div, divi, p, dm;

  obsprob_cope.resize(nobs);
  for(i = 0; i < nobs; i++) obsprob_cope[i].resize(ncomp+1);

  ncope = 0;
  for(i = 0; i < nind; i++){
    for(j = 0; j < indobs[i].size(); j++){
      ob = indobs[i][j];
      for(c = 0; c <= ncomp; c++){
        obsprob_cope[ob][c] = -1;
        if(c < ncomp){
          eq = obsprobeqn[ob][c];
          if(eq >= 0){
            jj = 0; while(jj < ncope && cope[jj] != eq) jj++;
            if(jj == ncope){ cope.push_back(eq); ncope++;}
            obsprob_cope[ob][c] = jj;
          }
        }
      }
    }
  }

  capprob_cpe.resize(ncap);
  for(i = 0; i < ncap; i++) capprob_cpe[i].resize(ncomp);

  ncpe = 0;    // makes a list of equations relating to capture probability
  for(cap = 0; cap < ncap; cap++){
    for(c = 0; c < ncomp; c++){
      eq = capprobeqn[cap][c];
      capprob_cpe[cap][c] = -1;
      if(eq >= 0){
        i = 0; while(i < ncpe && cpe[i] != eq) i++;
        if(i == ncpe){ cpe.push_back(eq); ncpe++;}
        capprob_cpe[cap][c] = i;
      }
    }
  }

  compcap.resize(ncomp); ncompcap.resize(ncomp);
  for(c = 0; c < ncomp; c++){         // lists captures that happen in each compartment
    //compcap.push_back(vector<long>());
    //compcapty.push_back(vector<long>());
    //compcapt.push_back(vector<double>());
    ncompcap[c] = 0;
    for(j = 0; j < ncap; j++){
      if(capprobeqn[j][c] >= 0){
        CC cc; cc.val = j; cc.ty = CAP; cc.t = capt[j]; compcap[c].push_back(cc); ncompcap[c]++;
        /*
        compcap[c].push_back(j);
        compcapty[c].push_back(CAP);
        compcapt[c].push_back(capt[j]);
        */
      }
    }

    for(p = 0; p < npopm; p++){
      if(popmcomp[p][c] == 1){
        CC cc; cc.val = p; cc.ty = POP; cc.t = popmt[p]; compcap[c].push_back(cc); ncompcap[c]++;
      }
    }

    for(dm = 0; dm < nderm; dm++){
      if(dermcomp[dm][c] == 1){
        CC cc; cc.val = dm; cc.ty = DER; cc.t = dermt[dm]; compcap[c].push_back(cc); ncompcap[c]++;
      }
    }

    sort(compcap[c].begin(),compcap[c].end(),comparecc);

    CC cc; cc.t = tmax; compcap[c].push_back(cc);
  }

  compcapdiv.resize(ncomp);
  for(c = 0; c < ncomp; c++){                    // creates a lookup table to calculate fast
    compcapdiv[c].resize(ncompcapdiv); 
    for(div = 0; div < ncompcapdiv; div++) compcapdiv[c][div] = ncompcap[c];
    divi = 0;
    for(j = 0; j < ncompcap[c]; j++){
      div = long(ncompcapdiv*compcap[c][j].t/tmax);
      while(divi <= div){ compcapdiv[c][divi] = j; divi++;}
    }
  }

  for(cap = 0; cap < ncap; cap++){
    capindob.push_back(vector<long>());
    for(i = 0; i < nind; i++) capindob[cap].push_back(-1);
  }

  for(i = 0; i < nind; i++){
    for(j = 0; j < nindobs[i]; j++){
      ob = indobs[i][j];
      c = obscap[ob];
      if(c >= 0) capindob[c][i] = ob;
    }
  }

  for(p = 0; p < npopm; p++) popmap.push_back(0);

  for(p = 0; p < nderm; p++) dermap.push_back(0);
}

void Chain::chainobsinit()                       // Initialises quantities relating to the observation probability
{
  long j, ob, fi;
  double pd;

  for(j = 0; j < ncope; j++){
    cope_num.push_back(0);
    cope_op.push_back(calcobsprob(cope[j]));
    cope_val.push_back(log(calcobsprob(cope[j])));
  }

  for(j = 0; j < ncpe; j++){
    cpe_num.push_back(0); cpe_oneminus.push_back(0);
    cpe_pd.push_back(calcobsprob(cpe[j]));
    cpe_val.push_back(log(calcobsprob(cpe[j])));
    cpe_valoneminus.push_back(log(1-calcobsprob(cpe[j])));
  }

  for(j = 0; j < ncapev; j++){
    capev_num.push_back(0); capev_oneminus.push_back(0);
    pd = calcobsprob(capevprobeqn[j]);
    capev_val.push_back(log(pd));
    capev_valoneminus.push_back(log(1-pd));
  }
}

double Chain::calcobsprob(long eq)       // Calculates the observation probaability
{
  double pr;

  if(eq < 0) emsg("calc obsprob");
  pr = calculatenotdep(eq,param);
  if(pr < obssmall) pr = obssmall; if(pr > 1-obssmall) pr = 1-obssmall;

  return pr;
}

double popcalc(long p, double val)          // Calculates the observation probability for a population
{
  return gammaprob(val+0.0001,popmal[p],popmbe[p]);
}

double dercalc(short j, double *popnum, double *param)        // Calculates the observation probability for a derived quantity
{
  long d;
  d = derm[j];
  return normalprob(calculate(derive[d],popnum,param),dermav[j],dermvar[j]);
}

long Chain::numwrong(long i, vector<EV> &vec)
{
  long numwrong, c, e, nev, ob, j;
  double t;

  nwrong = 0;

  c = NOTALIVE; e = 0; nev = vec.size();
  for(j = 0; j < nindobs[i]; j++){
    ob = indobs[i][j];

    t = obst[ob];
    while(e < nev && vec[e].t < t){ c = tra[vec[e].tr].cf; e++;}
    if(c == NOTALIVE) nwrong++;
    else{
      if(obsprobeqn[ob][c] == -1) nwrong++;
    }
  }
  return nwrong;
}

void altertemperature()
{
  long p, j;
  double dd, f;

  if(ch[0]->invT_pop > 1){  // Looks to relax population measurements
    if(ch[0]->invT_pop > 5) f = 0.5; else f = 2;

    for(p = 0; p < npopm; p++){
      dd = ch[0]->pop[p] - popmal[p]/popmbe[p];
      if(sqrt(dd*dd/(popmal[p]/(popmbe[p]*popmbe[p]))) > 0.5) break;
    }
    if(p == npopm){ ch[0]->invT_pop--; if(ch[0]->invT_pop < 1) ch[0]->invT_pop = 1;}

    if(samp > burnin/2) ch[0]->invT_pop = 1;
  }

  if(ch[0]->invT_der > 1){  // Looks to relax derived measurements
    for(j = 0; j < nderm; j++){
      dd = calculate(derive[derm[j]],ch[0]->dermpopnum[j],ch[0]->param) - dermav[j]; 
      if(sqrt(dd*dd/dermvar[j]) > 1) break;
    }
    if(j == nderm || samp > burnin) ch[0]->invT_der = 1;
  }
}


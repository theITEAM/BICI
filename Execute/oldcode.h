/*
void Chain::sim_prop(long i)                                                                 // Makes proposal by simulating a new sequence 
{
  short e, c;
  double probif, probfi, Ltoti, Ltotf, al;

  ist = i;
  evrev = indev[i];

  getlifespan(i);

  part_init(i,indtbirth[i],2);
  if(nevpsec != 1) emsg("should be 1");

  insertexisting();

  probif = part_simsec(i,0,nevp-1,part[1]);
  psel.resize(1); psel[0] = 1;

  generatenew();

  //cout << i << " i\n"; outputevents("J",evnew);

  probfi = part_simprob(i,0,nevp-1,part[0]);

  //if(i == 0) outputevents("Y",evnew);

  Ltoti = Lwnm();
  indchange(i);
  Ltotf = Lwnm();

  if(samp < 100) al = exp(Ltotf - Ltoti);
  else al = exp(Ltotf - Ltoti + probfi - probif); 
  //if(i == 0)  cout << al << " " << Ltoti << " " << Ltotf << " " << probif << " " << probfi << "al\n";

  ntr_sim[i]++;
  if(ran() < al){
    nac_sim[i]++;
   // cout <<part[pi].cstart << " " << part[pf].cstart << " " << tra[evnew[0].tr].cf <<  " c\n";

    if(samp < burnin){
      c = tra[indev[i][0].tr].cf; if(c == NOTALIVE) emsg("pr");
      if(indev[i][0].t == 0) sim_initc[i][c%ncomps]++; else sim_enterc[i][c%ncomps]++;
}
}
  else{
    evnew = evrev;
    indchange(i);
}
  
}
*/



void Chain::generatenew()
{
  short c, li, p, e, nev;

  evnew.clear();
  c = part[psel[0]].cstart;
  EV evbeg; evbeg.t = tent; if(tent == 0) evbeg.tr = trabeg+c; else{ evbeg.tr = compiftra[NOTALIVE][c]; if(evbeg.tr < 0) emsg("warnbeg");} evnew.push_back(evbeg);

  for(li = 0; li < nevpsec; li++){ p = psel[li]; nev = part[p].ev.size(); for(e = 0; e < nev; e++) evnew.push_back(part[p].ev[e]);}

  c = part[psel[nevpsec-1]].cend;
  EV evend; evend.t = tlea; if(tlea == tmax) evend.tr = traend+c; else{ evend.tr = compiftra[c][NOTALIVE]; if(evend.tr < 0) emsg("warnend");}
  e = evnew.size(); while(e > 0 && evnew[e-1].t > tmax) e--; evnew.insert(evnew.begin()+e,evend);
}




/*
double Chain::probsimple_sim(long i)                                   // Simulates using particles
{
  short c, e;
  double prob;

  part_init(0,i,indtbirth[i],1);
  if(nevpsec != 1) emsg("should be 1");

  ist = i;
  prob = part_simsec(i,0,nevp-1,part[0]);

  psel.resize(1); psel[0] = 0; generatenew();
  return prob;
}
    */
        
double Chain::part_simprob(long i, short evpbeg, short evpend, PART &pa)     // simulates new sequence
{
  short e, k, co, j, ch, c, cf, p, m, n, fl, ee, nev, a;
  double t, tt, tst, tnext, tnextold, z, R, r, val, kshape, futnext, tbeg, prob = 0, tmid;
  long tr;
  vector <EV> ev;

  tbeg = evp[evpbeg].t; t = tbeg;

  if(evpbeg == 0){
    c = pa.cstart;
    prob += probstart(i,tbeg,c%ncomps);
  }
  else c = part[pa.pback].cend;

  k = evpref[evpbeg]; if(k >= 0){ if(k == 0) co = NOTALIVE; else co = cexist[k-1]; tnextold = texist[k];}else{ tnextold = tmax; co = NOTALIVE;}
  ev = pa.ev; ee = 0; nev = ev.size();

  for(e = evpbeg; e < evpend; e++){
    if(evp[e].ty == CHC){ 
      cf = c + evp[e].ch; 
      if(ee == nev || ev[ee].tr != compiftra[c][cf]) emsg("event not agree");
      c = cf; ee++;
      //getnotdep(compleavesimnotdep[c]);
    }

    tnext = evp[e+1].t;
    if(t < tnext){
      getnotdep(compleavesimnotdep[c]);
      do{
        if(compleavesimdep[c].size() > 0){
          tt = t+partdt; if(tt > tnext) tt = tnext;
          tmid = 0.5*(t+tt); while(tnextold <= tmid){ co = cexist[k]; k++; if(k < nexist) tnextold = texist[k]; else tnextold = tmax;}
          getdep(compleavesimdep[c],co,tmid);
        }
        else{ Rdep = 0; tt = tnext;}

        R = Rdep + Rnotdep; 
        //if(nmfl == 1){ if(nfuteref > 0){ futnext = futev[futeref[nfuteref-1]].t; if(tt > futnext) tt = futnext;} else futnext = large;} else futnext = large;

        tst = t;
        if(ee < nev && ev[ee].t < tt){ // Actual event occurs
          tr = ev[ee].tr; t = ev[ee].t; ee++;
          switch(tra[tr].type){
            case EXP_TR:
              if(transdep[tra[tr].eq] == 1){
                j = 0; while(j < nRdep && compleavesimdep[c][j] != tr) j++; if(j == nRdep) emsg("rev prob");
                r = rdep[j];
              }
              else{
                j = 0; while(j < nRnotdep && compleavesimnotdep[c][j] != tr) j++; if(j == nRnotdep) emsg("rev prob");
                r = rnotdep[j];
              }
              prob += log(r);//cout << log(r) << " compr\n";
              c = tra[tr].cf;
              break;

              case END_TR: return prob;

            case GAMMAEND_TR: case WEIEND_TR: case FIXEDEND_TR:
              break;
              default: c = tra[tr].cf; break;
          }
          //if(nmfl == 1){ for(n = 0; n < nnonexptra[tr]; n++) addfute(nonexptra[tr][n],t);}    // adds on non-Markovian future events

          getnotdep(compleavesimnotdep[c]);
        }
        else t = tt;
        prob -= (t-tst)*R; 

  
        /*
        if(t == futnext){   // Future event
        n = futeref[nfuteref-1]; futeref.pop_back(); nfuteref--;
        tr = futev[n].tr;
        if(inside(c,n,tbeg,pa.pback) == 0){  // Adds the future event
        c += tra[tr].dc;
        EV ev; ev.tr = tr; ev.t = t; ev.tbeg = futev[n].tnow; evnew.push_back(ev);
        for(n = 0; n < nnonexptra[tr]; n++) addfute(nonexptra[tr][n],t);
      }
        else{
        EV ev; ev.tr = tra[tr].traother; ev.t = t; ev.tbeg = futev[n].tnow; evnew.push_back(ev);
      }
      }
        */
      }while(t < tnext);
    }
  }

  /*
  if(evpend == nevp-1){  // adds in any future events
  for(j = 0; j < nfuteref; j++){
  n = futeref[nfuteref-1]; tr = futev[n].tr;
  EV ev; ev.tr = tra[tr].traother; ev.t = t; ev.tbeg = futev[n].tnow; evnew.push_back(ev);
}
}

  pa.cend = c; pa.ev = evnew;

  if(nmfl == 1 && evpend < nevp-1) pa.futerefst = futeref;
  //tmp = 0;
  */
  return prob;
}

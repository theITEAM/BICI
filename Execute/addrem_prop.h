
void Chain::addrem_prop()                                                   // Adds and remove indiduals from the population
{
  long a, ag, c;
  long i, tr;
  double probif, probfi, z, tbirth, Ltoti, Ltotf, al, tleast, dLpri, pr;

  partcl = clall;

  //cout << "add rem\n";
  if(ran() < 0.5){     // adds an individual
    //cout << "add\n";
    ntr_add++;

    i = nindtot;

    probfi = log(1.0/(nindtot-nind+1));
    probif = 0;

    // First selects a time
    z = ran()*addremdivsumtot; a = 0; while(a <= naddremdiv && z > addremdivsum[a]) a++; if(a > naddremdiv) emsg("addremdic");
    probif += logaddremdiv[a];
 
    if(a == naddremdiv) tent = 0;
    else{ tent = (a+ran())*tmax/naddremdiv; probif += log(naddremdiv/tmax);}

    // Then selects the age
    z = ran()*induoagesumtot[a]; ag = 0; while(ag < nclassval[agecl] && z > induoagesum[a][ag]) ag++; if(ag == nclassval[agecl]) emsg("pro");
    probif += loginduoage[a][ag];

    if(tbirthfl == 1){
      if(ag == 0){ tbirth = tent - ran()*age[ag]; probif += log(1.0/age[ag]);}
      else{ tbirth = tent - (age[ag-1] + ran()*(age[ag]-age[ag-1])); probif += log(1.0/(age[ag]-age[ag-1]));}
      dLpri = agprior[ag];
    }
    else{ tbirth = tent; dLpri = 0;}

    if(sinkfl == 1){
      tlea = tent - log(ran())*induolifetime[a]; tleast = tlea;
      if(tlea > tmax){ tlea = tmax; probif += -(tlea-tent)/induolifetime[a];}
      else probif += -log(induolifetime[a]) - (tlea-tent)/induolifetime[a];
    }
    else tlea = tmax;

    addemptyind(tbirth);

    pr = indsimulate(i); if(pr == -large){ nfa_add++; rememptyind(i); return;}
    if(sourcefixfl == 1){
      tr = evnew[0].tr; if(tra[tr].capev >= 0 && capevall[tra[tr].capev] == 1){ nfa_add++; rememptyind(i); return;}
    }
    if(sinkfixfl == 1){
      tr = evnew[evnew.size()-1].tr; if(tra[tr].capev >= 0 && capevall[tra[tr].capev] == 1){ nfa_add++; rememptyind(i); return;}
    }

    probif += pr;

    if(checkevdt() == 1){ rememptyind(i); return;}

    Ltoti = Lwnm(); if(tempp == 1) Lout();
    indcha(i);
    Ltotf = Lwnm(); if(tempp == 1) Lout();

    al = exp(Ltotf - Ltoti + probfi - probif + dLpri); if(tempp == 1) cout << tempp << " " << al << " " << Ltotf << " "<< Ltoti << " " << probfi << " " << probif <<  " al add\n";
    //cout << al <<  "al\n";
    if(ran() < al){
      nac_add++;

      Lpri += dLpri;
      if(samp < burnin){    //  tunes proposals
        if(tent == 0) a = naddremdiv; else a = long(naddremdiv*tent/tmax);
        addremdiv[a]++;

        induoage[a][ag]++;

        c = tra[indev[i][0].tr].cf%ncomps; if(c == NOTALIVE) emsg("pr");
        if(indev[i][0].t == 0) simuo_initc[a][c]++; else simuo_enterc[a][c]++;

        if(sinkfl == 1) induolifetime[a] = 0.98*induolifetime[a] + (tleast - tent)*0.02;
      }
    }
    else{ indrev(i); rememptyind(i);}
  }
  else{     // removes an individual
    ntr_rem++;
    if(nindtot == nind) return;

    probif = log(1.0/(nindtot-nind));
    probfi = 0;

    i = nind+long(ran()*(nindtot-nind));

    getlifespan(i);
    if(tent == 0) a = naddremdiv; else{ a = long(naddremdiv*tent/tmax); probfi += log(naddremdiv/tmax);}

    probfi += logaddremdiv[a];

    c = tra[indev[i][0].tr].cf;

    ag = compval[c][agecl];
    probfi += loginduoage[a][ag];

    if(tbirthfl == 1){
      if(ag == 0) probfi += log(1.0/age[ag]);
      else probfi += log(1.0/(age[ag]-age[ag-1]));
      dLpri = -agprior[ag];
    }
    else dLpri = 0;

    if(sinkfl == 1){
      if(tlea == tmax) probfi += -(tlea-tent)/induolifetime[a];
      else probfi += -log(induolifetime[a]) - (tlea-tent)/induolifetime[a];
    }

    probfi += probindsimulate(i);

    evnew.clear();
    EV evbeg; evbeg.t = 0; evbeg.tr = tranull; evnew.push_back(evbeg);
    EV evend; evend.t = tmax; evend.tr = tranull; evnew.push_back(evend);

    Ltoti = Lwnm();  //Lout();
    indcha(i);
    Ltotf = Lwnm(); //Lout();

    al = exp(Ltotf - Ltoti + probfi - probif + dLpri); //cout << samp << " " << al << " " << Ltotf << " "<< Ltoti << " " << probfi << " " << probif <<  " al rem\n";

    if(ran() < al){ Lpri += dLpri; rememptyind(i); nac_rem++;}
    else indrev(i);
  }
}

void Chain::addreminit()                                     // Initialising proposals for adding and removing individuals
{
  long a, ag;

  induolifetime.resize(naddremdiv+1);
  addremdiv.resize(naddremdiv+1); addremdivsum.resize(naddremdiv+1); logaddremdiv.resize(naddremdiv+1);
  induoagesumtot.resize(naddremdiv+1); induoage.resize(naddremdiv+1); induoagesum.resize(naddremdiv+1); loginduoage.resize(naddremdiv+1);

  for(a = 0; a <= naddremdiv; a++){
    induolifetime[a] = tmax/10;

    addremdiv[a] = 50.0/naddremdiv;
    //if(a == naddremdiv) addremdiv[a] = 0;

    induoage[a].resize(nclassval[agecl]); induoagesum[a].resize(nclassval[agecl]); loginduoage[a].resize(nclassval[agecl]);

    for(ag = 0; ag < nclassval[agecl]; ag++) induoage[a][ag] = 50.0/nclassval[agecl];
  }
}

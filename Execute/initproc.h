// Initialialisation procedure for parameters 
// this disn't work because summary statistics not sufficiently good to descriminate parameters 


void Chain::initproc()
{
  long c, ce, i, j, f, faci, loop, loopmax = 4, num, n;
  double fac, Li, Lf, Lfmax = -large;
  vector<double> list;

  alinit.resize(ncompswa); for(c = 0; c < ncompswa; c++) alinit[c] = 1;
  nindinit = nind;

  for(ce = 0; ce < ncapev; ce++){
    list.clear();
    list.push_back(capevtmin[ce]);
    switch(capevtype[ce]){
      case CETRANS:
        for(f = 0; f < nfixev; f++){ if(fixev[f].capev == ce) list.push_back(fixev[f].t);}
        break;
    }
    list.push_back(capevtmax[ce]);
    sort(list.begin(),list.end());
    for(i = 0; i < list.size(); i++) cout << list[i] << ", "; cout<< "or\n";
    num = list.size();
    cout << num << " num\n";

    capevdivt.push_back(vector <double> ());
    capevdivnum.push_back(vector <long> ());

    for(loop = 0; loop <= loopmax; loop++){
      fac = (double(loop)/loopmax)*(num-1);
      faci = long(fac); fac -= faci; cout << faci << " " << fac << " ff\n";
      if(faci == 0) capevdivt[ce].push_back(list[faci]);
      else capevdivt[ce].push_back(list[faci]*(1-fac) + list[faci+1]*fac);
    }
    ncapevdiv.push_back(capevdivt[ce].size()-1);

    for(i = 0; i < ncapevdiv[ce]; i++){
      n = 0; for(j = 1; j < num-1; j++){ if(list[j] >= capevdivt[ce][i] && list[j] < capevdivt[ce][i+1]) n++;}
      capevdivnum[ce].push_back(n);
    }

    for(i = 0; i < ncapevdiv[ce]; i++){
      cout <<  ce << " " <<  capevdivt[ce][i] << " -> " << capevdivt[ce][1] << " " << capevdivnum[ce][i] << "  ran\n";
    }

    //cout << capevtmin[ce] << " " <<  capevtmax[ce] << " " << capevtype[ce] << "  tm\n";
    //for(c = 0; c < ncomp; c++) cout << capevfilt[ce][c] << ", "; cout << " filt\n";
  }

  initparamsamp();

  simon = 2;

  sim(0); Li = fit();

  traceinit();
  for(samp = 0; samp < 100; samp++){
    sim(0); Lf = fit();
    if(Lf > Lfmax) Lfmax = Lf;
    cout << Lf << " " << Lfmax << " Lf\n";
  }

 // for(i = 0; i < nindtot_sim; i++) cout << i << " " <<  nindev_sim[i] << " nev\n";
  simon = 0;
}

double Chain::fit()    // Estimates how well the simulated output agrees with the data
{
  long ce, c, i, e, j, ndiv;
  double tmi, tma, t, L, f;
  vector<long> num;

  L = 0;
  for(ce = 0; ce < ncapev; ce++){
    tmi = capevtmin[ce]; tma = capevtmax[ce]; ndiv = ncapevdiv[ce];

    num.resize(ndiv); for(j = 0; j < ndiv; j++) num[j] = 0;

    switch(capevtype[ce]){
      case CETRANS:
        j = 0;
        for(i = 0; i < nindtot_sim; i++){
          for(e = 1; e < nindev_sim[i]; e++){
            t = indev_sim[i][e].t;
            if(t > tmi){
              if(t > tma) break;
              while(t > capevdivt[ce][j+1]) j++;
              num[j]++;
            }
          }
        }
        break;
    }

    f = 0; for(j = 0; j < ndiv; j++) f += double(capevdivnum[ce][j]+num[j])/(2*ndiv);

    for(j = 0; j < ndiv; j++) L += -(capevdivnum[ce][j]-num[j])*(capevdivnum[ce][j]-num[j])/(f*f);

    cout << ce << ": "; 
    for(j = 0; j < ndiv; j++) cout << capevdivnum[ce][j] << "," << num[j] << "   ";
    cout << "num\n";

    /*
      cout << capevtmin[ce] << " " <<  capevtmax[ce] << " " << capevtype[ce] << "  tm\n";
      for(c = 0; c < ncomp; c++) cout << capevfilt[ce][c] << ", "; cout << " filt\n";
    }
    */
  }

  return L;
}


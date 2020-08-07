/*
const double shmin = 1, shmax = 100;
const short nsh = 500;
const double gmax = 20;
const short ng = 500;
*/

const double shmin = 1, shmax = 20;
const short nsh = 400;
const double gmax = 10;
const short ng = 400;

double tgammaprob(double x, double a, double b);

//vector <vector <double> > gammast, gammaupst;
vector <vector <double> > gammaupst;

/*
double Chain::likegammatot()
{
  short i;
  double Ligam = 0;
  for(i = 0; i < nindtot; i++) Ligam += likegamma(i);
  return Ligam;
}


void Chain::setgammaeq()
{
  short k;
  for(k = 0; k < ngammaeq; k++){
    gammaeqmean_val[k] = calculatenotdep(gammaeqmean[k],param);
    gammaeqshape_val[k] = calculatenotdep(gammaeqshape[k],param);
    if(gammaeqshape_val[k] > shmax) gammaeqshape_val[k] = shmax-tiny;
  }
}


double Chain::likegamma(short i)   // works out the likelihood of all the gamma shaped transitions
{
  short cl, cl2, clnext, trnext, stat[nclass], index[nclass], geq, ot, tr, ntra, tra, e;
  double Ligam, tenter, tleave, tnext, mean, kshape;

  Ligam = 0;
  nev = nindev[i][agecl];
  if(nev == 0){
    if(indinit[i][agecl] == -1) return 0;
    tenter = 0; 
    tleave = tmax;
  }
  else{
    if(transi[agecl][indev[i][agecl][0]] == -1) tenter = indevt[i][agecl][0]; else tenter = 0;
    if(transf[agecl][indev[i][agecl][nev-1]] == -1) tleave = indevt[i][agecl][nev-1]; else tleave = tmax;
  }

  for(cl = 0; cl < nclass; cl++){ stat[cl] = indinit[i][cl]; if(cl < nclass-2) index[cl] = -1; else index[cl] = 0;}

  do{
    tnext = tleave;
    for(cl = 0; cl < nclass; cl++){
      if(index[cl] == -1){ clnext = cl; tnext = tenter; trnext = -1; break;}
      if(index[cl] < nindev[i][cl] && indevt[i][cl][index[cl]] < tnext){
        tnext = indevt[i][cl][index[cl]]; clnext = cl; trnext = indev[i][cl][index[cl]];;
      }
    }
    if(tnext == tleave) break;

    cl = clnext;
    if(trnext >= 0){ stat[cl] = transf[cl][trnext];}
    index[cl]++;

    ntra = ntransgamma[cl][stat[cl]];
    for(tra = 0; tra < ntra; tra++){
      tr = transgamma[cl][stat[cl]][tra];

      ot = 0; for(cl2 = 0; cl2 < nclass; cl2++) ot += stat[cl2]*classtransmult[cl][cl2];
      geq = transeq[cl][tr][ot];
      if(geq >= 0){
        mean = gammaeqmean_val[geq];
        kshape = gammaeqshape_val[geq];
        e = index[cl];
        if(e < nindev[i][cl]){
          if(indev[i][cl][e] == tr) Ligam += gammaprob(indevt[i][cl][e] - tnext,kshape,kshape/mean);
          else Ligam += gammaup(indevt[i][cl][e] - tnext,kshape,kshape/mean);
        }
        else Ligam += gammaup(tleave - tnext,kshape,kshape/mean);
      }
    }
  }while(1 == 1);

  return Ligam;
}

double Chain::likegammasingle(short i, short cl) // works out the likelihood for a single classification
{
  short nev, j, trnext, geq, tra, ntra, eot, ot, e, tr;
  double Ligam, t, tnext, mean, kshape, tenter, tleave;

  nev = nindev[i][agecl];
  if(nev == 0){ tenter = 0; tleave = tmax;}
  else{
    if(transi[cl][indev[i][agecl][0]] == -1) tenter = indevt[i][agecl][0]; else tenter = 0;
    if(transf[cl][indev[i][agecl][nev-1]] == -1) tleave = indevt[i][agecl][nev-1]; else tleave = tmax;
  }

  j = indinit[i][cl];
  nev = nindev[i][cl];
  t = tenter;

  Ligam = 0;

  eot = 0; ot = evot[eot];
  for(e = 0; e <= nev; e++){
    if(e < nev){ tnext = indevt[i][cl][e];  trnext = indev[i][cl][e];} else{ tnext = tleave; trnext = -1;}

    while(eot+1 < nevot && evott[eot+1] < t){ eot++; ot = evot[eot];}

    ntra = ntransgamma[cl][j];
    for(tra = 0; tra < ntra; tra++){
      tr = transgamma[cl][j][tra];

      geq = transeq[cl][tr][ot];
      if(geq >= 0){
        mean = gammaeqmean_val[geq];
        kshape = gammaeqshape_val[geq];
        if(tr == trnext) Ligam += gammaprob(tnext - t,kshape,kshape/mean);
        else Ligam += gammaup(tnext - t,kshape,kshape/mean);
      }
    }
    j = transf[cl][trnext]; t = tnext;
  }

  return Ligam;
}
*/

void gammainit()
{
  short shi, i, lo, nlo = 400;
  double shape, x, dx, sum, beta, sumst[ng];

  gammaupst.resize(nsh); for(shi = 0; shi < nsh; shi++) gammaupst[shi].resize(ng);

  if(1 == 0){
    dx = gmax/ng;

    for(shi = 0; shi < nsh; shi++){
      shape = shmin + shi*(shmax-shmin)/nsh;
      beta = shape;

      sum = 0;
      for(i = 0; i < ng; i++){
        sumst[i] = sum;
        for(lo = 0; lo < nlo; lo++){ x = (i+double(lo)/nlo)*dx; sum += tgammaprob(x,shape,beta)*dx/nlo;}
      }
  
      for(i = 0; i < ng; i++){
        //gammast[shi][i] = tgammaprob(i*dx,shape,beta);
        gammaupst[shi][i] = log(1.0000000001-sumst[i]/sum);
      }
      cout << shi << "  "<< sum << "sum\n";
    }
    ofstream gam("GP/gamtab3");
    for(shi = 0; shi < nsh; shi++){
      for(i = 0; i < ng; i++){
        //gam <<  gammast[shi][i] << " " << gammaupst[shi][i] << "\n";
        gam << gammaupst[shi][i] << "\n";
      }
    }
  }
  else{
    //ifstream gam("GP/gamtab");
    ifstream gam("GP/gamtab3");
    for(shi = 0; shi < nsh; shi++){
      for(i = 0; i < ng; i++){
        //gam >> gammast[shi][i] >> gammaupst[shi][i];
        gam >> gammaupst[shi][i];
      }
    }
  }
}

double gammaup(double x, double a, double b)     // the integral of the gamma function from x up to infinity
{
  short shi, i;
  double sh, ifl, dsh, di;

  if(x < 0 || a < 0 ||  b < 0) emsg("gam neg prob");

  sh = (a - shmin)*nsh/(shmax-shmin);
  shi = short(sh); if(shi < 0 || shi >= nsh) emsg("gamma out of range");

  ifl = (ng/gmax)*x/(a/b);
  i = short(ifl); if(i < 0) emsg("gamma out of range");

  if(i >= ng) return gammaupst[shi][ng-1] - (i-ng);

  dsh = sh-shi;
  di = ifl-i;
  // cout << shi << " "<< i << " " << dsh << " "<< di << " ii\n";
  return gammaupst[shi][i]*(1-dsh)*(1-di) + gammaupst[shi+1][i]*dsh*(1-di) +gammaupst[shi][i+1]*(1-dsh)*di +gammaupst[shi+1][i+1]*dsh*di;
}

double gammasamp(double a, double b)   // draws a sample from    x^(a-1)*exp(-b*x)
{
  if(a < 1){
    double u = ran();
    return gammasamp(1.0 + a, b) * pow (u, 1.0 / a);
  }
  else{
    double x, v, u;
    double d = a - 1.0 / 3.0;
    double c = (1.0 / 3.0) / sqrt (d);
 
    while(1 == 1){
      do{
        x = sqrt(-2*log(ran()))*cos(2*3.141592654*ran());
        v = 1.0 + c * x;
      }while (v < 0);

      v = v*v*v;
      u = ran();

     // if (u < 1 - 0.0331*x*x*x*x) break;

      if (log(u) < 0.5*x*x + d*(1 - v + log(v))) break;
    }

    return d*v/b;
  }
}


double gammaprob(double x, double a, double b)                                      // The log of the probability from the gamma distribution
{
  if(x < 0 || a < 0 || b < 0) emsg("neg gam prob");
  return (a-1)*log(x) - b*x + a*log(b) - lgamma(a);
}

double tgammaprob(double x, double a, double b)
{
  return pow(x,(a-1))*exp(- b*x)*pow(a,b)/tgamma(a);
}

/*
double betaprob(double x, double a, double b)
{
  return (a-1)*log(x) + (b-1)*log(1-x) + qgammafunc(a+b) - qgammafunc(a) - qgammafunc(b);
}
*/

double weibullsamp(double lam, double kk)                                           // Samples from the Weibull distribution
{
  return -lam*pow(log(1-ran()),1.0/kk);
}

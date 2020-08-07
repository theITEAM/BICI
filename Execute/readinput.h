const long EXPONENTIAL = 0, GAM = 1, WEI = 4, SOURCE = 2, SINK = 3; 

void sortparam();
void sortcap();
void sortobs();
void sortpopm();
void sortderm();

//void sortfixev();

void createtra(long type, XMLNode* child3, long cl, long i, long f);
void findparam(string par, string s);      // Find the paramters with a particular name and attribute
vector<string> allotherlist;
vector<long> plist;

vector<string> getallpos(string name);

vector<long> state;
vector<string> staterate;

void readinput(string file)
{
  long fi, cl, i, f, fl, tr, j, eq, c, cc, ci, cf, fmi, fma, jj, k, p, pp, leng, ob, cap, trai, traf, d, dm, ty, ce, a;
  long s, nsmoothst;
  string s1, s2, s3, s4, name, id, st, na, type, capna, capevna, par;
  double t, val, sd, ti, tmin, dttiny;
  vector<string> oprob;
  vector< vector<string> > oprobvec;
  vector<long> emp;
  vector<long> list;

  XMLNode* child;
  XMLNode* child2;
  XMLNode* child3;
  XMLNode* child4;
  XMLDocument doc;

  ifstream inFile(file.c_str());

  stringstream strStream;
  strStream << inFile.rdbuf();
  string str = strStream.str();
  if(str.size() == 0) emsg2("No input file");

  char *xml=new char[long(str.size())+1];
  xml[str.size()]=0;
  memcpy(xml,str.c_str(),str.size());
  doc.Parse(xml);

  XMLElement* root = doc.FirstChildElement();

  for(child = root->FirstChild(); child; child = child->NextSibling()){                      // Reads in the classification names
    s1 = child->Value();
    if(s1.compare("model") == 0){
      for(child2 = child->FirstChild(); child2; child2 = child2->NextSibling()){
        s2 = child2->Value();
        if(s2 == "classification"){
          name = get(child2,"name");
          classname.push_back(name);
          classval.push_back(getcommasep(get(child2,"value")));
        }
      }
    }
  }

  ncapevtrange = 0;                                 // divides the timeline into regions with different capture event rules
  tmin = large; tmax = -large;
  for(child = root->FirstChild(); child; child = child->NextSibling()){
    s1 = child->Value();
    if(s1.compare("data") == 0){  // model
      for(child2 = child->FirstChild(); child2; child2 = child2->NextSibling()){
        s2 = child2->Value();
        if(s2.compare("unobserved")==0) addremfl = 1;

        if(s2.compare("capev") == 0 || s2.compare("capev") == 0 || s2.compare("capev") == 0){
          t = getnum(child2,"tmin"); if(t < tmin) tmin  = t; if(t > tmax) tmax = t;
          fi = 0; while(fi < ncapevtrange && t > capevtrange[fi]) fi++;
          if(!(fi < ncapevtrange && t == capevtrange[fi])){ capevtrange.insert(capevtrange.begin()+fi,t); ncapevtrange++;}

          t = getnum(child2,"tmax"); if(t < tmin) tmin  = t; if(t > tmax) tmax = t;
          fi = 0; while(fi < ncapevtrange && t > capevtrange[fi]) fi++;
          if(!(fi < ncapevtrange && t == capevtrange[fi])){ capevtrange.insert(capevtrange.begin()+fi,t); ncapevtrange++;}
        }

        if(s2.compare("capture") == 0 || s2.compare("population") == 0 || s2.compare("derived") == 0 || s2.compare("observation") == 0 || 
           s2.compare("transition") == 0 || s2.compare("move") == 0){
          t = getnum(child2,"time"); if(t < tmin) tmin  = t; if(t > tmax) tmax = t;
        }
      }
    }

    if(s1.compare("simulation") == 0) simon = 1;
    if(s1.compare("simulation") == 0 || s1.compare("inference") == 0){ 
      const char *a = child->ToElement()->Attribute("nsamp");
      if(a) nsamp = atoi(a);
      const char *a2 = child->ToElement()->Attribute("tmin");
      if(a2){
        tminactual = atof(a2);
        if(tminactual > tmin) emsg("inference time must be before data time h");
        if(tminactual < tmin) tmin = tminactual;
      }
      const char *a3 = child->ToElement()->Attribute("tmax");
      if(a3){
        if(simon == 1) tmax = tmin;
        tmaxactual = double(atof(a3));
        if(tmaxactual < tmax) emsg("inference time must be before data time");
        tmax2 = tmaxactual;//+dttiny;
      }

      indmax = long(getnum(child,"indmax"));
    }
  }

  if(tmin == large || tmax == -large) emsg("No data!");

  if(ncapevtrange > 0){ if(tmin == capevtrange[0]){ capevtrange.erase(capevtrange.begin()); ncapevtrange--;}}
  if(ncapevtrange > 0){ if(tmax == capevtrange[ncapevtrange-1]){ capevtrange.pop_back(); ncapevtrange--;}}
  capevtrange.push_back(tmax);

  if(tmax == tmin){ notimerange = 1; tmax += tiny*10;}
  dttiny = (tmax-tmin)*0.0001;
  tmin -= dttiny; tmax += dttiny;

  if(tmax2 < tmax) tmax2 = tmax;

  classname.push_back("Fix");
  vector<string> vals; for(fi = 0; fi <= ncapevtrange; fi++){ stringstream ss; ss << "Fix" << fi; vals.push_back(ss.str());}
  classval.push_back(vals);

  nclass = classname.size();
  for(cl = 0; cl < nclass; cl++) nclassval.push_back(classval[cl].size());

  agecl = nclass-3;
  settimecl = nclass-2;
  capevcl = nclass-1;

  clall = nclass-3;

  if(nclassval[agecl] > 1) tbirthfl = 1;

  ncomp = 1;
  for(cl = 0; cl < nclass; cl++){
    classmult.push_back(ncomp);
    ncomp *= nclassval[cl];
  }
  ncomps = classmult[agecl];
  ncompswa = classmult[settimecl];
  ncompswf = classmult[capevcl];

  NOTALIVE = ncomp;

  ntra = 0;
  //ngammaeq = 0;
  //age.push_back(0);
  for(child = root->FirstChild(); child; child = child->NextSibling()){
    s1 = child->Value();
    if(s1.compare("model") == 0){  // model
      for(child2 = child->FirstChild(); child2; child2 = child2->NextSibling()){
        s2 = child2->Value();
        if(s2 == "classification"){
          name = get(child2,"name");
          for(cl = 0; cl < nclass; cl++) if(classname[cl] == name) break;
          if(cl == nclass) emsg("classification wrong");

          for(child3 = child2->FirstChild(); child3; child3 = child3->NextSibling()){
            s3 = child3->Value();
            if(s3 == "transition"){
              i = getclassval(cl,get(child3,"from"));
              f =  getclassval(cl,get(child3,"to"));

              fl = 0;
              name = get(child3,"type");
              if(name == "Exponential"){ createtra(EXPONENTIAL,child3,cl,i,f); fl = 1;}

              if(name == "Gamma"){ createtra(GAM,child3,cl,i,f); fl = 1; gammafl=1; nmfl = 1;}

              if(name == "Weibull"){ createtra(WEI,child3,cl,i,f); fl = 1; nmfl = 1;}

              if(name == "Fixed"){
                getalltrans(child3,child3,cl,i,get(child3,"time"));
                for(j = 0; j < state.size(); j++){
                  c = state[j];
                  TT tt; tt.type = FIXED_TR; tt.ci = c; tt.cf = c+(f-i)*classmult[cl];
                  tt.cl = cl; tt.i = i; tt.f = f; tt.capev = -1; tt.like = 1; tt.eq = addequation(staterate[j]);
                  tra.push_back(tt);
                  ntra++;
                }
                fl = 1;
              }

              if(name == "Grow"){
                if(cl != agecl) emsg("grow prob");
                age.push_back(getnum(child3,"age"));

                getalltrans(child3,child3,cl,i,get(child3,"age"));
                for(j = 0; j < state.size(); j++){
                  c = state[j];
                  TT tt; tt.type = GROW_TR; tt.ci = c; tt.cf = c+(f-i)*classmult[cl]; tt.cl = cl; 
                  tt.i = i; tt.f = f;  tt.eq = -1; tt.capev = -1; tt.like = 1;
                  tra.push_back(tt);
                  ntra++;
                }
                fl = 1;
              }

              if(name == "Set"){
                if(cl != settimecl) emsg("grow prob");
                settime.push_back(getnum(child3,"time"));

                getalltrans(child3,child3,cl,i,get(child3,"time"));
                for(j = 0; j < state.size(); j++){
                  c = state[j];
                  TT tt; tt.type = SETTIME_TR; tt.ci = c; tt.cf = c+(f-i)*classmult[cl]; tt.cl = cl; 
                  tt.i = i; tt.f = f; tt.eq = -1; tt.capev = -1; tt.like = 1;
                  tra.push_back(tt);
                  ntra++;
                }
                fl = 1;
              }

              if(fl == 0) emsg("type not known");
            }
          }
        }
      }
    }
  }

  for(child = root->FirstChild(); child; child = child->NextSibling()){             // Source and sink information
    s1 = child->Value();
    if(s1.compare("model") == 0){
      for(child2 = child->FirstChild(); child2; child2 = child2->NextSibling()){
        s2 = child2->Value();
        if(s2 == "source"){ createtra(SOURCE,child2,-2,-1,-1); sourcefl = 1;}
        if(s2 == "sink"){ createtra(SINK,child2,-2,-1,-1); sinkfl = 1;}
      }
    }
  }

  nderive = 0;
  for(child = root->FirstChild(); child; child = child->NextSibling()){
    s1 = child->Value();
    if(s1.compare("model") == 0){  // model
      for(child2 = child->FirstChild(); child2; child2 = child2->NextSibling()){
        s2 = child2->Value();

        if(s2.compare("derived") == 0){
          derivename.push_back(get(child2,"name"));
          eq = addequation(get(child2,"value")); if(eq == -1) emsg("prob derive");
          derive.push_back(eq);
          nderive++;
        }
      }
    }
  }
  //long d; for(d = 0; d < nderive; d++) cout << derivename[d] << " " << derive[d] << "uu\n";

  for(f = 0; f < ncapevtrange; f++){  // adds in the capevtrange transitions
    for(c = 0; c < ncompswf; c++){
      TT tt; tt.type = CAPEVTRANGE_TR; tt.ci = c+f*classmult[capevcl]; tt.cf =  c+(f+1)*classmult[capevcl]; 
      tt.cl = capevcl; tt.i = f; tt.f = f+1;  tt.eq = -1; tt.capev = -1;
      tra.push_back(tt);
      ntra++;
    }
  }

  nage = age.size(); nsettime = settime.size();
  if(nage == 0) age.push_back(1); else age.push_back(age[nage-1]*(1+1.0/nage)); // arbitratily sets an upper bound on the age

  compiftra.resize(ncomp+1); for(c = 0; c <= ncomp; c++){ compiftra[c].resize(ncomp+1); for(cc = 0; cc <=  ncomp; cc++) compiftra[c][cc] = -1;}
  for(tr = 0; tr < ntra; tr++) compiftra[tra[tr].ci][tra[tr].cf] = tr;

  compval.resize(ncomp); for(c = 0; c < ncomp; c++){ compval[c].resize(nclass); for(cl = 0; cl < nclass; cl++) compval[c][cl] = (c/classmult[cl])%nclassval[cl];}

  compname.resize(ncomp+1);
  for(c = 0; c < ncomp; c++){
    stringstream ss;
    for(cl = 0; cl < nclass; cl++){ if(nclassval[cl] > 1){ if(cl > 0) ss << "_"; ss << classval[cl][compval[c][cl]];}}
    compname[c] = ss.str();
  }
  compname[ncomp] = "Not alive";

  /*
  for(p = 0; p < popnumname.size(); p++){
    cout << popnumname[p] << "new\n";
    for(j = 0; j < popnumterm[p].size(); j++){
      cout << compname[popnumterm[p][j]] << " " << popnumtermweight[p][j] << "  ";
    }
    cout << " vec\n\n";
 
  }
  emsg("P");
  */

  for(c = 0; c < ncomp; c++) emp.push_back(-1);

   // reads in capture information
  ncap = 0;
  for(child = root->FirstChild(); child; child = child->NextSibling()){
    s1 = child->Value();
    if(s1.compare("data") == 0){  // model
      for(child2 = child->FirstChild(); child2; child2 = child2->NextSibling()){
        s2 = child2->Value();
        if(s2.compare("capture") == 0){
          capfl = 1;
          capname.push_back(get(child2,"name"));
          capt.push_back(getnum(child2,"time"));
          capprobeqn.push_back(emp);
          for(child3 = child2->FirstChild(); child3; child3 = child3->NextSibling()){
            s3 = child3->Value();
            if(s3.compare("probability") == 0){
              getalltrans(child2,child3,-1,-1,get(child3,"value"));
              for(j = 0; j < state.size(); j++){
                eq = addequation(staterate[j]); if(eq == -1) emsg("Equation problem");
                capprobeqn[ncap][state[j]] = eq;
              }
            }
          }
          ncap++;
        }
      }
    }
  }
  sortcap();

  nobs = 0; nind = 0; npopm = 0; nderm = 0; nfixev = 0; ncapev = 0;
  oprob.resize(nclass); oprobvec.resize(nclass);
  for(child = root->FirstChild(); child; child = child->NextSibling()){
    s1 = child->Value();
    if(s1.compare("data") == 0){  // model
      for(child2 = child->FirstChild(); child2; child2 = child2->NextSibling()){
        s2 = child2->Value();
        if(s2.compare("population") == 0){      // a population measurement
          popfl = 1;
          popmt.push_back(getnum(child2,"time"));
          val = double(getint(child2,"value")); if(val == 0) val = 0.001;
          sd = getnum(child2,"sd"); if(sd == 0) sd = 0.001;
          popmal.push_back(val*val/(sd*sd));
          popmbe.push_back(val/(sd*sd));
          getalltrans(child2,child2,-1,-1,"");
          popmcomp.push_back(vector <long> ());
          popmcomp[npopm].resize(ncomp+1); for(c = 0; c <= ncomp; c++) popmcomp[npopm][c] = 0;
          for(j = 0; j < state.size(); j++) popmcomp[npopm][state[j]] = 1;
          npopm++;
        }

        if(s2.compare("derived") == 0){      // a derived measurement
          derfl = 1;
          name = get(child2,"name");
          d = 0; while(d < nderive && name != derivename[d]) d++;
          if(d == nderive) emsg("Cannot find derive name");
          derm.push_back(d);
          dermt.push_back(getnum(child2,"time"));
          val = double(getint(child2,"value"));
          sd = getnum(child2,"sd"); if(sd == 0) sd = 0.001;
          dermav.push_back(val);
          dermvar.push_back(sd*sd);
          nderm++;
        }

        if(s2.compare("capev") == 0){
          capevprobeqn.push_back(addequation(get(child2,"pd")));
          if(get(child2,"pd") == "1") capevall.push_back(1);
          else capevall.push_back(0);
          capevfl = 1;

          if(get(child2,"type") == "source") ty = CESOURCE;
          else{
            if(get(child2,"type") == "sink") ty = CESINK;
            else{
              if(get(child2,"type") == "trans") ty = CETRANS;
              else emsg("type not set");
            }
          }

          t = getnum(child2,"tmin"); fmi = 0; while(fmi <= ncapevtrange && capevtrange[fmi] <= t) fmi++;
          capevtmin.push_back(t);
          t = getnum(child2,"tmax"); fma = 0; while(fma <= ncapevtrange && capevtrange[fma] <= t) fma++;
          capevtmax.push_back(t);

          capevfilt.push_back(vector<long> ());   // This filter is used in initproc
          capevfilt[ncapev].resize(ncomp); for(c = 0; c < ncomp; c++) capevfilt[ncapev][c] = 0;
          capevtype.push_back(ty);

          capevname.push_back(get(child2,"name"));

          switch(ty){
            case CESOURCE:  // if source events are observed then it removes this possibility tra
              getalltrans(child2,child2,-1,-1,"");
              ci = NOTALIVE;
              for(jj = 0; jj < state.size(); jj++){
                cf = state[jj];
                fi = compval[cf][capevcl];
                if(fi >= fmi && fi < fma){ tr = compiftra[ci][cf]; if(tr != -1) tra[tr].capev = ncapev;}
              }
              break;

            case CESINK:         // if sink events are observed then it removes this possibility tra
              getalltrans(child2,child2,-1,-1,"");
              cf = NOTALIVE;
              for(jj = 0; jj < state.size(); jj++){
                ci = state[jj];
                fi = compval[ci][capevcl];
                if(fi >= fmi && fi < fma){ tr = compiftra[ci][cf]; if(tr != -1) tra[tr].capev = ncapev;}
              }
              break;

            case CETRANS:       // if transition events are observed then it removes this possibility tra
              for(cl = 0; cl < nclass; cl++){
                i = getclassval(cl,get(child2,"from"));
                f = getclassval(cl,get(child2,"to"));
                if(i >= 0 && f >= 0) break;
              }
              if(cl == nclass) emsg("Cannot find fixed transition");

              getalltrans(child2,child2,cl,i,"");

              for(jj = 0; jj < state.size(); jj++){
                ci = state[jj]; cf = state[jj]+(f-i)*classmult[cl];
                fi = compval[ci][capevcl];
                if(fi >= fmi && fi < fma){ 
                  tr = compiftra[ci][cf]; if(tr == -1) emsg("compif3"); tra[tr].capev = ncapev;
                }
              }
              break;
          }

          for(jj = 0; jj < state.size(); jj++) capevfilt[ncapev][state[jj]] = 1;

          ncapev++;
        }

        if(s2.compare("observation") == 0 || s2.compare("transition") == 0 ||  s2.compare("move") == 0){
          id = get(child2,"id");
          ti = getnum(child2,"time");
          i = 0; while(i < indid.size() && id != indid[i]) i++;
          if(i == indid.size()){
            indid.push_back(id); indobs.push_back(vector<long>()); indobst.push_back(vector<double>());
            indfixtenter.push_back(large); indfixtleave.push_back(large);
            indlikeenter.push_back(1); indlikeleave.push_back(1);
            indfixentercapev.push_back(-1); indfixleavecapev.push_back(-1);
            nind++; 
            if(nind > indmax){
              stringstream ss; ss << "Maximum number of " << indmax << " individuals exceeded.";
              emsg(ss.str()); 
            }
          }
        }

        if(s2.compare("observation") == 0){
          obst.push_back(ti);
          obsi.push_back(i);

          for(cl = 0; cl < nclass; cl++){
            if(exist(child2,classname[cl].c_str()) == 1) oprob[cl] = get(child2,classname[cl].c_str());
            else{
              stringstream ss;
              if(nclassval[cl] == 1) ss << classval[cl][0];
              else{
                for(j = 0; j < nclassval[cl]; j++){ ss << classval[cl][j] << ":1"; if(j < nclassval[cl]-1) ss << "|";}
              }
              oprob[cl] = ss.str();
            }
            oprobvec[cl] = getlinesep(oprob[cl]);
          }

          obsprobeqn.push_back(emp);
          //cout << i << ":";
          for(c = 0; c < ncomp; c++){
            st = "1";
            for(cl = 0; cl < nclass; cl++){
              na = classval[cl][compval[c][cl]];
              if(oprobvec[cl].size() == 1){ if(oprob[cl] != na && oprob[cl] != na+":1") break;}
              else{
                for(k = 0; k < oprobvec[cl].size(); k++){
                  if(oprobvec[cl][k].substr(0,len(na)+1) == na+":"){
                    if(oprobvec[cl][k].substr(len(na)+1) == "0") st = "0";
                    else{
                      if(st == "1") st = oprobvec[cl][k].substr(len(na)+1);
                      else st += "*("+oprobvec[cl][k].substr(len(na)+1)+")";
                    }
                    break;
                  }
                }
                if(k == oprobvec[cl].size() || st == "0") break;
              }
            }
            if(cl == nclass){
              eq = addequation(st); if(eq == -1) emsg("Equation problem");
              //cout << st << ",";
              obsprobeqn[nobs][c] = eq;
            }
          }
          //cout << " eq\n";

          obscap.push_back(-1);
          if(exist(child2,"capture")){ 
            capna = get(child2,"capture");
            cap = 0; while(cap < ncap && capname[cap] != capna) cap++;
            if(cap < ncap) obscap[nobs] = cap; else{ emsg("Cannot find capture");}
          }
          nobs++;
        }

        if(s2.compare("transition") == 0 || s2.compare("move") == 0){
          if(s2.compare("transition") == 0){
            capevna = get(child2,"capev");
            ce = 0; while(ce < ncapev && capevname[ce] != capevna) ce++; if(ce == ncapev) emsg("capev");
            ty = capevtype[ce];
          }
          else{
            ce = -1;
            if(get(child2,"type") == "source") ty = CESOURCE;
            else{
              if(get(child2,"type") == "sink") ty = CESINK;
              else{
                if(get(child2,"type") == "trans") ty = CETRANS;
                else emsg("type not set");
              }
            }
          }

          switch(ty){
            case CESOURCE:   // reads in fixed source transitions
              sourcefixfl = 1;
              indfixtenter[i] = ti-ran()*tiny;
              indfixentercapev[i] = ce;
              if(s2.compare("move") == 0){ indlikeenter[i] = 0; movefl = 1;}
             break;

            case CESINK:     // reads in fixed sink transitions
              sinkfixfl = 1;
              indfixtleave[i] = ti+ran()*tiny;
              indfixleavecapev[i] = ce;
              if(s2.compare("move") == 0){ indlikeleave[i] = 0; movefl = 1;}
              break;

            case CETRANS:  // reads in fixed transitions / moves
              fixfl = 1;

              for(cl = 0; cl < nclass; cl++){
                trai = getclassval(cl,get(child2,"from"));
                traf = getclassval(cl,get(child2,"to"));
                if(trai >= 0 && traf >= 0) break;
              }
              if(cl == nclass) emsg("Cannot find fixed transition");

              for(tr = 0; tr < ntra; tr++){ if(tra[tr].cl == cl && tra[tr].i == trai && tra[tr].f == traf) break;}

              if(tr == ntra) emsg("cannot find");

              FEV fev; 
              fev.i = i; fev.t = ti+ran()*tiny; fev.cl = cl; fev.trai = trai; fev.traf = traf; fev.capev = ce;
              fev.dc = (traf-trai)*classmult[cl]; fev.type = tra[tr].type;
              if(s2.compare("transition") == 0) fev.like = 1; else{ fev.like = 0; movefl = 1;}
              fixev.push_back(fev); nfixev++;
              break;
          }
        }
      }
    }
  }

  //for(i = 0; i < nind; i++){ if(nindfixev[i] > 0) sort(indfixev[i].begin(),indfixev[i].end(),compareFEV);}
  sort(fixev.begin(),fixev.end(),compareFEV);
  sortpopm(); sortderm(); sortobs(); //sortfixev();

  dirchal.resize(ncompswa); for(c = 0; c < ncompswa; c++) dirchal[c] = 1;

  // reads in prior information

  for(child = root->FirstChild(); child; child = child->NextSibling()){
    s1 = child->Value();
    if(s1.compare("model") == 0){  // model
      for(child2 = child->FirstChild(); child2; child2 = child2->NextSibling()){
        s2 = child2->Value();
        if(s2.compare("prior") == 0  || s2.compare("set") == 0){
          if(s2.compare("prior") == 0 && get(child2,"type") == "Dirichlet"){
            val = getnum(child2,"val1");
            vector<string> vec;
            vector<long> map;
            vec = getcommasep(get(child2,"parameter").substr(3));
            map.resize(ncompswa); for(c = 0; c < ncompswa; c++) map[c] = 1;

            for(j = 0; j < vec.size(); j++){
              st = vec[j];
              for(cl = 0; cl < nclass-2; cl++) if(classname[cl] == st) break;
              if(cl == nclass-2){
                for(cl = 0; cl < nclass-2; cl++){
                  for(k = 0; k < nclassval[cl]; k++){ if(classval[cl][k] == st) break;}
                  if(k < nclassval[cl]){
                    for(c = 0; c < ncompswa; c++){ if(compval[c][cl] != k) map[c] = 0;}
                    break;
                  }
                }
                if(cl == nclass-2) emsg("problem with dirichlet");
              }
            }
            for(c = 0; c < ncompswa; c++){ if(map[c] == 1) dirchal[c] = val;}
          }
          else{
            vector<string> pos;
            pos = getallpos(get(child2,"parameter"));

            for(pp = 0; pp < pos.size(); pp++){
              p = getparam(pos[pp]);
              paramprior[p] = priortype.size();
              paramclass[p] = get(child2,"class");

              if(s2.compare("set") == 0){
                priortype.push_back(FLAT);
                priorparam.push_back(p);
                priorminval.push_back(getnum(child2,"value"));
                priormaxval.push_back(getnum(child2,"value"));

                prioreq1.push_back(0);
                prioreq2.push_back(0);
              }
              else{
                type = get(child2,"type");

                fl = 0;
                if(type == "flat" || type == "Flat" || type == "FLAT"){
                  priortype.push_back(FLAT);
                  priorparam.push_back(p);
                  priorminval.push_back(getnum(child2,"val1"));
                  priormaxval.push_back(getnum(child2,"val2"));
                  prioreq1.push_back(-1);
                  prioreq2.push_back(-1);
                  fl = 1;
                }

                if(type == "unbounded" || type == "Unbounded" || type == "UNBOUNDED"){
                  priortype.push_back(FLAT);
                  priorparam.push_back(p);
                  priorminval.push_back(-large);
                  priormaxval.push_back(large);
                  prioreq1.push_back(-1);
                  prioreq2.push_back(-1);
                  fl = 1;
                }

                if(type == "Fix"){
                  priortype.push_back(FIX);
                  priorparam.push_back(p);
                  priorminval.push_back(getnum(child2,"val1"));
                  priormaxval.push_back(getnum(child2,"val1"));
                  prioreq1.push_back(-1);
                  prioreq2.push_back(-1);
                  fl = 1;
                }
  
                if(type == "gamma" || type == "Gamma"){
                  priortype.push_back(GAMMA);
                  priorparam.push_back(p);
                  priorminval.push_back(0);
                  priormaxval.push_back(0);
                  eq = addequation(get(child2,"val1")); if(eq == -1) emsg("prior prob1");
                  prioreq1.push_back(eq);
                  eq = addequation(get(child2,"val2")); if(eq == -1) emsg("prior prob2");
                  prioreq2.push_back(eq);
                  fl = 1;
                }
  
                if(type == "Normal"){
                  priortype.push_back(NORMAL);
                  priorparam.push_back(p);
                  priorminval.push_back(0);
                  priormaxval.push_back(0);
                  eq = addequation(get(child2,"val1")); if(eq == -1) emsg("prior prob3");
                  prioreq1.push_back(eq);
                  eq = addequation(get(child2,"val2")); if(eq == -1) emsg("prior prob4");
                  prioreq2.push_back(eq);
                  fl = 1;
                }

                if(type == "Log-Normal"){
                  priortype.push_back(LOGNORMAL);
                  priorparam.push_back(p);
                  priorminval.push_back(0);
                  priormaxval.push_back(0);
                  eq = addequation(get(child2,"val1")); if(eq == -1) emsg("prior prob5");
                  prioreq1.push_back(eq);
                  eq = addequation(get(child2,"val2")); if(eq == -1) emsg("prior prob6");
                  prioreq2.push_back(eq);
                  fl = 1;
                }
  
                if(type == "Exponential"){
                  priortype.push_back(EXPO);
                  priorparam.push_back(p);
                  priorminval.push_back(0);
                  priormaxval.push_back(0);
                  eq = addequation(get(child2,"val1")); if(eq == -1) emsg("prior prob7");
                  prioreq1.push_back(eq);
                  prioreq2.push_back(-1);
                  fl = 1;
                }
  
                if(type == "Beta"){
                  priortype.push_back(BETA);
                  priorparam.push_back(p);
                  priorminval.push_back(0);
                  priormaxval.push_back(0);
                  eq = addequation(get(child2,"val1")); if(eq == -1) emsg("prior prob9");
                  prioreq1.push_back(eq);
                  eq = addequation(get(child2,"val2")); if(eq == -1) emsg("prior prob10");
                  prioreq2.push_back(eq);
                  fl = 1;
                }
  
                if(type == "Weibull"){
                  priortype.push_back(WEIBULL);
                  priorparam.push_back(p);
                  priorminval.push_back(0);
                  priormaxval.push_back(0);
                  eq = addequation(get(child2,"val1")); if(eq == -1) emsg("prior prob11");
                  prioreq1.push_back(eq);
                  eq = addequation(get(child2,"val2")); if(eq == -1) emsg("prior prob12");
                  prioreq2.push_back(eq);
                  fl = 1;
                }
  
                if(fl == 0) emsg("prior type prob");
              }
            }
          }
        }
      }
    }
  }

  nprior = priorparam.size();
  nparam = paramname.size();
  npopnum = popnumname.size();

  /*
  for(p = 0; p < nparam; p++){
    cout << paramname[p] << " na\n";
    for(j = 0; j < nparampriordep[p]; j++){
      cout << parampriordep
    }
    
    nparampriordep;                        // the priors which depend on that parameter
    vector< vector<long> > parampriordep;
    
  }
  */
  sortparam();

  // smoothing priors
  nsmooth = 0;
  for(child = root->FirstChild(); child; child = child->NextSibling()){
    s1 = child->Value();
    if(s1.compare("model") == 0){  // model
      for(child2 = child->FirstChild(); child2; child2 = child2->NextSibling()){
        s2 = child2->Value();
        if(s2.compare("agesmooth") == 0 || s2.compare("timesmooth") == 0 ){
          par = get(child2,"parameter");
          if(get(child2,"type") == "Smooth") ty = SMOOTH;
          else{
            if(get(child2,"type") == "Log Smooth") ty = LOGSMOOTH;
            else emsg("Smooth prior problem");
          }
          val = getnum(child2,"val");

          nsmoothst = nsmooth;
          if(s2.compare("agesmooth") == 0){
            for(a = 0; a <= nage; a++){
              findparam(par,classval[agecl][a]);
              for(j = 0; j < plist.size(); j++){
                s = nsmoothst; while(s < nsmooth && smoothref[s] != allotherlist[j]) s++;
                if(s == nsmooth){
                  smoothname.push_back(par+" Age Smooth");
                  smoothtype.push_back(ty);
                  smoothref.push_back(allotherlist[j]);
                  smoothval.push_back(vector<double> ());
                  smoothparam.push_back(vector<long> ());
                  nsmooth++;
                }
                if(a == nage) smoothval[s].push_back(0);
                else smoothval[s].push_back(val);
                smoothparam[s].push_back(plist[j]);
              }
            }
            for(s = nsmoothst; s < nsmooth; s++){
              if(smoothparam[s].size() != nage+1) emsg("Smooth problem");
            }
          }
          else{
            for(a = 0; a <= nsettime; a++){
              findparam(par,classval[settimecl][a]);
              for(j = 0; j < plist.size(); j++){
                s = nsmoothst; while(s < nsmooth && smoothref[s] != allotherlist[j]) s++;
                if(s == nsmooth){
                  smoothname.push_back(par+" Time Smooth");
                  smoothtype.push_back(ty);
                  smoothref.push_back(allotherlist[j]);
                  smoothval.push_back(vector<double> ());
                  smoothparam.push_back(vector<long> ());
                  nsmooth++;
                }
                if(a == nsettime) smoothval[s].push_back(0);
                else smoothval[s].push_back(val);
                smoothparam[s].push_back(plist[j]);
              }
            }

            for(s = nsmoothst; s < nsmooth; s++){
              if(smoothparam[s].size() != nsettime+1) emsg("Smooth problem");
            }
          }

          for(s = nsmoothst; s < nsmooth; s++){ nsmoothparam.push_back(smoothparam[s].size());}
        }
      }
    }
  }

  paramsmooth.resize(nparam); paramsmoothi.resize(nparam);
  for(s = 0; s < nsmooth; s++){
    for(j = 0; j < nsmoothparam[s]; j++){
      paramsmooth[smoothparam[s][j]].push_back(s);
      paramsmoothi[smoothparam[s][j]].push_back(j);
    }
  }

  
  for(s = 0; s < nsmooth; s++){
    cout << smoothname[s] << " " << smoothtype[s] << " " << smoothref[s] << "  :";
    for(j = 0; j < smoothparam[s].size(); j++) cout << paramname[smoothparam[s][j]] << "," <<  smoothval[s][j] << ",   ";
    cout << "smooth\n";
  }
  

  //for(p = 0; p < nparam; p++) cout << paramname[p] << " param\n";
  //for(p = 0; p < nprior; p++) cout << paramname[priorparam] << " priorparam\n";
  /*
  ncompleave.resize(nclass); compleave.resize(nclass);
  for(cl = 0; cl < nclass; cl++){
    if(cl == agecl) jmax = nclassval[cl]+1; else jmax = nclassval[cl];
    ncompleave[cl].resize(jmax); compleave[cl].resize(jmax);
    for(tr = 0; tr < ntrans[cl]; tr++){
      i = transi[cl][tr]; if(i == -1) i = nclassval[cl];
      compleave[cl][i].push_back(tr);
    }

    for(j = 0; j < jmax; j++) ncompleave[cl][j] = compleave[cl][j].size();
  }

  // shifts all the times
  //tshift = tmin;
  */

  //for(j = 0; j < ncapevtrange; j++) cout << j << "  "<< capevtrange[j] << " rran\n"; 

  for(ob = 0; ob < nobs; ob++) obst[ob] -= tmin;
  for(fi = 0; fi < nfixev; fi++) fixev[fi].t -= tmin;
  for(i = 0; i < nind; i++){ if(indfixtenter[i] != large) indfixtenter[i] -= tmin;}
  for(i = 0; i < nind; i++){ if(indfixtleave[i] != large) indfixtleave[i] -= tmin;}

  for(cap = 0; cap < ncap; cap++) capt[cap] -= tmin;
  tmax -= tmin; tmax2 -= tmin;
  for(j = 0; j < nsettime; j++) settime[j] -= tmin;
  for(j = 0; j < ncapevtrange; j++) capevtrange[j] -= tmin;
  for(p = 0; p < npopm; p++) popmt[p] -= tmin;
  for(dm = 0; dm < nderm; dm++) dermt[dm] -= tmin;
}

void createtra(long type, XMLNode* child3, long cl, long i, long f)                    // gets all the exponential transition from the file
{
  long c, j, eq;
  string s4;
  vector <long> eqmap, eqmap2;
  XMLNode* child4;

  eqmap.resize(ncomp); eqmap2.resize(ncomp); for(c = 0; c < ncomp; c++){ eqmap[c] = -1; eqmap2[c] = -1;}

  for(child4 = child3->FirstChild(); child4; child4 = child4->NextSibling()){
    s4 = child4->Value();
    if(s4.compare("value") == 0){
      switch(type){
        case EXPONENTIAL:
          getalltrans(child3,child4,cl,i,get(child4,"rate"));
          for(j = 0; j < state.size(); j++) eqmap[state[j]] = addequation(staterate[j]);
          break;

        case SOURCE: case SINK:
          getalltrans(child3,child4,-1,-1,get(child4,"rate"));
          for(j = 0; j < state.size(); j++) eqmap[state[j]] = addequation(staterate[j]);
          break;

        case GAM:
          getalltrans(child3,child4,cl,i,get(child4,"mean"));
          for(j = 0; j < state.size(); j++) eqmap[state[j]] = addequation(staterate[j]);

          getalltrans(child3,child4,cl,i,get(child4,"shape"));
          for(j = 0; j < state.size(); j++) eqmap2[state[j]] = addequation(staterate[j]);
          break;

        case WEI:
          getalltrans(child3,child4,cl,i,get(child4,"λ"));
          for(j = 0; j < state.size(); j++) eqmap[state[j]] = addequation(staterate[j]);

          getalltrans(child3,child4,cl,i,get(child4,"k"));
          for(j = 0; j < state.size(); j++) eqmap2[state[j]] = addequation(staterate[j]);
          break;
      }
    }
    else emsg("rate not there");
  }

  for(c = 0; c < ncomp; c++){
    eq = eqmap[c];
    if(eq >= 0){
      TT tt;
      switch(type){
        case EXPONENTIAL: tt.type = EXP_TR; tt.ci = c; tt.cf = c+(f-i)*classmult[cl]; break;
        case SOURCE: tt.type = EXP_TR; tt.ci = NOTALIVE; tt.cf = c; break;
        case SINK: tt.type = EXP_TR; tt.ci = c; tt.cf = NOTALIVE; break;
        case GAM: tt.type = GAMMA_TR; tt.ci = c; tt.cf = c+(f-i)*classmult[cl]; tt.eqshape = eqmap2[c]; break;
        case WEI: tt.type = WEI_TR; tt.ci = c; tt.cf = c+(f-i)*classmult[cl]; tt.eqshape = eqmap2[c]; break;
      }
      tt.cl = cl; tt.i = i; tt.f = f; tt.eq = eq; tt.capev = -1; tt.like = 1;
      tra.push_back(tt);
      ntra++;
    }
  }
}

vector<string> getallpos(string name)    // If class specifiers are given then splits up
{
  long cl, p, j, jj, k, leng;
  vector<string> pos;
  string st;

  j = 0; while(j < len(name) && name.substr(j,1) != "_") j++;

  pos.push_back(name);
  if(j < len(name)){
    for(cl = 0; cl < nclass; cl++){
      vector<string> posnew;
      for(p = 0; p < pos.size(); p++){
        st = pos[p];
        leng = len(classname[cl]);
        jj = j+1;
        while(jj <= len(st)-leng && !(st.substr(jj,leng) == classname[cl] && (st.substr(jj-1,1) == "," || st.substr(jj-1,1) == "_") && (jj+leng == len(st)
              || st.substr(jj+leng,1) == ","))) jj++;
        if(jj <= len(st)-leng){
          for(k = 0; k < nclassval[cl]; k++) posnew.push_back(repla(st,classname[cl],classval[cl][k]));
        }
        else posnew.push_back(st);
      }
      pos = posnew;
    }
  }
  return pos;
}

string get(XMLNode* node, string attr)                                  // This gets an attribute
{
  const char *a = node->ToElement()->Attribute(attr.c_str());
  if(a){
    string s(a); return s;
  }
  else{
    stringstream ss; ss << "Cannot find attribute get " << attr;
    emsg(ss.str());
  }
}

double getnum(XMLNode* node, string attr)                          // Gets a number fron an XML attribute
{
  double num;
  const char *a = node->ToElement()->Attribute(attr.c_str());
  if(a){
    num = atof(a);
    stringstream ss; ss << a << "Not a number in attribute" << attr;
    if(isnan(num)) emsg(ss.str());
    return num;
  }
  else{
    stringstream ss; ss << "Cannot find attribute " << attr;
    emsg(ss.str());
  }
}


long exist(XMLNode* node, string attr)                            // Determines if an attribute exists
{
  const char *a = node->ToElement()->Attribute(attr.c_str());
  if(a) return 1;
  else return 0;
}

long getint(XMLNode* node, string attr)                           // Gets an integer XML attribute
{
  long num;
  const char *a = node->ToElement()->Attribute(attr.c_str());
  if(a){
    num = atoi(a);
    stringstream ss; ss << a << "Not a number in attribute" << attr;
    if(isnan(num)) emsg(ss.str());
    return num;
  }
  else{
    stringstream ss; ss << "Cannot find attribute " << attr;
    emsg(ss.str());
  }
}

vector<string> getlinesep(string a)           // Seperates using the | mark
{
  long j, jst;
  vector<string> linesep;

  j = 0;
  while(j < a.length()){
    jst = j;
    while(j < a.length() && a.substr(j,1) != "|") j++;
    linesep.push_back(a.substr(jst,j-jst));
    j++;
  }

  return linesep;
}

long getparam(string na)                                                      // Gets the parametere number fron the name
{
  long p;

  p = 0; while(p < paramname.size() && na != paramname[p]) p++;
  if(p == paramname.size()){
    paramname.push_back(na);
    paramprior.push_back(-1);
    paramclass.push_back("");
    parampriordep.push_back(vector<long>());

    /*
    cout << "beg" <<  na<< " cannot find\n";
    stringstream ss; ss << na << ": cannot find paramerter\n";
    for(p = 0; p < paramname.size(); p++) ss << paramname[p] << " jj\n";
    emsg(ss.str());
    */
  }
  return p;
}

string repla(string st, string sub1, string sub2)                              // Replaces one substring with another
{
  long i;

  i = 0;
  while(i <= len(st)-len(sub1)){
    if(st.substr(i,len(sub1)) == sub1){
      stringstream ss;
      ss << st.substr(0,i) << sub2 << st.substr(i+len(sub1));
      st = ss.str();
      i += len(sub2);
    }
    else i++;
  }

  return st;
}

void getalltrans(XMLNode* node, XMLNode* node2, long cl, long in, string rate)   // works out all the compartments consistent with a transition
{
  long cl2, filt[nclass], count[nclass], c, j, flag, i, ii;
  string st, s, name;

  state.clear(); staterate.clear();

  for(cl2 = 0; cl2 < nclass; cl2++){
    if(cl2 == cl) filt[cl2] = in;
    else{
      if(exist(node,classname[cl2].c_str()) == 1){
        st = get(node,classname[cl2].c_str());
        j = 0; while(j < nclassval[cl2] && classval[cl2][j] != st) j++;
        if(j == nclassval[cl2]){ cout << st << "st\n"; emsg("yy2");}
        filt[cl2] = j;
      }
      else{
        if(exist(node2,classname[cl2].c_str()) == 1){
          st = get(node2,classname[cl2].c_str());
          j = 0; while(j < nclassval[cl2] && classval[cl2][j] != st) j++;
          if(j == nclassval[cl2]){ cout << st << "st\n"; emsg("yy4");}
          filt[cl2] = j;
        }
        else filt[cl2] = -1;
      }
    }
  }
  //getalltrans(child2,child2,-1,-1,"");

  //if(cl == -1){ for(cl2 = 0; cl2 < nclass; cl2++) cout << filt[cl2] << ","; cout << " filt\n";}

  for(cl2 = 0; cl2 < nclass; cl2++){ if(filt[cl2] == -1) count[cl2] = 0; else count[cl2] = filt[cl2];}

  do{
    c = 0; for(cl2 = 0; cl2 < nclass; cl2++) c += classmult[cl2]*count[cl2];

    state.push_back(c);

    s = rate;                                                         // modifies expression to replace individual based classifiers with their actual values
    for(cl2 = 0; cl2 < nclass; cl2++){
      name = classname[cl2];

      ii = 0;
      do{
        if(s.substr(ii,1) == "{"){
          do{
            if(s.substr(ii,name.length()) == name && (s.substr(ii+len(name),1) == "," || s.substr(ii+len(name),1) == "}")){
              s = s.substr(0,ii)+classval[cl2][(c/classmult[cl2])%nclassval[cl2]]+s.substr(ii+len(name));
              ii += 2;
            }
            if(s.substr(ii,1) == "}") break;
            ii++;
          }while(ii < s.length());
        }
        else{
          if(s.substr(ii,1) == "["){
            ii++;
            do{
              if(s.substr(ii,name.length()) == name && (s.substr(ii+len(name),1) == "," || s.substr(ii+len(name),1) == "]")
                 && (s.substr(ii-1,1) == "," || s.substr(ii-1,1) == "_")){
                s = s.substr(0,ii)+classval[cl2][(c/classmult[cl2])%nclassval[cl2]]+s.substr(ii+len(name));
                ii += 2;
                 }
                 if(s.substr(ii,1) == "]") break;
                 ii++;
            }while(ii < s.length());
          }
          else ii++;
        }
      }while(ii < s.length());
    }
    //for(cl2 = 0; cl2 < nclass; cl2++) cout << count[cl2] << ","; cout << "   " <<s <<  "  count\n";

    staterate.push_back(s);
    j = 0;
    do{
      flag = 0;
      if(filt[j] != -1){ j++;flag = 1;}
      else{
        count[j]++; if(count[j] >= nclassval[j]){ count[j] = 0; j++; flag = 1;}
      }
    }while(flag == 1 && j < nclass);
  }while(j < nclass);
}

void findparam(string par, string s)      // Find the paramters with a particular name and attribute
{
  long p, j, jst;
  string name;

  allotherlist.clear(); plist.clear();
  for(p = 0; p < nparam; p++){
    name = paramname[p];
    j = 0; while(j < len(name) && name.substr(j,1) != "_") j++;
    //cout << name.substr(0,j) <<
    if(name.substr(0,j) == par){
      j++;
      while(j < len(name)){
        jst = j; while(j < len(name) && name.substr(j,1) != ",") j++;
        if(s == name.substr(jst,j-jst)){
          allotherlist.push_back(name.substr(0,jst)+name.substr(j,len(name)-j));
          plist.push_back(p);
        }
        j++;
      }
    }
  }
  if(plist.size() == 0) emsg("Could not set up smoothing prior");
}

struct PORD{ long param; long prior;};
bool pordcomp(PORD lhs, PORD rhs) { return lhs.prior < rhs.prior; }
void sortparam()
{
  long p;
  PORD paramord[nparam];

  for(p = 0; p < nparam; p++){ paramord[p].param = p; paramord[p].prior = paramprior[p];}
  sort(paramord,paramord+nparam,pordcomp);
  for(p = 0; p < nparam; p++) paramorder.push_back(paramord[p].param);
}

struct CA { string name; double time; vector <long> eqs;};
bool compare(CA lhs, CA rhs) { return lhs.time < rhs.time; }
void sortcap()
{

  long cap, c;
  CA caps[ncap];

  for(cap = 0; cap < ncap; cap++){
    caps[cap].name = capname[cap];
    caps[cap].time = capt[cap];
    caps[cap].eqs.resize(ncomp); for(c = 0; c < ncomp; c++) caps[cap].eqs[c] = capprobeqn[cap][c];
  }

  sort(caps,caps+ncap,compare);

  for(cap = 0; cap < ncap; cap++){
    capname[cap] = caps[cap].name;
    capt[cap] = caps[cap].time;
    for(c = 0; c < ncomp; c++) capprobeqn[cap][c] = caps[cap].eqs[c];
  }
}

struct OB { long i; double time; long cap; vector <long> eqs;};
bool compare2(OB lhs, OB rhs) { return lhs.time < rhs.time; }
void sortobs()
{
  long ob, c;
  OB obss[nobs];

  for(ob = 0; ob < nobs; ob++){
    obss[ob].i = obsi[ob];
    obss[ob].time = obst[ob];
    obss[ob].cap = obscap[ob];
    obss[ob].eqs.resize(ncomp); for(c = 0; c < ncomp; c++)  obss[ob].eqs[c] = obsprobeqn[ob][c];
  }

  sort(obss,obss+nobs,compare2);

  for(ob = 0; ob < nobs; ob++){
    obsi[ob] = obss[ob].i;
    obst[ob] = obss[ob].time;
    obscap[ob] = obss[ob].cap;
    for(c = 0; c < ncomp; c++) obsprobeqn[ob][c] = obss[ob].eqs[c];
  }
}

struct PO { double t; double al; double be; vector <long> comp;};
bool compare4(PO lhs, PO rhs) { return lhs.t < rhs.t; }
void sortpopm()
{
  long p, c;
  PO pops[npopm];

  for(p = 0; p < npopm; p++){
    pops[p].t = popmt[p];
    pops[p].al = popmal[p];
    pops[p].be = popmbe[p];
    pops[p].comp.resize(ncomp); for(c = 0; c < ncomp; c++)  pops[p].comp[c] = popmcomp[p][c];
  }

  sort(pops,pops+npopm,compare4);

  for(p = 0; p < npopm; p++){
    popmt[p] = pops[p].t;
    popmal[p] = pops[p].al;
    popmbe[p] = pops[p].be;
    for(c = 0; c < ncomp; c++) popmcomp[p][c] = pops[p].comp[c];
  }
}

struct DO { double t; double av; double var;};
bool compare5(DO lhs, DO rhs) { return lhs.t < rhs.t; }
void sortderm()
{
  long dm;
  DO ders[nderm];

  for(dm = 0; dm < nderm; dm++){
    ders[dm].t = dermt[dm];
    ders[dm].av = dermav[dm];
    ders[dm].var = dermvar[dm];
  }

  sort(ders,ders+nderm,compare5);

  for(dm = 0; dm < nderm; dm++){
    dermt[dm] = ders[dm].t;
    dermav[dm] = ders[dm].av;
    dermvar[dm] = ders[dm].var;
  }
}

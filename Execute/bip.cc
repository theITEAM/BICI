// Compile using:    g++ bip.cc tinyxml2.cc -O3 

// g++ bip.cc tinyxml2.cc  -g    use degugger
// gdb ./a.out
// run
// bt
//tr += moventra;

long noout = 0;            // Suppresses the output
const long checkon = 0;    // Determines if checking is done
long simnum = 0;

#include "tinyxml2.h"

using namespace tinyxml2;

using namespace std;

long ninput;
double input[10];

#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
#include <sstream>
#include <vector>
#include <algorithm>

//#include <cerrno>

#include "var.h"
#include "equation.h"
#include "readinput.h"
#include "incompgamma.h"
#include "dist.h"
#include "sim.h"
#include "output.h"
#include "init.h"
//#include "initproc.h"
#include "startseq.h"
#include "prior.h"
#include "likelihoodinit.h"
#include "likelihood.h"
#include "check.h"
#include "ratecalc.h"
#include "observation.h"
#include "param_prop.h"
#include "part_propnew.h"
#include "local_prop.h"
#include "move_prop.h"
#include "life_prop.h"
#include "addrem_prop.h"
#include "indsim_prop.h"
#include "derive.h"

int main (int argc, char *argv[])
{
  long i;
  string file;
  vector <double> v;

  if(argc == 1){ noout = 1; file = "init.xml"; simnum = 0;}
  if(argc > 1){ file = argv[1]; simnum = atoi(argv[2]);}                      // Creates a different see for different chains

  for(i = 0; i < 5555*simnum; i++) ran();

  cout << "3|0|\n"; cout.flush();

  readinput(file); cout << "3|10|\n"; cout.flush();

  init(); cout << "3|15|\n"; cout.flush();

  likelihoodinit(); cout << "3|20|\n"; cout.flush();

  if(derfl == 1) deriveinit();
  priorinit(); cout << "3|30|\n"; cout.flush();
  observationinit(); cout << "3|40|\n"; cout.flush();
  paraminit();
  partinit();
  localpropinit();
  lifepropinit(); cout << "3|50|\n"; cout.flush();
  moveinit();
  if(simon == 0) startseqinit(); cout << "3|60|\n"; cout.flush();

  outputmodel();

  ch[0] = new Chain;
  nchain = 1;

  ch[0]->invT_pd = 1;
  ch[0]->invTLi = 1;
  ch[0]->invT_st = 1;
  ch[0]->invT_pop = 20;
  ch[0]->invT_der = 10;

  //checkdist(); emsg("don");
  siminit();

  //ch[0]->initproc();
  //emsg("do");
  ch[0]->initparamsamp();

  ch[0]->initchain(); 
 
  if(simon == 0) ch[0]->start();
  
  cout << "3|90|\n"; cout.flush();
  
  cout << "3|95|\n"; cout.flush();

  if(nderive > 0) ch[0]->deriveplotinit();

  //siminit();
  cout << "3|100|\n"; cout.flush();

  traceinit();
  if(simon == 1){
    cout << "8|\n"; cout << "3|0|\n"; cout.flush();
    for(samp = 0; samp < nsamp; samp++){ 
      ch[0]->eventplot();
      ch[0]->traceplot();
    }
    //if(nderive > 0){cout << "HHHH\n"; ch[0] -> deriveplotcalc(); ch[0] -> derivepl();}
    return 0;
  }

  ch[0]->check(-1);

  //short p; for(p = 0; p < nparam; p++) cout << paramname[p] << ":" << ch[0]->param[p] << ",  "; cout << "\n"; return 0;
  for(samp = 0; samp < nsamp; samp++){
    altertemperature();

    //cout << samp << " "<< ch[0]->invT_pop << " temp\n";
    //if(samp < 100){ ch[0]->invT_pop = 10; ch[0]->invT_der = 10;}else{ ch[0]->invT_pop = 1; ch[0]->invT_der = 1;}

    if(noout == 1 && samp%100 == 0) cout << samp << " samp\n";

    if(samp%50 == 0) ch[0]->eventplot();

    if(samp%1 == 0) ch[0]->traceplot();

    ch[0]->update();

    if(samp%50 == 0) ch[0]->check(0);

    if(samp != 0 && samp%100 == 0) ch[0]->diagnosticschain();
  }
  ch[0]->diagnosticschain();
}

void Chain::update()                                                // Performs an MCMC update
{
  long chon = 0, paramloop, loop;
  long i, imax;
  
  //if(samp == 99) chon = 1;

  if(samp < burnin && samp%5 == 0) optimiseproposals();
  if(samp < burnin && samp%10 == 0) simsumcalc();

  if(totaltime == 0) paramloop = 4;
  else{ paramloop = long(5-8*timeprop[PARAM_PROP]/totaltime); if(paramloop < 1) paramloop = 1;}
 
  totaltime -= clock();

  timeprop[PARAM_PROP] -= clock();
  param_prop();
  for(loop = 0; loop < paramloop; loop++) paramnorm_prop();
  timeprop[PARAM_PROP] += clock();
  ntimeprop[PARAM_PROP]++;

  if(chon == 1) check(1);

  if(ran() < part_prob){
    timeprop[PART_PROP] -= clock();
    for(i = 0; i < nind; i++) part_prop(i);
    timeprop[PART_PROP] += clock();
    ntimeprop[PART_PROP]++;
  }

  if(chon == 1) check(2);

  if(nindtot > nind){
    if(ran() < induo_prob){
      timeprop[INDSIMUO_PROP] -= clock();
      for(i = nind; i < nindtot; i++) indsim_prop(i);
      timeprop[INDSIMUO_PROP] += clock();
      ntimeprop[INDSIMUO_PROP]++;
    }
  }

  if(chon == 1) check(3);

  if(ran() < move_prob){
    timeprop[MOVE_PROP] -= clock();
    for(i = 0; i < nind; i++) move_prop(i);
    timeprop[MOVE_PROP] += clock();
    ntimeprop[MOVE_PROP]++;
  }

  if(chon == 1) check(4);

  if(singeventfl == 1){
    if(ran() < sing_prob){
      timeprop[SINGEVENT_PROP] -= clock();
      for(i = 0; i < nind; i++) sing_prop(i);
      timeprop[SINGEVENT_PROP] += clock();
      ntimeprop[SINGEVENT_PROP]++;
    }
  }

  if(chon == 1) check(5);

  if(twothreefl == 1){
    if(ran() < twothree_prob){
      timeprop[TWOTHREE_PROP] -= clock();
      for(i = 0; i < nind; i++) twothree_prop(i);
      timeprop[TWOTHREE_PROP] += clock();
      ntimeprop[TWOTHREE_PROP]++;
    }
  }

  if(chon == 1) check(6);

  if(pairfl == 1){
    if(ran() < pair_prob){
      timeprop[PAIR_PROP] -= clock();
      for(i = 0; i < nind; i++) pair_prop(i);
      timeprop[PAIR_PROP] += clock();
      ntimeprop[PAIR_PROP]++;
    }
  }

  if(chon == 1) check(7);

  if(gapfl == 1){
    if(ran() < gap_prob){
      timeprop[GAP_PROP] -= clock();
      for(i = 0; i < nind; i++) gap_prop(i);
      timeprop[GAP_PROP] += clock();
      ntimeprop[GAP_PROP]++;
    }
  }

  if(chon == 1) check(8);

  if(tbirthfl == 1 || sourcefl == 1 || sinkfl == 1){
    if(ran() < life_prob){
      timeprop[LIFE_PROP] -= clock();
      for(i = 0; i < nindtot; i++) life_prop(i);
      timeprop[LIFE_PROP] += clock();
      ntimeprop[LIFE_PROP]++;
    }
  }

  if(chon == 1) check(9);

  if(addremfl == 1){
    if(ran() < addrem_prob){
      timeprop[ADDREM_PROP] -= clock();

      imax = (nind+nindtot)/4; if(imax == 0) imax = 1;
      for(i = 0; i < imax; i++) addrem_prop();
      timeprop[ADDREM_PROP] += clock();
      ntimeprop[ADDREM_PROP]++;
    }
  }

  if(chon == 1) check(10);

  totaltime += clock();
}

void Chain::optimiseproposals()                // optimises the frequency of different proposals
{
  long i;
  double fav, max;

  part_prob = 0.3;  induo_prob = 0.2; move_prob = 0.5; sing_prob = 0.2; twothree_prob = 0.2; pair_prob = 0.2; gap_prob = 0.2; life_prob = 1; addrem_prob = 1;

  
  //ind_prob = 0.2; 
  //part_prob = 1; ind_prob = 0;  induo_prob = 0.4; move_prob = 0; sing_prob = 0; twothree_prob = 0; pair_prob = 0; gap_prob = 0; life_prob = 1; addrem_prob = 1;
  /*
  return;

  if(samp == 0){
    //param_prob = 2; 
    part_prob = 1;; move_prob = 1; sing_prob = 1; twothree_prob = 1; pair_prob = 1; gap_prob = 1; life_prob = 1; addrem_prob = 1;
  }
  else{
    //part_prob = 1;

    fav = 0; for(i = 0; i < nind; i++) fav += (nac_part[i]+0.03)/(ntr_part[i]+0.1)/nind;
    part_prob = 1*fav*(ntimeprop[PART_PROP]+1)/(timeprop[PART_PROP]+1);

    move_prob = 1*0.5*(ntimeprop[MOVE_PROP]+1)/(timeprop[MOVE_PROP]+1);

    if(singeventfl == 0) sing_prob = 0;
    else{
      fav = 0; for(i = 0; i < nind; i++) fav += (nac_sing[i]+0.03)/(ntr_sing[i]+0.1)/nind;
      sing_prob = 0.3*fav*(ntimeprop[SINGEVENT_PROP]+1)/(timeprop[SINGEVENT_PROP]+1);
    }

    if(twothreefl == 0) twothree_prob = 0;
    else{
      fav = 0; for(i = 0; i < nind; i++) fav += (nac_twothree[i]+0.03)/(ntr_twothree[i]+0.1)/nind;
      twothree_prob = 0.3*fav*(ntimeprop[TWOTHREE_PROP]+1)/(timeprop[TWOTHREE_PROP]+1);
    }

    if(pairfl == 0) pair_prob = 0;
    else{
      fav = 0; for(i = 0; i < nind; i++) fav += (nac_pair[i]+0.03)/(ntr_pair[i]+0.1)/nind;
      pair_prob = 0.3*fav*(ntimeprop[PAIR_PROP]+1)/(timeprop[PAIR_PROP]+1);
    }

    if(gapfl == 0) gap_prob = 0;
    else{
      fav = 0; for(i = 0; i < nind; i++) fav += (nac_gap[i]+0.03)/(ntr_gap[i]+0.1)/nind;
      gap_prob = 0.3*fav*(ntimeprop[GAP_PROP]+1)/(timeprop[GAP_PROP]+1);
    }

    if(!(tbirthfl == 1 || sourcefl == 1 || sinkfl == 1)) life_prob = 0;
    else life_prob = 0.5*0.333*(ntimeprop[LIFE_PROP]+1)/(timeprop[LIFE_PROP]+1);

    if(!(capfl == 1 || popfl == 1 || derfl == 1)) addrem_prob = 0;
    else addrem_prob = 0.5*0.333*(ntimeprop[ADDREM_PROP]+1)/(timeprop[ADDREM_PROP]+1);
    cout << addrem_prob << " " << ntimeprop[ADDREM_PROP] << " " << timeprop[ADDREM_PROP] << "Pr\n";
    max = 0;

    //if(param_prob > max) max = param_prob;
    if(part_prob > max) max = part_prob;
    if(move_prob > max) max = move_prob;
    if(sing_prob > max) max = sing_prob;
    if(twothree_prob > max) max = twothree_prob;
    if(pair_prob > max) max = pair_prob;
    if(gap_prob > max) max = gap_prob;
    if(life_prob > max) max = life_prob;
    if(addrem_prob > max) max = addrem_prob;

    //param_prob /= max; 
    part_prob /= max; move_prob /= max; 
    sing_prob /= max; twothree_prob /= max; pair_prob /= max; gap_prob /= max; life_prob /= max; addrem_prob /= max;
  }
  */
}

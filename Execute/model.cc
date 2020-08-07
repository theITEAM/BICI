
using namespace std;

#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
#include <sstream>
#include <vector>
#include <algorithm>

double ran(){ double v = (double(rand())*29999+rand())/(30000.0*RAND_MAX); if(v == 0) return 0.1; else return v;}

void diseasetrans()
{
  ofstream init("Initpop1.txt");
  ofstream model("Model.txt");
  long t, indtot = 1, i;
  double z, p=0.4;

 
  for(t = 0; t < 64; t++){
    p = ran();

    model << "compartment classification='Group' name = 'G"<< (t+1) << "' x='"<< t%8 <<"' y='" << t/8 << "' color='#aaaaaa'\n";

    for(i = 0; i < 10; i++){
      if(i == 0) init << "Ind. " << indtot << "\t" << "G" << (t+1) << "\t" << "I\t";
      else init << "Ind. " << indtot << "\t" << "G" << (t+1) << "\t" << "S\t";

      z = ran(); 
      if(z < p*p) init << "AA\n";
      else{
        if(z-p*p < (1-p)*(1-p)) init << "BB\n";
        else  init << "AB\n";
      }

      /*
      if(t < 3){
      if(i == 0) init << "Ind. " << indtot << "\t" << "G" << (t+1) << "\t" << "I\tV\n";
      else init << "Ind. " << indtot << "\t" << "G" << (t+1) << "\t" << "S\tV\n";
    }
      else{
      if(i == 0) init << "Ind. " << indtot << "\t" << "G" << (t+1) << "\t" << "I\tNV\n";
      else init << "Ind. " << indtot << "\t" << "G" << (t+1) << "\t" << "S\tNV\n";
    }
  */
        indtot++;
    }
  }
}

void spatialdiffusion()
{
  const short GX = 11, GY = 11, range = 3;
  short i, j, ii, jj;
  double d;
  
  ofstream model("GP/Model 2.txt");
  model << "# This model uses a grid of compartments to represent a spatially distributed set of locations\n";

  for(i = 0; i < GX; i++){
    for(j = 0; j < GY; j++){
      
      model << "compartment classification='Loc' name = 'C" << i << "-" << j << "' x='" << i << "' y='" << j;
      
      if((i == 1 && j == 1) || (i == 1 && j == 5) || (i == 1 && j == 9) ||
          (i == 5 && j == 1) || (i == 5 && j == 5) || (i == 5 && j == 9) ||
          (i == 9 && j == 1) || (i == 9 && j == 5) || (i == 9 && j == 9)) model <<"' color='#0000ff'\n";
      else  model << "' color='#ff0000'\n";
      
  
      //model << "' color='#ff0000'\n";
    }
  }

  model << "\n";

  model << "# Individuals can move from one location to a neighbouring location\n";

  for(i = 0; i < GX; i++){
    for(j = 0; j < GY; j++){
      for(ii = i-1; ii <= i+1; ii++){
        for(jj = j-1; jj <= j+1; jj++){
          if(ii >= 0 && ii < GX && jj >= 0 && jj < GY && !(i == ii && j == jj)){
            model << "transition  from='C" << i << "-" << j << "' to='C" << ii << "-" << jj << "' type='exponential' ";
            if(i == ii || j == jj) model << "rate='[w]'\n";
            else model << "rate='[w]/2'\n";
          }
        }
      }
    }
  }

  model << "\n";

  model << "# For simulation initially 5 individuals start at each location\n"; 
  model << "siminitpopulation ";
  for(i = 0; i < GX; i++){
    for(j = 0; j < GY; j++){
      model << "C" << i << "-" << j << "='5' ";
    }
  }
  model << "\n";

  //model << "siminitpercent S='100' I='0'\n"; 
}

void spatialdiffusiondisease()
{
  const short GX = 11, GY = 11, range = 3;
  short i, j, ii, jj;
  double d;

  ofstream model("GP/Model 2.txt");
  model << "# This model uses a grid of compartments to represent a spatially distributed set of locations\n";

  for(i = 0; i < GX; i++){
    for(j = 0; j < GY; j++){
      model << "compartment classification='Loc' name = 'C" << i << "-" << j << "' x='" << i << "' y='" << j;
      /*
      if((i == 1 && j == 1) || (i == 1 && j == 5) || (i == 1 && j == 9) ||
          (i == 5 && j == 1) || (i == 5 && j == 5) || (i == 5 && j == 9) ||
          (i == 9 && j == 1) || (i == 9 && j == 5) || (i == 9 && j == 9)) model <<"' color='#0000ff'\n";
      else  model << "' color='#ff0000'\n";
      */
      model << "' color='#ff0000'\n";
    }
  }

  model << "\n";

  model << "# Individuals can move from one location to a neighbouring location\n";

  for(i = 0; i < GX; i++){
    for(j = 0; j < GY; j++){
      for(ii = i-1; ii <= i+1; ii++){
        for(jj = j-1; jj <= j+1; jj++){
          if(ii >= 0 && ii < GX && jj >= 0 && jj < GY && !(i == ii && j == jj)){
            model << "transition  from='C" << i << "-" << j << "' to='C" << ii << "-" << jj << "' type='exponential' ";
            if(i == ii || j == jj) model << "rate='[w]'\n";
            else model << "rate='[w]/2'\n";
          }
        }
      }
    }
  }

  model << "\n";

  model << "# Disease progression is represented an SI model\n";

  model << "compartment classification='DS' name = 'S' x='0' y='0'  color='#00ff00'\n";;
  model << "compartment classification='DS' name = 'I' x='1' y='0'  color='#ff0000'\n";;

  model << "\n";
  model << "transition from='S' to='I' type='exponential' rate='[a]+[β]*{I,Loc}'\n";
  model << "\n";

  model << "# For simulation initially 5 individuals start at each location all uninfected\n"; 
  model << "siminitpopulation ";
  for(i = 0; i < GX; i++){
    for(j = 0; j < GY; j++){
      model << "C" << i << "-" << j << "='5' ";
    }
  }
  model << "\n";

  model << "siminitpercent S='100'\n";
}

void spatialfarm()
{
  const double range = 10, distmin = 3.5;
  const short nfarm = 500;
  short i, j, fl;
  double d, fx[nfarm], fy[nfarm], dist;

  ofstream model("GP/Model 5.txt");
  model << "# This model uses points to represent a spatially distributed set of farms\n";

  for(i = 0; i < nfarm; i++){
    cout << i << " i\n";
    do{
      fx[i] = 100*ran(); fy[i] = 100*ran();
      for(j = 0; j < i; j++){
        dist = sqrt((fx[i]-fx[j])*(fx[i]-fx[j]) + (fy[i]-fy[j])*(fy[i]-fy[j]));
        if(dist < distmin) break;
      }
    }while(j < i);

    model << "compartment classification='Loc' name = 'F" << (i+1) << "' x='" << fx[i]<< "' y='" << fy[i];
    model << "' color='#aaaaaa'\n";
  }

  model << "\n";

  model << "# Each farm can be susceptible, infected or recovered\n";

  model << "compartment classification='DS' name = 'S' x='0' y='0' color='#00ff00'\n";
  model << "compartment classification='DS' name = 'I' x='1' y='0' color='#ff0000'\n";
  model << "compartment classification='DS' name = 'R' x='2' y='0' color='#0000ff'\n";
  model << "\n";

  model << "# The rate of becoming infected depends on the force of infection coming from nearby farms.\n";
  model << "transition  from='S' to='I' type='exponential' rate='\n";
  for(i = 0; i < nfarm; i++){
    fl = 0;
   
    switch(2){
      case 0: // Power
	model << "F" << i+1 << ":[a]+[w]*(";
	for(j = 0; j < nfarm; j++){
	  dist = sqrt((fx[i]-fx[j])*(fx[i]-fx[j]) + (fy[i]-fy[j])*(fy[i]-fy[j]));
	  if(dist > 0 && dist < range){
	    if(fl != 0) model << "+";
	    model << "{I,F" << j+1 << "}*(" << dist << "^[β])";
	    fl = 1;
	  }
	}
	model << ")";
	break;
	
      case 1:   // expoenential
	model << "F" << i+1 << ":[a]+[q]*(";
	for(j = 0; j < nfarm; j++){
	  dist = sqrt((fx[i]-fx[j])*(fx[i]-fx[j]) + (fy[i]-fy[j])*(fy[i]-fy[j]));
	  if(dist > 0 && dist < range){
	    if(fl != 0) model << "+";
	    model << "{I,F" << j+1 << "}*exp(-" << dist << "/[d0])";
	    fl = 1;
	  }
	}
	model << ")";
	break;
	
      case 2: //both
	model << "F" << i+1 << ":[a]+[s]*[w]*(";
	for(j = 0; j < nfarm; j++){
	  dist = sqrt((fx[i]-fx[j])*(fx[i]-fx[j]) + (fy[i]-fy[j])*(fy[i]-fy[j]));
	  if(dist > 0 && dist < range){
	    if(fl != 0) model << "+";
	    model << "{I,F" << j+1 << "}*(" << dist << "^[β])";
	    fl = 1;
	  }
	}
	model << ")";
	model << "+(1-[s])*[q]*(";
	fl = 0;
	for(j = 0; j < nfarm; j++){
	  dist = sqrt((fx[i]-fx[j])*(fx[i]-fx[j]) + (fy[i]-fy[j])*(fy[i]-fy[j]));
	  if(dist > 0 && dist < range){
	    if(fl != 0) model << "+";
	    model << "{I,F" << j+1 << "}*exp(-" << dist << "/[d0])";
	    fl = 1;
	  }
	}
	model << ")";
	break;
    }
    model << "\n";
  }
  model << "'\n";
  model << "\n";

  model << "# Farms recover at a given recovery rate γ.\n";
  model << "transition from='I' to='R' type='exponential' rate='[γ]'\n";
  model << "\n";

  model << "# For simulation initially 1 individual starts at each farm.\n"; 
  model << "siminitpopulation ";
  for(i = 0; i < nfarm; i++) model << "F" << (i+1) << "='1' ";
  model << "\n";

  model << "# All individuals are initially infected\n";
  model << "siminitpercent  S = '100' I = '0'\n";
  model << "\n";

  // model << "# For simulation initially 90 individuals start at teh central location\n"; 
  //model << "siminitpopulation C5-5='90'\n"; 
}

void badger()
{
  const double range = 32, distmin = 17;
  const short nfarm = 30;
  short i, j, fl;
  double d, fx[nfarm], fy[nfarm], dist;

  ofstream model("GP/Model 6.txt");
  model << "# This model uses points to represent a spatially distributed set of badger setts \n";

  for(i = 0; i < nfarm; i++){
    cout << i << " i\n";
    do{
      fx[i] = 100*ran(); fy[i] = 100*ran();
      for(j = 0; j < i; j++){
        dist = sqrt((fx[i]-fx[j])*(fx[i]-fx[j]) + (fy[i]-fy[j])*(fy[i]-fy[j]));
        if(dist < distmin) break;
      }
    }while(j < i);

    model << "compartment classification='Loc' name = 'L" << (i+1) << "' x='" << fx[i]<< "' y='" << fy[i];
    model << "' color='#aaaaaa'\n";
  }
  for(i = 0; i < nfarm; i++){
    for(j = 0; j < nfarm; j++){
      if(i != j){
        dist = sqrt((fx[i]-fx[j])*(fx[i]-fx[j]) + (fy[i]-fy[j])*(fy[i]-fy[j]));
        if(dist < range){
          model << "transition from='L" << (i+1)  << "' to='L" << (j+1) << "' type='exponential' rate='[ω]'\n";
        }
      }
    }
  }
  model << "\n";

  model << "# Disease progression is represented an SI model\n";

  model << "compartment classification='DS' name = 'S' x='0' y='0'  color='#00ff00'\n";;
  model << "compartment classification='DS' name = 'I' x='1' y='0'  color='#ff0000'\n";;
  model << "transition to='S' x='0' y ='-0.2' type='source' rate='[ν]*{Loc}*(1-{Loc}/[κ])'\n";
  model << "transition from='S' x='0' y ='0.2' type='sink' rate='[μ]'\n";
  model << "transition from='I' x='1' y ='0.2' type='sink' rate='[μ]*[m]'\n";
  model << "transition from='S' to='I' type='exponential' rate='[a]+[β]*{I,Loc}'\n";
  model << "\n";

  model << "# For simulation initially 5 badgers per sett with 30%  prevalence of disease.\n"; 
  model << "siminitpopulation ";
  for(i = 0; i < nfarm; i++) model << "L" << (i+1) << "='5' ";
  model << "\n";
  model << "siminitpercent  I = '30'\n";
  model << "\n";

  // model << "# For simulation initially 90 individuals start at teh central location\n"; 
  //model << "siminitpopulation C5-5='90'\n"; 
}

  
  // model << "# For simulation initially 90 individuals start at teh central location\n"; 
  //model << "siminitpopulation C5-5='90'\n";



  
  /*
  ofstream siminit("siminit.txt");
  for(i = 0; i < 60; i++){
    siminit << "ind" << i << "\t" << "T" << (1+i/9) << "\t";
    if(i%9 < 2) siminit << "I"; else siminit << "S";
    siminit << "\n";
}
  ofstream siminit("siminit.txt");
  long num = 1, loop;
  for(i = 0; i < GX; i++){
    for(j = 0; j < GY; j++){
      for(loop = 0; loop < 5; loop++){
        siminit << "Ind. " << num << "\t" <<  "C" << i << "-" << j << "\t";
        if(i == 5 && j == 5 && loop == 0) siminit << "I\n";
        else siminit << "S\n";
        num++;
}
}
}
  */

int main (int argc, char *argv[])
{
  //diseasetrans();
  //badger();
  
  //spatialdiffusiondisease();
  //diseasefarm();
  //spatialdiffusion();
  spatialfarm();
}


const double gammalim = 0.05;  // Limit on gamma / beta distribution

long samp, nsamp;
const long burnin = 1000;

long  num = 0;

short tempp=0;  //tempor

double av[100], nav[100];

long plfl = 0;

//double timeparam[10],timegradparam[10];

const long CHECK = 13, INDCHANGE = 14, TRACEPLOT = 15, EVENTPLOT = 16, PARTSIM = 17, PARTINIT = 18, INITBA = 19, INIT = 20;

const long PARAM_PROP = 0, PART_PROP = 1, MOVE_PROP = 2, LIFE_PROP = 5, ADDREM_PROP = 6, SINGEVENT_PROP = 7, TWOTHREE_PROP = 8, INDSIMUO_PROP = 10;
const long PAIR_PROP = 11, GAP_PROP = 12;

const long EXP_TR = 0, FIXED_TR = 1, GAMMA_TR = 2, WEI_TR = 3, GROW_TR = 4, SETTIME_TR = 5, CAPEVTRANGE_TR = 6, BEG_TR = 7, END_TR = 8, NULL_TR = 9;
const long FLAT = 0, GAMMA = 1, NORMAL = 2, LOGNORMAL = 3, EXPO = 4, BETA = 5, WEIBULL = 6, FIX = 7;
const long SMOOTH = 0, LOGSMOOTH = 1;

long DERIVEX = 400;

const long CAP = 0, POP = 1, DER = 2;
long indmax = 100000;
const double large = 10000000, mlarge = -large;
const double tiny = 0.000000001;
//const double evdtmin = 0.00000001; // smallest seperation between events
const double evdtmin = 0.00000001; // smallest seperation between events
const long logsummax = 30000;
const long depeqdivmax = 1000;
const long ncompcapdiv = 1000;
const double ratesmall = 0.00000001;  // the smallest rate size
const double obssmall = 0.00001;      // give a minimum limit on the obsevation probability
const double notobsdL = -1000;           // Give the penalty for a non-observation
const double onefac = 0.99999999;
const long npartmax = 40;            // The maximum particle number
const long npartstart = 10;           // The number of particles to start with

//ofstream nevpl("GP/nev3");

long  errorcode = -1;

double capprobest, muest;

long loopstart, loopstartmax = 1000; // Loop used when finding a valid initial starting event seqeucne 

double totaltime = 0, timeprop[30], ntimeprop[30];

long fixfl=0;                                       // This is set to one if any events are obseved in the data
long sinkfixfl=0;
long sourcefixfl=0;
long movefl=0;                                      // Set to one if the data contains any move, sourcemove, or sinkmove
long sinkfl=0, sourcefl=0;                          // Set to one if the model contains sinks or sources
long gammafl=0;                                     // Set to one if any gamma transitions in model
long nmfl=0;                                        // Set to one if any non-Markovian transitions in model
long capevfl = 0;                                   // Set to one if there is a probability associated with captures
long capfl=0;                                       // This is set to one if there are any captures
long popfl=0;                                       // This is set to one if there are any population measurements
long derfl=0;                                       // This is set to one if there are any derived measurements
long tbirthfl=0;                                    // This is set to one of there are more than one age classification
long pairfl=0;                                      // Set to when event propodals whihc go from zero events to two are used
long singeventfl=0;                                 // Set to one if single events can be inserted and deleted
long twothreefl=0;                                  // Switches between two and three events
long gapfl=0;                                       // Set to one if single events can be inserted and deleted
long addremfl=0;                                    // Set to one if there are unobserved individuals 

long simon=0;                                       // Set to one if a simulation in run

// Model parameters

long agecl, nage;                              // Characterises the age time line
vector<double> age;

long settimecl, nsettime;                      // Characterises the Time time line
vector<double> settime;

long nclass;                                        // # of classifications
vector<string> classname;                            // Names of classifications
vector<long> nclassval;                             // # of compartments in a classification
vector< vector<string> > classval;                   // Names of compartments
vector<long> classmult;

vector<string> compname;                            // Overall name 

vector<long> paramorder;                            // Orders the parameters for output

struct TT { 
  long type; long ci; long cf; long dc; long cl; long i; long f; long eq; long eqshape;
  long capev; long like; long trans; long ntraend; vector<long> traend; long nm;
  long nnonexptra; vector<long> nonexptra;   // gives the non exponential transition activated when individual undergoes transition
};

long ntra;                                          // The total number of transitions going between compartments
vector<TT> tra;
long moventra=0;                                    // Denotes which transitions are used for move events
long trabeg;                                        // Denotes which transitions are use to begin the timeline
long traend;                                        // Denotes which transitions are use to end the timeline
long tranull;                                       // Denotes a null transition

long ncomp;                                         // The total number of compartents
long ncomps;                                        // Does not include age, time and fix timelines
long ncompswa;                                      // Does not include time and fix timelines
long ncompswf;                                      // Does not include fix timelines
vector< vector <long> > compiftra;                  // [ci][cf] Gets the tra number
vector< vector <long> > compval;                    // [c][cl] Gets the class value for a particular compartment

vector <long> ncompleave;                           // Gives the tra leaving a compartment
vector < vector <long> > compleave;

vector< vector <long> > ncompclleave;               // Gives the tra leaving a compartment with transition in a particular cl
vector< vector < vector <long> > > compclleave;

vector < vector < vector <long> > > compleavesimcldep;
vector < vector < vector <long> > > compleavesimclnotdep;

vector <long> ncompleavedep;                        // Gives the tra leaving a compartment which are dependent (without taking into account fixed events)
vector < vector <long> > compleavedep;
vector <long> ncompleavenotdep;                     // Gives the tra leaving a compartment which are not dependent (without taking into account fixed events)
vector < vector <long> > compleavenotdep;


//vector <long> nnonexptra;                           // [tra][#] gives the non exponential transition activated when individual undergoes transition
//vector< vector <long> > nonexptra;

// Prior

long nprior;
vector<long> priortype;
vector<long> priorparam;
vector<double> priorminval;
vector<double> priormaxval;
vector<long> prioreq1;
vector<long> prioreq2;

vector<double> agprior;                              // The prior used for the birth time

vector <double> dirchal;                             // The Dirchlet prior placed on the initial composition

// Data

double tmax;                                         // The range over which inference is performed
double tmax2;                                        // Additional time used for forward prediction
//double tdatamin;                                     // The time when data starts
double tminactual, tmaxactual;                       // The total range over which outputs are given


long nind;                                          // The total number of observed individuals
vector<string> indid;                                // Individual ids
vector<long> nindobs;                               // Individual observations
vector< vector<long> > indobs;
vector< vector<double> > indobst;
vector<long> nindfixev;                             // lists all the fixed events associated with an individual
vector< vector<long> > indfixev;
vector<double> indfixtenter;                        // If the entry time of the individual is fixed then it is put here (else set to large)
vector<double> indfixtleave;                        // If the leave time of the individual is fixed then it is put here (else set to large)
vector<long> indlikeenter;                          // Determines if the enter likelihood is evaluated (for move events)
vector<long> indlikeleave;                          // Determines if the leave likelihood is evaluated (for move events)
vector<long> indfixentercapev;                      // The capev reference
vector<long> indfixleavecapev;                      // The capev reference

long nobs;                                          // The state observations
vector<double> obst;
vector<long> obsi;
vector<long> obscap;                                // stores the capture corrspending to a particular observation
vector< vector<long> > obsprobeqn;                  // [#obs][c]   gives the probability for observation in each of the model compartments

long ncope;                                         // This gives a list of all the euqation used in obsprobeqn
vector<long> cope;
vector< vector<long> > obsprob_cope;

struct FEV { long i; double t; long cl; long trai; long traf; long dc; long capev; long type; long like;};
bool compareFEV(FEV lhs, FEV rhs) { return lhs.t < rhs.t; }

long nfixev;                        // Stores all the fixed events given in the data
vector<FEV> fixev;

long ncapev;                                        // Event captures
vector<long> capevall;                              // Set if all all events are captured
vector<long> capevprobeqn;                          // The equation for the capture event probability

long capevcl, ncapevtrange;                         // Characterises the capture events time line
vector<double> capevtrange;

vector<string> capevname;
vector<double> capevtmin;
vector<double> capevtmax;
vector<long> capevtype;
vector< vector<long> > capevfilt;
const short CESOURCE = 0, CESINK = 1, CETRANS = 2;

long ncap;                                          // The capture campaigns
vector<string> capname;
vector<double> capt;
vector< vector<long> > capprobeqn;                  // [#cap][c]   gives the probability for capturing in each of the model compartments

long ncpe;                                          // Gives a list of equations used in capprobeqn
vector<long> cpe;
vector< vector<long> > capprob_cpe;

struct CC { long val; long ty; double t;};
bool comparecc(CC lhs, CC rhs) { return lhs.t < rhs.t; }
vector<long> ncompcap;                              // This gives the captures which occur in a given compartment
vector< vector<CC> > compcap;
vector< vector<long> > compcapdiv;

/*
vector< vector<long> > compcap;
vector< vector<long> > compcapty;
vector< vector<double> > compcapt;
*/

vector< vector<long> > capindob;                    // give the observation number for a given individual in a given capture

long npopm;                                         // Population measurements
vector<double> popmt;                                // The times of the measurements
vector<double> popmal;                               // The measurements are assumed to be gamma distributed with an alpha and beta
vector<double> popmbe;
vector< vector<long> > popmcomp;                    // The compartments that make up the population
vector <long> popmap;                               // Used to make changes in the observation probability for the population measurements
vector <long> poplist;

long nderm;                                         // Derived quantity measurements
vector<long> derm;                                  // Which quantity
vector<double> dermt;                               // The times of the measurements
vector<double> dermav;                              // The measurements are assumed to be gamma distributed with an alpha and beta
vector<double> dermvar;
vector< vector<long> > dermcomp;                    // The compartments which affect derivation
vector <long> dermap;                               // Used to make changes in the observation probability for the derived measurements
vector <long> derlist;

long NOTALIVE;                                      // Represents the not alive compartment

long nparam;                                        // Parameters in the model
vector<string> paramname;                           // The name of the parameter
vector<string> paramclass;                          // A reference names used in BIP to classify traces
vector<long> paramprior;                            // The prior used for the prior
vector<long> nparampriordep;                        // the priors which depend on that parameter
vector< vector<long> > parampriordep;

long nderive;
vector<string> derivename;
vector<long> derive;
vector< vector< long> > nderive_dpop;
vector< vector< vector<long> > > derive_dpop;       // Works out how populations change as compartments change
vector< vector< vector<double> > > derive_dpopweight;

long nsmooth;
vector<string> smoothname;
vector<long> smoothtype;
vector<string> smoothref;
vector<long> nsmoothparam;
vector< vector<double> > smoothval;
vector< vector<long> > smoothparam;

vector< vector<long> > paramsmooth;
vector< vector<long> > paramsmoothi;

long nderivetemporal;                               // Lists all erived quantities that are temporal
vector <long> derivetemporal;

vector<long> nparam_dep;                            // These are use to keep track of what changes a parameter makes (used in parameter proposals)
vector< vector<long> > param_dep;

vector<long> nparam_notdep;
vector< vector<long> > param_notdep;

vector<long> nparam_cope;
vector< vector<long> > param_cope;

vector<long> nparam_cpe;
vector< vector<long> > param_cpe;

vector<long> nparam_capev;
vector< vector<long> > param_capev;

vector<long> nparam_derm;
vector< vector<long> > param_derm;

vector<long> param_nmfl;                            // Determines if a change in a parameter requires the non-Markovian likelihood to be updated

vector <double> firstobst, lastobst;                 // The first and last time an  individual is observed

// Uses in MCMC proposals

double problifeprop[7] = {0.25,0.25,0.1,0.1,0.1,0.1,0.1};
double problifepropsum[7];

long nclnm;
vector <long> clnm;
vector <long> cstnm;
vector <double> tstnm;

struct EV { double t; long tr;};
bool compareev(EV lhs, EV rhs) { return lhs.t < rhs.t; }

struct PART { vector <EV> ev; long cend; double w; double wsum; long pback; long cstart; vector <long> futerefst;};
vector<PART> part;

vector<long> neq_popnum;                            // Gives a list of all the populations involved in an equation
vector< vector<long> > eq_popnum; 
vector<long> neq_param;                             // Gives a list of all the parameters involved in an equation
vector< vector<long> > eq_param;

long npopnum;                                       // The number of populations which are kept track of
vector<string> popnumname;                          // Name of population
vector< vector<long> > popnumterm;                  // The compartment populations which are added 
vector< vector<float> > popnumtermweight;           // The weightings given to compartments
long popnummax;                                     // The maximum number of populatins an equation can depend on

long nnmeq;                                         // List of equation used in non-Markovian 
vector<long> nmeq;

vector<long> transdep;                              // [eqn]   -1 = not an exponential transitions   0 = does not depend on population number   1 = depends on population number
vector<long> transdepref;                           // [eqn]  references which transdepeq
long ntransdepeq, ntransnotdepeq;
vector<long> transdepeq;                            // a list of all the equations dependent on population
vector<long> transnotdepeq;                         // a list of all the equations independent of population

vector< vector <long> > transdepeqrecalc;           // Determines if an equation needs to be recalculated based on removing the exiting individual
vector< vector <double> > popnumtake;

vector< vector< vector<long> > > transchref;

long ist; // temp
double Linmst;
long nwrong;

struct EQCH { long d; long n; long valch; double* popnum;};


vector< vector<EQCH> > transdepeqch;

struct NDEQCH { long d; long n;};

vector< vector<NDEQCH> > transnotdepeqch;

struct TRANS{ long cl; long i; long f;};           // Store all unique transitions in move (used for jump move)
long ntrans;
vector<TRANS> trans;

// Used for initproc


/*
vector <double> alinit;
long nindinit;
vector<long> ncapevdiv;
vector< vector<double> > capevdivt;
vector< vector<long> > capevdivnum;
*/

// Used for particles

const long CHC = 0, BEG = 2, END = 3, OBS = 4, FIXEV = 5, OTHCL = 6;
struct EVPART { double t; long ty; long ch;};

long nevp;                                                      // when doing simulation from particles this shows all the fixed times
vector <EVPART> evp;
vector<long> evpref;                                            // Gives a reference to cexist

long irev;                                                       // Used to keep track of event on the exisiting sequence
long erev, nevrev;
vector <EV> evrev;

vector <long> psel;

double wsum;
vector <double> wsumst;

long pbegst, pendst;                                            // Temporarily stores variables

long nexist;
vector<long> cexist;                                            // Keeps track of changes to c in the existing event sequence
vector<double> texist;

long nRdep, nRnotdep;                                           // Used when doing simulations
double Rdep,  Rnotdep;
vector <double>  Rdepst, Rnotdepst, rdep, rnotdep;

bool evpcomp(EVPART lhs, EVPART rhs) { return lhs.t < rhs.t; }

long nfuteref;
vector <long> futeref;
//vector <long> evc;
//vector <double> evt;

long npart;
long partcl;  // Determines if transitions are performed along a particular classification or all
long clall;   // The classification used to represent all

struct EVPSEC{ long evpbeg; long evpend; double tbeg; double tend;};

long nevpsec;
vector<EVPSEC> evpsec ;                                         // Divides up into sections which are simulated

//double part_sum;                                                 // Used for sampling the initial compartment for the particle
//vector<double> part_sumst;

//double wsum_part;

double partdt;                                                   // The Minimum time step used when doing simulation

double tent, tlea;                                               // stores the entry and leaving time for an individual

long ctop;                                                       // This represents age time and fix part of the compartment on individual entry

vector <vector <vector <vector< long> > > > fixnmmap;

// Used in adding and removing individuals

const long naddremdiv = 10;                                   // The number of time divisions used when inserting individual

vector <double> induoagesumtot;
vector <vector <double> > loginduoage, induoage;                                      // The age distribution for added individuals
vector <vector <double> > induoagesum;                                      // The age distribution for added individuals

double addremdivsumtot;
vector <double> logaddremdiv, addremdiv;
vector <double> addremdivsum;

vector <double> induolifetime;                                      // The expected lifetime for added individuals

// temporary variables used
double logsum[logsummax];
string capwarn;
double uocapmin, uocapmax;
long n_popnum;
long sampill;
//long tmp;
vector <EV> evnew;
//double temp;

short notimerange = 0;

double pmean, pvar, pgrad;
vector<double> pvel;

class Chain                           // all the quantities for each chain
{
  public:
    double Lir;
    double Liexp;
    double Linm;
    double Liinit;
    double Lev_pd;
    double Lob_pd;
    double Lob_st;
    double Lpop;      // The likelihood for the population measurements
    double Lder;      // The likelihood for the derived measurements
    double Lpri;

    vector <long> pop;                // current population at an observe population point
    vector <double> popL;             // the likelihood of this

    vector<double*> dermpopnum;       // Populations used to calculate derived values
    vector <double> derL;             // the likelihood of derived values

    double invT_pd;
    double invT_st;
    double invTLi;
    double invT_pop;
    double invT_der;

    //vector <double> obs_part;
    //vector <double> fix_part;
    vector< vector <double> > nindpart;            // The number of particles used for an individual

    long nindtot;                                  // Information about the individuals
    vector <double> indtbirth;
    vector<long> nindev;
    vector <vector <EV> > indev;

    double *param;

    long nderivetemporalplot;
    vector< vector< vector<double> > > derivetemporalplot;
 
    vector< vector <double> > pairac;   // Used to optimise zero two pairs

    vector <double> ntr_smooth, nac_smooth, jump_smooth;

    vector <double> ntr_param, nac_param, jump_param;
    vector <double> ntr_paramsamp, nac_paramsamp;
    vector <double> ntr_paramnorm, nac_paramnorm;
    vector < vector <double> > ntr_part, nac_part, nfa_part, nde_part;

    vector <double> ntr_indsim, nac_indsim, nfa_indsim;
    double ntr_indsimuo, nac_indsimuo, nfa_indsimuo;

    vector <double> ntr_tbirth, nac_tbirth, jump_tbirth;    // Diagnostics for proposals on life of observed individuals
    vector <double> ntr_tbirthent, nac_tbirthent, jump_tbirthent;
    vector <double> ntr_tent, nac_tent, jump_tent;
    vector <double> ntr_tlea, nac_tlea, jump_tlea;
    vector <double> ntr_entswitch, nac_entswitch;
    vector <double> ntr_birthentswitch, nac_birthentswitch;
    vector <double> ntr_leaswitch, nac_leaswitch;

    double ntr_tbirthuo, nac_tbirthuo, jump_tbirthuo;    // Diagnostics for proposals on life of unobserved individuals
    double ntr_tbirthentuo, nac_tbirthentuo, jump_tbirthentuo;
    double ntr_tentuo, nac_tentuo, jump_tentuo;
    double ntr_tleauo, nac_tleauo, jump_tleauo;
    double ntr_entswitchuo, nac_entswitchuo;
    double ntr_birthentswitchuo, nac_birthentswitchuo;
    double ntr_leaswitchuo, nac_leaswitchuo;

    double ntr_add, nac_add, nfa_add;
    double ntr_rem, nac_rem;

    vector< vector <double> > ntr_move, nac_move;
    vector <vector <double> > jump_move;

    vector< vector <double> > ntr_sing, nac_sing, nfa_sing;
    vector< vector <double> > ntr_twothree, nac_twothree, nfa_twothree;
    vector< vector <double> > ntr_pair, nac_pair, nfa_pair;
    vector< vector <double> > ntr_gap, nac_gap, nfa_gap;

    long nindinit;                      // The number of initial individuals
    vector <long> ncompinit;            // stores the number of individuals entering the poulation in the various compartments
    vector <double> probinit;           // The probability of an initial individual having a given compartment

    vector< vector<long> > depeqdiv;     // hash table in time

    long ndepeq_ev;                      // these give an unorder list describing the time variation in the dependent equations
    vector<long> depeq_evn;             // the number of times this equation appears in R
    vector< double* > depeq_evpopnum;    // the population numbers for that particular
    vector<double> depeq_evt;            // The time of the change in the 
    vector<long> depeq_evnext;           // The next member on the list
    vector<long> depeq_evback;           // The last member on the list
    vector<long> ndepeq_ev_ev;          // The cases where the actual event happens
    vector< vector<double> > depeq_ev_evt;
    vector<double> depeq_evval;          // The value for the equation

    long ndepeq_evstack;
    vector<long> depeq_evstack;

    vector<long> transnotdepeq_num;
    vector<double> transnotdepeq_dt;
    vector<double> transnotdepeq_val;

    vector<double> nmeq_val;

    vector<long> cope_num;
    vector<double> cope_val;
    vector<double> cope_op;

    vector<long> cpe_num;
    vector<long> cpe_oneminus;
    vector<double> cpe_pd;
    vector<double> cpe_val;
    vector<double> cpe_valoneminus;

    vector<double> capev_val;                           // The equation for the capture event probability
    vector<double> capev_valoneminus;                      // The equation for the capture event probability
    vector<long> capev_num;                             // The number of events observed
    vector<long> capev_oneminus;                          // The number of events not observed

    vector< vector <double> > simcl_initcsumtot;
    vector< vector< vector <double> > > logsimcl_initc, simcl_initc, simcl_initcsum;     // Used for sampling the initial compartment within a classification

    vector< vector <double> > simcl_entercsumtot;
    vector< vector< vector <double> > > logsimcl_enterc, simcl_enterc, simcl_entercsum;  // Used for sampling the enter compartment within a classification

    vector <double> sim_initcsumtot;
    vector< vector <double> > logsim_initc, sim_initc, sim_initcsum;           // Used for sampling the initial compartment when simulating
    vector <double> sim_entercsumtot;
    vector< vector <double> > logsim_enterc,sim_enterc, sim_entercsum;         // Used for sampling the enter compartment when simulating
    vector <double> simuo_initcsumtot;
    vector< vector <double> > logsimuo_initc, simuo_initc, simuo_initcsum;           // Used for sampling the initial compartment when simulating for the unobserved individuals
    vector <double> simuo_entercsumtot;
    vector< vector <double> > logsimuo_enterc, simuo_enterc, simuo_entercsum;         // Used for sampling the birth compartment when simulating for the unobserved individuals

    double part_prob, induo_prob, move_prob, sing_prob, twothree_prob, pair_prob, gap_prob, life_prob, addrem_prob;

  public:
  void initchain();
  long start();                                     // Works out the starting event sequence
  void addemptyind(double tbirth);                   // Adds a new empty individual (i.e. it is in the notalive category)
  void rememptyind(long i);                          // Removes an individual
  double priorprob(long pr);                        // Calculates the prior probability
  double priorsmooth(long s, long j);               // Calculates the smoothing prior
  double priorcalc();                                // The total prior probability
  void indchange(long i, long estart, double t, double tend);                           // Makes a change to the event sequence on individual i
  void indcha(long i);
  void indrev(long i);
  long eventmove(long i, long e, double tnew);       // Does fast movement of events
  void secchange(long i, long ci, long cf, double ti, double tf);  // Changes the state of an individual for a particular region
  void eqtimelineinit();                             // Initialises the timelines use for the dependent and not dependent equations
  long insertevent(long d, long e, double t, long eq, long dn, double *dpopnum, long sign);     // inserts event onto equation timeline
  void trydeleteevent(long d, long e);               // Tries to delete event on the equation time line
  void check(long num);                              // Used to chech the algorithm is working
  void checkparamgrad();                             // Checks that paramgrad is working correctly
  void checklike(long num);                          // Checks that the likelihood is correctly calculated
  void checkderive();                                // Checks that the derived likelihood is correct
  double calcobsprob(long eq);                       // Calculates the observation probaability
  void indremevent(long i, long e);                  // Removes an event for an indivual
  void indaddevent(long i, long k, EV ev);           // Adds an event for an indivual
  void chainobsinit();                               // Initialises quantities relating to the observation probability
  void param_prop();                                 // Makes proposals on model parameters
  void paramnorm_prop();                             // Makes a proposal which is normal fitted to posterio
  void paramgrad(long p);
  void param_gradprop();                                        // Makes proposals to increase posterior probability
  void param_sampprop();                                        // Makes proposals on model parameters by sampling from the prior
  void changeparam(long dir, long p, double param_new);       // Makes changes to parameters
  void diagnosticschain();                           // Outputs the success of MCMC proposals
  void part_prop(long i);                            // Makes a particle proposal on an individual
  void insertexisting(vector<EV> &ev);               // Populates particles with the existing sequence
  void addexistingfute(PART &pa, long tr, double tbeg, double tend);

  void getdep(vector <long> &vec, long co, double t);// Gets the total rate of dependent transitions
  void getnotdep(vector <long> &vec);               // Gets the total rate of not dependent transition
  void chainpartinit();                              // Initialises quantities used when sampling from partivles
  double transra(long eq, double t, long cnow);    // Calculates the transition rate but removes the individual at cnow

  double transradep(long eq, double t, long cnow);    // Calculates the transition rate but removes the individual at cnow
  void move_prop(long i);                            // moves the events in individual i
  void addfuture(long i, long c,long cl, double tstart);    // When simulating adds future events 
  void addfuturetra(long i, long tr, double tstart, double tadd);
  void sim(double tmin);                             // simulates from the model starting at time tmin
  void sim_indchange(long i, long ci, long cf);   // When an individual changes in the simulation
  void eventplot();                                  // Outputs an events sample
  void sim_check(double R, double t);                // Checks that the simulation is being performed correctly
  double likenm(long i);                             // The likelikelihood for the non Markovian transitions in individual i
  //double linknmtra(long i, long tr, long est);      // Calculates the non-Markovian probability for transition tr;
  //void addfute(long tr, double t);                   // Adds a future event
  //double addfuterev(long tr, double t);                // Adds the future event consistent with exisiting event sequence
  void part_simsec(long i, long evpbeg, long evpend, PART &pa);    // simulates new sequence
 // double part_revsec(long i, long evpbeg, long evpend, PART &pa);    // probability of simulating the existing event sequence
  void part_sim(long i);                          // Simulates using particles
  //double simple_sim(long i);
  double indsimulate(long i);                                   // Simulates event sequence for an individual
  double indsim(long i, long c);     // simulates a new sequence for an individual
  double probindsimulate(long i);                                   // Probability event sequence for an individual
  double probindsim(vector<EV> &ev);     // probability of simulating a new sequence for an individual
  long selectstart(long i, double t);                            // Selects the starting compartment for a simulation
  double probstart(long i, double t, long c);              // The probability of the starting compartment in a simulation
  void indsim_init(long i, double tbirth, long incobs);                                      // Initialises dividsion used in simualation
  
  double part_simprob(long i, long evpbeg, long evpend, PART &pa);
  long getlifespan(long i);                          // Gets the time over which an individual is alive

  //double part_rev(long i);                          // Calculates the probability of create the exisiting event sequence
  void setnmeq();                                    // Calculates all the equations needed for the likelihood of non-Markovian events
  void update();                                     // Performs an MCMC update
  void traceplot();                                  // Plots traces for the parameter values
  void nmchange(long tr, double dt, double dt2);     // Calculates the chnage in Linm from moving an evvent
  void nmupchange(long tr, double dt, double dt2);   // Calculates the chnage in Linm from moving an evvent

  double L();                                        // The posterior probability
  double Lwnm();                                     // The posterior probability without the non-Markovian part
  //long inside(long cl, double tsta, double tbeg, long pback);   // Checks if there has been some intervening event in that classification
  void life_prop(long i);                           // Makes changes to entry and leaving times as well as time of birth
  long lifechange(long i, double tentnew, double tleanre, double tbirthnew, double dprob);    // Makes a change to the life
  void addreminit();                                 // Initialising proposals for adding and removing individuals
  void addrem_prop();                                // Adds and remove indiduals from the population
  void simsumcalc();                                 // Used to initialise sampling of initial compartment
  double Lout();                                     // Outputs elements of the likelihood
  void plotparamlike(long p, double min, double max);// Plots the marginal likelihood distribution for a parameter
  long numwrong(long i, vector<EV> &vec);            // Gives the number of wrong observations
  void sing_prop(long i);
  void twothree_prop(long i);                        // Interchanges between two and three events
  void pair_prop(long i);
  void gap_prop(long i);                             // Resample event sequences between the gaps in obsrvations
  void initparamsamp();                              // Randomly samples parameter values
  void initproc();                                   // Finds a reasonable set of initial conditions
  double fit();                                      // An approximate fit function to the data

  long localresamp(long i, long cl, long estart, double t, double tend, long fixend, long evfac);
  void indsim_prop(long i);                          // Does a proposal which simple resimulates the trajectory for the individual
  double priorsamp(long pr);

  double fixnmtrans(vector<EV> &ev);                               // Accounts for the situation in which there are fixed non-Markovian transitions
  void addfixnmfute(long tr, double t, vector<EV> &ev, long e);
  long addevnew(long tr, double t);
  double startseq(long i);
  void nmcha(vector <EV> &ev, long e, long dir);
  void deriveplotinit();                            // Initialises plots for derived values
  void deriveplotcalc();                            // Calculates data for plots of derived values
  void derivepl();                                  // Plots derived values c
  void optimiseproposals();                            // optimises the frequency of different proposals
  double probuocap(long i, vector<EV> &vec);          // Gets the probability of captures which are not observed
};

const long chainmax = 1;
Chain* ch[chainmax];
long nchain;

vector <string> warning;                             // Keeps a list of any warning messages

// used in equation.h
vector<string> getcommasep(string a);                // Takes a string and splits it up by the commas
long getclassval(long cl, string st);               // Returns the value of a classification from a string
vector<long> getallcomp(string st);                 // Gets all the comprtments from a string

// used in readinput.h
string get(XMLNode* node, string attr);              // This gets an attribute
double getnum(XMLNode* node, string attr);           // Gets a number fron an XML attribute
void getalltrans(XMLNode* node, XMLNode* node2, long cl, long in, string rate);   // works out all the compartments consistent with a transition
long exist(XMLNode* node, string attr);             // Determines if an attribute exists
long getint(XMLNode* node, string attr);            // Gets an integer XML attribute
vector<string> getlinesep(string a);                 // Seperates using the | mark
string repla(string st, string sub1, string sub2);   // Replaces one substring with another
long getparam(string na);                           // Gets the parametere number fron the name

// using in output.h
void outputmodel();                                  // Outputs the model (used in debugging)
void outcomp(long c);                               // Outputs the state of a compartment
void traceinit();                                    // Initialises the trace plot

// udes in prior
void priorinit();                                    // Initialises the prior

// used in likelihoodinit.h
void likelihoodinit();                               // Initialises the quantities need to make fast likelihood calculations

// used in ratecalc.h
double ratecalc(long eq, double *popnum, double *param);      // Calculates the rate of an event
double ratecalcdep(long d, double *popnum, double *param);    // uses ncompind and popnum from the reduced set for a link
double ratecalcnotdep(long eq, double *paramval);    // uses ncompind and popnum from the reduced set for a link

// used in observation.h
void observationinit();                              // Initialises the observation model
double popcalc(long p, double val);                 // Calculates the observation probability for a population
double dercalc(short j, double *popnum, double *param);    // Calculates the observation probability for a derived quantity

// param_prop.h
void paraminit();                                    // Initialisation for parameter proposals
long notinsideprior(long p, double val);             // Determines if a value is within the prior

// part_prop.h
void partinit();                                     // Initialises quantities used when sampling from partivles

// ratecalc.h
double ratecalcdeptakeoff(long d, double *popnum, double *param, long c);    // Removes the individual at c when the quantity is calculated

// sim,h
void siminit();                                       // Initialises variables for use in simuations

// dist.h
double gammaprob(double x, double a, double b);       // The log of the probability from the gamma distribution
double dgammaprob(double x, double x2, double a, double b);   // The difference in the log of the probability from the gamma distribution

double gammaup(double x, double a, double b);         // The log of the integral of the gamma function from x up to infinity
double gammasamp(double a, double b);                 // Draws a sample from x^(a-1)*exp(-b*x)
void gammainit();                                     // Initialises the gamma distribution

double weibullprob(double x, double lam, double kk);  // The log of the probability of sampling from the Weibull distribution
double weibullup(double x, double lam, double kk);    // The log of the integral of the Weibull function from x up to infinity
double weibullsamp(double lam, double kk);            // Samples from the Weibull distribution
double dweibullprob(double x, double x2, double lam, double kk);   // The differnce log of the probability of sampling from the Weibull distribution
double dweibullup(double x, double x2, double lam, double kk);     // The difference in the log of the integral of the Weibull function from x up to infinity

// startseq.h
void startseqinit();                                  // Initialises quantitites used in setingup the initial sequence
long getctop(double t, double tbirth);                // Returns the compartment part related to the top three classifications
long checkforfixedinpath(long c, long cc, long ctopp);// Considers if fixed events stop a path from occuring

// check.h
void checkevseq(vector<EV> &vec);                     // CHecks that an event sequnece is consistent and correct
void checkevseqsimp(vector<EV> &vec);

// common function
double ran(){ double v = (double(rand())*29999+rand())/(30000.0*RAND_MAX); if(v == 0 || v == 1) return 0.1; else return v;}

double normal(float mu, float sd){ return mu + sd*sqrt(-2*log(ran()))*cos(2*3.141592654*ran());}
long len(string s){ return s.length();}

void emsg(string msg)
{
  if(checkon == 1){
    cout << "e|ERROR MSG: " << msg << " ";
    if(errorcode >= 0) cout << "EC:" << errorcode ;
    cout << "    samp:"<< samp <<  "\n";
  }
  else{
    cout << "e|" << msg << "\n";
  }
  exit (EXIT_FAILURE); 
}

void emsg(string msg, double v1){ stringstream ss; ss << msg << ":" << v1; emsg(ss.str());}
void emsg(string msg, double v1, double v2){ stringstream ss; ss << msg << ":" << v1 << "," << v2; emsg(ss.str());}
void emsg(string msg, double v1, double v2, double v3){ stringstream ss; ss << msg << ":" << v1 << "," << v2 << "," << v3; emsg(ss.str());}

void emsg2(string msg)
{
  warning.push_back(msg);
}

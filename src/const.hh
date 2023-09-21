/// Stores constant values used and also enumarated types

#ifndef BICI__CONST_HH
#define BICI__CONST_HH

#include <string>

using namespace std;

//#define WINDOWS

const bool testing = true;//false;                           // Set to true for additional testing functions

#ifdef WINDOWS
const bool com_op = true;                            // Set to true for command line output
#else
const bool com_op = false; 
#endif

/************************** Enumerated values ******************************/

// Enumerates different types of algorithm for simulation / inference
enum Algorithm { GILLESPIE, TAU, DA_MCMC, MFA, ABC_ALG, ABC_SMC_ALG, ABC_MBP, ABC_PAS, PMCMC, HMC };

// Enumerates different feature type in geoJSON files
enum FeatureType { POLYGON, MULTI_POLYGON, NO_FEATURE};

// Enumerates the different types of equation
enum EqnType { SE, SP, TRAP_PROB, COMP_PROB, BP, SOURCE_RATE, TRANS_MEAN, TRANS_RATE, TRANS_SHAPE, TRANS_SCALE, TRANS_CV, REPARAM, DIST, DERIVE_PARAM, DERIVE_EQN };

// Enumerates restrictions places on equarions
enum EqnMode { PARAM_ONLY, PARAM_WITH_DEP, ALL, DERIVE_PARAM_MODE, DERIVE_MODE };

// Different types of species model
enum SpeciesType { POPULATION, INDIVIDUAL };  

// Different types of transition distribution
enum TransType { EXP_RATE, EXP_MEAN, GAMMA, ERLANG, LOG_NORMAL, WEIBULL };

// Different types of markov transition distribution
//enum MarkovTransType { MARKOV_EXP_RATE, MARKOV_EXP_MEAN, MARKOV_UNSET}; 
	
// Different varieties of transition
enum TransVariety { NORMAL, SOURCE_TRANS, SINK_TRANS };

// Different modes of operation
enum Operation { SIM, INF, MODE_UNSET };

// Different prior possibilities
enum PriorPos { UNIFORM_PR, EXP_PR, NORMAL_PR, GAMMA_PR, LOG_NORMAL_PR, BETA_PR, BERNOULLI_PR, FIX_PR, FLAT_PR, DIRICHLET_PR };

// Ennumerates all possible command types
enum Command {SPECIES, CLASS, SET, CAMERA, COMP, TRANS, SOURCE_CD, SINK_CD, CLONE, DATA_DIR, DESC, LABEL, BOX, PARAM, DERIVED, IND_EFFECT, FIXED_EFFECT, INIT_POP, ADD_IND, REMOVE_IND, MOVE_IND, INIT_POP_SIM, ADD_IND_SIM, REMOVE_IND_SIM, MOVE_IND_SIM, INIT_POP_PRIOR, COMP_DATA, TRANS_DATA, SOURCE_DATA, SINK_DATA, TEST_DATA, POP_DATA, POP_TRANS_DATA, SET_TRAPS_DATA, IND_TRAPPED_DATA, GENETIC_DATA, SIMULATION, INFERENCE, SIM_PARAM, SIM_STATE, INF_PARAM, INF_STATE, 

// These are not commands but varient used when loading data
TRANS_TIMERANGE_DATA, SOURCE_TIMERANGE_DATA, SINK_TIMERANGE_DATA,
EMPTY };

// Different types of time range used in transition data specification
enum TimeRange { ALL_TIME, SPEC_TIME, FILE_TIME };

// Types of filtering which can be done on a classification
enum FilterClaType { ALL_FILT, FILE_FILT, COMP_FILT};
	
// Types of observation model done in data
enum ObsModelType { PERCENT_OBSMOD, SD_OBSMOD, FILE_OBSMOD };

enum Coord { CARTESIAN, LATLNG };

// Determines the result of a function
enum Result { SUCCESS, FAIL};

// The different types of column which can be loaded from data files
enum ColumnType { ID_COL, T_COL, SNP_COL, CL_PROB_COL, CL_ALL_COL, FROM_COL, TO_COL, POP_COL, ALPHA_COL, START_COL, END_COL, RESULT_COL, FILT_OBSPOP_COL, FILT_OBSPOPTRANS_COL, FILT_SETTRAPS_COL, SETTRAP_NAME_COL,COMP_NAME_COL, START_COMP_COL, END_COMP_COL};

// The different types of element in a table (this help to spot errors in data tables)
enum ElementType { TEXT_EL, TIME_FLOAT_EL, FLOAT_EL, COMP_PROB_EL, COMP_EL, POS_FLOAT_PERCENT_EL, POS_INT_EL, POS_FLOAT_EL, MULTI_COMP_EL, POS_ZERO_FLOAT_EL, COMP_SOURCE_SINK_EL, COMP_ALL_EL, DIAG_TEST_EL};
	
// Different types of warning	
enum ErrorType { ERROR_FATAL, ERROR_WARNING };

// Different types of smoothing splines
enum SmoothType { NORMAL_SMOOTH, LOG_NORMAL_SMOOTH };

// Different varieties of parameter
enum ParamVariety { PRIOR_PARAM, REPARAM_PARAM, DIST_PARAM, CONST_PARAM, UNSET_PARAM };

// Different types of text
enum TextType { NORMAL_TEXT, SUB_TEXT, SUP_TEXT };

// Different types of initial condition
enum InitCondType { POP_INIT, ZEROPOP_INIT, UNSET_INIT };

// Different types of event
enum EventType { ENTER_EV, LEAVE_EV, MOVE_EV, NM_TRANS_EV, M_TRANS_EV};

// Different types of observation
enum ObsType { OBS_TRANS_EV, OBS_COMP_EV, OBS_TEST_EV, OBS_POP};
	

// Different ways a parameter can affect likelihoods
enum AffectType { SPLINE_PRIOR_AFFECT, PRIOR_AFFECT, DIST_AFFECT, SPLINE_AFFECT, EXP_FE_AFFECT, DIV_VALUE_AFFECT, MARKOV_LIKE_AFFECT, POP_AFFECT, NM_TRANS_AFFECT, OMEGA_AFFECT, EXP_IE_AFFECT, LIKE_IE_AFFECT, INDFAC_INT_AFFECT, MARKOV_POP_AFFECT, LIKE_OBS_AFFECT};
	
/// Determines how affect proposals are acted on
enum AllTimeStatus { ALL_DIV, SEL_DIV};

/// Different proposal types
enum PropType { PARAM_PROP, IND_PROP, INDSAMP_PROP, MBP_PROP, MBPII_PROP, IE_PROP, IE_VAR_PROP};

/// Different types of MBP 
enum MBPType {TYPE_I, TYPE_II};

/// Different types of individual change
//enum IndChaType {REMOVE_IND_TRANS, CHANGE_LI_MARK_IND_TRANS};

/************************** Numeric constants ******************************/

const auto ERR_MSG_MAX = 10u;                     // The maximum number of error messages

const auto CODE = 99999989u;                      // THe threshold below which not code
const auto SOURCE = 99999990u;
const auto SINK = 99999991u;

const auto UNSET = 99999992u;
const auto TIME_VAR = 99999993u;

const auto IND_SAMP_THRESH = 0.9;                 // The threshold for individual event sampler

const auto LI_WRONG = -10;                        // Likelihood penalty for wrong (e.g. transition)

const double TINY = 0.00000001;                   // Used to represent a tiny number
const double TWO_TINY = 2*TINY;                   // Used to represent a tiny number
const double SMALL = 0.00001;                     // Used to represent a small number
const double LARGE = 1000000000;                  // Used to represent a big number
const double ALMOST_ONE = 0.99999999999999;          // Almost one
const double OVER_ONE = 1.0000000000001;
const double LOG_THRESH = 0.00001;                // The threshold below which logs not calculated

# define M_PI  3.14159265358979323846

const unsigned int NEXT_IND_EVENT_DIV = 100;      // Number of divisions next_event divided into

/************************** Text constants ******************************/

const vector <string> alphabet = {"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"};

const string chnotallowed = "\"!%$&=;@~#\\?";       // Character not allowed in equations
const string notparam_list = "+-*×/0123456789.{<>}() %$";   // Sets not a parameter character
const string paramend_list = "+-*×/0123456789.{<>}[] %$";	    // Sets characters which end parameter

const string err_mess = "Please correct these errors and rerun.";

// Different posibilities for transition type
//const vector <TransType> dist_pos = {  EXP_RATE, EXP_MEAN, GAMMA, ERLANG, LOG_NORMAL, WEIBULL };
//const vector <TransType> exp_dist_pos = {  EXP_RATE, EXP_MEAN, ERLANG}; 


/************************** Vector constants ******************************/
//const vector <PriorPos> prior_pos = {UNIFORM_PR,EXP_PR,NORMAL_PR,GAMMA_PR,LOG_NORMAL_PR,BETA_PR,BERNOULLI_PR,FLAT_PR,DIRICHLET_PR};
//const prior_pos_all = ["uniform","exp","normal","gamma","log-normal","beta","bernoulli","fix","flat","dirichlet"];
//const bp_prior_pos = ["flat","dirichlet"];
//const zeroone_prior_pos = ["flat","beta"];

/************************** Object constants ******************************/

struct EqnInfo {
	EqnInfo(EqnType ty, EqnMode mo, bool obs) { type = ty; mode = mo; obs_model = obs;}
	EqnType type;
	EqnMode mode;
	bool obs_model;
};

const vector <EqnInfo> eqn_types = { 
EqnInfo(SE,PARAM_ONLY,true),
EqnInfo(SP,PARAM_ONLY,true),
EqnInfo(TRAP_PROB,PARAM_ONLY,true),
EqnInfo(COMP_PROB,PARAM_ONLY,true),
EqnInfo(BP,PARAM_WITH_DEP,false),
EqnInfo(TRANS_MEAN,ALL,false),
EqnInfo(TRANS_RATE,ALL,false),
EqnInfo(TRANS_SHAPE,ALL,false),
EqnInfo(TRANS_SCALE,ALL,false),
EqnInfo(TRANS_CV,ALL,false),
EqnInfo(REPARAM,PARAM_ONLY,false),
EqnInfo(DIST,PARAM_ONLY,false),
//EqnInfo(DERIVE_PARAM,DERIVE_PARAM_MODE,false),
EqnInfo(DERIVE_EQN,DERIVE_MODE,false)
};

#endif
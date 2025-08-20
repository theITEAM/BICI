/// Stores constant values used and enumerated types

#pragma once

#include <string>
#include <cmath>

using namespace std;

// This is set if compilation is done on windows
#ifdef _WIN64
  #define WINDOWS
#elif __linux
  #define LINUX
#elif __APPLE__
  #define MAC
#else
# error "Unknown compiler"
#endif

#ifdef MAC
const string default_file = "/tmp/init.bici";        // This is used for Mac
#else
const string default_file = "Execute/init.bici";     // This is used for windows / linux
#endif

#define USE_MPI                                    // Sets if code can run in parallel

const bool jamie_code = true;                        // This turns in code for Jamie
const bool testing = true;                           // Set to true for additional testing
const bool slow_check = false;                       // Additional checks which are slow
const bool check_lin = true;                         // Checks linearisation of equations 
const bool adapt_prop_prob = false;//true;           // Determines if proposals prob adapted
const bool check_swap = false;                       // Used to check if new swap is working
const bool use_ind_key = true;                       // Determines if indindividual key
const bool print_diag_on = false;                    // Prints diagnostic statements to terminal
const bool cum_diag = true;                          // Cumulative probability diagnostic

const bool linearise_speedup = true;                 // Linearisation speed-up of likelihood
// This speeds up by taking account of linear population terms when calculation div values

const bool linearise_speedup2 = true;                // Linearisation speed-up likelihood
// This speeds up by considering changes to non-population terms when calculating div 

const bool linearise_speedup3 = true;                 // Linearisation speed-up of likelihood
// This speeds up update_ind by updating value based on changes in population (good for sums)

const bool removeparamvec_speedup = true;            // Removes zeros from paramter vector
const bool simplify_eqn = true;                      // Simplifies equations (to speed up)

/************************** Enumerated values ******************************/

// Information from tags when BICI is run
enum TagType { CHAIN, OP, SAMP_TYPE, SEED, TAG_UNSET };


// Different types of algorithm for simulation / inference
enum Algorithm { GILLESPIE, TAU, DA_MCMC, PAS_MCMC, MFA_ALG, ABC_ALG, ABC_SMC_ALG, ABC_MBP,  PMCMC, HMC };

// Different feature types in geoJSON files
enum FeatureType { POLYGON, MULTI_POLYGON, NO_FEATURE };

// Different types of equation
enum EqnType { SE, SP, COMP_PROB, TRANS_PROB, BP, SOURCE_RATE, SOURCE_MEAN, TRANS_RATE, TRANS_MEAN, TRANS_NM_RATE, TRANS_NM_MEAN, TRANS_SHAPE, TRANS_SCALE, TRANS_CV, REPARAM, REPARAM_EQN, DIST, DERIVE_PARAM, DERIVE_EQN, SEQ_VAR,CONST_EQN , MUT_RATE };

// Differetn restrictions placeed on equarions
enum EqnMode { PARAM_ONLY, PARAM_WITH_DEP, ALL, DERIVE_PARAM_MODE, DERIVE_MODE };

// Different types of species model
enum SpeciesType { POPULATION, INDIVIDUAL };  

// Determines if a compartment is infected or not
enum CompInfected { COMP_INFECTED, COMP_UNINFECTED, COMP_INFECTED_UNSET };

// Determines how infection status changes
enum TransInfectionType { TRANS_INFECTION, TRANS_RECOVERY, TRANS_INF_UNCHANGE };

// Different types of individuals
enum IndType { UNOBSERVED_IND, OBSERVED_IND };

// Different types of transition distribution
enum TransType { EXP_RATE, EXP_MEAN, EXP_RATE_NM, EXP_MEAN_NM, GAMMA, ERLANG, LOG_NORMAL, WEIBULL, PERIOD };

// Different ways individual factor can be changed
enum IndFacRateType { IND_EFF_MULT, IND_EFF_DIV, FIX_EFF_MULT, FIX_EFF_DIV };

// Different varieties of transition
enum TransVariety { NORMAL, SOURCE_TRANS, SINK_TRANS };

// Different modes of operation
enum Operation { SIM, INF, PPC, MODE_UNSET };

// Different prior possibilities
enum PriorPos { UNIFORM_PR, EXP_PR, NORMAL_PR, GAMMA_PR, LOG_NORMAL_PR, BETA_PR, BERNOULLI_PR, FIX_PR, DIRICHLET_PR, MDIR_PR, UNSET_PR };
 
// Different possible command types
enum Command {SPECIES, CLASS, SET, CAMERA, COMP, COMP_ALL, TRANS, TRANS_ALL, CLONE, DATA_DIR, DESC, LABEL, BOX, PARAM, DERIVED, IND_EFFECT, FIXED_EFFECT, INIT_POP, ADD_POP, REMOVE_POP, ADD_IND, REMOVE_IND, MOVE_IND, INIT_POP_SIM, ADD_POP_SIM, REMOVE_POP_SIM, ADD_IND_SIM, REMOVE_IND_SIM, MOVE_IND_SIM, ADD_POP_POST_SIM, REMOVE_POP_POST_SIM, ADD_IND_POST_SIM, REMOVE_IND_POST_SIM, MOVE_IND_POST_SIM, COMP_DATA, TRANS_DATA,  TEST_DATA, POP_DATA, POP_TRANS_DATA, IND_EFFECT_DATA, IND_GROUP_DATA, GENETIC_DATA, SIMULATION, INFERENCE, POST_SIM, SIM_PARAM, SIM_STATE, INF_PARAM, INF_PARAM_STATS, INF_STATE, POST_SIM_PARAM, POST_SIM_STATE, INF_DIAGNOSTICS, INF_GEN, MAP, PARAM_MULT, TRANS_DIAG, SIM_WARNING, INF_WARNING, PPC_WARNING,

// These are not commands but varient used when loading data
TRANS_TIMERANGE_DATA,
EMPTY };

// Different types of time range used in transition data specification
enum TimeRange { ALL_TIME, SPEC_TIME, FILE_TIME };

// Types of filtering which can be done on a classification
enum FilterClaType { FILE_FILT, COMP_FILT };

// Types of observation model done in data
enum ObsModelType { NORMAL_PERCENT_OBSMOD, NORMAL_SD_OBSMOD, NORMAL_FILE_OBSMOD, POISSON_OBSMOD, NEGBIN_OBSMOD, NEGBIN_FILE_OBSMOD };

enum ObsModelVariety { NORMAL_OBS, POISSON_OBS, NEGBIN_OBS};

// Different coordinate type
enum Coord { CARTESIAN, LATLNG };

// Determines the result of a function
enum Result { SUCCESS, FAIL };

// Different types of column which can be loaded from data files
enum ColumnType { ID_COL, T_COL, TSTART_COL, TEND_COL, SNP_COL, GENOBS_COL, CL_PROB_COL, CL_ALL_COL, CL_ALL_PROB_COL, FROM_COL, TO_COL, INIT_POP_COL, ADD_POP_COL, REM_POP_COL, START_COL, END_COL, RESULT_COL, FILT_OBSPOP_COL, FILT_OBSPOPTRANS_COL, COMP_NAME_COL, START_COMP_COL, END_COMP_COL };

// Different types of element in a table (this helps to spot errors in data tables)
enum ElementType { TEXT_EL, TIME_FLOAT_EL, FLOAT_EL, COMP_PROB_EL, COMP_PROB_NA_EL, COMP_EL, POS_FLOAT_PERCENT_EL, POS_INT_EL, POS_FLOAT_EL, ZERO_ONE_EL, POS_ZERO_FLOAT_EL, COMP_SOURCE_SINK_EL, COMP_ALL_EL, DIAG_TEST_EL, PRIOR_DIR_EL };

// Different types of warning
enum ErrorType { ERROR_FATAL, ERROR_WARNING };

// Different types of smoothing splines
enum SmoothType { NORMAL_SMOOTH, LOG_NORMAL_SMOOTH };

// Different varieties of parameter
enum ParamVariety { PRIOR_PARAM, REPARAM_PARAM, DIST_PARAM, CONST_PARAM, UNSET_PARAM };

// Different types of text
enum TextType { NORMAL_TEXT, SUB_TEXT, SUP_TEXT };

// Different types of event
enum EventType { ENTER_EV, LEAVE_EV, MOVE_EV, NM_TRANS_EV, M_TRANS_EV };

// Different types of trigger events
enum SimTrigEventType { ENTER_SIM_EV, LEAVE_SIM_EV, MOVE_SIM_EV, NM_TRANS_SIM_EV, M_TRANS_SIM_EV, DATA_TRANS_SIM_EV, SOURCE_SIM_EV };

// Differnt types of event when doing simulation
enum TrigEventType { TRIG_OBS_TRANS_EV, TRIG_MOVE_EV, TRIG_LEAVE_EV };

// Different types of observation
enum ObsType { OBS_SOURCE_EV, OBS_TRANS_EV, OBS_SINK_EV, OBS_COMP_EV, OBS_TEST_EV };

// Different ways a parameter can affect likelihoods
enum AffectType { SPLINE_PRIOR_AFFECT, PRIOR_AFFECT, DIST_AFFECT, SPLINE_AFFECT, EXP_FE_AFFECT, DIV_VALUE_AFFECT, DIV_VALUE_FAST_AFFECT, DIV_VALUE_LINEAR_AFFECT, MARKOV_LIKE_AFFECT, POP_AFFECT, NM_TRANS_AFFECT, NM_TRANS_BP_AFFECT, NM_TRANS_INCOMP_AFFECT, OMEGA_AFFECT, EXP_IE_AFFECT, LIKE_IE_AFFECT, INDFAC_INT_AFFECT, MARKOV_POP_AFFECT, LIKE_OBS_IND_AFFECT, LIKE_OBS_POP_AFFECT, LIKE_OBS_POP_TRANS_AFFECT, OBS_EQN_AFFECT, LIKE_UNOBS_TRANS_AFFECT, POP_DATA_CGL_TGL_AFFECT, LIKE_INIT_COND_AFFECT, PRIOR_INIT_COND_AFFECT, LIKE_GENETIC_PROCESS_AFFECT, GENETIC_VALUE_AFFECT, LIKE_GENETIC_OBS_AFFECT, IIF_W_AFFECT, POPNUM_IND_W_AFFECT, AFFECT_MAX };

// Different proposal types
enum PropType { PARAM_PROP, IND_EVENT_TIME_PROP, IND_MULTI_EVENT_PROP, IND_EVENT_ALL_PROP, IND_OBS_SAMP_PROP, IND_OBS_RESIM_PROP, IND_OBS_RESIM_SINGLE_PROP, IND_UNOBS_RESIM_PROP, IND_ADD_REM_PROP, IND_ADD_REM_TT_PROP, MBP_PROP, MBPII_PROP, MBP_IC_POP_PROP, MBP_IC_POPTOTAL_PROP, MBP_IC_RESAMP_PROP, INIT_COND_FRAC_PROP, IE_PROP, IE_VAR_PROP, IE_COVAR_PROP, IE_VAR_CV_PROP, TRANS_TREE_PROP,  TRANS_TREE_SWAP_INF_PROP, TRANS_TREE_MUT_PROP, TRANS_TREE_MUT_LOCAL_PROP, POP_ADD_REM_LOCAL_PROP, POP_MOVE_LOCAL_PROP, POP_IC_LOCAL_PROP, POP_END_LOCAL_PROP, POP_SINGLE_LOCAL_PROP, POP_IC_PROP, POP_IC_SWAP_PROP, PAR_EVENT_FORWARD_PROP, PAR_EVENT_FORWARD_SQ_PROP,PAR_EVENT_BACKWARD_SQ_PROP, IND_LOCAL_PROP, CORRECT_OBS_TRANS_PROP, IND_OBS_SWITCH_ENTER_SOURCE_PROP, IND_OBS_SWITCH_LEAVE_SINK_PROP };

// Different types of annealing stratregy
enum AnnealType { ANNEAL_NONE, ANNEAL_SCAN, ANNEAL_POWERAUTO, ANNEAL_LOGAUTO, ANNEAL_POWER };

// Different state timers
enum Timer { IND_TIMER, IND_POP_UPDATE_TIMER, UPDATE_SAMPLER_TIMER, CHECK_TIMER, TIMER_MAX };

// Different output timers
enum OutTimer { PARAM_OUTPUT, STATE_OUTPUT, OUTTIMER_MAX };

// Different check timers 
enum CheckTimer { CHECK_TRANS_NUM, CHECK_DEP_PARAM, CHECK_REF, CHECK_MARKOV, CHECK_NM, CHECK_LIKE, CHECK_SPLINE, CHECK_PRIOR, CHECK_POP, CHECK_POP2, CHECK_CPOP, CHECK_IE, CHECK_POP_LIKE, CHECK_MAPS, CHECK_EV_OBS, CHECK_OBS_LIKE, CHECK_IC, CHECK_LIN, CHECK_GEN,CHECK_POP_IND, CHECK_ADD_REM, CHECK_RANGE, CHECK_SIMP, CHECK_POP_IND_GENTIC, CHECK_MAX};

const vector <string> check_name = { "CHECK_TRANS_NUM", "CHECK_DEP_PARAM", "CHECK_REF", "CHECK_MARKOV", "CHECK_NM", "CHECK_LIKE", "CHECK_SPLINE", "CHECK_PRIOR", "CHECK_POP", "CHECK_POP2", "CHECK_CPOP", "CHECK_IE", "CHECK_POP_LIKE", "CHECK_MAPS", "CHECK_EV_OBS", "CHECK_OBS_LIKE", "CHECK_IC", "CHECK_LIN", "CHECK_GEN", "CHECK_POP_IND", "CHECK_ADD_REM", "CHECK_RANGE", "CHECK_SIMP", "CHECK_POP_IND_GENTIC", "CHECK_MAX"};

// Different proposal timers
enum PropTimer { PROP_TIMER, PROPTIMER_MAX };

// Different state species timers
enum StateSpeciesTimer { STSP_TIMER_MAX };

// INIT_POP_FIXED means init pop has been specified in the data file
// INIT_POP_DIST means a distribution has been specified in the data file
// INIT_POP_NONE means that no prior had been added to the initial population
// Only those individuals mentioned in the data enter
// (It doesn't mean there is no population)
enum InitPopType { INIT_POP_FIXED, INIT_POP_DIST, INIT_POP_NONE };

// Different types of result when changing initial conditions
enum ICResult { IC_SUCCESS, IC_FAIL, IC_ZERO };

// Determines how updates are done for non-markovian events
enum NmUpdate { UPDATE_NO_TIME, UPDATE_TIME, UPDATE_IND, NO_UPDATE };
 
// Different genetic characters
enum GenChar { A_CH, C_CH, T_CH, G_CH };
 
// Different infection node event type
enum InfEventType { INFECT_OTHER, GENETIC_OBS };
 
// Different types of genetic change
enum GenChaType { GENCHA_FAIL, NO_GENETIC_CHANGE, REGRAFT_NODE, ADD_NODE, REM_NODE };

// Determines how from and add are connected to the inf_node network
enum ConType { FROM_CON_ADD_CON, FROM_CON_ADD_UNCON, FROM_UNCON_ADD_CON, FROM_UNCON_ADD_UNCON };

// Different types of tree regrafting
enum RegraftType { NORMAL_REGRAFT, ORIGIN_REGRAFT, UNCHANGED_REGRAFT };

// Different local changes to a population model
enum LocalChangeType { REMOVE_EVENT, ADD_EVENT, REMOVE_IC, ADD_IC, ADD_C_IC, REMOVE_C_IC };

// Different type of non-Markovian update
enum NMUpdateType { NM_EVENT_STORE, NM_INCOMP_STORE };

// Determines is a branching probability is unset, set or derived
enum BpType { BP_UNSET, BP_SET, BP_DERIVED };
 
// Different directions of update for joint transition mean / event time proposals
// SQ refers to proposals which squash other transitions
enum Direction { FORWARD, FORWARD_SQ, BACKWARD_SQ };
 
// Different ways to restore state after an individual change
enum BackType { REMOVE_LI_MARKOV, ADD_LI_MARKOV, POP_DATA_NUM, LI_OBS_POP, POP_DATA_CGL, POP_TRANS_DATA_NUM, LI_OBS_POP_TRANS, LI_OBS_IND, TI_INDFAC, DLI_INDFAC, DIF_INDFAC, LI_MARKOV, OBS_TRANS_EQN_NUM, VALUE_MARKOV, IND_COND_VAL };

// Different ways restore population after an individual change
enum BackPopType { POP_NUM_T };

// Different ways of doing individual update
enum UpdateType { UP_SINGLE, UP_MULTI, UP_NOLIKE };

// Used in locally swapping events sequences for individuals
enum TrSwapType { ENTER_SW, SOURCE_SW, SINK_SW, TRANS_SW, LEAVE_SW, MOVE_SW };

// Used in locally swapping events sequences for individuals
enum EventLoc { DO_NOTHING_LOC, NEW_EVENT_LOC, OLD_EVENT_LOC, OBS_LOC };

// Used in locally swapping events sequences for individuals
enum CompRUpType { SINGLE_BRANCH_UP, SINGLE_CALCALL_UP, MARKOV_COMPR_UP, MARKOV_CALCALL_UP, NM_CALCALL_UP, NO_UP };

// Used in locally swapping events sequences for individuals
enum Modify { MULT_MOD, DIV_MOD };

// Different selections for how  individual updates are done (ALL_SAMP for everything)
enum SampType { LOCAL_SAMP, SAMP_SAMP, SIM_SAMP, SIM_CL_SAMP, ALL_SAMP };

// Different types of genetic data
enum GenDataType { SNP_DATA, MATRIX_DATA };

// Different ways in which local mutation are update on tree
enum GenBranchUpdate { OBS_GEN_UPDATE, SINGLE_BRANCH, TWO_BRANCH1, TWO_BRANCH2, TWO_BRANCH3, ORIGIN_GEN_UPDATE, GEN_UPDATE_MAX };
   
// Different modes when reading in state samples (for PPC)
enum InfStateMode { MODE_NODE, MODE_PARAM, MODE_PHYLO, MODE_DERIVE, MODE_INITIAL, MODE_POPCHANGE, MODE_TIMEPOINT, MODE_TRANS, MODE_INDIVIDUAL, MODE_TRANSDISTPROB };

// Sets direction of local event proposals
enum LocalDir { LOCAL_FORWARD, LOCAL_REVERSE };

// Sets different bound types on observation probabilities
enum BoundType { LOWER_BOUND, LOWER_UPPER_BOUND };

// Sets how markov equation is recalculated
enum RecalcMarkovType { RECALC, USE_POP_DIF, USE_POP_DIF_TIME};

// Sets diffent distribution for displaying text
enum DistText { NORM_TE, LOGNORM_TE, WEIBULL_TE, GAMMA_TE, BETA_TE, NEGBINO_TE, BERN_TE, EXP_RATE_TE, EXP_MEAN_TE, POIS_TE, PERIOD_TE};

// Different quantities
enum DistQuant { SD_QU, CV_QU, MEAN_QU, NORM_MEAN_QU, SHAPE_QU, SCALE_QU, ALPHA_QU, BETA_QU, P_QU, BERNP_QU, RATE_QU, EXP_MEAN_QU, POIS_QU, TIME_QU};

// Different percentage measures
enum PercentType { LOAD_PER, INIT_PER, RUN_PER, RUN_GEN_PER, ANNEAL_PER, OUTPUT_PER};

// Different types of derived function
enum DerFuncType { RN, RNE, RNC, GT, GTE, GTC, DF_UNSET};
 
/************************** Numeric constants ******************************/

const auto ERR_MSG_MAX = 5u;                     // The maximum number of error messages

const auto CODE = 99999989u;                      // Value above which is a non-numeric number
const auto UNSET = 99999990u;                     // Indicates an unset number
const auto USINT_MAX = 4294967295u;               // Largest number of unsigned int
const auto UNSET_LIST = USINT_MAX-1;              // Used in operator lists 
const auto CUT_LIST = USINT_MAX-2;                // Used to represent cut out of list
const auto OP_MAX = USINT_MAX-10;                 // Maximum number of operators
const auto UNSET_WILD = 99999980.0;               // Indicates an unset wildcard "*"
const auto UNSET_F = 99999990.0;                  // Floating point unset
const auto TIME_VAR = 99999991u;                  // Indicates time variable t (used in eqnations)
const auto OUTSIDE_INF = 99999992u;               // Stands for outside infection
const auto ENTER_INF = 99999993u;                 // Stands for entering infection
const auto BP_FROM_OTHERS = 99999994u;            // Branch probability is calculated from others
const auto GEN_PLOT = 99999995u;                  // Used for a generation plot
const auto LOG_HALF = log(0.5);

const auto LINK_CV_THRESH = 0.3;                  // CV threshold for linking for multi-event
const auto LINK_WEIBULL_SHAPE = 3.7;              // Shape thresh for linking for multi-event

const auto ESS_THRESH = 200;                      // Threshold for ESS statistic
const auto GR_THRESH = 1.01;                      // Threshold for GR statistic

const auto MEM_FRAC_MAX = 0.8;                    // Maximum fraction of memory usable
const auto MEM_FREE_MIN = 0.1;                    // Minimum fraction of memory available
const auto SIM_TRY_MAX = 10000u;                  // Maximum tries simulating nm trans
const auto TFORM_MAX = 100000u;                   // Maximum number for tform lookup table
const auto CGC_MAX = 100000u;                     // Maximum number for comp_global_convert
const auto EQ_POP_REF = 20u;                      // Store reference for population in eq
const auto IC_NUM_OPTIMUM = 100u;                 // Optimum number of events for IC local range
const auto RATE_RECOMMEND = 0.25;                 // Sets recommented wait (for warning)
const auto ADDREM_NUM_OPTIMUM = 5u;               // The optimum number for add/rem local
const auto WIN_UPDATE_ADAPT = 0.99;               // Rate of adaptation for win (local proposals)
const auto END_LOCAL_ADAPT = 0.99;                // Controls the rate of adaptation for end local
const auto LOOP_ADD_REM_LOCAL = 10u;              // Number of add/rem local changes
const auto LOOP_MOVE_LOCAL = 10u;                 // Number of move local changes
const auto LOOP_IC_LOCAL = 10u;                   // Number of initial condition local changes
const auto LOOP_IC = 10u;                         // Number of initial condition changes
const auto LOOP_IC_SWAP = 10u;                    // Number of initial condition swap changes
const auto LOOP_END_LOCAL = 10u;                  // Number of end local changes
const auto LOCAL_MOVE_RANGE = 0.1;                // The range over which local moves are made
const auto LOCAL_SINGLE_DIV = 10;                 // Number of divisions for local single updates
const auto IND_SAMP_THRESH = 0.9;                 // Threshold for individual event sampler
const auto ADD_REM_IND_MAX = 20u;                 // Max number of inds to add/rem per proposal
const auto UPDATE_RATE_FAC = 0.9;                 // Determines how quickly rates are estimated  
//const auto IND_LOCAL_TRANS_MAX = 3u;              // Max number of transition local change
const auto IND_LOCAL_TRANS_MAX = 2u;              // Max number of transition local change
const auto LI_WRONG = -100.0;                     // Likelihood penalty for wrong observation
const auto GEN_OBS_WRONG = -1000.0;               // Likelihood penalty genetic obs not made
const auto LOCAL_SMOOTH = 0.1;                    // Smoothing of transitions rate in local props
const auto PROP_SIM_NAC_INIT = 10;                // Initial value for nac for ind sim prob 
const auto PROP_SIM_PROB_FADE = 0.98;             // Used for time fading of sim prob
const auto PROP_SIM_PROB_MIN = 0.05;              // Minimum probability for ind sim proposal

const double TINY = 0.000000001;                  // Used to represent a tiny number
const double VTINY = 0.00000000000001;            // Used to represent a tiny number
const double EFFECT_MAX = 10000000;               // Sets maximum value for effect
const double EFFECT_MIN = 0.0000001;              // Sets minimum value for effect
const double CLIP_MIN = 0.99*log(EFFECT_MIN);     // Clips minimum 
const double CLIP_MAX = 0.99*log(EFFECT_MAX);     // Clips maximum 
const double DIF_THRESH = 0.00001;                // The threshold for a difference
const double THRESH_EXPAND = 100;                 // Expand threshold for certain quantities
const double DIF_THRESH_BURNIN = 0.1;             // The threshold during burnin
const unsigned int RANGE_MIN = 50;                // Minimum number to calc M

const unsigned int H_BIN = 10;                    // Used for distributions in cumulative prob

const double TWO_TINY = 2*TINY;                   // Used to represent a tiny number
const double SMALLISH = 0.01;                     // Used to represent quite a small number
const double SMALL = 0.00001;                     // Used to represent a small number
const double LARGE = 1000000000;                  // Used to represent a big number
const double LARGISH = 10000.0;                   // Used to represent a pretty large number 
const double ALPHA_ZERO = -0.000001;              // Used to represent a zero value for alpha in dirichlet
const int LARGE_INT = 1000000000;                 // Used to represent a big integer
const double INFY = 1000000001;                   // Used to represent infinity
const double UNDEF = 1000000002;                  // Used to represent undefined
const double ALMOST_ONE = 0.9999999999999;        // Almost one   
const double LOG_THRESH = 0.00001;                // The threshold below which logs not calculated
const double PROB_MOD = 0.001;                    // Avoids zero probability in local event props
const double LOW_BOUND = 0;                       // The lower bound for observation probability
const double UP_BOUND = 1;                        // The lower bound for observation probability
const double OBS_COMP_MIN = 0.00001;              // Minimum value for observed compartment
const double PROP_JOIN_COR_MIN = 0.8;             // The threshold corrlelation above which proposals join
const unsigned int SEED_MAX = 10000;              // The maximum seed number

const unsigned int TRUNC_MAX = 40;                // Truncates strings longer than this

const auto PHI_LIM = 1.0;                         // The final inverse temperature

const auto CHECK_THIN = 20u;                      // Determines how often MCMC is checked
const auto SD_MIN = TINY;                         // The minimum value for sd
const auto SD_MAX = LARGE;                        // The maximum value for sd
const auto CV_MIN = 0.01;                         // The minimum value for cv
const auto CV_MAX = 10.0;                         // The maximum value for cv
const auto MEAN_MIN = VTINY;                      // The minimum for dist mean
const auto MEAN_MAX = LARGE;                      // The maximim for dist mean
const auto NORM_MEAN_MIN = -LARGE;                // The minimum for normal mean
const auto NORM_MEAN_MAX = LARGE;                 // The maximum for normal mean
const auto SCALE_MIN = TINY;                      // The minimum for dist scale
const auto SCALE_MAX = LARGE;                     // The minimum for dist scale
const auto SHAPE_MIN = 0.1;                       // The minimum for dist shape
const auto SHAPE_MAX = 100.0;                     // The maximum for dist shape
const auto ALPBETA_MIN = 0.01;                    // The minumum for alpha or beta
const auto ALPBETA_MAX = 100;                     // The minumum for alpha or beta
const auto LAM_MIN = 0;                           // The minimum value for poisson
const auto LAM_MAX = LARGE;                       // The maximum value for poisson
const auto EXP_MEAN_MIN = TINY;                   // Minimum mean for exponential
const auto TIME_MIN = TINY;                       // The minimum time
const auto RATE_MIN = TINY;                       // The minimum rate
const auto P_MIN = TINY;                          // The minimum probability
const auto P_MAX = 1-TINY;                        // The maximum probability
 
const unsigned int SIM_NUM_DEFAULT = 1;           // The default simulation number
const unsigned int PPC_NUM_DEFAULT = 200;         // The default number of ppc simulations
const auto ALG_DEFAULT = DA_MCMC;                 // Default inference algorithm
const unsigned int SEED_DEFAULT = 0;              // The default seed
const double ANNEAL_POWER_DEFAULT = 4;            // Default annealing power
const double ANNEAL_RATE_DEFAULT = 0.01;          // Rate at which annealing is done
const unsigned int TI_DIV_MAX = 10000;            // Maximum number of time divisions

const unsigned int PARAM_OUTPUT_MAX_DEFAULT = 1000;// The default maximum number of tensor elements to be output
const auto INDMAX_DEFAULT = 20000u;               // The default maximum number of individuals
const auto BURNIN_FRAC_DEFAULT = 20.0;            // The default percentage burnin
const auto MCMC_SAMPLE_DEFAULT = 5000u;           // The default number of MCMC samples
const auto MCMC_OP_PARAM_DEFAULT = 1000u;         // The default number of output parameters
const auto MCMC_OP_STATE_DEFAULT = 200u;          // The default number of output states  
const auto MCMC_CHAIN_PER_CORE_DEFAULT = 1u;      // The default chains per core
const auto PAS_PART_PER_CORE_DEFAULT = 1u;        // The default particles per core
const auto PAS_GEN_UPDATE_DEFAULT = 1.0;          // The default updates per generation percent  

const auto ABC_SAMPLE_DEFAULT = 1000u;            // The default number of ABC samples
const auto ABC_ACFRAC_DEFAULT = 0.1;              // The default acceptance fraction for ABC
const auto ABCSMC_ACFRAC_DEFAULT = 0.5;           // The default acceptance fraction for ABCSMC
const auto ABCSMC_KERNEL_DEFAULT = 0.5;           // The default kernel size for ABCSMC 
const auto ABCSMC_GEN_DEFAULT = 5u;               // The number of generations for ABCSMC

const auto UPDATE_PROP_PROB = 10u;                // Sets how often the propsal prob is updated
const auto PROPOSAL_PROB_LIM = 0.05;              // The minimum probability for proposal
const auto PROP_LOCAL_W = 0.33;                   // The weight given to local proposals
const auto PROP_MBPII_W = 0.5;                    // The weight given to MBPII
const auto PROP_MBPIC_W = 0.5;                    // The weight give to MBP IC 
const auto POP_SI_LIM = 0.5;                      // Minimum proposal size when changing integer populations

const unsigned int NEXT_IND_EVENT_DIV = 100;      // Number of divisions next_event divided into

# define MM_PI  3.14159265358979323846             // Value for pi

const double PERCENT_STEP = 2.5;                  // Percentage per dot in running bar

const vector <unsigned int> factori = {1,1,2,6,24};// Fatorial N! for first 5 numbers

/************************** Text constants ******************************/

const vector <string> alphabet = {"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"};


const string chnotallowed = "\"$&=~#\\→";     // Character not allowed in equations
const string notparam_list = "+-*×/0123456789.{}<>〈〉() $′→;";   // Sets not a parameter character
const string paramend_list = "+-*×/|.{<>〈〉}[] $"; // Sets characters which end parameter
const string name_notallow = "|\"*×_ {<>〈〉}[]()=~′→;$&#\\";// Character not allowed in names
const string fe_char = "μ";                       // Character used for fixed effect parameter
const auto name_ch_max = 40;                      // Maximum number of characters allowed
const string endli = "\n";
const string dist_matrix_name = "D";              // The name of a distance matrix
const string iden_matrix_name = "δ";              // The name of identity matrix
const string iden_matrix_name2 = "\\delta";       // The name of identity matrix
const string RN_name = "RN";                      // The name of rep. number function
const string RNE_name = "RNE";                    // The effective rep number function
const string RNC_name = "RNC";                    // The computational rep number function
const string GT_name = "GT";                      // The basic generation time
const string GTE_name = "GTE";                    // The effective generation time function
const string GTC_name = "GTC";                    // The computational generation time function
const string ALG_WARN = "ALGORITHM WARNINGS";

// Checking for Σ and ∫ is done differently in c++ because it uses two characters
const string sigma = "Σ";  
const string tint = "∫";  

const string unobserved_str = "Unobserved";       // Name for unobserved 

const string missing_str = ".";                   // Determines a missing quantity
const string not_alive_str = "!";                 // Denotes individul not alive

const string err_mess = "Please correct these errors and rerun.";
const string err_mess_sing = "Please correct this error and rerun.";
const vector <string> must_term_str = { "simulation","sim","inference","inf","posterior-simulation","post-sim"};

/************************** Object constants ******************************/

struct EqnInfo {
	EqnInfo(EqnType ty, EqnMode mo, bool obs) { type = ty; mode = mo; obs_model = obs;}
	EqnType type;
	EqnMode mode;
	bool obs_model;
};

const vector <EqnInfo> eqn_types = { 
EqnInfo(MUT_RATE,PARAM_ONLY,true),
EqnInfo(SEQ_VAR,PARAM_ONLY,true),
EqnInfo(SE,PARAM_ONLY,true),
EqnInfo(SP,PARAM_ONLY,true),
EqnInfo(COMP_PROB,PARAM_ONLY,true),
EqnInfo(TRANS_PROB,PARAM_ONLY,true),
EqnInfo(BP,PARAM_WITH_DEP,false),
EqnInfo(SOURCE_RATE,ALL,false),
EqnInfo(SOURCE_MEAN,ALL,false),
EqnInfo(TRANS_RATE,ALL,false),
EqnInfo(TRANS_MEAN,ALL,false),
EqnInfo(TRANS_NM_RATE,ALL,false),
EqnInfo(TRANS_NM_MEAN,ALL,false),
EqnInfo(TRANS_SHAPE,ALL,false),
EqnInfo(TRANS_SCALE,ALL,false),
EqnInfo(TRANS_CV,ALL,false),
EqnInfo(REPARAM,PARAM_ONLY,false),
EqnInfo(DIST,PARAM_ONLY,false),
EqnInfo(DERIVE_EQN,DERIVE_MODE,false),
EqnInfo(REPARAM_EQN,ALL,false)
};

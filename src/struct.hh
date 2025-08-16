/// Stores commonly used data structures

#pragma once

#include <vector>

using namespace std;

#include "const.hh"
#include "hash.hh"

struct AlgWarn {                   // Stores an algorithm warning
	string te;                       // Text for warning
	unsigned int core;
	unsigned int sample;
	unsigned int num;
};

struct InfPeriod {                 // Defines an infectious period
 double start;                     // Start and end time for period
 double end;
 double end_inf;                   // The end time of when stops being infectious
 double num_inf;                   // Number of infections caused by individual
 double num_inf_t;                 // Number inf multiplied by time from infection
};

struct FRef {                      // Equation used in F for rep num
	unsigned int e;                  // Equation number
	unsigned int i;                  // popref number
};

struct SRef {                      // Equation used in S for rep num
	unsigned int e;                  // Equation number
	unsigned int m;                  // If branching this gives the nmtrans number
	TransType type;                  // Whether rate or mean
};

struct EqSign {                    // Includes equation ref with sign
	unsigned int k;
	int sign;
};

struct Fele {                      // Stores element of F matrix
	unsigned int k;
	unsigned int c_from; 
};

struct DerFunc {
	bool on;
	DerFuncType type;                // The type derived function
	string name;                     // The name
	string cont;                     // The content for specifying 
	unsigned int cl;                 // The classification
	vector <unsigned int> ref;       // References infected state
	vector < vector <unsigned int> > tau_eq; // Gives equations for calculating GT
	
	vector <FRef> Feq_ref;           // References equations used to calculate F
	vector <SRef> Seq_ref;           // References equations used to calculate S
	vector < vector < vector <Fele> > > F_eq;
	vector < vector < vector <EqSign> > > S_eq;
	string warn;
};

struct Warn {                      // Stores warning information from ESS and GR
	unsigned int th;                 // The parameter affected
	double num;                      // The number 
};

struct Sum {                       // Stores information about an equation sum
	vector <string> dep;
	unsigned int i_start;
	unsigned int i_end;
};

struct ParamRef {                  // Used to reference parameter
	ParamRef(){ th = UNSET; index = UNSET;}
	unsigned int th;                 // The number of the parameter
	unsigned int index;              // The index where to find 
};

struct DeriveRef {                 // Used to reference derived quantities
	DeriveRef(){ i = UNSET; index = UNSET; ti = UNSET;}
	unsigned int i;                  // The number of the parameter
	unsigned int index;              // The index where to find 
	unsigned int ti;                 // The time step
};

struct GradRef {                   // References population gradient
	unsigned int ref;                // References a Markov equation
	unsigned index;                  // References a population in Markov equation
};

struct PopAffect {                 // Lists all populations which affect DIV_VALUE_LINEAR_AFFECT
	unsigned int po;                 // Population number
	vector <GradRef> pop_grad_ref;   // References population in Markov equation 
};

struct LinearProp {                // Store information about DIV_VALUE_LINEAR_AFFECT proposal
	vector <unsigned int> me;        // Markov eqn ref
	vector <unsigned int> me_no_pop; // Which no_pop calculations need to be done
	vector <unsigned int> no_pop_ref;// References no_pop which needs to be updated
	vector <PopAffect> pop_affect;   // References populations which affect
};

struct AffectLike {                // Determines how a parameter affects likelihoods
	AffectType type;                 // The type of likelihood being affected
	unsigned int num;                // The reference number
	unsigned int num2;               // Another reference number
	vector <unsigned int> me_list;   // Stores a list of Markov equantions (used in DIV_VALUE_FAST_AFFECT)
	vector <bool> map;               // Maps all time elements which are affected
	vector <unsigned int> list;      // Lists all map which is true
	unsigned int order_num;          // Used to order terms by list priority
	LinearProp linear_prop;          // Used for DIV_VALUE_LINEAR_AFFECT
};

struct PopulationTerm {            // A term in the population
	unsigned int c;                  // Compartment
	double w;                        // Weight
};

struct PopRef {                    // References a population
	unsigned int po;                 // The population number
	unsigned int index;              // The index
};

struct PopMarkovEqnRef {           // References a Markov equation from a population
	unsigned int p;                  // The species
	unsigned int e;                  // The Markov eqn number
};

struct PopTransRef {               // References transition from a population
	unsigned int p;                  // The species
	unsigned int tr;                 // The transition number
};

struct Poss {                      // A possiblity used for sampling
	unsigned int i;                  // Reference number
	double weight;                   // Weight of possibility
};

struct IndPop {
	unsigned int c;                  // The compartment for individual
	vector <PopRef> pop_ref;         // Stores when on the population index
};

struct PopChange {                 // Stores the change in population
	unsigned int po;                 // The population
	double num;                      // The change
};

struct Population {                // Stores a population (used in an equation)
	string name;                     // The name of the population
	unsigned int sp_p;               // The species number
	bool ind_variation;              // Set if ind_eff_mult or fix_eff_mult are set
	vector <unsigned int> ind_eff_mult; // Individual effect which multiplies weight (if appropriate) 
	vector <unsigned int> fix_eff_mult; // Any fixed effects which modify weight
	vector <PopulationTerm> term;    // The global compartments which contribute to the popualtion
	vector <PopMarkovEqnRef> markov_eqn_ref;// References Markov equations which include population
	vector <PopTransRef> trans_ref;  // References transitions which include population
	//PopChangeSpeedup speedup;        // Used for speeding up population change in update_ind
};

struct ParamProp {                 // Stores properties of a parameter
	string name;                     // The name (without the 
	vector <string> dep;             // String form of dependency
	vector <string> dep_with_prime;  // String form of dependency (with primes)
	bool time_dep;                   // Time dependency
	bool prime;                      // Set to true of anu dependencies are primes
};

struct SplineInfo {                // Stores information about a parameter spline
	bool on;                         // Detemines if spline is switched on
	vector <double> knot_tdiv;      // Times of knots (if time variation)
	bool smooth;                     // Determines if smoothing is added
	SmoothType smooth_type;          // The type of smoothing (normal or log-normal)
	double smooth_value;             // The value used for the smoothing
};

struct SplineDiv {                 // A div in spline is calculated from a weighting of params
	unsigned int th1;                // The reference for the parameter 1
	unsigned int th2;                // The reference for the parameter 2
	double f;                        // The factor 1=th1  0=th2
};

struct Spline {                    // Stores an individual spline
	string name;                     // Name of the spline
	bool constant;                   // Determines if spline is constant
	unsigned int th;                 // The parameter
	unsigned int index;              // The index giving the dependency
	vector <unsigned int> param_ref; // References the parameters along the spline
	SplineInfo info;                 // Stores information about spline (copied from par)
	vector <SplineDiv> div;          // A division within the spline
	vector <double> const_val;       // If spline is constant then works out what value it should have
};

struct SplineValue {               // Stores values for a spline
	vector <double> val;             // The value down the spline
	vector <double> store;           // The value down the spline
};

struct CompPos {                   // Used to go through comparmtental possibilities
	vector <string> list;            // Stores a list of compartmental possibilities
	unsigned int index;              // Stores which of these possibilities currently looked at
}; 

struct DepInfo {                   // Stores dependency info froman equation
	unsigned int iend;               // The end position for the dependency
	vector <string> spl;             // Splits the dependency into individual elements
	vector <unsigned int> ipos;      // The positions for the different depedencies
	string warn;                     // Stores any warning message
};

struct SwapTemp {                  // Constructs a swap template
	string te;                       // Text fragment
	unsigned int num;                // Number to substitute
};

struct SwapResult {                // Stores the result of doing a substitution of indeces for values
	string warn;                     // Stores the warning message
	vector <bool> done;              // Registers true if a swap has been performed
	vector <SwapTemp> swap_temp;     // The template for swapping
};

struct DepConv {                   // Used to convert dependencies from index to value
	string before;                   // The string before change
	string after;                    // The string after change
};

struct EquationInfo {              // Stores information about an equation (prior to equation being evaluated)
	EquationInfo(){ te = ""; p = UNSET; cl = UNSET; infection_trans = false; error=false; emsg="";}
	string te_raw;                   // The rate text value (as seen in the input file)
	string te;                       // The text value after processing
	unsigned int p, cl;              // The species and classification for equation (if appropriate)
	//unsigned int c;                  // The global compartment (used for sum max function)
	EqnType type;                    // The type of the equation
	bool infection_trans;            // Set if infection transition (used to calculate transtree likelihood)
	double value;                    // The numerical value (if appropriate)
	unsigned int eq_ref;             // The equation (if appropriate)
	unsigned int line_num;           // Stores the line number the equation comes from
	bool error;                      // Set is there is an error in the equation
	string emsg;                     // Stores the error message
};

struct IEname {                    // Stores individual effect name
	string name;                     // The name of the individual effect
	unsigned int index;              // Indexes the ie in the species
};

struct Amatrix {                   // Stores an A matrix
	bool set;                        // Determines if the A matrix has been set
	vector < vector <double> > value;// The values in the matrix
	vector <string> ind_list;        // The names of individuals in the system
	Hash hash_ind_list;              // Hash table for individual list
};

struct IndEffect {                 // Store infortmation about an individual effect (in species)
	string name;                     // The name of the individual effect
	unsigned int index;              // References ind_eff_group
	unsigned int num;                // The number within ind_eff_group
	vector <unsigned int> markov_eqn_ref; // References any Markov equations the ie appears in
	vector <unsigned int> nm_trans_ref;   // Indexes nm_trans involving ies
	vector <unsigned int> nm_trans_incomp_ref;// Indexes nm_trans_incomp involving ies
	vector <unsigned int> pop_ref;   // Populations that ie appear in
	unsigned int line_num;           // The line number in input file 
};

struct IEgroup {                   // Stores information about a group of individual effects
	vector <IEname> list;            // The list of individual effect names
	Amatrix A_matrix;                // The A matrix
	vector < vector <unsigned int> > omega;// Matrix of parameter values 
	unsigned int line_num;           // Stores the line number the equation comes from
	bool ppc_resample;               // Determines if individual effects resampled 
	vector <unsigned int> markov_eqn_ref;  // Indexes any Markov equations involving ies
	vector <unsigned int> nm_trans_ref;    // Indexes nm_trans involving ies
};

struct NZElement {                 // Stores non-zero elements of a matrix
	unsigned int i;                  // x position
	double value;                    // Value
};

struct IEGroupSampler {            // Stores a sample for an individual effect group
	vector < vector <double> > omega;// The covariance matrix omega
	vector < vector <double> > omega_Z;  // The cholesky matrix for omega (if doing simulation)
	vector < vector <double> > omega_inv;// The inverse matrix for omega (if doing inference)
	vector < vector <double> > A;    // The A relationship matrix
	vector < vector <double> > A_Z;  // The cholesky matrix for A
	vector < vector <double> > A_inv;// The inverse relationship matrix
	vector <double> A_inv_diag;      // The diagonal elements of Ainv
	vector < vector <NZElement> > A_inv_nonzero;// The non-zero  elements
};

struct Dependency {                // A dependency in the model
	string index;                    // The index in the model
	string index_with_prime;         // The index in the model including any prime
	vector <string> list;            // List of possible values (either compartments or times along spline)
	Hash hash_list;                  // Generates a hash table for list
	unsigned int mult;               // Multiplicative factor when accessing an element
};

struct Prior {                     // Defines a parameter prior
	Prior(){ error = ""; type = UNSET_PR;}
	string name;                     // Stores text name
	PriorPos type;                   // The type of prior
	vector <EquationInfo> dist_param;// A list of equations to specify prior distribution parameters
	string in;                       // Descript of what prior is in (e.g. name of parameter)
	string error;                    // Stores error message if declared incorrectly
};

struct ParamVecEle                 // Stores information about an element in param vec
{
	string name;                     // The name of the parameter
	unsigned int th;                 // The number of the parameter
	unsigned int index;              // The index where to find 
	unsigned int prior_ref;          // References the prior
	ParamVariety variety;            // The parameter variety (copied from param)
	bool ppc_resample;               // Sets if parameter gets resampled for ppc
	bool prop_pos;                   // Set if it is possible to do a proposal on this parameter
	bool omega_fl;                   // Detemines if in omega
	vector <AffectLike> affect_like; // Determines how parameter affects likelihoods
	unsigned int ref;                // Reference param_vec_prop
};


struct ErrorMess {                 // Stores error messages when importing 
	unsigned int line_num;           // The number of the line in the input file
	string error;                    // The error message
	ErrorType type;                  // The type of the error (fatal or just a warning)
};

struct LatLng {                    // Denote a point with latitude and longitude
	double lat;                      // Latitude
	double lng;                      // Longitude
};

struct Pgon {                      // Denotes a polygon for a boundary
	vector <LatLng> point;           // Points on boundary
};

struct Boundary {
	vector <string> name;            // Potential names for boundary (taken from "properties"
	vector <Pgon> polygon;           // Boundary made up of polygons
};

struct GeoJSON {                   // Stores geoJSON information
	string file;                     // File name
	vector <Boundary> boundary;      // Boundary information
};

struct Details {                   // Stores details about simulation/inference/PPC
	double t_start;                  // The start time
	double t_end;                    // The end time
	double ppc_t_start;              // The start time for ppc
	double ppc_t_end;                // The end time for ppc
	double inf_t_end;                // The end time for inference
	double dt;                       // The time step 
	unsigned int T;                  // The total number of divisions
 	unsigned int number;             // Number of simulations to be performed
	unsigned int output_param;       // The number of output parameter samples (MCMC)
	unsigned int output_state;       // The number of output state samples (MCMC)
	unsigned int sample;             // Number of samples to be done (MCMC)
	double accfrac;                  // The acceptance fraction (ABC)
	unsigned int numgen;             // The number of generations (ABC-SMC)
	double kernelsize;               // The size of the kernel distribution (ABC-SMC)
	unsigned int individual_max;     // The maximum number of individuals (individual-based models)  
 	Algorithm algorithm;             // The algorithm used for simulation / inference
	bool stochastic;                 // Detemines if stochastic
	double burnin_frac;              // The fraction of samples used for burnin
	AnnealType anneal_type;          // The type of annealing performed
	double anneal_power;             // The power (if power annealing is used)
	double anneal_rate;              // The rate of annealing (if used)
	unsigned int gen_update;         // Number of updates per generation (PAS-MCMC)
	string ppc_resample;             // String determining which parameters must be resampled    
	unsigned int param_output_max;   // The maximum size of tensor which can be output
	unsigned int seed;               // Seed used to initialise RNG
	unsigned int nchain;             // The number of chains / particles
	unsigned int num_per_core;       // The number of samples / chains / particles per core
	bool diagnostics_on;             // Determines if MCMC diagnostics printed
};

struct FilterCla {                 // Used to filter a classification
	vector <string> comp_prob_str;   // Sring giving the probability of being in compartment
	FilterClaType type;              // Type of filter e.g. all comps, n file or specified by comp 
};

struct Filter {                    // Used to filter populations etc.
	vector <FilterCla> cla;          // Stores filters for each classification
};

struct Table {                     // Loads a table
	string file;                     // The file from which the tables was loaded
	unsigned int ncol;               // The number of columns
	unsigned int nrow;               // The number if rows
	vector <string> heading;         // The headings for the columns
	vector <vector <string> > ele;   // The elements of the table
	bool error;                      // Set if there is a problem loading the table
};

struct DiagTestSens {              // Determines compartments test sensitive to 
	unsigned int cl;                 // The classification
	vector <bool> comp;              // The compartments the test is sensitve to are set to true
};

struct ObsModel {                  // Stores the observation model assocaited with data
	ObsModelType type;               // The type of obsevation model
	string Se_str;                   // Gives the sensitivity of the test
	string Sp_str;                   // Gives the specigivity of the test
	string diag_pos;                 // The string for a positive diagnostic test (for diagnostic test data)
	string diag_neg;                 // The string for a negative diagnostic test (for diagnostic test data)
	DiagTestSens diag_test_sens;     // Determines compartments test sensitive to (for diagnostic test data)
	double percent;                  // percent value when 'error' normal obs model used 
	double sd;                       // sd value when 'error' normal obs model used 
	double p;                        // p value when 'error' neg-binomial obs model used 
};

struct LoadCol {                   // Used to store information about loading data columns
	LoadCol(string he, string de, ElementType ty, unsigned int cla = UNSET){
		heading = he; desc = de; type = ty; cl = cla;
	}
	string heading;                  // The heading for the column
	string desc;                     // A description
	ElementType type;                // The type of element
	unsigned int cl;                 // Stores classification number (if appropriate)
};

struct CompRef {                   // References a particular compartment
	unsigned int p;                  // Species number
	unsigned int cl;                 // Classification number
	unsigned int c;                  // Compartment number
	string error;                    // Error message
};

struct DataSource {                // Stores information about a data source
	string name;                     // The name of the data source 
	unsigned int index;              // Index within vector
	Command cname;                   // The name of the command
	Table table;                     // The table which stores the data
	unsigned int focal_cl;           // The focal classification (for init_pop data only)
	unsigned int p;                  // The species number
	unsigned int cl;                 // The classification number (for MOVE_IND, TRANS_DATA)
	InitPopType init_pop_type;       // Used to specifed the type of initial population
	Prior pop_prior;                 // The total population prior (INIT_POP) 
	TimeRange time_range;            // Determines what time range is used (used in TRANS_DATA)
	double time_start, time_end;     // Stores the time range (used in TRANS_DATA)
	string filter_str;               // String which filters compartment
	string filter_trans_str;         // String which filters transitions 
	ObsModel obs_model;              // The observation model (used in POP_DATA, POP_TRANS_DATA)
	string SNP_root;                 // Shows the root for the SNP data (used in GENETIC_DATA)
	GenDataType gen_data_type;       // Show type of genetic data (used in GENETIC_DATA)
	string mut_rate_str;             // The genetic mutation rate (used in GENETIC_DATA)
	string seq_var_str;              // The initial sequence variation (used in GENETIC_DATA)
	vector <LoadCol> load_col;       // Stores information about the loaded columns
	unsigned int line_num;           // The line number which loads up the data source
};

struct IslandRef {                 // References an island (a group of connected compartments)
	unsigned int i;                  // The island number
	unsigned int c;                  // The compartment within the island
};

struct Compartment {               // Information about a compartment
	string name;                     // The name of the compartment
	bool markov_branch;              // Determines if compartment is a branching if Markovian
	bool branch;                     // Determines if branch prob exists
	bool branching;                  // Determines if branching
	bool all_branches;               // Determines if all branches are specified
	bool erlang_hidden;              // Determines if hidden Erlang compartment
	string erlang_source;            // Determines Erlang source compartment
	unsigned int erlang_c_start;     // References the c compartment at the start of erlang
	double x, y;                     // The position of the compartment (optional, used if distance is in the model)
	double lng, lat;                 // The position of the compartment (optional, used if distance is in the model)
	IslandRef island_ref;            // References which island the compartment is in
	CompInfected infected;           // Determines if compartment is infected or not (or unset)
	vector <unsigned int> tr_leave;  // Transitions leaving this compartment
	vector <unsigned int> tr_enter;  // Transitions entering this compartment
};

struct ParamElement {              // Stores information about a parameter element
	ParamElement(){ param_vec_ref = UNSET;}

	EquationInfo value;
	bool used;                       // Determines if specific parameter is used	
	unsigned int prior_ref;          // References prior
	vector <ParamRef> parent;        // Parameters dependent on
	vector <ParamRef> child;         // Parameters on which dependent
	unsigned int param_vec_ref;      // References the value on param_vec_ref
};

struct Param {                     // Stores a model parameter
	string full_name;                // The full name (as defined in the input file)
	string name;                     // The name after removing dependencies
	vector <Dependency> dep;         // Dependency for parameter
	bool time_dep;                   // Time dependency
	SplineInfo spline_info;          // Stores information about spline (if used)
	bool sim_sample;                 // Set to true if value sampled from the distribution (simulation)
	ParamVariety variety;            // Determines the variety of parameter	
	unsigned int N;                  // The number of elements in the parameter
	bool trace_output;               // Determines if parameter is output (if not too big and not constant)
	
	unsigned int default_prior_ref;  // Default prior reference
		
	vector <unsigned int> element_ref; // This references a given element
	vector <ParamElement> element;   // This stores information about elements in parameter
	
	vector <double> weight;          // Different weights given to categories for factor

	vector <double> cons;            // Stores constant values
	string reparam_eqn;              // Used to store reparameterisation eqn
	
	unsigned int line_num;           // The line num in input file when parameter is defined
	bool factor;                     // Set if parameter is a param_mult
	unsigned int param_mult;         // If the parameter is multiplied by a factor this is th
	bool used;                       // Determines if param used in model
	bool not_set;                    // Detemines if not set

	bool cat_factor;                 // Set if parameter is a categorical factor
	bool cat_factor_weight_on;       // Determines if there is a weight factor
	
	unsigned int get_param_vec(unsigned int i) const; 
	const vector <ParamRef>& get_child(unsigned int i) const;
	const vector <ParamRef>& get_parent(unsigned int i) const;
	unsigned int get_eq_ref(unsigned int i) const;
	void all_elements();
	void add_element(unsigned int i); 
	void set_prior(unsigned int i, unsigned int prior_ref);
	void add_parent(unsigned int i, const ParamRef &pr);
	void add_child(unsigned int i, const ParamRef &pr);
	double get_value(unsigned int i) const;
	void set_used(unsigned int i);
	void set_value_te(unsigned int i, string te);
	string get_value_te(unsigned int i) const;
	void set_value_eqn(unsigned int i, const EquationInfo &val);
	//void set_weight(unsigned int i, double w);
	//double get_weight(unsigned int i);
	unsigned int get_prior_ref(unsigned int i) const;
	bool exist(unsigned int i) const;
	unsigned int add_cons(double val);
	void set_cons(unsigned int i, double val);
};

struct CompGlTransGroup {          // Information about transitions from a compartment (in a cl)
	bool branch;                     // Whether branching happens
	bool all_branches;               // Determines if all branches are specified
	bool markov;                     // Set if all transitions leaving are rates
	vector <unsigned int> tr_list;   // References global transitions leaving compartment group
};

struct CompGlobal {                // A global compartment (combines all classifications)
	string name;                     // The name of the global compartment
	unsigned int infected;           // Determines if an infected state (for trans-tree)
	bool erlang_hidden;              // Determines if added to make Erlang distribution
	unsigned int erlang_c_start;     // References the c compartment at the start of erlang
	vector <unsigned int> cla_comp;  // References the classification compartment
	vector <PopRef> pop_ref;         // References population and point on list
	vector <unsigned int> pop_ref_simp;       // References just population 
	vector <unsigned int> tr_enter;           // Transitions entering compartment
	vector <unsigned int> tr_leave;           // All the transitions leaving compartment 
	vector <unsigned int> tr_leave_markov;    // All Markov transitions leaving compartment 
	vector <CompGlTransGroup> tra_leave_group;// Information about transitions leaving [cl] 
	vector <unsigned int> me_ref;             // References markov equations leaving compartment
	vector <unsigned int> nmtransincomp_ref;  // References incomplete nm trans
};

struct Transition {                // Information about a transition in the model
	string name;                     // The name of the transitions
	unsigned int i;                  // The initial compartment (set to UNSET for source)
	unsigned int f;                  // The final compartment (set to UNSET for sink)
	TransType type;                  // The type of the transition
	TransVariety variety;            // Determines if normal, source or sink
	bool erlang_hidden;              // Determines if hidden Erlang transitions
	BpType bp_set;                   // Determines if branching probabilty has been set
	bool branch;                     // Determines if this is one of multiple branches leaving a compartment
	bool all_branches;               // Determines if all branches are specified
	EquationInfo bp;                 // References branching probability equation
	vector <EquationInfo> dist_param;// References a list of equations to specify transition distribution parameters
	unsigned int line_num;           // Stores the import line (for diagnostic error messages)
};
  
struct IndFacRate {                // Information about individual factors affecting rate
	IndFacRateType type;             // The type of individual effect change    
	unsigned int e;                  // Index                       
};
  
struct NMTrans {                   // Stores a non-Markovian transition
	string name;                     // The name of the transition
	TransType type;                  // The type of the transition
	unsigned int bp_eq;              // References branching probability equation (if applicable)
	vector <unsigned int> bp_other_eq;// If bp_eq=BP_FROM_OTHERS this gives the list of other equations
	bool all_branches;               // Determines if all branches are specified
	vector <unsigned int> bp_all_eq; // Gives equations for all bps 
	bool precalc_nm_rate;            // Set if it is possible to precalculate rate
	vector <unsigned int> dist_param_eq_ref;// Equations for distribution
	vector <IndFacRate> ind_fac_rate;// Stores individual factor for rate
	unsigned int trans_incomp_ref;   // References incomplete NM transition
};

struct NMTransIncomp {             // Stores an incomplete non-Markovian transition
	string name;                     // The name
	unsigned int cl;                 // The classification which undergoes the incomplete trans
	bool all_branches;               // Determines if all branches are specified
	vector <unsigned int> nmtrans_ref;// References nmtrans which come out
};

struct TransInfection {            // Information about transition infection (for trans trees)   
	TransInfectionType type;         // Type of transition (infection, recover or no change)
	vector <unsigned int> pop_ref;   // Stores populations
};

struct InfSourceSampler {          // Used to sample the source of infection (for trans trees)
	vector <double> val_store;       // The value for different sources
	vector <double> val_sum_store;   // The cumulative value
	double sum;                      // Total value for sum
	unsigned int sample_inf_source() const;
	double prob_inf_source(unsigned int j) const;
};

struct TransGlobal {               // A global transition (combines all classifications)
	string name;                     // The name of the transition
	unsigned int cl;                 // The classification
	unsigned int tr;                 // The transition number in the classifciation
	unsigned int i;                  // The initial compartment (set to UNSET for source)
	unsigned int f;                  // The final compartment (set to UNSET for sink)
	TransInfection infection;        // Under trans-tree provides informatransition is associates with infection
	bool erlang_hidden;              // Determines if hidden Erlang transitions
	TransType type;                  // The type of the transition
	EventType ev_type;               // NM_TRANS_EV or M_TRANS_EV
	TransVariety variety;            // Determines if normal, source or sink
	bool branch;                     // Determines if this is one of multiple branches leaving a compartment
	bool all_branches;               // Determines if all branches are specified
	EquationInfo bp;                 // References branching probability equation
	vector <unsigned int> bp_other_eq;// A list of other equations used to calculate bp = (1-sum other eqns) 
	BpType bp_set;                   // Determines if branching probabilty has been set
	vector <EquationInfo> dist_param;// References a list of equations to specify transition distribution parameters
	unsigned int nm_trans_ref;       // References nm_trans (if non-Markovian)
	unsigned int markov_eqn_ref;     // References Markovian equation (if Markovian)
	vector <unsigned int> tform;     // Transforms transition based on a shift in c
	bool time_vari;                  // Determines if time varying
	unsigned int line_num;           // Stores the import line (for diagnostic error messages)
};

struct MEref {                     // References an equation and time in markov_eqn
	unsigned int e;                  // Markov_eqn ref
	unsigned int ti;                 // Division time along Markov equation 
};

struct IslandTrans {               // Transitions within an island
	unsigned int tr;                 // References the transion
	unsigned int cf;                 // References the final compartment (in island units)
	bool nm_trans;                   // Determines if a non-Markovian transition
	vector < vector <MEref> > markov_eqn_ref; // References Markov equation [ti][cgl] 
	vector <unsigned int> nm_trans_ref;       // References non-Markovian transion [cgl]
};

struct IslandComp {                // Transitions within an island
	unsigned int c;                  // References the compartment
	bool nm_trans;                   // Set if non-Markovian transitions leaving compartment
	vector <IslandTrans> leave;      // Transitions which leave the compartment 
};

struct Island {                    // An island is a group of compartments not connected to rest
	vector <IslandComp> comp;        // Compartment in island
};

struct TrSwap {                    // Information about locally swapping transition events
	TrSwapType type;                 // The type of the swap
	unsigned int trc;                // The transition or compartment
};

struct Swap {                      // Allows for local swap of individual events
	string name;                     // The name of the swap
	unsigned int c;                  // The initial compartment (in classification)
	vector <TrSwap> start;           // The starting seqence of events
	vector < vector <TrSwap> > swap; // Possibilities for swapping
	vector <unsigned int> swap_rep_ref;// References swap_rep
};

struct SwapRep {                   // Swaps from one possiblity to another
	string name;                     // The name of the swap replacement   
	unsigned int swap_ref;           // References a swap
	unsigned int index;              // Gives the index of that swap
};

struct SwapProp {                  // Stores information about swap proposals
	unsigned int ntr;                // Number of times proposal has been performed
	unsigned int ntr2;               // Number of times after filter
	unsigned int nac;                // Number of times proposal has been accepted
	unsigned int nfa;                // Number of times proposal fails
	unsigned int nzero;              // Number of times proposal has zero acceptance prob
	unsigned int nfilt;              // Number of times proposal filter
};

struct GeneticProp {               // Used for local genetic proposals to mutation number
	GeneticProp(){ si = 1; ntr = 0; nac = 0;}
	double si;                       // The size of the proposals
	unsigned int ntr;                // Number of times proposal has been performed
	unsigned int nac;                // Number of times proposal has been accepted
};

struct Classification {            // Stores details of a species classification
	string name;                     // The name of the classification
	vector <Compartment> comp;       // Stores details of compartments
	Hash hash_comp;                  // Stores comp in a hash table
	unsigned int ncomp;              // The number of compartments
	vector <Transition> tra;         // Stores details of transitions
	unsigned int ntra;               // The number of transitions
	string index;                    // The mathematical index used to represent different compartments
	vector <Island> island;          // Compartments are split into "island" which are unconnected
	
	// Alternative sequences between two classications [ci][cf][seqnum]
	// ci=C and cf=C are used to represent UNSET
	
	// Determines if swapping between alternative local sequences is possible
	// ci=C and cf=C are used to represent UNSET
	vector < vector <bool> > tr_swap_possible;          // [ci][cf]
	
	//vector < vector < vector < vector <TrSwap> > > > tr_swap; 
	
	Hash hash_swap;                  // Allows to find local individual proposals
	vector <Swap> swap;              // Allows for local individual proposals
	vector <SwapRep> swap_rep;       // Allows for local individual proposals
	Coord coord;                     // The coordinate system used (cartesian on latlng)
};

struct InitCondValue {             // Sampled initial condition used in the state
	vector <unsigned int> cnum;      // Initial pop in different compartments (rounded)

	// With focal classification
	vector <unsigned int> N_focal;               // The number if different compartments
	vector <unsigned int> N_focal_unobs;         // The numnber of unobserved individuals
	vector < vector <unsigned int> > cnum_reduce;// Population in reduced compartments
 	vector < vector <double> > frac_focal;       // The fraction in other classifications
	vector <double> frac_comb;                   // The fraction combing all non-focal compartments

	// Without focal classification
	unsigned int N_total;            // Total number of individuals entering
	unsigned int N_total_unobs;      // The number of unobserved individuals entering
	vector <double> frac;            // Initial fraction in different compartments
};

struct CompRedRef {                // Obtains reduced reference from cgl
	unsigned int c;                  // Compartment in focal cl
	unsigned int cred;               // c in reduced compartments (comp_reduce)
};

struct InitCond {                  // Stores information about initial conditions in data
	InitPopType type;                // The type of initial condition (fixed or dist)
	
	// Used in INIT_POP_FIXED
	vector <unsigned int> cnum;      // The initial populations
	
	// Used in INIT_POP_DIST with focal
	unsigned int focal_cl;                 // The focal classification
	vector <Prior> comp_prior;
	vector < vector <double> > alpha_focal;// Gives the alpha values 
	unsigned int mult;                     // Multiplicative factor
	unsigned int N_reduce;                 // Number of elements 
	vector <unsigned int> comp_reduce;     // The global compartments
	vector <CompRedRef> comp_reduce_ref;   // Gets reduces comparmtnets (combining non-focal)
	
	// Used in INIT_POP_DIST without focal
	Prior pop_prior;                 // Gives the overall population prior
	vector <double> alpha;           // Gives the alpha values
};

struct Xvector {                   // Stores covariate in fixed effect
	vector <string> ind_list;        // The individuals
	Hash hash_ind_list;              // Hash table for the individual list
	vector <double> value;           // The covariate
};

struct FixedEffect {               // Stores an individual fixed IndEffMult
	string name;                     // The name of the parameter
	unsigned int th;                 // The number of the parameter
	Xvector X_vector;                // The covariate as loaded from the file
	vector <unsigned int> markov_eqn_ref;  // References any Markov equations the fe appears in
	vector <unsigned int> nm_trans_ref;    // Indexes nm_trans involving fes
	vector <unsigned int> nm_trans_incomp_ref;// Indexes nm_trans involving fes
	vector <unsigned int> pop_ref;         // Populations that fe appear in
	unsigned int line_num;                 // The line number
};

struct MarkovEqnRef {              // References when the individul is in the Markov equation
	unsigned int e;                  // The Markov equation
	unsigned int index;              // The index in the time division
};

struct IndEventRef {               // References an individual event
	unsigned int i;                  // Individual number
	unsigned int e;                  // Event number
};

struct NMEventRef {                // References an individual event
	unsigned int i;                  // Individual number
	unsigned int e_end;              // Event number at end of transition (UNSET for incomplete)
};

struct IndInfFrom {                // Stores information about the indivividual which caused an infection
	IndInfFrom(){ p = UNSET; i = UNSET; pref = UNSET; po = UNSET; w = UNSET;};
	unsigned int p;                  // The infecting population
	unsigned int i;                  // The individual doing the infecting
	unsigned int pref;               // Stores which pop_ref in equation is causing infection
	unsigned int po;                 // Stores which population is causing infection
	double w;                        // The weight associated with population
};

struct IIFChange {                 // Stores changes to iif
	unsigned int p_ev;               // Species of infection event
	unsigned int i_ev;               // Individual of infection event
	unsigned int e_ev;               // Event number of infection event
	IndInfFrom ind_inf_from;         // New information
	unsigned int ti;                 // Time division
	double Li_markov_bef;            // The likelihood before
	double Li_markov_aft;            // The likelihood after
};

struct Event {                     // A transition event
	EventType type;                  // The type of the event
	unsigned int tr_gl;              // The global transition which occurs 
	unsigned int cl;                 // The classification (M_TRANS_EV, NM_TRANS_EV, MOVE_EV)
	unsigned int move_c;             // The compartment (for entry events) (MOVE_EV)
	unsigned int c_after;            // The compartment after the transition
	double tdiv;                     // The time at which transition occurs
	double Li;                       // Likelihood for transition (for NM trans)
	double Li_bp;                    // Likelihood for branch prob (for NM trans)
	bool observed;                   // Determines if event is observed in data
	IndInfFrom ind_inf_from;         // The indivividual which caused an infection
	// If ENTER_EV then ind_inf_from.p is set to ENTER_INF
	// If infection comes from outside the system then ind_inf_from.p is set to OUTSIDE_INF
	
	unsigned int inf_node_ref;       // References an infection node
	
	// These variables are automatically dealt with in update_ind functions
	unsigned int m;                  // For NM m number, for M the Markov eqn
	unsigned int ti;                 // The time division for the origin
	unsigned int index;              // References nm_trans_ind (if appropriate)
	unsigned int e_origin;           // If NM this gives the origin event
};

struct GlobalEvent {               // Global time event (used in trans_tree proposals)
	double tdiv;                     // Event time                    
	unsigned int p;                  // Species
	unsigned int i;                  // Individual
	unsigned int e;                  // Event number
};

struct PopIndRef {                 // References an individual within a population
	unsigned int i;                  // The individual
	double w;                        // The weight
	unsigned int index;              // The index
};

struct PopnumIndRef {              // References popnum_ind
	unsigned int po;                 // The population number
	unsigned int ti;                 // The time division
	unsigned int index;              // The index
};

struct IncompNMTransRef {          // References an incomplete transition
	IncompNMTransRef(){ on = false;}
	bool on;                         // Determine if set
	unsigned int n;                  // References nm_trans_incomp
	unsigned int ti;                 // The time division
	double dtdiv;                    // The time span
	unsigned int e_begin;            // Event number at begining of transition
	double tdiv_end;                 // The end time of the transition
	unsigned int index;              // The index in nm_trans_incomp_ref
	double Li;                       // Likelihood of transition
};

struct NMIncompVal {               // Stores values for inclu   
	vector < vector <double> > ref_val;// Transition values [branching number][transition parameter]
	vector <double> bp_val;          // Branching probabilities [branching number]
};

struct Individual {                // Stores information about an individual
	IndType type;                    // Whether the individual is observed or not
	string name;                     // The name of the individual
	vector <Event> ev;               // Events which move individual between compartments
	vector <double> ie;              // Individual effect values
	vector <double> exp_ie;          // Expoential of individual effect values e^(ie)
	vector <double> X;               // The X value for the fixed effects
	vector <double> exp_fe;          // Expoential of individual fixed effect values e^(X*fe)
	vector <PopnumIndRef> popnum_ind_ref; // References popnum_ind (for transtree)
	vector <IncompNMTransRef> incomp_ref; // Stores incomplete transitions [cl]
	bool init_c_set;                 // Determines if the individual compartment set
	
	// Only used during simulation
	vector <MarkovEqnRef> markov_eqn_ref;  // References Markov equations that the individual is on
};


struct SimTrigEvent {              // An event trigger (e.g. ind enters, leaves, NM transition)
	SimTrigEventType type;           // The type of event
	unsigned int c;                  // The state (if entering)
	unsigned int i;                  // The individual on which the event happens
	unsigned int trg;                // The global transtion
	double tdiv;                     // The time the event happens
};

struct SimTrigEventDiv {           // Stores a division of non-Markovian next events (simulation) 
	vector <SimTrigEvent> ev;        // Events
};

struct TrigEventRef {              // Used when simulating new event sequence for individual
	TrigEventType type;              // The type of event
	unsigned int ref;                // If nm_trans then gives tr_gl, if trans data given obs num
	double tdiv;                     // The time the event happens
};

struct FutureNMEvent {             // Future  non-Markovian event (in simulation proposal)
	unsigned int trg;                // The global transtion
	double tdiv;                     // The time the event happens
};

struct MarkovNode {                // Stores Markov node (used during simulation)
	vector <unsigned int> child;     // Child nodes
	unsigned int parent;             // Parent node
};

struct MarkovTree {                // Stores a tree to be able to sample Markov events
	vector <MarkovNode> node;        // The nodes in the tree
};

struct IndTransRef {               // Stores individual transition associated with Markov eqn
	unsigned int i;                  // The individual number 
	unsigned int index;              // The index in ind.markov_eqn
	unsigned int tr_gl;              // The global transition number
};

struct DivIndRef {                 // References an individual event on a Markov eqn
	unsigned int i;                  // Indiviudal number
	unsigned int index;              // Index the event number
};

struct MarkovEqnDiv {              // Division in Markov timeline
	double value;                    // The value of the equation for the division
	double indfac_int;               // The integral of individual factors across div
	vector <DivIndRef> ind_trans;    // The individuals which undego transitions
};

struct MarkovEqn {                 // Stores information about the Markov equations
	unsigned int eqn_ref;            // The equation number
	bool time_vari;                  // Determines if equation has time variation
	bool source;                     // Set to true if equation is for a source
	bool rate;                       // Set if rate, unset it mean
	vector <unsigned> source_tr_gl;  // Global transition for source (if applicable)
	bool ind_variation;              // Set to true if there is individual variation
	bool infection_trans;            // Determines if relates to an infection transition (used for trans tree likelihood)
	vector <unsigned int> ind_eff_mult; // References an individual effect which multiplies equation
	vector <unsigned int> fix_eff_mult; // References fixed effects whihc multiply equation
};

struct MarkovEqnVariation {        // Stores variation in Markov equations for a species
	bool time_vari;                  // Determines if equation has time variation
	vector <IndTransRef> ind_tra;    // The individual transtions associated with equation 
	double dt;                      // The timestep used for divisions (set to UNSET for no time)
	
	// Used in simulation
	double value;                    // The value
	double indfac_sum;               // The sum of the potential individual effect acting on equation
	
	// Used in inference
	vector <MarkovEqnDiv> div;       // Information for each time division
};

struct EventData {                 // Stores individual data
	EventData(){ move_c = UNSET; cl = UNSET; tr = UNSET; tdiv = UNSET;}
	EventType type;                  // The type of the event
	unsigned int move_c;             // The compartment in classification (when moving)
	unsigned int cl;                 // The classification 
	unsigned int tr;                 // The transition which occurs
	double tdiv;                     // The time at which transition occurs
};

struct ObsData {
	ObsData(){ so = UNSET; ref = UNSET; cl = UNSET; test_res = UNSET; tdiv = UNSET;}
	
	ObsType type;                    // The type of the observation
	unsigned int so;                 // References the data source
	unsigned int ref;                // References population or trans pop data
	unsigned int c_exact;            // The value of c (if known precisely)
	vector <EquationInfo> c_obs_prob_eqn; // Equation giving the probability of being in compartment
	vector <bool> eqn_zero;          // Set to true if the equation is zero
	bool not_alive;                  // Set to true if compartmental observation is not alive
	vector <unsigned int> obs_eqn_ref;// References the observation equation
	unsigned int cl;                 // The classification
	EquationInfo Se_eqn;             // Stores the Se observation equation
	EquationInfo Sp_eqn;             // Stores the Sp observation equation
	unsigned int Se_obs_eqn_ref;     // References Se in obs_eqn
	unsigned int Sp_obs_eqn_ref;     // References Sp in obs_eqn
	bool test_res;                   // The result of a diagnostic test
	double tdiv;                     // The time at which observation occurs
	bool time_vari;                  // Determines if time variation
};

struct EnterCla {                  // Sets information about individuals entering a cla
	unsigned int c_set;              // The value of c (if set)
	vector <EquationInfo> eqn;       // Provides probability of different global c 
	vector <unsigned int> obs_eqn_ref;// References the observation equation
};

struct Enter {                     // Stores a probability distribution for entering ind
	string name;                     // Reference name
	double tdiv;                     // The time the individual enters
	unsigned int c_set;              // If the global c value is set (otherwise eqn is used)
	vector <EnterCla> cla;           // Information about probability entering different cla
	bool set;                        // Determines if set
};

struct IndData {                   // Stores data on individuals
	string name;                     // The name of the individual
	vector <EventData> ev;           // Any events on individual (enter, leave, move)
	vector <ObsData> obs;            // Any observed events on individual
	vector < vector <TrigEventRef> > cl_trig_ev_ref; // Combines ev and obs for classification (used for sampling proposals)
	vector < vector <TrigEventRef> > fixed_trans_ev; // Fixed transition events (observed and move)
	vector <TrigEventRef> trig_ev_ref;// Combines ev and obs (used for simulation proposals)
	vector <bool> sample_needed;     // Determines if individual sampling needed [cl]
	bool simulation_needed;          // Determines if individual resimulation
	bool move_needed;                // Determines if individual move events
	unsigned int enter_ref;          // References enter compartment probability
	bool init_c_set;                 // Determines if the init individual compartment set
	double tdivmin;                     // The minimum time for data
	double tdivmax;                     // The maximim time for data
};

struct WarnData {                  // Stores warning messages
	string te;                       // The text for the warning
	unsigned int line_num;           // The line number in the input file
};

struct SpeciesSimp {               // A simplified version of species to pass to equation
	SpeciesSimp(string name_, const vector <Classification> &cla, vector <IndEffect> &ind_effect, vector <FixedEffect> &fix_effect, const vector <unsigned int> &comp_mult, const vector <CompGlobal> &comp_gl, const vector <TransGlobal> &tra_gl, bool trans_tree_) : cla(cla), ind_effect(ind_effect), fix_effect(fix_effect), comp_mult(comp_mult), comp_gl(comp_gl), tra_gl(tra_gl)
	{
		name = name_;
		trans_tree = trans_tree_;
	}

	string name;                     // The name of the species
	const vector <Classification> &cla;// Stores classification within the species
	vector <IndEffect> &ind_effect;  // Any individual effects within species
	vector <FixedEffect> &fix_effect;// Any fixed effects within species
	const vector <unsigned int> &comp_mult; // Used to convert from classification definition to global
	const vector <CompGlobal> &comp_gl; // The global compartments
	const vector <TransGlobal> &tra_gl; // The global transition
	bool trans_tree;                 // Determines if the transmission tree is turned on
};

struct ParamTag {                  // Used to check tags on parameters specified correctly
	string val;                      // The value on the tag
	string tag;                      // The tag name
};


struct Derive {                    // Stores derived parameters
	string name;                     // The name of the derived quantity
	string full_name;                // Full name including dependency        
	bool time_dep;                   // Time dependency
	vector <Dependency> dep;         // Dependency for derived quantity (excluding time)
	vector <EquationInfo> eq;        // The equations for each quantity
	DerFunc func;                    // Stores information about derived function
	unsigned int line_num;           // Stores the line number the equation comes from
};

struct DeriveOutput {              // Calculates any derived outputs
	vector <string> value_str;       // Value  
};

struct Change {                    // Stores changes to state
	double num;                      // Value of change
	vector <double> vec;             // Vector of value changes
};

struct Like {                      // Stores a likelihood (or change in likelihood)
	Like(){ 
		init_cond = 0;
		init_cond_prior = 0;
		obs = 0;
		prior = 0; 
		spline_prior = 0; 
		dist = 0;
		markov = 0; 
		nm_trans = 0;
		genetic_process = 0;
		genetic_obs = 0;
		ie = 0;
	}
	double init_cond;                // The likelhood for initial conditions
	double init_cond_prior;          // The prior for indivdual conditions
	double obs;                      // The observation probability
	double prior;                    // The prior
	double spline_prior;             // The spline prior
	double dist;                     // The probability of internal distributions
	double markov;                   // The likelihood of Markovian transitions
	double nm_trans;                 // The likelihood of non-Markovian transitions
	double genetic_process;          // The likelihood for genetic processing
	double genetic_obs;              // The likelihood for genetic observations
	double ie;                       // The likelihood of individual effects
};

struct IndRef {                    // An individual reference
	unsigned int p;                  // Species number
	unsigned int i;                  // Individual number
};

struct IndSimProb {                // User to workout if simulation proposal should be done
	double ntr;                      // Number of proposals tried (with time fading)
	double nac;                      // Number of proposals accepted (with time fading)
	double prob;                     // Probability of trying 
	bool done; 
};

struct IndSampler {                // Used as an individual sampler
	unsigned int ntr;                // Number of tries
	unsigned int nac;                // Number accepted
	unsigned int nfa;                // Number failed
	double si;                       // Size of change
};

struct PopFilter {                 // Stores filter used for population data
	string name;                     // The name of the filter
	vector <EquationInfo> comp_prob_eqn; // Provides probability of in different global compartments
	vector <bool> eqn_zero;          // Set to true if the equation is zero
	vector <unsigned int> comp_obs_mod_ref; // References the observation model equation
	vector <unsigned int> c_nonzero; // Lists all nonzero values
	bool time_vari;                  // Determines if time varying
};

struct PopData {                   // Stores information about a population measurement
	unsigned int so;                 // References the source 
	double tdiv;                     // The time of the population measurement
	ObsModelVariety type;            // The type of observation
	double value;                    // The value of the measurement
	double obs_mod_val;              // The sd or p in the observation model
	unsigned int ref;                // References pop_filter
	vector <unsigned int> comp_obs_mod_ref; // References the obs model eqn (for time vari)
	bool time_vari;                  // Determines if time varying
};

struct PopTransFilter {            // Stores filter used for population transition data
	string name;                     // The name of the filter
	vector <EquationInfo> trans_prob_eqn; // Gives prob eqn for global transitions 
	vector <bool> eqn_zero;          // Set to true if the equation is zero
	vector <unsigned int> trans_obs_mod_ref; // References the observation model equation
	vector <unsigned int> tr_nonzero; // Lists all nonzero values
	bool time_vari;                  // Determines if time varying
};

struct PopTransData {              // Stores information about a population transitions
	unsigned int so;                 // References the source 
	double tdivmin;                  // The start time
	double tdivmax;                  // The end time
	ObsModelVariety type;            // The type of observation
	double value;                    // The value of the measurement
	double obs_mod_val;              // The SD or p in the observation model
	unsigned int ref;                // References pop_trans_filter
	vector <unsigned int> trans_obs_mod_ref;  // References the observation model equation
	bool time_vari;                  // Determines if time varying
};

struct TransDiagSpecies {          // Stores expected number of transitions based on model
	vector < vector <double> > exp_num;// [tgl][ti]
	unsigned int n;                  // The number which have been added
};

struct ParticleSpecies {           // Store the species state is a particle     
	InitCondValue init_cond_val;     // The initial state
	vector < vector <double> > trans_num; // Transition numbers (for population-based)
	vector <Individual> individual;  // Individuals (for individual-based)
	unsigned int nindividual;        // The number of individuals (used for trace plots)
	vector < vector <double> > exp_num;// Expected number of transitions [trg][ti]
	vector < vector <unsigned int> > cum_prob_dist;// Cumulative probability distribution [trg][b]
};

struct TransTreeStats {            // Stores statistics about transmission tree
	TransTreeStats(){ N_origin = UNSET; N_inf = UNSET; N_mut_tree = UNSET; N_mut_origin = UNSET; N_unobs = UNSET; t_root = UNSET;}
	
	unsigned int N_origin;           // The number of origins of infection entering
	unsigned int N_inf;              // Number of infected individuals
	unsigned int N_mut_tree;         // Number of mutations on the tree
	unsigned int N_mut_origin;       // Number of mutations for origins
	unsigned int N_unobs;            // Number of unobserved genetic observations
	unsigned int t_root;             // The time infection first enters 
};

struct SourceTransRate {           // Stores rates at which individuals enter system
	string rate;                     // The rate of entry
	vector <unsigned int> tr_gl;     // The global transitions to which source applies
};

struct FracSort {                  // Sorts individuals based on fraction (used for MVN sampling)
	unsigned int i;                  // Index
	double frac;                     // Fraction
};

struct ListMap {                   // Lists used for adding and removing unobserved individuals
	vector <unsigned int> list;      // List of individuals
	vector <unsigned int> map;       // References list
	unsigned int min;                // The minimum number on list
	
	ListMap(unsigned int min, unsigned int max);
	void add(unsigned int i);
	void remove(unsigned int i);
	void add_hole(unsigned int i);
	void remove_hole(unsigned int i);
	void check(unsigned int num);
};

struct ObsTrans {                  // Information about observed individual transitions
	string name;                     // The name of the observed transition
	ObsType type;                    // Determines if trans, source or sink
	vector <EquationInfo> tra_prob_eqn; // Gives transition probability [trg]
	vector <unsigned int> obs_eqn_ref;// References the observation equation (if no time vari)
	bool time_vari;                   // Determines if time variation
	unsigned int ti_min;             // The minimum time division for observation
	unsigned int ti_max;             // The maximum time division for observation
	vector <bool> is_one;            // See if the transition is observed perfectly
	bool single_trans;               // Sets if just a single transition is observed
};

struct DataTransPos {              // A transition possibility
	unsigned int trg;                // Global transition
	double prob;                     // Probability
	double probsum;                  // Cumulative probability
};

struct Bracket {                   // Stores bracket information
	string ch;                       // Character
	unsigned int i;                  // Index
};

struct NodeRef {                   // References an infection node
	unsigned int node;               // Node number
	unsigned int index;              // Index of event
};

struct InfNodeAlter {              // Determines an alternation in an infnode
	InfNodeAlter(){ possible = true; unchanged = false; e_add = UNSET;}

	unsigned e_add;                  // The event number for infection on new sequence
	double tdiv_start;               // The time of the addition
	double tdiv_rec;                 // The recovery time of the added individual 
	vector <unsigned int> obs_add_begin; // Genetic observations added at the beginning
	vector <unsigned int> obs_add_end;   // Genetic observations added at the beginning
	bool possible;                   // Set if change is possible
	bool unchanged;                  // Set if alteration leaves unchanged
};

struct InfEvent {                  // Infection event (references genetic observation or infection node)
	double tdiv;                     // The event time
	InfEventType type;               // INFECT_OTHER or GENETIC_OBS 
	unsigned int index;              // Node infected or number of observation
	unsigned int mut_num;            // The number of mutations immediately preceeding
};

struct InfOrigin {                 // An origin for a genetic strain 
	unsigned int node;               // Reference infection node
	unsigned int mut_num;            // The number of mutation for the source 
};

struct GenDifCha {                 // Stores a change in the genetic difference matrix 
	unsigned int ob1;                // First observation
	unsigned int ob2;                // Second observation
	int ch;                          // Change in matrix
};

struct ObsMut {                    // Used to store genetic observation (used in genetic props)
	unsigned int obs;                // Observation number
	unsigned int mut_num;            // Mutation number
};

struct InfNode {                   // Stores information about an infection events
	double tdiv_start;                  // The time the infection start
	double tdiv_rec;                    // The time the individual recovers
	unsigned int p;                  // The species of infected individual
	unsigned int i;                  // Individual number of infected individual
	unsigned int e;                  // The event number
	NodeRef from;                    // References node generating infection 
	                                 // (if entering from.node = ENTER_INF, if from origin from.node = OUTSIDE_INF, from.index gives inf_origin number)
	vector <InfEvent> inf_ev;        // References nodes which are infected + genetic obs
};

struct GenChange {                 // Stores information about a change to  
	GenChange(GenChaType type_) {
		probif = 0; probfi = 0; mut_num = UNSET; 
		dlike_genetic_process = 0; dlike_genetic_obs = 0;
		type = type_; 
	}
	void update_like_ch(Like &like_ch, double &dprob);
	
	GenChaType type;                 // The type of change
	unsigned int n;                  // The node being changed
	IndInfFrom iif_add;              // Stores iif_add
	NodeRef nr_add;                  // The node reference to which the node is being added
	InfNodeAlter node_alter;         // Store information about the alternation to the node
	RegraftType regraft_type;        // Determiones how the regrafting is done
	InfNode node_add;                // Stores information about node being added
 	double probif;                   // Probability of proposing
	double probfi;                   // Reverse probability of proposing
	unsigned int mut_num;            // Stores the mutation number of newly inserted connection
	vector <GenDifCha> gen_dif_cha;  // Change to observation matrix
	double dlike_genetic_process;    // Change in likelihood for genetic process
	double dlike_genetic_obs;        // Change in likelijood for genetic observation
};

struct ObsGeneticData {            // Store individual genetic data
	string name;                     // Name (Used for matrix data)
	unsigned int p;                  // The species
	unsigned int i;                  // The individual
	unsigned int r;                  // Allows for reordering
	double tdiv;                     // The time of observation
	vector <GenChar> snp;            // The snps (used for snp data)
};

struct ObsGenRef{                  // References a genetic observation             
	double tdiv;                     // The time of observation
	unsigned int m;                  // The observation number
};

struct GeneticData {               // Stores genetic data
	bool on;                         // Determines if genetic data is turned on
	EquationInfo mut_rate;           // The equation for the mutation rate
	EquationInfo seq_var;            // The equation for the sequence variation
	unsigned int nSNP;               // The number of SNPs
	double obs_scale;                // Scales the likelihood function
	vector <ObsGeneticData> obs;     // Observed individual genetic data 
	vector < vector < vector <ObsGenRef> > > ind_gen_obs; // The genetic observations on an individual [p][i][#]
	vector < vector <unsigned int> > gen_dif; // The genetic difference between sequences 
};

struct GeneticDataValue {          // Stores state values for genetic data
	double mut_rate;                 // The mutation rate
	double seq_var;                  // The sequence variation
	unsigned int nobs_not_infected;  // Number of genetic observations when not infected
	vector <NodeRef> obs_node_ref;   // Reference what nodes the observations are on
	vector <InfOrigin> inf_origin;   // The number of origins of virus
	vector < vector <unsigned int> > gen_dif; // The state genetic difference 
	vector <InfNode> inf_node;       // Infection nodes
};

struct Particle {                  // Stores information from state
	vector <double> param_val_prop;  // The parameter value (without const or reparameterisation)
	vector <ParticleSpecies> species;// Species state data
	
	vector <DeriveOutput> dir_out;   // Derived outputs
	Like like;                       // Stores the overall likelihood and prior
	TransTreeStats trans_tree_stats; // Statistics about the transmission tree
	vector <InfOrigin> inf_origin;   // The number of origins of virus
	vector <InfNode> inf_node;       // The infection nodes
	
	double w;                        // Particle weight (used in ABC-SMC)
	
	unsigned int s;                  // The sample number
	
	unsigned int chain;              // The chain number
};

struct BurnInfo {                  // Information about burnin phase                       
	BurnInfo(){ phi = 0; L_samp_sum = 0; L_samp_sum2 = 0; start = 0;} 

	bool on;                         // Determines if burning in
	double fac;                      // Allows for proposals to change quickly in initial burnin 
	vector <double> L_samp;          // Samples of the likelihood
	double L_samp_sum, L_samp_sum2;  // Keeps a running sum
	unsigned int start;              // First sample used in sum
	unsigned int range;              // The range over which varaince is estimated
	double alpha, epsilon;           // Parameters used for log-auto model
	double power;                    // The power of the annealing (dynamically updated)
	double phi;                      // The overall control parameter
	double init_cond;                // The inv temperature for likelihoods
	double init_cond_prior;               
	double obs; 
	double prior;
	double spline_prior;
	double dist;
	double markov;
	double nm_trans;
	double genetic_process;
	double genetic_obs;
	double ie;
	bool dprob_suppress;             // Supresses dprob term in MH proposal for initial samples 
	
	unsigned prop_join_step;         // How often to check if proposals should be jointed 
	
	void add_L(double L);
	void setup(unsigned int s, unsigned int &nburnin,  unsigned int &nsample, const Details &details);
	void pas_setup(unsigned int s, unsigned int g, unsigned int gen_update, double phi_);
	void pas_setup_run(unsigned int s, unsigned int &nburnin);
	void set_phi();
};

struct Fraction {                  // Used to store a fraction (used in initial state)
	unsigned int c;                  // The compartment number
	double fr;                       // The fraction
};

struct LocalChange {               // Stores a local event change
	LocalChangeType type;            // Type of event change
	unsigned int tr;                 // Transition number
	int ti;                          // Time division
};

struct AddRemLocal {               // List of transitions which leaves populations unchanged
	vector <unsigned int> tr_list;   // List of global transition
};

struct KernelAdd {                 // Kernel (used for local population proposals)
	int dti;                         // Change in time division
	double w;                        // The weight
};

struct Samp {                      // Used to sample from a 2D array (used in local pop props)
	unsigned int i;                  // i position in 2D array 
	unsigned int j;                  // i position in 2D array 
};

struct Sampler {                   // Sampler to sample from a vector of possibilities
	unsigned int N;                  // The number of compartments
	vector <double> num;             // Probability of sampling transitions at different times [tr][ti]
	vector <double> num_sum;         // Cumulative probability
 	double S;                        // Total sum
	
	void init(unsigned int _N);
	void setup();
	unsigned int sample(double &dprob) const;
	void update(unsigned int i, int sign);
	void check();
};

struct Sampler2D {                 // Used to sample from a two dimensional quantity
	unsigned int N;                  // The number in first dimension 
	unsigned int N2;                 // The number in second dimension
	vector <KernelAdd> kernel;       // Sets up a kernel when adding to sampler
	vector < vector <double> > num;  // Probability of sampling transitions at different times [tr][ti]
	vector < vector <double> > num_sum; // Cululative sum  
	vector <double> marg_sum;        // Marginal sum
	double S;                        // Total sum
	
	void init(unsigned int _N, unsigned int _N2);
	void setup();
	Samp sample(double &dprob) const;
	void update(unsigned int i, unsigned int j, int sign);
	void check();
	void print() const;
};


struct LocalSampler {              // Local sampler (for population-based models)
	vector <unsigned int> tr_list;   // The transitions used for add/remove local events 
	Sampler tr_ic;                   // The probability of sampling transitions to change IC
	Sampler comp_ic;                 // The probability of sampling compartment to change IC
	Sampler2D tr_samp;               // The transition sampler for add/rem single events
	double win;                      // The window size used for update
};

struct PropSpeed {                 // Stores proposals speeds
	unsigned int i;                  // The proposal number
	double time_per_prop;            // The time of a single proposal
};

struct IndEffGroupRef {            // References ind eff group element (for ie covar propoasls)
	unsigned int ieg;                // References an individual effect group
	unsigned int i;                  // i position in the omega matrix
	unsigned int j;                  // j position in 2D array 
};

struct NMupdate {                  // Stores information about a change in non-Markovian trans
	NMUpdateType type;               // Type of change
	unsigned int i;                  // Individual number
	unsigned int e;                  // Event number
	double Li;                       // nm transition likelihood
	double Li_bp;                    // Branching likelihood
};

struct ParEventJointInfo {         // Stores information about joint pa
	ParEventJointInfo(){ possible = true; p = UNSET;}
	bool possible;                   // Determines if possible
	unsigned int p;                  // Stores the species
	vector <unsigned int> tr_list;   // The transitions which parameter affects
};

struct ParEventJointProp {         // Stores information about a joint para event prop
	unsigned int th;                 // The param being changed 
	vector <unsigned int> tr_list;   // The corresponding transitions which change
	Direction dir;                   // The direction the change is done in
};

struct BackInd {                   // Stores information so ind. can be restored (if rejected)
	unsigned int i;                  // Individual undergoing change
	vector <Event> ev;               // The event sequence
	vector <IncompNMTransRef> incomp_ref; // Incomplete transition 
	vector <unsigned int> identical_old;  // Store events which are identical (e on old sequence)
	vector <unsigned int> identical_new;  // Store events which are identical (e on new sequence)
	vector <unsigned int> remove_old;// Events removed from old sequence
	vector <unsigned int> add_new;   // Events added to new sequence
};

struct Back {                      // Stores information so state can be restore (if rejected)
	Back(BackType type_, unsigned int i_, unsigned int index_, double value_){
		type = type_;
		i = i_;
		index = index_;
		value = value_;
	}

	Back(BackType type_, unsigned int i_, unsigned int index_, double value_, const vector<double> &vec_){
		type = type_;
		i = i_;
		index = index_;
		value = value_;
		vec = vec_;
	}

	Back(BackType type_, unsigned int i_, unsigned int index_, const vector<double> &vec_){
		type = type_;
		i = i_;
		index = index_;
		value = UNSET;
		vec = vec_;
	}

	Back(BackType type_, unsigned int i_, const vector<double> &vec_){
		type = type_;
		i = i_;
		vec = vec_;
		value = UNSET;
		index = UNSET;
	}

	Back(BackType type_, const vector<double> &vec_){
		type = type_;
		vec = vec_;
		i = UNSET;
		value = UNSET;
		index = UNSET;
	}

	Back(BackType type_, unsigned int i_, double value_){
		type = type_;
		i = i_;
		value = value_;
		index = UNSET;
	}

	BackType type;                   // The type of change
	unsigned int i;                  // Individual / index
	unsigned int index;              // Index
 	double value;                    // Value of change
 	vector <double> vec;             // Vector of values
};

struct BackPop {                   // Stores information so pop can be restore (if rejected)
	BackPop(BackPopType type_, unsigned int i_, unsigned int j_, unsigned int k_, double value_){
		type = type_;
		i = i_;
		j = j_;
		k = k_;
		value = value_;
	}

	BackPopType type;                // The type of change
	unsigned int i;                  // Index i
	unsigned int j;                  // Index j
	unsigned int k;                  // Index k
	double value;                    // Value for change
	vector <double> vec;             // Vector of values
};

struct PopUpdate {                 // Stores update in population
	unsigned int ti;                 // The initial time division 
	unsigned int ti_next;            // The final time division
	unsigned int c_old;              // The old compartment
	unsigned int c_new;              // The new compartment
};

struct EventCl {                   // Event for constructing timeline in a classification
 	TrSwap trs;                      // Potential swap
	unsigned int c_bef;              // The compartment before event
	unsigned int c_aft;              // The compartment after event
	unsigned int e;                  // References event
	double tdiv;                     // Time of event
};

struct LocalIndChange {            // Stores potential local changes to individual
	unsigned int i;                  // Individual
	unsigned int cl;                 // Classification
	double tdivmin;                  // Gives the minimum time for the local proposal
	double tdivmax;                  // Gives the maximum time for the local proposal
	vector <unsigned int> remove_e;  // Which events to remove
	vector <TrSwap> tr_swap;         // Transitions which need to be added    
	unsigned int swap_rep_ref;       // References swap_rep
	unsigned ref;                    // References lc_ref 
};

struct TraRate {                   // Stores rates for nm and Markovian transition 
 	bool markov;                     // Set if transition is Markovian
	bool bp_on;                      // Set if there is a branching prob
	bool ind_variation;              // Determines if there is individual variation
	vector <IndFacRate> ind_fac_rate;// Determines how rate is modified by ind effects
	vector <double> value_raw;       // Mean rate for transition multiplied by bp [ti]
	vector <double> bp_raw;          // Mean rate for transition [ti]
	vector <double> value;           // Mean rate for transition multiplied by bp [ti]
	vector <double> integral;        // Integrated rate
	vector <double> bp;              // Mean rate for transition [ti]
};

struct CompR {                     // Compartment rate
	CompRUpType update;              // The type of update
	vector <unsigned int> list;      // Stores Markovian or NM reference
	vector <double> integral;           // Sum of rate leaving compartment
};

struct RatePosteriorMean {         // Stores posterior mean of rates (for local proposals)
	RatePosteriorMean();

	bool first;                      // Set if first time to be updated
	vector <unsigned int> tramean_ref;         // References where to find transition rate [trg]
	vector <vector <unsigned int> > compR_ref; // References transition linked to comp [cl][cgl]
 	
	// Here markovian and non-Markovian transitions are stacked for TraMean
	vector <TraRate> tra_rate;       // Stores rates going down transitions
	vector <CompR> compR;            // Different combination of tra leaving
	
	double get_value(const vector <double> &vec, double v) const;
	void update_prob_try_Rint(double t1, double t2, unsigned int tr_gl, unsigned int cl, unsigned int c, double &prob_try, double &Rint, const Individual &ind, Modify mod) const;
};

struct SourceSamp {                // Used to generate a sample for the source
	unsigned int tr_gl;              // The golbal transition
	double prob_sum;                 // The cumulative probability
};

struct AllInd {                    // References all observed individuals
	unsigned int p;                  // Species number
	unsigned int i;                  // Individual number
};

struct SampleSpecies {             // Used to store information about a posterior state sample
	Table ind_tab;                   // Individuals
	Table trans_num_tab;             // Transition numbers 
	Table cpop_init_tab;             // Initial population
};

struct Sample {                    // Stores posterior sample (for PPC)
	unsigned int num;                // The sample number    
	vector < vector <double> > param_value; // The parameter sample
	vector <SampleSpecies> species;  // The state sample
};

struct TransRef {                  // Referenes transitions
	unsigned int cl;                 // The classificatiomn
	unsigned int tr;                 // The transition number
};

struct InfCause {                  // Reference how one indididual affect another
	InfCause(){ eq_ref=UNSET; pref = UNSET; po = UNSET;};
	unsigned int eq_ref;             // References the equation
	unsigned int pref;               // The pref
	unsigned int po;                 // The population
};

struct BICITag {                   // Gets a tag when running bici
	TagType type;                    // The type of tag
	unsigned int value;              // The value of the tag
	string te;                       // Text value
	bool processed;                  // Determines if the tag has been processed
};

struct PartReorder {               // Used for reorder particles based on likelihood
	double L;                        // Likelihood
	unsigned int i;                  // Particle number
};

struct Diagnostic {                // Stores diagnostic information about chain
	unsigned int ch;                 // Chain number
	string te;                       // Text
};

struct Stat {                      // Statistics
	double mean;                     // The mean
	double sd;                       // The standard deviation            
	double CImin;                    // Credible interval
	double CImax;
};

struct MeanSD {                    // Stores the mean and stadard deviation for a distribution
	double mean;
	double sd;
};

struct TRange {                    // Stores a time range
	double tdivmin;
	double tdivmax;
};

struct ListMove {                  // Stores a list and a corresponding move
	vector <unsigned int> list;
	vector <bool> move;
};

struct EventMove {                 // Used for reordering event for multi-event proposal
	double tdiv;
	unsigned int e;
	bool move;
};	

struct PropTime {                  // Used for orthering proposal times
	unsigned int i;
	long time;
};


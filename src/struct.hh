/// Stores commonly used data structures

#ifndef BICI__STRUCT_HH
#define BICI__STRUCT_HH

#include <vector>

using namespace std;

#include "const.hh"
//#include "equation.hh"

struct ParamRef                    // Used to reference parameter
{
	bool found;                      // Determines if parameter has been found
	unsigned int th;                 // The number of the parameter
	unsigned int index;              // The index where to find 
};
	
struct AffectLike {                // Determines how a parameter affects likelihoods
	AffectType type;                 // The type of likelihood being affected
	unsigned int num;                // The reference number
	unsigned int num2;               // Another reference number
	vector <bool> map;               // Maps all time elements which are affected
	vector <unsigned int> list;      // Lists all map which is true
	unsigned int order_num;          // Used to order terms by list priority
};

struct PopulationTerm {            // A term in the population
	unsigned int c;                  // Compartment
	double w;                        // Weight
};

struct PopRef {                    // References a population
	unsigned int po;                 // The population number
	unsigned int index;              // The index
};

struct PopMarkovEqnRef {           // References a markov equation from a population
	unsigned int p;                  // The species
	unsigned int e;                  // The Markov eqn number
};

struct PopTransRef {               // References transition from a population
	unsigned int p;                  // The species
	unsigned int tr;                 // The transition number
};
	
struct Population {                // Stores a population (used in an equation)
	string name;                     // The name of the population
	unsigned int sp_p;               // The species number
	vector <unsigned int> ind_eff_mult; // Individual effect which multiplies weight (if appropriate) 
	vector <unsigned int> fix_eff_mult; // Any fixed effects which modify weight
	vector <PopulationTerm> term;    // The global compartments which contribute to the popualtion
	vector <PopMarkovEqnRef> markov_eqn_ref;// References markov equations which include population
	vector <PopTransRef> trans_ref;  // References transitions which include population
};

struct ParamProp {                       // Stores properties of a parameter
	string name;                           // The name (without the 
	vector <string> dep;                   // String form of dependency
	vector <string> dep_with_prime;        // String form of dependency (with primes)
	bool time_dep;                         // Time dependency
	bool age_dep;                          // Age dependency
	bool prime;                            // Set to true of anu dependencies are primes
};
	
struct SplineInfo {                      // Stores information about a parameter spline
	bool on;                               // Detemines if spline is switched on
	vector <double> knot_times;            // Times of knots (if time or age variation)
	bool smooth;                           // Determines if smoothing is added
	SmoothType smooth_type;                // The type of smoothing (normal or log-normal)
	double smooth_value;                   // The value used for the smoothing
};

struct SplineDiv {                       // A div in spline is calculated from a weighting of params
	unsigned int th1;                      // The reference for the parameter 1
	unsigned int th2;                      // The reference for the parameter 2
	double f;                              // The factor 1=th1  0=th2
};

struct Spline {                          // Stores an individual spline
	string name;
	bool constant;                         // Determines if spline is constant
	unsigned int th;                       // The parameter
	unsigned int index;                    // The index giving the dependency
	vector <unsigned int> param_ref;       // References the parameters along the spline
	SplineInfo info;                       // Stores information about spline (copied from par)
	vector <SplineDiv> div;                // A division within the spline
	vector <double> const_val;             // If spline is constant then works out what value it should have
};

struct SplineValue {                     // Stores values for a spline
	vector <double> val;                   // The value down the spline
};


struct CompPos {                         // Used to go through comparmtental possibilities
	vector <string> list;
	unsigned int index;
}; 

struct DepInfo {                         // Stores dependency info froman equation
	unsigned int iend;                     // The end position for the dependency
	vector <string> spl;                   // Splits the dependency into individual elements
	Result result;                         // Stores any warning message
};

struct SwapResult {                      // Stores the result of doing a substitution of indeces for values
	string warn;                           // Stores the warning message0
	vector <bool> done;                    // Registers true if a swap has been performed
};

struct DepConv {                         // Used to convert dependencies from index to value
	string before;
	string after;
};

struct EquationInfo {                    // Stores information about an equation (prior to equation being evaluated)
	EquationInfo(){ te = ""; p = UNSET; cl = UNSET;}
	string te_raw;                         // The rate text value (as seen in the input file)
	string te;                             // The text value after processing
	unsigned int p, cl;                    // The species and classification for equation (if appropriate)
	EqnType type;                          // The type of the equation
	double value;                          // The numerical value (if appropriate)
	unsigned int eq_ref;                   // The equation (if appropriate)
	unsigned int line_num;                 // Stores the line number the equation comes from
};

struct IEname {                          // Stores individual effect name
	string name;                           // The name of the individual effect
	unsigned int index;                    // Indexes the ie in the species
};

struct Amatrix {                         // Stores an A matrix
	bool set;                              // Determines if the A matrix has been set
	vector < vector <double> > value;      // The values in the matrix
	vector <string> ind_list;              // The names of individuals in the system
};

struct IndEffect {                       // Store infortmation about an individual effect (in species)
	string name;                           // The name of the individual effect
	unsigned int index;                    // References ind_eff_group
	unsigned int num;
	vector <unsigned int> markov_eqn_ref;  // References any markov equations the ie appears in
	vector <unsigned int> nm_trans_ref;    // Indexes nm_trans involving ies
	vector <unsigned int> pop_ref;         // Populations that ie appear in
	unsigned int line_num;                 // The line number in input file 
};

struct IEgroup {                         // Stores information about a group of individual effects
	vector <IEname> list;                  // The list of individual effect names
	Amatrix A_matrix;                      // The A matrix
	
	vector < vector <unsigned int> > omega;// Matrix of parameter values 
	unsigned int line_num;                 // Stores the line number the equation comes from
	
	vector <unsigned int> markov_eqn_ref;  // Indexes any markov equations involving ies
	vector <unsigned int> nm_trans_ref;    // Indexes nm_trans involving ies
};

struct NZElement {                       // Stores non-zero elements of a matrix
	unsigned int i;                        // x position
	double value;                          // Value
};

struct IEGroupSampler {                  // Stores a sample for an individual effect group
	vector < vector <double> > omega;      // The covariance matrix omega
	vector < vector <double> > omega_Z;    // The cholesky matrix for omega (if doing simulation)
	vector < vector <double> > omega_inv;  // The inverse matrix for omega (if doing inference)
	vector < vector <double> > A;          // The A relationship matrix
	vector < vector <double> > A_Z;        // The cholesky matrix for A
	vector < vector <double> > A_inv;      // The inverse relationship matrix
	vector <double> A_inv_diag;             // The diagonal elements of Ainv
	vector < vector <NZElement> > A_inv_nonzero;// The non-zero  elements
};

struct Dependency {                      // A dependency in the model
	string index;                          // The index in the model
	string index_with_prime;               // The index in the model including any prime
	vector <string> list;                  // List of possible values (either compartments or times along spline)
	unsigned int mult;                     // Multiplicative factor when accessing an element
};

struct Prior {                           // Defines a parameter prior
	Prior(){ error = "";}
	PriorPos type;                         // The type of prior
	vector <EquationInfo> dist_param;      // A list of equations to specify prior distribution parameters
	string error;													 // Stores error message if declared incorrectly
};

struct ParamVecEle                 // Stores information about an element in param vec
{
	string name;                     // The name of the parameter
	unsigned int th;                 // The number of the parameter
	unsigned int index;              // The index where to find 
	Prior prior;                     // The prior distribution (copied from param)
	ParamVariety variety;            // The parameter variety (copied from param)
	
	vector <AffectLike> affect_like; // Determines how parameter affects likelihoods
};


struct Param {                           // Stores a model parameter
	string full_name;                      // The full name (as defined in the input file)
	string name;                           // The name after removing dependencies
	vector <Dependency> dep;               // Dependency for parameter
	bool time_dep;                         // Time dependency
	bool age_dep;                          // Age dependency
	SplineInfo spline_info;                // Stores information about spline (if used)
	bool auto_value;                       // Detemines if value is automatically calculated (e.g. for branching probabilities)
	
	bool sim_sample;                       // Set to true if value sampled from the distribution (simulation)
	
	ParamVariety variety;                  // Determines the variety of parameter	

	unsigned int N;                        // The number of elements in the parameter
	
	vector <EquationInfo> value;           // Different values (if parameter has dependency)
	
	vector <Prior> prior;                  // Different priors for different elements of the parameter
	
	vector < vector <ParamRef> > parent;   // Parameters which parameter is dependent on
	vector < vector <ParamRef> > child;    // Parameters which parameter is dependent on
	
	vector <unsigned int> param_vec_ref;   // References the value on param_vec_ref
	
	unsigned int line_num;                 // The line num in input file when parameter is defined
	bool used;                             // Determines if param used in model
};

struct ErrorMess {                       // Stores error messages when importing 
	unsigned int line_num;                 // The number of the line in the input file
	string error;                          // The error message
	ErrorType type;                        // The type of the error (fatal or just a warning)
};

struct LatLng {                          // Denote a point with latitude and longitude
	double lat;
	double lng;
};

struct Polygon {                         // Denotes a polygon for a boundary
	vector <LatLng> point;
};

struct Boundary {
	vector <string> name;                  // Potential names for boundary (taken from "properties"
	vector <Polygon> polygon;              // Boundary made up of polygons
};

struct GeoJSON {
	string file;
	vector <Boundary> boundary;
};

struct Details {                         // Stores details
	double t_start;                        // The start and end time for simulation
	double t_end;
	double dt;                             // The time step (used in tau leaping)
 	unsigned int number;                   // Number of simulations to be performed
	unsigned int thinparam;                // The thinning factor for parameters (MCMC)
	unsigned int thinstate;                // The thinning factor for parameters (MCMC)
	unsigned int sample;                   // Number of samples to be done (MCMC)
	double accfrac;                        // The acceptance fraction (ABC)
	unsigned int numgen;                   // The number of generations (ABC-SMC)
	double kernelsize;                     // The size of the kernel distribution (ABC-SMC)
	unsigned individual_max;               // The maximum number of individuals (individual-based models)  
 	Algorithm algorithm;                   // The algorithm used to do simulation / inference
	bool stochastic;                       // Detemines if stochastic
};

struct FilterCla {                       // Used to filter a classification
	vector <bool> comp;                    // Set to true if compartment included in filter
	FilterClaType type;                    // Type of filter e.g. all compartments, specified in file or sepcified by comp 
};

struct Filter {                          // Used to filter populations etc.
	vector <FilterCla> cla;                // Stores filters for each classification
};

struct Table {                           // Loads a table
	string file; 												   // The file from which the tables was loaded
	unsigned int ncol;                     // The number of columns
	unsigned int nrow;                     // The number if rows
	vector <string> heading;               // The headings for the columns
	vector <vector <string> > ele;         // The elements of the table
	bool error;                            // Set if there is a problem loading the table
};

struct DiagTestSens {                    // Determines compartments test sensitive to 
	unsigned int cl;                       // The classification
	vector <bool> comp;                    // The compartments the test is sensitve to are set to true
};

struct ObsModel {                        // Stores the observation model assocaited with data
 	ObsModelType type;                     // The type of obsevation model
	EquationInfo Se;                            // Reference to an equation giving sensitivity (for diagnostic test data)
	EquationInfo Sp;                            // Reference to an equation giving specificity (for diagnostic test data)
	EquationInfo trapprob;                      // Reference to an equation giving trapping probability (for trapping data)
	string diag_pos;                       // The string for a positive diagnostic test (for diagnostic test data)
	string diag_neg;                       // The string for a negative diagnostic test (for diagnostic test data)
	DiagTestSens diag_test_sens;           // Determines compartments test sensitive to (for diagnostic test data)
	double percent;                        // percent value when 'error' obs model used 
	double sd;                             // sd value when 'error' obs model used 
};

struct LoadCol {                         // Used to store information about loading data columns
	LoadCol(string he, string de, ElementType ty, unsigned int cla = UNSET){ 
		heading = he; desc = de; type = ty; cl = cla;
	}
	string heading;                        // The heading for the column
	string desc;                           // A description
	ElementType type;                      // The type of element
	unsigned int cl;                       // Stores classification number (if appropriate)
};

struct CompRef {                         // References a particular compartment
	unsigned int p;
	unsigned int cl;
	unsigned int c;
	string error;
};

struct DataSource {                      // Stores information about a data source
	unsigned int index;
	
	Command cname;                         // The name of the command
	Table table;                           // The table which stores the data
	unsigned int focal_cl;                 // The focal classification (for init_pop data only)
	unsigned int p;
	unsigned int cl;                       // Classification (used in the case of MOVE_IND)
	unsigned int tr;                       // Transition number (used in TRANS_DATA, SOURCE_DATA, SINK_DATA)
	PriorPos prior_type;                   // The type of prior (used in INIT_POP_PRIOR) 
	TimeRange time_range;                  // Determines what time range is used (used in TRANS_DATA, SOURCE_DATA, SINK_DATA)
	double time_start, time_end;           // Stores the time range (used in TRANS_DATA, SOURCE_DATA, SINK_DATA)
	Filter comp_filt;                      // FIlters the population (used in POP_DATA, SETTRAP_DATA, POP_TRANS_DATA)
	ObsModel obs_model;                    // The observation model (used in POP_DATA)
	vector <bool> trans_filt;              // Shows which transitions are filtered (used in POP_TRANS_DATA)
	string SNP_root;                       // Shows the root for the SNP data (used in GENETIC_DATA)
	vector <LoadCol> load_col;             // Stores information about the loaded columns
	unsigned int line_num;                 // The line number which loads up the data source
};

struct IslandRef {                       // References an island
	unsigned int i;                        // The island
	unsigned int c;                        // The compartment within the island
};

struct Compartment {
	string name;                           // The name of the compartment
	bool markov_branch;                    // Determines if compartment is a branching if markovian
	bool branch;                           // Determines if multiple transitions leave the compartment
	bool erlang_hidden;                    // Determines if hidden Erlang compartment
	string erlang_source;                  // Determines if hidden Erlang compartment
	double x, y;                           // The position of the compartment (optional, used if distance is in the model)
	double lng, lat;                       // The position of the compartment (optional, used if distance is in the model)
	IslandRef island_ref;                  // References which island the compartment is in
};

struct TransBPRef {                      // References a transition and a branching probability
	unsigned int tr;                       // A global transition reference
	unsigned int bp_th;                    // The parameter for branching
};

struct CompGlTrans {                     // Information about transitions from a compartment
	bool branch;                           // Whether branching happens
	vector <TransBPRef> tra_ref;           // References global transitions
};

struct CompGlobal {                      // A global compartment (combines all classifications)
	string name;                           // The name of the global compartment
	bool erlang_hidden;                    // Determines if added to make Erlang distribution
	vector <unsigned int> cla_comp;        // References the classification compartment
	vector <PopRef> pop_ref;               // References population and point on list
	vector <unsigned int> pop_ref_simp;    // References just population 
	vector <CompGlTrans> trainfo;          // Information about tranitions from branching from compartment  
	vector <unsigned int> me_ref;          // References markov equations leaving compartment
};

struct Transition{
	string name;                           // The name of the transitions
	unsigned int i;                        // The initial compartment (set to SOURCE for source)
	unsigned int f;                        // The final compartment (set to SINK for sink)
	
	TransType type;                        // The type of the transition
	TransVariety variety;                  // Determines if normal, source or sink

	bool erlang_hidden;                    // Determines if hidden Erlang transitions

	bool bp_set;                           // Determines if branching probabilty has been set
	bool branch;                           // Determines if this is one of multiple branches leaving a compartment
	EquationInfo bp;                       // References branching probability equation
	
	vector <EquationInfo> dist_param;      // References a list of equations to specify transition distribution parameters
	
	unsigned int line_num;                 // Stores the import line (for diagnostic error messages)
};

struct NMTransInd {                      // Information about a NM transition on an individual
	double t;                              // The time of the transition
	unsigned int i;                        // The individual
	unsigned int cl;                       // The classification
	unsigned int index;                    // Indexes ind.cl_ev.ref
};

struct NMTrans {                         // Stores a non-Markovian transition
	TransType type;                        // The type of the transition
	unsigned int bp_eq_ref;                // The branching probability equation (if applicable)
	vector <unsigned int> dist_param_eq_ref;// Equations for distribution
};

struct TransGlobal {                     // A global transition (combines all classifications)
	string name;
	unsigned int cl;                       // The classification
	unsigned int tr;                       // The transition number in the classifciation

	unsigned int i;                        // The initial compartment (set to SOURCE for source)
	unsigned int f;                        // The final compartment (set to SINK for sink)
	
	bool erlang_hidden;                    // Determines if hidden Erlang transitions

	TransType type;                        // The type of the transition
	TransVariety variety;                  // Determines if normal, source or sink

	bool bp_set;                           // Determines if branching probabilty has been set
	bool branch;                           // Determines if this is one of multiple branches leaving a compartment
	EquationInfo bp;                       // References branching probability equation
	
	vector <EquationInfo> dist_param;      // References a list of equations to specify transition distribution parameters
	
	unsigned int nm_ref;                   // References nm_trans (if non-Markovian)
	
	unsigned int markov_eqn_ref;           // References markovian equation (if Markovian)
	
	vector <unsigned int> transform;       // Transforms transition based on a shift in c
	
	unsigned int line_num;                 // Stores the import line (for diagnostic error messages)
};

/*
struct BranchPoint {                     // Information about a branching point
	unsigned int c;                        // The compartment where the branching occurs
	vector <unsigned int> tra;             // The transitions (the last is automatically calcualted
};
*/

struct MEref {                           // References an equation and time in markov_eqn
	unsigned int e;                        // Markov_eqn ref
	unsigned int ti;                       // Time
};

struct IslandTrans {                     // Transitions within an island
	unsigned int tr;                       // References the transion
	unsigned int cf;                       // References the final compartment (in island units)
	vector < vector <MEref> > markov_eqn_ref;  // References markov equation [ti][cgl] 
};

struct IslandComp {                     // Transitions within an island
	unsigned int c;                       // References the transion
	vector <IslandTrans> leave;
};

struct Island {                          // An island is a group of compartments not connected to rest
	vector <IslandComp> comp;              // Compartment in island
};

struct Classification {                  // Stores details of a species classification
	string name;                           // The name of the classification
	
	vector <Compartment> comp;             // Stores details of compartments
	unsigned int ncomp;                    // The number of compartments
	
	vector <Transition> tra;               // Stores details of transitions
	unsigned int ntra;                     // The number of transitions
	
	string index;                          // The mathematical index used to represent different compartments

	vector <Island> island;                // Compartments are split into "island" which are not connected
	
	Coord coord;                           // The coordinate system used (cartesian on latlng)
};

struct InitCondValue {                   // Sampled initial condition used in the state
	InitCondType type;                     // The type of initial condition
	vector <double> cpop;   
};

struct InitCond {                        // Stores information about initial conditions in data
	InitCondType type;                     // The type of initial condition
	vector <double> cpop;                   // The initial populations (if type POP_INIT)
};

struct Xvector {                         // Stores covariate in fixed effect
	vector <string> ind_list;              // The individuals
	vector <double> value;                 // The covariate
};

struct FixedEffect {                     // Stores an individual fixed IndEffMult
	string name;                           // The name of the parameter
	unsigned int th;                       // The number of the parameter
	Xvector X_vector;                      // The covariate as loaded from the file
	unsigned int line_num;                 // The line number
};

/*
struct MarkovEqnIndRef {                 // References when the event in in a Markov equation
	unsigned int e;                        // The Markov equation
	unsigned int ti;                       // The time division
	unsigned int index;                    // The index in the time division
};
*/

struct MarkovEqnRef {                    // References when the individul is in the Markov equation
	unsigned int e;                        // The Markov equation
	unsigned int index;                    // The index in the time division
};


struct IndDivRef {                       // Reference markov_eqn div from ind
	unsigned int e;                        // Markov equation
	unsigned int ti;
	unsigned int index;
};

struct Event {                           // A transition event
	EventType type;                        // The type of the event
	unsigned int tr_gl;                    // The global transition which occurs (for move used for cl)
	unsigned int c_gl;                     // The compartment (for entry events) (for move used for c)
	unsigned int c_after;                  // The compartment after the transition
	double t;                              // The time at which transition occurs
	IndDivRef markov_eq_ref;               // References the markov eqn the transition occurs in
	bool fixed;                            // Detemines if event is fixed
};

struct ClEvent {                         // References within individual classification
	unsigned int index_ev;                   // References the event timeline      
	unsigned int index_nmtrans;            // References nm_trans (if appropriate)
	unsigned int nm_ref;                   // Reference which nm_trans
};

struct Individual {                      // Stores information about an individual
	string name;                           // The name of the individual
	unsigned int cinit;                    // The initial compartment for individual
	unsigned int c;                        // The current compartment for individual
	vector <Event> ev;                     // Events which move individual between compartments
	vector < vector <ClEvent> > cl_ev;     // Classification event (which reference ev and nm_trans)
	vector <double> ie;                    // Individual effect values
	vector <double> exp_ie;                // Expoential of individual effect values e^(ie)
	vector <double> X;                     // The X value for the fixed effects
	vector <double> exp_fe;                // Expoential of individual fixed effect values e^(X*fe)

	// Used during simulation? 
	vector <MarkovEqnRef> markov_eqn_ref;  // References Markov equations that the individual is on
};


struct TrigEvent {                       // An event trigger (e.g. ind enters, leaves, NM transition)
	EventType type;                        // The type of event
	unsigned int c;                        // The state (if entering)
	unsigned int i;                        // The individual on which the event happens
	unsigned int trg;                      // The global transtion
	double t;                              // The time the event happens
};

struct TrigEventDiv {                    // Stores a division of non-Markovian next events (simulation) 
	vector <TrigEvent> ev;                   // Events
};

struct MarkovNode {
	vector <unsigned int> child;
	unsigned int parent;
};

struct MarkovTree {                      // Stores a tree to be able to sample Markov events
	vector <MarkovNode> node;
};

struct IndTransRef {
	unsigned int i;                        // The individual number 
	unsigned int index;                    // The index in ind.markov_eqn
	unsigned int tr_gl;                    // The global transition number
};

struct DivIndRef {                       // References an individual event on a markov eqn
	unsigned int i;
	unsigned int index;
};

struct MarkovEqnDiv {                    // Division in markov timeline
	double value;                          // The value of the equation	for the division
	double indfac_int;                     // The integral of individual factors across div
	vector <DivIndRef> ind_trans;          // The individuals which undego transitions
	
	//double indfac_sum;                     // Sum of individual factors acting (at the begining)
	//vector <MarkovEqnEvent> ev;            // Stores an event on the Markov equation timeline
};

struct MarkovEqn {                       // Stores information about the Markov equations
	unsigned int eqn_ref;                  // The equation number
	bool time_vari;                        // Determines if equation has time variation
	bool source;                           // Set to true if equation is for a source
	vector <unsigned> source_tr_gl;        // Global transition for source (if applicable)
	bool ind_variation;                    // Set to true if there is individual variation
	vector <unsigned int> ind_eff_mult;    // References an individual effect which multiplies equation
	vector <unsigned int> fix_eff_mult;    // References fixed effects whihc multiply equation
};

struct MarkovEqnVariation {              // Stores variation in markov equations for a species
	// These are used in simulation
	bool time_vari;                     	 // Determines if equation has time variation
	
	double value;                          // The value
	double indfac_sum;                     // The sum of the potential individual effect acting on equation
	vector <IndTransRef> ind_tra;          // The individual transtions associated with equation 
	
	// These are used in inference
	//vector <IndTransRef> ind_tra_init;     // The initial individuals attached to equations
	double dt;                             // The timestep used for divisions (set to UNSET for no time)
	vector <MarkovEqnDiv> div;             // Divisions
};

struct EventData {
	EventData(){ c = UNSET; cl = UNSET; tr = UNSET; t = UNSET;}
	EventType type;                        // The type of the event
	unsigned int c;                        // The compartment (if entering or moving)
	unsigned int cl;                       // The classification
	unsigned int tr;                       // The transition which occurs
	double t;                              // The time at which transition occurs
};

struct ObsData {
	ObsData(){ so = UNSET; ref = UNSET; c = UNSET; cl = UNSET; tr = UNSET; test_res = UNSET; t = UNSET;}
	ObsType type;                          // The type of the observation
	unsigned int so;                       // References the data source
	unsigned int ref;                      // References population or trans pop data
	unsigned int c;                        // The compartment (if entering or moving)
	unsigned int cl;                       // The classification
	unsigned int tr;                       // The transition which occurs
	unsigned int Se_obs_eqn;               // References the Se observation equation
	unsigned int Sp_obs_eqn;               // References the Sp observation equation
	bool test_res;                         // The result of a diagnostic test
	double t;                              // The time at which transition occurs
	//vector <bool> filt;                    // The filter used for the observation
};

struct IndData {                         // Stores data on individuals
	string name;                           // The name of the individual
	vector <EventData> ev;                 // Any events on individual (enter, leave, move)
	vector <ObsData> obs;                  // Any observed events on individual
};

struct WarnData {                        // Stores warning messages
	string te;                             // The text for the warning
	unsigned int line_num;                 // The line number in the input file
};

struct SpeciesSimp {                     // A simplified version of species to pass to equation
		SpeciesSimp(string name_, const vector <Classification> &cla, vector <IndEffect> &ind_effect, vector <FixedEffect> &fix_effect, const vector <unsigned int> &comp_mult, const vector <CompGlobal> &comp_gl) : cla(cla), ind_effect(ind_effect), fix_effect(fix_effect), comp_mult(comp_mult), comp_gl(comp_gl)
		{
			name = name_;
		}
 
		string name;                            // The name of the species
		const vector <Classification> &cla;     // Stores classification within the species
		vector <IndEffect> &ind_effect;         // Any individual effects within species
		vector <FixedEffect> &fix_effect;       // Any fixed effects within species	
		const vector <unsigned int> &comp_mult; // Used to convert from classification definition to global
		const vector <CompGlobal> &comp_gl;      // The global compartments
};


struct BranchParam {                        // Groups together branching parameters
	vector <unsigned int> group;
};

struct ParamTag {                           // Used to check tags on parameters specified correctly
	string val;
	string tag;
};


struct Derive {                             // Stores derived parameters
	string name;                              // The name of the derived quantity
	bool time_dep;                            // Time dependency
	vector <Dependency> dep;                  // Dependency for derived quantity (excluding time)
	vector <EquationInfo> eq;                 // The equations for each quantity
};

struct DeriveOutput {                       // Calculates any derived outputs
	vector <string> value_str;                // Value  
};

struct Change {                             // Stores changes to state
	double num;
	vector <double> vec;
};
	
struct Like {                               // Stores a likelihood (or change in likelihood
	Like(){ 
		obs = 0;
		prior = 0; 
		spline_prior = 0; 
		dist = 0;
		markov = 0; 
		nm_trans = 0;
		ie = 0;
	}
	double obs;                               // The observation probability
	double prior;                             // The prior
	double spline_prior;                      // The spline prior
	double dist;                              // The probability of internal distributions
	double markov;                            // The likelihood of Markovian transitions
	double nm_trans;                          // The likelihood of non-Markovian transitions
	double ie;                                // The likelihood of individual effects
};
	
	
struct IndRef {                             // An individual reference
	unsigned int p;                           // Species
	unsigned int i;                           // Individual
};

struct IndSampler {                         // Useds as an individual sampler
	unsigned int ntr;
	unsigned int nac;
	double si;
};

struct PopData {                            // Stores information about a population measurement
	unsigned int so;                          // References the source 
	double t;                                 // The time of the population measurement
	double value;                             // The value of the measurement
	double sd;                                // The standard deviation in the observation model
	vector <bool> filt;                       // Filters for all the global compartments  
};


struct PopTransData {                       // Stores information about a population transitions
	unsigned int so;                          // References the source 
	double tmin;                              // The start time
	double tmax;                              // The end time
	double value;                             // The value of the measurement
	double sd;                                // The standard deviation in the observation model
	vector <bool> filt;                       // Filters for all the global transitions 
};

struct ParticleSpecies {
	vector <double> cpop_init;             // The initial population in the compartments 
	vector < vector <double> > trans_num; 
	vector <Individual> individual;
};

struct Particle {                                     // Stores information from state
	vector <double> param_val;                          // The parameter value
	vector <ParticleSpecies> species;                   // Species state data
	vector <DeriveOutput> dir_out;                      // Derived outputs
	Like like;                                          // Stores the overall likelihood and prior
	double w;                                           // Particle weight (used in ABC-SMC)
};

#endif

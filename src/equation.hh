#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "const.hh"

enum EqItemType { LEFTBRACKET, RIGHTBRACKET, FUNCDIVIDE, ADD, TAKE, MULTIPLY, DIVIDE, REG, PARAMETER, SPLINE, POPNUM, TIME, IE, ONE, FE, NUMERIC, EXPFUNC, SINFUNC, COSFUNC, LOGFUNC, POWERFUNC, THRESHFUNC, UBOUNDFUNC, STEPFUNC, MAXFUNC, MINFUNC, ABSFUNC, SQRTFUNC, NOOP};

struct EqItem {                            // An individual operation in a calculation
	EqItem(){ num = UNSET; index = UNSET; constant = UNSET;}
	
	EqItemType type;
	unsigned int num;
	unsigned int index;
	double constant; 
};

struct Calculation                         // Stores a calculation (made up of operations)
{
	vector <EqItem> item;                    // Items used to make the calculation
	EqItemType op;                           // The operator used in the calculation
	unsigned int reg_store;                  // The register where to store the results
};

struct PopCalculation                      // Used to describe a population calculation
{
	vector <Calculation> calc;               // A calculation
	unsigned int po;                         // The population which multiplies calculation 
};

struct LinearCalculation                   // Stores a linear calculation
{
	PopCalculation no_pop_calc;              // Calculates non-population part
	vector <PopCalculation> pop_calc;        // Calculates parts multiplying populations
};

struct Linearise {                         // Information about linearising an equation 
	bool on;                                 // Determines if linearisation is possible
	vector <Calculation> no_pop_calc;        // The calculation not involving populations
	bool no_pop_calc_time_dep;               // Determines if time dependent
	
	vector < vector <Calculation> > pop_grad_calc; // Calculations involving populations (ordered by pop_ref)
	bool pop_grad_calc_time_dep;             // Determines if time dependent
	
	bool multi_source;                       // Detemines if comes from multiple sources (used in transmission trees)
};

class Equation                             // Stores information about the model
{
	public:
		vector <double> timer;
		
		bool plfl;                             // Determines if calculation is printed (for diagnostic purposes)
	
		vector <Calculation> calc;             // Stores calculation

		EqItem ans;                            // Stores the answer 
		unsigned int nreg;                     // The number of registers used
	
		vector <ParamRef> param_ref;           // Stores the parameters used in the equation
		
		vector <unsigned int> pop_ref;         // Stores the populations used in the equation
		
		vector <unsigned> source_tr_gl;        // Stores any global transitions from source (if eqn SOURCE_RATE)
		
		vector <unsigned int> ind_eff_mult;    // Potential individual effect multipier
		
		vector <bool> ind_eff_exist;           // Determines if individuals effect exists
		
		vector <unsigned int> fix_eff_mult;    // Potential fixed effect multipier
		
		string te_init;                        // The raw text used to generate the equation
		string te;                             // Text for the equation (modified during operation)
		string te_raw;                         // The text used to print to terminal 
		
		bool time_vari;                        // Determines if equation has time variation
		
		bool infection_trans;                  // Determines if relates to an infection transition
		
		vector < vector <unsigned int> > comp_pref_convert; // Converts from c to pref (used for infection transition) [p][cgl]		
		
		bool ind_eff;                          // Set if an individual effect exists
		
		EqnType type;                          // The type of the equation
		
		unsigned int sp_p, sp_cl;              // Species number and classification (if appropriate)
		unsigned int sp_c;                     // The compartment number (used for sum max function)

		unsigned int ti_fix;                   // The fixed time point (e.g. for pop measurements)
		
		string warn;                           // Stores error warning message

		unsigned int markov_eqn_ref;           // References Markov equation
		
		unsigned int line_num;                 // The line number in the input file
		
		vector <double> time;                  // Stores the times of the time points (such that t can be used)
		
		Linearise linearise;                   // Used for accelerated likelihood calculation
		
		Equation(string tex, EqnType ty, unsigned int p, unsigned int cl, unsigned int c, bool inf_trans, unsigned int tif, unsigned int li_num, vector <SpeciesSimp> &species, vector <Param> &param, vector <Spline> &spline, vector <ParamVecEle> &param_vec, vector <Population> &pop, Hash &hash_pop, const vector <double> &timepoint);
		double calculate_param_only(const vector <double> &param_val) const;
		double calculate_param_only_ti_fix(const vector <double> &param_val, const vector <SplineValue> &spline_val) const;
		double calculate_no_popnum(unsigned int ti, const vector <double> &param_val, const vector <SplineValue> &spline_val) const; 
		double calculate(unsigned int ti, const vector <double> &popnum, const vector <double> &param_val, const vector <SplineValue> &spline_val) const;
		void print_calculation() const;
		void print_ca(const Calculation &ca) const;
		
		double calculate_indfac(const Individual &ind, unsigned int ti, const vector <double> &popnum, const vector <double> &param_val, const vector <SplineValue> &spline_val) const;
		double indfac(const Individual &ind) const;
		bool param_linear(unsigned int th) const;
		bool is_zero() const;
		bool is_one() const;
		void simplify();
		vector < vector <unsigned int> > get_reg_used(const vector <Calculation> &calc) const;
		void check_reg_used(const vector <Calculation> &calc, vector < vector <unsigned int> > &reg_used) const;
		void setup_references();
		
	private:
		vector <unsigned int> getallcomp(string st);
		void print_operations(vector <EqItem> op) const;
		
		vector <string> find_list_from_index(string ind, double dist_max) const;
		vector <unsigned int> get_all_comp(unsigned int p, string st);
		double get_float(unsigned int i, unsigned int &raend) const;
		ParamRef get_param_name(unsigned int i, double &dist, unsigned int &raend);
		unsigned int get_pop(unsigned int i, unsigned int &raend);
		unsigned int get_ie(unsigned int i, unsigned int &raend);
		unsigned int get_fe(unsigned int i, unsigned int &raend);
		bool pop_combine(const vector <EqItem> &op, unsigned int i,	unsigned int &popcomb, double &popcombnum) const;
		bool quant(const vector <EqItem> &op, int i) const;
		bool is_func(const vector <EqItem> &op, int i) const;
		bool optype(const vector <EqItem> &op, int i, EqItemType type) const;
		int prio(const vector <EqItem> &op, int i) const;
		void unravel_sum();
		void minus_sign_adjust();
		double calculate_operation(EqItemType op, vector <double> &num) const;
		
		vector <EqItem> extract_operations();
		void combine_populations(vector <EqItem> &op) const;
		void create_calculation(vector <EqItem> &op);
		void extract_ind_eff();
		void extract_fix_eff();
		CompRef find_comp_from_name(unsigned int p, string te) const;
		void replace_reg(unsigned int num, EqItem it_change);
		void setup_comp_pref_convert();
		void check();
		void remove_unused();
		double add_const(EqItem item1, EqItem item2) const;
		double mult_const(EqItem item1, EqItem item2) const;
		void set_time_vari();
		double find_dist(unsigned int c, unsigned int cc, const vector <Compartment> &comp, Coord coord) const;
		double get_distance(const ParamProp &pp);
				
		vector <SpeciesSimp> &species;       // References the species from the model
		unsigned int nspecies;
		vector <Param> &param;               // References the parameters from the model
		vector <Spline> &spline;             // References splines from the model
		vector <ParamVecEle> &param_vec;     // References the param_vec from the model
		vector <Population> &pop;            // References the populations from the model
		Hash &hash_pop;                      // Hash table for pop
		const vector <double> &timepoint;    // References time points
		
	// In equation_linearise.cc
	public:
		void calculate_linearise();
		void setup_param_ref();
		void calculate_pop_ref();
		LinearCalculation convert_to_linear_calculation(const EqItem &it, EqItemType op, const vector <LinearCalculation> &lin_calc) const;
		void calc_mult(vector <Calculation> &calc, const vector <Calculation> &calc2) const;
		void calc_add(vector <Calculation> &calc, const vector <Calculation> &calc2) const;
		bool ca_is_one(const Calculation &ca) const;
		void calc_div(vector <Calculation> &calc, const vector <Calculation> &calc2) const;
		void print_linear_calc(string te, const LinearCalculation &lin) const;
		void print_calc(string st, const vector <Calculation> &calc) const;
		void single_param_func(Calculation ca, LinearCalculation &lin, const vector <LinearCalculation> &lin_calc) const;
		void two_param_func(Calculation ca, LinearCalculation &lin, const vector <LinearCalculation> &lin_calc) const;
		double calculate_calculation(const vector <Calculation> &calc, unsigned int ti, const vector <double> &param_val, const vector <SplineValue> &spline_val) const;
		double calculate_calculation_spline_store(const vector <Calculation> &calc, unsigned int ti, const vector <double> &param_val, const vector <SplineValue> &spline_val) const;
		bool calc_time_dep(const vector <Calculation> &calc) const;
		double calculate_calculation_notime(const vector <Calculation> &calc, const vector <double> &param_val) const;
		double calculate_linearise_check(unsigned int ti, const vector <double> &popnum, const vector <double> &param_val, const vector <SplineValue> &spline_val) const;
		bool equal_calc(const vector <Calculation> &calc1, const vector <Calculation> &calc2) const;
		vector <double> calculate_popnum_gradient(const vector <double> &param_val) const;
		double calculate_no_pop(unsigned int ti, const vector <double> &param_val, const vector <SplineValue> &spline_val) const;
		InfSourceSampler setup_source_sampler(unsigned int ti, const vector <double> &popnum, const vector <double> &param_val, const vector <SplineValue> &spline_val) const;
};


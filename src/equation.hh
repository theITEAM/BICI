#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "const.hh"

enum EqItemType { LEFTBRACKET, RIGHTBRACKET, FUNCDIVIDE, ADD, TAKE, MULTIPLY, DIVIDE, REG, PARAMETER, PARAMVEC, SPLINE, SPLINEREF, POPNUM, TIME, IE, ONE, FE, NUMERIC, EXPFUNC, SINFUNC, COSFUNC, LOGFUNC, POWERFUNC, THRESHFUNC, UBOUNDFUNC, STEPFUNC, MAXFUNC, MINFUNC, ABSFUNC, SQRTFUNC, SIGFUNC, TINT, INTEGRAL, DERIVE, NOOP};

struct EqItemList {                         // Allows for a list of equation items
	EqItemType type;
	unsigned int num;
	unsigned int prev;
	unsigned int next;
};

struct EqItem {                            // An individual operation in a calculation
	EqItem(){ num = UNSET;}
	
	EqItemType type;
	unsigned int num;
};

struct PopRefFromPo {                      // Used to get pop_ref from population number
	unsigned int i;                          // Index on pop_ref
	unsigned int po;                         // Population number
};

struct Calculation {                       // Stores a calculation (made up of operations)
	vector <EqItem> item;                    // Items used to make the calculation
	EqItemType op;                           // The operator used in the calculation
};

struct Integral {                          // Stores the calculation within an integral
	vector <Calculation> calc;    
	unsigned int ti_min;
	unsigned int ti_max;
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

	vector < vector <PopRefFromPo> > pop_ref_from_po; // Gets pop_ref from population	
		
	bool multi_source;                       // Detemines if comes from multiple sources (used in transmission trees)
	
	void init_pop_ref_from_po(const vector <unsigned int> &pop_ref);
	unsigned int get_pop_ref(unsigned int po) const;	
};

class Equation                             // Stores information about the model
{
	public:
		vector <double> timer;
		
		bool plfl;                             // Determines if calculation is printed (for diagnostic purposes)
	
		vector <Calculation> calcu;            // Stores calculation

		vector <Integral> integral;            // Calculation for integral

		vector <ParamRef> param_ref;           // Stores the parameters used in the equation
		
		vector <DeriveRef> derive_ref;         // Stores derived used in the equation
		
		vector <double> cons;                  // Stores any constant values
		
		vector <unsigned int> pop_ref;         // Stores the populations used in the equation
		
		vector <unsigned> source_tr_gl;        // Stores any global transitions from source (if eqn SOURCE_RATE)
		
		vector <unsigned int> ind_eff_mult;    // Potential individual effect multipier
		
		vector <bool> ind_eff_exist;           // Determines if individuals effect exists
		
		vector <unsigned int> fix_eff_mult;    // Potential fixed effect multipier
		
		string te;                             // Text for the equation (modified during operation)
		string te_raw;                         // The text used to print to terminal 
		
		bool time_vari;                        // Determines if equation has time variation
		
		bool infection_trans;                  // Determines if relates to an infection transition
		
		vector < vector <unsigned int> > comp_pref_convert; // Converts from c to pref (used for infection transition) [p][cgl]		
		
		bool ind_eff;                          // Set if an individual effect exists
		
		EqnType type;                          // The type of the equation
		
		unsigned int sp_p, sp_cl;              // Species number and classification (if appropriate)
		
		//unsigned int sp_c;                     // The compartment number (used for sum max function)

		unsigned int ti_fix;                   // The fixed time point (e.g. for pop measurements)
		
		string warn;                           // Stores error warning message

		unsigned int markov_eqn_ref;           // References Markov equation
		
		unsigned int line_num;                 // The line number in the input file
		
		vector <double> time;                  // Stores the times of the time points (such that t can be used)
		
		Linearise linearise;                   // Used for accelerated likelihood calculation
		
		Equation(string tex, EqnType ty, unsigned int p, unsigned int cl, bool inf_trans, unsigned int tif, unsigned int li_num, vector <SpeciesSimp> &species, vector <Param> &param, const vector <Derive> &derive, vector <Spline> &spline, vector <ParamVecEle> &param_vec, vector <Population> &pop, Hash &hash_pop, const vector <double> &timepoint);
		double calculate_param_only(const vector <double> &param_val) const;
		double calculate_param_only_ti_fix(const vector <double> &param_val, const vector <SplineValue> &spline_val) const;
		double calculate_no_popnum(unsigned int ti, const vector <double> &param_val, const vector <SplineValue> &spline_val) const; 
		double calculate(unsigned int ti, const vector <double> &popnum, const vector <double> &param_val, const vector <SplineValue> &spline_val) const;
		double calculate_integral(unsigned int i, const vector < vector <double> > &popnum, const vector <double> &param_val, const vector <SplineValue> &spline_val, const vector < vector < vector <double> > > &derive_val) const;
		double calculate_derive(unsigned int ti, const vector < vector <double> > &popnum, const vector <double> &param_val, const vector <SplineValue> &spline_val, const vector < vector < vector <double> > > &derive_val) const;
		void print_calculation() const;
		void print_ca(unsigned int i, const Calculation &ca) const;
		
		double calculate_indfac(const Individual &ind, unsigned int ti, const vector <double> &popnum, const vector <double> &param_val, const vector <SplineValue> &spline_val) const;
		double indfac(const Individual &ind) const;
		bool param_linear(unsigned int th) const;
		double is_num() const;
		bool is_zero() const;
		bool is_one() const;
		void insert_reg(vector <Calculation> &calc, vector <bool> &calc_on);
		void simplify(vector <Calculation> &calc);
		vector < vector <unsigned int> > get_reg_used(const vector <Calculation> &calc) const;
		void check_reg_used(const vector <Calculation> &calc, vector < vector <unsigned int> > &reg_used) const;
		void setup_references();
		double get_calc_hash_num(const vector <Calculation> &calc) const;
		unsigned int add_cons(double val);
		
	private:
		vector <unsigned int> getallcomp(string st);
		void print_operations(const vector <EqItem> &op) const;
		
		vector <string> find_list_from_index(string ind, double dist_max, string comp_max) const;
		vector <unsigned int> get_all_comp(unsigned int p, string te);
		double get_float(unsigned int i, unsigned int &raend) const;
		ParamRef get_param_name(unsigned int i, double &dist, unsigned int &raend);
		DeriveRef get_derive_name(unsigned int i, unsigned int &raend);
		unsigned int get_pop(unsigned int i, unsigned int &raend);
		unsigned int get_ie(unsigned int i, unsigned int &raend);
		unsigned int get_fe(unsigned int i, unsigned int &raend);
		bool pop_combine(const vector <EqItem> &op, unsigned int i,	unsigned int &popcomb, double &popcombnum) const;
		bool quant(const vector <EqItem> &op, int i) const;
		bool quantl(const vector <EqItemList> &opl, unsigned int i) const;
		bool is_func(const vector <EqItem> &op, int i) const;
		bool is_funcl(const vector <EqItemList> &opl, unsigned int i) const;
		bool optype(const vector <EqItem> &op, int i, EqItemType type) const;
		bool optypel(const vector <EqItemList> &opl, unsigned int i, EqItemType type) const;
		int prio(const vector <EqItem> &op, int i) const;
		int priol(const vector <EqItemList> &opl, unsigned int i) const;
		void unravel_sum();
		void minus_sign_adjust();
		double calculate_operation(EqItemType op, vector <double> &num) const;
		
		vector <EqItem> extract_operations();
		unsigned int get_integral_bound(string st);
		unsigned int extract_integral(const string &te, unsigned int i, vector <EqItem> &op);
		void combine_populations(vector <EqItem> &op) const;
		vector <EqItemList> create_opl(const vector <EqItem> &op) const;
		vector <EqItem> create_op_from_opl(const vector <EqItemList> &opl) const;
		void erase_opl(unsigned int i, unsigned int n, vector <EqItemList> &opl) const;
		EqItem get_op(const vector <EqItemList> &opl, unsigned int i) const;
		unsigned int get_opl_start(const vector <EqItemList> &opl) const;
		unsigned int get_opl_num(const vector <EqItemList> &opl) const;
		vector <Calculation> create_calculation(vector <EqItem> &op);
		void extract_ind_eff();
		void extract_fix_eff();
		CompRef find_comp_from_name(unsigned int p, string te) const;
		void replace_reg(const vector <EqItem> &reg_replace, vector <Calculation> &calc, const vector <bool> &calc_on);
		void setup_comp_pref_convert();
		void check();
		void check_opl(const vector <EqItemList> &opl) const;
		void remove_unused(vector <Calculation> &calc, vector <bool> &calc_on);
		unsigned int add_const(EqItem item1, EqItem item2);
		unsigned int mult_const(EqItem item1, EqItem item2);
		void set_time_vari();
		double find_dist(unsigned int c, unsigned int cc, const vector <Compartment> &comp, Coord coord) const;
		double get_distance(const ParamProp &pp);
		string op_name(EqItemType type) const;
		void check_repeated_operator(const vector <EqItem> &op);
		void time_integral(vector <EqItem> &op);
		void replace_minus(vector <EqItem> &op);
		unsigned int add_param_ref(const ParamRef &pref);
				
		vector <SpeciesSimp> &species;       // References the species from the model
		unsigned int nspecies;
		vector <Param> &param;               // References the parameters from the model
		const vector <Derive> &derive;       // Reference derived quantities in the model
		vector <Spline> &spline;             // References splines from the model
		vector <ParamVecEle> &param_vec;     // References the param_vec from the model
		vector <Population> &pop;            // References the populations from the model
		Hash &hash_pop;                      // Hash table for pop
		const vector <double> &timepoint;    // References time points
		
	// In equation_linearise.cc
	public:
		void calculate_linearise();
		void print_linear_final() const;
		void remove_unused_param_ref();
		void calculate_pop_ref();
		
		LinearCalculation convert_to_linear_calculation(const EqItem &it, EqItemType op, const vector <LinearCalculation> &lin_calc);
		void calc_mult(vector <Calculation> &calc, const vector <Calculation> &calc2);
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
		bool equal_calc(const vector <Calculation> &calc1, const vector <double> &cons, const vector <Calculation> &calc2, const vector <double> &cons2) const;
		vector <double> calculate_popnum_gradient(const vector <double> &param_val) const;
		double calculate_no_pop(unsigned int ti, const vector <double> &param_val, const vector <SplineValue> &spline_val) const;
		InfSourceSampler setup_source_sampler(unsigned int ti, const vector <double> &popnum, const vector <double> &param_val, const vector <SplineValue> &spline_val) const;
};


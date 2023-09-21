#ifndef BICI__EQUATION_HH
#define BICI__EQUATION_HH

#include <vector>

using namespace std;

//#include "species.hh"
#include "struct.hh"
#include "const.hh"


enum EqItemType { LEFTBRACKET, RIGHTBRACKET, FUNCDIVIDE, ADD, TAKE, MULTIPLY, DIVIDE, REG,
     PARAMETER, SPLINE, POPNUM, TIME, IE, ONE, FE, NUMERIC, EXPFUNC, SINFUNC, COSFUNC, LOGFUNC, POWERFUNC, THRESHFUNC, STEPFUNC, MAXFUNC, MINFUNC, ABSFUNC, SQRTFUNC};

struct EqItem{
	EqItem(){ num = UNSET; index = UNSET; constant = UNSET;}
	
	EqItemType type;
	unsigned int num;
	unsigned int index;
	double constant; 
};


struct Calculation
{
	vector <EqItem> item;           // Items used to make the calculation
	EqItemType op;                  // The operator used in the calculation
	unsigned int reg_store;         // The register where to store the results
};

class Equation                                            // Stores information about the model
{
	public:
		bool plfl;                      // Determines if calculation is printed (for diagnostic purposes)
	
		vector <Calculation> calc;

		EqItem ans;                     // Stores the answer 
		unsigned int nreg;              // The number of registers used
	
		vector <ParamRef> param_ref;    // Stores the parameters used in the equation
		
		vector <unsigned int> pop_ref;  // Stores the populations used in the equation
		
		vector <unsigned> source_tr_gl; // Stores any global transitions from source (if eqn SOURCE_RATE)
		
		vector <unsigned int> ind_eff_mult;// Potential individual effect multipier
		
		vector <unsigned int> fix_eff_mult;// Potential fixed effect multipier
		
		string te_init;                 // The raw text used to generate the equation
		string te;                      // The text for the equation (modified during operation)
		string te_raw;                  // The text used to print to terminal 
		
		EqnType type;                   // The type of the equation
		
		unsigned int sp_p, sp_cl;       // The species number and classification (if appropriate)
		
		string warn;                    // Stores error warning message

		unsigned int markov_eqn_ref;    // References the set of equations used for Markov transitions
		
		unsigned int line_num;
		
		vector <double> time;           // Stores the times of the time points (such that t can be used)
		
		Equation(string tex, EqnType ty, unsigned int p, unsigned int cl, unsigned int li_num, vector <SpeciesSimp> &species, vector <Param> &param, vector <Spline> &spline, vector <ParamVecEle> &param_vec, vector <Population> &pop, const vector <double> &timepoint);
		double calculate_param_only(const vector <double> &param_val) const;
		double calculate(unsigned int ti, const vector <double> &popnum, const vector <double> &param_val, const vector <SplineValue> &spline_val) const;
		void print_calculation() const;
		
		double calculate_indfac(const Individual &ind, unsigned int ti, const vector <double> &popnum, const vector <double> &param_val, const vector <SplineValue> &spline_val) const;
		double indfac(const Individual &ind) const;
	
	private:
		vector <unsigned int> getallcomp(string st);
		void print_operations(vector <EqItem> op) const;
		
		vector <string> find_list_from_index(string ind) const;
		vector <unsigned int> get_all_comp(unsigned int p, string st);
		double get_float(unsigned int i, unsigned int &raend) const;
		ParamRef get_param_name(unsigned int i, unsigned int &raend);
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
		void check();

		vector <SpeciesSimp> &species;       // References the species from the model
		unsigned int nspecies;
		vector <Param> &param;               // References the parameters from the model
		vector <Spline> &spline;             // References splines from the model
		vector <ParamVecEle> &param_vec;     // References the param_vec from the model
		vector <Population> &pop;            // References the populations from the model
		const vector <double> &timepoint;    // References time points
};
#endif

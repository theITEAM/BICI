#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "precalc.hh"
#include "hash_simp.hh"
#include "const.hh"

class Equation                             // Stores information about an equation
{
	public:
		vector <double> timer;
		
		bool plfl;                             // Determines if calculation is printed (for diagnostic purposes)
	
		vector <Calculation> calcu;            // Stores calculation

		vector <Integral> integral;            // Calculation for integral

		vector <SumInfo> sum_info;             // Information about sums

		vector <string> name_store;            // Used to store names of parameters and populations
		
		vector <ParamIndex> param_index;       // Information about parameter in index form
		
		vector <PopIndex> pop_index;           // Information about popualtion in index form
		
		vector <ParamRef> param_ref;           // Stores the parameters used in the equation
		
		vector <DeriveRef> derive_ref;         // Stores derived used in the equation
		
		vector <PopTimeRef> pop_time_ref;      // References a population time
		
		bool contain_population;               // Determines if contains a population
		
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
		
		unsigned int ti_fix;                   // The fixed time point (e.g. for pop measurements)
		
		string warn;                           // Stores error warning message

		unsigned int markov_eqn_ref;           // References Markov equation
		
		unsigned int line_num;                 // The line number in the input file
		
		vector <double> time;                  // Stores the times of the time points (such that t can be used)
		
		bool precalc_done;                     // Stores whether precalculation has been done
		
		Linearise lin;                         // Used for accelerated likelihood calculation
		
		Equation(EquationInfo &eqi, unsigned int tif, const vector <SpeciesSimp> &species, vector <Param> &param, vector <Prior> &prior, const vector <Derive> &derive, const vector <Spline> &spline, const vector <ParamVecEle> &param_vec, vector <Density> &density, vector <Population> &pop, Hash &hash_pop, vector < vector <PopComb> > &popcomb, vector <PopCombWeight> &popcombw, Constant &constant, const vector <double> &timepoint, const Details &details, vector <Define> &define);
		
		
		DefineOpStore get_define_op_store();
		void print_info() const;
		void print_calculation() const;
		void print_item(const EqItem &it) const;
		void print_ca(unsigned int i, const Calculation &ca) const;
		void print_mult(string name, const vector <EqItem> &item) const;
		
		double indfac(const Individual &ind) const;
		bool param_linear(unsigned int th) const;
		double is_num() const;
		bool is_zero() const;
		bool is_one() const;
		
		vector < vector <unsigned int> > get_reg_used(const vector <Calculation> &calc) const;
		void check_reg_used(const vector <Calculation> &calc, vector < vector <unsigned int> > &reg_used) const;
		void setup_references();
		double get_calc_hash_num(const vector <Calculation> &calc) const;
		
	private:
		vector <unsigned int> getallcomp(string st);
		void print_operations(const vector <EqItem> &op) const;
		
		void substitute_define(vector <EqItem> &op, bool &pass_again);
		void set_dep(vector <unsigned int> &dep, vector <IndexNotSet> &index_not_set, const vector <Substitution> &sub, const vector <string> &sum_index) const;
		void set_pop_dep(PopIndex &pi, const vector <Substitution> &sub, const vector <string> &sum_index) const;
		void copy_op(vector <EqItem> &op_new, const	vector <EqItem> &op_st, const vector <SumInfo> &sum_info_st, const vector <ParamIndex> &param_index_st, const vector <PopIndex> &pop_index_st, const vector <Substitution> &sub);
		CompPos find_list_from_index(string ind, double dist_max, string comp_max) const;
		vector <unsigned int> get_all_comp(const PopIndex &pind);
		double get_float(unsigned int i, unsigned int &raend) const;
		ParamIndex get_param_name(unsigned int i, unsigned int &raend);
		DeriveRef get_derive_name(unsigned int i, unsigned int &raend);
		PopIndex get_pop(unsigned int i, unsigned int &raend);
		bool pind_cl_exist(unsigned int cl, const PopIndex &pind) const;
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
		void unravel_sum(vector <EqItem> &op);
		void minus_sign_adjust();
		vector <EqItem> extract_operations();

		unsigned int get_other_bracket(unsigned int i, const vector <EqItem> &op) const;
		unsigned int get_integral_bound(string st);
		void simplify_operations(vector <EqItem> &op);
		void convert_param_index(vector <EqItem> &op);
		void convert_pop_index(vector <EqItem> &op);
		vector <unsigned int>	get_pop_hash_vec(const PopIndex &pind) const;
		string get_pop_name(const PopIndex &pind) const;
		unsigned int extract_sum(const string &te, unsigned int i, vector <EqItem> &op);
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
	
		void bracket_check(string st, const vector <EqItem> &op) const;
		void setup_comp_pref_convert();
		void check();
		void check_opl(const vector <EqItemList> &opl) const;
		void remove_unused(vector <Calculation> &calc, vector <bool> &calc_on);
		
		
		void set_time_vari();
		double find_dist(unsigned int c, unsigned int cc, const vector <Compartment> &comp, Coord coord) const;
		double geo_dist(double lat1, double lng1, double lat2, double lng2) const;
		double get_distance(const ParamIndex &pind);
		double get_identity(const ParamIndex &pind);
		vector <double> set_density(unsigned int p, unsigned int cl, double r, bool rel_den) const;
		double get_density(const ParamIndex &pind, bool rel_den);
		string op_name(EqItemType type) const;
		void check_repeated_operator(const vector <EqItem> &op);
		void time_integral(vector <EqItem> &op);
		void replace_minus(vector <EqItem> &op);
		unsigned int add_param_ref(const ParamRef &pref);
		bool model_type() const;
		bool not_time_vary(const EqItem &it) const;
		
		const vector <SpeciesSimp> &species;       // References the species from the model
		unsigned int nspecies;
		vector <Param> &param;               // References the parameters from the model
		vector <Prior> &prior;               // References the priors in the model
		const vector <Derive> &derive;       // Reference derived quantities in the model
		vector <Define> &define;             // References defined quantities
		const vector <Spline> &spline;             // References splines from the model
		const vector <ParamVecEle> &param_vec;     // References the param_vec from the model
		vector <Density> &density;           // Used for the DEN and RDEN functions 
		vector <Population> &pop;            // References the populations from the model
		Hash &hash_pop;                      // Hash table for pop
		vector < vector <PopComb> > &popcomb; // References the populations from the model
		vector <PopCombWeight> &popcombw;    // References the popcomb weights
		Constant &constant;                  // Stores constants in the model
		const vector <double> &timepoint;    // References time points
		const Details &details;              // Details
		
	// In equation_linearise.cc
	public:
		void calculate_linearise(Precalc &precalc_eqn, const vector <unsigned int> &param_vec_ref, const vector <unsigned int> &spline_ref, SpecPrecalc &spec_precalc);
		bool zero_eqn(const vector <Calculation> &calc) const;
		void print_linear_final(const vector <Calculation> &no_pop_calc, const vector < vector <Calculation> > &popcomb_grad_calc, const vector <unsigned int> &popcomb_list) const;
		void remove_unused_param_ref();
		void calculate_pop_ref();
		
		LinearCalculation convert_to_linear_calculation(const EqItem &it, EqItemType op, const vector <LinearCalculation> &lin_calc);
		void calc_mult(vector <Calculation> &calc, const vector <Calculation> &calc2);
		void calc_add(vector <Calculation> &calc, const vector <Calculation> &calc2) const;
		bool ca_is_one(const Calculation &ca) const;
		void calc_div(vector <Calculation> &calc, const vector <Calculation> &calc2) const;
		void print_linear_calc(string te, const LinearCalculation &lin) const;
		void print_calc(string st, const vector <Calculation> &calc) const;
		bool single_param_func(Calculation ca, LinearCalculation &lin, const vector <LinearCalculation> &lin_calc) const;
		bool two_param_func(Calculation ca, LinearCalculation &lin, const vector <LinearCalculation> &lin_calc) const;
		
		//double calculate_calculation_spline_store(const vector <Calculation> &calc, unsigned int ti, const vector <double> &precalc) const;
		bool calc_time_dep(const vector <Calculation> &calc) const;
		bool item_time_dep(const vector <EqItem> item) const;
		bool it_time_dep(const EqItem &it) const;
		
		double calculate_linearise_check(unsigned int ti, const vector <double> &popcomb, const vector <double> &precalc) const;
		bool equal_calc(const vector <Calculation> &calc1, const vector <Calculation> &calc2) const;
		
		InfSourceSampler setup_source_sampler(unsigned int ti, const vector <double> &popnum, const PV &param_val) const;
		//void get_pop_grad_calc_factorise();
		//void set_precalc();
		EqItem get_precalc(const vector <Calculation> &calc) const;
		
	private:
	
	// In equation_calculate.cc
	public:
		double calculate_constant() const;
		double calculate_param_ti_fix(const vector <double> &precalc) const;
		double calculate(unsigned int ti, const vector <double> &popcomb, const vector <double> &precalc) const;
		double calculate_all_time(unsigned int ti, const vector < vector <double> > &popcomb_t, const vector <double> &precalc) const;
		double calculate_calc(const vector <Calculation> &calc, unsigned int ti, const vector <double> &popcomb, const vector <double> &precalc, const vector < vector < vector <double> > > &derive_val) const;
		double print_working(const vector <Calculation> &calc, unsigned int ti, const vector <double> &popcomb, const vector <double> &precalc, const vector < vector < vector <double> > > &derive_val) const;
		void test_calculate_para(const vector <Calculation> &calc, const vector <unsigned int> &list, const vector < vector <double> > &popcomb_t, const vector <double> &precalc, const vector < vector < vector <double> > > &derive_val, string ref) const;
		vector <double> calculate_para(const vector <Calculation> &calc, const vector <unsigned int> &list, const vector < vector <double> > &popcomb_t, const vector <double> &precalc, const vector < vector < vector <double> > > &derive_val) const;
		double calculate_param(const vector <double> &precalc) const;
		double calculate_no_popcomb(unsigned int ti, const vector <double> &precalc) const; 
		double calculate_integral(unsigned int i, const vector < vector <double> > &popcomb_t, const vector <double> &precalc, const vector < vector < vector <double> > > &derive_val) const;
		double calculate_derive(unsigned int ti, const vector < vector <double> > &popcomb_t, const vector <double> &precalc, const vector < vector < vector <double> > > &derive_val) const;
		double calculate_indfac(const Individual &ind, unsigned int ti, const vector <double> &popcomb, const vector <double> &precalc) const;
		double calculate_pop_grad(unsigned int pref, unsigned int ti, const vector <double> &precalc) const;
		double calculate_pop_grad_no_time(unsigned int pref, const vector <double> &precalc) const;
		double calculate_pop_grad_without_factor_no_time(unsigned int pref, const vector <double> &precalc) const;
		double calculate_pop_grad_without_factor(unsigned int pref, unsigned int ti, const vector <double> &precalc) const;
		double calculate_factor(unsigned int ti, const vector <double> &precalc) const;
		vector <double> calculate_popnum_gradient_without_factor(const vector <double> &precalc) const;
		double calculate_no_pop(unsigned int ti, const vector <double> &precalc) const;
		double calculate_item(const EqItem &it, unsigned int ti, const vector <double> &precalc) const;
		double calculate_item_no_time(const EqItem &it, const vector <double> &precalc) const;
		double calculate_item_old(const EqItem &it, unsigned int ti, const vector <double> &precalc, const vector <double> &precalc_old) const;
		double calculate_item_old_no_time(const EqItem &it, const vector <double> &precalc, const vector <double> &precalc_old) const;
	private:
	
	// In equation_simplify.cc
	public:
		void simplify(vector <Calculation> &calc);
		void extract_popcomb(vector <Calculation> &calcu, Hash &hashw, Hash &hashpc, bool derive);
		
	private:
		void add_it(EqItemType type, vector <EqItem> &op) const;
		EqnRange get_element_from_end(unsigned int end, const vector <EqItem> &op) const;
		EqnRange get_element_from_start(unsigned int start, const vector <EqItem> &op) const;
		void print_eqn_range(const EqnRange &er, const vector <EqItem> &op) const;
		bool contain_time(const EqnRange &er, const vector <EqItem> &op) const;
		bool contain_ie_fe(const EqnRange &er, const vector <EqItem> &op) const;
		bool contain_dep(const EqnRange &er, const vector <string> &dep, const vector <EqItem> &op) const;
		bool contain_linear_pop(unsigned int start, unsigned int end, const vector <EqItem> &op) const;
		vector <EqnRange> eqn_split_add(unsigned int start, unsigned int end, const vector <EqItem> &op) const;	
		vector <EqnRange> eqn_split_mult(unsigned int start, unsigned int end, const vector <EqItem> &op) const;
		void replace_reg(const vector <EqItem> &reg_replace, vector <Calculation> &calc, const vector <bool> &calc_on, bool pl);
		void insert_reg(vector <Calculation> &calc, vector <bool> &calc_on);
		double numeric_one_function(double con, EqItemType ty);
		double numeric_two_function(double con1, double con2, EqItemType ty);
		unsigned int add_const(EqItem item1, EqItem item2);
		unsigned int mult_const(EqItem item1, EqItem item2);
		void add_popcomb(unsigned int i, const vector <PopCombTemp> &pop_co_temp, vector <bool> &remove, const vector < vector <unsigned int> > &used, vector <Calculation> &calcu, Hash &hashw, Hash &hashpc, bool derive);
		void add_popcomb_single(EqItem &item, Hash &hashw, Hash &hashpc, bool derive);
		string popcomb_name(unsigned int i) const;
		
};

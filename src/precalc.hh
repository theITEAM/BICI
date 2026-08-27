#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "const.hh"

class Precalc                             // Stores information about the model
{
	public:
		unsigned long long pcsize;            // The total size for the precalculation
	
		vector <PreCalc> pcalcu;              // Stores precalculation
		
		unsigned int num;
		
		vector <double> timer;                // Used for profiling
		
		Precalc(const vector <SpeciesSimp> &species, const vector <Spline> &spline, const vector <Param> &param, const vector <ParamVecEle> &param_vec, const vector <Population> &pop, Constant &constant, const vector <double> &timepoint, const Details &details);
		
		void clear_timer();
		void print_spec_precalc(string namne, const SpecPrecalc &spec) const;
		string print_spec_precalc_str(string name, const SpecPrecalc &spec) const;
		double get_splineval(const ElementRef &er, const vector <double> &value, const vector <double> &cval) const;
		vector <double> calculate_precalc_init(const SpecPrecalc &spec_precalc) const;
		const vector <unsigned int>& get_list_time(unsigned int cit, const SpecPrecalc &spec) const;
		void calculate(const SpecPrecalc &spec_calc, PV &param_val, bool store) const;
		void add_eqn(vector <Calculation> &calc, const vector <unsigned int> &param_vec_refq, const vector <unsigned int> &spline_refq, SpecPrecalc &spec_precalc);
		void add_eqn_simp(vector <Calculation> &calc, const vector <unsigned int> &param_vec_refq, const vector <unsigned int> &spline_refq);
		unsigned int add_param(unsigned int th);
		void calcu_add(PreCalc &ca);
		unsigned int add_spline(unsigned int th, SpecPrecalc &spec_precalc);
		EqItem add(Calculation &ca, SpecPrecalc &spec_precalc);
		unsigned int add_list_time(SpecPrecalc &sprec, const vector <unsigned int> &ltime) const;
		vector <unsigned int> combine_list_time(const vector <unsigned int> &lt1,const vector <unsigned int> &lt2) const;
		void map_remove_spec(vector <bool> &map, const SpecPrecalc &spec) const;
		vector <bool> map_empty() const;
		vector <bool> map_from_spec(const SpecPrecalc &spec) const;
		void map_add_spec(vector <bool> &map, const SpecPrecalc &spec) const;
		const SpecPrecalc spec_from_map(const vector <bool> &map) const;
		const SpecPrecalc spec_from_ti_ref_list(vector <unsigned int> list, const vector < vector <unsigned int> > &ti_ref_list, const vector < vector <unsigned int> > &list_time) const;
		//SpecPrecalc sp_add_map_PC(const SpecPrecalc &sprec_old, vector <bool> map_PC) const;
		void sp_add_test(SpecPrecalc &sprec, unsigned int i, const vector <unsigned int> &ltime, vector <double> &timer) const;
		void sp_add(SpecPrecalc &sprec, unsigned int q, const vector <unsigned int> &ltime) const;
		//void sp_add(SpecPrecalc &sprec, const vector <unsigned int> &q_list, const vector <unsigned int> &ltime) const;
		void print_calc() const;
		//void print_calc_section(unsigned int i_start) const;
		SpecPrecalc combine_pv_spec_precalc(const vector <unsigned int> &param_list_tot) const;
		SpecPrecalc combine_pv_spec_precalc_remove(const vector <unsigned int> &param_list_tot, const SpecPrecalc &remove) const;
		SpecPrecalc combine_spec_precalc(const vector <SpecPrecalc> &spec_list) const;
		SpecPrecalc shrink_sprec(unsigned int ti, const SpecPrecalc &spre_old) const;
		void set_param(SpecPrecalc &set_param_spec_precalc, SpecPrecalc &spec_precalc_after, bool spl_fl) const;
		//vector <bool> calculate_map_reparam_time_dep() const;
		SpecPrecalc calculate_spec_precalc_sample(const SpecPrecalc &spec_precalc) const;
		void set_all_time();
		//SpecPrecalc calculate_spec_precalc_all(const SpecPrecalc &spec_precalc) const;
		double hash_ca_mem() const;
		void hash_off();
		string print_ca(unsigned int i, const PreCalc &ca) const;
		bool combine_multiply_add(vector <Calculation> &calc, vector <unsigned int> &become_Rrecalc) const;
		void adjust_num(vector <Calculation> &calc) const;
		void adjust_it(EqItem &it) const;
	
	private:
		vector <unsigned int> get_vec(const PreCalc &ca) const;
		void calc_time_varying(unsigned int i, const PreCalc &ca, const vector <unsigned int> &list_time, vector <double> &precalc, CalcError &calc_err) const;
		double get_val(const PreEqItem &it, const vector <double> &precalc) const;
		void calc_time_varying_fast(unsigned int i, const PreCalc &ca, const vector <unsigned int> &list_time, vector <double> &precalc, CalcError &calc_err) const;
		void eqn_type_error(PreEqItemType type, unsigned int ref) const;
		
		Hash hash_ca;
	
		vector <unsigned int> all_time;
	
		const vector <SpeciesSimp> &species;       // References the species from the model
		const vector <Spline> &spline;             // References splines from the model
		const vector <Param> &param;               // References the parameters from the model
		const vector <ParamVecEle> &param_vec;     // References the param_vec from the model
		const vector <Population> &pop;            // References the populations from the model
		Constant &constant;                  // Stores constants in the model
		const vector <double> &timepoint;    // References time points
		const Details &details;              // Details
};

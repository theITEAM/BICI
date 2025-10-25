#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "const.hh"

class Precalc                             // Stores information about the model
{
	public:
		vector <PreCalc> calcu;             // Stores precalculation
		
		Precalc(const vector <SpeciesSimp> &species, const vector <Spline> &spline, const vector <ParamVecEle> &param_vec, const vector <Population> &pop, Constant &constant, const vector <double> &timepoint, const Details &details);
		
		//void calculate_all(const vector <unsigned int> &list_recalc, PV &param_val) const;
		double get_splineval(const ElementRef &er, const vector <double> &value, const vector <double> &cval) const;
		//void calc_spline_const(PV &param_val, const vector <unsigned int> &spline_ref) const;
		vector <double> calculate_precalc_init(const SpecPrecalc &spec_precalc) const;
		void calculate(const SpecPrecalc &spec_calc, PV &param_val, bool store) const;
		bool add_eqn(vector <Calculation> &calc, const vector <unsigned int> &param_vec_ref, const vector <unsigned int> &spline_ref, SpecPrecalc &spec_precalc, PrecalcAddType add_type=PRECALC_ALL);
		unsigned int add_param(unsigned int th);
		unsigned int add_spline(unsigned int th, SpecPrecalc &spec_precalc);
		EqItem add(PreCalc &pcalc, SpecPrecalc &spec_precalc);
		void sp_add(SpecPrecalc &sprec, unsigned int i, const vector <unsigned int> &ltime) const;
		void sp_add(SpecPrecalc &sprec, const vector <unsigned int> &i_list, const vector <unsigned int> &ltime) const;
		void print_calc() const;
		SpecPrecalc combine_spec_precalc(const vector <unsigned int> &param_list_tot) const;
		SpecPrecalc shrink_sprec(const vector <unsigned int> &lt, SpecPrecalc spre) const;
		void set_param(SpecPrecalc &set_param_spec_precalc, SpecPrecalc &spec_precalc_after, bool spl_fl) const;
		SpecPrecalc calculate_spec_precalc_sample(const SpecPrecalc &spec_precalc) const;
		void set_all_time();
		SpecPrecalc calculate_spec_precalc_all(const SpecPrecalc &spec_precalc, const vector <SpecPrecalcTime> &spec_precalc_time) const;
		
	private:
		vector <unsigned int> get_vec(const PreCalc &ca) const;
		void print_ca(unsigned int i, const PreCalc &ca) const;
	
		Hash hash_ca;
	
		vector <unsigned int> all_time;
	
		const vector <SpeciesSimp> &species;       // References the species from the model
		const vector <Spline> &spline;             // References splines from the model
		const vector <ParamVecEle> &param_vec;     // References the param_vec from the model
		const vector <Population> &pop;            // References the populations from the model
		Constant &constant;                  // Stores constants in the model
		const vector <double> &timepoint;    // References time points
		const Details &details;              // Details
};

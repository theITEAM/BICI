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
		
		void calculate_all(const vector <unsigned int> &list_recalc, PV &param_val) const;
		void calculate(const vector <unsigned int> &list_recalc, const vector <unsigned int> &list_time, PV &param_val, bool store) const;
		bool add_eqn(vector <Calculation> &calc, const vector <double> &param_vec_ref, const vector <double> &spline_ref, vector <unsigned int> &list_precalc, PrecalcAddType add_type=PRECALC_ALL);
		unsigned int add_param(unsigned int th);
		unsigned int add_spline(unsigned int th, vector <unsigned int> &list_precalc);
		EqItem add(PreCalc &pcalc, vector <unsigned int> &list_precalc);
		void print_calc() const;
		
	private:
		vector <unsigned int> get_vec(const PreCalc &ca) const;
		void print_ca(unsigned int i, const PreCalc &ca) const;
	
		Hash hash_ca;
	
		const vector <SpeciesSimp> &species;       // References the species from the model
		const vector <Spline> &spline;             // References splines from the model
		const vector <ParamVecEle> &param_vec;     // References the param_vec from the model
		const vector <Population> &pop;            // References the populations from the model
		Constant &constant;                  // Stores constants in the model
		const vector <double> &timepoint;    // References time points
		const Details &details;              // Details
};

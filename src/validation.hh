#pragma once

#include <string>

using namespace std;

#include "model.hh"
#include "output.hh"
#include "mpi.hh"

class Validation                                // Stores information about the model
{
	public:
		Validation();
		void tornado_setup(string add_info, Operation mode, ExtFactor ext_factor, string file, const vector <string> &data_sim_lines);
		void tornado_result(Operation mode, ExtFactor ext_factor, bool no_question, string file);
		void scan_setup(string scan_info, Operation mode, ExtFactor ext_factor, string file, const vector <string> &data_sim_lines, bool test);
		void scan_result(string scan_info, Operation mode, ExtFactor ext_factor, bool no_question, string file);
		
	private:
		string gp_header() const;
		string gp_header_multi() const;
		string convert_label(string name) const;
		double get_param_value(string na, const Model &model_base) const;
		ParamRef get_param_ref(string na, const Model &model) const;
		ScanInfo get_scan_info(string scan_info, const Model &model) const;
		vector <string> get_file_rep(unsigned int n, string run_dir) const;
		void run_sim_data(string name, string op_dir, bool test, vector <string> file_rep, unsigned int ncore_inf, ExtFactor ext_factor, const vector <string> &data_sim_lines, string file, string com) const;
};

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
		SpecInfo load_spec_info(SpecInfoType type, string spec_file);
		void tornado_setup(string spec_file, ExtFactor ext_factor, string file);
		vector <Table> consistent_row(const vector <Table> &list) const;
		string us(double val) const;
		void tornado_result(string spec_file, ExtFactor ext_factor, bool no_question, string file);
		void scan_setup(string spec_file, ExtFactor ext_factor, string file);
		bool is_log(const SpecInfo &si, vector <Stat> &stat) const;
		string get_si_name(const SpecInfo &si) const;
		void scan_result(string spec_file, ExtFactor ext_factor, bool no_question, string file);
		
	private:
		string gp_header() const;
		string gp_header_multi() const;
		string convert_label(string name) const;
		bool is_tt_stat(string na) const;
		double get_param_value(string na, const Model &model_base) const;
		void delete_files(const Model &model, Output &output, string tor_dir, string run_dir, string op_dir, string spec_file);
		
		ParamRef get_param_ref(string na, const Model &model) const;
		vector <ValidRun> get_valid_run(unsigned int n, string run_dir, const vector <string> &data_sim_lines) const;
		void run_sim_data(string name, string op_dir, bool para_run, const vector <ValidRun> &valid_run, unsigned int ncore_inf, ExtFactor ext_factor, string file, string com, string spec_file, const ExtinctInfo &extinct_info) const;
		ExtinctInfo load_extinction(string te) const;
};

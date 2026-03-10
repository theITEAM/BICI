#pragma once

#include <vector>

using namespace std;

#include "struct.hh"
#include "model.hh"
#include "output.hh"

class DataSim                              // Simulates data from state
{
	public:
		DataSim(Model &model, Output &output, bool sup_=false);
		void run(string data_sim_line);
		void del(string data_sim_line);
		void show() const;
		void clear(string data_sim_line);
		
	private:
		CommandLine line;
		string tag_find;
		string data_sim_line_store;
		
		string get_tag_value(string st);
		void cannot_find_tag() const;
		void error(string msg) const;
		unsigned int option_error(string na, string te, const vector <string> &pos, const vector <unsigned int> &conv) const;
		void init_pop_sim(unsigned int p, const State &state, Table &tab);
		string output_table(const Table &tab) const;
		string source_description(const Species &sp, const DataSource &so) const;
		vector <DataRef> get_data_list() const;
		vector <double> get_times();
		double get_frac();
		ObsModel load_obs_model();
		Stratify get_stratify(const Filter &cf, unsigned int p) const;
		vector <double> get_comp_prob(const Filter &cf, double tdiv, unsigned int cl_sel, unsigned int p) const;
		void set_filt_strat(unsigned int com, Filter &cf, Stratify strat, unsigned int p) const;
		double obs_model_sample(double mean, const ObsModel &om) const;
		vector <double> get_tr_prob(const vector <string> &trans_filt, double tdiv) const;
		void add_ind_sim(unsigned int p, const State &state, Table &tab);
		void copy_from_sim(unsigned int p, string st, Command cname, Table &tab) const;
		void comp_data_sim(unsigned int p, const State &state, Table &tab);
		void trans_data_sim(unsigned int p, const State &state, Table &tab);
		void test_data_sim(unsigned int p, const State &state, Table &tab);
		void pop_data_sim(unsigned int p, const State &state, Table &tab);
		void pop_trans_data_sim(unsigned int p, const State &state, Table &tab);
		void ind_effect_data_sim(unsigned int p, const State &state, Table &tab);
		void ind_group_data_sim(unsigned int p, Table &tab);
		void genetic_data_sim(const State &state, Table &tab);
		void simulate_genetic_snp(double mut_rate, double seq_var, string SNP_root, const State &state, Table &tab) const;
		void simulate_genetic_matrix(double mut_rate, double seq_var, const State &state, Table &tab);
		bool wildcard_match(const string &name, string filt) const;
		bool wildcard_match2(const string &name, const vector <string> &spl, bool root) const;
		vector < vector < vector <double> > > get_ind_times(const vector <double> &tims, const vector <IndList> &ind_list, const State &state) const;
		void print_inf_node(const vector <InfNode> &inf_node, const State &state) const;
		
		bool sup;                              // Determines if output is supressed
		
		Model &model;
		Output &output;
};

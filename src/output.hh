#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

#include "input.hh"
#include "model.hh"
#include "state.hh"
#include "struct.hh"
#include "mpi.hh"

class Output                               // Stores information about the data
{
	public:
		string diagdir;                        // The diagnostics directory (within Data if available)
		string sampledir;                      // The samples directory (within Data if available)
		
		vector <double> timer;                 // General purpose timers
		
		Output(const Model &model, const Input &input, Mpi &mpi);
		void check_open(ofstream &fout, string file) const;
		void summary(const Model &model) const;
		void prop_summary(string te) const;
		string print_affect_like(const AffectLike &al) const;
		void data_summary(const Model &model) const;
		void ensure_directory(const string &path) const;
		void print_initc(const vector <InitCondValue> &initc) const;
		//void updated_file(string file);
		void constants(string file);
		string to_str(double num) const;
		void param_sample(unsigned int s, unsigned int chain, State &state);
		void param_sample(const Particle &part);
		string param_output(const Particle &part, const vector < vector <double> > &value) const;
		void state_sample(unsigned int s, unsigned int chain, State &state);
		void state_sample(const Particle &part);
		string state_output(const Particle &part,	vector <string> &ind_key, Hash &hash_ind) const;
		string trace_init() const;
		string ic_output_head() const;
		string ic_output(const Particle &part) const;
		void set_output_burnin(double burnin_frac);
		void set_diagnostics(unsigned int ch, string diag);
		void end(string file, unsigned int total_cpu) const;
		vector <Particle> get_part_chain(unsigned int chain, const vector <Particle> &part) const;
		string print_obs_data(unsigned int p, const vector <ObsData> &obs) const;
		void print_individuals(unsigned int N, unsigned int p, const State &state) const;
		void final_time(long sec, long op_sec) const;
		void final_memory_usage() const;
	
	private:
		vector < vector <double> > param_value_from_vec(const vector <double> &param_val_prop) const;
		string generate_state_head(const vector <string> &ind_key) const;
		string output_param(const vector < vector <double> > &value, const Particle &part) const;
		vector < vector <double> > param_value_from_value() const;
		string transtype_text(TransType type) const;
		void number_part(vector <Particle> &part) const;
		void output_trace(unsigned int ch, const vector <Particle> &part, vector < vector < vector < vector <double> > > > &param_samp, ofstream &fout) const;
		string get_effective_sample_size(vector <double> vec) const;
		string get_Gelman_Rubin_statistic(const vector < vector <double> > &cha) const;
		string param_stat(unsigned int th, unsigned int i, const vector < vector < vector < vector <double> > > > &param_samp, vector <Warn> &ESS_warn, vector <Warn> &GR_warn) const;
		void output_param_statistics(const vector < vector < vector < vector <double> > > > &param_samp, ofstream &fout, vector <string> &final_warning) const;
		void output_state(unsigned int ch, const vector <Particle> &part, ofstream &fout) const;
		void output_trans_diag(const vector <TransDiagSpecies> &trans_diag, ofstream &fout) const;
		void add_warning(string err_msg, ofstream &fout) const;
		void output_add_ind_warning(vector <string> &final_warning) const;
		void output_rate_warning(unsigned int total_cpu, unsigned int per_start, unsigned int per_end, vector <string> &final_warning) const;
		void output_generation(const vector <Particle> &part, ofstream &fout) const;
		void output_diagnostic(const vector <Diagnostic> &diagnostic, bool &alg_warn_flag,ofstream &fout) const;
		vector <TransDiagSpecies> trans_diag_init() const;
		void trans_diag_add(vector <TransDiagSpecies> &trans_diag, const vector <Particle> &part) const;
		Stat get_statistic(vector <double> &vec) const;
		string generation_average(string head, string content) const;
		
		vector <string> lines_raw;    
	
		vector <Particle> param_store;
		vector <Particle> state_store;
		vector <Diagnostic> diagnostic_store;
		
		const Model &model;	
		Mpi &mpi;	
};


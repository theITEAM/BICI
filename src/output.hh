#ifndef BICI__OUTPUT_HH
#define BICI__OUTPUT_HH

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
		string outputdir;
		
		unsigned int chain;                    // Stores the chain number
		
		vector <long> timer;                   // General purpose timers
		
		Output(unsigned int _chain, const Model &model, const Input &input, Mpi &mpi);
		void check_open(ofstream &fout, string file) const;
		void summary(const Model &model) const;
		void prop_summary(string te) const;
		string print_affect_like(const AffectLike &al) const;
		void data_summary(const Model &model) const;
		void ensure_directory(const string &path) const;
		void print_param(const vector <double> &vec) const;
		void print_initc(const vector <InitCondValue> &initc) const;
		void updated_file(string file);
		void constants(string file);
		string to_str(double num) const;
		void param_sample(unsigned int s, const State &state);
		void param_sample(unsigned int s, const Particle &particle);
		void state_sample(unsigned int s, const State &state);
		void state_sample(unsigned int s, const Particle &part);
		void trace_init();
		void ic_output_head();
		void ic_output(const Particle &part);
		void set_output_burnin(double burnin_frac);
		void generate_files();
		string print_obs_data(unsigned int p, const vector <ObsData> &obs) const;
		void print_individuals(unsigned int N, unsigned int p, const State &state) const;
	
	private:
		vector < vector <double> > param_value_from_vec(const vector <double> &param_val) const;
		string output_param(const vector < vector <double> > &value, const Particle &part) const;
		vector < vector <double> > param_value_from_value() const;
		string transtype_text(TransType type) const;
		
		vector <string> lines_raw;    
	
		stringstream param_out;
		stringstream state_out;
	
		vector <string> ind_key;
		Hash hash_ind;
		
		const Model &model;	
		Mpi &mpi;	
};
#endif

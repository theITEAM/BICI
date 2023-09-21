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

class Output                                             // Stores information about the data
{
	public:
		string outputdir;
		
		Output(unsigned int chain, const Model &model, const Input &input);
		void check_open(ofstream &fout, string file) const;
		void summary(const Model &model) const;
		string print_affect_like(const AffectLike &al) const;
		void data_summary(const Model &model) const;
		void ensure_directory(const string &path) const;
		void print_param(const vector <double> &vec) const;
		void print_initc(const vector <InitCondValue> &initc) const;
		void sample(const State &state, bool include_state, string file);
		void updated_file(string file) const;
		void constants(string file);
		void param_sample(unsigned int s, const State &state);
		void param_sample(unsigned int s, const Particle &particle);
		void state_sample(unsigned int s, const State &state);
		void state_sample(unsigned int s, const Particle &part);
		void trace_init(unsigned int chain);
		void generate_files();
		
	private:
		vector < vector <double> > param_value_from_vec(const vector <double> &param_val) const;
		string output_param(const vector < vector <double> > &value, const Particle &part) const;
		vector < vector <double> > param_value_from_value() const;
		string transtype_text(TransType type) const;
		
		vector <string> lines_raw;    
	
		string param_out_file;
		string state_out_file;
		stringstream param_out;
		stringstream state_out;
		//ofstream param_out;
		//ofstream state_out;
		
		const Model &model;	
};
#endif

		
#pragma once

#include "struct.hh"
#ifdef USE_MPI
#include "mpi.h"
#endif
#include "model.hh"

#include "utils.hh"
#include <fstream>

struct Mpi {
	Mpi(unsigned int core_spec, const Model &model);
	
public:
	unsigned int ncore;                                          // The number of cores that MPI is using
	unsigned int core;                                           // The core of the current process
	bool core_spec_on;                                           // Determines if core has been specified
	
public:
	void send_particle(unsigned int co, const Particle &part);
	void get_particle(unsigned int co, Particle &part);
	void transfer_particle(vector <Particle> &part);
	void transfer_prop_info(vector <PropInfo> &prop_info);
	void transfer_samp(vector < vector <double> > &samp);
	void transfer_diagnostic(vector <Diagnostic> &diag);
	void transfer_terminal_info(vector <TerminalInfo> &term_info);
	void share_particle(vector <Particle> &part);
	void transfer_lines_raw(vector <string> &lines_raw);
	vector <double> gather(const vector <double> &vec);
	vector <unsigned int> gather(const unsigned int val);
	void bcast(unsigned int &val);
	void bcast(double &val);
	void bcast(vector <unsigned int> &vec);
	void bcast(vector <double> &vec);
	void bcast(vector <string> &vec);
	void distribute(vector < vector <double> > &samp);
	void distribute(vector <PropInfo> &prop_info);
	void barrier() const;
	void sample_barrier(unsigned int s, unsigned int nsample) const;
	void mess(string te) const;
	void sum(vector < vector <double> > &array);
	void sum(double &val);
	
private:
	vector<double> buffer;                                       // Stores packed up information to be sent between cores
	unsigned int k;                                              // Indexes the buffer
	
	void pack_initialise();                                      // Initilaises sending information
	void pack_particle(const vector <Particle> &part);
	void unpack_particle(vector <Particle> &part);
	
	void pack_num(unsigned int t);
	unsigned int unpack_num();
	
	template <class T>                                           // Generic pack/unpack commands
	void pack_item(T t);
	template <class T>
	void pack_item(const vector<T>& vec);	
	void pack_item(const string& vec);
	template<class T>
	void unpack_item(T &num);
	template <class T>
	void unpack_item(vector<T>& vec);
	void unpack_item(string &vec);
	
	void pack(const vector <string> &vec_str);                    // Specific pack/unpack commands
	void unpack(vector <string> &vec_str);
	void pack(const Particle &pa);
	void unpack(Particle &pa);
	void pack(const PropInfo &ps);
	void unpack(PropInfo &ps);
	
	size_t packsize();
	double *packbuffer();
	void unpack_check();
	void pack_send(const unsigned int co);
	void pack_recv(const unsigned int co);
	void pack_initialise(const size_t size);
	
	const Model &model;
};


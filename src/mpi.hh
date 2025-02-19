#ifndef BICI__MPI_HH
#define BICI__MPI_HH

#include "struct.hh"
#ifdef USE_MPI
#include "mpi.h"
#endif
#include "model.hh"
#include "utils.hh"
#include <fstream>

struct Mpi {
	Mpi(const Model &model);
	
public:
	unsigned int ncore;                                          // The number of cores that MPI is using
	unsigned int core;                                           // The core of the current process
	
public:
	void transfer_lines_raw(vector <string> &lines_raw);

private:
	vector<double> buffer;                                       // Stores packed up information to be sent between cores
	unsigned int k;                                              // Indexes the buffer
	
	void pack_initialise();                                      // Initilaises sending information
	
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
	
	size_t packsize();
	double *packbuffer();
	void unpack_check();
	void pack_send(const unsigned int co);
	void pack_recv(const unsigned int co);
	void pack_initialise(const size_t size);
	
	const Model &model;
};
#endif

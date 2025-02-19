/// Information and routines for transferring data between cores using MPI

#include <sstream>

using namespace std;

#include "mpi.hh"

Mpi::Mpi(const Model &model): model(model)
{
#ifdef USE_MPI
	int num;
	MPI_Comm_size(MPI_COMM_WORLD,&num); ncore = (unsigned int) num;
  MPI_Comm_rank(MPI_COMM_WORLD,&num); core = (unsigned int) num;
#else
	ncore = 1;
	core = 0;
#endif
}


#ifdef USE_MPI
/// Transfers lines_raw from other core to zero
void Mpi::transfer_lines_raw(vector <string> &lines_raw)
{
	if(ncore == 1) return;
	
	if(core != 0){
		pack_initialise();
		pack(lines_raw);
		pack_send(0);
	}
	else{
		for(auto co = 1u; co < ncore; co++){
			pack_recv(co);		
			vector <string> lines_raw_copy;
			unpack(lines_raw_copy);
			
			lines_raw.push_back("");
			for(const auto &li : lines_raw_copy){
				lines_raw.push_back(li);
			}
		}
	}
}


/// Initialises the buffer
void Mpi::pack_initialise()
{
	k = 0;
	buffer.resize(0);
}


/// Functions to pack up information

template <class T>
void Mpi::pack_item(T t)
{
	buffer.push_back(t); k++;
}

// Provide overloads for cases where default behaviour needs to be modified
void Mpi::pack_item(const string &vec)
{
	auto jmax = vec.length();
	pack_item(jmax);
	for(auto j = 0u; j < jmax; j++){
		pack_item(vec.at(j));
	}
}


/// Use template to share implementation details between cases
template <class T>
void Mpi::pack_item(const vector<T> &vec)
{
	pack_item(vec.size());
	for (auto &item : vec) {
		pack_item(item);
	}
}


/// Functions for unpacking informationtemplate<class T>
template <class T>
void Mpi::unpack_item(T &num)
{
	num = buffer[k]; k++;
}

void Mpi::unpack_item(string &vec)
{
	unsigned int jmax;
	
	unpack_item(jmax);
	stringstream ss; for(auto j = 0u; j < jmax; j++){ ss << (char) buffer[k]; k++;}
	vec = ss.str();
}


template <class T>
void Mpi::unpack_item(vector<T> &vec)
{
	unsigned int size;
	unpack_item(size);
	vec.resize(size);
	for (auto &item : vec) {
		unpack_item(item);
	}
}



void Mpi::pack(const vector <string> &vec_str)
{
	pack_item(vec_str);
}


void Mpi::unpack(vector <string> &vec_str)
{
	unpack_item(vec_str);
}

/// Checks all the values have been read from the buffer
void Mpi::unpack_check()
{
	if(k != buffer.size()) emsg("Mpi");
}


/// Returns the size of the buffer
size_t Mpi::packsize()
{
	return k;
}


/// The pointer to the buffer
double* Mpi::packbuffer()
{
	return &buffer[0];
}



/// Sends the buffer to core co
void Mpi::pack_send(const unsigned int co)
{
	unsigned int si = packsize();
	MPI_Send(&si,1,MPI_UNSIGNED,co,0,MPI_COMM_WORLD);
	MPI_Send(packbuffer(),si,MPI_DOUBLE,co,0,MPI_COMM_WORLD);
}


/// Recieves a message from core co and places it into the buffer
void Mpi::pack_recv(const unsigned int co)
{
	unsigned int si;
	MPI_Recv(&si,1,MPI_UNSIGNED,co,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	pack_initialise(si);
	MPI_Recv(packbuffer(),si,MPI_DOUBLE,co,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}	


/// Initialises the buffer
void Mpi::pack_initialise(const size_t size)
{
	k = 0;
	buffer.resize(size);
}
#endif


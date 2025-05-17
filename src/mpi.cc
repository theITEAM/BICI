/// Information and routines for transferring data between cores using MPI

#include <sstream>
#include <iostream>

using namespace std;

#include "mpi.hh"

Mpi::Mpi(unsigned int core_spec, const Model &model): model(model)
{
	core_spec_on = false;
	
#ifdef USE_MPI
	if(core_spec != UNSET) emsg("Core should not be specified");
	int num;
	MPI_Comm_size(MPI_COMM_WORLD,&num); ncore = (unsigned int) num;
  MPI_Comm_rank(MPI_COMM_WORLD,&num); core = (unsigned int) num;
#else
	ncore = 1;
	core = 0;
	if(core_spec != UNSET){
		core = core_spec;
		core_spec_on = true;
	}
#endif
}


#ifdef USE_MPI

/// Sends a particle to another core
void Mpi::send_particle(unsigned int co, const Particle &part)
{
	pack_initialise();
	pack(part);
	pack_send(co);
}


/// Gets a particle from another core
void Mpi::get_particle(unsigned int co, Particle &part)
{
	pack_recv(co);
	unpack(part);
	unpack_check();
}


/// Transfers particles from other cores to zero
void Mpi::transfer_particle(vector <Particle> &part)
{
	if(ncore == 1) return;
	
	if(core != 0){
		pack_initialise();
		pack_particle(part);
		pack_send(0);
	}
	else{
		for(auto co = 1u; co < ncore; co++){
			pack_recv(co);	
			unpack_particle(part);
		}
	}
}


/// Transfers diagnostic description from other cores to zero
void Mpi::transfer_diagnostic(vector <Diagnostic> &diag)
{
	if(ncore == 1) return;
	
	if(core != 0){
		pack_initialise();
		pack_num(diag.size());
		for(auto i = 0u; i < diag.size(); i++){
			pack_item(diag[i].ch);
			pack_item(diag[i].te);
		}
		pack_send(0);
	}
	else{
		for(auto co = 1u; co < ncore; co++){
			pack_recv(co);	
			auto N = unpack_num();
			for(auto i = 0u; i < N; i++){
				Diagnostic di; 
				unpack_item(di.ch);
				unpack_item(di.te);
				diag.push_back(di);
			}
			unpack_check();
		}
	}
}


/// Packs up a vector of particles
void Mpi::pack_particle(const vector <Particle> &part)
{
	pack_num(part.size());	
	for(auto i = 0u; i < part.size(); i++) pack(part[i]);
}


/// Unpacks a vector of particles and adds to part
void Mpi::unpack_particle(vector <Particle> &part)
{
	auto npart = unpack_num();
	for(auto i = 0u; i < npart; i++){
		Particle pa;
		unpack(pa);
		part.push_back(pa);
	}
	unpack_check();
}
	
	
/// Transfers particles from other cores to zero and then shares with all cores
void Mpi::share_particle(vector <Particle> &part)
{
	transfer_particle(part);
	
	pack_initialise();
	if(core == 0) pack_particle(part);
	
	unsigned int N = buffer.size();
	bcast(N);
	
	if(core != 0) buffer.resize(N);
	
	bcast(buffer);
	
	if(core != 0){
		part.clear();
		unpack_particle(part);
	}		
}


/// Packs up a particle
void Mpi::pack(const Particle &pa)
{
	pack_item(pa.param_val);
	
	// Species
	pack_num(pa.species.size());
	for(auto sp = 0u; sp < pa.species.size(); sp++){
		const auto &sp_val = pa.species[sp];
		const auto &icv = sp_val.init_cond_val;
		pack_item(icv.cnum);
		pack_item(icv.N_focal);
		pack_item(icv.N_focal_unobs);
		pack_item(icv.cnum_reduce);
		pack_item(icv.frac_focal);
		pack_item(icv.frac_comb);
		pack_item(icv.N_total);
		pack_item(icv.N_total_unobs);
		pack_item(icv.frac);
	
		pack_item(sp_val.trans_num);
		
		pack_item(sp_val.nindividual);
		
		pack_num(sp_val.individual.size());
		for(auto i = 0u; i < sp_val.individual.size(); i++){
			const auto &ind = sp_val.individual[i];
			pack_item(ind.type);
			pack_item(ind.name);
			pack_num(ind.ev.size());
			for(auto e = 0u; e < ind.ev.size(); e++){
				const auto &eve = ind.ev[e];
				pack_item(eve.type);
				pack_item(eve.tr_gl);
				pack_item(eve.cl);
				pack_item(eve.move_c);
				pack_item(eve.c_after);
				pack_item(eve.t);
				pack_item(eve.Li);
				pack_item(eve.Li_bp);
				pack_item(eve.observed);
				pack_item(eve.ind_inf_from.p);
				pack_item(eve.ind_inf_from.i);
				pack_item(eve.ind_inf_from.pref);
				pack_item(eve.ind_inf_from.po);
				pack_item(eve.ind_inf_from.w);
				pack_item(eve.inf_node_ref);
				pack_item(eve.m);
				pack_item(eve.ti);
				pack_item(eve.index);
				pack_item(eve.e_origin);
			}
			pack_item(ind.ie);
			pack_item(ind.exp_ie);
			pack_item(ind.X);
			pack_item(ind.exp_fe);
			
			pack_num(ind.popnum_ind_ref.size());
			for(auto j = 0u; j < ind.popnum_ind_ref.size(); j++){
				const auto &pir = ind.popnum_ind_ref[j];
				pack_item(pir.po);
				pack_item(pir.ti);
				pack_item(pir.index);
			}
		
			pack_num(ind.incomp_ref.size());
			for(auto j = 0u; j < ind.incomp_ref.size(); j++){
				const auto &ir = ind.incomp_ref[j];
				pack_item(ir.on);
				pack_item(ir.n);
				pack_item(ir.ti);
				pack_item(ir.dt);
				pack_item(ir.e_begin);
				pack_item(ir.t_end);
				pack_item(ir.index);
				pack_item(ir.Li);
			}

			pack_item(ind.init_c_set);
		}

		pack_item(sp_val.exp_num);
		pack_item(sp_val.cum_prob_dist);
	}
	
	pack_num(pa.dir_out.size());
	for(auto i = 0u; i < pa.dir_out.size(); i++){
		pack_item(pa.dir_out[i].value_str);
	}
	
	const auto &li = pa.like;
	pack_item(li.init_cond);
	pack_item(li.init_cond_prior);
	pack_item(li.obs);
	pack_item(li.prior);
	pack_item(li.spline_prior);
	pack_item(li.dist);
	pack_item(li.markov);
	pack_item(li.nm_trans);
	pack_item(li.genetic_process);
	pack_item(li.genetic_obs);
	pack_item(li.ie);

	const auto &tts = pa.trans_tree_stats;
	pack_item(tts.N_origin);
	pack_item(tts.N_inf);
	pack_item(tts.N_mut_tree);
	pack_item(tts.N_mut_origin);
	pack_item(tts.N_unobs);
	pack_item(tts.t_root);

	pack_num(pa.inf_origin.size());
	for(auto i = 0u; i < pa.inf_origin.size(); i++){
		const auto &io = pa.inf_origin[i];
		pack_item(io.node);
		pack_item(io.mut_num);
	}

	pack_num(pa.inf_node.size());
	for(auto i = 0u; i < pa.inf_node.size(); i++){
		const auto &in = pa.inf_node[i];
		pack_item(in.t_start);
		pack_item(in.t_rec);
		pack_item(in.p);
		pack_item(in.i);
		pack_item(in.e);
		pack_item(in.from.node);
		pack_item(in.from.index);
		
		pack_num(in.inf_ev.size());
		for(auto j = 0u; j < in.inf_ev.size(); j++){
			const auto &ie = in.inf_ev[j];
			pack_item(ie.t);
			pack_item(ie.type);
			pack_item(ie.index);
			pack_item(ie.mut_num);
		}
	}
	
	pack_item(pa.w);
	pack_item(pa.s);
	pack_item(pa.chain);
}


/// Packs up a particle
void Mpi::unpack(Particle &pa)
{
	unpack_item(pa.param_val);

	auto S = unpack_num();
	for(auto sp = 0u; sp < S; sp++){
		ParticleSpecies sp_val;
		auto &icv = sp_val.init_cond_val;
		unpack_item(icv.cnum);
		unpack_item(icv.N_focal);
		unpack_item(icv.N_focal_unobs);
		unpack_item(icv.cnum_reduce);
		unpack_item(icv.frac_focal);
		unpack_item(icv.frac_comb);
		unpack_item(icv.N_total);
		unpack_item(icv.N_total_unobs);
		unpack_item(icv.frac);
	
		unpack_item(sp_val.trans_num);
		
		unpack_item(sp_val.nindividual);
		
		auto N = unpack_num();
		for(auto i = 0u; i < N; i++){
			Individual ind;
			ind.type = (IndType)(unpack_num());
			unpack_item(ind.name);
			auto E = unpack_num();
			for(auto e = 0u; e < E; e++){
				Event eve;
				eve.type = (EventType) (unpack_num());
				unpack_item(eve.tr_gl);
				unpack_item(eve.cl);
				unpack_item(eve.move_c);
				unpack_item(eve.c_after);
				unpack_item(eve.t);
				unpack_item(eve.Li);
				unpack_item(eve.Li_bp);
				unpack_item(eve.observed);
				unpack_item(eve.ind_inf_from.p);
				unpack_item(eve.ind_inf_from.i);
				unpack_item(eve.ind_inf_from.pref);
				unpack_item(eve.ind_inf_from.po);
				unpack_item(eve.ind_inf_from.w);
				unpack_item(eve.inf_node_ref);
				unpack_item(eve.m);
				unpack_item(eve.ti);
				unpack_item(eve.index);
				unpack_item(eve.e_origin);
				
				ind.ev.push_back(eve);
			}
			unpack_item(ind.ie);
			unpack_item(ind.exp_ie);
			unpack_item(ind.X);
			unpack_item(ind.exp_fe);
			
			auto J = unpack_num();
			for(auto j = 0u; j < J; j++){
				PopnumIndRef pir;
				unpack_item(pir.po);
				unpack_item(pir.ti);
				unpack_item(pir.index);	
				ind.popnum_ind_ref.push_back(pir);
			}
		
			auto P = unpack_num();
			for(auto j = 0u; j < P; j++){
				IncompNMTransRef ir;
				unpack_item(ir.on);
				unpack_item(ir.n);
				unpack_item(ir.ti);
				unpack_item(ir.dt);
				unpack_item(ir.e_begin);
				unpack_item(ir.t_end);
				unpack_item(ir.index);
				unpack_item(ir.Li);
				ind.incomp_ref.push_back(ir);
			}
			
			unpack_item(ind.init_c_set);
			
			sp_val.individual.push_back(ind);
		}
		
		unpack_item(sp_val.exp_num);
		unpack_item(sp_val.cum_prob_dist);

		pa.species.push_back(sp_val);
	}
	
	auto N = unpack_num();
	for(auto i = 0u; i < N; i++){
		DeriveOutput dout;
		unpack_item(dout.value_str);
		pa.dir_out.push_back(dout);
	}
	
	auto &li = pa.like;
	unpack_item(li.init_cond);
	unpack_item(li.init_cond_prior);
	unpack_item(li.obs);
	unpack_item(li.prior);
	unpack_item(li.spline_prior);
	unpack_item(li.dist);
	unpack_item(li.markov);
	unpack_item(li.nm_trans);
	unpack_item(li.genetic_process);
	unpack_item(li.genetic_obs);
	unpack_item(li.ie);

	auto &tts = pa.trans_tree_stats;
	unpack_item(tts.N_origin);
	unpack_item(tts.N_inf);
	unpack_item(tts.N_mut_tree);
	unpack_item(tts.N_mut_origin);
	unpack_item(tts.N_unobs);
	unpack_item(tts.t_root);

	auto IO = unpack_num();
	for(auto i = 0u; i < IO; i++){
		InfOrigin io;
		unpack_item(io.node);
		unpack_item(io.mut_num);
		pa.inf_origin.push_back(io);
	}

	auto IN = unpack_num();
	for(auto i = 0u; i < IN; i++){
		InfNode in;
		unpack_item(in.t_start);
		unpack_item(in.t_rec);
		unpack_item(in.p);
		unpack_item(in.i);
		unpack_item(in.e);
		unpack_item(in.from.node);
		unpack_item(in.from.index);
		
		auto IEV = unpack_num();
		for(auto j = 0u; j < IEV; j++){
			InfEvent ie;
			unpack_item(ie.t);
			ie.type = (InfEventType)(unpack_num());
			unpack_item(ie.index);
			unpack_item(ie.mut_num);
			in.inf_ev.push_back(ie);
		}
		
		pa.inf_node.push_back(in);
	}
	
	unpack_item(pa.w);
	unpack_item(pa.s);
	unpack_item(pa.chain);
}


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

/// Packs a number of elements
void Mpi::pack_num(unsigned int t)
{
	buffer.push_back(t); k++;
}


/// Unpacks a number of elements
unsigned int Mpi::unpack_num()
{
	unsigned int t = buffer[k]; k++;
	return t;
}


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
	if(k != buffer.size()){
		emsg("Mpi unpack problem");
	}
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


/// Gathers a double vector across all cores and returns the combined vector to core 0
vector <double> Mpi::gather(const vector <double> &vec)
{
	vector <double> vectot(vec.size()*ncore);
	
	MPI_Gather(&vec[0],vec.size(),MPI_DOUBLE,&vectot[0],vec.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	return vectot;
}


/// Copies a variable in core 0 to all the other cores
void Mpi::bcast(unsigned int &val)
{
	MPI_Bcast(&val,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
}


/// Copies a variable in core 0 to all the other cores
void Mpi::bcast(double &val)
{
	MPI_Bcast(&val,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
}


/// Copies a vector in core 0 to all the other cores
void Mpi::bcast(vector <unsigned int> &vec)
{
	MPI_Bcast(&vec[0],vec.size(),MPI_UNSIGNED,0,MPI_COMM_WORLD);
}


/// Copies a vector in core 0 to all the other cores
void Mpi::bcast(vector <double> &vec)
{
	MPI_Bcast(&vec[0],vec.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
}

/// Copies a vector in core 0 to all the other cores
void Mpi::bcast(vector <string> &vec)
{
	pack_initialise();
	if(core == 0){
		pack_num(vec.size());
		for(auto i = 0u; i < vec.size(); i++) pack_item(vec[i]);
	}
	
	unsigned int B = buffer.size();
	bcast(B);

	if(core != 0) buffer.resize(B);
		
	bcast(buffer);
	
	if(core != 0){
		auto N = unpack_num();
		vec.resize(N);
		for(auto i = 0u; i < N; i++) unpack_item(vec[i]);
		unpack_check();
	}		
}


/// Waits for all processes
void Mpi::barrier() const
{
	MPI_Barrier(MPI_COMM_WORLD);  
}


/// Waits for all processes
void Mpi::mess(string te) const 
{
	MPI_Barrier(MPI_COMM_WORLD);  
	if(core == 0) cout << te << endl;
}
#endif

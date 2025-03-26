// Implements the PAS algorithm

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm> 
 
using namespace std;

#include "pas.hh"
#include "state.hh"
#include "utils.hh"
#include "matrix.hh"

PAS::PAS(const Model &model, Output &output, Mpi &mpi) : model(model), output(output),  mpi(mpi)
{	
	const auto &de = model.details;
	
	auto f = de.burnin_frac/100.0;
	nburnin = f*de.sample;
	nsample = de.sample;
			
	gen_update = de.gen_update;
	
	num_per_core = de.num_per_core;
	output_param = de.output_param;
	output_state = de.output_state;

	nparam_op = 0;
	nstate_op = 0;
	
	for(auto ch = 0u; ch < num_per_core; ch++){
		chain.push_back(Chain(nburnin,nsample,model,output));
	}
	
	phi = SMALL;
	phi_final = 1;
	
	quench_factor = 1;
	
	// Aim for 100 parameter samples per generation
	step_gen_param = (gen_update*de.nchain)/100;
	if(step_gen_param > gen_update) step_gen_param = gen_update;
	if(step_gen_param == 0) step_gen_param = 1;
	
	for(auto &ch : chain) ch.state.init();
}
 
 
/// Runs PAS algorithm
void PAS::run()
{
	//output.percentage(0,nsample);
	
	if(op()) cout << "Annealing..." << endl;
	
	for(auto &ch : chain) ch.state.dif_thresh = DIF_THRESH_BURNIN;
		
	for(auto &ch : chain) ch.init();

	auto time_anneal_start = clock();

	auto step_Lobs = gen_update/10;
	if(step_Lobs == 0) step_Lobs = 1;

	if(com_op) cout << "<ANNEALING>" << endl;
	
	auto g = 0;
	do{
		progress(phi,1);
			
		if(op()) cout << "Generation " << g << ":  φ=" << phi << endl;

		for(auto &ch : chain){ ch.Lobs_av = 0; ch.nLobs_av = 0;}
	 
		for(auto s = 0u; s < gen_update; s++) {   // Updates mcmc
			for(auto &ch : chain){
				ch.pas_burn_update(s,g,gen_update,phi);
				ch.update(s);
			}
			
			if(s%step_Lobs == step_Lobs-1){
				for(auto &ch : chain){
					ch.Lobs_av += ch.like_total_obs(); 
					ch.nLobs_av++;
				}
			}
			
			if(s%step_gen_param == step_gen_param-1){ // Stores samples for generation plot
				for(const auto &ch : chain) output.param_sample(g,GEN_PLOT,ch.state);
			}
		}

		if(phi == phi_final) break;
		
		bootstrap();
		g++;
	}while(true);
	
	auto time_anneal_end = clock();
	
	// Running phase
	if(op()) cout << "Running..." << endl;
	if(com_op) cout << "<RUNNING>" << endl;
	
	long time_start = clock();

	for(auto s = 0u; s < nsample; s++){
		//if(mpi.core == 14) cout << s << " samp" << endl;
	
		output.percentage(s,nsample);
		progress(s,nsample);
	
		for(auto &ch : chain){
			ch.burn_update(s);
			ch.update(s);
		}		
		sample_op(s);
	}
	output.percentage(nsample,nsample);
	
	output.set_output_burnin(double(100.0*nburnin)/nsample);
	
	for(auto ch = 0u; ch < num_per_core; ch++){
		auto diag = chain[ch].diagnostics(clock()-time_start,time_anneal_end-time_anneal_start);
		output.set_diagnostics(mpi.core*num_per_core+ch,diag);
	}
}


/// Determines if a state output is made
void PAS::sample_op(unsigned int s)
{
	auto base = mpi.core*num_per_core;
	
	if(nparam_op <= double(s*output_param)/nsample){
		for(auto ch = 0u; ch < num_per_core; ch++){
			output.param_sample(s,base+ch,chain[ch].state);
		}	
		nparam_op++;
	}
	
	if(s >= nburnin){
		if(nstate_op <= double((s-nburnin)*output_state)/(nsample-nburnin)){
			for(auto ch = 0u; ch < num_per_core; ch++){
				output.state_sample(s,base+ch,chain[ch].state);
			}
			nstate_op++;
		}
	}
}


// Used to order particles by EF
bool PartRO_ord (PartReorder p1,PartReorder p2)                      
{ return (p1.L < p2.L); };  


/// Works out which particles should be copied and sets new inverse temperature
void PAS::bootstrap()
{
	vector <double> L_store;
	
	for(const auto &ch : chain){
		L_store.push_back(ch.Lobs_av/ch.nLobs_av);
	}
	
#ifdef USE_MPI 
	auto L_store_tot = mpi.gather(L_store);
	unsigned int Ntot = L_store_tot.size();
	mpi.bcast(Ntot);
#else
	auto L_store_tot = L_store;
	unsigned int Ntot = L_store_tot.size();
#endif
	
	vector <unsigned int> partcopy(Ntot);
	
	if(op()){
		vector <PartReorder> part_ro(Ntot);
		for(auto i = 0u; i < Ntot; i++){ part_ro[i].i = i; part_ro[i].L = L_store_tot[i];}
		
		sort(part_ro.begin(),part_ro.end(),PartRO_ord);  // Sorts EFs
		
		auto i1 = int(0.5*Ntot);
		auto i2 = int(0.75*Ntot);
		if(i1 == i2){
			if(i1 > 0) i1--;
			else i2++;
		}
		auto L_range = part_ro[i2].L - part_ro[i1].L;
	
		if(L_range == 0) L_range = 1;
	
	  auto dphi = quench_factor/L_range;   // Sets the increase in inverse temperature
		
		if(phi+dphi > phi_final) dphi = phi_final-phi;  // Limits phi to phifinal
		
		auto Lav = 0.0;                           // Calculates the average in EF across runs
		for(auto val : L_store_tot) Lav += val;
		Lav /= Ntot;
		
		vector <double> wsum(Ntot);                             // Weights particles
		vector <unsigned int> num(Ntot,0);   
		auto sum = 0.0;
		for(auto j = 0u; j < Ntot; j++){
			sum += exp(0.5*dphi*(L_store_tot[j]-Lav));
			wsum[j] = sum;
		}
		if(sum == 0) emsg("Zero probability");
		
		for(auto j = 0u; j < Ntot; j++){             // Samples in proportion to particle weights
			auto z = ran()*sum;
			auto sel=0u; while(sel < Ntot && z > wsum[sel]) sel++;
			if(sel == Ntot) emsg("PAS");
			num[sel]++;
		}
		
		vector <unsigned int> list;                  // Works out which particles to copy and which to discard
		for(auto j = 0u; j < Ntot; j++){     
			if(num[j] > 1){ for(auto k = 1u; k < num[j]; k++) list.push_back(j);}
		}
		
		if(false) cout << Ntot - list.size() << " / "<< Ntot << " Particles kept" << endl;
		
		for(auto j = 0u; j < Ntot; j++){  
			if(num[j] > 0) partcopy[j] = UNSET;
			else{
				if(list.size() == 0) emsg("PAS");
				partcopy[j] = list[list.size()-1];
				list.pop_back();
			}
		}
		if(list.size() != 0) emsg("PAS");
		
		phi += dphi;
	}
	
#ifdef USE_MPI 
	mpi.bcast(partcopy);
	mpi.bcast(phi);
#endif
	
	for(auto j = 0u; j < Ntot; j++){  	
		auto from = partcopy[j];
		if(from != UNSET){        	
			auto ch = j%num_per_core;
			if(j/num_per_core == mpi.core){     // Receive 
				if(from/num_per_core == mpi.core){ // Comes from the same core	
					auto ch_from = from%num_per_core;

					auto part = chain[ch_from].state.generate_particle(0,ch,true);
					chain[ch].state.set_particle(part);
				}
				else{
#ifdef USE_MPI
					Particle part;
					mpi.get_particle(from/num_per_core,part);
					chain[ch].state.set_particle(part);
#else
					emsg("Should be only one core");
#endif
				}	
			}
			else{  
				if(from/num_per_core == mpi.core){     // Send
#ifdef USE_MPI
					auto ch_from = from%num_per_core;
					auto part = chain[ch_from].state.generate_particle(0,ch,true);
					mpi.send_particle(j/num_per_core,part);
#else
					emsg("SHould be only one core");
#endif
				}
			}
		}
	}
}

// Implements the PAS algorithm

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm> 
#include <unistd.h>

using namespace std;

#include "pas.hh"
#include "state.hh"
#include "utils.hh"
#include "matrix.hh"
#include "synchronise.hh"

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
	for(auto &ch : chain) ch.state.dif_thresh = DIF_THRESH_BURNIN;
		
	percentage_start(INIT_PER);	
	for(auto ch = 0u; ch < chain.size(); ch++){
		chain[ch].init(ch,chain.size());
	}
	percentage_end();
	
	percentage_start(ANNEAL_PER);	
	
	double time_anneal_start = clock();

	auto step_Lobs = gen_update/10;
	if(step_Lobs == 0) step_Lobs = 1;

	auto g = 0;
	do{
		estimate_percentage();
			
		//if(op()) cout << "Generation " << g << ":  phi=" << phi << endl;

		for(auto &ch : chain){ ch.Lobs_av = 0; ch.Lobs_av2 = 0; ch.nLobs_av = 0;}
	 
		for(auto s = 0u; s < gen_update; s++) {   // Updates mcmc
			//if(op()) cout << s << endl; 
			
			for(auto &ch : chain){
				ch.pas_burn_update(s,gen_update,phi);
				ch.update(s);
			}
			
			if(s%step_Lobs == step_Lobs-1){
				for(auto &ch : chain){
					auto L = ch.like_total_obs(); 
					ch.Lobs_av += L; ch.Lobs_av2 += L*L; 
					ch.nLobs_av++;
				}
			}
			
			if(s%step_gen_param == step_gen_param-1){ // Stores samples for generation plot
				for(auto &ch : chain) output.param_sample(g,GEN_PLOT,ch.state);
			}
		}

		if(phi == phi_final) break;
		
		bootstrap();
		
		if(model.sync_on) synchronise_proposal(UNSET,chain,mpi);
		
		g++;
	}while(true);
	percentage_end();
		
	double time_anneal_end = clock();
	
	// Running phase
	percentage_start(RUN_PER);	
	
	for(auto s = 0u; s < nsample; s++){
		percentage(s,nsample);
	
		if(model.sync_on) synchronise_proposal(s,chain,mpi);
	
		for(auto &ch : chain){
			ch.pas_burn_update_run(s);
			ch.update(s);
		}		
		sample_op(s);
		
#ifdef USE_MPI
		mpi.sample_barrier(s,nsample);
#endif
	}

#ifdef USE_MPI
	mpi.barrier();
#endif
	
	output.set_inference_prop(double(100.0*nburnin)/nsample,"burnin-percent",BURNIN_FRAC_DEFAULT);
	
	for(auto ch = 0u; ch < num_per_core; ch++){
		auto ch_tot = mpi.core*num_per_core+ch;
		const auto &cha = chain[ch];
		
		auto diag = cha.diagnostics(clock()-time_anneal_start,time_anneal_end-time_anneal_start);
		output.set_diagnostics(ch_tot,diag);
		output.terminal_info.push_back(cha.get_terminal_info(ch_tot));
	}
	
	percentage_end();
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
bool PartRO_ord (const PartReorder &p1, const PartReorder &p2)                      
{ return (p1.L < p2.L); };  


/// Works out which particles should be copied and sets new inverse temperature
void PAS::bootstrap()
{
	vector <double> L_av, L_var;
	for(const auto &ch : chain){
		auto av = ch.Lobs_av, av2 = ch.Lobs_av2;
		auto n = ch.nLobs_av;
		L_av.push_back(av/n);
		L_var.push_back(av2/n - (av/n)*(av/n)+TINY);
	}
	
#ifdef USE_MPI 
	auto L_av_tot = mpi.gather(L_av);
	auto L_var_tot = mpi.gather(L_var);
	unsigned int Ntot = L_av_tot.size();
	mpi.bcast(Ntot);
#else
	auto L_av_tot = L_av;
	auto L_var_tot = L_var;
	unsigned int Ntot = L_av_tot.size();
#endif
	
	vector <unsigned int> partcopy(Ntot);
	
	if(op()){
		// Variance within chain
		auto L_var_mean = mean(L_var_tot);
		
		// Variance between chain
		auto av = 0.0, av2 = 0.0, nav = 0.0; 
		for(auto va : L_av){ av += va; av2 += va*va; nav++;}
		auto L_var_bet = (av2/nav) - (av/nav)*(av/nav);
		
		auto L_sd_mean = sqrt(L_var_mean+L_var_bet+TINY);
		if(L_sd_mean == 0) L_sd_mean = 1;
	
	  auto dphi = quench_factor/L_sd_mean;   // Sets the increase in inverse temperature
		
		if(phi+dphi > phi_final) dphi = phi_final-phi;  // Limits phi to phifinal
		
		auto Lav = mean(L_av_tot);         // Calculates the average in EF across runs
	
		vector <double> wsum(Ntot);                             // Weights particles
		vector <unsigned int> num(Ntot,0);   
		auto sum = 0.0;
		for(auto j = 0u; j < Ntot; j++){
			sum += exp(0.5*dphi*(L_av_tot[j]-Lav));
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
		
		phi += dphi; //cout << phi << "phi" << endl;
	}
	
#ifdef USE_MPI 
	mpi.bcast(partcopy);
	mpi.bcast(phi);
#endif
	
	phi_list.push_back(phi);
	
	for(auto j = 0u; j < Ntot; j++){  	
		auto from = partcopy[j];
		if(from != UNSET){        	
			auto ch = j%num_per_core;
			if(j/num_per_core == mpi.core){     // Receive 
				if(from/num_per_core == mpi.core){ // Comes from the same core	
					auto ch_from = from%num_per_core;

					auto part = chain[ch_from].state.generate_particle(0,ch,true,false);
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
					auto part = chain[ch_from].state.generate_particle(0,ch,true,false);
					mpi.send_particle(j/num_per_core,part);
#else
					emsg("Should be only one core");
#endif
				}
			}
		}
	}
}


/// Estimates percentage based on list of phi
void PAS::estimate_percentage() const 
{
	auto N = phi_list.size();
	if(N < 2){
		percentage(0.01*N,1);
	}
	else{ // Performs linear regression to estimate when annealing will finish
		vector <double> x, y;
		for(auto j = 0u; j < N; j++){
			x.push_back(j+1); y.push_back(log(phi_list[j]));
		}
		
		vector <double> w(N);
		auto v = 1.0;
		for(int j = N-1; j >= 0; j--){
			w[j] = v; v *= 0.8;
		}	

		auto x_av = 0.0, y_av = 0.0, nav = 0.0;
		for(auto j = 0u; j < N; j++){
			x_av += w[j]*x[j];
			y_av += w[j]*y[j];
			nav += w[j];
		}
		x_av /= nav; y_av /= nav;
		
		auto sum = 0.0, sum2 = 0.0;
		for(auto j = 0u; j < N; j++){
			auto dx = x[j]-x_av, dy = y[j]-y_av;
			sum += w[j]*dx*dy; sum2 += w[j]*dx*dx;
		}
		auto beta = sum/sum2;
		auto alpha = y_av - beta*x_av;
		
		auto x_final = (unsigned int)(1-alpha/beta);  // The estimates number of generations

		if(false){
			ofstream fout("point.txt");
			for(auto j = 0u; j < N; j++){
				fout << x[j] << " " << y[j] << " " << w[j] << endl;
			}
		
			ofstream fout2("line.txt");
			fout2 << 0 << " " << alpha << endl;
			fout2 << 20 << " " << (alpha + 20*beta) << endl;
		}
		
		auto frac_done = double(N)/x_final;
		if(frac_done > 1) frac_done = 1;
		if(frac_done < 0.6){
			auto Nstart = 10u;
			if(N < Nstart){
				frac_done *= double(N)/Nstart;	
			}
		}
		percentage(frac_done,1);
	}
}
		
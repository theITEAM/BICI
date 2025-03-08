// Implements the basic ABC_SMC rejection algorithm

#include <iostream>
#include <fstream>
#include <algorithm> 
#include <sstream> 
#include <string>
#include <cmath>

using namespace std;

#include "abc_smc.hh"
#include "utils.hh"
#include "matrix.hh"

ABC_SMC::ABC_SMC(const Model &model, Output &output, Mpi &mpi) : model(model), output(output), mpi(mpi), state(model)
{	
	const auto &de = model.details;
	nsample = de.sample;
	acc_frac = de.accfrac;
	G = de.numgen;
	si = de.kernelsize;
	state.init();
}
 

/// Implements a version of simple ABC rejection algorithm
void ABC_SMC::run()
{
	auto prop = intialise_proposal();

	vector <Particle> particle_store;                // Stores the previous generation
	
	double obs_cutoff = -LARGE; 

	auto smax = model.details.num_per_core;
	for(auto g = 0u; g < G; g++){
		vector <Particle> particle;
		
		auto ntr = 0u, nac = 0u;
		
		output.percent_done = 0;
		
		auto last = false; if(g == G-1) last = true;
			
		if(g == 0){
			for(auto s = 0u; s < smax; s++){
				output.percentage(s,smax*G);
				progress(s,smax*G);
				
				auto param_val = model.param_sample();
				auto initc_val = model.initc_sample(param_val);

				state.simulate(param_val,initc_val);
				
				output.param_sample(g,GEN_PLOT,state);
				auto part = state.generate_particle(UNSET,0,last);
				part.w = 1;
				particle.push_back(part);
				ntr++; nac++;
			}
		}
		else{
			if(g == 1) prior_ref = model.prior_total(particle_store[0].param_val);
			
			cor_matrix.set_mvn_from_particle(particle_store);
			
			prop.set_mvn(si,cor_matrix);
		
			setup_particle_sampler(particle_store);         // Sets up sampler for particles
		
			auto s = 0u;
			do{
				output.percentage(s,smax);
				progress(smax*g+s,smax*G);
				
				ntr++;
					
				auto p = particle_sampler();                  // Samples from a particle in last gen
				const auto &part = particle_store[p];
				
				auto param_prop = prop.sample(part.param_val);// Proposes a new parameter set using MVN kernal
				
				if(model.inbounds(param_prop) == true){       // Checks if parameters within bounds
					auto initc_val = model.initc_sample(param_prop);

					state.simulate(param_prop,initc_val);       // Simulates a new state
																									  	// Calculates the weight for the sample
					auto w = calculate_particle_weight(param_prop,particle_store,prop); 

					if(state.like.obs > obs_cutoff){            // Checks if error function less than cutoff		
						output.param_sample(g,GEN_PLOT,state);
						auto part = state.generate_particle(UNSET,0,last);
						part.w = w;
						particle.push_back(part);
						nac++;
						s++;
					}
				}
			}while(s < smax);
		}
		output.percentage(smax,smax);
		
		if(!last) obs_cutoff = implement_cutoff_frac(particle); 
		
		particle_store = particle;
		
#ifdef USE_MPI 
		mpi.share_particle(particle_store);
#endif
		
		if(op()){
			cout << "Generation " << g << "   Acceptance: " << 100*double(nac)/ntr << "%  ";
			cout << "Cut-off: " << obs_cutoff << "   # Particle: " <<  particle_store.size() << endl << endl;
		}
	}

	setup_particle_sampler(particle_store);    
	for(auto s = 0u; s < smax; s++){
		auto p = particle_sampler();   
		const auto &part = particle_store[p];
	
		output.param_sample(part);
		output.state_sample(part);
	}
}
	
	
/// Intialises the proposal distribution
Proposal ABC_SMC::intialise_proposal() const 
{
	vector <unsigned int> vec;
	
	// Sets the vect to all the non-derived and non-fixed variables
	for(auto i = 0u; i < model.nparam_vec; i++){ // Univariate distributions
		const auto &pv = model.param_vec[i];
		const auto &par = model.param[pv.th];
		if(par.variety == PRIOR_PARAM || par.variety == DIST_PARAM){
			if(par.prior[pv.index].type != FIX_PR){
				vec.push_back(i);
			}
		}
	}
	
	Proposal pp(PARAM_PROP,vec,model,output,1,burn_info);
	
	return pp;
}


/// Cuts off the samples which do not agree as well
double ABC_SMC::implement_cutoff_frac(vector <Particle> &particle) const
{
	vector <double> obs_prob;                              // Calculates the cut-off
	for(auto &part : particle){
		obs_prob.push_back(part.like.obs);
	}

#ifdef USE_MPI 
	obs_prob = mpi.gather(obs_prob);
#endif

	sort(obs_prob.begin(),obs_prob.end());
	
	auto cutoff = obs_prob[(obs_prob.size()-1)*(1-acc_frac)];

#ifdef USE_MPI 
	mpi.bcast(cutoff);
#endif

	auto p = 0u;
	while(p < particle.size()){
		if(particle[p].like.obs < cutoff) particle.erase(particle.begin()+p);
		else p++;
	}	
	
	return cutoff;
}


/// Calculates the weights for the different particles
double ABC_SMC::calculate_particle_weight(const vector <double> &param_prop, const vector <Particle> &particle, const Proposal &prop) const 
{	
	auto sum = 0.0;
	for(const auto &part : particle){	
		sum += part.w*exp(prop.mvn_probability(param_prop,part.param_val));
	}

	return exp(model.prior_total(param_prop) - prior_ref)/sum;
}	


/// Generates a sampler to select a particle
void ABC_SMC::setup_particle_sampler(const vector <Particle> &particle)
{
	auto NN = particle.size();
	
	wsum.resize(NN);
	auto sum = 0.0;
	for(auto i = 0u; i < NN; i++){ 
		sum += particle[i].w;
		wsum[i] = sum;
	}
}


/// Samples the next particle
unsigned int ABC_SMC::particle_sampler() const
{
	auto NN = wsum.size();
	
	double z = ran()*wsum[NN-1]; 
	auto p = 0u; while(p < NN && z > wsum[p]) p++;
	if(p == NN) emsg("ABCSMC");
	
	return p;
}

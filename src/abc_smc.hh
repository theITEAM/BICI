#ifndef BICI__ABC_SMC_HH
#define BICI__ABC_SMC_HH

#include "struct.hh"
#include "model.hh"
#include "proposal.hh"
#include "state.hh"
#include "output.hh"

class ABC_SMC
{
	public:
		ABC_SMC(const Model &model, Output &output, Mpi &mpi);
		void run();
	
	private:
		Proposal intialise_proposal() const;
		double implement_cutoff_frac(vector <Particle> &particle) const;
		
		double calculate_particle_weight(const vector <double> &param_val, const vector <Particle> &particle, const Proposal &prop) const;
		void setup_particle_sampler(const vector <Particle> &particle);
		unsigned int particle_sampler() const;
	
		unsigned int nsample;                  // The total number of samples
		
		unsigned int G;                        // The number of generation
		
		double si;                             // The size of the kernel updates
		
		double acc_frac;                       // The acceptance fraction
		
		vector <double> wsum;                  // Used to sample particles
		
		double prior_ref;                      // This is used as a reference for the prior to avoid large numbers
		
		BurnInfo burn_info;                    // This is redundant (used to set up proposal)

		CorMatrix cor_matrix;                  // Stores information about the correlation matrix

		const Model &model;
		Output &output;
		Mpi &mpi;
		
		State state;
};
#endif

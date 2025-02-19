#ifndef BICI__POST_SIM_HH
#define BICI__POST_SIM_HH

#include <vector>

using namespace std;

#include "struct.hh"
#include "model.hh"
#include "state.hh"
#include "output.hh"

class PostSim                              // Simulates from the model
{
	public:
		PostSim(const Model &model, Output &output);
		void run();
		
	private:
		const Model &model;
		Output &output;
		
		State state;
};
#endif

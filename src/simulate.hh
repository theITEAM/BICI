#ifndef BICI__SIMULATE_HH
#define BICI__SIMULATE_HH

#include <vector>

using namespace std;

#include "struct.hh"
#include "model.hh"
#include "state.hh"
#include "output.hh"

class Simulate                             // Simulates from the model
{
	public:
		Simulate(const Model &model, Output &output);
		void run();
		
	private:
		const Model &model;
		Output &output;
		
		State state;
};
#endif

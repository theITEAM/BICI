
#CXX = g++
#CXX = mpicxx
CXX = mpic++
#CXXFLAGS = -O3 -std=c++11  -fmax-errors=3 -g -Wall -Wextra 
CXXFLAGS = -O3 -std=c++14  -fmax-errors=3 -g  -Wall
#CXXFLAGS = -O3 -std=c++14  -fmax-errors=3 -Wall
#CXXFLAGS = -O3 -std=c++11  -fmax-errors=3 -Wall -Wextra 

#CXX = mpicxx
#CXXFLAGS = -O3 -std=c++11 -march=native -MMD -MP   -fmax-errors=3 -g
#CXXLFAGS += -Wall -Wextra 
LIBRARIES = 
# src/update.cc
SOURCES = src/cor_matrix.cc src/chain.cc src/pas.cc src/mpi.cc src/like_integrate.cc src/mfa.cc src/post_sim.cc src/state_local.cc src/proposal_local.cc src/proposal_utils.cc src/state_local.cc src/state_genetic.cc  src/utils_eqn.cc src/source_sampler.cc src/input.cc src/input_commands.cc src/input_modelupdate.cc src/input_utils.cc src/input_check.cc  src/species_data.cc src/state_species_like.cc  src/state_species_local.cc src/state_species_like_obs.cc    src/hash.cc src/abc_smc.cc src/abc.cc src/ind_ev_simulate.cc src/ind_ev_simulate_single.cc src/ind_ev_sampler.cc src/ind_ev_sampler_noobs.cc  src/state_update_ind.cc src/species_update_ind.cc src/state_check.cc src/proposal.cc src/mcmc.cc src/species.cc src/state_species_ind.cc  src/state_species.cc src/matrix.cc src/state_species_check.cc src/output.cc src/simulate.cc src/state.cc src/equation.cc src/equation_linearise.cc src/bici.cc src/model.cc src/utils_check.cc src/utils.cc

# 
 #src/input_data.cc 
 

OBJECTS = $(SOURCES:.cc=.o)

.PHONY: clean all
.DEFAULT_GOAL := all

all: bici-para

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(LIBRARIES) $< -o $@ -c

bici-para: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LIBRARIES) $^ -o $@
	cp bici-para ..

clean:
	rm src/*.o

-include $(OBJECTS:.o=.d)



CXX = g++
#CXX = mpicxx
#CXXFLAGS = -O3 -std=c++14 -fmax-errors=3  -Wall -Wextra -static-libgcc -static-libstdc++
#CXXFLAGS = -O3 -std=c++14 -fmax-errors=3  -Wall -Wextra -g -static-libgcc -static-libstdc++
CXXFLAGS = -O3 -std=c++14 -fmax-errors=3  -Wall -Wextra -static-libgcc -static-libstdc++
#-static-libgcc -static-libstdc++
#CXXFLAGS = -O3

#CXX = mpicxx
#CXXFLAGS = -O3 -std=c++14 -march=native -MMD -MP   -fmax-errors=3 -g
#CXXLFAGS += -Wall -Wextra 
LIBRARIES = 


SOURCES = state_species_ind.cc   consistent.cc  model.cc equation_simplify.cc validation.cc state_profiler.cc lzw.cc input_profiler.cc data_sim.cc extend.cc synchronise.cc hash_simp.cc state_species_like.cc proposal_utils.cc state_species.cc  precalc.cc  equation_calculate.cc equation.cc equation_linearise.cc cor_matrix.cc  chain.cc pas.cc mpi.cc like_integrate.cc mfa.cc post_sim.cc state_local.cc proposal_local.cc state_local.cc state_genetic.cc  utils_eqn.cc source_sampler.cc input.cc input_commands.cc input_modelupdate.cc input_utils.cc input_check.cc  species_data.cc   state_species_local.cc state_species_like_obs.cc    hash.cc abc_smc.cc abc.cc ind_ev_simulate.cc ind_ev_simulate_single.cc ind_ev_sampler.cc ind_ev_sampler_noobs.cc  state_update_ind.cc species_update_ind.cc state_check.cc proposal.cc mcmc.cc species.cc matrix.cc state_species_check.cc output.cc simulate.cc state.cc  bici.cc utils_check.cc utils.cc


OBJECTS = $(SOURCES:.cc=.o)

.PHONY: clean all
.DEFAULT_GOAL := all

all: bici-core

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(LIBRARIES) $< -o $@ -c

bici-core: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LIBRARIES) $^ -o $@
#	mv bici-core.exe ..

clean:
	rm *.o

-include $(OBJECTS:.o=.d)


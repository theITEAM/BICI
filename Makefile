
CXX = g++
#CXXFLAGS = -O0 -std=c++11  -march=native -MMD -MP  -fmax-errors=3 -g
CXXFLAGS = -std=c++11  -fmax-errors=3 -g -Wall -Wextra 

#CXX = mpicxx
#CXXFLAGS = -O3 -std=c++11 -march=native -MMD -MP   -fmax-errors=3 -g
#CXXLFAGS += -Wall -Wextra 
LIBRARIES = 


SOURCES = src/abc_smc.cc src/abc.cc src/ind_ev_sampler.cc  src/state_update_ind.cc src/state_check.cc src/proposal.cc src/mcmc.cc src/state_species_like.cc src/species_data.cc  src/species.cc src/state_species_ind.cc  src/state_species.cc src/matrix.cc src/state_species_check.cc src/model_data.cc src/output.cc src/simulate.cc src/state.cc src/equation.cc src/bici.cc src/model.cc src/utils_check.cc src/utils.cc src/input.cc src/input_commands.cc src/input_modelupdate.cc src/input_utils.cc src/input_check.cc 

 #src/input_data.cc 
 

OBJECTS = $(SOURCES:.cc=.o)

.PHONY: clean all
.DEFAULT_GOAL := all

all: bici

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(LIBRARIES) $< -o $@ -c

bici: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LIBRARIES) $^ -o $@
	cp bici ..

clean:
	rm *.o bici *.d

-include $(OBJECTS:.o=.d)


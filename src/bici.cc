#include <iostream>
#include <sstream>
#include <math.h>
#include <map>
#include <algorithm>
#include <vector>
#include <iterator>
#include <signal.h>

using namespace std;

#include "model.hh"
#include "output.hh"
#include "input.hh"
#include "const.hh"
#include "utils.hh"
#include "simulate.hh"
#include "mcmc.hh"
#include "abc.hh"
#include "abc_smc.hh"

// Compile with make
// Run with ./bici 
// valgrind --leak-check=yes ./bici

// Compiling for windows
// .\make.bat from the directory D:\C_complier_BICI2.0\bin

int main(int argc, char** argv)
{
	//test_distribution(); return 0;
	
	string file = "Execute/init.txt";
	
	auto chain = 0;

	if(argc > 2){ cout << "Too many arguments" << endl; return 0;}
	if(argc == 2) chain = atoi(argv[1]);


	//if(argc == 2) file = argv[1];
	
	Model model;                            // Used to store the model structure 

	Input input(model,file);                // Imports information from input file into model

	Output output(chain,model,input);       // Sets up the class for model outputs

	output.summary(model);                  // Outputs a text file sumarising the model
	
	output.data_summary(model);             // Outputs a text file sumarising the data
	
	//output.constants("constants.txt");      // Creates a file with all model constants
	//return 0;
	
	if(model.mode == SIM) chain += 10;
	srand(chain+11);

	switch(model.mode){
	case SIM:                               // Simulates from the model
		{
			Simulate simu(model,output);
			simu.run();
		}
		break;
		
	case INF:                               // Performs inference on the model
		{
			switch(model.details.algorithm){
			case DA_MCMC:
				{
					MCMC mcmc(model,output);
					mcmc.run();
				}
				break;
				
			case ABC_ALG:
				{
					ABC abc(model,output);
					abc.run();
				}
				break;
				
			case ABC_SMC_ALG:
				{
					ABC_SMC abcsmc(model,output);
					abcsmc.run();
				}
				break;
				
			default: 
				emsg("This algorithm has not been implemented yet");
				return 0;
			}
		}
		break;
		
	case MODE_UNSET:
		emsg("The mode is unset");
		break;
	}

	output.updated_file(file);
	
	output.generate_files();
}

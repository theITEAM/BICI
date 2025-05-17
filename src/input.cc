/// Loads up model fron the input ./bici file

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath> 
 
using namespace std;

#include "input.hh"
#include "output.hh"
#include "utils.hh"

/// Initialises the model 
Input::Input(Model &model, string file, unsigned int seed, Mpi &mpi) : model(model), mpi(mpi)
{
	datadir = "";

	terminate = false;
	
	input_file = file;
	
	percentage_start(LOAD_PER);
	
	percentage(0,100);
	
	if(op()){
		ifstream fin(file);
		if(!fin) emsg("File '"+file+"' could not be loaded");
	
		do{
			string line;

			getline(fin,line);
		
			if(fin.eof()) break;
			
			remove_cr(line);
		
			lines_raw.push_back(line);
			
			line = remove_escape_char(line);
		
			lines.push_back(line);
		}while(true);
	}
	
#ifdef USE_MPI 
	mpi.bcast(lines);
	mpi.bcast(lines_raw);
#endif

	auto command_line = extract_command_line(lines); // Converts from text lines to command lines

	load_data_files(command_line);

	// Import happens in three stages:
	// (0) Determines if simulation, inference or post_sim
	// (1) Simulation or inference details
	// (2) Load species, classification and compartment information
	// (3) Load parameter information
	// (4) Adds inidividuals to the system
	// (5) Everything else

	for(auto loop = 1u; loop < 6; loop++){ 
		print_diag("loop="+to_string(loop));

		// Keeps track of the current species and classification 
		p_current = UNSET; cl_current = UNSET; 
	
		if(model.nspecies == 1){
			p_current = 0;
			if(model.species[p_current].ncla == 1) cl_current = 0;
		}
	
		for(const auto &line : command_line){
			line_num = line.line_num;
			auto cname = line.command;
			
			auto process = true;
			switch(loop){
			case 1:  // Gets simulation or inference details
				switch(cname){
				case SIMULATION: 
					if(model.mode != SIM) process = false; 
					break;
				
				case INFERENCE:
					if(model.mode != INF && model.mode != PPC) process = false; 
					break;
					
				case POST_SIM:
					if(model.mode != PPC) process = false; 
					break;
				
				case DATA_DIR: break;
				
				default: process = false; break;
				}
				break;
				
			case 2:   // In the first pass create species, classifications and compartments
				switch(cname){
				case SPECIES: case CLASS: case COMP: case COMP_ALL: break;
				default: process = false; break;
				}
				break;
				
			case 3:   // In the next pass the parameters 
				switch(cname){
				case PARAM: case DERIVED: break;
				default: process = false; break;
				}
				break;
			
			case 4:   // Adding individuals to the system
				switch(cname){
				case SPECIES: break;
				
				case PARAM_MULT:
					if(model.mode != PPC) process = false;
					break;
				
				case ADD_POP: case REMOVE_POP: case ADD_IND: 
					if(model.mode == SIM) process = false;
					break;
					
				case ADD_POP_SIM: case REMOVE_POP_SIM: case ADD_IND_SIM:
				case ADD_POP_POST_SIM: case REMOVE_POP_POST_SIM: case ADD_IND_POST_SIM:
					if(model.mode == INF) process = false;
					break;
					
				default: process = false; break;
				}
				break;
				
			case 5:   // In the last pass everything else
				switch(cname){
				case DATA_DIR: 
				case SIMULATION: case INFERENCE: case POST_SIM:
				case PARAM: case DERIVED: 
				case ADD_POP: case ADD_POP_SIM: case ADD_POP_POST_SIM:
				case REMOVE_POP: case REMOVE_POP_SIM:  case REMOVE_POP_POST_SIM:
				case ADD_IND: case ADD_IND_SIM: case ADD_IND_POST_SIM: 
				case PARAM_MULT:
				case COMP: case COMP_ALL: 
					process = false;
					break;
					
				case INIT_POP_SIM: case REMOVE_IND_SIM: case MOVE_IND_SIM:
					if(model.mode == INF) process = false;
					break;
					
				case INIT_POP: case REMOVE_IND: case MOVE_IND:
				case COMP_DATA: case TRANS_DATA: case TEST_DATA: case POP_DATA: case POP_TRANS_DATA:
					if(model.mode == SIM) process = false;
					break;
				
				default: break;
				}
			}
			
			if(model.mode == PPC){ // If PPC then do not load any data
				switch(cname){
				case INIT_POP_SIM: 
				case ADD_POP_SIM: case REMOVE_POP_SIM: case ADD_IND_SIM:
				case REMOVE_IND_SIM: case MOVE_IND_SIM:
				case COMP_DATA: case TRANS_DATA: case TEST_DATA:
				case POP_DATA: case POP_TRANS_DATA: case IND_EFFECT_DATA:
				case IND_GROUP_DATA: case GENETIC_DATA:
					process = false;
					break;
					
				default: break;
				}	
			}
			
			if(process == true){
				process_command(line,loop);
			
				if(terminate == true){
					output_error_messages(err_mess);
					return;
				}
			}
		}
		line_num = UNSET;
			
		switch(loop){
		case 1: calculate_timepoint(); break;
		case 2: check_comp_structure(); break; 		
		}
	}
	
	output_error_messages(err_mess);
	
	percentage(10,100);
	
	auto inf = false; if(model.mode == INF) inf = true;
	
	print_diag("loaded");
	
	set_contains_source_sink();        // Sets flag to determine if source or sink
	
	set_trans_tree();                  // Sets if there is a transition tree
	
	create_population_erlang();        // If a population-based model convert erlang to rate
	
	percentage(20,100);
	
	print_diag("loaded1");
	
	determine_branching();             // Determines which transitions branch from compartments
	
	print_diag("loaded2");
	
	bp_create_unset();                 // Substitutes "*" with the expression form branching prob
	
	print_diag("loaded3");
	
	markov_bp_convert();               // If all branches are rates then removes branching

	print_diag("loaded4");
	percentage(30,100);
	
	//population_bp_rate_combine();      // Combines branching probabilities with rates in population model
	
	print_diag("loaded5");
	
	check_initial_pop_error(true);     // Determines if initial population is specified correctly

	print_diag("loaded6");
	
	check_import_correct();            // Checks import has been successfully acheived

	print_diag("loaded7");
	percentage(40,100);
	
	global_comp_trans_init();          // Creates global compartments and transitions
	
	print_diag("loaded8");

	for(auto &sp : model.species){	
		sp.initialise_data();            // Extracts data structures from data sources
	
		for(const auto &wa : sp.warn) alert_line(wa.te,wa.line_num);    
	}
	
	percentage(50,100);

	print_diag("loaded9");
	
	model.set_hash_all_ind();           // Sets a hash table for all individuals
	
	add_genetic_data();                // Adds genetic data (potentially multi-species)
		
	source_rate_divide();              // Divides equations for source rates
	
	print_diag("h1");

	percentage(60,100);
	
	create_equations();                // Creates equation calculations
	
	percentage(70,100);
	
	print_diag("h1a");
	
	simplify_equations();              // Simplifies equations as much as possible

	percentage(80,100);
	
	print_diag("h1b");
		
	for(auto &sp : model.species){     // Classifies observed transitons as trans, source, or sink
		sp.set_ob_trans_ev(model.eqn);
	}

	print_diag("h1c");
	
	for(auto &eqn : model.eqn){        // Sets up reference (pop_ref, param_ref) in equations
		eqn.setup_references();
	}	
	
	check_reparam_time();              // Checks if reparameterised equations involve time
	
	print_diag("h1d");
	map_ind_effect();                  // Maps individual effects with groups

	set_param_use();                   // Sets which parameters are used in the model
	
	print_diag("h1e");
	
	set_param_parent_child();          // Sets parent child relationships for parameters
	
	print_diag("h1f");
	
	//exp_nm_convert();                  // If doesn't contain a population then EXP_RATE -> EXP_RATE_NM
	
	print_diag("h1g");
	
	setup_obs_trans();                 // Sets up obs_trans
	
	print_diag("h1h");
	
	setup_obs_trans_is_one();          // Sets up obs_trans is_one
	
	print_diag("h1i");
	
	create_nm_trans();                 // Creates a list of possible non-Markovian transitions

	print_diag("h2");
	
	source_equation_comp();            // Attaches compartments to source equations

	print_diag("h3");
	
	create_markov_eqn();               // Works out which equations are Markovian
	
	print_diag("h4");
	
	create_markov_comp_gl();           // Works out which markov_eqn applied to a compartment 
	
	print_diag("h5");
	
	create_markov_eqn_pop_ref();       // Works out which populations affect which markov eqns (and vice versa)
	
	print_diag("h6");

	create_pop_ref();                  // Create reference to pop in comp_gl
	
	print_diag("h8");
	
	create_param_vector();             // Creates an overall parameter vector
	
	print_diag("h9");
	
	create_spline();                   // Creates any splines in the model

	print_diag("h10");
	
	ind_fix_eff_group_trans_ref();     // References markov eqns and nm_trans ind_eff_group and fix_effect
	
	print_diag("h11");
	
	ind_fix_eff_pop_ref();             // References populations in ind_effect and fix_effect
	
	if(model.mode == INF || model.mode == PPC){
		create_trg_from_tr();      // Works out a convertion from tr to trg
	}
	
	print_diag("h12");
	
	set_tra_ie_fe_affect();            // Sets which ie affect transitions
	
	set_tr_leave_markov();             // Sets Markov transitions leaving compartment    
	
	print_diag("h13");
	
	set_multi_trans_cl();              // Sets if transitions in multiple classification
	
	print_diag("h14");
		
	set_tr_connected();                // Determines if global transitions are connected
	
	for(auto &sp : model.species){     // Determines if individual sampler needed
		if(inf) sp.set_ind_samp_needed(model.eqn);     
	}
	
	set_comp_period();                 // Determines if compartemnt begins or ends a transition period
	
	set_comp_terminal();               // Determines if a compartment is terminal
	
	set_sink_exist();                  // Determines if a sink exists
	
	set_eqn_ind_eff_exist();           // Sets if ind_eff exist in equations
	
	//set_pop_speedup();                 // Sets speedup for update_ind when population is changed
 	
	print_diag("h15");
	
	percentage(90,100);
	
	for(auto &eq : model.eqn){
		eq.calculate_linearise();        // Tries to linearise equations in terms of populations
	}
	
	if(inf) param_affect_likelihood(); // Works out how changes to parameters affect likelihoods
	
	print_diag("h16");
	
	check_param_used();                // Checks all defined parameters used in the model (and vice versa)
	
	print_diag("h18");
	
	for(auto &sp : model.species){
		sp.create_markov_tree();         // Creates sampler used to sample Markov events
		sp.init_pop_trans_ref();         // Initialises reference from global transition to obs         
		sp.init_pop_data_ref();          // Initialises reference from global compartment to obs         
		for(const auto &wa : sp.warn) alert_line(wa.te,wa.line_num);    
	}
	
	print_diag("h19");
	
	if(inf) create_island();                   // Works out how interconnected compartments are
	
	print_diag("h20");
	
	setup_trans_infection();                   // Sets up information about infection down transitions (for trans-tree)
	
	for(auto &sp : model.species){
		sp.create_ind_noobs();                   // Creates an "individual" with no observation
	}
	
	create_trig_ev_ref();                      // Generates trig_ev_ref (events on individual timelines)
	create_cl_trig_ev_ref();                   // Generates cl_trig_ev_ref (used in sampling)
	
	set_precalc_nm_rate();                     // Determines if nm_rate can be precalculated
	
	set_tr_enter();                            // Sets transitions entering compartment
 	               
	if(inf) set_joint_param_event();           // Sets any joint parameter and event proposals
	
	set_init_c_set();                          // Determines if initial compartment is set
	
	if(inf){ 
		set_local_ind_init();                    // Initialises local individual changes
		set_comp_global_convert();               // Sets comp_global_convert		
	}
	
	set_tra_global_time_vari();                // Sets the time varying parameter for transitions
	
	set_ind_variation();                       // Sets flag determining if individual variation
	
	set_cgl_begin_nm();                        // Determines if glob comp starts nm 
	
	set_cgl_tr_source_sink();                  // Sets trg reference from compartments
	
	set_omega_fl();                            // Sets a flag if param in in omega
	
	if(model.trans_tree){
		set_inf_cause();                         // Reference to see the population which causes infection 
	}
	
	if(model.mode == PPC && model.sample.size() == 0){
		alert_import("Posterior samples (from the commands 'inf-param' and 'inf-state') must be specified for 'post-sim' to be run.");
	}		

	if(model.mode == PPC) set_ppc_resample();
		
	check_markov_or_nm();                      // Checks that transition are either markovian or non-markovian
		
	check_nm_pop();                            // Checks no population in species for nm trans	
	
	set_seed(mpi.core,model.details,seed);     // Sets the psuedo random nunber generator 
	
	percentage(100,100);
	
	output_error_messages(err_mess,true);

	print_diag("Finish");
}


/// Converts from text lines to command lines
vector <CommandLine> Input::extract_command_line(vector <string> lines)
{
	vector <CommandLine> command_line;
	
	for(auto j = 0u; j < lines.size(); j++){
		line_num = j;
		
		auto trr = trim(lines[j]);
	
		auto flag = false;                    // Ignores line if empty or commented out
		if(trr.length() == 0) flag = true;
		if(trr.length() >= 2 && trr.substr(0,2) == "//") flag = true;
		if(trr.length() >= 1 && trr.substr(0,1) == "#") flag = true;
		
		if(flag == false){
			while(trr.length() > 3 && trr.substr(trr.length()-3,3) == "\"[["){
				j++;
				vector <string> file_lines;
				string file_name = "$FILE"+to_string(files.size())+"&";
				while(j < lines.size()){
					auto te = trim(lines[j]);
					if(te.length() >= 3 && te.substr(0,3) == "]]\""){
						trr += file_name+te;
						break;
					}
					file_lines.push_back(te);
					j++;
				}
				if(j == lines.size()) alert_import("Cannot find line with ']]\"' to specify end of file");
			
				FileStore file_st;
				file_st.name = "[["+file_name+"]]";
				file_st.lines = file_lines;
				file_st.sep = ',';
				files.push_back(file_st);
			}
			
			auto line = get_command_tags(trr,line_num);
			if(line.command != EMPTY){
				command_line.push_back(line);
			}				
		}
	}
	
	return command_line;
}
	
	
/// Loads up any data files
void Input::load_data_files(vector <CommandLine> &command_line)
{
	string data_dir = "";
	for(auto j = 0u; j < command_line.size(); j++){
		if(command_line[j].command == DATA_DIR){
			for(auto k = 0u; k < command_line[j].tags.size(); k++){
				const auto &tag = command_line[j].tags[k];
				if(tag.name == "folder") data_dir = tag.value;
			}
		}
	}
	
	convert_folder(data_dir);
	
	for(auto j = 0u; j < command_line.size(); j++){
		const auto &cl = command_line[j];
			
		for(auto k = 0u; k < cl.tags.size(); k++){
			const auto &tag = cl.tags[k];
			
			auto na = tag.name;
				
			if(na == "value" || na == "boundary" || na == "constant" || na == "reparam" || na == "prior-split" || na == "dist-split" || 
				na == "A" || na == "A-sparse" || na == "pedigree" || na == "X" || na == "file" || na == "text" || na == "ind-list" || na == "factor-weight"){
				auto file = tag.value;
				
				if(is_file(file)){
					auto k = 0u; while(k < files.size() && files[k].name != file) k++;
					if(k == files.size()){
						if(data_dir == "") alert_import("The 'data-dir' command must be set");
							
						auto full_name = data_dir+"/"+file; 
						ifstream fin(full_name);
						if(!fin) emsg_input("File '"+full_name+"' could not be loaded");
						
						vector <string> lines;
						
						string line;
						do{
							getline(fin,line);
							if(fin.eof()) break;
							remove_cr(line);
							lines.push_back(line);
						}while(true);
	
						FileStore file_st;
						file_st.name = file;
						file_st.lines = lines;
						file_st.sep = ',';
						auto end = file.substr(file.length()-4); if(end == ".tsv") file_st.sep = '\t';
						files.push_back(file_st);
					}
				}
			}
		}
	}
}


/// Adds text escape characters
string Input::remove_escape_char(string te)
{
	auto escape_char = get_escape_char();

	te = replace(te,"〈","<");
	te = replace(te,"〉",">");
	
	auto i = 0u;
	while(i < te.length()){
		string sq = "′";
		if(te.substr(i,sq.length()) == sq){   // Converts prime to single quote
			te = te.substr(0,i)+"'"+te.substr(i+sq.length());
		}
			
		if(te.substr(i,1) == "\\"){
			auto j = 0u; 
			while(j < escape_char.size() && te.substr(i,escape_char[j][0].length()) != escape_char[j][0]) j++;
			if(j < escape_char.size()){
				te = te.substr(0,i)+escape_char[j][1]+te.substr(i+escape_char[j][0].length());
			}
			else i++;
		}
		else i++;
	}
	
	return te;
}


/// Converts folder to linux format
void Input::convert_folder(string &data_dir) const
{
	if(begin_str(data_dir,"M:/Github")){
		data_dir = "/nfs/home/cpooley/Github"+data_dir.substr(9);
	}
	
/*
#ifdef WINDOWS
#else

	
	//data_dir = replace(data_dir,"\\","/");
#endif
*/
}


/// Gets all the command tags from a given line
CommandLine Input::get_command_tags(string trr, unsigned int line_num)
{
	auto num_quote = 0u;
	vector <unsigned int> quote_pos;
			
	for(auto i = 0u; i < trr.length(); i++){
		if(trr.substr(i,1) == "\""){ num_quote++; quote_pos.push_back(i);}
	}
	
	if(num_quote%2 != 0){ alert_import("Syntax error: Quotes do not match up."); return syntax_error();}
	
	for(auto i = 0u; i < num_quote; i += 2){
		if(quote_pos[i]+1 == quote_pos[i+1]){
			alert_import("Syntax error: No content within the quotation marks."); 
			return syntax_error();
		}
	}

	vector <Fragment> frag; 
	
	auto i = 0u; 
	auto quote = 0u;
	auto after_eq = false;
	
	do{
		auto ist = i; i++;
		
		while(i < trr.length() && !(trr.substr(i,1) == "=" && quote == 0) && 
				!(quote == 0 && trr.substr(i,1) == "\"")&& 
				!(quote == 1 && trr.substr(i,1) == "\"") && 
				!(quote == 0 && trr.substr(i,1) == " ") &&
				!(quote == 0 && after_eq == true)){
			i++;
		}

		after_eq = false;
		if(trr.substr(i,1) == "=") after_eq = true;
		
		auto te = trr.substr(ist,i-ist);
		while(te.substr(0,1) == "\n") te = te.substr(1);
		while(te.substr(te.length()-1,1) == "\n") te = te.substr(0,te.length()-1);

		te = trim(te); if(quote == 0) te = toLower(te);

		if(te != ""){
			auto pos = ist; while(pos < i && trr.substr(pos,1) == " ") pos++;
			Fragment fr; fr.text = trim(trr.substr(ist,i-ist)); fr.pos = pos; fr.pos_end = i; fr.quote = quote;
			frag.push_back(fr);
		}
		
		while(i < trr.length() && quote == 0 && trr.substr(i,1) == " ") i++;
		if(i < trr.length()){
			if(trr.substr(i,1) == "\""){ 
				quote = 1-quote; i++;
			}
		}
	}while(i < trr.length());
	
	if(frag.size() == 0){ alert_import("Does not contain any content"); return syntax_error();}
	
	if(frag[0].quote == 1){ alert_import("Should not start with quotes"); return syntax_error();}
	auto type = frag[0].text;
	
	Command com = EMPTY;
	if(type == "species") com = SPECIES;
	if(type == "classification" || type == "class") com = CLASS;
	if(type == "set") com = SET;
	if(type == "camera" || type == "view") com = CAMERA;
	if(type == "compartment" || type == "comp") com = COMP;
	if(type == "compartment-all" || type == "comp-all") com = COMP_ALL;
	if(type == "transition" || type == "trans") com = TRANS;
	if(type == "transition-all" || type == "trans-all") com = TRANS_ALL;
	if(type == "data-dir") com = DATA_DIR;
	if(type == "description" || type == "desc") com = DESC;
	if(type == "label") com = LABEL;
	if(type == "box") com = BOX;
	if(type == "parameter" || type == "param") com = PARAM;
	if(type == "derived" || type == "der") com = DERIVED;
	if(type == "init-pop") com = INIT_POP;
	if(type == "add-pop") com = ADD_POP;
	if(type == "remove-pop") com = REMOVE_POP;
	if(type == "add-ind") com = ADD_IND;
	if(type == "remove-ind") com = REMOVE_IND;
	if(type == "move-ind") com = MOVE_IND;
	if(type == "init-pop-sim") com = INIT_POP_SIM;
	if(type == "add-pop-sim") com = ADD_POP_SIM;
	if(type == "remove-pop-sim") com = REMOVE_POP_SIM;
	if(type == "add-ind-sim") com = ADD_IND_SIM;
	if(type == "remove-ind-sim") com = REMOVE_IND_SIM;
	if(type == "move-ind-sim") com = MOVE_IND_SIM;
	if(type == "add-pop-post-sim") com = ADD_POP_POST_SIM;
	if(type == "remove-pop-post-sim") com = REMOVE_POP_POST_SIM;
	if(type == "add-ind-post-sim") com = ADD_IND_POST_SIM;
	if(type == "remove-ind-post-sim") com = REMOVE_IND_POST_SIM;
	if(type == "move-ind-post-sim") com = MOVE_IND_POST_SIM;
	if(type == "comp-data") com = COMP_DATA;
	if(type == "trans-data") com = TRANS_DATA;
	if(type == "test-data") com = TEST_DATA;
	if(type == "pop-data") com = POP_DATA;
	if(type == "pop-trans-data") com = POP_TRANS_DATA;
	if(type == "ind-effect-data") com = IND_EFFECT_DATA;
	if(type == "ind-group-data") com = IND_GROUP_DATA;
	if(type == "genetic-data") com = GENETIC_DATA;
	if(type == "simulation" || type == "sim") com = SIMULATION;
	if(type == "inference" || type == "inf") com = INFERENCE;
	if(type == "posterior-simulation" || type == "post-sim") com = POST_SIM;
	if(type == "ind-effect") com = IND_EFFECT;
	if(type == "fixed-effect") com = FIXED_EFFECT;
	if(type == "sim-param") com = SIM_PARAM;
	if(type == "sim-state") com = SIM_STATE;
	if(type == "inf-param") com = INF_PARAM;
	if(type == "inf-generation") com = INF_GEN;
	if(type == "inf-state") com = INF_STATE;
	if(type == "post-sim-param") com = POST_SIM_PARAM;
	if(type == "post-sim-state") com = POST_SIM_STATE;
	if(type == "inf-diagnostics") com = INF_DIAGNOSTICS;
	if(type == "map") com = MAP;
	if(type == "post-sim" || type == "post-simulation" ) com = POST_SIM;
	
	if(type == "param-mult") com = PARAM_MULT;
	
	if(type == "trans-diag") com = TRANS_DIAG;
	
	if(com == EMPTY){ alert_import("Command '"+type+"' not recognised."); return syntax_error();}
	
	auto num = double(frag.size()-1)/3;
	auto numi = (unsigned int)(num);

	for(auto n = 0u; n < num; n++){
		auto ii = 1+n*3;

		if(frag[ii].text == "="){
			alert_import("An equal sign '=' is misplaced.");
			return syntax_error();
		}
			
		if(ii+2 >= frag.size()){
			if(ii+1 < frag.size() && frag[ii+1].text == "="){
				alert_import("The property '"+frag[ii].text+"' is unset");
			}
			else{
				alert_import("The text '"+frag[ii].text+"' cannot be understood");
			}
			return syntax_error();
		}
		
		if(frag[ii+1].text != "="){
			alert_import("The property '"+frag[ii].text+"' is missing an equals sign");
			return syntax_error();
		}
		
		if(ii+2 < frag.size() && frag[ii+2].text == "="){
			alert_import("The property '"+frag[ii].text+"' cannot be followed by '=='");
			return syntax_error();
		}
		
		if(ii+3 < frag.size() && frag[ii+3].text == "="){
			alert_import("The property '"+frag[ii].text+"' is unset");
			return syntax_error();
		}
	}

	if(num != numi){ alert_import("Syntax error"); return syntax_error();}
	
	vector <Tag> tags;

	for(auto n = 0u; n < num; n++){
		if(frag[1+3*n+0].quote != 0){ alert_import("Syntax error type 2"); return syntax_error();}
		if(frag[1+3*n+1].text != "="){ alert_import("Syntax error type 3"); return syntax_error();}
		if(frag[1+3*n+0].text == ""){ alert_import("Syntax error type 4"); return syntax_error();}
		if(frag[1+3*n+2].text == ""){ alert_import("Property "+frag[1+3*n+0].text+" must have content"); return syntax_error();}
		
		const auto &fr = frag[1+3*n+0];
		
		Tag tag; 
		tag.name = fr.text; tag.pos = fr.pos; tag.pos_end = fr.pos_end; 
		tag.value = frag[1+3*n+2].text; tag.done = 0;
		tags.push_back(tag);
	}
	
	if(tags.size() > 0){
		for(auto n = 0u; n < tags.size()-1; n++){
			for(auto nn = n+1; nn < tags.size(); nn++){
				if(tags[n].name == tags[nn].name){
					alert_import("The tag '"+tags[n].name+"' is set more than once");
					return syntax_error();
				}
			}
		}
	}
		
	CommandLine cline; 
	cline.command = com; cline.command_name = type; cline.type_pos = frag[0].pos; cline.tags = tags;
	cline.line_num = line_num;
	
	return cline;
}


/// Returns a syntax error for a line
CommandLine Input::syntax_error() const
{
	CommandLine cline; cline.command = EMPTY; return cline;
}


/// Error message for imported file 
void Input::alert(string st)                               
{
	add_error_mess(UNSET,st,ERROR_FATAL);
	
	if(fatal_error() == true && error_mess.size() >= ERR_MSG_MAX){
		output_error_messages("Total error limit exceeded");
	}
}


/// Adds a new error message to the list
void Input::add_error_mess(unsigned int line_num, string st, ErrorType type)
{
	if(st.length() > 0 && st.substr(st.length()-1,1) != ".") st += ".";
	ErrorMess em;
	em.line_num = line_num;
	em.error = st;
	em.type = type;
	
	// Makes sure message is not repeated
	auto i = 0u; 
	while(i < error_mess.size() && !(error_mess[i].type == type && 
	  (error_mess[i].error == st || (error_mess[i].line_num == line_num && line_num != UNSET)))){
		i++;
	}
	
	if(i == error_mess.size()) error_mess.push_back(em);	
}

 
/// Error message for imported file 
void Input::alert_import(string st, bool fatal)   
{
	add_error_mess(line_num,st,ERROR_FATAL);
	
	if(fatal){
		output_error_messages(err_mess); 
	}
	
	if(fatal_error() == true && error_mess.size() >= ERR_MSG_MAX){
		output_error_messages("Total error limit exceeded");
	}
}


/// Error message for imported file (for specific line)
void Input::alert_line(string st, unsigned int line)                               
{
	if(st.length() > 0 && st.substr(st.length()-1,1) != ".") st += ".";
	
	auto i = 0u; while(i < error_mess.size() && error_mess[i].line_num != line) i++;
	if(i == error_mess.size() || line == UNSET){
		add_error_mess(line,st,ERROR_FATAL);
	}

	if(fatal_error() == true && error_mess.size() >= ERR_MSG_MAX){
		output_error_messages("Total error limit exceeded");
	}
}


/// Alerts that an equation is incorectly specified
void Input::alert_equation(const EquationInfo &eqi, const string &warn)
{
	alert_line("For equation '"+eqi.te_raw+"': "+warn,eqi.line_num); 
}


/// Error message for imported file (for specific line)
void Input::alert_warning(string st)                               
{
	add_error_mess(line_num,st,ERROR_WARNING);
}


/// Determines if there has been a fatal error
bool Input::fatal_error() const 
{
	for(auto i = 0u; i < error_mess.size(); i++){
		if(error_mess[i].type == ERROR_FATAL) return true;
	}
	return false;
}


/// Outputs all the error messages
void Input::output_error_messages(string te, bool end) const 
{
	if(op() && (fatal_error() || end)){
		if(!com_op && error_mess.size() > 0) cout << endl;
		for(const auto &em : error_mess){
			if(!com_op && em.line_num < lines.size()){ 	
				cout << "\033[32m";
				cout <<  "Line " << em.line_num+1 << ": ";
				cout << "\033[0m";
				cout << add_escape_char(lines[em.line_num]) << endl;
			}

			switch(em.type){
			case ERROR_FATAL: display_error(em.error);	break;
			case ERROR_WARNING: display_warning(em.error); break;
			}
			cout << endl;
		}
	}
	
	if(fatal_error() == true){
		if(op()){
			if(te == err_mess && error_mess.size() == 1) te = err_mess_sing;
			cout << te << endl << endl;
		}
		end_code();
	}
}


/// Processes a given command
void Input::process_command(const CommandLine &cline, unsigned int loop)
{
	cline_store = cline;
	all_row = UNSET;
	
	auto nem = error_mess.size();
	
	auto cname = cline.command;
	
	switch(cname){
	case SPECIES: species_command(loop); break;
	case CLASS: classification_command(loop); break;
	case SET: set_command(); break;
	case CAMERA: camera_command(); break;
	case COMP: compartment_command(); break;
	case COMP_ALL: compartment_all_command(); break;
	case TRANS: transition_command(); break;
	case TRANS_ALL: transition_all_command(); break;
	case DATA_DIR: datadir_command(); break;
	case DESC: description_command(); break;
	case LABEL: label_command(); break;
	case BOX: box_command(); break;
	case PARAM: param_command(); break;
	case DERIVED: derived_command(); break;
	case SIMULATION: simulation_command(); break;
	case INFERENCE: inference_command(); break;
	case POST_SIM: post_sim_command(); break;
	
	case IND_EFFECT: ind_effect_command(); break;
	case FIXED_EFFECT: fixed_effect_command(); break;
	case PARAM_MULT: param_mult_command(); break;
	
	case INIT_POP_SIM:
	case ADD_POP_SIM: case REMOVE_POP_SIM: 
	case ADD_IND_SIM: case REMOVE_IND_SIM: case MOVE_IND_SIM: 
		if(model.mode == SIM) import_data_table_command(cname);
		else alert_warning("Line ignored because not need for simulation");
		break;
	
	case ADD_POP_POST_SIM: case REMOVE_POP_POST_SIM: 
	case ADD_IND_POST_SIM: case REMOVE_IND_POST_SIM: case MOVE_IND_POST_SIM: 
		if(model.mode == PPC) import_data_table_command(cname);
		else alert_warning("Line ignored because not need for simulation");
		break;
		
	case INIT_POP: 
	case ADD_POP: case REMOVE_POP:
	case ADD_IND: case REMOVE_IND: case MOVE_IND: 
	case COMP_DATA: case TRANS_DATA:
	case TEST_DATA: case POP_DATA: case POP_TRANS_DATA: 
	case GENETIC_DATA:
		if(model.mode != SIM) import_data_table_command(cname);
		else{
			alert_warning("Line ignored because not need for inference");
			for(auto &ta : cline_store.tags) ta.done = 1;
		}
		break;
	
	case IND_EFFECT_DATA:
		get_tag_value("name"); 
		get_tag_value("file"); 
		break;
		
	case IND_GROUP_DATA:
		get_tag_value("name"); 
		get_tag_value("file"); 
		break;
	
	case SIM_PARAM: case SIM_STATE: dummy_file_command(); break;
	case POST_SIM_PARAM: case POST_SIM_STATE: dummy_file_command(); break;
	case INF_PARAM: dummy_file_command(); break;
	case INF_DIAGNOSTICS: dummy_file_command(); break;
	case INF_GEN: dummy_file_command(); break;
	
	case INF_STATE: 
		if(model.mode == PPC) inf_state_command();
		else dummy_file_command();
		break;
		
	case MAP: map_command(); break;
		
	case TRANS_DIAG: 
		dummy_file_command();
		break;
		
	default: alert_import("Command not recognised"); return;
	}
	
	
	// Checks if all tags are used
	auto error = false;
	for(auto e = nem; e <  error_mess.size(); e++){
		if(error_mess[e].type == ERROR_FATAL) error = true;
	}

	if(error == false){
		const auto &tags = cline_store.tags;
		for(auto n = 0u; n < tags.size(); n++){
			if(tags[n].done != 1){ 
				alert_warning("Tag '"+tags[n].name+"' not used"); return;
			}
		}
	}
}


/// Gets the value of a tag
string Input::get_tag_value(string st)                                   
{
	auto &tags = cline_store.tags;

	tag_find = st;
	for(auto i = 0u; i < tags.size(); i++){
		if(tags[i].name == st){ tags[i].done = 1; return tags[i].value;}
	}
	return "";
}


/// Gets the value of a tag from a specified list
string Input::get_tag_val(string st, vector <Tag> &tags) 
{
	tag_find = st;
	for(auto i = 0u; i < tags.size(); i++){
		if(tags[i].name == st){ tags[i].done = 1; return tags[i].value;}
	}
	return "";
}


/// Error massage if a tag cannot be found
void Input::cannot_find_tag(bool fatal)                                   
{
	string te = "Cannot find the '"+tag_find+"' tag for '"+cline_store.command_name+"'";
	if(all_row != UNSET) te += "( line "+to_string(all_row+2)+" in file)";
	
	alert_import(te,fatal);
}


/// Checks is a particular value is one of a list of possible options
unsigned int Input::option_error(string na, string te, const vector <string> &pos, const vector <unsigned int> &conv)
{
	if(pos.size() != conv.size()){ 
		alert_import("Options do not match up"); return UNSET;
	}
	
	auto k = find_in(pos,te);
	if(k != UNSET) return conv[k];
	
	auto st = "'"+na+"' has a value '"+te+"' but should be chosen from one of the following options: ";
	for(auto i = 0u; i < pos.size(); i++){
		if(i != 0){
			if(i+1 == pos.size()) st += " or ";
			else st += ", ";
		}
		st += "'"+pos[i]+"'";
	}
	
	alert_import(st);
	return UNSET;
}


/// Checks if lat and lng are correctly specified
bool Input::check_latlng_error(double lat, double lng, bool fatal)
{
	if(lat > 90 || lat < -90){ 
		alert_import("For 'lat' the value '"+tstr(lat)+"' must be in the range -90 to 90",fatal); 
		return true;
	}
	
	if(lng > 180 || lng < -180){ 
		alert_import("For 'lng' the value '"+tstr(lng)+"' must be in the range -180 to 180",fatal); 
		return true;
	}

	return false;
}


/// Checks if a string is a number
bool Input::is_number(string num, string tag, bool fatal)
{
	auto x = number(num);
	if(x == UNSET){
		alert_import("For '"+tag+"' the value '"+num+"' must be a number",fatal); 
		return false;
	}
	return true;
}


/// Checks if a string is a positive number
bool Input::is_positive(string num, string tag, bool fatal)
{
	auto val = number(num);
	if(val == UNSET || val <= 0){
		alert_import("For '"+tag+"' the value '"+num+"' must be a positive number",fatal); 
		return false;
	}
	return true;
}


/// Checks if a string is between zero and one
bool Input::is_zeroone(string num, string tag)
{
	auto val = number(num);
	if(val == UNSET || val < 0 || val > 1){
		alert_import("For '"+tag+"' thel value '"+num+"' must be between 0 and 1"); 
		return false;
	}
	return true;
}


/// Maps from individual effects on species to those stored in IEgroups
void Input::map_ind_effect()
{
	for(auto &sp : model.species){
		// Works out parameter references for omega
		for(auto i = 0u; i < sp.ind_eff_group.size(); i++){
			auto &ieg = sp.ind_eff_group[i];
			
			auto N = ieg.list.size();
			
			ieg.omega.resize(N); 	for(auto j = 0u; j < N; j++) ieg.omega[j].resize(N);
			
			for(auto j = 0u; j < N; j++){
				for(auto jj = j; jj < N; jj++){
					auto ie1 = ieg.list[j].name, ie2 = ieg.list[jj].name;
					auto sup = "^"+ie1+","+ie2;
					
					string name;
					if(j == jj) name = "Ω"+sup;
					else name = "ω"+sup;
					
					auto th = 0u; while(th < model.param.size() && model.param[th].name != name) th++;
					
					if(th == model.param.size()){
						alert_line("Parameter '"+name+"' must be specified for these individual effects",ieg.line_num);
					}
					else{
						ieg.omega[j][jj] = th;
						ieg.omega[jj][j] = th;
						model.param[th].used = true;
					}
				}
			}		
		}
		
		// Links from ind_effect to ind_effect_group
		for(auto i = 0u; i < sp.ind_effect.size(); i++){
			auto &ie = sp.ind_effect[i];
			auto name = ie.name;
				
			auto flag = false;
			for(auto j = 0u; j < sp.ind_eff_group.size(); j++){
				auto &ieg = sp.ind_eff_group[j];
				for(auto k = 0u; k < ieg.list.size(); k++){
					if(ieg.list[k].name == name){
						ieg.list[k].index = i;
						ie.index = j;
						ie.num = k;
						flag = true;
						break;
					}
				}
				if(flag == true) break;
			}
			
			if(flag == false){
				alert_line("Individual effect '"+name+"' is not defined using the 'ind-effect' command.",ie.line_num);
			}
		}
		
		for(auto j = 0u; j < sp.ind_eff_group.size(); j++){
			const auto &ieg = sp.ind_eff_group[j];
			for(auto k = 0u; k < ieg.list.size(); k++){
				const auto &li = ieg.list[k];
				if(li.index == UNSET){
					alert_line("Individual effect '"+li.name+"' not used in the model.",ieg.line_num);
				}
			}
		}
	}
}


/// Checks if a parameter has been used or not
void Input::check_param_used()
{
	for(const auto &par : model.param){
		if(par.used == false){
			if(model.mode != PPC){
				alert_line("Parameter '"+par.full_name+"' is not used in the model",par.line_num);
			}
		}
	}
}

	
/// Orders parameters in a vector such that only dependent on an params with a smaller index
void Input::create_param_vector()
{
	auto N = model.param.size();
	
	// This has been moved to "use" in parameter
	/*
	// Finds which parameters are actually used in the equations
	vector < vector <bool> > eqn_used;
	eqn_used.resize(N);
	for(auto th = 0u; th < N; th++){
		const auto &par = model.param[th];
		eqn_used[th].resize(par.N,false);
	}
	
	for(auto &eq : model.eqn){
		for(auto &ca : eq.calc){
			for(auto &it : ca.item){
				if(it.type == PARAMETER){
					if(it.num >= model.param.size()) emsg_input("Out of range1");
					if(it.index >= model.param[it.num].N) emsg_input("Out of range2");
					eqn_used[it.num][it.index] = true;
				}
			}
		}
	
		if(eq.ans.type == PARAMETER){
			if(eq.ans.num == UNSET) emsg_input("done");
			eqn_used[eq.ans.num][eq.ans.index] = true;
		}
	}
	
	for(const auto &sp : model.species){
		for(const auto &ieg : sp.ind_eff_group){
			auto N = ieg.list.size(); 
			for(auto j = 0u; j < N; j++){
				for(auto i = 0u; i < N; i++){
					eqn_used[ieg.omega[j][i]][0] = true;
				}
			}
		}
		
		for(const auto &fe : sp.fix_effect){
			eqn_used[fe.th][0] = true;
		}
	}
	
	if(false){
		for(auto th = 0u; th < N; th++){
			const auto &par = model.param[th];
			cout << par.name << endl;
			auto imax = eqn_used[th].size(); if(imax > 100) imax = 100;
			for(auto i = 0u; i < imax; i++) cout << eqn_used[th][i] << " ";
			cout << endl;
		}
		emsg_input("Shows parameters used");
	}
	*/
	
	vector < vector <unsigned int> > needed;
	
	needed.resize(N);

	for(auto th = 0u; th < N; th++){		
		const auto &par = model.param[th];
		for(auto i = 0u; i < par.parent.size(); i++){
			for(auto j = 0u; j < par.parent[i].size(); j++){
				auto th2 = par.parent[i][j].th;
			
				if(th2 == th){
					alert_import("In parameter '"+par.name+"' error with one element dependent on another element");
					return;
				}
				else{
					add_to_vec(needed[th],th2);
				}
			}				
		}
	}	
	
	vector <bool> defined(N,false);
	
	vector <unsigned int> list;
	
	bool flag;
	do{
		flag = false;
		
		for(auto th = 0u; th < N; th++){
			if(defined[th] == false){
				auto j = 0u; while(j < needed[th].size() && defined[needed[th][j]] == true) j++;
			
				if(j == needed[th].size()){
					defined[th] = true; 	
					list.push_back(th);
					flag = true;
				}
			}
		}
	}while(flag == true);
	
	if(list.size() != N){
		alert_import("Could not order parameters (e.g. cyclic dependencies exist such as A dependent on B and B dependent on A)");
	}
	
	for(auto i = 0u; i < N; i++){
		auto th = list[i];
		auto &par = model.param[th];
		
		par.param_vec_ref.resize(par.N,UNSET);
		
		if(par.variety != CONST_PARAM || par.factor){	
			for(auto j = 0u; j < par.N; j++){
				if(removeparamvec_speedup == false || 
				   //eqn_used[th][j] == true || par.spline_info.on == true){
				   par.use[j] == true || par.spline_info.on == true){
					par.param_vec_ref[j] = model.param_vec.size();
					
					ParamVecEle pr; 
					pr.name = get_param_name_with_dep(par,par.dep,j);
					pr.th = th; 
					pr.index = j;
					pr.prior = par.prior[j];
					pr.variety = par.variety;
					pr.ppc_resample = false;
				
					pr.prop_pos = false;
					if(pr.variety == PRIOR_PARAM || pr.variety == DIST_PARAM){
						if(pr.prior.type != FIX_PR){
							pr.prop_pos = true;
						}
					}
					
					model.param_vec.push_back(pr);
				}
			}		 
		}
	}
	
	// Converts parameter references in equation to param_vec
	for(auto &eq : model.eqn){
		for(auto &ca : eq.calc){
			for(auto &it : ca.item){
				if(it.type == PARAMETER){
					if(it.num >= model.param.size()) emsg_input("Out of range1");
					if(it.index >= model.param[it.num].param_vec_ref.size()) emsg_input("Out of range2");
					it.num = model.param[it.num].param_vec_ref[it.index];
					it.index = UNSET;
				}
			}
		}
	
		if(eq.ans.type == PARAMETER){
			eq.ans.num = model.param[eq.ans.num].param_vec_ref[eq.ans.index];
			eq.ans.index = UNSET;
		}
	}
	model.nparam_vec = model.param_vec.size();
	
	if(false){
		auto imax = model.nparam_vec;
		if(imax > 100) imax = 100;
		for(auto i = 0u; i < imax; i++){
			cout << model.param_vec[i].name << " vec" << endl;
		}
		
		emsg_input("param vec");
	}
}


/// Creates a reference to pop from comp_gl
void Input::create_pop_ref() 
{
	for(auto po = 0u; po < model.pop.size(); po++){
		const auto &pop = model.pop[po];
		for(auto i = 0u; i < pop.term.size(); i++){
			const auto &te = pop.term[i];
			PopRef pr; pr.po = po; pr.index = i;
			model.species[pop.sp_p].comp_gl[te.c].pop_ref.push_back(pr);
			model.species[pop.sp_p].comp_gl[te.c].pop_ref_simp.push_back(po);
		}
	}
}


/// Simplifies equations as much as possible
void Input::simplify_equations()
{
	if(simplify_eqn == false) return; 
	
	bool flag_global;
	do{
		flag_global = false;
			
		// If a reparameterised parameter is used and it is constant than substitutes 
		for(auto &eq : model.eqn){
			auto flag = false;
			for(auto &ca : eq.calc){
				for(auto &it : ca.item){
					if(it.type == PARAMETER){
						
						const auto &par = model.param[it.num];
						if(par.variety == REPARAM_PARAM){
							const auto &eqn = model.eqn[par.value[it.index].eq_ref];
							
							if(eqn.calc.size() == 0){
								if(eqn.ans.type == NUMERIC){
									it.type = NUMERIC;
									it.constant = eqn.ans.constant;
									flag = true;
								}
							}
						}
					}
				}
			}
			
			if(flag == true){
				eq.simplify();
				flag_global = true;
			}
		}
	}while(flag_global == true);
}

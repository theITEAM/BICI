/// Loads up model and data fron input file

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
 
//#include "math.h"

using namespace std;

#include "input.hh"
#include "output.hh"
#include "utils.hh"

/// Initialises the model 
Input::Input(Model &model, string file) : model(model)
{
	//generate_data();
	terminate = false;
	
	input_file = file;
	
	ifstream fin(file);
	if(!fin) emsg("File '"+file+"' could not be loaded");
	
	do{
		string line;
		getline(fin,line);
		
		if(fin.eof()) break;
		
		remove_cr(line);
	
		lines_raw.push_back(line);
		
		line = add_escape_char(line);
	
		lines.push_back(line);
	}while(true);

	auto command_line = extract_command_line(lines); // Converts from text lines to command lines
	
	// Import happens in three stages:
	// (0) Determines if simulation or inference
	// (1) Load species, classification and compartment information
	// (2) Load parameter information
	// (3) Everything else

	for(auto loop = 0u; loop < 4; loop++){ 
		//cout << loop << "LOOP\n";

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
			case 0:  // The initial pass determines if simulation or inference is being performed
				switch(cname){
				case SIMULATION: case INFERENCE: case DATA_DIR: break;
				default: process = false; break;
				}
				break;
				
			case 1:   // In the first pass create species, classifications and compartments
				switch(cname){
				case SPECIES: case CLASS: case COMP: break;
				default: process = false; break;
				}
				break;
				
			case 2:   // In the next pass the parameters
				switch(cname){
				case PARAM: case DERIVED: break;
				default: process = false; break;
				}
				break;
				
			case 3:   // In the last pass everything else
				switch(cname){
				case SIMULATION: case INFERENCE: case DATA_DIR: case PARAM: case DERIVED: case COMP:
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
		
		if(loop == 1){
			check_comp_structure(); 	
		}
	}

	//raw_process();
	
	create_population_erlang();        // If a population-based model convert elang to rate
	
	determine_branching();             // Determines which transitions branch from compartments

	population_bp_rate_combine();      // Combines branching probabilities with rates in population model
	
	check_initial_pop_error(true);     // Determines if initial population is specified correctly

	if(datadir == "") alert_import("'datadir' must be set");
		
	check_import_correct();            // Checks import has been successfully acheived

	global_comp_trans_init();          // Creates global compartments and transitions

	create_equations();                // Creates equation calculations

	source_equation_comp();            // Attaches compartments to source equations

	create_markov_eqn();               // Works out which equations are Markovian

	create_markov_comp_gl();           // Works out which markov_eqn applie to a compartment 
	
	create_markov_eqn_pop_ref();       // Works out which populations affect which markov eqns (and vice versa)

	create_nm_trans();                 // Creates a list of possible non-Markovian transitions
	
	create_pop_ref();                  // Create reference to pop in comp_gl
	
	create_param_vector();             // Creates an overall parameter vector
	
	calculate_timepoint();             // Caclualte the timepoints used to solve the equations
	
	create_spline();                   // Creates any splines in the model

	branch_param_group();              // Determines groups of parameters which determine branching
	
	map_ind_effect();                  // Maps individual effects with groups
	
	ind_eff_group_trans_ref();         // References markov eqns and nm_trans ind_eff_group
	
	ind_eff_pop_ref();                 // References populations in ind_effect
	
	if(model.mode == INF){ 	
		param_affect_likelihood();       // Works out how changes to parameters affect likelihoods
	}
	
	if(model.mode == INF){
		create_island();                 // Works out how interconnected compartments are
	}
	
	create_trg_from_tr();              // Works out a convertion from tr to trg
	
	check_param_used();                // Checks all defined parameters used in the model (and vice versa)
	
	for(auto &sp : model.species){
		sp.create_markov_tree();         // Creates sampler used to sample Markov events
		sp.initialise_data();            // Extracts data structures from data sources
		for(const auto &wa : sp.warn) alert_line(wa.te,wa.line_num);    
	}
	
	output_error_messages(err_mess);
}


/// Converts from text lines to command lines
vector <CommandLine> Input::extract_command_line(vector <string> lines)
{
	vector <CommandLine> command_line;
	
	for(line_num = 0; line_num < lines.size(); line_num++){
		auto trr = trim(lines[line_num]);
	
		auto flag = false;                    // Ignores line if empty or commented out
		if(trr.length() == 0) flag = true;
		if(trr.length() >= 2 && trr.substr(0,2) == "//") flag = true;
		if(trr.length() >= 1 && trr.substr(0,1) == "#") flag = true;
		
		if(flag == false){
			auto line = get_command_tags(trr,line_num);
			if(line.command != EMPTY){
				command_line.push_back(line);
			}				
		}
	}
	
	return command_line;
}
	

/// Adds text escape characters
string Input::add_escape_char(string te)
{
	vector< vector <string> > escape_char;
	
	escape_char.push_back({"\\alpha","α"});
	escape_char.push_back({"\\beta","β"});
	escape_char.push_back({"\\gamma","γ"});
	escape_char.push_back({"\\delta","δ"});
	escape_char.push_back({"\\epsilon","ε"});
	escape_char.push_back({"\\zeta","ζ"});
	escape_char.push_back({"\\eta ","η"});
	escape_char.push_back({"\\theta","θ"});
	escape_char.push_back({"\\iota","ι"});
	escape_char.push_back({"\\kappa","κ"});
	escape_char.push_back({"\\lambda ","λ"});
	escape_char.push_back({"\\mu","μ"});
	escape_char.push_back({"\\nu","ν"});
	escape_char.push_back({"\\xi","ξ"});
	escape_char.push_back({"\\pi","π"});
	escape_char.push_back({"\\rho","ρ"});
	escape_char.push_back({"\\sigma","σ"});
	escape_char.push_back({"\\tau","τ"});
	escape_char.push_back({"\\upsilon","υ"});
	escape_char.push_back({"\\phi","φ"});
	escape_char.push_back({"\\chi","χ"});
	escape_char.push_back({"\\psi","ψ"});
	escape_char.push_back({"\\omega","ω"});	
	escape_char.push_back({"\\sum","Σ"});		
	
	te = replace(te,"〈","<");
	te = replace(te,"〉",">");
	//te = replace(te,"<","〈");
	//te = replace(te,">","〉");
	
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


/// Gets all the command tags from a given line
CommandLine Input::get_command_tags(string trr, unsigned int line_num)
{
	auto num_quote = 0u;
	vector <unsigned int> quote_pos;
	
	for(auto i = 0u; i < trr.length(); i++){
		if(trr.substr(i,1) == "\""){ num_quote++; quote_pos.push_back(i);}
	}
	
	if(num_quote%2 != 0){ alert_import("Syntax error: Quotes do not match up."); return syntax_error();}
	
	for(auto i = 0u; i < num_quote; i+= 2){
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
	if(type == "camera") com = CAMERA;
	if(type == "compartment" || type == "comp") com = COMP;
	if(type == "transition" || type == "trans") com = TRANS;
	if(type == "source") com = SOURCE_CD;
	if(type == "sink") com = SINK_CD;
	if(type == "data-dir") com = DATA_DIR;
	if(type == "description" || type == "desc") com = DESC;
	if(type == "label") com = LABEL;
	if(type == "box") com = BOX;
	if(type == "parameter" || type == "param") com = PARAM;
	if(type == "derived" || type == "der") com = DERIVED;
	if(type == "init-pop") com = INIT_POP;
	if(type == "add-ind") com = ADD_IND;
	if(type == "remove-ind") com = REMOVE_IND;
	if(type == "move-ind") com = MOVE_IND;
	if(type == "init-pop-sim") com = INIT_POP_SIM;
	if(type == "add-ind-sim") com = ADD_IND_SIM;
	if(type == "remove-ind-sim") com = REMOVE_IND_SIM;
	if(type == "move-ind-sim") com = MOVE_IND_SIM;
	if(type == "init-pop-prior") com = INIT_POP_PRIOR;
	if(type == "comp-data") com = COMP_DATA;
	if(type == "trans-data") com = TRANS_DATA;
	if(type == "source-data") com = SOURCE_DATA;
	if(type == "sink-data") com = SINK_DATA;
	if(type == "test-data") com = TEST_DATA;
	if(type == "pop-data") com = POP_DATA;
	if(type == "pop-trans-data") com = POP_TRANS_DATA;
	if(type == "set-traps-data") com = SET_TRAPS_DATA;
	if(type == "ind-trapped-data") com = IND_TRAPPED_DATA;
	if(type == "genetic-data") com = GENETIC_DATA;
	if(type == "simulation" || type == "sim") com = SIMULATION;
	if(type == "inference" || type == "inf") com = INFERENCE;
	if(type == "ind-effect") com = IND_EFFECT;
	if(type == "fixed-effect") com = FIXED_EFFECT;
	if(type == "sim-param") com = SIM_PARAM;
	if(type == "sim-state") com = SIM_STATE;
	if(type == "inf-param") com = INF_PARAM;
	if(type == "inf-state") com = INF_STATE;
	
	if(com == EMPTY){ alert_import("Command '"+type+"' not recognised."); return syntax_error();}
	
	auto num = double(frag.size()-1)/3;

	if(num != int(num)){ alert_import("Syntax error"); return syntax_error();}
	
	vector <Tag> tags;

	for(auto n = 0u; n < num; n++){
		if(frag[1+3*n+0].quote != 0){ alert_import("Syntax error type 2"); return syntax_error();}
		if(frag[1+3*n+1].text != "="){ alert_import("Syntax error type 3"); return syntax_error();}
		if(frag[1+3*n+0].text == ""){ alert_import("Syntax error type 4"); return syntax_error();}
		if(frag[1+3*n+2].text == ""){ alert_import("Property "+frag[1+3*n+0].text+" must have content"); return syntax_error();}
		
		const auto &fr = frag[1+3*n+0];
		
		Tag tag; tag.name = fr.text; tag.pos = fr.pos; tag.pos_end = fr.pos_end; tag.value = frag[1+3*n+2].text; tag.done = 0;
		tags.push_back(tag);
	}
	
	for(auto n = 0u; n < tags.size()-1; n++){
		for(auto nn = n+1; nn < tags.size(); nn++){
			if(tags[n].name == tags[nn].name){
				alert_import("The tag '"+tags[n].name+"' is set more than once");
				return syntax_error();
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
	if(st.length() > 0 && st.substr(st.length()-1,1) != ".") st += ".";
	
	ErrorMess em;
	em.line_num = UNSET;
	em.error = st;
	em.type = ERROR_FATAL;
	
	error_mess.push_back(em);	
	if(fatal_error() == true && error_mess.size() >= ERR_MSG_MAX) output_error_messages("Total error limit exceeded");
}

 
/// Error message for imported file 
void Input::alert_import(string st)                               
{
	if(st.length() > 0 && st.substr(st.length()-1,1) != ".") st += ".";
	
	ErrorMess em;
	em.line_num = line_num;
	em.error = st;
	em.type = ERROR_FATAL;
	
	error_mess.push_back(em);	
	if(fatal_error() == true && error_mess.size() >= ERR_MSG_MAX) output_error_messages("Total error limit exceeded");
}


/// Error message for imported file (for specific line)
void Input::alert_line(string st, unsigned int line)                               
{
	if(st.length() > 0 && st.substr(st.length()-1,1) != ".") st += ".";
	
	auto i = 0u; while(i < error_mess.size() && error_mess[i].line_num != line) i++;

	if(i == error_mess.size()){
		ErrorMess em;
		em.line_num = line;
		em.error = st;
		em.type = ERROR_FATAL;
	
		error_mess.push_back(em);	
	}

	if(fatal_error() == true &&  error_mess.size() >= ERR_MSG_MAX){
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
	if(st.length() > 0 && st.substr(st.length()-1,1) != ".") st += ".";
	
	ErrorMess em;
	em.line_num = line_num;
	em.error = st;
	em.type = ERROR_WARNING;
	
	error_mess.push_back(em);	
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
void Input::output_error_messages(string te) const 
{
	for(const auto &em : error_mess){
		if(em.line_num < lines.size()){ 	
			cout << "\033[32m";
			cout <<  "Line " << em.line_num+1 << ": ";
			cout << "\033[0m";
			cout << lines[em.line_num] << endl;
		}

		switch(em.type){
		case ERROR_FATAL: 
			cout << "\033[31m";
			cout << "ERROR: ";
			cout << "\033[0m";
			break;
			
		case ERROR_WARNING:
			cout << "\033[35m";
			cout << "WARNING: ";
			cout << "\033[0m";
			break;
		}
		
		cout << em.error << endl;
		cout << endl;
	}
	
	if(fatal_error() == true){
		cout << te << endl;
		exit (EXIT_FAILURE);
	}
}


/// Processes a given command
void Input::process_command(const CommandLine &cline, unsigned int loop)
{
	cline_store = cline;
	
	auto cname = cline.command;
	
	switch(cname){
	case SPECIES: species_command(loop); break;
	case CLASS: classification_command(loop); break;
	case SET: set_command(); break;
	case CAMERA: camera_command(); break;
	case COMP: compartment_command(); break;
	case TRANS: transition_command(cname); break;
	case SOURCE_CD: transition_command(cname); break;
	case SINK_CD: transition_command(cname); break;
	case DATA_DIR: datadir_command(); break;
	case DESC: description_command(); break;
	case LABEL: label_command(); break;
	case BOX: box_command(); break;
	case PARAM: param_command(); break;
	case DERIVED: derived_command(); break;
	case SIMULATION: simulation_command(); break;
	case INFERENCE: inference_command(); break;
	case IND_EFFECT: ind_effect_command(); break;
	case FIXED_EFFECT: fixed_effect_command(); break;
	
	case INIT_POP_SIM: case ADD_IND_SIM: case REMOVE_IND_SIM: case MOVE_IND_SIM: 
		if(model.mode == SIM) import_data_table_command(cname);
		else alert_warning("Line ignored because not need for simulation");
		break;
	
	case INIT_POP: case ADD_IND: case REMOVE_IND: case MOVE_IND: 
	case INIT_POP_PRIOR: case COMP_DATA: case TRANS_DATA:
	case SOURCE_DATA: case SINK_DATA:
	case TEST_DATA: case POP_DATA: case POP_TRANS_DATA: 
	case SET_TRAPS_DATA: case IND_TRAPPED_DATA:
	case GENETIC_DATA:
		if(model.mode == INF) import_data_table_command(cname);
		else alert_warning("Line ignored because not need for inference");
		break;
	
	case SIM_PARAM: case SIM_STATE: dummy_file_command(); break;
	case INF_PARAM: case INF_STATE: dummy_file_command(); break;
	
	//case SIM_CONST: case INF_CONST: sample_command(); break;
	
	default: alert_import("Command not recognised"); return;
	}
	
	const auto &tags = cline_store.tags;
	for(auto n = 0u; n < tags.size(); n++){
		if(tags[n].done != 1){ 
			alert_warning("Tag '"+tags[n].name+"' not used"); return;
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


/// Error massage if a tag cannot be found
void Input::cannot_find_tag()                                   
{
	alert_import("Cannot find the '"+tag_find+"' tag for '"+cline_store.command_name+"'");
}


/// Checks is a particular value is one of a list of possible options
unsigned int Input::option_error(string na, string te, const vector <string> &pos, const vector <unsigned int> &conv)
{
	if(pos.size() != conv.size()){ alert_import("Options do not match up"); return UNSET;}
	
	auto k = find_in(pos,te);
	if(k != UNSET) return conv[k];
	
	auto st = "'"+na+"' has a value '"+te+"' but should be chosen from one of the following options: ";
	for(auto i = 0u; i < pos.size(); i++){
		if(i != 0) st += ", ";
		st += "'"+pos[i]+"'";
	}
	
	alert_import(st);
	return UNSET;
}


/// Checks if lat and lng are correctly specified
bool Input::check_latlng_error(double lat, double lng)
{
	if(lat > 90 || lat < -90){ alert_import("'lat' must be in the range -90 to 90"); return true;}
	
	if(lng > 180 || lng < -180){ alert_import("'lng' must be in the range -180 to 180"); return true;}

	return false;
}


/// Checks if a string is a number
bool Input::is_number(string num, string tag)
{
	auto x = number(num);
	if(x == UNSET){
		alert_import("'"+tag+"' must be a number"); 
		return false;
	}
	return true;
}


/// Checks if a string is a positive number
bool Input::is_positive(string num, string tag)
{
	auto val = number(num);
	if(val == UNSET || val <= 0){
		alert_import("'"+tag+"' must be a positive number"); 
		return false;
	}
	return true;
}


/// Checks if a string is between zero and one
bool Input::is_zeroone(string num, string tag)
{
	auto val = number(num);
	if(val == UNSET || val < 0 || val > 1){
		alert_import("'"+tag+"' must be between 0 and 1"); 
		return false;
	}
	return true;
}


/// Converts a text string to a prior specification
Prior Input::convert_text_to_prior(string te)
{
	Prior pri;
	
	te = trim(te);
	
	auto spl = split(te,'(');
	if(spl.size() != 2){ pri.error = true; return pri;}
	
	if(spl[1].size() == 0){ pri.error = "Prior has a syntax error"; return pri;}
	if(spl[1].substr(spl[1].length()-1,1) != ")"){ pri.error = "Prior has a syntax error"; return pri;}
	
	auto type = toLower(spl[0]);
	if(type == "dir") type = "dirichlet";
	if(type == "bern") type = "bernoulli";
	
	auto fl = false;
	if(type == "uniform"){ pri.type = UNIFORM_PR; fl = true;}
	if(type == "exp"){ pri.type = EXP_PR; fl = true;}
	if(type == "normal"){ pri.type = NORMAL_PR; fl = true;}
	if(type == "gamma"){ pri.type = GAMMA_PR; fl = true;}
	if(type == "log-normal"){ pri.type = LOG_NORMAL_PR; fl = true;}
	if(type == "beta"){ pri.type = BETA_PR; fl = true;}
	if(type == "bernoulli"){ pri.type = BERNOULLI_PR; fl = true;}
	if(type == "fix"){ pri.type = FIX_PR; fl = true;}
	if(type == "flat"){ pri.type = FLAT_PR; fl = true;}
	if(type == "dirichlet"){ pri.type = DIRICHLET_PR; fl = true;}
	
	if(fl == false){
		pri.error = "Distribution '"+type+"' not recognised"; return pri;
	}
	
	auto bra = spl[1].substr(0,spl[1].length()-1);
	auto spl2 = split(bra,',');

	switch(pri.type){
	case FIX_PR:
		{
			if(spl2.size() != 1){ pri.error = "Expected one value in the brackets"; return pri;}
		}
		break;
	
	case UNIFORM_PR:
		{
			if(spl2.size() != 2){ pri.error = "Expected two values in the brackets"; return pri;}
		}
		break;
	
	case EXP_PR:
		{
			if(spl2.size() != 1){ pri.error = "Expected one value in the brackets"; return pri;}
		}
		break;
		
	case NORMAL_PR:
		{
			if(spl2.size() != 2){ pri.error = "Expected two values in the brackets"; return pri;}
		}
		break;
		
	case GAMMA_PR:
		{
			if(spl2.size() != 2){ pri.error = "Expected two values in the brackets"; return pri;}
		}
		break;
		
	case LOG_NORMAL_PR:
		{
			if(spl2.size() != 2){ pri.error = "Expected two values in the brackets"; return pri;}
		}
		break;
		
	case BETA_PR:
		{
			if(spl2.size() != 2){ pri.error = "Expected two values in the brackets"; return pri;}
		}
		break;

	case BERNOULLI_PR:
		{
			if(spl2.size() != 1){ pri.error = "Expected one value in the brackets"; return pri;}
		}
		break;
		
	case FLAT_PR:
		break;
		
	case DIRICHLET_PR:
		{
			if(spl2.size() != 1){ pri.error = "Expected one value in the brackets"; return pri;}
		}
		break;
	}

	for(auto i = 0u; i < spl2.size(); i++){
		pri.dist_param.push_back(add_equation_info(spl2[i],DIST));
	}
	
	return pri;
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
					auto sup = "^"+ie1+ie2;
					
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
			alert_line("Parameter '"+par.full_name+"' is not used in the model",par.line_num);
		}
	}
}

	
/// Orders parameters in a vector such that only dependent on an params with a smnaller index
void Input::create_param_vector()
{
	auto N = model.param.size();
	
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
			
		if(par.variety != CONST_PARAM){
			for(auto j = 0u; j < par.N; j++){
				par.param_vec_ref[j] = model.param_vec.size();
				
				ParamVecEle pr; 
				pr.name = get_param_name_with_dep(par,par.dep,j);
				pr.th = th; 
				pr.index = j;
				pr.prior = par.prior[j];
				pr.variety = par.variety;
				model.param_vec.push_back(pr);
			}		 
		}
	}
	
	/// Converts parameter references in equation to param_vec
	for(auto &eq : model.eqn){
		for(auto &ca : eq.calc){
			for(auto &it : ca.item){
				if(it.type == PARAMETER){
					if(it.num >= model.param.size()) emsg("Out of range1");
					if(it.index >= model.param[it.num].param_vec_ref.size()) emsg("Out of range2");
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


/// Performs raw processing of data
void Input::raw_process()
{
	vector <string> vecr;
	
	vecr.push_back("Clackmannanshire");
	vecr.push_back("Dumfries and Galloway");
	vecr.push_back("East Ayrshire");
	vecr.push_back("East Lothian");
	vecr.push_back("East Renfrewshire");
	vecr.push_back("Na h-Eileanan Siar");
	vecr.push_back("Falkirk");
	vecr.push_back("Highland");
	vecr.push_back("Inverclyde");
	vecr.push_back("Midlothian");
	vecr.push_back("Moray");
	vecr.push_back("North Ayrshire");
	vecr.push_back("Orkney Islands");
	vecr.push_back("Scottish Borders");
	vecr.push_back("Shetland Islands");
	vecr.push_back("South Ayrshire");
	vecr.push_back("South Lanarkshire");
	vecr.push_back("Stirling");
	vecr.push_back("Aberdeen City");
	vecr.push_back("Aberdeenshire");
	vecr.push_back("Argyll and Bute");
	vecr.push_back("City of Edinburgh");
	vecr.push_back("Renfrewshire");
	vecr.push_back("West Dunbartonshire");
	vecr.push_back("West Lothian");
	vecr.push_back("Angus");
	vecr.push_back("Dundee City");
	vecr.push_back("East Dunbartonshire");
	vecr.push_back("Fife");
	vecr.push_back("Perth and Kinross");
	vecr.push_back("Glasgow City");
	vecr.push_back("North Lanarkshire");


	ofstream fout("Testing/Map data/matrix_col.csv");
	fout << "r',r,Value" << endl;
	
	string file = "Testing/Map data/M_Scotland_LA.csv";
	ifstream in(file);
	if(!in){ return;}
	
	string line;
	
	auto j = 0u;
	do{
		getline(in,line);
		
		if(in.eof()) break;

		remove_cr(line);

		auto vec = split(line,',');
	
		for(auto i = 0u; i < vec.size(); i++){
			fout << "\"" << vecr[j] << "\",";
			fout << "\"" << vecr[i] << "\",";
			fout << vec[i] << endl;
		}
		j++;
	}while(true);

}


/*
/// Finds the possible set of priors for a parameter type
vector <PriorPos> Input::set_pri_pos(type)
{
	let pri_pos = prior_pos;           
	switch(type){
	case "trans_bp": pri_pos = bp_prior_pos; break;
	case "trap_prob": case "Se": case "Sp": pri_pos = zeroone_prior_pos; break;
	}
	
	return pri_pos;
}
*/
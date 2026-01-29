// Simulates data from state

// ./bici-para Execute/init.bici data-sim 'init-pop-inf'
// ./bici-para Execute/init.bici data-show
// ./bici-para Execute/init.bici data-sim ic
// ./bici-para Execute/init.bici data-sim 'init-pop-inf'
// ./bici-para Execute/init.bici data-sim 'remove-pop-inf'
// ./bici-para Execute/init.bici data-sim 'add-pop-inf'
// ./bici-para Execute/init.bici data-sim 'remove-ind-inf'
// ./bici-para Execute/init.bici data-sim 'remove-ind-inf'
// ./bici-para Execute/init.bici data-sim 'comp-data name="comp" class="DS" dt="10"'
// ./bici-para Execute/init.bici data-sim 'trans-data name="trans" trans="I->R"'
// ./bici-para Execute/init.bici data-sim 'test-data name="test" Se=1 Sp=1 comp="I" dt="10"'
// ./bici-para Execute/init.bici data-sim 'pop-data name="pop" filter="DS=file" dt=10 error="poisson"'
// ./bici-para Execute/init.bici data-sim 'pop-trans-data name="poptrans"  trans="I->R"  filter="Sex=file" dt=10 error="normal:1"'
// ./bici-para Execute/init.bici data-sim 'ind-effect-data name="iedata" ie="g"'
// ./bici-para Execute/init.bici data-sim 'pop-trans-data name="poptrans"  trans="I->R"  filter="DS=file" dt=10 error="poisson"'
// gdb ./bici-para Execute/init.bici data-sim 'comp-data'
// run  Execute/init.bici data-sim 'comp-data'
// valgrind --exit-on-first-error=yes --error-exitcode=1 --leak-check=yes -s  ./bici-para Execute/init.bici data-sim 'ic'


#include <iostream>
#include <fstream>
#include <algorithm>
#include <assert.h>
#include <math.h>

using namespace std;

#include "data_sim.hh"
#include "state.hh"
#include "utils.hh"

/// Initilaises the simulation
DataSim::DataSim(Model &model, Output &output) : model(model), output(output)
{	
}


/// Performs a simulation of a data source
void DataSim::run(string data_sim_line)
{
	data_sim_line_store = data_sim_line;
	
	State state(model);
	state.init();
	
	line = get_command_tags(data_sim_line,UNSET);
	if(line.emsg != "") error(line.emsg);
	
	auto cname = line.command;
	
	switch(cname){
	case IC_DATA:
	case INIT_POP: case ADD_POP: case REMOVE_POP: 
	case ADD_IND: case REMOVE_IND: case MOVE_IND:
	case COMP_DATA: case TRANS_DATA: case TEST_DATA:
	case POP_DATA: case POP_TRANS_DATA:
	case IND_EFFECT_DATA: case IND_GROUP_DATA: case GENETIC_DATA:
		break;
		
	default:
		error("BICI command '"+line.command_name+"' cannot be simulated");
		break;
	}
		
	if(model.sample.size() == 0) error("No simulated samples exist (run BICI using 'sim')");
		
	const auto &samp = model.sample[0];	
	auto param_val = model.post_param(samp);	
		
	state.load_samp(param_val,samp);
		
	auto p = UNSET;
	if(model.species.size() == 1) p = 0;
	else{
		auto spe = get_tag_value("species"); 	if(spe == "") cannot_find_tag(); 
		p = model.find_p(spe);
		if(p == UNSET) error("Species '"+spe+"' is not reconginsed");
	}
		
	auto file = get_tag_value("file"); 	
	if(file == ""){
		if(output.datadir != ""){
			error("A 'file' tag must be added to give the filename used to store the data.");
		}		
	}
	else{	
		if(output.datadir == ""){
			error("A 'file' tag is not required because results are embedded into the BICI-script.");
		}
	}
	
	const auto &sp = model.species[p];
	
	//if(sp.type == POPULATION || sp.type == DETERMINISTIC) state.species[p].set_cpop_st();
	
	if(sp.type != INDIVIDUAL){
		switch(cname){
		case ADD_IND: case REMOVE_IND: case MOVE_IND:
		case COMP_DATA: case TRANS_DATA: case TEST_DATA:
		case IND_EFFECT_DATA: case IND_GROUP_DATA: case GENETIC_DATA:
			error("Cannot simulate individual-based data from a population-based species");
			break;
			
		default: break;
		}
	}
	
	if(cname == IC_DATA){ // Copies across initial condition data
		// First removes any existing
		vector <unsigned int> existing;
		for(auto k = 0u; k < sp.source.size(); k++){
			const auto &so = sp.source[k];
			switch(so.cname_raw){
			case INIT_POP: case ADD_POP: case REMOVE_POP:
			case ADD_IND: case REMOVE_IND: case MOVE_IND:
				existing.push_back(k);
				break;
			default: break;
			}
		}
		
		if(existing.size() > 0){
			auto res = model.question("This will remove existing data. Is that OK?");
			
			if(res == false){
				cout << "Terminated" << endl; 
				exit (EXIT_FAILURE);
			}
		}
		
		// Deletes any existing definitions
		for(auto k : existing){
			const auto &so = sp.source[k];
			output.delete_command(so.name,so.line_num,false);
		}
		
		for(int k = sp.source.size()-1; k >= 0; k--){
			const auto &so = sp.source[k];
			switch(so.cname_raw){
			case INIT_POP_SIM: case ADD_POP_SIM: case REMOVE_POP_SIM:
			case ADD_IND_SIM: case REMOVE_IND_SIM: case MOVE_IND_SIM:
				{
					string name;
					switch(so.cname_raw){
					case INIT_POP_SIM: name = "init-pop-inf"; break;
					case ADD_POP_SIM: name = "add-pop-inf"; break;
					case REMOVE_POP_SIM: name = "remove-pop-inf"; break;
					case ADD_IND_SIM: name = "add-ind-inf"; break;
					case REMOVE_IND_SIM: name = "remove-ind-inf"; break;
					case MOVE_IND_SIM: name = "move-ind-inf"; break;
					default: break;
					}
			
					auto st = name;
					for(const auto &ta : so.tags){
						if(ta.name != "file"){
							st += " "+ta.name+"=\""+ta.value+"\"";
						}
					}
					
					output.insert_command("",p,INF_BANNER,PARAM_BANNER,st,output_table(so.table),file); 
				}
				break;
			default: break;
			}
		}
	}
	else{
		auto name = get_tag_value("name"); 	
		if(name == "") get_default_name(name,cname);
		if(name == "") cannot_find_tag(); 
			
		string te="";
		
		Table tab;
		
		switch(cname){
		case INIT_POP: init_pop_sim(p,state,tab); break;
		case ADD_IND: add_ind_sim(p,state,tab); break;
		case ADD_POP: copy_from_sim(p,"add-pop-sim",ADD_POP_SIM,tab); break;
		case REMOVE_POP: copy_from_sim(p,"remove-pop-sim",REMOVE_POP_SIM,tab); break;
		case REMOVE_IND: copy_from_sim(p,"remove-ind-sim",REMOVE_IND_SIM,tab); break;
		case MOVE_IND: copy_from_sim(p,"move-ind-sim",MOVE_IND_SIM,tab); break;
		
		case COMP_DATA: comp_data_sim(p,state,tab); break;
		case TRANS_DATA: trans_data_sim(p,state,tab); break;
		case TEST_DATA: test_data_sim(p,state,tab); break;
		case POP_DATA: pop_data_sim(p,state,tab); break; 
		case POP_TRANS_DATA: pop_trans_data_sim(p,state,tab); break; 
		case IND_EFFECT_DATA: ind_effect_data_sim(p,state,tab); break; 
		case IND_GROUP_DATA: ind_group_data_sim(p,state,tab); break; 
		case GENETIC_DATA: genetic_data_sim(p,state,tab); break;
		default: emsg("Simulation prob"); break;
		}

		auto st = line.command_name;
		auto &tags = line.tags;
		
		auto i = 0u;
		while(i < tags.size()){
			const auto &ta = tags[i];
			if(!ta.done){
				cout << "\033[35m";
				cout << "WARNING: ";
				cout << "\033[0m";
				cout << "In 'data-sim' tag '" << ta.name << "' has not been used." << endl;
				
				tags.erase(tags.begin()+i);
			}
			else{
				if(ta.name != "dt" && ta.name != "times" && ta.name != "frac" && ta.name != "ac"){
					st += " "+ta.name+"=\""+ta.value+"\"";
				}
				i++;
			}
		}
		
		if(tab.ele.size() > 0){
			output.insert_command(name,p,INF_BANNER,PARAM_BANNER,st,output_table(tab),file); 
		}
	}
}


/// Gets the value of a tag
string DataSim::get_tag_value(string st)                                  
{
	auto &tags = line.tags;

	tag_find = st;
	for(auto i = 0u; i < tags.size(); i++){
		if(tags[i].name == st){ tags[i].done = 1; return tags[i].value;}
	}
	
	return "";
}

/// Error massage if a tag cannot be found
void DataSim::cannot_find_tag() const                                 
{
	string te = "Cannot find the '"+tag_find+"' tag for '"+line.command_name+"'";
	error(te);
}


/// Dislays an error message
void DataSim::error(string msg) const
{
	add_full_stop(msg);
	
	cout << endl;
	cout << "\033[32m";
	cout <<  "'";
	switch(model.mode){
	case DATA_SIM: cout << "data-sim"; break;
	case DATA_SHOW: cout << "data-show"; break;
	case DATA_DEL: cout << "data-del"; break;
	case DATA_CLEAR: cout << "clear"; break;
	default: emsg("Cannot be here"); break;
	}
	cout << "': ";
	cout << "\033[0m";
	cout << data_sim_line_store << endl;
				
	cout << "\033[31m" << "ERROR: ";
	cout << "\033[0m";
	cout << msg;
	cout << endl;
	
	end_code();
}


/// Checks is a particular value is one of a list of possible options
unsigned int DataSim::option_error(string na, string te, const vector <string> &pos, const vector <unsigned int> &conv) const
{
	if(pos.size() != conv.size()){ 
		error("Options do not match up"); 
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
	
	error(st);
	return UNSET;
}


/// Outputs a simple table 
string DataSim::output_table(const Table &tab) const
{
	string st = "";
	for(auto i = 0u; i < tab.heading.size(); i++){
		if(i != 0) st += ",";
		st += "\""+tab.heading[i]+"\"";
	}
	st += endli;

	for(auto r = 0u; r < tab.ele.size(); r++){
		auto ele = tab.ele[r];
		if(ele.size() != tab.heading.size()) emsg("Output columns to not match"); 

		for(auto i = 0u; i < ele.size(); i++){
			if(i != 0) st += ',';
			if(number(ele[i]) == UNSET) st += "\""+ele[i]+"\"";
			else st += ele[i];
		}
		st += endli;
	}
	
	return st;
}


/// Shows any loaded data
void DataSim::show() const
{
	auto list = get_data_list();
	
	if(list.size() == 0) cout << "No data." << endl;
	
	for(auto ty = 0u; ty < 3; ty++){
		auto type = DataSourceType(ty);
	
		auto k = 0u;
		while(k < list.size() && list[k].type != type) k++;
		if(k < list.size()){
			switch(ty){
			case SIM_DATA: cout << "SIMULATION DATA" << endl; break;
			case INF_DATA: cout << "INFERENCE DATA" << endl; break;
			case POST_SIM_DATA: cout << "POST-SIM DATA" << endl; break;
			}
			
			for(auto p = 0u; p < model.nspecies; p++){
				const auto &sp = model.species[p];
				
				auto k = 0u; 
				while(k < list.size() && !(list[k].type == type && list[k].p == p)) k++;
				if(k < list.size()){
					if(model.nspecies > 1){
						cout << "For species '" << sp.name << "':" << endl;
					}
					for(auto k = 0u; k < list.size(); k++){
						const auto &dr = list[k];
						if(dr.type == type && dr.p == p){
							const auto &so = sp.source[dr.i];
							auto desc = source_description(sp,so);
							cout << "(" << k+1 << ") " << so.name;
							if(desc != "") cout << " (" << desc << ")";
							cout << " active: " << so.active;
							cout << endl;
						}
					}
				}
			}
			cout << endl;
		}
	}
}


/// Organises data into sections
vector <DataRef> DataSim::get_data_list() const
{
	vector <DataRef> list;
	
	for(auto ty = 0u; ty < 3; ty++){
		auto type = DataSourceType(ty);
		
		for(auto p = 0u; p < model.species.size(); p++){
			const auto &sp = model.species[p];
			for(auto i = 0u; i < sp.source.size(); i++){
				auto add = false;
				
				const auto &so = sp.source[i];
				switch(so.cname_raw){
				case INIT_POP_SIM:
				case ADD_POP_SIM: case REMOVE_POP_SIM: 
				case ADD_IND_SIM: case REMOVE_IND_SIM: case MOVE_IND_SIM: 
					if(type == SIM_DATA) add = true;
					break;
					
				case INIT_POP: 
				case ADD_POP: case REMOVE_POP:
				case ADD_IND: case REMOVE_IND: case MOVE_IND: 
				case COMP_DATA: case TRANS_DATA:
				case TEST_DATA: case POP_DATA: case POP_TRANS_DATA: 
				case GENETIC_DATA:
					if(type == INF_DATA) add = true;
					break;
	
				case ADD_POP_POST_SIM: case REMOVE_POP_POST_SIM: 
				case ADD_IND_POST_SIM: case REMOVE_IND_POST_SIM: case MOVE_IND_POST_SIM: 
					if(type == POST_SIM_DATA) add = true;
					break;
					
				default: emsg("Data type not recognised"); break;
				}
				
				if(add){
					DataRef dr; dr.type = type; dr.p = p; dr.i = i;
					list.push_back(dr);
				}
			}
		}			
	}
	
	return list;
}


/// Deletes a data source
void DataSim::del(string data_sim_line)
{	
	data_sim_line_store = data_sim_line;

	auto dsl = data_sim_line;
	
	if(dsl == "all"){
		vector <IndRef> list;
		for(auto p = 0u; p < model.nspecies; p++){
			const auto &sp = model.species[p];
			for(auto i = 0u; i < sp.source.size(); i++){
				const auto &so = sp.source[i];
				
				switch(so.cname_raw){
				case INIT_POP: case ADD_POP: case REMOVE_POP: 
				case ADD_IND: case REMOVE_IND: case MOVE_IND:
				case COMP_DATA: case TRANS_DATA: case TEST_DATA:
				case POP_DATA: case POP_TRANS_DATA:
				case IND_EFFECT_DATA: case IND_GROUP_DATA: case GENETIC_DATA:
					{
						IndRef po; po.p = p; po.i = i;
						list.push_back(po);
					}
					break;
				
				default: break;
				}	
			}
		}
		
		if(list.size() == 0) error("There was no inference data to delete");
		else{
			auto res = model.question("TAre you sure you want to delete all inference data?");
			
			if(res == false){
				cout << "Terminated" << endl; 
				exit (EXIT_FAILURE);
			}
			
			for(const auto &po : list){
				const auto &so = model.species[po.p].source[po.i];
				output.delete_command(so.name,so.line_num,false);
			}
		}
	}
	else{
		auto p_sel = UNSET;
		if(model.nspecies == 1) p_sel = 0;
		
		auto i = 0u; 
		while(i < data_sim_line.length() && data_sim_line.substr(i,1) != ":" && data_sim_line.substr(i,1) != "\"") i++;
		
		if(i < data_sim_line.length() && data_sim_line.substr(i,1) == ":"){
			auto sp_na = trim(data_sim_line.substr(0,i));
			auto p = model.find_p(sp_na);
			if(p != UNSET) p_sel = p;
			else error("Species name '"+sp_na+"' not recognised");
			
			data_sim_line = data_sim_line.substr(i);
		}

		vector <IndRef> pos;
		
		auto num = number(data_sim_line);
		if(num != UNSET){
			auto list = get_data_list();
			
			auto numi = (unsigned int) num;
			if(num < 1 || num > list.size() || numi != num){
				error("'"+data_sim_line+"' is not a valid number");
			}
			
			const auto &dr = list[numi-1];
			IndRef po; po.p = dr.p; po.i = dr.i;
			pos.push_back(po);
		}
		else{
			auto p_min = 0u, p_max = model.nspecies;
			if(p_sel != UNSET){ p_min = p_sel; p_max = p_sel+1;}
				
			for(auto p = p_min; p < p_max; p++){
				const auto &sp = model.species[p];
				
				for(auto i = 0u; i < sp.source.size(); i++){
					const auto &so = sp.source[i];
					if(begin_str(so.name,data_sim_line)){
						IndRef po; po.p = p; po.i = i;
						pos.push_back(po);
					}
				}
			}
		}
		
		if(pos.size() == 0) error("No data matches for '"+dsl+"'");
		
		if(pos.size() > 1) error("Multiple data sources match '"+dsl+"'");
		
		{
			const auto &so = model.species[pos[0].p].source[pos[0].i];
			output.delete_command(so.name,so.line_num,true);
		}
	}
}


/// Clears output data from sim/inf/post-sim
void DataSim::clear(string data_sim_line)
{
	data_sim_line_store = data_sim_line;
	
	auto fl = false;

	vector <string> list, banner;
	if(data_sim_line == "sim" || data_sim_line == "all"){ 
		fl = true;
		list.push_back("param-sim"); 
		list.push_back("state-sim"); 
		list.push_back("warning-sim"); 
		banner.push_back(OUT_SIM_BANNER);
	}
	
	if(data_sim_line == "inf" || data_sim_line == "all"){
		fl = true;
		list.push_back("param-inf"); 
		list.push_back("state-inf"); 
		list.push_back("generation-inf");
		list.push_back("trans-diag-inf"); 
		list.push_back("proposal-inf"); 
		list.push_back("diagnostics-inf"); 
		list.push_back("param-stats-inf"); 
		list.push_back("warning-inf"); 
		banner.push_back(OUT_INF_BANNER);
	}
	
	if(data_sim_line == "post-sim" || data_sim_line == "inf" || data_sim_line == "all"){
		fl = true;
		list.push_back("param-post-sim"); 
		list.push_back("state-post-sim"); 
		list.push_back("warning-post-sim"); 
		banner.push_back(OUT_POST_SIM_BANNER);
	}
	
	if(!fl) error("The value after 'clear' must be 'sim', 'inf', 'post-sim' or 'all'");
	
	output.delete_commands(data_sim_line,list);
	output.delete_lines(banner);
	output.remove_double_space();
}


/// Gives a description of the data source
string DataSim::source_description(const Species &sp, const DataSource &so) const 
{
	string desc;

	const auto &tab = so.table;
	auto c = tab.ncol-1;
	
	auto num=tab.nrow;
	
	switch(so.cname){
	case INIT_POP:
		switch(so.init_pop_type){
		case INIT_POP_FIXED:		
			if(so.focal_cl != UNSET){
				auto cl = so.focal_cl;
				
				num = 0;
				for(auto r = 0u; r < tab.nrow; r++){
					auto val = tab.ele[r][c];
					if(!is_percent(val)){
						auto nu = number(val);
						if(nu != UNSET) num += nu;		
					}						
				}
				if(sp.ncla > 1) desc += "focal '"+sp.cla[cl].name+"'";
			}
			else{
				num = 0;
				for(auto r = 0u; r < tab.nrow; r++){
					auto nu = number(tab.ele[r][c]);
					if(nu != UNSET) num += nu;						
				}
			}
			break;
			
		case INIT_POP_DIST:
			if(so.focal_cl != UNSET){
				if(sp.ncla > 1) desc += "focal '"+sp.cla[so.focal_cl].name+"'";
				if(desc != "") desc += " ";
				desc += "distribution";
			}
			else{
				desc += "N ~ "+ get_prior_string(so.pop_prior);
			}
			num = UNSET;
			break;

			case INIT_POP_NONE: break;		
		}
		break;

	case ADD_IND: case REMOVE_IND: 
		break;
	
	case ADD_POP: case REMOVE_POP: 
		num = 0;
		for(auto r = 0u; r < tab.nrow; r++){
			auto nu = number(tab.ele[r][c]);
			if(nu != UNSET) num += nu;						
		}
		break;
		
	case MOVE_IND:
		break;

	case COMP_DATA: 
		desc += so.filter_str;
		break;
		
	case TEST_DATA:
		{
			const auto &om = so.obs_model;
			const auto &dts = om.diag_test_sens;
		
			string sens;
			for(auto i = 0u; i < dts.comp.size(); i++){
				if(i != 0) sens += ",";
				sens += sp.cla[dts.cl].comp[dts.comp[i]].name;
			}
			desc = sp.cla[dts.cl].name+" "+sens+" Se="+om.Se_str+" Sp="+om.Sp_str;
		}
		break;
	
	case IND_EFFECT_DATA:
		break;
		
	case IND_GROUP_DATA:
		break;
	
	case POP_DATA:	
		desc = so.filter_str;
		num = tab.nrow;
		break;
		
	case TRANS_DATA: case POP_TRANS_DATA:
		desc = so.filter_trans_str;
		break;

	case GENETIC_DATA:
		switch(so.gen_data_type){
		case MATRIX_DATA: desc = "Matrix"; break;
		case SNP_DATA: desc = "# SNP "+(tab.ncol-2); break;
		}
		break;
	
	default: desc = "no description"; break;
	}

	if(desc != "") desc = "details: "+desc;
	
	if(desc != "") desc = ", "+desc;
	desc = "command: "+so.command_name+desc;

	if(num != UNSET){
		if(desc != "") desc += ", ";
		desc += "number:"+tstr(num);
	}
	return desc;
}


/// Loads up the observation model
ObsModel DataSim::load_obs_model()
{
	ObsModel om;
	
	auto error_str = get_tag_value("error"); if(error_str == "") cannot_find_tag(); 
			
	auto spl = split(error_str,':');

	auto str = toLower(spl[0]);
	if(str == "normal"){
		if(spl.size() != 2) error("For 'error' problem with expression '"+error_str+"'");
	
		auto val = trim(spl[1]);
		
		if(toLower(val) == "file"){
			om.type = NORMAL_FILE_OBSMOD;
		}
		else{
			if(is_percent(val) == true){
				om.type = NORMAL_PERCENT_OBSMOD;
				om.percent = number(val.substr(0,val.length()-1));
				if(om.percent == UNSET){
					error("For 'error' the expression '"+tstr(om.percent)+"' is not a percentage");
				}
				
				if(om.percent <= 0){
					error("For 'error' the value '"+tstr(om.percent)+"%' must be a positive percentage");
				}
			}	
			else{
				om.type = NORMAL_SD_OBSMOD;
				om.sd = number(val);
				if(om.sd == UNSET){
					error("For 'error' the expression '"+val+"' must be a number or a percentage");
				}
				
				if(om.sd <= 0){
					error("For 'error' the value '"+tstr(om.sd)+"' must be positive");
				}
			}
		}
	}
	else{
		if(str == "poisson"){
			if(spl.size() != 1){
				error("For 'error' there is a problem with expression '"+error_str+"'");
			}
			
			om.type = POISSON_OBSMOD;
		}
		else{
			if(str == "neg-binomial"){
				if(spl.size() != 2){
					error("For 'error' there is a problem with expression '"+error_str+"'");
				}
				
				auto val = trim(spl[1]);
			
				if(toLower(val) == "file"){
					om.type = NEGBIN_FILE_OBSMOD;
				}
				else{
					om.type = NEGBIN_OBSMOD;
					om.p = number(val);
						
					if(om.p == UNSET){
						error("For 'error' the expression '"+val+"' must be a number");
					}
					
					if(om.p <= 0 || om.p >= 1){
						error("For 'error' the value '"+val+"' must be between zero and one");
					}
				}
			}
			else{
				error("For 'error' the observation error '"+str+"' is not recognised. It must be chosen from the following options: 'normal', 'poisson' or 'neg-binomial'.");
			}
		}
	}
	
	return om;
}


/// Work out how filter is stratified
Stratify DataSim::get_stratify(const Filter &cf, unsigned int p) const
{
	const auto &sp = model.species[p];
	
	Stratify strat;
	
	vector <unsigned int> max;
	auto &cl_list = strat.cl_list;
	for(auto cl = 0u; cl < sp.ncla; cl++){
		if(cf.cla[cl].type == FILE_FILT){
			cl_list.push_back(cl);
			max.push_back(sp.cla[cl].comp.size());
		}
	}
	
	auto N = cl_list.size();
	vector <unsigned int> index(N,0);

	do{
		strat.comb.push_back(index);
		
		auto i = 0u;
		auto fl = true;
		while(fl == true && i < N){
			fl = false;
			index[i]++; if(index[i] == max[i]){ index[i] = 0; i++; fl = true;}
		}
		if(i == N) break;
	}while(true);
	
	return strat;
}


/// For a given stratification set the filter
void DataSim::set_filt_strat(unsigned int com, Filter &cf, Stratify strat, unsigned int p) const
{
	const auto &sp = model.species[p];
	
	for(auto i = 0u; i < strat.cl_list.size(); i++){
		auto cl = strat.cl_list[i];
		auto &cf_cl = cf.cla[cl];
		
		auto c_sel = strat.comb[com][i];
		
		cf_cl.type = COMP_FILT;
		const auto &comp = sp.cla[cl].comp;
		cf_cl.comp_prob_str.resize(comp.size());
		for(auto c = 0u; c < sp.cla[cl].comp.size(); c++){
			if(c == c_sel) cf_cl.comp_prob_str[c] = "1";
			else cf_cl.comp_prob_str[c] = "0";
		}
	}
}


/// Gets to probability of being in different global compartments
vector <double> DataSim::get_comp_prob(const Filter &cf, double tdiv, unsigned int cl_sel, unsigned int p) const
{
	const auto &sp = model.species[p];
	
	vector < vector <double> > prob;
	prob.resize(sp.ncla);
	for(auto cl = 0u; cl < sp.ncla; cl++){
		if(cl != cl_sel){
			for(auto c = 0u; c < sp.cla[cl].comp.size(); c++){
				auto probc_str = cf.cla[cl].comp_prob_str[c];
								
				string err;
				auto probc = model.calculate_equation(probc_str,tdiv,err);
				if(err != ""){
					error("An error was generated for filter '"+probc_str+"': "+err);
				}
				
				if(probc < 0){
					error("An error was generated for filter '"+probc_str+"': The value is negative");
				}

				prob[cl].push_back(probc);								
			}
		}
	}

	vector <double> res;
	for(auto c = 0u; c < sp.comp_gl.size(); c++){
		const auto &cgl = sp.comp_gl[c];
		
		auto pr = 1.0;
		for(auto cl = 0u; cl < sp.ncla; cl++){
			if(cl != cl_sel) pr *= prob[cl][cgl.cla_comp[cl]];
		}
		res.push_back(pr);
	}
	
	return res;
}
			
			
/// Gets a vector of transition probabilities
vector <double> DataSim::get_tr_prob(const vector <string> &trans_filt, double tdiv) const
{
	vector <double> res;
	
	for(auto tr = 0u; tr < trans_filt.size(); tr++){
		auto str = trans_filt[tr];
		
		string err;
		auto prob = model.calculate_equation(str,tdiv,err);
		if(err != ""){
			error("An error was generated for filter '"+str+"': "+err);
		}
		if(prob < 0){
			error("An error was generated for filter '"+str+"': Number negative");
		}
		
		res.push_back(prob);
	}
	
	return res;
}

		
/// Samples from the observation model
double DataSim::obs_model_sample(double mean, const ObsModel &om) const
{
	if(mean < 0) mean = 0;
	
	string warn;
	double val;
	val = UNSET;
	
	switch(om.type){
	case NORMAL_PERCENT_OBSMOD: case NORMAL_SD_OBSMOD:
		{
			double sd;
			if(om.type == NORMAL_PERCENT_OBSMOD) sd = mean*om.percent/100;
			else sd = om.sd;
			
			auto loop = 0u, loopmax = 1000u;
			do{ val = normal_sample(mean,sd,warn); loop++;}while(val < 0 && loop < loopmax);
			if(loop == loopmax) warn = "Negative sampling problem";
		}
		break;
	
	case POISSON_OBSMOD:
		val = poisson_sample(mean,warn);
		break;
		
	case NEGBIN_OBSMOD:
		val = neg_binommial_sample(mean,om.p,warn);
		break;
		
	case NORMAL_FILE_OBSMOD: case NEGBIN_FILE_OBSMOD:
		error("'file' cannot be used in observation model");
		break;
	}
	
	if(warn != "") error("Problem with observation model: "+warn);
	
	return val;	
}


/// Simulates initial population
void DataSim::init_pop_sim(unsigned int p, const State &state, Table &tab)
{
	const auto &sp = model.species[p];
	const auto &ssp = state.species[p];
	
	auto focal_cl = UNSET;
	
	const auto &source = sp.source;
	
	auto i = 0u; 
	while(i < source.size() && source[i].cname_raw != INIT_POP_SIM) i++;
	if(i < source.size()){
		const auto &so = source[i];
		if(so.focal_cl != UNSET) focal_cl = so.focal_cl;
	}
	
	auto focal = get_tag_value("focal");
	if(focal != ""){
		error("When simulating 'focal' is not defined (because it is automatically determined).");
	}
	
	if(focal_cl != UNSET){
		Tag ta; ta.name = "focal"; ta.value = sp.cla[focal_cl].name;
		line.tags.insert(line.tags.begin(),ta);
	}
	
	auto type = get_tag_value("type");
	if(type != ""){
		error("When simulating 'type' is not defined (because it is automatically 'fixed')");
	}
	
	Tag ta; ta.name = "type"; ta.value = "fixed";
	line.tags.insert(line.tags.begin(),ta);
	
	if(focal_cl == UNSET){
		for(auto cl = 0u; cl < sp.ncla; cl++) tab.heading.push_back(sp.cla[cl].name);
		tab.heading.push_back("Population");

		for(auto c = 0u; c < sp.comp_gl.size(); c++){
			const auto &cgl = sp.comp_gl[c];
			
			vector <string> row;		
			for(auto cl = 0u; cl < sp.ncla; cl++){
			 row.push_back(sp.cla[cl].comp[cgl.cla_comp[cl]].name);
			}
			row.push_back(tstr(ssp.init_cond_val.cnum[c]));
			tab.ele.push_back(row);
		}
	}
	else{
		tab.heading.push_back("Compartment");
		tab.heading.push_back("Population");
		
		vector < vector <unsigned int> > comp_split;
		comp_split.resize(sp.ncla);
		for(auto cl = 0u; cl < sp.ncla; cl++){
			const auto &claa = sp.cla[focal_cl];
			comp_split[cl].resize(claa.ncomp,0);
		}
		
		auto num_tot = 0u;
		for(auto c = 0u; c < sp.comp_gl.size(); c++){
			const auto &cgl = sp.comp_gl[c];
			
			auto num = ssp.init_cond_val.cnum[c];
			
			for(auto cl = 0u; cl < sp.ncla; cl++){
				comp_split[cl][cgl.cla_comp[cl]] += num; 
			}
			num_tot += num;
		}
		
		const auto &claa = sp.cla[focal_cl];
		for(auto c = 0u; c < claa.ncomp; c++){
			vector <string> row;		
			row.push_back(claa.comp[c].name);
			row.push_back(tstr(comp_split[focal_cl][c]));
			tab.ele.push_back(row);
		}
		
		for(auto cl = 0u; cl < sp.ncla; cl++){
			if(cl != focal_cl){
				const auto &claa = sp.cla[cl];
			
				for(auto c = 0u; c < claa.ncomp-1; c++){
					auto fr = double(comp_split[cl][c])/num_tot;
					
					vector <string> row;		
					row.push_back(claa.comp[c].name);
					row.push_back(tstr(fr*100)+"%");
					tab.ele.push_back(row);
				}
			}
		}
	}
}


/// Simulates initial population
void DataSim::add_ind_sim(unsigned int p, const State &state, Table &tab)
{
	const auto &sp = model.species[p];
	const auto &ssp = state.species[p];
	
	tab.heading.push_back("ID");
	tab.heading.push_back("t");
	for(auto cl = 0u; cl < sp.ncla; cl++) tab.heading.push_back(sp.cla[cl].name);
	
	for(const auto &ind : ssp.individual){
		if(ind.ev.size() > 0){
			const auto &ev = ind.ev[0];
			if(ev.type == ENTER_EV){
				vector <string> row;
				row.push_back(ind.name);
			
				row.push_back(tstr(model.calc_t(ev.tdiv)));
				
				const auto &cgl = sp.comp_gl[ev.c_after];
				for(auto cl = 0u; cl < sp.ncla; cl++){
					row.push_back(sp.cla[cl].comp[cgl.cla_comp[cl]].name);
				}
				tab.ele.push_back(row);
			}
		}
	}
}


/// Copies a data source from simulation
void DataSim::copy_from_sim(unsigned int p, string st, Command cname, Table &tab) const 
{
	const auto &sp = model.species[p];
	
	const auto &source = sp.source;
	
	auto i = 0u; 
	while(i < source.size() && source[i].cname_raw != cname) i++;
	if(i < source.size()){
		tab = source[i].table;
	}
	else{
		error("Could not find '"+st+"' data source to copy");
	}
}


/// Gets times for observations
vector <double> DataSim::get_times()
{
	auto dt = get_tag_value("dt"); 
	
	auto times = get_tag_value("times"); 
	
	if(dt == "" && times == "") error("Either 'dt' or 'times' must be set");
	if(dt != "" && times != "") error("'dt' and 'times' cannot both be set");
	
	auto t_start = model.details.t_start;
	auto t_end = model.details.t_end;
	
	vector <double> tims;
	
	if(dt != ""){
		auto ddt = number(dt);
		if(ddt == UNSET || ddt <= 0){
			error("In 'dt' the value '"+dt+"' must be a positive number");
		}

		for(auto t = t_start; t <= t_end; t += ddt){
			tims.push_back(t);
		}
	}
	else{
		auto spl = split(times,',');
		for(auto i = 0u; i < spl.size(); i++){
			auto val = spl[i];
			auto num = number(val);
			if(num == UNSET){
				error("In 'times' the value '"+val+"' is not a number");
			}
			else{
				if(num < t_start || num > t_end){
					error("In 'times' the value "+val+" is not between the start and end times");
				}
				
				tims.push_back(num);
			}
		}
	}
	
	return tims;
}


/// Gets the fraction of observations
double DataSim::get_frac()
{
	auto frac = get_tag_value("frac"); 
	auto fr = 1.0;
	if(frac != ""){
		fr = number(frac);
		if(fr == UNSET){
			error("In 'frac' the value '"+frac+"' is not a number");
		}
		
		if(fr <= 0 || fr > 1){
			error("In 'frac' the value '"+frac+"' must be a number between 0 and 1");
		}
	}
	
	return fr;
}


/// Simulates compartmental data
void DataSim::comp_data_sim(unsigned int p, const State &state, Table &tab)
{
	const auto &sp = model.species[p];			
	
	auto name = get_tag_value("class");
	if(name == ""){
		if(sp.ncla == 1) name = sp.cla[0].name;
		else{ cannot_find_tag(); return;}
	}
			
	auto cl = model.find_cl(p,name);
	if(cl == UNSET) error("In 'class' the classification '"+name+"' is not found");
		
	auto accu = get_tag_value("ac"); 
	auto acc = 1.0;
	if(accu != ""){
		acc = number(accu);
		if(acc == UNSET) error("In 'ac' the value '"+accu+"' is not a number");
		if(acc <= 0 || acc > 1) error("In 'ac' the value '"+accu+"' must be a number between 0 and 1");
	}
	
	auto fr = get_frac();
	
	auto tims = get_times();
		
	tab.heading.push_back("ID");
	tab.heading.push_back("t");
	tab.heading.push_back(name);
	
	const auto &ssp = state.species[p];
	
	for(const auto &ind : ssp.individual){
		for(auto k = 0u; k < tims.size(); k++){
			if(ran() < fr){	
				auto t = tims[k];
				auto tdiv = model.calc_tdiv(t);
			 
				auto c = UNSET;
				auto e = 0u; 
				while(e < ind.ev.size() && ind.ev[e].tdiv <= tdiv){
					c = ind.ev[e].c_after;
					e++;
				}
			
				string te;
			
				vector <string> row;	
				row.push_back(ind.name);
				row.push_back(tstr(t));
				
				if(c == UNSET){
					te = "!";
				}
				else{
					auto cgl = sp.comp_gl[c];
					const auto &comp = sp.cla[cl].comp;
					auto ci = cgl.cla_comp[cl];
				
					if(acc == 1){
						te = comp[ci].name;
					}
					else{
						te = comp[ci].name+":"+tstr(acc);
						if(comp.size() > COMP_NOISY_MAX){  // If many comparmtents then randomly chooses some 
							vector <unsigned int> list;
							for(auto cii = 0u; cii < comp.size(); cii++){
								if(cii != ci) list.push_back(cii);
							}
							
							auto ac_other = precision(((1-acc)/(COMP_NOISY_MAX-1)),3);
							for(auto k = 0u; k < COMP_NOISY_MAX-1; k++){
								auto m = (unsigned int)(ran()*list.size());
								auto cii = list[m];
								te += "|"+comp[cii].name+":"+ac_other;
								list[m] = list[list.size()-1];
								list.pop_back();
							}
						}
						else{
							auto ac_other = precision(((1-acc)/(comp.size()-1)),3);
							for(auto cii = 0u; cii < comp.size(); cii++){
								if(cii != ci){
									te += "|"+comp[cii].name+":"+ac_other;
								}
							}
						}
					}
				}
				row.push_back(te);
				
				tab.ele.push_back(row);
			}
		}
	}
}

/// Simulates transition data
void DataSim::trans_data_sim(unsigned int p, const State &state, Table &tab)
{
	const auto &sp = model.species[p];
	
	auto tmin = model.details.t_start;
	auto tmax = model.details.t_end;
	
	auto filter_str = get_tag_value("filter");
	
	auto trans = get_tag_value("trans"); if(trans == "") cannot_find_tag(); 
		
	auto cl_sel = model.get_cl_from_trans(trans,p);
	if(cl_sel == UNSET){
		error("For 'trans' the value '"+trans+"' is not recognised");
		return;
	}
	
	auto obsran = toLower(get_tag_value("obsrange"));
	if(obsran == "") obsran = "all";
	
	auto time_range = TimeRange(option_error("obsrange",obsran,{"all","specify","file"},{ ALL_TIME, SPEC_TIME, FILE_TIME }));
	if(time_range == UNSET) return;
		
	switch(time_range){
	case ALL_TIME: case FILE_TIME: 
		break;
		
	case SPEC_TIME:
		{
			auto start = get_tag_value("start"); if(start == "") cannot_find_tag();
			tmin = number(start);
			if(tmin == UNSET) error("In 'start' the value '"+start+"' is not a number");
			
			auto end = get_tag_value("end"); if(end == "") cannot_find_tag();
			tmax = number(end);
			if(tmax == UNSET) error("In 'end' the value '"+end+"' is not a number");
			
			if(tmin >= tmax){
				error("The start time '"+start+"' must be before the end time '"+end+"'"); 
			}
		}
		break;
	}
	
	string errmsg;
	auto cf = sp.set_comp_filt_raw(filter_str,UNSET,LOWER_UPPER_BOUND,errmsg);
	if(errmsg != "") error("In 'filter' problem with value '"+filter_str+"'");
	
	auto trans_filt = sp.set_trans_filt(cl_sel,trans,LOWER_UPPER_BOUND,errmsg);
	if(errmsg != "") error("In 'trans' problem with value '"+trans+"'");
	
	tab.heading.push_back("ID");
	tab.heading.push_back("t");
		
	const auto &ssp = state.species[p];
	
	auto tmin_div = model.calc_tdiv(tmin);
	auto tmax_div = model.calc_tdiv(tmax);
			 
	for(const auto &ind : ssp.individual){
		auto c = UNSET;
		for(const auto &ev : ind.ev){
			if(ev.type == NM_TRANS_EV || ev.type == M_TRANS_EV){
				auto tdiv = ev.tdiv;
				if(tdiv > tmin_div && tdiv < tmax_div){
					const auto &tra = sp.tra_gl[ev.tr_gl];
					auto tr = tra.tr;
					auto prob_str = trans_filt[tr];
					
					string err;
					auto prob = model.calculate_equation_zero_one(prob_str,tdiv,err);
					if(err != ""){
						error("An error was generated for filter '"+prob_str+"': "+err);
					}
					
					if(prob > 0){
						auto cc = c; if(cc == UNSET) cc = ev.c_after;
						
						const auto &cgl = sp.comp_gl[cc];
						
						// Implements filter
						for(auto cl = 0u; cl < sp.ncla; cl++){
							if(cl != cl_sel){
								auto ci = cgl.cla_comp[cl];
								auto probc_str = cf.cla[cl].comp_prob_str[ci];
								
								string err;
								auto probc = model.calculate_equation_zero_one(probc_str,tdiv,err);
								if(err != ""){
									error("An error was generated for filter '"+probc_str+"': "+err);
								}
								prob *= probc;								
							}
						}
						
						if(ran() < prob){
							vector <string> row;
							row.push_back(ind.name);
							row.push_back(tstr(model.calc_t(tdiv)));
							tab.ele.push_back(row);
						}
					}
				}
				
				c = ev.c_after;
			}
		}
	}
}


/// Simulates test data
void DataSim::test_data_sim(unsigned int p, const State &state, Table &tab)
{
	const auto &sp = model.species[p];
	
	auto Se_str = get_tag_value("Se"); if(Se_str == "") cannot_find_tag();	
	auto Sp_str = get_tag_value("Sp"); if(Sp_str == "") cannot_find_tag(); 
		
	auto pos = get_tag_value("pos");
	if(pos == "") pos = "1";
			
	auto neg = get_tag_value("neg");
	if(neg == "") neg = "0";
			
	if(pos == neg){ 
		error("'pos' and 'neg' cannot both have the same value"); 
	}
					
	auto comp = get_tag_value("comp"); if(comp == "") cannot_find_tag(); 
		
	auto spl = split(comp,',');
			
	auto cl = model.get_cl_from_comp(spl[0],p);
	if(cl == UNSET){ 
		error("For 'comp' value '"+spl[0]+"' is not a compartment");
	}
			
	const auto &claa = model.species[p].cla[cl];
	vector <bool> diag_test_sens(claa.ncomp,false);
			
	for(auto c = 0u; c < claa.ncomp; c++){
		auto name = claa.comp[c].name;
		if(find_in(spl,name) != UNSET) diag_test_sens[c] = true;
	}
			
	auto fr = get_frac();
	
	auto tims = get_times();
	
	tab.heading.push_back("ID");
	tab.heading.push_back("t");
	tab.heading.push_back("Result");

	const auto &ssp = state.species[p];
	
	for(const auto &ind : ssp.individual){
		for(auto k = 0u; k < tims.size(); k++){
			if(ran() < fr){	
				auto t = tims[k];
				auto tdiv = model.calc_tdiv(t);
				
				auto c = UNSET;
				auto e = 0u; 
				while(e < ind.ev.size() && ind.ev[e].tdiv <= tdiv){
					c = ind.ev[e].c_after;
					e++;
				}
			
				if(c != UNSET){
					const auto &cgl = sp.comp_gl[c];
					auto ci = cgl.cla_comp[cl];

					string res;
					
					if(diag_test_sens[ci]){
						string err;
						auto Se = model.calculate_equation(Se_str,tdiv,err);
						if(err != "") error("Problem with sensitivity '"+Se_str+"': "+err);
							
						if(ran() < Se) res = pos; else res = neg;
					}
					else{
						string err;
						auto Sp = model.calculate_equation(Sp_str,tdiv,err);
						if(err != "") error("Problem with specificity '"+Sp_str+"': "+err);
						
						if(ran() < Sp) res = neg; else res = pos;
					}
					
					vector <string> row;
					row.push_back(ind.name);
					row.push_back(tstr(model.calc_t(tdiv)));
					row.push_back(res);
					tab.ele.push_back(row);
				}
			}
		}
	}
}


/// Simulates population data
void DataSim::pop_data_sim(unsigned int p, const State &state, Table &tab)
{
	const auto &sp = model.species[p];
	
	auto filter_str = get_tag_value("filter");
	
	string errmsg;
	auto cf = sp.set_comp_filt_raw(filter_str,UNSET,LOWER_BOUND,errmsg);
	if(errmsg != "") error("In 'filter' problem with value '"+filter_str+"'");
	
	auto tims = get_times();
	
	auto om = load_obs_model();
	
	auto strat = get_stratify(cf,p);
	
	tab.heading.push_back("t");
	for(auto i = 0u; i < strat.cl_list.size(); i++){
		tab.heading.push_back(sp.cla[strat.cl_list[i]].name);
	}
	tab.heading.push_back("Population");
		
	auto &ssp = state.species[p];
	
	for(auto com = 0u; com < strat.comb.size(); com++){
		set_filt_strat(com,cf,strat,p);
			
		for(auto k = 0u; k < tims.size(); k++){
			auto t = tims[k];
			auto tdiv = model.calc_tdiv(t);
			auto ti = get_ti(tdiv);
			if(ti == state.T) ti = state.T-1;
			
			auto posg = get_comp_prob(cf,tdiv,UNSET,p);
		
			auto N = sp.comp_gl.size();
			vector <double> popc(N,0);
			
			switch(sp.type){
			case INDIVIDUAL:
				{
					for(const auto &ind : ssp.individual){
						auto c = UNSET;
						auto e = 0u; 
						while(e < ind.ev.size() && ind.ev[e].tdiv <= tdiv){
							c = ind.ev[e].c_after;
							e++;
						}
						
						if(c != UNSET) popc[c]++;
					}
				}
				break;
				
			case POPULATION: case DETERMINISTIC:
				for(auto c = 0u; c < N; c++){
					popc[c] = ssp.cpop_st[ti][c];
				}
				break;
			}
	
			vector <string> row;
			row.push_back(tstr(t));
			
			for(auto i = 0u; i < strat.cl_list.size(); i++){
				auto cl = strat.cl_list[i];
				row.push_back(sp.cla[cl].comp[strat.comb[com][i]].name);
			}
			
			auto sum = 0.0;
			for(auto c = 0u; c < N; c++){
				sum += posg[c]*popc[c];
			}
			
			auto result = obs_model_sample(sum,om);
			
			row.push_back(precision(result,5));

			tab.ele.push_back(row);
		}
	}
}


/// Simulates population transition data
void DataSim::pop_trans_data_sim(unsigned int p, const State &state, Table &tab)
{
	const auto &sp = model.species[p];
	
	auto filter_str = get_tag_value("filter");
	
	string errmsg;
	auto cf = sp.set_comp_filt_raw(filter_str,UNSET,LOWER_BOUND,errmsg);
	if(errmsg != "") error("In 'filter' problem with value '"+filter_str+"'");
	
	auto tims = get_times();
	
	auto om = load_obs_model();
	
	auto trans = get_tag_value("trans"); if(trans == "") cannot_find_tag(); 
		
	auto cl_sel = model.get_cl_from_trans(trans,p);
	if(cl_sel == UNSET){
		error("For 'trans' the value '"+trans+"' is not recognised");
		return;
	}
	
	auto strat = get_stratify(cf,p);
	
	if(find_in(strat.cl_list,cl_sel) != UNSET){
		error("'filter'	cannot specify transition classification as 'file'");
	}
	
	auto trans_filt = sp.set_trans_filt(cl_sel,trans,LOWER_UPPER_BOUND,errmsg);
	if(errmsg != "") error("In 'trans' problem with value '"+trans+"'");
	
	auto ndiv = tims.size()-1;
		
	const auto &tp = model.timepoint;

	auto T = state.T;

	vector <unsigned int> div_ref(T);

	auto k = 0u;
	for(auto j = 0u; j < ndiv; j++){
		auto t = tims[j+1];
		while(k < T && tp[k] < t-TINY){ div_ref[k] = j; k++;}		
		
		if(dif(t,tp[k],TINY)) error("Times must coincide with time points used to simulate.");
	}
	
	tab.heading.push_back("Start");
	tab.heading.push_back("End");
	for(auto i = 0u; i < strat.cl_list.size(); i++){
		tab.heading.push_back(sp.cla[strat.cl_list[i]].name);
	}
	tab.heading.push_back("Number");
	
	const auto &ssp = state.species[p];
	
	for(auto com = 0u; com < strat.comb.size(); com++){
		set_filt_strat(com,cf,strat,p);
		
		vector < vector <double> > res;
		res.resize(ndiv);
		for(auto d = 0u; d < ndiv; d++){
			res[d].resize(sp.tra_gl.size(),0);
		}

		switch(sp.type){
		case INDIVIDUAL:
			{
				for(const auto &ind : ssp.individual){
					for(auto e = 0u; e < ind.ev.size(); e++){
						const auto &ev = ind.ev[e];
						if(ev.type == NM_TRANS_EV || ev.type == M_TRANS_EV){
							auto ti = get_ti(ev.tdiv);
							res[div_ref[ti]][ev.tr_gl]++;
						}
					}
				}
			}
			break;
			
		case POPULATION: case DETERMINISTIC:
			{
				for(auto ti = 0u; ti < T; ti++){
					auto d = div_ref[ti];
					for(auto trg = 0u; trg < sp.tra_gl.size(); trg++){
						res[d][trg] += ssp.trans_num[trg][ti];
					}
				}
			}
			break;
		}
		
		for(auto d = 0u; d < ndiv; d++){
			vector <string> row;
			row.push_back(tstr(tims[d]));
			row.push_back(tstr(tims[d+1]));

			for(auto i = 0u; i < strat.cl_list.size(); i++){
				auto cl = strat.cl_list[i];
				row.push_back(sp.cla[cl].comp[strat.comb[com][i]].name);
			}
			
			auto t = (tims[d]+tims[d+1])/2;
			auto tdiv = model.calc_tdiv(t);

			auto posg = get_comp_prob(cf,tdiv,cl_sel,p); 

			auto tr_prob = get_tr_prob(trans_filt,tdiv);

			auto sum = 0.0;
			for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
				auto re = res[d][tr];
				if(re != 0){
					const auto &tra = sp.tra_gl[tr];
					auto c = tra.i; if(c == UNSET) c = tra.f;
					sum += tr_prob[tra.tr]*posg[c]*re;		
				}
			}
		
			auto result = obs_model_sample(sum,om);
			row.push_back(precision(result,5));
			tab.ele.push_back(row);
		}
	}
}
	
	
/// Simulates individual effect data
void DataSim::ind_effect_data_sim(unsigned int p, const State &state, Table &tab)
{
	const auto &sp = model.species[p];
	
	auto ie = get_tag_value("ie"); if(ie == "") cannot_find_tag(); 

	auto e = 0u; 
	while(e < sp.ind_effect.size() && sp.ind_effect[e].name != ie) e++;
	if(e == sp.ind_effect.size()) error("In 'ie' cannot find individual effect '"+ie+"'");
	
	const auto &ssp = state.species[p];
	
	tab.heading.push_back("ID");
	tab.heading.push_back("Value");
	for(const auto &ind : ssp.individual){
		vector <string> row;
		row.push_back(ind.name); 
		row.push_back(tstr(ind.ie[e]));
		tab.ele.push_back(row);
	}
}


/// Simulates individual effect data
void DataSim::ind_group_data_sim(unsigned int p, const State &state, Table &tab)
{
}


/// Simulates genetic data
void DataSim::genetic_data_sim(unsigned int p, const State &state, Table &tab)
{
	error("Simulation of genetic data hasn't been implemneted in BICI yet");
	auto type = get_tag_value("type"); if(type == ""){ cannot_find_tag(); return;}
			
			/*
			ds.gen_data_type = GenDataType(option_error("type",type,{"snp","matrix"},{ SNP_DATA, MATRIX_DATA}));
			if(ds.gen_data_type == UNSET) return;
			
			auto root = get_tag_value("root");

			if(type == "snp"){
				if(root == ""){ cannot_find_tag(); return;}
			}
			ds.SNP_root = root;
			
			auto mut_rate = get_tag_value("mut-rate"); if(mut_rate == ""){ cannot_find_tag(); return;}
			ds.mut_rate_str = mut_rate;
			
			auto seq_var = get_tag_value("seq-var"); if(seq_var == ""){ cannot_find_tag(); return;}
			ds.seq_var_str = seq_var;
			*/
}

/// Processes commands from the input file 

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sys/stat.h>
#include <iomanip>  

using namespace std;

#include "input.hh"
#include "utils.hh"
#include "matrix.hh"

/// Imports a data-table 
void Input::import_data_table_command(Command cname)
{
	if(model.mode == PPC){                       // If doing PPC then ignore data
		switch(cname){
		case COMP_DATA: case TRANS_DATA: case TEST_DATA: case POP_DATA:
		case POP_TRANS_DATA: case IND_EFFECT_DATA: case IND_GROUP_DATA: case GENETIC_DATA:
			return;
		default: break;
		}
	}
	
	switch(cname){
	case ADD_POP_SIM: case ADD_POP_POST_SIM: cname = ADD_POP; break;
	case REMOVE_POP_SIM: case REMOVE_POP_POST_SIM: cname = REMOVE_POP; break;
	case ADD_IND_SIM: case ADD_IND_POST_SIM: cname = ADD_IND; break;
	case REMOVE_IND_SIM: case REMOVE_IND_POST_SIM: cname = REMOVE_IND; break;
	case MOVE_IND_SIM: case MOVE_IND_POST_SIM: cname = MOVE_IND; break;
	case INIT_POP_SIM: cname = INIT_POP; break;
	default: break;
	}
	
	auto file = get_tag_value("file"); if(file == ""){ cannot_find_tag(); return;}

	auto tab = load_table(file); if(tab.error == true) return;

	auto cols = get_tag_value("cols");
	
	auto p = p_current;
	if(p == UNSET){ 
		alert_import("To load the data file the species must be set"); 
		return;
	}
	
	DataSource ds;
	ds.line_num = line_num;
	ds.p = p; ds.cl = UNSET;
	ds.cname = cname; ds.focal_cl = UNSET;
	ds.time_range = ALL_TIME; ds.time_start = UNSET; ds.time_end = UNSET;
	ds.focal_cl = UNSET;

	auto &om = ds.obs_model;

	switch(cname){  // Modifies specification dependent on data source
	case INIT_POP:
		{
			auto focal = get_tag_value("focal");
		
			if(focal != ""){
				ds.focal_cl = find_cl(p,focal); 
				if(ds.focal_cl == UNSET){
					alert_import("The focal classification '"+focal+"' not recognised");
					return;
				}
			}
			
			auto type = get_tag_value("type");
			ds.init_pop_type = InitPopType(option_error("type",type,{"fixed","dist"},{ INIT_POP_FIXED, INIT_POP_DIST}));
			if(ds.init_pop_type == UNSET) return;
			
			if(ds.init_pop_type == INIT_POP_DIST){
				if(ds.focal_cl == UNSET){
					auto prior = get_tag_value("prior"); 
					auto pri = convert_text_to_prior(prior,line_num,false);
					if(pri.error != ""){
						alert_import("For 'prior' error with expression '"+prior+"': "+pri.error);
						return;
					}
					ds.pop_prior = pri;
				}
			}
		}	
		break;
	
	case ADD_POP: case REMOVE_POP:
		break;
		
	case ADD_IND:
		break;
		
	case REMOVE_IND:
		break;
		
	case MOVE_IND:
		{
			auto name = get_tag_value("class"); if(name == ""){ cannot_find_tag(); return;}
			
			ds.cl = find_cl(p,name);
			if(ds.cl == UNSET){ 
				alert_import("For 'class' the value '"+name+"' is not a classification");
				return;
			}
		}
		break;
	
	case COMP_DATA:
		{
			auto name = get_tag_value("class"); if(name == ""){ cannot_find_tag(); return;}
			
			ds.cl = find_cl(p,name);
		
			if(ds.cl == UNSET){
				alert_import("For 'class' the value '"+name+"' is not a classification");
				return;
			}
		}	
		break;
	
	case TEST_DATA:
		{
			auto Se = get_tag_value("Se"); if(Se == ""){ cannot_find_tag(); return;}
			om.Se_str = Se;
			
			auto Sp = get_tag_value("Sp"); if(Sp == ""){ cannot_find_tag(); return;}
			om.Sp_str = Sp;
			
			auto pos = get_tag_value("pos");
			if(pos == "") pos = "1";
			
			auto neg = get_tag_value("neg");
			if(neg == "") pos = "0";
			
			if(pos == neg){ 
				alert_import("'pos' and 'neg' cannot both have the same value"); 
				return;
			}
				
			om.diag_pos = pos; om.diag_neg = neg;
				
			auto comp = get_tag_value("comp"); if(comp == ""){ cannot_find_tag(); return;}
		
			auto spl = split(comp,',');
			
			auto cl = get_cl_from_comp(spl[0],p);
			if(cl == UNSET){ 
				alert_import("Value '"+spl[0]+"' is not a compartment5");
				return;
			}
			
			om.diag_test_sens.cl = cl;
			ds.cl = cl;
			
			const auto &claa = model.species[p].cla[cl];
			om.diag_test_sens.comp.resize(claa.ncomp,false);
			
			for(auto c = 0u; c < claa.ncomp; c++){
				auto name = claa.comp[c].name;
				if(find_in(spl,name) != UNSET) om.diag_test_sens.comp[c] = true;
			}
		}
		break;
	
	case POP_DATA: 
		{	
			ds.filter_str = get_tag_value("filter");

			load_obs_model(om);
		}
		break;
	
	case TRANS_DATA:
		{
			ds.filter_str = get_tag_value("filter");
			
			auto name = get_tag_value("name");
			ds.filter_trans_str = name;
		
			auto cl_sel = get_cl_from_trans(name,p);
		
			if(cl_sel == UNSET){
				alert_import("For 'name' the value '"+name+"' is not recognised");
				return;
			}
		
			ds.cl = cl_sel;
			
			auto obsran = toLower(get_tag_value("obsrange"));
			if(obsran == ""){ cannot_find_tag(); return;}
	
			ds.time_range = TimeRange(option_error("obsrange",obsran,{"all","specify","file"},{ ALL_TIME, SPEC_TIME, FILE_TIME }));
			if(ds.time_range == UNSET) return;
	
			switch(ds.time_range){
			case ALL_TIME: case FILE_TIME: 
				break;
				
			case SPEC_TIME:
				auto start = get_tag_value("start"); if(start == ""){ cannot_find_tag(); return;} 
				if(!is_number(start,"start")) return;
				ds.time_start = number(start);
			
				auto end = get_tag_value("end"); if(end == ""){ cannot_find_tag(); return;} 
				if(!is_number(end,"end")) return;
				ds.time_end = number(end);
			
				if(ds.time_start >= ds.time_end){
					alert_import("The start time '"+start+"' must be before the end time '"+end+"'"); 
					return;
				}
				break;
			}
		}
		break;
		
	case POP_TRANS_DATA:
		{
			ds.filter_str = get_tag_value("filter");
			
			auto name = get_tag_value("name");
			ds.filter_trans_str = name;
			
			auto cl_sel = get_cl_from_trans(name,p);
		
			if(cl_sel == UNSET){
				alert_import("For 'name' the value '"+name+"' is not recognised"); return;
			}
		
			ds.cl = cl_sel;
		
			load_obs_model(om);
		}
		break;
		
	case GENETIC_DATA:
		{
			auto type = get_tag_value("type"); if(type == ""){ cannot_find_tag(); return;}
			
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
		}
		break;
		
	default: alert_import("Should not be default2"); return;
	}
	
	if(set_loadcol(cname,ds) == false) return;
	
	if(cname == GENETIC_DATA){
		switch(ds.gen_data_type){
		case SNP_DATA: set_SNP_columns(tab,ds); break;
		case MATRIX_DATA: set_genetic_matrix_columns(tab,ds); break;
		}
	}
	
	vector <string> col_name;
	for(auto c = 0u; c < ds.load_col.size(); c++){
		col_name.push_back(ds.load_col[c].heading);
	}
	
	if(cols != ""){
		auto spl = split(cols,',');

		if(spl.size() != col_name.size()){
			alert_import("'cols' does not have the correct number of entries (expected something in the order '"+stringify(col_name)+"')"); 
			return;
		}
		
		for(auto i = 0u; i < spl.size(); i++){
			for(auto j = 0u; j < col_name.size(); j++){
				if(spl[i] == col_name[j] && i != j){
					alert_import("'cols' does not have the correct order (expected something in the order '"+stringify(col_name)+"')"); 
					return;
				}
			}
		}
		
		col_name = spl;
	}
	
	ds.table = get_subtable(tab,col_name); if(ds.table.error == true) return;
	
	if(cname == TRANS_DATA){ // Removes "no" entries
		auto &tab = ds.table;
		auto r = 0u;
		while(r < tab.nrow){
			if(tab.ele[r][1] == "no"){
				tab.ele.erase(tab.ele.begin()+r);
				tab.nrow--;
			}
			else r++;
		}
	}
	
	if(false){
		print_table(ds.table);
		for(auto col : ds.load_col){
			cout << col.heading << " " << col.cl <<" " << col.type << " Data columns" << endl;
		}
	}
	
	data_source_check_error(ds);
	
	ds.index = model.species[p].source.size();
	
	model.species[p].source.push_back(ds);
}


/// Adds a species to the model
void Input::species_command(unsigned int loop)
{
	auto name = get_tag_value("name");
	if(name == ""){ cannot_find_tag(); return;}
	
	auto type = toLower(get_tag_value("type"));
	if(type == ""){ cannot_find_tag(); return;}
	
	auto sp_type = SpeciesType(option_error("type",type,{"population","individual"},{ POPULATION, INDIVIDUAL}));
	if(sp_type == UNSET) return;

	auto trans_tree = false;
	if(sp_type == INDIVIDUAL){
		auto trans_tree_str = toLower(get_tag_value("trans-tree"));
		if(trans_tree_str != ""){
			if(trans_tree_str == "on") trans_tree = true;
			else{
				if(trans_tree_str != "off"){ 
					alert_import("'trans-tree' must be either 'off' or 'on'");
					return;
				}
			}
		}
	}

	if(loop == 2){
		auto nsp = model.nspecies;
		add_species(name,sp_type,trans_tree);
		if(model.nspecies != nsp) p_current = model.nspecies-1;		
	}
	else{
		p_current = find_p(name);
		if(p_current == UNSET){
			alert_import("Problem finding species '"+name+"'"); 
			return;
		}
	}
	cl_current = UNSET;
}



/// Adds a classification to the model
void Input::classification_command(unsigned int loop)
{
	auto p = p_current;
	if(p == UNSET){ 
		alert_import("A species must be defined before a classification can be added",true); 
		return;
	}
		
	auto &sp = model.species[p];
	
	auto name = get_tag_value("name");
	if(name == ""){ cannot_find_tag(true); return;}
	
	auto clone = get_tag_value("clone");
	if(clone != ""){ // Clones a classification from another species
		if(loop == 2){
			auto p2 = find_p(clone);
			if(p2 == UNSET){
				alert_import("For 'clone' cannot find the species '"+clone+"'");
				return;
			}
			
			auto cl2 = find_cl(p2,name);
			
			if(cl2 == UNSET){
				alert_import("In clone species '"+clone+"' cannot find the classification '"+name+"'");
				return;
			} 
				
			auto ncl = sp.ncla;
			clone_class(p,p2,cl2);
			if(ncl != sp.ncla) cl_current = sp.ncla-1;		
		}
		else{
			cl_current = find_cl(p,name);
			if(cl_current == UNSET){
				alert_import("Problem finding classification '"+name+"'"); 
				return;
			}
		}
	}
	else{
		auto index = get_tag_value("index");
		if(index == ""){ cannot_find_tag(true); return;}
		if(index.length() != 1){ 
			alert_import("The index '"+index+"' must be just a single character",true); 
			return;
		}
		
		if(index == "t"){ 
			alert_import("The index 't' cannot be used because it is reserved for time variation",true); 
			return;
		}
		
		if(find_in(alphabet,index) == UNSET){ 
			alert_import("Index '"+index+"' must be from the lowercase alphabet",true); 
			return;
		}
		
		auto coord = toLower(get_tag_value("coord"));
		
		auto coord_type = CARTESIAN;
		
		if(coord != ""){
			coord_type = Coord(option_error("coord",coord,{"cartesian", "latlng"},{ CARTESIAN, LATLNG}));
			if(coord_type == UNSET) return;
		}
		
		if(coord == "latlng"){
			auto def_map = get_tag_value("map");
			if(def_map != ""){
				if(def_map != "load"){ 
					alert_import("'map' must be set to 'load'");
					return;
				}
			}
		}
		
		if(loop == 2){
			auto ncl = sp.ncla;
			add_classification(p,name,index,coord_type);
			if(sp.ncla != ncl) cl_current = sp.ncla-1;		
		}
		else{
			cl_current = find_cl(p,name);
			if(cl_current == UNSET){ 
				alert_import("Problem finding classification '"+name+"'"); 
				return;
			}
		}
	}
}


/// Adds a multiplier to a parameter (used in PPC)
void Input::param_mult_command()
{
	auto full_name = get_tag_value("name"); if(full_name == ""){ cannot_find_tag(); return;}
			
	auto th = 0u; 
	while(th < model.param.size() && model.param[th].full_name != full_name) th++;
	if(th == model.param.size()){
		alert_import("Could not find parameter '"+full_name+"'");
		return;
	}
		
	auto &par_orig = model.param[th];

	par_orig.param_mult = model.param.size();

	Param par; 
	par.variety = CONST_PARAM;
	par.time_dep = true;
	par.spline_info.on = true;
	par.used = true;
	par.param_mult = UNSET;
	par.factor = true;
	par.line_num = line_num;
	par.cat_factor = false;
	par.cat_factor_weight_on = false;
	
	auto knot_times_str = get_tag_value("knot-times"); 
	if(knot_times_str == ""){ cannot_find_tag(); return;}
	
	vector <string> knot_times;
	set_spline(knot_times_str,"",knot_times,false,par);
	
	par.dep = par_orig.dep;

	if(par_orig.time_dep == false){
		par.name = "f~"+par_orig.name;
		par.full_name = "f~"+par_orig.full_name+"(t)";
	}
	else{
		par.name = "f~"+par_orig.name;
		par.full_name = "f~"+par_orig.full_name;
	}
	auto pp = get_param_prop(par.full_name);
		
	auto mult = get_dependency(par.dep,pp,knot_times); if(mult == UNSET) return; 
	
	par.element_ref.resize(mult,UNSET);
	
	par.N = mult;
	par.trace_output = false;
	
	auto file = get_tag_value("constant"); if(file == ""){ cannot_find_tag(); return;}
	load_param_value(pp,file,par,"In 'file'");

	model.param.push_back(par);
};


/// Sets the species/classification which is being worked on
void Input::set_command()
{
	auto p = p_current;
	auto sp = get_tag_value("species");
	if(sp != ""){
		p_current = find_p(sp);
		if(p_current == UNSET){
			alert_import("Cannot set the species '"+sp+"' as this does not exist"); 
			return;
		}
	}
	
	auto name = get_tag_value("classification");
	if(name != ""){
		if(p == UNSET){
			alert_import("A species must be set before a 'classification' can be set"); 
			return;
		}
	
		cl_current = find_cl(p,name);
		if(cl_current == UNSET){
			alert_import("Cannot set the classification '"+name+"'"); 
			return;
		}
	}
}


/// Sets a camera in a classification (not used in c++ version)
void Input::camera_command()
{
	const auto &claa = get_claa();
	
	auto grid = toLower(get_tag_value("grid"));
	if(grid != ""){
		if(grid != "on" && grid != "off"){
			alert_import("For 'grid' the value '"+grid+"' must be 'on' or 'off'"); 
		}
	}
	
	auto sc = toLower(get_tag_value("comp-scale"));
	if(sc != ""){
		is_positive(sc,"sc");	
		auto num = number(sc);
		if(num < 0.1 || num > 10){
			alert_import("For 'comp-scale' the value '"+sc+"' must be between 0.1 and 10"); 
		}
	}
	
	auto scale = get_tag_value("scale");
	if(scale != "") is_positive(scale,"scale");
	
	auto x = get_tag_value("x");

	switch(claa.coord){
	case CARTESIAN:
		{
			auto x = get_tag_value("x"); if(x == "") cannot_find_tag();
			auto y = get_tag_value("y"); if(y == "") cannot_find_tag();
		
			is_number(x,"x");	
			is_number(y,"y");
		}
		break;
		
	case LATLNG:
		{
			auto lat_str = get_tag_value("lat"); if(lat_str == "") cannot_find_tag();
			auto lng_str = get_tag_value("lng"); if(lng_str == "") cannot_find_tag();
			
			if(!is_number(lat_str,"lat")) return;
			auto lat = number(lat_str);
			
			if(!is_number(lng_str,"lng")) return;
			auto lng = number(lng_str);
				
			check_latlng_error(lat,lng);
		}
		break;
		
	default:
		alert_import("'coord' must be 'cartesian' or 'latlng'"); 
		break;
	}
}


/// Gets a list of tags from a file 
vector < vector <Tag> > Input::get_tags_list(string file) 
{
	vector < vector <Tag> > tags_list;
	
	auto tab = load_table(file);
	
	auto N = tab.ncol;
	
	const auto &head = tab.heading;
	
	for(auto r = 0u; r < tab.nrow; r++){
		const auto &row = tab.ele[r];
		auto tags = cline_store.tags;
		for(auto i = 0u; i < N; i++){
			auto te = row[i];
			if(te != ""){
				Tag tag; tag.name = head[i]; tag.value = te; tag.done = 0;
				tags.push_back(tag);
			}
		}
		
		tags_list.push_back(tags);	

		if(false){
			for(auto ta:tags){
				cout << ta.name << " " << ta.value << "  "; 
			}
			cout << endl;
		}
	}

	return tags_list;
}


/// Checks all tags used 
void Input::check_tags_used(unsigned int r, const vector <Tag> &tags)
{
	for(auto i = 0u; i < tags.size(); i++){
		const auto &tag = tags[i];
		if(tag.done == 0){
			alert_import("The value '"+tag.value+"' for property '"+tag.name+"' is not used (line "+to_string(r+2)+" in file).");
			return;
		}			
	}
}


/// Adds a compartment to the model
void Input::compartment_command()
{
	compartment_command2(cline_store.tags);
}


/// Adds all compartments to the model
void Input::compartment_all_command()
{
	auto file = get_tag_value("file"); if(file == ""){ cannot_find_tag(true); return;}

	auto tags_list = get_tags_list(file);
	
	for(auto i = 0u; i < tags_list.size(); i++){
		auto &tags = tags_list[i];
		all_row = i;
		compartment_command2(tags);
		check_tags_used(i,tags);
	}
	
	for(auto &ta : cline_store.tags) ta.done = 1;
}


/// Adds a compartment to the model
void Input::compartment_command2(vector <Tag> &tags)
{
	if(check_claa_error() == true) return;
	
	const auto &claa = get_claa();
	auto p = p_current;
	auto cl = cl_current;
	
	auto name = get_tag_val("name",tags); if(name == ""){ cannot_find_tag(true); return;}
	 
	// Checks that a compartment name is valid
	check_name_input(name,"Compartment name");
	
	auto color = get_tag_val("color",tags); 
	if(color != ""){
		if(!is_Color(color)) alert_import("'"+color+"' is not a valid color.",true); 
	}
	
	auto fix_str = get_tag_val("fix",tags);
	
	double x = UNSET, y = UNSET, lat = UNSET, lng = UNSET;
	
	switch(claa.coord){
	case CARTESIAN:
		{
			auto x_str = get_tag_val("x",tags);
			auto y_str = get_tag_val("y",tags);
			
			if(x_str == "" && y_str == ""){
			}
			else{ 
				x_str = get_tag_val("x",tags);
				if(x_str == ""){ cannot_find_tag(true); return;}
				if(!is_number(x_str,"x",true)) return;
				x = number(x_str); 
			
				y_str = get_tag_val("y",tags);
				if(y_str == ""){ cannot_find_tag(true); return;}
				if(!is_number(y_str,"y",true)) return;
				y = number(y_str);
				
				if(x != UNSET && y != UNSET){
					for(auto i = 0u; i < claa.ncomp; i++){
						if(claa.comp[i].x == x && claa.comp[i].y == y){		
							alert_import("Compartments '"+name+"' and '"+claa.comp[i].name+"' cannot have the same position",true);	
							return;
						}
					}
				}
			}
		}
		break;
	
	case LATLNG:
		{
			auto bound_file = get_tag_val("boundary",tags);
			if(bound_file != ""){
				auto i = import_geojson(bound_file); if(i == UNSET) return;
				
				auto po = boundary_mean_latlng(i,name);
			
				lat = po.lat; lng = po.lng;
			}
			else{
				auto lat_str = get_tag_val("lat",tags);
				if(!is_number(lat_str,"lat",true)) return;
				lat = number(lat_str);
				
				auto lng_str = get_tag_val("lng",tags);
				if(!is_number(lng_str,"lng",true)) return;
				lng = number(lng_str);
			
				if(lat != UNSET && lng != UNSET){
					if(check_latlng_error(lat,lng,true) == true) return;
				}					
			}
		}
		break;
	}

	CompInfected infected = COMP_INFECTED_UNSET;
	auto infected_str = get_tag_val("infected",tags); 
	if(infected_str != ""){
		if(infected_str == "true") infected = COMP_INFECTED;
		else{
			if(infected_str == "false") infected = COMP_UNINFECTED;
			else{
				alert_import("'infected' must be 'true' or 'false'"); 
				return;
			}
		}
	}
	
	auto markov_branch = false;
	auto markov_branch_str = get_tag_val("branch-prob",tags); 
	if(markov_branch_str != ""){
		if(markov_branch_str == "true") markov_branch = true;
		else{
			if(markov_branch_str == "false") markov_branch = false;
			else{
				alert_import("'branch-prob' must be 'true' or 'false'"); 
				return;
			}
		}
	}
	
	add_compartment(name,p,cl,x,y,lat,lng,markov_branch,infected,"");
}


/// Adds a transition / source / sink to the model
void Input::transition_command()
{
	transition_command2(cline_store.tags);
}


/// Adds a transition / source / sink to the model
void Input::transition_all_command()
{
	auto file = get_tag_value("file"); if(file == ""){ cannot_find_tag(); return;}

	auto tags_list = get_tags_list(file);
	
	for(auto i = 0u; i < tags_list.size(); i++){
		auto &tags = tags_list[i];
		all_row = i;
		transition_command2(tags);
		check_tags_used(i,tags);
	}
	
	for(auto &ta : cline_store.tags) ta.done = 1;
}


/// Adds a transition / source / sink to the model
void Input::transition_command2(vector <Tag> &tags)
{
	if(check_claa_error() == true) return;
	auto &claa = get_claa();

	auto p = p_current;
	auto cl = cl_current;
	
	auto te = get_tag_val("name",tags); if(te == ""){ cannot_find_tag(); return;}
	
	if(te.length() < 2){ 
		alert_import("Transition name '"+te+"' is not understood. Expected format: 'Initial compartment'->'final compartment'."); 
		return;
	}
	
	auto i = 0u; while(i < te.length()-1 && te.substr(i,2) != "->") i++;
	
	if(i == te.length()-1){ 
		alert_import("Transition name '"+te+"' is not understood. Expected format: 'Initial compartment'->'final compartment'."); 
		return;
	}
	
	string fr = trim(te.substr(0,i));
	if(fr == "+") fr = "Source";
	
	string to = trim(te.substr(i+2));
	if(to == "-") to = "Sink";
	
	auto value_str = get_tag_val("value",tags); if(value_str == ""){ cannot_find_tag(); return;}

	auto trans_def = extract_trans_def(value_str);

	if(trans_def.set == false){ 
		alert_import("In the expression '"+value_str+"' the transition distribution is incorrectly specified. Expected format: 'exp(rate:...)', 'gamma(mean:...,cv:...)', 'erlang(mean:...,shape:)', 'log-normal(mean:...,cv:...)', 'weibull(scale:...,shape:...)' or 'period(time:...)'.");
		return;
	}
	
	unsigned int ci, cf;
	
	if(fr == "Source") ci = UNSET;
	else{
		ci = find_c(p,cl,fr);
		if(ci == UNSET){ cannot_find_trans_comp(te,p,cl,fr); return;}
	}
	
	if(to == "Sink") cf = UNSET;
	else{
		cf = find_c(p,cl,to);
		if(cf == UNSET){
			cannot_find_trans_comp(te,p,cl,to);
			return;
		}
	}

	for(auto tr = 0u; tr < claa.ntra; tr++){
		if(claa.tra[tr].i == ci && claa.tra[tr].f == cf){
			if(ci == UNSET){ 
				alert_import("Source to '"+to+"' already exists"); 
				return;
			}
			else{
				if(cf == UNSET){
					alert_import("Sink from '"+fr+"' already exists"); 
					return;
				}
				else{
					alert_import("Transition from '"+fr+"' to '"+to+"' already exists");
					return;
				}
			}
		}
	}

	auto x = get_tag_val("x",tags);
	auto y = get_tag_val("y",tags);
				
	auto xmid = get_tag_val("mid-x",tags);
	auto ymid = get_tag_val("mid-y",tags);
			
	auto lng = get_tag_val("lng",tags);
	auto lat = get_tag_val("lat",tags);
	auto lngmid = get_tag_val("mid-lng",tags);
	auto latmid = get_tag_val("mid-lat",tags);
			
	add_transition(p,cl,ci,cf,trans_def.type);

	auto &tra = claa.tra[claa.ntra-1];
	
	string sup, sup2;
	if(ci == UNSET){ sup = claa.comp[cf].name; sup2 = sup;}
	else{
		if(cf == UNSET){ sup = "("+claa.comp[ci].name+"→Sink)"; sup2 = claa.comp[ci].name;}
		else{
			sup = "("+claa.comp[ci].name+"→"+claa.comp[cf].name+")"; sup2 = claa.comp[ci].name;
		}
	}
		
	auto bp = get_tag_val("bp",tags);
	
	if(bp == "" || bp == "*"){
		tra.bp_set = BP_UNSET;
		tra.bp.eq_ref = UNSET;
	}
	else{
		tra.bp_set = BP_SET;
		tra.bp = he(add_equation_info(bp,BP,p,cl));
	}
	
	switch(trans_def.type){
	case EXP_RATE:
		{	
			auto eqtype = TRANS_RATE; if(tra.variety == SOURCE_TRANS) eqtype = SOURCE_RATE;
			
			tra.dist_param.push_back(he(add_equation_info(trans_def.rate,eqtype,p,cl)));
		}
		break;
		
	case EXP_MEAN:
		{	
			auto eqtype = TRANS_MEAN; if(tra.variety == SOURCE_TRANS) eqtype = SOURCE_MEAN;
			
			tra.dist_param.push_back(he(add_equation_info(trans_def.mean,eqtype,p,cl)));
		}
		break;
		
	case EXP_RATE_NM:
		{	
			auto eqtype = TRANS_NM_RATE; 
			
			tra.dist_param.push_back(he(add_equation_info(trans_def.rate,eqtype,p,cl)));
		}
		break;
		
	case EXP_MEAN_NM:
		{	
			auto eqtype = TRANS_NM_MEAN; 
			
			tra.dist_param.push_back(he(add_equation_info(trans_def.mean,eqtype,p,cl)));
		}
		break;
	
	case GAMMA:
		{
			tra.dist_param.push_back(he(add_equation_info(trans_def.mean,TRANS_MEAN,p,cl)));
			tra.dist_param.push_back(he(add_equation_info(trans_def.cv,TRANS_CV,p,cl)));
		}
		break;
		
	case ERLANG:
		{
			tra.dist_param.push_back(he(add_equation_info(trans_def.mean,TRANS_MEAN,p,cl)));
			
			
			auto num = number(trans_def.shape);
			
			if(num == UNSET || num <= 0 || int(num) != num){
				alert_import("For an Erlang distribution the shape parameter must be a positive integer");
			}
			
			tra.dist_param.push_back(he(add_equation_info(trans_def.shape,TRANS_SHAPE,p,cl)));
		}
		break;
		
	case LOG_NORMAL:
		{
			tra.dist_param.push_back(he(add_equation_info(trans_def.mean,TRANS_MEAN,p,cl)));
			tra.dist_param.push_back(he(add_equation_info(trans_def.cv,TRANS_CV,p,cl)));
		}
		break;
		
	case WEIBULL:
		{
			tra.dist_param.push_back(he(add_equation_info(trans_def.scale,TRANS_SCALE,p,cl)));
			tra.dist_param.push_back(he(add_equation_info(trans_def.shape,TRANS_SHAPE,p,cl)));
		}
		break;
	
	case PERIOD:
		{
			auto period = number(trans_def.time);
			if(period == UNSET || period <= 0){
				alert_import("The period must be a positive number");		
			}
		
			period += NEAR_DIV_THRESH*model.details.dt*2;
			stringstream ss;
			ss << fixed << std::setprecision(16);
			ss << period;
			 
			tra.dist_param.push_back(he(add_equation_info(ss.str(),TRANS_MEAN,p,cl)));
		}
		break;
	}
}


/// Sets the data directory
void Input::datadir_command()
{
	datadir = get_tag_value("folder"); if(datadir == ""){ cannot_find_tag(); return;}
	
	datadir = get_data_dir(datadir);
	
	if(check_char_allowed(datadir,"<>\"|?*") == false) return;

	struct stat st;
	if (stat(datadir.c_str(), &st) == -1){  
		alert_import("The data folder '"+datadir+"' does not exist");
		return;
	}
}


/// Sets a description of the analysis (not used in c++ version)
void Input::description_command()
{
	auto te = get_tag_value("text"); if(te == ""){ cannot_find_tag(); return;}
}


/// Adds a label to the model (not used in c++ version)
void Input::label_command()
{
	auto te = get_tag_value("text"); if(te == ""){ cannot_find_tag(); return;}
	auto tesize = get_tag_value("text-size"); 
	auto x = get_tag_value("x"); 
	auto y = get_tag_value("y");
	auto lat = get_tag_value("lat"); 
	auto lng = get_tag_value("lng"); 
	auto color = get_tag_value("color"); 
}


/// Adds a box around compartments in the model (not used in c++ version)
void Input::box_command()
{
	auto te = get_tag_value("text");
	
	auto tesize = get_tag_value("textsize"); 
	if(tesize != ""){ 
		if(!is_number(tesize,"textsize")) return;
	}
	
	auto comps = get_tag_value("comps"); if(comps == ""){ cannot_find_tag(); return;}
	auto color = get_tag_value("color"); 
}


/// Sets the value for a parameter in the model
void Input::param_command()
{
	auto full_name = get_tag_value("name"); if(full_name == ""){ cannot_find_tag(); return;}
	
	auto pp = get_param_prop(full_name);

	Param par; 
	par.variety = UNSET_PARAM;
	par.name = pp.name;
	par.full_name = full_name;
	par.time_dep = pp.time_dep;
	par.spline_info.on = false;
	par.used = false;
	par.param_mult = UNSET;
	par.factor = false;
	par.line_num = line_num;
	
	vector <string> knot_times;
	if(pp.time_dep == true){
		auto knot_times_str = get_tag_value("knot-times"); if(knot_times_str == ""){ cannot_find_tag(); return;}
		
		auto smooth = trim(toLower(get_tag_value("smooth")));
		
		set_spline(knot_times_str,smooth,knot_times,true,par);
	}
	
	auto cons = get_tag_value("constant"); 
	auto dist = get_tag_value("dist");
	auto dist_split = get_tag_value("dist-split");
	auto value = get_tag_value("value");
	auto reparam = get_tag_value("reparam"); 
	auto prior = get_tag_value("prior"); 
	auto prior_split = get_tag_value("prior-split"); 

	auto mode = model.mode;

	switch(mode){
	case INF: case PPC:
		value = "";
		break;
		
	case SIM: 
		prior = ""; prior_split = "";
		break;
	default: break;
	}
	
	vector <ParamTag> param_tag;
	ParamTag pt; 
	pt.val = cons; pt.tag = "constant"; param_tag.push_back(pt);
	pt.val = dist; pt.tag = "dist"; param_tag.push_back(pt);
	pt.val = dist_split; pt.tag = "dist-split"; param_tag.push_back(pt);
	pt.val = value; pt.tag = "value"; param_tag.push_back(pt);
	pt.val = reparam; pt.tag = "reparam"; param_tag.push_back(pt);
	pt.val = prior; pt.tag = "prior"; param_tag.push_back(pt);
	pt.val = prior_split; pt.tag = "prior-split"; param_tag.push_back(pt);
	
	for(auto j = 0u; j < param_tag.size(); j++){
		for(auto i = j+1; i < param_tag.size(); i++){
			if(param_tag[j].val != "" && param_tag[i].val != ""){
				alert_import("'"+param_tag[j].tag+"' and '"+param_tag[i].tag+"' cannot both be set"); 
				return;
			}
		}
	}
	
	par.not_set = false;
	auto j = 0u; while(j < param_tag.size() && param_tag[j].val == "") j++;
	if(j == param_tag.size()){
		par.not_set = true;
	}
	
	auto mult = get_dependency(par.dep,pp,knot_times); if(mult == UNSET) return; 
	
	par.element_ref.resize(mult,UNSET);
	
	par.N = mult;
	
	// Specifies a parameter as a categorical factor
	{
		par.cat_factor = false;
		par.cat_factor_weight_on = false;
		auto cat_factor = toLower(get_tag_value("factor")); 
		if(cat_factor != ""){
			if(cat_factor == "true"){
				par.cat_factor = true;
				if(par.time_dep){
					alert_import("'factor' cannot be set to 'true' for a time varying parameter");
				}
				
				if(reparam != ""){
					alert_import("'factor' cannot be set to 'true' if 'reparam' is set");
				}
				if(dist != ""){
					alert_import("'factor' cannot be set to 'true' if 'dist' is set");
				}
				if(prior_split != ""){
					alert_import("'factor' cannot be set to 'true' if 'prior-split' is set");
				}
				
				if(par.dep.size() == 0){
					alert_import("'factor' cannot be set to 'true' unless the parameter has some dependency");
				}								
			}
			else{
				if(cat_factor != "false"){
					alert_import("'factor' can only take the values 'true' or 'false'");
				}
			}
		}
	}
	
	if(par.cat_factor){
		//for(auto i = 0u; i < mult; i++) par.add_element(i);
		par.weight.resize(par.N,1);
		
		auto weight = get_tag_value("factor-weight"); 
		if(weight != ""){
			par.cat_factor_weight_on = true;
			if(is_file(weight) == false){
				alert_import("'factor-weight' must be a data table");
			}		
			else{
				load_weight_value(pp,weight,par,"In 'file'");
			}
		}
	}
	
	// Sets default value to zero
	//for(auto i = 0u; i < mult; i++) par.value[i].te = "0";
	
	par.trace_output = true;
	if(cons != "" || mult > model.details.param_output_max) par.trace_output = false;

	if(par.name == dist_matrix_name){
		alert_import("The distance matrix '"+par.full_name+"' must not be set");
	}
	
	auto pre = "Parameter '"+par.full_name+"': ";
	
	if(value == "auto") emsg("auto no longer supported");
		
	if(pp.dep.size() == 0){
		if(value != ""){
			auto val = number(value);
			if(val == UNSET){
				alert_import(pre+"The value '"+value+"' must be a number"); 
			}
			
			par.variety = CONST_PARAM;
			par.set_cons(0,val);
		}
		
		if(reparam != ""){
			par.variety = REPARAM_PARAM;
			par.reparam_eqn = reparam;
			//par.set_value_eqn(0,reparam);
		}
		
		if(cons != ""){
			auto val = number(cons);
			if(val == UNSET){
				alert_import(pre+"The constant '"+cons+"' must be a number"); 
			}
			
			par.variety = CONST_PARAM;
			par.set_cons(0,val);
		}
	}
	else{
		if(value != "" || cons != "" || reparam != ""){
			par.variety = CONST_PARAM;		
			
			string desc = pre+"For 'value'";
			auto valu = value; 
			if(valu == ""){
				if(cons != ""){
					valu = cons; desc = pre+"For 'const'";
				}
				else{
					if(reparam != ""){
						valu = reparam; desc = pre+"For 'reparam'";
						par.variety = REPARAM_PARAM;		
					}
					else{ 
						alert_import("Problem importing"); 
						return;
					}
				}
			}
			
			if(reparam != "" && is_file(valu) == false){
				par.reparam_eqn = reparam;
				//load_reparam_eqn(reparam,par);
			}
			else{
				if(is_file(valu) == false){
					double val = number(valu);
					
					if(par.variety != CONST_PARAM) emsg("Should be const");
				
					if(val == UNSET){
						alert_import(desc+" '"+valu+"' is not a number");
						return;
					}
						
					auto ref = par.add_cons(val);
					for(auto k = 0u; k < mult; k++) par.element_ref[k] = ref;
				}
				else{
					load_param_value(pp,valu,par,desc);
				}
			}
		}
	}

	if(prior != ""){
		par.variety = PRIOR_PARAM;
		
		auto pri = convert_text_to_prior(prior,line_num,false);
		if(pri.error != ""){
			alert_import("For 'prior' error with expression '"+prior+"': "+pri.error);
			return;
		}
		
		if(par.cat_factor){
			if(pri.type != MDIR_PR){
				alert_import("If 'factor' is set to 'true' the prior must be 'mdir(...)'");
				return;
			}
		}
		else{
			if(pri.type == MDIR_PR){
				alert_import("Prior 'mdir(...)' can only be used when 'factor' is set to 'true'");
				return;
			}
		}
		
		auto ref = model.prior.size();
		par.default_prior_ref = ref;
		model.prior.push_back(pri);
		for(auto i = 0u; i < mult; i++) par.set_prior(i,ref);
	}
	
	if(prior_split != ""){
		//par.prior.resize(mult);
		par.variety = PRIOR_PARAM;
		
		if(par.dep.size() == 0){
			alert_import("'prior-split' can only be used if the parameter has a dependency."); 
			return;
		}
	
		auto tab = load_table(prior_split);
		if(tab.error == true) return;	
		
		auto col_name = pp.dep;
		col_name.push_back("Prior");
		
		auto subtab = get_subtable(tab,col_name); if(subtab.error == true) return;
		
		auto ncol = subtab.ncol;
		
		for(auto r = 0u; r < subtab.nrow; r++){
			vector <unsigned int> ind(ncol-1);
			for(auto i = 0u; i < ncol-1; i++){
				ind[i] = par.dep[i].hash_list.find(subtab.ele[r][i]);
				if(ind[i] == UNSET){ 
					alert_import("The table element '"+subtab.ele[r][i]+"' is not valid (column '"+subtab.heading[i]+"', row "+tstr(r+2)+")");
					return;
				}
			}
			
			auto pri = convert_text_to_prior(subtab.ele[r][ncol-1],line_num,false);
			
			if(pri.error != ""){
				alert_import("The table element '"+subtab.ele[r][ncol-1]+"' is not a valid prior specification: "+pri.error+" (col '"+subtab.heading[ncol-1]+"', row "+tstr(r+2)+").");
				return;
			}
			
			set_prior_element(par,ind,pri);
		}
	}
	
	if(dist != ""){
		par.variety = DIST_PARAM;
		
		auto pri = convert_text_to_prior(dist,line_num,true);
		if(pri.error != ""){
			alert_import("For 'dist' error with expression '"+dist+"': "+pri.error);
			return;
		}
		
		par.default_prior_ref = model.prior.size();
		model.prior.push_back(pri);
	}

	if(dist_split != ""){
		par.variety = DIST_PARAM;
	
		if(par.dep.size() == 0){
			alert_import("'dist-split' can only be used if the parameter has a dependency."); 
			return;
		}
		
		auto tab = load_table(dist_split);
		if(tab.error == true) return;	
		
		auto col_name = pp.dep;
		col_name.push_back("Dist");
		
		auto subtab = get_subtable(tab,col_name); if(subtab.error == true) return;
		
		auto ncol = subtab.ncol;
		
		for(auto r = 0u; r < subtab.nrow; r++){
			vector <unsigned int> ind(ncol-1);
			for(auto i = 0u; i < ncol-1; i++){
				ind[i] = par.dep[i].hash_list.find(subtab.ele[r][i]);
	
				if(ind[i] == UNSET){ 
					alert_import("The table element '"+subtab.ele[r][i]+"' is not valid (column '"+subtab.heading[i]+"', row "+tstr(r+2)+")");
					return;
				}
			}
			
			auto pri = convert_text_to_prior(subtab.ele[r][ncol-1],line_num,true);
			
			if(pri.error != ""){
				alert_import("The table element '"+subtab.ele[r][ncol-1]+"' is not a valid distribution specification: "+pri.error+" (col '"+subtab.heading[ncol-1]+"', row "+tstr(r+2)+").");
				return;
			}
			
			set_prior_element(par,ind,pri);
		}
	}
			
	auto sim_sample = toLower(get_tag_value("sim-sample")); 
	if(sim_sample == "") sim_sample = "true";
	else{
		if(dist == ""){
			alert_import("'sim-sample' can only be set if a distribution is set through 'dist'");
			return;
		}
	}
	
	if(sim_sample == "false") par.sim_sample = false;
	else{
		if(sim_sample == "true") par.sim_sample = true;
		else{
			alert_import("'sim-sample' must be either 'true' or 'false'."); 
			return;
		}
	}
	
	for(const auto &par2 : model.param){
		if(par2.name == par.name){
			alert_import("Parameter already defined"); 
			return;
		}
	}
	
	if(par.variety == UNSET_PARAM){
		if(par.not_set){
			switch(model.mode){	
			case INF: case PPC:
				alert_import("A prior should be set for parameter '"+par.full_name+"'");
				break;
	
			case SIM: 
				alert_warning("A value has not been set for parameter '"+par.full_name+"'");
				break;
			default: break;
			}
		}
		else{
			alert_import("Parameter variety for '"+par.name+"' is unset");
		}
		return;
	}
	else{	
		if(par.cat_factor && model.mode == INF){ // Adds an extra parameter to model categorical factor
			add_param_cat_factor(par);
		}
		
		model.param.push_back(par);
	}
	
	if(par.variety == REPARAM_PARAM){
		for(const auto &ele : par.element) he(ele.value);
	}
}


/// Definition for derived quantities
void Input::derived_command()
{
	auto full_name = get_tag_value("name"); if(full_name == ""){ cannot_find_tag(); return;}
	auto eqn_name = get_tag_value("eqn"); if(eqn_name == ""){ cannot_find_tag(); return;}
	
	auto pp = get_param_prop(full_name);
	if(pp.time_dep == true) pp.dep.pop_back();
	
	if(pp.name == "D") alert_import("Name 'D' is reserved for the distance matrix");
	if(pp.name == "t") alert_import("Name 't' is reserved for time");
	
	Derive der;
	der.name = pp.name;
	der.time_dep = pp.time_dep;

	auto mult = get_dependency(der.dep,pp,vector <string> ()); if(mult == UNSET) return; 
	
	const auto &depend = der.dep;
	
	auto der_eqn_raw = he(add_equation_info(eqn_name,DERIVE_EQN));
	
	vector <DepConv> dep_conv;
	for(auto d = 0u; d < depend.size(); d++){
		const auto &dep = depend[d];
		DepConv dc; 
		dc.before = dep.index_with_prime;
		dep_conv.push_back(dc);
	}
	
	auto swap_temp = swap_template(der_eqn_raw.te,dep_conv);
	if(swap_temp.warn != ""){ 
		alert_import(swap_temp.warn); 
		return;
	}
	
	for(auto i = 0u; i < mult; i++){
		for(auto d = 0u; d < depend.size(); d++){
			const auto &dep = depend[d];
			dep_conv[d].after = dep.list[(i/dep.mult)%dep.list.size()];
		}
		
		auto der_eqn = der_eqn_raw;
		
		der_eqn.te = swap_index_temp(dep_conv,swap_temp);
		
		if(check_swap){
			auto te_ch = der_eqn_raw.te;
			auto res = swap_index(te_ch,dep_conv);
			if(res.warn != ""){ 
				alert_import(res.warn); 
				return;
			}
			
			if(te_ch != der_eqn.te){
				cout << der_eqn.te << " " << te_ch << " compare" << endl; 
				emsg_input("Swap index dif res");
			}
		}
		
		der.eq.push_back(der_eqn);
	}	
						
	model.derive.push_back(der);
}


/// Applies the simulation command
void Input::simulation_command()
{
	auto &details = model.details;
	
	auto start = get_tag_value("start"); 
	if(start == ""){ terminate = true; cannot_find_tag(); return;}
	if(!is_number(start,"start")){ terminate = true; return;}
	details.t_start = number(start);	
	
	auto end = get_tag_value("end");
	if(end == ""){ terminate = true; cannot_find_tag(); return;}

	if(!is_number(end,"end")){ terminate = true; return;}
	auto end_num = number(end);
	
	details.t_end = end_num;
	
	if(details.t_start >= details.t_end){
		terminate = true;
		alert_import("The start time '"+start+"' must be before the end time '"+end+"'");
		return;
	}

	details.number = check_pos_integer("number",SIM_NUM_DEFAULT);

	if(details.number%mpi.ncore != 0 && model.mode == SIM){
		alert_import("For 'number' the value '"+tstr(details.number)+"' must be a multiple of the number of cores");
	}		 
	details.num_per_core = details.number/mpi.ncore;

	details.nchain = 1;
	
	details.seed = get_seed();
	
	auto alg = toLower(get_tag_value("algorithm"));
	if(alg == "") alg = "gillespie";
	
	details.algorithm = Algorithm(option_error("algorithm",alg,{"gillespie","tau"},{ GILLESPIE, TAU }));
	if(details.algorithm == UNSET) return;
	
	auto dt_str = get_tag_value("timestep"); 
	if(dt_str == ""){ terminate = true; cannot_find_tag(); return;}
	
	if(!is_positive(dt_str,"timestep")){ terminate = true; return;}
	
	details.dt = number(dt_str);
	details.stochastic = true;
	details.individual_max = check_pos_integer("ind-max",INDMAX_DEFAULT);
	details.param_output_max = check_pos_integer("param-output-max",PARAM_OUTPUT_MAX_DEFAULT);
	details.anneal_type = ANNEAL_NONE;
	details.anneal_rate = UNSET;
	details.anneal_power = UNSET;
	details.diagnostics_on = false;
	
	check_dt(details);
}


/// Applies the inference command
void Input::inference_command()
{
	auto &details = model.details;

	auto start = get_tag_value("start"); 
	if(start == ""){ terminate = true; cannot_find_tag(); return;}
	if(!is_number(start,"start")){ terminate = true; return;}
	details.t_start = number(start);
	
	auto end = get_tag_value("end"); 
	if(end == ""){ terminate = true; cannot_find_tag(); return;}
	if(!is_number(end,"end")){ terminate = true; return;}
	details.t_end = number(end);
	
	if(details.t_start >= details.t_end){
		terminate = true;
		alert_import("The start time '"+start+"' must be before the end time '"+end+"'");
		return;
	}
	
	details.seed = get_seed();
	
	details.algorithm = DA_MCMC;
	
	auto alg = toUpper(get_tag_value("algorithm"));
	
	if(alg != ""){
		details.algorithm = Algorithm(option_error("algorithm",alg,{"DA-MCMC","PAS-MCMC","MFA","ABC","ABC-SMC","ABC-MBP","PMCMC","HMC"},{ DA_MCMC, PAS_MCMC, MFA_ALG, ABC_ALG, ABC_SMC_ALG, ABC_MBP, PMCMC, HMC }));
		if(details.algorithm == UNSET) return;
	}
	
	auto dt_str = get_tag_value("timestep"); 
	if(dt_str == ""){ terminate = true; cannot_find_tag(); return;}

	if(!is_positive(dt_str,"timestep")){ terminate = true; return;}
	
	details.dt = number(dt_str);
	details.sample = MCMC_SAMPLE_DEFAULT;
	details.output_param = MCMC_OP_PARAM_DEFAULT;
	details.output_state = MCMC_OP_STATE_DEFAULT;
	
	details.nchain = 1;
	details.gen_update = UNSET;
	
	details.individual_max = check_pos_integer("ind-max",INDMAX_DEFAULT);
	details.param_output_max = check_pos_integer("param-output-max",PARAM_OUTPUT_MAX_DEFAULT);

	auto algo = details.algorithm;
	
	switch(algo){
	case DA_MCMC:
		details.sample = check_pos_integer("update",MCMC_SAMPLE_DEFAULT);
		details.nchain = check_pos_integer("nchain");
		if(details.nchain%mpi.ncore != 0 && model.mode == INF){
			alert_import("For 'nchain' the value '"+tstr(details.nchain)+"' must be a multiple of the number of cores");
		}
		details.num_per_core = check_pos_integer("chain-per-core",MCMC_CHAIN_PER_CORE_DEFAULT);
		break;
		
	case ABC_ALG:
		{
			details.sample = check_pos_integer("sample",ABC_SAMPLE_DEFAULT);
			details.accfrac = check_zero_one("acc-frac",ABC_ACFRAC_DEFAULT);
			auto num = (unsigned int)(details.sample/details.accfrac);
			
			if(num%mpi.ncore != 0 && model.mode == INF){
				alert_import("'sample' divided by 'acc-frac' must be a multiple of the number of cores"); 
			}			
			details.num_per_core = num/mpi.ncore;
		}
		break;
		
	case ABC_SMC_ALG:
		details.sample = check_pos_integer("sample",ABC_SAMPLE_DEFAULT);
		details.accfrac = check_zero_one("acc-frac",ABCSMC_ACFRAC_DEFAULT);
		details.numgen = check_pos_integer("gen",ABCSMC_GEN_DEFAULT);
		details.kernelsize = check_pos_number("kernel-size",ABCSMC_KERNEL_DEFAULT);
		
		if(details.sample%mpi.ncore != 0 && model.mode == INF){
			alert_import("For 'sample' the value '"+tstr(details.sample)+"' must be a multiple of the number of cores");
		}	
		details.num_per_core = details.sample/mpi.ncore;
		break;
	
	case PAS_MCMC:
		details.sample = check_pos_integer("update",MCMC_SAMPLE_DEFAULT);
		details.nchain = check_pos_integer("npart");
		if(details.nchain <= 1){
			alert_import("For 'npart' the value '"+tstr(details.nchain)+"' must be 2 or above");
		}
		
		if(details.nchain%mpi.ncore != 0 && model.mode == INF){
			alert_import("For 'npart' the value '"+tstr(details.nchain)+"' must be a multiple of the number of cores");
		}
		details.num_per_core = check_pos_integer("part-per-core",PAS_PART_PER_CORE_DEFAULT);
		break;
		
	default: break;
	}
	
	details.diagnostics_on = true;
	if(algo == DA_MCMC || algo == PAS_MCMC){
		details.output_param = check_pos_integer("param-output",MCMC_OP_PARAM_DEFAULT);
		details.output_state = check_pos_integer("state-output",MCMC_OP_STATE_DEFAULT);
		details.output_param /= details.nchain;
		details.output_state /= details.nchain;	
		if(details.output_param < 1) details.output_param = 1;
		if(details.output_state < 1) details.output_state = 1;	
		
		auto ncore = details.nchain/details.num_per_core;
		if(ncore != mpi.ncore && model.mode == INF && mpi.core_spec_on == false){
			alert_import("The number of cores is '"+tstr(mpi.ncore)+"' and should be '"+tstr(ncore)+"'");
		}
	
		auto diag = get_tag_value("diagnostics");
		
		if(diag != ""){
			if(diag == "on") details.diagnostics_on = true;
			else{
				if(diag != "off"){
					alert_import("'diagnostics' has a value '"+diag+"' but must be 'off' or 'on'");
					return;
				}
			}
		}
	}
	
	model.details.stochastic = true;

	details.burnin_frac = 20;
	details.anneal_type = ANNEAL_NONE;
	details.anneal_rate = ANNEAL_RATE_DEFAULT;
	details.anneal_power = ANNEAL_POWER_DEFAULT;
	
	if(algo == DA_MCMC || algo == PAS_MCMC){
		auto burnin_str = get_tag_value("burnin-percent"); 
		
		if(burnin_str != ""){
			auto burnin = number(burnin_str);
			if(burnin == UNSET || burnin < 1 || burnin > 90){
				terminate = true;
				alert_import("For 'burnin-percent' the value '"+burnin_str+"' must be a number between 1 and 90");
				return;		
			}
			details.burnin_frac = burnin;
		}	
		else{
			details.burnin_frac = BURNIN_FRAC_DEFAULT;
		}
	}
	
	if(algo == PAS_MCMC){
		auto genper_str = get_tag_value("gen-percent"); 
		
		auto genper = PAS_GEN_UPDATE_DEFAULT;
		if(genper_str != ""){
			genper = number(genper_str);
			if(genper == UNSET || genper <= 0 || genper >= 100){
				terminate = true;
				alert_import("For 'gen-percent' the value '"+genper_str+"' must be a number between 0 and 100, exclusive");
				return;		
			}
		}	
		
		auto update = (unsigned int)(genper*details.sample/100.0);
		if(update < 10){
			update = 10;
			alert_warning("'gen-percent' set to "+tstr(update*100/details.sample)+"% because the value "+tstr(genper)+"% was too small.");
		}
		details.gen_update = update;
	}
	
	if(algo == DA_MCMC){
		auto anneal_str = get_tag_value("anneal"); 
		if(anneal_str != ""){
			details.anneal_type = AnnealType(option_error("anneal",anneal_str,{"none","scan","power-auto","log-auto","power"},{ ANNEAL_NONE, ANNEAL_SCAN, ANNEAL_POWERAUTO, ANNEAL_LOGAUTO, ANNEAL_POWER}));
			if(details.anneal_type == UNSET) return;
		
			switch(details.anneal_type){
			case ANNEAL_SCAN:
				{
					auto rate_str = get_tag_value("rate");
					if(rate_str != ""){ 
						if(!is_positive(rate_str,"rate")) return;
						details.anneal_rate = number(rate_str);
					}
				}
				break;
				
			case ANNEAL_POWER:	
				{
					auto power_str = get_tag_value("power");
					if(power_str != ""){ 
						if(!is_positive(power_str,"power")) return;
						details.anneal_power = number(power_str);
					}
				}
				break;
				
			default: break;
			}
		}
	}
}


/// Applies the post-sim command
void Input::post_sim_command()
{
	auto &details = model.details;
	
	auto start = get_tag_value("start"); 
	if(start == ""){ terminate = true; cannot_find_tag(); return;}
	if(!is_number(start,"start")){ terminate = true; return;}
	details.ppc_t_start = number(start);
	
	auto end = get_tag_value("end"); 
	if(end == ""){ terminate = true; cannot_find_tag(); return;}
	if(!is_number(end,"end")){ terminate = true; return;}
	auto end_num = number(end);
	
	if(details.ppc_t_start >= end_num){
		alert_import("The start time '"+start+"' must be before the end time '"+end+"'.");
		return;
	}
	
	if(details.ppc_t_start < details.t_start){
		alert_import("The start time '"+start+"' cannot be before the inference start time '"+tstr(details.t_start)+"'.");
		return;
	}
	
	if(details.ppc_t_start > details.t_end){
		alert_import("The start time '"+start+"' cannot be after the inference end time '"+tstr(details.t_end)+"'.");
		return;
	}
	
	if(end_num < details.t_end){
		alert_import("The end time '"+end+"' cannot be before the inference end time '"+tstr(details.t_end)+"'.");
		return;
	}
	
	details.ppc_t_end = end_num;
	details.inf_t_end = details.t_end;
	if(end_num > details.t_end) details.t_end = end_num;
	
	if(details.t_start >= details.t_end){
		terminate = true;
		alert_import("The start time '"+start+"' must be before the end time '"+end+"'");
		return;
	}

	details.number = check_pos_integer("number",PPC_NUM_DEFAULT);
	
	if(details.number%mpi.ncore != 0 && model.mode == PPC){
		alert_import("For 'number' the value '"+tstr(details.number)+"' must be a multiple of the number of cores");
	}		 
	details.num_per_core = details.number/mpi.ncore;

	details.individual_max = check_pos_integer("ind-max",INDMAX_DEFAULT);
	details.param_output_max = check_pos_integer("param-output-max",PARAM_OUTPUT_MAX_DEFAULT);
	
	details.seed = get_seed();
	
	details.ppc_resample = get_tag_value("resample");
	details.diagnostics_on = false;
}


/// Checks that a number is a positive integer
unsigned int Input::check_pos_integer(string te, unsigned int def)
{
	auto value = get_tag_value(te);
	
	if(value == "" && def != UNSET) return def;
	
	auto num = number(value);
	if(num == UNSET || num <= 0 || num != (int)num){
		alert_import("For '"+te+"' the value '"+value+"' must be a positive integer");
		num = 1;
	}
	return (unsigned int) num;
}


/// Checks that a number is positive 
double Input::check_pos_number(string te, unsigned int def)
{
	auto value = get_tag_value(te);
	
	if(value == "") return def;
	
	auto num = number(value);
	if(num == UNSET || num <= 0){
		alert_import("For '"+te+"' the value '"+value+"' must be a positive number");
	}
	return num;
}


/// Checks that a number is between zero and one
double Input::check_zero_one(string te, double def)
{
	auto value = get_tag_value(te);
	
	if(value == "") return def;
	
	auto num = number(value);
	if(num == UNSET || num <= 0 || num > 1){
		alert_import("For '"+te+"' the value '"+value+"' must be between zero and one");
	}
	return num;
}


/// Loads individual effects into the model
void Input::ind_effect_command()
{
	auto p = p_current;
	if(p == UNSET){ 
		terminate = true;
		alert_import("To add individual effects the species must be set");
		return;
	}
	
	auto name = get_tag_value("name"); if(name == ""){ cannot_find_tag(); return;}
	
	auto spl = split(name,',');
	vector <IEname> list; 
	for(auto i = 0u; i < spl.size(); i++){
		auto nam = spl[i];
		check_name_input(nam,"Individual effect name",true);
	
		IEname na; na.name = nam; na.index = UNSET;
		list.push_back(na);
	}
	
	auto A = get_tag_value("A"); 
	auto A_sparse = get_tag_value("A-sparse"); 
	auto pedigree = get_tag_value("pedigree");
	
	auto num = 0u;
	if(A != "") num++;
	if(A_sparse != "") num++;
	if(pedigree != "") num++;
	if(num > 1){
		alert_import("Cannot specify more than one of 'A', 'A-sparse' and 'pedigree'."); 
		return;
	}
	
	Amatrix A_matrix;
	A_matrix.set = false;
	
	if(num == 1){
		if(pedigree != ""){
			auto tab_ped = load_table(pedigree);
				
			if(tab_ped.ncol != 3){
				alert_import("'pedigree' file should have 3 columns"); return;
			}
				
			if(tab_ped.heading[0] != "ID" || tab_ped.heading[1] != "sire" || 
					tab_ped.heading[2] != "dam"){
				alert_import("'pedigree' file should have columns with headings 'ID', 'sire' and 'dam'"); return;
			}
			
			for(auto r = 0u; r < tab_ped.nrow; r++){
				A_matrix.ind_list.push_back(tab_ped.ele[r][0]);
			}
			
			A_matrix.hash_ind_list.create(A_matrix.ind_list);
			
			auto N = A_matrix.ind_list.size();
		
			vector < vector <double> > Ainv;    	// Creates an inverse
			Ainv.resize(N);
			for(auto r = 0u; r < N; r++) Ainv[r].resize(N,0);
			for(auto r = 0u; r < N; r++) Ainv[r][r] = 1;
			
			for(auto r = 0u; r < tab_ped.nrow; r++){
				vector <unsigned int> par;
				
				auto te = tab_ped.ele[r][1];
				if(te != "."){
					auto i = A_matrix.hash_ind_list.find(te);
					if(i == UNSET){
						alert_import("Sire '"+te+"' cannot be found in individual list");
					}
					else par.push_back(i);
				}					
				
				te = tab_ped.ele[r][2];
				if(te != "."){
					auto i = A_matrix.hash_ind_list.find(te);
					if(i == UNSET){
						alert_import("Dam '"+te+"' cannot be found in individual list");
					}
					else par.push_back(i);
				}			

				switch(par.size()){
				case 0: break;   // Both parents unknown
				case 1:          // One parent known
					{
						auto p = par[0];
						Ainv[p][p] += 1.0/3;
						Ainv[r][r] += 1.0/3;
						Ainv[p][r] -= 2.0/3;
						Ainv[r][p] -= 2.0/3;
					}
					break;
					
				case 2:          // Both parents known
					{
						auto par1 = par[0];
						auto par2 = par[1];
						
						Ainv[par1][par1] += 1.0/2;
						Ainv[par2][par2] += 1.0/2;
						Ainv[r][r] += 1.0;
						Ainv[par1][par2] += 1.0/2;
						Ainv[par2][par1] += 1.0/2;
						Ainv[par1][r] -= 1.0;
						Ainv[r][par1] -= 1.0;
						Ainv[par2][r] -= 1.0;
						Ainv[r][par2] -= 1.0;
					}
					break;
					
				default: emsg("Not right"); break;
				}				
			}
		
			auto A = invert_matrix(Ainv);
			
			tidy(A);
		
			A_matrix.value = A;
		}
		
		if(A_sparse != ""){
			A_matrix.set = true;
			
			auto ind_list = get_tag_value("ind-list"); 
			
			auto tab_ind = load_table(ind_list);
			if(tab_ind.error == true) return;
			vector <string> list;
			for(auto r = 0u; r < tab_ind.nrow; r++){ 
				list.push_back(tab_ind.ele[r][0]);
			}			
		
			A_matrix.ind_list = list;
		
			auto N = list.size();
			vector < vector <double> > val;
			val.resize(N);
			for(auto j = 0u; j < N; j++) val[j].resize(N,0);
		
			auto tab = load_table(A_sparse);
			
			if(tab.error == true) return;
		
			for(auto r = 0u; r < tab.nrow; r++){
				auto j = integer(tab.ele[r][0]);
				auto i = integer(tab.ele[r][1]);
				auto value = number(tab.ele[r][2]);
				if(j == UNSET || i == UNSET || value == UNSET){
					alert_import("Problem loading table on line "+to_string(r+1));
					return;
				}
				
				val[j][i] = value;
				val[i][j] = value;
			}
			A_matrix.value = val;
			A_matrix.hash_ind_list.create(A_matrix.ind_list);
		}
	
		if(A != ""){
			A_matrix.set = true;
			auto tab = load_table(A);
			if(tab.error == true) return;
		
			A_matrix.ind_list = tab.heading;
			if(tab.nrow != tab.ncol){
				alert_import("The file '"+tab.file+"' must contain an equal number of columns and rows."); 
				return;
			}
			
			vector < vector <double> > val;
			val.resize(tab.nrow);
			for(auto r = 0u; r < tab.nrow; r++){
				val[r].resize(tab.ncol);
				for(auto c = 0u; c < tab.ncol; c++){
					auto ele = number(tab.ele[r][c]);
					if(ele == UNSET){
						alert_import(in_file_text(tab.file)+ " the element '"+tstr(ele)+"' is not a number2 (row "+tstr(r+2)+", col "+tstr(c+1)+")");
						return;
					}
					val[r][c] = ele;
				}
			}
			A_matrix.value = val;
			A_matrix.hash_ind_list.create(A_matrix.ind_list);
		}
	}
	
	// Checks that ind effects do not already exist
	for(const auto &sp : model.species){
		for(const auto &ieg : sp.ind_eff_group){
			for(const auto &li : ieg.list){
				for(const auto &li2 : list){
					if(li.name == li2.name){
						alert_import("Individual effect '"+li2.name+"' already exists in the model");
						return;
					}
				}
			}
		}
	}
	
	IEgroup ieg; 
	ieg.list = list; 
	ieg.A_matrix = A_matrix;
	ieg.line_num = line_num;
	ieg.ppc_resample = false;
	model.species[p].ind_eff_group.push_back(ieg);
}


/// Loads fixed effects into the model
void Input::fixed_effect_command()
{
	auto p = p_current;
	if(p == UNSET){
		terminate = true;
		alert_import("To add individual effects the species must be set"); 
		return;
	}
	
	auto name = get_tag_value("name"); if(name == ""){ cannot_find_tag(); return;}

	check_name_input(name,"Individual fixed effect name",true);

	auto X = get_tag_value("X"); if(X == ""){ cannot_find_tag(); return;}
	
	Xvector X_vector;
	
	auto tab = load_table(X);
	if(tab.error == true) return;
	
	vector <string> col_name; 
	col_name.push_back("ID");
	col_name.push_back("value"); 
	
	auto subtab = get_subtable(tab,col_name); if(subtab.error == true) return;
	
	for(auto r = 0u; r < subtab.nrow; r++){
		auto id = subtab.ele[r][0];
		X_vector.ind_list.push_back(id);
	
		auto val = number(tab.ele[r][1]);
		if(val == UNSET){
			alert_import(in_file_text(tab.file)+" the element '"+tstr(val)+"' is not a number3 (row "+tstr(r+2)+")");
			return;
		}
		X_vector.value.push_back(val);
	}
	
	X_vector.hash_ind_list.create(X_vector.ind_list);

	auto fix_name = fe_char+"^"+name;
	
	auto th = 0u; while(th < model.param.size() && model.param[th].name != fix_name) th++;
	if(th == model.param.size()){
		alert_import("Parameter '"+fix_name+"' is not specified by 'param'");
		return;
	}
	model.param[th].used = true;
							
	FixedEffect fe;
	fe.name = name;    
	fe.th = th;	
	fe.X_vector = X_vector;
	fe.line_num = line_num;       
	model.species[p].fix_effect.push_back(fe);
}

	
/// Ignores map command 
void Input::map_command()
{
	auto file = get_tag_value("file"); 
}


/// Loads inference states into the model (for PPC)
void Input::inf_state_command()
{
	auto file = get_tag_value("file");
	auto chain = get_tag_value("chain"); 
	
	auto i = 0u; while(i < files.size() && files[i].name != file) i++;
	if(i == files.size()){
		alert_import("Could not find '"+file+"'");
		return;
	}

	const auto &flines = files[i].lines;

	auto li = 0u;
	
	// Reads in individual key
	vector <string> ind_key;
	{
		string warn = "Problem loading state file";
		
		while(li < flines.size() && trim(flines[li]) != "{") li++;
		if(li == flines.size()) alert_import(warn+"1");
		while(li < flines.size() && !begin_str(flines[li],"timepoint")) li++;
		if(li == flines.size()) alert_import(warn+"2");
		li++;
		while(li < flines.size()){
			auto line = trim(flines[li]);
			if(line == "}") break;
			if(line != ""){
				auto spl = split(flines[li],':');
				if(spl.size() != 2) alert_import(warn+"3");
				auto num = number(spl[0]);
				auto name = spl[1];
				if(num >= ind_key.size()) ind_key.resize(num+1);
				ind_key[num] = name;
			}
			li++;
		}
		
		while(li < flines.size() && !begin_str(flines[li],"<<")) li++;
		if(li == flines.size()) alert_import(warn+"4");
	}
	
	vector <string> lines; 
	
	while(li < flines.size()){	
		auto lin = trim(flines[li]);
		if(lin.length() > 2 && lin.substr(0,2) == "<<"){
			if(lines.size() > 0) read_state_sample(lines,ind_key);
			lines.clear();
		}
		lines.push_back(lin);
		li++;
	}

	if(lines.size() > 0) read_state_sample(lines,ind_key);
}

/// Warning command
void Input::warning_command()
{
	auto text = get_tag_value("text"); 
}


/// Ignores sample command 
void Input::dummy_file_command()
{
	auto file = get_tag_value("file"); 
	auto chain = get_tag_value("chain"); 
}

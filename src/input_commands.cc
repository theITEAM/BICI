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
	
	auto name = get_tag_value("name"); 
	if(name == ""){ 
		switch(cname){
		case INIT_POP: name = "Initial population"; break;
		case ADD_IND: name = "Added individuals"; break;
		case REMOVE_IND: name = "Removed individuals"; break;
		case ADD_POP: name = "Added populations"; break;
		case REMOVE_POP: name = "Removed populations"; break;
		default: cannot_find_tag(); return;
		}
	}
	
	DataSource ds;
	ds.name = name;
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
					auto pri = convert_text_to_prior(prior,line_num,"Total population",false);
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
				alert_import("Value '"+spl[0]+"' is not a compartment");
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
			
			auto trans = get_tag_value("trans");
			ds.filter_trans_str = trans;
		
			auto cl_sel = get_cl_from_trans(trans,p);
		
			if(cl_sel == UNSET){
				alert_import("For 'trans' the value '"+trans+"' is not recognised");
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
				{
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
				}
				break;
			}
		}
		break;
		
	case POP_TRANS_DATA:
		{
			ds.filter_str = get_tag_value("filter");
			
			auto trans = get_tag_value("trans");
			ds.filter_trans_str = trans;
			
			auto cl_sel = get_cl_from_trans(trans,p);
		
			if(cl_sel == UNSET){
				alert_import("For 'trans' the value '"+trans+"' is not recognised");
				return;
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
		
	default: alert_emsg_input("Should not be default2"); return;
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
			alert_import("'cols' does not have the correct number of entries (expected a table with headings in the order '"+stringify(col_name)+"')"); 
			return;
		}
		
		for(auto i = 0u; i < spl.size(); i++){
			for(auto j = 0u; j < col_name.size(); j++){
				if(spl[i] == col_name[j] && i != j){
					alert_import("'cols' does not have the correct order (expected a table with headings in the order '"+stringify(col_name)+"')"); 
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
bool Input::species_command(unsigned int loop)
{
	auto name = get_tag_value("name");
	if(name == ""){ cannot_find_tag(); return false;}
	
	auto type = toLower(get_tag_value("type"));
	if(type == ""){ cannot_find_tag(); return false;}
	
	auto sp_type = SpeciesType(option_error("type",type,{"population","individual"},{ POPULATION, INDIVIDUAL}));
	if(sp_type == UNSET) return false;

	auto trans_tree = false;
	if(sp_type == INDIVIDUAL){
		auto trans_tree_str = toLower(get_tag_value("trans-tree"));
		if(trans_tree_str != ""){
			if(trans_tree_str == "on") trans_tree = true;
			else{
				if(trans_tree_str != "off"){ 
					alert_import("'trans-tree' must be either 'off' or 'on'");
					return false;
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
			return false;
		}
	}
	cl_current = UNSET;
	
	return true;
}



/// Adds a classification to the model
bool Input::classification_command(unsigned int loop)
{
	auto p = p_current;
	if(p == UNSET){ 
		alert_import("A species must be defined before a classification can be added"); 
		return false;
	}
		
	auto &sp = model.species[p];
	
	auto name = get_tag_value("name");
	if(name == ""){ cannot_find_tag(true); return false;}
	
	auto clone = get_tag_value("clone");
	if(clone != ""){ // Clones a classification from another species
		if(loop == 2){
			auto p2 = find_p(clone);
			if(p2 == UNSET){
				alert_import("For 'clone' cannot find the species '"+clone+"'");
				return false;
			}
			
			auto cl2 = find_cl(p2,name);
			
			if(cl2 == UNSET){
				alert_import("In clone species '"+clone+"' cannot find the classification '"+name+"'");
				return false;
			} 
				
			auto ncl = sp.ncla;
			clone_class(p,p2,cl2);
			if(ncl != sp.ncla) cl_current = sp.ncla-1;		
		}
		else{
			cl_current = find_cl(p,name);
			if(cl_current == UNSET){
				alert_import("Problem finding classification '"+name+"'"); 
				return false;
			}
		}
	}
	else{
		auto index = get_tag_value("index");
		if(index == ""){ cannot_find_tag(true); return false;}
		if(index.length() != 1){ 
			alert_import("The index '"+index+"' must be just a single character"); 
			return false;
		}
		
		if(index == "t"){ 
			alert_import("The index 't' cannot be used because it is reserved for time variation"); 
			return false;
		}
		
		if(index == "z"){ 
			alert_import("The index 'z' cannot be used because it is reserved for covariance matrices"); 
			return false;
		}
		
		if(find_in(alphabet,index) == UNSET){ 
			alert_import("Index '"+index+"' must be from the lowercase alphabet"); 
			return false;
		}
		
		auto coord = toLower(get_tag_value("coord"));
		
		auto coord_type = CARTESIAN;
		
		if(coord != ""){
			coord_type = Coord(option_error("coord",coord,{"cartesian", "latlng"},{ CARTESIAN, LATLNG}));
			if(coord_type == UNSET) return false;
		}
		
		if(coord == "latlng"){
			auto def_map = get_tag_value("map");
			if(def_map != ""){
				if(def_map != "load"){ 
					alert_import("'map' must be set to 'load'");
					return false;
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
				return false;
			}
		}
	}
	
	return true;
}


/// Loads up proposal information
void Input::proposal_info_command()
{
	auto file = get_tag_value("file"); if(file == ""){ cannot_find_tag(); return;}
	
	auto i = 0u; while(i < files.size() && files[i].name != file) i++;
	if(i == files.size()){
		alert_import("Could not find the file '"+file+"'");
		return;
	}

	auto ch = get_chain();

	// For ext in -core mode only need to load states for that core
	if(model.mode == EXT && mpi.core != ch/model.details.num_per_core) return;
	
	auto err = model.load_prop_info(ch,files[i].lines);
	
	if(err != "") alert_import("Proposal information could not be loaded: "+err);
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

	Param par(model.constant); 
	par.variety = CONST_PARAM;
	par.time_dep = true;
	par.spline_out.on = false;
	par.spline_info.on = true;
	par.spline_outside = false;
	par.spline_outside_par = false;
	par.used = true;
	par.param_mult = UNSET;
	par.factor = true;
	par.line_num = line_num;
	par.cat_factor = false;
	par.cat_factor_weight_on = false;
	par.default_prior_ref = UNSET;
	
	auto knot_times_str = get_tag_value("knot-times"); 
	if(knot_times_str == ""){ cannot_find_tag(); return;}
	
	auto type = LINEAR_SPL;
		
	auto sptype = trim(toLower(get_tag_value("spline-type")));
	if(sptype != ""){
		type = SplineType(option_error("spline-type",sptype,{"linear","square","cubic +ve","cubic"},{ LINEAR_SPL, SQUARE_SPL, CUBICPOS_SPL, CUBIC_SPL}));
	}
		
	vector <string> knot_times, knot_times_out;
	set_spline(type,knot_times_str,"",knot_times,knot_times_out,false,par);
	
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
		
	auto mult = get_dependency(par.dep,pp,knot_times,knot_times_out); if(mult == UNSET) return; 
	
	par.element_ref.resize(mult);
	for(auto &er : par.element_ref) er.index = UNSET;
	
	par.N = mult;
	par.trace_output = false;
	
	auto file = get_tag_value("constant"); if(file == ""){ cannot_find_tag(); return;}
	load_param_value(pp,file,par,"In 'file'",VALUE_LOAD);

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
			alert_import("A species must be specified before a 'classification' can be set"); 
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
		vector <Tag> tags; 
		for(const auto &ta : cline_store.tags){
			if(ta.name != "file"){
				tags.push_back(ta);
			}
		}
		
		for(auto i = 0u; i < N; i++){
			auto te = row[i];
			te = trim(te);
			if(te != ""){
				Tag tag; tag.name = trim(head[i]); tag.value = te; tag.done = 0;
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
			alert_import("The value '"+tag.value+"' for tag '"+tag.name+"' is not used (line "+to_string(r+2)+" in file).");
			return;
		}			
	}
}


/// Adds a compartment to the model
bool Input::compartment_command()
{
	return compartment_command2(cline_store.tags);
}


/// Adds all compartments to the model
bool Input::compartment_all_command()
{
	auto file = get_tag_value("file"); if(file == ""){ cannot_find_tag(true); return false;}

	auto tags_list = get_tags_list(file);
	
	for(auto i = 0u; i < tags_list.size(); i++){
		error_line = in_file_text(file)+" on line "+tstr(i+2)+": ";
		
		auto &tags = tags_list[i];
	
		if(!compartment_command2(tags)) return false;
		
		check_tags_used(i,tags);
	}
	
	for(auto &ta : cline_store.tags) ta.done = 1;
	
	error_line = "";
	
	return true;
}


/// Adds a compartment to the model
bool Input::compartment_command2(vector <Tag> &tags)
{
	if(check_claa_error() == true) return false;

	const auto &claa = get_claa();
	auto p = p_current;
	auto cl = cl_current;
	
	auto name = get_tag_val("name",tags); if(name == ""){ cannot_find_tag(true); return false;}

	// Checks that a compartment name is valid
	check_name_input(name,"Compartment name");

	auto color = get_tag_val("color",tags); 
	if(color != ""){
		if(!is_Color(color)){
			alert_import("'"+color+"' is not a valid color.",true); 
			return false;
		}
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
				if(x_str == ""){ cannot_find_tag(true); return false;}
				if(!is_number(x_str,"x",true)) return false;
				x = number(x_str); 
			
				y_str = get_tag_val("y",tags);
				if(y_str == ""){ cannot_find_tag(true); return false;}
				if(!is_number(y_str,"y",true)) return false;
				y = number(y_str);
				
				if(x != UNSET && y != UNSET){
					for(auto i = 0u; i < claa.ncomp; i++){
						if(claa.comp[i].x == x && claa.comp[i].y == y){		
							alert_import("Compartments '"+name+"' and '"+claa.comp[i].name+"' cannot have the same position");	
							return false;
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
				auto i = import_geojson(bound_file); if(i == UNSET) return false;
				
				auto po = boundary_mean_latlng(i,name);
			
				lat = po.lat; lng = po.lng;
			}
			else{
				auto lat_str = get_tag_val("lat",tags);
				if(!is_number(lat_str,"lat",true)) return false;
				lat = number(lat_str);
				
				auto lng_str = get_tag_val("lng",tags);
				if(!is_number(lng_str,"lng",true)) return false;
				lng = number(lng_str);
			
				if(lat != UNSET && lng != UNSET){
					if(check_latlng_error(lat,lng,true) == true) return false;
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
				return false;
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
				return false;
			}
		}
	}
	
	if(!add_compartment(name,p,cl,x,y,lat,lng,markov_branch,infected,"")) return false;

	return true;
}


/// Adds a transition / source / sink to the model
bool Input::transition_command()
{
	return transition_command2(cline_store.tags);
}


/// Adds a transition / source / sink to the model
bool Input::transition_all_command()
{
	auto file = get_tag_value("file"); if(file == ""){ cannot_find_tag(); return false;}

	auto tags_list = get_tags_list(file);

	for(auto i = 0u; i < tags_list.size(); i++){
		error_line = in_file_text(file)+" on line "+tstr(i+2)+": ";
		
		auto &tags = tags_list[i];

		if(!transition_command2(tags)) return false;
		check_tags_used(i,tags);
	}
	
	for(auto &ta : cline_store.tags) ta.done = 1;
	
	error_line = "";
	
	return true;
}


/// Adds a transition / source / sink to the model
bool Input::transition_command2(vector <Tag> &tags)
{
	if(check_claa_error() == true) return false;
	auto &claa = get_claa();

	auto p = p_current;
	auto cl = cl_current;
	
	auto te = get_tag_val("name",tags); if(te == ""){ cannot_find_tag(); return false;}
	
	if(te.length() < 2){ 
		alert_import("Transition name '"+te+"' is not understood. Expected format: 'Initial compartment'->'final compartment'."); 
		return false;
	}
	
	auto i = 0u; while(i < te.length()-1 && te.substr(i,2) != "->") i++;
	
	if(i == te.length()-1){ 
		alert_import("Transition name '"+te+"' is not understood. Expected format: 'Initial compartment'->'final compartment'."); 
		return false;
	}
	
	string fr = trim(te.substr(0,i));
	if(fr == "+") fr = "Source";
	
	string to = trim(te.substr(i+2));
	if(to == "-") to = "Sink";
	
	auto value_str = get_tag_val("value",tags); 
	if(value_str == ""){ cannot_find_tag(); return false;}

	auto trans_def = extract_trans_def(value_str);

	if(trans_def.set == false){ 
		alert_import("In the expression '"+value_str+"' the transition distribution is incorrectly specified. Expected format: 'exp(rate:...)', 'gamma(mean:...,cv:...)', 'erlang(mean:...,shape:)', 'log-normal(mean:...,cv:...)', 'weibull(scale:...,shape:...)' or 'period(time:...)'.");
		return false;
	}
	
	unsigned int ci, cf;
	
	if(fr == "Source") ci = UNSET;
	else{
		ci = find_c(p,cl,fr);
		if(ci == UNSET){ cannot_find_trans_comp(te,p,cl,fr); return false;}
	}
	
	if(to == "Sink") cf = UNSET;
	else{
		cf = find_c(p,cl,to);
		if(cf == UNSET){
			cannot_find_trans_comp(te,p,cl,to);
			return false;
		}
	}

	for(auto tr = 0u; tr < claa.ntra; tr++){
		if(claa.tra[tr].i == ci && claa.tra[tr].f == cf){
			if(ci == UNSET){ 
				alert_import("Source to '"+to+"' already exists"); 
				return false;
			}
			else{
				if(cf == UNSET){
					alert_import("Sink from '"+fr+"' already exists"); 
					return false;
				}
				else{
					alert_import("Transition from '"+fr+"' to '"+to+"' already exists");
					return false;
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
			
	if(!add_transition(p,cl,ci,cf,trans_def.type)) return false;

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
		
			stringstream ss;
			ss << fixed << std::setprecision(16);
			ss << period;
			 
			tra.dist_param.push_back(he(add_equation_info(ss.str(),TRANS_MEAN,p,cl)));
		}
		break;
	}
	
	return true;
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


/// Sets any defined equations
void Input::define_command()
{
	auto full_name = get_tag_value("name"); if(full_name == ""){ cannot_find_tag(); return;}
	auto value = get_tag_value("value"); if(value == ""){ cannot_find_tag(); return;}
	
	auto pp = get_param_prop(full_name);
	
	Define def;
	def.name = pp.name;

	def.line_num = line_num;
	def.full_name = full_name;
	def.time_dep = pp.time_dep;
	
	auto mult = get_dependency(def.dep,pp,vector <string> (),vector <string> ()); 
	if(mult == UNSET) return; 
	
	auto def_eqn_raw = he(add_equation_info(value,DEFINE_EQN));
	def.value = def_eqn_raw;
	
	string warn = "";
	auto eq_dep = model.equation_dep(def_eqn_raw.te,warn);
	if(warn != ""){ alert_import(warn,true); return;}
	
	auto res = dep_agree(full_name,def.dep,eq_dep);
	if(res != ""){ alert_import(res,true); return;}
		
	const auto &depend = def.dep;
		
	vector <DepConv> dep_conv;
	for(auto d = 0u; d < depend.size(); d++){
		const auto &dep = depend[d];
		DepConv dc; 
		dc.before = dep.index_with_prime;
		dep_conv.push_back(dc);
	}
	
	auto swap_temp = swap_template(def_eqn_raw.te,dep_conv);
	if(swap_temp.warn != ""){ 
		alert_import(swap_temp.warn); 
		return;
	}
	
	def.swap_temp = swap_temp;
	
	model.define.push_back(def);
}

	
/// Sets the value for a parameter in the model
void Input::param_command()
{
	auto full_name = get_tag_value("name"); if(full_name == ""){ cannot_find_tag(); return;}
	
	auto pp = get_param_prop(full_name);

	Param par(model.constant); 
	par.variety = UNSET_PARAM;
	par.name = pp.name;
	par.full_name = full_name;
	par.time_dep = pp.time_dep;
	par.spline_out.on = false;
	par.spline_info.on = false;
	par.spline_outside = false;
	par.spline_outside_par = false;
	par.used = false;
	par.param_mult = UNSET;
	par.factor = false;
	par.line_num = line_num;
	par.default_prior_ref = UNSET;
		
	vector <string> knot_times, knot_times_out;
	if(pp.time_dep == true){
		auto knot_times_str = get_tag_value("knot-times"); if(knot_times_str == ""){ cannot_find_tag(); return;}
		
		auto smooth = trim(toLower(get_tag_value("smooth")));
		
		auto type = LINEAR_SPL;
		
		auto sptype = trim(toLower(get_tag_value("spline-type")));
		if(sptype != ""){
			type = SplineType(option_error("spline-type",sptype,{"linear","square","cubic +ve","cubic"},{ LINEAR_SPL, SQUARE_SPL, CUBICPOS_SPL, CUBIC_SPL}));
		}
		
		set_spline(type,knot_times_str,smooth,knot_times,knot_times_out,true,par);
	}
	
	auto cons = get_tag_value("constant"); 
	auto dist = get_tag_value("dist");
	auto dist_split = get_tag_value("dist-split");
	auto value = get_tag_value("value");
	auto reparam = get_tag_value("reparam"); 
	auto prior = get_tag_value("prior"); 
	auto prior_split = get_tag_value("prior-split"); 

	auto mode = model.mode;

	auto sim_sample = toLower(get_tag_value("sim-sample")); 
	if(sim_sample == "") sim_sample = "true";
	else{
		if(dist == "" && dist_split == ""){
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
	
	switch(mode){
	case INF: case PPC: case EXT:
		value = "";
		break;
		
	case SIM: 
		prior = ""; prior_split = "";
		if(par.sim_sample == false){ dist = ""; dist_split = "";}
		break;
	default: break;
	}
	
	vector <ParamTag> param_tag;
	ParamTag pt; 
	pt.val = cons; pt.tag = "constant"; param_tag.push_back(pt);
	pt.val = value; pt.tag = "value"; param_tag.push_back(pt);
	pt.val = dist; pt.tag = "dist"; param_tag.push_back(pt);
	pt.val = dist_split; pt.tag = "dist-split"; param_tag.push_back(pt);
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
	
	auto mult = get_dependency(par.dep,pp,knot_times,knot_times_out); if(mult == UNSET) return; 
	
	par.element_ref.resize(mult);
	for(auto &er : par.element_ref) er.index = UNSET; 
	
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
		par.weight.resize(par.N,1);
		
		auto weight = get_tag_value("factor-weight"); 
		if(weight != ""){
			par.cat_factor_weight_on = true;
			if(is_file(weight) == false){
				alert_import("'factor-weight' must be defined through a data table");
			}		
			else{
				string desc = "For 'factor-weight'";
				load_param_value(pp,weight,par,desc,FACW_LOAD);
			}
		}
	}
	
	par.trace_output = true;
	if(cons != "" || mult > model.details.param_output_max) par.trace_output = false;

	auto warn = check_reserved_name(par.name);	
	if(warn != "") alert_import(warn);
	
	auto pre = "Parameter '"+par.full_name+"': ";
	
	if(value == "auto") alert_emsg_input("auto no longer supported");
		
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
			auto eqn = he(add_equation_info(reparam,REPARAM_EQN));
			par.reparam_eqn = eqn;
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
						alert_emsg_input("Problem importing"); 
						return;
					}
				}
			}
			
			if(reparam != "" && is_file(valu) == false){
				auto eqn = he(add_equation_info(reparam,REPARAM_EQN));

				par.reparam_eqn = eqn;
			}
			else{
				if(is_file(valu) == false){
					double val = number(valu);
					
					if(par.variety != CONST_PARAM) alert_emsg_input("Should be const");
				
					if(val == UNSET){
						alert_import(desc+" '"+valu+"' is not a number");
						return;
					}
						
					ElementRef er; er.index = par.add_cons(val); er.cons = true;
					for(auto k = 0u; k < mult; k++){
						par.element_ref[k] = er;
					}
				}
				else{
					load_param_value(pp,valu,par,desc,VALUE_LOAD);
				}
			}
		}
	}

	if(par.reparam_eqn.te != ""){
		string warn = "";
		auto eq_dep = model.equation_dep(par.reparam_eqn.te,warn);
		if(warn != ""){ alert_import(warn,true); return;}
			
		auto res = dep_agree(full_name,par.dep,eq_dep);
		if(res != ""){ alert_import(res,true); return;}
	}

	if(prior != ""){
		par.variety = PRIOR_PARAM;
		
		auto pri = convert_text_to_prior(prior,line_num,par.full_name,false);
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
		
		if(model.is_symmetric(par)){                             // Deals with covariance matrix 
			model.prior.push_back(pri);
			auto L = par.dep[0].list.size();
			for(auto i = 0u; i < L; i++) par.set_prior(i*L+i,ref); // Places prior along diagonal
			
			if(L > 1){                                             // Places correlations off diagonal
				Prior pri_off;
				pri_off.name = "Correlation";
				pri_off.type = MVN_COR_PR;
				auto ref_cor = model.prior.size();
				for(auto j = 0u; j < L; j++){
					for(auto i = 0u; i < L; i++){
						if(i != j) par.set_prior(j*L+i,ref_cor); 
					}
				}
				model.prior.push_back(pri_off);
			}
			
			if(false){
				for(auto j = 0u; j < L; j++){
					for(auto i = 0u; i < L; i++){
						auto pri = model.prior[par.element[par.element_ref[j*L+i].index].prior_ref];
						cout << pri.type << ",";
					}
					cout << " pr" << endl;
				}
			}				
		}
		else{
			model.prior.push_back(pri);
			for(auto i = 0u; i < mult; i++) par.set_prior(i,ref);
		}
	}
	
	if(prior_split != ""){
		par.variety = PRIOR_PARAM;
		
		if(par.dep.size() == 0){
			alert_import("'prior-split' can only be used if the parameter has a dependency."); 
			return;
		}
	
		if(par.cat_factor){
			alert_import("'prior-split' cannot be used for a factor."); 
			return;
		}
		
		if(model.is_symmetric(par)){
			alert_import("'prior-split' cannot be used for a covariance matrix."); 
			return;
		}
		
		string desc = "For 'prior-split'";
		load_param_value(pp,prior_split,par,desc,PRIOR_SPLIT_LOAD);
	}
	
	if(dist != ""){
		par.variety = DIST_PARAM;
		
		auto pri = convert_text_to_prior(dist,line_num,par.full_name,true);
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
		
		string desc = "For 'dist-split'";
		load_param_value(pp,dist_split,par,desc,DIST_SPLIT_LOAD);
	}
		
	for(const auto &par2 : model.param){
		if(par2.name == par.name){
			alert_import("Parameter '"+par.name+"' is already defined"); 
			return;
		}
	}
	
	if(par.variety == UNSET_PARAM){
		if(par.not_set){
			switch(model.mode){	
			case INF: case PPC: case EXT:
				alert_import("A prior should be set for parameter '"+par.full_name+"'");
				break;
	
			case SIM: 
				alert_warning("A value has not been set for parameter '"+par.full_name+"'. Ignore this warning if this parameter is used in the observation process.");
				break;
				
			default: emsg("option he"); break;
			}
		}
		else{
			alert_import("Parameter variety for '"+par.name+"' is unset");
		}
		return;
	}
	else{	
		// Adds an extra parameter to model categorical factor
		if(par.cat_factor && (model.mode == INF || model.mode == EXT)){ 
			add_param_cat_factor(par);
		}
		
		model.param.push_back(par);
	}

	if(par.variety == REPARAM_PARAM){	
		for(const auto &ele : par.element){		
			he(ele.value);
		}
	}
}


/// Definition for derived quantities
void Input::derived_command()
{
	auto full_name = get_tag_value("name"); if(full_name == ""){ cannot_find_tag(); return;}
	auto eqn_name = get_tag_value("eqn"); if(eqn_name == ""){ cannot_find_tag(); return;}
	
	auto pp = get_param_prop(full_name);
	auto dep_eqn = pp.dep_with_prime;
	if(pp.time_dep == true) pp.dep.pop_back();
	
	auto warn = check_reserved_name(pp.name);
	if(warn != "") alert_import(warn);
	
	Derive der;
	der.name = pp.name;
	der.line_num = line_num;
	der.full_name = full_name;
	der.time_dep = pp.time_dep;
	der.func.on = false;
	
	auto te = trim(eqn_name);

	auto df_type = DF_UNSET;
	if(begin_str(te,RN_name+"(")) df_type = RN;
	if(begin_str(te,RNE_name+"(")) df_type = RNE;
	if(begin_str(te,RNC_name+"(")) df_type = RNC;
	if(begin_str(te,GT_name+"(")) df_type = GT;
	if(begin_str(te,GTE_name+"(")) df_type = GTE;
	if(begin_str(te,GTC_name+"(")) df_type = GTC;
	
	if(df_type != DF_UNSET){
		setup_der_func(df_type,eqn_name,der.func);
	}
	else{			
		auto mult = get_dependency(der.dep,pp,vector <string> (),vector <string> ()); 
		if(mult == UNSET) return; 
		
		const auto &depend = der.dep;
		
		auto der_eqn_raw = he(add_equation_info(eqn_name,DERIVE_EQN));
		
		der.eq_raw = der_eqn_raw;
		
		string warn = "";
		auto eq_dep = model.equation_dep(der_eqn_raw.te,warn);
		if(warn != ""){ alert_import(warn,true); return;}
			
		auto res = dep_agree(full_name,dep_eqn,eq_dep);
		if(res != ""){ alert_import(res,true); return;}
						
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
					alert_emsg_input("Swap index dif res");
				}
			}
			
			der.eq.push_back(der_eqn);
		}	
	}
									
	model.derive.push_back(der);
}


/// Applies the simulation command
bool Input::simulation_command()
{
	auto &details = model.details;
	
	auto start = get_tag_value("start"); 
	if(start == ""){ cannot_find_tag(); return false;}
	if(!is_number(start,"start")) return false;
	details.t_start = number(start);	
	
	auto end = get_tag_value("end");
	if(end == ""){ cannot_find_tag(); return false;}

	if(!is_number(end,"end")) return false;
	auto end_num = number(end);
	
	details.t_end = end_num;
	
	if(details.t_start >= details.t_end){
		alert_import("The start time '"+start+"' must be before the end time '"+end+"'");
		return false;
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
	if(details.algorithm == UNSET) return false;
	
	auto dt_str = get_tag_value("timestep"); 
	if(dt_str == ""){ cannot_find_tag(); return false;}
	
	if(!is_positive(dt_str,"timestep")) return false;
	
	details.dt = number(dt_str);
	details.stochastic = true;
	details.individual_max = check_pos_integer("ind-max",INDMAX_DEFAULT);
	details.param_output_max = check_pos_integer("param-output-max",PARAM_OUTPUT_MAX_DEFAULT);
	details.anneal_type = ANNEAL_NONE;
	details.anneal_rate = UNSET;
	details.anneal_power = UNSET;
	details.diagnostics_on = false;
	
	if(check_dt(details) == false) return false;
	
	return true;
}


/// Applies the inference command
bool Input::inference_command()
{
	auto &details = model.details;

	auto start = get_tag_value("start"); 
	if(start == ""){ cannot_find_tag(); return false;}
	if(!is_number(start,"start")) return false;
	details.t_start = number(start);
	
	auto end = get_tag_value("end"); 
	if(end == ""){ cannot_find_tag(); return false;}
	if(!is_number(end,"end")) return false;
	details.t_end = number(end);
	
	if(details.t_start >= details.t_end){
		alert_import("The start time '"+start+"' must be before the end time '"+end+"'");
		return false;
	}
	
	details.seed = get_seed();
	
	details.algorithm = DA_MCMC;
	
	auto alg = toUpper(get_tag_value("algorithm"));
	
	if(alg != ""){
		details.algorithm = Algorithm(option_error("algorithm",alg,{"DA-MCMC","PAS-MCMC","MFA","ABC","ABC-SMC","ABC-MBP","PMCMC","HMC"},{ DA_MCMC, PAS_MCMC, MFA_ALG, ABC_ALG, ABC_SMC_ALG, ABC_MBP, PMCMC, HMC }));
		if(details.algorithm == UNSET) return false;
	}
	
	auto dt_str = get_tag_value("timestep"); 
	if(dt_str == ""){ cannot_find_tag(); return false;}

	if(!is_positive(dt_str,"timestep")) return false;
	
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
		if(details.nchain%mpi.ncore != 0 && (model.mode == INF || model.mode == EXT)){
			alert_import("For 'nchain' the value '"+tstr(details.nchain)+"' must be a multiple of the number of cores");
		}
		details.num_per_core = check_pos_integer("chain-per-core",MCMC_CHAIN_PER_CORE_DEFAULT);
		break;
		
	case ABC_ALG:
		{
			details.sample = check_pos_integer("sample",ABC_SAMPLE_DEFAULT);
			details.accfrac = check_zero_one("acc-frac",ABC_ACFRAC_DEFAULT);
			auto num = (unsigned int)(details.sample/details.accfrac);
			
			if(num%mpi.ncore != 0 && (model.mode == INF || model.mode == EXT)){
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
		
		if(details.sample%mpi.ncore != 0 && (model.mode == INF || model.mode == EXT)){
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
		
		if(details.nchain%mpi.ncore != 0 && (model.mode == INF || model.mode == EXT)){
			alert_import("For 'npart' the value '"+tstr(details.nchain)+"' must be a multiple of the number of cores");
		}
		details.num_per_core = check_pos_integer("part-per-core",PAS_PART_PER_CORE_DEFAULT);
		break;
		
	default: break;
	}
	
	model.sync_on = true;
	auto sync = get_tag_value("sync");
	auto sync_lo = toLower(sync);
	
	if(sync_lo != ""){
		if(sync_lo == "on") model.sync_on = true;
		else{
			if(sync_lo == "off") model.sync_on = false;
			else{
				alert_import("For 'sync' the value '"+sync+"' must be a multiple of the number of cores");
			}
		}
	}
		
	if(com_op){
		if(details.num_per_core != details.nchain){
			model.sync_on = false;
		}
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
		if(ncore != mpi.ncore && (model.mode == INF || model.mode == EXT) && mpi.core_spec_on == false){
			alert_import("The number of cores is '"+tstr(mpi.ncore)+"' and should be '"+tstr(ncore)+"'");
		}
	
		auto diag = get_tag_value("diagnostics");
		
		if(diag != ""){
			if(diag == "on") details.diagnostics_on = true;
			else{
				if(diag != "off"){
					alert_import("'diagnostics' has a value '"+diag+"' but must be 'off' or 'on'");
					return false;
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
				alert_import("For 'burnin-percent' the value '"+burnin_str+"' must be a number between 1 and 90");
				return false;		
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
				alert_import("For 'gen-percent' the value '"+genper_str+"' must be a number between 0 and 100, exclusive");
				return false;		
			}
		}	
		
		auto update = (unsigned int)(genper*details.sample/100.0);
		if(update < 5){
			update = 5;
			alert_warning("'gen-percent' set to "+tstr(update*100/details.sample)+"% because the value "+tstr(genper)+"% was too small.");
		}
		details.gen_update = update;
	}
	
	if(algo == DA_MCMC){
		auto anneal_str = get_tag_value("anneal"); 
		if(anneal_str != ""){
			details.anneal_type = AnnealType(option_error("anneal",anneal_str,{"none","scan","power-auto","log-auto","power"},{ ANNEAL_NONE, ANNEAL_SCAN, ANNEAL_POWERAUTO, ANNEAL_LOGAUTO, ANNEAL_POWER}));
			if(details.anneal_type == UNSET) return false;
		
			switch(details.anneal_type){
			case ANNEAL_SCAN:
				{
					auto rate_str = get_tag_value("rate");
					if(rate_str != ""){ 
						if(!is_positive(rate_str,"rate")) return false;
						details.anneal_rate = number(rate_str);
					}
				}
				break;
				
			case ANNEAL_POWER:	
				{
					auto power_str = get_tag_value("power");
					if(power_str != ""){ 
						if(!is_positive(power_str,"power")) return false;
						details.anneal_power = number(power_str);
					}
				}
				break;
				
			default: break;
			}
		}
	}
	
	if(check_dt(details) == false) return false;
	
	switch(algo){
	case DA_MCMC: case PAS_MCMC:
		if(model.mode == EXT && !model.ext_factor.percent){
			if(model.ext_factor.value <= details.sample){
				alert_import("The extended number of updates "+tstr(model.ext_factor.value)+" must be larger than the current number "+tstr(details.sample)+".");
			}
		}
		break;
		
	default: break;
	}
	
	return true;
}


/// Applies the post-sim command
bool Input::post_sim_command()
{
	auto &details = model.details;
	
	auto start = get_tag_value("start"); 
	if(start == ""){ cannot_find_tag(); return false;}
	if(!is_number(start,"start")) return false;
	details.ppc_t_start = number(start);
	
	auto end = get_tag_value("end"); 
	if(end == ""){ cannot_find_tag(); return false;}
	if(!is_number(end,"end")) return false;
	auto end_num = number(end);
	
	auto dt = get_tag_value("timestep"); 
	if(dt != ""){
		alert_import("'time-step' should not be set for 'post-sim'");
	}
	
	if(details.ppc_t_start >= end_num){
		alert_import("The start time '"+start+"' must be before the end time '"+end+"'.");
		return false;
	}
	
	if(details.ppc_t_start < details.t_start){
		alert_import("The start time '"+start+"' cannot be before the inference start time '"+tstr(details.t_start)+"'.");
		return false;
	}
	
	if(details.ppc_t_start > details.t_end){
		alert_import("The start time '"+start+"' cannot be after the inference end time '"+tstr(details.t_end)+"'.");
		return false;
	}
	
	if(end_num < details.t_end){
		alert_import("The end time '"+end+"' cannot be before the inference end time '"+tstr(details.t_end)+"'.");
		return false;
	}
	
	details.ppc_t_end = end_num;
	details.inf_t_end = details.t_end;
	if(end_num > details.t_end) details.t_end = end_num;
	
	if(details.t_start >= details.t_end){
		alert_import("The start time '"+start+"' must be before the end time '"+end+"'");
		return false;
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
	
	if(check_dt(details) == false) return false;
	
	return true;
}


/// Checks that a number is a positive integer
unsigned int Input::check_pos_integer(string te, unsigned int def)
{
	auto value = get_tag_value(te);
	
	if(value == ""){
		if(def == UNSET) cannot_find_tag();
		return def;
	}
	
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
	
	if(value == ""){
		if(def == UNSET) cannot_find_tag();
		return def;
	}
	
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
	
	if(value == ""){
		if(def == UNSET) cannot_find_tag();
		return def;
	}
	
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
	
	auto ie = get_tag_value("ie"); if(ie == ""){ cannot_find_tag(); return;}
	
	auto spl = split(ie,',');
	vector <IEname> list; 
	for(auto i = 0u; i < spl.size(); i++){
		auto nam = spl[i];
		check_name_input(nam,"Individual effect name",true);
	
		IEname na; na.name = nam; na.index = UNSET;
		list.push_back(na);
	}
	
	auto A = get_tag_value("A"); 
	auto A_sparse = get_tag_value("A-sparse"); 
	auto A_inv = get_tag_value("Ainv"); 
	auto pedigree = get_tag_value("pedigree");
	
	auto num = 0u;
	if(A != "") num++;
	if(A_inv != "") num++;
	if(A_sparse != "") num++;
	if(pedigree != "") num++;
	if(num > 1){
		alert_import("Cannot specify more than one of 'A', 'Ainv', 'A-sparse' or 'pedigree'."); 
		return;
	}
	
	Amatrix A_matrix;
	A_matrix.set = false;

	if(num == 1){
		if(pedigree != ""){
			A_matrix.set = true;
			
			auto tab_ped = load_table(pedigree);
			if(tab_ped.error == true) return;
				
			if(tab_ped.ncol != 3){
				alert_import("The 'pedigree' file should have three columns"); return;
			}
				
			if(tab_ped.heading[0] != "ID" || tab_ped.heading[1] != "sire" || 
					tab_ped.heading[2] != "dam"){
				alert_import("The 'pedigree' file should have columns with headings 'ID', 'sire' and 'dam'"); return;
			}
			
			for(auto r = 0u; r < tab_ped.nrow; r++){
				auto id = tab_ped.ele[r][0];
				auto k = A_matrix.hash_ind_list.find(id);
				if(k != UNSET){
					alert_import("Individual '"+id+"' is referenced more than once"); 
					return;
				}
				
				A_matrix.hash_ind_list.add(A_matrix.ind_list.size(),id);
				A_matrix.ind_list.push_back(id);
			}
			
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
						alert_import("Dam '"+te+"' cannot be found on individual list");
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
					
				default: alert_emsg_input("Not right"); break;
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
			
			for(auto r = 0u; r < tab_ind.nrow; r++){ 
				auto id = tab_ind.ele[r][0];
				auto k = A_matrix.hash_ind_list.find(id);
				if(k != UNSET){
					alert_import("Individual '"+id+"' is referenced more than once"); 
					return;
				}
				
				A_matrix.hash_ind_list.add(A_matrix.ind_list.size(),id);
				A_matrix.ind_list.push_back(id);
			}			
		
			const auto &list = A_matrix.ind_list;
		
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
		}
	
		if(A != ""){
			A_matrix.set = true;
			auto tab = load_table(A);
			if(tab.error == true) return;
		
			for(auto i = 0u; i < tab.heading.size(); i++){
				auto id = tab.heading[i];
				auto k = A_matrix.hash_ind_list.find(id);
				if(k != UNSET){
					alert_import("Individual '"+id+"' is referenced more than once"); 
					return;
				}
				
				A_matrix.hash_ind_list.add(A_matrix.ind_list.size(),id);
				A_matrix.ind_list.push_back(id);
			}
			
			if(tab.nrow != tab.ncol){
				alert_import("The file '"+tab.file+"' must contain an equal number of columns and rows "+print_row_col(tab)+"."); 
				return;
			}
			
			vector < vector <double> > val;
			val.resize(tab.nrow);
			for(auto r = 0u; r < tab.nrow; r++){
				val[r].resize(tab.ncol);
				for(auto c = 0u; c < tab.ncol; c++){
					auto ele = number(tab.ele[r][c]);
					if(ele == UNSET){
						alert_import(in_file_text(tab.file)+ " the element '"+tstr(ele)+"' is not a number (row "+tstr(r+2)+", col "+tstr(c+1)+")");
						return;
					}
					val[r][c] = ele;
				}
			}
			A_matrix.value = val;
		}
		
		if(A_inv != ""){
			A_matrix.set = true;
			auto tab = load_table(A_inv);
			if(tab.error == true) return;
		
			for(auto i = 0u; i < tab.heading.size(); i++){
				auto id = tab.heading[i];
				auto k = A_matrix.hash_ind_list.find(id);
				if(k != UNSET){
					alert_import("Individual '"+id+"' is referenced more than once"); 
					return;
				}
				
				A_matrix.hash_ind_list.add(A_matrix.ind_list.size(),id);
				A_matrix.ind_list.push_back(id);
			}
	
			if(tab.nrow != tab.ncol){
				alert_import("The file '"+tab.file+"' must contain an equal number of columns and rows "+print_row_col(tab)+"."); 
				return;
			}
			
			vector < vector <double> > val;
			val.resize(tab.nrow);
			for(auto r = 0u; r < tab.nrow; r++){
				val[r].resize(tab.ncol);
				for(auto c = 0u; c < tab.ncol; c++){
					auto ele = number(tab.ele[r][c]);
					if(ele == UNSET){
						alert_import(in_file_text(tab.file)+ " the element '"+tstr(ele)+"' is not a number (row "+tstr(r+2)+", col "+tstr(c+1)+")");
						return;
					}
					val[r][c] = ele;
				}
			}
			
			auto val_inv = invert_matrix(val);
			A_matrix.value = val_inv;
		}
		
		if(A_matrix.set != true) emsg("A-matrix should be set");
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

	if(A_matrix.set){
		const auto &A = A_matrix.value;
		auto M = A.size();
		auto num = 0.0;		
		for(auto j = 0u; j < M; j++){
			for(auto i = 0u; i < M; i++){
				if(A[j][i] != 0) num++;
			}
		}
		
		A_matrix.sparsity = num/(M*M);
	}
	else A_matrix.sparsity = UNSET; 
	
	IEgroup ieg; 
	ieg.name = name;
	ieg.list = list; 
	ieg.A_matrix = A_matrix;
	ieg.line_num = line_num;
	ieg.ppc_resample = false;
	ieg.prior_ref = UNSET;
	ieg.log_det_min = UNSET;
	ieg.var_max = UNSET;
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
		
		auto j = X_vector.hash_ind_list.find(id);
		if(j != UNSET){
			alert_import(in_file_text(tab.file)+" individual '"+id+"' has more than one entry (row "+tstr(r+2)+")");
		}
		else{
			X_vector.hash_ind_list.add(X_vector.ind_list.size(),id);
			X_vector.ind_list.push_back(id);
		
			auto val = number(tab.ele[r][1]);
			if(val == UNSET){
				alert_import(in_file_text(tab.file)+" the element '"+tstr(val)+"' is not a number (row "+tstr(r+2)+")");
				return;
			}
			X_vector.value.push_back(val);
		}
	}
	
	auto fix_name = fe_char+"^"+name;
	
	auto th = 0u; while(th < model.param.size() && model.param[th].name != fix_name) th++;
	if(th == model.param.size()){
		alert_import("Parameter '"+fix_name+"' is not specified by a 'param' command");
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


/// Loads inference parameter into the model (for EXT)
void Input::inf_param_command()
{
	auto file = get_tag_value("file");
	auto ch = get_chain();
	
	if(find_in(inf_state_list,ch) != UNSET){
		alert_import("There are multiple 'state-inf' commands for chain "+tstr(ch+1)+"'.");
	}
	inf_param_list.push_back(ch);
	
	// Only loads parameter samples for the output core
	if(!op() || model.mode != EXT) return; 
	
	auto tab = load_table(file); if(tab.error == true) return;

	if(model.param_samp_store.size() == 0){
		model.param_samp_store.resize(model.details.nchain);
	}
	
	model.param_samp_store[ch] = tab;
}


/// Loads inference states into the model (for PPC/EXT)
void Input::inf_state_command()
{
	auto file = get_tag_value("file");
	auto ch = get_chain();

	if(find_in(inf_state_list,ch) != UNSET){
		alert_import("There are multiple 'param-inf' commands for chain "+tstr(ch+1)+"'.");
	}
	inf_state_list.push_back(ch);

	if(model.mode != PPC && model.mode != EXT) return;
		
	// For ext in -core mode only need to load states for that core
	if(model.mode == EXT && mpi.core != ch/model.details.num_per_core) return;
	
	//term_out("load chain"+tstr(ch));
		
	auto i = 0u; while(i < files.size() && files[i].name != file) i++;
	if(i == files.size()){
		alert_import("Could not find the file '"+file+"'");
		return;
	}

	const auto &flines = files[i].lines;	
	
	auto li = 0u;
	
	// Reads in individual key
	vector <string> ind_key;
	{
		string warn = "Problem loading state file";
		
		while(li < flines.size() && trim(flines[li]) != "{") li++;
		if(li == flines.size()) alert_import(warn);
		li++;
		
		auto li_st = li;
		
		while(li < flines.size() && trim(flines[li]) != "}") li++;
		if(li == flines.size()) alert_import(warn);
		
		for(auto i = li_st; i < li; i++){
			auto st = trim(flines[i]);
			if(st.length() > 0 && !begin_str(st,"#")){
				if(begin_str(st,"timepoint")){
				}
				else{
					if(begin_str(st,"spline-out")){
						auto spl = split(st,' ');
						if(spl.size() != 3){ alert_import(warn); return;}
						
						auto name = spl[1];
						auto th = 0u;
						while(th < model.param.size() && model.param[th].name != name) th++;
						if(th == model.param.size()){ alert_import(warn); return;}
						
						auto &par = model.param[th];
						if(!par.time_dep){ alert_import(warn); return;}
						
						if(!par.spline_out.on){
							par.spline_out.on = true;
							par.spline_out.list = split(spl[2],',');
						}
					}
					else{
						auto spl = split(st,':');
						if(spl.size() != 2) alert_import(warn);
						auto num = number(spl[0]);
						auto name = spl[1];
						if(num >= ind_key.size()) ind_key.resize(num+1);
						ind_key[num] = name;
					}
				}
			}
		}
		
		while(li < flines.size() && !begin_str(flines[li],"<<")) li++;
		if(li == flines.size()) alert_import(warn);
	}
	
	vector <string> lines; 
	
	while(li < flines.size()){	
		auto lin = trim(flines[li]);
		if(lin.length() > 2 && lin.substr(0,2) == "<<"){
			if(lines.size() > 0) read_state_sample(ch,lines,ind_key);
			lines.clear();
		}
		lines.push_back(lin);
		li++;
	}

	if(lines.size() > 0) read_state_sample(ch,lines,ind_key);
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
	if(false) cout << file << " " << chain << endl;
}


/// Processes commands from the input file 

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sys/stat.h>

using namespace std;

#include "input.hh"
#include "utils.hh"

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
	if(p == UNSET){ alert_import("To load the data file the species must be set"); return;}
	
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
				if(ds.focal_cl == UNSET){ alert_import("The focal classification '"+focal+"' not recognised"); return;}
			}
			
			auto type = get_tag_value("type");
			ds.init_pop_type = InitPopType(option_error("type",type,{"fixed","dist"},{ INIT_POP_FIXED, INIT_POP_DIST}));
			
			if(ds.init_pop_type == INIT_POP_DIST){
				if(ds.focal_cl == UNSET){
					auto prior = get_tag_value("prior"); 
	
					ds.pop_prior = convert_text_to_prior(prior,line_num);
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
			if(ds.cl == UNSET){ alert_import("In 'class' the value '"+name+"' is not a classification"); return;}
		}
		break;
	
	case COMP_DATA:
		{
			auto name = get_tag_value("class"); if(name == ""){ cannot_find_tag(); return;}
			
			ds.cl = find_cl(p,name);
		
			if(ds.cl == UNSET){ alert_import("In 'class' the value '"+name+"' is not a classification"); return;}
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
			
			if(pos == neg){ alert_import("'pos' and 'neg' cannot both have the same value"); return;}
				
			om.diag_pos = pos; om.diag_neg = neg;
				
			auto comp = get_tag_value("comp"); if(comp == ""){ cannot_find_tag(); return;}
		
			auto spl = split(comp,',');
			
			auto cl = get_cl_from_comp(spl[0],p);
			if(cl == UNSET){ alert_import("Value '"+spl[0]+"' is not a compartment5"); return;}
			
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
			//max(min(a+b×cos(0.06×t)|0.999)|0.001)
			
			auto cl_sel = get_cl_from_trans(name,p);
		
			if(cl_sel == UNSET){
				alert_import("Error with expression '"+name+"'"); return;
			}
		
			ds.cl = cl_sel;
			
			auto obsran = toLower(get_tag_value("obsrange"));
			if(obsran == ""){ cannot_find_tag(); return;}
	
			ds.time_range = TimeRange(option_error("obsrange",obsran,{"all","specify","file"},{ ALL_TIME, SPEC_TIME, FILE_TIME }));
	
			switch(ds.time_range){
			case ALL_TIME: case FILE_TIME: 
				break;
				
			case SPEC_TIME:
				auto start = get_tag_value("start"); if(start == ""){ cannot_find_tag(); return;} 
				ds.time_start = number(start);
				if(ds.time_start == UNSET){ alert_import("'start' must be a number"); return;}
				
				auto end = get_tag_value("end"); if(end == ""){ cannot_find_tag(); return;} 
				ds.time_end = number(end);
				if(ds.time_end == UNSET){ alert_import("'end' must be a number"); return;}
				
				if(ds.time_start >= ds.time_end){ alert_import("'start' must be before 'end'"); return;}
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
				alert_import("Error with expression '"+name+"'"); return;
			}
		
			ds.cl = cl_sel;
		
			load_obs_model(om);
		}
		break;
		
	case GENETIC_DATA:
		{
			auto type = get_tag_value("type"); if(type == "") cannot_find_tag();
			
			ds.gen_data_type = GenDataType(option_error("type",type,{"snp","matrix"},{ SNP_DATA, MATRIX_DATA}));
			
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
	if(type == "") cannot_find_tag(); 
	
	auto sp_type = SpeciesType(option_error("type",type,{"population","individual"},{ POPULATION, INDIVIDUAL}));

	auto trans_tree = false;
	if(sp_type == INDIVIDUAL){
		auto trans_tree_str = toLower(get_tag_value("trans-tree"));
		if(trans_tree_str != ""){
			if(trans_tree_str == "on") trans_tree = true;
			else{
				if(trans_tree_str != "off"){ alert_import("'trans-tree' must be either 'off' or 'on'"); return;}
			}
		}
	}

	if(loop == 2){
		add_species(name,sp_type,trans_tree);
		p_current = model.nspecies-1;		
	}
	else{
		p_current = find_p(name);
		if(p_current == UNSET){ alert_import("Problem finding species '"+name+"'"); return;}
	}
	cl_current = UNSET;
}


/// Adds a classification to the model
void Input::classification_command(unsigned int loop)
{
	auto p = p_current;

	if(p == UNSET){ alert_import("A species must be defined before a classification can be added"); return;}
		
	auto &sp = model.species[p];
	
	auto name = get_tag_value("name");
	if(name == ""){ cannot_find_tag(); return;}
	
	auto clone = get_tag_value("clone");
	if(clone != ""){ // Clones a classification from another species
		if(loop == 2){
			auto p2 = find_p(clone);
			if(p2 == UNSET){
				alert_import("In 'clone' cannot find the species '"+clone+"'");
				return;
			}
			
			auto cl2 = find_cl(p2,name);
			
			if(cl2 == UNSET){
				alert_import("In clone species '"+clone+"' cannot find the classification '"+name+"'");
				return;
			} 
				
			clone_class(p,p2,cl2);
			cl_current = sp.ncla-1;		
		}
		else{
			cl_current = find_cl(p,name);
			if(cl_current == UNSET){ alert_import("Problem finding classification '"+name+"'"); return;}
		}
	}
	else{
		auto index = get_tag_value("index");
		if(index == ""){ cannot_find_tag(); return;}
		if(index.length() != 1){ alert_import("The index '"+index+"' must be just a single character"); return;}
		
		if(index == "t"){ alert_import("The index 't' cannot be used because it is reserved for time variation"); return;}
		
		if(find_in(alphabet,index) == UNSET){ alert_import("Index '"+index+"' must be from the lower case alphabet"); return;}
		
		auto coord = toLower(get_tag_value("coord"));
		
		auto coord_type = CARTESIAN;
		
		if(coord != ""){
			coord_type = Coord(option_error("coord",coord,{"cartesian", "latlng"},{ CARTESIAN, LATLNG}));
		}
		
		if(coord == "latlng"){
			auto def_map = get_tag_value("map");
			if(def_map != ""){
				if(def_map != "load"){ alert_import("'map' must be set to 'load'"); return;}
			}
		}
		
		if(loop == 2){
			add_classification(p,name,index,coord_type);
			cl_current = sp.ncla-1;		
		}
		else{
			cl_current = find_cl(p,name);
			if(cl_current == UNSET){ alert_import("Problem finding classification '"+name+"'"); return;}
		}
	}
}


/// Adds a multiplier to a parameter (used in PPC)
void Input::param_mult_command()
{
	auto full_name = get_tag_value("name"); if(full_name == "") cannot_find_tag();
			
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
	par.auto_value = false;
	par.line_num = line_num;
	
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
	
	par.value.resize(mult);
	par.prior.resize(mult);
	par.parent.resize(mult);
	par.child.resize(mult);
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
		if(p_current == UNSET){ alert_import("Cannot set the species '"+sp+"' as this does not exist"); return;}
	}
	
	auto name = get_tag_value("classification");
	if(name != ""){
		if(p == UNSET){ alert_import("A species must be set before a 'classification' can be set"); return;}
	
		cl_current = find_cl(p,name);
		if(cl_current == UNSET){ alert_import("Cannot set the classification '"+name+"'"); return;}
	}
}


/// Sets a camera in a classification (not used in c++ version)
void Input::camera_command()
{
	get_tag_value("scale");
	get_tag_value("x");
	get_tag_value("y");
	get_tag_value("lat"); 
	get_tag_value("lng");
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
	auto file = get_tag_value("file"); if(file == "") cannot_find_tag();

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
	
	auto name = get_tag_val("name",tags); if(name == ""){ cannot_find_tag(); return;}
	
	/// Checks that a compartment name is valid
	for(auto i = 0u; i < name.length(); i++){
		auto ch = name.substr(i,1);
		if(includes(compnotallow,ch)){
			alert_import("Compartment name '"+name+"' cannot use character '"+ch+"'");
		}
		
		if(str_eq(name,i,sigma)){
			alert_import("Compartment name '"+name+"' cannot use character '"+sigma+"'");
		}
	}
	
	auto color = get_tag_val("color",tags); 
	
	auto fix_str = get_tag_val("fix",tags);
	
	double x = UNSET, y = UNSET, lat = UNSET, lng = UNSET;
	
	switch(claa.coord){
	case CARTESIAN:
		{
			auto x_str = get_tag_val("x",tags);
			if(x_str != ""){
				x = number(x_str); if(x == UNSET){ alert_import("'x' must be a number"); return;}
			}
			
			auto y_str = get_tag_val("y",tags);
			if(y_str != ""){
				y = number(y_str); if(y == UNSET){ alert_import("'y' must be a number"); return;}
			}
			
			if(x != UNSET && y != UNSET){
				for(auto i = 0u; i < claa.ncomp; i++){
					if(claa.comp[i].x == x && claa.comp[i].y == y){		
						alert_import("Compartments '"+name+"' and '"+claa.comp[i].name+"' cannot have the same position");	
						return;
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
				if(lat_str != ""){
					lat = number(lat_str);
					if(lat == UNSET){ alert_import("'lat' must be a number"); return;}
				}
				
				auto lng_str = get_tag_val("lng",tags);
				if(lng_str != ""){
					lng = number(lng_str);
					if(lng == UNSET){ alert_import("'lng' must be a number"); return;}
				}
			
				if(lat != UNSET && lng != UNSET){
					if(check_latlng_error(lat,lng) == true) return;
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
	auto file = get_tag_value("file"); if(file == "") cannot_find_tag();

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
	
	auto te = get_tag_val("name",tags); if(te == "") cannot_find_tag();
	
	if(te.length() < 2){ alert_import("Expression '"+te+"' is not understood"); return;}
	
	auto i = 0u; while(i < te.length()-1 && te.substr(i,2) != "->") i++;
	
	if(i == te.length()-1){ alert_import("Expression '"+te+"' is not understood"); return;}
	
	string fr = trim(te.substr(0,i));
	if(fr == "+") fr = "Source";
	
	string to = trim(te.substr(i+2));
	if(to == "-") to = "Sink";
	
	auto value_str = get_tag_val("value",tags); if(value_str == "") cannot_find_tag();
	
	auto trans_def = extract_trans_def(value_str);
	if(trans_def.set == false){ 
		alert_import("There is a syntax error in value '"+value_str+"'"); return;
	}
	
	unsigned int ci, cf;
	
	if(fr == "Source") ci = UNSET;
	else{
		ci = find_c(p,cl,fr);
		if(ci == UNSET){ alert_import("Cannot find compartment '"+fr+"'"); return;}
	}
	
	if(to == "Sink") cf = UNSET;
	else{
		cf = find_c(p,cl,to);
		if(cf == UNSET){ alert_import("Cannot find compartment '"+to+"'"); return;}
	}

	for(auto tr = 0u; tr < claa.ntra; tr++){
		if(claa.tra[tr].i == ci && claa.tra[tr].f == cf){
			if(ci == UNSET){ alert_import("Source to '"+to+"' already exists"); return;}
			else{
				if(cf == UNSET){ alert_import("Sink from '"+fr+"' already exists"); return;}
				else{ alert_import("Transition from '"+fr+"' to '"+to+"' already exists"); return;}
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
		
	case EXP_RATE_NM:
		{	
			auto eqtype = TRANS_NM_RATE; 
			
			tra.dist_param.push_back(he(add_equation_info(trans_def.rate,eqtype,p,cl)));
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
				return;
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
			tra.dist_param.push_back(he(add_equation_info(trans_def.time,TRANS_MEAN,p,cl)));
		}
		break;
	}
}


/// Sets the data directory
void Input::datadir_command()
{
	datadir = get_tag_value("folder"); if(datadir == ""){ cannot_find_tag(); return;}

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
		if(number(tesize) == UNSET){ alert_import("'textsize' must be a number"); return;}
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
	//auto dist_mat = get_tag_value("distance-matrix"); 

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
	//pt.val = dist_mat; pt.tag = "distance-matrix"; param_tag.push_back(pt);
	
	for(auto j = 0u; j < param_tag.size(); j++){
		for(auto i = j+1; i < param_tag.size(); i++){
			if(param_tag[j].val != "" && param_tag[i].val != ""){
				alert_import("'"+param_tag[j].tag+"' and '"+param_tag[i].tag+"' cannot both be set"); 
				return;
			}
		}
	}
	
	auto j = 0u; while(j < param_tag.size() && param_tag[j].val == "") j++;
	
	auto mult = get_dependency(par.dep,pp,knot_times); if(mult == UNSET) return; 
	
	par.value.resize(mult);
	par.prior.resize(mult);
	par.parent.resize(mult);
	par.child.resize(mult);
	
	par.N = mult;
	
	// Sets default value to zero
	for(auto i = 0u; i < mult; i++) par.value[i].te = "0";
	
	par.trace_output = true;
	if(mult > model.details.param_output_max) par.trace_output = false;
	
	par.auto_value = false;
	
	if(par.name == dist_matrix_name) alert_import("The distance matrix '"+par.full_name+"' must not be set");

	if(value == "auto"){
		par.auto_value = true;
		for(auto &val : par.value) val.te = "auto";
		par.variety = CONST_PARAM;
	}
	else{
		if(pp.dep.size() == 0){
			if(value != ""){
				par.value[0].te = value;
				par.variety = CONST_PARAM;
			}
			
			if(reparam != ""){
				par.value[0] = he(add_equation_info(reparam,REPARAM));
				par.variety = REPARAM_PARAM;
			}
			
			if(cons != ""){
				par.value[0].te = cons;
				par.variety = CONST_PARAM;
			}
		}
		else{
			if(value != "" || cons != "" || reparam != ""){
				par.variety = CONST_PARAM;		
				
				string desc = "In 'value'";
				auto valu = value; 
				if(valu == ""){
					if(cons != ""){
						valu = cons; desc = "In 'const'";
					}
					else{
						if(reparam != ""){
							valu = reparam; desc = "In 'reparam'";
							par.variety = REPARAM_PARAM;		
						}
						else{ alert_import("Problem importing"); return;}
					}
				}
				
				if(reparam != "" && is_file(valu) == false){
					load_reparam_eqn(reparam,par);
				}
				else{
					if(is_file(valu) == false){
						double val = number(valu);
						
						switch(par.variety){
						case CONST_PARAM:
							if(val == UNSET){
								alert_import(desc+" '"+valu+"' is not a number");
								return;
							}
							
							for(auto k = 0u; k < mult; k++){
								par.value[k].te = valu;
							}
							break;
						
						case REPARAM_PARAM:
							{
								if(val == UNSET){
									if(check_eqn_valid(valu) != SUCCESS){
										alert_import(desc+" '"+value+"' is not a valid equation.");
										return;
									}
								}
								
								auto ei = he(add_equation_info(valu,REPARAM));
								for(auto k = 0u; k < mult; k++){
									par.value[k] = ei;
								}
							}
							break;
							
						default: alert_import("Should not be default3"); return;
						}
					}
					else{
						load_param_value(pp,valu,par,desc);
					}
				}
			}
		}
	}

	if(prior != ""){
		par.variety = PRIOR_PARAM;
		
		auto pri = convert_text_to_prior(prior,line_num);
			
		if(pri.error != ""){
			alert_import("Prior syntax error");
			return;
		}
		
		for(auto i = 0u; i < mult; i++) par.prior[i] = pri;
	}
	
	if(prior_split != ""){
		par.variety = PRIOR_PARAM;
		
		if(par.dep.size() == 0){
			alert_import("'prior-split' can only be used if the parameter has a dependency.");  
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
				//const auto &hash_list = par.dep[i].hash_list;			
				//auto vec = hash_list.get_vec_string(subtab.ele[r][i]);
				//ind[i] = hash_list.existing(vec);
				ind[i] = par.dep[i].hash_list.find(subtab.ele[r][i]);
	
				//ind[i] = find_in(par.dep[i].list,subtab.ele[r][i]);
				if(ind[i] == UNSET){ 
					alert_import("The table element '"+subtab.ele[r][i]+"' is not valid (column '"+subtab.heading[i]+"', row "+tstr(r+2)+")");
					return;
				}
			}
			
			auto pri = convert_text_to_prior(subtab.ele[r][ncol-1],line_num);
			
			if(pri.error != ""){
				alert_import("For the table element '"+subtab.ele[r][ncol-1]+"': "+pri.error+" (col '"+subtab.heading[ncol-1]+"', row "+tstr(r+2)+").");
				return;
			}
			
			set_prior_element(par.prior,par.dep,ind,pri);
		}
	}
	
	if(dist != ""){
		par.variety = DIST_PARAM;
		
		auto pri = convert_text_to_prior(dist,line_num);
		
		for(auto i = 0u; i < mult; i++) par.prior[i] = pri;
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
				//const auto &hash_list = par.dep[i].hash_list;			
				//auto vec = hash_list.get_vec_string(subtab.ele[r][i]);
				//ind[i] = hash_list.existing(vec);
				ind[i] = par.dep[i].hash_list.find(subtab.ele[r][i]);
	
				//ind[i] = find_in(par.dep[i].list,subtab.ele[r][i]);
				if(ind[i] == UNSET){ 
					alert_import("The table element '"+subtab.ele[r][i]+"' is not valid (column '"+subtab.heading[i]+"', row "+tstr(r+2)+")");
					return;
				}
			}
			
			auto pri = convert_text_to_prior(subtab.ele[r][ncol-1],line_num);
			
			if(pri.error != ""){
				alert_import("For the table element '"+subtab.ele[r][ncol-1]+"': "+pri.error+" (col '"+subtab.heading[ncol-1]+"', row "+tstr(r+2)+").");
				return;
			}
			
			set_prior_element(par.prior,par.dep,ind,pri);
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
			alert_import("'sim-sample' must be either 'true' or 'false'."); return;
		}
	}
	
	for(const auto &par2 : model.param){
		if(par2.name == par.name){
			alert_import("Parameter already defined"); return;
		}
	}
	
	if(par.variety == UNSET_PARAM){
		if(false){ alert_import("Parameter variety for '"+par.name+"' is unset"); return;}
	}
	else{
		model.param.push_back(par);
	}
}


/// Loads up a reparameterisation
void Input::load_reparam_eqn(string te, Param &par)
{
	const auto &depend = par.dep;
	
	auto eqn_raw = he(add_equation_info(te,REPARAM_EQN));
		
	vector <DepConv> dep_conv;
	for(auto d = 0u; d < depend.size(); d++){
		const auto &dep = depend[d];
		DepConv dc; 
		dc.before = dep.index_with_prime;
		dep_conv.push_back(dc);
	}
	
	auto swap_temp = swap_template(eqn_raw.te,dep_conv);
	if(swap_temp.warn != ""){ alert_import(swap_temp.warn); return;}
		
	for(auto i = 0u; i < par.N; i++){
		for(auto d = 0u; d < depend.size(); d++){
			const auto &dep = depend[d];
			dep_conv[d].after = dep.list[(i/dep.mult)%dep.list.size()];
		}
		
		auto eqn = eqn_raw;
		eqn.te = swap_index_temp(dep_conv,swap_temp);
		
		if(check_swap){
			auto te_st = eqn.te;
			auto res = swap_index(eqn.te,dep_conv);
			if(res.warn != ""){ alert_import(res.warn); return;}
			if(eqn.te != te_st){
				cout << eqn.te << " " << te_st << " compare" << endl; 
				emsg("Swap index dif res");
			}
		}
		
		eqn.type = REPARAM;
		
		auto ind = find_index(i,depend);
		set_reparam_element(par.value,par.dep,ind,eqn);
	}	
}


/// Definition for derived quantities
void Input::derived_command()
{
	auto full_name = get_tag_value("name"); if(full_name == ""){ cannot_find_tag(); return;}
	auto eqn_name = get_tag_value("eqn"); if(eqn_name == ""){ cannot_find_tag(); return;}
	
	auto pp = get_param_prop(full_name);
	if(pp.time_dep == true) pp.dep.pop_back();
	
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
	if(swap_temp.warn != ""){ alert_import(swap_temp.warn); return;}
	
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
			if(res.warn != ""){ alert_import(res.warn); return;}
			if(te_ch != der_eqn.te){
				cout << der_eqn.te << " " << te_ch << " compare" << endl; 
				emsg("Swap index dif res");
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
	
	auto start_num = number(start);
	if(start_num == UNSET){ terminate = true; alert_import("'start' must be a number"); return;}
	details.t_start = start_num;
	
	auto end = get_tag_value("end");
	if(end == ""){ terminate = true; cannot_find_tag(); return;}

	auto end_num = number(end);
	if(end_num == UNSET){ terminate = true; alert_import("'end' must be a number"); return;}
	details.t_end = end_num;
	
	if(details.t_start >= details.t_end){
		alert_import("'start' must before 'end'");
	}

	auto num_str = get_tag_value("number");
	double num = SIM_NUM_DEFAULT;
	if(num_str != "") num = number(num_str);
	if(num == UNSET || num != int(num) || num <= 0){ 
		alert_import("'number' must be a positive integer"); 
		return;
	}
	details.number = (unsigned int)num;
	
	details.seed = get_seed();
	
	auto alg = toLower(get_tag_value("algorithm"));
	if(alg == "") alg = "gillespie";
	
	details.algorithm = Algorithm(option_error("algorithm",alg,{"gillespie","tau"},{ GILLESPIE, TAU }));
	
	auto dt_str = get_tag_value("timestep"); 
	if(dt_str == ""){ terminate = true; cannot_find_tag(); return;}
	
	auto dt = number(dt_str);
	if(dt == UNSET || dt <= 0){
		terminate = true;
		alert_import("'timestep' must be a positive number");
		return;
	}
	
	details.dt = dt;
	details.stochastic = true;
	details.individual_max = check_pos_integer("ind-max",INDMAX_DEFAULT);
	details.param_output_max = check_pos_integer("param-output-max",PARAM_OUTPUT_MAX_DEFAULT);
	details.anneal_type = ANNEAL_NONE;
	details.anneal_rate = UNSET;
	details.anneal_power = UNSET;
}


/// Applies the inference command
void Input::inference_command()
{
	auto &details = model.details;

	auto start = get_tag_value("start"); 
	if(start == ""){ terminate = true; cannot_find_tag(); return;}
	
	auto start_num = number(start);
	if(start_num == UNSET){ terminate = true; alert_import("'start' must be a number"); return;}
	
	details.t_start = start_num;
	
	auto end = get_tag_value("end"); 
	if(end == ""){ terminate = true; cannot_find_tag(); return;}
	
	auto end_num = number(end);
	if(end_num == UNSET){ terminate = true; alert_import("'end' must be a number"); return;}
	
	details.t_end = end_num;
	
	if(details.t_start >= details.t_end){
		alert_import("'start' must before 'end'");
	}
	
	details.seed = get_seed();
	
	details.algorithm = DA_MCMC;
	
	auto alg = toUpper(get_tag_value("algorithm"));
	
	if(alg != ""){
		details.algorithm = Algorithm(option_error("algorithm",alg,{"DA-MCMC","PAS-MCMC","MFA","ABC","ABC-SMC","ABC-MBP","PMCMC","HMC"},{ DA_MCMC, PAS_MCMC, MFA_ALG, ABC_ALG, ABC_SMC_ALG, ABC_MBP, PMCMC, HMC }));
	}
	
	auto dt_str = get_tag_value("timestep"); if(dt_str == ""){ terminate = true; cannot_find_tag(); return;}

	auto dt = number(dt_str);
	if(dt == UNSET || dt <= 0){
		terminate = true; alert_import("'timestep' must be a positive number");
		return;
	}
	
	get_tag_value("nchain"); // nchain is not needed
	
	model.details.dt = dt;
	model.details.sample = MCMC_SAMPLE_DEFAULT;
	model.details.output_param = MCMC_OP_PARAM_DEFAULT;
	model.details.output_state = MCMC_OP_STATE_DEFAULT;

	details.individual_max = check_pos_integer("ind-max",INDMAX_DEFAULT);
	details.param_output_max = check_pos_integer("param-output-max",PARAM_OUTPUT_MAX_DEFAULT);

	switch(details.algorithm){
	case DA_MCMC:
		details.sample = check_pos_integer("sample",MCMC_SAMPLE_DEFAULT);
		details.output_param = check_pos_integer("param-output",MCMC_OP_PARAM_DEFAULT);
		details.output_state = check_pos_integer("state-output",MCMC_OP_STATE_DEFAULT);
		if(num_core() != 1){
			details.output_param /= num_core();
			details.output_state /= num_core();
		}
		else{
			if(nchain != UNSET){
				details.output_param /= nchain;
				details.output_state /= nchain;
			}
		}
		
		if(details.output_param < 1) details.output_param = 1;
		if(details.output_state < 1) details.output_state = 1;
		break;
		
	case ABC_ALG:
		details.sample = check_pos_integer("sample",ABC_SAMPLE_DEFAULT);
		details.accfrac = check_zero_one("acc-frac",ABC_ACFRAC_DEFAULT);
		break;
		
	case ABC_SMC_ALG:
		details.sample = check_pos_integer("sample",ABC_SAMPLE_DEFAULT);
		details.accfrac = check_zero_one("acc-frac",ABCSMC_ACFRAC_DEFAULT);
		details.numgen = check_pos_integer("gen",ABCSMC_GEN_DEFAULT);
		details.kernelsize = check_pos_number("kernel-size",ABCSMC_KERNEL_DEFAULT);
		break;
		
	default: break;
	}
	
	model.details.stochastic = true;

	details.burnin_frac = 20;
	details.anneal_type = ANNEAL_NONE;
	details.anneal_rate = ANNEAL_RATE_DEFAULT;
	details.anneal_power = ANNEAL_POWER_DEFAULT;
	
	if(details.algorithm == DA_MCMC){
		auto burnin_str = get_tag_value("burnin-frac"); 
		
		if(burnin_str != ""){
			auto burnin = number(burnin_str);
			if(burnin == UNSET || burnin < 1 || burnin > 90){
				terminate = true;
				alert_import("'burnin-frac' must be a number between 1 and 90");
				return;		
			}
			details.burnin_frac = burnin;
		}	
		else{
			details.burnin_frac = BURNIN_FRAC_DEFAULT;
		}
	
		auto anneal_str = get_tag_value("anneal"); 
		if(anneal_str != ""){
			details.anneal_type = AnnealType(option_error("anneal",anneal_str,{"none","scan","power-auto","log-auto","power"},{ ANNEAL_NONE, ANNEAL_SCAN, ANNEAL_POWERAUTO, ANNEAL_LOGAUTO, ANNEAL_POWER}));
		
			switch(details.anneal_type){
			case ANNEAL_SCAN:
				{
					auto rate_str = get_tag_value("rate");
					if(rate_str != ""){ 
						auto rate = number(rate_str);
						if(rate == UNSET || rate <= 0){
							terminate = true;
							alert_import("'rate' must be a positive number");
							return;		
						}
					
						details.anneal_rate = rate;
					}
				}
				break;
				
			case ANNEAL_POWER:	
				{
					auto power_str = get_tag_value("power");
					if(power_str != ""){ 
						auto power = number(power_str);
						if(power == UNSET || power <= 0){
							terminate = true;
							alert_import("'power' must be a positive number");
							return;		
						}
					
						details.anneal_power = power;
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
	
	auto start_num = number(start);
	if(start_num == UNSET){ terminate = true; alert_import("'start' must be a number"); return;}
	details.ppc_t_start = start_num;
	
	auto end = get_tag_value("end"); 
	if(end == ""){ terminate = true; cannot_find_tag(); return;}
	
	auto end_num = number(end);
	if(end_num == UNSET){ terminate = true; alert_import("'end' must be a number"); return;}
	details.ppc_t_end = end_num;
	details.inf_t_end = details.t_end;
	if(end_num > details.t_end) details.t_end = end_num;
	
	if(details.t_start >= details.t_end){
		alert_import("'start' must before 'end'");
	}

	auto num_str = get_tag_value("number");
	double num = PPC_NUM_DEFAULT;
	if(num_str != "") num = number(num_str);
	if(num == UNSET || num != int(num) || num <= 0){ 
		alert_import("'number' must be a positive integer"); 
		return;
	}
	
	details.individual_max = check_pos_integer("ind-max",INDMAX_DEFAULT);
	details.param_output_max = check_pos_integer("param-output-max",PARAM_OUTPUT_MAX_DEFAULT);
	
	details.seed = get_seed();

	details.number = (unsigned int) num;
	details.ppc_resample = get_tag_value("resample");
}


/// Checks that a number is a positive integer
unsigned int Input::check_pos_integer(string te, unsigned int def)
{
	auto value = get_tag_value(te);
	
	if(value == "") return def;
	
	auto num = number(value);
	if(num == UNSET || num <= 0 || num != (int)num){
		alert_import("'"+te+"' must be a positive integer");
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
		alert_import("'"+te+"' must be a positive number");
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
		alert_import("'"+te+"' must be between zero and one");
	}
	return num;
}


/// Loads individual effects into the model
void Input::ind_effect_command()
{
	auto p = p_current;
	if(p == UNSET){ alert_import("To add individual effects the species must be set"); return;}
	
	auto name = get_tag_value("name"); if(name == ""){ cannot_find_tag(); return;}
	auto spl = split(name,',');
	vector <IEname> list; 
	for(auto i = 0u; i < spl.size(); i++){
		IEname na; na.name = spl[i]; na.index = UNSET;
		list.push_back(na);
	}
	
	auto A = get_tag_value("A"); 
	auto A_sparse = get_tag_value("A-sparse"); 

	Amatrix A_matrix;
	if(A == ""){
		if(A_sparse == "") A_matrix.set = false;
		else{
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
					alert_import("Problem loading table on line "+to_string(r+1)); return;
				}
				val[j][i] = value;
				val[i][j] = value;
			}
			A_matrix.value = val;
		}
	}
	else{
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
					alert_import(in_file_text(tab.file)+ "the element '"+tstr(ele)+"' is not a number2 (row "+tstr(r+2)+", col "+tstr(c+1)+")");
					return;
				}
				val[r][c] = ele;
			}
		}
		A_matrix.value = val;
	}
	A_matrix.hash_ind_list.create(A_matrix.ind_list);
	
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
	if(p == UNSET){ alert_import("To add individual effects the species must be set"); return;}
	
	auto name = get_tag_value("name"); if(name == ""){ cannot_find_tag(); return;}

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
			alert_import(in_file_text(tab.file)+"the element '"+tstr(val)+"' is not a number3 (row "+tstr(r+2)+")");
			return;
		}
		X_vector.value.push_back(val);
	}
	
	X_vector.hash_ind_list.create(X_vector.ind_list);

	auto fix_name = "ν^"+name;
	
	auto th = 0u; while(th < model.param.size() && model.param[th].name != fix_name) th++;
	if(th == model.param.size()){
		alert_import("Parameter '"+name+"' is not specified by 'param'");
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
	if(i == files.size()) alert_import("Could not find '"+file+"'");

	const auto &flines = files[i].lines;

	auto li = 0u;
	
	// Reads in individual key
	vector <string> ind_key;
	{
		auto warn = "Problem loading state file";
		
		while(li < flines.size() && trim(flines[li]) != "{") li++;
		if(li == flines.size()) alert_import(warn);
		while(li < lines.size() && !begin_str(flines[li],"timepoint")) li++;
		if(li == flines.size()) alert_import(warn);
		li++;
		while(li < flines.size()){
			auto line = trim(flines[li]);
			if(line == "}") break;
			if(line != ""){
				auto spl = split(flines[li],':');
				if(spl.size() != 2) alert_import(warn);
				auto num = number(spl[0]);
				auto name = spl[1];
				if(num >= ind_key.size()) ind_key.resize(num+1);
				ind_key[num] = name;
			}
			li++;
		}
		
		while(li < lines.size() && !begin_str(flines[li],"<<")) li++;
		if(li == flines.size()) alert_import(warn);
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


/// Ignores sample command 
void Input::dummy_file_command()
{
	auto file = get_tag_value("file"); 
	auto chain = get_tag_value("chain"); 
}

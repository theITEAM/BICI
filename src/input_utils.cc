/// Utility functions used in input

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include "json.hpp" // This is a JSON parser from https://github.com/nlohmann/json

using namespace std;

using json = nlohmann::json;

#include "input.hh"
#include "utils.hh"

/// Returns the species number from its name
unsigned int Input::find_p(string name) const 
{
	for(auto p = 0u; p < model.nspecies; p++){
		if(toLower(model.species[p].name) == toLower(name)) return p;
	}
	return UNSET;
}


/// Returns the classification number from its name
unsigned int Input::find_cl(unsigned int p, string name) const 
{
	const auto &sp = model.species[p];
	for(auto cl = 0u; cl < sp.ncla; cl++){
		if(toLower(sp.cla[cl].name) == toLower(name)) return cl;
	}
	return UNSET;
}


/// Returns the classification number from its index
unsigned int Input::find_cl_index(unsigned int p, string index) const
{
	const auto &sp = model.species[p];
	for(auto cl = 0u; cl < sp.ncla; cl++){
		if(sp.cla[cl].index == index) return cl;
	}
	return UNSET;
}


/// Returns the compartment number from its name
unsigned int Input::find_c(unsigned int p, unsigned int cl, string name) const 
{
	//const auto &claa = model.species[p].cla[cl];
	//auto vec = claa.hash_comp.get_vec_string(name);
	//return claa.hash_comp.existing(vec);
	return model.species[p].cla[cl].hash_comp.find(name);
}


/// Checks that it is OK to access claa
bool Input::check_claa_error()
{
	if(p_current == UNSET){
		alert_import("A species needs to be specified before the '"+cline_store.command_name+"' command can be added",true);
		return true;
	}
			
	if(cl_current == UNSET){
		alert_import("A classification need to be specified before the '"+cline_store.command_name+"' command can be added",true);
		return true;
	}
	
	return false;
}


/// Gets the classification object based on store values for p and cl
Classification& Input::get_claa()
{
	if(p_current == UNSET){
		alert_import("A species needs to be specified before the '"+cline_store.command_name+"' command can be added",true);
	}
	
	if(cl_current == UNSET){
		alert_import("A classification needs to be specified before the '"+cline_store.command_name+"' command can be added",true);
	}
	
	return model.species[p_current].cla[cl_current];
}


/// Checks for disallowed characters
bool Input::check_char_allowed(string st, string not_allowed)
{
	for(auto i = 0u; i < st.length(); i++){
		auto ch = st.substr(i,1);
		if(find_in(not_allowed,ch) != UNSET){
			alert_import("In '"+st+"' the character '"+ch+"' is not allowed");
			return false;
		}
	}
	
	return true;
}


/// Checks if a compartment already exists 
bool Input::check_comp_exist(string name, unsigned int p)
{
	if(name == "") return false;
	
	const auto &sp = model.species[p];
	for(auto cl = 0u; cl < sp.ncla; cl++){
		const auto &claa = sp.cla[cl];
		if(claa.name == name){
			alert_import("Cannot have a classification and compartment both with the name '"+name+"'",true);
			return true;
		}
		
		if(claa.index == name){ 
			alert_import("Cannot have a classification index and compartment both with the name '"+name+"'",true);
			return true;
		}
		
		auto c = find_c(p,cl,name);
		if(c != UNSET){
			alert_import("There is already a compartment with the name '"+name+"'",true);
			return true;
		}
	}
	
	return false;
}


/// Displays text referencing a file/table
string Input::in_data_source(const DataSource &ds) const 
{
	return "In table for data source '"+ds.name+"'";
}


/// Displays text referencing a file/table
string Input::in_file_text(string te) const 
{
	if(te.substr(0,7) == "[[$FILE") return "In table '[[...]]'";
	return "In file '"+te+"'";
}


/// Loads a table from a file
Table Input::load_table(const string file)
{
	Table tab; tab.error = false;
	
	auto i = 0u; while(i < files.size() && files[i].name != file) i++;
	
	if(i == files.size()){
		alert_import("Could not find file '"+file+"'");
		tab.error = true; 
		return tab;
	}
	
	const auto &fs = files[i];

	auto sep = fs.sep;
	
	tab.file = fs.name;
	
	auto lines = fs.lines;

	auto j = 0u;
	while(j < lines.size() && (begin_str(trim(lines[j]),"#") || trim(lines[j]) == "")) j++;

	remove_cr(lines[j]);
	
	tab.heading = split(lines[j],sep);

	for(auto &val : tab.heading){
		if(sep == ',') remove_tab(val);
		val = trim(remove_quote(val));
	}
	j++;

	tab.ncol = tab.heading.size();
	while(j < lines.size()){
		remove_cr(lines[j]);
		if(lines[j] != ""){
			auto vec = split(lines[j],sep);
	
			if(vec.size() != tab.ncol){
				alert_import(in_file_text(file)+" there are "+tstr(tab.ncol)+" headings but line "+tstr(j+1)+" ("+trunc(lines[j],30)+") contains "+tstr(vec.size())+" columns.");
				tab.error = true; 
				return tab;
			}

			for(auto &val : vec){
				if(sep == ',') remove_tab(val);
				val = trim(remove_quote(val));
			}
			
			tab.ele.push_back(vec);
		}
		j++;
	}
	
	tab.nrow = tab.ele.size();
	
	return tab;
}


/// Gets a sub-table based on a series of column heading names
Table Input::get_subtable(const Table &tab, const vector <string> &col_name)
{
	Table table;
	table.error = false;
	table.file = tab.file;
	table.ncol = col_name.size();
	table.nrow = tab.nrow;
	table.heading = col_name;
	table.ele.resize(tab.nrow);
	
	vector <unsigned int> col;
	for(auto i = 0u; i < col_name.size(); i++){
		auto cname = col_name[i];
		
		auto c = find_string_in(tab.heading,cname);
		if(c == UNSET){
			alert_import(in_file_text(tab.file)+" cannot find heading '"+cname+"'");
			table.error = true;
			return table;
		}
		col.push_back(c);
	}
	
	for(auto r = 0u; r < table.nrow; r++){
		for(auto i = 0u; i < table.ncol; i++){
			table.ele[r].push_back(tab.ele[r][col[i]]);
		}
	}		
	
	return table;
}
 

/// Gets the classification from the compartment name (or if formated using S|E or S:0.5|E:1) 
unsigned int Input::get_cl_from_comp(string name, unsigned int p) const
{
	const auto &sp = model.species[p];

	auto spl = split_with_bracket(name,'|');
	
	if(spl.size() == 1){	
		auto spl2 = split(name,':');
		for(auto cl = 0u; cl < sp.ncla; cl++){
			const auto &claa = sp.cla[cl];
			for(auto c = 0u; c < claa.ncomp; c++){
				if(claa.comp[c].name == spl2[0]) return cl;
			}
		}
	}
	else{
		auto cl_st = UNSET;
		for(auto i = 0u; i < spl.size(); i++){
			auto cl = get_cl_from_comp(spl[i],p);
			if(cl_st == UNSET) cl_st = cl;
			else{
				if(cl_st != cl) return UNSET;
			}
		}
		return cl_st;
	}
	
	return UNSET;
}


/// Gets the classification from the trans name (or if formated using S->E|E->I or S->E:0.5|E->I:1) 
unsigned int Input::get_cl_from_trans(string name, unsigned int p) const
{
	const auto &sp = model.species[p];
	
	name = replace(name,"->","→");

	auto spl = split_with_bracket(name,'|');
	
	if(spl.size() == 1){	
		auto spl2 = split(name,':');
		for(auto cl = 0u; cl < sp.ncla; cl++){
			const auto &claa = sp.cla[cl];
			for(auto c = 0u; c < claa.ntra; c++){
				if(claa.tra[c].name == spl2[0]) return cl;
			}
		}
	}
	else{
		auto cl_st = UNSET;
		for(auto i = 0u; i < spl.size(); i++){
			auto cl = get_cl_from_trans(spl[i],p);
			if(cl_st == UNSET) cl_st = cl;
			else{
				if(cl_st != cl) return UNSET;
			}
		}
		return cl_st;
	}
	
	return UNSET;
}


/// Finds a value in list of strings (case-independent)
unsigned int Input::find_string_in(const vector <string> &arr, string val) const 
{
	auto i = 0u; while(i < arr.size() && toLower(arr[i]) != toLower(val)) i++;
	if(i == arr.size())	return UNSET;
	return i;
}


/// Sets the columns which need to be loaded 
bool Input::set_loadcol(Command cname, DataSource &ds)
{
	auto &load_col = ds.load_col;

	// If time range of transition observation is needed then this adds extra columns to data
	if(cname == TRANS_DATA){
		if(ds.time_range == FILE_TIME) cname = TRANS_TIMERANGE_DATA; 
	}
	
	auto i = 0u; while(i < data_template.size() && data_template[i].cname != cname) i++;
	
	if(i == data_template.size()){ emsg_input("Cannot find command"); return false;}
	
	auto cols = data_template[i].cols;
	
	for(auto c = 0u; c < cols.size(); c++){
		vector <string> spl = split(cols[c],',');
		
		auto type = ColumnType(option_error("cols",spl[0],{"ID","t","tstart","tend","snp","cl_prob","cl_all","cl_all_prob","from","to","init_pop","add_pop","rem_pop","start","end","result","filt_obspop","filt_obspoptrans","comp_name","start_comp","end_comp",},{ ID_COL, T_COL, TSTART_COL, TEND_COL, SNP_COL, CL_PROB_COL, CL_ALL_COL, CL_ALL_PROB_COL, FROM_COL, TO_COL, INIT_POP_COL, ADD_POP_COL, REM_POP_COL, START_COL, END_COL, RESULT_COL, FILT_OBSPOP_COL, FILT_OBSPOPTRANS_COL, COMP_NAME_COL, START_COMP_COL, END_COMP_COL,  }));
		if(type == UNSET) return false;
		
		if(type == SNP_COL && ds.gen_data_type == MATRIX_DATA) type = GENOBS_COL;
		
		switch(type){
			case ID_COL:
				load_col.push_back(LoadCol("ID","individual ID",TEXT_EL));
				break;

			case T_COL:
				load_col.push_back(LoadCol("t",spl[1],TIME_FLOAT_EL));
				break;
				
			case TSTART_COL:
				load_col.push_back(LoadCol("Start","measurement start time",TIME_FLOAT_EL));
				break;
			
			case TEND_COL:
				load_col.push_back(LoadCol("End","measurement end time",TIME_FLOAT_EL));
				break;

			case SNP_COL:
				load_col.push_back(LoadCol("SNP",spl[1],TEXT_EL));
				break;
				
			case GENOBS_COL:
				load_col.push_back(LoadCol("Obs","Observation ID",TEXT_EL));
				break;

			case CL_PROB_COL:
				{
					auto p = ds.p; if(p == UNSET){ emsg_input("p should be defined"); return false;}
					auto cl = ds.cl; if(cl == UNSET){ emsg_input("cl should be defined"); return false;}
					
					auto name = model.species[p].cla[cl].name; 
						
					load_col.push_back(LoadCol(name,"the classification in '"+name+"'",COMP_PROB_NA_EL,cl));
				}
				break;

			case CL_ALL_COL:
				{
					auto p = ds.p; if(p == UNSET){ emsg_input("p should be defined"); return false;}
						
					const auto &sp = model.species[p];
					
					switch(cname){
					case INIT_POP: 
						if(ds.focal_cl != UNSET){
							load_col.push_back(LoadCol("Compartment","comparment",COMP_ALL_EL));
						}
						else{
							for(auto cl = 0u; cl < sp.ncla; cl++){
								load_col.push_back(LoadCol(sp.cla[cl].name,"'"+sp.cla[cl].name+"'",COMP_EL,cl));
							}
						}
						break;
						
					default:
						for(auto cl = 0u; cl < sp.ncla; cl++){
							load_col.push_back(LoadCol(sp.cla[cl].name,"'"+sp.cla[cl].name+"' into which individuals are added",COMP_EL,cl));
						}
						break;
					}
				}
				break;
				
			case CL_ALL_PROB_COL:
				{
					auto p = ds.p; if(p == UNSET){ emsg_input("p should be defined"); return false;}
						
					const auto &sp = model.species[p];
					for(auto cl = 0u; cl < sp.ncla; cl++){
						load_col.push_back(LoadCol(sp.cla[cl].name,"'"+sp.cla[cl].name+"' into which individuals are added",COMP_PROB_EL,cl));
					}
				}
				break;

			case FROM_COL:
				{
					auto p = ds.p; if(p == UNSET){ emsg_input("p should be defined"); return false;}
					auto cl = ds.cl; if(cl == UNSET){ emsg_input("cl should be defined"); return false;}
					
					auto cl_name = model.species[p].cla[cl].name; 
					load_col.push_back(LoadCol("From","the compartment in '"+cl_name+"' from which individuals come",COMP_SOURCE_SINK_EL,cl));
				}
				break;

			case TO_COL:
				{
					auto p = ds.p; if(p == UNSET){ emsg_input("p should be defined"); return false;}
					auto cl = ds.cl; if(cl == UNSET){ emsg_input("cl should be defined"); return false;}
					
					auto cl_name = model.species[p].cla[cl].name; 
					load_col.push_back(LoadCol(cl_name,"the compartment in '"+cl_name+"' to which individuals go",COMP_SOURCE_SINK_EL,cl));
				}
				break;

			case INIT_POP_COL:
				if(ds.cname != INIT_POP) emsg_input("Should be init pop");
			
				switch(ds.init_pop_type){
				case INIT_POP_FIXED: 
					if(ds.focal_cl != UNSET){
						load_col.push_back(LoadCol("Population","the population",POS_FLOAT_PERCENT_EL));
					}
					else{
						load_col.push_back(LoadCol("Population","the population",POS_ZERO_FLOAT_EL));
					}
					break;
				
				case INIT_POP_DIST:
					if(ds.focal_cl != UNSET){
						load_col.push_back(LoadCol("Distribution","the distribution/α value",PRIOR_DIR_EL));
					}
					else{
						load_col.push_back(LoadCol("Alpha","the alpha value",POS_FLOAT_EL));
					}
					break;
					
				default: emsg_input("Def op wrong"); break;
				}
				break;
		
			case ADD_POP_COL:
				load_col.push_back(LoadCol("Population","the population added",POS_FLOAT_EL));
				break;
				
			case REM_POP_COL:
				load_col.push_back(LoadCol("Population","the population removed",POS_FLOAT_EL));
				break;
			
			case START_COL:
				load_col.push_back(LoadCol("Start","the observation start time",FLOAT_EL));
				break;

			case END_COL:
				load_col.push_back(LoadCol("End","the observation end time",FLOAT_EL));
				break;

			case RESULT_COL:
				load_col.push_back(LoadCol("Result","the disease diagnostic result",DIAG_TEST_EL));
				break;
				
			case FILT_OBSPOP_COL:
				{
					auto p = ds.p; if(p == UNSET){ emsg_input("p should be defined"); return false;}
					auto &sp = model.species[p];
					
					const auto filt = sp.set_comp_filt(ds.filter_str,UNSET,LOWER_BOUND,ds);
					
					for(auto cl = 0u; cl < sp.ncla; cl++){
						if(filt.cla[cl].type == FILE_FILT){
							load_col.push_back(LoadCol(sp.cla[cl].name,"'"+sp.cla[cl].name+"' which specifies sub-population",COMP_PROB_EL,cl));
						}
					}
				
					load_col.push_back(LoadCol("Population","the observed population estimates",POS_ZERO_FLOAT_EL));
					
					switch(ds.obs_model.type){
					case NORMAL_FILE_OBSMOD:
						load_col.push_back(LoadCol("SD","the standard deviation in population estimates",POS_FLOAT_EL));
						break;
						
					case NEGBIN_FILE_OBSMOD:
						load_col.push_back(LoadCol("p","the probability p in negative binomial error",ZERO_ONE_EL));
						break;
						
					default: break;
					}
				}
				break;
				
			case FILT_OBSPOPTRANS_COL:
				{		
					auto p = ds.p; if(p == UNSET){ emsg_input("p should be defined"); return false;}
					auto cl = ds.cl; if(cl == UNSET){ emsg_input("cl should be defined"); return false;}
					auto &sp = model.species[p];
					
					const auto filt = sp.set_comp_filt(ds.filter_str,cl,LOWER_BOUND,ds);
			
					for(auto cl2 = 0u; cl2 < sp.ncla; cl2++){
						if(cl2 != cl){
							if(filt.cla[cl2].type == FILE_FILT){
								load_col.push_back(LoadCol(sp.cla[cl2].name,"'"+sp.cla[cl2].name+"' which specifies sub-population",COMP_PROB_EL,cl2));
							}
						}
					}
					
					load_col.push_back(LoadCol("Number","the observed number of transitions",POS_ZERO_FLOAT_EL));
					
					switch(ds.obs_model.type){
					case NORMAL_FILE_OBSMOD:
						load_col.push_back(LoadCol("SD","the standard deviation in population transition estimates",POS_FLOAT_EL));
						break;
						
					case NEGBIN_FILE_OBSMOD:
						load_col.push_back(LoadCol("p","the probability p in negative binomial error",ZERO_ONE_EL));
						break;
					
					default: break;
					}
				}
				break;
			
			case COMP_NAME_COL:
				load_col.push_back(LoadCol("Name","the compartment name",TEXT_EL));
				break;
				
			case START_COMP_COL:
				load_col.push_back(LoadCol("From","the compartment from which individuals move",TEXT_EL));
				break;
				
			case END_COMP_COL:
				load_col.push_back(LoadCol("To","the compartment to which individuals move",TEXT_EL));
				break;
		}
	}
	
	return true;
}


/// Sets the columns which need to be loaded SNP data
void Input::set_SNP_columns(const Table &tab, DataSource &ds)
{
	auto root = ds.SNP_root;
	auto len = root.length();
	
	ds.load_col.pop_back();
	
	for(auto c = 0u; c < tab.ncol; c++){
		auto head = tab.heading[c];
		if(head.length() >= len && head.substr(0,len) == root){
			ds.load_col.push_back(LoadCol(head,"SNP value",TEXT_EL));
		}
	}
}


/// Sets the columns which need to be loaded in matrix data
void Input::set_genetic_matrix_columns(const Table &tab, DataSource &ds)
{
	auto c = ds.load_col.size()-1;
	for(auto r = 0u; r < tab.ele.size(); r++){
		ds.load_col.push_back(LoadCol(tab.ele[r][c],"Matrix data",POS_INT_EL));
	}
}


/// Loads up the observation model
void Input::load_obs_model(ObsModel &om)
{
	auto error = get_tag_value("error"); if(error == ""){ cannot_find_tag(); return;}
			
	auto spl = split(error,':');

	auto str = toLower(spl[0]);
	if(str == "normal"){
		if(spl.size() != 2) alert_import("For 'error' problem with expression '"+error+"'");
	
		auto val = trim(spl[1]);
		
		if(toLower(val) == "file"){
			om.type = NORMAL_FILE_OBSMOD;
		}
		else{
			if(is_percent(val) == true){
				om.type = NORMAL_PERCENT_OBSMOD;
				om.percent = number(val.substr(0,val.length()-1));
				if(om.percent == UNSET){
					alert_import("For 'error' the expression '"+tstr(om.percent)+"' is not a percentage");
					return;
				}
				
				if(om.percent <= 0){
					alert_import("For 'error' the value '"+tstr(om.percent)+"%' must be a positive percentage");
					return;
				}
			}	
			else{
				om.type = NORMAL_SD_OBSMOD;
				om.sd = number(val);
				if(om.sd == UNSET){
					alert_import("For 'error' the expression '"+val+"' must be a number or a percentage");
					return;
				}
				
				if(om.sd <= 0){
					alert_import("For 'error' the value '"+tstr(om.sd)+"' must be positive");
					return;
				}
			}
		}
	}
	else{
		if(str == "poisson"){
			if(spl.size() != 1) alert_import("For 'error' there is a problem with expression '"+error+"'");
			om.type = POISSON_OBSMOD;
		}
		else{
			if(str == "neg-binomial"){
				if(spl.size() != 2) alert_import("For 'error' there is a problem with expression '"+error+"'");
	
				auto val = trim(spl[1]);
			
				if(toLower(val) == "file"){
					om.type = NEGBIN_FILE_OBSMOD;
				}
				else{
					om.type = NEGBIN_OBSMOD;
					om.p = number(val);
						
					if(om.p == UNSET){
						alert_import("For 'error' the expression '"+val+"' must be a number");
					}
					
					if(om.p <= 0 || om.p >= 1){
						alert_import("For 'error' the value '"+val+"' must be between zero and one");
					}
				}
			}
			else{
				alert_import("For 'error' the observation error '"+str+"' is not recognised. It must be chosen from the following options: 'normal', 'poisson' or 'neg-binomial'.");
			}
		}
	}
}


/// Write an array as a string
string Input::stringify(const vector <string> &arr) const
{
	string te="";
	for(auto i = 0u; i < arr.size(); i++){
		if(i != 0) te += ","; 
		te += arr[i];
	}
	return te;
}


/// Finds a compartment based on a name
CompRef Input::find_comp_from_name(unsigned int p, string te) const 
{
	te = trim(te);
	if(te != ""){ 
		const auto &sp = model.species[p];
		for(auto cl = 0u; cl < sp.ncla; cl++){
			const auto &claa = sp.cla[cl];
			auto c = claa.hash_comp.find(te);
			if(c != UNSET){
				CompRef cr; cr.p = p; cr.cl = cl; cr.c = c; cr.error = "";
				return cr;
			}
		}
	}
	
	CompRef cr; cr.p = UNSET; cr.cl = UNSET; cr.c = UNSET; 
	if(te == "") cr.error = "No compartment specfied";
	else cr.error = "Compartment '"+te+"' not specfied";
	return cr;
}


/// Loads up a geoJSON file
unsigned int Input::import_geojson(string file)
{
	if(check_char_allowed(file,"/<>:\"\\|?*") == false) return UNSET;

	auto i = 0u; while(i < files.size() && files[i].name != file) i++;
	if(i == files.size()) alert_import("Could not find file '"+file+"'",true);
	
	const auto &fs = files[i];
	
	for(auto i = 0u; i < geo_json.size(); i++){
		if(geo_json[i].file == file) return i;
	}
	
	GeoJSON gj; gj.file = file;
	
	stringstream ss;
	for(const auto &li : fs.lines){
		ss << li << endl;
	}
	
	json data = json::parse(ss.str());
	
	for(const auto &ob : data.items()){
		if(ob.key() == "features"){
			auto val = ob.value();
			for(const auto &ob2 : val.items()){
				auto val2 = ob2.value();
				
				vector <string> name;
				for(const auto &ob3 : val2.items()){
					if(ob3.key() == "properties"){
						auto val3 = ob3.value();
						for(const auto &ob4 : val3.items()){
							name.push_back(remove_quote(ob4.value()));
						}	
					}
				}
				
				if(name.size() == 0){ alert_import("Problem loading file '"+file+"'",true); return UNSET;}
					
				for(const auto &ob3 : val2.items()){
					if(ob3.key() == "geometry"){
						auto type =  NO_FEATURE; 
						 
						auto val3 = ob3.value();
						for(const auto &ob4 : val3.items()){
							if(ob4.key() == "type"){
								auto val4 = remove_quote(ob4.value());
								if(val4 == "Polygon") type = POLYGON;
								if(val4 == "MultiPolygon") type = MULTI_POLYGON;
							}
						}
						
						if(type != NO_FEATURE){
							Boundary bound; bound.name = name;
							 
							for(const auto &ob4 : val3.items()){
								if(ob4.key() == "coordinates"){
									auto val4 = ob4.value();
									
									for(const auto &ob5 : val4){
										switch(type){
										case MULTI_POLYGON:
											for(const auto &ob6 : ob5){
												Pgon poly;
												for(const auto &ob7 : ob6){
													if(ob7.size() != 2){ alert_import("Problem loading file '"+file+"'",true); return UNSET;}
													LatLng po; po.lng = ob7[0]; po.lat = ob7[1];
													poly.point.push_back(po);
												}
												bound.polygon.push_back(poly);
											}
											break;
										
										case POLYGON:
											{
												Pgon poly;
												for(const auto &ob6 : ob5){
													if(ob6.size() != 2){ alert_import("Problem loading file '"+file+"'",true); return UNSET;}
													LatLng po; po.lat = ob6[0]; po.lng = ob6[1];
													poly.point.push_back(po);
												}
												bound.polygon.push_back(poly);
											}
											break;
										
										case NO_FEATURE: 
											break;
										}
									}
								}
							}
							
							gj.boundary.push_back(bound);
						}
					}			
				}
			}
		}
	}

	geo_json.push_back(gj);

	return geo_json.size()-1;
}


/// Returns the mean latitude and longitude of a boundary with a given name
LatLng Input::boundary_mean_latlng(unsigned int i, string name)
{
	for(const auto &bo : geo_json[i].boundary){
		if(find_in(bo.name,name) != UNSET){
			auto lat_av = 0.0, lng_av = 0.0, nav = 0.0;
			for(const auto &poly : bo.polygon){
				for(const auto &po : poly.point){
					lat_av += po.lat; lng_av += po.lng; nav++;
				}
			}
			
			LatLng po; po.lat = lat_av/nav; po.lng = lng_av/nav;			
			return po;
		}
	}
	
	alert_import("Problem finding '"+name+"' in file '"+geo_json[i].file+"'",true);
	LatLng po; po.lat = UNSET; po.lng = UNSET;			
	return po;
}


/// Sets the value of a constant based on a dependency and an index
void Input::set_const(Param &par, const vector <unsigned int> &ind, double val)
{
	const auto &dep = par.dep;
	
	auto sum = 0u;
	for(auto i = 0u; i < dep.size(); i++) sum += dep[i].mult*ind[i];
	if(sum >= par.N) emsg_input("Problem setting element");
 
	par.set_cons(sum,val);
}
			
			
/// Sets the text value of an element based on a dependency and an index
void Input::set_element(Param &par, const vector <unsigned int> &ind, string te)
{
	const auto &dep = par.dep;
	
	auto sum = 0u;
	for(auto i = 0u; i < dep.size(); i++) sum += dep[i].mult*ind[i];
	if(sum >= par.N) emsg_input("Problem setting element");
 
	par.set_value_te(sum,te);
}


/// Sets the text value of an element based on a dependency and an index
void Input::set_reparam_element(Param &par, const vector <unsigned int> &ind, const EquationInfo &val)
{
	const auto &dep = par.dep;
	
	auto sum = 0u;
	for(auto i = 0u; i < dep.size(); i++) sum += dep[i].mult*ind[i];
	if(sum >= par.N) emsg_input("Problem setting element");
 
	par.set_value_eqn(sum,val);
}


/// Sets the text value of an element based on a dependency and an index
void Input::set_prior_element(Param &par, const vector <unsigned int> &ind, Prior pri)
{
	const auto &dep = par.dep;
	
	auto sum = 0u;
	for(auto i = 0u; i < dep.size(); i++) sum += dep[i].mult*ind[i];
	if(sum >= par.N) emsg_input("Problem setting element");

	par.set_prior(sum,model.prior.size());
	model.prior.push_back(pri);
}


/// Prints a table to the console
void Input::print_table(const Table &tab) const 
{
	for(auto i = 0u; i < tab.ncol; i++) cout << tab.heading[i] << "\t"; 
	cout << endl;
	auto row_max = tab.nrow; 
	if(row_max > 5) row_max = 5;
	for(auto r = 0u; r < row_max; r++){
		for(auto i = 0u; i < tab.ncol; i++) cout << tab.ele[r][i] << "\t";
		cout << endl;
	}
}


/// Adds a parameter reference to a list (if it doesn't already exist)
void Input::add_to_list(vector <ParamRef> &list, const ParamRef &pr) const 
{
	auto i = 0u; while(i < list.size() && !(list[i].th == pr.th && list[i].index == pr.index)) i++;
	if(i == list.size()) list.push_back(pr);
}


/// Creates a dependency from an vector of indices
unsigned int Input::get_dependency(vector <Dependency> &dep_alter, const ParamProp &pp, const vector <string> &knot_times, const vector <string> &knot_times_out)
{
	dep_alter.clear();
	
	for(auto i = 0u; i < pp.dep.size(); i++){
		auto index = pp.dep[i];
		
		Dependency dep;
		dep.index = index;
		dep.index_with_prime = pp.dep_with_prime[i];
		
		if(index == "t"){
			if(!(pp.time_dep == true && index == "t")){
				emsg_input("Problem with dependency"); 
				return UNSET;
			}
		
			dep.list = knot_times;
			dep.list_out = knot_times_out;
		}
		else{
			auto flag = false;
			for(auto p = 0u; p < model.nspecies; p++){
				const auto &sp = model.species[p];
			
				for(auto cl = 0u; cl < sp.ncla; cl++){
					const auto &claa = sp.cla[cl];
				
					if(claa.index == index){
						for(auto c = 0u; c < claa.ncomp; c++) dep.list.push_back(claa.comp[c].name);
						flag = true;
						break;
					}
				}
				if(flag == true) break;
			}
			
			if(flag != true){ emsg_input("Cannot find the index '"+index+"'"); return UNSET;}
		}
		
		dep.hash_list.create(dep.list);
		dep.hash_list_out.create(dep.list_out);
		/*
		for(auto k = 0u; k < dep.list.size(); k++){
			auto te = dep.list[k];
			auto vec = dep.hash_list.get_vec_string(te);
			dep.hash_list.add(k,vec);
		}
		*/
	
		dep_alter.push_back(dep);
	}
		
	auto mult = 1u;
	for(int i = dep_alter.size()-1; i >= 0; i--){
		dep_alter[i].mult = mult;
		mult *= dep_alter[i].list.size();
	}
	
	return mult;
}


/// Handles any error message
EquationInfo Input::he(EquationInfo eqn_inf, unsigned int lnum)
{
	if(lnum != UNSET) eqn_inf.line_num = lnum;
	else eqn_inf.line_num = line_num;
	if(eqn_inf.error) alert_line(eqn_inf.emsg,eqn_inf.line_num);	
	return eqn_inf;
}	


/// Determines if text string is for a file
bool Input::is_file(string te) const
{
	if(te.length() < 4) return false;
	
	if(te.substr(0,3) == "[[$" && te.substr(te.length()-3,3) == "&]]")	return true;
	
	int k = te.length()-1;
	int kmin = k-10;
	
	while(k > 0 && k > kmin && te.substr(k,1) != ".") k--;
	if(k == 0 || k == kmin) return false;
		
	auto end = toLower(te.substr(k));
	
	if(end != ".csv" && end != ".txt" && end != ".tsv" && end != ".geojson") return false;

	return true;
}


/// Based on i and a dependency works out the index
vector <unsigned int> Input::find_index(unsigned int i, const vector <Dependency> &depend) const
{ 
	vector <unsigned int> ind;
	for(auto d = 0u; d < depend.size(); d++){
		const auto &dep = depend[d];
		ind.push_back((i/dep.mult)%dep.list.size());
	}
	return ind;
}


/// Takes a transition definition and extracts information from it
TransDef Input::extract_trans_def(string value) const
{
	TransDef tdef; tdef.set = false;
	
	value = trim(value);
	if(value == "") return tdef;
	auto spl = split(value,'(');
	if(spl.size() == 0) return tdef; 
	
	if(value.substr(value.length()-1,1) != ")") return tdef;

	auto ty = toLower(spl[0]);
	
	if(ty == "exp"){
		auto rate_temp = get_prop(value,"rate:","end");
		if(rate_temp != ""){
			tdef.rate = rate_temp; tdef.type = EXP_RATE; tdef.set = true;
		}
		else{
			auto mean_temp = get_prop(value,"mean:","end");
			if(mean_temp != ""){
				tdef.mean = mean_temp; tdef.type = EXP_MEAN; tdef.set = true;
			}
		}
	}
	
	if(ty == "gamma"){
		tdef.mean = get_prop(value,"mean:","cv:");
		tdef.cv = get_prop(value,"cv:","end");
		if(tdef.mean != "" && tdef.cv != ""){ tdef.type = GAMMA; tdef.set = true;}
	}
	
	if(ty == "erlang"){
		tdef.mean = get_prop(value,"mean:","shape:");
		tdef.shape = get_prop(value,"shape:","end");
		if(tdef.mean != "" && tdef.shape != ""){ tdef.type = ERLANG; tdef.set = true;}
	}
		
	if(ty == "log-normal"){
		tdef.mean = get_prop(value,"mean:","cv:");
		tdef.cv = get_prop(value,"cv:","end");
		if(tdef.mean != "" && tdef.cv != ""){ tdef.type = LOG_NORMAL; tdef.set = true;}
	}
		
	if(ty == "weibull"){
		tdef.scale = get_prop(value,"scale:","shape:");
		tdef.shape = get_prop(value,"shape:","end");
		if(tdef.scale != "" && tdef.shape != ""){ tdef.type = WEIBULL; tdef.set = true;}
	}
		
	if(ty == "period"){
		tdef.time = get_prop(value,"time:","end");
		if(tdef.time != ""){ tdef.type = PERIOD; tdef.set = true;}
	}
	
	return tdef;
}


/// Tries to get a property from a string
string Input::get_prop(string value, string prop, string end) const
{
	auto i = value.find(prop,0);

	if(i != string::npos){
		i += prop.length();
		if(end == "end"){
			return trim(value.substr(i,value.length()-1-i));
		}
		else{
			auto i2 = value.find(end,0);
			if(i2 != string::npos){
				while(i2 > 0 && value.substr(i2,1) != ",") i2--;
				if(i2 > 0){
					return trim(value.substr(i,i2-i));
				}
			}
		}
	}
	
	return "";
}


/// Generates a warning if the sample cannot be read
void Input::alert_sample(string warn, unsigned int num)
{
	alert_import(warn+" "+to_string(num)); 
}


/// Reads a state file from text
void Input::read_state_sample(const vector <string> &lines, const vector <string> &ind_key)
{
	Sample samp;
	samp.param_value.resize(model.param.size());
	
	samp.species.resize(model.species.size());
	
	auto mode = MODE_NODE;
	string warn = "Problem loading inferred state";
	
	auto p = UNSET;
	for(auto i = 0u; i < lines.size(); i++){
		auto li = lines[i];
		if(li != ""){
			if(li.substr(0,1) == "<"){
				if(li.substr(li.length()-1,1) != ">") alert_sample(warn,1); 
				li = li.substr(1,li.length()-2);
			
				auto spl = split(li,' ');
			
				auto va = spl[0];
				
				if(va == "<STATE"){
					auto spl2 = split(spl[1],'>');
					samp.num = number(spl2[0]);
				}
				
				if(va == "PARAMETERS"){
					mode = MODE_PARAM; 
				}
					
				if(va == "SPECIES"){
					if(spl.size() != 2) alert_sample(warn,2);
					auto sp_name = remove_quote(spl[1]);
					
					p = 0; 
					while(p < model.species.size() && model.species[p].name != sp_name) p++;
					if(p == model.species.size()) alert_sample(warn,3);
				}
				
				if(va == "INITIAL"){
					if(p == UNSET) alert_sample(warn,4);
					if(model.species[p].type == INDIVIDUAL) alert_sample(warn,5);
					mode = MODE_INITIAL;
				}
					
				if(va == "POPULATION"){
					if(p == UNSET) alert_sample(warn,4);
					if(model.species[p].type == INDIVIDUAL) alert_sample(warn,5);
					mode = MODE_POPCHANGE;
				}
					
				if(va == "TRANSITIONS"){ 
					if(p == UNSET) alert_sample(warn,6);
					if(model.species[p].type == INDIVIDUAL) alert_sample(warn,7);
					mode = MODE_TRANS; 
				}
				
				if(va == "TIMEPOINT"){
					mode = MODE_TIMEPOINT;
				}
					
				if(va == "INDIVIDUALS"){
					if(p == UNSET) alert_sample(warn,10);
					if(model.species[p].type == POPULATION) alert_sample(warn,11);
					mode = MODE_INDIVIDUAL;
				}
					
				if(va == "DERIVED"){
					mode = MODE_DERIVE;
				}
					
				if(va == "TRANSTREE"){
					mode = MODE_PHYLO;
				}
				
				if(va == "TRANSDISTPROB"){
					mode = MODE_TRANSDISTPROB;
				}
			}
			else{
				switch(mode){
				case MODE_PHYLO:
					break;
					
				case MODE_PARAM: 
					i = get_param_value(samp.param_value,i,lines,warn);
					break;
					
				case MODE_DERIVE: break;
				
				case MODE_TRANSDISTPROB: break;
				
				case MODE_INITIAL:
					{
						if(p == UNSET) alert_sample(warn,26);
						auto spl = split(lines[i],',');
						
						auto &tab = samp.species[p].cpop_init_tab;
						if(tab.heading.size() == 0){
							tab.heading = spl;
							tab.ncol = spl.size();
						}
						else{
							tab.ele.push_back(spl);
							tab.nrow++;
						}
					}
					break;
					
				case MODE_POPCHANGE:
					break;
					
				case MODE_TRANS:
					{
						if(p == UNSET) alert_sample(warn,26);
						
						auto &tab = samp.species[p].trans_num_tab;
						if(tab.heading.size() == 0){
							tab.heading.push_back("Transitions");
							tab.heading.push_back("Timeline");
							tab.ncol = 2;
						}
						
						auto spl = split(lines[i],':');
						if(spl.size() != 2) alert_sample(warn,27);
						
						spl[0] = replace(spl[0],"->","→");
						tab.ele.push_back(spl);
						tab.nrow++;
					}
					break;
				
				case MODE_TIMEPOINT:
					break;
					
				case MODE_INDIVIDUAL:
					{
						if(p == UNSET) alert_sample(warn,26);
						auto spl = split(lines[i],',');
						
						auto &ind_tab = samp.species[p].ind_tab;
						if(spl.size() > 2 && spl[0] == "index" && spl[1] == "source"){
							ind_tab.heading = spl;
							ind_tab.ncol = spl.size();
						}
						else{
							auto num = number(spl[0]);
							if(num == UNSET) alert_sample(warn,26);
							spl[0] = ind_key[num];
							ind_tab.ele.push_back(spl);
							ind_tab.nrow++;
						}
					}
					break;
					
				default:
					alert_sample(warn,40); 
					return;
				}
			}
		}
	}
	
	if(false){
		for(auto th = 0u; th < model.param.size(); th++){
			cout << model.param[th].name << " ";
			for(auto va : samp.param_value[th]) cout << va << ",";
			cout << "  value" << endl;
		}
	}
	
	model.sample.push_back(samp);
}


/// Gets parameter value from line in the file
unsigned int Input::get_param_value(vector < vector <double> > &param_value, unsigned int i, const vector <string> &lines, string warn)
{
	string not_param = "L^markov,L^non-markov,L^ie,L^dist,L^obs,L^init,Prior";
	auto not_param_list = split(not_param,',');
	
	auto spl = comma_split(lines[i]);
	auto name = remove_quote(replace(spl[0],"->","→"));// remove_quote(spl[0].replace(/->/g,"→"));
	name = remove_escape_char(name);

	if(find_in(not_param_list,name) != UNSET){
		return i;
	}
	
	auto th = 0u;
	while(th < model.param.size() && model.param[th].name != name) th++;
	if(th == model.param.size()){
		alert_sample(warn,13);
		return i;
	}

	const auto &par = model.param[th];
	param_value[th].resize(par.N,UNSET);

	//if(par.used == false) return i;

	switch(spl.size()){
	case 1:           // Set parameter value with dependency
		{
			i++;
			auto title = comma_split(lines[i]);
	
			if(title.size() != par.dep.size()+1) alert_sample(warn,14);
			
			for(auto k = 0u; k < par.dep.size(); k++){
				if(title[k] != par.dep[k].index_with_prime) alert_sample(warn,15);
			}
			if(title[par.dep.size()] != "Value") alert_sample(warn,16);
			
			if(par.spline_out.on){ // Spline has been restricted
				const auto &dp_last = par.dep[par.dep.size()-1];
				const auto &list = dp_last.list;
				const auto &solist = par.spline_out.list;
				
				if(solist.size() > list.size()) alert_sample(warn,17);
				for(auto k = 0u; k < solist.size(); k++){
					if(solist[k] != list[k]) alert_sample(warn,18);
				}

				for(auto j = 0u; j < par.N; j++){
					auto mm = (unsigned int)(j/dp_last.mult)%list.size();
					if(mm < solist.size()){
						i++;
						auto li = comma_split(lines[i]);
						for(auto k = 0u; k < par.dep.size(); k++){
							const auto &dp = par.dep[k];
							auto m = (unsigned int)(j/dp.mult)%dp.list.size();
							if(li[k] != dp.list[m]) alert_sample(warn,17);
						}
						auto value = number(li[par.dep.size()]);
						if(value == UNSET) alert_sample(warn,19);
						param_value[th][j] = value;
					}
					else{
						param_value[th][j] = param_value[th][j-1]; 
					}
				}
			}
			else{
				for(auto j = 0u; j < par.N; j++){
					i++;
					auto li = comma_split(lines[i]);
					for(auto k = 0u; k < par.dep.size(); k++){
						const auto &dp = par.dep[k];
						auto m = (unsigned int)(j/dp.mult)%dp.list.size();
						if(li[k] != dp.list[m]) alert_sample(warn,17);
					}
					auto value = number(li[par.dep.size()]);
					if(value == UNSET) alert_sample(warn,19);
					param_value[th][j] = value;
				}
			}
		}
		break;
		
	case 2:    // Directly set parameter value (univariate)
		{
			auto value = number(spl[1]);
			if(value == UNSET) alert_sample(warn,19);
			param_value[th][0] = value;
		}
		break;
		
	default: alert_sample(warn,20); return i;
	}
	
	return i;
}


/// Loads parameter values from a file 
void Input::load_param_value(const ParamProp &pp, string valu, Param &par, string desc)
{
	auto tab = load_table(valu); if(tab.error == true) return;

	auto col_name = pp.dep_with_prime;

	col_name.push_back("Value");
	
	auto subtab = get_subtable(tab,col_name); if(subtab.error == true) return;
	
	auto ncol = subtab.ncol;
	
	auto ndep = pp.dep.size();
	
	for(auto r = 0u; r < subtab.nrow; r++){
		auto fl = false;
		
		vector <unsigned int> ind(ndep);
		for(auto i = 0u; i < ndep; i++){
			auto ele = subtab.ele[r][i];
			ind[i] = par.dep[i].hash_list.find(ele);
	
			if(ind[i] == UNSET){ 
				fl = true;
				
				if(par.dep[i].hash_list_out.find(ele) == UNSET){	
					alert_import(desc+" the element '"+ele+"' is not valid (column '"+subtab.heading[i]+"', row "+tstr(r+2)+")");
					return;
				}
			}
		}
		
		if(fl == false){
			auto ele = subtab.ele[r][ncol-1];

			double val = number(ele);
			switch(par.variety){
			case CONST_PARAM:
				if(val == UNSET){
					if(par.cat_factor && ele == "*"){
					}
					else{
						alert_import(desc+" the element '"+ele+"' is not a number (column '"+subtab.heading[ncol-1]+"', row "+tstr(r+2)+")");
						return;
					}
				}
							
				set_const(par,ind,val);
				break;
			
			case REPARAM_PARAM:
				if(val == UNSET){
					if(check_eqn_valid(ele) != SUCCESS){
						alert_import(desc+" the element '"+ele+"' is not a valid equation (column '"+subtab.heading[ncol-1]+"', row "+tstr(r+2)+")");
						return;
					}
				}
				set_reparam_element(par,ind,he(add_equation_info(ele,REPARAM)));
				break;
				
			default: emsg_input("Should not be default3"); return;
			}
		}
	}
	
	// Sets default values for any missing
	switch(par.variety){
	case CONST_PARAM:
		{
			auto ref = par.add_cons(0);
			for(auto i = 0u; i < par.N; i++){
				if(par.element_ref[i] == UNSET) par.element_ref[i] = ref;
			}
		}
		break;
		
	case REPARAM_PARAM:
		{
			for(auto i = 0u; i < par.N; i++){
				if(par.element_ref[i] == UNSET) par.set_value_te(i,"0"); 
			}
		}
		break;
		
	default: emsg_input("Should not be default3"); return;
	}
	
	if(par.cat_factor){
		auto num = 0u, numi = 0u;
		auto sum = 0.0;
		auto wsum = 0.0, wi = 0.0;
		for(auto i = 0u; i < par.N; i++){
			auto ref = par.element_ref[i];
			if(ref == UNSET) emsg_input("ref should not be unset"); 
		
			auto val = par.cons[ref];
			auto w = par.weight[i];
			
			if(val == UNSET){
				num++; numi = i; wi = w;
			}
			else{
				if(val < 0){
					emsg_input(desc+" the value '"+tstr(val)+"' must be positive for a factor.");
				}
				sum += w*val;
			}
			wsum += w;
		}
		
		if(num > 1){
			alert_import(desc+" this shouldn't contain more than one element with value '*'");
		}
		
		if(num == 1){
			auto val_new = (wsum-sum)/wi;
			if(val_new < 0){
				alert_import(desc+" the calculated value for '*' is '"+tstr(val_new)+"', which must be positive");
			}
			
			par.set_cons(numi,val_new);
		}
		
		if(num == 0){
			auto mean = sum/wsum;
			if(dif(mean,1,DIF_THRESH)){
				if(!par.cat_factor_weight_on){ 	
					alert_import(desc+" the mean of elements is "+tstr(mean)+" and it should be 1");
				}
				else{
					alert_import(desc+" the weighted mean of elements is "+tstr(mean)+" and it should be 1");
				}
			}
		}
		
		if(false){
			for(auto i = 0u; i < par.N; i++){
				auto ref = par.element_ref[i];
				cout << par.cons[ref] << "  param fact" << endl;
			}
		}
	}
		
	if(false){
		for(auto i = 0u; i < par.N; i++){
			auto ref = par.element_ref[i];
			if(ref == UNSET) emsg_input("ref should not be unset"); 
			
			auto &ele = par.element[ref];
			cout << i << " " << ele.value.te << " val" << endl;
		}
	}
}


/// Loads parameter weight from a file 
void Input::load_weight_value(const ParamProp &pp, string valu, Param &par, string desc)
{
	auto tab = load_table(valu); if(tab.error == true) return;

	auto col_name = pp.dep_with_prime;

	col_name.push_back("Value");
	
	auto subtab = get_subtable(tab,col_name); if(subtab.error == true) return;
	
	auto ncol = subtab.ncol;
	
	auto ndep = pp.dep.size();
	
	par.weight.resize(par.N,1);
	
	for(auto r = 0u; r < subtab.nrow; r++){
		vector <unsigned int> ind(ndep);
		
		auto fl = false;
		for(auto i = 0u; i < ndep; i++){
			auto ele = subtab.ele[r][i];
			
			ind[i] = par.dep[i].hash_list.find(ele);
	
			if(ind[i] == UNSET){ 
				fl = true;
				
				if(par.dep[i].hash_list_out.find(ele) == UNSET){	
					alert_import(desc+" the element '"+ele+"' is not valid (column '"+subtab.heading[i]+"', row "+tstr(r+2)+")");
					return;
				}
			}
		}
		
		if(fl == false){
			auto ele = subtab.ele[r][ncol-1];

			double val = number(ele);
		
			if(val == UNSET){
				alert_import(desc+" the element '"+ele+"' is not a number (column '"+subtab.heading[ncol-1]+"', row "+tstr(r+2)+")");
				return;
			}
			
			auto sum = 0u;
			for(auto i = 0u; i < par.dep.size(); i++) sum += par.dep[i].mult*ind[i];
			if(sum >= par.N) emsg_input("Problem setting element");
	 
			par.weight[sum] = val;
		}
		//par.set_weight(sum,val);
	}
}

	
/// Sets spline based on knot-times and smooth specification
void Input::set_spline(SplineType type, string knot_times_str, string smooth, vector <string> &knot_times, vector <string> &knot_times_out, bool use_inf_time, Param &par)
{
	par.spline_info.on = true;

	par.spline_info.type = type;
	
	vector <double> times;
	
	auto t_start = model.details.t_start;
	if(model.mode == PPC && !use_inf_time) t_start = model.details.ppc_t_start; 
	
	auto t_end = model.details.t_end;
	if(model.mode == PPC && use_inf_time) t_end = model.details.inf_t_end;
		
	knot_times = split(knot_times_str,',');
	for(auto j = 0u; j < knot_times.size(); j++){
		double num;
		
		auto te = knot_times[j];
		if(te == "start"){
			if(j != 0){ alert_import("For 'knot_times' the value 'start' should only occur at the start of the knot definition"); return;}
			num = t_start;
		}
		else{
			if(te == "end"){
				if(j != knot_times.size()-1){ alert_import("For 'knot_times' the value 'end' should only occur at the end of the knot definition"); return;}
				num = t_end;
			}
			else{
				num = number(te);
				if(num == UNSET){
					alert_import("For 'knot_times' the value '"+te+"' must be a number"); return;
				}
				else{
					if(num < t_start || num > t_end){
						par.spline_outside_par = true;
					}
					
					if(num == t_start && j != 0){
						alert_import("For 'knot_times' the end time is set twice"); return;
					}
					
					if(num == t_end && j+1 != knot_times.size()){
						alert_import("For 'knot_times' the end time is set twice"); return;
					}
				}
			}
		}
		
		if(j > 0 && num < times[times.size()-1]){
			alert_import("'knot_times' must be time ordered"); return;
		}
			
		times.push_back(num);
	}

	// Converts to tdiv
	for(auto &ti : times) ti = model.calc_tdiv(ti);
	
	// Removes times which are outside of range
	auto &info = par.spline_info;
	auto T = model.details.T;
	
	auto kmax = times.size()-1;
	while(kmax > 0 && times[kmax] >= T) kmax--;
	if(info.type != SQUARE_SPL) kmax++;
	kmax++;
	
	auto kmin = 1u; while(kmin < times.size() && times[kmin] <= 0) kmin++;
	kmin--;
	
	//tdiv.resize(k+1);
	if(kmin > 0 || kmax < times.size()){
		par.spline_outside = true;
		
		auto times_old = times;
		auto knot_times_old = knot_times;
		times.clear();
		knot_times.clear();
		for(auto k = 0u; k < times_old.size(); k++){
			if(k >= kmin && k < kmax){
				times.push_back(times_old[k]);
				knot_times.push_back(knot_times_old[k]);
			}
			else{
				knot_times_out.push_back(knot_times_old[k]);
			}
		}			
	}
	
	par.spline_info.knot_tdiv = times;
	
	if(smooth == ""){
		par.spline_info.smooth = false;
	}
	else{
		par.spline_info.smooth = true;
		
		auto spl = split(smooth,'(');
		
		auto type = SmoothType(option_error("smooth",spl[0],{"normal","log-normal"},{ NORMAL_SMOOTH, LOG_NORMAL_SMOOTH }));
		if(type == UNSET) return;
			
		par.spline_info.smooth_type = type;

		if(spl.size() != 2){ alert_import("There is syntax error in 'smooth'"); return;}
		
		if(spl[1].substr(spl[1].length()-1,1) != ")"){ alert_import("There is syntax error in 'smooth'"); return;}
		
		auto val_str = spl[1].substr(0,spl[1].length()-1);
		auto val = number(val_str);
		if(val == UNSET){ alert_import("For 'smooth' the value '"+val_str+"' is not a number"); return;}
		if(val <= 0){ alert_import("For 'smooth' the value '"+val_str+"' is not positive"); return;}
		par.spline_info.smooth_value = val;
	}
}


/// Gets the seed value from the input file
unsigned int Input::get_seed()
{
	auto seed_str = get_tag_value("seed");
	if(seed_str != ""){
		auto num = number(seed_str);
		if(num == UNSET){ terminate = true; alert_import("Seed '"+seed_str+"' must be a number"); return UNSET;}
		if(num != int(num)){ terminate = true; alert_import("Seed '"+seed_str+"' must be an integer"); return UNSET;}
		if(num < 0 || num > SEED_MAX){ terminate = true; alert_import("Seed '"+seed_str+"' must be in the range 0-"+to_string(SEED_MAX)); return UNSET;}
		
		return num;
	}
	return UNSET;
}


/// Checks the timestep fits between the start and end times
bool Input::check_dt(const Details &details)
{
	auto num = (details.t_end-details.t_start)/details.dt;
	if(dif(num,double(round_int(num)),DIF_THRESH)){
		alert_import("The difference between the start '"+tstr(details.t_start)+"' and end '"+tstr(details.t_end)+"' times is not a multiple of the time-step '"+tstr(details.dt)+"'.");
		return false;
	}
	
	if(model.calc_tdiv(details.t_end) > TI_DIV_MAX){
		alert_import("The number of time divisions must be fewer than "+tstr(TI_DIV_MAX));
		return false;
	}
	
	return true;
}


/// Adds an extra parameter to model categorical factor
void Input::add_param_cat_factor(Param &par)
{
	Param par_basic = par;
	
	string pre = "basic";
	par_basic.full_name = pre+par_basic.full_name;  
	par_basic.name = pre+par_basic.name;  
	par_basic.trace_output = false;

	auto N = par_basic.N;
	const auto &pri = model.prior[par_basic.get_prior_ref(0)];
	
	auto sigma = number(pri.dist_param[0].te);
	
	//const auto &weight = par.weight;
	
	auto wsum = 0.0;
	for(auto i = 0u; i < N; i++) wsum += par.weight[i];
	
	for(auto i = 0u; i < N; i++){
		auto pri = model.prior[par_basic.get_prior_ref(i)];
		pri.type = DIRICHLET_PR;
		pri.dist_param[0].te = tstr((((N-1)/(sigma*sigma))-1)*(par.weight[i]/wsum));
		par_basic.set_prior(i,model.prior.size());
		model.prior.push_back(pri);
	}
	
	const auto &depend = par.dep;
	par.variety = REPARAM_PARAM;
	
	string denom = "";

	auto name = par_basic.name;
	for(auto i = 0u; i < par.N; i++){
		if(i != 0) denom += "+";
		denom += name+"_";
		for(auto d = 0u; d < depend.size(); d++){
			const auto &dep = depend[d];
			if(d != 0) denom += ",";
			denom += dep.list[(i/dep.mult)%dep.list.size()];
		}
	}
	
	for(auto i = 0u; i < par.N; i++){
		stringstream ss;
		ss << (wsum/par.weight[i]) << "*";
		//ss << "%" << name << "_";
			ss << name << "_";
		for(auto d = 0u; d < depend.size(); d++){
			const auto &dep = depend[d];
			if(d != 0) ss << ",";
			ss << dep.list[(i/dep.mult)%dep.list.size()];
		}
		ss << "/(" << denom << ")";
		par.set_value_eqn(i,he(add_equation_info(ss.str(),REPARAM)));
	}
	
	model.param.push_back(par_basic);
}


/// Loads up a reparameterisation
void Input::add_reparam_eqn(Param &par, Hash &hash_eqn)
{
	const auto &depend = par.dep;
	auto te = par.reparam_eqn;
	par.reparam_eqn = "";
	
	if(par.time_dep){
		if(par.spline_info.type != SQUARE_SPL){
			alert_line("A square spline must be used for time-varying reparameterised parameter '"+par.full_name+"'.",par.line_num);
		}			
	}
	
	auto eqn_raw = he(add_equation_info(te,REPARAM_EQN),par.line_num);
		
	vector <DepConv> dep_conv;
	for(auto d = 0u; d < depend.size(); d++){
		const auto &dep = depend[d];
		DepConv dc; 
		dc.before = dep.index_with_prime;
		dep_conv.push_back(dc);
	}
	
	auto swap_temp = swap_template(eqn_raw.te,dep_conv);
	if(swap_temp.warn != ""){ 
		alert_import(swap_temp.warn); 
		return;
	}
		
	for(auto i = 0u; i < par.N; i++){
		auto ref = par.element_ref[i];
		if(ref != UNSET){
			//auto &ele = par.element[ref];
			
			for(auto d = 0u; d < depend.size(); d++){
				const auto &dep = depend[d];
				dep_conv[d].after = dep.list[(i/dep.mult)%dep.list.size()];
			}
			
			auto eqn = eqn_raw;
			eqn.te = swap_index_temp(dep_conv,swap_temp);
			
			if(check_swap){
				auto te_st = eqn.te;
				auto res = swap_index(eqn.te,dep_conv);
				if(res.warn != ""){
					alert_import(res.warn); 
					return;
				}
				
				if(eqn.te != te_st){
					cout << eqn.te << " " << te_st << " compare" << endl; 
					emsg_input("Swap index dif res");
				}
			}
			
			eqn.type = REPARAM;
			
			model.add_eq_ref(eqn,hash_eqn);
			
			const auto &eqq = model.eqn[eqn.eq_ref];
			if(eqq.time_vari == true){
				if(par.spline_info.type != SQUARE_SPL){
					alert_line("A square spline must be used for time-varying reparameterised expression '"+eqq.te_raw+"'.",par.line_num);	
				}
			}
		
			par.element[ref].value = eqn;
		}
	}
}


/// Gets the data directory
string Input::get_data_dir(string data_dir)
{
	if(!includes(data_dir,":")){
		data_dir = find_root(input_file)+data_dir;
	}
	return data_dir;
}
				
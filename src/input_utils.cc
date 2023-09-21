/// Utility functions used in input

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
//#include <sys/stat.h>
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
	const auto &comp = model.species[p].cla[cl].comp;
	for(auto c = 0u; c < comp.size(); c++){
		if(comp[c].name == name) return c;
	}
	return UNSET;
}


/// Returns the transition number from its name
unsigned int Input::find_tr(unsigned int p, unsigned int cl, string name) const 
{
	const auto &tra = model.species[p].cla[cl].tra;
	for(auto tr = 0u; tr < tra.size(); tr++){
		if(tra[tr].name == name) return tr;
	}
	return UNSET;
}


/*
/// Clones the compatments from one classification to another
void Input::clone_compartments(unsigned int p_from, unsigned int cl_from, unsigned int p_to, unsigned int cl_to)
{
	auto &cla_from = model.species[p_from].cla[cl_from];
	auto &cla_to = model.species[p_to].cla[cl_to];
	cla_to.comp = cla_from.comp;
	cla_to.ncomp = cla_from.ncomp;
	cla_to.tra.resize(0);
	cla_to.ntra = 0;
}
*/


/// Checks that it is OK to access claa
bool Input::check_claa_error()
{
	if(p_current == UNSET){
		alert_import("Population needs to be specified before "+cline_store.command_name);
		return true;
	}
			
	if(cl_current == UNSET){
		alert_import("Classification need to be specified before "+cline_store.command_name);
		return true;
	}
	
	return false;
}


/// Gets the classification object based on store values for p and cl
Classification& Input::get_claa()
{
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
			alert_import("Cannot have a classification and compartment both '"+name+"'");
			return true;
		}
		
		if(claa.index == name){ 
			alert_import("Cannot have a classification index and compartment both '"+name+"'");
			return true;
		}
		
		auto c = find_c(p,cl,name);
		if(c != UNSET){
			alert_import("There is already a compartment with the name '"+name+"'");
			return true;
		}
	}
	
	return false;
}


/// Loads a table from a file
Table Input::load_table(const string file)
{
	Table tab; tab.error = false;
		
	char sep = ',';
	if(file.length() > 4){
		auto end = file.substr(file.length()-4);
		if(end == ".txt" || end == ".tsv") sep = '\t';
	}
	
	if(datadir == ""){ alert_import("'datadir' must be first set"); tab.error = true; return tab;}
	
	auto filefull = datadir+"/"+file;
	ifstream in(filefull);
	if(!in){ alert_import("Cannot open the file '"+filefull+"'."); tab.error = true; return tab;}
	
	tab.file = file;
	
	string line;
	
	do{
		getline(in,line);
	}while(line.substr(0,1) == "#");
		
	if(com_op == false) cout << "Loaded table '" << file << "'." << endl;

	remove_cr(line);

	tab.heading = split(line,sep);
	
	for(auto &val : tab.heading) val = trim(remove_quote(val));

	tab.ncol = tab.heading.size();
	
	do{
		getline(in,line);
		if(in.eof()) break;

		remove_cr(line);

		auto vec = split(line,sep);
	
		if(vec.size() != tab.ncol){
			alert_import("Rows in the file '"+file+"' do not all share the same number of columns.");
			tab.error = true; 
			return tab;
		}

		for(auto &val : vec) val = trim(remove_quote(val));
	
		tab.ele.push_back(vec);
	}while(true);
	tab.nrow = tab.ele.size();
	
	return tab;
}


/// Gets a subtable based on a series of column heading names
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
		auto c = find_string_in(tab.heading,col_name[i]);
		if(c == UNSET){
			alert_import("Cannot find heading '"+col_name[i]+"' in the file '"+tab.file+"'");
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
 

/// Gets the classification from the compartment name
unsigned int Input::get_cl_from_comp(string name, unsigned int p) const
{
	const auto &sp = model.species[p];
	for(auto cl = 0u; cl < sp.ncla; cl++){
		const auto &claa = sp.cla[cl];
		for(auto c = 0u; c < claa.ncomp; c++){
			if(claa.comp[c].name == name) return cl;
		}
	}
	return UNSET;
}


/// Finds a value in list of string case-independent
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
	if(cname == TRANS_DATA || cname == SOURCE_DATA || cname == SINK_DATA){
		if(ds.time_range == FILE_TIME){
			switch(cname){
			case TRANS_DATA: cname = TRANS_TIMERANGE_DATA; break;
			case SOURCE_DATA: cname = SOURCE_TIMERANGE_DATA; break;
			case SINK_DATA: cname = SINK_TIMERANGE_DATA; break;
			default: break;
			}
		}
	}
	
	auto i = 0u; while(i < data_template.size() && data_template[i].cname != cname) i++;
	
	if(i == data_template.size()){ alert_import("Cannot find command"); return false;}
	
	auto cols = data_template[i].cols;
	
	for(auto c = 0u; c < cols.size(); c++){
		vector <string> spl = split(cols[c],',');
		
		auto type = ColumnType(option_error("cols",spl[0],{"ID","t","snp","cl_prob","cl_all","from","to","pop","alpha","start","end","result","filt_obspop","filt_obspoptrans","filt_settraps","setrap_name","comp_name","start_comp","end_comp"},{ ID_COL, T_COL, SNP_COL, CL_PROB_COL, CL_ALL_COL, FROM_COL, TO_COL, POP_COL, ALPHA_COL, START_COL, END_COL, RESULT_COL, FILT_OBSPOP_COL, FILT_OBSPOPTRANS_COL, FILT_SETTRAPS_COL, SETTRAP_NAME_COL,COMP_NAME_COL, START_COMP_COL, END_COMP_COL }));
	
		switch(type){
			case ID_COL:
				load_col.push_back(LoadCol("ID","individual ID",TEXT_EL));
				break;

			case T_COL:
				load_col.push_back(LoadCol("t",spl[1],TIME_FLOAT_EL));
				break;

			case SNP_COL:
				load_col.push_back(LoadCol("SNP",spl[1],TEXT_EL));
				break;

			case CL_PROB_COL:
				{
					auto p = ds.p; if(p == UNSET){ alert_import("p should be defined"); return false;}
					auto cl = ds.cl; if(cl == UNSET){ alert_import("cl should be defined"); return false;}
					
					auto name = model.species[p].cla[cl].name; 
						
					load_col.push_back(LoadCol(name,"the classification in '"+name+"'",COMP_PROB_EL,cl));
				}
				break;

			case CL_ALL_COL:
				{
					if(cname == INIT_POP && ds.focal_cl != UNSET){
						load_col.push_back(LoadCol("Compartment","Comparment",COMP_ALL_EL));
					}
					else{
						auto p = ds.p; if(p == UNSET){ alert_import("p should be defined"); return false;}
						
						const auto &sp = model.species[p];
						for(auto cl = 0u; cl < sp.ncla; cl++){
							load_col.push_back(LoadCol(sp.cla[cl].name,"'"+sp.cla[cl].name+"' into which individuals are added",COMP_EL,cl));
						}
					}
				}
				break;

			case FROM_COL:
				{
					auto p = ds.p; if(p == UNSET){ alert_import("p should be defined"); return false;}
					auto cl = ds.cl; if(cl == UNSET){ alert_import("cl should be defined"); return false;}
					
					auto cl_name = model.species[p].cla[cl].name; 
					load_col.push_back(LoadCol("From","the compartment in '"+cl_name+"' from which individuals come",COMP_SOURCE_SINK_EL,cl));
				}
				break;

			case TO_COL:
				{
					auto p = ds.p; if(p == UNSET){ alert_import("p should be defined"); return false;}
					auto cl = ds.cl; if(cl == UNSET){ alert_import("cl should be defined"); return false;}
					
					auto cl_name = model.species[p].cla[cl].name; 
					load_col.push_back(LoadCol(cl_name,"the compartment in '"+cl_name+"' to which individuals go",COMP_SOURCE_SINK_EL,cl));
				}
				break;

			case POP_COL:
				if(cname == INIT_POP && ds.focal_cl != UNSET){
					load_col.push_back(LoadCol("Population","the population",POS_FLOAT_PERCENT_EL));
				}
				else{
					load_col.push_back(LoadCol("Population","the population",POS_ZERO_FLOAT_EL));
				}
				break;

			case ALPHA_COL:
				load_col.push_back(LoadCol("Alpha","the alpha value",POS_FLOAT_EL));
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
					const auto &filt = ds.comp_filt;
					
					auto p = ds.p; if(p == UNSET){ alert_import("p should be defined"); return false;}

					const auto &sp = model.species[p];
					
					for(auto cl = 0u; cl < sp.ncla; cl++){
						if(filt.cla[cl].type == FILE_FILT){
							load_col.push_back(LoadCol(sp.cla[cl].name,"'"+sp.cla[cl].name+"' which specifies sub-population",MULTI_COMP_EL));
						}
					}
				
					load_col.push_back(LoadCol("Population","the observed population estimates",POS_ZERO_FLOAT_EL));
					
					if(ds.obs_model.type == FILE_OBSMOD){
						load_col.push_back(LoadCol("sd","the standard deviation in population estimates",POS_FLOAT_EL));
					}
				}
				break;
				
			case FILT_OBSPOPTRANS_COL:
				{	
					const auto &filt = ds.comp_filt;
						
					auto p = ds.p; if(p == UNSET){ alert_import("p should be defined"); return false;}
					auto cl = ds.cl; if(cl == UNSET){ alert_import("cl should be defined"); return false;}
					
					const auto &sp = model.species[p];
					for(auto cl2 = 0u; cl2 < sp.ncla; cl2++){
						if(cl2 != cl){
							if(filt.cla[cl2].type == FILE_FILT){
								load_col.push_back(LoadCol(sp.cla[cl2].name,"'"+sp.cla[cl2].name+"' which specifies sub-population",COMP_EL));
							}
						}
					}
					
					load_col.push_back(LoadCol("Number","the observed number of transitions",POS_ZERO_FLOAT_EL));
					
					if(ds.obs_model.type == FILE_OBSMOD){
						load_col.push_back(LoadCol("sd","the standard deviation in population estimates",POS_FLOAT_EL));
					}
				}
				break;
			
			case FILT_SETTRAPS_COL:
				{
					const auto &filt = ds.comp_filt;
					auto p = ds.p; if(p == UNSET){ alert_import("p should be defined"); return false;}
					
					const auto &sp = model.species[p];
					for(auto cl = 0u; cl < sp.ncla; cl++){
						if(filt.cla[cl].type == FILE_FILT){
							load_col.push_back(LoadCol(sp.cla[cl].name,"'"+sp.cla[cl].name+"' which specifies sub-population",COMP_EL));
						}
					}
				}
				break;

			case SETTRAP_NAME_COL:
				load_col.push_back(LoadCol("Name","the unique name for the trapping event",TEXT_EL));
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


/// Loads up the observation model
void Input::load_obs_model(ObsModel &om)
{
	auto error = get_tag_value("error"); if(error == ""){ cannot_find_tag(); return;}
			
	auto spl = split(error,':');
	if(spl.size() != 2){ alert_import("Problem with expression '"+error+"'"); return;}
	
	if(toLower(spl[0]) != "normal"){ alert_import("'"+spl[0]+"' should be 'normal'"); return;}
	
	auto val = trim(spl[1]);
	
	om.percent = 10;
	om.sd = UNSET;
			
	if(toLower(val) == "file"){
		om.type = FILE_OBSMOD;
	}
	else{
		if(is_percent(val) == true){
			om.type = PERCENT_OBSMOD;
			om.percent = number(val.substr(0,val.length()-1));
			if(om.percent == UNSET){
				alert_import("The expression '"+tstr(om.percent)+"' is not a percentage");
				return;
			}
			
			if(om.percent <= 0){
				alert_import("The value '"+tstr(om.percent)+"%' must be a positive percentage");
				return;
			}
		}	
		else{
			om.type = SD_OBSMOD;
			om.sd = number(val);
			if(om.sd == UNSET){
				alert_import("The expression '"+val+"' must be a number or a percentage");
				return;
			}
			
			if(om.sd <= 0){
				alert_import("The value '"+tstr(om.sd)+"' must be positive");
				return;
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
			for(auto c = 0u; c < claa.ncomp; c++){
				if(claa.comp[c].name == te){
					CompRef cr; cr.p = p; cr.cl = cl; cr.c = c; cr.error = "";
					return cr;
				}
			}
		}
	}
	
	CompRef cr; cr.p = UNSET; cr.cl = UNSET; cr.c = UNSET; cr.error = "Compartment '"+te+"' not specfied";
	return cr;
}


/// Loads up a geoJSON file
unsigned int Input::import_geojson(string file)
{
	if(check_char_allowed(file,"/<>:\"\\|?*") == false) return UNSET;
	
	if(datadir == ""){ alert_import("'datadir' must be first set"); return UNSET;}
	
	for(auto i = 0u; i < geo_json.size(); i++){
		if(geo_json[i].file == file) return i;
	}
	
	GeoJSON gj; gj.file = file;
	
	file = datadir+"/"+file;
	
	ifstream f(file);
	if(!f) emsg("File '"+file+"' could not be loaded");
	
	json data = json::parse(f);
	
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
				
				if(name.size() == 0){ alert_import("Problem loading '"+file+"'"); return UNSET;}
					
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
												Polygon poly;
												for(const auto &ob7 : ob6){
													if(ob7.size() != 2){ alert_import("Problem loading '"+file+"'"); return UNSET;}
													LatLng po; po.lng = ob7[0]; po.lat = ob7[1];
													poly.point.push_back(po);
												}
												bound.polygon.push_back(poly);
											}
											break;
										
										case POLYGON:
											{
												Polygon poly;
												for(const auto &ob6 : ob5){
													if(ob6.size() != 2){ alert_import("Problem loading '"+file+"'"); return UNSET;}
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
	
	alert_import("Problem finding '"+name+"' in '"+geo_json[i].file+"'");
	LatLng po; po.lat = UNSET; po.lng = UNSET;			
	return po;
}


/// Sets the text value of an element based on a dependency and an index
void Input::set_element(vector <EquationInfo> &value, const vector <Dependency> &dep, const vector <unsigned int> &ind, string te)
{
	auto sum = 0u;
	for(auto i = 0u; i < dep.size(); i++) sum += dep[i].mult*ind[i];
	
	if(sum >= value.size()) emsg("Problem setting element");
 
	value[sum].te = te; 
}


/// Sets the text value of an element based on a dependency and an index
void Input::set_reparam_element(vector <EquationInfo> &value, const vector <Dependency> &dep, const vector <unsigned int> &ind, const EquationInfo &val)
{
	auto sum = 0u;
	for(auto i = 0u; i < dep.size(); i++) sum += dep[i].mult*ind[i];
	
	if(sum >= value.size()) emsg("Problem setting element");
 
	value[sum] = val; 
}



/// Sets the text value of an element based on a dependency and an index
void Input::set_prior_element(vector <Prior> &prior, const vector <Dependency> &dep, const vector <unsigned int> &ind, Prior pri)
{
	auto sum = 0u;
	for(auto i = 0u; i < dep.size(); i++) sum += dep[i].mult*ind[i];
	
	if(sum >= prior.size()) emsg("Problem setting element");
 
	prior[sum] = pri; 
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
unsigned int Input::get_dependency(vector <Dependency> &dep_alter, const ParamProp &pp, const vector <string> &knot_times)
{
	for(auto i = 0u; i < pp.dep.size(); i++){
		auto index = pp.dep[i];
		
		Dependency dep;
		dep.index = index;
		dep.index_with_prime = pp.dep_with_prime[i];
		
		if(index == "t" || index == "a"){
			if(!(pp.time_dep == true && index == "t") && !(pp.age_dep == true && index == "a")){
				alert_import("Problem with dependency"); 
				return UNSET;
			}
		
			dep.list = knot_times;
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
			
			if(flag != true){ alert_import("Cannot find the index '"+index+"'"); return UNSET;}
		}
		
		dep_alter.push_back(dep);
	}
		
	auto mult = 1u;
	for(int i = dep_alter.size()-1; i >= 0; i--){
		dep_alter[i].mult = mult;
		mult *= dep_alter[i].list.size();
	}
	
	return mult;
}

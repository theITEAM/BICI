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


/// Adds a species to the model
void Input::species_command(unsigned int loop)
{
	auto name = get_tag_value("name");
	if(name == ""){ cannot_find_tag(); return;}
	
	auto type = toLower(get_tag_value("type"));
	if(type == ""){ cannot_find_tag(); return;}
	
	auto sp_type = SpeciesType(option_error("type",type,{"population","individual"},{ POPULATION, INDIVIDUAL}));

	if(loop == 1){
		add_species(name,sp_type);
	
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
		if(loop == 1){
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
		
		if(index == "a"){ alert_import("The index 'a' cannot be used because it is reserved for age variation"); return;}
		
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
		
		if(loop == 1){
			add_classification(p,name,index,coord_type);
			cl_current = sp.ncla-1;		
		}
		else{
			cl_current = find_cl(p,name);
			if(cl_current == UNSET){ alert_import("Problem finding classification '"+name+"'"); return;}
		}
	}
}


/// Sets the classification which is being worked on
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


/// Adds a compartment to the model
void Input::compartment_command()
{
	if(check_claa_error() == true) return;
	const auto &claa = get_claa();
	auto p = p_current;
	auto cl = cl_current;
	
	auto name = get_tag_value("name"); if(name == ""){ cannot_find_tag(); return;}
	
	auto color = get_tag_value("color"); 
	
	auto fix_str = get_tag_value("fix");
	
	double x = UNSET, y = UNSET, lat = UNSET, lng = UNSET;
	
	switch(claa.coord){
	case CARTESIAN:
		{
			auto x_str = get_tag_value("x");
			if(x_str != ""){
				x = number(x_str); if(x == UNSET){ alert_import("'x' must be a number"); return;}
			}
			
			auto y_str = get_tag_value("y");
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
			auto bound_file = get_tag_value("boundary");
			if(bound_file != ""){
				auto i = import_geojson(bound_file); if(i == UNSET) return;
				
				auto po = boundary_mean_latlng(i,name);
			
				lat = po.lat; lng = po.lng;
			}
			else{
				auto lat_str = get_tag_value("lat");
				if(lat_str != ""){
					lat = number(lat_str);
					if(lat == UNSET){ alert_import("'lat' must be a number"); return;}
				}
				
				auto lng_str = get_tag_value("lng");
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

	auto markov_branch = false;
	auto markov_branch_str = get_tag_value("branch-prob"); 
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
	
	add_compartment(name,p,cl,x,y,lat,lng,markov_branch,"");
}


/// Adds a transition / source / sink to the model
void Input::transition_command(const Command cname)
{
	if(check_claa_error() == true) return;
	auto &claa = get_claa();

	auto p = p_current;
	auto cl = cl_current;
	
	auto type = EXP_RATE;
	
	auto type_str = get_tag_value("type");
	if(type_str != ""){
		if(cname == SOURCE_CD){
			type = TransType(option_error("type",type_str,{"exp(rate)","exp(mean)","erlang"},{ EXP_RATE, EXP_MEAN, ERLANG}));
		}
		else{
			type = TransType(option_error("type",type_str,{"exp(rate)","exp(mean)","gamma","erlang","log-normal","weibull"},{ EXP_RATE, EXP_MEAN, GAMMA, ERLANG, LOG_NORMAL, WEIBULL}));
		}
	}
	
	string fr, to;
	if(cname ==  SOURCE_CD){
		fr = "Source"; 
		auto te = get_tag_value("from"); 
		if(te != ""){ alert_import("A source transition should not have a 'from' tag"); return;}
	}
	else{ 
		fr = get_tag_value("from");
		if(fr == ""){ cannot_find_tag(); return;} 
	}
	
	if(cname == SINK_CD){
		to = "Sink"; 
		auto te = get_tag_value("to"); 
		if(te != ""){ alert_import("A sink transition should not have a 'to' tag"); return;}
	}
	else{ 
		to = get_tag_value("to"); 
		if(to == ""){ cannot_find_tag(); return;} 
	}
	
	unsigned int ci, cf;
	
	if(fr == "Source") ci = SOURCE;
	else{
		ci = find_c(p,cl,fr);
		if(ci == UNSET){ alert_import("Cannot find compartment '"+fr+"'"); return;}
	}
	
	if(to == "Sink") cf = SINK;
	else{
		cf = find_c(p,cl,to);
		if(cf == UNSET){ alert_import("Cannot find compartment '"+to+"'"); return;}
	}
	
	for(auto tr = 0u; tr < claa.ntra; tr++){
		if(claa.tra[tr].i == ci && claa.tra[tr].f == cf){
			if(ci == SOURCE){ alert_import("Source to '"+to+"' already exists"); return;}
			else{
				if(cf == SINK){ alert_import("Sink from '"+fr+"' already exists"); return;}
				else{ alert_import("Transition from '"+fr+"' to '"+to+"' already exists"); return;}
			}
		}
	}

	auto x = get_tag_value("x");
	auto y = get_tag_value("y");
				
	auto xmid = get_tag_value("mid-x");
	auto ymid = get_tag_value("mid-y");
			
	auto lng = get_tag_value("lng");
	auto lat = get_tag_value("lat");
	auto lngmid = get_tag_value("mid-lng");
	auto latmid = get_tag_value("mid-lat");
			
	add_transition(p,cl,ci,cf,type);

	auto &tra = claa.tra[claa.ntra-1];
	
	string sup, sup2;
	if(ci ==	SOURCE){ sup = claa.comp[cf].name; sup2 = sup;}
	else{
		if(cf == SINK){ sup = "("+claa.comp[ci].name+"→Sink)"; sup2 = claa.comp[ci].name;}
		else{
			sup = "("+claa.comp[ci].name+"→"+claa.comp[cf].name+")"; sup2 = claa.comp[ci].name;
		}
	}
		
	auto bp = get_tag_value("bp");
	if(bp != "") tra.bp_set = true;
	if(bp == "") bp = "b^"+sup;
	tra.bp = add_equation_info(bp,BP,p,cl);
	
	switch(type){
	case EXP_RATE:
		{
			auto rate = get_tag_value("rate"); 
			if(rate == ""){ if(tra.variety == SOURCE_TRANS) rate = "s^"+sup2; else rate = "r^"+sup;}
			
			auto eqtype = TRANS_RATE; if(tra.variety == SOURCE_TRANS) eqtype = SOURCE_RATE;
			
			tra.dist_param.push_back(add_equation_info(rate,eqtype,p,cl));
		}
		break;
		
	case EXP_MEAN:
		{
			auto mean = get_tag_value("mean"); 
			if(mean == ""){ if(tra.variety == SOURCE_TRANS) mean = "μ^"+sup2; else mean = "μ^"+sup;}
			auto rate = "1/("+mean+")";
			
			tra.type = EXP_RATE;
			auto eqtype = TRANS_RATE; if(tra.variety == SOURCE_TRANS) eqtype = SOURCE_RATE;
			
			tra.dist_param.push_back(add_equation_info(rate,eqtype,p,cl));
		}
		break;
		
	case GAMMA:
		{
			auto mean = get_tag_value("mean"); if(mean == "") mean = "μ^"+sup;
			tra.dist_param.push_back(add_equation_info(mean,TRANS_MEAN,p,cl));
			
			auto cv = get_tag_value("cv"); if(cv == "") cv = "cv^"+sup2;
			tra.dist_param.push_back(add_equation_info(cv,TRANS_CV,p,cl));
		}
		break;
		
	case ERLANG:
		{
			auto mean = get_tag_value("mean"); if(mean == "") mean = "μ^"+sup;
			tra.dist_param.push_back(add_equation_info(mean,TRANS_MEAN,p,cl));
			
			auto shape = get_tag_value("shape"); if(shape == ""){ cannot_find_tag(); return;} 
			
			auto num = number(shape);
			
			if(num == UNSET || num <= 0 || int(num) != num){
				alert_import("For an Erlang distribution the shape parameter must be a positive integer");
				return;
			}
			
			tra.dist_param.push_back(add_equation_info(shape,TRANS_SHAPE,p,cl));
		}
		break;
		
	case LOG_NORMAL:
		{
			auto mean = get_tag_value("mean"); if(mean == "") mean = "μ^"+sup;
			tra.dist_param.push_back(add_equation_info(mean,TRANS_MEAN,p,cl));
			
			auto cv = get_tag_value("cv"); if(cv == "") cv = "cv^"+sup2;
			tra.dist_param.push_back(add_equation_info(cv,TRANS_CV,p,cl));
		}
		break;
		
	case WEIBULL:
		{
			auto scale = get_tag_value("scale"); if(scale == "") scale = "ν^"+sup2;
			tra.dist_param.push_back(add_equation_info(scale,TRANS_SCALE,p,cl));
			
			auto shape = get_tag_value("shape"); if(shape == "") shape = "k^"+sup2;
			tra.dist_param.push_back(add_equation_info(shape,TRANS_SHAPE,p,cl));
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
	
	/// Sets the default output directory 
	if(outputdir == ""){	
		auto i = input_file.length()-1;
		while(i > 0 && input_file.substr(i,1) != ".") i--;
		//outputdir = datadir+"/Output";
		outputdir = "Output";
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
	auto tesize = get_tag_value("textsize"); if(tesize == ""){ cannot_find_tag(); return;}
	if(number(tesize) == UNSET){ alert_import("'textsize' must be a number"); return;}
	
	auto x = get_tag_value("x"); if(x == ""){ cannot_find_tag(); return;}
	auto y = get_tag_value("y"); if(y == ""){ cannot_find_tag(); return;}
	if(number(x) == UNSET){ alert_import("'x' must be a number"); return;}
	if(number(y) == UNSET){ alert_import("'y' must be a number"); return;}
	auto color = get_tag_value("color"); if(color == ""){ cannot_find_tag(); return;}
}


/// Adds a box around compartments in the model (not used in c++ version)
void Input::box_command()
{
	auto te = get_tag_value("text"); if(te == ""){ cannot_find_tag(); return;}
	auto tesize = get_tag_value("textsize"); if(tesize == ""){ cannot_find_tag(); return;}
	if(number(tesize) == UNSET){ alert_import("'textsize' must be a number"); return;}
	
	auto comps = get_tag_value("compartments"); if(comps == ""){ cannot_find_tag(); return;}
	auto color = get_tag_value("color"); if(color == ""){ cannot_find_tag(); return;}
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
	par.age_dep = pp.age_dep;
	par.spline_info.on = false;
	par.used = false;
	par.line_num = line_num;
	
	vector <string> knot_times;
	if(pp.time_dep == true || pp.age_dep == true){
		par.spline_info.on = true;

		if(pp.time_dep == true && pp.age_dep == true){ alert_import("Cannot have time and age dependency"); return;}
		
		auto knot_times_str = get_tag_value("knot-times"); if(knot_times_str == ""){ cannot_find_tag(); return;}
		
		vector <double> times;
		
		knot_times = split(knot_times_str,',');
		for(auto j = 0u; j < knot_times.size(); j++){
			double num;
			
			auto te = knot_times[j];
			if(te == "start"){
				if(j != 0){ alert_import("In 'knot_times' 'start' should only occur at the start of the knot definition"); return;}
				
				num = model.details.t_start;
			}
			else{
				if(te == "end"){
					if(j != knot_times.size()-1){ alert_import("In 'knot_times' 'end' should only occur at the end of the knot definition"); return;}
				
					num = model.details.t_end;
				}
				else{
					num = number(te);
					if(num == UNSET){
						alert_import("In 'knot_times' the value '"+te+"' must be a number"); return;
					}
				}
			}
			
			if(j > 0 && num < times[times.size()-1]){
				alert_import("'knot_times' must be time ordered"); return;
			}
			
			times.push_back(num);
		}
				
		par.spline_info.knot_times = times;
		
		auto smooth = trim(toLower(get_tag_value("smooth")));
		if(smooth == ""){
			par.spline_info.smooth = false;
		}
		else{
			par.spline_info.smooth = true;
			
			auto spl = split(smooth,'(');
			
			auto type = SmoothType(option_error("smooth",spl[0],{"normal","log-normal"},{ NORMAL_SMOOTH, LOG_NORMAL_SMOOTH }));
			par.spline_info.smooth_type = type;
	
			if(spl.size() != 2){ alert_import("There is syntax error in 'smooth'"); return;}
			
			if(spl[1].substr(spl[1].length()-1,1) != ")"){ alert_import("There is syntax error in 'smooth'"); return;}
			
			auto val_str = spl[1].substr(0,spl[1].length()-1);
			auto val = number(val_str);
			if(val == UNSET){ alert_import("In 'smooth' the value '"+val_str+"' is not a number"); return;}
			if(val <= 0){ alert_import("In 'smooth' the value '"+val_str+"' is not positive"); return;}
			par.spline_info.smooth_value = val;
		}
	}
	
	auto cons = get_tag_value("constant"); 
	auto dist = get_tag_value("dist");
	auto dist_split = get_tag_value("dist-split");
	auto value = get_tag_value("value");
	auto reparam = get_tag_value("reparam"); 
	auto prior = get_tag_value("prior"); 
	auto prior_split = get_tag_value("prior-split"); 
	auto dist_mat = get_tag_value("distance-matrix"); 

	vector <ParamTag> param_tag;
	ParamTag pt; 
	pt.val = cons; pt.tag = "constant"; param_tag.push_back(pt);
	pt.val = dist; pt.tag = "dist"; param_tag.push_back(pt);
	pt.val = dist_split; pt.tag = "dist-split"; param_tag.push_back(pt);
	pt.val = value; pt.tag = "value"; param_tag.push_back(pt);
	pt.val = reparam; pt.tag = "reparam"; param_tag.push_back(pt);
	pt.val = prior; pt.tag = "prior"; param_tag.push_back(pt);
	pt.val = prior_split; pt.tag = "prior-split"; param_tag.push_back(pt);
	pt.val = dist_mat; pt.tag = "distance-matrix"; param_tag.push_back(pt);
	
	for(auto j = 0u; j < param_tag.size(); j++){
		for(auto i = j+1; i < param_tag.size(); i++){
			if(param_tag[j].val != "" && param_tag[i].val != ""){
				alert_import("'"+param_tag[j].tag+"' and '"+param_tag[i].tag+"' cannot both be set"); 
				return;
			}
		}
	}
	
	auto j = 0u; while(j < param_tag.size() && param_tag[j].val == "") j++;
	
	if(j == param_tag.size()){
		string te = "One of these possibilities must be set: ";
		for(auto j = 0u; j < param_tag.size(); j++){
			if(j != 0) te += ", ";
			te += "'"+param_tag[j].tag+"'";
		}
		alert_import(te); 
		return;
	}

	auto mult = get_dependency(par.dep,pp,knot_times); if(mult == UNSET) return; 
	
	par.value.resize(mult);
	par.prior.resize(mult);
	par.parent.resize(mult);
	par.child.resize(mult);
	
	par.N = mult;
	
	par.auto_value = false;
	
	if(dist_mat != ""){
		if(dist_mat != "true") alert_import("'distance-matrix' must be set to 'true'");
		
		set_dist(par);
		
		model.param.push_back(par);
		return;
	}
	
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
				par.value[0] = add_equation_info(reparam,REPARAM,UNSET,UNSET);
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
				
				auto tab = load_table(valu); if(tab.error == true) return;
			
				auto col_name = pp.dep_with_prime;
			
				col_name.push_back("Value");
				
				auto subtab = get_subtable(tab,col_name); if(subtab.error == true) return;
			
				print("info",subtab);
				auto ncol = subtab.ncol;
				
				auto ndep = pp.dep.size();
				
				for(auto r = 0u; r < subtab.nrow; r++){
					vector <unsigned int> ind(ndep);
					for(auto i = 0u; i < ndep; i++){
						ind[i] = find_in(par.dep[i].list,subtab.ele[r][i]);
						if(ind[i] == UNSET){ 
							alert_import(desc+" the element '"+subtab.ele[r][i]+"' is not valid (column '"+subtab.heading[i]+"', row "+tstr(r+2)+")");
							return;
						}
					}
					
					auto ele = subtab.ele[r][ncol-1];

					double val = number(ele);
				
					switch(par.variety){
					case CONST_PARAM:
						if(val == UNSET){
							alert_import(desc+" the element '"+ele+"' is not a number (column '"+subtab.heading[ncol-1]+"', row "+tstr(r+2)+")");
							return;
						}
						set_element(par.value,par.dep,ind,ele);
						break;
					
					case REPARAM_PARAM:
						if(val == UNSET){
							if(check_eqn_valid(ele) != SUCCESS){
								alert_import(desc+" the element '"+ele+"' is not a valid equation (column '"+subtab.heading[ncol-1]+"', row "+tstr(r+2)+")");
								return;
							}
						}
						set_reparam_element(par.value,par.dep,ind,add_equation_info(ele,REPARAM,UNSET,UNSET));
						break;
						
					default: alert_import("Should not be default3"); return;
					}
				}
			}
		}
	}

	if(prior != ""){
		par.variety = PRIOR_PARAM;
		Prior pri;
	
		pri.type = PriorPos(option_error("prior",prior,{"uniform","exp","normal","gamma","log-normal","beta","bernoulli","fix","flat","dirichlet"},{UNIFORM_PR,EXP_PR,NORMAL_PR,GAMMA_PR,LOG_NORMAL_PR,BETA_PR,BERNOULLI_PR,FIX_PR,FLAT_PR,DIRICHLET_PR}));
		
		switch(pri.type){
		case FIX_PR:
			{
				auto mean = get_tag_value("fixed"); if(mean == ""){ cannot_find_tag(); return;}
				if(!is_number(mean,"fixed")) return;
				pri.dist_param.push_back(add_equation_info(mean,DIST));
			}
			break;
			
		case UNIFORM_PR:
			{	
				auto min = get_tag_value("min"); if(min == ""){ cannot_find_tag(); return;}
				if(!is_number(min,"min")) return;
				pri.dist_param.push_back(add_equation_info(min,DIST));
				
				auto max = get_tag_value("max"); if(max == ""){ cannot_find_tag(); return;}
				if(!is_number(max,"max")) return;
				pri.dist_param.push_back(add_equation_info(max,DIST));
				
				if(number(min) >= number(max)){
					alert_import("'min' must be smaller than 'max'"); 
					return;
				}
			}
			break;
			
		case NORMAL_PR:
			{
				auto mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				if(!is_number(mean,"mean")) return;
				pri.dist_param.push_back(add_equation_info(mean,DIST));
					
				auto sd = get_tag_value("sd"); if(sd == ""){ cannot_find_tag(); return;}
				if(!is_positive(sd,"sd")) return;
				pri.dist_param.push_back(add_equation_info(sd,DIST));
			}
			break;
			
		case LOG_NORMAL_PR:
			{
				auto mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				if(!is_positive(mean,"mean")) return;
				pri.dist_param.push_back(add_equation_info(mean,DIST));
					
				auto cv = get_tag_value("cv"); if(cv == ""){ cannot_find_tag(); return;}
				if(!is_positive(cv,"cv")) return;
				pri.dist_param.push_back(add_equation_info(cv,DIST));
			}
			break;
			
		case EXP_PR:
			{
				auto mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				if(!is_positive(mean,"mean")) return;
				pri.dist_param.push_back(add_equation_info(mean,DIST));
			}
			break;
			
		case GAMMA_PR:
			{
				auto mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				if(!is_positive(mean,"mean")) return;
				pri.dist_param.push_back(add_equation_info(mean,DIST));
				
				auto cv = get_tag_value("cv"); if(cv == ""){ cannot_find_tag(); return;}
				if(!is_positive(cv,"cv")) return;
				pri.dist_param.push_back(add_equation_info(cv,DIST));
			}
			break;
		
		case BETA_PR:
			{
				auto alpha = get_tag_value("alpha"); if(alpha == ""){ cannot_find_tag(); return;}
				if(!is_positive(alpha,"alpha")) return;
				pri.dist_param.push_back(add_equation_info(alpha,DIST));
				
				auto beta = get_tag_value("beta"); if(beta == ""){ cannot_find_tag(); return;}
				if(!is_positive(beta,"beta")) return;
				pri.dist_param.push_back(add_equation_info(beta,DIST));
			}
			break;
			
		case BERNOULLI_PR:
			{
				auto mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				if(!is_zeroone(mean,"mean")) return;
				pri.dist_param.push_back(add_equation_info(mean,DIST));
			}
			break;
	
		case FLAT_PR:
			break;
		
		case DIRICHLET_PR:
			{
				auto alpha = get_tag_value("alpha"); if(alpha == ""){ cannot_find_tag(); return;}
				if(!is_positive(alpha,"alpha")) return;
				pri.dist_param.push_back(add_equation_info(alpha,DIST));
			}
			break;
	
		default: 
			alert_import("Prior type '"+prior+"' not recognised");
			break;
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
				ind[i] = find_in(par.dep[i].list,subtab.ele[r][i]);
				if(ind[i] == UNSET){ 
					alert_import("The table element '"+subtab.ele[r][i]+"' is not valid (column '"+subtab.heading[i]+"', row "+tstr(r+2)+")");
					return;
				}
			}
			
			auto pri = convert_text_to_prior(subtab.ele[r][ncol-1]);
			
			if(pri.error != ""){
				alert_import("For the table element '"+subtab.ele[r][ncol-1]+"': "+pri.error+" (col '"+subtab.heading[ncol-1]+"', row "+tstr(r+2)+").");
				return;
			}
			
			set_prior_element(par.prior,par.dep,ind,pri);
		}
	}
	
	if(dist != ""){
		par.variety = DIST_PARAM;
		
		Prior pri;
		pri.type = PriorPos(option_error("dist",dist,{"uniform","exp","normal","gamma","log-normal","beta","bernoulli","fix","flat","dirichlet"},{UNIFORM_PR,EXP_PR,NORMAL_PR,GAMMA_PR,LOG_NORMAL_PR,BETA_PR,BERNOULLI_PR,FIX_PR,FLAT_PR,DIRICHLET_PR}));
		
		switch(pri.type){
		case FIX_PR: alert_import("'dist' cannot be fixed"); return;
		
		case UNIFORM_PR:
			{	
				auto min = get_tag_value("min"); if(min == ""){ cannot_find_tag(); return;}
				pri.dist_param.push_back(add_equation_info(min,DIST));
				
				auto max = get_tag_value("max"); if(max == ""){ cannot_find_tag(); return;}
				pri.dist_param.push_back(add_equation_info(max,DIST));
			}
			break;
			
		case NORMAL_PR:
			{
				auto mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				pri.dist_param.push_back(add_equation_info(mean,DIST));
					
				auto sd = get_tag_value("sd"); if(sd == ""){ cannot_find_tag(); return;}
				pri.dist_param.push_back(add_equation_info(sd,DIST));
			}
			break;
			
		case LOG_NORMAL_PR:
			{
				auto mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				pri.dist_param.push_back(add_equation_info(mean,DIST));
					
				auto cv = get_tag_value("cv"); if(cv == ""){ cannot_find_tag(); return;}
				pri.dist_param.push_back(add_equation_info(cv,DIST));
			}
			break;
			
		case EXP_PR:
			{
				auto mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				pri.dist_param.push_back(add_equation_info(mean,DIST));
			}
			break;
			
		case GAMMA_PR:
			{
				auto mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				pri.dist_param.push_back(add_equation_info(mean,DIST));
				
				auto cv = get_tag_value("cv"); if(cv == ""){ cannot_find_tag(); return;}
				pri.dist_param.push_back(add_equation_info(cv,DIST));
			}
			break;
		
		case BETA_PR:
			{
				auto alpha = get_tag_value("alpha"); if(alpha == ""){ cannot_find_tag(); return;}
				pri.dist_param.push_back(add_equation_info(alpha,DIST));
				
				auto beta = get_tag_value("beta"); if(beta == ""){ cannot_find_tag(); return;}
				pri.dist_param.push_back(add_equation_info(beta,DIST));
			}
			break;
			
		case BERNOULLI_PR:
			{
				auto mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				pri.dist_param.push_back(add_equation_info(mean,DIST));
			}
			break;
			
		case FLAT_PR:
			break;
		
		case DIRICHLET_PR:
			auto alpha = get_tag_value("alpha"); if(alpha == ""){ cannot_find_tag(); return;}
			pri.dist_param.push_back(add_equation_info(alpha,DIST));
			break;
		}
		
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
				ind[i] = find_in(par.dep[i].list,subtab.ele[r][i]);
				if(ind[i] == UNSET){ 
					alert_import("The table element '"+subtab.ele[r][i]+"' is not valid (column '"+subtab.heading[i]+"', row "+tstr(r+2)+")");
					return;
				}
			}
			
			auto pri = convert_text_to_prior(subtab.ele[r][ncol-1]);
			
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
		alert_import("Parameter variety for '"+par.name+"' is unset"); return;
	}
	
	model.param.push_back(par);
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
	
	for(auto i = 0u; i < mult; i++){
		vector <DepConv> dep_conv;
		for(auto d = 0u; d < depend.size(); d++){
			const auto &dep = depend[d];
			
			DepConv dc; 
			dc.before = dep.index_with_prime;
			dc.after = dep.list[(i/dep.mult)%dep.list.size()];
			dep_conv.push_back(dc);
		}
		
		
		auto der_eqn = add_equation_info(eqn_name,DERIVE_EQN);
		auto res = swap_index(der_eqn.te,dep_conv);
		if(res.warn != ""){ alert_import(res.warn); return;}
		
		der.eq.push_back(der_eqn);
	}	
						
	model.derive.push_back(der);
}


/// Applies the simulation command
void Input::simulation_command()
{
	if(model.mode != MODE_UNSET){ alert_import("Cannot set 'simulation' and 'inference'"); return;}
	model.mode = SIM;

	auto &details = model.details;
	
	auto start = get_tag_value("start"); if(start == ""){ terminate = true; cannot_find_tag(); return;}
	auto start_num = number(start);
	if(start_num == UNSET){ terminate = true; alert_import("'start' must be a number"); return;}
	details.t_start = start_num;
	
	auto end = get_tag_value("end"); if(end == ""){ terminate = true; cannot_find_tag(); return;}
	auto end_num = number(end);
	if(end_num == UNSET){ terminate = true; alert_import("'end' must be a number"); return;}
	details.t_end = end_num;
	
	if(details.t_start >= details.t_end){
		alert_import("'start' must before 'end'");
	}
	
	auto num_str = get_tag_value("number");
	if(num_str == "") num_str = "1";
	auto num = number(num_str);
	if(num == UNSET || num != int(num) || num <= 0){ 
		alert_import("'number' must be a positive integer"); 
		return;
	}
	details.number = (unsigned int)num;
	
	auto alg = toLower(get_tag_value("algorithm"));
	if(alg == "") alg = "gillespie";
	
	details.algorithm = Algorithm(option_error("algorithm",alg,{"gillespie","tau"},{ GILLESPIE, TAU }));
	
	auto dt_str = get_tag_value("timestep"); if(dt_str == ""){ terminate = true; cannot_find_tag(); return;}
	auto dt = number(dt_str);
	if(dt == UNSET || dt <= 0){
		terminate = true;
		alert_import("'timestep' must be a positive number");
		return;
	}
	
	details.dt = dt;
	
	details.stochastic = true;
	
	details.individual_max = check_pos_integer("ind-max",10000);
}


/// Applies the simulation command
void Input::inference_command()
{
	if(model.mode != MODE_UNSET){ alert_import("Cannot set 'simulation' and 'inference'"); return;}
	model.mode = INF;
	
	auto &details = model.details;
	
	auto start = get_tag_value("start"); if(start == ""){ terminate = true; cannot_find_tag(); return;}
	auto start_num = number(start);
	if(start_num == UNSET){ terminate = true; alert_import("'start' must be a number"); return;}
	details.t_start = start_num;
	
	auto end = get_tag_value("end"); if(end == ""){ terminate = true; cannot_find_tag(); return;}
	auto end_num = number(end);
	if(end_num == UNSET){ terminate = true; alert_import("'end' must be a number"); return;}
	details.t_end = end_num;
	
	if(details.t_start >= details.t_end){
		alert_import("'start' must before 'end'");
	}
	
	auto alg = toUpper(get_tag_value("algorithm"));
	
	details.algorithm = Algorithm(option_error("algorithm",alg,{"DA-MCMC","MFA","ABC","ABC-SMC","ABC-MBP","ABC_PAS","PMCMC","HMC"},{ DA_MCMC, MFA, ABC_ALG, ABC_SMC_ALG, ABC_MBP, ABC_PAS, PMCMC, HMC }));
	
	auto dt_str = get_tag_value("timestep"); if(dt_str == ""){ terminate = true; cannot_find_tag(); return;}

	auto dt = number(dt_str);
	if(dt == UNSET || dt <= 0){
		terminate = true; alert_import("'timestep' must be a positive number");
		return;
	}
	model.details.dt = dt;
	
	model.details.sample = 10000;
	model.details.thinparam = 10;
	model.details.thinstate = 100;

	details.individual_max = check_pos_integer("ind-max",10000);

	switch(details.algorithm){
	case DA_MCMC:
		details.sample = check_pos_integer("sample",10000);
		details.thinparam = check_pos_integer("param-thin",10);
		details.thinstate = check_pos_integer("state-thin",100);
		break;
		
	case ABC_ALG:
		details.sample = check_pos_integer("sample",1000);
		details.accfrac = check_zero_one("acc-frac",0.1);
		break;
		
	case ABC_SMC_ALG:
		details.sample = check_pos_integer("sample",1000);
		details.accfrac = check_zero_one("acc-frac",0.5);
		details.numgen = check_pos_integer("gen",5);
		details.kernelsize = check_pos_number("kernel-size",0.5);
		break;
		
	default: break;
	}
	
	model.details.stochastic = true;
}


/// Check that a number is a positive integer
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


/// Check that a number is positive 
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


/// Check that a number is between zero and one
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
	
	auto  A = get_tag_value("A"); 
	
	Amatrix A_matrix;
	if(A == ""){
		A_matrix.set = false;
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
					alert_import("In file '"+tab.file+"' the element '"+tstr(ele)+"' is not a number (row "+tstr(r+2)+", col "+tstr(c+1)+")");
					return;
				}
				val[r][c] = ele;
			}
		}
		A_matrix.value = val;
	}
	
	IEgroup ieg; 
	ieg.list = list; 
	ieg.A_matrix = A_matrix; 
	ieg.line_num = line_num;
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
		X_vector.ind_list.push_back(subtab.ele[r][0]);
		auto val = number(tab.ele[r][1]);
		if(val == UNSET){
			alert_import("In file '"+tab.file+"' the element '"+tstr(val)+"' is not a number (row "+tstr(r+2)+")");
			return;
		}
		X_vector.value.push_back(val);
	}

	auto th = 0u; while(th < model.param.size() && model.param[th].name != name) th++;
	if(th == model.param.size()){
		alert_import("Parmeter '"+name+"' is not specified by 'param'");
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


/// Imports a data-table 
void Input::import_data_table_command(Command cname)
{
	switch(cname){
	case ADD_IND_SIM: cname = ADD_IND; break;
	case REMOVE_IND_SIM: cname = REMOVE_IND; break;
	case MOVE_IND_SIM: cname = MOVE_IND; break;
	case INIT_POP_SIM: cname = INIT_POP; break;
	default: break;
	}
	
	auto file = get_tag_value("file"); if(file == ""){ cannot_find_tag(); return;}
	
	auto tab = load_table(file); if(tab.error == true) return;

	auto cols = get_tag_value("cols");
	
	auto p = p_current;
	if(p == UNSET){ alert_import("To load the data file the species must be set"); return;}
	
	auto &sp = model.species[p];
	
	DataSource ds;
	ds.line_num = line_num;
	ds.p = p; ds.cl = UNSET;
	ds.tr = UNSET;
	ds.cname = cname; ds.focal_cl = UNSET;
	ds.time_range = ALL_TIME; ds.time_start = UNSET; ds.time_end = UNSET;

	auto &om = ds.obs_model;

	switch(cname){  // Modifies specification dependent on data source
	case INIT_POP:
		{
			auto focal = get_tag_value("focal");
			if(focal != ""){
				ds.focal_cl = find_cl(p,focal); 
				if(ds.focal_cl == UNSET){ alert_import("The focal classification '"+focal+"' not recognised"); return;}
			}
		}	
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
	
	case INIT_POP_PRIOR: 
		{
			auto ty = toLower(get_tag_value("type"));
			ds.prior_type = PriorPos(option_error("type",ty,{"flat","dirichlet"},{ FLAT_PR, DIRICHLET_PR}));
		}
		break;
	
	case COMP_DATA:
		{
			auto name = get_tag_value("class"); if(name == ""){ cannot_find_tag(); return;}
			
			ds.cl = find_cl(p,name);
			if(ds.cl == UNSET){ alert_import("In 'class' the value '"+name+"' is not a classification"); return;}
		}	
		break;
		
	case TRANS_DATA: case SOURCE_DATA: case SINK_DATA:
		{
			switch(cname){
			case TRANS_DATA:
				{
					auto from = get_tag_value("from"); if(from == ""){ cannot_find_tag(); return;}
			
					auto cl_fr = get_cl_from_comp(from,p); 
					if(cl_fr == UNSET){
						alert_import("'"+from+"' is not a valid compartment"); 
						return;
					}
				
					auto to = get_tag_value("to"); if(to == ""){ cannot_find_tag(); return;}
			
					auto cl_to = get_cl_from_comp(to,p);
					if(cl_to == UNSET){ 
						alert_import("'"+to+"' is not a valid compartment"); 
						return;
					}
				
					if(cl_fr != cl_to){ 
						alert_import("'"+from+"' and '"+to+"' cannot be in different classifications");
						return;
					}	
	
					auto c_fr = find_c(p,cl_fr,from);
					if(c_fr == UNSET){
						alert_import("Compartment '"+from+"' not recognised");
						return;
					}
					
					auto c_to = find_c(p,cl_fr,to);
					if(c_to == UNSET){
						alert_import("Compartment '"+to+"' not recognised");
						return;
					}
					
					const auto &claa = sp.cla[cl_fr];
					
					auto tr = 0u; 
					while(tr < claa.tra.size() && !(claa.tra[tr].i == c_fr && claa.tra[tr].f == c_to)) tr++;
					
					if(tr == claa.tra.size()){
						alert_import("Transition not recognised");
						return;
					}
					
					ds.tr = tr;
					ds.cl = cl_fr;
				}
				break;
				
			case SOURCE_DATA:
				{
					auto to = get_tag_value("to"); if(to == ""){ cannot_find_tag(); return;}
			
					auto cl_to = get_cl_from_comp(to,p);
					if(cl_to == UNSET){ 
						alert_import("'"+to+"' is not a valid compartment"); 
						return;
					}
					
					auto c_to = find_c(p,cl_to,to);
					if(c_to == UNSET){
						alert_import("Compartment '"+to+"' not recognised");
						return;
					}
					
					const auto &claa = sp.cla[cl_to];
					
					auto tr = 0u; 
					while(tr < claa.tra.size() && !(claa.tra[tr].i == SOURCE && claa.tra[tr].f == c_to)) tr++;
					
					if(tr == claa.tra.size()){
						alert_import("Transition not recognised");
						return;
					}
					
					ds.tr = tr;
					ds.cl = cl_to;
				}
				break;
				
			case SINK_DATA:
				{
					auto from = get_tag_value("from"); if(from == ""){ cannot_find_tag(); return;}

					auto cl_fr = get_cl_from_comp(from,p);
					if(cl_fr == UNSET){
						alert_import("'"+from+"' is not a valid compartment"); 
						return;
					}
					
					auto c_fr = find_c(p,cl_fr,from);
					if(c_fr == UNSET){
						alert_import("Compartment '"+from+"' not recognised");
						return;
					}
					
					const auto &claa = sp.cla[cl_fr];
					
					auto tr = 0u; 
					while(tr < claa.tra.size() && !(claa.tra[tr].i == c_fr && claa.tra[tr].f == SINK)) tr++;
					
					if(tr == claa.tra.size()){
						alert_import("Transition not recognised");
						return;
					}
					
					ds.tr = tr;
					ds.cl = cl_fr;
				}
				break;
				
			default: alert_import("Should not be default1"); return;
			}
		
			
			
			auto obsran = toLower(get_tag_value("obsrange")); if(obsran == ""){ cannot_find_tag(); return;}
			
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
	
	case TEST_DATA:
		{
			auto Se = get_tag_value("Se"); if(Se == ""){ cannot_find_tag(); return;}
			om.Se= add_equation_info(Se,SE,p);
			
			auto Sp = get_tag_value("Sp"); if(Sp == ""){ cannot_find_tag(); return;}
			om.Sp = add_equation_info(Sp,SP,p);
			
			auto pos = get_tag_value("pos");
			if(pos == "") pos = "1";
			
			auto neg = get_tag_value("neg");
			if(neg == "") pos = "0";
			
			if(pos == neg){ alert_import("'pos' and 'neg' cannot both have the same value"); return;}
				
			om.diag_pos = pos; om.diag_neg = neg;
				
			auto comp = get_tag_value("comp"); if(comp == ""){ cannot_find_tag(); return;}
		
			auto spl = split(comp,',');
			
			auto cl = get_cl_from_comp(spl[0],p);
			if(cl == UNSET){ alert_import("Value '"+spl[0]+"' is not a compartment"); return;}
			
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
	
	case POP_DATA: case SET_TRAPS_DATA:
		{	
			auto &filt = ds.comp_filt;
			
			filt.cla.resize(sp.ncla);
			for(auto cl = 0u; cl < sp.ncla; cl++){
				const auto &claa = sp.cla[cl];
				
				filt.cla[cl].type = ALL_FILT;
				filt.cla[cl].comp.resize(claa.ncomp,false);
			}
			
			auto filter = get_tag_value("filter");
			
			if(filter != ""){
				auto spl = split(filter,',');
				for(auto j = 0u; j < spl.size(); j++){
					auto spl2 = split(spl[j],':');
					if(spl2.size() != 2){ 
						alert_import("In 'filter' error understanding '"+filter+"'"); 
						return;
					}
					
					if(spl2[0] == ""){
						alert_import("In 'filter' the value '"+filter+"' does not specify a classification");
						return;
					}
					
					auto cl = find_cl(p,spl2[0]);
					if(cl == UNSET){ 
						alert_import("In 'filter' the value '"+spl2[0]+"' is not a classification"); 
						return;
					}
					
					if(toLower(spl2[1]) == "file"){
						filt.cla[cl].type = FILE_FILT;
					}
					else{
						filt.cla[cl].type = COMP_FILT;
						
						auto spl3 = split(spl2[1],'|');
						for(auto k = 0u; k < spl3.size(); k++){
							auto c = find_c(p,cl,spl3[k]);
							
							if(spl3[k] == ""){
								alert_import("In 'filter' the classification '"+spl2[0]+"' is not set");
								return;
							}
							
							if(c == UNSET){ 
								alert_import("In 'filter' the value '"+spl3[k]+"' is not a compartment in '"+spl2[0]+"'"); 
								return;
							}
							
							filt.cla[cl].comp[c] = true;
						}
					}
				}
			}
			
			if(cname == SET_TRAPS_DATA){
				auto prob = get_tag_value("prob"); if(prob == ""){ cannot_find_tag(); return;}
				om.trapprob = add_equation_info(prob,TRAP_PROB,p);
			}

			if(cname == POP_DATA){
				load_obs_model(om);
			}
		}
		break;
	
	case POP_TRANS_DATA:
		{
			auto from = get_tag_value("from"); if(from == ""){ cannot_find_tag(); return;}
			auto from_spl = split(from,',');
			
			auto to = get_tag_value("to"); if(to == ""){ cannot_find_tag(); return;}
			auto to_spl = split(to,',');

			if(from_spl.size() != to_spl.size()){
				alert_import("The 'from' and 'to' must have the same size list");
				return;
			}
			
			auto cl_sel = UNSET;
			for(auto i = 0u; i < from_spl.size(); i++){
				auto cl = get_cl_from_comp(from_spl[i],p);
				if(cl == UNSET){
					alert_import("The value '"+from_spl[i]+"' is not a compartment");
					return;
				}
				
				if(cl_sel == UNSET) cl_sel = cl;
				else{
					if(cl != cl_sel){
						alert_import("In 'from' the value '"+from+"' has compartments from different classifications");
						return;
					}
				}
			}
			
			if(cl_sel == UNSET){
				alert_import("No transitions have been selected");
				return;
			}
			ds.cl = cl_sel;
			
			const auto &claa_sel = sp.cla[cl_sel];
			
			ds.trans_filt.resize(claa_sel.ntra,false);
		
			for(auto j = 0u; j < from_spl.size(); j++){
				auto name = from_spl[j]+"→"+to_spl[j];
				
				auto tr = find_tr(p,cl_sel,name);
				if(tr == UNSET){ alert_import("The transition '"+name+"' does not exist"); return;}
				
				ds.trans_filt[tr] = true;
			}
		
			auto &filt = ds.comp_filt;
			
			filt.cla.resize(sp.ncla);
			for(auto cl = 0u; cl < sp.ncla; cl++){
				const auto &claa = sp.cla[cl];
				
				filt.cla[cl].type = ALL_FILT;
				filt.cla[cl].comp.resize(claa.ncomp,false);
			}
		
			auto filter = get_tag_value("filter");
			if(filter != ""){
				auto spl = split(filter,',');
				for(auto j = 0u; j < spl.size(); j++){
					auto spl2 = split(spl[j],':');
					if(spl2.size() != 2){ 
						alert_import("In 'filter' error understanding '"+filter+"'"); 
						return;
					}
					
					if(spl2[0] == ""){
						alert_import("In 'filter' the value '"+filter+"' does not specify a classification");
						return;
					}
					
					auto cl = find_cl(p,spl2[0]);
					if(cl == UNSET){
						alert_import("In 'filter' the value '"+spl2[0]+"' is not a classification"); 
						return;
					}
					
					if(cl == cl_sel){ 
						alert_import("In 'filter' the value '"+spl2[0]+"' cannot be the same as the transition classification"); 
						return;
					}
						
					if(toLower(spl2[1]) == "file"){
						filt.cla[cl].type = FILE_FILT;
					}
					else{
						filt.cla[cl].type = COMP_FILT;
				
						auto spl3 = split(spl2[1],'|');
						for(auto k = 0u; k < spl3.size(); k++){
							if(spl3[k] == ""){
								alert_import("In 'filter' the classification '"+spl2[0]+"' is not set");
								return;
							}
							
							auto c = find_c(p,cl,spl3[k]);
							if(c == UNSET){
								alert_import("In 'filter' the value '"+spl3[k]+"' is not a compartment in '"+spl2[0]+"'");
								return;
							}
							
							filt.cla[cl].comp[c] = true;
						}
					}
				}
			}
		
			load_obs_model(om);
		}
		break;
		
	case GENETIC_DATA:
		{
			auto root = get_tag_value("root"); if(root == ""){ cannot_find_tag(); return;}
			ds.SNP_root = root;
		}
		break;
		
	default: alert_import("Should not be default2"); return;
	}
	
	if(set_loadcol(cname,ds) == false) return;
	
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
	
	//print_table(ds.table);
		
	data_source_check_error(ds);
	
	ds.index = model.species[p].source.size();
	
	model.species[p].source.push_back(ds);
}


/// Ignores sample command 
void Input::dummy_file_command()
{
	auto file = get_tag_value("file"); 
	auto chain = get_tag_value("chain"); 
}


"use strict";
// Function which implement import commands

/// Imports a data-table 
function import_data_table_command(cname)
{
	let siminf = "inf";

	switch(cname){
	case "add-pop-sim": cname = "add-pop"; siminf = "sim"; break;
	case "remove-pop-sim": cname = "remove-pop"; siminf = "sim"; break;	
	case "add-ind-sim": cname = "add-ind"; siminf = "sim"; break;
	case "remove-ind-sim": cname = "remove-ind"; siminf = "sim"; break;
	case "move-ind-sim": cname = "move-ind"; siminf = "sim"; break;
	case "init-pop-sim": cname = "init-pop"; siminf = "sim"; break;
	case "add-pop-post-sim": cname = "add-pop"; siminf = "ppc"; break;
	case "remove-pop-post-sim": cname = "remove-pop"; siminf = "ppc"; break;	
	case "add-ind-post-sim": cname = "add-ind"; siminf = "ppc"; break;
	case "remove-ind-post-sim": cname = "remove-ind"; siminf = "ppc"; break;
	case "move-ind-post-sim": cname = "move-ind"; siminf = "ppc"; break;
	}
	
	let i = find(convert,"command",cname);
	if(i == undefined){
		error("Counld not not import data table:"+cname);
		return;
	}
	
	let type = convert[i].type;
	
	let file = get_tag_value("file"); if(file == "") cannot_find_tag();
	
	let tab = load_table(file.te,true,file.sep,file.name);
	
	if(typeof tab == 'string'){
		alert_import(tab);
	}
	
	let cols = get_tag_value("cols");
	let cols_split = cols.split(",");
	
	let p = imp.p;
	if(p == undefined){
		alert_import("To load the data file the species must be set");
	}
	
	let sp = model.species[p];
	
	let spec;
	let info = {p:p, siminf:siminf};
	
	switch(type){                                    // Modifies specification dependent on data source
	case "Init. Pop.":
		{
			let focal = get_tag_value("focal");
		
			if(focal == "") spec = {radio:{value:"File"},radio2:{value:"All"}, focal:{te:""}};
			else{
				let cl = find_nocase(sp.cla,"name",focal);
				if(cl == undefined) alert_import("The focal classification '"+focal+"' is not recognised");
				else focal = sp.cla[cl].name;
		
				spec = {radio:{value:"File"},radio2:{value:"Focal"}, focal:{te:focal}};
			}
			
			let type = get_tag_value("type"); if(type == "") type = "fixed";
			if(option_error("type",type,["fixed","dist"]) == true) return;
				
			if(type == "fixed") spec.radio_dist = {value:"Fixed"};
			else spec.radio_dist = {value:"Dist"};
		}	
		break;
	
	case "Add Pop.": case "Remove Pop.":
		break
		
	case "Add Ind.":
		break;
		
	case "Remove Ind.":
		break;
	
	case "Move Ind.":
		{
			let name = get_tag_value("class"); if(name == "") cannot_find_tag();
		
			let cl = find_nocase(sp.cla,"name",name);
			if(cl == undefined){
				alert_import("For 'class' the value '"+name+"' is not a classification"); 
			}
			
			spec = {cl_drop:{te:sp.cla[cl].name}};
		}
		break;
		
	case "Ind. Eff.":
		{
			let name = get_tag_value("name"); if(name == "") cannot_find_tag();
		
			spec = {drop:{te:name}};
		}
		break;
		
	case "Ind. Group":
		{
			let name = get_tag_value("name"); if(name == "") cannot_find_tag();
			spec = { gname:name};
		}
		break;
	
	case "Compartment":
		{
			let name = get_tag_value("class");
			if(name == "") cannot_find_tag();
			
			let cl = find_nocase(sp.cla,"name",name);
			if(cl == undefined){
				alert_import("For 'class' the value '"+name+"' is not a classification"); 
			}
		
			spec = {cl_drop:{te:sp.cla[cl].name}};
		}	
		break;
		
	case "Transition": 
		{
			let name = get_tag_value("name"); if(name == ""){ cannot_find_tag(); return;}
			name = name.replace(/->/g,"→")
			
			let cl_sel = get_cl_from_trans(name,p);
		
			if(cl_sel == undefined){
				alert_import("For 'name' the value '"+name+"' is not recognised"); return; 
			}

			let filter = get_tag_value("filter");
			let filt = get_filt_from_str(filter,p);
			if(filt == undefined) return;

			let claa_sel = sp.cla[cl_sel];
		
			set_get_tra_filt_from_str(claa_sel,filt,name);

			spec = {cl_drop:{te:claa_sel.name}, filter:filt};

			let or = get_tag_value("obsrange").toLowerCase();
			if(or == "") cannot_find_tag();
			
			if(option_error("obsrange",or,["all","specify"]) == true) return;
			
			spec.time_start = "";
			spec.time_end = "";
			
			switch(or){
			case "all": 
				spec.time_radio = {value:"All"};
				break;
				
			case "specify":
				spec.time_radio = {value:"Spec"}; 
				
				let start = get_tag_value("start");
				if(start == "") cannot_find_tag();
				is_number(start,"start");
				spec.time_start = start;
				
				let end = get_tag_value("end");
				if(end == "") cannot_find_tag();
				is_number(end,"end");
				spec.time_end = end;
				
				if(Number(start) >= Number(end)){
					alert_import("The start time '"+start+"' must be before the end time '"+end+"'");
				}
				break;
			
			default: 
				alert_import("The value '"+or+"' not recognised"); 
				break;
			}
		}
		break;
		
	case "Diag. Test":
		{
			let Se = get_tag_value("Se");
			if(Se == "") cannot_find_tag();
			
			output_error(is_eqn(Se,"Se",{zero_one_range:true}));
			
			let Sp = get_tag_value("Sp");
			if(Sp == "") cannot_find_tag();
			
			output_error(is_eqn(Sp,"Sp",{zero_one_range:true}));
			
			let pos = get_tag_value("pos");
			if(pos == "") pos = "1";
			
			let neg = get_tag_value("neg");
			if(neg == "") neg = "0";
			
			if(pos == neg){
				alert_import("'pos' and 'neg' cannot both have the same value");
			}
			
			let comp = get_tag_value("comp");
			if(comp == "") cannot_find_tag();
		
			let spl = comp.split(",");
			
			let cl = get_cl_from_comp(spl[0],p);
			if(cl == undefined){
				alert_import("Value '"+spl[0]+"' is not a compartment");
			}
			
			let claa = model.species[p].cla[cl];
			
			let cb = {name:claa.name, value:[]};
			
			for(let c = 0; c < claa.ncomp; c++){
				let name = claa.comp[c].name;
				if(find_in(spl,name) != undefined) cb.value[c] = {comp_name_store:name, check:true};
				else cb.value[c] = {comp_name_store:name, check:false};
			}
	
			let Se_eqn = create_equation(Se,"Se");
			let Sp_eqn = create_equation(Sp,"Sp");
			
			spec = {cl_drop:{te:claa.name}, Se_eqn:Se_eqn, Sp_eqn:Sp_eqn, pos_result:pos, neg_result:neg, check_box:cb};
		}
		break;
		
	case "Population":
		{		
			let filter = get_tag_value("filter");
			
			let filt = get_filt_from_str(filter,p);
			if(filt == undefined) return;
		
			spec = {filter:filt};
				
			load_obs_model(spec);
		}
		break;
		
	case "Pop. Trans.":
		{
			let name = get_tag_value("name"); if(name == ""){ cannot_find_tag(); return;}
			name = name.replace(/->/g,"→")
			
			let cl_sel = get_cl_from_trans(name,p);
		
			if(cl_sel == undefined){
				alert_import("For 'name' the value '"+name+"' is not recognised"); return;
			}

			let filter = get_tag_value("filter");
			let filt = get_filt_from_str(filter,p);
			if(filt == undefined) return;

			let claa_sel = sp.cla[cl_sel];
		
			set_get_tra_filt_from_str(claa_sel,filt,name);

			spec = {cl_drop:{te:claa_sel.name}, filter:filt};

			load_obs_model(spec);
		}
		break;
		
	case "Genetic":
		{
			let type = get_tag_value("type"); if(type == "") cannot_find_tag();
			if(option_error("type",type,["matrix","snp"]) == true) return;
			
			let snp_root = get_tag_value("root");
			if(type == "snp"){
				if(snp_root == "") cannot_find_tag();
			}
			
			let mut_rate = get_tag_value("mut-rate"); if(mut_rate == "") cannot_find_tag();
			let seq_var = get_tag_value("seq-var"); if(seq_var == "") cannot_find_tag();
			spec = {type_radio:{value:type}, snp_root:snp_root, mut_rate_eqn:create_equation(mut_rate,"mut_rate"), seq_var_eqn:create_equation(seq_var,"seq_var"), };
		}
		break;
	}
	
	start_data_source(type,spec,info);
	
	let so = edit_source;
	
	if(type == "Genetic"){
		switch(so.spec.type_radio.value){
		case "matrix": set_genetic_matrix_columns(tab.ele,so); break;
		case "snp": set_SNP_columns(tab,so); break;
		}
	}
	
	// Determines if dates need to be converted to times
	so.time_format = "float";
	for(let c = 0; c < so.load_col.length; c++){
		if(so.load_col[c].heading == "t"){
			let date = get_tag_value("date");
				
			if(date != ""){
				if(option_error("date",date,["dd/mm/yyyy","dd/mm/yy","dd.mm.yyyy","yyyy-mm-dd","mm/dd/yyyy","mm/dd/yy"]) == true) return;
				
				so.time_format = date;
			}
		}
	}
	
	so.data_table_use = "not loaded";
	
	let col_name = [];
	for(let c = 0; c < so.load_col.length; c++){
		col_name.push(so.load_col[c].heading);
	}
	
	if(cols != ""){
		let spl = cols.split(",");

		if(spl.length != col_name.length){
			alert_import("'cols' does not have the correct number of entries (expected something in the order '"+stringify(col_name)+"')"); 
		}
		
		for(let i = 0; i < spl.length; i++){
			for(let j = 0; j < col_name.length; j++){
				if(spl[i] == col_name[j] && i != j){
					alert_import("'cols' does not have the correct order (expected something in the order '"+stringify(col_name)+"')"); 
				}
			}
		}
		
		col_name = spl;
	}

	so.table = get_subtable(tab,col_name);	
	if(so.table.error != "") alert_import("File problem – "+so.table.error);

	convert_date_time(so);

	so.table_loaded = true;
	so.table.edit = true;

	switch(type){                                    // Adds prior for information for total population
	case "Init. Pop.":
		{
			init_pop_convert_to_graphical(so,true);

			if(so.spec.radio_dist.value == "Dist" && so.spec.radio2.value == "All"){
				so.pop_dist = unset_prior("pop_prior");
				
				let prior = get_tag_value("prior").toLowerCase(); 
				let pri = convert_text_to_prior(prior,prior_pos);
				if(pri.err == true){
					alert_import("For 'prior' error with expression '"+prior+"': "+pri.msg);
					return;
				}
				so.pop_dist = pri;
			}
		}
		break;
	} 
	
	data_source_check_error("import",so)
	
	edit_source.info.imp = true;
	data_source("Add",edit_source);
}


/// Adds a species to the model
function species_command(loop)
{
	let name = get_tag_value("name");
	if(name == "") cannot_find_tag();
	
	check_name_input(name,"Species name",true);
	
	let type = get_tag_value("type").toLowerCase();
	if(type == "") cannot_find_tag();
	
	if(option_error("type",type,["individual","population"]) == true) return;

	let trans_tree = false; 
	if(type == "individual"){
		let trans_tree_str = get_tag_value("trans-tree").toLowerCase();
		if(trans_tree_str != ""){
			if(trans_tree_str == "on") trans_tree = true;
			else{
				if(trans_tree_str != "off"){ alert_import("'trans-tree' must be either 'off' or 'on'"); return;}
			}
		}
	}
	
	if(loop == 0){
		let type_conv = type.substr(0,1).toUpperCase()+type.substr(1);
		output_error(model.add_species(name,type_conv,trans_tree));
		imp.p = model.species.length-1;		
	}
	else{
		imp.p = find(model.species,"name",name);
		if(imp.p == undefined) alert_import("Problem finding species '"+name+"'");
	}

	imp.cl = undefined;	
}


/// Adds a classification to the model
function classification_command(loop)
{
	let p = imp.p;

	if(p == undefined){
		alert_import("A species must be defined before a classification can be added",);
	}
		
	let sp = model.species[p];
	
	let name = get_tag_value("name");
	if(name == "") cannot_find_tag();
	
	check_name_input(name,"Classification name",true);
	
	let clone = get_tag_value("clone");
	if(clone != ""){                                 // Clones a classification from another species
		if(loop == 0){
			let p2 = find(model.species,"name",clone);
			if(p2 == undefined){
				alert_import("For 'clone' cannot find the species '"+clone+"'");
			}
			
			let sp2 = model.species[p2];
			
			let cl2 = find(sp2.cla,"name",name);
			
			if(cl2 == undefined){
				alert_import("In clone species '"+clone+"' cannot find the classification '"+name+"'");
			} 
				
			model.clone_class(p,p2,cl2)	
			sp.cla[sp.ncla-1].import_line = imp.line;
		}
		else{
			imp.cl = find(model.species[p].cla,"name",name);
			if(imp.cl == undefined) alert_import("Problem finding classification '"+name+"'");
		}
	}
	else{
		let index = get_tag_value("index");
		if(index == "") cannot_find_tag();
		if(index.length != 1){
			alert_import("The index '"+index+"' must be just a single character");
		}
		
		if(index == "t"){
			alert_import("The index 't' cannot be used because it is reserved for time variation");
		}	
		
		if(find_in(alphabet,index) == undefined){
			alert_import("Index '"+index+"' must be from the lowercase alphabet");
		}
		
		let op = {};
		
		let coord = get_tag_value("coord");
		if(coord != ""){
			coord = coord.toLowerCase();
			if(coord != "latlng" && coord != "cartesian"){
				alert_import("'coord' must be either 'latlng' or 'cartesian'");
			}					
			op.coord = coord;
		}
		
		if(coord == "latlng"){
			let def_map = get_tag_value("map");
			if(def_map != ""){
				if(def_map != "load"){
					alert_import("'map' must be set to 'load'");
				}
				op.default_map = true;
			}
		}
		
		if(loop == 0){
			let warn_name = model.check_name(name,[]);
			if(warn_name) output_error(err(warn_name));
				
			output_error(model.add_classification(p,name,index,op));
		
			sp.cla[sp.ncla-1].import_line = imp.line;
		
			imp.cl = sp.ncla-1;		
		}
		else{
			imp.cl = find(model.species[p].cla,"name",name);
			if(imp.cl == undefined) alert_import("Problem finding classification '"+name+"'"); 
		}
	}
}


/// Sets the classification which is being worked on
function set_command()
{
	let p = imp.p;
	let sp = get_tag_value("species");
	if(sp != ""){
		imp.p = find(model.species,"name",sp);
		if(imp.p == undefined){
			alert_import("Cannot set the species '"+sp+"' as this does not exist");
		}
	}
	
	let name = get_tag_value("classification");
	
	if(name != ""){
		let p = imp.p;
		if(p == undefined){
			alert_import("A species must be set before a 'classification' can be set");
		};
	
		imp.cl = find(model.species[p].cla,"name",name);
		if(imp.cl == undefined){
			alert_import("Cannot set the classification '"+name+"'");
		}
	}
}


/// Sets a camera in a classification
function camera_command()
{
	let claa = get_claa(); if(claa == undefined) return;

	let grid = get_tag_value("grid").toLowerCase(); 
	if(grid == "") grid = "off";
	if(option_error("grid",grid,["off","on"]) == true) return;
	
	claa.camera.grid = grid;
	
	let sc = get_tag_value("comp-scale"); 
	if(sc == "") sc = "1";
	is_positive(sc,"comp-scale");
	sc = Number(sc);
	if(sc < 0.1 || sc > 10) alert_import("For 'comp-scale' the value '"+sc+"' must be between 0.1 and 10"); 
	claa.camera.slider.value = Math.log(sc);
		
	if(option_error("grid",grid,["off","on"]) == true) return;
	claa.camera.grid = grid;
	
	let scale = get_tag_value("scale"); if(scale == "") cannot_find_tag();
	is_positive(scale,"scale");
	
	claa.camera.set = true;
	
	claa.camera.scale = scale;

	switch(claa.camera.coord){
	case "cartesian":
		let x = get_tag_value("x"); if(x == "") cannot_find_tag();
		let y = get_tag_value("y"); if(y == "") cannot_find_tag();
		
		is_number(x,"x");
		x = Number(x)*import_scale_factor; 
		
		is_number(y,"y");
		y = Number(y)*import_scale_factor; 
		
		claa.camera.x = x;
		claa.camera.y = y;
		break;
		
	case "latlng":
		let lat = get_tag_value("lat"); if(lat == "") cannot_find_tag();
		let lng = get_tag_value("lng"); if(lng == "") cannot_find_tag();
		check_latlng(lat,"lat",lng,"lng");
		
		let p = transform_latlng(lng,lat);
		claa.camera.x = p.x;
		claa.camera.y = p.y;
		break;
		
	default:
		alert_import("'coord' must be 'cartesian' or 'latlng'"); 
		break;
	}
}


/// Sets a camera in a classification
function warning_command()
{
	let warn = get_tag_value("text"); if(warn == "") cannot_find_tag();
	
	model.warning = {show:false, te:warn};
}


/// Gets a list of tags from a file 
function get_tags_list(file)
{
	let tags_list=[];
	let tab = load_table(file.te,true,file.sep,file.name);
	let N = tab.ncol;
	
	let head = tab.heading;
	let tags_root=[];
	for(let i = 0; i < imp.tags.length; i++){
		let ta = copy(imp.tags[i]);
		if(ta.name == "file") ta.value = "";
		tags_root.push(ta);
	}

	for(let r = 0; r < tab.nrow; r++){
		let row = tab.ele[r];
		let tags = copy(tags_root);
		for(let i = 0; i < N; i++){
			let te = row[i];
			if(te != "") tags.push({name:head[i],value:te,done:0});
		}
		
		tags_list.push(tags);		
	}
	
	return tags_list;
}
	

/// Checks all tags used 
function check_tags_used(r,tags)
{
	for(let i = 0; i < tags.length; i++){
		let tag = tags[i];
		if(tag.done == 0){
			alert_import("The value '"+tag.value+"' for property '"+tag.name+"' is not used (line "+(r+2)+" in file).");
			return;
		}			
	}
}


/// Adds a compartment to the model
function compartment_command()
{
	compartment_command2(imp.tags)
}


/// Adds compartments to the model
function compartment_all_command()
{
	let file = get_tag_value("file"); if(file == "") cannot_find_tag();

	let tags_list = get_tags_list(file);
	
	for(let i = 0; i < tags_list.length; i++){
		let tags = tags_list[i];
		imp.all_row = "(line "+(i+2)+" in file)";
		compartment_command2(tags);
		check_tags_used(i,tags);
	}
	
	for(let m = 0; m < imp.tags.length; m++) imp.tags[m].done = 1; 
}

	
/// Adds a compartment to the model
function compartment_command2(tags)
{
	let claa = get_claa(); if(claa == undefined) return;
	let p = imp.p;
	let cl = imp.cl;
	
	let name = get_tag_val("name",tags); if(name == "") cannot_find_tag();
	let warn_comp = model.check_comp_name(name,"full");
	if(warn_comp != undefined) alert_import(warn_comp);
	
	let cam = claa.camera;
	
	let color = get_tag_val("color",tags); 
	if(color == ""){ 
		let check_col;
		switch(name){
		case "S": check_col = DGREEN; break;
		case "E": check_col = ORANGE; break;
		case "I": check_col = RED; break;
		case "R": check_col = BLUE; break;
		case "D": check_col = BLACK; break;
		}
			
		if(check_col != undefined){
			if(find(claa.comp,"col",check_col) == undefined) color = check_col;
		}
		
		if(color == ""){
			let j = 0; 
			while(j < auto_color.length && find(claa.comp,"col",auto_color[j]) != undefined) j++;
			
			if(j == auto_color.length) j = Math.floor(Math.random()*auto_color.length);
				
			color = auto_color[j];
		}	
	}
	else{
		if(is_Color(color) == false) alert_import("'"+color+"' is not a valid color."); 
	}
	
	let fix = false;
	let fix_str = get_tag_val("fix",tags).toLowerCase(); 
	if(fix_str != ""){
		if(fix_str == "true") fix = true;
		else{
			if(fix_str == "false") fix = false;
			else{
				alert_import("'fix' must be 'true' or 'false'"); 
			}
		}
	}
	
	let infected = false;
	let infected_str = get_tag_val("infected",tags).toLowerCase(); 
	if(infected_str != ""){
		if(infected_str == "true") infected = true;
		else{
			if(infected_str == "false") infected = false;
			else{
				alert_import("'infected' must be 'true' or 'false'"); 
			}
		}
	
		let sp = model.species[p];
		let claname = sp.cla[cl].name;
		
		if(sp.infection_cl.te == select_str) sp.infection_cl.te = claname;
		else{
			if(sp.infection_cl.te != claname){
				alert_import("Cannot have 'infected' specified in more than one classification for species '"+sp.name+"'.");
			}
		}
	}

	switch(cam.coord){
	case "cartesian":
		{
			let x = get_tag_val("x",tags);
			let y = get_tag_val("y",tags);
			if(x == "" && y == ""){
				x = "auto"; y = "auto";
			}
			else{
				x = get_tag_val("x",tags); if(x == ""){ cannot_find_tag(true); return;}
				y = get_tag_val("y",tags); if(y == ""){ cannot_find_tag(true); return;}
				
				is_number(x,"x");
				x = Number(x)*import_scale_factor; 
				
				is_number(y,"y");
				y = Number(y)*import_scale_factor; 
			
				for(let i = 0; i < claa.ncomp; i++){
					if(claa.comp[i].x == x && claa.comp[i].y == y){		
						alert_import("Compartments '"+name+"' and '"+claa.comp[i].name+"' cannot have the same position");
					}	
				}
			}
			
			output_error(model.add_compartment(name,p,cl,x,y,color,fix,infected,false));
			claa.comp[claa.ncomp-1].import_line = imp.line;
		}
		break;
	
	case "latlng":
		{
			let bound_file = get_tag_val("boundary",tags);
			if(bound_file != ""){
				let claa = get_claa(); if(claa == undefined) return;

				output_error(add_individual_compartment_boundary(name,imp.p,imp.cl,bound_file,color,infected));
			}
			else{
				let lat = get_tag_val("lat",tags); if(lat == "") cannot_find_tag();
				let lng	= get_tag_val("lng",tags); if(lng == "") cannot_find_tag();
				
				check_latlng(lat,"lat",lng,"lng");
				
				let pt = transform_latlng(lng,lat);
				
				output_error(model.add_latlng_compartment(name,imp.p,imp.cl,pt.x,pt.y,color,fix,infected));		
			}
		}
		break;

	default: error("Option not recognised 50"); break;
	}
	
	let markov_branch = false;
	let markov_branch_str = get_tag_val("branch-prob",tags); 
	if(markov_branch_str != ""){
		if(markov_branch_str == "true") markov_branch = true;
		else{
			if(markov_branch_str == "false") markov_branch = false;
			else{
				alert_import("'branch-prob' must be 'true' or 'false'"); 
			}
		}
	}
	let co = claa.comp[claa.ncomp-1];

	co.markov_branch = markov_branch;
}


/// Sets the data directory
function datadir_command()
{
	let dir = get_tag_value("folder"); 
	if(dir == "") cannot_find_tag();

	check_char_allowed(dir,"<>\"|?*");
}


/// Determines if a file exists or not
function check_file_exists(file)
{
}


/// Checks for disallowed characters
function check_char_allowed(st,not_allowed)
{
	for(let i = 0; i < st.length; i++){
		let ch = st.substr(i,1);
		if(find_in(not_allowed,ch) != undefined){
			alert_import("In '"+st+"' the character '"+ch+"' is not allowed");
		}
	}
}


/// Sets a description of the analysis
function description_command()
{
	let te = get_tag_value("text"); if(te == "") cannot_find_tag();
	if(typeof te == "object") te = te.te;
	else te = te.replace(/\|/g,"\n");

	model.description.te = te;
}
	
/// Adds a transition / source / sink to the model
function transition_command()
{
	transition_command2(imp.tags);
}


/// Adds compartments to the model
function transition_all_command()
{
	let file = get_tag_value("file"); if(file == "") cannot_find_tag();

	let tags_list = get_tags_list(file);

	for(let i = 0; i < tags_list.length; i++){
		let tags = tags_list[i];
		transition_command2(tags);
		check_tags_used(i,tags);
	}
	
	for(let m = 0; m < imp.tags.length; m++) imp.tags[m].done = 1; 
}


/// Adds a transition / source / sink to the model
function transition_command2(tags)
{
	let claa = get_claa(); if(claa == undefined) return;
	let cam = claa.camera;

	let p = imp.p;
	let cl = imp.cl;
	
	let value = get_tag_val("value",tags); if(value == "") cannot_find_tag();
	
	let te = get_tag_val("name",tags); if(te == "") cannot_find_tag();
	
	let spl = te.split("->");
	if(spl.length != 2){ alert_import("Transition name '"+te+"' is not understood. Expected format: 'Initial compartment'->'final compartment'."); return;}

	let fr = spl[0].trim();
	if(fr == "+") fr = SOURCE; 
	
	let to = spl[1].trim();
	if(to == "-") to = SINK; 
	
	let ci, cf;
	
	if(fr == SOURCE) ci = SOURCE;
	else{
		ci = hash_find(claa.hash_comp,fr);
		if(ci == undefined) cannot_find_trans_comp(te,p,cl,fr);
	}
	
	if(to == SINK) cf = SINK;
	else{
		cf = hash_find(claa.hash_comp,to);
		if(cf == undefined) cannot_find_trans_comp(te,p,cl,to);
	}
	
	let pos = dist_pos;	if(ci == SOURCE) pos = source_dist_pos;

	for(let tr = 0; tr < claa.ntra; tr++){
		if(claa.tra[tr].i == ci && claa.tra[tr].f == cf){
			if(ci == SOURCE) alert_import("Source to '"+to+"' already exists");
			else{
				if(cf == SINK) alert_import("Sink from '"+fr+"' already exists");
				else alert_import("Transition from '"+fr+"' to '"+to+"' already exists");
			}
		}
	}

	let midp=[];

	switch(cam.coord){
	case "cartesian": 
		{
			if(ci == SOURCE){ 
				let x = get_tag_val("x",tags);
				let y = get_tag_val("y",tags);
				if(x == "" && y == ""){
					x = "auto"; y = "auto";
				}
				else{
					if(x == "" || y == ""){
						alert_import("The 'x' and 'y' positions must be set for the source");
					}

					is_number(x,"x");
					x = Number(x)*import_scale_factor; 
					
					is_number(y,"y");
					y = Number(y)*import_scale_factor; 
				}
				
				midp.push({x:x, y:y});
			}

			let xmid = get_tag_val("mid-x",tags);
			let ymid = get_tag_val("mid-y",tags);
			if(xmid != ""){
				var xlist = xmid.split(",");
				if(ymid == "") alert_import("Must set 'mid-x' and 'mid-y' points"); 
				var ylist = ymid.split(",");
				if(xlist.length != ylist.length){
					alert_import("Must have the same number of 'x-mid' and 'y-mid' points");
				}

				for(let n = 0; n < xlist.length; n++){
					is_number(xlist[n],"mid-x");
					xlist[n] = Number(xlist[n])*import_scale_factor; 
					
					is_number(ylist[n],"mid-y");
					ylist[n] = Number(ylist[n])*import_scale_factor; 
				}
				
				for(let n = 0; n < xlist.length; n++){
					midp.push({x:xlist[n], y:ylist[n]});
				}
			}
			else{
				if(ymid != "") alert_import("Must set 'mid-x' and 'mid-y' points");
			}

			if(cf == SINK){ 
				let x = get_tag_val("x",tags);
				let y = get_tag_val("y",tags);
				
				if(x == "" && y == ""){
					x = "auto"; y = "auto";
				}
				else{
					if(x == "" || y == ""){
						alert_import("The 'x' and 'y' positions must be set for the sink");
					}

					is_number(x,"x");
					x = Number(x)*import_scale_factor; 
				
					is_number(y,"y");
					y = Number(y)*import_scale_factor; 
				}
				
				midp.push({x:x, y:y});
			}
		}
		break;

	case "latlng":
		{
			if(ci == SOURCE){
				let lng = get_tag_val("lng",tags);
				let lat = get_tag_val("lat",tags);
				check_latlng(lat,"lat",lng,"lng");
				
				let p = transform_latlng(lng,lat);
				midp.push(p);
			}

			let lngmid = get_tag_val("mid-lng",tags);
			let latmid = get_tag_val("mid-lat",tags);
			if(lngmid != ""){
				var lnglist = lngmid.split(",");
				if(latmid == ""){ 
					alert_import("Must set 'mid-lat' and 'mid-lng' points");
				}

				var latlist = latmid.split(",");
				if(lnglist.length != latlist.length){ 
					alert_import("Must have the same number of 'lat-mid' and 'lng-mid' points");
				}

				for(let n = 0; n < lnglist.length; n++){
					check_latlng(latlist[n],"lat-mid",lnglist[n],"lng-mid");
					
					let p = transform_latlng(lnglist[n],latlist[n]);
					midp.push(p);
				}
			}
			else{
				if(latmid != "") alert_import("Must set 'mid-lat' and 'mid-lng' points");
			}

			if(cf == SINK){ 
				let lng = get_tag_val("lng",tags);
				let lat = get_tag_val("lat",tags);
				check_latlng(lat,"lat",lng,"lng");
				 
				let p = transform_latlng(lng,lat);
				midp.push(p);
			}
		}
		break;

	default: error("Option not recognised 130"); break;
	}

	let trans_def = extract_trans_def(value);
	if(trans_def == undefined){
		alert_import("In the expression '"+value+"' the transition distribution is incorrectly specified. Expected format: 'exp(rate:...)', 'gamma(mean:...,cv:...)', 'erlang(mean:...,shape:)', 'log-normal(mean:...,cv:...)', 'weibull(scale:...,shape:...)' or 'period(time:...)'.");
	}
	
	option_error("type",trans_def.type,pos);

	if(output_error(model.add_transition(p,cl,ci,cf,midp,trans_def.type,"simple"))) return;

	let traa = claa.tra[claa.ntra-1];
	traa.import_line = imp.line;
	
	traa.bp_set = false;
	
	let val = traa.value;
	
	switch(trans_def.type){
	case "exp(rate)":
		import_eqn_value(val.rate_eqn,trans_def.rate);
		break;
		
	case "exp(mean)":
		import_eqn_value(val.mean_eqn,trans_def.mean);
		break;
		
	case "gamma":
		import_eqn_value(val.mean_eqn,trans_def.mean);
		import_eqn_value(val.cv_eqn,trans_def.cv);
		break;
		
	case "erlang":
		{
			import_eqn_value(val.mean_eqn,trans_def.mean);
			
			let shape = trans_def.shape;
			let num = Number(shape);
			if(isNaN(shape) || num <= 0 || Math.round(num) != num){
				alert_import("For an Erlang distribution the shape parameter '"+shape+"' must be a positive integer");				
			}
			
			import_eqn_value(val.shape_erlang,shape);
		}
		break;
		
	case "log-normal":
		import_eqn_value(val.mean_eqn,trans_def.mean);
		import_eqn_value(val.cv_eqn,trans_def.cv);
		break;
		
	case "weibull":
		import_eqn_value(val.scale_eqn,trans_def.scale);
		import_eqn_value(val.shape_eqn,trans_def.shape);
		break;

	case "period":
		{
			let period = Number(trans_def.time);
			if(isNaN(trans_def.time) || period <= 0){
				alert_import("The period must be a positive number");				
			}
		
			import_eqn_value(val.mean_eqn,trans_def.time);
		}
		break;
		
	default: error("Option not recognised 51"); break;
	}
	
	let bp = get_tag_val("bp",tags);
	if(bp != ""){
		traa.bp_set = true;
		import_eqn_value(val.bp_eqn,bp);
		
		if(val.bp_eqn.te == "×"){
			traa.branch_select = false;
			val.bp_eqn.te = "";
		}
		else traa.branch_select = true;
	}
}


/// Adds a label to the model
function label_command()
{
	let claa = get_claa(); if(claa == undefined) return;
	
	let te = get_tag_value("text"); if(te == "") cannot_find_tag();
	let size = get_tag_value("text-size");
	if(size == "") size = size_annotation_default;
	else{
		is_number(size,"text-size");
		size = Number(size); 
	}
	
	let x, y;
	switch(claa.camera.coord){
	case "cartesian":
		{
			x = get_tag_value("x"); if(x == "") cannot_find_tag();
			y = get_tag_value("y"); if(y == "") cannot_find_tag();
			
			is_number(x,"x");
			x = Number(x)*import_scale_factor; 
			
			is_number(y,"y");
			y = Number(y)*import_scale_factor;
		}
		break;
		
	case "latlng":
		{
			let lat = get_tag_value("lat"); if(lat == "") cannot_find_tag();
			let lng = get_tag_value("lng"); if(lng == "") cannot_find_tag();
			check_latlng(lat,"lat",lng,"lng");
		
			let p = transform_latlng(lng,lat);
			x = p.x;
			y = p.y;
		}
		break;
	
	default:
		alert_import("'coord' must be 'cartesian' or 'latlng'"); 
		break;
	}
	
	let color = get_tag_value("color"); if(color == "") color = BLACK;

	claa.annotation.push({type:"text", te:te, size:size, x:x, y:y, color:color});
}


/// Adds a box around compartments in the model
function box_command()
{
	let claa = get_claa(); if(claa == undefined) return;
	
	let te = get_tag_value("text"); if(te == "") cannot_find_tag();
	let size = get_tag_value("text-size"); if(size == "") size = size_annotation_default;
	else{
		is_number(size,"text-size");
		size = Number(size); 
	}
	
	let comps = get_tag_value("comps"); if(comps == "") cannot_find_tag();
	let color = get_tag_value("color"); if(color == "") color = BLACK;
			
	let list=[];
	let spl = comps.split(",");
	if(spl.length == 0) alert_import("Error 'box' contains no compartments.");
	for(let k = 0; k < spl.length; k++){
		let c = 0; while(c < claa.ncomp && claa.comp[c].name != spl[k]) c++;
		if(c == claa.ncomp){ 
			alert_import("Cannot find compartment '"+spl[k]+"'.");
		}
		list.push(spl[k]);
	}
	
	claa.annotation.push({type:"box", te:te, size:size, comps:list, color:color});
}


/// Adds an annotation map
function map_command()
{
	get_claa();
	let p = imp.p;
	let cl = imp.cl;
	
	let file = get_tag_value("file"); if(filev == "") cannot_find_tag();
	
	let res = load_annotation_map(file.te);
	
	let name = "file"+Math.random();
	map_store.push({name:name, feature:res.feature, box:res.box});
				
	let claa = model.species[p].cla[cl];
	claa.annotation.push({type:"map", map_ref:name});
}

	
/// Sets param_mult
function param_mult_command(per_start,per_end)
{
	let full_name = get_tag_value("name"); if(full_name == "") cannot_find_tag();
	
	let par = get_param_prop(full_name);
	par.full_name = param_name(par);
	
	let f_par = get_param_mult(par);
	model.param_factor.push({f_param:f_par, param:{name:par.name,dep:copy(par.dep), full_name:par.full_name}});

	param_command2(remove_eq_quote(f_par.full_name),per_start,per_end,"mult");
}


/// Sets the value for a parameter in the model
function param_command(per_start,per_end)
{
	let full_name = get_tag_value("name"); if(full_name == "") cannot_find_tag();
	param_command2(full_name,per_start,per_end);
}


/// Processes param / param-mult command
function param_command2(full_name,per_start,per_end,op)
{
	let pp = get_param_prop(full_name);

	let par = create_new_param(pp,"normal");

	par.import_line = imp.line;

	if(par.time_dep == true){
		let knot_times = get_tag_value("knot-times"); if(knot_times == "") cannot_find_tag();
		let warn = check_knot_times(knot_times);
		if(warn != "") alert_import(warn);
		
		change_spline_knot(par,knot_times.split(","));
		par_set_default(par);
		
		let smooth = get_tag_value("smooth").toLowerCase().trim();
		
		par.spline.smooth = default_smooth();
		
		if(smooth != ""){
			par.spline.smooth.check = true;
			
			let spl = smooth.split("(");
			
			if(option_error("'smooth'",spl[0],["normal","log-normal"]) == true) return;
		
			switch(spl[0]){
			case "normal": par.spline.smooth.type = {value:"Normal"}; break;
			case "log-normal": par.spline.smooth.type = {value:"Log-Normal"}; break;
			}
	
			if(spl.length != 2) alert_import("There is syntax error in 'smooth'");
			
			if(spl[1].substr(spl[1].length-1,1) != ")"){
				alert_import("There is syntax error in 'smooth'");
			}
			
			let val = spl[1].substr(0,spl[1].length-1);
			if(isNaN(val)) alert_import("For 'smooth' the value '"+val+"' is not a number");
			if(Number(val) <= 0) alert_import("For 'smooth' the value '"+val+"' is not positive");
			par.spline.smooth.value = val;
		}
	}

	if(par.dep.length > 0){
		par.list = par_find_list(par);
	}
	
	par.factor = false;
	let fact = get_tag_value("factor");
	if(fact != ""){
		if(option_error("factor",fact,["true","false"]) == true) return;
		if(fact == "true"){
			par.factor = true;
			par.pri_pos = prior_factor_pos;
		}
	}
	
	par.factor_weight_on = {check:false};
	let fw = get_tag_value("factor-weight");
	if(fw != ""){
		if(par.factor == false){
			alert_import("'factor-weight' can only be set if 'factor' is set to 'true'"); 
		}
	
		par.factor_weight_on = {check:true};
		
		set_default_factor_weight(par);
		
		let tab = load_table(fw.te,true,fw.sep,fw.name);
			
		if(tab == undefined) return;
		if(typeof tab == 'string') alert_import(tab);
	
		let col_name = copy(par.dep);
		col_name.push("Value");
		
		let subtab = get_subtable(tab,col_name);
		if(subtab.error != "") alert_import(subtab.error);

		let ncol = subtab.ncol;
	
		for(let r = 0; r < subtab.nrow; r++){
			let ind = [];
			for(let i = 0; i < ncol-1; i++){
				ind[i] = find_in(par.list[i],subtab.ele[r][i]);
				if(ind[i] == undefined){
					alert_import("For 'factor-weight' the element '"+subtab.ele[r][i]+"' is not valid (column '"+subtab.heading[i]+"', row "+(r+2)+")");
				}
			}
			let ele = subtab.ele[r][ncol-1];
			let val = Number(ele);
			
			if(isNaN(val)){
				alert_import("For 'factor-weight' the element '"+ele+"' is not a number (column '"+subtab.heading[ncol-1]+"', row "+(r+2)+")");
			}
				
			set_element(par.factor_weight,ind,val);
		}	
	}
	
	let cons = get_tag_value("constant"); 
	let dist = get_tag_value("dist");	
	let dist_split = get_tag_value("dist-split");
	let value = get_tag_value("value");
	let reparam = get_tag_value("reparam"); 
	let prior = get_tag_value("prior"); 
	let prior_split = get_tag_value("prior-split"); 
	
	//let dist_mat = get_tag_value("distance-matrix").toLowerCase(); 
	//let dist_mat = "";

	let param_tag = [];
	param_tag.push({val:cons, tag:"constant"});
	param_tag.push({val:dist, tag:"dist"});
	param_tag.push({val:dist_split, tag:"dist-split"});
	param_tag.push({val:value, tag:"value"});
	param_tag.push({val:reparam, tag:"reparam"});
	param_tag.push({val:prior, tag:"prior"});
	param_tag.push({val:prior_split, tag:"prior-split"});
	//param_tag.push({val:dist_mat, tag:"distance-matrix"});

	for(let j = 0; j < param_tag.length; j++){
		for(let i = j+1; i < param_tag.length; i++){
			if(param_tag[j].val != "" && param_tag[i].val != ""){
				if(!(param_tag[j].tag == "value" && param_tag[j].tag != "prior") && 
					 !(param_tag[j].tag == "value" && param_tag[j].tag != "prior_split")){
					alert_import("'"+param_tag[j].tag+"' and '"+param_tag[i].tag+"' cannot both be set"); 
				}
			}
		}
	}
	
	let j = 0; while(j < param_tag.length && param_tag[j].val == "") j++;
	
	if(par.dist_mat){
		alert_import("The distance matrix '"+dist_matrix_name+"' cannot be set"); 
	}
	
	par.reparam_eqn = "";
	par.reparam_eqn_on = false;

	if(par.dep.length == 0){
		if(value != ""){
			par.value = value;
			par.variety = "normal";
		}
		
		if(reparam != ""){
			par.value = reparam;
			par.variety = "reparam";
		}
		
		if(cons != ""){
			par.value = cons;
			par.variety = "const";
		}
	}
	else{
		if(value != "" || cons != "" || reparam != ""){
			set_default_value(par);
			
			let desc = "For 'value'";
			let valu = value; 
			if(valu == ""){
				if(cons != ""){
					valu = cons; desc = "For 'constant'";
					par.variety = "const";		
				}
				else{
					if(reparam != ""){
						valu = reparam; desc = "For 'reparam'";
						par.variety = "reparam";		
					}
					else error("Problem importing");
				}
			}
		
			if(reparam != "" && is_file(valu) == false){
				par.reparam_eqn = reparam;
				par.reparam_eqn_on = true;
			}
			else{		
				par.value = param_blank(par);
			
				let dim = get_dimensions(par.value);
				let ele_list = get_element_list(par.value,dim);
		
				if(is_file(valu) == false){ // Sets all elements to the same	
					for(let k = 0; k < ele_list.length; k++){
						set_element(par.value,ele_list[k],valu);
					}
				}
				else{
					// Sets a default value of zero
					for(let k = 0; k < ele_list.length; k++){
						set_element(par.value,ele_list[k],"0");
					}
					
					let tab = load_table(valu.te,true,valu.sep,valu.name);
				
					if(tab == undefined) return;
					if(typeof tab == 'string') alert_import(tab);
				
					let col_name = copy(par.dep);
					col_name.push("Value");
					
					let subtab = get_subtable(tab,col_name);
					if(subtab.error != "") alert_import(subtab.error);
			
					let ncol = subtab.ncol;
				
					for(let r = 0; r < subtab.nrow; r++){
						let ind = [];
						for(let i = 0; i < ncol-1; i++){
							ind[i] = find_in(par.list[i],subtab.ele[r][i]);
							if(ind[i] == undefined){
								alert_import(desc+" the element '"+subtab.ele[r][i]+"' is not valid (column '"+subtab.heading[i]+"', row "+(r+2)+")");
							}
						}
						let ele = subtab.ele[r][ncol-1];
						let val;
						
						switch(par.variety){
						case "normal": case "const":
							if(par.factor && ele == "*"){
								val = ele;
							}
							else{
								val = Number(ele);
						
								if(isNaN(val)){
									alert_import(desc+" the element '"+ele+"' is not a number (column '"+subtab.heading[ncol-1]+"', row "+(r+2)+")");
								}
							}
							break;
						
						case "reparam":
							if(!isNaN(ele)) val = Number(ele);
							else{
								let conv_res = detect_greek(ele,0);
								ele = conv_res.te;
								subtab.ele[r][ncol-1] = ele;

								val = ele;
								
								let res = is_eqn(val,"Table element",{});
								if(res.err == true){
									alert_import(desc+" the element '"+ele+"' is not a valid equation (column '"+subtab.heading[ncol-1]+"', row "+(r+2)+")");
								}
							}
							break;
							
						default: error("option not recognised1"); break;
						}
								
						set_element(par.value,ind,val);
					}
				}
				
				let err = check_param_value("Set Param",par,par.value);
				if(typeof err == 'string') alert_import(desc+": "+err);
			}
			par.set = true;
		}
	}

	if(prior != ""){
		let pri = convert_text_to_prior(prior,par.pri_pos);
		if(pri.err == true){
			alert_import("For 'prior' error with expression '"+prior+"': "+pri.msg);
			return;
		}
		par.prior = pri;
	}
	
	if(prior_split != ""){
		par.prior_split_check = {check:true};
		par.prior_split_set = true;
		par.prior_split = param_blank(par);
		if(par.dep.length == 0){
			alert_import("'prior-split' can only be used if the parameter has a dependency.");  
		}
	
		let tab = load_table(prior_split.te,true,prior_split.sep,prior_split.name);
		if(tab == undefined) return;	
		if(typeof tab == 'string') alert_import(tab);
		
		let col_name = copy(par.dep);
		col_name.push("Prior");
		
		let subtab = get_subtable(tab,col_name);
		if(subtab.error != "") alert_import(subtab.error);
	
		let ncol = subtab.ncol;
		
		for(let r = 0; r < subtab.nrow; r++){
			let ind = [];
			for(let i = 0; i < ncol-1; i++){
				ind[i] = find_in(par.list[i],subtab.ele[r][i]);
				if(ind[i] == undefined){ 
					alert_import("The table element '"+subtab.ele[r][i]+"' is not valid (column '"+subtab.heading[i]+"', row "+(r+2)+")");
				}
			}
			
			let pri = convert_text_to_prior(subtab.ele[r][ncol-1],par.pri_pos);
			if(pri.err == true){
				alert_import("The table element '"+subtab.ele[r][ncol-1]+"' is not a valid prior specification (col '"+subtab.heading[ncol-1]+"', row "+(r+2)+"). "+pri);
			}
			
			set_element(par.prior_split,ind,pri);
		}
	}
	
	if(dist != ""){
		let pri = convert_text_to_prior(dist,par.pri_pos,true);
		if(pri.err == true){
			alert_import("For 'dist' error with expression '"+dist+"': "+pri.msg);
			return;
		}
		par.prior = pri;
		par.variety = "dist";
	}
	
			
	if(dist_split != ""){
		par.prior_split_check = {check:true};
		par.variety = "dist";
		par.prior_split_set = true;
		
		par.prior_split_check = {check:true};
		if(par.dep.length == 0){
			alert_import("'dist-split' can only be used if the parameter has a dependency.");  
		}
		
		let tab = load_table(dist_split.te,true,dist_split.sep,dist_split.name);
		//let tab = load_table_from_file(dist_split);
		
		if(tab == undefined) return;	
		if(typeof tab == 'string') alert_import(tab);
		
		let col_name = copy(par.dep);
		col_name.push("Dist");
		
		let subtab = get_subtable(tab,col_name);
		if(subtab.error != "") alert_import(subtab.error);
	
		let ncol = subtab.ncol;
		
		for(let r = 0; r < subtab.nrow; r++){
			let ind = [];
			for(let i = 0; i < ncol-1; i++){
				ind[i] = find_in(par.list[i],subtab.ele[r][i]);
				if(ind[i] == undefined){ 
					alert_import("The table element '"+subtab.ele[r][i]+"' is not valid (column '"+subtab.heading[i]+"', row "+(r+2)+")");
				}
			}
			
			let pri = convert_text_to_prior(subtab.ele[r][ncol-1],par.pri_pos,true);
			if(pri.err == true){
				alert_import("The table element '"+subtab.ele[r][ncol-1]+"' is not a valid distribution specification (col '"+subtab.heading[ncol-1]+"', row "+(r+2)+"). "+pri);
			}
			
			set_element(par.prior_split,ind,pri);
		}
	}
	
	let sim_sample = get_tag_value("sim-sample").toLowerCase(); 
	if(sim_sample == "") sim_sample = "true";
	else{
		if(dist == ""){
			alert_import("'sim-sample' can only be set if a distribution is set through 'dist'");
		}
	}
	if(option_error("sim-sample",sim_sample,["true","false"]) == true) return;
	
	if(dist != ""){
		switch(sim_sample){
		case "false": par.sim_sample = {check:false}; break;
		case "true": par.sim_sample = {check:true}; break;
		}
	}
	
	par.full_name = param_name(par);

	get_reparam_param_list(par);
	
	get_prior_param_list(par);
	
	model.param.push(par);
}


/// Definition for derived quantities
function derived_command()
{
	let full_name = get_tag_value("name"); if(full_name == "") cannot_find_tag();
	let eqn_name = get_tag_value("eqn"); if(eqn_name == "") cannot_find_tag();
	
	let eqn1 = create_equation(full_name,"derive_param");
	if(eqn1.warn.length > 0) alert_import("For 'name': "+eqn1.warn[0].te);
	if(eqn1.param.length != 1){
		alert_import("'name' must contain a single paramter");
	}
	else{
		let name = eqn1.param[0].name;
		if(name == "D") alert_import("Name 'D' is reserved for the distance matrix");
		if(name == "t") alert_import("Name 't' is reserved for time");
	}

	let eqn2 = create_equation(eqn_name,"derive_eqn");
	if(eqn2.warn.length > 0) alert_import("In 'eqn': "+eqn2.warn[0].te);
	
	model.derive.push({eqn1:eqn1, eqn2:eqn2});
}


/// Applies the simulation command
function simulation_command()
{
	let details = model.sim_details;
	
	let start = get_tag_value("start");
	if(start != "") is_number(start,"start");
	details.t_start = start;
	
	let end = get_tag_value("end"); 
	if(end != "") is_number(end,"end");
	details.t_end = end;
	
	if(details.t_start != "" && details.t_end != ""){
		if(Number(details.t_start) >= Number(details.t_end)){
			alert_import("The start time '"+start+"' must be before the end time '"+end+"'");
		}
	}
	
	let num = get_tag_value("number").toLowerCase();
	if(num == "") num = SIM_NUM_DEFAULT;
	is_pos_int(num,"number")
	details.number = Number(num);
	
	details.seed_on = {value:"No"};
	details.seed = SEED_DEFAULT;
	let seed = get_tag_value("seed"); 
	if(seed != ""){ details.seed = seed; details.seed_on.value = "Yes";}
	
	let alg = get_tag_value("algorithm").toLowerCase();
	if(alg == "") alg = "gillespie";
	
	if(option_error("algorithm",alg,sim_alg_list) == true) return;
	details.algorithm.value = alg;		
	
	let dt_str = get_tag_value("timestep");
	if(dt_str == "") details.timestep = "";
	else{
		let dt = Number(dt_str);
		if(isNaN(dt_str) || dt <= 0){
			alert_import("For 'timestep' the value '"+dt_str+"' must be a positive number");
		}
		details.timestep = dt;
	}
	
	details.indmax = check_pos_integer("ind-max",INDMAX_DEFAULT);
	
	details.param_output_max = check_pos_integer("param-output-max",PARAM_OUTPUT_MAX_DEFAULT);
}


/// Applies the simulation command
function inference_command()
{
	let details = model.inf_details;
	
	let start = get_tag_value("start"); 
	if(start != "") is_number(start,"start");
	details.t_start = start;
	
	let end = get_tag_value("end"); 
	if(end != "") is_number(end,"end");
	details.t_end = end;
	
	if(details.t_start != "" && details.t_end != ""){
		if(Number(details.t_start) >= Number(details.t_end)){
			alert_import("The start time '"+start+"' must be before the end time '"+end+"'");
		}
	}
	
	details.seed_on = {value:"No"};
	details.seed = SEED_DEFAULT;
	let seed = get_tag_value("seed"); 
	if(seed != ""){ details.seed = seed; details.seed_on.value = "Yes";}
	
	let alg = get_tag_value("algorithm").toUpperCase();
	if(alg == "" || set_default) alg = ALG_DEFAULT;
	
	if(option_error("algorithm",alg,inf_alg_list) == true) return;
	
	details.algorithm.value = alg;	

	let dt_str = get_tag_value("timestep");
	if(dt_str == "") details.timestep = "";
	else{
		is_positive(dt_str,"timestep");
		details.timestep = Number(dt_str);
	}
	
	details.indmax = check_pos_integer("ind-max",INDMAX_DEFAULT);
	
	details.param_output_max = check_pos_integer("param-output-max",PARAM_OUTPUT_MAX_DEFAULT);
	
	switch(alg){
	case "DA-MCMC":
		details.nchain = check_pos_integer("nchain",MCMC_CHAIN_DEFAULT);
		details.sample = check_pos_integer("update",MCMC_SAMPLE_DEFAULT);
		details.output_param = check_pos_integer("param-output",MCMC_OP_PARAM_DEFAULT);
		details.output_state = check_pos_integer("state-output",MCMC_OP_STATE_DEFAULT);
		details.cha_per_core = check_pos_integer("chain-per-core",MCMC_CHAIN_PER_CORE_DEFAULT);
		break;
	
	case "PAS-MCMC":
		details.npart = check_pos_integer("npart",PAS_PART_DEFAULT);
		details.gen_update = check_percent("gen-percent",PAS_GEN_UPDATE_DEFAULT);
		if(details.gen_update == 0) alert_import("'gen-percent' cannot be 0");
		details.sample = check_pos_integer("update",MCMC_SAMPLE_DEFAULT);
		details.output_param = check_pos_integer("param-output",MCMC_OP_PARAM_DEFAULT);
		details.output_state = check_pos_integer("state-output",MCMC_OP_STATE_DEFAULT);
		details.part_per_core = check_pos_integer("part-per-core",PAS_PART_PER_CORE_DEFAULT);
		break;
		
	case "ABC":
		details.abcsample = check_pos_integer("sample",ABC_SAMPLE_DEFAULT);
		details.accfrac = check_zero_one("acc-frac",ABC_ACFRAC_DEFAULT);
		break;
		
	case "ABC-SMC":
		details.abcsample = check_pos_integer("sample",ABCSMC_SAMPLE_DEFAULT);
		details.accfrac = check_zero_one("acc-frac",ABCSMC_ACFRAC_DEFAULT);
		details.numgen = check_pos_integer("gen",ABCSMC_GEN_DEFAULT);
		details.kernelsize = check_pos_number("kernel-size",ABCSMC_KERNEL_DEFAULT);
		break;
		
	default: break;
		break;
	}
	
	details.burnin_frac = BURNIN_FRAC_DEFAULT;
	details.anneal_type = {te:ANNEAL_DEFAULT};
	details.anneal_rate = ANNEAL_RATE_DEFAULT;
	details.anneal_power = ANNEAL_POWER_DEFAULT;
	
	if(alg == "DA-MCMC" || alg == "PAS-MCMC"){
		let burnin_str = get_tag_value("burnin-percent"); 
		if(burnin_str != ""){
			let burnin = Number(burnin_str);
			if(isNaN(burnin) || burnin < 1 || burnin > 90){
				alert_import("For 'burnin-percent' the value '"+burnin_str+"' must be a number between 1 and 90");
				return;		
			}
			details.burnin_frac = burnin_str;
		}	
	}
	
	if(alg == "DA-MCMC"){
		let anneal_str = get_tag_value("anneal"); 
		if(anneal_str != ""){
			if(option_error("anneal",anneal_str,["none","scan","power-auto","log-auto","power"]) == true) return;
			
			details.anneal_type = {te:anneal_str};
			
			switch(anneal_str){
			case "scan":
				{
					let rate_str = get_tag_value("rate");
					if(rate_str != ""){		
						is_positive(rate_str,"rate");
						details.anneal_rate = rate_str;
					}
				}
				break;
				
			case "power":	
				{
					let power_str = get_tag_value("power");
					if(power_str != ""){
						is_positive(power_str,"power");
						details.anneal_power = power_str;
					}
				}
				break;
				
			default: break;
			}
		}
	}
	
	let ppc = model.ppc_details;

	if(ppc.ppc_t_start == ""){
		ppc.ppc_t_start = details.t_start;
		ppc.ppc_t_end = details.t_end;
	}
	
	let diag = get_tag_value("diagnostics").toLowerCase();
	if(diag != ""){
		if(option_error("diagnostics",diag,["on","off"]) == true) return;
	}
}


/// Applies the post-sim command
function post_sim_command()
{
	let details = model.ppc_details;
	
	// Copies information from inference
	let inf_de = model.inf_details;
	details.timestep = inf_de.timestep;
	details.t_start = inf_de.t_start;
	details.t_end = inf_de.t_end;
	details.inf_t_end = inf_de.t_end;
	
	let start = get_tag_value("start"); 
	if(start != "") details.ppc_t_start = start;
	else details.ppc_t_start = details.t_start;
		
	let end = get_tag_value("end"); 
	if(end != "") details.ppc_t_end = end;
	else details.ppc_t_end = details.t_end;
		
	if(details.ppc_t_start != "" && details.ppc_t_end != ""){
		if(Number(details.ppc_t_start) >= Number(details.ppc_t_end)){
			alert_import("The start time '"+start+"' must be before the end time '"+end+"'");
		}
		if(Number(end) > Number(details.t_end)) details.t_end = end;
	}
	
	let num = get_tag_value("number").toLowerCase();
	if(num == "") num = PPC_NUM_DEFAULT;
	if(set_default) num = PPC_NUM_DEFAULT;
	
	is_pos_int(num,"number");
	details.number = num;
	
	details.seed_on = {value:"No"};
	details.seed = SEED_DEFAULT;
	let seed = get_tag_value("seed"); 
	if(seed != ""){ details.seed = seed; details.seed_on.value = "Yes";}
	
	let alg = "gillespie";
	
	if(option_error("algorithm",alg,sim_alg_list) == true) return;
	details.algorithm.value = alg;	

	let resample = get_tag_value("resample");
	if(resample != ""){
		details.resample = resample.split(","); 		
		details.line = imp.line;	
	}
}


/// Check that a number is a positive integer
function check_pos_integer(te,def)
{
	let value = get_tag_value(te);
	
	if(value == "" || set_default) return def;
	
	let num = Number(value);
	if(isNaN(value) || num <= 0 || num != Math.floor(num)){
		alert_import("For '"+te+"' the value '"+value+"' must be a positive integer");
	}
	return num;
}


/// Check that a number is a percentage
function check_percent(te,def)
{
	let value = get_tag_value(te);
	
	if(value == "" || set_default) return def;
	
	let num = Number(value);
	if(isNaN(value) || num < 0 || num > 100){
		alert_import("For '"+te+"' the value '"+value+"' must be between 0 and 100");
	}
	return num;
}

/// Check that a number is positive 
function check_pos_number(te,def)
{
	let value = get_tag_value(te);
	
	if(value == "" || set_default) return def;
	
	let num = Number(value);
	if(isNaN(value) || num <= 0 || num >= 100){
		alert_import("For '"+te+"' the value '"+value+"' must be between 0 and 100, exclusive");
	}
	return num;
}


/// Check that a number is between zero and one
function check_zero_one(te,def)
{
	let value = get_tag_value(te);
	
	if(value == "" || set_default) return def;
	
	let num = Number(value);
	if(isNaN(value) || num <= 0 || num > 1){
		alert_import("For '"+te+"' the value '"+value+"' must be between zero and one");
	}
	return num;
}


/// Loads individual effects into the model
function ind_effect_command()
{
	let p = imp.p;

	if(p == undefined){
		alert_import("A species must be defined before an individual effect can be added");
	}
		
	let name = get_tag_value("name"); if(name == "") cannot_find_tag();
	
	let spl = name.split(",");
	let ie_list = []; 
	for(let i = 0; i < spl.length; i++){
		let nam = spl[i];
		check_name_input(nam,"Individual effect name",true);
		ie_list.push({name:nam});
	}
	
	let A = get_tag_value("A"); 
	let A_sparse = get_tag_value("A-sparse"); 
	let pedigree = get_tag_value("pedigree"); 
	
	let num = 0;
	if(A != "") num++;
	if(A_sparse != "") num++;
	if(pedigree != "") num++;
	
	if(num > 1){
		alert_import("Cannot specify more than one of 'A', 'A-sparse' and 'pedigree'."); return;
	}
	
	let A_matrix = {check:false, loaded:false, A_value:[], ind_list:[], pedigree:false, sire_list:[], dam_list:[]};
	
	if(num == 1){
		if(pedigree != ""){
			A_matrix = {check:true, loaded:true, A_value:[], ind_list:[], pedigree:true, sire_list:[], dam_list:[]};
	
			let tab_ped = load_table(pedigree.te,true,pedigree.sep,pedigree.name);
				
			if(tab_ped.ncol != 3){
				alert_import("'pedigree' file should have 3 columns"); return;
			}
				
			if(tab_ped.heading[0] != "ID" || tab_ped.heading[1] != "sire" || 
					tab_ped.heading[2] != "dam"){
				alert_import("'pedigree' file should have columns with headings 'ID', 'sire' and 'dam'"); return;
			}
			
			for(let r = 0; r < tab_ped.nrow; r++){
				A_matrix.ind_list.push(tab_ped.ele[r][0]);
				A_matrix.sire_list.push(tab_ped.ele[r][1]);
				A_matrix.dam_list.push(tab_ped.ele[r][2]);
			}					
		}
		
		if(A_sparse != ""){
			A_matrix = {check:true, loaded:true, pedigree:false, sire_list:[], dam_list:[]};
			
			let ind_list = get_tag_value("ind-list"); 
		
			let tab_ind = load_table(ind_list.te,true,ind_list.sep,ind_list.name);
			if(tab_ind == undefined) return;
			
			let list=[];
			for(let r = 0; r < tab_ind.nrow; r++){ 
				list.push(tab_ind.ele[r][0]);
			}			
		
			A_matrix.ind_list = list;
		
			let N = list.length;
			let val=[];
			for(let r = 0; r < N; r++){
				val[r]=[];
				for(let c = 0; c < N; c++) val[r][c] = 0;
			}
		
			let tab = load_table(A_sparse.te,true,A_sparse.sep,A_sparse.name);
			if(tab == undefined) return;
		
			for(let r = 0; r < tab.nrow; r++){
				let fl = false;
				let j = Number(tab.ele[r][0]); if(isNaN(j) || j != Math.floor(j)) fl = true;
				let i = Number(tab.ele[r][1]); if(isNaN(i) || i != Math.floor(i)) fl = true;
				let value = Number(tab.ele[r][2]); if(isNaN(value)) fl = true;
				if(fl == true){
					alert_import("Problem loading table on line "+to_string(r+1)); return;
				}
				val[j][i] = value;
				val[i][j] = value;
			}
			A_matrix.A_value = val;
		}
	
		if(A != ""){
			A_matrix = {check:true, loaded:true, pedigree:false, sire_list:[], dam_list:[]};
			
			let tab = load_table(A.te,true,A.sep,A.name);
			//let tab = load_table_from_file(A);
			if(tab == undefined) return;
			if(typeof tab == 'string') alert_import(tab); 

			A_matrix.ind_list = tab.heading;
			if(tab.nrow != tab.ncol){
				alert_import("The file '"+tab.filename+"' must contain an equal number of columns and rows."); 
			}
			
			let val = [];
			for(let r = 0; r < tab.nrow; r++){
				val[r] = [];
				for(let c = 0; c < tab.ncol; c++){
					let ele = tab.ele[r][c];
					if(isNaN(ele)){
						alert_import(in_file_text(tab.filename)+" the element '"+ele+"' is not a number (row "+(r+2)+", col "+(c+1)+")");
					}
					val[r][c] = Number(ele);
				}
			}
			A_matrix.A_value = val;
		}
	}
	
	// Checks that ind effects do not already exist in the model
	for(let pp = 0; pp < model.species.length; pp++){
		let iegroup = model.species[pp].ind_eff_group;
		for(let k = 0; k < iegroup.length; k++){
			let list = iegroup[k].ie_list;
			for(let ii = 0; ii < list.length; ii++){				
				for(let j = 0; j < ie_list.length; j++){	
					if(find(list,"name",ie_list[j].name) != undefined){
						alert_import("Individual effect '"+ie_list[j].name+"' already in the model");
					}
				}
			}
		}
	}
	
	model.species[p].ind_eff_group.push({ie_list:ie_list, A_matrix:A_matrix});
}


/// Displays text referencing a file/table
function in_file_text(te)
{
	if(te == "inline") return "In table '[[...]]'";
	return "In file '"+te+"'";
}


/// Loads fixed effects into the model
function fixed_effect_command()
{
	let p = imp.p;

	if(p == undefined){
		alert_import("A species must be defined before a fixed effect can be added");
	}
		
	let name = get_tag_value("name"); if(name == "") cannot_find_tag();
	check_name_input(name,"Individual fixed effect",true);

	let X_vector = { loaded:false, ind_list:[], X_value:[]};

	let X = get_tag_value("X"); 
	if(X != ""){
		X_vector.loaded = true;		
		let tab = load_table(X.te,true,X.sep,X.name);		
		if(tab == undefined) return;
		
		let col_name = ["ID","value"]; 
		let subtab = get_subtable(tab,col_name);
		if(subtab.error != "") alert_import(subtab.error);
		
		for(let r = 0; r < subtab.nrow; r++){
			X_vector.ind_list.push(subtab.ele[r][0]);
			let val = tab.ele[r][1];
			if(isNaN(val)){
				alert_import(in_file_text(tab.filename)+" the element '"+val+"' is not a number (row "+(r+2)+")");
			}
			X_vector.X_value.push(val);
		}
	}
	
	model.species[p].fix_eff.push({name:name, X_vector:X_vector, defined:true});
}


/// Converts dates to times
function convert_date_time(so)
{
	if(so.time_format != "float"){
		let sep;
		if(so.time_format.includes("/")) sep = "/";
		if(so.time_format.includes(".")) sep = ".";
		if(so.time_format.includes("-")) sep = "-";
		let temp = so.time_format.split(sep);
		
		for(let c = 0; c < so.load_col.length; c++){
			if(so.load_col[c].heading == "t"){
				let tab = so.table;
				for(let r = 0; r < tab.nrow; r++){
					let el = tab.ele[r][c];
					let val = el.split(sep);
					
					let warn = false;
					
					if(temp.length != val.length) warn = true;
					
					let year, month, day;
					for(let k = 0; k < temp.length; k++){
						if(isNaN(val[k])) warn = true; 
						else{
							if(temp[k].substr(0,1) == "d") day = Number(val[k]);
							if(temp[k].substr(0,1) == "m") month = Number(val[k]);
							if(temp[k].substr(0,1) == "y"){
								year = Number(val[k]);
							}
						}
					}
					
					if(warn == true){
						alert_import("For '"+tab.filename+"': '"+el+"' does not follow format '"+so.time_format+"' (col '"+tab.heading[c]+"', row "+(r+2)+")");
					}
				
					let time = (new Date(year, month, day)).getTime();
					let time_beg = (new Date(year, 1, 1)).getTime();
					let time_end = (new Date(year+1, 1, 1)).getTime();
					let frac = (time - time_beg)/(time_end-time_beg);
					
					tab.ele[r][c] = String((year+frac).toFixed(4));
				}
			}
		}
	}
}


/// Sets the columns which need to be loaded SNP data
function set_SNP_columns(tab,ds)
{
	let snp_root = ds.spec.snp_root;
	let len = snp_root.length;	
	ds.load_col.pop();

	for(let c = 0; c < tab.ncol; c++){
		let head = tab.heading[c];
		if(head.length >= len && head.substr(0,len) == snp_root){
			ds.load_col.push({heading:head, desc:"SNP data",type:"text"});
		}
	}
}


/// Sets the columns which need to be loaded matrix data
function set_genetic_matrix_columns(ele,ds)
{
	//let obs_col = ds.load_col[ds.load_col.length-1];
	let c = ds.load_col.length-1;
	for(let r = 0; r < ele.length; r++){
		ds.load_col.push({heading:ele[r][c], desc:"Matrix data",type:"positive int"});
	}
}


/// Loads up the observation model
function load_obs_model(spec)
{
	let error = get_tag_value("error");
	if(error == "") cannot_find_tag();
			
	let spl = error.split(":");

	switch(spl[0].toLowerCase()){
	case "normal": spec.obs_error_dist = {value:"Normal"}; break;
	case "poisson": spec.obs_error_dist = {value:"Poisson"}; break;
	case "neg-binomial": spec.obs_error_dist = {value:"Negative binomial"}; break;
	default:
		alert_import("For 'error' the observation error '"+spl[0]+"' not recognised. It must be chosen from the following options: 'normal', 'poisson' or 'neg-binomial'.");
		break;
	}
	
	spec.p = "0.5";
	spec.percent = "10";
	spec.sd = "";
	spec.obs_error = {value:"percent"};	
	spec.obs_error_p = {value:"p"};	

	switch(spec.obs_error_dist.value){
	case "Normal":
		{
			if(spl.length != 2) alert_import("For 'error' problem with expression '"+error+"'");
	
			let val = spl[1].trim();
		
			if(val.toLowerCase() == "file"){
				spec.obs_error = {value:"file"};	
			}
			else{
				if(is_percent(val) == true){
					spec.obs_error = {value:"percent"};
					spec.percent = val.substr(0,val.length-1);
					if(isNaN(spec.percent)){
						alert_import("For 'error' the expression '"+spec.percent+"' is not a percentage");
					}
					
					if(Number(spec.percent) <= 0){
						alert_import("For 'error' the value '"+spec.percent+"%' must be a positive percentage");
					}
				}	
				else{
					spec.obs_error = {value:"sd"};
					spec.sd = val;
					
					if(isNaN(spec.sd)){
						alert_import("For 'error' the expression '"+spec.sd+"' must be a number or a percentage");
					}
					
					if(Number(spec.sd) <= 0){
						alert_import("For 'error' the value '"+spec.sd+"' must be positive");
					}
				}
			}
		}
		break;
		
	case "Poisson":
		if(spl.length != 1) alert_import("For 'error' problem with expression '"+error+"'");
		break;
	
	case "Negative binomial":
		{
			if(spl.length != 2) alert_import("For 'error' problem with expression '"+error+"'");
	
			let val = spl[1].trim();
		
			if(val.toLowerCase() == "file"){
				spec.obs_error_p = {value:"file"};
			}
			else{
				spec.obs_error_p = {value:"p"};
				spec.p = val;
					
				if(isNaN(spec.p)){
					alert_import("For 'error' the expression '"+spec.p+"' must be a number");
				}
				
				if(Number(spec.p) <= 0 || Number(spec.p) >= 1){
					alert_import("For 'error' the value '"+spec.p+"' must be between zero and one");
				}
			}
		}
		break;
	}
}


/// Gets a filter from a string
function get_filt_from_str(filter,p)
{
	let filt = load_default_filt(p);
	
	if(filter != ""){
		let sp = model.species[p];
		
		let spl = filter.split(",");
		for(let j = 0; j < spl.length; j++){
			let spl2 = spl[j].split("=");
		
			if(spl2.length != 2){ 
				alert_import("For 'filter' error understanding '"+filter+"'"); 
				return;
			}
			
			if(spl2[0] == ""){
				alert_import("For 'filter' the value '"+filter+"' does not specify a classification");
				return;
			}
			
			let cl = find_nocase(sp.cla,"name",spl2[0]);
			if(cl == undefined){ 
				alert_import("For 'filter' the value '"+spl2[0]+"' is not a classification"); 
				return;
			}
			let claa = sp.cla[cl];
			
			if(spl2[1].includes(":")){                  // Load up an observation model
				filt.cla[cl].radio.value = "ObsMod";
				
				//let spl3 = spl2[1].split("|");
				let spl3 = split_with_bracket(spl2[1],"|");
				for(let k = 0; k < spl3.length; k++){
					let spl4 = spl3[k].split(":");
				
					if(spl4[0] == ""){
						alert_import("In 'filter' there was an error");
						return;
					}
						
					let c = hash_find(claa.hash_comp,spl4[0]);
					if(c == undefined){ 
						alert_import("For 'filter' the value '"+spl3[k]+"' is not a compartment in '"+spl2[0]+"'"); 
						return;
					}
					
					switch(spl4.length){
					case 1: filt.cla[cl].comp[c] = true; break;
					case 2: filt.cla[cl].comp[c].prob_eqn.te = spl4[1]; break;
					default: alert_import("In 'filter' there was an error"); return;
					}
				}
			}
			else{
				if(spl2[1].toLowerCase() == "file"){
					filt.cla[cl].radio.value = "File";
				}
				else{
					filt.cla[cl].radio.value = "Comp";
				
					//let spl3 = spl2[1].split("|");
					let spl3 = split_with_bracket(spl2[1],"|");
				
					for(let k = 0; k < spl3.length; k++){
						if(spl3[k] == ""){
							alert_import("For 'filter' the classification '"+spl2[0]+"' is not set");
							return;
						}
						
						let c = hash_find(claa.hash_comp,spl3[k]);
						if(c == undefined){ 
							alert_import("In 'filter' the value '"+spl3[k]+"' is not a compartment in '"+spl2[0]+"'"); 
							return;
						}
						
						filt.cla[cl].comp[c].check = true;
					}
				}
			}
		}
	}
	
	return filt;
}


/// Sets the transition filter from a string (e.g. S->E or S->E:0.5|E->I:x)
function set_get_tra_filt_from_str(claa_sel,filt,name)
{
	filt.te = claa_sel.name;

	let tra_fi = [];
	for(let i = 0; i < claa_sel.ntra; i++){
		let tra = claa_sel.tra[i];
		tra_fi[i] = { tra_name_store:tra.name, check:false, prob_eqn:create_equation("0","trans_prob")};
	}
	
	filt.trans_obs_model = {value: "off"};
	if(name.includes(":")) filt.trans_obs_model.value = "on";
	
	//let spl = name.split("|");
	let spl = split_with_bracket(name,"|");
	for(let j = 0; j < spl.length; j++){
		let spl2 = spl[j].split(":");
		
		switch(spl2.length){
		case 1:
			{
				let i = hash_find(claa_sel.hash_tra,spl2[0]);
				if(i == undefined) alert_import("The transition '"+spl2[0]+"' does not exist");
				tra_fi[i].prob_eqn = create_equation("1","trans_prob");
				tra_fi[i].check = true;
			}
			break;
		
		case 2:
			{
				let i = hash_find(claa_sel.hash_tra,spl2[0]);
				if(i == undefined) alert_import("The transition '"+spl2[0]+"' does not exist");
				tra_fi[i].prob_eqn = create_equation(spl2[1],"trans_prob");
				tra_fi[i].check = true;
			}
			break;
		
		default:
			alert_import("Expression '"+spl[j]+"' not understood");
			break;
		}
	}

	filt.tra = tra_fi;
}


/// Reads in parameter samples into sim_results
function sim_param_command()
{
	let file = get_tag_value("file");
	if(file == "") cannot_find_tag();
	
	read_param_samples_file(0,file,sim_result);
}


/// Gets the chain value (increments for plotting)
function get_chain_value()
{
	let chain = get_tag_value("chain");
	if(chain == "") chain = "0";
	return String(Number(chain)+1);
}


/// Reads in parameter samples into inf_results
function inf_param_command()
{
	let file = get_tag_value("file");
	if(file == "") cannot_find_tag();
	
	let chain = get_chain_value();

	read_param_samples_file(chain,file,inf_result);
}


/// Reads in parameter samples into inf_results
function inf_generation_command()
{
	let file = get_tag_value("file");
	if(file == "") cannot_find_tag();
	
	read_param_samples_file("gen-plot",file,inf_result);
}


/// Reads in diagnostics into inf_results
function inf_diagnostics_command()
{
	let file = get_tag_value("file");
	if(file == "") cannot_find_tag();
	
	let chain = get_chain_value();
	
	inf_result.diagnostics[chain] = file.te;
	inf_result.prop_diag_on = true;
}


/// Reads in parameter samples into ppc_results
function post_sim_param_command()
{
	let file = get_tag_value("file");
	if(file == "") cannot_find_tag();
	
	read_param_samples_file(0,file,ppc_result);
}


/// Reads in state samples into sim_results
function sim_state_command()
{
	let file = get_tag_value("file");
	if(file == "") cannot_find_tag();
	
	read_state_samples_file(0,file,sim_result);
}


/// Reads in state samples into inf_results
function inf_state_command()
{
	let file = get_tag_value("file");
	if(file == "") cannot_find_tag();
	
	let chain = get_chain_value();
	
	read_state_samples_file(chain,file,inf_result);
}


/// Reads in state samples into sim_results
function post_sim_state_command()
{
	let file = get_tag_value("file");
	if(file == "") cannot_find_tag();
	
	read_state_samples_file(0,file,ppc_result);
}


/// Reads in transition diagnostic information
function trans_diag_command()
{
	let file = get_tag_value("file");
	if(file == "") cannot_find_tag();
	
	read_trans_diag(file,inf_result);
}

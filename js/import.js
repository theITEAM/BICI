// Imports information into BICI

"use strict";

/// Processes a given command
function process_command(cname,tags,loop)
{
	inter.imp.tags = tags;
	inter.imp.cname = cname;
	//if(cname != "param") inter.imp.need_init_param = true;
 
	switch(cname){
	case "species": species_command(loop); break;
	case "classification": classification_command(loop); break;
	case "set": set_command(); break;
	case "camera": camera_command(); break;
	case "compartment": compartment_command(); break;
	case "transition": transition_command(cname); break;
	case "source": transition_command(cname); break;
	case "sink": transition_command(cname); break;
	case "data-dir": datadir_command(); break;
	case "description": description_command(); break;
	case "label": label_command(); break;
	case "box": box_command(); break;
	case "param": param_command(); break;
	case "derived": derived_command(); break;
	case "simulation": simulation_command(); break;
	case "inference": inference_command(); break;
	case "ind-effect": ind_effect_command(); break;
	case "fixed-effect": fixed_effect_command(); break;
		
	case "init-pop": case "add-ind": case "remove-ind": case "move-ind": 
	case "init-pop-sim": case "add-ind-sim": case "remove-ind-sim": case "move-ind-sim": 
	case "init-pop-prior": case "comp-data": case "trans-data":
	case "source-data": case "sink-data":
	case "test-data": case "pop-data": case "pop-trans-data": 
	case "set-traps-data": case "ind-trapped-data":
	case "genetic-data":
		import_data_table_command(cname);
		break;
	
	case "sim-param": sim_param_command(); break;
	case "sim-state": sim_state_command(); break;
	case "inf-param": inf_param_command(); break;
	case "inf-state": inf_state_command(); break;
	
	//case "sim-const": sim_const_command(); break;
	//case "inf-const": inf_const_command(); break;
	
	default: 
		error("'"+cname+"' not recognised");
		alert_import("Command '"+cname+"' not recognised"); 
		return;
	}
	
	let warn = check_initial_pop_error(false);
	if(warn != undefined){ alert_import(warn); return;}
}


/// Loads up an imported file (for testing purposes only)
function load_import_test()
{	
	let fs = require('fs');	
	//let te = fs.readFileSync("Testing/init.txt", 'utf8');
	let te = fs.readFileSync("Execute/init.txt", 'utf8');

	import_file(te);
}


/// Import a script to define (or partially define) model and data 
function import_file(te)                                
{
	let store = save_BICI_file();     // Saves so that can be reset if import unsuccessful
	
	model.start_new();
	
	te = add_escape_char(te);
	
	let	lines = te.split('\n');

	let formatted = get_formatted(lines);

	close_help();

	let p, cl;
	if(model.species.length == 1){
		p = 0;
		if(model.species[p].ncla == 1) cl = 0;
	}  
	
	// Import happens in four stages:
	// (0) Load species and classification information
	// (1) Load parameter information
	// (2) Everything else apart from output samples
	// (3) Output samples
	
	// Keeps track of the current species and classification 
	inter.imp = { datadir:undefined, script:formatted, previous_loaded_table:[]};   
		
	for(let loop = 0; loop < 4; loop++){
		inter.imp.p = p; inter.imp.cl = cl;
		
		if(loop == 3){
			init_param();

			results_add_model(sim_result,model.sim_details,"sim");
			results_add_model(inf_result,model.inf_details,"inf");
		}
			
		for(let j = 0; j < lines.length; j++){
			let trr = lines[j].trim();
		
			let flag = false;                    // Ignores line if empty or commented out
			if(trr.length == 0) flag = true;
			if(trr.length >= 2 && trr.substr(0,2) == "//") flag = true;
			if(trr.length >= 1 && trr.substr(0,1) == "#") flag = true;
			
			if(flag == false){
				inter.imp.line = j;
			
				let line = get_command_tags(trr);
			
				if(check_import_error(store) == true) return; 
				
				let cname = line.type;
				
				switch(cname){                   // Accounts for shortened versions of commands
				case "comp": cname = "compartment"; break;
				case "trans": cname = "transition"; break;
				case "class": cname = "classification"; break;
				case "desc": cname = "description"; break;
				case "parameter": cname = "param"; break;
				case "der": cname = "derived"; break;
				case "sim": cname = "simulation"; break;
				case "inf": cname = "inference"; break;
				}
	
				if(check_import_error(store) == true) return; 
				
				let process = true;
				switch(loop){
				case 0:   // In the first pass create species, classifications, compartments and loads data directory
					switch(cname){
					case "species": case "classification": case "comp": case "compartment": case "data-dir": break;
					default: process = false; break;
					}
					break;
					
				case 1:   // The next pass does loads the parameters 
					process = false; 
					switch(cname){
					case "species": case "classification": case "param": case "derived": 
						process = true; 
						break;
					}
					break;
					
				case 2:   // The last pass does everything except the parameters and compartments
					switch(cname){
					case "comp": case "compartment": case "param": case "derived": case "data-dir":
					case "sim-param": case "sim-state":
					case "inf-param": case "inf-state": 
					//case "sim-const": case "inf-const": 
						process = false; break;
					}
					break
					
				case 3:
					process = false; 
					switch(cname){
					case "sim-param": case "sim-state":
					case "inf-param": case "inf-state": 
						process = true; 
						break;
					}
					break;
				}
				
				if(process == true){
					process_command(cname,line.tags,loop);
					if(check_import_error(store) == true) return; 
		
					for(let n = 0; n < line.tags.length; n++){
						if(line.tags[n].done != 1){ 
							alert_import("Tag '"+line.tags[n].name+"' not used"); break;
						}
					}
				
					if(check_import_error(store) == true) return; 
				}
			}
		}
	}
	let warn = check_initial_pop_error(true);
	if(warn != undefined) alert_import(warn);
	
	if(check_import_error(store) == true) return; 

	check_param_complete(); // Checks information about all parameters have been loaded correctly
 
	if(check_import_error(store) == true) return; 

	model_set_auto_positions();
	
	if(check_import_error(store) == true) return; 
	
	check_import_correct();
	
	if(check_import_error(store) == true) return; 
	
	setup_camera();
	
	model.update_pline_all();

	model.load_default_map();

	if(sim_result.on == true){
		intialise_plot_filters(sim_result);
	}
	
	if(inf_result.on == true){
		intialise_plot_filters(inf_result);
	}
	
	initialise_pages();
	
	change_page({pa:"Model"});
}


/// Imports a geojson file
function import_geojson(file)  
{
	if(check_char_allowed(file,"/<>:\"\\|?*") != "success") return; 
	
	file = inter.imp.datadir+"/"+file;
	
	if(check_file_exists(file) != "success") return;
		
	let tab = previously_loaded(file);

	if(tab != undefined) return tab;

	let fs = require('fs');	

	let te = fs.readFileSync(file, 'utf8');
	
	let da = JSON.parse(te);
	
	let feature = get_feature(da);
	if(feature.length == 0){
		return "There are no features in the file '"+file+"'";
	}	
		
	tab = get_feature_table(da,feature);
	add_previously_loaded(file,tab);

	return tab;
}


/// Loads up a table in table format
function import_table(file)                                // Imports a table of data
{
	if(check_char_allowed(file,"/<>:\"\\|?*") != "success") return; 
	
	file = inter.imp.datadir+"/"+file;
	
	if(check_file_exists(file) != "success") return;
	
	let fs = require('fs');
	
	let te = fs.readFileSync(file, 'utf8');
	
	let lines = te.split('\n');

	let heading = lines[0].split(',');
	let ncol = heading.length;
	
	let ele = [];
	for(let j = 1; j < lines.length; j++){
		if(lines[j] != ""){
			ele[j-1] = lines[j].split(',');
		}
	}
	
	let nrow = ele.length;
	if(nrow > 0){
		for(let j = 0; j < nrow; j++){
			if(ele[j].length != ncol){
				alert_import("The rows do not all have the same size:"+ele[j].length+" "+ncol);
				return;
			}
		}
	}
	else{ 
		alert_import("Table does not contain any data"); return;
	}

	let tab = { heading:heading, ele:ele, nrow:nrow, ncol:ncol};
	add_previously_loaded(file,tab)

	return tab;	
}


/// If a table has been previously loaded then use that
function previously_loaded(file)
{
	let i = find(inter.imp.previous_loaded_table,"file",file);
	if(i != undefined) return inter.imp.previous_loaded_table[i].tab;
}


/// Adds a table to those previously loaded
function add_previously_loaded(file,tab)
{
	inter.imp.previous_loaded_table.push({file:file, tab:tab});
}


/// Checks if lat and lng are correctly specified
function check_latlng(lat,lng)
{
	if(lng == "" || lat == ""){
		alert_import("The 'lat' and 'lng' must be set"); return;
	}
	if(isNaN(lat)){
		alert_import("'lat' must be a number"); return;
	}

	if(isNaN(lng)){
		alert_import("'lng' must be a number"); return;
	}

	if(lat > 90 || lat < -90){
		alert_import("'lat' must be in the range -90 to 90"); return;
	}
	
	if(lng > 180 || lng < -180){
		alert_import("'lng' must be in the range -180 to 180"); return;
	}
	
	return "success";
}


/// Loads up data
function load_data()
{
	if(inter.imp.datadir == undefined){
		alert_import("The data directory must be set for files can be loaded.");
		return;
	}
	
	let type = get_tag_value("type"); if(type == ""){ cannot_find_tag(); return;}
	let file = get_tag_value("file"); if(file == ""){ cannot_find_tag(); return;}
	
	let p = inter.imp.p;
	if(p == undefined){
		alert_import("Population needs to be specified before "+cname); return;
	}
	
	let name = get_tag_value("name"); if(name == ""){ cannot_find_tag(); return;}

	let tab = import_table(file);
	
	if(tab == 1){ alert_import(errormsg); return;}
	
	let comp_flag = false;
	
	let k = 0; while(k < data_template.length && data_template[k].type != type) k++;
	if(k == data_template.length){
		alert_import("The data type '"+type+"' is not recognised"); return;
	}
		
	let dt = data_template[k];
	
	let collist=[];
	for(let i = 0; i < dt.cols.length; i++){
		let tag = get_tag_value(dt.cols[i]); if(tag == ""){ cannot_find_tag(); return;}
		collist.push({name:dt.cols[i], c:find_column(tab,tag)});
	}
	
	if(dt.comp_flag == true){
		let sp = model.species[p]; 
		
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			var col_name = get_tag_value(claa.name); 
			if(col_name == ""){ cannot_find_tag(); return;}
			collist.push({name:claa.name, c:find_column(tab,col_name)});
		}
	}

	let source = {type:type, name:name, values:[]};

	if(dt.from_to_flag == true){
		let from = get_tag_value("from"); if(from == ""){ cannot_find_tag(); return;}
		let to = get_tag_value("to"); if(from == ""){ cannot_find_tag(); return;}

		let flag = false;
		
		let sp = model.species[p];
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			for(let tr = 0; tr < claa.ntra; tr++){	
				let traa = claa.tra[tr];
				if((traa.i >= 0 && claa.comp[traa.i].name == from) &&
					 (traa.f >= 0 && claa.comp[traa.f].name == to)) flag = true;
			}
		}
		if(flag == false){ alert_import("Could not find transition '"+from+"→"+to+"'"); return;}
		
		source.from = from;
		source.to = to;
	}
	
	if(dt.class_flag == true){
		let clas = get_tag_value("classification"); if(clas == ""){ cannot_find_tag(); return;}
		source.classification = clas;
	}
	
	if(dt.min_max_flag == true){
		let min = get_tag_value("min"); if(min == ""){ cannot_find_tag(); return;}
		let minnum = Number(min);
		if(isNaN(minnum)){ alert_import("For 'min' the value '"+min+"' is not a number"); return;} 
		source.min = minnum;
		
		let max = get_tag_value("max"); if(max == ""){ cannot_find_tag(); return;}
		let maxnum = Number(max);
		if(isNaN(maxnum)){ alert_import("For 'max' the value '"+max+"' is not a number"); return;} 
		source.max = maxnum;
	}

	for(let r = 0; r < tab.nrow; r++){
		let row = {};
		for(let i = 0; i < collist.length; i++){
			row[collist[i].name] = tab.ele[r][collist[i].c];
		}
		source.values.push(row);
	}

	data.source.push(source);
}


/// Gets all the command tags from a given line
function get_command_tags(trr)
{
	let num_quote = 0;
	let quote_pos = [];
	
	for(let i = 0; i < trr.length; i++){
		if(trr.substr(i,1) == '"'){ num_quote++; quote_pos.push(i);}
	}
	if(num_quote%2 != 0){ alert_import("Syntax error: Quotes do not match up."); return;}
	
	for(let i = 0; i < num_quote; i+= 2){
		if(quote_pos[i]+1 == quote_pos[i+1]){
			alert_import("Syntax error: No content within the quotation marks."); return;
		}
	}
	
	let frag=[]; 
	
	let i = 0; 
	let quote = 0;
	let after_eq = false;
	
	do{
		let ist = i; i++;
		
		while(i < trr.length && !(trr.substr(i,1) == "=" && quote == 0) && 
				!(quote == 0 && trr.substr(i,1) == '"')&& 
				!(quote == 1 && trr.substr(i,1) == '"') && 
				!(quote == 0 && trr.substr(i,1) == " ") &&
				!(quote == 0 && after_eq == true)){
			i++;
		}
	
		after_eq = false;
		if(trr.substr(i,1) == "=") after_eq = true;
		
		let te = trr.substr(ist,i-ist);
		while(te.substr(0,1) == "\n") te = te.substring(1);
		while(te.substr(te.length-1,1) == "\n") te = te.substring(0,te.length-1);
		
		te = te.trim(); if(quote == 0) te = te.toLowerCase();
		if(te != ""){
			let pos = ist; while(pos < i && trr.substr(pos,1) == " ") pos++;
			frag.push({text:trr.substr(ist,i-ist).trim(), pos:pos, pos_end:i, quote:quote});
		}
		while(i < trr.length && quote == 0 && trr.substr(i,1) == " ") i++;
		if(i < trr.length){
			if(trr.substr(i,1) == '"'){ 
				quote = 1-quote; i++;
			}
		}
	}while(i < trr.length);
	
	if(frag.length == 0){ alert_import("Does not contain any content"); return;}
	
	if(frag[0].quote == 1){ alert_import("Should not start with quotes"); return;}
	let type = frag[0].text, type_pos = frag[0].pos; 
	
	if(find_in(command_list,type) == undefined){ 
		alert_import("Command '"+type+"' not recognised."); return;
	}
	
	inter.imp.typest = type;
	
	let num = (frag.length-1)/3;

	if(num != Math.floor(num)){ alert_import("Syntax error"); return;}
	let tags=[];

	for(let n = 0; n < num; n++){
		if(frag[1+3*n+0].quote != 0){ alert_import("Syntax error type 2"); return;}
		if(frag[1+3*n+1].text != "="){ alert_import("Syntax error type 3"); return;}
		if(frag[1+3*n+0].text == ""){ alert_import("Syntax error type 4"); return;}
		if(frag[1+3*n+2].text == ""){ alert_import("Property "+frag[1+3*n+0].text+" must have content"); return;}
		let fr = frag[1+3*n+0];
		tags.push({name:fr.text, pos:fr.pos, pos_end:fr.pos_end, value:frag[1+3*n+2].text, done:0});
	}
	
	for(let n = 0; n < tags.length-1; n++){
		for(let nn = n+1; nn < tags.length; nn++){
			if(tags[n].name == tags[nn].name){ alert_import("The tag '"+tags[n].name+"' is set more than once"); return;}
		}
	}
	
	return {type:type, type_pos, tags:tags};
}


/// Finds a column in a table
function find_column(tab,colname)
{
	let i = 0; while(i < tab.ncol && tab.heading[i] != colname) i++;
	if(i == tab.ncol)	alert_import("The table does not have the column '"+colname+"'");
	
	return i;
}


/// Gets the value of a tag
function get_tag_value(st)                                   
{
	let tags = inter.imp.tags;
	inter.imp.tagst = st;
	for(let i = 0; i < tags.length; i++){
		if(tags[i].name == st){ tags[i].done = 1; return tags[i].value;}
	}
	return "";
}


/// Error massage if a tag cannot be found
function cannot_find_tag()                                   
{
	alert_import("Cannot find the '"+inter.imp.tagst+"' tag for '"+inter.imp.typest+"'");
}


/// Error message for imported file (for specific line)
function alert_import(st)                                     
{
	let te = "On line "+(inter.imp.line+1)+": "+st
	if(te.length > 0 && te.substr(te.length-1,1) != ".") te += ".";
	
	inter.help = { title:"Error importing file!", te:te, st:st, line:inter.imp.line, scroll_to_line:true, script:inter.imp.script};
	
	generate_screen();
}


/// Alerts a specific line in the import file
function alert_line(st,line)                                     
{
	let te = "On line "+(line+1)+": "+st
	if(te.length > 0 && te.substr(te.length-1,1) != ".") te += ".";
	
	inter.help = { title:"Error importing file!", te:te, line:line, scroll_to_line:true, script:inter.imp.script};
	
	generate_screen();
}


/// Alerts a general error in the import file
function alert_noline(st)                                     
{
	let te = st;
	if(te.length > 0 && te.substr(te.length-1,1) != ".") te += ".";
	
	inter.help = { title:"Error importing file!", te:te, scroll_to_line:false, script:inter.imp.script};
	
	generate_screen();
}


/// Checks if an error occurs 
function check_for_error()
{
	let te = inter.help.te;
	if(te != undefined && te.substr(0,7) != "On line"){ alert_import(te); return true;}

	return false;  
}


/// Corrects captialisation deviation
function correct_caps(te,pos)
{
	for(let i = 0; i < pos.length; i++){
		if(te.toLowerCase() == pos[i].toLowerCase()) te = pos[i];
	}
	return te;
}
 
 
/// Checks is a particular value is one of a list of possible options
function option_error(na,te,pos)
{
	for(let i = 0; i < pos.length; i++) if(pos[i] == te) return false;
	
	let st = "'"+na+"' has a value '"+te+"' but should be chosen from one of the following options: ";
	for(let i = 0; i < pos.length; i++){
		if(i != 0) st += ", ";
		st += "'"+pos[i]+"'";
	}
	
	alert_import(st);
	return true;
}


/// Gets the classification object based on store values for p and cl
function get_claa()
{
	let p = inter.imp.p;
	if(p == undefined){
		alert_import("Population needs to be specified before "+inter.imp.cname); return;
	}
			
	let cl = inter.imp.cl;
	if(cl == undefined){
		alert_import("Classification need to be specified before "+inter.imp.cname); return;
	}
	return model.species[p].cla[cl];
}		


/// Adds colour formatting to the script
function get_formatted(lines)
{
	let formatted = "";
	for(let j = 0; j < lines.length; j++){
		let trr = lines[j].trim();
	
		let flag = false;                    // Ignores line if empty or commented out
		if(trr.length == 0) flag = true;
		if(trr.length >= 2 && trr.substr(0,2) == "//") flag = true;
		if(trr.length >= 1 && trr.substr(0,1) == "#") flag = true;
			
		if(flag == true) trr = "<g>"+trr+"</g>";
		else{
			let line = get_command_tags(trr);

			if(trr.substr(0,1) != "#" && trr.substr(0,2) != "//" && line != undefined){
				for(let j = line.tags.length-1; j >= 0; j--){
					let k = line.tags[j].pos;
					let k2 = line.tags[j].pos_end;
				
					trr = trr.substr(0,k2)+"</t>"+trr.substr(k2);
					trr = trr.substr(0,k)+"<t>"+trr.substr(k);
				}
				
			
				let kk = 1;  // Pads space around certain characters
				for(let kk = 1; kk < trr.length-1; kk++){
					let ch = trr.substr(kk,1);
					/*
					if(ch == "="){
						if(trr.substr(kk-1,1) != " ") trr = trr.substr(0,kk)+" "+trr.substr(kk);
						if(trr.substr(kk+1,1) != " ") trr = trr.substr(0,kk+1)+" "+trr.substr(kk+1);
						kk++
					}
					*/
				}
		
				let k = line.type_pos;
				trr = trr.substr(0,k)+"<c>"+line.type+"</c>" + trr.substr(k+line.type.length);
			}
		}
		
		formatted += trr+"\n";
	}
	return formatted;
}
	

/// Sets the camera if not specified in the input file
function setup_camera()
{
	for(let p = 0; p < model.species.length; p++){ 	
		for(let cl = 0; cl < model.species[p].ncla; cl++){
			let cam = model.species[p].cla[cl].camera;
			if(cam.scale == "") set_camera(p,cl);
			else set_ruler(p,cl);
		}
	}
}


/// Adds text escape characters
function add_escape_char(te)
{
	const escape_char = [["\\alpha","α"],["\\beta","β"],["\\gamma","γ"],["\\delta","δ"],["\\epsilon","ε"],["\\zeta","ζ"],["\\eta ","η"],["\\theta","θ"],["\\iota","ι"],["\\kappa","κ"],["\\lambda ","λ"],["\\mu","μ"],["\\nu","ν"],["\\xi","ξ"],["\\pi","π"],["\\rho","ρ"],["\\sigma","σ"],["\\tau","τ"],["\\upsilon","υ"],["\\phi","φ"],["\\chi","χ"],["\\psi","ψ"],["\\omega","ω"],["\\sum","Σ"]];
	let i = 0;

	while(i < te.length){
		if(te.substr(i,1) == "\\"){
			let j = 0; 
			while(j < escape_char.length && te.substr(i,escape_char[j][0].length) != escape_char[j][0]) j++;
			if(j < escape_char.length){
				te = te.substr(0,i)+escape_char[j][1]+te.substr(i+escape_char[j][0].length);
			}
			else i++;
		}
		else i++;
	}
	
	return te;
}


/*
/// Gets the classification from elements in the table
function get_cl_from_table(tab,p,col)
{
	let c = find_string_in(tab.heading,col);
	
	if(c == undefined) return "Could not find the column '"+col+"' in file '"+tab.filename+"'";
	
	if(tab.nrow == 0) return "The file '"+tab.filename+"' contains no rows";
	
	let name = tab.ele[0][c];
	
	let cl = get_cl_from_comp(name,p);
	if(cl == undefined) return "In file, '"+name+"' is not a compartment in species '"+model.species[p].name+"'";
	
	return cl;
}
*/

/// Gets the classification from the compartment name
function get_cl_from_comp(name,p)
{
	let sp = model.species[p];
	for(let cl = 0; cl < sp.ncla; cl++){
		let claa = sp.cla[cl];
		for(let c = 0; c < claa.ncomp; c++){
			if(claa.comp[c].name == name) return cl;
		}
	}
}


/// Determines if a valid colour 
function is_Color(color)
{
  const s = new Option().style;
  s.color = color;
  if(s.color == '') return false;
	return true;
}


/// Automatically sets the positions of compartments
function model_set_auto_positions()
{
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
	
			if(claa.camera.coord == "cartesian"){	
				let ncomp = claa.ncomp;
				let comp = claa.comp;
				
				let ntra = claa.ntra;
				let tra = claa.tra;
				
				if(ncomp == 0){
					inter.imp.line = claa.import_line;
					alert_import("The classification '"+claa.name+"' does not contain any compartments.") 
					return;
				}
		
				if(find(comp,"x","auto") != undefined){
					for(let c = 0; c < ncomp; c++){
						if(comp[c].x != "auto"){
							inter.imp.line = comp[c].import_line;
							alert_import("'x' and 'y' should not be specified (because classification '"+claa.name+"' is set to have positions calculated automatically)."); 
							return;
						}
					}
					
					for(let tr = 0; tr < ntra; tr++){
						let traa = tra[tr];
						
						for(let j = 0; j < traa.midp.length; j++){
							if(traa.midp[j].x != "auto"){
								inter.imp.line = traa.import_line;
								
								if((traa.variety == "Source" && j == 0) || 
								   (traa.variety == "Sink" && j == traa.midp.length-1)){
									alert_import("'x' and 'y' should not be specified (because classification '"+claa.name+"' is set to have positions calculated automatically)."); 
								}
								else{
									alert_import("'mid-x' and 'mid-y' should not be specified (because classification '"+claa.name+"' is set to have positions calculated automatically)."); 
								}
								return;
							}
						}
					}
					
					/// Picks the compartment with the fewest entering it as the first 
					let num_in = [];
					for(let c = 0; c < ncomp; c++) num_in[c] = 0;
					
					for(let tr = 0; tr < ntra; tr++){
						if(tra[tr].variety == "Normal"){
							num_in[tra[tr].f]++; 
						}
					}

					let max = LARGE;
					let c_sel;
					for(let c = 0; c < ncomp; c++){
						if(num_in[c] < max){ max = num_in[c]; c_sel = c;}
					}
					
					let list = [];
					for(let c = 0; c < ncomp; c++){
						if(c != c_sel) list.push(c);
					}
					
					let col = [];
					col[0] = {comp:[]};
					col[0].comp.push(c_sel);
					
					let co = 1;
					while(list.length > 0){
						col[co] = {comp:[]};
						
						for(let i = 0; i < col[co-1].comp.length; i++){
							let c = col[co-1].comp[i];

							for(let tr = 0; tr < ntra; tr++){
								if(tra[tr].i == c){
									let cc = tra[tr].f;
									let j = find_in(list,cc);
									if(j != undefined){
										col[co].comp.push(cc);
										list.splice(j,1);
									}
								}
							}							
						}
						
						
						
						if(col[co].comp.length == 0){
							col[co].comp.push(list[0]);
							list.splice(0,1);
						}
						
						co++;
					}
					
					// Places compartments in classification
					
					let dx = import_scale_factor; if(col.length <= 3) dx *= 2;
					
					for(let c = 0; c < comp.length; c++){
						if(comp[c].w > dx/2) dx = Math.floor(2*comp[c].w+1); 
					}
					
					let dy = import_scale_factor;
					
					for(let co = 0; co < col.length; co++){
						let x = co*dx;
						for(let j = 0; j < col[co].comp.length; j++){
							let c = col[co].comp[j];
							
							let y = (j-(col[co].comp.length-1)/2)*dy;
				
							claa.comp[c].x = x; claa.comp[c].y = y; 
						}
					}
					
					for(let tr = 0; tr < ntra; tr++){
						let traa = tra[tr];
						
						switch(traa.variety){
						case "Normal": // Works out if intermediate points need to be added
							{
								let pi = {x:comp[traa.i].x, y:comp[traa.i].y};
								let pf = {x:comp[traa.f].x, y:comp[traa.f].y};
								
								if(comp_line_colide(comp,pi,pf) == true){
									if(pi.y == pf.y){
										let y = pf.y-3;
										
										while(comp_line_colide(comp,{x:pi.x, y:y},{x:pf.x, y:y}) == true) y -= 0.5;
										
										traa.midp.push({x:pi.x, y:y});
										traa.midp.push({x:pf.x, y:y});
									}
									else{
										if(pi.x == pf.x){
											let x = pf.x+dx/2;
										
											while(comp_line_colide(comp,{x:x, y:pi.y},{x:x, y:pf.y}) == true) x += 0.5;
										
											traa.midp.push({x:x, y:pi.y});
											traa.midp.push({x:x, y:pf.y});
										}
									}
								}
							}
							break;
						}
					}						
							
					for(let tr = 0; tr < ntra; tr++){
						let traa = tra[tr];
						
						switch(traa.variety){
						case "Source": case "Sink":
							{
								let i = 0; 
								let co = claa.comp[traa.f];
								if(traa.variety == "Sink"){
									i = traa.midp.length-1;
									co = claa.comp[traa.i];
								}							
								
								let x = co.x-dx, y = co.y; 
								
								if(check_pos_exist(claa,x,y) == false){
									traa.midp[i].x = x; traa.midp[i].y = y;
								}
								else{
									x = co.x; y = co.y+dy; 
									if(check_pos_exist(claa,x,y) == false){
										traa.midp[i].x = x; traa.midp[i].y = y;
									}
									else{
										x = co.x; y = co.y-dy; 
										if(check_pos_exist(claa,x,y) == false){
											traa.midp[i].x = x; traa.midp[i].y = y;
										}
										else{
											x = co.x+dx; y = co.y; 
											if(check_pos_exist(claa,x,y) == false){
												traa.midp[i].x = x; traa.midp[i].y = y;
											}
											else{
												traa.midp[i].x = co.x; traa.midp[i].y = co.y+dy/2;
											}
										}
									}
								}
							}
							break;
						}
					}
				}
			}
		} 
	}
}


/// Determines if a line collides with a compartment
function comp_line_colide(comp,pi,pf)
{
	for(let c = 0; c < comp.length; c++){
		let co = comp[c];
		let d = (co.x-pi.x)*(co.x-pi.x) + (co.y-pi.y)*(co.y-pi.y);
		if(d > TINY){
			let dx = co.x - pi.x, dy = co.y - pi.y;
			let ax = pf.x-pi.x, ay = pf.y-pi.y;
			
			let dist = (ax*dy-ay*dx)/Math.sqrt(ax*ax+ay*ay);
			let alpha = (dx*ax+dy*ay)/(ax*ax+ay*ay);
		
			if(alpha > 0.01 && alpha < 0.99 && dist > -co.h/2 && dist < co.h/2){
				return true;
			}
		}
	}
	
	return false;
}


/// Checks if there is already a compartment at a location
function check_pos_exist(claa,x,y)
{
	for(let c = 0; c < claa.ncomp; c++){
		let co = claa.comp[c];
		if((co.x-x)*(co.x-x) + (co.y-y)*(co.y-y) < TINY) return true;
	}
	
	for(let tr = 0; tr < claa.ntra; tr++){
		let tra = claa.tra[tr];
		
		for(let k = 0; k < tra.midp.length; k++){
			let p = tra.midp[k];
			if((p.x-x)*(p.x-x) + (p.y-y)*(p.y-y) < TINY) return true;
		}
	}
	
	return false;
}


/// Checks that the imported model is correct
function check_import_correct()
{
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		
		// Checks that branching probability has not been specified in error 
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			for(let i = 0; i < claa.ntra; i++){
				let tr = claa.tra[i];
				if(tr.branch == false && tr.bp_set == true){
					inter.imp.line = tr.import_line;
					alert_import("'bp' should not be set for this transition.") 
					return;
				}
			}
		}
		
		/// Checks that non-Markovian transitions not used for population models
		if(sp.type == "population"){
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				for(let i = 0; i < claa.ntra; i++){
					let tr = claa.tra[i];
					if(tr.type != "exp(rate)" && tr.type != "exp(mean)" && tr.type != "erlang"){
						inter.imp.line = tr.import_line;
						alert_import("Because species '"+sp.name+"' has a population-based model this transition must have a type 'exp(rate)', 'exp(mean)' or 'erlang'.");
						return;
					}						
				}
			}
		}
	}
}


/// Checks if a string is a number
function is_number(num,tag)
{
	if(isNaN(num)){
		alert_import("'"+tag+"' must be a number"); 
		return false;
	}
	return true;
}


/// Checks if a string is a positive number
function is_positive(num,tag)
{
	let val = Number(num);
	if(isNaN(num) || val <= 0){
		alert_import("'"+tag+"' must be a positive number"); 
		return false;
	}
	return true;
}


// Checks if a string is between zero and one
function is_zeroone(num,tag)
{
	let val = Number(num);
	if(isNaN(num) || val < 0 || val > 1){
		alert_import("'"+tag+"' must be between 0 and 1"); 
		return false;
	}
	return true;
}


/// Write an array as a string
function stringify(arr)
{
	let te="";
	for(let i = 0; i < arr.length; i++){
		if(i != 0) te += ","; 
		te += arr[i];
	}
	return te;
}

"use strict";
// Functions for imports information from a .bici file

/* // List of commands in alphabetical order
	add-ind
	add-ind-post-sim
	add-ind-sim
	add-pop
	add-pop-post-sim
	add-pop-sim
	box
	camera
	classification / class
	compartment / comp
	compartment-all / comp-all
	comp-data
	data-dir
	derived / der
	description /desc
	fixed-effect
	genetic-data
	ind-effect
	ind-effect-data
	ind-group-data
	inference / inf
	inf-diagnostics
	inf-generation
	inf-param
	inf-state
	inf-diagnostics
	init-pop
	init-pop-sim
	label
	map
	move-ind
	move-ind-post-sim
	move-ind-sim
	param-mult
	parameter
	pop-data
	pop-trans-data
	posterior-simulation/post-sim
	post-sim-param
	post-sim-state
	remove-ind
	remove-ind-post-sim
	remove-ind-sim
	remove-pop
	remove-pop-post-sim
	remove-pop-sim
	species
	set
	simulation / sim
	sim-param
	sim-state
	test-data
	trans-data
	transition / trans
	transition-all / trans-all
	*/
	
let imp = {};                                      // Stores information as import is done

/// Import a script to define (or partially define) model and data 
function import_file(te,clear_results)                                
{		
	import_te = te;
	
	percent(5);
			
	te = add_escape_char(te);
	
	let	lines = te.split('\n');
	
	imp.script = [];
	for(let j = 0; j < lines.length; j++) imp.script.push({line:j, te:lines[j]});
	
	let pro = process_lines(lines);
	
	load_data_files(pro);

	// Keeps track of the current species and classification 
	imp = { lines:lines, script:pro.formatted, previous_loaded_table:[], warn:false}; 

	model = new Model();	 
	model.start_new();
	
	init_result(pro,clear_results);

	percent(10);
		
	for(let loop = 0; loop < 4; loop++){ 
		// Import happens in four stages:
		// (0) Load species and classification information
		// (1) Load parameter information
		// (2) Everything else apart from output samples
		// (3) Output samples
			
		let p, cl;
		if(model.species.length == 1){
			p = 0;
			if(model.species[p].ncla == 1) cl = 0;
		}  
		imp.p = p; imp.cl = cl;
		
		if(loop == 3){
			imp.line = undefined;
			
			percent(70);

			model.determine_branching("set all_branches");

			// Checks information about all parameters have been loaded correctly
			check_param_complete(); 

			let warn = check_initial_pop_error(true);
			if(warn != undefined) alert_import(warn);

			percent(75);
		
			check_import_correct();
			model_set_auto_positions();
			setup_camera();
	
			hash_compgl = undefined;
			hash_tragl = undefined;	

			percent(80);
			
			model.update_pline_all();
		
			if(sim_result.load) results_add_model(sim_result,model.sim_details,"sim");
			if(inf_result.load) results_add_model(inf_result,model.inf_details,"inf");
			if(ppc_result.load) results_add_model(ppc_result,model.ppc_details,"ppc");
			
			percent(85);
		}
		
		for(let m = 0; m < pro.processed.length; m++){
			let line = pro.processed[m];
		
			imp.line = line.line;
			imp.typest = line.type;
		
			let cname = line.type;
			//pr(loop+" "+cname+" loop command");
			switch(cname){                               // Accounts for shortened versions of commands
			case "comp": cname = "compartment"; break;
			case "comp-all": cname = "compartment-all"; break;
			case "trans": cname = "transition"; break;
			case "trans-all": cname = "transition-all"; break;
			case "class": cname = "classification"; break;
			case "desc": cname = "description"; break;
			case "parameter": cname = "param"; break;
			case "der": cname = "derived"; break;
			case "sim": cname = "simulation"; break;
			case "inf": cname = "inference"; break;
			case "posterior-simulation": cname = "post-sim"; break;
			}

			let process = true;
			switch(loop){
			case 0: // In the first pass create species, classifications, compartments and loads data directory
				switch(cname){
				case "species": case "classification": case "camera": 
				case "compartment": case "compartment-all": 
				case "data-dir": break;
				default: process = false; break;
				}
				break;
				
			case 1: // The next pass does loads the parameters 
				process = false; 
				switch(cname){
				case "species": case "classification": case "camera": case "param": case "derived": 
					process = true; 
					break;
				}
				break;
				
			case 2: // The last pass does everything except the parameters and compartments
				switch(cname){
				case "compartment": case "compartment-all":
				case "param": case "derived": case "data-dir":
				case "sim-param": case "sim-state":
				case "inf-param": case "inf-state": 
				case "post-sim-param": case "post-sim-state": 
				case "inf-generation":
				case "inf-diagnostics":
					process = false; break;
				}
				break
				
			case 3:
				process = false; 
				switch(cname){
				case "sim-param": case "sim-state":
				case "inf-param": case "inf-state": 
				case "post-sim-param": case "post-sim-state": 
				case "inf-generation":
				case "inf-diagnostics":
					process = true; 
					break;
				}
				break;
			}
			
			if(process == true){
				process_command(cname,line.tags,loop);
				
				for(let n = 0; n < line.tags.length; n++){
					if(line.tags[n].done != 1){ 
						alert_import("Tag '"+line.tags[n].name+"' not used");
					}
				} 
			}
			
			if(loop < 3){
				percent(10+80*(loop*0.25 + 0.25*(m/pro.processed.length)));
			}
			else{
				percent(80 + 10*(m/pro.processed.length));
			}
			
			if(imp.warn == true) return;
		}
	}
	
	if(sim_result.load) results_finalise(sim_result);
	if(inf_result.load) results_finalise(inf_result);
	if(ppc_result.load) results_finalise(ppc_result);
	
	percent(90);
	
	load_default_map();
	
	percent(95);
	
	//ppc_add_ind();
	
	percent(99);
	
	initialise_filters();

	let in2; if(input.info) in2 = input.info.type;

	if(input.type == "Load Example") model.example = input.info;

	let ans = { model:strip_heavy(model), info:in2, map_store:map_store};
	if(sim_result.load) ans.sim_load = true;
	if(inf_result.load) ans.inf_load = true;
	if(ppc_result.load) ans.ppc_load = true;
	
	import_post_mess(ans);
	
	map_store = [];
}


/// Initialises results based on what information is loaded in file
function init_result(pro,clear_results)
{
	if(clear_results){
		sim_result = {siminf:"sim"};                  // Stores results from simulation
		inf_result = {siminf:"inf"};                  // Stores results from inference
		ppc_result = {siminf:"ppc"};
	}
	
	sim_result.load = false;
	inf_result.load = false;
	ppc_result.load = false;

	for(let m = 0; m < pro.processed.length; m++){
		let line = pro.processed[m];
		let cname = line.type;
		
		switch(cname){
		case "sim-param": case "sim-state": turn_result_on(sim_result,"sim"); break;
		case "inf-param": case "inf-state": turn_result_on(inf_result,"inf"); break;
		case "post-sim-param": case "post-sim-state": turn_result_on(ppc_result,"ppc"); break;
		default: break;
		}
	}
}
	
	
/// Turns a result on
function turn_result_on(result,siminf)
{
	if(result.load == false){
		for(let ele in result){
			delete result[ele];
		}
	
		result.siminf = siminf;
		result.load = true;
		result.on = true;
	}		
}
	
	
/// Adds add_ind, move_ind and remove_ind results to PPC
function ppc_add_ind()
{
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let j = 0; j < sp.inf_source.length; j++){
			let so = sp.inf_source[j];
			switch(so.type){
			case "Add Ind.": case "Remove Ind.": case "Move Ind.":
				{
					let so_new = copy(so);
					so_new.check_info.i = sp.ppc_source.length;
					sp.ppc_source.push(so_new);
				}
				break;
			}
		}
	}
}


/// Loads up default map (if neccesary);
function load_default_map()
{
	let map;
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = model.species[p].cla[cl];			
			if(claa.default_map == "loading"){
				if(map == undefined){			
					let mfi = "../Maps/World.json";
					if(load_map_fast) mfi = "D:/BICI_nolongerused/Maps/World.json";
				
					let te = load_file_http(mfi,"import");
					if(te == "") alert_import("The file '"+filename(mfi)+"' does not exist or is empty.");
						 
					map = load_annotation_map(te);
				}
			
				let name = "file"+Math.random();
				map_store.push({name:name, feature:map.feature, box:map.box});
				
				claa.annotation.push({type:"map", map_ref:name, default_map:true});
	
				claa.default_map = true;	
			}
		}
	}
}
	

/// Processes a given command
function process_command(cname,tags,loop)
{
	imp.tags = tags;
	imp.cname = cname;
	imp.all_row = "";

	switch(cname){
	case "species": species_command(loop); break;
	case "classification": classification_command(loop); break;
	case "set": set_command(); break;
	case "camera": camera_command(); break;
	case "compartment": compartment_command(); break;
	case "compartment-all": compartment_all_command(); break;
	case "transition": transition_command(); break;
	case "transition-all": transition_all_command(); break;
	case "data-dir": datadir_command(); break;
	case "description": description_command(); break;
	case "label": label_command(); break;
	case "box": box_command(); break;
	case "map": map_command(); break;
	case "param": param_command(); break;
	case "param-mult": param_mult_command(); break;
	case "derived": derived_command(); break;
	case "simulation": simulation_command(); break;
	case "inference": inference_command(); break;
	case "post-sim": post_sim_command(); break;
	case "ind-effect": ind_effect_command(); break;
	case "fixed-effect": fixed_effect_command(); break;
	case "sim-param": sim_param_command(); break;
	case "sim-state": sim_state_command(); break;
	case "inf-param": inf_param_command(); break;
	case "inf-state": inf_state_command(); break;
	case "post-sim-param": post_sim_param_command(); break;
	case "post-sim-state": post_sim_state_command(); break;
	case "inf-diagnostics": inf_diagnostics_command(); break;
	case "inf-generation": inf_generation_command(); break;
	
	default: 
		if(find_in(data_command_list,cname) != undefined){
			import_data_table_command(cname);
		}
		else{
			error("'"+cname+"' not recognised");
			alert_import("Command '"+cname+"' not recognised"); 
		}
		break;
	}
	
	let warn = check_initial_pop_error(false);
	if(warn != undefined) alert_import(warn);
}


/// Imports a geojson file
function import_geojson(file)  
{
	let te = file.te;
	
	let da = JSON.parse(te);
	
	let feature = get_feature(da);
	if(feature.length == 0){
		return "There are no features in the file '"+file+"'";
	}	
		
	let tab = get_feature_table(da,feature);
	
	return tab;
}


/// Loads up a table in table format
function import_table(te)
{	
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
			}
		}
	}
	else{ 
		alert_import("Table does not contain any data");
	}

	let tab = { heading:heading, ele:ele, nrow:nrow, ncol:ncol};
	
	return tab;	
}


/// Checks if lat and lng are correctly specified
function check_latlng(lat,lng)
{
	if(lng == "" || lat == ""){
		alert_import("The 'lat' and 'lng' must be set");
	}
	if(isNaN(lat)){
		alert_import("'lat' must be a number");
	}

	if(isNaN(lng)){
		alert_import("'lng' must be a number");
	}

	if(lat > 90 || lat < -90){
		alert_import("'lat' must be in the range -90 to 90");
	}
	
	if(lng > 180 || lng < -180){
		alert_import("'lng' must be in the range -180 to 180");
	}
}


/// Gets all the command tags from a given line
function get_command_tags(trr)
{
	let num_quote = 0;
	let quote_pos = [];
	
	for(let i = 0; i < trr.length; i++){
		if(trr.substr(i,1) == '"'){ num_quote++; quote_pos.push(i);}
	}
	if(num_quote%2 != 0) alert_import("Syntax error: Quotes do not match up."); 
	

	for(let i = 0; i < num_quote; i+= 2){
		if(quote_pos[i]+1 == quote_pos[i+1]){
			alert_import("Syntax error: No content within the quotation marks."); 
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
	
	if(frag.length == 0) alert_import("Does not contain any content"); 
	
	if(frag[0].quote == 1) alert_import("Should not start with quotes");
	let type = frag[0].text, type_pos = frag[0].pos; 
	
	if(find_in(command_list,type) == undefined){ 
		alert_import("Command '"+type+"' not recognised.");
	}
	
	let num = (frag.length-1)/3;

	if(num != Math.floor(num)){ alert_import("Syntax error"); return;}
	let tags=[];

	for(let n = 0; n < num; n++){
		if(frag[1+3*n+0].quote != 0) alert_import("Syntax error type 2"); 
		if(frag[1+3*n+1].text != "=") alert_import("Syntax error type 3");
		if(frag[1+3*n+0].text == "") alert_import("Syntax error type 4");
		if(frag[1+3*n+2].text == "") alert_import("Property "+frag[1+3*n+0].text+" must have content");
		let fr = frag[1+3*n+0];
		tags.push({name:fr.text, pos:fr.pos, pos_end:fr.pos_end, value:frag[1+3*n+2].text, done:0});
	}
	
	for(let n = 0; n < tags.length-1; n++){
		for(let nn = n+1; nn < tags.length; nn++){
			if(tags[n].name == tags[nn].name) alert_import("The tag '"+tags[n].name+"' is set more than once");
		}
	}
	
	return {type:type, type_pos:type_pos, tags:tags};
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
	let tags = imp.tags; 
	imp.tagst = st; 
	for(let i = 0; i < tags.length; i++){
		if(tags[i].name == st){ tags[i].done = 1; return tags[i].value;}
	}
	return "";
}


/// Gets the value of a tag from a specified list
function get_tag_val(st,tags)                                   
{
	imp.tagst = st;
	for(let i = 0; i < tags.length; i++){
		if(tags[i].name == st){ tags[i].done = 1; return tags[i].value;}
	}
	return "";
}


/// Error massage if a tag cannot be found
function cannot_find_tag()                                   
{
	let te = "Cannot find the '"+imp.tagst+"' tag for '"+imp.typest+"'";
	if(imp.all_row != "") te += " "+imp.all_row;
	alert_import(te);
}


/// Error message for imported file (for specific line)
function alert_import(st,line)                                     
{
	if(!line) line = imp.line;
	
	let te;
	if(line) te = "On line "+(line+1)+": "+st;
	else te = st;
	
	if(te.length > 0 && te.substr(te.length-1,1) != ".") te += ".";
	
	get_formatted_line_width();

	throw({ type:"Alert Import", title:"Error importing file!", te:te, st:st, line:line, scroll_to_line:true, script:imp.script});
}


/// Gets the line width 
function get_formatted_line_width()
{
	for(let i = 0; i < imp.script.length; i++){
		let line = imp.script[i];
		let te = line.te;
		te = te.replace(/<c>/g,"");
		te = te.replace(/<t>/g,"");
		te = te.replace(/<p>/g,"");
		te = te.replace(/<g>/g,"");
		te = te.replace(/<\/c>/g,"");
		te = te.replace(/<\/t>/g,"");
		te = te.replace(/<\/p>/g,"");
		te = te.replace(/<\/g>/g,"");
		line.w = 1.05*text_width_para(te);
	}
}


/// Transfers help messages to import messages
function alert_help(title,te)
{		
	throw({type:"AlertHelp", title:title, te:te});
}


/// If the result of a function throws an error then outputs to help
function output_help(res)
{
	if(res.err == true){
		throw({type:"AlertHelp", title:"Sorry an error occured!", te:res.msg});
	}
}

/// Transfers alertp messages to window
function alertp(te)
{
	if(false){ pr("throw"); pr({type:"AlertP", te:te});}
	throw({type:"AlertP", te:te});
}


/// Alerts a specific line in the import file
function alert_line(st,line)                                     
{
	let te = "On line "+(line+1)+": "+st
	if(te.length > 0 && te.substr(te.length-1,1) != ".") te += ".";
	
	throw({type:"Alert Import", title:"Error importing file!", te:te, line:line, scroll_to_line:true, script:imp.script});
}


/// Alerts a general error in the import file
function alert_noline(st)                                     
{
	let te = st;
	if(te.length > 0 && te.substr(te.length-1,1) != ".") te += ".";
	
	throw({ type:"Alert Import", title:"Error importing file!", te:te, scroll_to_line:false, script:imp.script});
}


/// Checks if an error occurs 
function output_error(res)
{
	if(res.err == true){ alert_import(res.msg); return true;}
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
function option_error(na,te,pos,line)
{
	for(let i = 0; i < pos.length; i++) if(pos[i] == te) return false;
	
	let st = "'"+na+"' has a value '"+te+"' but should be chosen from one of the following options: ";
	for(let i = 0; i < pos.length; i++){
		if(i != 0) st += ", ";
		st += "'"+pos[i]+"'";
	}
	
	alert_import(st,line);
	return true;
}


/// Gets the classification object based on store values for p and cl
function get_claa()
{
	let p = imp.p;
	if(p == undefined){
		alert_import("Species needs to be specified before "+imp.cname);
	}
			
	let cl = imp.cl;
	if(cl == undefined){
		alert_import("Classification needs to be specified before "+imp.cname);
	}
	return model.species[p].cla[cl];
}		


/// Processes the lines to extract file information and formats commands
function process_lines(lines)
{
	let formatted=[];
	let processed=[];
	
	let data_dir;
	
	let j = 0;
	while(j < lines.length){
		imp.line = j;
		
		let trr = lines[j].trim();
	
		let flag = false;                              // Ignores line if empty or commented out
		if(trr.length == 0) flag = true;
		if(trr.length >= 2 && trr.substr(0,2) == "//") flag = true;
		if(trr.length >= 1 && trr.substr(0,1) == "#") flag = true;
			
		let command_line; 
		let flag_file = false;
			
		if(flag == true) trr = "<g>"+trr+"</g>";
		else{	
			let files = [];
		
			while(trr.length > 3 && trr.substr(trr.length-3,3) == '"[['){
				j++;
				let file_te = "";
				let file_formatted = "";
				let file_name = "$FILE"+files.length+"$";
				while(j < lines.length){
					let te = lines[j].trim();
		
					if(te.length >= 3 && te.substr(0,3) == ']]"'){
						trr += file_name+te;
						break;
					}
					
					file_formatted += te+"|";
					
					file_te += te+"\n";
					j++;
				}
				if(j == lines.length) alert_import("Cannot find line with ']]\"' to specify end of file");
				
				if(file_formatted.length > 0){
					file_formatted = file_formatted.substr(0,file_formatted.length-1);
				}
				if(file_formatted.length > DISPLAY_FILE_LEN_MAX) file_formatted = file_formatted.substr(0,DISPLAY_FILE_LEN_MAX)+"..."
				
				files.push({name:file_name, te:file_te, formatted:file_formatted});
			}
	
			command_line = get_command_tags(trr);

			let tags = command_line.tags;
			for(let i = 0; i < files.length; i++){
				let k = 0; while(k < tags.length && tags[k].value != "[["+files[i].name+"]]") k++;
				if(k == tags.length) alertp("Could not find file");
				
				tags[k].value = {name:"inline", sep:"comma", te:files[i].te};
			}
		
			if(command_line != undefined){
				for(let m = command_line.tags.length-1; m >= 0; m--){
					let k = command_line.tags[m].pos;
					let k2 = command_line.tags[m].pos_end;
				
					trr = trr.substr(0,k2)+"</t>"+trr.substr(k2);
					trr = trr.substr(0,k)+"<t>"+trr.substr(k);
				}
		
				let k = command_line.type_pos;
				trr = trr.substr(0,k)+"<c>"+command_line.type+"</c>" + trr.substr(k+command_line.type.length);
			}
			
			for(let j = files.length-1; j >= 0; j--){
				let name = files[j].name;
				let i = trr.indexOf(name);
				if(i == -1) alertp("Problem loading");
				trr = trr.substr(0,i-2)+"<p>[["+files[j].formatted+"]]</p>"+trr.substr(i+name.length+2);
			}
		}
		
		if(command_line != undefined){
			if(command_line.type == "data-dir") data_dir = command_line.tags[0].value;
		}
		
		formatted.push({line:imp.line, te:trr});
	
		if(command_line != undefined){
			command_line.line = imp.line;
			processed.push(command_line);
		}		
			
		j++;
	}
	
	return {formatted:formatted, processed:processed, data_dir:data_dir};
}


/// Loads up any data files
function load_data_files(pro)
{
	let previous_loaded = [];

	for(let j = 0; j < pro.processed.length; j++){
		let cl = pro.processed[j];
		
		for(let k = 0; k < cl.tags.length; k++){
			let tag = cl.tags[k];
			
			let file = tag.value;
			if(typeof file == 'string'){
				switch(tag.name){
				case "value": case "boundary": case "constant": case "reparam": case "prior-split": case "dist-split": 
				case "A": case "A-sparse": case "pedigree": case "X": case "file": case "text": case "ind-list":
					{
						if(is_file(file)){
							let k = 0; while(k < previous_loaded.length && file != previous_loaded[k]) k++;
							
							if(k < previous_loaded.length){
								tag.value =  previous_loaded[k];
							}
							else{
								check_char_allowed(file,"<>:\"|?*");

								if(pro.data_dir == undefined) alert_import("The 'data-dir' command must be set");
								let full_name = pro.data_dir+"/"+file;
							
								let te = load_file_http(full_name,"import");
								if(te == "") alert_import("The file '"+filename(file)+"' does not exist or is empty.",cl.line);
						 
								let sep = "comma"; if(file.substr(file.length-4,4) == ".tsv") sep = "tab";
								tag.value = {name:file, full_name:full_name, sep:sep, te:te};
								previous_loaded.push(tag.value);
							}
						}
					}
					break;
				}
			}
		}
	}
}
	

/// Sets the camera if not specified in the input file
function setup_camera()
{
	for(let p = 0; p < model.species.length; p++){ 	
		for(let cl = 0; cl < model.species[p].ncla; cl++){
			let cam = model.species[p].cla[cl].camera;

			if(cam.set == false) set_camera(p,cl);
			else set_ruler(p,cl);
		}
	}
}


/// Adds text escape characters
function add_escape_char(te)
{
	let escape_char = [];
	for(let i = 0; i < greek_latex.length; i++){
		escape_char.push(["\\"+greek_latex[i][0],greek_latex[i][1]]);
	}
	escape_char.push(["\\sum","Σ"]);
	
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


/// Gets the classification from the compartment name (or if formated using S|E or S:0.5|E:1) 
function get_cl_from_comp(name,p)
{
	let sp = model.species[p];
	let spl = split_with_bracket(name,"|");
	
	if(spl.length == 1){	
		let spl2 = name.split(":");
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			for(let c = 0; c < claa.ncomp; c++){
				if(claa.comp[c].name == spl2[0]) return cl;
			}
		}
	}
	else{
		let cl_st;
		for(let i = 0; i < spl.length; i++){
			let cl = get_cl_from_comp(spl[i],p);
			if(cl_st == undefined) cl_st = cl;
			else{
				if(cl_st != cl) return;
			}
		}
		return cl_st;
	}
}


/// Gets the classification from the transition name (or if formated using S|E or S:0.5|E:1) 
function get_cl_from_trans(name,p)
{
	let sp = model.species[p];
	
	let spl = split_with_bracket(name,"|");
	
	if(spl.length == 1){	
		let spl2 = name.split(":");
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			for(let c = 0; c < claa.ntra; c++){
				if(claa.tra[c].name == spl2[0]) return cl;
			}
		}
	}
	else{
		let cl_st;
		for(let i = 0; i < spl.length; i++){
			let cl = get_cl_from_trans(spl[i],p);
			if(cl_st == undefined) cl_st = cl;
			else{
				if(cl_st != cl) return;
			}
		}
		return cl_st;
	}
}




/// Determines if a valid colour 
function is_Color(color)//zz
{
	let allow = "0123456789abcdefABCDEF"
	color = color.trim();
	if(color.substr(0,1) == "#"){
		if(color.length != 7) return false;
		for(let i = 1; i < 7; i++){
			let ch = color.substr(i,1);
			if(!allow.includes(ch)) return false;
		}
	}
	else{
		if(color.substr(0,4) != "rgb(") return false;
		if(color.substr(color.length-1,1) != ")") return false;
		
		let cont = color.substr(4,color.length-5);
		let spl = cont.split(",");
		if(spl.length != 3) return false;
		for(let i = 0; i < 3; i++){
			if(isNaN(spl[i])) return false;
			let num = Number(spl[i]);
			if(num != Math.floor(num)) return false;
			if(num < 0 || num > 255) return false;
		}
	}
	
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
				
				if(find(comp,"x","auto") != undefined){
					for(let c = 0; c < ncomp; c++){
						if(comp[c].x != "auto"){
							imp.line = comp[c].import_line;
							alert_import("'x' and 'y' should not be specified (because classification '"+claa.name+"' is set to have positions calculated automatically)."); 
						}
					}
					
					for(let tr = 0; tr < ntra; tr++){
						let traa = tra[tr];
						
						for(let j = 0; j < traa.midp.length; j++){
							if(traa.midp[j].x != "auto"){
								imp.line = traa.import_line;
								
								if((traa.variety == SOURCE && j == 0) || 
								   (traa.variety == SINK && j == traa.midp.length-1)){
									alert_import("'x' and 'y' should not be specified (because classification '"+claa.name+"' is set to have positions calculated automatically)."); 
								}
								else{
									alert_import("'mid-x' and 'mid-y' should not be specified (because classification '"+claa.name+"' is set to have positions calculated automatically)."); 
								}
							}
						}
					}
					
					// Picks the compartment with the fewest entering it as the first 
					let num_in = [];
					for(let c = 0; c < ncomp; c++) num_in[c] = 0;
					
					for(let tr = 0; tr < ntra; tr++){
						if(tra[tr].variety == NORMAL){
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
						case NORMAL:                         // Works out if intermediate points need to be added
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
						case SOURCE: case SINK:
							{
								let i = 0; 
								let co = claa.comp[traa.f];
								if(traa.variety == SINK){
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
				
				// If sources or sinks set to auto then set position
				for(let tr = 0; tr < ntra; tr++){
					let traa = tra[tr];
					if(traa.midp.length == 1 && traa.midp[0].x == "auto"){
						switch(traa.variety){
							case SINK:
								{
									let co = claa.comp[traa.i];  
									traa.midp[0].x = co.x; traa.midp[0].y = co.y+4; 
								}
								break;
								
							case SOURCE:
								{
									let co = claa.comp[traa.f];  
									traa.midp[0].x = co.x; traa.midp[0].y = co.y-4; 
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
					imp.line = tr.import_line;
					alert_import("'bp' should not be set for this transition.") 
				}
			}
		}
		
		// Checks that non-Markovian transitions not used for population models
		if(sp.type == "population"){
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				for(let i = 0; i < claa.ntra; i++){
					let tr = claa.tra[i];
					//if(tr.type != "exp(rate)" && tr.type != "exp(mean)" && tr.type != "erlang"){
					if(tr.type != "exponential" && tr.type != "erlang"){
						imp.line = tr.import_line;
						alert_import("Because species '"+sp.name+"' has a population-based model this transition must have a type 'exponential' or 'erlang'.");
					}						
				}
			}
		}
	}
}


/// Checks if a string is a number
function is_number(num,tag)
{
	if(isNaN(num)) alert_import("'"+tag+"' must be a number"); 
}


/// Checks if a string is a positive number
function is_positive(num,tag)
{
	let val = Number(num);
	if(isNaN(num) || val <= 0){
		alert_import("'"+tag+"' must be a positive number"); 
	}
}


/// Checks if a string is between zero and one
function is_zeroone(num,tag)
{
	let val = Number(num);
	if(isNaN(num) || val < 0 || val > 1){
		alert_import("'"+tag+"' must be between 0 and 1"); 
	}
}


/// Posts a message (but only one);
function import_post_mess(mess)
{
	if(imp.warn == false) post(mess);
	imp.warn = true;
}

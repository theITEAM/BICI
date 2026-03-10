"use strict";
// Functions for imports information from a .bici file

/* // List of commands in alphabetical order
	add-ind-inf
	add-ind-post-sim
	add-ind-sim
	add-pop-inf
	add-pop-post-sim
	add-pop-sim
	box
	classification / class
	compartment / comp
	compartment-all / comp-all
	comp-data
	data-dir
	define
	derived / der
	description /desc
	diagnostics-inf
	fixed-effect
	generation-inf
	genetic-data
	ind-effect
	ind-effect-data
	ind-group-data
	inference / inf
	init-pop-inf
	init-pop-sim
	label
	map
	move-ind-inf
	move-ind-post-sim
	move-ind-sim
	param / parameter
	param-inf
	param-mult
	param-sim
	param-stats-inf
	pop-data
	pop-trans-data
	posterior-simulation/post-sim
	param-post-sim
	proposal-inf
	remove-ind-inf
	remove-ind-post-sim
	remove-ind-sim
	remove-pop-inf
	remove-pop-post-sim
	remove-pop-sim
	set
	simulation / sim
	species
	state-inf
	state-post-sim
	state-sim
	test-data
	trans-data
	trans-diag-inf
	transition / trans
	transition-all / trans-all
	view
	warning-inf
	warning-post-sim
	warning-sim
	*/
	
let imp = {};                                      // Stores information as import is done

/// Import a script to define (or partially define) model and data 
function import_file(te,file,clear_results)                                
{	
	loading_mess("Processing...");

	percent(0);
	
	import_te = te;

	//te = te.replace(/\r/g,"");
	te = remove_escape_char(te);

	percent(1);
	
	let	lines = te.split('\n');
	
	imp.script = [];
	let jmax = lines.length;
	let fl = false;
	if(jmax > SCRIPT_LINE_MAX){ fl = true; jmax = SCRIPT_LINE_MAX;}
	for(let j = 0; j < jmax; j++) imp.script.push({line:j, te:lines[j]});
	if(fl) imp.script.push({line:jmax, te:(lines.length-jmax)+" more lines (too long to show)"});
		
	let pro = process_lines(lines,file,2,15);
	
	percent(10);

	let data_file_list = get_data_file_list(pro,15,20);

	if(model != undefined){
		model_store = copy(model);
	}
	
	init_result(pro,clear_results);

	// Keeps track of the current species and classification 
	imp = { pro:pro, lines:lines, clear_results:clear_results, script:pro.formatted, previous_loaded_table:[], warn:false}; 

	if(data_file_list.length > 0){
		if(begin_str(data_file_list[0].full_name,"..")) load_local(data_file_list,0,100);
		else post({type:"Import model files", data_file_list:data_file_list});
	}
	else{
		import_file2(data_file_list);
	}
}


/// Loads up files locally
function load_local(data_file_list,per_start,per_end)
{
	loading_mess("Loading...");

	for(let i = 0; i < data_file_list.length; i++){
		percent(per_start+((i+0.5)/data_file_list.length)*(per_end-per_start));
		let dfl = data_file_list[i];
		dfl.te = load_file_local(dfl.full_name);
	}
	
	import_file2(data_file_list);
}


/// Continues importing after all files have been loaded
function import_file2(data_file_list)
{

	let ts = clock();
	loading_mess("Processing...");
	
	let pro = imp.pro;
	imp.data_file_list = data_file_list;

	model = new Model();	 
	model.start_new();
	
	let total_pt = assign_processing_time();    // Estimates how long each command will take
	
	for(let loop = 0; loop <= 4; loop++){ 
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
			
			model.determine_branching("set all_branches");

			// Checks information about all parameters have been loaded correctly			
			check_param_complete(); 

			add_non_proc_time(0.1);
			
			let warn = check_initial_pop_error(true);
			if(warn != undefined) alert_import(warn);

			add_non_proc_time(0.1);
			
			check_import_correct();
			model_set_auto_positions();
			setup_camera();
	
			hash_compgl = undefined;
			hash_tragl = undefined;	
			
			model.update_pline_all();

			add_non_proc_time(0.1);
			
			if(sim_result_import.load) results_add_model(sim_result_import,model.sim_details,"sim");
			if(inf_result_import.load) results_add_model(inf_result_import,model.inf_details,"inf");
			if(ppc_result_import.load) results_add_model(ppc_result_import,model.ppc_details,"ppc");
		
			if(inf_result_import.load) inf_result_import.diagnostics_on = true;
			
			add_non_proc_time(0.5);
			
			initialise_filters_setup();
			
			add_non_proc_time(0.2);
		}
		
		for(let m = 0; m < pro.processed.length; m++){
			let line = pro.processed[m];
		
			imp.line = line.line;
			imp.typest = line.type;
		
			let cname = line.type;
			//prr(loop+" "+cname+" loop command");
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
				case "species": case "classification": case "view":
				case "compartment": case "compartment-all": 
				case "data-dir": 
				case "ind-effect":
				case "define":
					break;
				default:
					process = false; 
					break;
				}
				break;
				
			case 1: // The next pass does loads the parameters 
				process = false; 
				switch(cname){
				case "species": case "classification": case "view": 
				case "param":
					process = true; 
					break;
				}
				break;
				
			case 2: // The last pass does everything except the parameters and compartments
				switch(cname){
				case "compartment": case "compartment-all":
				case "param": case "data-dir":  case "define":
				case "ind-effect":
				case "param-sim": case "state-sim":
				case "param-inf": case "state-inf": 
				case "param-post-sim": case "state-post-sim": 
				case "generation-inf":
				case "diagnostics-inf":
				case "trans-diag-inf": 
				case "add-ind-post-sim": case "remove-ind-post-sim": case "move-ind-post-sim":
				case "add-pop-post-sim": case "remove-pop-post-sim":
					process = false; 
					break;
				}
				break

			case 3:
				process = false; 
				switch(cname){
				case "state-sim": case "state-inf": case "state-post-sim": 
					process = true; 
					break;
				}
				break;
				
			case 4:
				process = false; 
				switch(cname){
				case "param-sim": 
				case "param-inf":
				case "param-post-sim": 
				case "generation-inf":
				case "diagnostics-inf":
				case "trans-diag-inf":
				case "add-ind-post-sim": case "remove-ind-post-sim": case "move-ind-post-sim":
				case "add-pop-post-sim": case "remove-pop-post-sim":
					process = true; 
					break;
				}
				break;
			}
			
			if(cname == "param-stats-inf") process = false;
			
			if(process == true){	
				process_command(cname,loop,line);
				
				for(let n = 0; n < line.tags.length; n++){
					if(line.tags[n].done != 1){ 
						alert_import("Tag '"+line.tags[n].name+"' not used");
					}
				} 
			}
		
			if(imp.warn == true) return;
		}
	}
		
	if(dif(imp.pro.proc_time.pt_sum,1)){
		prr("PTSUM NOT MATCH"+imp.pro.proc_time.pt_sum);
	}
	
	if(sim_result_import.load) results_finalise(sim_result_import);
	if(inf_result_import.load) results_finalise(inf_result_import);
	if(ppc_result_import.load) results_finalise(ppc_result_import);

	percent(90);
	
	load_default_map();
	
	percent(95);
	
	//ppc_add_ind();
	
	percent(99);
	
	initialise_filters();

	let in2; if(input.info) in2 = input.info.type;

	if(input.type == "Load Example") model.example = input.info;

	let ans = { model:strip_heavy(model), info:in2, map_store:map_store};
	
	if(!imp.warn){
		if(sim_result_import.load){ sim_result = sim_result_import; ans.sim_load = true;}
		if(inf_result_import.load){ inf_result = inf_result_import; ans.inf_load = true;}
		if(ppc_result_import.load){ ppc_result = ppc_result_import; ans.ppc_load = true;}
	}
	
	import_post_mess(ans);
	
	sim_result_import = {siminf:"sim"};  
	inf_result_import = {siminf:"inf"};
	ppc_result_import = {siminf:"ppc"};
	
	map_store = [];
	
	imp = {};  // Removes to save memory	
	
	//profiling();  // Looks at where memory is being used
	
	//prr("Processing time: "+(clock()-ts)/1000);
}


/// Estimates how long each command will take
function assign_processing_time()
{
	let pro = imp.pro;
	
	let total_pt = 0;
	/// Estimates time taken to process each command (based on file size) 
	for(let i = 0; i < pro.processed.length; i++){
		let pc = pro.processed[i];
		pc.pt = 1;
		for(let j = 0; j < pc.tags.length; j++){
			let tag = pc.tags[j];
			if(tag.need_pt){
				if(tag.value != undefined && tag.value.te != undefined){
					let fac = 1;
					
					switch(pc.type){
					case "comp-all": fac = 10; break;
					case "param-stats-inf": pc.pt = 0; fac = 0; break;
					}
					
					pc.pt += 0.1*tag.value.te.length*fac;
				}
			}
		}
		total_pt += pc.pt;
	}
	
	for(let i = 0; i < pro.processed.length; i++){
		let pc = pro.processed[i];
		pc.pt = frac_proc*(pc.pt/total_pt);
	}
	
	pro.proc_time = {per_start:20, per_end:90, per:0, pt_sum:0, frac_no_proc:1-frac_proc};
}


/// Adds non_processing 
function add_non_proc_time(fr)
{
	let pro = imp.pro;
	pro.proc_time.pt_sum += fr*pro.proc_time.frac_no_proc;
	show_proc_percent(pro);
}


/// Adds non_processing 
function add_proc_time(num)
{
	let pro = imp.pro;
	pro.proc_time.pt_sum += num;
	show_proc_percent(pro);
}


/// Gets current percentage
function get_percent()
{
	let pro = imp.pro;
	let pt = pro.proc_time;
	return pt.per_start+(pt.per_end-pt.per_start)*(pt.pt_sum);
}


/// Shows percent during the initialisation phase 
function show_proc_percent()
{
	let pro = imp.pro;
	let pt = pro.proc_time;
	let per_new = Math.floor(get_percent());
	if(per_new != pt.per){ percent(per_new); pt.per = per_new;}
}


/// Initialises results based on what information is loaded in file
function init_result(pro,clear_results)
{
	if(clear_results){
		sim_result_import = {siminf:"sim"};                    // Stores results from simulation
		inf_result_import = {siminf:"inf"};                  // Stores results from inference
		ppc_result_import = {siminf:"ppc"};
	}
		
	sim_result_import.load = false;
	inf_result_import.load = false;
	ppc_result_import.load = false;

	for(let m = 0; m < pro.processed.length; m++){
		let line = pro.processed[m];
		let cname = line.type;
		
		switch(cname){
		case "param-sim": case "state-sim": turn_result_on(sim_result_import,"sim"); break;
		case "param-inf": case "state-inf": turn_result_on(inf_result_import,"inf"); break;
		case "param-post-sim": case "state-post-sim": turn_result_on(ppc_result_import,"ppc"); break;
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
		let sp_ppc = model.inf_res.plot_filter.species[p];
		for(let j = 0; j < sp.inf_source.length; j++){
			let so = sp.inf_source[j];
			switch(so.type){
			case "Add Ind.": case "Remove Ind.": case "Move Ind.":
				{
					let so_new = copy(so);
					so_new.check_info.i = sp_ppc.ppc_source.length;
					sp_ppc.ppc_source.push(so_new);
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
				
					let te = load_file_local(mfi,"import");
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
function process_command(cname,loop,line)
{
	let tags = line.tags;
	//prr("PROCESS "+line.pt+" "+cname); 	
	imp.tags = tags;
	imp.cname = cname;
	imp.all_row = "";

	/*
	switch(cname){
	case "param-sim": case "state-sim": case "param-inf": case "state-inf":
	case "param-post-sim": case "state-post-sim": case "diagnostics-inf":
	//case "generation-inf":
	//case "trans-diag-inf": 
	case "proposal-inf": 

		add_proc_time(line.pt);
		line.pt = 0;
		break;
	}
	*/
	
	switch(cname){
	case "species": species_command(loop); break;
	case "classification": classification_command(loop); break;
	case "set": set_command(); break;
	case "view": camera_command(); break;
	case "warning-sim": warning_command("sim"); break;
	case "warning-inf": warning_command("inf"); break;
	case "warning-post-sim": warning_command("ppc"); break;
	case "compartment": compartment_command(); break;
	case "compartment-all": compartment_all_command(line); break;
	case "transition": transition_command(); break;
	case "transition-all": transition_all_command(); break;
	case "data-dir": datadir_command(); break;
	case "description": description_command(); break;
	case "label": label_command(); break;
	case "box": box_command(); break;
	case "map": map_command(); break;
	case "define": define_command(); break;
	case "param": param_command(line); break;
	case "param-mult": param_mult_command(line); break;
	case "derived": derived_command(); break;
	case "simulation": simulation_command(); break;
	case "inference": inference_command(); break;
	case "post-sim": post_sim_command(); break;
	case "ind-effect": ind_effect_command(); break;
	case "fixed-effect": fixed_effect_command(); break;
	case "param-sim": sim_param_command(line); break;
	case "state-sim": sim_state_command(line); break;
	case "param-inf": inf_param_command(line); break;
	case "state-inf": inf_state_command(line); break;
	case "param-post-sim": post_sim_param_command(line); break;
	case "state-post-sim": post_sim_state_command(line); break;
	case "diagnostics-inf": inf_diagnostics_command(); break;
	case "generation-inf": inf_generation_command(line); break;
	case "trans-diag-inf": trans_diag_command(line); break;
	case "proposal-inf": proposal_inf_command(); break;
	
	default: 
		if(find_in(data_command_list,cname) != undefined){
			import_data_table_command(cname);
		}
		else{
			error("'"+cname+"' not recognised");
			alert_import("Command '"+cname+"' not recognised."); 
		}
		break;
	}
	
	let warn = check_initial_pop_error(false);
	if(warn != undefined) alert_import(warn);
	
	if(line.pt != 0){
		add_proc_time(line.pt);
		line.pt = 0;
	}
}


/// Imports a geojson file
function import_geojson(file)  
{
	let te = file.te;
	
	let da = JSON.parse(te);
	
	let feature = get_feature(da);
	if(feature.length == 0){
		return in_file_text(file)+" there are no features";
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
				alert_import("The rows do not all have the same size: "+ele[j].length+" "+ncol);
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
function check_latlng(lat,lat_tag,lng,lng_tag)
{
	if(lng == "" || lat == ""){
		alert_import("'"+lat_tag+"' and '"+lng_tag+"' must both be set");
	}
	
	is_number(lat,lat_tag);
	
	is_number(lng,lng_tag);

	if(lat > 90 || lat < -90){
		alert_import("For '"+lat_tag+"' the value '"+lat+"' must be in the range -90 to 90");
	}
	
	if(lng > 180 || lng < -180){
		alert_import("For '"+lng_tag+"' the value '"+lng+"' must be in the range -180 to 180");
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
	
	switch(type){
	case "sim-param": type = "param-sim"; break;
	case "inf-param": type = "param-inf"; break;
	case "post-sim-param": type = "param-post-sim"; break;
	case "sim-state": type = "state-sim"; break;
	case "inf-state": type = "state-inf"; break;
	case "post-sim-state": type = "state-post-sim"; break;
	case "sim-warning": type = "warning-sim"; break;
	case "inf-warning": type = "warning-inf"; break;
	case "post-sim-warning": type = "warning-post-sim"; break;
	}
	
	let k = find(command_list,"na",type);
	if(k == undefined) alert_import("Command '"+type+"' not recognised.");
		
	let num = (frag.length-1)/3;
	let numi = Math.floor(num);

	for(let n = 0; n < num; n++){
		let ii = 1+n*3;

		if(frag[ii].text == "="){
			alert_import("An equals sign '=' is misplaced.");
			return syntax_error();
		}
	}
	
	for(let n = 0; n < num; n++){
		let ii = 1+n*3;
			
		if(ii+2 >= frag.length){
			if(ii+1 < frag.length && frag[ii+1].text == "="){
				alert_import("The tag '"+frag[ii].text+"' is unset");
			}
			else{
				alert_import("The tag '"+frag[ii].text+"' cannot be understood");
			}
			return syntax_error();
		}
		
		if(frag[ii+1].text != "="){
			alert_import("The tag '"+frag[ii].text+"' is missing an equals sign");
			return syntax_error();
		}
		
		if(ii+2 < frag.length && frag[ii+2].text == "="){
			alert_import("The tag '"+frag[ii].text+"' cannot be followed by '=='");
			return syntax_error();
		}
		
		if(ii+3 < frag.length && frag[ii+3].text == "="){
			alert_import("The tag '"+frag[ii].text+"' is unset");
			return syntax_error();
		}
	}

	if(num != numi){ alert_import("Syntax error"); return;}
	let tags=[];

	for(let n = 0; n < num; n++){
		if(frag[1+3*n+0].quote != 0) alert_import("Syntax error"); 
		if(frag[1+3*n+1].text != "=") alert_import("Syntax error");
		if(frag[1+3*n+0].text == "") alert_import("Syntax error: Tag name not specified");
		if(frag[1+3*n+2].text == "") alert_import("Tag "+frag[1+3*n+0].text+" must have content");
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
	if(i == tab.ncol)	alert_import("The table does not have a column with heading '"+colname+"'");
	
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


/// Gets the process time used for a given tag
function get_tag_pt(st)                                   
{
	let tags = imp.tags; 
	imp.tagst = st; 
	for(let i = 0; i < tags.length; i++){
		if(tags[i].name == st){ tags[i].done = 1; return tags[i].pt;}
	}
	return 0;
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
	alert_import(te);
}


/// Error message for imported file (for specific line)
function alert_import(st,line)                                     
{
	if(!line) line = imp.line;
	
	let te;
	if(line) te = "On line "+(line+1)+": "+st;
	if(imp.all_row != undefined && imp.all_row != "") te += " "+imp.all_row;
	
	else te = st;
	
	te = add_full_stop(te);
	
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
	if(false){ error("throw"); error({type:"AlertP", te:te});}
	throw({type:"AlertP", te:te});
}


/// Alerts a specific line in the import file
function alert_line(st,line)                                     
{
	let te = "On line "+(line+1)+": "+st
	te = add_full_stop(te);
	
	throw({type:"Alert Import", title:"Error importing file!", te:te, line:line, scroll_to_line:true, script:imp.script});
}


/// Alerts a general error in the import file
function alert_noline(st)                                     
{
	let te = add_full_stop(st);
	
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
		if(i != 0){
			if(i+1 == pos.length) st += " or ";
			else st += ", ";
		}
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
		alert_import("A species needs to be specified before the '"+imp.cname+"' command can be added");
	}
			
	let cl = imp.cl;
	if(cl == undefined){
		alert_import("A classification needs to be specified before the '"+imp.cname+"' command can be added");
	}
	
	return model.species[p].cla[cl];
}		


/// Processes the lines to extract file information and formats commands
function process_lines(lines,file,per_start,per_end)
{
	let formatted=[];
	let processed=[];
	
	let data_dir;
	
	let per = per_start;
	
	let j = 0;
	while(j < lines.length){
		let per_new = per_start + Math.floor((per_end-per_start)*j/lines.length);
		if(per_new != per){
			percent(per_new);
			per = per_new;
		}
		
		imp.line = j;
		
		// Removes brackets for individual effects
		let trr = lines[j].trim();
		
		for(let i = 0; i < trr.length; i++){
			let ch = trr.substr(i,1);
			if(ch == "<") trr = trr.substr(0,i)+"〈"+trr.substr(i+1);
			if(ch == ">" && !(i > 0 && trr.substr(i-1,1)== "-")){
				trr = trr.substr(0,i)+"〉"+trr.substr(i+1);
			}			
		}
	
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
			if(command_line.type == "data-dir"){
				data_dir = command_line.tags[0].value;
				if(!data_dir.includes(":") && file != ""){ // This is a relative path					
					let root = find_root(file);
					data_dir = root+data_dir;
					data_dir = data_dir.replace(/\\/g,"/");
				}
			}
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
function get_data_file_list(pro,per_start,per_end)
{
	let load_list=[];
	for(let j = 0; j < pro.processed.length; j++){
		let cl = pro.processed[j];
		for(let k = 0; k < cl.tags.length; k++){		
			let tag = cl.tags[k];
			
			let file = tag.value;
			switch(tag.name){
			case "value": case "constant": case "reparam": case "text":  // May or may not be a file
				tag.need_pt = true;
				if(typeof file == 'string' && is_file(file)){
					load_list.push({j:j, k:k});
				}
				break;
					
			case "file": case "boundary": 
			case "prior-split": case "dist-split": 
			case "A": case "Ainv": case "A-sparse": case "pedigree": case "X": 
			case "ind-list": case"factor-weight":
				{
					tag.need_pt = true;
					if(typeof file == 'string' && is_file(file)){
						load_list.push({j:j, k:k});
					}
					else{
						if(typeof(file) != "object"){
							alert_line("A data table was expected for '"+tag.name+"'.",cl.line);
						}
					}
				}
				break;
				
			default:
				if(typeof(file) == "object" || (typeof file == 'string' && is_file(file))){
					let extra="";
					switch(tag.name){
					case "prior": 
						extra = " The tag 'prior-split' can be used to specify priors separately for each parameter element."; 
						break;
					case "dist": 
						extra = " The tag 'dist-split' can be used to specify distributions separately for each parameter element."; 
						break;
					}
					alert_line("A data table was not expected for '"+tag.name+"'."+extra,cl.line);
				}
				break;
			}
		}	
	}

	let data_file_list = [];
	let per = per_start;
	for(let i = 0; i < load_list.length; i++){
		let per_new = per_start + Math.floor((per_end-per_start)*i/load_list.length);
		if(per_new != per){
			percent(per_new);
			per_new = per;
		}
		
		let j = load_list[i].j, k = load_list[i].k;
	
		let cl = pro.processed[j];
		let tag = cl.tags[k];
			
		let file = tag.value;
		
		let m = 0; while(m < data_file_list.length && file != data_file_list[m].name) m++;
		
		tag.value = {ref:m};
		
		if(m == data_file_list.length){
			check_char_allowed(file,"<>:\"|?*");

			if(pro.data_dir == undefined) alert_import("The 'data-dir' command must be set");
		
			let full_name = pro.data_dir+"/"+file;
			
			let sep = "comma"; if(file.substr(file.length-4,4) == ".tsv") sep = "tab";
			
			data_file_list.push({name:file, full_name:full_name, sep:sep});
		}
	}
	
	percent(per_end);
	
	return data_file_list;
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
function is_Color(color)
{
	let allow = "0123456789abcdefABCDEF";
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
		switch(sp.type){
		case "Population": case "Deterministic":
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				for(let i = 0; i < claa.ntra; i++){
					let tr = claa.tra[i];
					if(tr.type != "exp(rate)" && tr.type != "exp(mean)" && tr.type != "erlang"){
						imp.line = tr.import_line;
						alert_import("Because species '"+sp.name+"' has a population-based model this transition must have a type 'exponential' or 'erlang'.");
					}						
				}
			}
			break;
		}
	}
}


/// Checks if a string is a number
function is_number(num,tag)
{
	if(isNaN(num)) alert_import("For '"+tag+"' the value '"+num+"' must be a number"); 
}


/// Checks if a string is a positive number
function is_positive(num,tag)
{
	let val = Number(num);
	if(isNaN(num) || val <= 0){
		alert_import("For '"+tag+"' the value '"+num+"' must be a positive number"); 
	}
}


/// Checks if a string is a positive integet
function is_pos_int(num,tag)
{
	let val = Number(num);
	if(isNaN(num) || val != Math.floor(val) || val <= 0){
		alert_import("For '"+tag+"' the value '"+num+"' must be a positive integer"); 
	}
}


/// Checks if a string is between zero and one
function is_zeroone(num,tag)
{
	let val = Number(num);
	if(isNaN(num) || val < 0 || val > 1){
		alert_import("For '"+tag+"' the value '"+num+"' must be between 0 and 1"); 
	}
}


/// Posts a message (but only one);
function import_post_mess(mess)
{
	if(imp.warn == false) post(mess);
	imp.warn = true;
}


/// Gets file text from a given tag
function get_fi(file)
{
	let fi;
	if(file.te != undefined) fi = file;
	else{
		if(file.ref != undefined) fi = imp.data_file_list[file.ref];
	}
	
	if(fi == undefined) error("Problem getting file");
	else{
		let enc = get_tag_value("compress");
		if(enc.toLowerCase() == "true") fi.encode = true;
	}
	
	return fi;
}


/// If text from file then extracts
function text_from_file(te)
{
	if(is_file(te)){
		let fi = get_fi(te);
		return fi.te;
	}
	else{
		return te.replace(/\|/g,"\n");
	}
}


/// Loads up parameter values from a file
function load_param_value(par,value,head_col,valu,desc,dpt)
{
	let dim = get_dimensions(value);
	
	let ele_list = get_element_list(value,dim);
	
	add_proc_time(0.2*dpt);
		
	// Sets a default value of zero
	for(let k = 0; k < ele_list.length; k++){
		set_element(value,ele_list[k],"0");
	}
	
	add_proc_time(0.2*dpt);

	let fi = get_fi(valu);
	
	let tab = load_table(fi.te,true,fi.sep,fi.name);
	
	add_proc_time(0.1*dpt);
	
	if(typeof tab == 'string') alert_import(tab);

	let load_type = "normal";
	
	let col_name = copy(par.dep);
	col_name.push(head_col);
	
	let subtab = get_subtable(tab,col_name);
	if(subtab.error != ""){
		if(is_matrix(par)){ // Tries to load a matrix
			col_name = par.list[0];
			subtab = get_subtable(tab,col_name);
			if(subtab.error != ""){ alert_import(subtab.error); return;}
			load_type = "matrix";
		}
		else alert_import(subtab.error);
	}
	
	switch(load_type){
	case "matrix":    // Loads using a matrix format
		{
			let sym = is_symmetric(par);
			let fpt = 0.5*dpt/par.list[0].length;
			for(let j = 0; j < par.list[0].length; j++){
				add_proc_time(fpt);
				for(let i = 0; i < par.list[1].length; i++){
					let ind=[]; ind[0] = j; ind[1] = i;
					let ele = subtab.ele[j][i];
					
					if(sym == true && i < j){
						if(ele != "."){
							let name = "<e>"+par.name+"_"+par.list[0][j]+","+par.list[0][i]+"</e>";
							alert_import(desc+" the element '"+ele+"' for "+name+" should be set to '.' because it is in the bottom left-hand corneter of the covariance matrix definition");
						}
						set_element(value,ind,ele);
					}
					else{
						let val = get_val_from_ele(ele,par,head_col,desc,j,subtab.heading[i]);

						set_element(value,ind,val);
					}
				}
			}
		}
		break;
		
	case "normal":
		{
			let ncol = subtab.ncol;

			let hash_vec = [];
			for(let i = 0; i < ncol-1; i++){ 
				hash_vec[i] = new Hash();
				
				hash_vec[i].create_vec(par.list[i]);
			}
			
			let fpt = 0.5*dpt/subtab.nrow;
			for(let r = 0; r < subtab.nrow; r++){
				add_proc_time(fpt);
				
				let eler = subtab.ele[r];
				
				let ind = [];
				for(let i = 0; i < ncol-1; i++){
					ind[i] = hash_vec[i].find(eler[i]);
					if(ind[i] == undefined){
						alert_import(desc+" the element '"+eler[i]+"' is not valid (column '"+subtab.heading[i]+"', row "+(r+2)+")");
					}
				}
				
				let ele = eler[ncol-1];
				
				let val = get_val_from_ele(ele,par,head_col,desc,r,subtab.heading[ncol-1]);
				
				set_element(value,ind,val);
			}
		}
		break;
		
	default: error("op problem"); break;
	}
			
	if(head_col == "Value" && is_symmetric(par)){  // Checks that matrix is symmtric 
		for(let j = 0; j < par.list[0].length; j++){
			for(let i	= 0; i < j; i++){
				let ind1=[]; ind1[0] = j; ind1[1] = i;
				let el1 = get_element(value,ind1);
				
				if(el1 != "."){
					let name = "<e>"+par.name+"_"+par.list[0][j]+","+par.list[0][i]+"</e>";
					alert_import(desc+" the element '"+el1+"' for "+name+" should be set to '.' because it is in the bottom left-hand corneter of the covariance matrix definition");
				}
			}
		}
	}
}


/// Gets a value from an element
function get_val_from_ele(ele,par,head_col,desc,r,col)
{
	let val;
	
	ele = ele.trim();
	if(ele == ""){
		alert_import(desc+" no content in table element (column '"+col+"', row "+(r+2)+").");
		return ele;
	}

	switch(head_col){
	case "Prior":
		val = convert_text_to_prior(ele,par.pri_pos);
		if(val.err == true){
			alert_import(desc+" the table element '"+ele+"' is not a valid prior specification (column '"+col+"', row "+(r+2)+"). "+val.msg);
		}
		break;
		
	case "Dist":
		val = convert_text_to_prior(ele,par.pri_pos,true);
		if(val.err == true){
			alert_import(desc+" the table element '"+ele+"' is not a valid distribution specification (column '"+col+"', row "+(r+2)+"). "+val.msg);
		}
		break;

	default:
		switch(par.variety){
		case "normal": case "const":
			if(par.factor && ele == "*"){
				val = ele;
			}
			else{
				val = Number(ele);
				if(isNaN(val)){
					alert_import(desc+" the element '"+ele+"' is not a number (column '"+col+"', row "+(r+2)+")");
				}
			}
			break;
		
		case "reparam":
			if(!isNaN(ele)) val = Number(ele);
			else{
				let conv_res = detect_greek(ele,0);
				ele = conv_res.te;
			
				val = ele;
				
				let res = is_eqn(val,"Table element",{});
				if(res.err == true){
					alert_import(desc+" the element '"+ele+"' is not a valid equation (column '"+col+"', row "+(r+2)+")");
				}
			}
			break;
			
		default: error("option not recognised1"); break;
		}
	}
	
	return val;
}


/// Sets up filters 
function initialise_filters_setup()
{
	if(sim_result_import.on == true) initialise_plot_filters_setup(sim_result_import,model.sim_res);
	if(inf_result_import.on == true) initialise_plot_filters_setup(inf_result_import,model.inf_res);
	if(ppc_result_import.on == true) initialise_plot_filters_setup(ppc_result_import,model.ppc_res);
}
	
	
/// Initialises filters for output plots 
function initialise_filters()
{
	if(sim_result_import.on == true) initialise_plot_filters(sim_result_import,model.sim_res);
	if(inf_result_import.on == true) initialise_plot_filters(inf_result_import,model.inf_res);
	if(ppc_result_import.on == true) initialise_plot_filters(ppc_result_import,model.ppc_res);
}


/// Compresses a string using the LZW algorithm
function decode(fi)
{
	if(fi.encode != true) return fi.te;

	let output=[];
	
	let te = fi.te;
	for(let i = 0; i < te.length; i++){ 	
		let c = te.charCodeAt(i);
		if(c >= COMPRESS_NUM_MAX){
			i++;
			let c2 = te.charCodeAt(i);
			if(c2 > 92) c2--;
			c += COMPRESS_NUM_MAX*(c2-35);
		}
		
		if(c != 10){
			if(c > 92) c--;
			output.push(c-35);
		}
	}		
	
	let imax = output.length;

	let dic=[];
	
	for(let i = 0; i < dic_list.length; i++){
		dic.push(dic_list[i]);
	}
	
	te="";

	let prev="";
	for(let i = 0; i < imax; i++){
		let val = output[i];
		let W;
		if(val < dic.length) W = dic[val];
		else W = prev +prev[0];
		
		te += W;
		if(prev != ""){
			let st = prev + W[0];
			
			if(dic.length < COMPRESS_DIC_MAX) dic.push(st);
		}
		prev = W;
	}
	
	// Adds in any repeated lines
	{
		let fl = false;
		let lines = te.split('\n');
		for(let li = 0; li < lines.length; li++){
			if(lines[li].length > 1){
				if(lines[li].substr(0,1) == "@"){
					let num = Number(lines[li].substr(1));
					if(!isNaN(num)){
						lines[li] = lines[li-num];
						fl = true;
					}
				}
			}
		}		
		
		if(fl == true){
			te = "";
			for(let li = 0; li < lines.length; li++){
				te += lines[li];
				if(li+1 < lines.length) te += endl;
			}
		}
	}
	return te;
}


/// Sets the optimise options
function set_optimise_option(details)
{
	let optimise = get_tag_value("optimise").toLowerCase();
	if(optimise == "") optimise = "auto";
	
	if(option_error("optimise",optimise,["speed","memory","auto"]) == true) return;
	
	details.optimise.value = optimise;
}


/// Sets the optimise options
function set_compress_option(details)
{
	let compress = get_tag_value("compress").toLowerCase();
	if(compress == "") compress = "auto";
	
	if(option_error("compress",compress,["always","never","auto"]) == true) return;
	
	details.compress.value = compress;
}

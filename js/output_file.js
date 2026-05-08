"use strict";
// Functions which generate output .bici file

/// Creates an output file
// save_type has the values:
// sim - Saves a simulation (this performs additional simulation checks)
// inf - Saves an inference (this performs additional inference checks)
// save - Saves files (somes things are not checked)

// one_file deterines if a a single file is generated or a file and directory are made

/// Creates the ouytput file from the model
function create_output_file(save_type,one_file,map_store,bscript_add)
{	
	let bscript = create_bscript(save_type,map_store);
	
	if(bscript_add != undefined){
		for(let i = 0; i < bscript_add.length; i++) bscript.push(bscript_add[i]);
	}
	
	let file_list=[];
	let te = generate_text_from_bscript(bscript,file_list,one_file);
	
	if(input.type == "View Code"){
		te = remove_escape_char(te);	
		let	lines = te.split('\n');
		
		let pro = process_lines(lines,"");
		
		post({formatted:pro.formatted, param_factor:strip_heavy(model.param_factor), param:strip_heavy(model.param), species:strip_heavy(model.species)});
	}
	else{
		percent(90);
		write_file_store(te,"bicifile",file_list,"bicifile");
		percent(100);
		post({save_type:save_type, file_list:file_list, param_factor:strip_heavy(model.param_factor), param:strip_heavy(model.param), species:strip_heavy(model.species)});
	}
}


/// Creates a bscript from the model
function create_bscript(save_type,map_store)
{
	percent(0);

	update_model();

	percent(30);
	
	let bscript=[];
	
	create_output_details(save_type,bscript);

	percent(40);

	check_data_valid_all(save_type);
	
	check_data_time(save_type);
	
	check_param_valid(save_type);
	
	check_integral_bounds(save_type);
	
	model.check_det_species(save_type);
	
	percent(45);
	
	for(let p = 0; p < model.species.length; p++){
		create_output_species(p,save_type,map_store,bscript);
	}

	percent(50);

	create_output_param(save_type,bscript);
	
	percent(60);
	
	create_output_derived(bscript);
	
	if(model.warn.length > 0) err_warning();
	
	percent(70);
	
	percent(80);

	return bscript;
}


/// Generates the final text file from the bscript
function generate_text_from_bscript(bscript,file_list,one_file)
{
	let sampledir_rel_sim = "output-sim/";
	let sampledir_rel_inf = "output-inf/";
	let sampledir_rel_post_sim = "output-post-sim/";
	
	let te = "";
	
	let species_list=[];
	let chain_list=[];
	for(let i = 0; i < bscript.length; i++){ // Works out all the species 
		let com = bscript[i];
		if(com.type == "species"){	
			let name = find_tag(com,"name");
			if(find_in(species_list,name) == undefined) species_list.push(name);
		}
		if(tag_exist(com,"chain")){
			let name = find_tag(com,"chain");
			if(find_in(chain_list,name) == undefined) chain_list.push(name);
		}			
	}
	let multi_species = false;
	if(species_list.length > 1) multi_species = true;
	let multi_chain = false;
	if(chain_list.length > 1) multi_chain = true;
	
	let imin = 0;               // Puts process comment at beginning
	if(bscript.length > 0){
		let com = bscript[0];
		if(com.type == "comment" && begin_str(com.te,"# PROCESSED USING BICI")){
			te += com.te+endl+endl;
			imin++;
		}		
	}
	
	if(!one_file){
		te += banner_text("DATA DIRECTORY")+endl;
		te += endl;
		te += 'data-dir folder="."'+endl;
		te += endl;
	}
	
	let species_sel = "";
	let cla_sel = "";
	
	for(let i = imin; i < bscript.length; i++){
		let com = bscript[i];
	
		switch(com.type){  // Keeps track of which which
		case "species": species_sel = find_tag(com,"name"); break;
		case "class": cla_sel = find_tag(com,"name"); break;	
		case "set":
			{
				if(tag_exist(com,"species")) species_sel = find_tag(com,"species");
				if(tag_exist(com,"classification")) cla_sel = find_tag(com,"classification");
			}
			break;
		}
	
		if(com.type == "comment") te += com.te;
		else{
			te += com.type;
			for(let j = 0; j < com.tags.length; j++){
				let ta = com.tags[j];
				
				te += " "+ta.name+"=";
				let val = ta.value;
				if(typeof val == 'object'){
					if(val.name == "inline") val = val.te;
					else error("problem gen bscript");
				}

				let is_file = false;
				
				switch(is_tag_name_file(ta.name)){
				case "yes": 
					is_file = true; 
					break;
					
				case "maybe": 
					if(val.indexOf("\n") != -1) is_file = true; 
					break;
				}
				
				if(is_file){
					if(one_file || begin_str(com.type,"warning")){
						te += '"[['+endl+val+']]"';
					}
					else{
						let fi_name = com.type;
						let end = '.csv';
						switch(com.type){
						case "param-sim": 
							fi_name = sampledir_rel_sim+"param"+chain_num(com,multi_chain);
							break;
							
						case "state-sim": 
							fi_name = sampledir_rel_sim+"state"+chain_num(com,multi_chain);
							end = ".txt";
							break;
							
						case "param-inf": 
							fi_name = sampledir_rel_inf+"param"+chain_num(com,multi_chain);
							break;
							
						case "state-inf": 
							fi_name = sampledir_rel_inf+"state"+chain_num(com,multi_chain);
							end = ".txt";
							break;
							
						case "trans-diag-inf":
							fi_name = sampledir_rel_inf+"diagnostic/trans_diag"+chain_num(com,multi_chain);
							end = ".txt";
							break;
							
						case "proposal-inf":
							fi_name = sampledir_rel_inf+"prop-info"+chain_num(com,multi_chain);
							end = ".txt";
							break;
							
						case "diagnostics-inf":
							fi_name = sampledir_rel_inf+"diagnostic/info"+chain_num(com,multi_chain);
							end = ".txt";
							break;
						
						case "param-stats-inf":
							fi_name = sampledir_rel_inf+"param-stats"+chain_num(com,multi_chain);
							break;
							
						case "generation-inf":
							fi_name = sampledir_rel_inf+"diagnostic/generation";
							break;
						
						case "param-post-sim": 
							fi_name = sampledir_rel_sim+"param"+chain_num(com,multi_chain);
							break;
							
						case "state-post-sim": 
							fi_name = sampledir_rel_sim+"state"+chain_num(com,multi_chain);
							end = ".txt";
							break;
							
						case "map": end = '.geojson'; break;
						
						case "comp":
							if(ta.name == "boundary"){
								let name = find_tag(com,"name");	
								fi_name = name+"_boundary"; end = '.geojson';
							}
							break;
						
						case "comp-all":
							fi_name = "comp_"+cla_sel;
							break;
							
						case "param": 
							{
								let name = find_tag(com,"name");
								let spl = name.split("_");
								fi_name = spl[0]+"_"+ta.name;
							}
							break;
						
						case "ind-effect":
							fi_name = find_tag(com,"name")+"_"+ta.name;
							break;
							
						case "fixed-effect":
							fi_name = find_tag(com,"name")+"_"+ta.name;
							break;
							
						case "description": 
							end = ".txt";
							break;
						}		

						let file = get_unique_file(fi_name,file_list,end)
						write_file_store(val,file,file_list);
						
						te += '"'+file+'"';
					}
				}
				else{
					if(isNaN(val)) te += '"'+val+'"';
					else te += val;
				}
			}
		}
		te += endl+endl;
	}

	return te;
}


/// Gets the chain number if needed
function chain_num(com,multi_chain)
{
	if(multi_chain == true && tag_exist(com,"chain") != undefined){
		return "_"+find_tag(com,"chain");
	}
	return "";
}


/// Determines if a tag exists or not
function tag_exist(com,na)
{
	if(find(com.tags,"name",na) != undefined) return true;
	return false;
}


/// Finds the value of a tag
function find_tag(com,na)
{
	let k = find(com.tags,"name",na);
	if(k == undefined) return "file";
	
	return com.tags[k].value;
}

/// If there are error then display
function err_warning()
{
	if(model.warn.length > 0){
		if(loading_from_file()){
			let wa = model.warn[0];
			let msg = wa.mess+": "+wa.mess2;
			if(wa.line) alert_line(msg,wa.line)
			else alert_import(msg);
			return;
		}
	}
	
	let full_warn = false;
	
	switch(input.type){
	case "View Code": case "Export BICI": case "Start": case "Save BICI": 
	case "StartClusterSave": case "StartClusterExport": case "StartPPC": case "StartEXT":
		full_warn = true;
		break;
	}

	let wa = {type:"Model warning", command_type:input.type, species:strip_heavy(model.species), warn:model.warn, full_warn:full_warn};
	if(try_on == false) prr(wa);
	else throw(wa);
}


/// Creates a banner for the output file
function banner(name,bscript)
{
	bscript.push({type:"comment", tags:[], te:banner_text(name)});
}


/// Creates a banner for the output file
function mini_banner(name,bscript)
{
	bscript.push({type:"comment", tags:[], te:"# "+name});
}


/// Generates details
function create_output_details(save_type,bscript)
{
	banner("DESCRIPTION",bscript);
	add_bscript("description",[{name:"text", value:model.description.te}],bscript);

	banner("DETAILS",bscript);
	create_output_siminf(save_type,bscript);
}
	
	
/// Adds a command to the bscript
function add_bscript(name,tags,bscript)
{
	bscript.push({type:name, tags:tags});
}
	

/// Adds a tag to a vector of tags
function add_tag(name,value,tags)
{
	tags.push({name:name, value:value});
}

	
/// Outputs all the results for a given species
function create_output_species(p,save_type,map_store,bscript)
{
	let sp = model.species[p];
		
	banner("DEFINE MODEL AND DATA FOR SPECIES "+sp.name.toUpperCase(),bscript);
	
	{
		let ty = sp.type.toLowerCase();
		
		let tags=[];
		
		add_tag("name",sp.name,tags);
		add_tag("type",ty,tags);
		
		if(ty == "individual"){
			if(sp.trans_tree.check == true) add_tag("trans-tree","on",tags);
		}
		
		add_bscript("species",tags,bscript);
	}
	
	create_output_compartments(p,map_store,bscript);

	create_output_ind_eff(p,bscript);

	create_output_fix_eff(p,save_type,bscript);

	create_output_sim_inf_source(p,"sim",bscript);
	if(save_type == "sim") output_check("sim");	

	create_output_sim_inf_source(p,"inf",bscript);
	if(save_type == "inf") output_check("inf");	
	
	if(model.inf_res.on == true){
		create_output_sim_inf_source(p,"ppc",bscript);
		if(save_type == "ppc") output_check("ppc");	
	}
}


/// Creates output for all the compartments and transitions in the model
function create_output_compartments(p,map_store,bscript)
{
	mini_banner("SPECIES MODEL",bscript);
		
	let sp = model.species[p];
	
	if(sp.type == "Individual" && sp.trans_tree.check == true){
		if(sp.infection_cl.te == select_str){
			add_warning({mess:"Transition tree problem", mess2:"The infection classification for species '"+sp.name+"' must be set.", p:p, warn_type:"TransTreeInf"});
		}
		else{
			let cl = find(sp.cla,"name",sp.infection_cl.te);
			if(cl == undefined) error("Problem with classification name");
			else{
				let claa = sp.cla[cl];
				
				let c = 0; while(c < claa.ncomp && claa.comp[c].infected.check != true) c++;
				if(c == claa.ncomp){
					add_warning({mess:"Transition tree problem", mess2:"The infection state of at least one compartment in classification '"+claa.name+"' must be set.", p:p, cl:cl, warn_type:"ModelClass"});
				}
			}
		}
	}
	
	for(let cl = 0; cl < sp.ncla; cl++){
		let claa = sp.cla[cl];
		let cam = claa.camera;
		
		{
			let tags=[];
			add_tag("name",claa.name,tags);
		
			if(claa.default_map == true) add_tag("map","load",tags);
			
			let clp = clone_p(p,cl);

			if(clp != undefined) add_tag("clone",model.species[clp].name,tags);
			else{
				add_tag("index",claa.index,tags);
				if(cam.coord != COORD_DEFAULT) add_tag("coord",cam.coord,tags);
			}
			add_bscript("class",tags,bscript);
		
			if(clp == undefined){
				if(cam.grid == "on" || cam.slider.value != 0 || cam.coord == "latlng"){
					let tags=[];
					output_coords(cam.x,cam.y,cam,tags);
							
					if(cam.grid == "on") add_tag("grid","on",tags);
					
					if(cam.slider.value != 0) add_tag("comp-scale",Math.exp(cam.slider.value).toPrecision(2),tags);
					
					add_tag("scale",precision(cam.scale),tags);
						
					add_bscript("view",tags,bscript);			
				}
	
				let store=[];
				for(let c = 0; c < claa.ncomp; c++){
					let labs=[];
					
					let co = claa.comp[c];
					
					labs.push({tag:"name", val:co.name});
					
					labs.push({tag:"color", val:co.col});
					
					labs.push({tag:"fix", val:String(co.fixed.check)});
			
					if(sp.type == "Individual" && sp.trans_tree.check == true){
						if(sp.infection_cl.te == sp.cla[cl].name && co.infected.check != undefined){
							labs.push({tag:"infected", val:String(co.infected.check)});
						}
					}
						
					if(co.choose_branch == true){
						if(co.branch == true) labs.push({tag:"branch-prob", val:"true"});
						else labs.push({tag:"branch-prob", val:"false"});
					}
					
					switch(cam.coord){
					case "cartesian": 
						labs.push({tag:"x", val:precision(co.x/import_scale_factor)});
						labs.push({tag:"y", val:precision(co.y/import_scale_factor)});
						break;
						
					case "latlng":
						if(co.type == "boundary"){
							let k = find(map_store,"name",co.map_ref); 
							if(k == undefined) error("Cannot find feature");
		
							let map_da = { crs:{},features:[],type:"FeatureCollection"};
							map_da.features.push(generate_JSON_feature(co.name,map_store[k].feature));
							
							let file = JSON.stringify(map_da)+"\n";
							
							labs.push({tag:"boundary", val:file});
						}
						else{
							let p = transform_latlng_inv(co.x,co.y);
							labs.push({tag:"lat", val:precision(p.lat)});
							labs.push({tag:"lng", val:precision(p.lng)});
						}
						break;

					default: error("Option not recognised 111"); break;
					}
					
					store.push(labs);
				}
				
				if(store.length > 0){
					output_command_list("comp",claa.name,store,bscript);
				}
			}
		}
			
		{
			let store=[];
			for(let j = 0; j < claa.ntra; j++){
				let labs=[];
				
				let tr = claa.tra[j];
				
				let mid_start = 0, mid_end = tr.midp.length;
				let i_sel;
				
				let name = tr.name.replace(/→/g,"->")
			
				let splna = name.split("->");
				if(splna[0] == "+"){ i_sel = 0; mid_start++;}
				if(splna[1] == "-"){ i_sel = mid_end-1; mid_end--;}
				
				labs.push({tag:"name", val:name});
				
				if(i_sel != undefined){
					switch(cam.coord){
					case "cartesian":
						labs.push({tag:"x", val:precision(tr.midp[i_sel].x/import_scale_factor)});
						labs.push({tag:"y", val:precision(tr.midp[i_sel].y/import_scale_factor)});
						break;
					
					case "latlng":
						{
							let p = transform_latlng_inv(tr.midp[i_sel].x,tr.midp[i_sel].y);
							labs.push({tag:"lat", val:precision(p.lat)});
							labs.push({tag:"lng", val:precision(p.lng)});
						}
						break;

					default: error("Option not recognised 113"); break;
					}
				}
			
				let value = model.get_trans_value(tr);
		
				let tra_type = value.tra_type;
				let nm_flag = value.nm_flag;
				let err_fl = value.err_gl;
				let te = value.te;
				
				labs.push({tag:"value", val:te});
				
				if(err_fl != undefined){
					add_warning({mess:"Transition problem.", mess2:"The value for '"+err_fl+"' is not specified", p:p, cl:cl, i:j, warn_type:"TransMistake"});
				}
				
				if(nm_flag){
					switch(sp.type){
					case "Population": case "Deterministic":
						add_warning({mess:"Transition problem.", mess2:"Transition '"+tr.name+"' cannot have a '"+tra_type+"' distribution because '"+sp.name+"' uses a population-based model.", p:p, cl:cl, i:j, warn_type:"TransMistake"});
						break;
					}
				}
				
				let bp_te;
				
				if(tr.branch == true){
					if(tr.all_branches == true) bp_te = esc(tr.value.bp_eqn.te);
					else{
						if(tr.branch_select == true) bp_te = esc(tr.value.bp_eqn.te);
						else bp_te = "*";
					}
				}
				
				if(bp_te != undefined){
					if(bp_te == ""){
						add_warning({mess:"Branching problem", mess2:"Transition '"+tr.name+"' must have a branching probability assocaited with it.", p:p, cl:cl, i:j, warn_type:"TransMistake"});
					}
					else{
						labs.push({tag:"bp", val:bp_te});
					}
				}
				
				switch(cam.coord){
				case "cartesian":
					let xlist = '', ylist = '';
					
					for(let i = mid_start; i < mid_end; i++){
						if(i != mid_start){ xlist += ','; ylist += ',';}
						
						xlist += precision(tr.midp[i].x/import_scale_factor);
						ylist += precision(tr.midp[i].y/import_scale_factor);
					}
					
					if(xlist != ''){
						labs.push({tag:"mid-x", val:xlist});
						labs.push({tag:"mid-y", val:ylist});
					}
					break;
					
				case "latlng":
					let latlist = '', lnglist = '';
					
					for(let i = mid_start; i < mid_end; i++){
						if(i != mid_start){ latlist += ','; lnglist += ',';}
						
						let p = transform_latlng_inv(tr.midp[i].x,tr.midp[i].y);

						latlist += precision(p.lat);
						lnglist += precision(p.lng);
					}
					
					if(latlist != ''){
						labs.push({tag:"mid-lat", val:latlist});
						labs.push({tag:"mid-lng", val:lnglist});
					}
					break;

				default: error("Option not recognised 115"); break;
				}
				
				store.push(labs);
			}
			
			if(store.length > 0){
				output_command_list("trans",claa.name,store,bscript);
			}
		}
		
		if(claa.annotation.length > 0){
			for(let i = 0; i < claa.annotation.length; i++){
				let anno = claa.annotation[i];
				switch(anno.type){
				case "box":
					{
						let tags=[];
						
						let st = "";
						for(let k = 0; k < anno.comps.length; k++){
							if(k != 0) st += ",";
							st += anno.comps[k];
						}						
						add_tag("comps",st,tags);
						
						if(anno.color != annotation_col_default) add_tag("color",anno.color,tags); 
						if(anno.size != size_annotation_default) add_tag("text-size",anno.size,tags);  
						
						add_bscript("box",tags,bscript);
					}
					break;
					
				case "text":
					{
						let tags=[];
						add_tag("text",anno.te,tags);
						output_coords(anno.x,anno.y,cam,tags);
						if(anno.color != annotation_col_default) add_tag("color",anno.color,tags);
						if(anno.size != size_annotation_default) add_tag("text-size",anno.size,tags);
					
						add_bscript("label",tags,bscript);
					}
					break;
					
				case "map":
					if(anno.default_map != true){
						let k = find(map_store,"name",anno.map_ref);
						if(k == undefined) error("map undefined");
						else{
							let map_da = { crs:{},features:[],type:"FeatureCollection"};
							for(let j = 0; j < map_store[k].feature.length; j++){
								map_da.features.push(generate_JSON_feature("Feature"+k,map_store[k].feature[j]));
							}
						
							let file = JSON.stringify(map_da)+"\n";
							
							let tags=[];
							add_tag("file",file,tags);
							add_bscript("map",tags,bscript);
						}
					}
					break;
					
				default:
					break;
				}
			}
		}
	}
}

	
/// Gets string for coordinates
function output_coords(x,y,cam,tags)
{
	switch(cam.coord){
	case "cartesian": 
		add_tag("x",precision(x),tags);
		add_tag("y",precision(y),tags);
		break;
		
	case "latlng": 
		let p = transform_latlng_inv(x,y);
		add_tag("lat",precision(p.lat),tags);
		add_tag("lng",precision(p.lng),tags);
		break;

	default: error("Option not recognised 110"); break;
	}
}

	
/// Outputs a list of commands either seperately or in a file (for comps-all and trans-all)
function output_command_list(type,claa_name,store,bscript)
{
	let fl = false;
	if(store.length >= 10){
		fl = true;
		for(let j = 0; j < store.length; j++){
			let labs = store[j];
			for(let i = 0; i < labs.length; i++){
				if(labs[i].tag == "boundary") fl = false;
			}
		}
	}
	
	if(fl == true){	
		let head=[];
		for(let j = 0; j < store.length; j++){
			let labs = store[j];
			for(let i = 0; i < labs.length; i++){
				let ta = labs[i].tag;
				let j = find_in(head,ta);
				if(j == undefined){ j = head.length; head.push(ta);}
				labs[i].col = j;
			}
		}
		
		let row=[];
		let N = head.length;
		for(let j = 0; j < store.length; j++){
			let ro=[];
			for(let i = 0; i < N; i++) ro[i] = "";
			
			let labs = store[j];
			for(let i = 0; i < labs.length; i++){
				ro[labs[i].col] = labs[i].val;
			}
			
			row.push(ro);
		}
		
		let tags=[];
		
		// Looks to remove tags which are all the same
		let i = 0;
		while(i < N){
			let val = row[0][i];
			
			let j = 1; while(j < store.length && row[j][i] == row[0][i]) j++;
			
			if(j == store.length){
				add_tag(head[i],val,tags);
				
				head.splice(i,1);
				for(let j = 0; j < store.length; j++){
					row[j].splice(i,1);
				}
				N--;
			}
			else i++;
		}
		
		add_tag("file",output_table_simp_file(head,row),tags);
	
		add_bscript(type+"-all",tags,bscript);
	}	
	else{		
		for(let j = 0; j < store.length; j++){
			let labs = store[j];

			let tags=[];
			for(let k = 0; k < labs.length; k++){
				if(!(labs[k].tag == "fix" && labs[k].val == "false")){
					add_tag(labs[k].tag,labs[k].val,tags);
				}
			}
			
			add_bscript(type,tags,bscript);
		}
	}
}


/// Determines if not a number
function is_nan(val)
{
	if(isNaN(val)) return true;
	if(typeof(val) == "boolean") return true;
	return false
}


/// Determines if cloning is possible because classification compartments already given previously
function clone_p(p,cl)
{
	let claa = model.species[p].cla[cl];
	let cam = claa.camera;
	for(let p2 = 0; p2 < p; p2++){
		let sp2 = model.species[p2];
		for(let cl2 = 0; cl2 < sp2.ncla; cl2++){
			let claa2 = sp2.cla[cl2];
			if(claa2.name == claa.name){
				if(claa2.ncomp == claa.ncomp){
					let prob = 0;
					for(let c = 0; c < claa.ncomp; c++){
						let comp = claa.comp[c];
						let comp2 = claa2.comp[c];
						if(comp.name != comp2.name) prob = 1;
						
						if(cam.coord == "latlng" && comp.type == "boundary"){}
						else{
							if(comp.x != comp2.x || comp.y != comp2.y) prob = 2;
						}
						
						if(comp.col != comp2.col) prob = 3;
						
						if(comp.fixed.check != comp2.fixed.check) prob = 4;
					}
					if(prob == 0) return p2;
				}
			}
		}			
	}
}
 

/// Create output of parameters
function create_output_param(save_type,bscript)
{
	if(model.param.length == 0) return "";
	
	banner("PARAMETERS",bscript);
	
	for(let th = 0; th < model.param.length; th++){
		let par = model.param[th];
		if(par.variety == "define") output_define(par,bscript)
	}
	
	for(let th = 0; th < model.param.length; th++){
		let par = model.param[th];
		if(par.variety != "define"){
			output_param(model.param[th],save_type,bscript);
		}
	}
}


/// Outputs any parameter definitions
function output_define(par,bscript)
{
	if(par.variety != "define") error("Equation should be define");
	
	let tags=[];
	add_tag("name",remove_eq_quote(par.full_name),tags);
	add_tag("value",par.define_eqn.trim(),tags);

	add_bscript("define",tags,bscript);
}
		
		
/// Outputs a single parameter
function output_param(par,save_type,bscript)
{
	let dist_done = false;

	let tags=[];
	
	if(par.type != "derive_param"){		
		let num_warn = model.warn.length;
	
		let com = "param";
		if(par.type == "param factor"){
			let pfac = model.param_factor;
		
			let k = 0, kmax = pfac.length; 

			while(k < kmax && pfac[k].f_param.full_name != par.full_name) k++;
			if(k == kmax) alert_help("Problem loading param factor");
	
			com += "-mult";
			add_tag("name",remove_eq_quote(pfac[k].param.full_name),tags);
		}
		else{
			add_tag("name",remove_eq_quote(par.full_name),tags);
		}
		
		let display = false;
		
		if(par.dist_mat || par.iden_mat || par.den_vec){
			display = true;
		}
		else{
			// Constants and simulation parameters
			if(par.variety == "reparam" && par.reparam_eqn_on){	
				if(par.reparam_eqn.trim() != ""){
					add_tag("reparam",par.reparam_eqn,tags);
				}
				display = true;
			}
			else{
				if((par.dep.length == 0 && par.value == set_str) || (par.set == false && !(par.variety == "dist" && par.sim_sample.check == true))){
					if(save_type == "sim"){
						if(find_in(sim_param_not_needed,par.type) == undefined){
							add_warning({mess:"A value for parameter "+par.full_name+" must be set.", mess2:"Parameter values must be set before simulation can be perfomed", warn_type:"MissingSimValue", name:par.name});
						}
					}
					
					if(save_type == "ppc"){
						add_warning({mess:"A value for parameter multiplier "+par.full_name+" must be set.", mess2:"Parameter multiplier values must be set before posterior simulation can be perfomed", warn_type:"MissingParamMultValue", name:par.name});
					}
				}
				else{
					if(par.variety == "dist" && par.sim_sample.check == true){
						output_add_prior_distribution(save_type,par,"dist",tags);
						display = true;
						dist_done = true;
					}
					else{
						let exp = "value";
						if(par.variety == "const") exp = 'constant';
						else{
							if(par.variety == "reparam") exp = 'reparam';
						}
									
						let value;
						if(par.dep.length > 0){
							value = output_value_table(par,par.value,"Value","value",save_type);
						}
						else{
							value = par.value;
						}
					
						add_tag(exp,value,tags);
						display = true;
					}	
				}
			}
			
			if(par.prior_const_on == true && par.prior_const_set == true){
				output_add_prior_const(save_type,par,tags);
				display = true;
			}
			else{
				if(par.variety != "const" && par.variety != "reparam"){
					switch(par.variety){
					case "normal":
						output_add_prior_distribution(save_type,par,"prior",tags);
						display = true;
						break;
						
					case "dist":
						if(dist_done == false){
							output_add_prior_distribution(save_type,par,"dist",tags);
							display = true;
						}
						break;
						
					default: error("Option not recognised 1B"+par.variety); break;
					}
				}
			}
			
			if(par.spline.on == true){
				if(par.spline.knot.length == 0) error("There should be knot times specified");
				
				// Checks to see if spline is within range
				if(save_type == "sim" || save_type == "inf"){	
					let details = model.sim_details;
					if(save_type == "inf") details =  model.inf_details;
				
					let t_start = Number(details.t_start);
					let t_end = Number(details.t_end);
					for(let j = 0; j < par.spline.knot.length; j++){
						let num = par.spline.knot[j];
						if(num != "start" && num != "end"){
							/*
							num = Number(num);
							if(num < t_start || num > t_end){
								add_warning({mess:"Knot problem", mess2:"Knot time '"+num+"' must be between the start and end times", warn_type:"KnotProblem", name:par.name});
							}
							*/
						}
					}
				}
					
				add_tag("knot-times",stringify(par.spline.knot),tags);
				
				let sm = par.spline.smooth;
				if(sm.check == true){
					add_tag("smooth",sm.type.value.toLowerCase()+'('+sm.value+')',tags);
				}
				
				let sel = par.spline.spline_radio.value;

				if(sel != "Linear"){
					add_tag("spline-type",sel.toLowerCase(),tags); 
				}
			}
		}
		
		if(par.factor){
			add_tag("factor","true",tags);
		
			if(par.factor_weight_on.check == true){
				let file = output_value_table(par,par.factor_weight,"Value","value",save_type);
				add_tag("factor-weight",file,tags);
			}
		}
		
		if(par.sim_sample && par.sim_sample.check == false){
			add_tag("sim-sample","false",tags);
		}
				
		if(!par.dist_mat && !par.iden_mat && !par.den_vec){
			add_bscript(com,tags,bscript);
		}
		
		if(display == false && (save_type == "sim" || save_type == "inf")){
			if(num_warn == model.warn.length){
				add_warning({mess:"Parameter "+par.full_name+" could not be exported.", mess2:"Unknown error", warn_type:"", name:par.name});
			}
		}
	}
}
	
	
/// Outputs the value of a table
function output_value_table(par,value,head_col,key,save_type)
{
	let list = par.list;
	
	if(value == undefined){
		let wt, ty;
		switch(head_col){
		case "Dist": 
			wt = "MissingDistValue"; 
			if(save_type != "sim" && save_type != "inf") return; 	
			break;
			
		case "Prior":
			wt = "MissingPriorSplitValue"; 
			if(save_type != "inf") return; 	
			break;
			
		case "Value":
			wt = "MissingSimValue";
			if(save_type != "sim" && save_type != "inf") return; 	
			break;
			
		case "ValueInf":
			wt = "MissingPriorConst";
			if(save_type != "sim" && save_type != "inf") return; 	
			break;
		}
		
		add_warning({mess:"Missing "+key, mess2:"The "+key+" for "+par.full_name+" must be set.", warn_type:wt, name:par.name});
		return;
	}
						
	let data = '';
	if(is_matrix(par)){   // Matrix format
		for(let i = 0; i < list[1].length; i++){
			if(i != 0) data += ",";
			data += '"'+list[1][i]+'"';
		}
		data += endl;
		
		for(let j = 0; j < list[0].length; j++){
			for(let i = 0; i < list[1].length; i++){
				if(i != 0) data += ",";
				
				let index = []; index[0] = j; index[1] = i;
				let el = get_element(value,index);
				if(el == undefined){ param_undefined(par,index); return;}
					
				let val = get_op_val(el,head_col,par);
				if(isNaN(val)) data += '"'+val+'"';
				else data += val;
			}
			data += endl;
		}
	}
	else{
		for(let d = 0; d < par.dep.length; d++){	
			data += '"'+par.dep[d]+'",';
		}
		data += '"'+head_col+'"'+endl;

		let co_list = generate_co_list(par.list);
		for(let i = 0; i < co_list.length; i++){
			let comb = co_list[i];
			let el = get_element(value,comb.index);
			if(el == undefined){ param_undefined(par,comb.index); return;}
			
			if(el != 0){	
				for(let d = 0; d < par.dep.length; d++){
					data += '"'+list[d][comb.index[d]]+'",';
				}
				
				let val = get_op_val(el,head_col,par);
				if(isNaN(val)) data += '"'+val+'"';
				else data += val;
				data += endl;
			}
		}
	}

	return data;
}


/// For output convert an element to output text
function get_op_val(el,head_col,par)
{
	if(head_col == "Dist" || head_col == "Prior"){
		if(el.type.te == select_str){
			add_warning({mess:"Missing "+key, mess2:"The "+key+" for "+par.full_name+" must be set.", warn_type:"MissingPriorSplitValue", name:par.name});
			return "";
		}
	
		return get_prior_string(el);
	}
	
	return el;
}

					
/// Outputs a warning message because parameter undefined
function param_undefined(par,index)
{
	let me = "A value for parameter "+par.full_name+" is undefined (";
	for(let j = 0; j < par.dep.length; j++){
		if(j != 0) me += ",";
		me += list[j][index[j]];
	}
	me += ").";
				
	add_warning({mess:"Parameter value is undefined", mess2:me, warn_type:"ParamValueProblem", par:par});
}


/// Adds text associated with adding a distribution
function output_add_prior_distribution(save_type,par,variety,tags)
{
	let key;
	switch(variety){
	case "dist": key = "distribution"; break;
	case "prior": key = "prior"; break;
	}
		
	if(par.dep.length > 0 && par.prior_split_check.check == true && !par.factor){  /// Prior split up
		let warn_type;
		switch(variety){
		case "dist": key = "distribution"; warn_type = "MissingDistValue"; break;
		case "prior": key = "prior"; warn_type = "MissingPriorValue";  break;
		default: error("option prob"); break;
		}
			
		if(par.prior_split_set == false){
			add_warning({mess:"Missing "+key+" distribution information", mess2:"The "+key+" for parameter "+par.full_name+" must be set.", warn_type:warn_type, name:par.name});
		}
		else{
			let ta, head_col;
			switch(variety){
			case "dist": ta = "dist-split"; head_col = "Dist"; break;
			case "prior": ta = "prior-split"; head_col = "Prior"; break;
			default: error("option prob"); break;
			}
			
			let file = output_value_table(par,par.prior_split,head_col,key,save_type);
			if(file == undefined) return "";
			
			add_tag(ta,file,tags);
		}
	}
	else{
		let pri = par.prior;
		
		let prior_string = get_prior_output(save_type,par.prior,variety,par);
		if(prior_string == "") return "";
	
		let ta;
		switch(variety){ 
		case "dist": ta = "dist"; break;
		case "prior": ta = "prior"; break;
		default: error("option prob"); break;
		}
		
		add_tag(ta,prior_string,tags);
	}
}


/// Adds a prior constant
function output_add_prior_const(save_type,par,tags)
{
	let te = ' prior-const=';
		
	let value;
	if(par.dep.length > 0){
		let err = check_param_value("Set Param",par,par.prior_const);
		if(err){
			let wt = "MissingSimValue"; if(par.variety == "reparam") wt = "ReparamValue";
			
			add_warning({mess:"Problem with "+par.full_name+" value", mess2:err, warn_type:wt, name:par.name});
			return;
		}
	
		value = output_value_table(par,par.prior_const,"Value","value",save_type);
	}
	else{
		value = JSON.stringify(par.prior_const); 
	}
	
	add_tag("prior-const",value,tags);
}


/// Gets the output string from a prior definition
function get_prior_output(save_type,pri,variety,par)
{
	let tex = pri.type.te;
	if(tex == select_str){
		if(save_type == "inf"){	
			add_warning({mess:"Missing "+variety, mess2:"The "+variety+" for parameter "+par.full_name+" must be set.", warn_type:"MissingPriorValue", name:par.name});
		}
		return "";
	}
	
	let st = get_prior_string(pri);
	if(save_type == "inf" && par){	
		let dist = true; if(variety == "prior") dist = false;
		let pri = convert_text_to_prior(st,par.pri_pos,dist);
	
		if(pri.err == true){
			add_warning({mess:"Error for "+variety, mess2:"The "+variety+" for parameter "+par.full_name+" has an error: "+pri.msg, warn_type:"ErrorPriorValue", name:par.name});
		}
	}
	
	return st;
}


/// Checks that an output for a 
function check_prior_val(val,variety,par,op)
{
	let key, warn_type, warn_text;

	switch(variety){
	case "prior":
		key = "prior"; 
		warn_type = "MissingPriorValue"; 
	
		if(isNaN(val)) warn_text = "Not a value";
		else{
			if(op.positive == true){
				if(Number(val) <= 0) warn_text = "Must be positive2";
			}
		}
		break;

	case "dist":
		{
			key = "distribution"; warn_type = "MissingDistValue"; 
			let res = is_eqn(val,"min",op);
			if(res.err == true) warn_text = "'"+val+"' is not valid.";
		}
		break;
		
	default: error("not op"); break;
	}

	if(warn_text != undefined){
		add_warning({mess:"Error in "+key+" for "+par.full_name, mess2:warn_text, warn_type:warn_type, name:par.name});
	}
}

	
/// Create output of derived quantities
function create_output_derived(bscript)
{
	let derive = model.derive;

	if(derive.length == 0) return "";

	mini_banner("DERIVED",bscript);

	for(let i = 0; i < derive.length; i++){
		let der = derive[i];
		
		if(der.eqn1.param.length != 1) error("problem with derived");
		else{
			let name = remove_eq_quote(der.eqn1.param[0].full_name);	
			
			let tags=[];
			add_tag("name",name,tags);
			add_tag("eqn",der.eqn2.te,tags);
			add_bscript("derived",tags,bscript);
		}
	}
}


/// Create output of individual effects
function create_output_ind_eff(p,bscript)
{
	let iegs = model.species[p].ind_eff_group;
		
	if(iegs.length == 0) return "";
	
	mini_banner("INDIVIDUAL EFFECTS",bscript);

	for(let i = 0; i < iegs.length; i++){
		let ieg = iegs[i];
		let ie = "";
		for(let j = 0; j < ieg.ie_list.length; j++){
			if(j > 0) ie += ",";
			ie += ieg.ie_list[j].name;
		}
		
		let name = ieg.name;
		
		let tags=[];
		
		add_tag("name",ieg.name,tags);
		add_tag("ie",ie,tags);
		
		let mat = ieg.A_matrix;
		
		if(mat.check == true){
			if(mat.loaded == false){
				add_warning({mess:"Error in individual effects", mess2:"The A matrix for individual effects '"+name+"' must be loaded", warn_type:"A matrix"});
			}
			else{
				switch(mat.type){
				case "pedigree":
					{
						let data = 'ID,sire,dam'+endl;
						for(let c = 0; c < mat.ind_list.length; c++){ 
							data += '"'+mat.ind_list[c]+'"'+",";
							data += '"'+mat.sire_list[c]+'"'+",";
							data += '"'+mat.dam_list[c]+'"'+endl;
						}
							
						add_tag("pedigree",data,tags);
					}
					break;
					
				case "A": 
					{	
						if(true){                                   // Sparse representation
							let data_ind = 'Individual'+endl;
							for(let c = 0; c < mat.ind_list.length; c++){ 
								data_ind += '"'+mat.ind_list[c]+'"'+endl;
							}
						
							add_tag("ind-list",data_ind,tags);
							
							let data = 'j,i,value'+endl;
							for(let r = 0; r <  mat.ind_list.length; r++){
								for(let c = 0; c <= r; c++){
									let val = mat.A_value[r][c];
									
									if(val != 0) data += r+","+c+","+val+endl;
								}
							}
				
							add_tag("A-sparse",data,tags);
						}
						else{	
							let data = '';
							for(let c = 0; c < mat.ind_list.length; c++){ 
								if(c != 0) data += ',';
								data += '"'+mat.ind_list[c]+'"';
							}
							data += endl;

							for(let r = 0; r <  mat.ind_list.length; r++){
								for(let c = 0; c < mat.ind_list.length; c++){ 
									if(c != 0) data += ',';
									data += mat.A_value[r][c];
								}
								data += endl;
							}
				
							add_tag("A",data,tags);
						}
					}
					break;
					
				case "Ainv": 
					{
						let data = '';
						for(let c = 0; c < mat.ind_list.length; c++){ 
							if(c != 0) data += ',';
							data += '"'+mat.ind_list[c]+'"';
						}
						data += endl;

						for(let r = 0; r <  mat.ind_list.length; r++){
							for(let c = 0; c < mat.ind_list.length; c++){ 
								if(c != 0) data += ',';
								data += mat.A_value[r][c];
							}
							data += endl;
						}
				
						add_tag("Ainv",data,tags);
					}
					break;
				}
			}
		}	

		add_bscript("ind-effect",tags,bscript);
	}
}


/// Create output of fixed effects
function create_output_fix_eff(p,save_type,bscript)
{
	let fegs = model.species[p].fix_eff;
	
	if(fegs.length == 0) return "";
	
	mini_banner("FIXED EFFECTS",bscript);

	for(let i = 0; i < fegs.length; i++){
		let feg = fegs[i];
		let name = feg.name;
		
		let tags=[];
		
		add_tag("name",feg.name,tags);
		
		let vec = feg.X_vector;
		if(vec.loaded == false){
			if(save_type == "sim" || save_type == "inf"){
				add_warning({mess:"Error in fixed effects", mess2:"The vector <e>X^"+name+"</e> for fixed effect '"+name+"' must be loaded", warn_type:"X vector"});
			}
		}
		else{
			let data = 'ID,value'+endl;
			for(let c = 0; c < vec.ind_list.length; c++){ 
				data += '"'+vec.ind_list[c]+'",'+vec.X_value[c]+endl;
			}
			add_tag("X",data,tags);
		}
		
		add_bscript("fixed-effect",tags,bscript);
	}
}

	
/// Creates information about simulation, inference or posterior simulation
function create_output_siminf(save_type,bscript)
{
	{                                                // Deals with sim details
		let details = model.sim_details;
		
		if(save_type ==  "sim"){
			if(details.t_start == "" || isNaN(details.t_start)){
				add_warning({mess:"Error in simulation details", mess2:"The start time must be set", warn_type:"SimDetails"});
			}
		
			if(details.t_end == "" || isNaN(details.t_end)){
				add_warning({mess:"Error in simulation details", mess2:"The end time must be set", warn_type:"SimDetails"});
			}
		
			let num = Number(details.number);
			if(details.number == "" || isNaN(num) || num != Math.floor(num) || num <= 0){
				add_warning({mess:"Error in simulation details", mess2:"The simulation number must be set to a positive integer", warn_type:"SimDetails"});
			}
		
			let timestep = Number(details.timestep);
			if(details.timestep == "" || isNaN(timestep) || timestep <= 0){
				add_warning({mess:"Error in simulation details", mess2:"The simulation timestep must be set to a positive number", warn_type:"SimDetails"});
			}
		}
		
		let tags=[];
		if(details.t_start != "") add_tag("start",details.t_start,tags);
		if(details.t_end != "") add_tag("end",details.t_end,tags);
		
		if(details.number != ""){
			let num = Number(details.number);
			if(num != SIM_NUM_DEFAULT) add_tag("number",num,tags);
		}
		
		if(details.timestep != "") add_tag("timestep",details.timestep,tags); 
		
		if(details.seed_on.value == "Yes") add_tag("seed",details.seed,tags); 
		 
		if(details.indmax != INDMAX_DEFAULT) add_tag("ind-max",details.indmax,tags); 
		
		if(details.param_output_max != PARAM_OUTPUT_MAX_DEFAULT){
			add_tag("param-output-max",details.param_output_max,tags); 
		}
		
		if(details.optimise.value != "auto") add_tag("optimise",details.optimise.value,tags); 
		
		if(details.compress.value != "auto") add_tag("compress",details.compress.value,tags); 
		
		add_bscript("simulation",tags,bscript);
	}

	{                                                // Deals with inf details
		let details = model.inf_details;
		
		if(save_type == "inf"){
			if(details.t_start == "" || isNaN(details.t_start)){
				add_warning({mess:"Error in inference details", mess2:"The start time must be set", warn_type:"InfDetails"});
			}
		
			if(details.t_end == "" || isNaN(details.t_end)){
				add_warning({mess:"Error in inference details", mess2:"The end time must be set", warn_type:"InfDetails"});
			}

			let timestep = Number(details.timestep);
			if(details.timestep == "" || isNaN(timestep) || timestep <= 0){
				add_warning({mess:"Error in inference details", mess2:"The timestep must be set to a positive number", warn_type:"InfDetails"});
			}
		}
		
		let tags=[];
		
		if(details.t_start != "") add_tag("start",details.t_start,tags); 
		if(details.t_end != "") add_tag("end",details.t_end,tags);
		if(details.timestep != "") add_tag("timestep",details.timestep,tags); 
		
		if(details.seed_on.value == "Yes") add_tag("seed",details.seed,tags);
		
		let alg = details.algorithm.value;
		
		switch(alg){
		case "DA-MCMC": case "PAS-MCMC":
			{
				let num_samp = Number(details.sample);
				if(num_samp != MCMC_SAMPLE_DEFAULT) add_tag("update",num_samp,tags); 
				
				let num_param = Number(details.output_param);
				if(num_param != MCMC_OP_PARAM_DEFAULT) add_tag("param-output",num_param,tags); 
			
				let num_state = Number(details.output_state);
				if(num_state != MCMC_OP_STATE_DEFAULT) add_tag("state-output",num_state,tags); 
			
				if(alg == "DA-MCMC"){
					let cha_per_core = Number(details.cha_per_core);
					if(cha_per_core != MCMC_CHAIN_PER_CORE_DEFAULT) add_tag("chain-per-core",cha_per_core,tags); 
				
					if(Number(details.nchain)%cha_per_core != 0){
						add_warning({mess:"Error with chain number", mess2:"The number of chains '"+details.nchain+"' must be a multiple of the number of chains per core '"+cha_per_core+"'", warn_type:"InfFurtherDetails"});
					}
				}
				
				if(alg == "PAS-MCMC"){
					let part_per_core = Number(details.part_per_core);
					if(part_per_core != PAS_PART_PER_CORE_DEFAULT) add_tag("part-per-core",part_per_core,tags); 
					
					if(Number(details.npart)%part_per_core != 0){
						add_warning({mess:"Error with particle number", mess2:"The number of particles '"+details.npart+"' must be a multiple of the number of particles per core '"+part_per_core+"'", warn_type:"InfFurtherDetails"});
					}
				}
			}
			break;
		
		case "ABC":
			add_tag("sample",details.abcsample,tags); 
			add_tag("acc-frac",details.accfrac,tags); 
			break;
			
		case "ABC-SMC":
			add_tag("sample",details.abcsample,tags); 
			add_tag("acc-frac",details.accfracsmc,tags); 
			add_tag("gen",details.numgen,tags); 
			add_tag("kernel-size",details.kernelsize,tags); 
			break;
		}
		
		if(details.indmax != INDMAX_DEFAULT) add_tag("ind-max",details.indmax,tags);
			
		let poutmax = Number(details.param_output_max);
		if(poutmax != PARAM_OUTPUT_MAX_DEFAULT) add_tag("param-output-max",poutmax,tags); 
		
		if(details.algorithm.value != ALG_DEFAULT) add_tag("algorithm",details.algorithm.value,tags); 
		
		switch(details.algorithm.value){
		case "DA-MCMC":
			{	
				add_tag("nchain",details.nchain,tags); 
				
				let type = details.anneal_type.te;
				
				if(type != ANNEAL_DEFAULT){
					add_tag("anneal",type,tags); 
					switch(type){
					case "scan": add_tag("rate",details.anneal_rate,tags); break;	
					case "power": add_tag("power",details.anneal_power,tags); break;
					}
				}
				
				if(type != "scan"){
					let nbf = Number(details.burnin_frac);
					if(nbf != BURNIN_FRAC_DEFAULT) add_tag("burnin-percent",details.burnin_frac,tags);
				}
			}
			break;
				
		case "PAS-MCMC":
			{
				add_tag("npart",details.npart,tags); 
				add_tag("gen-percent",details.gen_update,tags); 
				let nbf = Number(details.burnin_frac);
				if(nbf != BURNIN_FRAC_DEFAULT) add_tag("burnin-percent",details.burnin_frac,tags); 
			}
			break;
		}
	
		switch(details.algorithm.value){
		case "DA-MCMC": case "PAS-MCMC":
			if(details.sync_on.value == "Off") add_tag("sync","Off",tags);  
			break;
		}
		
		if(details.optimise.value != "auto") add_tag("optimise",details.optimise.value,tags); 
		
		if(details.compress.value != "auto") add_tag("compress",details.compress.value,tags); 
		
		add_bscript("inference",tags,bscript);
	}
	
	{                                               // Deals with ppc details
		let tags = get_ppc_tags(save_type);
		add_bscript("post-sim",tags,bscript);
	}
}


/// Creates tags for the ppc command
function get_ppc_tags(save_type)
{
	let details = model.ppc_details;
		
	if(save_type == "ppc"){
		if(details.ppc_t_start == "" || isNaN(details.ppc_t_start)){
			add_warning({mess:"Error in simulation details", mess2:"The start time must be set", warn_type:"SimDetails"});
		}
	
		if(details.ppc_t_end == "" || isNaN(details.ppc_t_end)){
			add_warning({mess:"Error in simulation details", mess2:"The end time must be set", warn_type:"SimDetails"});
		}
	
		let num = Number(details.number);
		if(details.number == "" || isNaN(num) || num != Math.floor(num) || num <= 0){
			add_warning({mess:"Error in simulation details", mess2:"The simulation number must be set to a positive integer", warn_type:"SimDetails"});
		}
	}
	
	let tags=[];
	
	if(details.ppc_t_start != "") add_tag("start",details.ppc_t_start,tags);
	if(details.ppc_t_end != "") add_tag("end",details.ppc_t_end,tags);
	
	if(details.number != ""){
		let num = Number(details.number);
		if(num != PPC_NUM_DEFAULT) add_tag("number",details.number,tags);
	}
	
	if(details.seed_on.value == "Yes") add_tag("seed",details.seed,tags);
	
	if(details.optimise.value != "auto") add_tag("optimise",details.optimise.value,tags);
		
	if(details.compress.value != "auto") add_tag("compress",details.compress.value,tags);
	
	if(details.run_inf_model.value == "No") add_tag("param-only","true",tags);
	
	switch(details.run_post.value){
	case "postsample": break;
	case "postmean": add_tag("param-sample","mean",tags); break;
	case "postmedian": add_tag("param-sample","median",tags); break;
	default: error("Problem with option"); break;
	}
	
	let cbl = details.check_box_list;

	if(cbl){
		let st = "";
		let fl = false;
		for(let k = 0; k < cbl.length; k++){
			if(cbl[k].checkb.check == true){
				if(fl) st += ",";
				st += cbl[k].name_simp;
				fl = true;
			}
		}		
		
		if(fl) add_tag("resample",st,tags);
	}
	
	return tags;
}


/// Adds information about simulation and inference
function create_output_sim_inf_source(p,sim_or_inf,bscript)
{
	let sp = model.species[p];
	if(sim_or_inf == "ppc"){
		if(!model.inf_res.plot_filter) return;
		sp = model.inf_res.plot_filter.species[p];
		if(!sp) return;
	}
	
	let source, pa;
	
	switch(sim_or_inf){
	case "sim": 
		pa = "Simulation"; 
		source = sp.sim_source; if(source.length == 0) return "";
		mini_banner("SIMULATION INITIAL CONDITIONS",bscript); 
		break;
		
	case "inf": 
		pa = "Inference"; 
		source = sp.inf_source; if(source.length == 0) return "";
		mini_banner("INFERENCE DATA",bscript); 
		break;
		
	case "ppc": 
		pa = "Population Mod."; 
		source = sp.ppc_source; if(source.length == 0) return ""; 
		mini_banner("MODIFICATION DATA",bscript); 
		break;
		
	default: error("Option not recognised 116"+sim_or_inf); break;
	}

	for(let i = 0; i < source.length; i++){
		let so = source[i];
	
		so.check_info = {siminf:sim_or_inf, p:p, i:i}; 

		data_source_check_error("warn",so);

		if(so.error == true){
		}
		else{
			let add_flag = false;
			
			switch(so.type){
			case "Init. Pop.":
				{
					let com = "init-pop-inf"; if(sim_or_inf == "sim") com = "init-pop-sim";
		
					let ic_type = so.spec.radio_dist.value;
					
					let head = [];
					let row = [];
					
					let tags=[];
					
					switch(ic_type){
					case "Fixed": add_tag("type","fixed",tags); break;
					case "Dist": add_tag("type","dist",tags); break;
					}
					
					switch(so.spec.radio2.value){
					case "Focal":
						if(so.cla.length != sp.ncla){
							so.error = true;
							add_warning({mess:"Initial population problem", mess2:"The number of classifications in the initial population for '"+sp.name+"' is not right",warn_type:"Init_pop", pa:pa, p:p, i:i});
						}
						else{
							for(let cl = 0; cl < sp.ncla; cl++){
								let claa = sp.cla[cl];

								let ip = so.cla[cl].comp_init_pop;
								if(ic_type == "Dist") ip = so.cla[cl].comp_init_dist;
								
								if(ip.length != claa.ncomp){
									so.error = true;
									add_warning({mess:"Initial population problem", mess2:"The number of compartments in classification '"+claa.name+"' is not right",warn_type:"Init_pop", pa:pa, p:p, i:i});
								}
								else{
									for(let c = 0; c < ip.length; c++){
										if(ip[c].comp_name_store != claa.comp[c].name){
											if(ip[c].name == undefined){
												so.error = true;
												add_warning({mess:"Initial population problem", mess2:"The data source has become invalid",warn_type:"Init_pop", pa:pa, p:p, i:i});
											}
											else{
												so.error = true;
												add_warning({mess:"Initial population problem", mess2:"The names for the compartments '"+ip[c].name+"' and '"+claa.comp[c].name+"' do not agree",warn_type:"Init_pop", pa:pa, p:p, i:i});
											}
										}
									}
								}
							}
						}
						
						if(model.warn.length == 0){
							head.push("Compartment");
							switch(ic_type){
							case "Fixed": head.push("Population"); break;
							case "Dist": head.push("Distribution"); break;
							}
						
							for(let cl = 0; cl < sp.ncla; cl++){
								let claa = sp.cla[cl];
								let cmax = claa.ncomp;
								if(ic_type == "Fixed" && claa.name != so.spec.focal.te) cmax--;
								
								for(let c = 0; c < cmax; c++){
									let ele = [];
									ele.push(claa.comp[c].name);
									
									switch(ic_type){
									case "Fixed":
										if(claa.name == so.spec.focal.te){
											ele.push(so.cla[cl].comp_init_pop[c].pop);
										}
										else{
											ele.push(so.cla[cl].comp_init_pop[c].pop_per+"%");
										}
										break;
										
									case "Dist":										
										if(claa.name == so.spec.focal.te){
											let val = get_prior_string(so.cla[cl].comp_init_dist[c].dist);
											ele.push(val);
										}
										else{
											ele.push(so.cla[cl].comp_init_dist[c].alpha);
										}
										break;
									}
									
									row.push(ele);
								}
							}
							
							add_tag("focal",so.spec.focal.te,tags);
							
							let tab = output_table_simp_file(head,row);
							add_tag("file",tab,tags);
						}
						break;
							
					case "All":
						if(so.spec.radio_dist.value == "Dist"){
							add_tag("prior",get_prior_output(sim_or_inf,so.pop_dist,"prior",undefined),tags);
						}
						
						for(let cl = 0; cl < sp.ncla; cl++) head.push(sp.cla[cl].name);
						switch(ic_type){
						case "Fixed": head.push("Population"); break;
						case "Dist": head.push("Alpha"); break;
						}
						
						for(let j = 0; j < so.glob_comp.length; j++){
							let gc = so.glob_comp[j];
							if(gc.cla.length != sp.ncla){
								so.error = true;
								add_warning({mess:"Initial population problem", mess2:"The number of classifications in the initial population for '"+sp.name+"' is not right",warn_type:"Init_pop", pa:pa, p:p, i:i});
								return;
							}
					
							let ele=[];
							for(let cl = 0; cl < sp.ncla; cl++) ele.push(gc.cla[cl]);
							
							switch(ic_type){
							case "Fixed": ele.push(gc.pop); break;
							case "Dist": ele.push(gc.alpha); break;
							}
							row.push(ele);
						}
					
						let tab = output_table_simp_file(head,row);
						add_tag("file",tab,tags);
						break;
					}
					
					add_bscript(com,tags,bscript);
				}
				break;

			case "Add Pop.": case "Remove Pop.":
				{
					let tab = so.table;
					for(let cl = 0; cl < sp.cla.length; cl++){
						let name = sp.cla[cl].name;
						if(find_in(tab.heading,name) == undefined){
							so.error = true;
							add_warning({mess:so.type+" problem", mess2:"The table does not contain the column '"+name+"'",warn_type:"Init_pop", pa:pa, p:p, i:i});
							return;
						}
					}
				}
				add_flag = true;
				break;
				
			case "Add Ind.":                            // Makes sure Add Ind contains all columns
				{
					let tab = so.table;
					for(let cl = 0; cl < sp.cla.length; cl++){
						let name = sp.cla[cl].name;
						if(find_in(tab.heading,name) == undefined){
							so.error = true;
							add_warning({mess:"Add Ind. problem", mess2:"The table does not contain the column '"+name+"'",warn_type:"Init_pop", pa:pa, p:p, i:i});
							return;
						}
					}
				}
				add_flag = true;
				break;
			
			default:
				add_flag = true;
				break;
			}
			
			if(add_flag == true){	
				add_data_table(so.type,so.table,so,p,i,sim_or_inf,bscript);
			}
		}
	}
}


/// Generates a data-table and command from a data source
function add_data_table(type,tab,so,p,index,sim_or_inf,bscript)
{
	let i = find(convert,"type",type);
	if(i == undefined){
		error("Counld not not create data table:"+type);
		return;
	}
	
	let com = convert[i].command;
	
	if(sim_or_inf == "sim"){
		switch(com){
		case "add-pop-inf": com = "add-pop-sim"; break;
		case "remove-pop-inf": com = "remove-pop-sim"; break;
		case "add-ind-inf": com = "add-ind-sim"; break;
		case "remove-ind-inf": com = "remove-ind-sim"; break;
		case "move-ind-inf": com = "move-ind-sim"; break;
		case "init-pop-inf": com = "init-pop-sim"; break;
		case "test-and-cull": com = "test-and-cull-sim"; break;
		}
	}
	
	if(sim_or_inf == "ppc"){
		switch(com){
		case "add-pop-inf": com = "add-pop-post-sim"; break;
		case "remove-pop-inf": com = "remove-pop-post-sim"; break;
		case "add-ind-inf": com = "add-ind-post-sim"; break;
		case "remove-ind-inf": com = "remove-ind-post-sim"; break;
		case "move-ind-inf": com = "move-ind-post-sim"; break;
		}
	}
	
	let file = com;
	
	let tags=[];
	
	add_tag("name",so.name,tags);
	
	let spec = so.spec;
	
	switch(so.type){                                // Adds information dependent on data source
	case "Init. Pop.":
		if(spec.radio2.value == "Focal") add_tag("focal",spec.focal.te,tags);
		break;
	
	case "Move Ind.":
		add_tag("class",spec.cl_drop.te,tags);
		break;
		
	case "Compartment":
		add_tag("class",spec.cl_drop.te,tags); 
		break;
	
	case "Transition": 
		{
			let sp = model.species[p];
		
			let cl = find(sp.cla,"name",spec.cl_drop.te);

			let trans_fi = output_add_trans_filt(spec.filter,p,cl,so,index);
			if(trans_fi != "") add_tag("trans",trans_fi,tags);
			
			let filt = output_add_comp_filt(spec,sp,cl,p,so,index);
			if(filt != "") add_tag("filter",filt,tags);

			let obsrange = spec.time_radio.value.toLowerCase();
			
			switch(spec.time_radio.value){
			case "All": obsrange = "all"; break;
			case "Spec": obsrange = "specify"; break;
			case "File": obsrange = "file"; break;
			default: error("Option not recognised 1C"); break;
			}
			
			let ty = so.type;
			add_tag("obsrange",obsrange,tags); 
		
			switch(obsrange){
			case "all": break;
			case "specify":
				{
					let start = spec.time_start;
					if(isNaN(start)){	
						add_warning({mess:ty+" data", mess2:"The start time for the "+ty.toLowerCase()+" data must be set to a number", warn_type:"TransDataSpecProb", ind:index});
					}
					
					let end = spec.time_end;
					if(isNaN(end)){	
						add_warning({mess:ty+" data", mess2:"The end time for the "+ty.toLowerCase()+" data must be set to a number", warn_type:"TransDataSpecProb", ind:index});
					}
					
					add_tag("start",start,tags);
					add_tag("end",end,tags);
				}
				break;

			case "file":
				break;

			default: error("Option not recognised 121"+obsrange); break;
			}
		}
		break;
		
	case "Pop. Trans.":
		{
			let sp = model.species[p];
		
			let cl = find(sp.cla,"name",spec.cl_drop.te);

			let trans_fi = output_add_trans_filt(spec.filter,p,cl,so,index);
			if(trans_fi != "") add_tag("trans",trans_fi,tags);
			
			let filt = output_add_comp_filt(spec,sp,cl,p,so,index);
			if(filt != "") add_tag("filter",filt,tags);

			add_obs_model(spec,tags);
		}
		break;

	case "Diag. Test": case "Test-and-cull":
		{
			let err_fl;
			
			let flag = false;
			
			let cb = spec.check_box;
			
			let pcl = get_p_cl_from_claname(cb.name);
			let claa = model.species[pcl.p].cla[pcl.cl];
			
			let st = "";
			for(let i = 0; i < cb.value.length; i++){
				if(cb.value[i].check == true){
					if(flag == true) st += ',';
					st += claa.comp[i].name;
					st += "[Se:"+cb.value[i].Se_eqn.te+"]";
					flag = true;
				}
			}
			add_tag("comp",st,tags);
			
			let Sp = esc(spec.Sp_eqn.te); if(Sp == "") err_fl = "Sp";
			add_tag("Sp",Sp,tags);
			
			let pos = esc(spec.pos_result); if(pos == "") err_fl = "positive value";
			add_tag("pos",pos,tags);
			
			let neg = esc(spec.neg_result); if(neg == "") err_fl = "negative value";
			add_tag("neg",neg,tags);
			
			if(so.type == "Test-and-cull"){
				add_tag("time-gap",spec.time_gap,tags);
				
				let fr = so.frac_obs;
				if(Number(fr) != 1) add_tag("frac",fr,tags);
				
				switch(so.time_radio.value){
				case "Periodic": add_tag("dt",so.time_gen,tags); break;
				case "Specified": add_tag("times",so.time_gen,tags); break;
				default: error("time_radio val"); break;
				}				
			}
			
			if(err_fl != undefined){
				add_warning({mess:"Problem with diagnostic test data", mess2:"The value for '"+eff_fl+"' is not set", warn_type:"SourceProb", p:p, ind:index});
			}
		}
		break;
		
	case "Ind. Eff.":
		add_tag("ie",spec.drop.te,tags);
		break;
		
	case "Ind. Group":
		break;
		
	case "Population": 
		{
			let sp = model.species[p];
		
			let filt = output_add_comp_filt(spec,sp,undefined,p,so,index);
			if(filt != "") add_tag("filter",filt,tags);
			add_obs_model(spec,tags);
		}
		break;
		
	case "Genetic":
		switch(spec.type_radio.value){
		case "matrix": 
			add_tag("type","matrix",tags);
			break;
			
		case "snp":	
			add_tag("type","snp",tags);
			add_tag("root",spec.snp_root,tags);
			break;
		}
		add_tag("mut-rate",spec.mut_rate_eqn.te,tags);
		add_tag("seq-var",spec.seq_var_eqn.te,tags);
		break;

	default: break;	
	}

	if(so.type != "Test-and-cull"){
		add_tag("file",output_table_file(tab),tags);
	}
	
	add_bscript(com,tags,bscript);
}


/// Adds the text for 'filter'
function output_add_comp_filt(spec,sp,cl_sel,p,so,index)
{
	if(spec.filter.cla.length != sp.ncla){
		so.error = true;
		add_warning({mess:so.type, mess2:"The population filter does not agree with the model classifications.", warn_type:"SourceProb", p:p, ind:index});
		return;
	}

	let filt='';
	for(let cl = 0; cl < sp.ncla; cl++){
		if(cl != cl_sel){
			let clz = spec.filter.cla[cl];
			let claa = sp.cla[cl];
			switch(clz.radio.value){
			case "All": break;
			case "Comp": 
				{
					if(filt != "") filt += ',';
					filt += claa.name+'=';
					let flag = false;
					for(let c = 0; c < claa.ncomp; c++){
						if(clz.comp[c].check == true){
							if(flag == true) filt += '|'; flag = true;
							filt += claa.comp[c].name;
						}
					}
				}
				break;
			
			case "ObsMod":
				{
					if(filt != "") filt += ',';
					filt += claa.name+'=';
					let flag = false;
					for(let c = 0; c < claa.ncomp; c++){
						let prob = clz.comp[c].prob_eqn.te.trim();
						if(prob != "0"){
							if(flag == true) filt += '|'; flag = true;
							filt += claa.comp[c].name+":"+prob;
						}
					}
				}
				break;
				
			case "File":
				if(filt != "") filt += ',';
				filt += claa.name+'=file';
				break;
			}
		}
	}

	return filt;
}
				

/// Adds the text for 'name' in pop-trans-data
function output_add_trans_filt(tfilt,p,cl,so,index)
{
	let sp = model.species[p];
	let claa = sp.cla[cl];
	let tra = claa.tra;	
	let trz = tfilt.tra;

	if(tra.length != trz.length){ 
		add_warning({mess:so.type, mess2:"The transition filter does not agree with the model.", warn_type:"SourceProb", p:p, ind:index});
		return;
	}

	if(tfilt.cla.length != sp.ncla){
		add_warning({mess:so.type, mess2:"The transition population filter does not agree with with the model.", warn_type:"SourceProb", p:p, ind:index});
		return;
	}
	
	let name_fi = "";
	
	for(let j = 0; j < trz.length; j++){
		let name = tra[j].name.replace(/→/g,"->")
		
		switch(tfilt.trans_obs_model.value){
		case "on":
			{
				let prob = esc(trz[j].prob_eqn.te);
				if(prob != "0"){
					if(name_fi != "") name_fi += "|";
					name_fi += name+":"+prob;
				}
			}
			break;
			
		case "off":
			if(trz[j].check == true){
				if(name_fi != "") name_fi += "|";
				name_fi += name;
			}
			break;
		
		default: error("Should bnot be def"); break;
		}
	}

	return name_fi;	
}

				
// Ensures that all file names are unique 
function get_unique_file(file,file_list,end)
{	
	file = add_escape_char(file);
	file = file.replace(/ /g,"-");	
	file = file.replace(/\\/g,"");

	let file_not_same;
	let num = 1;
	do{
		file_not_same = file; if(num > 1) file_not_same += "-"+num;
		file_not_same += end;

		if(find(file_list,"file",file_not_same) == undefined) break;
		num++;
	}while(true);
	
	return file_not_same;
}


/// Saves the files to memory for writing later
function write_file_store(data,file,file_list,type)
{
	file_list.push({file:file, data:data, type:type});
}


/// Adds on an observation model to tags
function add_obs_model(spec,tags)
{
	let val;
	
	switch(spec.obs_error_dist.value){
	case "Normal":		
		switch(spec.obs_error.value){
		case "percent": val = "normal:"+spec.percent+"%"; break;
		case "sd": val = "normal:"+spec.sd; break;
		case "file": val = "normal:file"; break;
		default: error("Option not recognised 123"); break;
		}
		break;
	
	case "Poisson":		
		val = "poisson";
		break;
		
	case "Negative binomial":	 
		switch(spec.obs_error_p.value){
		case "p": val = "neg-binomial:"+spec.p;	break;
		case "file": val = "neg-binomial:file";
		default: error("Option not recognised 123"); break;
		}
		break;
		
	default: error("obs dist not found"); break;
	}
	
	add_tag("error",val,tags);
}


/// Outputs a table 
function output_table_file(tab)
{
	var st = '';

	for(let i = 0; i < tab.heading.length; i++){
		if(i != 0) st += ',';
		st += '"'+tab.heading[i]+'"';
	}
	st += endl;

	for(let r = 0; r < tab.ele.length; r++){
		let el = tab.ele[r];
		if(el.length != tab.heading.length){
			error("Could not write file...columns do not match"); return;
		}

		for(let i = 0; i < el.length; i++){
			if(i != 0) st += ',';
			
			let val = el[i];
			if(isNaN(val)) st += '"'+val+'"';
			else st += val;
		}
		st += endl;
	}
	
	return st;
}


/// Outputs a simple table 
function output_table_simp_file(head,row) 
{
	var st = '';
	for(let i = 0; i < head.length; i++){
		if(i != 0) st += ',';
		st += '"'+head[i]+'"';
	}
	st += endl;

	for(let r = 0; r < row.length; r++){
		let ele = row[r];
		if(ele.length != head.length){
			error("Could not write file...columns to not match"); return;
		}

		for(let i = 0; i < ele.length; i++){
			if(i != 0) st += ',';
			let val = ele[i];
			if(isNaN(val)) st += '"'+val+'"';
			else st += val;
		}
		st += endl;
	}
	
	return st;
}
							

/// Checks the model is correctly specified
function output_check(save_type)
{
	// Checks to make sure initial population is set correctly
	let spec;
	if(save_type == "ppc"){
		let pf = model.inf_res.plot_filter;
		if(pf){
			spec = pf.species;
		}
	}
	else{
		spec = model.species;
	}
	
	if(!spec) return;
	
	for(let p = 0; p < spec.length; p++){
		let sp = spec[p];
		
		let source, pa;
		
		switch(save_type){
		case "sim": pa = "Simulation"; source = sp.sim_source; break;
		case "inf": pa = "Inference"; source = sp.inf_source; break;
		case "ppc": pa = "Modification"; source = sp.ppc_source; break;
		default: error("Option not recognised 125"); break;
		}
		
		let ninitpop = 0, ninitpop_fix = 0, naddind = 0, naddpop = 0;
		//if(so.init_pop_type == INIT_POP_FIXED) ninitpop_fix++; 

		
		for(let i = 0; i < source.length; i++){
			let so = source[i];
	
			switch(so.type){
			case "Init. Pop.": ninitpop++; if(so.spec.radio_dist.value == 'Fixed') ninitpop_fix++; break;
			case "Add Ind.": naddind++; break;
			case "Add Pop.": naddpop++; break;
			}
		}

		switch(save_type){
		case "sim":
			if(ninitpop > 0 && sp.fix_eff.length > 0){
				add_warning({mess:"Cannot add individuals with fixed effect", mess2:"The individuals from 'Init. Pop.' cannot be added because of individual fixed effects in the model", warn_type:"SimPopulationProb", p:p});
			}
			
			if(ninitpop > 1){
				add_warning({mess:"Initial population", mess2:"Only one initial population should be set for species '"+sp.name+"'", warn_type:"SimPopulationProb", p:p});
			}
			break;
		
		case "inf":
			if(ninitpop_fix > 0 && sp.type == "Individual"){
				add_warning({mess:"Initial population", mess2:"Fixed initial populations cannot be used with individual-based models for species '"+sp.name+"'. Either a population distribution is used, or the initial state of individuals is specified using 'add-ind' (which can be used to capture the any initial uncertainty).", warn_type:"InfPopulationProb", p:p});
			}
		
			if(ninitpop > 1){
				add_warning({mess:"Initial population", mess2:"Only one initial population should be set for species '"+sp.name+"'", warn_type:"InfPopulationProb", p:p});
			}

			switch(sp.type){
			case "Population": case "Deterministic":
				if(ninitpop == 0 && naddpop == 0 && model.source_exist(p) == false){
					add_warning({mess:"No individuals", mess2:"Either 'Init. Pop.' or 'Add Pop.' must be set for species '"+sp.name+"'", warn_type:"InfPopulationProb", p:p});
				}
				break;
			
			case "Individual":
				if(ninitpop == 0 && naddind == 0 && model.source_exist(p) == false){
					add_warning({mess:"No individuals", mess2:"Either 'Init. Pop.' or 'Add Ind' must be set for species '"+sp.name+"'", warn_type:"InfPopulationProb", p:p});
				}
				break;
			
			default: error("op prob"); break;	
			}
			
			if(ninitpop > 0 && sp.fix_eff.length > 0){
				add_warning({mess:"Cannot add individuals with fixed effect", mess2:"The individuals from 'Init. Pop.' cannot be added because of individual fixed effects in the model", warn_type:"InfPopulationProb", p:p});
			}
			break;

		case "ppc":
			break;

		default: error("Option not recognised 126"); break;
		}
	}
}


/// Checks compartmental structure (looks for repeated names etc...) 
function check_comp_structure()
{
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			
			for(let p2 = p; p2 < model.species.length; p2++){
				let sp2 = model.species[p2];
				let cl2min = 0; if(p == p2) cl2min = cl+1;
				for(let cl2 = cl2min; cl2 < sp2.ncla; cl2++){
					let claa2 = sp2.cla[cl2];
					
					if(!(p == p2 && cl == cl2)){
						if(claa.name == claa2.name){
							if(p == p2){
								add_warning({mess:"Repeated classification", mess2:"The classification '"+claa.name+"' is repeated twice in '"+sp.name+"'", warn_type:"ModelClass", p:p, cl:cl});
								return;
							}
							else{										
								if(claa.index != claa2.index){
									add_warning({mess:"Different indices", mess2:"The classification '"+claa.name+"' on species '"+sp.name+"' and '"+sp2.name+"' have different indices.", warn_type:"ModelClass", p:p, cl:cl});
									return;
								}
								else{
									let vec1=[], vec2=[];
									for(let c = 0; c < claa.ncomp; c++) vec1.push(claa.comp[c].name);
									for(let c = 0; c < claa2.ncomp; c++) vec2.push(claa2.comp[c].name);
									
									if(equal_vec(vec1,vec2) == false){
										add_warning({mess:"Compartments different", mess2:"For classification '"+claa.name+"' the compartments do not agree between species '"+sp.name+"' and '"+sp2.name+"'", warn_type:"ModelClass", p:p, cl:cl});
										return;
									}
								}
							}
						}
						else{
							if(claa.index == claa2.index){
								if(p == p2){
									add_warning({mess:"Identical indices", mess2:"In species '"+sp.name+"' the classifications '"+claa.name+"' and '"+claa2.name+"' should not share the same index '"+claa.index+"'", warn_type:"ModelClass", p:p, cl:cl});
									return;
								}
								else{
									add_warning({mess:"Identical indices", mess2:"Species '"+sp.name+" classifications '"+claa.name+"' and species '"+sp2.name+"' classification '"+claa2.name+"' should not share the same index '"+claa.index+"'", warn_type:"ModelClass", p:p, cl:cl});
									return;
								}
							}
						}
					}
				}
			}
		}
	}
}


/// Removes specified commands from a file
function remove_command(bscript,com) 
{
	let pos = com.split(",");
	
	let i = 0;
	while(i < bscript.length){
		let com = bscript[i];
		if(find_in(pos,com.type) != undefined || com.type == "comment" && find_in(pos,com.te)){
			bscript.splice(i,1);
		}
		else i++;
	}
}


/// Extracts a particular command type from the 
function extract_command(bscript,bscript_add,com) // TO DO is this needed?
{
	let pos = com.split(",");
	
	for(let i = 0; i < bscript.length; i++){
		let com = bscript[i];
		if(find_in(pos,com.type) != undefined) bscript_add.push(com);
	}
}	

 
/// Creates a file to start a posterior simulation
function create_ppc_file(info,one_file)
{
	let bscript = copy(inf_result.bscript);
	
	if(info.use_inf == undefined) error("Use inf should be defined");
	
	if(info.use_inf == false){
		let bscript_add=[];
		mini_banner("OUTPUT INFERENCE",bscript_add);
		
		extract_command(bscript,bscript_add,"param-inf,state-inf,param-stats-inf,diagnostics-inf,generation-inf,trans-diag-inf,proposal-inf,warning-inf");
		
		create_output_file(info.save_type,one_file,info.map_store,bscript_add);
		return;
	}		

	remove_ppc_command(bscript,true);
	
	// Alters post-sim command
	let i = find(bscript,"type","post-sim");
	if(i == undefined) i = find(bscript,"type","posterior-simulation");
	
	let tags = get_ppc_tags(info.save_type);
	
	if(i == undefined){
		add_bscript("post-sim",tags,bscript);
	}
	else{
		bscript[i].tags = tags;
	}
	
	for(let p = 0; p < inf_result.species.length; p++){
		let sp = inf_result.species[p];
		if(inf_result.species.length > 1){
			let tags=[];
			add_tag("species","sp.name",tags);
			add_bscript("set",tags,bscript);
		}
		create_output_sim_inf_source(p,"ppc",bscript);
	}
	
	
	let pfac_list=[];
	for(let th = 0; th < model.param.length; th++){
		let par = model.param[th];
		if(par.type == "param factor") pfac_list.push(th);
	}
	
	if(pfac_list.length > 0){
		mini_banner("PARAMETER MULTIPLIERS",bscript);
		for(let k = 0; k < pfac_list.length; k++){
			let par = model.param[pfac_list[k]];
			if(par.type == "param factor"){
				output_param(par,"ppc",bscript);
			}
		}
	}
 
	if(model.warn.length > 0) err_warning();
	
	let file_list=[];
	let te = generate_text_from_bscript(bscript,file_list,one_file);

	write_file_store(te,"bicifile",file_list,"bicifile");

	post({save_type:"ppc", file_list:file_list, param_factor:strip_heavy(model.param_factor), param:strip_heavy(model.param), species:strip_heavy(model.species)});
}


/// Removes all commands associated with ppc
function remove_ppc_command(bscript,rem_data)
{
	
	remove_command(bscript,"param-post-sim,state-post-sim,warning-post-sim,# OUTPUT POSTERIOR SIMULATION"); 
	
	if(rem_data){
		remove_command(bscript,"param-mult,add-pop-post-sim,remove-pop-post-sim,add-ind-post-sim,remove-ind-post-sim,move-ind-post-sim,# MODIFICATION DATA,# PARAMETER MULTIPLIERS");
	}
}


/// Creates a file to start extension
function create_ext_file(one_file)
{
	let bscript = copy(inf_result.bscript);
	
	remove_ppc_command(bscript,true);
	
	let file_list=[];
	let te = generate_text_from_bscript(bscript,file_list,one_file);
	
	write_file_store(te,"bicifile",file_list,"bicifile");

	post({save_type:"ext", file_list:file_list, param_factor:strip_heavy(model.param_factor), param:strip_heavy(model.param), species:strip_heavy(model.species)});
}

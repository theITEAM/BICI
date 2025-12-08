"use strict";
// Functions which generate output .bici file

/// Creates an output file
// save_type has the values:
// sim - Saves a simulation (this performs additional simulation checks)
// inf - Saves an inference (this performs additional inference checks)
// save - Saves files (somes things are not checked)

// one_file deterines if a a single file is generated or a file and directory are made

/// Creates the ouytput file from the model
function create_output_file(save_type,one_file,map_store)
{
	let file_list=[];  // Used to make sure that files to not share the same name
	
	percent(0);

	update_model();

	percent(30);
	
	let te = create_output_details(save_type,file_list,one_file);

	percent(40);

	check_data_valid_all(save_type);
	
	check_data_time(save_type);
	
	check_param_valid(save_type);
	
	check_integral_bounds(save_type);
	
	percent(45);
	
	for(let p = 0; p < model.species.length; p++){
		te += create_output_species(p,file_list,save_type,map_store,one_file);
	}

	percent(50);

	te += create_output_param(save_type,file_list,one_file);
	
	percent(60);
	
	te += create_output_derived();
	
	te += endl;
		
	if(model.warn.length > 0) err_warning();
	
	percent(70);
	
	//check_comp_structure();

	percent(80);

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
	case "StartClusterSave": case "StartClusterExport": case "StartPPC":
		full_warn = true;
		break;
	}

	let wa = {type:"Model warning", command_type:input.type, species:strip_heavy(model.species), warn:model.warn, full_warn:full_warn};
	if(try_on == false) prr(wa);
	else throw(wa);
}


/// Creates a banner for the output file
function banner(name)
{
	let num = 58;
	let num2 = num - name.length;
	let num3 = Math.floor(num2/2);
	
	let te = endl;
	for(let i = 0; i < num3; i++) te += "#";
	te += " "+name+" ";
	for(let i = 0; i < num2-num3; i++) te += "#";
	te += endl+endl;
	
	return te;
}


/// Creates a banner for the output file
function mini_banner(name)
{
	return "# "+name+endl+endl;
}


/// Generates details
function create_output_details(save_type,file_list,one_file)
{
	let te = "";

	if(!one_file){
		te += banner("DATA DIRECTORY");
		te += 'data-dir folder="."'+endl;
		te += endl;
	}
	
	te += banner("DESCRIPTION");
	
	let desc = model.description.te;
	let j = desc.length-1;
	while(j >= 0 && (desc.substr(j,1) == "\n" || desc.substr(j,1) == "\r")) j--;
	desc = desc.substr(0,j+1);
	desc += endl;
	
	let file = get_unique_file("description",file_list,'.txt');
	write_file_store(desc,file,file_list);
	
	//output_table_simp_file(file,head,row,file_list);
	if(one_file) file = get_one_file(file_list);
							
	//te += 'description text="[['+endl+desc+endl+']]"'+endl;
	te += 'description text="'+file+'"'+endl;
	te += endl;
	
	te += banner("DETAILS");
	
	te += create_output_siminf(save_type);

	return te;
}
	
	
/// Outputs all the results for a given species
function create_output_species(p,file_list,save_type,map_store,one_file)
{
	let sp = model.species[p];
		
	let te = banner("DEFINE MODEL AND DATA FOR SPECIES "+sp.name.toUpperCase());
	
	let ty = sp.type.toLowerCase();
	
	te += 'species name="'+sp.name+'" type="'+ty+'"';

	if(ty == "individual"){
		if(sp.trans_tree.check == true){
			te += ' trans-tree="';
			if(sp.trans_tree.check == true) te += 'on';
			else te += 'off';
			te += '"';
		}
	}
	
	te += endl+endl;
	
	te += create_output_compartments(p,file_list,map_store,one_file);
	
	te += create_output_ind_eff(p,file_list,one_file);

	te += create_output_fix_eff(p,file_list,save_type,one_file);

	te += create_output_sim_inf_source(p,"sim",file_list,one_file);
	if(save_type == "sim") output_check("sim");	

	te += create_output_sim_inf_source(p,"inf",file_list,one_file);
	if(save_type == "inf") output_check("inf");	
	
	if(model.inf_res.on == true){
		te += create_output_sim_inf_source(p,"ppc",file_list,one_file);
		if(save_type == "ppc") output_check("ppc");	
	}

	return te;
}


/// Creates output for all the compartments and transitions in the model
function create_output_compartments(p,file_list,map_store,one_file)
{
	let te = mini_banner("SPECIES MODEL");
		
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
		
		te += 'class name="'+claa.name+'"';
		
		if(claa.default_map == true) te += ' map="load"'; 
		
		let clp = clone_p(p,cl);

		if(clp != undefined) te += ' clone="'+model.species[clp].name+'"'+endl+endl;
		else{
			te += ' index="'+claa.index+'"';
			if(cam.coord != COORD_DEFAULT) te += ' coord="'+cam.coord+'"';
			
			te += endl+endl;
		
			if(one_file || cam.coord == "latlng"){
				te += 'view '+output_coords(cam.x,cam.y,cam);
				
				if(cam.grid == "on"){
					te += ' grid="on"'; 
				}
				
				if(cam.slider.value != 0){
					te += ' comp-scale='+Math.exp(cam.slider.value).toPrecision(2);
				}
				
				te += ' scale='+precision(cam.scale);
				te += endl+endl;
			}
			
			let store=[];
			for(let c = 0; c < claa.ncomp; c++){
				let labs=[];
				
				let co = claa.comp[c];
				
				labs.push({tag:"name", val:co.name});
				
				labs.push({tag:"color", val:co.col});
				
				labs.push({tag:"fix", val:co.fixed.check});
			
				if(sp.type == "Individual" && sp.trans_tree.check == true){
					if(sp.infection_cl.te == sp.cla[cl].name && co.infected.check != undefined){
						labs.push({tag:"infected", val:co.infected.check});
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
						let file = get_unique_file(co.name+"_boundary",file_list,'.geojson');
					
						let k = find(map_store,"name",co.map_ref); 
						if(k == undefined) error("Cannot find feature");
	
						let map_da = { crs:{},features:[],type:"FeatureCollection"};
						map_da.features.push(generate_JSON_feature(co.name,map_store[k].feature));
						
						write_file_store(JSON.stringify(map_da)+"\n",file,file_list);
						
						if(one_file) file = get_one_file(file_list);
						
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
				te += output_command_list("comp",claa.name,store,file_list,one_file);
				te += endl;	
			}
		}

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
		
			let tra_type = tr.type.toLowerCase();
			let nm_flag = false;
			let err_fl;
			
			let te;
			switch(tra_type){
			case "exp(rate)": 
				{
					let rate = esc(tr.value.rate_eqn.te); if(rate == "") err_fl = "rate";
					te = "exp(rate:"+rate+")";
				}
				break;
				
			case "exp(mean)": 
				{
					let mean = esc(tr.value.mean_eqn.te); if(mean == "") err_fl = "mean";
					te = "exp(mean:"+mean+")";
				}
				break;
				
			case "gamma":
				{
					let mean = esc(tr.value.mean_eqn.te); if(mean == "") err_fl = "mean";
					let cv = esc(tr.value.cv_eqn.te); if(cv == "") err_fl = "cv";
					te = "gamma(mean:"+mean+", cv:"+cv+")";
					nm_flag = true;
				}
				break;
				
			case "erlang":
				{
					let mean = esc(tr.value.mean_eqn.te); if(mean == "") err_fl = "mean";
					let sh = esc(tr.value.shape_erlang.te); if(sh == "") err_fl = "shape";
					te = "erlang(mean:"+mean+", shape:"+sh+")";
				}
				break;
				
			case "log-normal": 
				{
					let mean = esc(tr.value.mean_eqn.te); if(mean == "") err_fl = "mean";
					let cv = esc(tr.value.cv_eqn.te); if(cv == "") err_fl = "cv";
					te = "log-normal(mean:"+mean+", cv:"+cv+")";
					nm_flag = true;
				}
				break;
				
			case "weibull":
				{
					let sc = esc(tr.value.scale_eqn.te); if(sc == "") err_fl = "scale";
					let sh = esc(tr.value.shape_eqn.te); if(sh == "") err_fl = "shape";
					te = "weibull(scale:"+sc+", shape:"+sh+")";
					nm_flag = true;
				}
				break;

			case "period": 
				{
					let time = esc(tr.value.mean_eqn.te); if(time == "") err_fl = "time";
					te = "period(time:"+time+")";
					nm_flag = true;
				}
				break;
			
			default: error("Option not recognised 114"); break;
			}
			
			labs.push({tag:"value", val:te});
			
			if(err_fl != undefined){
				add_warning({mess:"Transition problem.", mess2:"The value for '"+err_fl+"' is not specified", p:p, cl:cl, i:j, warn_type:"TransMistake"});
			}
			
			if(sp.type == "Population" && nm_flag){ 		
				add_warning({mess:"Transition problem.", mess2:"Transition '"+tr.name+"' cannot have a '"+tra_type+"' distribution because '"+sp.name+"' uses a population-based model.", p:p, cl:cl, i:j, warn_type:"TransMistake"});
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
			te += output_command_list("trans",claa.name,store,file_list,one_file);
		
			te += endl;
		}
		
		if(claa.annotation.length > 0){
			for(let i = 0; i < claa.annotation.length; i++){
				let anno = claa.annotation[i];
				switch(anno.type){
				case "box":
					te += 'box text="'+anno.te+'" comps="';
					for(let k = 0; k < anno.comps.length; k++){
						if(k != 0) te += ",";
						te += anno.comps[k];
					}						
					te += '"'
					if(anno.color != annotation_col_default) te += ' color="'+anno.color+'"';
					if(anno.size != size_annotation_default) te += ' text-size='+anno.size;
					te += endl;
					break;
					
				case "text":
					te += 'label text="'+anno.te+'" '+output_coords(anno.x,anno.y,cam);
					if(anno.color != annotation_col_default) te += ' color="'+anno.color+'"';
					if(anno.size != size_annotation_default) te += ' text-size='+anno.size;
					te += endl;
					break;
					
				case "map":
					if(anno.default_map != true){
						let file = get_unique_file("map",file_list,'.geojson');
						let k = find(map_store,"name",anno.map_ref);
						if(k == undefined) error("map undefined");
						else{
							let map_da = { crs:{},features:[],type:"FeatureCollection"};
							for(let j = 0; j < map_store[k].feature.length; j++){
								map_da.features.push(generate_JSON_feature("Feature"+k,map_store[k].feature[j]));
							}
						
							write_file_store(JSON.stringify(map_da)+"\n",file,file_list);
							if(one_file) file = get_one_file(file_list);
							te += 'map file="' +file+'"' + endl;
						}
					}
					break;
					
				default:
					break;
				}
			}
			te += endl;
		}
	}

	return te;
}


/// Gets string for coordinates
function output_coords(x,y,cam)
{
	let te = "";
	switch(cam.coord){
	case "cartesian": 
		te += 'x='+precision(x)+' y='+precision(y); 
		break;
		
	case "latlng": 
		let p = transform_latlng_inv(x,y);
		te += 'lat='+precision(p.lat)+' lng='+precision(p.lng); 
		break;

	default: error("Option not recognised 110"); break;
	}
	
	return te;
}

	
/// Outputs a list of commands either seperately or in a file (for comps-all and trans-all)
function output_command_list(type,claa_name,store,file_list,one_file)
{
	let te = "";
	
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
		
		te += type+'-all';
		
		// Looks to remove tags which are all the same
		let i = 0;
		while(i < N){
			let val = row[0][i];
			
			let j = 1; while(j < store.length && row[j][i] == row[0][i]) j++;
			
			if(j == store.length){
				te += " "+head[i]+"=";
				if(is_nan(val)) te += '"'+val+'"';
				else te += val;
				
				head.splice(i,1);
				for(let j = 0; j < store.length; j++){
					row[j].splice(i,1);
				}
				N--;
			}
			else i++;
		}
		
		let file = get_unique_file(type+"-"+claa_name,file_list,'.csv');

		output_table_simp_file(file,head,row,file_list);
		if(one_file) file = get_one_file(file_list);
		
		te += ' file="'+file+'"'+endl;
	}	
	else{		
		for(let j = 0; j < store.length; j++){
			let labs = store[j];
			te += type;
			for(let k = 0; k < labs.length; k++){
				if(!(labs[k].tag == "fix" && labs[k].val == false)){
					te += " "+labs[k].tag+"=";
					let val = labs[k].val;
					if(is_nan(val)) te += '"'+val+'"';
					else te += val;
				}
			}
			te += endl;
		}
	}
	
	return te;
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
function create_output_param(save_type,file_list,one_file)
{
	if(model.param.length == 0) return "";
	
	let te = banner("PARAMETERS");
	
	for(let th = 0; th < model.param.length; th++){
		te += output_param(model.param[th],save_type,file_list,one_file);
	}
	
	te += endl;

	return te;
}


/// Outputs a single parameter
function output_param(par,save_type,file_list,one_file)
{
	let te = "";
	let dist_done = false;
	
	if(par.type != "derive_param"){	 
		let num_warn = model.warn.length;
	
		let te1 = "param";
		if(par.type == "param factor"){
			let pfac = model.param_factor;
		
			let k = 0, kmax = pfac.length; 

			while(k < kmax && pfac[k].f_param.full_name != par.full_name) k++;
			if(k == kmax) alert_help("Problem loading param factor");
	
			te1 += '-mult name="'+remove_eq_quote(pfac[k].param.full_name)+'"';
		}
		else{
			te1 += ' name="'+remove_eq_quote(par.full_name)+'"';
		}
		
		let te2 = "";
		let display = false;
		
		if(par.dist_mat || par.iden_mat || par.den_vec){
			display = true;
		}
		else{
			// Constants and simulation parameters
			if(par.variety == "reparam" && par.reparam_eqn_on){
				if(par.reparam_eqn.trim() != ""){
					te1 += ' reparam="'+par.reparam_eqn+'"';
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
						te2 += output_add_prior_distribution(save_type,par,"dist",file_list,one_file);
						dist_done = true;
					}
					else{
						let exp = "value";
						if(par.variety == "const"){
							exp = 'constant';
						}
						else{
							if(par.variety == "reparam") exp = 'reparam';
						}
									
						te2 += ' '+exp+'=';
					
						let file;				
						if(par.dep.length > 0){
							file = get_unique_file(exp+"-"+par.name,file_list,'.csv');
						}
						
						if(file != undefined){
							let err = check_param_value("Set Param",par,par.value);
							if(err){
								let wt = "MissingSimValue"; if(par.variety == "reparam") wt = "ReparamValue";
								
								add_warning({mess:"Problem with "+par.full_name+" value", mess2:err, warn_type:wt, name:par.name});
								return;
							}
						
							file = output_value_table(par,par.value,"Value",file,file_list,one_file,"value",save_type);
							te2 += '"'+file+'"';
						}
						else{
							te2 += JSON.stringify(par.value); 
						}
					}	
				}
			}
			
			if(par.variety != "const" && par.variety != "reparam"){
				switch(par.variety){
				case "normal":
					te2 += output_add_prior_distribution(save_type,par,"prior",file_list,one_file);
					break;
					
				case "dist":
					if(dist_done == false){
						te2 += output_add_prior_distribution(save_type,par,"dist",file_list,one_file);
					}
					break;
					
				default: error("Option not recognised 1B"+par.variety); break;
				}
			}
			
			if(te2 != "") display = true;
			
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
					
				te2 += ' knot-times="'+stringify(par.spline.knot)+'"';
				
				let sm = par.spline.smooth;
				if(sm.check == true){
					te2 += ' smooth="'+sm.type.value.toLowerCase()+'('+sm.value+')"';
				}
				
				let sel = par.spline.spline_radio.value;

				if(sel != "Linear"){
					switch(sel){
					case "Linear": te2 += ' spline-type="linear"'; break;
					case "Square": te2 += ' spline-type="square"'; break;
					case "Cubic +ve": te2 += ' spline-type="cubic +ve"'; break;
					case "Cubic": te2 += ' spline-type="cubic"'; break;
					}
				}
			}
		}
		
		if(par.factor){
			te1 += ' factor="true"';
			
			if(par.factor_weight_on.check == true){
				let file = get_unique_file("weight-"+par.name,file_list,'.csv');
				file = output_value_table(par,par.factor_weight,"Value",file,file_list,one_file,"value",save_type);
				te1 += ' factor-weight="'+file+'"';
			}
		}
		
		
		if(par.sim_sample && par.sim_sample.check == false){
			te1 += ' sim-sample="false"';
		}
				
		if(!par.dist_mat && !par.iden_mat && !par.den_vec){
			te += te1+te2+endl+endl;
		}
		
		if(display == false && (save_type == "sim" || save_type == "inf")){
			if(num_warn == model.warn.length){
				add_warning({mess:"Parameter "+par.full_name+" could not be exported.", mess2:"Unknown error", warn_type:"", name:par.name});
			}
		}
	}
	
	return te;
}
	
	
/// Outputs the value of a table
function output_value_table(par,value,head_col,file,file_list,one_file,key,save_type)
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
		}
		
		add_warning({mess:"Missing "+key, mess2:"The "+key+" for "+par.full_name+" must be set.", warn_type:wt, name:par.name});
		return;
	}
						
	let data = '';
	if(is_matrix(par)){   // Matrix format
		for(let i = 0; i < list[1].length; i++){
			if(i != 0) data += ",";
			data += list[1][i];
		}
		data += endl;
		
		for(let j = 0; j < list[0].length; j++){
			for(let i = 0; i < list[1].length; i++){
				if(i != 0) data += ",";
				
				let index = []; index[0] = j; index[1] = i;
				let el = get_element(value,index);
				if(el == undefined){ param_undefined(par,index); return;}
					
				data += '"'+get_op_val(el,head_col,par)+'"';
			}
			data += endl;
		}
	}
	else{
		for(let d = 0; d < par.dep.length; d++){	
			data += '"'+par.dep[d]+'",';
		}
		data += head_col + endl;

		let co_list = generate_co_list(par.list);
		for(let i = 0; i < co_list.length; i++){
			let comb = co_list[i];
			let el = get_element(value,comb.index);
			if(el == undefined){ param_undefined(par,comb.index); return;}
			
			if(el != 0){	
				for(let d = 0; d < par.dep.length; d++){
					data += '"'+list[d][comb.index[d]]+'",';
				}
				data += '"'+get_op_val(el,head_col,par)+'"'+endl;
			}
		}
	}

	write_file_store(data,file,file_list);
	if(one_file) file = get_one_file(file_list);
	
	return file;
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
function output_add_prior_distribution(save_type,par,variety,file_list,one_file)
{
	let st="";

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
			let file = get_unique_file(key+"-"+par.name,file_list,'.csv');
			
			let head_col;
			switch(variety){
			case "dist": st += ' dist-split'; head_col = 'Dist'; break;
			case "prior": st += ' prior-split'; head_col = 'Prior'; break;
			default: error("option prob"); break;
			}
			
			file = output_value_table(par,par.prior_split,head_col,file,file_list,one_file,key,save_type);
			if(file == undefined) return "";
			st += '="'+file+'"';
		}
	}
	else{
		let pri = par.prior;
		
		let prior_string = get_prior_output(save_type,par.prior,variety,par);
		if(prior_string == "") return "";
	
		switch(variety){ 
		case "dist": st += ' dist'; break;
		case "prior": st += ' prior'; break;
		default: error("option prob"); break;
		}
		
		st += '='+prior_string;
	}
		
	return st;
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
	
	return '"'+st+'"';
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
function create_output_derived()
{
	let derive = model.derive;

	if(derive.length == 0) return "";

	let te = mini_banner("DERIVED");

	for(let i = 0; i < derive.length; i++){
		let der = derive[i];
		
		if(der.eqn1.param.length != 1) error("problem with derived");
		else{
			let name = remove_eq_quote(der.eqn1.param[0].full_name);	
			te += 'derived name="'+name+'" eqn="'+der.eqn2.te+'"'+endl;
		}
	}
	te += endl;

	return te;
}


/// Create output of individual effects
function create_output_ind_eff(p,file_list,one_file)
{
	let iegs = model.species[p].ind_eff_group;
		
	if(iegs.length == 0) return "";
	
	let te = mini_banner("INDIVIDUAL EFFECTS");

	for(let i = 0; i < iegs.length; i++){
		let ieg = iegs[i];
		let ie = "";
		for(let j = 0; j < ieg.ie_list.length; j++){
			if(j > 0) ie += ",";
			ie += ieg.ie_list[j].name;
		}
			
		te += 'ind-effect name="'+ieg.name+'" ie="'+ie+'"';
		
		let mat = ieg.A_matrix;
		
		if(mat.check == true){
			if(mat.loaded == false){
				add_warning({mess:"Error in individual effects", mess2:"The A matrix for individual effects '"+name+"' must be loaded", warn_type:"A matrix"});
			}
			else{
				switch(mat.type){
				case "pedigree":
					{
						let file = get_unique_file("pedigree-"+name,file_list,'.csv');
						
						let data = 'ID,sire,dam'+endl;
						
						for(let c = 0; c < mat.ind_list.length; c++){ 
							data += '"'+mat.ind_list[c]+'"'+",";
							data += '"'+mat.sire_list[c]+'"'+",";
							data += '"'+mat.dam_list[c]+'"'+endl;
						}
							
						write_file_store(data,file,file_list);
						if(one_file) file = get_one_file(file_list);
						te += ' pedigree="'+file+'"';
					}
					break;
					
				case "A": 
					{
						let file = get_unique_file("A-matrix-"+name,file_list,'.csv');
						
						if(true){                                   // Sparse representation
							let file_ind = get_unique_file("ind-list-"+name,file_list,'.csv');
						
							let data_ind = 'Individual'+endl;
							for(let c = 0; c < mat.ind_list.length; c++){ 
								data_ind += '"'+mat.ind_list[c]+'"'+endl;
							}
						
							write_file_store(data_ind,file_ind,file_list);
							if(one_file) file_ind = get_one_file(file_list);
							te += ' ind-list="'+file_ind+'"';
							
							let data = 'j,i,value'+endl;
							for(let r = 0; r <  mat.ind_list.length; r++){
								for(let c = 0; c <= r; c++){
									let val = mat.A_value[r][c];
									
									if(val != 0) data += r+","+c+","+val+endl;
								}
							}
				
							write_file_store(data,file,file_list);
							if(one_file) file = get_one_file(file_list);
							te += ' A-sparse="'+file+'"';
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
				
							write_file_store(data,file,file_list);
							if(one_file) file = get_one_file(file_list);
							te += ' A="'+file+'"';
						}
					}
					break;
					
				case "Ainv": 
					{
						let file = get_unique_file("Ainv-matrix-"+name,file_list,'.csv');
					
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
				
						write_file_store(data,file,file_list);
						if(one_file) file = get_one_file(file_list);
						te += ' Ainv="'+file+'"';
					}
					break;
				}
			}
		}		
		te += endl;
	}

	te += endl;

	return te;
}


/// Create output of fixed effects
function create_output_fix_eff(p,file_list,save_type,one_file)
{
	let fegs = model.species[p].fix_eff;
	
	if(fegs.length == 0) return "";
	
	let te = mini_banner("FIXED EFFECTS");

	for(let i = 0; i < fegs.length; i++){
		let feg = fegs[i];
		let name = feg.name;
		
		te += 'fixed-effect name="'+feg.name+'"';
		
		let vec = feg.X_vector;
		if(vec.loaded == false){
			if(save_type == "sim" || save_type == "inf"){
				add_warning({mess:"Error in fixed effects", mess2:"The vector <e>X^"+name+"</e> for fixed effect '"+name+"' must be loaded", warn_type:"X vector"});
			}
		}
		else{
			let file = get_unique_file("X-vector-"+name,file_list,'.csv');

			let data = 'ID,value'+endl;
			for(let c = 0; c < vec.ind_list.length; c++){ 
				data += '"'+vec.ind_list[c]+'",'+vec.X_value[c]+endl;
			}
			write_file_store(data,file,file_list);
			if(one_file) file = get_one_file(file_list);
			te += ' X="'+file+'"';
		}
		te += endl;
	}

	te += endl;

	return te;
}

	
/// Creates information about simulation or inference
function create_output_siminf(save_type)
{
	let te = "";

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
		
		te += "simulation";
		if(details.t_start != "") te +=" start="+details.t_start;
		if(details.t_end != "") te += " end="+details.t_end;
		
		if(details.number != ""){
			let num = Number(details.number);
			if(num != SIM_NUM_DEFAULT) te += " number="+num;
		}
		
		if(details.timestep != "") te += " timestep="+details.timestep; 
		
		if(details.seed_on.value == "Yes") te += " seed="+Number(details.seed);
		 
		if(details.indmax != INDMAX_DEFAULT){
			te += " ind-max="+details.indmax;
		}
		
		if(details.param_output_max != PARAM_OUTPUT_MAX_DEFAULT){
			te += " param-output-max="+details.param_output_max;
		}
		
		te += endl;
		te += endl;
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
		
		te += "inference";
		if(details.t_start != "") te +=" start="+details.t_start;
		if(details.t_end != "") te += " end="+details.t_end;
		if(details.timestep != "") te += " timestep="+details.timestep; 
		
		if(details.seed_on.value == "Yes") te += " seed="+Number(details.seed);
		
		let alg = details.algorithm.value;
		
		switch(alg){
		case "DA-MCMC": case "PAS-MCMC":
			{
				let num_samp = Number(details.sample);
				if(num_samp != MCMC_SAMPLE_DEFAULT){
					te += " update="+num_samp;
				}
				
				let num_param = Number(details.output_param);
				if(num_param != MCMC_OP_PARAM_DEFAULT){
					te += " param-output="+num_param;
				}
				
				let num_state = Number(details.output_state);
				if(num_state != MCMC_OP_STATE_DEFAULT){
					te += " state-output="+num_state;
				}
				
				if(alg == "DA-MCMC"){
					let cha_per_core = Number(details.cha_per_core);
					if(cha_per_core != MCMC_CHAIN_PER_CORE_DEFAULT){
						te += " chain-per-core="+cha_per_core;
					}
					
					if(Number(details.nchain)%cha_per_core != 0){
						add_warning({mess:"Error with chain number", mess2:"The number of chains '"+details.nchain+"' must be a multiple of the number of chains per core '"+cha_per_core+"'", warn_type:"InfFurtherDetails"});
					}
				}
				
				if(alg == "PAS-MCMC"){
					let part_per_core = Number(details.part_per_core);
					if(part_per_core != PAS_PART_PER_CORE_DEFAULT){
						te += " part-per-core="+part_per_core;
					}
					
					if(Number(details.npart)%part_per_core != 0){
						add_warning({mess:"Error with particle number", mess2:"The number of particles '"+details.npart+"' must be a multiple of the number of particles per core '"+part_per_core+"'", warn_type:"InfFurtherDetails"});
					}
				}
				
				//if(debug){
					//te += ' diagnostics="on"';
				//}
			}
			break;
		
		case "ABC":
			te += " sample="+details.abcsample;
			te += " acc-frac="+details.accfrac;
			break;
			
		case "ABC-SMC":
			te += " sample="+details.abcsample;
			te += " acc-frac="+details.accfracsmc;
			te += " gen="+details.numgen;
			te += " kernel-size="+details.kernelsize;
			break;
		}
		
		if(details.indmax != INDMAX_DEFAULT){
			te += " ind-max="+details.indmax;
		}
		
		if(details.param_output_max != PARAM_OUTPUT_MAX_DEFAULT){
			te += " param-output-max="+details.param_output_max;
		}
		
		if(details.algorithm.value != ALG_DEFAULT){
			te += ' algorithm="'+details.algorithm.value+'"';
		}
		
		switch(details.algorithm.value){
		case "DA-MCMC":
			{	
				let nchain = Number(details.nchain);
				te += " nchain="+nchain;
				
				let type = details.anneal_type.te;
				
				if(type != ANNEAL_DEFAULT){
					te += ' anneal="' +type+'"';
					switch(type){
					case "scan":
						te += ' rate="' +details.anneal_rate+'"';
						break;
						
					case "power":
						te += ' power="' +details.anneal_power+'"';
						break;
					}
				}
				
				if(type != "scan"){
					if(Number(details.burnin_frac) != BURNIN_FRAC_DEFAULT){
						te += ' burnin-percent="' +details.burnin_frac+'"';
					}
				}
			}
			break;
				
		case "PAS-MCMC":
			te += " npart="+details.npart;
			te += ' gen-percent=' +details.gen_update;
			if(Number(details.burnin_frac) != BURNIN_FRAC_DEFAULT){
				te += ' burnin-percent="' +details.burnin_frac+'"';
			}
			break;
		}
		
		te += endl;
		te += endl;
	}
	
	{                                               // Deals with ppc details
		te += get_ppc_command(save_type)
		te += endl;
		te += endl;
	}

	return te;
}


/// Creates text for the ppc command
function get_ppc_command(save_type)
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
	
	let te = "post-sim";
	if(details.ppc_t_start != "") te +=" start="+details.ppc_t_start;
	if(details.ppc_t_end != "") te += " end="+details.ppc_t_end;
	
	if(details.number != ""){
		let num = Number(details.number);
		if(num != PPC_NUM_DEFAULT){
			te += " number="+details.number;
		}
	}
	
	if(details.seed_on.value == "Yes") te += " seed="+Number(details.seed);
	
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
		
		if(fl){
			te += ' resample="'+st+'"';	
		}
	}
	
	return te;
}


/// Adds information about simulation and inference
function create_output_sim_inf_source(p,sim_or_inf,file_list,one_file)
{
	let te = "";
	
	let sp = model.species[p];
	if(sim_or_inf == "ppc"){
		if(!model.inf_res.plot_filter) return te;
		sp = model.inf_res.plot_filter.species[p];
		if(!sp) return te;
	}
	
	let source, pa;
	
	switch(sim_or_inf){
	case "sim": pa = "Simulation"; source = sp.sim_source; te += mini_banner("SIMULATION INITIAL CONDITIONS"); break;
	case "inf": pa = "Inference"; source = sp.inf_source; te += mini_banner("INFERENCE DATA"); break;
	case "ppc": pa = "Population Mod."; source = sp.ppc_source; te += mini_banner("MODIFICATION DATA"); break;
	default: error("Option not recognised 116"+sim_or_inf); break;
	}
	
	if(source.length == 0) return "";

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
					
					te += com;
					
					switch(ic_type){
					case "Fixed": te += ' type="fixed"'; break;
					case "Dist": te += ' type="dist"'; break;
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
							
							let file = get_unique_file(com+"-"+sp.name,file_list,'.csv');
							output_table_simp_file(file,head,row,file_list);
							if(one_file) file = get_one_file(file_list);
							te += ' focal="'+so.spec.focal.te+'" file="'+file+'"'+endl;
							
						}
						break;
							
					case "All":
						if(so.spec.radio_dist.value == "Dist"){
							te += " prior="+get_prior_output(sim_or_inf,so.pop_dist,"prior",undefined);
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
					
						let file = get_unique_file(com+"-"+sp.name,file_list,'.csv');

						output_table_simp_file(file,head,row,file_list);
						if(one_file) file = get_one_file(file_list);
						te += ' file="'+file+'"'+endl;
						break;
					}
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
				te += add_data_table(so.type,so.table,so,p,i,file_list,sim_or_inf,one_file);
			}
		}
		te += endl;
	}
	te += endl;

	return te;
}


/// Gets that last thing added to file_list and outputs directly into the file
function get_one_file(file_list)
{
	if(file_list.length == 0) alert_help("Problem with file");
	let te = file_list[file_list.length-1].data;
	file_list.pop();
	return "[[\n"+te+"]]";
}


/// Generates a data-table and command from a data source
function add_data_table(type,tab,so,p,index,file_list,sim_or_inf,one_file)
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
	
	let te = com+' name=\"'+so.name+'\" ';
	
	let spec = so.spec;
	
	switch(so.type){                                // Adds information dependent on data source
	case "Init. Pop.":
		if(so.spec.radio2.value == "Focal"){
			te +=  'focal="'+so.spec.focal.te+'" ';
		}
		break;
	
	case "Move Ind.":
		{
			let cl_name = so.spec.cl_drop.te;
			te += 'class="'+cl_name+'" ';
		}
		break;
		
	case "Compartment":
		{
			let cl_name = so.spec.cl_drop.te;
			te += 'class="'+cl_name+'" ';
			file += '-'+cl_name.replace(/ /g,"-");
		}
		break;
		
	
	case "Transition": 
		{
			let sp = model.species[p];
		
			let cl = find(sp.cla,"name",spec.cl_drop.te);

			let trans_fi = output_add_trans_filt(spec.filter,p,cl,so,index);
			if(trans_fi != "") te += 'trans="'+trans_fi+'" ';
			
			let filt = output_add_comp_filt(spec,sp,cl,p,so,index);
			if(filt != "") te += 'filter="'+filt+'" ';

			let obsrange = spec.time_radio.value.toLowerCase();
			
			switch(spec.time_radio.value){
			case "All": obsrange = "all"; break;
			case "Spec": obsrange = "specify"; break;
			case "File": obsrange = "file"; break;
			default: error("Option not recognised 1C"); break;
			}
			
			let ty = so.type;
			te += 'obsrange="'+obsrange+'" ';
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
					
					te += 'start="'+start+'" ';
					te += 'end="'+end+'" ';
				}
				break;

			case "file": break;

			default: error("Option not recognised 121"+obsrange); break;
			}
		}
		break;
		
	case "Pop. Trans.":
		{
			let sp = model.species[p];
		
			let cl = find(sp.cla,"name",spec.cl_drop.te);

			let trans_fi = output_add_trans_filt(spec.filter,p,cl,so,index);
			if(trans_fi != "") te += 'trans="'+trans_fi+'" ';
			
			let filt = output_add_comp_filt(spec,sp,cl,p,so,index);
			if(filt != "") te += 'filter="'+filt+'" ';

			te += add_obs_model(spec);
		}
		break;

	case "Diag. Test":
		{
			let err_fl;
			
			let Se = esc(so.spec.Se_eqn.te); if(Se == "") err_fl = "Se";
			te += 'Se="'+Se+'" ';
			
			let Sp = esc(so.spec.Sp_eqn.te); if(Sp == "") err_fl = "Sp";
			te += 'Sp="'+Sp+'" ';
			
			let pos = esc(so.spec.pos_result); if(pos == "") err_fl = "positive value";
			te += 'pos="'+pos+'" ';
			
			let neg = esc(so.spec.neg_result); if(neg == "") err_fl = "negative value";
			te += 'neg="'+neg+'" ';
			
			if(err_fl != undefined){
				add_warning({mess:"Problem with diagnostic test data", mess2:"The value for '"+eff_fl+"' is not set", warn_type:"SourceProb", p:p, ind:index});
		
			
			}
			
			te += 'comp="';
			let flag = false;
			
			let cb = so.spec.check_box;
			
			let pcl = get_p_cl_from_claname(cb.name);
			let claa = model.species[pcl.p].cla[pcl.cl];
			
			for(let i = 0; i < cb.value.length; i++){
				if(cb.value[i].check == true){
					if(flag == true) te += ',';
					te += claa.comp[i].name;
					flag = true;
				}
			}
			te += '" ';
		}
		break;
		
	case "Ind. Eff.":
		{
			te += 'ie="'+so.spec.drop.te+'" ';
		}
		break;
		
	case "Ind. Group":
		{
			//te += 'name="'+so.spec.gname+'" ';
		}
		break;
		
	case "Population": 
		{
			let sp = model.species[p];
		
			let filt = output_add_comp_filt(spec,sp,undefined,p,so,index);
			if(filt != "") te += 'filter="'+filt+'" ';

			te += add_obs_model(spec);
		}
		break;
		
	case "Genetic":
		switch(so.spec.type_radio.value){
		case "matrix": 
			te += 'type="matrix" '; 
			break;
			
		case "snp":	
			te += 'type="snp" ';
			te += 'root="'+spec.snp_root+'" ';
			break;
		}
		te += 'mut-rate="'+spec.mut_rate_eqn.te+'" ';
		te += 'seq-var="'+spec.seq_var_eqn.te+'" ';
		break;

	default: break;	
	}

	file = get_unique_file(file,file_list,'.csv');
	
	output_table_file(file,tab,file_list);
	if(one_file) file = get_one_file(file_list);
	
	te += 'file="'+file+'"'+endl;
	
	return te;
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
	// Converts any greek letters to word version
	for(let i = 0; i < file.length; i++){
		let ch = file.substr(i,1);
		for(let j = 0; j < greek_latex.length; j++){
			if(greek_latex[j][1] == ch){
				file = file.substr(0,i)+greek_latex[j][0]+file.substr(i+1);
			}
		}
		
		if(ch == "→") file = file.substr(0,i)+"-"+file.substr(i+1);
	}
	
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


function add_obs_model(spec)
{
	switch(spec.obs_error_dist.value){
	case "Normal":		
		switch(spec.obs_error.value){
		case "percent": return 'error="normal:'+spec.percent+'%" ';	
		case "sd": return 'error="normal:'+spec.sd+'" ';
		case "file": return 'error="normal:file" ';
		default: error("Option not recognised 123"); break;
		}
		break;
	
	case "Poisson":		
		return 'error="poisson" ';
		
	case "Negative binomial":	 
		switch(spec.obs_error_p.value){
		case "p": return 'error="neg-binomial:'+spec.p+'" ';	
		case "file": return 'error="neg-binomial:file" ';
		default: error("Option not recognised 123"); break;
		}
		break;
		
	default: error("obs dist not found"); break;
	}
}


/// Outputs a table 
function output_table_file(file,tab,file_list)
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
			
			if(isNaN(el[i])) st += '"'+el[i]+'"';
			else st += el[i];
		}
		st += endl;
	}
	
	write_file_store(st,file,file_list);
}


/// Outputs a simple table 
function output_table_simp_file(file,head,row,file_list)
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
			if(isNaN(ele[i])) st += '"'+ele[i]+'"';
			else st += ele[i];
		}
		st += endl;
	}
	
	write_file_store(st,file,file_list);
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
		
		let ninitpop = 0, naddind = 0, naddpop = 0;
		
		for(let i = 0; i < source.length; i++){
			let so = source[i];
			
			switch(so.type){
			case "Init. Pop.": ninitpop++; break;
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
			if(ninitpop > 1){
				add_warning({mess:"Initial population", mess2:"Only one initial population should be set for species '"+sp.name+"'", warn_type:"InfPopulationProb", p:p});
			}

			switch(sp.type){
			case "Population":
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


/// Ensures that line breaks are escaped 
function esc(st)
{
	st = st.replace(/\n/g,"");
	st = st.trim();
	return st;
}


/// Removes specified commands from a file
function remove_command(te,com)
{
	let pos = com.split(",");
	
	let lines = te.split("\n");
	
	let te_new = "";
	
	for(let li = 0; li < lines.length; li++){
		let line = lines[li];
		
		let fl = false;
			
		let i = 0; while(i < line.length && line.substr(i,1) == " ") i++;
	
		for(let j = 0; j < pos.length; j++){
			if(i+pos[j].length <= line.length && line.substr(i,pos[j].length) == pos[j]){
				if(i+pos[j].length < line.length){
					if(line.substr(i+pos[j].length,1) == " ") fl = true;
				}
				else fl = true;
			}
		}
	
		if(fl == false) te_new += line+endl;
		else{
			let tri = line.trim();
			if(tri.length > 3 && tri.substr(tri.length-3,3) == '"[['){
				while(li < lines.length){
					let tri = lines[li].trim();
				
					if(tri.length >= 3 && tri.substr(0,3) == ']]"'){
						if(tri.substr(tri.length-3,3) != '"[[') break;
					}
					li++;
				}
			}
		}
	}
	
	return te_new;
}


/// Tidys the code
function tidy_code(te)
{
	let spl = te.split("\n");
	
	// Trims lines
	for(let li = 0; li < spl.length; li++){
		spl[li]= spl[li].trim();
	}
	
	// Gets rid of double spaces
	let li = 0;
	while(li < spl.length-1){
		if(spl[li] == "" && spl[li+1] == ""){
			spl.splice(li,1);
		}
		else li++;
	}
	
	let te_new="";
	for(let li = 0; li < spl.length; li++){
		te_new += spl[li]+endl;
	}
	
	return te_new;
}

 
/// Creates a file to start a posterior simulation
function create_ppc_file()
{
	let te = inf_result.import_te;
	
	let file_list=[];
	
	te = remove_command(te,"post-sim,posterior-simulation,add-pop-post-sim,remove-pop-post-sim,add-ind-post-sim,remove-ind-post-sim,move-ind-post-sim,param-mult,# MODIFICATION DATA");
	te += endl+get_ppc_command("ppc")+endl+endl;

	for(let p = 0; p < inf_result.species.length; p++){
		let sp = inf_result.species[p];
		if(inf_result.species.length > 1) te += 'set species="'+sp.name+'"'+endl;
		te += create_output_sim_inf_source(p,"ppc",file_list,true);
	}
	
	let pfac_list=[];
	for(let th = 0; th < model.param.length; th++){
		let par = model.param[th];
		if(par.type == "param factor") pfac_list.push(th);
	}
	
	if(pfac_list.length > 0){
		te += mini_banner("PARAMETER MULTIPLIERS");
		for(let k = 0; k < pfac_list.length; k++){
			let par = model.param[pfac_list[k]];
			if(par.type == "param factor"){
				te += output_param(par,"ppc",file_list,true);
			}
		}
		te += endl;
	}
 
	te = tidy_code(te);
	
	if(model.warn.length > 0) err_warning();
	
	write_file_store(te,"bicifile",file_list,"bicifile");

	post({save_type:"ppc", file_list:file_list, param_factor:strip_heavy(model.param_factor), param:strip_heavy(model.param), species:strip_heavy(model.species)});
}

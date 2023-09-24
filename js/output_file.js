"use strict";

/// Creates an output file
function create_output_file(type,filename)
{
	let file_list=[];  // Used to make sure that files to not share the same name
	
	let dir = "data-files";
	if(filename != undefined){
		let spl = filename.split(".");
		dir = spl[0]+"-data-files";
	}
	if(ver == "mac") dir = "/tmp/"+dir;
		
	const fs = require('fs');  // Creates the data directory
	if (!fs.existsSync(dir)) {
		fs.mkdir(dir, function(err) {
			if(err) { 
				alertp("There was a problem creating the directoy '"+dir+"'");
			}
		});
	}
	
	init_param();
	if(model.warn.length > 0){ generate_screen(); return;}
	
	let map_data = { crs:{},features:[],type:"FeatureCollection"};
	
	let exporting = false;
	if(type == "siminf") exporting = true;
	
	let te = create_output_details(type,dir);

	for(let p = 0; p < model.species.length; p++){
		te += create_output_species(p,map_data,file_list,type,dir,exporting);
	}

	te += create_output_param(type,dir,file_list,exporting);
		
	te += create_output_derived();

	check_data_valid(type);
	
	if(model.warn.length > 0){ generate_screen(); return;}
	
	check_comp_structure();

	turn_off_import_warning();
	
	if(model.warn.length > 0){ generate_screen(); return;}
		
	if(filename == "NO_SAVE") return te;
	
	if(map_data.features.length > 0){
		write_file(JSON.stringify(map_data),dir+"/map-data.json");
	}

	if(ver == "mac") filename = "/tmp/"+filename;
	
	write_file(te,filename);
	
	return "success";
}


/// Turns off any additional import warnings which get set during exporting
function turn_off_import_warning()
{
	if(inter.help.title == "Error importing file!") close_help();
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
function create_output_details(type,dir)
{
	let te = banner("DETAILS");
		
	let desc = model.description.te.replace(/\n/g,"|");
		
	te += 'description text="'+desc+'"'+endl;
	te += endl;
	
	te += 'data-dir folder="'+dir+'"'+endl;
	te += endl;
	
	te += create_output_siminf(type);
	return te;
}
	
	
/// Outputs all the results for a given species
function create_output_species(p,map_data,file_list,type,dir,exporting)
{
	let sp = model.species[p];
		
	let te = banner("DEFINE MODEL AND DATA FOR SPECIES "+sp.name.toUpperCase());
	
	te += 'species name="'+sp.name+'" type="'+sp.type.toLowerCase()+'"'+endl+endl;
	
	te += create_output_compartments(p,map_data,file_list,dir);
	
	te += create_output_ind_eff(p,dir,file_list);

	te += create_output_fix_eff(p,dir,file_list);

	if(type == "sim" || type == "siminf"){
		te += create_output_sim_inf_source(p,"sim",dir,file_list,exporting);
		output_check("sim",exporting);	
	}

	if(type == "inf" || type == "siminf"){
		te += create_output_sim_inf_source(p,"inf",dir,file_list,exporting);
		output_check("inf",exporting);	
	}
	
	return te;
}


/// Creates output for all the compartments in the model
function create_output_compartments(p,map_data,file_list,dir)
{
	let te = mini_banner("COMPARTMENTS");
	
	let sp = model.species[p];
	for(let cl = 0; cl < sp.ncla; cl++){
		let claa = sp.cla[cl];
		let cam = claa.camera;
		
		te += 'class name="'+claa.name+'"';
		
		if(claa.default_map == true) te += ' map="load"'; 
		
		let clp = clone_p(p,cl);

		if(clp != undefined) te += ' clone="'+model.species[clp].name+'"';
		else{
			te += ' index="'+claa.index+'" coord="'+cam.coord+'"';
			
			te += endl+endl;
		
			te += 'camera';
			switch(cam.coord){
			case "cartesian": 
				te += ' x='+precision(cam.x)+' y='+precision(cam.y); 
				break;
				
			case "latlng": 
				let p = transform_latlng_inv(cam.x,cam.y);
				te += ' lat='+precision(p.lat)+' lng='+precision(p.lng); 
				break;

			default: error("Option not recognised 110"); break;
			}
			
			te += ' scale='+precision(cam.scale)+endl+endl;
			
			for(let c = 0; c < claa.ncomp; c++){
				let co = claa.comp[c];
				
				te += 'comp name="'+co.name+'"';
				switch(cam.coord){
				case "cartesian": 
					te += ' x='+precision(co.x/import_scale_factor);
					te += ' y='+precision(co.y/import_scale_factor);
					break;
					
				case "latlng":
					if(co.type == "boundary"){
						te += ' boundary="map-data.json"';
						
						let fea = generate_JSON_feature(co.name,co.feature);

						map_data.features.push(fea);
					}
					else{
						let p = transform_latlng_inv(co.x,co.y);
						te += ' lat='+precision(p.lat)+' lng='+precision(p.lng); 
					}
					break;

				default: error("Option not recognised 111"); break;
				}
		
				te += ' color="'+co.col+'"';
		
				te += ' fix="'+co.fixed.check+'"';

				if(co.choose_branch == true){
					if(co.branch == true) te += ' branch-prob="true"';
					else te += ' branch-prob="false"';
				}
				te += endl;
			}
			te += endl;
		}
	
		for(let j = 0; j < claa.ntra; j++){
			let tr = claa.tra[j];
			
			let mid_start = 0, mid_end = tr.midp.length;
			let i_sel;
			
			switch(tr.variety){
			case "Normal":
				te += 'trans from="'+claa.comp[tr.i].name+'"';
				te += ' to="'+claa.comp[tr.f].name+'"';
				break;
				
			case "Source":
				te += 'source to="'+claa.comp[tr.f].name+'"';
				mid_start = 1; i_sel = 0;
				break;
			
			case "Sink":
				te += 'sink from="'+claa.comp[tr.i].name+'"';
				mid_end--; i_sel = mid_end;
				break;

			default: error("Option not recognised 112"); break;
			}
			
			if(i_sel != undefined){
				switch(cam.coord){
				case "cartesian":
					te += ' x="'+precision(tr.midp[i_sel].x/import_scale_factor)+'"';
					te += ' y="'+precision(tr.midp[i_sel].y/import_scale_factor)+'"';
					break;
				
				case "latlng":
					let p = transform_latlng_inv(tr.midp[i_sel].x,tr.midp[i_sel].y);
					te += ' lat="'+precision(p.lat)+'" lng="'+precision(p.lng)+'"';
					break;

				default: error("Option not recognised 113"); break;
				}
			}
			
			te += ' type="'+tr.type+'"';
			switch(tr.type.toLowerCase()){
			case "exp(rate)": 
				te += ' rate="'+esc(tr.value.rate_eqn.te)+'"';
				break;
			
			case "exp(mean)": 
				te += ' mean="'+esc(tr.value.mean_eqn.te)+'"';
				break;
				
			case "gamma":
				te += ' mean="'+esc(tr.value.mean_eqn.te)+'"';
				te += ' cv="'+esc(tr.value.cv_eqn.te)+'"';
				break;
				
			case "erlang":
				te += ' mean="'+esc(tr.value.mean_eqn.te)+'"';
				te += ' shape="'+esc(tr.value.shape_erlang.te)+'"';
				break;
				
			case "log-normal": 
				te += ' mean="'+esc(tr.value.mean_eqn.te)+'"';
				te += ' cv="'+esc(tr.value.cv_eqn.te)+'"';
				break;
				
			case "weibull":
				te += ' scale="'+esc(tr.value.scale_eqn.te)+'"';
				te += ' shape="'+esc(tr.value.shape_eqn.te)+'"';
				break;

			default: error("Option not recognised 114"); break;
			}

			if(tr.branch == true){
				te += ' bp="'+esc(tr.value.bp_eqn.te)+'"';
			}
			
			switch(cam.coord){
			case "cartesian":
				let xlist = '', ylist = '';
				
				for(let i = mid_start; i < mid_end; i++){
					if(i != mid_start){ xlist += ','; ylist += ',';}
					
					xlist += precision(tr.midp[i].x/import_scale_factor);
					ylist += precision(tr.midp[i].y/import_scale_factor);
				}
				
				if(xlist != '') te += ' mid-x="'+xlist+'" mid-y="'+ylist+'"';
				break;
				
			case "latlng":
				let latlist = '', lnglist = '';
				
				for(let i = mid_start; i < mid_end; i++){
					if(i != mid_start){ latlist += ','; lnglist += ',';}
					
					let p = transform_latlng_inv(tr.midp[i].x,tr.midp[i].y);

					latlist += precision(p.lat);
					lnglist += precision(p.lng);
				}
				
				if(latlist != '') te += ' mid-lat="'+latlist+'" mid-lng="'+lnglist+'"';
				break;

			default: error("Option not recognised 115"); break;
			}

			te += endl;
		}
		te += endl;
	}
	
	return te;
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
						
						//if(comp.choose_branch != comp2.choose_branch) prob = 5;
						
						//if(comp.branch != comp2.branch) prob = 6;
					}
					if(prob == 0) return p2;
				}
			}
		}			
	}
}


/// Sets the precision of the output	
function precision(x,digit)
{
	if(isNaN(x)) error(x+" is not a number");
	x = Number(x);

	if(digit == undefined) digit = 6;
	
	let mod = x; if(mod < 0) mod = -mod; 
	
	let digit_min = Math.floor(1+Math.log10(mod));
	if(digit_min < 1) digit_min = 1;

	if(digit < digit_min) digit = digit_min;

	let val = Number(x).toPrecision(digit);
	for(let j = digit-1; j >= digit_min; j--){
		let val2 = x.toPrecision(j);
		if(Number(val) == Number(val2)) val = val2; else break;
	}
 
	return val;
}


/// Create output of parameters
function create_output_param(type,dir,file_list,exporting)
{
	let te = banner("PARAMETERS");
	
	for(let th = 0; th < model.param.length; th++){
		let par = model.param[th];
		
		let dist_done = false;
		
		if(!(type == "sim" && find_in(sim_param_not_needed,par.type) != undefined) &&
			 !(type == "inf" && find_in(inf_param_not_needed,par.type) != undefined)){
			let num_warn = model.warn.length;
		
			let te1 = 'param name="'+remove_eq_quote(par.full_name); te1 += '"';
			
			let te2 = "";
			let display = false;
				
			if(par.dist_matrix.check == true){
				te1 += ' distance-matrix="true"';
				display = true;
			}
			else{
				// Constants and simulation parameters
				if(par.variety == "const" || par.variety == "reparam" || type == "sim" || type == "siminf"){ 
					if((par.dep.length == 0 && par.value == set_str) || (par.set == false && !(par.variety == "dist" && par.sim_sample.check == true))){
						if(exporting == false){
							model.warn.push({mess:"A value for "+par.full_name+" must be set.", mess2:"Parameter values must be set before simulation can be perfomed", warn_type:"MissingSimValue", name:par.name});
						}
					}
					else{
						if(par.variety == "dist" && par.sim_sample.check == true){
							te2 += output_add_prior_distribution(par,"dist",file_list,dir);
							dist_done = true;
						}
						else{
							let exp = "value";
							if(par.variety == "const") exp = 'constant';
							else{
								if(par.variety == "reparam") exp = 'reparam';
							}
							te2 += ' '+exp+'=';
							
							if(exp == "value" && par.auto_value == true){
								te2 += '"auto"';
							}
							else{
								let file;				
								if(par.dep.length > 0){
									file = get_unique_file(exp+"-"+par.name,file_list,'.csv');
								}
								
								if(file != undefined){
									te2 += '"'+file+'"';
									
									let data = '';
									for(let d = 0; d < par.dep.length; d++){	
										data += '"'+par.dep[d]+'",';
									}
									data += 'Value' + endl;
								
									for(let i = 0; i < par.comb_list.length; i++){
										let comb = par.comb_list[i];
										for(let d = 0; d < par.dep.length; d++){
											data += '"'+comb.combination[d]+'",';
										}
										data += '"'+get_element(par.value,comb.index)+'"'+endl;
									}
								
									write_file(data,dir+"/"+file);
								}
								else{
									te2 += JSON.stringify(par.value); 
								}
							}
						}	
					}
				}
				
				if(par.variety != "const" && par.variety != "reparam" &&
					(type == "inf" || type == "siminf")){
					
					switch(par.variety){
					case "normal":
						te2 += output_add_prior_distribution(par,"prior",file_list,dir);
						break;
						
					case "dist":
						if(dist_done == false){
							te2 += output_add_prior_distribution(par,"dist",file_list,dir);
						}
						break;
						
					default: error("Option not recognised 1B"+par.variety); break;
					}
				}
				
				if(te2 != "") display = true;
				
				if(par.spline.on == true){
					if(par.spline.knot.length == 0) error("There should be knot times specified");
					
					te2 += ' knot-times="'+stringify(par.spline.knot)+'"';
					
					let sm = par.spline.smooth;
					if(sm.check == true){
						te2 += ' smooth="'+sm.type.value.toLowerCase()+'('+sm.value+')"';
					}
				}
			}
				
			if(display == true){
				te += te1+te2+endl;
			}
			else{
				if(num_warn == model.warn.length){
					model.warn.push({mess:"Parameter "+par.full_name+" could not be exported.", mess2:"Unknown error", warn_type:"", name:par.name});
				}
			}
		}
	}
	te += endl;

	return te;
}


/// Adds text associated with adding a distribution
function output_add_prior_distribution(par,variety,file_list,dir,)
{
	let st="";

	let key;
	switch(variety){
	case "dist": key = "distribution"; break;
	case "prior": key = "prior"; break;
	}
		
	if(par.dep.length > 0 && par.prior_split_check.check == true){  /// Prior split up
		let warn_type;
		switch(variety){
		case "dist": key = "distribution"; warn_type = "MissingDistValue"; break;
		case "prior": key = "prior"; warn_type = "MissingPriorValue";  break;
		default: error("option prob"); break;
		}
			
		if(par.prior_split_set == false){
			model.warn.push({mess:"Missing "+key+" distribution information", mess2:"The "+key+" for "+par.full_name+" must be set.", warn_type:warn_type, name:par.name});
		}
		else{
			let file = get_unique_file(key+"-"+par.name,file_list,'.csv');
			
			let col;
			
			switch(variety){
			case "dist": st += ' dist-split'; col = 'Dist'; break;
			case "prior": st += ' prior-split'; col = 'Prior'; break;
			default: error("option prob"); break;
			}
			
			st += '="'+file+'"';
								
			let data = '';
			for(let d = 0; d < par.dep.length; d++){
				data += '"'+par.dep[d]+'",';
			}
			data += col+endl;
			
			//let list = generate_list_from_dep(par);
			
			//for(let i = 0; i < list.length; i++){
			for(let i = 0; i < par.comb_list.length; i++){
				let comb = par.comb_list[i];
				for(let d = 0; d < par.dep.length; d++){
					data += '"'+comb.combination[d]+'",';
				}
				
				let pri = get_element(par.prior_split,comb.index);
				
				if(pri.type.te == select_str){
					model.warn.push({mess:"Missing "+key, mess2:"The "+key+" for "+par.full_name+" must be set.", warn_type:"MissingPriorSplitValue", name:par.name});
					return "";
				}
				
				data += '"'+get_prior_string(pri)+'"'+endl;
			}
		
			write_file(data,dir+"/"+file);
		}
	}
	else{
		let pri = par.prior;
	
		let tex = pri.type.te;
		if(tex == select_str){
			model.warn.push({mess:"Missing "+key, mess2:"The "+key+" for "+par.full_name+" must be set.", warn_type:"MissingPriorValue", name:par.name});
			return "";
		}
		
		let te = '"'+tex.toLowerCase()+'"';
	
		switch(pri.type.te){
		case "fix":
			{
				let mean = esc(pri.value.mean_eqn.te);
				check_prior_val(mean,variety,par,{positive:true});
				te += ' fixed="'+mean+'"';
			}
			break;
			
		case "uniform":
			{
				let min = esc(pri.value.min_eqn.te);
				check_prior_val(min,variety,par,{});
				te += ' min="'+min+'"';
				
				let max = esc(pri.value.max_eqn.te);
				check_prior_val(max,variety,par,{});
				te += ' max="'+max+'"';
			}
			break;

		case "exp":
			{
				let mean = esc(pri.value.mean_eqn.te);
				check_prior_val(mean,variety,par,{positive:true});
				te += ' mean="'+mean+'"';
			}
			break;

		case "normal":
			{
				let mean = esc(pri.value.mean_eqn.te);
				check_prior_val(mean,variety,par,{});
				te += ' mean="'+mean+'"';
				
				let sd = esc(pri.value.sd_eqn.te);
				check_prior_val(sd,variety,par,{positive:true});
				te += ' sd="'+sd+'"';
			}
			break;

		case "log-normal":
			{
				let mean = esc(pri.value.mean_eqn.te);
				check_prior_val(mean,variety,par,{positive:true});
				te += ' mean="'+mean+'"';
				
				let cv = esc(pri.value.cv_eqn.te);
				check_prior_val(cv,variety,par,{positive:true});
				te += ' cv="'+cv+'"';
			}
			break;
			
		case "gamma":
			{
				let mean = esc(pri.value.mean_eqn.te);
				check_prior_val(mean,variety,par,{positive:true});
				te += ' mean="'+mean+'"';
				
				let cv = esc(pri.value.cv_eqn.te);
				check_prior_val(cv,variety,par,{positive:true});
				te += ' cv="'+cv+'"';
			}
			break;
		
		case "beta":
			{
				let alpha = esc(pri.value.alpha_eqn.te);
				check_prior_val(alpha,variety,par,{positive:true});
				te += ' alpha="'+alpha+'"';
				
				let beta = esc(pri.value.beta_eqn.te);
				check_prior_val(beta,variety,par,{positive:true});
				te += ' beta="'+beta+'"';
			}
			break;
			
		case "bernoulli":
			{
				let mean = esc(pri.value.mean_eqn.te);
				check_prior_val(mean,variety,par,{positive:true});
				te += ' mean="'+mean+'"';
			}
			break;
			
		case "flat":
			break;
			
		case "dirichlet":
			{
				let alpha = esc(pri.value.alpha_eqn.te);
				check_prior_val(alpha,variety,par,{positive:true});
				te += ' alpha="'+alpha+ '"';
			}
			break;
		
		default: error("op error"); break;
		}
	
		switch(variety){ 
		case "dist": st += ' dist'; break;
		case "prior": st += ' prior'; break;
		default: error("option prob"); break;
		}
		
		st += '='+te;
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
				if(Number(val) <= 0) warn_text = "Must be positive";
			}
		}
		break;

	case "dist":
		key = "distribution"; warn_type = "MissingDistValue"; 
		if(!is_eqn(val,"min",op)) warn_text = "'"+val+"' is not valid.";
		break;
		
	default: error("not op"); break;
	}

	if(warn_text != undefined){
		model.warn.push({mess:"Error in "+key+" for "+par.full_name, mess2:warn_text, warn_type:warn_type, name:par.name});
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
			let name = remove_eq_quote(get_full_parameter_name(der.eqn1.param[0].name));
			te += 'derived name="'+name+'" eqn="'+der.eqn2.te+'"'+endl;
		}
	}
	te += endl;

	return te;
}


/// Create output of individual effects
function create_output_ind_eff(p,dir,file_list)
{
	let iegs = model.species[p].ind_eff_group;
		
	if(iegs.length == 0) return "";
	
	let te = mini_banner("INDIVIDUAL EFFECTS");

	for(let i = 0; i < iegs.length; i++){
		let ieg = iegs[i];
		let name = "";
		for(let j = 0; j < ieg.list.length; j++){
			if(j > 0) name += ",";
			name += ieg.list[j].name;
		}
			
		te += 'ind-effect name="'+name+'"';
		
		let mat = ieg.A_matrix;
		if(mat.check == true){
			if(mat.loaded == false){
				model.warn.push({mess:"Error in individual effects", mess2:"The A matrix for individual effects '"+name+"' must be loaded", warn_type:"A matrix"});
			}
			else{
				let file = get_unique_file("A-matrix-"+name,file_list,'.csv');
	
				te += ' A="'+file+'"';
				
				let data = '';
				for(let c = 0; c < mat.ind_list.length; c++){ 
					if(c != 0) data += ',';
					data += '"'+mat.ind_list[c]+'"';
				}
				data += endl;

				for(let r = 0; r <  mat.ind_list.length; r++){
					for(let c = 0; c < mat.ind_list.length; c++){ 
						if(c != 0) data += ',';
						data += mat.value[r][c];
					}
					data += endl;
				}
	
				write_file(data,dir+"/"+file);
			}
		}		
		te += endl;
	}

	te += endl;

	return te;
}


/// Create output of fixed effects
function create_output_fix_eff(p,dir,file_list)
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
			model.warn.push({mess:"Error in fixed effects", mess2:"The vector <e>X^"+name+"</e> for fixed effect '"+name+"' must be loaded", warn_type:"X vector"});
		}
		else{
			let file = get_unique_file("X-vector-"+name,file_list,'.csv');

			te += ' X="'+file+'"';
			
			let data = 'ID,value'+endl;
			for(let c = 0; c < vec.ind_list.length; c++){ 
				data += '"'+vec.ind_list[c]+'",'+vec.value[c]+endl;
			}
			write_file(data,dir+"/"+file);
		}
		te += endl;
	}

	te += endl;

	return te;
}


/// Creates information about simulation or inference
function create_output_siminf(type)
{
	let te = "";
	
	if(type == "sim" || type == "siminf"){
		let flag = false;
		
		let details = model.sim_details;
		
		if(details.t_start == "" || isNaN(details.t_start)){
			if(type == "siminf") flag = true;
			else{
				model.warn.push({mess:"Error in simulation details", mess2:"The start time must be set", warn_type:"SimDetails"});
			}
		}
		
		if(details.t_end == "" || isNaN(details.t_end)){
			if(type == "siminf") flag = true;
			else{
				model.warn.push({mess:"Error in simulation details", mess2:"The end time must be set", warn_type:"SimDetails"});
			}
		}
		
		let num = Number(details.number);
		if(details.number == "" || isNaN(num) || num != Math.floor(num) || num <= 0){
			if(type == "siminf") flag = true;
			else{
				model.warn.push({mess:"Error in simulation details", mess2:"The simulation number must be set to a positive integer", warn_type:"SimDetails"});
			}
		}
		
		let timestep = Number(details.timestep);
		if(details.timestep == "" || isNaN(timestep) || timestep <= 0){
			if(type == "siminf") flag = true;
			else{
				model.warn.push({mess:"Error in simulation details", mess2:"The simulation timestep must be set to a positive number", warn_type:"SimDetails"});
			}
		}
		
		if(flag == false){
			te += "simulation start="+details.t_start+" end="+details.t_end+" number="+num;
			//if(timestep != "") 
			te += " timestep="+timestep; 
			
			te += " ind-max="+details.indmax;
			
			te += endl;
			te += endl;
		}
	}	

	if(type == "inf" || type == "siminf"){ 
		let flag = false;
		
		let details = model.inf_details;
		
		if(details.t_start == "" || isNaN(details.t_start)){
			if(type == "siminf") flag = true;
			else{
				model.warn.push({mess:"Error in inference details", mess2:"The start time must be set", warn_type:"InfDetails"});
			}
		}
		
		if(details.t_end == "" || isNaN(details.t_end)){
			if(type == "siminf") flag = true;
			else{
				model.warn.push({mess:"Error in inference details", mess2:"The end time must be set", warn_type:"InfDetails"});
			}
		}
	
		let timestep = Number(details.timestep);
		if(details.timestep == "" || isNaN(timestep) || timestep <= 0){
			if(type == "siminf") flag = true;
			else{
				model.warn.push({mess:"Error in inference details", mess2:"The timestep must be set to a positive number", warn_type:"InfDetails"});
			}
		}
		
		if(flag == false){
			te += "inference start="+details.t_start+" end="+details.t_end;
			if(timestep != "") te += " timestep="+timestep; 
			
			
			switch(details.algorithm.value){
			case "DA-MCMC":
				te += " sample="+details.sample;
				te += " param-thin="+details.thinparam;
				te += " state-thin="+details.thinstate;
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
			
			te += " ind-max="+details.indmax;
				
			te += ' algorithm="'+details.algorithm.value+'"'+endl;
			te += endl;
		}
	}	

	return te;
}

/*
/// given a list of dependencies this provides a list going thorough each
function generate_list_from_dep(par)
{
	let dep = par.dep;
	
	let ndep = dep.length;
	let list = [];
	
	let pos = [];
	let index=[];
	for(let d = 0; d < ndep; d++){
		if(dep[d] == "t" || dep[d] == "a"){
			pos[d] = par.spline.knot;
		}
		else{
			pos[d] = find_comp_from_index(dep[d]);
		}
		index[d] = 0;
	}

	let flag;
	do{
		let comb = [];
		for(let d = 0; d < ndep; d++) comb.push(pos[d][index[d]]);
		
		list.push({combination:comb, index:copy(index)});
		
		let i = 0;
		do{
			flag = false;
			
			index[i]++; 
			if(index[i] >= pos[i].length){
				index[i] = 0; 
				i++; 
				flag = true;
			}
		}while(i < ndep && flag == true);
	}while(flag == false);
	
	return list;
}
*/


function create_output_sim_inf_source(p,type,dir,file_list,exporting)
{
	let sp = model.species[p];
	
	let source, pa;
		
	let te = "";
	
	switch(type){
	case "sim": pa = "Simulation"; source = sp.sim_source; te += mini_banner("SIMULATION DATA"); break;
	case "inf": pa = "Inference"; source = sp.inf_source; te += mini_banner("DATA"); break;
	default: error("Option not recognised 116"+type); break;
	}
	
	if(source.length == 0) return "";
	
	for(let i = 0; i < source.length; i++){
		let so = source[i];
	
		if(so.error == true){
			if(exporting == false){
				let te;
				if(so.table.filename == "") te = "Please delete data and reload.";
				else te = "Data table filename: "+so.table.filename;
				model.warn.push({mess:"Data source is invalid", mess2:te, warn_type:"SourceInvalid", siminf:type, p:p, ind:i});
			}
		}
		else{
			let add_flag = false;
			
			switch(so.type){
			case "Init. Pop.":
				let com = "init-pop"; if(type == "sim") com = "init-pop-sim";
	
				switch(so.spec.radio.value){
				case "Graphical":				
					let head = [];
					let row = [];
					
					switch(so.spec.radio2.value){
					case "Focal":
						if(so.cla.length != sp.ncla){
							model.warn.push({mess:"Initial population problem", mess2:"The number of classifications in the initial population for '"+sp.name+"' is not right",warn_type:"Init_pop", pa:pa, p:p, i:i});
						}
						else{
							for(let cl = 0; cl < sp.ncla; cl++){
								let claa = sp.cla[cl];

								let ip = so.cla[cl].comp_init_pop;
							
								if(ip.length != claa.ncomp){
									model.warn.push({mess:"Initial population problem", mess2:"The number of compartments in classification '"+claa.name+"' is not right",warn_type:"Init_pop", pa:pa, p:p, i:i});
								}
								else{
									for(let c = 0; c < ip.length; c++){
										if(ip[c].comp_name_store != claa.comp[c].name){
											if(ip[c].name == undefined){
												model.warn.push({mess:"Initial population problem", mess2:"The data source has become invalid",warn_type:"Init_pop", pa:pa, p:p, i:i});
											}
											else{
												model.warn.push({mess:"Initial population problem", mess2:"The names for the compartments '"+ip[c].name+"' and '"+claa.comp[c].name+"' do not agree",warn_type:"Init_pop", pa:pa, p:p, i:i});
											}
										}
									}
								}
							}
						}
						
						if(model.warn.length == 0){
							head.push("Compartment");
							head.push("Population");
						
							for(let cl = 0; cl < sp.ncla; cl++){
								let claa = sp.cla[cl];
								let cmax = claa.ncomp;
								if(claa.name != so.spec.focal.te) cmax--;
								
								for(let c = 0; c < cmax; c++){
									let ele = [];
									ele.push(claa.comp[c].name);
									if(claa.name == so.spec.focal.te){
										ele.push(so.cla[cl].comp_init_pop[c].pop);
									}
									else{
										ele.push(so.cla[cl].comp_init_pop[c].pop_per+"%");
									}
									row.push(ele);
								}
							}
							
							let file = get_unique_file(com+"-"+sp.name,file_list,'.csv');

							te += com+' focal="'+so.spec.focal.te+'" file="'+file+'"'+endl;
							output_table_simp_file(dir+"/"+file,head,row);
						}
						break;
							
					case "All":
						for(let cl = 0; cl < sp.ncla; cl++) head.push(sp.cla[cl].name);
						head.push("Population");
						
						for(let j = 0; j < so.glob_comp.length; j++){
							let gc = so.glob_comp[j];
							if(gc.cla.length != sp.ncla){
								model.warn.push({mess:"Initial population problem", mess2:"The number of classifications in the initial population for '"+sp.name+"' is not right",warn_type:"Init_pop", pa:pa, p:p, i:i});
								return;
							}
					
							let ele=[];
							for(let cl = 0; cl < sp.ncla; cl++) ele.push(gc.cla[cl]);
							ele.push(gc.pop);
							
							row.push(ele);
						}
					
						let file = get_unique_file(com+"-"+sp.name,file_list,'.csv');

						//let file = dir+"/"+get_unique_file("init-pop-"+sp.name,file_list,'.csv');

						te += com + ' file="'+file+'"'+endl;
						output_table_simp_file(dir+"/"+file,head,row);
						break;
					}
					break;
					
				case "File":
					add_flag = true;
					break;

				default: error("Option not recognised 120"); break;
				}
				break;
			
			case "Init. Pop. Prior":
				switch(so.variety){
				case "Flat": te += 'init-pop-prior type="flat"'+endl; break;
				case "Dirichlet": add_flag = true; break;
				}
				break;
				
			default:
				add_flag = true;
				break;
			}
			
			if(add_flag == true){
				te += add_data_table(so.type,so.table,so,p,i,file_list,type,dir);
			}
		}
	}
	te += endl;

	return te;
}


/// Generates a data-table and command from a data source
function add_data_table(type,tab,so,p,index,file_list,siminf,dir)
{
	let i = find(convert,"type",type);
	if(i == undefined){
		error("Counld not not create data table:"+type);
		return;
	}
	
	let com = convert[i].command;
	
	if(siminf == "sim"){
		switch(com){
		case "add-ind": com = "add-ind-sim"; break;
		case "remove-ind": com = "remove-ind-sim"; break;
		case "move-ind": com = "move-ind-sim"; break;
		case "init-pop": com = "init-pop-sim"; break;
		}
	}
	
	let file = com;
	
	let te = com+' ';
	
	let spec = so.spec;
	
	switch(so.type){  // Adds in extra information dependent on data source
	case "Init. Pop.":
		if(so.spec.radio2.value == "Focal"){
			te +=  'focal="'+so.spec.focal.te+'" ';
		}
		break;
	
	case "Init. Pop. Prior": 
		te += 'type="'+so.variety.toLowerCase()+'" ';
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
			file += '-'+cl_name.replace(" ","-");
		}
		break;
		
	case "Transition": case "Source": case "Sink":
		let name = spec.tr_drop.te;

		switch(so.type){
		case "Transition":
			let spl = name.split("→");
			if(spl.length != 2) error("Split not correct");
			te += 'from="'+spl[0]+'" to="'+spl[1]+'" ';
			break;
		case "Source": te += 'to="'+name+'" '; break;
		case "Sink": te += 'from="'+name+'" '; break;
		}
		
		file += '-'+name;
	
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
			let start = spec.time_start;
			if(isNaN(start)){	
				model.warn.push({mess:ty+" data", mess2:"The start time for the "+ty.toLowerCase()+" data must be set to a number", warn_type:"TransDataSpecProb", ind:index});
			}
			
			let end = spec.time_end;
			if(isNaN(end)){	
				model.warn.push({mess:ty+" data", mess2:"The end time for the "+ty.toLowerCase()+" data must be set to a number", warn_type:"TransDataSpecProb", ind:index});
			}
			
			te += 'start="'+start+'" ';
			te += 'end="'+end+'" ';
			break;

		case "file": break;

		default: error("Option not recognised 121"+obsrange); break;
		}
		break;
		
	case "Diag. Test":
		te += 'Se="'+esc(so.spec.Se_eqn.te)+'" ';
		te += 'Sp="'+esc(so.spec.Sp_eqn.te)+'" ';
		te += 'pos="'+esc(so.spec.pos_result)+'" ';
		te += 'neg="'+esc(so.spec.neg_result)+'" ';
		
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
		break;
		
	case "Population": case "Set Traps":
		{
			if(so.type == "Set Traps"){
				te += 'prob="'+esc(spec.trap_prob_eqn.te)+'" ';
			}
			
			let sp = model.species[p];
			
			if(spec.filter.cla.length != sp.ncla){
				model.warn.push({mess:so.type, mess2:"The population filter does not agree with the model classifications.", warn_type:"SourceProb", p:p, ind:index});
			}
			else{
				let filt='';
				for(let cl = 0; cl < sp.ncla; cl++){
					let clz = spec.filter.cla[cl];
					let claa = sp.cla[cl];
					switch(clz.radio.value){
					case "All": break;
					case "Comp": 
						filt += claa.name+':';
						let flag = false;
						for(let c = 0; c < claa.ncomp; c++){
							if(clz.comp[c].check == true){
								if(flag == true) filt += '|'; flag = true;
								filt += claa.comp[c].name;
							}
						}
						filt += ',';
						break;
						
					case "File":
						filt += claa.name+':file,';
						break;
					}
				}
				
				if(filt != "") te += 'filter="'+filt.substr(0,filt.length-1)+'" ';

				if(so.type == "Population"){
					te += add_obs_model(spec);
				}
			}
		}
		break;
		
	case "Pop. Trans.":
		{
			let sp = model.species[p];
		
			let cl = find(sp.cla,"name",spec.cl_drop.te);
			let claa = sp.cla[cl];
			let tra = claa.tra;
		
			let trz = spec.filter.tra;

			if(tra.length != trz.length){ 
				model.warn.push({mess:so.type, mess2:"The transition filter does not agree with the model.", warn_type:"SourceProb", p:p, ind:index});
			}
			else{
				if(spec.filter.cla.length != sp.ncla){
					model.warn.push({mess:so.type, mess2:"The transition population filter does not agree with with the model.", warn_type:"SourceProb", p:p, ind:index});
				}
				else{
					let fr = "", to = "";
					for(let j = 0; j < trz.length; j++){
						if(trz[j].check == true){
							file += '-'+tra[j].name;
							
							let fro, too;
							if(tra[j].i == "Source") fro = "+"; else fro = claa.comp[tra[j].i].name;
							if(tra[j].f == "Sink") too = "-"; else too = claa.comp[tra[j].f].name;
							
							if(fr != ""){ fr += ","; to += ",";}
							fr += fro; to += too;
						}
					}
				
					te += 'from="'+fr+'" to="'+to+'" ';
					
					let filt='';
					for(let cl2 = 0; cl2 < sp.ncla; cl2++){
						if(cl2 != cl){
							let clz = spec.filter.cla[cl2];
							let claa = sp.cla[cl2];
							switch(clz.radio.value){
							case "All": break;
							case "Comp": 
								filt += claa.name+':';
								let flag = false;
								for(let c = 0; c < claa.ncomp; c++){
									if(clz.comp[c].check == true){
										if(flag == true) filt += '|'; flag = true;
										filt += claa.comp[c].name;
									}
								}
								filt += ',';
								break;
								
							case "File":
								filt += claa.name+':file,';
								break;
								
							default: error("Problem with output option"); break;
							}
						}
					}
					
					if(filt != "") te += 'filter="'+filt.substr(0,filt.length-1)+'" ';
				
					te += add_obs_model(spec);
				}
			}
		}
		break;

	case "Genetic":
		te += 'root="SNP" ';
		break;

	default: break;	
	}

	file = get_unique_file(file,file_list,'.csv');
	
	output_table_file(dir+"/"+file,tab);
	
	te += 'file="'+file+'"'+endl;
	return te;
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
		if(find_in(file_list,file_not_same) == undefined) break;
		num++;
	}while(true);
	
	file_list.push(file_not_same);
	
	return file_not_same+end;
}


function add_obs_model(spec)
{
	switch(spec.obs_error.value){
	case "percent": return 'error="normal:'+spec.percent+'%" ';	
	case "sd": return 'error="normal:'+spec.sd+'" ';
	case "file": return 'error="normal:file" ';
	default: error("Option not recognised 123"); break;
	}
}


/// Outputs a table 
function output_table_file(file,tab)
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
	
	write_file(st,file);
}


/// Outputs a simple table 
function output_table_simp_file(file,head,row)
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
	
	write_file(st,file);
}
							

/// Checks the models is correctly specified
function output_check(type,exporting)
{
	// Checks to make sure initial population is set correctly
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		let source, pa;
		
		switch(type){
		case "sim": pa = "Simulation"; source = sp.sim_source; break;
		case "inf": pa = "Inference"; source = sp.inf_source; break;
		default: error("Option not recognised 125"); break;
		}
		
		let ninitpopprior = 0, ninitpop = 0, naddind = 0;
		
		for(let i = 0; i < source.length; i++){
			let so = source[i];
			
			switch(so.type){
			case "Init. Pop. Prior": ninitpopprior++; break;
			case "Init. Pop.": ninitpop++; break;
			case "Add Ind.": naddind++; break;
			}
		}

		switch(type){
		case "sim":
			if(ninitpopprior > 0){
				model.warn.push({mess:"Initial population prior", mess2:"The initial population prior should not be set for species '"+sp.name+"'", warn_type:"SimPopulationProb", p:p});
			}
			
			if(ninitpop > 0 && sp.fix_eff.length > 0){
				model.warn.push({mess:"Cannot add individuals with fixed effect", mess2:"The individuals from 'Init. Pop.' cannot be added because of individual fixed effects in the model", warn_type:"SimPopulationProb", p:p});
			}
	
			/*
			if(ninitpop == 0 && exporting == false){
				model.warn.push({mess:"Initial population", mess2:"The initial population should be set for species '"+sp.name+"'", warn_type:"SimPopulationProb", p:p});
			}
			*/
			
			if(ninitpop > 1){
				model.warn.push({mess:"Initial population", mess2:"Only one initial population should be set for species '"+sp.name+"'", warn_type:"SimPopulationProb", p:p});
			}
			break;
		
		case "inf":
			if(ninitpopprior > 1){
				model.warn.push({mess:"Initial population prior", mess2:"Only one initial population prior should be set for species '"+sp.name+"'", warn_type:"InfPopulationProb", p:p});
			}
			
			if(ninitpop > 1){
				model.warn.push({mess:"Initial population", mess2:"Only one initial population should be set for species '"+sp.name+"'", warn_type:"InfPopulationProb", p:p});
			}

			if(ninitpopprior == 0 && ninitpop == 0 && naddind == 0 && exporting == false){
				model.warn.push({mess:"No individuals", mess2:"Either 'Init. Pop.', 'Init. Pop. Prior' or 'Add Ind' must be set for species '"+sp.name+"'", warn_type:"InfPopulationProb", p:p});
			}
			
			if(ninitpopprior == 1 && ninitpop == 1){
				model.warn.push({mess:"Initial population", mess2:"'Init. Pop.' and 'Init. Pop. Prior' cannot both be set for species '"+sp.name+"'", warn_type:"InfPopulationProb", p:p});
			}
			
			if(ninitpop > 0 && sp.fix_eff.length > 0){
				model.warn.push({mess:"Cannot add individuals with fixed effect", mess2:"The individuals from 'Init. Pop.' cannot be added because of individual fixed effects in the model", warn_type:"InfPopulationProb", p:p});
			}
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
								model.warn.push({mess:"Repeated classification", mess2:"The classification '"+claa.name+"' is repeated twice in '"+sp.name+"'", warn_type:"ModelClass", p:p, cl:cl});
								return;
							}
							else{										
								if(claa.index != claa2.index){
									model.warn.push({mess:"Different indices", mess2:"The classification '"+claa.name+"' on species '"+sp.name+"' and '"+sp2.name+"' have different indices.", warn_type:"ModelClass", p:p, cl:cl});
									return;
								}
								else{
									let vec1=[], vec2=[];
									for(let c = 0; c < claa.ncomp; c++) vec1.push(claa.comp[c].name);
									for(let c = 0; c < claa2.ncomp; c++) vec2.push(claa2.comp[c].name);
									
									if(equal_vec(vec1,vec2) == false){
										model.warn.push({mess:"Compartments different", mess2:"For classification '"+claa.name+"' the compartments do not agree between species '"+sp.name+"' and '"+sp2.name+"'", warn_type:"ModelClass", p:p, cl:cl});
										return;
									}
								}
							}
						}
						else{
							if(claa.index == claa2.index){
								if(p == p2){
									model.warn.push({mess:"Identical indices", mess2:"In species '"+sp.name+"' the classifications '"+claa.name+"' and '"+claa2.name+"' should not share the same index '"+claa.index+"'", warn_type:"ModelClass", p:p, cl:cl});
									return;
								}
								else{
									model.warn.push({mess:"Identical indices", mess2:"Species '"+sp.name+" classifications '"+claa.name+"' and species '"+sp2.name+"' classification '"+claa2.name+"' should not share the same index '"+claa.index+"'", warn_type:"ModelClass", p:p, cl:cl});
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
	return st;
}

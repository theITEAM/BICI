"use strict";
// Function related to worker input and output

/// Adds an object to be saved
function save_add(te,ob)
{
	return te+":"+JSON.stringify(ob)+"\n";
}


/// Creates a new model
function new_model()                                          
{
	model = new Model();	 
	
	sim_result = {};    
	inf_result = {}; 
	ppc_result = {}; 

	model.start_new();
	
	post({ type:input.type, model:strip_heavy(model), info:input.info});
}


/// Sets up filters 
function initialise_filters_setup()
{
	if(sim_result.on == true) initialise_plot_filters_setup(sim_result,model.sim_res);
	if(inf_result.on == true) initialise_plot_filters_setup(inf_result,model.inf_res);
	if(ppc_result.on == true) initialise_plot_filters_setup(ppc_result,model.ppc_res);
}


/// Initialises filters for output plots 
function initialise_filters()
{
	if(sim_result.on == true) initialise_plot_filters(sim_result,model.sim_res);
	if(inf_result.on == true) initialise_plot_filters(inf_result,model.inf_res);
	if(ppc_result.on == true) initialise_plot_filters(ppc_result,model.ppc_res);
}


/// Loads up a file
function load_fi(file,file_type,heading,format,op)
{
	var fileReader = new FileReaderSync();
	let te = fileReader.readAsText(file, "UTF-8");

	te = te.replace(/\r/g, "");

	let error = false;
		
	switch(file_type){
	case "Observation_File": 
		break;
	
	case "GEOJSON_annotation":
		{
			let map = load_annotation_map(te);	
		
			post({type:"Load map", source:te, feature:map.feature, box:map.box, info:op});
		}
		break;
		
	case "GEOJSON_compartment":
		load_compartment_map(te);
		break;
	
	case "BICI_file": 
		import_file(te,true); 
		break;
		
	case "Import file": 
		import_file(te,true); 
		break;
	
	case "A matrix":
		{
			let sep = "comma";
			if(format == "tsv") sep = "tab"; 
			let tab = load_table(te,true,sep,file.name);
		
			if(typeof tab == 'string'){ alertp(tab); return;}
			else{
				A_matrix_loaded(tab,op.p,op.i);
				post({species:strip_heavy(model.species), info:file_type});
			}
		}
		break;
		
	case "Data file":
		{
			let head = true;
			let sep = "comma";
			if(heading == "No") head = false;
			if(format == "tsv") sep = "tab"; 

			let tab = load_table(te,head,sep,file.name);
	
			if(typeof tab == 'string'){
				alertp(tab);
				return;
			}
			else{
				post({ type:input.type, info:file_type, tab:tab});
			}
		}
		break;
		
	default: alertp("File type not recognised"); break;
	}
}


/// Loads vector X from table
function inf_eff_load_X(so)
{
	let i = so.info.i;
	let X = model.species[so.info.p].fix_eff[i].X_vector;
	
	X.loaded = true;
	X.X_value = [];
	X.ind_list = [];
	let tab = so.table;
	for(let j = 0; j < tab.nrow; j++){
		X.ind_list.push(tab.ele[j][0]);
		X.X_value.push(Number(tab.ele[j][1]));
	}
}


/// The function is called when the A pedigree is loaded
function add_A_pedigree(tab,p,i)
{
	let sp = model.species[p];
	let Amat = sp.ind_eff_group[i].A_matrix;

	Amat.pedigree = true;
	Amat.loaded = true;
	Amat.ind_list = [];
	Amat.sire_list = [];
	Amat.dam_list = [];
	
	for(let r = 0; r < tab.nrow; r++){
		Amat.ind_list[r] = tab.ele[r][0];
		Amat.sire_list[r] = tab.ele[r][1];
		Amat.dam_list[r] = tab.ele[r][2];
	}
}


/// The function is called when the A matrix is loaded
function A_matrix_loaded(tab,p,i)
{
	if(tab.ncol == 0){ alertp("The data table does not have any columns"); return;}
	if(tab.nrow == 0){ alertp("The data table does not have any rows"); return;}
	if(tab.nrow != tab.ncol){
		alertp("The data table must be square and have headings which give the IDs of the individuals"); 
		return;
	}
	
	for(let r = 0; r < tab.nrow; r++){
		for(let c = 0; c < tab.ncol; c++){
			if(isNaN(tab.ele[r][c])){
				alertp("The element '"+tab.ele[r][c]+"' is not a number (col "+(c+1)+", row "+(r+2)+")"); 
				return;
			}
		}
	}
	
	let sp = model.species[p];
	
	let Amat = sp.ind_eff_group[i].A_matrix;
	Amat.ind_list = tab.heading;
	Amat.loaded = true;
	Amat.pedigree = false;
	
	Amat.A_value = [];
	for(let r = 0; r < tab.nrow; r++){
		Amat.A_value[r] = [];
		for(let c = 0; c < tab.ncol; c++){
			Amat.A_value[r][c] = Number(tab.ele[r][c]);
		}
	}
}
				


/// Loads a table from text
function load_table(te,head,sep,filename)
{
	let file_simp = cut_path(filename);
	
	te.replace(/\r/g,"");
	let lines = te.split("\n");
	
	while(lines.length > 0 && (begin_str(lines[0].trim(),"#") || lines[0].trim() == "")){
		lines.splice(0,1);
	}
	
	let ele = [];
	for(let r = 0; r < lines.length; r++){
		let st = lines[r];
		if(st.length > 0){
			let li=[];
			
			switch(sep){
			case "tab":
				li = st.split("\t");
				break;
			
			case "comma":
				{
					let quote = "";
					let ist = 0;
					for(let i = 0; i < st.length; i++){
						let ch = st.substr(i,1);
						switch(ch){
						case "\"":
							if(quote == "\"") quote = "";
							else{
								if(quote == "") quote = "\"";
							}
							break;
						
						case ",":
							if(quote == ""){
								li.push(st.substr(ist,i-ist));
								ist = i+1;
							}
							break;
						}
					}
					
					if(ist < st.length) li.push(st.substr(ist,st.length-ist));
				}
				break;

			default: error("Option not recognised 62"); break;
			}
		
			if(head == true && ele.length > 0 && ele[0].length != li.length){
				return in_file_text(file_simp)+" there are "+ele[0].length+" headings but line "+(r+1)+" ("+trunc(st,30)+") contains "+li.length+" columns.";
			}

			ele.push(li);
		}
	}
	
	// Removes spaces and quotation marks
	for(let r = 0; r < ele.length; r++){
		for(let c = 0; c < ele[r].length; c++){ 	
			let te = ele[r][c].trim();
			if(te.length >= 2){
				if((te.substr(0,1) == "'" && te.substr(te.length-1,1) == "'") ||
				   (te.substr(0,1) == "\"" && te.substr(te.length-1,1) == "\"")){
					te = te.substr(1,te.length-2).trim();
				}						 
			}
			ele[r][c] = te;
		}
	}
	
	if(ele.length == 0) return in_file_text(file_simp)+" no content was found";
		
	let ncol = ele[0].length;
	
	for(let r = 0; r < ele.length; r++){
		if(ele[r].length != ncol){
			return in_file_text(file_simp)+" not all the columns have the same size (e.g. row "+(r+1)+" has "+ele[r].length+" columns).";
		}
	}
	
	let heading=[];
	
	if(head == true){
		heading = ele[0];
		ele.splice(0,1);
	}
	else{
		for(let c = 0; c < ncol; c++) heading.push("Column "+(c+1));
	}
	
	return { filename:file_simp, heading:heading, ele:ele, ncol:ncol, nrow:ele.length};
}
 
 
/// Reads parameter samples from a file 
function read_param_samples_file(chain,file,result)
{
	let te= file.te;

	let warn = "Problem loading file '"+file.name+"'";
	if(file.name == "inline") warn = "Problem loading parameter sample";
	
	read_param_samples(chain,te,result,warn)
}


/// Reads in a parameter samples from a file
function read_param_samples(chain,te,result,warn)
{
	let lines = te.split('\n');
	
	let spl = comma_split(lines[0].trim());
	
	if(spl[0] != "State"){ alert_param_sample(warn,-1); return;}

	let J = spl.length;

	let hash_ref = result.hash_ref;
	
	let ref = [];
	ref.push({});
	for(let i = 1; i < J; i++){
		let name = remove_quote(spl[i]).replace(/->/g,"→");
		name = name.replace(/\|/g,",");
	
		let ref_add = hash_ref.find(name)
		if(ref_add == undefined){
			if(find_in(like_name,name) != undefined){
				let th = result.param.length;
				ref_add = {th:th, index:[]};
				
				result.param.push({name:name, full_name:name, dep:[], time_dep:false, variety: "likelihood", selop:"like"});
				
			}
			else{
				if(find_in(trans_tree_name,name) != undefined){
					let th = result.param.length;
					ref_add = {th:th, index:[]};
					result.param.push({name:name, full_name:name, dep:[], time_dep:false, variety: "trans_tree", selop:"tree"});
				}
				else{
					let pp = get_param_prop(name);
			
					let th = find(result.param,"name",pp.name);
			
					if(th == undefined && begin(name,"N^")){
						for(let p = 0; p < result.species.length; p++){
							if(name == "N^"+result.species[p].name){
								th = result.param.length;
								let par = {name:pp.name, full_name:"<e>"+pp.name+"</e>", dep:[], time_dep:false, selop:"ic"};
			
								result.param.push(par);
								break;
							}			
							
							if(name == "N^"+result.species[p].name || name == "N^"+result.species[p].name+"-total"){
								th = result.param.length;
								let par = {name:pp.name, full_name:"<e>"+pp.name+"</e>", dep:[], time_dep:false, selop:"indpop"};
			
								result.param.push(par);
								break;
							}			
						}
					}
									
					if(th == undefined){
						// Adds parameters for initial condition and factor parameters
						let name = pp.name;
						let fl = false;
						let time_dep = false;
						if(name == "N^init") fl = true;
						if(name == "f^init") fl = true;
					
						if(fl == true){
							th = result.param.length;
						
							let dep = copy(pp.dep);
							
							let imax = dep.length; if(time_dep) imax--;
							for(let i = 0; i < imax; i++){
								let val = dep[i]; 
						
								dep[i] = "";
								for(let p = 0; p < result.species.length; p++){
									let sp = result.species[p];
									for(let cl = 0; cl < sp.ncla; cl++){
										let claa = sp.cla[cl];
										for(let c = 0; c < claa.comp.length; c++){
											let co = claa.comp[c];
											if(val == co.name) dep[i] = claa.index;
										}
									}									
								}
								if(dep[i] == "") alert_param_sample(warn,-40); 
							}
							
							if(time_dep) dep[imax] = "t";
						
							let full = pp.name+"_";
							for(let i = 0; i < imax; i++){
								if(i != 0) full += ",";
								full += dep[i];
							}
							if(time_dep) full += "(t)";
							
							let par = {name:pp.name, full_name:"<e>"+full+"</e>", dep:dep, time_dep:time_dep};
				
							/*
							if(false){                             // Turned off because no time variables added 
								if(time_dep){
									let knot;
									for(let p = 0; p < result.species.length; p++){
										let ppc_so = result.species[p].ppc_source;
										for(let k = 0; k < ppc_so.length; k++){
											let ppcs = ppc_so[k];
											if(ppcs.type == "Parameter Mult."){
												let na = "f~"+ppcs.spec.full_name;
												if(!(na.length > 3 && na.substr(na.length-3,3) == "(t)")) na += "(t)";
												if(na == full){
													knot = ppcs.spec.knot_times.split(",");
													break;
												}
											}
										}
									}
								
									if(knot == undefined) alert_param_sample(warn,-5); 
									par.spline = {on:true, smooth:{check:false}, time_dep:true, knot:knot};
								}
							}
							*/
						
							par.list = par_find_list(par);
							par.comb_list = generate_comb_list(par.list);
							result.param.push(par);
						}
					}
				
					if(th == undefined){
						alert_param_sample(warn,-2);
						return;
					}
					
					let par = result.param[th];
					if(pp.dep.length != par.dep.length){ alert_param_sample(warn,-3); return;}
					
					let ind = [];
					for(let j = 0; j < par.dep.length; j++){
						ind[j] = find_in(par.list[j],pp.dep[j]);	
						if(ind[j] == undefined){ alert_param_sample(warn,41); return;}
					}
					
					//if(par.output == false){ alert_param_sample(warn,5); return;}
					
					ref_add = {th:th, index:ind};
				}
			}
			hash_ref.add(name,ref_add);
		}
		ref.push(ref_add);
	}

	for(let i = 1; i < ref.length; i++){
		result.param[ref[i].th].output = true;
	}

	let temp = [];
	
	for(let th = 0; th < result.param.length; th++){
		let par = result.param[th];

		let value; 
		if(par.output == false) value = "No output";
		else{
			if(par.variety == "const"){
				value = "const";
			}
			else{
				if(par.dep.length == 0){
					value = "unset";
				}
				else{
					if(par.type == "derive_param" && par.time_dep) value = "derived";
					else{
						if(par.output){
							let list = par_find_list(par);
							value = par_find_template(list);
						}
					}
				}
			}
		}
		temp[th] = value;
	}
	
	for(let li = 1; li < lines.length; li++){
		if(lines[li] != ""){
			let spl = lines[li].split(",");
		
			if(spl.length != J){ alert_param_sample(warn,-4); return;}
			let param = copy(temp);
	
			let sa;
			
			if(chain == "gen-plot"){
				let par_mean = copy(temp);
				let par_CImin = copy(temp);
				let par_CImax = copy(temp);
				
				for(let i = 1; i < J; i++){
					let div = spl[i].split("|");
					if(div.length != 3){ alert_param_sample(warn,-5); return;}
					
					let mean = Number(div[0]); if(isNaN(mean)){ alert_param_sample(warn,-5); return;}
					let CImin = Number(div[1]); if(isNaN(CImin)){ alert_param_sample(warn,-5); return;}
					let CImax = Number(div[2]); if(isNaN(CImax)){ alert_param_sample(warn,-5); return;}
					
					let th = ref[i].th;
					if(ref[i].index.length == 0){
						par_mean[th] = mean;
						par_CImin[th] = CImin;
						par_CImax[th] = CImax;
					}
					else{
						set_element(par_mean[th],ref[i].index,mean);
						set_element(par_CImin[th],ref[i].index,CImin);
						set_element(par_CImax[th],ref[i].index,CImax);
					}
				}

				sa = { num:Number(spl[0]), mean:par_mean, CImin:par_CImin, CImax:par_CImax};
				result.generation.push(sa);
			}
			else{
				for(let i = 1; i < J; i++){
					let num = Number(spl[i]); if(isNaN(num)){ alert_param_sample(warn,-5); return;}
					
					let th = ref[i].th;
					if(ref[i].index.length == 0) param[th] = num;
					else{
						set_element(param[th],ref[i].index,num);
					}
				}
			
				sa = { chain:chain, num:Number(spl[0]), param:param};
				
				result.par_sample.push(sa);
			}
			
			result.param_memory += obj_memory(sa);
			if(result.param_memory > mem_param_sample_max*bytes_in_GB) out_of_param_memory(result);	
		}
	}
	
	// Checks no values for distance matrix or reparameterisation equations
	for(let th = 0; th < result.param.length; th++){
		let par = result.param[th];
		
		if(par.dist_mat || (par.variety == "reparam" && par.reparam_eqn_on)){
			if(par.comb_list.length != 0) error("Param error");
			if(par.value.length != 0) error("Param error");
			if(par.prior_split.length != 0) error("Param error");
		}
	}
	
	//prr("res"); prr(result);
	
	result.on = true;
}


/// Reads in state samples from a file
function read_state_samples_file(chain,file,result)
{	
	let te = file.te;
	
	let warn = "Problem loading file '"+file.name+"'";
	if(file.name == "inline") warn = "Problem loading state sample";
	
	read_state_samples(chain,te,result,warn);
}
	

/// Reads in state samples from text
function read_state_samples(chain,te,result,warn)
{
	let i = 0;
	while(i < te.length && te.substr(i,1) != "{") i++;
	let i_st = i;
	while(i < te.length && te.substr(i,1) != "}") i++;
	if(i == te.length){
		alert_sample("Error loading sample",1);
		return;
	}
	
	load_state_sample_header(te.substr(i_st,i),result,warn);
	i++;
	
	te = te.substr(i);

	let sec = te.split("<<");

	for(let s = 1; s < sec.length; s++){
		read_state_sample("<<"+sec[s],chain,result,warn);
	}	
}


/// Loads up information from the state sample header (time points and individual key)
function load_state_sample_header(te,result,warn)
{
	let spl = te.split("\n");
	let ind_key=[];
	for(let j = 0; j < spl.length; j++){
		let st = spl[j].trim();
		if(st != ""){
			if(st.length > 10 && st.substr(0,10) == "timepoint "){
				let spl = st.split(" ");
				let spl2 = spl[1].split(":");
				
				if(spl2.length != 3) alert_sample(warn,29);
				let t_start = Number(spl2[0]); if(isNaN(t_start)) alert_sample(warn,30);
				let dt = Number(spl2[1]); if(isNaN(dt)) alert_sample(warn,31);
				let t_end = Number(spl2[2]); if(isNaN(t_end)) alert_sample(warn,32);
				
				let timepoint=[];
				for(let t = t_start; t <= t_end+TINY; t += dt){
					timepoint.push(t);
				}
				
				result.timepoint = timepoint;
			}
			else{
				let spl = st.split(":");
				if(spl.length == 2){
					let i = Number(spl[0]); if(isNaN(i)) alert_sample(warn,33);
					ind_key[i] = spl[1];
				}
			}
		}		
	}		

	result.ind_key = ind_key;
	if(result.timepoint == undefined) alert_sample(warn,29);
}


/// Creates a hash table for compartments and transitions 
function create_compartment_hash(result)
{
	if(hash_compgl != undefined) return;
		
	hash_compgl=[]; hash_tragl=[];
	for(let p = 0; p < result.species.length; p++){
		let sp = result.species[p];
	
		hash_compgl[p] = new Hash();
		for(let c = 0; c < sp.comp_gl.length; c++){
			hash_compgl[p].add(sp.comp_gl[c].name,c);
		}

		hash_tragl[p] = new Hash();
		for(let tr = 0; tr < sp.tra_gl.length; tr++){
			hash_tragl[p].add(sp.tra_gl[tr].name,tr);
		}
	}
}


/// Reads a state file from text
function read_state_sample(te,chain,result,warning)
{
	create_compartment_hash(result);
	
	let lines = te.split('\n');
	
	for(let i = 0; i < lines.length; i++) lines[i] = lines[i].trim();
			
	let sample = { chain:chain, param:[], species:[], phytree:false, inf_node:[]};
	
	for(let th = 0; th < result.param.length; th++) sample.param[th] = undefined;
	
	let t_start = Number(result.details.t_start);
	
	let mode = "none";
	
	let sample_num;
	
	let get_inf_from = [];
	
	let ind_key = result.ind_key;
	
	let timepoint = result.timepoint;
	
	let p;
	for(let i = 0; i < lines.length; i++){
		let li = lines[i];
	
		if(li != ""){
			let warn = warning+": Line "+(i+1)+"   "+li;
	
			if(li.substr(0,1) == "<"){
				if(li.substr(li.length-1,1) != ">") alert_sample(warn,100);
				li = li.substr(1,li.length-2);
			
				let spl = li.split(" ");
			
				switch(spl[0]){
				case "<STATE":
					{
						let spl2 = spl[1].split(">");
						sample.num = spl2[0];
					}
					break;
					
				case "TRANSDISTPROB":
					mode = "transdiag"; 
					break;
					
				case "PARAMETERS":
					mode = "param"; 
					break;
					
				case "SPECIES":
					{
						mode = "species";
						if(spl.length != 2) alert_sample(warn,2);
						let sp_name = remove_quote(spl[1]);
						
						p = find(result.species,"name",sp_name);
						if(p == undefined){
							alert_sample(warn,"Species name '"+sp_name+"' not found");
						}
						
						let sp = result.species[p];
						
						switch(sp.type){
						case "Population": 
							sample.species[p] = {type:"Population", cpop_init:[], transnum_tl:[]}; 
							break;
							
						case "Individual": 
							sample.species[p] = {type:"Individual", individual:[], cpop_init:[]};
							break;
							
						default: error("PROBLEM"); break;
						}
					
						if(true){
							sample.species[p].trans_hbin = [];
							for(let tr = 0; tr < sp.tra_gl.length; tr++){
								sample.species[p].trans_hbin[tr] = [];
							}
						}
						
						let dpop=[];
						for(let c = 0; c < sp.comp_gl.length; c++){
							dpop[c]=[];
							for(let ti = 0; ti < timepoint.length-1; ti++) dpop[c][ti] = 0;
						}	
						sample.species[p].dpop = dpop;
					}						
					break;
				
				case "INITIAL":
					if(p == undefined) alert_sample(warn,4);
					if(sample.species[p].type == "Individual") alert_sample(warn,5);
					mode = "initial";
					break;
					
				case "POPULATION":
					if(p == undefined) alert_sample(warn,4);
					if(sample.species[p].type == "Individual") alert_sample(warn,5);
					mode = "popchange";
					break;
					
				case "TRANSITIONS": 
					if(p == undefined) alert_sample(warn,6);
					if(sample.species[p].type == "Individual") alert_sample(warn,7);
					mode = "trans"; 
					break;
				
				case "INDIVIDUALS":
					if(p == undefined) alert_sample(warn,10);
					if(sample.species[p].type == "Population") alert_sample(warn,11);
					mode = "individual";
					break;
					
				case "DERIVED":
					mode = "derive";
					break;
					
				case "PHYLOTREE":
					mode = "phylo";
					break;
					
				default: alert_sample(warn,12); return;
				}
			}
			else{
				switch(mode){
				case "phylo":
					{
						let spl = comma_split(lines[i]);
						let spl2 = spl[0].split("|");
						if(spl2.length != 4) alert_sample(warn,129); 
						let fr = spl2[2];
						if(fr != "ENT" && fr != "OUT") fr = Number(fr);
						let num = spl2[3]; if(num == "NA") num = UNSET; 
						let ino = {ind:spl2[0], t:Number(spl2[1]), from:fr, num:Number(num), inf_ev:[]}
						
						for(let k = 1; k < spl.length; k++){
							let spl2 = spl[k].split("|");
							
							if(spl2.length != 4) alert_sample(warn,130); 
							let num = spl2[3]; if(num == "NA") num = UNSET; 
							let iev = {ty:spl2[0], i:Number(spl2[1]), t:Number(spl2[2]), num:Number(num)};
							ino.inf_ev.push(iev);
						}
						sample.inf_node.push(ino);
						sample.phytree = true;
					}
					break;
					
				case "param": case "derive":		
					i = get_param_value(i,sample,lines,result,warn,mode);
					break;
					
				case "initial":
					{
						let spl = comma_split(lines[i]);
						if(spl.length != 2) alert_sample(warn,21);
						
						let name = spl[0], pop = spl[1];
				
						if(!(name == "compartment" && pop == "population")){
							if(isNaN(pop)) alert_sample(warn,23);
						
							let sp = result.species[p];
							
							if(name.includes("_erlang")){
								let k = 1;
								while(k <= name.length-7){
									if(name.substr(k,7) == "_erlang"){
										let kk = k-1;
										while(kk > 0 && name.substr(kk,1) != "_") kk--;
										name = name.substr(0,kk)+name.substr(k+7);
										k = kk;
									}
									else k++;
								}
								let c = hash_compgl[p].find(name);
								if(c == undefined) alert_sample(warn,122);
							
								sample.species[p].cpop_init[c] += Number(pop);
							}
							else{
								let c = hash_compgl[p].find(name);
								if(c == undefined) alert_sample(warn,123);
							
								sample.species[p].cpop_init[c] = Number(pop);
							}
						}
					}
					break;
					
				case "popchange":
					{
						let spl = comma_split(lines[i]);
						if(spl.length != 3) alert_sample(warn,21);
						
						let tp = spl[0], name = spl[1], change = spl[2];
						
						if(!(tp == "tp" && name == "compartment" && change == "change")){
							let sp = result.species[p];
							
							let c = hash_compgl[p].find(name);
							if(c == undefined) alert_sample(warn,124);
							
							if(isNaN(change)) alert_sample(warn,23);
							if(isNaN(tp)) alert_sample(warn,23);
							
							let samp_sp = sample.species[p];
							let ti = Number(tp);
							samp_sp.dpop[c][ti] += Number(change);
							if(!samp_sp.add_rem_pop_change) samp_sp.add_rem_pop_change=[];
							if(!samp_sp.add_rem_pop_change[ti]) samp_sp.add_rem_pop_change[ti] = [];
							samp_sp.add_rem_pop_change[ti].push({c:c,dpop:Number(change)});
						}
					}
					break;
				
				case "transdiag":
					{
						let sp = result.species[p];
							
						let ssp = sample.species[p];
						
						let spl = lines[i].split(":");
						if(spl.length != 2) alert_sample(warn,131);
							
						let name = spl[0].replace(/->/g,"→");
						if(!name.includes("_erlang")){
							let trg = hash_tragl[p].find(name);
							if(trg == undefined){
								alert_sample(warn,132);
							}
					
							let vec=[];
							let spl2 = spl[1].split("|");
							if(spl2.length != H_BIN) alert_sample(warn,133);
							
							let hbin = ssp.trans_hbin[trg];
							for(let b = 0; b < H_BIN; b++){
								if(isNaN(spl2[b])) alert_sample(warn,134);
								hbin[b] = Number(spl2[b]);
							}
						}
					}
					break;
					
				case "trans":
					{
						let sp = result.species[p];
						let ssp = sample.species[p];
						
						let spl = lines[i].split(":");
						if(spl.length != 2) alert_sample(warn,135);
							
						let name = spl[0].replace(/->/g,"→");
						if(!name.includes("_erlang")){
							let trg = hash_tragl[p].find(name);
							if(trg == undefined) alert_sample(warn,127);
					
							let vec=[];
							let spl2 = spl[1].split(",");
							for(let k = 0; k < spl2.length; k++){
								let val = spl2[k];
								if(val.substr(0,1) == "Z"){
									let num = Number(val.substr(1));
									for(let j = 0; j < num; j++) vec.push(0);
								}
								else{
									let num = Number(val);
									vec.push(num);
								}
							}
							
							if(vec.length != timepoint.length-1){
								alert_sample(warn,128);
							}
							
							ssp.transnum_tl[trg] = get_timeline(vec,result);
						}
					}
					break;
					
				case "individual":
					{
						let sp = result.species[p];
						let ssp = sample.species[p];
			
						let spl = lines[i].split(",");
				
						if(spl.length > 2 && spl[0] == "name" && spl[1] == "source"){
							sp.ind_effect = [];
							for(let k = 2; k < spl.length-1; k++) sp.ind_effect.push(spl[k]);	
						}
						else{	
							if(spl.length != 3+sp.ind_effect.length) alert_sample(warn,30);
							
							let name = spl[0].trim();
							name = ind_key[name];
							let so_str = spl[1].trim();
							
							let cinit;
							switch(so_str){
							case "yes": cinit = SOURCE; break;
							case "no": cinit = OUT; break;
							default: error("So string prob"); break;
							}
						
							let c = cinit;
							
							let ie = [];
							for(let k = 0; k < sp.ind_effect.length; k++){
								let val = Number(spl[k+2]);
								if(isNaN(val)) alert_sample(warn,32);
								ie.push(val);
							}
							
							let ev=[];
							
							let ev_text = spl[spl.length-1].trim();
							let obs = true;
							if(ev_text == "" || ev_text == "unobserved"){	
								if(ev_text == "unobserved") obs = false;
							}
							else{							
								let col = ev_text.split(" ");		
								
								for(let e = 0; e < col.length; e++){
									let tspl = col[e].split(":");
									if(tspl.length != 2) alert_sample(warn,33);
									let t = Number(tspl[1]); if(isNaN(t)) alert_sample(warn,34);
									
									if(c == OUT){                    // Individual enters the system
										c = hash_compgl[p].find(tspl[0]);
										if(c == undefined) alert_sample(warn,35);
										
										if(t == t_start) cinit = c;
										else{
											ev.push({type:EV_ENTER, t:t, c:c});
										}
									}
									else{
										let tr = tspl[0];
										if(c == SOURCE){                // Individual enters from source
											tr = extract_infection_info(tr,get_inf_from,p,ev,ssp);
										
											tr = tr.replace(/->/g,"→");
										
											let trg = hash_tragl[p].find(tr);
											if(trg == undefined) alert_sample(warn,36);
											
											let trgl = sp.tra_gl[trg];
											if(trgl.i != c) alert_sample(warn,37);
											ev.push({type:EV_TRANS, t:t, trg:trg});
											c = trgl.f;
										}
										else{
											if(tr == "-"){               // Individual leaves the system
												ev.push({type:EV_LEAVE, t:t});
											}
											else{			
												if(tr.length > 6 && tr.substr(0,6) == "move->"){
													let nam = tr.substr(6);
													let co = find_comp_from_name(nam,sp.name);
													if(co.warn != undefined) alert_sample(warn,38);
													
													let cnew = update_c_comp(sp,c,co.cl,co.c);
													if(c != cnew){
														ev.push({type:EV_MOVE, t:t, cl:co.cl, ci:c, cf:cnew});
														c = cnew;
													}
												}
												else{
													tr = extract_infection_info(tr,get_inf_from,p,ev,ssp);
				
													// POTENTIALLY SLOW
													let k = find(sp.comp_gl[c].tra_gl,"name",tr);
													if(k == undefined) alert_sample(warn);
													let trg = sp.comp_gl[c].tra_gl[k].num;
													let trgl = sp.tra_gl[trg];
													if(trgl.i != c) alert_sample(warn,39);
													
													ev.push({type:EV_TRANS, t:t, trg:trg});
													c = trgl.f;
												}
											}
										}
									}
								}
							}
							
							let ind = {name:name, cinit:cinit, ev:ev, ie:ie, obs:obs};
							
							ssp.individual.push(ind);
						}
					}
					break;
					
				default: alert_sample(warn,40); return;
				}
			}
		}
	}
	
	generate_trans_tree(get_inf_from,result,sample);

	generate_cpop_init_from_ind(result,sample);

	generate_transnum_from_ind(result,sample);

	generate_cpop_from_transnum(result,sample);

	update_average(result,sample);
	
	remove_dpop(sample);
	
	result.state_memory += obj_memory(sample);
	if(result.state_memory > mem_state_sample_max*bytes_in_GB) out_of_state_memory(result);

	result.sample.push(sample);
	if(result.on != true) result.on = true;
}


/// Reads in  information about transition diagnostics
function read_trans_diag(file,result)
{
	create_compartment_hash(result);
	
	let warning = "Problem loading file '"+file.name+"'";
	if(file.name == "inline") warning = "Problem loading state sample";
	
	let line = file.te.split("\n");
	
	let T = result.timepoint.length-1;
					
	let p = 0;
	
	for(let i = 0; i < line.length; i++){
		let li = line[i];
		
		if(li != ""){
			let warn = warning+": Line "+(i+1)+"   "+li;
	
			if(li.substr(0,1) == "<"){
				if(li.substr(li.length-1,1) != ">") alert_sample(warn,100);
				li = li.substr(1,li.length-2);
			
				let spl = li.split(" ");
			
				switch(spl[0]){
				case "SPECIES":
					{
						if(spl.length != 2) alert_sample(warn,2);
						let sp_name = remove_quote(spl[1]);
						
						p = find(result.species,"name",sp_name);
						if(p == undefined) alert_sample(warn,3);
					}						
					break;
				
				default :
					alert_sample(warn,100);
					break;
				}
			}
			else{
				let spl = li.split(":");
				if(spl.length != 2) alert_sample(warn,129);
					
				let name = spl[0].replace(/->/g,"→");
				if(!name.includes("_erlang")){
					let trg = hash_tragl[p].find(name);
					if(trg == undefined) alert_sample(warn,130);
					
					let spl2 = spl[1].split(",");
					if(spl2.length != T) alert_sample(warn,28);
					let exp_num=[];
					for(let ti = 0; ti < T; ti++){
						if(isNaN(spl2[ti])) alert_sample(warn,28);
						exp_num[ti] = Number(spl2[ti]);
					}
					
					if(!result.species[p].exp_num) result.species[p].exp_num = [];
					result.species[p].exp_num[trg] = exp_num;					
				}
			}				
		}
	}
	
	//prr("res");
	//prr(result);
}


/// This is called if parameter memory has run out
function out_of_param_memory(result)
{
	delete_results(result);	
	throw({ type:"Alert Param Memory", siminf:result.siminf});
}
		
		
/// This is called if state sample memory has run out
function out_of_state_memory(result)
{
	delete_results(result);
	throw({ type:"Alert State Memory", siminf:result.siminf});
}
		

/// Initialises average values (used for fast plotting and saving memory)
function average_init(result)
{
	result.average = [];
	result.average[0] = { species:[], num:LARGE, numav:0};
}
	
	
/// Updates average values (used for fast plotting and saving memory)
function update_average(result,sample)
{
	let T = result.timepoint.length-1;

	for(let p = 0; p < result.species.length; p++){
		let ssp = sample.species[p];
		
		let transnum = reconstruct_timeline_vec(ssp.transnum_tl,result);
		
		let av = result.average[0];
		if(av.numav == 0){
			let transnum = reconstruct_timeline_vec(ssp.transnum_tl,result);
			let cpop = reconstruct_timeline_vec(ssp.cpop_tl,result);
			av.species[p] = {transnum:transnum, cpop:cpop};
		}
		else{
			let transnum = av.species[p].transnum;
			
			for(let i = 0; i < ssp.transnum_tl.length; i++){
				add_tl(transnum[i],ssp.transnum_tl[i],result);
			}
			
			let cpop = av.species[p].cpop;
			for(let i = 0; i < ssp.cpop_tl.length; i++){
				add_tl(cpop[i],ssp.cpop_tl[i],result);
			}
		}
	}
	
	result.average[0].numav++;
}


/// Removes dpop information (because no longer needed)
function remove_dpop(result)
{
	for(let p = 0; p < result.species.length; p++){
		delete result.species[p].dpop;
	}
}


/// Rescales average values (used for fast plotting and saving memory)
function average_rescale(result)
{
	let T = result.timepoint.length-1;
	
	let numav = result.average[0].numav;
	for(let p = 0; p < result.species.length; p++){
		let as = result.average[0].species[p];
		let transnum = as.transnum;
		as.transnum_tl = [];
		for(let k = 0; k < transnum.length; k++){
			for(let ti = 0; ti < T; ti++) transnum[k][ti] /= numav;
			as.transnum_tl[k] = get_timeline(transnum[k],result);
		}
		delete as.transnum;
		
		let cpop = as.cpop;
		as.cpop_tl = [];
		for(let k = 0; k < cpop.length; k++){
			for(let ti = 0; ti <= T; ti++) cpop[k][ti] /= numav;
			as.cpop_tl[k] = get_timeline(cpop[k],result);
		}
		delete as.cpop;
	}
}
	
	
/// Extracts any infection information
function extract_infection_info(tr,get_inf_from,p,ev,ssp)
{
	if(tr.substr(tr.length-1,1) == "]"){
		let spl = tr.split("[");
		if(spl.length != 2) alert_sample("Problem reading file");
		tr = spl[0];
		get_inf_from.push({ te:spl[1].substr(0,spl[1].length-1), p:p, i:ssp.individual.length, e:ev.length});
	}
	return tr;
}


/// Finds which individuals infect which other ones
function generate_trans_tree(get_inf_from,result,sample)
{
	let trans_tree = [];
	
	let fl = false;
	for(let p = 0; p < result.species.length; p++){
		if(result.species[p].trans_tree.check == true){ trans_tree[p] = true; fl = true;}
		else trans_tree[p] = false;
	}
	
	if(fl == false) return;
	
	// Creates a hash table of inidividuals
	let hash_vec = []; 
	for(let p = 0; p < sample.species.length; p++){
		hash_vec[p] = new Hash();
		let sp = sample.species[p];
		for(let i = 0; i < sp.individual.length; i++){
			hash_vec[p].add(sp.individual[i].name,i);
		}
	}
	
	for(let k = 0; k < get_inf_from.length; k++){
		let gif = get_inf_from[k];
		let te = gif.te;
		
		let p_from, i_from;
		
		if(te == "OUT_INF"){
			p_from = OUTSIDE_INF;
			i_from = OUTSIDE_INF;
		}
		else{
			let spl = te.split("|");
			if(spl.length > 2) error("Problem with split");
			p_from = gif.p;
			if(spl.length == 2){
				p_from = find(result.species,"name",spl[0]);			
				if(p_from == undefined) alert_import("Could not find species '"+spl[0]+"'");
				te = spl[1];
			}
			
			i_from = hash_vec[p_from].find(te);
			if(i_from == undefined) alert_import("Could not find individual '"+te+"'");
		}
		
		sample.species[gif.p].individual[gif.i].ev[gif.e].infection = {p:p_from, i:i_from};
	}
}
	

/// Removes path from filename
function cut_path(filename)
{
	let j = filename.length-1;
	while(j >= 0 && filename.substr(j,1) != "/" && filename.substr(j,1) != "\\") j--;
	return filename.substr(j+1);
}


/// Gets a subtable based on a series of column heading names
function get_subtable(tab,col_name)
{
	let col = [];
	for(let i = 0; i < col_name.length; i++){
		let c = find_string_in(tab.heading,col_name[i]);
		if(c == undefined){
			return {error:in_file_text(tab.filename)+" cannot find heading '"+col_name[i]+"'"};
		}
		col.push(c);
	}
	
	let ele = [];
	for(let r = 0; r < tab.nrow; r++){
		ele[r]=[];
		for(let i = 0; i < col.length; i++){
			ele[r][i] = tab.ele[r][col[i]];
		}
	}		
	
	return {error:"", filename:tab.filename, heading:col_name, nrow:tab.nrow, ncol:col.length, ele:ele};
}
 

/// Loads a file using a http request
function load_file_http(file,type)
{
	let root = "Examples\\";
	if(begin(file,root)) file = "..\\"+file;
	
	const xhr = new XMLHttpRequest();
	
	if(ver=="windows") file = file.replace(/\//g,"\\");
	
  xhr.open("GET",file,false); // synchronous request
	
	let fl = false;
	
	try{
		xhr.send();
	}catch(e){
		fl = true;
	}
	
	if(xhr.responseText == "") fl = true;

	if(fl == true && type != "import"){
		alertp("The file '"+file+"' does not exist or is empty.");
	}
	
	return xhr.responseText;
}


/// Strips file from path
function filename(te)
{
	let spl = te.split("\\");
	return spl[spl.length-1];
}



/// Works out the width of text in a table
function text_width_worker(te,si,bold)
{
	cv.font = "30px Times";
	if(bold == true) cv.font = "bold 30px Times";
	return si*cv.measureText(te).width/30;
}	


/// Reduces the size of a tensor returned from the worker
function reduce_size(info,par)
{
	let scale = ELE_REDUCE_FAC*Math.pow(ELEMENT_MAX/num_element(par),1.0/par.dep.length);
			
	let list_shrink=[];
	for(let d = 0; d < par.dep.length; d++){
		let li=[];
		let imax = Math.floor(scale*par.list[d].length);
		if(imax < 0) imax = 1;
		
		for(let i = 0; i < imax; i++){
			li[i] = par.list[d][i];
		}
		list_shrink.push(li);
	}

	let temp = par_find_template(list_shrink);
	let comb_list = generate_comb_list(list_shrink);
	
	let value = info.value;
	
	for(let i = 0; i < comb_list.length; i++){
		let ind = comb_list[i].index;
		set_element(temp,ind,get_element(value,ind));
	}
			
	info.value = temp;
	info.list = list_shrink;
	info.too_big = true;
}


/// Gets results from menu name
function get_result(name)
{
	switch(name){
	case "Simulation": return sim_result;
	case "Inference": return inf_result;
	case "Post. Simulation": return ppc_result;
	}
}


/// Extracts information from the input text
function extract_text_samples(siminf,type,result)
{
	let com = siminf+"-"+type;
	
	let nchain = 1; if(result.chains) nchain = result.chains.length;
	
	let lines = result.import_te.split("\n");
	let pro = process_lines(lines);
	
	let te="";
	for(let i = 0; i < pro.processed.length; i++){
		let command = pro.processed[i];
		if(command.type == com){
			let ch = 0;
			for(let k = 0; k < command.tags.length; k++){
				let tag = command.tags[k];
				if(tag.name == "chain") ch = Number(tag.value);
				if(tag.name == "file"){
					switch(type){
					case "state":
						if(nchain > 1){
							te += "CHAIN "+(ch+1)+endl;
						}
						te += tag.value.te;
						break;
						
					case "param":
						{
							let spl = tag.value.te.split("\n");
							
							if(te == "") te += spl[0];
							if(nchain > 1) te += ",Chain";
							te += endl;
							
							for(let i = 1; i < spl.length; i++){
								if(spl[i] != ""){
									te += spl[i];
									if(nchain > 1) te += ","+(ch+1);
									te += endl;
								}
							}								
						}
						break;
					}
				}
			}
		}
	}
	
	if(te == "") alert_help("Export problem","There was a problem exporting this data");
	
	te = add_escape_char(te);
	
	return te;
}

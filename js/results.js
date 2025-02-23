"use strict";
// Functions associated with loading up results from BICI


/// Finalises after the import file has been loaded
function results_finalise(result)
{
	average_rescale(result);
	
	delete result.hash_tl;
}


/// Adds information about the model into the results (such that results can be plotted)
function results_add_model(result,details,siminf)
{
	result.sample = [];
	result.par_sample = [];
	result.plot_average = false;
	
	result.state_memory = 0;               // Keeps track of state memory
	result.param_memory = 0;               // Keeps track of parameter memory
	
	result.tl_store=[];                    // Stores timelines
	result.hash_tl = new Hash();           // Hash table for timelines
	result.hash_tl.create_ran();

	result.siminf = siminf;

	result.import_te = import_te;

	result.details = copy(details);
	
	result.species = copy(model.species);
	
	average_init(result);

	global_comp_trans_init(result.species);
	
	result.param = copy(model.param);

	create_param_const(result,siminf);            // Creates constant values for non-variable parameters 

	set_genetic_data(result);
	
	if(siminf == "inf") set_ppc_checkbox(model);
}


/// Deletes results (because memory has run out
function delete_results(result)
{
	delete result.sample;
	delete result.par_sample;
	delete result.tl_store;
	delete result.hash_tl;
	result.on = false;
}


/// Sets potential checkboxes for ppc
function set_ppc_checkbox(param)
{
	let check_box_list = [];
	for(let th = 0; th < model.param.length; th++){
		let par = model.param[th];
		
		if(par.variety == "dist"){
			let name = remove_eq_quote(par.full_name);
			if(find(check_box_list,"name",name) == undefined){		
				check_box_list.push({name:name, name_simp:par.name, checkb:{check:false, param:true}});
			}
		}
	}

	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let j = 0; j < sp.ind_eff_group.length; j++){
			let ieg = sp.ind_eff_group[j];
			let name = "";
			for(let k = 0; k < ieg.ie_list.length; k++){
				if(k != 0) name += "-";
				name += ieg.ie_list[k].name;
			}
			
			if(find(check_box_list,"name",name) == undefined){
				check_box_list.push({name:name, name_simp:name, checkb:{check:false, param:true}});
			}
		}
	}

	// Makes use of the resample 
	let resample = model.ppc_details.resample;
	if(resample){
		for(let k = 0; k < resample.length; k++){
			let na = resample[k];
			let j = find(check_box_list,"name_simp",na);
			if(j == undefined) alert_import("In resample could not find '"+na+"'",model.ppc_details.line);
			else check_box_list[k].checkb.check = true;
		}
	}
	
	model.ppc_details.check_box_list = check_box_list;
}


/// If there is genetic data then the table elements are put into result
function set_genetic_data(result)
{
	if(result.siminf == "inf"){
		for(let p = 0; p < result.species.length; p++){
			let sp = result.species[p];
			for(let i = 0; i < sp.inf_source.length; i++){
				let is = sp.inf_source[i];
				if(is.type == "Genetic"){				
					result.genetic_data = copy(is.table.ele);
				}
			}
		}
	}
}


/// Creates constant values for non-variable parameters 
function create_param_const(result,siminf)
{
	// For simulation removes parameter associated with observation process
	if(siminf == "sim"){ 
		let th = 0;
		while(th < result.param.length){
			if(find_in(sim_param_not_needed,result.param[th].type) == undefined ||
				result.param[th].type =="derive_param"){
				th++;
			}
			else result.param.splice(th,1);
		}
	}
	
	for(let th = 0; th < result.param.length; th++){
		let par = result.param[th];
		
		par.kind = "variable";
	
		if(par.variety == "const" || siminf == "sim"){ 
			if(!(par.variety == "dist" && par.sim_sample.check == true) 
				&& par.type != "derive_param"
				&& par.variety != "reparam"
				){
				par.kind = "const";
			}
		}
	}		
}


/// Sends an error that the sample could not be loaded
function alert_sample(warn,code)
{
	pr(warn+" CODE: "+code);
	alert_import(warn+" CODE: "+code);
}


/// Sends an error that the sample could not be loaded
function alert_param_sample(warn,code)
{
	alert_import(warn+" CODE: "+code);
}


/// Comma splits text but allows for quotation marks
function comma_split(te)
{
	let vec = [];
	vec[0] = "";
	let quote = 0;
	let num = 0;
	for(let i = 0; i < te.length; i++){
		let cha = te.substr(i,1);
		switch(cha){
		case ',': 
			if(quote == 0){ num++; vec[num] = "";}
			else vec[num] += cha;
			break;
			
		case '"': 
			quote = 1-quote;
			break;
			
		default: 
			vec[num] += cha;
			break;
		}
	}

	return vec;
}


/// Calculate the burn-in time
function calculate_burnin(result,source)
{
	if(source.plot_filter.details.algorithm.value != "DA-MCMC"){ source.burnin = 0; return;}
	
	// Parameter burnin
	let cha = result.chains;
	let samp_max = [];
	for(let ch = 0; ch < cha.length; ch++) samp_max[ch] = 0;
	
	for(let i = 0; i < result.par_sample.length; i++){
		let samp = result.par_sample[i];
		let ch = find_in(cha,samp.chain);
		if(ch == undefined) error("Should not be undefined");
		if(samp.num > samp_max[ch]) samp_max[ch] = samp.num;
	}
	
	let max = 0; for(let ch = 0; ch < cha.length; ch++) if(samp_max[ch] > max) max = samp_max[ch];
	source.samp_max = max;
	source.burnin = Math.floor(source.burnin_factor*source.samp_max);
}


/// Gets parameter value from line in the file
function get_param_value(i,source,lines,result,warn,mode)
{
	let param = result.param;

	let spl = comma_split(lines[i]);
	let name = remove_quote(spl[0].replace(/->/g,"→"));
	
	let th = find(result.param,"name",name);
	if(th == undefined) alert_sample(warn,13);
	
	let par = result.param[th];
	
	let der_time_var = false;
	if(par.time_dep == true && mode == "derive"){
		der_time_var = true;
	}

	if(der_time_var) par.dep.pop();
			
	let value;
	switch(spl.length){
	case 1:                                          // Set 	
		{		
			let list = par_find_list(par);
				
			value = par_find_template(list);
			let comb_list = generate_comb_list(list);
			par.list = list;
			par.comb_list = comb_list;
		
			i++;
			let spl_head = comma_split(lines[i]);
			
			let dep = par.dep;
			
			let ndep = dep.length;
		
			if(spl_head.length != ndep+1) alert_sample(warn,14);
			
			for(let d = 0; d < ndep; d++){
				if(dep[d] != spl_head[d]) alert_sample(warn,15); 
			}

			if(spl_head[ndep] != "Value") alert_sample(warn,16); 
				
			i++;
			for(let j = 0; j < comb_list.length; j++){
				let spl_row = lines[i].split(",");
			
				if(spl_row.length != ndep+1) alert_sample(warn,17);
				let index=[];
				for(let d = 0; d < ndep; d++){
					let k = find_in(list[d], spl_row[d]);
					if(k == undefined) alert_sample(warn,18);
					index.push(k);
				}
				let val = spl_row[ndep];
				
				if(der_time_var){
					val = split_spline(val);
				}
				else{
					if(isNaN(val)) alert_sample(warn,100);
					val = Number(val);
				}

				set_element(value,index,val);		
				i++;
			}
			i--;
		}
		break;
		
	case 2:                                          // Directly set parameter value (univariate)
		value = spl[1];
		
		if(der_time_var){
			value = split_spline(value);
		}
		else{
			value = Number(value);
			if(isNaN(value)) alert_sample(warn,19);
		}
		break;
		
	default: alert_sample(warn,20); return;
	}
	
	if(der_time_var){
		par.dep.push("t");
		let tp = get_time_points(result.details);
		tp.pop();
		par.spline.knot = tp;
		
		let list = par_find_list(par);
		
		par.list = list;
		par.comb_list = generate_comb_list(list);
	}
	
	if(par.comb_list){
		for(let j = 0; j < par.comb_list.length; j++){
			if(get_element(value,par.comb_list[j].index) == undefined) alert_sample(warn,101);
		}
	}
			
	source.param[th] = value;
	
	return i;
}
	

/// Splits text "1|3.." into a numeric vector
function split_spline(val)
{
	val = val.split("|");
	for(let k = 0; k < val.length; k++){
		if(isNaN(val[k])) alert_sample(warn,101);
		val[k] = Number(val[k]);
	}
	return val;		
}


/// Updates the global compartment based on a classification compartment being set to a given value
function update_c_comp(sp,c,cl,c_comp)
{
	let cco = sp.comp_gl[c].cla_comp[cl];
	let cnew = c +(c_comp - cco)*sp.comp_mult[cl];
	return cnew;
}


/// Generates cpop_init from individuals
function generate_cpop_init_from_ind(result,sample)
{
	for(let p = 0; p < sample.species.length; p++){
		let sp = result.species[p];
		let ssp = sample.species[p];
		
		if(ssp.type == "Individual"){
			let cpop = [];
			for(let c = 0; c < sp.comp_gl.length; c++) cpop[c] = 0;
			
			for(let i = 0; i < ssp.individual.length; i++){
				let ind = ssp.individual[i];
				let c = ind.cinit; 
				if(c != OUT && c != SOURCE){
					if(c < 0 || c >= sp.comp_gl.length) error("Out of range"+c);
					cpop[c]++;
				}
			}
			ssp.cpop_init = cpop;
		}
	}
}


/// Generates transnum from individuals
function generate_transnum_from_ind(result,sample)
{
	for(let p = 0; p < sample.species.length; p++){
		let sp = result.species[p];
		let ssp = sample.species[p];
		
		if(ssp.type == "Individual"){
			let timepoint = result.timepoint;
			let T = timepoint.length-1;
		
			let transnum=[];
			for(let k = 0; k < sp.tra_gl.length; k++){
				transnum[k] = [];
				for(let ti = 0; ti < T; ti++) transnum[k][ti] = 0;
			}
			
			let dpop = ssp.dpop;
			
			let tim = result.timepoint;
			let dt = 1;
			if(tim.length > 1) dt = tim[1]-tim[0];
			let t_start = tim[0];
	
			for(let i = 0; i < ssp.individual.length; i++){
				let ind = ssp.individual[i];
				
				let c = ind.cinit;
				
				for(let e = 0; e < ind.ev.length; e++){
					let ev = ind.ev[e];
					let t = ev.t;
				
					let ti = Math.floor(ALMOST_ONE*(t-t_start)/dt);
					if(ti < 0 || ti >= T) error("Out of range");
				
					switch(ev.type){
					case EV_TRANS:
						let tra = sp.tra_gl[ev.trg];
						if(tra.i != c){
							error("c is not consistent");
						}
						c = tra.f;
						transnum[ev.trg][ti]++;
						break;
						
					case EV_ENTER:
						dpop[ev.c][ti]++;
						if(c != OUT) error("Should be out");
						c = ev.c;
						break;
						
					case EV_LEAVE:
						if(e != ind.ev.length-1) error("Leave should be last event");
						dpop[c][ti]--;
						break;
						
					case EV_MOVE:
						if(ev.ci != c) error("move c is not consistent");
						dpop[ev.ci][ti]--;
						dpop[ev.cf][ti]++;
						c = ev.cf;
						break;
						
					default: error("event type not recognised:"+ev.type); break;
					}
				}
			}
			
			sample.species[p].transnum_tl = [];
			for(let k = 0; k < sp.tra_gl.length; k++){
				sample.species[p].transnum_tl[k] = get_timeline(transnum[k],result);
			}
		}
	}
}


/// Reconstructs timelines from timeline references 
function reconstruct_timeline_vec(vec,result)
{
	let tl_store = result.tl_store;

	let mat=[]
	for(let i = 0; i < vec.length; i++){
		mat.push(reconstruct_timeline(vec[i],result));
	}
	
	return mat;
}


/// Reconstructs a single timeline
function reconstruct_timeline(ref,result)
{
 	let tl = result.tl_store[ref];
	let v = tl.v, n = tl.n;
	let vec=[];
	for(let i = 0; i < v.length; i++){
		let va = v[i], num = n[i];
		for(let k = 0; k < num; k++) vec.push(va);
	}
	return vec;
}


/// Adds a compactified timeline onto an existing vector
function add_tl(vec,ref,result,op)
{
	let tl = result.tl_store[ref];
	let v = tl.v, n = tl.n;

	let j = 0;
	for(let i = 0; i < v.length; i++){
		let va = v[i], num = n[i];
		if(op){ if(op == "positive" && va < 0) va = 0;}
		for(let k = 0; k < num; k++){ vec[j] += va; j++;}
	}
	if(j != vec.length) pr("vec not right length");
}


/// Gets a value from a timeline at a given time point
function get_val_tl(ti,ref,result,op)
{
	let tl = result.tl_store[ref];
	let v = tl.v, n = tl.n;

	let i = 0;
	let ti_sum = 0;
	while(i < v.length && ti_sum+n[i] < ti){
		ti_sum += n[i];
		i++;
	}
	if(i == v.length) error("Problem with get val ti");
	
	let va = v[i];
	if(op){ if(op == "positive" && va < 0) va = 0;}
	
	return va;
}


/// Adds a time vector to the timeline
function get_timeline(vec,result)
{
	// First compatfies 
	let comp=[];
	let val=vec[0], num=1;
	let v=[], n=[];
	for(let i = 1; i < vec.length; i++){
		let va = vec[i];
		if(va != val){ 
			v.push(val); n.push(num);
			val = va; num = 1;
		}
		else num++;
	}
	v.push(val); n.push(num);
	
	// Looks in hash table
	let hash_tl = result.hash_tl;
	let tl_store = result.tl_store;
	
	let H = hash_tl.create_H(v,n);
	
	let j = hash_tl.find(H);
	if(j != undefined){
		let tl = tl_store[j];
		if(!equal_vec(tl.v,v)){ 
			j = undefined; pr("ERROR vec not equal1");
		}
		if(!equal_vec(tl.n,n)){ j = undefined; pr("ERROR vec not equal2");}
	}
	
	if(j == undefined){
		j = tl_store.length;
		let tls = {v:v, n:n};
		tl_store.push(tls);
	
		result.state_memory += obj_memory(tls);
		if(result.state_memory > mem_state_sample_max*bytes_in_GB) out_of_state_memory(result);
	
		hash_tl.add(H,j);
	}

	return j;
}


/// Generates cpop using global transitions specified in transnum
function generate_cpop_from_transnum(result,sample)
{
	for(let p = 0; p < sample.species.length; p++){
		let sp = result.species[p];
		let ssp = sample.species[p];

		//let dpop = reconstruct_timeline_vec(ssp.dpop_tl,result);
		let dpop = ssp.dpop;
		let transnum = reconstruct_timeline_vec(ssp.transnum_tl,result);
		
		let cpop=[];
		for(let c = 0; c < sp.comp_gl.length; c++){
			cpop[c]=[];
			cpop[c].push(ssp.cpop_init[c]);
		}

		for(let ti = 0; ti < result.timepoint.length-1; ti++){
			let cp=[];
			for(let c = 0; c < sp.comp_gl.length; c++){
				cp[c] = cpop[c][ti] + dpop[c][ti];
			}
				
			for(let j = 0; j < sp.tra_gl.length; j++){
				let tra = sp.tra_gl[j];
				let num = transnum[j][ti];
				let i = tra.i, f = tra.f;
				if(i != SOURCE) cp[i] -= num;
				if(f != SINK) cp[f] += num;
			}
			
			for(let c = 0; c < sp.comp_gl.length; c++) cpop[c].push(cp[c]);
		}
		
		ssp.cpop_tl=[];
		for(let c = 0; c < sp.comp_gl.length; c++){
			ssp.cpop_tl[c] = get_timeline(cpop[c],result);
		}
	}
}


/// Initialises global compartments and transitions
function global_comp_trans_init(species)
{
	for(let p = 0; p < species.length; p++){
		let sp = species[p];
		let inf_cl = find(sp.cla,"name",sp.infection_cl.te);
		
		sp.comp_mult=[];
		
		let N = 1;                                     // Works out global conversion
		for(let cl = 0; cl < sp.ncla; cl++){
			sp.comp_mult.push(N);
			N *= sp.cla[cl].comp.length;
		}
		
		sp.N = N;
	
		sp.comp_gl = [];                               // Constructs global compartments
		
		for(let c = 0; c < N; c++){                
			let co = {cla_comp:[], tra_gl:[]};
			
			let na = "";
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				let cc = Math.floor(c/sp.comp_mult[cl])%claa.comp.length;
				co.cla_comp.push(cc);

				if(na != "") na += "|";
				na += claa.comp[cc].name;
				
				if(cl == inf_cl && claa.comp[cc].infected.check == true){
					co.infected = true;
				}
			}
			co.name = na;
				
			sp.comp_gl.push(co);
		}

		sp.tra_gl = [];
		for(let cl = 0; cl < sp.ncla; cl++){           // Constructs global transitions
			let cla = sp.cla[cl];
			for(let j = 0; j < cla.ntra; j++){
				let tr = cla.tra[j];

				for(let c = 0; c < N; c++){
					let co = sp.comp_gl[c];		

					let cc = co.cla_comp[cl];
					let ii, ff;

					if(tr.i == SOURCE){
						if(tr.f == cc){ ii = SOURCE; ff = c;}
					}
					else{
						if(tr.f == SINK){
							if(tr.i == cc){ ii = c; ff = SINK;}
						}
						else{
							if(tr.i == cc){
								let dif = Number(tr.f) - Number(tr.i);
								ii = c; ff = c + dif*sp.comp_mult[cl];
							}
						}
					}

					if(ii != undefined){
						if(ii != SOURCE && ii != SINK){
							sp.comp_gl[ii].tra_gl.push({name:tr.name.replace(/→/g,"->"), num:sp.tra_gl.length});
						}
						
						let name = "";
						for(let cl2 = 0; cl2 < sp.ncla; cl2++){
							let claa2 = sp.cla[cl2];
							if(name != "") name += "|";
							if(cl2 == cl) name += tr.name;
							else name += claa2.comp[co.cla_comp[cl2]].name;
						}
			
						sp.tra_gl.push({name:name, cl:cl, tr:j, i:ii, f:ff});
					}
				}
			}
		}
	}
}


/// Works out the variety of plot based on the view type
function plot_variety(view)
{
	switch(view){
	case "Graph": case "Graph (all)": case "Graph (split)": 
	case "Graph (lines)": case "Graph (CI)": case "Data": return "Line plot";
	case "Scatter": return "Scatter plot";
	case "Compartment": return "Comp plot";
	case "Density": return "Density plot";
	case "Individual": return "Individual plot";
	case "TransTree": return "TransTree plot";
	case "PhyloTree": return "PhyloTree plot";
	case "Histogram": return "Histogram plot";
	case "HistoAnim": return "HistoAnim plot";
	case "Matrix": return "Matrix plot";
	case "MatrixAnim": return "MatrixAnim plot";
	case "CompMatrix": return "CompMatrix plot";
	case "CompMatrixAnim": return "CompMatrixAnim plot";
	case "CompVector": return "CompVector plot";
	case "Stat Table": return "Stat table plot";
	case "No Graph": return "No graph plot";
	default: error("Error variety:"+view); break;
	}
}


/// Initialises filters used to plot the results
function initialise_plot_filters(result,source)
{
	source.on = true;
	source.plot_filter = {};
	let rpf = source.plot_filter;
	
	rpf.plot_average = result.plot_average;
	rpf.siminf = result.siminf;
	rpf.details = result.details;
	
	// Initialises distribution plots
	let pos_bin = [];
	pos_bin.push({te:10}); pos_bin.push({te:20}); pos_bin.push({te:50}); 
	pos_bin.push({te:100}); pos_bin.push({te:200}); pos_bin.push({te:500});
		
	let pos_h = [];
	pos_h.push({te:0.01}); pos_h.push({te:0.02}); pos_h.push({te:0.05}); 
	pos_h.push({te:0.1}); pos_h.push({te:0.2}); pos_h.push({te:0.5});
	
	rpf.dist_settings = { show_mean:{check:true}, show_prior:{check:true}, radio:{value:"kde"}, sel_h:{te:0.05}, pos_h:pos_h, sel_bin:{te:50}, pos_bin:pos_bin};
	
	rpf.sim_val = {check:true};
	
	rpf.scatter_settings = {show_eb:{check:true}};

	// Possibilities for chains 
	result.chains = find_chains(result);
	rpf.pos_chain = [];
	rpf.pos_chain.push({te:"All", i:"All"});
	if(result.chains.length > 1){
		for(let i = 0; i < result.chains.length; i++){
			rpf.pos_chain.push({te:result.chains[i], i:i});
		}
	}	
	rpf.sel_chain = copy(rpf.pos_chain[0]);	

	// Possibilities for chain distribution 
	rpf.pos_distchain = [];
	rpf.pos_distchain.push({te:"Combine", i:"Combine"});
	if(result.chains.length > 1){
		rpf.pos_distchain.push({te:"All", i:"All"});
		for(let i = 0; i < result.chains.length; i++){
			rpf.pos_distchain.push({te:"Chain "+result.chains[i], i:i});
		}
	}	
	rpf.sel_distchain = copy(rpf.pos_distchain[0]);	
	
	// Possibilies for parameter view
	rpf.pos_paramview = get_pos_paramview(result);
	
	if(rpf.pos_paramview.length > 0){
		rpf.sel_paramview = copy(rpf.pos_paramview[0]);
	}
	
	// Possiblities for parameter checkbox
	rpf.param_check = get_param_check(result);
	
	// Sets up radio button for axes 
	rpf.xaxis_radio = get_tot_param_radio(result);
	rpf.yaxis_radio = get_tot_param_radio(result);
	if(rpf.yaxis_radio.param_list.length > 0 && rpf.yaxis_radio.radio.value+1 < rpf.yaxis_radio.param_list.length) rpf.yaxis_radio.radio.value++;
	
	// Possibilies for individual view
	rpf.pos_indview = [];
	
	rpf.pos_indview.push({te:"Timeline"});
	
	{
		for(let p = 0; p < result.species.length; p++){
			if(result.species[p].trans_tree.check == true){
				rpf.pos_indview.push({te:"Trans. Tree"});
				break;
			}
		}
	}
	
	{
		if(result.sample.length > 0){
			let sa = result.sample[0];
			if(sa.phytree){ 				
				rpf.pos_indview.push({te:"Phylo. Tree"});
			}
		}
	}
	
	if(rpf.siminf == "inf" || rpf.siminf == "ppc"){
		rpf.sel_ie_view = [];
		
		let fl = false;
		for(let p = 0; p < result.species.length; p++){
			let list = [];
			let sp = result.species[p];
			for(let k = 0; k < sp.inf_source.length; k++){
				let so = sp.inf_source[k];
				
				if(so.type == "Ind. Eff."){
					list.push({name:so.spec.drop.te});
					fl = true;
				}
			}
			rpf.sel_ie_view.push({ radio: {value: 0}, list:list});
		}
		
		if(fl == true) rpf.pos_indview.push({te:"Ind. Eff."});
	}
					
	rpf.pos_indview.push({te:"Table"});
	
	rpf.sel_indview = copy(rpf.pos_indview[0]);
	
	// Possiblities for time-step
	rpf.pos_timestep = get_possible_timestep(result.details);

	if(rpf.pos_timestep.length > 0){
		rpf.sel_timestep = copy(rpf.pos_timestep[0]);
	}
	
	// Possibilities for classification
	rpf.species = [];
	for(let p = 0; p < result.species.length; p++){
		let sp = result.species[p];
		
		let pos_class = [];
		for(let cl = 0; cl < sp.ncla; cl++){
			pos_class.push({te:sp.cla[cl].name, siminf:rpf.siminf, p:p, cl:cl, check_filter:true});
		}

		let cla_red = [];
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			cla_red[cl] = { name:claa.name};
			if(claa.ncomp > COMP_FILTER_MAX) cla_red[cl].comp_too_large = true;
			else{ cla_red[cl].comp = claa.comp; cla_red[cl].ncomp = claa.ncomp;} 
			
			if(claa.ntra > TRANS_FILTER_MAX) cla_red[cl].tra_too_large = true;
			else{ cla_red[cl].tra = claa.tra; cla_red[cl].ntra = claa.ntra;} 
		}
		
		rpf.species[p] = { sel_class:copy(pos_class[0]), pos_class:pos_class, cla:cla_red, ncla:sp.ncla, name:sp.name, type:sp.type, gen_source:[], filter:[]};
	}
	
	// Possibilities for sample
	let pos_sample = [];
	if(result.sample.length > 1) pos_sample.push({te:"All", i:"All"});
	for(let i = 0; i < result.sample.length; i++){
		let te=""; if(result.chains.length > 1) te = "[Chain "+result.sample[i].chain+"] ";
		te += result.sample[i].num;
		pos_sample.push({te:te, i:i});
	}
	
	if(pos_sample.length == 0) error("No samples");
	else rpf.sel_sample = copy(pos_sample[0]);
	rpf.pos_sample = pos_sample;

	// Possibilities for sample in different chains
	rpf.pos_sample_chain = [];
	rpf.sel_sample_chain = [];
	for(let ch = 0; ch < result.chains.length; ch++){
		rpf.pos_sample_chain[ch] = [];
		rpf.pos_sample_chain[ch].push({te:"All", i:"All"});
		for(let i = 0; i < result.sample.length; i++){
			if(result.sample[i].chain == result.chains[ch]){
				rpf.pos_sample_chain[ch].push({te:result.sample[i].num, i:i});
			}
		}
		
		rpf.sel_sample_chain[ch] = copy(rpf.pos_sample_chain[ch][0]);
	}
	
	// Possibilities for population and poptrans data
	let pos_popdata = [];
	let pos_poptransdata = [];
	if(result.siminf == "inf" || result.siminf == "ppc"){
		for(let p = 0; p < result.species.length; p++){
			let sp = result.species[p];
			for(let i = 0; i < sp.inf_source.length; i++){
				let so = sp.inf_source[i];
			
				switch(so.type){
				case "Population": case "Pop. Trans.": 
					{
						let desc_base ;
						if(so.type == "Population") desc_base = "Pop. "; 
						else{
							let spl = so.desc.split(" ");
							desc_base = spl[0]+" ";
						}		
						
						let filt = so.spec.filter;				
						
						for(let cl = 0; cl < filt.cla.length; cl++){
							let claa = filt.cla[cl];
							switch(claa.radio.value){
							case "Comp":
								for(let c = 0; c < claa.comp.length; c++){
									if(claa.comp[c].check == true) desc_base += claa.comp[c].comp_name_store+" ";
								}
								break;
								
							case "ObsMod":
								for(let c = 0; c < claa.comp.length; c++){
									if(claa.comp[c].prob_eqn.te != "0") desc_base += claa.comp[c].comp_name_store+":"+claa.comp[c].prob_eqn.te+" ";
								}
								break;
							}
						}
						
						
						let strat = get_stratify_data(filt);
			
						for(let k = 0; k < strat.comb.length; k++){
							let ind = strat.comb[k].index;
							let desc = desc_base;
							let strat_filt=[];
							for(let m = 0; m < ind.length; m++){
								if(m != 0) desc += ",";
								desc += strat.list[m][ind[m]];
								strat_filt[strat.cl_list[m]] = ind[m];
							}
						
							if(so.type == "Population") pos_popdata.push({te:desc, p:p, i:i, strat_filt:strat_filt});
							else pos_poptransdata.push({te:desc, p:p, i:i, strat_filt:strat_filt});
						}
					}
					break;
				}
			}
		}
	}			
	
	if(pos_popdata.length > 0) rpf.sel_popdata = copy(pos_popdata[0]);
	rpf.pos_popdata = pos_popdata;
	
	if(pos_poptransdata.length > 0) rpf.sel_poptransdata = copy(pos_poptransdata[0]);
	rpf.pos_poptransdata = pos_poptransdata;
	
	// Possibilities for viewing populations
	let pos_view = [];
	if(pos_sample.length == 1){
		pos_view.push({te:"Graph"});
	}
	else{
		pos_view.push({te:"Graph (CI)"});
		pos_view.push({te:"Graph (lines)"});	
	}
	rpf.slice_time = result.details.t_start;
	pos_view.push({te:"Slice"});
	pos_view.push({te:"Compartment"});
	pos_view.push({te:"Density"});
	if(pos_popdata.length > 0) pos_view.push({te:"Data"});
	
	for(let k = 0; k < pos_view.length; k++) pos_view[k].bub_off = true;
	
	rpf.sel_view = copy(pos_view[0]);
	rpf.pos_view = pos_view;
	
	// Possibilities for viewing individual 
	
	let pos_ind_sel_view=[];
	pos_ind_sel_view.push({te:"Graph"});
	pos_ind_sel_view.push({te:"Compartment"});
	pos_ind_sel_view.push({te:"Density"});
	rpf.sel_ind_sel_view = copy(pos_ind_sel_view[0]);
	rpf.pos_ind_sel_view = pos_ind_sel_view;
	
	// Possibilities for viewing transitions
	let pos_trans_view = [];
	if(pos_sample.length == 1){
		pos_trans_view.push({te:"Graph"});
	}
	else{
		pos_trans_view.push({te:"Graph (CI)"});
		pos_trans_view.push({te:"Graph (lines)"});	
	}
	if(pos_poptransdata.length) pos_trans_view.push({te:"Data"});
	
	rpf.sel_trans_view = copy(pos_trans_view[0]);
	rpf.pos_trans_view = pos_trans_view;
	
	// Sets up potential groups of individuals
	
	for(let p = 0; p < result.species.length; p++){
		let pos = [];
		pos.push({te:"All"});
		let inf_so = result.species[p].inf_source;
		for(let i = 0; i < inf_so.length; i++){
			let is = inf_so[i];
			if(is.type == "Ind. Group"){
				let ele = is.table.ele;
				let name_list = [];
				for(let j = 0; j < ele.length; j++){
					name_list.push(ele[j][0]);
				}
				pos.push({te:is.spec.gname, name_list:name_list});
			}
		}
		
		rpf.species[p].pos_indgroup = pos;
		rpf.species[p].sel_indgroup = copy(pos[0]);
	}

	source.burnin_factor = 0.2;
	if(rpf.details.burnin_frac){
		source.burnin_factor = Number(rpf.details.burnin_frac)/100;
	}
	
	calculate_burnin(result,source);
	
	source.ind_flag = false;

	if(find(result.species,"type","Individual") != undefined) source.ind_flag = true;
	
	source.species_name = [];
	for(let p = 0; p < result.species.length; p++) source.species_name.push(result.species[p].name);
	
	source.param = result.param;
}


/// Works out how data is stratified
function get_stratify_data(filt)
{
	let cl_list = [];
	let list = [];
	for(let cl = 0; cl < filt.cla.length; cl++){
		let fcl = filt.cla[cl];
		if(fcl.radio.value == "File"){
			cl_list.push(cl);
			let comp_list = [];
			
			for(let c = 0; c < fcl.comp.length; c++){
				comp_list.push(fcl.comp[c].comp_name_store);
			}
			list.push(comp_list);
		}
	}
	
	let comb = generate_comb_list(list);

	return { cl_list:cl_list, comb:comb, list:list};
}


/// Gets possible views for a parameter
function get_par_pos_view(par,so)
{
	let pos_view = [];

	let name = par.full_name;
	
	if(par.spline && par.spline.on == true){ // Possible views for splines
		if(par.dep.length > 1){
			pos_view.push({te:"Graph (all)", title:"Spline profile "+name, type:"Timevary"});
			pos_view.push({te:"Graph (split)", title:"Spline profile "+name, type:"Timevary"});
		}
		else{
			pos_view.push({te:"Graph", title:"Spline profile "+name, type:"Timevary"});
		}
		
		let title = "Time variation in "+name;
		
		switch(par.dep.length){
		case 2:
			{
				let index = remove_prime(par.dep[0]);
				let pcl = find_p_cl_index(so,index);
				if(pcl != undefined){
					pos_view.push({te:"Compartment", title:title, type:"Timevary", p:pcl.p, cl:pcl.cl});
					pos_view.push({te:"Density", title:title, type:"Timevary", p:pcl.p, cl:pcl.cl});
					pos_view.push({te:"Histogram", title:title, type:"Timevary", p:pcl.p, cl:pcl.cl});
				}
			}
			break;
			
		case 3:
			{
				pos_view.push({te:"Matrix", title:title, type:"Timevary"});
				
				let index = remove_prime(par.dep[0]);
				if(index == remove_prime(par.dep[1])){
					let pcl = find_p_cl_index(so,index);
					if(pcl != undefined){
						pos_view.push({te:"Compartment", title:title, type:"Timevary", p:pcl.p, cl:pcl.cl});
					}
				}
			}
			break;
		}
	}
	else{
		if(par.dep.length == 1){                       // Possible views for vectors
			let title = "Vector "+name;
		
			pos_view.push({te:"Histogram", title:title, type:"Vector"});
		
			let index = remove_prime(par.dep[0]);
			let pcl = find_p_cl_index(so,index);
			if(pcl != undefined){
				pos_view.push({te:"Compartment", title:title, type:"Vector", p:pcl.p, cl:pcl.cl});
				pos_view.push({te:"Density", title:title, type:"Vector", p:pcl.p, cl:pcl.cl});
			}
		}
		
		if(par.dep.length == 2){                       // Possible views for matrices
			let title = "Matrix "+name;
			
			pos_view.push({te:"Matrix", title:title, type:"Matrix"});
		
			let index = remove_prime(par.dep[0]);
			if(index == remove_prime(par.dep[1])){
				let pcl = find_p_cl_index(so,index);
				if(pcl != undefined){
					pos_view.push({te:"Compartment", title:title, type:"Matrix", p:pcl.p, cl:pcl.cl});
				}
			}
		}
	}

	return pos_view;
}


/// Given an index finds which species and classification
function find_p_cl_index(so,index)
{
	for(let p = 0; p < so.species.length; p++){
		let sp = so.species[p];
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			if(claa.index == index) return {p:p, cl:cl};
		}						
	}
}

				
/// Finds a list of all chains
function find_chains(result)
{
	let list = [];
	for(let i = 0; i < result.par_sample.length; i++){
		let samp = result.par_sample[i];
		if(find_in(list,samp.chain) == undefined) list.push(samp.chain);
	}
	return list;
}
	
	
/// When chaning classification removes filters (if neccesary)
function change_remove_filter(pos)
{
	let rpf;
	switch(pos.siminf){
	case "sim":	rpf = model.sim_res.plot_filter; break;
	case "inf":	rpf = model.inf_res.plot_filter; break;
	case "ppc":	rpf = model.ppc_res.plot_filter; break;
	default: error("crf prob"); break;
	}

	let cl_now = rpf.species[pos.p].sel_class.cl;

	let filter = rpf.species[pos.p].filter;
	let j = 0;
	while(j < filter.length){
		if(filter[j].cl == pos.cl || filter[j].cl == cl_now){ filter.splice(j,1); close_bubble();}
		else j++;
	}
}
	

/// Gets a list of all possible parameter views
function get_pos_paramview(result)
{
	// Determines which parameters to output
	for(let th = 0; th < result.param.length; th++){
		let par = result.param[th];
		
		par.output = true;
		if(par.list && par.type != "derive_param"){
			let N = 1; for(let i = 0; i < par.list.length; i++) N *= par.list[i].length;
			if(N > result.details.param_output_max) par.output = false;
		}
	}	
	
	let total_param_list = [];
	
	let pos_paramview = [];
	
	{                                               // Univariate parameters
		let list = [];
		
		for(let th = 0; th < result.param.length; th++){
			let par = result.param[th];

			if(par.output == true && par.dep.length == 0 && par.selop == undefined){
				list.push({th:th, index:[], name:param_name_index(par,[])});
				if(par.kind != "const"){
					total_param_list.push({th:th, index:[], name:param_name_index(par,[])});
				}
			}
		}
		
		if(list.length > 0){
			pos_paramview.push({te:"Univariate", param:true, radio:{value:0, param:true}, list:list});	
		}
	}
	
	for(let th = 0; th < result.param.length; th++){
		let par = result.param[th];
		if(par.dep.length > 0 && par.type != "param factor"){
			let list = [];
			let list_split = [];
			let para = true;
			if(par.output == true){
				for(let j = 0; j < par.comb_list.length; j++){	
					let ind = par.comb_list[j].index;
					
					let name = param_name_index(par,ind);
					list.push({th:th, index:ind, name:name});
					
					if(par.time_dep){
						if(ind[ind.length-1] == 0){
							let ind2 = copy(ind); ind2.pop();
							let spl = name.split("(");
							list_split.push({th:th, index:ind2, name:spl[0]});
						}
					}
					
					if(par.kind != "const" && par.type != "derive_param"){
						total_param_list.push({th:th, index:ind, name:param_name_index(par,ind)});
					}
				}
			}
			else para = "too big";
			
			let ppv = {te:remove_eq_quote(par.full_name), param:para, radio:{value:0, param:true}, list:list, radio_split:{value:0, param:true}, list_split:list_split};
			
			pos_paramview.push(ppv);
		}
	}

	result.total_param_list = total_param_list;

	if(total_param_list.length > 1 && (result.siminf == "inf" || result.siminf == "ppc")){
		pos_paramview.push({te:"Scatter", param:true, radio:{value:0, param:true}, radioy:{value:1, param:true}, list:total_param_list});
		
		pos_paramview.push({te:"Correlation", param:true, radio:{value:0, param:true}, radioy:{value:1, param:true}, list:total_param_list});
	}
	
	let list_trans = [], list_like = [], list_ic = [], list_indpop = [];
	for(let th = 0; th < result.param.length; th++){
		let par = result.param[th];
		switch(par.selop){
		case "tree": list_trans.push({th:th, index:[], name:par.name}); break;
		case "like": list_like.push({th:th, index:[], name:par.name}); break;
		case "ic": list_ic.push({th:th, index:[], name:par.name}); break;
		case "indpop": list_indpop.push({th:th, index:[], name:par.name}); break;
		}
	}
	
	if(list_trans.length > 0){
		pos_paramview.push({te:"Trans. Tree", param:true, radio:{value:0, param:true},  list:list_trans});	
	}
	
	if(result.siminf == "inf" || result.siminf == "ppc"){
		if(list_like.length > 0){
			pos_paramview.push({te:"Likelihoods", param:true, radio:{value:0, param:true},  list:list_like});	
		}
	
		if(list_ic.length > 0){
			pos_paramview.push({te:"Init. Pop.", param:true, radio:{value:0, param:true}, list:list_ic});	
		}
	
		if(list_indpop.length > 0){
			pos_paramview.push({te:"Ind. Pop.", param:true, radio:{value:0, param:true}, list:list_indpop});	
		}
	}

	return pos_paramview;
}


/// Checkboxes which allow
function get_param_check(result)
{
	let param_list = [];
	
	for(let i = 0; i < result.total_param_list.length; i++){
		let tpl = result.total_param_list[i];
		if(!param_list[tpl.th]){
			if(tpl.index.length == 0){
				param_list[tpl.th] = { name:tpl.name, checkb:{check:true, param:true, update:true}, tpl_ref:i};
			}
			else{
				let par = result.param[tpl.th];
				param_list[tpl.th] = { name:remove_eq_quote(par.full_name), checkb:{check:false, param:true, update:true}, list:[]};
			}
		}
		
		if(tpl.index.length > 0){
			param_list[tpl.th].list.push({ name:tpl.name, checkb:{check:false, param:true, update:true}, tpl_ref:i}); 
		}
	}
	
	for(let th = 0; th < param_list.length; th++){
		let pl = param_list[th];
		if(pl && pl.list){
			if(pl.list.length > 4){
				pl.all_checkb = { check:false, list:pl.list, all_check:true};
			}
			else{ // Switchs on if 4 or fewer elements
				pl.checkb.check = true;
				for(let k = 0; k < pl.list.length; k++){
					pl.list[k].checkb.check = true;
				}
			}
		}
	}
	
	remove_undefined_param_list(param_list)
	
	return param_list;
}


/// Constructs a radio button for all parameters
function get_tot_param_radio(result)
{
	let param_list = [];
	for(let i = 0; i < result.total_param_list.length; i++){
		let tpl = result.total_param_list[i];
		if(!param_list[tpl.th]){
			if(tpl.index.length == 0){
				param_list[tpl.th] = { name:tpl.name, tpl_ref:i};
			}
			else{
				let par = result.param[tpl.th];
				param_list[tpl.th] = { name:remove_eq_quote(par.full_name), radio:{value:0, param:true}, list:[]};
			}
		}
		
		if(tpl.index.length > 0){
			param_list[tpl.th].list.push({ name:tpl.name, tpl_ref:i}); 
		}
	}
	
	remove_undefined_param_list(param_list);
	
	let param_radio = { radio:{value:0, param:true}, param_list:param_list};
	
	return param_radio;
}
 

/// Removes any undefined parameters from list
function remove_undefined_param_list(param_list)
{	
	let i = 0;
	while(i < param_list.length){
		if(!param_list[i]) param_list.splice(i,1);
		else i++;
	}
}
	
	
	
/// Gets possible timesteps based on dt
function get_possible_timestep(details)
{
	let dt = details.timestep;
	
	let range = Number(details.t_end) - Number(details.t_start);
	let dt_max = range;

	let pos_timestep = [];
	let pos = [1,2,4,7,10,20,40,70,100];

	if(dt < 1){
		let num = 1.0/dt;
		let d = num - Math.round(num);
		if(d*d < TINY){
			let j = 1;
			while(j < num){
				add_timestep(pos_timestep,dt,dt_max,j);
				j *= 2;
			}
			
			for(let i = 0; i < pos.length; i++){
				add_timestep(pos_timestep,dt,dt_max,num*pos[i]);
			}
		}
		else{
			for(let i = 0; i < pos.length; i++){
				add_timestep(pos_timestep,dt,dt_max,pos[i]);
			}
		}
	}
	else{
		for(let i = 0; i < pos.length; i++){
			add_timestep(pos_timestep,dt,dt_max,pos[i]);
		}
	}
	
	return pos_timestep;
}


/// Adds a timestep to the possibility
function add_timestep(pos_timestep,dt,dt_max,j)
{
	if(j*dt <= dt_max) pos_timestep.push({te:dt*j, i:j});
}

				
/// Adds bubble content which allows for a filter to be added
function filter_bubble(bu,cont)
{
	let filt=[];

	let rpf = bu.rpf;
	
	let p = model.get_p();
	
	let rpf2 = rpf.species[p];
	let cl = rpf2.sel_class.cl;
	
	for(let i = 0; i < rpf2.filter.length; i++) filt[rpf2.filter[i].cl] = "not";
	
	let sp = rpf.species[p];
	
	let pos=[];
	for(let cl2 = 0; cl2 < sp.ncla; cl2++){
		if(filt[cl2] != "not") pos.push({te:sp.cla[cl2].name, cl:cl2});	
	}
	
	if(pos.length == 0){
		cont.dx = 9;
		bubble_addtitle(cont,"Add Filter",{});
	
		bubble_addparagraph(cont,"There are no potential classifications which can be filtered.",0,cont.dx);
		add_end_button(cont,"OK","CloseBubble",{});	
	}
	else{
		cont.dx = 7;
		bubble_addtitle(cont,"Add Filter",{});
	
		for(let i = 0; i < pos.length; i++){
			cont.lay.add_button({te:pos[i].te, x:0, y:cont.y, dx:cont.dx, dy:dropdown_height, rpf:rpf, cl:pos[i].cl, ac:"AddFilter2", type:"PopFilt"});
			cont.y += dropdown_height+0.2;
		}
	
		add_end_button(cont,"OK","Done",{});
	}
}


/// Adds a filter to a population plot
function add_pop_filter(rpf,cl)
{
	let p = model.get_p();

	let sp = rpf.species[p];

	let rpf2 = rpf.species[p];
	
	let claa = sp.cla[cl];
		
	if(claa.comp_too_large){ alert_help("Too many compartments to filter"); return;}
	
	let comp_filt = [];
	let comp_name = [];
	for(let c = 0; c < claa.comp.length; c++){
		comp_filt.push({check:true});
		comp_name.push(sp.cla[cl].comp[c].name);
	}
	
	let tra_filt = [];
	let tra_name = [];
	for(let j = 0; j < claa.ntra; j++){
		let tr = claa.tra[j];
		tra_filt.push({check:true});
		tra_name.push(claa.tra[j].name);
	}
	
	let filter = {name:sp.cla[cl].name, p:p, cl:cl, fraction:{check:false, update_fraction:true, p:p, cl:cl}, all_check:{check:true}, radio_sel:{value:comp_name[0], noupdate:true}, comp_name:comp_name, comp_filt:comp_filt, tra_filt:tra_filt, tra_name:tra_name, single:{}, radio:{value:"select", noupdate:true}};
	
	rpf2.filter.push(filter);

	select_bubble_but_prop("RightMenu","PopFilt",["num"],rpf2.filter.length-1);
	inter.bubble.popfilt = {rpf:rpf, filter:rpf2.filter[rpf2.filter.length-1]};		
}


/// Makes sure that fraction can only be specified for one filter
function update_popfilt_fraction(p,cl)
{	
	let rpf = get_inf_res().plot_filter;
	let filter = rpf.species[p].filter;
	for(let i = 0; i < filter.length; i++){
		let fi = filter[i];
		if(fi.cl != cl){
			if(fi.fraction.check == true){
				fi.fraction.check = false;
			}
		}	
	}
}


/// Activates when editing filter is done
function done_filter()
{
	if(inter.bubble.check_radio_press == true) inter.graph.init = undefined;
	close_bubble();	
}
	

/// Converts the knot text entries into numerical value
function get_times(knot,details)
{
	let t_start = details.t_start;
	let t_end = details.t_end;
	if(details.inf_t_end != undefined) t_end = details.inf_t_end;
	
	if(details.ppc_t_start != undefined) t_start = details.ppc_t_start;
	if(details.ppc_t_end != undefined) t_end = details.ppc_t_end;

	let times = [];
	for(let k = 0; k < knot.length; k++){
		let t;
		let val = knot[k];
		if(val == "start"){		
			if((typeof t_start) == "string" && (t_start == "" || isNaN(t_start))){
				return err("The simulation start time must be specified");
			}
			val = Number(t_start);
		}
		else{
			if(val == "end"){
				if((typeof t_end) == "string" && (t_end == "" || isNaN(t_end))){
					return err("The simulation end time must be specified");
				}
				val = Number(t_end);
			}
			else{
				if(isNaN(val)){
					return err("Knot value not a number");
				}
				val = Number(val);
			}
		}
		times.push(val);
	}
	return { err:false, times:times};
}


/// Given a set of values knots and timelines 
function set_spline_value(value,times,tp)
{
	let i = 0; // Indexes times
	let div = [];
	for(let ti = 0; ti < tp.length-1; ti++){
		let tmid = (tp[ti]+tp[ti+1])/2;
		while(i+1 < times.length && times[i+1] < tmid) i++;
		if(i+1 >= times.length) error("Range problem");
		let f = (times[i+1]-tmid)/(times[i+1]-times[i]);
		div.push({t:tmid, value:f*value[i]+(1-f)*value[i+1]});
	}
	return div;
}


/// The scrollable box which allows the used to specify population filter
function filterpos_scrollable(lay)
{
	let cx = 0.5, cy = 0;
	let dx = lay.dx;
	let op = lay.op;
	
	let popfilt = inter.bubble.popfilt;
	
	let filt = popfilt.filter;
	let rpf = popfilt.rpf;
	let claa = rpf.species[filt.p].cla[filt.cl];
	
	if(rpf.species[filt.p].sel_class.cl == filt.cl){  // Filtering on observered cl
		switch(subsubtab_name()){
		case "Populations": case "Individuals":
			cy = add_checkbox_list(lay,claa.comp,"name",filt.comp_filt,cy,dx); 
			break;

		case "Transitions": 
			cy = add_checkbox_list(lay,claa.tra,"name",filt.tra_filt,cy,dx); 
			break;
		default: error("tree op not recognised"); break;
		}
	}
	else{
		if(filt.radio.value == "select"){
			cy = add_checkbox_list(lay,claa.comp,"name",filt.comp_filt,cy,dx);
		}
		else{
			cy = add_radio_list(lay,claa.comp,"name",filt.radio_sel,cy,dx);
		}
	}
	
	cy = lay.add_gap(cy,0.1);
	
	return cy;
}

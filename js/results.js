// Functions associated with loading up results from BICI

"use strict";

/// Adds information about the model into the results (such that results can be plotted)
function results_add_model(result,details,siminf)
{
	result.siminf = siminf;

	result.details = copy(details);
	
	result.species = copy(model.species);
	
	global_comp_trans_init(result.species);
	
	result.param = copy(model.param);
	
	create_param_const(result,siminf);            // Creates constant values for non-variable parameters 
	
	create_spline();

	result.spline = copy(model.spline);
}


/// Creates constant values for non-variable parameters 
function create_param_const(result,siminf)
{
	if(siminf == "sim"){  // For simulation removes parameter associated with observation process
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
			if(!(par.variety == "dist" && par.sim_sample.check == true) && par.type != "derive_param"){
				par.kind = "const";
			}
		}
	}		
}


/// Sends an error that the sample could not be loaded
function alert_sample(warn,code)
{
	alert_import(warn+" CODE: "+code);
}


/// Sends an error that the sample could not be loaded
function alert_param_sample(warn,code)
{
	alert_import(warn+" CODE: "+code);
}

function read_param_samples_file(chain,file,result)
{
	let fs = require('fs');	
	
	let file_full = inter.imp.outputdir+"/"+file;
	
	if(check_file_exists(file_full) != "success") return;
	
	let te = fs.readFileSync(file_full, 'utf8');
	
	let warn = "Problem loading file '"+file+"'";
	
	read_param_samples(chain,te,result,warn)
}

/// Reads in a parameter samples from a file
function read_param_samples(chain,te,result,warn)
{
	let lines = te.split('\n');
	
	let spl = comma_split(lines[0]);
	
	if(spl[0] != "State"){ alert_param_sample(warn,-1); return;}
	
	let J = spl.length;
	
	let ref = [];
	ref.push({});
	for(let i = 1; i < J; i++){
		let name = remove_quote(spl[i]).replace(/->/g,"→");
	
		if(find_in(like_name,name) != undefined){
			ref.push({th:result.param.length, index:[]});
			result.param.push({name:name, full_name:name, dep:[], time_dep:false, age_dep:false, kind:"variable", variety: "likelihood", like:true});
		}
		else{
			let pp = get_param_prop(name);
			let th = find(result.param,"name",pp.name);
			if(th == undefined){ alert_param_sample(warn,-2); return;}
			
			let par = result.param[th];
			if(pp.dep.length != par.dep.length){ alert_param_sample(warn,-3); return;}
			
			let ind = [];
			for(let j = 0; j < par.dep.length; j++){
				ind[j] = find_in(par.list[j],pp.dep[j]);
				if(ind[j] == undefined){ alert_param_sample(warn,4); return;}
			}
			
			ref.push({th:th, index:ind});
		}
	}

	let temp = [];
	
	for(let th = 0; th < result.param.length; th++){
		let par = result.param[th];
		
		if(par.kind == "variable"){
			let value; 
			if(par.dep.length == 0){
				value = "unset";
			}
			else{
				let list = par_find_list(par);
				value = par_find_template(list);
			}
			temp[th] = value;
		}
	}
	
	for(let li = 1; li < lines.length; li++){
		if(lines[li] != ""){
			let spl = lines[li].split(",");
		
			if(spl.length != J){ alert_param_sample(warn,-4); return;}
			let param = copy(temp);
			for(let i = 1; i < J; i++){
				let num = Number(spl[i]); if(isNaN(num)){ alert_param_sample(warn,-5); return;}
				
				let th = ref[i].th;
				if(ref[i].index.length == 0) param[th] = num;
				else set_element(param[th],ref[i].index,num);
			}
		
			result.par_sample.push({ chain:chain, num:Number(spl[0]), param:param});
		}
	}
		
	if(result.on != true) result.on = true;
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


/// Returns the name of the algorithm
function alg(result)
{
	return result.details.algorithm.value;
}


/// Calculate the burn-in time
function calculate_burnin(result)
{
	if(alg(result) != "DA-MCMC"){ result.burnin = 0; return;}
	
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

	result.burnin = Math.floor(result.burnin_factor*max);
}


/// Reads in state samples from a file
function read_state_samples_file(chain,file,result)
{
	let fs = require('fs');	
	
	let file_full = inter.imp.outputdir+"/"+file;
	
	if(check_file_exists(file_full) != "success") return;
	
	let te = fs.readFileSync(file_full, 'utf8');

	let warn = "Problem loading file '"+file+"'";
	
	read_state_samples(chain,te,result,warn);
}
	

/// Reads in state samples from text
function read_state_samples(chain,te,result,warn)
{
	let sec = te.split("<<");
	
	for(let s = 1; s < sec.length; s++){
		read_state_sample("<<"+sec[s],chain,result,warn);
	}		
}


/// Reads a state file from text
function read_state_sample(te,chain,result,warning)
{
	let lines = te.split('\n');
	for(let i = 0; i < lines.length; i++) lines[i] = lines[i].trim();
				
	let sample = { chain:chain, timepoint:[], param:[], spline:[], species:[]};
	
	for(let th = 0; th < result.param.length; th++) sample.param[th] = undefined;
	
	let mode = "none";
	
	let sample_num;
	
	let p;
	for(let i = 0; i < lines.length; i++){
		let li = lines[i];
		if(li != ""){
			let warn = warning+": Line "+(i+1)+"   "+li;
	
			if(li.substr(0,1) == "<"){
				if(li.substr(li.length-1,1) != ">"){ alert_sample(warn,1); return;}	
				li = li.substr(1,li.length-2);
			
				let spl = li.split(" ");
			
				switch(spl[0]){
				case "<STATE":
					{
						let spl2 = spl[1].split(">");
						sample.num = spl2[0];
					}
					break;
					
				case "PARAMETERS":
					mode = "param"; 
					break;
					
				case "SPECIES":
					{
						mode = "species";
						if(spl.length != 2){ alert_sample(warn,2); return;}	
						let sp_name = remove_quote(spl[1]);
						
						p = find(result.species,"name",sp_name);
						if(p == undefined){ alert_sample(warn,3); return;}	
						
						switch(result.species[p].type){
						case "Population": sample.species[p] = {type:"Population", cpop_init:[], timeline:[]}; break;
						case "Individual": sample.species[p] = {type:"Individual", individual:[], cpop_init:[], timeline:[]}; break;
						default: error("PROBLEM"); break;
						}
						
						for(let ti = 0; ti < sample.timepoint.length; ti++){
							sample.species[p].timeline[ti] = {};
						}
					}
					break;
				
				case "INITIAL":
					if(p == undefined){ alert_sample(warn,4); return;}	
					if(sample.species[p].type == "Individual"){ alert_sample(warn,5); return;}	
					mode = "initial";
					break;
					
				case "TRANSITIONS": 
					if(p == undefined){ alert_sample(warn,6); return;}	
					if(sample.species[p].type == "Individual"){ alert_sample(warn,7); return;}	
					mode = "trans"; 
					break;
				
				case "TIMEPOINT":
					mode = "timepoint";
					break;
					
				case "INDIVIDUALS":
					if(p == undefined){ alert_sample(warn,10); return;}	
					if(sample.species[p].type == "Population"){ alert_sample(warn,11); return;}	
					mode = "individual";
					break;
					
				case "DERIVED":
					mode = "derive";
					break;
					
				default: alert_sample(warn,12); return;
				}
			}
			else{
				switch(mode){
				case "param": case "derive":
					i = get_param_value(i,sample,lines,result,warn,mode);
					break;
					
				case "initial":
					{
						let spl = comma_split(lines[i]);
						if(spl.length != 2){ alert_sample(warn,21); return;}
						
						let name = spl[0], pop = spl[1];
						if(!(name == "compartment" && pop == "population")){
							let sp = result.species[p];
							
							let c = find(sp.comp_gl,"name",name);
							if(c == undefined){ alert_sample(warn,22); return;}
							
							if(isNaN(pop)){ alert_sample(warn,23); return;}
							
							sample.species[p].cpop_init[c] = Number(pop);
						}
					}
					break;
					
				case "trans":
					{
						let sp = result.species[p];
						
						let ssp = sample.species[p];
						
						let spl = lines[i].split(",");
						if(spl.length != sp.tra_gl.length){ alert_sample(warn,25); return;}
						i++;
						
						for(let j = 0; j < sp.tra_gl.length; j++){
							if(spl[j].replace(/->/g,"→") != sp.tra_gl[j].name){ alert_sample(warn,24); return;}
						}
						
						for(let ti = 0; ti < sample.timepoint.length-1; ti++){
							let spl = lines[i].split(",");
							if(spl.length != sp.tra_gl.length){ alert_sample(warn,25); return;}
							i++;
						
							let transnum = [];
						
							for(let j = 0; j < sp.tra_gl.length; j++){
								let num = Number(spl[j]);
								if(isNaN(num)){ alert_sample(warn,27); return;}
								transnum[j] = num;
							}
							
							ssp.timeline[ti].transnum = transnum;
						}
					}
					break;
				
				case "timepoint":
					{
						let spl = lines[i].split(",");
						
						for(let i = 0; i < spl.length; i++){
							let t = Number(spl[i]);
							if(isNaN(t)){ alert_sample(warn,29); return;}
							sample.timepoint.push(t);
						}
					}
					break;
					
				case "individual":
					{
						let sp = result.species[p];
						let ssp = sample.species[p];
				
						let spl = lines[i].split(",");
				
						if(spl.length > 2 && spl[0] == "name" && spl[1] == "init"){
							sp.ind_effect = [];
							for(let k = 2; k < spl.length-1; k++) sp.ind_effect.push(spl[k]);	
						}
						else{	
							if(spl.length != 3+sp.ind_effect.length){ alert_sample(warn,30); return;}
							
							let name = spl[0].trim();
							
							let c;
								
							let cinit_str = spl[1].trim();
							if(cinit_str == "") c = "out";
							else{
								if(cinit_str == "+") c = "Source";
								else{								
									c = find(sp.comp_gl,"name",spl[1]);
									if(c == undefined){ alert_sample(warn,31); return;}
								}
							}
							
							let cinit = c;
								
							let ie = [];
							for(let k = 0; k < sp.ind_effect.length; k++){
								let val = Number(spl[k+2]);
								if(isNaN(val)){ alert_sample(warn,32); return;}
								ie.push(val);
							}
							
							let ev=[];
							
							let ev_text = spl[spl.length-1].trim();
							if(ev_text != ""){		
								let col = ev_text.split(" ");		
								
								for(let e = 0; e < col.length; e++){
									let tspl = col[e].split(":");
									if(tspl.length != 2){ alert_sample(warn,33); return;}
									let t = Number(tspl[1]); if(isNaN(t)){ alert_sample(warn,34); return;}
									
									if(c == "out"){   // Individual enters the system
										c = find(sp.comp_gl,"name",tspl[0]);
										if(c == undefined){ alert_sample(warn,35); return;}
										
										ev.push({type:"enter", t:t, c:c});
									}
									else{
										let tr = tspl[0];
										
										if(c == "Source"){   // Individual enters from source
											tr = tr.replace(/->/g,"→");
										
											let trg = find(sp.tra_gl,"name",tr);
											if(trg == undefined){ alert_sample(warn,36); return;}
											
											let trgl = sp.tra_gl[trg];
											if(trgl.i != c){ alert_sample(warn,37); return;}
											ev.push({type:"trans", t:t, trg:trg});
											c = trgl.f;
										}
										else{
											if(tr == "-"){ // Individual leaves the system
												ev.push({type:"leave", t:t});
											}
											else{			
												if(tr.length > 6 && tr.substr(0,6) == "move->"){
													let nam = tr.substr(6);
													let co = find_comp_from_name(nam,sp.name);
													if(co.warn != undefined){ alert_sample(warn,38); return;}
													
													let cnew = update_c_comp(sp,c,co.cl,co.c);
													if(c != cnew){
														ev.push({type:"move", t:t, cl:co.cl, ci:c, cf:cnew});
														c = cnew;
													}
												}
												else{
													let k = find(sp.comp_gl[c].tra_gl,"name",tr);
													if(k == undefined){ alert_sample(warn); return;}
													let trg = sp.comp_gl[c].tra_gl[k].num;
													let trgl = sp.tra_gl[trg];
													if(trgl.i != c){ alert_sample(warn,39); return;}
													ev.push({type:"trans", t:t, trg:trg});
													c = trgl.f;
												}
											}
										}
									}
								}
							}
							
							let ind = {name:name, cinit:cinit, ev:ev, ie:ie};
							
							ssp.individual.push(ind);
						}
					}
					break;
					
				default: alert_sample(warn,40); return;
				}
			}
		}
	}

	generate_spline(result,sample);
	
	generate_cpop_init_from_ind(result,sample);

	generate_transnum_from_ind(result,sample);

	generate_cpop_from_transnum(result,sample);

	for(let th = 0; th < result.param.length; th++){
		let par = result.param[th];
		if(par.kind == "variable" && sample.param[th] == undefined){
			if(!(par.type == "derive_param" && par.time_dep == true) && !(par.variety == "likelihood")){
				alert_sample(warning+": Parameter '"+result.param[th].name+"' not set",41); return;
			}
		}
	}		

	result.sample.push(sample);
	if(result.on != true) result.on = true;
}


/// Gets parameter value from line in the file
function get_param_value(i,source,lines,result,warn,mode)
{
	let param = result.param;

	let spl = comma_split(lines[i]);
	let name = remove_quote(spl[0].replace(/->/g,"→"));
	
	let th = find(result.param,"name",name);
	if(th == undefined){pr(name+" name"); alert_sample(warn,13); return;}	
	
	let par = result.param[th];
	
	let time_var = false;
	if(par.time_dep == true && mode == "derive"){
		time_var = true;
		par.dep.pop();
	}
	
	let value;
	switch(spl.length){
	case 1:           // Set 
		let list = par_find_list(par);
		value = par_find_template(list);
		let comb_list = generate_comb_list(list);
		par.list = list;
		par.comb_list = comb_list;
	
		i++;
		let spl_head = comma_split(lines[i]);
		
		let ndep = par.dep.length;
		
		if(spl_head.length != ndep+1){ alert_sample(warn,14); return;}	
		
		for(let d = 0; d < ndep; d++){
			if(par.dep[d] != spl_head[d]){ alert_sample(warn,15); return;}	
		}

		if(spl_head[ndep] != "Value"){ alert_sample(warn,16); return;}	
			
		i++;
		for(let j = 0; j < comb_list.length; j++){
			let spl_row = lines[i].split(",");
			if(spl_row.length != ndep+1){ alert_sample(warn,17); return;}	
			let index=[];
			for(let d = 0; d < ndep; d++){
				let k = find_in(list[d], spl_row[d]);
				if(k == undefined){ alert_sample(warn,18); return;}	
				index.push(k);
			}
			let val = spl_row[ndep];
			
			if(!time_var){
				val = Number(val);
				if(isNaN(val)){ alert_sample(warn,100); return;}
			}

			set_element(value,index,val);		
			i++;
		}
		i--;
		
		for(let j = 0; j < comb_list.length; j++){
			if(get_element(value,comb_list[j].index) == undefined){ alert_sample(warn,101); return;}	
		}
		break;
		
	case 2:    // Directly set parameter value (univariate)
		value = spl[1];
		
		if(!time_var){
			value = Number(value);
			
			if(isNaN(value)){ alert_sample(warn,19); return;}	
		}
		break;
		
	default: alert_sample(warn,20); return;
	}
	
	source.param[th] = value;
	
	return i;
}
	

/// Converts a string in | format into spline data
function get_spline_from_string(st,sample)
{
	let spl = st.split("|");
									
	let tim = sample.timepoint;
	if(spl.length != tim.length-1){
		alert_sample(warn,180); return;
	}	
	let val_new = [];
	
	for(let ti = 0; ti < tim.length-1; ti++){
		let value = Number(spl[ti]);
		if(isNaN(value)){ alert_sample(warn,181); return;}
		val_new.push({t:(tim[ti]+tim[ti+1])/2, value:value});
	}
	
	return val_new;
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
				if(c != "out"){
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
			for(let ti = 0; ti < sample.timepoint.length-1; ti++){
				let tim = ssp.timeline[ti];
				tim.transnum = [];
				for(let k = 0; k < sp.tra_gl.length; k++) tim.transnum[k] = 0;
				tim.dpop = [];
				for(let c = 0; c < sp.comp_gl.length; c++) tim.dpop[c] = 0;
			}
			
			let tim = sample.timepoint;
			let dt = 1;
			if(tim.length > 1) dt = tim[1]-tim[0];
			let t_start = tim[0];
	
			for(let i = 0; i < ssp.individual.length; i++){
				let ind = ssp.individual[i];
				
				let c = ind.cinit;
				
				for(let e = 0; e < ind.ev.length; e++){
					let ev = ind.ev[e];
					let t = ev.t;
				
					let ti = Math.floor((t-t_start)/dt);
					if(ti < 0 || ti >= ssp.timeline.length-1)	error("Out of range");
					
					switch(ev.type){
					case "trans":
						let tra = sp.tra_gl[ev.trg];
						if(tra.i != c) error("c is not consistent");
						c = tra.f;
						ssp.timeline[ti].transnum[ev.trg]++;
						break;
						
					case "enter":
						ssp.timeline[ti].dpop[ev.c]++;
						if(c != "out") error("Should be out");
						c = ev.c;
						break;
						
					case "leave":
						if(e != ind.ev.length-1) error("Leave should be last event");
						ssp.timeline[ti].dpop[c]--;
						break;
						
					case "move":
						ssp.timeline[ti].dpop[ev.ci]--;
						ssp.timeline[ti].dpop[ev.cf]++;
						break;
						
					default: error("event type not recognised:"+ev.type); break;
					}
				}
			}
		}
	}
}


/// Generates cpop using global transitions specified in transnum
function generate_cpop_from_transnum(result,sample)
{
	for(let p = 0; p < sample.species.length; p++){
		let sp = result.species[p];
		let ssp = sample.species[p];

		for(let t = 0; t < ssp.timeline.length; t++){
			let tline = ssp.timeline[t];
			if(t == 0){
				tline.cpop = ssp.cpop_init;
			}
			else{
				let last = ssp.timeline[t-1];
				tline.cpop = [];
				for(let c = 0; c < sp.comp_gl.length; c++){
					tline.cpop[c] = last.cpop[c];
					if(sp.type == "Individual") tline.cpop[c] += last.dpop[c];
				}
				
				for(let j = 0; j < sp.tra_gl.length; j++){
					let tra = sp.tra_gl[j];
					let num = last.transnum[j];
					let i = tra.i, f = tra.f;
					if(i != "Source") tline.cpop[i] -= num;
					if(f != "Sink") tline.cpop[f] += num;
				}
			}
		}					
	}
}


/// Intialises global compartments and transitions
function global_comp_trans_init(species)
{
	for(let p = 0; p < species.length; p++){
		let sp = species[p];
		
		sp.comp_mult=[];
		
		let N = 1;                                                  // Works out global conversion
		for(let cl = 0; cl < sp.ncla; cl++){
			sp.comp_mult.push(N);
			N *= sp.cla[cl].comp.length;
		}
		
		sp.N = N;
	
		sp.comp_gl = [];                                                 // Constructs global compartments
		
		for(let c = 0; c < N; c++){                
			let co = {cla_comp:[], tra_gl:[]};
			
			let na = "";
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				let cc = Math.floor(c/sp.comp_mult[cl])%claa.comp.length;
				co.cla_comp.push(cc);

				if(na != "") na += "|";
				na += claa.comp[cc].name;
			}
			co.name = na;
			
			sp.comp_gl.push(co);
		}

		sp.tra_gl = [];
		for(let cl = 0; cl < sp.ncla; cl++){                              // Constructs global transitions
			let cla = sp.cla[cl];
			for(let j = 0; j < cla.ntra; j++){
				let tr = cla.tra[j];

				for(let c = 0; c < N; c++){
					let co = sp.comp_gl[c];		

					let cc = co.cla_comp[cl];
					let ii, ff;

					if(tr.i == "Source"){
						if(tr.f == cc){ ii = "Source"; ff = c;}
					}
					else{
						if(tr.f == "Sink"){
							if(tr.i == cc){ ii = c; ff = "Sink";}
						}
						else{
							if(tr.i == cc){
								let dif = Number(tr.f) - Number(tr.i);
								ii = c; ff = c + dif*sp.comp_mult[cl];
							}
						}
					}

					if(ii != undefined){
						if(ii != "Source" && ii != "Sink"){
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
	case "Graph": case "Graph (lines)": case "Graph (CI)": case "Data": return "Line plot";
	case "Scatter": return "Scatter plot";
	case "Compartment": return "Comp plot";
	case "Individual": return "Individual plot";
	case "Histogram": return "Histogram plot";
	case "Matrix": return "Matrix plot";
	case "Stat Table": return "Stat table plot";
	case "No Graph": return "No graph plot";
	default: error("Error variety:"+view); break;
	}
}


/// Initialises filters used to plot the results
function intialise_plot_filters(result)
{
	result.plot_filter = {};
	let rpf = result.plot_filter;
	
	// Initialises distribution plots
	let pos_bin = [];
	pos_bin.push({te:10}); pos_bin.push({te:20}); pos_bin.push({te:50}); 
	pos_bin.push({te:100}); pos_bin.push({te:200}); pos_bin.push({te:500});
	
	let pos_h = [];
	pos_h.push({te:0.01}); pos_h.push({te:0.02}); pos_h.push({te:0.05}); 
	pos_h.push({te:0.1}); pos_h.push({te:0.2}); pos_h.push({te:0.5});
	
	rpf.dist_settings = { show_mean:{check:true}, show_prior:{check:true}, radio:{value:"kde"}, sel_h:{te:0.05}, pos_h:pos_h, sel_bin:{te:50}, pos_bin:pos_bin};
	
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

	/// Possibilities for chain distribution 
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
	rpf.sel_paramview = copy(rpf.pos_paramview[0]);
	
	// Possibilies for individual view
	rpf.pos_indview = [];
	rpf.pos_indview.push({te:"Timeline"});
	rpf.pos_indview.push({te:"Table"});
	rpf.sel_indview = copy(rpf.pos_indview[0]);
	
	// Possiblities for time-step
	rpf.pos_timestep = get_possible_timestep(result.details);
	rpf.sel_timestep = copy(rpf.pos_timestep[0]);
	
	// Possibilities for species
	let pos_species = [];
	for(let p = 0; p < result.species.length; p++){
		pos_species.push({te:result.species[p].name, p:p});
	}
	rpf.sel_species = copy(pos_species[0]);
	rpf.pos_species = pos_species;
	
	// Possibilities for classification
	rpf.species = [];
	for(let p = 0; p < result.species.length; p++){
		let sp = result.species[p];
		
		let pos_class = [];
		for(let cl = 0; cl < sp.ncla; cl++){
			pos_class.push({te:sp.cla[cl].name, cl:cl});
		}
	
		rpf.species[p] = {};
		rpf.species[p].sel_class = copy(pos_class[0]);
		rpf.species[p].pos_class = pos_class;
		
		rpf.species[p].filter=[];
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
	
	rpf.sel_sample = copy(pos_sample[0]);
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
	if(result.siminf == "inf"){
		for(let p = 0; p < result.species.length; p++){
			let sp = result.species[p];
			for(let i = 0; i < sp.inf_source.length; i++){
				let so = sp.inf_source[i];
				let spl = so.desc.split("#");
				let desc = spl[0].trim();
				switch(so.type){
				case "Population": pos_popdata.push({te:desc, p:p, i:i}); break;
				case "Pop. Trans.": pos_poptransdata.push({te:desc, p:p, i:i}); break;
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
	pos_view.push({te:"Compartment"});
	if(pos_popdata.length) pos_view.push({te:"Data"});
	
	rpf.sel_view = copy(pos_view[0]);
	rpf.pos_view = pos_view;
	
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
	
	// Possibilities for splines
	
	let pos_spline = [];
	for(let spl = 0; spl < result.spline.length; spl++){
		pos_spline.push({te:result.spline[spl].name, i:spl});
	}
	if(result.spline.length > 0) rpf.sel_spline = copy(pos_spline[0]);
	rpf.pos_spline = pos_spline;
	
	// Possibilities for viewing splines
	let pos_spline_view = [];
	if(pos_sample.length == 1){
		pos_spline_view.push({te:"Graph"});
	}
	else{
		pos_spline_view.push({te:"Graph (CI)"});
		pos_spline_view.push({te:"Graph (lines)"});	
	}
	
	rpf.sel_spline_view = copy(pos_spline_view[0]);
	rpf.pos_spline_view = pos_spline_view;
	
	result.burnin_factor = 0.2;
	calculate_burnin(result);
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
	

/// Gets a list of all possible parameter views
function get_pos_paramview(result)
{
	let total_param_list = [];
	
	let pos_paramview = [];
	let list = [];
	
	for(let th = 0; th < result.param.length; th++){
		let par = result.param[th];
		if(par.dep.length == 0 && par.like != true){
			list.push({th:th, index:[]});
			if(par.kind != "const"){
				total_param_list.push({th:th, index:[], name:param_name_index(par,[])});
			}
		}
	}
	
	if(list.length > 0){
		pos_paramview.push({te:"Univariate", param:true, radio:{value:0, param:true},  list:list});	
	}

	for(let th = 0; th < result.param.length; th++){
		let par = result.param[th];
		if(par.dep.length > 0){
			let list = [];
			for(let j = 0; j < par.comb_list.length; j++){	
				let ind = par.comb_list[j].index;
				
				list.push({th:th, index:ind});
				
				if(par.kind != "const"){
					total_param_list.push({th:th, index:ind, name:param_name_index(par,ind)});
				}
			}
			pos_paramview.push({te:remove_eq_quote(par.full_name), param:true, radio:{value:0, param:true}, list:list});
		}
	}

	result.total_param_list = total_param_list;

	if(total_param_list.length > 1){
		pos_paramview.push({te:"Scatter", param:true, radio:{value:0, param:true}, radioy:{value:1, param:true}, list:total_param_list});
		
		pos_paramview.push({te:"Correlation", param:true, radio:{value:0, param:true}, radioy:{value:1, param:true}, list:total_param_list});
	}
	
	let list_like = [];
	for(let th = 0; th < result.param.length; th++){
		let par = result.param[th];
		if(par.like == true){
			list_like.push({th:th, index:[]});
		}
	}
	
	if(list_like.length > 0){
		pos_paramview.push({te:"Likelihoods", param:true, radio:{value:0, param:true},  list:list_like});	
	}
	
	return pos_paramview;
}


/// Gets possible timesteps based on dt
function get_possible_timestep(details)
{
	let dt = details.timestep;
	
	let range = Number(details.t_end) - Number(details.t_start);
	let dt_max = range/4;

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
	if(j*dt < dt_max) pos_timestep.push({te:dt*j, i:j});
}

				
/// Adds bubble content which allows for a filter to be added
function filter_bubble(bu,cont)
{
	let filt=[];
	
	let result = bu.result;
	let rpf = result.plot_filter;
	let p = rpf.sel_species.p;
	
	let rpf2 = rpf.species[p];
	let cl = rpf2.sel_class.cl;
	filt[cl] = "not";
	
	for(let i = 0; i < rpf2.filter.length; i++) filt[rpf2.filter[i].cl] = "not";
	
	let sp = result.species[p];
	
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
	
		bubble_addtext(cont,"Classification:");
		
		cont.y -= 0.4;
		bubble_adddropdown(cont,0,7,inter.bubble.sel_filter,pos);
		cont.y += 0.3;
		
		let ac; if(inter.bubble.sel_filter.te != select_str) ac = "AddFilter2";
			
		add_end_button(cont,"OK",ac,{result:result});
	}
}


/// Adds a filter to a population plot
function add_pop_filter(result)
{
	let rpf = result.plot_filter;
	let p = rpf.sel_species.p;

	let sp = result.species[p];
	
	let cl = inter.bubble.sel_filter.cl;
	
	let rpf2 = rpf.species[p];
	
	// Possibilities for filter
	let pos_comp = [];
	for(let c = 0; c < sp.cla[cl].comp.length; c++){
		pos_comp.push({te:sp.cla[cl].comp[c].name, c:c});
	}
	pos_comp.push({te:"All", c:"All"});
	pos_comp.push({te:"Remove", ac:"Remove filter", p:p, filter:rpf2.filter, num:rpf2.filter.length, c:"Remove"});
	
	let filter = {name:sp.cla[cl].name, cl:cl, sel_comp:{te:select_str}, pos_comp:pos_comp};
	
	rpf2.filter.push(filter);

	close_bubble();
}


/// Generates values for the spline based on the parameter samples
function generate_spline(result,sample)
{
	let tp = sample.timepoint;
	
	sample.spline = [];

	for(let i = 0; i < result.spline.length; i++){
		let th = result.spline[i].th;
		let par = result.param[th];
		
		if(par.type == "derive_param"){
			let value;
			if(par.dep.length == 0){
				value = get_spline_from_string(sample.param[th].value,sample);	
			}
			else{
				let list = par_find_list(par);
			
				value = par_find_template(list);
				let comb_list = generate_comb_list(list);
				
				for(let i = 0; i < comb_list.length; i++){
					let index = comb_list[i].index;
					let str = get_element(sample.param[th].value,index);
				
					let val_new = get_spline_from_string(str,sample);
					if(list.length == 0) value = val_new;
					else set_element(value,index,val_new);
				}
			}
			sample.spline[i] = value;
		}
		else{
			let ps = sample.param[th];
			
			let list = par_find_list(par);
		
			let knot = par.spline.knot;
			let times = get_times(knot,tp[0],tp[tp.length-1]);
			if(typeof times == 'string') error(times);
			
			if(list[list.length-1].length != knot.length) error("Length does not match");
			
			let list_new = copy(list);
			list_new.pop();
			
			let value;
			if(par.kind == "const") value = copy(result.param[th].value);
			else value = ps;
		
			let comb_list = generate_comb_list(list_new);
		
			for(let j = 0; j < comb_list.length; j++){
				let index = comb_list[j].index;
				let val_new = set_spline_value(get_element(value,index),times,tp);
				if(list_new.length == 0) value = val_new;
				else set_element(value,index,val_new);
			}

			sample.spline[i] = value;
		}
	}
}


/// Converts the knot text entries into numerical value
function get_times(knot,t_start,t_end)
{
	let times = [];
	for(let k = 0; k < knot.length; k++){
		let t;
		let val = knot[k];
		if(val == "start"){		
			if((typeof t_start) == "string" && (t_start == "" || isNaN(t_start))) return "The simulation start time must be specified";
			val = Number(t_start);
		}
		else{
			if(val == "end"){
				if((typeof t_end) == "string" && (t_end == "" || isNaN(t_end))) return "The simulation end time must be specified";
				val = Number(t_end);
			}
			else{
				if(isNaN(val)) return "Knot value not a number";
				val = Number(val);
			}
		}
		times.push(val);
	}
	return times;
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

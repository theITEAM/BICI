"use strict";
// Functions which check that the data being added into the model is valid

/// imports an equation value
function import_eqn_value(eqn,te,op)
{
	te = char_replace(te);
 
	eqn.te = te;
	extract_equation_properties(eqn);
	if(eqn.warn.length > 0){
		alert_import("Equation '"+te+"' – "+ eqn.warn[0].te);
	}
}


/// Check for initial pop
function check_initial_pop_error(end)
{ 
	for(let p = 0; p < model.species.length; p++){
		for(let loop = 0; loop < 2; loop++){
			let sp = model.species[p];
		
			let type = get_siminf(loop);
			let source = get_source(type,p); 
		
			let ninitpop = 0;
			for(let i = 0; i < source.length; i++){
				let so = source[i];
				
				if(so.type ==  "Init. Pop.") ninitpop++;
			}
			
			switch(type){
			case "sim":
				if(ninitpop > 1){
					return "Only one initial population should be set for species '"+sp.name+"'";
				}
				break;
			
			case "inf":
				if(ninitpop > 1){
					return "Only one initial population should be set for species '"+sp.name+"'";
				}
				break;
			}
		}
	}
}


/// Checks that parameter specification are complete
function check_param_complete()
{
	let param = copy(strip_heavy(model.param));

	update_model();

	let flag=[];
	for(let th = 0; th < param.length; th++){
		let par = param[th];
		
		let th2 = find(model.param,"name",par.name);
		if(th2 == undefined){
			alert_line("Parameter "+par.full_name+" is not used in the model",par.import_line);
			return;
		}
		else{
			let par2 = model.param[th2];
			if(equal_vec(par.dep,par2.dep) == false){
				alert_line("Parameter "+par.full_name+" has a different dependency in the model "+par2.full_name,par.import_line);
			}
		
			// Checks the prior are within those possible
			if(par2.prior.type.te != select_str){
				if(option_error("prior",par2.prior.type.te,par2.pri_pos,par.import_line) == true) return;
			}
		
			flag[th2] = true;
		}
	}
	
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let i = 0; i < sp.fix_eff.length; i++){
			let fe = sp.fix_eff[i];
			if(fe.defined != true){
				alert_noline("Fixed effect 〈"+fe.name+"〉 is unspecified. It must be defined through the 'fixed-effect' command.");
			}
		}
	}
	
	for(let th = 0; th < model.param.length; th++){
		let par = model.param[th];
		if(flag[th] != true && !is_in_obsmodel(par) && par.derive != true){
			if(!par.dist_mat && !par.iden_mat && !par.den_vec){
				alert_noline("Parameter "+par.full_name+" is not specified as a 'param'");
			}
		}
	}
}


/// Determines if text string is for a file
function is_file(te)
{
	if(te == "") return false;
	
	if(te.te != undefined || te.ref != undefined) return true;
	
	let k = te.length-1;
	let kmin = k-10;
	while(k > 0 && k > kmin && te.substr(k,1) != ".") k--;
	if(k == 0 || k == kmin) return false;
	let end = te.substr(k).toLowerCase();
	if(end != ".csv" && end != ".txt" && end != ".tsv" && end != ".geojson") return false;

	return true;
}


/// Calls when cannot find a compartment in 
function cannot_find_trans_comp(te,p,cl,co)
{
	let st = "In transition '"+te+"' cannot find compartment '"+co+"'";

	if(co == ""){
		st = "In transition '"+te+"' the initial and final compartment aren't specified";
	}
	else{
		let sp = model.species[p];
		
		if(sp.ncla > 1) st += " in classification '"+sp.cla[cl].name+"'";
		st += ".";
		
		for(let cl2 = 0; cl2 < sp.ncla; cl2++){
			if(cl2 != cl){
				let c = hash_find(sp.cla[cl2].hash_comp,co);
				if(c != undefined){
					st += " Perhaps this transition should be placed in classification '"+sp.cla[cl2].name+"'?";
				}
			}
		}
	}
		
	alert_import(st);
}
			

/// Prints the number of columns and rows in a table
function print_row_col(tab)
{
	let te = "(there ";
	if(tab.ncol == 1) te += "is "+tab.ncol+" column"; 
	else te += "are "+tab.ncol+" columns";
	te += " and ";
	if(tab.nrow == 1) te += tab.nrow+" row"; 
	else te += tab.nrow+" rows)";
	
	return te;
}


/// Looks at how memory is being used
function profiling()
{
	prr(mem(model)+" model");
	prr(mem(model_store)+" model_store");
	prr(mem(sim_result)+" sim_result");
	prr(mem(inf_result)+" inf_result");
	prr(mem(ppc_result)+" ppc_result");
	prr(mem(sim_result_import)+" sim_result_import");
	prr(mem(inf_result_import)+" inf_result_import");
	prr(mem(ppc_result_import)+" ppc_result_import");
	//prr(mem(import_bscript)+" import_bscript");
	prr(mem(edit_source)+" edit_source");
	prr(mem(input)+" input");
	
	
	mem_split(inf_result,"inf_result");
	
	mem_split(inf_result.bscript,"inf_result->bscript");
	
	mem_split(inf_result.sample[0],"inf_result->sample[0]");
	
	mem_split(inf_result.sample[0].species[0],"inf_result->sample[0]->species[0]");
	
	//mem_split(sim_result,"sim_result");
	
	//mem_split(sim_result.sample[0],"sim_result -> sample");
	
	//mem_split(sim_result.sample[0].species[0],"sim_result -> sample -> species");
	
	//mem_split(sim_result.sample[0].species[0].individual[0],"sim_result -> sample -> species -> individual");
	//mem_split(sim_result.sample[0].species[0].individual[0].ev[0],"sim_result -> sample -> species -> individual");
}


/// Estimates memory used for an object
function mem(ob)
{
	if(ob == undefined) return 0;
	
	//return roughSizeOfObject(ob);
	//obj_memory(ob)

	return JSON.stringify(ob).length;	
}


/// Estimates memory used for an object
function mem_split(ob,te)
{
	prr("\nSPLIT OBJECT: "+te);
	prr(ob);
	let mem_list=[];
	for (var key in ob){
		mem_list.push({name:key, mem:mem(ob[key]), type:typeof(ob[key])});
	}
	
	mem_list.sort( function(a, b){ return b.mem-a.mem});
	
	for(let j = 0; j < mem_list.length; j++){
		let ml = mem_list[j];
		prr(ml.name+": "+ml.mem+" "+ml.type);
	}

	prr("SPLIT DONE");
	prr("");
}


/// Checks that bscript is properly defined
function check_bscript(pro,map_store,data_file_list)
{
	prr("CHECK");

	let bsc = create_bscript("sim",map_store);
	let bsc2 = copy(pro.bscript);
	bsc2.splice(0,2);
	
	if(bsc.length != bsc2.length){
		prr("Not same length of commands");
		
		prr("bsc");
		prr(bsc);
		prr(bsc2);
		for(let i = 0; i < bsc.length; i++){
			if(bsc[i].type != bsc2[i].type){ prr(i+" wrong"); break;}
		}
		return;
	}
	
	for(let i = 0; i < bsc.length; i++){
		let sc = bsc[i];
	
		if(i>= bsc2.length){
			error("length wrong");
			return;
		}
		
		let sc2 = bsc2[i];
	
		if(sc.type != sc2.type){
			prr("Type is wrong:"+i+" "+sc.type+" "+sc2.type);
			prr("bsc");
			prr(bsc);
			prr(bsc2);
		}
		
		if(sc.tags.length != sc2.tags.length){
			prr("# tags different");
			prr(sc);
			prr(sc2);
		}
		else{
			for(let j = 0; j < sc.tags.length; j++){
				let ta = sc.tags[j];
				
				let k = 0; while(k < sc2.tags.length && sc2.tags[k].name != ta.name) k++;
				if(k == sc2.tags.length){
					prr("Tag could not be found");
					prr(ta);
					prr( sc2.tags[j]);
				}
				else{
					let ta2 = sc2.tags[k];
				
					let va = ta.value;
					let va2 = ta2.value;
					if(va2.ref != undefined){
						va2 = data_file_list[va2.ref].te;
					}
					
					if(va != va2){
						if(!(!isNaN(va) && !isNaN(va2) && Number(va) == Number(va2))){ 
							prr("tag value different");
							prr(va);
							prr(va2);
							prr(ta);
							prr(ta2);
							prr(sc);
							prr(sc2);
							prr("bsc");
							prr(bsc);
							prr(bsc2);
							if(typeof va == "string" && typeof va2 == "string"){
								for(let i = 0; i < va.length; i++){
									if(va.substr(i,1) != va2.substr(i,1)){
										prr(i+" wrong");
										prr(va.substr(0,i));
										break;
									}
								}
							}
						}
					}
				}
			}
		}
	}
}


/// Removes all but one chain
function leave_only_one_chain(pro)
{
	let bnew=[];
	for(let i = 0; i < pro.bscript.length; i++){
		let fl = false;
		let li = pro.bscript[i];
		for(let j = 0; j < li.tags.length; j++){
			//if(li.tags[j].name == "chain" && li.tags[j].value != "1") fl = true;
			if(li.tags[j].name == "chain" && li.tags[j].value != "1" && Number(li.tags[j].value) > 6) fl = true;
		}
		if(fl == false) bnew.push(li);
	}
	pro.bscript = bnew;
}

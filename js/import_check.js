"use strict";
// Functions which check that the data being added into the model is valid

/// imports an equation value
function import_eqn_value(eqn,te,op)
{
	te = te.replace(/\*/g,"×");
	
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
		let sp = model.species[p];
		
		for(let loop = 0; loop < 3; loop++){
			let source, type;
			switch(loop){
			case 0: source = sp.sim_source; type = "sim"; break;
			case 1: source = sp.inf_source; type = "inf"; break;
			case 2: source = sp.ppc_source; type = "ppc"; break;
			}
			
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
	
	for(let th = 0; th < model.param.length; th++){
		let par = model.param[th];
		if(flag[th] != true && !is_in_obsmodel(par) && par.type != "derive_param"){
			if(par.name != dist_matrix_name){
				alert_noline("Parameter "+par.full_name+" is not specified as a 'param'");
			}
		}
	
		if(par.type == unset_type){
			alert_noline("Parameter "+par.full_name+" does not have a type")   
		}
	}
}


/// Determines if text string is for a file
function is_file(te)
{
	if(te == "") return false;
	
	if(te.te != undefined) return true;
	
	let k = te.length-1;
	let kmin = k-10;
	while(k > 0 && k > kmin && te.substr(k,1) != ".") k--;
	if(k == 0 || k == kmin) return false;
	let end = te.substr(k).toLowerCase();
	if(end != ".csv" && end != ".txt" && end != ".tsv" && end != ".geojson") return false;

	return true;
}

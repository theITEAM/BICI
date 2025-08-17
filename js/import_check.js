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
		
			let source, type;
			switch(loop){
			case 0: source = sp.sim_source; type = "sim"; break;
			case 1: source = sp.inf_source; type = "inf"; break;
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
		if(flag[th] != true && !is_in_obsmodel(par) && par.type != "derive_param"){
			if(!par.dist_mat && !par.iden_mat){
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
			

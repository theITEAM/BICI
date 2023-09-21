/// Checks that the data being added into the model is valid

"use strict";

/// Checks that 
function edit_param_check_error(error_mess_only)
{
	return false;
}

/// Before adding a data source check that everything is correct
function data_source_check_error(error_mess_only)
{
	let so = inter.edit_source;

	let tab = so.table;
	if(tab != undefined){
		// Checks that elements are correctly specified
		for(let c = 0; c < tab.ncol; c++){
			for(let r = 0; r < tab.nrow; r++){
				let warn = check_element(tab.ele[r][c],c);
			
				if(warn != ""){
					if(error_mess_only == true){
						if(warn.substr(warn.length-1,1) == ".") warn = warn.substr(0,warn.length-1);
						
						inter.help = {title:"Error with data table", te:warn+" (col '"+tab.heading[c]+"', row "+(r+2)+")."};
					}
					else{
						select_table_elelent(r,c);
						inter.bubble.error_warning = warn;
					}
					return true;
				}
			}
		}

		// Checks not adding individual twice
		if(so.type == "Add Ind." || so.type == "Remove Ind."){ 
			let list=[];
			for(let r = 0; r < tab.nrow; r++) list.push(tab.ele[r][0]);
			list.sort();
			for(let r = 0; r < list.length-1; r++){
				if(list[r] == list[r+1]){
					let te; if(so.type == "Add Ind.") te = "added"; else te = "removed";
					
					inter.help = {title:"Error with data table", te:"Individual '"+list[r]+"' is "+te+" more than once."};
					return true;
				}
			}
		}
		
		if((so.type == "Transition" || so.type == "Source" || so.type == "Sink") && 
				so.spec.time_radio.value == "File"){
			for(let r = 0; r < tab.nrow; r++){
				if(Number(tab.ele[r][2]) >= Number(tab.ele[r][3])){				
					inter.help = {title:"Error with data table", te:"For the observations the 'start' time must be before the 'end' time (row "+(r+2)+")."};
					return true;
				}
			}
		}
		
		if(so.type == "Diag. Test"){
			for(let r = 0; r < tab.nrow; r++){
				let res = tab.ele[r][2];
				
				let neg_res = so.spec.neg_result;
				let pos_res = so.spec.pos_result;
				
				if(res != neg_res && res != pos_res){
					inter.help = {title:"Error with data table", te:"The test result '"+res+"' is neither negative '"+neg_res+"' or positive '"+pos_res+"' (col '"+tab.heading[2]+"', row "+(r+2)+")."};
					return true;
				}
			}
		}
		
		if(so.type == "Set Traps"){
			let list=[];
			for(let r = 0; r < tab.nrow; r++) list.push(tab.ele[r][0]);
			list.sort();
			for(let r = 0; r < list.length-1; r++){
				if(list[r] == list[r+1]){	
					inter.help = {title:"Error with data table", te:"In column '"+tab.heading[0]+"' the trapping event '"+list[r]+"' is set more than once."};
					return true;
				}
			}
		}
	}
	
	let desc = check_data_source_valid(so);
	if(so.error == true){
		if(desc == undefined) desc = so.error_mess;
		inter.help = {title:"Error with data table", te:desc};
		return true;
	}
	
	return false;
}


/// Checks the value of an element is correct
function check_element(te,c)
{
	let col = inter.edit_source.load_col[c];
	switch(col.type){
	case "float":
		{
			let num = Number(te);
			if(isNaN(num)) return "Must be a number";
		}
		break;
		
	case "time_float":
		{
			if(te != "start" && te != "end"){
				let num = Number(te);
				if(isNaN(num)) return "Must be a number";
			}
		}
		break;
		
	case "time_trans_float":
		{
			if(te != "start" && te != "end" && te != "no"){
				let num = Number(te);
				if(isNaN(num)) return "Must be a number";
			}
		}
		break;
	
	case "pos_float":
		{
			let num = Number(te);
			if(isNaN(num)) return "Must be a number";
			if(num <= 0) return "Must be positive";
		}
		break;
		
	case "pos_zero_float":
		{
			let num = Number(te);
			if(isNaN(num)) return "Must be a number";
			if(num < 0) return "Must be positive or zero";
		}
		break;
		
	case "positive int":
		{
			let num = Number(te);
			if(isNaN(num)) return "The Must be a number";
			if(num != Math.round(num)) return "Must be an integer";	
			if(num < 0) return "Must be non-negative";
		}
		break;

	case "positive float/percent":
		{
			if(te.length > 0 && te.substr(te.length-1,1) == "%"){
				te = te.substr(0,te.length-1);
			}
			
			let num = Number(te);
			if(isNaN(num)) return "Must be a number";
			if(num < 0) return "Must be non-negative";
		}
		break;
	
	case "text":
		if(te == "") return "Element empty";
	
		if(col.cl != undefined) error("SHOULD BE 'compartment'");
		break;
		
	case "compartment": 
		{  // This allows 'S', 'E' etc...
			let claa =  model.species[col.p].cla[col.cl];
			let comp = claa.comp;
		
			if(find(comp,"name",te) == undefined){
				return "Compartment '"+te+"' not in classification '"+claa.name+"'";
			}
		}
		break;
		
	case "comp_all":
		{
			let result = find_comp_from_name(te,model.species[col.p].name);
			if(typeof result == 'string') return result;
		}
		break;
		
	case "comp_source_sink": 
		{  // Allows for compartment as well as potentially source or sink
			let claa = model.species[col.p].cla[col.cl];
			let comp = claa.comp;
		
			if(te != "+" && te != "-"){
				if(find(comp,"name",te) == undefined){
					return "Compartment '"+te+"' not in classification '"+claa.name+"'";
				}
			}
		}
		break;
		
	case "compartment_prob": 
		{ // This allows 'S', 'E', 'S|E', 'S:0.4|E:0.6' etc...
			let claa = model.species[col.p].cla[col.cl];
			let comp = claa.comp;
		
			let spl = te.split("|");
		
			let syntax_error = false;
				
			let colon_flag = false, normal_flag = false;
			for(let i = 0; i < spl.length; i++){
				if(spl[i] == ""){ syntax_error = true; break;}
				else{
					let spl2 = spl[i].split(":");
					
					if(spl2[0] == ""){ syntax_error = true; break;}
					
					if(find(comp,"name",spl2[0]) == undefined){
						return "Compartment '"+spl2[0]+"' not in classification '"+claa.name+"'";
					}
					
					if(spl2.length == 1){
						normal_flag = true;
					}							
					else{
						colon_flag = true; 
						if(spl2[1] == ""){ syntax_error = true; break;}
						
						let valid = check_eqn_valid(spl2[1]);
						if(valid != "success"){
							return "Error with equation'"+spl2[1]+"': "+valid;
						}
					}
					
					if(spl2.length > 2){ syntax_error = true; break;}
				}
			}
				
			if(normal_flag == true && colon_flag == true) syntax_error = true;

			if(syntax_error == true){
				return "There is a syntax error in table element.";
			}
		}
		break;
		
	case "multi_comp":
		{
			let claa = model.species[col.p].cla[col.cl];
			let comp = claa.comp;
			let spl = te.split("|");
			
			for(let j = 0; j < spl.length; j++){
				if(find(comp,"name",spl[j]) == undefined){
					return "Value '"+spl[j]+"' not in classification '"+claa.name+"'";
				}
			}
		}
		break;
		
	case "diag test":
		if(te != col.pos_res && te != col.neg_res){
			return "Value '"+te+"' is neither '"+col.pos_res+"' or '"+col.neg_res+"'";
		}
		break;
	
	case "colour":
		if(isColor(te) == false) return "Value '"+te+"' is not a colour";
		break;
		
	case "eqn":
		{
			let p = col.p;
			let cl = col.cl;
			let eqn = create_equation(te,"comp_prob",p,cl);
			if(eqn.warn.length > 0) return eqn.warn[0].te;
		}
		break;
	
	case "lng":
		{
			let num = Number(te);
			if(isNaN(num)) return "Must be a number";
			if(num < -180 || num > 180) return "Must be in range -180 - 180";
		}
		break;
		
	case "lat":
		{
			let num = Number(te);
			if(isNaN(num)) return "Must be a number";
			if(num < -90 || num > 90) return "Must be in range -90 - 90";
		}
		break;
		
	case "prior":
		let pri = convert_text_to_prior(te,col.pri_pos);
		if(typeof pri == 'string') return pri;
		break;
		
	default: error("Option not recognised 61"+col.type); break;
	}
	
	return "";
}

/// Checks is an equation is valid
function check_eqn_valid(te)
{
	let eqn = create_equation(te,"test");		
	if(eqn.warn.length > 0) return eqn.warn[0].te;
	
	return "success";
}


/// Determines if a string is a valid colour
function isColor(strColor)
{
	var s = new Option().style;
	s.color = strColor;
	return s.color.length > 0;
}


/// imports an equation value
function import_eqn_value(eqn,te,op)
{
	eqn.te = te;
	extract_equation_properties(eqn);
	if(eqn.warn.length > 0){
		alert_import("Equation '"+te+"' – "+ eqn.warn[0].te); return;
	}
	
	return "success";
}


/// Checks if a string is an equation with a single parameter or a number
function is_eqn(te,tag,op)
{
	let eqn = create_equation(te,"test");
	
	if(eqn.warn.length > 0){
		alert_import("For '"+tag+"' — "+eqn.warn[0].te); 
		return false;
	}
	
	if(eqn.mode == "param_only"){
		if(eqn.ind_eff.length > 0){
			alert_import("'"+tag+"' cannot contain individual effects"); 
			return false;
		}
		
		if(eqn.param.length > 1){
			alert_import("'"+tag+"' can only be constant or contain one parameter"); 
			return false;
		}
	}
	
	if(eqn.param.length == 0){
		/*
		if(isNaN(te)){
			alert_import("'"+tag+"' must be a number or a parameter"); 
			return false;
		}
		*/
		
		if(op.positive == true){
			if(Number(te) <= 0){
				alert_import("'"+tag+"' must be a positive number"); 
				return false;
			}
		}
		
		if(op.zero_one_range == true){
			if(Number(te) < 0 || Number(te) > 1){
				alert_import("'"+tag+"' must be in the range 0 - 1"); 
				return false;
			}
		}
	}

	return true;
}


/// Check for initial pop
function check_initial_pop_error(end)
{
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		
		for(let loop = 0; loop < 2; loop++){
			let source, type;
			if(loop == 0){ source = sp.sim_source; type = "sim";}
			else{ source = sp.inf_source; type = "inf";}
				
			let ninitpopprior = 0, ninitpop = 0;
			for(let i = 0; i < source.length; i++){
				let so = source[i];
				
				switch(so.type){
				case "Init. Pop. Prior": ninitpopprior++; break;
				case "Init. Pop.": ninitpop++; break;
				}
			}
			
			switch(type){
			case "sim":
				if(ninitpopprior > 0){
					return "The initial population prior should not be set for species '"+sp.name+"'";
				}
				
				/*
				if(end == true){
					if(ninitpop == 0){
						return "The initial population should be set for species '"+sp.name+"'";
					}
				}
				*/
				
				if(ninitpop > 1){
					return "Only one initial population should be set for species '"+sp.name+"'";
				}
				break;
			
			case "inf":
				if(ninitpopprior > 1){
					return "Only one initial population prior should be set for species '"+sp.name+"'";
				}
				
				if(ninitpop > 1){
					return "Only one initial population should be set for species '"+sp.name+"'";
				}
				
				/*
				if(end == true){
					if(ninitpopprior == 0 && ninitpop == 0){
						return "Either 'Init. Pop.' or 'Init. Pop. Prior' must be set for species '"+sp.name+"'";
					}
				}
				*/
				
				if(ninitpopprior == 1 && ninitpop == 1){
					return "'Init. Pop.' and 'Init. Pop. Prior' cannot both be set for species '"+sp.name+"'";
				}
				break;
			}
		}
	}
}


/// Checks that parameter specification are complete
function check_param_complete()
{
	let param = copy(model.param);
	
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
				return;
			}
		
			/// Checks the prior are within those possible
			if(par2.prior.type.te != select_str){
				if(option_error("prior",par2.prior.type.te,par2.pri_pos) == true){
					alert_line(inter.help.st,par.import_line);
					return;
				}
			}
		
			flag[th2] = true;
		}
	}
	
	for(let th = 0; th < model.param.length; th++){
		let par = model.param[th];
		if(flag[th] != true && !is_in_obsmodel(par) && par.type != "derive_param"){
			alert_noline("Parameter "+par.full_name+" is not specitified as a 'param'")   
			return;
		}
	
		if(par.type == unset_type){
			alert_noline("Parameter "+par.full_name+" does not have a type")   
			return;
		}
	}
}


/// Checks if there is an error in importing
function check_import_error(store) 
{			
	if(inter.help.te != undefined){
		load_BICI_file(store);
		generate_screen();
		return true;
	}
	
	return false;
}
		
		


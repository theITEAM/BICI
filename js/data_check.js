"use strict";
// Functions which check the data inputs

/// Checks if an individual data source is valid
function data_source_check_error(out_type,so)
{	
	let p = so.info.p;
	
	let sp = model.species[p];

	if(so.info.siminf == "gen") sp = model.sim_res.plot_filter.species[p];
	so.error = false;
	
	let tab = so.table;
	
	if(tab != undefined){                            // Checks that elements are correctly specified
		/* // This was turned off becase it lead to error e.g. "E | I"
		// If a compartment then corrects name (e.g. space -> "-") 
		if(so.load_col){
			for(let c = 0; c < so.load_col.length; c++){
				switch(so.load_col[c].type){
				case "compartment": case "comp_all": 
				case "compartment_prob": case "compartment_prob_na":
				case "comp_source_sink": case "comptext":
					for(let r = 0; r < tab.nrow; r++){	
						tab.ele[r][c] = make_comp_name(tab.ele[r][c]);
					}
					break;
				}
			}
		}
		*/
	
		for(let c = 0; c < tab.ncol; c++){
			for(let r = 0; r < tab.nrow; r++){
				let warn = check_element(tab.ele[r][c],c,so);			
				if(warn != "") data_error(warn,out_type,so,{r:r,c:c});
			}
		}
	}

	
	switch(so.type){
	case "Pop. Trans.":
		for(let r = 0; r < tab.nrow; r++){
			let sta_str = tab.ele[r][0], end_str = tab.ele[r][1];
			let sta = Number(sta_str), end = Number(end_str);
			if(sta >= end){
				data_error("The start time '"+sta_str+"' must be less than the end time '"+end_str+"'.",out_type,so,{r:r});
			}
		}	
		break;
	}

	// Checks for repeated elements
	switch(so.type){
	case "Comp File": case "Comp File Colour": case "Comp File Pos": case "Comp File Pos Colour":
	case "CompMap":
		{	 
			let c = 0; if(so.type == "CompMap") c = 1;
		
			// Checks to see if compartments have the same name
			let list = []
			for(let r = 0; r < tab.nrow; r++){	
				list.push(tab.ele[r][c]);
			}
			
			let repeat = check_repeat(list);
		
			if(repeat != ""){
				so.warning = "Note, repeated compartments are ignored: "+repeat+".";
			}
		}
		break;
		
	case "Trans File": case "Trans File Pos":
		break;

	case "Add Pop.": case "Remove Pop.":
		break;
		
	case "Add Ind.": case "Remove Ind.":
		{
			// Checks not adding individual twice
			let list=[];
			for(let r = 0; r < tab.nrow; r++) list.push(tab.ele[r][0]);
			
			let repeat = check_repeat(list);
		
			if(repeat != ""){
				let te; if(so.type == "Add Ind.") te = "added"; else te = "removed";
				
				data_error("The following individuals are "+te+" more than once:"+repeat,out_type,so);
			}
		}
		break;
		
	case "Transition": 
		if(so.spec.time_radio.value == "File"){
			for(let r = 0; r < tab.nrow; r++){
				if(Number(tab.ele[r][2]) >= Number(tab.ele[r][3])){		
					data_error("For the observations the 'start' time must be before the 'end' time (row "+(r+2)+").",out_type,so);
				}
			}
		}
		break;
		
	case "Diag. Test":
		for(let r = 0; r < tab.nrow; r++){
			let res = tab.ele[r][2];
			
			let neg_res = so.spec.neg_result;
			let pos_res = so.spec.pos_result;
			
			if(res != neg_res && res != pos_res && res != "NA"){
				data_error("The test result '"+res+"' is neither negative '"+neg_res+"' or positive '"+pos_res+"'",out_type,so,{r:r,c:2});
			}
		}
		break;
	}
	
	let mess_def = data_invalid;

	switch(so.type){
	case "Init. Pop.":
		if(so.spec.radio_dist.value == "Dist" && so.spec.radio.value == "Graphical"){
			switch(so.spec.radio2.value){
			case "Focal":
				{
					let p = so.info.p;
					let sp = model.species[p];
					
					let foc_cl = so.spec.focal.te;
					let cl = find(sp.cla,"name",foc_cl);
					if(cl == undefined) error("problem finding cl");
					else{
						for(let c = 0; c < sp.cla[cl].ncomp; c++){
							if(so.cla[cl] && so.cla[cl].comp_init_dist[c]){
								if(so.cla[cl].comp_init_dist[c].dist.type.te == select_str){
									data_error("The population distribution for <e>"+sp.cla[cl].comp[c].name+"</e> must be set.",out_type,so);
								}
							}								
						}
					}
				}
				break;
				
			case "All":
				if(so.pop_dist.type.te == select_str){
					data_error("The population distribution must be set",out_type,so);
				}
				break;
			}
		}
		break;
	}

	if(sp){
		if(sp.type == "Population"){
			switch(so.type){
			case "Add Ind.":case "Move Ind.":  case "Remove Ind.": case "Compartment":
			case "Transition": case "Diag. Test": case "Genetic": case "Ind. Eff.": case "Ind. Group":
				data_error("Individual-based data '"+so.type+"' cannot be added to a population-based model",out_type,so);
				break;
			}
		}
		
		if(sp.type == "Individual"){
			switch(so.type){
			case "Add Pop.":case "Remove Pop.":
				data_error("Population-based data '"+so.type+"' cannot be added to an individual-based model",out_type,so);
				break;
			}
		}
	}
	
	if(so.is_data){
		switch(so.type){
		case "Init. Pop.":
			{
				let di_type = so.spec.radio_dist.value;
				
				if(so.spec.radio.value ==  "Graphical"){
					switch(so.spec.radio2.value){
					case "Focal":
						let focal = so.spec.focal.te;
						let cl = find(sp.cla,"name",focal);

						if(cl == undefined){
							data_error("The focal value '"+focal+"' is not a classification",out_type,so);
						}
						else{
							if(so.cla.length != sp.ncla) data_error(mess_def,out_type,so,{num:1});
							else{
								for(let cl = 0; cl < sp.ncla; cl++){
									let claa = sp.cla[cl];
									let clz = so.cla[cl];

									if(claa.name != clz.cl_name_store) data_error(mess_def,out_type,so,{num:2});
									else{
										switch(di_type){
										case "Fixed":
											if(clz.comp_init_pop.length != claa.ncomp) data_error(mess_def,out_type,so,{num:3});
											else{
												for(let c = 0; c < claa.ncomp; c++){
													if(clz.comp_init_pop[c].comp_name_store != claa.comp[c].name){
														data_error(mess_def,out_type,so,{num:4});
													}
												}
											}
											break;
											
										case "Dist":
											if(clz.comp_init_dist.length != claa.ncomp) data_error(mess_def,out_type,so,{num:5});
											else{
												for(let c = 0; c < claa.ncomp; c++){
													if(clz.comp_init_dist[c].comp_name_store != claa.comp[c].name){
														data_error(mess_def,out_type,so,{num:6});;
													}
												}
											}
											break;
										}
									}
								}
							}
						}
						break;

					case "All":	
						{
							let glob_comp = get_glob_comp(p);
							if(glob_comp.length != so.glob_comp.length) data_error(mess_def,out_type,so,{num:7});
							else{
								for(let j = 0; j < glob_comp.length; j++){
									let gc = glob_comp[j].cla;
									let gc2 = so.glob_comp[j].cla;
									if(gc.length != gc2.length) data_error(mess_def,out_type,so,{num:8});
									else{
										for(let cl = 0; cl < gc.length; cl++){
											if(gc[cl] != gc2[cl]) data_error(mess_def,out_type,so,{num:9});
										}
									}
								}
							}
						}
						break;
					
					default: error("Option not recognised 200"); break;
					}
				}
			}	
			break;
	
		case "Add Pop.": case "Remove Pop.":
			{
				let tab = so.table;
				if(tab.ncol != sp.ncla+2) data_error(mess_def,out_type,so,{num:10});
				for(let cl = 0; cl < sp.ncla; cl++){
					let claa = sp.cla[cl];

					if(tab.heading[cl+1] != claa.name) data_error(mess_def,out_type,so,{num:11});
					else{
					
					}
				}
			}
			break;
			
		case "Add Ind.":
			{
				let tab = so.table;
				if(tab.ncol != sp.ncla+2) data_error(mess_def,out_type,so,{num:10});
				else{
					for(let cl = 0; cl < sp.ncla; cl++){
						let claa = sp.cla[cl];

						if(tab.heading[cl+2] != claa.name) data_error(mess_def,out_type,so,{num:11});
						else{
						
						}
					}
				}
			}
			break;

		case "Move Ind.":
			{
				let cl = find(sp.cla,"name",so.spec.cl_drop.te);
				if(cl == undefined) data_error(mess_def,out_type,so,{num:12});
				else{
					let claa = sp.cla[cl];

					let tab = so.table;
					for(let r = 0; r < tab.nrow; r++){
						let to = tab.ele[r][2];
						if(hash_find(claa.hash_comp,to) == undefined){
							data_error("The value '"+to+"' is not a compartment in classification '"+claa.name+"'",out_type,so,{r:r,c:3,num:13});
						}
					}
				}
			}
			break;

		case "Compartment":
			{
				let cl = find(sp.cla,"name",so.spec.cl_drop.te);
				if(cl == undefined) data_error(mess_def,out_type,so,{num:14});
				else{
					let claa = sp.cla[cl];

					let tab = so.table;
					if(tab.heading[2] != claa.name) data_error(mess_def,out_type,so,{num:15});
					else{
					}
				}
			}
			break;
			
		case "Diag. Test":
			{
				let cb = so.spec.check_box;
				let cl = find(sp.cla,"name",cb.name);
				if(cl == undefined) data_error(mess_def,out_type,so,{num:16});
				else{
					let claa = sp.cla[cl];
					if(cb.value.length != claa.ncomp) data_error(mess_def,out_type,so,{num:17});
					else{
						for(let c = 0; c < claa.ncomp; c++){
							if(cb.value[c].comp_name_store != claa.comp[c].name) data_error(mess_def,out_type,so,{num:18});
						}
					}
				}
			}
			break;
			
		case "Ind. Eff.":
			{
				let tab = so.table;
				for(let r = 0; r < tab.nrow; r++){
					let ie = tab.ele[r][1];
					if(isNaN(ie)){
						data_error("The value '"+ie+"' is not a number",mess_def,out_type,so,{r:r,c:3,num:19});
					}
				}
			}
			break;

		case "Ind. Group":
			{
			}
			break;
		
		case "Population":
			{
				let filt = so.spec.filter;
		
				if(filt.cla.length != sp.ncla) data_error(mess_def,out_type,so,{num:21});
				else{
					for(let cl = 0; cl < sp.ncla; cl++){
						let claa = sp.cla[cl];
						let clz = filt.cla[cl];
						if(claa.name != clz.cl_name_store) data_error(mess_def,out_type,so,{num:22});
						else{
							if(clz.comp.length != claa.ncomp) data_error(mess_def,out_type,so,{num:122});
							else{
								for(let c = 0; c < claa.ncomp; c++){
									if(clz.comp[c].comp_name_store != claa.comp[c].name) data_error(mess_def,out_type,so,{num:23});
								}
							}
						}
					}
				}
			}
			break;


		case "Transition": 
		case "Pop. Trans.":
			{
				let cl = find(sp.cla,"name",so.spec.cl_drop.te);
				if(cl == undefined) data_error(mess_def,out_type,so,{num:24});
				else{
					let claa = sp.cla[cl];

					let filt = so.spec.filter;

					if(filt.tra.length != claa.ntra) data_error(mess_def,out_type,so,{num:25}); 
					else{
						let flag = false;
						for(let j = 0; j < claa.ntra; j++){
							let tra = filt.tra[j];
							if(tra.tra_name_store != claa.tra[j].name) data_error(mess_def,out_type,so,{num:26});
							
							switch(filt.trans_obs_model.value){
							case "off": if(tra.check == true) flag = true; break;
							case "on": 
								if(!tra.prob_eqn) data_error(mess_def,out_type,so,{num:26});
								else{
									if(tra.prob_eqn.te != "0") flag = true; 
								}
								break;
							default: data_error(mess_def,out_type,so,{num:27}); break;
							}
						}
						if(flag == false) data_error(mess_def,out_type,so,{num:28});
					}

					if(filt.cla.length != sp.ncla){
						data_error(mess_def,out_type,so,{num:29});
					}
					else{
						for(let cl2 = 0; cl2 < sp.ncla; cl2++){
							if(cl2 != cl){
								let claa2 = sp.cla[cl2];
								let clz = filt.cla[cl2];
								if(clz.cl_name_store != claa2.name) data_error(mess_def,out_type,so,{num:30});
								else{	
									if(claa2.ncomp != clz.comp.length) data_error(mess_def,out_type,so,{num:31});
									else{
										for(let c = 0; c < claa2.ncomp; c++){
											if(clz.comp[c].comp_name_store != claa2.comp[c].name){
												data_error(mess_def,out_type,so,{num:32});
											}	
										}		
									}
								}
							}
						}
					}
				}
			}
			break;

		case "Remove Ind.": break;
		case "Genetic": break;	

		default: error("Source not checked error: "+so.type); break;
		}
	}
}


/// Checks that the time column is within the range
function check_data_time(out_type)
{
	if(out_type != "sim" && out_type != "inf") return;

	let details;
	switch(out_type){
	case "sim": details = model.sim_details; break;
	case "inf": details = model.inf_details; break;
	}
		
	let t_start = Number(details.t_start);
	let t_end = Number(details.t_end);
	
	let mess_def = "Data source time is invalid";

	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		
		let source;
		switch(out_type){
		case "sim": source = sp.sim_source; break;
		case "inf": source = sp.inf_source; break;
		}
		
		for(let k = 0; k < source.length; k++){
			let so = source[k];
			let tab = so.table;
			if(tab){
				for(let c = 0; c < so.load_col.length; c++){
					if(so.load_col[c].type == "time_float"){
						for(let r = 0; r < tab.nrow; r++){
							let te = tab.ele[r][c];
							if(te != "start" && te != "end"){
								if(isNaN(te)){
									data_error(mess_def,"warn",so,{c:c, r:r, num:130});
									break;
								}
								else{
									let t = Number(te);
									if(t < t_start || t > t_end){
										data_error("The time '"+t+"' must be between the start and end times","warn",so,{minor:true, c:c, r:r, num:130});
										break;
									}	

									if(so.type == "Add Ind." && t == t_end){
										data_error("When adding individuals the time '"+t+"' must be before the end time","warn",so,{minor:true, c:c, r:r, num:130});
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


/// Checks the value of an element is correct
function check_element(te,c,so)
{
	let col = so.load_col[c];	
	
	switch(col.type){
	case "float":
		{
			let num = Number(te);
			if(isNaN(num)) return "The value '"+te+"' must be a number";
		}
		break;
		
	case "time_float":
		{
			if(te != "start" && te != "end"){
				let num = Number(te);
				if(isNaN(num)) return "The value '"+te+"' must be a number";
			}
		}
		break;
		
	case "time_trans_float":
		{
			if(te != "start" && te != "end" && te != "no"){
				let num = Number(te);
				if(isNaN(num)) return "The value '"+te+"' must be a number";
			}
		}
		break;
	
	case "pos_float":
		{
			let num = Number(te);
			if(isNaN(num)) return "The value '"+te+"' must be a number";
			if(num <= 0) return "The value '"+te+"' must be positive";
		}
		break;
		
	case "pos_zero_float":
		{
			let num = Number(te);
			if(isNaN(num)) return "The value '"+te+"' must be a number";
			if(num < 0) return "The value '"+te+"' must be positive or zero";
		}
		break;

	case "zero_one_float":
		{
			let num = Number(te);
			if(isNaN(num)) return "The value '"+te+"' must be a number";
			if(num <= 0 || num >= 1) return "The value '"+te+"' must be between zero and one";
		}
		break;
		
	case "positive int":
		{
			let num = Number(te);
			if(isNaN(num)) return "The value '"+te+"' must be a number";
			if(num != Math.round(num)) return "The value '"+te+"' must be an integer";	
			if(num < 0) return "The value '"+te+"' must be non-negative";
		}
		break;
		
	case "positive nonzero int":
		{
			let num = Number(te);
			if(isNaN(num)) return "The value '"+te+"' must be a number";
			if(num != Math.round(num)) return "The value '"+te+"' must be an integer";	
			if(num <= 0) return "The value '"+te+"' must be positive";
		}
		break;

	case "positive float/percent":
		{
			if(te.length > 0 && te.substr(te.length-1,1) == "%"){
				te = te.substr(0,te.length-1);
			}
			
			let num = Number(te);
			if(isNaN(num)) return "The value '"+te+"' must be a number";
			if(num < 0) return "The value '"+te+"' must be non-negative";
		}
		break;
	
	case "text":
		if(te == ""){
			return "Element empty";
		}
		if(te.trim().includes(" ")){ return "The value '"+te+"' should not include a space. Names in BICI are not allowed spaces.";}
		if(col.cl != undefined) error("SHOULD BE 'compartment'");
		break;
	
	case "comptext":
		if(te == "") return "Element empty";	
		let compwarn = model.check_comp_name(te);
		if(compwarn != undefined) return compwarn;
		break;
		
	case "knot":
		if(find_in(col.knot,te) == undefined){
			return "The value '"+te+"' is not a knot time";
		}
		break;
		
	case "compartment": 
		{  // This allows 'S', 'E' etc...
			let claa = model.species[col.p].cla[col.cl];
			if(hash_find(claa.hash_comp,te) == undefined){
				return "Compartment '"+te+"' is not in classification '"+claa.name+"'";
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
			te = te.trim();
			if(te != "+" && te != "-"){
				if(hash_find(claa.hash_comp,te) == undefined){
					return "Compartment '"+te+"' is not in classification '"+claa.name+"'";
				}
			}
		}
		break;
		
	case "compartment_prob": case "compartment_prob_na": 
		{ 
			// This allows 'S', 'E', 'S|E', 'S:0.4|E:0.6' etc...
			// "compartment_prob_na" also allow '!' for "not alive"
			if(te != missing_str && !(col.type == "compartment_prob_na" && te == not_alive_str)){
				let claa = model.species[col.p].cla[col.cl];
				let hash_comp = claa.hash_comp;
			
				let spl = split_with_bracket(te,"|");
			
				let syntax_error = false;
					
				let colon_flag = false, normal_flag = false;
				let list=[];
				for(let i = 0; i < spl.length; i++){
					if(spl[i] == ""){ syntax_error = true; break;}
					else{
						let spl2 = spl[i].split(":");
						
						if(spl2[0] == ""){ syntax_error = true; break;}
						
						let k = hash_find(hash_comp,spl2[0]);
						if(k == undefined){
							return "In table element '"+te+"', compartment '"+spl2[0]+"' is not in classification '"+claa.name+"'";
						}
						if(find_in(list,k)){
							return "In table element '"+te+"', compartment '"+spl2[0]+"' is repeated.";
						}
						
						list.push(k);
						
						if(spl2.length == 1){
							normal_flag = true;
						}							
						else{
							colon_flag = true; 
							if(spl2[1] == ""){ syntax_error = true; break;}
							
							let valid = check_eqn_valid(spl2[1]);
							if(valid.err == true){
								return "In table element '"+te+"', error with equation '"+spl2[1]+"': "+valid.msg;
							}
						}
						
						if(spl2.length > 2){ syntax_error = true; break;}
					}
				}
					
				if(normal_flag == true && colon_flag == true){
					return "In table element '"+te+"', either all compartments should use a probability or all not (they cannot be mixed). For example, 'S:0.5|I:0.5' uses probabilites and 'S|I' does not.";
				}
				
				if(syntax_error == true){
					return "There is a syntax error in table element with value '"+te+"'.";
				}
			}
		}
		break;
		
	case "diag test":
		if(te != col.pos_res && te != col.neg_res && te != "NA"){
			return "The value '"+te+"' is neither '"+col.pos_res+"' or '"+col.neg_res+"'";
		}
		break;
	
	case "colour":
		if(isColor(te) == false) return "The value '"+te+"' is not a colour";
		break;
		
	case "trans_value":
		{
			let trans_def = extract_trans_def(te);
			if(trans_def == undefined) return "Problem extracing transition information from '"+te+"' (e.g. should have a a format similar to 'exp(rate:r)')";
		}
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
			if(isNaN(num)) return "The value '"+te+"' must be a number";
			if(num < -180 || num > 180) return "The value '"+te+"' must be in the range -180 - 180";
		}
		break;
		
	case "lat":
		{
			let num = Number(te);
			if(isNaN(num)) return "The value '"+te+"' must be a number";
			if(num < -90 || num > 90) return "The value '"+te+"' must be in the range -90 - 90";
		}
		break;
		
	case "prior":
		{
			let pri = convert_text_to_prior(te,col.pri_pos);
			if(pri.err == true) return pri.msg;
		}
		break;
		
	case "dist":
		{
			let pri = convert_text_to_prior(te,col.pri_pos,true);
			if(pri.err == true) return pri.msg;
		}
		break;
	
	case "boundary":
		break;
		
	case "prior/dirichlet":
		if(isNaN(te)){
			let pri = convert_text_to_prior(te,col.pri_pos);
			if(pri.err == true) return pri.msg;
		}
		break;
		
	default: error("Option not recognised 61"+col.type); break;
	}
	
	return "";
}


/// Checks to see if elements are repeated
function check_repeat(list)
{
	list.sort();
	let repeat="";
	for(let i = 0; i < list.length-1; i++){
		if(list[i] == list[i+1]){
			if(repeat != "") repeat += ", ";
			repeat += list[i];
			while(i < list.length-1 && list[i] == list[i+1]) i++;
		}
	}
	return repeat;
}


/// Adds a warning to the list
function add_warning(wa)
{
	let i = 0; while(i < model.warn.length && !(model.warn[i].mess == wa.mess && model.warn[i].mess2 == wa.mess2)) i++;
	
	if(i == model.warn.length){
		model.warn.push(wa);
	}
}	
	
	
/// Output a data error warning
function data_error(warn,out_type,so,op)
{
	if(false){ error("ERR"); error(so);}

	if(so.error == true) return;
	
	let minor_warn = false;
	
	if(op != undefined && op.minor) minor_warn = true;
	else so.error = true;
	
	if(out_type == "warn"){
		let de = so.type;
		if(so.desc != undefined) de += " "+so.desc;
	
		if(de.length > 20) de = de.substr(0,20)+"...";
	
		warn = "For '"+de+"': "+warn;
	}
	else{
		if(so.table != undefined){
			warn = in_file_text(so.table.filename)+": "+warn;
		}
	}
	
	if(warn.substr(warn.length-1,1) == ".") warn = warn.substr(0,warn.length-1);	
	if(op != undefined && op.r != undefined){
		let tab = so.table;
		warn += " (";
		if(op.c) warn += "col '"+tab.heading[op.c]+"', ";
		warn += "row "+(op.r+2)+").";
	}
			
	switch(out_type){
	case "warn":
		{
			let ci = so.check_info;	
			if(minor_warn){
				add_warning({mess:"Data error", mess2:warn, siminf:ci.siminf, p:ci.p, warn_type:"DataProblem"});
			}
		}
		break;
				
	case "import":
		alert_import(warn);
		break;
	
	case "alert warn":
		throw({type:"AlertP", te:warn});
		break;
		
	case "worker":
		if(!op || !op.r) throw({type:"AlertP", te:warn});
		throw({type:"Select Element", r:op.r, c:op.c, warn:warn});
		break;
	
	default:
		error("Should not be def");
		break;
	}
}


/// Goes through all data soruces and check if valid given the current model
function check_data_valid_all(type)
{
	let list = [];
	if(type == "sim" || type == "siminf") list.push("sim");
	if(type == "inf" || type == "siminf") list.push("inf");
	if(type == "ppc" || type == "siminf") list.push("ppc");

	for(let loop = 0; loop < list.length; loop++){
		let siminf = list[loop];
		
		let spec;
		if(siminf == "ppc"){
			let pf = model.inf_res.plot_filter;
			if(pf) spec = pf.species;
		}
		else{
			spec = model.species;
		}
		
		if(spec){
			for(let p = 0; p < spec.length; p++){
				let sp = spec[p];
		
				let source;
				switch(siminf){
				case "sim": source = sp.sim_source; break;
				case "inf": source = sp.inf_source; break;
				case "ppc": source = sp.ppc_source; break;
				}

				for(let i = 0; i < source.length; i++){
					let so = source[i];
					so.check_info = {siminf:list[loop], p:p, i:i}; 
					data_source_check_error("warn",so);
					if(so.error != true) add_source_description(so);
				}
			}
		}
	}
}

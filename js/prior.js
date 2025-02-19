"use strict";
// Functions related to the prior

/// Adds the button for the page in which parameter priors are input
function add_param_prior_buts(lay)
{
	let cx = corner.x;
	let cy = corner.y;
	
	cy = lay.add_title("Parameter priors",cx,cy,{te:editprior_text2});
	
	cy += 0.5;
		
	add_layer("ParamPriorContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-2,{});	
}


/// Adds a screen allowing parameters priors to be edited
function add_param_prior_content(lay)
{
	let y = 0;
	let dy = 2.0;
	
	let x = maximim_label_width() + 4;
		
	let param_cat = get_param_cat(inf_param_not_needed,"only normal");
	
	if(no_param(param_cat)) center_message("No priors need to be set.",lay);
	
	for(let cat = 0; cat < param_cat.length; cat++){
		let pc = param_cat[cat];
		if(pc.list.length > 0){
			y = lay.add_subtitle(pc.name,1,y+0.2,WHITE,{te:pc.inf_te});
	
			for(let cati = 0; cati < pc.list.length; cati++){
				let i = pc.list[cati];
				let par = model.param[i];
			
				let w = wright2;
				
				if(par.dep.length > 0){
					w -= 4.5;
					lay.add_checkbox(w,y+0.4,"Split","Split",par.prior_split_check,WHITE);
				}
				
				if(par.dep.length == 0 || par.prior_split_check.check == false){
					display_distribution(i,x,y,lay,true,true,w);
				}
				else{
					if(par.prior_split_desc == no_elements){
						display_no_element(par,x,y,lay,w);
					}
					else{
						display_distribution_split(i,x,y,lay,true,true,"prior",w);
					}
				}
				
				y += dy;
			}
		}
	}
}


/// Displays a row allowing a variable distribution to be set
function display_distribution(i,x,y,lay,allow_eqn,allow_edit,w)
{		
	let par = model.param[i];
	let prior = par.prior;

	lay.display_param(x-par.label_info.dx-0.7,y-0.1,par.label_info);
	
	let si = 1.5;
	lay.add_button({te:"~", x:x, y:y, dy:si, type:"Text", font:get_font(si), si:si, col:BLACK});

	let te = get_prior_string(prior);
	
	let fo = get_font(1.1,"","times");
	
	let ac; if(allow_edit != false) ac = "EditPriorElement";
	
	let ty = "PriorElement"; if(par.variety == "dist") ty = "DistElement";
	
	lay.add_button({te:te, x:x+1.6, y:y, dx:w-x-1.6, dy:1.6, type:ty, font:fo, ac:ac, i:i, name:par.name, allow_eqn:allow_eqn});
}


/// Converts a prior to a string
function get_prior_string(prior)
{
	switch(prior.type.te){
	case select_str:
		return select_str;
		
	case "fix":
		return "fix("+prior.value.mean_eqn.te+")";
		
	case "uniform":
		return "uniform("+prior.value.min_eqn.te+","+prior.value.max_eqn.te+")";
		
	case "exp":
		return "exp("+prior.value.mean_eqn.te+")";
		
	case "normal":
		return "normal("+prior.value.mean_eqn.te+","+prior.value.sd_eqn.te+")";
		
	case "gamma":
		return "gamma("+prior.value.mean_eqn.te+","+prior.value.cv_eqn.te+")";
		
	case "log-normal":
		return "log-normal("+prior.value.mean_eqn.te+","+prior.value.cv_eqn.te+")";
		
	case "beta":
		return "beta("+prior.value.alpha_eqn.te+","+prior.value.beta_eqn.te+")";
	
	case "bernoulli":
		return "bernoulli("+prior.value.mean_eqn.te+")";
	
	case "flat":
		return "flat";
	
	case "dirichlet":
		return "dir("+prior.value.alpha_eqn.te+")";
		
	default: error("Option not recognised 105"); break;
	}
}


/// Gets a description of value to send back to inferface
function get_prior_split_desc(par)
{
	let te;
	if(par.prior_split != undefined){
		if(par.comb_list.length == 0){
			return no_elements;
		}
		
		let dim = get_dimensions(par.prior_split);
		switch(dim.length){
		case 1: te = "Vector("; break;
		case 2: te = "Matrix("; break;
		default: te = "Tensor("; break;
		}
		
		let tot = 1; 
		for(let k = 0; k < dim.length; k++){
			if(k != 0) te += ",";
			te += dim[k];
			tot *= dim[k];
		}
		te += ")";
		
		if(par.prior_split_set == true && tot < 10){
			let tensor = get_prior_split_tensor(par)
			te = JSON.stringify(tensor);
			te = te.replace(/\"/g," ");
			te = te.replace(/null/g,"Unset");
		}
	}
	
	return te;
}


/// Displays a row allowing a variable distribution to be set
function display_distribution_split(i,x,y,lay,allow_eqn,allow_edit,variety,w)
{
	let par = model.param[i];
	
	lay.display_param(x-par.label_info.dx-0.7,y-0.1,par.label_info);

	let si = 1.5;
	lay.add_button({te:"~", x:x, y:y, dy:si, type:"Text", font:get_font(si), si:si, col:BLACK});
	
	if(par.dep.length == 0){
		error("Should not be here");
		return;
	}
	
	let te="";
	
	if(par.prior_split_set != true) te += "Unset ";
	
	te += par.prior_split_desc;
	
	let fo = get_font(1.1,"","times");
		
	let ac;
	if(allow_edit != false){
		ac = "EditDistSplitValue"; 
		if(variety == "prior") ac = "EditPriorSplitValue";
	}
	
	let ty = "ParamSimElement"; if(par.variety == "dist") ty = "DistSimElement";
	if(te == no_elements) ac = undefined;
	
	lay.add_button({te:te, x:x+1.6, y:y+0., dx:w-x-1.6, dy:1.6, type:ty, font:fo, ac:ac, i:i, name:par.name, label_info:par.label_info, allow_eqn:allow_eqn});
}



/// Creates a text tensor from prior_split
function get_prior_split_tensor(par)
{
	let pri = par.prior_split;

	let tensor = copy(pri);
		
	let dim = get_dimensions(tensor);
	let ele_list = get_element_list(tensor,dim);
	
	for(let k = 0; k < ele_list.length; k++){
		let pindex = ele_list[k];
			
		let te = get_prior_string(get_element(pri,pindex)); 
			
		set_element(tensor,pindex,te);
	}
	
	return tensor;
}


/// Converts a text string to a prior specification
function convert_text_to_prior(te,pri_pos,dist)
{
	let pri = unset_prior();
	
	let spl = te.trim().split("(");
	if(spl.length != 2) return err("Prior has a syntax error");
	
	if(spl[1].length == 0) return err("Prior has a syntax error");
	if(spl[1].substr(spl[1].length-1,1) != ")") return err("Prior has a syntax error");
	
	let type = spl[0].toLowerCase();
	if(type == "dir") type = "dirichlet";
	if(type == "bern") type = "bernoulli";
	
	if(!pri_pos) error("prior pos is not set!");

	if(find_in(pri_pos,type) == undefined){
		if(find_in(prior_pos,type) == undefined){
			return err("'"+type+"' is not a valid prior type");
		}
		else{
			return err("The prior type '"+type+"' cannot be used for this parameter");
		}
	}
	
	let bra = spl[1].substr(0,spl[1].length-1);
	let spl2 = bra.split(",");
		
	pri.type.te = type;
	
	switch(type){
	case "fix":
	{
		if(spl2.length != 1) return err("Expected one value in the brackets");
			if(dist == true){
				if(!is_eqn(spl2[0],"value",{})){
					return err("'"+spl2[0]+"' is not valid");
				}
			}
			else{
				if(isNaN(spl2[0])) return err("'"+spl2[0]+"' is not a number");
				//if(Number(spl2[0]) <= 0) return err("Value must be positive");
			}
			pri.value.mean_eqn.te = spl2[0];
		}
		break;
		
	case "uniform":
		{
			if(spl2.length != 2) return err("Expected two values in the brackets");
			if(dist == true){
				if(!is_eqn(spl2[0],"min",{})) return err("'"+spl2[0]+"' is not valid");
				if(!is_eqn(spl2[1],"max",{})) return err("'"+spl2[1]+"' is not valid");
				if(!isNaN(spl2[0]) && !isNaN(spl2[1]) && Number(spl2[0]) >= Number(spl2[1])){
					return err("Minimum must be smaller than maximum");
				}
			}
			else{
				if(isNaN(spl2[0])) return err("'"+spl2[0]+"' is not a number");
				if(isNaN(spl2[1])) return err("'"+spl2[1]+"' is not a number");
				if(Number(spl2[0]) >= Number(spl2[1])) return err("Minimum must be smaller than maximum");
			}
			pri.value.min_eqn.te = spl2[0];
			pri.value.max_eqn.te = spl2[1];
		}
		break;
	
	case "exp":
		{
			if(spl2.length != 1) return err("Expected one value in the brackets");
			if(dist == true){
				if(!is_eqn(spl2[0],"mean",{positive:true})) return err("'"+spl2[0]+"' is not valid");
			}
			else{
				if(isNaN(spl2[0])) return err("'"+spl2[0]+"' is not a number");
				if(Number(spl2[0]) <= 0) return err("Mean must be positive");
			}
			pri.value.mean_eqn.te = spl2[0];
		}
		break;
		
	case "normal":
		{
			if(spl2.length != 2) return err("Expected two values in the brackets");
			if(dist == true){
				if(!is_eqn(spl2[0],"mean",{})) return err("'"+spl2[0]+"' is not valid");
				if(!is_eqn(spl2[1],"sd",{positive:true})) return err("'"+spl2[1]+"' is not valid");
			}
			else{
				if(isNaN(spl2[0])) return err("'"+spl2[0]+"' is not a number");
				if(isNaN(spl2[1])) return err("'"+spl2[1]+"' is not a number");
				if(Number(spl2[1]) <= 0) return err("Standard deviation must be positive");
			}
			
			pri.value.mean_eqn.te = spl2[0];
			pri.value.sd_eqn.te = spl2[1];
		}
		break;
		
	case "gamma":
		{
			if(spl2.length != 2) return err("Expected two values in the brackets");
			if(dist == true){
				if(!is_eqn(spl2[0],"mean",{positive:true})) return err("'"+spl2[0]+"' is not valid");
				if(!is_eqn(spl2[1],"cv",{positive:true})) return err("'"+spl2[1]+"' is not valid");
			}
			else{
				if(isNaN(spl2[0])) return err("'"+spl2[0]+"' is not a number");
				if(isNaN(spl2[1])) return err("'"+spl2[1]+"' is not a number");
				if(Number(spl2[0]) <= 0) return err("Mean must be positive");
				if(Number(spl2[1]) <= 0) return err("CV must be positive");
			}
			
			pri.value.mean_eqn.te = spl2[0];
			pri.value.cv_eqn.te = spl2[1];
		}
		break;
		
	case "log-normal":
		{
			if(spl2.length != 2) return err("Expected two values in the brackets");
			if(dist == true){
				if(!is_eqn(spl2[0],"mean",{})) return err("'"+spl2[0]+"' is not valid");
				if(!is_eqn(spl2[1],"cv",{positive:true})) return err("'"+spl2[1]+"' is not valid");
			}
			else{
				if(isNaN(spl2[0])) return err("'"+spl2[0]+"' is not a number");
				if(isNaN(spl2[1])) return err("'"+spl2[1]+"' is not a number");
				if(Number(spl2[0]) <= 0) return err("Mean must be positive");
				if(Number(spl2[1]) <= 0) return err("CV must be positive");
			}
			pri.value.mean_eqn.te = spl2[0];
			pri.value.cv_eqn.te = spl2[1];
		}
		break;
		
	case "beta":
		{
			if(spl2.length != 2) return err("Expected two values in the brackets");
			if(dist == true){
				if(!is_eqn(spl2[0],"alpha",{})) return err("'"+spl2[0]+"' is not valid");
				if(!is_eqn(spl2[1],"beta",{positive:true})) return err("'"+spl2[1]+"' is not valid");
			}
			else{
				if(isNaN(spl2[0])) return err("'"+spl2[0]+"' is not a number");
				if(isNaN(spl2[1])) return err("'"+spl2[1]+"' is not a number");
				if(Number(spl2[0]) <= 0) return err("Alpha must be positive");
				if(Number(spl2[1]) <= 0) return err("Beta must be positive");
			}
			
			pri.value.alpha_eqn.te = spl2[0];
			pri.value.beta_eqn.te = spl2[1];
		}
		break;

	case "bernoulli":
		{
			if(spl2.length != 1) return err("Expected one value in the brackets");
			if(dist == true){
				if(!is_eqn(spl2[0],"mean",{zero_one_range:true})){
					return err("'"+spl2[0]+"' is not valid");
				}
			}
			else{
				if(isNaN(spl2[0])) return err("'"+spl2[0]+"' is not a number");
				if(Number(spl2[0]) < 0 || Number(spl2[0]) > 1) return err("Mean must be between 0 and 1");
			}
			
			pri.value.mean_eqn.te = spl2[0];
		}
		break;
		
	case "flat":
		break;
		
	case "dirichlet":
		{
			if(spl2.length != 1) return err("Expected one value in the brackets");
			if(dist == true){
				if(!is_eqn(spl2[0],"alpha",{positive:true})) return err("'"+spl2[0]+"' is not valid");
			}
			else{
				if(isNaN(spl2[0])) return err("'"+spl2[0]+"' is not a number");
				if(Number(spl2[0]) > 0) return err("Alpha must be positive");
			}
			
			pri.value.alpha_eqn.te = spl2[0];
		}
		break;
	
	default: error("Option problem"); return err("Unknown problem");
	}

	return pri;
}


/// This is activated when the 'done' button is pressed on the edit prior bubble
function done_prior()
{
	if(bubble_check_error() == false){
		copy_back_to_source();
		
		if(invalid_prior() == false){
			let bubpri = inter.bubble.prior;
		
			let pri = model.param[bubpri.i].prior;
			pri.type = bubpri.type;
			pri.value = bubpri.value;
			update_param();
			close_bubble();
		}
	}
}


/// This is activated when the 'done' button is pressed on the edit prior bubble
function done_comp_prior()
{
	if(bubble_check_error() == false){
		copy_back_to_source();
		
		if(invalid_prior() == false){
			let bubpri = inter.bubble.prior;
			let bu = inter.bubble.bu;
			let pri = bu.prior;
			pri.type = bubpri.type;
			pri.value = bubpri.value;
			update_param();
			close_bubble();
		}
	}
}


/// Checks to see if the prior is invalid
function invalid_prior()
{
	let bubpri = inter.bubble.prior;
		
	if(bubpri.type.te == "uniform"){
		if(Number(bubpri.value.max_eqn.te) <= (Number(bubpri.value.min_eqn.te))){
			set_warning("Must be larger than minimum value",["prior_max","prior_dist_max"]);
			return true;
		}
	}
	
	return false;
}


/// Sets a warning 
function set_warning(te,type_arr)
{
	let sto = inter.textbox_store;
	
	for(let i = 0; i < sto.length; i++){
		let ty = sto[i].source.type;
		for(let j = 0; j < type_arr.length; j++){
			let ty2 = type_arr[j];
			if(ty.length >= ty2.length){
				if(ty.substr(0,ty2.length) == ty2){
					sto[i].warning = te; return;
				}
			}
		}
	}
	
	error("Could not set warning");
}

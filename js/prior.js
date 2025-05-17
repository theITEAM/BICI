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
				
				if(par.dep.length > 0 && !par.factor){
					w -= 4.5;
					lay.add_checkbox(w,y+0.4,"Split","Split",par.prior_split_check,WHITE);
				}
				
				if(par.dep.length == 0 || par.prior_split_check.check == false || par.factor){
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
	
	case "dirichlet":
		return "dir("+prior.value.alpha_eqn.te+")";
		
	case "mdir":
		return "mdir("+prior.value.sigma_eqn.te+")";
	
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
	if(spl.length != 2) return err("Syntax error");
	
	if(spl[1].length == 0) return err("Syntax error");
	if(spl[1].substr(spl[1].length-1,1) != ")") return err("Syntax syntax error");
	
	let type = spl[0].toLowerCase();
	if(type == "dir") type = "dirichlet";
	if(type == "bern") type = "bernoulli";
	
	if(!pri_pos) error("prior pos is not set!");

	if(find_in(pri_pos,type) == undefined){
		if(find_in(prior_pos,type) == undefined){
			return err("'"+type+"' is not a valid type");
		}
		else{
			return err("The type '"+type+"' cannot be used for this parameter");
		}
	}
	
	let bra = spl[1].substr(0,spl[1].length-1);
	let spl2 = bra.split(",");
		
	pri.type.te = type;
	
	let er;
	
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
				if(isNaN(spl2[0])) return err("'"+spl2[0]+"' must be a number");
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
			
			er = check_thresh(EXP_MEAN_TE,EXP_MEAN_QU,spl2[0]); if(er) return er;
			
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
		
			er = check_thresh(NORM_TE,SD_QU,spl2[1]); if(er) return er;
			er = check_thresh(NORM_TE,NORM_MEAN_QU,spl2[0]); if(er) return er;

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
			
			er = check_thresh(GAMMA_TE,CV_QU,spl2[1]); if(er) return er;
			er = check_thresh(GAMMA_TE,MEAN_QU,spl2[0]); if(er) return er;

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
			
			er = check_thresh(LOGNORM_TE,CV_QU,spl2[1]); if(er) return er;
			er = check_thresh(LOGNORM_TE,MEAN_QU,spl2[0]); if(er) return er;

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
			
			er = check_thresh(BETA_TE,BETA_QU,spl2[1]); if(er) return er;
			er = check_thresh(BETA_TE,ALPHA_QU,spl2[0]); if(er) return er;

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
			
			er = check_thresh(BERN_TE,BERNP_QU,spl2[0]); if(er) return er;
		
			pri.value.mean_eqn.te = spl2[0];
		}
		break;
		
	case "dirichlet":
		{
			if(spl2.length != 1) return err("Expected one value in the brackets");
			if(dist == true){
				if(!is_eqn(spl2[0],"alpha",{positive:true})) return err("'"+spl2[0]+"' is not valid");
			}
			else{
				if(isNaN(spl2[0])) return err("'"+spl2[0]+"' is not a number");
				if(Number(spl2[0]) <= 0) return err("Alpha must be positive");
			}
			
			pri.value.alpha_eqn.te = spl2[0];
		}
		break;
	
	case "mdir":
		if(spl2.length != 1) return err("Expected one value in the brackets");
			if(dist == true){
				return err("mdir should not be applied to a distribution");
			}
			else{
				if(isNaN(spl2[0])) return err("'"+spl2[0]+"' is not a number");
				if(Number(spl2[0]) <= 0) return err("Sigma must be positive");
			}
			
			pri.value.sigma_eqn.te = spl2[0];
		break;
		
	default:  err("Distribution '"+type+"' not recognised");
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


/// Checks if values are within thresholds
function check_thresh(dist,dq,val)
{
	if(isNaN(val)) return;
	val = Number(val);
	
	let min_th, max_th;
	
	switch(dq){
	case MEAN_QU: 
		if(val < MEAN_MIN) min_th = MEAN_MIN;
		if(val > MEAN_MAX) max_th = MEAN_MAX; 
		break;
		
	case NORM_MEAN_QU: 
		if(val < NORM_MEAN_MIN) min_th = NORM_MEAN_MIN;
		if(val > NORM_MEAN_MAX) max_th = NORM_MEAN_MAX; 
		break;
	
	case SD_QU: 
		if(val < SD_MIN) min_th = SD_MIN;
		if(val > SD_MAX) max_th = SD_MAX; 
		break;
	
	case CV_QU: 
		if(val < CV_MIN) min_th = CV_MIN;
		if(val > CV_MAX) max_th = CV_MAX; 
		break;
		
	case SHAPE_QU: 
		if(val < SHAPE_MIN) min_th = SHAPE_MIN;
		if(val > SHAPE_MAX) max_th = SHAPE_MAX; 
		break;
		
	case SCALE_QU: 
		if(val < SCALE_MIN) min_th = SCALE_MIN;
		if(val > SCALE_MAX) max_th = SCALE_MAX; 
		break;
	
	case ALPHA_QU: case BETA_QU:
		if(val < ALPBETA_MIN) min_th = ALPBETA_MIN;
		if(val > ALPBETA_MAX) max_th = ALPBETA_MAX; 
		break;
		
	case P_QU: 
		if(val < P_MIN) min_th = P_MIN;
		if(val > P_MAX) max_th = P_MAX; 
		break;
		
	case BERNP_QU: 
		if(val < 0) min_th = 0;
		if(val > 1) max_th = 1; 
		break;
	
	case RATE_QU: 
		if(val < RATE_MIN) min_th = RATE_MIN;
		break;
		
	case EXP_MEAN_QU: 
		if(val < EXP_MEAN_MIN) min_th = EXP_MEAN_MIN;
		break;
		
	case TIME_QU: 
		if(val < TIME_MIN) min_th = TIME_MIN;
		break;
		
	case POIS_QU: 
		if(val < LAM_MIN) min_th = LAM_MIN;
		if(val > LAM_MAX) max_th = LAM_MAX; 
		break;
	}

	let dist_te;
	if(min_th != undefined || max_th != undefined){
		switch(dist){
		case LOGNORM_TE: dist_te = "Log-normal"; break;
		case NORM_TE: dist_te = "Normal"; break;
		case WEIBULL_TE: dist_te = "Weibull"; break;
		case GAMMA_TE: dist_te = "Gamma"; break;
		case BETA_TE: dist_te = "Beta"; break;
		case NEGBINO_TE: dist_te = "Negitive-binomial"; break;
		case BERN_TE: dist_te = "Bernoulli"; break;
		case EXP_RATE_TE: dist_te = "Exponential"; break;
		case EXP_MEAN_TE: dist_te = "Exponential"; break;
		case POIS_TE: dist_te = "Poisson"; break;
		case PERIOD_TE: dist_te = "Period"; break;
		}
		
		let quant;
		switch(dq){
		case MEAN_QU: quant = "mean"; break;
		case NORM_MEAN_QU: quant = "mean"; break;
		case SD_QU: quant = "standard deviation"; break;
		case CV_QU: quant = "coefficient of variation"; break;
		case SHAPE_QU: quant = "shape"; break;
		case SCALE_QU: quant = "scale"; break;
		case ALPHA_QU: quant = "alpha"; break;
		case BETA_QU: quant = "beta"; break;
		case P_QU: quant = "probability"; break;
		case BERNP_QU: quant = "probability"; break;
		case RATE_QU: quant = "rate"; break;
		case EXP_MEAN_QU: quant = "mean"; break;
		case POIS_QU: quant = "mean"; break;
		case TIME_QU: quant = "time"; break;
		}		
		
		if(min_th != UNSET){
			return err(dist_te+" "+quant+" has value '"+val+"' which should not be below threshold '"+min_th+"'");
		}
		
		if(max_th != UNSET){
			return err(dist_te+" "+quant+" has value '"+val+"' which should not be above threshold '"+max_th+"'");
		}
	}
}


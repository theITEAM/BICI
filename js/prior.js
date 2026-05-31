"use strict";
// Functions related to the prior

/// Adds the button for the page in which parameter priors are input
function add_param_prior_buts(lay)
{
	let cx = corner.x;
	let cy = corner.y;
	
	cy = lay.add_title("Parameter priors",cx,cy,{te:editprior_text2});
	
	cy += 0.5;
		
	add_layer("ParamPriorContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-3.5,{});	
	
	let active = true;
	if(model.param.length == 0) active = false;
	
	let x = 1.2, y = lay.dy-1.6;
	let gap = 3.5;
	
	let w = model.add_object_button(lay,"Const. Prior",x,y,"SetPriorConst",{ back:WHITE, active:active, info:{}, title:"Constant under inference", te:prior_const_text}); 
	x += w+gap;
}


/// Adds a screen allowing parameters priors to be edited
function add_param_prior_content(lay)
{
	let y = 0;
	let dy = 2.0;
	
	let x = maximim_label_width() + 4;
		
	let param_cat = get_param_cat("inf","only normal");
	
	if(no_param(param_cat)) center_message("No priors need to be set.",lay);
	
	for(let cat = 0; cat < param_cat.length; cat++){
		let pc = param_cat[cat];
		if(pc.list.length > 0){
			y = lay.add_subtitle(pc.name,1,y+0.2,WHITE,{te:pc.inf_te});
	
			for(let cati = 0; cati < pc.list.length; cati++){
				let i = pc.list[cati];
				let par = model.param[i];
			
				let w = wright2;
				
				if(par.prior_const_on == true){
					w -= 1.5;
					
					if(par.ndep_cont == 0){
						display_prior_const(i,x,y,lay,w);
					}
					else{
						if(par.prior_const_desc == no_elements){
							display_no_element(par,x,y,lay,w);
						}
						else{
							if(add_prior_const_view_button(par,w-4,y,i,lay,model)) w -= 4.5;
							
							display_prior_const(i,x,y,lay,w);
						}
					}
					
					let del_x = lay.dx - 2.5;
					let del_dx = 1.3;
					lay.add_button({x:del_x, y:y+0.05, dx:del_dx, dy:del_dx, type:"Delete", i:i, ac:"DeleteParamPriorConst"});
				}
				else{
					if(par.ndep_cont > 0 && !par.factor && !is_symmetric(par)){
						w -= 4.5;
						lay.add_checkbox(w,y+0.4,"Split","Split",par.prior_split_check,WHITE);
					}
					
					if(par.ndep_cont == 0 || par.prior_split_check.check == false || par.factor){
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
	
	case "covar-default":
		return "covar-default()";
	
	case "covar-normal-lkj":
		return "covar-normal-lkj("+prior.value.sd_eqn.te+","+prior.value.eta_eqn.te+")";
	
	case "covar-uniform-lkj":
		return "covar-uniform-lkj("+prior.value.min_eqn.te+","+prior.value.max_eqn.te+","+prior.value.eta_eqn.te+")";
	
	case "covar-inv-wishart":
		return "covar-inv-wishart("+prior.value.S_eqn.te+","+prior.value.nu_eqn.te+")";
	
	case "covar-jeffreys":
		return "covar-jeffreys("+prior.value.min_eqn.te+","+prior.value.max_eqn.te+")";
	
	case "covar-uniform":
		return "covar-uniform("+prior.value.min_eqn.te+","+prior.value.max_eqn.te+")";
	
	case "inverse":
		return "inverse("+prior.value.min_eqn.te+","+prior.value.max_eqn.te+")";
	
	case "power":
		return "power("+prior.value.min_eqn.te+","+prior.value.max_eqn.te+","+prior.value.power_eqn.te+")";
	
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
	
	let dim = [];
	for(let d = 0; d < par.list.length; d++) dim.push(par.list[d].length);
	
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
		
	if(par.prior_split != undefined){
		if(tot < 10){
			let tensor = get_prior_split_tensor(par.prior_split)
			te = JSON.stringify(tensor);
			te = te.replace(/\"/g," ");
			te = te.replace(/null/g,"Unset");
		}
	}
	else te = "Unset "+te;
	
	return te;
}


/// Displays a row allowing a variable distribution to be set
function display_distribution_split(i,x,y,lay,allow_eqn,allow_edit,variety,w)
{
	let par = model.param[i];
	
	lay.display_param(x-par.label_info.dx-0.7,y-0.1,par.label_info);

	let si = 1.5;
	lay.add_button({te:"~", x:x, y:y, dy:si, type:"Text", font:get_font(si), si:si, col:BLACK});
	
	if(par.ndep_cont == 0){
		error("Should not be here");
		return;
	}
	
	let te="";
	
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
function get_prior_split_tensor(pri)
{
	let tensor = copy(pri);
	
	let dim = get_dimensions(tensor);
	let ele_list = get_element_list(tensor,dim);
	
	for(let k = 0; k < ele_list.length; k++){
		let pindex = ele_list[k];
		let pri = get_element(tensor,pindex);
		
		let te;
		if(pri == undefined) te = select_str
		else te = get_prior_string(pri); 
		set_element(tensor,pindex,te);
	}
	
	return tensor;
}


/// Determines if prior is correctly specified
function check_prior(pri,dist)
{
	let er;

	switch(pri.type.te){
	case "fix":
		{
			let mean = pri.value.mean_eqn.te;
			if(dist == true){
				er = is_eqn(mean,"value","norm"); if(er.err == true) return er.msg;
			}
			else{
				if(isNaN(mean)) return "'"+mean+"' must be a number";
			}
		}
		break;
		
	case "uniform":
		{
			let min = pri.value.min_eqn.te;
			let max = pri.value.max_eqn.te;
			
			if(dist == true){
				er = is_eqn(min,"min","norm"); if(er.err == true) return er.msg;
				er = is_eqn(max,"max","norm"); if(er.err == true) return er.msg;
				if(!isNaN(min) && !isNaN(max) && Number(min) >= Number(max)){
					return "minimum must be smaller than maximum";
				}
			}
			else{
				if(isNaN(min)) return "'"+min+"' is not a number";
				if(isNaN(max)) return "'"+max+"' is not a number";
				if(Number(min) >= Number(max)) return "minimum must be smaller than maximum";
			}
		}
		break;
	
	case "covar-default":
		break;
		
	case "covar-normal-lkj":
		{
			let sd = pri.value.sd_eqn.te;
			let eta = pri.value.eta_eqn.te;
		
			if(isNaN(sd)) return "'"+sd+"' is not a number";
			if(isNaN(eta)) return "'"+eta+"' is not a number";
			if(Number(sd) < SD_VAR_MIN) return "sd cannot be less than "+SD_VAR_MIN;
			if(Number(sd) > SD_VAR_MAX) return "sd cannot be more than "+SD_VAR_MAX;
			if(Number(eta) < 1) return "eta cannot be less than one";
			if(Number(eta) > ETA_MAX) return "eta cannot be more than "+ETA_MAX;
		}
		break;
		
	case "covar-uniform-lkj":
		{
			let min = pri.value.min_eqn.te;
			let max = pri.value.max_eqn.te;	
			let eta = pri.value.eta_eqn.te;
				
			if(isNaN(min)) return "'"+min+"' is not a number";
			if(isNaN(max)) return "'"+max+"' is not a number";
			if(isNaN(eta)) return "'"+eta+"' is not a number";
			if(Number(min) < VAR_MIN) return "minimum cannot be less than "+VAR_MIN;
			if(Number(max) > VAR_MAX) return "maximum cannot be more than "+VAR_MAX;
			if(Number(min) >= Number(max)) return "minimum must be smaller than maximum";
			if(Number(eta) < 1) return "eta cannot be less than one";
			if(Number(eta) > ETA_MAX) return "eta cannot be more than "+ETA_MAX;
		}
		break;
		
	case "covar-inv-wishart":
		{
			let S = pri.value.S_eqn.te;
			let nu = pri.value.nu_eqn.te;
				
			if(isNaN(S)) return "'"+S+"' is not a number";
			if(isNaN(nu)) return "'"+nu+"' is not a number";
			if(Number(S) <= 0) return "S must be positive";
			
			let num = Number(nu);
			if(num <= 0) return "nu must be positive";
			if(num != Math.floor(num)) "nu must be an integer";
		}
		break;
		
	case "covar-jeffreys":
		{
			let min = pri.value.min_eqn.te;
			let max = pri.value.max_eqn.te;
		
			if(dist == true){
				er = is_eqn(min,"min","norm"); if(er.err == true) return er.msg;
				er = is_eqn(max,"max","norm"); if(er.err == true) return er.msg;
				if(!isNaN(min) && Number(min) <= 0){
					return "minimum must be positive";
				}
				if(!isNaN(max) && Number(max) <= 0){
					return "maximum must be positive";
				}
				if(!isNaN(min) && !isNaN(max) && Number(min) >= Number(max)){
					return "minimum must be smaller than maximum";
				}
			}
			else{
				if(isNaN(min)) return "'"+min+"' is not a number";
				if(isNaN(max)) return "'"+max+"' is not a number";
				if(Number(min) <= 0) return "minimum must be positive";
				if(Number(max) <= 0) return "maximum must be positive";
				if(Number(min) >= Number(max)) return "minimum must be smaller than maximum";
			}
		}
		break;
		
	case "covar-uniform":
		{
			let min = pri.value.min_eqn.te;
			let max = pri.value.max_eqn.te;
			
			if(dist == true){
				er = is_eqn(min,"min","norm"); if(er.err == true) return er.msg;
				er = is_eqn(max,"max","norm"); if(er.err == true) return er.msg;
				if(!isNaN(min) && Number(min) <= 0) return "minimum must be positive";
				if(!isNaN(max) && Number(max) <= 0) return "maximum must be positive";
				if(!isNaN(min) && !isNaN(max) && Number(min) >= Number(max)){
					return "minimum must be smaller than maximum";
				}
			}
			else{
				if(isNaN(min)) return "'"+min+"' is not a number";
				if(isNaN(max)) return "'"+max+"' is not a number";
				if(Number(min) <= 0) return "minimum must be positive";
				if(Number(max) <= 0) return "maximum must be positive";
				if(Number(min) >= Number(max)) return "minimum must be smaller than maximum";
			}
		}
		break;
		
	case "inverse":
		{
			let min = pri.value.min_eqn.te;
			let max = pri.value.max_eqn.te;
			
			if(dist == true){
				er = is_eqn(min,"min","norm"); if(er.err == true) return er.msg;
				er = is_eqn(max,"max","norm"); if(er.err == true) return er.msg;
				if(!isNaN(min) && Number(min) <= 0) return "minimum must be positive";
				if(!isNaN(max) && Number(max) <= 0) return "maximum must be positive";
				if(!isNaN(min) && !isNaN(max) && Number(min) >= Number(max)){
					return "minimum must be smaller than maximum";
				}
			}
			else{
				if(isNaN(min)) return "'"+min+"' is not a number";
				if(isNaN(max)) return "'"+max+"' is not a number";
				if(Number(min) <= 0) return "minimum must be positive";
				if(Number(max) <= 0) return "maximum must be positive";
				if(Number(min) >= Number(max)) return "minimum must be smaller than maximum";
			}
		}
		break;
		
	case "power":
		{
			let min = pri.value.min_eqn.te;
			let max = pri.value.max_eqn.te;
			let pow = pri.value.power_eqn.te;
			
			if(dist == true){
				er = is_eqn(min,"min","norm"); if(er.err == true) return er.msg;
				er = is_eqn(max,"max","norm"); if(er.err == true) return er.msg;
				er = is_eqn(pow,"power","norm"); if(er.err == true) return er.msg;
				if(!isNaN(min) && Number(min) <= 0) return "minimum must be positive";
				if(!isNaN(max) && Number(max) <= 0) return "maximum must be positive";
				if(!isNaN(min) && !isNaN(max) && Number(min) >= Number(max)){
					return "minimum must be smaller than maximum";
				}
			}
			else{
				if(isNaN(min)) return "'"+min+"' is not a number";
				if(isNaN(max)) return "'"+max+"' is not a number";
				if(isNaN(pow)) return "'"+pow+"' is not a number";
				if(Number(min) <= 0) return "minimum must be positive";
				if(Number(max) <= 0) return "maximum must be positive";
				if(Number(min) >= Number(max)) return "minimum must be smaller than maximum";
			}
		}
		break;
	
	case "exp":
		{
			let mean = pri.value.mean_eqn.te;
		
			if(dist == true){
				er = is_eqn(mean,"mean","pos"); if(er.err == true) return er.msg;
			}
			else{
				if(isNaN(mean)) return "'"+mean+"' is not a number";
				if(Number(mean) <= 0) return "mean must be positive";
			}
			
			er = check_thresh(EXP_MEAN_TE,EXP_MEAN_QU,mean); if(er.err == true) return er.msg;
		}
		break;
		
	case "normal":
		{
			let mean = pri.value.mean_eqn.te;
			let sd = pri.value.sd_eqn.te;
			
			if(dist == true){
				er = is_eqn(mean,"mean","norm"); if(er.err == true) return er.msg;
				er = is_eqn(sd,"sd","pos"); if(er.err == true) return er.msg;
			}
			else{
				if(isNaN(mean)) return "'"+mean+"' is not a number";
				if(isNaN(sd)) return "'"+sd+"' is not a number";
				if(Number(sd) <= 0) return "standard deviation must be positive";
			}
			pr("sd");
		pr(sd);
		pr(er);
			er = check_thresh(NORM_TE,SD_QU,sd); if(er.err == true) return er.msg;
			er = check_thresh(NORM_TE,NORM_MEAN_QU,mean); if(er.err == true) return er.msg;
		}
		break;
		
	case "gamma":
		{
			let mean = pri.value.mean_eqn.te;
			let cv = pri.value.cv_eqn.te;
			
			if(dist == true){
				er = is_eqn(mean,"mean","pos"); if(er.err == true) return er.msg;
				er = is_eqn(cv,"cv","pos"); if(er.err == true) return er.msg;
			}
			else{
				if(isNaN(mean)) return "'"+mean+"' is not a number";
				if(isNaN(cv)) return "'"+cv+"' is not a number";
				if(Number(mean) <= 0) return "mean must be positive";
				if(Number(cv) <= 0) return "cv must be positive";
			}
			
			er = check_thresh(GAMMA_TE,CV_QU,cv); if(er.err == true) return er.msg;
			er = check_thresh(GAMMA_TE,MEAN_QU,mean); if(er.err == true) return er.msg;
		}
		break;
		
	case "log-normal":
		{
			let mean = pri.value.mean_eqn.te;
			let cv = pri.value.cv_eqn.te;
		
			if(dist == true){
				er = is_eqn(mean,"mean","norm"); if(er.err == true) return er.msg;
				er = is_eqn(cv,"cv","pos"); if(er.err == true) return er.msg;
			}
			else{
				if(isNaN(mean)) return "'"+mean+"' is not a number";
				if(isNaN(cv)) return "'"+cv+"' is not a number";
				if(Number(mean) <= 0) return "mean must be positive";
				if(Number(cv) <= 0) return "cv must be positive";
			}
			
			er = check_thresh(LOGNORM_TE,CV_QU,cv); if(er.err == true) return er.msg;
			er = check_thresh(LOGNORM_TE,MEAN_QU,mean); if(er.err == true) return er.msg;
		}
		break;
		
	case "beta":
		{
			let al = pri.value.alpha_eqn.te;
			let be = pri.value.beta_eqn.te;
			
			if(dist == true){
				er = is_eqn(al,"alpha","pos"); if(er.err == true) return er.msg;
				er = is_eqn(be,"beta","pos"); if(er.err == true) return er.msg;
			}
			else{
				if(isNaN(al)) return "'"+al+"' is not a number";
				if(isNaN(be)) return "'"+be+"' is not a number";
				if(Number(al) <= 0) return "alpha value must be positive";
				if(Number(bl) <= 0) return "beta must be positive";
			}
			
			er = check_thresh(BETA_TE,BETA_QU,be); if(er.err == true) return er.msg;
			er = check_thresh(BETA_TE,ALPHA_QU,al); if(er.err == true) return er.msg;
		}
		break;

	case "bernoulli":
		{
			let mean = pri.value.mean_eqn.te;
			
			if(dist == true){
				er = is_eqn(mean,"mean","zeroone"); if(er.err == true) return er.msg;
			}
			else{
				if(isNaN(mean)) return "'"+mean+"' is not a number";
				if(Number(mean) < 0 || Number(mean) > 1) return "mean must be between 0 and 1";
			}
		
			er = check_thresh(BERN_TE,BERNP_QU,mean); 
			if(er.err == true) return er.msg;
		}
		break;
		
	case "dirichlet":
		{
			let al = pri.value.alpha_eqn.te;
		
			if(dist == true){
				er = is_eqn(al,"alpha","pos"); if(er.err == true) return er.msg;
			}
			else{
				if(isNaN(al)) return "'"+al+"' is not a number";
				if(Number(al) <= 0) return "alpha must be positive";
			}
		}
		break;
	
	case "mdir":
		{
			let sig = pri.value.sigma_eqn.te;
			
			if(dist == true){
				return "mdir should not be applied to a distribution";
			}
			else{
				if(isNaN(sig)) return "'"+sig+"' is not a number";
				if(Number(sig) <= 0) return "sigma must be positive";
			}
		}
		break;
		
	default: return "type '"+pri.type.te+"' not recognised";
	}
}


/// Converts a text string to a prior specification
function convert_text_to_prior(te,pri_pos,dist)
{
	let pri = unset_prior();

	let start = "In '"+te+"' ";
	
	let spl = te.trim().split("(");
	if(spl.length != 2) return err(start+"there is a syntax error");
	
	if(spl[1].length == 0) return err(start+"there is a syntax error");
	if(spl[1].substr(spl[1].length-1,1) != ")") return err(start+"there is a syntax error");
	
	let type = spl[0].trim().toLowerCase();
	if(type == "dir") type = "dirichlet";
	if(type == "bern") type = "bernoulli";
	
	if(!pri_pos) error("prior pos is not set!");
	
	if(find_in(pri_pos,type) == undefined){
		if(find_in(prior_pos,type) == undefined && find_in(prior_factor_pos,type) == undefined ){
			return err(start+"the distribution '"+type+"' is not recognised");
		}
		else{
			return err(start+"the prior specification '"+type+"' cannot be used for this parameter");
		}
	}
	
	let bra = spl[1].substr(0,spl[1].length-1);
	let spl2 = bra.split(",");
		
	pri.type.te = type;
	
	let er, wa;

	switch(type){
	case "fix":
		{
			wa = ". The fixed prior should have the format 'fix(value)'";
			if(spl2.length != 1) return err(start+"expected one value in the brackets"+wa);	
			pri.value.mean_eqn.te = spl2[0];
		}
		break;
		
	case "uniform":
		{
			wa = ". The uniform prior should have the format 'uniform(min,max)'";
			if(spl2.length != 2) return err(start+"expected two values in the brackets"+wa);
			pri.value.min_eqn.te = spl2[0];
			pri.value.max_eqn.te = spl2[1];
		}
		break;
	
	case "covar-default":
		{
			wa = ". This prior should have the format 'covar-default()'";
			if(!(spl2.length == 1 && spl2[0] == "")) return err(start+"expected no values in the brackets"+wa);
		}
		break;
		
	case "covar-normal-lkj":
		{
			wa = ". This should have the format 'covar-normal-lkj(sd,eta)'";	
			if(spl2.length != 2) return err(start+"expected two values in the brackets"+wa);
			pri.value.sd_eqn.te = spl2[0];
			pri.value.eta_eqn.te = spl2[1];
		}
		break;
		
	case "covar-uniform-lkj":
		{
			wa = ". The multivariate uniform LKJ prior should have the format 'covar-uniform-lkj(min,max,eta)'";
			if(spl2.length != 3) return err(start+"expected three values in the brackets"+wa);
			pri.value.min_eqn.te = spl2[0];
			pri.value.max_eqn.te = spl2[1];
			pri.value.eta_eqn.te = spl2[2];
		}
		break;
		
	case "covar-inv-wishart":
		{
			wa = ". This should have the format 'covar-inv-wishart(sd,eta)'";	
			if(spl2.length != 2) return err(start+"expected two values in the brackets"+wa);
			pri.value.S_eqn.te = spl2[0];
			pri.value.nu_eqn.te = spl2[1];
		}
		break;
		
	case "covar-jeffreys":
		{
			wa = ". The multivariate Jeffreys prior should have the format 'covar-jeffreys(min,max)'";	
			if(spl2.length != 2) return err(start+"expected two values in the brackets"+wa);
			pri.value.min_eqn.te = spl2[0];
			pri.value.max_eqn.te = spl2[1];
		}
		break;
		
	case "covar-uniform":
		{
			wa = ". The multivariate uniform prior should have the format 'covar-uniform(min,max)'";
			if(spl2.length != 2) return err(start+"expected two values in the brackets"+wa);
			pri.value.min_eqn.te = spl2[0];
			pri.value.max_eqn.te = spl2[1];
		}
		break;
		
	case "inverse":
		{
			wa = ". The inverse prior should have the format 'inverse(min,max)'";	
			if(spl2.length != 2) return err(start+"expected two values in the brackets"+wa);
			pri.value.min_eqn.te = spl2[0];
			pri.value.max_eqn.te = spl2[1];
		}
		break;
		
	case "power":
		{
			wa = ". The power prior should have the format 'power(min,max,power)'";	
			if(spl2.length != 3) return err(start+"expected three values in the brackets"+wa);	
			pri.value.min_eqn.te = spl2[0];
			pri.value.max_eqn.te = spl2[1];
			pri.value.power_eqn.te = spl2[2];
		}
		break;
	
	case "exp":
		{
			wa = ". The exponential prior should have the format 'exp(mean)'";	
			if(spl2.length != 1) return err(start+"expected one value in the brackets"+wa);
			pri.value.mean_eqn.te = spl2[0];
		}
		break;
		
	case "normal":
		{
			wa = ". The normal prior should have the format 'normal(mean,sd)'";	
			if(spl2.length != 2) return err("Expected two values in the brackets");
			pri.value.mean_eqn.te = spl2[0];
			pri.value.sd_eqn.te = spl2[1];
		}
		break;
		
	case "gamma":
		{
			wa = ". The gamma prior should have the format 'gamma(mean,cv)'";
			if(spl2.length != 2) return err(start+"expected two values in the brackets"+wa);
			pri.value.mean_eqn.te = spl2[0];
			pri.value.cv_eqn.te = spl2[1];
		}
		break;
		
	case "log-normal":
		{
			wa = ". The log-normal prior should have the format 'log-normal(mean,cv)'";
			if(spl2.length != 2) return err(start+"expected two values in the brackets"+wa);
			pri.value.mean_eqn.te = spl2[0];
			pri.value.cv_eqn.te = spl2[1];
		}
		break;
		
	case "beta":
		{
			wa = ". The beta prior should have the format 'beta(alpha,beta)'";
			if(spl2.length != 2) return err(start+"expected two values in the brackets"+wa);
			pri.value.alpha_eqn.te = spl2[0];
			pri.value.beta_eqn.te = spl2[1];
		}
		break;

	case "bernoulli":
		{
			wa = ". The Bernoulli prior should have the format 'bern(mean)'";
			if(spl2.length != 1) return err(start+"expected one value in the brackets"+wa);
			pri.value.mean_eqn.te = spl2[0];
		}
		break;
		
	case "dirichlet":
		{
			wa = ". The Dirichlet prior should have the format 'dir(alpha)'";	
			if(spl2.length != 1) return err(start+"expected one value in the brackets"+wa);
			pri.value.alpha_eqn.te = spl2[0];
		}
		break;
	
	case "mdir":
		{
			wa = ". The modified Dirichlet prior should have the format 'mdir(sigma)'";
			if(spl2.length != 1) return err(start+"expected one value in the brackets"+wa);
			pri.value.sigma_eqn.te = spl2[0];
		}
		break;
		
	default: err(start+"distribution '"+type+"' not recognised");
	}

	let priw = check_prior(pri,dist);
	if(priw != undefined) return err(start+priw+wa);
 
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

	switch(bubpri.type.te){
	case "uniform": case "inverse": case "power": 
	case "covar-jeffreys": case "covar-uniform-lkj": case "covar-uniform":
		if(Number(bubpri.value.max_eqn.te) <= (Number(bubpri.value.min_eqn.te))){
			set_warning("Must be larger than minimum value",["prior_max","prior_dist_max"]);
			return true;
		}
		break;
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
	if(isNaN(val)) return success();
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
	
	return success();
}


"use strict";
// Function related to editing parameters


/// Allows the user to edit parameter values
function create_edit_param(lay)
{
	lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, type:"Rect", val:WHITE});
		
	let cx = corner.x;
	let cy = corner.y;

	let title = "Edit tensor values";
	let help = edittensor_text;
	let doneac = "EditParamDone";
	
	let	load_title = "Load tensor values", load_te = load_tensor_text;
	let load_ac = "LoadTensor";
	
	switch(inter.edit_param.type){
	case "PriorConst":
		title = "Edit constant prior values";
		help = editpriorconst_text;
		doneac = "EditPriorConstDone";
		load_title = "Load const values";
		load_te = load_prior_const_text;
		break;
		
	case "weight":
		title = "Edit tensor weights";
		help = edittensorweight_text;
		doneac = "EditWeightDone";
		load_title = "Load weight values";
		load_te = load_weight_text;
		break;
	}
	
	let par = model.param[inter.edit_param.i];
	
	let too_big = false;
	if(inter.edit_param.too_big == true) too_big = true;
	
	switch(lay.op.type){
	case "Const":
		title = "Edit constant tensor values";
		load_te = load_const_tensor_text;
		break;
		
	case "Value":
		title = "Edit tensor values";
		break;
		
	case "Reparam":
		title = "Edit reparameterisation";
		help = editreparam_text;
		doneac = "EditReparamDone";
		load_title = "Load reparameterisation", load_te = load_reparam_text; 
		load_ac = "LoadReparam";
		break;
	
	case "Define":
		title = "Edit definition";
		help = editdefine_text;
		doneac = "EditDefineDone";
		load_title = "Load definition", load_te = load_define_text; 
		load_ac = "LoadDefine";
		break;
		
	case "PriorSplit":
		title = "Edit priors";
		help = editpriorsplit_text; 
		doneac = "EditPriorSplitDone";
		load_title = "Load priors", load_te = load_priorsplit_text; 
		load_ac = "LoadPriorSplit";
		break;
		
	case "DistSplit":
		title = "Edit distributions";
		help = editdistsplit_text;
		doneac = "EditDistSplitDone";
		load_title = "Load distributions", load_te = load_distsplit_text; 
		load_ac = "LoadDistSplit";
		break;
	}

	let dim = par.ndep_cont;
	
	title = convert_tensor_text(title,dim);
	help = convert_tensor_text(help,dim);
	load_title = convert_tensor_text(load_title,dim);
	load_te = convert_tensor_text(load_te,dim);

	if(dim > 0){
		switch(lay.op.type){
		case "Value":	load_te += tab_par("Value",par); break;
		case "Reparam": load_te += tab_par("Equation",par); break;
		case "PriorSplit": load_te += tab_par("Prior",par); break;
		case "DistSplit": load_te += tab_par("Distribution",par); break;
		}
	}
	
	title += " for "+par.full_name;

	let fl = false;
	if(par.dist_mat){ title = "Distance matrix"; help = dist_mat_text; fl = true;}
	if(par.iden_mat){ title = "Identity matrix"; help = iden_mat_text; fl = true;}
	if(par.den_vec){
		if(relative_den(par.name)){ title = "Relative density vector"; help = rden_vec_text; fl = true;}
		else{ title = "Density vector"; help = den_vec_text; fl = true;}
	}
	
	if(is_symmetric(par)) title += " (variances on diagonal, otherwise correlations)";

	if(too_big){
		if(fl) title += " (too large to show all)";
		else title += " (too large to edit)";
	}
	
	cy = lay.add_title(title,cx,cy,{te:help});
	
	add_layer("CreateEditParamContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-3.5,{type:lay.op.type});
	
	if(too_big == true || par.dist_mat || par.iden_mat || par.den_vec){
		lay.add_corner_button([["Back","Grey","CancelEditParam"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	}
	else{
		lay.add_corner_button([["Cancel","Grey","CancelEditParam"],["Update","Grey",doneac]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	}
	
	let x = 1.2, y = lay.dy-1.6;
	let gap = 3.5;
	
	if(!par.dist_mat && !par.iden_mat && !par.den_vec){
		let w = model.add_object_button(lay,"Load",x,y,load_ac,{ back:WHITE, active:true, info:{}, title:load_title, te:load_te}); 
		x += w+gap;
	}
}


/// Determines if relative density
function relative_den(name)
{
	if(begin_str(name,rdensity_name)) return true;
	return false;
}


/// Plots vector/matrix to allow for values to be editted			
function add_create_edit_param_buts(lay)
{
	let eparam = inter.edit_param;
	let i = eparam.i;

	let par = eparam.par_st; if(par == undefined) par = model.param[i];
	
	let action = "EditParamElement";
		
	let dx_table_param = 7;
	
	let ele_type = "ParamElement";
	
	switch(lay.op.type){
	case "weight":
		ele_type = "ParamWeightConst";
		break;
		
	case "Const": case "Value": case "PriorConst":
		ele_type = "ParamElementConst";
		if(par.factor) ele_type = "ParamFactorConst";
		break;
		
	case "Reparam":
		ele_type = "ReparamTableElement"; 
		dx_table_param = 12;
		break;
		
 case "Define":
		ele_type = "DefineTableElement"; 
		dx_table_param = 12;
		break;
		
	case "PriorSplit": 
		ele_type = "PriorSplitElement";
		action = "EditPriorSplitElement";
		dx_table_param = 12;
		break;
		
	case "DistSplit":
		ele_type = "DistSplitElement";
		action = "EditDistSplitElement";
		dx_table_param = 12;
		break;
		
	default: error(lay.op.type); error("option not recog"); break;
	}
	
	if(eparam.too_big == true || par.dist_mat || par.iden_mat || par.den_vec){ 
		ele_type = "TooBigElement"; action = undefined;
	}
	
	let value = eparam.value;
	
	let list = eparam.list;
	
	let dep = par.dep;
	let ndep = par.ndep_cont;
	
	let si_mar = 1;
	let fo_mar = get_font(si_mar);
	let fo_head = get_font(si_mar,"bold italic","times");
	let fo_table = get_font(si_table);
	let mar = 0.4;
	let gap = 0.7;
	let cy = 0;
	let mar_col = EDIT_MARGIN_COL;
	
	let sym = is_symmetric(par);
	
	let w_dep = [];
	for(let i = 0; i < ndep; i++){
		let wmax = 0;
		
		for(let j = 0; j < list[i].length; j++){
			let te = list[i][j];
			let w = text_width(te,fo_mar);
			if(w > wmax) wmax = w;
		}
		if(wmax > width_table_max) wmax = width_table_max ;
		
		w_dep[i] = wmax+1;
	}
		
	if(ndep == 2){                             // Case of matrix	 
		let ch = 0;
		let longest;
		for(let j = 0; j < list[0].length; j++){
			for(let i = 0; i < list[1].length; i++){
				let va = value[j][i];
				let val;
				
				if(sym && i > j){
					val = ".";
				}
				else{
					if(va == undefined) val = "Unset";
					else val = String(value[j][i]);
					if(val.length > ch){
						ch = val.length;
						longest = val;
					}
				}
			}
		}
		
		let welemax = text_width(longest,fo_table);
		welemax++;
		
		if(welemax > width_table_max) welemax = width_table_max ;
		
		let dx = w_dep[1];
		if(welemax > dx) dx = welemax;
  	if(dx < 5) dx = 5;
		
		let cx = 2;
		
		lay.add_button({te:dep[1], x:cx, y:cy, dx:w_dep[0]-0.3, dy:dy_table, type:"RightText", si:si_mar, font:fo_head, col:mar_col});
		
		lay.add_button({te:dep[0], x:cx, y:cy+1.1, dx:w_dep[0]-0.3, dy:dy_table, type:"RightText", si:si_mar, font:fo_head, col:mar_col});
		
		cx += w_dep[0]+gap;
			
		for(let i = 0; i < list[1].length; i++){
			lay.add_button({te:list[1][i], x:cx, y:cy, dx:dx, dy:dy_table_param, type:"Text", si:si_mar, font:fo_mar, col:mar_col});
			cx += dx;
		}
		
		cy += dy_table_param+0.4;
		
		let out_dx = dx*list[1].length+2*mar;
		let out_dy = dy_table_param*list[0].length+2*mar;
		
		if(eparam.too_big == true){
			if(eparam.shrunk[1]) out_dx += 2; 
			if(eparam.shrunk[0]) out_dy += dy_table_param;
		}
		
		cx = 2;
		lay.add_button({x:cx+w_dep[0]+gap-mar, y:cy-mar, dx:out_dx, dy:out_dy, type:"Outline", col:BLACK});
		
		for(let j = 0; j < list[0].length; j++){
			let cx = 2;
			lay.add_button({te:list[0][j], x:cx, y:cy, dx:w_dep[0]-0.3, dy:dy_table_param, type:"RightText", si:si_mar, font:fo_mar, col:mar_col});
			cx += w_dep[0]+gap;
			
			for(let i = 0; i < list[1].length; i++){
				let pindex = [j,i];
				let val = get_element(value,pindex);
				if(val == undefined) val = "Unset";
				
				if(sym && i < j){
					lay.add_button({te:"", x:cx, y:cy, dx:dx, dy:dy_table_param, type:"SymmetricParam", font:fo_table});
				}
				else{
					let sym2 = sym; if(i == j) sym2 = false;
					lay.add_button({te:val, x:cx, y:cy, dx:dx, dy:dy_table_param, type:ele_type, font:fo_table, i:i, pindex:pindex, sym:sym2, ac:action});
				}
				cx += dx;
			}
			
			if(eparam.too_big == true && eparam.shrunk[1]){
				lay.add_button({te:"...", x:cx, y:cy, dx:2, dy:dy_table_param, type:"Text", si:si_mar, font:fo_mar, col:mar_col});
			}
			cy += dy_table_param;
		}
		
		if(eparam.too_big == true && eparam.shrunk[0]){
			let cx = 2+w_dep[0]+gap;
			for(let i = 0; i < list[1].length; i++){
				lay.add_button({te:"⋮", x:cx, y:cy, dx:2, dy:dy_table_param, type:"RightText", si:si_mar, font:fo_mar, col:mar_col});
				cx += dx;
			}
		}
	}
	else{	                                           // Case of tensor
		let dim = get_dimensions(value);
		let ele_list = get_element_list(value,dim);
	
		let cx = 2;
		for(let j = 0; j < ndep; j++){
			lay.add_button({te:dep[j], x:cx, y:cy, dx:w_dep[j], dy:dy_table_param, type:"Text", si:si_mar, font:fo_head, col:mar_col});
			cx += w_dep[j];
		} 
		cy += dy_table_param;
			
		let out_dx = dx_table_param+2*mar;
		let out_dy = ele_list.length*dy_table_param+2*mar;
		
		if(eparam.too_big == true) out_dy += dy_table_param;
		
		lay.add_button({x:cx+gap-mar, y:cy-mar, dx:out_dx, dy:out_dy, type:"Outline", col:BLACK});
		
		for(let k = 0; k < ele_list.length; k++){
			let cx = 2;
			for(let j = 0; j < ndep; j++){
				lay.add_button({te:list[j][ele_list[k][j]], x:cx, y:cy, dx:w_dep[j], dy:dy_table_param, type:"Text", si:si_mar, font:fo_mar, col:mar_col});
				cx += w_dep[j];
			}
		
			cx += gap;
			
			let pindex = ele_list[k];
			
			let val = get_element(value,pindex);
			if(val == undefined) val = "Unset";
				
			lay.add_button({te:val, x:cx, y:cy, dx:dx_table_param, dy:dy_table_param, type:ele_type, font:fo_table, i:i, pindex:pindex, ac:action});
			
			cy += dy_table_param;
		}
		
		if(eparam.too_big == true){
			let cx = 2+gap;
			for(let j = 0; j < ndep; j++) cx += w_dep[j];
			lay.add_button({te:"⋮", x:cx+dx_table_param/2-0.3, y:cy, dx:dx_table_param, dy:dy_table_param, type:ele_type, font:fo_table, i:i, ac:action});
		}
	}
	
	lay.add_button({x:0, y:cy, dx:0, dy:0.5, type:"Nothing"});
}


/// Determines if a parameter is a covariance matrix
function is_symmetric(par)
{
	if(is_matrix(par) && begin(par.name,"Ω")) return true;
	return false;
}


/// Loads tensor values
function load_tensor(ep,source)
{
	let par = model.param[ep.i];
	
	let dep = par.dep;
	let ndep = par.ndep_cont;
	
	let list = par.list;
	if(ep.too_big){
		ep.list = copy(list);
		ep.value = par_find_template(list);
	}
	let ep_value = ep.value;
	
	set_zero(ep_value);
	
	let tab = source.table;

	let hash_list = calc_hash_list(list);

	for(let r = 0; r < tab.nrow; r++){
		let ind=[];
		for(let i = 0; i < ndep; i++){
			let te = tab.ele[r][i].trim();
		
			let te2 = te;
			if(dep[i] == "t") te2 = convert_knottime(te,list[i]);
			
			let k = hash_list[i].find(te2);
			if(k == undefined) alertp("The value '"+te+"' is not found (col "+(i+1)+", row "+(r+1)+")");
			else ind[i] = k;
		}
		
		let va = tab.ele[r][ndep];
		
		if(par.factor == true && va == "*") set_element(ep_value,ind,va);
		else{
			if(isNaN(va)) alertp("Problem loading. The value '"+va+"' on line "+(r+1)+" is not a number");
			let val = Number(va);
			set_element(ep_value,ind,val);
		}
	}
	
	if(ep.too_big){ 
		if(ep.type == "weight"){
			par.factor_weight = copy(ep_value);
			reduce_size(ep,par);
			ep.weight_desc = get_weight_desc(par);
			par.weight_desc = get_weight_desc(par);
		}			
		else{
			par.value = copy(ep_value);
			par.set = true;

			reduce_size(ep,par);
		
			ep.value_desc = get_value_desc(par);
			par.value_desc = get_value_desc(par);
		}
	}
	if(ep.type != "weight") ep.set = par.set;
}


/// Sets all values to zero
function set_zero(value)
{
	let dim = get_dimensions(value);
	let ele_list = get_element_list(value,dim);

	for(let k = 0; k < ele_list.length; k++){
		set_element(value,ele_list[k],0);
	}
}


/// Loads reparameterised values
function load_reparam(ep,source)
{
	let par = ep.par_st;
	
	let dep = par.dep;
	let ndep = par.ndep_cont;
	
	let list = par.list;
	
	if(ep.too_big){
		ep.list = copy(list);
		ep.value = par_find_template(list);
	}
	
	let ep_value = ep.value;
	let tab = source.table;
	
	set_zero(ep_value);
	
	let hash_list = calc_hash_list(list);
	
	for(let r = 0; r < tab.nrow; r++){
		let ind=[];
		for(let i = 0; i < ndep; i++){
			let te = tab.ele[r][i].trim();
			
			let te2 = te;
			if(dep[i] == "t") te2 = convert_knottime(te,list[i]);
			
			let k = hash_list[i].find(te2);
			if(k == undefined) alertp("The value '"+te+"' is not found (col "+(i+1)+", row "+(r+1)+")");
			else ind[i] = k;
		}
		
		let val = tab.ele[r][ndep];
		if(!isNaN(val)) val = Number(val);
		set_element(ep_value,ind,val);
	}
	
	if(ep.too_big){ 
		par.value = copy(ep_value);
		reduce_size(ep,par);
		par.set = true;
		get_defrep_param_list(par);
	
		ep.value_desc = get_value_desc(par);
		par.value_desc = get_value_desc(par);
	}
	ep.set = par.set;
}


/// Creates a hash table for the list
function calc_hash_list(list)
{
	let hash_list = [];

	for(let d = 0; d < list.length; d++){
		hash_list[d] = new Hash();
		
		let li = list[d];
		for(let i = 0; i < li.length; i++){
			hash_list[d].add(li[i],i);
		}
	}
	
	return hash_list;
}


/// Loads prior split values
function load_priorsplit(ep,source,dist)
{
	let par = model.param[ep.i];
	let dep = par.dep;
	let ndep = par.ndep_cont;
	
	let list = par.list;
	
	if(ep.too_big){
		ep.list = copy(list);
		ep.value = par_find_template(list);
		ep.prior_split = par_find_template(list);
	}
	
	let prior_split = ep.prior_split;
	let value = ep.value;
	let tab = source.table;
	
	let hash_list = calc_hash_list(list);
	
	for(let r = 0; r < tab.nrow; r++){
		let ind=[];
		for(let i = 0; i < ndep; i++){
			let te = tab.ele[r][i].trim();
			
			let te2 = te;
			if(dep[i] == "t") te2 = convert_knottime(te,list[i]);
			
			let k = hash_list[i].find(te2);
			if(k == undefined) alertp("The value '"+te+"' is not found (col "+(i+1)+", row "+(r+1)+")");
			else ind[i] = k;
		}
	
		let ele = tab.ele[r][ndep];
	
		let pri = convert_text_to_prior(ele,par.pri_pos,dist);
		if(pri.err == true){
			alertp("Problem loading the element '"+ele+"' (col "+(ndep+1)+", row "+(r+1)+"): "+pri.msg+" .");
		}
		
		set_element(prior_split,ind,pri);
		set_element(value,ind,get_prior_string(pri));
	}
	
	if(ep.too_big){ 
		let dim = get_dimensions(prior_split);
		let ele_list = get_element_list(prior_split,dim);

		for(let k = 0; k < ele_list.length; k++){
			let el = get_element(prior_split,ele_list[k]);
			if(el == undefined){
				let te = "Not all elements in the tensor could be loaded (e.g. '";
				for(let i = 0; i < ndep; i++){
					if(i != 0) te += ", ";
					te += "'"+list[i][ele_list[k][i]]+"'";
				}
				te += "). For split priors all elements must be specified.";
				alertp(te);
			}
		}
	
		par.prior_split = copy(prior_split);
		par.prior_split_set = true;
		reduce_size(ep,par);
		ep.prior_split_desc = get_prior_split_desc(par);
		par.prior_split_desc = get_prior_split_desc(par);
	}
	ep.prior_split_set = par.prior_split_set;
}


/// Updates prior_split after it has been editted
function update_prior_split()
{
	let bubpr = inter.bubble.prior;
	let pri = get_element(inter.edit_param.prior_split,bubpr.pindex);

	pri.type = bubpr.type;
	pri.value = bubpr.value;
			
	set_element(inter.edit_param.value,bubpr.pindex,get_prior_string(pri));
}


/// Finds the lists for the dependencies
function par_find_list(par,mod)
{
	let dep = par.dep;
	let ndep = dep.length;
	
	if(ndep == 0) return;
	
	var list = [];
	
	for(let i = 0; i < ndep; i++){	
		if(dep[i] == "t"){
			list[i]	= par.spline.knot;
		}
		else{
			if(remove_prime(dep[i]) == "z"){
				list[i] = find_ieg_list(par,mod);
				if(list.length == 0) error("Cannot find ie list");
			}
			else{
				list[i] = find_comp_from_index(dep[i],mod);
				if(list.length == 0) error("Cannot find compartments from index");
			}
		}
	}
	
	// In the case of the distance matrix truncates if too large
	if(par.dist_mat || par.iden_mat){
		let list_max = Math.floor(Math.sqrt(ELEMENT_MAX))+1;
		
		for(let k = list_max; k < list[0].length; k++){
			list[0][k] = undefined;
			list[1][k] = undefined;
		}		
	}
	
	// In the case of the distance matrix truncates if too large
	if(par.dist_mat || par.iden_mat || par.den_vec){
		let list_max = ELEMENT_MAX;
		for(let k = list_max; k < list[0].length; k++){
			list[0][k] = undefined;
		}		
	}
	
	return list;
}


/// Checks that knot times are set correctly 
function check_knot_times(te)
{		
	let spl = te.split(",");
	for(let j = 0; j < spl.length; j++){
		let te = spl[j].trim();
		if(!(te == "start" && j == 0) && !(te == "end" && j == spl.length-1)){
			if(isNaN(te)){
				return "For '"+spl+"' the value '"+te+"' must be a number";
				break;
			}
			else{
				if(j > 0 && !isNaN(spl[j-1]) && Number(te) <= Number(spl[j-1])){
					return "'"+spl+"' must be time ordered";
				}
			}
		}
	}
	
	return "";
}


/// Changes the knots used in a spline
function change_spline_knot(par,knot)
{
	for(let i = 0; i < knot.length; i++) knot[i] = knot[i].trim();
	
	let dif = false;
	if(par.spline.knot.length != knot.length) dif = true;
	else{
		for(let i = 0; i < knot.length; i++) if(knot[i] != par.spline.knot[i]) dif = true;
	}

	if(dif == true){
		par.spline.knot = knot;
		update_param();
	}
}


/// Checks that parameters are valid
function check_param_valid(type)
{
	for(let i = 0; i < model.param.length; i++){
		let par = model.param[i];
			
		let warn = check_reserved_name(par.name,"par_allow");
		if(warn != ""){
			add_warning({mess:"Parameter name error", mess2:warn, warn_type:"ParamPage", siminf:type, name:par.name});
		}
	
		if(par.variety == "reparam"){
			if(par.time_dep){
				if(par.spline.spline_radio.value != "Square"){
					add_warning({mess:"Reparameterisation error", mess2:"A square spline must be used for time-varying reparameterised parameter "+par.full_name+".", warn_type:"ReparamSquareSpline", name:par.name});
				}
			}
		}
	}
}


/// This ensures that a parameter is in view when it is altered on model->param page (e.g. reparamerisation)
function par_in_view(name)
{
	generate_screen();
	
	let l = find(inter.layer,"name","ModelParamContent");
	if(l == undefined) return;

	let th = find(model.param,"name",name);
	if(th == undefined){ error("Could not find par view"); return;}
	
	let lay = inter.layer[l];
	for(let i = 0; i < lay.but.length; i++){
		let bu = lay.but[i];
		if(bu.i == th){
			switch(bu.type){
			case "ParamSimElement": case "FactorElement": 
			case "ReparamElement":
			case "ReparamEqn": case "DefineEqn":
				shift_button_in_view(l,i);
				return;
			}
		}
	}
	
	prr("Could not find to shift");
	prr(th);
	prr(lay.but);
}

/// Determines if a value needs to be set for a parameter
function sim_value_required(par)
{
	if(par.variety != "const" && par.variety != "reparam" &&  par.variety != "define"){
		if(param_needed(par,"sim")) return true;
	}
	
	return false;
}


/// Replaces parameter values with posterior mean
function sim_param_post_mean()
{
	let param = model.param;
	let inf_param = inf_result.param;
	let details = inf_result.details;
	let burnin = model.inf_res.burnin

	let ps = inf_result.par_sample;
	
	let not_found = [];
	for(let th = 0; th < param.length; th++){
		let par = param[th];
	
		if(sim_value_required(par)){
			let th2 = find(inf_param,"name",par.name);
		
			let fl = false;
			
			if(th2 == undefined) fl = true;
			else{
				let inf_par = inf_param[th2];
				if(!equal_vec(par.dep,inf_par.dep)) fl = true;
				else{
					if(par.ndep_cont == 0){
						let av = 0, nav = 0;
						for(let i = 0; i < ps.length; i++){
							if(ps[i].num >= burnin){
								av += Number(ps[i].param[th2]);
								nav++;
							}
						}							
						if(nav == 0) fl = true;
						else par.value = precision(av/nav,5);
					}
					else{
						for(let d = 0; d < par.dep.length; d++){
							if(!equal_vec(par.list[d],inf_par.list[d])) fl = true;
						}

						if(fl == false){
							let co_list = generate_co_list(par.list);
							
							par.value = par_find_template(par.list);
							for(let k = 0; k < co_list.length; k++){
								let ind = co_list[k].index;
								
								let av = 0, nav = 0;
								for(let i = 0; i < ps.length; i++){
									if(ps[i].num >= burnin){
										let val = Number(get_element(ps[i].param[th2],ind));
										av += val;
										nav++;
									}
								}
							
								if(nav == 0) fl = true;
								else set_element(par.value,ind,precision(av/nav,5));
							}
						}							
					
						par.value_desc = get_value_desc(par);
						par.set = true;
					}
				}				
			}
			
			if(fl) not_found.push(par.full_name);
		}
	}
	
	let st;
	if(not_found.length > 0){
		st = "Parameters not set: ";
		for(let i = 0; i < not_found.length; i++){
			if(i != 0) st += ", ";
			st += not_found[i]; 
		}			
	}
	
	post({ param:strip_heavy(model.param), comment:st});
}


/// On the model parameter page ensures that a section is viewable
function view_section(sec)
{
	let l = find(inter.layer,"name","ModelParamContent");
	if(l == undefined) return;
	
	let but = inter.layer[l].but;
		
	let i = 0; 
	while(i < but.length && !(but[i].te == sec && but[i].type == "CurvedOutline")) i++;
	if(i < but.length) shift_button_in_view(l,i);
	else prr("prob view sec");
}


/// Sets ndep_cont property
function set_ndep_cont(par)
{
	par.ndep_cont = par.dep.length;
	
	if(par.variety == "define" && par.time_dep){
		if(par.dep.length > 0 && par.dep[par.dep.length-1] == "t"){
			par.ndep_cont--;
		}
	}
}

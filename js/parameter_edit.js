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
	
	let too_big = inter.edit_param.too_big;
	
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
		else title += " (too large to edit, please load below)";
		if(title.substr(0,4) == "Edit") title = title.substr(5,1).toUpperCase()+title.substr(6);
	}
	
	cy = lay.add_title(title,cx,cy,{te:help});
	
	add_layer("CreateEditParamContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-3.5,{type:lay.op.type});
	
	//if(too_big == true || par.dist_mat || par.iden_mat || par.den_vec){
	if(par.dist_mat || par.iden_mat || par.den_vec){
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


/// Determines if a parameter is being editted
function editting_param()
{
	let ep = inter.edit_param;
	if(ep.i == undefined) return false;
	
	let par = model.param[inter.edit_param.i];
	if(par.dist_mat || par.iden_mat || par.den_vec) return false;
	return true;
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
	let ep = inter.edit_param;
	
	let i = ep.i;

	let par = ep.par_st; if(par == undefined) par = model.param[i];
	
	let vari = par.variety; if(ep.vari_new != undefined) vari = ep.vari_new;
	
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
	
	let is_cov = is_covar(par);
	
	if(ep.too_big == true || par.dist_mat || par.iden_mat || par.den_vec){ 
		ele_type = "TooBigElement"; action = undefined;
	}
	
	let value = ep.value;
	
	let list = ep.list;
	if(ep.list_shrink != undefined) list = ep.list_shrink;
	
	let dep = par.dep;
	let ndep = get_ndep_cont(par,vari); 
	
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
		
		if(ep.too_big == true){
			if(ep.shrunk[1]) out_dx += 2; 
			if(ep.shrunk[0]) out_dy += dy_table_param;
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
					lay.add_button({te:val, x:cx, y:cy, dx:dx, dy:dy_table_param, type:ele_type, font:fo_table, i:i, pindex:pindex, is_cov:is_cov, sym:sym2, ac:action});
				}
				cx += dx;
			}
			
			if(ep.too_big == true && ep.shrunk[1]){
				lay.add_button({te:"...", x:cx, y:cy, dx:2, dy:dy_table_param, type:"Text", si:si_mar, font:fo_mar, col:mar_col});
			}
			cy += dy_table_param;
		}
		
		if(ep.too_big == true && ep.shrunk[0]){
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
	
		let out_dx = dx_table_param+2*mar;
		let out_dy = ele_list.length*dy_table_param+2*mar;
		
		let too_big_gap=[]
		if(ep.too_big == true){
			out_dy += dy_table_param;
		
			for(let k = 0; k < ele_list.length-1; k++){
				let el1 = ele_list[k];
				let el2 = ele_list[k+1];
				for(let j = 0; j < el1.length; j++){
					if(ep.shrunk[j] && el2[j] < el1[j]){
						too_big_gap[k] = true; 
						out_dy += dy_table_param;
						break;
					}
				}
			}
		}
			
		let cx = 2;
		for(let j = 0; j < ndep; j++){
			lay.add_button({te:dep[j], x:cx, y:cy, dx:w_dep[j], dy:dy_table_param, type:"Text", si:si_mar, font:fo_head, col:mar_col});
			cx += w_dep[j];
		} 
		cy += dy_table_param;
			
		lay.add_button({x:cx+gap-mar, y:cy-mar, dx:out_dx, dy:out_dy, type:"Outline", col:BLACK});
		
		let cx_dots = 2+gap;
		for(let j = 0; j < ndep; j++) cx_dots += w_dep[j];
			
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
			
			if(too_big_gap[k]){
				lay.add_button({te:"⋮", x:cx_dots+dx_table_param/2-0.3, y:cy, dx:dx_table_param, dy:dy_table_param, type:ele_type, font:fo_table, i:i, ac:action});
				cy += dy_table_param;
			}
		}
		
		if(ep.too_big == true){
			lay.add_button({te:"⋮", x:cx_dots+dx_table_param/2-0.3, y:cy, dx:dx_table_param, dy:dy_table_param, type:ele_type, font:fo_table, i:i, ac:action});
		}
	}
	
	lay.add_button({x:0, y:cy, dx:0, dy:0.5, type:"Nothing"});
}


/// Determines if a parameter is a covariance matrix
function is_symmetric(par)
{
	if(is_matrix(par) && is_covar(par)) return true;
	return false;
}


/// Loads tensor values
function load_tensor(ep,source)
{
	let par = model.param[ep.i];
	
	let dim = get_dimensions(ep.value);
	
	let dep = par.dep;
	let ndep = dim.length;
	
	let list = [];
	for(let d = 0; d < ndep; d++) list.push(par.list[d]);
	
	ep.list = list;
	ep.value = par_find_template(list,par.ndep_cont);
	
	let ep_value = ep.value;

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
			
			let el = get_element(ep_value,ind);
			if(el != undefined) ep.multi_set = true;
			set_element(ep_value,ind,val);
		}
	}
	
	ep.set_zero = set_zero(ep_value);
	
	if(ep.too_big){ 
		too_big_value_store = ep.value;
		reduce_size(ep,par);
	}
	
	//if(ep.type != "weight") ep.set = par.set;
	if(ep.type != "weight") ep.set = true;
}


/// Sets any undefined values to zero
function set_zero(value)
{
	let dim = get_dimensions(value);
	let ele_list = get_element_list(value,dim);

	let fl = false;
	for(let k = 0; k < ele_list.length; k++){
		let el = get_element(value,ele_list[k]);
		if(el == undefined){
			set_element(value,ele_list[k],0);
			fl = true;
		}
	}
	
	return fl;
}


/// Loads reparameterised values
function load_reparam(ep,source)
{
	let par = ep.par_st;
	
	let dim = get_dimensions(ep.value);
	
	let dep = par.dep;
	let ndep = dim.length;
	
	let list = [];
	for(let d = 0; d < ndep; d++) list.push(par.list[d]);
	
	ep.list = list;
	ep.value = par_find_template(list,par.ndep_cont);
	
	let ep_value = ep.value;
	let tab = source.table;
	
	let hash_list = calc_hash_list(list);
	
	ep.multi_set = false;
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
		
		let el = get_element(ep_value,ind);
		if(el != undefined) ep.multi_set = true;
		set_element(ep_value,ind,val);
	}
	
	ep.set_zero = set_zero(ep_value);
	
	if(ep.too_big){ 
		too_big_value_store = ep.value;
		reduce_size(ep,par);
	}
	
	ep.set = true;
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
		ep.value = par_find_template(list,par.ndep_cont);
		ep.prior_split = par_find_template(list,par.ndep_cont);
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
		too_big_prior_split_store = ep.prior_split;
		too_big_value_store = ep.value;
		ep.prior_split = undefined;
		reduce_size(ep,par);
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
	
	var list = [];
	
	//if(ndep == 0) return;
	
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
	
	/*
	if(par.too_big && (par.dist_mat || par.iden_mat || par.den_vec)){  // If too big then restricts the size of the list
		for(let i = 0; i < ndep; i++){	
			let list_max = par.list_shrink[i];
			for(let k = list_max; k < list[i].length; k++) list[i][k] = undefined;	
		}
	}
	*/

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
				
				if(par.reparam_eqn_on){				
					let eqn = create_equation(par.reparam_eqn,"reparam_eqn");
					for(let j = 0; j < eqn.param.length; j++){
						if(eqn.param[j].time_dep){
							add_warning({mess:"Reparameterisation error", mess2:"Parameter '"+par.full_name+"' cannot depend on time-varying parameter '"+eqn.param[j].full_name+"'.", warn_type:"RepEqValue", name:par.name});			
						}
					}
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
	if(par.variety != "const" && par.variety != "reparam" && par.variety != "define" && par.variety != "dynamic"){
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
							
							par.value = par_find_template(par.list,par.ndep_cont);
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
					}
					par.set = true;
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


/// Displays any warning after parameter information has been loaded
function load_param_warning(ep)
{
	if(ep.multi_set){
		alert_help("Parameter warning!","Parameter values set multiple times (the last allocation is the one used).");
		ep.multi_set = false;
	}
	else{
		if(ep.set_zero){
			alert_help("Parameter warning!","Unset parameter values set to zero.");
			ep.set_zero = false;
		}
	}
}

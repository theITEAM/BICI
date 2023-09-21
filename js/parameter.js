"use strict";
/// Give all functions relating to parameters

/// Genrates a list of all model parameters
function generate_parameter_list()
{
	//coment("Gather list of parameters");
	let eq_list = model.find_equation_list();
	
	//comment("eq list");comment(eq_list);
	
	let list = [];
	let ind_eff_list = [];
	let fix_eff_list = [];
	
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			for(let c = 0; c < claa.ncomp; c++){
				let co = claa.comp[c];
				if(co.name == ""){
					model.warn.push({mess:"Missing compartment name", mess2:"Name has not been specified in classification '"+claa.name+"'", p:p, cl:cl, c:c, warn_type:"MissingComp"});
				}
				else{
					if(co.col == WHITE){
						model.warn.push({mess:"Compartment colour not set", mess2:"In compartment '"+co.name+"' classification '"+claa.name+"'", p:p, cl:cl, c:c, warn_type:"MissingColour"});
					}
				}
			}
		}
	}
	
	for(let eq = 0; eq < eq_list.length; eq++){
		let eqn = eq_list[eq];

		extract_equation_properties(eqn);

		let info = eqn.eqn_info;
		
		let mess;
		switch(eqn.type){
		case "trans_bp":
			mess = "For the branching probability for transition '"+model.species[info.p].cla[info.cl].tra[info.i].name+"'"
			break;
			
		case "trans_mean": case "trans_rate": case "trans_shape":
		case "trans_scale": case "trans_cv":
			mess = "For transition '"+model.species[info.p].cla[info.cl].tra[info.i].name+"'";
			break;
		
		case "derived":
			{
				let pd = model.derive[info.i];
				mess = "For a derived expression for parameter <e>"+pd.eqn1.te+"</e>"
			}
			break;
			
		case "derive_param":
			{	
				let pd = model.derive[info.i];
				mess = "For derived parameter <e>"+pd.eqn1.te+"</e>"
			}
			break;
			
		case "Se": 
			mess = "For a disease diagnostic test sensitivity"; 
			break;

		case "Sp": 
			mess = "For a disease diagnostic test specificity";
			break;
			
		case "trap_prob":
			mess = "For a trapping probability";
			break;
			
		case "comp_prob":
			mess = "For a compartmental probability";
			break;
						
		case "reparam":
			mess = "For reparameterisation of <e>"+get_full_parameter_name(info.par_name)+"</e>";
			break;
			
		case "prior":
			mess = "For a distribution";
			break;
			
		default: 
			error("type not included: "+eqn.type); 
			break;
		}
			
		if(eqn.warn.length > 0){		
			model.warn.push({mess:mess, mess2:eqn.warn[0].te, eqn_info:eqn.eqn_info, eqn_type:eqn.type, warn_type:"Equation"});
		}
		
		for(let i = 0; i < eqn.param.length; i++){
			let par = copy(eqn.param[i]);
			add_par_to_list(par,list,copy(eqn),mess);
		}
		
		for(let i = 0; i < eqn.fix_eff.length; i++){
			let par = { name:eqn.fix_eff[i].name, age_dep: false, time_dep: false, dep:[], eqn_appear:[], type: "fixed effect"};
			add_par_to_list(par,list,copy(eqn),mess);
		}
		
		for(let i = 0; i < eqn.ind_eff.length; i++){
			let ief = eqn.ind_eff[i]; 
			
			let p_name = eqn.p_name;
			if(ief.p_name_filter != "") p_name = ief.p_name_filter;
			
			let name = eqn.ind_eff[i].name;
			
			let j;
			for(j = 0; j < ind_eff_list.length; j++){
				let ielist = ind_eff_list[j];
				if(ielist.name == name){
					if(ielist.p_name != p_name){
						model.warn.push({mess:mess, mess2:"The individual effect '"+name+"' appears in the transitions of two species", eqn_info:eqn.eqn_info, eqn_type:eqn.type, warn_type:"Equation"});
						break;
					}
				}
			}
			
			if(j == ind_eff_list.length){
				ind_eff_list.push({name:name, p_name:p_name, eqn_appear:[]});
			}
			ind_eff_list[j].eqn_appear.push(copy(eqn));
		}
		
		for(let i = 0; i < eqn.fix_eff.length; i++){
			let fe = eqn.fix_eff[i];
			let name = fe.name;
			
			let p_name = eqn.p_name;
			if(fe.p_name_filter != "") p_name = fe.p_name_filter;
			
			let j;
			for(j = 0; j < fix_eff_list.length; j++){
				let felist = fix_eff_list[j];
				if(felist.name == name){
					if(felist.p_name != p_name){
						model.warn.push({mess:mess, mess2:"The fixed effect '"+name+"' appears in the transitions of two species", eqn_info:eqn.eqn_info, eqn_type:eqn.type, warn_type:"Equation"});
						break;
					}
				}
			}
			
			if(j == fix_eff_list.length){
				fix_eff_list.push({name:name, p_name:p_name, eqn_appear:[]});
			}
			fix_eff_list[j].eqn_appear.push(copy(eqn));
		}
	}
	
	return {par_list:list, ind_eff_list:ind_eff_list, fix_eff_list:fix_eff_list};
}


/// Adds a new parameter to the list
function add_par_to_list(par,list,eqn,mess)
{
	for(let j = 0; j < par.dep.length; j++){
		par.dep[j] = remove_prime(par.dep[j]);
	}
	
	// Adds primes back in for repeated indices
	for(let j = par.dep.length-1; j > 0; j--){
		for(let jj = 0; jj < j; jj++){
			if(par.dep[jj] == par.dep[j]) par.dep[j] += "'";
		}
	}	
	
	let j;
	for(j = 0; j < list.length; j++){
		let parlist = list[j];
		
		if(par.name == parlist.name){
			if(par.type != "derived" && parlist.type != "derived"){					
				if(par.type != parlist.type){
					model.warn.push({mess:"The parameter '"+par.name+"' has contrasting types in different parts of the model", mess2:"This parameter is used "+param_eqn_desc(eqn), eqn_info:eqn.eqn_info, eqn_type:eqn.type, warn_type:"Equation"});
					
					let eqn2 = parlist.eqn_appear[parlist.eqn_appear.length-1];
				
					model.warn.push({mess:"", mess2:"And also "+param_eqn_desc(eqn2), eqn_info:eqn2.eqn_info, eqn_type:eqn2.type, warn_type:"Equation"});
					break;
				}
		
				if(par.p_name != parlist.p_name){
					model.warn.push({mess:mess, mess2:"The parameter '"+par.name+"' appears in the transitions of two species", eqn_info:eqn.eqn_info, eqn_type:eqn.type, warn_type:"Equation"});
					break;
				}
			
				if(par.cl_name != parlist.cl_name){
					model.warn.push({mess:mess, mess2:"The parameter '"+par.name+"' appears in two different classifications", eqn_info:eqn.eqn_info, eqn_type:eqn.type, warn_type:"Equation"});
					break;
				}
			}
			
			if(par.dep.length != parlist.dep.length){	
				model.warn.push({mess:mess, mess2:"The parameter '"+par.name+"' has different dependencies in different parts of the model", warn_type:"Equation", eqn:eqn, eqn_type:eqn.type});
				break;
			}
			break;		
		}
	}
	if(j == list.length){
		par.eqn_appear=[];
		list.push(par);
	}
	
	if(eqn != undefined){
		list[j].eqn_appear.push(eqn);
	}
}


/// Provides a description of the equation that the parameter appears in 
function param_eqn_desc(eqn)
{
	switch(eqn.type){
	case "Se": return "as a diagnostic test sensitivity";
	case "Sp": return "as a diagnostic test specificity";
	case "trap_prob": return "as a trapping probability";
	case "trans_bp": return "as a branching probability";
	case "trans_mean": return "in a transition mean";
	case "trans_rate": return "in a transition rate";
	case "trans_shape": return "in a transition shape";
	case "trans_scale": return "in a transition scale";
	case "trans_cv": return "in a transition coefficient of variation";
	case "comp_prob": return "in a compartmental probability";
	case "reparam": return "in a reparametersation";
	default: error("option not pos"); break;	
	}
}


/// Finds the lists for the dependencies
function par_find_list(par)
{
	let dep = par.dep;
	let ndep = dep.length;
	
	if(ndep == 0) return;
	if(par.type == "ind_eff") return;
	
	var list = [];
	
	for(let i = 0; i < ndep; i++){	
		if(dep[i] == "t" || dep[i] == "a"){
			list[i]	= par.spline.knot;
		}
		else{
			list[i] = find_comp_from_index(dep[i]);
			if(list.length == 0) error("Cannot find compartments from index");
		}
	}
	
	return list;
}


/// Finds the tensor template for a parameter (this is a multi-dimensional tensor with unset elements)
function par_find_template(list)
{
	let ndep = list.length;
	
	let index = [];
	for(let i = 0; i < ndep; i++) index[i] = 0;

	let template = [];
	
	let flag;
	do{
		let ob = template;
		let res = true;
		for(let i = 0; i < ndep; i++){
			let v = index[i];
			
			if(i < ndep-1){ 
				if(ob[v] == undefined) ob[v] = [];
				ob = ob[v];
			}
			else{
				ob[v] = undefined;
			}
		}
	
		let i = 0;
		do{
			flag = false;
			index[i]++; 
			if(index[i] >= list[i].length){
				index[i] = 0; 
				i++; 
				flag = true;
			}
		}while(i < ndep && flag == true);
	}while(flag == false);
	
	return template;
}


/// Gets a list of every possible element (dim_used is optional
function get_element_list(ob,dim,dim_used)
{
	var ele_list =[];
	
	let index = [];	
	let ndim = dim.length;
	if(dim_used != undefined) ndim = dim_used;
	
	for(let i = 0; i < ndim; i++) index[i] = 0;
	
	if(ndim > 0){
		let flag;
		do{
			ele_list.push(copy(index));
			
			let i = 0; 
			do{
				flag = false;
				index[i]++; if(index[i] >= dim[i]){ index[i] = 0; i++; flag = true;}
			}while(i < ndim && flag == true);
		}while(flag == false);
	}
	
	return ele_list;
}


/*
/// Performs an OR opertaion on the needed objects
function add_needed(needed,add)
{
	let dim = get_dimentions(needed);
	
	let ndim = dim.length;
	let index = [];	
	for(let i = 0; i < ndim; i++) index[i] = 0;
	
	let flag;
	do{
		if(get_element(add,index) == true){
			set_element(needed,index,true);
		}
		
		let i = 0; 
		do{
			flag = false;
			index[i]++; if(index[i] >= dim[i]){ index[i] = 0; i++; flag = true;}
		}while(i < ndim && flag == true);
	}while(flag == false);
}
*/


/// Gets the value of the object based on the index
function get_element(ob,ind)
{
	switch(ind.length){
	case 0: return ob;
	case 1: return ob[ind[0]];
	case 2: return ob[ind[0]][ind[1]];
	case 3: return ob[ind[0]][ind[1]][ind[2]];
	case 4: return ob[ind[0]][ind[1]][ind[2]][ind[3]];
	case 5: return ob[ind[0]][ind[1]][ind[2]][ind[3]][ind[4]];
	case 6: return ob[ind[0]][ind[1]][ind[2]][ind[3]][ind[4]][ind[5]];
	default: error("Dimention out of range");
	}
}


/// Sets the value of the object based on the index
function set_element(ob,ind,value)
{
	switch(ind.length){
	case 0: ob = value; break;
	case 1: ob[ind[0]] = value; break;
	case 2: ob[ind[0]][ind[1]] = value; break;
	case 3: ob[ind[0]][ind[1]][ind[2]] = value; break;
	case 4: ob[ind[0]][ind[1]][ind[2]][ind[3]] = value; break;
	case 5: ob[ind[0]][ind[1]][ind[2]][ind[3]][ind[4]] = value; break;
	case 6: ob[ind[0]][ind[1]][ind[2]][ind[3]][ind[4]][ind[5]] = value; break;
	default: error("Dimention out of range");
	}
}


/// Gets the dimentions of the object
function get_dimentions(ob)
{
	let dim=[];
	while(typeof(ob) == "object"){
		if(ob == undefined || ob.length == undefined) break;
		dim.push(ob.length);
		ob = ob[0];
	}
	return dim;
}


/// Given the name of an index this gets a list of the possible compartments
function find_comp_from_index(index)
{
	let list = [];
	index = remove_prime(index);
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			if(claa.index == index){
				for(let c = 0; c < claa.ncomp; c++){
					list.push(claa.comp[c].name);
				}
				return list;
			}				
		}
	}
	return list;
}


/// Finds the classification from the index name
function find_cla_from_index(name)
{
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			if(claa.index == name) return { p:p, cl:cl};
		}
	}
}


/// Finds the index associated with a given compartment
function find_index_from_comp(p_name,name)
{
	let p = find(model.species,"name",p_name);
	if(p == undefined) error("Problem finding index");

	let sp = model.species[p];
	for(let cl = 0; cl < sp.ncla; cl++){
		let claa = sp.cla[cl];
		if(find(claa.comp,"name",name) != undefined){
			return claa.index;
		}
	}
	return "";
}
	

/// Closes the parameter source 
function close_param_source()
{
	inter.edit_param = {};
}

	
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
	
	let par = model.param[inter.edit_param.i];
	
	switch(lay.op.type){
	case "Value":
		title = "Edit constant tensor values";
		break;
		
	case "Reparam":
		title = "Edit parameter reparameterisation";
		help = editreparam_text;
		doneac = "EditReparamDone";
		load_title = "Load reparameterisation", load_te = load_reparam_text; 
		load_ac = "LoadReparam";
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

	let dim = par.dep.length;
	
	title = convert_tensor_text(title,dim);
	help = convert_tensor_text(help,dim);
	load_title = convert_tensor_text(load_title,dim);
	load_te = convert_tensor_text(load_te,dim);
 
	title += " for "+par.full_name;
	cy = lay.add_title(title,cx,cy,{te:help});
	
	//let xx = cx + text_width(te,get_font(1))+0.4;
	
	//lay.display_param(xx,cy-0.15,inter.edit_param.label_info);
	
	//cy = lay.add_paragraph(te,lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
	//cy += 0.5;
	//cy -=1;
	
	add_layer("CreateEditParamContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-3.5,{type:lay.op.type});
	
	lay.add_corner_button([["Cancel","Grey","CancelEditParam"],["Done","Grey",doneac]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	 
	let x = 1.2, y = lay.dy-1.6;
	let gap = 3.5;
	
	let w = model.add_object_button(lay,"Load",x,y,load_ac,{ back:WHITE, active:true, info:{}, title:load_title, te:load_te}); 
	x += w+gap;
}


/// Plots vector/matrix to allow for values to be editted			
function add_create_edit_param_buts(lay)
{
	let i = inter.edit_param.i;
	//source =
	
	let par = model.param[i];

	let action = "EditParamElement";
	let dx_table_param = 7;
	
	let ele_type = "ParamElement";
	
	
	switch(lay.op.type){
	case "Value":
		ele_type = "ParamElementConst";
		break;
		
	case "Reparam":
		ele_type = "ReparamTableElement"; 
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
		
	default: error("option not recog"); break;
	}
	
	let value = inter.edit_param.value;
	//if(
	
	//let need = inter.edit_param.need;
	
	let list = par.list;
	
	let dep = par.dep;
	
	let si_mar = 1;
	let fo_mar = get_font(si_mar);
	let fo_head = get_font(si_mar,"bold italic","times");
	let fo_table = get_font(si_table);
	let mar = 0.4;
	let gap = 0.7;
	let cy = 0;//1;
	let mar_col = DRED;
	
	let w_dep = [];
	for(let i = 0; i < dep.length; i++){
		let wmax = 0;
		
		for(let j = 0; j < list[i].length; j++){
			let te = list[i][j];
			let w = text_width(te,fo_mar);
			if(w > wmax) wmax = w;
		}
		w_dep[i] = wmax+1;
	}
		
	if(dep.length == 2){ // matrix
		let welemax = 0;
		for(let j = 0; j < list[0].length; j++){
			for(let i = 0; i < list[1].length; i++){
				let w = text_width(value[j][i],fo_table);
				if(w > welemax) welemax  = w;
			}
		}
		welemax++;
		
		let dx = w_dep[1];
		if(welemax > dx) dx = welemax;
  	
		let cx = 2;
		
		lay.add_button({te:dep[1], x:cx, y:cy, dx:w_dep[0]-0.3, dy:dy_table, type:"RightText", si:si_mar, font:fo_head, col:mar_col});
		
		lay.add_button({te:dep[0], x:cx, y:cy+1.1, dx:w_dep[0]-0.3, dy:dy_table, type:"RightText", si:si_mar, font:fo_head, col:mar_col});
		
		cx += w_dep[0]+gap;
			
		for(let i = 0; i < list[1].length; i++){
			lay.add_button({te:list[1][i], x:cx, y:cy, dx:dx, dy:dy_table_param, type:"Text", si:si_mar, font:fo_mar, col:mar_col});
			cx += dx;
		}
		cy += dy_table_param+0.4;
		
		cx = 2;
		lay.add_button({x:cx+w_dep[0]+gap-mar, y:cy-mar, dx:dx*list[1].length+2*mar, dy:dy_table_param*list[0].length+2*mar, type:"Outline", col:BLACK});
		
		for(let j = 0; j < list[0].length; j++){
			let cx = 2;
			lay.add_button({te:list[0][j], x:cx, y:cy, dx:w_dep[0]-0.3, dy:dy_table_param, type:"RightText", si:si_mar, font:fo_mar, col:mar_col});
			cx += w_dep[0]+gap;
			
			for(let i = 0; i < list[1].length; i++){
				let pindex = [j,i];
				let val = get_element(value,pindex);
		
				lay.add_button({te:val, x:cx, y:cy, dx:dx, dy:dy_table_param, type:ele_type, font:fo_table, i:i, pindex:pindex, ac:action});
				cx += dx;
			}
			cy += dy_table_param;
		}
	}
	else{	
		let dim = get_dimentions(value);
		let ele_list = get_element_list(value,dim);
	
		let cx = 2;
		for(let j = 0; j < dep.length; j++){
			lay.add_button({te:dep[j], x:cx, y:cy, dx:w_dep[j], dy:dy_table_param, type:"Text", si:si_mar, font:fo_head, col:mar_col});
			cx += w_dep[j];
		} 
		cy += dy_table_param;
			
		lay.add_button({x:cx+gap-mar, y:cy-mar, dx:dx_table_param+2*mar, dy:ele_list.length*dy_table_param+2*mar, type:"Outline", col:BLACK});
		
		for(let k = 0; k < ele_list.length; k++){
			let cx = 2;
			for(let j = 0; j < dep.length; j++){
				lay.add_button({te:list[j][ele_list[k][j]], x:cx, y:cy, dx:w_dep[j], dy:dy_table_param, type:"Text", si:si_mar, font:fo_mar, col:mar_col});
				cx += w_dep[j];
			}
			
			cx += gap;
			
			let pindex = ele_list[k];
			
			let val = get_element(value,pindex);
				
			//let ac; if(get_element(need,pindex) == true) ac = action;
				
			lay.add_button({te:val, x:cx, y:cy, dx:dx_table_param, dy:dy_table_param, type:ele_type, font:fo_table, i:i, pindex:pindex, ac:action});
			
			cy += dy_table_param;
		}
	}
	
	lay.add_button({x:0, y:cy, dx:0, dy:0.5, type:"Nothing"});	
}


/// Loads tensor values
function load_tensor()
{
	let ep = inter.edit_param;
	let par = model.param[ep.i];
	let dep = par.dep;
	let list = par.list;
	let value = ep.value;
	//let need = ep.need;
	let tab = inter.edit_source.table;
	
	for(let r = 0; r < tab.nrow; r++){
		let ind=[];
		let flag = false;
		for(let i = 0; i < dep.length; i++){
			let te = tab.ele[r][i].trim();
			let k = 0; while(k < list[i].length && te != list[i][k]) k++;
			if(k < list[i].length) ind[i] = k;
			else{ flag = true; error("Problem loading");}
		}
		
		if(flag == false){
			let val = Number(tab.ele[r][dep.length]);
			//if(get_element(need,ind) == true){
			set_element(value,ind,val);
			//}
		}
	}
	
	close_data_source();
}


/// Loads reparameterised values
function load_reparam()
{
	let ep = inter.edit_param;
	let par = model.param[ep.i];
	let dep = par.dep;
	let list = par.list;
	let value = ep.value;
	//let need = ep.need;
	let tab = inter.edit_source.table;
	
	for(let r = 0; r < tab.nrow; r++){
		let ind=[];
		let flag = false;
		for(let i = 0; i < dep.length; i++){
			let te = tab.ele[r][i].trim();
			let k = 0; while(k < list[i].length && te != list[i][k]) k++;
			if(k < list[i].length) ind[i] = k;
			else{ flag = true; error("Problem loading");}
		}
		
		if(flag == false){
			let val = tab.ele[r][dep.length];
			if(!isNaN(val)) val = Number(val);
			
			//if(get_element(need,ind) == true){
			set_element(value,ind,val);
			//}
		}
	}
	
	close_data_source();
}


/// Loads prior split values
function load_priorsplit()
{
	let ep = inter.edit_param;
	let par = model.param[ep.i];
	
	let dep = par.dep;
	let list = par.list;
	let prior_split = ep.prior_split;
	let value = ep.value;
	//let need = ep.need;
	let tab = inter.edit_source.table;
	
	for(let r = 0; r < tab.nrow; r++){
		let ind=[];
		let flag = false;
		for(let i = 0; i < dep.length; i++){
			let te = tab.ele[r][i].trim();
			let k = 0; while(k < list[i].length && te != list[i][k]) k++;
			if(k < list[i].length) ind[i] = k;
			else{ flag = true; error("Problem loading");}
		}
		
		if(flag == false){
			//if(get_element(need,ind) == true){
			let ele = tab.ele[r][dep.length];
		
			let pri = convert_text_to_prior(ele,par.pri_pos);
			if(typeof pri == 'string'){
				alertp("Problem loading the element '"+ele+"' (col "+(dep.length+1)+", row "+(r+1)+"): "+pri+".");
				return;
			}
			set_element(prior_split,ind,pri);
			set_element(value,ind,get_prior_string(pri));
			//}
		}
	}
	
	close_data_source();
}


/// Generates a list of all model parameters
/// goes through the exisiting list and copies any info on values/priors
function init_param()
{
	model.warn = [];
	
	let lists = generate_parameter_list();   // Generates a list of all model parameters

	let ielist = lists.ind_eff_list;       
  let felist = lists.fix_eff_list;       
	let par_list = lists.par_list;
	
	update_model_ind_eff(ielist);            // Updates individual effect groupings in model
	
	if(model.warn > 0) return;	
	
	update_model_fix_eff(felist);            // Updates fixed effects in model 
	
	if(model.warn > 0) return;	

	add_ind_eff_param(par_list);             // Adds individual effect parameters (i.e. variances) to model
 	
	update_model_param(par_list);            // Updates the model parameters
	
	create_spline();                         // Creates splines based on parameters

	if(model.warn > 0) return;	

	for(let th = 0; th < model.param.length; th++){  // Sets the distance matrix
		let par = model.param[th];	
		if(par.dist_matrix.check == true){
			set_dist(par);                            
		}
	}
}


// Updates individual effect groups stored in model
function update_model_ind_eff(ielist)
{
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		
		for(let i = 0; i < sp.ind_eff_group.length; i++){
			let ieg = sp.ind_eff_group[i];
			for(let j = 0; j < ieg.list.length; j++) ieg.list[j].used = false;
		}
	
		for(let k = 0; k < ielist.length; k++){
			let ie = ielist[k];
			let flag = false;
			
			if(ie.p_name == sp.name){
				for(let i = 0; i < sp.ind_eff_group.length; i++){
					let ieg = sp.ind_eff_group[i];
					let list = ieg.list;
					for(let j = 0; j < list.length; j++){
						if(ie.name == list[j].name){
							list[j].used = true;
							flag = true;
						}
					}
				}
			
				if(flag == false){
					sp.ind_eff_group.push({p_name:ie.p_name, A_matrix:{check:false, loaded:false, value:[], ind_list:[]}, list:[{name:ie.name}]});
				}
			}
		}
		
		// Deletes unused individual effects
		let i = 0;
		while(i < sp.ind_eff_group.length){
			let ieg = sp.ind_eff_group[i];
			let list = ieg.list;
				
			let j = 0;
			while(j < list.length){
				if(list[j].used == false) list.splice(j,1);
				else j++;
			}
			
			if(list.length == 0) sp.ind_eff_group.splice(i,1);
			else i++;
		}
	}
	
	// Checks that there are no individual effects in population-based models
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
	
		if(sp.type == "Population" && sp.ind_eff_group.length > 0){
			let name = sp.ind_eff_group[0].list[0].name;
			model.warn.push({mess:"Cannot have individual effects", mess2:"A population-based model cannot have individual effects such as '"+name+"'", warn_type:"Equation"});
		}
	}
}


// Updates fixed effects stored in model
function update_model_fix_eff(felist)
{
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		
		for(let i = 0; i < sp.fix_eff.length; i++){
			sp.fix_eff[i].used = false;
		}
	
		for(let k = 0; k < felist.length; k++){
			let fe = felist[k];
		
			let flag = false;
			if(fe.p_name == sp.name){
				for(let i = 0; i < sp.fix_eff.length; i++){
					let feg = sp.fix_eff[i];
					if(feg.name == fe.name){
						feg.used = true;
						flag = true;
					}
				}
				
				if(flag == false){
					sp.fix_eff.push({name:fe.name, X_vector:{loaded:false, value:[], ind_list:[]}});
				}
			}
		}
		
		// Deletes unused fixed effects
		let i = 0;
		while(i < sp.fix_eff.length){
			if(sp.fix_eff[i].used == false) sp.fix_eff.splice(i,1);
			else i++;
		}
	}
	
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
	
		if(sp.type == "Population" && sp.fix_eff.length > 0){
			let name = sp.eqn.fix_eff[0].name;
			model.warn.push({mess:"Cannot have individual effects", mess2:"A population-based model cannot have fixed effects such as '"+name+"'", warn_type:"Equation"});
		}
	}
}


/// Adds individual effect parameters (i.e. variances) to model
function add_ind_eff_param(par_list)
{
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		
		for(let i = 0; i < sp.ind_eff_group.length; i++){
			let ieg = sp.ind_eff_group[i];
			for(let j = 0; j < ieg.list.length; j++){
				for(let jj = j; jj < ieg.list.length; jj++){
					let ie1 = ieg.list[j].name, ie2 = ieg.list[jj].name;
					let sup = "^"+ie1+ie2;
					
					let par;
					if(j == jj){
						par = { name:"Ω"+sup, age_dep: false, time_dep: false, dep:[], eqn_appear:[], type: "variance"};
					}
					else{
						par = { name:"ω"+sup, age_dep: false, time_dep: false, dep:[], eqn_appear:[], type: "correlation"};
					}
				
					add_par_to_list(par,par_list,{p:p, ieg:i, ie1:ie1, ie2:ie2},"Individual effect variance");
				}
			}		
		}
	}
} 	


/// Updates the model parameters
function update_model_param(par_list)
{
	let param_old = model.param;
				
	let param=[];
	for(let i = 0; i < par_list.length; i++){
		let par = par_list[i];
		
		let j;                                 // Determines if parameter already exists
		for(j = 0; j < param_old.length; j++){
			let par_old = param_old[j];
			if(par.name == par_old.name){	
				if(par.dep.length == par_old.dep.length){
					let k = 0; 
					while(k < par.dep.length && par.dep[k] == par_old.dep[k]) k++;
					if(k == par.dep.length) break;
				}
			}
		}
		
		if(j < param_old.length){              // Copy across information from existing parameter definition
			param.push(copy_param_info(par,param_old[j]));
			add_dependent_param(param,par_list);
		}
		else{                                  // Create new parameter
			param.push(create_new_param(par,"normal"));
		}
	}

	param.sort( function(a, b){              // Alp[habetically sorts parameters
		if(a.dep.length > b.dep.length) return 1; 
		if(a.dep.length == b.dep.length){
			if(a.name.toLowerCase() > b.name.toLowerCase()) return 1;
			else return -1;	
		}
		return -1;});
	
	let branch_list = set_model_branch_list(param);
	
	check_derived_param_all(param);
	
	if(model.warn > 0) return;
	
	model.param = param;

	model.branch_list = branch_list;

	branch_auto_calculate();
}


/// Goes through reparameterisation and distribution and adds to parameter list
function add_dependent_param(param,par_list)
{
	let par = param[param.length-1];
	
	let vari = par.variety;
	
	switch(vari){
	case "reparam":
		{
			if(par.dep.length == 0){
				add_ele_param(par.value,vari,par_list);
			}
			else{
				for(let i = 0; i < par.comb_list.length; i++){
					let comb = par.comb_list[i];
					let ele = get_element(par.value,comb.index);
					add_ele_param(ele,vari,par_list);
				}
			}
		}
		break;
		
	case "dist":
		{
			if(par.prior_split_check.check == false){	
				add_distribution_eqn(par.prior,vari,par_list);	
			}
			else{
				for(let i = 0; i < par.comb_list.length; i++){
					let comb = par.comb_list[i];
					let ele = get_element(par.prior_split,comb.index);
				
					add_distribution_eqn(ele,vari,par_list);	
				}
			}
		}
		break;
	}
	
	par.full_name = param_name(par);
}


/// Adds any parameters from a distribution to the list of equations
function add_distribution_eqn(dist,vari,par_list)
{
	let val = dist.value;
	
	switch(dist.type.te){
	case "fix": 
		add_ele_param(val.mean_eqn.te,vari,par_list);
		break;
		
	case "uniform":
		add_ele_param(val.min_eqn.te,vari,par_list);
		add_ele_param(val.max_eqn.te,vari,par_list);
		break;
		
	case "exp": 
		add_ele_param(val.mean_eqn.te,vari,par_list);
		break;
		
	case "normal": 
		add_ele_param(val.mean_eqn.te,vari,par_list);
		add_ele_param(val.sd_eqn.te,vari,par_list);
		break;
	
	case "log-normal": 
		add_ele_param(val.mean_eqn.te,vari,par_list);
		add_ele_param(val.cv_eqn.te,vari,par_list);
		break;
		
	case "gamma":
		add_ele_param(val.mean_eqn.te,vari,par_list);
		add_ele_param(val.cv_eqn.te,vari,par_list);
		break;
		
	case "bernoulli":
		add_ele_param(val.mean_eqn.te,vari,par_list);
		break;
					
	case "beta":
		add_ele_param(val.alpha_eqn.te,vari,par_list);
		add_ele_param(val.beta_eqn.te,vari,par_list);
		break;
	
	case "flat":
		break;
		
	case "dirichlet":
		add_ele_param(val.alpha_eqn.te,vari,par_list);
		break;
	
	default: error("not op:"+dist.type.te); break;
	}
}


/// Adds any parameters from a text element to the parameter list
function add_ele_param(ele,variety,par_list)
{
	if(isNaN(ele)){
		let eqn = create_equation(ele,variety);
	
		for(let i = 0; i < eqn.param.length; i++){
			add_par_to_list(eqn.param[i],par_list,{},"Reparameterisation");
		}
	}
}


/// Creates a new parameter
// par must have the properties: name, dep, type, time_dep, age_dep, (optionally eqn_appear)
// On top of these new properties are added: 
// list: shows the potential values for the dependencies
// comb_list: this is a list giving every possible combination of dependency  
function create_new_param(par,variety)
{
	if(par.name == undefined) error("Parameter name must be set");
	if(par.dep == undefined) error("Parameter dependency must be set");
	if(par.type == undefined) error("Parameter type must be set");
	if(par.time_dep == undefined) error("Parameter time dependency must be set");
	if(par.age_dep == undefined) error("Parameter age dependency must be set");

	par.variety = variety;
	
	par.pri_pos = set_pri_pos(par.type);

	par.label_info = get_label_info(par);

	let knot = ["start","end"];           // Intialises a spline
	
	let on = false;
	if(par.time_dep == true || par.age_dep == true) on = true;
	
	par.spline = { on:on, smooth:default_smooth(), time_dep:par.time_dep, age_dep:par.age_dep, knot:knot}; 
		
	par.sim_sample = {check:true};
		 
	par.prior = unset_prior(par.type);

	par_set_default(par);

	par.dist_matrix = false; 
	
	par.full_name = param_name(par);

	return par;
}


/// The default value for smoothing
function default_smooth()
{
	return {check:false, type:{value:"Log-Normal"}, value:0.1};
}


/// Sets default values for value and prior_split
function par_set_default(par)
{
	if(par.dep.length == 0){
		par.value = set_str;
		if(par.list){
			delete par.list; delete par.comb_list;
			delete par.prior_split; delete par.prior_split_set; delete par.prior_split_check;
			delete par.set;
		}
	}
	else{
		let list = par_find_list(par);
		let temp = par_find_template(list);
		let comb_list = generate_comb_list(list);
		
		let value = copy(temp);
		let prior_split = copy(temp);
		
		let dim = get_dimentions(temp);
		let ele_list = get_element_list(temp,dim);
		
		for(let k = 0; k < ele_list.length; k++){	
			let	val = 0;
			let val_eq = create_equation("","reparam");
			let pr_sp = unset_prior(par.type);
			
			set_element(value,ele_list[k],val);
			set_element(prior_split,ele_list[k],pr_sp);
		}
		
		par.list = list;
		par.comb_list = comb_list;
		par.value = value;
		par.prior_split = prior_split;
		par.prior_split_check = {check:false};
		par.set = false;
	}
}


/// Copies infomation from old parameter definition to new 		
function copy_param_info(par,old)
{
	par.variety = old.variety;
	par.pri_pos = set_pri_pos(par.type);
	par.spline = old.spline;
	par.prior = old.prior;
	par.sim_sample = old.sim_sample;
	par.auto_value = old.auto_value;
	par.label_info = old.label_info;
	par.set = old.set;
	par.dist_matrix = old.dist_matrix;
	
	if(par.dep.length == 0){	
		par.value = old.value;
	}
	else{
		let list = par_find_list(par);
		let temp = par_find_template(list);
		let comb_list = generate_comb_list(list);
		
		let map=[];
		for(let k = 0; k < par.dep.length; k++){
			map[k]=[];
			
			let list_new = list[k];
			let list_old = old.list[k];
			let used=[];
			for(let ij = 0; ij < list_new.length; ij++){
				let val = list_new[ij];
				let m = 0; while(m < list_old.length && !(list_old[m] == val && used[m] == undefined)) m++;
				if(m < list_old.length){ map[k][ij] = m; used[m] = true;}
			}
		}
	
		let value = copy(temp);
		let prior_split = copy(temp);
		
		let dim = get_dimentions(temp);
		let ele_list = get_element_list(temp,dim);
		let prior = copy(old.prior);
		
		for(let k = 0; k < ele_list.length; k++){
			let ind = ele_list[k];
			
			let val, val_eqn, pr_sp;
			
			let ind_old=[];
			let ii;
			for(ii = 0; ii < par.dep.length; ii++){
				ind_old[ii] = map[ii][ele_list[k][ii]];
				if(ind_old[ii] == undefined) break;
			}
			
			if(ii == par.dep.length){
				val = get_element(old.value,ind_old);
				pr_sp = get_element(old.prior_split,ind_old);
			}
			else{
				val = 0;
				val_eqn = create_equation("","reparam");
				pr_sp = unset_prior(par.type);
			}
			
			set_element(value,ind,val);
			set_element(prior_split,ele_list[k],pr_sp);
		}

		par.list = list;
		par.comb_list = comb_list;
		par.value = value;
		par.prior_split = prior_split;
		par.prior_split_set = old.prior_split_set;
		par.prior_split_check = old.prior_split_check;
	}

	return par;
}


function get_label_info(par)
{	
	let fo_big = get_font(si_big,undefined,"Times");
	let fo_sup = get_font(si_sub,undefined,"Times");
	let fo_sub = get_font(si_sub,"italic","Times");

	let te = par.name;
	let spl = te.split("^");
	par.name_raw = spl[0];
	par.sup = ""; if(spl.length == 2)	par.sup = remove_bracket(spl[1]); 
	
	par.name_raw_w = text_width(par.name_raw,fo_big);
	par.sup_w = text_width(par.sup,fo_sup);
	
	let sub = "";
	for(let j = 0; j < par.dep.length; j++){
		if(par.dep[j] != "t" && par.dep[j] != "a"){
			if(sub != "") sub += ",";
			sub += par.dep[j];
		}
	}
	par.sub = sub;
	par.sub_w = text_width(par.sub,fo_sub);
	
	let w = par.name_raw_w;
	let w2 = par.sup_w; if(par.sub_w > w2) w2 = par.sub_w;
	let w_tot = w+w2+0.1;
	
	let age_time = "";
	if(par.time_dep == true) age_time = "(t)";
	if(par.age_dep == true) age_time = "(a)";
	
	let age_time_pos = w_tot;
	if(age_time != ""){
		w_tot += text_width(age_time,fo_big);
	}
	
	return {name_full:par.name, name:par.name_raw, age_time_pos:age_time_pos, age_time:age_time, sup:par.sup, sub:par.sub, fo_big:fo_big, fo_sup:fo_sup, fo_sub:fo_sub, dx:w_tot};
}


/// Finds the possible set of priors for a parameter type
function set_pri_pos(type)
{
	let pri_pos = prior_pos;           
	switch(type){
	case "trans_bp": pri_pos = bp_prior_pos; break;
	case "trap_prob": case "Se": case "Sp": pri_pos = zeroone_prior_pos; break;
	}
	
	return pri_pos;
}


/// Given a list of dependencies this provides a list of cominations going through each
function generate_comb_list(list)
{
	let ndep = list.length;
	let comb_list = [];
	if(ndep == 0) return [{combination:[], index:[]}];
	
	let index=[];
	for(let d = 0; d < ndep; d++) index[d] = 0;

	let flag;
	do{
		let comb = [];
		for(let d = 0; d < ndep; d++) comb.push(list[d][index[d]]);
		
		comb_list.push({combination:comb, index:copy(index)});
		
		let i = 0;
		do{
			flag = false;
			
			index[i]++; 
			if(index[i] >= list[i].length){
				index[i] = 0; 
				i++; 
				flag = true;
			}
		}while(i < ndep && flag == true);
	}while(flag == false);
	
	return comb_list;
}

				
/// Sets up a description of the branching in the model
function set_model_branch_list(param)
{
	// Ensures that one branching prob is set to auto_set
	let branch_list = model.find_branch_list();

	let gg = copy(branch_list);
	
	for(let i = 0; i < branch_list.length; i++){
		let bl = branch_list[i];
		
		let par_list = [];
		for(let j = 0; j < bl.tra.length; j++){
			let fl = false;
			for(let k = 0; k < param.length; k++){
				let par = param[k];
			
				if(par.type == "trans_bp"){
					for(let m = 0; m < par.eqn_appear.length; m++){
						let info = par.eqn_appear[m].eqn_info;
						if(info.p == bl.p && info.cl == bl.cl && info.i == bl.tra[j]){
							if(fl == true) error("branching parameter appears more than once");
							
							par_list.push(k);
							
							fl = true;
						}
					}
					if(fl == true) break;
				}
			}
			if(fl == false){ error("Could not get bp"); break;}
		}
		
		let N = par_list.length;
		if(N != bl.tra.length){
			error("Could not set branching auto");
		}
		else{
			/// Checks that the dependencies on the different parameters are all the same
			let dep = param[par_list[0]].dep;
			let ill = false;
			
			for(let j = 1; j < N; j++){
				if(equal_vec(dep,param[par_list[j]].dep) != true) ill = true;
			}		
		
			if(ill == true){
				for(let j = 0; j < N; j++){
					let par = param[par_list[j]];
					let eqn_info = par.eqn_appear[0].eqn_info;
					
					if(j == 0){
						model.warn.push({mess:"The dependency on the branching probabilties is not the same", mess2:"View parameter '"+par.name+"'", eqn_info:eqn_info, eqn_type:"trans_bp", warn_type:"Equation"});
					}
					else{
							model.warn.push({mess:"", mess2:"View parameter '"+par.name+"'", eqn_info:eqn_info, eqn_type:"trans_bp", warn_type:"Equation"});
					}
				}
			}
			else{
				bl.dep = dep;
			
				let j = 0; while(j < N && param[par_list[j]].auto_value != true) j++;
				if(j == N) j = N-1;
				
				for(let k = 0; k < N; k++){
					if(k == j) param[par_list[k]].auto_value = true;
					else param[par_list[k]].auto_value = false;
				}
			}
			
			branch_list[i].par_list = par_list;
		}
	}
	
	return branch_list;
}
	

/// Sets the parameter which is automatically calculated during a branching event
function set_branch_auto(th)
{
	let fl = false;
	
	let param = model.param;
	for(let i = 0; i < model.branch_list.length; i++){
		let par_list = model.branch_list[i].par_list;
		let j = find_in(par_list,th);
		if(j != undefined){
			let N = par_list.length;
			for(let k = 0; k < N; k++){
				if(k == j) param[par_list[k]].auto_value = true;
				else param[par_list[k]].auto_value = false;
			}
			fl = true;
			break;
		}
	}
	
	if(fl == false) error("Could not change auto branch");
}


/// Automatically calculates branching probability
function branch_auto_calculate()
{
	let param = model.param;

	for(let i = 0; i < model.branch_list.length; i++){
		let bl = model.branch_list[i];
		let par_list = bl.par_list;
		let N = par_list.length;
		
		let dep = bl.dep;
		if(dep != undefined){
			if(dep.length == 0){
				let sum = 0;
				let par_st;
				let fl = false;
				for(let j = 0; j < N; j++){
					let par = param[par_list[j]];
					if(par.auto_value == true) par_st = par;
					else{
						if(isNaN(par.value)) fl = true;
						else sum += Number(par.value);
					}
				}				
				
				if(par_st == undefined) error("Problem calculating branching");
				else{
					if(fl == false){
						let val = Number((1-sum).toPrecision(8));
						if(val > -TINY && val < TINY) val = 0;
						par_st.value = val;
					}
				}
			}
			else{
				let flag = false;   // Makes sure all non-auto values are set
				let jst;
				for(let j = 0; j < N; j++){
					let par = param[par_list[j]];
					if(par.auto_value == false){
						if(par.set != true) flag = true;
					}
					else jst = j;
				}
			
				if(flag == false){
					let par = param[par_list[0]];
					
					let value = par.value;
					let dim = get_dimentions(value);
					let ele_list = get_element_list(value,dim);
					
					let par_set = param[par_list[jst]];
					par_set.set = true;
					
					for(let k = 0; k < ele_list.length; k++){
						let sum = 0;
						let fl = false;
						for(let j = 0; j < N; j++){
							if(j != jst){
								let par = param[par_list[j]];
								sum += Number(get_element(par.value,ele_list[k]));
							}
						}
					
						let val = Number((1-sum).toPrecision(8));
						if(val > -TINY && val < TINY) val = 0;
						
						set_element(par_set.value,ele_list[k],val)
					}
				}
			}
		}
	}
}


/// Sets an unspecified prior
function unset_prior(par_type)
{
	var eb = {te:"", type:"prior"};
	var value = {min_eqn:copy(eb), max_eqn:copy(eb), mean_eqn:copy(eb), shape_eqn:copy(eb), cv_eqn:copy(eb), sd_eqn:copy(eb), alpha_eqn:copy(eb), beta_eqn:copy(eb)};
	
	let te = select_str;
	
	switch(par_type){  // Loads up standard priors for some parameters
	case "trans_bp": te = "flat"; break;
	case "correlation": te = "uniform"; value.min_eqn.te = "-0.9"; value.max_eqn.te = "0.9"; break;
	}	
				
	return {type:{te:te}, value:value};
}


/// Checks that knot times are set correctly 
function check_knot_times(te)
{		
	let spl = te.split(",");
	for(let j = 0; j < spl.length; j++){
		let te = spl[j];
		if(!(te == "start" && j == 0) && !(te == "end" && j == spl.length-1)){
			if(isNaN(te)){
				return "In '"+spl+"' the value '"+te+"' must be a number";
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
	let par_copy = copy(par);
	par_copy.spline.knot = knot;
	copy_param_info(par,par_copy);
}


/// CHecks all the derived parameters
function check_derived_param_all(param)
{
	for(let i = 0; i < model.derive.length; i++){
		let eqn1 = model.derive[i].eqn1;
		let eqn2 = model.derive[i].eqn2;
		
		let wa = check_derived_param(eqn2,eqn1,param);
		if(wa != "success"){
			model.warn.push({mess:"Error in derived expression", mess2:wa, warn_type:"Equation", eqn_info:{i:i}, eqn_type:eqn2.type});
		}
	}	
}


/// Checks that derived equation only contains only parameters from the model 
function check_derived_param(eqn,eqn_param,param)
{
	if(param == undefined) param = model.param;
	
	extract_equation_properties(eqn);
	
	if(eqn.warn.length > 0) return eqn.warn[0].te;

	for(let i = 0; i < eqn.param.length; i++){
		let par = eqn.param[i];
		
		let j = 0;
		while(j < param.length && par_same(par,param[j]) == false) j++;
		
		if(j == param.length){
			if(find(param,"name",par.name) != undefined){
				return "Parameter '"+par.full_name+"' has a different dependency in the model";
			}
			else return "Parameter '"+par.full_name+"' is not found in the model";
		}
	}
	
	// Checks that the dependencies for a derived quantity agree with equation
	if(eqn_param == undefined){ error("should be defined"); return;}
	
	extract_equation_properties(eqn_param);

	if(eqn_param.param.length == 1){
		let par_param = eqn_param.param[0];
		for(let i = 0; i < param.length; i++){
			let par = param[i];
			if(par.type != "derive_param"){
				if(par_param.name == param[i].name){
					return "The model already contains a parameter with name '"+par_param.name+"'";
				}
			}
		}
	
		let dep = equation_dep(eqn).sort();
		let dep_param = equation_dep(eqn_param).sort();
		
		if(find_in(dep,"a") != undefined){
			return "Derived quantities cannot have an age dependency";
		}
		
		if(dep.length != dep_param.length){
			let par = eqn_param.param[0];
			let par2 = {name:par.name, dep:dep};
			let full_name = param_name(par2);
			return "Problem with dependency - Expect '"+full_name+"' not '"+par.full_name+"'.";
		}
	}
	
	return "success";
}


/// Returns true if two parameters the same
function par_same(par1,par2)
{
	if(par1.name != par2.name) return false;
	
	if(par1.dep.length != par2.dep.length) return false;
	
	for(let i = 0; i < par1.dep.length; i++){
		if(par1.dep[i] != par2.dep[i]) return false;
	}
	
	return true;
}


/// Returns text with the name of a parameter
function param_name(par)
{
	let name = par.name;
	if(par.dep.length > 0){
		let flag = false;
		for(let i = 0; i < par.dep.length; i++){
			let dep = par.dep[i];
			if(dep == "t") name += "(t)";
			else{
				if(dep == "a") name += "(a)";
				else{
					if(flag == false){ name += "_"; flag = true;}
					if(i != 0) name += ",";
					name += dep;
				}
			}
		}
	}
	
	return "<e>"+name+"</e>";
}


/// Returns text with the name of a parameter with the indexes substituted by values
function param_name_index(par,index)
{
	let name = par.name;
	
	if(par.dep.length > 0){
		let flag = false;
		for(let i = 0; i < par.dep.length; i++){
			let dep = par.dep[i];
			let val = par.list[i][index[i]];
			if(dep == "t") name += "("+val+")";
			else{
				if(dep == "a") name += "("+")";
				else{
					if(flag == false){ name += "_"; flag = true;}
					if(i != 0) name += ",";
					name += val;
				}
			}
		}
	}
	
	return name;
}


/// Updates prior_split after it has been editted
function update_prior_split()
{
	let bubpr = inter.bubble.prior;
	let pri = get_element(inter.edit_param.prior_split,bubpr.pindex);

	pri.type = bubpr.type;
	pri.value = bubpr.value;
			
	set_element(inter.edit_param.value,bubpr.pindex,get_prior_string(pri));
			
	//init_param();
}

	
/// Gets the full parameter name 
function get_full_parameter_name(name)
{
	let i = find(model.param,"name",name);
	if(i == undefined){
		error("could not find parameter '"+name+"'");
		return name;
	}
	
	return param_name(model.param[i]);
}


/// Determines if two objects are the same or not
function equal_ob(ob1, ob2)
{
	if(JSON.stringify(ob1) != JSON.stringify(ob2)) return false;
	return true;
}


/// Checks that parameters are correctly specified
function check_param()
{
	for(let i = 0; i < model.param.length; i++){
		let par = model.param[i];
		
		if(par.dep.length > 0){
			if(equal_ob(par.list,par_find_list(par)) == false) error("List not correct");
			if(equal_ob(par.comb_list,generate_comb_list(par.list)) == false) error("Comb list not correct");
		}
		
		if(par.dep.length == 0){
			if(par.list != undefined){ error("List should not be specified");}
			if(par.comb_list != undefined) error("comb_list should not be specified");
			if(par.prior_split != undefined) error("prior_split should not be specified");
			if(par.prior_split_check != undefined) error("prior_split_check should not be specified");
			if(par.set != undefined) error("Set should not be specified");
		}
		
		if(par.dep.length > 0){
			let temp = par_find_template(par.list);
		
			let dim = get_dimentions(temp);
	
			if(equal_vec(dim,get_dimentions(par.value)) == false){
				error("Problem with value dimension");
			}
		
			if(equal_vec(dim,get_dimentions(par.prior_split)) == false) error("Problem with prior_split dimension");
		}
	}
}



/// Creates splines in the model 
function create_spline()
{
	model.spline = [];
	
	for(let th = 0; th < model.param.length; th++){
		let par = model.param[th];

		if(par.time_dep == true){
			let list = par_find_list(par);
			list.pop();
			let comb_list = generate_comb_list(list);
			
			for(let k = 0; k < comb_list.length; k++){
				comb_list[k].name = get_spline_name(par,comb_list[k].index,list);
			}
	
			model.spline.push({name:par.name+"(t)", th:th, comb_list:comb_list});
		}			
	}
}


/// Gets the name of a spline from a parameter and a set of indices
function get_spline_name(par,ind,list)
{
	let name = par.name;
	if(ind.length > 0){
		name += "_";
		for(let d = 0; d < ind.length; d++){
			if(d != 0) name += ","; name += list[d][ind[d]];
		}
	}
	name += "(t)";

	return name;
}


/*
/// Determines if the parameter is "dist" giving the distance between compartments
function is_dist(par)
{
	if(par.name.toLowerCase() == "dist" && par.dep.length == 2){
		let index = remove_prime(par.dep[0]);
		if(index == remove_prime(par.dep[1])) return true;
	}
	return false;
}
*/


/// If any parameters are called "dist" then put values for the distance matrix into them
function set_dist(par)
{
	let claa = get_cla_from_index(model,par.dep[0]);
	if(claa == undefined){ error("Classification is undefined"); return;}
	let value = par.value;
	
	let N = claa.comp.length;

	switch(claa.camera.coord){
	case "cartesian":
		{
			let px=[], py=[];
	
			for(let i = 0; i < N; i++){
				px.push(claa.comp[i].x); py.push(claa.comp[i].y); 
			}
			
			for(let j = 0; j < N; j++){
				for(let i = j; i < N; i++){
					if(i == j) value[j][i] = 0;
					else{
						let dx = px[i]-px[j], dy = py[i]-py[j];
						value[j][i] = Math.sqrt(dx*dx + dy*dy);
					}
					
					value[i][j] = value[j][i];
				}
			}
		}
		break;
		
	case "latlng":
		{
			let lat=[], lng=[];
			
			for(let i = 0; i < N; i++){
				let p = transform_latlng_inv(claa.comp[i].x,claa.comp[i].y);
				lat.push(p.lat*Math.PI/180);
				lng.push(p.lng*Math.PI/180);
			}
			
			let r = 6371;

			for(let j = 0; j < N; j++){
				for(let i = j; i < N; i++){
					if(i == j) value[j][i] = 0;
					else{
						let si = Math.sin(0.5*(lat[j]-lat[i]));
						let si2 = Math.sin(0.5*(lng[j]-lng[i]))
						value[j][i] = 2*r*Math.asin(Math.sqrt(si*si + Math.cos(lat[i])*Math.cos(lat[j])*si2*si2));
					}
					
					value[i][j] = value[j][i];
				}
			}
		}
		break;
	}
	
	par.variety = "const";	
	par.set = true;
}
	

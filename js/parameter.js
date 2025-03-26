"use strict";
// Functions relating to parameters

/// Generates a list of all model parameters
// Goes through the exisiting list and copies any info on values/priors
function update_model()
{
	model.warn = [];
	
	check_data_valid_all("siminf");

	let lists = generate_parameter_list();           // Generates a list of all model parameters

	if(false){ pr("lists"); pr(lists);}

	let ielist = lists.ind_eff_list;       
  let felist = lists.fix_eff_list;       
	let par_list = lists.par_list;
	
	update_model_ind_eff(ielist);                    // Updates individual effect groupings in model
	
	update_model_fix_eff(felist);                    // Updates fixed effects in model 
	
	add_param_ppc_factor(par_list);                  // Adds any parameter for ppc factors

	add_ind_eff_param(par_list);                     // Adds individual effect parameters (i.e. variances) to model
	
	update_model_param(par_list);                    // Updates the model parameters

	for(let th = 0; th < model.param.length; th++){  // Sets the distance matrix
		let par = model.param[th];	
		if(par.name == dist_matrix_name){
			set_dist(par);                            
		}
	}
	
	set_generate_pos();
	
	if(model.warn.length > 0) err_warning();
}


/// Sets if it is possible to generate data from results (by seeing if the model had been altered)
function set_generate_pos()
{
	let sp = model.species;
	let sp2 = sim_result.species;
	
	for(let p = 0; p < sp.length; p++){
		sp[p].generate_pos = false;
	}	
	
	if(sp2 && sp2.length == sp.length){
		for(let p = 0; p < sp.length; p++){
			let fl = 0;
			
			if(sp[p].type != sp2[p].type) fl = 1;
			
			if(sp[p].cla.length != sp2[p].cla.length) fl = 2;
			else{
				for(let cl = 0; cl < sp[p].cla.length; cl++){
					let claa = sp[p].cla[cl], claa2 = sp2[p].cla[cl];
					
					if(claa.name != claa2.name || claa.index != claa2.index) fl = 3;
					else{
						if(claa.comp.length != claa2.comp.length) fl = 4;
						else{
							for(let c = 0; c < claa.comp.length; c++){
								let co = claa.comp[c], co2 = claa2.comp[c];
								if(co.name != co2.name || co.branch != co2.branch) fl = 5;
								if(co.infected.check != co2.infected.check) fl = 6;
							}
						}
		
						if(claa.tra.length != claa2.tra.length) fl = 7;
						else{
							for(let tr = 0; tr < claa.tra.length; tr++){
								let tra = claa.tra[tr], tra2 = claa2.tra[tr];
								if(tra.i != tra2.i) fl = 8;
								if(tra.f != tra2.f) fl = 9;
								if(tra.type != tra2.type) fl = 10;
								if(tra.bp_set != tra2.bp_set) fl = 11;
								
								let br = false; if(tra.branch == true) br = true;
								let br2 = false; if(tra2.branch == true) br2 = true;
								if(br != br2)	fl = 12;
							}				
						}
					}
				}
			}
			
			if(fl == 0) sp[p].generate_pos = true;
			else sp[p].generate_pos = false;
		}
	}		
}


/// Genrates a list of all model parameters
function generate_parameter_list()
{
	let eq_list = find_equation_list();
	
	if(false){ pr("eq_list"); pr(eq_list);}
	
	let list = [];
	let ind_eff_list = [];
	let fix_eff_list = [];
	
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			if(claa.ncomp == 0){
				model.warn.push({mess:"Classification problem", mess2:"Classification '"+claa.name+"' does not contain any compartments", p:p, cl:cl, warn_type:"NoClaComp"});
			}
			for(let c = 0; c < claa.ncomp; c++){
				let co = claa.comp[c];
				if(co.name == ""){
					model.warn.push({mess:"Missing compartment name", mess2:"Name has not been specified in classification '"+claa.name+"'", p:p, cl:cl, c:c, warn_type:"MissingComp"});
				}
				else{
					if(co.col == "#FFFFFF"){
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
			
		case "mut_rate": 
			mess = "For a genetic mutation rate";
			break;
			
		case "seq_var": 
			mess = "For a genetic sequence variation";
			break;
			
		case "comp_prob":
			mess = "For a compartmental probability";
			break;
		
		case "trans_prob":
			mess = "For a transition probability";
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
			//let par = { name:eqn.fix_eff[i].name, time_dep: false, dep:[], eqn_appear:[], type: "fixed effect"};
			let par = { name:"ν^"+eqn.fix_eff[i].name, time_dep: false, dep:[], eqn_appear:[], type: "fixed effect"};
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
	case "mut_rate": return "as a genetic mutation rate";
	case "seq_var": return "as a genetic sequence variation";
	case "trans_bp": return "as a branching probability";
	case "trans_mean": return "in a transition mean";
	case "trans_rate": return "in a transition rate";
	case "trans_shape": return "in a transition shape";
	case "trans_scale": return "in a transition scale";
	case "trans_cv": return "in a transition coefficient of variation";
	case "comp_prob": return "in a compartmental probability";
	case "reparam": return "in a reparametersation";
	default: error("option not pos:"+eqn.type); break;	
	}
}


/// Finds the tensor template for a parameter (this is a multi-dimensional tensor with unset elements)
function par_find_template(list)
{
	let ndep = list.length;
	
	let list_len = [];
	for(let d = 0; d < ndep; d++){
		let k = 0;
		while(k < list[d].length && list[d][k] != undefined) k++;
		list_len[d] = k;
	}
	
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
			if(index[i] >= list_len[i]){
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


/// Gets the dimentions
function get_value_dim(par)
{
	let dim = [];
	if(Array.isArray(par.value) && par.comb_list.length > 0){
		for(let d = 0; d < par.list.length; d++) dim.push(par.list[d].length);
	}
	return dim;
}


/// Gets a description of value to send back to inferface
function get_value_desc(par)
{
	let te;
	
	if(Array.isArray(par.value)){
		if(par.comb_list.length == 0){
			return no_elements;
		}
		
		let dim = [];
		for(let d = 0; d < par.list.length; d++) dim.push(par.list[d].length);
		
		//let dim = get_list_len(par.list);//  get_dimensions(par.value);
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
		
		if(par.set == true && tot < 10){
			te = JSON.stringify(par.value);
			te = te.replace(/\"/g," ")
			te = te.replace(/null/g,"Unset");
		}
	}
	
	return te;
}


/// Gets the dimensions of the object
function get_dimensions(ob)
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
		if(hash_find(claa.hash_comp,name) != undefined){
			return claa.index;
		}
	}
	return "";
}
	
	
/// Updates individual effect groups stored in model
function update_model_ind_eff(ielist)
{
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		
		for(let i = 0; i < sp.ind_eff_group.length; i++){
			let ieg = sp.ind_eff_group[i];
			for(let j = 0; j < ieg.ie_list.length; j++) ieg.ie_list[j].used = false;
		}
	
		for(let k = 0; k < ielist.length; k++){
			let ie = ielist[k];
			let flag = false;
			
			if(ie.p_name == sp.name){
				for(let i = 0; i < sp.ind_eff_group.length; i++){
					let ieg = sp.ind_eff_group[i];
					let ie_list = ieg.ie_list;
					for(let j = 0; j < ie_list.length; j++){
						if(ie.name == ie_list[j].name){
							ie_list[j].used = true;
							flag = true;
						}
					}
				}
			
				if(flag == false){
					sp.ind_eff_group.push({p_name:ie.p_name, A_matrix:{check:false, loaded:false, value:[], ind_list:[], pedigree:false, sire_list:[], dam_list:[]}, ie_list:[{name:ie.name}]});
				}
			}
		}
		
		// Deletes unused individual effects
		let i = 0;
		while(i < sp.ind_eff_group.length){
			let ieg = sp.ind_eff_group[i];
			let ie_list = ieg.ie_list;
				
			let j = 0;
			while(j < ie_list.length){
				if(ie_list[j].used == false) ie_list.splice(j,1);
				else j++;
			}
			
			if(ie_list.length == 0) sp.ind_eff_group.splice(i,1);
			else i++;
		}
	}
	
	// Checks that there are no individual effects in population-based models
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
	
		if(sp.type == "Population" && sp.ind_eff_group.length > 0){
			let name = sp.ind_eff_group[0].ie_list[0].name;
			model.warn.push({mess:"Cannot have individual effects", mess2:"A population-based model cannot have individual effects such as '"+name+"'", warn_type:"Equation"});
		}
	}
}


/// Updates fixed effects stored in model
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
			let name = sp.fix_eff[0].name;
			model.warn.push({mess:"Cannot have fixed effects", mess2:"A population-based model cannot have fixed effects such as '"+name+"'", warn_type:"Equation"});
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
			for(let j = 0; j < ieg.ie_list.length; j++){
				for(let jj = j; jj < ieg.ie_list.length; jj++){
					let ie1 = ieg.ie_list[j].name, ie2 = ieg.ie_list[jj].name;
					let sup = "^"+ie1+","+ie2;
					
					let par;
					if(j == jj){
						par = { name:"Ω"+sup, time_dep: false, dep:[], eqn_appear:[], type: "variance"};
					}
					else{
						par = { name:"ω"+sup, time_dep: false, dep:[], eqn_appear:[], type: "correlation"};
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
	check_derived_param_all(model.param);

	if(model.warn.length > 0) err_warning();
	
	let param_old = model.param;
				
	let param=[];
	for(let i = 0; i < par_list.length; i++){
		let par = par_list[i];
		
		let j;                                         // Determines if parameter already exists
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
		
		if(j < param_old.length){                      // Copies information from existing parameter definition	
			param.push(copy_param_info(par,param_old[j]));
			add_dependent_param(param,par_list);
		}
		else{                                          // Create new parameter
			param.push(create_new_param(par,"normal"));
		}
	}

	param.sort( function(a, b){                      // Alphabetically sorts parameters
		if(a.dep.length > b.dep.length) return 1; 
		if(a.dep.length == b.dep.length){
			if(a.name.toLowerCase() > b.name.toLowerCase()) return 1;
			else return -1;	
		}
		return -1;});
	
	if(model.warn.length > 0) err_warning();
	
	model.param = param;
}


/// Goes through reparameterisation and distribution and adds to parameter list
function add_dependent_param(param,par_list)
{
	let par = param[param.length-1];
	
	let vari = par.variety;
	
	switch(vari){
	case "reparam":
		if(par.dep.length == 0) add_ele_param(par.value,vari,par_list,{par_name:par.name});
		else{
			if(par.reparam_eqn_on) add_ele_param(par.reparam_eqn,"reparam_eqn",par_list,{par_name:par.name});
			else add_to_par_list(par.reparam_param_list,par_list);
		}
		break;
		
	case "dist":
		{
			if(par.dep.length == 0 || par.prior_split_check.check == false){	
				add_distribution_eqn(par.prior,vari,par_list,{par_name:par.name});	
			}
			else{
				add_to_par_list(par.prior_param_list,par_list);
			}
		}
		break;
	}
	
	par.full_name = param_name(par);
}


/// Adds parameters defined in reparameterisation or priors
function add_to_par_list(add_param_list,par_list)
{
	for(let i = 0; i < add_param_list.length; i++){
		let par_add = add_param_list[i]; 

		let j = 0; while(j < par_list.length && par_list[j].name != par_add.name) j++;
		if(j == par_list.length) par_list.push(par_add);
		else{
			for(let k = 0; k < par_add.eqn_appear.length; k++){
				par_list[j].eqn_appear.push(par_add.eqn_appear[k]);
			}
		}
	}
}
		

/// Adds any parameters from a distribution to the list of equations
function add_distribution_eqn(dist,vari,par_list,eqn_info,op)
{
	let val = dist.value;
	
	switch(dist.type.te){
	case "fix": 
		add_ele_param(val.mean_eqn.te,vari,par_list,eqn_info,op);
		break;
		
	case "uniform":
		add_ele_param(val.min_eqn.te,vari,par_list,eqn_info,op);
		add_ele_param(val.max_eqn.te,vari,par_list,eqn_info,op);
		break;
		
	case "exp": 
		add_ele_param(val.mean_eqn.te,vari,par_list,eqn_info,op);
		break;
		
	case "normal": 
		add_ele_param(val.mean_eqn.te,vari,par_list,eqn_info,op);
		add_ele_param(val.sd_eqn.te,vari,par_list,eqn_info,op);
		break;
	
	case "log-normal": 
		add_ele_param(val.mean_eqn.te,vari,par_list,eqn_info,op);
		add_ele_param(val.cv_eqn.te,vari,par_list,eqn_info,op);
		break;
		
	case "gamma":
		add_ele_param(val.mean_eqn.te,vari,par_list,eqn_info,op);
		add_ele_param(val.cv_eqn.te,vari,par_list,eqn_info,op);
		break;
		
	case "bernoulli":
		add_ele_param(val.mean_eqn.te,vari,par_list,eqn_info,op);
		break;
					
	case "beta":
		add_ele_param(val.alpha_eqn.te,vari,par_list,eqn_info,op);
		add_ele_param(val.beta_eqn.te,vari,par_list,eqn_info,op);
		break;
	
	case select_str:
		break;
		
	default: error("not op:"+dist.type.te); break;
	}
}


/// Adds any parameters from a text element to the parameter list
function add_ele_param(ele,variety,par_list,eqn_info,op)
{
	if(isNaN(ele)){
		let eqn = create_equation(ele,variety);
		eqn.eqn_info = eqn_info;
		
		if(op == "single eqn_appear"){  // In this case only want one eqn_appear 
			for(let i = 0; i < eqn.param.length; i++){
				let par = eqn.param[i];
				let j = 0; while(j < par_list.length && par_list[j].name != par.name) j++;
				
				if(j == par_list.length){
					add_par_to_list(par,par_list,copy(eqn),"Reparameterisation");
				}
			}
		}
		else{
			for(let i = 0; i < eqn.param.length; i++){
				add_par_to_list(eqn.param[i],par_list,copy(eqn),"Reparameterisation");
			}
		}
	}
}


/// Gets a list of all parameters in a reparameterisation
function get_reparam_param_list(par) 
{
	par.reparam_param_list=[];
	if(par.variety == "reparam" && par.dep.length > 0 && par.reparam_eqn_on != true){
		let eqn_info = {par_name:par.name};
		let par_list = [];
		let vari = par.variety;
		
		for(let i = 0; i < par.comb_list.length; i++){
			let comb = par.comb_list[i];
			let ele = get_element(par.value,comb.index);
			eqn_info.index = comb.index;
			add_ele_param(ele,vari,par_list,eqn_info,"single eqn_appear");
		}
	
		par.reparam_param_list = par_list;
	}
}


/// Gets a list of all parameters in prior
function get_prior_param_list(par) 
{
	par.prior_param_list=[];
	
	if(par.dep.length > 0 && par.prior_split_check.check == true){
		let eqn_info = {par_name:par.name};
		let par_list = [];
		let vari = par.variety;
		
		for(let i = 0; i < par.comb_list.length; i++){
			let comb = par.comb_list[i];
			let ele = get_element(par.prior_split,comb.index);
				
			add_distribution_eqn(ele,vari,par_list,eqn_info,"single eqn_appear");	
		}

		par.prior_param_list = par_list;
	}
}


/// Creates a new parameter
// par must have the properties: name, dep, type, time_dep, (optionally eqn_appear)
// On top of these new properties are added: 
// list: shows the potential values for the dependencies
// comb_list: this is a list giving every possible combination of dependency  
function create_new_param(par,variety)
{
	if(par.name == undefined) error("Parameter name must be set");
	if(par.dep == undefined) error("Parameter dependency must be set");
	if(par.type == undefined) error("Parameter type must be set");
	if(par.time_dep == undefined) error("Parameter time dependency must be set");

	if(par.type == "param factor") variety = "const"; 

	par.variety = variety;
	
	par.pri_pos = set_pri_pos(par.type);

	let knot = ["start","end"];                      // Initialises a spline
	
	let on = false;
	if(par.time_dep == true) on = true;
	
	par.spline = { on:on, smooth:default_smooth(), time_dep:par.time_dep, knot:knot}; 
		
	par.sim_sample = {check:true};
		 
	par.prior = unset_prior(par.type);
	
	par.reparam_eqn_on = false;
	par.reparam_eqn = "";
	
	par_set_default(par);

	//par.dist_matrix = {check:false}; 
	
	par.full_name = param_name(par);

	par.reparam_param_list = [];
	par.prior_param_list = [];

	return par;
}


/// The default value for smoothing
function default_smooth()
{
	return {check:false, type:{value:"Log-Normal"}, value:0.5};
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
		
		par.list = list;
		par.comb_list = comb_list;
		par.value = value;
		par.prior_split = prior_split;
		par.prior_split_check = {check:false};
		par.set = false;
	}
}


/// Sets default values for prior
function set_default_prior_split(par)
{
	let prior_split = par.prior_split;
	let dim = get_dimensions(prior_split);
	let ele_list = get_element_list(prior_split,dim);
		
	for(let k = 0; k < ele_list.length; k++){	
		let pr_sp = unset_prior(par.type);	
		set_element(prior_split,ele_list[k],pr_sp);
	}
}


/// Sets default values for tensor
function set_default_value(par)
{
	let value = par.value;
	let dim = get_dimensions(value);
	let ele_list = get_element_list(value,dim);
		
	for(let k = 0; k < ele_list.length; k++){	
		set_element(value,ele_list[k],0);
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
	par.label_info = old.label_info;
	par.set = old.set;
	//par.dist_matrix = old.dist_matrix;
	par.import_line = old.import_line;
	par.reparam_param_list = old.reparam_param_list;
	par.prior_param_list = old.prior_param_list;
	par.reparam_eqn_on = old.reparam_eqn_on;
	par.reparam_eqn = old.reparam_eqn;
	
	if(par.dep.length == 0){	
		par.value = old.value;
	}
	else{
		let list = par_find_list(par);
		let list_old = old.list;
		
		// Checks if list is the same as the 
		let dif = false;
		
		if(!list_old) dif = "unset";
		else{
			if(list.length != list_old.length) dif = true;
			else{
				for(let k = 0; k < list.length; k++){
					let lik = list[k];
					let lik_old = list_old[k];
					if(lik.length != lik_old.length){ dif = true; break;}
					for(let j = 0; j < lik.length; j++){
						if(lik[j] != lik_old[j]){ dif = true; break;}
					}
				}
			}
		}
		
		if(dif != false){
			let temp = par_find_template(list);
			let comb_list = generate_comb_list(list);
		
			let value = copy(temp);
			let prior_split = copy(temp);
		
			let prior = copy(old.prior);
			
			let old_fl = false;
			
			// Attempts to copy across existing values
			if(dif == true){
				let map=[];
				for(let k = 0; k < par.dep.length; k++){
					map[k]=[];
					
					let list_new = list[k];
					if(k < old.list.length){
						let list_old = old.list[k];
						
						let used=[];
						for(let ij = 0; ij < list_new.length; ij++){
							let val = list_new[ij];
							let m = 0; while(m < list_old.length && !(list_old[m] == val && used[m] == undefined)) m++;
							if(m < list_old.length){ map[k][ij] = m; used[m] = true;}
						}
					}
				}
				
				let dim = get_dimensions(temp);
				let ele_list = get_element_list(temp,dim);
				
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
						old_fl = true;
						val = get_element(old.value,ind_old);
						pr_sp = get_element(old.prior_split,ind_old);
					}
					else{
						val = undefined;
						val_eqn = create_equation("","reparam");
						pr_sp = unset_prior(par.type);
					}
					
					set_element(value,ind,val);
					set_element(prior_split,ele_list[k],pr_sp);
				}
			}
			
			par.list = list;
			par.comb_list = comb_list;
			par.value = value;
			par.prior_split = prior_split;
			if(old_fl == false) par.set = false;
		}
		else{
			par.list = old.list;
			par.comb_list = old.comb_list;
			par.value = old.value;
			par.prior_split = old.prior_split;
		}
		
		par.prior_split_set = old.prior_split_set;
		par.prior_split_check = old.prior_split_check;
	}
	
	return par;
}


/// Finds the possible set of priors for a parameter type
function set_pri_pos(type)
{
	let pri_pos = prior_pos;  

	return pri_pos;
}


/// Given a list of dependencies this provides a list of cominations going through each
function generate_comb_list(list)
{
	let ndep = list.length;
	let comb_list = [];
	
	let list_len = [];
	for(let d = 0; d < ndep; d++){
		let k = 0;
		while(k < list[d].length && list[d][k] != undefined) k++;
		list_len[d] = k;
	}
	
	if(ndep == 0) return [{index:[]}];
	
	for(let d = 0; d < ndep; d++){
		if(list_len[d] == 0){
			return comb_list;
		}
	}
	
	let index=[];
	for(let d = 0; d < ndep; d++) index[d] = 0;

	let flag;
	do{
		comb_list.push({index:copy(index)});
		
		let i = 0;
		do{
			flag = false;
			
			index[i]++; 
			if(index[i] >= list_len[i]){
				index[i] = 0; 
				i++; 
				flag = true;
			}
		}while(i < ndep && flag == true);
	}while(flag == false);
	
	return comb_list;
}


/// Sets an unspecified prior
function unset_prior(par_type)
{
	var eb = {te:"", type:"prior", mode:"param only"};
	var value = {min_eqn:copy(eb), max_eqn:copy(eb), mean_eqn:copy(eb), shape_eqn:copy(eb), cv_eqn:copy(eb), sd_eqn:copy(eb), alpha_eqn:copy(eb), beta_eqn:copy(eb)};
	
	let te = select_str;
	
	switch(par_type){  // Loads up standard priors for some parameters
	case "correlation": te = "uniform"; value.min_eqn.te = "-0.9"; value.max_eqn.te = "0.9"; break;
	}	
				
	return {type:{te:te}, value:value};
}


/// Checks all the derived parameters
function check_derived_param_all(param)
{
	for(let i = 0; i < model.derive.length; i++){
		let eqn1 = model.derive[i].eqn1;
		let eqn2 = model.derive[i].eqn2;
		
		let res = check_derived_param(eqn2,eqn1,param);
		
		if(res.err == true){
			model.warn.push({mess:"Error in derived expression", mess2:res.msg, warn_type:"Equation", eqn_info:{i:i}, eqn_type:eqn2.type});
		}
	}	
}


/// Checks that derived equation only contains parameters from the model 
function check_derived_param(eqn,eqn_param,param)
{
	if(param == undefined) param = model.param;
	
	extract_equation_properties(eqn);
	
	if(eqn.warn.length > 0) return err(eqn.warn[0].te);

	for(let i = 0; i < eqn.param.length; i++){
		let par = eqn.param[i];
		
		let j = 0;
		while(j < param.length && par_same(par,param[j]) == false) j++;
		
		if(j == param.length){
			if(find(param,"name",par.name) != undefined){
				return err("Parameter '"+par.full_name+"' has a different dependency in the model");
			}
			else return err("Parameter '"+par.full_name+"' is not found in the model");
		}
		else{
			if(param[j].type == "derive_param"){
				return err("Derived parameter '<e>"+eqn_param.te+"</e>' cannot depend on on another derived parameter '"+par.full_name+"'.");
			}
		}
	}
	
	// Checks that the dependencies for a derived quantity agree with equation
	if(eqn_param == undefined) return err("Should be defined");
	
	extract_equation_properties(eqn_param);

	if(eqn_param.param.length == 1){
		let par_param = eqn_param.param[0];
		for(let i = 0; i < param.length; i++){
			let par = param[i];
			if(par.type != "derive_param"){
				if(par_param.name == param[i].name){
					return err("The model already contains a parameter with name '"+par_param.name+"'");
				}
			}
		}
	
		let dep = equation_dep(eqn).sort();
		let dep_param = equation_dep(eqn_param).sort();
		
		if(dep.length != dep_param.length){
			let par = eqn_param.param[0];
			let par2 = {name:par.name, dep:dep};
			let full_name = param_name(par2);
			return err("Problem with dependency - Expect '"+full_name+"' not '"+par.full_name+"'.");
		}
	}
	
	return success();
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
				if(flag == false){ name += "_"; flag = true;}
				if(i != 0) name += ",";
				name += dep;
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
			if(dep == "t") name += "(t="+val+")";
			else{
				if(flag == false){ name += "_"; flag = true;}
				if(i != 0) name += ",";
				name += val;
			}
		}
	}
	
	return name;
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
		
			let dim = get_dimensions(temp);
	
			if(equal_vec(dim,get_dimensions(par.value)) == false){
				error("Problem with value dimension");
			}
		
			if(equal_vec(dim,get_dimensions(par.prior_split)) == false) error("Problem with prior_split dimension");
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


/// If any parameters are called "dist" then put values for the distance matrix into them
function set_dist(par)
{
	let claa = get_cla_from_index(model,par.dep[0]);
	if(claa == undefined){ error("Classification is undefined"); return;}
	let value = par.value;
	
	let N = claa.comp.length;

	// Works out if we need to calculate the full matrix
	if(N > MATRIX_COMP_MAX && N*N > MATRIX_PLOT_MAX){
		if(N*N > ELEMENT_MAX){
			N = Math.floor(ELE_REDUCE_FAC*Math.sqrt(ELEMENT_MAX));
		}
	}
	if(N > value.length) N = value.length;
	
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
				let co = claa.comp[i];
				
				let x, y;
				if(co.type == "boundary"){ x = co.xmid; y = co.ymid;}
				else{ x = co.x; y = co.y;}
				
				let p = transform_latlng_inv(x,y);
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
						value[j][i] = precision(2*r*Math.asin(Math.sqrt(si*si + Math.cos(lat[i])*Math.cos(lat[j])*si2*si2)),5);
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
	

/// Finds a list of all the equations in the model
function find_equation_list(all_param)
{
	let eqn_list = [];
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			for(let i = 0; i < claa.ntra; i++){
				let eqn_info = {p:p,cl:cl,i:i};
				
				let traa = claa.tra[i];
			
				let val = traa.value;
			
				if(val.bp_eqn != undefined){
					if((traa.branch_select == true || traa.all_branches) && (all_param == true || traa.branch == true)){
						add_equation_to_list(eqn_list,val.bp_eqn,eqn_info);
					}
				}
				
				if(val.mean_eqn != undefined){
					if(all_param == true || traa.type == "gamma" || traa.type == "erlang" || traa.type == "log-normal" ||  traa.type == "period"){
						add_equation_to_list(eqn_list,val.mean_eqn,eqn_info);
					}
				}
				
				if(val.rate_eqn != undefined){
					if(all_param == true || traa.type == "exponential"){
						add_equation_to_list(eqn_list,val.rate_eqn,eqn_info);
					}
				}
				
				if(val.scale_eqn != undefined){
					if(all_param == true || traa.type == "weibull"){
						add_equation_to_list(eqn_list,val.scale_eqn,eqn_info);
					}
				}
				
				if(val.shape_eqn != undefined){
					if(all_param == true || traa.type == "weibull"){
						add_equation_to_list(eqn_list,val.shape_eqn,eqn_info);
					}
				}
					
				if(val.cv_eqn != undefined){
					if(all_param == true || traa.type == "gamma" || traa.type == "log-normal"){
						add_equation_to_list(eqn_list,val.cv_eqn,eqn_info);
					}
				}
			}
		}
		
		// Adds parameters associated with the observation process
		for(let i = 0; i < sp.inf_source.length; i++){
			let so = sp.inf_source[i];
	
			if(so.error != true){
				let eqn_info = {p:p, i:i};
				
				switch(so.type){
				case "Diag. Test":
					add_equation_to_list(eqn_list,so.spec.Se_eqn,eqn_info);
					add_equation_to_list(eqn_list,so.spec.Sp_eqn,eqn_info);
					break;
				
				case "Population":
					add_comp_filt_eqn_list(so.spec.filter,eqn_list,i,p,so);
					break;
				
				case "Pop. Trans.": case "Transition":
					{
						let filt = so.spec.filter;
						let cl = find(sp.cla,"name",filt.te);
						if(cl == undefined){ error("Cannot find classification"); return;}
						add_trans_filt_eqn_list(filt,eqn_list,i,p,cl);
						add_comp_filt_eqn_list(filt,eqn_list,i,p,so,cl);
					}
					break;
				
				case "Add Ind.":
					{
						let tab = so.table;
						for(let r = 0; r < tab.nrow; r++){
							for(let c = 2; c < tab.ncol; c++){
								add_eqn_filter(tab.ele[r][c],"comp_prob",p,i,2,r,eqn_list);
							}
						}	
					}
					break;

				case "Compartment":
					{
						let tab = so.table;
						for(let r = 0; r < tab.nrow; r++){
							add_eqn_filter(tab.ele[r][2],"comp_prob",p,i,2,r,eqn_list);
						}	
					}
					break;
					
				case "Genetic":
					add_equation_to_list(eqn_list,so.spec.mut_rate_eqn,eqn_info);
					add_equation_to_list(eqn_list,so.spec.seq_var_eqn,eqn_info);
					break;
				}
			}
		}		
	}
	
	// Goes through derived parameters
	for(let i = 0; i < model.derive.length; i++){
		let der = model.derive[i];
		
		let eqn_info = {i:i};
		
		add_equation_to_list(eqn_list,der.eqn1,eqn_info);		
	}
	
	return eqn_list;
}


/// Adds any equations from the compartmental filter to the equation list
function add_comp_filt_eqn_list(filt,eqn_list,i,p,so,cl_sel)
{
	let cla = filt.cla;
	for(let cl = 0; cl < cla.length; cl++){
		if(cl != cl_sel){
			switch(cla[cl].radio.value){
			case "ObsMod":
				{
					let comp = cla[cl].comp;
					for(let c = 0; c < comp.length; c++){
						let eqn_info2 = {p:p,i:i,cl:cl,c:c};
						add_equation_to_list(eqn_list,comp[c].prob_eqn,eqn_info2);
					}
				}
				break;

			case "File":
				{
					let tab = so.table;
					let c = find_in(tab.heading,model.species[p].cla[cl].name);
					if(c == undefined) error("c should not be undef");
					for(let r = 0; r < tab.nrow; r++){
						add_eqn_filter(tab.ele[r][c],"comp_prob",p,i,c,r,eqn_list);
					}
				}			
				break;
			}
		}
	}
}


/// Adds any equations from the transition filter to the equation list
function add_trans_filt_eqn_list(filt,eqn_list,i,p,cl)
{
	if(filt.trans_obs_model.value == "on"){
		for(let tr = 0; tr < filt.tra.length; tr++){
			let eqn_info2 = {p:p,i:i,cl:cl,tr:tr};
			if(filt.tra[tr].prob_eqn){
				add_equation_to_list(eqn_list,filt.tra[tr].prob_eqn,eqn_info2);
			}
		}
	}
}


/// Adds equations based on a filter such as
function add_eqn_filter(te,type,p,i,c,r,eqn_list)
{
	//let spl = te.split("|");
	let spl = split_with_bracket(te,"|");
	for(let k = 0; k < spl.length; k++){
		let spl2 = spl[k].split(":");
		if(spl2.length == 2){
			let te2 = spl2[1];
			
			if(isNaN(te2)){
				let eqn_info2 = {p:p,i:i,r:r,c:c};
				
				let eqn = create_equation(te2,type,p,undefined);  
				add_equation_to_list(eqn_list,eqn,eqn_info2);
			}
		}
	}
}

	
/// Updates compartment name in parameters
function param_update_rename_compartment(p,cl,c,old_name,new_name)
{
	let index = model.species[p].cla[cl].index;
	for(let i = 0; i < model.param.length; i++){
		let par = model.param[i];
	
		if(par.list){
			let flag = false;
			for(let j = 0; j < par.dep.length; j++){
				
				if(index == remove_prime(par.dep[j])){
					if(c < par.list[j].length){
						if(par.list[j][c] != old_name) error("Old name not correct");
						par.list[j][c] = new_name;
						flag = true;
					}
				}
			}
			
			if(flag == true){
				par.comb_list = generate_comb_list(par.list);
			}
		}
	}
}


/// Adds equation to a list and appends information about source 
function add_equation_to_list(eqn_list,source,eqn_info)
{
	source.eqn_info = eqn_info;
	eqn_list.push(source);
}


/// Gets the specification for param_mult based on par
function get_param_mult(par)
{
	let dep = copy(par.dep);
	let name = "f~"+par.name;
	let full_name = "<e>"+name;
	if(dep.length > 0){
		full_name += "_";
		for(let j = 0; j < dep.length; j++){
			if(j > 0) full_name += ",";
			full_name += dep[j];
		}
	}
	
	if(par.time_dep == false) dep.push("t");
	full_name += "(t)</e>";
	return {name:name, full_name:full_name, dep:dep, time_dep:true, type:"param factor"};
}


/// Adds factors into the model
function add_param_ppc_factor(par_list)
{
	if(!model.param_factor) return;
	for(let j = 0; j < model.param_factor.length; j++){
		let pf = model.param_factor[j];
	
		// Makes sure parameter is in the model
		let fl = false;
		for(let th = 0; th < par_list.length; th++){
			if(pf.param.name == par_list[th].name){
				if(equal_vec(pf.param.dep,par_list[th].dep)){
					fl = true; break;
				}
			}
		}
		
		if(fl == true){
			par_list.push(pf.f_param);
		}
	}
}

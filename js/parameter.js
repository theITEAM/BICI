"use strict";
// Functions relating to parameters

/// Generates a list of all model parameters
// Goes through the exisiting list and copies any info on values/priors
function update_model()
{
	model.warn = [];

	check_data_valid_all("siminf","warn");

	let lists = generate_parameter_list();           // Generates a list of all model parameters

	if(false){ prr("lists"); prr(lists);}
	
	update_model_param(lists);                    // Updates the model parameters

	set_generate_pos();                              // Sets if possible to generate results
	
	update_desc();                                   // Updates descriptions of parameter

	if(model.warn.length > 0) err_warning();
	
	for(let th = 0; th < model.param.length; th++){
		let par = model.param[th];
		if(par.type != undefined) error("par type should not be set");
		if(par.ndep_cont == undefined) error("ndep_cont should be set"); 
	}
}


/// Updates descriptions of parameter
function update_desc()
{
	for(let i = 0; i < model.param.length; i++){
		let par = model.param[i];
		if(par.ndep_cont > 0){
			par.value_desc = get_value_desc(par);
			par.weight_desc = get_weight_desc(par);
			par.prior_split_desc = get_prior_split_desc(par);
			par.prior_const_desc = get_prior_const_desc(par);
			par.dim = get_value_dim(par);
		}
	}
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


/// Adds information about where parameter comes from 
function add_within(par,eqn)
{
	let wi = {type:eqn.type, eqn_info:eqn.eqn_info};
	add_wi(par,wi);
}


/// Adds within 
function add_wi(par,wi)
{
	if(par.within == undefined) return;
	if(par.within.length == WITHIN_MAX) return;

	let i;
	for(i = 0; i < par.within.length; i++){
		let wi2 = par.within[i];
		if(wi2.type == wi.type && wi2.p == wi2.p && wi2.cl == wi2.cl && wi2.i == wi2.i && 
		  wi2.r == wi2.r && wi2.c == wi2.c && wi2.eso == wi2.eso && wi2.ieg == wi2.ieg && 
			wi2.ie == wi2.ie) break;
	}
	
	if(i == par.within.length) par.within.push(wi);
}


/// Genrates a list of all model parameters
function generate_parameter_list()
{
	let eq_list = find_equation_list();
	
	if(false){ prr("eq_list"); prr(eq_list);}
	
	let list = [];
	let ind_eff_list = [];
	let fix_eff_list = [];
	
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			/*
			if(claa.ncomp == 0){
				add_warning({mess:"Classification problem", mess2:"Classification '"+claa.name+"' does not contain any compartments", p:p, cl:cl, warn_type:"NoClaComp"});
			}
			*/
			
			for(let c = 0; c < claa.ncomp; c++){
				let co = claa.comp[c];
				if(co.name == ""){
					add_warning({mess:"Missing compartment name", mess2:"Name has not been specified in classification '"+claa.name+"'", p:p, cl:cl, c:c, warn_type:"MissingComp"});
				}
				else{
					if(co.col == "#FFFFFF"){
						add_warning({mess:"Compartment colour not set", mess2:"In compartment '"+co.name+"' classification '"+claa.name+"'", p:p, cl:cl, c:c, warn_type:"MissingColour"});
					}
				}
			}
				
			for(let k = 0; k < claa.annotation.length; k++){
				let anno = claa.annotation[k];
				switch(anno.type){
				case "box":
					if(anno.te == ""){
						add_warning({mess:"Missing bounding box name", mess2:"Bounding box name has not been specified in classification '"+claa.name+"'", p:p, cl:cl, k:k, warn_type:"MissingBox"});
					}
					break;
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
			
		case "SeTC": 
			mess = "For test-and-cull diagnostic test sensitivity"; 
			break;

		case "SpTC": 
			mess = "For test-and-cull diagnostic test specificity";
			break;
			
		case "mut_rate": 
			mess = "For a genetic mutation rate";
			break;
			
		case "seq_var": 
			mess = "For a genetic sequence variation";
			break;
			
		case "comp_prob": case "sim_comp_prob":
			mess = "For a compartmental probability";
			break;
		
		case "trans_prob":
			mess = "For a transition probability";
			break;
						
		case "reparam_ele":
			mess = "For reparameterisation of <e>"+get_full_parameter_name(info.par_name)+"</e>";
			break;
			
		case "define_eqn":
			mess = "For definition of <e>"+get_full_parameter_name(info.par_name)+"</e>";
			break;
			
		case "prior":
			mess = "For a distribution";
			break;
			
		default: 
			error("type not included: "+eqn.type); 
			break;
		}
			
		if(eqn.warn.length > 0){		
			add_warning({mess:mess, mess2:eqn.warn[0].te, eqn_info:eqn.eqn_info, eqn_type:eqn.type, warn_type:"Equation"});
		}
		
		for(let i = 0; i < eqn.param.length; i++){
			let par = copy(eqn.param[i]);
			if(eqn.type == "derive_param") par.derive = true;
			add_par_to_list(par,list,copy(eqn),mess);
		}
		
		for(let i = 0; i < eqn.ind_eff.length; i++){
			let ief = eqn.ind_eff[i]; 
			
			let p_name = eqn.p_name;
			
			let name = eqn.ind_eff[i].name;
			
			let j;
			for(j = 0; j < ind_eff_list.length; j++){
				let ie_list = ind_eff_list[j];
				if(ie_list.name == name){
					if(ie_list.p_name != p_name){
						add_warning({mess:mess, mess2:"The individual effect '"+name+"' appears in the transitions of two species", eqn_info:eqn.eqn_info, eqn_type:eqn.type, warn_type:"Equation"});
						break;
					}
				}
			}
			
			if(j == ind_eff_list.length){
				ind_eff_list.push({name:name, p_name:p_name, within:[]});
			}
			add_within(ind_eff_list[j],eqn);
		}
		
		for(let i = 0; i < eqn.fix_eff.length; i++){
			let fe = eqn.fix_eff[i];
			let name = fe.name;
			
			let p_name = eqn.p_name;
			
			let j;
			for(j = 0; j < fix_eff_list.length; j++){
				let fe_list = fix_eff_list[j];
				if(fe_list.name == name){
					if(fe_list.p_name != p_name){
						add_warning({mess:mess, mess2:"The fixed effect '"+name+"' appears in the transitions of two species", eqn_info:eqn.eqn_info, eqn_type:eqn.type, warn_type:"Equation"});
						break;
					}
				}
			}
			
			if(j == fix_eff_list.length){
				fix_eff_list.push({name:name, p_name:p_name, within:[]});
			}
			add_within(fix_eff_list[j],eqn);
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
			if(par.p_name != undefined){
				if(parlist.p_name == undefined) parlist.p_name = par.p_name;
				else{
					if(parlist.p_name != par.p_name){
						parlist.p_name = mixed_p_name;
					}
				}
			}
			
			check_if_param_same(par,parlist);
			break;		
		}
	}

	if(j == list.length){
		par.within=[];
		list.push(par);
	}
	
	if(eqn != undefined){
		add_within(list[j],eqn);
	}
}


/// Provides a description of the equation that the parameter appears in 
function param_eqn_desc(eqn)
{
	switch(eqn.type){
	case "Se": return "as a diagnostic test sensitivity";
	case "Sp": return "as a diagnostic test specificity";
	case "SeTC": return "as a test-and-cull test sensitivity";
	case "SpTC": return "as a test-and-cull test specificity";
	case "mut_rate": return "as a genetic mutation rate";
	case "seq_var": return "as a genetic sequence variation";
	case "trans_bp": return "as a branching probability";
	case "trans_mean": return "in a transition mean";
	case "trans_rate": return "in a transition rate";
	case "trans_shape": return "in a transition shape";
	case "trans_scale": return "in a transition scale";
	case "trans_cv": return "in a transition coefficient of variation";
	case "comp_prob": case "sim_comp_prob": return "in a compartmental probability";
	case "reparam_eqn": return "in a reparameterised equation";
	case "reparam_ele": return "in a reparametersation";
	case "define_eqn": return "in a definition equation";
	case "define_ele": return "in a definition";
	case "derive_param": return "in a derived parameter";
	case "derive_eqn": return "in a derived equation";
	default: error("option not pos:"+eqn.type); break;	
	}
}


/// Finds the tensor template for a parameter (this is a multi-dimensional tensor with unset elements)
function par_find_template(list,ndep_cont)
{
	let ndep = list.length;
	if(ndep_cont != undefined) ndep = ndep_cont;
	
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
			let ind = [];
			for(let i = 0; i < ndim; i++) ind.push(index[i]);
			ele_list.push(ind);
			
			let i = ndim-1; 
			do{
				flag = false;
				index[i]++; if(index[i] >= dim[i]){ index[i] = 0; i--; flag = true;}
			}while(i >= 0 && flag == true);
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
	if(Array.isArray(par.value)){
		for(let d = 0; d < par.list.length; d++) dim.push(par.list[d].length);
	}
	return dim;
}


/// Gets a description of value to send back to inferface
function get_value_desc(par)
{
	let te;
	
	let dim = [];
	for(let d = 0; d < par.ndep_cont; d++) dim.push(par.list[d].length);
		
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
	
	if(par.value != undefined){
		if(tot < 10){
			te = JSON.stringify(par.value);
			te = te.replace(/\"/g," ")
			te = te.replace(/null/g,"Unset");
		}
	}
	else te = "Unset "+te;
	
	return te;
}


/// Gets a description of prior_const to send back to inferface
function get_prior_const_desc(par)
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
	
	if(par.prior_const != undefined){
		if(tot < 10){
			te = JSON.stringify(par.prior_const);
			te = te.replace(/\"/g," ")
			te = te.replace(/null/g,"Unset");
		}
	}
	else te = "Unset "+te;
	
	return te;
}


/// Gets a description of weight factor to send back to inferface
function get_weight_desc(par)
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
	
	if(par.factor_weight != undefined){
		if(tot < 10){
			te = JSON.stringify(par.factor_weight);
			te = te.replace(/\"/g," ")
			te = te.replace(/null/g,"Unset");
		}
	}
	else te = "Unset " + te;
	
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
function find_comp_from_index(index,mod)
{
	if(mod == undefined) mod = model;
	
	let list = [];
	index = remove_prime(index);
	for(let p = 0; p < mod.species.length; p++){
		let sp = mod.species[p];
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


/// Finds the individual effect in a ind effect group from the parameter make
function find_ieg_list(par,mod)
{
	let list = [];
	
	if(mod == undefined) mod = model;
	
	let spl = par.name.split("^");
	if(spl.length != 2){ error("Spl should be 2"); return;}
	if(spl[0] != "Ω"){ error("Should be Ω"); return;}
	
	let name = spl[1];
	for(let p = 0; p < mod.species.length; p++){
		let sp = mod.species[p];
		for(let j = 0; j < sp.ind_eff_group.length; j++){
			let ieg = sp.ind_eff_group[j];
			if(ieg.name == name){
				for(let k = 0; k < ieg.ie_list.length; k++) list.push(ieg.ie_list[k].name);
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
function update_model_ind_eff(ie_list)
{
	for(let i = 0; i < ie_list.length; i++){
		let ie = ie_list[i];
		
		// Checks that all the species are correctly identified
		if(ie.p_name == mixed_p_name){
			//add_warning({mess:"Parameter issue", mess2:"Parameter features in multiple species and individual effect cannot be specified", eqn_info:eqn.eqn_info, eqn_type:eqn.type, warn_type:"Equation"});
			
			add_warning({mess:"Individual effect problem", mess2:"Cannot allocate the species for individual effect '["+ie.name+"]'", warn_type:"Equation"});
		}
	}
	
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		
		for(let i = 0; i < sp.ind_eff_group.length; i++){
			let ieg = sp.ind_eff_group[i];
			for(let j = 0; j < ieg.ie_list.length; j++) ieg.ie_list[j].used = false;
		}
	
		for(let k = 0; k < ie_list.length; k++){
			let ie = ie_list[k];
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
					sp.ind_eff_group.push({name:generate_iegroup_name(), p_name:ie.p_name, A_matrix:{check:false, loaded:false, value:[], ind_list:[], pedigree:false, sire_list:[], dam_list:[]}, ie_list:[{name:ie.name}]});
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

		if(sp.ind_eff_group.length > 0){
			switch(sp.type){
			case "Population": case "Deterministic":
				{
					let name = sp.ind_eff_group[0].ie_list[0].name;
					add_warning({mess:"Cannot have individual effects", mess2:"A population-based model cannot have individual effects such as '["+name+"]'", warn_type:"Equation"});
				}
				break;
			}
		}
	}
}


/// Generates a unique iegroup name
function generate_iegroup_name()
{
	let list=[];
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let i = 0; i < sp.ind_eff_group.length; i++){
			list.push(sp.ind_eff_group[i].name);
		}
	}

	let k = 1;
	while(true){
		let name = "CM"; if(k > 1) name += k;
		if(find_in(list,name) == undefined) return name;
		k++;
	}	
}


/// Updates fixed effects stored in model
function update_model_fix_eff(fe_list)
{
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		
		for(let i = 0; i < sp.fix_eff.length; i++){
			sp.fix_eff[i].used = false;
		}
	
		for(let k = 0; k < fe_list.length; k++){
			let fe = fe_list[k];
	
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
	
		if(sp.fix_eff.length > 0){
			switch(sp.type){
			case "Population": case "Deterministic":
				{
					let name = sp.fix_eff[0].name;
					add_warning({mess:"Cannot have fixed effects", mess2:"A population-based model cannot have fixed effects such as '〈"+name+"〉'", warn_type:"Equation"});
				}
				break;
			}
		}
	}
}


/// Adds individual effect parameters (i.e. variances) to model
function add_ind_eff_param(par_list,param)
{
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		
		for(let i = 0; i < sp.ind_eff_group.length; i++){
			let ieg = sp.ind_eff_group[i];
			
			let name = "Ω^"+ieg.name;
			
			if(find(param,"name",name) == undefined){
				if(ieg.ie_list.length == 1){
					let par = { name:name, time_dep:false, dep:[], within:[], variance:true};
					par.full_name = param_name(par);	
				
					add_par_to_list(par,par_list,{eqn_info:{p:p, ieg:i, eso:"model"}},"Individual effect variance");
				}
				else{
					let dep = ["z","z'"];
					let par = { name:name, time_dep:false, dep:dep, within:[], variance:true};
					par.full_name = param_name(par);			
					par_find_list(par);
					add_par_to_list(par,par_list,{eqn_info:{p:p, ieg:i, ie:ieg.ie_list[0].name, eso:"model"}},"Individual effect variance");
				}
			}
		}
	}
} 	


/// Adds individual effect parameters (i.e. variances) to model
function add_fix_eff_param(par_list,param)
{	
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
	
		for(let i = 0; i < sp.fix_eff.length; i++){
			let fe = sp.fix_eff[i];
			let name = fe_char+"^"+fe.name;
			if(find(param,"name",name) == undefined){
				let par = { name:name, time_dep: false, dep:[], within:[], fixed_effect:true};
				add_par_to_list(par,par_list,{eqn_info:{p:p, i:i, eso:"model"}},"Fixed effect parameter");
			}
		}
	}
}


/// Updates the model parameters
function update_model_param(lists)
{ 
	let par_list = lists.par_list;
	let ie_list = lists.ind_eff_list;       
  let fe_list = lists.fix_eff_list;      
	
	add_param_ppc_factor(par_list);                  // Adds any parameter for ppc factors
	
	if(model.warn.length > 0) err_warning();
	
	let param_old = model.param;
	
	// Copies across any definitions
	for(let th = 0; th < param_old.length; th++){
		let par = param_old[th];
		if(par.variety == "define"){
			if(find(par_list,"name",par.name) == undefined) par_list.push(par);
		}
	}
	
	let param=[];

	param_add_par_list(param,param_old,par_list,ie_list,fe_list);

	update_model_ind_eff(ie_list);                    // Updates individual effect groupings in model
	
	update_model_fix_eff(fe_list);                    // Updates fixed effects in model 

	{
		var var_par_list = [];
		add_ind_eff_param(var_par_list,param);          // Adds ind. effect parameters to model
		add_fix_eff_param(var_par_list,param);          // Adds fixed effect parameters to model
	
		param_add_par_list(param,param_old,var_par_list,ie_list,fe_list);
	}

	param.sort( function(a, b){                      // Alphabetically sorts parameters
		if(a.ndep_cont > b.ndep_cont) return 1; 
		if(a.ndep_cont == b.ndep_cont){
			if(a.name.toLowerCase() > b.name.toLowerCase()) return 1;
			else return -1;	
		}
		return -1;});
	
	if(model.warn.length > 0) err_warning();
	
	model.param = param;
	
	check_derived_param_all(model.param);
	check_reparam_all(model.param);
	check_define_all(model.param);
}


/// Adds the elements of par_list to param
function param_add_par_list(param,param_old,par_list,ie_list,fe_list)
{
	for(let i = 0; i < par_list.length; i++){
		let par = par_list[i];
		
		if(find(param,"name",par.name) != undefined) error("param should not exist"); 
		
		let j;                                         // Determines if parameter already exists
		for(j = 0; j < param_old.length; j++){
			let par_old = param_old[j];
			if(par.name == par_old.name){	
				if(par.dep.length == par_old.dep.length){
					let k = 0; 
					while(k < par.dep.length && par.dep[k] == par_old.dep[k]) k++;
					if(k == par.dep.length) break;
				}
				
				if(loading_from_file()){
					alert_line("The parameter definition '"+par_old.full_name+"' doesn't agree with that in the model '"+par.full_name+"'.",par_old.import_line);
				}
			}
		}
		
		if(j < param_old.length){                      // Copies information from existing parameter definition	
			param.push(copy_param_info(par,param_old[j]));
		
			add_dependent_param(param,par_list,ie_list,fe_list);
		}
		else{                                          // Create new parameter
			param.push(create_new_param(par,"normal"));
		}
	}
}


/// Goes through reparameterisation, definition and distribution and adds to parameter list
function add_dependent_param(param,par_list,ie_list,fe_list)
{
	let par = param[param.length-1];
	let p_name = par.p_name;

	let eqn_info = {par_name:par.name, eso:"dep"};
	
	let vari = par.variety;
	
	let wtype;
	let warn =[];
	
	switch(vari){
	case "reparam":
		if(par.reparam_eqn_on){
			wtype = "RepEqValue";
			add_ele_param_ie_fe(par.name,par.reparam_eqn,"reparam_eqn",par_list,ie_list,fe_list,eqn_info,warn,p_name);
		}
		else{
			wtype = "RepValue";
			for(let k = 0; k < par.defrep_warn.length; k++) warn.push(par.defrep_warn[k]);

			add_to_par_list(par.name,par.defrep_param_list,par_list);
			add_to_ie_list(par.defrep_ie_list,ie_list,p_name);
			add_to_fe_list(par.defrep_fe_list,fe_list,p_name);
		}

		if(warn.length > 0){
			add_warning({mess:"Reparameterisation error", mess2:"The reparameterisation of '"+par.full_name+"': "+warn[0], warn_type:wtype, name:par.name});
		}
		break;
	
	case "define":
		if(par.define_eqn_on){
			wtype = "DefEqValue";
			add_ele_param_ie_fe(par.name,par.define_eqn,"define_eqn",par_list,ie_list,fe_list,eqn_info,warn,p_name);
		}
		else{
			wtype = "DefValue";
			for(let k = 0; k < par.defrep_warn.length; k++) warn.push(par.defrep_warn[k]);
			
			add_to_par_list(par.name,par.defrep_param_list,par_list);
			add_to_ie_list(par.defrep_ie_list,ie_list,p_name);
			add_to_fe_list(par.defrep_fe_list,fe_list,p_name);
		}
		
		if(warn.length > 0){
			add_warning({mess:"Definition error", mess2:"The definition of '"+par.full_name+"': "+warn[0], warn_type:wtype, name:par.name});
		}
		break;
		
	case "dist":
		{
			if(par.ndep_cont == 0 || par.prior_split_check.check == false){	
				wtype = "DistValue";
				add_distribution_eqn(par.prior,vari,par_list,eqn_info,warn);	
			}
			else{
				wtype = "DistSplitValue";
				for(let k = 0; k < par.prior_warn.length; k++) warn.push(par.prior_warn[k]);
			
				add_to_par_list(par.name,par.prior_param_list,par_list);
			}
			
			if(warn.length > 0){
				add_warning({mess:"Distribution error", mess2:"The distribution of '"+par.full_name+"': "+warn[0], warn_type:wtype, name:par.name});
			}
		}
		break;
	}
	
	par.full_name = param_name(par);
}


/// Adds parameters defined in reparameterisation or priors
function add_to_par_list(par_name,add_param_list,par_list)
{
	for(let i = 0; i < add_param_list.length; i++){
		let par_add = add_param_list[i]; 

		let j = find(par_list,"name",par_add.name)
		
		if(j == undefined){
			par_list.push(par_add);
		}
		else{
			check_if_param_same(par_add,par_list[j]);

			for(let k = 0; k < par_add.within.length; k++){
				let wi = par_add.within[k];
				add_wi(par_list[j],wi);
			}
		}
		
		check_nested(par_name,par_add);
	}
}


/// Throws a warning if a parameter is within itself
function check_nested(par_name,par)
{
	if(par.name != par_name) return;
	
	let te = "The parameter "+par.full_name+" cannot be contained within itself.";
	
	if(par.within != undefined && par.within.length > 0){
		let wi = par.within[0];
		add_warning({mess:"Parameter problem", mess2:te, eqn_type: wi.type, eqn_info:wi.eqn_info, warn_type:"Equation"});
	}
	else{
		add_warning({mess:"Parameter problem", mess2:te, warn_type:"ParamPage"});
	}
}


/// Throws a warning if two parameters are different 
function check_if_param_same(par1,par2)
{
	if(par_same(par1,par2)) return;
	
	let line;
	if(par1.import_line != undefined) line = par1.import_line;
	if(par2.import_line != undefined) line = par2.import_line;

	let te = "The parameter "+par1.name+" has different dependencies (with "+par1.full_name+" in one place and "+par2.full_name+" in another).";
	
	if(par1.within != undefined && par1.within.length > 0){
		let wi = par1.within[0];
		add_warning({mess:"Model inconsistency", mess2:te, eqn_type:wi.type, eqn_info:wi.eqn_info, line:line, warn_type:"Equation"});
	}
	else{
		add_warning({mess:"Model inconsistency", mess2:te, line:line, warn_type:"ParamPage"});
	}
}
	

/// Adds individual effects onto list
function add_to_ie_list(add_ie_list,ie_list,p_name)
{
	for(let i = 0; i < add_ie_list.length; i++){
		let ie = add_ie_list[i];
		let j = find(ie_list,"name",ie.name);
		if(j == undefined){ j = ie_list.length; ie_list.push(ie);}
		
		if(ie_list[j].p_name == undefined) ie_list[j].p_name = p_name;
	}
}


/// Adds individual effects onto list
function add_to_fe_list(add_fe_list,fe_list,p_name)
{
	for(let i = 0; i < add_fe_list.length; i++){
		let fe = add_fe_list[i];
		let j = find(fe_list,"name",fe.name);
		if(j == undefined){ j = fe_list.length; fe_list.push(fe);}
		
		if(fe_list[j].p_name == undefined) fe_list[j].p_name = p_name;
	}
}
	

/// Adds any parameters from a distribution to the list of equations
function add_distribution_eqn(dist,vari,par_list,eqn_info,warn)
{
	let val = dist.value;
	
	switch(dist.type.te){
	case "fix": 
		add_ele_param(val.mean_eqn.te,vari,par_list,eqn_info,warn);
		break;
		
	case "uniform":
		add_ele_param(val.min_eqn.te,vari,par_list,eqn_info,warn);
		add_ele_param(val.max_eqn.te,vari,par_list,eqn_info,warn);
		break;
		
	case "covar-normal-lkj": 
		add_ele_param(val.sd_eqn.te,vari,par_list,eqn_info,warn);
		add_ele_param(val.eta_eqn.te,vari,par_list,eqn_info,warn);
		break;
		
	case "covar-uniform-lkj": 
		add_ele_param(val.min_eqn.te,vari,par_list,eqn_info,warn);
		add_ele_param(val.max_eqn.te,vari,par_list,eqn_info,warn);
		add_ele_param(val.eta_eqn.te,vari,par_list,eqn_info,warn);
		break;
		
	case "covar-inv-wishart":
		add_ele_param(val.S_eqn.te,vari,par_list,eqn_info,warn);
		add_ele_param(val.nu_eqn.te,vari,par_list,eqn_info,warn);
		break;
		
	case "covar-jeffreys": case "covar-uniform":
		add_ele_param(val.min_eqn.te,vari,par_list,eqn_info,warn);
		add_ele_param(val.max_eqn.te,vari,par_list,eqn_info,warn);
		break;
		
	case "inverse":
		add_ele_param(val.min_eqn.te,vari,par_list,eqn_info,warn);
		add_ele_param(val.max_eqn.te,vari,par_list,eqn_info,warn);
		break;
		
	case "power":
		add_ele_param(val.min_eqn.te,vari,par_list,eqn_info,warn);
		add_ele_param(val.max_eqn.te,vari,par_list,eqn_info,warn);
		add_ele_param(val.power_eqn.te,vari,par_list,eqn_info,warn);
		break;
		
	case "exp": 
		add_ele_param(val.mean_eqn.te,vari,par_list,eqn_info,warn);
		break;
		
	case "normal": 
		add_ele_param(val.mean_eqn.te,vari,par_list,eqn_info,warn);
		add_ele_param(val.sd_eqn.te,vari,par_list,eqn_info,warn);
		break;
	
	case "log-normal": 
		add_ele_param(val.mean_eqn.te,vari,par_list,eqn_info,warn);
		add_ele_param(val.cv_eqn.te,vari,par_list,eqn_info,warn);
		break;
		
	case "gamma":
		add_ele_param(val.mean_eqn.te,vari,par_list,eqn_info,warn);
		add_ele_param(val.cv_eqn.te,vari,par_list,eqn_info,warn);
		break;
		
	case "bernoulli":
		add_ele_param(val.mean_eqn.te,vari,par_list,eqn_info,warn);
		break;
					
	case "beta":
		add_ele_param(val.alpha_eqn.te,vari,par_list,eqn_info,warn);
		add_ele_param(val.beta_eqn.te,vari,par_list,eqn_info,warn);
		break;
	
	case "mdir":
		add_ele_param(val.sigma_eqn.te,vari,par_list,eqn_info,warn);
		break;
		
	case select_str:
		break;
		
	default: error("not op:"+dist.type.te); break;
	}
}


/// Gets the type of equation fro 
function get_type_from_vari(vari,time_dep)
{
	switch(vari){
	case "define": 
		if(time_dep == true) return "define_ele";
		return "define_ele_notime";
		
	case "reparam": 
		if(time_dep == true) return "reparam_ele";
		return "reparam_ele_notime";
	
	case "dist": 
		return "dist";
	
	default: error("Which ty from vari:"+vari); break;
	}
}

		
/// Adds any parameters from a text element to the parameter list
function add_ele_param(ele,vari,par_list,eqn_info,warn,p_name)
{
	if(isNaN(ele)){
		let p;
		if(p_name != undefined && p_name != mixed_p_name){
			p = find(model.species,"name",p_name);
			if(p == undefined) error("Still can't find species");
		}
		
		let eqn = create_equation(ele,get_type_from_vari(vari,false),p);
		eqn.eqn_info = eqn_info;
		
		if(eqn.warn.length > 0){
			warn.push(eqn.warn[0].te);
		}
	
		for(let i = 0; i < eqn.param.length; i++){
			let par = eqn.param[i];
			par.p_name = p_name;
			add_par_to_list(par,par_list,copy(eqn),"Reparameterisation");
		}
	}
}


/// Adds any parameters, ie or fe from a text element to lists
function add_ele_param_ie_fe(par_name,ele,eq_type,par_list,ie_list,fe_list,eqn_info,warn,p_name)
{
	if(isNaN(ele)){
		let p;
		if(p_name != undefined && p_name != mixed_p_name){
			p = find(model.species,"name",p_name);
			if(p == undefined) error("Still can't find species");
		}
	
		let eqn = create_equation(ele,eq_type,p);
		eqn.eqn_info = eqn_info;
		
		if(eqn.warn.length > 0){
			warn.push(eqn.warn[0].te);
		}
		
		for(let i = 0; i < eqn.ind_eff.length; i++){
			let ie = eqn.ind_eff[i];
			if(find(ie_list,"name",ie.name) == undefined) ie_list.push(ie);
		}
		
		for(let i = 0; i < eqn.fix_eff.length; i++){
			let fe = eqn.fix_eff[i];
			if(find(fe_list,"name",fe.name) == undefined) fe_list.push(fe);
		}
	
		for(let i = 0; i < eqn.param.length; i++){
			let par = eqn.param[i];	
			par.p_name = p_name;
			add_par_to_list(par,par_list,copy(eqn),"Reparameterisation");
			check_nested(par_name,par);	
		}
	}
}


/// Gets a list of all parameters in a reparameterisation
function get_defrep_param_list(par) 
{
	par.defrep_param_list=[];
	par.defrep_ie_list=[];
	par.defrep_fe_list=[];
	par.defrep_warn=[];
	if((par.variety == "reparam" && par.reparam_eqn_on != true) || 
		 (par.variety == "define" && par.define_eqn_on != true)){
		let eq_type = get_type_from_vari(par.variety,par.time_dep);
			 
		let eqn_info = {type:eq_type, par_name:par.name, eso:"rep"};
		let par_list = [];
		let ie_list = [];
		let fe_list = [];
		
		let ndep = par.ndep_cont;
		
		if(ndep == 0) error("Should not be here");
		
		let co_list = generate_co_list(par.list,ndep);
		for(let i = 0; i < co_list.length; i++){
			let comb = co_list[i];
			let ele = get_element(par.value,comb.index);
			eqn_info.index = comb.index;	
			add_ele_param_ie_fe(par.name,ele,eq_type,par_list,ie_list,fe_list,eqn_info,par.defrep_warn);
		}
	
		par.defrep_param_list = par_list;
		
		par.defrep_ie_list = ie_list;
		par.defrep_fe_list = fe_list;
	}
}


/// Gets a list of all parameters in prior
function get_prior_param_list(par) 
{
	par.prior_param_list=[];
	par.prior_warn=[];
	if(par.ndep_cont > 0 && par.prior_split_check.check == true){
		let eqn_info = {par_name:par.name, eso:"pri"};
		let par_list = [];
		let vari = par.variety;
		
		let co_list = generate_co_list(par.list);
		for(let i = 0; i < co_list.length; i++){
			let comb = co_list[i];
			let ele = get_element(par.prior_split,comb.index);
				
			add_distribution_eqn(ele,vari,par_list,eqn_info,par.prior_warn);	
		}

		par.prior_param_list = par_list;
	}
}


/// Creates a new parameter
// par must have the properties: name, dep, type, time_dep, (optionally within)
// On top of these new properties are added: 
// list: shows the potential values for the dependencies
function create_new_param(par,variety)
{
	if(par.name == undefined) error("Parameter name must be set");
	if(par.dep == undefined) error("Parameter dependency must be set");
	if(par.time_dep == undefined) error("Parameter time dependency must be set");
	
	par.variety = variety;
	
	if(par.variance == undefined) par.variance = false;
	if(par.derive == undefined) par.derive = false;
	if(par.param_fac == undefined) par.param_fac = false;
	if(par.fixed_effect == undefined) par.fixed_effect = false;
	
	par.pri_pos = set_pri_pos(par.variance,false);

	set_ndep_cont(par);
	
	let knot = ["start","end"];                      // Initialises a spline
	
	let on = false;
	if(par.time_dep == true) on = true;
	
	par.spline = { on:on, smooth:default_smooth(), spline_radio:spline_radio_pos[0], time_dep:par.time_dep, knot:knot}; 
		
	par.sim_sample = {check:true};
		 
	par.prior = unset_prior();
	
	par.reparam_eqn_on = false; if(par.ndep_cont == 0) par.reparam_eqn_on = true;
	par.reparam_eqn = "";
	
	par.define_eqn_on = false; if(par.ndep_cont == 0) par.define_eqn_on = true;
	par.define_eqn = "";
	
	par_set_default(par);

	par.full_name = param_name(par);

	if(par.variance == true) par.prior.type.te = "covar-default";
	
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
	par.dist_mat = false;
	par.iden_mat = false;
	par.den_vec = false;
	
	if(par.ndep_cont == 0){
		par.value = set_str;
		par.prior_const = set_str;
		if(par.list){
			delete par.list;
			delete par.prior_split; delete par.prior_split_set; delete par.prior_split_check;
			delete par.factor_weight;
			delete par.set;
			delete par.prior_const_set;
		}
	}
	else{
		if(is_density_name(par.name)){ density_set_default(par); return;}
		
		if(par.name == dist_matrix_name){ dist_set_default(par); return;}
		
		if(par.name == iden_matrix_name || par.name == iden_matrix_name2){ iden_set_default(par); return;}
	
		if((par.variety == "reparam" && par.reparam_eqn_on) || 
			(par.variety == "define" && par.define_eqn_on)){
			reparam_eqn_set_default(par); return;
		}
	
		par.list = par_find_list(par);
		par.prior_split_check = {check:false};
		par.set = false;
		par.prior_const_set = false;
	}
	par.prior_const_on = false;
	par.defrep_param_list = [];
	par.defrep_ie_list=[];
	par.defrep_fe_list=[];
	par.defrep_warn = [];
	par.prior_param_list = [];
	par.prior_warn = [];
	par.factor = false;
	par.factor_weight_on = {check:false};
}


/// Sets default values for distance matrix
function dist_set_default(par)
{
	par.list = par_find_list(par);
	if(par.prior_split) delete par.prior_split;
	if(par.value) delete par.value;
	par.prior_split_check = {check:false};
	par.variety = "const";	
	par.dist_mat = true;
	par.iden_mat = false;
	par.den_vec = false;
	par.prior_const_set = false;
	par.set = false;
}


/// Sets default values for distance matrix
function density_set_default(par)
{
	par.list = par_find_list(par);
	if(par.prior_split) delete par.prior_split;
	if(par.value) delete par.value;
	par.prior_split_check = {check:false};
	par.variety = "const";	
	par.den_vec = true;
	par.dist_mat = false;
	par.iden_mat = false;
	par.prior_const_set = false;
	par.set = false;
}


/// Sets default values for identity matrix
function iden_set_default(par)
{
	par.list = par_find_list(par);
	//par.value = [];
	//par.prior_split = [];
	if(par.prior_split) delete par.prior_split;
	if(par.value) delete par.value;
	par.prior_split_check = {check:false};
	par.variety = "const";
	par.dist_mat = false;	
	par.iden_mat = true;
	par.den_vec = false;
	par.prior_const_set = false;
	par.set = false;
}


/// Sets default values for distance matrix
function reparam_eqn_set_default(par)
{
	par.list = par_find_list(par);
	if(par.prior_split) delete par.prior_split;
	if(par.value) delete par.value;
	par.prior_split_check = {check:false};
	par.set = false;
}


/// Sets default values for tensor
function set_default_value(par)
{
	if(par.dist_mat || par.iden_mat || par.den_vec) return;
	
	let value = par.value;
	let dim = get_dimensions(value);
	let ele_list = get_element_list(value,dim);
		
	for(let k = 0; k < ele_list.length; k++){	
		set_element(value,ele_list[k],0);
	}
}


/// Sets default values for tensor
function set_default_factor_weight(par)
{
	par.factor_weight = param_blank(par);
		
	let value = par.factor_weight;
	let dim = get_dimensions(value);
	let ele_list = get_element_list(value,dim);
		
	for(let k = 0; k < ele_list.length; k++){	
		set_element(value,ele_list[k],1);
	}
}


/// Copies infomation from old parameter definition to new 		
function copy_param_info(par,old)
{
	par.variety = old.variety;
	
	/*
	par.variance = old.variance;
	par.derive = old.derive;
	par.param_fac = old.param_fac;
	par.fixed_effect = old.fixed_effect;
	*/
	set_ndep_cont(par);
	
	par.spline = old.spline;
	par.prior = old.prior;
	par.sim_sample = old.sim_sample;
	par.label_info = old.label_info;
	par.set = old.set;
	par.prior_const_set = old.prior_const_set;
	par.prior_const_on = old.prior_const_on;
	par.import_line = old.import_line;
	par.defrep_param_list = old.defrep_param_list;
	par.defrep_ie_list = old.defrep_ie_list;
	par.defrep_fe_list = old.defrep_fe_list;
	par.defrep_warn = old.defrep_warn;
	par.prior_warn = old.prior_warn;
	par.prior_param_list = old.prior_param_list;
	par.reparam_eqn_on = old.reparam_eqn_on;
	par.reparam_eqn = old.reparam_eqn;
	par.define_eqn_on = old.define_eqn_on;
	par.define_eqn = old.define_eqn;
	par.dist_mat = false;
	par.iden_mat = false;
	par.den_vec = false;
	par.factor = old.factor;
	par.factor_weight_on = old.factor_weight_on;
	
	par.pri_pos = set_pri_pos(par.variance,par.factor);
	
	if(is_density_name(par.name)){ density_set_default(par); return par;}
		
	if(par.name == dist_matrix_name){ dist_set_default(par); return par;}
	
	if(par.name == iden_matrix_name || par.name == iden_matrix_name2){ iden_set_default(par); return par;}
	
	if((par.variety == "reparam" && par.reparam_eqn_on) || 
			(par.variety == "define" && par.define_eqn_on)){
		reparam_eqn_set_default(par); 
	}
	
	if(par.ndep_cont == 0){	
		par.value = old.value;
		par.prior_const = old.prior_const;
	}
	else{
		let list = par_find_list(par);
	
		let list_old = old.list;
	
		// Checks if list is the same as the old version
		let dif = false;
		
		let val_set = false; if(old.value != undefined) val_set = true;
		let prc_set = false; if(old.prior_const != undefined) prc_set = true;
		let pri_set = false; if(old.prior_split != undefined) pri_set = true;
		let wei_set = false; if(old.factor_weight != undefined) wei_set = true;
		
		if(val_set || prc_set || pri_set || wei_set){
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
				let co_list = generate_co_list(list);
			
				let value; if(val_set) value = copy(temp);
				let prior_const; if(prc_set) prior_const = copy(temp);
				let prior_split; if(pri_set) prior_split = copy(temp);
				let weight; if(wei_set) weight = copy(temp);
			
				//let prior = copy(old.prior);
				
				let old_fl = false;
				
				// Attempts to copy across existing values
				if(dif == true){
					let map=[];
					for(let k = 0; k < par.ndep_cont; k++){
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
						
						let val, pr_co, pr_sp, wei;
						
						let ind_old=[];
						let ii;
						for(ii = 0; ii < par.ndep_cont; ii++){
							ind_old[ii] = map[ii][ele_list[k][ii]];
							if(ind_old[ii] == undefined) break;
						}
						
						if(ii == par.ndep_cont){
							old_fl = true;
							if(val_set) val = get_element(old.value,ind_old);
							if(prc_set) pr_co = get_element(old.prior_const,ind_old);
							if(pri_set) pr_sp = get_element(old.prior_split,ind_old);
							if(wei_set) wei = get_element(old.factor_weight,ind_old);
						}
						else{
							if(val_set) val = undefined;
							if(prc_set) pr_co = undefined;
							if(pri_set) pr_sp = unset_prior();
							if(wei_set) wei = 1;
						}
						
						if(val_set) set_element(value,ind,val);
						if(prc_set) set_element(prior_const,ind,pr_co);
						if(pri_set) set_element(prior_split,ind,pr_sp);
						if(wei_set) set_element(weight,ind,wei);
					}
				}
				
				par.list = list;
				if(val_set) par.value = value;
				if(prc_set) par.prior_const = prior_const;
				if(pri_set) par.prior_split = prior_split;
				if(wei_set) par.factor_weight = weight;
				if(old_fl == false){
					par.set = false;
					par.prior_const_set = false;
				}
			}
			else{
				par.list = old.list;
				if(val_set) par.value = old.value;
				if(prc_set) par.prior_const = old.prior_const;
				if(pri_set) par.prior_split = old.prior_split;
				if(wei_set) par.factor_weight = old.factor_weight;
			}
		}
		else{
			par.list = list;
		}
			
		par.prior_split_set = old.prior_split_set;
		par.prior_split_check = old.prior_split_check;
	}
	
	return par;
}


/// Determines if the parameter is a density
function is_density_name(name)
{
	let spl = name.split("^");
	if(spl[0] == density_name || spl[0] == rdensity_name) return true;
	return false;
}


/// Finds the possible set of priors for a parameter type
function set_pri_pos(variance,factor)
{
	let pri_pos = prior_pos;  
	if(factor) pri_pos = prior_factor_pos;
	if(variance == true) pri_pos = prior_cv_pos;
	
	return pri_pos;
}


/// Given a list of dependencies this provides a list of cominations going through each
function generate_co_list(list,ndep_cont)
{
	let ndep = list.length;
	if(ndep_cont != undefined) ndep = ndep_cont;
	
	let co_list = [];
	
	let list_len = [];
	for(let d = 0; d < ndep; d++){
		let k = 0;
		while(k < list[d].length && list[d][k] != undefined) k++;
		list_len[d] = k;
	}
	
	if(ndep == 0) return [{index:[]}];
	
	for(let d = 0; d < ndep; d++){
		if(list_len[d] == 0){
			return co_list;
		}
	}
	
	let index=[];
	for(let d = 0; d < ndep; d++) index[d] = 0;

	let flag;
	do{
		co_list.push({index:copy(index)});
		
		let i = ndep-1;
		do{
			flag = false;
			
			index[i]++; 
			if(index[i] >= list_len[i]){
				index[i] = 0; 
				i--; 
				flag = true;
			}
		}while(i >= 0 && flag == true);
	}while(flag == false);
	
	return co_list;
}


/// Sets an unspecified prior
function unset_prior()
{
	var eb = {te:"", type:"prior", mode:"param only"};
	var value = {min_eqn:copy(eb), max_eqn:copy(eb), mean_eqn:copy(eb), power_eqn:copy(eb), shape_eqn:copy(eb), cv_eqn:copy(eb), sd_eqn:copy(eb), alpha_eqn:copy(eb), beta_eqn:copy(eb), sigma_eqn:copy(eb), eta_eqn:{te:"1", type:"prior",  mode:"param only"}, S_eqn:{te:"1", type:"prior", mode:"param only"}, nu_eqn:copy(eb)};
	
	let te = select_str;
	
	value.sigma_eqn.te = "0.5";
				
	return {type:{te:te}, value:value};
}


/// Checks all the derived parameters
function check_derived_param_all(param)
{
	for(let i = 0; i < model.derive.length; i++){
		let der = model.derive[i];
		let eqn1 = der.eqn1;
		let eqn2 = der.eqn2;
		
		let res = check_derived_param(eqn2,eqn1,param);
		
		if(res.err == true){
			add_warning({mess:"Error in derived expression", mess2:res.msg, warn_type:"Equation", eqn_info:{i:i}, eqn_type:eqn2.type, line:der.line});
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
		
		let j = find(param,"name",par.name);
	
		if(j == undefined){
			return err("Parameter '"+par.full_name+"' is not found in the model");
		}
		else{
			if(par_same(par,param[j]) == false){
				let par_exist = param[j];
				return err("Parameter '"+par.full_name+"' should have the dependency '"+par_exist.full_name+"' from the model");
			}
		}
	}
	
	// Checks that the dependencies for a derived quantity agree with equation
	if(eqn_param == undefined) return err("eqn_param Should be defined");
	
	extract_equation_properties(eqn_param);

	if(eqn_param.param.length == 1){
		let par_param = eqn_param.param[0];
		for(let i = 0; i < param.length; i++){
			let par = param[i];
			if(par.derive != true){
				if(par_param.name == param[i].name){
					return err("The model already contains a parameter with name '"+par_param.name+"'");
				}
			}
		}

		let dep = copy(eqn.dep);
		dep = dep.sort();
		
		let dep_param = copy(eqn_param.dep);
		dep_param = dep_param.sort();
		
		let par = eqn_param.param[0];
		let par2 = {name:par.name, dep:dep};
		let full_name = param_name(par2);
			
		if(dep.length != dep_param.length){
			return err("Problem with dependency - Expected '"+full_name+"' not '"+par.full_name+"'.");
		}
		else{
			for(let i = 0; i < dep.length; i++){
				if(remove_prime(dep[i]) != remove_prime(dep_param[i])){
					return err("Problem with dependency - Expected '"+full_name+"' not '"+par.full_name+"'.");
				}
			}
			
			for(let i = 0; i < dep.length; i++){
				if(dep[i] != dep_param[i]){
					return err("Problem with dependency - Primes do not match up. Expected '"+full_name+"' not '"+par.full_name+"'.");
				}
			}
		}
	}
	
	return success();
}


/// Checks all the reparameterised parameters
function check_reparam_all(param)
{
	for(let th = 0; th < param.length; th++){
		let par = param[th];
		if(par.variety == "reparam" && par.reparam_eqn_on){
			let res = define_reparam_eqn_check(par.reparam_eqn,th);
			
			if(res.err == true){
				add_warning({mess:"Error in reparameterised expression", mess2:res.msg, warn_type:"Equation", eqn_info:{i:th}, eqn_type:"reparam_eqn", line:par.import_line});
			}
		}
	}
}


/// Checks all the defined parameters
function check_define_all(param)
{
	for(let th = 0; th < param.length; th++){
		let par = param[th];
		if(par.variety == "define" && par.define_eqn_on){
			let res = define_reparam_eqn_check(par.define_eqn,th);
			
			if(res.err == true){
				add_warning({mess:"Error in defined expression", mess2:res.msg, warn_type:"Equation", eqn_info:{i:th}, eqn_type:"reparam_eqn", line:par.import_line});
			}
		}
	}
}


/// Checks reparameterisation is correct
function define_reparam_eqn_check(te,th)
{
	let eqn = create_equation(te,"reparam_eqn");
	
	let par = model.param[th];
	
	// Makes sure that all dependencies in equation are in the parameter
	let fl = false;
	for(let d = 0; d < eqn.dep.length; d++){
		let de = eqn.dep[d];
		if(find_in(par.dep,de) == undefined){
			//if(!(par.variety == "define" && de == "t")){
			let te = "Problem with dependency - The expression depends on ";
			if(de == "t") te += "time";
			else te += "index "+de;
			te += ", but this isn't contained in parameter "+par.name;
			return err(te);
			//}
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
			//if(equal_ob(par.co_list,generate_co_list(par.list)) == false) error("Co list not correct");
		}
		
		if(par.dep.length == 0){
			if(par.list != undefined){ error("List should not be specified");}
			//if(par.co_list != undefined) error("co_list should not be specified");
			if(par.prior_split != undefined) error("prior_split should not be specified");
			if(par.prior_split_check != undefined) error("prior_split_check should not be specified");
			if(par.factor_weight != undefined) error("factor_weight should not be specified");
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


/// If any parameters are called "D" then put values for the distance matrix into them
function set_dist(info,par)
{
	let claa = get_cla_from_index(model,par.dep[0]);
	if(claa == undefined){ error("Classification is undefined"); return;}
	
	let N = claa.comp.length;

	// Works out if we need to calculate the full matrix
	if(N*N > ELEMENT_MAX){
		N = Math.floor(ELE_REDUCE_FAC*Math.sqrt(ELEMENT_MAX));
		info.too_big = true;
	}
	
	let li=[];
	for(let i = 0; i < N; i++){
		li.push(claa.comp[i].name);
	}
	
	info.list[0] = li;
	info.list[1] = li;
	
	let value = [];
	for(let j = 0; j < N; j++) value[j]=[];
				
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
						value[j][i] =  precision(Math.sqrt(dx*dx + dy*dy),5);
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
		
			for(let j = 0; j < N; j++){
				for(let i = j; i < N; i++){
					if(i == j) value[j][i] = 0;
					else{
						value[j][i] = precision(geo_dist(lat[i],lng[i],lat[j],lng[j]),5);
					}
					
					value[i][j] = value[j][i];
				}
			}
		}
		break;
	}
	
	info.value = value;
	let shrunk=[];
	for(let d = 0; d < par.ndep_cont; d++){
		if(info.list.length < par.list[d].length) shrunk[d] = true;
		else shrunk[d] = false;
	}
	info.shrunk = shrunk;
}


/// Calculates geographical distance	
function geo_dist(lat1,lng1,lat2,lng2)
{
	let r = 6371;
	let si = Math.sin(0.5*(lat2-lat1));
	let si2 = Math.sin(0.5*(lng2-lng1))
	return 2*r*Math.asin(Math.sqrt(si*si + Math.cos(lat1)*Math.cos(lat2)*si2*si2));
}


/// If any parameters are called "D" then put values for the distance matrix into them
function set_iden(info,par)
{
	let claa = get_cla_from_index(model,par.dep[0]);
	if(claa == undefined){ error("Classification is undefined"); return;}
	
	let N = claa.comp.length;

	// Works out if we need to calculate the full matrix
	if(N*N > ELEMENT_MAX){
		N = Math.floor(ELE_REDUCE_FAC*Math.sqrt(ELEMENT_MAX));
		info.too_big = true;
	}
	
	let li=[];
	for(let i = 0; i < N; i++){
		li.push(claa.comp[i].name);
	}
	
	info.list[0] = li;
	info.list[1] = li;
	
	let value = [];
	for(let j = 0; j < N; j++){
		value[j]=[];
		for(let i = 0; i < N; i++){
			if(i == j) value[j][i] = 1;
			else value[j][i] = 0;
		}
	}
	
	info.value = value;
	let shrunk=[];
	for(let d = 0; d < par.dep.length; d++){
		if(info.list.length < par.list[d].length) shrunk[d] = true;
		else shrunk[d] = false;
	}
	info.shrunk = shrunk;
}
	

/// If any parameters are called "DEN" then put values for the density vector into them
function set_density(info,par,limit_on)
{
	let spl = par.name.split("^");
	let kernel_r = Number(spl[1]);
	
	let claa = get_cla_from_index(model,par.dep[0]);
	if(claa == undefined){ error("Classification is undefined"); return;}
	
	let N = claa.comp.length;
	let N_tot = N;
	
	/*
	// Works out if we need to calculate the full matrix
	if(N > ELEMENT_MAX && limit_on != false){
		N = ELEMENT_MAX;
		info.too_big = true;
	}
	*/
	
	let li=[];
	for(let i = 0; i < N; i++){
		li.push(claa.comp[i].name);
	}
	
	info.list=[];
	info.list[0] = li;
	
	let value = [];
				
	switch(claa.camera.coord){
	case "cartesian":
		{
			let px=[], py=[];
	
			let xmin = LARGE, xmax = -LARGE, ymin = LARGE, ymax = -LARGE;
		
			for(let i = 0; i < N_tot; i++){
				let x = claa.comp[i].x;
				let y = claa.comp[i].y;
				
				px.push(x);
				py.push(y);
				
				if(x > xmax) xmax = x;
				if(x < xmin) xmin = x;
				if(y > ymax) ymax = y;
				if(y < ymin) ymin = y;
			}
		
			let xdist = xmax-xmin;
			let ydist = ymax-ymin;
			
			let LX = Math.floor(xdist/kernel_r); if(LX > 100) LX = 100;
			let LY = Math.floor(ydist/kernel_r); if(LY > 100) LY = 100;
			
			let dx = xdist/(LX*ALMOST_ONE);
			let dy = ydist/(LY*ALMOST_ONE);
			
			let grid=[];
			for(let j = 0; j < LY; j++){
				grid[j]=[];
				for(let i = 0; i < LX; i++){
					grid[j][i]=[];
				}
			}
			
			for(let k = 0; k < N_tot; k++){
				let i = Math.floor((px[k]-xmin)/dx); if(i < 0 || i >= LX) error("Out of range");
				let j = Math.floor((py[k]-ymin)/dy); if(j < 0 || j >= LY) error("Out of range");
				grid[j][i].push(k);
			}

			for(let k = 0; k < N; k++){
				let imid = Math.floor((px[k]-xmin)/dx); if(imid < 0 || imid >= LX) error("Out of range");
				let jmid = Math.floor((py[k]-ymin)/dy); if(jmid < 0 || jmid >= LY) error("Out of range");
		
				let sum = 0.0;
				for(let j = jmid-density_kernel_max; j <= jmid+density_kernel_max; j++){
					if(j >= 0 && j < LY){
						for(let i = imid-density_kernel_max; i <= imid+density_kernel_max; i++){
							if(i >= 0 && i < LX){
								for(let m = 0; m < grid[j][i].length; m++){
									let kk = grid[j][i][m];
									let ddx = px[kk]-px[k];
									let ddy = py[kk]-py[k];
									let dd = ddx*ddx + ddy*ddy;
						
									let dd_sc = dd/(kernel_r*kernel_r);
									if(dd_sc < density_kernel_max*density_kernel_max){
										sum += Math.exp(-0.5*dd_sc);
									}
								}
							}
						}		
					}
				}
				value[k] = sum/(2*Math.PI*kernel_r*kernel_r);
			
				if(false){  // Used for checking
					let sum_ch = 0.0;
					for(let i = 0; i < N_tot; i++){	
						let ddx = px[i]-px[k];
						let ddy = py[i]-py[k];
						let dd = ddx*ddx + ddy*ddy;
			
						let dd_sc = dd/(kernel_r*kernel_r);
						if(dd_sc < density_kernel_max*density_kernel_max){
							sum_ch += Math.exp(-0.5*dd_sc);
						}
					}
					
					if(dif(sum,sum_ch)) error("Sum not agree");
				}
			}
		}
		break;
		
	case "latlng":
		{
			let lat=[], lng=[];
			
			let lat_min = LARGE, lat_max = -LARGE, lng_min = LARGE, lng_max = -LARGE;
			
			for(let i = 0; i < N_tot; i++){
				let co = claa.comp[i];
				
				let x, y;
				if(co.type == "boundary"){ x = co.xmid; y = co.ymid;}
				else{ x = co.x; y = co.y;}
				
				let p = transform_latlng_inv(x,y);
				
				let la = p.lat*Math.PI/180;
				let ln = p.lng*Math.PI/180;
				
				lat.push(la);
				lng.push(ln);
				
				if(la > lat_max) lat_max = la;
				if(la < lat_min) lat_min = la;
				if(ln > lng_max) lng_max = ln;
				if(ln < lng_min) lng_min = ln;
			}
		
			let lat_far = lat_max; if(lat_min*lat_min > lat_max*lat_max) lat_far = lat_min;
			
			let lng_dist = geo_dist(lat_far,lng_min,lat_far,lng_max);
			let lat_dist = geo_dist(lat_min,lng_max,lat_max,lng_max);
			
			let LX = Math.floor(lng_dist/kernel_r); if(LX > 100) LX = 100;
			let LY = Math.floor(lat_dist/kernel_r); if(LY > 100) LY = 100;
			
			let dlng = (lng_max-lng_min)/(LX*ALMOST_ONE);
			let dlat = (lat_max-lat_min)/(LY*ALMOST_ONE);
			
			let grid=[];
			for(let j = 0; j < LY; j++){
				grid[j]=[];
				for(let i = 0; i < LX; i++){
					grid[j][i]=[];
				}
			}
			
			for(let k = 0; k < N_tot; k++){
				let i = Math.floor((lng[k]-lng_min)/dlng); if(i < 0 || i >= LX) error("Out of range");
				let j = Math.floor((lat[k]-lat_min)/dlat); if(j < 0 || j >= LY) error("Out of range");
				grid[j][i].push(k);
			}

			for(let k = 0; k < N; k++){
				let imid = Math.floor((lng[k]-lng_min)/dlng); if(imid < 0 || imid >= LX) error("Out of range");
				let jmid = Math.floor((lat[k]-lat_min)/dlat); if(jmid < 0 || jmid >= LY) error("Out of range");
		
				let sum = 0.0;
				for(let j = jmid-density_kernel_max; j <= jmid+density_kernel_max; j++){
					if(j >= 0 && j < LY){
						for(let i = imid-density_kernel_max; i <= imid+density_kernel_max; i++){
							if(i >= 0 && i < LX){
								for(let m = 0; m < grid[j][i].length; m++){
									let kk = grid[j][i][m];
									let d = geo_dist(lat[kk],lng[kk],lat[k],lng[k]);
						
									let d_sc = d/kernel_r;
									if(d_sc < density_kernel_max){
										sum += Math.exp(-0.5*d_sc*d_sc);
									}
								}
							}
						}		
					}
				}
				value[k] = sum/(2*Math.PI*kernel_r*kernel_r);
			
				if(false){  // Used for checking
					let sum_ch = 0.0;
					for(let i = 0; i < N_tot; i++){
						let d = geo_dist(lat[i],lng[i],lat[k],lng[k]);
						
						let d_sc = d/kernel_r;
						if(d_sc < density_kernel_max){
							sum_ch += Math.exp(-0.5*d_sc*d_sc);
						}
					}
					
					if(dif(sum,sum_ch)) error("Sum not agree");
				}
			}
		}
		break;
	}
	
	if(relative_den(par.name)){ // Converts to the relative density 
		let sum = 0.0;
		for(let k = 0; k < N; k++) sum += value[k];
		sum /= N;
		for(let k = 0; k < N; k++) value[k] /= sum;
	}
	
	// Sets to a certain precision
	for(let k = 0; k < N; k++) value[k] = precision(value[k],5);
		
	info.value = value;
	let shrunk=[];
	for(let d = 0; d < par.ndep_cont; d++){
		if(info.list.length < par.list[d].length) shrunk[d] = true;
		else shrunk[d] = false;
	}
	info.shrunk = shrunk;
}


/// Finds a list of all the equations in the model
function find_equation_list()
{
	let eqn_list = [];
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			for(let i = 0; i < claa.ntra; i++){
				let eqn_info = {p:p, cl:cl, i:i, eso:"model"};
				
				let traa = claa.tra[i];
			
				let val = traa.value;
			
				if(val.bp_eqn != undefined){
					if((traa.branch_select == true || traa.all_branches) && traa.branch == true){
						add_equation_to_list(eqn_list,val.bp_eqn,eqn_info);
					}
				}
				
				if(val.mean_eqn != undefined){
					if(traa.type == "exp(mean)" || traa.type == "gamma" || traa.type == "erlang" || traa.type == "log-normal" ||  traa.type == "period"){
						add_equation_to_list(eqn_list,val.mean_eqn,eqn_info);
					}
				}
				
				if(val.rate_eqn != undefined){
					if(traa.type == "exp(rate)"){
						add_equation_to_list(eqn_list,val.rate_eqn,eqn_info);
					}
				}
				
				if(val.scale_eqn != undefined){
					if(traa.type == "weibull"){
						add_equation_to_list(eqn_list,val.scale_eqn,eqn_info);
					}
				}
				
				if(val.shape_eqn != undefined){
					if(traa.type == "weibull"){
						add_equation_to_list(eqn_list,val.shape_eqn,eqn_info);
					}
				}
					
				if(val.cv_eqn != undefined){
					if(traa.type == "gamma" || traa.type == "log-normal"){
						add_equation_to_list(eqn_list,val.cv_eqn,eqn_info);
					}
				}
			}
		}
		
		// Adds parameters associated with simulation
		for(let i = 0; i < sp.sim_source.length; i++){
			let so = sp.sim_source[i];
	
			if(so.error != true){
				let eqn_info = {p:p, i:i, eso:"sim"};
			
				switch(so.type){
				case "Add Ind.":
					{
						let tab = so.table;
						for(let r = 0; r < tab.nrow; r++){
							for(let c = 2; c < tab.ncol; c++){
								add_eqn_filter(tab.ele[r][c],"sim_comp_prob",p,i,2,r,eqn_list);
							}
						}	
					}
					break;
					
				case "Test-and-cull":
					{
						let cb = so.spec.check_box.value; 
						for(let ci = 0; ci < cb.length; ci++){
							if(cb[ci].check){
								let Se = cb[ci].Se_eqn.te;
								add_equation_to_list(eqn_list,cb[ci].Se_eqn,eqn_info);
							}
						}
		
						add_equation_to_list(eqn_list,so.spec.Sp_eqn,eqn_info);
					}
					break;
				}
			}
		}		
		
		
		// Adds parameters associated with the observation process
		for(let i = 0; i < sp.inf_source.length; i++){
			let so = sp.inf_source[i];
	
			if(so.error != true){
				let eqn_info = {p:p, i:i, eso:"inf"};
				
				switch(so.type){
				case "Diag. Test": case "Test-and-cull":
					{
						let cb = so.spec.check_box.value; 
						for(let ci = 0; ci < cb.length; ci++){
							if(cb[ci].check){
								let Se = cb[ci].Se_eqn.te;
								add_equation_to_list(eqn_list,cb[ci].Se_eqn,eqn_info);
							}
						}
		
						add_equation_to_list(eqn_list,so.spec.Sp_eqn,eqn_info);
					}
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
		
		let eqn_info = {i:i, eso:"der"};
		
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
						let eqn_info2 = {p:p,i:i,cl:cl,c:c, eso:"inf"};
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
			let eqn_info2 = {p:p,i:i,cl:cl,tr:tr, eso:"inf"};
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
			//let flag = false;
			for(let j = 0; j < par.ndep_cont; j++){
				
				if(index == remove_prime(par.dep[j])){
					if(c < par.list[j].length){
						if(par.list[j][c] != old_name) error("Old name not correct");
						par.list[j][c] = new_name;
					}
				}
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
	if(par.time_dep) dep.pop();
	
	let name = "f~"+par.name;
	let full_name = "<e>"+name;
	if(dep.length > 0){
		full_name += "_";
		for(let j = 0; j < dep.length; j++){
			if(j > 0) full_name += ",";
			full_name += dep[j];
		}
	}
	
	dep.push("t");
	full_name += "(t)";
	
	full_name += "</e>";
	return {name:name, full_name:full_name, dep:dep, time_dep:true, param_fac:true};
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


/// Sets up a blank parameter specification
function param_blank(par)
{
	return par_find_template(par.list,par.ndep_cont);
}


/// Gets the number of elements in a parameter value
function param_num_element(par)
{
	let num = 1; 
	if(par.ndep_cont == 0) return num;
	for(let i = 0; i < par.list.length; i++){
		num *= par.list[i].length;
	}
	return num;
}


/// Determines if param is needed for simulation/inference
function param_needed(par,siminf)
{
	if(siminf == "all") return true;

	if(par.within == undefined) return false; // This can happen if defined quantities not used
	
	for(let i = 0; i < par.within.length; i++){
		let wi = par.within[i];
		
		let eso = wi.eqn_info.eso;
		switch(eso){
		case "model": case "dep": case "rep":
			return true;
		
		case "der":
			return false;
			
		case "sim":
			if(siminf == "sim") return true; 
			break;
			
		case "inf":
			if(siminf == "inf") return true; 
			break;
			
		default: 
			error("eso not found:"+eso);
			break;
		}
	}
	
	return false;
}

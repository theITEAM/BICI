"use strict";
// Functions which take an equation and extract information from the text

/// Extracts parameters from an equation
function extract_equation_properties(eqn)
{
	eqn.param = [];                                  // Parameters in the equation
	eqn.ind_eff = [];                                // Individual effect in the equation
	eqn.fix_eff = [];                                // Fixed effects in the equation
	eqn.pop = [];                                    // Populations in the equation
	eqn.sum = [];                                    // Sums in the equation
	eqn.warn = [];                                   // Warnings
	
	// These are used such that names can be changed
	eqn.sp_name_list = [];   
	eqn.index_name_list = [];
	eqn.comp_name_list = [];

	let icur = 0;

	let conv_res = detect_greek(eqn.te,0);
	eqn.te = conv_res.te;
	
	let lin = eqn.te;
	
	let res = check_brackets_match(lin);

	if(res.err == true){
		eqn.warn.push({te:res.msg, cur:res.op, len:1});
		return;
	}
		
	res = check_chnotallowed(lin);
	if(res.err == true){
		eqn.warn.push({te:res.msg, cur:res.op, len:1});
		return;
	}
	
	let lin2 = lin;
	let icur2 = icur;

	let two_variable_flag = false;

	let i = 0; 
	while(i < lin2.length){
		let type, right;
		while(i < lin2.length){
			if(i < lin2.length-4){
				if(lin2.substr(i,4) == "exp(" || lin2.substr(i,4) == "cos(" ||
					lin2.substr(i,4) == "sin(" || lin2.substr(i,4) == "log(" || lin2.substr(i,4) == "pow("){
					i += 4;
				}
			}
			
			if(i < lin2.length-5 && lin2.substr(i,5) == "step(") i += 5;			
			
			if(i < lin2.length-4 && lin2.substr(i,4) == "pow("){ two_variable_flag = true; i += 4;}			
			
			if(i < lin2.length-7 && lin2.substr(i,7) == "thresh("){ two_variable_flag = true; i += 7;}
			
			if(i < lin2.length-7 && lin2.substr(i,7) == "ubound("){ two_variable_flag = true; i += 7;}
			
			if(i < lin2.length-4 && lin2.substr(i,4) == "max("){ two_variable_flag = true; i += 4;}
	
			if(i < lin2.length-4 && lin2.substr(i,4) == "min("){ two_variable_flag = true; i += 4;}
			
			if(i < lin2.length-4 && lin2.substr(i,4) == "abs(") i += 4;			
			
			if(i < lin2.length-5 && lin2.substr(i,5) == "sqrt(") i += 5;			
			
			if(lin2.substr(i,1) == "|" && two_variable_flag == true) i++;
			
			if(is_sigma(lin2,i)){
				let ist = i;
				i = param_end(lin2,i+1,"("); 
				
				if(typeof i == 'string'){
					eqn.warn.push({te:i.replace("parameter","sum"), cur:icur2+ist, len:1});
					return;
				}
				
				let tex = lin2.substr(ist,i-ist);
				check_sum(tex,icur2+ist,eqn);
			}
			
			let ch = lin2.substr(i,1);
			
			if(ch == "{"){ type = "pop"; right = "}"; break;}
			if(ch == "["){ type = "ie"; right = "]"; break;}
			if(ch == "〈"){ type = "fe"; right = "〉"; break;}
			
			if(ch == "}"){
				eqn.warn.push({te:"There should not be a right bracket '}'", cur:icur2+i, len:1});
				return;
			}
			if(ch == "]"){
				eqn.warn.push({te:"There should not be a right bracket ']'", cur:icur2+i, len:1});
				return;
			}
		
			if(ch == "〉"){
				eqn.warn.push({te:"There should not be a right bracket '〉'", cur:icur2+i, len:1});
				return;
			}
			
			if(notparam_list.includes(ch) == false){
				type = "param"; break;
			}
			i++;
		}
			
		if(type != undefined){
			let ist = i;
			
			if(type == "param"){
				i = param_end(lin2,i); 
				
				if(typeof i == 'string'){
					eqn.warn.push({te:i, cur:icur2+ist, len:1});
					return;
				}

				let tex = lin2.substr(ist,i-ist);
			
				if(tex != "t"){
					let icur3 = icur2+ist;
				
					check_parameter(tex,icur3,eqn);
				}
			}
			else{
				while(i < lin2.length && lin2.substr(i,1) != right) i++;
				if(i == lin2.length){
					eqn.warn.push({te:"There should be a right bracket '"+right+"'", cur:icur2+ist, len:lin2.length-(ist)}); 
				}
				else{
					let tex = lin2.substr(ist+1,i-ist-1);
					let icur3 = icur2+ist+1;
					switch(type){
					case "pop": check_population(tex,icur3,eqn); break; 
					case "ie": check_ie(tex,"",icur3,eqn); break;  
					case "fe": check_fe(tex,"",icur3,eqn); break;  
					default: error("Option not recognised 40"); break;
					}
				}
				i++;
			}
		}
		else i++;

		icur = lin.length+1;
	}
	
	check_sum_bracket(eqn);
	
	check_indexes_match(eqn);
}


/// Checks that sums are encaptulated by brackets
function check_sum_bracket(eqn)
{
	let lin = eqn.te;
	
	let sigma = "Σ";
	
	let i = 0;
	while(i < lin.length-1){
		while(i < lin.length-1 && !is_sigma(lin,i)) i++;
		
		if(i < lin.length-1){
			let i_st = i;
			i += 1;
			if(lin.substr(i,1) == "^"){
				while(i < lin.length && lin.substr(i,1) != "_") i++;
				if(i == lin.length){
					eqn.warn.push({te:"The sum should contain a dependency '_'", cur:i_st, len:1}); 
					return;
				}
			}
			
			if(lin.substr(i,1) != "_"){
				eqn.warn.push({te:"The sum should be followed by a dependency '_'", cur:i, len:1}); 
			}
			else{
				i++;
				let ist = i;
				while(i < lin.length && lin.substr(i,1) != " " && lin.substr(i,1) != "(") i++;
				while(i < lin.length && lin.substr(i,1) == " ") i++;
				
				if(i == lin.length || lin.substr(i,1) != "("){
					eqn.warn.push({te:"Content in sums must be surrounded by round brackets (...)", cur:i, len:1}); 
				}
				
				// Checks to see if primed indexes appear in brackets
				let dep = lin.substr(ist,i-ist).split(",");
				for(let k = 0; k < dep.length; k++){
					dep[k] = dep[k].trim();
				}
				
				i++;
				let ibeg = i;
				let num = 1;
				
				while(i < lin.length){
					if(lin.substr(i,1) == "(") num++;
					if(lin.substr(i,1) == ")") num--;
					if(num == 0) break;
					i++;
				}					
				if(i == lin.length){ 
					eqn.warn.push({te:"The sum bracket does not match up", cur:ibeg, len:1}); 
				}
				
				let inside = lin.substr(ibeg,i-ibeg);
			
				for(let k = 0; k < dep.length; k++){
					if(!inside.includes(dep[k])){
						eqn.warn.push({te:"The sum does not include the index '"+dep[k]+"'", cur:ibeg-1, len:i-ibeg+2}); 
					}
				}
			}
		}
	}
}


/// check that primed indexes match with those on populations or parameter
function check_indexes_match(eqn)
{
	let list = [];
	let used = [];
	for(let i = 0; i < eqn.sum.length; i++){
		let su = eqn.sum[i];
		for(let j = 0; j < su.dep.length; j++){
			let de = su.dep[j];
			if(find_in(list,de) == undefined){ list.push(de); used.push(false);}
		}
	}
	
	for(let i = 0; i < eqn.param.length; i++){
		let par = eqn.param[i];
		for(let j = 0; j < par.dep.length; j++){
			let ind = par.dep[j];
			if(is_prime(ind)){
				let k = find_in(list,ind);
				if(k == undefined){
					if(eqn.type != "reparam_eqn"){
						eqn.warn.push({te:"Primed index '"+ind+"' for parameter '"+par.name+"' should be summed over"});
					}
				}					
				else used[k] = true;
			}
		}
	}

	for(let i = 0; i < eqn.pop.length; i++){
		let index = eqn.pop[i].index;
		for(let j = 0; j < index.length; j++){
			let ind = index[j];
			if(is_prime(ind)){
				let k = find_in(list,ind);
				if(k == undefined){
					eqn.warn.push({te:"Population primed index '"+ind+"' must be summed over", cur:0, len:0});
				}					
				else used[k] = true;
			}
		}
	}
	
	for(let i = 0; i < list.length; i++){
		if(used[i] == false){
			eqn.warn.push({te:"Index sum '"+list[i]+"' does not appear in either parameters or populations"});
		}
	}
}


/// Determines if sigma is at a particular point in a string
function is_sigma(st,i)
{
	let ch = "Σ";
	let len = ch.length;
  if(i <= st.length-len && st.substr(i,len) == ch) return true;
	return false;
}


/// Determines if a sum is defined correctly
function check_sum(tex,icur,eqn)
{
	let ch = "Σ";
	let i = ch.length;
	
	if(tex.substr(i,1) == "^"){
		i++;
		let ist = i;
		while(i < tex.length && tex.substr(i,1) != "_") i++;
		if(i == tex.length){
			eqn.warn.push({te:"The sum should contain a dependency2 '_'", cur:ist, len:1}); 
			return;
		}
		
		let sup = tex.substr(ist,i-ist);
		
		let fl = false;
		if(sup.length < 4) fl = true;
		else{
			if(sup.substr(0,4) != "max:") fl = true;
			let num = sup.substr(4);
			if(num == "" || isNaN(num)){
				eqn.warn.push({te:"Sum superscript number '"+num+"' not recognised.", cur:ist+6, len:sup.length-4}); 
				return;
			}
		}
		
		if(fl == true){
			eqn.warn.push({te:"Sum superscript '"+sup+"' not recognised.", cur:ist+2, len:sup.length}); 
			return;
		}
	}
			
	if(tex.substr(i,1) != "_"){
		eqn.warn.push({te:"The sum should be followed by '_'", cur:icur, len:i});
	}
	i++;
	
	let dep_str = tex.substr(i);
	if(dep_str == ""){
		eqn.warn.push({te:"The sum must have a dependency", cur:icur+i-1, len:1});
	}
	
	let dep = dep_str.split(",");
	
	eqn.sum.push({dep:dep});
	
	for(let j = 0; j < dep.length; j++){
		let de = dep[j];
		if(!is_prime(de)){
			eqn.warn.push({te:"The index '"+de+"' must have a prime if it is summed over", cur:icur+i, len:de.length});
		}
		
		let de_raw = remove_prime(de);
		let fl = false;
		for(let p = 0; p < model.species.length; p++){
			let sp = model.species[p];
			for(let cl = 0; cl < sp.cla.length; cl++){
				if(sp.cla[cl].index == de_raw) fl = true;
			}
		}
		if(fl == false){
			eqn.warn.push({te:"The index '"+de+"' is not in the model", cur:icur+i, len:de.length});
		}
		
		eqn.index_name_list.push({ index_name:remove_prime(de), icur:icur+i});
		i += de.length+1;
	}
}


/// Determines when a parameter definition ends
function param_end(st,i,op)
{
	let sub_brack = false;
	let sup_brack = false;
	let type = "normal";
	
	do{
		let ch = st.substr(i,1);
		if(ch == op) break;
		
		if(ch == "^"){
			if(type != "normal") return "Problem understanding parameter definition.";
			type = "sup";
		}
			
		if(ch == "_"){
			if(type != "normal" && type != "sup"){
				return "Problem understanding parameter definition.";
			}
			type = "sub";
		}
		
		if(ch == "("){
			if(i > 0 && st.substr(i-1,1) == "^") sup_brack = true;
			else{
				if(i > 0 && st.substr(i-1,1) == "_") sub_brack = true;
				else break;
			}
		}
		
		if(ch == ")"){
			if(!(type == "sup" && sup_brack == true) && !(type == "sub" && sub_brack == true)) break;
			sub_brack = false; sup_brack = false;
		}
		
		if(paramend_list.includes(ch)) break;
		i++;
	}while(i < st.length);
	
	if(i <= st.length - 3){
		if(st.substr(i,3) == "(t)") i += 3;
	}

	return i;
}			

		
/// Checks that parameter is correctly defined			
function check_parameter(te,icur,eqn) 
{
	let p_name = eqn.p_name;
	let p; if(p_name != undefined) p = find(model.species,"name",p_name);
	let sp; if(p != undefined) sp = model.species[p];
	let cl; if(eqn.cl_name != undefined) cl = find(sp.cla,"name",eqn.cl_name);

	let time_dep = false;
	
	if(te.length > 3){
		switch(te.substr(te.length-3,3)){
			case "(t)": time_dep = true; te = te.substr(0,te.length-3); break;
		}
	}
	
	let dep=[]; 
	let spl = te.split("_");
	if(spl.length > 2){
		eqn.warn.push({te:"Did not expect more than one underscore", cur:icur, len:te.length});
	}
	
	let name = spl[0];
	
	check_valid_name(name,icur,eqn);
	
	// Checks if compartment names in parameter name 
	if(eqn.type == "trans"){
		let spl = name.split("^");
		if(spl.length == 2){
			let icur2 = icur+spl[0].length+1;
			let tex = spl[1];
			let tex2 = remove_bracket(tex);
			if(tex2.length < tex.length) icur2++; 
			let spl2 = tex2.split("→");
			
			let p = find(model.species,"name",eqn.p_name);
			let sp = model.species[p];
			let cl = find(sp.cla,"name",eqn.cl_name);
			let claa = sp.cla[cl];
	
			for(let i = 0; i < spl2.length; i++){
				let c = hash_find(claa.hash_comp,spl2[i]);
				if(c != undefined){
					eqn.comp_name_list.push({ p_name:eqn.p_name, cl_name:eqn.cl_name, comp_name:spl2[i], icur:icur2});
				}
				icur2 += spl2[i].length+1;
			}
		}
	}
	
	if(spl.length == 2){
		let icur2 = icur+name.length+1;
		let tee = spl[1];
		
		if(eqn.mode == "param only"){
			eqn.warn.push({te:"This parameter cannot have any dependency", cur:icur2, len:tee.length});
		}
		else{
			let tee2 = remove_bracket(tee);
			if(tee2.length < tee.length) icur2++; 

			dep = tee2.split(",");
		
			for(let i = 0; i < dep.length; i++){
				let na = dep[i];
				let index = remove_prime(na);
			
				if(eqn.type != "derive_eqn" && eqn.type != "derive_param" && eqn.type != "reparam_eqn" && eqn.type != "test"){
					if(sp == undefined) error("SHould be defined");
					let cl2; cl2 = find(sp.cla,"index",index);	
					if(cl2 == undefined){
						/// Looks for potential indices in other species
						let cl3;
						for(let p2 = 0; p2 < model.species.length; p2++){
							if(p2 != p){
								cl3 = find(model.species[p2].cla,"index",index);
								if(cl3 != undefined) break;
							}
						}
						if(cl3 == undefined){
							eqn.warn.push({te:"'"+index+"' is not a classification index name", cur:icur2, len:index.length});
						}
						else{
							if(index == na){
								eqn.warn.push({te:"The index '"+index+"' must be primed because it is in another species", cur:icur2, len:index.length});
							}
						}
					}
					else{
						if(cl2 == cl && index == na){
							eqn.warn.push({te:"The index '"+na+"' must be primed", cur:icur2, len:na.length});
						}
					}
				}
				
				eqn.index_name_list.push({ p_name:p_name, index_name:index, icur:icur2});

				icur2 += na.length+1;
			}
		}
	}
	
	if(time_dep == true) dep.push("t");
	
	let par = { name:name, p_name:eqn.p_name, cl_name:eqn.cl_name, dep:dep, type:eqn.type, time_dep:time_dep};

	par.full_name = param_name(par);
	
	eqn.param.push(par);
}


/// Determines if two vectors are the same
function equal_vec(vec1,vec2)
{
	if(vec1.length != vec2.length) return false;
	for(let i = 0; i < vec1.length; i++){
		if(vec1[i] != vec2[i]) return false;
	}
	
	return true;
}


/// Determines if a parameter is primed
function is_prime(te)
{
	if(te.length > 0 && te.substr(te.length-1,1) == "'") return true;
	return false;
}


/// Removes any primes
function remove_prime(te)
{
	let i = te.length-1;
	while(i >= 0 && te.substr(i,1) == "'") i--;
	return te.substr(0,i+1);
}


/// Removes round brakets from a string
function remove_bracket(te)
{
	if(te.length > 1){
		if(te.substr(0,1) == "(" && te.substr(te.length-1,1) == ")"){
			return te.substr(1,te.length-2);
		}
	}
	return te;	
}
				

/// Checks that parameter is correctly defined			
function check_population(te,icur,eqn) 
{
	let p_name = eqn.p_name;
	
	let p; if(p_name != undefined) p = find(model.species,"name",p_name);
	
	// Removes species filter
	let i = 0; while(i < te.length && te.substr(i,1) != ":") i++;
	if(i < te.length){
		p_name = te.substr(0,i);
		eqn.sp_name_list.push({ p_name:p_name, icur:icur});
		
		p = find(model.species,"name",p_name);
		if(p == undefined){
			eqn.warn.push({te:"The species '"+p_name+"' is not recognised", cur:icur, len:i});
			return;
		}
		icur += i+1;
		te = te.substr(i+1);
	}
	
	if(p == undefined){
		if(model.species.length == 1) p = 0;
		else{
			eqn.warn.push({te:"Species must be defined.", cur:icur, len:te.length});
			return;
		}
	}
	
	// Removes any individual effects
	do{
		let j = 0; while(j < te.length && te.substr(j,1) != "[") j++;
		if(j < te.length){
			let jst = j;
			while(j < te.length && te.substr(j,1) != "]") j++;
			if(j == te.length){
				eqn.warn.push({te:"Individual effect should have a right bracket ']'", cur:icur+jst, len:j-jst});
			}
			else{
				let ie = te.substr(jst+1,j-(jst+1));
				check_ie(ie,p_name,icur+jst+1,eqn);
				te = te.substr(0,jst)+te.substr(j+1);
			}
		}
		else break;
	}while(true);
	
	// Removes any fixed effects
	{
		let j = 0; 
		while(j < te.length){
			while(j < te.length && te.substr(j,1) != "〈") j++;
			if(j < te.length){
				let jst = j;
				while(j < te.length && te.substr(j,1) != "〉") j++;
				if(j == te.length){
					eqn.warn.push({te:"Fixed effect should have a right bracket '〉'", cur:icur+jst, len:j-jst});
					return;
				}
				else{
					let fe = te.substr(jst+1,j-(jst+1));
					check_fe(fe,p_name,icur+jst+1,eqn);
					te = te.substr(0,jst)+te.substr(j+1);
					j = jst;
				}
			}
		}
	}
	
	// Checks that the population being set is consistent 
	let sp = model.species[p];
	
	let te_comp = te;
	
	let cl_ref=[];

	let pop = {index:[]};
	
	// Goes through compartment names
	let icur3 = icur;
	if(te_comp != "" && te_comp != "All"){
		let spl = te_comp.split(",");
		
		for(let i = 0; i < spl.length; i++){
			let tex = spl[i];
			
			let index = remove_prime(tex);
		
			let cl_sp = find(sp.cla,"index",index);      // Detects if an index
			
			if(cl_sp != undefined){
				pop.index.push(tex);
				
				eqn.index_name_list.push({ index_name:index, icur:icur3});
			
				if(cl_ref[cl_sp] != undefined){
					eqn.warn.push({te:"Cannot have multiple compartments in the same classification", cur:cl_ref[cl_sp],len:icur3+te.length-cl_ref[cl_sp]});
				}
				cl_ref[cl_sp] = icur3;
			}
			else{
				let warn;
				let cl_sp;
				
				let spl2 = tex.split("|");
				
				let icur4 = icur3;
				for(let k = 0; k < spl2.length; k++){
					let te2 = spl2[k];
					let co = find_comp_from_name(te2,p_name);
					if(co == undefined){
						warn = "Compartment '"+te2+"' not found";
						if(index != tex) warn = "Index '"+te2+"' not found";
					}
					else{
						if(co.warn != undefined) warn = co.warn;
						else{
							let cl_name = sp.cla[co.cl].name;
								
							let comp_name = sp.cla[co.cl].comp[co.c].name;
							eqn.comp_name_list.push({ p_name:p_name, cl_name:cl_name, comp_name:comp_name, icur:icur4});
							
							if(cl_sp == undefined) cl_sp = co.cl;
							else{
								if(cl_sp != co.cl){
									warn = "The compartments '"+spl2[k-1].trim()+"' and '"+spl2[k].trim()+"' are on different classifications";
								}
							}
						}
					}
					icur4 += te2.length+1;
				}
				if(warn != undefined){
					eqn.warn.push({te:warn, cur:icur3, len:tex.length});
				}
				
				if(cl_ref[cl_sp] != undefined){
					eqn.warn.push({te:"Cannot have multiple compartments in the same classification", cur:cl_ref[cl_sp],len:icur3+te.length-cl_ref[cl_sp]});
				}
				cl_ref[cl_sp] = icur3;
			}
			
			icur3 += te_comp.length+1;
		}
	}
	
	eqn.pop.push(pop);
}


/// Checks that individual effect is correctly defined			
function check_ie(te,p_name_filter,icur,eqn) 
{
	check_valid_name(te,icur,eqn);
	
	let i = find(eqn.ind_eff,"name",te);
	if(i == undefined){
		eqn.ind_eff.push({name:te, p_name_filter:p_name_filter});
	}
	else{
		if(eqn.ind_eff[i].p_name_filter != p_name_filter){
			eqn.warn.push({te:"Species filters for individual effect "+te+" do not match up"});
		}
	}
}


/// Checks that individual effect is correctly defined			
function check_fe(te,p_name_filter,icur,eqn) 
{
	check_valid_name(te,icur,eqn);
	
	let i = find(eqn.fix_eff,"name",te);
	if(i == undefined){
		eqn.fix_eff.push({name:te, p_name_filter:p_name_filter});
	}
	else{
		if(eqn.ind_eff[i].p_name_filter != p_name_filter){
			eqn.warn.push({te:"Species filters for fixed effect "+te+" do not match up"});
		}
	}
}


/// Finds a compartment based on a name
function find_comp_from_name(te,p_name)
{
	te = te.trim();
	if(te == "") return {warn:"Compartment not specfied"};
	
	if(p_name == undefined){
		for(let p = 0; p < model.species.length; p++){
			let sp = model.species[p];
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				for(let c = 0; c < claa.ncomp; c++){
					if(claa.comp[c].name == te){
						return {p:p, cl:cl, c:c};
					}
				}
			}
		}
	}
	else{
		let p = find(model.species,"name",p_name);
		if(p == undefined){
			return {warn:"Problem with expression"};
		}
		
		for(let cl = 0; cl < model.species[p].ncla; cl++){
			let claa = model.species[p].cla[cl];
			for(let c = 0; c < claa.ncomp; c++){
				if(claa.comp[c].name == te){
					return {p:p, cl:cl, c:c};
				}
			}
		}
	}
	return;
	
	//return { warn:"Could not find compartment '"+te+"'"};
}


/// Checks that a string is a valid name for a parameter
function check_valid_name(name,icur,eqn)
{
	for(let j = 0; j < name.length; j++){
		let ch = name.substr(j,1);
		let k = 0; while(k < char_not_allowed.length && ch != char_not_allowed[k]) k++;
		if(k < char_not_allowed.length){
			if(ch == "_") ch = "underscore";
			eqn.warn.push({te:"In '"+name+"' the character '"+ch+"' is not expected", cur:icur+j, len:1});
		}
	}
}


/// Gets the properties of a parameter from the string
function get_param_prop(st)
{
	st = st.trim();
	let time_dep = false;
	
	if(st.length > 3){
		let end = st.substr(st.length-3,3);
		if(end == "(t)"){ time_dep = true; st = st.substr(0,st.length-3);}
	}
	
	let j = 0;
	while(j < st.length && st.substr(j,1) != "_") j++;

	let name = st.substr(0,j);
	let name_raw = name;
	let sup = "";
	
	let jj = 0;
	while(jj < name.length && name.substr(jj,1) != "^") jj++;
	if(jj < name.length){
		sup = name.substr(jj+1);
		sup = remove_bracket(sup);
		name_raw = name.substr(0,jj);
	}

	let dep = [];
	let sub = "";
	
	if(st.substr(j,1) == "_"){
		j++;
		let jst = j; 
		while(j < st.length && st.substr(j,1) != "[") j++;
		sub = st.substr(jst,j-jst);
		sub = remove_bracket(sub);
		dep = sub.split(",");
	}
					
	if(time_dep == true) dep.push("t");

	return {name:name, name_raw:name_raw, dep:dep, sub:sub, sup:sup, time_dep:time_dep, type:unset_type};
}


/// Turns a dependency array into a string
function string_dep(dep)
{
	if(dep.length == 0) return "none";
	
	let te = "";
	for(let i = 0; i < dep.length; i++){
		if(i != 0) te += ",";
		te += dep[i];
	}
	
	return te;
}


/// Returns the overall equation dependency
function equation_dep(eqn)
{
	extract_equation_properties(eqn);
	
	let time_fl = false;
	
	let dep=[];
	for(let i = 0; i < eqn.param.length; i++){
		let par = eqn.param[i];
		for(let j = 0; j < par.dep.length; j++){
			let de = par.dep[j];
			if(de == remove_prime(de)){
				if(de == "t") time_fl = true;
				else{
					if(find_in(dep,de) == undefined){
						dep.push(de);
					}
				}
			}
		}
	}
	
	for(let i = 0; i < eqn.pop.length; i++){
		let po = eqn.pop[i];
		for(let j = 0; j < po.index.length; j++){
			let de = po.index[j];
			if(find_in(dep,de) == undefined){
				dep.push(de);
			}
		}
	}
	
	if(eqn.te.includes("{") == true) time_fl = true;
	
	if(time_fl == true) dep.push("t");

	return dep;
}


/// Checks that all brackets match in the equation
function check_brackets_match(te)
{
	let brack_list = [];
	for(let i = 0; i < te.length; i++){
		let ch = te.substr(i,1);
		if(ch == "(" || ch == "[" || ch == "{") brack_list.push({i:i, ch:ch});
		if(ch == ")" || ch == "]" || ch == "}"){
			if(brack_list.length == 0){
				return err("The bracket '"+ch+"' does not match up",i);
			}
			
			let chlast = brack_list[brack_list.length-1].ch;
			if((chlast == "(" && ch == ")") || (chlast == "[" && ch == "]") ||(chlast == "{" && ch == "}")){
				brack_list.pop();
			}
			else{
				return err("The bracket '"+chlast+"' does not match up",brack_list[brack_list.length-1].i);
			}
		}		
	}
	
	if(brack_list.length != 0){
		let last = brack_list[brack_list.length-1];
		return err("The bracket '"+last.ch+"' does not match up",last.i);
	}
	
	return success();
}


/// Checks that the expression does not include ceratin characters which are not allowed
function check_chnotallowed(te)
{
	for(let i = 0; i < te.length; i++){
		let ch = te.substr(i,1);
		if(chnotallowed.includes(ch)){
			return err("The character '"+ch+"' is not allowed",i);
		}
	}
	
	return success();
}

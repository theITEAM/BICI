"use strict";
// General purpose functions

/// Prints to the console
function pr(te)
{ 
	if(te == undefined){ prr("UNDEFINED"); return;}
	console.log(copy(te));
}


/// Prints to the console
function prr(te)
{ 
	console.log(te);
}


/// Prints to the console
function comment(te)
{ 
	console.log(te);
}


/// Prints a text error to the console
function error(te)
{
	prr("ERROR:");
	console.log(te);
}


/// Copies an object
function copy(inp){
	return JSON.parse(JSON.stringify(inp));
}


/// Finds a value in list of objects
function find(arr,prop,val)
{
	let i = 0; while(i < arr.length && arr[i][prop] != val) i++;
	if(i == arr.length) return;
	return i;
}



/// Finds a value in list of objects (irrespective of case)
function find_nocase(arr,prop,val)
{
	let i = 0; while(i < arr.length && arr[i][prop].toLowerCase() != val.toLowerCase()) i++;
	if(i == arr.length){
		//error("Cannot find array value. prop:"+prop+" val:"+val+"  arr:"+arr);
		return;
	}
	return i;
}


/// Finds a value in list of objects
function find_in(arr,val)
{
	let i = 0; while(i < arr.length && arr[i] != val) i++;
	if(i == arr.length)	return;
	return i;
}


/// Finds a value in list of string case-independent
function find_string_in(arr,val)
{
	let i = 0; while(i < arr.length && arr[i].toLowerCase() != val.toLowerCase()) i++;
	if(i == arr.length)	return;
	return i;
}


/// When typing this convert \alpha to greek alpha
function detect_greek(te,curi)
{
	let i = 0;
	while(i < te.length){
		while(i < te.length && te.substr(i,1) != "\\") i++;
	
		if(i < te.length){
			i++;
			for(let j = 0; j < greek_latex.length; j++){
				let wo = greek_latex[j][0];
				let len = wo.length;
				if(i <= te.length-len){
					if(te.substr(i,len) == wo){
						te = te.substr(0,i-1)+greek_latex[j][1]+te.substr(i+len);
						if(curi > i) curi -= len;
					}
				}
			}
		}
	}
	
	return {te:te,curi:curi};
}


/// Returns an error from a function
function err(msg,op)
{
	return { err:true, msg:msg, op:op};
}


/// Returns success from a function
function success()
{
	return { err:false};
}


/// Gets the font based on size (and potenitally if bold ot italic)
function get_font(fac,te,type)
{
	if(type == undefined) type = "arial";
	let fo = Math.round(fac*inter.sca)+"px "+type;
	if(te) return te+" "+fo;
	return fo;
}


/// Rounds to nearest number, e.g. converts 0.199900009999 to 0.2
function round_small(val)
{
	//return Math.round(val*10000000000)/10000000000;
	return Number(val.toPrecision(6));
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


/// Write an array as a string
function stringify(arr)
{
	let te="";
	for(let i = 0; i < arr.length; i++){
		if(i != 0) te += ","; 
		te += arr[i];
	}
	return te;
}


/// Removes <e> </e> equation quotes
function remove_eq_quote(te)
{
	if(typeof te == "string"){
		if(te.length >= 3 && te.substr(0,3) == "<e>") te = te.substr(3);
		if(te.length >= 4 && te.substr(te.length-4,4) == "</e>") te = te.substr(0,te.length-4);
	}
	return te;
}


/// Removes quotation marks from a string
function remove_quote(te)
{
	if(te.substr(0,1) == "\"" && te.substr(te.length-1,1) == "\""){
		te = te.substr(1,te.length-2);
	}
	return te;
}


/// Returns a darker version of a colour
function dark_colour(col)
{ 
	if(col == BLACK) return DDGREY;

	const frac = 0.7;
	let bigint = parseInt(col.substring(1), 16);	
	let r = (bigint >> 16) & 255;
	let g = (bigint >> 8) & 255;
	let b = bigint & 255;
	
	return "#" + hex(frac*r) + hex(frac*g) + hex(frac*b);
}


/// Returns a lighter version of a colour
function light_colour(col,frac) 
{ 
	if(frac == undefined) frac = 0.13;

	let bigint = parseInt(col.substring(1), 16);
	let r = 255-(255-((bigint >> 16) & 255))*frac;
	let g = 255-(255-((bigint >> 8) & 255))*frac;
	let b = 255-(255-(bigint & 255))*frac;
	return "#" + hex(r) + hex(g) + hex(b);
}


/// If dark, makes light, and vice-versa
function shift_colour(col,fr)
{ 
	let bigint = parseInt(col.substring(1), 16);	
	let r = (bigint >> 16) & 255;
	let g = (bigint >> 8) & 255;
	let b = bigint & 255;
	
	if(r+g+b < 150){
		let frac = 0.8; if(fr != undefined) frac = fr;
		r = 255-(255-r)*frac;
		g = 255-(255-g)*frac;
		b = 255-(255-b)*frac;
	}
	else{
		let frac = 0.8; if(fr != undefined) frac = fr;
		r *= frac;
		g *= frac;
		b *= frac;
	}
	
	return "rgb(" + Math.floor(r) +","+Math.floor(g)+","+Math.floor(b)+")";
}


/// Converts decimal to hexidecimal
function hex(c) 
{
	let hex = (Math.floor(c)).toString(16);
	return hex.length == 1 ? "0" + hex : hex;
}


/// Based on the individual classifications, this generates all the compartments in the system
function get_glob_comp(p)
{
	let sp = model.species[p];
	
	let glob_comp=[];
	
	if(sp.ncla == 0) return glob_comp;
	
	let index=[];	
	for(let cl = 0; cl < sp.ncla; cl++){
		if(sp.cla[cl].ncomp == 0) return glob_comp;
		index[cl] = 0;
	}
	
	do{
		let co = { cla:[]};
		
		for(let cl = 0; cl < sp.ncla; cl++){
			co.cla[cl] = sp.cla[cl].comp[index[cl]].name;
		}
		glob_comp.push(co);
			
		let cl = 0; 
		let flag;
		do{
			flag = false;
			index[cl]++; if(index[cl] >= sp.cla[cl].ncomp){ index[cl] = 0; cl++; flag = true;}
		}while(flag == true && cl < sp.ncla);
		if(cl == sp.ncla) break;
	}while(true);
		
	return glob_comp;
}


/// Given a classification name return species
function get_p_cl_from_claname(name)
{
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let cl = 0; cl < sp.ncla; cl++){
			if(sp.cla[cl].name == name) return {p:p, cl:cl};
		}
	}
}


/// Sets the precision of the output	
function precision(x,digit)
{
	if(isNaN(x)) error(x+" is not a number");
	x = Number(x);

	if(digit == undefined) digit = 5;
	
	let mod = x; if(mod < 0) mod = -mod; 
	
	let digit_min = Math.floor(1+Math.log10(mod));
	if(digit_min < 1) digit_min = 1;

	if(digit < digit_min) digit = digit_min;
	
	let val = Number(x).toPrecision(digit);
	for(let j = digit-1; j >= digit_min; j--){
		let val2 = x.toPrecision(j);
		if(Number(val) == Number(val2)) val = val2; else break;
	}
 
	return val;
}


/// Updates a model property and sets update_model_needed if changed
function update(ob,te)
{
	if(ob.te == undefined){ error("Should not be undefined"); return;}
	if(ob.te != te){
		ob.te = te;
		update_param();
	}
}


/// Converts text to italic
function italic(te)
{
	return "<e>"+te+"</e>";
}


/// Sets for the parameters to be updated
function update_param()
{
	model.update_model_needed = true;
}


/// Samples Poisson. This version should not fail when mean becomes large
function poisson_sample(mean)
{
	let left = mean, k = 0, p = 1, STEP=500;
 
	do{
		k++;
		p *= Math.random();
		while(p < 1 && left > 0){
			if(left > STEP){
				p *= Math.exp(STEP);
				left -= STEP;
      }
			else{
				p *= Math.exp(left);
        left = 0;
			}
		}
	}while(p > 1);
	return k-1;
}


/// Draws a sample from a Poisson distribution
function poisson_sample_old(mean)
{
	var L = Math.exp(-mean);
	var p = 1.0;
	var k = 0;

	do{
    k++;
    p *= Math.random();
	}while (p > L);
	return k-1;
}


/// Checks the Poisson sampling is being done correctly
function check_poisson()
{
	let m = 10000;
	let max = 100000;
	let av = 0, av2 = 0;
	for(let l = 0; l < max; l++){
		let sa = poisson_sample(m);
		av += sa; av2 += sa*sa;
	}	
	av /= max; av2 /= max;
	
	prr(m+" "+av+" "+(av2-av*av)+" Check Poisson");
}

			
/// Draws a sample from a normal distribution
function normal_sample(mean,sd)
{
	let u = 1 - Math.random(); // Converting [0,1) to (0,1]
  let v = Math.random();
  let z = Math.sqrt(-2.0*Math.log(u))*Math.cos(2.0*Math.PI*v);
  return mean + z*sd;
}


/// Draws a sample from a negative binomial distribution
function neg_binommial_sample(mean,p)
{
	let r = Math.round(mean*p/(1-p));
	if(r < 1) r = 1;
	
	let nfail = 0, nsuc = 0; 
	do{
		if(Math.random() < p) nsuc++;
		else nfail++;
	}while(nsuc < r);
	return nfail;
}


/// Copies a source to a destination but omits certain data 
function copy_strip(source,dest)
{
	if(typeof source != 'object'){
		error("SHOULD BE OBJECT");
		prr(source);
	}
	
	let is_arr = Array.isArray(source);
	
	if(is_arr){
		let num = dest.length - source.length;	
		for(let i = 0; i < num; i++) dest.pop(); 
	}
	for(let ele in source){
		let cont = false;
		
		switch(ele){
		case "factor_weight":
			dest.weight_desc = get_weight_desc(source);
			//dest.dim = get_value_dim(source);
			break;
			
		case "value": 
			if(source.ind_list == undefined){
				dest.value_desc = get_value_desc(source);
				dest.dim = get_value_dim(source);
			}
			
			if(!Array.isArray(source.value)) dest.value = source.value;
			break;
			
		case "prior_split": 
			dest.prior_split_desc = get_prior_split_desc(source);
			break;
			
		case "table": case "A_value": case "X_value": case "ind_list":
		case "list": case "comb_list": case "reparam_param_list": case "prior_param_list":
			break;
		
		case "param": case "species": 
		case "inf_source": case "sim_source":  case "ppc_source": case "gen_source":
		case "ind_eff_group": case "A_matrix": case "fix_eff": case "X_vector":
		case "sim_res": case "inf_res": case "ppc_res": case "plot_filter":
			cont = true;
			break;
			
		default:		
			if(is_arr) cont = true;
			else dest[ele] = source[ele];
			break;
		}
		
		if(cont == true){
			if(dest[ele] == undefined){
				if(Array.isArray(source[ele])) dest[ele]=[];
				else dest[ele]={};
			}
			
			copy_strip(source[ele],dest[ele]);
		}
	}
}


/// Strips out parts with large memory footprint
function strip_heavy(source)
{
	let source_strip;
	if(Array.isArray(source)) source_strip=[];
	else source_strip={};
	copy_strip(source,source_strip);
	return source_strip;	
}


/// Determines if parameter is in the observation model
function is_in_obsmodel(par)
{
	let i = find(eqn_types,"name",par.type);
	if(i != undefined && eqn_types[i].obs_model == true) return true;
	return false;
}


/// Determines if a string is a valid colour
function isColor(strColor)
{
	return true;
	var s = new Option().style;
	s.color = strColor;
	return s.color.length > 0;
}


/// Takes a transition definition and extracts information from it
// exp(rate:r)
function extract_trans_def(value)
{
	value = value.trim();
	if(value == "") return;
	let spl = value.split("(");
	if(spl.length == 0) return; 
	
	if(value.substr(value.length-1,1) != ")") return;

	let type, rate, mean, cv, shape, scale, time;

	switch(spl[0].toLowerCase()){
	case "exp":
		{
			rate = get_prop(value,"rate:","end");
			if(rate != undefined) type = "exponential";
		}
		break;
		
	case "gamma":
		mean = get_prop(value,"mean:","cv:");
		cv = get_prop(value,"cv:","end");
		if(mean != undefined && cv != undefined) type = "gamma";
		break;
	
	case "erlang":
		mean = get_prop(value,"mean:","shape:");
		shape = get_prop(value,"shape:","end");
		if(mean != undefined && shape != undefined) type = "erlang";
		break;
		
	case "log-normal":
		mean = get_prop(value,"mean:","cv:");
		cv = get_prop(value,"cv:","end");
		if(mean != undefined && cv != undefined) type = "log-normal";
		break;
		
	case "weibull":
		scale = get_prop(value,"scale:","shape:");
		shape = get_prop(value,"shape:","end");
		if(scale != undefined && shape != undefined) type = "weibull";
		break;
		
	case "period":
		time = get_prop(value,"time:","end");
		if(time != undefined) type = "period";
		break;
	}
	
	if(type == undefined) return;
	
	return { type:type, rate:rate, mean:mean, cv:cv, shape:shape, scale:scale, time:time};
}


/// Tries to get a property from a string
function get_prop(value,prop,end)
{
	let i = value.indexOf(prop);
	if(i != -1){
		i += prop.length;
		if(end == "end"){
			let val = value.substr(i,value.length-1-i).trim();
			return val;
		}
		else{
			let i2 = value.indexOf(end);
			if(i2 != -1){
				while(i2 > 0 && value.substr(i2,1) != ",") i2--;
				if(i2 > 0){
					let val = value.substr(i,i2-i).trim();
					return val;
				}
			}
		}
	}
}


/// Determines the total number of elements for parameter
function num_element(par)
{
	let num = 1;
	for(let k = 0; k < par.dep.length; k++) num *= par.list[k].length;
	
	return num;
}


/// Gets the clock time (in milliseconds)
function clock()
{
	return (new Date()).getTime();
}


/// Makes sure the camera is up-to-date for any clones
function clone_camera(p,cl)
{
	let list = model.find_clones(p,cl);	
	if(list.length > 0){
		let cam = model.species[p].cla[cl].camera;
	
		for(let li of list){
			let cam2 = model.species[li.p].cla[li.cl].camera;
			cam2.x = cam.x; cam2.y = cam.y; cam2.scale = cam.scale; cam2.ruler = cam.ruler;
		}
		model.check_clones();
	}
}


/// Calculates if two numbers is different (subject to numerical noise) 
function dif(a,b)
{
	let dif = a-b; if(dif < 0) dif = -dif;
	if(dif > TINY) return true;
	return false;
}


/// Determines if string beg is at the begining of te
function begin(te,beg)
{
	if(te.length < beg.length) return false;
	if(te.substr(0,beg.length) == beg) return true;
	return false;
}


/// Gets inf_res
function get_inf_res()
{
	if(inter.page.length > 0){
		let pag = inter.page[inter.pa];
		switch(pag.name){
		case "Simulation": return model.sim_res; 
		case "Inference": return model.inf_res; 
		case "Post. Simulation": return model.ppc_res;
		}
	}
}


/// Split up a string accounting for brackets and quotation marks
function split_with_bracket(s,delimiter)
{                              
 let splits=[];                       
 
  let quoteon = false;
	let num_brac = 0;
	
	let j = 0;
	for(let i = 0; i < s.length; i++){
		let letter = s.substr(i,1);
		
		if(letter == "\""){
			if(quoteon == false) quoteon = true; else quoteon = false;
		}
	
		if(letter == "(") num_brac++;
		if(letter == ")") num_brac--;
		
		if(letter == delimiter && num_brac == 0 && quoteon == false){
			splits.push(s.substr(j,i-j).trim()); j = i+1; 
		}
	}
	splits.push(s.substr(j,s.length-j).trim());

	return splits;                                           
}


/// CHecks if two hash tables are the same
function check_hash_same(store1,store2)
{
	if(store1.length != store2.length) error("Hash length wrong");
	else{
		for(let i = 0; i < store1.length; i++){
			let vec1 = store1[i];
			let vec2 = store2[i];
			if(vec1 == undefined || vec2 == undefined){
				if(vec1 != undefined){
					if(vec1.length > 0) error("Hash defined wrong1");
				}
				
				if(vec2 != undefined){
					if(vec2.length > 0){
						prr("vecs");
						prr(vec1);
						prr(vec2);
						error("Hash defined wrong2");
					}
				}
			}
			else{
				if(vec1.length != vec2.length) error("Hash length wrong2");
				else{
					vec1.sort( function(a, b){ return a.num-b.num});
					vec2.sort( function(a, b){ return a.num-b.num});
				
					for(let j = 0; j < vec1.length; j++){
						if(vec1[j].num != vec2[j].num){
							error("Hash num wrong");
						}
						if(vec1[j].ref != vec2[j].ref) error("Hash ref wrong");
					}
				}
			}
		}
	}
}


/// Checks that the hash tables for comp and tra are correctly specified
function check_hash_comp_tra()
{
	prr("CHECK HASH");
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let cl = 0; cl < sp.cla.length; cl++){
			let claa = sp.cla[cl];
			
			{
				let hash = [];
				for(let c = 0; c < claa.comp.length; c++){
					hash_add(hash,claa.comp[c].name,c);
				}
				check_hash_same(hash,claa.hash_comp);
			}
		
			{
				let hash = [];
				for(let c = 0; c < claa.tra.length; c++){
					hash_add(hash,claa.tra[c].name,c);
				}
				check_hash_same(hash,claa.hash_tra);
			}
		}
	}
}


/// Checks that nsource and nsink are correctly specified
function check_nsource_nsink()
{
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let cl = 0; cl < sp.cla.length; cl++){
			let claa = sp.cla[cl];
			
			let nsource = 0, nsink = 0;
			for(let j = 0; j < claa.tra.length; j++){
				let tr = claa.tra[j];
				if(tr.i == SOURCE) nsource++;
				if(tr.f == SINK) nsink++;
			}
			
			if(nsource != claa.nsource) error("nsource problem");
			if(nsink != claa.nsink) error("nsink problem");
		}
	}
}


/// Returns the size of an object
function obj_memory(obj)
{
	var bytes = 0;
	if(obj !== null && obj !== undefined) {
		switch (typeof obj) {
			case "number":
				bytes += 8;
				break;
			case "string":
				bytes += obj.length * 2;
				break;
			case "boolean":
				bytes += 4;
				break;
			case "object":
				var objClass = Object.prototype.toString.call(obj).slice(8, -1);
				if (objClass === "Object" || objClass === "Array") {
					for (var key in obj) {
						if (!obj.hasOwnProperty(key)) continue;
						//prr("insde");
						//prr(obj[key]);
						bytes += obj_memory(obj[key]);
					}
				} else bytes += obj.toString().length * 2;
				break;
		}
	}
  return bytes;
}


/// Truncates a text string to a certain length
function trunc(te,len)
{
	if(te.length < len) return te;
	return te.substr(0,len)+"...";
}


/// Determines if a string begins with another string
function begin_str(st,st2)
{
	let tr = st.trim();
	if(st.length < st2.length) return false;
	if(st.substr(0,st2.length) == st2) return true;
	return false;
}


/// Checks if a name is invalid
function check_invalid_name(name)
{
	for(let i = 0; i < invalid_name.length; i++){
		if(name.toLowerCase() == invalid_name[i].toLowerCase()) return true;
	}
	return false;
}

	
/// Checks the name is valid
function check_name_warn(name,te) 
{
	if(check_invalid_name(name)){
		return te+" cannot use reserved word '"+name+"'";
	}
	
	for(let i = 0; i < name.length; i++){
		let ch = name.substr(i,1);
		if(name_notallow.includes(ch)){
			return te+" '"+name+"' cannot use character '"+ch+"'";
		}
	}
	
	return "";
}


/// Checks the name is valid
function check_name_input(name,te) 
{
	let warn = check_name_warn(name,te);
	if(warn != "") alert_import(warn);
}


/// Replaces charceters in a string
function char_replace(te)
{
	te = te.replace(/\*/g,"×");
	te = te.replace(/</g,"〈");
	return te.replace(/>/g,"〉");
}


/// Remove text escape characters
function remove_escape_char(te)
{
	let escape_char = [];
	for(let i = 0; i < greek_latex.length; i++){
		escape_char.push(["\\"+greek_latex[i][0],greek_latex[i][1]]);
	}
	escape_char.push(["\\sum","Σ"]);
	
	let i = 0;

	while(i < te.length){
		if(te.substr(i,1) == "\\"){
			let j = 0; 
			while(j < escape_char.length && te.substr(i,escape_char[j][0].length) != escape_char[j][0]) j++;
			if(j < escape_char.length){
				te = te.substr(0,i)+escape_char[j][1]+te.substr(i+escape_char[j][0].length);
			}
			else i++;
		}
		else i++;
	}
	
	return te;
}


/// Replaces greek letters with escape characters
function add_escape_char(te)
{
	for(let i = 0; i < greek_latex.length; i++){
		let st1 = "\\"+greek_latex[i][0];
		let st2 = greek_latex[i][1];
		let regex = new RegExp(st2,"g"); 
		te = te.replace(regex,st1);
	}
	
	te = te.replace("×","*");
	te = te.replace("→","->");
	
	return te;
}


/// Determines if a filter should be applied (there are differences between population and transition plots
function apply_filter(rpf2,i,name)
{
	let rpf3 = rpf2.filter[i];
	
	let ty;
	switch(name){
	case "Transitions": ty = "trans"; break;
	case "Populations": case "Individuals":  ty = "pop"; break;
	default: error("Name not recognised"+name); break;
	}
	
	if(ty == "pop" && rpf3.type == "trans_filt") return false;

	if(ty == "trans" && rpf3.type == "pop_filt" && rpf2.sel_class.cl == rpf3.cl) return false;
	return true;
}


/// Finds the root for a file name
function find_root(file)
{
	let i = file.length-1; 
	while(i >= 0 && file.substr(i,1) != "\\" && file.substr(i,1) != "/") i--;

	return file.substr(0,i+1);
}


/// Finds the local file for a file name
function find_file(file)
{
	let i = file.length-1; 
	while(i >= 0 && file.substr(i,1) != "\\" && file.substr(i,1) != "/") i--;
	return file.substr(i+1);
}


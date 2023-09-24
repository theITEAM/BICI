// Implements an import command

"use strict";

/// Adds a species to the model
function species_command(loop)
{
	let name = get_tag_value("name");
	if(name == ""){ cannot_find_tag(); return;}
	
	let type = get_tag_value("type").toLowerCase();
	if(type == ""){ cannot_find_tag(); return;}
	
	if(option_error("type",type,["individual","population"]) == true) return;

	if(loop == 0){
		let type_conv = type.substr(0,1).toUpperCase()+type.substr(1);
	
		model.add_species(name,type_conv);
	
		if(check_for_error() == true) return;	
		inter.imp.p = model.species.length-1;		
	}
	else{
		inter.imp.p = find(model.species,"name",name);
		if(inter.imp.p == undefined){ alert_import("Problem finding species '"+name+"'"); return;}
	}
	inter.imp.cl = undefined;	

	initialise_pages();
}


/// Adds a classification to the model
function classification_command(loop)
{
	let p = inter.imp.p;

	if(p == undefined){
		alert_import("A species must be defined before a classification can be added");
		return;
	}
		
	let sp = model.species[p];
	
	let name = get_tag_value("name");
	if(name == ""){ cannot_find_tag(); return;}
	
	let clone = get_tag_value("clone");
	if(clone != ""){ // Clones a classification from another species
		if(loop == 0){
			let p2 = find(model.species,"name",clone);
			if(p2 == undefined){
				alert_import("In 'clone' cannot find the species '"+clone+"'");
				return;
			}
			
			let sp2 = model.species[p2];
			
			let cl2 = find(sp2.cla,"name",name);
			
			if(cl2 == undefined){
				alert_import("In clone species '"+clone+"' cannot find the classification '"+name+"'");
				return;
			} 
				
			model.clone_class(p,p2,cl2)	
			sp.cla[sp.ncla-1].import_line = inter.imp.line;
		}
		else{
			inter.imp.cl = find(model.species[p].cla,"name",name);
			if(inter.imp.cl == undefined){ alert_import("Problem finding classification '"+name+"'"); return;}
		}
	}
	else{
		let index = get_tag_value("index");
		if(index == ""){ cannot_find_tag(); return;}
		if(index.length != 1){
			alert_import("The index '"+index+"' must be just a single character");
			return;
		}
		
		if(index == "t"){
			alert_import("The index 't' cannot be used because it is reserved for time variation");
			return;
		}	
		
		if(index == "a"){
			alert_import("The index 'a' cannot be used because it is reserved for age variation");
			return;
		}
		
		if(find_in(alphabet,index) == undefined){
			alert_import("Index '"+index+"' must be from the lower case alphabet");
		}
		
		let op = {};
		
		let coord = get_tag_value("coord");
		if(coord != ""){
			coord = coord.toLowerCase();
			if(coord != "latlng" && coord != "cartesian"){
				alert_import("'coord' must be either 'latlng' or 'cartesian'");
				return;
			}					
			op.coord = coord;
		}
		
		if(coord == "latlng"){
			let def_map = get_tag_value("map");
			if(def_map != ""){
				if(def_map != "load"){
					alert_import("'map' must be set to 'load'");
					return;
				}
				op.default_map = true;
			}
		}
		
		if(loop == 0){
			model.add_classification(p,name,index,op);
		
			sp.cla[sp.ncla-1].import_line = inter.imp.line;
		
			if(check_for_error() == true) return;	
			inter.imp.cl = sp.ncla-1;		
		}
		else{
			inter.imp.cl = find(model.species[p].cla,"name",name);
			if(inter.imp.cl == undefined){ alert_import("Problem finding classification '"+name+"'"); return;}
		}
	}
	
	initialise_pages();
}


/// Sets the classification which is being worked on
function set_command()
{
	let p = inter.imp.p;
	let sp = get_tag_value("species");
	if(sp != ""){
		inter.imp.p = find(model.species,"name",sp);
		if(inter.imp.p == undefined){
			alert_import("Cannot set the species '"+sp+"' as this does not exist");
			return;
		}
	}
	
	let name = get_tag_value("classification");
	
	if(name != ""){
		let p = inter.imp.p;
		if(p == undefined){
			alert_import("A species must be set before a 'classification' can be set");
			return;
		};
	
		inter.imp.cl = find(model.species[p].cla,"name",name);
		if(inter.imp.cl == undefined){
			alert_import("Cannot set the classification '"+name+"'");
			return;
		}
	}
}


/// Sets a camera in a classification
function camera_command()
{
	let claa = get_claa(); if(claa == undefined) return;

	let scale = get_tag_value("scale"); if(scale == ""){ cannot_find_tag(); return;}
	if(isNaN(scale)){ alert_import("'scale' must be a number"); return;}
	if(scale <= 0){ alert_import("'scale' must be positive"); return;}
	
	claa.camera.scale = scale;

	switch(claa.camera.coord){
	case "cartesian":
		let x = get_tag_value("x"); if(x == ""){ cannot_find_tag(); return;}
		let y = get_tag_value("y"); if(y == ""){ cannot_find_tag(); return;}
		
		x = Number(x)*import_scale_factor; if(isNaN(x)){ alert_import("'x' must be a number"); return;}
		
		y = Number(y)*import_scale_factor; if(isNaN(y)){ alert_import("'y' must be a number"); return;}
			
		claa.camera.x = x;
		claa.camera.y = y;
		break;
		
	case "latlng":
		let lat = get_tag_value("lat"); if(lat == ""){ cannot_find_tag(); return;}
		let lng = get_tag_value("lng"); if(lng == ""){ cannot_find_tag(); return;}
		if(check_latlng(lat,lng) != "success") return;
		
		let p = transform_latlng(lng,lat);
		claa.camera.x = p.x;
		claa.camera.y = p.y;
		break;
		
	default:
		alert_import("'coord' must be 'cartesian' or 'latlng'"); 
		return;
	}
}


/// Adds a compartment to the model
function compartment_command()
{
	let claa = get_claa(); if(claa == undefined) return;
	let p = inter.imp.p;
	let cl = inter.imp.cl;
	
	let name = get_tag_value("name"); if(name == ""){ cannot_find_tag(); return;}
	
	let cam = claa.camera;
	
	let color = get_tag_value("color"); 
	if(color == ""){ 
		let check_col;
		switch(name){
		case "S": check_col = DGREEN; break;
		case "E": check_col = ORANGE; break;
		case "I": check_col = RED; break;
		case "R": check_col = BLUE; break;
		case "D": check_col = BLACK; break;
		}
			
		if(check_col != undefined){
			if(find(claa.comp,"col",check_col) == undefined) color = check_col;
		}
		
		if(color == ""){
			let j = 0; 
			while(j < auto_color.length && find(claa.comp,"col",auto_color[j]) != undefined) j++;
			
			if(j == auto_color.length) j = Math.floor(Math.random()*auto_color.length);
				
			color = auto_color[j];
		}	
	}
	else{
		if(is_Color(color) == false) alert_import("'"+color+"' is not a valid color."); 
	}
	
	let fix = false;
	let fix_str = get_tag_value("fix").toLowerCase(); 
	if(fix_str != ""){
		if(fix_str == "true") fix = true;
		else{
			if(fix_str == "false") fix = false;
			else{
				alert_import("'fix' must be 'true' or 'false'"); 
			}
		}
	}

	switch(cam.coord){
	case "cartesian":
		{
			let x = get_tag_value("x");
			let y = get_tag_value("y");
			if(x == "" && y == ""){
				x = "auto"; y = "auto";
			}
			else{
				if(x == ""){ cannot_find_tag(); return;}
				if(y == ""){ cannot_find_tag(); return;}
				
				x = Number(x)*import_scale_factor; if(isNaN(x)){ alert_import("'x' must be a number"); return;}
				y = Number(y)*import_scale_factor; if(isNaN(y)){ alert_import("'y' must be a number"); return;}
			
				for(let i = 0; i < claa.ncomp; i++){
					if(claa.comp[i].x == x && claa.comp[i].y == y){		
						alert_import("Compartments '"+name+"' and '"+claa.comp[i].name+"' cannot have the same position"); return;
					}	
				}
			}
			
			model.add_compartment(name,p,cl,x,y,color,fix,false);
			claa.comp[claa.ncomp-1].import_line = inter.imp.line;
		}
		break;
	
	case "latlng":
		{
			let bound_file = get_tag_value("boundary");
			if(bound_file != ""){
				let claa = get_claa(); if(claa == undefined) return;

				let res = add_individual_compartment_boundary(name,inter.imp.p,inter.imp.cl,bound_file,color);
				if(res != "success"){
					if(inter.help.te == undefined) alert_import(res);
					return;
				}
			}
			else{
				let lat = get_tag_value("lat"); if(lat == ""){ cannot_find_tag(); return;}
				let lng	= get_tag_value("lng"); if(lng == ""){ cannot_find_tag(); return;}
				if(check_latlng(lat,lng) != "success") return;
				
				let pt = transform_latlng(lng,lat);
			
				for(let i = 0; i < claa.ncomp; i++){
					if(claa.comp[i].x == pt.x && claa.comp[i].y == pt.y){		
						alert_import("Compartments '"+name+"' and '"+claa.comp[i].name+"' cannot have the same position"); return;
					}	
				}
				
				model.add_latlng_compartment(name,inter.imp.p,inter.imp.cl,pt.x,pt.y,color,fix);		
			}
		}
		break;

	default: error("Option not recognised 50"); break;
	}

	if(check_for_error() == true) return;	
	
	let markov_branch = false;
	let markov_branch_str = get_tag_value("branch-prob"); 
	if(markov_branch_str != ""){
		if(markov_branch_str == "true") markov_branch = true;
		else{
			if(markov_branch_str == "false") markov_branch = false;
			else{
				alert_import("'branch-prob' must be 'true' or 'false'"); 
			}
		}
	}
	let co = claa.comp[claa.ncomp-1];

	co.markov_branch = markov_branch;
}


/// Sets the data directory
function datadir_command()
{
	let dir = get_tag_value("folder"); 
	
	if(dir == ""){ cannot_find_tag(); return;}

	if(check_char_allowed(dir,"<>\"|?*") != "success") return;
	
	const fs = require('fs');
	
	if(!fs.existsSync(dir)){ 
		alert_import("The directory '"+dir+"' does not exist"); return;
	}
	
	inter.imp.datadir = dir;
	inter.imp.outputdir = dir+"/Output";
}


/// Determines if a file exists or not
function check_file_exists(file)
{
	let fs = require('fs');	
	if(!fs.existsSync(file)){
		alert_import("The file '"+file+"' does not exist"); return;
	}
	return "success";
}


/// Checks for disallowed characters
function check_char_allowed(st,not_allowed)
{
	for(let i = 0; i < st.length; i++){
		let ch = st.substr(i,1);
		if(find_in(not_allowed,ch) != undefined){
			alert_import("In '"+st+"' the character '"+ch+"' is not allowed"); return;
		}
	}
	
	return "success";
}


/// Sets a description of the analysis
function description_command()
{
	let te = get_tag_value("text"); if(te == ""){ cannot_find_tag(); return;}

	model.description.te = te.replace(/\|/g,"\n");
}
	
	
/// Adds a transition / source / sink to the model
function transition_command(cname)
{
	let pos = dist_pos;	if(cname == "source") pos = exp_dist_pos;

	let claa = get_claa(); if(claa == undefined) return;
	let cam = claa.camera;

	let p = inter.imp.p;
	let cl = inter.imp.cl;
	
	let type = get_tag_value("type");
	if(type == "") type = "exp(rate)";
	
	let fr, to;
	if(cname == "source"){
		fr = "Source"; 
		let te = get_tag_value("from"); 
		if(te != ""){ alert_import("A source transition should not have a 'from' tag"); return;}
	}
	else{ 
		fr = get_tag_value("from");
		if(fr == ""){ cannot_find_tag(); return;}
	}
	
	if(cname == "sink"){
		to = "Sink"; 
		let te = get_tag_value("to"); 
		if(te != ""){ alert_import("A sink transition should not have a 'to' tag"); return;}
	}
	else{ 
		to = get_tag_value("to"); 
		if(to == ""){ cannot_find_tag(); return;}
	}
	
	let ci, cf;
	
	if(fr == "Source") ci = "Source";
	else{
		ci = find(claa.comp,"name",fr);
		if(ci == undefined){ alert_import("Cannot find compartment '"+fr+"'"); return;}
	}
	
	if(to == "Sink") cf = "Sink";
	else{
		cf = find(claa.comp,"name",to);
		if(cf == undefined){ alert_import("Cannot find compartment '"+to+"'"); return;}
	}
	
	for(let tr = 0; tr < claa.ntra; tr++){
		if(claa.tra[tr].i == ci && claa.tra[tr].f == cf){
			if(ci == "Source") alert_import("Source to '"+to+"' already exists");
			else{
				if(cf == "Sink") alert_import("Sink from '"+fr+"' already exists");
				else alert_import("Transition from '"+fr+"' to '"+to+"' already exists");
			}
			return;
		}
	}

	let midp=[];

	switch(cam.coord){
	case "cartesian": 
		{
			if(cname == "source"){ 
				let x = get_tag_value("x");
				let y = get_tag_value("y");
				if(x == "" && y == ""){
					x = "auto"; y = "auto";
				}
				else{
					if(x == "" || y == ""){
						alert_import("The 'x' and 'y' positions must be set for the source"); return;
					}

					x = Number(x)*import_scale_factor; 
					if(isNaN(x)){ alert_import("'x' must be a number"); return;}
					y = Number(y)*import_scale_factor; 
					if(isNaN(y)){ alert_import("'y' must be a number"); return;}
				}
				
				midp.push({x:x, y:y});
			}

			let xmid = get_tag_value("mid-x");
			let ymid = get_tag_value("mid-y");
			if(xmid != ""){
				var xlist = xmid.split(",");
				if(ymid == ""){ alert_import("Must set 'mid-x' and 'mid-y' points"); return;}
				var ylist = ymid.split(",");
				if(xlist.length != ylist.length){
					alert_import("Must have the same number of 'x-mid' and 'y-mid' points"); return;
				}

				for(let n = 0; n < xlist.length; n++){
					xlist[n] = Number(xlist[n])*import_scale_factor; 
					if(isNaN(xlist[n])){ alert_import("'mid-x' must be a number"); return;}
					ylist[n] = Number(ylist[n])*import_scale_factor; 
					if(isNaN(ylist[n])){ alert_import("'mid-y' must be a number"); return;}
				}
				
				for(let n = 0; n < xlist.length; n++){
					midp.push({x:xlist[n], y:ylist[n]});
				}
			}
			else{
				if(ymid != ""){ alert_import("Must set 'mid-x' and 'mid-y' points"); return;}
			}

			if(cname == "sink"){ 
				let x = get_tag_value("x");
				let y = get_tag_value("y");
				
				if(x == "" && y == ""){
					x = "auto"; y = "auto";
				}
				else{
					if(x == "" || y == ""){
						alert_import("The 'x' and 'y' positions must be set for the sink"); return;
					}

					x = Number(x)*import_scale_factor; 
					if(isNaN(x)){ alert_import("'x' must be a number"); return;}
					y = Number(y)*import_scale_factor; 
					if(isNaN(y)){ alert_import("'y' must be a number"); return;}
				}
				
				midp.push({x:x, y:y});
			}
		}
		break;

	case "latlng":
		{
			if(cname == "source"){ 
				let lng = get_tag_value("lng");
				let lat = get_tag_value("lat");
				if(check_latlng(lat,lng) != "success") return;
				
				let p = transform_latlng(lng,lat);
				midp.push(p);
			}

			let lngmid = get_tag_value("mid-lng");
			let latmid = get_tag_value("mid-lat");
			if(lngmid != ""){
				var lnglist = lngmid.split(",");
				if(latmid == ""){ 
					alert_import("Must set 'mid-lat' and 'mid-lng' points"); return;
				}

				var latlist = latmid.split(",");
				if(lnglist.length != latlist.length){ 
					alert_import("Must have the same number of 'lat-mid' and 'lng-mid' points"); return;
				}

				for(let n = 0; n < lnglist.length; n++){
					if(check_latlng(latlist[n],lnglist[n]) != "success") return;
					let p = transform_latlng(lnglist[n],latlist[n]);
					midp.push(p);
				}
			}
			else{
				if(latmid != ""){ alert_import("Must set 'mid-lat' and 'mid-lng' points"); return;}
			}

			if(cname == "sink"){
				let lng = get_tag_value("lng");
				let lat = get_tag_value("lat");
				if(check_latlng(lat,lng) != "success") return;
				
				let p = transform_latlng(lng,lat);
				midp.push(p);
			}
		}
		break;

	default: error("Option not recognised 130"); break;
	}

	type = correct_caps(type,pos);

	if(option_error("type",type,pos) == true) return;

	model.add_transition(p,cl,ci,cf,midp,type);

	if(check_for_error() == true) return;	
	
	let traa = claa.tra[claa.ntra-1];
	traa.import_line = inter.imp.line;
	
	traa.bp_set = false;
	
	let val = traa.value;
	
	switch(type){
	case "exp(rate)":
		{
			let rate = get_tag_value("rate");
			if(rate != ""){
				if(import_eqn_value(val.rate_eqn,rate) != "success") return;
			}
		}
		break;
		
	case "exp(mean)":
		{
			let mean = get_tag_value("mean"); 
			if(mean != ""){ 
				if(import_eqn_value(val.mean_eqn,mean) != "success") return;
			}
		}
		break;
		
	case "gamma":
		{
			let mean = get_tag_value("mean"); 
			if(mean != ""){ 
				if(import_eqn_value(val.mean_eqn,mean) != "success") return;
			}
			
			let cv = get_tag_value("cv"); 
			if(cv != ""){ 
				if(import_eqn_value(val.cv_eqn,cv) != "success") return;
			}
		}
		break;
		
	case "erlang":
		{
			let mean = get_tag_value("mean"); 
			if(mean != ""){ 
				if(import_eqn_value(val.mean_eqn,mean) != "success") return;
			}
			
			let shape = get_tag_value("shape"); 
			if(shape == ""){ cannot_find_tag(); return;}
			
			let num = Number(shape);
			if(isNaN(shape) || num <= 0 || Math.round(num) != num){
				alert_import("For an Erlang distribution the shape parameter must be a positive integer");				
				return;
			}
			
			if(import_eqn_value(val.shape_erlang,shape) != "success") return;
		}
		break;
		
	case "log-normal":
		{
			let mean = get_tag_value("mean"); 
			if(mean != ""){ 
				if(import_eqn_value(val.mean_eqn,mean) != "success") return;
			}
			
			let cv = get_tag_value("cv"); 
			if(cv != ""){ 
				if(import_eqn_value(val.cv_eqn,cv) != "success") return;
			}
		}
		break;
		
	case "weibull":
		{
			let scale = get_tag_value("scale"); 
			if(scale != ""){ 
				if(import_eqn_value(val.scale_eqn,scale) != "success") return;
			}
			
			let shape = get_tag_value("shape"); 
			if(shape != ""){ 
				if(import_eqn_value(val.shape_eqn,shape) != "success") return;
			}
		}
		break;

	default: error("Option not recognised 51"); break;
	}
	
	let bp = get_tag_value("bp");
	if(bp != ""){
		traa.bp_set = true;
		import_eqn_value(val.bp_eqn,bp);
	}
	
	model.update_pline(p,cl);
}


/// Adds a label to the model
function label_command()
{
	let claa = get_claa(); if(claa == undefined) return;
	
	let te = get_tag_value("text"); if(te == ""){ cannot_find_tag(); return;}
	let tesize = get_tag_value("textsize"); if(tesize == ""){ cannot_find_tag(); return;}
	tesize = Number(tesize); if(isNaN(tesize)){ alert_import("'textsize' must be a number"); return;}
	
	let x = get_tag_value("x"); if(x == ""){ cannot_find_tag(); return;}
	let y = get_tag_value("y"); if(y == ""){ cannot_find_tag(); return;}
	x = Number(x)*import_scale_factor; if(isNaN(x)){ alert_import("'x' must be a number"); return;}
	y = Number(y)*import_scale_factor; if(isNaN(y)){ alert_import("'y' must be a number"); return;}
	let color = get_tag_value("color"); if(color == ""){ cannot_find_tag(); return;}

	claa.annotation.push({type:"text", te:te, tesize:tesize, x:x, y:y, color:color});
}


/// Adds a box around compartments in the model
function box_command()
{
	let claa = get_claa(); if(claa == undefined) return;
	
	let te = get_tag_value("text"); if(te == ""){ cannot_find_tag(); return;}
	let tesize = get_tag_value("textsize"); if(tesize == ""){ cannot_find_tag(); return;}
	tesize = Number(tesize); if(isNaN(tesize)){ alert_import("'textsize' must be a number"); return;}
	
	let comps = get_tag_value("compartments"); if(comps == ""){ cannot_find_tag(); return;}
	let color = get_tag_value("color"); if(color == ""){ cannot_find_tag(); return;}
			
	let list=[];
	let spl = comps.split(",");
	if(spl.length == 0){ alert_import("Error 'box' contains no compartments."); return;}
	for(let k = 0; k < spl.length; k++){
		let c = 0; while(c < claa.ncomp && claa.comp[c].name != spl[k]) c++;
		if(c == claa.ncomp){ 
			alert_import("Cannot find compartment '"+spl[k]+"'."); return;
		}
		list.push(spl[k]);
	}
	
	claa.annotation.push({type:"box", te:te, tesize:tesize, comps:list, color:color});
}


/// Loads data inito the model
function data_command()
{
	load_data();
	//if(inter.help.te != undefined) break;
}
		
		
/// Sets the value for a parameter in the model
function param_command()
{
	let full_name = get_tag_value("name"); if(full_name == ""){ cannot_find_tag(); return;}
	
	let pp = get_param_prop(full_name);

	let par = create_new_param(pp,"normal");

	par.import_line = inter.imp.line;
	
	if(par.time_dep == true || par.age_dep == true){
		let knot_times = get_tag_value("knot-times"); if(knot_times == ""){ cannot_find_tag(); return;}
		let warn = check_knot_times(knot_times);
		if(warn != ""){ alert_import(warn); return;}
		
		change_spline_knot(par,knot_times.split(","));
	
		let smooth = get_tag_value("smooth").toLowerCase().trim();
		
		par.spline.smooth = default_smooth();
		
		if(smooth != ""){
			par.spline.smooth.check = true;
			
			let spl = smooth.split("(");
			
			if(option_error("'smooth'",spl[0],["normal","log-normal"]) == true) return;
		
			switch(spl[0]){
			case "normal": par.spline.smooth.type = {value:"Normal"}; break;
			case "log-normal": par.spline.smooth.type = {value:"Log-Normal"}; break;
			}
	
			if(spl.length != 2){ alert_import("There is syntax error in 'smooth'"); return;}
			
			if(spl[1].substr(spl[1].length-1,1) != ")"){ alert_import("There is syntax error in 'smooth'"); return;}
			let val = spl[1].substr(0,spl[1].length-1);
			if(isNaN(val)){ alert_import("In 'smooth' the value '"+val+"' is not a number"); return;}
			if(Number(val) <= 0){ alert_import("In 'smooth' the value '"+val+"' is not positive"); return;}
			par.spline.smooth.value = val;
		}
	}
	
	let cons = get_tag_value("constant"); 
	let dist = get_tag_value("dist").toLowerCase();
	let dist_split = get_tag_value("dist-split").toLowerCase();
	let value = get_tag_value("value");
	let reparam = get_tag_value("reparam"); 
	let prior = get_tag_value("prior").toLowerCase(); 
	let prior_split = get_tag_value("prior-split").toLowerCase(); 
	let dist_mat = get_tag_value("distance-matrix").toLowerCase(); 

	let param_tag = [];
	param_tag.push({val:cons, tag:"constant"});
	param_tag.push({val:dist, tag:"dist"});
	param_tag.push({val:dist_split, tag:"dist-split"});
	param_tag.push({val:value, tag:"value"});
	param_tag.push({val:reparam, tag:"reparam"});
	param_tag.push({val:prior, tag:"prior"});
	param_tag.push({val:prior_split, tag:"prior-split"});
	param_tag.push({val:dist_mat, tag:"distance-matrix"});

	for(let j = 0; j < param_tag.length; j++){
		for(let i = j+1; i < param_tag.length; i++){
			if(param_tag[j].val != "" && param_tag[i].val != ""){
				alert_import("'"+param_tag[j].tag+"' and '"+param_tag[i].tag+"' cannot both be set"); 
				return;
			}
		}
	}
	
	let j = 0; while(j < param_tag.length && param_tag[j].val == "") j++;
	
	if(j == param_tag.length){
		let te = "One of these possibilities must be set: ";
		for(let j = 0; j < param_tag.length; j++){
			if(j != 0) te += ", ";
			te += "'"+param_tag[j].tag+"'";
		}
		alert_import(te); 
		return;
	}
		
	par.dist_matrix = { check:false}; 
	if(dist_mat != ""){
		if(dist_mat != "true") alert_import("'distance-matrix' must be set to 'true'");
		par.dist_matrix.check = true;
	}
	
	par.auto_value = false;
	
	if(value == "auto"){
		par.auto_value = true;
		par.set = true;
	}
	else{
		if(par.dep.length == 0){
			if(value != ""){
				par.value = value;
				par.variety = "normal";
			}
			
			if(reparam != ""){
				par.value = reparam;
				par.variety = "reparam";
			}
			
			if(cons != ""){
				par.value = cons;
				par.variety = "const";
			}
		}
		else{
			if(value != "" || cons != "" || reparam != ""){
				let desc = "In 'value'";
				let valu = value; 
				if(valu == ""){
					if(cons != ""){
						valu = cons; desc = "In 'constant'";
						par.variety = "const";		
					}
					else{
						if(reparam != ""){
							valu = reparam; desc = "In 'reparam'";
							par.variety = "reparam";		
						}
						else error("Problem importing");
					}
				}
				
				let tab = load_table_from_file(valu);
				if(tab == undefined) return;
				if(typeof tab == 'string'){ alert_import(tab); return;}
			
				let col_name = copy(par.dep);
				col_name.push("Value");
				
				let subtab = get_subtable(tab,col_name);
				if(subtab.error != ""){	alert_import(subtab.error); return;}
			
				let ncol = subtab.ncol;
			
				for(let r = 0; r < subtab.nrow; r++){
					let ind = [];
					for(let i = 0; i < ncol-1; i++){
						ind[i] = find_in(par.list[i],subtab.ele[r][i]);
						if(ind[i] == undefined){ 
							alert_import(desc+" the element '"+subtab.ele[r][i]+"' is not valid (column '"+subtab.heading[i]+"', row "+(r+2)+")");
							return;
						}
					}
					let ele = subtab.ele[r][ncol-1];
					let val;
					
					switch(par.variety){
					case "normal": case "const":
						val = Number(ele);
					
						if(isNaN(val)){
							alert_import(desc+" the element '"+ele+"' is not a number (column '"+subtab.heading[ncol-1]+"', row "+(r+2)+")");
							return;
						}
						break;
					
					case "reparam":
						if(!isNaN(ele)) val = Number(ele);
						else{
							val = ele;
							if(is_eqn(val,"Table element",{}) == false){
								alert_import(desc+" the element '"+ele+"' is not a valid equation (column '"+subtab.heading[ncol-1]+"', row "+(r+2)+")");
								return;
							}
						}
						break;
						
					default: error("option not recognised"); break;
					}
					
					set_element(par.value,ind,val);
				}
				par.set = true;
			}
		}
	}
	
	if(prior != ""){
		let pri = par.prior;
			
		let val;
		
		pri.type.te = prior;
			
		switch(prior){
		case "fix":
			{
				let mean = get_tag_value("fixed"); if(mean == ""){ cannot_find_tag(); return;}
				if(!is_number(mean,"fixed")) return;
				pri.value.mean_eqn.te = mean;
			}
			break;
			
		case "uniform":
			{	
				let min = get_tag_value("min"); if(min == ""){ cannot_find_tag(); return;}
				if(!is_number(min,"min")) return;
				pri.value.min_eqn.te = min;
				
				let max = get_tag_value("max"); if(max == ""){ cannot_find_tag(); return;}
				if(!is_number(max,"max")) return;
				pri.value.max_eqn.te = max;
				
				if(Number(min) >= Number(max)){
					alert_import("'min' must be smaller than 'max'"); 
					return;
				}
			}
			break;
			
		case "normal":
			{
				let mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				if(!is_number(mean,"mean")) return;
				pri.value.mean_eqn.te = mean;
					
				let sd = get_tag_value("sd"); if(sd == ""){ cannot_find_tag(); return;}
				if(!is_positive(sd,"sd")) return;
				pri.value.sd_eqn.te = sd;
			}
			break;
			
		case "log-normal":
			{
				let mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				if(!is_positive(mean,"mean")) return;
				pri.value.mean_eqn.te = mean;
					
				let cv = get_tag_value("cv"); if(cv == ""){ cannot_find_tag(); return;}
				if(!is_positive(cv,"cv")) return;
				pri.value.cv_eqn.te = cv;
			}
			break;
			
		case "exp":
			{
				let mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				if(!is_positive(mean,"mean")) return;
				pri.value.mean_eqn.te = mean;
			}
			break;
			
		case "gamma":
			{
				let mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				if(!is_positive(mean,"mean")) return;
				pri.value.mean_eqn.te = mean;
				
				let cv = get_tag_value("cv"); if(cv == ""){ cannot_find_tag(); return;}
				if(!is_positive(cv,"cv")) return;
				pri.value.cv_eqn.te = cv;
			}
			break;
		
		case "beta":
			{
				let alpha = get_tag_value("alpha"); if(alpha == ""){ cannot_find_tag(); return;}
				if(!is_positive(alpha,"alpha")) return;
				pri.value.alpha_eqn.te = alpha;
				
				let beta = get_tag_value("beta"); if(beta == ""){ cannot_find_tag(); return;}
				if(!is_positive(beta,"beta")) return;
				pri.value.beta_eqn.te = beta;
			}
			break;
			
		case "bernoulli":
			{
				let mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				if(!is_zeroone(mean,"mean")) return;
				pri.value.mean_eqn.te = mean;
			}
			break;
			
		case "flat":
			break;
		
		case "dirichlet":
			let alpha = get_tag_value("alpha"); if(alpha == ""){ cannot_find_tag(); return;}
			if(!is_positive(alpha,"alpha")) return;
			pri.value.alpha_eqn.te = alpha;
			break;
	
		default: 
			error("Prior type '"+prior+"' not recognised");
			break;
		}
	}
	
	if(prior_split != ""){
		par.prior_split_check = {check:true};
		par.prior_split_set = true;
		if(par.dep.length == 0){
			alert_import("'prior-split' can only be used if the parameter has a dependency.");  
		}
	
		let tab = load_table_from_file(prior_split);
		if(tab == undefined) return;	
		if(typeof tab == 'string'){ alert_import(tab); return;}
		
		let col_name = copy(par.dep);
		col_name.push("Prior");
		
		let subtab = get_subtable(tab,col_name);
		if(subtab.error != ""){ alert_import(subtab.error); return;}
	
		let ncol = subtab.ncol;
		
		for(let r = 0; r < subtab.nrow; r++){
			let ind = [];
			for(let i = 0; i < ncol-1; i++){
				ind[i] = find_in(par.list[i],subtab.ele[r][i]);
				if(ind[i] == undefined){ 
					alert_import("The table element '"+subtab.ele[r][i]+"' is not valid (column '"+subtab.heading[i]+"', row "+(r+2)+")");
					return;
				}
			}
			
			let pri = convert_text_to_prior(subtab.ele[r][ncol-1],par.pri_pos)
			if(typeof pri == 'string'){
				alert_import("The table element '"+subtab.ele[r][ncol-1]+"' is not valid (col '"+subtab.heading[ncol-1]+"', row "+(r+2)+"). "+pri);
				return;
			}
			
			set_element(par.prior_split,ind,pri);
		}
	}
	
	if(dist != ""){
		let pri = par.prior;
		par.variety = "dist";
		
		pri.type.te = dist;
				
		if(option_error("type",dist,prior_pos) == true) return;
		
		switch(dist){
		case "fix":
			{
				let mean = get_tag_value("value"); if(mean == ""){ cannot_find_tag(); return;}
				if(!is_eqn(mean,"value",{})) return;
				pri.value.mean_eqn.te = mean;
			}
			break;
			
		case "uniform":
			{	
				let min = get_tag_value("min"); if(min == ""){ cannot_find_tag(); return;}
				if(!is_eqn(min,"min",{})) return;
				pri.value.min_eqn.te = min;
				
				let max = get_tag_value("max"); if(max == ""){ cannot_find_tag(); return;}
				if(!is_eqn(max,"max",{})) return;
				pri.value.max_eqn.te = max;
			}
			break;
			
		case "normal":
			{
				let mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				if(!is_eqn(mean,"mean",{})) return;
				pri.value.mean_eqn.te = mean;
					
				let sd = get_tag_value("sd"); if(sd == ""){ cannot_find_tag(); return;}
				if(!is_eqn(sd,"sd",{positive:true})) return;
				pri.value.sd_eqn.te = sd;
			}
			break;
			
		case "log-normal":
			{
				let mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				if(!is_eqn(mean,"mean",{positive:true})) return;
				pri.value.mean_eqn.te = mean;
					
				let sd = get_tag_value("sd"); if(sd == ""){ cannot_find_tag(); return;}
				if(!is_eqn(sd,"sd",{positive:true})) return;
				pri.value.sd_eqn.te = sd;
			}
			break;
			
		case "exp":
			{
				let mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				if(!is_eqn(mean,"mean",{positive:true})) return;
				pri.value.mean_eqn.te = mean;
			}
			break;
			
		case "gamma":
			{
				let mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				if(!is_eqn(mean,"mean",{positive:true})) return;
				pri.value.mean_eqn.te = mean;
				
				let cv = get_tag_value("cv"); if(cv == ""){ cannot_find_tag(); return;}
				if(!is_eqn(cv,"cv",{positive:true})) return;
				pri.value.cv_eqn.te = cv;
			}
			break;
		
		case "beta":
			{
				let alpha = get_tag_value("alpha"); if(alpha == ""){ cannot_find_tag(); return;}
				if(!is_eqn(alpha,"alpha",{positive:true})) return;
				pri.value.alpha_eqn.te = alpha;
				
				let beta = get_tag_value("beta"); if(beta == ""){ cannot_find_tag(); return;}
				if(!is_eqn(beta,"beta",{positive:true})) return;
				pri.value.beta_eqn.te = beta;
			}
			break;
			
		case "bernoulli":
			{
				let mean = get_tag_value("mean"); if(mean == ""){ cannot_find_tag(); return;}
				if(!is_eqn(mean,"mean",{zero_one_range:true})) return;
				pri.value.mean_eqn.te = mean;
			}
			break;
			
		case "flat":
			break;
		
		case "dirichlet":
			let alpha = get_tag_value("alpha"); if(alpha == ""){ cannot_find_tag(); return;}
			if(!is_eqn(alpha,"alpha",{positive:true})) return;
			pri.value.alpha_eqn.te = alpha;
			break;
		
		default: 
			error("Distribution type '"+dist+"' not recognised");
			break;
		}
	}
	
			
	if(dist_split != ""){
		par.prior_split_check = {check:true};
		par.variety = "dist";
		par.prior_split_set = true;
		
		par.prior_split_check = {check:true};
		if(par.dep.length == 0){
			alert_import("'dist-split' can only be used if the parameter has a dependency.");  
		}
		
		let tab = load_table_from_file(dist_split);
		
		if(tab == undefined) return;	
		if(typeof tab == 'string'){ alert_import(tab); return;}
		
		let col_name = copy(par.dep);
		col_name.push("Dist");
		
		let subtab = get_subtable(tab,col_name);
		if(subtab.error != ""){ alert_import(subtab.error); return;}
	
		let ncol = subtab.ncol;
		
		for(let r = 0; r < subtab.nrow; r++){
			let ind = [];
			for(let i = 0; i < ncol-1; i++){
				ind[i] = find_in(par.list[i],subtab.ele[r][i]);
				if(ind[i] == undefined){ 
					alert_import("The table element '"+subtab.ele[r][i]+"' is not valid (column '"+subtab.heading[i]+"', row "+(r+2)+")");
					return;
				}
			}
			
			let pri = convert_text_to_prior(subtab.ele[r][ncol-1],par.pri_pos,true);
			
			if(typeof pri == 'string'){
				alert_import("The table element '"+subtab.ele[r][ncol-1]+"' is not valid (col '"+subtab.heading[ncol-1]+"', row "+(r+2)+"). "+pri);
				return;
			}
			
			set_element(par.prior_split,ind,pri);
		}
	}
	
	let sim_sample = get_tag_value("sim-sample").toLowerCase(); 
	if(sim_sample == "") sim_sample = "true";
	else{
		if(dist == ""){
			alert_import("'sim-sample' can only be set if a distribution is set through 'dist'");
			return;
		}
	}
	if(option_error("sim-sample",sim_sample,["true","false"]) == true) return;
	
	if(dist != ""){
		switch(sim_sample){
		case "false": par.sim_sample = {check:false}; break;
		case "true": par.sim_sample = {check:true}; break;
		}
	}
	
	par.full_name = param_name(par);

	model.param.push(par);
}


/// Definition for derived quantities
function derived_command()
{
	let full_name = get_tag_value("name"); if(full_name == ""){ cannot_find_tag(); return;}
	let eqn_name = get_tag_value("eqn"); if(eqn_name == ""){ cannot_find_tag(); return;}
	
	let eqn1 = create_equation(full_name,"derive_param");
	if(eqn1.warn.length > 0){ alert_import("In 'name': "+eqn1.warn[0].te); return;}
	if(eqn1.param.length != 1){
		pr("PP");
		pr(eqn1);
		pr(full_name);
		alert_import("'name' must contain a single paramter"); return;
	}

	let eqn2 = create_equation(eqn_name,"derive_eqn");
	if(eqn2.warn.length > 0){ alert_import("In 'eqn': "+eqn2.warn[0].te); return;}
	
	model.derive.push({eqn1:eqn1, eqn2:eqn2});
}


/// Applies the simulation command
function simulation_command()
{
	let details = model.sim_details;
	
	let start = get_tag_value("start"); if(start == ""){ cannot_find_tag(); return;}
	if(isNaN(start)){ alert_import("'start' must be a number"); return;}
	details.t_start = start;
	
	let end = get_tag_value("end"); if(end == ""){ cannot_find_tag(); return;}
	if(isNaN(end)){ alert_import("'end' must be a number"); return;}
	details.t_end = end;
	
	if(Number(details.t_start) >= Number(details.t_end)){
		alert_import("'start' must before 'end'"); return;
	}
	
	let num = get_tag_value("number").toLowerCase();
	if(num == "") num = 1;
	let numf = Number(num);
	
	if(isNaN(num) || numf != Math.floor(numf) || numf <= 0){ 
		alert_import("'number' must be a positive integer"); 
		return;
	}
	details.number = num;
	
	let alg = get_tag_value("algorithm").toLowerCase();
	if(alg == "") alg = "gillespie";
	
	if(option_error("algorithm",alg,sim_alg_list) == true) return;
	details.algorithm.value = alg;		
	
	let dt_str = get_tag_value("timestep"); if(dt_str == ""){ cannot_find_tag(); return;}
	let dt = Number(dt_str);
	if(isNaN(dt_str) || dt <= 0){
		alert_import("'timestep' must be a positive number"); return;
	}
	details.timestep = dt;
	
	details.indmax = check_pos_integer("ind-max",10000);
}


/// Applies the simulation command
function inference_command()
{
	let details = model.inf_details;
	
	let start = get_tag_value("start"); if(start == ""){ cannot_find_tag(); return;}
	if(isNaN(start)){ alert_import("'start' must be a number"); return;}
	details.t_start = start;
	
	let end = get_tag_value("end"); if(end == ""){ cannot_find_tag(); return;}
	if(isNaN(end)){ alert_import("'end' must be a number"); return;}
	details.t_end = end;
	
	if(Number(details.t_start) >= Number(details.t_end)){
		alert_import("'start' must before 'end'"); return;
	}
	
	let alg = get_tag_value("algorithm").toUpperCase();
	if(alg == "") alg = "DA-MCMC";
	
	if(option_error("algorithm",alg,inf_alg_list) == true) return;
	
	details.algorithm.value = alg;	

	let dt_str = get_tag_value("timestep"); if(dt_str == ""){ cannot_find_tag(); return;}
	let dt = Number(dt_str);
	if(isNaN(dt_str) || dt <= 0){
		alert_import("'timestep' must be a positive number"); return;
	}
	details.timestep = dt;
	
	details.indmax = check_pos_integer("ind-max",10000);

	switch(alg){
	case "DA-MCMC":
		details.sample = check_pos_integer("sample",10000);
		details.thinparam = check_pos_integer("param-thin",10);
		details.thinstate = check_pos_integer("state-thin",100);
		break;
		
	case "ABC":
		details.abcsample = check_pos_integer("sample",1000);
		details.accfrac = check_zero_one("acc-frac",0.1);
		break;
		
	case "ABC-SMC":
		details.abcsample = check_pos_integer("sample",1000);
		details.accfrac = check_zero_one("acc-frac",0.5);
		details.numgen = check_pos_integer("gen",5);
		details.kernelsize = check_pos_number("kernel-size",0.5);
		break;
		
	default: break;
	
		break;
	}
	
}


/// Check that a number is a positive integer
function check_pos_integer(te,def)
{
	let value = get_tag_value(te);
	
	if(value == "") return def;
	
	let num = Number(value);
	if(isNaN(value) || num <= 0 || num != Math.floor(num)){
		alert_import("'"+te+"' must be a positive integer");
	}
	return num;
}

/// Check that a number is positive 
function check_pos_number(te,def)
{
	let value = get_tag_value(te);
	
	if(value == "") return def;
	
	let num = Number(value);
	if(isNaN(value) || num <= 0){
		alert_import("'"+te+"' must be a positive number");
	}
	return num;
}


/// Check that a number is between zero and one
function check_zero_one(te,def)
{
	let value = get_tag_value(te);
	
	if(value == "") return def;
	
	let num = Number(value);
	if(isNaN(value) || num <= 0 || num > 1){
		alert_import("'"+te+"' must be between zero and one");
	}
	return num;
}


/// Loads individual effects into the model
function ind_effect_command()
{
	let p = inter.imp.p;

	if(p == undefined){
		alert_import("A species must be defined before an individual effect can be added");
		return;
	}
		
	let name = get_tag_value("name"); if(name == ""){ cannot_find_tag(); return;}
	let spl = name.split(",");
	let list = []; for(let i = 0; i < spl.length; i++) list.push({name:spl[i]});
	
	let A = get_tag_value("A"); 
	
	let A_matrix;
	if(A == ""){
		A_matrix = {check:false, loaded:false, value:[], ind_list:[]};
	}
	else{
		A_matrix = {check:true, loaded:true};
		
		let tab = load_table_from_file(A);
		if(tab == undefined) return;
		if(typeof tab == 'string'){ alert_import(tab); return;}

		A_matrix.ind_list = tab.heading;
		if(tab.nrow != tab.ncol){
			alert_import("The file '"+tab.filename+"' must contain an equal number of columns and rows."); 
			return;
		}
		
		let val = [];
		for(let r = 0; r < tab.nrow; r++){
			val[r] = [];
			for(let c = 0; c < tab.ncol; c++){
				let ele = tab.ele[r][c];
				if(isNaN(ele)){
					alert_import("In file '"+tab.filename+"' the element '"+ele+"' is not a number (row "+(r+2)+", col "+(c+1)+")");
					return;
				}
				val[r][c] = Number(ele);
			}
		}
		A_matrix.value = val;
	}
	
	model.species[p].ind_eff_group.push({list:list, A_matrix:A_matrix});
}


/// Loads fixed effects into the model
function fixed_effect_command()
{
	let p = inter.imp.p;

	if(p == undefined){
		alert_import("A species must be defined before a fixed effect can be added");
		return;
	}
		
	let name = get_tag_value("name"); if(name == ""){ cannot_find_tag(); return;}

	let X = get_tag_value("X"); if(X == ""){ cannot_find_tag(); return;}
	
	let X_vector = { loaded:true, ind_list:[], value:[]};
		
	let tab = load_table_from_file(X);
	if(tab == undefined) return;
	
	let col_name = ["ID","value"]; 
	let subtab = get_subtable(tab,col_name);
	if(subtab.error != ""){	alert_import(subtab.error); return;}
	
	for(let r = 0; r < subtab.nrow; r++){
		X_vector.ind_list.push(subtab.ele[r][0]);
		let val = tab.ele[r][1];
		if(isNaN(val)){
			alert_import("In file '"+tab.filename+"' the element '"+val+"' is not a number (row "+(r+2)+")");
			return;
		}
		X_vector.value.push(val);
	}

	model.species[p].fix_eff.push({name:name, X_vector:X_vector});
}


/// Imports a data-table 
function import_data_table_command(cname)
{
	let siminf = "inf";
	switch(cname){
	case "add-ind-sim": cname = "add-ind"; siminf = "sim"; break;
	case "remove-ind-sim": cname = "remove-ind"; siminf = "sim"; break;
	case "move-ind-sim": cname = "move-ind"; siminf = "sim"; break;
	case "init-pop-sim": cname = "init-pop"; siminf = "sim"; break;
	}
	
	let i = find(convert,"command",cname);
	if(i == undefined){
		error("Counld not not import data table:"+cname);
		return;
	}
	
	let type = convert[i].type;
	
	let file = get_tag_value("file");
	let tab; 
	if(file != ""){
		tab = load_table_from_file(file);
		if(tab == undefined) return;
		if(typeof tab == 'string'){ alert_import(tab); return;}
	}
	
	let cols = get_tag_value("cols");
	let cols_split = cols.split(",");
	
	let p = inter.imp.p;
	if(p == undefined){
		alert_import("To load the data file the species must be set"); return;
	}
	
	let sp = model.species[p];
	
	let spec;
	let info = {p:p, siminf:siminf};
	
	switch(type){  // Modifies specification dependent on data source
	case "Init. Pop.":
		{
			let focal = get_tag_value("focal");
		
			if(focal == "") spec = {radio:{value:"File"},radio2:{value:"All"}, focal:{te:""}};
			else{
				let cl = find_nocase(sp.cla,"name",focal);
				if(cl == undefined){ alert_import("'"+focal+"' is not a classification name"); return;}
				else focal = sp.cla[cl].name;
		
				spec = {radio:{value:"File"},radio2:{value:"Focal"}, focal:{te:focal}};
			}
		}	
		break;
		
	case "Add Ind.":
		break;
		
	case "Remove Ind.":
		break;
	
	case "Move Ind.":
		{
			let name = get_tag_value("class");
			if(name == ""){ cannot_find_tag(); return;}
		
			let cl = find_nocase(sp.cla,"name",name);
			if(cl == undefined){
				alert_import("In 'class' the value '"+name+"' is not a classification"); 
				return;
			}
			
			spec = {cl_drop:{te:sp.cla[cl].name}};
		}
		break;
	
	case "Init. Pop. Prior": 
		{
			let ty = get_tag_value("type").toLowerCase();
	
			if(option_error("type",ty,["flat","dirichlet"]) == true) return;
	
			switch(ty){
			case "flat": spec = {radio:{value:"Flat"}}; break;
			case "dirichlet": spec = {radio:{value:"Dirichlet"}}; break;
			default: error("Op not recognised"); break;
			}
		}
		break;
	
	case "Compartment":
		{
			let name = get_tag_value("class");
			if(name == ""){ cannot_find_tag(); return;}
			
			let cl = find_nocase(sp.cla,"name",name);
			if(cl == undefined){
				alert_import("In 'class' the value '"+name+"' is not a classification"); 
				return;
			}
		
			spec = {cl_drop:{te:sp.cla[cl].name}};
		}	
		break;
		
		
	case "Transition": case "Source": case "Sink":
		{
			switch(type){
			case "Transition":
				{
					let from = get_tag_value("from");
					if(from == ""){ cannot_find_tag(); return;}
			
					let cl_fr = get_cl_from_comp(from,p);
					if(cl_fr == undefined){
						alert_import("'"+from+"' is not a valid compartment"); 
						return;
					}
				
					let to = get_tag_value("to");
					if(to == ""){ cannot_find_tag(); return;}
			
					let cl_to = get_cl_from_comp(to,p);
					if(cl_to == undefined){ 
						alert_import("'"+to+"' is not a valid compartment"); 
						return;
					}
				
					if(cl_fr != cl_to){ 
						alert_import("'"+from+"' and '"+to+"' cannot be in different classifications");
						return;
					}
						
					spec = {tr_drop:{te:from+"→"+to}};				
				}
				break;
				
			case "Source":
				{
					let to = get_tag_value("to");
					if(to == ""){ cannot_find_tag(); return;}
			
					let cl_to = get_cl_from_comp(to,p);
					if(cl_to == undefined){ 
						alert_import("'"+to+"' is not a valid compartment"); 
						return;
					}
					
					spec = {tr_drop:{te:to}};
				}
				break;
				
			case "Sink":
				{
					let from = get_tag_value("from");
					if(from == ""){ cannot_find_tag(); return;}

					let cl_fr = get_cl_from_comp(from,p);
					if(cl_fr == undefined){
						alert_import("'"+from+"' is not a valid compartment"); 
						return;
					}
					
					spec = {tr_drop:{te:from}};
				}
				break;
			}
		
			let or = get_tag_value("obsrange").toLowerCase();
			if(or == ""){ cannot_find_tag(); return;}
			
			if(option_error("obsrange",or,["all","specify","file"]) == true) return;
			
			spec.time_start = "";
			spec.time_end = "";
			
			switch(or){
			case "all": 
				spec.time_radio = {value:"All"};
				break;
				
			case "specify":
				spec.time_radio = {value:"Spec"}; 
				
				let start = get_tag_value("start");
				if(start == ""){ cannot_find_tag(); return;}
				if(isNaN(start)){ alert_import("'start' must be a number"); return;}
				spec.time_start = start;
				
				let end = get_tag_value("end");
				if(end == ""){ cannot_find_tag(); return;}
				if(isNaN(end)){ alert_import("'end' must be a number"); return;}
				spec.time_end = end;
				
				if(Number(start) >= Number(end)){
					alert_import("'start' must be before 'end'"); return;
				}
				break;
				
			case "file": 
				spec.time_radio = {value:"File"}; 
				break;
			
			default: 
				alert_import("The value '"+or+"' not recognised"); 
				return;
			}
		}
		break;
		
	case "Diag. Test":
		{
			let Se = get_tag_value("Se");
			if(Se == ""){ cannot_find_tag(); return;}
			
			if(!is_eqn(Se,"Se",{zero_one_range:true})) return;
			
			let Sp = get_tag_value("Sp");
			if(Sp == ""){ cannot_find_tag(); return;}
			
			if(!is_eqn(Sp,"Sp",{zero_one_range:true})) return;
			
			let pos = get_tag_value("pos");
			if(pos == "") pos = "1";
			
			let neg = get_tag_value("neg");
			if(neg == "") pos = "0";
			
			if(pos == neg){
				alert_import("'pos' and 'neg' cannot both have the same value");
			}
			
			let comp = get_tag_value("comp");
			if(comp == ""){ cannot_find_tag(); return;}
		
			let spl = comp.split(",");
			
			let cl = get_cl_from_comp(spl[0],p);
			if(cl == undefined){
				alert_import("Value '"+spl[0]+"' is not a compartment");
				return;
			}
			
			let claa = model.species[p].cla[cl];
			
			let cb = {name:claa.name, value:[]};
			
			for(let c = 0; c < claa.ncomp; c++){
				let name = claa.comp[c].name;
				if(find_in(spl,name) != undefined) cb.value[c] = {comp_name_store:name, check:true};
				else cb.value[c] = {comp_name_store:name, check:false};
			}
	
			let Se_eqn = create_equation(Se,"Se");
			let Sp_eqn = create_equation(Sp,"Sp");
			
			spec = {cl_drop:{te:claa.name}, Se_eqn:Se_eqn, Sp_eqn:Sp_eqn, pos_result:pos, neg_result:neg, check_box:cb};
		}
		break;
		
	case "Population": case "Set Traps":
		{		
			let filt = [];
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				filt[cl] = {comp:[], radio:{value:"All"}};
				for(let c = 0; c < claa.ncomp; c++) filt[cl].comp[c] = false;
			}
			
			let filter = get_tag_value("filter");
			
			if(filter != ""){
				let spl = filter.split(",");
				for(let j = 0; j < spl.length; j++){
					let spl2 = spl[j].split(":");
					if(spl2.length != 2){ 
						alert_import("In 'filter' error understanding '"+filter+"'"); 
						return;
					}
					
					if(spl2[0] == ""){
						alert_import("In 'filter' the value '"+filter+"' does not specify a classification");
						return;	
					}
					
					let cl = find_nocase(sp.cla,"name",spl2[0]);
					if(cl == undefined){ 
						alert_import("In 'filter' the value '"+spl2[0]+"' is not a classification"); 
						return;
					}
					
					let claa = sp.cla[cl];
					
					if(spl2[1].toLowerCase() == "file"){
						filt[cl].radio.value = "File";
					}
					else{
						filt[cl].radio.value = "Comp";
						
						let spl3 = spl2[1].split("|");
						for(let k = 0; k < spl3.length; k++){
							let c = find(claa.comp,"name",spl3[k]);
							
							if(spl3[k] == ""){
								alert_import("In 'filter' the classification '"+spl2[0]+"' is not set");
								return;
							}
							
							if(c == undefined){ 
								alert_import("In 'filter' the value '"+spl3[k]+"' is not a compartment in '"+spl2[0]+"'"); 
								return;
							}
							
							filt[cl].comp[c] = true;
						}
					}
				}
			}

			let cla_fi = [];
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				
				let comp = [];
				
				for(let c = 0; c < claa.ncomp; c++){
					comp[c] = {check: filt[cl].comp[c], comp_name_store:claa.comp[c].name};
				}
						
				cla_fi[cl] = {cl_name_store:claa.name, radio:filt[cl].radio, comp:comp};
			}
	
			spec = {filter:{cla:cla_fi}};
			
			if(type == "Set Traps"){
				let prob = get_tag_value("prob");
				if(prob == ""){ cannot_find_tag(); return;}
			
				if(!is_eqn(prob,"prob",{zero_one_range:true})) return;
				
				spec.trap_prob_eqn = create_equation(prob,"trap_prob",p,undefined);
			}

			if(type == "Population"){
				if(load_obs_model(spec) != "success") return;
			}
		}
		break;
		
	case "Pop. Trans.":
		{
			let from = get_tag_value("from");
			if(from == ""){ cannot_find_tag(); return;}
			let from_spl = from.split(",");
			
			let to = get_tag_value("to");
			if(to == ""){ cannot_find_tag(); return;}
			let to_spl = to.split(",");

			if(from_spl.length != to_spl.length){
				alert_import("The 'from' and 'to' must have the same size list");
				return;
			}
			
			let cl_sel;
			for(let i = 0; i < from_spl.length; i++){
				let cl = get_cl_from_comp(from_spl[i],p);
				if(cl == undefined){
					alert_import("The value '"+from_spl[i]+"' is not a compartment");
					return;
				}
				if(cl_sel == undefined) cl_sel = cl;
				else{
					if(cl != cl_sel){
						alert_import("In 'from' the value '"+from+"' has compartments from different classifications");
						return;
					}
				}
			}
			
			if(cl_sel == undefined){
				alert_import("No transitions have been selected");
				return;
			}
					
			let claa_sel = sp.cla[cl_sel];
			
			let tra_fi = [];
			for(let i = 0; i < claa_sel.ntra; i++){
				let tra = claa_sel.tra[i];
				tra_fi[i] = { tra_name_store:tra.name, check:false};
			}
			
			for(let j = 0; j < from_spl.length; j++){
				let name = from_spl[j]+"→"+to_spl[j];
				let i = find(claa_sel.tra,"name",name);
				if(i == undefined){
					alert_import("The transition '"+name+"' does not exist");
					return;
				}
				
				tra_fi[i].check = true;
			}
		
			let filt = [];
			for(let cl = 0; cl < sp.ncla; cl++){
				if(cl == cl_sel){
					filt[cl] = {};
				}
				else{
					let claa = sp.cla[cl];
					filt[cl] = {comp:[], radio:{value:"All"}};
					for(let c = 0; c < claa.ncomp; c++) filt[cl].comp[c] = false;
				}
			}
		
			let filter = get_tag_value("filter");
			if(filter != ""){
				let spl = filter.split(",");
				for(let j = 0; j < spl.length; j++){
					let spl2 = spl[j].split(":");
					if(spl2.length != 2){ 
						alert_import("In 'filter' error understanding '"+filter+"'"); 
						return;
					}
					
					if(spl2[0] == ""){
						alert_import("In 'filter' the value '"+filter+"' does not specify a classification");
						return;	
					}
					
					let cl = find(sp.cla,"name",spl2[0]);
					if(cl == undefined){
						alert_import("In 'filter' the value '"+spl2[0]+"' is not a classification"); 
						return;
					}
					
					if(cl == cl_sel){ 
						alert_import("In 'filter' the value '"+spl2[0]+"' cannot be the same as the transition classification"); 
						return;
					}
					
					let claa = sp.cla[cl];
					
					if(spl2[1].toLowerCase() == "file"){
						filt[cl].radio.value = "File";
					}					
					else{	
						filt[cl].radio.value = "Comp";
					
						let spl3 = spl2[1].split("|");
						for(let k = 0; k < spl3.length; k++){
							if(spl3[k] == ""){
								alert_import("In 'filter' the classification '"+spl2[0]+"' is not set");
								return;
							}
							
							let c = find(claa.comp,"name",spl3[k]);
							if(c == undefined){
								alert_import("In 'filter' the value '"+spl3[k]+"' is not a compartment in '"+spl2[0]+"'"); 
								return;
							}
							
							filt[cl].comp[c] = true;
						}
					}
				}
			}
			
			let cla_fi = [];
			for(let cl = 0; cl < sp.ncla; cl++){
				if(cl == cl_sel){
					cla_fi[cl] = {};
				}
				else{
					let claa = sp.cla[cl];
					
					let comp = [];
					
					for(let c = 0; c < claa.ncomp; c++){
						comp[c] = {check: filt[cl].comp[c], comp_name_store:claa.comp[c].name};
					}
							
					cla_fi[cl] = {cl_name_store:claa.name, radio:filt[cl].radio, comp:comp};
				}
			}
	
			spec = {cl_drop:{te:claa_sel.name}, filter:{te:claa_sel.name, cla:cla_fi, tra:tra_fi}};
			
			if(load_obs_model(spec) != "success") return;
		}
		break;
		
	case "Genetic":
		{
			let root = get_tag_value("root");
			if(root == ""){ cannot_find_tag(); return;}
			
			spec = {root:root};
		}
		break;
	}
		
	start_data_source(type,spec,info);

	let so = inter.edit_source;

	switch(type){  // Adds in extra information dependent on data source
	case "Init. Pop. Prior": 
		{
			so.variety = so.spec.radio.value;
			if(so.variety == "Flat"){
				data_source("Add");
				return;
			}
		}
		break;
	}
	//if(type == "Init. Pop. Prior": 
	
	so.data_table_use = "not loaded";
	
	let col_name = [];
	for(let c = 0; c < so.load_col.length; c++){
		col_name.push(so.load_col[c].heading);
	}
	
	if(cols != ""){
		let spl = cols.split(",");

		if(spl.length != col_name.length){
			alert_import("'cols' does not have the correct number of entries (expected something in the order '"+stringify(col_name)+"')"); 
			return;
		}
		
		for(let i = 0; i < spl.length; i++){
			for(let j = 0; j < col_name.length; j++){
				if(spl[i] == col_name[j] && i != j){
					alert_import("'cols' does not have the correct order (expected something in the order '"+stringify(col_name)+"')"); 
					return;
				}
			}
		}
		
		col_name = spl;
	}
	
	so.table = get_subtable(tab,col_name);	
	if(so.table.error != ""){ alert_import("File problem – "+so.table.error); return;}

	so.table_loaded = true;
	so.table.edit = true;

	if(data_source_check_error(true) == false){
		data_source("Add");
	}
	else{
		if(inter.help.te != undefined){
			close_data_source();
			alert_import("File problem – "+inter.help.te);
			return;
		}
	}
}


/// Loads up the observation model
function load_obs_model(spec)
{
	let error = get_tag_value("error");
	if(error == ""){ cannot_find_tag(); return;}
			
	let spl = error.split(":");
	if(spl.length != 2){ alert_import("Problem with expression '"+error+"'"); return;}
	
	if(spl[0].toLowerCase() != "normal"){ alert_import("'"+spl[0]+"' should be 'normal'"); return;}
	
	let val = spl[1].trim();
	
	spec.percent = "10";
	spec.sd = "";
			
	if(val.toLowerCase() == "file"){
		spec.obs_error = {value:"file"};	
	}
	else{
		if(is_percent(val) == true){
			spec.obs_error = {value:"percent"};
			spec.percent = val.substr(0,val.length-1);
			if(isNaN(spec.percent)){
				alert_import("The expression '"+spec.percent+"' is not a percentage");
				return;
			}
			
			if(Number(spec.percent) <= 0){
				alert_import("The value '"+spec.percent+"%' must be a positive percentage");
				return;
			}
		}	
		else{
			spec.obs_error = {value:"sd"};
			spec.sd = val;
			
			if(isNaN(spec.sd)){
				alert_import("The expression '"+spec.sd+"' must be a number or a percentage");
				return;
			}
			
			if(Number(spec.sd) <= 0){
				alert_import("The value '"+spec.sd+"' must be positive");
				return;
			}
		}
	}
	
	return "success";
}


/// Reads in parameter samples into sim_results
function sim_param_command()
{
	let file = get_tag_value("file");
	if(file == ""){ cannot_find_tag(); return;}
	
	read_param_samples_file(0,file,sim_result);
}


/// Reads in parameter samples into inf_results
function inf_param_command()
{
	let file = get_tag_value("file");
	if(file == ""){ cannot_find_tag(); return;}
	
	let chain = get_tag_value("chain");
	if(chain == ""){ cannot_find_tag(); return;}
	
	read_param_samples_file(chain,file,inf_result);
}


/// Reads in state samples into sim_results
function sim_state_command()
{
	let file = get_tag_value("file");
	if(file == ""){ cannot_find_tag(); return;}
	
	read_state_samples_file(0,file,sim_result);
}


/// Reads in state samples into inf_results
function inf_state_command()
{
	let file = get_tag_value("file");
	if(file == ""){ cannot_find_tag(); return;}
	
	let chain = get_tag_value("chain");
	if(chain == ""){ cannot_find_tag(); return;}
	
	read_state_samples_file(chain,file,inf_result);
}


/*
/// Load in constant parameters
function sim_const_command()
{
	let file = get_tag_value("file");
	if(file == ""){ cannot_find_tag(); return;}
	
	read_constant(file,sim_result);
}
*/


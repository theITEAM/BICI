/// Bubble which allows user to select a constant parameter
function set_param_dynamic_bubble(cont)
{	
	if(check_param_free("param_dynamic") == true){
		cont.dx = 20;
		bubble_addtitle(cont,"Set dynamic parameter");

		bubble_addparagraph(cont,"Select parameter to be dynamically updated:",0,cont.dx);
		cont.y += 0.2;
	
		bubble_addscrollable(cont,{type:"param dynamic sel", ymax:10, ac:"AddParamDynamic"}); 
		
		add_bubble_end(cont);
	}
	else{
		cont.dx = 15;
		bubble_addtitle(cont,"Set dynamic parameter");

		bubble_addparagraph(cont,"There are currently no free parameters to be dynamically updated.",0,cont.dx); 
	}
}


/// Allows for bubble to allocate param dynamic behaviour
function param_dynamic_bubble(bu,cont)
{
	cont.dx = 11.5;
	bubble_addtitle(cont,"Dynamic update",{te:dynamicupdate_text});
	
	let par = model.param[bu.i];
	
	let dynamic_info = inter.bubble.dynamic_info;
	
	let wmax = 9;
	
	{
		let pos = get_param_dynamic_pos(par);
		bubble_adddropdown(cont,2.5,wmax,dynamic_info.type,pos); cont.y -= 1.1;
		bubble_add_minititle(cont,"Type:");
	}
	
	let ac;
	
	let type = dynamic_info.type.te;
	if(type != select_drop_str){
		let sp_str;
		{  // Checks species is selected
			let sp_pos = [];
			for(let p = 0; p < model.species.length; p++){
				sp_pos.push({te:model.species[p].name, p:p});
			}

			if(sp_pos.length > 1){
				bubble_adddropdown(cont,3.5,wmax-1,dynamic_info.sp_drop,sp_pos); cont.y -= 1.1;
				bubble_add_minititle(cont,"Species:");
			}
			
			let te = dynamic_info.sp_drop.te;
			if(te != select_drop_str) sp_str = te;
		}
		
		let index_sel;
		
		if(sp_str != undefined){ // Select from an index
			switch(type){
			case "bin-thresh-dist": case "bin-thresh-region": case "bin-min-max-dist": case "bin-min-max-region":
				{	
					let pos = get_index_pos(par);
					if(pos.length == 1){
						dynamic_info.index_drop.te = pos[0].te;
					}
					bubble_adddropdown(cont,2.5,wmax,dynamic_info.index_drop,pos); cont.y -= 1.1;
					bubble_add_minititle(cont,"Index:");
					
					let te = dynamic_info.index_drop.te;
					if(te != select_drop_str) index_sel = te;
					cont.y += 0.2;
				}
				break;
			
			default: index_sel = "none"; break;
			}
		}
		
		if(index_sel != undefined && sp_str != undefined){
			ac = "DoneParamDynamic";
			
			let p = find(model.species,"name",sp_str);
		
			switch(type){
			case "bin-thresh-eqn": case "bin-min-max-eqn":
				bubble_input(cont,"Equation:",{type:"dynamic_eqn", eqn:true});
				break;
			
			default:
				{
					let disable = false;
					bubble_addscrollable(cont,{type:"param_dynamic", p:p, index_sel, dynamic_info:dynamic_info, disable:disable, ymax:bubblescroll_dymax}); 
		
					bubble_addcheckbox(cont,0,"Fraction",dynamic_info.fraction);
				}
				break;
			}
			
			switch(type){
			case "bin-thresh-dist": case "bin-thresh-region": case "bin-min-max-dist": case "bin-min-max-region":
				{	
					bubble_addcheckbox(cont,0,"Weight",dynamic_info.weight);
					
					if(dynamic_info.weight.check == true){
						bubble_input(cont,"Weight:",{type:"dynamic_weight", eqn:true});
					}
				}
				break;
			}
			
			switch(type){
			case "bin-thresh-dist": case "bin-min-max-dist":
				bubble_input(cont,"Distance:",{type:"dist"});
				break;
			}
			
			// Select threshold
			switch(type){
			case "bin-thresh": case "bin-thresh-dist": case "bin-thresh-region": case "bin-thresh-eqn":
				bubble_input(cont,"Threshold:",{type:"thresh"});
				break;
				
			case "bin-min-max": case "bin-min-max-dist": case "bin-min-max-region": case "bin-min-max-eqn":
				bubble_input(cont,"Min:",{type:"threshmin"});	
				bubble_input(cont,"Max:",{type:"threshmax"});
				break;
			}
		}
	}

	add_end_button(cont,"Done",ac);	
}


/// Sets p_name and cl_name for dynamic_info weight_eqn
function set_di_weight_pcl(di)
{
	switch(di.type){
	case "bin-thresh-dist": case "bin-min-max-dist": case "bin-thresh-region": case "bin-min-max-region": break;
	default: return;
	}
	
	let p_name = di.sp_drop.te;
	di.weight_eqn.p_name = p_name;
	let p = find(model.species,"name",p_name);
	if(p == undefined){ error("Problem setting weight_pcl"); return;}
	
	let sp = model.species[p];
	let cl = find(sp.cla,"index",di.index_drop.te);
	if(cl == undefined){ error("Problem setting weight_pcl"); return;}
	di.weight_eqn.cl_name = sp.cla[cl].name;
}
	

/// Get the dynamic-sim possibilities for a parameter
function get_param_dynamic_pos(par)
{
	if(par.dep.length > 1) return tensor_dynamic_pos;
	return param_dynamic_pos; 
}


/// Gets the posibilities for the index
function get_index_pos(par)
{
	let pos =[];
	for(let i = 0; i < par.dep.length; i++){
		let te = remove_prime(par.dep[i]);
		if(te != "t") pos.push({te:te});
	}
	
	return pos;
}


/// The scrollable box which allows the used to specify compartments that make up population
function param_dynamic_scrollable(lay)
{
	let cy = 0;
	let dx = lay.dx;
	let w = dx;
	let op = lay.op;
	op.back_col = BUBBLE_COL;

	let mar = 0.2;
	let si = 0.7;
	let fo = get_font(si,"bold");
	let fosi = get_font(si);

	let p = lay.op.p;
	let disable = lay.op.disable;
	
	let dynamic_info = lay.op.dynamic_info;
	let filt = dynamic_info.filter;

	cy += 0.4;

	let sp = model.species[p];
	
	let cl_not = undefined;
	let index_sel = lay.op.index_sel;

	for(let cl2 = 0; cl2 < sp.ncla; cl2++){
		let claa = sp.cla[cl2];
		if(claa.index != index_sel){
			let te = claa.name+" filter:";
			
			lay.add_button({te:te, x:mar, y:cy, dx:w, dy:1, type:"Text", si:si, font:fo, col:BLACK, back_col:BUBBLE_COL});
			cy += 1.1;
		
			lay.add_radio(0,cy,"All","All",filt.cla[cl2].radio,{disable:disable, back_col:BUBBLE_COL});
			lay.add_radio(3.5,cy,"Comp","Comp.",filt.cla[cl2].radio,{disable:disable, back_col:BUBBLE_COL});
				
			cy += 1.6;
			
			switch(filt.cla[cl2].radio.value){
			case "Comp":
				cy = add_checkbox_list(lay,claa.comp,"name",filt.cla[cl2].comp,cy,w);
				break;
			}
		}
	}

	return cy;
}


/// Displays information about dynamic control
function display_dynamic(i,x,y,lay,w)
{
	let source = model;
	let par = source.param[i];
	
	lay.display_param(x-par.label_info.dx-0.7,y-0.1,par.label_info);

	let si = 1.5;
	lay.add_button({te:"=", x:x, y:y, dy:si, type:"Text", font:get_font(si), si:si, col:BLACK});

	let te = dynamic_definition(par.dynamic_info,par);
	
	let ac = "EditParamDynamic"; //if(allow_edit != false) ac = "EditSimValue";
	
	let fo = get_font(1.1,"","times");
	
	lay.add_button({te:te, x:x+1.6, y:y+0., dx:w-x-1.6, dy:1.6, type:"ParamDynamic", source:source, font:fo, ac:ac, i:i, name:par.name, label_info:par.label_info});
}


/// Sets dynamic definition
function dynamic_definition(di,par,check)
{	
	let type = di.type.te;
	
	if(type == select_drop_str) return "Unset dynamic update";
	
	let p; 
	if(model.species == 1) p = 0;
	else p = find(model.species,"name",di.sp_drop.te);

	let vec=[];
	
	let sp = model.species[p];
	let filt = di.filter;

	let filt_te="";  // Generates the population filter text
	
	switch(type){
	case "bin-thresh-eqn": case "bin-min-max-eqn":
		break;
	
	default:
		{	 
			if(model.species.length > 1){
				filt_te += sp_drop.te+":";
			}
			
			let fll = false;
			for(let cl = 0; cl < filt.cla.length; cl++){
				let fcl = filt.cla[cl];
				if(fcl.radio.value == "Comp"){
					let comp = fcl.comp;
					let st = "";
					for(let c = 0; c < comp.length; c++){
						if(comp[c].check == true){
							if(st != "") st += "|";
							st += comp[c].comp_name_store;
						}
					}
					if(fll == true) filt_te += ",";
					filt_te += st;
					fll = true;
				}
			}

			if(di.fraction.check == true) vec.push("frac");
		}
		break;
	}
	
	let warn;
	
	if(par.time_dep != true) warn = "Parameter should be time dependent";
	if(par.spline.on != true) warn = "Parameter should have a spline";
	if(par.spline.spline_radio.value != "Square") warn = "A square spline must be used for dynamic parameter"; 

	switch(type){
	case "bin-thresh-dist": case "bin-min-max-dist": case "bin-thresh-region": case "bin-min-max-region":
		{
			let pos = get_index_pos(par);
			let index = "index:"+di.index_drop.te;
			if(warn == undefined && index == select_drop_str) warn = "Please select index";
			vec.push("index:"+di.index_drop.te);	
			
			if(di.weight.check){
				vec.push("weight:"+di.weight_eqn.te);	
			}
		}
		break;
	}
	
	switch(type){
	case "bin-thresh-dist": case "bin-min-max-dist": 
		if(warn == undefined) warn = check_thresh_value(di.dist,"Distance",false,type);
		vec.push("dist:"+di.dist);	
		break;
	}
	
	switch(type){
	case "bin-thresh": case "bin-thresh-dist": case "bin-thresh-region": case "bin-thresh-eqn":
		if(warn == undefined) warn = check_thresh_value(di.thresh,"Threshold",di.fraction.check,type);
		vec.push("thresh:"+di.thresh);	
		break;
		
	case "bin-min-max": case "bin-min-max-dist": 	case "bin-min-max-region": case "bin-min-max-eqn":
		if(warn == undefined) warn = check_thresh_value(di.threshmin,"Minimum",di.fraction.check,type);
		if(warn == undefined) warn = check_thresh_value(di.threshmax,"Maximum",di.fraction.check,type);
		if(Number(di.threshmin) > Number(di.threshmax)) warn = "Minimum must be less than maximum";
		vec.push("min:"+di.threshmin);	
		vec.push("max:"+di.threshmax);	
		break;
	
	default: error("dynamic-sim type problem"); break;
	}
	
	let te = type+"(";
	
	switch(type){
	case "bin-thresh-eqn": case "bin-min-max-eqn": te += di.eqn.te; break;
	default: te += "{"+filt_te+"}";
	}
	
	for(let i = 0; i < vec.length; i++) te += ", "+vec[i];
	te += ")";
	
	if(check && warn != undefined) return err(warn);
	
	return te;
}


/// Checks that threshold value is suitable set
function check_thresh_value(te,tag,frac,type)
{
	if(type == "bin-thresh-eqn" || type == "bin-min-max-eqn") return;
	
	if(te.trim() == "") return tag+" must be set"; 
	let num = Number(te);
	if(isNaN(num)) return tag+" must be a number";
	else{
		if(num < 0) return tag+" must be non-negative";
		if(frac){
			if(num > 1) return tag+" must be between zero and one";
		}
	}
	return;
}


/// Converts from text to dynamic info
function dynamic_info_from_text(te,par)
{
	let fl = false;
	
	let di = set_init_dynamic_info();
	
	let type;
	
	let vec;
	
	let i = 0; while(i < te.length && te.substr(i,1) != "(") i++;
	if(i == te.length) fl = true;
	else{
		type = te.substr(0,i);
		
		let pos = get_param_dynamic_pos(par);
		let k = find(pos,"te",type);
		if(k == undefined) alert_import("In 'dynamic-sim' the type '"+type+"' is not recognised"); 
		else{
			di.type.te = type;
			
			if(!end_str(te,")")) fl = true;
			else{
				vec = split_with_curly_bracket(te.substr(i+1,te.length-2-i),",");
			}
		}	
	}
	
	if(fl){
		alert_import("Syntax error in '"+te+"' for 'dynamic-sim'"); 
		return;
	}
	
	for(let i = 0; i < vec.length; i++){
		let spl = vec[i].split(":");
		switch(spl.length){
		case 1:
			{
				let val = spl[0].trim();
				if(val == "frac"){
					switch(type){
					case "bin-thresh-eqn": case "bin-min-max-eqn": 
						alert_import("Syntax error in 'dynamic-sim': 'frac' should not be set");
						break;
					default:
						di.fraction.check = true;
						break;
					}
				}
				else{
					switch(type){
					case "bin-thresh-eqn": case "bin-min-max-eqn": 
						di.eqn.te = val;
						break;
					
					default:
						if(begin_str(val,"{")){
							let res = get_pop_filt_from_text(val,di);
							if(res.err) alert_import(res.msg); 
						}
						else{
							alert_import("In '"+te+"' for 'dynamic-sim' the value '"+val+"' is not recognised"); 
						}
						break;
					}
				}
			}
			break;
		default: break;
		}
	}
		
	for(let i = 0; i < vec.length; i++){
		let spl = vec[i].split(":");
		switch(spl.length){
		case 1: break;
		
		case 2:
			{
				let tag = spl[0].trim();
				let val = spl[1].trim();
				switch(tag){	
				case "thresh":
					{
						if(type != "bin-thresh" && type != "bin-thresh-dist" && type != "bin-thresh-region" && type != "bin-thresh-eqn"){
							alert_import("In '"+te+"' for 'dynamic-sim' a value for '"+tag+"' should not be set"); 
						}
						
						{
							let warn = check_thresh_value(val,"Threshold",di.fraction.check,type);
							if(warn != undefined) alert_import("In '"+te+"' for 'dynamic-sim': "+warn); 
						}
						
						di.thresh = val;
					}
					break;
					
				case "min":
					{
						if(type != "bin-min-max" && type != "bin-min-max-dist" && type != "bin-min-max-region" && type != "bin-min-max-eqn"){
							alert_import("In '"+te+"' for 'dynamic-sim' a value for '"+tag+"' should not be set"); 
						}
						
						{
							let warn = check_thresh_value(val,"Minimum",di.fraction.check,type);
							if(warn != undefined) alert_import("In '"+te+"' for 'dynamic-sim': "+warn); 
						}
						
						di.threshmin = val;
					}
					break;
					
				case "max":
					{
						if(type != "bin-min-max" && type != "bin-min-max-dist" && type != "bin-min-max-region"){
							alert_import("In '"+te+"' for 'dynamic-sim' a value for '"+tag+"' should not be set"); 
						}
						
						{
							let warn = check_thresh_value(val,"Maximum",di.fraction.check,type);
							if(warn != undefined) alert_import("In '"+te+"' for 'dynamic-sim': "+warn); 
						}
						
						di.threshmax = val;
					}
					break;
					
				case "dist":
					{
						if(type != "bin-thresh-dist" && type != "bin-min-max-dist"){
							alert_import("In '"+te+"' for 'dynamic-sim' a value for '"+tag+"' should not be set"); 
						}
						
						{
							let warn = check_thresh_value(val,"Distance",false,type);
							if(warn != undefined) alert_import("In '"+te+"' for 'dynamic-sim': "+warn); 
						}
					
						di.dist = val;
					}
					break;
					
				case "index":
					{
						if(type != "bin-thresh-dist" && type != "bin-min-max-dist" && type != "bin-thresh-region" && type != "bin-min-max-region"){
							alert_import("In '"+te+"' for 'dynamic-sim' a value for '"+tag+"' should not be set"); 
						}
						
						di.index_drop.te = val;
					}
					break;
				
				case "weight":
					{
						if(type != "bin-thresh-dist" && type != "bin-min-max-dist" && type != "bin-thresh-region" && type != "bin-min-max-region"){
							alert_import("In '"+te+"' for 'dynamic-sim' a value for '"+tag+"' should not be set"); 
						}
						
						di.weight.check = true;
						di.weight_eqn.te = val;
					}
					break;
					
				default: 
					alert_import("In '"+te+"' for 'dynamic-sim' the tag '"+tag+"' is not recognised"); 
					break; 
				}
			}
			break;
		
		default: 	
			alert_import("Syntax error in '"+te+"' for 'dynamic-sim'"); 
			break;
		}
	}
	
	switch(di.type.te){
	case "bin-min-max": case "bin-min-max-dist": case "bin-min-max-region":
		if(Number(di.threshmin) > Number(di.threshmax)){
			alert_import("For 'dynamic-sim' the minimum must be less than or equal to maximum value");
		}
		break;
	}
	
	set_di_weight_pcl(di);
	
	return di;
}


/// Gets a filter from a string
function get_pop_filt_from_text(te,di)
{
	if(!begin_str(te,"{") || !end_str(te,"}")) return err("Syntax error in expression '"+te+"'");

	let te_st = te;

	te = te.substr(1,te.length-2).trim();
	
	if(te == "") return success();
	
	let p; // Finds the species
	if(model.species.length == 1) p = 0;
	else{
		let spl = te.split(":");
		switch(spl.length){
		case 1: return err("The species must be set in expression '"+te_st+"'");
		case 2:
			p = find(model.species,"name",spl[0]);
			if(p == undefined) return err("The species '"+spl[0]+"' in expression '"+te_st+"' is not recognised");
			te = spl[1];
			break;
		default: return err("Syntax error in expression '"+te_st+"'");
		}
	}
	
	let sp = model.species[p];
	let p_name = sp.name;
	di.sp_drop.te = p_name;
	
	di.filter = load_default_filt(p);	
	 
	let filt = di.filter;
	
	let spl = te.split(",");
	
	for(let i = 0; i < spl.length; i++){
		let cl_sp;
	
		let spl2 = spl[i].split("|");
		for(let k = 0; k < spl2.length; k++){
			let te2 = spl2[k];
			let co = find_comp_from_name(te2,p_name);
		
			if(co == undefined) return err("Syntax error in expression '"+te_st+"', compartment '"+te2+"' not found");
			else{
				if(co.warn != undefined) return err("Syntax error in expression '"+te_st+"': "+co.warn);
				else{
					let cl_name = sp.cla[co.cl].name;
						
					let comp_name = sp.cla[co.cl].comp[co.c].name;
					
					if(cl_sp == undefined) cl_sp = co.cl;
					else{
						if(cl_sp != co.cl){
							return err("Syntax error in expression '"+te_st+"': The compartments '"+spl2[k-1].trim()+"' and '"+spl2[k].trim()+"' are on different classifications");
						}
					}
					
					filt.cla[co.cl].comp[co.c].check = true;
				}
			}
		}
		
		if(filt.cla[cl_sp].radio.value == "Comp"){
			return err("Syntax error in expression '"+te_st+"': Specification on classification '"+filt.cla[cl_sp].cl_name_store+"' are made more than once");
		}
		filt.cla[cl_sp].radio.value = "Comp";
	}
	
	return success();
}
	

/// Sets the initial dynamic parameter information
function set_init_dynamic_info()
{
	let filt;
	let sp_drop = {te:select_drop_str};
	
	if(model.species.length == 1){
		filt = load_default_filt(0);	
		sp_drop.te = model.species[0].name;
	}
			
	return {type:{te:select_drop_str}, filter:filt, thresh:"", threshmin:"", threshmax:"", index_drop:{te:select_drop_str}, dist:"", fraction:{check:false}, sp_drop:sp_drop, region:{loaded:false, source:undefined}, weight:{check:false}, weight_eqn:create_equation("1","dynamic_weight"), eqn:create_equation("","dynamic_eqn")};
}
			

/// This is activated when the 'done' button is pressed on the edit prior bubble
function done_param_dynamic()
{
	if(bubble_check_error() == false){
		copy_back_to_source();
		let bub = inter.bubble;
		
		let di = bub.dynamic_info;
		
		let fl = false;
		switch(di.type.te){
		case "bin-min-max": case "bin-min-max-dist": case "bin-min-max-region":
			if(Number(di.threshmin) > Number(di.threshmax)){
				set_warning("Must be less than or equal to maximum value",["threshmin","threshmax"]);
				fl = true;
			}
			break;
		}
		
		if(!fl){
			model.param[bub.bu.i].dynamic_info = di;
			
			update_param();
			close_bubble();
		}
	}
}

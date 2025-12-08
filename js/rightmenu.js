"use strict";
// Functions to set up the right menu bar

/// Generates content for the right bar
function right_menu_buts(lay)
{	
	let type = lay.op.type;

	if(type != "Diagnostics" && inter.graph.init != true && inter.graph.init != "no data") return;

	let dygap = 0;//0.2;
	
	let rpf; if(type != "GraphView") rpf = get_inf_res().plot_filter;

	let y = 0.1;

	if(inter.graph.init == "no data"){
		y = stand_filt(y,rpf,lay);
		return;
	}

	switch(type){
	case "Populations":
		{
			if(rpf.pos_view.length > 1){
				y = add_filter("View",y,rpf.sel_view,rpf.pos_view,lay);
				
				if(rpf.sel_view.te == "Slice"){
					y = add_slice_time("Time",rpf,y,lay);
				}
				
				y += dygap;		
			}
			
			if(rpf.sel_view.te == "Data"){		
				y = add_filter("Data sources",y,rpf.sel_popdata,rpf.pos_popdata,lay);
				if(!rpf.plot_average) y = chain_sample_filt(y,rpf,lay);
			}
			else{
				y = stand_filt(y,rpf,lay,"pop");
				y = additional_filt(y,rpf,lay);
			}
		}
		break;

	case "Transitions":	
		{
			if(rpf.pos_trans_view.length > 1){
				y = add_filter("View",y,rpf.sel_trans_view,rpf.pos_trans_view,lay);
				y += dygap;	
			}
			
			if(rpf.sel_trans_view.te == "Data"){		
				y = add_filter("Data sources",y,rpf.sel_poptransdata,rpf.pos_poptransdata,lay);
				if(!rpf.plot_average) y = chain_sample_filt(y,rpf,lay);
			}
			else{
				y = stand_filt(y,rpf,lay,"trans");
				y = additional_filt(y,rpf,lay);
			}
			
			if(rpf.pos_timestep.length > 1){
				y = add_filter("Time-step",y,rpf.sel_timestep,rpf.pos_timestep,lay);
			}
		}
		break;
		
	case "Individuals":
		{
			let inf_eff_pl = false;
			
			if(rpf.pos_indview.length > 1 && inter.graph.ind_sel == undefined){
				y = add_filter("View",y,rpf.sel_indview,rpf.pos_indview,lay);
				y += dygap;	
				
				if(rpf.sel_indview.te == "Ind. Eff."){
					if(rpf.pos_indeffview.length > 1){
						y = add_filter("Graph",y,rpf.sel_indeffview,rpf.pos_indeffview,lay);
						y += dygap;
					}
					inf_eff_pl = true;
				}
			}
			
			if(inter.graph.ind_sel){
				y = add_filter("View",y,rpf.sel_ind_sel_view,rpf.pos_ind_sel_view,lay);
				y += dygap;	
			}
	
			if(inf_eff_pl){			
				y = chain_sample_filt(y,rpf,lay);
			}
			else{
				y = stand_filt(y,rpf,lay);		
			}
			
			if(rpf.siminf == "inf" || rpf.siminf == "ppc"){ // Individual group filter
				let p = model.get_p();
				let rpf2 = rpf.species[p];
				
				if(rpf2.pos_indgroup.length > 1){
					y = add_filter("Ind. Group",y,rpf2.sel_indgroup,rpf2.pos_indgroup,lay);
				}
				
				if(rpf.sel_indview.te == "Timeline" || rpf.sel_indview.te == "Table"){
					if(rpf2.pos_sort.length > 1){
						y = add_filter("Sort",y,rpf2.sel_sort,rpf2.pos_sort,lay);
					}
				}
			}
	
			if(false){
				let dx = right_menu_width-1.6;
				lay.add_button({te:"Create SIRE", x:0.1, y:y, dx:dx, dy:dropdown_height, ac:"CreateSIRE", type:"SliceTime"}); y += 2;
			}
			
			if(inf_eff_pl && rpf.sel_indeffview.te == "Scatter"){
				let p = model.get_p();
				let rpf2 = rpf.species[p];
				y = add_colour_filter(y,rpf,rpf2.col_filt,lay);
			}
				
			if(inter.graph.ind_sel){
				y = additional_filt(y,rpf,lay);
			}
			else{
				y = additional_filt(y,rpf,lay);
			}
			
			if(inf_eff_pl && rpf.sel_indeffview.te == "Scatter"){
				let p = model.get_p();
				let rpf2 = rpf.species[p];
				y = add_filter("X-axis",y,rpf2.sel_ie_xaxis,rpf2.pos_ie_axis,lay);
				y = add_filter("Y-axis",y,rpf2.sel_ie_yaxis,rpf2.pos_ie_axis,lay);		
			}
		}
		break;
		
	case "Derived":
		{
			if(rpf.pos_spline.length > 1){
				y = add_filter("Quantity",y,rpf.sel_spline,rpf.pos_spline,lay);
				y += dygap;		
			}
			
			y = chain_sample_filt(y,rpf,lay);
		}
		break;
	
	case "Parameters":
		{
			if(rpf.sel_paramview == undefined) return;
			
			y = add_filter("View",y,rpf.sel_paramview,rpf.pos_paramview,lay);
			y += dygap;				
						
			switch(rpf.sel_paramview.te){
			case "Scatter": case "Correlation": break;
			default:
				if(rpf.pos_paramviewtype){
					if(rpf.pos_paramviewtype.length > 1){
						y = add_filter("Graph",y,rpf.sel_paramviewtype,rpf.pos_paramviewtype,lay);
						y += dygap;
					}
				
					if(rpf.sel_paramviewtype.te == "Distribution"){
						if(rpf.pos_distchain.length > 1){
							y = add_filter("Chain",y,rpf.sel_distchain,rpf.pos_distchain,lay);
						}
					}
					else{
						if(rpf.pos_chain.length > 1){
							y = add_filter("Chain",y,rpf.sel_chain,rpf.pos_chain,lay);
						}
					}
				}
				break;
			}
		}
		break;
	
	case "Generations":
		{
			if(rpf.sel_genview == undefined) return;
			
			y = add_filter("View",y,rpf.sel_genview,rpf.pos_genview,lay);
			y += dygap;		
		}
		break;
	
	case "GraphView":
		{
			let vg = inter.view_graph;
			if(vg.pos_view.length > 1){
				add_filter("View",y,vg.sel_view,vg.pos_view,lay);
			}
		}
		break;
		
	case "Diagnostics":
		{
			y = add_filter("View",y,rpf.sel_diag_view,rpf.pos_diag_view,lay);
			
			let sdv = rpf.sel_diag_view;
			
			switch(sdv.te){
			case "Trans. (exp.)":
				y = stand_filt(y,rpf,lay,"trans");
				y = additional_filt(y,rpf,lay);
				if(rpf.pos_timestep.length > 1){
					y = add_filter("Time-step",y,rpf.sel_timestep,rpf.pos_timestep,lay);
				}
				break;
				
			case "Trans. (dist.)":
				y = stand_filt(y,rpf,lay,"trans");	
				y = additional_filt(y,rpf,lay);
				break;
				
			case "Trans. (bias)": case "Trans. (p-val.)":
				y = add_filter("Graph",y,rpf.sel_bias_view,rpf.pos_bias_view,lay);		
				y = stand_filt(y,rpf,lay,"trans");	
				y = additional_filt(y,rpf,lay);
				break;
			
			case "Proposals":
				if(rpf.pos_diag_chain.length > 1){
					y = add_filter("Chain",y,rpf.sel_diag_chain,rpf.pos_diag_chain,lay);	
				}
				break;
			}
		}
		break;
	}
}


/// A standard set of filters (e.g. as used on population)
function stand_filt(y,rpf,lay,op)
{	
	let p = model.get_p();
	
	let rpf2 = rpf.species[p];
	
	if(rpf2.pos_class.length > 1){
		y = add_filter("Classification",y,rpf2.sel_class,rpf2.pos_class,lay);
	}
	
	if(!(rpf.plot_average && (op == "pop" || op == "trans")) && op != "trans-diag-inf"){
		y = chain_sample_filt(y,rpf,lay);
	}
	
	return y;
}


/// Filters chain and sample
function chain_sample_filt(y,rpf,lay)
{
	if(rpf.pos_chain.length > 1){
		y = add_filter("Chain",y,rpf.sel_chain,rpf.pos_chain,lay);
	}

	if(rpf.sel_chain.te == "All"){
		if(rpf.pos_sample.length > 1){
			y = add_filter("Sample",y,rpf.sel_sample,rpf.pos_sample,lay);
		}
	}
	else{
		let ch = rpf.sel_chain.i;
		if(rpf.pos_sample_chain[ch].length > 1){
			y = add_filter("Sample",y,rpf.sel_sample_chain[ch],rpf.pos_sample_chain[ch],lay);
		}
	}
	
	return y;
}


/// All for filters to be applied to other classification
function additional_filt(y,rpf,lay)
{
	let p = model.get_p();
	let name = subsubtab_name();
	
	let rpf2 = rpf.species[p];
	
	for(let i = 0; i < rpf2.filter.length; i++){
		let te = "";

		if(apply_filter(rpf2,i,name)){
			let rpf3 = rpf2.filter[i];
				
			switch(rpf3.type){
			case "trans_filt":
				{
					let num = 0;
					for(let k = 0; k < rpf3.tra_name.length; k++){
						if(rpf3.tra_filt[k].check == true){
							if(te != "") te += ",";
							te += rpf3.tra_name[k];
							num++;
						}
					}
					if(num == rpf3.tra_name.length) te = "All";
					if(num == 0) te = "None";
				}
				break;
				
			case "pop_filt":
				{
					switch(rpf3.radio.value){
					case "select":
						let num = 0;
						for(let c = 0; c < rpf3.comp_name.length; c++){
							if(rpf3.comp_filt[c].check == true){
								if(te != "") te += ",";
								te += rpf3.comp_name[c];
								num++;
							}
						}
						if(num == rpf3.comp_name.length) te = "All";
						if(num == 0) te = "None";
						break;
						
					case "single":
						te = rpf3.radio_sel.value;
						break;
					
					default: error("Option prob"); break;
					}
					if(rpf3.fraction.check) te = "Frac. "+te;
				}
				break;
				
			default: error("Filt type not recognised"); break;
			}
			
			y = add_population_filter(rpf3.name,y,te,rpf,rpf3,i,lay);
		}
	}
	
	
	let te = "Add filter";

	let si = 0.8;
	let font = get_font(si);
	let w = text_width(te,font);
	lay.add_button({te:te, x:0.1, y:y, dx:w+1.2, dy:1.0, type:"AddFilter", ac:"AddFilter", rpf:rpf});
		
	lay.add_help_button(w+1.5,y+1.4,{title:"Filter", te:filter_text, back_col:WHITE});
		
	y += RIGHT_MENU_DY;
	
	return y;
}


/// Used to add a filter to the model
function add_filter(te,y,source,pos,lay)
{
	let dx = right_menu_width-1.6;
	if(te != undefined){
		lay.add_button({te:te+":", x:0.1, y:y, dx:dx, dy:0.8, type:"InputBoxName", back_col:WHITE});
		y += 0.9;
	}
	
	lay.add_dropdown(0.1,y,dx,10,source,pos);
	
	return y+RIGHT_MENU_DY;
}


/// Used to add a filter to the model
function add_colour_filter(y,rpf,filter,lay)
{
	let dx = right_menu_width-1.6;

	lay.add_button({te:"Colour:", x:0.1, y:y, dx:dx, dy:0.8, type:"InputBoxName"});
	y += 0.9;
	
	let desc = "";
	
	let cs = filter.cla_split;
	for(let k = 0; k < cs.length; k++){
		let fi = cs[k];
		if(fi.checkb.check == true){
			if(desc != "") desc += ",";
			desc += fi.name;
		}
	}
	if(desc == "") desc = "None";

	lay.add_button({te:desc, x:0.1, y:y, dx:dx, dy:dropdown_height, rpf:rpf, filter:filter, ac:"ColourFilt", type:"ColourFilt"});

	return y+RIGHT_MENU_DY;
}


/// Used to add a filter to the model
function add_population_filter(te,y,desc,rpf,filter,num,lay)
{
	let dx = right_menu_width-2.6;
	if(te != undefined){
		lay.add_button({te:te+":", x:0.1, y:y, dx:dx, dy:0.8, type:"InputBoxName"});
		y += 0.9;
	}
	
	lay.add_button({te:desc, x:0.1, y:y, dx:dx, dy:dropdown_height, rpf:rpf, filter:filter, num:num, ac:"PopFilt", type:"PopFilt"});

	let si = 1;
	lay.add_button({x:0.1+dx+0.2, y:y+0.1, dx:si, dy:si, op:{p:filter.p, num:num, rpf:rpf}, ac:"RemoveFilter", type:"Delete"});
	
	return y+RIGHT_MENU_DY;
}


/// Used to add slice time
function add_slice_time(te,rpf,y,lay)
{
	let dx = right_menu_width-1.6;
	lay.add_button({te:te+":", x:0.1, y:y, dx:dx, dy:0.8, type:"InputBoxName"});
	y += 0.9;
	
	lay.add_button({te:rpf.slice_time, x:0.1, y:y, dx:dx, dy:dropdown_height, rpf:rpf, ac:"SliceTime", type:"SliceTime"});
	
	return y+RIGHT_MENU_DY;
}


/// Generates slider for right menu
function right_menu_slider_buts(lay)
{
	lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, type:"TableSlider2", ac:"TableSlider2"});
}


/// Generates content for the right bottom menu
function rightbot_menu_buts(lay)
{
	if(inter.graph.type == undefined) return;

	let gir = get_inf_res();
	if(gir){
		let rpf = gir.plot_filter;
	
		if(rpf){
			switch(subsubtab_name()){
			case "Parameters":
				{
					let selpv = rpf.sel_paramview;
					if(selpv && selpv.te == "Scatter"){
						tot_param_list("",rpf.yaxis_radio,lay);
						return;
					}
				}
				break;
				
			/*
			case "Individuals":
				{
					let sel_ind = rpf.sel_indeffview;
					if(sel_ind && sel_ind.te == "Scatter"){
						let p = model.get_p();
						let seliev = rpf.sel_ie_view2[p]; 
						param_list("",seliev.radio,seliev.list,lay);
						return;
					}
				}
				break;
				*/
			}
		}
	}
	
	let key = inter.graph.op.key;
	
	if(key != undefined && inter.graph.init == true){
		let dy = 1.2;

    let y = lay.dy-key.length*dy-4; 
		if(y < 0) y = 0;
		
		let imax = key.length;
		if(imax > KEY_MAX) imax = KEY_MAX;	
		for(let i = 0; i < imax; i++){
			let ke = key[i];
			
			switch(ke.type){
			case "Line": 
				lay.add_button({te:ke.te, x:0, y:y, dx:lay.dx-0.8, dy:dy, dash:ke.dash, thick:ke.thick, col:ke.col, type:"Key"}); 
				break;
			
			case "Rect": 
				lay.add_button({te:ke.te, x:0, y:y, dx:lay.dx-0.8, dy:dy, col:ke.col, type:"KeyRect"}); 
				break;
			
			case "ErrBar":
				lay.add_button({te:ke.te, x:0, y:y, dx:lay.dx-0.8, dy:dy, col:ke.col, type:"KeyErrBar"}); 
				break;
			
			case "Cross":
				lay.add_button({te:ke.te, x:0, y:y, dx:lay.dx-0.8, dy:dy, col:ke.col, type:"KeyCross"}); 
				break;
				
			case "AddObs": case "RemObs": case "MoveObs": case "TransObs": 
			case "CompObs": case "DiagObs": case "GeneticObs": 
				lay.add_button({te:ke.te, x:0, y:y, dx:lay.dx-0.8, dy:dy, col:ke.col, type:"KeySymb", symb:ke.type});
				break;
			}
			y += dy; 
		}
		
		if(imax != key.length){
			lay.add_button({te:"⋮", x:0, y:y, dx:lay.dx-0.8, dy:dy, type:"KeyNone"});
		}			
	}
	
	let mkey = inter.graph.colour_key;
		
	if(mkey && inter.graph.init == true)
	{
		let dy = 8;
		let y = lay.dy-dy-key_gap; if(y < 0){ dy += y; y = 0;}
		
		lay.add_button({x:1, y:y, dx:lay.dx-1, dy:dy, mkey:mkey, type:"MatrixKey"});	
	}		
}


/// Generates content for the right middle menu
function rightmid_menu_buts(lay)
{
	let rpf = get_inf_res().plot_filter;

	if(rpf){
		switch(subsubtab_name()){
		case "Diagnostics":
			{
				let seldv = rpf.sel_diag_view;
				let p = model.get_p();
				let rpf2 = rpf.species[p];
				let cl = rpf2.sel_class.cl;
				let sel = rpf2.cla[cl].sel;
	
				switch(seldv.te){
				case "Trans. (dist.)": param_list("",sel.radio,sel.list,lay); break;
				}
			}
			break;
		
		case "Parameters":
			{
				let selpv = rpf.sel_paramview;
				
				let title = ""; 
				
				switch(selpv.te){
				case "Correlation":
					param_check_box(rpf.param_check,lay);
					break;
					
				case "Scatter":
					tot_param_list("",rpf.xaxis_radio,lay);
					break;
				
				default:
					switch(rpf.sel_paramviewtype.te){
					case "Graph (split)": param_list(title,selpv.radio_split,selpv.list_split,lay); break;
					case "Histogram": break;
					default: param_list(title,selpv.radio,selpv.list,lay); break;
					}
					break;
				}
			}
			break;
			
		case "Generations":
			{
				let selpv = rpf.sel_genview;
				let title = ""; 
				param_list(title,selpv.radio,selpv.list,lay);
			}
			break;
			
		case "Individuals":
			{
				let p = model.get_p();
				
				let seliev;
				switch(rpf.sel_indeffview.te){
				case "Pred. Acc.": seliev = rpf.sel_ie_data_view[p]; break;
				case "Dist. (log)": seliev = rpf.sel_ie_view[p]; break;
				case "Dist. (norm)": seliev = rpf.sel_ie_view[p]; break;
				case "Scatter": seliev = rpf.sel_ie_view[p]; break;
				default: error("Default option problem"); break;
				}
			
				param_list("",seliev.radio,seliev.list,lay);
			}
			break;
		}
	}
}


/// Generates a list of parameters
function param_list(title,rad,list,lay)
{
	let cy = 0;
	
	if(title != ""){
		let si = 0.8;
		lay.add_button({te:title, x:2, y:cy, dx:5, dy:0.8, si:si, font:get_font(si,"bold"), type:"Text", col:BLACK});
		cy += 1.3;
	}
	
	let too_big = false;
	let imax = list.length; if(imax > PARAM_LIST_MAX){ imax = PARAM_LIST_MAX; too_big = true;}
	for(let i = 0; i < imax; i++){
		let li = list[i];
		lay.add_radio(0,cy,i,list[i].name,rad,{});
		cy += 1.3;	
	}
	if(too_big){ lay.add_text("... too many",2,cy+0.4,BLACK,0.8); cy += 1.3;}	
}
	
	
/// Generates checkboxes to select parameters
function param_check_box(param_check,lay)
{
	let cy = 0;
	
	for(let th = 0; th < param_check.length; th++){
		let pc = param_check[th];
		if(pc){
			lay.add_checkbox(0,cy,pc.name,pc.name,pc.checkb,WHITE);
			cy += 1.3;

			if(pc.checkb.check && pc.list){
				let too_big = false;
				let kmax = pc.list.length;
				if(kmax > PARAM_LIST_MAX){ kmax = PARAM_LIST_MAX; too_big = true;}
				
				if(pc.all_checkb && !too_big){
					lay.add_checkbox(0.5,cy,"All","All",pc.all_checkb,WHITE,{});
					cy += 1.3;	
				}

				for(let k = 0; k < kmax; k++){
					let pcl = pc.list[k];
					lay.add_checkbox(0.5,cy,pcl.name,pcl.name,pcl.checkb,WHITE,{});
					cy += 1.3;	
				}				
				if(too_big){ lay.add_text("... too many",2.5,cy+0.4,BLACK,0.8); cy += 1.3;}	
			
				cy += 0.2;
			}
		}
	}
}


/// Constructs a searchable radio list for parameters
function tot_param_list(title,param_radio,lay)
{
	let cy = 0;
	
	if(title != ""){
		let si = 0.8;
		lay.add_button({te:title, x:2, y:cy, dx:5, dy:0.8, si:si, font:get_font(si,"bold"), type:"Text", col:BLACK});
		cy += 1.3;
	}
	
	for(let i = 0; i < param_radio.param_list.length; i++){
		let prpl = param_radio.param_list[i];
		lay.add_radio(0,cy,i,prpl.name,param_radio.radio,{});
		cy += 1.3;	
		if(param_radio.radio.value == i && prpl.list){
			let too_big = false;
			let jmax = prpl.list.length;
			if(jmax > PARAM_LIST_MAX){ jmax = PARAM_LIST_MAX; too_big = true;}
			for(let j = 0; j < jmax; j++){
				let prpll = prpl.list[j];
				lay.add_radio(0.5,cy,j,prpll.name,prpl.radio,{});
				cy += 1.3;	
			}
			if(too_big){ lay.add_text("... too many",2.5,cy+0.4,BLACK,0.8); cy += 1.3;}	
			
			cy += 0.2;	
		}
	}
}

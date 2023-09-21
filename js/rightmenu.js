"use strict";

/// Generates content for the right bar
function right_menu_buts(lay)
{
	let type = lay.op.type;
	
	let result, rpf;
	if(type != "GraphView"){
		let tree = inter.page_name.split("->");
		
		switch(tree[0]){
		case "Simulation": result = sim_result; break;
		case "Inference": result = inf_result; break;
		default: error("Problem tree"); break;
		}
		
		rpf = result.plot_filter;
	}
	
	let y = 0.1;

	switch(type){
	case "Populations":
		{
			if(rpf.pos_view.length > 1){
				y = add_filter("View",y,rpf.sel_view,rpf.pos_view,lay);
				y += 2;		
			}
			
			if(rpf.sel_view.te == "Data"){		
				y = add_filter("Data sources",y,rpf.sel_popdata,rpf.pos_popdata,lay);
				y = chain_sample_filt(y,rpf,lay,result);
			}
			else{
				y = stand_filt(y,rpf,lay,result);
				y = addutional_filt(y,rpf,lay,result);
			}
		}
		break;

	case "Transitions":	
		{
			if(rpf.pos_trans_view.length > 1){
				y = add_filter("View",y,rpf.sel_trans_view,rpf.pos_trans_view,lay);
				y += 2;	
			}
			
			if(rpf.sel_trans_view.te == "Data"){		
				y = add_filter("Data sources",y,rpf.sel_poptransdata,rpf.pos_poptransdata,lay);
				y = chain_sample_filt(y,rpf,lay,result);
			}
			else{
				y = stand_filt(y,rpf,lay,result);
				y = addutional_filt(y,rpf,lay,result);
			}
			
			if(rpf.pos_timestep.length > 1){
				y = add_filter("Time-step",y,rpf.sel_timestep,rpf.pos_timestep,lay);
			}
		}
		break;
		
	case "Individuals":
		{
			if(rpf.pos_indview.length > 1){
				y = add_filter("View",y,rpf.sel_indview,rpf.pos_indview,lay);
				y += 2;	
			}
			
			y = stand_filt(y,rpf,lay,result);
		}
		break;
		
	case "Splines":
		{
			if(rpf.pos_spline_view.length > 1){
				y = add_filter("View",y,rpf.sel_spline_view,rpf.pos_spline_view,lay);	
				y += 2;		
			}
		
			if(rpf.pos_spline.length > 1){
				y = add_filter("Spline",y,rpf.sel_spline,rpf.pos_spline,lay);
			}
			
			y = chain_sample_filt(y,rpf,lay,result)
		}
		break;
	
	case "Parameters":
		{
			if(rpf.pos_paramview.length > 1){
				y = add_filter("View",y,rpf.sel_paramview,rpf.pos_paramview,lay);
				y += 2;				
			}			
			
			switch(rpf.sel_paramview.te){
			case "Scatter": case "Correlation": break;
			default:
				{
					if(rpf.pos_paramviewtype.length > 1){
						y = add_filter("Graph",y,rpf.sel_paramviewtype,rpf.pos_paramviewtype,lay);
						y += 2;
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
			}
		}
		break;
	
	case "GraphView":
		break;
	}
}


/// A standard set of filters (e.g. as used on population)
function stand_filt(y,rpf,lay,result)
{
	if(rpf.pos_species.length > 1){
		y = add_filter("Species",y,rpf.sel_species,rpf.pos_species,lay);
	}
	
	let p = rpf.sel_species.p;
	
	let rpf2 = rpf.species[p];
	if(rpf2.pos_class.length > 1){
		y = add_filter("Classification",y,rpf2.sel_class,rpf2.pos_class,lay);
	}
	
	y = chain_sample_filt(y,rpf,lay,result);
	
	return y;
}


/// Filters chain and sample
function chain_sample_filt(y,rpf,lay,result)
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
function addutional_filt(y,rpf,lay,result)
{
	let p = rpf.sel_species.p;
	
	let rpf2 = rpf.species[p];
	
	for(let i = 0; i < rpf2.filter.length; i++){
		let rpf3 = rpf2.filter[i];
		y = add_filter(rpf3.name,y,rpf3.sel_comp,rpf3.pos_comp,lay);
	}
	
	lay.add_button({te:"Add filter", x:0.1, y:y, dx:4.5, dy:1.0, type:"AddFilter", ac:"AddFilter", result:result});
	y += 2;
	
	return y;
}

	
function add_filter(te,y,source,pos,lay)
{
	let dx = right_menu_width-1.6;
	if(te != undefined){
		lay.add_button({te:te+":", x:0.1, y:y, dx:dx, dy:0.8, type:"InputBoxName"});
		y += 0.9;
	}
	
	lay.add_dropdown(0.1,y,dx,10,source,pos);
	
	return y+2;
}


/// Generates content for the right bottom menu
function rightbot_menu_buts(lay)
{
	if(inter.graph.type == undefined) return;
	
	if(tab_name() == "Inference"){
		let result = inf_result;
		let rpf = result.plot_filter;

		if(rpf.sel_paramview && subsubtab_name() == "Parameters"){
			let selpv = rpf.sel_paramview;
		
			if(selpv.te == "Scatter"){
				param_list("Y AXIS",selpv.radioy,selpv.list,result,lay);
				return;
			}
		}
	}
	
	let key = inter.graph.op.key;

	if(key != undefined){
		let dy = 1.5;

		let y = lay.dy-key.length*dy-4; 
		if(y < 0) y = 0;
		
		for(let i = 0; i < key.length; i++){
			let ke = key[i];
			
			switch(ke.type){
			case "Line": 
				lay.add_button({te:ke.te, x:0, y:y, dx:lay.dx-0.8, dy:1.3, dash:ke.dash, col:ke.col, type:"Key"}); 
				break;
			
			case "Rect": 
				lay.add_button({te:ke.te, x:0, y:y, dx:lay.dx-0.8, dy:1.3, col:ke.col, type:"KeyRect"}); 
				break;
			
			case "ErrBar":
				lay.add_button({te:ke.te, x:0, y:y, dx:lay.dx-0.8, dy:1.3, col:ke.col, type:"KeyErrBar"}); 
				break;
			
			case "ObsComp": 
				lay.add_button({te:ke.te, x:0, y:y, dx:lay.dx-0.8, dy:1.3, col:ke.col, type:"KeyObsComp"});
				break;
			}
			y += dy; 
		}
	}
	
	if(inter.graph.variety == "Matrix")
	{
		let mkey = inter.graph.matrix_key;
		
		let dy = 12;
		let y = lay.dy-dy-5; 
		
		lay.add_button({x:1, y:y, dx:lay.dx-1, dy:dy, mkey:mkey, type:"MatrixKey"});	
	}
}


/// Generates a list of parameters
function param_list(title,rad,list,result,lay)
{
	let cy = 0;
	
	if(title != ""){
		let si = 0.8;
		lay.add_button({te:title, x:2, y:cy, dx:5, dy:0.8, si:si, font:get_font(si,"bold"), type:"Text", col:BLACK});
		cy += 1.3;
	}
	
	for(let i = 0; i < list.length; i++){
		let th = list[i].th;
		let ind = list[i].index;
		
		let par = result.param[th];
		
		let name;
		if(ind) name = param_name_index(par,ind);
		else name = par.name;
		
		lay.add_radio(0,cy,i,name,rad,{});
		cy += 1.3;	
	}
}
	

/// Generates content for the right middle menu
function rightmid_menu_buts(lay)
{
	if(tab_name() == "Inference"){
		let result = inf_result;
		let rpf = result.plot_filter;
		let selpv = rpf.sel_paramview;
		let title = ""; if(selpv.te == "Scatter") title = "X AXIS";
	
		param_list(title,selpv.radio,selpv.list,result,lay);
	}
}

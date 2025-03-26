"use strict";
// Functions related to data input

class Data
{
	source_width = [];                               // Stores widths of columns for source table (allows for editing)
	
	table_width = [];                                // Stores widths of columns for datatables (in help pages) 
	
	table = [];                                      // Stores data table
	
	filename = "";                                   // Store filename of recently loaded file
	
	constructor()
	{
		this.source_width = [8,23,7,3.5,3,1.5]; 
		this.table_width = [1,17.5,5,5,1.5];
	}
}

	
/// Starts loading a data source
function start_data_source(type,spec,info)
{
	if(info.p == undefined && type != "LoadTensor" && 
		type != "LoadPriorSplit" && type != "LoadDistSplit" && type != "LoadReparam" && type != "LoadAmatrix" &&
		type != "Fixed Effect"){
		info.p = model.get_p();
	}
	
	edit_source = { type:type, table_loaded:false, load_datatable:true, table:{filename:"",heading:[], col_used:[], ele:[], ncol:0, nrow:0, edit:false}, spec:spec, info:info};

	set_loadcol();
}


/// Closes the data source 
function close_data_source()
{
	edit_source = {};
}


/// Generates a page giving initial conditions 
function add_init_cont_buts(lay,siminf)
{
	let cx = corner.x;
	let cy = corner.y;

	let so = edit_source;
	
	let type = so.type;
	
	let p	= model.get_p();
	
	let sp = model.species[p];
	//if(siminf == "gen") sp = model.sim_res.plot_filter.species[p];
	
	let title, te;
	switch(siminf){
	case "sim": title = "Initial conditions"; te = initcond_text; break;
	//case "inf": title = "Data"; te = data_text; break;
	//case "ppc": title = "Population modification"; te = modification_text; break;
	//case "gen": title = "Simulated data"; te = data_text; break;
	}
	
	if(model.species.length > 1) title += " for species '"+sp.name+"'";
		
	cy = lay.add_title(title,cx,cy,{te:te});
	
	let mar = 0.5, mar2 = 1;
	//let cx = mar;
	//let cy = 0.2;

	let dy = 1.3;
	let si = dy*si_comp_text_frac;
	let fo = get_font(si,"","times");
	let si_sup = si*0.6;
	
	pr("sp");
	pr(sp);
	
	let ic = sp.sim_init_cond;
	pr(ic);
	
	cy = lay.add_paragraph(":",lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
	
	lay.add_radio(4.7,cy,"Fixed","Fixed",ic.type,{});
			
	//bubble_addradio(cont,0,"Fixed","Fixed",spec.radio_dist,disable);
	//bubble_addradio(cont,0,"Dist","Distribution",spec.radio_dist,disable);
		return;
		
	switch(ic.type.value){
	case "Fixed":
		switch(ic.focal.value){
		case "Focal":
			{
				if(cx != mar) cy += 1.5;
				cx = mar;
				
				cy = lay.add_paragraph("Set the populations in each of the compartments:",lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
			
				let cl = ic.focal_cl;
				pr("cl");
				pr(cl);
				
				//let te = edit_source.spec.focal.te;
				//let cl = find(sp.cla,"name",te);
				
				let claa = sp.cla[cl];
				
				lay.add_button({te:claa.name, x:mar, y:cy+0.06, dy:dy, type:"SubTitle", col:BLACK}); 
				cy += 2;
				
				let source = edit_source.cla[cl];
				
				for(let c = 0; c < claa.ncomp; c++){
					let co = claa.comp[c];
					let te = co.name+"=";	
					let ac = "EditCompPop";
					let val = source.comp_init_pop[c].pop;
					te += val;
					
					let w = text_width(te,fo);
					let dx = w+2;
					
					if(cx+dx > lay.dx-mar2){ cx = mar; cy += dy+1;}	
					
					lay.add_button({te:te, x:cx, y:cy, dx:dx, dy:dy, p:p, cl:cl, c:c, type:"CompPop", val:val, ac:ac, col:co.col, col_dark:dark_colour(co.col)}); 
					cx += dx+1;
				}
				cx = mar; cy += dy+1;
				
				if(sp.ncla > 1){
					cy += 0.5;
					cy = lay.add_paragraph("Select the percentage of the population in different compartments:",lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
				
					for(let cl = 0; cl < sp.ncla; cl++){
						let claa = sp.cla[cl];
						if(edit_source.spec.focal.te != claa.name){
							lay.add_button({te:claa.name, x:mar, y:cy+0.06, dy:dy, type:"SubTitle", col:BLACK}); 
							cy += 2;
							
							let source = edit_source.cla[cl];
							
							for(let c = 0; c < claa.ncomp; c++){
								let co = claa.comp[c];
								let te = co.name+"=";
								
								let ac = "EditCompPopPercent"; 
								if(c == claa.ncomp-1) ac = "NotEditCompPopPercent";

								let val = source.comp_init_pop[c].pop_per;
								te += precision(val,4)+"%";
								
								let w = text_width(te,fo);
								let dx = w+2;
								
								if(cx+dx > lay.dx-mar2){ cx = mar; cy += dy+1;}	
								
								lay.add_button({te:te, x:cx, y:cy, dx:dx, dy:dy, p:p, cl:cl, c:c, type:"CompPop", val:val, ac:ac, col:co.col, col_dark:dark_colour(co.col)}); 
								cx += dx+1;
							}
							
							cx = mar; cy += dy+1;
						}
					}
				}
			}
			break;
			
		case "All":
			{
				let te = "Specify the population in each of the compartment combinations:";
				if(sp.ncla == 1) te = "Specify the population in each of the compartments:"
				
				cy = lay.add_paragraph(te,lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
				cy += 0.5;
				
				let glob = edit_source.glob_comp;
				
				for(let c = 0; c < glob.length; c++){
					let claa = glob[c].cla;
					
					let col_list=[];
					
					let te = "";
					for(let cl = 0; cl < claa.length; cl++){
						if(cl > 0) te += ",";
						te += claa[cl];
						
						let c = hash_find(sp.cla[cl].hash_comp,claa[cl]);
						col_list.push(sp.cla[cl].comp[c].col);
					}					
					te += "="+glob[c].pop;
					let tsa = text_sup_anno(te,si);
					let w = tsa.w;
					
					let dx = w+2;
								
					if(cx+dx > lay.dx-mar2){ cx = mar; cy += dy+1;}	
								
					lay.add_button({tsa:tsa, x:cx, y:cy, dx:dx, dy:dy, c:c, type:"MultiComp", ac:"EditGlobCompPop", col_list:col_list}); 
					cx += dx+1;
				}
			
				cy += 1.3;
			}
			break;	
		
		default: error("Option not recognised 104"); break;
		}
		break;
		
	case "Dist":
		switch(edit_source.spec.radio2.value){
		case "Focal":
			{
				if(cx != mar) cy += 1.5;
				cx = mar;
				
				cy = lay.add_paragraph("Set the population distribution in each of the compartments:",lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
			
				let te = edit_source.spec.focal.te;
				let cl = find(sp.cla,"name",te);
				
				let claa = sp.cla[cl];
				
				lay.add_button({te:claa.name, x:mar, y:cy+0.06, dy:dy, type:"SubTitle", col:BLACK}); 
				cy += 2;
				
				let source = edit_source.cla[cl];
				
				for(let c = 0; c < claa.ncomp; c++){
					let co = claa.comp[c];
					let val = get_prior_string(source.comp_init_dist[c].dist);
					
					let te = "N", sup = co.name, te2 = " ~ "+val;
					
					let sup_tsa = text_sup_anno(sup,si_sup,LARGE);
						
					let w = text_width(te,fo)+sup_tsa.w+text_width(te2,fo);
									
					let dx = w+2;
					
					if(cx+dx > lay.dx-mar2){ cx = mar; cy += dy+1;}	
					
					lay.add_button({te:te, sup_tsa:sup_tsa, te2:te2, x:cx, y:cy, dx:dx, dy:dy, p:p, cl:cl, c:c, prior:source.comp_init_dist[c].dist, type:"CompPrior", val:val, ac:"EditPopPrior", fo:fo, col:co.col}); 
					cx += dx+1;
				}
				cx = mar; cy += dy+1;
				
				if(sp.ncla > 1){
					cy += 0.5;
					cy = lay.add_paragraph("Select the α value for the Dirichlet distribution:",lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
				
					for(let cl = 0; cl < sp.ncla; cl++){
						let claa = sp.cla[cl];
						if(edit_source.spec.focal.te != claa.name){
							lay.add_button({te:claa.name, x:mar, y:cy+0.06, dy:dy, type:"SubTitle", col:BLACK}); 
							cy += 2;
							
							let source = edit_source.cla[cl];
							
							for(let c = 0; c < claa.ncomp; c++){
								let co = claa.comp[c];
								let val = source.comp_init_dist[c].alpha;
								
								let te = "α", sup = co.name, te2 = " = "+val;
								
								let sup_tsa = text_sup_anno(sup,si_sup,LARGE);
									
								let w = text_width(te,fo)+sup_tsa.w+text_width(te2,fo);
								
								let dx = w+2;
								
								if(cx+dx > lay.dx-mar2){ cx = mar; cy += dy+1;}	
								
								lay.add_button({te:te, sup_tsa:sup_tsa, te2:te2, x:cx, y:cy, dx:dx, dy:dy, source:source.comp_init_dist[c], type:"CompAlpha", val:val, ac:"EditAlpha", fo:fo, col:co.col}); 
								cx += dx+1;
							}
							
							cx = mar; cy += dy+1;
						}
					}
				}
			}
			break;
			
		case "All":
			{
				cy = lay.add_paragraph("Set the total population prior:",lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
				cy += 0.2;
				
				let val = get_prior_string(edit_source.pop_dist);
				let te = "N", sup = "total", te2 = " ~ "+val;
				let sup_tsa = text_sup_anno(sup,si_sup,LARGE);
				let w = text_width(te,fo)+sup_tsa.w+text_width(te2,fo);
				let dx = w+2;
			
				lay.add_button({te:te, sup_tsa:sup_tsa, te2:te2, x:cx, y:cy, dx:dx, dy:dy, prior:edit_source.pop_dist, type:"CompPrior", ac:"EditPopPrior", fo:fo, col:WHITE}); 
			
				cx = mar; cy += dy+1;
			
				cy = lay.add_paragraph("Select the α values for the Dirichlet prior in each of the compartment combinations:",lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
				cy += 0.5;
				
				let glob = edit_source.glob_comp;
				
				for(let c = 0; c < glob.length; c++){
					let claa = glob[c].cla;
					
					let col_list=[];
					
					let te = "α";
					let sup = "";
					for(let cl = 0; cl < claa.length; cl++){
						if(cl > 0) sup += ",";
						sup += claa[cl];
						
						let c = hash_find(sp.cla[cl].hash_comp,claa[cl]);
						col_list.push(sp.cla[cl].comp[c].col);
					}					
					let te2 = " = "+glob[c].alpha;
								
					let sup_tsa = text_sup_anno(sup,si_sup,LARGE);
								
					let w = text_width(te,fo)+sup_tsa.w+text_width(te2,fo);
					let dx = w+2;
								
					if(cx+dx > lay.dx-mar2){ cx = mar; cy += dy+1;}	
								
					lay.add_button({te:te, sup_tsa:sup_tsa, te2:te2, x:cx, y:cy, dx:dx, dy:dy, c:c, source:glob[c], type:"MultiCompAlpha", ac:"EditAlpha", col_list:col_list, fo:fo}); 
					cx += dx+1;
				}
			
				cy += 1.3;
			}
			break;	
		
		default: error("Option not recognised 104"); break;
		}
		break;
	}
	
	lay.add_space(cy,0.3);
	
	
	let info = {siminf:siminf};
		let x = 1.2, y = lay.dy-1.6;
	let gap = 3.4;
		model.add_object_button(lay,"Init. Pop.",x,y,"AddInitialPopulation",{back:WHITE, active:true, info:info, title:"Initial population", te:init_pop_text, siminf:siminf}); 
	
}


/// Generate the page givng the table of data sources
function add_data_buts(lay,siminf)
{
	let cx = corner.x;
	let cy = corner.y;

	let so = edit_source;
	
	let type = so.type;
	
	let p	= model.get_p();
	
	let sp = model.species[p];
	if(siminf == "gen") sp = model.sim_res.plot_filter.species[p];
	
	let title, te;
	switch(siminf){
	case "sim": title = "Initial conditions"; te = initcond_text; break;
	case "infic": title = "Initial conditions"; te = initcond_text; break;
	case "inf": title = "Data"; te = data_text; break;
	case "ppc": title = "Population modification"; te = modification_text; break;
	case "gen": title = "Simulated data"; te = data_text; break;
	}

	if(model.species.length > 1) title += " for species '"+sp.name+"'";
		
	cy = lay.add_title(title,cx,cy,{te:te});
	
	let table;
	if(model.species.length == 0){
		center_message("A species needs to be added to the model.",lay);
		return;
	}
	
	if(siminf == "gen" && model.species[p].generate_pos == false){
		center_message("The model has been changed so data cannot be simulated.",lay);
		return;
	}
	
	let source;
	switch(siminf){
	case "sim": source = sp.sim_source; break;
	case "inf": case "infic": source = sp.inf_source; break;
	case "ppc": source = sp.ppc_source; break;
	case "gen": source = sp.gen_source; break;
	}

	if(source.length == 0){
		table = "There are currently no data sources added.";
		if(siminf == "sim") table = "There is currently no initial condition information.";
	}
	else{
		table = { width:data.source_width, heading:[{name:"Type"},{name:"Details"},{name:"Number"},{name:"Spec."},{name:"Table"},{name:""}], content:[]};
	
		for(let i = 0; i < source.length; i++){
			let so = source[i];
			
			let fn = "N/A"; 
			if(so.table != undefined){
				if(so.table.filename != "") fn = so.table.filename;
			}

			let info = {siminf:siminf, p:p, i:i};
			let view_ob = {te:"View",ac:"ViewSource",info:info};
			let edit_ob = {te:"Edit",ac:"EditSource",info:info};
			
			switch(so.type){
			case "Add Ind.": case "Remove Ind.": case "Move Ind.":
			case "Add Pop.": case "Remove Pop.": 
			case "Ind. Group": //case "Parameter Mult.":
				view_ob.ac = undefined;
				break;
			}
			
			if(so.error == true){ 
				edit_ob = {te:""};
				view_ob = {te:"Invalid", ac:"Invalid", info:info};
			}

			let desc_ob = {te:so.desc};
			table.content.push([{te:so.type},desc_ob,{te:so.num},view_ob,edit_ob,{te:"Delete",ac:"DeleteSource",error:so.error,siminf:siminf,p:p,i:i}]);
		}
	} 
	
	let x = 1.2, y = lay.dy-1.6;
	let gap = 3.4;

	let active = true; if(model.species.length == 0) active = false;

	let info = {siminf:siminf};
	
	let w;

	if(inter.data_type == undefined) inter.data_type = {te:data_types[1]};
	let pos = []; for(let i = 0; i < data_types.length; i++) pos.push({te:data_types[i]});

	let data_ty;

	switch(siminf){
	case "sim": case "ppc": case "infic":
		data_ty = "Init. Cond."; 
		break;

	case "inf": case "gen":
		data_ty = inter.data_type.te;
		lay.add_dropdown(x,y-0.15,6,10,inter.data_type,pos);
		x += 7;
		break;
	}
	
	if(siminf == "gen"){
		lay.add_corner_button([["Copy","Grey","CopyData"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy},{p:p, p_name:sp.name});
	}
	
	if(false){
		if(siminf == "inf"){
			lay.add_corner_button([["SIRE","Grey","CreateSIRE"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy},{p:p, p_name:sp.name});
		}
	}

	switch(data_ty){
	case "Init. Cond.":
		{
			active = true; 
			if(siminf != "ppc"){
				w = model.add_object_button(lay,"Init. Pop.",x,y,"AddInitialPopulation",{back:WHITE, active:active, info:info, title:"Initial population", te:init_pop_text, siminf:siminf}); x += w+gap;
			}
			
			/*
			{
				w = model.add_object_button(lay,"Add Pop.",x,y,"AddPop",{back:WHITE, active:active, info:info, title:"Add Population", te:add_ind_text, siminf:siminf}); x += w+gap;
			
				w = model.add_object_button(lay,"Remove Pop.",x,y,"RemPop",{back:WHITE, active:active, info:info, title:"Remove Population", te:add_ind_text, siminf:siminf}); x += w+gap;
			}
			
			{
				active = true;
				if(sp.type != "Individual") active = false;
				x += 3;

				w = model.add_object_button(lay,"Add Ind.",x,y,"AddIndividuals",{back:WHITE, active:active, info:info, title:"Add Individuals", te:add_ind_text, siminf:siminf}); x += w+gap;
				
				w = model.add_object_button(lay,"Remove Ind.",x,y,"RemIndividuals",{back:WHITE, active:active, info:info, title:"Remove Individuals", te:rem_ind_text, siminf:siminf}); x += w+gap;
				
				w = model.add_object_button(lay,"Move Ind.",x,y,"MoveIndividuals",{back:WHITE, active:active, info:info, title:"Move Individuals", te:move_ind_text, siminf:siminf}); x += w+gap;
			}
			*/
		
			switch(sp.type){
			case "Population":
				w = model.add_object_button(lay,"Add Pop.",x,y,"AddPop",{back:WHITE, active:active, info:info, title:"Add Population", te:add_ind_text, siminf:siminf}); x += w+gap;
			
				w = model.add_object_button(lay,"Remove Pop.",x,y,"RemPop",{back:WHITE, active:active, info:info, title:"Remove Population", te:add_ind_text, siminf:siminf}); x += w+gap;
				break;
			
			case "Individual":
				if(sp.type != "Individual") active = false;

				w = model.add_object_button(lay,"Add Ind.",x,y,"AddIndividuals",{back:WHITE, active:active, info:info, title:"Add Individuals", te:add_ind_text, siminf:siminf}); x += w+gap;
				
				w = model.add_object_button(lay,"Remove Ind.",x,y,"RemIndividuals",{back:WHITE, active:active, info:info, title:"Remove Individuals", te:rem_ind_text, siminf:siminf}); x += w+gap;
				
				w = model.add_object_button(lay,"Move Ind.",x,y,"MoveIndividuals",{back:WHITE, active:active, info:info, title:"Move Individuals", te:move_ind_text, siminf:siminf}); x += w+gap;
				break;
			}
		}
		break;
		
	case "Individual":
		{
			active = true; if(sp.type != "Individual") active = false;
			
			w = model.add_object_button(lay,"Compartment",x,y,"CompData",{ back:WHITE, active:active, info:info, title:"Compartmental data", te:comp_data_text, siminf:siminf}); x += w+gap;
	
			w = model.add_object_button(lay,"Transition",x,y,"TransData",{ back:WHITE, active:active, info:info, title:"Transition data", te:trans_data_text, siminf:siminf}); x += w+gap;
			
			w = model.add_object_button(lay,"Diag. Test",x,y,"DiagTestData",{ back:WHITE, active:active, info:info, title:"Diagnostic test data", te:diag_test_data_text, siminf:siminf}); x += w+gap;
			
			w = model.add_object_button(lay,"Genetic data",x,y,"SeqData",{ back:WHITE, active:active, info:info, title:"Genetic data", te:seq_data_text, siminf:siminf}); x += w+gap;	
		}
		break;

	case "Population":
		{	
			active = true;
			w = model.add_object_button(lay,"Population",x,y,"PopulationData",{ back:WHITE, active:active, info:info, title:"Population data", te:pop_data_text, siminf:siminf}); x += w+gap;
	
			w = model.add_object_button(lay,"Pop. Transition",x,y,"PopTransData",{ back:WHITE, active:active, info:info, title:"Aggregated transition data", te:poptrans_data_text, siminf:siminf}); x += w+gap;
		}
		break;

	case "Additional":
		active = true; if(sp.type != "Individual") active = false;	
		w = model.add_object_button(lay,"Ind. Effect",x,y,"IndEffData",{ back:WHITE, active:active, info:info, title:"Individual effect data", te:ind_eff_data_text, siminf:siminf}); x += w+gap;
	
		w = model.add_object_button(lay,"Ind. Group",x,y,"IndGroupData",{ back:WHITE, active:active, info:info, title:"Individual group data", te:ind_group_data_text, siminf:siminf}); x += w+gap;
		break;
		
	default: error("Option not recognised 20"); break;		
	}

	add_layer("TableContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-3.5,{table:table});
}


/// Based on the worker reply updates which data sources have an error
function update_data_update_source_error(siminf,warn)
{
	for(let p = 0; p < model.species; p++){
		let sp = model.species[p];

		let source;
		switch(siminf){
		case "sim": source = sp.sim_source; break;
		case "inf": source = sp.inf_source; break;
		case "ppc": source = sp.ppc_source; break;
		}
		
		for(let i = 0; i < source.length; i++) source[i].error = false;
	}
	
	for(let j = 0; j < warn.length; j++){
		let wa = warn[j];
		if(wa.mess == "Data error"){
			let sp = model.species[wa.p];
			switch(wa.siminf){
			case "sim": sp.sim_source[wa.i].error = true; break;
			case "inf": sp.inf_source[wa.i].error = true; break;
			case "ppc": sp.ppc_source[wa.i].error = true; break;
			}
		}
	}		
}


/// Detemines if the value is a percentage
function is_percent(val) 
{
	val = String(val);
	if(val.length == 0) return false;
	if(val.substr(val.length-1,1) == "%") return true;
	return false;
}


/// Start adding a datatable
function file_add_datatable()
{
	if(subtab_name() == "Generate Data") sim_data(); 
	else{
		load_data_table_initialise();	
		close_bubble();
	}
}


/// Sets the columns which need to be loaded 
function set_loadcol()
{
	let so = edit_source;

	let x_str = "x position", x_head = "x", x_type = "float";
	let y_str = "y position", y_head = "y", y_type = "float";
	
	if(so.info.p != undefined && so.info.cl != undefined){
		if(model.species[so.info.p].cla[so.info.cl].camera.coord == "latlng"){
			x_str = "latitude"; x_head = "Lat"; x_type = "lat";
			y_str = "longitude"; y_head = "Lng"; y_type = "lng";
		}
	}

	let load_col=[];

	let type = so.type;
	if(type == "Transition"){
		let time_radio = edit_source.spec.time_radio.value;
			
		switch(time_radio){
		case "All": case "Spec": type = "Transition"; break;
		case "File":  type = "Trans. in time range"; break;
		default: error("Can't find option"); break;
		}
	}
	
	let i = find(data_template,"type",type);

	if(i == undefined){ pr(type+"tye"); error("Cannot find data template"); return;}

	let temp = copy(data_template[i]);

	let cols = temp.cols;
	so.title = temp.title;
	so.help = temp.help;
	so.is_data = temp.is_data;
	
	for(let c = 0; c < cols.length; c++){
		let spl = cols[c].split(",");

		if(spl[0] == "snp"){
			if(so.spec.type_radio.value == "matrix") spl[0] = "genobs";
		}
	
		switch(spl[0]){
			case "ID":
				load_col.push({heading:"ID", desc:"individual ID",type:"text"});
				break;

			case "sire":
				load_col.push({heading:"sire", desc:"sire ID",type:"text"});
				break;
				
			case "dam":
				load_col.push({heading:"dam", desc:"dam ID",type:"text"});
				break;

			case "t":
				load_col.push({heading:"t", desc:spl[1], type:"time_float"});
				break;
			
			case "tstart":
				load_col.push({heading:"Start", desc:"measurement start time", type:"time_float"});
				break;

			case "tend":
				load_col.push({heading:"End", desc:"measurement end time", type:"time_float"});
				break;

			case "ttrans":
				load_col.push({heading:"t", desc:spl[1], type:"time_trans_float"});
				break;
				
			case "snp":
				load_col.push({heading:"SNP", desc:spl[1], type:"text"});
				break;
	
			case "genobs":
				load_col.push({heading:"Obs", desc:"ID for genetic observation", type:"text"});
				break;

			case "cl_prob":
				{
					let p = so.info.p;
					if(p == undefined){ error("p should be defined"); return;}
						
					let sp = model.species[p]; 
					
					let cl_drop = so.spec.cl_drop;
					if(cl_drop.te != select_drop_str){
						let cl = find_nocase(sp.cla,"name",cl_drop.te);
						if(cl == undefined){ error("cl should be defined"); return;}
						
						load_col.push({heading:cl_drop.te, desc:"the classification in '"+cl_drop.te+"'",type:"compartment_prob", p:p, cl:cl});
					}
				}
				break;

			case "cl_all":
				{
					let p = so.info.p;
					if(p == undefined){ error("p should be defined"); return;}
							
					let sp = model.species[p];
					
					switch(so.type){
					case "Init. Pop.":
						if(so.spec.radio2 != undefined && so.spec.radio2.value == "Focal"){
							load_col.push({heading:"Compartment", desc:"compartment", type:"comp_all", p:p});
						}
						else{
							for(let cl = 0; cl < sp.ncla; cl++){
								load_col.push({heading:sp.cla[cl].name, desc:"compartment in '"+sp.cla[cl].name+"'",type:"compartment", p:p, cl:cl});
							}
						}
						break;
						
					default:
						for(let cl = 0; cl < sp.ncla; cl++){
							load_col.push({heading:sp.cla[cl].name, desc:"compartment in '"+sp.cla[cl].name+"' individuals are added",type:"compartment", p:p, cl:cl});
						}
						break;
					}
				}
				break;

			case "cl_all_prob":
				{
					let p = so.info.p;
					if(p == undefined){ error("p should be defined"); return;}
						
					let sp = model.species[p];
					for(let cl = 0; cl < sp.ncla; cl++){
						load_col.push({heading:sp.cla[cl].name, desc:"compartment individuals are added",type:"compartment_prob", p:p, cl:cl});
					}
				}
				break;

			case "from":
				{
					let p = so.info.p;
					if(p == undefined){ error("p should be defined"); return;}
					let cl_name = so.spec.cl_drop.te;
					let cl = find_nocase(model.species[p].cla,"name",cl_name);
					
					load_col.push({heading:"From", desc:"the compartment from which individuals come", type:"comp_source_sink", p:p, cl:cl});
				}
				break;

			case "to":
				{
					let p = so.info.p;
					if(p == undefined){ error("p should be defined"); return;}
					let cl_name = so.spec.cl_drop.te;
					let cl = find_nocase(model.species[p].cla,"name",cl_name);
					
					load_col.push({heading:cl_name, desc:"the compartment in to which individuals go", type:"comp_source_sink", p:p, cl:cl});
				}
				break;

			case "add_pop":
				if(so.type != "Add Pop.") error("Should be add pop");
				load_col.push({heading:"Population", desc:"the population added", type:"pos_zero_float"});
				break;
				
			case "rem_pop":
				if(so.type != "Remove Pop.") error("Should be rem pop");
				load_col.push({heading:"Population", desc:"the population removed", type:"pos_zero_float"});
				break;
					
			case "init_pop":
				if(so.type != "Init. Pop.") error("Should be init pop");
				
				switch(so.spec.radio_dist.value){
				case "Fixed":
					if(so.spec.radio2.value == "Focal"){
						load_col.push({heading:"Population", desc:"the population", type:"positive float/percent"});
					}
					else{
						load_col.push({heading:"Population", desc:"the population", type:"pos_zero_float"});
					}
					break;
					
				case "Dist":
					if(so.spec.radio2.value == "Focal"){
						load_col.push({heading:"Distribution", desc:"the distriution/α value", pri_pos:prior_pos_positive, type:"prior/dirichlet"});
					}
					else{
						load_col.push({heading:"Alpha", desc:"the Dirchlet α value", type:"pos_float"});
					}
					break;
				}
				break;
				
			case "start":
				load_col.push({heading:"Start", desc:"the observation start time", type:"float"});
				break;

			case "end":
				load_col.push({heading:"End", desc:"the observation end time", type:"float"});
				break;

			case "result":
				load_col.push({heading:"Result", desc:"the disease diagnostic result", pos_res:so.spec.pos_result, neg_res:so.spec.neg_result, type:"diag test"});
				break;
				
			case "filt_obspop":
				{
					let spec = edit_source.spec;
			
					let filt = spec.filter;
					let p = so.info.p;

					let sp = model.species[p];
					
					for(let cl = 0; cl < sp.ncla; cl++){
						if(filt.cla[cl].radio.value == "File"){
							load_col.push({heading:sp.cla[cl].name, desc:"'"+sp.cla[cl].name+"' which specifies sub-population",type:"compartment_prob", p:p, cl:cl});
						}
					}
				
					load_col.push({heading:"Population", desc:"the observed population estimates", type:"pos_zero_float"});
					
					switch(spec.obs_error_dist.value){
					case "Normal":
						if(spec.obs_error.value == "file"){
							load_col.push({heading:"SD", desc:"the standard deviation in population estimates", type:"pos_zero_float"});
						}
						break;
					case "Negative binomial":
						if(spec.obs_error_p.value == "file"){
							load_col.push({heading:"p", desc:"the probability p in negative binomial error", type:"zero_one_float"});
						}
						break;
					}
				}
				break;
				
			case "filt_obspoptrans":
				{	
					let spec = edit_source.spec;
					let cl_sel = spec.cl_drop.te;
					if(cl_sel != "Select"){
						let filt = spec.filter;
						
						let p = so.info.p;
						let sp = model.species[p];
						let cl = find(sp.cla,"name",cl_sel);

						if(cl == undefined) error("Cannot find cl");
						else{
							for(let cl2 = 0; cl2 < sp.ncla; cl2++){
								if(cl2 != cl){
									if(filt.cla[cl2].radio.value == "File"){
										load_col.push({heading:sp.cla[cl2].name, desc:"'"+sp.cla[cl2].name+"' which specifies sub-population",type:"compartment_prob", p:p, cl:cl2});
									}
								}
							}
							
							load_col.push({heading:"Number", desc:"the observed number of transitions", type:"pos_zero_float"});
							
							switch(spec.obs_error_dist.value){
							case "Normal":
								if(spec.obs_error.value == "file"){
									load_col.push({heading:"SD", desc:"the standard deviation in population-level transition estimates", type:"pos_zero_float"});
								}
								break;
							case "Negative binomial":
								if(spec.obs_error_p.value == "file"){
									load_col.push({heading:"p", desc:"the probability p in negative binomial error", type:"zero_one_float"});
								}
								break;
							}
						}
					}
				}
				break;

			case "comp_name":
				load_col.push({heading:"Name", desc:"the compartment name", type:"comptext"});
				break;
				
			case "start_comp":
				load_col.push({heading:"From", desc:"the compartment from which individuals move", p:so.info.p, cl:so.info.cl, type:"comp_source_sink"});
				break;
				
			case "end_comp":
				load_col.push({heading:"To", desc:"the compartment to which individuals move", p:so.info.p, cl:so.info.cl, type:"comp_source_sink"});
				break;
				
			case "enter_comp":
				load_col.push({heading:"Enter", desc:"the compartment individuals enter", p:so.info.p, cl:so.info.cl, type:"compartment"});
				break;
				
			case "leave_comp":
				load_col.push({heading:"Leave", desc:"the compartment individuals leave", p:so.info.p, cl:so.info.cl, type:"compartment"});
				break;
				
			case "comp_x":
				load_col.push({heading:x_head, desc:"the compartment "+x_str, type:x_type});
				break;
				
			case "comp_y":
				load_col.push({heading:y_head, desc:"the compartment "+y_str, type:y_type});
				break;
			
			case "source_x":
				load_col.push({heading:x_head, desc:"the source "+x_str, type:x_type});
				break;
				
			case "source_y":
				load_col.push({heading:y_head, desc:"the source "+y_str, type:y_type});
				break;
				
			case "sink_x":
				load_col.push({heading:x_head, desc:"the sink "+x_str, type:x_type});
				break;
				
			case "sink_y":
				load_col.push({heading:y_head, desc:"the sink "+y_str, type:y_type});
				break;
				
			case "colour":
				load_col.push({heading:"Color", desc:"the compartment colour", type:"colour"});
				break;
			
			case "trans_value":
				load_col.push({heading:"Value", desc:"the transition distribution", p:so.info.p, cl:so.info.cl, type:"trans_value"});
				break;
				
			case "rate":
				load_col.push({heading:"Rate", desc:"the transition rate", p:so.info.p, cl:so.info.cl, type:"eqn"});
				break;
			
			case "mean":
				load_col.push({heading:"Mean", desc:"the transition mean time", p:so.info.p, cl:so.info.cl, type:"eqn"});
				break;
				
			case "period":
				load_col.push({heading:"Time", desc:"the transition time", p:so.info.p, cl:so.info.cl, type:"eqn"});
				break;
		
			case "cv":
				load_col.push({heading:"CV", desc:"the transition coefficient of variation", p:so.info.p, cl:so.info.cl, type:"eqn"});
				break;
				
			case "sd":
				load_col.push({heading:"SD", desc:"the transition standard deviation", p:so.info.p, cl:so.info.cl, type:"eqn"});
				break;
			
			case "shape":
				load_col.push({heading:"Shape", desc:"the transition shape parameter", p:so.info.p, cl:so.info.cl, type:"eqn"});
				break;
				
			case "shape (integer)":
				load_col.push({heading:"Shape (integer)", desc:"the transition shape parameter", p:so.info.p, cl:so.info.cl, type:"positive nonzero int"});
				break;
			
			case "scale":
				load_col.push({heading:"Scale", desc:"the transition scale parameter", p:so.info.p, cl:so.info.cl, type:"eqn"});
				break;
				
			case "boundary":
				load_col.push({heading:"Boundary", desc:"", type:"boundary"});
				break;
				
			case "dep":
				{
					let par = model.param[inter.edit_param.i];	
					let dep = par.dep;
			
					let desc = convert_tensor_text("the tensor");
					desc += " index giving '";
					for(let k = 0; k < dep.length; k++){
						let de = remove_prime(dep[k]);
						if(de == "t"){
							load_col.push({heading:dep[k], desc:"knot times", knot:par.spline.knot, type:"knot"});
						}
						else{
							let res = find_cla_from_index(de);
						
							if(res != undefined){
								let claa = model.species[res.p].cla[res.cl];
								load_col.push({heading:dep[k], desc:desc+claa.name+"'", p:res.p, cl:res.cl, type:"compartment"});
							}
						}
					}
				}
				break;
				
			case "dep_index":
				{
					let name = so.spec.full_name;
					let pp = get_param_prop(name);
					let dep = pp.dep;
					let desc = " '";
					for(let k = 0; k < dep.length; k++){
						let de = remove_prime(dep[k]);
						if(de != "t"){
							let res = find_cla_from_index(de);
						
							if(res != undefined){
								let claa = model.species[res.p].cla[res.cl];
								load_col.push({heading:dep[k], desc:desc+claa.name+"'", p:res.p, cl:res.cl, type:"compartment"});
							}
						}
					}
				}
				break;
				
			case "value":
				load_col.push({heading:"Value", desc:"the element value", type:"float"});
				break;
			
			case "prior":
				{
					let par = model.param[inter.edit_param.i];
					load_col.push({heading:"Prior", desc:"the element prior", pri_pos:par.pri_pos, type:"prior"});
				}
				break;
				
			case "dist":
				{
					let par = model.param[inter.edit_param.i];
					load_col.push({heading:"Distribution", desc:"the element distribution", pri_pos:par.pri_pos, type:"dist"});
				}
				break;
				
			case "eqn":
				load_col.push({heading:"Equation", desc:"the equation", type:"eqn"});
				break;
				
			case "A":
				break;
		
			default: error(spl[0]+" in data column not recognised"); break;
		}
	}

	so.load_col = load_col;
}

	
/// Allows data sources to be added and editted
function data_source(type,edit_source,info)
{
	if(info == undefined) info = edit_source.info;
	let p = info.p; 
	if(p == undefined){ error("Data source does not have p"); return}

	let sp = model.species[p];
	if(info.siminf == "gen") sp = model.sim_res.plot_filter.species[p];
	
	let so;
	switch(type){
	case "Add":
		so = edit_source;
		add_source_description(so);
		break;
		
	case "Edit":
		break;
	
	case "ReplaceView":
		switch(info.siminf){
		case "sim": so = sp.sim_source[info.i]; break;
		case "inf": so = sp.inf_source[info.i]; break;
		case "ppc": so = sp.ppc_source[info.i]; break;
		case "gen": so = sp.gen_source[info.i]; break;
		}
		so.spec = edit_source.spec;
		add_source_description(so);
		break;
		
	case "Replace": 
		so = edit_source;
		so.edit_spec = false;
		if(so.table == undefined){ // When changing spec copies back table
			switch(info.siminf){
			case "sim": so.table = sp.sim_source[info.i].table; break;
			case "inf": so.table = sp.inf_source[info.i].table; break;
			case "ppc": so.table = sp.ppc_source[info.i].table; break;
			case "gen": so.table = sp.gen_source[info.i].table; break;
			}
		}
		add_source_description(so);
		break;
		
	case "Delete":
		break;
		
	default: error("Error option"); break;
	}

	switch(info.siminf){
	case "sim": 
		switch(type){
		case "Add": sp.sim_source.push(so); break;
		case "Edit": so = sp.sim_source[info.i]; break;
		case "Replace": sp.sim_source[info.i] = so; break;
		case "ReplaceView": break;
		case "Delete": 
			sp.sim_source.splice(info.i,1); 
			break; 
		default: error("Error option"); break;
		}
		break;
		
	case "inf": 
		switch(type){
		case "Add": sp.inf_source.push(so); break;
		case "Edit": so = sp.inf_source[info.i]; break;
		case "Replace": sp.inf_source[info.i] = so; break;
		case "ReplaceView": break;
		case "Delete":
			sp.inf_source.splice(info.i,1);
			break; 
		default: error("Error option"); break;
		}
		break;
		
	case "ppc": 
		switch(type){
		case "Add": sp.ppc_source.push(so); break;
		case "Edit": so = sp.ppc_source[info.i]; break;
		case "Replace": sp.ppc_source[info.i] = so; break;
		case "ReplaceView": break;
		case "Delete":
			sp.ppc_source.splice(info.i,1);
			break; 
		default: error("Error option"); break;
		}
		break;
		
	case "gen": 
		switch(type){
		case "Add": sp.gen_source.push(so); break;
		case "Edit": so = sp.gen_source[info.i]; break;
		case "Replace": sp.gen_source[info.i] = so; break;
		case "ReplaceView": break;
		case "Delete": 
			sp.gen_source.splice(info.i,1);
			break; 
		default: error("Error option"); break;
		}
		break;
		
	default: error("Option not recognised 21b"+info.siminf); break;
	}
	
	if(info.imp == true) return;
	
	if(type == "Edit"){
		so.info = info;
		data_source_check_error("worker",so);
		
		post({ type:"Edit Data Source", edit_source:so});
		return;
	}
	
	switch(info.siminf){
	case "sim": post({p:p, siminf:info.siminf, sim_source:strip_heavy(sp.sim_source)}); break;
	case "inf": post({p:p, siminf:info.siminf, inf_source:strip_heavy(sp.inf_source)}); break;
	case "ppc": post({p:p, siminf:info.siminf, ppc_source:strip_heavy(sp.ppc_source)}); break;
	case "gen": post({p:p, siminf:info.siminf, gen_source:strip_heavy(sp.gen_source)}); break;
	}
}


/// Gets reply from data_source
function data_source_reply(ans)
{
	close_data_source();
	
	switch(ans.siminf){
	case "sim": model.species[ans.p].sim_source = ans.sim_source; break;
	case "inf": model.species[ans.p].inf_source = ans.inf_source; break;
	case "ppc": model.species[ans.p].ppc_source = ans.ppc_source; break;
	case "gen": model.sim_res.plot_filter.species[ans.p].gen_source = ans.gen_source; break;
	}
}


/// Adds a description of the data source
function add_source_description(so)
{
	if(so.error == true) return;
	
	let desc="";
	let num="-";
	
	let sp = model.species[so.info.p];
	let spec = so.spec;
	
	switch(so.type){
	case "Init. Pop.":
		switch(spec.radio_dist.value){
		case "Fixed":		
			switch(spec.radio2.value){
			case "Focal":
				{	
					let focal = spec.focal.te;
					let cl = find(sp.cla,"name",focal);
					if(cl == undefined){
						desc += "'"+focal+"' no longer exists";
					}
					else{
						let cip = so.cla[cl].comp_init_pop;
						let sum = 0;
						for(let c = 0; c < cip.length; c++) sum += Number(cip[c].pop);
						num = sum;
						if(sp.ncla > 1) desc += "Focal '"+sp.cla[cl].name+"'";
					}
				}
				break;
				
			case "All": 
				{	
					let sum = 0;
					for(let c = 0; c < so.glob_comp.length; c++){
						sum += Number(so.glob_comp[c].pop);
					}
					num = sum;
				}
				break;
			}
			break;
			
		case "Dist":
			switch(spec.radio2.value){
			case "Focal":
				{	
					let foc_cl = spec.focal.te;
					let cl = find(sp.cla,"name",foc_cl);
					if(cl == undefined) error("Cannot find classification");
					else{
						if(sp.cla[cl]){
							for(let c = 0; c < sp.cla[cl].ncomp; c++){
								if(c != 0) desc += ", ";
								if(sp.cla[cl].comp && so.cla[cl] && so.cla[cl].comp_init_dist[c]){
									desc += "{"+sp.cla[cl].comp[c].name+"} ~ "+get_prior_string(so.cla[cl].comp_init_dist[c].dist);
								}
							}
						}
					}
				}
				break;
				
			case "All":
				{	
					desc += "N ~ "+get_prior_string(so.pop_dist);
				}
			}
			break;	
		}
		break;

	case "Add Ind.": case "Remove Ind.": 
		num = so.table.nrow;
		break;
	
	case "Add Pop.": case "Remove Pop.": 
		{
			let sum = 0;
			let tab = so.table;
			for(let r = 0; r < tab.nrow; r++){
				sum += Number(tab.ele[r][tab.ncol-1]);
			}
			num = sum;
		}
		break;
		
	case "Move Ind.":
		num = so.table.nrow;
		break;

	case "Compartment": 
		{
			let cl_se = find(sp.cla,"name",spec.cl_drop.te);
			desc = spec.cl_drop.te;
			num = so.table.nrow;
		}
		break;
		
	case "Diag. Test":
		{
			let cl_se = find(sp.cla,"name",spec.cl_drop.te);
			
			let sen = "";
			let val = spec.check_box.value;
			for(let c = 0; c < val.length; c++){
				if(val[c].check == true){
					if(sen != "") sen += ",";
					sen += val[c].comp_name_store;
				}
			}
					
			desc = spec.cl_drop.te+" "+sen+" Se="+spec.Se_eqn.te+" Sp="+spec.Sp_eqn.te;
			
			num = so.table.nrow;
		}
		break;
		
	case "Ind. Eff.":
		{
			desc = spec.drop.te;
			num = so.table.nrow;
		}
		break;
		
	case "Ind. Group":
		{
			desc = spec.gname;
			num = so.table.nrow;
		}
		break;
	
	case "Population":	
		desc = comp_filter_desc(spec.filter.cla,sp);
		num = so.table.nrow;
		break;
		
	case "Transition": case "Pop. Trans.":
		{
			let filt = spec.filter;
			
			let trz = filt.tra;
			
			let cl_se = find(sp.cla,"name",spec.cl_drop.te);
			if(cl_se == undefined) desc = "'"+spec.cl_drop.te+"' no longer exists";
			else{
				if(cl_se >= sp.cla.length){ so.error = true; break;}
				
				let tra = sp.cla[cl_se].tra;
				if(tra.length != trz.length){ so.error = true; break;}
				
				let st = "";
				for(let j = 0; j < trz.length; j++){
					switch(filt.trans_obs_model.value){
					case "on":
						{
							let pro = trz[j].prob_eqn.te.trim();
							
							if(pro != "0"){
								if(st != "") st += "|";
								st += tra[j].name+":"+pro;
							}
						}
						break;
					
					case "off":
						if(trz[j].check == true){
							if(st != "") st += "|";
							st += tra[j].name;
						}
						break;
					}
				}
				if(st != "") desc += st+"  ";
				
				let te = comp_filter_desc(spec.filter.cla,sp);

				if(te != "") desc += te+" ";
				
				if(spec.time_radio && spec.time_radio.value == "Spec"){
					desc += "(t="+spec.time_start+" - "+spec.time_end+")";
				}
	
				num = so.table.nrow;
			}
		}
		break;

	case "Genetic":
		switch(spec.type_radio.value){
		case "matrix": desc = "Matrix"; break;
		case "snp": desc = "# SNP "+(so.table.ncol-2); break;
		}
		num = so.table.nrow;
		break;
		
	case "APed":
		num = so.table.nrow;
		break;
	
	
	default: error(so.type+" DESC NOT ADDED"); break;
	}
	
	if(so.error == true) so.desc = "Invalid data source";
	else{
		if(desc == "") desc = "-";
		so.desc = desc;
	}
	so.num = num;
}


/// Gives a description of a compartment filter
function comp_filter_desc(clz,sp)
{
	let te = "", st = "";
	for(let cl = 0; cl < clz.length; cl++){
		let clzz = clz[cl];
		if(clzz != undefined && clzz.radio != undefined){
			switch(clzz.radio.value){
			case "Comp":
				{
					let te2 = "";
					for(let j = 0; j < clzz.comp.length; j++){
						if(clzz.comp[j].check == true){
							if(te2 != "") te2 += "|";
							te2 += sp.cla[cl].comp[j].name;
						}
					}
					te += te2+" ";
				}
				break;
			
			case "ObsMod":
				{
					let te2 = "";
					for(let j = 0; j < clzz.comp.length; j++){
						let pro = clzz.comp[j].prob_eqn.te.trim();
						if(pro != "0"){
							if(te2 != "") te2 += "|";
							te2 += sp.cla[cl].comp[j].name+":"+pro;
						}
					}
					te += te2+" ";
				}
				break;

			case "File":
				if(st != "") st += ",";
				st += sp.cla[cl].name;
				break;

			case "All":
				break;

			default: error("Option not recognised 23"+clzz.radio.value); break;
			}
		}
	}
	if(te != "") te += "  ";
	if(st != "") te += st+"  ";
	
	return te;
}


/// Gets a list of all possible classifications
function get_cl_pos()
{
	let sp = model.get_sp();

	let cl_pos = [];
	for(let cl = 0; cl < sp.ncla; cl++){
		if(find(cl_pos,"te",sp.cla[cl].name) == undefined){
			cl_pos.push({te:sp.cla[cl].name});
		}
	}

	return cl_pos;
}

				
/// Comverts from a name to a compartment name
function make_comp_name(name)
{
	name = name.replace(/ /g,"-");
	name = name.replace(/_/g,"-");
	name = name.replace(/\(/g,"|");
	name = name.replace(/\)/g,"|");
	name = name.replace(/</g,"|");
	name = name.replace(/>/g,"|");
	name = name.replace(/{/g,"|");
	name = name.replace(/}/g,"|");
	
	return name;
}


/// Removes repeated entries (used in testing)
function remove_repeated(tab)
{
	let c = 0;
	for(let j = 0; j < tab.nrow; j++){
		pr(j+"j");
		let i = j+1; while(i < tab.nrow && tab.ele[i][c] != tab.ele[j][c]) i++;
		if(i < tab.nrow){
			tab.ele.splice(j,1); j--;
			tab.nrow--;
		}
	}
	
	let s = "";
	for(let i = 0; i < tab.ncol; i++){
		if(i != 0) s += ",";
		s += tab.heading[i];
	}
	s += "\n";
	
	for(let j = 0; j < tab.nrow; j++){
		for(let i = 0; i < tab.ncol; i++){
			if(i != 0) s += ",";
			s += '"'+tab.ele[j][i]+'"';
		}
		s += "\n";
	}	

	write_file_async(s,"D:/AI temp/pop.csv");
}


/// Used to generate a bit transition table (for testing);
function create_trans_data(tab)
{
	let te = "From,To,rate\n";
	for(let i = 0; i < 40000; i++){
		let ci = Math.floor(Math.random()*tab.nrow);
		let cf = Math.floor(Math.random()*tab.nrow);
		te += '"'+tab.ele[ci][0]+'"';
		te += ",";
		te += '"'+tab.ele[cf][0]+'"';
		te += ",0.01\n";
	}
	
	write_file_async(te,"D:/AI temp/trans_big.csv");
}


/// Sets a default filter
function load_default_filt(p)
{
	let sp = model.species[p];
	let filt={};

	filt.cla=[];	
	for(let cl = 0; cl < sp.ncla; cl++){
		let claa = sp.cla[cl];
		filt.cla[cl] = { cl_name_store:claa.name, radio:{value:"All"}, comp:[]};
		for(let c = 0; c < claa.ncomp; c++){
			filt.cla[cl].comp[c] = {check:false, prob_eqn:create_equation("0","comp_prob",p,cl), comp_name_store:claa.comp[c].name};
		}
	}
	
	return filt;
}


/// Checks that paramater is correctly specified
function edit_param_check_error(error_mess_only)
{
	return false;
}


/// Checks is an equation is valid
function check_eqn_valid(te)
{
	let eqn = create_equation(te,"test");		
	if(eqn.warn.length > 0) return err(eqn.warn[0].te);
	
	return success();
}

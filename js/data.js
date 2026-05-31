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
		this.source_width = [10,8,15,5,3.5,3,1.5]; 
		this.table_width = [1,17.5,5,5,1.5];
	}
}

	
/// Starts loading a data source
function start_data_source(type,spec,info)
{
	if(info.p == undefined && type != "LoadTensor" &&
		type != "LoadPriorSplit" && type != "LoadDistSplit" && type != "LoadReparam" && type != "LoadDefine" && type != "LoadAmatrix" &&
		type != "Fixed Effect"){
		info.p = model.get_p();
	}
	
	// Gets the default name for the data source
	let name;
	switch(type){
	case "Init. Pop.": name = "Initial population"; break;
	case "Move Ind.": name = "Individual move data"; break;
	case "Compartment": name = "Compartment data"; break;
	case "Transition": name = "Transition data"; break;
	case "Diag. Test": name = "Diagnostic test data"; break;
	case "Test-and-cull": name = "Test and cull"; break;
	case "Population": name = "Pop. data"; break;
	case "Pop. Trans.": name = "Pop. trans. data"; break;
	case "Genetic": name = "Genetic data"; break;
	case "Ind. Eff.": name = "True individual effect"; break;
	case "Set IE": name = "Set individual effect"; break;
	case "Ind. Group": name = "Group"; break;
	case "Add Pop.": name = "Added populations"; break;
	case "Remove Pop.": name = "Removed populations"; break;
	case "Add Ind.": name = "Added individuals"; break;
	case "Remove Ind.": name = "Removed individuals"; break;
	case "KnotTimes": name = "Knot times"; break;
	case "APed": name = "Pedigree data"; break;
	case "Fixed Effect": name = "Fixed effect data"; break;
	default: name = "Data"; break;
	}
	
	{ // Adds number if multiple sources with the same name
		let siminf = info.siminf;
		if(siminf != undefined){
			let p = info.p;
			
			let source = get_source(siminf,p);
			
			let num = 1;
			do{
				let name2 = name; if(num > 1) name2 += num;
				if(find(source,"name",name2) == undefined){
					name = name2; break;
				}
				num++;
			}while(true);
		}
	}
	
	if(type == 'Test-and-cull' && spec.time_gap == undefined){
		spec.time_gap = "";
	}
	
	edit_source = { name:name, type:type, table_loaded:false, load_datatable:true, table:{filename:"",heading:[], col_used:[], ele:[], ncol:0, nrow:0, edit:false}, spec:spec, info:info};
	
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
	case "sim": title = "Simulation setup"; te = initcond_text; break;
	}
	
	if(model.species.length > 1) title += " for species '"+sp.name+"'";
		
	cy = lay.add_title(title,cx,cy,{te:te});
	
	let mar = 0.5, mar2 = 1;
	
	let dy = 1.3;
	let si = dy*si_comp_text_frac;
	let fo = get_font(si,"","times");
	let si_sup = si*0.6;
	
	let ic = sp.sim_init_cond;
	
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
	let nsp = get_so_nsp(siminf);

	let title, te;
	switch(siminf){
	case "sim": title = "Simulation setup"; te = initcond_text; break;
	case "inf": title = "Data"; te = data_text; break;
	case "ppc": title = "Population modification"; te = modification_text; break;
	case "gen": title = "Simulated data"; te = data_sim_text; break;
	}
	
	cy = lay.add_title(title,cx,cy,{te:te});
	
	let table;
	if(nsp == 0){
		center_message("A species needs to be added to the model.",lay);
		return;
	}
	
	let sp = get_so_sp(siminf,p);
	let source = get_source(siminf,p);
	
	if(nsp > 1) title += " for species '"+sp.name+"'";
	
	if(siminf == "gen" && (p >= model.species.length || model.species[p].generate_pos == false)){
		center_message("The model has changed, so data cannot be simulated.",lay);
		return;
	}
		
	if(source.length == 0){
		table = "No data sources added.";
		if(siminf == "sim") table = "No setup information.";
	}
	else{
		table = { width:data.source_width, heading:[{name:"Name"},{name:"Type"},{name:"Details"},{name:"Number"},{name:"Spec."},{name:"Table"},{name:""}], content:[]};
	
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
			case "Ind. Group":
				view_ob.ac = undefined;
				break;
			
			case "Test-and-cull":
				edit_ob.ac = undefined;
				break;
			}
			
			if(so.error == true){ 
				edit_ob = {te:""};
				view_ob = {te:"Invalid", ac:"Invalid", info:info};
			}

			let desc_ob = {te:so.desc};
			table.content.push([{te:"Source Name",name:so.name,ac:"EditSourceName",so:so},{te:so.type},desc_ob,{te:so.num},view_ob,edit_ob,{te:"Delete",ac:"DeleteSource",error:so.error,siminf:siminf,p:p,i:i}]);
		}
	} 
	
	let x = 1.2, y = lay.dy-1.6;
	let gap = 3.4;

	let active = true; if(nsp == 0) active = false;

	let info = {siminf:siminf};
	
	let w;
	
	let data_ty;

	switch(siminf){
	case "sim": case "ppc": 
		{
			if(inter.data_type_sim == undefined) inter.data_type_sim = {te:data_types_sim[0]}; // Default data type
			let pos = []; for(let i = 0; i < data_types_sim.length; i++) pos.push({te:data_types_sim[i]});		

			data_ty = inter.data_type_sim.te;
			lay.add_dropdown(x,y-0.15,6.4,10,inter.data_type_sim,pos);
			x += 7.4;
		}
		break;

	case "inf": case "gen":
		{
			if(inter.data_type == undefined) inter.data_type = {te:data_types[0]}; // Default data type
			let pos = []; for(let i = 0; i < data_types.length; i++) pos.push({te:data_types[i]});
		
			data_ty = inter.data_type.te;
			lay.add_dropdown(x,y-0.15,6.4,10,inter.data_type,pos);
			x += 7.4;
		}
		break;
	}
	
	if(siminf == "gen"){
		let rpf = model.sim_res.plot_filter;
		if(rpf.pos_sim.length > 1){
			lay.add_dropdown(x,y-0.15,6,10,rpf.sel_sim,rpf.pos_sim);
			x += 7;
		}
	}
	
	if(siminf == "gen"){
		lay.add_corner_button([["Copy","Grey","CopyData"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy},{p:p, p_name:sp.name});
	}
	
	if(false){
		if(siminf == "inf"){
			lay.add_corner_button([["SIRE","Grey","CreateSIRE"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy},{p:p, p_name:sp.name});
		}
	}

	//data_ty = "Individual"; prr("Enfore data type");

	switch(data_ty){
	case "Dynamic Int.":
		{
			let te = test_and_cull_text, ti = "Implement test-and-cull";
						
			w = model.add_object_button(lay,"Test-and-cull",x,y,"TestAndCull",{ back:WHITE, active:active, info:info, title:ti, te:te, siminf:siminf}); 
			x += w+gap;
		}
		break;
	
	case "Init. Cond.":
		{
			active = true; 
			if(siminf != "ppc"){
				let te = init_pop_text, ti = "Initial population";
				if(siminf == "gen"){ te = sim_init_pop_text; ti = "Generate initial population data";}
				w = model.add_object_button(lay,"Init. Pop.",x,y,"AddInitialPopulation",{back:WHITE, active:active, info:info, title:ti, te:te, siminf:siminf}); 
				x += w+gap;
			}
			
			switch(sp.type){
			case "Population": case "Deterministic":
				{
					let te = add_pop_text, ti = "Add Population";
					if(siminf == "gen"){ te = sim_add_pop_text; ti = "Generate added population data";}
					w = model.add_object_button(lay,"Add Pop.",x,y,"AddPop",{back:WHITE, active:active, info:info, title:ti, te:te, siminf:siminf}); 
					x += w+gap;
				}
				
				{
					let te = remove_pop_text, ti = "Remove Population";
					if(siminf == "gen"){ te = sim_remove_pop_text; ti = "Generate remnoved population data";}
					w = model.add_object_button(lay,"Remove Pop.",x,y,"RemPop",{back:WHITE, active:active, info:info, title:ti, te:te, siminf:siminf}); 
					x += w+gap;
				}
				break;
			
			case "Individual":
				if(sp.type != "Individual") active = false;

				{
					let te = add_ind_text, ti = "Add Individuals";
					if(siminf == "gen"){ te = sim_add_ind_text; ti = "Generate added individual data";}		
					w = model.add_object_button(lay,"Add Ind.",x,y,"AddIndividuals",{back:WHITE, active:active, info:info, title:ti, te:te, siminf:siminf}); x += w+gap;
				}
				
				{
					let te = rem_ind_text, ti = "Remove Individuals";
					if(siminf == "gen"){ te = sim_rem_ind_text; ti = "Generate removed individual data";}			
					w = model.add_object_button(lay,"Remove Ind.",x,y,"RemIndividuals",{back:WHITE, active:active, info:info, title:ti, te:te, siminf:siminf}); x += w+gap;
				}
				
				{
					let te = move_ind_text, ti = "Move Individuals";
					if(siminf == "gen"){ te = sim_move_ind_text; ti = "Generate move individual data";}			
					w = model.add_object_button(lay,"Move Ind.",x,y,"MoveIndividuals",{back:WHITE, active:active, info:info, title:ti, te:te, siminf:siminf}); x += w+gap;
				}
				
				if(siminf == "sim"){
					let active2 = active;
					if(sp.ind_eff_group.length == 0) active2 = false;
					let te = fix_IE_text, ti = "Set individual effect";
					w = model.add_object_button(lay,"Set IE",x,y,"SetIE",{back:WHITE, active:active2, info:info, title:ti, te:te, siminf:siminf}); x += w+gap;
				}
				break;
			}
		}
		break;
		
	case "Individual":
		{
			active = true; if(sp.type != "Individual") active = false;
			
			{
				let te = comp_data_text, ti = "Compartmental data";
				if(siminf == "gen"){ te = sim_comp_data_text; ti = "Generate individual compartmental data";}			
				w = model.add_object_button(lay,"Compartment",x,y,"CompData",{ back:WHITE, active:active, info:info, title:ti, te:te, siminf:siminf});
				x += w+gap;
			}
			
			{
				let te = trans_data_text, ti = "Transition data";
				if(siminf == "gen"){ te = sim_trans_data_text; ti = "Generate individual transition data";}			
				w = model.add_object_button(lay,"Transition",x,y,"TransData",{ back:WHITE, active:active, info:info, title:ti, te:te, siminf:siminf});
				x += w+gap;
			}
			
			{
				let te = diag_test_data_text, ti = "Diagnostic test data";
				if(siminf == "gen"){ te = sim_diag_test_data_text; ti = "Generate individual diagnostic test data";}			
				w = model.add_object_button(lay,"Diag. Test",x,y,"DiagTestData",{ back:WHITE, active:active, info:info, title:ti, te:te, siminf:siminf}); 
				x += w+gap;
			}
			
			if(siminf == "gen"){
				if(sp.pos_test_and_cull != undefined){
					let te = test_and_cull_data_text, ti = "Test-and cull diagnostic test data";
					w = model.add_object_button(lay,"Test-and-cull",x,y,"TestAndCullData",{ back:WHITE, active:true, info:info, title:ti, te:te, siminf:siminf}); 
					x += w+gap;
				}
			}
			
			{
				let te = seq_data_text, ti = "Genetic data";
				if(siminf == "gen"){ te = sim_seq_data_text; ti = "Generate pathogen genetic data";}			
				w = model.add_object_button(lay,"Genetic data",x,y,"SeqData",{ back:WHITE, active:active, info:info, title:ti, te:te, siminf:siminf}); 
				x += w+gap;	
			}
		}
		break;

	case "Population":
		{	
			active = true;
			{
				let te = pop_data_text, ti = "Population data";
				if(siminf == "gen"){ te = sim_pop_data_text; ti = "Generate population data";}			
				w = model.add_object_button(lay,"Population",x,y,"PopulationData",{ back:WHITE, active:active, info:info, title:ti, te:te, siminf:siminf}); 
				x += w+gap;
			}
			
			{
				let te = poptrans_data_text, ti = "Aggregated transition data";
				if(siminf == "gen"){ te = sim_poptrans_data_text; ti = "Generate aggregated transition data";}		
				w = model.add_object_button(lay,"Pop. Transition",x,y,"PopTransData",{ back:WHITE, active:active, info:info, title:ti, te:te, siminf:siminf});
				x += w+gap;
			}
		}
		break;

	case "Additional":
		if(siminf == "gen"){	
			let param = model.sim_res.param; 
			let active = false;
			for(let i = 0; i < model.param.length; i++){
				let par = model.param[i];
				if(param_pos(par,"priorconst")) active = true;
			}
	
			let te = const_prior_text, ti = "Create constant prior";
			w = model.add_object_button(lay,"Const. Prior",x,y,"ConstPriorData",{ back:WHITE, active:active, info:info, title:ti, te:te, siminf:siminf}); 
			x += w+gap;
		}
		
		active = true; if(sp.type != "Individual") active = false;	
		
		{
			let te = ind_eff_data_text, ti = "Individual effect data";
			if(siminf == "gen"){ te = sim_ind_eff_data_text; ti = "Generate individual effect data";}		
			w = model.add_object_button(lay,"Ind. Effect",x,y,"IndEffData",{ back:WHITE, active:active, info:info, title:ti, te:te, siminf:siminf}); 
			x += w+gap;
		}
		
		{
			let te = ind_group_data_text, ti = "Individual group data";
			if(siminf == "gen"){ te = sim_ind_group_data_text; ti = "Generate individual group data";}		
			w = model.add_object_button(lay,"Ind. Group",x,y,"IndGroupData",{ back:WHITE, active:active, info:info, title:ti, te:te, siminf:siminf}); 
			x += w+gap;
		}
		break;
		
	default: error("Option not recognised 20"); break;		
	}

	add_layer("TableContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-3.5,{table:table});
}


/// Gets the source of data
function get_source(siminf,p)
{
	switch(siminf){
	case "sim": return model.species[p].sim_source; 
	case "inf": return model.species[p].inf_source; 
	case "gen": return model.sim_res.plot_filter.species[p].gen_source;
	case "ppc": return model.inf_res.plot_filter.species[p].ppc_source;
	}
	error("Cannot get source");
}


/// Gets the source of data
function set_source(siminf,p,source)
{
	switch(siminf){
	case "sim": model.species[p].sim_source = source; break;
	case "inf": model.species[p].inf_source = source; break; 
	case "gen": model.sim_res.plot_filter.species[p].gen_source = source; break;
	case "ppc": model.inf_res.plot_filter.species[p].ppc_source = source; break;
	default: error("Cannot set source:"+siminf); break;
	}
}


/// Gets number of species for source
function get_so_nsp(siminf)
{
	switch(siminf){
	case "sim": return model.species.length; 
	case "inf": return model.species.length; 
	case "gen":
		{
			let pf = model.sim_res.plot_filter;
			if(pf == undefined) return;
			return pf.species.length;
		}
	case "ppc":
		{
			let pf = model.inf_res.plot_filter;
			if(pf == undefined) return;
			return pf.species.length;
		}
	}
	error("Cannot get source species");
}


/// Gets model species for source
function get_so_sp(siminf,p)
{
	if(siminf == undefined) return model.species[p]; 
		
	switch(siminf){
	case "sim": return model.species[p]; 
	case "inf": return model.species[p]; 
	case "gen": 
		if(is_interface || sim_result.species == undefined) return model.sim_res.plot_filter.species[p];
		return sim_result.species[p];
	case "ppc": 
		if(is_interface || inf_result.species == undefined) return model.inf_res.plot_filter.species[p];
		return inf_result.species[p];
	}
	error("Cannot get source species:"+siminf);
}


/// Based on the worker reply updates which data sources have an error
function update_data_update_source_error(siminf,warn)
{
	for(let p = 0; p < model.species; p++){
		let source = get_source(siminf,p);
		for(let i = 0; i < source.length; i++) source[i].error = false;
	}
	
	for(let j = 0; j < warn.length; j++){
		let wa = warn[j];
		if(wa.mess == "Data error"){
			let source = get_source(wa.siminf,wa.p);
			source[wa.i].error = true; 
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
		let sp = get_so_sp(so.info.siminf,so.info.p);
		if(sp.cla[so.info.cl].camera.coord == "latlng"){
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

	if(i == undefined){ error("Cannot find data template"); return;}

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
	
		let cty = spl[0];
		switch(cty){
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
						
					let sp = get_so_sp(so.info.siminf,p);
					
					let cl_drop = so.spec.cl_drop;
					if(cl_drop.te != select_drop_str){
						let cl = find_nocase(sp.cla,"name",cl_drop.te);
						if(cl == undefined){ error("cl should be defined"); return;}
						
						load_col.push({heading:cl_drop.te, desc:"the classification in '"+cl_drop.te+"'",type:"compartment_prob_na", p:p, cl:cl});
					}
				}
				break;

			case "cl_all":
				{
					let p = so.info.p;
					if(p == undefined){ error("p should be defined"); return;}
				
					let sp = get_so_sp(so.info.siminf,p);
					
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
						
					let sp = get_so_sp(so.info.siminf,p);
				
					for(let cl = 0; cl < sp.ncla; cl++){
						load_col.push({heading:sp.cla[cl].name, desc:"compartment individuals are added",type:"compartment_prob", p:p, cl:cl});
					}
				}
				break;

			case "from":
				{
					let p = so.info.p;
					if(p == undefined){ error("p should be defined"); return;}
					let sp = get_so_sp(so.info.siminf,p);
					
					let cl_name = so.spec.cl_drop.te;
					let cl = find_nocase(sp.cla,"name",cl_name);
					
					load_col.push({heading:"From", desc:"the compartment from which individuals come", type:"comp_source_sink", p:p, cl:cl});
				}
				break;

			case "to":
				{
					let p = so.info.p;
					if(p == undefined){ error("p should be defined"); return;}
					let sp = get_so_sp(so.info.siminf,p);
					
					let cl_name = so.spec.cl_drop.te;
					let cl = find_nocase(sp.cla,"name",cl_name);
					
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

					let sp = get_so_sp(so.info.siminf,p);
					
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
						let sp = get_so_sp(so.info.siminf,p);
			
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
					let ep = inter.edit_param;
					let par = ep.par_st; if(par == undefined) par = model.param[ep.i];
					let ndep = par.ndep_cont;
					let dep = par.dep;
			
					let desc = convert_tensor_text("the tensor");
					desc += " index giving '";
					for(let k = 0; k < ndep; k++){
						let de = remove_prime(dep[k]);
						if(de == "t"){
							load_col.push({heading:dep[k], desc:"knot times", knot:par.spline.knot, type:"knot"});
						}
						else{
							let res = find_cla_from_index(de);
						
							if(res != undefined){
								let sp = get_so_sp(so.info.siminf,res.p);
								
								let claa = sp.cla[res.cl];
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
								let sp = get_so_sp(so.info.siminf,res.p);
								let claa = sp.cla[res.cl];
								load_col.push({heading:dep[k], desc:desc+claa.name+"'", p:res.p, cl:res.cl, type:"compartment"});
							}
						}
					}
				}
				break;
				
			case "value":
				load_col.push({heading:"Value", desc:"the element value", type:"float"});
				break;
				
			case "pos_value":
				load_col.push({heading:"Value", desc:"the element value", type:"pos_float"});
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
				
			case "eqn": case "define_ele": case "reparam_ele": 
				load_col.push({heading:"Equation", desc:"the equation", type:cty});
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

	let source = get_source(info.siminf,p);
	
	let so;
	switch(type){
	case "Add":
		so = edit_source;
		add_source_description(so);
		break;
		
	case "Edit":
		break;
	
	case "ReplaceView":
		so = source[info.i];
		so.spec = edit_source.spec;
		add_source_description(so);
		break;
		
	case "Replace": 
		so = edit_source;
		so.edit_spec = false;
		if(so.table == undefined){ // When changing spec copies back table
			so.table = source[info.i].table;
		}
		add_source_description(so);
		break;
		
	case "Delete":
		break;
		
	default: error("Error option"); break;
	}

	switch(type){
	case "Add": source.push(so); break;
	case "Edit": so = source[info.i]; break;
	case "Replace": source[info.i] = so; break;
	case "ReplaceView": break;
	case "Delete": source.splice(info.i,1); break; 
	default: error("Error option"); break;
	}
	
	/*
	switch(info.siminf){
	case "sim": 
		
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
	*/

	if(info.imp == true) return;
	
	if(type == "Edit"){
		so.info = info;
		data_source_check_error("worker",so);
		
		post({ type:"Edit Data Source", edit_source:so});
		return;
	}
	
	post({p:p, siminf:info.siminf, source:strip_heavy(source)});
}


/// Gets reply from data_source
function data_source_reply(ans)
{
	close_data_source();
	
	let source = set_source(ans.siminf,ans.p,ans.source);
	/*
	switch(ans.siminf){
	case "sim": source = ans.sim_source; break;
	case "inf": source = ans.inf_source; break;
	case "ppc": source = ans.ppc_source; break;
	case "gen": source = ans.gen_source; break;
	}
	*/
}


/// Adds a description of the data source
function add_source_description(so)
{
	if(so.error == true) return;
	
	let desc="";
	let num="-";
	
	let sp = get_so_sp(so.info.siminf,so.info.p);
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
		
	case "Diag. Test": case "Test-and-cull":
		{
			let cl_se = find(sp.cla,"name",spec.cl_drop.te);
			
			let sen = "";
			let val = spec.check_box.value;
			for(let c = 0; c < val.length; c++){
				if(val[c].check == true){
					if(sen != "") sen += ",";
					sen += val[c].comp_name_store;
					
					sen += "[Se:"+val[c].Se_eqn.te+"]";
				}
			}
					
			desc = spec.cl_drop.te+" "+sen+" Sp:"+spec.Sp_eqn.te;
			
			if(so.type == "Test-and-cull");// desc += " (gap="+spec.time_gap+")";
			else num = so.table.nrow;
		}
		break;
		
	case "Ind. Eff.": case "Set IE":
		{
			desc = "["+spec.drop.te+"]";
			num = so.table.nrow;
		}
		break;
		
	case "Ind. Group":
		{
			desc = "Ind.: ";
			let tab = so.table;
			for(let i = 0; i < tab.nrow; i++){
				if(i != 0) desc += ",";
				desc += tab.ele[i][0];
				if(desc.length > 100) break;
			}
			num = tab.nrow;
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
								if(st != "") st += " | ";
								st += tra[j].name+" : "+pro;
							}
						}
						break;
					
					case "off":
						if(trz[j].check == true){
							if(st != "") st += " | ";
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
		let clz_ = clz[cl];
		if(clz_ != undefined && clz_.radio != undefined){
			switch(clz_.radio.value){
			case "Comp":
				{
					let te2 = "";
					for(let j = 0; j < clz_.comp.length; j++){
						if(clz_.comp[j].check == true){
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
					for(let j = 0; j < clz_.comp.length; j++){
						let pro = clz_.comp[j].prob_eqn.te.trim();
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

			default: error("Option not recognised 23"+clz_.radio.value); break;
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


/// Checks to see if prior is correctly specified
function check_prior_split(type,par,prior_split)
{
	let dist = false; if(par.variety == "dist") dist = true;
		
	let co_list = generate_co_list(par.list);
	for(let i = 0; i < co_list.length; i++){
		let comb = co_list[i];
		let el = get_element(prior_split,comb.index);
		
		if(!el) return "This contains unset elements";
		if(!el.type) return "This contains unset elements";
		if(el.type.te == select_str) return "This contains unset elements";
		
		let warn = check_prior(el,dist);
		
		if(warn != undefined){
			let type = "prior"; if(dist) type = "distribution";
			return element_error(par,comb.index,type,warn);
		}
	}
}


/// Generates an error for an element
function element_error(par,index,type,warn)
{
	let te = "Problem with ";
	
	te += " for parameter "+par.full_name+": ";
	te += warn;
	te += " (element ";

	for(let i = 0; i < index.length; i++){
		if(i != 0) te += ",";
		te += par.list[i][index[i]];
	}
	te += ")";

	return te;
}

	
/// Checks a define or reparam element
function define_reparam_element_check(el,is_define,time_dep)
{
	let ty;
	if(time_dep){
		if(is_define) ty = "define_ele"; else ty = "reparam_ele";
	}
	else{
		if(is_define) ty = "define_ele_notime"; else ty = "reparam_ele_notime";
	}
	
	let eqn = create_equation(String(el),ty);
	if(eqn.warn.length > 0) return eqn.warn[0].te;
	return "";
}


/// Checks parameter values
function unescape_param_value(type,par,value)
{
	let co_list = generate_co_list(par.list,par.ndep_cont);
	for(let i = 0; i < co_list.length; i++){
		let comb = co_list[i];
		let el = get_element(value,comb.index);
		if(isNaN(el)){
			el = remove_escape_char(el);
			el = char_replace(el);
			set_element(value,comb.index,el);
		}
	}
}

		
/// Checks parameter values
function check_param_value(type,par,value)
{
	let co_list = generate_co_list(par.list,par.ndep_cont);
	for(let i = 0; i < co_list.length; i++){
		let comb = co_list[i];
		let el = get_element(value,comb.index);
		
		if(el == undefined){
			return "This contains unset elements";
		}
		
		switch(par.variety){
		case "define":
			{
				let warn = define_reparam_element_check(el,true,par.time_dep);
				if(warn != "") return element_error(par,comb.index,"definition",warn);
			}
			break;
			
		case "reparam":
			{
				let warn = define_reparam_element_check(el,false,par.time_dep);
				if(warn != "") return element_error(par,comb.index,"definition",warn);
			}
			break;
		}
	}
	
	switch(type){
	case "Set Param":
		if(par.factor){
			let weight;
			if(par.factor_weight_on.check == true) weight = par.factor_weight;
			let num = 0;
			let numi;
			let sum = 0;
			let wsum = 0, wi;
			for(let i = 0; i < co_list.length; i++){
				let comb = co_list[i];
				let el = get_element(value,comb.index);
				let val = Number(el);
				if(val < 0) return "The value '"+el+"' must be positive for a factor.";
					
				let w = 1;
				if(weight) w = Number(get_element(weight,comb.index));
					
				if(el == "*"){ num++; numi = i; wi = w;}
				else{
					sum += w*Number(el);
				}
				wsum += w;
			}
			
			if(num > 1){
				return "Parameter factors shouldn't contain more than one element with value '*'";
			}
			
			if(num == 1){
				let val_new = (wsum-sum)/wi;
				if(val_new < 0){
					return "The calculated value for '*' is '"+round_small(val_new)+"' which must be positive";
				}
				
				let comb = co_list[numi];
				set_element(value,comb.index,round_small(val_new));
			}
			
			if(num == 0){
				let mean = round_small(sum/wsum);
				if(dif(mean,1)){
					if(!weight){ 	
						return "The mean of elements is "+mean+" and it should be 1";
					}
					else{
						return "The weighted mean of elements is "+mean+" and it should be 1";
					}
				}
			}
		}
		
		if(par.time_dep && par.variety == "reparam"){
			if(par.spline.spline_radio.value != "Square"){	
				for(let i = 0; i < co_list.length; i++){
					let comb = co_list[i];
					let el = String(get_element(value,comb.index));
					if(contains_population(el)){
						return "The value '"+el+"' can only contain a population if the spline type is set to square.";
					}
				}
			}
		}
		break;
		
	case "Set Weight":
		{
			for(let i = 0; i < co_list.length; i++){
				let comb = co_list[i];
				let el = get_element(value,comb.index);
				if(isNaN(el)) return "The value '"+el+"' must be a number for a weight.";
				let val = Number(el);
				if(val <= 0) return "The value '"+el+"' must be positive for a weight.";
			}
		}
		break;
	}
}

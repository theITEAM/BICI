"use strict";

class Data
{
	source_width = [];             // Stores widths of columns for source table (allows for editing)
	
	table_width = [];              // Stores widths of columns for datatables (in help pages) 
	
	table = [];
	
	filename = "";                 // Store filename of recently loaded file
	
	constructor()
	{
		this.source_width = [8,12,18,3.5,3,1.5]; 
		this.table_width = [1,17.5,5,5,1.5];
		//this.popdata_width = [8,10,18.5,3,1.5];
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
	
	inter.edit_source = { type:type, table_loaded:false, load_datatable:true, table:{filename:"",heading:[], col_used:[], ele:[], ncol:0, nrow:0, edit:false}, spec:spec, info:info};

	set_loadcol();
}


/// Closes the data source 
function close_data_source()
{
	inter.edit_source = {};
}


/// Generate the page givng the table of data sources
function add_data_buts(lay,siminf)
{
	let cx = corner.x;
	let cy = corner.y;

	let so = inter.edit_source;
	
	let type = so.type;
	
	let p	= model.get_p();
	
	let sp = model.species[p];
	if(siminf == "gen") sp = sim_result.species[p];
	
	let title, te;
	switch(siminf){
	case "sim": title = "Population"; te = population_text; break;
	case "inf": title = "Data"; te = data_text; break;
	case "gen": title = "Simulated data"; te = data_text; break;
	}

	if(model.species.length > 1) title += " for species '"+sp.name+"'";
		
	cy = lay.add_title(title,cx,cy,{te:te});
	
	let table;
	if(model.species.length == 0){
		center_message("A species needs to be added to the model.",lay);
		return;
	}
	
	let source;
	switch(siminf){
	case "sim": source = sp.sim_source; break;
	case "inf": source = sp.inf_source; break;
	case "gen": source = sp.gen_source; break;
	}

	if(source.length == 0){
		table = "There are currently no data sources added.";
	}
	else{
		table = { width:data.source_width, heading:[{name:"Type"},{name:"File"},{name:"Details"},{name:"Spec."},{name:"Table"},{name:""}], content:[]};
	
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
			case "Init. Pop. Prior": 
				view_ob.ac = undefined;
				if(so.variety == "Flat") edit_ob.ac = undefined;
				break;
			
			case "Add Ind.": case "Remove Ind.": case "Move Ind.": 
			case "Ind. Captured": case "Ind. Trapped":
				view_ob.ac = undefined;
				break;
			}
			if(so.error == true){ 
				edit_ob = {te:""};
				view_ob = {te:"Invalid", ac:"Invalid", info:info};
			}

			table.content.push([{te:so.type},{te:fn},{te:so.desc},view_ob,edit_ob,{te:"Delete",ac:"DeleteSource",source:source,i:i}]);
		}
	} 
	
	let x = 1.2, y = lay.dy-1.6;
	let gap = 3.4;

	let active = true; if(model.species.length == 0) active = false;

	let info = {siminf:siminf};
	
	let w;

	if(inter.data_type == undefined) inter.data_type = {te:data_types[0]};
	let pos = []; for(let i = 0; i < data_types.length; i++) pos.push({te:data_types[i]});

	let data_ty;

	switch(siminf){
	case "sim": 
		data_ty = "Init. Pop."; 
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

	switch(data_ty){
	case "Init. Pop.":
		{
			active = true; 
			let tree = inter.page_name.split("->");

			if(tree[0] == "Inference"){
				w = model.add_object_button(lay,"Init. Pop. Prior",x,y,"AddInitPopPrior",{back:WHITE, active:active, info:info, title:"Initial population prior", te:init_pop_prior_text, siminf:siminf}); x += w+gap;
			}

			w = model.add_object_button(lay,"Init. Pop.",x,y,"AddInitialPopulation",{back:WHITE, active:active, info:info, title:"Initial population", te:init_pop_text, siminf:siminf}); x += w+gap;
		
			if(sp.type != "Individual") active = false;

			w = model.add_object_button(lay,"Add Ind.",x,y,"AddIndividuals",{back:WHITE, active:active, info:info, title:"Add Individuals", te:add_ind_text, siminf:siminf}); x += w+gap;
			
			w = model.add_object_button(lay,"Remove Ind.",x,y,"RemIndividuals",{back:WHITE, active:active, info:info, title:"Remove Individuals", te:rem_ind_text, siminf:siminf}); x += w+gap;
			
			w = model.add_object_button(lay,"Move Ind.",x,y,"MoveIndividuals",{back:WHITE, active:active, info:info, title:"Move Individuals", te:move_ind_text, siminf:siminf}); x += w+gap;
		}
		break;
		
	case "Individual":
		{
			active = true; if(sp.type != "Individual") active = false;
			
			w = model.add_object_button(lay,"Compartment",x,y,"CompData",{ back:WHITE, active:active, info:info, title:"Compartmental data", te:comp_data_text, siminf:siminf}); x += w+gap;
	
			w = model.add_object_button(lay,"Transition",x,y,"TransData",{ back:WHITE, active:active, info:info, title:"Transition data", te:trans_data_text, siminf:siminf}); x += w+gap;

			w = model.add_object_button(lay,"Source",x,y,"SourceData",{ back:WHITE, active:active, info:info, title:"Source data", te:source_data_text, siminf:siminf}); x += w+gap;

			w = model.add_object_button(lay,"Sink",x,y,"SinkData",{ back:WHITE, active:active, info:info, title:"Sink data", te:sink_data_text, siminf:siminf}); x += w+gap;

			w = model.add_object_button(lay,"Diag. Test",x,y,"DiagTestData",{ back:WHITE, active:active, info:info, title:"Diagnostic test data", te:diag_test_data_text, siminf:siminf}); x += w+gap;
		}
		break;

	case "Population":
		{	
			active = true;
			w = model.add_object_button(lay,"Population",x,y,"PopulationData",{ back:WHITE, active:active, info:info, title:"Population data", te:pop_data_text, siminf:siminf}); x += w+gap;
	
			w = model.add_object_button(lay,"Pop. Transition",x,y,"PopTransData",{ back:WHITE, active:active, info:info, title:"Aggregated transition data", te:poptrans_data_text, siminf:siminf}); x += w+gap;
		}
		break;

	case "Capture":
		{
			active = true; if(sp.type != "Individual") active = false;
			w = model.add_object_button(lay,"Set Traps",x,y,"SetTrapsData",{ back:WHITE, active:active, info:info, title:"Set Traps", te:traps_data_text, siminf:siminf}); x += w+gap;
			w = model.add_object_button(lay,"Ind. Trapped",x,y,"IndTrappedData",{ back:WHITE, active:active, info:info, title:"Individuals trapped", te:ind_trap_data_text, siminf:siminf}); x += w+gap;
		}
		break;

	case "Genetic":
		active = true; if(sp.type != "Individual") active = false;
		w = model.add_object_button(lay,"Sequence data",x,y,"SeqData",{ back:WHITE, active:active, info:info, title:"Sequence data", te:seq_data_text, siminf:siminf}); x += w+gap;
		break;

	default: error("Option not recognised 20"); break;		
	}

	add_layer("TableContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-3.5,{table:table});
}


/// Works out if data is still valid given the current model
function check_data_valid(type)
{
	let list = [];
	if(type == "sim" || type == "siminf") list.push("sim");
	if(type == "inf" || type == "siminf") list.push("inf");
	
	
	for(let loop = 0; loop < list.length; loop++){
		
		for(let p = 0; p < model.species.length; p++){
			let sp = model.species[p];

			let source;
			if(list[loop] == "sim") source = sp.sim_source; 
			else source = sp.inf_source; 

			for(let i = 0; i < source.length; i++){
				let so = source[i];
				check_data_source_valid(so);
				
				if(so.error == true){
					model.warn.push({mess:"Data error", mess2:so.error_mess, warn_type:"DataProblem", siminf:type, p:p});
				}
			}
		}
	}
	
	if(type == "inf"){  // If individual data makes sure that add ind is set
		for(let p = 0; p < model.species.length; p++){
			let sp = model.species[p];
			if(find(sp.inf_source,"type","Add Ind.") == undefined){
				for(let i = 0; i < sp.inf_source.length; i++){
					let so = sp.inf_source[i];

					switch(so.type){
					case "Remove Ind.": case "Move Ind.": case "Fixed Eff.": 
					case "Compartment": case "Transition": case "Source": case "Sink":
					case "Diag. Test": case "Ind. Trapped": case "Genetic":
						model.warn.push({mess:"Data error", mess2:"For '"+so.type+"' data to be added individuals must be added to the system using 'Add Ind.'", warn_type:"DataProblem", p:p});
						break;
					default: break;
					}
				}
			}
		}
	}
}


/// Checks if an individual data source is valid
function check_data_source_valid(so)
{
	let p = so.info.p;
	let sp = model.species[p];

	//if(so.info.siminf == "gen"){ so.info.siminf = "inf"; pr("conv");}
	if(so.info.siminf == "gen") sp = sim_result.species[p];
			
	so.error = false;

	let mess_def = "Data source is invalid";

	if(sp.type == "Population"){
		switch(so.type){
		case "Add Ind.":case "Move Ind.":  case "Remove Ind.": case "Compartment":
		case "Transition":case "Source": case "Sink": case "Diag. Test": 
		case "Ind. Trapped": case "Genetic": 
			source_error(so,"Individual-based data '"+so.type+"' cannot be added to a population-based model",300);
			return;
		}
	}
			
	switch(so.type){
	case "Init. Pop.":
		{
			switch(so.spec.radio.value){
			case "Graphical":
				switch(so.spec.radio2.value){
				case "Focal":
					let focal = so.spec.focal.te;
					let cl = find(sp.cla,"name",focal);

					if(cl == undefined){
						source_error(so,"The focal value '"+focal+"' is not a classification",300);
					}
					else{
						if(so.cla.length != sp.ncla) source_error(so,mess_def,1);
						else{
							for(let cl = 0; cl < sp.ncla; cl++){
								let claa = sp.cla[cl];
								let clz = so.cla[cl];

								if(claa.name != clz.cl_name_store) source_error(so,mess_def,2);
								else{
									if(clz.comp_init_pop.length != claa.ncomp) source_error(so,mess_def,3);
									else{
										for(let c = 0; c < claa.ncomp; c++){
											if(clz.comp_init_pop[c].comp_name_store != claa.comp[c].name){
												source_error(so,mess_def,4);
											}
										}
									}
								}
							}
						}
					}
					break;

				case "All":	
					let glob_comp = model.get_glob_comp(p);
					if(glob_comp.length != so.glob_comp.length) source_error(so,mess_def,5);
					else{
						for(let j = 0; j < glob_comp.length; j++){
							let gc = glob_comp[j].cla;
							let gc2 = so.glob_comp[j].cla;
							if(gc.length != gc2.length) source_error(so,mess_def,6);
							else{
								for(let cl = 0; cl < gc.length; cl++){
									if(gc[cl] != gc2[cl]) source_error(so,mess_def,7);
								}
							}
						}
					}
					break;
				
				default: error("Option not recognised 200"); break;
				}
				break;

			case "File":
				let tab = so.table;
				switch(so.spec.radio2.value){
				case "Focal":	
					let focal = so.spec.focal.te;
					let cl2 = find(sp.cla,"name",focal);
					if(cl2 == undefined){
						source_error(so,mess_def,8);
						return "The value '"+focal+"' is not a classification.";
					}
					else{
						let list=[];
						for(let cl = 0; cl < sp.ncla; cl++) list[cl]=[];

						for(let r = 0; r < tab.nrow; r++){
							let val = tab.ele[r][0];
							
							if(val == ""){
								source_error(so,mess_def,8);
								return "The element is empty (col '"+tab.heading[0]+"', row "+(r+2)+").";
							}
							
							let pop = tab.ele[r][1];

							let cl = 0; 
							while(cl < sp.ncla && find(sp.cla[cl].comp,"name",val) == undefined) cl++;
							if(cl == sp.ncla){
								source_error(so,mess_def,9);
								return "The compartment '"+val+"' does not exist (col '"+tab.heading[0]+"', row "+(r+2)+").";
							}
							else{
								if(cl == cl2){
									if(is_percent(pop) == true){
										source_error(so,mess_def,10);
										return "The value '"+pop+"' should not be a percentage (col '"+tab.heading[1]+"', row "+(r+2)+").";
									}
								}
								else{
									if(is_percent(pop) == false){
										source_error(so,mess_def,11);
										return "The value '"+pop+"' should be a percentage (col '"+tab.heading[1]+"', row "+(r+2)+").";	
									}
								}

								if(find_in(list[cl],val) != undefined) source_error(so,mess_def,12);
								list[cl].push(val);
							}
						}

						for(let cl = 0; cl < sp.ncla; cl++){
							if(cl == cl2){
								if(list[cl].length != sp.cla[cl].ncomp){
									source_error(so,mess_def,13);
									return "Not all the compartmental populations are specified in '"+focal+"'"; 
								}
							}	
							else{
								if(list[cl].length < sp.cla[cl].ncomp-1){
									source_error(so,mess_def,14);
									return "All but one compartmental percentages/fractions should be specified in '"+sp.cla[cl].name+"'"; 
								}

								if(list[cl].length == sp.cla[cl].ncomp){
									source_error(so,mess_def,14);
									return "Every compartmental percentages/fractions should not be specified in '"+sp.cla[cl].name+"'"; 
								}
							}
						}
					}
					break;

				case "All":	
					if(tab.ncol != sp.ncla+1) source_error(so,mess_def,15);
					else{
						for(let cl = 0; cl < sp.ncla; cl++){
							let claa = sp.cla[cl];
							if(tab.heading[cl] != claa.name) source_error(so,mess_def,1011);
						}
						let glob_comp = model.get_glob_comp(p);
					
						let flag = [];

						for(let r = 0; r < tab.nrow; r++){
							let j = 0;
							while(j < glob_comp.length){
								let cl = 0; while(cl < sp.ncla && tab.ele[r][cl] == glob_comp[j].cla[cl]) cl++;
								if(cl == sp.ncla) break;
								j++;
							}
							if(j == glob_comp.length){
								source_error(so,mess_def,1012);
								return "The combination of compartments on line "+(r+1)+" of the table is not valid";
							}

							if(flag[j] != undefined){
								source_error(so,mess_def,117);
								return "Line "+(r+1)+" has the same combination of compartments as line "+(flag[j]+1);
							}	
							flag[j] = r;
						}
	
						for(let j = 0; j < glob_comp.length; j++){
							if(flag[j] == undefined){
								source_error(so,mess_def,1006);
								let te = "The combination ";
								for(let cl = 0; cl < sp.ncla; cl++){
									if(cl > 0) te += ",";
									te += "'"+glob_comp[j].cla[cl]+"' ";
								}
								te += "is not set";
								return te;
							}
						}
					}
					break;

				default: error("Option not recognised 201"); break;
				}
				break;

			default: error("Option not recognised 202"); break;	
			}
		}	
		break;

	case "Init. Pop. Prior":
		switch(so.spec.radio.value){
		case "Flat":
			break;

		case "Dirichlet":
			let tab = so.table;
			if(tab.ncol != sp.ncla+1) source_error(so,mess_def,1500);
			else{
				for(let cl = 0; cl < sp.ncla; cl++){
					let claa = sp.cla[cl];
					if(tab.heading[cl] != claa.name) source_error(so,mess_def,1007);	
				}
				let glob_comp = model.get_glob_comp(p);
			
				let flag = [];

				for(let r = 0; r < tab.nrow; r++){
					let j = 0;
					while(j < glob_comp.length){
						let cl = 0; while(cl < sp.ncla && tab.ele[r][cl] == glob_comp[j].cla[cl]) cl++;
						if(cl == sp.ncla) break;
						j++;
					}
					if(j == glob_comp.length){
						source_error(so,mess_def,1008);
						return "The combination of compartments on line "+(r+1)+" of the table is not valid";
					}

					if(flag[j] != undefined){
						source_error(so,mess_def,1009);
						return "Line "+(r+1)+" has the same combination of compartments as line "+(flag[j]+1);
					}	
					flag[j] = r;
				}

				for(let j = 0; j < glob_comp.length; j++){
					if(flag[j] == undefined){
						source_error(so,mess_def,1010);
						let te = "The combination ";
						for(let cl = 0; cl < sp.ncla; cl++){
							if(cl > 0) te += ",";
							te += "'"+glob_comp[j].cla[cl]+"' ";
						}
						te += "is not set";
						return te;
					}
				}
			}
			break;
		}
		break;

	case "Population":
		{
			let filt = so.spec.filter;
			if(filt.cla.length != sp.ncla) source_error(so,mess_def,21);
			else{
				for(let cl = 0; cl < sp.ncla; cl++){
					let claa = sp.cla[cl];
					let clz = filt.cla[cl];
					if(claa.name != clz.cl_name_store) source_error(so,mess_def,22);
					else{
						if(clz.comp.length != claa.ncomp) source_error(so,mess_def,23);
						else{
							for(let c = 0; c < claa.ncomp; c++){
								if(clz.comp[c].comp_name_store != claa.comp[c].name) source_error(so,mess_def,24);
							}
						}
					}
				}
			}
		}
		break;

	case "Add Ind.":
		{
			let tab = so.table;
			if(tab.ncol != sp.ncla+2) source_error(so,mess_def,25);
			else{
				for(let cl = 0; cl < sp.ncla; cl++){
					let claa = sp.cla[cl];

					if(tab.heading[cl+2] != claa.name) source_error(so,mess_def,26);
					else{
						for(let r = 0; r < tab.nrow; r++){
							let val = tab.ele[r][cl+2];
							if(find(claa.comp,"name",val) == undefined) source_error(so,mess_def,27);
						}
					}
				}
			}
		}
		break;

	case "Move Ind.":
		{
			let cl = find(sp.cla,"name",so.spec.cl_drop.te);
			if(cl == undefined) source_error(so,mess_def,28);
			else{
				let claa = sp.cla[cl];

				let tab = so.table;
				for(let r = 0; r < tab.nrow; r++){
					let to = tab.ele[r][2];
					if(find(claa.comp,"name",to) == undefined){
						source_error(so,mess_def,30);
						return "The value '"+to+"' is not a compartment in classification '"+claa.name+"' (col '"+tab.heading[3]+"', row "+(r+2)+").";
					}
				}
			}
		}
		break;

	case "Compartment":
		{
			let cl = find(sp.cla,"name",so.spec.cl_drop.te);
			if(cl == undefined) source_error(so,mess_def,31);
			else{
				let claa = sp.cla[cl];

				let tab = so.table;
				if(tab.heading[2] != claa.name) source_error(so,mess_def,32);
				else{
				}
			}
		}
		break;

	case "Transition":
		{
			let tr_name = so.spec.tr_drop.te;

			let cl = 0; 
			while(cl < sp.ncla && find(sp.cla[cl].tra,"name",tr_name) == undefined){
				cl++;
			}
			if(cl == sp.ncla){	
				source_error(so,mess_def,33);
				return "Transition '"+tr_name+"' is not found";
			}
		}
		break;

	case "Source":
		{
			let tr_name = so.spec.tr_drop.te;

			let flag = false;
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				for(let j = 0; j < claa.tra.length; j++){
					let tra = claa.tra[j];
					if(tra.i == "Source"){
						if(tr_name == claa.comp[tra.f].name) flag = true;
					}
				}
			}

			if(flag == false) source_error(so,mess_def,34);
		}
		break;

	case "Sink":
		{
			let tr_name = so.spec.tr_drop.te;

			let flag = false;
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				for(let j = 0; j < claa.tra.length; j++){
					let tra = claa.tra[j];
					if(tra.f == "Sink"){
						if(tr_name == claa.comp[tra.i].name) flag = true;
					}
				}
			}

			if(flag == false) source_error(so,mess_def,35);
		}
		break;

	case "Diag. Test":
		{
			let cb = so.spec.check_box;
			let cl = find(sp.cla,"name",cb.name);
			if(cl == undefined) source_error(so,mess_def,36);
			else{
				let claa = sp.cla[cl];
				if(cb.value.length != claa.ncomp) source_error(so,mess_def,37);
				else{
					for(let c = 0; c < claa.ncomp; c++){
						if(cb.value[c].comp_name_store != claa.comp[c].name) source_error(so,mess_def,38);
					}
				}
			}
		}
		break;

	case "Pop. Trans.":
		{
			let cl = find(sp.cla,"name",so.spec.cl_drop.te);
			if(cl == undefined) source_error(so,mess_def,39);
			else{
				let claa = sp.cla[cl];

				let filt = so.spec.filter;

				if(filt.tra.length != claa.ntra) source_error(so,mess_def,40); 
				else{
					let flag = false;
					for(let j = 0; j < claa.ntra; j++){
						let tra = filt.tra[j];
					
						if(tra.tra_name_store != claa.tra[j].name) source_error(so,mess_def,41);
						if(tra.check == true) flag = true;
					}
					if(flag == false) source_error(so,mess_def,142);
				}

				if(filt.cla.length != sp.ncla){
					pr(filt.cla.length +" "+ sp.ncla+" filt");
					source_error(so,mess_def,42); 
				}
				else{
					for(let cl2 = 0; cl2 < sp.ncla; cl2++){
						if(cl2 != cl){
							let claa2 = sp.cla[cl2];
							let clz = filt.cla[cl2];
							if(clz.cl_name_store != claa2.name) source_error(so,mess_def,43);
							else{	
								if(claa2.ncomp != clz.comp.length) source_error(so,mess_def,44);
								else{
									for(let c = 0; c < claa2.ncomp; c++){
										if(clz.comp[c].comp_name_store != claa2.comp[c].name){
											source_error(so,mess_def,45);
										}	
									}		
								}
							}
						}
					}
				}
			}
		}
		break;

	case "Set Traps":
		{
			let filt = so.spec.filter;
			if(filt.cla.length != sp.ncla) source_error(so,mess_def,46);
			else{
				for(let cl = 0; cl < sp.ncla; cl++){
					let claa = sp.cla[cl];
					let clz = filt.cla[cl];
					if(claa.name != clz.cl_name_store) source_error(so,mess_def,47);
					else{
						if(clz.comp.length != claa.ncomp) source_error(so,mess_def,48);
						else{
							for(let c = 0; c < claa.ncomp; c++){
								if(clz.comp[c].comp_name_store != claa.comp[c].name) source_error(so,mess_def,49);
							}
						}
					}
				}
			}
		}
		break;

	case "Remove Ind.": break;
	case "Ind. Trapped": break;
	case "Genetic": break;

	case "LoadTensor": break;
	case "LoadReparam": break;
	case "LoadPriorSplit": break;
	case "LoadDistSplit": break;
	case "LoadXvector2": break;
	case "Fixed Effect": break;
	case "Trans exp(rate)": break;
	
	default: error("Source not checked error: "+so.type); break;
	}
}


/// Sets an error in a source
function source_error(so,te,num)
{
	so.error = true;
	if(true){ error("Source error: "+num);}
	if(te != undefined) so.error_mess = te;
	else so.error_mess = "An unspecified error occurred";
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
	close_bubble();
	
	let tree = inter.page_name.split("->");
	if(tree[1] == "Generate Data") sim_data(); 
	else load_data_table_initialise();
}


/// Sets the columns which need to be loaded 
function set_loadcol()
{
	let so = inter.edit_source;

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
	if(type == "Transition" || type == "Source" || type == "Sink"){
		let time_radio = inter.edit_source.spec.time_radio.value;
			
		switch(time_radio){
		case "All": case "Spec": type = "Transition"; break;
		case "File": 
			switch(type){
			case "Transition": type = "Trans. in time range"; break;
			case "Source": type = "Source in time range"; break;
			case "Sink": type = "Sink in time range"; break;
			}
			break;
		default: error("Can't find option"); break;
		}
	}
	
	let i = find(data_template,"type",type);
	
	if(i == undefined){ error("Cannot find data template"); return;}

	let temp = copy(data_template[i]);

	let cols = temp.cols;
	for(let c = 0; c < cols.length; c++){
		let spl = cols[c].split(",");

		switch(spl[0]){
			case "ID":
				load_col.push({heading:"ID", desc:"individual ID",type:"text"});
				break;

			case "t":
				load_col.push({heading:"t", desc:spl[1], type:"time_float"});
				break;

			case "ttrans":
				load_col.push({heading:"t", desc:spl[1], type:"time_trans_float"});
				break;
				
			case "snp":
				load_col.push({heading:"SNP", desc:spl[1], type:"text"});
				break;

			case "cl_prob":
				{
					let p = so.info.p;
					if(p == undefined){ error("p should be defined"); return;}
						
					let sp = model.species[p]; 
					
					let cl_drop = so.spec.cl_drop;
					let cl = find_nocase(sp.cla,"name",cl_drop.te);
					if(cl == undefined){ error("cl should be defined"); return;}
						
					load_col.push({heading:cl_drop.te, desc:"the classification in '"+cl_drop.te+"'",type:"compartment_prob", p:p, cl:cl});
				}
				break;

			case "cl_all":
				{
					let p = so.info.p;
					if(p == undefined){ error("p should be defined"); return;}
						
					if(so.type == "Init. Pop." && so.spec.radio2 != undefined && 
							so.spec.radio2.value == "Focal"){
						load_col.push({heading:"Compartment", desc:"Compartment", type:"comp_all", p:p});
					}
					else{
						let sp = model.species[p];
						for(let cl = 0; cl < sp.ncla; cl++){
							load_col.push({heading:sp.cla[cl].name, desc:"compartment in '"+sp.cla[cl].name+"' individuals are added",type:"compartment", p:p, cl:cl});
						}
					}
				}
				break;

			case "from":
				{
					let p = so.info.p;
					if(p == undefined){ error("p should be defined"); return;}
					let cl_name = so.spec.cl_drop.te;
					let cl = find_nocase(model.species[p].cla,"name",cl_name);
					
					load_col.push({heading:"From", desc:"the compartment in '"+cl_name+"' from which individuals come", type:"comp_source_sink", p:p, cl:cl});
				}
				break;

			case "to":
				{
					let p = so.info.p;
					if(p == undefined){ error("p should be defined"); return;}
					let cl_name = so.spec.cl_drop.te;
					let cl = find_nocase(model.species[p].cla,"name",cl_name);
					
					load_col.push({heading:cl_name, desc:"the compartment in '"+cl_name+"' to which individuals go", type:"comp_source_sink", p:p, cl:cl});
				}
				break;

			case "pop":
				if(so.type == "Init. Pop." && so.spec.radio2.value == "Focal"){
					load_col.push({heading:"Population", desc:"the population", type:"positive float/percent"});
				}
				else{
					load_col.push({heading:"Population", desc:"the population", type:"pos_zero_float"});
				}
				break;

			case "alpha":
				load_col.push({heading:"Alpha", desc:"the alpha value", type:"pos_float"});
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
					let spec = inter.edit_source.spec;
			
					let filt = spec.filter;
					let p = so.info.p;

					let sp = model.species[p];
					
					for(let cl = 0; cl < sp.ncla; cl++){
						if(filt.cla[cl].radio.value == "File"){
							load_col.push({heading:sp.cla[cl].name, desc:"'"+sp.cla[cl].name+"' which specifies sub-population",type:"multi_comp", p:p, cl:cl});
						}
					}
				
					load_col.push({heading:"Population", desc:"the observed population estimates", type:"pos_zero_float"});
					
					if(spec.obs_error.value == "file"){
						load_col.push({heading:"sd", desc:"the standard deviation in population estimates", type:"pos_zero_float"});
					}
				}
				break;
				
			case "filt_obspoptrans":
				{	
					let spec = inter.edit_source.spec;
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
										load_col.push({heading:sp.cla[cl2].name, desc:"'"+sp.cla[cl2].name+"' which specifies sub-population",type:"compartment", p:p, cl:cl2});
									}
								}
							}
							
							load_col.push({heading:"Number", desc:"the observed number of transitions", type:"pos_zero_float"});
							
							if(spec.obs_error.value == "file"){
								load_col.push({heading:"sd", desc:"the standard deviation in population estimates", type:"pos_float"});
							}
						}
					}
				}
				break;
			
			case "filt_settraps":
				{
					let spec = inter.edit_source.spec;
			
					let filt = spec.filter;
					let p = so.info.p;
					let sp = model.species[p];
					
					for(let cl = 0; cl < sp.ncla; cl++){
						if(filt.cla[cl].radio.value == "File"){
							load_col.push({heading:sp.cla[cl].name, desc:"'"+sp.cla[cl].name+"' which specifies sub-population",type:"compartment", p:p, cl:cl});
						}
					}
				}
				break;

			case "setrap_name":
				load_col.push({heading:"Name", desc:"the unique name for the trapping event", type:"text"});
				break;

			case "comp_name":
				load_col.push({heading:"Name", desc:"the compartment name", type:"text"});
				break;
				
			case "start_comp":
				load_col.push({heading:"From", desc:"the compartment from which individuals move", p:so.info.p, cl:so.info.cl, type:"compartment"});
				break;
				
			case "end_comp":
				load_col.push({heading:"To", desc:"the compartment to which individuals move", p:so.info.p, cl:so.info.cl, type:"compartment"});
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
				load_col.push({heading:"Colour", desc:"the compartment colour", type:"colour"});
				break;
				
			case "rate":
				load_col.push({heading:"Rate", desc:"the transition rate", p:so.info.p, cl:so.info.cl, type:"eqn"});
				break;
			
			case "mean":
				load_col.push({heading:"Mean", desc:"the transition mean time", p:so.info.p, cl:so.info.cl, type:"eqn"});
				break;
		
			case "sd":
				load_col.push({heading:"SD", desc:"the transition standard deviation", p:so.info.p, cl:so.info.cl, type:"eqn"});
				break;
			
			case "shape":
				load_col.push({heading:"Shape", desc:"the transition shape parameter", p:so.info.p, cl:so.info.cl, type:"eqn"});
				break;
			
			case "scale":
				load_col.push({heading:"Scale", desc:"the transition scale parameter", p:so.info.p, cl:so.info.cl, type:"eqn"});
				break;
				
			case "boundary":
				load_col.push({heading:"Boundary", desc:"", type:"boundary"});
				break;
				
			case "dep":
				let dep = model.param[inter.edit_param.i].dep;
				let desc = convert_tensor_text("the tensor");
				desc += " index giving '";
				for(let k = 0; k < dep.length; k++){
					let res = find_cla_from_index(remove_prime(dep[k]));;
					if(res != undefined){
						let claa = model.species[res.p].cla[res.cl];
						load_col.push({heading:dep[k], desc:desc+claa.name+"'", p:res.p, cl:res.cl, type:"compartment"});
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
					load_col.push({heading:"Dist", desc:"the element distribution", pri_pos:par.pri_pos, type:"prior"});
				}
				break;
				
			case "eqn":
				load_col.push({heading:"Eqn", desc:"the equation", type:"eqn"});
				break;
				
			case "A":
				//load_col.push({heading:"Eqn", desc:"the equation", type:"eqn"});
				break;
				
			default: error(spl[0]+" in data column not recognised"); break;
		}
	}
	
	so.load_col = load_col;
}

function data_source_add_and_show()
{
	data_source("Add");
	generate_screen();
	scroll_view_last_button("TableContent");
}
	
	
/// Allows data sources to be added and editted
function data_source(type,info)
{
	if(info == undefined) info = inter.edit_source.info;
	let p = info.p; 
	if(p == undefined){ error("Data source does not have p"); return}
	
	let sp = model.species[p];
	if(info.siminf == "gen") sp = sim_result.species[p];
	
	let so;
	switch(type){
	case "Add":
		so = copy(inter.edit_source);
		add_source_description(so);
		break;
		
	case "Edit": 
		break;
	
	case "Replace": 
		so = copy(inter.edit_source);
		so.edit_spec = false;
		add_source_description(so);
		info = so.info;
		break;
		
	default: error("Error option"); break;
	}

	switch(info.siminf){
	case "sim": 
		switch(type){
		case "Add": sp.sim_source.push(so); break;
		case "Edit": so = sp.sim_source[info.i]; break;
		case "Replace": sp.sim_source[info.i] = so; break;
		default: error("Error option"); break;
		}
		break;
		
	case "inf": 
		switch(type){
		case "Add": sp.inf_source.push(so); break;
		case "Edit": so = sp.inf_source[info.i]; break;
		case "Replace": sp.inf_source[info.i] = so; break;
		default: error("Error option"); break;
		}
		break;
		
	case "gen": 
		switch(type){
		case "Add": sp.gen_source.push(so); break;
		case "Edit": so = sp.gen_source[info.i]; break;
		case "Replace": sp.gen_source[info.i] = so; break;
		default: error("Error option"); break;
		}
		break;
		
	default: error("Option not recognised 21"); break;
	}
		
	switch(type){
	case "Add": 
		close_data_source();
		break;

	case "Replace":
		close_data_source();
		break;
		
	case "Edit":
		inter.edit_source = copy(so);
		inter.edit_source.info = info;
		break;

	default: error("Option not recognised 22"); break;
	}
}


/// Adds a description of the data source
function add_source_description(so)
{
	let desc="";

	let sp = model.species[so.info.p];
	let spec = so.spec;
	
	//pr(so.type+"type");
	switch(so.type){
	case "Init. Pop.":
		switch(spec.radio.value){
		case "File":
			let tab = so.table;
			let num = 0;
			
			switch(spec.radio2.value){
			case "Focal":
				let focal = so.spec.focal.te;
				let cl = find(sp.cla,"name",focal);
				if(cl == undefined){
					if(so.error != true) error("Focal not recognised");
				}
				else{
					let claa = sp.cla[cl];
				
					desc += focal+"  ";

					for(let r = 0; r < tab.nrow; r++){
						if(find(claa.comp,"name",tab.ele[r][0]) != undefined){
							num += Number(tab.ele[r][1]);
						}
					}
				}
				break;
			
			case "All":
				let col = tab.ncol-1;
				for(let r = 0; r < tab.nrow; r++){
					num += Number(tab.ele[r][col]);
				}
				break;
				
			default: error("Option not reco"); break;
			}
				
			desc += "# Ind. "+num;
			break;
			
		case "Graphical":
			switch(so.spec.radio2.value){
			case "Focal":
				{	
					let focal = so.spec.focal.te;
					let cl = find(sp.cla,"name",focal);
					if(cl == undefined){
						desc += "'"+focal+"' no longer exists";
					}
					else{
						let cip = so.cla[cl].comp_init_pop;
						let sum = 0;
						for(let c = 0; c < cip.length; c++) sum += cip[c].pop;
						desc += "# Ind. "+sum+" ";
						if(sp.ncla > 1) desc += "using '"+sp.cla[cl].name+"'";
					}
				}
				break;
				
			case "All":
				let sum = 0;
				for(let c = 0; c < so.glob_comp.length; c++){
					sum += so.glob_comp[c].pop;
				}
				desc += "# Ind. "+sum;
				//if(pop.ncla > 1) desc += "All compartments specified.";
				break;
			}
		}
		break;

	case "Init. Pop. Prior":
		switch(so.variety){
		case "Flat": desc += "Flat prior"; break;
		case "Dirichlet": desc += "Dirichlet prior"; break;
		}	
		break;
		
	case "Add Ind.": case "Remove Ind.": 
		desc = "# Ind. "+so.table.nrow;
		break;

	case "Move Ind.":
		desc = "# Ind. moves: "+so.table.nrow;
		break;

	case "Compartment": case "Diag. Test":
		{
			let cl_se = find(sp.cla,"name",spec.cl_drop.te);
			desc = so.spec.cl_drop.te+"  # Obs. "+so.table.nrow;
		}
		break;
		
	case "Transition": case "Source": case "Sink":
		desc = so.spec.tr_drop.te+"  # Obs. "+so.table.nrow;
		break;
		
	case "Population":	
		desc = comp_filter_desc(so.spec.filter.cla,sp)+"# Obs. "+so.table.nrow;
		break;
		
	case "Pop. Trans.":
		{
			let trz = so.spec.filter.tra;
			
			let cl_se = find(sp.cla,"name",spec.cl_drop.te);
			if(cl_se == undefined) desc = "'"+spec.cl_drop.te+"' no longer exists";
			else{
				let tra = sp.cla[cl_se].tra;
				let st = "";
				for(let j = 0; j < trz.length; j++){
					if(trz[j].check == true){
						if(st != "") st += "|";
						st += tra[j].name;
					}
				}
				if(st != "") desc += st+"  ";
				
				let te = comp_filter_desc(so.spec.filter.cla,sp);

				if(te != "") desc += te+"  ";
				
				desc += "# Obs. "+so.table.nrow;	
			}
		}
		break;
	
	case "Set Traps":
		desc = comp_filter_desc(so.spec.filter.cla,sp)+"# Trap events "+so.table.nrow;
		break;

	case "Ind. Trapped":
		desc = "# Trapped "+so.table.nrow;
		break;

	case "Genetic":
		desc = "# SNP "+(so.table.ncol-2)+"  # Seq. "+so.table.nrow;
		break;

	default: error(so.type+" DESC NOT ADDED"); break;
	}
	
	so.desc = desc;
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
				let te2 = "";
				for(let j = 0; j < clzz.comp.length; j++){
					if(clzz.comp[j].check == true){
						if(te2 != "") te2 += "|";
						te2 += sp.cla[cl].comp[j].name;
					}
				}
				te += te2+" ";
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

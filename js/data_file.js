"use strict";
// Functions relating to the loading of data files
 
//////////////////////// INIT POP

function initpop_data_bubble(cont,type)
{
	let sp = model.get_sp();
	
	let spec = edit_source.spec;
	
	cont.dx = 10;
	
	bubble_addtitle(cont,"Initial population",{te:init_pop_text});
	
	let warn = "";
	if(sp.ncla == 0) warn = "First a classification must be added to '"+sp.name+"'.";
	else{
		for(let cl = 0; cl < sp.ncla; cl++){
			if(sp.cla[cl].ncomp == 0){
				warn = "First compartments need to be added to classification '"+sp.cla[cl].name+"'.";
			}
		}
	}
	
	if(warn != 0){
		bubble_addparagraph(cont,warn,0,10); cont.y += 0.2;
		add_end_button(cont,"Close","CloseBubble");
	}
	else{
		let disable = false; if(type == "view") disable = true;
		
		bubble_addradio(cont,0,"Fixed","Fixed",spec.radio_dist,disable);
		bubble_addradio(cont,0,"Dist","Distribution",spec.radio_dist,disable);
		cont.y += 0.2;
	
		let ac = "AddInitialPopulation2";
		
		if(!(sp.ncla ==	1 && spec.radio_dist.value == "Fixed")){
			cont.y += 0.3;
			bubble_addparagraph(cont,"Defined by:",0,10); cont.y += 0.2;
				
			if(sp.ncla ==	1){
				bubble_addradio(cont,0,"Focal","Compartments",spec.radio2,disable);
				bubble_addradio(cont,0,"All","Total population",spec.radio2,disable);
			}
			else{
				bubble_addradio(cont,0,"Focal","Focal classification",spec.radio2,disable);
				if(spec.radio2.value == "Focal"){
					let cl_pos = get_pos(sp.cla,"name");
					let wmax = get_w_drop(cl_pos);
					if(type == "view") cl_pos = undefined;
					if(wmax+2.2 > cont.dx) cont.dx = wmax+2.2;
					bubble_adddropdown(cont,2.2,wmax,spec.focal,cl_pos);
					if(spec.focal.te == "Select") ac = undefined;
					cont.y += 0.6;
				}
				bubble_addradio(cont,0,"All","All compartments",spec.radio2,disable);
			}
		}
		
		cont.y += 0.3;
	
		if(disable == false){
			bubble_addparagraph(cont,"Set using:",0,10); cont.y += 0.2;
			bubble_addradio(cont,0,"Graphical","Graphical interface",spec.radio,disable);
			bubble_addradio(cont,0,"File","File",spec.radio,disable);
		}
		
		if(type == "view") add_end_button(cont,"Done","DataTableBack");
		else add_end_button(cont,"Next",ac);
	}		
}


/// Initialises init_pop data
function initpop_data(bu)
{
	start_data_source("Init. Pop.",{radio_dist:{value:"Fixed"}, radio:{value:"Graphical"},radio2:{value:"Focal"}, focal:{te:"Select"}},bu.op.info);
	
	if(bu.op.siminf == "gen"){
		let spec = edit_source.spec;
		spec.radio_dist.value = "Fixed";
		spec.radio.value = "File";
		spec.radio2.value = "All";
		set_loadcol();
		sim_data(); 
	}
	else select_bubble_over();
}


/// Sets up after click
function initpop_data2()
{
	let spec = edit_source.spec;
	
	let sp = model.get_sp(); 
	
	if(sp.ncla ==	1){
		spec.focal.te = sp.cla[0].name;
		if(spec.radio_dist.value == "Fixed") spec.radio2.value = "Focal";
	}
	
	set_loadcol();
	
	switch(spec.radio.value){
	case "Graphical":
		start_add_initial_population(edit_source);
		delete edit_source.load_datatable;
		delete edit_source.table_loaded;
		break;

	case "File":
		file_add_datatable();
		break;

	default: error("Option not recognised 24"); break;
	}
	close_bubble();	
}


/// Start editing of initial population
function start_add_initial_population(so)
{
	let p = so.info.p;
	let sp = model.species[p]; 
	
	switch(so.spec.radio2.value){
	case "Focal":
		so.cla = [];
	
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			
			switch(so.spec.radio_dist.value){
			case "Fixed":
				{
					so.cla[cl] = { cl_name_store:claa.name, comp_init_pop:[]};

					let per = 100/claa.ncomp;
					for(let c = 0; c < claa.ncomp; c++){
						let co = claa.comp[c];
						so.cla[cl].comp_init_pop[c] = { comp_name_store:co.name, pop:0, pop_per:per};
					}
				}
				break;
				
			case "Dist":
				{
					so.cla[cl] = { cl_name_store:claa.name, comp_init_dist:[]};

					for(let c = 0; c < claa.ncomp; c++){
						let co = claa.comp[c];
						so.cla[cl].comp_init_dist[c] = { comp_name_store:co.name, dist:unset_prior("pop_prior"), alpha:1};
					}
				}
				break;
			}
		}
		break;
		
	case "All":	
		so.glob_comp = get_glob_comp(p);
	
		switch(so.spec.radio_dist.value){
		case "Fixed":
			for(let c = 0; c < so.glob_comp.length; c++){
				so.glob_comp[c].pop = 0;
			}
			break;
		
		case "Dist":
			so.pop_dist = unset_prior("pop_prior");
			for(let c = 0; c < so.glob_comp.length; c++){
				so.glob_comp[c].alpha = 1;
			}
			break;
		}
		break;

	default: error("Option not recognised 25"); break;
	}
}


/// Sets value to undefined (these are to be loaded from a file)
function reset_add_initial_population(so)
{
	let p = so.info.p;
	let sp = model.species[p]; 
	
	switch(so.spec.radio2.value){
	case "Focal":
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			let clas = so.cla[cl];
		
			switch(so.spec.radio_dist.value){
			case "Fixed":
				{
					for(let c = 0; c < claa.ncomp; c++){
						clas.comp_init_pop[c].pop = undefined;
						clas.comp_init_pop[c].pop_per = undefined;
					}
				}
				break;
				
			case "Dist":
				{
					for(let c = 0; c < claa.ncomp; c++){
						clas.comp_init_dist[c].alpha = undefined; 
					}
				}
				break;
			}
		}
		break;
		
	case "All":	
		so.glob_comp = get_glob_comp(p);
	
		switch(so.spec.radio_dist.value){
		case "Fixed":
			for(let c = 0; c < so.glob_comp.length; c++){
				so.glob_comp[c].pop = undefined;
			}
			break;
		
		case "Dist":
			for(let c = 0; c < so.glob_comp.length; c++){
				so.glob_comp[c].alpha = undefined;
			}
			break;
		}
		break;

	default: error("Option not recognised 25"); break;
	}
}


/// Fills in unspecified percentages
function fill_in_percentages(so)
{
	let p = so.info.p;
	let sp = model.species[p]; 

	switch(so.spec.radio2.value){
	case "Focal":
		{
			let foc_cl = so.spec.focal.te;
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				let clas = so.cla[cl];
				
				switch(so.spec.radio_dist.value){
				case "Fixed":
					{
						if(claa.name == foc_cl){
							for(let c = 0; c < claa.ncomp; c++){
								if(clas.comp_init_pop[c].pop == undefined){
									clas.comp_init_pop[c].pop = 0;
									//alert_import("Initial population for '"+clas.comp_init_pop[c].comp_name_store+"' not set");
								}
							}
						}
						else{
							let sum = 0;
							let list = [];
							for(let c = 0; c < claa.ncomp; c++){
								if(clas.comp_init_pop[c].pop_per != undefined) sum += clas.comp_init_pop[c].pop_per;
								else{
									list.push(c);
								}
							}
							if(list.length != 1){
								alert_import("In classification '"+claa.name+"' all but one percentage must be specified");
							}
							if(sum > 100) alert_import("In classification '"+claa.name+"' the percentages add up to more than 100.");
							clas.comp_init_pop[list[0]].pop_per = 100-sum;
						}
					}
					break;
					
				case "Dist":
					{
						if(claa.name == foc_cl){
							for(let c = 0; c < claa.ncomp; c++){
								if(clas.comp_init_dist[c].dist.type.te == select_str){
									alert_import("The prior for the population in '"+clas.comp_init_dist[c].comp_name_store+"' must be set");
								}
							}
						}
						else{
							for(let c = 0; c < claa.ncomp; c++){
								let co = claa.comp[c];
								if(clas.comp_init_dist[c].alpha == undefined){
									alert_import("Alpha value for '"+clas.comp_init_dist[c].comp_name_store+"' not set");
								}
							}
						}
					}
					break;
				}
			}
		}
		break;
		
	case "All":	
		switch(so.spec.radio_dist.value){
		case "Fixed":
			for(let c = 0; c < so.glob_comp.length; c++){
				if(so.glob_comp[c].pop == undefined){
					so.glob_comp[c].pop = 0;
					//alert_import("Population for '"+get_glob_comp_name(so.glob_comp[c])+"' not set");
				}
			}
			break;
		
		case "Dist":
			for(let c = 0; c < so.glob_comp.length; c++){
				if(so.glob_comp[c].alpha == undefined){
					alert_import("Population for '"+get_glob_comp_name(so.glob_comp[c])+"' not set");
				}
			}
			break;
		}
		break;

	default: error("Option not recognised 25"); break;
	}
}


/// Returns the name of the global compartment
function get_glob_comp_name(gc)
{
	let te = "";
	let list = gc.cla;
	for(let i = 0; i < list.length; i++){
		if(te != "") te += "|"; 
		te += list[i];
	}
	return te;
}


/// Converts from table to graphical interface 
function init_pop_convert_to_graphical(so,imp)
{
	start_add_initial_population(so);
	
	reset_add_initial_population(so);
			
	let tab = so.table;

	let p = so.info.p;
	let sp = model.species[p]; 

	so.spec.radio = {value: "Graphical"};
	so.table = undefined;

	switch(so.spec.radio2.value){
	case "Focal":
		{	
			for(let r = 0; r < tab.nrow; r++){
				let name = tab.ele[r][0];
				let val = tab.ele[r][1];
				
				let flag = false;
				for(let cl = 0; cl < sp.ncla; cl++){
					let claa = sp.cla[cl];
					for(let c = 0; c < claa.ncomp; c++){
						switch(so.spec.radio_dist.value){
						case "Fixed":
							{
								let cip = so.cla[cl].comp_init_pop[c]
								if(cip.comp_name_store == name){
									if(so.spec.focal.te == claa.name) cip.pop = Number(val);
									else{
										if(is_percent(val) == true) val = Number(val.substr(0,val.length-1));
										cip.pop_per = Number(val);
									}
									flag = true;
								}
							}
							break;
								
						case "Dist":
							{
								let cip = so.cla[cl].comp_init_dist[c];
								if(cip.comp_name_store == name){
									if(so.spec.focal.te == claa.name){
										let pri = convert_text_to_prior(val,prior_pos_positive,false);
										if(pri.err == true){
											let te = "Problem converting prior string '"+val+"':"+pri.msg;
											if(imp) alert_import(te); 
											else alert_help("Prior problem",te);
										}
										cip.dist = pri;
									}
									else{
										cip.alpha = Number(val);
									}
									flag = true;
								}
							}
							break;
						}
					}
				}
				
				if(flag == false){ alert_help("Could not find"); return;}
			}
		}
		break;
		
	case "All":	
		{
			let gc = so.glob_comp;
			
			for(let r = 0; r < tab.nrow; r++){
				let row = tab.ele[r];
				let c_gl;
				for(c_gl = 0; c_gl < gc.length; c_gl++){
					let cl = 0; while(cl < sp.cla.length && gc[c_gl].cla[cl] == row[cl]) cl++;
					if(cl == sp.cla.length){
						switch(so.spec.radio_dist.value){
						case "Fixed": gc[c_gl].pop = Number(row[cl]); break;
						case "Dist": gc[c_gl].alpha = Number(row[cl]); break;
						}
						break;
					}
				}
				if(c_gl == gc.length){ alert_help("Could not find"); return;}
			}
		}
		break;
	}

	fill_in_percentages(so);
	
	delete so.load_datatable;
	delete so.table_loaded;
}


/////////////////////// MOVE INDIVIUDUALS

function move_data_bubble(cont,type)
{
	let sp = model.get_sp();
	let cl_pos;
	if(type != "view") cl_pos = get_pos(sp.cla,"name");

	let wmax = get_w_drop(cl_pos);
	let cl_drop = edit_source.spec.cl_drop;

	cont.dx = 5.2+wmax;

	bubble_addtitle(cont,"Move individuals",{te:move_ind_text});
	bubble_addparagraph(cont,"The classification in which individuals move compartment must be identified.",0,cont.dx); cont.y += 0.2;
	
	bubble_adddropdown(cont,5.2,wmax,cl_drop,cl_pos); cont.y -= 1.0;
	bubble_addparagraph(cont,"Classification:",0,10); cont.y += 0.4;
	cont.y += 0.3;
	
	if(type != "view"){
		let but_name = "Next"; if(type == "view") but_name = "Done";

		let ac; 
		if(cl_drop.te != select_drop_str){
			ac = "AddDataTable"; if(type == "view") ac = "DataTableBack";
		}	
		
		add_end_button(cont,but_name,ac);	
	}
}


/// Loads up the move data bubble
function move_data(bu)
{
	start_data_source("Move Ind.",{p:model.get_p(), cl_drop:{te:select_drop_str}},bu.op.info);
	select_bubble_over();
}


//////////////////////// COMP Data

function comp_data_bubble(cont,type)
{
	let cl_pos = get_cl_pos();
	
	let wmax = get_w_drop(cl_pos);
	let cl_drop = edit_source.spec.cl_drop;

	cont.dx = 10; if(cont.dx < wmax+5) cont.dx = wmax+5;

	if(type == "view") cl_pos = undefined;
	
	bubble_addtitle(cont,"Compartmental data",{te:comp_data_text});
	
	if(cl_pos && cl_pos.length == 1){
		cl_drop.te = cl_pos[0].te;
	}
		
	bubble_addparagraph(cont,"Select the classification on which measurements are made.",0,cont.dx+2); 
	cont.y += 0.2;

	bubble_adddropdown_text(cont,"Classification:",5.5,wmax,cl_drop,cl_pos); 
	cont.y += 0.5;
	
	if(sim_options()){
		add_timepoint_options(cont);
		
		cont.y += 0.6;

		let so = edit_source;
		bubble_add_minititle(cont,"Observations:"); 
		cont.y -= 0.4;

		bubble_addradio(cont,0,"Exact","Exact",so.exact); cont.y -= 1.4;
		bubble_addradio(cont,6,"Noisy","Noisy",so.exact);

		if(so.exact.value == "Noisy"){
			bubble_input(cont,"Accuracy",{x:1, w:cont.dx-1, type:"comp_acc"});
		}
	}
	
	if(type == "view") add_end_button(cont,"Done","DataTableBack");
	else{
		let ac; if(cl_drop.te != select_drop_str) ac = "AddDataTable";
		add_end_button(cont,"Next",ac);	
	}
}


/// Adds options for time points (used to generate simulated data)
function add_timepoint_options(cont)
{
	cont.y += 0.6;

	let so = edit_source;
	bubble_add_minititle(cont,"Times of observations:"); 
	cont.y -= 0.4;

	bubble_addradio(cont,0,"Periodic","Periodic",so.time_radio); cont.y -= 1.4;
	bubble_addradio(cont,6,"Specified","Specified",so.time_radio);

	let name;
	switch(so.time_radio.value){
	case "Periodic": name = "Time-step"; break;
	case "Specified": name = "Times (comma separated)"; break;
	}
	
	bubble_input(cont,name,{x:1, w:cont.dx-1, type:"time_gen"});
}
			
			
/// Loads up the compartment data bubble
function comp_data(bu)
{
	start_data_source("Compartment",{cl_drop:{te:select_drop_str}},bu.op.info);
	select_bubble_over();
	
	edit_source.time_radio = {value:"Periodic"};
	edit_source.exact = {value:"Exact"};
	edit_source.comp_acc = "1";
	edit_source.time_gen = "";
}


///////////////////// Transition data

function trans_data_bubble(cont,type)
{
	inter.bubble.check = "trans_checkbox";
	
	let disable = false; if(type == "view") disable = true;

	cont.dx = 17;
	
	bubble_addtitle(cont,"Transition data",{te:poptrans_data_text});
	
	let res = bubble_add_trans_info(cont,disable);
	
	if(res.err){
		bubble_add_minititle(cont,res.msg); 
		add_end_button(cont,"Next");	
		return;
	}
	
	let spec = edit_source.spec;
	let cl_sel = spec.cl_drop.te;
		
	let but_name = "Next"; if(type == "view") but_name = "Done";

	let ac; 
	if(cl_sel != select_drop_str){
		ac = "AddDataTable"; if(type == "view") ac = "DataTableBack";
	}	
	
	add_end_button(cont,but_name,ac);	
}


/// Loads up the transition data bubble
function trans_data(bu)
{
	let p = model.get_p();
	let filt = load_default_filt(p);

	start_data_source("Transition",{time_radio:{value:"All"},time_start:"",time_end:"",cl_drop:{te:select_drop_str},filter:filt},bu.op.info);
	select_bubble_over();
	
	edit_source.time_step = "";
}


///////////////////// Disease diagnostic test

function diagtest_data_bubble(cont,type)
{
	inter.bubble.check = "checkbox";
	
	let cl_pos = get_cl_pos();

	let wmax = get_w_drop(cl_pos);

	cont.dx = 16.5;

	let spec = edit_source.spec;
	let cl_drop = spec.cl_drop;
	let cl_sel = cl_drop.te;
	
	if(type == "view") cl_pos = undefined;
	
	bubble_addtitle(cont,"Diagnostic tests",{te:diag_test_data_text});
	
	bubble_addparagraph(cont,"Select the classification the test is sensitive to:",0,cont.dx);
	cont.y += 0.2;
	bubble_adddropdown_text(cont,"Classification:",5.5,wmax,cl_drop,cl_pos); 
	
	cont.y += 0.2;
	
	if(cl_sel != select_drop_str){
		let pcl = get_p_cl_from_claname(cl_sel);
		let p = pcl.p, cl = pcl.cl;
		
		cont.y += 0.4;
		bubble_addparagraph(cont,"Select the compartments the test is positive to:",0,cont.dx);
	
		bubble_addscrollable(cont,{type:"comp list", p:p, cl:cl, ymax:bubblescroll_dymax}); 
		
		cont.y += 0.3;
		
		let ch_warn = inter.bubble.check_warning;

		if(ch_warn != undefined){
			cont.lay.add_button({te:ch_warn, x:0, y:cont.y-0.2, dx:cont.dx, dy:0.8, type:"Warning"});
			cont.y += 1;
		}
		
		let eqn_on = true;
		if(sim_options()) eqn_on = false;
			
		bubble_addparagraph(cont,"Sensitivity and specificity of the diagnostic test:",0,cont.dx); 
		bubble_double_input(cont,"Sensitivity:",{type:"Se", eqn:eqn_on },
														 "Specificity:",{type:"Sp", eqn:eqn_on });
		
		cont.y += 0.4;
		bubble_addparagraph(cont,"Set the text used to represent test results:",0,cont.dx); 
		bubble_double_input(cont,"Positive text:",{type:"pos_result"},
														 "Negative text:",{type:"neg_result"});
														 
		if(sim_options()) add_timepoint_options(cont);
	}	
	else{ 
		bubble_input(cont,"Sensitivity:",{hidden:true, type:"Se", eqn:true});
		bubble_input(cont,"Specificity:",{hidden:true, type:"Sp", eqn:true});
		bubble_input(cont,"Positive text:",{hidden:true, type:"pos_result"});
		bubble_input(cont,"Negative text:",{hidden:true, type:"neg_result"});
		cont.y += 0.4;
	}

	let but_text = "Next"; if(type == "view") but_text = "Done";
	let ac; 
	if(cl_sel != select_drop_str){
		ac = "AddDataTable"; if(type == "view") ac = "DataTableBack";
	}	
	
	add_end_button(cont,but_text,ac);
}


/// The scrollable box which allows the used to specify compartments
function diagtest_scrollable(lay)
{
	let cx = 0.5, cy = 0;
	let dx = lay.dx;
	let op = lay.op;
	
	let p = op.p, cl = op.cl;
	let claa = model.species[p].cla[cl];
	let cb = edit_source.spec.check_box;
	if(cb.name != claa.name){
		cb.name = claa.name;
		for(let c = 0; c < claa.comp.length; c++) cb.value[c] = { comp_name_store:claa.comp[c].name, check:false};
	}

	cy = add_checkbox_list(lay,claa.comp,"name",cb.value,cy,dx);
	
	return cy;
}
		
	
/// Initialises the diagtest bubble	
function diagtest_data(bu)
{
	let Se_eqn = create_equation("1","Se");
	let Sp_eqn = create_equation("1","Sp");
	
	start_data_source("Diag. Test",{cl_drop:{te:select_drop_str},check_box:{ name:"",value:[]},Se_eqn:Se_eqn,Sp_eqn:Sp_eqn,pos_result:"+",neg_result:"-"},bu.op.info);
	select_bubble_over();
	
	edit_source.time_radio = {value:"Periodic"};
	edit_source.time_gen = "";
}


/////////////////////////// POPULATION DATA

function population_data_bubble(cont,type)
{
	let spec = edit_source.spec;
	
	cont.dx = 20.5;
	
	bubble_addtitle(cont,"Population data",{te:pop_data_text});
	
	let p = model.get_p();
	
	let disable = false; if(type == "view") disable = true;

	bubble_addscrollable(cont,{type:"pop list", disable:disable, p:p, ymax:bubblescroll_dymax, viewtype:type}); 
	
	let ch_warn = inter.bubble.check_warning;

	if(ch_warn != undefined){
		cont.lay.add_button({te:ch_warn, x:0, y:cont.y-0.4, dx:cont.dx, dy:0.8, type:"Warning"});
		cont.y += 0.8;
	}
	cont.y -= 0.2;
		
	bub_observation_error(cont,spec,disable);
	
	if(sim_options()) add_timepoint_options(cont);
	
	let but_name = "Next"; if(type == "view") but_name = "Done";

	let ac; 
	ac = "AddDataTable"; if(type == "view") ac = "DataTableBack";
	
	add_end_button(cont,but_name,ac);	
}


/// Buttons related to observation error
function bub_observation_error(cont,spec,disable)
{	
	let shup = 0.5;
	let shdo = 0.3;
	
	bubble_add_minititle(cont,"Observation error:"); 
	cont.y -=  shup;
	
	bubble_addradio(cont,0,"Poisson","Poisson",spec.obs_error_dist,disable); cont.y -= 1.4;
	bubble_addradio(cont,6,"Normal","Normal",spec.obs_error_dist,disable); cont.y -= 1.4;
	bubble_addradio(cont,11.5,"Negative binomial","Negative binomial",spec.obs_error_dist,disable);
	cont.y += shdo;

	let odv = spec.obs_error_dist.value;

	switch(odv){
	case "Normal":
		bubble_add_minititle(cont,"Normal definition:",1); 
		cont.y -=  shup;
		bubble_addradio(cont,1,"percent","Percentage",spec.obs_error,disable); cont.y -= 1.4;
		bubble_addradio(cont,8.5,"sd","SD",spec.obs_error,disable);
		if(edit_source.info.siminf != "gen"){
			cont.y -= 1.4;
			bubble_addradio(cont,13,"file","SD in file",spec.obs_error,disable);
		}
		cont.y += shdo;
		break;
	
	case "Poisson":
		cont.y += shdo;
		break;
		
	case "Negative binomial":
		bubble_add_minititle(cont,"Negative binomial definition:",1); 
		cont.y -= shup;
		bubble_addradio(cont,1,"p","Value of p",spec.obs_error_p,disable); 
		if(edit_source.info.siminf != "gen"){
			cont.y -= 1.4;
			bubble_addradio(cont,7,"file","p in file",spec.obs_error_p,disable);
		}
		cont.y += shdo;
		break;
	}
	
	let hid = false; if(!(odv == "Normal" && spec.obs_error.value == "percent")) hid = true;		
	bubble_input(cont,"",{hidden:hid, type:"percent", x:2});
	
	hid = false; if(!(odv == "Normal" && spec.obs_error.value == "sd")) hid = true;		
	bubble_input(cont,"",{hidden:hid, type:"sd", x:2});
	
	hid = false; if(!(odv == "Negative binomial" && spec.obs_error_p.value == "p")) hid = true;		
	bubble_input(cont,"",{hidden:hid, type:"p", x:2}); 
}


/// The scrollable box which allows the used to specify compartments that make up population
function population_scrollable(lay)
{
	let mar = 0.2;
	inter.bubble.check = "pop_checkbox";
	
	let spec = edit_source.spec;
	let filt = spec.filter;
	let sp = model.get_sp();
	
	let cy = 0;
	let dx = lay.dx;
	let op = lay.op;

	let disable = lay.op.disable;
	cy = add_comp_filter(cy,op.p,undefined,mar,dx,filt,disable,lay);
	
	return cy;
}


/// Loads up the population data bubble
function population_data(bu)
{
	let p = model.get_p();
	
	let filt = load_default_filt(p);
	
	start_data_source("Population",{filter:filt,obs_error_dist:{value:"Poisson"},obs_error:{value:"percent"},percent:"10",sd:"",obs_error_p:{value:"p"},p:"0.5"},bu.op.info);
	select_bubble_over();
	
	edit_source.time_radio = {value:"Periodic"};
	edit_source.time_gen = "";
}


////////////////////////// POPULATION TRANSITION DATA

function poptrans_data_bubble(cont,type)
{
	inter.bubble.check = "trans_checkbox";
	
	let disable = false; if(type == "view") disable = true;

	cont.dx = 20;
	
	bubble_addtitle(cont,"Aggregated transition data",{te:poptrans_data_text});
	
	let res = bubble_add_trans_info(cont,disable);
	
	if(res.err){
		bubble_add_minititle(cont,res.msg); 
		add_end_button(cont,"Next");	
		return;
	}
	
	let spec = edit_source.spec;
	let cl_sel = spec.cl_drop.te;
	
	if(cl_sel != select_drop_str){
		bub_observation_error(cont,spec,disable);
		if(sim_options()) add_timepoint_options(cont);
	}	
	
	let but_name = "Next"; if(type == "view") but_name = "Done";

	let ac; 
	if(cl_sel != select_drop_str){
		ac = "AddDataTable"; if(type == "view") ac = "DataTableBack";
	}	
	
	add_end_button(cont,but_name,ac);	
}


/// Adds information which specifies transition and filter
function bubble_add_trans_info(cont,disable)
{		
	let cl_pos = get_cl_pos();
	let wmax = get_w_drop(cl_pos);

	let spec = edit_source.spec;
	
	let cl_drop = spec.cl_drop;

	if(cl_pos.length == 1) cl_drop.te = cl_pos[0].te;
	let cl_sel = cl_drop.te;
	
	bubble_adddropdown(cont,5.5,wmax,cl_drop,cl_pos); cont.y -= 1.1;
	bubble_add_minititle(cont,"Classification:");

	if(cl_sel != select_drop_str){
		let p = model.get_p();
		let sp = model.get_sp();
		let cl = find(sp.cla,"name",cl_sel);
	
		if(sp.cla[cl].ntra == 0) return err("No transitions found for this classification")
	
		bubble_addscrollable(cont,{type:"poptrans list", p:p, cl:cl, disable:disable, ymax:bubblescroll_dymax}); 
		
		let ch_warn = inter.bubble.check_warning;

		if(ch_warn != undefined){
			cont.lay.add_button({te:ch_warn, x:0, y:cont.y-0.4, dx:cont.dx, dy:0.8, type:"Warning"});
			cont.y += 0.8;
		}
		
		switch(edit_source.type){
		case "Transition":
			{
				bubble_add_minititle(cont,"Observation time period:"); 
				cont.y -= 0.4;
				bubble_addradio(cont,0,"All","All time",spec.time_radio,disable);
				cont.y -= 1.4;
				bubble_addradio(cont,10,"Spec","Specified",spec.time_radio,disable);
				
				if(spec.time_radio.value == "Spec"){
					bubble_input(cont,"Start:",{x:1, w:7, no_down:true, type:"time_start"});
					bubble_input(cont,"End:",{x:9, w:7, type:"time_end"});
					cont.y -= 0.2;
				}
				else{
					bubble_input(cont,"",{hidden:true, type:"time_start"});
					bubble_input(cont,"",{hidden:true, type:"time_end"});
				}

				cont.y += 0.3;
			}
			break;
			
		case "Pop. Trans.":
			break;
			
		default: error("Wrong source"); break;
		}
	}
	
	return success();
}


/// The scrollable box which allows the used to specify compartments that make up population
function poptrans_scrollable(lay)
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

	lay.add_button({te:"Transitions:", x:mar, y:cy, dx:w, dy:1, type:"Text", si:si, font:fo, col:BLACK, back_col:BUBBLE_COL});
	cy += 0.8;
			
	let spec = edit_source.spec;

	let filt = spec.filter;
	
	let p = lay.op.p;
	let cl = lay.op.cl;
	let disable = lay.op.disable;

	let sp = model.species[p];
	let claa = sp.cla[cl];
	
	let te = spec.cl_drop.te;

	if(filt.te != te){	
		filt.trans_obs_model = {value:"off"};
	
		filt.te = te;
		filt.tra=[];
		for(let i = 0; i < claa.ntra; i++){
			let tra = claa.tra[i];
			filt.tra[i] = { tra_name_store:tra.name, check:false, prob_eqn:create_equation("0","trans_prob",p,cl), i:tra.i, f:tra.f};
		}
	}
	
	if(claa.ntra == 0){
		let fo2 = get_font(si);
		lay.add_button({te:"No transitions in this classification.", x:mar, y:cy, dx:dx, dy:1, type:"Text", si:si, font:fo2, col:BLACK, back_col:BUBBLE_COL});
		cy += 1;
	}
	else{
		lay.add_button({te:"Obs. model:", x:mar+0.5, y:cy+0.2, dx:dx, dy:1, type:"Text", si:si, font:fosi, col:BLACK, back_col:BUBBLE_COL});
		
		lay.add_radio(4.7,cy,"off","Off",filt.trans_obs_model,{back_col:BUBBLE_COL});
		lay.add_radio(8.7,cy,"on","On",filt.trans_obs_model,{back_col:BUBBLE_COL});
		cy += 1.2;
		
		switch(filt.trans_obs_model.value){
		case "on":
			cy = add_prob_list(lay,claa.tra,"name",filt.tra,cy,dx,p,cl);
			break;
			
		case "off":			
			cy += 0.3;
			cy = add_checkbox_list(lay,claa.tra,"name",filt.tra,cy,dx);
			break;
		}
	}

	cy += 0.4;
	
	cy = add_comp_filter(cy,p,cl,mar,w,filt,disable,lay);
	
	return cy;
}


/// Adds comartment filter 
function add_comp_filter(cy,p,cl_not,mar,w,filt,disable,lay)
{
	let sp = model.species[p];
	
	let si = 0.7;
	let fo = get_font(si,"bold");
	let fosi = get_font(si);
	
	for(let cl2 = 0; cl2 < sp.ncla; cl2++){
		if(cl2 != cl_not){
			let claa = sp.cla[cl2];
			let te = claa.name+" filter:";
			
			lay.add_button({te:te, x:mar, y:cy, dx:w, dy:1, type:"Text", si:si, font:fo, col:BLACK, back_col:BUBBLE_COL});
			cy += 1.1;
		
			if(edit_source.type == "Transition"){
				if(filt.cla[cl2].comp.length == 1){
					lay.add_radio(0,cy,"All","All",filt.cla[cl2].radio,{disable:disable, back_col:BUBBLE_COL});
					lay.add_radio(5,cy,"ObsMod","Obs. Model",filt.cla[cl2].radio,{disable:disable, back_col:BUBBLE_COL});
				}
				else{
					lay.add_radio(0,cy,"All","All",filt.cla[cl2].radio,{disable:disable, back_col:BUBBLE_COL});
					lay.add_radio(4.5,cy,"Comp","Comp.",filt.cla[cl2].radio,{disable:disable, back_col:BUBBLE_COL});
					lay.add_radio(10,cy,"ObsMod","Obs. Model",filt.cla[cl2].radio,{disable:disable, back_col:BUBBLE_COL});
				}
			}
			else{
				if(filt.cla[cl2].comp.length == 1){
					lay.add_radio(0,cy,"All","All",filt.cla[cl2].radio,{disable:disable, back_col:BUBBLE_COL});
					lay.add_radio(4,cy,"ObsMod","Obs. Model",filt.cla[cl2].radio,{disable:disable, back_col:BUBBLE_COL});
				}
				else{
					lay.add_radio(0,cy,"All","All",filt.cla[cl2].radio,{disable:disable, back_col:BUBBLE_COL});
					lay.add_radio(3.5,cy,"Comp","Comp.",filt.cla[cl2].radio,{disable:disable, back_col:BUBBLE_COL});
					lay.add_radio(8,cy,"ObsMod","Obs. Model",filt.cla[cl2].radio,{disable:disable, back_col:BUBBLE_COL});
					lay.add_radio(15,cy,"File","Stratify",filt.cla[cl2].radio,{disable:disable, back_col:BUBBLE_COL});
				}
			}
			cy += 1.6;
			
			switch(filt.cla[cl2].radio.value){
			case "Comp":
				cy = add_checkbox_list(lay,claa.comp,"name",filt.cla[cl2].comp,cy,w);
				break;
				
			case "ObsMod":
				cy = add_prob_list(lay,claa.comp,"name",filt.cla[cl2].comp,cy,w,p,cl2);
				break;
			}
		}
	}

	return cy;
}


/// Loads up the population data bubble
function poptrans_data(bu)
{
	let p = model.get_p();
	let filt = load_default_filt(p);

	start_data_source("Pop. Trans.",{cl_drop:{te:select_drop_str},filter:filt,obs_error_dist:{value:"Poisson"},obs_error:{value:"percent"},percent:"10",sd:"",obs_error_p:{value:"p"},p:"0.5"},bu.op.info);
	select_bubble_over();
	
	edit_source.time_radio = {value:"Periodic"};
	edit_source.time_gen = "";
}


/// Plots a list of checkbox items
function add_checkbox_list(lay,list,prop,so,cy,dx)
{
	inter.bubble.all_check = {check:false, all_check_flag:true, so:so};
	let c = 0; while(c < so.length && so[c].check == true) c++;
	if(c == so.length) inter.bubble.all_check.check = true;

	if(list.length == 0){
		cy = lay.add_paragraph("There are no options.",dx,0.5,cy,BLACK,para_si,para_lh);
		return cy;
	}
	
	if(list.length > 1){
		lay.add_checkbox(0.5,cy,"All","All",inter.bubble.all_check,BUBBLE_COL);
		cy += 1.4;
	}
	
	let w = [];
	let wav = 0;
	for(let c = 0; c < list.length; c++){
		let te = list[c][prop];
		w[c] = 2+text_width(te,get_font(si_radio));
		wav += w[c];
	}
	wav /= list.length;

	let mar = 0.5;
	for(let ncol = 5; ncol >= 1; ncol--){
		if(ncol <= list.length){
			let colw=[];
			let colwsum_st=[], colwsum = 0;
			
			for(let c = 0; c < ncol; c++){
				let wmax = 0;
				for(let i = c; i < list.length; i += ncol){
					if(w[i] > wmax) wmax = w[i];
				}
				colwsum_st[c] = colwsum;
				colw[c] = wmax+1;
				colwsum += colw[c];
			}
			if(ncol == 1 || colwsum < dx-mar){
				let pad = 0; if(ncol > 1) pad = ((dx-2*mar)-colwsum)/(ncol-1);
				for(let c = 0; c < list.length; c++){
					let col = c%ncol;
					let row = Math.floor(c/ncol);
					let te = list[c][prop];
					lay.add_checkbox(mar+colwsum_st[col]+pad*col,cy+row*1.2,te,te,so[c],BUBBLE_COL);
				}
			
				return cy + Math.floor((list.length+ncol-1)/ncol)*1.2+0.2;
			}
		}
	}
	
	return cy;
}


/// Plots a list of checkbox items
function add_prob_list(lay,list,prop,so,cy,dx,p,cl)
{
	inter.bubble.all_check = {check:false, all_check_flag:true, so:so};
	
	let te_list=[]
	for(let c = 0; c < list.length; c++){
		te_list[c] = "P("+list[c][prop]+") = "+so[c].prob_eqn.te;
	}
	
	let w = [];
	let wav = 0;
	for(let c = 0; c < list.length; c++){
		w[c] = 2+text_width(te_list[c],get_font(si_radio));
		wav += w[c];
	}
	wav /= list.length;

	let mar = 0.5;
	for(let ncol = 5; ncol >= 1; ncol--){
		if(ncol <= list.length){
			let colw=[];
			let colwsum_st=[], colwsum = 0;
			
			for(let c = 0; c < ncol; c++){
				let wmax = 0;
				for(let i = c; i < list.length; i += ncol){
					if(w[i] > wmax) wmax = w[i];
				}
				colwsum_st[c] = colwsum;
				colw[c] = wmax+1;
				colwsum += colw[c];
			}
			
			if(ncol == 1 || colwsum < dx-mar){
				let dy = 1.6;
				let pad = 0; if(ncol > 1) pad = ((dx-2*mar)-colwsum)/(ncol-1);
				for(let c = 0; c < list.length; c++){
					let col = c%ncol;
					let row = Math.floor(c/ncol);
				
					let ww = (dx-mar)/ncol;
					let fo = get_font(1.1,"","times");
					
					lay.add_button({te:te_list[c], x:mar+colwsum_st[col]+pad*col, y:0.1+cy+row*dy, dx:ww, dy:dy, eqn:so[c].prob_eqn, ac:"ProbEqn", type:"ProbEqn", si:si_radio,  font:fo, col:BLACK});
				}
			
				return cy + Math.floor((list.length+ncol-1)/ncol)*dy+0.2;
			}
		}
	}
	
	return cy;
}


/// Plots a list of radio items
function add_radio_list(lay,list,prop,so,cy,dx)
{
	let w = [];
	let wav = 0;
	for(let c = 0; c < list.length; c++){
		let te = list[c][prop];
		w[c] = 2+text_width(te,get_font(si_radio));
		wav += w[c];
	}
	wav /= list.length;

	let mar = 0.5;
	for(let ncol = 5; ncol >= 1; ncol--){
		if(ncol <= list.length){
			let colw=[];
			let colwsum_st=[], colwsum = 0;
			
			for(let c = 0; c < ncol; c++){
				let wmax = 0;
				for(let i = c; i < list.length; i += ncol){
					if(w[i] > wmax) wmax = w[i];
				}
				colwsum_st[c] = colwsum;
				colw[c] = wmax+1;
				colwsum += colw[c];
			}
			if(ncol == 1 || colwsum < dx-mar){
				let pad = 0; if(ncol > 1) pad = ((dx-2*mar)-colwsum)/(ncol-1);
				for(let c = 0; c < list.length; c++){
					let col = c%ncol;
					let row = Math.floor(c/ncol);
					let te = list[c][prop];
					lay.add_radio(mar+colwsum_st[col]+pad*col,cy+row*1.2,te,te,so,BUBBLE_COL);
				}
			
				return cy + Math.floor((list.length+ncol-1)/ncol)*1.2+0.2;
			}
		}
	}
}


/////////////////////////// SEQUENCE DATA

function sequence_data_bubble(cont,type)
{
	cont.dx = 11;
	
	bubble_addtitle(cont,"Genetic data",{te:seq_data_text});
	
	let eqn_on = true;
	if(sim_options()) eqn_on = false;
			
	let bub = inter.bubble;
	bubble_add_minititle(cont,"Data type:");
	cont.y -= 0.3;
	let spec = edit_source.spec;
	bubble_addradio(cont,0,"matrix","Matrix",spec.type_radio); cont.y -= 1.4;
	bubble_addradio(cont,6,"snp","SNP",spec.type_radio);
	cont.y += 0.3;
	bubble_input(cont,"Mutation rate:",{type:"mut_rate", eqn:eqn_on });
	
	bubble_input(cont,"Sequence variation:",{type:"seq_var", eqn:eqn_on });
	
	if(spec.type_radio.value == "snp"){
		if(edit_source.edit_spec != true){
			bubble_input(cont,"SNP root:",{type:"snp_root" });
		}
	}
	
	if(sim_options()){
		if(spec.type_radio.value == "snp"){
			bubble_input(cont,"Base-pairs:",{type:"num_basep"});
		}
		
		bubble_input(cont,"Fraction observed:",{type:"frac_obs"});
		
		add_timepoint_options(cont);
	}
	
	if(type == "view") add_end_button(cont,"Done","DataTableBack");
	else add_end_button(cont,"Next","AddDataTable");
}


/// Determines if simulation options should appead in bubble
function sim_options()
{
	if(edit_source.info.siminf == "gen" && edit_source.edit_spec != true) return true;
	return false;
}


/// Loads up the genetic data bubble
function sequence_data(bu)
{
	let mut_rate_eqn = create_equation("1","mut_rate");
	let seq_var_eqn = create_equation("1","seq_var");
	
	start_data_source("Genetic",{mut_rate_eqn:mut_rate_eqn, seq_var_eqn:seq_var_eqn,snp_root:"SNP",type_radio:{value:"matrix"}},bu.op.info);
	
	select_bubble_over();
	
	edit_source.time_radio = {value:"Periodic"};
	edit_source.time_gen = "";
	
	edit_source.numbp = 1000;
	edit_source.frac_obs = 1;
}


/////////////////////// IND EFFECT DATA

function ind_eff_data_bubble(cont,type)
{
	let sp = model.get_sp();
	
	let pos = [];
	for(let g = 0; g < sp.ind_eff_group.length; g++){
		let ieg = sp.ind_eff_group[g];
		for(let j = 0; j < ieg.ie_list.length; j++){
			pos.push({te:ieg.ie_list[j].name});
		}
	}
	
	let wmax = 7;

	cont.dx = 4.2+wmax;

	bubble_addtitle(cont,"Individual effect",{te:move_ind_text});
	
	let drop = edit_source.spec.drop;
	
	if(type != "view"){
		if(pos.length == 0){
			bubble_addparagraph(cont,"There are no individual effect in the model.",0,cont.dx); cont.y += 0.2;
		}
		else{
			bubble_addparagraph(cont,"The individual effect for which values are known:",0,cont.dx);cont.y += 0.2;
			
			bubble_adddropdown(cont,0.5,wmax,drop,pos);
			cont.y += 0.2;
			
			let ac; if(drop.te != select_drop_str) ac = "AddDataTable";
		
			add_end_button(cont,"Next",ac);	
		}
	}
	else{
		bubble_addparagraph(cont,"The individual effect for which values are known: <e>"+drop.te+"</e>",0,cont.dx); cont.y += 0.2;
	}
}


/// Loads up the ind effect data bubble
function ind_eff_data(bu)
{
	start_data_source("Ind. Eff.",{p:model.get_p(), drop:{te:select_drop_str}},bu.op.info);
	select_bubble_over();
}


/////////////////////// IND GROUP DATA

function ind_group_data_bubble(cont,type)
{
	let spec = edit_source.spec;
	
	cont.dx = 20;
	bubble_addtitle(cont,"Individual group",{te:move_ind_text});
	
	bubble_input(cont,"Group name",{x:0, w:cont.dx-1, type:"group_name"});
	
	let bub = inter.bubble;
	
	bubble_addradio(cont,0,"Specify","Specify",bub.ind_group_radio); cont.y -= 1.4;
	bubble_addradio(cont,6,"Wildcard","Use wildcard",bub.ind_group_radio);
	if(subtab_name() != "Generate Data"){
		cont.y -= 1.4;
		bubble_addradio(cont,13,"Load","Load",bub.ind_group_radio);
	}
	
	cont.y += 0.2;
	
	switch(bub.ind_group_radio.value){
	case "Specify":
		cont.y += 0.2;
		bubble_add_minititle(cont,"Individuals in group:"); 
		
		bubble_addscrollable(cont,{type:"ind list", p:spec.p, name_list:inter.bubble.name_list, ymax:bubblescroll_dymax}); 
			
		cont.y += 0.4;
		add_end_button(cont,"Done","AddIndGroup");	
		break;
		
	case "Wildcard":
		bubble_input(cont,"Wildcard",{x:0, w:cont.dx-1, type:"wild_card"});
		add_end_button(cont,"Apply","IndWildcard");	
		break;
		
	case "Load":
		add_end_button(cont,"Load","AddDataTable");	
		break;
	}
}


/// Loads up the individual group data bubble
function ind_group_data(bu)
{
	select_bubble_over();
	
	let bub = inter.bubble;
	bub.name_list = [];
	bub.check_box = {value:[]};
	bub.ind_group_radio = {value:"Specify"};
	bub.wildcard = "";
	
	start_data_source("Ind. Group",{p:model.get_p(),gname:""},bu.op.info);
	
	let ty = "generate"; if(tab_name() == "Inference") ty = "data"; 
		
	start_worker("Get Ind List",{p:model.get_p(), ty:ty});	
}


/// The scrollable box which allows the user to specify individuals
function individual_scrollable(lay)
{
	let cx = 0.5, cy = 0;
	let dx = lay.dx;
	let op = lay.op;
	
	let p = op.p;
	let cb = inter.bubble.check_box;

	cy = add_checkbox_list(lay,op.name_list,"name",cb.value,cy,dx);
	
	return cy;
}
		
		
/// Loads up the individual group data bubble
function ind_group_data2(name_list)
{
	let bub = inter.bubble;
	bub.name_list = name_list;
	
	let cb = bub.check_box;
	for(let i = 0; i < name_list.length; i++){
		cb.value[i] = {check:false};
	}
	
	generate_screen();
}


/// Selects individuals using the wildcard
function select_ind_wildcard()
{
	let bub = inter.bubble;
	
	let filt = bub.wildcard;
	
	bub.ind_group_radio.value = "Specify";
	
	let name_list = bub.name_list;
	
	let cb = bub.check_box;
	for(let i = 0; i < name_list.length; i++){
		cb.value[i].check = false;
		if(wildcard_match(name_list[i].name,filt) == true){
			cb.value[i].check = true;
		}
	}
}


/// Determines if a string matchs a wildcard
function wildcard_match(name,filt)
{
	let spl = filt.split("×");

	return wildcard_match2(name,spl);
}


/// In text name looks for spl[0]
function wildcard_match2(name,spl)
{
	let te = spl[0].toLowerCase(), len = te.length;

	if(name.length >= len){
		if(spl.length == 1){
			if(len == 0) return true;
			if(name.toLowerCase() == te) return true;
		}
		else{
			let imax = name.length-len;
			for(let i = 0; i <= imax; i++){
				if(name.substr(i,len).toLowerCase() == te){
					if(spl.length == 1) return true;
					else{
						let name_next = name.substr(i+len);
						let spl_next=[];
						for(let k = 1; k < spl.length; k++) spl_next.push(spl[k]);
						return wildcard_match2(name_next,spl_next);
					}
				}
			}
		}
	}
	
	return false;
}

function add_ind_group_data()
{
	let head = [];
	head.push("ID");
		
	let bub = inter.bubble;

	let cb = bub.check_box.value;
	let name_list = bub.name_list;
	
	let ele = [];
	for(let i = 0; i < name_list.length; i++){
		if(cb[i].check){
			let row = [];
			row.push(name_list[i].name); 
			ele.push(row);
		}
	}
	
	add_gen_data(head,ele);		
	generate_screen();
}

"use strict";
// Gives all the functions relating to the loading of data files

//////////////////////// INIT POP

function initpop_data_bubble(cont,type)
{
	let sp = model.get_sp();
	
	let spec = inter.edit_source.spec;
	
	cont.dx = 7;
	
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
		
		bubble_addparagraph(cont,"Set using:",0,10); cont.y += 0.2;
		bubble_addradio(cont,0,"Graphical","Graphical interface",spec.radio,disable);
		bubble_addradio(cont,0,"File","File",spec.radio,disable);
	
		let ac = "AddInitialPopulation2";
		
		//if(spec.radio.value == "Graphical" && model.species[p].ncla > 1){
		if(sp.ncla > 1){
			cont.y += 0.3;
			bubble_addparagraph(cont,"Define:",0,10); cont.y += 0.2;
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
			
		cont.y += 0.3;

		if(type == "view") add_end_button(cont,"Done","DataTableBack");
		else add_end_button(cont,"Next",ac);
	}		
}

				
function initpop_data(bu)
{
	start_data_source("Init. Pop.",{radio:{value:"Graphical"},radio2:{value:"Focal"}, focal:{te:"Select"}},bu.op.info);
	
	if(bu.op.siminf == "gen"){
		let spec = inter.edit_source.spec;
		spec.radio.value = "File";
		spec.radio2.value = "All";
		set_loadcol();
		sim_data(); 
	}
	else select_bubble_over();
}

function initpop_data2()
{
	let spec = inter.edit_source.spec;
	
	let sp = model.get_sp(); 
	if(sp.ncla == 1){
		spec.radio2.value = "Focal";
		spec.focal.te = sp.cla[0].name;
	}
	
	set_loadcol();
	
	switch(spec.radio.value){
	case "Graphical":
		start_add_initial_population();
		break;

	case "File":
		file_add_datatable();
		break;

	default: error("Option not recognised 24"); break;
	}
	close_bubble();	
}


/// Start editing of initial population
function start_add_initial_population()
{
	let sp = model.get_sp();
	
	let so = inter.edit_source;

	delete so.load_datatable;
	delete so.table_loaded;

	switch(so.spec.radio2.value){
	case "Focal":
		so.cla = [];
	
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			so.cla[cl] = { cl_name_store:claa.name, comp_init_pop:[]};

			let per = 100/claa.ncomp;
			for(let c = 0; c < claa.ncomp; c++){
				let co = claa.comp[c];
				so.cla[cl].comp_init_pop[c] = { comp_name_store:co.name, pop:0, pop_per:per};
			}
		}
		break;
		
	case "All":	
		let p = model.get_p();
		so.glob_comp = model.get_glob_comp(p);
		
		for(let c = 0; c < inter.edit_source.glob_comp.length; c++){
			so.glob_comp[c].pop = 0;
		}
		break;

	default: error("Option not recognised 25"); break;
	}
}


//////////////////////// INIT POP PRIOR

function initpopprior_data_bubble(cont,type)
{
	let spec = inter.edit_source.spec;
	
	cont.dx = 7;
	bubble_addtitle(cont,"Population prior",{title:"Initial population prior", te:init_pop_prior_text});
	
	//let p = model.get_p();
	let sp = model.get_sp();  //species[p];
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
		bubble_addparagraph(cont,"Prior type:",0,10);
		bubble_addradio(cont,0,"Flat","Flat prior",spec.radio);
		bubble_addradio(cont,0,"Dirichlet","Dirichlet prior",spec.radio);

		cont.y += 0.3;

		add_end_button(cont,"Add","AddInitPopPrior2");
	}		
}

				
function initpopprior_data(bu)
{
	start_data_source("Init. Pop. Prior",{radio:{value:"Flat"},radio2:{value:"Focal"}},bu.op.info);
	select_bubble_over();
}

function initpopprior_data2()
{
	let spec = inter.edit_source.spec;
	
	inter.edit_source.variety = spec.radio.value;
	switch(spec.radio.value){
	case "Flat":
		if(data_source_check_error() == false){
			data_source_add_and_show();
		}
		break;

	case "Dirichlet":
		file_add_datatable();
		break;

	default: error("Option not recognised 26"); break;
	}

	close_bubble();	
}


/////////////////////// MOVE INDIVIUDUALS

function move_data_bubble(cont,type)
{
	let sp = model.get_sp();
	let cl_pos;
	if(type != "view") cl_pos = get_pos(sp.cla,"name");

	let wmax = get_w_drop(cl_pos);
	let cl_drop = inter.edit_source.spec.cl_drop;

	cont.dx = 5.2+wmax;

	bubble_addtitle(cont,"Move individuals",{te:move_ind_text});
	bubble_addparagraph(cont,"The classification in which individuals move compartment must be identified.",0,cont.dx); cont.y += 0.2;
	
	bubble_adddropdown(cont,5.2,wmax,cl_drop,cl_pos); cont.y -= 1.2;
	bubble_addparagraph(cont,"Classification:",0,10); cont.y += 0.2;
	cont.y += 0.3;
	
	if(type != "view"){
		let but_name = "Next"; if(type == "view") but_name = "Done";

		let ac; 
		if(cl_drop.te != "Select"){
			ac = "AddDataTable"; if(type == "view") ac = "DataTableBack";
		}	
		
		add_end_button(cont,but_name,ac);	
	}
}

/// Loads up the move data bubble
function move_data(bu)
{
	start_data_source("Move Ind.",{p:model.get_p(), cl_drop:{te:"Select"}},bu.op.info);
	select_bubble_over();
}


//////////////////////// COMP Data

function comp_data_bubble(cont,type)
{
	let cl_pos = get_cl_pos();
	
	let wmax = get_w_drop(cl_pos);
	let cl_drop = inter.edit_source.spec.cl_drop;

	cont.dx = 10; if(cont.dx < wmax+5) cont.dx = wmax+5;

	if(type == "view") cl_pos = undefined;
	
	bubble_addtitle(cont,"Compartmental data",{te:comp_data_text});
	bubble_addparagraph(cont,"Select the classification on which measurements are made.",0,cont.dx+2); 
	cont.y += 0.2;

	bubble_adddropdown_text(cont,"Classification:",5.5,wmax,cl_drop,cl_pos); 
	cont.y += 0.5;
	
	if(inter.edit_source.info.siminf == "gen"){
		add_timepoint_options(cont);
	}
	
	if(type == "view") add_end_button(cont,"Done","DataTableBack");
	else{
		let ac; if(cl_drop.te != "Select") ac = "AddDataTable";
		add_end_button(cont,"Next",ac);	
	}
}


/// Adds options for time points (used to generate simulated data)
function add_timepoint_options(cont)
{
	let so = inter.edit_source;
	bubble_add_minititle(cont,"Times of observations:"); 

	if(cont.dx > 15.5){ // Sideways
		bubble_addradio(cont,0,"Periodic","Periodic",so.time_radio); cont.y -= 1.4;
		bubble_addradio(cont,5.5,"Fixed time","Fixed time",so.time_radio); cont.y -= 1.4;
		bubble_addradio(cont,11,"Specified","Specified",so.time_radio);
	}
	else{
		bubble_addradio(cont,0,"Periodic","Periodic",so.time_radio);
		bubble_addradio(cont,0,"Fixed time","Fixed time",so.time_radio);
		bubble_addradio(cont,0,"Specified","Specified",so.time_radio);
	}
	
	let name;
	switch(so.time_radio.value){
	case "Periodic": name = "Time-step:"; break;
	case "Fixed time": name = "Time:"; break;
	case "Specified": name = "Times (comma separated):"; break;
	}
	
	bubble_input(cont,name,{type:"time_gen"});
}
			
			
/// Loads up the compartment data bubble
function comp_data(bu)
{
	start_data_source("Compartment",{cl_drop:{te:"Select"}},bu.op.info);
	select_bubble_over();
	
	inter.edit_source.time_radio = {value:"Periodic"};
	inter.edit_source.time_gen = "";
}

///////////////////// Transition data

function trans_data_bubble(cont,type,variety)
{
	let pos=[];
	
	let sp = model.get_sp();
	
	for(let cl = 0; cl < sp.ncla; cl++){
		let claa = sp.cla[cl];
		for(let i = 0; i < claa.ntra; i++){
			let tr = claa.tra[i];
			let vari = tr.variety;
			
			let fl = false;
			switch(variety){
			case "TransVar": 
				if(vari == "Normal") pos.push({te:tr.name, cl:cl, i:i}); 
				break;
				
			case "SourceVar": 
				if(vari == "Source") pos.push({te:claa.comp[tr.f].name, cl:cl, i:i}); 
				break;
			
			case "SinkVar": 
				if(vari == "Sink") pos.push({te:claa.comp[tr.i].name, cl:cl, i:i}); 
				break;

			default: error("Option not recognised 26"); break;
			}
		} 
	}
	let npos = pos.length;
	
	let wmax = get_w_drop(pos);
		
	let spec = inter.edit_source.spec;
	let tr_drop = spec.tr_drop;
	let tr_sel = tr_drop.te;

	cont.dx = 5.2+wmax;

	if(type == "view") pos = undefined;
	
	let title, desc, mini, te;
	switch(variety){
	case "TransVar":
		title = "Transition data";
		desc = "Select the transition on which measurements are made."; 
		if(npos == 0) desc = "There are no transitions in the model.";
		mini = "Transition:";
		te = trans_data_text;
		break;
		
	case "SourceVar": 
		title = "Source data";
		desc = "Select the compartment into which individuals enter.";
		if(npos == 0) desc = "There are no sources in the model.";
		mini = "Compartment:";
		te = source_data_text;
		break;
		
	case "SinkVar": 
		title = "Sink data"; 
		desc = "Select the compartment from which individuals leave."; 
		if(npos == 0) desc = "There are no sinks in the model.";
		mini = "Compartmemt:";
		te = sink_data_text;
		break;

	default: error("Option not recognised 27"); break;
	}
	
	bubble_addtitle(cont,title,{te:te});

	bubble_addparagraph(cont,desc,0,cont.dx); cont.y += 0.2;

	if(npos == 0){
		add_end_button(cont,"Done","CloseBubble");	
		return;
	}
	
	bubble_adddropdown_text(cont,mini,5.2,wmax,tr_drop,pos); 
			
	cont.y += 0.5;

	if(tr_sel != "Select"){
		let disable = false; if(type == "view") disable = true;
		bubble_addparagraph(cont,"Observation time period:",0,10); cont.y += 0.2;
		bubble_addradio(cont,0,"All","All time",spec.time_radio,disable);
		bubble_addradio(cont,0,"Spec","Specified time range",spec.time_radio,disable);
		
	
		if(inter.edit_source.info.siminf != "gen"){
			bubble_addradio(cont,0,"File","Time range in file",spec.time_radio,disable);
		}
		
		if(spec.time_radio.value == "Spec"){
			cont.y += 0.3;
			bubble_input(cont,"Time start observation:",{type:"time_start"});
			bubble_input(cont,"Time end observation:",{type:"time_end"});
		}
		else{
			bubble_input(cont,"",{hidden:true, type:"time_start"});
			bubble_input(cont,"",{hidden:true, type:"time_end"});
		}

		cont.y += 0.3;
	}
	
	let but_name = "Next"; if(type == "view") but_name = "Done";

	let ac; 
	if(tr_sel != "Select"){
		ac = "AddDataTable"; if(type == "view") ac = "DataTableBack";
	}	

	add_end_button(cont,but_name,ac);	
}


/// Loads up the transition data bubble
function trans_data(bu,type)
{
	start_data_source(type,{tr_drop:{te:"Select"},time_radio:{value:"All"},time_start:"",time_end:""},bu.op.info);
	select_bubble_over();
}

///////////////////// Disease diagnostic test

function diagtest_data_bubble(cont,type)
{
	inter.bubble.check = "checkbox";
	
	let cl_pos = get_cl_pos();

	let wmax = get_w_drop(cl_pos);

	cont.dx = 16.5;

	let spec = inter.edit_source.spec;
	let cl_drop = spec.cl_drop;
	let cl_sel = cl_drop.te;
	
	if(type == "view") cl_pos = undefined;
	
	bubble_addtitle(cont,"Diagnostic tests",{te:diag_test_data_text});
	
	bubble_addparagraph(cont,"Select the classification the test is sensitive to:",0,cont.dx);
	cont.y += 0.2;
	bubble_adddropdown_text(cont,"Classification:",5.5,wmax,cl_drop,cl_pos); 
	
	cont.y += 0.2;
	
	if(cl_sel != "Select"){
		let pcl = get_p_cl_from_claname(cl_sel);
		let p = pcl.p, cl = pcl.cl;
		
		cont.y += 0.4;
		bubble_addparagraph(cont,"Select the compartments the test is positive to:",0,cont.dx);
	
		bubble_addscrollable(cont,{type:"comp list", p:p, cl:cl, ymax:15}); 
		
		cont.y += 0.3;
		
		let ch_warn = inter.bubble.check_warning;

		if(ch_warn != undefined){
			cont.lay.add_button({te:ch_warn, x:0, y:cont.y-0.2, dx:cont.dx, dy:0.8, type:"Warning"});
			cont.y += 1;
		}
		
		let eqn_on = true;
		if(inter.edit_source.info.siminf == "gen") eqn_on = false;
			
		bubble_addparagraph(cont,"Sensitivity and specificity of the diagnostic test:",0,cont.dx); 
		bubble_double_input(cont,"Sensitivity:",{type:"Se", eqn:eqn_on },
														 "Specificity:",{type:"Sp", eqn:eqn_on });
		
		cont.y += 0.4;
		bubble_addparagraph(cont,"Set the text used to represent test results:",0,cont.dx); 
		bubble_double_input(cont,"Positive text:",{type:"pos_result"},
														 "Negative text:",{type:"neg_result"});
														 
		if(inter.edit_source.info.siminf == "gen"){
			add_timepoint_options(cont);
		}
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
	if(cl_sel != "Select"){
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
	let cb = inter.edit_source.spec.check_box;
	if(cb.name != claa.name){
		cb.name = claa.name;
		for(let c = 0; c < claa.comp.length; c++) cb.value[c] = { comp_name_store:claa.comp[c].name, check:false};
	}

	cy = add_checkbox_list(lay,claa.comp,"name",cb.value,cy,dx);
	
	cy = lay.add_gap(cy,0.5);
	
	return cy;
}
		
	
/// Initialises the diagtest bubble	
function diagtest_data(bu)
{
	let Se_eqn = create_equation("1","Se");
	let Sp_eqn = create_equation("1","Sp");
	
	start_data_source("Diag. Test",{cl_drop:{te:"Select"},check_box:{ name:"",value:[]},Se_eqn:Se_eqn,Sp_eqn:Sp_eqn,pos_result:"+",neg_result:"-"},bu.op.info);
	select_bubble_over();
	
	inter.edit_source.time_radio = {value:"Periodic"};
	inter.edit_source.time_gen = "";
}

/////////////////////////// POPULATION DATA

function population_data_bubble(cont,type)
{
	let spec = inter.edit_source.spec;
	
	cont.dx = 16.5;
	
	bubble_addtitle(cont,"Population data",{te:pop_data_text});
	
	let p = model.get_p();
		
	cont.y += 0.4;
	
	bubble_add_minititle(cont,"Compartments:"); 
	cont.y -= 0.4;
	bubble_addparagraph(cont,"Select compartments which define population:",0,cont.dx);
	cont.y += 0.2;
	bubble_addscrollable(cont,{type:"pop list", p:p, ymax:15, viewtype:type}); 
	
	let ch_warn = inter.bubble.check_warning;

	if(ch_warn != undefined){
		cont.lay.add_button({te:ch_warn, x:0, y:cont.y-0.4, dx:cont.dx, dy:0.8, type:"Warning"});
		cont.y += 0.8;
	}
		
	let disable = false; if(type == "view") disable = true;
		
	bubble_add_minititle(cont,"Observation error:"); 
	cont.y -= 0.4;
	bubble_addparagraph(cont,"Select the observation error :",0,cont.dx);
	cont.y += 0.2;
	bubble_addradio(cont,0,"percent","Normal (specified percentage)",spec.obs_error,disable);
	
	let hid = false; if(spec.obs_error.value != "percent") hid = true;		
	bubble_input(cont,"",{hidden:hid, type:"percent", x:2,w:8});
	
	bubble_addradio(cont,0,"sd","Normal (specified SD)",spec.obs_error,disable);
	
	hid = false; if(spec.obs_error.value != "sd") hid = true;		
	bubble_input(cont,"",{hidden:hid, type:"sd", x:2,w:8});
	
	if(inter.edit_source.info.siminf != "gen"){
		bubble_addradio(cont,0,"file","Normal (SD specified in file)",spec.obs_error,disable);
	}
	
	if(inter.edit_source.info.siminf == "gen"){
		add_timepoint_options(cont);
	}
	
	let but_name = "Next"; if(type == "view") but_name = "Done";

	let ac; 
	ac = "AddDataTable"; if(type == "view") ac = "DataTableBack";
	
	add_end_button(cont,but_name,ac);	
}


// The scrollable box which allows the used to specify compartments that make up population
function population_scrollable(lay)
{
	let mar = 0.2;
	inter.bubble.check = "pop_checkbox";
	
	let spec = inter.edit_source.spec;
	let filt = spec.filter;
	let sp = model.get_sp();
	
	let cy = 0;
	let dx = lay.dx;
	let op = lay.op;

	for(let cl = 0; cl < sp.ncla; cl++){
		let claa = sp.cla[cl];
		
		let si = 0.7;
		let fo = get_font(si,"bold");
		let te = (claa.name+":").toUpperCase();
		let w = dx;
		
		lay.add_button({te:te, x:mar, y:cy, dx:w, dy:1, type:"Text", si:si, font:fo, col:BLACK, back_col:BUBBLE_COL});
		cy += 1.1;
	
		let op = {back_col:BUBBLE_COL}; if(lay.op.viewtype == "view") op.disable = true;
		
		if(inter.edit_source.info.siminf == "gen"){
			lay.add_radio(0,cy,"All","All",filt.cla[cl].radio,op);
			lay.add_radio(6,cy,"Comp","Compartment",filt.cla[cl].radio,op);
		}
		else{
			lay.add_radio(0,cy,"All","All",filt.cla[cl].radio,op);
			lay.add_radio(4,cy,"Comp","Compartment",filt.cla[cl].radio,op);
			lay.add_radio(12,cy,"File","File",filt.cla[cl].radio,op);
		}
		cy += 1.4;
		
		if(filt.cla[cl].radio.value == "Comp"){
			cy = add_checkbox_list(lay,claa.comp,"name",filt.cla[cl].comp,cy,dx);
		}
	}
	
	cy = lay.add_gap(cy,0.8);
		
	return cy;
}


/// Loads up the population data bubble
function population_data(bu)
{
	let filt={};
	let sp = model.get_sp();
	
	filt.cla=[];	
	for(let cl = 0; cl < sp.ncla; cl++){
		let claa = sp.cla[cl];
		filt.cla[cl] = { cl_name_store:claa.name, radio:{value:"All"}, comp:[]};
		for(let c = 0; c < claa.ncomp; c++){
			filt.cla[cl].comp[c] = {check:false, comp_name_store:claa.comp[c].name};
		}
	}

	start_data_source("Population",{filter:filt,obs_error:{value:"percent"},percent:"10",sd:""},bu.op.info);
	select_bubble_over();
	
	inter.edit_source.time_radio = {value:"Periodic"};
	inter.edit_source.time_gen = "";
}

//////////////////////////POPULATION TRANSITION DATA

function poptrans_data_bubble(cont,type)
{
	inter.bubble.check = "trans_checkbox";
	
	let disable = false; if(type == "view") disable = true;
		
	let cl_pos = get_cl_pos();
	let wmax = get_w_drop(cl_pos);

	let spec = inter.edit_source.spec;
	let cl_drop = spec.cl_drop;
	let cl_sel = cl_drop.te;
	
	cont.dx = 16.5;
	
	bubble_addtitle(cont,"Aggregated transition data",{te:poptrans_data_text});
	
	bubble_addparagraph(cont,"Select classification in which the transition lies:",0,cont.dx); cont.y += 0.2;

	bubble_adddropdown(cont,5.2,wmax,cl_drop,cl_pos); cont.y -= 1.2;
	bubble_addparagraph(cont,"Classification:",0,10); cont.y += 0.2;
	cont.y += 0.3;

	if(cl_sel != "Select"){
		let p = model.get_p();
		let sp = model.get_sp();
		let cl = find(sp.cla,"name",cl_sel);
	
		cont.y += 0.4;
		
		if(sp.cla[cl].ntra == 0){
			bubble_add_minititle(cont,"No transitions found for this classification"); 
			add_end_button(cont,"Next");	
			return;
		}
		
		bubble_add_minititle(cont,"Transitions:"); 
		cont.y -= 0.4;
		bubble_addparagraph(cont,"Select observed transition(s) (with optional compartmental filters in other classifications):",0,cont.dx);
		cont.y += 0.2;
		
		bubble_addscrollable(cont,{type:"poptrans list", p:p, cl:cl, disable:disable, ymax:11}); 
		let ch_warn = inter.bubble.check_warning;

		if(ch_warn != undefined){
			cont.lay.add_button({te:ch_warn, x:0, y:cont.y-0.4, dx:cont.dx, dy:0.8, type:"Warning"});
			cont.y += 0.8;
		}
		
		cont.y += 0.2;
		bubble_add_minititle(cont,"Observation error:"); 
		cont.y -= 0.4;
		bubble_addparagraph(cont,"Select the observation error :",0,cont.dx);
		cont.y += 0.2;
	
		bubble_addradio(cont,0,"percent","Normal (specified percentage)",spec.obs_error,disable);
		
		let hid = false; 
		if(spec.obs_error.value != "percent") hid = true;		
		else cont.y += 0.1;
		bubble_input(cont,"",{hidden:hid, type:"percent", x:2,w:8});
		
		bubble_addradio(cont,0,"sd","Normal (specified SD)",spec.obs_error,disable);
		
		hid = false; if(spec.obs_error.value != "sd") hid = true;		
		bubble_input(cont,"",{hidden:hid, type:"sd", x:2,w:8});
		
		if(inter.edit_source.info.siminf != "gen"){
			bubble_addradio(cont,0,"file","Normal (SD specified in file)",spec.obs_error,disable);
		}
	
		if(inter.edit_source.info.siminf == "gen"){
			bubble_input(cont,"Observation time-step",{type:"time_step"});
		}
	}

	let but_name = "Next"; if(type == "view") but_name = "Done";

	let ac; 
	if(cl_sel != "Select"){
		ac = "AddDataTable"; if(type == "view") ac = "DataTableBack";
	}	
	
	add_end_button(cont,but_name,ac);	
}


// The scrollable box which allows the used to specify compartments that make up population
function poptrans_scrollable(lay)
{
	let mar = 0.2;
	let spec = inter.edit_source.spec;

	let filt = spec.filter;
	let p = lay.op.p;
	let cl = lay.op.cl;
	let disable = lay.op.disable;

	let sp = model.species[p];
	let claa = sp.cla[cl];
	
	let te = spec.cl_drop.te;
		
	if(filt.te != te){	
		filt.te = te;
		filt.tra=[];
		for(let i = 0; i < claa.ntra; i++){
			let tra = claa.tra[i];
			filt.tra[i] = { tra_name_store:tra.name, check:false, i:tra.i, f:tra.f};
		}
		
		filt.cla=[];	
		for(let cl2 = 0; cl2 < sp.ncla; cl2++){
			let claa = sp.cla[cl2];
			if(cl2 != cl){
				filt.cla[cl2] = { cl_name_store:claa.name, radio:{value:"All"}, comp:[]};
				for(let c = 0; c < claa.ncomp; c++){
					filt.cla[cl2].comp[c] = {check:false, comp_name_store:claa.comp[c].name};
				}
			}
			else filt.cla[cl2] = {};
		}
	}
	
	let cy = 0;
	let dx = lay.dx;
	let op = lay.op;
	op.back_col = BUBBLE_COL;
	
	let si = 0.7;
	let fo = get_font(si,"bold");

	lay.add_button({te:"TRANSITION:", x:mar, y:cy, dx:dx, dy:1, type:"Text", si:si, font:fo, col:BLACK, back_col:BUBBLE_COL});
	cy += 1.1;
	
	if(claa.ntra == 0){
		let fo2 = get_font(si);
		lay.add_button({te:"No transitions in this classification.", x:mar, y:cy, dx:dx, dy:1, type:"Text", si:si, font:fo2, col:BLACK, back_col:BUBBLE_COL});
		cy += 1;
	}
	else{
		cy = add_checkbox_list(lay,claa.tra,"name",filt.tra,cy,dx);
	}

	cy += 0.4;
	
	for(let cl2 = 0; cl2 < sp.ncla; cl2++){
		if(cl2 != cl){
			let claa = sp.cla[cl2];
			let te = (claa.name+" filter:").toUpperCase();
			let w = dx;
			
			lay.add_button({te:te, x:mar, y:cy, dx:w, dy:1, type:"Text", si:si, font:fo, col:BLACK, back_col:BUBBLE_COL});
			cy += 1.1;
		
			if(inter.edit_source.info.siminf == "gen"){
				lay.add_radio(0,cy,"All","All",filt.cla[cl2].radio,{disable:disable, back_col:BUBBLE_COL});
				lay.add_radio(6,cy,"Comp","Compartment",filt.cla[cl2].radio,{disable:disable, back_col:BUBBLE_COL});
			}
			else{
				lay.add_radio(0,cy,"All","All",filt.cla[cl2].radio,{disable:disable, back_col:BUBBLE_COL});
				lay.add_radio(4,cy,"Comp","Compartment",filt.cla[cl2].radio,{disable:disable, back_col:BUBBLE_COL});
				lay.add_radio(12,cy,"File","File",filt.cla[cl2].radio,{disable:disable, back_col:BUBBLE_COL});
			}
			cy += 1.4;
			
			if(filt.cla[cl2].radio.value == "Comp"){
				cy = add_checkbox_list(lay,claa.comp,"name",filt.cla[cl2].comp,cy,dx);
			}
		}
	}
	
	cy = lay.add_gap(cy,0.8);
	
	return cy;
}


/// Loads up the population data bubble
function poptrans_data(bu)
{
	start_data_source("Pop. Trans.",{cl_drop:{te:"Select"},filter:{},obs_error:{value:"percent"},percent:"10",sd:""},bu.op.info);
	select_bubble_over();
	
	inter.edit_source.time_step = "";
}


/// Plots a list of checkbox items
function add_checkbox_list(lay,list,prop,so,cy,dx)
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
				let pad = 0; if(ncol > 1) pad = ((dx-mar)-colwsum)/(ncol-1);
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


/////////////////////////// SET TRAPS DATA

function settraps_data_bubble(cont,type)
{
	let spec = inter.edit_source.spec;
	
	cont.dx = 16.5;
	
	bubble_addtitle(cont,"Set traps data",{te:traps_data_text});
	
	let p = model.get_p();
		
	cont.y += 0.4;
	
	bubble_add_minititle(cont,"Compartments:"); 
	cont.y -= 0.4;
	bubble_addparagraph(cont,"Select compartments which define population:",0,cont.dx);
	cont.y += 0.2;
	bubble_addscrollable(cont,{type:"pop list", p:p, ymax:15, viewtype:type}); 
	
	let ch_warn = inter.bubble.check_warning;

	if(ch_warn != undefined){
		cont.lay.add_button({te:ch_warn, x:0, y:cont.y-0.4, dx:cont.dx, dy:0.8, type:"Warning"});
		cont.y += 0.8;
	}
		
	bubble_add_minititle(cont,"Trapping probability:"); 
	cont.y -= 0.4;
	bubble_addparagraph(cont,"Define the probability of trapping an individual:",0,cont.dx);
	cont.y += 0.4;
	bubble_input(cont,"",{type:"trap_prob", eqn:true});

	let but_name = "Next"; if(type == "view") but_name = "Done";

	let ac; 
	ac = "AddDataTable"; if(type == "view") ac = "DataTableBack";
	
	add_end_button(cont,but_name,ac);	
}


// The scrollable box which allows the used to specify compartments that make up population
function settraps_scrollable(lay)
{
	let mar = 0.2;
	inter.bubble.check = "pop_checkbox";
	
	let spec = inter.edit_source.spec;
	let filt = spec.filter;
	let sp = model.get_sp();
	
	let cy = 0;
	let dx = lay.dx;
	let op = lay.op;
	op.back_col = BUBBLE_COL;
	
	for(let cl = 0; cl < sp.ncla; cl++){
		let claa = sp.cla[cl];
		
		let si = 0.7;
		let fo = get_font(si,"bold");
		let te = (claa.name+":").toUpperCase();
		let w = dx;
		
		lay.add_button({te:te, x:mar, y:cy, dx:w, dy:1, type:"Text", si:si, font:fo, col:BLACK, back_col:BUBBLE_COL});
		cy += 1.1;
	
		let op = {}; if(lay.op.viewtype == "view") op.disable = true;
		
		lay.add_radio(0,cy,"All","All",filt.cla[cl].radio,op);
		lay.add_radio(4,cy,"Comp","Compartment",filt.cla[cl].radio,op);
		lay.add_radio(12,cy,"File","File",filt.cla[cl].radio,op);
		cy += 1.4;
		
		if(filt.cla[cl].radio.value == "Comp"){
			cy = add_checkbox_list(lay,claa.comp,"name",filt.cla[cl].comp,cy,dx);
		}
	}
	
	cy = lay.add_gap(cy,0.5);
	
	return cy;
}


/// Loads up the population data bubble
function settraps_data(bu)
{
	let filt={};
	let sp = model.get_sp();
	
	filt.cla=[];	
	for(let cl = 0; cl < sp.ncla; cl++){
		let claa = sp.cla[cl];
		filt.cla[cl] = { cl_name_store:claa.name, radio:{value:"All"}, comp:[]};
		for(let c = 0; c < claa.ncomp; c++){
			filt.cla[cl].comp[c] = {check:false, comp_name_store:claa.comp[c].name};
		}
	}
	
	let trap_prob_eqn = create_equation("","trap_prob");
	
	start_data_source("Set Traps",{filter:filt, trap_prob_eqn},bu.op.info);
	select_bubble_over();
}

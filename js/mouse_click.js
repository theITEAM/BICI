"use strict";
// Functions for when a mouse click is performed

/// Fires when the mouse is clicked 
function mouse_click(x,y)                                 
{
	inter.mx = x; inter.my = y;
	
	let l = inter.over.layer;
	let i = inter.over.i;
	
	if(l != undefined){
		let but = inter.layer[l].but[i];
	
		if(inter.loading_symbol.on && but.ac != "Stop") return;
	
		button_action(but,"click");
	}
	else{
		close_dropdown();

		if(inter.bubble.lay_name != undefined){
			close_bubble();
			generate_screen();
		}
	}	
}


/// Actives a button and makes sure in view
function activate_button(lay,i)
{
	inter.over = {layer:lay.index, i:i};
	shift_button_in_view(lay.index,i);	
	button_action(lay.but[i],"click");		
}


/// The action generated when a button is clicked
function button_action(bu,action_type)                                  
{
	let x = Math.round(bu.x), y = Math.round(bu.y);
	let dx = Math.round(bu.dx), dy = Math.round(bu.dy);
	let te = bu.te;
	let val = bu.val, val2 = bu.val2, ac = bu.ac;

	let double_click = false;
	if(inter.mouse_time_down - inter.mouse_time_down_prev < 350){
		let ov1 = inter.mouse_down_over_prev, ov2 = inter.mouse_down_over;
		if(ov1.l == ov2.l && ov1.i == ov2.i) double_click = true;
	}
	
	if(inter.bubble.lay_name != undefined && action_type == "click"){
		let na = inter.layer[inter.over.layer].name;

		switch(na){
		case "Bubble": case "BubbleBack": case "BubbleScrollable": 
		case"Dropdown": case "DropdownOptions":
		case "Expand": case "Equation": 
		case "EquationBackground": case "EquationAddQuantity":
		case "Help": case "HelpBackground": case "Yscroll": case "Xscroll":
		case "Input": case "TextBox":
			break;
		
		default: close_bubble(); break;
		}	
	}
	
	if(ac != "Dropdown") close_dropdown();

	//turn_off_cursor();
	
	pr(ac+" ac");

	switch(ac){
	case "Tab":		
		change_page({pa:val});
		break;
		
	case "PageSub":			
		change_page({su:val});
		break;
		
	case "PageSubSub":	
		let pag = inter.page[inter.pa];
		if(pag.sub[pag.index].index == val){
			if(tab_name() == "Model" && subtab_name() == "Compartments"){
				select_bubble_over();
				let p = model.get_p();
				let sp = model.species[p];
				let bub = inter.bubble;
				bub.radio = {value:sp.type};
				bub.checkbox = {check:sp.trans_tree.check};
				
				bub.infection_cl = sp.infection_cl;		
				let pos=[];
				for(let cl = 0; cl < sp.ncla; cl++) pos.push({te:sp.cla[cl].name});
				bub.infection_cl_pos = pos;		
			}
		}
		else{
			change_page({susu:val});
		}
		break;
		
	case "PageSubSubSub":	
		change_page({sususu:val});
		break;
		
	case "ClassTab":
		if(bu.cl == bu.cl_sel){
			let ind = model.species[bu.p].cla[bu.cl].index;
			let pos = model.get_index_possibility(ind);
			select_bubble_over();
			inter.bubble.drop  = {te:ind, pos:pos}; 
		}
		else change_page({sususu:bu.cl});
		break;
		
	case "HelpIcon":
		inter.help = bu.help;
		break;
		
	case "CloseHelp":
		close_help();
		break;		
	
	case "CloseEquation":
		turn_off_cursor();
		inter.equation={};
		break;
			
	case "ScrollUp":
		change_scroll(-1,bu);
		break;
		
	case "ScrollDown":
		change_scroll(1,bu);
		break;
	
	case "WarningStore":
		model.warn_view = true;
		if(model.update_model_needed == true){
			model.update_model_needed = false;
			start_worker("UpdateModel");
		}
		break;
		
	case "ExpandIcon":  // Turns on the equation editor	
		{
			let l = inter.over.layer; while(l >= 0 && inter.layer[l].name != "Input") l--;
			if(l < 0) error("Cannot find");

			let i = inter.cursor.i;
			
			let lay = inter.layer[l];

			let bs = lay.get_text_box_store();
			
			start_equation(bs.te,bs.eqn,lay.op.source,i);
		}
		break;

	case "PositionCursor":
		break;

	case "EditDescription":
		model.description.edit = true;
		reset_text_box();
		generate_screen();
		set_focus_first();
		break;
		
	case "CancelEditDescription":
		close_description();
		break;
		
	case "DoneEditDescription":
		close_description();
		copy_back_to_source(inter.over.layer);
		break;
		
	case "AddCompartment":
		set_ruler(model.get_p(),model.get_cl());
		inter.mode={ type:"Add_Compartment"};
		
		close_bubble();
		break;
		
	case "ImportCompMap2":
		if(bubble_check_error() == false){
			inter.selected_col = inter.bubble.radio.col;
			loading_dialogue("",".geojson,.json","GEOJSON_compartment"); 
			close_bubble();
		}
		break;
	
		
	case "AddTransition":
		inter.mode={ type:"Add_Transition", i:undefined, f:undefined, midp:[]};
		close_bubble();
		break;
		
	case "ImportTrans2":
		{
			let p = model.get_p();
			let cl = model.get_cl();
			let ty = "Trans File";
			if(inter.bubble.pos_check.check == true) ty += " Pos";
			
			start_data_source(ty,{},{load_file:"Trans", p:p, cl:cl});
			file_add_datatable(); 
			close_bubble();
		}
		break;
		
	case "AddSource":
		inter.mode={ type:"Add_Source", midp:[]};
		close_bubble();
		break;
	
	case "AddSink":
		inter.mode={ type:"Add_Sink", midp:[]};
		close_bubble();
		break;
	
	case "IntermediateAdd":
		{
			let p = model.get_mouse_to_desktop();
			
			switch(inter.mode.type){
			case "Add_Transition":
				if(inter.mode.i != undefined) inter.mode.midp.push(p);
				break;
			
			case "Add_Source":
				inter.mode.midp.push(p);
				break;
				
			case "Add_Sink":
				inter.mode.midp.splice(0,0,p);
				break;
			}
		}
		break;
		
	case "ZoomIn":
		zoom_factor(bu.p,bu.cl,1.3);
		break;
	
	case "ZoomOut":
		zoom_factor(bu.p,bu.cl,1/1.3);
		break;
		
	case "ZoomInTimeline":
		inter.graph.zoom_timeline_factor(1.3);
		break;
	
	case "ZoomOutTimeline":
		inter.graph.zoom_timeline_factor(1.0/1.3);
		break;
		
	case "ZoomInGraph":
		inter.graph.zoom_graph_factor(1.3);
		break;
	
	case "ZoomOutGraph":
		inter.graph.zoom_graph_factor(1.0/1.3);
		break;
		
	case "TimelineGrab":
		if(double_click == true){
			let mx = inter.mx;
			let lay = get_lay("GraphContent");
			inter.graph.zoom_timeline_factor(1.3,(mx-lay.x)/lay.dx);
		}			
		break;
		
	case "GraphGrab":
		if(double_click == true){
			let mx = inter.mx, my = inter.my;
			let lay = get_lay("GraphContent");
			inter.graph.zoom_graph_factor(1.3,(mx-lay.x)/lay.dx,(my-lay.y)/lay.dy);
		}			
		break;
	
	case "Compartment":
		{
			let mo = inter.mode;
			
			let add_trans = {};
			
			switch(mo.type){
			case "Add_Transition":		
				if(mo.i == undefined) mo.i = bu.i;
				else add_trans = { i:mo.i, f:bu.i};
				break;
			
			case "Add_Source":
				add_trans = { i:SOURCE, f:bu.i};	
				break;
				
			case "Add_Sink":
				add_trans = { i:bu.i, f:SINK};	
				break;
				
			default:
				select_bubble_over();
				
				inter.comp_select.p = bu.p;
				inter.comp_select.cl = bu.cl;
				inter.comp_select.list = [bu.i];
				
				let co = model.species[bu.p].cla[bu.cl].comp[bu.i];
				inter.bubble.checkbox = {check:co.markov_branch};
				inter.bubble.checkbox2 = {check:co.all_branches};
				break;
			}
			
			if(add_trans.i != undefined){
				mode_off();
				let res = model.add_transition(bu.p,bu.cl,add_trans.i,add_trans.f,mo.midp,dist_pos[0]);
				output_help(res);
				if(res.err == false){
					model.select_transition_bubble(bu.p,bu.cl,model.species[bu.p].cla[bu.cl].ntra-1);
				}			
			}
		}
		break;
	
	case "Transition":
		select_bubble_over();
		break;
		
	case "TransitionPoint":
		if(bu.variety == MIDPOINT){
			let mo = inter.mode;
			if(mo.type == "Add_Transition" || mo.type == "Add_Source" || mo.type == "Add_Sink"){
				mo.midp.splice(bu.index,1);
			}
			else model.delete_transitionpoint(bu.p,bu.cl,bu.i,bu.index);
		}
		break;
		
	case "CloseBubble":
		let bub = inter.bubble;
		if(bub.bu && bub.bu.type == "PopFilt") done_filter();
		close_bubble();
		break;
	
	case "ColourSelect":
		{
			let sel_bu = bu.sel_bu;
		
			switch(sel_bu.type){
			case "AddButton": 
				inter.bubble.radio.col = bu.col; 
				break;
				
			case "LabelText": case "Box":
				model.species[sel_bu.p].cla[sel_bu.cl].annotation[sel_bu.i].color = bu.col;
				break;
			
			default:
				model.species[sel_bu.p].cla[sel_bu.cl].comp[sel_bu.i].col = bu.col;
				break;
			}
		}
		break;
		
	case "ChangeComp":
		if(bubble_check_error() == false){
			copy_back_to_source();
			
			let te = inter.rename_compartment;
			let bu = inter.bubble.bu;
			let co = model.species[bu.p].cla[bu.cl].comp[bu.i];
		
			
			if(co.choose_branch == true){
				co.markov_branch = inter.bubble.checkbox.check;
				model.determine_branching();
				model.update_pline(bu.p,bu.cl);
			}
			
			if(co.branch == true){
				co.all_branches = inter.bubble.checkbox2.check;
				model.determine_branching();
				model.update_pline(bu.p,bu.cl);
			}
			
			if(co.name == te) close_bubble();
			else{
				start_worker("Rename Compartment",{p:bu.p, cl:bu.cl, c:bu.i, old_name:co.name, new_name:te});	
			}
		}
		break;
		
	case "DeleteComp":
		inter.help = {title: "Delete Compartment", te: "Are you sure you want to delete this compartment?", p:bu.op.p, cl:bu.op.cl, i:bu.op.i, ok:"DeleteCompConfirm"};
		break;

	case "DeleteCompConfirm":
		{
			let he = inter.help;
			start_worker("Delete Comp",{p:he.p, cl:he.cl, c:he.i});
		}
		break;
		
	case "DeleteTrans":
		inter.help = {title: "Delete Transition", te: "Are you sure you want to delete this transition?", p:bu.op.p, cl:bu.op.cl, i:bu.op.i, ok:"DeleteTransConfirm"};
		break;
	
	case "DeleteTransConfirm":
		{
			let he = inter.help;
			start_worker("Delete Trans",{p:he.p, cl:he.cl, i:he.i});
		}
		break;
		
	case "ClassificationBack":
		if(double_click == true) zoom_double_click(bu.p,bu.cl);
		break;
		
	case "MessageClose":
		mode_off();
		break;
	
	case "NothingMode":
		{
			let mo = inter.mode;
			switch(mo.type){
			case "Add_Compartment":
				if(mo.p != undefined){
					let name = model.get_new_comp_name();
					
					let ty;
					if(mo.coord == "latlng"){
						ty = "CompLatLng";
						let res = model.add_latlng_compartment(name,mo.p,mo.cl,mo.mp.x,mo.mp.y,BLACK,false);
						output_help(res);
						if(res.err == false){
							set_ruler(mo.p,mo.cl);
						}
					}
					else{
						ty = "Compartment";
						let res = model.add_compartment(name,mo.p,mo.cl,mo.mp.x,mo.mp.y,WHITE,false,false,true);
						output_help(res);
					}
					
					mode_off();
					model.select_button_bubble("Compartment",ty,mo.p,mo.cl,model.species[mo.p].cla[mo.cl].ncomp-1);
				}
				break;
			
			case "Add_Label":
				model.add_label(mo.te,mo.size,model.get_p(),model.get_cl());
				mode_off();
				break;
			}
		}
		break;
		
	case "SliceTime":
		select_bubble_over();
		inter.bubble.slice_time = bu.rpf.slice_time;
		break;
	
	case "DoneSliceTime":
		if(bubble_check_error() == false){	
			copy_back_to_source();	
			inter.bubble.bu.rpf.slice_time = inter.bubble.slice_time;
			close_bubble();
			inter.graph.init = undefined;
		}
		break;
		
	case "CreateSIRE":
		start_worker("Create SIRE");	
		break;
		
	case "AddSpecies":
		select_bubble_over();
		inter.bubble.name = "";
		inter.bubble.radio = {value:"Population"};
		inter.bubble.checkbox = {check:false};
		break;
		
	case "AddSpeciesOK":
		if(bubble_check_error() == false){	
			copy_back_to_source();		
			let res = model.add_species(inter.bubble.species_name,inter.bubble.radio.value,inter.bubble.checkbox.check);
			output_help(res);
			if(res.err == false){
				initialise_pages();
				change_page({susu:model.species.length-1});
			};
		}
		break;
	
	case "AddClassification":
		{
			let pos = model.get_index_possibility();
			if(pos.length == 0) alertp("There are no more indexes available.");
			else{
				select_bubble_over();
		
				inter.bubble.drop  = {te:pos[0].te, pos:pos}; 
				inter.bubble.radio = {value:"cartesian"};
				inter.bubble.checkbox = {check:true};
			}
		}
		break;
		
	case "AddClassificationOK":
		if(bubble_check_error() == false){
			copy_back_to_source();		
			
			let pag = inter.page[inter.pa];
			let p = pag.sub[pag.index].index;
		
			let bub = inter.bubble;
			
			start_worker("Add Cla",{p:p, name:bub.classification_name, index:bub.drop.te, op:{coord:bub.radio.value, default_map:bub.checkbox.check}});
		}
		break;
		
	case "RadioButton":
		if(inter.bubble.check_radio_press == false){
			inter.bubble.check_radio_press = true;
		}
		
		bu.source.value = bu.value;
		if(bu.source.noupdate != true) inter.graph.init = undefined;
		break;
		
	case "DeleteSpecies":
		inter.help = { title:"Delete", te:"Are you sure you want to delete this species?", ok:"DeleteSpecies2"};
		break;
	
	case "DeleteSpecies2":
		{
			let p = model.get_p();
			start_worker("Delete Species",{p:p});
		}
		break;

	case "DeleteClass":
		inter.help = { title:"Delete", te:"Are you sure you want to delete this classification?", ok:"DeleteClass2"};
		break;

	case "DeleteClass2":
		{
			let p = model.get_p();
			let cl = model.get_cl();
			start_worker("Delete Cla",{p:p, cl:cl});
		}
		break;

	case "ChangeClassName": 
		if(bubble_check_error() == false){
			copy_back_to_source();
			let index_new = inter.bubble.drop.te;
			
			let bu = inter.bubble.bu;
			let claa = model.species[bu.p].cla[bu.cl];
				
			let te = inter.rename_classification;
			if(claa.name != te){
				start_worker("Rename Classification",{p:bu.p, cl:bu.cl, old_name:claa.name, new_name:te});	
			}
			
			if(index_new != claa.index){
				start_worker("Rename Index",{p:bu.p, cl:bu.cl, index_old:claa.index, index_new:index_new});	
			}
			else close_bubble();
		}
		break;
		
	case "ChangeSpeciesName":
		if(bubble_check_error() == false){
			copy_back_to_source();
			
			let bu = inter.bubble.bu;
		
			let sp = model.species[bu.p];
			sp.type = inter.bubble.radio.value;
			sp.trans_tree.check = inter.bubble.checkbox.check;
	
			update_param();
	
			if(inter.rename_species != sp.name){
				start_worker("Rename Species",{p:bu.p, species_old:sp.name, species_new:inter.rename_species});	
			}
			else{
				close_bubble();
			}
		}
		break;
	
	case "Dropdown":
		if(inter.dropdown.source == bu.source) inter.dropdown = {};
		else{
			if(inter.bubble && inter.bubble.lay_name == "Frame") close_bubble();
			inter.dropdown = {source:bu.source, style:bu.style};
		}
		break;
		
	case "DropdownOption":
		{
			let source = bu.source
			let pos = bu.pos;
			source.te = pos.te;
		
			if(pos.check_filter == true) change_remove_filter(pos);
		
			if(pos.i != undefined) source.i = pos.i;
			if(pos.p != undefined) source.p = pos.p;
			if(pos.cl != undefined) source.cl = pos.cl;
			if(pos.tr != undefined) source.tr = pos.cl;
			if(pos.c != undefined) source.c = pos.c;
			if(pos.param != undefined) source.param = pos.param;
			if(pos.list != undefined) source.list = pos.list;
			if(pos.radio != undefined) source.radio = pos.radio;
			if(pos.list_split != undefined) source.list_split = pos.list_split;
			if(pos.radio_split != undefined) source.radio_split = pos.radio_split;
			if(pos.radioy != undefined) source.radioy = pos.radioy;
			if(pos.title != undefined) source.title = pos.title;
			if(pos.type != undefined) source.type = pos.type;
			if(pos.name_list != undefined) source.name_list = pos.name_list;
			if(pos.strat_filt != undefined) source.strat_filt = pos.strat_filt;
			if(pos.createviewflag != undefined) source.createviewflag = pos.createviewflag;
			else source.createviewflag = undefined;
			
			if(source == inter.data_type) close_bubble();			
			inter.dropdown = {};
				
			if(pos.trans_dist){	
				model.determine_branching();
				inter.bubble.find_focus = true;			
			}
		
			if(pos.noupdate != true) inter.graph.init = undefined;
			
			if(bu.source.bub_off) close_bubble();
		}
		break;
		
	case "ChangeTransValue":
		if(bubble_check_error() == false){	
			model.species[bu.op.p].cla[bu.op.cl].tra[bu.op.i].type = inter.bubble.trans_type.te;
			update_param();
			copy_back_to_source();
			model.determine_branching();
			model.update_pline(bu.op.p,bu.op.cl);    
			close_bubble();
		}
		break;
		
	case "Calculator_But":
		{
			if(inter.cursor.i != undefined){
				let paste = bu.te; 
				let shift = 0;
				if(find_in(functi,paste) != undefined){
					if(paste == "pow" || paste == "max" || paste == "min" || paste == "thresh" || paste == "ubound"){
						paste += "(|)";
						shift = 2;
					}
					else{
						paste += "()";
						shift = 1;
					}
				}
				
				if(paste == "Σ"){ paste += "_()"; shift = 2;}
		
				cursor_paste(paste);
				
				inter.cursor.i -= shift;
			}
		}
		break;
	
	case "Letter_But":
		inter.equation.letter = bu.te;
		break;
	
	case "AddToolbar":
		initialise_toolbar(bu.te2);
		break;

	case "CloseToolbar":
		delete inter.equation.toolbar;
		break;
		
	case "CheckboxButton":
		if(inter.bubble.check_radio_press == false){
			inter.bubble.check_radio_press = true;
		}
	
		let bso = bu.source;
		
		if(bso.check == false){
			bso.check = true;
		}
		else bso.check = false;
	
		if(bso.all_check_flag == true){
			let so = bso.so;
			for(let c = 0; c < so.length; c++) so[c].check = bso.check;
		}
		
		if(bso.update == true) inter.graph.init = undefined;
		
		if(bso.update_fraction && bso.check == true) update_popfilt_fraction(bso.p,bso.cl);
	
		if(bso.all_check){
			for(let j = 0; j < bso.list.length; j++){
				bso.list[j].checkb.check = bso.check;
			}				
			inter.graph.init = undefined;
		}
		break;
		
	case "AddParameterDone":
		equation_add_parameter();
		break;
		
	case "AddTensorDone":
		equation_add_tensor();
		break;
		
	case "AddSumDone":
		equation_add_sum();
		break;
		
	case "AddPopDone":
		equation_add_population();
		break;
		
	case "AddIndEffectDone":
		equation_add_indeffect();
		break;
		
	case "AddFixEffectDone":
		equation_add_fixeffect();
		break;
		
	case "AddModParamEqn":
		equation_add_modelparam(bu.i);
		break;
	
	case "AddTimeDone":
		equation_add_time(bu.i);
		break;
		
	case "AddDistanceDone":
		equation_add_distance(bu.i);
		break;
		
	case "Menu":
		select_bubble_over();
		break;	
		
	case "Load":	
		close_bubble();
		loading_dialogue("",".bici","BICI_file");
		break;
		
	case "Save":
		close_bubble();
		inter.file_store.type = "Save";
		start_worker("Save BICI",{ save_type:"save", map_store:map_store});
		break;
		
	case "SaveAs":
		close_bubble();
		inter.file_store.type = "SaveAs";
		start_worker("Save BICI",{save_type:"save", map_store:map_store});
		break;
		
	case "SaveSim":
		close_bubble();
		start_worker("Save Sim");
		break;
	
	case "SaveInf":
		close_bubble();
		start_worker("Save Inf");
		break;
		
	case "SavePPC":
		close_bubble();
		start_worker("Save PPC");
		break;
	
	case "Import":
		change_bubble_mode("Import Menu");
		break;
		
	case "ImportScript":
		close_bubble();
		loading_dialogue("",".bici","Import file");
		break;
		
	case "NewModel":
		start_worker("New model");
		break;
		
	case "EquationDone":
		equation_done();
		break;

	case "CancelEditSource":
		close_bubble();
		close_data_source();
		break;
	
	case "EditCompPop": case "EditCompPopPercent": case "EditGlobCompPop":
		select_bubble_over();
		break;

	case "NotEditCompPopPercent":
		inter.help = { title:"Compartmental percentage", te:"This cannot be editted because it is set by requiring that all percentages add up to 100%"};
		break;

	case "Done": 
		if(bubble_check_error() == false){
			copy_back_to_source();
			close_bubble();
		}
		break;
		
	case "DoneReparamEquation": 
		if(bubble_check_error() == false){
			copy_back_to_source();
			close_bubble();
			update_param();
		}
		break;
		
	case "DeleteSource":	
		if(bu.error == true){  // If error is source then delete immediately
			start_worker("Data Source",{type:"Delete", edit_source:{}, info:{siminf:bu.siminf,p:bu.p,i:bu.i}});
		}
		else{	
			inter.help = { title:"Delete", te:"Are you sure you want to delete this?", ok:"DeleteSource2",  siminf:bu.siminf, p:bu.p, i:bu.i};
		}
		break;
	
	case "DeleteSource2":
		{
			let he = inter.help;
			start_worker("Data Source",{type:"Delete", edit_source:{}, info:{siminf:he.siminf,p:he.p,i:he.i}});
			close_help();
		}
		break;
		
	case "DeleteSource3":
		bu.source.splice(bu.i,1);
		break;
	
	case "AddDataSource":
		start_worker("Data Source",{type:"Add", edit_source:edit_source});
		break;
		
	case "InitPopAdd":
		start_worker("Init pop convert",edit_source);
		break;
		
	case "EditSource":
		scroll_to_top();
	
		start_worker("Data Source",{type:"Edit", info:bu.info});
		break;
		
	case "ViewSource":
		{
			let info = bu.info;
			let p = info.p; 
			if(p == undefined){ error("Data source does not have p"); return}

			let sp = model.species[p];
			switch(bu.info.siminf){
			case "sim": edit_source = sp.sim_source[info.i]; break;
			case "inf": edit_source = sp.inf_source[info.i]; break;
			case "ppc": edit_source = sp.ppc_source[info.i]; break;
			case "gen": edit_source = model.sim_res.plot_filter.species[p].gen_source[info.i]; break;
			}
			edit_source.info_store = info;
		}
		edit_source.edit_spec = true;
		select_bubble_over();
		break;
		
	case "PlaceSourceBack":
		start_worker("Data Source",{type:"Replace", edit_source:edit_source, info:edit_source.info_store});
		break;
		
	case "UploadDataTable":
		{
			if(edit_source.type == "LoadAmatrix" && inter.radio_heading.value == "No"){
				inter.help = {title:"Problem loading", te:"The A matrix must be loaded with a heading line which gives the names of individuals."};
			}
			else{
				let pos;
				switch(inter.radio_format.value){
				case "csv": pos = ".csv"; break;
				case "tsv": pos = ".tsv,.txt"; break;
				default: error("Format wrong"); break;
				}
				
				let file_type = "Data file";
				if(edit_source.type == "LoadAmatrix") file_type = "A matrix";
				
				loading_dialogue("",pos,file_type);
			}
		}
		break;
		
	case "UsePreloadedTable":
		edit_source.data_table_use = bu.i;
		edit_source.table_loaded = true;
		close_help();
		break;
	
	case "SelectColumn":
		transfer_column(bu.c);
		break;
		
	case "BackEditSource":
		transfer_column_back();
		break;
		
	case "AddParamMult":
		select_bubble_over();
		break;
		
	case "AddParamMultParam":
		{	
			let par = model.param[bu.i];
			
			model.param_factor.push({f_param:get_param_mult(par), param:{name:par.name,dep:copy(par.dep), full_name:par.full_name}});
			close_bubble();
			start_worker("UpdateModel");
		}
		break;
		
	case "AddPop":
		start_data_source("Add Pop.",{},bu.op.info);
		file_add_datatable();
		break;
		
	case "RemPop":
		start_data_source("Remove Pop.",{},bu.op.info);
		file_add_datatable();
		break;
		
	case "AddIndividuals":
		start_data_source("Add Ind.",{},bu.op.info);
		file_add_datatable();
		break;

	case "RemIndividuals":
		start_data_source("Remove Ind.",{},bu.op.info);
		file_add_datatable();
		break;

	case "MoveIndividuals": 
		move_data(bu); 
		break;
	
	case "FixedEffect": 
		fixed_effect_data(bu); 
		break;
	
	case "AddInitialPopulation":
		initpop_data(bu); 
		break;
	
	case "AddInitialPopulation2": 
		initpop_data2(); 
		break;

	case "EditTableElement": case "EditTableHead":
		select_bubble_over();
		break;
	
	case "EditParamElement": 
		select_bubble_over();
		break;
	
	case "EditSplineKnots":
		{
			select_bubble_over();
			inter.bubble.knot_radio = {value:"knot_time"};
			
			let t_start = "", t_end = "";
			if(model.sim_details.t_start != "" && model.sim_details.t_end != ""){
				t_start = model.sim_details.t_start;
				t_end = model.sim_details.t_end;
			}
			else{
				if(model.inf_details.t_start != "" && model.inf_details.t_end != ""){
					t_start = model.inf_details.t_start;
					t_end = model.inf_details.t_end;
				}
			}
			
			inter.bubble.knot_times = stringify(model.param[bu.i].spline.knot);
			inter.bubble.knot_t_start = t_start;
			inter.bubble.knot_t_end = t_end;
			inter.bubble.knot_dt = "";
		}
		break;
		
	case "DoneKnots":
		if(bubble_check_error() == false){
			copy_back_to_source();
				
			let bub = inter.bubble;
			if(bub.knot_radio.value == "time_step"){
				let t_start = Number(bub.knot_t_start);
				let t_end = Number(bub.knot_t_end);
				let dt = Number(bub.knot_dt);
				if(t_start >= t_end) alertp("The start time must be before the end time");
				else{
					let vec = [];
					for(let t = t_start; t < t_end-TINY; t+= dt){
						if(t == t_start) vec.push("start");
						else vec.push(String(t));
					}						
					vec.push("end");
				
					change_spline_knot(model.param[bub.bu.i],vec);
					update_param();					
					close_bubble();
				}
			}
			else{
				change_spline_knot(model.param[bub.bu.i],bub.knot_times.split(",")); 
				update_param();
				close_bubble();
			}
		}
		break;
	
	case "UploadKnots":
		start_data_source("KnotTimes",{},{p:"not needed", i:bu.op.i});
		
		file_add_datatable(); 
		break;
		
	case "AddKnotTimes":
		{
			let tab = edit_source.table;
		
			let vec = [];
			for(let r = 0; r < tab.nrow; r++) vec.push(tab.ele[r][0]);
			
			change_spline_knot(model.param[edit_source.info.i],vec); 
			update_param();
			close_data_source();
		}
		break;
	
	case "Search":
		inter.bubble.find = "";
		change_bubble_mode("Search");
		break;
		
	case "DoneSearch":
		if(bubble_check_error() == false){
			copy_back_to_source();
			done_search();
		}
		break;
	
	case "SearchBack": search_back(); break;
	
	case "SearchNext": search_next(); break;
	
	case "Replace": 
		inter.bubble.find = "";
		inter.bubble.replace = "";
		change_bubble_mode("Replace");
		break;
		
	case "DoneReplace":
		if(bubble_check_error() == false){
			copy_back_to_source();
			done_replace();
		}
		break;

	case "Order":
		inter.bubble.order_radio = {value:"A-Z"};
		change_bubble_mode("Order");
		break;

	case "DoneOrder":
		done_order();
		break;

	case "DeleteRows":
		inter.bubble.find = "";
		change_bubble_mode("Delete");
		break;
		
	case "TableEditOn":
		edit_source.table.edit = true;
		break;
		
	case "DoneDeleteRows":
		if(bubble_check_error() == false){
			copy_back_to_source();
			done_delete_rows();
		}
		break;

	case "RowNumber":
		select_bubble_over();
		break;

	case "DeleteRow":
		done_delete_row();
		break;

	case "AddInitPopPrior": 
		initpopprior_data(bu);
		break;
		
	case "AddInitPopPrior2":
		initpopprior_data2(); 
		break;

	case "AddDataTable": 
		if(bubble_check_error() == false){
			scroll_to_top();
			copy_back_to_source();
			set_loadcol();
			file_add_datatable(); 
		}
		break;
	
	case "AddIndGroup":
		if(bubble_check_error() == false){
			copy_back_to_source();
			add_ind_group_data();
		}
		break;
		
	case "DataTableBack":
		if(bubble_check_error() == false){
			copy_back_to_source();
			close_bubble();	
			start_worker("Data Source",{type:"ReplaceView", edit_source:edit_source, info:edit_source.info_store});
		}
		break;
		
	case "DataViewBack":
		if(bubble_check_error() == false){
			copy_back_to_source();
			close_bubble();	
			start_worker("Data Source",{type:"Replace", edit_source:edit_source, info:edit_source.info_store});
		}
		break;

	case "IndWildcard": 
		if(bubble_check_error() == false){
			copy_back_to_source();
			select_ind_wildcard(); 
		}
		break;

	case "CompData": comp_data(bu);	break;
	case "CompView": comp_data_view(bu); break;
	case "TransData": trans_data(bu,"Transition"); break;	
	case "SourceData": trans_data(bu,"Source"); break;	
	case "SinkData": trans_data(bu,"Sink"); break;	
	case "DiagTestData": diagtest_data(bu); break;
	case "SeqData": sequence_data(bu); break;
	case "IndEffData": ind_eff_data(bu); break;
	case "IndGroupData": ind_group_data(bu); break;
	case "PopulationData": population_data(bu); break;
	case "PopTransData": poptrans_data(bu); break;
	case "TableSlider": case "TableSlider2": break;
	case "CloneClassInit": clone_class_init(bu.op.p);	break;
	
	case "CloneClass":
		{
			let drop = inter.bubble.drop;
			model.clone_class(drop.p_sel,drop.p,drop.cl);			
			change_page({sususu:model.species[model.get_p()].ncla-1});
			close_bubble();
		}
		break;

	case "ViewModelError":
		view_warning(bu.i);			
		break;
		
	case "AddAnnotation":
		select_bubble_over();
		inter.bubble.radio = {value:"Label"};
		break;
		
	case "ImportModel":
		select_bubble_over();
		inter.bubble.radio = {value:"Compartments"};
		break;
		
	case "ImportModelNext":
		switch(inter.bubble.radio.value){
		case "Compartments":
			inter.bubble.radio = {value:"SelectColour", col:WHITE};
			inter.bubble.pos_check = {check:true};
			change_bubble_mode("LoadComps");
			break;
			
		case "Transitions":
			inter.bubble.pos_check = {check:false};
			change_bubble_mode("LoadTrans");
			break;
			
		case "Comp. Map":
			inter.p_cl_store = {p:model.get_p(), cl:model.get_cl()};
			inter.bubble.radio = {value:"SelectColour", col:WHITE};
			change_bubble_mode("LoadCompMaps");
			break;
		}
		break;
		
	case "AddAnnotationOK":
		inter.bubble.radio.col = annotation_col_default;
	
		switch(inter.bubble.radio.value){
		case "Label": 
			change_bubble_mode("add_label");
			inter.bubble.label = { te:"", size:String(size_annotation_default)};
			break;
			
		case "Box":
			{
				let op = inter.bubble.bu.op;
				if(model.species[op.p].cla[op.cl].ncomp == 0){
					inter.help = { title:"Bounding box", te:"A bounding box can only be added if there are compartments"};
				}
				else{
					inter.mode={ type:"Add_Box", p:op.p, cl:op.cl};
					close_bubble();
				}
			}
			break;
			
		case "Map":
			inter.p_cl_store = {p:model.get_p(), cl:model.get_cl()};
		
			loading_dialogue("",".geojson,.json","GEOJSON_annotation"); 
			close_bubble();
			break;
		
		case "Remove":
			remove_annotation_init();
			break;
		}
		break;
		
	case "AddLabelOK":
		if(bubble_check_error() == false){
			copy_back_to_source();
		
			inter.mode={ type:"Add_Label", te:inter.bubble.label.te, size:inter.bubble.label.size, col:inter.bubble.radio.col};
			close_bubble();
		}
		break;
		
	case "DeriveOK":
		if(bubble_check_error() == false){
			copy_back_to_source();
			
			let res = model.add_modify_derived(inter.bubble.derived,"modify");
			if(res.err == false) close_bubble();
			else output_help(res);
		}
		break;
		
	case "DeleteAnnotation":
		delete_annotation();
		close_bubble();	
		break;
		
	case "LabelText":
		select_bubble_over();
		inter.bubble.label = {te:bu.te, size:String(size_annotation_default)};
		break;
	
	case "LabelOK":
		if(bubble_check_error() == false){
			copy_back_to_source();
			close_bubble();
		}
		break;
		
	case "DeleteAnno":
		model.species[bu.op.p].cla[bu.op.cl].annotation.splice(bu.op.i,1);
		close_bubble();
		break;
	
	case "Select": break;
	
	case "LoadCompartmentNext":
		if(bubble_check_error() == false){
			let bub = inter.bubble;
			let p = model.get_p();
			let cl = model.get_cl();
			let ty = "Comp File";
			if(bub.pos_check.check == true) ty += " Pos";
			if(bub.radio.value == "FileColour") ty += " Colour";
			
			start_data_source(ty,{},{load_file:"Comp", p:p, cl:cl, col:bub.radio.col});
		
			file_add_datatable(); 
			close_bubble();
		}	
		break;
	
	case "AddFileComp":
		start_worker("Add file comp",edit_source);
		break;
		
	case "AddFileTrans":
		start_worker("Add file trans",edit_source);
		break;
		
	case "Box":
		select_bubble_over();
		break;
		
	case "CloseBubble":
		close_bubble();
		break;	
		
	case "AddCompMap":
		edit_source.col = inter.selected_col;
		start_worker("Add comp map",edit_source);
		break;
	
	case "CancelEditParam":
		close_bubble();
		close_param_source();
		break;
	
	case "CancelAmatrix":
		close_bubble();
		close_Amatrix_source();
		break;
		
	case "CancelXvector":
		close_bubble();
		close_Xvector_source();
		break;
		
	case "EditParamDone":
		if(edit_param_check_error() == false){
			let ep = inter.edit_param;
			
			start_worker("Set Param",{i:ep.i, value:ep.value});
			
			close_bubble();
			close_param_source();
		}
		break;
		
	case "EditPriorSplitDone": case "EditDistSplitDone": 
		if(edit_param_check_error() == false){
			let ep = inter.edit_param;
			start_worker("Set PriorSplit",{i:ep.i, prior_split:ep.prior_split});	
		}
		break;
		
	case "EditReparamDone":
		if(edit_param_check_error() == false){
			let ep = inter.edit_param;
			start_worker("Set Reparam",{i:ep.i, value:ep.value});
			close_bubble();
			close_param_source();
		}
		break;
		
	case "SmoothValueDone":
		if(bubble_check_error() == false){
			copy_back_to_source();
			close_bubble();
		}
		break;
		
	case "LoadTensor": 
		start_data_source(ac,{},{load_file:"LoadTensor"});
		file_add_datatable(); 	
		break;
		
	case "LoadReparam": 
		start_data_source(ac,{},{load_file:"LoadReparam"});
		file_add_datatable(); 	
		break;
	
	case "LoadPriorSplit":
		start_data_source(ac,{},{load_file:"LoadPriorSplit"});
		file_add_datatable(); 	
		break;
		
	case "LoadDistSplit":
		start_data_source(ac,{},{load_file:"LoadDistSplit"});
		file_add_datatable(); 	
		break;
		
	case "LoadTensorDone": 
		start_worker("Load Tensor",{ep:inter.edit_param, source:edit_source});	
		break;
	
	case "LoadReparamDone":
		start_worker("Load Reparam",{ep:inter.edit_param, source:edit_source});	
		break;
		
	case "LoadPriorSplitDone": 
		start_worker("Load PriorSplit",{ep:inter.edit_param, source:edit_source});	
		break;
		
	case "LoadDistSplitDone":
		start_worker("Load PriorSplit",{ep:inter.edit_param, source:edit_source, dist:true});	
		break;
		
	case "SetConstant": case "SetReparam":
		select_bubble_over();
		inter.bubble.set_reparam_type = false;
		break;
		
	case "AddConstParam":
		{
			model.param[bu.i].variety = "const";
			close_bubble();
		}
		break;
		
	case "AddReparamParam":
		{
			let par = model.param[bu.i];
			if(par.dep.length > 0){
				inter.bubble.set_reparam_type = true;
				inter.bubble.k = bu.i;
				inter.bubble.radio = {value:"equation"};
			}
			else{
				par.variety = "reparam";
				par.reparam_eqn_on = false;
				press_button_prop("ModelParamContent","ReparamElement",["name"],par.name);
			}
		}
		break;
		
	case "AddReparamParam2":
		{
			let par = model.param[inter.bubble.k];
			par.variety = "reparam";
			par.reparam_eqn_on = false;
			if(inter.bubble.radio.value == "equation") par.reparam_eqn_on = true;
			close_bubble();
		}
		break;
		
	case "SetDistribution":
		select_bubble_over();
		break;
		
	case "AddDistParam":
		{
			let i = bu.i;
			model.param[i].variety = "dist";
			select_bubble_but_prop("ModelParamContent","DistElement",["i"],i);

			inter.bubble.prior = copy(model.param[i].prior);
			inter.bubble.prior.i = i;
		}
		break;
		
	case "DeleteParamConst":
		inter.help = {title: "Delete constant", te: "Are you sure you want to delete this constant definition?", i:bu.i, ok:"DeleteParamDistConfirm"};
		break;
		
	case "DeleteParamMult":
		inter.help = {title: "Delete multiplier", te: "Are you sure you want to delete this parameter multiplier?", i:bu.i, ok:"DeleteParamMultConfirm"};
		break;
		
	case "DeleteParamReparam":
		inter.help = {title: "Delete reparameterisation", te: "Are you sure you want to delete this reparameterisation?", i:bu.i, ok:"DeleteParamDistConfirm"};
		break;
		
	case "DeleteParamDist":
		inter.help = {title: "Delete parameter distribution", te: "Are you sure you want to delete this parameter distribution?", i:bu.i, ok:"DeleteParamDistConfirm"};
		break;
	
	case "DeleteParamDistConfirm":
		{
			let he = inter.help;
			start_worker("DeleteParamDist",{ i:he.i});
			close_help();
		}
		break;
		
	case "DeleteParamMultConfirm":
		{
			let he = inter.help;
			start_worker("DeleteParamMult",{ i:he.i});
			close_help();
		}
		break;
	
	case "SetDerived":
		{
			select_bubble_over();
			let eqn1 = create_equation("","derive_param");
			let eqn2 = create_equation("","derive_eqn");
			inter.bubble.derived = { eqn1:eqn1, eqn2:eqn2};
			update_param();
		}
		break;
		
	case "AddDerived":
		if(bubble_check_error() == false){
			copy_back_to_source();
			let res = model.add_modify_derived(inter.bubble.derived);
			if(res.err == false){
				close_bubble();
				update_param();
			}
			else output_help(res);
		}
		break;
	
	case "LoadAMarix":
		start_data_source("LoadAmatrix",{},{p:bu.p,i:bu.i});
		file_add_datatable(); 
		break;
	
	case "EditAmatrix":
		start_worker("Edit A matrix",{ i:bu.i, p:bu.p});
		break;
		
	case "EditSimValue":
		edit_sim_value(bu.i,inter.layer[inter.over.layer].name,inter.over.i,bu.source);
		break;
		
	case "EditPriorSplitValue":
		{
			let th = bu.i;
			let par = model.param[th];
			start_worker("Edit PriorSplit",{type:"PriorSplit", label_info:par.label_info, i:th});
		}
		break;
		
	case "EditPriorSplitDone":
		break;
		
	case "EditDistSplitValue":
		{
			let th = bu.i;
			let par = model.param[th];		
		
			start_worker("Edit DistSplit",{type:"DistSplit", source:bu.source, label_info:par.label_info, i:th})
		}
		break;
		
	case "EditReparamValue":
		edit_reparam_value(bu.i,inter.layer[inter.over.layer].name,inter.over.i,bu.source);
		break;
		
	case "EditAlpha":
		select_bubble_over();
		inter.bubble.alpha_val = bu.source.alpha;
		break;
		
	case "DoneAlpha":
		if(bubble_check_error() == false){
			copy_back_to_source();	
			inter.bubble.bu.source.alpha = inter.bubble.alpha_val;
			close_bubble();
		}
		break;
		
	case "EditPriorElement":
		select_bubble_over();
		inter.bubble.prior = copy(model.param[bu.i].prior);
		inter.bubble.prior.i = bu.i;
		break;
	
	case "EditPopPrior":
		select_bubble_over();
		inter.bubble.prior = copy(bu.prior);
		break;
		
	case "DoneCompPrior":
		done_comp_prior();
		break;
		
	case "DonePrior": case "DoneDist":
		done_prior();
		break;
	
	case "EditPriorSplitElement": case "EditDistSplitElement": 
		select_bubble_over();
	
		let prior = get_element(inter.edit_param.prior_split,bu.pindex);
		inter.bubble.prior = copy(prior);
		inter.bubble.prior.i = bu.i;
		inter.bubble.prior.pindex = bu.pindex;
		break;
		
	case "DonePriorSplit": case "DoneDistSplit":
		if(bubble_check_error() == false){
			copy_back_to_source();
			if(invalid_prior() == false){
				update_prior_split();
				close_bubble();
			}
		}
		break;
		
	case "MenuBack":
		inter.bubble.radio = {value:"Compartments"};
		change_bubble_mode(undefined);
		break;
		
	case "EditDerive":
		{
			select_bubble_over();
			let der = model.derive[bu.val];
			inter.bubble.derived = { eqn1:copy(der.eqn1), eqn2:copy(der.eqn2), val:bu.val};
		}
		break;
	
	case "DeleteDerive":
		inter.help = {title: "Delete derived quantity", te: "Are you sure you want to delete this derived quantity?", val:bu.val, ok:"DeleteDeriveConfirm"};
		break;
	
	case "DeleteDeriveConfirm":
		{
			let he = inter.help;
			model.derive.splice(he.val,1);
			update_param();
			close_bubble();
			close_help();
		}
		break;
	
	case "StartSimulation":
		if(check_error_textbox() == false){
			copy_back_to_source();
			let de = model.sim_details;
			if(Number(de.t_start) >= Number(de.t_end)){
				alertp("The start time must be before the end time.");
			}
			else{
				if(check_time_error() == false){	
					switch(model.sim_details.run_local.value){
					case "Yes":
						start_worker("Start",{save_type:"sim", map_store:map_store, ver:ver});
						break;
						
					case "No":
						run_cluster("sim");
						break;
					}
				}
			}
		}
		break;
		
	case "StartPPC":
		if(check_error_textbox() == false){
			copy_back_to_source();
		
			let de = model.ppc_details;
			
			if(Number(de.ppc_t_start) >= Number(de.ppc_t_end)){
				alertp("The start time must be before the end time.");
			}
			else{
				if(Number(de.ppc_t_start) < Number(de.t_start)){
					alertp("The start time cannot be before the inference start time");
				}
				else{
					if(Number(de.ppc_t_start) > Number(de.inf_t_end)){
						alertp("The start time cannot be after the inference end time");
					}
					else{
						if(check_time_error() == false){	
		
							switch(model.ppc_details.run_local.value){
							case "Yes":
								start_worker("StartPPC",{save_type:"ppc", map_store:map_store, ver:ver});
								break;
								
							case "No":
								run_cluster("ppc");
								break;
							}
						}
					}
				}
			}
		}
		break;
	
	case "LoadInference":
		{
			start_worker("Import output");
		}
		break;
		
	case "StartInference":
		if(check_disable() == false){
			if(check_error_textbox() == false){
				copy_back_to_source();
				
				if(check_time_error() == false){
					switch(model.inf_details.run_local.value){
					case "Yes":
						if(model.inf_details.algorithm.value == "PAS-MCMC"){
							alert_help("Problem starting...","The PAS-MCMC algorithm cannot be run on the local machine (because it requires communication between cores).");
						}
						else{
							start_worker("Start",{save_type:"inf", map_store:map_store, ver:ver});
						}
						break;
					
					case "No":
						run_cluster("inf");
						break;
					}
				}
			}
		}
		break;
		
	case "StartCluster":
		start_worker("StartCluster",{save_type:inter.help.siminf, map_store:map_store, ver:ver});		
		break;

	case "CopyText":
		navigator.clipboard.writeText(bu.info.substr(2));
		break;

	case "BayesFactor":
		select_bubble_over();
		inter.bubble.BF_val = "0";
		break;
		
	case "CalculateBF":
		if(check_error_textbox() == false){
			copy_back_to_source();
			inter.graph.calculate_BF();
		}
		break;
		
	case "Invalid":
		select_bubble_over();
		break;
	
	case "Nothing":
		break;
	
	case "ExportScript":
		close_bubble();
		inter.file_store.type = "SaveAs";
		start_worker("Save BICI",{save_type:"export", map_store:map_store});
		break;
		
	case "ParamInfo":
		select_bubble_over();
		break;
	
	case "SelectParam":
		{
			let par = bu.op.par;
		
			switch(par.variety){
			case "normal":
				if(par.type == "derive_param"){
					let der = model.derive;
					
					let fl = false;
					for(let i = 0; i < der.length; i++){
						if(der[i].eqn1.param.length == 1){
							if(der[i].eqn1.param[0].name == par.name){
								fl = true;
								select_bubble_derived(i);
								break;
							}
						}
					}
					
					if(fl == false) error("Could not select");
				}
				break;
				
			case "const": case "reparam": case "dist":
				select_bubble_param(par);
				break;
					
			default: error("sel op error"); break;
			}
		}
		break;
		
	case "SelectTransEqn":
		{
			let info = bu.eqn_info;
			select_bubble_transition(info.p,info.cl,info.i);
		}
		break;
		
	case "SelectDataSpec":
		{
			let info = bu.eqn_info;
			select_bubble_data_spec(info.p,info.i);
		}
		break;
		
	case "SelectDataElement":
		{
			let info = bu.eqn_info;
			select_bubble_data_element(info.p,info.i,info.r,info.c);
		}
		break;
	
	case "SelectReparam":
		{
			let info = bu.eqn_info;
			select_reparam_element(info.par_name,info.index);
		}
		break;
		
	case "SelectDist":
		{
			let info = bu.eqn_info;
			select_dist_element(info.par_name,info.index);
		}
		break;
		
	case "SelectDerived":
		{
			let info = bu.eqn_info;
			select_bubble_derived(info.i);
		}
		break;
	
	case "X Tick":
		select_bubble_over();
		inter.bubble.min = inter.graph.range.xmin;
		inter.bubble.max = inter.graph.range.xmax;
		break;
		
	case "xaxisOK":
		if(bubble_check_error() == false){
			copy_back_to_source();
			if(inter.bubble.min >= inter.bubble.max){
				inter.help = { title:"Ordering problem", te:"Minumum must be smaller than maximum"};
			}
			else{
				inter.graph.range.xmin = inter.bubble.min;
				inter.graph.range.xmax = inter.bubble.max;
				inter.graph.set_axes();
				close_bubble();
			}
		}
		break;
		
	case "Y Tick":
		select_bubble_over();
		inter.bubble.min = inter.graph.range.ymin;
		inter.bubble.max = inter.graph.range.ymax;
		break;
		
	case "yaxisOK":
		if(bubble_check_error() == false){
			copy_back_to_source();
			if(inter.bubble.min >= inter.bubble.max){
				inter.help = { title:"Ordering problem", te:"Minumum must be smaller than maximum"};
			}
			else{
				inter.graph.range.ymin = inter.bubble.min;
				inter.graph.range.ymax = inter.bubble.max;
				inter.graph.set_axes();
				close_bubble();
			}
		}
		break;
		
	case "CombineIEBubble":
		select_bubble_over();
		inter.bubble.ie_p = bu.p;
		inter.bubble.ie_combine = bu.i;
		break;
		
	case "CombineIE":
		start_worker("Combine IE",{p:inter.bubble.ie_p, e:inter.bubble.ie_combine, i:bu.i});
		close_bubble();		
		break;
		
	case "SplitIE":
		inter.help = { title:"Split individual effect group", te:"Are you sure you want to split this group of individual effects?", ok:"SplitIE2", p:bu.p, i:bu.i};
		break;
		
	case "SplitIE2":
		start_worker("Split IE",{p:inter.help.p, i:inter.help.i});	
		close_help();
		break;
		
	case "EditXvectorElement":
		select_bubble_over();
		break;
	
	case "EditAmatrixElement":
		select_bubble_over();
		break;
	
	case "DoneAmatrix":
		{
			let ea = inter.edit_Amatrix;
			start_worker("Set A matrix",{p:ea.p, i:ea.i, value:ea.A_value });
		}
		break;
		
	case "EditXvector":
		start_worker("Edit Xvector",{p:bu.p, i:bu.i});
		break;
		
	case "DoneXvector":
		{
			let ea = inter.edit_Xvector;
			start_worker("Set Xvector",{p:ea.p, i:ea.i, value:ea.X_value });
		}
		break;
		
	case "LoadXvector":
		start_data_source("Fixed Effect",{},{load_file:"Comp", p:bu.p, i:bu.i});
		file_add_datatable(); 	
		break;
	
	case "LoadXvector2":
		start_worker("Load Xvector",{edit_source:edit_source});
		break;
		
	case "AddFilter":
		select_bubble_over();
		break;
	
	case "AddFilter2":
		add_pop_filter(bu.rpf,bu.cl);
		break;
		
	case "PlayButton":
		inter.graph.press_play();
		break;
		
	case "PlayForward":
		inter.graph.change_frame(1);
		break;
		
	case "PlayBackward":
		inter.graph.change_frame(-1);
		break;
			
	case "Timebar":
		inter.graph.press_timebar(bu);
		break;
		
	case "Settings":
		select_bubble_over();
		break;
		
	case "Defaultaxis":
		inter.graph.range = copy(inter.graph.range_default);
		inter.graph.initialise_axes();
		close_bubble();
		break;

	case "AddConstantTime":
		select_bubble_over();
		inter.bubble.radio = {value:"start"};
		inter.bubble.fixed_time = "";
		break;
		
	case "AddFixedtimeColumn":
		if(bubble_check_error() == false){	
			copy_back_to_source();	
			if(inter.bubble.radio.value == "start") add_constant_column("start");
			else add_constant_column(inter.bubble.fixed_time);
			close_bubble();
		}	
		break;
	
	case "AddConstantColumn":
		add_constant_column(bu.info);
		break;
	
	case "AddConstantValue":
		select_bubble_over();
		inter.bubble.const_column = "";
		break;
		
	case "AddConstantValue2":
		if(bubble_check_error() == false){	
			copy_back_to_source();	
			add_constant_column(inter.bubble.const_column);
			close_bubble();
		}	
		break;
		
	case "ViewParam":
		{
			let par = bu.source.param[bu.i];
			start_worker("View Param",{i:bu.i, par:par, source:bu.source, pos_view:bu.pos_view, sel_view:copy(bu.pos_view[0])});
		}
		break;
		
	case "CancelViewGraph":
		close_view_graph();
		break;
		
	case "CancelViewInd":
		close_view_ind();
		break;
		
	case "ViewCode":
		close_bubble();
		start_worker("View Code",{save_type:"export", map_store:map_store});
		break;
	
	case "GraphSettings":
		select_bubble_over();
		break;
		
	case "TraceSettings":
		select_bubble_over();
		inter.bubble.burnin = Math.floor(model.inf_res.burnin_factor*100);
		break;
		
	case "ScatterSettings":
		select_bubble_over();
		break;
		
	case "BurninOK":
		if(bubble_check_error() == false){			
			copy_back_to_source();
			let source = model.inf_res;
			source.burnin_factor = Number(inter.bubble.burnin)/100;
			source.burnin = Math.floor(source.burnin_factor*source.samp_max);
			close_bubble();
			inter.graph.init = undefined;
		}
		break;
	
	case "ScatterOK": case "DistributionPlotOK":
		if(bubble_check_error() == false){			
			copy_back_to_source();
			close_bubble();
			inter.graph.init = undefined;
		}	
		break;
	
	case "ExportGraph":
		close_bubble();
		saving_dialogue("",".png","Export graph");
		break;
		
	case "PrintGraph":
		close_bubble();
		inter.graph.export_image();
		break;
	
	case "ExportTable":
		close_bubble();
		saving_dialogue("",".csv","Export table");
		break;

	case "ExportTableContent":
		close_bubble();
		saving_dialogue("",".csv","Export table content");
		break;
		
	case "ExportMatrixTable":
		close_bubble();
		saving_dialogue("",".csv","Export matrix table");
		break;
		
	case "CopyData":
		if(model.sim_res.plot_filter.species[bu.op.p].gen_source.length == 0){
			inter.help = {title: "Copy Data", te: "This copies the simulated data as real data, such that inference can be performed on it. To do this, however, data first has to be added."};
		}
		else{	
			inter.help = {title: "Copy Data", te: "This copies the simulated data as real data, such that inference can be performed on it.\nAre you sure you want to do this?", op:bu.op, ok:"CopyData2"};
		}
		break;
		
	case "CopyData2":
		{
			start_worker("Copy simulated");
			close_bubble();
		}
		break;
		
	case "Options":
		inter.options = true;
		break;
		
	case "OptionsDone":
		if(check_error_textbox() == false){
			copy_back_to_source();
			inter.options = false;
		}
		break;
	
	case "Stop":
		if(inter.running_status == true) terminate();
		else stop_worker();
		break;
		
	case "MatrixEleBut":	
		select_bubble_over();
		break;
	
	case "PopFilt":
		select_bubble_over();
		inter.bubble.popfilt = {rpf:bu.rpf, filter:bu.filter};
		break;
		
	case "RemoveFilter":
		{
			let op = bu.op;
			remove_filter(op.p,op.num,op.rpf);
		}
		break;
		
	case "DoneFilter":
		done_filter();
		break;
		
	case "Slider":
		inter.graph.press_slider(bu,true);
		break;
	
	case "CompVectorClick":
		select_bubble_over();
		break;
	
	case "ProbEqn":
		start_equation(bu.eqn.te,bu.eqn,{ref:"ProbEqn", te_store:bu.eqn.te, eqn:bu.eqn},undefined);
		break;
	
	case "DeleteBranchProb":
		model.delete_branch_prob(bu.p,bu.cl,bu.i);
		break;
	
	case "SelectInd":
		inter.graph.ind_sel = bu.te;
		inter.graph.init = undefined;
		break;
		
	case "CopyPic":
		select_bubble_over();
		inter.bubble.suffix = "";
		break;

	case "DeleteSel":
		inter.help = { title:"Delete", te:"Are you sure you want to delete the selected compartemnts?", ok:"DeleteSel2"};
		break;
		
	case "DeleteSel2":
		model.delete_selected();
		close_bubble();
		close_help();
		break;

	case "DoneCopyPic":
		if(check_error_textbox() == false){
			copy_back_to_source();
			
			model.copy_selected(inter.bubble.suffix);
			close_bubble();
		}
		break;
	
	case "Selected":
		break;
		
	case "CopyComp":
		model.copy_selected("Copy");
		close_bubble();
		generate_screen();
		model.select_single();
		break;
		
	case "ExampleModel":
		if(model.start == true){ 
			inter.help = { title:"Load example", te:"Are you sure you want to discard the currently loaded model?", ok:"LoadExample", file:bu.file};
		}
		else{
			start_worker("Load Example",bu.file);
		}
		break;
			
	case "Run Anyway":
		start_spawn(inter.help.file_list);
		close_help();
		break;
			
	case "LoadExample":
		start_worker("Load Example",inter.help.file);
		close_help();
		break;
		
	case "ModType": case "ShowModel":
		select_bubble_over();
		break;
	
	case "OutputStates":
		start_worker("Output States",{name:tab_name()});	
		break;
		
	case "OutputParams":
		start_worker("Output Params",{name:tab_name()});	
		break;
		
	default:
		error(ac+" Action not recognised"); 
		break;
	}
	
	if(false) check_param();                       // This is used for diagnostics
	
	generate_screen();
}

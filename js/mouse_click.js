"use strict";

/// Fires when the mouse is clicked 
function mouse_click(x,y)                                 
{
	inter.mx = x; inter.my = y;
	
	let l = inter.over.layer;
	let i = inter.over.i;
		
	if(l != undefined){
		button_action(inter.layer[l].but[i],"click");
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
		case "Help": case "HelpBackground": case "Yscroll": case "Xscroll": break;
		
		default: close_bubble(); break;
		}	
	}
	
	if(ac != "Dropdown") close_dropdown();
	
	//comment(ac+" action  "+bu.type+" type");	
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
			let tree = inter.page_name.split("->");
			if(tree[0] == "Model" && tree[1] == "Compartments"){
				select_bubble_over();
				let p = model.get_p();
				inter.bubble.radio = {value:model.species[p].type}
			}
		}
		else{
			change_page({susu:val});
		}
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

	case "Example": 
		if(model.start == true){ 
			inter.help = { title:"Load example", te:"Are you sure you want to discard the currently loaded model?", ok:"LoadExample", file:bu.file};
		}
		else{
			load_example(bu.file);
		}
		break;

	case "LoadExample":
		load_example(inter.help.file);
		close_help();
		break;

	case "EditDescription":
		model.description.edit = true;
		reset_text_box();
		generate_screen();
		set_focus_first();
		break;
		
	case "CancelEditDescription":
		model.description.edit = false;
		break;
		
	case "DoneEditDescription":
		model.description.edit = false;
		copy_back_to_source(inter.over.layer);
		break;
		
	case "AddCompartment":
		inter.mode={ type:"Add_Compartment"};
		close_bubble();
		break;
	
	case "ImportCompMap":
		loading_dialogue("",".geojson,.json","GEOJSON_compartment"); 
		close_bubble();
		break;
	
	case "ImportComp":
		inter.bubble.radio = {value:"SelectColour", col:WHITE};
		inter.bubble.pos_check = {check:true};
		change_bubble_mode("LoadComps");
		break;
		
	case "AddTransition":
		inter.mode={ type:"Add_Transition", i:undefined, f:undefined, midp:[]};
		close_bubble();
		break;
	
	case "ImportTrans":
		inter.bubble.trans_type = {te:dist_pos[0]};
		change_bubble_mode("LoadTrans");
		break;
		
	case "ImportTrans2":
		{
			let p = model.get_p();
			let cl = model.get_cl();
			start_data_source("Trans "+inter.bubble.trans_type.te,{},{load_file:"Trans", p:p, cl:cl});
			file_add_datatable(); 
			close_bubble();
		}
		break;
		
	case "AddSource":
		inter.mode={ type:"Add_Source", midp:[]};
		close_bubble();
		break;
		
	case "ImportSource":
		inter.bubble.trans_type = {te:exp_dist_pos[0]};
		inter.bubble.pos_check = {check:false};
		change_bubble_mode("LoadSource");
		break;
		
	case "ImportSource2":
		{
			let p = model.get_p();
			let cl = model.get_cl();
			let ty = "Source ";
			if(inter.bubble.pos_check.check == true) ty += "Pos ";
			ty += inter.bubble.trans_type.te;
			start_data_source(ty,{},{load_file:"Source", p:p, cl:cl});
			file_add_datatable(); 
			close_bubble();
		}
		break;
	
	case "AddSink":
		inter.mode={ type:"Add_Sink", midp:[]};
		close_bubble();
		break;
	
	case "ImportSink":
		inter.bubble.trans_type = {te:dist_pos[0]};
		inter.bubble.pos_check = {check:false};
		change_bubble_mode("LoadSink");
		break;
	
	case "ImportSink2":
		{
			let p = model.get_p();
			let cl = model.get_cl();
			let ty = "Sink ";
			if(inter.bubble.pos_check.check == true) ty += "Pos ";
			ty += inter.bubble.trans_type.te;		
			start_data_source(ty,{},{load_file:"Sink", p:p, cl:cl});
			file_add_datatable(); 
			close_bubble();
		}
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
				add_trans = { i:"Source", f:bu.i};	
				break;
				
			case "Add_Sink":
				add_trans = { i:bu.i, f:"Sink"};	
				break;
				
			default:
				select_bubble_over();
				
				let co = model.species[bu.p].cla[bu.cl].comp[bu.i];
				inter.bubble.checkbox = {check:co.markov_branch};
				break;
			}
			
			if(add_trans.i != undefined){
				mode_off();
				if(model.add_transition(bu.p,bu.cl,add_trans.i,add_trans.f,mo.midp,dist_pos[0]) == true){
					model.select_transition_bubble(bu.p,bu.cl,model.species[bu.p].cla[bu.cl].ntra-1);
				}			
			}
		}
		break;
	
	case "Transition":
		select_bubble_over();
		break;
		
	case "TransitionPoint":
		if(bu.variety == "Midpoint"){
			let mo = inter.mode;
			if(mo.type == "Add_Transition" || mo.type == "Add_Source" || mo.type == "Add_Sink"){
				mo.midp.splice(bu.index,1);
			}
			else model.delete_transitionpoint(bu.p,bu.cl,bu.i,bu.index);
		}
		break;
		
	case "CloseBubble":
		close_bubble();
		break;
	
	case "ColourSelect":
		let sel_bu = bu.sel_bu;
		if(sel_bu.type == "Menu") inter.bubble.radio.col = bu.col;
		else model.species[sel_bu.p].cla[sel_bu.cl].comp[sel_bu.i].col = bu.col;
		break;
		
	case "ChangeComp":
		if(bubble_check_error() == false){
			copy_back_to_source();
			let bu = inter.bubble.bu;
			let co = model.species[bu.p].cla[bu.cl].comp[bu.i];
			if(co.choose_branch == true){
				co.markov_branch = inter.bubble.checkbox.check;
				model.determine_branching();
				model.update_pline(bu.p,bu.cl);
			}
			
			// Make sure all the clones are up-to-date
			let list = model.find_clones(bu.p,bu.cl);
		
			if(list.length > 0){
				for(let li of list){
					let co2 = model.species[li.p].cla[li.cl].comp[bu.i];
					
					co2.col = co.col;
					co2.fixed.check = co.fixed.check;
					
					if(co2.name != co.name){
						model.rename_compartment(p,cl,bu.i,co.name);
					}
				}
				model.check_clones();
			}
			
			close_bubble();
		}
		break;
		
	case "DeleteComp":
		inter.help = {title: "Delete Compartment", te: "Are you sure you want to delete this compartment?", p:bu.op.p, cl:bu.op.cl, i:bu.op.i, ok:"DeleteCompConfirm"};
		break;

	case "DeleteCompConfirm":
		{
			let he = inter.help;
			model.delete_compartment(he.p,he.cl,he.i,true);
			close_help();
			close_bubble();
		}
		break;
		
	case "DeleteTrans":
		inter.help = {title: "Delete Transition", te: "Are you sure you want to delete this transition?", 
		              p:bu.op.p, cl:bu.op.cl, i:bu.op.i, ok:"DeleteTransConfirm"};
		break;
	
	case "DeleteTransConfirm":
		{
			let he = inter.help;
			model.delete_transition(he.p,he.cl,he.i);
			close_help();
			close_bubble();
		}
		break;
		
	case "ClassificationBack":
		if(double_click == true) zoom_double_click(bu.p,bu.cl);
		break;
		
	case "MessageClose":
		mode_off();
		//model.determine_branching();
		//model.update_pline_all();
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
						model.add_latlng_compartment(name,mo.p,mo.cl,mo.mp.x,mo.mp.y,BLACK);
						set_ruler(mo.p,mo.cl);
					}
					else{
						ty = "Compartment";
						model.add_compartment(name,mo.p,mo.cl,mo.mp.x,mo.mp.y,WHITE,false,true);
					}
					
					mode_off();
					model.select_button_bubble("Compartment",ty,mo.p,mo.cl,model.species[mo.p].cla[mo.cl].ncomp-1);
				}
				break;
			
			case "Add_Label":
				model.add_label(mo.te,model.get_p(),model.get_cl());
				mode_off();
				break;
			}
		}
		break;
		
	case "AddSpecies":
		select_bubble_over();
		inter.bubble.name = "";
		inter.bubble.radio = {value:"Individual"};
		break;
		
	case "AddSpeciesOK":
		if(bubble_check_error() == false){	
			copy_back_to_source();		
			let result = model.add_species(inter.bubble.species_name,inter.bubble.radio.value);
			
			if(result == "success"){
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
				inter.bubble.checkbox = {check:false};
			}
		}
		break;
		
	case "AddClassificationOK":
		if(bubble_check_error() == false){
			copy_back_to_source();		
			
			let pag = inter.page[inter.pa];
			let p = pag.sub[pag.index].index;
		
			let bub = inter.bubble;
			let result = model.add_classification(p,bub.classification_name,bub.drop.te,{coord:bub.radio.value, default_map:bub.checkbox.check});
			
			if(result == "success"){
				model.load_default_map();
				
				close_bubble();
				initialise_pages();
			
				change_page({sususu:model.species[p].ncla-1});
			}
		}
		break;
		
	case "RadioButton":
		bu.source.value = bu.value;
		inter.graph.init = undefined;
		break;
		
	case "DeleteSpecies":
		inter.help = { title:"Delete", te:"Are you sure you want to delete this species?", ok:"DeleteSpecies2"};
		break;
	
	case "DeleteSpecies2":
		{
			let p = model.get_p();
			model.delete_species(p);
			change_page({});
		}
		break;

	case "DeleteClass":
		inter.help = { title:"Delete", te:"Are you sure you want to delete this classification?", ok:"DeleteClass2"};
		break;

	case "DeleteClass2":
		{
			let p = model.get_p();
			let cl = model.get_cl();
			model.delete_classification(p,cl);
			change_page({});
		}
		break;

	case "ChangeClassName": 
		if(bubble_check_error() == false){
			copy_back_to_source();
			let index_new = inter.bubble.drop.te;
			let claa = model.get_cla();
			
			if(index_new != claa.index){
				model.rename_index(index_new,model.get_p(),model.get_cl(),true);
			}

			close_bubble();
			initialise_pages();
		}

	case "ChangeSpeciesName":
		if(bubble_check_error() == false){
			model.species[model.get_p()].type = inter.bubble.radio.value;
			
			copy_back_to_source();
			close_bubble();
			initialise_pages();
		}
		break;
	
	case "Dropdown":
		if(inter.dropdown.source == bu.source) inter.dropdown = {};
		else inter.dropdown = {source:bu.source, style:bu.style};
		break;
		
	case "DropdownOption":
		{
			bu.source.te = bu.pos.te;
			let pos = bu.pos;
			
			if(pos.i != undefined) bu.source.i = pos.i;
			if(pos.p != undefined) bu.source.p = pos.p;
			if(pos.cl != undefined) bu.source.cl = pos.cl;
			if(pos.tr != undefined) bu.source.tr = pos.cl;
			if(pos.c != undefined) bu.source.c = pos.c;
			if(pos.param != undefined) bu.source.param = pos.param;
			if(pos.list != undefined) bu.source.list = pos.list;
			if(pos.radio != undefined) bu.source.radio = pos.radio;
			if(pos.radioy != undefined) bu.source.radioy = pos.radioy;
				
			if(bu.source == inter.data_type) close_bubble();			
			inter.dropdown = {};
				
			if(find_in(prior_pos,bu.pos.te) != undefined){
				inter.bubble.find_focus = true;
			}
			
			if(pos.ac == "Remove filter") pos.filter.splice(pos.num,1);
		
			inter.graph.init = undefined;
		}
		break;
		
	case "ChangeTransValue":
		if(bubble_check_error() == false){
			model.species[bu.op.p].cla[bu.op.cl].tra[bu.op.i].type = inter.bubble.trans_type.te;

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
				let func = false;
				if(paste == "exp" || paste == "log" || paste == "cos" || paste == "sin"){ func = true; paste += "()";}
				cursor_paste(paste);
				if(func == true) inter.cursor.i--;
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
		if(bu.source.check == false){
			bu.source.check = true;
		}
		else bu.source.check = false;
	
		if(bu.value && bu.value.dist_matrix){
			if(bu.source.check == true)	set_dist(bu.value);
			else bu.value.set = false;
		}
		
		inter.graph.init = undefined;
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
		
	case "AddModParamEqn":
		equation_add_modelparam(bu.i);
		break;
		
	case "Menu":
		select_bubble_over();
		break;	
		
	case "Load":	
		close_bubble();
		loading_dialogue("",".bici","BICI_file");
		break;
		
	case "Save":
		save_file(model.filename,"BICI_file"); 
		break;
		
	case "SaveAs":
		close_bubble();
		saving_dialogue("",".bici","BICI_file");
		break;
		
	case "Import":
		change_bubble_mode("Import Menu");
		break;
		
	case "ImportScript":
		close_bubble();
		loading_dialogue("",".txt","Import file");
		break;
		
	case "NewModel":
		model.start_new();
		change_page({pa:"Model"});
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
		
	case "DoneEquation": 
		if(bubble_check_error() == false){
			copy_back_to_source();
			close_bubble();
			init_param();
		}
		break;
		
	case "DeleteSource":	
		if(bu.source[bu.i].error == true){  // If error is source then delete immediately
			bu.source.splice(bu.i,1);
		}
		else{
			inter.help = { title:"Delete", te:"Are you sure you want to delete this?", ok:"DeleteSource2", source:bu.source, i:bu.i};
		}
		break;
	
	case "DeleteSource2":
		inter.help.source.splice(inter.help.i,1);
		close_help();
		break;
		
	case "DeleteSource3":
		bu.source.splice(bu.i,1);
		break;
	
	case "AddDataSource":
		if(data_source_check_error() == false){
			data_source_add_and_show();
		}
		break;
		
	case "EditSource":
		data_source("Edit",bu.info);
		break;
		
	case "ViewSource":
		data_source("Edit",bu.info);
		inter.edit_source.edit_spec = true;
		select_bubble_over();
		break;
		
	case "PlaceSourceBack":
		if(data_source_check_error() == false){
			data_source("Replace");
		}
		break;
		
	case "UploadDataTable":
		{
			if(inter.edit_source.type == "LoadAmatrix" && inter.radio_heading.value == "No"){
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
				if(inter.edit_source.type == "LoadAmatrix") file_type = "A matrix";
				
				loading_dialogue("",pos,file_type);
			}
		}
		break;
		
	case "UsePreloadedTable":
		inter.edit_source.data_table_use = bu.i;
		inter.edit_source.table_loaded = true;
		close_help();
		break;
	
	case "SelectColumn":
		transfer_column(bu.c);
		break;
		
	case "BackEditSource":
		transfer_column_back();
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
		disactivated();
		//move_data(bu); 
		break;
	
	case "FixedEffect": 
		fixed_effect_data(bu); 
		break;
	
	case "AddInitialPopulation":
		initpop_data(bu); 
		break;
	
	case "AddInitialPopulation2": initpop_data2(); break;

	case "IndTrappedData":
		disactivated();
		//start_data_source("Ind. Trapped",{},bu.op.info);
		//file_add_datatable();
		break;

	case "EditTableElement": case "EditTableHead":
		select_bubble_over();
		break;
	
	case "EditParamElement": 
		select_bubble_over();
		break;
	
	case "EditSplineKnots":
		select_bubble_over();
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
			let c = inter.bubble.bu.c;
			let tab = inter.edit_source.table;
			let num = 0;
			
			let find = inter.bubble.find;
			let replace = inter.bubble.replace;
		
			for(let r = 0; r < tab.nrow; r++){
				if(tab.ele[r][c] == find){
					tab.ele[r][c] = replace;
					num++;
				}
			}
			
			inter.bubble.num = num;
			change_bubble_mode("ReplaceResult");
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
		inter.edit_source.table.edit = true;
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
		disactivated();
		//initpopprior_data(bu);
		break;
		
	case "AddInitPopPrior2":
		initpopprior_data2(); 
		break;

	case "AddDataTable":
		if(bubble_check_error() == false){
			copy_back_to_source();
			set_loadcol();
			file_add_datatable(); 
		}
		break;
		
	case "DataTableBack":
		if(bubble_check_error() == false){
			copy_back_to_source();
			close_bubble();
			data_source("Replace");
		}
		break;

	case "CompData": comp_data(bu);	break;
	case "CompView": comp_data_view(bu); break;
	case "TransData": trans_data(bu,"Transition"); break;	
	case "SourceData": trans_data(bu,"Source"); break;	
	case "SinkData": trans_data(bu,"Sink"); break;	
	case "DiagTestData": diagtest_data(bu); break;
	case "PopulationData": population_data(bu); break;
	case "PopTransData": poptrans_data(bu); break;
	case "SetTrapsData": 
		disactivated();
		//settraps_data(bu); 
		break;
		
	case "SeqData":
		disactivated();
		//start_data_source("Genetic",{},bu.op.info);
		//file_add_datatable();
		break;

	case "TableSlider": break;
		
	case "CloneClassInit":
		clone_class_init(bu.op.p);
		break;
	
	case "CloneClass":
		{
			let drop = inter.bubble.drop;
			model.clone_class(drop.p_sel,drop.p,drop.cl);
			//change_page({sususu:cl});
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
		
	case "AddAnnotationOK":
		switch(inter.bubble.radio.value){
		case "Label": 
			change_bubble_mode("add_label");
			inter.bubble.label = { te:""};
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
			inter.mode={ type:"Add_Label", te:inter.bubble.label.te};
			close_bubble();
		}
		break;
		
	case "DeriveOK":
		if(bubble_check_error() == false){
			copy_back_to_source();
			if(model.add_modify_derived(inter.bubble.derived,"modify") == "success"){
				close_bubble();
			}
		}
		/*
		if(bubble_check_error() == false){
			copy_back_to_source();
			close_bubble();
		}
		*/
		break;
		
	case "DeleteAnnotation":
		delete_annotation();
		close_bubble();	
		break;
		
	case "LabelText":
		select_bubble_over();
		inter.bubble.label = {te:bu.te};
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
		{
			if(data_source_check_error() == false){
				let so = inter.edit_source;
				let tab = so.table;
				let col = so.info.col;
				model.add_file_compartment(so.info.p,so.info.cl,tab,col);
				close_data_source();
			}
		}
		break;
		
	case "AddFileTrans":
		{
			if(data_source_check_error() == false){
				let so = inter.edit_source;
				let tab = so.table;
				let type = so.type.substr(6);
				model.add_file_transition(so.info.p,so.info.cl,tab,type);
				close_data_source();
			}
		}
		break;
		
	case "AddFileSource":
		{
			if(data_source_check_error() == false){
				let so = inter.edit_source;
				let tab = so.table;
				let type = so.type.substr(7);
				model.add_file_source(so.info.p,so.info.cl,tab,type);
				close_data_source();
			}
		}
		break;
	
	case "AddFileSink":
		{
			if(data_source_check_error() == false){
				let so = inter.edit_source;
				let tab = so.table;
				let type = so.type.substr(5);
				model.add_file_sink(so.info.p,so.info.cl,tab,type);
				close_data_source();
			}
		}
		break;
		
	case "Box":
		select_bubble_over();
		break;
		
	case "CloseBubble":
		close_bubble();
		break;	
		
	case "AddCompMap":
		add_compartment_map();
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
			model.param[ep.i].value = ep.value;
			model.param[ep.i].set = true;
			
			branch_auto_calculate();
			
			close_bubble();
			close_param_source();
		}
		break;
		
	case "EditPriorSplitDone": case "EditDistSplitDone": 
		if(edit_param_check_error() == false){
			let ep = inter.edit_param;
			model.param[ep.i].prior_split = ep.prior_split;
	
			model.param[ep.i].prior_split_set = true;
			
			close_bubble();
			close_param_source();
			if(ac == "EditDistSplitDone") init_param();
		}
		break;
		
	case "EditReparamDone":
		if(edit_param_check_error() == false){
			let ep = inter.edit_param;
			model.param[ep.i].value = ep.value;
			model.param[ep.i].set = true;
			close_bubble();
			close_param_source();
			init_param();
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
		if(data_source_check_error() == false){
			load_tensor();
		}
		break;
	
	case "LoadReparamDone":
		if(data_source_check_error() == false){
			load_reparam();
		}
		break;
		
	case "LoadPriorSplitDone": case "LoadDistSplitDone":
		if(data_source_check_error() == false){
			load_priorsplit();
		}
		break;
		
	case "SetConstant": case "SetReparam":
		select_bubble_over();
		break;
		
	case "AddConstParam":
		{
			generate_screen();
			let par = model.param[bu.i];
			par.variety = "const";
			press_button_prop("ModelParamContent","ParamSimElement",["name"],par.name);
		}
		break;
		
	case "AddReparamParam":
		{
			let par = model.param[bu.i];
			par.variety = "reparam";
			press_button_prop("ModelParamContent","ReparamElement",["name"],par.name);
		}
		break;
		
	case "SetDistribution":
		select_bubble_over();
		//inter.bubble.sel_param = "";
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
		
	case "DeleteParamReparam":
		inter.help = {title: "Delete reparameterisation", te: "Are you sure you want to delete this reparameterisation?", i:bu.i, ok:"DeleteParamDistConfirm"};
		break;
		
	case "DeleteParamDist":
		inter.help = {title: "Delete parameter distribution", te: "Are you sure you want to delete this parameter distribution?", i:bu.i, ok:"DeleteParamDistConfirm"};
		break;
	
	case "DeleteParamDistConfirm":
		{
			let he = inter.help;
			model.param[he.i].variety = "normal";
			init_param();
			close_help();
		}
		break;
	
	case "SetDerived":
		{
			select_bubble_over();
			let eqn1 = create_equation("","derive_param");
			let eqn2 = create_equation("","derive_eqn");
			inter.bubble.derived = { eqn1:eqn1, eqn2:eqn2};
			//init_param();
		}
		break;
		
	case "AddDerived":
		if(bubble_check_error() == false){
			copy_back_to_source();
			if(model.add_modify_derived(inter.bubble.derived) == "success"){
				close_bubble();
			}
		}
		break;
	
	case "EditSimValue":
		edit_sim_value(bu.i,inter.layer[inter.over.layer].name,inter.over.i,bu.source);
		break;
		
	case "EditPriorSplitValue":
		{
			let th = bu.i;
			let par = model.param[th];		
			inter.edit_param = {type:"PriorSplit", prior_split:copy(par.prior_split), value:get_prior_split_tensor(par), label_info:par.label_info, i:th};
		}
		break;
		
	case "EditDistSplitValue":
		{
			let th = bu.i;
			let par = model.param[th];		
			inter.edit_param = {type:"DistSplit", prior_split:copy(par.prior_split), value:get_prior_split_tensor(par), label_info:par.label_info, i:th};
		}
		break;
		
	case "EditReparamValue":
		edit_reparam_value(bu.i,inter.layer[inter.over.layer].name,inter.over.i);
		break;
		
	case "EditPriorElement":
		select_bubble_over();
		inter.bubble.prior = copy(model.param[bu.i].prior);
		inter.bubble.prior.i = bu.i;
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
		
	case "DoneKnots":
		if(bubble_check_error() == false){
			copy_back_to_source();
			close_bubble();
		}
		break;
		
	case "MenuBack":
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
			init_param();
			close_help();
		}
		break;
	
	case "StartSimulation":
		if(check_error_textbox() == false){
			copy_back_to_source();
			let de = model.sim_details;
			if(Number(de.t_start) >= Number(de.t_end)) alertp("The start time must be before the end time.");
			else 
			if(check_time_error() == false){		
				start_simulation();
			}
		}
		break;
		
	case "StartInference":
		if(check_disable() == false){
			if(check_error_textbox() == false){
				copy_back_to_source();
				
				if(check_time_error() == false){
					start_inference();
				}
			}
		}
		break;

	case "Invalid":
		select_bubble_over();
		//inter.help = { title:"Invalid data source", te:"Because of changes to the model this data source is not longer valid. Please delete and reload."};
		break;
	
	case "Nothing":
		break;
	
	case "ExportScript":
		saving_dialogue("",".txt","Export script");
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
	
	case "Auto":
		set_branch_auto(bu.i);
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
		add_ie_to_group(inter.bubble.ie_p,inter.bubble.ie_combine,bu.i);
		close_bubble();
		break;
		
	case "SplitIE":
		inter.help = { title:"Split individual effect group", te:"Are you sure you want to split this group of individual effects?", ok:"SplitIE2", p:bu.p, i:bu.i};
		break;
		
	case "SplitIE2":
		split_ie(inter.help.p,inter.help.i);
		close_help();
		break;
		
	case "LoadAMarix":
		start_data_source("LoadAmatrix",{},{p:bu.p,i:bu.i});
		file_add_datatable(); 
		break;
	
	case "EditAmatrix":
		{
			let ieg = model.species[bu.p].ind_eff_group[bu.i];
			inter.edit_Amatrix = {value:copy(ieg.A_matrix.value), i:bu.i, p:bu.p};	
		}
		break;
		
	case "EditXvector":
		{
			let feg = model.species[bu.p].fix_eff[bu.i];
			inter.edit_Xvector = {value:copy(feg.X_vector.value), p:bu.p, i:bu.i};	
		}
		break;
		
	case "EditXvectorElement":
		select_bubble_over();
		break;
	
	case "EditAmatrixElement":
		select_bubble_over();
		break;
	
	case "DoneAmatrix":
		done_Amatrix();
		break;
		
	case "DoneXvector":
		done_Xvector();
		break;
		
	case "LoadXvector":
		start_data_source("Fixed Effect",{},{load_file:"Comp", p:bu.p, i:bu.i});
		file_add_datatable(); 	
		break;
	
	case "LoadXvector2":
		if(data_source_check_error() == false){
			inf_eff_load_X();
		}
		break;
		
	case "AddFilter":
		select_bubble_over();
		inter.bubble.sel_filter = {te:select_str};
		break;
	
	case "AddFilter2":
		add_pop_filter(bu.op.result);
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
		//add_constant_column(bu.info);
		break;
		
	case "AddConstantValue2":
		if(bubble_check_error() == false){	
			copy_back_to_source();	
			add_constant_column(inter.bubble.const_column);
		}	
		break;
		
	case "ViewSpline":
		inter.view_graph = {type:"ViewSpline", source:bu.source, title:"Spline profile "+bu.source.param[bu.i].full_name, i:bu.i};
		break;
		
	case "ViewVector":
		inter.view_graph = {type:"ViewVector", source:bu.source, title:"Parameter vector "+bu.source.param[bu.i].full_name, i:bu.i};
		break;
		
	case "ViewMatrix":
		inter.view_graph = {type:"ViewMatrix", source:bu.source, title:"Parameter matrix "+bu.source.param[bu.i].full_name, i:bu.i};
		break;
		
	case "CancelViewGraph":
		close_view_graph();
		break;
		
	case "ViewCode":
		{
			close_bubble();
			
			let type = "sim"; if(inter.page[inter.pa].name == "Inference") type = "inf";
			
			let te = create_output_file(type,"NO_SAVE");
			
			if(te){
				te = add_escape_char(te);	
		
				let	lines = te.split('\n');

				let formatted = get_formatted(lines);
				inter.help = { title:"Source code", te:"This is the source code when exported (it includes information about simulation and inference).", st:"", scroll_to_line:false, script:formatted};
			}
		}
		break;
	
	case "GraphSettings":
		select_bubble_over();
		break;
		
	case "TraceSettings":
		select_bubble_over();
		inter.bubble.burnin = Math.floor(inf_result.burnin_factor*100);
		inter.bubble.find = "hello";
		break;
		
	case "BurninOK":
		if(bubble_check_error() == false){
			copy_back_to_source();
			inf_result.burnin_factor = Number(inter.bubble.burnin)/100;
			calculate_burnin(inf_result);
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
		
	case "ExportMatrixTable":
		close_bubble();
		saving_dialogue("",".csv","Export matrix table");
		break;
		
	case "CopyData":
		if(sim_result.species[bu.op.p].gen_source.length == 0){
			inter.help = {title: "Copy Data", te: "This copies the simulated data as real data, such that inference can be performed on it. To do this, however, data first has to be added."};
		}
		else{
			inter.help = {title: "Copy Data", te: "This copies the simulated data as real data, such that inference can be performed on it.\nAre you sure you want to do this?", op:bu.op, ok:"CopyData2"};
		}
		break;
		
	case "CopyData2":
		{
			copy_simulated_data();
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
		break;
		
	default:
		error(ac+" Action not recognised"); 
		break;
	}
	
	if(false) check_param(); // This is used for diagnostics
	
	generate_screen();
}


"use strict";
// Functions related to the Javascript interface

/// Change the page being viewed
function change_page(page_new)  
{
	if(model && model.description && model.description.edit == true){
		alert_help("Close editing","Please finish editing the model description.");
		return;
	}
	
	if(inter.edit_source == "editted"){
		alert_help("Data source changed","Click 'Update' to revise data source otherwise 'Cancel'.");
		return;
	}
	
	if(inter.edit_source == "adding"){
		alert_help("Add data source","Click 'Done' to add data source otherwise 'Cancel'.");
		return;
	}
	
	let pa_new = page_new.pa;
	if(pa_new != undefined){
		if(isNaN(pa_new)){
			pa_new = find(inter.page,"name",pa_new);
			if(pa_new == undefined){ alertp("Problem changing page"); return;}
		}
		inter.pa = pa_new;
	}
	
	let su_new = page_new.su;
	if(su_new != undefined){
		if(isNaN(su_new)){
			su_new = find(inter.page[inter.pa].sub,"name",su_new);
			if(su_new == undefined){ alertp("Problem changing page2"); return;}
		}
	
		if(false){
			let max = inter.page[inter.pa].sub.length;
			if(su_new >= max) su_new = max-1;
		}
		
		inter.page[inter.pa].index = su_new;
	}
	
	let susu_new = page_new.susu;
	if(susu_new != undefined){
		let pag = inter.page[inter.pa];
		
		if(isNaN(susu_new)){
			susu_new = find(pag.sub[pag.index].sub,"name",susu_new);
			if(susu_new == undefined){ alertp("Problem changing page3"); return;}
		}
		
		if(false){
			let max = pag.sub[pag.index].sub.length;
			if(susu_new >= max) susu_new = max-1;
		}
		
		pag.sub[pag.index].index = susu_new;
	}
	
	let sususu_new = page_new.sususu;
	if(sususu_new != undefined){
		let pag = inter.page[inter.pa];
		let p = pag.sub[pag.index].index;
		let susu = pag.sub[pag.index].sub[p];
		if(false){
			let max = susu.sub.length;
			if(sususu_new >= max) sususu_new = max-1;
		}
		
		susu.index = sususu_new;
	}
	
	let pag = inter.page[inter.pa];
	let na = pag.name;
	if(pag.sub.length > 0){
		let sub = pag.sub[pag.index];
		na += "->"+sub.name;
		if(sub.sub.length > 0){
			let subsub = sub.sub[sub.index];
			na += "->"+subsub.name;
		}
	}

	inter.page_name = na;

	inter.options = false;
	clear_comp_select();
	close_help();
	close_bubble();
	close_data_source();
	close_param_source();
	close_Amatrix_source();
	close_Xvector_source();
	close_view_graph();
	close_view_ind();
	inter.edit_source = false;
	mode_off();
	close_description();
	if(make_fig) inter.figure=[];
		
	if(inter.popup_timeout != undefined) clearTimeout(inter.popup_timeout);
	
	let tree = na.split("->");
	
	generate_screen();
	
	if(tree[0] == "Model" && tree[1] == "Compartments"){
		let p = model.get_p();
		if(model.species.length == 0){
			inter.popup_timeout = setTimeout(function(){ model.add_species_popup();},500);
		}
		else{
			if(model.species[p].ncla == 0){
				inter.popup_timeout = setTimeout(function(){ model.add_classification_popup();},500);
			}
			else{
				let cl = model.get_cl();
				if(model.species[p].cla[cl].ncomp == 0){
					inter.popup_timeout = setTimeout(function(){ model.add_compartment_popup();},500);
				}
			}
		}
	}
	
	model.check_show_warning();
}


/// Updates edit_source property
function update_edit_source(bu)
{
	if(bu.type == "BubbleEndBut"){
		if(inter.edit_source == true) inter.edit_source = "editted";
		if(inter.edit_source == "add") inter.edit_source = "adding";
	}
}


/// Clears any selected compartments
function clear_comp_select()
{
	inter.comp_select={list:[]};
}


/// Closes a dropdown menu if it is on
function close_dropdown()
{
	if(inter.dropdown.source != undefined){
		inter.dropdown = {};
		generate_screen();
	}
}


/// Closes the parameter source 
function close_param_source()
{
	inter.edit_param = {};
}


/// Fires when a key is pressed
function key_press(e)
{
	if(false) error(e.keyCode+" code");
	
	let code = e.keyCode;	
	//prr(code+"code");
	
	if(code == 223 && tab_name() == "Inference" && testing){
		start_worker("CopyValue");
	}
	
	if(make_fig){
		switch(code){
		case 8: // Deletes later letter
			if(inter.figure.length > 0){ inter.figure.pop(); generate_screen();}
			break;
			
		case 36:
			saving_dialogue("",".png","Export figure");
			break;
		}
		
		if(code >= 65 && code < 76){
			let letter = String.fromCharCode(code);
			inter.figure.push({te:letter, x:inter.mx-letter_size/2, y:inter.my-letter_size/2});
			generate_screen();
		}
	}
	
	if(test_comment){
		if(code == 38){ 
			inter.help.index--; 
			if(inter.help.index < 0) inter.help.index = 0;
		}
		if(code == 40){ 
			inter.help.index++; 
			if(inter.help.index >= all_comments_name.length) inter.help.index = all_comments_name.length-1;
		}
		generate_screen();
	}		

	if(code == 16 || code == 17) return;
	
	if(code == 27 && testing){ // Escape key
		console.clear();
		location.reload(true);
	}
	
	let cur = inter.cursor;
	if(code != 38 && code != 40) cur.xstore = undefined;

	if(cur.i != undefined && cur.l >= inter.layer.length){
		cur.i = undefined; cur.l = undefined;
	}

	let i = cur.i;

	if(i != undefined){	
		if(code == 9){                                // Tab key
			let lay_name = inter.bubble.lay_name;
			
			let l = cur.l;
			let ll = l+1; 
			while(ll < inter.layer.length && !(inter.layer[ll].name == "Input" && inter.layer[ll].op.source.hidden != true)) ll++;
		
			if(ll < inter.layer.length){
				unfocus_check();		
				
				cur.l = ll;
				cur.i = 0;
				cur.on =1;
				cur.click_pos = 0;
				cur.select_pos = undefined;

				reset_cursor();
				return;
			}
			else{
				if(bubble_on() == true && 
					(lay_name == "ParamValueContent" || lay_name == "ModelParamContent" || lay_name =="CreateEditParamContent" || lay_name == "CreateEditAmatrixContent" || lay_name == "InitialPopulationContent")){
					if(bubble_check_error() == false){
						let bub = inter.bubble;
					
						let i_now = inter.bubble.bu.ind;
					
						copy_back_to_source();
						
						if(bub.bu.type == "DistSplitElement" || bub.bu.type == "PriorSplitElement") update_prior_split();
				
						close_bubble();
					
						generate_screen(false);	

						let lay = get_lay(lay_name);
					
						if(lay != undefined){
							for(let i = i_now+1; i < lay.but.length; i++){
								let bu = lay.but[i];
								
								if(bu.ac != undefined){
									if(bu.type == "ParamSimElement" || bu.type == "DistSimElement" ||  bu.type == "PriorElement" ||
										bu.type == "DistElement" || bu.ac == "EditDerive" ||
										bu.type == "ParamElement" || bu.type == "ReparamElement" ||
										bu.type == "ParamElementConst" || 
										bu.type == "ParamFactorConst" || 
										bu.type == "ParamWeightConst" || 
										bu.type == "PriorSplitElement" || 
										bu.type == "DistSplitElement" || 
										bu.type == "CompPop"){
										activate_button(lay,i);
										return;
									}
								}
							}
						}
						return;
					}
				}
			}
		}

		let text_lay = inter.layer[cur.l];
		
		let st = text_lay.get_text_from_source();

		if(code == 37 || code == 38 || code == 39 || code == 40){
			if(e.shiftKey == true){ 	
				if(cur.select_pos == undefined) cur.select_pos = i;
			}
			else{
				cur.select_pos = undefined;
			}
			
			switch(code){
			case 37:                                     // Left arrow
				if(i > 0) i--; 
				break;
		
			case 39:                                     // Right arrow
				if(i < st.length) i++;
				break;

			case 38:                                     // Up arrow
				i = text_lay.arrow_up_down(i,"up");
				break;

			case 40:                                     // Down arrow
				i = text_lay.arrow_up_down(i,"down");
				break;
			}
			cur.i = i; if(cur.select_pos == i) cur.select_pos = undefined;
			text_lay.ensure_cursor_in_view();
			
			cur.on = 1;
		
			text_lay.plot_buttons();
			plot_screen();
		}
		else{
			let remember = { on:true, i:cur.i, select_pos:cur.select_pos};

			switch(code){
			case 8:                                      // Backspace
				if(cur.select_pos != undefined){ st = remove_selection(st); i = cur.i;}
				else{
					if(i > 0){
						st = st.substr(0,i-1)+st.substr(i);
						i--;
					}
				}
				break;

			case 46:                                     // Delete
				if(cur.select_pos != undefined){ st = remove_selection(st); i = cur.i;}
				else{
					if(i < st.length){
						st = st.substr(0,i)+st.substr(i+1);
					}
				}
				break;

			case 13:                                     // Enter
				if(text_lay.name == "TextBox"){
					if(inter.equation.te != undefined){
						equation_done();
						generate_screen();
						return;
					}
					else{
						st = st.substr(0,i)+"\n"+st.substr(i);
						i++;
					}
				}
				else{
					unfocus_check();		
				
					let flag = false;
					let l = inter.cursor.l;
					if(l != undefined){
						let sto = inter.textbox_store;
					}
				
					let fin_but = inter.bubble.final_button;
					if(fin_but != undefined && fin_but.ac != undefined && flag == false){	
						button_action(fin_but,"enter");  
						return;
					}
					else{
						turn_off_cursor();
						generate_screen();
					}
				}
				break;
			}

			if(e.key.length == 1 && e.ctrlKey == false){
				st = remove_selection(st); i = cur.i;
				let key = e.key;
				if(key == "*"){
					let fl = false;
					if(text_lay){
						switch(text_lay.op.source.type){
						case "find": case "replace": case "description": case "element_factor_const":
						case "element_weight_const": case "wild_card":
							fl = true; 
							break;				
						}
					}
					
					if(fl == false) key = "×";
				}
				if(key == "<") key = "〈";
				if(key == ">") key = "〉";
				st = st.substr(0,i)+key+st.substr(i);
				i++; 
				if(st.includes("\\")){
					let result = detect_greek(st,i);
					st = result.te; i = result.curi;
				}
			}

			if(e.keyCode == 90 && e.ctrlKey == true){
				let cur = inter.cursor;
				if(cur.i != undefined){
					let tbs = text_lay.get_text_box_store();
					if(tbs.store_text.length > 0){
						let store = tbs.store_text[tbs.store_text.length-1];
						st = store.text;
						i = store.i;
						cur.select_pos = store.select_pos;
						tbs.store_text.pop();
						remember.on = false;
					}
				}
			}
			
			text_lay.put_text_in_source(st,remember);
			
			cur.i = i; if(cur.select_pos == i) cur.select_pos = undefined;

			if(inter.equation.warning != undefined){
				delete inter.equation.warning.cur;
			}

			let lay = inter.layer[cur.l];
			lay.initialise_buttons();
			lay.plot_buttons();	
			if(lay.need_add_scroll()) generate_screen();
			plot_screen();
		
			inter.layer[cur.l].ensure_cursor_in_view();
		}
	}
}


/// Gets the mouse position
function mouse_position(canvas, evt)                         
{
	var rect = canvas.getBoundingClientRect();
	return {
		x: (evt.clientX - rect.left)/inter.sca,
		y: (evt.clientY - rect.top)/inter.sca
	};
}

	
/// Loads a data table into the data source
function load_data_table_initialise()
{
	inter.radio_heading = {value:"Yes"};
	inter.radio_format = {value:"csv"};
	
	inter.help = { title:"Data tables", te:"", show_datatable:true, upload:"UploadDataTable"};
}


/// Creates a new tables by transfering columns from a stored data table
function create_edit_table(lay)
{
	lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, type:"Rect", val:WHITE});
	
	let cx = corner.x;
	let cy = corner.y;
	
	let so = edit_source;
		
	let tab_source = data.table[so.data_table_use];
	let tab_dest = so.table

	let title = so.type, help;
	
	if(so.title) title = so.title;
	if(so.help) help = so.help;
		
	switch(title){
	case "Init. Pop.": 
		if(so.spec.radio_dist.value == "Dist"){
			title = "Initial population distribution";
			help = load_initpopdist_text;
		}
		break;

	case "LoadTensor":
		title = convert_tensor_text("Load tensor values",so.load_col.length-1); 
		help = convert_tensor_text(load_tensor_text2,so.load_col.length-1);
		break;
	}
	
	cy = lay.add_title(title,cx,cy,{title:title, te:help});
	
	let te;
	let ac;

	if(tab_dest.ncol >= so.load_col.length){
		if(tab_dest.edit == false){
			if(so.info.load_file){
				let type;
				switch(so.info.load_file){
				case "Comp": type = "compartment"; break;
				case "Trans": type = "transition"; break;
				case "Source": type = "source transition"; break;
				case "Sink": type = "sink transition"; break;
				case "LoadTensor": type = "constant"; break;
				case "LoadReparam": type = "reparameterisation"; break;
				case "LoadPriorSplit": type = "prior"; break;
				case "LoadDistSplit": type = "distribution"; break;
				default: error("Option not recognised 60"+so.info.load_file); break;
				}
				te = "This sets the "+type+" information. Click 'Next' to continue";
			}
			else{
				te = "This sets the data table used for analysis. Click 'Next' to continue."; 
			}
			ac = "TableEditOn";
		}
		else{
			let nrow;
			if(tab_source != undefined) nrow = tab_source.nrow;
			if(tab_dest != undefined) nrow = tab_dest.nrow;
			
			if(nrow != undefined && nrow > TABLE_ROW_MAX){
				te = "This table is too large to edit. Click 'Done' to complete."; 
			}
			else{
				te = "Edit this table (optional) and click 'Done' to complete."; 
			}
			
			switch(so.type){
			case "KnotTimes": ac = "AddKnotTimes"; break;
			case "CompMap": ac = "AddCompMap"; break;
			case "LoadTensor": ac = "LoadTensorDone"; break;
			case "LoadPriorSplit": ac = "LoadPriorSplitDone"; break;
			case "LoadDistSplit": ac = "LoadDistSplitDone"; break;
			case "LoadReparam": ac = "LoadReparamDone"; break;	
			case "Fixed Effect": ac = "LoadXvector2"; break;
			case "Init. Pop.": ac = "InitPopAdd"; break;
			case "APed": ac = "APedAdd"; break;
			default:
				if(so.info.load_file) ac = "AddFile" + so.info.load_file; 
				else{		
					if(so.info.i != undefined)	ac = "PlaceSourceBack";
					else ac = "AddDataSource";
				}
				break;
			}
		}
	}
	else{
		let col = so.load_col[tab_dest.ncol];
		te = "Select column '<b>"+col.heading+"</b>' for "+col.desc+":";
		
		let teri, teri2;
	
		if(col.heading == "t"){
			teri = "(click ['here','AddConstantTime'] to add a constant)";
			teri2 = "(click here to add a constant)";
		}		
		else{
			if(col.type == "compartment" && model.species[col.p].cla[col.cl].ncomp == 1){
				let claa = model.species[col.p].cla[col.cl];
				let name = claa.comp[0].name;
				teri = "(click ['here','AddConstantColumn|"+name+"'] to add '"+name+"')";
				teri2 = "(click here to add '"+name+"')";
			}
			else{
				switch(col.heading){
				case "ID": case "snp": case "Name": case "From": case "To": 
				case "x": case "y": case "Lat": case "Lng": break;
				default:
					teri = "(click ['here','AddConstantValue'] to add a constant value)";
					teri2 = "(click here to add a constant value)";
					break;
				}
			}
		}
		
		if(teri){
			let w = text_width(teri2,get_font(para_si));
			lay.add_paragraph(teri,w+1,lay.dx-w-3,cy,BLACK,para_si,para_lh);
		}
	}
	
	cy = lay.add_paragraph(te,lay.dx-4,cx,cy,BLACK,para_si,para_lh);
	cy += 0.4;
	
	add_layer("CreateEditTableContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-3.5,{});

	if(tab_dest.ncol == 0 || (tab_dest.ncol == 1 && so.type == "CompMap")){
		lay.add_corner_button([["Cancel","Grey","CancelEditSource"],["Next","GreyNext",ac]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	}		
	else{
		if(tab_dest.edit == false){
			lay.add_corner_button([["Back","GreyBack","BackEditSource"],["Cancel","Grey","CancelEditSource"],["Next","GreyNext",ac]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
		}
		else{
			lay.add_corner_button([["Cancel","Grey","CancelEditSource"],["Done","Grey",ac]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
		}
	}
}


/// Changes text such that tensor becomes vector or matrix
function convert_tensor_text(te,dim)
{
	switch(dim){
	case 1: 
		te = te.replace(/tensor/g,"vector"); 
		te = te.replace(/Tensor/g,"Vector"); 
		te = te.replace(/columns with headings given by the indices/g,"a column with heading given by the index"); 
		break;
	case 2: 
		te = te.replace(/tensor/g,"matrix"); 
		te = te.replace(/Tensor/g,"Matrix"); 
		break;
	}
	return te;
}


/// Displays the creation of a new data table					
function add_create_edit_table_buts(lay)
{
	let tab_source = data.table[edit_source.data_table_use];
	
	let tab_dest = edit_source.table;

	lay.op.table = edit_source.table;

	let x = 0;
	
	x = draw_data_table(lay,x,0,tab_dest,{name:"Data table", edit:tab_dest.edit});
	x += 1;
	
	if(tab_dest.ncol < edit_source.load_col.length){
		lay.add_button({x:x, y:1.3, dx:2, dy:1.3, type:"LeftArrow", col:BLACK});
		x += 2;
		
		draw_data_table(lay,x,0,tab_source,{name:"Loaded table", col_used:tab_dest.col_used, ac:"SelectColumn"});
	}
}


/// Draws a specified data table at a given location
function draw_data_table(lay,x,y,tab,op)
{
	let nrow = tab.nrow;
	let too_big = false;
	if(nrow > TABLE_ROW_MAX){ nrow = TABLE_ROW_MAX; too_big = true;}
	let nr = (nrow+1);
	
	let pad = 1;
	let mar = 0.4;
	let name = op.name;

	let x_start = x;
	let si_num = 0.8;
	let fo_num = get_font(si_num);
	let wnum = text_width(String(nrow),fo_num);
	x += wnum+0.6;
	
	let wmin = 0;
	if(name != undefined){
		let fo = get_font(0.8,"bold");
		wmin = text_width(name,fo)+0.5;
		lay.add_button({te:name, x:x+0.3, y:y, dx:wmin, dy:1, type:"Text", si:0.8, font:fo, col:BLACK});
		y += 1;
	}
	
	let fo_head = get_font(si_table,"bold");
	let fo = get_font(si_table);

	let cx = x+mar;
	
	for(let c = 0; c < tab.ncol; c++){
		let flag = false;
		if(op.col_used != undefined){
			for(let i = 0; i < op.col_used.length; i++) if(c == op.col_used[i]) flag = true;
		}
		
		if(flag == false){
			let wmax = 4;
			let head = tab.heading[c];
			let w = text_width(head,fo_head);
			if(w > wmax) wmax = w;
			
			// Works out if a compartment or index
			let comp_col = false;
			for(let p = 0; p < model.species.length; p++){
				let sp = model.species[p];
				for(let cl = 0; cl < sp.ncla; cl++){
					let claa = sp.cla[cl];
					if(head == claa.name || remove_prime(head) == claa.index) comp_col = true;
				}
			}
			
			let lenmax = 0, longest="";
			for(let r = 0; r < nrow; r++){
				let te = get_table_ele_text(tab,r,c);
				let l = te.length;
				if(l > lenmax){ lenmax = l; longest = te;}
			}
			
			w = 1.1*text_width(longest,fo)+1; if(w > 10) w = 10;
			if(w > wmax) wmax = w;
			if(wmax < 5) wmax = 5;
	
			wmax += pad;
			
			let cy = y+mar;

			if(op.edit == true && too_big == false){
				lay.add_button({te:head, x:cx, y:cy, dx:wmax, dy:dy_table, type:"Element", c:c, ac:"EditTableHead", head:true, font:get_font(si_table,"bold")}); cy += dy_table;
				for(let r = 0; r < nrow; r++){
					let te = get_table_ele_text(tab,r,c); 
					lay.add_button({te:te, x:cx, y:cy, dx:wmax, comp_col:comp_col, dy:dy_table, type:"Element", c:c, r:r, ac:"EditTableElement", font:get_font(si_table)}); cy += dy_table;
				}
			}
			else{
				let column=[];
				column.push(head);
				for(let r = 0; r < nrow; r++){
					let te = get_table_ele_text(tab,r,c); 
					column.push(te);
				}
				
				lay.add_button({te:column, x:cx, y:cy, dx:wmax, dy:nr*dy_table, type:"Column", space:dy_table, c:c, ac:op.ac});
			}
				
			if(too_big){
				let si_mar = 1;
				let fo_mar = get_font(si_mar);
				let mar_col = EDIT_MARGIN_COL;
				lay.add_button({te:"⋮", x:cx, y:cy+nr*dy_table, dx:2, dy:dy_table_param, type:"RightText", si:si_mar, font:fo_mar, col:mar_col});
			}
			cx += wmax;
		}
	}
	
	if(cx == x+mar){
		let column=[];
		let te = "Empty"; nr = 1;
		column.push(te);
		let dx = text_width(te,fo_head)+1.5;
		lay.add_button({te:column, x:cx, y:y+mar, dx:dx, dy:column.length*dy_table, type:"Column", space:dy_table});
		cx += dx;
	}
	
	let ww = cx+mar-x;
	for(let r = 0; r < nr-1; r++){
		let ac; if(op.edit == true && too_big == false) ac = "RowNumber"
		lay.add_button({te:String(r+1), x:x_start, y:y+mar+(1+r)*dy_table, dx:x-x_start-0.1, dy:dy_table, type:"RowNumber", r:r, si:si_num, font:fo_num, col:RED, ac:ac});
	}
	let dyout = nr*dy_table+2*mar; if(too_big) dyout += dy_table_param;
	
	lay.add_button({x:x, y:y, dx:ww, dy:dyout, type:"Outline", col:BLACK});

	lay.add_button({x:x, y:y+nr*dy_table+2*mar, dx:ww+0.1, dy:0.1, type:"Nothing"});

	return cx+mar;
}


/// Gets the text for a table element
function get_table_ele_text(tab,r,c)
{
	 let te = tab.ele[r][c];
	 if(typeof(te) == "object"){
		 return "# Polygons: "+te.polygon.length;
	 }		 
	 return te;
}


/// Converts from a text string to a filter
function find_filt(find)
{
	let filt = [];
	let spl = find.split("*");
	for(let i = 0; i < spl.length; i++){
		let len_after = 0;
		for(let j = i; j < spl.length; j++){
			len_after += spl[j].length;
		}
		filt.push({len:spl[i].length, te:spl[i], len_after:len_after});
	}
	
	return filt;
}


/// Works out if a string agees with a filter
function filt_agree(st,filt)
{
	st = String(st);
	if(filt.length == 1){
		if(st.toLowerCase() == filt[0].te.toLowerCase()) return true;
		return false;
	}
	return filt_agree2(st,0,filt,0);
}

function filt_agree2(st,i,filt,num)
{
	let fi = filt[num];
	if(i > st.length-fi.len_after) return false;
	
	let jmin, jmax;
	if(num == 0){ jmin = 0; jmax = 1;}
	else{
		if(num == filt.length-1){ jmin = st.length-fi.len; jmax = jmin+1;}
		else{
			jmin = i; jmax = st.length-fi.len_after;
		}
	}

	for(let j = jmin; j < jmax; j++){
		if(st.substr(j,fi.len).toLowerCase() == fi.te.toLowerCase()){
			if(num == filt.length-1) return true;
			else{
				if(filt_agree2(st,j+fi.len,filt,num+1) == true) return true;
			}
		}
	}
	return false;
}


/// Performs a search 
function done_search()
{
	let c = inter.bubble.bu.c;
	let tab = edit_source.table;
	
	let find = inter.bubble.find;
	
	let filt = find_filt(find);

	let row_find = [];
	for(let r = 0; r < tab.nrow; r++){
		if(filt_agree(tab.ele[r][c],filt)){
			row_find.push(r);
		}
	}
	
	if(row_find.length == 0){
		alert_help("Search problem","No search results were found");
	}
	else{
		if(row_find.length > 0){
			select_table_elelent(row_find[0],c);
		}
		inter.bubble.row_find = row_find;
		inter.bubble.search_select = 0;
		
		change_bubble_mode("SearchResult");
	}
}


/// Works out if a string agees with a filter and replaces
function filt_replace(st,filt,replace)
{
	st = String(st);
	if(filt.length == 1){
		if(st.toLowerCase() == filt[0].te.toLowerCase()) return replace[0];
		return st;
	}
	
	let res = filt_replace2(st,0,filt,0,replace);
	if(res == false) return st;
	return res;
}


/// Works out if a string agees with a filter and replaces
function filt_replace2(st,i,filt,num,replace)
{
	let fi = filt[num];
	if(i > st.length-fi.len_after) return false;
	
	let jmin, jmax;
	if(num == 0){ jmin = 0; jmax = 1;}
	else{
		if(num == filt.length-1){ jmin = st.length-fi.len; jmax = jmin+1;}
		else{
			jmin = i; jmax = st.length-fi.len_after;
		}
	}

	for(let j = jmin; j < jmax; j++){
		if(st.substr(j,fi.len).toLowerCase() == fi.te.toLowerCase()){
			if(num == filt.length-1){
				let te = ""; if(jmin > i) te += st.substr(i,jmin-i);
				return te+replace[num];
			}
			else{
				let res = filt_replace2(st,j+fi.len,filt,num+1,replace);
				if(res != false){
					let te = ""; if(jmin > i) te += st.substr(i,jmin-i);
					return te+replace[num]+res;
				}
			}
		}
	}
	return false;
}


/// Performs a replace 
function done_replace()
{
	let c = inter.bubble.bu.c;
	let tab = edit_source.table;
	let num = 0;
	
	let filt = find_filt(inter.bubble.find);
	let replace = inter.bubble.replace.split("*");
	if(filt.length != replace.length){
		alert_help("Wildcard error", "The number of wildcards '*' must match between 'find' and 'replace'.");
	}
	else{
		for(let r = 0; r < tab.nrow; r++){
			let te = String(tab.ele[r][c]);
			let te_new = filt_replace(te,filt,replace);
			
			if(te != te_new){
				tab.ele[r][c] = te_new;
				num++;
			}
		}
		inter.bubble.num = num;
		change_bubble_mode("ReplaceResult");
	}
}


/// Delete rows based on finding a string
function done_delete_rows()
{
	let c = inter.bubble.bu.c;
	let tab = edit_source.table;
	
	let find = inter.bubble.find;
	
	let filt = find_filt(find);

	let num = 0;
	let r = 0;
	while(r < tab.nrow){
		if(filt_agree(tab.ele[r][c],filt)){
			tab.ele.splice(r,1);
			tab.nrow--;
			num++;
		}
		else r++;
	}

	close_bubble();
	generate_screen();
	select_table_elelent(undefined,c);

	inter.bubble.num = num;
	change_bubble_mode("DeleteRowsResult");
}


/// Delete selected row
function done_delete_row()
{
	let r = inter.bubble.bu.r;
	let tab = edit_source.table;
	
	tab.ele.splice(r,1);
	tab.nrow--;
	close_bubble();
}

		
/// Goes back on a search result
function search_back()
{
	let bub = copy(inter.bubble);
	bub.search_select--;
	if(bub.search_select < 0) bub.search_select = bub.row_find.length-1;
	select_table_elelent(bub.row_find[bub.search_select],inter.bubble.bu.c);
	inter.bubble.mode = bub.mode;
	inter.bubble.search_select = bub.search_select;
	inter.bubble.row_find = bub.row_find;
}


/// Goes next on a search result
function search_next()
{
	let bub = copy(inter.bubble);

	bub.search_select++;
	if(bub.search_select == bub.row_find.length) bub.search_select = 0;

	select_table_elelent(bub.row_find[bub.search_select],inter.bubble.bu.c);
	inter.bubble.mode = bub.mode;
	inter.bubble.search_select = bub.search_select;
	inter.bubble.row_find = bub.row_find;
}


/// Performs an orderinf of a column 
function done_order()
{
	let c = inter.bubble.bu.c;
	let tab = edit_source.table;
	
	let list=[];
	for(let r = 0; r < tab.nrow; r++){
		list.push({r:r, te:tab.ele[r][c]});
	}

	switch(inter.bubble.order_radio.value){
	case "A-Z": list.sort( function(a, b){ return sort_string(a.te,b.te,1)}); break;
	case "Z-A": list.sort( function(a, b){ return sort_string(a.te,b.te,-1)}); break;
	case "Low-High": list.sort( function(a, b){ return sort_number(a.te,b.te,1)}); break;
	case "High-Low": list.sort( function(a, b){ return sort_number(a.te,b.te,-1)}); break;
	default: error("Cannot find ordering"); break;
	}

	let ele_copy = copy(tab.ele);
	tab.ele = [];
	for(let r = 0; r < tab.nrow; r++){
		tab.ele.push(ele_copy[list[r].r]);
	}
	
	close_bubble();
}


/// Selects a specified table element
function select_table_elelent(r,c)
{
	close_bubble();
	
	let lay_name = "CreateEditTableContent";
	
	let l = find(inter.layer,"name",lay_name);
	if(l == undefined){ error("Cannot find layer1 "+lay_name); return;}
	
	let lay = inter.layer[l];
	let but = lay.but;
	
	let i = 0; while(i < but.length && !(but[i].r == r && but[i].c == c)) i++;
	if(i == but.length){ error("Button not found"); return;}
	
	inter.bubble = { lay_name:lay_name, bu:but[i], i:i, show_warning:false, warning:false, op:{}, find_focus:true};
	
	reset_text_box();
	
	shift_button_in_view(l,i);
}


/// Selects a specified table element
function select_param_element(index)
{
	close_bubble();
	
	let lay_name = "CreateEditParamContent";
		
	let l = 0; while(l < inter.layer.length && inter.layer[l].name != lay_name) l++;
	if(l == inter.layer.length){ error("Cannot find layer2 "+lay_name); return;}

	let lay = inter.layer[l];
	let but = lay.but;
	
	let i = 0; 
	while(i < but.length && (!(but[i].ac == "EditParamElement" 
	   || but[i].ac == "EditDistSplitElement") || equal_vec(but[i].pindex,index) == false)) i++;
	if(i == but.length){ error("Button not found"); return;}

	activate_button(lay,i);
	
	shift_button_in_view(l,i);
}


/// Makes sure that button i can be viewed
function shift_button_in_view(l,i)
{
	let lay = inter.layer[l];
	let but = lay.but;
	
	if(l < inter.layer.length-1 && inter.layer[l+1].name == "Yscroll"){	
		let but_sc = inter.layer[l+1].but[0];
		let pos = inter.scroll_position[but_sc.scroll_ref];
		let sh = pos.shift;
	
		let view_range = 5;
		if(but[i].y < sh+view_range){
			change_scroll(but[i].y-(sh+view_range),but_sc,"page_shift");
		}
		
		if(but[i].y+but[i].dy > sh+lay.dy-view_range){
			change_scroll((but[i].y+but[i].dy)-(sh+lay.dy-view_range),but_sc,"page_shift");
		}
	}
}


/// Transfers a column from a loaded data table to a data source 
function transfer_column(c,all_snp_flag)
{
	let so = edit_source;

	let tab_source = data.table[so.data_table_use];
	let tab_dest = so.table;
	
	let replace_space = false;
	if(so.load_col[tab_dest.ncol].type == "comptext") replace_space = true;
	
	tab_dest.filename = tab_source.filename;
	for(let r = 0; r < tab_source.nrow; r++){
		if(tab_dest.nrow == 0) tab_dest.ele[r]=[];
		
		let te = tab_source.ele[r][c];
		if(replace_space) te = te.replace(/ /g,"-");
		tab_dest.ele[r].push(te);
	}

	tab_dest.heading.push(so.load_col[tab_dest.ncol].heading);
	
	tab_dest.col_used.push(c);
	tab_dest.ncol++; 
	tab_dest.nrow = tab_source.nrow;

	if(so.type == "Genetic" && tab_dest.ncol == 3 && all_snp_flag != false){
		let list_extra = [];
		switch(so.spec.type_radio.value){
		case "snp":
			{
				let snp_root = so.spec.snp_root;

				for(let c = 0; c < tab_source.ncol; c++){
					if(find_in(tab_dest.col_used,c) == undefined){
						if(tab_source.heading[c].substr(0,snp_root.length) == snp_root){
							for(let r = 0; r < tab_source.nrow; r++){
								tab_dest.ele[r].push(tab_source.ele[r][c]);
							}
			
							let head = "SNP"+(list_extra.length+2);
							tab_dest.heading.push(head);
							list_extra.push(c);
							tab_dest.ncol++;

							so.load_col.push({heading:head, desc:"SNP column",type: "text"});
						}
					}
				}
			}
			break;
			
		case "matrix":
			{	
				for(let r = 0; r < tab_dest.ele.length; r++){
					let val = tab_dest.ele[r][2];
					let c = find_in(tab_source.heading,val);
					if(c != undefined){
						for(let r = 0; r < tab_source.nrow; r++){
							tab_dest.ele[r].push(tab_source.ele[r][c]);
						}
							
						tab_dest.heading.push(val);
						list_extra.push(c);
						tab_dest.ncol++;

						so.load_col.push({heading:val, desc:"SNP column",type: "text"});
					}
				}
			}
			break;
		}
		tab_dest.list_extra = list_extra;
	}
}


/// Transfers a column back				
function transfer_column_back()
{
	let so = edit_source;
	let tab_dest = so.table;

	if(tab_dest.list_extra != undefined){  // For genetic data removes extra columns
		for(let i = 0; i < tab_dest.list_extra.length; i++){
			let c = tab_dest.list_extra[i];
			for(let r = 0; r < tab_dest.nrow; r++) tab_dest.ele[r].pop();
	
			tab_dest.heading.pop();
			tab_dest.ncol--;
			so.load_col.pop();
		} 
		delete tab_dest.list_extra;
	}

	for(let r = 0; r < tab_dest.nrow; r++) tab_dest.ele[r].pop();
	tab_dest.ncol--;
	tab_dest.heading.pop();
	tab_dest.col_used.pop();
}


/// Adds a constant column to the table
function add_constant_column(val)
{
	let so = edit_source;
		
	let tab_source = data.table[so.data_table_use];
	tab_source.heading.push("Constant");
	for(let r = 0; r < tab_source.nrow; r++){
		tab_source.ele[r].push(val);
	}
	tab_source.ncol++;
	transfer_column(tab_source.ncol-1);
}


/// Adds the contents of a table to a layer
function add_table_content(lay,table,back_col,dy)
{
	let mar = 0.5;
	let cx = mar, cy = 0.0;
	let si = 1;
	
	let col = BLACK;
	if(back_col != WHITE) col = WHITE;
	
	
	if(typeof table == 'string'){
		lay.add_button({te:table, x:cx, y:cy, dx:lay.dx-5, dy:1.3, type:"Text", si:si, font:get_font(si), col:col}); 
		return;
	}
	
	let head = table.heading;
	let ncol = head.length;
	
	for(let c = 0; c < ncol; c++){
		let w = table.width[c];
		let fo = get_font(si,"bold");

		lay.add_button({te:head[c].name, x:cx, y:cy, dx:w, dy:1.3, type:"Text", back_col:back_col, si:si, font:fo, col:col}); 
		cx += w;
	}
	
	cy += dy;
	
	for(let r = 0; r < table.content.length; r++){
		let row = table.content[r];
		cx = mar;
		for(let c = 0; c < ncol; c++){
			let ele = row[c];
			let w = table.width[c];
			
			switch(ele.te){
			case "":
				break;

			case "Invalid":
				lay.add_button({te:ele.te, x:cx, y:cy+0.1, dx:6.5, dy:1.1, type:"RedInvalid", back_col:back_col, info:ele.info, ac:ele.ac}); 
				break;

			case "Edit":
				lay.add_button({te:ele.te, x:cx, y:cy+0.1, dx:w, dy:1.2, type:"GreyEdit", back_col:back_col, info:ele.info, ac:ele.ac}); 
				break;
				
			case "View":
				lay.add_button({te:ele.te, x:cx, y:cy+0.1, dx:w-0.5, dy:1.2, type:"GreyView", back_col:back_col, info:ele.info, ac:ele.ac}); 
				break;
			
			case "Source Name":
				lay.add_button({te:ele.name, x:cx, y:cy+0.1, dx:w-0.5, dy:1.2, type:"NameLink", ac:ele.ac, si:si, font:get_font(si), so:ele.so}); 
				break;
				
			case "Delete":
				{
					let si = 1.3;
					lay.add_button({te:ele.te, x:cx+w-si/2, y:cy, dx:si, dy:si, type:"Delete", back_col:back_col, error:ele.error, siminf:ele.siminf, source:ele.source, p:ele.p, i:ele.i, ac:ele.ac}); 
				}
				break;
			
			default:
				{
					let fo = get_font(si);
				
					if(ele.ac != undefined){
						let ww = text_width(ele.te,fo)+0.1;
						if(ww > w) ww = w;
						lay.add_button({te:ele.te, x:cx, y:cy, dx:ww, dy:1.3, type:"Link", i:ele.i, ac:ele.ac, si:si, font:fo});
					}
					else{
						if(ele.pname == true){
							lay.add_button({te:ele.te, x:cx, y:cy, dx:w, dy:1.3, type:"ParamText", back_col:back_col, si:si, font:fo, col:col});
						}
						else{
							let comp_col;
							if(head[c].name == "Details") comp_col = true; 
								
							lay.add_button({te:ele.te, comp_col:comp_col, x:cx, y:cy, dx:w, dy:1.3, type:"Text", back_col:back_col, si:si, font:fo, col:col});
						}
					}					
				}
				break;	
			}
			
			cx += w;
		}
		cy += dy;
	}
	
	cx = mar+table.width[0];
	for(let c = 1; c < ncol; c++){
		let d = 0.4;
		
		if(head[c-1].name != "" && head[c].name != "" && head[c-1].name != "Spec." && head[c].name != "Spec."){
			lay.add_button({x:cx-d, y:0, dx:d, dy:cy, col:c, width:table.width, type:"TableSlider", ac:"TableSlider"});
		}
		cx += table.width[c];		
	}
	
}	


/// Drags the tableslider (assuming a mouse move even)
function drag_table_slider(x)
{
	let wmin = 2;
	let mod = inter.mode;
	let dx = x-mod.x_st;
	let c = mod.col;
	let w1 = mod.width_store[c-1]+dx;
	let w2 = mod.width_store[c]-dx;
	if(w1 < wmin) dx += wmin-w1;
	if(w2 < wmin) dx -= wmin-w2; 
	w1 = mod.width_store[c-1]+dx;
	w2 = mod.width_store[c]-dx;
	
	inter.mode.width[c-1] = w1;
	inter.mode.width[c] = w2;
	generate_screen();
}


/// Drags the right menu slider 
function drag_rightmenu_slider(x)
{
	let mod = inter.mode;
	
	right_menu_width = mod.width_store-(x-mod.x_st);
	if(right_menu_width < 7) right_menu_width = 7;
	if(right_menu_width > 30) right_menu_width = 30;
	
	generate_screen();
}


/// Rounds a number to reduce the number of decimal places
function round(num)
{
	let n = Number(num);
	let i;
	for(i = 0; i < 10; i++){
		let dif = n- Math.round(n);
		if(dif > -TINY && dif < TINY) break;	
		n *= 10; 
	}
	return num.toFixed(i);
}


/// Gets an element in DOM
function ById(a){
	return document.getElementById(a);
}     


/// Switches any mode being operater in off (e.g. placing a compartment)
function mode_off()
{
	inter.mode={};	
}


/// Prints a text error to the screen
function alertp(te)
{
	te = add_full_stop(te);
	
	inter.help = { title:"Sorry an error occured!", te:te};
	generate_screen();
}


/// If the result of a function throws an error then outputs to help
function output_help(res)
{
	if(res.err == true){
		inter.help = { title:"Sorry an error occured!", te:res.msg};
		generate_screen();
	}
}


/// Counts the number of a value in list of objects
function count(arr,prop,val)
{
	let num = 0;
	for(let i = 0; i < arr.length; i++){
		if(arr[i][prop] == val) num++;
	}
	return num;
}


/// Gets the possibilities from an array with a given property
function get_pos(ar,ob)
{
	let pos = [];
	for(let i = 0; i < ar.length; i++) pos.push({te:ar[i][ob], i:i});
	return pos;
}


/// Gets the width a dropdown menu need to have
function get_w_drop(pos)
{
	let fo = get_font(si_drop);
	let wmax = 0;
	if(pos != undefined){
		for(let i = 0; i < pos.length; i++){
			let w = text_width(pos[i].te,fo)+2.5;
			if(w > wmax) wmax = w;
		}
	}
	
	let w = text_width("Select",fo)+2.5;
	if(w > wmax) wmax = w;

	return wmax;
}


/// Views a warning
function view_warning(i)
{
	let warn = model.warn[i];

	model.warn_view = false;

	//error(warn.warn_type+" warn");
	switch(warn.warn_type){
	case "TransTreeInf":
		change_page({pa:"Model", su:"Compartments", susu:warn.p});
		press_button_prop("Menu","PageSubSub",["val"],warn.p);	
		break;
	
	case "ModelClass":
		change_page({pa:"Model", su:"Compartments", susu:warn.p, sususu:warn.cl});
		break;
		
	case "SimDetails":
		change_page({pa:"Simulation", su:"Run"});
		break;
		
	case "InfDetails":
		change_page({pa:"Inference", su:"Run"});
		break;
		
	case "InfFurtherDetails":
		change_page({pa:"Inference", su:"Run"});
		inter.options = true;
		break;
		
	case "Init_pop": 
		switch(inter.page[warn.pa].name){
		case "Simulation":
			change_page({pa:"Simulation", su:"Initial Conditions", susu:warn.p}); 
			break;
		case "Inference": 
			change_page({pa:"Inference", su:"Data", susu:warn.p});
			break;
		}
		
		press_button_prop("TableContent","GreyEdit",["info","i"],warn.i);	
		break;
		
	case "NoClaComp":
		change_page({pa:"Model", su:"Compartments", susu:warn.p, sususu:warn.cl});
		break;
		
	case "TransTreeComp":
		change_page({pa:"Model", su:"Compartments", susu:warn.p, sususu:warn.cl});
		break;
	
	case "MissingComp": case "MissingColour":
		select_bubble_compartment(warn.p,warn.cl,warn.c);
		break;
	
	case "MissingBox": 
		select_bubble_box(warn.p,warn.cl,warn.k);
		break;
		
	case "TransMistake":
		select_bubble_transition(warn.p,warn.cl,warn.i);
		break;
		
	case "MissingParamMultValue":
		change_page({pa:"Post. Simulation", su:"Parameter Mult."});
		break;
		
	case "MissingSimValue":
		change_page({pa:"Simulation", su:"Parameters"});
		press_button_prop("ParamValueContent","ParamSimElement",["name"],warn.name);
		break;
		
	case "ReparamValue":
		change_page({pa:"Model", su:"Parameters"});
		press_button_prop("ModelParamContent","ReparamElement",["name"],warn.name);
		break;
		
	case "MissingPriorValue": case "ErrorPriorValue":
		change_page({pa:"Inference", su:"Prior"});
		press_button_prop("ParamPriorContent","PriorElement",["name"],warn.name);
		break;
		
	case "MissingPriorSplitValue":
		change_page({pa:"Inference", su:"Prior"});
		press_button_prop("ParamPriorContent","ParamSimElement",["name"],warn.name);
		break;
		
	case "MissingDistValue":
		change_page({pa:"Model", su:"Parameters"});
		press_button_prop("ModelParamContent","DistElement",["name"],warn.name);
		break;
		
	case "reparam":
		change_page({pa:"Model", su:"Parameters"});
		press_button_prop("ModelParamContent","ReparamElement",["name"],warn.name);
		break;
	
	case "ReparamSquareSpline":
		change_page({pa:"Model", su:"Parameters"});
		break;
		
	case "SimPopulationProb":
		change_page({pa:"Simulation", su:"Initial Conditions", susu:warn.p});
		break;

	case "SourceProb": 
		change_page({pa:"Inference", su:"Data", susu:warn.p});
		press_button_prop("TableContent","GreyView",["info","i"],warn.ind);	
		break;

	case "SourceInvalid": 
		{
			let paname = "Inference", supage = "Data";
			if(warn.siminf == "sim"){ paname = "Simulation"; supage = "Population";}
			
			change_page({pa:paname, su:supage, susu:warn.p});
			press_button_prop("TableContent","RedInvalid",["info","i"],warn.ind);	
		}
		break;
		
	case "TransDataSpecProb":
		change_page({pa:"Inference", su:"Data", susu:warn.p});
		press_button_prop("TableContent","GreyView",["info","i"],warn.ind);	
		break;
		
	case "InfPopulationProb":
		change_page({pa:"Inference", su:"Data", susu:warn.p});
		break;
	
	case "Equation":
		{
			let info = warn.eqn_info;
	
			switch(warn.eqn_type){
			case "comp_prob": case "sim_comp_prob":
				select_bubble_data_element(info.p,info.i,info.r,info.c);
				break;
			
			case "reparam":
				select_reparam_element(info.par_name, info.index);
				break;
			
			case "Se": case "Sp":
				select_bubble_data_spec(info.p,info.i);
				break;
			
			case "trans_mean": case "trans_rate":
			case "trans_scale": case "trans_shape":
			case "trans_cv":
				select_bubble_transition(info.p,info.cl,info.i);;
				break;
			
			case "trans_bp": 
				select_bubble_transition(info.p,info.cl,info.i);
				break;
				
			case "derive_param": case "derive_eqn":
				select_bubble_derived(info.i);
				break;
			
			default: 
				error("This equation type is not recognised: "+warn.eqn_type);
				break;
			}
		}
		break;
	
	case "A matrix":
		change_page({pa:"Model", su:"Parameters"});
		break;
		
	case "X vector":
		change_page({pa:"Model", su:"Parameters"});
		break;
		
	case "DataProblem": 
		if(warn.siminf == "sim") change_page({pa:"Simulation", su:"Initial Conditions", susu:warn.p});
		else change_page({pa:"Inference", su:"Data", susu:warn.p});
		break;
	
	case "DeriveProblem":
		select_bubble_derived(warn.i);
		break;
		
	case "KnotProblem":
		change_page({pa:"Model", su:"Parameters"});
		break;
		
	case "ParamValueProblem":
		if(warn.par.type == "param factor"){
			change_page({pa:"Post. Simulation", su:"Parameter Mult."});
		}
		else{				
			switch(warn.par.variety){
			case "normal": change_page({pa:"Simulation", su:"Parameters"}); break;
			default: change_page({pa:"Model", su:"Parameters"}); break;
			}
		}
		break;
		
	case "KnotTimeProblem":
		if(warn.par.type == "param factor"){
			change_page({pa:"Post. Simulation", su:"Parameter Mult."});
		}
		else{				
			change_page({pa:"Model", su:"Parameters"}); 
		}	
		break;
	}
}


/// Gets the minimum width needed for a dropdown menu with a given set of posibilities
function get_drop_width(wmin,pos)
{
	let fo = get_font(si_drop);
	let w = 0;
	for(let i = 0; i < pos.length; i++){
		let ww = text_width(pos.te,fo);
		if(ww > w) w = ww;
	}
	w += 1;
	if(w < wmin) w = wmin;
	
	return w;
}


/// Removes a filter
function remove_filter(p,i,rpf)
{
	rpf.species[p].filter.splice(i,1);
	inter.graph.init = undefined;
	close_bubble();
}


/// Sets up the bubble for cloning a classification
function clone_class_init(p_sel)
{
	let pos=[];   // Works out all the classifications which could be cloned
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		if(p != p_sel){
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				if(find(model.species[p_sel].cla,"name",claa.name) == undefined){
					pos.push({te:sp.name+": "+claa.name, p:p, cl:cl});
				}
			}
		}
	}		

	let w = get_drop_width(12,pos);
	
	inter.bubble.drop = { te:select_str, p_sel:p_sel, w:w, pos:pos};

	change_bubble_mode("clone");
}


/// Displays a message at the center of the screen
function center_message(te,lay)
{
	let dx = 40;
	
	let si = para_si;
	
	let y = page_char_hei/2-lay.y;
	let spl = te.split("\n");
	for(let j = 0; j < spl.length; j++){			
		lay.add_button({te:spl[j], x:page_char_wid/2-lay.x-dx/2+menu_width/2, y:y, dx:dx, dy:si, type:"CenterText", font:get_font(si)}); 
		y += 1.5;
	}
}


/// Displays a message at the center of a box
function center_message_box(x,y,dx,dy,te,lay)
{
	let si = para_si;
	let ddy = 1.5;
	
	let spl = te.split("\n");
	y += dy/2 - 2 - (spl.length-1)*ddy/2;
	
	for(let j = 0; j < spl.length; j++){			
		lay.add_button({te:spl[j], x:x, y:y, dx:dx, dy:si, type:"CenterText", font:get_font(si)}); 
		y += 1.5;
	}
}


/// Gets the main tab page name
function tab_name()
{
	let tree = inter.page_name.split("->");
	return tree[0];
}


/// Gets the sub tab  name
function subtab_name()
{
	if(!inter.page_name) return "";
	let tree = inter.page_name.split("->");
	return tree[1];
}
	

/// Gets the sub sub tab  name
function subsubtab_name()
{
	let tree = inter.page_name.split("->");
	return tree[2];
}


/// Removes quotation marks from a string
function add_quote(te)
{
	if(te == "") te = "\"\"";
	else{
		if(te.substr(0,1) != "\"") te = "\""+te;
		if(te.substr(te.length-1,1) != "\"") te += "\"";
	}
	return te;
}


/// A message showing that a feature is currently not avaiable
function disactivated(te)
{
	if(te == undefined) te = ""; else te += " ";
	te += "This feature is currently under development, but unfortunately is unavailable yet.";
	
	alert_help("Feature not available",te);
}


/// Starts the worker doing a job
function start_worker(type,info)
{
	if(false) error("start work: "+type+" "+info);
	
	start_loading_symbol(0,type);
	
	generate_screen();
	
	inter.worker_mess.next = {type:type, info:info};
	
	if(false){ // This is used to diagnose which part of info is using a lot of memory
		if(info && info.plot_filter){
			for(let ele in info.plot_filter){
				let el = info.plot_filter[ele]
				let te = JSON.stringify(el);
				prr(ele+":"+te.length);		
			}
		}
	}
	
	if(subtab_name() != "Results") inter.worker_mess.next.model = model;

	setTimeout(worker_mess_process,10);
}


/// Processes any web worker messages 
function worker_mess_process()
{
	let wm = inter.worker_mess;
	if(wm.next && wm.active == false){
		wm.active = true;
		worker.postMessage(wm.next);
		delete wm.next;
	}
}


/// For some scroll bars make sure they go to the top
function scroll_to_top()
{
	let list = ["CreateEditTableContent"];
	for(let i = 0; i < inter.scroll_position.length; i++){
		let sp = inter.scroll_position[i];
	
		let j = 0; while(j < list.length && !sp.page.includes(list[j])) j++;
		if(j < list.length) sp.shift = 0;
	}
}


/// Stops the current message
function stop_worker()
{
	stop_loading_symbol();
	worker.postMessage({type:"Stop Load"});
	inter.worker_mess.active = false;
}


/// Determines a match between the two textstore objects
function find_match(ts,ts2)
{
	if(ts.length != ts2.length) return;
	let match = [];
	let flag = [];
	for(let i = 0; i < ts.length; i++){
		let type = ts[i].source.type;
		for(let j = 0; j < ts2.length; j++){
			if(flag[j] != true){
				if(ts2[j].source.type == type){
					flag[j] = true;
					match[i] = j;
					break;
				}
			}
		}
		if(match[i] == undefined) return;
	}
	return match;
}


/// Pastes content at the cursor location
function cursor_paste(paste)
{
	if(paste.length > 7 && paste.substr(0,7) == "OBJECT-"){
		if(bubble_on()){
			let bub_sum = get_bubble_summary();
			let bub_sum2 = JSON.parse(paste.substr(7));

			let match = find_match(bub_sum.ts,bub_sum2.ts);
			
			if(bub_sum.type != bub_sum2.type || 
				bub_sum.drop.length != bub_sum2.drop.length || 
				bub_sum.check.length != bub_sum2.check.length || 
				bub_sum.radio.length != bub_sum2.radio.length || 
				match == undefined){
				alert_help("Error pasting", "The copied object cannot be placed into this element.");
				generate_screen();
			}
			else{
				for(let i = 0; i < match.length; i++){
					inter.textbox_store[i].te = bub_sum2.ts[match[i]].te;
				}
				
				let lay = get_lay("Bubble");
				for(let i = 0; i < bub_sum.check.length; i++){
					let ch = bub_sum.check[i];
					lay.but[ch.i].source.check = ch.check;
				}
				
				for(let i = 0; i < bub_sum.radio.length; i++){
					let ra = bub_sum.radio[i];
					lay.but[ra.i].source.value = ra.value;
				}
				
				for(let i = 0; i < bub_sum.drop.length; i++){
					inter.layer[bub_sum.drop[i].l].op.source.te = bub_sum2.drop[i].te;
				}
				
				button_action(inter.bubble.final_button,"enter");  
			}
		}
	}
	else{
		paste = paste.replace(/</g,"〈");
		paste = paste.replace(/>/g,"〉");
		if(subtab_name() != "Description") paste = paste.replace(/\*/g,"×");
				
		let cur = inter.cursor;
		if(cur.i == undefined) return;

		let i = 0; 
		while(i < paste.length){
			let ch = paste.charCodeAt(i);
			if(ch == 13){
				if(inter.layer[cur.l].name == "TextBox"){
					paste = paste.substr(0,i)+"\n"+paste.substr(i+1);
				}
			}
			else{
				if(ch >= 0 && ch < 32){
					paste = paste.substr(0,i)+paste.substr(i+1);
					i--;
				}
			}
			i++;
		}

		let text_lay = inter.layer[cur.l];
		
		let te = text_lay.get_text_from_source();

		let remember = { on:true, i:cur.i, select_pos:cur.select_pos};

		te = remove_selection(te);
		
		te = te.substr(0,cur.i)+paste+te.substr(cur.i);
		text_lay.put_text_in_source(te,remember);
		cur.i += paste.length; 
		reset_cursor();
	}
}


/// Gets a summary of all the properties of a bubble (so it can be copied)
function get_bubble_summary()
{
	let bu = inter.bubble.bu;
	
	let sum = {type:bu.type, ts:inter.textbox_store, check:[], drop:[], radio:[]};
	
	let lay = get_lay("Bubble");
	for(let j = 0; j < lay.but.length; j++){
		let bu = lay.but[j];
		switch(bu.type){
		case "CheckboxButton":
			sum.check.push({ i:j, check:bu.source.check});
			break;
			
		case "RadioButton":
			sum.radio.push({ i:j, value:bu.source.value});
			break;	
		}
	}
	
	for(let l = lay.index+1; l < inter.layer.length; l++){
		let la = inter.layer[l];
		if(la.name == "Dropdown"){
			sum.drop.push({l:l, te:la.op.source.te});
		}
	}
	
	return sum;
}


/// Copies selected text to clipboard
function cursor_copy(cut)
{
	let te = get_selected_text();
	if(te.length == 0){
		if(bubble_on()){	
			let bub_sum = get_bubble_summary();
			
			te = "OBJECT-"+JSON.stringify(bub_sum);
		}
	}
	else{
		if(cut == true){
			let cur = inter.cursor;
			if(cur.l != undefined){
				let text_lay = inter.layer[cur.l];
			
				let st = text_lay.get_text_from_source();

				if(cur.select_pos != undefined){
					let remember = { on:true, i:cur.i, select_pos:cur.select_pos};
					
					st = remove_selection(st);
					text_lay.put_text_in_source(st,remember);
					reset_cursor();
				}
			}
		}
	
		let i = 0;
		while(i < te.length){
			if(te.substr(i,1) == "\n") te = te.substr(0,i)+"\n"+te.substr(i+1);
			i++;
		}
	}
	
	navigator.clipboard.writeText(te);
}


/// Specifies a selection box for selecting compartments
function select_box()
{
	let x1 = inter.mode.mx, y1 = inter.mode.my;
	let x2 = inter.mx, y2 = inter.my;
	if(x1 > x2){ let xx = x1; x1 = x2; x2 = xx;}
	if(y1 > y2){ let yy = y1; y1 = y2; y2 = yy;}

	let box = {x1:x1,y1:y1,x2:x2,y2:y2};
	if((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) > 0.5){ 
		box.on = true;
	}
	
	return box;
}

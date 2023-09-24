"use strict";

/// Change the page being viewed
function change_page(page_new)  
{
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
		
		inter.page[inter.pa].index = su_new;
	}
	
	let susu_new = page_new.susu;
	if(susu_new != undefined){
		let pag = inter.page[inter.pa];
		
		if(isNaN(susu_new)){
			susu_new = find(pag.sub[pag.index].sub,"name",susu_new);
			if(susu_new == undefined){ alertp("Problem changing page3"); return;}
		}
		
		pag.sub[pag.index].index = susu_new;
	}
	
	let sususu_new = page_new.sususu;
	if(sususu_new != undefined){
		let pag = inter.page[inter.pa];
		let p = pag.sub[pag.index].index;
		pag.sub[pag.index].sub[p].index = sususu_new;
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
	
	model.warn.length = 0;
	
	inter.options = false;
	
	close_help();
	close_bubble();
	close_data_source();
	close_param_source();
	close_Amatrix_source();
	close_Xvector_source();
	close_view_graph();
	
	mode_off();
	
	if(inter.popup_timeout != undefined) clearInterval(inter.popup_timeout);
	
	let tree = na.split("->");
	if((tree[0] == "Model" && tree[1] == "Parameters") ||
		 (tree[0] == "Simulation" && tree[1] == "Parameters") ||
		 (tree[0] == "Inference" && tree[1] == "Prior")){
		init_param();
	}
	
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
}


/// Closes a dropdown menu if it is on
function close_dropdown()
{
	if(inter.dropdown.source != undefined){
		inter.dropdown = {};
		generate_screen();
	}
}


/// Fires when a key is pressed
function key_press(e)
{
	//error(e.keyCode+" code");
	let code = e.keyCode;	

	if(code == 16 || code == 17) return;
	
	if(code == 27) location.reload(true); // Escape key
	
	let cur = inter.cursor;
	if(code != 38 && code != 40) cur.xstore = undefined;

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
						generate_screen();		
						
						let lay = get_lay(lay_name);
			
						for(let i = i_now+1; i < lay.but.length; i++){
							let bu = lay.but[i];
							
							if(bu.ac != undefined){
								if(bu.type == "ParamSimElement" || bu.type == "PriorElement" ||
									bu.type == "DistElement" || bu.ac == "EditDerive" ||
									bu.type == "ParamElement" || bu.type == "ReparamElement" ||
									bu.type == "ParamElementConst" || 
									bu.type == "PriorSplitElement" || 
									bu.type == "DistSplitElement" || 
									bu.type == "CompPop"){
									activate_button(lay,i);
									return;
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
			case 37:                             // Left arrow
				if(i > 0) i--; 
				break;
		
			case 39:                             // Right arrow
				if(i < st.length) i++;
				break;

			case 38:                             // Up arro
				i = text_lay.arrow_up_down(i,"up");
				break;

			case 40:                             // Down arrow
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
			case 8:                              // Backspace
				if(cur.select_pos != undefined){ st = remove_selection(st); i = cur.i;}
				else{
					if(i > 0){
						st = st.substr(0,i-1)+st.substr(i);
						i--;
					}
				}
				break;

			case 46:                             // Delete
				if(cur.select_pos != undefined){ st = remove_selection(st); i = cur.i;}
				else{
					if(i < st.length){
						st = st.substr(0,i)+st.substr(i+1);
					}
				}
				break;

			case 13:                             // Enter
				if(text_lay.name == "TextBox"){
					st = st.substr(0,i)+"\n"+st.substr(i);
					i++;
				}
				else{
					unfocus_check();		
				
					let flag = false;
					let l = inter.cursor.l;
					if(l != undefined){
						let sto = inter.textbox_store;
					}
				
					if(inter.bubble.final_button != undefined && flag == false){
						button_action(inter.bubble.final_button,"enter");  
						return;
					}
					else{
						turn_off_cursor();
					}
				}
				break;
			}

			if(e.key.length == 1 && e.ctrlKey == false){
				st = remove_selection(st); i = cur.i;
				let key = e.key;
				if(key == "*") key = "×";
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

			reset_cursor();
			
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
	
	let so = inter.edit_source;
		
	let tab_source = data.table[so.data_table_use];
	let tab_dest = so.table

	let title = so.type, help;
	switch(title){
	case "Init. Pop.": 
		title = "Initial population"; 
		help = load_initpop_text;
		break;
		
	case "Init. Pop. Prior":
		title = "Initial population prior";
		help = load_initpopprior_text;
		break;
		
	case "Add Ind.": 
		title = "Add individuals"; 
		help = load_add_ind_text;
		break;
		
	case "Remove Ind.":
		title = "Remove individuals";
		help = load_rem_ind_text;
		break;
		
	case "Move Ind.": 
		title = "Move individuals"; 
		help = load_move_ind_text;
		break;
	
	case "Comp File":
		title = "Load compartments"; 
		help = load_comp_text;
		break;
		
	case "Comp File Colour":
		title = "Load compartments"; 
		help = load_compcol_text;
		break;
		
	case "Comp File Pos": 
		title = "Load compartments"; 
		help = load_comppos_text;
		break;
		
	case "Comp File Pos Colour":
		title = "Load compartments"; 
		help = load_compposcol_text;
		break;
		
	case "CompMap":
		title = "Compartment map"; 
		help = load_compmap_text;
		break;
		
	case "LoadTensor":
		title = convert_tensor_text("Load tensor values",so.load_col.length-1); 
		help = convert_tensor_text(load_tensor_text2,so.load_col.length-1);
		break;
			
	case "LoadReparam":
		title = "Load reparameterisation";
		help = load_reparam_text2;
		break;
		
	case "LoadPriorSplit": 
		title = "Load priors"; 
		help = load_priorsplit_text2;
		break;
		
	case "LoadDistSplit": 
		title = "Load distributions"; 
		help = load_distsplit_text2;
		break;
	
	case "Population":
		title = "Load population data"; 
		help = load_population_text;
		break;
	
	case "Pop. Trans.":
		title = "Load population transition data"; 
		help = load_poptrans_text;
		break;
		
	case "Fixed Effect":
		title = "Load covariate vector for fixed effect"; 
		help = fixed_eff_text;
		break;
	
	case "Compartment":
		title = "Load compartmental data"; 
		help = load_compartment_text;
		break;
		
	case "Transition":
		title = "Load transition data"; 
		help = load_compartment_text;
		break;
		
	default:
		//error("Option error:"+title);
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
			te = "Edit this table (optional) and click 'Done' to complete."; 

			switch(so.type){
			case "CompMap": ac = "AddCompMap"; break;
			case "LoadTensor": ac = "LoadTensorDone"; break;
			case "LoadPriorSplit": ac = "LoadPriorSplitDone"; break;
			case "LoadDistSplit": ac = "LoadDistSplitDone"; break;
			case "LoadReparam": ac = "LoadReparamDone"; break;	
			case "Fixed Effect": ac = "LoadXvector2"; break;
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
		//te = "Select column '<b>"+col.heading.toLowerCase()+"</b>' for "+col.desc+":";
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
		//indices of the tensor
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
	let tab_source = data.table[inter.edit_source.data_table_use];
	let tab_dest = inter.edit_source.table

	let x = 0;
	
	x = draw_data_table(lay,x,0,tab_dest,{name:"Data table", edit:tab_dest.edit});
	x += 1;
	
	if(tab_dest.ncol < inter.edit_source.load_col.length){
		lay.add_button({x:x, y:1.3, dx:2, dy:1.3, type:"LeftArrow", col:BLACK});
		x += 2;
		
		draw_data_table(lay,x,0,tab_source,{name:"Loaded table", col_used:tab_dest.col_used, ac:"SelectColumn"});
	}
}


/// Draws a specified data table at a given location
function draw_data_table(lay,x,y,tab,op)
{
	let pad = 1;
	let mar = 0.4;
	let name = op.name;

	let x_start = x;
	let si_num = 0.8;
	let fo_num = get_font(si_num);
	let wnum = text_width(String(tab.nrow-1),fo_num);
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
			let wmax = 0;
			let head = tab.heading[c];
			let w = text_width(head,fo_head);
			if(w > wmax) wmax = w;
			
			for(let r = 0; r < tab.nrow; r++){
				let te = get_table_ele_text(tab,r,c); 
				let w = text_width(te,fo);
				if(w > wmax) wmax = w;
			}
			wmax += pad;
			if(wmax < 4) wmax = 4;
			if(c == tab.ncol-1 && cx+wmax < x+wmin) wmax = x+wmin-cx;
			
			let cy = y+mar;

			if(op.edit == true){
				lay.add_button({te:head, x:cx, y:cy, dx:wmax, dy:dy_table, type:"Element", c:c, ac:"EditTableHead", font:get_font(si_table,"bold")}); cy += dy_table;
				for(let r = 0; r < tab.nrow; r++){
					let te = get_table_ele_text(tab,r,c); 
					lay.add_button({te:te, x:cx, y:cy, dx:wmax, dy:dy_table, type:"Element", c:c, r:r, ac:"EditTableElement", font:get_font(si_table)}); cy += dy_table;
				}
			}
			else{
				let column=[];
				column.push(head);
				for(let r = 0; r < tab.nrow; r++){
					let te = get_table_ele_text(tab,r,c); 
					column.push(te);
				}
				
				lay.add_button({te:column, x:cx, y:cy, dx:wmax, dy:column.length*dy_table, type:"Column", space:dy_table, c:c, ac:op.ac});
			}
			
			cx += wmax;
		}
	}
	
	let nr = (tab.nrow+1);
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
		let ac; if(op.edit == true) ac = "RowNumber"
		lay.add_button({te:String(r+1), x:x_start, y:y+mar+(1+r)*dy_table, dx:x-x_start-0.1, dy:dy_table, type:"RowNumber", r:r, si:si_num, font:fo_num, col:RED, ac:ac});
	}
	lay.add_button({x:x, y:y, dx:ww, dy:nr*dy_table+2*mar, type:"Outline", col:BLACK});

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


/// Performs a search 
function done_search()
{
	let c = inter.bubble.bu.c;
	let tab = inter.edit_source.table;
	
	let find = inter.bubble.find;
	
	let row_find = [];
	for(let r = 0; r < tab.nrow; r++){
		if(tab.ele[r][c].toLowerCase() == find.toLowerCase()){
			row_find.push(r);
		}
	}
	
	if(row_find.length > 0){
		select_table_elelent(row_find[0],c);
	}
	inter.bubble.row_find = row_find;
	inter.bubble.search_select = 0;
	
	change_bubble_mode("SearchResult");
}


/// Delete rows based on finding a string
function done_delete_rows()
{
	let c = inter.bubble.bu.c;
	let tab = inter.edit_source.table;
	
	let find = inter.bubble.find;
	
	let num = 0;
	let r = 0;
	while(r < tab.nrow){
		if(tab.ele[r][c].toLowerCase() == find.toLowerCase()){
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
	let tab = inter.edit_source.table;
	
	tab.ele.splice(r,1);
	tab.nrow--;
	close_bubble();
}

		
/// Goes back on a search result
function search_back()
{
	let bub = inter.bubble;
	bub.search_select--;
	if(bub.search_select < 0) bub.search_select = bub.row_find.length-1;
	select_table_elelent(bub.row_find[bub.search_select],inter.bubble.bu.c);
}


/// Goes next on a search result
function search_next()
{
	let bub = inter.bubble;
	bub.search_select++;
	if(bub.search_select ==  bub.row_find.length) bub.search_select = 0;
	select_table_elelent(bub.row_find[bub.search_select],inter.bubble.bu.c);
}


/// Sorts string
function sort_string(a,b,sign)
{
	const A = String(a).toUpperCase(); 
  const B = String(b).toUpperCase();
  if(A < B) return -1*sign; 
	if(A > B) return 1*sign;
  return 0;
}


/// Sorts string
function sort_number(a,b,sign)
{
	const A = Number(a); 
  const B = Number(b);
	if(isNaN(A) || isNaN(B)) return 0;
  if(A < B) return -1*sign; 
	if(A > B) return 1*sign;
  return 0;
}



/// Performs an orderinf of a column 
function done_order()
{
	let c = inter.bubble.bu.c;
	let tab = inter.edit_source.table;
	
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
		
	let l = 0; while(l < inter.layer.length && inter.layer[l].name != lay_name) l++;
	if(l == inter.layer.length){
		error("Cannot find layer1 "+lay_name); return;
	}

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
	while(i < but.length && (!(but[i].ac == "EditParamElement" || but[i].ac == "EditDistSplitElement") || equal_vec(but[i].pindex,index) == false)) i++;
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
function transfer_column(c)
{
	let so = inter.edit_source;

	let tab_source = data.table[so.data_table_use];
	let tab_dest = so.table;
	
	tab_dest.filename = tab_source.filename;
	
	for(let r = 0; r < tab_source.nrow; r++){
		if(tab_dest.nrow == 0) tab_dest.ele[r]=[];
		tab_dest.ele[r].push(tab_source.ele[r][c]);
	}
	
	tab_dest.heading.push(so.load_col[tab_dest.ncol].heading);
	tab_dest.col_used.push(c);
	tab_dest.ncol++; 
	tab_dest.nrow = tab_source.nrow;

	if(so.type == "Genetic" && tab_dest.ncol == 3){
		let head = tab_source.heading[c];
		let i = 0; while(i < head.length && isNaN(head.substr(i,1))) i++;
		let root = head.substr(0,i);

		let list_extra = [];
		for(let c = 0; c < tab_source.ncol; c++){
			if(find_in(tab_dest.col_used,c) == undefined){
				if(tab_source.heading[c].substr(0,root.length) == root){
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
		tab_dest.list_extra = list_extra;
	}
}


/// Transfers a column back				
function transfer_column_back()
{
	let so = inter.edit_source;
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
	let so = inter.edit_source;
		
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
	
	/*
	if(outline != undefined){
		let sum = 0;
		for(let c = 0; c < ncol; c++) sum += head[c].width;
		lay.add_button({x:mar, y:0.5, dx:sum, dy:(table.content.length+1)*dy, type:"Outline", back_col:WHITE, col:WHITE}); 
	}
	*/
	
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
				lay.add_button({te:ele.te, x:cx, y:cy+0.1, dx:w, dy:1.1, type:"GreyEdit", back_col:back_col, info:ele.info, ac:ele.ac}); 
				break;
				
			case "View":
				lay.add_button({te:ele.te, x:cx, y:cy+0.1, dx:w-0.5, dy:1.1, type:"GreyView", back_col:back_col, info:ele.info, ac:ele.ac}); 
				break;
				
			case "Delete":
				{
					let si = 1.3;
					lay.add_button({te:ele.te, x:cx+w-si/2, y:cy, dx:si, dy:si, type:"Delete", back_col:back_col, source:ele.source, i:ele.i, ac:ele.ac}); 
				}
				break;
			
			default:
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
						lay.add_button({te:ele.te, x:cx, y:cy, dx:w, dy:1.3, type:"Text", back_col:back_col, si:si, font:fo, col:col});
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
		if(head[c-1].name != "" && head[c].name != ""){
			lay.add_button({x:cx-d, y:0, dx:d, dy:cy, col:c, width:table.width, type:"TableSlider", ac:"TableSlider"});
		}
		cx += table.width[c];		
	}
}	


/// Given a classification name return species
function get_p_cl_from_claname(name)
{
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let cl = 0; cl < sp.ncla; cl++){
			if(sp.cla[cl].name == name) return {p:p, cl:cl};
		}
	}
}

/*
/// Gets p and cl from classification name
function get_p_cl_from_name(te)
{
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let cl = 0; cl < sp.ncla; cl++){
			if(sp.cla[cl].name == te) return {p:p, cl:cl};
		}
	}
}
*/


/// Drags the tableslide (assuming a mouse move even)
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


/// Prints to the console
function pr(te)
{ 
	console.log(te);
}

/// Prints to the console
function comment(te)
{ 
	console.log(te);
}

/// Prints a text error to the console
function error(te)
{
	console.log("ERROR: "+te);
}


/// Prints a text error to the screen
function alertp(te)
{
	if(te.length > 0 && te.substr(te.length-1,1) != ".") te += ".";
	
	inter.help = { title:"Sorry an error occured!", te:te};
	generate_screen();
}


/// Copies an object
function copy(inp){ 
	return JSON.parse(JSON.stringify(inp));
}


/// Finds a value in list of objects
function find(arr,prop,val)
{
	let i = 0; while(i < arr.length && arr[i][prop] != val) i++;
	if(i == arr.length){
		//error("Cannot find array value. prop:"+prop+" val:"+val+"  arr:"+arr);
		return;
	}
	return i;
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


/// Finds a value in list of objects (irrespective of case)
function find_nocase(arr,prop,val)
{
	let i = 0; while(i < arr.length && arr[i][prop].toLowerCase() != val.toLowerCase()) i++;
	if(i == arr.length){
		//error("Cannot find array value. prop:"+prop+" val:"+val+"  arr:"+arr);
		return;
	}
	return i;
}


/// Finds a value in list of objects
function find_in(arr,val)
{
	let i = 0; while(i < arr.length && arr[i] != val) i++;
	if(i == arr.length)	return;
	return i;
}


/// Finds a value in list of string case-independent
function find_string_in(arr,val)
{
	let i = 0; while(i < arr.length && arr[i].toLowerCase() != val.toLowerCase()) i++;
	if(i == arr.length)	return;
	return i;
}


/// When typing this convert \alpha to greek alpha
function detect_greek(te,curi)
{
	let i = 0;
	while(i < te.length){
		while(i < te.length && te.substr(i,1) != "\\") i++;
	
		if(i < te.length){
			i++;
			for(let j = 0; j < greek_latex.length; j++){
				let wo = greek_latex[j][0];
				let len = wo.length;
				if(i <= te.length-len){
					if(te.substr(i,len) == wo){
						te = te.substr(0,i-1)+greek_latex[j][1]+te.substr(i+len);
						if(curi > i) curi -= len;
					}
				}
			}
		}
	}
	
	return {te:te,curi:curi};
}


/*
/// Finds a value in list of objects assuming lower case
function find_lc(arr,prop,val)
{
	let i = 0; 
	while(i < arr.length && arr[i][prop].toLowerCase() != val.toLowerCase()) i++;
	if(i == arr.length) return;
	return i;
}
*/

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
	let warn = copy(model.warn[i]);
	model.warn.length = 0;
	
	let info = warn.eqn_info;
			
	//error(warn.warn_type+" warn type"); error(warn);

	switch(warn.warn_type){
	case "ModelClass":
		change_page({pa:"Model", su:"Compartments", susu:warn.p, sususu:warn.cl});
		break;
		
	case "SimDetails":
		change_page({pa:"Simulation", su:"Run"});
		break;
		
	case "InfDetails":
		change_page({pa:"Inference", su:"Run"});
		break;
		
	case "Init_pop": 
		change_page({pa:warn.pa, su:"Population", susu:warn.p});
		press_button_prop("TableContent","GreyEdit",["info","i"],warn.i);	
		break;
		
	case "MissingComp": case "MissingColour":
		select_bubble_compartment(warn.p,warn.cl,warn.c);
		break;
	
	case "MissingSimValue":
		change_page({pa:"Simulation", su:"Parameters"});
		press_button_prop("ParamValueContent","ParamSimElement",["name"],warn.name);
		break;
		
	case "MissingPriorValue":
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
		press_button_prop("ModelParamContent","DistElement",["name"],warn.name);
		error("Need to do");
		break;
	
	case "SimPopulationProb":
		change_page({pa:"Simulation", su:"Population", susu:warn.p});
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
		switch(warn.eqn_type){
		case "comp_prob":
			select_bubble_data_element(info.p,info.i,info.r,info.c);
			break;
		
		case "reparam":
			select_reparam_element(info.par_name, info.index);
			break;
		
		case "trap_prob":	case "Se": case "Sp":
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
			error("This equation type is not recognised:"+warn.eqn_type);
			break;
		}
		break;
	
	case "A matrix":
		change_page({pa:"Model", su:"Parameters"});
		break;
		
	case "X vector":
		change_page({pa:"Model", su:"Parameters"});
		break;
		
	case "DataProblem": 
		if(warn.siminf == "sim") change_page({pa:"Simulation", su:"Population", susu:warn.p});
		else change_page({pa:"Inference", su:"Data", susu:warn.p});
		break;
			
	default:
		error("This warning type is not recognised:"+warn.warn_type);
		break;
	}
	
	model.warn.length = 0;
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


/// Makes sure the camera is up-to-date for any clones
function clone_camera(p,cl)
{
	let list = model.find_clones(p,cl);	
	if(list.length > 0){
		let cam = model.species[p].cla[cl].camera;
	
		for(let li of list){
			let cam2 = model.species[li.p].cla[li.cl].camera;
			cam2.x = cam.x; cam2.y = cam.y; cam2.scale = cam.scale; cam2.ruler = cam.ruler;
		}
		model.check_clones();
	}
}


/// Displays a message at the center of the screen
function center_message(te,lay)
{
	let dx = 40;
	
	let si = para_si;
	lay.add_button({te:te, x:page_char_wid/2-lay.x-dx/2+menu_width/2, y:page_char_hei/2-lay.y, dx:dx, dy:si, type:"CenterText", font:get_font(si)}); 
}


// Gets the main tab page name
function tab_name()
{
	let tree = inter.page_name.split("->");
	return tree[0];
}


/// Gets the sub tab  name
function subtab_name()
{
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
function remove_quote(te)
{
	if(te.substr(0,1) == "\"" && te.substr(te.length-1,1) == "\""){
		te = te.substr(1,te.length-2);
	}
	return te;
}


/// Removes <e> </e> equation quotes
function remove_eq_quote(te)
{
	if(typeof te == "string"){
		if(te.length >= 3 && te.substr(0,3) == "<e>") te = te.substr(3);
		if(te.length >= 4 && te.substr(te.length-4,4) == "</e>") te = te.substr(0,te.length-4);
	}
	return te;
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


/// Calculates if two numbers is different (subject to numerical noise) 
function dif(a,b)
{
	let dif = a-b; if(dif < 0) dif = -dif;
	if(dif > TINY) return true;
	return false;
}


/// A message showing that a feature is currently not avaiable
function disactivated(te)
{
	if(te == undefined) te = ""; else te += " ";
	te += "This feature is currently under development, but unfortunately is unavailable yet.";
	inter.help = {title:"Feature not available", te:te};
}
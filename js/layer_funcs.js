"use strict";

/// Generates all the layers on the screen
function generate_screen()
{
	let w = window.innerWidth;
	let h = window.innerHeight;
	let width = w, height = h;
	inter.canw = width, inter.canh = height;

	if(height < 538) height = 538; 
	if(width < 917) width = 917;
	
	inter.canvas.width = width;
  inter.canvas.height = height;
		
	
	let ratio = page_char_wid/page_char_hei;
	let wid_show = Math.round(height*ratio);
	let hei_show = height;
	
	if(wid_show > width){
		wid_show = width; hei_show = Math.round(wid_show/ratio);
	}

	inter.sca = wid_show/page_char_wid;

	inter.over = {};
	
	inter.layer = [];

	turn_off_cursor_flash();

	//check_error_textbox();

	reset_textbox_store();
	
	add_layer("Screen",0,0,width/inter.sca,height/inter.sca,{});
	
	purge_textbox_store();
	
	inter.over = get_button_over(inter.mx,inter.my);
	
	for(let l = 0; l < inter.layer.length; l++){
		inter.layer[l].plot_buttons() 
	}
	
	start_cursor_flash();
	
	plot_screen();
}


/// Determines if the mouse is within a specied layer
function within_layer(name)
{
	let lay = get_lay(name);
	
	if(inter.mx > lay.x && inter.mx < lay.x+lay.dx && inter.my > lay.y && inter.my < lay.y+lay.dy){
		return true;
	}		
	return false;
}


/// Prints the names of all the layers (for diagnostic purposes)
function print_layer()
{
	error("LAYER NAMES:");
	for(let l = 0; l < inter.layer.length; l++){
		error(inter.layer[l].name);
	}
}


/// Add a layer to the screen
function add_layer(na,x,y,dx,dy,op)
{
	x = nearest_pixel(x); y = nearest_pixel(y); 
	
	switch(na){
	case "Dropdown":
		dx = nearest_pixel(dx); dy = nearest_pixel(dy);
		break;
	}
	//dx = nearest_pixel(dx); dy = nearest_pixel(dy);
	
	inter.layer.push(new Layer(na,x,y,dx,dy,op));

	
	let lay = inter.layer[inter.layer.length-1];
	
	let opnew = {lay:lay, style:"Grey", position:"right"};
	let wid = scroll_width;

	switch(na){
	case "HelpContent": 
		if(inter.help.show_datatable == true || inter.help.script != undefined){
			wid = scroll_width*0.8;
		}
		else{
			opnew.style = "White"; wid = scroll_width*0.8;
		}
		break;
	case "EquationAddQuantity": opnew.style = "White"; wid = scroll_width*0.8; break;
	case "BubbleScrollable": opnew.style = "Bubble"; wid = scroll_width*0.8; break;
	case "Menu": opnew.style = "White"; wid = scroll_width*0.8; opnew.position = "left"; break;
	case "RightBotMenu": opnew.style = "TextBox"; wid = scroll_width*0.8; break;
	case "DropdownOptions": opnew.style = "TextBox"; wid = scroll_width*0.8; break;
	case "TextBox": opnew.style = "TextBox"; wid = scroll_width*0.8; break;
	case "UpperMenu": opnew.style = "TabSlider"; wid = scroll_width*0.7; break;
	case "Input": wid = 0; break;
	}

	lay.inner_dx = lay.dx;
	
	lay.initialise_buttons();

	lay.inner_dy = lay.dy;

	if(lay.invisible == true) return;
	
	switch(lay.name){
	case "Main": case "Compartment": case "Transition": 
	case "Annotation": case "AnnotationMap":
	case "GraphCompartments": case "GraphTransitions":
		let tree = inter.page_name.split("->");
		if(tree[0] == "Model" && tree[1] == "Compartments") na = "NoScroll";
		break;
	}

	switch(na){
	case "Yscroll": case "Xscroll": case "CornerScroll": case "NoScroll": case "Help": case "Menu": case "Expand":
		break;

	default:
		let box = get_but_box(lay.but);

		if(box.ymax > lay.dy+0.05){
			lay.y_scroll = true;
			if(opnew.position == "right") lay.inner_dx -= wid;
			lay.initialise_buttons();

			box = get_but_box(lay.but);

			if(box.xmax > lay.inner_dx){
				lay.x_scroll = true;
			
				lay.inner_dy -= wid;
				lay.initialise_buttons();
			}
		}
		else{
			if(box.xmax > lay.dx+TINY){
				lay.x_scroll = true;
			
				lay.inner_dy -= wid;
				lay.initialise_buttons();
			}
		}

		if(lay.y_scroll == true){
			let xpos = x+lay.dx-wid;

			//if(lay.name == "TextBox") xpos -= 0.2;
			//else{
				if(opnew.position != "right") xpos = x-wid;
			//}

			add_layer("Yscroll",xpos,y,wid,lay.dy,opnew);
		
/*
				if(lay.name == "TextBox") add_layer("Yscroll",x+lay.dx-wid-0.2,y,wid,lay.dy,opnew);
				else add_layer("Yscroll",x+lay.dx-wid,y,wid,lay.dy,opnew);
			}	
			else{
				add_layer("Yscroll",x-wid,y,wid,lay.dy,opnew);
			}
*/
		}

		if(lay.x_scroll == true) add_layer("Xscroll",x,y+dy-wid,dx,wid,opnew);
		
		if(lay.y_scroll == true && lay.x_scroll == true){
			add_layer("CornerScroll",x+dx-wid,y+dy-wid,wid,wid,opnew);
		}
	}

	//lay.plot_buttons();
}


/// Shifts a layer to a given location
function shift_layer(name,x,y)
{
	let l = find(inter.layer,"name",name);
	let lay_sc = inter.layer[l];
	let dx = x-lay_sc.x, dy = y-lay_sc.y;

	lay_sc.x += dx;
	lay_sc.y += dy;
	if(l < inter.layer.length-1){
		let lay = inter.layer[l+1];
		if(lay.name == "Yscroll"){
			lay.x += dx;
			lay.y += dy;
		}
	}
}


/// Scrolls the layer such that the last button can be viewed
function scroll_view_last_button(lay_name)
{
	let lay = get_lay(lay_name);

	let i = lay.but.length-1; while(i > 0 && lay.but[i].type == "TableSlider") i--;
	
	scroll_to_view(lay,lay.but[i]);
}


/// Shifts the scroll location of a layer such that a button can be observed
function scroll_to_view(lay,bu)
{
	let mar = 1;

	if(lay.x_scroll == true){
		let dx = 0;
		if(bu.x-lay.x_shift < mar) dx = bu.x-lay.x_shift - mar;
		if(bu.x+bu.dx-lay.x_shift > lay.inner_dx-mar){
			dx = bu.x+bu.dx-lay.x_shift - (lay.inner_dx-mar);
		}
		
		if(dx != 0){
			let l = lay.index;
			let ll = l+1;
			if(inter.layer[ll].name != "Xscroll") ll++;
			if(inter.layer[ll].name != "Xscroll") error("Cannot find scroll layer");
			else{
				change_scroll(dx,inter.layer[ll].but[0],"page_shift");
			}
		}
	}
	
	if(lay.y_scroll == true){
		let dy = 0;
		if(bu.y-lay.y_shift < mar) dy = bu.y-lay.y_shift - mar;
		if(bu.y+bu.dy-lay.y_shift > lay.inner_dy-mar){
			dy = bu.y+bu.dy-lay.y_shift - (lay.inner_dy-mar);
		}
		
		if(dy != 0){
			let l = lay.index;
			let ll = l+1;
			if(inter.layer[ll].name != "Yscroll") ll++;
			if(inter.layer[ll].name != "Yscroll") error("Cannot find scroll layer");
			else{
				change_scroll(dy,inter.layer[ll].but[0],"page_shift");
			}
		}
	}
}


/// Operates the up and down scroll buttons
function change_scroll(sign,bu,op)
{
	let pos = inter.scroll_position[bu.scroll_ref];

	let sh = pos.shift;
	
	switch(op){
	case "wheel":
		sh += 0.04*sign;
		break;
		
	case "page_shift":
		sh += sign;
		break;
		
	case "page_set":
		sh = sign;
		break;
		
	default:
		if(pos.discretise == true) sh += sign*textbox_linesi;
		else sh += 0.8*sign*pos.max*bu.frac;
		break;
	}
	
	if(sh < 0) sh = 0;
	if(sh > pos.max*(1-pos.frac)) sh = pos.max*(1-pos.frac);

	pos.shift = sh;
}


/// Removes selected text
function remove_selection(te)
{
	let cur = inter.cursor;
	if(cur.select_pos != undefined){
		if(cur.i > cur.select_pos){
			te = te.substr(0,cur.select_pos)+te.substr(cur.i)
			cur.i = cur.select_pos;
		}
		else{
			te = te.substr(0,cur.i)+te.substr(cur.select_pos);
		}
	}
	cur.select_pos = undefined;

	return te;
}


/// Actiaved when the mouse moves and the cursor is active
function drag_cursor(x,y)
{
	let cur = inter.cursor;
	let cur_over;
	
	if(cur.l != undefined){		
		let lay = inter.layer[cur.l];
		
		let over_new = get_button_over(x,y);
	
		let l = over_new.layer;
		
		if(l == cur.l){ // Selects based on being over a character
			let i = over_new.i;
		
			let lay = inter.layer[l];
			let bu = lay.but[i];
			
			cur_over = bu.i;
		}
		else{   // Selects based on mouse being outside area
			let sh_x = lay.x-lay.x_shift;
			let sh_y = lay.y-lay.y_shift;

			cur_over = 0;
			let dist = LARGE;
			let last;
			for(let i = 0; i < lay.but.length; i++){
				let bu = lay.but[i];
				if(bu.ac == "PositionCursor" && bu.ac){
					if(y >= sh_y+bu.y && y <= sh_y+bu.y+bu.dy){
						let d = (x-(sh_x+bu.x))*(x-(sh_x+bu.x));
						if(d < dist){
							dist = d; cur_over = bu.i;
						}
					}
					last = i;
				}
			}
			
			if(cur_over == 0 && last != undefined){	
				let bu = lay.but[last];
				if(y >= sh_y+bu.y+bu.dy) cur_over = bu.i;
			}
			
		
			
			//if(y < lay.y){ // Above upper edge
			//cur_over = 0;
			//	}
			//if(x >= lay.x && x <= lay.x+lay.dx && y >= lay.y && y <= lay.y+lay.dy){
			/*
			for(let i = lay.but.length-1; i >= 0; i--){
				let bu = lay.but[i];

				if(bu.ac && x >= sh_x+bu.x && x <= sh_x+bu.x+bu.dx && y >= sh_y+bu.y && y <= sh_y+bu.y+bu.dy){
			*/

			
		
		}

		if(cur_over != undefined){
			if(cur_over != cur.select_pos){
				cur.i = cur_over;
				cur.select_pos = cur.click_pos;
				lay.ensure_cursor_in_view();
				
				lay.plot_buttons();
				plot_screen();
			}
		}
	}
}		
				
				
/// Pastes content at the cursor location
function cursor_paste(paste)
{
	paste = paste.replace(/</g,"〈");
	paste = paste.replace(/>/g,"〉");
	paste = paste.replace(/\*/g,"×");
			
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


/// Copies selected text to clipboard
function cursor_copy(cut)
{
	let te = get_selected_text();

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
	
	navigator.clipboard.writeText(te);
}


/// Gets selected text
function get_selected_text()
{
	let cur = inter.cursor;
	if(cur.i == undefined) return "";

	let te = inter.layer[cur.l].get_text_from_source();
	if(cur.select_pos < cur.i){
		return te.substr(cur.select_pos, cur.i-cur.select_pos);
	}
	else{
		return te.substr(cur.i,cur.select_pos-cur.i);
	}
}	


/// When unfocusing cursor it checks to find errors in box
function unfocus_check()
{
	let l = inter.cursor.l;
	if(l != undefined){
		let lay = inter.layer[l];
		let source = lay.op.source;
			
		if(source != undefined && input_error(source) == true){
			let ref = lay.op.source.ref;
			let sto = inter.textbox_store;
	
			let i = 0; while(i < sto.length && sto[i].ref != ref) i++;
			if(i < sto.length){
				let warn_before = sto[i].warning;
				check_error_textbox2(sto[i]);
				if(sto[i].warning != warn_before) return true;
			}
		}
	}
	return false;
}	
		
		
/// Given a source works out if there is already a warning
function input_error(source)
{
	let k = find(inter.textbox_store,"ref",source.ref);
	if(k != undefined){
		let sto = inter.textbox_store[k];
		if(sto.warning != undefined && sto.warning != "") return true;
	}
		
	return false;
}

				
/// Switches off the cursor if it is on
function turn_off_cursor()
{
	turn_off_cursor_flash(); 
	inter.cursor = {};
}


/// Turns off the regular interval uded to initiate a cursor flash
function turn_off_cursor_flash()
{
	if(inter.cursor.interval != undefined){
		clearInterval(inter.cursor.interval);
		inter.cursor.interval = undefined;
	}
}


/// Starts the cursor flashing
function start_cursor_flash()
{
	if(inter.cursor.i == undefined) return;
	
	turn_off_cursor_flash();

	inter.cursor.interval = setInterval(
		function (){ 
			let cur = inter.cursor;
			if(cur.i != undefined){
				cur.on = 1-cur.on;
				inter.layer[cur.l].plot_buttons();
				plot_screen();
			}
		}, 500);
}


/// Resets the curson (so it is displayed)
function reset_cursor()
{
	let cur = inter.cursor;
	cur.on = 1;
	
	generate_screen();
}


/// Constructs the different layers on the screen
function add_screen_buts(lay)
{
	let x1 = nearest_pixel(lay.dx/2-page_char_wid/2);
	let x2 = x1 + nearest_pixel(menu_width);
	let x3 = nearest_pixel(lay.dx/2+page_char_wid/2);
	let y1 = nearest_pixel(0.5);
	
	let y2 = lay.dy-y1;
	let y3 = lay.dy-y1;
	let y4 = lay.dy-y1;

	let ya = 1.5;
	let yb = 5;
	let yc = y2-2;

	let show_model = model.get_show_model();
	
	let so = inter.edit_source;
	let type = so.type;

	if(inter.edit_param.type != undefined || inter.edit_Amatrix.i != undefined ||
		inter.edit_Xvector.i != undefined || inter.view_graph.type != undefined ||
		(type != undefined && so.table_loaded != false && so.edit_spec != true)){
	
		if(so.type != undefined && so.table_loaded == undefined){
			switch(so.type){
			case "Init. Pop.": add_layer("EditInitPop",x2,y1,x3-x2,y2-y1,{}); break;
			case "Init. Pop. Prior": add_layer("EditInitPopPrior",x2,y1,x3-x2,y2-y1,{}); break;
			default: error("Should not be here"); break;
			}
		}
		else{
			if(so.load_datatable == true && so.data_table_use != undefined){
				add_layer("EditTable",x2,y1,x3-x2,y2-y1,{});
			}	
			else{
				let type = inter.edit_param.type;
				if(type != undefined){
					add_layer("EditParam",x2,y1,x3-x2,y2-y1,{type:type});
				}
				
				if(inter.edit_Amatrix.i != undefined){
					add_layer("EditAmatrix",x2,y1,x3-x2,y2-y1,{});
				}
				
				if(inter.edit_Xvector.i != undefined){
					add_layer("EditXvector",x2,y1,x3-x2,y2-y1,{});
				}
				
				if(inter.view_graph.type != undefined){
					add_layer("ViewGraph",x2,y1,x3-x2,y2-y1,{});
				}
			}
		}
	}
	else{
		add_layer("Main",x2,y1,x3-x2,y2-y1,{});
	
		if(show_model == true){
			add_layer("AnnotationMap",x2,y1,x3-x2,y2-y1,{});
			add_layer("Annotation",x2,y1,x3-x2,y2-y1,{});
			add_layer("Compartment",x2,y1,x3-x2,y2-y1,{});
			add_layer("Transition",x2,y1,x3-x2,y2-y1,{});
		}
	
		if(inter.mode.type == "Add_Box"){
			add_layer("Select",x2,y1,x3-x2,y2-y1,{});
		}
		
		switch(inter.mode.type){
		case "Add_Compartment":
			{
				let w, h;
				let scale = model.get_scale();
				let coord = model.get_coord();
				let ruler = model.get_ruler();
				switch(coord){
				case "cartesian": w = 4*grid*scale; h = compartment_height*scale; break;
				case "latlng": w = 3*latlng_radius*scale*ruler; h = w; break;
				default: error("Option not recognised 67"); break;
				}
				add_layer("AddCompartment",1000,10,w,h,{coord:coord});
			}
			break;
		
		case "Add_Label":
			{
				let mo = inter.mode;
				mo.si = si_annotation*model.get_scale()*model.get_ruler();
				mo.fo = get_font(mo.si);
				mo.w = text_width(mo.te,mo.fo)+1;
				mo.h = mo.si;
				add_layer("AddLabel",1000,10,mo.w,mo.h,{});
			}
			break;
			
		}
		
		if(show_model == true){
			add_layer("LowerMenu",x2,y2-2.4,x3-x2,2.4,{});
			add_layer("UpperMenu",x2+1.2,y1,x3-x2-7,2,{});
		}
		
		//add_layer("LowerMenu",x2,y2-2.4,x3-x2,2.4,{});
	}
	
	let showbubble = false;
	
	let right_menu = "";
	let right_mid;
	
	let tree = inter.page_name.split("->");
	switch(tree[0]){
	case "Simulation": case "Inference":
		switch(tree[1]){
		case "Results":
			switch(tree[2]){
			case "Populations": case "Transitions": case "Individuals": case "Splines": case "Parameters":
				right_menu = tree[2];
				if(tree[0] == "Inference" && tree[2] == "Parameters"){
					let rpf = inf_result.plot_filter;
					
					switch(rpf.sel_paramview.te){
					case "Scatter": right_mid = "split"; break;
					case "Correlation": break;
					default:
						let pl = rpf.sel_paramviewtype.te;
						if(pl == "Trace" || pl =="Samples" || pl == "Distribution") right_mid = "with key";
						break;
					}
				}
				break;
			}
			break;
		}
		break;
	}
	
	if(inter.view_graph.type != undefined) right_menu = "GraphView";
	
	if(right_menu != ""){
		let xr = x3-right_menu_width;
		let yr = y1 + right_menu_top;
	
		let l = inter.layer.length;
		add_layer("RightMenu",xr,yr,x3-xr,y2-yr,{type:right_menu});
		
		let box = get_but_box(inter.layer[l].but);
		
		let ymax = box.ymax; if(ymax == -LARGE) ymax = 0; 
		
		let yr2 = yr+ymax+3;	
	
		if(right_mid){
			let hh = y2-yr2-0.5;
			switch(right_mid){
			case "with key":
				add_layer("RightMidMenu",xr,yr2,x3-xr-0.5,0.63*hh,{});
				add_layer("RightBotMenu",xr,yr2+0.65*hh,x3-xr-0.5,0.35*hh,{});
				break;
				
			case "split":
				add_layer("RightMidMenu",xr,yr2,x3-xr-0.5,0.48*hh,{});
				add_layer("RightBotMenu",xr,yr2+0.50*hh,x3-xr-0.5,0.48*hh,{});
				break;
			}
		}
		else{
			add_layer("RightBotMenu",xr,yr2,x3-xr-0.5,y2-yr2-0.5,{});
		}
	}

	/*
	if(inter.bubble.i != undefined && inter.bubble.l < inter.layer.length){
		add_layer("BubbleBack",x2,y1,x3-x2,y2-y1,{});
		add_layer("Bubble",10,10,10,10,{});
		showbubble = true;
	}
	*/
	
	add_layer("Frame",0,0,lay.dx,lay.dy,{x:x2, y:y1, dx:x3-x2, dy:y2-y1});

	add_layer("Menu",x1,yb,x2-x1,yc-yb,{});

	if(inter.bubble.i != undefined && showbubble == false){
		add_layer("BubbleBack",x2,y1,x3-x2,y2-y1,{});
		add_layer("Bubble",10,10,10,10,{});
	}
	
	add_layer("Logo",x1,ya,x2-x1,yb-ya,{});

	if(inter.equation.te != undefined){
		add_layer("EquationBackground",0,0,lay.dx,lay.dy,{});
		let wid = 0.54*(x3-x1);
		
		add_layer("Equation",lay.dx/2-wid/2+menu_width/2,0,wid,lay.dy,{});
	}
	
	add_dropdown_select();
	
	if(inter.help.title != undefined){
		add_layer("HelpBackground",0,0,lay.dx,lay.dy,{});
		let wid = 0.6*(x3-x1); 
		if(inter.help.scroll_to_line != undefined) wid = 0.8*(x3-x1);
		add_layer("Help",lay.dx/2-wid/2+menu_width/2,0,wid,lay.dy,{});
	}
	
	if(inter.loading_symbol.on){
		let six = 4.7, siy = 3.7;
		add_layer("LoadingSymbol",(x2+x3)/2-six/2,(y1+y2)/2-siy/2,six,siy,{});
	}
} 


// Places the selected dropdown layer on top of all the others
function add_dropdown_select()
{
	let source = inter.dropdown.source;
	let style = inter.dropdown.style;

	if(source != undefined){
		for(let l = 0; l < inter.layer.length; l++){
			let lay = inter.layer[l];
		
			if(lay.name == "Dropdown" && lay.op.source == source){
				let pos = lay.op.pos;

				let nl = pos.length;
				let nlmax = Math.floor((dropdown_dymax-dropdown_height)/dropdown_opheight);
				if(nlmax < 2) nlmax = 2;
				if(nl > nlmax) nl = nlmax;
		
				let dy = dropdown_height+dropdown_opheight*nl+0.7;
				dy = nearest_pixel(dy);
				
				if(lay.y+dy > page_char_hei -1){
					add_layer("Dropdown",lay.x,lay.y-dy+dropdown_height,lay.dx,dy,{source:source, pos:pos, direction:"Up", style:style});
					add_layer("DropdownOptions",lay.x+0.3,lay.y-dy+dropdown_height+0.3,lay.dx-0.6,dropdown_opheight*nlmax,{source:source, pos:pos});
				}
				else{
					add_layer("Dropdown",lay.x,lay.y,lay.dx,dy,{source:source, pos:pos, style:style});
					add_layer("DropdownOptions",lay.x+0.3,lay.y+dropdown_height+0.3,lay.dx-0.6,dropdown_opheight*nlmax,{source:source, pos:pos});
				}
				return;
			}
		}

		error("problem here");
	}
}


/// Starts the loading symbol
function start_loading_symbol(per)                                  
{
	if(!inter.loading_symbol.on){
		let interval = setInterval(function(){ inter.loading_symbol.offset++; replot_loading_symbol(); plot_screen();},110);
		inter.loading_symbol = {on:true, offset:0, interval:interval, processing:false, percent:per};
		generate_screen();
	}
}


/// Sets the loading percentage
function set_loading_percent(per)
{
	inter.loading_symbol.percent = per;
}


/// Stops the loading symbol
function stop_loading_symbol()                                   
{
	if(inter.loading_symbol.on){
		clearInterval(inter.loading_symbol.interval);
		inter.loading_symbol = {};
	}
}


/// Plots all the buttons on the current page
function plot_screen()                                    
{
	cv = inter.canvas_cv;
	
	clear_rectangle(0,0,inter.canw,inter.canh);   
	
	for(let i = 0; i < inter.layer.length; i++){
		let lay = inter.layer[i];
	
		if(lay.invisible != true && lay.can && lay.can.width > 0 && lay.can.height > 0){
			cv = inter.canvas_cv;

			if(lay.outline){	
				let out = lay.outline;
				if(out.style == "round"){
					let marx = out.marx, mary = out.mary;

					let extend = 0;
					if(lay.name == "Input" && lay.op.source.eqn == true) extend = expand_dx; 
					draw_round_rectangle(lay.x-marx,lay.y-mary,lay.dx+2*marx+extend,lay.dy+2*mary,0.5,WHITE,out.col);  
				}
				else{
					draw_rectangle(lay.x,lay.y,lay.dx,lay.dy,out.col,out.line);
				}
			}
		
			cv.drawImage(lay.can,Math.round(lay.x*inter.sca),Math.round(lay.y*inter.sca));
		}
	}
	
	if(inter.bubble.out_of_range == true){	
		close_bubble();
		generate_screen();
		return;  
	}
}


/// Gets the font based on size (and potenitally if bold ot italic)
function get_font(fac,te,type)
{
	if(type == undefined) type = "arial";
	let fo = Math.round(fac*inter.sca)+"px "+type;
	if(te) return te+" "+fo;
	return fo;
}


/// Returns a layer based on the layer name
function get_lay(te)
{
	for(let l = 0; l < inter.layer.length; l++){
		if(inter.layer[l].name == te) return inter.layer[l];
	}
	error("Cannot find layer3 "+te);
}


/// Finds a bounding box around all the button on a layer
function get_but_box(but)
{
	let xmin = LARGE, xmax = -LARGE;
	let ymin = LARGE, ymax = -LARGE;
	
	for(let i = 0; i < but.length; i++){
		let bu = but[i];
		if(bu.x < xmin) xmin = bu.x; if(bu.x+bu.dx > xmax) xmax = bu.x+bu.dx; 
		if(bu.y < ymin) ymin = bu.y; if(bu.y+bu.dy > ymax) ymax = bu.y+bu.dy; 
	}
	
	return {xmin:xmin, xmax:xmax, ymin:ymin, ymax:ymax};
}


/// Stores all the information in text boxes when they are created (such that editing the box does not directly effect the contents) 
function reset_text_box()
{
	inter.textbox_store = [];
}


/// Constructs a unique reference for a text box
function add_ref(source)
{
	let ref = source.type;
	if(source.val != undefined) ref += "_"+source.val;
	if(source.val2 != undefined) ref += "_"+source.val2;
	source.ref = ref;
}


function replot_loading_symbol()
{
	let lay = get_lay("LoadingSymbol");
	if(lay){
		lay.plot_button(lay.but[0]);
	}
}

/// Replots a given layer
function replot_layer(name)
{
	inter.replotting = true;
	
	let lay = get_lay(name);
	
	if(lay == undefined) return;

	lay.initialise_buttons();
	lay.plot_buttons();
}


/// Sets focus on 
function set_focus_first()
{
	let sto = inter.textbox_store;
	let num = 0; while(num < sto.length && sto[num].source.hidden == true) num++;
	if(num < sto.length) set_focus(num,"end");
}


/// Sets focus on a text box
function set_focus(num,i)
{
	if(i == undefined) i = 0;
	let sto = inter.textbox_store;
	if(sto.length == 0) return;
	
	if(sto.length <= num){ error("set focus problem"); return;}

	let ref = sto[num].ref;
	
	let l = 0;
	while(l < inter.layer.length && 
	   !((inter.layer[l].name == "Input" || inter.layer[l].name == "TextBox") && ref == inter.layer[l].op.source.ref)) l++;
	
	if(l == inter.layer.length){ error("Could not find layer"); return;}
	
	if(i == "end"){
		let te = inter.layer[l].get_text_from_source();
		i = te.length;
	}
	inter.cursor = {l:l, i:i, on:0, interval:undefined, click_pos:0, time:0};	
}
	

/// Resets textbox_store such that only those used are kept	
function reset_textbox_store()
{
	let st = inter.textbox_store;
	for(let i = 0; i < st.length; i++) st[i].used = false;
}


/// Purges those which have not been used
function purge_textbox_store()
{
	let st = inter.textbox_store;
	let i = 0; 
	while(i < st.length){
		if(st[i].used == false) st.splice(i,1);
		else i++;
	}
}

	
/// Copies the altered contents of text boxes to their source
function copy_back_to_source()
{
	let sto = inter.textbox_store;
	
	for(let i = 0; i < sto.length; i++){
		copy_back_to_source2(sto[i]);
	}
}

function copy_back_to_source2(tbs)
{
	let so = tbs.source;
	let te = tbs.te;
	
	switch(so.type){
	case "equation": inter.equation.te = te; break;
	case "description": model.description.te = te; break;
	case "compartment": model.rename_compartment(so.p,so.cl,so.i,te); break;
	case "add_species_name": inter.bubble.species_name = te; break;
	case "add_classification_name": inter.bubble.classification_name = te; break;
	case "classification_name": model.rename_classification(te,so.p,so.cl); break;
	case "species_name": model.rename_species(te,so.p); break;
	case "init_population": inter.edit_source.cla[so.cl].comp_init_pop[so.c].pop = Number(te); break;
		
	case "init_per":
		let cip = inter.edit_source.cla[so.cl].comp_init_pop;
		cip[so.c].pop_per = Number(te);
	
		let sum = 0;
		for(let c = 0; c < cip.length-1; c++) sum += cip[c].pop_per
		if(sum > 100 && sum < 100+TINY) sum = 100;
		cip[cip.length-1].pop_per = 100-sum;
		break;
		
	case "init_globpopulation": inter.edit_source.glob_comp[so.c].pop = Number(te); break;
	case "element": inter.edit_source.table.ele[so.r][so.c] = te;	break;
	
	case "element_param": 
		{
			let val = Number(te); if(isNaN(val)) val = te; 
			set_element(inter.edit_param.value,so.pindex,val);
		}
		break;
		
	case "element_eqn": 
		{
			let val = Number(te); if(isNaN(val)) val = te; 
			set_element(inter.edit_param.value,so.pindex,val);
		}
		break;
		
	case "reparam_eqn": 
		{
			let val = Number(te); if(isNaN(val)) val = te; 
			model.param[inter.bubble.th].value = val;
		}
		break;
		
	case "element_param_const": 
		{
			let val = Number(te);
			set_element(inter.edit_param.value,so.pindex,val);
		}
		break;
		
	case "element_Amatrix": inter.edit_Amatrix.value[so.j][so.i] = Number(te); break;	
	case "element_Xvector": inter.edit_Xvector.value[so.i] = Number(te); break;	
	case "find": inter.bubble.find = te; break;
	case "burnin": inter.bubble.burnin = Number(te); break;
	case "replace": inter.bubble.replace = te; break;
	case "time_start": inter.edit_source.spec.time_start = te; break;
	case "time_end": inter.edit_source.spec.time_end = te; break;
	case "pos_result": inter.edit_source.spec.pos_result = te; break;
	case "neg_result": inter.edit_source.spec.neg_result = te; break;
	case "percent": inter.edit_source.spec.percent = te; break;
	case "sd": inter.edit_source.spec.sd = te; break;	
	
	// equations
	case "trans_mean": 
		model.species[so.p].cla[so.cl].tra[so.i].value.mean_eqn.te = te;
		model.find_trans_pline(so.p,so.cl,so.i);
		break;
		
	case "trans_rate":
		model.species[so.p].cla[so.cl].tra[so.i].value.rate_eqn.te = te;
		model.find_trans_pline(so.p,so.cl,so.i);
		break;			
		
	case "trans_shape": 
		model.species[so.p].cla[so.cl].tra[so.i].value.shape_eqn.te = te; 
		model.find_trans_pline(so.p,so.cl,so.i);
		break;	
		
	case "trans_scale": 
		model.species[so.p].cla[so.cl].tra[so.i].value.scale_eqn.te = te;
		model.find_trans_pline(so.p,so.cl,so.i);
		break;
		
	case "trans_cv": 
		model.species[so.p].cla[so.cl].tra[so.i].value.cv_eqn.te = te; 
		model.find_trans_pline(so.p,so.cl,so.i);
		break;
		
	case "trans_bp": 
		model.species[so.p].cla[so.cl].tra[so.i].value.bp_eqn.te = te; 
		model.find_trans_pline(so.p,so.cl,so.i);
		break;

	case "trans_shape_erlang":
		model.species[so.p].cla[so.cl].tra[so.i].value.shape_erlang.te = te; 
		model.find_trans_pline(so.p,so.cl,so.i);
		break;

	case "Se": inter.edit_source.spec.Se_eqn.te = te; break;
	case "Sp": inter.edit_source.spec.Sp_eqn.te = te; break;
	case "trap_prob": inter.edit_source.spec.trap_prob_eqn.te = te; break;
	case "param_val":	
		model.param[so.val].value = te; 
		if(model.param[so.val].type == "trans_bp") branch_auto_calculate();
		break; 
	case "label": inter.bubble.label.te = te; break;
	case "label_anno": model.species[so.p].cla[so.cl].annotation[so.i].te = te; break;
	case "prior_min": case "prior_dist_min": 
		inter.bubble.prior.value.min_eqn.te = te; 
		break;
	case "prior_max": case "prior_dist_max": 
		inter.bubble.prior.value.max_eqn.te = te; 
		break;
	case "prior_mean": case "prior_dist_mean": 
		inter.bubble.prior.value.mean_eqn.te = te; 
		break;
	case "prior_shape": case "prior_dist_shape": 
		inter.bubble.prior.value.shape_eqn.te = te; 
		break;
	case "prior_sd": case "prior_dist_sd": 
		inter.bubble.prior.value.sd_eqn.te = te; 
		break;
	case "prior_cv": case "prior_dist_cv": 
		inter.bubble.prior.value.cv_eqn.te = te; 
		break;
	case "prior_alpha": case "prior_dist_alpha":
		inter.bubble.prior.value.alpha_eqn.te = te; 
		break;
	case "prior_beta": case "prior_dist_beta":
		inter.bubble.prior.value.beta_eqn.te = te; 
		break;
	case "derive_eqn1": inter.bubble.derived.eqn1.te = te; break;
	case "derive_eqn2": inter.bubble.derived.eqn2.te = te; break;
	case "derive_eqn": 
		model.derive[so.val].eqn2.te = te; 
		init_param();
		break;
	case "deriveparam_eqn": 
		model.derive[so.val].eqn1.te = te; 
		init_param();
		break;
	case "sim_t_start": model.sim_details.t_start = te; break;
	case "sim_t_end": model.sim_details.t_end = te; break;
	case "sim_number": model.sim_details.number = te; break;
	case "sim_timestep": model.sim_details.timestep = te; break;
		
	case "inf_t_start": model.inf_details.t_start = te; break;
	case "inf_t_end": model.inf_details.t_end = te; break;
	case "inf_timestep": model.inf_details.timestep = te; break;
	case "inf_sample": model.inf_details.sample = te; break;
	case "inf_chain": model.inf_details.nchain = te; break;
	case "inf_abcsample": model.inf_details.abcsample = te; break;
	case "inf_thinparam": model.inf_details.thinparam = te; break;
	case "inf_thinstate": model.inf_details.thinstate = te; break;
	case "inf_accfrac": model.inf_details.accfrac = te; break;
	case "inf_accfracsmc": model.inf_details.accfracsmc = te; break;
	case "inf_numgen": model.inf_details.numgen = te; break;
	case "inf_kernelsize": model.inf_details.kernelsize = te; break;
	case "inf_indmax": model.inf_details.indmax = te; break;
	case "sim_indmax": model.sim_details.indmax = te; break;


	case "knot_times": change_spline_knot(model.param[so.val],te.split(",")); break;
	case "smooth_value": model.param[so.i].spline.smooth.value = te; break; 
	case "min": inter.bubble.min = Number(te); break;
	case "max": inter.bubble.max = Number(te); break;
	case "fixed_time": inter.bubble.fixed_time = Number(te); break;
	case "const_column": inter.bubble.const_column = te; break;
	case "time_gen": inter.edit_source.time_gen = te; break;
	case "time_step": inter.edit_source.time_step = te; break;
	
	default: error("SOURCE PROBLEM: "+so.type); break;
	}
}


/// Checks that a string is a number
function check_number(te)
{
	let num = Number(te);
	if(isNaN(num)) return "Must be a number";
	return "";
}
				
/// Checks that a string is a non-negative number
function check_nonnegative(te)
{
	let num = Number(te);
	if(isNaN(num)) return "Must be a number";
	else{
		if(num < 0) return "Must be non-negative";
	}
	return "";
}


/// Checks that a string is a positive number
function check_posnumber(te)
{
	let num = Number(te);
	if(isNaN(num)) return "Must be a number";
	else{
		if(num <= 0) return "Must be positive";
	}
	return "";
}

/// Checks that a string is a number between zero and one
function check_zeroone(te)
{
	let num = Number(te);
	if(isNaN(num)) return "Must be a number";
	else{
		if(num < 0 || num > 1) return "Must be between 0 and 1 inclusive";
	}
	return "";
}


/// Checks that a string is a percentage
function check_percent(te)
{
	let num = Number(te);
	if(isNaN(num)) return "Must be a number";
	else{
		if(num <= 0) return "Number be positive";
		if(num > 100) return "Cannot be greater than 100%";
	}
	return "";
}


/// Checks that a string is either a single parameter or value in range 0-1
function check_param_or_number(te,type)
{
	let eqn = create_equation(te,"test");
	if(eqn.warn.length > 0) return eqn.warn[0].te;
	
	if(eqn.ind_eff.length > 0) return "Cannot contain individual effects"; 
	if(eqn.param.length > 1) return "Can only be constant or contain one parameter"; 
		
	if(!isNaN(te)){
		let val = Number(te);
		if(type == "zeroone"){
			if(val < 0 || val > 1) return "Value must be in the range 0 to 1";
		}
		
		if(type == "pos"){
			if(val <= 0) return "Value must be positive";
		}
	}
	return "";
}


/// Checks that a string is a positive integer
function check_posinteger(te)
{
	let num = Number(te);
	if(isNaN(num) || num <= 0 || num != Math.floor(num)){
		return "Must be a positive integer";
	}
	return "";
}


/// Checks for any error in text boxes
function check_error_textbox()
{
	let sto = inter.textbox_store;
	
	let flag = false;

	for(let i = 0; i < sto.length; i++){
		if(check_error_textbox2(sto[i]) == true) flag = true;
	}
	
	return flag;
}
	
function check_error_textbox2(tbs)
{
	let warn = "";
	
	let bub = inter.bubble;
	
	let source = tbs.source;

	if(source.hidden != true){
		let te = tbs.te;
		
		if(te == "") warn = "Cannot be empty!";
		
		let char_lim = 20;

		if(warn == ""){
			switch(source.type){ 			
			case "compartment":
				{
					let flag = false;
					
					let list = model.find_clones(source.p,source.cl);	
					
					for(let p = 0; p < model.species.length; p++){
						let sp = model.species[p];
						for(let cl = 0; cl < sp.ncla; cl++){
							for(let c = 0; c < sp.cla[cl].ncomp; c++){
								if(sp.cla[cl].comp[c].name == te){
									if(!(p == source.p && cl == source.cl && c == source.i)){
										let k = 0; 
										while(k < list.length && !(p == list[k].p && cl == list[k].cl && c == source.i)) k++;
									
										if(k == list.length){
											warn = "Name must be unique";
										}
									}
								}
							}
						}
					}
				}
				break;
				
			case "add_classification_name":
				{
					//let list = model.find_clones(bu.p,bu.cl);
					for(let p = 0; p < model.species.length; p++){
						let claa = model.species[p].cla;	
						let cl = 0; while(cl < claa.length && claa[cl].name.toLowerCase() != te.toLowerCase()) cl++;
						if(cl < claa.length) warn = "Classifications must be different";
					}
					char_lim = 20;
				}
				break;
				
			case "classification_name":
				{
					let list = model.find_clones(source.p,source.cl);	
					
					for(let p = 0; p < model.species.length; p++){
						let claa = model.species[p].cla;
						let cl = 0; 
						while(cl < claa.length && !(claa[cl].name.toLowerCase() == te.toLowerCase())) cl++;
						if(cl < claa.length){
							if(!(p == source.p && cl == source.cl)){
								let k = 0; 
								while(k < list.length && !(list[k].p == p && list[k].cl == cl)) k++; 
	
								if(k == list.length){
									warn = "Classifications must be different";
								}
							}
						}
					}
					char_lim = 20;
				}
				break;
				
			case "species_name":
				break;
			
			case "init_population": case "init_per": case "init_globpopulation":
				{
					warn = check_nonnegative(te);
				
					if(warn == "" && source.type == "init_per"){
						let cip = inter.edit_source.cla[source.cl].comp_init_pop
						let sum = 0;
						for(let c = 0; c < cip.length-1; c++){
							sum += cip[c].pop_per;
						}
						if(sum > 100+TINY) warn = "This is too large";
					}
				}
				break;
				
			case "element":
				warn = check_element(te,source.c);
				break;
			
			case "find":
				{
					let c = inter.bubble.bu.c;
					let tab = inter.edit_source.table;
					let r = 0;
					while(r < tab.nrow && String(tab.ele[r][c]).toLowerCase() != te.toLowerCase()) r++;
					if(r == tab.nrow) warn = "Cannot find";
				}
				break;
			
			case "replace":
				break;

			case "time_start": case "time_end":
				{
					let num = Number(te);
					if(isNaN(num)) warn = "Must be a number";
				}
				break;
				
			case "sim_t_start": case "sim_t_end": case "inf_t_start": case "inf_t_end":
				warn = check_number(te);
				break;
				
			case "sim_timestep": case "inf_timestep":
				warn = check_posnumber(te);
				break;

			case "inf_sample": case "inf_abcsample": 
				warn = check_posinteger(te);
				if(warn == "" && Number(te) < 100) warn = "Must be 100 or more";
				break;
			
			case "inf_chain":
				warn = check_posinteger(te);
				if(warn == "" && Number(te) > 10) warn = "Cannot be greater than 10";
				break;
				
			case "inf_thinparam": 
				warn = check_posinteger(te);
				break;
			
			case "inf_thinstate":
				warn = check_posinteger(te);
				break;
				
			case "inf_accfrac": case "inf_accfracsmc": 
				warn = check_zeroone(te)
				break;
		
			case "inf_numgen":
				warn = check_posinteger(te);
				break;
				
			case "inf_kernelsize":
				warn = check_posnumber(te);
				break;
				
			case "inf_indmax":
				warn = check_posinteger(te);
				break;
				
			case "sim_indmax":
				warn = check_posinteger(te);
				break;
			
			case "Se": case "Sp":
				warn = check_param_or_number(te,"zeroone");
				break;
				
			case "pos_result": case "neg_result":
				break;
				
			case "percent": 
				warn = check_percent(te);
				break;
				
			case "sd":
				warn = check_posnumber(te);
				break;
				
			case "label_anno":
				char_lim = 100; 
				break;
			
			case "param_val":
				warn = check_number(te);
				break;
			
			case "trap_prob":
				warn = check_param_or_number(te,"zeroone");
				break;

			case "element_param":
				break;

			case "element_eqn":
				break;
			
			case "reparam_eqn":
				break;
				
			case "element_param_const": 
				{
					let num = Number(te);
					if(isNaN(num)) warn = "Must be a number";
				}
				break;
				
			case "element_Amatrix":
				{	
					let num = Number(te);
					if(isNaN(num)) warn = "Must be a number";
				}
				break;
				
			case "trans_rate":
			case "trans_mean":
			case "trans_shape":
			case "trans_bp":
			case "trans_cv":
				//warn = check_valid_equation(te);
				break;
			
			case "trans_shape_erlang": 
				warn = check_posinteger(te);
				break;
			
			case "prior_min": case "prior_max": case "prior_mean":
				switch(bub.prior.type.te){
				case "exp": case "gamma": case "log-normal": warn = check_posnumber(te); break;
				case "normal": warn = check_number(te); break;
				case "uniform": warn = check_number(te); break;
				case "fix": warn = check_number(te); break;
				case "bernoulli": warn = check_zeroone(te); break;
				default: error("option problem"); break;
				}
				break;
			
			case "prior_alpha": warn = check_posnumber(te); break;

			case "prior_beta": warn = check_posnumber(te); break;
						
			case "prior_shape": case "prior_sd": case "prior_cv": 
				warn = check_posnumber(te);
				break;
				
			case "prior_dist_min": case "prior_dist_max": 
				warn = check_param_or_number(te);			
				break;
				
			case "prior_dist_mean":
				warn = check_param_or_number(te);			
				break;
				
			case "prior_dist_shape": case "prior_dist_sd": case "prior_dist_cv": 
				warn = check_param_or_number(te,"pos");		
				break;
					
			case "knot_times":
				warn = check_knot_times(te);
				char_lim = 1000;
				break;
				
			case "smooth_value":
				warn = check_posnumber(te);
				break;
		
			case "add_species_name":
				break;
				
			case "derive_eqn":
				//warn = check_valid_equation(te);
				break;
			
			case "derive_eqn1": case "derive_eqn2": break;
			case "derive_eqn": case "deriveparam_eqn": break;
			
			case "sim_number": 
				if(isNaN(te)) warn = "Must be a number";
				let num = Number(te);
				if(num <= 0 || num != Math.floor(num)) warn = "Must be a positive integer";
				break;
			
			case "min": case "max":
				if(isNaN(te)) warn = "Must be a number";
				break;
				
			case "fixed_time":
				if(isNaN(te)) warn = "Must be a number";
				break;
			
			case "const_column":
				break;
			
			case "element_Xvector":
				if(isNaN(te)) warn = "Must be a number";
				break;
				
			case "burnin": 
				if(isNaN(te)) warn = "Must be a number";
				else{
					let num = Number(te);
					if(num < 0 || num > 90)  warn = "Must be between 0% and 90%";
				}
				break;
			
			case "time_gen":
				{
					let t_start = Number(sim_result.details.t_start);
					let t_end = Number(sim_result.details.t_end);
	
					switch(inter.edit_source.time_radio.value){
					case "Periodic": case "Fixed time": 
						let t = Number(te);
						if(isNaN(te)) warn = "Must be a number"; 
						else{
							if(t < t_start || t > t_end) warn = "Must be in the simulation time range"; 
						}
						break;
					
					case "Specified":
						let spl = te.split(",");
						for(let j = 0; j < spl.length; j++){
							let t = Number(spl[j]);
							if(isNaN(spl[j])) warn = "Not valid";
							else{
								if(t < t_start || t > t_end) warn = "Must be in the simulation time range"; 
							}
						}
						break;
					}
				}
				break;
				
			case "time_step":
				warn = check_posnumber(te);	
				break;
			
			case "":
				break;
				
			default: error("Option not recognised 80"+source.type); break;
			}
		}
				
		if(source.eqn != true && te.length > char_lim && warn == ""){
			warn = "Not more than "+char_lim+" characters";
		}
		
		if(source.eqn == true && warn == ""){
			let eqn = tbs.eqn;
			
			eqn.te = te;
			
			extract_equation_properties(eqn);
			
			if(eqn.warn.length > 0){
				warn = eqn.warn[0].te;
			}
			else{
				if(source.type == "derive_eqn"){
					if(source.val != undefined){
						let eqn1 = model.derive[source.val].eqn1;
						if(eqn1 != undefined){
							let wa = check_derived_param(eqn,eqn1);
							if(wa != "success") warn = wa;
						}
					}
				}
				
				if(eqn.mode == "derive_param"){
					if(eqn.param.length == 0) warn = "Must contain a parameter";
					if(eqn.param.length > 1) warn = "Must contain just a single parameter";
					if(eqn.ind_eff.length > 0) warn = "Cannot contain individual effects";
				}		
				
				if(eqn.mode == "param_only"){
					if(eqn.ind_eff.length > 0) warn = "Cannot contain individual effects";
				}
			}
		}
	}
	
	tbs.warning = warn;
	
	if(warn != "") return true;
	return false;
}

"use strict";
// Functions when the mouse button goes up or down

/// Fires if mouse button is pressed down
function mouse_down(xx,yy,evt)                                     
{
	if(false){ prr("evt"); prr(evt);}

  let d = new Date(); 
	let over = inter.over;

	inter.mouse_time_down_prev = inter.mouse_time_down;
	inter.mouse_down_over_prev = inter.mouse_down_over; 
	
	inter.mouse_time_down = d.getTime();
	inter.mouse_down_over = copy(over);
	inter.mouse_down = {x:xx, y:yy, which:evt.which};

	let l = over.layer;
	let i = over.i;
	
	let keep_cursor_flag = false;
	
	if(l != undefined){
		let bu = inter.layer[l].but[i];
		
		switch(evt.which){	
		case 3:                                        // Right click
			{	
				switch(bu.ac){	
				case "Compartment": case "Transition":
				case "TransitionPoint":
					{
						inter.mode={ type:"Drag_Classification", p:bu.p, cl:bu.cl, l:l};
						set_arrow_icon();
					}
					break;
					
				case "ClassificationBack":
					{
						inter.mode={ type:"Drag_Selector", mx:xx, my:yy,  p:bu.p, cl:bu.cl, l:l};
						set_arrow_icon();
					}
					break;
					
				case "ClassGraphBack":
					{
						clear_comp_select();
						inter.mode={ type:"Drag_Classification", p:bu.p, cl:bu.cl, l:l};
						set_arrow_icon();
						generate_screen();
					}
					break;
				}
			}
			break;
			
		case 1:                                        // Left button
			{
				switch(bu.type){
				case "Letter":
					inter.mode={ type:"Drag_Letter", index:bu.index};
					break;
				
				case "Timebar":
					{
						let lay = get_lay("AnimControls");
						let anim = inter.graph.animation;
			
						let x = inter.mx-bu.x-lay.x;
						let xnow = bu.dy/2 + (bu.dx-bu.dy)*anim.playframe/anim.playframe_max;
						if(x >= xnow - bu.dy/2 && x <= xnow + bu.dy/2){
							inter.mode={ type:"Drag_Timebar", drag:true, bu:bu};
						}
					}
					break;
					
					case "Slider":
					{
						let info = bu.info;
						
						let lay = get_lay(info.lay);
						let anim = inter.graph.animation;
			
						let x = inter.mx-bu.x-lay.x;
						let fr = (info.value-info.min)/(info.max-info.min);
						let xnow = slider_dx/2 + (bu.dx-slider_dx)*fr;
						
						if(x >= xnow - slider_dx/2 && x <= xnow + slider_dx/2){
							inter.mode={ type:"Drag_Slider", drag:true, bu:bu};
						}
					}
					break;
					
				case "VerticalSlider":
					{
						let pos = inter.scroll_position[bu.scroll_ref];
						inter.mode={ type:"Drag_VerticalSlider", drag:true, scroll_ref:bu.scroll_ref, scale:bu.scale, pn:bu.pn, sh:pos.shift, y:inter.my};
					}
					break;
					
				case "HorizontalSlider":
					{
						let pos = inter.scroll_position[bu.scroll_ref];
						inter.mode={ type:"Drag_HorizontalSlider", drag:true, scroll_ref:bu.scroll_ref, scale:bu.scale, pn:bu.pn, sh:pos.shift, x:inter.mx};
					}
					break;
				
				case "TableSlider":
					inter.mode={ type:"Drag_TableSlider", drag:true, width_store:copy(bu.width), width:bu.width, x_st:inter.mx,  col:bu.col};
					break;
					
				case "TableSlider2":
					inter.mode={ type:"Drag_TableSlider2", drag:true, width_store:right_menu_width, x_st:inter.mx, col:bu.col};
					break;
				
				case "Compartment": case "CompLatLng": case "LabelText": 
					if(model.get_show_model() == true && inter.mode.type != "Add_Transition" &&
   					inter.mode.type != "Add_Source" && inter.mode.type != "Add_Sink"){
						
						if(bu.type == "Compartment" && evt.ctrlKey){
							inter.comp_select.list.push(bu.i);
						}
				
						let fix = false;
						if(bu.type == "Compartment" || bu.type == "CompLatLng"){
							fix = model.species[bu.p].cla[bu.cl].comp[bu.i].fixed.check;
						}
						
						if(fix == false){
							let fr = get_button_frac_over();
							inter.mode={ type:"Drag_"+bu.type, drag:true, p:bu.p, cl:bu.cl, i:bu.i, shx:(fr.x-0.5)*bu.dx, shy:(fr.y-0.5)*bu.dy}; 
						}
					}
					break;
					
				case "Selected":
					{
						// Stores positions of selected
						let cs = inter.comp_select;
						let claa = model.species[cs.p].cla[cs.cl];
						
						let fixed = false;
						
						let pos=[];
						for(let i = 0; i < cs.list.length; i++){
						
							let c = cs.list[i];
							let co = claa.comp[c];
							if(co.fixed.check == true) fixed = true;
							
							pos.push({c:c, x:co.x, y:co.y});
						}
						
						if(pos.length > 0 && fixed == false){
							inter.mode={ type:"Drag_Selected", drag:true, p:bu.p, cl:bu.cl, i:bu.i, pos:pos, mx:inter.mx, my:inter.my}; 
						}
					}
					break;
					
				case "Box":
					if(inter.mode.type != "Add_Transition" && inter.mode.type != "Add_Source" && inter.mode.type != "Add_Sink"){
						let claa = model.species[bu.p].cla[bu.cl];
						if(claa.camera.coord == "cartesian"){
							let comps = claa.annotation[bu.i].comps;
							
							let pos=[];
							for(let i = 0; i < comps.length; i++){
								let c = hash_find(claa.hash_comp,comps[i]);
								pos.push({c:c, x:claa.comp[c].x, y:claa.comp[c].y});
							}
							inter.mode={ type:"Drag_Box", p:bu.p, cl:bu.cl, i:bu.i, pos:pos, mx:inter.mx, my:inter.my}; 
						}
					}
					break;
				
					
				case "Transition":
					if(inter.mode.type != "Add_Transition" && inter.mode.type != "Add_Source" && inter.mode.type != "Add_Sink"){
						inter.mode={ type:"Split_Transition", p:bu.p, cl:bu.cl, i: bu.i, index:inter.over.index, mx_store:inter.mx, my_store:inter.my};
					}
					break;
					
				case "TransitionPoint":
					let mo = inter.mode;
					if(mo.type != "Add_Transition" && mo.type != "Add_Source" && mo.type != "Add_Sink"){
						inter.mode = { type:"Drag_TransitionPoint", p:bu.p, cl:bu.cl, i: bu.i, index: bu.index};
					}
					break;
				}

				switch(bu.ac){
				case "TimelineGrab":
					{
						inter.mode={ type:"Drag_Timeline", l:l};
						set_arrow_icon();
					}
					break;
					
				case "GraphGrab":
					{
						inter.mode={ type:"Drag_Graph", l:l};
						set_arrow_icon();
					}
					break;
				
				case "ClassificationBack": case "ClassGraphBack":
					{
						clear_comp_select();
						inter.mode={ type:"Drag_Classification", p:bu.p, cl:bu.cl, l:l};
						set_arrow_icon();
						generate_screen();
					}
					break;
				
				case "PositionCursor":
					{	
						inter.mode={ type:"Drag_Cursor"};
						let cur = inter.cursor; 	
						
						if(cur.l != l){
							if(unfocus_check() == true) generate_screen();
						}
						
						turn_off_cursor_flash();
						
						if(evt.shiftKey == true && cur.click_pos != undefined && cur.l == l){
							cur.i = bu.i;
							cur.on = 0;
						}
						else{
							let time = (new Date()).getTime();
						
							if(cur.time != undefined && time-cur.time < 400 && cur.l == l){  // Double click on text box
								let te = inter.layer[cur.l].get_text_from_source();
								
								let i_start = cur.i; while(i_start > 0 && te.substr(i_start-1,1) != " ") i_start--;
								let i_end = cur.i; while(i_end < te.length && te.substr(i_end,1) != " ") i_end++;
								
								if(cur.select_pos == i_start && cur.i == i_end){
									cur.select_pos = 0;
									cur.i = te.length;
								}
								else{
									cur.select_pos = i_start;
									cur.i = i_end;
								}	
							}
							else{
								let l_old = cur.l;
							
								inter.cursor = {l:l, i:bu.i, on:0, interval:undefined, click_pos:bu.i, time:time};	
							
								if(l_old != undefined && l_old != l){
									inter.layer[l_old].plot_buttons();
								}
							}
						}
					
						close_dropdown();
						
						inter.layer[l].plot_buttons();
						plot_screen();
						
						keep_cursor_flag = true;
					}
					break;
			
				case "ExpandIcon": keep_cursor_flag = true; break;
				
				case "Select":
					inter.mode.mx = inter.mx;
					inter.mode.my = inter.my;
					break;
					
				case "CloseBubble": case "RadioButton": case "CheckboxButton": 
				case "Dropdown": case "DropdownOption": 
				case "VerticalSlider": case "HorizontalSlider":
					keep_cursor_flag = true; 
					break;
				}
			}
			break;
		}
	}

	if(inter.equation.te != undefined) keep_cursor_flag = true;

	if(keep_cursor_flag == false){
		if(inter.cursor.i != undefined && inter.mode.type != "Drag_Cursor"){
			unfocus_check();			
			turn_off_cursor();
			
			// If button not active then replot screen
			let fl = false;
			if(l != undefined){
				let bu = inter.layer[l].but[i];
				if(bu.ac != undefined) fl = true;
			}	
			if(fl == false) generate_screen();
		}
	}
}

 
/// Fires when mouse button is released
function mouse_up(xx,yy)                                       
{
	switch(inter.mode.type){	
	case "Drag_Letter":
		mode_off();
		generate_screen();
		break;

	case "Drag_Selector":
		model.select_compartments();
		mode_off();
		generate_screen();
		break;
		
	case "Drag_Timebar": 
		mode_off();
		break;
		
	case "Drag_Slider":
		inter.graph.press_slider(inter.mode.bu,true);
		mode_off();
		generate_screen();
		break;
		
	case "Drag_VerticalSlider": case "Drag_HorizontalSlider": 
		mode_off();
		generate_screen();
		return;

	case "Drag_TableSlider": case "Drag_TableSlider2": 
		mode_off();
		return;
		
	case "Drag_Cursor":
		mode_off();
		let cur = inter.cursor; 

		if(cur.i != undefined){
			if(cur.i == cur.select_pos) cur.select_pos = undefined;
			
			start_cursor_flash();
			cur.on = 1;
			inter.layer[cur.l].plot_buttons();
			plot_screen();
		}
		return;
		
	case "Drag_Graph": 	case "Drag_Timeline": 
		mode_off();
		break;
		
	case "Drag_Classification": case "Drag_Compartment": case "Drag_CompLatLng":
	case "Drag_LabelText": case "Drag_Box": case "Drag_Selected":
	case "Drag_TransitionPoint": case "Split_Transition":
		if(inter.mode.moved == true){ mode_off(); return;}
		mode_off();
		break;
		
	case "Add_Box":
		{
			let mo = inter.mode;
			let p = mo.p, cl = mo.cl;
			add_box();
			mode_off();
			
			model.select_button_bubble("Annotation","Box",p,cl,model.species[p].cla[cl].annotation.length-1);
		}
		break;
	}
	
	let timeup = (new Date()).getTime();
	let timeclick  = timeup - inter.mouse_time_down; 
	
	if(timeclick < 350 && inter.mouse_down.which == 1){ 
		mouse_click(xx,yy);
	}
}

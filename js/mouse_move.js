"use strict";
// Functions when the mouse moves

/// Fires when the mouse moves
function mouse_move(x,y)                                   
{
	let xst = inter.mx, yst = inter.my
	inter.mx = x; inter.my = y;

	//prr("ty"); prr(inter.mode.type);
	switch(inter.mode.type){
	case "Drag_Letter":
		inter.figure[inter.mode.index].x = x-letter_size/2;
		inter.figure[inter.mode.index].y = y-letter_size/2;
		replot_layer("Figure");
		plot_screen();
		break;
			
	case "Drag_Selector":
		replot_layer("Frame");
		plot_screen();
		break;
		
	case "Drag_Timebar":
		{
			let anim = inter.graph.animation;
			let st = anim.playframe;
			inter.graph.press_timebar(inter.mode.bu);
			
			if(st != anim.playframe) inter.graph.replot_anim();
		}
		break;
		
	case "Drag_Slider":
		{
			inter.graph.press_slider(inter.mode.bu);	
			replot_layer(inter.mode.bu.info.lay);
			plot_screen();
		}
		break;
		
	case "Drag_VerticalSlider": case "Drag_HorizontalSlider":
		{	
			let pos = inter.scroll_position[inter.mode.scroll_ref];
		
			let sh;
			if(inter.mode.type == "Drag_VerticalSlider") sh	= inter.mode.sh + (y-inter.mode.y)*pos.scale;
			else sh	= inter.mode.sh + (x-inter.mode.x)*pos.scale;
			
			/*
			let sh = pos.shift;
			if(inter.mode.type == "Drag_VerticalSlider") sh += (y-yst)*pos.scale;
			else sh += (x-xst)*pos.scale;
			*/
			
			if(sh < 0) sh = 0;
			if(sh > pos.max*(1-pos.frac)) sh = pos.max*(1-pos.frac);
			pos.shift = sh;

			let lay = inter.layer[pos.scroll_layer];
		
			if(bubble_on() == true) generate_screen();
			else{
				lay.initialise_buttons();
				lay.plot_buttons();
				lay.op.lay.plot_buttons();
					
				plot_screen();
			}
		}
		break;

	case "Drag_Cursor":
		drag_cursor(x,y);
		break;
	
	case "Drag_TableSlider":
		drag_table_slider(x);
		break;
		
	case "Drag_TableSlider2":
		drag_rightmenu_slider(x);
		break;
		
	case "Drag_Classification":
		{
			let show = model.get_show_model();
			
			let cam;
			if(show == true) cam = model.species[inter.mode.p].cla[inter.mode.cl].camera;
			else cam = inter.graph.get_cla(inter.mode.p,inter.mode.cl).camera;
		
			let lay = {dx:0, dy:0};
			let pst = trans_point_rev(xst,yst,cam,lay);
			let p = trans_point_rev(x,y,cam,lay);

			let sh = {x:x-xst, y:y-yst};
			cam.x -= p.x-pst.x;
			cam.y -= p.y-pst.y;
		
			if(inter.bubble.lay_name != undefined) generate_screen();
			else model.replot(sh);
			
			model.ensure_lng(cam);
			
			if(show == true) clone_camera(inter.mode.p,inter.mode.cl);
		}
		break;
		
	case "Drag_Graph":
		{
			let mo = inter.mode;
			let shx = inter.mx-xst 
			let shy = inter.my-yst;
			let la = inter.layer[inter.mode.l];
			let w = la.dx, h = la.dy;
			let ra = inter.graph.range;
			let dx = -(ra.xmax-ra.xmin)*shx/w;
			let dy = (ra.ymax-ra.ymin)*shy/h;
			ra.xmin += dx; ra.xmax += dx;
			if(!(inter.graph.op.yaxis == false)){
				ra.ymin += dy; ra.ymax += dy;
			}
			inter.graph.replot();
		}
		break;
		
	case "Drag_Timeline":
		{
			let mo = inter.mode;
			let shx = inter.mx-xst;
			let la = inter.layer[inter.mode.l];
			let w = la.dx;
			let ra = inter.graph.range;
			let dx = -(ra.xmax-ra.xmin)*shx/w;
			ra.xmin += dx; ra.xmax += dx;
			
			let shy = inter.my-yst;
			
			if(shy != 0) scroll_layer("GraphContent",-shy);
			inter.graph.replot();
		}
		break;
		
	case "Drag_Box": case "Drag_Selected":
		{
			let mo = inter.mode;
			let p = mo.p;
			let cl = mo.cl;
			let i = mo.i;
			let pos = mo.pos; 
			
			let shx = inter.mx-mo.mx;
			let shy = -(inter.my-mo.my);
		
			let claa = model.species[p].cla[cl];
			let cam = claa.camera;
			if(cam.coord == "latlng"){
				let sc = Number(cam.scale);
				shx /= sc; shy /= sc;
				
				// Ensures that does not go out of range
				let max = -LARGE, min = LARGE;
				for(let i = 0; i < pos.length; i++){
					let xx = pos[i].x + shx;
					if(xx > max) max = xx;
					if(xx < min) min = xx;
				}
				
				if(max > Math.PI) shx -= max-Math.PI + TINY;
				if(min < -Math.PI) shx += -Math.PI-min + TINY;
			}
		
			let mgroup=[];	
			for(let i = 0; i < pos.length; i++){
				mgroup.push({c:pos[i].c, x:pos[i].x + shx, y:pos[i].y + shy});
			}				
			
			model.move_group(mgroup,cl,p);
		}
		break;
	
	case "Drag_LabelText":
		{
			let p = inter.mode.p;
			let cl = inter.mode.cl;
			let i = inter.mode.i;
			
			let claa = model.species[p].cla[cl];
			let cam = claa.camera;
			
			let mo = inter.mode;
	
			let lay = get_lay("Compartment");
			let po = trans_point_rev(x-lay.x-mo.shx,y-lay.y-mo.shy,cam,lay);
			let ob = claa.annotation[i];
			let cx = ob.x, cy = ob.y;
			ob.x = po.x;
			ob.y = po.y;
			if(cam.coord == "cartesian") model.snap_comp_to_grid(ob);
			
			if(cx != ob.x || cy != ob.y){
				model.replot();
			}
		}
		break;
		
	case "Drag_Compartment": case "Drag_CompLatLng": 
		{
			let p = inter.mode.p;
			let cl = inter.mode.cl;
			let i = inter.mode.i;
			
			let claa = model.species[p].cla[cl];
			let cam = claa.camera;
			
			let mo = inter.mode;
	
			let lay = get_lay("Compartment");
			let po = trans_point_rev(x-lay.x-mo.shx,y-lay.y-mo.shy,cam,lay);
			
			let mgroup=[];	
			mgroup.push({c:i, x:po.x, y:po.y});
			model.move_group(mgroup,cl,p);
		}
		break;
	
	case "Drag_TransitionPoint":
		{
			let mp = model.get_mouse_to_desktop();
			let mo = inter.mode;
			
			let pos = model.species[mo.p].cla[mo.cl].tra[mo.i].midp[mo.index];
			
			if(pos.x != mp.x || pos.y != mp.y){
				let claa = model.species[mo.p].cla[mo.cl];
				claa.tra[mo.i].midp[mo.index] = mp;
				model.find_trans_pline(claa,mo.i);   
			
				if(inter.bubble.lay_name != undefined){
					generate_screen();
				}
				else{
					replot_layer("Transition");
					plot_screen();
				}
			}
		}
		break;
		
	case "Split_Transition":
		{
			let dx = inter.mx - inter.mode.mx_store, dy = inter.my - inter.mode.my_store;
			if(dx*dx + dy*dy > SPLIT_SIZE*SPLIT_SIZE){
				model.split_transitionpoint();
				replot_layer("Transition");
				plot_screen();
			}
		}
		break;
	
	case "Add_Box":
		replot_layer("Select");
		plot_screen();
		break;
		
	default:
		{	
			switch(inter.mode.type){
			case "Add_Compartment": model.place_add_compartment(inter.mode.type); break;
			case "Add_Label": model.place_add_compartment(inter.mode.type); break;
			}
			
			let over_new = get_button_over(x,y);
			
			let replot_flag = false; 
		
			if(inter.mode.type == "Add_Transition" || inter.mode.type == "Add_Source" || inter.mode.type == "Add_Sink") replot_flag = true;
				
			if(over_new.layer != inter.over.layer || over_new.i != inter.over.i){
				let l = inter.over.layer;
				let i = inter.over.i;
				if(l){
					let lay = inter.layer[l];
			
					if(lay.name == "Transition") replot_flag = true;
					else lay.plot_button(lay.but[i],false);
				}
	
				l = over_new.layer;
				i = over_new.i;
				if(l){
					let lay = inter.layer[l];
					if(lay.name == "Transition") replot_flag = true;
					else lay.plot_button(lay.but[i],true);
				}
			
				inter.over = over_new;
				
				if(replot_flag == true) replot_layer("Transition");
				
				plot_screen();
			}
			else{
				if(replot_flag == true){
					replot_layer("Transition");	 
					plot_screen();
				}
			}
		}
		break;
	}

	set_arrow_icon();
}


/// Sets the arrow icon depending on what the arrow is over
function set_arrow_icon()
{
	let arrow_icon = "";
	let l = inter.over.layer;

	if(l != undefined && inter.mode.type != "Add_Compartment"){
		switch(inter.layer[l].but[inter.over.i].ac){
			case "TableSlider": case "TableSlider2": arrow_icon = "ColShift"; break;
		}
		
		switch(inter.mode.type){
		case "Drag_Classification": case "Drag_Timeline": case "Drag_Graph":
			arrow_icon = "Grab";
			break;
		}
	}
	
	// Sets the cursor (i.e. a hand is used for grabbing)
	if(arrow_icon != inter.arrow_icon){        
		inter.arrow_icon = arrow_icon;
		let val = "";
		switch(arrow_icon){
		case "Grab": val = "-webkit-grab"; break;
		case "ColShift": val = "ew-resize"; break;
		}
		
		document.getElementById("bod").style.cursor = val;
	}
}

	
/// Retuns which button the mouse is over
function get_button_over(x,y)
{		
	for(let l = inter.layer.length-1; l >= 0; l--){
		let lay = inter.layer[l];
		
		if(x >= lay.x && x <= lay.x+lay.dx && y >= lay.y && y <= lay.y+lay.dy){
			let sh_x = lay.x-lay.x_shift;
			let sh_y = lay.y-lay.y_shift;
		
			for(let i = lay.but.length-1; i >= 0; i--){
				let bu = lay.but[i];
				if(bu.ac && x >= sh_x+bu.x && x <= sh_x+bu.x+bu.dx 
				         && y >= sh_y+bu.y && y <= sh_y+bu.y+bu.dy){
					
					switch(bu.type){
					case "Transition":
						{
							let index = model.mouse_over_transition(bu.points,lay);
							if(index != undefined) return { layer:l, i:i, index:index};
						}
						break;
						
					case "CompMap": case "CompMapGraph":
						{
							let frx = (x-(sh_x+bu.x))/bu.dx;
							let fry = 1-(y-(sh_y+bu.y))/bu.dy;
							if(mouse_over_comp_map(frx,fry,bu.mask) == true){
								return { layer:l, i:i};
							}
						}
						break;
						
					default:
						{
							let bub = inter.bubble;
							
							if((bub.op && bub.op.type == "Popup") || 
								!(bub.lay_name == lay.name && bub.i == i)){
								return { layer:l, i:i};
							}
						}
						break;
					}
				}
			}
		}
	}		
	
	return {};
}


/// Gets the fraction the mouse is across the button
function get_button_frac_over()
{
	let l = inter.over.layer;
	let i = inter.over.i;
	
	let lay = inter.layer[l];
	let bu = lay.but[i];
	
	let sh_x = lay.x-lay.x_shift;
	let sh_y = lay.y-lay.y_shift;
	
	return { x:(inter.mx-(sh_x+bu.x))/bu.dx, y:(inter.my-(sh_y+bu.y))/bu.dy};
}

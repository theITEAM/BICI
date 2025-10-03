"use strict";
// Functions for individual effects

/// Display an individual effect group
function display_ind_eff_group(p,i,x,y,lay,wmax)
{		
	let ieg = model.species[p].ind_eff_group[i];
	let ie_list = ieg.ie_list;
	
	//ieg.A_matrix.name = get_A_matrix_name(ie_list);
	
	lay.add_checkbox(wmax-4,y+0.3,"Ind. Cor.","Ind. Cor.",ieg.A_matrix,WHITE,{title:"Individual correlation", te:corrlation_text});
	wmax -= 4;
	
	if(ie_list.length == 1){	
		if(combine_possible(i,p) == true){
			lay.add_button({te:"Combine", x:wmax-4.8, y:y+0.2, dx:4.4, dy:1.2, ac:"CombineIEBubble", type:"GreyView", p:p, i:i});
			wmax -= 5.3;
		}
	}
	else{
		lay.add_button({te:"Split", x:wmax-4, y:y+0.2, dx:3.5, dy:1.2, ac:"SplitIE", type:"GreyView", p:p, i:i});
		wmax -= 4.5;
	}
	
	let fo = get_font(1.2,"","times");
	
	var si = 1.5;
	
	if(ie_list.length == 0) error("Should contain IE");
	
	let te_si = 0.8, te_fo = get_font(te_si,"","times");
	
	lay.add_button({te:"NAME:", x:x, y:y+0.45, dy:te_si, si:te_si, font:te_fo, type:"Text", col:BLACK});

	let fo_name = get_font(1.1,"","times");
	lay.add_button({te:ieg.name, x:x+2.7, y:y+0., dx:4, dy:1.6, type:"IEGroupName", font:fo_name, ac:"IEGroupName", source:ieg });
	
	let xx = x + 8;
	
	if(ie_list.length == 1){
		let st = ie_list[0].name;
		let st2 = "["+st+"]";
		let w = text_width(st2,fo);
	
		lay.add_button({te:st2, x:xx, y:y, dy:si_big, type:"Text", font:fo, si:si_big, col:BLACK});
	
		lay.add_button({te:"~", x:xx+w+0.6, y:y, dy:si, type:"Text", font:get_font(si), si:si, col:BLACK});
		
		let eq = "Ω^"+ieg.name; 
		let te = "<e>Shifted LN("+eq;
		
		if(ieg.A_matrix.check == true){
			te += "</e><e>⊗</e><e>A^"+ieg.name+"</e>";
		}
		te += ")</e>";
		
		lay.add_paragraph(te,wmax-x-w-0.6-1.6,xx+w+0.6+1.6,y+0.2,BLACK,para_eq_si,para_lh);
	}
	else{
		let st = "";
		for(let j = 0; j < ie_list.length; j++){
			if(st != "") st += ", ";
			st += "["+ie_list[j].name+"]";
		}
		st = "("+st+")";
		
		let w = text_width(st,fo);
		
		//if(x-w-0.6 < 2) x = 2+w+0.6;
		lay.add_button({te:st, x:xx, y:y, dy:si_big, type:"Text", font:fo, si:si_big, col:BLACK});
		
		lay.add_button({te:"~", x:xx+w+0.6, y:y, dy:si, type:"Text", font:get_font(si), si:si, col:BLACK});
		
		let eq = "Ω^"+ieg.name; 
		let te = "<e>Shifted MVLN("+eq;
		
		if(ieg.A_matrix.check == true) te += "</e><e>⊗</e><e>A^"+ieg.name;
		te += ")</e>";
		
		lay.add_paragraph(te,wmax-x-w-0.6-1.6,xx+w+0.6+1.6,y+0.2,BLACK,para_eq_si,para_lh);
	}
}


/// Returns the name of the A matrix
function get_A_matrix_name(ie_list)
{
	
	let st = "";
	for(let j = 0; j < ie_list.length; j++){
		if(st != "") st += ",";
		st += ie_list[j].name;
	}
	return "A^"+st;
}


/// Bubble for selecting which IE to combine with
function combine_ie_bubble(cont)
{
	cont.dx = 15;
	bubble_addtitle(cont,"Combine individual effects",{});
	
	let p = inter.bubble.ie_p;
	let ie = inter.bubble.ie_combine;
	let name = model.species[p].ind_eff_group[ie].ie_list[0].name;

	bubble_addparagraph(cont,"Select individual effect to combine <e>"+name+"</e> with:",0,cont.dx);
	
	cont.y += 0.2;

	bubble_addscrollable(cont,{type:"combine IE", p:p, ymax:10});
	add_bubble_end(cont);
}
	
	
/// Allows individual effects to be combined
function combineIE_scrollable(lay)
{
	let cx = 0, cy = 0;
	let dy = 2.4;
	let gap = 0.5;
	
	let si = si_big;
	let fo = get_font(si,"","times");
	
	let sp = model.species[lay.op.p];
	for(let i = 0; i < sp.ind_eff_group.length; i++){
		let ieg = sp.ind_eff_group[i];
		if(i != inter.bubble.ie_combine){
			let te = "";
	
			for(let j = 0; j < ieg.ie_list.length; j++){
				if(te != "") te += ",";
				te += ieg.ie_list[j].name;
			}
		
			let dx = 1+text_width(te,fo);
			if(cx+dx+0.3 >= lay.inner_dx){ cx = 0; cy += dy;}
			
			lay.add_button({te:te, x:cx, y:cy, dx:dx, dy:1.4*si, type:"IEcombine", back_col:BUBBLE_COL, font:fo, i:i, ac:"CombineIE"});
			
			cx += dx+0.3;
		}
	}
	
	cy += dy;
	
	return cy;
}


/// Determines if it is possible to combine indiv
function combine_possible(isel,psel)
{
	let sp = model.species[psel];
	for(let i = 0; i < sp.ind_eff_group.length; i++){
		let ieg = sp.ind_eff_group[i];
		if(i != isel){
			if(ieg.ie_list.length > 0) return true;
		}
	}
	
	return false;
}


/// Adds one individual effect group to another
function add_ie_to_group(p,i1,i2)
{
	let sp = model.species[p];
	let ieg1 = sp.ind_eff_group[i1];
	let ieg2 = sp.ind_eff_group[i2];
	for(let i = 0; i  < ieg1.ie_list.length; i++){
		ieg2.ie_list.push({name:ieg1.ie_list[i].name});
	}
	
	sp.ind_eff_group.splice(i1,1);
}


/// Splits a group of individual effects into individual elements
function split_ie(p,i)
{
	let sp = model.species[p];
	let ieg = sp.ind_eff_group[i];
	for(let j = 0; j < ieg.ie_list.length; j++){
		sp.ind_eff_group.push(copy(ieg));
		let ieg2 = sp.ind_eff_group[sp.ind_eff_group.length-1];
		ieg2.ie_list.length = 0;
		ieg2.ie_list.push(copy(ieg.ie_list[j]));
		if(j > 0) ieg2.name = generate_iegroup_name();
	}
	
	sp.ind_eff_group.splice(i,1);
	
	update_param();
}


/// Allows the user to edit Amatrix values
function create_edit_Amatrix(lay)
{
	lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, type:"Rect", val:WHITE});
	
	let cx = corner.x;
	let cy = corner.y;
	
	let Amat = inter.edit_Amatrix;
	
	let too_big = false; if(Amat.too_big == true) too_big = true;
	let title = "Edit A matrix";
	if(Amat.type == "Ainv") title = "Edit inverse A matrix";
	
	if(too_big) title += " (too large to edit)";
	
	cy = lay.add_title(title,cx,cy,{te:Amatrix_text});
	
	add_layer("CreateEditAmatrixContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-3.5,{type:lay.op.type, too_big:too_big});
	
	if(too_big == true){
		lay.add_corner_button([["Back","Grey","CancelAmatrix"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	}
	else{
		lay.add_corner_button([["Cancel","Grey","CancelAmatrix"],["Done","Grey","DoneAmatrix"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	}
}


/// Closes the Amatrix being editted
function close_Amatrix_source()
{
	inter.edit_Amatrix = {};
}


/// Plots A matrix to allow for values to be editted			
function add_Amatrix_buts(lay)
{
	let edit_A = inter.edit_Amatrix;
	let ind_list = edit_A.ind_list;
	
	let sp = model.species[edit_A.p];
	let Amat = sp.ind_eff_group[edit_A.i].A_matrix
	
	let si_mar = 1;
	let fo_mar = get_font(si_mar);
	let fo_table = get_font(si_table);

	let w_max = 0;
	for(let i = 0; i < ind_list.length; i++){
		let w = text_width(ind_list[i],fo_mar);
		if(w > w_max) w_max = w;		
	}
	w_max += 1;
	
	let mar = 0.4;
	let gap = 0.7;
	let cy = 0;
	let dx = 3;
	
	// Finds the longest element
	let str = "";
	for(let j = 0; j < ind_list.length; j++){
		for(let i = 0; i < ind_list.length; i++){
			let st = String(inter.edit_Amatrix.A_value[j][i]);
			if(st.length > str.length) str = st;
		}
	}
	let w_ele = 1.2*text_width(str,fo_table);
	if(w_ele > dx) dx = w_ele;
	if(dx > 8) dx = 8;
	
	let mar_col = DRED;
	
	let cx = 2+w_max+gap;
	for(let i = 0; i < ind_list.length; i++){
		lay.add_button({te: ind_list[i], x:cx+dx/2-0.6, y:cy, dx:1.2, dy:w_max, type:"VertText", font:fo_mar, col:mar_col});
		cx += dx;
	}
	cy += w_max+gap;
	
	let too_big = edit_A.too_big;
	
	cx = 2;
	let out_dx = dx*ind_list.length+2*mar;
	let out_dy = dy_table_param*ind_list.length+2*mar;
		
	if(too_big){ out_dx += 2; out_dy += dy_table_param;}
	
	lay.add_button({x:cx+w_max+gap-mar, y:cy-mar, dx:out_dx, dy:out_dy, type:"Outline", col:BLACK});
	
	let ele_type = "AmatrixElement";
	let action = "EditAmatrixElement";

	if(too_big){ ele_type = "TooBigElement"; action = undefined;}
	
	for(let j = 0; j < ind_list.length; j++){
		let cx = 2;
		lay.add_button({te: ind_list[j], x:cx, y:cy, dx:w_max, dy:dy_table_param, type:"RightText", si:si_mar, font:fo_mar, col:mar_col});
	
		cx += w_max+gap;
	
		for(let i = 0; i < ind_list.length; i++){
			let val = inter.edit_Amatrix.A_value[j][i];
			lay.add_button({te:val, x:cx, y:cy, dx:dx, dy:dy_table_param, type:ele_type, font:fo_table, i:i, j:j, ac:action});
			cx += dx;
		}
		if(too_big){
			lay.add_button({te:"...", x:cx, y:cy, dx:2, dy:dy_table_param, type:"Text", si:si_mar, font:fo_mar, col:mar_col});
		}
		cy += dy_table_param;
	}
	
	if(too_big == true){
		let cx = 2+w_max+gap;
		for(let i = 0; i < ind_list.length; i++){
			lay.add_button({te:"⋮", x:cx, y:cy, dx:2, dy:dy_table_param, type:"RightText", si:si_mar, font:fo_mar, col:mar_col});
			cx += dx;
		}
	}
	
	lay.add_button({x:0, y:cy, dx:0, dy:0.5, type:"Nothing"});	
}

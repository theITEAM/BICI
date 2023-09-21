"use strict";
/// Functions for individual effects

/// Display an individual effect group
function display_ind_eff_group(p,i,x,y,lay,wmax)
{		
	let ieg = model.species[p].ind_eff_group[i];
	let list = ieg.list;
	
	lay.add_checkbox(wmax-4,y+0.3,"Ind. Cor.","Ind. Cor.",ieg.A_matrix,WHITE,{title:"Individual correlation", te:corrlation_text});
	wmax -= 4;
	
	if(list.length == 1){	
		if(combine_possible(i,p) == true){
			lay.add_button({te:"Combine", x:wmax-4, y:y+0.3, dx:3.7, dy:1, ac:"CombineIEBubble", type:"CombineIE", p:p, i:i});
			wmax -= 4.5;
		}
	}
	else{
		lay.add_button({te:"Split", x:wmax-4, y:y+0.3, dx:3.7, dy:1, ac:"SplitIE", type:"CombineIE", p:p, i:i});
		wmax -= 4.5;
	}
	
	let fo = get_font(1.2,"","times");
	
	var si = 1.5;
	
	if(list.length == 0) error("Should contain IE");
	
	if(list.length == 1){
		let st = list[0].name;
		let st2 = "["+st+"]";
		let w = text_width(st2,fo);

		lay.add_button({te:st2, x:x-w-0.6, y:y, dy:si_big, type:"Text", font:fo, si:si_big, col:BLACK});
	
		lay.add_button({te:"~", x:x, y:y, dy:si, type:"Text", font:get_font(si), si:si, col:BLACK});
		
		let eq = "Ω^"+st+st; 
		let te = "<e>Shifted LN("+eq;
		
		if(ieg.A_matrix.check == true) te += "</e><e>×</e><e>"+ieg.A_matrix.name+"</e>";
		te += ")</e>";
		
		lay.add_paragraph(te,wmax-x-1.6,x+1.6,y+0.2,BLACK,para_eq_si,para_lh);
	}
	else{
		let st = "", st2 = "";
		for(let j = 0; j < list.length; j++){
			if(st += "") st += ", ";
			st += "["+list[j].name+"]";	st2 += list[j].name;
		}
		st = "("+st+")";
		
		let w = text_width(st,fo);
		
		if(x-w-0.6 < 2) x = 2+w+0.6;
		lay.add_button({te:st, x:x-w-0.6, y:y, dy:si_big, type:"Text", font:fo, si:si_big, col:BLACK});
		
		lay.add_button({te:"~", x:x, y:y, dy:si, type:"Text", font:get_font(si), si:si, col:BLACK});
		
		let eq = "Ω^"+st2; 
		let te = "<e>Shifted MVLN("+eq;
		
		if(ieg.A_matrix.check == true) te += "</e><e>⊗</e><e>"+ieg.A_matrix.name;
		te += ")</e>";
		
		lay.add_paragraph(te,wmax-x-1.6,x+1.6,y+0.2,BLACK,para_eq_si,para_lh);
	}
	
	ieg.A_matrix.name = get_A_matrix_name(list);
}


/// Returns the name of the A matrix
function get_A_matrix_name(list)
{
	if(list.length == 1) return "A^"+list[0].name;
	else{
		let st = "";
		for(let j = 0; j < list.length; j++) st += list[j].name;
		return "A^"+st;
	}
}


/// Bubble for selecting which IE to combine with
function combine_ie_bubble(cont)
{
	cont.dx = 15;
	bubble_addtitle(cont,"Combine individual effects",{});
	
	let p = inter.bubble.ie_p;
	let ie = inter.bubble.ie_combine;
	let name = model.species[p].ind_eff_group[ie].list[0].name;

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
	
			for(let j = 0; j < ieg.list.length; j++){
				if(te != "") te += ",";
				te += ieg.list[j].name;
			}
		
			let dx = 1+text_width(te,fo);
			if(cx+dx+1 >= lay.inner_dx){ cx = 0; cy += dy;}
			
			lay.add_button({te:te, x:cx, y:cy, dx:dx, dy:1.4*si, type:"IEcombine", back_col:BUBBLE_COL, font:fo, i:i, ac:"CombineIE"});
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
			if(ieg.list.length > 0) return true;
		}
	}
	
	return false;
}


// Adds one individual effect group to another
function add_ie_to_group(p,i1,i2)
{
	let sp = model.species[p];
	let ieg1 = sp.ind_eff_group[i1];
	let ieg2 = sp.ind_eff_group[i2];
	for(let i = 0; i  < ieg1.list.length; i++){
		ieg2.list.push({name:ieg1.list[i].name});
	}
	
	sp.ind_eff_group.splice(i1,1);
	
	init_param();
}


/// Splits a group of individual effects into individual elements
function split_ie(p,i)
{
	let sp = model.species[p];
	let ieg = sp.ind_eff_group[i];
	for(let j = 0; j < ieg.list.length; j++){
		sp.ind_eff_group.push(copy(ieg));
		let ieg2 = sp.ind_eff_group[sp.ind_eff_group.length-1];
		ieg2.list.length = 0;
		ieg2.list.push(copy(ieg.list[j]));
	}
	
	sp.ind_eff_group.splice(i,1);
	
	init_param();
}


/// The function is called when the A matrix is loaded
function A_matrix_loaded(tab)
{
	if(tab.ncol == 0){ alertp("The data table does not have any columns"); return;}
	if(tab.nrow == 0){ alertp("The data table does not have any rows"); return;}
	if(tab.nrow != tab.ncol){ alertp("The data table must be square and have headings which give the IDs of the individuals"); return;}
	
	for(let r = 0; r < tab.nrow; r++){
		for(let c = 0; c < tab.ncol; c++){
			if(isNaN(tab.ele[r][c])){
				alertp("The element '"+tab.ele[r][c]+"' is not a number (col "+(c+1)+", row "+(r+2)+")"); 
				return;
			}
		}
	}
	
	let sp = model.species[inter.edit_source.info.p];
	
	let Amat = sp.ind_eff_group[inter.edit_source.info.i].A_matrix;
	Amat.ind_list = tab.heading;
	Amat.loaded = true;
	
	Amat.value = [];
	for(let r = 0; r < tab.nrow; r++){
		Amat.value[r] = [];
		for(let c = 0; c < tab.ncol; c++){
			Amat.value[r][c] = Number(tab.ele[r][c]);
		}
	}
	
	close_help();
}


/// Allows the user to edit Amatrix values
function create_edit_Amatrix(lay)
{
	lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, type:"Rect", val:WHITE});
	
	let cx = corner.x;
	let cy = corner.y;
	
	cy = lay.add_title("Edit A matrix",cx,cy,{te:Amatrix_text});
	
	add_layer("CreateEditAmatrixContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-3.5,{type:lay.op.type});
	
	lay.add_corner_button([["Cancel","Grey","CancelAmatrix"],["Done","Grey","DoneAmatrix"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	 
	//let x = 1.2, y = lay.dy-1.6;
	//let gap = 3.5;
	
//	let w = model.add_object_button(lay,"Load",x,y,load_ac,{ back:WHITE, active:true, info:{}, title:load_title, te:load_te}); 
	//x += w+gap;
}


/// Closes the Amatrix being editted
function close_Amatrix_source()
{
	inter.edit_Amatrix = {};
}


/// Plots A matrix to allow for values to be editted			
function add_Amatrix_buts(lay)
{
	let sp = model.species[inter.edit_Amatrix.p];
	let Amat = sp.ind_eff_group[inter.edit_Amatrix.i].A_matrix
	
	let si_mar = 1;
	let fo_mar = get_font(si_mar);
	let fo_table = get_font(si_table);

	let w_max = 0;
	for(let i = 0; i < Amat.ind_list.length; i++){
		let w = text_width(Amat.ind_list[i],fo_mar);
		if(w > w_max) w_max = w;		
	}
	w_max += 1;
	
	let mar = 0.4;
	let gap = 0.7;
	let cy = 0;
	let dx = 4;
	let mar_col = DRED;
	
	let cx = 2+w_max+gap;
	for(let i = 0; i < Amat.ind_list.length; i++){
		lay.add_button({te: Amat.ind_list[i], x:cx+dx/2-0.6, y:cy, dx:1.2, dy:w_max, type:"VertText", font:fo_mar, col:mar_col});
		cx += dx;
	}
	cy += w_max+gap;
	
	cx = 2;
	lay.add_button({x:cx+w_max+gap-mar, y:cy-mar, dx:dx*Amat.ind_list.length+2*mar, dy:dy_table_param*Amat.ind_list.length+2*mar, type:"Outline", col:BLACK});
		
	for(let j = 0; j < Amat.ind_list.length; j++){
		let cx = 2;
		lay.add_button({te: Amat.ind_list[j], x:cx, y:cy, dx:w_max, dy:dy_table_param, type:"RightText", si:si_mar, font:fo_mar, col:mar_col});
	
		cx += w_max+gap;
	
		for(let i = 0; i < Amat.ind_list.length; i++){
			let val = inter.edit_Amatrix.value[j][i];
			lay.add_button({te:val, x:cx, y:cy, dx:dx, dy:dy_table_param, type:"AmatrixElement", font:fo_table, i:i, j:j, ac:"EditAmatrixElement"});
			cx += dx;
		}
		
		cy += dy_table_param;
	}
	
	lay.add_button({x:0, y:cy, dx:0, dy:0.5, type:"Nothing"});	
}


/// Runs when editting of the A matrix is complete
function done_Amatrix()
{
	let ea = inter.edit_Amatrix;
	
	model.species[ea.p].ind_eff_group[ea.i].A_matrix.value = ea.value;
	close_bubble();
	close_Amatrix_source();
}


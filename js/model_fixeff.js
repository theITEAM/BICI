"use strict";
/// Functions for fixed effects

/// Display an individual effect 
function display_fix_eff(p,i,x,y,lay,wmax)
{		
	let feg = model.species[p].fix_eff[i];

	let fo = get_font(si_big,"","times");
	
	var si = 1.2;
	
	let st = "〈"+feg.name+"〉";
	
	let w = text_width(st,fo);

	if(x-w-0.6 < 2) x = 2+w+0.6;

	lay.add_button({te:st, x:x-w-0.6, y:y, dy:si_big, type:"Text", font:fo, si:si_big, col:BLACK});
	
	lay.add_button({te:"∝", x:x, y:y+0.2, dy:si, type:"Text", font:get_font(si), si:si, col:BLACK});

	let te = "<e>exp(X^r_i</e><e>×</e><e>"+feg.name+")</e>";
	lay.add_paragraph(te,wmax-x,x+1.6,y+0.2,BLACK,para_eq_si,para_lh);// zzz
}


/// Loads vector X from table
function inf_eff_load_X()
{
	let so = inter.edit_source;

	let i = so.info.i;
	let X = model.species[so.info.p].fix_eff[i].X_vector;
	
	X.loaded = true;
	X.value = [];
	X.ind_list = [];
	let tab = so.table;
	for(let j = 0; j < tab.nrow; j++){
		X.ind_list.push(tab.ele[j][0]);
		X.value.push(Number(tab.ele[j][1]));
	}
	
	close_data_source();
}


/// Allows the user to edit X vector values
function create_edit_Xvector(lay)
{
	lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, type:"Rect", val:WHITE});
	
	let cx = corner.x;
	let cy = corner.y;
	
	cy = lay.add_title("Edit X vector",cx,cy,{te:editxvector_text});
	
	add_layer("CreateEditXvectorContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-3.5,{type:lay.op.type});
	
	lay.add_corner_button([["Cancel","Grey","CancelXvector"],["Done","Grey","DoneXvector"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
}


/// Closes the X vector being editted
function close_Xvector_source()
{
	inter.edit_Xvector = {};
}


/// Plots X vector to allow for values to be editted			
function add_Xvector_buts(lay)
{
	let Xvec = model.species[inter.edit_Xvector.p].fix_eff[inter.edit_Xvector.i].X_vector;
	
	let si_mar = 1;
	let fo_mar = get_font(si_mar);
	let fo_table = get_font(si_table);

	let w_max = 0;
	for(let i = 0; i < Xvec.ind_list.length; i++){
		let w = text_width(Xvec.ind_list[i],fo_mar);
		if(w > w_max) w_max = w;		
	}
	w_max += 1;
	
	let mar = 0.4;
	let gap = 0.7;
	let cy = 1;
	let dx = 10;
	let mar_col = DRED;
	
	let cx = 2+gap;
	
	cx = 2;
	lay.add_button({x:cx+w_max+gap-mar, y:cy-mar, dx:dx+2*mar, dy:dy_table_param*Xvec.ind_list.length+2*mar, type:"Outline", col:BLACK});
		
	for(let j = 0; j < Xvec.ind_list.length; j++){
		let cx = 2;
		lay.add_button({te: Xvec.ind_list[j], x:cx, y:cy, dx:w_max, dy:dy_table_param, type:"RightText", si:si_mar, font:fo_mar, col:mar_col});
	
		cx += w_max+gap;
	
		let val = inter.edit_Xvector.value[j];
		lay.add_button({te:val, x:cx, y:cy, dx:dx, dy:dy_table_param, type:"XvectorElement", font:fo_table, i:j, ac:"EditXvectorElement"});
	
	/*
		for(let i = 0; i < Amat.ind_list.length; i++){
			
			lay.add_button({te:val, x:cx, y:cy, dx:dx, dy:dy_table_param, type:"AmatrixElement", font:fo_table, i:i, j:j, ac:"EditAmatrixElement"});
			cx += dx;
		}
		*/
		cy += dy_table_param;
	}
	
	lay.add_button({x:0, y:cy, dx:0, dy:0.5, type:"Nothing"});	
}


/// Runs when editing of the X vector is complete
function done_Xvector()
{
	let ea = inter.edit_Xvector;
	model.species[ea.p].fix_eff[ea.i].X_vector.value = ea.value;
	close_bubble();
	close_Xvector_source();
}


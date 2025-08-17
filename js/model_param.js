"use strict";
// Functions for model parameters

/// Adds the button for the page in which simulation parameters are input
function add_model_param_buts(lay)
{
	let cx = corner.x;
	let cy = corner.y;
	
	
	cy = lay.add_title("Model parameters",cx,cy,{te:modelparam_text});

	let active = true;
	if(model.param.length == 0){
		active = false;
	}
	
	cy += 0.5;

	add_layer("ModelParamContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-3.5,{});	
	
	let x = 1.2, y = lay.dy-1.6;
	let gap = 3.5;
	
	let w = model.add_object_button(lay,"Constant",x,y,"SetConstant",{ back:WHITE, active:active, info:{}, title:"Constant", te:constant_text}); 
	x += w+gap;
	
	w = model.add_object_button(lay,"Reparam.",x,y,"SetReparam",{ back:WHITE, active:active, info:{}, title:"Reparameterisation", te:reparam_text}); 
	x += w+gap;
	
	w = model.add_object_button(lay,"Distribution",x,y,"SetDistribution",{ back:WHITE, active:active, info:{}, title:"Distribution", te:distribution_text}); 
	x += w+gap;
	
	w = model.add_object_button(lay,"Derived",x,y,"SetDerived",{ back:WHITE, active:active, info:{}, title:"Derived", te:derived_text}); 
	x += w+gap;
	
	w = model.add_object_button(lay,"Factor",x,y,"SetFactor",{ back:WHITE, active:active, info:{}, title:"Facto", te:factor_text}); 
	x += w+gap;
}


/// Adds the button for the parameter mult page
function add_param_mult_buts(lay)
{
	let cx = corner.x;
	let cy = corner.y;
	
	cy = lay.add_title("Model parameter multipliers",cx,cy,{te:modelparam_text});

	let active = true;
	if(model.param.length == 0){
		active = false;
	}
	
	//cy += 0.5;

	add_layer("ParamMultContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-3.5,{});	
	
	let x = 1.2, y = lay.dy-1.6;
	let gap = 3.5;
	
	let w = model.add_object_button(lay,"Multiplier",x,y,"AddParamMult",{ back:WHITE, active:active, info:{}, title:"Multiplier", te:add_param_mult_text }); 
}


/// The content within the param mult info
function add_param_mult_content(lay)
{
	let y = 0.5;
	let mar = 2;
	let dy = 4;
	let col_round = LLLBLUE, col_text = DDBLUE;
	
	let x = maximim_label_width() + 4;
	let del_x = lay.dx - 4.5;
	let del_dx = 1.3;

	for(let i = 0; i < model.param.length; i++){
		let par = model.param[i];
		if(par.type == "param factor"){
			lay.add_button({te:"Parameter factor", x:1, y:y, dx:lay.dx-3, dy:dy+1.5, col:col_round, col2:col_text, type:"CurvedOutline"});
			y += 1.4;
			
			let si = 1.0, fo = get_font(si,"","times");
			let xx = x-par.label_info.dx-0.7-1.8;
			lay.add_button({te:"Knots =", x:xx+1.8, y:y+0.3, dy:si, si:si, font:fo, type:"Text", col:BLACK});
			
			let spl = par.spline;
			let te = stringify(spl.knot);
					
			let w = wright;
				
			lay.add_button({te:te, x:xx+5.6, y:y+0., dx:w-xx-5.6, dy:1.6, type:"SplineKnots", font: get_font(1.1,"","times"), ac:"EditSplineKnots", i:i});
			
			y += 2;
		
			if(add_view_button(par,w-4,y,i,lay,model)) w -= 4.5;
		
			display_constant(i,x,y,lay,w);
			
			lay.add_button({x:del_x, y:y+0.2, dx:del_dx, dy:del_dx, type:"Delete", i:i, ac:"DeleteParamMult"});
			y += dy;
		}
	}
	
	if(y == 0.5){
		let si = 1;
		let col = BLACK;
		lay.add_button({te:"There are currently no parameter factors added to the model.", x:0.5, y:0, dx:lay.dx-5, dy:1.3, type:"Text", si:si, font:get_font(si), col:col}); 
	}		
}


/// Adds a screen allowing simulation parameters to be edited
function add_model_param_content(lay)
{
	let y = 0;
	let mar = 2;
	let x = mar;
	let dy_param = 2.5;
	let gap_param = 0.5;
	
	let depnum = -1;
	
	let param_cat = get_param_cat([],"");
	
	if(no_param(param_cat)) center_message("No parameters are in the model.",lay);
	
	for(let cat = 0; cat < param_cat.length; cat++){
		let pc = param_cat[cat];
		if(pc.list.length > 0){
			y = lay.add_subtitle(pc.name,1,y+0.2,WHITE,{te:pc.sim_te});
		
			x = mar;
			for(let k = 0; k < pc.list.length; k++){
				let i = pc.list[k];
				let par = model.param[i];
				let info = par.label_info;
			
				if(x+info.dx+1 >= lay.inner_dx){ x = mar; y += dy_param;}
		
				lay.add_button({x:x, y:y, dx:info.dx+1, dy:1.4*si_big, type:"ParamLabel2", info:info, col:BLACK, i:i, eqn_appear:par.eqn_appear, ac:"ParamInfo", back_col:WHITE});
				
				x += info.dx+1+gap_param;
			}
			x = mar; y += dy_param;
		}
	}		
	
	y += 1;
	
	let si_small = 0.8;
	let fo_small = get_font(si_small,"Bold");
	
	x = maximim_label_width() + 4;

	let del_x = lay.dx - 4.5;
	let del_dx = 1.3;
	
	let col_round = LLLBLUE, col_text = DDBLUE;
	
	let dy = 2;
	
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		
		let num = sp.ind_eff_group.length;
		if(num > 0){
			let h = 0;
			for(let i = 0; i < num; i++){
				h += dy; if(sp.ind_eff_group[i].A_matrix.check == true) h += dy; 
			}
			
			let tex = "Individual effects";
			if(model.species.length > 1) tex += " for "+sp.name;
			lay.add_button({te:tex, x:1, y:y, dx:lay.dx-3, dy:h+1.7, col:col_round, col2:col_text, type:"CurvedOutline"});
			y += 1.4;
			
			for(let i = 0; i < sp.ind_eff_group.length; i++){
				display_ind_eff_group(p,i,x,y,lay,wright);
				y += dy;
				
				let ieg = sp.ind_eff_group[i];
				if(ieg.A_matrix.check == true){
					
					let si = 0.9;
					
					let te_mat = "Load <e>"+ieg.A_matrix.name+"</e> matrix"; 
				
					if(ieg.A_matrix.loaded == true){ 
						te_mat = "Reload <e>"+ieg.A_matrix.name+"</e> matrix";
					}
					
					let xload = 4;
					if(te_mat != ""){
						let text_anno = text_convert_annotation(te_mat,si,si,100,"",BLACK);

						let dxload = text_anno.wmax+0.5;
					
						lay.add_button({word:text_anno.word, x:xload, y:y+0.1, dx:dxload, dy:1.1, type:"LinkPara", ac:"LoadAMatrix", p:p, i:i});
						xload += dxload+0.3;
		
						if(ieg.A_matrix.loaded == true){
							lay.add_button({te:"Edit", x:xload, y:y-0.05, dx:3.75, dy:1.2, ac:"EditAmatrix", type:"GreyView", p:p, i:i});
							xload += 3.75;
						}
						xload += 3;
					}
					
					y += dy;
				}
			}
			y += 1;
		}
	}
	
	// Fixed effects
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		
		if(sp.fix_eff.length > 0){
			let h = 2*dy*sp.fix_eff.length;
		
			let tex = "Fixed effects";
			if(model.species.length > 1) tex += " for "+sp.name;
			
			lay.add_button({te:tex, x:1, y:y, dx:lay.dx-3, dy:h+1.7, col:col_round, col2:col_text, type:"CurvedOutline"});
			y += 1.4;
			
			for(let i = 0; i < sp.fix_eff.length; i++){
				display_fix_eff(p,i,x,y,lay,wright);
				y += dy;
				
				let feg = sp.fix_eff[i];
				
				let si = 0.9;
				let te = "Load <e>X^"+feg.name+"</e> vector", dx = 7; 
				if(feg.X_vector.loaded == true){ te = "Reload <e>X^"+feg.name+"</e> vector"; dx = 8;}
					
				
				let text_anno = text_convert_annotation(te,si,si,100,"",BLACK);

				let xload = 4;
				let dxload = text_anno.wmax+0.5;
					
				lay.add_button({word:text_anno.word, x:xload, y:y+0.1, dx:dxload, dy:1.1, type:"LinkPara", ac:"LoadXvector", p:p, i:i});

				if(feg.X_vector.loaded == true){
					lay.add_button({te:"Edit", x:xload+dxload+0.3, y:y-0.1, dx:3.5, dy:1.2, ac:"EditXvector", type:"GreyView", p:p, i:i});
				}
				y += dy;
			}
			y += 1;
		}
	}

	// SPLINES
	{
		let num = 0;
		for(let i = 0; i < model.param.length; i++){
			let par = model.param[i];
			
			if(par.spline.on == true && par.type != "derive_param" && par.type != "param factor"){
				num++;
				if(model.param[i].spline.smooth.check == true) num++;
			}
		}

		if(num > 0){
			lay.add_button({te:"Splines", x:1, y:y, dx:lay.dx-3, dy:dy*num+1.5, col:col_round, col2:col_text, type:"CurvedOutline"});
			y += 1.4;
			
			for(let i = 0; i < model.param.length; i++){
				let par = model.param[i];
				let spl = par.spline;
				if(spl.on == true && par.type != "derive_param" && par.type != "param factor"){
					lay.display_param(x-par.label_info.dx-0.7,y-0.1,par.label_info);
				
					let si = 1.0, fo = get_font(si,"","times");
					lay.add_button({te:"Knots =", x:x+1.8, y:y+0.3, dy:si, si:si, font:fo, type:"Text", col:BLACK});
					
					let te = stringify(spl.knot);
					
					let w = wright;
					
					lay.add_checkbox(w-4,y+0.3,"Smooth","Smooth",spl.smooth,WHITE,{title:"Smoothing", te:smoothing_text});
					w -= 4;
				
					lay.add_button({te:te, x:x+5.6, y:y+0., dx:w-x-5.6, dy:1.6, type:"SplineKnots", font: get_font(1.1,"","times"), ac:"EditSplineKnots", i:i});
				
					y += dy;
					
					if(spl.smooth.check == true){
						lay.add_radio(x+1.6,y,"Log-Normal","Log-Normal",spl.smooth.type);
						lay.add_radio(x+8.6,y,"Normal","Normal",spl.smooth.type);
						
						let si = 0.9, fo = get_font(si,"","times");
						let sh = 0.15;
						lay.add_button({te:"Value =", x:x+16.6, y:y+0.03+sh, dy:si, si:si, font:fo, type:"Text", col:BLACK});
					
						lay.add_button({te:spl.smooth.value, x:x+19.4, y:y+sh, dx:8, dy:1.2, type:"SmoothValue", font: get_font(0.9), ac:"EditSplineKnots", i:i});
						y += dy;
					}
				}
			}
			y += 1;
		}
	}
	
	{ // CONSTS
		let num = 0;
		for(let i = 0; i < model.param.length; i++){
			let par = model.param[i];
			if(par.variety == "const" && par.type != "param factor") num++;
		}
				
		if(num > 0){
			lay.add_button({te:"Constants", x:1, y:y, dx:lay.dx-3, dy:dy*num+1.5, col:col_round, col2:col_text, type:"CurvedOutline"});
			y += 1.4;
		
			for(let i = 0; i < model.param.length; i++){
				let par = model.param[i];
				if(par.variety == "const" && par.type != "param factor"){
					let w = wright;
					
					if(par.value_desc == no_elements){
						display_no_element(par,x,y,lay,w);
					}
					else{
						if(add_view_button(par,w-4,y,i,lay,model)) w -= 4.5;
						
						display_constant(i,x,y,lay,w);
					}
					
					if(!par.dist_mat && !par.iden_mat){
						lay.add_button({x:del_x, y:y+0.2, dx:del_dx, dy:del_dx, type:"Delete", i:i, ac:"DeleteParamConst"});
					}
					
					y += dy;
				}
			}
			y += 1;
		}
	}
	
	{ // Factors
		let num = 0;
		for(let i = 0; i < model.param.length; i++){
			let par = model.param[i];
			if(par.factor) num++;
		}
				
		if(num > 0){
			lay.add_button({te:"Factors", x:1, y:y, dx:lay.dx-3, dy:dy*num+1.5, col:col_round, col2:col_text, type:"CurvedOutline"});
			y += 1.4;
		
			for(let i = 0; i < model.param.length; i++){
				let par = model.param[i];
				if(par.factor){
					let w = wright;
					display_factor(i,x,y,lay,w);
	
					lay.add_button({x:del_x, y:y+0.2, dx:del_dx, dy:del_dx, type:"Delete", i:i, ac:"DeleteParamFactor"});				
					y += dy;
				}
			}
			y += 1;
		}
	}
	
	{ // REPARAM
		let num = count(model.param,"variety","reparam");
		if(num > 0){
			lay.add_button({te:"Reparameterised", x:1, y:y, dx:lay.dx-3, dy:dy*num+1.5, col:col_round, col2:col_text, type:"CurvedOutline"});
			y += 1.4;

			for(let i = 0; i < model.param.length; i++){
				let par = model.param[i];
				if(par.variety == "reparam"){
					let w = wright;
					
					if(par.value_desc == no_elements && par.reparam_eqn_on != true){
						display_no_element(par,x,y,lay,w);
					}
					else{
						display_reparam(i,x,y,lay,w);
					}
					lay.add_button({x:del_x, y:y+0.2, dx:del_dx, dy:del_dx, type:"Delete", i:i, ac:"DeleteParamReparam"});
					y += dy;
				}
			}
			y += 1;
		}
	}
	
	{ // DIST
		let num = count(model.param,"variety","dist");
		if(num > 0){
			lay.add_button({te:"Distributions", x:1, y:y, dx:lay.dx-3, dy:dy*num+1.5, col:col_round, col2:col_text, type:"CurvedOutline"});
			y += 1.4;

			for(let i = 0; i < model.param.length; i++){
				let par = model.param[i];
				if(par.variety == "dist"){
					let w = wright;
					
					if(par.dep.length > 0){
						lay.add_checkbox(w-4,y+0.4,"Split","Split",par.prior_split_check,WHITE);
						w -= 4;
					}
		
					if(par.dep.length == 0 || par.prior_split_check.check == false){
						display_distribution(i,x,y,lay,true,true,w);
					}
					else{
						if(par.prior_split_desc == no_elements){
							display_no_element(par,x,y,lay,w);
						}
						else{
							display_distribution_split(i,x,y,lay,true,true,"dist",w);
						}
					}
					
					lay.add_button({x:del_x, y:y+0.2, dx:del_dx, dy:del_dx, type:"Delete", i:i, ac:"DeleteParamDist"});
					y += dy;
				}
			}
			y += 1;
		}
	}
	
	{ // DERIVE	
		let der = model.derive;
		let num = der.length;
		if(num > 0){
			lay.add_button({te:"Derived", x:1, y:y, dx:lay.dx-3, dy:dy*num+1.7, col:col_round, col2:col_text, type:"CurvedOutline"});
			y += 1.5;

			for(let i = 0; i < der.length; i++){
				let w = wright;
				
				display_derive(i,x,y,lay,true,w);
					
				lay.add_button({x:del_x, y:y+0.2, dx:del_dx, dy:del_dx, type:"Delete", val:i, ac:"DeleteDerive"});
				y += dy;
			}
			y += 1;
		}
	}
}


/// Bubble which allows user to select a constant parameter
function set_constant_bubble(cont)
{	
	if(check_param_free("const") == true){
		cont.dx = 20;
		bubble_addtitle(cont,"Set constant");

		bubble_addparagraph(cont,"Select parameter which is going to be set as a constant:",0,cont.dx);
		cont.y += 0.2;
	
		bubble_addscrollable(cont,{type:"const param sel", ymax:10, ac:"AddConstParam"}); 
		
		add_bubble_end(cont);
	}
	else{
		cont.dx = 10;
		bubble_addtitle(cont,"Set constant");

		bubble_addparagraph(cont,"There are currently no free parameters to set.",0,cont.dx); 
	}
}


/// Bubble which allows user to select a constant parameter
function set_factor_bubble(cont)
{	
	if(check_param_free("fac") == true){
		cont.dx = 20;
		bubble_addtitle(cont,"Set factor");

		bubble_addparagraph(cont,"Select parameter which is going to be set as a factor:",0,cont.dx);
		cont.y += 0.2;
	
		bubble_addscrollable(cont,{type:"fac param sel", ymax:10, ac:"AddFactorParam"}); 
		
		add_bubble_end(cont);
	}
	else{
		cont.dx = 10;
		bubble_addtitle(cont,"Set factor");

		bubble_addparagraph(cont,"There are currently no parameters to set as factors.",0,cont.dx); 
	}
}


/// Bubble which allows user to select a reparam parameter
function set_reparam_bubble(cont)
{	
	let bub = inter.bubble;
	if(bub.set_reparam_type){
		cont.dx = 15;
		bubble_addtitle(cont,"Set reparameterisation");
		bubble_addparagraph(cont,"Set via an equation or split into elements:",0,cont.dx); 
		bubble_addradio(cont,0,"equation","Equation",bub.radio);
		bubble_addradio(cont,0,"split","Split",bub.radio);	
		add_end_button(cont,"Done","AddReparamParam2");		
		return;
	}
	
	if(check_param_free("reparam") == true){
		cont.dx = 20;
		bubble_addtitle(cont,"Set reparameterisation");

		bubble_addparagraph(cont,"Select parameter which is going to be reparameterised:",0,cont.dx);
		cont.y += 0.2;
	
		bubble_addscrollable(cont,{type:"reparam param sel", ymax:10, ac:"AddReparamParam"});
		add_bubble_end(cont);	
	}
	else{
		cont.dx = 17;
		bubble_addtitle(cont,"Set reparameterisation");

		bubble_addparagraph(cont,"There are currently no free parameters to set.",0,cont.dx); 
	}
}


/// Determines if reparameterisation is done via an equation  
function set_reparam_bubble2(cont)
{	
	cont.dx = 20;
	bubble_addtitle(cont,"Set reparameterisation");
	bubble_addparagraph(cont,"Set via an equation .",0,cont.dx); 
}


/// Bubble which allows user to select a parameter for distribution
function set_distribution_bubble(cont)
{
	if(check_param_free("dist") == true){
		cont.dx = 18;
		bubble_addtitle(cont,"Set distribution");
		bubble_addparagraph(cont,"Select parameter to be sampled from a distribution:",0,cont.dx); cont.y += 0.2;
		
		bubble_addscrollable(cont,{type:"dist param sel", ymax:10, ac:"AddDistParam"}); 
		
		add_bubble_end(cont);
		//cont.y += 0.3;
		//let ac; if(inter.bubble.sel_param != "") ac = "AddDistParam";
		//add_end_button(cont,"Set",ac);
	}
	else{
		cont.dx = 10;
		bubble_addtitle(cont,"Set distribution");

		bubble_addparagraph(cont,"There are currently no free parameters to set.",0,cont.dx); 
	}
}


/// Bubble which allows user to select a parameter for modification 
function add_param_mult_bubble(cont)
{
	cont.dx = 20;
	bubble_addtitle(cont,"Set parameter");

	bubble_addparagraph(cont,"Select parameter which is going to be modified:",0,cont.dx);
	cont.y += 0.2;

	bubble_addscrollable(cont,{type:"trans param sel", ymax:10, ac:"AddParamMultParam"});
	add_bubble_end(cont);	
}


/// Bubble which allows user to select a parameter for distribution
function set_derived_bubble(cont)
{
	cont.dx = 25;
	bubble_addtitle(cont,"Set derived parameter",{te:setderived_text});
	bubble_addparagraph(cont,"Define the new parameter and how it relates to others:",0,cont.dx); 
	cont.y += 0.7;
	
	let x = 7, gap = 2;
	
	let y_st = cont.y;
	bubble_input(cont,"",{type:"derive_eqn1", eqn:true, w:x});
	let y_end = cont.y;
	
	cont.y = y_st;
	bubble_input(cont,"",{type:"derive_eqn2", eqn:true, x:x+gap, w:cont.dx-x-gap});
	let y_end2 = cont.y;
	
	cont.y = y_st;
	bubble_simple_text(cont,"=",1.2,BLACK,x+0.7); cont.y -= 1.2+0.3;
	
	if(y_end > y_end2) cont.y = y_end; else cont.y = y_end2;
	
	cont.y += 0.2;
	
	add_end_button(cont,"Add","AddDerived");
}


/// Bubble which allows user to select a parameter for distribution
function alter_derived_bubble(cont,bu)
{
	cont.dx = 25;
	bubble_addtitle(cont,"Modify derived parameter",{te:setderived_text});
	bubble_addparagraph(cont,"Change how the derived parameter relates to others:",0,cont.dx); 
	cont.y += 0.7;
	
	let x = 7, gap = 2;
	
	let y_st = cont.y;
	
	bubble_input(cont,"",{type:"derive_eqn1", eqn:true, w:x});
	let y_end = cont.y;
	
	cont.y = y_st;
	bubble_input(cont,"",{type:"derive_eqn2", eqn:true, x:x+gap, w:cont.dx-x-gap});
	let y_end2 = cont.y;
	
	cont.y = y_st;
	bubble_simple_text(cont,"=",1.2,BLACK,x+0.7); cont.y -= 1.2+0.3;
	
	if(y_end > y_end2) cont.y = y_end; else cont.y = y_end2;
	
	cont.y += 0.2;
	
	add_end_button(cont,"OK","DeriveOK");
}


/// Works out if there is a free parameter which can be selected
function check_param_free(op)
{
	for(let i = 0; i < model.param.length; i++){
		let par = model.param[i];
		if(param_pos(par,op)) return true;
	}
	return false;
}


/// Determines if possible to select parameter
function param_pos(par,op)
{
	if(par.type == "derive_param") return false;
		 
	switch(op){
	case "const":
		if(par.variety != "normal") return false;
		break;
		
	case "reparam":
		if(par.variety != "normal") return false;
		if(par.time_dep) return false;
		if(par.factor) return false;
		break;
		
	case "dist":
		if(par.variety != "normal") return false;
		if(par.factor) return false;
		if(par.dep.length == 0) return false;
		break;
	
	case "fac":
		if(par.factor) return false;
		if(par.time_dep) return false;
		if(par.variety != "normal" && par.variety != "const") return false;
		if(par.dep.length == 0) return false;
		break;
		
	case "trans":
		if(par.variety != "normal") return false;
		if(par.type != "trans_rate") return false;
		break;
	}

	return true;
}


/// Allows for a parameter to be selected
function param_sel_scrollable(lay,op)
{
	let cx = 0, cy = 0;
	let dy = 2.4;
	let gap = 0.5;
	
	for(let i = 0; i < model.param.length; i++){
		let par = model.param[i];

		if(param_pos(par,op)){	
			let info = par.label_info;
			if(cx+info.dx+1 >= lay.inner_dx){ cx = 0; cy += dy;}
		
			lay.add_button({x:cx, y:cy, dx:info.dx+1, dy:1.4*si_big, type:"ParamLabel", info:info, col:WHITE, back_col:BUBBLE_COL, i:i, ac:lay.op.ac});
		
			cx += info.dx+1+gap;
		}
	}

	cy += dy;
	
	return cy;
}


/// Displays a row allowing a derived variable to be set
function display_derive(j,x,y,lay,edit,w)
{
	let der = model.derive[j];
	let si = 1.5;
	let fo = get_font(si);

	let si_text = 1.1;
	let fo_text = get_font(si_text);

	let right_mar = lay.dx-5;
	let right_mar2 = lay.dx-8;
	let set_dx = 0.8;

	/*
	let par;
	if(der.eqn1.param.length == 1){
		let param = der.eqn1.param[0];

		let i = 0; while(i < model.param.length && par_same(param,model.param[i]) == false) i++;
		if(i < model.param.length) par = model.param[i];
	}
	*/
	let par;
	if(der.eqn1.param.length == 1){
		par = der.eqn1.param[0];
		model.get_label_info(par);
	}
	
	let ac_lab; if(edit == true) ac_lab = "EditDerive";

	lay.add_button({te_eqn1:der.eqn1.te, te_eqn2:der.eqn2.te, val:j, x:2, y:y-0.05, x_eq:x, dx:w-2, dy:1.7, si:1.4, par:par, type:"Derive", col:BLUE, ac:ac_lab});
}


/// Displays bubble which is selected on model parameter page
function set_paramselect_bubble(cont,eqn_appear)
{
	cont.dx = 17;
		
	bubble_addtitle(cont,"Parameter",{te:derived_text});
	
	let par = model.param[cont.bu.i];

	let te;
	let fl = false;

	switch(par.variety){
	case "normal":
		if(par.type == "fixed effect"){ te = "This parameter is a fixed effect."; fl = true;}
		if(par.type == "derive_param") te = "This parameter is derived."; 
		if(par.type == "variance" || par.type == "correlation"){
			let info = par.eqn_appear[0];
			let ieg = model.species[info.p].ind_eff_group[info.ieg];
			fl = true;

			if(par.type == "variance"){
				te = "This is the variance for individual effect <e>"+info.ie1+"</e>.";
			}
			else{
				te = "Correlation between individual effects <e>"+info.ie1+"</e> and <e>"+info.ie2+"</e>.";
			}
		}
		
		if(is_in_obsmodel(par)){
			te = "This parameter is in the observation model.";
			fl = true;
		}
		break;
		
	case "reparam": te = "This parameter is reparameterised."; break;
	case "const": te = "This parameter is set to a constant."; break;
	case "dist": te = "This parameter is sampled."; break;
	default: error("var option error"); break;
	}

	let te2 = "This parameter appears in:";
	if(te != undefined){
		cont.y += 0.2;
		if(fl == false){
			bubble_addbutton(cont,"View",cont.dx-3.6,cont.y-0.2,3.5,1.2,"GreyView","SelectParam",{par:par});
		}
		
		bubble_addparagraph(cont,te,0,cont.dx);
	
		cont.y += 0.5;
		te2 = "It appears in:";
	}

	if(par.type != "derive_param" && par.type != "variance" && par.type != "correlation"){
		bubble_addparagraph(cont,te2,0,cont.dx);
	
		cont.y += 0.2;
		bubble_addscrollable(cont,{type:"param details", eqn_appear:eqn_appear, ymax:10});
	}
	
	add_bubble_end(cont);	
}

		
/// Allows for a parameter
function param_details_scrollable(lay)
{
	let cx = 0, cy = 0.2;
	let dy = 2.4;
	let gap = 0.5;

	let dx = lay.inner_dx;
	
	let store = [];
	
	let eq_app = lay.op.eqn_appear;
	
	let too_big = false;
	
	for(let i = 0; i < eq_app.length; i++){
		let eq = eq_app[i];

		let te, ac;
	
		switch(eq.type){
		case "trans_mean": case "trans_rate": case "trans_shape":
		case "trans_scale": case "trans_cv":
			{
				let tr = model.species[eq.eqn_info.p].cla[eq.eqn_info.cl].tra[eq.eqn_info.i];
				let spl = tr.name.split("→");
				
				te = "Transition "+tr.name;		
				if(tr.variety == SOURCE) te = "Source for "+spl[1];
				if(tr.variety == SINK) te = "Sink for "+spl[0];
				ac = "SelectTransEqn";
			}
			break;
			
		case "trans_bp":
			{
				let tr = model.species[eq.eqn_info.p].cla[eq.eqn_info.cl].tra[eq.eqn_info.i];
			
				te = "Branching prob. for "+tr.name;
				ac = "SelectTransEqn";
			}
			break;
			
		case "Se": case "Sp":
			te = "Diagnostic test data";
			ac = "SelectDataSpec";
			break;
			
		case "mut_rate": case "seq_var":
			te = "Genetic data";
			ac = "SelectDataSpec";
			break;
			
		case "derived":
			te = "Derived quantity";
			ac = "SelectDerived";
			break;
		
		case "comp_prob":
			te = "Compartmental probability";
			ac = "SelectDataElement";
			break;
			
		case "trans_prob":
			te = "Transition observation prob.";
			ac = "SelectDataSpec";
			break;
		
			
		case "reparam":
			{
				te = "Reparameterisation";
				if(eq.eqn_info != undefined){
					let name = get_full_parameter_name(eq.eqn_info.par_name);
					te = "Reparameterisation of <e>"+name+"</e>";
					ac = "SelectReparam";
				}
			}
			break;
			
		case "prior": case "dist":
			{
				let name = get_full_parameter_name(eq.eqn_info.par_name);
				te = "Distribution for <e>"+name+"</e>";
				ac = "SelectDist";
			}
			break;
			
		case "derive_param":
			break;
		
		default: error("Option could not be recognised: "+eq.type); break;
		}
		
		if(find_in(store,te) == undefined){
			lay.add_button({te:"View", x:dx-3.8, y:cy-0.2, dx:3.5, dy:1.2, ac:ac, type:"GreyView", eqn_info:eq.eqn_info});
			cy = lay.add_paragraph("<b>"+te+"</b>",dx-4,cx,cy,BLACK,warn_si,warn_lh);
		
			cy += 0.5;
		
			store.push(te);
			if(store.length > 10){
				cy = lay.add_paragraph("+ more",dx-4,cx+0.3,cy-0.3,BLACK,warn_si,warn_lh);
				break;
			}
		}
	}
	
	return cy;
}


/// Allows for pedigree to be editted
function edit_A_pedigree(Amat,p,i)
{
	start_data_source("APed",{},{p:p,i:i});
	
	let so = edit_source;
	so.table_loaded = true;
	so.data_table_use = UNSET;

	let tab = so.table;
	
	for(let r = 0; r < Amat.ind_list.length; r++){
		let row = [];
		row[0] = Amat.ind_list[r];
		row[1] = Amat.sire_list[r];
		row[2] = Amat.dam_list[r];
		tab.ele.push(row);
	}
	tab.nrow = tab.ele.length;
	tab.ncol = 3;
	tab.edit = true;
	for(let j = 0; j < so.load_col.length; j++){
		tab.col_used.push(j);
		tab.heading.push(so.load_col[j].heading);
	}
	
	generate_screen();
}


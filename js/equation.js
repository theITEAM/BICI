"use strict";
// Function which incorporate equations into the model

// Different equation modes determine what properties the equation can have.
// These are now automatically generated when an equation is created
// "all" - Used generally for transition. No restrictions
// "param_only" - Allows only for a constant or single parameter 
// "derive" - Used in a derived formula. Can only involve existing parameters in the model or populations 
// "derive_param" - A derived parameter (cannot already be in the model). Must just be a parameter (with potentially some dependencies which agree with 'derive' expression).

// Setting positive:true means if equation is a number it must be positive
// Setting zero_one_range:true means if equation is a number it must be in range 0 - 1

/// Adds background to edit equation page
function add_equation_background_buts(lay)
{
	lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, type:"EquationBackground"});
	lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, ac:"Nothing", type:"Nothing"});
}


/// Adds buttons to edit equation page
function add_equation_buts(lay)
{
	let marx = 1;

	let dx = lay.dx;
	
	let martop = 2.9, marbot = 2.8;
	let ytop = 5-martop;
	let dy = martop+marbot+11.4;
	
	let warn = inter.equation.warning;

	if(warn != undefined){
		warn.dx = dx-3; 
		let text_anno = text_convert_annotation(warn.te,warn_lh,warn_lh,warn.dx,"",BLACK);
		dy += text_anno.height;
	}

	lay.add_button({x:0, y:ytop, dx:dx, dy:dy, type:"EquationBack"});

	lay.add_title("Equation editor",1.3,ytop + 0.6,{col:WHITE, col_line:WHITE,  te:eqn_text, back_col:GREY_BACK});
	
	let mode = inter.equation.mode;
	
	equation_calulator(lay,1,ytop+3,lay.inner_dx-2,inter.equation.source,warn,mode);
	
	let marrx = 1.3,  marry = 0.6, busi = 1.2;
	lay.add_button({x:dx-marx-marrx, y:ytop+marry, dx:busi, dy:busi, ac:"CloseEquation", type:"WhiteCross"});

	let widbut = but_width, heibut = but_height;

	lay.add_corner_button([["Done","White","EquationDone"]],{x:dx-marx+0.3, y:ytop+dy-0.2});

	let toolbar = inter.equation.toolbar;

	if(toolbar != undefined) add_toolbar(toolbar,ytop+dy+0.7,lay);
}


/// Adds the tool bar to the edit equartion 
function add_toolbar(te,ytop,lay)
{
	let marx = 1;
	let dx = lay.dx;
	
	let max = 11;
	let shift;
	switch(te){
		case "Add model parameter": shift = 2.7; max += 2; break;
		case "Add parameter": shift = 4.7; break;
		case "Add tensor": shift = 4.7; break;
		case "Add sum": shift = 2.7; break;
		case "Add population": shift = 2.7; break;
		case "Add individual effect": shift = 4.7; break;
		case "Add fixed effect": shift = 4.7; break;
		case "Add time": shift = 2.7; break;
		case "Add distance": shift = 2.7; break;
		default: error("Option not recognised 30"); break;
	}
	
	let l_store = inter.layer.length;
	add_layer("EquationAddQuantity",lay.x+marx,ytop+shift,dx-2*marx-0.2,0,{te:te, ymax:max-shift});	
	let hei = inter.layer[l_store].dy;

	let dy = shift+hei+2.3;

	lay.add_button({te:te, x:0, y:ytop, dx:dx, dy:dy, type:"EquationAddBack"});

	let marrx = 0.9,  marry = 0.8, busi = 1;
	lay.add_button({x:dx-marx-marrx, y:ytop+marry, dx:busi, dy:busi, ac:"CloseToolbar", type:"WhiteCross"});

	let ac;
	switch(te){
	case "Add model parameter":
		{
		}
		break;
		
	case "Add parameter":
		{
			let bdx = 1.59, bdy = 1.8;
			
			let x = 0.8, y = ytop+2.6;
			for(let i = 0; i < greek.length; i++){
				lay.add_button({te:greek[i], x:x, y:y, dx:bdx, dy:bdy, ac:"Letter_But", type:"Letter_But", col:BLACK, colov:DDGREY});
				x += bdx+0.2;
			}		
			
			ac = "AddParameterDone";	
		}
		break;
	
	case "Add tensor":
		{
			let bdx = 1.59, bdy = 1.8;
			
			let x = 0.8, y = ytop+2.6;
			for(let i = 0; i < greek_capital.length; i++){
				lay.add_button({te:greek_capital[i], x:x, y:y, dx:bdx, dy:bdy, ac:"Letter_But", type:"Letter_But", col:BLACK, colov:DDGREY});
				x += bdx+0.2;
			}		
		
			ac = "AddTensorDone";	
		}
		break;
		
	case "Add sum":
		ac = "AddSumDone";	
		break;
		
	case "Add population":
		ac = "AddPopDone";	
		break;
		
	case "Add individual effect":
		{
			let x = 0.5, y = ytop+2.6;
			let bdx = 1.47, bdy = 1.8;
			for(let i = 0; i < alphabet.length; i++){
				lay.add_button({te:alphabet[i], x:x, y:y, dx:bdx, dy:bdy, ac:"Letter_But", type:"Letter_But", col:BLACK, colov:DDGREY});
				x += bdx+0.2;
			}			
			
			ac = "AddIndEffectDone";	
		}
		break;
		
	case "Add fixed effect":
		{
			let x = 0.5, y = ytop+2.6;
			let bdx = 1.47, bdy = 1.8;
			for(let i = 0; i < alphabet.length; i++){
				lay.add_button({te:alphabet[i], x:x, y:y, dx:bdx, dy:bdy, ac:"Letter_But", type:"Letter_But", col:BLACK, colov:DDGREY});
				x += bdx+0.2;
			}			
			
			ac = "AddFixEffectDone";	
		}
		break;
		
		case "Add time":
		{
			let x = 1, y = ytop+2.6;
			lay.add_paragraph("Adds the time variable to the equation.",lay.dx,x,y,WHITE,para_si,para_lh);
			ac = "AddTimeDone";	
		}
		break;
		
	case "Add distance":
		ac = "AddDistanceDone";	
		break;

	default: error("Option not recognised 31"); break;
	}
	
	let warn = inter.equation.toolbar_warning;
	if(warn != undefined && warn != ""){
		lay.add_paragraph(warn,dx-marx-4,1,ytop+dy-1.8,LRED,para_si,para_lh);
	}
	
	lay.add_button({te:"Add", x:dx-marx-3.5, y:ytop+dy-2, dx:but_width_small, dy:but_height_small, ac:ac, type:"SmallButton"});
}
	

/// Buttons for adding a quantity (e.g. parameter, variable etc..)	
function add_quantity_content_buts(lay)
{	
	let eqn = inter.equation;
	
	let x = 0, y = 0, dy = 1;

	switch(lay.op.te){
	case "Add model parameter":
		{
			let dy = 2.4;
			let gap = 0.5;
	
			for(let i = 0; i < model.param.length; i++){
				let par = model.param[i];
				let info = par.label_info;
				if(x+info.dx+0.5 >= lay.inner_dx){ x = 0; y += dy;}

				lay.add_button({x:x, y:y, dx:info.dx+0.5, dy:1.4*si_big, type:"ParamLabel", info:info, col:WHITE, back_col:BLACK, i:i, ac:"AddModParamEqn"});
				
				x += info.dx+0.5+gap;
			}
			x = 0; y += dy;
		}
		break;
		
	case "Add parameter":
		{		
			let marleft = 6;
			
			lay.add_button({te:"VARIATION:", x:0, y:y, dx: marleft, dy:dy, type:"Text", col:WHITE, si:0.8, font:get_font(0.8)});
				
			lay.add_radio_white(marleft,y,"None","None",eqn.dep_radio);
			lay.add_radio_white(marleft+5,y,"Time","Time",eqn.dep_radio);
			
			y += 1.6;
			
			switch(eqn.mode){
			case "derive_param":	
				{	
					lay.add_button({te:"INDICES:", x:0, y:y, dx: marleft, dy:dy, type:"Text", col:WHITE, si:0.8, font:get_font(0.8)});
					let gap = 3.4;
					let ml = 4;
					
					for(let i = -1; i < index_max; i++){
						lay.add_radio_white(ml+(i+1)*gap,y,i+1,i+1,eqn.num_radio);
					}
					y += 1.6;
					
					for(let i = 0; i < eqn.num_radio.value; i++){
						lay.add_dropdown(ml+(i+1)*gap,y,2.9,10,eqn.index_drop[i],eqn.pos,"white_style");	
					}
					y += 1;
				}
				break;
				
			case "all":
				{
					let p = eqn.origin.p;
					let cl_view = eqn.origin.cl;
			
					let sp = model.species[p];
					let ncla = sp.cla.length;
					if(ncla > 1){
						lay.add_button({te:"DEPENDENCY:", x:x, y:y, dx: marleft, dy:dy, type:"Text", col:WHITE, si:0.8, font:get_font(0.8)});
					
						x = marleft;
						
						for(let cl = 0; cl < ncla; cl++){
							if(cl != cl_view){
								let name = sp.cla[cl].name;
						
								let w = 1.6+text_width(name,get_font(si_radio))+1;
								if(x+w > lay.inner_dx){ x = marleft; y += 1.4;}
						
								lay.add_checkbox_white(x,y,name,name,eqn.cla_check[cl]);
								x += w;
							}
						}
						
						y += 1.2;
					}
				}
				break;
			}
		}
		break;
		
		
	case "Add tensor": 
		{	
			let p = eqn.origin.p;
			let cl_view = eqn.origin.cl;
			
			let marleft = 6;
			
			lay.add_button({te:"VARIATION:", x:0, y:y, dx: marleft, dy:dy, type:"Text", col:WHITE, si:0.8, font:get_font(0.8)});
			
			lay.add_radio_white(marleft,y,"None","None",eqn.dep_radio);
			lay.add_radio_white(marleft+5,y,"Time","Time",eqn.dep_radio);
			y += 1.6;
			
			let ml = 4;
			
			lay.add_button({te:"INDICES:", x:0, y:y, dx: marleft, dy:dy, type:"Text", col:WHITE, si:0.8, font:get_font(0.8)});
			let gap = 3.4;
			for(let i = 0; i < index_max; i++){
				lay.add_radio_white(ml+i*gap,y,i+1,i+1,eqn.num_radio);
			}
			y += 1.6;
			
			for(let i = 0; i < eqn.num_radio.value; i++){
				lay.add_dropdown(ml+i*gap,y,2.9,10,eqn.index_drop[i],eqn.pos,"white_style");	
			}
			y += 1;			
		}
		break;
		
	case "Add sum": 
		{	
			let p = eqn.origin.p;
			let cl_view = eqn.origin.cl;
			
			let marleft = 6;
			
			let ml = 4;
			
			lay.add_button({te:"INDICES:", x:0, y:y, dx: marleft, dy:dy, type:"Text", col:WHITE, si:0.8, font:get_font(0.8)});
			let gap = 3.4;
			for(let i = 0; i < index_max; i++){
				lay.add_radio_white(ml+i*gap,y,i+1,i+1,eqn.num_radio);
			}
			y += 1.6;
			
			for(let i = 0; i < eqn.num_radio.value; i++){
				lay.add_dropdown(ml+i*gap,y,2.9,10,eqn.index_drop[i],eqn.pos,"white_style");	
			}
			y += 1;			
		}
		break;
		
	case "Add population":
		{
			x = 0;
			lay.add_button({te:"Select species:", x:x, y:y, dx:lay.dx-2, dy:dy, type:"Text", col:WHITE, si:0.8, font:get_font(0.8,"bold")});
			y += 1.4; 
			
			for(let p = 0; p < model.species.length; p++){	
				let name = model.species[p].name;
				
				let w = 1.6+text_width(name,get_font(si_radio))+1;
				if(x+w > lay.inner_dx){ x = marleft; y += 1.4;}
				lay.add_radio_white(x,y,name,name,eqn.sp_radio);
				x += w;
			}
			y += 2.1; 
			
			
			let p = find(model.species,"name",eqn.sp_radio.value);
			
			x = 0;
			let sp = model.species[p];
			
			for(let cl = 0; cl < sp.ncla; cl++){
				let tree = eqn.button_tree[p][cl];
				let claa = sp.cla[cl];
				let ind = claa.index;
				
				let name = claa.name;
				let te = name.toUpperCase()+":";
				let font = get_font(0.8);
				let left = 0.;
				
				lay.add_button({te:te, x:0, y:y, dx:lay.dx-2, dy:dy, type:"Text", col:WHITE, si:0.8, font:font});

				x = left; y += 1.4;
				lay.add_radio_white(x,y,"All","All",tree.sel_radio);

				lay.add_radio_white(x+5,y,"Compartment","Compartment",tree.sel_radio);
			
				lay.add_radio_white(x+14,y,"Index","Index",tree.sel_radio);
				y += 1.4;
				
				let marleft = left+2;
				x = marleft;
				
				switch(tree.sel_radio.value){
				case "Compartment":
					for(let c = 0; c < claa.ncomp; c++){
						let name = claa.comp[c].name;
					
						let w = 1.6+text_width(name,get_font(si_radio))+1;
						if(x+w > lay.inner_dx){ x = marleft; y += 1.4;}
						lay.add_checkbox_white(x,y,name,name,tree.comp_check[c]);
						x += w;
					}
					y += 1.4;
					break;
				
				case "Index":
					let fo = get_font(si_radio,"italic","times"); 
					lay.add_radio_white(marleft,y,ind,ind,tree.index_radio,{fo:fo});
					lay.add_radio_white(marleft+3,y,ind+"'",ind+"'",tree.index_radio,{fo:fo});
					lay.add_radio_white(marleft+6,y,ind+"''",ind+"''",tree.index_radio,{fo:fo});
					y += 1.4;
					break;
				}
			
				y += 0.3;
			}
		}
		break;

	case "Add individual effect": break;
	case "Add fixed effect": break;
	case "Add time": break;
	case "Add distance":
		{
			
			let marleft = 7;
			
			let p = eqn.origin.p;
			let sp = model.species[p];
			let ncla = sp.cla.length;
			lay.add_button({te:"CLASSIFICATION:", x:x, y:y+0.1, dx: marleft, dy:dy, type:"Text", col:WHITE, si:0.8, font:get_font(0.8)});
					
			x = marleft;
						
			for(let cl = 0; cl < ncla; cl++){
				let name = sp.cla[cl].name;
		
				let w = 1.6+text_width(name,get_font(si_radio))+1;
				if(x+w > lay.inner_dx){ x = marleft; y += 1.4;}
		
				lay.add_radio_white(x,y,name,name,eqn.cla_radio);
				x += w;
			}
			
			y += 1.6;
		}
		break;

	default: error("Option not recognised 32"); break;
	}
	
	if(y > lay.op.ymax) y = lay.op.ymax;
	lay.dy = y;
}


/// Adds buttons for the equation editor
function equation_calulator(lay,cx,cy,width,source,warn,mode)
{
	let nrow = 3;

	source.warn = warn;
	lay.add_textbox(cx+0.5,cy,width-1,nrow,source);
	
	cy += nrow*equation_linesi+0.8;
	
	if(warn != undefined){
		cy = lay.add_paragraph(warn.te,warn.dx,1.5,cy-0.4,LRED,warn_si,warn_lh);
		cy += 0.1;
	}

	let dx = 2.7, dy=1.5, gap = 0.5;
	let dx2 = 1.52, gap2 = 0.5;
	let dx3 = 1.52, gap3 = 0.5;
	
	let x = cx+0.2;
	for(let i = 0; i < functi.length; i++){
		let ddx = 2.95; if(functi_dx[i] != "def") ddx = functi_dx[i];
		lay.add_button({te:functi[i], x:x, y:cy, dx:ddx, dy:dy, ac:"Calculator_But", type:"Calculator_But", col:GREY, colov:LGREY});
		x += ddx+gap;
	}			
	
	cy += dy+0.4;
		
	x = cx+0.2;
	for(let i = 0; i < numbut.length; i++){
		lay.add_button({te:numbut[i], x:x, y:cy, dx:dx3, dy:dy, ac:"Calculator_But", type:"Calculator_But", col:"#666666", colov:DGREY});
		x += dx3+gap3;
	}

	for(let i = 0; i < opbut.length; i++){
		lay.add_button({te:opbut[i], x:x, y:cy, dx:dx2, dy:dy, ac:"Calculator_But", type:"Calculator_But", col:DGREY, colov:GREY});
		x += dx2+gap2;
	}
	cy += dy+0.4;
	
	cy += 0.2;
	x = cx+0.2;
	let font = get_font(si_toolbut);	
	let pad = 1.1;
	let gap4 = 0.5;
	
	let model_param_fl = false, param_fl = false, pop_fl = false;
	let indeff_fl = false, sum_fl = false, tensor_fl = false;
	
	let eqn = inter.equation;

	switch(eqn.mode){
	case "all": 
		param_fl = true; pop_fl = true; indeff_fl = true; sum_fl = true, tensor_fl = true;
		break;
		
	case "derive_param": 
		param_fl = true; 
		break;
	
	case "derive": 
		model_param_fl = true; pop_fl = true; sum_fl = true; 
		break;
	
	case "param only": 
		param_fl = true; 
		break; 
		
	case "param with dep":
		param_fl = true; 
		break;
	
	default:
		error("option prob 100"+eqn.mode); 
		break;
	}
	
	let ac = "AddToolbar";	
	
	if(model_param_fl == true){
		let te = "Model parameters", w = text_width(te,font)+pad;
		lay.add_button({te:te, te2:"Add model parameter", x:x, y:cy, dx:w, dy:dy, ac:ac, type:"Toolbar_But"});
		x += w+gap4;
	}
	
	if(param_fl == true){
		let te = "Parameter", w = text_width(te,font)+pad;
		lay.add_button({te:te, te2:"Add parameter", x:x, y:cy, dx:w, dy:dy, ac:ac, type:"Toolbar_But"});
		x += w+gap4;
	}
	
	if(pop_fl == true){
		let te = "Population", w = text_width(te,font)+pad;
		lay.add_button({te:te, te2:"Add population", x:x, y:cy, dx:w, dy:dy, ac:ac, type:"Toolbar_But"});
		x += w+gap4;
	}
	
	if(param_fl == true){
		let te = "t", w = text_width(te,font)+pad;
		lay.add_button({te:te, te2:"Add time", x:x, y:cy, dx:w, dy:dy, ac:ac, type:"Toolbar_But"});
		x += w+gap4;
	}
	
	if(param_fl == true){
		let te = "D", w = text_width(te,font)+pad;
		lay.add_button({te:te, te2:"Add distance", x:x, y:cy, dx:w, dy:dy, ac:ac, type:"Toolbar_But"});
		x += w+gap4;
	}
	
	if(indeff_fl == true){
		let te = "Ind. Effect", w = text_width(te,font)+pad;
		lay.add_button({te:te, te2:"Add individual effect", x:x, y:cy, dx:w, dy:dy, ac:ac, type:"Toolbar_But"});
		x += w+gap4;
	}
	
	if(indeff_fl == true){
		let te = "Fixed Effect", w = text_width(te,font)+pad;
		lay.add_button({te:te, te2:"Add fixed effect", x:x, y:cy, dx:w, dy:dy, ac:ac, type:"Toolbar_But"});
		x += w+gap4;
	}
	
	if(sum_fl == true){
		let te = "Sum", w = text_width(te,font)+pad;
		lay.add_button({te:te, te2:"Add sum", x:x, y:cy, dx:w, dy:dy, ac:ac, type:"Toolbar_But"});
		x += w+gap4;
	}
	
	if(tensor_fl == true){
		let te = "Tensor", w = text_width(te,font)+pad;
		lay.add_button({te:te, te2:"Add tensor", x:x, y:cy, dx:w, dy:dy, ac:ac, type:"Toolbar_But"});
		x += w+gap4;
	}
	
	cy += dy+0.6;
	
	return cy;
}


/// Initialises all the variables on a tool bar
function initialise_toolbar(type)
{
	let eqn = inter.equation;
	eqn.letter="";
	
	eqn.toolbar = type;
	
	switch(type){
	case "Add model parameter":
		break;
		
	case "Add parameter":
		{
			eqn.dep_radio = { value:"None"};

			if(eqn.mode == "derive_param"){
				eqn.num_radio = { value:0};
				let pos = [];
				for(let p = 0; p < model.species.length; p++){
					let sp = model.species[p];
					for(let cl = 0; cl < sp.ncla; cl++){
						let claa = sp.cla[cl];
						pos.push({te:claa.index});
					}
				}
				eqn.pos = pos;
				eqn.index_drop=[];
				for(let loop = 0; loop < index_max; loop++){
					eqn.index_drop[loop] = copy(pos[0]);
				}
			}
			
			if(eqn.mode == "all"){
				let p = eqn.origin.p;
				
				eqn.cla_check = [];
				let sp = model.species[p];
				for(let cl = 0; cl < sp.ncla; cl++){
					eqn.cla_check.push({check:false});
				}
			}
		}
		break;
		
	case "Add tensor":
		{
			eqn.dep_radio = { value:"None"};
			
			eqn.num_radio = { value:1};
			let pos = [];
			for(let p = 0; p < model.species.length; p++){
				let sp = model.species[p];
				for(let cl = 0; cl < sp.ncla; cl++){
					let claa = sp.cla[cl];
					let val = claa.index;
					for(let loop = 0; loop < 3; loop++){
						pos.push({te:val});
						val += "'";
					}
				}
			}
			eqn.pos = pos;
			eqn.index_drop=[];
			for(let loop = 0; loop < index_max; loop++){
				eqn.index_drop[loop] = copy(pos[0]);
			}
		}
		break;
		
	case "Add sum":
		{
			eqn.num_radio = { value:1};
			let pos = [];
			for(let p = 0; p < model.species.length; p++){
				let sp = model.species[p];
				for(let cl = 0; cl < sp.ncla; cl++){
					let claa = sp.cla[cl];
					let val = claa.index+"'";
					for(let loop = 0; loop < 2; loop++){
						pos.push({te:val});
						val += "'";
					}
				}
			}
			eqn.pos = pos;
			eqn.index_drop=[];
			for(let loop = 0; loop < index_max; loop++){
				eqn.index_drop[loop] = copy(pos[0]);
			}
		}
		break;
		
	case "Add population":		
		{
			eqn.sp_radio = {value:model.species[0].name};
			
			eqn.button_tree = [];
			let tree = eqn.button_tree;
			for(let p = 0; p < model.species.length; p++){
				let sp = model.species[p];
				tree[p] = [];
				for(let cl = 0; cl < sp.ncla; cl++){
					let claa = sp.cla[cl];
					tree[p][cl] = { sel_radio:{value:"All"}, index_radio:{value:claa.index}, comp_check:[]};
					for(let c = 0; c < claa.ncomp; c++){
						tree[p][cl].comp_check.push({check:false});
					}
				}
			}
		}
		break;
		
	case "Add individual effect": break;
	case "Add fixed effect": break;
	case "Add time": break;
	case "Add distance": 
		{
			let p = eqn.origin.p;	
			let sp = model.species[p];
			eqn.cla_radio = {value:sp.cla[0].name};
		}
		break;

	default: error("Option not recognised 33"+type); break;
	}
}


/// Starts the editing of an equation
function start_equation(te,eqn,origin_source,i)
{
	turn_off_cursor();
	
	inter.equation = { te:te, mode:eqn.mode, origin:origin_source, source:{ref:"equation", type:"equation"}};
					
	let j = 0; 
	while(j < inter.textbox_store.length){
		if(inter.textbox_store[j].source.type == "equation") inter.textbox_store.splice(j,1);
		else j++;
	}

	inter.textbox_store.push({ref: "equation", source: {type:"equation"}, te:te, eqn:eqn, store_text:[]});	
	
	generate_screen();
	set_focus(inter.textbox_store.length-1,i); 
	inter.over = {};
}


/// Add parameters to an equation
function equation_add_parameter()
{
	let eqn = inter.equation;
	eqn.toolbar_warning = "";
	
	if(eqn.letter == undefined || eqn.letter == ""){
		eqn.toolbar_warning = "A greek letter must be selected";
		return;
	}
	
	if(inter.cursor.i == undefined){
		eqn.toolbar_warning = "The cursor location must be selected"; 
		return;
	}
	
	let paste = eqn.letter;
	let dep = "";
	
	switch(eqn.mode){
	case "derive_param":
		{
			for(let i = 0; i < eqn.num_radio.value; i++){
				if(dep != "") dep += ",";
				dep += eqn.index_drop[i].te;
			}
		}
		break;
		
	case "all":
		{
			let p = eqn.origin.p;
			if(p != undefined){
				let cl_view = eqn.origin.cl;
					
				let sp = model.species[p];
				for(let cl = 0; cl < sp.ncla; cl++){
					if(cl != cl_view && eqn.cla_check[cl].check == true){
						if(dep != "") dep += ",";
						dep += sp.cla[cl].index;
					}
				}
			}
		}
		break;
	}
	
	if(dep != "") paste += "_"+dep;
		
	switch(inter.equation.dep_radio.value){
	case "Time": paste += "(t)"; break;
	case "None": break;
	default: error("Option not recognised 35"); break;
	}

	cursor_paste_with_space(paste);
	delete inter.equation.toolbar;
}


/// Add model parameter to an equation
function equation_add_modelparam(i)
{
	let par = model.param[i];
	let paste = remove_eq_quote(par.full_name);
	cursor_paste_with_space(paste);
	delete inter.equation.toolbar;
}


/// Add tensor to an equation
function equation_add_tensor()
{
	let eqn = inter.equation;
	eqn.toolbar_warning = "";
	
	if(eqn.letter == undefined || eqn.letter == ""){
		eqn.toolbar_warning = "A greek letter must be selected";
		return;
	}
	
	if(inter.cursor.i == undefined){
		eqn.toolbar_warning = "The cursor location must be selected"; 
		return;
	}
	
	for(let i = 0; i < eqn.num_radio.value; i++){
		for(let j = i+1; j < eqn.num_radio.value; j++){
			if(eqn.index_drop[i].te == eqn.index_drop[j].te){
				eqn.toolbar_warning = "The index '"+eqn.index_drop[i].te+"' is repeated more than once."; 
				return;
			}
		}
	}

	let paste = eqn.letter;
	
	let dep = "";
	for(let i = 0; i < eqn.num_radio.value; i++){
		if(dep != "") dep += ",";
		dep += eqn.index_drop[i].te;
	}
	if(dep != "") paste += "_"+dep;
		
	switch(inter.equation.dep_radio.value){
	case "Time": paste += "(t)"; break;
	default: error("Option not recognised 36"); break;
	}

	cursor_paste_with_space(paste);
	delete inter.equation.toolbar;
}


/// Add sum to an equation
function equation_add_sum()
{
	let eqn = inter.equation;
	eqn.toolbar_warning = "";
	
	if(inter.cursor.i == undefined){
		eqn.toolbar_warning = "The cursor location must be selected"; 
		return;
	}
	
	for(let i = 0; i < eqn.num_radio.value; i++){
		for(let j = i+1; j < eqn.num_radio.value; j++){
			if(eqn.index_drop[i].te == eqn.index_drop[j].te){
				eqn.toolbar_warning = "The index '"+eqn.index_drop[i].te+"' is repeated more than once."; 
				return;
			}
		}
	}
	
	let paste = "Σ"+eqn.letter;
	
	let dep = "";
	for(let i = 0; i < eqn.num_radio.value; i++){
		if(dep != "") dep += ",";
		dep += eqn.index_drop[i].te;
	}
	if(dep != "") paste += "_"+dep;
		
	paste += "()";
	
	cursor_paste_with_space(paste,-1);
	
	delete inter.equation.toolbar;
}


/// Pastes with some space padding
function cursor_paste_with_space(paste,shift)
{
	let i = inter.cursor.i;
	
	let st = inter.layer[inter.cursor.l].get_text_from_source();
	if(i > 0 && st.substr(i-1,1) != " ") paste = " "+paste;
	if(i < st.length && st.substr(i,1) != " "){ paste = paste + " "; if(shift != undefined) shift--;}

	cursor_paste(paste);
	if(shift != undefined) inter.cursor.i += shift;
}


/// Add population to an equation
function equation_add_population()
{
	let eqn = inter.equation;
	eqn.toolbar_warning = "";
	
	if(inter.cursor.i == undefined){
		eqn.toolbar_warning = "The cursor location must be selected"; 
		return;
	}
	
	let paste = "{";
	
	let p = eqn.origin.p;
	let psel = find(model.species,"name",eqn.sp_radio.value);
	if(p != psel) paste += model.species[psel].name+":";
	
	let start = paste;

	let sp = model.species[psel];
	for(let cl = 0; cl < sp.ncla; cl++){
		let claa = sp.cla[cl];
		
		let tree = eqn.button_tree[psel][cl];
	
		switch(tree.sel_radio.value){
		case "All":
			break;
			
		case "Compartment":
			if(paste != start) paste += ",";
			
			let co = "";
			for(let c = 0; c < claa.ncomp; c++){
				if(tree.comp_check[c].check == true){
					if(co != "") co += "|";
					co += claa.comp[c].name;
				}
			}
			paste += co; 
			
			if(co == ""){
				eqn.toolbar_warning = "'"+ claa.name+"' must have at least one compartment set"; 
				return;
			}
			break;
			
		case "Index":
			if(paste != start) paste += ",";
			paste += tree.index_radio.value;
			break;

		default: error("Option not recognised 37"); break;
		}
	}
	if(paste == start) paste += "All";
	
	paste += "}";
	cursor_paste(paste);
	delete inter.equation.toolbar;
}


/// Add individual effect to an equation
function equation_add_indeffect()
{
	let eqn = inter.equation;
	eqn.toolbar_warning = "";
	
	if(eqn.letter == undefined){
		eqn.toolbar_warning = "A letter must be selected";
		return;
	}
	
	if(inter.cursor.i == undefined){
		eqn.toolbar_warning = "The cursor location must be selected"; 
		return;
	}
		
	let paste = "["+eqn.letter+"]";
	cursor_paste(paste);
	delete inter.equation.toolbar;
}


/// Add fixed effect to an equation
function equation_add_fixeffect()
{
	let eqn = inter.equation;
	eqn.toolbar_warning = "";
	
	if(eqn.letter == undefined){
		eqn.toolbar_warning = "A letter must be selected";
		return;
	}
	
	if(inter.cursor.i == undefined){
		eqn.toolbar_warning = "The cursor location must be selected"; 
		return;
	}
	
	let paste = "<"+eqn.letter+">";
	//cursor_paste_with_space(paste);
	cursor_paste(paste);
	delete inter.equation.toolbar;
}


/// Add time to an equation
function equation_add_time()
{
	let eqn = inter.equation;
	eqn.toolbar_warning = "";

	if(inter.cursor.i == undefined){
		eqn.toolbar_warning = "The cursor location must be selected"; 
		return;
	}
	
	let paste = "t";
	cursor_paste_with_space(paste);
	delete inter.equation.toolbar;
}


/// Add distance to an equation
function equation_add_distance()
{
	let eqn = inter.equation;
	eqn.toolbar_warning = "";

	if(inter.cursor.i == undefined){
		eqn.toolbar_warning = "The cursor location must be selected"; 
		return;
	}
	
	let val = eqn.cla_radio.value;
	let p = eqn.origin.p;
	let sp = model.species[p];
	let cl = 0; while(cl < sp.cla.length && sp.cla[cl].name != val) cl++;
	if(cl == sp.cla.length) error("Could not find index");
	let index = sp.cla[cl].index;
	
	let paste = "D_"+index+","+index+"′";
	
	cursor_paste_with_space(paste);
	delete inter.equation.toolbar;
}


/// Run when the equation done button is clicked
function equation_done()
{
	let sto = inter.textbox_store[inter.textbox_store.length-1];

	let te;
	if(sto.ref != "equation"){ error("Should be equation"); return;}
	else te = sto.te;
	
	let eqn = sto.eqn;
	
	eqn.te = te;
	extract_equation_properties(eqn);

	// Looks for warnings associated with calculating equations for simulating data
	if(eqn.warn.length == 0 && subtab_name() == "Generate Data"){  
		switch(eqn.type){
		case "trans_prob": case "comp_prob":
			{
				let eq_calc = setup_eqn(te);
				if(eq_calc.err){ 
					eqn.warn.push({te:eq_calc.msg, cur:0, len:0});
				}
			}
			break;
		}
	}
	
	if(eqn.warn.length > 0){
		inter.equation.warning = eqn.warn[0];
		return;
	}
	
	turn_off_cursor();
	
	let orig = inter.equation.origin;
	let ref = orig.ref;
	
	if(ref == "ProbEqn"){
		if(te != orig.te_store) update_param();
		inter.equation = {};
		inter.textbox_store.pop();
		inter.bubble.check_warning = undefined;
	}
	else{
		let i = 0; while(i < inter.textbox_store.length && inter.textbox_store[i].ref != ref) i++;
		if(i == inter.textbox_store.length) error("Could not find");
		else{
			let tbs = inter.textbox_store[i];
			
			tbs.te = te;
			tbs.eqn = eqn;
			
			let origin = inter.equation.origin;
			if(origin.p != undefined && origin.cl != undefined) model.update_pline(origin.p,origin.cl);    
		
			inter.equation = {};
			inter.textbox_store.pop();
			
			bubble_check_error();
			
			press_bubble_OK();
		}
	}
}


/// Creates a new equation with a text string and a type. p and cl are available for transitions
function create_equation(te,type,p,cl)
{
	te = te.replace(/\*/g,"×");
	
	let i = find(eqn_types,"name",type);
	if(i == undefined){
		error("Equation type not found: "+type); 
		return;
	}		
	
	let eqn = {te:te, type:type, mode:eqn_types[i].mode};
	
	if(p != undefined) eqn.p_name = model.species[p].name;
	if(cl != undefined) eqn.cl_name = model.species[p].cla[cl].name;

	extract_equation_properties(eqn);
	
	return eqn;
}


/// Checks if a string is an equation with a single parameter or a number
function is_eqn(te,tag,op)
{
	let eqn = create_equation(te,"test");
	
	if(eqn.warn.length > 0) return err("For '"+tag+"' — "+eqn.warn[0].te); 
	
	if(eqn.mode == "param_only"){
		if(eqn.ind_eff.length > 0){
			return err("'"+tag+"' cannot contain individual effects"); 
		}
		
		if(eqn.param.length > 1){
			return err("'"+tag+"' can only be constant or contain one parameter"); 
		}
	}
	
	if(eqn.param.length == 0){
		if(op.positive == true){
			if(Number(te) <= 0){
				return err("'"+tag+"' must be a positive number"); 
			}
		}
		
		if(op.zero_positive == true){
			if(Number(te) < 0){
				return err("'"+tag+"' must be a non-negative number"); 
			}
		}
		
		if(op.zero_one_range == true){
			if(Number(te) < 0 || Number(te) > 1){
				return err("'"+tag+"' must be in the range 0 - 1"); 
			}
		}
	}
	
	return success();
}


/// Renames compartment in equations
function equation_rename_compartment(p,cl,old_name,new_name)
{
	let sp = model.species[p];
	let cl_name = sp.cla[cl].name;
		
	let eq_list = model.find_equation_list();
	
	let len = old_name.length;
	let dif = new_name.length - old_name.length;
	
	for(let i = 0; i < eq_list.length; i++){
		let eqn = eq_list[i];
		extract_equation_properties(eqn);
	
		let te = eqn.te;
		for(let j = 0; j < eqn.comp_name_list.length; j++){
			let ch = eqn.comp_name_list[j];
		
			if(ch.cl_name == cl_name && ch.comp_name == old_name){
				te = te.substr(0,ch.icur)+new_name+te.substr(ch.icur+len);
				for(let jj = j+1; jj < eqn.comp_name_list.length; jj++){
					eqn.comp_name_list[jj].icur += dif;
				}
			}
		}
		eqn.te = te;
	}	
}

"use strict";

/// Creates graphical interface for initial population
function create_initial_pop(lay)
{
	let so = inter.edit_source;

	lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, type:"Rect", val:WHITE});
	
	let cx = corner.x;
	let cy = corner.y;
	
	let ac = "AddDataSource";				
	let te = "Initial population";
		
	if(so.info.i != undefined){
		ac = "PlaceSourceBack"; te = "Edit initial population";
	}
			
	cy = lay.add_title(te,cx,cy,{te:init_pop_text2});

	add_layer("InitialPopulationContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-3.5,{});

	lay.add_corner_button([["Cancel","Grey","CancelEditSource"],["Done","Grey",ac]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
}


/// Creates graphical interface for initial population prior
function create_initial_pop_prior(lay)
{
	let so = inter.edit_source;
	
	lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, type:"Rect", val:WHITE});
	
	let cx = corner.x;
	let cy = corner.y;
	
	let ac = "AddDataSource";				
	let te = "Prior on initial population";

	if(so.info.i != undefined){
		ac = "PlaceSourceBack"; te = "Edit prior on initial population";
	}
	
	cy = lay.add_title(te,cx,cy,{te:initpopprior_text});

	add_layer("InitPopPriorContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-3.5,{});

	lay.add_corner_button([["Cancel","Grey","CancelEditSource"],["Done","Grey",ac]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
}


/// Generate the content for the add initial population page (graphical version) 
function add_initial_population_content(lay)
{
	let p = model.get_p();
	let sp = model.species[p];
	
	let mar = 0.5, mar2 = 1;
	let cx = mar;
	let cy = 0.2;

	let dy = 1.3;
	let si = dy*si_comp_text_frac;
	let fo = get_font(si,"","times");
	
	switch(inter.edit_source.spec.radio2.value){
	case "Focal":
		{
			if(cx != mar) cy += 1.5;
			cx = mar;
			
			cy = lay.add_paragraph("Set the populations in each of the compartments:",lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
		
			let te = inter.edit_source.spec.focal.te;
			let cl = find(sp.cla,"name",te);
			
			let claa = sp.cla[cl];
			
			lay.add_button({te:claa.name, x:mar, y:cy+0.06, dy:dy, type:"SubTitle", col:BLACK}); 
			cy += 2;
			
			let source = inter.edit_source.cla[cl];
			
			for(let c = 0; c < claa.ncomp; c++){
				let co = claa.comp[c];
				let te = co.name+"=";	
				let ac = "EditCompPop";
				let val = source.comp_init_pop[c].pop;
				te += val;
				
				let w = text_width(te,fo);
				let dx = w+2;
				
				if(cx+dx > lay.dx-mar2){ cx = mar; cy += dy+1;}	
				
				lay.add_button({te:te, x:cx, y:cy, dx:dx, dy:dy, p:p, cl:cl, c:c, type:"CompPop", val:val, ac:ac, col:co.col}); 
				cx += dx+1;
			}
			cx = mar; cy += dy+1;
			
			if(sp.ncla > 1){
				cy += 0.5;
				cy = lay.add_paragraph("Select the percentage of the population in different compartments:",lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
			
				for(let cl = 0; cl < sp.ncla; cl++){
					let claa = sp.cla[cl];
					if(inter.edit_source.spec.focal.te != claa.name){
						lay.add_button({te:claa.name, x:mar, y:cy+0.06, dy:dy, type:"SubTitle", col:BLACK}); 
						cy += 2;
						
						let source = inter.edit_source.cla[cl];
						
						for(let c = 0; c < claa.ncomp; c++){
							let co = claa.comp[c];
							let te = co.name+"=";
							
							let ac = "EditCompPopPercent"; 
							if(c == claa.ncomp-1) ac = "NotEditCompPopPercent";

							let val = source.comp_init_pop[c].pop_per;
							te += precision(val,4)+"%";
							
							let w = text_width(te,fo);
							let dx = w+2;
							
							if(cx+dx > lay.dx-mar2){ cx = mar; cy += dy+1;}	
							
							lay.add_button({te:te, x:cx, y:cy, dx:dx, dy:dy, p:p, cl:cl, c:c, type:"CompPop", val:val, ac:ac, col:co.col}); 
							cx += dx+1;
						}
						
						cx = mar; cy += dy+1;
					}
				}
			}
		}
		break;
		
	case "All":
		{
			let te = "Specify the population in each of the compartment combinations:";
			if(sp.ncla == 1) te = "Specify the population in each of the compartments:"
			
			cy = lay.add_paragraph(te,lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
			cy += 0.5;
			
			let glob = inter.edit_source.glob_comp;
			
			for(let c = 0; c < glob.length; c++){
				let claa = glob[c].cla;
				
				let col_list=[];
				
				let te = "";
				for(let cl = 0; cl < claa.length; cl++){
					if(cl > 0) te += ",";
					te += claa[cl];
					
					let c = find(sp.cla[cl].comp,"name",claa[cl]);
					col_list.push(sp.cla[cl].comp[c].col);
				}					
				te += "="+glob[c].pop;
							
				let w = text_width(te,fo);
				let dx = w+2;
							
				if(cx+dx > lay.dx-mar2){ cx = mar; cy += dy+1;}	
							
				lay.add_button({te:te, x:cx, y:cy, dx:dx, dy:dy, c:c, type:"MultiComp", ac:"EditGlobCompPop", col_list:col_list}); 
				cx += dx+1;
			}
		
			cy += 1.3;
		}
		break;	
	
	default: error("Option not recognised 104"); break;
	}
	
	lay.add_space(cy,0.3);
}


/// Allows the user to edit the prior for the initial population
function add_init_pop_prior_content(lay)
{
	let mar = 0.5;
	let cx = mar;
	let cy = 0.2;
	cy = lay.add_paragraph("Select which classification specifies population:",lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
}


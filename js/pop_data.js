"use strict";
// Functions relating to population data

/// Creates graphical interface for initial population
function create_initial_pop(lay)
{
	let so = edit_source;

	lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, type:"Rect", val:WHITE});
	
	let cx = corner.x;
	let cy = corner.y;
	
	let ac = "AddDataSource";				
	let title = "Initial population";
	let help = init_pop_text2;
	
	if(so.spec.radio_dist.value == "Dist"){
		title = "Initial population distribution";
		help = initpopdist_text;
	}
	
	let bute = "Done";
	if(so.info.i != undefined){
		ac = "PlaceSourceBack";
		title = "Edit initial population";
		if(so.spec.radio_dist.value == "Dist"){
			title = "Edit initial population distribution";	
		}
		bute = "Update";
	}
			
	cy = lay.add_title(title,cx,cy,{te:help});

	add_layer("InitialPopulationContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-3.5,{});

	lay.add_corner_button([["Cancel","Grey","CancelEditSource"],[bute,"Grey",ac]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
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
	let si_sup = si*0.6;

	switch(edit_source.spec.radio_dist.value){
	case "Fixed":
		switch(edit_source.spec.radio2.value){
		case "Focal":
			{
				if(cx != mar) cy += 1.5;
				cx = mar;
				
				cy = lay.add_paragraph("Set the populations in each of the compartments:",lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
			
				let te = edit_source.spec.focal.te;
				let cl = find(sp.cla,"name",te);
				
				let claa = sp.cla[cl];
				
				lay.add_button({te:claa.name, x:mar, y:cy+0.06, dy:dy, type:"SubTitle", col:BLACK}); 
				cy += 2;
				
				let source = edit_source.cla[cl];
				
				for(let c = 0; c < claa.ncomp; c++){
					let co = claa.comp[c];
					let te = co.name+"=";	
					let ac = "EditCompPop";
					let val = source.comp_init_pop[c].pop;
					te += val;
					
					let w = text_width(te,fo);
					let dx = w+2;
					
					if(cx+dx > lay.dx-mar2){ cx = mar; cy += dy+1;}	
					
					lay.add_button({te:te, x:cx, y:cy, dx:dx, dy:dy, p:p, cl:cl, c:c, type:"CompPop", val:val, ac:ac, col:co.col, col_dark:dark_colour(co.col)}); 
					cx += dx+1;
				}
				cx = mar; cy += dy+1;
				
				if(sp.ncla > 1){
					cy += 0.5;
					cy = lay.add_paragraph("Select the percentage of the population in different compartments:",lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
				
					for(let cl = 0; cl < sp.ncla; cl++){
						let claa = sp.cla[cl];
						if(edit_source.spec.focal.te != claa.name){
							lay.add_button({te:claa.name, x:mar, y:cy+0.06, dy:dy, type:"SubTitle", col:BLACK}); 
							cy += 2;
							
							let source = edit_source.cla[cl];
							
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
								
								lay.add_button({te:te, x:cx, y:cy, dx:dx, dy:dy, p:p, cl:cl, c:c, type:"CompPop", val:val, ac:ac, col:co.col, col_dark:dark_colour(co.col)}); 
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
				
				let glob = edit_source.glob_comp;
				
				for(let c = 0; c < glob.length; c++){
					let claa = glob[c].cla;
					
					let col_list=[];
					
					let te = "";
					for(let cl = 0; cl < claa.length; cl++){
						if(cl > 0) te += ",";
						te += claa[cl];
						
						let c = hash_find(sp.cla[cl].hash_comp,claa[cl]);
						col_list.push(sp.cla[cl].comp[c].col);
					}					
					te += "="+glob[c].pop;
					let tsa = text_sup_anno(te,si);
					let w = tsa.w;
					
					let dx = w+2;
								
					if(cx+dx > lay.dx-mar2){ cx = mar; cy += dy+1;}	
								
					lay.add_button({tsa:tsa, x:cx, y:cy, dx:dx, dy:dy, c:c, type:"MultiComp", ac:"EditGlobCompPop", col_list:col_list}); 
					cx += dx+1;
				}
			
				cy += 1.3;
			}
			break;	
		
		default: error("Option not recognised 104"); break;
		}
		break;
		
	case "Dist":
		switch(edit_source.spec.radio2.value){
		case "Focal":
			{
				if(cx != mar) cy += 1.5;
				cx = mar;
				
				cy = lay.add_paragraph("Set the population distribution in each of the compartments:",lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
			
				let te = edit_source.spec.focal.te;
				let cl = find(sp.cla,"name",te);
				
				let claa = sp.cla[cl];
				
				lay.add_button({te:claa.name, x:mar, y:cy+0.06, dy:dy, type:"SubTitle", col:BLACK}); 
				cy += 2;
				
				let source = edit_source.cla[cl];
				
				for(let c = 0; c < claa.ncomp; c++){
					let co = claa.comp[c];
					let val = get_prior_string(source.comp_init_dist[c].dist);
					
					let te = "N", sup = co.name, te2 = " ~ "+val;
					
					let sup_tsa = text_sup_anno(sup,si_sup,LARGE);
						
					let w = text_width(te,fo)+sup_tsa.w+text_width(te2,fo);
									
					let dx = w+2;
					
					if(cx+dx > lay.dx-mar2){ cx = mar; cy += dy+1;}	
					
					lay.add_button({te:te, sup_tsa:sup_tsa, te2:te2, x:cx, y:cy, dx:dx, dy:dy, p:p, cl:cl, c:c, prior:source.comp_init_dist[c].dist, type:"CompPrior", val:val, ac:"EditPopPrior", fo:fo, col:co.col}); 
					cx += dx+1;
				}
				cx = mar; cy += dy+1;
				
				if(sp.ncla > 1){
					cy += 0.5;
					cy = lay.add_paragraph("Select the α value for the Dirichlet distribution:",lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
				
					for(let cl = 0; cl < sp.ncla; cl++){
						let claa = sp.cla[cl];
						if(edit_source.spec.focal.te != claa.name){
							lay.add_button({te:claa.name, x:mar, y:cy+0.06, dy:dy, type:"SubTitle", col:BLACK}); 
							cy += 2;
							
							let source = edit_source.cla[cl];
							
							for(let c = 0; c < claa.ncomp; c++){
								let co = claa.comp[c];
								let val = source.comp_init_dist[c].alpha;
								
								let te = "α", sup = co.name, te2 = " = "+val;
								
								let sup_tsa = text_sup_anno(sup,si_sup,LARGE);
									
								let w = text_width(te,fo)+sup_tsa.w+text_width(te2,fo);
								
								let dx = w+2;
								
								if(cx+dx > lay.dx-mar2){ cx = mar; cy += dy+1;}	
								
								lay.add_button({te:te, sup_tsa:sup_tsa, te2:te2, x:cx, y:cy, dx:dx, dy:dy, source:source.comp_init_dist[c], type:"CompAlpha", val:val, ac:"EditAlpha", fo:fo, col:co.col}); 
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
				cy = lay.add_paragraph("Set the total population prior:",lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
				cy += 0.2;
				
				let val = get_prior_string(edit_source.pop_dist);
				let te = "N", sup = "total", te2 = " ~ "+val;
				let sup_tsa = text_sup_anno(sup,si_sup,LARGE);
				let w = text_width(te,fo)+sup_tsa.w+text_width(te2,fo);
				let dx = w+2;
			
				lay.add_button({te:te, sup_tsa:sup_tsa, te2:te2, x:cx, y:cy, dx:dx, dy:dy, prior:edit_source.pop_dist, type:"CompPrior", ac:"EditPopPrior", fo:fo, col:WHITE}); 
			
				cx = mar; cy += dy+1;
			
				cy = lay.add_paragraph("Select the α values for the Dirichlet prior in each of the compartment combinations:",lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
				cy += 0.5;
				
				let glob = edit_source.glob_comp;
				
				for(let c = 0; c < glob.length; c++){
					let claa = glob[c].cla;
					
					let col_list=[];
					
					let te = "α";
					let sup = "";
					for(let cl = 0; cl < claa.length; cl++){
						if(cl > 0) sup += ",";
						sup += claa[cl];
						
						let c = hash_find(sp.cla[cl].hash_comp,claa[cl]);
						col_list.push(sp.cla[cl].comp[c].col);
					}					
					let te2 = " = "+glob[c].alpha;
								
					let sup_tsa = text_sup_anno(sup,si_sup,LARGE);
								
					let w = text_width(te,fo)+sup_tsa.w+text_width(te2,fo);
					let dx = w+2;
								
					if(cx+dx > lay.dx-mar2){ cx = mar; cy += dy+1;}	
								
					lay.add_button({te:te, sup_tsa:sup_tsa, te2:te2, x:cx, y:cy, dx:dx, dy:dy, c:c, source:glob[c], type:"MultiCompAlpha", ac:"EditAlpha", col_list:col_list, fo:fo}); 
					cx += dx+1;
				}
			
				cy += 1.3;
			}
			break;	
		
		default: error("Option not recognised 104"); break;
		}
		break;
	}
	
	lay.add_space(cy,0.3);
}

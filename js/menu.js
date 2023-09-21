"use strict";

/// Initialises all the pages and scrollbars
function initialise_pages()                                     
{	
	let page=[];
	
	page.push({name:"Home", sub:[]});

	page.push({name:"Model", sub:[]});

	{
		let sub = page[page.length-1].sub;
		sub.push({name:"Compartments", sub:[]});
		
		{
			let subsub = sub[sub.length-1].sub;
			
			for(let p = 0; p < model.species.length; p++){
				subsub.push({name:model.species[p].name, sub:[]});
			
				let subsubsub = subsub[subsub.length-1].sub;
				for(let cl = 0; cl < model.species[p].ncla; cl++){	 
					subsubsub.push({ name:model.species[p].cla[cl].name});
				}
			}
		
			subsub.push({name:"Add species"});
		}
		
		sub.push({name:"Parameters", sub:[]});
		
		sub.push({name:"Description", sub:[]});
	}
	
	page.push({name:"Simulation", sub:[]});
	
	{
		let sub = page[page.length-1].sub;
	
		sub.push({name:"Population", sub:[]});
		
		{
			let subsub = sub[sub.length-1].sub;
		
			if(model.species.length > 1){
				for(let p = 0; p < model.species.length; p++){
					subsub.push({name:model.species[p].name, sub:[]});
				}
			}
		}
		
		sub.push({name:"Parameters", sub:[]});
		sub.push({name:"Run", sub:[]});
		
		if(sim_result.on == true){
			sub.push({name:"Results", sub:[]});
		
			let subsub = sub[sub.length-1].sub;
			subsub.push({name:"Populations"});
			subsub.push({name:"Transitions"});
			
			if(sim_result.species){
				if(find(sim_result.species,"type","Individual") != undefined){
					subsub.push({name:"Individuals"});
				}
			}
			
			if(sim_result.spline.length > 0) subsub.push({name:"Splines"});
			
			subsub.push({name:"Parameters"});
			
			sub.push({name:"Generate Data", sub:[]});
			
			{
				let subsub = sub[sub.length-1].sub;
			
				if(sim_result.species.length > 1){
					for(let p = 0; p < sim_result.species.length; p++){
						subsub.push({name:sim_result.species[p].name, sub:[]});
					}
				}
			}
		
		}
	}
	
	page.push({name:"Inference", sub:[]});
		
	{
		let sub = page[page.length-1].sub;
	
		sub.push({name:"Data", sub:[]});
		
		{
			let subsub = sub[sub.length-1].sub;
		
			if( model.species.length > 1){
				for(let p = 0; p < model.species.length; p++){
					subsub.push({name:model.species[p].name, sub:[]});
				}
			}
		}
		
		sub.push({name:"Prior", sub:[]});
		
		sub.push({name:"Run", sub:[]});
		
		if(inf_result.on == true){
			sub.push({name:"Results", sub:[]});
		
			let subsub = sub[sub.length-1].sub;
			subsub.push({name:"Parameters"});
			subsub.push({name:"Populations"});
			subsub.push({name:"Transitions"});
			
			if(inf_result.species){
				if(find(inf_result.species,"type","Individual") != undefined){
					subsub.push({name:"Individuals"});
				}
			}
			
			if(inf_result.spline){
				if(inf_result.spline.length > 0) subsub.push({name:"Splines"});
			}
		}
	}
	
	set_page_index(page);
	
	inter.page = page;
}


/// Sets the indexes for the pages
function set_page_index(page)
{
	// Copies page numbers from previous version
	let pold = inter.page;
	for(let p = 0; p < pold.length; p++){
		let pmap = 0; while(pmap < page.length && pold[p].name != page[pmap].name) pmap++;
		if(pmap < page.length){
			let max = page[pmap].sub.length; 
			if(max > 0){
				page[pmap].index = pold[p].index; 
				if(page[pmap].index >= max) page[pmap].index = max-1;
			}
			
			let pag_old = pold[p];
			let pag = page[pmap];
			for(let s = 0; s < pag_old.sub.length; s++){
				let smap = 0; while(smap < pag.sub.length && pag_old.sub[s].name != pag.sub[smap].name) smap++;
				if(smap < pag.sub.length){
					let max = pag.sub[smap].sub.length; 
					if(max > 0){
						pag.sub[smap].index = pag_old.sub[s].index;
						if(pag.sub[smap].index >= max) pag.sub[smap].index = max-1;
					}
					
					for(let ss = 0; ss < pag_old.sub[s].sub.length; ss++){
						let ssmap = 0; while(ssmap < pag.sub[smap].sub.length && pag_old.sub[s].sub[ss].name != pag.sub[smap].sub[ssmap].name) ssmap++;
						if(ssmap < pag.sub[smap].sub.length){
							if(pag.sub[smap].sub[ssmap].sub != undefined){
								let max = pag.sub[smap].sub[ssmap].sub.length; 
								if(max > 0){
									pag.sub[smap].sub[ssmap].index = pag_old.sub[s].sub[ss].index;
									if(pag.sub[smap].sub[ssmap].index >= max) pag.sub[smap].sub[ssmap].index = max-1;			
								}
							}
						}
					}
				}
			}
		}
	}
	
	// Sets all undefined page indexes to zero
	for(let p = 0; p < page.length; p++){
		let pag = page[p];
		if(pag.index == undefined) pag.index = 0;
		for(let s = 0; s < pag.sub.length; s++){
			if(pag.sub[s].index == undefined) pag.sub[s].index = 0;
			for(let ss = 0; ss < pag.sub[s].sub.length; ss++){
				if(pag.sub[s].sub[ss].sub != undefined){
					if(pag.sub[s].sub[ss].index == undefined) pag.sub[s].sub[ss].index = 0;
				}
			}
		}
	}
}


/// Sets all the page indexes to zero
function zero_page_index()
{
	let page = inter.page;
	for(let p = 0; p < page.length; p++){
		let pag = page[p];
		pag.index = 0;
		for(let s = 0; s < pag.sub.length; s++){
			pag.sub[s].index = 0;
			for(let ss = 0; ss < pag.sub[s].sub.length; ss++){
				pag.sub[s].sub[ss].index = 0;
			}
		}
	}
}


/// Adds the logo button
function add_logo_buts(lay) 
{
	let mar = 0.05*lay.dx;
	let ddx = lay.dx-2*mar;
	lay.add_button({x:mar, y:0, dx:ddx, dy:(370/927)*ddx, type:"Logo"});
}


/// Adds buttons assocaited with the main menu on the left
function add_menu_buts(lay)
{
	let dx = lay.dx;
	
	let ddy = nearest_pixel(1.7);
	let dby = nearest_pixel(0.15);
	let dby2 = nearest_pixel(1.5);
	let dysub = nearest_pixel(1.6);
	let dysubsub = nearest_pixel(1.3);
	let xsub = nearest_pixel(0.5);
	let xsubsub = nearest_pixel(1.0);
	
	let y = 0;

	for(let p = 0; p < inter.page.length; p++){
		let pag = inter.page[p];

		let ac;
		if(p == 0 || model.start) ac = "Tab";
	
		lay.add_button({te:pag.name, x:0, y:y, dx:dx, dy:ddy, ac:ac, type:"Tab", val:p}); 
		y += ddy; 
		y += dby;
		if(p == inter.pa && pag.sub.length > 0){	
			y += 0.5*dby;		
			for(let s = 0; s < pag.sub.length; s++){
				let sub = pag.sub[s];
			
				let xx = 0.7;
				
				lay.add_button({te:sub.name, x:xsub, y:y, dx:dx-xsub, dy:dysub, ac:"PageSub", type:"PageSub", val:s});
				y += dysub;
				
				if(s == inter.page[inter.pa].index && sub.sub.length > 0){
					y += dby;
					for(let ss = 0; ss < sub.sub.length; ss++){
						let subsub = sub.sub[ss];
						
						if(subsub.name == "Add species"){
							lay.add_button({te:"Add species", x:xx, y:y, dx:dx-xx, dy:dysubsub, ac:"AddSpecies", type:"AddSpecies", val:s});
						}
						else{
							lay.add_button({te:subsub.name, x:xsubsub, y:y, dx:dx-xsubsub, dy:dysubsub, ac:"PageSubSub", type:"PageSubSub", val:ss});
						}
						y += dysubsub;
					}
					y += nearest_pixel(0.5);
				}
				y += dby;
			}
			y += nearest_pixel(1);
		}
		y += dby2;
	}
}


/// Used to add a sub-menu
function add_submenu(submenu,dx,y,shift)                                     
{
	const xx = 15, dddy = dysub;
		
	for(let i = 0; i < submenu.length; i++){
		let ii = i; if(shift) ii += shift;
		//this.add_button({x:xx, y:y, dx:dx-xx, dy:inter.menu_submenu_size, type:"PageSubBack"});
		this.add_button({te:submenu[i], x:xx, y:y, dx:dx-xx, dy:inter.menu_submenu_size, ac:"PageSub", type:"PageSub", val:ii});
		y += inter.menu_submenu_size;
	}
	
	return y;
}


/// Used to add a sub-sub-menu
function add_subsubmenu(submenu,y,dx)                                     
{
	const xx = 28; dddy = dysubsub;

	for(let i = 0; i < submenu.length; i++){
		this.add_button({te:out[i], x:xx, y:y, dx:dx-xx, dy:dddy, ac:"PageSubSub", type:"PageSubSub", val:i}); y += dddy;
	}
	y += 5;
	return y;
}


/// Adds buttons with block any underlying content around the main paged (i.e. generate the rounded corners on the main area)
function add_frame_buts(lay) 
{
	let op = lay.op;
	/*
	lay.add_button({x:0, y:0, dx:op.x, dy:lay.dy, type:"Rect", col:BACKGROUND});
	lay.add_button({x:op.x, y:0, dx:op.dx, dy:op.y, type:"Rect", col:BACKGROUND});
	lay.add_button({x:op.x+op.dx, y:0, dx:lay.dx-(op.x+op.dx), dy:lay.dy, type:"Rect", col:BACKGROUND});
	lay.add_button({x:op.x, y:op.y+op.dy, dx:op.dx, dy:lay.dy-(op.y+op.dy), type:"Rect", col:BACKGROUND});
	*/
	lay.add_button({x:op.x, y:op.y, dx:op.dx, dy:op.dy, type:"Corners"});
	
	//lay.add_button({te:"File", x:op.x+op.dx-5, y:op.y, dx:3, dy:1.3, ac:"File", type:"File"});		
	
	//lay.add_button({x:op.x+op.dx-3, y:op.y+0.6, dx:1.8, dy:1.8, ac:"Menu", type:"Menu"});	
	
	let x = op.x+op.dx-4, y = op.y+1;
	let tree = inter.page_name.split("->");
	//if(tree[0] == "Model" && tree[1] == "Compartments"){ x = op.x+op.dx-3; y = op.y+0.6;}
		
	let show_model = model.get_show_model();
	if(show_model == true){
		if(model.species.length == 0 || model.species[model.get_p()].ncla == 0) show_model = false;
	}
	
	let coord; if(show_model == true) coord = model.get_coord();
		
	lay.add_button({x:x, y:y, dx:1.8, dy:1.8, ac:"Menu", type:"Menu", show_model:show_model, coord:coord});		
}


/// Sets the main dropdown menu
function drop_menu(cont,lay)
{
	let bub = inter.bubble;
	switch(bub.mode){
	case "LoadComps":
		{
			cont.dx = 10;
			bubble_addtitle(cont,"Load compartments",{te:loadcomp_text});
			
			if(bub.bu.coord == "cartesian"){
				bubble_addcheckbox(cont,0,"Load position data",bub.pos_check);
			}
			cont.y += 0.3;

			bubble_addradio(cont,0,"SelectColour","Select colour",bub.radio);
			bubble_colour(cont,2);
			cont.y += 0.3;

			bubble_addradio(cont,0,"FileColour","Load colour from file",bub.radio);
			
			add_end_button(cont,"Back","MenuBack");		
			add_end_button(cont,"Next","LoadCompartmentNext");		
		}
		break;
	
	case "LoadTrans":
		{
			cont.dx = 12;
			bubble_addtitle(cont,"Load transitions",{te:loadtrans_text});
				
			bubble_addparagraph(cont,"Select distribution for the loaded transitions:",0,cont.dx);
			cont.y += 0.3;
			bubble_addtext(cont,"Distribution:",1);
			cont.y -= 1.7;

			let dist = get_dist_pos();
			let pos=[]; for(let i = 0; i < dist.length; i++) pos.push({te:dist[i]});

			bubble_adddropdown(cont,5.6,6.5,inter.bubble.trans_type,pos);
			cont.y += 0.4;
			
			add_end_button(cont,"Back","MenuBack");		
			add_end_button(cont,"Load","ImportTrans2");		
		}
		break;
	
	case "LoadSource":
		{
			cont.dx = 11;
			bubble_addtitle(cont,"Load sources",{te:loadsource_text});
				
			bubble_addcheckbox(cont,0,"Load position data",bub.pos_check);
			
			bubble_addparagraph(cont,"Select distribution for the loaded sources:",0,cont.dx);
			cont.y += 0.5;
			bubble_addtext(cont,"Distribution:",0);
			cont.y -= 1.7;

			let dist = exp_dist_pos;
			let pos=[]; for(let i = 0; i < dist.length; i++) pos.push({te:dist[i]});
			
			bubble_adddropdown(cont,4.6,6.5,inter.bubble.trans_type,pos);
			cont.y += 0.6;
			
			add_end_button(cont,"Back","MenuBack");		
			add_end_button(cont,"Load","ImportSource2");		
		}
		break;
		
	case "LoadSink":
		{
			cont.dx = 11;
			bubble_addtitle(cont,"Load sinks",{te:loadsink_text});
				
			bubble_addcheckbox(cont,0,"Load position data",bub.pos_check);
			
			bubble_addparagraph(cont,"Select distribution for the loaded sinks:",0,cont.dx);
			cont.y += 0.3;
			bubble_addtext(cont,"Distribution:",0);
			cont.y -= 1.7;

			let dist = get_dist_pos();
			let pos=[]; for(let i = 0; i < dist.length; i++) pos.push({te:dist[i]});

			bubble_adddropdown(cont,4.6,6.5,inter.bubble.trans_type,pos);
			cont.y += 0.6;
			
			add_end_button(cont,"Back","MenuBack");		
			add_end_button(cont,"Load","ImportSink2");		
		}
		break;

	default:
		{				
			let dx = 6.5, dy = 1.2;
			cont.dx = dx;
			
			let menu = [];
			{
				let sub = [];
				sub.push({te:"Load", ac:"Load"});
				let active; if(model.filename == "") active = false;
				sub.push({te:"Save", ac:"Save", active:active});
				sub.push({te:"Save As...", ac:"SaveAs"});
				
				let active2 = true; 
				if(inter.page[inter.pa].name != "Simulation" && inter.page[inter.pa].name != "Inference"){
					active2 = false;
				}
				
				sub.push({te:"View Code", ac:"ViewCode", active:active2});
				
				menu.push({te:"File",sub:sub, help:file_text});
			}
			{
				let sub = [];
				sub.push({te:"Script", ac:"ImportScript"});
				let active = cont.bu.show_model;
				sub.push({te:"Compartments", ac:"ImportComp", active:active});
				let active2 = active; if(active2 == true && cont.bu.coord != "latlng") active2 = false;
				sub.push({te:"Comp. Map", ac:"ImportCompMap", active:active2});			
							
				sub.push({te:"Transitions", ac:"ImportTrans", active:active});
				sub.push({te:"Sources", ac:"ImportSource", active:active});
				sub.push({te:"Sinks", ac:"ImportSink", active:active});
				menu.push({te:"Import",sub:sub, help:import_text});
			}
			{
				let sub = [];
				let active = false; 
				
				sub.push({te:"Script", ac:"ExportScript", active:true});
				
				if(inter.graph.type) active = true;
				if(model.get_show_model()) active = true;
				
				sub.push({te:"Print", ac:"PrintGraph", active:active});
				sub.push({te:"Image (.png)", ac:"ExportGraph", active:active});
				
				let ac = "ExportTable";
				active = false;
				if(find(inter.layer,"name","TableContent")) active = true;
				
				if(inter.graph.type == "Matrix"){
					active = true;
					ac = "ExportMatrixTable";
				}
				
				sub.push({te:"Table (.csv)", ac:ac, active:active});
				
				menu.push({te:"Export",sub:sub, help:export_text});
			}
			
			cont.y -= 0.4;
			for(let i = 0; i < menu.length; i++){
				let te = menu[i].te;
				lay.add_button({te:te, x:0, y:cont.y, dx:dx, dy:dy, type:"MiniMenu"});

				let w = text_width(te,get_font(0.8));
				cont.lay.add_help_button(w+0.4,cont.y+1.7,{title:te, te:menu[i].help, back_col:BUBBLE_COL});
					
				cont.y += 1.6;
				for(let j = 0; j < menu[i].sub.length; j++){
					let sub = menu[i].sub[j];
					let ac; if(sub.active != false) ac = sub.ac;
					lay.add_button({te:sub.te, x:0, y:cont.y, dx:dx, dy:dy, ac:ac, type:"MenuBut"});
					
					cont.y += 1.2;
				}
			}
		}
		break;
	}
}

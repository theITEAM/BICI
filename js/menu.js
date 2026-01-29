"use strict";
// Functions which generate the menu

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
	
		sub.push({name:"Initial Conditions", sub:[]});
		add_species_options(sub,model.species);
		
		sub.push({name:"Parameters", sub:[]});
		sub.push({name:"Run", sub:[]});
		
		if(model.sim_res.on == true){
			sub.push({name:"Results", sub:[]});
		
			let subsub = sub[sub.length-1].sub;
			subsub.push({name:"Populations", sub:[]});
			add_species_options(subsub,model.sim_res.species_name);
			
			subsub.push({name:"Transitions", sub:[]});
			add_species_options(subsub,model.sim_res.species_name);
			
			if(model.sim_res.ind_flag == true){
				subsub.push({name:"Individuals", sub:[]});
				add_species_options(subsub,model.sim_res.species_name);
			}
			
			subsub.push({name:"Parameters"});
			
			sub.push({name:"Generate Data", sub:[]});
			add_species_options(sub,model.sim_res.species_name);
		}
	}
	
	page.push({name:"Inference", sub:[]});
		
	{
		let sub = page[page.length-1].sub;
	
		sub.push({name:"Data", sub:[]});
		add_species_options(sub,model.species);
		
		sub.push({name:"Prior", sub:[]});
		
		sub.push({name:"Run", sub:[]});
		
		if(model.inf_res.on == true){
			sub.push({name:"Results", sub:[]});
		
			let subsub = sub[sub.length-1].sub;
			subsub.push({name:"Parameters"});
	
			subsub.push({name:"Populations", sub:[]});
			add_species_options(subsub,model.inf_res.species_name);
		
			subsub.push({name:"Transitions", sub:[]});
			add_species_options(subsub,model.inf_res.species_name);
			
			if(model.inf_res.ind_flag == true){
				subsub.push({name:"Individuals", sub:[]});
				add_species_options(subsub,model.inf_res.species_name);
			}
			
			let alg = model.inf_res.plot_filter.details.algorithm.value;
			if(alg == "PAS-MCMC" || alg == "ABC-SMC") subsub.push({name:"Generations"});
	
			if(model.inf_res.plot_filter.diagnostics_on){
				subsub.push({name:"Diagnostics", sub:[]});
				add_species_options(subsub,model.inf_res.species_name);
			}	
			
			if(alg == "DA-MCMC" || alg == "PAS-MCMC"){
				subsub.push({name:"Extend", sub:[]});
			}
		}
	}
	
	if(model.inf_res.on == true){
		page.push({name:"Post. Simulation", sub:[]});
		
		let sub = page[page.length-1].sub;
		sub.push({name:"Population Mod.", sub:[]});
		add_species_options(sub,model.species);
		
		sub.push({name:"Parameter Mult.", sub:[]});
		
		sub.push({name:"Run", sub:[]});

		if(model.ppc_res.on == true){
			sub.push({name:"Results", sub:[]});
		
			let subsub = sub[sub.length-1].sub;
			subsub.push({name:"Populations", sub:[]});
			add_species_options(subsub,model.ppc_res.species_name);
			
			subsub.push({name:"Transitions", sub:[]});
			add_species_options(subsub,model.ppc_res.species_name);
			
			
			if(model.ppc_res.ind_flag == true){
				subsub.push({name:"Individuals", sub:[]});
				add_species_options(subsub,model.inf_res.species_name);
			}
			
			subsub.push({name:"Parameters"});
		}
	}
	
	set_page_index(page);
	
	inter.page = page;
}


/// Adds species options onto the menu
function add_species_options(sub,list)
{
	let subsub = sub[sub.length-1].sub;
		
	if(list.length > 1){
		for(let p = 0; p < list.length; p++){
			let te = list[p]; if(list[p].name) te = list[p].name;
			subsub.push({name:te, sub:[]});
		}
	}
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
function zero_page_index(not_zero_page)
{
	let page = inter.page;
	if(not_zero_page == true){ // Just sets defined and makes sure in range
		for(let p = 0; p < page.length; p++){
			let pag = page[p];
			if(pag.index == undefined || pag.index >= pag.sub.length) pag.index = 0;
			for(let s = 0; s < pag.sub.length; s++){
				let pag_sub = pag.sub[s];
				if(pag_sub.index == undefined || pag_sub.index >= pag_sub.sub.length) pag_sub.index = 0;
				for(let ss = 0; ss < pag_sub.sub.length; ss++){
					let pag_sub_sub = pag.sub[s].sub[ss];
					if(pag_sub_sub.index == undefined || pag_sub_sub.index >= pag_sub_sub.sub.length){ 		
						pag_sub_sub.index = 0;
					}
				}
			}
		}
	}
	else{
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
}


/// Adds the logo button
function add_logo_buts(lay) 
{
	let mar = 0.05*lay.dx;
	let ddx = lay.dx-2*mar;
	
	let ac; if(debug == true && testing == true) ac = "LoadInference"; 
	
	lay.add_button({x:mar, y:0, dx:ddx, dy:(370/927)*ddx, ac:ac, type:"Logo"});
}


/// Adds buttons assocaited with the main menu on the left
function add_menu_buts(lay)
{
	let dx = lay.dx;
	
	let ddy = nearest_pixel(1.7);
	let dby = nearest_pixel(0.15);
	let dby2 = nearest_pixel(1.1);
	let dysub = nearest_pixel(1.6);
	let dysubsub = nearest_pixel(1.3);
	let xsub = nearest_pixel(0.5);
	let xsubsub = nearest_pixel(1);
	let xsubsubsub = nearest_pixel(1);
	
	let tree = inter.page_name.split("->");
	if(tab_name() == "Simulation" && subtab_name() == "Reults" &&
  	model.sim_res.on && model.inf_res.on){
		dysubsub = nearest_pixel(1);
	}
		
	let y = 0;

	for(let p = 0; p < inter.page.length; p++){
		let pag = inter.page[p];

		let ac;
		if((p == 0 || model.start) && p != inter.pa) ac = "Tab";
	
		lay.add_button({te:pag.name, x:0, y:y, dx:dx, dy:ddy, ac:ac, type:"Tab", val:p}); 
		y += ddy; 
		y += dby;
		if(p == inter.pa && pag.sub.length > 0){	
			y += 0.5*dby;		
			for(let s = 0; s < pag.sub.length; s++){
				let sub = pag.sub[s];
			
				let xx = 0.7;
			
				let acsub; if(s != pag.index) acsub = "PageSub";
				lay.add_button({te:sub.name, x:xsub, y:y, dx:dx-xsub, dy:dysub, ac:acsub, type:"PageSub", val:s});
				y += dysub;
				
				if(s == inter.page[inter.pa].index && sub.sub.length > 0){
					y += dby;
					for(let ss = 0; ss < sub.sub.length; ss++){
						let subsub = sub.sub[ss];
						
						if(subsub.name == "Add species"){
							lay.add_button({te:"Add species", x:xx, y:y, dx:dx-xx, dy:dysubsub, ac:"AddSpecies", type:"AddSpecies", val:s});
						}
						else{
							let acsubsub = "PageSubSub";
				
							lay.add_button({p:ss, te:subsub.name, x:xsubsub, y:y, dx:dx-xsubsub, dy:dysubsub, ac:acsubsub, type:"PageSubSub", val:ss});
						}
						
						y += dysubsub;
						
						if(sub.name != "Compartments" && subsub.sub && ss == sub.index && subsub.sub.length > 1){
							for(let sss = 0; sss < subsub.sub.length; sss++){
								let subsubsub = subsub.sub[sss];
								lay.add_button({p:sss, te:subsubsub.name, x:xsubsubsub, y:y, dx:dx-xsubsubsub, dy:dysubsub, ac:"PageSubSubSub", type:"PageSubSubSub", val:sss});
								y += dysubsub;
							}
						}
					}
						
					y += nearest_pixel(0.5);
				}
				y += dby;
			}
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
	
	lay.add_button({x:op.x, y:op.y, dx:op.dx, dy:op.dy, type:"Corners"});
	
	let x = op.x+op.dx-4, y = op.y+1;
		
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
	
	let graph_show = false;
	if(find(inter.layer,"name","GraphContent") != undefined) graph_show = true; 
	
	if(bub.mode == "SelectLine"){
		let dx = 8, dy = 1.2;
		cont.dx = dx;
		bubble_addtitle(cont,"Export data",{title:"Select line", te:exportline_text});
	
		let dg = inter.graph.data_group;
	
		cont.y += 0.2;
		
		for(let i = 0; i < dg.length; i++){
			bubble_addradio(cont,0,i,dg[i].name,inter.bubble.line_radio); 
		}			
		cont.y += 0.2;
	
		add_end_button(cont,"Export","ExportLine");
	
		return;
	}
		
	if(bub.mode == "Video"){
		let dx = 8, dy = 1.2;
		cont.dx = dx;
		
		bubble_addtitle(cont,"Export Video",{title:"Export video", te:exportvideo_text});
	
		bubble_add_minititle(cont,"Quality:");
		cont.y -= 0.4;
		bubble_addradio(cont,0,"Low","Low quality",inter.mp4quality_radio);
		bubble_addradio(cont,0,"Medium","Medium quality",inter.mp4quality_radio); 
		bubble_addradio(cont,0,"High","High quality",inter.mp4quality_radio);
		cont.y += 0.3;
		bubble_input(cont,"Frames per second:",{type:"fps"});
		add_end_button(cont,"Export","ExportVideo2");
		return;
	}
		
	let dx = 6.5, dy = 1.2;
	cont.dx = dx;
	
	let active_mod = true; if(model.start == false) active_mod = false;	
		
	let menu = [];
	{
		let sub = [];
		sub.push({te:"Load", ac:"Load"});
		let active; if(inter.file_store.filename == "") active = false;
		sub.push({te:"Save", ac:"Save", active:active});
		sub.push({te:"Save As...", ac:"SaveAs", active:active_mod});
		
		{
			let ac; if(model.sim_res.on) ac = "SaveSim";
			sub.push({te:"Save Simulation", ac:ac});
		}
		
		{
			let ac; if(model.inf_res.on) ac = "SaveInf";
			sub.push({te:"Save Inference", ac:ac});
		}
		
		{
			let ac; if(model.ppc_res.on) ac = "SavePPC";
			sub.push({te:"Save Post. Sim.", ac:ac});
		}

		{
			sub.push({te:"View Code", ac:"ViewCode", active:active_mod});
		}
		
		menu.push({te:"File",sub:sub, help:file_text});
	}
	
	{
		let sub = [];
		
		{
			sub.push({te:"Script", ac:"ExportScript", active:active_mod});
		}
		
		let active = false; 
		if(graph_show) active = true;
		if(model.get_show_model()) active = true;
		
		sub.push({te:"Print", ac:"PrintGraph", active:active});
		sub.push({te:"Image (.png)", ac:"ExportGraph", active:active});
		
		let active_video = false;	
		if(find(inter.layer,"name","AnimControls") != undefined) active_video = true;
	
		sub.push({te:"Video (.mp4)", ac:"ExportVideo", active:active_video});
		
		let ac = "ExportTable";

		active = false;
		for(let l = 0; l < inter.layer.length; l++){
			switch(inter.layer[l].name){
			case "TableContent": 
			case "CreateEditXvectorContent":
			case "CreateEditTableContent":
			case "CreateEditParamContent":
			case "CreateEditAmatrixContent":
				active = true;
				ac = "ExportTableContent"
				break;
			}
		}
	
		if(graph_show){
			if(inter.graph.type == "Matrix"){
				active = true;
				ac = "ExportMatrixTable";
			}
					
			switch(inter.graph.data_group.length){
			case 0: break;
			case 1: active = true; ac = "ExportDataLine"; break;
			default: active = true; ac = "SelectDataLine"; break;
			}				
		}
		
		sub.push({te:"Table (.csv)", ac:ac, active:active});
		
		active = false;
		if(subtab_name() == "Results") active = true;
		sub.push({te:"Params (.csv)", ac:"OutputParams", active:active});
		sub.push({te:"States (.txt)", ac:"OutputStates", active:active});
		
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

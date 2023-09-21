"use strict";

class Model
{
	start;                         // Set to true when a new model is started 
	
	description;                   // A description of the analysis
	filename;                      // The current filename 
	species = [];                  // Stores information about species
	param = [];                    // Stores parameter information
	spline = [];                   // Stores spline information
	warn = [];                     // Stores model warnings
	derive = [];                   // Sets derived quantities
	
	sim_details = {};              // Details of simulation
	inf_details = {};              // Details of inference
	
	constructor()
	{
		this.start = false;
	}

	///	Starts a new model
	start_new(op)
	{
		this.description = {te:"Model and analysis description."};
		this.species = [];
		this.param = [];
		this.derive = [];   
		this.sim_details = { t_start:"", t_end:"", timestep:"", indmax:10000, algorithm:{value:"gillespie"}, number:1};    
		this.inf_details = { t_start:"", t_end:"", timestep:"", abcsample:"1000", sample:"10000", thinparam:"10", thinstate:"50", accfrac:"0.1", accfracsmc:"0.5", numgen:"5", kernelsize:"0.5", indmax:10000, nchain:"3", algorithm:{value:"DA-MCMC"}}; 
		
		this.start = true;
		this.filename = "";
		initialise_pages();
		zero_page_index();
		
		if(op != "without results"){
			sim_result = clear_result();
			inf_result = clear_result();
		}
		
		initialise_pages();
	}
	
	
	/// Gets the species number currently being viewed
	get_p()
	{
		let pag = inter.page[inter.pa];
		let sub = pag.sub[pag.index];
		if(!(pag.name =="Simulation" && sub.name == "Population") &&
		  !(pag.name =="Simulation" && sub.name == "Generate Data") &&
			!(pag.name =="Inference" && sub.name == "Data") &&
			!(pag.name =="Model" && sub.name == "Compartments")){
			error("Should not be getting p");
		}
		return sub.index;
	}
	
	
	/// Gets the species currently being viewed
	get_sp()
	{
		return this.species[this.get_p()];
	}
	
	
	/// Gets the classification being viewed
	get_cl()
	{
		let pag = inter.page[inter.pa];
		let p = pag.sub[pag.index].index;
		return pag.sub[pag.index].sub[p].index;
	}

	get_cla()
	{
		return this.species[this.get_p()].cla[this.get_cl()];
	}

	///	Returns true if the compartmental model is being viewed 
	get_show_model()
	{
		if(this.warn.length > 0) return false;
		
		let tree = inter.page_name.split("->");
		if(tree[0] == "Model" && tree[1] == "Compartments"){
			return true;
		}
		return false;
	}
	
	
	/// Adds all the compartment buttons
	add_compartment_buts(lay)
	{
		let p = this.get_p();
		let cl = this.get_cl();
		
		if(this.species.length == 0 || cl == this.species[p].ncla) return;

		let claa = this.species[p].cla[cl];
	
		let mo = inter.mode;
		if(mo.type == "Add_Transition" || mo.type == "Add_Source" || mo.type == "Add_Sink"){
			lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, ac:"IntermediateAdd", type:"Nothing", p:p, cl:cl});
		}
		
		let cam = claa.camera;                    // Loads up the camera used to view model
		let bub = inter.bubble;

		for(let k = 0; k < claa.ncomp; k++){
			let c = claa.comp[k];
		
			let ac = "Compartment";
			if((inter.mode.type == "Add_Source" || inter.mode.type == "Add_Sink") && inter.mode.midp.length == 0) ac = undefined;	

			if(bub.lay_name != undefined && bub.bu.type == "CompMap" && bub.bu.i == k){
				ac = undefined;
			}
			
			switch(c.type){
			case "box":
				{
					let pt = trans_point(c.x,c.y,cam,lay);
					let w =	c.w*cam.scale;
					let h =	c.h*cam.scale; 
					let x = pt.x-w/2, y = pt.y-h/2;
					lay.add_button({te:claa.comp[k].name, x:x, y:y, dx:w, dy:h, ac:ac, type:"Compartment", col:c.col, p:p, cl:cl, i:k});
				}
				break;
			
			case "latlng":
				{
					let pt = trans_point(c.x,c.y,cam,lay);
					let w =	2*latlng_radius*cam.scale*cam.ruler, h = w; 
					let x = pt.x-w/2, y = pt.y-h/2;
					lay.add_button({te:claa.comp[k].name, x:x, y:y, dx:w, dy:h, ac:ac, type:"CompLatLng", col:c.col, p:p, cl:cl, i:k});
				}
				break;
				
			case "boundary":
				{
					let box = c.feature.box;
				
					let pmin = trans_point(box.xmin,box.ymin,cam,lay);
					let pmax = trans_point(box.xmax,box.ymax,cam,lay);
					
					lay.add_button({x:pmin.x, y:pmin.y, dx:pmax.x-pmin.x, dy:pmax.y-pmin.y, ac:ac, type:"CompMap", polygon:c.feature.polygon, mask:c.mask, col:c.col, p:p, cl:cl, i:k});
				}
				break;

			default: error("Option not recognised 82"); break;
			}
		}
		
		if(mo.type == "Add_Compartment" || mo.type == "Add_Label"){
			lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, ac:"NothingMode", type:"Nothing"});	
		}
	}


	/// Adds all the buttons associated with transitions
	add_transition_buts(lay)
	{
		let pag = inter.page[inter.pa];
		let p = pag.sub[pag.index].index;
		let cl = pag.sub[pag.index].sub[p].index;
		
		if(this.species.length == 0 || cl == this.species[p].ncla) return;
				
		let claa = this.species[p].cla[cl];
		
		let mo = inter.mode;
		
		let trans_temp_add = false;
		let add_mouse = false;
		if(mo.type == "Add_Transition" && within_layer("Compartment") == true){
			if(mo.i != undefined){
				let f;
				let midp=[]; for(let k = 0; k < mo.midp.length; k++) midp.push(mo.midp[k]);
				
				let flag = false;
				
				let l = inter.over.layer;
				if(l != undefined){
					let bu = inter.layer[l].but[inter.over.i];
					if(bu != undefined){
						switch(bu.ac){
						case "Compartment": f = bu.i; break;
						case "IntermediateAdd": break;
						default: flag = true; break;
						}
					}
				}
				
				if(f == undefined){ midp.push(this.get_mouse_to_desktop()); add_mouse = true;}
				
				if(flag == false){
					this.add_transition(p,cl,mo.i,f,midp,"");
				
					trans_temp_add = true;
				}
			}
		}
		
		if(mo.type == "Add_Source" && within_layer("Compartment") == true){		
			let f;
			let midp=[]; for(let k = 0; k < mo.midp.length; k++) midp.push(mo.midp[k]);
			
			let flag = false;
			let l = inter.over.layer;
			if(l != undefined){
				let bu = inter.layer[l].but[inter.over.i];
				if(bu != undefined){
					switch(bu.ac){
					case "Compartment": if(inter.mode.midp.length > 0) f = bu.i; break;
					case "IntermediateAdd": break;
					default: flag = true; break;
					}
				}
			}
				
			if(f == undefined){ midp.push(this.get_mouse_to_desktop()); add_mouse = true;}
			
			if(flag == false){
				this.add_transition(p,cl,"Source",f,midp,"");
					
				trans_temp_add = true;
			}
		}
		
		if(mo.type == "Add_Sink" && within_layer("Compartment") == true){		
			let i;
			let midp=[]; for(let k = 0; k < mo.midp.length; k++)	midp.push(mo.midp[k]);
			
			let flag = false, flag_na = false;
			let l = inter.over.layer;
			if(l != undefined){
				let bu = inter.layer[l].but[inter.over.i];
				if(bu != undefined){
					switch(bu.ac){
					case "Compartment": if(inter.mode.midp.length > 0) i = bu.i; break;
					case "IntermediateAdd": break;
					case "TransitionPoint": flag_na = true; break;
					default: flag = true; break;
					}
				}
				else flag_na = true;
			}
					
			if(i == undefined && flag_na == false){
				midp.splice(0,0,this.get_mouse_to_desktop()); add_mouse = true;
			}
				
			if(flag == false){
				this.add_transition(p,cl,i,"Sink",midp,"");
						
				trans_temp_add = true;
			}
		}
		
		this.add_transition_buts2(p,cl,inter.mode.type,claa,add_mouse,lay);
	
		if(trans_temp_add == true){
			claa.tra.pop(); claa.ntra--;
		}
	}


	/// Adds the basic transitions
	add_transition_buts2(p,cl,mode_type,claa,add_mouse,lay)
	{
		let cam = claa.camera;

		for(let k = 0; k < claa.ntra; k++){
			let tr = claa.tra[k];
		
			let points=[];
			for(let i = 0; i < tr.pline.length; i++){
				let pl = tr.pline[i];
				points[i] = {x:lay.dx/2 + (pl.x-cam.x)*cam.scale, y:lay.dy/2 + (pl.y-cam.y)*cam.scale, index:pl.index};
			}
			
			let center;
			if(tr.center != undefined){
				let cen = tr.center;
				center = {x:lay.dx/2 + (cen.x-cam.x)*cam.scale, y:lay.dy/2 + (cen.y-cam.y)*cam.scale, nx:cen.nx, ny:cen.ny, w:cen.w};
			}
			
			let d = TRANS_OVER_RANGE;
		
			let ac = "Transition";
			if(mode_type == "Add_Source" || mode_type == "Add_Sink" || mode_type == "No Click") ac = undefined;
			
			lay.add_button({x:lay.dx/2 + (tr.box.x-cam.x)*cam.scale - d, y:lay.dy/2 + (tr.box.y-cam.y)*cam.scale-d, dx:tr.box.dx*cam.scale+2*d, dy:tr.box.dy*cam.scale+2*d, ac:ac, type:"Transition", cam:cam, points:points, center:center, tr:tr, p:p, cl:cl, i:k});	
				
			for(let j = 0; j < points.length; j++){
				let po = points[j];
				if(po.index >= 0 && po.index < tr.midp.length){
					let r = TRANS_POINT_R; 
					let variety = "Midpoint";
					let ac = "TransitionPoint";
					if(tr.i == "Source" && j == 0){ 
						r = 0.7*cam.scale*cam.ruler;
						variety = "Source"; ac = "Transition";
					}
					if(tr.f == "Sink" && j == points.length-1){ 
						r = 0.7*cam.scale*cam.ruler;
						variety = "Sink"; ac = "Transition";
					}

					if((variety == "Source" || variety == "Sink") && cam.coord == "latlng") r *= 0.6;
					
					if(mode_type == "Add_Transition" || mode_type == "Add_Source" || mode_type == "Add_Sink"){
						if(k < claa.ntra-1) ac = undefined;
						else{
							if(add_mouse == true){
								if(mode_type == "Add_Sink"){
									if(po.index == 0) ac = undefined;
								}
								else{
									if(po.index == tr.midp.length-1) ac = undefined;
								}
							}
							if(variety == "Source" || variety == "Sink") ac = undefined;
						}
					}

					if(ac != undefined || variety == "Source" || variety == "Sink"){
						lay.add_button({x:po.x - r, y:po.y - r, dx:2*r, dy:2*r, ac:ac, type:"TransitionPoint", tr:tr, p:p, cl:cl, i:k, index:po.index, variety:variety, center:center, tr:tr, p:p, cl:cl, i:k});
					}
				}
			}
		}
	}
/*
	/// Sets camera based on the location of compartments, transitions and annotations
	set_camera(p,cl)
	{
		let scw = page_char_wid-menu_width;
		let sch = page_char_hei;
		
		let xmin = LARGE, xmax = -LARGE, ymin = LARGE, ymax = -LARGE;
		
		let claa = this.species[p].cla[cl];
		for(let i = 0; i < claa.ncomp; i++){
			let co = claa.comp[i];
			let x = co.x; if(x < xmin) xmin = x; if(x > xmax) xmax = x;
			let y = co.y; if(y < ymin) ymin = y; if(y > ymax) ymax = y;
			
			x = co.x+co.w; if(x < xmin) xmin = x; if(x > xmax) xmax = x;
		  y = co.y+co.h; if(y < ymin) ymin = y; if(y > ymax) ymax = y;
		}
		
		for(let i = 0; i < claa.ntra; i++){
			let traa = claa.tra[i];
			for(let i = 0; i < traa.midp.length; i++){
				let x = traa.midp[i].x; if(x < xmin) xmin = x; if(x > xmax) xmax = x;
				let y = traa.midp[i].y; if(y < ymin) ymin = y; if(y > ymax) ymax = y;
			}
		}
	
		for(let i = 0; i < claa.annotation.length; i++){
			let an = claa.annotation[i];
			switch(an.type){
			case "text":
		
				let x = an.x; if(x < xmin) xmin = x; if(x > xmax) xmax = x;
				let y = an.y; if(y < ymin) ymin = y; if(y > ymax) ymax = y;
				break;
			}
		}
		
		let sca = scw/(xmax-xmin);
		if(0.9*sch/(ymax-ymin) < sca) sca = 0.9*sch/(ymax-ymin);
		sca *= 0.9;
		if(sca > 1) sca = 1;
		
		claa.camera = {x:(xmin+xmax)/2, y:(ymin+ymax)/2, scale:sca};
	}
	*/
	
	
	/// Adds a compartment button
	add_addcompartment_buts(lay)
	{
		inter.mode.l = lay.index;

		lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, type:"AddCompartmentMode"});
	}
	
	
	/// Adds a label button
	add_addlabel_buts(lay)
	{
		let mo = inter.mode;
		mo.l = lay.index;
		lay.add_button({te:mo.te, x:0, y:0, dx:lay.dx, dy:lay.dy, si:mo.si, font:mo.fo, type:"Text", col:GREY});
	}
	

	/// Adda a label annotation to the model
	add_label(te,p,cl)
	{
		let mo = inter.mode;
		let claa = this.species[p].cla[cl];
		claa.annotation.push({type:"text", te:te, tesize:si_annotation, x:mo.mp.x, y:mo.mp.y, color:BLACK});
	}

	
	/// Gets the mouse position on the desktop (i.e. the position in pixels is converted to model units accounting for camera)
	get_mouse_to_desktop()
	{
		let lay_main = get_lay("Main");
			
		let pag = inter.page[inter.pa];
		let p = pag.sub[pag.index].index;
		let cl = pag.sub[pag.index].sub[p].index;
		let cam = this.species[p].cla[cl].camera;

		let xc = ((inter.mx - lay_main.x - lay_main.dx/2)/cam.scale)+cam.x;
		let yc = ((inter.my - lay_main.y - lay_main.dy/2)/cam.scale)+cam.y;
			
		if(cam.coord == "cartesian"){
			return {x:Math.round(xc/grid)*grid, y:Math.round(yc/grid)*grid};
		}
	
		return {x:xc, y:yc};	
	//	return {x:Math.round(xc/grid)*grid, y:Math.round(yc/grid)*grid};
	}


	/// Gets the scale 
	get_scale()
	{
		return this.species[this.get_p()].cla[this.get_cl()].camera.scale;
	}
		
	
	/// Gets the coordinate system
	get_coord()
	{
		return this.species[this.get_p()].cla[this.get_cl()].camera.coord;
	}
	
	/*
	/// Gets the coord 
	get_coord()
	{
		let pag = inter.page[inter.pa];
		let p = pag.sub[pag.index].index;
		let cl = pag.sub[pag.index].sub[p].index;
		let cam = this.species[p].cla[cl].camera;
		return cam.coord;
	}
	*/
	
	/// Gets the ruler
	get_ruler()
	{
		return this.species[this.get_p()].cla[this.get_cl()].camera.ruler;
	}
		
		
	/// Converts from posiiton on desktop (i.e. model) to position on screen 
	get_desktop_to_image(mp)
	{
		let lay_main = get_lay("Main");
			
		let pag = inter.page[inter.pa];
		let p = pag.sub[pag.index].index;
		let cl = pag.sub[pag.index].sub[p].index;
		let cam = this.species[p].cla[cl].camera;

		let x = lay_main.dx/2 + (mp.x-cam.x)*cam.scale; 
		let y = lay_main.dy/2 + (mp.y-cam.y)*cam.scale;

		return {x:x, y:y};
	}


	/// In the "add compartment" mode this adds the new compartment
	place_add_compartment()
	{
		let lay_comp = inter.layer[inter.mode.l];
					
		lay_comp.x = 1000;
		
		let l = inter.over.layer;
		if(l != undefined){
			if(inter.layer[l].but[inter.over.i].ac == "NothingMode"){
				let mp = this.get_mouse_to_desktop();
				let ip = this.get_desktop_to_image(mp);
				
				let lay_main = get_lay("Main");
			
				let pag = inter.page[inter.pa];
				let p = pag.sub[pag.index].index;
				let cl = pag.sub[pag.index].sub[p].index;

				lay_comp.x = lay_main.x + ip.x-lay_comp.dx/2;
				lay_comp.y = lay_main.y + ip.y-lay_comp.dy/2; 
			
				let mo = inter.mode;
				mo.mp = mp;
				mo.p = p;
				mo.cl = cl;
				mo.coord = this.species[p].cla[cl].camera.coord;
			}
		}
		
		plot_screen();
	}


	/// Adds a new species to the model
	add_species(name,type)
	{
		let p = find(this.species,"name",name);
		if(p != undefined){
			alertp("There is already a species with the name '"+name+"'");
			return;
		}
		
		for(let p = 0; p < this.species.length; p++){
			let sp = this.species[p];
			
			let i = find(sp.cla,"name",name);
			if(i != undefined){
				alertp("Cannot have the same name as classification '"+sp.cla[i].name+"' ");
				return;
			}
			
			let j = find(sp.cla,"index",name);
			if(j != undefined){
				alertp("Cannot have the same name as index for classification '"+sp.cla[j].name+"' ");
				return;
			}
		
			for(let cl = 0; cl < sp.ncla; cl++){
				let c = find(sp.cla[cl].comp,"name",name);
				if(c != undefined){
					alertp("Cannot have the same name as compartment '"+sp.cla[cl].comp[i].name+"' in classification '"+sp.cla[cl].name+"'");
					return;
				}
			}
		}
		
		this.species.push({ name:name, sim_source:[], inf_source:[], gen_source:[], data_type:{te:data_types[0]
		}, cla:[], ncla:0, fix_eff:[], ind_eff_group:[], type:type});		
			
		return "success";
	}
	
	/// Adds a new classification to the model
	add_classification(p,name,index,op)
	{
		if(p == undefined){
			alertp("Cannot add a classification when the species is not set");
			return;
		}
		
		let cl = find(this.species[p].cla,"name",name);
		if(cl != undefined){
			alertp("There is already a classification with the name '"+name+"'");
			return;
		}
		
		if(name == index){
			alertp("Name and index must be different");
			return;
		}
		
		let coord = op.coord; if(coord == undefined) coord = "cartesian";
		
	
		let cam = {x:0, y:0, scale:1, ruler:1, coord:coord};
		
		let dmap = false; if(coord == "latlng" && op.default_map == true){ dmap = "loading"; cam.ruler = 0.05;}
		
		this.species[p].cla.push({ name:name, index:index, tra:[], ntra:0, comp:[], ncomp:0, annotation:[], default_map:dmap, camera:cam});
		this.species[p].ncla++;
		
		data_update_add_classification(p);
		check_data_valid();
		
		return "success";
	}
	
	
	/// Loads up default map (if neccesary);
	load_default_map()
	{
		for(let p = 0; p < this.species.length; p++){
			let sp = this.species[p];
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = this.species[p].cla[cl];
				if(claa.default_map == "loading"){
					claa.default_map = true;
					load_map("world",p,cl);
				}
			}
		}
	}
	

	/// Clones a classification from another species
	clone_class(p_to,p_from,cl_from)
	{
		let cla_from = this.species[p_from].cla[cl_from];
	
		this.add_classification(p_to,cla_from.name,cla_from.index,{coord:cla_from.camera.coord});
		
		let cl = this.species[p_to].ncla-1;
		
		this.species[p_to].cla[cl] = copy(cla_from);
		
		let cla_to = this.species[p_to].cla[cl];
		cla_to.tra = [];
		cla_to.ntra = 0;
	}


	/// Adds a new compartment to the model
	add_compartment(name,p,cl,x,y,col,fix,include_clone)
	{
		if(fix == undefined) fix = false;
		
		if(this.check_comp_exist(name,p) == true) return;
	
		let claa = this.species[p].cla[cl];
		let comp = claa.comp;
			
		comp.push({name:name, type:"box", x:x, y:y, col:col, fixed:{check:fix}, markov_branch:false});
		claa.ncomp++;
			
		let i = comp.length-1;
		this.set_compartment_size(comp[i]);
		
		this.snap_comp_to_grid(comp[i]);
		this.determine_branching();
		
		data_update_add_compartment(p,cl);
		
		if(include_clone == true) this.update_param_index(claa.index); // Updates param
		
		check_data_valid();

		if(include_clone == true){  // Adds to clones 
			let list = this.find_clones(p,cl);
		
			for(let li of list){
				this.add_compartment(name,li.p,li.cl,x,y,col,fix,false);
			}
		}
		
		return "success";
	}
	
		
	/// Updates parameter if the is a compartmental change in an index
	update_param_index(index)
	{
		for(let th = 0; th < this.param.length; th++){
			let par = this.param[th];
		
			let j = find_in(par.dep,index);
			if(j != undefined){
				let par_copy = copy(par);
				copy_param_info(par,par_copy);
			}
		}
	}
	
	
	/// Adds a new compartment to the model
	add_latlng_compartment(name,p,cl,x,y,col,fix)
	{
		if(fix == undefined) fix = false;
	
		if(this.check_comp_exist(name,p) == true) return;
		
		let claa = this.species[p].cla[cl];
		let comp = claa.comp;
		comp.push({name:name, type:"latlng", x:x, y:y, col:col, fixed:{check:fix}, markov_branch:false});
		claa.ncomp++;
		
		let i = comp.length-1;
		
		this.determine_branching();

		data_update_add_compartment(p,cl);
		check_data_valid();
		
		return "success";
	}
	
	
	/// Finds all the clones of a given classification
	find_clones(p,cl)
	{
		let list=[];
		let na =  this.species[p].cla[cl].name;
		for(let p2 = 0; p2 < this.species.length; p2++){
			if(p2 != p){
				let sp2 = this.species[p2];
				for(let cl2 = 0; cl2 < sp2.ncla; cl2++){
					let claa2 = sp2.cla[cl2];
					if(claa2.name == na){
						list.push({p:p2, cl:cl2});
					}
				}
			}
		}
		
		return list;
	}
	
	
	/// Gets a new unique compartment name
	get_new_comp_name()
	{
		let num = 0;
		do{
			let name = "";
			if(num > 0) name += "C"+num;
			
			let flag = false;
			
			for(let p = 0; p < this.species.length; p++){
				let sp = this.species[p];
				if(sp.name == name) flag = true;
				for(let cl = 0; cl < sp.ncla; cl++){
					let claa = sp.cla[cl];
					if(claa.name == name) flag = true;
					if(claa.index == name) flag = true;
					for(let c = 0; c < claa.ncomp; c++){
						let co = claa.comp[c];
						if(co.name == name) flag = true;
					}
				}
			}
			if(flag == false) return name;
			num++;
		}while(true);
	}
	
	
	/// Adds a new map compartment to the model
	add_map_compartment(name,p,cl,feature,mask,col)
	{
		if(this.check_comp_exist(name,p) == true) return;
		
		/// Calculates the midpoint of the region as the middle of the mask
		let claa = this.species[p].cla[cl];
		let comp = claa.comp;
		let xav = 0, yav = 0, nav = 0;
		for(let j = 0; j < mask_size; j++){
			for(let i = 0; i < mask_size; i++){
				if(mask[j*mask_size+i] == 1){
					xav += i; yav += j; nav++;
				}
			}
		}
		let box = feature.box;
		let xmid = box.xmin+((xav/nav)/mask_size)*(box.xmax-box.xmin);
		let ymid = box.ymin+((yav/nav)/mask_size)*(box.ymax-box.ymin);
		
		comp.push({name:name, type:"boundary", xmid:xmid, ymid:ymid, feature:feature, fixed:{check:true}, mask:mask, col:col, markov_branch:false});
		claa.ncomp++;
		this.determine_branching();
	}

	
	/// Checks if a compartment already exists 
	check_comp_exist(name,p)
	{
		if(name != ""){
			let sp = this.species[p];
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				if(claa.name == name){
					alertp("Cannot have a classification and compartment both '"+name+"'");
				}					
				
				if(claa.index == name){
					alertp("Cannot have a classification index and compartment both '"+name+"'");
				}		
				
				let c = find(claa.comp,"name",name);
				if(c != undefined){
					alertp("There is already a compartment with the name '"+name+"'");
					return true;
				}
			}
		}
		return false;
	}
	
	/// Selects a compartment (i.e. generate a bubble)
	select_button_bubble(lay_name,type,p,cl,i)
	{
		generate_screen();
		let lay = get_lay(lay_name);
		
		for(let b = 0; b < lay.but.length; b++){
			let bu = lay.but[b];
			
			if(bu.type == type && bu.p == p && bu.cl == cl && bu.i == i){
				select_bubble(lay.name,b);
				generate_screen();
				return;
			}
		}
		error("Could not select");
	}
	

	/// Selects a transition (i.e. generate a bubble)
	select_transition_bubble(p,cl,i)
	{
		generate_screen();
		let lay = get_lay("Transition");
		
		for(let b = 0; b < lay.but.length; b++){
			let bu = lay.but[b];
			
			if(bu.type == "Transition" && bu.p == p && bu.cl == cl && bu.i == i){
				select_bubble(lay.name,b);
				return;
			}
		}
		error("Could not select");
	}
	
	
	/// Adds a transition to the model
	add_transition(p,cl,i,f,midp,ty)
	{
		let claa = this.species[p].cla[cl];
		let tra = claa.tra;
	
		if(ty != ""){
			for(let j = 0; j < tra.length; j++){
				if(tra[j].i == i && tra[j].f == f){
					if(i == "Source") alertp("Cannot have two sources into the same compartment");
					else{
						if(f == "Sink") alertp("Cannot have two sinks leaving the same compartment");
						else{
							alertp("Cannot have two transitions with the same source and destination compartments");
						}
					}
					return false;
				}
			}
			
			if(i == f){
				alertp("The 'from' and 'to' compartments must be different");
				return false;
			}
		}
		
		let value = {};
		if(ty != ""){
			if(i == "Source"){
				let sup2 = claa.comp[f].name;
				value = this.new_trans_equations(sup2,sup2,true,p,cl);
			}
			else{	
				if(f == "Sink"){
					let sup = "("+claa.comp[i].name+"→Sink)";
					let sup2 = claa.comp[i].name;
					value = this.new_trans_equations(sup,sup2,true,p,cl);
				}
				else{
					let sup = "("+claa.comp[i].name+"→"+claa.comp[f].name+")";
					let sup2 = claa.comp[i].name;
					value = this.new_trans_equations(sup,sup2,true,p,cl);
				}
			}
		}

		tra.push({type:ty, i:i, f:f, pline:[], midp:midp, value:value});	
		this.set_transition_name(p,cl,claa.ntra);
		claa.ntra++;
		
		if(ty != "") this.determine_branching();
		
		this.update_pline(p,cl);  

		if(ty != ""){
			data_update_add_transition(p,cl);
			check_data_valid();
		}
	
		return true;
	}
	
	
	/// Creates a new set of transition equations
	new_trans_equations(sup,sup2,flag,p,cl)
	{
		let value = {};
		if(flag == false){
			value.mean_eqn = create_equation("μ^"+sup2,"trans_mean",p,cl);
			value.rate_eqn = create_equation("s^"+sup2,"trans_rate",p,cl);
		}
		else{
			value.mean_eqn = create_equation("μ^"+sup,"trans_mean",p,cl);
			value.rate_eqn = create_equation("r^"+sup,"trans_rate",p,cl);
			value.bp_eqn = create_equation("b^"+sup,"trans_bp",p,cl);
			value.shape_eqn = create_equation("k^"+sup2,"trans_shape",p,cl);
			value.scale_eqn = create_equation("ν^"+sup2,"trans_scale",p,cl);
			value.cv_eqn = create_equation("cv^"+sup2,"trans_cv",p,cl);
			value.shape_erlang = {te:"1"};
		}
		
		return value;
	}
	
	
	/// Adds a derived quantity to the model
	add_modify_derived(ob,mod)
	{
		extract_equation_properties(ob.eqn1);
		extract_equation_properties(ob.eqn2);
		
		let par = ob.eqn1.param;
		if(par.length != 1){ 
			if(par.length == 0){ alertp("Does not contain a parameter"); return;}
			if(par.length > 1){	alertp("Should not contain more than one parameter"); return;}
		}
		
		let warn = check_derived_param(ob.eqn2,ob.eqn1);
		if(warn != "success"){ alertp(warn); return;}
		
		let name = par[0].name;
		let i = find(this.param,"name",name);
		if(i != undefined && !(mod == "modify" && par[0].full_name == this.derive[ob.val].eqn1.param[0].full_name)){
			alertp("The parameter '"+name+"' is already used in the model"); 
			return;
		}
	
		if(mod == "modify") this.derive[ob.val] = ob;
		else this.derive.push(ob);
		init_param();
		return "success";
	}
	
	
	/// Determines if compartments/transition are branching
	determine_branching()
	{
		for(let p = 0; p < this.species.length; p++){
			let sp = this.species[p];
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				for(let tr = 0; tr < claa.ntra; tr++){
					claa.tra[tr].branch = false;
				}

				for(let c = 0; c < claa.ncomp; c++){
					let co = claa.comp[c];
					
					let list=[]; 
					let flag = false;
					for(let tr = 0; tr < claa.ntra; tr++){
						let trr = claa.tra[tr];
						if(trr.i == c){
							list.push(tr);
							if(trr.type != "exp(rate)" && trr.type != "exp(mean)" && tr.type != "erlang") flag = true;
						}
					}
					
					if(list.length > 1 && flag == false) co.choose_branch = true;
					else co.choose_branch = false;
				
					if(list.length > 1 && (flag == true || co.markov_branch == true)){
						co.branch = true;
					}
					else co.branch = false;
					
					for(let i = 0; i < list.length; i++){
						claa.tra[list[i]].branch = co.branch;
					}
				}					
			}
		}
	}
	

	/// Sets the name of a transition
	set_transition_name(p,cl,tr)
	{
		let claa = this.species[p].cla[cl];
		let tra = claa.tra[tr];
		
		if(tra.type != ""){
			let variety = "Normal", name;
			let i = tra.i, f = tra.f;
			if(i == "Source"){ variety = "Source"; name = "+→"+claa.comp[f].name;}
			else{
				if(f == "Sink"){ variety = "Sink"; name = claa.comp[i].name+"→-";}
				else name = claa.comp[i].name+"→"+claa.comp[f].name;
			}
			tra.variety = variety; tra.name = name;
		}
	}


	/// Sets the size of a compartment (depending on the length of the text which defined it)
	set_compartment_size(comp)
	{
		comp.h  = compartment_height
		
		let font = get_font(compartment_height*0.85,"","times");
		let w = text_width(comp.name,font)+1;
		if(w < 3) w = 3;
		
		w = (Math.floor(w/grid)+1)*grid;
		comp.w = w;
	}


	/// Adds the lower buttons (e.g. for adding compartments etc...) 
	add_lower_menu_buts(lay)
	{
		let message = "";
		switch(inter.mode.type){
		case "Add_Compartment": message = "Place compartment on screen"; break;
	
		case "Add_Transition":
			if(inter.mode.i == undefined) message = "Select initial compartment"; 
			else message = "Select final compartment (with optional intermediate points)";
			break;
		
		case "Add_Source":
			if(inter.mode.midp.length == 0) message = "Select source position"; 
			else message = "Select compartment individuals entering (with optional intermediate points)";
			break;
			
		case "Add_Sink":
			if(inter.mode.midp.length == 0) message = "Select sink position"; 
			else message = "Select compartment individuals leaving (with optional intermediate points)";		
			break;
			
		case "Add_Label": message = "Place label on screen"; break;
		
		case "Add_Box": message = "Drag to select a box containing compartments"; break;
		}
		
		if(message != ""){ 		
			lay.add_button({te:message, x:0, y:0, dx:lay.dx, dy:lay.dy, type:"Message"});
			lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, ac:"Nothing", type:"Nothing"});
			lay.add_button({te:message, x:lay.dx-1.5, y:0.2, dx:1, dy:1, ac:"MessageClose", type:"MessageClose"});
			return;
		}
		
		let x = 1.2, y = lay.dy-1.6;
		
		let pag = inter.page[inter.pa];
		let p = pag.sub[pag.index].index;

		let n = 0; if(this.species.length > 0) n = this.species[p].ncla;

		let cl = pag.sub[pag.index].sub[p].index;
		
		let claa; if(cl < n) claa = this.species[p].cla[cl];
		
		let gap = 4;
		let w;		
		
		let active;
		if(cl < n) active = true; else active = false;
		let coord; if(cl < n) coord = claa.camera.coord;
		w = this.add_object_button(lay,"Compartment",x,y,"AddCompartment",{active:active, title:"Add a compartment", te:compartment_text, back:WHITE}); x += w+gap;
	
		if(cl < n && claa.ncomp > 1) active = true; else active = false;
		w = this.add_object_button(lay,"Transition",x,y,"AddTransition",{active:active, title:"Add a transition", te:trans_text, back:WHITE}); x += w+gap;
		
		if(cl < n && claa.ncomp > 0) active = true; else active = false;
		w = this.add_object_button(lay,"Source",x,y,"AddSource",{active:active, title:"Add a source", te:source_text, back:WHITE}); x += w+gap;
		
		w = this.add_object_button(lay,"Sink",x,y,"AddSink",{active:active, title:"Add a sink", te:sink_text, back:WHITE}); x += w+gap;
	
		x += gap;
		if(cl < n) active = true; else active = false;
		w = this.add_object_button(lay,"Annotation",x,y,"AddAnnotation",{active:active, p:p, cl:cl, title:"Add an annotation", te:annotation_text, back:WHITE});
		x += w+gap;
	
		let ac_in, ac_out;
	
		if(cl < n && (claa.ncomp > 0 || claa.annotation.length > 0)){ 
			ac_in = "ZoomIn"; ac_out = "ZoomOut";
		}
		
		lay.add_button({x:lay.dx-5, y:0, dx:1.8, dy:2, ac:ac_in, type:"ZoomIn", p:p, cl:cl});

		lay.add_button({x:lay.dx-3, y:0, dx:1.8, dy:1.8, ac:ac_out, type:"ZoomOut", p:p, cl:cl});
	}


	/// Generates the popup bubble which appears on the add classification button 
	add_classification_popup()
	{
		let i = get_lay("UpperMenu").search_button_ac("AddClassification");

		select_bubble("UpperMenu",i,{type:"Popup"});
		generate_screen();
	}
	
	
	/// Generates the popup bubble which appears on the add compartment button 
	add_compartment_popup()
	{
		let i = get_lay("LowerMenu").search_button_ac("AddCompartment");

		select_bubble("LowerMenu",i,{type:"Popup"});
		generate_screen();
	}
	

	/// Adds the popup bubble which appears next to the add species button 
	add_species_popup()
	{
		let i = get_lay("Menu").search_button_ac("AddSpecies");
		select_bubble("Menu",i,{type:"Popup"});
		generate_screen();
	}
	
	
	/// Adds button on the top menu (e.g. tabs for the different classifications)
	add_upper_menu_buts(lay)
	{			
		let pag = inter.page[inter.pa];
		let p = pag.sub[pag.index].index;
		let cl_sel = pag.sub[pag.index].sub[p].index;
		
		if(this.species == 0) return;
		
		let cla = this.species[p].cla;
		let x = 0.8, y = 0;

		for(let cl = 0; cl < cla.length; cl++){
			let si = si_clatab;
			let font = get_font(si);
			let te = cla[cl].name;
			let font_it = get_font(si_claindextab,"italic","times");
			let index = cla[cl].index;
			let w = text_width(te,font)+text_width(index,font_it)+2.5;
			
			lay.add_button({te:te, index:index, x:x, y:y, dx:w, dy:1.3, ac:"ClassTab", type:"ClassTab", si:si, font:font,font_it:font_it, p:p, cl:cl, cl_sel:cl_sel});
			
			x += w+1.5;
		}
	
		x -= 0.5;
		let yy = y+0.1;
		if(cla.length == 0){ x = 0; yy = y+0.5;}
		this.add_object_button(lay,"Classification",x,yy,"AddClassification",{active:true, p:p, title:"Add a classification", te:class_text, back:WHITE});
	}


	/// Adds an add object button (e.g. add compartment) 
	add_object_button(lay,te,x,y,ac,op)
	{
		if(op == undefined) op = {};
		if(op.active == false) ac = undefined;
		
		let si = 0.8;
		let font = get_font(si);
		let w = text_width(te,font);
		lay.add_button({te:te, x:x, y:y, dx:w+1.2, dy:1.1, ac:ac, type:"AddButton", si:si, font:font, op:op});	
		
		if(ac && op && op.title){
			lay.add_help_button(x+w+1.2,y+1.4,{title:op.title, te:op.te, back_col:op.back_col});
		}
		
		return w;
	}


	/// Ensures the positions of the midpoint of compartments an on a grid
	//snap_comp_to_grid(p,cl,i)
	snap_comp_to_grid(ob)
	{
		if(ob.x == "auto") return;
		
		ob.x = Math.round(ob.x/grid)*grid;
		ob.y = Math.round(ob.y/grid)*grid;
	}


	/// Deletes a transition  point
	delete_transitionpoint(p,cl,i,index)
	{
		this.species[p].cla[cl].tra[i].midp.splice(index,1);
		this.update_pline(p,cl); 
	}
	
	
	/// Creates a new transition point from dragging transition
	split_transitionpoint()
	{
		let mo = inter.mode;
		this.species[mo.p].cla[mo.cl].tra[mo.i].midp.splice(mo.index+1,0,this.get_mouse_to_desktop());
		this.update_pline(mo.p,mo.cl);
		
		inter.mode={ type:"Drag_TransitionPoint", p:mo.p, cl:mo.cl, i: mo.i, index: mo.index+1};
	}


	/// Updates all the plines
	update_pline_all()
	{
		for(let p = 0; p < this.species.length; p++){
			let sp = this.species[p];
			for(let cl = 0; cl < sp.ncla; cl++){
				this.update_pline(p,cl);
			}
		}
	}
	
	
	/// Updates all plines in a classification
	update_pline(p,cl)
	{
		for(let i = 0; i < this.species[p].cla[cl].ntra; i++){
			this.find_trans_pline(p,cl,i);
		}
	}
	
	
	/// Converts points in a transition to points on the model
	find_trans_pline(p,cl,i)                                 
	{
		let claa = this.species[p].cla[cl];
		let cam = claa.camera;
		let r_latlng = latlng_radius*cam.ruler;
		
		let tr = claa.tra[i];

		tr.pair = false;                        // Finds if a transition is pair with its reverse
		for(let j = 0; j < claa.ntra; j++){
			if(j != i){
				let tr2 = claa.tra[j];
				if(tr.i == tr2.f && tr.f == tr2.i){
					if(tr.midp.length == 0 && tr2.midp.length == 0){
						tr.pair = true;
					}
				}
			}
		}
		
		tr.pline=[];
		
		let pair_shift = {x:0, y:0};
		let co = claa.comp[tr.i];
		
		if(tr.pair == true){
			let pi = this.comp_center(claa.comp[tr.i]);
			let pf = this.comp_center(claa.comp[tr.f]);
			
			let dx = pf.x-pi.x, dy = pf.y-pi.y;
			let r = Math.sqrt(dx*dx+dy*dy);
			
			let shift_gap = 0.3;
			
			pair_shift.x = shift_gap*dy/r;
			pair_shift.y = -shift_gap*dx/r;
		}			
		
		if(tr.i != undefined && tr.i != "Source"){
			let co = claa.comp[tr.i];
			
			let p1 = this.comp_center(co);
		
			p1.x += pair_shift.x; p1.y += pair_shift.y;
			
			let	p2;
			if(tr.midp.length > 0) p2 = tr.midp[0];
			else{
				if(tr.f != "Sink" && tr.f != undefined) p2 = this.comp_center(claa.comp[tr.f]);
			}
			
			if(p2 != undefined){
				let pos = this.intersect_comp_line(co,p1,p2,r_latlng);
				tr.pline.push({x:pos.x, y:pos.y, index:-1});
			}
		}
		
		for(let i = 0; i < tr.midp.length; i++){
			tr.pline.push({x:tr.midp[i].x, y:tr.midp[i].y, index:i});
		}
		
		if(tr.f != undefined && tr.f != "Sink"){
			let co = claa.comp[tr.f];
			
			let p1 = this.comp_center(co);
		
			p1.x += pair_shift.x; p1.y += pair_shift.y;
			
			let p2;
			if(tr.midp.length > 0) p2 = tr.midp[tr.midp.length-1];
			else{ if(tr.i != "Source") p2 = this.comp_center(claa.comp[tr.i]);}
		
			let pos = this.intersect_comp_line(co,p1,p2,r_latlng);
			tr.pline.push({x:pos.x, y:pos.y, index:tr.midp.length});
		}
		

		let np = tr.pline.length;
		let l = 0;
		let dist=[];
		for(let p = 0; p < np-1; p++){
			let dx = tr.pline[p+1].x - tr.pline[p].x, dy = tr.pline[p+1].y - tr.pline[p].y;
			dist[p] = Math.sqrt(dx*dx+dy*dy)+TINY; 
			l += dist[p];
		}

		let l2 = 0; 
		for(let p = 0; p < np-1; p++){
			let dx = tr.pline[p+1].x - tr.pline[p].x, dy = tr.pline[p+1].y - tr.pline[p].y;
			let dl = dist[p];
			if(l2+dl > l/2){
				let f = (l/2-l2)/dl;
				let nx = dx/dist[p];
				let ny = dy/dist[p];
				tr.center = { x:tr.pline[p].x+f*dx, y:tr.pline[p].y+f*dy, nx:nx, ny:ny};
				break;
			}
			l2 += dl;
		}
		tr.l = l;
		
		if(l > 0){
			let wmin = LARGE, wmax = -LARGE; // Finds the potential width for the label
			for(let p = 0; p < np; p++){
				let dx = tr.pline[p].x - tr.center.x;
				let dy = tr.pline[p].y - tr.center.y;
				
				let d = dx*tr.center.nx + dy*tr.center.ny;
				if(d < wmin) wmin = d; if(d > wmax) wmax = d;
			}
			tr.center.w = wmax-wmin;
		}
		
		let xmax = -LARGE, xmin = LARGE;   // Finds a box around the transition points
		let ymax = -LARGE, ymin = LARGE;
			
		for(let j = 0; j < tr.pline.length; j++){
			let po = tr.pline[j];
			if(po.x < xmin) xmin = po.x; if(po.x > xmax) xmax = po.x;
			if(po.y < ymin) ymin = po.y; if(po.y > ymax) ymax = po.y;	
		}
		tr.box = {x:xmin, y:ymin, dx:xmax-xmin, dy:ymax-ymin};
		
		let lab="";
	
		let val = tr.value;

		switch(tr.type){
		case "exp(mean)": 
			if(tr.branch == true) lab = val.bp_eqn.te+"/";
			else lab = "1/";
			lab += "("+val.mean_eqn.te+")"; 
			break;
			
		case "exp(rate)": 
			if(tr.branch == true) lab = val.bp_eqn.te+"×"; 
			lab += val.rate_eqn.te; 
			break;
			
		case "gamma": 
			if(tr.branch == true) lab = val.bp_eqn.te+","; 
			lab += "Γ("+val.mean_eqn.te+" , "+val.cv_eqn.te+")";
			break;
			
		case "erlang": 
			if(tr.branch == true) lab = val.bp_eqn.te+","; 
			lab += "EL("+val.mean_eqn.te+" , "+val.shape_erlang.te+")";
			break;
			
		case "log-normal":
			if(tr.branch == true) lab = val.bp_eqn.te+","; 
			lab += "LN("+val.mean_eqn.te+" , "+val.cv_eqn.te+")"; 
			break;
			
		case "weibull":
			if(tr.branch == true) lab = val.bp_eqn.te+","; 
			lab += "W("+val.scale_eqn.te+" , "+val.shape_eqn.te+")";
			break;
			
		default: delete tr.label; break;
		}
		
		this.set_transition_name(p,cl,i);
	
		tr.label = lab;
	}
	

	/// 	Returns the center of a compartment
	comp_center(co)
	{
		switch(co.type){
		case "box": case "latlng": return {x:co.x, y:co.y};
		case "boundary": return {x:co.xmid, y:co.ymid};
		default: error("Option not recognised 85"); break;
		}
	}
	
	
	/// Finds the intersection between a rounded compartment and a line 
	intersect_comp_line(co,p1,p2,r_latlng)
	{
		switch(co.type){
		case "box":
			{
				let lx=[], ly=[], jmax = 4;
							
				let dx = co.w, dy = co.h;
				let xx = co.x-dx/2, yy = co.y-dy/2;
		
				let nl = 0;
				let r = COMPARTMENT_CURVE*dy;
				for(let j = 0; j <= jmax; j++){ let th = j*Math.PI/(2*jmax); lx[nl] = xx+dx-r+r*Math.sin(th); ly[nl] = yy+r-r*Math.cos(th); nl++;}
				for(let j = 0; j <= jmax; j++){ let th = j*Math.PI/(2*jmax); lx[nl] = xx+dx-r+r*Math.cos(th); ly[nl] = yy+dy-r+r*Math.sin(th); nl++;}
				for(let j = 0; j <= jmax; j++){ let th = j*Math.PI/(2*jmax); lx[nl] = xx+r-r*Math.sin(th); ly[nl] = yy+dy-r+r*Math.cos(th); nl++;}
				for(let j = 0; j <= jmax; j++){ let th = j*Math.PI/(2*jmax); lx[nl] = xx+r-r*Math.cos(th); ly[nl] = yy+r-r*Math.sin(th); nl++;}

				let nx = p2.x - p1.x, ny = p2.y - p1.y;

				for(let ii = 0; ii < nl; ii++){
					let x2 = lx[ii], y2 = ly[ii];
					let nx2 = lx[(ii+1)%nl]-x2, ny2 = ly[(ii+1)%nl]-y2;

					r = nx*ny2-ny*nx2;
					if(r != 0){
						let al = ((x2-p1.x)*ny2-(y2-p1.y)*nx2)/r;
						let be = ((x2-p1.x)*ny-(y2-p1.y)*nx)/r;
						if(al >= 0 && al <= 1 && be >= 0 && be <= 1){
							nx *= al; ny *= al;
						}
					}
				}
				return {x:p1.x+nx, y:p1.y+ny};
			}
		
		case "boundary":
			return {x:co.xmid, y:co.ymid};
			
		case "latlng":
			{
				let po;
				if(p1.x == co.x && p1.y == co.y) po = p2;
				if(p2.x == co.x && p2.y == co.y) po = p1;
				let ddx = po.x-co.x, ddy = po.y-co.y;
				
				let r = Math.sqrt(ddx*ddx + ddy*ddy);
		
				if(r < r_latlng) return {x:po.x, y:po.y};
				else{
					let fr = r_latlng/r;
					return {x:co.x+fr*ddx, y:co.y+fr*ddy};
				}
			}
			break;

		default: error("Option not recognised 87"); break;
		}
	}
	
	
	/// Deletes a specified compartment
	delete_compartment(p,cl,i,include_clone)
	{
		let claa = this.species[p].cla[cl];
		
		// First removes any transitions connected with compartment
		for(let	tr = claa.ntra-1; tr >= 0; tr--){
			let traa = claa.tra[tr];
			if(traa.i == i || traa.f == i) this.delete_transition(p,cl,tr);
		}

		// Removes compartment from any box annotations
		let name = claa.comp[i].name;
		let k = 0;
		while(k < claa.annotation.length){
			let an = claa.annotation[k];
			if(an.type == "box"){
				let j = 0; 
				while(j < an.comps.length){
					if(an.comps[j] == name) an.comps.splice(j,1);
					else j++;
				}
			}
			
			if(an.type == "box" && an.comps.length == 0) claa.annotation.splice(k,1);
			else k++;
		}
	
		claa.comp.splice(i,1);
		claa.ncomp--;

		for(let	tr = 0; tr < claa.ntra; tr++){
			let traa = claa.tra[tr];
			if(traa.i > i) traa.i--;
			if(traa.f > i) traa.f--;
		}

		this.determine_branching();

		this.update_pline(p,cl);

		data_update_delete_compartment(p,cl,i,name);
		
		if(include_clone == true) this.update_param_index(claa.index); // Updates param
		
		if(include_clone == true){  // Adds to clones 
			let list = this.find_clones(p,cl);
		
			for(let li of list){
				this.delete_compartment(li.p,li.cl,i,false);
			}
			
			check_data_valid();
		}
	}


  /// Deletes a specified transition
	delete_transition(p,cl,i)                             
	{
		let claa = this.species[p].cla[cl];
		
		claa.tra.splice(i,1);
		claa.ntra--;
		this.determine_branching();
		this.update_pline(p,cl); 

		data_update_delete_transition(p,cl,i);
		check_data_valid();
	}


	/// Changes the name of a compartment
	rename_compartment(p,cl,c,new_name)
	{
		let sp = this.species[p];
		let p_name = sp.name;
		let claa = sp.cla[cl];
		let cl_name = claa.name;
		let old_name = claa.comp[c].name;
		
		let eq_list = this.find_equation_list();
		
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
		clone_camera(p,cl);
	
		// this converts all classification on all species
		for(let p = 0; p < this.species.length; p++){
			let sp = this.species[p];
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				if(claa.name == cl_name){
					for(let c = 0; c < claa.ncomp; c++){
						let co = claa.comp[c];
						if(co.name == old_name){
							co.name = new_name;
							this.set_compartment_size(co);
							this.update_pline(p,cl);
						}
					}
				}
			}
		}

		// Changes name on param
		for(let i = 0; i < this.param.length; i++){
			let par = this.param[i];
			
			let flag = false;
			for(let j = 0; j < par.dep.length; j++){
				if(par.dep[j] == claa.index){
					if(c < par.list[j].length){
						if(par.list[j][c] != old_name) error("Old name not correct");
						par.list[j][c] = new_name;
						flag = true;
					}
				}
			}
			
			if(flag == true){
				par.comb_list = generate_comb_list(par.list);
			}
		}
		
		// Changes any names on box annotations
		for(let k = 0; k < claa.annotation.length; k++){
			let an = claa.annotation[k];
			if(an.type == "box"){
				for(let k = 0; k < an.comps.length; k++){
					if(an.comps[k] == old_name) an.comps[k] = new_name;
				}
			}
		}
	
		data_update_rename_compartment(p,cl,old_name,new_name);
		
		this.rename_ob(model,"model","comp_name",old_name,new_name);
		this.rename_ob(model,"model","comp_name_store",old_name,new_name);
	
		if(false) this.check_ob_string_exist(model,"model",old_name);
			
		this.update_pline_all();

		check_data_valid();
	}


	/// Rename a classification
	rename_classification(new_name,p,cl)
	{
		let p_name = this.species[p].name;
		
		let old_name = this.species[p].cla[cl].name
		if(new_name == old_name) return;

		let eq_list = this.find_equation_list();
		
		let len = old_name.length;
		let dif = new_name.length - old_name.length;
		
		for(let i = 0; i < eq_list.length; i++){
			let eqn = eq_list[i];
			extract_equation_properties(eqn);
		
			if(eqn.cl_name == old_name) eqn.cl_name = new_name;
		}
	
		// this converts all classification on all species
		for(let p = 0; p < this.species.length; p++){
			let sp = this.species[p];
			for(let cl = 0; cl < sp.ncla; cl++){
				if(this.species[p].cla[cl].name == old_name){
					this.species[p].cla[cl].name = new_name;
				}
			}
		}

		data_update_rename_classification(p,cl,old_name,new_name);
		
		this.rename_ob(model,"model","cl_name",old_name,new_name);
		this.rename_ob(model,"model","cl_name_store",old_name,new_name);
		
		if(false) this.check_ob_string_exist(model,"model",old_name);
		
		this.update_pline_all();
		
		check_data_valid();
	}
	
	/// Rename a species
	rename_species(new_name,p)
	{	
		let old_name = this.species[p].name;
		if(new_name == old_name) return;

		let eq_list = this.find_equation_list();
		
		let len = old_name.length;
		let dif = new_name.length - old_name.length;
		
		for(let i = 0; i < eq_list.length; i++){
			let eqn = eq_list[i];
			extract_equation_properties(eqn);
			if(eqn.p_name == old_name) eqn.p_name = new_name;
			
			let te = eqn.te;
			for(let j = 0; j < eqn.sp_name_list.length; j++){
				let ch = eqn.sp_name_list[j];
				if(ch.p_name == old_name){
					te = te.substr(0,ch.icur)+new_name+te.substr(ch.icur+len);
					for(let jj = j+1; jj < eqn.sp_name_list.length; jj++){
						eqn.sp_name_list[jj].icur += dif;
					}
				}
			}
			eqn.te = te;
		}
	
		this.species[p].name = new_name;
		
		this.rename_ob(model,"model","p_name",old_name,new_name);
	
		if(false) this.check_ob_string_exist(model,"model",old_name);
		
		this.update_pline_all();
		
		check_data_valid();
	}
	
	
	/// Renames an index on a classification
	rename_index(new_name,p,cl,include_clone)
	{
		let claa = this.species[p].cla[cl];
		
		let old_name = claa.index;
		if(new_name == old_name) return "success";
		
		claa.index = new_name;
		
		if(include_clone == true){ 
			for(let i = 0; i < this.param.length; i++){
				let par = this.param[i];
				for(let d = 0; d < par.dep.length; d++){
					if(par.dep[d] == old_name){
						par.dep[d] = new_name;
						par.label_info = get_label_info(par);
					}
				}
			}
		}
		
		if(include_clone == true){ 
			let eq_list = this.find_equation_list();
			
			for(let i = 0; i < eq_list.length; i++){
				let eqn = eq_list[i];
				extract_equation_properties(eqn);
				
				let te = eqn.te;
				for(let j = 0; j < eqn.index_name_list.length; j++){
					let ind = eqn.index_name_list[j];
					if(ind.index_name == old_name){
						ind.index_name = new_name;
						te = te.substr(0,ind.icur)+new_name+te.substr(ind.icur+1);
					}
				}
				eqn.te = te;
			}
		}
		
		if(include_clone == true){  // Adds to clones 
			let list = this.find_clones(p,cl);
		
			for(let li of list){
				this.rename_index(new_name,li.p,li.cl,false);
			}
			
			if(false) this.check_ob_string_exist(model,"model",old_name);
			
			this.update_pline_all();
		}
	}


	/// Renames a property within the model
	rename_ob(ob,root,prop,old_name,new_name)
	{
		if(old_name == "") return;
		
		for(let ele in ob){
			if(ele == prop){ 
				if(ob[ele] == old_name) ob[ele] = new_name; 
			}
			else{
				if(ele != "eqn_appear"){ 
					let ob2 = ob[ele];
					
					if(typeof ob2 == 'object'){
						this.rename_ob(ob2,root+"->"+ele,prop,old_name,new_name);
					}
				}
			}
		}
	}
	
	
	/// Checks if a string exists in an of the properties of an object (diagnostic)
	check_ob_string_exist(ob,root,name)
	{
		if(name == "") return;
		
		for(let ele in ob){
			if(ele != "eqn_appear"){ 
				let ob2 = ob[ele];
				if(ob2 == name) error("'"+name+"' found here: "+root);
				
				if(typeof ob2 == 'object'){
					this.check_ob_string_exist(ob2,root+"->"+ele,name);
				}
			}
		}
	}
	
	
	/// Checks that the name no longer exist is the model (a diagnostic check)
	check_rename(name)
	{
		let te = JSON.stringify(model);
		if(te.includes('"'+name+'"')){
			error("problem still there");
		}
	}
		
	
	/// Finds a list of branching points
	find_branch_list()
	{
		let branch_list = [];
		
		for(let p = 0; p < this.species.length; p++){
			let sp = this.species[p];
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				for(let c = 0; c < claa.ncomp; c++){
					let bl = { p:p, cl:cl, c:c, tra:[]};
					for(let i = 0; i < claa.ntra; i++){
						let tr = claa.tra[i];
						if(tr.branch == true && tr.i == c) bl.tra.push(i);
					}
					if(bl.tra.length > 0){
						if(bl.tra.length == 1) error("Cannot have one branch");
						branch_list.push(bl); 
					}
				}
			}
		}
		
		return branch_list; 
	}
	
	
	/// Finds a list of all the equations in the model
	find_equation_list(all_param)
	{
		let eqn_list = [];
		for(let p = 0; p < this.species.length; p++){
			let sp = this.species[p];
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				for(let i = 0; i < claa.ntra; i++){
					let eqn_info = {p:p,cl:cl,i:i};
					
					let traa = claa.tra[i];
				
					let val = traa.value;
					if(val.bp_eqn != undefined){
						if(all_param == true || traa.branch == true){
							this.add_equation_to_list(eqn_list,val.bp_eqn,eqn_info);
						}
					}
					
					if(val.mean_eqn != undefined){
						if(all_param == true || traa.type == "exp(mean)" || traa.type == "gamma" || traa.type == "erlang" || traa.type == "log-normal"){
							this.add_equation_to_list(eqn_list,val.mean_eqn,eqn_info);
						}
					}
					
					if(val.rate_eqn != undefined){
						if(all_param == true || traa.type == "exp(rate)"){
							this.add_equation_to_list(eqn_list,val.rate_eqn,eqn_info);
						}
					}
					
					if(val.scale_eqn != undefined){
						if(all_param == true || traa.type == "weibull"){
							this.add_equation_to_list(eqn_list,val.scale_eqn,eqn_info);
						}
					}
					
					if(val.shape_eqn != undefined){
						if(all_param == true || traa.type == "weibull"){
							this.add_equation_to_list(eqn_list,val.shape_eqn,eqn_info);
						}
					}
						
					if(val.cv_eqn != undefined){
						if(all_param == true || traa.type == "gamma" || traa.type == "log-normal"){
							this.add_equation_to_list(eqn_list,val.cv_eqn,eqn_info);
						}
					}
				}
			}
			
			for(let i = 0; i < sp.inf_source.length; i++){
				let so = sp.inf_source[i];
				
				let eqn_info = {p:p, i:i};
					
				if(so.type=="Diag. Test"){
					this.add_equation_to_list(eqn_list,so.spec.Se_eqn,eqn_info);
					this.add_equation_to_list(eqn_list,so.spec.Sp_eqn,eqn_info);
				}
				
				if(so.type=="Set Traps"){
					this.add_equation_to_list(eqn_list,so.spec.trap_prob_eqn,eqn_info);
				}
				
				if(so.type == "Compartment"){
					let tab = so.table;
					for(let r = 0; r < tab.nrow; r++){
						let te = tab.ele[r][2];
						let spl = te.split("|");
						for(let k = 0; k < spl.length; k++){
							let spl2 = spl[k].split(":");
							if(spl2.length == 2){
								let te2 = spl2[1];
								
								if(isNaN(te2)){
									let eqn_info2 = {p:p,i:i,r:r,c:2};
									
									let eqn = create_equation(te2,"comp_prob",p,undefined);  
									this.add_equation_to_list(eqn_list,eqn,eqn_info2);
								}
							}
						}
					}
				}
			}		
		}
		
		/// Goes through derived parameters
		for(let i = 0; i < this.derive.length; i++){
			let der = this.derive[i];
			
			let eqn_info = {i:i};
			
			this.add_equation_to_list(eqn_list,der.eqn1,eqn_info);		
		}
		
		/*
		/// Goes through reparameterisation
		for(let th = 0; th < this.param.length; th++){  
			let par = this.param[th];
			if(par.variety == "reparam"){		
				if(par.dep.length == 0){
					let ele = par.value;
					if(isNaN(ele)){
						let eqn = create_equation(ele,par.variety);
						let eqn_info = {par_name:par.name};
						this.add_equation_to_list(eqn_list,eqn,eqn_info);		
					}
				}
				else{
					for(let i = 0; i < par.comb_list.length; i++){
						let comb = par.comb_list[i];
						let ele = get_element(par.value,comb.index);
						if(isNaN(ele)){
							let eqn = create_equation(ele,par.variety);
							let eqn_info = {par_name:par.name, index:comb.index};
							this.add_equation_to_list(eqn_list,eqn,eqn_info);		
						}
					}
				}
			}
			
			if(par.variety == "dist"){
				if(par.prior_split_check.check == false){
					let eqn_info = {par_name:par.name};
					
					this.add_distribution_eqn(eqn_list,par.prior,eqn_info);	
				}
				else{
					for(let i = 0; i < par.comb_list.length; i++){
						let comb = par.comb_list[i];
						let ele = get_element(par.prior_split,comb.index);
					
						let eqn_info = {par_name:par.name, index:comb.index};
					
						this.add_distribution_eqn(eqn_list,ele,eqn_info);	
					}
				}
			}
		}
		*/
		
		return eqn_list;
	}
			
	
	
	
	
	/*
	/// Adds a distribution equation to equation list
	add_dist_eqn(eqn_list,eqn,eqn_info)	
	{
		if(isNaN(eqn.te)){
			this.add_equation_to_list(eqn_list,eqn,eqn_info);
		}			
	}
	*/


	/// Adds equation to a list and appends information about source 
	add_equation_to_list(eqn_list,source,eqn_info)
	{
		source.eqn_info = eqn_info;
		eqn_list.push(source);
	}
	
	
	/// Deletes a species
	delete_species(p)
	{
		this.species.splice(p,1);
		
		if(p == this.species.length && p > 0) change_page({pa:"Model",su:"Compartments",susu:p-1});
		initialise_pages();
		reset_info_p();
	}


	/// Deletes a classification
	delete_classification(p,cl)
	{
		let sp = this.species[p];
		let claa = sp.cla[cl];
		let name = claa.name;
		let index =claa.index;
		
		let comp_list = [];
		for(let i = 0; i < claa.ncomp; i++) comp_list.push(claa.comp[i].name);

		sp.cla.splice(cl,1);
		sp.ncla--;
	
		initialise_pages();

		data_update_delete_classification(p,cl,name,comp_list);
		
		for(let i = 0; i < this.param.length; i++){
			let par = this.param[i];
			
			let j = find_in(par.dep,index);
			if(j != undefined){
				par.dep.splice(j,1);
				par_set_default(par);
			}
		}
		
		check_data_valid();
	}


	/// Determines if mouse is over a transition
	mouse_over_transition(points,lay)                         
	{		
		if(inter.mode.type == "Add_Transition") return UNSET;
		
		let mx = inter.mx-lay.x, my = inter.my - lay.y;
		let np = points.length;
		for(let p = 0; p < np-1; p++){
			let x = points[p].x, y = points[p].y;
			let xx = points[p+1].x, yy = points[p+1].y;
			
			let dx = xx-x, dy = yy-y;
			let dmx = mx-x, dmy = my-y;
			let l = (dx*dmx + dy*dmy)/(dx*dx+dy*dy);
			if(l > 0 && l < 1){
				let para = (dx*dmy - dy*dmx)/Math.sqrt(dx*dx+dy*dy);
				if(para < TRANS_OVER_RANGE && para > -TRANS_OVER_RANGE) return points[p].index;
			}
		}

		return UNSET;
	}


	/// Based on the individual classifications, this generates all the compartments in the system
	get_glob_comp(p)
	{
		let sp = this.species[p];
		
		let glob_comp=[];
		
		let index=[];	
		for(let cl = 0; cl < sp.ncla; cl++){
			if(sp.cla[cl].ncomp == 0) return glob_comp;
			index[cl] = 0;
		}
		
		do{
			let co = { cla:[]};
			
			for(let cl = 0; cl < sp.ncla; cl++){
				co.cla[cl] = sp.cla[cl].comp[index[cl]].name;
			}
			glob_comp.push(co);
				
			let cl = 0; 
			let flag;
			do{
				flag = false;
				index[cl]++; if(index[cl] >= sp.cla[cl].ncomp){ index[cl] = 0; cl++; flag = true;}
			}while(flag == true && cl < sp.ncla);
			if(cl == sp.ncla) break;
		}while(true);
			
		return glob_comp;
	}
	
	
	add_warning_buts(lay)
	{
		let cx = corner.x;
		let cy = corner.y;
	
		cy = lay.add_title("Some things need fixing!",cx,cy,{te:fix_text});
		
		cy = lay.add_paragraph("The following errors need to be corrected:",lay.inner_dx-cx-1,cx,cy,BLACK,para_si,para_lh);
		cy += 0.5;
		
		add_layer("WarningContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-3.5,{});
	}
	
	
	/// Displays a screen giving warning messages
	add_warning_content(lay)
	{
		let x = 0.5;
		let y = 0.5;
		let dx = lay.inner_dx;
	
		let warn = this.warn;
		for(let i = 0; i < warn.length; i++){
			if(warn[i].mess != ""){ 
				lay.add_button({x:x, y:y+0.5, dx:2.3, dy:2, type:"WarnPic"});
			}
			else y -= 1;
				
			lay.add_button({te:"View", x:dx-5, y:y+0.9, dx:4, dy:1.5, ac:"ViewModelError", type:"GreyView", i:i});
	
			y = lay.add_paragraph("<b>"+warn[i].mess+"</b>",dx-8,x+3,y+0.4,BLACK,para_si,para_lh);
			
			let yy = y; if(warn[i].mess == "") yy -= 0.4;
			
			y = lay.add_paragraph(warn[i].mess2,dx-8,x+3,yy,BLACK,warn_si,warn_lh);
			
			y += 0.5;
		}
	}
	

	/// Adds compartments from a table
	add_file_compartment(p,cl,tab,col)
	{
		let claa = this.species[p].cla[cl];	
		let cam = claa.camera;

		let pos_x=[], pos_y=[];
	
		if(tab.ncol <= 2){ // Makes up a set of points on a grid
			let font = get_font(compartment_height*0.85,"","times");
		
			let w_max = 0;
			for(let r = 0; r < tab.nrow; r++){
				let w = text_width(tab.ele[r][0],font)+1;
				if(w > w_max) w_max = w;
			}
			if(w_max < 3) w_max = 3;
			
			let gapx = 4, gapy = 4;
			
			let ratio = (page_char_wid-menu_width)/page_char_hei;
			
			ratio *= (compartment_height+gapy)/(w_max+gapx);
			
			let d = 0.7*Math.sqrt(tab.nrow/ratio);
			let LX, LY;
			do{
				d *= 1.01;
				LX = Math.floor(ratio*d); LY = Math.floor(d); 
			}while(LX*LY < tab.nrow);
			
			for(let j = 0; j < LY; j++){
				for(let i = 0; i < LX; i++){
					let r = j*LX+i;
					if(r < tab.nrow){
						pos_x[r] = i*(w_max+gapx);
						pos_y[r] = j*(compartment_height+gapy);
					}
				}
			}
		}
		else{
			for(let r = 0; r < tab.nrow; r++){
				switch(cam.coord){
				case "cartesian":								
					pos_x[r] = Number(tab.ele[r][1])*import_scale_factor;
					pos_y[r] = Number(tab.ele[r][2])*import_scale_factor;
					break;
					
				case "latlng":
					let pt = transform_latlng(tab.ele[r][2],tab.ele[r][1]);
					pos_x[r] = pt.x; pos_y[r] = pt.y; 
					break;

				default: error("Option not recognised 88"); break;
				}
			}
		}
		
		for(let r = 0; r < tab.nrow; r++){	
			let col2 = col; if(tab.ncol == 2 || tab.ncol == 4) col2 = tab.ele[r][tab.ncol-1];
			let name = tab.ele[r][0];
			
			let x = pos_x[r], y = pos_y[r];
			
			switch(cam.coord){
			case "cartesian":		
				{		
					let c = find(claa.comp,"name",name);
					if(c == undefined) this.add_compartment(name,p,cl,x,y,col2,true,true);		
					else{
						let co = claa.comp[c];
						co.x = x; co.y = y; co.col = col2;
						this.set_compartment_size(co);
					}
				}
				break;
			
			case "latlng":
				{
					let c = find(claa.comp,"name",name);
					if(c == undefined) this.add_latlng_compartment(name,p,cl,x,y,col2,true);		
					else{
						let co = claa.comp[c];
						co.x = x; co.y = y; co.col = col2;
					}
				}
				break;

			default: error("Option not recognised 89"); break;
			}
		}

		set_camera(p,cl);
	}
	
	
	/// Adds compartments from a table
	add_file_transition(p,cl,tab,type)
	{
		let claa = this.species[p].cla[cl];	
		for(let r = 0; r < tab.nrow; r++){	
			let i	= find(claa.comp,"name",tab.ele[r][0]);
			let f	= find(claa.comp,"name",tab.ele[r][1]);
			
			if(i != undefined && f != undefined && i != f){
				let k = 0; 
				while(k < claa.tra.length && !(claa.tra[k].i == i && claa.tra[k].f == f)) k++;
				
				let value;
				if(k < claa.tra.length){
					let traa = claa.tra[k];
					traa.type = type;
					traa.midp = [];
					value = traa.value;
				}
				else{
					this.add_transition(p,cl,i,f,[],type);
					value = claa.tra[claa.tra.length-1].value;
				}
				
				switch(type){
				case "exp(rate)":
					value.rate_eqn.te = tab.ele[r][2];
					break;
					
				case "exp(mean)": 
					value.mean_eqn.te = tab.ele[r][2];
					break;
					
				case "erlang": 
					value.mean_eqn.te = tab.ele[r][2];
					value.shape_eqn.te = tab.ele[r][3];
					break;
					
				case "gamma":
					value.mean_eqn.te = tab.ele[r][2];
					value.cv_eqn.te = tab.ele[r][3];
					break;
					
				case "log-normal": 
					value.mean_eqn.te = tab.ele[r][2];
					value.cv_eqn.te = tab.ele[r][3];
					break;
					
				case "weibull": 
					value.scale_eqn.te = tab.ele[r][2];
					value.shape_eqn.te = tab.ele[r][3];
					break;

				default: error("Option not recognised 90"); break;
				}
			}
		}
		
		this.determine_branching();
		this.update_pline(p,cl);   
			
		set_camera(p,cl);
	}
	
	
	/// Adds sources from a table
	add_file_source(p,cl,tab,type)
	{
		let claa = this.species[p].cla[cl];	
		let cam = claa.camera;
		
		let midp_list=[];
	
		if(tab.ncol == 2){ 
			for(let r = 0; r < tab.nrow; r++){	
				let f	= find(claa.comp,"name",tab.ele[r][0]);
				if(f != undefined){
					let co = claa.comp[f];
					let x, y;
					switch(cam.coord){
					case "cartesian":	
						x = co.x;
						y = co.y+compartment_height+1;
						break;
						
					case "latlng":
						x = co.x;
						y = co.y;
						break;

					default: error("Option not recognised 91"); break;
					}
					midp_list[r] = [{x:x,y:y}];
				}
			}
		}
		else{
			for(let r = 0; r < tab.nrow; r++){	
				let x,y;
				switch(cam.coord){
				case "cartesian":								
					x = Number(tab.ele[r][1])*import_scale_factor;
					y = Number(tab.ele[r][2])*import_scale_factor;
				break;
				
				case "latlng":
					let pt = transform_latlng(tab.ele[r][2],tab.ele[r][1]);
					x = pt.x; y = pt.y; 
					break;

				default: error("Option not recognised 92"); break;
				}
				
				midp_list[r] = [{x:x,y:y}];
			}
		}
		
		for(let r = 0; r < tab.nrow; r++){	
			let f	= find(claa.comp,"name",tab.ele[r][0]);
			
			if(f != undefined){
				let k = 0; 
				while(k < claa.tra.length && !(claa.tra[k].i == "Source" && claa.tra[k].f == f)) k++;
				
				let value;
				if(k < claa.tra.length){
					let traa = claa.tra[k];
					traa.type = type;
					traa.midp = midp_list[r];
					value = traa.value;
				}
				else{
					this.add_transition(p,cl,"Source",f,midp_list[r],type);
					value = claa.tra[claa.tra.length-1].value;
				}
				
				switch(type){
				case "exp(rate)":
					value.rate_eqn.te = tab.ele[r][tab.ncol-1];
					break;
					
				case "exp(mean)": 
					value.mean_eqn.te = tab.ele[r][tab.ncol-1];
					break;

				default: error("Option not recognised 93"); break;
				}
			}
		}
		
		this.determine_branching();
		this.update_pline(p,cl);   
			
		set_camera(p,cl);
	}
	
	
	/// Adds sinks from a table
	add_file_sink(p,cl,tab,type)
	{
		let claa = this.species[p].cla[cl];	
		let cam = claa.camera;
		
		let midp_list=[];
	
		if(tab.ncol == 2){ 
			for(let r = 0; r < tab.nrow; r++){	
				let i	= find(claa.comp,"name",tab.ele[r][0]);
				if(i != undefined){
					let co = claa.comp[i];
					let x, y;
					
					switch(cam.coord){
					case "cartesian":	
						x = co.x;
						y = co.y+compartment_height+1;
						break;
						
					case "latlng":
						x = co.x;
						y = co.y;
						break;

					default: error("Option not recognised 95"); break;
					}
					midp_list[r] = [{x:x,y:y}];
				}
			}
		}
		else{
			for(let r = 0; r < tab.nrow; r++){	
				let x,y;
				switch(cam.coord){
				case "cartesian":								
					x = Number(tab.ele[r][1])*import_scale_factor;
					y = Number(tab.ele[r][2])*import_scale_factor;
				break;
				
				case "latlng":
					let pt = transform_latlng(tab.ele[r][2],tab.ele[r][1]);
					x = pt.x; y = pt.y; 
					break;

				default: error("Option not recognised 96"); break;
				}
				
				midp_list[r] = [{x:x,y:y}];
			}
		}
		
		for(let r = 0; r < tab.nrow; r++){	
			let i	= find(claa.comp,"name",tab.ele[r][0]);
			
			if(i != undefined){
				let k = 0; 
				while(k < claa.tra.length && !(claa.tra[k].i == i && claa.tra[k].f == "Sink")) k++;
				
				let value;
				if(k < claa.tra.length){
					let traa = claa.tra[k];
					traa.type = type;
					traa.midp = midp_list[r];
					value = traa.value;
				}
				else{
					this.add_transition(p,cl,i,"Sink",midp_list[r],type);
					value = claa.tra[claa.tra.length-1].value;
				}
				
				switch(type){
				case "exp(rate)":
					value.rate_eqn.te = tab.ele[r][tab.ncol-1];
					break;
					
				case "exp(mean)": 
					value.mean_eqn.te = tab.ele[r][tab.ncol-1];
					break;
					
				case "gamma": 
					value.mean_eqn.te = tab.ele[r][tab.ncol-2];
					value.cv_eqn.te = tab.ele[r][tab.ncol-1];
					break;
					
				case "erlang": 
					value.mean_eqn.te = tab.ele[r][tab.ncol-2];
					value.shape_eqn.te = tab.ele[r][tab.ncol-1];
					break;
					
				case "log-normal": 
					value.mean_eqn.te = tab.ele[r][tab.ncol-2];
					value.cv_eqn.te = tab.ele[r][tab.ncol-1];
					break;
					
				case "weibull": 
					value.scale_eqn.te = tab.ele[r][tab.ncol-2];
					value.shape_eqn.te = tab.ele[r][tab.ncol-1];
					break;

				default: error("Option not recognised 98"); break;
				}
			}
		}
		
		this.determine_branching();
		this.update_pline(p,cl);   
			
		set_camera(p,cl);
	}

	
	/// Adds a classification background buttons (this is activated when clicked)
	add_classback_buts(lay)
	{
		if(this.species.length == 0) return;
		let p = this.get_p();
		if(this.species[p].ncla == 0) return;
		let cl = this.get_cl();
		lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, ac:"ClassificationBack", type:"Nothing", p:p, cl:cl});
	}
	
	
	/// Replots the model
	replot()
	{
		if(find(inter.layer,"name","Annotation") == undefined){
			replot_layer("GraphCompartments");
			replot_layer("GraphTransitions");
		}			
		else{
			replot_layer("Annotation");
			replot_layer("AnnotationMap");
			replot_layer("Compartment");
			replot_layer("Transition");
		}
		
		plot_screen();
	}
	
	
	/// Gets all the values an index cound take
	get_index_possibility(def_pos)
	{
		let pos = [];
		for(let i = 0; i < alphabet.length; i++){
			let ch = alphabet[i];
			let flag = false;
			if(ch == "a" || ch == "t") flag = true;
			
			if(ch != def_pos){
				for(let p = 0; p < this.species.length; p++){
					let sp = this.species[p];
					for(let cl = 0; cl < sp.ncla; cl++){
						let claa = sp.cla[cl];
						
						if(claa.index == ch) flag = true;
						for(let c = 0; c < claa.ncomp; c++){
							if(claa.comp[c].name == ch) flag = true;
						}
					}
				}
			}
			if(flag == false) pos.push({te:ch});
		}
		
		return pos;
	}

	/// Check clones all have the same specifications for compartments
	check_clones()
	{
		if(check_clone_on != true) return;
	
		for(let p = 0; p < this.species.length; p++){
			let sp = this.species[p];
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				
				for(let p2 = p+1; p2 < this.species.length; p2++){
					let sp2 = this.species[p2];
					for(let cl2 = 0; cl2 < sp2.ncla; cl2++){
						let claa2 = sp2.cla[cl2];
				
						if(claa.name == claa2.name){
							if(claa.index != claa2.index) error("CLONE index not right");
							
							let cam = claa.camera;
							let cam2 = claa2.camera;
				
							if(cam.coord != cam2.coord) error("CLONE cam coord not right");
							if(cam.ruler != cam2.ruler) error("CLONE cam ruler not right");
							if(cam.scale != cam2.scale) error("CLONE cam scale not right");
							if(cam.x != cam2.x) error("CLONE cam x not right");
							if(cam.y != cam2.y) error("CLONE cam y not right");
							
							if(claa.ncomp != claa2.ncomp) error("CLONE ncomp not right");
							else{
								for(let c = 0; c < claa.ncomp; c++){
									let co = claa.comp[c];
									let co2 = claa2.comp[c];
						
									if(cam.coord == "latlng" && co.type == "boundary"){}
									else{
										if(co.x != co2.x || co.y != co2.y) error("CLONE x,y not right");
									}
									
									if(co.col != co2.col) error("CLONE col not right");
						
									if(co.fixed.check != co2.fixed.check) error("CLONE fixed not right");
								}
							}
						}
					}
				}
			}
		}
	}
}

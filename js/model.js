"use strict";
// Function that define and control the model

class Model
{
	start;                                           // Set to true when a new model is started 
	
	description;                                     // A description of the analysis
	filename;                                        // The current filename 
	species = [];                                    // Stores information about species
	param = [];                                      // Stores parameter information
	param_factor = [];                               // Stores factors which multiply parameter (for ppc)
	spline = [];                                     // Stores spline information
	warn = [];                                       // Stores model warnings
	derive = [];                                     // Sets derived quantities
	
	sim_details = {};                                // Details of simulation
	inf_details = {};                                // Details of inference
	ppc_details = {};                                // Details of posterior predictive check
	
	sim_res = { siminf:"sim"};                       // Filter used to plot sim results
	inf_res = { siminf:"inf"};                       // Filter used to plot inf results
	ppc_res = { siminf:"ppc"};                       // Filter used to plot post-sim results
	
	update_model_needed = false;                     // Determines if init parameter is needed
	warn_view = true;                                // Determines if warning are shown on main page
	
	constructor()
	{
		this.start = false;
		this.update_model_needed = true;
	}


	///	Starts a new model
	start_new(op)
	{
		this.description = {te:"# Model and analysis description."};
		this.example = undefined;
		this.species = [];
		this.param = [];
		this.param_factor = [];
		this.derive = [];   
		
		this.sim_details = { t_start:"", t_end:"", timestep:"", indmax:INDMAX_DEFAULT, param_output_max:PARAM_OUTPUT_MAX_DEFAULT, algorithm:{value:"gillespie"}, number:SIM_NUM_DEFAULT, run_local:{value:"Yes"}, run_save_type:{value:"Save"}, seed_on:{value:"No"}, seed:SEED_DEFAULT };    
		
		this.inf_details = { t_start:"", t_end:"", timestep:"", abcsample:String(ABC_SAMPLE_DEFAULT), sample:String(MCMC_SAMPLE_DEFAULT), output_param :String(MCMC_OP_PARAM_DEFAULT), output_state:String(MCMC_OP_STATE_DEFAULT), accfrac:String(ABC_ACFRAC_DEFAULT), accfracsmc:String(ABCSMC_ACFRAC_DEFAULT), numgen:String(ABCSMC_GEN_DEFAULT), kernelsize:String(ABCSMC_KERNEL_DEFAULT), indmax:INDMAX_DEFAULT, param_output_max:PARAM_OUTPUT_MAX_DEFAULT, nchain:String(MCMC_CHAIN_DEFAULT), algorithm:{value:ALG_DEFAULT}, run_local:{value:"Yes"}, run_save_type:{value:"Save"}, seed_on:{value:"No"}, seed:SEED_DEFAULT, burnin_frac:BURNIN_FRAC_DEFAULT, anneal_type:{te:ANNEAL_DEFAULT}, anneal_rate:ANNEAL_RATE_DEFAULT, anneal_power:ANNEAL_POWER_DEFAULT, npart:String(PAS_PART_DEFAULT), gen_update:String(PAS_GEN_UPDATE_DEFAULT), cha_per_core:String(MCMC_CHAIN_PER_CORE_DEFAULT), part_per_core:String(PAS_PART_PER_CORE_DEFAULT) };
		
		this.ppc_details = {  ppc_t_start:"", ppc_t_end:"", t_start:"", t_end:"", algorithm:{value:"gillespie"}, number:PPC_NUM_DEFAULT, run_local:{value:"Yes"}, run_save_type:{value:"Save"}, seed_on:{value:"No"}, seed:SEED_DEFAULT};    
			
		this.start = true;
		this.filename = "";
		
		this.sim_res.run_warning=[];
		this.sim_res.run_warning_show = false; 
		this.inf_res.run_warning=[];
		this.inf_res.run_warning_show = false; 
		this.ppc_res.run_warning=[];
		this.ppc_res.run_warning_show = false; 
	}
	
	
	/// Gets the species number currently being viewed
	get_p(op)
	{
		let pag = inter.page[inter.pa];
		let sub = pag.sub[pag.index];
		
		if(sub.name == "Results"){
			if(!sub.sub) return 0;
			
			let subsub = sub.sub[sub.index];
		
			switch(subsub.name){
			case "Populations":
			case "Transitions":
			case "Individuals":
			case "Diagnostics":
				return subsub.index;
				
			default:
				if(op == "no error") return;
				error("Should not be getting p");
				return;
			}
		}
		else{
			if(!(pag.name =="Simulation" && sub.name == "Initial Conditions") &&
				!(pag.name =="Simulation" && sub.name == "Generate Data") &&
				!(pag.name =="Inference" && sub.name == "Initial Conditions") &&
				!(pag.name =="Inference" && sub.name == "Data") &&
				!(pag.name =="Post. Simulation" && sub.name == "Population Mod.") &&
				!(pag.name =="Model" && sub.name == "Compartments")){
				if(op == "no error") return;
				error("Should not be getting p");
			}
			
			return sub.index;
		}
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


	/// Gets the classification
	get_cla()
	{
		let p = this.get_p(); if(p == undefined) return;
		let cl = this.get_cl(); if(cl == undefined) return;
		if(p >= this.species.length) return;
		if(cl >= this.species[p].cla.length) return;
		return this.species[p].cla[cl];
	}


	///	Returns true if the compartmental model is being viewed 
	get_show_model()
	{
		if(this.show_warning() == true) return false;
		
		if(tab_name() == "Model" && subtab_name() == "Compartments"){
			return true;
		}
		return false;
	}
	
	
	/// Updates a sbox region
	update_sbox(sbox,x,y,dx,dy,mar)
	{
		if(mar < 1.5) mar = 1.5;
		let x1 = x-mar, y1 = y-mar, x2 = x+dx+mar, y2 = y+dy+mar; 
		if(!sbox.x1){
			sbox = {on:true, x1:x1, y1:y1, x2:x2, y2:y2};
		}
		else{
			if(x1 < sbox.x1) sbox.x1 = x1;
			if(x2 > sbox.x2) sbox.x2 = x2;
			if(y1 < sbox.y1) sbox.y1 = y1;
			if(y2 > sbox.y2) sbox.y2 = y2;
		}
	
		return sbox;
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
		
		let cam = claa.camera;                         // Loads up the camera used to view model
		
		let bub = inter.bubble;

		let sbox = {};
		
		for(let k = 0; k < claa.ncomp; k++){
			let sel = false; if(find_in(inter.comp_select.list,k) != undefined) sel = true;
			
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
					let w =	c.w*cam.scale*Math.exp(cam.slider.value);
					let h =	c.h*cam.scale*Math.exp(cam.slider.value); 
					let x = pt.x-w/2, y = pt.y-h/2;
					if(sel){
						sbox = this.update_sbox(sbox,x,y,w,h,h/2);			
					}
					
					lay.add_button({te:claa.comp[k].name, x:x, y:y, dx:w, dy:h, ac:ac, type:"Compartment", col:c.col, col_dark:dark_colour(c.col), p:p, cl:cl, i:k});
				}
				break;
			
			case "latlng":
				{
					let si = cam.scale*cam.ruler*Math.exp(cam.slider.value);
					if(si > si_limit_circle){
						let w =	2*latlng_radius*si, h = w; 
						let pt = trans_point(c.x,c.y,cam,lay);
						let x = pt.x-w/2, y = pt.y-h/2;
						if(sel){
							sbox = this.update_sbox(sbox,x,y,w,h,h/2);			
						}
						lay.add_button({te:claa.comp[k].name, x:x, y:y, dx:w, dy:h, ac:ac, type:"CompLatLng", col:c.col, col_dark:dark_colour(c.col), p:p, cl:cl, i:k});
					}
				}
				break;
				
			case "boundary":
				{
					let ms = find_map_store(c.map_ref);
					let box = ms.feature.box;
				
					let pmin = trans_point(box.xmin,box.ymax,cam,lay);
					let pmax = trans_point(box.xmax,box.ymin,cam,lay);
					
					if(sel){
						sbox = this.update_sbox(sbox,pmin.x,pmin.y,pmax.x-pmin.x,pmax.y-pmin.y,0);
					}
						
					lay.add_button({x:pmin.x, y:pmin.y, dx:pmax.x-pmin.x, dy:pmax.y-pmin.y, ac:ac, type:"CompMap", polygon:ms.feature.polygon, mask:ms.mask, col:c.col, p:p, cl:cl, i:k});
				}
				break;

			default: error("Option not recognised 82"); break;
			}
		}
		
		if(mo.type == "Add_Compartment" || mo.type == "Add_Label"){
			lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, ac:"NothingMode", type:"Nothing"});	
		}
		
		if(sbox.on){
			inter.comp_select.sbox = sbox;
		}
	}


	/// Adds a button for selected compartments
	add_comp_select_buts(lay)
	{
		let sbox = inter.comp_select.sbox;
		let x1 = sbox.x1, y1 = sbox.y1, x2 = sbox.x2, y2 = sbox.y2;
		lay.add_button({x:x1, y:y1, dx:x2-x1, dy:y2-y1, ac:"Selected", type:"Selected", p:this.get_p(), cl:this.get_cl()});	
	}
	
	
	/// Adds a button for selected compartments
	add_comp_select_button_buts(lay)
	{
		let cs = inter.comp_select;
		
		let sbox = cs.sbox;
		let x1 = sbox.x1, y1 = sbox.y1, x2 = sbox.x2, y2 = sbox.y2;

		let ddx = 1.1, ddy = ddx*231/218, mar = 0.4;
			
		let claa = this.species[cs.p].cla[cs.cl];
		let no_copy = false;
		for(let i = 0; i < cs.list.length; i++){
			let co = claa.comp[cs.list[i]];
			if(co.type == "boundary") no_copy = true;
		}
		
		if(no_copy == false){
			lay.add_button({x:x2-ddx-mar, y:y2-ddy-mar, dx:ddx, dy:ddy, ac:"CopyPic", type:"CopyPic"});	
		}
		
		ddx = 1.1, ddy = 1.1;
		lay.add_button({x:x2-ddx-mar, y:y1+mar, dx:ddx, dy:ddy, ac:"DeleteSel", type:"Delete"});	
	}


	/// Copies the selected compartments 
	copy_selected(suf)
	{
		let p = inter.comp_select.p, cl = inter.comp_select.cl;
		let clone = model.find_clones(p,cl);	
		clone.push({p:p, cl:cl});
		
		let list = inter.comp_select.list;
	
		for(let loop = 0; loop < clone.length; loop++){
			let clo = clone[loop];
			let p = clo.p, cl = clo.cl;
		
			let claa = this.species[p].cla[cl];	
		
			let list_new = [];
			for(let i = 0; i < list.length; i++){
				let c = list[i];
				list_new.push(claa.ncomp);
				let co = copy(claa.comp[c]);
				co.fixed = {check: false}
				co.name += suf;
				switch(co.type){
				case "box": 
					this.set_compartment_size(co); 
					co.x += 2; co.y += 2; 
					break;
				case "latlng":
					{
						let d = 2*claa.camera.ruler;
						co.x += d; co.y += d; 
					}
					break;
				}
				
				hash_add(claa.hash_comp,co.name,claa.ncomp);
			
				claa.comp.push(co);
				claa.ncomp++;
			}
			
			let map_c = [];
			for(let i = 0; i < list.length; i++){
				map_c[list[i]] = list_new[i];
			}
		
			// Copies any transitions
			for(let i = 0; i < claa.ntra; i++){
				let tr = claa.tra[i];
				
				let i_in; if(tr.i != UNSET && map_c[tr.i] != undefined) i_in = map_c[tr.i];
				let f_in; if(tr.f != UNSET && map_c[tr.f] != undefined) f_in = map_c[tr.f];
				
				if(i_in != undefined || f_in != undefined){
					let tr_new = copy(tr);
					if(i_in != undefined) tr_new.i = i_in;
					if(f_in != undefined) tr_new.f = f_in;
					
					claa.tra.push(tr_new);
					this.set_transition_name(claa,claa.ntra);
		
					if(tr_new.i == SOURCE) claa.nsource++;
					if(tr_new.f == SINK) claa.nsink++;
					claa.ntra++;
				}
			}
	
			this.update_pline(claa);
			inter.comp_select.list = list_new;
		}
	}
	
	
	/// Deletes selected compartments
	delete_selected()
	{
		let p = inter.comp_select.p, cl = inter.comp_select.cl;
		let clone = model.find_clones(p,cl);	
		clone.push({p:p, cl:cl});
			
		let list = inter.comp_select.list;
	
		let map_c = [];
		for(let i = 0; i < list.length; i++) map_c[list[i]] = true;
	
		for(let loop = 0; loop < clone.length; loop++){
			let clo = clone[loop];
			let p = clo.p, cl = clo.cl;
		
			let claa = this.species[p].cla[cl];	
		
			let i = 0;
			while(i < claa.ntra){
				let tr = claa.tra[i];
				if((tr.i != UNSET && map_c[tr.i] == true) ||
					 (tr.f != UNSET && map_c[tr.f] == true)){
				 	if(tr.i == SOURCE) claa.nsource--;
					if(tr.f == SINK) claa.nsink--;
					
					claa.tra.splice(i,1);
					claa.ntra--;
				}					 
				else i++;
			}
			hash_redo(claa.hash_tra,claa.tra);

			for(let c = claa.ncomp-1; c >= 0; c--){
				if(map_c[c] == true){
					claa.comp.splice(c,1);
					claa.ncomp--;
				}
			}	
			hash_redo(claa.hash_comp,claa.comp);
			
			this.update_pline(claa);
		}
		
		clear_comp_select();
	}
	
	
	/// New added compartment
	add_compartment_new_buts(lay) 
	{
		let p = this.get_p();
		let cl = this.get_cl();
		
		let cam = this.species[p].cla[cl].camera;
	
		switch(cam.coord){
		case "cartesian": lay.add_button({te:"", x:0, y:0, dx:lay.dx-marnew, dy:lay.dy-marnew, type:"Compartment", col:WHITE, col_dark:dark_colour(WHITE)}); break;
		case "latlng": lay.add_button({te:"", x:0, y:0, dx:lay.dx-marnew, dy:lay.dy-marnew, type:"CompLatLng", col:BLACK, col_dark:DDGREY}); break;
		default: error("SHoudl not be def"); break;
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
					let res = this.add_transition(p,cl,mo.i,f,midp,"");
					output_help(res);
					
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
				let res = this.add_transition(p,cl,SOURCE,f,midp,"");
				output_help(res);
				
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
				let res = this.add_transition(p,cl,i,SINK,midp,"");
				output_help(res);
				
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
				points[i] = {x:lay.dx/2 + (pl.x-cam.x)*cam.scale, y:lay.dy/2 - (pl.y-cam.y)*cam.scale, index:pl.index};
			}
			
			let center;
			if(tr.center != undefined){
				let cen = tr.center;
				center = {x:lay.dx/2 + (cen.x-cam.x)*cam.scale, y:lay.dy/2 - (cen.y-cam.y)*cam.scale, nx:cen.nx, ny:-cen.ny, w:cen.w};
			}
			
			let d = TRANS_OVER_RANGE;
		
			let ac = "Transition";
			if(mode_type == "Add_Source" || mode_type == "Add_Sink" || mode_type == "No Click") ac = undefined;
			
			lay.add_button({x:lay.dx/2 + (tr.box.x-cam.x)*cam.scale - d, y:lay.dy/2 - ((tr.box.y+tr.box.dy)-cam.y)*cam.scale-d, dx:tr.box.dx*cam.scale+2*d, dy:tr.box.dy*cam.scale+2*d, ac:ac, type:"Transition", cam:cam, points:points, center:center, tr:tr, p:p, cl:cl, i:k});	
				
			for(let j = 0; j < points.length; j++){
				let po = points[j];
				if(po.index >= 0 && po.index < tr.midp.length){
					let r = TRANS_POINT_R; 
					let variety = MIDPOINT;
					let ac = "TransitionPoint";
					if(tr.i == SOURCE && j == 0){ 
						r = 0.7*cam.scale*cam.ruler;
						variety = SOURCE; ac = "Transition";
					}
					if(tr.f == SINK && j == points.length-1){ 
						r = 0.7*cam.scale*cam.ruler;
						variety = SINK; ac = "Transition";
					}

					if((variety == SOURCE || variety == SINK) && cam.coord == "latlng") r *= 0.6;
					
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
							if(variety == SOURCE || variety == SINK) ac = undefined;
						}
					}

					if(ac != undefined || variety == SOURCE	|| variety == SINK){
						lay.add_button({x:po.x - r, y:po.y - r, dx:2*r, dy:2*r, ac:ac, type:"TransitionPoint", tr:tr, p:p, cl:cl, i:k, index:po.index, variety:variety, center:center, tr:tr, p:p, cl:cl, i:k});
					}
				}
			}
		}
	}

	
	/// Determines if a source exists in the model
	source_exist(p)
	{
		let sp = this.species[p];
		
		for(let cl = 0; cl < sp.cla.length; cl++){
			let claa = sp.cla[cl];
			
			for(let tr = 0; tr < claa.tra.length; tr++){
				let tra = claa.tra[tr];
				if(claa.tra[tr].i == SOURCE) return true;
			}
		}
	
		return false;
	}
	
	
	/// Adds a label button
	add_addlabel_buts(lay)
	{
		let mo = inter.mode;
		mo.l = lay.index;
		lay.add_button({te:mo.te, x:0, y:0, dx:lay.dx, dy:lay.dy, si:mo.si, font:mo.fo, type:"Text", col:mo.col});
	}
	

	/// Adds a label annotation to the model
	add_label(te,size,p,cl)
	{
		let mo = inter.mode;
		let claa = this.species[p].cla[cl];
		claa.annotation.push({type:"text", te:te, size:size, x:mo.mp.x, y:mo.mp.y, color:mo.col});
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
		let yc = -((inter.my - lay_main.y - lay_main.dy/2)/cam.scale)+cam.y;
			
		if(cam.coord == "cartesian"){
			return {x:Math.round(xc/grid)*grid, y:Math.round(yc/grid)*grid};
		}
	
		return {x:xc, y:yc};	
	}


	/// Gets the scale 
	get_scale()
	{
		return this.species[this.get_p()].cla[this.get_cl()].camera.scale;
	}
		
		
	/// Gets the camera 
	get_cam()
	{
		return this.species[this.get_p()].cla[this.get_cl()].camera;
	}
	
	
	/// Gets the coordinate system
	get_coord()
	{
		return this.species[this.get_p()].cla[this.get_cl()].camera.coord;
	}
	
	
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
		let y = lay_main.dy/2 - (mp.y-cam.y)*cam.scale;

		return {x:x, y:y};
	}


	/// In the "add compartment" mode this adds the new compartment
	place_add_compartment(type)
	{
		let lay_comp;
		switch(type){
		case "Add_Compartment": lay_comp = get_lay("AddCompartment"); break;
		case "Add_Label": lay_comp = get_lay("AddLabel"); break;
		default: error("Add op problem"); break;
		}
		
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
		
				if(lay_comp.x < lay_main.x || lay_comp.x+lay_comp.dx > lay_main.x+lay_main.dx ||
					 lay_comp.y < lay_main.y || lay_comp.y+lay_comp.dy > lay_main.y+lay_main.dy) lay_comp.x = 1000;
			
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
	add_species(name,type,trans_tree)
	{
		let p = find(this.species,"name",name);
	
		if(p != undefined){
			return err("There is already a species with the name '"+name+"'");
		}
		
		for(let p = 0; p < this.species.length; p++){
			let sp = this.species[p];
			
			let i = find(sp.cla,"name",name);
			if(i != undefined){
				return err("Cannot have the same name as classification '"+sp.cla[i].name+"' ");
			}
			
			let j = find(sp.cla,"index",name);
			if(j != undefined){
				return err("Cannot have the same name as index for classification '"+sp.cla[j].name+"' ");
			}
		
			for(let cl = 0; cl < sp.ncla; cl++){
				let c = find(sp.cla[cl].comp,"name",name);
				if(c != undefined){
					return err("Cannot have the same name as compartment '"+sp.cla[cl].comp[i].name+"' in classification '"+sp.cla[cl].name+"'");
				}
			}
		}
		
		this.species.push({ name:name, sim_source:[], inf_source:[], data_type:{te:data_types[0]
		}, cla:[], ncla:0, fix_eff:[], ind_eff_group:[], type:type, trans_tree:{check:trans_tree}, infection_cl:{te:select_str}});		
		//, sim_init_cond:{ type:{value:"Fixed"}, focal:{value:"Focal"}, focal_cl:0}}
		
		if(debug == true) this.check_consistent();  
		
		update_param();
		
		return success();
	}
	
	
	/// Adds a new classification to the model
	add_classification(p,name,index,op)
	{
		if(p == undefined){
			return err("Cannot add a classification when the species is not set");
		}
	
		let cl = find(this.species[p].cla,"name",name);
		if(cl != undefined){
			return err("There is already a classification with the name '"+name+"'");
		}
		
		if(name == index){
			return err("Name and index must be different");
		}
		
		let coord = op.coord; if(coord == undefined) coord = COORD_DEFAULT;
	
		let slider = {update:"scale",lay:"LowerMenu", min:-Math.log(10),max:Math.log(10), value:0};
		
		let cam = {x:0, y:0, scale:1, ruler:1, coord:coord, grid:"off", set:false, slider:slider};
		
		let dmap = false; if(coord == "latlng" && op.default_map == true) dmap = "loading"; 
		if(coord == "latlng") cam.ruler = default_ruler;
		
		this.species[p].cla.push({ name:name, index:index, tra:[], ntra:0, hash_tra:[], comp:[], ncomp:0, hash_comp:[], nsource:0, nsink:0, annotation:[], default_map:dmap, camera:cam});
		this.species[p].ncla++;
		
		data_update_add_classification(p);
		
		if(debug == true) this.check_consistent();  
		
		return success();
	}
	

	/// Clones a classification from another species
	clone_class(p_to,p_from,cl_from)
	{
		let cla_from = this.species[p_from].cla[cl_from];
	
		this.add_classification(p_to,cla_from.name,cla_from.index,{coord:cla_from.camera.coord});
		
		let cl = this.species[p_to].ncla-1;
		
		this.species[p_to].cla[cl] = copy(cla_from);
		
		let cla_to = this.species[p_to].cla[cl];
		cla_to.tra=[]; cla_to.ntra = 0; cla_to.tra_hash=[];
		cla_to.nsource = 0; cla_to.nsink = 0;
	}


	/// Adds a new compartment to the model
	add_compartment(name,p,cl,x,y,col,fix,infected,include_clone)
	{
		if(fix == undefined) fix = false;
		
		let res = this.check_comp_exist(name,p); if(res.err == true) return res;
	
		let sp = this.species[p];
		let claa = sp.cla[cl];
		let comp = claa.comp;
		
		if(infected == true) sp.infection_cl.te = claa.name;
		
		hash_add(claa.hash_comp,name,claa.ncomp);
		
		comp.push({name:name, type:"box", x:x, y:y, col:col, fixed:{check:fix}, infected:{check:infected}, markov_branch:false, all_branches:false});
		claa.ncomp++;
			
		let i = comp.length-1;
		
		this.set_compartment_size(comp[i]);
		
		this.snap_comp_to_grid(comp[i]);
		this.determine_branching();
		
		if(include_clone == true){                     // Adds to clones 
			let list = this.find_clones(p,cl);
		
			for(let li of list){
				let res = this.add_compartment(name,li.p,li.cl,x,y,col,fix,infected,false);
				if(res.err == true) return res;
			}
		}
		
		if(debug == true) this.check_consistent();  
		
		update_param();
		
		return success();
	}
	

	/// Adds a new compartment to the model
	add_latlng_compartment(name,p,cl,x,y,col,fix,infected)
	{
		if(fix == undefined) fix = false;
	
		let res = this.check_comp_exist(name,p); if(res.err == true) return res;
		
		let sp = this.species[p];
		let claa = sp.cla[cl];
		let comp = claa.comp;
		
		if(infected == true) sp.infection_cl.te = claa.name;
		
		hash_add(claa.hash_comp,name,claa.ncomp);
		
		comp.push({name:name, type:"latlng", x:x, y:y, col:col, fixed:{check:fix}, infected:{check:infected}, markov_branch:false, all_branches:false});
		claa.ncomp++;
		
		let i = comp.length-1;
		
		this.determine_branching();

		if(debug == true) this.check_consistent();  
		
		update_param();
		
		return success();
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
	add_map_compartment(name,p,cl,feature,mask,col,infected)
	{
		let res = this.check_comp_exist(name,p); if(res.err == true) return res;
		
		// Calculates the midpoint of the region as the middle of the mask
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
		
		let map_name = "file"+Math.random();
		map_store.push({name:map_name, feature:feature, mask:mask});
		
		hash_add(claa.hash_comp,name,claa.ncomp);
		
		comp.push({name:name, type:"boundary", xmid:xmid, ymid:ymid, map_ref:map_name, col:col, fixed:{check:true}, infected:{check:infected}, markov_branch:false, all_branches:false});
		claa.ncomp++;
		this.determine_branching();
		
		update_param();
		
		return success();
	}

	
	/// Checks if a compartment already exists 
	check_comp_exist(name,p)
	{
		if(name != ""){
			let sp = this.species[p];
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				if(claa.name == name){
					return err("Cannot have a classification and compartment both '"+name+"'");
				}					
				
				if(claa.index == name){
					return err("Cannot have a classification index and compartment both '"+name+"'");
				}		
				
				let c = hash_find(claa.hash_comp,name);
				if(c != undefined){
					return err("There is already a compartment with the name '"+name+"'");
				}
			}
		}
		return success();
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
	add_transition(p,cl,i,f,midp,ty,op)
	{
		let sp = this.species[p];
		let claa = sp.cla[cl];
		let tra = claa.tra;
	
		if(ty != ""){
			if(i == SOURCE){
				for(let cl2 = 0; cl2 < cl; cl2++){
					if(cl2 != cl && sp.cla[cl2].nsource > 0){
						return err("Cannot have sources in multiple classifications (here sources are in '"+sp.cla[cl].name+"' and '"+sp.cla[cl2].name+"')");
					}
				}
			}

			if(f == SINK){
				for(let cl2 = 0; cl2 < cl; cl2++){
					if(cl2 != cl && sp.cla[cl2].nsink > 0){
						return err("Cannot have sinks in multiple classifications (here sinks are in '"+sp.cla[cl].name+"' and '"+sp.cla[cl2].name+"')");
					}
				}
			}

			let exist_fl = false;
			
			let name = this.get_trans_name(claa,i,f);
			
			if(hash_find(claa.hash_tra,name) != undefined){
				if(i == SOURCE){
					return err("Cannot have two sources into the same compartment");
				}
				else{
					if(f == SINK){
						return err("Cannot have two sinks leaving the same compartment");
					}
					else{
						return err("Cannot have two transitions with the same source and destination compartments");
					}
				}
			}
	
			if(i == f){
				return err("The 'from' and 'to' compartments must be different");
			}
		}
		
		let value = {};
		if(ty != ""){
			if(i == SOURCE){
				let sup2 = claa.comp[f].name;
				value = this.new_trans_equations(sup2,sup2,true,p,cl);
			}
			else{	
				if(f == SINK){
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
		
		tra.push({type:ty, i:i, f:f, pline:[], branch_select:undefined, midp:midp, value:value});	
		this.set_transition_name(claa,claa.ntra);
		
		if(ty != ""){
			if(i == SOURCE) claa.nsource++;
			if(f == SINK) claa.nsink++;	
		}
		
		claa.ntra++;
		
		if(op != "simple"){
			if(ty != "") this.determine_branching();
		
			this.update_pline(claa);  

			if(debug == true) this.check_consistent();  
		}
		
		update_param();
	
		return success();
	}
	
	
	/// Updates the model
	update_check(p,cl)
	{
		this.determine_branching();
		let claa = this.species[p].cla[cl];		
		this.update_pline(claa);  
		if(debug == true) this.check_consistent();  
	}
	
	
	/// Creates a new set of transition equations
	new_trans_equations(sup,sup2,flag,p,cl)
	{
		let value = {};
		if(true){                                     // Turned off default transition nales
			if(flag == false){
				value.mean_eqn = create_equation("","trans_mean",p,cl);
				value.rate_eqn = create_equation("","trans_rate",p,cl);
			}
			else{
				value.mean_eqn = create_equation("","trans_mean",p,cl);
				value.rate_eqn = create_equation("","trans_rate",p,cl);
				value.bp_eqn = create_equation("","trans_bp",p,cl);
				value.shape_eqn = create_equation("","trans_shape",p,cl);
				value.scale_eqn = create_equation("","trans_scale",p,cl);
				value.cv_eqn = create_equation("","trans_cv",p,cl);
				value.shape_erlang = {te:"2"};
			}
		}
		else{
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
			if(par.length == 0)	return err("Does not contain a parameter");
			if(par.length > 1) return err("Should not contain more than one parameter");
		}
		else{
			if(par[0].name == "D") return err("Name 'D' is reserved for the distance matrix");
			if(par[0].name == "t") return err("Name 't' is reserved for time");
		}
		
		let res = check_derived_param(ob.eqn2,ob.eqn1);
		if(res.err == true) return res;
		
		let name = par[0].name;

		let i = find(this.param,"name",name);

		if(i != undefined && mod != "modify"){
			return err("The parameter '"+name+"' is already used in the model"); 
		}
	
		if(mod == "modify") this.derive[ob.val] = ob;
		else this.derive.push(ob);
		
		update_param();
	
		return success();
	}
	
	
	/// Determines if compartments/transition are branching
	determine_branching(op)
	{
		for(let p = 0; p < this.species.length; p++){
			let sp = this.species[p];
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				for(let tr = 0; tr < claa.ntra; tr++){
					claa.tra[tr].branch = false;
					claa.tra[tr].all_branches = false;
				}

				for(let c = 0; c < claa.ncomp; c++){
					let co = claa.comp[c];
					
					let list=[]; 
					let flag = false;
					let all_set = true;
					for(let tr = 0; tr < claa.ntra; tr++){
						let trr = claa.tra[tr];
						if(trr.i == c){
							list.push(tr);
										
							if(trr.type != "exp(rate)" && trr.type != "exp(mean)" && tr.type != "erlang") flag = true;
							
							if(trr.branch_select != true) all_set = false;
						}
					}
					
					
					if(op == "set all_branches"){           // Sets the all_branches flag on compartments
						if(list.length > 1 && all_set == true){
							co.all_branches = true;
						}
					}
					
					if(list.length > 1 && all_set == true){
						claa.tra[list[list.length-1]].branch_select = false;
					}
					
					for(let i = 0; i < list.length; i++){
						let tra = claa.tra[list[i]];
						tra.all_branches = co.all_branches;
					}
					
					if(list.length > 1 && flag == false) co.choose_branch = true;
					else co.choose_branch = false;
				
					if(list.length > 1 && (flag == true || co.markov_branch == true)){
						co.branch = true;
					}
					else co.branch = false;
					
					for(let i = 0; i < list.length; i++){
						let tra = claa.tra[list[i]];
						tra.branch = co.branch;
					
						if(list.length > 1 && tra.branch_select == undefined){
							if(i == 0) tra.branch_select = false;
							else tra.branch_select = true;
						}
					}
				}					
			}
		}
	}


	/// Deletes the branching probability (such that it is defined for another transition)
	delete_branch_prob(p,cl,i)
	{
		let claa = this.species[p].cla[cl];
		let tra_sel = claa.tra[i];
		
		for(let tr = 0; tr < claa.ntra; tr++){
			let trr = claa.tra[tr];
			if(trr.i == tra_sel.i) trr.branch_select = true;
		}
		
		tra_sel.branch_select = false;
		
		this.determine_branching();
		this.update_pline(claa);
	}
	
	get_trans_name(claa,i,f)
	{
		if(i == SOURCE) return "+→"+claa.comp[f].name;
		else{
			if(f == SINK) return claa.comp[i].name+"→-";
			else return claa.comp[i].name+"→"+claa.comp[f].name;
		}	
	}
	
	/// Sets the name of a transition
	set_transition_name(claa,tr)
	{
		let tra = claa.tra[tr];
		
		if(tra.type != ""){
			let name = this.get_trans_name(claa,tra.i,tra.f);
			
			let variety = NORMAL;
			if(tra.i == SOURCE) variety = SOURCE;
			else{
				if(tra.f == SINK) variety = SINK; 
			}
	
			if(name != tra.name){
				if(tra.name != undefined && tra.name != "") hash_remove(claa.hash_tra,tra.name);
				hash_add(claa.hash_tra,name,tr);
				tra.name = name;
			}
			
			tra.variety = variety; 
		}
	}


	/// Sets the size of a compartment (depending on the height of the text which defined it)
	set_compartment_size(comp)
	{
		comp.h = compartment_height;

		cv.font = "30px Times";
		let name = comp.name;
		let w = 2;
		if(name.includes("^")){
			let spl = name.split("^");
			w += comp.h*si_comp_text_frac*(cv.measureText(spl[0]).width+cv.measureText(spl[1]).width*si_sup_fac)/30;
		}
		else{
			w += comp.h*si_comp_text_frac*(cv.measureText(name).width/30);
		}
		
		if(w < compartment_width_min) w = compartment_width_min;
		
		w = (Math.floor(w/grid))*grid;
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
			lay.add_button({te:message, x:0, y:0, dx:lay.dx, dy:lay.dy, type:"MessageMenu"});
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
		
		let gap = 3.2;
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
	
		x += 2;
		if(cl < n) active = true; else active = false;
		w = this.add_object_button(lay,"Annotation",x,y,"AddAnnotation",{active:active, p:p, cl:cl, title:"Add an annotation", te:annotation_text, back:WHITE});
		x += w+gap;
		
		w = this.add_object_button(lay,"Import",x,y,"ImportModel",{active:active, p:p, cl:cl, title:"Import", te:import_text, back:WHITE});
	
		let ac_in, ac_out;
	
		if(cl < n){ 
			ac_in = "ZoomIn"; ac_out = "ZoomOut";
		}
		
		lay.add_button({x:lay.dx-5, y:0, dx:1.8, dy:2, ac:ac_in, type:"ZoomIn", p:p, cl:cl});

		lay.add_button({x:lay.dx-3, y:0, dx:1.8, dy:1.8, ac:ac_out, type:"ZoomOut", p:p, cl:cl});
		
		let grid_ac; if(cl < n) grid_ac = "GridIcon";
		lay.add_button({x:lay.dx-7, y:0.4, dx:1.5, dy:1.5, ac:grid_ac, type:"GridIcon", p:p, cl:cl});
		
		if(claa){
			let cam = claa.camera;
			lay.add_button({x:lay.dx-12.5, y:0.7, dx:5, dy:1, p:p, cl:cl, info:cam.slider, ac:"Slider", type:"Slider"});
		}
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
		if(loading_sym()) return;
		
		let lay = get_lay("LowerMenu"); if(lay == undefined) return;
		
		let i = lay.search_button_ac("AddCompartment");

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
		
		
	/// Adds button to view warning
	add_view_warning_buts(lay)
	{		
		lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, ac:"WarningStore", type:"WarningStore"});
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
	snap_comp_to_grid(ob)
	{
		if(ob.x == "auto") return;
		
		ob.x = Math.round(ob.x/grid)*grid;
		ob.y = Math.round(ob.y/grid)*grid;
	}


	/// Deletes a transition  point
	delete_transitionpoint(p,cl,i,index)
	{
		let claa = this.species[p].cla[cl];
		claa.tra[i].midp.splice(index,1);
		this.update_pline(claa); 
	}
	
	
	/// Creates a new transition point from dragging transition
	split_transitionpoint()
	{
		let mo = inter.mode;
		let claa = this.species[mo.p].cla[mo.cl];
		claa.tra[mo.i].midp.splice(mo.index+1,0,this.get_mouse_to_desktop());
		this.update_pline(claa);
		
		inter.mode={ type:"Drag_TransitionPoint", p:mo.p, cl:mo.cl, i: mo.i, index: mo.index+1};
	}


	/// Updates all the plines
	update_pline_all()
	{
		for(let p = 0; p < this.species.length; p++){
			let sp = this.species[p];
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				this.update_pline(claa);
			}
		}
	}
	
	
	/// Updates all plines in a classification
	update_pline(claa)
	{
		for(let i = 0; i < claa.ntra; i++){
			this.find_trans_pline(claa,i);
		}
	}
	
	
	/// Converts points in a transition to points on the model
	find_trans_pline(claa,i)                                 
	{
		let cam = claa.camera;
		let r_latlng = latlng_radius*cam.ruler;
		let slid = Math.exp(cam.slider.value);
		
		let tr = claa.tra[i];

		let midp = tr.midp;
		
		tr.pair = false;                              // Finds if a transition is pair with its reverse
		for(let j = 0; j < claa.ntra; j++){
			if(j != i){
				let tr2 = claa.tra[j];
				if(tr.i == tr2.f && tr.f == tr2.i){
					if(midp.length == 0 && midp.length == 0){
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
			
			pair_shift.x = -shift_gap*dy/r;
			pair_shift.y = shift_gap*dx/r;
		}			
		
		if(tr.i != undefined && tr.i != SOURCE){
			let co = claa.comp[tr.i];
			let p1 = this.comp_center(co);
		
			p1.x += pair_shift.x; p1.y += pair_shift.y;
			
			let	p2;
			if(midp.length > 0) p2 = midp[0];
			else{
				if(tr.f != SINK && tr.f != undefined) p2 = this.comp_center(claa.comp[tr.f]);
			}
			
			if(p2 != undefined){
				let pos = this.intersect_comp_line(co,p1,p2,r_latlng,slid);
				tr.pline.push({x:pos.x, y:pos.y, index:-1});
			}
		}
		
		for(let i = 0; i < midp.length; i++){
			tr.pline.push({x:midp[i].x, y:midp[i].y, index:i});
		}

		if(tr.f != undefined && tr.f != SINK){
			let co = claa.comp[tr.f];
			
			let p1 = this.comp_center(co);
		
			p1.x += pair_shift.x; p1.y += pair_shift.y;
			
			let p2;
			if(midp.length > 0) p2 = midp[midp.length-1];
			else{ if(tr.i != SOURCE) p2 = this.comp_center(claa.comp[tr.i]);}
		
			let pos = this.intersect_comp_line(co,p1,p2,r_latlng,slid);
			tr.pline.push({x:pos.x, y:pos.y, index:midp.length});
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
			let wmin = LARGE, wmax = -LARGE;            // Finds the potential width for the label
			for(let p = 0; p < np; p++){
				let dx = tr.pline[p].x - tr.center.x;
				let dy = tr.pline[p].y - tr.center.y;
				
				let d = dx*tr.center.nx + dy*tr.center.ny;
				if(d < wmin) wmin = d; if(d > wmax) wmax = d;
			}
			tr.center.w = wmax-wmin;
		}
		
		let xmax = -LARGE, xmin = LARGE;              // Finds a box around the transition points
		let ymax = -LARGE, ymin = LARGE;
			
		for(let j = 0; j < tr.pline.length; j++){
			let po = tr.pline[j];
			if(po.x < xmin) xmin = po.x; if(po.x > xmax) xmax = po.x;
			if(po.y < ymin) ymin = po.y; if(po.y > ymax) ymax = po.y;	
		}
		
		tr.box = {x:xmin, y:ymin, dx:xmax-xmin, dy:ymax-ymin};
		
		let val = tr.value;
		
		let lab="";
	
		if(tr.branch == true){
			if(tr.branch_select == true || tr.all_branches == true){
				let te = val.bp_eqn.te; if(te == "") te = "?";
				lab += te;
			}
			else lab += "*";
			
			lab += " | ";
		}
	
		switch(tr.type){
		case "exp(rate)": 
			lab += val.rate_eqn.te.trim(); 
			break;
			
		case "exp(mean)": 
			lab += val.mean_eqn.te.trim(); 
			break;
			
		case "gamma": 
			lab += "Γ("+val.mean_eqn.te.trim()+" , "+val.cv_eqn.te+")";
			break;
			
		case "erlang": 
			lab += "EL("+val.mean_eqn.te.trim()+" , "+val.shape_erlang.te.trim()+")";
			break;
			
		case "log-normal":
			lab += "LN("+val.mean_eqn.te.trim()+" , "+val.cv_eqn.te.trim()+")"; 
			break;
			
		case "weibull":
			lab += "W("+val.scale_eqn.te.trim()+" , "+val.shape_eqn.te.trim()+")";
			break;
			
		case "period": 
			lab += "P("+val.mean_eqn.te.trim()+")"; 
			break;
			
		default: delete tr.label; break;
		}
		
		this.set_transition_name(claa,i);
	
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
	intersect_comp_line(co,p1,p2,r_latlng,slid)
	{
		switch(co.type){
		case "box":
			{
				let lx=[], ly=[], jmax = 4;
					
				let dx = co.w*slid, dy = co.h*slid;
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
				r_latlng *= slid;
				
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

		let co = claa.comp[i];
		if(co == undefined) return;
		
		// First removes any transitions connected with compartment
		for(let	tr = claa.ntra-1; tr >= 0; tr--){
			let traa = claa.tra[tr];
			if(traa.i == i || traa.f == i) this.delete_transition(p,cl,tr);
		}

		// Removes compartment from any box annotations
		let name = co.name;
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

		// Removes any boundary data
		if(co.map_ref != undefined) delete_map_store(co.map_ref);
		
		claa.comp.splice(i,1);
		claa.ncomp--;

		hash_redo(claa.hash_comp,claa.comp);
	
		for(let	tr = 0; tr < claa.ntra; tr++){
			let traa = claa.tra[tr];
			if(traa.i > i) traa.i--;
			if(traa.f > i) traa.f--;
		}

		this.determine_branching();

		this.update_pline(claa);

		if(include_clone == true){  // Adds to clones 
			let list = this.find_clones(p,cl);
		
			for(let li of list){
				this.delete_compartment(li.p,li.cl,i,false);
			}
		}
		
		update_param();

		if(debug == true) this.check_consistent();  
	}


  /// Deletes a specified transition
	delete_transition(p,cl,i)                             
	{
		let claa = this.species[p].cla[cl];
		
		let tr = claa.tra[i];
		if(tr.i == SOURCE) claa.nsource--;
		if(tr.f == SINK) claa.nsink--;
					
		claa.tra.splice(i,1);
		claa.ntra--;
		
		hash_redo(claa.hash_tra,claa.tra);
		
		this.determine_branching();
		this.update_pline(claa); 

		update_param();
		
		if(debug == true) this.check_consistent();  
	}


	/// Checks that a compartment name is valid
	check_comp_name(te,op)
	{
		if(check_invalid_name(te)) return "Cannot use reserved word '"+te+"'";
	
		for(let i = 0; i < te.length; i++){
			let ch = te.substr(i,1);
		
			if(name_notallow.includes(ch)){
				if(ch == " ") return "Cannot contain any spaces";
				else{
					if(op == "full") return "Compartment name '"+te+"' cannot use the character '"+ch+"'";
					else return "Cannot use the character '"+ch+"'";
				}
			}
		}
	}
	

	/// Changes the name of a compartment
	rename_compartment(p,cl,c,new_name)
	{
		let sp = this.species[p];
		let p_name = sp.name;
		let claa = sp.cla[cl];
		let cl_name = claa.name;
		let old_name = claa.comp[c].name;
		
		clone_camera(p,cl);
	
		// This converts all classification on all species
		for(let p = 0; p < this.species.length; p++){
			let sp = this.species[p];
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				if(claa.name == cl_name){
					for(let c = 0; c < claa.ncomp; c++){
						let co = claa.comp[c];
						if(co.name == old_name){
							hash_remove(claa.hash_comp,old_name);
							hash_add(claa.hash_comp,new_name,c);
						
							co.name = new_name;
							this.set_compartment_size(co);
							this.update_pline(claa);
						}
					}
				}
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
		
		this.rename_ob(model,"model","comp_name",old_name,new_name);
		this.rename_ob(model,"model","comp_name_store",old_name,new_name);
		
		this.update_pline_all();

		if(debug == true) this.check_consistent();  
	}


	/// Rename a classification
	rename_classification(new_name,p,cl)
	{
		let sp = this.species[p];
		let p_name = sp.name;
		let claa = sp.cla[cl]
		
		let old_name = claa.name
		if(new_name == old_name) return;

		let eq_list = this.find_equation_list();
		
		let len = old_name.length;
		let dif = new_name.length - old_name.length;
		
		for(let i = 0; i < eq_list.length; i++){
			let eqn = eq_list[i];
			extract_equation_properties(eqn);
		
			if(eqn.cl_name == old_name) eqn.cl_name = new_name;
		}
	
		// This converts all classification on all species
		for(let p = 0; p < this.species.length; p++){
			let sp = this.species[p];
			for(let cl = 0; cl < sp.ncla; cl++){
				if(this.species[p].cla[cl].name == old_name){
					this.species[p].cla[cl].name = new_name;
				}
			}
		}

		this.rename_ob(model,"model","cl_name",old_name,new_name);
		this.rename_ob(model,"model","cl_name_store",old_name,new_name);
		
		if(sp.infection_cl.te == old_name) sp.infection_cl.te = new_name;
		
		this.update_pline_all();
		
		if(debug == true) this.check_consistent();  
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
	
		this.update_pline_all();
		
		if(debug == true) this.check_consistent();  
	}
	
	
	/// Changes a parameter property
	change_param(th,prop,val)
	{
		if(val != this.param[th][prop]){
			this.param[th][prop] = val;
			update_param();
		}
	}
	
	
	/// Renames an index on a classification
	rename_index(new_name,p,cl,include_clone)
	{
		let claa = this.species[p].cla[cl];
		
		let old_name = claa.index;
		if(new_name == old_name) return success();
		
		if(include_clone == true){ 
			for(let i = 0; i < this.param.length; i++){
				let par = this.param[i];
				for(let d = 0; d < par.dep.length; d++){
					let de = par.dep[d];
					let de_remove = remove_prime(de)
					if(de_remove == old_name){
						par.dep[d] = new_name+de.substr(de_remove.length);
						par.full_name = param_name(par);
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
				
				extract_equation_properties(eqn);
			}
	
			let list = this.find_clones(p,cl);
		
			for(let li of list){
				this.rename_index(new_name,li.p,li.cl,false);
			}
			
			claa.index = new_name;
	
			if(debug == true) this.check_consistent();  
		
			this.update_pline_all();
		}
		
		return success();
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
				switch(ele){
				case "eqn_appear": break;
				default:
					{
						let ob2 = ob[ele];
					
						if(typeof ob2 == 'object'){
							this.rename_ob(ob2,root+"->"+ele,prop,old_name,new_name);
						}
					}
					break;
				}
			}
		}
	}
	
	
	/// Checks if a string exists in an of the properties of an object (diagnostic)
	check_ob_string_exist(ob,root,name)
	{
		return;// TURN OFF IF NOT DEBUGGING
		
		if(root == "model") error("Checking for "+name);
		if(name == "") return;
		
		for(let ele in ob){
			if(ele != "eqn_appear"){ 
				let ob2 = ob[ele];
				if(ob2 == name){
					error("'"+name+"' found here: "+root);
				}
				
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
				
				if(so.type == "Compartment"){
					let tab = so.table;
					for(let r = 0; r < tab.nrow; r++){
						let te = tab.ele[r][2];
						//let spl = te.split("|");
						let spl = split_with_bracket(te,"|");
						
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
		
		// Goes through derived parameters
		for(let i = 0; i < this.derive.length; i++){
			let der = this.derive[i];
			
			let eqn_info = {i:i};
			
			this.add_equation_to_list(eqn_list,der.eqn1,eqn_info);		
		}
		
		return eqn_list;
	}
	
	
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
		reset_info_p();
		if(debug == true) this.check_consistent();  
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
	
		data_update_delete_classification(p,cl,name,comp_list);
		
		for(let i = 0; i < this.param.length; i++){
			let par = this.param[i];
			
			let j = find_in(par.dep,index);
			if(j != undefined){
				par.dep.splice(j,1);
				par_set_default(par);
			}
		}
		
		if(debug == true) this.check_consistent(); 
	}


	/// Checks that the model has been formed correctly
	check_consistent()
	{
		for(let p = 0; p < model.species.length; p++){
			let sp = model.species[p];
			if(sp.ncla != sp.cla.length) error(sp.ncla+" "+sp.cla.length+" ncla not consistent")
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				if(claa.ntra != claa.tra.length) error(claa.ntra+" "+claa.tra.length+" ntra ncla not consistent");
				if(claa.ncomp != claa.comp.length) error(claa.ncomp+" "+claa.comp.length+" ncomp ncla not consistent")
			}
		}
	}
	
	
	/// Determines if mouse is over a transition
	mouse_over_transition(points,lay)                         
	{		
		if(inter.mode.type == "Add_Transition") return;
		
		let mx = inter.mx-lay.x, my = inter.my - lay.y;
		let np = points.length;
		for(let p = 0; p < np-1; p++){
			let x = points[p].x, y = points[p].y;
			let xx = points[p+1].x, yy = points[p+1].y;
			
			let dx = xx-x, dy = yy-y;
			let dmx = mx-x, dmy = (my-y);
			let l = (dx*dmx + dy*dmy)/(dx*dx+dy*dy);
			if(l > 0 && l < 1){
				let para = (dx*dmy - dy*dmx)/Math.sqrt(dx*dx+dy*dy);
				if(para < TRANS_OVER_RANGE && para > -TRANS_OVER_RANGE) return points[p].index;
			}
		}
	}


	/// Adds warning buttons to the screen
	add_warning_buts(lay)
	{
		lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, type:"Rect", val:WHITE});
		
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
		percent(0);
		
		this.remove_repeated(tab,0);
		
		clear_classification(p,cl);
		
		let claa = this.species[p].cla[cl];	
		let cam = claa.camera;

		let pos_x=[], pos_y=[];
	
		if(tab.ncol <= 2){ // Makes up a set of points on a grid
			let si = compartment_height*si_comp_text_frac;
		
			let w_max = 0;
			for(let r = 0; r < tab.nrow; r++){
				let w = text_width_worker(tab.ele[r][0],si,false)+1;
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
				percent(100*r/tab.nrow);
				
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
					let c = hash_find(claa.hash_comp,name);
					if(c == undefined){
						let res = this.add_compartment(name,p,cl,x,y,col2,true,false,true);		
						output_help(res);
					}
					else{
						let co = claa.comp[c];
						co.x = x; co.y = y; co.col = col2;
						this.set_compartment_size(co);
					}
				}
				break;
			
			case "latlng":
				{
					let c = hash_find(claa.hash_comp,name);
					if(c == undefined){
						let res = this.add_latlng_compartment(name,p,cl,x,y,col2,true,false);
						output_help(res);
					}
					else{
						let co = claa.comp[c];
						co.x = x; co.y = y; co.col = col2;
					}
				}
				break;

			default: error("Option not recognised 89"); break;
			}
		}

		if(debug == true) this.check_consistent();  
	}
	
	
	/// Removes lines with repeated elements on 
	remove_repeated(tab,col)
	{
		let hash = new Hash();
		let r = 0;
		while(r < tab.nrow){
			let val = tab.ele[r][col];
			if(hash.find(val) != undefined){ tab.ele.splice(r,1); tab.nrow--;}
			else{ hash.add(val,r); r++;}
		}
	}
	
	
	/// Adds transitions from a table
	add_file_transition(p,cl,tab,so_type)
	{
		let claa = this.species[p].cla[cl];	
		claa.tra=[]; claa.ntra = 0; claa.hash_tra=[] 
		claa.nsource = 0; claa.nsink = 0;
		
		let cam = claa.camera;
		
		let equal_flag = false;
		for(let r = 0; r < tab.nrow; r++){	
			if(r%10) percent(100*r/tab.nrow);
		
			let i_str = tab.ele[r][0];
			let i;
			if(i_str == "+") i = SOURCE;
			else i = hash_find(claa.hash_comp,i_str);
				
			let f_str = tab.ele[r][1];
			let f;
			if(f_str == "-") f = SINK;
			else f = hash_find(claa.hash_comp,f_str);
			
			let trans_def = extract_trans_def(tab.ele[r][tab.ncol-1])
			
			let pos = dist_pos;	if(i == SOURCE) pos = source_dist_pos;

			option_error("type",trans_def.type,pos);
			
			if(i == SOURCE && f == SINK) alert_help("Cannot have '+->-'");
				
			if(i == f) equal_flag = true;
			
			if(i != undefined && f != undefined && i != f){	
				let midp=[];
	
				if(so_type == "Trans File Pos"){
					let x_spl = tab.ele[r][2].split(",");
					let y_spl = tab.ele[r][3].split(",");
					if(x_spl.length != y_spl.length){
						alert_help("On line "+(r+1)+" the x and y positions do not have the same number of points");
					}
					for(let k = 0; k < x_spl.length; k++){
						let numx = Number(x_spl[k]);
						let numy = Number(y_spl[k]);
						if(isNaN(numx)){
							alert_help("On line "+(r+1)+" the x positions is not a number");
						}
						
						if(isNaN(numy)){
							alert_help("On line "+(r+1)+" the y positions is not a number");
						}
						
						let x,y
						switch(cam.coord){
						case "cartesian":				
							midp.push({x:numx*import_scale_factor, y:numy*import_scale_factor});						
							break;
						
						case "latlng":
							let pt = transform_latlng(numx,numy);
							midp.push({x:pt.x, y:pt.y});	
							break;
						}
					}
				}
				else{
					if(i == SOURCE){
						let co = claa.comp[f];
						switch(cam.coord){
						case "cartesian":	midp = [{x:co.x,y:co.y+compartment_height+1}]; break;
						case "latlng": midp = [{x:co.x,y:co.y}]; break;
						default: error("Option not recognised 95"); break;
						}
					}
					
					if(f == SINK){
						let co = claa.comp[i];
						switch(cam.coord){
						case "cartesian":	midp = [{x:co.x,y:co.y-compartment_height-1}]; break;
						case "latlng": midp = [{x:co.x,y:co.y}]; break;
						default: error("Option not recognised 95"); break;
						}
					}
				}
				
				let k = 0; 
				while(k < claa.tra.length && !(claa.tra[k].i == i && claa.tra[k].f == f)) k++;
				
				let value;
				if(k < claa.tra.length){
					let traa = claa.tra[k];
					traa.type = trans_def.type;
					traa.midp = midp;
					value = traa.value;
				}
				else{
					let res = this.add_transition(p,cl,i,f,midp,trans_def.type,"simple");
					output_help(res);
					
					value = claa.tra[claa.tra.length-1].value;
				}
			
				switch(trans_def.type){
				case "exp(rate)":
					value.rate_eqn.te = trans_def.rate;
					break;
					
				case "exp(mean)":
					value.mean_eqn.te = trans_def.mean;
					break;
				
				case "erlang": 
					value.mean_eqn.te = trans_def.mean;
					value.shape_erlang.te = trans_def.shape;
					break;
					
				case "gamma":
					value.mean_eqn.te = trans_def.mean;
					value.cv_eqn.te = trans_def.cv;
					break;
					
				case "log-normal": 
					value.mean_eqn.te = trans_def.mena;
					value.cv_eqn.te = trans_def.cv;
					break;
					
				case "weibull": 
					value.scale_eqn.te = trans_def.scale;
					value.shape_eqn.te = trans_def.shape;
					break;

				case "period": 
					value.mean_eqn.te = trans_def.time;
					break;
					
				default: error("Option not recognised 90"); break;
				}
			}
		}
		
		this.update_check(p,cl);
	
		if(equal_flag) alert_help("Warning!","Some transitions are set with the same 'from' and 'to' compartments. These are ignored");
	}

	
	/// Adds a classification background buttons (this is activated when clicked)
	add_classback_buts(lay)
	{
		if(this.species.length == 0) return;
		let p = this.get_p();
		if(this.species[p].ncla == 0) return;
		let cl = this.get_cl();
		lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, ac:"ClassificationBack", type:"ClassBack", p:p, cl:cl});
	}
	
	
	/// Replots the model
	replot(sh)
	{	
		if(!layer_exist("Annotation")){	
			if(layer_exist("GraphAnnotations")) replot_layer("GraphAnnotations",sh);
			if(layer_exist("GraphCompartments")) replot_layer("GraphCompartments",sh);
			if(layer_exist("GraphTransitions")) replot_layer("GraphTransitions",sh);
		}			
		else{
			replot_layer("Annotation",sh);
			replot_layer("AnnotationMap",sh);
			if(layer_exist("PointLabel")) replot_layer("PointLabel",sh);
			replot_layer("Compartment",sh);
			replot_layer("Transition",sh);
			if(inter.comp_select.sbox && inter.comp_select.list.length > 1){
				replot_layer("CompSelect",sh);
				replot_layer("CompSelectButton",sh);
			}
		}
		
		plot_screen();
	}
	
	
	/// Gets all the values an index could take
	get_index_possibility(def_pos)
	{
		let pos = [];
		for(let i = 0; i < alphabet.length; i++){
			let ch = alphabet[i];
			let flag = false;
			if(ch == "t") flag = true;
			if(ch == "z") flag = true;
			
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


	/// Gets all the label info for the parameters
	get_label_info_all()
	{
		for(let th = 0; th < this.param.length; th++) this.get_label_info(this.param[th]);
	}


	/// Gets label information about a parameter
	get_label_info(par)
	{	
		let fo_big = get_font(si_big,undefined,"Times");
		let fo_sup = get_font(si_sub,undefined,"Times");
		let fo_sub = get_font(si_sub,"italic","Times");

		let te = par.name;
		let spl = te.split("^");
		par.name_raw = spl[0];
		par.sup = ""; if(spl.length == 2)	par.sup = remove_bracket(spl[1]); 
		
		par.name_raw_w = text_width(par.name_raw,fo_big);
		par.sup_w = text_width(par.sup,fo_sup);
		
		let sub = "";
		for(let j = 0; j < par.dep.length; j++){
			if(par.dep[j] != "t"){
				if(sub != "") sub += ",";
				sub += par.dep[j];
			}
		}
		par.sub = sub;
		par.sub_w = text_width(par.sub,fo_sub);
		
		let w = par.name_raw_w;
		let w2 = par.sup_w; if(par.sub_w > w2) w2 = par.sub_w;
		let w_tot = w+w2+0.1;
		
		let time = "";
		if(par.time_dep == true) time = "(t)";
		
		let time_pos = w_tot;
		if(time != ""){
			w_tot += text_width(time,fo_big);
		}
		
		par.label_info = {name_full:par.name, name:par.name_raw, time_pos:time_pos, time:time, sup:par.sup, sub:par.sub, fo_big:fo_big, fo_sup:fo_sup, fo_sub:fo_sub, dx:w_tot};
	}


	/// Loads up a new model based on information from a webworker
	load(ans,not_zero_page)
	{
		let exa_store;
		
		if(this.example && not_zero_page) exa_store = this.example;
	
		let mod = ans.model;
		
		for(let ele in mod) this[ele] = mod[ele];
		
		this.get_label_info_all();
		
		if(ans.map_store) map_store = ans.map_store;
		
		initialise_pages();
		zero_page_index(not_zero_page);
	
		if(ans.sim_load) model_sim = this.create_ans_model();
		if(ans.inf_load) model_inf = this.create_ans_model();
	
		// Switches page depeding on what model has been loaded
		if(ans.inf_load == true){
			if(ans.ppc_load == true){
				change_page({pa:"Post. Simulation", su:"Results"});
			}
			else{
				change_page({pa:"Inference", su:"Results"});
			}
			inter.options = true;	generate_screen();
		}
		else{
			if(ans.sim_load == true){
				change_page({pa:"Simulation", su:"Results"});
			}
			else{
				change_page({pa:"Model", su:"Compartments"});
			}
		}
		
		if(exa_store) this.example = exa_store;
	
		model.check_show_warning(); // Checks if we should show a model warning
	}
	
	
	/// Gets all individuals mentioned in the data
	get_all_data_individual(p)
	{
		let hash = new Hash();
		
		let list=[];
		
		let inf_source = model.species[p].inf_source;
		for(let i = 0; i < inf_source.length; i++){
			let so = inf_source[i];
			switch(so.type){
			case "Add Pop.": case "Remove Pop.": 
			case "Add Ind.": case "Remove Ind.": case "Move Ind.": 
			case "Compartment":case "Transition": case "Trans. in time range":
			case "Diag. Test": case "Genetic":
				{
					let tab = so.table;
					for(let r = 0; r < tab.nrow; r++){
						let name = tab.ele[r][0];
						if(hash.find(name) == undefined){
							hash.add(name,list.length);
							list.push({name:name});
						}
					}
				}
				break;
			}
		}
		
		let ind_eff_group = model.species[p].ind_eff_group;
		
		for(let i = 0; i < ind_eff_group.length; i++){
			let A = ind_eff_group[i].A_matrix;
			if(A.check == true && A.loaded == true){
				for(let j = 0; j < A.ind_list.length; j++){
					let name = A.ind_list[j];
					if(hash.find(name) == undefined){
						hash.add(name,list.length);
						list.push({name:name});
					}
				}					
			}
		}
	
		return list;
	}
	
	
	/// Determines if a warning should be shown
	show_warning()
	{
		if(tab_name() == "Home") return false;
		
		if(model.warn.length > 0){
			if(model.warn_view == true) return true;
			else return "button";
		}
		return false;
	}

	
	/// Checks that a name is not already existing in the model (list stores exceptions)
	check_name(te,list)
	{
		if(check_invalid_name(te)) return "Cannot use reserved word '"+te+"'";
	
		for(let i = 0; i < te.length; i++){
			let ch = te.substr(i,1);
		
			if(name_notallow.includes(ch)){
				if(ch == " ") return "Cannot contain any spaces";
				else return "Cannot use the character '"+ch+"'";
			}
		}
		
		for(let p = 0; p < this.species.length; p++){
			let sp = this.species[p];
			
			if(sp.name.toLowerCase() == te.toLowerCase()){	
				let k = 0; 
				while(k < list.length && !(list[k].p == p && list[k].cl == undefined && list[k].c == undefined)) k++; 
				if(k == list.length){
					return "Must have a different name to species '"+sp.name+"'";
				}
			}
			
			for(let cl = 0; cl < sp.cla.length; cl++){
				let claa = sp.cla[cl];	
				
				if(claa.name.toLowerCase() == te.toLowerCase()){
					let k = 0; 
					while(k < list.length && !(list[k].p == p && list[k].cl == cl && list[k].c == undefined)) k++; 
					if(k == list.length){
						let st = "Must have a different name to classification '";
						if(this.species.length > 1) st += sp.name+"→";
						st += claa.name+"'";
						return st;
					}
				}				

				for(let c = 0; c < claa.ncomp; c++){
					let co = claa.comp[c];
					if(co.name.toLowerCase() == te.toLowerCase()){
						let k = 0; 
						while(k < list.length && !(list[k].p == p && list[k].cl == cl && list[k].c == c)) k++; 
						if(k == list.length){
							let st = "Must have a different name to compartment '";
							if(this.species.length > 1) st += sp.name+"→";
							if(sp.cla.length > 1) st += claa.name+"→";
							st += co.name+"'";
							return st;
						}
					}
				}
			}
		}
	}
	
	
	/// Selects compartments based on a mouse selection
	select_compartments()
	{
		let box = select_box();
		box.x1 -= menu_width; box.x2 -= menu_width; 
		
		let im = inter.mode;
		let lay = inter.layer[im.l];
		
		inter.comp_select={p:im.p, cl:im.cl, list:[]};
		let cs_list = inter.comp_select.list;
	
		let claa = this.species[im.p].cla[im.cl];
		let cam = claa.camera;  
				
		for(let k = 0; k < claa.ncomp; k++){
			let c = claa.comp[k];
			let pt;
			switch(c.type){
			case "box": pt = trans_point(c.x,c.y,cam,lay); break;
			case "latlng": pt = trans_point(c.x,c.y,cam,lay); break;
			case "boundary":
				{
					let ms = find_map_store(c.map_ref);
					let box = ms.feature.box;
					pt = trans_point((box.xmin+box.xmax)/2,(box.ymin+box.ymax)/2,cam,lay);
				}
				break;
			}
			
			let x = pt.x, y = pt.y;
			if(x > box.x1 && x < box.x2 && y > box.y1 && y < box.y2){
				cs_list.push(k);
			}
		}
		
		if(cs_list.length == 1){                      // If only one selected the click on
			this.select_single();
		}
	}


	/// Determines if a given suffix is valid
	allow_suffix(suf)
	{
		let cs = inter.comp_select;
		
		let hash = new Hash();
		let claa = model.species[cs.p].cla[cs.cl];
		for(let i = 0; i < cs.list.length; i++){
			let co = claa.comp[cs.list[i]];
			hash.add(co.name+suf,i);
		}
		
		for(let p = 0; p < model.species.length; p++){
			let sp = model.species[p];
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				for(let c = 0; c < claa.ncomp; c++){
					if(hash.find(claa.comp[c].name) != undefined){
						error(claa.comp[c].name+" not allow");
						return false;
					}
				}
			}
		}
		
		return true;
	}
	
	
	/// Selects a single selected compartment
	select_single()
	{
		let cs = inter.comp_select;
		let cs_list = cs.list;
	
		if(cs_list.length != 1){ error("Should only be one"); return;}
		
		let lay = get_lay("Compartment");
		
		for(let i = 0; i < lay.but.length; i++){
			let bu = lay.but[i];
			
			switch(bu.type){
			case "Compartment": case "CompLatLng":
				if(bu.p == cs.p && bu.cl == cs.cl && bu.i == cs_list[0]){
					inter.over = {layer:lay.index, i:i};
					button_action(bu,"click");
					return;
				}
				break;
			}
		}
		error("Could not find single");
	}
	
	
	/// Creates a simplified model used for outputting simulation/inference/PPC
	create_ans_model()
	{
		let fac = 0.85;
		
		let model_new = {species:[]};
		for(let p = 0; p < this.species.length; p++){
			let sp = this.species[p];
			model_new.species[p] = {cla:[]};
			for(let cl = 0; cl < sp.ncla; cl++){
				let claa = sp.cla[cl];
				model_new.species[p].cla[cl] = copy(claa);
				let claa_new = model_new.species[p].cla[cl];
				claa_new.camera.scale = fac*Number(claa_new.camera.scale);
			}
		}
		return model_new;
	}
	

	///	Moves a group of compartments
	move_group(mgroup,cl,p)
	{
		let repl = false;
		for(let i = 0; i < mgroup.length; i++){
			let mg = mgroup[i];
			if(this.move_comp(mg.c,mg.x,mg.y,p,cl) == true){
				repl = true;
			}
		}

		if(repl){
			inter.mode.moved = true;
			let claa = this.species[p].cla[cl];
			this.update_pline(claa);  
			let clone = inter.mode.clone;
			if(clone){
				for(let k = 0; k < clone.length; k++){
					let clo = clone[k];
					let claa2 = this.species[clo.p].cla[clo.cl];
					this.update_pline(claa2);  
				}
			}
			
			if(inter.bubble.lay_name != undefined) generate_screen();
			else this.replot();
		}
	}
		
		
	/// Moves a compartment (with potentiall any sources or sinks also)
	move_comp(c,x,y,p,cl,op)
	{
		let claa = this.species[p].cla[cl];
		let co = claa.comp[c];
		let x_old = co.x;
		let y_old = co.y;
		co.x = x;
		co.y = y;
		
		if(claa.camera.coord == "cartesian") this.snap_comp_to_grid(co);
		
		let shx = co.x-x_old;
		let shy = co.y-y_old;
		if(shx != 0 || shy != 0){
			let mode = inter.mode;
			
			if(mode != undefined){			
				// Creates information which can be used to determine transition entering leaving
				if(!mode.comp_source_sink){
					mode.comp_source_sink = {species:[]};
					for(let p2 = 0; p2 < this.species.length; p2++){
						let sp2 = this.species[p2];
						mode.comp_source_sink.species[p2]={cla:[]};
					}
				}
				
				let css_sp = mode.comp_source_sink.species[p];
				if(css_sp.cla[cl] == undefined){
					css_sp.cla[cl] = {comp:[]};
					let css_cl = css_sp.cla[cl];
					for(let c = 0; c < claa.ncomp; c++){
						css_cl.comp[c] = {list:[]};
					}
					for(let tr = 0; tr < claa.ntra; tr++){
						let tra = claa.tra[tr];
						if(tra.i == SOURCE) css_cl.comp[tra.f].list.push(tr);
						if(tra.f == SINK) css_cl.comp[tra.i].list.push(tr);
					}
				}
			
				let list = css_sp.cla[cl].comp[c].list;
				for(let i = 0; i < list.length; i++){ 
					let tra = claa.tra[list[i]];
					for(let k = 0; k < tra.midp.length; k++){
						tra.midp[k].x += shx;
						tra.midp[k].y += shy;
					}
				}
		
				if(mode.clone == undefined){
					mode.clone = model.find_clones(p,cl);	
				}
			
				if(op != "no clone"){
					for(let loop = 0; loop < mode.clone.length; loop++){
						let clo = mode.clone[loop];
						this.move_comp(c,x,y,clo.p,clo.cl,"no clone");
					}
				}
			}
			
			return true;
		}				
		return false;
	}
	

	/// Ensures that camera stays within the lng range
	ensure_lng(cam,ww)
	{
		if(cam.coord == "latlng"){
			let no_plot = false;
			if(ww == undefined){
				let show = model.get_show_model();
				if(show) ww = get_lay("Main").dx;
				else ww = get_lay("GraphCompartments").dx;
			}
			else no_plot = true;
		
			let max = Math.PI - (ww/2)/cam.scale;
			let min = -Math.PI+ (ww/2)/cam.scale;
			
			let fl = false;
			if(min > max){
				cam.x = 0;
				cam.scale = ww/(2*Math.PI);
				fl = true;
			}
			else{
				if(cam.x < min){ cam.x = min; fl = true;}
				if(cam.x > max){ cam.x = max; fl = true;}
			}
			
			if(fl && no_plot == false){
				this.replot();
			}
		}
	}
	
	
	/// Checks if we should show a model warning
	check_show_warning()
	{
		let na = inter.page_name;
		let tree = na.split("->");
		
		let res;
		if(tree[1] == "Results"){
			switch(tree[0]){
			case "Simulation": res = this.sim_res; break;
			case "Inference": res = this.inf_res; break;
			case "Post. Simulation": res = this.ppc_res; break;
			default: break;
			}
		}			
		
		if(res == undefined) return;
	
		if(res.run_warning_show == true) return;

		let wa = res.run_warning;
		if(wa.length > 0){
			if(wa.length == 1){
				alert_help("Run time warning",wa[0]);
			}
			else{
				let te="";
				for(let i = 0; i < wa.length; i++){
					te += "• "+wa[i]+endl;
				}
				alert_help("Run time warnings",te);
			}
		
			res.run_warning_show = true;
		}
	}
}

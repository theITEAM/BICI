"use strict";
// Functions for setting up simulation

const spawn = require('child_process').spawn;

/// Starts creation of files for spawning
function start_spawn(file_list)
{
	let file = "Execute/init.bici"; if(ver == "mac") file = "/tmp/init.bici";
	let dir;

	create_files(file_list,file,dir,"spawn");
}

					
/// Saves all the files needed for analysis
function create_files(file_list,file,dir,type)
{
	if(dir != undefined){
		const fs = require('fs');  // Creates the data directory
		if (dir != undefined && !fs.existsSync(dir)) {
			fs.mkdir(dir, function(err) {
				if(err) { 
					alertp("There was a problem creating the directoy '"+dir+"'");
				}
			});
		}
	}
	
	inter.saving = { num:0, num_to_do:file_list.length, type:type}
		
	for(let i = 0; i < file_list.length; i++){
		let finfo = file_list[i]
		if(finfo.type == "bicifile"){
			write_file_async(finfo.data,file);
		}
		else{
			if(dir == undefined) error("dir should be defined");
			write_file_async(finfo.data,dir+"/"+finfo.file);
		}
	}

	generate_screen();
}


/// Saving files has finished
function saving_done()
{
	if(make_file){
		stop_loading_symbol();
		generate_screen();
		return;
	}
	
	switch(inter.saving.type){
	case "spawn":
		start(inter.save_type);
		break;
	default: break;
	}
}


/// Starts a simulation or inference 
function start(siminf)
{
	if(win_linux){
		setTimeout(function(){ 
			stop_loading_symbol();
			alert_help("Model saved!"); 
			generate_screen();
		}, 10);
		return;
	}
	
	start_loading_symbol(0,"Spawn");

	inter.running_status = true;

	inter.chain = [];
	
	let nchain = 1; 
	if(inter.save_type == "inf" && model.inf_details.algorithm.value =="DA-MCMC"){
		nchain = model.inf_details.nchain;
	}
	
	for(let ch = 0; ch < nchain; ch++){
		startspawn(ch,nchain,siminf);
	}
}


/// Gets the value for the seed
function get_seed_not_set(siminf)
{
	let details;
	switch(siminf){
	case "sim": details = model.sim_details; break;
	case "inf": details = model.inf_details; break;
	case "ppc": details = model.ppc_details; break;
	}

	if(!turn_off_random_seed){
		if(details.seed_on.value != "Yes") return Math.floor(Math.random()*seed_max); 
	}
}


/// Starts execution of C++ code
function startspawn(ch,nchain,siminf)                                    
{	
	let seed_not_set = get_seed_not_set(siminf);

	inter.chain[ch] = { done:false, siminf:siminf, term:false, prog:0, leftover:"", lines:[]};

	let do_com = siminf; if(do_com == "ppc") do_com="post-sim";

	let file = "bici-core.exe";
	if(ver == "mac") file = "./bici-core";
	
	if(seed_not_set != undefined){
		inter.child[ch] = spawn(file,["default.bici",do_com,"-chain="+ch,"-seed="+seed_not_set]);
	}
	else{
		inter.child[ch] = spawn(file,["default.bici",do_com,"-chain="+ch]);
	}		
			
	funct(inter.child[ch],ch);
}


/// This function executes the c++ code
function funct(chi,ch)                             // Gathers output of C++ file
{
	chi.stdout.on('data', function (data) {
		let cha = inter.chain[ch];
		let st = cha.leftover + data;

		let lines = st.split('\n');
		cha.leftover = lines[lines.length-1];
		
		for(let li = 0; li < lines.length-1; li++){
			let line = lines[li];
			
			let fl = false;
			if(line.length > 1 && line.substr(0,1) == "<"){
				if(begin(line,"<CREATING>")){
					loading_symbol_message("Creating...");
					fl = true;
				}
				
				if(begin(line,"<INIT>")){
					loading_symbol_message("Initialising...");
					fl = true;
				}
				
				if(begin(line,"<RUNNING>")){
					loading_symbol_message("Running...");
					fl = true;
				}
				
				if(begin(line,"<RUNGEN>")){
					let gen = Number(line.substr(8));
					loading_symbol_message("Generation "+gen);
					fl = true;
				}
				
				if(begin(line,"<ANNEALING>")){
					loading_symbol_message("Annealing...");
					fl = true;
				}
		
				if(begin(line,"<PROGRESS>")){
					cha.prog = Number(line.substr(10));
					let progmin = LARGE;
					for(let ch2 = 0; ch2 < inter.chain.length; ch2++){
						if(inter.chain[ch2].prog < progmin) progmin = inter.chain[ch2].prog;
					}
					set_loading_percent(progmin);
					fl = true;
				}
					
				if(begin(line,"<OUTPUTTING>")){
					loading_symbol_message("Outputting...");
					fl = true;
				}
			}
			
			if(fl == false){
				cha.lines.push(line);
			}
		}
	});

	chi.stderr.on('data', function (data) {
		terminate("The BICI core code terminated unexpectedly"); 
	});

	chi.on('close', function (code) {
		if(inter.running_status == true){
			let cha = inter.chain[ch];
			cha.done = true;
			
			// Extracts the file from the output
			let lines = cha.lines;
			for(let c = 0; c < lines.length; c++){ lines[c] = lines[c].replace(/\r/g, "");}
			
			let i = 0;
			while(i < lines.length && lines[i] != "<<ERROR>>" && lines[i] != "<<OUTPUT FILE>>") i++;
		
			if(i == lines.length){
				terminate("There was a problem running the code."); 
			}
			else{
				let li = lines[i];
				
				i++;
				let content="";
				while(i < lines.length && lines[i] != "<<END>>"){
					content += lines[i]+endl;
					i++;
				}
			
				switch(li){
				case "<<OUTPUT FILE>>": 
					cha.content = content;
					break;
					
				case "<<ERROR>>":
					terminate(content.trim());
					return;
				} 
			}
		
			let j = 0; 
			while(j < inter.chain.length && inter.chain[j].done == true && inter.chain[j].term == false) j++;
			if(j == inter.chain.length){
				setTimeout(function(){ process_all_chains()}, 10);
			}
		}
	});
}	


/// Processes information from all chains
function process_all_chains()
{
	loading_symbol_message("Processing...");
	
	inter.running_status = false;
	
	let content = inter.chain[0].content;
	for(let ch = 1; ch < inter.chain.length; ch++){
		content += endl+inter.chain[ch].content;
	}
	 
	set_loading_percent(1);
	start_worker("Spawn Output",{content:content});	
}


/// Terminates the code
function terminate(te)
{
	if(inter.running_status == false) return;
	
	for(let ch = 0; ch < inter.chain.length; ch++){
		let cha = inter.child[ch];
		inter.chain[ch].term = true;
		inter.child[ch].stdin.pause();
		inter.child[ch].kill();
	}
	
				
	inter.running_status = false;
	
	if(te){
		inter.help = { title:"Error running core code", te:"The following error was generated which caused BICI to terminate:\n'ERROR: "+te+"'"};
	}
	
	setTimeout(function(){ stop_loading_symbol(); generate_screen();}, 10);
}


/// Adds buttons to start a simulation		
function add_sim_start_buts(lay)
{
	//if(inter.running_status == true) return;
	
	let cx = corner.x;
	let cy = corner.y;
	
	if(inter.options == false){
		{
			cy = lay.add_title("Start simulation",cx,cy,{te:start_sim_text});
		
			cy += 2;
		
			cy = lay.add_subtitle("Time range",cx,cy,WHITE,{te:time_range_text});
			cy = lay.add_paragraph("Set the time range over which simulation is performed:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
			let yy = cy-2.5;
			add_right_input_field(yy,"Start time",{type:"sim_t_start",update:true},lay);
			add_right_input_field(yy+3.5,"End time",{type:"sim_t_end",update:true},lay);
		}
		
		cy += 5;
		
		{
			cy = lay.add_subtitle("Time-step",cx,cy,WHITE,{te:sim_timestep});

			if(model.sim_details.algorithm == undefined) model.sim_details.algorithm = {value:"gillespie"};
		
			let xx = 33;
		
			cy = lay.add_paragraph("Set the time-step used to simulate:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);

			let yy = cy-2.5;
			add_right_input_field(yy,"Time-step",{type:"sim_timestep",update:true},lay);
		}
		
		cy += 2;
		
		{
			cy = lay.add_subtitle("Simulation number",cx,cy,WHITE,{te:sim_num_text});
			cy = lay.add_paragraph("Set the number of simulations to be generated:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
			let yy = cy-2.5;
			add_right_input_field(yy,"Number",{type:"sim_number",update:true},lay);
		}
		
		cy += 2;
		
		cy = run_local(cx,cy,model.sim_details,lay);
		
		add_corner_link("> Further options","Options",lay);
		
		lay.add_corner_button([["Start","Grey","StartSimulation"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	}
	else{
		let te = "Further simulation options";
		cy = lay.add_title(te,cx,cy,{te:further_sim_text});
		
		cy += 2;
		
		let gap = 2;
		
		{
			cy = lay.add_subtitle("Individual max",cx,cy,WHITE,{te:indmax_text});
		
			cy = lay.add_paragraph("Maximum number of individuals (for individual-based models)",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
	
			let yy = cy-2.5;
			add_right_input_field(yy,"Maximum",{type:"sim_indmax",update:true},lay);
		}
		
		cy += gap;
		
		{
			cy = lay.add_subtitle("Parameter outputs",cx,cy,WHITE,{te:paramoutputmax_text});
				
			cy = lay.add_paragraph("Threshold number of tensor elements above which tensor not output.",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
	
			let yy = cy-2.5;
			add_right_input_field(yy,"Maximum",{type:"sim_paramout",update:true},lay);
		}
		
		cy += gap;
		
		cy = set_seed(cx,cy,"sim_seed",model.sim_details,lay);
	
		cy += gap;
		
		lay.add_corner_button([["Done","Grey","OptionsDone"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	}
}


/// Sets button if running locally
function run_local(cx,cy,details,lay)
{
	cy = lay.add_subtitle("Run on local machine",cx,cy,WHITE,{te:run_local_text});
		
	let xx = 3, gap = 0.5;
	xx = lay.add_radio(xx,cy+0.1,"Yes","Yes",details.run_local,{back_col:WHITE});
	xx += gap;
	xx = lay.add_radio(xx,cy+0.1,"No","No",details.run_local,{back_col:WHITE});		
	
	return cy;
}
		
		
/// Sets the seed option
function set_seed(cx,cy,ref,details,lay)
{
	cy = lay.add_subtitle("Set random seed",cx,cy,WHITE,{te:seed_text});
		
	let xx = 3, gap = 0.5;
	xx = lay.add_radio(xx,cy+0.1,"Yes","Yes",details.seed_on,{back_col:WHITE});
	xx += gap;
	xx = lay.add_radio(xx,cy+0.1,"No","No",details.seed_on,{back_col:WHITE});		
	
	if(details.seed_on.value == "Yes"){
		let yy = cy-1.5;
		add_right_input_field(yy,"Seed",{type:ref,update:true},lay);
	}
	
	return cy;
}


/// Adds buttons to start a posterior predictive check
function add_ppc_start_buts(lay)
{
	if(inter.running_status == true) return;
	
	let cx = corner.x;
	let cy = corner.y;
	inter.options = false;
	if(inter.options == false){
		cy = lay.add_title("Start posterior predictive check",cx,cy,{te:start_sim_text});
		
		cy += 2;
		
		cy = lay.add_subtitle("Time range",cx,cy,WHITE,{te:time_range_text});
		cy = lay.add_paragraph("Set the time range over which simulation is performed:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
		let yy = cy-2.5;
		add_right_input_field(yy,"Start time",{type:"ppc_t_start",update:true},lay);
		add_right_input_field(yy+4,"End time",{type:"ppc_t_end",update:true},lay);

		cy += 6;
		
		let xx = 33;
		
		cy += 2;
		
		cy = lay.add_subtitle("Simulation number",cx,cy,WHITE,{te:ppc_num_text});
		cy = lay.add_paragraph("Set the number of simulations to be generated:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
		yy = cy-2.5;
		add_right_input_field(yy,"Number",{type:"ppc_number",update:true},lay);
		
		cy += 2;
	
		cy = run_local(cx,cy,model.ppc_details,lay);
	
		if(model.ppc_details.check_box_list){
			if(model.ppc_details.check_box_list.length > 0){
				cy += 2;
				cy = lay.add_subtitle("Resampling",cx,cy,WHITE,{te:ppc_num_text});
				cy = lay.add_paragraph("Select which parameter or individual effects get resampled:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
				
				cy += 0.5;
				
				add_layer("SamplingCheckbox",lay.x+cx,lay.y+cy,lay.dx-2*cx-20,lay.dy-cy-3.5,{});
			}
		}
		
		lay.add_corner_button([["Start","Grey","StartPPC"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	}
	else{
		let te = "Further posterior simulation options";
		cy = lay.add_title(te,cx,cy,{te:further_post_sim_text});
		
		cy += 2;
		
		let gap = 2;
			
		cy = set_seed(cx,cy,"PPC_seed",model.ppc_details,lay);
	}
}


/// Shows check boxes for resampling
function sampling_checkbox(lay)
{
	let cbl = model.ppc_details.check_box_list;
	let cy = 0;
	for(let i = 0; i < cbl.length; i++){
		let cb = cbl[i];
		lay.add_checkbox(1,cy,cb.name,cb.name,cb.checkb);
		cy += 1.2;
	}	
}


/// Adds a text input field on the right hand menu
function add_right_input_field(yy,te,op,lay,shift)
{
	let dx = 8;
	let xx = lay.inner_dx-dx-3;
	if(shift) xx -= shift;
	add_input_field(xx,yy,dx,te,op,lay);
}


/// Adds a text input field on the right hand menu
function add_right_dropdown(yy,te,so,pos,lay)
{	
	let dx = 8;
	let xx = lay.inner_dx-dx-3;
	
	lay.add_button({te:te, x:xx, y:yy, dx:dx, dy:0.8, type:"InputBoxName", back_col:WHITE});
	yy += 1;
	
	lay.add_dropdown(xx,yy,dx,10,so,pos);
	
	yy += 2;
}


/// Adds an input box
function add_input_field(xx,yy,dx,te,op,lay)
{
	lay.add_button({te:te, x:xx, y:yy, dx:dx, dy:0.8, type:"InputBoxName", back_col:WHITE});
	yy += 1.2; if(op.narrow == true) yy -= 0.2;
	lay.add_input(xx,yy,dx,op);
	
	yy += 2;
	let sto = inter.textbox_store;
	let k = 0; while(k < sto.length && sto[k].ref != op.type) k++;
	
	if(k < sto.length){
		let warn = sto[k].warning;
		if(warn != undefined){
			lay.add_paragraph(warn,dx,xx,yy-0.4,RED,warn_si,warn_lh,undefined,"center");
		}
	}
}
		

/// Adds the button for the page in which simulation parameters are input
function add_param_value_buts(lay)
{
	let cx = corner.x;
	let cy = corner.y;
	
	cy = lay.add_title("Simulation parameter values",cx,cy, {te:sim_text});
		
	cy += 0.5;
		
	add_layer("ParamValueContent",lay.x+cx,lay.y+cy,lay.dx-2*cx,lay.dy-cy-2,{});	
}



/// Gets a list of categories from the parameter list
function get_param_cat(filt_list,filt_type)
{
	let param = model.param;
	
	let param_cat = [];
	param_cat.push({name:"Parameters", sim_te:sim_uni_text, inf_te:inf_uni_text, list:[]});
	param_cat.push({name:"Parameter vectors", sim_te:sim_vector_text, inf_te:inf_vector_text, list:[]});
	param_cat.push({name:"Parameter matrices", sim_te:sim_matrix_text, inf_te:inf_matrix_text, list:[]});
	param_cat.push({name:"Parameter tensors", sim_te:sim_tensor_text, inf_te:inf_tensor_text, list:[]});
	param_cat.push({name:"Individual effects", sim_te:sim_ie_variance_text, inf_te:inf_ie_variance_text, list:[]});
	param_cat.push({name:"Fixed effects", sim_te:fixed_eff_text, inf_te:fixed_eff_text, list:[]});
	
	for(let i = 0; i < param.length; i++){
		if(find_in(filt_list,param[i].type) == undefined && 
		!(filt_type == "only normal" && param[i].variety != "normal") &&
		!(filt_type == "for sim" && param[i].variety == "const") &&
		!(filt_type == "for sim" && param[i].variety == "reparam") && 
		param[i].type != "param factor"
		){
			if(param[i].type == "fixed effect"){
				 param_cat[5].list.push(i);
			}
			else{
				if(param[i].type == "variance" || param[i].type == "correlation") param_cat[4].list.push(i);
				else{
					let n = param[i].dep.length; if(n > 3) n = 3;
					param_cat[n].list.push(i);
				}
			}
		}
	}
	
	return param_cat;
}


/// Adds a screen allowing simulation parameters to be edited
function add_param_value_content(lay)
{
	let si_small = 0.8;
	let fo_small = get_font(si_small);
	
	let y = 0;
	let dy = 2.0;

	let vector_flag = false;
	
	let x = maximim_label_width() + 4;

	let param = model.param;
	
	let param_cat = get_param_cat(sim_param_not_needed,"for sim");
	
	if(no_param(param_cat)) center_message("No parameters need to be set.",lay);
		
	for(let cat = 0; cat < param_cat.length; cat++){
		let pc = param_cat[cat];
		if(pc.list.length > 0){
			y = lay.add_subtitle(pc.name,1,y+0.2,WHITE,{te:pc.sim_te});
		
			for(let cati = 0; cati < pc.list.length; cati++){
				let i = pc.list[cati];
				let par = param[i];
			
				if(par.variety != "const" && par.variety != "reparam"){
					switch(par.type){
					case "Se": case "Sp":// case "trap_prob": 
					case "mut_rate": case "seq_var":
					case "comp_prob": case "trans_prob": case "derive_param":
						break;
					
					default:
						{
							let w = wright2;
						
							if(par.value_desc == no_elements){
								display_no_element(par,x,y,lay,w);
							}
							else{
								if(add_view_button(par,w-4,y,i,lay,model)) w -= 4.5;
								
								if(par.variety == "dist"){
									w -= 4.5;
									lay.add_checkbox(w,y+0.4,"Sample","Sample",par.sim_sample,WHITE);
								}
								
								if(par.variety == "dist" && par.sim_sample.check == true){
									if(par.dep.length == 0 || par.prior_split_check.check == false){
										display_distribution(i,x,y,lay,true,false,w);
									}
									else{
										if(par.prior_split_desc == no_elements){
											display_no_element(par,x,y,lay,w);
										}
										else{
											display_distribution_split(i,x,y,lay,true,false,"dist",w);
										}
									}
								}
								else{
									display_constant(i,x,y,lay,w);
								}
							}
							y += dy;
						}
						break;
					}
				}
			}
		}
	}
}


/// Determines if there are no parameters in param_cat
function no_param(param_cat)
{
	for(let cat = 0; cat < param_cat.length; cat++){
		if(param_cat[cat].list.length > 0) return false;
	}
	return true;
}


/// Displays a row showing that a parameter down not have any elements
function display_no_element(par,x,y,lay,w)
{
	lay.display_param(x-par.label_info.dx-0.7,y-0.1,par.label_info);

	let si = 1.5;
	lay.add_button({te:"=", x:x, y:y, dy:si, type:"Text", font:get_font(si), si:si, col:BLACK});
	
	let fo = get_font(1.1,"","times");
		
	lay.add_button({te:no_elements, x:x+1.6, y:y+0., dx:w-x-1.6, dy:1.6, type:"No element",font:fo});
}


/// Displays a row allowing a variable constant to be set
function display_constant(i,x,y,lay,w,allow_edit,source)
{
	if(source == undefined) source = model;
	let par = source.param[i];
	
	lay.display_param(x-par.label_info.dx-0.7,y-0.1,par.label_info);

	let si = 1.5;
	lay.add_button({te:"=", x:x, y:y, dy:si, type:"Text", font:get_font(si), si:si, col:BLACK});
	
	let te="";
	if(par.dep.length == 0){
		te = par.value;
		if(allow_edit == false && isNaN(te)) te = "Auto calculate"; 
	}
	else{
		if(par.set == false) te = "Unset ";
		
		if(par.dist_mat) te = "Distance "; 
		te += par.value_desc;
	}
	
	let fo = get_font(1.1,"","times");
		
	let ac; if(allow_edit != false) ac = "EditSimValue";
	
	lay.add_button({te:te, x:x+1.6, y:y+0., dx:w-x-1.6, dy:1.6, type:"ParamSimElement", source:source, font:fo, ac:ac, i:i, name:par.name, label_info:par.label_info});
}


/// Displays a row allowing a variable factor to be set
function display_factor(i,x,y,lay,w,allow_edit,source)
{
	if(source == undefined) source = model;

	let par = model.param[i];
	
	lay.display_param(x-par.label_info.dx-0.7,y-0.1,par.label_info);

	if(par.factor_weight_on.check){
		let te = par.weight_desc;
	
		let si = 1.0, fo_la = get_font(si,"","times");			
		lay.add_button({te:"Weight =", x:x+1.8, y:y+0.3, dy:si, si:si, font:fo_la, type:"Text", col:BLACK});
		
		let fo = get_font(1.1,"","times");
		
		let ac = "EditWeightValue";
		lay.add_button({te:te, x:x+5.9, y:y+0., dx:w-x-9.6, dy:1.6, type:"ParamSimElement", source:source, font:fo, ac:ac, i:i, name:par.name, label_info:par.label_info});
	}
	
	lay.add_checkbox(w-4,y+0.3,"Weight","Weight",par.factor_weight_on,WHITE,{title:"Factor weight", te:factor_weight_text});
}


/// Displays a row allowing a variable constant to be set
function display_reparam(i,x,y,lay,w)
{
	let par = model.param[i];
	
	lay.display_param(x-par.label_info.dx-0.7,y-0.1,par.label_info);

	let si = 1.5;
	lay.add_button({te:"=", x:x, y:y, dy:si, type:"Text", font:get_font(si), si:si, col:BLACK});
	
	let te="", ac = "EditReparamValue", type = "ReparamElement";
	if(par.dep.length == 0){
		te = par.value;
	}
	else{
		if(par.reparam_eqn_on == true){
		 te = par.reparam_eqn;
		 if(te == "") te = "Unset";
		 type = "ReparamEqn";
		}
		else{
			if(par.set == false) te = "Unset ";
			te += par.value_desc;
		}
	}
	
	let fo = get_font(1.1,"","times");
	
	lay.add_button({te:te, x:x+1.6, y:y+0., dx:w-x-1.6, dy:1.6, type:type, font:fo, ac:ac, i:i, name:par.name, label_info:par.label_info});
}


/// Click to edit simulation value
function edit_sim_value(th,lay_name,i,source)
{
	let par = source.param[th];

	if(par.dep.length == 0){
		select_bubble(lay_name,i,{});
	}
	else{
		start_worker("Edit Param",{type:"Value", source:source, label_info:par.label_info, i:th});
	}
}


/// Click to edit reparameterised expressions
function edit_reparam_value(th,lay_name,i,source)
{
	let par = model.param[th];

	if(par.dep.length == 0){
		select_bubble_over();
		inter.bubble.th = th;
	}
	else{
		if(par.reparam_eqn_on == true){
			select_bubble_over();
			inter.bubble.th = th;
		}
		else{
			start_worker("Edit Reparam",{type:"Reparam", source:source, label_info:par.label_info, i:th});
		}
	}
}


/// Works out the maximum width of a parameter label
function maximim_label_width()
{
	let w_max = 0;
	for(let i = 0; i < model.param.length; i++){
		let par =  model.param[i];
		if(par.label_info.dx > w_max) w_max = par.label_info.dx;
	}
	
	return w_max;
}


/// Adds a "View" button to allow for the data to be visualised
function add_view_button(par,x,y,i,lay,source)
{
	if(par.set == false) return false;
		
	if(!par.dist_mat && (par.variety == "normal" || par.variety == "const")){
		let pos_view = get_par_pos_view(par,source);
		if(pos_view.length > 0){
			lay.add_button({te:"View", source:source, x:x, y:y+0.2, dx:3.5, dy:1.2, ac:"ViewParam", type:"GreyView", i:i, pos_view:pos_view});
			return true;
		}
	}
	
	return false;
}


/// Checks there is enough memory to perform task
function check_memory(ans)
{
	let type = ans.save_type;
	
	const RAM = navigator.deviceMemory;
		
	let C_total = 0, T_total = 0;
	let C_st=[], T_st=[];
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		let N=[];
		let C=1;
		for(let cl = 0; cl < sp.cla.length; cl++){
			let claa = sp.cla[cl];
			N[cl] = claa.comp.length;
			C *= N[cl];
		}
	
		let T = 0;
		for(let cl = 0; cl < sp.cla.length; cl++){
			let claa = sp.cla[cl];	
			T += claa.tra.length*C/N[cl];
		}
		
		C_st[p] = C; T_st[p] = T;
		
		C_total += C; T_total += T;
	}
	
	let S, P; // Number of state and parameter samples
	let nchain = 1;
	
	let details;

	switch(type){
	case "sim":
		details = model.sim_details;
		S = Number(details.number);
		P = S;
		break;
	
	case "inf":
		details = model.inf_details;
		S = details.output_state;
		P = details.output_param;
		nchain = Number(details.nchain);
		break;
		
	case "ppc":
		details = model.ppc_details;
		S = Number(details.number);
		P = S;
		break;
	}
	
	let t_start = Number(details.t_start);
	let t_end = Number(details.t_end);
	if(type == "ppc"){
		t_start = Number(details.ppc_t_start);
		t_end = Number(details.ppc_t_end);
	}
	
	let timestep = Number(details.timestep);
	let TS = (t_end-t_start)/timestep;
	
	let state_mem = (C_total + T_total)*TS*10*8/bytes_in_GB;
	
	let NE_max = details.param_output_max;
	
	let param_ele_total = 0;
	for(let k = 0; k < model.param.length; k++){
		let par = model.param[k];
		if(par.dep.length == 0) param_ele_total++;
		else{
			let NE = 1; 
			let dim = par.dim;
			for(let j = 0; j < dim.length; j++) NE *= dim[j];
			if(NE < NE_max) param_ele_total += NE;
		}
	}		
	
	let param_mem = param_ele_total*P*8/bytes_in_GB;
	
	//let interface_mem = state_mem+param_mem;
	
	let run_mem = 0;
	let temp_mem = 0;
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		
		let source;
		switch(type){
		case "sim": source = sp.sim_source; break;
		case "inf": case "ppc": source = sp.inf_source; break; 
		}
		
		let NI = 0, Npopdata = 0, Npoptransdata = 0;
		let gen_data = false;
		for(let k = 0; k < source.length; k++){
			let so = source[k];
			switch(so.type){
			case "Add Ind.": NI += so.num; break;
			case "Population": Npopdata += so.num; break;
			case "Pop. Trans.": Npoptransdata += so.num; break;
			case "Genetic": gen_data = true; break;
			}
		}
		
		run_mem += T_st[p]*C_st[p]*4;   // transform in TransGlobal
	
		if(type == "inf"){
			run_mem += C_st[p]*T_st[p]*4; // trg_from_tr
			run_mem += C_st[p]*C_st[p]*4; // comp_global_convert
			
			let tm = C_st[p]*C_st[p]*8; // tr_swap (temporarily) 
			if(tm > temp_mem) temp_mem = tm;
				
			run_mem += C_st[p]*C_st[p]*4; // comp_global_convert

			tm = T_st[p]*T_st[p];   // tr_connected (only temporary)
			if(tm > temp_mem) temp_mem = tm;
		}
	
		run_mem += T_st[p]*TS*4*10;   // obs_trans_eqn_ref
		run_mem += C_st[p]*TS*4*10;   // pop_data_ref
		run_mem += T_st[p]*TS*4*10;   // pop_trans_ref
		if(gen_data) run_mem += C_total*T_total*12; // inf_cause
			
		run_mem += Npopdata*C_st[p]*8; // pop_data_cgl
		run_mem += Npoptransdata*T_st[p]*8; //pop_trans_data_tgl
		
		switch(sp.type){
		case "Individual":
			run_mem += NI*100;            // individual
			break;
		
		case "Population":			
			run_mem += T_st[p]*TS*8;      // trans_num
			run_mem += T_st[p]*TS*8;      // tnum_mean_st
			run_mem += T_st[p]*TS*8;      // Li_markov_pop
			run_mem += C_st[p]*TS*8;      // cpop_st
			run_mem += T_st[p]*TS*8;      // trans_num_f
			run_mem += T_st[p]*TS*8;      // tnum_mean_st_f
			run_mem += C_st[p]*TS*8;      // cpop_st_f
			break;
		}
	}
	if(temp_mem > run_mem) run_mem == temp_mem;
	
	run_mem *= nchain;
	run_mem /= bytes_in_GB;

	let warn=[], sug=[];
	
	let error;
	if(state_mem > mem_state_sample_max){
		warn.push("Results may be hard to load into the interface due to insufficent memory to store state samples."); 
		sug = state_memory_error(type); 
	}
	
	if(param_mem > mem_param_sample_max){
		warn.push("Results may be hard to load into the interface due to insufficent memory to store parameter samples"); 
		sug = param_memory_error(type); 
	}
	
	if(run_mem > 0.5*RAM){
		warn.push("BICI may have difficulty running due to insufficent memory on this computer ("+precision(run_mem,2)+" GB are required and this computer has only "+RAM+" GB)."); 
		
		sug.push("Run the output file on a computer with more memory.");
		if(nchain > 1){
			sug.push("Reduce the number of chains run.");
		}
	}		

	if(warn.length > 0){
		inter.help = {title:"Limited computational resources", te:"", ok:"Run Anyway", sug:sug, file_list:ans.file_list, warn:warn};
		stop_loading_symbol();
		generate_screen();
		return false;
	}
	return true;
}


/// Generates a state memory warning
function state_memory_error(type)
{
	let sug=[];

	switch(type){
	case "sim":
		sug.push("The number of simulations could be decreased.");
		sug.push("The time-step could be increased.");
		break;
		
	case "inf":
		sug.push("Under 'Further options' the output number of posterior state samples could be decreased.");
		sug.push("The time-step could be increased.");
		break;
		
	case "ppc":
		sug.push("The number of posterior simulations could be decreased.");
		break;
	}

	return sug;
}


/// Generates a parameter memory warning
function param_memory_error(type)
{
	let sug=[];

	switch(type){
	case "sim":
		sug.push("The number of simulations could be decreased.");
		sug.push("Under 'Further options' the threshold number of tensor elements could be decreased.");
		break;
		
	case "inf":
		sug.push("Under 'Further options' the output number of posterior parameter samples could be decreased.");
		sug.push("Under 'Further options' the threshold number of tensor elements could be decreased.");
		break;
		
	case "ppc":
		sug.push("The number of posterior simulations could be decreased.");
		sug.push("Under 'Further options' the threshold number of tensor elements could be decreased.");
		break;
	}
	
	return sug;
}
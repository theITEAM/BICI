"use strict";

/// Starts a simulation
function start_simulation()
{
	start("sim",1);
}


/// Starts a simulation or inference 
function start(siminf,nchain)
{
	let stat;
	if(ver == "mac") stat = create_output_file(siminf,"init.txt"); 
	else stat = create_output_file(siminf,"Execute/init.txt"); 
		
	if(stat == "success"){
		inter.running_status = true;
		generate_screen();
		start_loading_symbol(1);
		
		inter.chain = [];
		
		if(siminf == "sim"){
			sim_result = clear_result();			
			results_add_model(sim_result,model.sim_details,siminf);
		}
		else{
			inf_result = clear_result();			
			results_add_model(inf_result,model.inf_details,siminf);
		}
		
		let nchain = 1; 
		if(siminf == "inf" && model.inf_details.algorithm.value =="DA-MCMC"){
			nchain = model.inf_details.nchain;
		}
		
		for(let ch = 0; ch < nchain; ch++){
			startspawn(ch,siminf);  
		}
	}
}


/// Clears any results
function clear_result()
{
	return {sample:[], par_sample:[]};
}


/// Starts execution of C++ code
function startspawn(ch,siminf)                                    
{		
	let num = Math.floor(Math.random()*1000); 
	
	inter.chain[ch] = { done:false, siminf:siminf, term:false, prog:0, leftover:"", lines:[]};
	
	switch(ver){
	case "windows": inter.child[ch] = spawn('Execute/bici.exe',["Execute/init.txt",num]); break;
	case "mac": inter.child[ch] = spawn('Execute/bici.out',["/tmp/init.txt",num]); break;
	}
	
	funct(inter.child[ch],ch);
}


/// This function executes the c++ code
function funct(chi,ch)                                      // Gathers output of C++ file
{
	chi.stdout.on('data', function (data) {
		//pr(data+" data");
		let cha = inter.chain[ch];
		let st = cha.leftover + data;
		
		let lines = st.split('\n');
		cha.leftover = lines[lines.length-1];
		
		for(let li = 0; li < lines.length-1; li++){
			let line = lines[li];
	
			if(line.length > 10 && line.substr(0,10) == "<PROGRESS>"){
				cha.prog = Number(line.substr(10));
				let progmin = LARGE;
				for(let ch2 = 0; ch2 < inter.chain.length; ch2++){
					if(inter.chain[ch2].prog < progmin) progmin = inter.chain[ch2].prog;
				}
				if(progmin > 0) set_loading_percent(progmin);
			}
			
			cha.lines.push(line);
		}
	});

	chi.stderr.on('data', function (data) {
		alertp("There was a problem running the code: "+endl+data); 
		terminate(); 
	});

	chi.on('close', function (code) {
		//pr("close");
		inter.chain[ch].done = true;
		
		let j = 0; 
		while(j < inter.chain.length && inter.chain[j].done == true && inter.chain[j].term == false) j++;
		if(j == inter.chain.length){
			setTimeout(function(){ process_all_chains()}, 10);
		}
		/*
		if(inter.chain[ch].term == false){
			setTimeout(function(){pr("ch"+ch); process_info(ch);}, 10);
		}
		*/
	});
}	


/// Processes information from all chains
function process_all_chains(ch)
{
	inter.loading_symbol.processing = true;
	inter.loading_symbol.percent = 1;
	
	replot_layer("LoadingSymbol");
	replot_loading_symbol(); plot_screen();

	setTimeout(function(){ process_chain(0)}, 10);
}


/// Processes a single chain	
function process_chain(ch)
{
	let siminf = inter.chain[ch].siminf;
	
	let result;
	switch(siminf){
	case "sim": result = sim_result; break;
	case "inf": result = inf_result; break;
	}
	
	process_info(ch,result);
	
	ch++;
	
	inter.loading_symbol.percent = Math.floor(100*ch/inter.chain.length);
	replot_loading_symbol(); plot_screen();
	
	if(ch < inter.chain.length){ setTimeout(function(){ process_chain(ch)}, 10);}
	else{		
		inter.running_status = false;
		stop_loading_symbol();

		intialise_plot_filters(result);
		initialise_pages();
	
		let newpage = "Simulation"; if(siminf == "inf") newpage = "Inference";
		change_page({pa:newpage, su:"Results"});
	}
}


/// Processes the output from the c++ code
function process_info(ch,result)
{
	let lines = inter.chain[ch].lines;
	
	for(let c = 0; c < lines.length; c++){ lines[c] = lines[c].replace(/\r/g, "");}

	// checks for any errors in the 
	let c = 0; 
	while(c < lines.length){
		
		inter.loading_symbol.percent = Math.floor(100*c/lines.length);
	
		let line = lines[c];
	
		if(line == "<<PARAMS>>" || line == "<<STATES>>" || line == "<<ERROR>>"){
			c++;
			let content="";
			while(c < lines.length && lines[c] != "<<END>>"){
				content += lines[c]+endl;
				c++;
			}
		
			if(c == lines.length){
				alertp("There was a problem running the code."); 
				terminate(); 
				return;
			}
			
			switch(line){
			case "<<PARAMS>>":
				{
					let warn = "Problem trasfering data";
					
					read_param_samples(ch+1,content,result,warn);
				}
				break;
			
			case "<<STATES>>":
				read_state_samples(ch+1,content,result);
				break;
				
			case "<<ERROR>>":
				{
					alertp("The following error was found:"+endl+content);
					terminate(); 
					return;
				}
				break;
			}
		}
	
		c++;
	}
}


/// Terminates the code
function terminate()
{
	if(inter.running_status == false) return;
	
	for(let ch = 0; ch < inter.chain.length; ch++){
		let cha = inter.child[ch];
		inter.chain[ch].term = true;
		inter.child[ch].stdin.pause();
		inter.child[ch].kill();
	}
	
	inter.running_status = false;
	stop_loading_symbol();
	generate_screen();
}


/// Adds buttons to start a simulation		
function add_sim_start_buts(lay)
{
	if(inter.running_status == true) return;
	
	let cx = corner.x;
	let cy = corner.y;
	
	if(inter.options == false){
		cy = lay.add_title("Start simulation",cx,cy,{te:start_sim_text});
		
		cy += 2;
		
		cy = lay.add_subtitle("Time range",cx,cy,WHITE,{te:time_range_text});
		cy = lay.add_paragraph("Set the time range over which simulation is performed:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
		let yy = cy-2.5;
		add_right_input_field(yy,"Start time",{type:"sim_t_start",update:true},lay);
		add_right_input_field(yy+4,"End time",{type:"sim_t_end",update:true},lay);

		cy += 6;
		
		cy = lay.add_subtitle("Simulation number",cx,cy,WHITE,{te:sim_num_text});
		cy = lay.add_paragraph("Set the number of simulations to be generated:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
		yy = cy-2.5;
		add_right_input_field(yy,"Number",{type:"sim_number",update:true},lay);
		
		cy += 2;
		
		cy = lay.add_subtitle("Time-step",cx,cy,WHITE,{te:sim_timestep});

		if(model.sim_details.algorithm == undefined) model.sim_details.algorithm = {value:"gillespie"};
		
		let xx = 33;
		
		cy = lay.add_paragraph("Set the time-step used to simulate",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);

		yy = cy-2.5;
		add_right_input_field(yy,"Time-step",{type:"sim_timestep",update:true},lay);
		
		add_corner_link("> Further options","Options",lay);
		
		lay.add_corner_button([["Start","Grey","StartSimulation"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	}
	else{
		let te = "Further simulation options";
		cy = lay.add_title(te,cx,cy,{te:further_sim_text});
		
		cy += 2;
		
		let gap = 2;
		
		cy = lay.add_subtitle("Individual max",cx,cy,WHITE,{te:indmax_text});
		
		cy = lay.add_paragraph("Maximum number of individuals (for individual-based models)",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
	
		let yy = cy-2.5;
		add_right_input_field(yy,"Maximum",{type:"sim_indmax",update:true},lay);
		
		cy += gap;
		
		/*
		cy = lay.add_subtitle("Seed",cx,cy,WHITE,{te:seed_text});
		
		cy = lay.add_paragraph("Set the random seed",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
	
		let yy = cy-2.5;
		add_right_input_field(yy,"Maximum",{type:"sim_indmax",update:true},lay);
		*/
		
		
		lay.add_corner_button([["Done","Grey","OptionsDone"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	}
}


function add_right_input_field(yy,te,op,lay)
{
	let dx = 8;
	let xx = lay.inner_dx-dx-3;
	add_input_field(xx,yy,dx,te,op,lay);
}
	
function add_input_field(xx,yy,dx,te,op,lay)
{
	lay.add_button({te:te, x:xx, y:yy, dx:dx, dy:0.8, type:"InputBoxName", back_col:WHITE});
	yy += 1.2;
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
		!(filt_type == "for sim" && param[i].variety == "reparam")
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
					case "Se": case "Sp": case "trap_prob": case "comp_prob": case "derive_param":
						break;
					
					default:
						{
							let w = wright;
							
							if(add_view_button(par,w-4,y,i,lay,model)) w -= 4.5;
							
							if(par.variety == "dist"){
								lay.add_checkbox(w-4.5,y+0.4,"Sample","Sample",par.sim_sample,WHITE);
								w -= 4.5;
							}
							
							if(par.variety == "dist" && par.sim_sample.check == true){
								if(par.dep.length == 0 || par.prior_split_check.check == false){
									display_distribution(i,x,y,lay,true,false,w);
								}
								else{
									display_distribution_split(i,x,y,lay,true,false,"dist",w);
								}
							}
							else{
								if(par.type == "trans_bp" && par.auto_value != undefined){
									let si = 1;
									let ac, allow_edit = false, te = "Auto";
									if(par.auto_value == false){
										ac = "Auto";
										te = "Auto";
										allow_edit = true;
										lay.add_button({te:te, x:38.5, y:y+0.25, dx:3, dy:1, ac:ac, type:"ViewSmall", i:i});
									}
									
									//if(is_dist(par) == true) 
									
									display_constant(i,x,y,lay,w-3,allow_edit);
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

	
/// Displays a row allowing a variable constant to be set
function display_constant(i,x,y,lay,w,allow_edit,source)
{
	if(source == undefined) source = model;
	let par = source.param[i];
	
	lay.display_param(x-par.label_info.dx-0.7,y-0.1,par.label_info);

	let si = 1.5;
	lay.add_button({te:"=", x:x, y:y, dy:si, type:"Text", font:get_font(si), si:si, col:BLACK});
	
	let te;
	if(par.dep.length == 0){
		te = par.value;
		if(allow_edit == false && isNaN(te)) te = "Auto calculate"; 
	}
	else{
		if(par.set == false){
			te = "Unset ";
			switch(par.dep.length){
			case 1: 
				te += "Vector("+par.list[0].length+")";
				break;

			case 2:
				te += "Matrix("+par.list[0].length+"×"+par.list[1].length+")";
				break;

			default:
				te += "Tensor("+par.list[0].length;
				for(let k = 1; k < par.dep.length; k++){
					te += "×"+par.list[k].length;
				}
				te += ")";
				break;
			}
		}
		else{
			te = JSON.stringify(par.value);
			te = te.replace(/\"/g," ");
		}
	}
	
	let fo = get_font(1.1,"","times");
		
	let ac; if(allow_edit != false) ac = "EditSimValue";
	
	lay.add_button({te:te, x:x+1.6, y:y+0., dx:w-x-1.6, dy:1.6, type:"ParamSimElement", source:source, font:fo, ac:ac, i:i, name:par.name, label_info:par.label_info});
}


/// Displays a row allowing a variable constant to be set
function display_reparam(i,x,y,lay,w)
{
	let par = model.param[i];
	
	lay.display_param(x-par.label_info.dx-0.7,y-0.1,par.label_info);

	let si = 1.5;
	lay.add_button({te:"=", x:x, y:y, dy:si, type:"Text", font:get_font(si), si:si, col:BLACK});
	
	let te, ac;
	if(par.dep.length == 0){
		te = par.value;
	}
	else{
		if(par.set == false){
			te = "Unset ";
			switch(par.dep.length){
			case 1: 
				te += "Vector("+par.list[0].length+")";
				break;

			case 2:
				te += "Matrix("+par.list[0].length+"×"+par.list[1].length+")";
				break;

			default:
				te += "Tensor("+par.list[0].length;
				for(let k = 1; k < par.dep.length; k++){
					te += "×"+par.list[k].length;
				}
				te += ")";
				break;
			}
		}
		else{
			te = JSON.stringify(par.value);
			te = te.replace(/\"/g," ");
		}
	}
	
	let fo = get_font(1.1,"","times");
	
	lay.add_button({te:te, x:x+1.6, y:y+0., dx:w-x-1.6, dy:1.6, type:"ReparamElement", font:fo, ac:"EditReparamValue", i:i, name:par.name, label_info:par.label_info});
}


/// Click to edit simulation value
function edit_sim_value(th,lay_name,i,source)
{
	let par = source.param[th];

	if(par.dep.length == 0){
		select_bubble(lay_name,i,{});
	}
	else{
		inter.edit_param = {type:"Value", source:source, value:copy(par.value), label_info:par.label_info, i:th};	
	}
}


/// Click to edit reparameterised expressions
function edit_reparam_value(th,lay_name,i)
{
	let par = model.param[th];

	if(par.dep.length == 0){
		select_bubble_over();
		inter.bubble.th = th;
	}
	else{
		inter.edit_param = {type:"Reparam", value:copy(par.value), label_info:par.label_info, i:th};	
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
	if(par.variety == "normal" || par.variety == "const"){
		if(par.spline.on == true){
			lay.add_button({te:"View", source:source, x:x, y:y+0.3, dx:3.7, dy:1, ac:"ViewSpline", type:"CombineIE", i:i});
			return true;
		}
		else{
			if(par.dep.length == 1){
				lay.add_button({te:"View", source:source, x:x, y:y+0.3, dx:3.7, dy:1, ac:"ViewVector", type:"CombineIE", i:i});
				return true;
			}
			
			if(par.dep.length == 2){
				lay.add_button({te:"View", source:source, x:x, y:y+0.3, dx:3.7, dy:1, ac:"ViewMatrix", type:"CombineIE", i:i});
				return true;
			}
		}
	}
	
	return false;
}


/// Adds a distance checkbox for a matrix
function add_distance_button(par,x,y,lay)
{
	if(par.dep.length == 2){
		let index = remove_prime(par.dep[0]);
		if(index == remove_prime(par.dep[1])){
			lay.add_checkbox(x+0.5,y+0.3,par,"Distance",par.dist_matrix,WHITE,{title:"Distance matrix", te:distmat_text});	
			return true;
		}
	}
	
	return false;
}
	

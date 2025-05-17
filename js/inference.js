"use strict";
// Functions related to inference

/// Starts inference on a linux cluster
function run_cluster()
{
	let siminf;
	switch(tab_name()){
	case "Simulation": siminf = "sim"; break;
	case "Inference": siminf = "inf"; break;
	case "Post. Simulation": siminf = "ppc"; break;
	}
	
	let ncore;
	if(siminf == "inf"){
		switch(model.inf_details.algorithm.value){
		case "DA-MCMC":
			ncore = Number(model.inf_details.nchain)/Number(model.inf_details.cha_per_core);
			break;
		case "PAS-MCMC":
			ncore = Number(model.inf_details.npart)/Number(model.inf_details.part_per_core);
			break;
		default: break;
		}
	}
	
	let line1 = "module load mpi/openmpi-x86_64";
	
	let line2="";
	if(ncore != 1){
		line2 += "mpirun -n ";
		if(ncore != undefined) line2 += ncore+" "; else line2 += "[cores] ";
	}		
	line2 += "./bici-para ";
	if(true) line2 += "[file.bici]";
	else line2 += "Execute/init.bici";

	switch(siminf){
	case "sim": line2 += " sim"; break;
	case "inf": line2 += " inf"; break;
	case "ppc": line2 += " post-sim"; break;
	}
	
	let te = "To run BICI on a Linux cluster the following steps must be followed:\n• <b>Create BICI file</b> – Click on the 'Save' button below to create and save the initialisation BICI file. Copy this to the cluster where you want BICI to run.\n• <b>Executable</b> – The executable 'bici-para' must also be copied from the BICI main folder to the cluster (<i>e.g.</i> this could be in the same directory as the BICI file).\n";
	
	//te += "• <b>Load MPI</b> – MPI will need to be installed on the cluster. Once installed , it can be loaded using:\n]>'<b>"+line1+"</b>' ["+butte+",'CopyText|BB"+line1+"'] \n";

	let butte = "'⟨⟨COPY⟩⟩'"; if(line2 == inter.copied) butte = "'⟨⟨COPIED⟩⟩'";
	te += "• <b>Run</b> – BICI is run using:\n]><b>'"+line2+"</b>' ["+butte+",'CopyText|BB"+line2+"']\n>>(where ";
	if(ncore == undefined) te += "'[core]' is the number of CPU cores and ";
	te += "'[file.bici]' is replaced by the name of the BICI file you created). ";

	butte = "'⟨⟨COPY⟩⟩'"; if(line1 == inter.copied) butte = "'⟨⟨COPIED⟩⟩'";
	te += "Note, if the error “mpirun: command not found...” is encountered, MPI can be loaded with the terminal command:\n]>'<b>"+line1+"</b>' ["+butte+",'CopyText|BB"+line1+"']\n• <b>Visualise</b> – Once BICI has run it puts its results into the script file. Copy this back to your local computer and load into the interface.";
	
	inter.help = { title:"Run BICI on Linux cluster", te:te, siminf:siminf, save:"StartCluster"};  
	inter.copied = undefined;	
}			


/// Checks if some features have been disabled
function check_disable()
{
	let alg = model.inf_details.algorithm.value;
	switch(alg){
	case "ABC": case "ABC-SMC":
		for(let p = 0; p < model.species.length; p++){
			let so = model.species[p].inf_source;
			for(let j = 0; j < so.length; j++){
				switch(so[j].type){
				case "Compartment": case "Transition": 			
				case "Diag. Test":
					alertp("Individual-level data cannot be handled by the "+alg+" algorithm");
					return true;
				}
			}
		}
		break;
	}
	
	switch(alg){
	case "PMCMC":
		disactivated("Inferences with the PMCMC algorithm cannot be performed.");
		return true;
	
	case "ABC-MBP":
		disactivated("Inferences with the ABC-MBP algorithm cannot be performed.");
		return true;
	
	//case "PAS-MCMC":
		//disactivated("Inferences with the PAS-MBP algorithm cannot be performed.");
		//return true;
	
	case "Emulation":
		disactivated("Inferences with the emulation algorithm cannot be performed.");
		return true;
	}
	
	return false;
}


/// Checks that the initial start and end times and time-step are all consistent
function check_time_error()
{	
	let de;
	let pag = inter.page[inter.pa];
	switch(pag.name){
	case "Simulation": de = model.sim_details; break;
	case "Inference": de = model.inf_details; break;
	case "Post. Simulation": de = model.ppc_details; break;
	default: error("option not recognised2"); return;
	}
		
	let t_start = Number(de.t_start);
	let t_end = Number(de.t_end);
	
	if(de.ppc_t_start) t_start = Number(de.ppc_t_start);
	if(de.ppc_t_end) t_end = Number(de.ppc_t_end);
	
	let dt = Number(de.timestep);

	if((t_end - t_start)/dt > TI_DIV_MAX){
		alertp("The number of time divisions must be fewer than "+TI_DIV_MAX);
		return true;
	}
	
	if(t_start >= t_end){
		alertp("The start time must be before the end time.");
		return true;
	}
	
	let num = (t_end-t_start)/dt;
	if(dif(num,Math.round(num))){
		alertp("The difference between the start '"+t_start+"' and end '"+t_end+"' times is not a multiple of the time-step '"+dt+"'.");
		return true;
	}
	
	// Check spline times
	for(let th = 0; th < model.param.length; th++){
		let par = model.param[th];
		if(par.spline.on){
			let gtimes = get_times(par.spline.knot,de);
			
			if(gtimes.err){
				return add_knot_warning(gtimes.msg,par);
			}
			
			let times = gtimes.times;
			
			if(times[0] != t_start){
				return add_knot_warning("The knot start time must be 'start' instead of '"+times[0]+"'.",par);
			}
			
			if(times[times.length-1] != t_end){
				return add_knot_warning("The knot end time must be 'end' instead of '"+times[times.length-1]+"'.",par);
			}
			
			for(let i = 1; i < times.length-1; i++){
				if(times[i] < t_start || times[i] > t_end){
					return add_knot_warning("The knot time '"+times[i]+"' is outside the time range.",par);
				}
				
				if(times[i] == t_start){
					return add_knot_warning("The start knot time '"+times[i]+"' is set twice.",par);
				}
				
				if(times[i] == t_end){
					return add_knot_warning("The end knot time '"+times[i]+"' is set twice.",par);
				}
			}
			
			for(let i = 0; i < times.length-1; i++){
				if(times[i] > times[i+1]){
					return add_knot_warning("The knots must be time ordered.",par);
				}
			}				
		}
	}
	
	return false;
}


/// Warning for knot timings
function add_knot_warning(te,par)
{
	add_warning({mess:"Problem with knot times", mess2:"For parameter '"+par.full_name+"': "+te, par:par, warn_type:"KnotTimeProblem"});
	return true;
}


/// Checks that the timestep fits neatly into 
function check_timestep_fit(t,de)
{
	let t_start = Number(de.t_start);
	let dt = Number(de.timestep);
	let num = (t-t_start)/dt;
	if(dif(num,Math.round(num))) return false;
	
	return true;
}


/// Sets buttons used to start inference
function add_inf_start_buts(lay)
{
	//if(inter.running_status == true) return;
	
	let cx = corner.x;
	let cy = corner.y;
	
	if(inter.options == false){
		cy = lay.add_title("Start inference",cx,cy,{te:start_inf_text});
		
		cy += 1;
		
		cy = lay.add_subtitle("Time range",cx,cy,WHITE,{te:inf_time_range_text});
		cy = lay.add_paragraph("Set time range over which inference is performed:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
		let yy = cy-2.5;
		add_right_input_field(yy,"Start time",{type:"inf_t_start",update:true},lay);
		add_right_input_field(yy+3.5,"End time",{type:"inf_t_end",update:true},lay);

		cy += 4;
			
		cy = lay.add_subtitle("Time-step",cx,cy,WHITE,{te:sim_timestep});

		cy = lay.add_paragraph("Set the time-step:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
		yy = cy-2.5;
		add_right_input_field(yy,"Time-step",{type:"inf_timestep",update:true},lay);

		cy += 1;
		
		cy = lay.add_subtitle("Algorithm",cx,cy,WHITE,{te:inf_algorithm});

		//cy = lay.add_paragraph("Set the algorithm used to perform inference:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
			
		//cy += 0.5;
		
		let xx = 3, gap = 0.5;
		xx = lay.add_radio(xx,cy+0.1,"DA-MCMC","DA-MCMC",model.inf_details.algorithm,{back_col:WHITE});
		xx += gap;
	
		xx = lay.add_radio(xx,cy+0.1,"PAS-MCMC","PAS-MCMC",model.inf_details.algorithm,{back_col:WHITE});
		xx += gap;
		
		xx = lay.add_radio(xx,cy+0.1,"ABC","ABC",model.inf_details.algorithm,{back_col:WHITE});
		xx += gap;
			
		xx = lay.add_radio(xx,cy+0.1,"ABC-SMC","ABC-SMC",model.inf_details.algorithm,{back_col:WHITE});
		xx += gap;
			
		xx = lay.add_radio(xx,cy+0.1,"PMCMC","PMCMC",model.inf_details.algorithm,{back_col:WHITE});
		xx += gap;
		
		//xx = lay.add_radio(xx,cy+0.1,"ABC-MBP","ABC-MBP",model.inf_details.algorithm,{back_col:WHITE});
		//xx += gap;
		
		//xx = lay.add_radio(xx,cy+0.1,"HMC","HMC",model.inf_details.algorithm,{back_col:WHITE});
		
		xx = lay.add_radio(xx,cy+0.1,"Emulation","Emulation",model.inf_details.algorithm,{back_col:WHITE});

		cy += 2.5;

		let gapop = 1;
			
		switch(model.inf_details.algorithm.value){
		case "DA-MCMC":
			{
				cy = lay.add_subtitle("Updates",cx,cy,WHITE,{te:sample_text});

				cy = lay.add_paragraph("Set the number of MCMC updates:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
				yy = cy-2.5;
				add_right_input_field(yy,"Updates",{type:"inf_sample",update:true},lay);
			}
			cy += gapop;	

			{			
				cy = lay.add_subtitle("Chains",cx,cy,WHITE,{te:sample_text});
			
				cy = lay.add_paragraph("Set the number of MCMC chains:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
				yy = cy-2.5;
				add_right_input_field(yy,"Chains",{type:"inf_chain",update:true},lay);
			}
			break;
		
		case "PAS-MCMC":
			{
				cy = lay.add_subtitle("Updates",cx,cy,WHITE,{te:pas_sample_text});

				cy = lay.add_paragraph("Set the number of MCMC updates:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
				yy = cy-2.5;
				add_right_input_field(yy,"Updates",{type:"inf_sample",update:true},lay);
			}
			
			cy += gapop;	
			
			{
				cy = lay.add_subtitle("Particles",cx,cy,WHITE,{te:particle_text});
			
				cy = lay.add_paragraph("Set the number of particles:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
				yy = cy-2.5;
				add_right_input_field(yy,"Particles",{type:"inf_npart",update:true},lay);
			}
			
			cy += gapop;	

			{
				cy = lay.add_subtitle("Updates per generation",cx,cy,WHITE,{te:gen_update_text});
				
				cy = lay.add_paragraph("Set the number of updates per generation as a percentage of MCMC updates:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
				yy = cy-2.5;
				add_right_input_field(yy,"Percentage",{type:"inf_gen_update",update:true},lay);
			}
			break;
			
		case "ABC":
			{
				cy = lay.add_subtitle("Samples",cx,cy,WHITE,{te:abc_sample_text});

				cy = lay.add_paragraph("Set the number of samples:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
				yy = cy-2.5;
				add_right_input_field(yy,"Samples",{type:"inf_abcsample",update:true},lay);
			}
			
			cy += gapop;
			
			{
				cy = lay.add_subtitle("Acceptance",cx,cy,WHITE,{te:acceptance_text});

				cy = lay.add_paragraph("Set the acceptance fraction:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
				yy = cy-2.5;
				add_right_input_field(yy,"Fraction",{type:"inf_accfrac",update:true},lay);
			}
			break;
			
		case "ABC-SMC":
			{
				cy = lay.add_subtitle("Samples",cx,cy,WHITE,{te:abc_sample_text});

				cy = lay.add_paragraph("Set the number of samples:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
				yy = cy-2.5;
				add_right_input_field(yy,"Samples",{type:"inf_abcsample",update:true},lay);
			}
			
			cy += gapop;

			{
				cy = lay.add_subtitle("Generations",cx,cy,WHITE,{te:abc_gen_text});

				cy = lay.add_paragraph("Set the number of generations:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
				yy = cy-2.5;
				add_right_input_field(yy,"Generations",{type:"inf_numgen",update:true},lay);
			}
			break;
		}

		cy += gapop;	
		
		cy = run_local(cx,cy,model.inf_details,lay);
	
		cy += gapop;	
		
		add_corner_link("> Further options","Options",lay);
		
		lay.add_corner_button([["Start","Grey","StartInference"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	}
	else{
		let alg = model.inf_details.algorithm.value;
		let te = "Further options for the "+alg+" algorithm";
		
		cy = lay.add_title(te,cx,cy,{te:further_inf_text});
		
		cy += 1;
		
		let gap = 1.8;
		if(alg =="DA-MCMC") gap = 0.9;
		if(alg =="PAS-MCMC") gap = 1.4;
		
		switch(alg){
		case "DA-MCMC": case "PAS-MCMC": 
			{
				cy = lay.add_subtitle("Output parameter sample number",cx,cy,WHITE,{te:param_sample_text});
				
				cy = lay.add_paragraph("Set the number of parameter samples to be generated:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
			
				let yy = cy-2.5;
				add_right_input_field(yy,"Param. samples output",{type:"inf_output_param ",update:true},lay);

				cy += gap;
			
				cy = lay.add_subtitle("Output state sample number",cx,cy,WHITE,{te:state_sample_text});
				
				cy = lay.add_paragraph("Set the number of state samples to be generated:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
			
				yy = cy-2.5;
				add_right_input_field(yy,"State samples output",{type:"inf_outputstate",update:true},lay);
				
				cy += gap;
				
				if(alg == "DA-MCMC"){
					cy = lay.add_subtitle("Chains per core",cx,cy,WHITE,{te:chain_per_core_text});
				
					cy = lay.add_paragraph("Set the number of chains to be run on each CPU core:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
			
					yy = cy-2.5;
					add_right_input_field(yy,"Chains per core",{type:"inf_cha_per_core",update:true},lay);
				}
				
				if(alg == "PAS-MCMC"){
					cy = lay.add_subtitle("Particles per core",cx,cy,WHITE,{te:part_per_core_text});
				
					cy = lay.add_paragraph("Set the number of particles to be run on each CPU core:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
			
					yy = cy-2.5;
					add_right_input_field(yy,"Particles per core",{type:"inf_part_per_core",update:true},lay);
				}
				
				cy += gap;
			}
			break;
			
		case "ABC-SMC":
			{
				cy = lay.add_subtitle("Acceptance",cx,cy,WHITE,{te:abcsmc_acceptance_text});
			
				cy = lay.add_paragraph("Set the acceptance fraction for each generation:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
			
				let yy = cy-2.5;
				add_right_input_field(yy,"Fraction",{type:"inf_accfracsmc",update:true},lay);
		
				cy += gap;
	
				cy = lay.add_subtitle("Kernel",cx,cy,WHITE,{te:kernel_text});
			
				cy = lay.add_paragraph("Set size of proposal kernel relative to the posterior:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
			
				yy = cy-2.5;
				add_right_input_field(yy,"Kernel size",{type:"inf_kernelsize",update:true},lay);
				
				cy += gap;
			}
			break;
		}
		
		cy = lay.add_subtitle("Individual maximum",cx,cy,WHITE,{te:indmax_text});
				
		cy = lay.add_paragraph("Maximum number of individuals (for individual-based models):",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
	
		let yy = cy-2.5;
		add_right_input_field(yy,"Maximum",{type:"inf_indmax",update:true},lay);
	
		cy += gap;
		
		cy = lay.add_subtitle("Parameter outputs",cx,cy,WHITE,{te:paramoutputmax_text});
				
		cy = lay.add_paragraph("Threshold number of tensor elements above which tensor not output:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
	
		yy = cy-2.5;
		add_right_input_field(yy,"Maximum",{type:"inf_paramout",update:true},lay);
	
		cy += gap;
	
		switch(alg){
		case "DA-MCMC":
			{
				cy = lay.add_subtitle("Burn-in",cx,cy,WHITE,{te:burnin_text});
				
				cy = lay.add_paragraph("Annealing used during the burn-in phase:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
				
				let pos = [{te:"none"},{te:"scan"},{te:"power-auto"},{te:"log-auto"},{te:"power"}];
				
				let yy = cy-2.5;
				add_right_dropdown(yy,"Type",model.inf_details.anneal_type,pos,lay);
				
				let type = model.inf_details.anneal_type.te;
				
				switch(type){
				case "scan":
					{
						//cy += 1.5;
						
						//cy = lay.add_paragraph("Rate of annealing",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
			
						let yy = cy-2.5;
						add_right_input_field(yy,"Rate",{type:"anneal_rate",update:true},lay,9);
					}
					break;
					
				case "power":
					{
						//cy += 1.5;
						
						//cy = lay.add_paragraph("Power in annealing model",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
			
						let yy = cy-2.5;
						add_right_input_field(yy,"Power",{type:"anneal_power",update:true},lay,9);
					}
					break;
				}
				
				if(type != "scan"){
					cy += 1.5;
						
					cy = lay.add_paragraph("Set the number of updates for burn-in as a percentage of MCMC updates:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
					let yy = cy-2.5;
					add_right_input_field(yy,"Percentage",{type:"burnin_frac",update:true},lay);
				}
				
				cy += gap;
			}
			break;
			
		case "PAS-MCMC":
			{
				cy = lay.add_subtitle("Burn-in",cx,cy,WHITE,{te:burnin_text});
				
				cy = lay.add_paragraph("Set the number of updates for burn-in as a percentage of MCMC updates:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
	
				let yy = cy-2.5;
				add_right_input_field(yy,"Percentage",{type:"burnin_frac",update:true},lay);
				cy += gap;
			}	
			break;
		}
		
		cy = set_seed(cx,cy,"inf_seed",model.inf_details,lay);
		
		lay.add_corner_button([["Done","Grey","OptionsDone"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	}
}


/// Adds a link in the corner of the page
function add_corner_link(te,ac,lay)
{
	let si = 0.8;
	let fo = get_font(si);
	let w = text_width(te,fo);
	lay.add_button({te:te, x:corner.x, y:lay.dy-2, dx:w+0.2, dy:1, type:"Link", ac:ac, si:si, font:fo});
}

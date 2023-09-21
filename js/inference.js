"use strict";

/// Starts inference
function start_inference()
{
	start("inf",1);
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
				case "Compartment": case "Transition": case "Source": case "Sink": case "Diag. Test":
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
	
	case "PAS-MBP":
		disactivated("Inferences with the PAS-MBP algorithm cannot be performed.");
		return true;
	
	case "Emulation":
		disactivated("Inferences with the emulation algorithm cannot be performed.");
		return true;
	}
	
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p]; 
		for(let cl = 0; cl < sp.ncla; cl++){
			let claa = sp.cla[cl];
			for(let tr = 0; tr < claa.tra.length; tr++){
				let tra = claa.tra[tr];
				if(tra.variety == "Source"){
					disactivated("Inferences with source transitions cannot be performed.");
					return true;
				}
				
				if(tra.variety == "Sink"){
					disactivated("Inferences with sink transitions cannot be performed.");
					return true;
				}
			
				if(tra.type == "gamma"){
					disactivated("Inferences with gamma distributed transitions cannot be performed.");
					return true;
				}
				
				if(tra.type == "log-normal"){
					disactivated("Inferences with log-normally distributed transitions cannot be performed.");
					return true;
				}
				
				if(tra.type == "weibull"){
					disactivated("Inferences with Weibull distributed transitions cannot be performed.");
					return true;
				}
			}
		}
	}
	
	return false;
}


/// Checks that the intial start and end times and time-step are all consistent
function check_time_error()
{
	let de;
	switch(inter.page[inter.pa].name){
	case "Simulation": de = model.sim_details; break;
	case "Inference": de = model.inf_details; break;
	default: error("option not recognised"); return;
	}
		
	let t_start = Number(de.t_start);
	let t_end = Number(de.t_end);
	let dt = Number(de.timestep);

	if(t_start >= t_end){
		alertp("The start time must be before the end time.");
		return true;
	}
	
	let num = (t_end-t_start)/dt;
	if(dif(num,Math.round(num))){
		alertp("The timestep must fit exactly as an integer quantity between the start and end times.");
		return true;
	}
		
	if(inter.page[inter.pa].name == "Inference"){
		for(let p = 0; p < model.species.length; p++){
			let sp = model.species[p];
			let so = sp.inf_source;
			let spec = so.spec;
			
			switch(so.type){
			case "Transition":
				// TO DO
				break;
			}
		}	
	}
	
	return false;
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
	if(inter.running_status == true) return;
	
	let cx = corner.x;
	let cy = corner.y;
	
	if(inter.options == false){
		cy = lay.add_title("Start inference",cx,cy,{te:start_inf_text});
		
		cy += 1;
		
		cy = lay.add_subtitle("Time range",cx,cy,WHITE,{te:inf_time_range_text});
		cy = lay.add_paragraph("Set time range over which inference is performed:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
		let yy = cy-2.5;
		add_right_input_field(yy,"Start time",{type:"inf_t_start",update:true},lay);
		add_right_input_field(yy+4,"End time",{type:"inf_t_end",update:true},lay);

		cy += 5;
			
		cy = lay.add_subtitle("Time-step",cx,cy,WHITE,{te:sim_timestep});

		cy = lay.add_paragraph("Set the time-step:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
		yy = cy-2.5;
		add_right_input_field(yy,"Time-step",{type:"inf_timestep",update:true},lay);

		cy += 2;
		
		cy = lay.add_subtitle("Algorithm",cx,cy,WHITE,{te:inf_algorithm});

		cy = lay.add_paragraph("Set the algorithm used to perform inference:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
			
		cy += 1;
		
		let xx = 3, gap = 0.5;
		xx = lay.add_radio(xx,cy+0.1,"DA-MCMC","DA-MCMC",model.inf_details.algorithm,{back_col:WHITE});
		xx += gap;
		
		xx = lay.add_radio(xx,cy+0.1,"ABC","ABC",model.inf_details.algorithm,{back_col:WHITE});
		xx += gap;
			
		xx = lay.add_radio(xx,cy+0.1,"ABC-SMC","ABC-SMC",model.inf_details.algorithm,{back_col:WHITE});
		xx += gap;
			
		xx = lay.add_radio(xx,cy+0.1,"PMCMC","PMCMC",model.inf_details.algorithm,{back_col:WHITE});
		xx += gap;
		
		xx = lay.add_radio(xx,cy+0.1,"ABC-MBP","ABC-MBP",model.inf_details.algorithm,{back_col:WHITE});
		xx += gap;
		
		xx = lay.add_radio(xx,cy+0.1,"PAS-MBP","PAS-MBP",model.inf_details.algorithm,{back_col:WHITE});
		xx += gap;
			
		//xx = lay.add_radio(xx,cy+0.1,"HMC","HMC",model.inf_details.algorithm,{back_col:WHITE});
		
		xx = lay.add_radio(xx,cy+0.1,"Emulation","Emulation",model.inf_details.algorithm,{back_col:WHITE});

		cy += 3;

		let gapop = 1;
			
		switch(model.inf_details.algorithm.value){
		case "DA-MCMC":
			cy = lay.add_subtitle("Samples",cx,cy,WHITE,{te:sample_text});

			cy = lay.add_paragraph("Set the number of samples:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
			yy = cy-2.5;
			add_right_input_field(yy,"Samples",{type:"inf_sample",update:true},lay);

			cy += gapop;	
			
			cy = lay.add_subtitle("Chains",cx,cy,WHITE,{te:sample_text});
			
			cy = lay.add_paragraph("Set the number of MCMC chains:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
			yy = cy-2.5;
			add_right_input_field(yy,"Chains",{type:"inf_chain",update:true},lay);

			cy += gapop;	
			break;
		
		case "ABC":
			cy = lay.add_subtitle("Samples",cx,cy,WHITE,{te:abc_sample_text});

			cy = lay.add_paragraph("Set the number of samples:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
			yy = cy-2.5;
			add_right_input_field(yy,"Samples",{type:"inf_abcsample",update:true},lay);

			cy += gapop;
			
			cy = lay.add_subtitle("Acceptance",cx,cy,WHITE,{te:acceptance_text});

			cy = lay.add_paragraph("Set the acceptance fraction:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
			yy = cy-2.5;
			add_right_input_field(yy,"Fraction",{type:"inf_accfrac",update:true},lay);
			break;
			
		case "ABC-SMC":
			cy = lay.add_subtitle("Samples",cx,cy,WHITE,{te:abc_sample_text});

			cy = lay.add_paragraph("Set the number of samples:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
			yy = cy-2.5;
			add_right_input_field(yy,"Samples",{type:"inf_abcsample",update:true},lay);

			cy += gapop;

			cy = lay.add_subtitle("Generations",cx,cy,WHITE,{te:abc_gen_text});

			cy = lay.add_paragraph("Set the number of generations:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
		
			yy = cy-2.5;
			add_right_input_field(yy,"Generations",{type:"inf_numgen",update:true},lay);
			
			cy += gapop;
			break;
		}

		add_corner_link("> Further options","Options",lay);
		
		lay.add_corner_button([["Start","Grey","StartInference"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	}
	else{
		let alg = model.inf_details.algorithm.value;
		let te = "Further options for the "+alg+" algorithm";
		
		cy = lay.add_title(te,cx,cy,{te:further_inf_text});
		
		cy += 2;
		
		let gap = 2;
		
		switch(model.inf_details.algorithm.value){
		case "DA-MCMC": 
			{
				cy = lay.add_subtitle("Parameter thinning",cx,cy,WHITE,{te:param_thin_text});
				
				cy = lay.add_paragraph("Set the parameter thinning factor",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
			
				let yy = cy-2.5;
				add_right_input_field(yy,"Parameter thinning",{type:"inf_thinparam",update:true},lay);

				cy += gap;
			
				cy = lay.add_subtitle("State thinning",cx,cy,WHITE,{te:state_thin_text});
				
				cy = lay.add_paragraph("Set the state thinning factor",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
			
				yy = cy-2.5;
				add_right_input_field(yy,"State thinning",{type:"inf_thinstate",update:true},lay);
				
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
			
				cy = lay.add_paragraph("Set size of proposal kernel relative to posterior:",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
			
				yy = cy-2.5;
				add_right_input_field(yy,"Kernel size",{type:"inf_kernelsize",update:true},lay);
		
				cy += gap;
			}
			break;
		}
		
		cy = lay.add_subtitle("Individual maximum",cx,cy,WHITE,{te:indmax_text});
				
		cy = lay.add_paragraph("Maximum number of individuals (for individual-based models)",lay.inner_dx-2*cx,cx,cy,BLACK,para_si,para_lh);
	
		let yy = cy-2.5;
		add_right_input_field(yy,"Maximum",{type:"inf_indmax",update:true},lay);
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

"use strict";
// Functions which simulated data

/// Simulate data
function sim_data()
{
	start_worker("Sim data",edit_source);	
}


/// Samples from the observation model
function obs_model_sample(mean,so)
{
	if(mean < 0) mean = 0;
	
	switch(so.spec.obs_error_dist.value){
	case "Normal":
		{
			let sd;
			switch(so.spec.obs_error.value){
			case "percent": sd = mean*Number(so.spec.percent)/100; break;
			case "sd": sd = Number(so.spec.sd); break;
			default: error("option not recognised"); break;
			}
			
			do{
				let res = normal_sample(mean,sd);
				if(res >= 0) return res;
			}while(true);
		}
		break;
	
	case "Poisson":
		return poisson_sample(mean);
		
	case "Negative binomial":
		return neg_binommial_sample(mean,Number(so.spec.p));
	}
}


/// Simulates data in the webworker
function simulate_data(so)
{	
	let p = so.info.p;

	let result = sim_result;
	
	let t_start = Number(result.details.t_start);
	let t_end = Number(result.details.t_end);
	let timestep = Number(result.details.timestep);
	
	let sp = result.species[p];

	let samp = result.sample[0];

	let tp = result.timepoint;
	
	let sim = samp.species[p];

	let transnum = reconstruct_timeline_vec(sim.transnum_tl,result);
	let cpop = reconstruct_timeline_vec(sim.cpop_tl,result);
	
	let spec = so.spec;
		
	let head = [];
	let ele = [];

	switch(so.type){
	case "Genetic":
		switch(spec.type_radio.value){
		case "snp": simulate_genetic_snp(t_start,t_end,ele,head,so); break;
		case "matrix": simulate_genetic_matrix(t_start,t_end,ele,head,so); break;
		}
		result.genetic_data = copy(ele);
		break;
		
	case "Init. Pop.":
		{
			for(let cl = 0; cl < sp.ncla; cl++) head.push(sp.cla[cl].name);
			head.push("Population");

			for(let c = 0; c < sp.comp_gl.length; c++){
				let cgl = sp.comp_gl[c];
					
				let row = [];
				for(let cl = 0; cl < sp.ncla; cl++){
					row.push(sp.cla[cl].comp[cgl.cla_comp[cl]].name);
				}
				row.push(sim.cpop_init[c]);
				ele.push(row);
			}
		}
		break;
	
	case "Add Pop.": case "Remove Pop.":
		{
			if(sim.add_rem_pop_change){
				let sign = 1; if(so.type == "Remove Pop.") sign = -1;
				let T = tp.length-1;
		
				head.push("t");
				for(let cl = 0; cl < sp.ncla; cl++) head.push(sp.cla[cl].name);
				head.push("Population");
				
				for(let ti = 0; ti < T; ti++){	
					if(sim.add_rem_pop_change[ti]){	
						let arpc = sim.add_rem_pop_change[ti];
						for(let k = 0; k < arpc.length; k++){
							let val = arpc[k].dpop*sign;
							if(val > 0){
								let c = arpc[k].c;
								let cgl = sp.comp_gl[c];
								
								let row = [];
								row.push(t_start+ti*timestep);
							
								for(let cl = 0; cl < sp.ncla; cl++){
									row.push(sp.cla[cl].comp[cgl.cla_comp[cl]].name);
								}
								
								row.push(val);
								ele.push(row);
							}
						}
					}
				}
			}
		}
		break;
	
	case "Add Ind.":
		{
			head.push("ID");
			head.push("t");
			for(let cl = 0; cl < sp.ncla; cl++) head.push(sp.cla[cl].name);
			
			for(let i = 0; i < sim.individual.length; i++){
				let ind = sim.individual[i];
				if(ind.obs){	
					let row = [];
					row.push(ind.name);
					
					let cinit = ind.cinit
					if(cinit != SOURCE){ 
						if(cinit == OUT){
							let ev = ind.ev[0]
							row.push(ev.t);
							cinit = ev.c;
						}
						else row.push(t_start);
						
						let cgl = sp.comp_gl[cinit];
						for(let cl = 0; cl < sp.ncla; cl++){
							row.push(sp.cla[cl].comp[cgl.cla_comp[cl]].name);
						}
						ele.push(row);
					}
				}
			}
		}
		break;
		
	case "Remove Ind.":
		{
			head.push("ID");
			head.push("t");
		
			for(let i = 0; i < sim.individual.length; i++){
				let ind = sim.individual[i];
				if(ind.obs){	
					let row = [];
					row.push(ind.name);
					
					if(ind.ev.length > 0){
						let ev = ind.ev[ind.ev.length-1];
						if(ev.type == EV_LEAVE){
							row.push(ev.t);	
							ele.push(row);
						}
					}
				}
			}
		}
		break;
		
	case "Move Ind.":
		{
			head.push("ID");
			head.push("t");
			let cl_drop = spec.cl_drop;
			head.push(cl_drop.te);
			let cl_sel = cl_drop.i;
			
			for(let i = 0; i < sim.individual.length; i++){
				let ind = sim.individual[i];
				if(ind.obs){	
					for(let e = 0; e < ind.ev.length; e++){
						let ev = ind.ev[e];
						if(ev.type == EV_MOVE && ev.cl == cl_sel){
							let row = [];
							row.push(ind.name);
							row.push(ev.t);
							let cgl = sp.comp_gl[ev.cf];
							row.push(sp.cla[cl_sel].comp[cgl.cla_comp[cl_sel]].name);
							ele.push(row);
						}
					}
				}
			}
		}
		break;
		
	case "Compartment":
		{
			let exact = false;
			let acc;
			if(so.exact.value == "Exact") exact = true;
			else acc = Number(so.comp_acc);
				
			let times = get_timepoints(t_start,t_end,so);
			if(typeof times == 'string'){ alert_help("Problem with times",times); return;}
	
			let f = so.frac_obs;
			
			head.push("ID");
			head.push("t");
			let cl_drop = spec.cl_drop;
			head.push(cl_drop.te);
			
			let cl_sel = find(sp.cla,"name",cl_drop.te);
			
			for(let i = 0; i < sim.individual.length; i++){
				let ind = sim.individual[i];
				if(ind.obs){	
					for(let k = 0; k < times.length; k++){
						if(Math.random() < f){	
							let t = times[k];
						
							let c = ind.cinit;
							let e = 0; 
							while(e < ind.ev.length && ind.ev[e].t < t){
								c = c_after_event(ind.ev[e],sp);
								e++;
							}
						
							let te;
						
							let row = [];	
							row.push(ind.name);
							row.push(t);
							
							if(c == SOURCE || c == SINK || c == OUT){
								te = "!";
							}
							else{
								let cgl = sp.comp_gl[c];
							
								let comp = sp.cla[cl_sel].comp;
								let ci = cgl.cla_comp[cl_sel];
							
								
								if(exact){
									te = comp[ci].name;
								}
								else{
									te = comp[ci].name+":"+acc;
									if(comp.length > COMP_NOISY_MAX){  // If many comparmtents then randomly chooses some 
										let list=[];
										for(let cii = 0; cii < comp.length; cii++){
											if(cii != ci) list.push(cii);
										}
										let ac_other = precision(((1-acc)/(COMP_NOISY_MAX-1)),3);
										for(let k = 0; k < COMP_NOISY_MAX-1; k++){
											let m = Math.floor(Math.random()*list.length);
											let cii = list[m];
											te += "|"+comp[cii].name+":"+ac_other;
											list[m] = list[list.length-1];
											list.pop();
										}
									}
									else{
										let ac_other = precision(((1-acc)/(comp.length-1)),3);
										for(let cii = 0; cii < comp.length; cii++){
											if(cii != ci){
												te += "|"+comp[cii].name+":"+ac_other;
											}
										}
									}
								}
							}
							row.push(te);
							ele.push(row);
						}
					}
				}
			}
		}
		break;
		
	case "Transition":
		{
			let tmin = t_start, tmax = t_end;
			if(spec.time_radio.value == "Spec"){
				
				tmin = Number(spec.time_start);
				tmax = Number(spec.time_end);
	
				if(tmin < t_start){ 
					alert_help("The start time '"+tmin+"' cannot be before the simulation start time '"+t_start+"'"); return;
				}
				
				if(tmax > t_end){ 
					alert_help("The end time '"+tmax+"' cannot be after the simulation end time '"+t_end+"'"); return;
				}
				
				if(tmin >= tmax){
					alert_help("The start time '"+tmin+"' must be before the end time '"+tmax+"'"); return;
				}					
			}
			
			let cl_drop = spec.cl_drop;
			let cl_sel = find(sp.cla,"name",cl_drop.te);
		
			let filt = copy(spec.filter);
		
			data_sim_setup_eqn(filt,sp,cl_sel);

			head.push("ID");
			head.push("t");
		
			for(let i = 0; i < sim.individual.length; i++){
				let ind = sim.individual[i];
				if(ind.obs){	
					for(let e = 0; e < ind.ev.length; e++){
						let ev = ind.ev[e];
						if(ev.type == EV_TRANS){
							if(ev.t > tmin && ev.t < tmax){
								let prob = get_tra_gl_prob(ev.trg,ev.t,filt,sp,cl_sel);

								if(Math.random() < prob){
									let row = [];
									row.push(ind.name);
									row.push(ev.t);
									ele.push(row);
								}
							}
						}
					}
				}
			}
		}
		break;
	
	case "Diag. Test":
		{
			let times = get_timepoints(t_start,t_end,so);
			if(typeof times == 'string'){ alert_help("Problem with times",times); return;}
	
			let cl_drop = spec.cl_drop;
			let cl_sel = find(sp.cla,"name",cl_drop.te);
			
			head.push("ID");
			head.push("t");
			head.push("Result");

			let Se = get_equation_value(spec.Se_eqn);
			let Sp = get_equation_value(spec.Sp_eqn);
			let neg = spec.neg_result;
			let pos = spec.pos_result;
		
			let f = so.frac_obs;
		
			for(let i = 0; i < sim.individual.length; i++){
				let ind = sim.individual[i];
				if(ind.obs){	
					for(let k = 0; k < times.length; k++){
						if(Math.random() < f){	
							let t = times[k];
						
							let c = ind.cinit;
							let e = 0; 
							while(e < ind.ev.length && ind.ev[e].t < t){
								c = c_after_event(ind.ev[e],sp);
								e++;
							}
						
							if(c != SOURCE && c != SINK && c != OUT){
								let cgl = sp.comp_gl[c];
								let ci = cgl.cla_comp[cl_sel];

								let res;
								if(spec.check_box.value[ci].check == true){
									if(Math.random() < Se) res = pos; else res = neg;
								}
								else{
									if(Math.random() < Sp) res = neg; else res = pos;
								}
								
								let row = [];
								row.push(ind.name);
								row.push(t);
								row.push(res);
								ele.push(row);
							}
						}
					}
				}
			}
		}
		break;
		
	case "Population":
		{
			let times = get_timepoints(t_start,t_end,so);
		
			if(typeof times == 'string'){ alert_help("Problem with times",times); return;}

			let filt = copy(spec.filter);
	
			let strat = get_stratify(filt,sp);
		
			data_sim_setup_eqn(filt,sp);
			
			head.push("t");
			for(let i = 0; i < strat.cl_list.length; i++) head.push(sp.cla[strat.cl_list[i]].name);
			head.push("Population");
			
			for(let com = 0; com < strat.comb.length; com++){
				set_filt_strat(com,strat,filt);
				
				for(let k = 0; k < times.length; k++){
					let t = times[k];
				
					let posg = get_comp_prob(t,filt,sp);
					for(let c = 0; c < sp.comp_gl.length; c++){
						if(posg[c] < 0){ alert_help("Problem with observation model","The observation model has become negative"); return;}
					}
				
					let popc = [];
					
					switch(sp.type){
					case "Individual":
						{
							for(let c = 0; c < sp.comp_gl.length; c++) popc[c] = 0;
						
							for(let i = 0; i < sim.individual.length; i++){
								let ind = sim.individual[i];
								if(ind.obs){	
									let c = ind.cinit;
									let e = 0; 
									while(e < ind.ev.length && ind.ev[e].t < t){
										c = c_after_event(ind.ev[e],sp);
										e++;
									}
									
									if(!isNaN(c)) popc[c]++;
								}
							}
						}
						break;
						
					case "Population":
						{
							let ti = Math.floor((t-t_start)/timestep);
							for(let c = 0; c < sp.comp_gl.length; c++) popc[c] = cpop[c][ti];
						}
						break;
					}
					
					let row = [];
					row.push(t);
					
					for(let i = 0; i < strat.cl_list.length; i++){
						row.push(strat.list[i][strat.comb[com].index[i]]);
					}
					
					let sum = 0;
					for(let c = 0; c < sp.comp_gl.length; c++){
						sum += posg[c]*popc[c];
					}
					
					let result = obs_model_sample(sum,so);

					row.push(precision(result,5));

					ele.push(row);
				}
			}
		}
		break;
		
	case "Pop. Trans.":
		{
			let times = get_timepoints(t_start,t_end,so);
			if(typeof times == 'string'){ alert_help("Problem with times",times); return;}
	
			if(times[0] != t_start) times.splice(0,0,t_start);
			if(times[times.length-1] != t_end) times.push(t_end);
		
			let ndiv = times.length-1;
			let div_ref = [];
			
			let k = 0;
			for(let j = 0; j < ndiv; j++){
				let t = times[j+1];
				while(k < tp.length && tp[k] < t-TINY){ div_ref[k] = j; k++;}		
				if(dif(t,tp[k])) alertp("Times must coincide with time points used to simulate.")
			}
			
			let cl_drop = spec.cl_drop;
			let cl_sel = find(sp.cla,"name",cl_drop.te);
		
			let filt = copy(spec.filter);
			
			let strat = get_stratify(filt,sp);
			
			data_sim_setup_eqn(filt,sp,cl_sel);

			head.push("Start");
			head.push("End");
			for(let i = 0; i < strat.cl_list.length; i++) head.push(sp.cla[strat.cl_list[i]].name);
			head.push("Number");
			
			for(let com = 0; com < strat.comb.length; com++){
				set_filt_strat(com,strat,filt);
				
				let res = [];
				for(let d = 0; d < ndiv; d++){
					res[d] = [];
					for(let tri = 0; tri < sp.tra_gl.length; tri++){
						res[d][tri] = 0;
					}
				}

				switch(sp.type){
				case "Individual":
					{
						for(let i = 0; i < sim.individual.length; i++){
							let ind = sim.individual[i];
							if(ind.obs){	
								for(let e = 0; e < ind.ev.length; e++){
									let ev = ind.ev[e];
									if(ev.type == EV_TRANS){
										let k = Math.floor((ev.t-t_start)/timestep);
										res[div_ref[k]][ev.trg]++;
									}
								}
							}
						}
					}
					break;
					
				case "Population":
					{
						for(let ti = 0; ti < tp.length-1; ti++){
							let d = div_ref[ti];
							for(let trg = 0; trg < sp.tra_gl.length; trg++){
								res[d][trg] += transnum[trg][ti];
							}
						}
					}
					break;
				}
					
				for(let d = 0; d < ndiv; d++){
					let row = [];
					row.push(times[d]);
					row.push(times[d+1]);

					for(let i = 0; i < strat.cl_list.length; i++){
						row.push(strat.list[i][strat.comb[com].index[i]]);
					}
					
					let t = (times[d]+times[d+1])/2;

					let posg = get_comp_prob(t,filt,sp,cl_sel);

					let tr_prob = get_tr_prob(t,filt,posg,sp,cl_sel);

					let sum = 0;
					for(let tri = 0; tri < sp.tra_gl.length; tri++){
						sum += tr_prob[tri]*res[d][tri];			
					}
					
					let result = obs_model_sample(sum,so);
					row.push(precision(result,5));
					ele.push(row);
				}
			}
		}
		break;
	
	case "Ind. Group":
		{
			head.push("ID");
		
			let cb = so.temp.check_box.value;
			let name_list = so.temp.name_list;
			
			for(let i = 0; i < name_list.length; i++){
				if(cb[i].check){
					let row = [];
					row.push(name_list[i].name); 
					ele.push(row);
				}
			}
		}
		break;
	
	case "Ind. Eff.":
		{
			let e = find_in(sp.ind_effect,so.spec.drop.te);
			if(e == undefined) alertp("Cannot find individual effect");
			head.push("ID");
			head.push("Value");
			for(let i = 0; i < sim.individual.length; i++){
				let ind = sim.individual[i];
				let row = [];
				row.push(ind.name); 
				row.push(ind.ie[e]);
				ele.push(row);
			}
		}
		break;
		
	default:
		error("Cannot simulate this data type");
		break;
	}

	if(ele.length == 0 && so.type != "Transition"){
		alert_help("No data","This data is empty.");
	}
	
	post({ type:input.type, head:head, ele:ele});
}


/// Works out the value of c after an event
function c_after_event(ev,sp)
{
	switch(ev.type){
	case EV_ENTER: return ev.c;
	case EV_MOVE: return ev.cf;
	case EV_LEAVE: return OUT; 
	case EV_TRANS: return sp.tra_gl[ev.trg].f;
	default: error("Cannot find event type"); break;
	}
}


/// Work out how filter is stratified
function get_stratify(filt,sp)
{
	let cl_list = [];
	let list = [];
	for(let cl = 0; cl < sp.cla.length; cl++){
		let fcl = filt.cla[cl];
		let claa = sp.cla[cl];
		if(fcl.radio.value == "File"){
			cl_list.push(cl);
			let comp_list = [];
			if(fcl.comp.length != claa.comp.length) error("comp size wrong");
			
			for(let c = 0; c < claa.comp.length; c++){
				comp_list.push(claa.comp[c].name);
			}
			list.push(comp_list);
			
			fcl.radio.value = "Comp";
		}
	}
	
	let comb = generate_co_list(list);
	return { cl_list:cl_list, comb:comb, list:list};
}


/// For a given stratification set the filter
function set_filt_strat(com,strat,filt)
{
	for(let i = 0; i < strat.cl_list.length; i++){
		let fcl = filt.cla[strat.cl_list[i]];
		let c_sel = strat.comb[com].index[i];
		for(let ci = 0; ci < fcl.comp.length; ci++){
			if(ci == c_sel) fcl.comp[ci].check = true; 
			else fcl.comp[ci].check = false; 
		}
	}
}


/// Sets up any equations used to simulate data
function data_sim_setup_eqn(filt,sp,cl_sel)
{		
	for(let cl = 0; cl < sp.ncla; cl++){
		if(cl == cl_sel){
			let tra = sp.cla[cl].tra;

			let tr_filt = filt.tra;

			switch(filt.trans_obs_model.value){
			case "on":
				for(let ti = 0; ti < tra.length; ti++){
					let prob_eqn = tr_filt[ti].prob_eqn;
					let te = prob_eqn.te;
					if(isNaN(te)){
						let eqn = setup_eqn(te);
						if(eqn.err == true){ alert_help("Error evaluating equation",eqn.msg); return;}
						prob_eqn.evaluate = eqn;
					}
				}
				break;
			}
		}
		else{
			let comp = sp.cla[cl].comp;
			let fcl = filt.cla[cl];
			switch(fcl.radio.value){
			case "ObsMod":
				for(let ci = 0; ci < sp.cla[cl].comp.length; ci++){
					let prob_eqn = fcl.comp[ci].prob_eqn;
					let te = prob_eqn.te;
					if(isNaN(te)){
						let eqn = setup_eqn(te);
						if(eqn.err == true){ alert_help("There was a problem evaluating this equation: "+te); return;}
						prob_eqn.evaluate = eqn;
					}
				}
				break;
			}
		}
	}
}
			

/// Gets the probability for a specific transition 
function get_tra_gl_prob(tri,t,filt,sp,cl_sel)
{
	let tr = sp.tra_gl[tri];
	if(tr.cl != cl_sel) return 0;
						
	let c = tr.i; if(c == SOURCE) c = tr.f;

	let prob = 1;
	let tr_filt = filt.tra;
		
	switch(filt.trans_obs_model.value){
	case "on":
		{	
			let prob_eqn = tr_filt[tr.tr].prob_eqn;
			
			let te = prob_eqn.te;
			if(isNaN(te)){
				let eqn = prob_eqn.evaluate;
				if(eqn == undefined) error("Equation not set"); 
				prob *= calculate(eqn,t);
			}
			else{
				prob *= Number(te);
			}
		}
		break;

	case "off":
		if(tr_filt[tr.tr].check != true) prob = 0;
		break;
	} 

	if(prob == 0) return 0;

	for(let cl = 0; cl < sp.ncla; cl++){
		let comp = sp.cla[cl].comp;
		
		if(cl != cl_sel){
			let ci = sp.comp_gl[c].cla_comp[cl];
		
			let fcl = filt.cla[cl];
	
			switch(fcl.radio.value){
			case "All":
				break;
				
			case "Comp":
				if(fcl.comp[ci].check != true) prob = 0;
				break;
				
			case "ObsMod":
				{
					let prob_eqn = fcl.comp[ci].prob_eqn;
					let te = prob_eqn.te;
				
					if(isNaN(te)){
						let eqn = prob_eqn.evaluate;
						if(eqn == undefined) error("Equation not set"); 
						prob *= calculate(eqn,t);
					}
					else{
						prob *= Number(te);
					}
				}
				break;
			
			default: error("get_comp_prob problem"); break;
			}	
		}
	}
	
	return prob;
}


/// Gets the probability of being in different global compartments
function get_comp_prob(t,filt,sp,cl_sel)
{
	// Gets the factor which each comartment is multiplied by
	let pos = [];
	for(let cl = 0; cl < sp.ncla; cl++){
		let comp = sp.cla[cl].comp;
		let row = [];
		
		if(cl != cl_sel){
			let fcl = filt.cla[cl];
			switch(fcl.radio.value){
			case "All":
				for(let ci = 0; ci < comp.length; ci++){
					row[ci] = 1;
				}
				break;
				
			case "Comp":
				for(let ci = 0; ci < comp.length; ci++){
					if(fcl.comp[ci].check == true) row[ci] = 1;
					else row[ci] = 0;
				}
				break;
				
			case "ObsMod":
				for(let ci = 0; ci < sp.cla[cl].comp.length; ci++){
					let prob_eqn = fcl.comp[ci].prob_eqn;
					let te = prob_eqn.te;
					if(isNaN(te)){
						let eqn = prob_eqn.evaluate;
						if(eqn == undefined) error("Equation not set"); 
						row[ci] = calculate(eqn,t);
					}
					else{
						row[ci] = Number(te);
					}
				}
				break;
			
			default: error("get_comp_prob problem"); break;
			}	
		}
		pos.push(row);
	}

	let posg = [];
	for(let c = 0; c < sp.comp_gl.length; c++){
		let cgl = sp.comp_gl[c];
		let fac = 1;
		for(let cl = 0; cl < sp.ncla; cl++){
			if(cl != cl_sel) fac *= pos[cl][cgl.cla_comp[cl]];
		}
		posg[c] = fac;
	}

	return posg;
}


/// Gets the transition probability
function get_tr_prob(t,filt,posg,sp,cl_sel)
{
	let tr_filt = filt.tra;

	let tr_prob = [];	
		
	for(let tri = 0; tri < sp.tra_gl.length; tri++){
		let tr = sp.tra_gl[tri];
						
		tr_prob[tri] = 0;
						
		if(tr.cl == cl_sel){
			let c = tr.i; if(c == SOURCE) c = tr.f;

			let prob = posg[c];
			if(prob != 0){
				switch(filt.trans_obs_model.value){
				case "on":
					{	
						let prob_eqn = tr_filt[tr.tr].prob_eqn;
						let te = prob_eqn.te;
						if(isNaN(te)){
							let eqn = prob_eqn.evaluate;
							if(eqn == undefined) error("Equation not set"); 
							tr_prob[tri] = prob*calculate(eqn,t);
						}
						else{
							tr_prob[tri] = prob*Number(te);
						}
					}
					break;
		
				case "off":
					if(tr_filt[tr.tr].check == true) tr_prob[tri] = prob;
					break;
				} 
			}
		}
	}

	return tr_prob;
}


/// Gets the equation value for observation model
function get_equation_value(eqn)
{
	if(isNaN(eqn.te)) post({type:"Error", te:"Cannot find equation value: "+eqn.te});
	
	return Number(eqn.te);
}


/// Gets timepoints specified in the bubble
function get_timepoints(t_start,t_end,so)
{
	let times = [];
	
	let val = so.time_gen;

	switch(so.time_radio.value){
	case "Periodic":
		{
			if(isNaN(val)) return "The value '"+val+"' must be a number";
			let dt = Number(val);
			if(dt <= 0) return "The value '"+val+"' must be positive";
			for(let t = t_start+dt; t <= t_end; t += dt){
				times.push(t);
			}
		}			
		break;
	
	case "Specified": 
		let spl = val.split(",");
		for(let k = 0; k < spl.length; k++){
			if(isNaN(spl[k])) return "The value '"+spl[k]+"' must be a number";
			let num = Number(spl[k]);
			if(num < t_start || num > t_end){
				return "The value '"+spl[k]+"' must be between the start '"+t_start+"' and end '"+t_end+"' times";
			}
			times.push(num);
		}
		break;
	
	case "Data":
		alertp("Should not be data");
		break;
	}
	
	return times;
}


/// Gets timepoints specified in the bubble
function get_ind_times(t_start,t_end,so,samp,hash,ind_list)
{
	let ind_times = [];
	
	for(let p = 0; p < samp.species.length; p++){
		ind_times[p]=[];
		let sosp = samp.species[p];
		for(let i = 0; i < sosp.individual.length; i++) ind_times[p][i] = [];
	}
		
	switch(so.time_radio.value){
	case "Periodic": case "Specified": 
		{	
			let times = get_timepoints(t_start,t_end,so);
			
			for(let p = 0; p < samp.species.length; p++){
				let sosp = samp.species[p];
				for(let i = 0; i < sosp.individual.length; i++){
					for(let j = 0; j < times.length; j++) ind_times[p][i].push(times[j]);
				}
			}
		}
		break;
	
	case "Data":
		{
			let sel = so.sel_data;
			let gen_so = model.sim_res.plot_filter.species[sel.p].gen_source[sel.i];
			
			switch(gen_so.type){
			case "Compartment": case "Transition": case "Diag. Test":
				{
					let ele = gen_so.table.ele;
					for(let r = 0; r < ele.length; r++){
						let ro = ele[r];
						let id = ro[0];
						let t = ro[1];
						
						let j = hash.find(id);
						if(j == undefined) alertp("Should not be undefined");
						else{
							let il = ind_list[j];
							let tf = Number(t);
							if(isNaN(tf)){
								if(t == "start") tf = t_start;
								else{
									if(t == "end") tf = t_end;
									else{
										alertp("Not a number");
									}
								}
							}
							if(gen_so.type == "Transition") tf += TINY;
							ind_times[il.p][il.i].push(tf);
						}
					}
				}
				break;
				
			default:
				alert_help("Cannot used this data source");
				break;
			}
		}			
		break;	
	}
	
	return ind_times;
}


/// Adds a data source onto the generate data
function add_gen_data(head,ele)
{
	let tab = { filename:"Simulated data", heading:head, ele:ele, ncol:head.length, nrow:ele.length};

	data.table.push(tab);
	edit_source.data_table_use = data.table.length-1;
	edit_source.table_loaded = true;

	if(edit_source.type == "Genetic"){
		switch(edit_source.spec.type_radio.value){
		case "matrix":
			{
				let c = edit_source.load_col.length-1;
				let ele = tab.ele;
				for(let r = 0; r < ele.length; r++){
					edit_source.load_col.push({heading:ele[r][c], desc:"Matrix data",type:"positive int"});	
				}
			}
			break;
		
		case "snp":
			{
				edit_source.load_col.pop();
				
				for(let i = 2; i < tab.ncol; i++){
					edit_source.load_col.push({heading:tab.heading[i], desc:"SNP column",type: "text"});
				}
			}
			break;
		}
	}

	for(let i = 0; i < head.length; i++){
		transfer_column(i,false);
	}
	
	edit_source.table.edit = true;
}


/// Copies simulated data into data
function copy_simulated_data()
{
	let rpf = model.sim_res.plot_filter;
	let sp = rpf.species;

	if(model.species.length != sp.length){
		alertp("The number of species in the model and simulated data do not agree");
		return;
	}
	
	for(let p = 0; p < model.species.length; p++){
		if(model.species[p].name != sp[p].name){
			alertp("The names of the species in the model and simulated data do not agree");
			return;
		}
	}
	
	for(let p = 0; p < model.species.length; p++){
		model.species[p].inf_source = copy(sp[p].gen_source);
		for(let j = 0; j < model.species[p].inf_source.length; j++){
			model.species[p].inf_source[j].info.siminf = "inf";
		}
	}
	
	model.inf_details.t_start = rpf.details.t_start;
	model.inf_details.t_end = rpf.details.t_end;
	model.inf_details.timestep = rpf.details.timestep;
}


/// Simulates genetic snp data 
function simulate_genetic_snp(t_start,t_end,ele,head,so)
{
	let result = sim_result;
	
	let mu = so.spec.mut_rate_eqn.te;
	if(isNaN(mu)){ alert_help("The mutation rate must be a number"); return;}
	mu = Number(mu);
	
	let seq_var = so.spec.seq_var_eqn.te;
	if(isNaN(seq_var) ){ alert_help("The sequence variation must be a non-negative number"); return;}
	seq_var = Number(seq_var);
	if(seq_var < 0){ alert_help("The sequence variation must be a non-negative number"); return;}
	
	let N = so.numbp;
	let f = so.frac_obs;

	head.push("ID");
	head.push("t");
	for(let i = 0; i < N; i++) head.push(so.spec.snp_root+(i+1));
	
	let samp = result.sample[0];

	// Makes a hash table of individuals
	let hash = new Hash();
	
	let ind_list = [];
	for(let p = 0; p < samp.species.length; p++){
		let sosp = samp.species[p];
		for(let i = 0; i < sosp.individual.length; i++){
			let ind = sosp.individual[i];
	
			let name = ind.name;
			let j = hash.find(name);
			if(j == undefined){
				j = ind_list.length;
				hash.add(name,ind_list.length);
				ind_list.push({name:name, p:p, i:i});
			}
		}
	}
	
	let ind_times = get_ind_times(t_start,t_end,so,samp,hash,ind_list);
	
	//let times = get_timepoints(t_start,t_end,so);
	//if(typeof times == 'string'){ alert_help("Problem with times",times); return;}
	
	let ref_seq = [];
	for(let j = 0; j < N; j++){
		ref_seq.push(Math.floor(4*Math.random()));
	}
	
	// Generates timelines in which an individual infects another 
	let ind_infect=[];
	for(let p = 0; p < samp.species.length; p++){
		ind_infect[p]=[];
		let indi = samp.species[p].individual;
		for(let i = 0; i < indi.length; i++){
			ind_infect[p][i]=[];
		}
	}
	
	// Makes a list of when infection enters from outside the system
	let inf_enter = [];
	for(let p = 0; p < samp.species.length; p++){
		let sp = result.species[p];
		let indi = samp.species[p].individual;
		for(let i = 0; i < indi.length; i++){
			let ind = indi[i];
			let c = ind.cinit;
			if(c != OUT && c != SOURCE && sp.comp_gl[c].infected == true){
				inf_enter.push({p:p, i:i, t:t_start});
			}
		}
	}

	for(let p = 0; p < samp.species.length; p++){	
		let indi = samp.species[p].individual;
		
		let sp = result.species[p];
		let tra_gl = sp.tra_gl;
		let comp_gl = sp.comp_gl;
		
		for(let i = 0; i < indi.length; i++){
			let ev = indi[i].ev;
			let c = indi[i].cinit; 
			let inf_flag = false;
			if(c != OUT && c != SOURCE && comp_gl[c].infected == true) inf_flag = true;
				
			for(let e = 0; e < ev.length; e++){
				c = c_after_event(ev[e],sp);
				let infe = ev[e].infection;
				if(infe){
					if(comp_gl[c].infected != true) error("Should be infected");
					if(inf_flag != false) error("Already infected");
					inf_flag = true;	
					
					if(infe.p == OUTSIDE_INF){
						inf_enter.push({p:p, i:i, t:ev[e].t});
					}
					else{
						ind_infect[p][i].push({type:"infected", t:ev[e].t});  
					
						ind_infect[infe.p][infe.i].push({type:"infect", t:ev[e].t, p:p, i:i}); 
					}
				}
				else{
					if(inf_flag == true){		
						if(c == OUT || c == SOURCE || c == SINK || comp_gl[c].infected != true){ 
							ind_infect[p][i].push({type:"recovery", t:ev[e].t});  
							inf_flag = false;
						}
					}
				}
			}
		}
	}
	
	// Orders ind_infect
	for(let p = 0; p < samp.species.length; p++){
		for(let i = 0; i < ind_infect[p].length; i++){
			ind_infect[p][i].sort( function(a, b){ return a.t-b.t});
		}
	}
	
	let node_list = [];
	
	// Adds new infections to the system
	for(let i = 0; i < inf_enter.length; i++){
		let ie = inf_enter[i];
		
		let new_seq = copy(ref_seq);
		let num = poisson_sample(seq_var+(ie.t-t_start)*mu);
		mutate_seq(num,new_seq);
		
		seq_sim_ind(ie.p,ie.i,new_seq,ie.t,ind_infect,ele,mu,f,ind_times,samp,node_list);
	}
	
	node_list.sort( function(a, b){ return a.i-b.i});
	
	let te = "";
	for(let r = 0; r < node_list.length; r++){
		let no = node_list[r];
		
		te += no.i+" "+no.t+"  ";
		for(let j = 0; j < no.inf_ev.length; j++){
			te += no.inf_ev[j].t+" "+no.inf_ev[j].type+" "
			if(no.inf_ev[j].type == "infect") te += no.inf_ev[j].infi+" ";
			te += no.inf_ev[j].mut_num+",  ";
		}	
		te += "\n";
	}
	
	ele.sort( function(a, b){ 
		let tea = a[0].toLowerCase();
		let teb = b[0].toLowerCase();
		if(tea > teb) return 1;
		else{
			if(tea < teb) return -1;
			else return a[1].t-b[1].t;
		}		
	});
}


/// Simulates matrix of genetic differences
function simulate_genetic_matrix(t_start,t_end,ele,head,so)
{
	let result = sim_result;
	
	let mu = so.spec.mut_rate_eqn.te;
	if(isNaN(mu)){ alert_help("The mutation rate must be a number"); return;}
	mu = Number(mu);
	
	let seq_var = so.spec.seq_var_eqn.te;
	if(isNaN(seq_var) ){ alert_help("The sequence variation must be a non-negative number"); return;}
	seq_var = Number(seq_var);
	if(seq_var < 0){ alert_help("The sequence variation must be a non-negative number"); return;}
	
	let f = so.frac_obs;

	let samp = result.sample[0];
	
	let inf_node = samp.inf_node;
	
	if(!samp.inf_node_store) samp.inf_node_store = copy(inf_node);
	else inf_node = copy(samp.inf_node_store);
	
	// Makes a hash table of individuals
	let hash = new Hash();
	
	let ind_list = [];
	for(let p = 0; p < samp.species.length; p++){
		let sosp = samp.species[p];
		for(let i = 0; i < sosp.individual.length; i++){
			let ind = sosp.individual[i];
	
			let name = ind.name;
			let j = hash.find(name);
			if(j == undefined){
				j = ind_list.length;
				hash.add(name,ind_list.length);
				ind_list.push({name:name, p:p, i:i});
			}
		}
	}
	
	let ind_times = get_ind_times(t_start,t_end,so,samp,hash,ind_list);
	
	// Adds in observations
	let obs = [];
	for(let n = 0; n < inf_node.length; n++){
		let infn = inf_node[n];
		
		let j = hash.find(infn.ind);
		let il = ind_list[j];
		let p = il.p;
		let ind = samp.species[p].individual[il.i];
		let ev = ind.ev;
		let t = infn.t;
		
		// Works out the time when the individual is no longer in the infectious state
		let e = 0;
	
		if(infn.from == "ENT" || infn.from == "OUT"){
			infn.num = poisson_sample(seq_var+(t-t_start)*mu);
		}
		else{
			while(e < ev.length && ev[e].t != t) e++;
			if(e == ev.length) alertp("Could not find event");
		}
		
		let sp = result.species[p];
		
		let tmax = t_end;
		while(e < ev.length){
			let c = c_after_event(ev[e],sp);
			if(c == OUT || c == SOURCE || c == SINK || sp.comp_gl[c].infected != true){
				tmax = ev[e].t;
				break;
			}
			e++;
		}

		let ti_list = ind_times[p][il.i];
		for(let k = 0; k < ti_list.length; k++){
			let tt = ti_list[k];
			if(tt > t && tt < tmax){
				if(Math.random() < f){
					infn.inf_ev.push({ty: "OBS", i:obs.length, t:tt, num:UNSET});
					obs.push({name:"Obs-"+(obs.length+1), n:n, ind:ind.name, t:tt});
				}
			}					
		}
		
		infn.inf_ev.sort( function(a, b){ return a.t-b.t});
		
		for(let k = 0; k < infn.inf_ev.length; k++){
			let tt = infn.inf_ev[k].t;
			infn.inf_ev[k].num = poisson_sample(mu*(tt-t));
			t = tt;
		}
	}	
	
	// Corrects inf_node.num (because observations have been added)
	for(let n = 0; n < inf_node.length; n++){
		let infn = inf_node[n];
		if(infn.from != "ENT" && infn.from != "OUT"){
			let nn = infn.from;
			let infe = inf_node[nn].inf_ev;
			let num = infn.num;
			let ninf = 0;
			let e = 0; 
			for(e = 0; e < infe.length; e++){
				if(infe[e].ty == "INF"){
					if(ninf == num){ infn.num = e; break;}
					ninf++;
				}
			}
			
			if(e == infe.length){
				alertp("Could not correct");
			}
		}
	}
	
	// Calculates matrix of observations
	let mat = [];
	for(let j = 0; j < obs.length; j++) mat[j]=[];
	
	let obs_visit=[];
	let obs_infev=[];
	for(let n = 0; n < inf_node.length; n++){
		obs_visit[n]=[];
		obs_infev[n]=[];
		let inf_ev = inf_node[n].inf_ev;
		for(let e = 0; e < inf_ev.length; e++){
			obs_infev[n][e]=[];
			if(inf_ev[e].ty == "OBS"){
				obs_infev[n][e].push({m:inf_ev[e].i, sum:0});
			}
		}
	}
		
	for(let n = 0; n < inf_node.length; n++){
		let inf_ev = inf_node[n].inf_ev;
		for(let k = 0; k < inf_ev.length; k++){
			let ie = inf_ev[k];
			if(ie.ty == "OBS"){
				let m = ie.i;
				let kk = k, nn = n;
				let sum = 0;
				do{
					let infn = inf_node[nn];
					
					let oie = obs_infev[nn][kk];
					for(let j = 0; j < oie.length; j++){
						//if(infn.inf_ev[kk].ty == "OBS"){
						//let mm = infn.inf_ev[kk].i;
						let mm = oie[j].m;
						if(mat[m][mm] == undefined){
							let summ = oie[j].sum;
							mat[m][mm] = sum+summ;
							mat[mm][m] = sum+summ;
						}
					}
					oie.push({m:m, sum:sum});
							
					sum += infn.inf_ev[kk].num;
				
					kk--;
					if(kk == -1){ // Goes up the node
						if(inf_node[nn].from == "ENT" || inf_node[nn].from == "OUT"){
							obs_visit[nn].push({m:m,sum:sum});
							break;
						}					
						kk = inf_node[nn].num;
						nn = inf_node[nn].from;
					}
				}while(true);
			}
		}
	}
	
	let enter_list = [];
	for(let n = 0; n < inf_node.length; n++){
		if(inf_node[n].from == "ENT" || inf_node[n].from == "OUT"){
			if(obs_visit[n].length > 0){
				enter_list.push(n);
			}
		}
	}
	
	// Accounts for differences between observations on different origins
	if(enter_list.length > 1){
		for(let j = 0; j < enter_list.length; j++){
			let n = enter_list[j];
			let ov = obs_visit[n];
			let num = inf_node[n].num;
			for(let k = 0; k < ov.length; k++){
				let m = ov[k].m;
				let sum_k = ov[k].sum;
				for(let jj = j+1; jj < enter_list.length; jj++){
					let nn = enter_list[jj];
					let ovv = obs_visit[nn];
					let numm = inf_node[nn].num;
					for(let kk = 0; kk < ovv.length; kk++){
						let mm = ovv[kk].m;
						let sum_kk = ovv[kk].sum;
					
						if(mat[m][mm] != undefined)	alertp("Should not be defined");
						mat[m][mm] = sum_k + sum_kk+ num + numm;
						mat[mm][m] = mat[m][mm];
					}
				}
			}
		}
	}
	
	head.push("ID");
	head.push("t");
	head.push("Obs");
	for(let i = 0; i < obs.length; i++){
		head.push(obs[i].name);
	}
	
	for(let r = 0; r < obs.length; r++){
		let ob = obs[r];
		let row=[];
		row.push(ob.ind);
		row.push(ob.t);
		row.push(ob.name);
		for(let i = 0; i < obs.length; i++){
			row.push(mat[r][i]);
		}
		ele.push(row);
	}

	set_genetic_matrix_columns(ele,so);
}


/// Sets sequence for individual and evolves in time
function seq_sim_ind(p,i,new_seq,t,ind_infect,ele,mu,f,ind_times,samp,node_list)
{
	let ind_inf = ind_infect[p][i];
	let node = {t:t, p:p, i:i, inf_ev:[]};
	
	let e = 0;
	while(e < ind_inf.length && ind_inf[e].t <= t) e++;
	
	let ti_list = ind_times[p][i];
	
	let m = 0; 
	while(m < ti_list.length && ti_list[m] <= t) m++;
	
	let inf_flag = true;
	
	do{
		let t_next = LARGE;
		
		if(m < ti_list.length && ti_list[m] < t_next) t_next = ti_list[m];
		
		if(e < ind_inf.length && ind_inf[e].t < t_next) t_next = ind_inf[e].t;
		
		if(t_next == LARGE) break;
		
		let num = poisson_sample(mu*(t_next-t));
		
		mutate_seq(num,new_seq);
		
		t = t_next;
		
		if(m < ti_list.length && ti_list[m] == t){
			if(inf_flag == true){
				node.inf_ev.push({type:"obs", t:t, mut_num:num});
			
				let row=[];
				row.push(samp.species[p].individual[i].name);
				row.push(t);
				for(let j = 0; j < new_seq.length; j++){
					let ch;
					switch(new_seq[j]){
					case 0: ch = "A"; break;
					case 1: ch = "C"; break;
					case 2: ch = "T"; break;
					case 3: ch = "G"; break;
					}
					row.push(ch);
				}
				ele.push(row);
			}
			m++;
		}
		
		if(e < ind_inf.length && ind_inf[e].t == t){
			let indinf = ind_inf[e];
			
			if(indinf.type == "infected") break;
			else{
				if(indinf.type == "infect"){
					node.inf_ev.push({type:"infect", t:t, infi:indinf.i, mut_num:num});
					seq_sim_ind(indinf.p,indinf.i,copy(new_seq),t,ind_infect,ele,mu,f,ind_times,samp,node_list);
				}
				else{
					if(indinf.type != "recovery") error("Option not recognised");
					inf_flag = false;
				}
			}
			e++;
		}
	}while(true);
	
	node_list.push(node);
}


/// Mutates an element of the sequence
function mutate_seq(num,seq)
{
	for(let k = 0; k < num; k++){
		let i = Math.floor(seq.length*Math.random());
		seq[i] = (seq[i]+1+Math.floor(3*Math.random()))%4;
	}
}

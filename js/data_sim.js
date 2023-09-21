"use strict";

/// Simulate data

function sim_data()
{
	let so = inter.edit_source;
	
	let result = sim_result;
	
	let t_start = Number(result.details.t_start);
	let t_end = Number(result.details.t_end);
	let timestep = Number(result.details.timestep);
	
	let p = model.get_p();
	let sp = result.species[p];
	
	let sim = result.sample[0].species[p];
	
	let spec = so.spec;
		
	let head = [];
	let ele = [];
	
	switch(so.type){
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
		
	case "Add Ind.":
		{
			head.push("ID");
			head.push("t");
			for(let cl = 0; cl < sp.ncla; cl++) head.push(sp.cla[cl].name);
			
			for(let i = 0; i < sim.individual.length; i++){
				let ind = sim.individual[i];
				let row = [];
				row.push(ind.name);
				
				let cinit = ind.cinit
				if(ind.cinit != "Source"){ 
					if(cinit == "out"){
						let ev = ind.ev[0]
						row.push(ev.t);
						cinit = ev.c;;
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
		break;
		
	case "Remove Ind.":
		{
			head.push("ID");
			head.push("t");
		
			for(let i = 0; i < sim.individual.length; i++){
				let ind = sim.individual[i];
				let row = [];
				row.push(ind.name);
				
				let ev = ind.ev[ind.ev.length-1];
				if(ev.type == "leave"){
					row.push(ev.t);	
					ele.push(row);
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
				
				for(let e = 0; e < ind.ev.length; e++){
					let ev = ind.ev[e];
					if(ev.type == "move" && ev.cl == cl_sel){
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
		break;
		
	case "Compartment":
		{
			let times = get_timepoints(t_start,t_end);
			head.push("ID");
			head.push("t");
			let cl_drop = spec.cl_drop;
			head.push(cl_drop.te);
			
			let cl_sel = find(sp.cla,"name",cl_drop.te);
			
			for(let i = 0; i < sim.individual.length; i++){
				let ind = sim.individual[i];
				for(let k = 0; k < times.length; k++){
					let t = times[k];
				
					let c = ind.cinit;
					let e = 0; 
					while(e < ind.ev.length && ind.ev[e].t < t){
						let ev = ind.ev[e];
						switch(ev.type){
						case "enter": c = ev.c; break;
						case "move": c = ev.cf; break;
						case "leave": c = "out"; break;
						case "trans": c = sp.tra_gl[ev.trg].f; break;
						}
						e++;
					}
				
					if(c != "Source" && c != "Sink" && c != "out"){
						let row = [];
						row.push(ind.name);
						row.push(t);
						let cgl = sp.comp_gl[c];
					
						row.push(sp.cla[cl_sel].comp[cgl.cla_comp[cl_sel]].name);
						ele.push(row);
					}
				}
			}
		}
		break;
		
	case "Transition": case "Source": case "Sink":
		{
			let tmin = t_start, tmax = t_end;
			if(spec.time_radio.value == "Spec"){
				tmin = Number(spec.time_start);
				tmax = Number(spec.time_end);
			}
			
			let cl_sel = spec.tr_drop.cl;
			let tri = spec.tr_drop.i;
			
			head.push("ID");
			head.push("t");
		
			for(let i = 0; i < sim.individual.length; i++){
				let ind = sim.individual[i];
				let row = [];
					
				let flag = false;
				for(let e = 0; e < ind.ev.length; e++){
					let ev = ind.ev[e];
					if(ev.type == "trans" && ev.t >= tmin && ev.t < tmax){
						let trg = sp.tra_gl[ev.trg];
						if(trg.cl == cl_sel && trg.tr == tri){
							row.push(ind.name);
							row.push(ev.t);
							flag = true;
						}
					}
				}
				
				if(flag == false){
					row.push(ind.name);
					row.push("no");
				}
				ele.push(row);
			}
		}
		break;
	
	case "Diag. Test":
		{
			let times = get_timepoints(t_start,t_end);
			
			let cl_drop = spec.cl_drop;
			let cl_sel = find(sp.cla,"name",cl_drop.te);
			
			head.push("ID");
			head.push("t");
			head.push("Result");

			let Se = Number(spec.Se_eqn.te);
			let Sp = Number(spec.Sp_eqn.te);
			let neg = spec.neg_result;
			let pos = spec.pos_result;
		
			for(let i = 0; i < sim.individual.length; i++){
				let ind = sim.individual[i];
				for(let k = 0; k < times.length; k++){
					let t = times[k];
				
					let c = ind.cinit;
					let e = 0; 
					while(e < ind.ev.length && ind.ev[e].t < t){
						let ev = ind.ev[e];
						switch(ev.type){
						case "enter": c = ev.c; break;
						case "move": c = ev.cf; break;
						case "leave": c = "out"; break;
						case "trans": c = sp.tra_gl[ev.trg].f; break;
						}
						e++;
					}
				
					if(c != "Source" && c != "Sink" && c != "out"){
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
		break;
		
	case "Population":
		{
			let times = get_timepoints(t_start,t_end);
			
			let filt = spec.filter;

			let pos = [];
			for(let cl = 0; cl < sp.ncla; cl++){
				let row = [];
				
				let fcl = filt.cla[cl];
				if(fcl.radio.value == "All"){
					for(let ci = 0; ci < sp.cla[cl].comp.length; ci++){
						row[ci] = true;
					}
				}
				else{
					for(let ci = 0; ci < sp.cla[cl].comp.length; ci++){
						if(fcl.comp[ci].check == true) row[ci] = true;
						else row[ci] = false;
					}
				}
				pos.push(row);
			}
			
			let posg = [];
			for(let c = 0; c < sp.comp_gl.length; c++){
				let cgl = sp.comp_gl[c];
				let cl = 0; while(cl < sp.ncla && pos[cl][cgl.cla_comp[cl]] == true) cl++;
				if(cl == sp.ncla) posg[c] = true;
				else posg[c] = false;
			}
			
			head.push("t");
			head.push("Population");
			
			for(let k = 0; k < times.length; k++){
				let t = times[k];
				
				let popc = [];
				
				switch(sp.type){
				case "Individual":
					{
						for(let c = 0; c < sp.comp_gl.length; c++) popc[c] = 0;
					
						for(let i = 0; i < sim.individual.length; i++){
							let ind = sim.individual[i];
					
							let c = ind.cinit;
							let e = 0; 
							while(e < ind.ev.length && ind.ev[e].t < t){
								let ev = ind.ev[e];
								switch(ev.type){
								case "enter": c = ev.c; break;
								case "move": c = ev.cf; break;
								case "leave": c = "out"; break;
								case "trans": c = sp.tra_gl[ev.trg].f; break;
								}
								e++;
							}
							
							if(!isNaN(c)) popc[c]++;
						}
					}
					break;
					
				case "Population":
					{
						let ti = Math.floor((t-t_start)/timestep);
						
						for(let c = 0; c < sp.comp_gl.length; c++) popc = sim.timeline[ti].cpop;
					}
					break;
				}
				
				let row = [];
				row.push(t);
				
				let sum = 0;
				for(let c = 0; c < sp.comp_gl.length; c++){
					if(posg[c] == true) sum += popc[c];
				}
				if(sum < 0) sum = 0;
				row.push(sum);
				ele.push(row);
			}
		}
		break;
		
	case "Pop. Trans.":
		{
			let dt = Number(so.time_step);
			
			let cl_drop = spec.cl_drop;
			let cl_sel = find(sp.cla,"name",cl_drop.te);
		
		
			let filt = spec.filter;
		
			let pos = [];  // Possible compartments
			for(let cl = 0; cl < sp.ncla; cl++){
				let row = [];
				if(cl == cl_sel){
					for(let ci = 0; ci < sp.cla[cl].comp.length; ci++){
						row[ci] = true;
					}
				}
				else{
					let fcl = filt.cla[cl];
					if(fcl.radio.value == "All"){
						for(let ci = 0; ci < sp.cla[cl].comp.length; ci++){
							row[ci] = true;
						}
					}
					else{
						for(let ci = 0; ci < sp.cla[cl].comp.length; ci++){
							if(fcl.comp[ci].check == true) row[ci] = true;
							else row[ci] = false;
						}
					}
				}
				pos.push(row);
			}
			
			let posg = [];
			for(let c = 0; c < sp.comp_gl.length; c++){
				let cgl = sp.comp_gl[c];
				let cl = 0; while(cl < sp.ncla && pos[cl][cgl.cla_comp[cl]] == true) cl++;
				if(cl == sp.ncla) posg[c] = true;
				else posg[c] = false;
			}
		
			let tr_filt = spec.filter.tra;
		
			let trpos = [];
			for(let tri = 0; tri < sp.tra_gl.length; tri++){
				let tr = sp.tra_gl[tri];
				
				trpos[tri] = false;
				
				if(tr.cl == cl_sel){
					if(tr_filt[tr.tr].check == true){
						if(posg[tr.i] == true || posg[tr.f] == true){
							trpos[tri] = true;
						}
					}
				} 
			}
		
			head.push("t");
			head.push("Population");
			
			let res = [];
			let ndiv = 1+Math.floor(ALMOST_ONE*(t_end-t_start)/dt);
			for(let d = 0; d < ndiv; d++) res[d] = 0;
			
			switch(sp.type){
			case "Individual":
				{
					for(let i = 0; i < sim.individual.length; i++){
						let ind = sim.individual[i];
					
						for(let e = 0; e < ind.ev.length; e++){
							let ev = ind.ev[e];
							if(ev.type == "trans"){
								if(trpos[ev.trg] == true){
									let d = Math.floor(ALMOST_ONE*((ev.t-t_start)/dt));
									res[d]++;
								}
							}
						}
					}
				}
				break;
				
			case "Population":
				{
					for(let ti = 0; ti <sim.timeline.length-1; ti++){
						let d = Math.floor(ti*timestep/dt);
						for(let trg = 0; trg < sp.tra_gl.length; trg++){
							if(trpos[trg] == true){
								res[d] += sim.timeline[ti].transnum[trg];
							}
						}
					}
				}
				break;
			}
				
			for(let d = 0; d < ndiv; d++){
				let row = [];
				row.push(t_start + d*dt);
				row.push(res[d]);
				ele.push(row);
			}
		}
		break;
	}

	add_gen_data(head,ele);
}


/// Gets timepoints specified in the bubble
function get_timepoints(t_start,t_end)
{
	let times = [];
	
	let so = inter.edit_source;
	
	let val = so.time_gen;
	
	switch(so.time_radio.value){
	case "Periodic":
		let dt = Number(val);
		for(let t = t_start+dt; t <= t_end; t += dt){
			times.push(t);
		}			
		break;
		
	case "Fixed time": 
		times.push(Number(val)); 
		break;
		
	case "Specified": 
		let spl = val.split(",");
		for(let k = 0; k < spl.length; k++){
			times.push(Number(spl[k]));
		}
		break;
	}
	
	return times;
}


/// Adds a data source onto the generate data
function add_gen_data(head,ele)
{
	let tab = { filename:"Simulated data", heading:head, ele:ele, ncol:head.length, nrow:ele.length};

	data.table.push(tab);
	inter.edit_source.data_table_use = data.table.length-1;
	inter.edit_source.table_loaded = true;

	for(let i = 0; i < head.length; i++){
		transfer_column(i);
	}
	
	inter.edit_source.table.edit = true;
}


/// Copies generated data into data
function copy_simulated_data(p_from,p_to)
{
	if(model.species.length != sim_result.species.length){
		alertp("The number of species in the model and simulated data do not agree");
		return;
	}
	
	for(let p = 0; p < model.species.length; p++){
		if(model.species[p].name != sim_result.species[p].name){
			alertp("The names of the species in the model and simulated data do not agree");
			return;
		}
	}
	
	for(let p = 0; p < model.species.length; p++){
		model.species[p].inf_source = copy(sim_result.species[p].gen_source);
		for(let j = 0; j < model.species[p].inf_source.length; j++){
			model.species[p].inf_source[j].info.siminf = "inf";
		}
	}
	
	model.inf_details.t_start = model.sim_details.t_start;
	model.inf_details.t_end = model.sim_details.t_end;
	model.inf_details.timestep = model.sim_details.timestep;

	change_page({pa:"Inference", su:"Data"});
}


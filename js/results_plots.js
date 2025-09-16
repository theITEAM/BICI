"use strict";
// Functions which set up the plots from results

/// Adds buttons associated with a population plot
function add_pop_buts(res,lay)
{
	if(inter.graph.init == "loading") return;
	
	let cx = corner.x;
	let cy = corner.y;
	
	let ti = "Populations";
	
	let rpf = res.plot_filter;	
	
	if(rpf.sel_view.te == "Data") ti = rpf.sel_popdata.te;
	
	rpf = get_inf_res().plot_filter;
	
	cy = lay.add_title(ti,cx,cy,{te:pop_plot_text});
	
	cy += 1;

	if(inter.graph.init != true){
		inter.graph.init = "loading"; 
		start_worker("Graph pop",res_worker(res));
		return;
	}

	switch(plot_variety(inter.graph.type)){
	case "Line plot": inter.graph.create(cx,cy,graph_width-right_menu_width,29,lay); break;
	case "Comp plot": inter.graph.create(1,cy-1,graph_width-right_menu_width+0.5,lay.dy-cy,lay); break;
	case "Density plot": inter.graph.create(1,cy-1,graph_width-right_menu_width+0.5,lay.dy-cy,lay); break;
	case "No graph plot": inter.graph.no_graph(cx,cy,graph_width-right_menu_width,29,lay); break;
	default: error("No graph option"); break;
	}
}


/// Outpus the percentage complete
function output_percent(stepp,i,N)
{
	if(i%stepp == 0) percent(100*i/N);
}


/// Calculates quantities needed for the graph
function graph_pop_calculate(result,rpf,burn,p)
{
	// Sets classification
	let rpf2 = rpf.species[p];
	
	let cl = rpf2.sel_class.cl;
	let filter = rpf2.filter;

	// Incorprates filters
	let popfilt = get_popfilt(rpf2,"Populations");
	
	// Sets the chain
	let chsel = rpf.sel_chain.te;
	
	let sp = result.species[p];
	let claa = sp.cla[cl];
	
	let sample = result.sample;
	
	// Sets samples
	let imin = 0, imax = sample.length;
	let isel = get_isel(rpf); if(isel != "All"){ imin = isel; imax = isel+1; burn = 0;}

	let view = rpf.sel_view.te;

	if(result.plot_average || ((view == "Compartment" || view == "Density") && isel == "All" && chsel == "All")){
		imin = 0; imax = 1;
		isel = "All";
		sample = result.average;
	}			
	
	let data = [];
	
	let key = [];
	
	let line_max = false;
	
	let timepoint = result.timepoint;
	let T = timepoint.length;
	
	if(result.siminf == "ppc"){
		let det = result.details;
		if(det.ppc_t_start != det.t_start){
			data.push({te:"Start", col:BLACK, thick:1, x:Number(det.ppc_t_start), type:"VertLine"});
		}
		
		if(det.ppc_t_end != det.t_end){
			data.push({te:"End", col:BLACK, thick:1, x:Number(det.ppc_t_end), type:"VertLine"});
		}
	}
	
	if(view == "Data"){                              // If data is viewed 
		let seldata = rpf.sel_popdata;
		
		let so = result.species[seldata.p].inf_source[seldata.i];
	
		let strat_filt = seldata.strat_filt;

		// Works out which transitions apply to data
		let gfilt = [];	
		for(let c = 0; c < sp.comp_gl.length; c++){
			let cgl = sp.comp_gl[c];
			
			let flag = false;
			
			let cl2;
			for(cl2 = 0; cl2 < sp.ncla; cl2++){
				let filt = so.spec.filter.cla[cl2];
				let cc = cgl.cla_comp[cl2];
			
				if(!(filt.radio.value == "All" || (filt.radio.value == "ObsMod" && filt.comp[cc].prob_eqn.te != "0") ||
        				filt.comp[cc].check == true || strat_filt[cl2] == cc)) break;
			}
			if(cl2 == sp.ncla) flag = true;
			gfilt.push(flag);
		}
	
		let line_post = [];
		{                                              // Posterior
			
			for(let i = imin; i < imax; i++){
				let samp = sample[i];
				if(samp.num >= burn && !(chsel != "All" && samp.chain != chsel)){
					let sa = samp.species[p];
					
					let point = [];
				
					let y=[];
					for(let t = 0; t < T; t++) y[t] = 0;
					
					for(let k = 0; k < sp.comp_gl.length; k++){
						if(gfilt[k] == true) add_tl(y,sa.cpop_tl[k],result);
					}
					
					for(let t = 0; t < T; t++){					
						point.push({x:timepoint[t], y:y[t]});
					}
					
					line_post.push(point);
				}
			}
		
			add_line(view,line_post,BLUE,posterior_name(result),data,key);
		}
		
		{                                              // Data
			let line_mean = get_line_mean(line_post);
			
			let line = [];
			let tab = so.table;
		
			let spec = so.spec;
			let dist = spec.obs_error_dist.value;
			let ty = spec.obs_error.value;
		
			let col = BLACK;

			let co = tab.ncol-1;

			for(let j = 0; j < tab.nrow; j++){
				let t = Number(tab.ele[j][0]);
				
				let dfilt=[];
				for(let k = 0; k < so.load_col.length; k++){
					let cl = so.load_col[k].cl;
					if(cl != undefined){
						let co = tab.ele[j][k];
						let c = hash_find(sp.cla[cl].hash_comp,co);
						if(c == undefined) error("Could not find "+co);
						dfilt[cl] = c;					
					}
				}
						
				let cl = 0; while(cl < sp.ncla && dfilt[cl] == strat_filt[cl]) cl++;
				
				if(cl == sp.ncla){
					let m = find(so.load_col,"heading","Population");
					if(m == undefined) error("Cannot find population");
					
					let val_data = Number(tab.ele[j][m]);
						
					// Generate error bar on posterior
					let ti = get_ti(t,timepoint);
					if(ti != undefined){
						let val = line_mean[ti]; 
						
						let sd;
						switch(dist){
						case "Normal": 
							switch(ty){
							case "percent": sd = Number(spec.percent)*val/100; break;
							case "sd": sd = Number(spec.sd); break;
							case "file":
								{
									let m = find(so.load_col,"heading","sd");
									if(m == undefined) error("Cannot find sd");
									sd = Number(tab.ele[j][m]); break;
								}
								break;
							}
							break;
							
						case "Poisson":
							sd = Math.sqrt(val);
							break;
		
						case "Negative binomial":
							{
								let p;
								switch(spec.obs_error_p.value){
								case "p": p = Number(spec.p); break;
								case "file":
									{
										let m = find(so.load_col,"heading","p");
										if(m == undefined) error("Cannot find p");
										p = Number(tab.ele[j][m]); break;
									}
									break;
								}
								sd = Math.sqrt(val/p);
							}
							break;
						}
						
						if(!isNaN(t) && !isNaN(val) && !isNaN(sd)){	
							let ymin = val-sd; if(ymin < 0) ymin = 0;
							data.push({type:"ErrorBar", x:t, ymin:ymin, y:val, ymax:val+sd, col:col});
						}
					}
					
					data.push({type:"Cross", x:t, y:val_data, col:RED});	
				}
			}
			key.push({type:"Cross", te:"Data", col:RED});
		}
	}
	else{
		let stepp = get_percent_step(claa.ncomp);
		
		let c_list_rel = get_clist_rel(cl,popfilt,sp);
	
		let fraction = get_fraction(filter);
		let fraction_cl = fraction_in_cl(cl,filter);
		
		let c_list_rel_denom;
		let pos_op;
		if(fraction){
			pos_op = "positive";
			let popfilt_denom = get_popfilt(rpf2,"Populations","denom");
		
			let op; if(fraction_cl) op = "all";
		
			c_list_rel_denom = get_clist_rel(cl,popfilt_denom,sp,op);
		}
		
		for(let c = 0; c < claa.ncomp; c++){
			output_percent(stepp,c,claa.ncomp);
			
			if(!(popfilt[cl] != undefined && popfilt[cl][c] == false)){
				let col = claa.comp[c].col;
				
				let sh = c*sp.comp_mult[cl];
			
				let c_list = [];
				for(let k = 0; k < c_list_rel.length; k++){
					c_list.push(c_list_rel[k]+sh);
				}
				
				let c_list_denom = [];
				if(fraction){
					if(fraction_cl) sh = 0;
					for(let k = 0; k < c_list_rel_denom.length; k++){
						c_list_denom.push(c_list_rel_denom[k]+sh);
					}
				}
				
				let line = [];
			
				if(view == "Slice"){
					let t = Number(rpf.slice_time);
					
					let vec = [];
					
					let ti = get_ti(t,timepoint);
					if(ti != undefined){
						for(let i = imin; i < imax; i++){
							let samp = sample[i];
							if(samp.num >= burn && !(chsel != "All" && samp.chain != chsel)){
								let sa = samp.species[p];
							
								let pop = 0;
								for(let ij = 0; ij < c_list.length; ij++){
									pop += get_val_tl(ti,sa.cpop_tl[c_list[ij]],result,pos_op);
								}
								
								if(fraction){
									let pop_denom = 0;
									for(let ij = 0; ij < c_list_denom.length; ij++){			
										pop_denom += get_val_tl(ti,sa.cpop_tl[c_list_denom[ij]],result,pos_op);
									}
						
									if(pop_denom == 0) pop = 1;
									else pop /= pop_denom;
								}
								
								vec.push(pop);
							}
						}

						//key.push({type:"Line", te:claa.comp[c].name, col:col, c:c});
					}
					
					let show_mean = rpf.dist_settings.show_mean.check;
					if(vec.length > 0){
						let int_req, max = LARGE; 
						if(!fraction) int_req = "integer";
						else max = 1;
						distribution_add_data_line(vec,claa.comp[c].name,col,data,key,rpf,show_mean,0,max,int_req);
					}
				}
				else{
					for(let i = imin; i < imax; i++){
						let samp = sample[i];
						
						if(samp.num >= burn && !(chsel != "All" && samp.chain != chsel)){
							let sa = samp.species[p];
							
							let point = [];
							
							let y=[];
							for(let t = 0; t < T; t++) y[t] = 0;
							
							for(let ij = 0; ij < c_list.length; ij++){
								add_tl(y,sa.cpop_tl[c_list[ij]],result,pos_op);
							}
							
							if(fraction){
								let y_denom=[];
								for(let t = 0; t < T; t++) y_denom[t] = 0;
						
								for(let ij = 0; ij < c_list_denom.length; ij++){			
									add_tl(y_denom,sa.cpop_tl[c_list_denom[ij]],result,pos_op);
								}
								
								for(let t = 0; t < T; t++){
									let st = y[t];
									if(y_denom[t] == 0) y[t] = 1;
									else y[t] /= y_denom[t];
								}
							}
							
							for(let t = 0; t < T; t++){ 
								point.push({x:timepoint[t], y:y[t]});
							}
							
							line.push(point);
						}
					}
				
					switch(view){
					case "Graph": case "Graph (all)": case "Graph (split)": case "Graph (lines)":
					
						for(let j = 0; j < line.length; j++){
							data.push({point:line[j], col:claa.comp[c].col, view:view, type:"Line"});
							
							if(view == "Graph (lines)" && data.length > GRAPH_LINE_MAX){
								line_max = true; break;
							}
						}
						
						key.push({type:"Line", te:claa.comp[c].name, col:col, c:c});
						break;
						
					case "Graph (CI)": case "Data":
						{
							let line_stats = get_line_stats(line);
							data.push({point:line_stats, col:col, type:"Line CI"});
							key.push({type:"Line", te:claa.comp[c].name, col:col, c:c});
						}
						break;
						
					case "Compartment": case "Density":
						{
							let line_stats = get_line_mean(line);
							data.push({y_vec:line_stats, col:col, c:c, type:"Comp Video"});
						}
						break;
						
					default: error("view option error1"); break;
					}
				}
			}
			
			if(key.length > KEY_LINE_MAX){ line_max = true; break;}
		}
	}
			
	switch(view){
	case "Graph": case "Graph (all)": case "Graph (split)": 
	case "Graph (lines)": case "Graph (CI)": case "Data":
		post({type:"Graph define", variety:"Population", view:view, data:data, op:{ def_xrange:time_range(result.details), x_label:"Time", y_label:"Population", key:key, line_max:line_max}});
		break;
		
	case "Compartment": case "Density": 
		post({type:"Graph define", variety:"Population", view:view, data:data, op:{p:p, cl:cl, timepoint:result.timepoint, number_comp:true}});
		break;
		
	case "Slice":
		post({type:"Graph define", variety:"Distribution", view:"Graph", data:data, op:{x_label:"Population", x_param:true, y_label:"Probability", key:key, line_max:line_max, yaxis:false}});
		break;
		
	default: error("Type not recognised"); break;
	}
}

	
/// Gets name which is posterior under inference
function posterior_name(result)
{
	switch(result.siminf){
	case "sim": return "Simulation";
	case "inf": return "Posterior";
	case "ppc": return "Post. Sim.";
	default: return "Unset";
	}
}
			

/// Gets ti from time t 
function get_ti(t,tp)
{
	let N = tp.length;
	if(t < tp[0] || t > tp[N-1]) return;
	if(tp.length == 1) return 0;
	
	let dt = tp[1]-tp[0];
	let i = Math.floor((t-tp[0])/dt); if(i == N) i = N-1;
	return i;
}


/// Used to get a default xrange for graphs
function time_range(details)
{
	return {min:Number(details.t_start), max:Number(details.t_end)};
}


/// Determines if fraction is set for filter
function get_fraction(filter)
{
	for(let i = 0; i < filter.length; i++){
		if(filter[i].fraction.check == true) return true;
	}
	return false;
}


/// Determines if fraction is in a specific cl
function fraction_in_cl(cl,filter)
{
	for(let i = 0; i < filter.length; i++){
		if(filter[i].fraction.check == true && filter[i].cl == cl) return true;
	}
	return false;
}


/// Gets population filter
function get_popfilt(rpf2,name,op)
{
	let popfilt = [];
	
	for(let i = 0; i < rpf2.filter.length; i++){
		if(apply_filter(rpf2,i,name)){
			let rpf3 = rpf2.filter[i];
			
			let vec = [];

			if(op == "denom" && rpf3.fraction.check){
				for(let j = 0; j < rpf3.comp_filt.length; j++) vec.push(true);
			}
			else{
				if(rpf3.radio.value == "single"){
					for(let c = 0; c < rpf3.comp_name.length; c++){
						if(rpf3.comp_name[c] == rpf3.radio_sel.value) vec.push(true);
						else vec.push(false);
					}
				}
				else{
					for(let j = 0; j < rpf3.comp_filt.length; j++){
						if(rpf3.comp_filt[j].check == true) vec.push(true);
						else vec.push(false);
					}
				}
			}
			popfilt[rpf3.cl] = vec;
		}
	}
	return popfilt;
}

	
/// Determines if two vectors are the same
function vec_equal(vec1,vec2)
{
	if(vec1.length != vec2.length) return false;
	for(let i = 0; i < vec1.length; i++) if(vec1[i] != vec2[i]) return false;
	return true;
}


/// Based on a population filter gets relative position of all comartments allowed 
function get_clist_rel(cl,popfilt,sp,op)
{
	let clist_rel = [];
	for(let k = 0; k < sp.comp_gl.length; k++){
		let co = sp.comp_gl[k];
		let cl2;
		for(cl2 = 0; cl2 < sp.ncla; cl2++){
			if(cl2 == cl){
				if(op != "all" && co.cla_comp[cl2] != 0) break;
			}
			else{
				if(popfilt[cl2] != undefined){
					if(popfilt[cl2][co.cla_comp[cl2]] == false) break;
				}
			}
		}
				
		if(cl2 == sp.ncla) clist_rel.push(k);
	}
	return clist_rel;
}

	
/// Converts from a filter in different classifications to a global filter
function global_convert(filt,sp)
{
	let gfilt = [];
	
	for(let c = 0; c < sp.comp_gl.length; c++){
		let cgl = sp.comp_gl[c];
	
		let cl = 0; while(cl < sp.ncla && filt.cla[cl].comp[cgl.cla_comp[cl]].check == true) cl++;
		
		if(cl == sp.ncla) gfilt.push(true);
		else gfilt.push(false);
	}
	
	return gfilt;
}


/// Gets a specific sample
function get_isel(rpf)
{
	if(rpf.sel_chain.te == "All") return rpf.sel_sample.i;
	else return rpf.sel_sample_chain[rpf.sel_chain.i].i;
}


/// Adds buttons associated with transition plots
function add_trans_buts(res,lay)
{
	if(inter.graph.init == "loading") return;
		
	let cx = corner.x;
	let cy = corner.y;

	cy = lay.add_title("Transition rates",cx,cy,{te:trans_plot_text});
	
	cy += 1;

	if(inter.graph.init != true){
		inter.graph.init = "loading"; 
		start_worker("Graph trans",res_worker(res));
		return;
	}
	
	inter.graph.create(cx,cy,graph_width-right_menu_width,29,lay); 
}


/// Creates a reference between tra and tra_gl (for fast plotting of transitions)
function set_tr_tra_gl_ref(sp)
{
	sp.tr_tra_gl_ref=[];
	for(let cl = 0; cl < sp.cla.length; cl++){
		sp.tr_tra_gl_ref[cl] = []; 
		let claa = sp.cla[cl];
		for(let tr = 0; tr < claa.tra.length; tr++){
			sp.tr_tra_gl_ref[cl][tr] = [];
		}
	}
		
	for(let trg = 0; trg < sp.tra_gl.length; trg++){
		let tra = sp.tra_gl[trg];
		sp.tr_tra_gl_ref[tra.cl][tra.tr].push(trg);
	}
}

	
/// Calculates quantities needed for the graph
function graph_trans_calculate(result,rpf,burn,p)
{
	// Sets classification
	let rpf2 = rpf.species[p];
	let cl = rpf2.sel_class.cl;

	// Incorprates filters
	let popfilt = get_popfilt(rpf2,"Transitions");

	// Sets the chain
	let chsel = rpf.sel_chain.te;
	
	let sample = result.sample;
	
	// Sets samples
	let imin = 0, imax = result.sample.length;
	let isel = get_isel(rpf); if(isel != "All"){ imin = isel; imax = isel+1; burn = 0;}

	let view = rpf.sel_trans_view.te;

	if(result.plot_average || ((view == "Compartment" || view == "Density") && isel == "All" && chsel == "All")){
		imin = 0; imax = 1;
		isel = "All";
		sample = result.average;
	}		
	
	let sp = result.species[p];
	let claa = sp.cla[cl];
	
	if(!sp.tr_tra_gl_ref) set_tr_tra_gl_ref(sp);

	let step = rpf.sel_timestep.i;

	let data = [];

	let key = [];
	
	let line_max = false;
	
	let timepoint = result.timepoint;
	
	let T = timepoint.length-1;
						
	if(view == "Data"){                              // If data is viewed 
		let seldata = rpf.sel_poptransdata;
		let so = result.species[seldata.p].inf_source[seldata.i];
	
		let strat_filt = seldata.strat_filt;

		let cl = find(sp.cla,"name",so.spec.cl_drop.te);
		
		// Works out which transitions apply to data
		let gfilt = [];	
		for(let tr = 0; tr < sp.tra_gl.length; tr++){
			let trg = sp.tra_gl[tr];
			
			let flag = false;
			if(trg.cl == cl){
				let c = trg.i; if(c == SOURCE) c = trg.f;
				
				let cgl = sp.comp_gl[c];
				
				let cl2;
				for(cl2 = 0; cl2 < sp.ncla; cl2++){
					if(cl == cl2){
						if(so.spec.filter.tra[trg.tr].check != true) break;
					}
					else{
						let filt = so.spec.filter.cla[cl2];
						let cc = cgl.cla_comp[cl2];
						
						if(!(filt.radio.value == "All" || (filt.radio.value == "ObsMod" && filt.comp[cc].prob_eqn.te != "0") ||
        				filt.comp[cc].check == true || strat_filt[cl2] == cc)) break;
					}
				}
				if(cl2 == sp.ncla) flag = true;
			}
			gfilt.push(flag);
		}
		
		let line_post = [];
		let line_mean = [];
		for(let t = 0; t < T; t++) line_mean[t] = 0;
		let nline_mean = 0;
		{                                              // Posterior
			let tr_list = [];
			for(let k = 0; k < sp.tra_gl.length; k++){
				if(gfilt[k] == true) tr_list.push(k);
			}
			
			for(let i = imin; i < imax; i++){
				let samp = sample[i];
				if(samp.num >= burn && !(chsel != "All" && samp.chain != chsel)){
					let sa = samp.species[p];
					
					let point = [];
					
					let y=[];
					for(let t = 0; t < T; t++) y[t] = 0;
					
					for(let ij = 0; ij < tr_list.length; ij++){
						add_tl(y,sa.transnum_tl[tr_list[ij]],result);
					}
					
					for(let t = 0; t < T; t++) line_mean[t] += y[t];
					nline_mean++;
					
					for(let t = 0; t < T; t += step){
						let step_max = t + step; if(step_max > T) step_max = T;
						let dt = timepoint[step_max] - timepoint[t];
						let num = 0;
						for(let tt = t; tt < step_max; tt++) num += y[tt];	
				
						point.push({x:timepoint[t] + dt/2, y:num/dt});
					}

					line_post.push(point);
				}
			}
			add_line(view,line_post,BLUE,posterior_name(result),data,key);
			
			for(let t = 0; t < T; t++) line_mean[t] /= nline_mean;
		}
		
		{                                              // Data
			let line = [];
			let tab = so.table;
						
			let spec = so.spec;
			let dist = spec.obs_error_dist.value;
			let ty = spec.obs_error.value;
		
			let col = BLACK;
			
			let point = [];		
			for(let j = 0; j < tab.nrow; j++){
				let t = Number(tab.ele[j][0]);
				let t2 = Number(tab.ele[j][1]);
				let dtdata = t2 - t;
			
				let dfilt=[];
				for(let k = 0; k < so.load_col.length; k++){
					let cl = so.load_col[k].cl;
					if(cl != undefined){
						let co = tab.ele[j][k];
						let c = hash_find(sp.cla[cl].hash_comp,co);
						if(c == undefined) error("Could not find "+co);
						dfilt[cl] = c;					
					}
				}
				
				let cl = 0; while(cl < sp.ncla && dfilt[cl] == strat_filt[cl]) cl++;
				
				if(cl == sp.ncla){	
					let m = find(so.load_col,"heading","Number");
					if(m == undefined) error("Cannot find Number");
					
					let val_data = Number(tab.ele[j][m]);
						
					// Generate error bar on posterior
					let ti = get_ti(t,timepoint);
					let ti2 = get_ti(t2,timepoint);
					if(ti != undefined && ti2 != undefined && ti != ti2){
						let val = ((timepoint[ti+1]-t)/(timepoint[ti+1]-timepoint[ti]))*line_mean[ti];
						for(let ii = ti+1; ii < ti2; ii++){
							if(ii < line_mean.length) val += line_mean[ii];
						}							
						if(ti2 < line_mean.length) val += ((t2-timepoint[ti2])/(timepoint[ti2+1]-timepoint[ti]))*line_mean[ti2];
					
						let sd;
						switch(dist){
						case "Normal": 
							switch(ty){
							case "percent": sd = Number(spec.percent)*val/100; break;
							case "sd": sd = Number(spec.sd); break;
							case "file":
								{
									let m = find(so.load_col,"heading","sd");
									if(m == undefined) error("Cannot find sd");
									sd = Number(tab.ele[j][m]); break;
								}
								break;
							}
							break;
							
						case "Poisson":
							sd = Math.sqrt(val);
							break;
		
						case "Negative binomial":
							{
								let p;
								switch(spec.obs_error_p.value){
								case "p": p = Number(spec.p); break;
								case "file":
									{
										let m = find(so.load_col,"heading","p");
										if(m == undefined) error("Cannot find p");
										p = Number(tab.ele[j][m]); break;
									}
									break;
								}
								sd = Math.sqrt(val/p);
							}
							break;
						}
						
						if(!isNaN(t) && !isNaN(val) && !isNaN(sd)){	
							let ymin = (val-sd)/dtdata; if(ymin < 0) ymin = 0;
							data.push({type:"ErrorBar", x:t+dtdata/2, ymin:ymin, y:val/dtdata, ymax:(val+sd)/dtdata, col:col});
						}
					}
					
					if(!isNaN(t) && !isNaN(val_data)){
						point.push({x:t, y:val_data/dtdata});
						point.push({x:t+dtdata, y:val_data/dtdata});
					}
				}
			}
			line.push(point);
			add_line(view,line,RED,"Data",data,key);
		}
	}
	else{                                            // General transitions are observed
		// Potentially only shows some lines
		let line_filt;
		for(let i = 0; i < rpf2.filter.length; i++){
			let rpf3 = rpf2.filter[i];
			if(rpf3.type == "trans_filt" && rpf3.cl == cl){ line_filt = rpf3.tra_filt; break;}
		}
		
		let col_line = get_col_trans(claa);
		
		let stepp = get_percent_step(claa.ntra);
		for(let j = 0; j < claa.ntra; j++){
			output_percent(stepp,j,claa.ntra);
				
			if(line_filt == undefined || line_filt[j].check == true){
				let tr = claa.tra[j];
				
				let col = col_line[j];
			
				let line = [];
			
				let ttgr = sp.tr_tra_gl_ref[cl][j];
				
				// Makes a list of a transitions to be included
				let tr_list=[];
				for(let ij = 0; ij < ttgr.length; ij++){
					let k = ttgr[ij];
					let tr = sp.tra_gl[k];
			
					let cc = tr.i; if(cc == SOURCE) cc = tr.f;
						
					let co = sp.comp_gl[cc];
						
					let cl2 = 0; 
					while(cl2 < sp.ncla && (popfilt[cl2] == undefined || popfilt[cl2][co.cla_comp[cl2]] == true)) cl2++;
				
					if(cl2 == sp.ncla) tr_list.push(k);
				}
				
				for(let i = imin; i < imax; i++){
					let samp = sample[i];
					if(samp.num >= burn && !(chsel != "All" && samp.chain != chsel)){
						let sa = samp.species[p];
						
						let point = [];
						
						let y=[];
						for(let t = 0; t < T; t++) y[t] = 0;
						
						for(let ij = 0; ij < tr_list.length; ij++){
							add_tl(y,sa.transnum_tl[tr_list[ij]],result);
						}
							
						for(let t = 0; t < T; t += step){
							let step_max = t + step; if(step_max > T) step_max = T;
						
							let dt = timepoint[step_max] - timepoint[t];
							
							let num = 0;
							for(let tt = t; tt < step_max; tt++) num += y[tt]
							
							if(T <= step){
								point.push({x:timepoint[t], y:num/dt});
								point.push({x:timepoint[t] + dt, y:num/dt});
							}
							else{
								point.push({x:timepoint[t] + dt/2, y:num/dt});
							}	
						}							
					
						line.push(point);
					}
				}
				add_line(view,line,col,tr.name,data,key);
				
				if(key.length > KEY_LINE_MAX){ line_max = true; break;}
			}
		}
	}
	
	post({type:"Graph define", variety:"Transition", view:view, data:data, op:{ def_xrange:time_range(result.details), x_label:"Time", y_label:"Transition rate", key:key, line_max:line_max}});
}
 

/// Gets colours for transitions
function get_col_trans(claa)
{
	let col_line = [];
	for(let j = 0; j < claa.ntra; j++){ 
		let tr = claa.tra[j];
		let col_pos;
		if(tr.i == SOURCE || tr.f == SINK) col_pos = BLACK;
		else col_pos = claa.comp[tr.i].col;
		
		if(find_in(col_line,col_pos) == undefined) col_line.push(col_pos);
		else{
			let k = 0; while(k < auto_color.length && find_in(col_line, auto_color[k]) != undefined) k++;
			if(k < auto_color.length) col_line.push(auto_color[k]);
			else col_line.push(auto_color[Math.floor(Math.random()*auto_color.length)]);
		}
	}
	
	return col_line;
}
		
		
/// Time variation in excess events
function graph_trans_expect_calculate(result,rpf,burn,p,type)
{
	// Sets classification
	let rpf2 = rpf.species[p];
	let cl = rpf2.sel_class.cl;

	// Incorprates filters
	let popfilt = get_popfilt(rpf2,"Transitions");

	// Sets the chain
	let chsel = rpf.sel_chain.te;
	
	let sample = result.sample;
	
	// Sets samples
	let imin = 0, imax = result.sample.length;
	let isel = get_isel(rpf); if(isel != "All"){ imin = isel; imax = isel+1; burn = 0;}

	let view = "Graph (CI)";
	
	if(result.plot_average || ((view == "Compartment" || view == "Density") && isel == "All" && chsel == "All")){
		imin = 0; imax = 1;
		isel = "All";
		sample = result.average;
	}		
	
	let sp = result.species[p];
	let claa = sp.cla[cl];
	
	if(!sp.tr_tra_gl_ref) set_tr_tra_gl_ref(sp);

	let step = rpf.sel_timestep.i;

	let data = [];

	let key = [];
	
	let line_max = false;
	
	let timepoint = result.timepoint;
	
	let T = timepoint.length-1;
						
	// Potentially only shows some lines
	let line_filt;
	for(let i = 0; i < rpf2.filter.length; i++){
		if(rpf2.filter[i].cl == cl){ line_filt = rpf2.filter[i].tra_filt; break;}
	}
	
	let col_line = get_col_trans(claa);
	
	let stepp = get_percent_step(claa.ntra);
	for(let j = 0; j < claa.ntra; j++){
		output_percent(stepp,j,claa.ntra);
			
		if(line_filt == undefined || line_filt[j].check == true){
			let tr = claa.tra[j];
			
			let col = col_line[j];
		
			let line = [];
		
			let ttgr = sp.tr_tra_gl_ref[cl][j];
			
			// Makes a list of a transitions to be included
			let tr_list=[];
			for(let ij = 0; ij < ttgr.length; ij++){
				let k = ttgr[ij];
				let tr = sp.tra_gl[k];
		
				let cc = tr.i; if(cc == SOURCE) cc = tr.f;
					
				let co = sp.comp_gl[cc];
					
				let cl2 = 0; 
				while(cl2 < sp.ncla && (popfilt[cl2] == undefined || popfilt[cl2][co.cla_comp[cl2]] == true)) cl2++;
			
				if(cl2 == sp.ncla) tr_list.push(k);
			}
			
			for(let i = imin; i < imax; i++){
				let samp = sample[i];
				if(samp.num >= burn && !(chsel != "All" && samp.chain != chsel)){
					let sa = samp.species[p];
					
					let point = [];
					
					let y=[];
					for(let t = 0; t < T; t++) y[t] = 0;
					
					for(let ij = 0; ij < tr_list.length; ij++){
						add_tl(y,sa.transnum_tl[tr_list[ij]],result);
					}
						
					for(let t = 0; t < T; t += step){
						let step_max = t + step; if(step_max > T) step_max = T;
					
						let dt = timepoint[step_max] - timepoint[t];
						
						let num = 0;
						for(let tt = t; tt < step_max; tt++) num += y[tt]
						
						if(T <= step){
							point.push({x:timepoint[t], y:num/dt});
							point.push({x:timepoint[t] + dt, y:num/dt});
						}
						else{
							point.push({x:timepoint[t] + dt/2, y:num/dt});
						}	
					}							
				
					line.push(point);
				}
			}
		
			add_line(view,line,col,tr.name,data,key);
			
			// Adds in line for expected transition number
			let point=[];
			for(let t = 0; t < T; t += step){
				let step_max = t + step; if(step_max > T) step_max = T;
					
				let dt = timepoint[step_max] - timepoint[t];
						
				let num = 0;
				for(let ij = 0; ij < tr_list.length; ij++){
					let exp_num = sp.exp_num[tr_list[ij]];
					for(let tt = t; tt < step_max; tt++) num += exp_num[tt];
				}
		
				if(T <= step){
					point.push({x:timepoint[t], y:num/dt});
					point.push({x:timepoint[t] + dt, y:num/dt});
				}
				else{
					point.push({x:timepoint[t] + dt/2, y:num/dt});
				}	
			}				
				
			let col_sh = shift_colour(col,0.6);
			data.push({point:point, col:col_sh, view:view, type:"Line", dash:TRANS_EXP_DASH, thick:TRANS_EXP_THICK});
			
		}
	}
	
	key.push({type:"Line", te:"Expected", col:BLACK, dash:TRANS_EXP_DASH, thick:TRANS_EXP_THICK});
	
	post({type:"Graph define", variety:"Transition", view:view, data:data, op:{ def_xrange:time_range(result.details), x_label:"Time", y_label:"Transition rate", key:key, line_max:line_max}});
}


/// Calculates quantities needed for the graph
function graph_trans_hbin_calculate(result,rpf,burn,p,type)
{
	// Sets classification
	let rpf2 = rpf.species[p];
	let cl = rpf2.sel_class.cl;

	let popfilt = get_popfilt(rpf2,"Transitions");
	
	let data = [];

	let key = [];
	
	let sp = result.species[p];
	let claa = sp.cla[cl];
	
	let sample = result.sample;
	
	// Sets the chain
	let chsel = rpf.sel_chain.te;

	// Sets samples
	let imin = 0, imax = sample.length;
	let isel = get_isel(rpf); if(isel != "All"){ imin = isel; imax = isel+1; burn = 0;}

	if(!sp.tr_tra_gl_ref) set_tr_tra_gl_ref(sp);

	switch(type){
	case "Distribution":
		{
			let sel = rpf2.cla[cl].sel;	
			let j = sel.radio.value;
				
			let ttgr = sp.tr_tra_gl_ref[cl][j];
					
			// Makes a list of a transitions to be included
			let tr_list=[];
			for(let ij = 0; ij < ttgr.length; ij++){
				let k = ttgr[ij];
				let tr = sp.tra_gl[k];
		
				let cc = tr.i; if(cc == SOURCE) cc = tr.f;
					
				let co = sp.comp_gl[cc];
					
				let cl2 = 0; 
				while(cl2 < sp.ncla && (popfilt[cl2] == undefined || popfilt[cl2][co.cla_comp[cl2]] == true)) cl2++;
			
				if(cl2 == sp.ncla) tr_list.push(k);
			}
			
		
			let hbin=[];
			for(let b = 0; b < H_BIN; b++) hbin[b]=[];
				
			for(let i = imin; i < imax; i++){
				let samp = sample[i];
				if(samp.num >= burn && !(chsel != "All" && samp.chain != chsel)){
					let sa = samp.species[p];
					for(let j = 0; j < tr_list.length; j++){
						let hbin2 = sa.trans_hbin[tr_list[j]];
						for(let b = 0; b < H_BIN; b++) hbin[b].push(hbin2[b]);
					}
				}
			}
		
			for(let b = 0; b < H_BIN; b++){
				let name = (b/H_BIN)+"-"+((b+1)/H_BIN);
				let stat = get_statistic(hbin[b]);
				data.push({type:"Bar", name:name, x:b+0.5, y:stat.mean, thick:bar_thick, col:GREY});
				data.push({type:"ErrorBar", x:b+0.5, ymin:stat.CImin, y:stat.mean, ymax:stat.CImax, col:BLACK});
			}
			
			post({type:"Graph define", variety:"Histogram", view:"Histogram", data:data, op:{x_label:"Probability distribution", x_param:false, y_label:"Frequency"}});
		}
		break;
	
	case "Histogram":
		{
			let sel = rpf2.cla[cl].sel;	
			let j = sel.radio.value;
				
			// Potentially only shows some lines
			let line_filt;
			for(let i = 0; i < rpf2.filter.length; i++){
				if(rpf2.filter[i].cl == cl){ line_filt = rpf2.filter[i].tra_filt; break;}
			}
		
			let col_line = get_col_trans(claa);
			
			let view = rpf.sel_diag_view.te;
			
			for(let j = 0; j < claa.ntra; j++){	
				if(line_filt == undefined || line_filt[j].check == true){		
					let ttgr = sp.tr_tra_gl_ref[cl][j];
							
					// Makes a list of a transitions to be included
					let tr_list=[];
					for(let ij = 0; ij < ttgr.length; ij++){
						let k = ttgr[ij];
						let tr = sp.tra_gl[k];
				
						let cc = tr.i; if(cc == SOURCE) cc = tr.f;
							
						let co = sp.comp_gl[cc];
							
						let cl2 = 0; 
						while(cl2 < sp.ncla && (popfilt[cl2] == undefined || popfilt[cl2][co.cla_comp[cl2]] == true)) cl2++;
					
						if(cl2 == sp.ncla) tr_list.push(k);
					}
					
					let vec = [];					
					for(let i = imin; i < imax; i++){
						let samp = sample[i];
						if(samp.num >= burn && !(chsel != "All" && samp.chain != chsel)){
							let sa = samp.species[p];
							for(let j = 0; j < tr_list.length; j++){
								let hbin = sa.trans_hbin[tr_list[j]];
								
								let sum = 0, nav = 0;
								
								switch(view){
								case "Trans. (bias)": 
									for(let b = 0; b < H_BIN; b++){
										sum += (2*(b+0.5)/H_BIN - 1)*hbin[b];
										nav += hbin[b];
									}
									if(nav != 0) sum /= nav;
									break;
								
								case "Trans. (p-val.)": 
									{
										for(let b = 0; b < H_BIN; b++) nav += hbin[b];
									
										if(nav != 0){
											let E = nav/H_BIN;
											let xis = 0;
											for(let b = 0; b < H_BIN; b++){
												xis += (hbin[b]-E)*(hbin[b]-E)/E;
											}
					
											sum = -Math.log(xi_sq_p_value(xis));
										}
									}
									break;
								}
								
								vec.push(sum);
							}
						}
					}
				
					let stat = get_statistic(vec);
				
					data.push({type:"Bar", name:claa.tra[j].name, x:j+0.5, y:stat.mean, thick:bar_thick, col:col_line[j]});
					data.push({type:"ErrorBar", x:j+0.5, ymin:stat.CImin, y:stat.mean, ymax:stat.CImax, col:BLACK});
				}
			}
			
			let y_label;
			switch(view){
			case "Trans. (bias)": 
				y_label = "Excess transitions";
				break;
				
			case "Trans. (p-val.)":
				{		
					y_label = "-log(p-value)";
					
					let point=[];
					point.push({x:0, y:-Math.log(0.05)});
					point.push({x:claa.ntra, y:-Math.log(0.05)});
			
					data.push({point:point, col:BLACK, type:"Line", thick:SIM_VALUE_THICK, dash:SIM_VALUE_DASH});
				}
				break;
			}
				
			post({type:"Graph define", variety:"Histogram", view:"Histogram", data:data, op:{x_label:"Transition", x_param:false, y_label:y_label}});
		}
		break;
	
	default:
		error("Type not found");
		break;
	}
}
	
	
/// Adds a line to the graph
function add_line(view,line,col,name,data,key)
{
	switch(view){
	case "Graph": case "Graph (all)": case "Graph (split)": case "Graph (lines)":
		for(let j = 0; j < line.length; j++){
			data.push({point:line[j], col:col, view:view, type:"Line"});
		}
		key.push({type:"Line", te:name, col:col});
		break;
		
	case "Graph (CI)": case "Data":
		{
			let line_stats = get_line_stats(line);
			data.push({point:line_stats, col:col, type:"Line CI"});
			key.push({type:"Line", te:name, col:col});
		}
		break;
	
	default: error("view option error2"); break;
	}
}

			
/// Adds buttons associated with an individual plot
function add_individual_buts(res,lay)
{
	if(inter.graph.init == "loading") return;
	
	let cx = corner.x;
	let cy = corner.y;
	
	let ind_sel = inter.graph.ind_sel;
	if(ind_sel != undefined){
		cy = lay.add_title(ind_sel,cx,cy,{te:ind_sing_plot_text});
		cy += 1;
		
		if(inter.graph.init != true){
			inter.graph.init = "loading"; 
			start_worker("Graph ind sing",{res:res_worker(res), ind_sel:ind_sel});
			return;
		}
		
		inter.graph.create(cx,cy,graph_width-right_menu_width,29,lay);
		
		lay.add_corner_button([["Back","Grey","CancelViewInd"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
		return;
	}
	
	cy = lay.add_title("Individuals",cx,cy,{te:ind_plot_text});
	

	// Sets species
	let rpf = res.plot_filter;	
	let p = model.get_p();
	
	if(rpf.species[p].type == "Population"){
		inter.graph.no_graph(cx,cy,graph_width-right_menu_width,29,lay,"No individual data for a population model.");
		inter.graph.init = "no data";
		return;
	}

	cy += 1;

	if(inter.graph.init != true){
		inter.graph.init = "loading"; 
		start_worker("Graph ind",res_worker(res));
		return;
	}
	
	let vari = plot_variety(inter.graph.type);
	
	switch(vari){
	case "Line plot": inter.graph.create(cx,cy,graph_width-right_menu_width,29,lay); break;
	case "Stat table plot": inter.graph.create(cx,cy,graph_width-right_menu_width,29,lay); break;
	case "Scatter plot": 
		{
			let x = cx, w = graph_width-right_menu_width, h = 29;
			if(false){                                   // Makes square
				let tw = inter.graph.tick.wmax;
				let ww = w-graph_mar.right-graph_mar.left-tw;
				let hh = h-graph_mar.top-graph_mar.bottom;
				if(ww > hh){ 
					let d = (ww-hh)/2;
					x += 0.5*d; w -= 2*d;
				}
			}
			inter.graph.create(x,cy,w,h,lay); 
		}
		break;
	case "Individual plot": inter.graph.create(0,cy,graph_width+2-right_menu_width,29,lay); break;
	case "TransTree plot": inter.graph.create(0,cy,graph_width+2-right_menu_width,29,lay); break;
	case "PhyloTree plot": inter.graph.create(0,cy,graph_width+2-right_menu_width,29,lay); break;
	case "Histogram plot": inter.graph.create(0,cy,graph_width+2-right_menu_width,29,lay); break;
	case "No graph plot": inter.graph.create(0,cy,graph_width+2-right_menu_width,29,lay); break;
	default: error("Graph not plotting"+vari); break;
	}
}


/// Displays diagnostic information
function add_diagnostic_buts(res,lay)
{
	if(inter.graph.init == "loading") return;
	
	let cx = corner.x;
	let cy = corner.y;

	let w = graph_width-right_menu_width;
	
	let rpf = res.plot_filter;
	
	let view = rpf.sel_diag_view.te;

	switch(view){
	case "Proposals":
		cy = lay.add_title("Diagnostics",cx,cy,{te:diagnostic_text});
		break;
		
	case "Trans. (exp.)": 
		cy = lay.add_title("Transition rate vs expected from model",cx,cy,{te:trans_exp_text});
		break;
		
	case "Trans. (dist.)": 
		cy = lay.add_title("Cumulative probability distribution for transitions",cx,cy,{te:trans_dist_text});
		break;
		
	case "Trans. (bias)": 
		cy = lay.add_title("Cumulative probability distribution bias for transitions",cx,cy,{te:trans_bias_text});
		break;
		
	case "Trans. (p-val.)": 
		cy = lay.add_title("Cumulative probability distribution p-value for transitions",cx,cy,{te:trans_pval_text});
		break;
	}
	
	if(view == "Proposals"){
		let x1 = lay.x+cx;
		let x2 = lay.x+lay.dx-cx;
		let y1 = lay.y+cy+0.5;
		let y2 = lay.y+lay.dy-1;

		add_layer("Diagnostics",x1,y1,w-0.5,y2-y1,{te:res.plot_filter.diagnostics[rpf.sel_diag_chain.i+1]});
		return;
	}
	if(inter.graph.init == "loading") return;
	
	cy += 1;

	if(inter.graph.init != true){
		inter.graph.init = "loading"; 
		start_worker("Graph "+view,res_worker(res));
		return;
	}
	
	inter.graph.create(cx,cy,graph_width-right_menu_width,29,lay); 
}


/// Page giving diagnistics information
function diagnostics_buts(lay)
{	
	let te = lay.op.te;
	let lines = te.split("\n");
	
	let cx = 0, cy = 0;
	let mar = 1;
	for(let i = 0; i < lines.length; i++){
		cy = lay.add_paragraph(lines[i],lay.dx-4*mar,cx,cy,BLACK,para_si,para_lh);
	}
}


/// Calculates quantities needed for the graph
function graph_ind_calculate(result,rpf,burn,p)
{
	percent(1);
	
	// Sets classification
	let rpf2 = rpf.species[p];
	let cl = rpf2.sel_class.cl;
	
	// Sets the chain
	let chsel = rpf.sel_chain.te;
	
	// Sets samples
	let imin = 0, imax = result.sample.length;
	let isel = get_isel(rpf); if(isel != "All"){ imin = isel; imax = isel+1; burn = 0;}
	
	let sp = result.species[p];
	
	let claa = sp.cla[cl];

	let tp = result.timepoint;
	let tmin = tp[0], tmax = tp[tp.length-1];

	if(rpf.sel_indview.te == "Trans. Tree"){
		add_trans_tree(imin,imax,tmin,tmax,chsel,result,rpf,burn);
		return;		
	}
	
	if(rpf.sel_indview.te == "Phylo. Tree"){
		add_phylo_tree(imin,imax,tmin,tmax,chsel,result,rpf,burn);
		return;		
	}
	
	if(rpf.sel_indview.te == "First Inf. Time"){
		add_first_inf_time(imin,imax,tmin,tmax,chsel,result,rpf,burn);
		return;		
	}
	
	if(rpf.sel_indview.te == "First Inf. Ind."){
		add_first_inf_ind(imin,imax,tmin,tmax,chsel,result,rpf,burn);
		return;		
	}
	
	let key = [];
	let data = [];

	// Sets up a potential individual filter
	let hash_filt;
	if(rpf2.sel_indgroup.te != "All"){
		let name_list = rpf2.sel_indgroup.name_list;
		hash_filt = new Hash();
		for(let k = 0; k < name_list.length; k++) hash_filt.add(name_list[k],k);
	}

	// Makes a list of individuals from all samples
	let hash = new Hash();

	let ind_max = false;
	
	let ind_list = [];
	for(let samp = imin; samp < imax; samp++){
		let sampl = result.sample[samp];
		if(sampl.num >= burn && !(chsel != "All" && sampl.chain != chsel)){		
			let sa = sampl.species[p];
			let imax = sa.individual.length;
			if(imax > IND_PLOT_MAX){
				imax = IND_PLOT_MAX;
				ind_max = true;
			}
				
			for(let i = 0; i < imax; i++){
				let ind = sa.individual[i];
				
				let name = ind.name;
				if(hash_filt == undefined || hash_filt.find(name) != undefined){
					let j = hash.find(name);
					if(j == undefined){
						j = ind_list.length;
						hash.add(name,ind_list.length);
						ind_list.push({name:name, sa_ref:[], obs:[]});
					}
					ind_list[j].sa_ref.push({sample:samp, i:i});
				}
			}
		}
	}

	let popfilt = get_popfilt(rpf2,"Populations");
		
	/// Filters individuals based on a population filter
	if(popfilt && popfilt.length > 0){
		let g_filt = get_g_filt(popfilt,sp);
		
		let ind_list_new=[];
		
		for(let k = 0; k < ind_list.length; k++){
			let ind_li = ind_list[k];
			
			let fl = false;
			for(let j = 0; j < ind_li.sa_ref.length; j++){
				let sr = ind_li.sa_ref[j];
				
				let ind = result.sample[sr.sample].species[p].individual[sr.i];
				
				let c = ind.cinit;
				if(g_filt[c] == true){ fl = true; break;}
				
				let eve = ind.ev;
						
				for(let e = 0; e < eve.length; e++){
					let ev = eve[e];
						
					switch(ev.type){
					case EV_TRANS: c = sp.tra_gl[ev.trg].f; break;
					case EV_ENTER: c = ev.c; break;	
					case EV_LEAVE: c = OUT; break;		
					case EV_MOVE: c = ev.cf; break;
					default: error("option not recognised"); break;
					}
					
					if(c != OUT && g_filt[c] == true){ fl = true; break;}
				}
			
				if(fl == true) break;
			}
			
			if(fl == true) ind_list_new.push(ind_li);
		}
		
		ind_list = ind_list_new;
		
		// Redos hash for reduced list
		hash = new Hash(); for(let j = 0; j < ind_list.length; j++) hash.add(ind_list[j].name,j);
	}

	percent(10);

	// Looks to add actual values for individual effects from data
	if(result.siminf == "inf" || result.siminf == "ppc"){   
		for(let k = 0; k < sp.inf_source.length; k++){
			let so = sp.inf_source[k];
			if(so.type == "Ind. Eff."){
				let tab = so.table;
				let ie = so.spec.drop.te;
				for(let r = 0; r < tab.nrow; r++){
					let j = hash.find(tab.ele[r][0]);
					if(j != undefined){
						ind_list[j][ie] = tab.ele[r][1];
					}
				}
			}
		}
	}
	
	switch(result.siminf){
	case "sim": 
		add_source_obs(sp.sim_source,claa,hash,result,ind_list);
		break;
	case "inf": 
		add_source_obs(sp.inf_source,claa,hash,result,ind_list);
		break;
	case "ppc": 
		add_source_obs(sp.inf_source,claa,hash,result,ind_list);
		add_source_obs(model.inf_res.plot_filter.species[p].ppc_source,claa,hash,result,ind_list);
		break;
	}
	
	// Genetic data could be places in any species
	if(result.siminf == "inf" || result.siminf == "ppc"){     
		for(let p2 = 0; p2 < result.species.length; p2++){
			let sp2 = result.species[p2];
		
			for(let i = 0; i < sp2.inf_source.length; i++){
				let so = sp2.inf_source[i];

				if(so.type == "Genetic"){
					let tab = so.table;
					
					for(let r = 0; r < tab.nrow; r++){
						let j = hash.find(tab.ele[r][0]);
						if(j != undefined){
							let tstr = get_time(tab.ele[r][1],result);
							let t = Number(tstr);	
							ind_list[j].obs.push({type:"GeneticObs", t:t, tstr:tstr, i:i});
						}
					}
				}
			}
		}
	}
	
	percent(20);
	
	// Gets statistics about individual effects
	for(let i = 0; i < ind_list.length; i++){
		if(i%100 == 0) percent(20+60*i/ind_list.length);
			
		let comb_ind = ind_list[i];
		
		let ie_vec = [], log_ie_vec = [];
		let nie = sp.ind_effect.length;
		for(let e = 0; e < nie; e++){ ie_vec[e] = []; log_ie_vec[e] = [];}
		
		let kmax = comb_ind.sa_ref.length;
		for(let k = 0; k < kmax; k++){
			let sr = comb_ind.sa_ref[k];
			let ind = result.sample[sr.sample].species[p].individual[sr.i];
			if(ind.ie.length != nie) error("ie not the right size");
			for(let e = 0; e < nie; e++){
				ie_vec[e].push(ind.ie[e]);
				log_ie_vec[e].push(Math.log(ind.ie[e]));
			}
		}
		
		comb_ind.ie_stat = []; 
		comb_ind.log_ie_stat = [];
		for(let e = 0; e < nie; e++){
			comb_ind.ie_stat[e] = get_statistic(ie_vec[e]);
			comb_ind.log_ie_stat[e] = get_statistic(log_ie_vec[e]);
		}
	}
	
	percent(80);
	
	// Goes through all individuals from all samples
	switch(rpf.sel_indview.te){
	case "Timeline": add_timelines(ind_list,tmin,tmax,p,cl,result,data,key,ind_max); break
	case "Ind. Eff.":
		switch(rpf.sel_indeffview.te){
		case "Scatter": add_ind_eff_scatter(ind_list,p,result,rpf,data,key,ind_max); break;
		case "Distribution": add_ind_eff_distribution(ind_list,p,result,rpf,data,key,ind_max); break;
		default: error("Option sel ind eff view not here");
		}
		break
	case "Table": add_individual_table(ind_list,p,result,data,ind_max); break;
	default: error("Option not here: "+rpf.sel_indview.te); break;
	}
}


/// Adds symbols from observed data
function add_source_obs(source,claa,hash,result,ind_list)
{
	// Puts observed data on the time line
	for(let i = 0; i < source.length; i++){
		let so = source[i];

		switch(so.type){
		case "Add Ind.":
			if(so.table_loaded == true){
				let tab = so.table;							
				let c = find_in(tab.heading,claa.name);
				if(c == undefined) error("Could not find classification1");
				else{
					for(let r = 0; r < tab.nrow; r++){
						let j = hash.find(tab.ele[r][0]);
						if(j != undefined){
							let tstr = get_time(tab.ele[r][1],result);
							
							let t = Number(tstr);
					
							let list = find_obs_colour(claa,tab.ele[r][c]);
							ind_list[j].obs.push({type:"AddObs", t:t, tstr:tstr, col_list:list, i:i});
						}
					}
				}
			}
			break;
		
		case "Remove Ind.":
			if(so.table_loaded == true){
				let tab = so.table;							
			
				for(let r = 0; r < tab.nrow; r++){
					let j = hash.find(tab.ele[r][0]);
					if(j != undefined){
						let tstr = get_time(tab.ele[r][1],result);
						let t = Number(tstr);
					
						ind_list[j].obs.push({type:"RemObs", t:t, tstr:tstr, i:i});
					}
				}
			}
			break;
			
		case "Move Ind.":
			if(so.table_loaded == true){
				let tab = so.table;						
				let spec = so.spec;
		
				if(claa.name == spec.cl_drop.te){
					for(let r = 0; r < tab.nrow; r++){
						let j = hash.find(tab.ele[r][0]);
						if(j != undefined){
							let tstr = get_time(tab.ele[r][1],result);
							let t = Number(tstr);
					
							let list = find_obs_colour(claa,tab.ele[r][2]);
							ind_list[j].obs.push({type:"MoveObs", t:t, tstr:tstr, col_list:list, i:i});
						}
					}
				}
			}
			break;
			
		case "Transition":
			if(so.table_loaded == true){
				let tab = so.table;
				let spec = so.spec;
		
				if(claa.name == spec.cl_drop.te){
					let pos_tra=[];
					let tra = spec.filter.tra;
				
					for(let k = 0; k < tra.length; k++){
						if(tra[k].check == true) pos_tra.push(tra[k].tra_name_store);
					}
					
					let col = []; 
					if(pos_tra.length == 1){
						let spl = pos_tra[0].split("→");
						for(let loop = 0; loop < 2; loop++){
							let na = spl[loop];
							col[loop] = WHITE;
					
							for(let c = 0; c < claa.ncomp; c++){
								if(claa.comp[c].name == na) col[loop] = claa.comp[c].col;
							}
						}
					}
			
					for(let r = 0; r < tab.nrow; r++){
						let j = hash.find(tab.ele[r][0]);
						if(j != undefined){
							let tstr = get_time(tab.ele[r][1],result);
							let t = Number(tstr);
						
							ind_list[j].obs.push({type:"TransObs", t:t, tstr:tstr, col:col, i:i});
						}
					}
				}
			}
			break;
			
		case "Compartment":
			if(so.table_loaded == true){
				let tab = so.table;
			
				let spec = so.spec;
				if(claa.name == spec.cl_drop.te){
					for(let r = 0; r < tab.nrow; r++){
						let j = hash.find(tab.ele[r][0]);
						if(j != undefined){
							let tstr = get_time(tab.ele[r][1],result);
							let t = Number(tstr);
							let cname = tab.ele[r][2];
							if(j != undefined){
								let list = find_obs_colour(claa,cname);
								ind_list[j].obs.push({type:"CompObs", t:t, tstr:tstr, col_list:list, i:i});
							}
						}
					}
				}
			}
			break;
		
		case "Diag. Test":
			if(so.table_loaded == true){
				let tab = so.table;
			
				let spec = so.spec;
				if(claa.name == spec.cl_drop.te){
					for(let r = 0; r < tab.nrow; r++){
						let j = hash.find(tab.ele[r][0]);
						if(j != undefined){
							let tstr = get_time(tab.ele[r][1],result);
							let t = Number(tstr);
							
							let res_str = tab.ele[r][2];
							let res;
							if(res_str == spec.neg_result) res = false;
							if(res_str == spec.pos_result) res = true;
							
							if(res != undefined && j != undefined){
								ind_list[j].obs.push({type:"DiagObs", t:t, tstr:tstr, res:res, i:i});
							}
						}
					}
				}
			}
			break;
		}
	}
}


/// Gets a global compartment filter from a population filter
function get_g_filt(popfilt,sp)
{
	let g_filt = [];
	for(let c = 0; c < sp.comp_gl.length; c++){	
		let cgl = sp.comp_gl[c];
			
		let cl = 0; 
		while(cl < sp.ncla && (popfilt[cl] == undefined || popfilt[cl][cgl.cla_comp[cl]] == true)) cl++;
		if(cl == sp.ncla) g_filt[c] = true;
		else g_filt[c] = false;
	}
	
	return g_filt;
}


/// Calculates quantities needed for the graph
function graph_ind_sing_calculate(ind_sel,result,rpf,burn,p)
{
	percent(1);

	// Sets classification
	let rpf2 = rpf.species[p];
	let cl = rpf2.sel_class.cl;

	// Sets the chain
	let chsel = rpf.sel_chain.te;
	
	let view = rpf.sel_ind_sel_view.te;
	
	// Sets samples
	let imin = 0, imax = result.sample.length;
	let isel = get_isel(rpf); if(isel != "All"){ imin = isel; imax = isel+1; burn = 0;}
	
	let sp = result.species[p];
	let claa = sp.cla[cl];

	let tp = result.timepoint;
	let tmin = tp[0], tmax = tp[tp.length-1];
	
	let key = [];
	let data = [];

	let line_max = false;
	
	let popfilt = get_popfilt(rpf2,"Populations");
	
	// Constructs a global filter 
	let g_filt = get_g_filt(popfilt,sp);
	
	// Gets the individual number
	let i_store = [];
	let i_ind;
	for(let samp = imin; samp < imax; samp++){
		let sampl = result.sample[samp];
		if(sampl.num >= burn && !(chsel != "All" && sampl.chain != chsel)){		
			let sa = sampl.species[p];
			
			if(i_ind == undefined || i_ind >= sa.individual.length || sa.individual[i_ind].name != ind_sel){ 	
				i_ind = find(sa.individual,"name",ind_sel);
			}
			
			i_store[samp] = i_ind;
		}
	}
		
	let stepp = get_percent_step(claa.ncomp);
	for(let c = 0; c < claa.ncomp; c++){
		output_percent(stepp,c,claa.ncomp);
		
		if(!(popfilt[cl] != undefined && popfilt[cl][c] == false)){
			let line = [];
			
			for(let samp = imin; samp < imax; samp++){
				let sampl = result.sample[samp];
				if(sampl.num >= burn && !(chsel != "All" && sampl.chain != chsel)){		
					let sa = sampl.species[p];
					
					let i = i_store[samp];
					if(i != undefined){
						let ind = sa.individual[i];
						
						let cc = ind.cinit;
					
						let point = [];
					
						let eve = ind.ev;
						
						let e = 0;
						for(let t = 0; t < tp.length; t++){
							let tl = tp[t];
							
							while(e < eve.length && eve[e].t < tl){
								let ev = eve[e];
						
								switch(ev.type){
								case EV_TRANS: cc = sp.tra_gl[ev.trg].f; break;
								case EV_ENTER: cc = ev.c; break;	
								case EV_LEAVE: cc = OUT; break;		
								case EV_MOVE: cc = ev.cf; break;
								default: error("option not recognised"); break;
								}
								e++;
							}
								
							let prob = 0;
							if(cc != SOURCE && cc != SINK && cc != OUT){
								if(sp.comp_gl[cc].cla_comp[cl] == c){
									if(g_filt[cc] == true) prob = 1;
								}
							}								
							point.push({x:tl, y:prob});
						}
						line.push(point);
					}
				}
			}
		
			let col = claa.comp[c].col;
			
			switch(view){	
			case "Graph":
				{
					let line_stats = get_line_stats(line);
		
					data.push({point:line_stats, col:col, type:"Line CI"});
					key.push({type:"Line", te:claa.comp[c].name, col:col, c:c});
					if(key.length > KEY_LINE_MAX){ line_max = true; break;}
				}
				break;
			
			case "Compartment": case "Density":
				{
					let line_stats = get_line_mean(line);
					data.push({y_vec:line_stats, col:col, c:c, type:"Comp Video"});
				}
				break;
			}
		}
		
		if(key.length > KEY_LINE_MAX){ line_max = true; break;}
	}
	
	percent(100);
	
	switch(view){
	case "Graph":
		post({type:"Graph define", variety:"Population", view:"Graph (CI)", data:data, op:{ def_xrange:time_range(result.details), x_label:"Time", y_label:"Probability", key:key, line_max:line_max}});
		break;
		
	case "Compartment": case "Density": 
		post({type:"Graph define", variety:"Population", view:view, data:data, op:{p:p, cl:cl, timepoint:result.timepoint, number_comp:true}});
		break;
	}
}


/// Gets the time from a string (and accounts for "start" and "end" key words)
function get_time(te,result)
{
	if(te == "start") return result.details.t_start;
	if(te == "end") return result.details.t_end;
	return te;	
}


/// Takes "S", or "S|I" or "S:0.7|I:0.5 and converts into a list of colours probabilities
function find_obs_colour(claa,st)
{
	let list = [];

	let flag = false;

	let frac_sum = 0;

	if(st == "."){
		list.push({col:WHITE, frac:1});
	}
	else{
		if(st == "!"){ 
			list.push({col:"notalive", frac:1});
		}
		else{
			let spl = split_with_bracket(st,"|");
			for(let i = 0; i < spl.length; i++){
				let spl2 = spl[i].split(":");
				
				let c = hash_find(claa.hash_comp,spl2[0]);	
				if(c == undefined){ error("Cannot find obs compartment"); return;}

				let frac = 1; 
				if(spl2.length == 2){
					if(isNaN(spl2[1])) flag = true;
					else frac = Number(spl2[1]);
				}
				frac_sum += frac;

				list.push({col:claa.comp[c].col, frac:frac});
			}

			if(flag == false){
				for(let i = 0; i < list.length; i++) list[i].frac /= frac_sum;
			}
			else{
				for(let i = 0; i < list.length; i++) list[i].frac = 1/list.length;	
			}	
		}
	}
	
	return list;
}


/// Adds individual table
function add_individual_table(ind_list,p,result,data,ind_max)
{
	let sp = result.species[p];
	
	let twid = [];
	
	let table = {width:twid, heading:[], content:[]};
	
	for(let i = 0; i < ind_list.length; i++) table.content.push([]);
		
	let vec = [];
	for(let i = 0; i < ind_list.length; i++) vec.push(ind_list[i].name);
	
	add_table_column("ID",vec,table);
	
	percent(85);
	
	let nie = sp.ind_effect.length;
	for(let e = 0; e < nie; e++){	
		let flag = false;
		
		{
			let name = "["+sp.ind_effect[e]+"]";
			
			let vec=[];
			for(let i = 0; i < ind_list.length; i++){
				let stat = ind_list[i].ie_stat[e];
				let te = stat.mean.toPrecision(pre);
				if(stat.CImin != stat.mean || stat.CImax != stat.mean){
					te += " ("+stat.CImin.toPrecision(pre)+" — "+stat.CImax.toPrecision(pre)+")";
				}
				vec.push(te);
			}
			
			add_table_column(name,vec,table);
		}
		
		{
			let name = "log(["+sp.ind_effect[e]+"])";
			let vec=[];
			for(let i = 0; i < ind_list.length; i++){
				let stat = ind_list[i].log_ie_stat[e];
				let te = stat.mean.toPrecision(pre);
				if(stat.CImin != stat.mean || stat.CImax != stat.mean){
					te += " ("+stat.CImin.toPrecision(pre)+" — "+stat.CImax.toPrecision(pre)+")";
				}
				if(ind_list[i][name] != undefined) flag = true;
				vec.push(te);
			}
			
			add_table_column(name,vec,table);
			
			if(flag == true){ // True ie values exist
				let vec=[];
				for(let i = 0; i < ind_list.length; i++){
					if(ind_list[i][name] != undefined) vec.push(ind_list[i][name]);
					else vec.push("-");
				}
				add_table_column("True "+name,vec,table);
			}
		}
	}	
	
	percent(90);
	
	let fe_lookup = create_fe_lookup(sp);
	
	let head = fe_lookup.head;
	for(let f = 0; f < head.length; f++){
		let hash = fe_lookup.hash_fe[f];
		
		let vec=[];
		for(let i = 0; i < ind_list.length; i++){
			let X = hash.find(ind_list[i].name);
			if(X == undefined) X = "-";
			vec.push(X);
		}
		
		add_table_column(head[f],vec,table);
	}
	
	// Outputs individual data observations
	let list = [];
	for(let i = 0; i < ind_list.length; i++){
		let obs = ind_list[i].obs;
		for(let k = 0; k < obs.length; k++){
			if(find_in(list,obs[k].i) == undefined) list.push(obs[k].i);
		}
	}
	
	percent(95);
	
	let mat=[];
	for(let j = 0; j < list.length; j++){
		mat[j] = [];
		for(let i = 0; i < ind_list.length; i++) mat[j][i] = "";
	}
	
	for(let i = 0; i < ind_list.length; i++){
		let obs = ind_list[i].obs;
		for(let k = 0; k < obs.length; k++){
			let ob = obs[k];
			let j = find_in(list,ob.i);
			if(mat[j][i] != "") mat[j][i] += "|"; 
			mat[j][i] += "t="+ob.tstr;
			if(ob.res != undefined){
				if(ob.res == true) mat[j][i] += ":+"; else mat[j][i] += ":-";
			}
		}
	}
	
	for(let j = 0; j < list.length; j++){
		let so;
		switch(result.siminf){
		case "sim": so = sp.sim_source[list[j]]; break;
		case "inf": so = sp.inf_source[list[j]]; break;
		default: break;
		}
		
		if(so != undefined){
			let te = so.name;
			add_table_column(te,mat[j],table);
		}
	}
	
	// Adds padding on the last column (if possible)
	let wid = 0;
	for(let i = 0; i < twid.length; i++) wid += twid[i];
	
	let wid_max = graph_width-right_menu_width-3;
	if(wid < wid_max) twid[twid.length-1] += wid_max-wid;
	
	data.push({type:"Table", table:table});
	
	percent(100);
	
	post({type:"Graph define", variety:"Stat Table", view:"Stat Table", data:data, op:{ind_max:ind_max}});
}


/// Adds a column to a table
function add_table_column(head,vec,table)
{
	let wmax = 0;
	table.heading.push({name:head});
	let w = text_width_table(head,true); if(w > wmax) wmax = w;
	
	for(let j = 0; j < vec.length; j++){
		table.content[j].push({te:vec[j]});
		let w = text_width_table(vec[j]); if(w > wmax) wmax = w;
	}
	wmax += 3;
	
	if(wmax < 5) wmax = 5;
	
	table.width.push(wmax);
}


/// Works out the width of text in a table
function text_width_table(te,bold)
{
	cv.font = "30px Times";
	if(bold == true) cv.font = "bold 30px Times";
	return si_table*cv.measureText(te).width/30;
}
	
	
/// Works out the width of text in a table
function text_width_para(te)
{
	cv.font = "30px arial";
	return para_si*cv.measureText(te).width/30;
}
	
	
/// Adds scatter plot showing true vs estimated individual effect
function add_ind_eff_scatter(ind_list,p,result,rpf,data,key,ind_max)
{		
	let seliev = rpf.sel_ie_data_view[p];
	
	if(seliev.list.length == 0){
		post(no_graph_msg("No individual effects"));
		return;
	}
	
	if(seliev.radio.value >= seliev.list.length) seliev.radio.value = 0;
	
	let name = seliev.list[seliev.radio.value].name;
	
	let ie = name.substr(1,name.length-2);
	
	let ind_effect = result.species[p].ind_effect;
	
	let e = find_in(ind_effect,ie);
	if(e == undefined) alertp("Problem getting ie");
	
	let namex = "True log("+name+")";
					
	let namey = "Estimated log("+name+")";
			
	let vecA=[], vecB=[];
	let point=[];
	for(let i = 0; i < ind_list.length; i++){
		if(i%100 == 0) percent(80+20*i/ind_list.length);
		
		if(ind_list[i][ie] != undefined){
			let val = Math.log(Number(ind_list[i][ie]));
			if(!isNaN(val)){
				let stat = ind_list[i].log_ie_stat[e];
	
				vecA.push(val); vecB.push(stat.mean);
				
				if(rpf.scatter_settings.show_eb.check == true){
					data.push({type:"ErrorBar", x:val, ymin:stat.CImin, y:stat.mean, ymax:stat.CImax, col:BLACK}); 
				}
				
				point.push({x:val, y:stat.mean});
			}
		}	
	}

	data.push({point:point, type:"Points", col:BLUE});
	
	percent(100);
	
	let cor = get_correlation(vecA,vecB);
	let title = "Prediction accuracy = "+cor.toPrecision(pre);
	
	post({type:"Graph define", variety:"Scatter", view:"Scatter", data:data, op:{x_label:namex, x_param:true, y_label:namey, y_param:true, key:[], title:title, ind_max:ind_max}});		
}


/// Shows distribution in individual effect
function add_ind_eff_distribution(ind_list,p,result,rpf,data,key,ind_max)
{	
	let seliev = rpf.sel_ie_view[p];
	
	if(seliev.list.length == 0){
		post(no_graph_msg("No individual effects"));
		return;
	}
	
	if(seliev.radio.value >= seliev.list.length) seliev.radio.value = 0;
	
	let name = seliev.list[seliev.radio.value].name;
	
	let ie = name.substr(1,name.length-2);
	
	let ind_effect = result.species[p].ind_effect;
	
	let e = find_in(ind_effect,ie);
	if(e == undefined) alertp("Problem getting ie");
		
	let vec=[];
	for(let i = 0; i < ind_list.length; i++){
		vec.push(ind_list[i].log_ie_stat[e].mean);
	}
	
	percent(90);
		
	distribution_add_data_line(vec,posterior_name(result),BLUE,data,key,rpf,true);
		
	percent(100);
	
	post({type:"Graph define", variety:"Distribution", view:"Graph", data:data, op:{x_label:"log("+name+")", x_param:true, y_label:"Probability", param:true, key:key, line_max:false, yaxis:false}});
}


/// Creates a look up table to find X fixed effect easily
function create_fe_lookup(sp)
{
	let hash_fe = []; 
	
	let nfe = sp.fix_eff.length; 
	let head = [];
	for(let f = 0; f < nfe; f++){
		let fe = sp.fix_eff[f];
	
		hash_fe[f] = new Hash();
		head.push("X("+fe.name+")");
		
		let X = fe.X_vector;
		let Xind_list = X.ind_list;
		let value = X.X_value;
		for(let k = 0; k < Xind_list.length; k++){
			hash_fe[f].add(Xind_list[k],value[k]);
		}
	}
	
	return { head:head, hash_fe:hash_fe};
}


/// Adds individual timelines
function add_timelines(ind_list,tmin,tmax,p,cl,result,data,key,ind_max)
{
	let sp = result.species[p];
	let claa = sp.cla[cl];
	
	// Works out colours for each of the compartments
	let comp_col = [];
	for(let c = 0; c < claa.comp.length; c++){
		let co = claa.comp[c];
		let col = co.col;
		let bigint = parseInt(col.substring(1), 16);	
		comp_col.push({r:(bigint >> 16) & 255, g:(bigint >> 8) & 255, b:bigint & 255});
		
		key.push({type:"Rect", te:co.name, col:co.col});
	}
	
	let addobs_fl, remobs_fl, moveobs_fl, transobs_fl, compobs_fl, diagobs_fl, geneticobs_fl;
	
	// Adds key for compartmental observations	
	for(let i = 0; i < ind_list.length; i++){
		for(let k = 0; k < ind_list[i].obs.length; k++){
			switch(ind_list[i].obs[k].type){
			case "AddObs": addobs_fl = true; break;
			case "RemObs": remobs_fl = true; break;
			case "MoveObs": moveobs_fl = true; break;
			case "TransObs": transobs_fl = true; break;
			case "CompObs": compobs_fl = true; break;
			case "DiagObs": diagobs_fl = true; break;
			case "GeneticObs": geneticobs_fl = true; break;
			default: error("Option obs not recognised"); break;
			}
		}
	}
	
	if(addobs_fl) key.push({type:"AddObs", te:"Add Ind."});
	if(remobs_fl) key.push({type:"RemObs", te:"Remove Ind."});
	if(moveobs_fl) key.push({type:"MoveObs", te:"Move Ind."});
	if(transobs_fl) key.push({type:"TransObs", te:"Trans. Data"});
	if(compobs_fl) key.push({type:"CompObs", te:"Comp. Data"});
	if(diagobs_fl) key.push({type:"DiagObs", te:"Diag. Tests"});
	if(geneticobs_fl) key.push({type:"GeneticObs", te:"Gen. Data"});
	
	let fe_lookup = create_fe_lookup(sp);

	for(let i = 0; i < ind_list.length; i++){
		if(i%100 == 0) percent(80+20*i/ind_list.length);
		
		let comb_ind = ind_list[i];
		
		// Gets information about transitions			
		let cpop = [];
		for(let c = 0; c < claa.ncomp; c++) cpop[c] = 0; 
		let nalive = 0;
	
		let nsamp = comb_ind.sa_ref.length;
		
		let col_timeline = [];
			
		let cla_ev = [];
		for(let k = 0; k < nsamp; k++){
			let sr = comb_ind.sa_ref[k];
			let ind = result.sample[sr.sample].species[p].individual[sr.i];
		
			let c = OUT;
			
			if(ind.cinit != OUT && ind.cinit != SOURCE){
				nalive++; 
				let cgl = sp.comp_gl[ind.cinit];
				c = cgl.cla_comp[cl];
				cpop[c]++;
			}
			
			for(let e = 0; e < ind.ev.length; e++){
				let ev = ind.ev[e];
			
				switch(ev.type){
				case EV_TRANS:
					let tr = sp.tra_gl[ev.trg];
					if(tr.cl == cl || tr.i == SOURCE || tr.f == SINK){
						let cf;
						if(tr.f == SINK) cf = OUT; 
						else cf = sp.comp_gl[tr.f].cla_comp[cl];	
						cla_ev.push({t:ev.t, ci:c, cf:cf});
						c = cf;
					}
					break;
					
				case EV_ENTER:
					let cf = sp.comp_gl[ev.c].cla_comp[cl];	
					cla_ev.push({t:ev.t, ci:c, cf:cf});
					c = cf;
					break;
					
				case EV_LEAVE:
					cla_ev.push({t:ev.t, ci:c, cf:OUT});
					c = OUT;
					break;
					
				case EV_MOVE:
					if(ev.cl == cl){
						let cf = sp.comp_gl[ev.cf].cla_comp[cl];	
						cla_ev.push({t:ev.t, ci:c, cf:cf});
						c = cf;
					}
					break;
				}
			}
		}
		
		cla_ev.sort( function(a, b){ return a.t-b.t});
		
		if(nalive > 0 || cla_ev.length > 0){
			col_timeline.push({t:tmin, col:get_mix_colour(comp_col,cpop), alpha:nalive/nsamp});
			
			for(let k = 0; k < cla_ev.length; k++){
				let clev = cla_ev[k];	
				let ci = clev.ci; if(ci == OUT) nalive++; else cpop[ci]--;
				let cf = clev.cf; if(cf == OUT) nalive--; else cpop[cf]++;
				
				col_timeline.push({t:clev.t, col:get_mix_colour(comp_col,cpop), alpha:nalive/nsamp});
			}
			col_timeline.push({t:tmax, col:get_mix_colour(comp_col,cpop), alpha:nalive/nsamp});
		}
		
		let info = "";
		
		for(let e = 0; e < comb_ind.ie_stat.length; e++){
			if(e != 0) info += "   ";
			let name = "["+sp.ind_effect[e]+"]";
			if(nsamp == 1) info += name;
			else info += "⟨"+name+"⟩";
			info += "= "+comb_ind.ie_stat[e].mean.toPrecision(pre);
			
			if(comb_ind[name] != undefined){
				info += " (True ="+Number(comb_ind[name]).toPrecision(pre)+")";
			}
		}
		
		let head = fe_lookup.head;
		for(let f = 0; f < head.length; f++){
			let hash = fe_lookup.hash_fe[f];
		
			let X = hash.find(ind_list[i].name);
			if(X != undefined){
				if(info != "") info += "   ";
				info += head[f]+"="+X;
			}
		}
		
		data.push({type:"Individual", name:comb_ind.name, info:info, tmin:tmin, tmax:tmax, col_timeline:col_timeline, obs:comb_ind.obs});	
	}

	percent(100);

	post({type:"Graph define", variety:"Individual", view:"Individual", data:data, op:{ def_xrange:time_range(result.details), x_label:"Time", key:key, ind_max:ind_max}});
}


/// Creates a plot for the transition tree
function add_trans_tree(imin,imax,tmin,tmax,chsel,result,rpf,burn)
{	
	percent(1);

	// Finds which classification defined infection transitions
	let inf_cl=[];
	for(let p = 0; p < result.species.length; p++){
		let sp = result.species[p];
		inf_cl[p] = find(sp.cla,"name",sp.infection_cl.te);
	}
	
	let ind_max = false;
	 
	let hash = new Hash();
	
	let ind_list = [];
	for(let samp = imin; samp < imax; samp++){
		percent(1+9*(samp-imin)/(imax-imin));
		
		let sampl = result.sample[samp];
		if(sampl.num >= burn && !(chsel != "All" && sampl.chain != chsel)){	
			for(let p = 0; p < sampl.species.length; p++){
				let sa = sampl.species[p];
				
				let imax = sa.individual.length;
				if(imax > IND_PLOT_MAX){ imax = IND_PLOT_MAX; ind_max = true;}
				
				for(let i = 0; i < imax; i++){
					let ind = sa.individual[i];
					let name = ind.name;		
					
					let j = hash.find(name);
					if(j == undefined){
						j = ind_list.length;
						hash.add(name,ind_list.length);
						ind_list.push({ name:name, sa_ref:[], p:p, inf_from:[], inft_av:0, ninft_av:0});
					}
					ind_list[j].sa_ref.push({s:samp,i:i});
				}
			}
		}
	}
	
	percent(10);

	let t_start = Number(result.details.t_start);

	// Calcualtes average time when individual is infected
	for(let k = 0; k < ind_list.length; k++){
		percent(10+30*(k/ind_list.length));
			
		let il = ind_list[k];
		let p = il.p;
		let sp = result.species[p];
		let cl = inf_cl[p];
		let comp_cl = sp.cla[cl].comp;
		
		for(let m = 0; m < il.sa_ref.length; m++){
			let sar = il.sa_ref[m];
			let ind = result.sample[sar.s].species[p].individual[sar.i];
				
			let c = ind.cinit;
			if(c != OUT && c != SOURCE && sp.comp_gl[c].infected == true){
				il.inft_av += t_start; il.ninft_av++;
			}
			
			let ev = ind.ev;
			
			for(let e = 0; e < ev.length; e++){
				let eve = ev[e];
				if(eve.infection){
					il.inft_av += eve.t; il.ninft_av++;
				}
			}
		}
	}		

	// Constructs index to order individuals
	for(let k = 0; k < ind_list.length; k++){
		let il = ind_list[k];
		if(il.ninft_av == 0) il.index = LARGE;
		else il.index = il.inft_av/il.ninft_av;
	}
	
	ind_list.sort( function(a, b){ return a.index-b.index});
		
	// Constructs ind_ref
	let ind_ref = [];
	for(let samp = imin; samp < imax; samp++){
		ind_ref[samp] = [];
		for(let p = 0; p < result.species.length; p++){
			ind_ref[samp][p] = [];
		}
	}
	
	for(let k = 0; k < ind_list.length; k++){
		let il = ind_list[k];
		for(let j = 0; j < il.sa_ref.length; j++){
			let sr = il.sa_ref[j];
			ind_ref[sr.s][il.p][sr.i] = k;
		}
	}
	
	if(false){ prr("Ind list"); prr(ind_list); prr(ind_ref);}
	
	// Works out where individual become infected from
	for(let k = 0; k < ind_list.length; k++){
		percent(30+20*(k/ind_list.length));
			
		let il = ind_list[k];
		il.num_inf = 0;
		il.num_out_inf = 0;
		
		let p = il.p;
		let sp = result.species[p];
		let cl = inf_cl[p];
		let comp_cl = sp.cla[cl].comp;
		
		for(let m = 0; m < il.sa_ref.length; m++){
			let sar = il.sa_ref[m];
			let ind = result.sample[sar.s].species[p].individual[sar.i];
		
			let ev = ind.ev;
			
			let inf_fl = false, out_inf_fl = false;
			for(let e = 0; e < ev.length; e++){
				let eve = ev[e];
				if(eve.infection){
					inf_fl = true;
					let pp = eve.infection.p;
					if(pp == OUTSIDE_INF) out_inf_fl = true;
					else{
						let kk = ind_ref[sar.s][pp][eve.infection.i];
						
						// POTENTIALLY SLOW
						let j =	find(il.inf_from,"k",kk);
					
						if(j == undefined) il.inf_from.push({k:kk, num:1});
						else il.inf_from[j].num++;
					}
				}
			}
			
			if(inf_fl == true) il.num_inf++;
			if(out_inf_fl == true) il.num_out_inf++;
		}
	}	
	
	let transmission = [];
	
	for(let k = 0; k < ind_list.length; k++){
		percent(50+20*(k/ind_list.length));
		
		let il = ind_list[k];
		
		if(il.index != LARGE){
			let kmin = k, kmax = k;
			let num = il.sa_ref.length;
			let inf_list=[];
			for(let j = 0; j < il.inf_from.length; j++){
				let kk = il.inf_from[j].k;
				if(kk < kmin) kmin = kk;
				if(kk > kmax) kmax = kk;
				inf_list.push({k:kk, frac:il.inf_from[j].num/num});
			}
			if(inf_list.length == 0) transmission.push({ k:k, frac:il.num_inf/num, frac_out:il.num_out_inf/num, t:il.index});
			else transmission.push({ k:k, frac:il.num_inf/num, frac_out:il.num_out_inf/num, t:il.index, min:kmin, max:kmax, inf_list:inf_list});
		}
	}
	
	let key = [];
	let data = [];
	
	percent(70);
	
	// Generates timelines for each individual showing when they are infectious
	let tp = result.timepoint;
	let T = tp.length-1;
	for(let k = 0; k < ind_list.length; k++){
		percent(70+30*(k/ind_list.length));
		
		let il = ind_list[k];
		
		let inf_line=[];
		for(let ti = 0; ti < T; ti++) inf_line[ti] = 0;
		let nav = 0;
		
		let sr = il.sa_ref;
		let p = il.p;
		let i = il.i;
		let sp = result.species[p];
		let cl = inf_cl[p];
		let tra_cl = sp.cla[cl].tra;
		let comp_cl = sp.cla[cl].comp;
				
		let tra_gl = result.species[p].tra_gl;
		
		for(let j = 0; j < sr.length; j++){
			let sar = sr[j];
			
			let ind = result.sample[sar.s].species[p].individual[sar.i];
		
			let ev = ind.ev;
			
			let infec = false;
				
			let c = ind.cinit;
			if(c != OUT && c != SOURCE && sp.comp_gl[c].infected == true) infec = true;
			
			if(infec == true || ev.length > 0){
				let e = 0;
			
				let inf_flag = false;
				for(let ti = 0; ti < T; ti++){
					let t = tmin+(tmax-tmin)*(ti+1)/T;
					while(e < ev.length && ev[e].t < t){
						let eve = ev[e];
						let infe = eve.infection;
						if(infe) infec = true;
						else{
							if(infec == true){
								let trg = sp.tra_gl[ev[e].trg];
								if(trg.cl == cl){
									if(trg.f == OUT || trg.f == SINK) infec = false;
									else{
										if(comp_cl[tra_cl[trg.tr].f].infected.check == false) infec = false;
									}
								}
							}
						}
						e++;
					}
					if(infec == true){
						inf_line[ti]++;
						inf_flag = true;
					}
				}	
				if(inf_flag == true) nav++;
			}
		}
		
		if(nav != 0){
			for(let ti = 0; ti < T; ti++) inf_line[ti] /= nav;
		}
		
		data.push({type:"InfBar", name:il.name, tmin:tmin, tmax:tmax, y:k, inf_line:inf_line});
	}
	
	data.push({type:"Transmission", transmission:transmission});		
	
	percent(100);
	
	post({type:"Graph define", variety:"TransTree", view:"TransTree", data:data, op:{ def_xrange:time_range(result.details), x_label:"Time", key:key, ind_max:ind_max}});
}


/// Gets a list of first infected individuals
function get_first_infected_time(sampl)
{
	let tmin = LARGE;
	for(let p = 0; p < sampl.species.length; p++){
		let sa = sampl.species[p];

		for(let i = 0; i < sa.individual.length; i++){
			let ind = sa.individual[i];
	
			let ev = ind.ev;
			for(let e = 0; e < ev.length; e++){
				let eve = ev[e];
				if(eve.infection){
					if(eve.t < tmin) tmin = eve.t;
				}
			}
		}
	}
	
	if(tmin == LARGE) return;
	
	return tmin;
}


function get_first_infected_list(tmin,sampl)	
{
	let fi = [];
	
	for(let p = 0; p < sampl.species.length; p++){
		let sa = sampl.species[p];

		for(let i = 0; i < sa.individual.length; i++){
			let ind = sa.individual[i];
	
			let ev = ind.ev;
			for(let e = 0; e < ev.length; e++){
				let eve = ev[e];
				if(eve.infection){
					if(eve.t == tmin){
						fi.push({p:p, i:i, e:e});
					}
				}
			}
		}
	}

	return fi;
}


/// Creates a distribution plot for the first infection time
function add_first_inf_time(imin,imax,tmin,tmax,chsel,result,rpf,burn)
{	
	let vec=[];
	for(let samp = imin; samp < imax; samp++){
		let sampl = result.sample[samp];
		if(sampl.num >= burn && !(chsel != "All" && sampl.chain != chsel)){	
			let tmin = get_first_infected_time(sampl);
			if(tmin != undefined) vec.push(tmin);
		}						
	}
	
	if(vec.length == 0){
		post(no_graph_msg("No first infection"));
		return;
	}
		
	let key = [];
	let data = [];
	
	let show_mean = rpf.dist_settings.show_mean.check;
	let line_max = false;
	
	distribution_add_data_line(vec,"",BLUE,data,key,rpf,show_mean);
	
	post({type:"Graph define", variety:"Distribution", view:"Graph", data:data, op:{x_label:"Time", x_param:true, y_label:"Probability", param:false, key:key, line_max:line_max, yaxis:false}});
}


/// Creates a distribution plot for the first infection time
function add_first_inf_ind(imin,imax,tmin,tmax,chsel,result,rpf,burn)
{	
	let fii_all = [];

	let sampl = result.sample[0];
		
	{	
		for(let p = 0; p < sampl.species.length; p++){
			fii_all[p]=[];
			
			let sa = sampl.species[p];
			for(let i = 0; i < sa.individual.length; i++){
				fii_all[p][i]=[];
			}			
		}
	}
	
	let num = 0;
	for(let samp = imin; samp < imax; samp++){
		let sampl = result.sample[samp];
		if(sampl.num >= burn && !(chsel != "All" && sampl.chain != chsel)){	
			let tmin = get_first_infected_time(sampl);
			if(tmin != undefined){
				let fii=[];
				for(let p = 0; p < sampl.species.length; p++){
					fii[p]=[];
			
					let sa = sampl.species[p];
					for(let i = 0; i < sa.individual.length; i++) fii[p][i] = 0;
				}
				
				let fil = get_first_infected_list(tmin,sampl);
				
				let fac = 1.0/fil.length;
				
				for(let k = 0; k < fil.length; k++){
					let val = fil[k];
					fii[val.p][val.i] += fac;
				}
				num++;
				
				for(let p = 0; p < sampl.species.length; p++){
					let sa = sampl.species[p];
					for(let i = 0; i < sa.individual.length; i++){
						fii_all[p][i].push(fii[p][i]);
					}
				}
			}
		}						
	}
	
	if(num == 0){
		post(no_graph_msg("No first infection"));
		return;
	}
	
	let data=[];
	
	let x = 0;
	for(let p = 0; p < sampl.species.length; p++){
		let sa = sampl.species[p];
		for(let i = 0; i < sa.individual.length; i++){
			let ind = sa.individual[i];
			
			let fa = fii_all[p][i];
			let j = 0; while(j < fa.length && fa[j] == 0) j++;
			if(j < fa.length){
				let stat = get_statistic(fa);
				
				data.push({type:"Bar", name:ind.name, x:x+0.5, y:stat.mean, thick:bar_thick, col:DGREY});
				data.push({type:"ErrorBar", x:x+0.5, ymin:stat.CImin, y:stat.mean, ymax:stat.CImax, col:BLACK});
				x++;
			}
		}
	}
		
	post({type:"Graph define", variety:"Histogram", view:"Histogram", data:data, op:{x_label:"Individual", x_param:false, y_label:"Probability"}});
}


/// Creates a plot for the phylogenetic tree
function add_phylo_tree(imin,imax,tmin,tmax,chsel,result,rpf,burn)
{	
	if(!result.genetic_data || result.genetic_data.length == 0){
		post(no_graph_msg("No genetic data"));
		return;
	}
	
	let Nobs = result.genetic_data.length;

	percent(1);

	let hash_ori = new Hash();
	hash_ori.set_phylo_num(Nobs); 
	
	let hash_bra = new Hash();
	hash_bra.set_phylo_num(Nobs); 
	
	let node_max = false;
	
	for(let samp = imin; samp < imax; samp++){
		percent(1+89*(samp-imin)/(imax-imin));
		
		let sampl = result.sample[samp];
	
		if(sampl.num >= burn && !(chsel != "All" && sampl.chain != chsel)){	
			let inf_node = copy(sampl.inf_node);
			
			if(false){ prr("ind_nod"); prr(inf_node);}
			
			let nmax = inf_node.length;
			if(nmax > NODE_PLOT_MAX){ no_graph_msg("Too many nodes to plot"); return;}
			
			// Back propagates observatation up the branches
			for(let n = 0; n < nmax; n++){
				let ino = inf_node[n];
				ino.obs_list_root = [];
				for(let k = 0; k < ino.inf_ev.length; k++){
					let iev = ino.inf_ev[k];
					iev.obs_list = [];
				}
			}				
			
			for(let n = 0; n < nmax; n++){
				let ino = inf_node[n];
				for(let k = 0; k < ino.inf_ev.length; k++){
					let iev = ino.inf_ev[k];
					if(iev.ty == "OBS"){
						let m = iev.i; 
					
						let nn = n, kk = k;
						do{
							let ino2 = inf_node[nn];
				
							ino2.inf_ev[kk].obs_list.push(m);
							
							if(kk > 0) kk--;
							else{
								ino2.obs_list_root.push(m);
								
								let fr = ino2.from;
								if(fr == "ENT" || fr == "OUT") break;
								else{
									nn = fr; kk = ino2.num;
								}
							}
						}while(true);
					}
				}
			}
			
			// Finds code for initial state
			let code_list=[];
			let arr = [];
			let t_vec = [];                             // Times when individuals enter
			let mnum_vec = []; 
			
			for(let n = 0; n < nmax; n++){
				let ino = inf_node[n];
				
				let fr = ino.from;
				if(fr == "ENT" || fr == "OUT"){
					code_list.push(hash_ori.get_phylo_code(ino.obs_list_root));
					arr.push(ino.obs_list_root);
					t_vec.push(ino.t);
					mnum_vec.push(ino.num);
				}
			}
			code_list.sort();
			
			hash_ori.add_phylo_ori(JSON.stringify(code_list),t_vec,mnum_vec,arr);
			
			// Adds branching point to hash table
			for(let n = 0; n < nmax; n++){
				let ino = inf_node[n];
				let ind = ino.ind;
				
				for(let k = 0; k < ino.inf_ev.length; k++){
					let iev = ino.inf_ev[k];
					let obs_list = iev.obs_list;
					let code_tot = hash_bra.get_phylo_code(obs_list);
					let ob;
						
					switch(iev.ty){
					case "OBS":
						{
							ob = {ty:"OBS", t:iev.t, m:iev.i, mnum:iev.num, ind:ind};
						}	
						break;
						
					case "INF":
						{
							let list_br = inf_node[iev.i].obs_list_root;
							let list_br2 = [];
							let map=[];
							for(let kk = 0; kk < list_br.length; kk++) map[list_br[kk]] = true;
							for(let kk = 0; kk < obs_list.length; kk++){
								let m = obs_list[kk];
								if(map[m] != true) list_br2.push(m);
							}
							
							let code_br = hash_bra.get_phylo_code(list_br);	
							ob = {ty:"SPLIT", t:iev.t, code_br:code_br, list_br:list_br, list_br2:list_br2, mnum:iev.num, ind:ind};
						}
						break;
					}				

					if(false){
						prr("");
						prr("ADD BRANCH");
						prr(String(code_tot));
						prr(obs_list)
						prr(ob);
					}
					
					hash_bra.add_phylo_br(String(code_tot),ob);
				}
			}
		}
	}

	// Samples the most frequent phylogentic tree 
	let freq_ori = hash_ori.get_phylo_ori()

	let pnode = [];
	
	if(false){
		prr("freq");
		prr(freq_ori);
		prr("HASH");
		prr(hash_ori);
		prr(hash_bra);
	}
	
	for(let i = 0; i < freq_ori.arr.length; i++){
		percent(90+10*(i)/(freq_ori.arr.length));
		
		let n = pnode.length;
		pnode.push({ty:"ROOT", t:freq_ori.t_vec[i], mnum:freq_ori.mnum_vec[i], obs:[], ind_list:[]});
		
		add_phylo_branch(n,freq_ori.arr[i],pnode,hash_bra);
	}
	
	// Sets individual labels
	for(let i = 0; i < pnode.length; i++){
		let ind_list = pnode[i].ind_list;
		
		let max = 0, tot = 0;
		let iname;
		for(let j = 0; j < ind_list.length; j++){
			if(ind_list[j].num > max){ iname = ind_list[j].ind; max = ind_list[j].num;}
			tot += ind_list[j].num;
		}		
		if(max != tot) iname += "("+Math.floor((max*100)/tot)+"%)";

		pnode[i].label = iname;
	}
	
	// Works out the y positions for the various nodes
	
	let y = 1;                                      // Sets Tips
	for(let i = 0; i < pnode.length; i++){
		let pn = pnode[i];
		if(pn.branch == undefined){ pn.y = y; y++;}
	}
	
	for(let i = pnode.length-1; i >= 0; i--){
		let pn = pnode[i];
		if(pn.branch != undefined){
			let yav = 0;
			for(let k = 0; k < pn.branch.br.length; k++){	
				let y = pnode[pn.branch.br[k]].y;
				if(y == undefined) error("y is undefined");
				yav += y;
			}
			pn.y = yav/pn.branch.br.length;
		}		
	}
	
	let key = [];
	key.push({type:"GeneticObs", te:"Gen. Data"});
		
	let data = [];
	data.push({type:"Phylo", pnode:pnode});		
	
	percent(100);
	
	post({type:"Graph define", variety:"PhyloTree", view:"PhyloTree", data:data, op:{ def_xrange:time_range(result.details), x_label:"Time", key:key}});
}


/// Adds information to a phylogentic branch
function add_phylo_branch(n,obs_list,pnode,hash_bra)
{
	let po = pnode[n];
	
	let tnow = po.t; if(po.obs.length > 0) tnow = po.obs[po.obs.length-1].t;
	
	let br = hash_bra.get_phylo_br(obs_list,tnow);

	if(br != undefined){
		for(let k = 0; k < br.ind_list.length; k++){
			let indl = br.ind_list[k];
		
			// POTENTIALLY SLOW
			let j = find(pnode[n].ind_list,"ind",indl.ind);
			if(j == undefined){
				pnode[n].ind_list.push(indl);
			}
			else{
				pnode[n].ind_list[j].num += indl.num;
			}	
		}
	
		switch(br.ty){
		case "OBS":
			{
				pnode[n].obs.push({ m:br.m, mnum:br.mnum_sum/br.num, t:br.t_sum/br.num});
		
				let ii = 0; while(ii < obs_list.length && obs_list[ii] != br.m) ii++;
				if(ii == obs_list.length) error("Cannot find ons");
				if(ii+1 < obs_list.length) obs_list[ii] = obs_list[obs_list.length-1];
				obs_list.pop();
			
				if(obs_list.length > 0){
					add_phylo_branch(n,obs_list,pnode,hash_bra);
				}
			}
			break;
		
		case "SPLIT":
			{
				let t = br.t_sum/br.num;
				
				pnode[n].branch = {mnum:br.mnum_sum/br.num, frac:br.frac, t:t, br:[]};
				
				if(br.list_br.length > 0){
					let nn = pnode.length;
					pnode[n].branch.br.push(nn);
					pnode.push({ty:"INT", t:t, back:n, obs:[], ind_list:[]});
					add_phylo_branch(nn,br.list_br,pnode,hash_bra);
				}
				
				if(br.list_br2.length > 0){
					let nn = pnode.length;
					pnode[n].branch.br.push(nn);
					pnode.push({ty:"INT", t:t, back:n, obs:[], ind_list:[]});
					add_phylo_branch(nn,br.list_br2,pnode,hash_bra);
				}
			}
			break;
		}
	}
	
	return br;
}


/// Mixes the colours from the compartments to get the 
function get_mix_colour(comp_col,cpop)
{
	let rsum = 0, gsum = 0, bsum = 0, sum = 0; 
	for(let c = 0; c < comp_col.length; c++){
		let num = cpop[c];
		if(num != 0){
			let cc = comp_col[c];
			rsum += cc.r*num; gsum += cc.g*num; bsum += cc.b*num;
			sum += num;
		}
	}
	if(sum == 0) return "255,255,255";
	
	return Math.floor(rsum/sum)+","+Math.floor(gsum/sum)+","+Math.floor(bsum/sum);
}


/*
/// Adds buttons associated with a spline plot
function add_spline_buts(res,lay)
{
	if(inter.graph.init == "loading") return;
	
	let cx = corner.x;
	let cy = corner.y;
	
	cy = lay.add_title("Time varying derived quantities",cx,cy,{te:derived_plot_text});
	
	cy += 1;

	if(inter.graph.init != true){
		inter.graph.init = "loading"; 
		start_worker("Graph spline",res_worker(res));
		return;
	}
	
	inter.graph.create(cx,cy,graph_width-right_menu_width,29,lay); 
}
*/


/// Calculates spline information
function graph_spline_calculate(result,rpf,burn)
{
	let spl = rpf.sel_spline.i; 

	// Sets samples
	let imin = 0, imax = result.sample.length;
	let isel = rpf.sel_sample.i;
	
	if(isel != "All"){ imin = isel; imax = isel+1; burn = 0;}

	let view = rpf.sel_spline_view.te;
	
	let data = [];
	
	let key = [];
	
	let spline = result.spline[spl];
	
	let par = result.param[spline.th];

	let col_list = get_col_list(par.dep,result,spline.co_list.length);
	
	let line_ref = [];
	
	for(let i = 0; i < spline.co_list.length; i++){  // There is potential for filtering
		line_ref.push({index:spline.co_list[i].index, col:col_list[i], name:spline.co_list[i].name});
	}
	
	let line_max = false;
	
	let stepp = get_percent_step(line_ref.length);
	for(let j = 0; j < line_ref.length; j++){
		output_percent(stepp,j,line_ref.length);
		
		let lref = line_ref[j];
		let col = lref.col;
		
		let line = [];
	
		for(let i = imin; i < imax; i++){
			let samp = result.sample[i];
			if(samp.num >= burn){				
				let sa = get_element(samp.spline[spl],lref.index);
			
				let point = [];
				
				for(let t = 0; t < sa.length; t++){
					let tl = sa[t];
					point.push({x:tl.t, y:tl.value});
				}
				
				line.push(point);
			}
		}
		
		add_line_to_data(view,line_ref[j].name,col,line,key,data);
		if(key.length > KEY_LINE_MAX){ line_max = true; break;}
	}

	post({type:"Graph define", variety:"Population", view:view, data:data, op:{ def_xrange:time_range(result.details), x_label:"Time", y_label:spline.name, y_param:true, key:key, line_max:line_max}});
}

		
/// Based on a given dependency calculate a colour list for lines
function get_col_list(dep,mod,N)
{
	let col_list = []; 

	// If only one index then uses the colour scheme from the compartments
	if(dep.length == 1 && dep[0] != "t"){    
		let indname = dep[0];
		
		let cla = get_cla_from_index(mod,dep[0])
		if(cla == undefined) error("Could not find classification2");
		else{
			let comp = cla.comp;
			if(comp.length != N) error("col_list problem");
			for(let k = 0; k < comp.length; k++) col_list[k] = comp[k].col;
		}
	}
	else{
		for(let k = 0; k < N; k++) col_list[k] = get_colour(k);
	}
	
	return col_list;
}


/// Returns a colour following a sequence
function get_colour(k)
{
	return auto_color[k%auto_color.length];
}

		
/// Gets the classification from a given index
function get_cla_from_index(mod,index)
{
	for(let p = 0; p < mod.species.length; p++){
		let sp = mod.species[p];
		let cl = find(sp.cla,"index",index);
		if(cl != undefined) return sp.cla[cl];
	}
	return;
}


/// Adds a line to the data
function add_line_to_data(view,name,col,line,key,data)
{
	switch(view){
	case "Graph": case "Graph (all)": case "Graph (split)": case "Graph (lines)": 
		for(let k = 0; k < line.length; k++){
			data.push({point:line[k], col:col, view:view, type:"Line"});
		}
		if(name != "") key.push({type:"Line", te:name, col:col});
		break;
		
	case "Graph (CI)": case "Data":
		{
			let line_stats = get_line_stats(line);
			data.push({point:line_stats, col:col, type:"Line CI"});
			if(name != "") key.push({type:"Line", te:name, col:col});
		}
		break;

	default: error("view option error3"); break;
	}
}


/// Used for testin different graph functions		
function test_graph(lay)
{
	let cx = corner.x;
	let cy = corner.y;
	
	if(false){
		cy = lay.add_title("Histogram",cx,cy,{te:""});
		if(inter.graph.init != true){
			inter.graph.init = true;
			let data = [];
			data.push({type:"Bar", name:"A", x:0.5, y:3, thick:0.7, col:RED});
			data.push({type:"ErrorBar", x:0.5, ymin:2, y:3, ymax:5, col:BLACK});
			data.push({type:"Bar", name:"B", x:1.5, y:10, thick:0.7, col:GREEN});
			
			inter.graph.define("Histogram","Histogram",data,{x_label:"Category", y_label:"Value"});
		}
		
		inter.graph.create(cx,cy,graph_width-right_menu_width,29,lay);
	}
	
	cy = lay.add_title("Matrix",cx,cy,{te:""});
	if(inter.graph.init != true){
		inter.graph.init = true;
		let data = [];
		
		let mat = []; mat[0] = []; mat[1] = []; 
		mat[0][0] = 1; mat[0][1] = 2;
		mat[1][0] = 5; mat[1][1] = 1;
		
		let xlab = ["Cat 1","Cat 2"];
		let ylab = ["Cadcsdcsdcsdt A","Casdcsdcsdcsdct B"];
		
		data.push({type:"Matrix", xlab:xlab, ylab:ylab, mat:mat});
		
		inter.graph.define("Matrix","Matrix",data,{x_label:"Category", y_label:"Value"});
	}
	
	inter.graph.create(cx,cy,graph_width-right_menu_width,29,lay);
}


/// Allows the user to view some graph specified in inter.view_graph
function create_view_graph(lay)
{
	let vg = inter.view_graph;
	let gr = inter.graph;
	
	let so = vg.source; 
	let par = vg.par;
	let sel_view = vg.sel_view;
	let value = vg.value;
	
	lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, type:"Rect", val:WHITE});
	
	let warn = "";
	
	let cx = corner.x;
	let cy = corner.y;
	
	cy = lay.add_title(sel_view.title,cx,cy,{});
	
	if(gr.init != true){
		inter.graph.init = "loading"; 
		
		start_worker("Graph view",{name:par.name, sel_view:sel_view, so:so});
		return;
	}
	
	inter.graph.create(cx,cy,graph_width-right_menu_width,30,lay); 
	
	lay.add_corner_button([["Back","Grey","CancelViewGraph"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
}


/// Calculates a graph view
function create_view_graph_calculate(name,sel_view,so)
{
	let i = 0; while(i < model.param.length && name != model.param[i].name) i++;
	
	if(i == model.param.length) alertp("Could not plot");
		
	let par = model.param[i];
	
	let det = so.sim_details;
	if(par.type == "param factor") det = so.ppc_details;
	
	post(define_parameter_plot("view_graph",par,par.value,undefined,undefined,sel_view,det,so));
}


/// Based on spline information constructs over time
function construct_spline_timevariation(par,value,details)
{
	if(value == undefined) return err("No value");

	let knot = par.spline.knot;

	let res = get_times(knot,details);
	if(res.err) return err(res.msg);
	
	let times = res.times;
		
	let tp = get_time_points(details);

	let val_spl = copy(value);
	let type = par.spline.spline_radio.value;

	switch(par.dep.length){
	case 1: 
		val_spl = get_time_variation(val_spl,times,tp,type); 
		break;
		
	case 2: 
		for(let j = 0; j < val_spl.length; j++){
			val_spl[j] = get_time_variation(val_spl[j],times,tp,type);
		}
		break;
		
	case 3:
		for(let j = 0; j < val_spl.length; j++){
			for(let i = 0; i < val_spl[j].length; i++){
				val_spl[j][i] = get_time_variation(val_spl[j][i],times,tp,type);
			}
		}
		break;
	}
	
	return { err:false, value:val_spl};
}


/// Given value this defines a parameter plot
function define_parameter_plot(from,par,value,CImin,CImax,sel_view,details,so,rpf) 
{
	let view = sel_view.te;

	if(graph_dia) prr("GRAPH DIA: define_param_plot: "+sel_view.type+",  "+view);

	let data = [];
	let key = [];
	
	switch(sel_view.type){
	case "Timevary":
		{			
			switch(view){
			case "Graph": case "Graph (all)": case "Graph (split)":
				{
					let knot = par.spline.knot;
					let res = get_times(knot,details);
					
					if(res.err) return no_graph_msg(res.msg);
						
					let times = res.times;
				
					let tp = get_time_points(details);
					
					let spl_type = par.spline.spline_radio.value;
					
					if(typeof times == 'string'){ alertp(times); close_view_graph(); return;}
					
					let ntimes = times.length;
				
					let list = par_find_list(par);
				
					if(list[list.length-1].length != ntimes) error("Should be same length");
					list.pop();
					
					let co_list = generate_co_list(list);
				
					let col_list = get_col_list(par.dep,so,co_list.length);

					let sim_flag = false;
				
					let imin = 0, imax = co_list.length;
					if(view == "Graph (split)" && rpf){
						imin = rpf.sel_paramview.radio_split.value; imax = imin+1;
					}
					
					let line_max = false;
					
					for(let i = imin; i < imax; i++){
						let line = [];
						let point = [];
						
						let index = co_list[i].index;
						
						let vec_val = get_element(value,index);
						
						for(let k = 0; k < vec_val.length; k++){
							if(vec_val[k] == undefined){
								return no_graph_msg("Not all spline values are specified");
							}
						}
						
						let type = "Line";
						if(CImin){
							type = "Line CI";
							let vec_CImin = get_element(CImin,index);
							let vec_CImax = get_element(CImax,index);
							
							for(let t = 0; t < tp.length; t++){
								point.push({x:tp[t], y:vec_val[t], CImin:vec_CImin[t], CImax:vec_CImax[t]});
							}
						}
						else{	
							switch(from){
							case "param_calc":
								for(let t = 0; t < tp.length; t++){
									point.push({x:tp[t], y:vec_val[t]});
								}
								break;
								
							case "view_graph":
								switch(spl_type){
								case "Square":
									{
										let tend;
										if(!isNaN(details.t_end)){
											let num =  Number(details.t_end);
											if(num > times[times.length-1]) tend = num;
										}
										for(let t = 0; t < ntimes; t++){
											point.push({x:times[t], y:vec_val[t]});
											if(t+1 < ntimes) point.push({x:times[t+1], y:vec_val[t]});
											else{
												if(tend != undefined) point.push({x:tend, y:vec_val[t]});
											}
										}
									}
									break;
									
								case "Linear":
									for(let t = 0; t < ntimes; t++){
										point.push({x:times[t], y:vec_val[t]});
									}
									break;
									
								case "Cubic": case "Cubic +ve":
									{
										let dt;
										let tsta = times[0], tend = times[times.length-1];
										if(!isNaN(details.timestep)) dt = Number(details.timestep);
										else dt = (tend-tsta)/100;
										
										let cspl = solve_cubic_spline(times,vec_val,spl_type);
										if(typeof cspl == "string"){
											return no_graph_msg(cspl);
										}
										
										for(let t = tsta; t <= tend; t += dt){
											point.push({x:t, y:calculate_cubic_spline(t,cspl)});
										}
									}
									break;
								}
								break;
								
							default: error("Option not here:"+from); break;
							}
						}
						
						let name = "";
						for(let d = 0; d < index.length; d++){
							if(d != 0) name += ","; name += list[d][index[d]];
						}
			
						let col = col_list[i];
						data.push({point:point, col:col, type:type});
						if(name != "" && view != "Graph (split)"){
							key.push({type:"Line", te:name, col:col});
							if(key.length > KEY_LINE_MAX){ line_max = true; break;}
						}
						
						// Adds in simuation values
						if(par.variety != "const" && rpf && rpf.sim_val.check == true && rpf.siminf != "sim" && par.type != "derive_param"){	
							let res = construct_spline_timevariation(par,par.value,so.details);
							if(!res.err){ 
								let vec_val = get_element(res.value,index);
								
								let po = [];
								for(let t = 0; t < tp.length; t++){
									po.push({x:tp[t], y:vec_val[t]});
								}
								
								data.push({point:po, col:SIM_VALUE_COL, type:"Line", thick:SIM_VALUE_THICK, dash:SIM_VALUE_DASH});
								
								sim_flag = true;
							}
						}
					}
					
					if(sim_flag == true){
						key.push({type:"Line", te:"Sim. Val.", dash:SIM_VALUE_DASH, thick:SIM_VALUE_THICK, col:SIM_VALUE_COL});
					}
					
					return {type:"Graph define", variety:"Population", view:"Graph", data:data, op:{ def_xrange:time_range(details), x_label:"Time", y_label:par.name+"(t)", y_param:true, key:key, line_max:line_max}};
				}
				break;
		
			case "Matrix":
				{
					let list = par_find_list(par);
					let vecy = list[0];
					let vecx = list[1];

					if(vecy.length*vecx.length > MATRIX_PLOT_MAX){
						return no_graph_msg("Matrix too large");
					}
					
					let val = value;
					if(!CImin){
						let res = construct_spline_timevariation(par,par.value,details);
						if(!res.err) val = res.value;
					}
					
					if(so.siminf == "sim" || par.variety == "const"){
						data.push({type:"Matrix", xlab:vecx, ylab:vecy, mat:val});
					}
					else{
						data.push({type:"Matrix", xlab:vecx, ylab:vecy, mat:val, CImin:CImin, CImax:CImax});
					}
					
					return {type:"Graph define", variety:"MatrixAnim", view:"MatrixAnim", data:data, op:{x_label:par.dep[1], x_param:false, y_label:par.dep[0], y_param:false, italic:true}};
				}
				break;
				
			case "Compartment": case "Density": case "Histogram":
				{
					let p = sel_view.p;
					let cl = sel_view.cl;
					let claa = so.species[p].cla[cl];
					let comp = claa.comp;
					
					let val = value;
					
					if(!CImin){
						let res = construct_spline_timevariation(par,par.value,details);
						if(!res.err) val = res.value;
					}
					
					if(par.dep.length == 3){
						if(comp.length > MATRIX_COMP_MAX){
							return no_graph_msg("Too many lines to plot");
						}
						
						data.push({type:"CompMatrixAnim", value:val});
		
						return {type:"Graph define", variety:"CompMatrixAnim", view:"CompMatrixAnim", data:data, op:{p:p, cl:cl}};
					}
					else{
						for(let c = 0; c < comp.length; c++){
							let point = [];
							for(let i = 0; i < val[c].length; i++){
								if(CImin){
									point.push({ x:i, y:val[c][i], CImin:CImin[c][i], CImax:CImax[c][i], col:comp[c].col});
								}
								else{
									point.push({ x:i, y:val[c][i], col:comp[c].col});
								}
							}
							data.push({type:"SplineAnim", c:c, point:point});
						}
						
						switch(view){
						case "Compartment":	case "Density":
							return {type:"Graph define", variety:"Population", view:view, data:data, op:{p:p, cl:cl}};
							
						case "Histogram":
							if(comp.length > HISTO_PLOT_MAX){
								return no_graph_msg("Too many bars to plot");
							}
						
							return {type:"Graph define", variety:"HistoAnim", view:"HistoAnim", data:data, op:{x_label:claa.name, x_param:false, y_label:"Value", p:p, cl:cl}};
						}
					}
				}
				break;
			}
		}
		break;
		
	case "Vector":
		{
			switch(view){
			case "Histogram":
				{
					let list = par_find_list(par);
					let vec = list[0];
					if(vec.length > HISTO_PLOT_MAX){
						return no_graph_msg("Too many bars to plot");
					}
					
					let col_list = get_col_list(par.dep,so,vec.length);
					
					let sim_val_on = false;
					if(par.variety != "const" && rpf && rpf.siminf != "sim" && rpf.sim_val.check == true && 
						par.value && par.value.length > 0){
						sim_val_on = true;
					}
					
					for(let i = 0; i < vec.length; i++){
						let ind = [i]; 
						
						data.push({type:"Bar", name:par.list[0][i], x:i+0.5, y:value[i], thick:bar_thick, col:col_list[i]});
						
						if(CImin != undefined){
							data.push({type:"ErrorBar", x:i+0.5, ymin:CImin[i], y:value[i], ymax:CImax[i], col:BLACK});
						}
						
						if(sim_val_on){
							let val = Number(get_element(par.value,ind));	
							
							let point = [];
							let mar = 0.1;
							point.push({x:i+mar, y:val}); point.push({x:i+1-mar, y:val});
								
							data.push({point:point, col:SIM_VALUE_COL, type:"Line", thick:SIM_VALUE_THICK, dash:SIM_VALUE_DASH});
						}
					}
					
					if(sim_val_on){
						key.push({type:"Line", te:"Sim. Val.", dash:SIM_VALUE_DASH, thick:SIM_VALUE_THICK, col:SIM_VALUE_COL});
					}
					
					let claa = get_cla_from_index(so,par.dep[0]);
				
					return {type:"Graph define", variety:"Histogram", view:"Histogram", data:data, op:{x_label:claa.name, x_param:false, y_label:"Value", key:key}};
				}
				break;
				
			case "Compartment": case "Density":
				{
					let p = sel_view.p;
					let cl = sel_view.cl;
					let comp = so.species[p].cla[cl].comp;
					
					for(let c = 0; c < comp.length; c++){
						let point = [];
						if(CImin) point.push({ x:0, y:value[c], CImin:CImin[c], CImax:CImax[c]});
						else point.push({ x:0, y:value[c]});
						data.push({type:"Vector", c:c, point:point});
					}
				
					return {type:"Graph define", variety:"CompVector", view:view, data:data, op:{p:p, cl:cl}};
				}
				break;
			}
		}
		break;
		
	case "Matrix":
		{
			switch(view){
			case "Matrix":
				{
					let list = par_find_list(par);
					let vecy = list[0];
					let vecx = list[1];
			
					if(vecy.length*vecx.length > MATRIX_PLOT_MAX){
						return no_graph_msg("Matrix too large");
					}
					
					data.push({type:"Matrix", xlab:vecx, ylab:vecy, mat:value});
		
					return {type:"Graph define", variety:"Matrix", view:"Matrix", data:data, op:{x_label:par.dep[1], x_param:false, y_label:par.dep[0], y_param:false, italic:true}};
				}
				break;
				
			case "Compartment":
				{
					let p = sel_view.p;
					let cl = sel_view.cl;
					let claa = so.species[p].cla[cl];
					let comp = claa.comp;
					
					if(comp.length > MATRIX_COMP_MAX){
						return no_graph_msg("Too many lines to plot");
					}
				
					data.push({type:"CompMatrix", value:value});
		
					return {type:"Graph define", variety:"CompMatrix", view:"CompMatrix", data:data, op:{p:p, cl:cl}};
				}
				break;
				
			default: error("View not recognised"); break;
			}
		}
		break;
		
	default: error("Graph view option not recognised"); break;
	}
}


/// Get an equally spaced timeline from start to end time 
function get_time_points(details)
{
	let ti = Number(details.t_start);
	let tf = Number(details.t_end);
	let dt = Number(details.timestep);
	let tp = [];
	
	for(let t = ti; t < tf; t += dt) tp.push(t);
	tp.push(tf);
	return tp;
}


/// Converts splines from times at know to time on global time line
function get_time_variation(val,times,tp,type)
{
	let value = [];
	
	if(type == "Cubic" || type == "Cubic +ve"){	
		let cspl = solve_cubic_spline(times,val,type);
		if(typeof cspl == "string"){
			return no_graph_msg(cspl);
		}
																			
		for(let i = 0; i < tp.length; i++){
			value.push(calculate_cubic_spline(tp[i],cspl));
		}
	}
	else{
		let j = 0;
		for(let i = 0; i < tp.length; i++){
			let t = tp[i];
			while(j < times.length-1 && times[j+1] < t) j++;
			
			if(j == times.length-1){                      // Extended beyond end time (for PCC)
				value.push(val[j]);
			}
			else{
				switch(type){
				case "Square": 
					value.push(val[j]);
					break;
				
				case "Linear":
					{
						let frac = (t-times[j])/(times[j+1]-times[j]);
						value.push(val[j]*(1-frac) + val[j+1]*frac);
					}
					break;
		
				default: error("Spline type to do"); break;
				}
			}
		}
	}
	return value;
}


/// Closes the view graph
function close_view_graph()
{
	inter.view_graph = {};
	inter.graph.init = undefined;
	inter.graph.type = undefined;
}


/// Closes the view individual
function close_view_ind()
{
	inter.graph.ind_sel = undefined;
	inter.graph.init = undefined;
	inter.graph.type = undefined;
}


/// Adds buttons for the parameters page
function add_parameter_buts(res,lay)
{
	if(inter.graph.init == "loading") return;
	
	let cx = corner.x;
	let cy = corner.y;
	
	cy = lay.add_title("Parameters",cx,cy,{te:param_plot_text});
	cy += 1;
	
	let rpf = res.plot_filter;
	let pview = rpf.sel_paramview;

	if(pview == undefined){
		center_message("No parameters are in the model.",lay);
		return;
	}
	
	if(inter.graph.init != true){
		inter.graph.init = "loading"; 
		
		if(pview.list.length > 0){
			let li = pview.list[0];
			
			let par = res.param[li.th];
			rpf.pos_paramviewtype = [];
			let pvt = rpf.pos_paramviewtype;
			
			let posv = get_par_pos_view(par,rpf);
			for(let i = 0; i < posv.length; i++){
				posv[i].createviewflag = true;
				pvt.push(posv[i]);
			}
					
			if(res.siminf == "inf"){
				switch(rpf.details.algorithm.value){
				case "DA-MCMC": case "PAS-MCMC": pvt.push({te:"Trace"}); break;
				default: pvt.push({te:"Samples"}); break;
				}
			}
			
			pvt.push({te:"Distribution"});	
				
			pvt.push({te:"Statistics"});
		
			// Works out if to keep the current view
			if(!rpf.sel_paramviewtype) rpf.sel_paramviewtype = copy(pvt[0]);
			else{
				let k = find(rpf.pos_paramviewtype,"te",rpf.sel_paramviewtype.te);
				if(k == undefined) rpf.sel_paramviewtype = copy(pvt[0]);
				else rpf.sel_paramviewtype = copy(rpf.pos_paramviewtype[k]);
			}
		}

		start_worker("Graph param",res_worker(res));
		return;
	}
	
	let x = cx, y = cy;
	let w = graph_width-right_menu_width;
	let h = 29;
	
	inter.graph.create(x,y,w,h,lay);
}


/// Adds buttons for the generations page
function add_generation_buts(res,lay)
{
	if(inter.graph.init == "loading") return;
	
	let cx = corner.x;
	let cy = corner.y;
	
	cy = lay.add_title("Generations",cx,cy,{te:generations_text});
	cy += 1;
	
	let rpf = res.plot_filter;
	let pview = rpf.sel_genview;

	if(pview == undefined){
		center_message("No parameters are in the model.",lay);
		return;
	}
	
	if(inter.graph.init != true){
		inter.graph.init = "loading"; 
		
		start_worker("Graph generation",res_worker(res));
		return;
	}
	
	let x = cx, y = cy;
	let w = graph_width-right_menu_width;
	let h = 29;
	
	inter.graph.create(x,y,w,h,lay);
}


/// Returns only those parts of res used in worker
function res_worker(res)
{
	return {siminf:res.siminf, plot_filter:res.plot_filter, burnin:res.burnin, p:model.get_p("no error")};
}


/// Gets a selected parameter
function get_radio_param(ax)
{
	let axpl = ax.param_list[ax.radio.value];
	if(axpl && axpl.tpl_ref != undefined) return axpl.tpl_ref;
	
	return axpl.list[axpl.radio.value].tpl_ref;
}

			
/// Initialises parameter plots
function graph_param_calculate(result,rpf,burn)
{	
	let data = [];
	let pview = rpf.sel_paramview;
	
	if(graph_dia) prr("GRAPH DIA  graph_param_calculate: pview:"+pview.te);
	
	if(pview.param == "too big"){
		post(no_graph_msg("This quantity is too large to output.\nThe threshold number of tensor elements can be altered under 'Further options'."));
		return;
	}

	let pviewtype = rpf.sel_paramviewtype;

	switch(pview.te){
	case "Scatter":
		{
			let sp = rpf.sel_paramview;
		
			let valx = get_radio_param(rpf.xaxis_radio);
			
			let thx = sp.list[valx].th;
			let indx = sp.list[valx].index;
			let parx = result.param[thx];
			let namex = param_name_index(parx,indx);
			
			let valy = get_radio_param(rpf.yaxis_radio);
			
			let thy = sp.list[valy].th;
			let indy = sp.list[valy].index;
			let pary = result.param[thy];
			let namey = param_name_index(pary,indy);
			
			let chmin = 0, chmax = result.chains.length;
			if(rpf.sel_chain.te != "All"){
				chmin = find_in(result.chains,rpf.sel_chain.te);
				chmax = chmin+1;
			}
			
			let point_max = false;
			
			let point=[];	
			for(let ch = chmin; ch < chmax; ch++){
				let cha = result.chains[ch];
				
				for(let i = 0; i < result.par_sample.length; i++){
					let samp = result.par_sample[i];
					if(samp.num >= burn && samp.chain == cha){
						let valx = get_param_val(indx,samp.param[thx],parx);
						let valy = get_param_val(indy,samp.param[thy],pary);
						point.push({x:valx, y:valy});
						if(point.length > POINT_PLOT_MAX){ point_max = true; break;};
					}
				}
				if(point_max) break;
			}
			data.push({point:point, type:"Points", col:BLUE});
		
			if(rpf.sim_val.check == true){
				let parx = result.param[thx];
				if(parx.value && parx.value.length > 0){
					let vax = get_element(parx.value,indx);
					if(vax != undefined && !isNaN(vax)){ 
						let valx = Number(vax);	
						data.push({x:valx, type:"SimX", col:SIM_VALUE_COL, dash:SIM_VALUE_DASH });
					}
				}
				
				let pary = result.param[thy];
				if(pary.value && pary.value.length > 0){
					let vay = get_element(pary.value,indy);
					
					if(vay != undefined && !isNaN(vay)){
						let valy = Number(vay);	
						data.push({y:valy, type:"SimY", col:SIM_VALUE_COL, dash:SIM_VALUE_DASH });
					}
				}
			}

			post({type:"Graph define", variety:"Scatter", view:"Scatter", data:data, op:{x_label:namex, x_param:true, y_label:namey, y_param:true, key:[], point_max:point_max}});
			return;
		}
		break;
		
	case "Correlation":
		{	
			// Works out which parameter to get correlation for
			let list_plot = [];
			for(let th = 0; th < rpf.param_check.length; th++){
				let pc = rpf.param_check[th];

				if(pc.checkb.check){
					if(pc.tpl_ref != undefined) list_plot.push(pc.tpl_ref);
					else{					
						for(let k = 0; k < pc.list.length; k++){
							let pcl = pc.list[k];
							if(pcl.checkb.check) list_plot.push(pcl.tpl_ref);
						}
					}
				}
			}
		
			let len = list_plot.length;
			
			if(len == 0){
				post(no_graph_msg("No parameters have been selected"));	
				return;
			}
			
			if(len*len > MATRIX_PLOT_MAX){
				post(no_graph_msg("Too many elements to be plotted"));	
				return;
			}
			
			let mat = get_correlation_matrix(list_plot,result,rpf,burn);
		
			let vec = [];
			for(let k = 0; k < list_plot.length; k++){
				vec.push(result.total_param_list[list_plot[k]].name);
			}
			
			data.push({type:"Correlation", xlab:vec, ylab:vec, x_param:true, y_param:true, mat:mat});

			post({type:"Graph define", variety:"Matrix", view:"Matrix", data:data, op:{x_label:"Parameters", y_label:"Parameters"}});
			return;
		}
		break;
		
	default:
		{
			if(pviewtype.createviewflag == true){       // Sets up graphs for a multivariate parameter
				multivariate_param_plot(result,rpf,burn);
				return;
			}
			
			let vte = rpf.sel_paramviewtype.te;
	
			switch(vte){
			case "Trace": case "Samples":
				{
					let val = rpf.sel_paramview.radio.value;
					let th = rpf.sel_paramview.list[val].th;
					let ind = rpf.sel_paramview.list[val].index;
				
					let par = result.param[th];
				
					let name = param_name_index(result.param[th],ind);
					
					let chmin = 0, chmax = result.chains.length;
					if(rpf.sel_chain.te != "All"){
						chmin = find_in(result.chains,rpf.sel_chain.te);
						chmax = chmin+1;
					}
					
					let key = [];
					for(let ch = chmin; ch < chmax; ch++){
						let cha = result.chains[ch];
						
						let point=[];
						
						let sample = result.par_sample;
						
						if(par.type == "derive_param" && par.time_dep) sample = result.sample;
						
						let imin = 0;
						while(imin < sample.length && sample[imin].chain != cha) imin++;
						
						let iburn = imin;
						
						if(vte == "Trace"){
							while(iburn < sample.length && sample[iburn].num < burn) iburn++;
							
							if(iburn > imin){
								let point=[];
								for(let i = imin; i <= iburn; i++){
									let samp = sample[i];
									if(samp.chain == cha){
										point.push({x:samp.num, y:get_param_val(ind,samp.param[th],par)});
									}
								}
						
								let col = light_colour(get_colour(ch-chmin));
						
								data.push({point:point, col:col, type:"BurninLine"});
							}
						}
						
						for(let i = iburn; i < sample.length; i++){
							let samp = sample[i];
							if(samp.chain == cha){
								point.push({x:samp.num, y:get_param_val(ind,samp.param[th],par)});
							}
						}
						
						let col = get_colour(ch-chmin);
						
						data.push({point:point, col:col, type:"Line"});
					
						let lab = "Sample"; if(result.chains.length > 1) lab = "Chain "+cha;
						key.push({type:"Line", te:lab, col:col});
					}
					
					if(vte == "Trace"){
						data.push({te:"Burn-in", col:BLACK, thick:1, x:burn, type:"VertLine"});
					}
					
					if(rpf.sim_val.check == true){
						let end = result.par_sample[result.par_sample.length-1].num;
						
						if(par.variety != "const" && par.value && par.value.length > 0){
							let va = get_element(par.value,ind);
							
							if(va != undefined && !isNaN(va)){
								let val = Number(va);	
							
								let point = [];
								point.push({x:0, y:val}); point.push({x:end, y:val});
								
								data.push({point:point, col:SIM_VALUE_COL, type:"Line", thick:SIM_VALUE_THICK, dash:SIM_VALUE_DASH});
								
								key.push({type:"Line", te:"Sim. Val.", dash:SIM_VALUE_DASH, thick:SIM_VALUE_THICK, col:SIM_VALUE_COL});
							}
						}
					}
					
					if(key.length == 1) key = [];
			
					let xax = "# Updates";
					if(vte == "Samples") xax = "Sample"; 
					
					post({type:"Graph define", variety:vte, view:"Graph", data:data, op:{x_label:xax, y_label:name, y_param:true, key:key}});
					return;
				}
				break;
			
			case "Distribution":
				setup_distribution(result,rpf,burn);
				return;
				
			case "Statistics":
				{	
					let table_width, table;
					
					let param_max = false;
					
					let imax = pview.list.length;
					if(imax > PARAM_STAT_MAX){ imax = PARAM_STAT_MAX; param_max = true;}
						
					switch(result.siminf){
					case "sim": case "ppc":
						table_width = [10,10,9,9];
						table = {width:table_width, heading:[{name:"Parameter"},{name:"Mean"},{name:"CI min"},{name:"CI max"}], content:[]};

						for(let i = 0; i < imax; i++){
							let li = pview.list[i];
						
							let stat = get_param_stats(li.th,li.index,result,rpf,burn);
							table.content.push([{te:stat.name, pname:true},{te:stat.mean},{te:stat.CImin},{te:stat.CImax}]);
						}
						break;
						
					case "inf":
						table_width = [8,8,7,7,5,5];
						table = {width:table_width, heading:[{name:"Parameter"},{name:"Mean"},{name:"CI min"},{name:"CI max"},{name:"ESS"},{name:"GR"}], content:[]};

						for(let i = 0; i < imax; i++){
							let li = pview.list[i];
						
							let stat = get_param_stats(li.th,li.index,result,rpf,burn);
							table.content.push([{te:stat.name, pname:true},{te:stat.mean},{te:stat.CImin},{te:stat.CImax},{te:stat.ESS},{te:stat.GR}]);
						}
						break;
					}
					
					data.push({type:"Table", table:table});
			
					post({type:"Graph define", variety:"Stat Table", view:"Stat Table", data:data, op:{param_max:param_max}});
					return;
				}
				break;
			}
		}
		break;
	}
	
	throw({type:"Error", te:"Cannot find graph2"});
}


/// Initialises generation plots
function graph_generation_calculate(result,rpf,burn)
{	
	let data = [];
	
	let pview = rpf.sel_genview;
	
	if(graph_dia) prr("GRAPH DIA  graph_generation_calculate: pview:"+pview.te);
	
	if(pview.param == "too big"){
		post(no_graph_msg("This quantity is too large to output.\nThe threshold number of tensor elements can be altered under 'Further options'."));
		return;
	}
	
	let val = pview.radio.value;
	let th = pview.list[val].th;
	let ind = pview.list[val].index;
				
	let par = result.param[th];
				
	let name = param_name_index(result.param[th],ind);
					
	let key = [];
		
	let pmean=[], pCImin=[], pCImax=[];
	for(let g = 0; g < result.generation.length; g++){
		let gen = result.generation[g];
	
		pmean.push({x:g+0.5, y:get_element(gen.mean[th],ind)});
		pCImin.push({x:g+0.5, y:get_element(gen.CImin[th],ind)});
		pCImax.push({x:g+0.5, y:get_element(gen.CImax[th],ind)});
	}
										
	data.push({point:pmean, col:GENRATION_COL, type:"Line"});
	data.push({point:pCImin, col:GENRATION_CI_COL, type:"Line", dash:GENRATION_DASH});
	data.push({point:pCImax, col:GENRATION_CI_COL, type:"Line", dash:GENRATION_DASH});			
	key.push({type:"Line", te:"Mean", thick:GENRATION_THICK, col:GENRATION_COL});
	key.push({type:"Line", te:"95% CI", thick:GENRATION_THICK, dash:GENRATION_DASH, col:GENRATION_CI_COL});
	
	let sample = result.par_sample;
	if(par.type == "derive_param" && par.time_dep) sample = result.sample;
	
	let vec = [];
	for(let i = 0; i < sample.length; i++){
		let samp = sample[i];
		if(samp.num >= burn){
			vec.push(get_param_val(ind,samp.param[th],par));
		}
	}
	
	let end = result.generation[result.generation.length-1].num+1;
		
	if(vec.length > 0){
		let stat = get_statistic(vec);
	
		data.push({type:"ErrorBar", x:end, ymin:stat.CImin, y:stat.mean, ymax:stat.CImax, col:BLACK});
		key.push({type:"ErrBar", te:posterior_name(result), thick:GENRATION_THICK, col:BLACK});
	}
	
	if(rpf.sim_val.check == true){
		if(par.variety != "const" && par.value){
			let va = get_element(par.value,ind);
			
			if(va != undefined && !isNaN(va)){
				let val = Number(va);	
			
				let point = [];
				point.push({x:0, y:val}); point.push({x:end, y:val});
				
				data.push({point:point, col:SIM_VALUE_COL, type:"Line", thick:SIM_VALUE_THICK, dash:SIM_VALUE_DASH});
				
				key.push({type:"Line", te:"Sim. Val.", dash:SIM_VALUE_DASH, thick:SIM_VALUE_THICK, col:SIM_VALUE_COL});
			}
		}
	}
	
	post({type:"Graph define", variety:"Samples", view:"Graph", data:data, op:{x_label:"Generations", y_label:name, y_param:true, key:key}});
}


/// Sets up graphs for a multivariate parameter
function multivariate_param_plot(result,rpf,burn)
{
	let pview = rpf.sel_paramview;

	let pviewtype = rpf.sel_paramviewtype;
	
	if(graph_dia) prr("GRAPH DIA: multivariate");
	
	let li = pview.list[0];
		
	let par = result.param[li.th];
	
	// Sets the chain
	let chsel = rpf.sel_chain.te;

	let value, CImin, CImax;

	let samp_val_list = []

	if(par.type == "derive_param" && par.time_dep){
		for(let s = 0; s < result.sample.length; s++){
			let samp = result.sample[s];
			if(samp.num >= burn && !(chsel != "All" && samp.chain != chsel)){
				let samp_val = samp.param[li.th];	
				samp_val_list.push(samp_val);
			}
		}
	}
	else{
		for(let s = 0; s < result.par_sample.length; s++){
			let samp = result.par_sample[s];
			if(samp.num >= burn && !(chsel != "All" && samp.chain != chsel)){
				let samp_val = samp.param[li.th];
				if(samp_val == "const") samp_val = result.param[li.th].value;
			
				if(pviewtype.type == "Timevary"){ 
					let res = construct_spline_timevariation(par,samp_val,result.details);
					
					if(res.err){ post(no_graph_msg(res.msg)); return;}
					samp_val = res.value; 
				}
				
				samp_val_list.push(samp_val);
			}
		}
	}
	
	let nsamp = samp_val_list.length;
	if(nsamp == 0){
		throw({type:"Error", te:"Cannot find graph1"});
	}
	
	let sam = samp_val_list[0];
	value = copy(sam);
	CImin = copy(sam);
	CImax = copy(sam);

	let dim = get_dimensions(sam);
	let list = get_element_list(sam,dim);
	for(let i = 0; i < list.length; i++){
		let index = list[i];
	
		let vec = []; 
		for(let s = 0; s < nsamp; s++){
			vec.push(get_element(samp_val_list[s],index));
		}
		
		let stat = get_statistic(vec);
		set_element(value,index,stat.mean);
		set_element(CImin,index,stat.CImin);
		set_element(CImax,index,stat.CImax);
	}
	//}

	post(define_parameter_plot("param_calc",par,value,CImin,CImax,pviewtype,result.details,result,rpf));
}
			
			
/// Gets a matrix giving the correlations between parameters
function get_correlation_matrix(list_plot,result,rpf,burn)
{
	let list = [];
	for(let j = 0; j < list_plot.length; j++){
		list.push(result.total_param_list[list_plot[j]]);
	}
	let N = list.length;
	
	let chmin = 0, chmax = result.chains.length;
	if(rpf.sel_chain.te != "All"){
		chmin = find_in(result.chains,rpf.sel_chain.te);
		chmax = chmin+1;
	}
	
	let value = [], trace = [];
	for(let j = 0; j < N; j++){
		value[j] = [];
		trace[j] = [];
	}

	for(let ch = chmin; ch < chmax; ch++){
		let cha = result.chains[ch];
		
		let point=[];
		for(let i = 0; i < result.par_sample.length; i++){
			let samp = result.par_sample[i];
			if(samp.num >= burn && samp.chain == cha){
				for(let j = 0; j < N; j++){
					let th = list[j].th;
					trace[j].push(get_param_val(list[j].index,samp.param[th],result.param[th]));
				}
			}
		}
	}

	// Transforms such that trace has mean zero and variance 1;
	let K = trace[0].length;
	
	let unvar=[];
	
	for(let j = 0; j < N; j++){
		let sum = 0.0, sum2 = 0.0;
		for(let k = 0; k < K; k++){
			sum += trace[j][k]; sum2 += trace[j][k]*trace[j][k];
		}
	
		let mean = sum/K;
		let vari = sum2/K - mean*mean; if(vari < 0) vari = 0;
		if(vari < TINY) unvar[j] = true;
		else{
			let sd = Math.sqrt(vari);
			for(let k = 0; k < K; k++) trace[j][k] = (trace[j][k]-mean)/sd;
		}
	}
	
	for(let j = 0; j < N; j++){
		for(let i = j; i < N; i++){
			if(unvar[j] == true || unvar[i] == true){
				value[j][i] = 0;
				value[i][j] = 0;
			}
			else{
				let cor = 0.0;
				for(let k = 0; k < K; k++) cor += trace[j][k]*trace[i][k];
				cor /= K;
				
				value[j][i] = cor.toPrecision(pre);			
				value[i][j] = cor.toPrecision(pre);
			}
		}
	}
	
	return value;
}


/// Gets the statistics for a given element in a given the parameter
function get_param_stats(th,index,result,rpf,burn)
{
	let mean, CImin = "-", CImax = "-", ESS = "-", GR = "-";
	let par = result.param[th];

	// Sets the chain
	let chsel = rpf.sel_chain.te;
			
	let name = param_name_index(par,index);

	let vec = []; 
	
	if(par.type == "derive_param" && par.time_dep){
		for(let s = 0; s < result.sample.length; s++){
			let samp = result.sample[s];
			if(samp.num >= burn && !(chsel != "All" && samp.chain != chsel)){
				vec.push(get_param_val(index,samp.param[th],par));
			}
		}
	}
	else{
		for(let s = 0; s < result.par_sample.length; s++){
			let samp = result.par_sample[s];
			if(samp.num >= burn && !(chsel != "All" && samp.chain != chsel)){
				let samp_val;
				if(index == undefined) samp_val = samp.param[th];
				else samp_val = get_param_val(index,samp.param[th],par);
				vec.push(Number(samp_val));
			}
		}
	}

	let stat = get_statistic(vec);

	mean = stat.mean.toPrecision(pre); 
	if(stat.CImin != stat.mean) CImin = stat.CImin.toPrecision(pre);
	if(stat.CImax != stat.mean) CImax = stat.CImax.toPrecision(pre);

	ESS = get_effective_sample_size(vec);
	
	if(result.chains.length > 1 && chsel == "All"){
		let cha = []; 
		for(let c = 0; c < result.chains.length; c++) cha[c] = [];
		
		for(let s = 0; s < result.par_sample.length; s++){
			let samp = result.par_sample[s];
	
			if(samp.num >= burn){
				let ch = find_in(result.chains,samp.chain);
				if(ch == undefined) error("Cannot find chain");
			
				if(index == undefined) cha[ch].push(samp.param[th]);
				else cha[ch].push(get_param_val(index,samp.param[th],par));
			}
		}
		
		GR = get_Gelman_Rubin_statistic(cha);
	}
	
	return { name:name, mean:mean, CImin:CImin, CImax:CImax, ESS:ESS, GR:GR};
}


/// Gets a parameter value
function get_param_val(index,spar,par)
{
	if(par.variety == "const") return get_element(par.value,index);
	return get_element(spar,index);
}
			
			
/// Displays no graph with a message
function no_graph_msg(msg)
{
	return {type:"Graph define", variety:"No Graph", view:"No Graph", data:[], op:{label:msg, key:[]}};
}


/// Sets up a distribution plot
function setup_distribution(result,rpf,burn)
{	
	let val = rpf.sel_paramview.radio.value;
	let th = rpf.sel_paramview.list[val].th;
	let ind = rpf.sel_paramview.list[val].index;
	
	let key = [];
	let data = [];
	
	let par = result.param[th];
		
	let name = par.name;
	
	let line = [];
	switch(rpf.sel_distchain.te){
	case "Combine":
		{
			line.push({ name:posterior_name(result), filt:"All"});
		}
		break;
	
	case "All":
		{
			let chmax = result.chains.length;
			for(let ch = 0; ch < chmax; ch++){
				let lab = "Sample"; if(chmax > 1) lab = "Chain "+result.chains[ch];
				line.push({ name:lab, filt:result.chains[ch]});
			}
		}
		break;
		
	default:
		{
			let ch = rpf.sel_distchain.i;
			line.push({ name:rpf.sel_distchain.te, filt:result.chains[ch]});
		}
		break;
	}
	
	let line_max = false;
	
	let stepp = get_percent_step(line.length);
	for(let li = 0; li < line.length; li++){
		output_percent(stepp,li,line.length);
		
		let col = get_colour(li);
		
		let vec = [];
		let sample = result.par_sample;
		if(par.type == "derive_param" && par.time_dep) sample = result.sample;
		
		for(let i = 0; i < sample.length; i++){
			let samp = sample[i];
			if(samp.num >= burn){
				if(line[li].filt == "All" || line[li].filt == samp.chain){
					vec.push(get_param_val(ind,samp.param[th],par));
				}
			}
		}
	
		let show_mean = rpf.dist_settings.show_mean.check;
		if(line.length > 1) show_mean = false;

		let clip = get_prior_clip(par,ind);
		
		distribution_add_data_line(vec,line[li].name,col,data,key,rpf,show_mean,clip.min,clip.max);
		if(key.length > KEY_LINE_MAX){ line_max = true; break;}
	}

	// Puts on simulation value
	if(par.variety != "const" && rpf.siminf != "sim" && rpf.sim_val.check == true && par.value && par.value.length > 0){
		let val = Number(get_element(par.value,ind));	
		if(!isNaN(val)){
			data.push({te:"Simulation value", col:BLACK, thick:1, x:val, type:"VertDashLine"});
		}
	}
	
	// Draws the prior
	if(result.siminf != "sim" && rpf.dist_settings.show_prior.check == true && par.variety != "reparam" && par.variety != "likelihood" && par.prior && par.type != "derive_param"){
		let N = 200;                                   // The number of points
		
		let pri;
		if(par.dep.length > 0 && par.prior_split_set == true) pri = get_element(par.prior_split,ind);
		else pri = par.prior;
		
		let point = [];
			
		switch(pri.type.te){
		case "uniform":
			{
				let min = Number(pri.value.min_eqn.te);
				let max = Number(pri.value.max_eqn.te);
				point.push({x:min, y:0});
				point.push({x:min, y:1.0/(max-min)});
				point.push({x:max, y:1.0/(max-min)});
				point.push({x:max, y:0});
			}
			break;
			
		case "exp":
			{	
				let mean = Number(pri.value.mean_eqn.te);
				let vec = [];
				for(let i = 0; i < result.par_sample.length; i++){
					let samp = result.par_sample[i]
					if(samp.num >= burn){
						vec.push(samp.param[th]);
					}
				}
		
				let range = extend_range(get_range(vec),0.4);
				if(range.min > 0) range.min = 0;
				if(range.max < 10*mean) range.max = 10*mean;
			
				let dx = range.max-range.max
				for(let i = 0; i <= N; i++){
					let x = range.min + i*(range.max-range.min)/N;
					let prob = Math.exp(-x/mean)/mean;
					point.push({x:x, y:prob});
				}
			}
			break;
			
		case "normal":
			{
				let mean = Number(pri.value.mean_eqn.te);
				let sd = Number(pri.value.sd_eqn.te);
				
				let min = mean - 5*sd;
				let max = mean + 5*sd;
				
				for(let i = 0; i <= N; i++){
					let x = min + i*(max-min)/N;
					let prob = Math.exp(-0.5*Math.log(2*Math.PI*sd*sd) - (x-mean)*(x-mean)/(2*sd*sd));
					point.push({x:x, y:prob});
				}
			}
			break;
			
		case "gamma":
			{
				let mean = Number(pri.value.mean_eqn.te);
				let cv = Number(pri.value.cv_eqn.te);
				let sd = mean*cv;
				let min = 0;
				let max = mean + 5*sd;
				
				let shape = 1.0/(cv*cv);
				let b = shape/mean;
 
				let sum = 0.0;
				let prob = [];
				for(let i = 0; i < N; i++){
					let x = min + (i+0.5)*(max-min)/N;
					
					prob[i] = Math.exp((shape-1)*Math.log(x) - b*x + shape*Math.log(b));
					sum += prob[i];
				}
				let fac = sum*(max-min)/N;
				
				for(let i = 0; i < N; i++) point.push({x:min + (i+0.5)*(max-min)/N, y:prob[i]/fac});
			}
			break;
			
		case "log-normal":
			{
				let mean = Number(pri.value.mean_eqn.te);
				let cv = Number(pri.value.cv_eqn.te);
				let sd = mean*cv;
					
				let min = 0;
				let max = mean + 5*sd;
				
				let vari = Math.log(1+cv*cv);             // Works out variables on log scale
				let mu = Math.log(mean)-vari/2;
		
				for(let i = 0; i < N; i++){
					let x = min + (i+0.5)*(max-min)/N;
					let logx = Math.log(x);
					let prob = Math.exp(-0.5*Math.log(2*Math.PI*vari*x*x) - (logx-mu)*(logx-mu)/(2*vari));
					point.push({x:x, y:prob});
				}
			}
			break;
			
		case "beta":
			{
				let alpha = Number(pri.value.alpha_eqn.te);
				let beta = Number(pri.value.beta_eqn.te);
			
				let min = 0, max = 1;
				
				let sum = 0.0;
				let prob = [];
				for(let i = 0; i < N; i++){
					let x = min + (i+0.5)*(max-min)/N;
					
					prob[i] = Math.exp((alpha-1)*Math.log(x) + (beta-1)*Math.log(1-x));
					sum += prob[i];
				}
				let fac = sum*(max-min)/N;
				
				for(let i = 0; i < N; i++) point.push({x:min + (i+0.5)*(max-min)/N, y:prob[i]/fac});	
			}
			break;
			
		case "bernoulli":
			break;
			
		case "fix":
			break;
		
		case "mdir":
			break;
			
		default: error("prior type not recognised: "+pri.type.te); break; 
		}
		
		if(point.length > 0){
			data.push({point:point, col:BLACK, thick:1, type:"PriorLine"});
		
			key.push({type:"Line", te:"Prior", thick:1, col:BLACK, dash:0});
		}
	}

	post({type:"Graph define", variety:"Distribution", view:"Graph", data:data, op:{x_label:param_name_index(par,ind), x_param:true, y_label:"Probability", param:true, key:key, line_max:line_max, yaxis:false}});
}


/// Gets clipped edges based on a parameter prior
function get_prior_clip(par,ind)
{
	let clip_min, clip_max;
		
	if(par.variety != "reparam" && par.variety != "likelihood" && par.prior){
		let pri;
		if(par.dep.length > 0 && par.prior_split_set == true) pri = get_element(par.prior_split,ind);
		else pri = par.prior;
		
		switch(pri.type.te){
		case "uniform": 
			return {min:Number(pri.value.min_eqn.te), max: Number(pri.value.max_eqn.te)};
			
		case "gamma": case "log-normal":
			return {min:0, max: LARGE};
		
		case "beta":
			return {min:0, max: 1};
		}
	}
	return {min:undefined, max:undefined};
}


/// Adds a distribution line on a distribution plot
function distribution_add_data_line(vec,name,col,data,key,rpf,show_mean,clip_min,clip_max,op)
{
	// Makes sure values are within the range
	if(clip_min != undefined){
		for(let i = 0; i < vec.length; i++){
			if(vec[i] < clip_min) vec[i] = clip_min;
		}
	}
	
	if(clip_max != undefined){
		for(let i = 0; i < vec.length; i++){
			if(vec[i] > clip_max) vec[i] = clip_max;
		}
	}
	
	let stat = get_statistic(vec);
	let range = get_range(vec);

	if(range.min == range.max){
		data.push({te:"Constant", col:col, thick:1, x:range.min, type:"VertLine2"});
		key.push({type:"Line", te:name, col:col});
		return;
	}
	
	let point = [];
			
	switch(rpf.dist_settings.radio.value){
	case "bin":
		{
			let nbin = rpf.dist_settings.sel_bin.te;
			
			let bin = [];
			for(let b = 0; b < nbin; b++) bin[b] = 0;
		
			let mi = range.min;
			let ma = range.max;
			let dx = (ma-mi)/nbin;
			if(op == "integer"){
				dx = Math.round(dx);
				if(dx == 0) dx = 1;
				ma = mi+dx*nbin;
			}
			
			let num = (1.0/vec.length)/dx;
			
			for(let i = 0; i < vec.length; i++){
				bin[Math.floor(ALMOST_ONE*nbin*(vec[i]-range.min)/(range.max-range.min))] += num;
			}
			
			point.push({x:range.min, y:0});
			for(let b = 0; b < nbin; b++){
				point.push({x:range.min+(b/nbin)*(range.max-range.min), y:bin[b]});
				point.push({x:range.min+((b+1)/nbin)*(range.max-range.min), y:bin[b]});
			}
			point.push({x:range.max, y:0});
		}
		break;

	case "kde":
		{
			let h = rpf.dist_settings.sel_h.te;
		
			let d = range.max-range.min;
			let ma = range.max, mi = range.min;
		
			mi -= d*h;
			ma += d*h;
			
			if(op == "integer"){
				mi = Math.floor(mi);
				ma = 1+Math.floor(ma);
			}
			
			if(clip_min != undefined && mi < clip_min) mi = clip_min;
			if(clip_max != undefined && ma > clip_max) ma = clip_max;
		
			let nbin = Math.floor(5/h); 
			if(nbin < 100) nbin = 100; if(nbin > 500) nbin = 500;
			
			if(op == "integer"){ // nbin must be some multiple
				let step = Math.floor((ma-mi)/nbin);
				if(step <= 1) nbin = ma-mi;
				else ma = mi*step;
			}
			
			let dd = ma-mi;	
			
			let db = dd/nbin;
			
			let hb = Math.floor(1+d*h/db);
			
			let bin=[];
			for(let b = 0; b <= nbin; b++) bin[b] = 0;
			
			for(let i = 0; i < vec.length; i++){
				let bf = nbin*(vec[i]-mi)/dd;
				let bmid = Math.floor(bf);
				
				let bmi = bmid-hb, bma = bmid+hb;
			
				for(let b = bmi; b <= bma; b++){
					let bsh = b;//+0.5;
					
					let d;
					if(bsh < bf) d = 1-(bf-bsh)/hb;
					else d = 1+(bf-bsh)/hb;
					if(d > 0){
						let bb = b;
						if(bb < 0) bb = -bb; 
						if(bb > nbin) bb = nbin -(bb-(nbin));
						if(bb >= 0 && bb <= nbin){
							bin[bb] += d;
							if(bb == 0 && clip_min != undefined) bin[bb] += d;
							if(bb == nbin && clip_max != undefined) bin[bb] += d;
						}
					}		
				}
			}
			
			let sum = 0.0; for(let b = 0; b <= nbin; b++) sum += bin[b];
			
			for(let b = 0; b <= nbin; b++) bin[b] /= (db*sum);
			
			point.push({x:mi, y:0});
			for(let b = 0; b <= nbin; b++){
				point.push({x:mi+(b/nbin)*dd, y:bin[b]});
			}
			point.push({x:ma, y:0});
		}
		break;
		
	default: error("Should not be default"); break;
	}
	
	if(show_mean){
		data.push({te:"Mean", col:BLACK, thick:1, x:stat.mean, type:"VertLine"});
	}
	
	data.push({point:point, col:col, CImin:stat.CImin, CImax:stat.CImax, type:"Distribution"});
	
	key.push({type:"Line", te:name, col:col});
}


/// Extends a range
function extend_range(ran,frac)
{
	let d = ran.max-ran.min;
	
	let min_new = ran.min - frac*d; if(min_new*ran.min < 0) min_new = 0;
	let max_new = ran.max + frac*d; if(max_new*ran.max < 0) max_new = 0;
	
	return {min:min_new, max:max_new};
}


/// Works out how often to send a percentage update
function get_percent_step(N)
{
	let step = Math.floor(N/20);
	if(step == 0) step = 1;
	return step;
}

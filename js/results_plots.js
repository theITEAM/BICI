// Functions which set up the plots for results

"use strict";


/// Adds buttons associated with a population plot
function add_pop_buts(result,lay)
{
	//test_graph(lay);
	
	let cx = corner.x;
	let cy = corner.y;
	
	cy = lay.add_title("Populations",cx,cy,{te:pop_plot_text});
	
	cy += 1;

	if(inter.graph.init != true){
		inter.graph.init = true;
		
		// Sets species
		let rpf = result.plot_filter;
		let p = rpf.sel_species.p;

		// Sets classification
		let rpf2 = rpf.species[p];
		let cl = rpf2.sel_class.cl;
	
		// Incorprates filter
		let filt = [];
		
		for(let i = 0; i < rpf2.filter.length; i++){
			let rpf3 = rpf2.filter[i];
			let c = rpf3.sel_comp.c;
			if(!isNaN(c)) filt[rpf3.cl] = c;
		}

		// Sets the chain
		let chsel = rpf.sel_chain.te;
		
		let burn = result.burnin;
		
		// Sets samples
		let imin = 0, imax = result.sample.length;
		let isel = get_isel(rpf); if(isel != "All"){ imin = isel; imax = isel+1; burn = 0;}
	
		let sp = result.species[p];
		let claa = sp.cla[cl];
		
		let view = rpf.sel_view.te;
		
		let data = [];
		
		let key = [];
		
		if(view == "Data"){                                                         // If data is viewed 
			let seldata = rpf.sel_popdata;
			let so = result.species[seldata.p].inf_source[seldata.i];
		
			// Works out which transitions apply to data
			let gfilt = [];	
			for(let c = 0; c < sp.comp_gl.length; c++){
				let cgl = sp.comp_gl[c];
				
				let flag = false;
				
				let cl2;
				for(cl2 = 0; cl2 < sp.ncla; cl2++){
					let filt = so.spec.filter.cla[cl2];
					let cc = cgl.cla_comp[cl2];
					
					if(!(filt.radio.value == "All" || filt.comp[cc].check == true)) break;
				}
				if(cl2 == sp.ncla) flag = true;
				gfilt.push(flag);
			}
			
			{ // Posterior
				let line = [];
			
				for(let i = imin; i < imax; i++){
					let samp = result.sample[i];
					if(samp.num >= burn && !(chsel != "All" && samp.chain != chsel)){
						let sa = samp.species[p];
						
						let point = [];
						
						for(let t = 0; t < samp.timepoint.length; t++){
							let tl = sa.timeline[t];
							
							let pop = 0;
							for(let k = 0; k < sp.comp_gl.length; k++){
								if(gfilt[k] == true) pop += tl.cpop[k];
							}
							point.push({x:samp.timepoint[t], y:pop});
						}
						
						line.push(point);
					}
				}
				add_line(view,line,BLUE,"Posterior",data,key);
			}
			
			{ // Data
				let line = [];
				let tab = so.table;
							
				let spec = so.spec;
				let ty = spec.obs_error.value;
			
				let col = BLACK;
				
				for(let j = 0; j < tab.nrow; j++){
					let t = Number(tab.ele[j][0]);
					let val = Number(tab.ele[j][1]);
					let sd;
					switch(ty){
					case "percent": sd = spec.percent*val/100; break;
					case "sd": sd = spec.sd; break;
					case "file": sd = Number(tab.ele[j][2]); break;
					}
					
					if(!isNaN(t) && !isNaN(val) && !isNaN(sd)){	
						data.push({type:"ErrorBar", x:t, ymin:val-sd, y:val, ymax:val+sd, col:col});
					}
				}
				key.push({type:"ErrBar", te:"Data", col:col});
			}
		}
		else{
			for(let c = 0; c < claa.ncomp; c++){
				let col = claa.comp[c].col;
				
				filt[cl] = c;
				
				let filtg = get_filtg_from_filt(filt,sp);
				
				let line = [];
			
				for(let i = imin; i < imax; i++){
					let samp = result.sample[i];
					if(samp.num >= burn && !(chsel != "All" && samp.chain != chsel)){
						let sa = samp.species[p];
						
						let point = [];
						
						for(let t = 0; t < samp.timepoint.length; t++){
							let tl = sa.timeline[t];
							
							let pop = 0;
							for(let k = 0; k < sp.comp_gl.length; k++){
								if(filtg[k] == true) pop += tl.cpop[k];
							}
							point.push({x:samp.timepoint[t], y:pop});
						}
						
						line.push(point);
					}
				}
				
				switch(view){
				case "Graph": case "Graph (lines)":
					for(let j = 0; j < line.length; j++){
						data.push({point:line[j], col:claa.comp[c].col, view:view, type:"Line"});
					}
					key.push({type:"Line", te:claa.comp[c].name, col:claa.comp[c].col});
					break;
					
				case "Graph (CI)": case "Data":
					{
						let line_stats = get_line_stats(line);
						data.push({point:line_stats, col:col, type:"Line CI"});
						key.push({type:"Line", te:claa.comp[c].name, col:col});
					}
					break;
					
				case "Compartment":
					let line_stats = get_line_stats(line);
					data.push({point:line_stats, col:col, type:"Comp Video"});
					break;
					
				default: error("view option error1"); break;
				}
			}
		}
		
		switch(view){
		case "Graph": case "Graph (lines)": case "Graph (CI)": case "Data":
			inter.graph.define("Population",view,data,{x_label:"Time", y_label:"Population", key:key});
			break;
			
		case "Compartment": 
			inter.graph.define("Population",view,data,{species:result.species, p:p, cl:cl});
			break;
			
		default: error("Type not recognised"); break;
		}
	}
	
	switch(plot_variety(inter.graph.type)){
	case "Line plot": inter.graph.create(cx,cy,graph_width-right_menu_width,29,lay); break;
	case "Comp plot": inter.graph.create(0,3,graph_width+2-right_menu_width,lay.dy-3,lay); break;
	}
}


/// Determines if two vectors are the same
function vec_equal(vec1,vec2)
{
	if(vec1.length != vec2.length) return false;
	for(let i = 0; i < vec1.length; i++) if(vec1[i] != vec2[i]) return false;
	return true;
}


/// Convert from filt to filtg
function get_filtg_from_filt(filt,sp)
{
	let filtg = [];
			
	for(let k = 0; k < sp.comp_gl.length; k++){
		let co = sp.comp_gl[k];
		let cl2 = 0; 
		while(cl2 < sp.ncla && (filt[cl2] == undefined || co.cla_comp[cl2] == filt[cl2])) cl2++;
					
		if(cl2 == sp.ncla) filtg.push(true);
		else filtg.push(false);
	}
	
	return filtg;
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
function add_trans_buts(result,lay)
{
	let cx = corner.x;
	let cy = corner.y;

	cy = lay.add_title("Transition rates",cx,cy,{te:trans_plot_text});
	
	cy += 1;

	if(inter.graph.init != true){
		inter.graph.init = true;
		
		// Sets species
		let rpf = result.plot_filter;
		let p = rpf.sel_species.p;

		// Sets classification
		let rpf2 = rpf.species[p];
		let cl = rpf2.sel_class.cl;
	
		// Incorprates filter
		let filt = [];
		for(let i = 0; i < rpf2.filter.length; i++){
			let rpf3 = rpf2.filter[i];
			let c = rpf3.sel_comp.c;
			if(!isNaN(c)) filt[rpf3.cl] = c;
		}

		// Sets the chain
		let chsel = rpf.sel_chain.te;
		
		let burn = result.burnin;
		
		// Sets samples
		let imin = 0, imax = result.sample.length;
		let isel = get_isel(rpf); if(isel != "All"){ imin = isel; imax = isel+1; burn = 0;}

		let sp = result.species[p];
		let claa = sp.cla[cl];
	
		let step = rpf.sel_timestep.i;
			
		let data = [];
	
		let key = [];
	
		let view = rpf.sel_trans_view.te;
		
		if(view == "Data"){                                                         // If data is viewed 
			let seldata = rpf.sel_poptransdata;
			let so = result.species[seldata.p].inf_source[seldata.i];
		
			let cl = find(sp.cla,"name",so.spec.cl_drop.te);
			
			// Works out which transitions apply to data
			let gfilt = [];	
			for(let tr = 0; tr < sp.tra_gl.length; tr++){
				let trg = sp.tra_gl[tr];
				
				let flag = false;
				if(trg.cl == cl){
					let c = trg.i; if(c == "Source") c = trg.f;
					
					let cgl = sp.comp_gl[c];
					
					let cl2;
					for(cl2 = 0; cl2 < sp.ncla; cl2++){
						if(cl == cl2){
							if(so.spec.filter.tra[trg.tr].check != true) break;
						}
						else{
							let filt = so.spec.filter.cla[cl2];
							let cc = cgl.cla_comp[cl2];
							
							if(!(filt.radio.value == "All" || filt.comp[cc].check == true)) break;
						}
					}
					if(cl2 == sp.ncla) flag = true;
				}
				gfilt.push(flag);
			}
			
			{ // Posterior
				let line = [];
				
				for(let i = imin; i < imax; i++){
					let samp = result.sample[i];
					if(samp.num >= burn && !(chsel != "All" && samp.chain != chsel)){
						let sa = samp.species[p];
						
						let point = [];
						
						let tmax = samp.timepoint.length;
						for(let t = 0; t < tmax-1; t += step){
							let step_max = t + step; if(step_max > tmax-1) step_max = tmax-1;
							
							let dt = samp.timepoint[step_max] - samp.timepoint[t];
							
							let num = 0;
							for(let k = 0; k < sp.tra_gl.length; k++){
								if(gfilt[k] == true){
									for(let tt = t; tt < step_max; tt++){
										num += sa.timeline[tt].transnum[k];
									}
								}
							}
							point.push({x:samp.timepoint[t] + dt/2, y:num/dt});
						}
						line.push(point);
					}
				}
				add_line(view,line,BLUE,"Posterior",data,key);
			}
			
			{ // Data
				let line = [];
				let tab = so.table;
							
				let spec = so.spec;
				let ty = spec.obs_error.value;
			
				let dtdata = Number(tab.ele[1][0]) - Number(tab.ele[0][0])
				
				let col = BLACK;
				
				let point = [];		
				for(let j = 0; j < tab.nrow; j++){
					let t = Number(tab.ele[j][0]);
					let val = Number(tab.ele[j][1]);
					let sd;
					switch(ty){
					case "percent": sd = spec.percent*val/100; break;
					case "sd": sd = spec.sd; break;
					case "file": sd = Number(tab.ele[j][2]); break;
					}
					
					if(!isNaN(t) && !isNaN(val) && !isNaN(sd)){
						point.push({x:t, y:val/dtdata});
						point.push({x:t+dtdata, y:val/dtdata});
						
						data.push({type:"ErrorBar", x:t+dtdata/2, ymin:(val-sd)/dtdata, y:val/dtdata, ymax:(val+sd)/dtdata, col:col});
					}
				}
				line.push(point);
				add_line(view,line,BLACK,"Data",data,key);
			}
		}
		else{                                                                  // General transitions are observed
			// Gets colours for the lines
			let col_line = [];
			for(let j = 0; j < claa.ntra; j++){
				let tr = claa.tra[j];
				let col_pos;
				if(tr.i == "Source" || tr.f == "Sink") col_pos = BLACK;
				else col_pos = claa.comp[tr.i].col;
				
				if(find_in(col_line,col_pos) == undefined) col_line.push(col_pos);
				else{
					let k = 0; while(k < auto_color.length && find_in(col_line, auto_color[k]) != undefined) k++;
					if(k < auto_color.length) col_line.push(auto_color[k]);
					else col_line.push(auto_color[Math.floor(Math.random()*auto_color.length)]);
				}
			}
			
			for(let j = 0; j < claa.ntra; j++){
				let tr = claa.tra[j];
				
				let col = col_line[j];
			
				let line = [];
			
				for(let i = imin; i < imax; i++){
					let samp = result.sample[i];
					if(samp.num >= burn && !(chsel != "All" && samp.chain != chsel)){
						let sa = samp.species[p];
						
						let point = [];
						
						let tmax = samp.timepoint.length;
						for(let t = 0; t < tmax-1; t += step){
							let step_max = t + step; if(step_max > tmax-1) step_max = tmax-1;
							
							let dt = samp.timepoint[step_max] - samp.timepoint[t];
							
							let num = 0;
							for(let k = 0; k < sp.tra_gl.length; k++){
								let tr = sp.tra_gl[k];
								if(tr.cl == cl && tr.tr == j){
									let cc = tr.i; if(cc == "Source") cc = tr.f;
									
									let co = sp.comp_gl[cc];
									let cl2 = 0; 
									while(cl2 < sp.ncla && (filt[cl2] == undefined || co.cla_comp[cl2] == filt[cl2])) cl2++;
								
									if(cl2 == sp.ncla){
										for(let tt = t; tt < step_max; tt++){
											num += sa.timeline[tt].transnum[k];
										}
									}
								}
							}
							point.push({x:samp.timepoint[t] + dt/2, y:num/dt});
						}
						line.push(point);
					}
				}
				add_line(view,line,col,tr.name,data,key);
			}
		}
		
		inter.graph.define("Transition",view,data,{x_label:"Time", y_label:"Transition rate", key:key});
	}
	
	inter.graph.create(cx,cy,graph_width-right_menu_width,29,lay); 
}


/// Adds a line to the graph
function add_line(view,line,col,name,data,key)
{
	switch(view){
	case "Graph": case "Graph (lines)":
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
function add_individual_buts(result,lay)
{
	let cx = corner.x;
	let cy = corner.y;
	
	cy = lay.add_title("Individuals",cx,cy,{te:ind_plot_text});
	
	// Sets species
	let rpf = result.plot_filter;	
	let p = rpf.sel_species.p;
	
	if(result.species[p].type == "Population"){
		let dx = 20;
		lay.add_paragraph("No individual data for a population model.",dx,lay.dx/2-dx/2,lay.dy/2,BLACK,para_si,para_lh);	
		return;
	}

	cy += 1;

	if(inter.graph.init != true){
		inter.graph.init = true;
		
		// Sets classification
		let rpf2 = rpf.species[p];
		let cl = rpf2.sel_class.cl;
	
		// Sets the chain
		let chsel = rpf.sel_chain.te;
		
		let burn = result.burnin;
		
		// Sets samples
		let imin = 0, imax = result.sample.length;
		let isel = get_isel(rpf); if(isel != "All"){ imin = isel; imax = isel+1; burn = 0;}
		
		let sp = result.species[p];
		let claa = sp.cla[cl];
	
		let tl = result.sample[0].species[p].timeline;
		
		let tp = result.sample[0].timepoint;
		let tmin = tp[0], tmax = tp[tp.length-1];

		let key = [];
		let data = [];
		
		/// First makes a list of individuals from all samples
		let ind_list = [];
		for(let samp = imin; samp < imax; samp++){
			let sampl = result.sample[samp];
			if(sampl.num >= burn && !(chsel != "All" && sampl.chain != chsel)){		
				let sa = sampl.species[p];
				for(let i = 0; i < sa.individual.length; i++){
					let ind = sa.individual[i];
					
					let name = ind.name;
					let j = find(ind_list,"name",name);
					if(j == undefined){
						j = ind_list.length;
						ind_list.push({name:name, sa_ref:[], obs:[]});
					}
					ind_list[j].sa_ref.push({sample:samp, i:i});
				}
			}
		}

		if(result.siminf == "inf"){
			// Puts observed data on the time line
			for(let i = 0; i < sp.inf_source.length; i++){
				let so = sp.inf_source[i];
		
				switch(so.type){
				case "Transition":
					if(so.table_loaded == true){
						let tab = so.table;

						let spec = so.spec;
						let name = spec.tr_drop.te;
						let spl = name.split("→");
						
						let col = []; 
						let fl = false;
						for(let loop = 0; loop < 2; loop++){
							let na = spl[loop];
							col[loop] = WHITE;
							
							for(let c = 0; c < claa.ncomp; c++){
								if(claa.comp[c].name == na){ col[loop] = claa.comp[c].col; fl = true;}
							}
						}
					
						if(fl == true){
							for(let r = 0; r < tab.nrow; r++){
								let j = find(ind_list,"name",tab.ele[r][0]);
								let t = Number(tab.ele[r][1]);
								if(j != undefined){
									ind_list[j].obs.push({type:"Transition", t:t, col:col});
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
								let j = find(ind_list,"name",tab.ele[r][0]);
								let t = Number(tab.ele[r][1]);
								let cname = tab.ele[r][2];
								if(j != undefined){
									let k = find(claa.comp,"name",cname);			
									ind_list[j].obs.push({type:"Compartment", t:t, col:claa.comp[k].col});
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
								let j = find(ind_list,"name",tab.ele[r][0]);
								let t = Number(tab.ele[r][1]);
								
								let res_str = tab.ele[r][2];
								let res;
								if(res_str == spec.neg_result) res = false;
								if(res_str == spec.pos_result) res = true;
								
								if(res != undefined && j != undefined){
									ind_list[j].obs.push({type:"Diag. Test", t:t, res:res});
								}
							}
						}
					}
					break;
				}
			}
		}
		
		// Gets statistics about individual effects
		for(let i = 0; i < ind_list.length; i++){
			let comb_ind = ind_list[i];
			
			let ie_vec = [];
			let nie = sp.ind_effect.length;
			for(let e = 0; e < nie; e++) ie_vec[e] = [];
			
			let kmax = comb_ind.sa_ref.length;
			for(let k = 0; k < kmax; k++){
				let sr = comb_ind.sa_ref[k];
				let ind = result.sample[sr.sample].species[p].individual[sr.i];
				if(ind.ie.length != nie) error("ie not the right size");
				for(let e = 0; e < nie; e++) ie_vec[e].push(ind.ie[e]);
			}
			
			comb_ind.ie_stat = [];
			for(let e = 0; e < nie; e++){
				comb_ind.ie_stat[e] = get_statistic(ie_vec[e]);
			}
		}
		
		// Goes through all individuals from all samples
		
		switch(rpf.sel_indview.te){
		case "Timeline": add_timelines(ind_list,tmin,tmax,p,cl,result,data,key); break
		case "Table": add_individual_table(ind_list,p,result,data); break;
		default: error("Option not here"); break;
		}
	}
	
	switch(plot_variety(inter.graph.type)){
	case "Stat table plot": inter.graph.create(cx,cy,graph_width-right_menu_width,29,lay); break;
	case "Individual plot": inter.graph.create(0,cy,graph_width+2-right_menu_width,29,lay); break;
	}
}


/// Adds individual table
function add_individual_table(ind_list,p,result,data)
{
	let sp = result.species[p];
	
	let twid = inter.graph.table_width;
	twid = [];
	
	let table = {width:twid, heading:[], content:[]};
	
	for(let i = 0; i < ind_list.length; i++) table.content.push([]);
		
	let vec = [];
	for(let i = 0; i < ind_list.length; i++) vec.push(ind_list[i].name);
	
	add_table_column("ID",vec,table);
	
	let nie = sp.ind_effect.length;
	for(let e = 0; e < nie; e++){	
		let vec=[];
		for(let i = 0; i < ind_list.length; i++){
			let stat = ind_list[i].ie_stat[e];
			let te = stat.mean.toPrecision(pre);
			if(stat.CImin != stat.mean || stat.CImax != stat.mean){
				te += " ("+stat.CImin.toPrecision(pre)+" — "+stat.CImax.toPrecision(pre)+")";
			}
			vec.push(te);
		}
		
		add_table_column(sp.ind_effect[e],vec,table);
	}	
	
	let nfe = sp.fix_eff.length; 
	for(let f = 0; f < nfe; f++){
		let fe = sp.fix_eff[f];
		
		let X = fe.X_vector
		let Xind_list = X.ind_list;
		let value = X.value;
		
		let vec=[];
		for(let i = 0; i < ind_list.length; i++){
			let val;
			if(i < Xind_list.length && ind_list[i].name ==  Xind_list[i]) val = value[i];
			else{
				let j = find(Xind_list,ind_list[i].name);
				if(j != undefined) val = value[j];
				else val = "-";
			}
			vec.push(val);
		}
		
		add_table_column("X for "+fe.name,vec,table);
	}
	
	/// Adds padding on the last column (if possible)
	let wid = 0;
	for(let i = 0; i < twid.length; i++) wid += twid[i];
	
	let wid_max = graph_width-right_menu_width-2;
	if(wid < wid_max) twid[twid.length-1] += wid_max-wid;
	
	data.push({type:"Table", table:table});
	
	inter.graph.define("Stat Table","Stat Table",data,{});
}


/// Adds a column to a table
function add_table_column(head,vec,table)
{
	let fo = get_font(si_table);
	
	let wmax = 0;
	table.heading.push({name:head});
	let w = text_width(head,fo); if(w > wmax) wmax = w;
	
	for(let j = 0; j < vec.length; j++){
		table.content[j].push({te:vec[j]});
		let w = text_width(vec[j],fo); if(w > wmax) wmax = w;
	}
	wmax += 2;
	
	if(wmax < 5) wmax = 5;
	
	table.width.push(wmax);
}


/// Adds individual timelines
function add_timelines(ind_list,tmin,tmax,p,cl,result,data,key)
{
	let sp = result.species[p];
	let claa = sp.cla[cl];
	
	/// Works out colours for each of the compartments
	let comp_col = [];
	for(let c = 0; c < claa.comp.length; c++){
		let co = claa.comp[c];
		let col = co.col;
		let bigint = parseInt(col.substring(1), 16);	
		comp_col.push({r:(bigint >> 16) & 255, g:(bigint >> 8) & 255, b:bigint & 255});
		
		key.push({type:"Rect", te:co.name, col:co.col});
	}
		
	for(let i = 0; i < ind_list.length; i++){
		if(find(ind_list[i].obs,"type","Compartment") != undefined){
			key.push({type:"ObsComp", te:"Data"});
			break;
		}
	}
	
	for(let i = 0; i < ind_list.length; i++){
		let comb_ind = ind_list[i];
		
		// Gets information about transitions			
		let cpop = [];
		for(let c = 0; c < claa.ncomp; c++) cpop[c] = 0; 
		let nalive = 0;
	
		let nsamp = comb_ind.sa_ref.length;
		
		let cla_ev = [];
		for(let k = 0; k < nsamp; k++){
			let sr = comb_ind.sa_ref[k];
			let ind = result.sample[sr.sample].species[p].individual[sr.i];
			
			let c = "out";
			
			if(ind.cinit != "out" && ind.cinit != "Source"){
				nalive++; 
				let cgl = sp.comp_gl[ind.cinit];
				c = cgl.cla_comp[cl];
				cpop[c]++;
			}
			
			for(let e = 0; e < ind.ev.length; e++){
				let ev = ind.ev[e];
			
				switch(ev.type){
				case "trans":
					let tr = sp.tra_gl[ev.trg];
					if(tr.cl == cl || tr.i == "Source"){
						let cf;
						if(tr.f == "Sink") cf = "out"; 
						else cf = sp.comp_gl[tr.f].cla_comp[cl];	
						cla_ev.push({t:ev.t, ci:c, cf:cf});
						c = cf;
					}
					break;
					
				case "enter":
					let cf = sp.comp_gl[ev.c].cla_comp[cl];	
					cla_ev.push({t:ev.t, ci:c, cf:cf});
					c = cf;
					break;
					
				case "leave":
					cla_ev.push({t:ev.t, ci:c, cf:"out"});
					c = "out";
					break;
					
				case "move":
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
		
		let col_timeline = [];
			
		col_timeline.push({t:tmin, col:get_mix_colour(comp_col,cpop), alpha:nalive/nsamp });
				
		for(let k = 0; k < cla_ev.length; k++){
			let clev = cla_ev[k];	
			let ci = clev.ci; if(ci == "out") nalive++; else cpop[ci]--;
			let cf = clev.cf; if(cf == "out") nalive--; else cpop[cf]++;
			
			col_timeline.push({t:clev.t, col:get_mix_colour(comp_col,cpop), alpha:nalive/nsamp});
		}
		
		col_timeline.push({t:tmax, col:get_mix_colour(comp_col,cpop), alpha:nalive/nsamp});
		
		let info = "";
		
		for(let e = 0; e < comb_ind.ie_stat.length; e++){
			if(e != 0) info += "   ";
			if(nsamp == 1) info += sp.ind_effect[e];
			else info += "〈"+sp.ind_effect[e]+"〉";
			info += "= "+comb_ind.ie_stat[e].mean.toPrecision(pre);
		}
	
		data.push({type:"Individual", name:comb_ind.name, info:info, tmin:tmin, tmax:tmax, col_timeline:col_timeline, obs:comb_ind.obs});	
	}

	inter.graph.define("Individual","Individual",data,{x_label:"Time", key:key});
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


/// Adds buttons associated with a spline plot
function add_spline_buts(result,lay)
{
	let cx = corner.x;
	let cy = corner.y;
	
	cy = lay.add_title("Spline variation",cx,cy,{te:spline_plot_text});
	
	cy += 1;

	if(inter.graph.init != true){
		inter.graph.init = true;
		
		// Sets spline parameter
		let rpf = result.plot_filter;
		let spl = rpf.sel_spline.i; 

		let burn = result.burnin;

		// Sets samples
		let imin = 0, imax = result.sample.length;
		let isel = rpf.sel_sample.i;
		
		if(isel != "All"){ imin = isel; imax = isel+1; burn = 0;}

		let view = rpf.sel_spline_view.te;
		
		let data = [];
		
		let key = [];
		
		let spline = result.spline[spl];
		
		let par = result.param[spline.th];
	
		let col_list = get_col_list(par.dep,result,spline.comb_list.length);
		
		let line_ref = [];
		
		for(let i = 0; i < spline.comb_list.length; i++){  // There is potential for filtering
			
			line_ref.push({index:spline.comb_list[i].index, col:col_list[i], name:spline.comb_list[i].name});
		}
		
		for(let j = 0; j < line_ref.length; j++){
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
		}
		
		inter.graph.define("Population",view,data,{x_label:"Time", y_label:spline.name, y_param:true, key:key});
	}
	
	inter.graph.create(cx,cy,graph_width-right_menu_width,29,lay); 
}

	
/// Based on a given dependency calculate a colour list for lines
function get_col_list(dep,mod,N)
{
	let col_list = [];    
	if(dep.length == 2){    // If only one index then uses the colour scheme from the compartments
		let indname = dep[0];
		
		let cla = get_cla_from_index(mod,dep[0])
		if(cla == undefined) error("Could not find classification");
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
		return sp.cla[cl];
	}
	return;
}

				

/// Adds a line to the data
function add_line_to_data(view,name,col,line,key,data)
{
	switch(view){
	case "Graph": case "Graph (lines)": 
		for(let k = 0; k < line.length; k++){
			data.push({point:line[k], col:col, view:view, type:"Line"});
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
	
	lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, type:"Rect", val:WHITE});
	
	let warn = "";
	
	let cx = corner.x;
	let cy = corner.y;
	
	cy = lay.add_title(vg.title,cx,cy,{});
	
	cy += 1;

	if(gr.init != true){
		gr.init = true;
		
		let data = [];
		
		let key = [];
		
		switch(vg.type){
		case "ViewSpline":
			{
				let par = so.param[vg.i];
				
				let value = par.value;
				
				let knot = par.spline.knot;
				let times = get_times(knot,so.sim_details.t_start,so.sim_details.t_end);
				
				if(typeof times == 'string'){ alertp(times); close_view_graph(); return;}
				
				let ntimes = times.length;
				
				let list = par_find_list(par);
				if(list[list.length-1].length != ntimes) error("Should be same length");
				list.pop();
				
				let comb_list = generate_comb_list(list);
				
				let col_list = get_col_list(par.dep,so,comb_list.length);

				for(let i = 0; i < comb_list.length; i++){
					let line = [];
					
					let vec = get_element(value,comb_list[i].index);
					
					let point = [];
					for(let t = 0; t < ntimes; t++) point.push({x:times[t], y:vec[t]});
				
					let name = get_spline_name(par,comb_list[i].index,list);
					
					line.push(point);
					add_line_to_data("Graph",name,col_list[i],line,key,data);
				}
		
				gr.define("Population","Graph",data,{x_label:"Time", y_label:par.name+"(t)", y_param:true, key:key});
			}
			break;
			
		case "ViewVector":
			{
				let par = so.param[vg.i];
				
				let value = par.value;
				
				let list = par_find_list(par);
				let vec = list[0];
			
				let col_list = get_col_list(par.dep,so,vec.length);

				for(let i = 0; i < vec.length; i++){
					let ind = [i]; 
					data.push({type:"Bar", name:param_name_index(par,ind), x:i+0.5, y:value[i], thick:0.7, col:col_list[i]});
				}
				
				let claa = get_cla_from_index(so,par.dep[0]);
				inter.graph.define("Histogram","Histogram",data,{x_label:claa.name, y_label:"Value"});
			}
			break;
			
		case "ViewMatrix":
			{
				let par = so.param[vg.i];

				let value = par.value;
				
				let list = par_find_list(par);
				let vecy = list[0];
				let vecx = list[1];

				data.push({type:"Matrix", xlab:vecx, ylab:vecy, mat:value});
		
				inter.graph.define("Matrix","Matrix",data,{x_label:par.dep[1], x_param:true, y_label:par.dep[0], y_param:true});
			}
			break;
			
		default: error("Graph view option not recognised"); break;
		}
	}
	
	inter.graph.create(cx,cy,graph_width-right_menu_width,29,lay); 
	
	lay.add_corner_button([["Back","Grey","CancelViewGraph"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
}


/// Closes the data source 
function close_view_graph()
{
	inter.view_graph = {};
	inter.graph.init = undefined;
	inter.graph.type = undefined;
}


/// Adds buttons for the parameters page
function add_parameter_buts(result,lay)
{
	let cx = corner.x;
	let cy = corner.y;
	
	cy = lay.add_title("Parameters",cx,cy,{te:param_plot_text});
	
	if(inter.graph.init != true){
		inter.graph.init = true;
		
		// Sets species
		let rpf = result.plot_filter;
	
		let pview = rpf.sel_paramview;
		
		if(pview.list.length > 0){
			let li = pview.list[0];
			let par = result.param[li.th];
			rpf.pos_paramviewtype = [];
			
			if(result.siminf == "inf"){
				if(alg(result) == "DA-MCMC") rpf.pos_paramviewtype.push({te:"Trace"});
				else rpf.pos_paramviewtype.push({te:"Samples"});
					
				rpf.pos_paramviewtype.push({te:"Distribution"});	
			}
			rpf.pos_paramviewtype.push({te:"Statistics"});
		
			if(par.dep.length == 1){
				rpf.pos_paramviewtype.push({te:"Histogram"});
			}
			if(par.dep.length == 2){
				rpf.pos_paramviewtype.push({te:"Matrix"});
			}
		
			if(!rpf.sel_paramviewtype || find(rpf.pos_paramviewtype,"te",rpf.sel_paramviewtype.te) == undefined){
				rpf.sel_paramviewtype = copy(rpf.pos_paramviewtype[0]);
			}
		}
		
		let data = [];
		
		switch(rpf.sel_paramview.te){
		case "Scatter":
			{
				let sp = rpf.sel_paramview;
				let valx = sp.radio.value;
				let thx = sp.list[valx].th;
				let indx = sp.list[valx].index;
				let namex = param_name_index(result.param[thx],indx);
						
				let valy = sp.radioy.value;
				let thy = sp.list[valy].th;
				let indy = sp.list[valy].index;
				let namey = param_name_index(result.param[thy],indy);
				
				let chmin = 0, chmax = result.chains.length;
				if(rpf.sel_chain.te != "All"){
					chmin = find_in(result.chains,rpf.sel_chain.te);
					chmax = chmin+1;
				}
				
				let point=[];	
				for(let ch = chmin; ch < chmax; ch++){
					let cha = result.chains[ch];
					
					for(let i = 0; i < result.par_sample.length; i++){
						let samp = result.par_sample[i];
						if(samp.num >= result.burnin && samp.chain == cha){
							let valx = get_element(samp.param[thx],indx);
							let valy = get_element(samp.param[thy],indy);
							point.push({x:valx, y:valy});
						}
					}
				}
				data.push({point:point, type:"Points", col:BLUE});
			
				inter.graph.define("Scatter","Scatter",data,{x_label:namex, x_param:true, y_label:namey, y_param:true, key:[]});
			}
			break;
			
		case "Correlation":
			{	
				let mat = get_correlation_matrix(result);
			
				let vec = [];
				for(let j = 0; j < result.total_param_list.length; j++){
					vec.push(result.total_param_list[j].name);
				}
				
				data.push({type:"Correlation", xlab:vec, ylab:vec, x_param:true, y_param:true, mat:mat});

				inter.graph.define("Matrix","Matrix",data,{x_label:"Parameter", y_label:"Parameter"});
			}
			break;
			
		default:
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
						let imin = 0;
						while(imin < result.par_sample.length && result.par_sample[imin].chain != cha) imin++;
						
						let iburn = imin;
						
						if(vte == "Trace"){
							while(iburn < result.par_sample.length && result.par_sample[iburn].num < result.burnin) iburn++;
							
							if(iburn > imin){
								let point=[];
								for(let i = imin; i <= iburn; i++){
									let samp = result.par_sample[i];
									if(samp.chain == cha){
										let val;
										if(par.kind == "const") val = get_element(par.value,ind)
										else val = get_element(samp.param[th],ind);
										point.push({x:samp.num, y:val});
									}
								}
						
								let col = light_colour(get_colour(ch-chmin));
						
								data.push({point:point, col:col, type:"BurninLine"});
							}
						}
						
						let point=[];
						for(let i = iburn; i < result.par_sample.length; i++){
							let samp = result.par_sample[i];
							if(samp.chain == cha){
								let val;
								if(par.kind == "const") val = get_element(par.value,ind)
								else val = get_element(samp.param[th],ind);
								point.push({x:samp.num, y:val});
							}
						}
						
						let col = get_colour(ch-chmin);
						
						data.push({point:point, col:col, type:"Line"});
					
						key.push({type:"Line", te:"Chain "+cha, col:col});
					}
					
					if(vte == "Trace"){
						data.push({te:"Burn-in", col:BLACK, thick:1, x:result.burnin, type:"VertDashLine"});
					}
					
					if(key.length == 1) key = [];
					inter.graph.define(vte,"Graph",data,{x_label:"Sample", y_label:name, y_param:true, key:key});
				}
				break;
			
			case "Distribution":
				{
					setup_distribution(result);
				}
				break;
				
			case "Statistics":
				{	
					inter.graph.table_width = [8,8,7,7,5,5];
					let table = {width:inter.graph.table_width, heading:[{name:"Parameter"},{name:"Mean"},{name:"CI min"},{name:"CI max"},{name:"ESS"},{name:"GR"}], content:[]};
		
					for(let i = 0; i < pview.list.length; i++){
						let li = pview.list[i];
						let stat = get_param_stats(li.th,li.index,result);
				
						table.content.push([{te:stat.name, pname:true},{te:stat.mean},{te:stat.CImin},{te:stat.CImax},{te:stat.ESS},{te:stat.GR}]);
					}
		
					data.push({type:"Table", table:table});
			
					inter.graph.define("Stat Table","Stat Table",data,{});
				}
				break;
				
			case "Histogram":
				{
					let li = pview.list[0];
				
					let par = result.param[li.th];
					
					let col_list = get_col_list(par.dep,result,pview.list.length);

					for(let i = 0; i < pview.list.length; i++){
						let li = pview.list[i];

						let stat = get_param_stats(li.th,li.index,result);

						data.push({type:"Bar", name:stat.name, x:i+0.5, y:stat.mean, thick:0.7, col:col_list[i]});
						data.push({type:"ErrorBar", x:i+0.5, ymin:stat.CImin, y:stat.mean, ymax:stat.CImax, col:BLACK});
					}
				
					let de = par.dep[0];
					let name;
					if(de == "t") name = "Time";
					else name = get_cla_from_index(result,de).name;
					inter.graph.define("Histogram","Histogram",data,{x_label:name, y_label:"Value"});
				}
				break;
				
			case "Matrix":
				{
					let li = pview.list[0];
					let par = result.param[li.th];

					let value = par_find_template(par.list);

					for(let i = 0; i < pview.list.length; i++){
						let li = pview.list[i];

						let stat = get_param_stats(li.th,li.index,result);
						set_element(value,li.index,stat.value);
					}
					
					let vecy = par.list[0];
					let vecx = par.list[1];

					data.push({type:"Matrix", xlab:vecx, ylab:vecy, mat:value});
			
					inter.graph.define("Matrix","Matrix",data,{x_label:par.dep[1], x_param:true, y_label:par.dep[0], y_param:true});
				}
			}
			break;
		}
	}
	
	
	let x = cx, y = cy;
	let w = graph_width-right_menu_width;
	let h = 29;
	
	if(inter.graph.variety == "Scatter"){   // Makes square
		let tw = inter.graph.tick.wmax;
		let ww = w-graph_mar.right-graph_mar.left-tw;
		let hh = h-graph_mar.top-graph_mar.bottom;
		if(ww > hh){ 
			let d = (ww-hh)/2;
			x += 0.6*d; w -= 2*d;
		}
	}
				
	inter.graph.create(x,y,w,h,lay);
}


/// Gets a matrix giving the correlations between parameters
function get_correlation_matrix(result)
{
	let rpf = result.plot_filter;
	
	let list = result.total_param_list;
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
			if(samp.num >= result.burnin && samp.chain == cha){
				for(let j = 0; j < N; j++){
					trace[j].push(get_element(samp.param[list[j].th],list[j].index));
				}
			}
		}
	}

	// Transforms such that trace has mean zero and variance 1;
	let K = trace[0].length;
	
	for(let j = 0; j < N; j++){
		let sum = 0.0, sum2 = 0.0;
		for(let k = 0; k < K; k++){
			sum += trace[j][k]; sum2 += trace[j][k]*trace[j][k];
		}
	
		let mean = sum/K;
		let vari = sum2/K - mean*mean; if(vari < 0) vari = 0;
		let sd = Math.sqrt(vari);
		for(let k = 0; k < K; k++) trace[j][k] = (trace[j][k]-mean)/sd;
	}
	
	for(let j = 0; j < N; j++){
		for(let i = j; i < N; i++){
			let cor = 0.0;
			for(let k = 0; k < K; k++) cor += trace[j][k]*trace[i][k];
			cor /= K;
			
			value[j][i] = cor.toPrecision(pre);			
			value[i][j] = cor.toPrecision(pre);
		}
	}
	
	return value;
}
				
				
/// Gets the statistics for a given element in a given the parameter
function get_param_stats(th,index,result)
{
	let mean, CImin = "-", CImax = "-", ESS = "-", GR = "-";
	let par = result.param[th];
	
	// Sets the chain
	let chsel = result.plot_filter.sel_chain.te;
			
	let name = param_name_index(par,index);

	if(par.kind == "const"){	
		if(index == undefined) mean = par.value;
		else mean = get_element(par.value,index);
	}
	else{
		let vec = []; 
		for(let s = 0; s < result.par_sample.length; s++){
			let samp = result.par_sample[s];
		
			if(samp.num >= result.burnin && !(chsel != "All" && samp.chain != chsel)){
				if(index == undefined) vec.push(samp.param[th]);
				else vec.push(get_element(samp.param[th],index));
			}
		}
	
		let stat = get_statistic(vec);
	
		mean = stat.mean.toPrecision(pre); 
		CImin = stat.CImin.toPrecision(pre);
		CImax = stat.CImax.toPrecision(pre);
		
		ESS = get_effective_sample_size(vec);
		
		if(result.chains.length > 1 && chsel == "All"){
			let cha = []; 
			for(let c = 0; c < result.chains.length; c++) cha[c] = [];
			
			for(let s = 0; s < result.par_sample.length; s++){
				let samp = result.par_sample[s];
		
				if(samp.num >= result.burnin){
					let ch = find_in(result.chains,samp.chain);
					if(ch == undefined) error("Cannot find chain");
				
					if(index == undefined) cha[ch].push(samp.param[th]);
					else cha[ch].push(get_element(samp.param[th],index));
				}
			}
			
			GR = get_Gelman_Rubin_statistic(cha).toPrecision(pre);
		}
	}
	
	return { name:name, mean:mean, CImin:CImin, CImax:CImax, ESS:ESS, GR:GR};
}


/// Sets up a distribution plot
function setup_distribution(result)
{
	let rpf = result.plot_filter;
		
	let val = rpf.sel_paramview.radio.value;
	let th = rpf.sel_paramview.list[val].th;
	let ind = rpf.sel_paramview.list[val].index;
	
	let key = [];
	
	let data = [];
	
	let par = result.param[th];
	if(par.kind == "const"){
		inter.graph.define("No Graph","No Graph",data,{label:"This parameter is constant", key:key});
		return;
	}
		
	let name = result.param[th].name;
	
	
	
	let line = [];
	switch(rpf.sel_distchain.te){
	case "Combine":
		{
			line.push({ name:"Posterior", filt:"All"});
		}
		break;
	
	case "All":
		{
			for(let ch = 0; ch < result.chains.length; ch++){
				line.push({ name:"Chain "+result.chains[ch], filt:result.chains[ch]});
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
	
	for(let li = 0; li < line.length; li++){
		let col = get_colour(li);
		
		let vec = [];
		for(let i = 0; i < result.par_sample.length; i++){
			let samp = result.par_sample[i];
			if(samp.num >= result.burnin){
				if(line[li].filt == "All" || line[li].filt == samp.chain){
					vec.push(get_element(samp.param[th],ind));
				}
			}
		}
	
		let stat = get_statistic(vec);
		let range = get_range(vec);

		let point = [];
				
		switch(rpf.dist_settings.radio.value){
		case "bin":
			{
				let nbin = rpf.dist_settings.sel_bin.te;
				
				let bin = [];
				for(let b = 0; b < nbin; b++) bin[b] = 0;
			
				let dx = (range.max-range.min)/nbin;
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
				range.min -= d*h;
				range.max += d*h;
				
				let ma = range.max, mi = range.min;
				let dd = ma-mi;
				
				let nbin = Math.floor(5/h); 
				if(nbin < 100) nbin = 100; if(nbin > 500) nbin = 500;
				
				let db = dd/nbin;
				
				let hb = Math.floor(1+d*h/db);
				
				let bin=[];
				for(let b = 0; b < nbin; b++) bin[b] = 0;
				
				for(let i = 0; i < vec.length; i++){
					let bf = nbin*(vec[i]-mi)/dd;
					let bmid = Math.floor(bf);
					
					let bmi = bmid-hb, bma = bmid+hb;
				
					for(let b = bmi; b <= bma; b++){
						let bsh = b+0.5;
						
						let d;
						if(bsh < bf) d = 1-(bf-bsh)/hb;
						else d = 1+(bf-bsh)/hb;
						if(d > 0){
							let bb = b;
							if(bb < 0) bb = -bb; 
							if(bb >= nbin) bb = nbin-1 -(bb-(nbin-1));
							bin[bb] += d;
						}		
					}
				}
				
				let sum = 0.0; for(let b = 0; b < nbin; b++) sum += bin[b];
				
				for(let b = 0; b < nbin; b++) bin[b] /= (db*sum);
				
				point.push({x:mi, y:0});
				for(let b = 0; b < nbin; b++){
					point.push({x:mi+((b+0.5)/nbin)*dd, y:bin[b]});
				}
				point.push({x:ma, y:0});
			}
			break;
			
		default: error("Should not be default"); break;
		}
		
		data.push({point:point, col:col, CImin:stat.CImin, CImax:stat.CImax, type:"Distribution"});
		
		key.push({type:"Line", te:line[li].name, col:col});
		
		if(line.length == 1 && rpf.dist_settings.show_mean.check == true){
			data.push({te:"Mean", col:BLACK, thick:1, x:stat.mean, type:"VertLine"});
		}
	}
	
	// Draws the prior
	if(rpf.dist_settings.show_prior.check == true && par.variety != "reparam" && par.variety != "likelihood" && par.kind != "const"){
		let N = 200; // The number of points
		
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
					if(samp.num >= result.burnin){
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
				
				let vari = Math.log(1+cv*cv);   // Works out variables on log scale
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
			
		default: error("prior type not recognised: "+pri.type.te); break; 
		}
		
		if(point.length > 0){
			data.push({point:point, col:BLACK, dash:2, type:"PriorLine"});
		
			key.push({type:"Line", te:"Prior", col:BLACK, dash:2});
		}
	}
	
	inter.graph.define("Distribution","Graph",data,{x_label:name, x_param:true, y_label:"Probability", key:key, yaxis:false});
}


/// Extends a range
function extend_range(ran,frac)
{
	let d = ran.max-ran.min;
	
	let min_new = ran.min - frac*d; if(min_new*ran.min < 0) min_new = 0;
	let max_new = ran.max + frac*d; if(max_new*ran.max < 0) max_new = 0;
	
	return {min:min_new, max:max_new};
}

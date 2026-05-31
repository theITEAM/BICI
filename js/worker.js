"use strict";
// Functions for the worker

importScripts("utils.js");
importScripts("help_text.js");
importScripts("const.js");
importScripts("hash.js");	
importScripts("equation.js");
importScripts("equation_prop.js");
importScripts("equation_calculate.js");
importScripts("parameter.js");
importScripts("parameter_edit.js");
importScripts("prior.js");
importScripts("worker_io.js");
importScripts("model.js");
importScripts("model_indeff.js");
importScripts("data.js");
importScripts("data_check.js");
importScripts("data_file.js");
importScripts("data_update.js");
importScripts("import.js");
importScripts("import_check.js");
importScripts("import_commands.js");
importScripts("annotation.js");
importScripts("output_file.js");
importScripts("results.js");
importScripts("results_plots.js");
importScripts("statistics.js");
importScripts("data_sim.js");
importScripts("create_sire.js");
		
let xi_sq =	calc_xi_sq();

let model;                                        // Stores the model

let model_store;                                  // Stores model in case loading fails

let sim_result = {siminf:"sim"};                  // Stores results from simulation
let inf_result = {siminf:"inf"};                  // Stores results from inference
let ppc_result = {siminf:"ppc"};                  // Stores results from ppc

let sim_result_import = {siminf:"sim"};           // Creates on import and copies once done
let inf_result_import = {siminf:"inf"};
let ppc_result_import = {siminf:"ppc"};

let edit_source;                                  // Used when editting a data source

let input;                                        // Stores the input into the worker
      
const offscreen = new OffscreenCanvas(10,10);     // Offscreen canvas (for measuring text width)
let cv = offscreen.getContext('2d');

let hash_compgl, hash_tragl;                      // Hash tables for comps and tras

let xmlhttp = new XMLHttpRequest();               // Used for reading in files

let map_store = [];                               // Stores annotation maps

let child=[];

let timer = new Timer();                          // Stores timers (for diagnostic purposes)
	
let run_chain = [];

let is_interface = false;                         // Shows not the interface

/*
let test_ch={testtest:[]};
for(let i = 0; i< 10000000; i++){
	test_ch.testtest[i] = Math.random();
}
*/
/*
for(let i = 0; i< 1000000; i++){
	let arr =[];
	for(let j = 0; j < 100; j++){
		arr[j] = Math.random();
	}
	//21 217 016
	 149 217 016

	test_ch.testtest[i] =arr;
}
*/

//let test_ch=
//5 216 972
//5 216 884

let fileReader = new FileReader();

onmessage = function(e) 
{
	if(false){ prr("Worker start"); prr(e.data);}
	//prr("START"+e.data.type);

	if(try_on){
		try {
			process(e);
		}catch(e){
			restore_model();
			post(e);
		}
	}
	else{
		process(e);
	}
}


/// Processes new message
function process(e)
{
	percent(1)
	
	input = e.data;
	
	if(input.type == "Stop Load"){
		fileReader.abort();
		return;
	}
	
	let info = input.info;

	let update_mod = false;

	if(input.type != "Import2") model_store = undefined;

	if(model && input.model){
		copy_strip(input.model,model);
		if(model.warn.length > 0) update_mod = true;
	}
	
	let itype = input.type;
	//prr(itype+" type");
	switch(itype){
	case "Factor reduce":
		model.factor_reduce(info.p,info.cl,info.fac);
		update_mod = true;
		break;
		
	case "Add A pedigree":
		{
			let so = info.edit_source;
			add_A_pedigree(so.table,so.info.p,so.info.i);
			update_mod = true;
		}
		break;
		
	case "EditAPed":
		{
			let p = info.p, i = info.i;
			post({Amat:model.species[p].ind_eff_group[i].A_matrix, p:p, i:i});
		}
		break;
		
	case "Output States": case "Output Params":
		{
			let result = get_result(input.info.name);
			
			let siminf;
			switch(input.info.name){
			case "Simulation": siminf = "sim"; break;
			case "Inference": siminf = "inf"; break;
			case "Post. Simulation": siminf = "ppc"; break;
			}
			
			let type;
			if(input.type == "Output States") type = "state";
			else type = "param";
			
			let te = extract_text_samples(siminf,type,result);
	
			let file_list=[];
			write_file_store(te,"bicifile",file_list,"bicifile");
			post({file_list:file_list});
		}	
		break;
		
	case "Create SIRE":
		create_sire();
		break;
		
	case "DeleteParamPriorConst":
		{
			let par = model.param[info.i];
			par.prior_const_on = false;
			update_mod = true;
		}
		break;
		
	case "DeleteParamDist":
		{
			let par = model.param[info.i];
			if(par.variety != "const" && par.dep.length != 0){
				let list = par_find_list(par);
				par.value = undefined;
				par.set = false;
			}
			par.variety = "normal";
			par.prior = unset_prior();
			par.prior_split_check = {check:false};
			par.defrep_param_list=[];
			par.defrep_ie_list=[];
			par.defrep_fe_list=[];
			par.defrep_warn=[];
			set_ndep_cont(par);
			update_mod = true;
		}
		break;
	
	case "DeleteParamFactor":
		{
			let par = model.param[info.i];
			par.factor = false;
			update_mod = true;
		}
		break;
		
	case "DeleteParamMult":
		{
			let name = model.param[info.i].name;
			
			model.param.splice(info.i,1);
			
			let pf = model.param_factor;
			
			let k = 0; while(k < pf.length && pf[k].f_param.name != name) k++;
			if(k == pf.length) error("Cannot find parameter factor");
			else pf.splice(k,1);
			update_mod = true;
		}
		break;
		
	case "Add Cla":
		{
			let p = info.p;
			
			let res = model.add_classification(p,info.name,info.index,info.op);
		
			if(res.err == false){			
				if(info.op.coord == "latlng" && info.op.default_map == true){
					let mfi = "../Maps/World.json"; if(load_map_fast) mfi = "D:/BICI_nolongerused/Maps/World.json";
					
					let te = load_file_local(mfi);
					let map = load_annotation_map(te);	
					
					post({species:strip_heavy(model.species), feature:map.feature, box:map.box, info:{p:p,cl:model.species[p].ncla-1,default_map:true}});
					return;
				}					 
			}
	
			post({ species:strip_heavy(model.species)});
		}
		break;
		
	case "Delete Cla":
		model.delete_classification(info.p,info.cl);	
		post({ species:strip_heavy(model.species)});
		break;
		
	case "Delete Species":
		model.delete_species(info.p);
		post({ species:strip_heavy(model.species)});
		break;
		
	case "Delete Comp":
		model.delete_compartment(info.p,info.cl,info.c,true);
		post({ species:strip_heavy(model.species)});
		break;
		
	case "Delete Trans":
		model.delete_transition(info.p,info.cl,info.i);
		post({ species:strip_heavy(model.species)});
		break;
		
	case "Init pop convert":
		{
			let so = info;
			data_source_check_error("worker",so);
			init_pop_convert_to_graphical(so,false);
			post({ so:so});
		}
		break;
		
	case "Add comp map":
		{
			let so = info;
			data_source_check_error("worker",so);
			add_compartment_map(so.col,so);
			post({ p:so.info.p, cl:so.info.cl, species:strip_heavy(model.species), map_store:map_store});
		}
		break;
		
	case "Add file trans":
		{
			let so = info;
			data_source_check_error("worker",so);
			model.add_file_transition(so.info.p,so.info.cl,so.table,so.type);
			post({ p:so.info.p, cl:so.info.cl, species:strip_heavy(model.species), warning:so.warning});
		}
		break;

	case "Add file comp":
		{
			let so = info;
			data_source_check_error("worker",so);
			model.add_file_compartment(so.info.p,so.info.cl,so.table,so.info.col);
			post({ p:so.info.p, cl:so.info.cl, species:strip_heavy(model.species), warning:so.warning});
		}
		break;
								
	case "View Code": case "Export BICI": case "Start": case "Save BICI": 
		{
			let one_file = true; if(input.type == "Export BICI") one_file = false;
			create_output_file(info.save_type,one_file,info.map_store);
		}
		break;
		
	case "StartClusterSave": case "StartClusterExport":
		{
			let one_file = true; if(input.type == "StartClusterExport") one_file = false;
		
			switch(info.save_type){
			case "ext": create_ext_file(one_file); break;
			case "ppc": create_ppc_file(info,one_file); break;
			default: create_output_file(info.save_type,one_file,info.map_store); break;
			}
		}
		break;
		
	case "Save Sim": case "Save Inf": case "Save PPC":
	case "Export Sim": case "Export Inf": case "Export PPC":
		{
			let spl = input.type.split(" ");
			let one_file = true; if(begin_str(input.type,"Export")) one_file = false;
			
			let res;
			switch(spl[1]){
			case "Sim": res = sim_result; break;
			case "Inf": res = inf_result; break;
			case "PPC": res = ppc_result; break;
			default: error("op problem"); break;
			}
		
			let bscr = copy(res.bscript);

			if(spl[1] == "Inf") remove_ppc_command(bscr,false);
		
			let file_list=[];
			let te = generate_text_from_bscript(bscr,file_list,one_file);
	
			write_file_store(te,"bicifile",file_list,"bicifile");
			
			post({file_list:file_list});
		}
		break;
		
	case "StartPPC":
		create_ppc_file(info,true);
		break;
		
	case "StartEXT":
		create_ext_file(true);
		break;

	case "Import2":
		import_file2(info.data_file_list);
		break;
		
	case "Spawn Output":
		import_file(info.content,"Execute/init.bici",false);
		break;

	case "Load Default": 
		load_bici("D:/AI temp/test.bici")
		break;
		
	case "Load Example": 
		if(model){
			model.clear_res();
		}
		load_bici("../Examples/"+info+".bici");
		break;
		
	case "Copy simulated":
		copy_simulated_data();
		post({ species:strip_heavy(model.species), inf_details:model.inf_details});
		break;
		
	case "Rename Species":
		model.rename_species(info.species_new,info.p); 
		model.check_ob_string_exist(model,"model",info.species_old);
		post({species:strip_heavy(model.species), param:strip_heavy(model.param),});
		break;
		
	case "Rename Classification":
		model.rename_classification(info.new_name,info.p,info.cl); 
		data_update_rename_classification(info.p,info.cl,info.old_name,info.new_name);
		model.check_ob_string_exist(model,"model",info.old_name);
		post({species:strip_heavy(model.species)});
		break;
		
	case "Rename Index":
		model.rename_index(info.index_new,info.p,info.cl);			
		post({index_old:info.index_old, param_factor:strip_heavy(model.param_factor), param:strip_heavy(model.param), species:strip_heavy(model.species)});
		break;
		
	case "Rename Compartment":
		model.rename_compartment(info.p,info.cl,info.c,info.new_name); 
		equation_rename_compartment(info.p,info.cl,info.old_name,info.new_name);
		param_update_rename_compartment(info.p,info.cl,info.c,info.old_name,info.new_name);
		data_update_rename_compartment(info.p,info.cl,info.old_name,info.new_name);
		model.check_ob_string_exist(model,"model",info.old_name); // Checks not in model
		post({species:strip_heavy(model.species)});
		break;
	
	case "Set Xvector":
		{
			let feg = model.species[info.p].fix_eff[info.i];
			feg.X_vector.X_value = info.value;
			post({});
		}
		break;
		
	case "Edit Xvector":
		{
			let feg = model.species[info.p].fix_eff[info.i];
			info.X_value = feg.X_vector.X_value;
			info.ind_list = feg.X_vector.ind_list;
			post({ info:info});
		}
		break;
		
	case "Load Xvector":
		{
			let so = info.edit_source;
			data_source_check_error("worker",so);
			inf_eff_load_X(so);
			post({p:so.info.p, i:so.info.i});
		}
		break;
		
	case "Set A matrix":
		{
			let ieg = model.species[info.p].ind_eff_group[info.i];
			ieg.A_matrix.A_value = info.value;
			post({});
		}
		break;
		
	case "Edit A matrix":
		{
			let ieg = model.species[info.p].ind_eff_group[info.i];	
			
			let ind_list = ieg.A_matrix.ind_list;
			let A_value = ieg.A_matrix.A_value;
			
			if(ind_list.length*ind_list.length > ELEMENT_MAX){
				let max = Math.floor(Math.sqrt(ELEMENT_MAX));
				
				let ind_list_shrink=[], A_value_shrink=[];
				for(let j = 0; j < max; j++){
					ind_list_shrink.push(ind_list[j]);
					
					A_value_shrink[j]=[];
					for(let i = 0; i < max; i++){
						A_value_shrink[j][i] = A_value[j][i];
					}
				}
				
				info.type = ieg.A_matrix.type;
				info.A_value = A_value_shrink;
				info.ind_list = ind_list_shrink;
				info.too_big = true;
			}
			else{
				info.A_value = A_value;
				info.ind_list = ind_list;
			}
			
			post({ info:info});
		}
		break;
		
	case "Set PriorSplit":
		{
			let i = info.i;
			let par = model.param[i];
			
			let err = check_prior_split(input.type,par,info.prior_split);
			if(typeof err == 'string') alert_help("Problem updating",err);
			
			par.prior_split = info.prior_split;
			par.prior_split_set = true;
			get_prior_param_list(par);
			post({ i:i, prior_split_desc: get_prior_split_desc(par)});
		}
		break;
		
	case "Load PriorSplit":
		{
			let ep = info.ep;
			data_source_check_error("worker",info.source);
			load_priorsplit(ep,info.source,info.dist);
			post({ type:"Load PriorSplit", ep:ep});
		}
		break;
		
	case "Edit PriorSplit": case "Edit DistSplit":
		{
			let par = model.param[info.i];
		
			if(par.prior_split == undefined){
				let prior_split = par_find_template(par.list);
				let dim = get_dimensions(prior_split);
				let ele_list = get_element_list(prior_split,dim);
		
				for(let k = 0; k < ele_list.length; k++){	
					let pr_sp = unset_prior();	
					set_element(prior_split,ele_list[k],pr_sp);
				}
				
				info.prior_split = prior_split;
			}
			else info.prior_split = par.prior_split;
			
			info.value = get_prior_split_tensor(info.prior_split);
			
			info.list = par.list;
			
			if(num_element(par) > ELEMENT_MAX) reduce_size(info,par);
		
			post({ info:info});
		}
		break;
		
	case "View Param":
		{	
			let par = model.param[info.i];
			
			if(par.den_vec){ set_density(info,par,false); par.value = info.value;}
			else info.value = par.value;
			
			let dim = get_dimensions(par.value);
		
			for(let i = 0; i < dim.length; i++){
				if(dim[i] != par.list[i].length){
					alert_help("Too large to plot");
				}
			}
			
			post({ info:info});
		}
		break;
		
	case "View Prior Const":
		{	
			let par = model.param[info.i];
			
			info.value = par.prior_const;
			
			let dim = get_dimensions(par.prior_const);
		
			for(let i = 0; i < dim.length; i++){
				if(dim[i] != par.list[i].length){
					alert_help("Too large to plot");
				}
			}
			post({ info:info});
		}
		break;
		
	case "Edit Param": case "Edit Reparam": case "Edit Define": 
		{
			let par = info.par_st;
			let parm = model.param[info.i];
			par.list = parm.list;
			
			info.value = parm.value;	
			if(info.value == undefined) info.value = param_blank(par);
			
			info.time_dep = par.time_dep;
			
			info.list = parm.list;
			if(par.dist_mat) set_dist(info,par);
			else{
				if(par.iden_mat) set_iden(info,par);
				else{
					if(par.den_vec) set_density(info,par);
					else{
						if(num_element(par) > ELEMENT_MAX) reduce_size(info,par);
					}
				}
			}
			
			post({ info:info});
		}
		break;
		
	case "Edit Prior Const":
		{
			let par = model.param[info.i];
		
			info.value = par.prior_const;	
			if(info.value == undefined) info.value = param_blank(par);
			
			info.list = par.list;
			
			if(num_element(par) > ELEMENT_MAX) reduce_size(info,par);
			
			post({ info:info});
		}
		break;
		
	case "Edit Weight":
		{
			let par = model.param[info.i];
			
			info.value = par.factor_weight;	
			if(info.value == undefined) info.value = param_blank(par);
			
			info.list = par.list;
			if(num_element(par) > ELEMENT_MAX) reduce_size(info,par);
			info.type = "weight";
			post({ info:info});
		}
		break;
		
	case "Set Param": case "Set Reparam": case "Set Define":
		{
			let i = info.i;
			let par = info.par_st;	
			
			if(info.vari_new != undefined) par.variety = info.vari_new;
			par.set = true;
			
			let err = check_param_value(input.type,par,info.value);
			if(typeof err == 'string') alert_help("Problem updating",err);
		
			unescape_param_value(input.type,par,info.value);
			
			par.value = info.value;
		
			get_defrep_param_list(par);
			model.param[i] = par;
			
			update_model();
			post({ par_name:par.name, param:strip_heavy(model.param), species:strip_heavy(model.species)});
		}
		break;
		
	case "Set Prior Const": 
		{
			let i = info.i;
			let par = model.param[i];
			let err = check_param_value(input.type,par,info.value);
			if(typeof err == 'string') alert_help("Problem updating",err);
			
			par.prior_const = info.value;
			par.prior_const_set = true;
	
			//get_defrep_param_list(par);
			
			post({ i:i, prior_const_desc:get_prior_const_desc(par)});
		}
		break;
		
	case "Set Weight": 
		{
			let i = info.i;
			let par = model.param[i];
			let err = check_param_value(input.type,par,info.value);
			if(typeof err == 'string') alert_help("Problem updating",err);
			
			par.factor_weight = info.value;
			post({ i:i, weight_desc:get_weight_desc(par)});
		}
		break;
		
	case "Load Reparam": case "Load Define":
		{
			let ep = info.ep;
			data_source_check_error("worker",info.source);
			load_reparam(ep,info.source,info.dist);
			post({ ep:ep});
		}
		break;
		
	case "Load Tensor":
		{
			let ep = info.ep;
			data_source_check_error("worker",info.source);
			load_tensor(ep,info.source);
		
			post({ ep:ep});
		}
		break;
		
	case "Data Source": 
		if(info.type == "Add" || info.type == "Replace"){
			data_source_check_error("alert warn",info.edit_source);
		}
	
		data_source(info.type,info.edit_source,info.info);
		break;
	
	case "New model":
		new_model();
		break;
		
	case "Sim data":
		simulate_data(info);
		break;
	
	case "Generate const prior":
		generate_const_prior(info.th);
		break;
	
	case "Add ind single":
		add_sing_ind_per_comp(info);
		break;
		
	case "Graph spline":
		{
			let res = info;
			switch(res.siminf){
			case "sim": graph_spline_calculate(sim_result,res.plot_filter,res.burnin); break;
			case "inf": graph_spline_calculate(inf_result,res.plot_filter,res.burnin); break;
			case "ppc": graph_spline_calculate(ppc_result,res.plot_filter,res.burnin); break;
			default: error("Problem sim inf"); break;
			}
		}
		break;
		
	case "Graph pop":
		{
			let res = info;
		
			switch(res.siminf){
			case "sim": graph_pop_calculate(sim_result,res.plot_filter,res.burnin,res.p); break;
			case "inf": graph_pop_calculate(inf_result,res.plot_filter,res.burnin,res.p); break;
			case "ppc": graph_pop_calculate(ppc_result,res.plot_filter,res.burnin,res.p); break;
			default: error("Problem sim inf"); break;
			}
		}
		break;
		
	case "Graph trans":
		{
			let res = info;
			switch(res.siminf){
			case "sim": graph_trans_calculate(sim_result,res.plot_filter,res.burnin,res.p); break;
			case "inf": graph_trans_calculate(inf_result,res.plot_filter,res.burnin,res.p); break;
			case "ppc": graph_trans_calculate(ppc_result,res.plot_filter,res.burnin,res.p); break;
			default: error("Problem sim inf"); break;
			}
		}
		break;
		
	case "Graph Trans. (dist.)": 
		graph_trans_hbin_calculate(inf_result,info.plot_filter,info.burnin,info.p,"Distribution"); 
		break;
	
	case "Graph Trans. (exp.)":
		graph_trans_expect_calculate(inf_result,info.plot_filter,info.burnin,info.p,input.type); 
		break;
	
	case "Graph Trans. (bias)": case "Graph Trans. (p-val.)":
		{
			let rpf = info.plot_filter;
			graph_trans_hbin_calculate(inf_result,rpf,info.burnin,info.p,rpf.sel_bias_view.te);
		}
		break;
		
	case "Graph ind":
		{
			let res = info;
			switch(res.siminf){
			case "sim": graph_ind_calculate(sim_result,res.plot_filter,res.burnin,res.p); break;
			case "inf": graph_ind_calculate(inf_result,res.plot_filter,res.burnin,res.p); break;
			case "ppc": graph_ind_calculate(ppc_result,res.plot_filter,res.burnin,res.p); break;
			default: error("Problem sim inf"); break;
			}
		}
		break;
		
	case "Graph ind sing":
		{
			let res = info.res;
			let ind_sel = info.ind_sel;
			
			switch(res.siminf){
			case "sim": graph_ind_sing_calculate(ind_sel,sim_result,res.plot_filter,res.burnin,res.p); break;
			case "inf": graph_ind_sing_calculate(ind_sel,inf_result,res.plot_filter,res.burnin,res.p); break;
			case "ppc": graph_ind_sing_calculate(ind_sel,ppc_result,res.plot_filter,res.burnin,res.p); break;
			default: error("Problem sim inf"); break;
			}
		}
		break;
		
	case "Graph generation":
		graph_generation_calculate(inf_result,info.plot_filter,info.burnin); 
		break;
		
	case "Graph param":
		{
			let res = info;
			switch(res.siminf){
			case "sim": graph_param_calculate(sim_result,res.plot_filter,res.burnin); break;
			case "inf": graph_param_calculate(inf_result,res.plot_filter,res.burnin); break;
			case "ppc": graph_param_calculate(ppc_result,res.plot_filter,res.burnin); break;
			default: error("Problem sim inf"); break;
			}
		}
		break;
	
	case "Graph view":
		create_view_graph_calculate(info.name,info.prior_const,info.sel_view,info.so);
		break;
	
	case "Load File":
		load_fi(info.file,info.type,info.heading,info.format,info);
		break;
	
	case "Import output":
		{
			let fi = "M:/Github/theITEAM/BICI/Execute/init.bici";
			if(ver =="windows"){
				fi = "C:/Users/cpooley/Desktop/BICI_release/BICI_v0.8_windows/Execute/init.bici";
			}
			if(ver =="mac"){
				fi = "/tmp/init.bici";
			}
			percent(2)
			load_bici(fi);
		}
		break;

	case "Import output2":
		load_bici("D:/Execute/init.bici");
		break;
		
	case "Combine IE":
		add_ie_to_group(info.p,info.i,info.e);
		update_mod = true;
		break;
		
	case "Split IE":
		split_ie(info.p,info.i);
		update_mod = true;
		break;
	
	case "Get Ind List":
		{
			let name_list =[];
			switch(info.ty){
			case "generate":                             // Gets potential individuals from simulation
				{
					let indi = sim_result.sample[0].species[info.p].individual;
					for(let i = 0; i < indi.length; i++){
						name_list.push({name:sim_result.all_ind_list[indi[i].all_ind_ref].name});
					}
				}
				break;
				
			case "data":                                 // Gets potential individuals from data
				name_list = model.get_all_data_individual(info.p);
				break;
			}
			
			post({ name_list:name_list});
		}
		break;
		
	case "Gen. test-and-cull":
		generate_test_and_cull(info.p,info.i,info.sel_sim);
		break;
		
	case "UpdateModel":
		update_mod = true;
		break;
		
	case "CopyValue":
		copy_model_value();
		update_mod = true;
		break;
		
	case "Stop":
		break;
		
	case "PosteriorMean":
		sim_param_post_mean();
		break;
	
	default:
		error("WORKER PROBLEM:"+input.type);
		break;
	}
	
	if(update_mod == true){
		update_model();

		let do_after;
		if(input.type == "UpdateModel") do_after = input.info;

		post({ type:"UpdateModel", do_after:do_after, param_factor:strip_heavy(model.param_factor), param:strip_heavy(model.param), species:strip_heavy(model.species)});
	}
}


/// If there are any invalid data sources then removes and generates a message
function create_invalid_message(mess)
{
	if(model && model.species){
		let list = [];
		
		for(let loop = 0; loop < 3; loop++){
			let siminf = get_siminf(loop);
			let nsp = get_so_nsp(siminf);
			if(nsp != undefined){
				for(let p = 0; p < nsp; p++){	
					let source = get_source(siminf,p);
					let type;
					switch(loop){
					case 0: type = "Population"; break;
					case 1: type = "Data"; break;
					case 2: type = "Population Mod."; break;
					}
					
					let j = 0; 
					while(j < source.length){
						let so = source[j];
						if(so.error == true){
							list.push("• <b>"+type+"</b> - "+so.name);
							source.splice(j,1);
						}
						else j++;
					}
				}
			}
		}
		
		if(list.length > 0){
			let te = "Due to changes in the model these data sources have been removed:\n";
			for(let j = 0; j < list.length; j++){
				te += list[j]+"\n";
			}
			mess.invalid_message = te;
			mess.species = strip_heavy(model.species);
		}
	}
}


/// Posts a message back to worker_reply
function post(mess)
{
	if(mess.type == undefined) mess.type = input.type;
	create_invalid_message(mess);
		
	postMessage(mess);
}


/// Outputs percentage complete
function percent(per)
{
	post({ type:"Percent", info:per});
}


/// Outputs percentage complete
function loading_mess(te)
{
	post({ type:"LoadMess", info:te});
}

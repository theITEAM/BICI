"use strict";
// Functions which receive reply from the worker

worker.onmessage = function (e) 
{
	let ans = e.data;
	
	//prr("Worker ans"); prr(ans);
	
	if(ans.type == "Percent"){
		let mem = performance.memory;
		set_loading_percent(ans.info);
		return;
	}
	
	if(ans.type == "LoadMess"){
		loading_symbol_message(ans.info);
		return;
	}

	if(false){ prr("Worker reply: "+ans.type); prr("ans"); prr(ans);}

	if(inter.worker_mess.active != "stop"){
		if(ans.type != "Start" && ans.type != "StartPPC" && ans.type != "Import model files") stop_loading_symbol();
		
		switch(ans.type){		
		case "EditAPed":
			edit_A_pedigree(ans.Amat,ans.p,ans.i);
			break;
		
		case "Create SIRE":
			write_file_async(ans.te,"../ZSIRE2.1/file.sire","export");
			write_file_async(ans.dt,"../ZSIRE2.1/table.csv","export");
			write_file_async(ans.Ate,"../ZSIRE2.1/A.csv","export");
			generate_screen();
			break;
		
		case "Delete Cla": case "Add Cla":
		case "Delete Species": case "Delete Comp": case "Delete Trans":
			model.species = ans.species;	
		
			update_param();
			close_bubble();
			close_help();	
			initialise_pages();
			
			switch(ans.type){
			case "Delete Species":
				{
					let p = model.get_p();
					if(p == model.species.length && p > 0){
						change_page({pa:"Model",su:"Compartments",susu:p-1});
					}
				}
				break;
				
			case "Add Cla":
				{
					if(ans.feature != undefined) load_map(ans);		
				
					let p = model.get_p();
					change_page({sususu:model.species[p].ncla-1});
				}
				break;
			}
			generate_screen();
			break;
			
		case "Select Element":
			generate_screen();
			select_table_elelent(ans.r,ans.c);
			inter.bubble.error_warning = ans.warn;
			generate_screen();
			break;
			
		case "Load Comp Map":
			{
				let p = inter.p_cl_store.p, cl = inter.p_cl_store.cl;
				start_data_source("CompMap",{},{p:p, cl:cl});
				data.table.push(ans.tab);

				let so = edit_source;
				so.table_loaded = true;
				so.p = p; so.cl = cl;

				so.data_table_use = data.table.length-1;

				transfer_column(0);
				generate_screen();
			}
			break;
			
		case "View Code":
			{
				inter.help = { title:"Source code", te:"This shows the source code when exported in .bici format.", st:"", scroll_to_line:false, script:ans.formatted};
				model_updated(ans);
			}
			break;
			
		case "Save BICI":	case "Export BICI": case "StartCluster":
			inter.file_store.file_list = ans.file_list;
			inter.file_store.save_type = ans.save_type;
			
			if(inter.file_store.type == "Save") save_bici(inter.file_store.filename)
			else saving_dialogue("",".bici","BICI_file");
			model_updated(ans);
			break;
		
		case "Save Sim": case "Save Inf": case "Save PPC":
			inter.file_store.file_list = ans.file_list;
			inter.file_store.save_type = "save";
			saving_dialogue("",".bici","BICI_file");
			break;
		
		case "Output States": case "Output Params":
			inter.file_store.file_list = ans.file_list;
			inter.file_store.save_type = "save";
			let ty = ".txt"; if(ans.type == "Output Params") ty = ".csv";
			saving_dialogue("",ty,"BICI_file");
			close_bubble();
			generate_screen();
			break;
			
		case "Copy simulated":
			model.inf_details = ans.inf_details;
			model.species = ans.species;
			change_page({pa:"Inference", su:"Data", susu:model.get_p()});
			update_param();
			break;
			
		case "Data Source":	
			data_source_reply(ans);
			update_param();
			generate_screen();
			break;
		
		case "Rename Species":
			model.species = ans.species;
			update_param();
			close_bubble();
			initialise_pages();
			generate_screen();
			break;
			
		case "Rename Classification":
			model.species = ans.species;
			update_param();
			close_bubble();
			initialise_pages();
			generate_screen();
			break;
			
		case "Rename Index":
			model.species = ans.species;
			model.param	= ans.param;
			model.get_label_info_all();
			model.check_ob_string_exist(model,"model",ans.index_old);// Checks not in model
			update_param();
			close_bubble();
			generate_screen();
			break;
			
		case "Rename Compartment":
			model.species = ans.species;
			update_param();
			close_bubble();
			generate_screen();
			break;
			
		case "Set Xvector":
			close_bubble();
			close_Xvector_source();
			generate_screen();
			break;
			
		case "Edit Xvector":
			inter.edit_Xvector = ans.info;
			generate_screen();
			break;
			
		case "Load Xvector":
			model.species[ans.p].fix_eff[ans.i].X_vector.loaded = true;
			close_data_source();
			generate_screen();
			break;
			
		case "Set A matrix":
			close_bubble();
			close_Amatrix_source();
			generate_screen();
			break;
			
		case "Edit A matrix":
			inter.edit_Amatrix = ans.info;
			generate_screen();
			break;
			
		case "Set PriorSplit":
			model.param[ans.i].prior_split_desc = ans.prior_split_desc;
			model.param[ans.i].prior_split_set = true;
			close_bubble();
			close_param_source();
			update_param();
			generate_screen();
			break;
		
	 case "Load Reparam": case "Load Tensor":
			{
				inter.edit_param = ans.ep;
				let par = model.param[ans.ep.i];
				par.value_desc = ans.ep.value_desc;
				par.set = ans.ep.set;
				close_data_source();
				if(ans.type == "Load Reparam") update_param();
				generate_screen();
			}
			break;
			
		case "Load PriorSplit":
			{
				inter.edit_param = ans.ep;
				let par = model.param[ans.ep.i];
				par.prior_split_desc = ans.ep.prior_split_desc;
				par.prior_split_set = ans.ep.prior_split_set;
				close_data_source();
				generate_screen();
			}
			break;
			
		case "Edit PriorSplit": case "Edit DistSplit":
			inter.edit_param = ans.info;
			generate_screen();
			break;
			
		case "View Param":
			{
				inter.view_graph = ans.info;	
				let pos_view = inter.view_graph.pos_view;
				let j = 0;
				while(j < pos_view.length){
					if(pos_view[j].te == "Graph (split)") pos_view.splice(j,1);
					else j++;
				}
				generate_screen();
			}
			break;
			
		case "Edit Param": case "Edit Reparam": case "Edit Weight":
			inter.edit_param = ans.info;
			generate_screen();
			break;
			
		case "Set Param": case "Set Reparam":
			inter.edit_source = false;
			close_bubble();
			close_param_source();
			
			model.param[ans.i].value_desc = ans.value_desc;
			model.param[ans.i].set = true;
			if(ans.type == "Set Reparam") update_param();
			generate_screen();
			break;
			
		case "Set Weight": 
			model.param[ans.i].weight_desc = ans.weight_desc;
			close_param_source();
			generate_screen();
			break;
			
		case "Close Data Source":
			close_data_source();
			generate_screen();
			break;
			
		case "Edit Data Source":
			edit_source = ans.edit_source;
			edit_source.edit_spec = false;
			generate_screen();
			break;
			
		case "Load map":
			load_map(ans);		
			generate_screen();
			break;
			
		case "New model":
			model.load(ans);
			inter.file_store.filename = "";
			map_store=[];
			break;
			
		case "Sim data":
			{	
				close_bubble();
				add_gen_data(ans.head,ans.ele);
				generate_screen();
			}
			break;
				
		case "Graph define":
			{
				inter.graph.define(ans.variety,ans.view,ans.data,ans.op);
				inter.graph.init = true;
				generate_screen();
			}
			break;
		
		case "Load Example":
			model.load(ans);
			break;
		
		case "Import model files":
			load_BICI_files(ans,15,30);
			break;
		
		case "Import2":
			model.load(ans);	
			break;
		
		case "Load File":
			switch(ans.info){
			case "BICI_file": case "Import file":
				model.load(ans);	
				break;
			
			case "Data file":
				data.table.push(ans.tab);
				edit_source.data_table_use = data.table.length-1;
				edit_source.table_loaded = true;
				close_help();				
				generate_screen();
				break;
				
			case "A matrix": case "Ainv matrix": 
				model.species = ans.species;
				close_help();
				generate_screen();
				break;
			}
			break;
		
		case "Spawn Output":
			model.load(ans,true);	
			break;
		
		case "Import output": case "Import output2": case "Load Default": 
			model.load(ans);
			break;
		
		case "Add comp map":
			{
				model.species = ans.species;
				for(let i = 0; i < ans.map_store.length; i++){
					map_store.push(ans.map_store[i]);
				}
			}
			set_camera(ans.p,ans.cl);
			update_param();
			close_data_source();
			generate_screen();
			break;
			
		case "Init pop convert":
			edit_source = ans.so;
			generate_screen();
			break;
			
		case "Add file trans":
			model.species = ans.species;
			set_camera(ans.p,ans.cl);
			update_param();
			close_data_source();
			generate_screen();
			break;
			
		case "Add file comp":
			model.species = ans.species;
			set_camera(ans.p,ans.cl);
			update_param();
			close_data_source();
			if(ans.warning && ans.warning != "") alert_help("Warning",ans.warning);
			generate_screen();
			break;
			
		case "UpdateModel":
			model_updated(ans);
			break;
			
		case "Start": case "StartPPC":
			{
				model_updated(ans);
		
				inter.save_type = ans.save_type;
				if(check_memory(ans)){
					start_spawn(ans.file_list);
				}
			}
			break;
			
		case "Model warning":
			{	
				if(ans.full_warn) model.warn_view = true;
				model.warn = ans.warn;
				model.species = ans.species;
				close_help();
				generate_screen();
			}
			break;
		
		case "Alert State Memory":
			{
				let warn=[];
				warn.push("Unfortunately the BICI file could not be loaded because the memory required to store state samples has run out (a maximuim of "+mem_state_sample_max+" GB is allocted for this)."); 
	
				inter.help = {title:"Memory has run out", te:"", sug:state_memory_error(ans.siminf), warn:warn};
				stop_loading_symbol();
				generate_screen();
			}
			break;
			
		case "Alert Param Memory":
			{
				let warn=[];
				warn.push("Unfortunately the BICI file could not be loaded because the memory required to store parameter samples has run out (a maximuim of "+mem_param_sample_max+" GB is allocted for this)."); 
	
				inter.help = {title:"Memory has run out", te:"", sug:param_memory_error(ans.siminf), warn:warn};
				stop_loading_symbol();
				generate_screen();
			}
			break;
			
		case "Alert Import":
			inter.help = { title:ans.title, te:ans.te, st:ans.st, line:ans.line, scroll_to_line:ans.scroll_to_line, script:ans.script};
			add_script_width(inter.help.script);
			generate_screen();
			break;
			
		case "AlertP":
			alertp(ans.te);
			break;
			
		case "AlertHelp":
			alert_help(ans.title,ans.te);
			generate_screen();
			break;
			
		case "Get Ind List":
			ind_group_data2(ans.name_list);
			break;
			
		case "Error":
			error("Error: "+ans.te);
			break;
			
		default:
			error("Option not recognisted: "+ans.type);
			break;
		}
	}
	else error("ignore reply");
	
	if(ans.invalid_message){ 
		model.species = ans.species;
		alert_help("Invalid data sources",ans.invalid_message);
		generate_screen();
	}
	
	inter.worker_mess.active = false;
	
	if(inter.worker_mess.next){ setTimeout(worker_mess_process,10);}
};


function model_updated(ans)
{
	model.param = ans.param;
	model.param_factor = ans.param_factor;
	model.species = ans.species;
	model.warn.length = 0;
	model.warn_view = true;
	model.get_label_info_all();	
	generate_screen();
}

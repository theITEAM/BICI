"use strict";

/// Saves the model/analysis in a special "BICI" format (this uses JSON to store the various objects)
function save_BICI_file(with_results) 
{
	let st = "";
	
	st += save_add("start",model.start);
	st += save_add("description",model.description);
	st += save_add("species",model.species);
	st += save_add("param",model.param);
	st += save_add("derive",model.derive);
	st += save_add("sim_details",model.sim_details);
	st += save_add("inf_details",model.inf_details);
	st += save_add("data_table",data.table);
	st += save_add("data_source",data.source);
	//st += save_add("sim_result",sim_result);
	//st += save_add("inf_result",inf_result);
	
	
	return st;
}


/// Adds an object to be saved
function save_add(te,ob)
{
	return te+":"+JSON.stringify(ob)+"\n";
}


/// Loads a BICI file
function load_BICI_file(te)                                          
{		
	model.start_new();
	
	te = te.replace(/′/g,"'");
	
	let lines = te.split('\n');

	for(let n = 0; n < lines.length; n++){  // Removes and '\r' characters
		if(lines[n].substr(lines[n].length-1,1) == "\r") lines[n] = lines[n].substr(0,lines[n].length-1);
	}
	
	for(let n = 0; n < lines.length; n++){
		let st = lines[n];
		if(st.length > 0){
			let j = 0; while(j < st.length && st.substr(j,1) != ":") j++; 
			if(j == st.length){ alertp("Problem loading");}
			
			let code = st.substr(0,j); st = st.substr(j+1); 
			if(st != "undefined"){
				switch(code){
				case "start": model.start = JSON.parse(st); break;
				case "description": model.description = JSON.parse(st); break;
				case "species": model.species = JSON.parse(st); break;
				case "param": model.param = JSON.parse(st); break;
				case "derive": model.derive = JSON.parse(st); break;
				case "sim_details": model.sim_details = JSON.parse(st); break;
				case "inf_details": model.inf_details = JSON.parse(st); break;
				case "data_table": data.table = JSON.parse(st); break;
				case "data_source": data.source = JSON.parse(st); break;
				case "sim_result":
					//sim_result = JSON.parse(st); 
					break;
				case "inf_result": 
					//inf_result = JSON.parse(st); 
					break;
				default: error("Loading '"+code+"' not recognised"); break; 
				}
			}
		}
	}
	
	/*
	{ // TO DO TURN OFF
			model.sim_details.indmax= 10000;
			model.sim_details.algorithm = {value:"gillespie"};
			model.sim_details.number = 1;
			
			model.inf_details.abcsample = "1000";
			model.inf_details.sample = "10000";
			model.inf_details.thinparam = "10";
			model.inf_details.thinstate ="50";
			model.inf_details.accfrac ="0.1";
			model.inf_details.accfracsmc = "0.5";
			model.inf_details.numgen = "5";
			model.inf_details.kernelsize = "0.5";
			model.inf_details.indmax = 10000;
			model.inf_details.algorithm = {value:"DA-MCMC"}; 
	}
	*/
	model.inf_details.nchain = "3";
	
	sim_result = clear_result();
	inf_result = clear_result();
	
	model.determine_branching();

	check_data_valid();
}


/// Loads one of the example files
function load_example(val)                                   
{	 
	start_loading_symbol(0);
	
	let file = "Examples/"+val+".bici";
	//let file = "Testing/"+val+".bici";
	model.filename = file;
	
	//let xmlhttp = new XMLHttpRequest();  
	xmlhttp.abort();
	
	xmlhttp.addEventListener("progress", function(ev) { 
		inter.loading_symbol.percent = Math.floor(100*ev.loaded/ev.total);
	});
	
	xmlhttp.addEventListener("loadend",function(){
		stop_loading_symbol();
	
		let te = xmlhttp.responseText;

		load_BICI_file(te);
		
		initialise_pages();
		
		change_page({pa:"Model", su:"Compartments"});
		//change_page({pa:"Simulation", su:"Run"});
		//change_page({pa:"Inference", su:"Run"});
	});
  
	xmlhttp.open("GET",file,true); xmlhttp.send();
}


/// Opens a load file dialogue
function loading_dialogue(value,filter,type)                                       
{
	ById("fileToLoad").value = value;
	ById("fileToLoad").accept = filter; 
	inter.file_type = type;
	ById("fileToLoad").click();	
}


/// Opens a save file dialogue
function saving_dialogue(value,filter,type)                                       
{
	ById("fileToSave").value = value;
	ById("fileToSave").accept = filter; 
	inter.file_type = type;
	ById("fileToSave").click();	
}


/// Loads a file
function load_file()
{
	close_bubble(); generate_screen();
	
	let value = ById("fileToLoad").value;	
	let fileToLoad = ById("fileToLoad").files[0];
	if(!fileToLoad){ alertp("Could not load this file."); return;}
	
	let filename = fileToLoad.name;

	data.filename = filename;
	
	var fileReader = new FileReader();
	fileReader.onload = function(fileLoadedEvent){	
		stop_loading_symbol();
	
		let te = fileLoadedEvent.target.result;
	
		te = te.replace(/\r/g, "");
	
		let error = false;
		
		switch(inter.file_type){
		case "Observation_File": 
			/*
				res = load_observation_file(); 
				if(res == 0){
					datashow = "table";
					dataselected = data.length;
					ncoldef = 0;
					initialise_datatemp();
					if(ncol < ncoldefmax){ helptype = 112; return 2;}
					addingdata = 1; 
				}
				*/
			break;
		
		//case "Import_File": res = import_file(textFromFile); break;
		
		//case "InitPopulation_File": res = load_simulation_init(textFromFile); break;

		case "GEOJSON_annotation":
			load_annotation_map(te,model.get_p(),model.get_cl());
			break;
			
		case "GEOJSON_compartment":
			load_compartment_map(te);
			break;
		
		case "BICI_file": 
			load_BICI_file(te); 
			initialise_pages();
			change_page({pa:"Model"});
			break;
			
		case "Import file": 
			import_file(te); 
			break;
		
		case "A matrix":
			{
				let sep = "comma";
				if(inter.radio_format.value == "tsv") sep = "tab"; 
				let tab = load_table(te,true,sep,data.filename);
			
				if(typeof tab == 'string'){
					alertp(tab);
					return;
				}
				else{
					A_matrix_loaded(tab);
				}
				generate_screen();
			}
			break;
			
		case "Data file":
			{
				let head = true;
				let sep = "comma";
				if(inter.radio_heading.value == "No") head = false;
				if(inter.radio_format.value == "tsv") sep = "tab"; 

				let tab = load_table(te,head,sep,data.filename);
		
				if(typeof tab == 'string'){
					alertp(tab);
					return;
				}
				else{
					data.table.push(tab);
					inter.edit_source.data_table_use = data.table.length-1;
					inter.edit_source.table_loaded = true;
					close_help();
				}
				generate_screen();
			}
			break;
			
		default: alertp("File type not recognised"); break;
		}
	
		//res = process_loaded_file();
	/*
		loading = 0; 
		if(res == 1){ 
			if(errormsg != "") helptype = 100; 
			else{ alertp("There was a problem loading this file");}
		}
		*/
		
		generate_screen();
		//button_initialise();
	};
	fileReader.readAsText(fileToLoad, "UTF-8");
	start_loading_symbol();
	
	ById("fileToLoad").value="";
}


/// Loads a table from a file
function load_table_from_file(file)
{
	if(check_char_allowed(file,"<>:\"|?*") != "success") return; 

	file = inter.imp.datadir+"/"+file;
	
	if(check_file_exists(file) != "success") return;
	
	let fs = require('fs');	
	let te = fs.readFileSync(file, 'utf8');

	let sep = "comma"; if(file.substr(file.length-4,4) == ".tsv") sep = "tab";
	
	return load_table(te,true,sep,file);
}


/// Loads a table from text
function load_table(te,head,sep,filename)
{
	te.replace(/\r/g,"");
	let lines = te.split("\n");
	
	while(lines.length > 0 && lines[0].length > 0 && lines[0].substr(0,1) == "#") lines.splice(0,1);

	let ele = [];
	for(let r = 0; r < lines.length; r++){
		let st = lines[r];
		if(st.length > 0){
			let li=[];
			
			switch(sep){
			case "tab":
				li = st.split("\t");
				break;
			
			case "comma":
				{
					let quote = "";
					let ist = 0;
					for(let i = 0; i < st.length; i++){
						let ch = st.substr(i,1);
						switch(ch){
						/*
						case "'":
							if(quote == "'") quote = "";
							else{
								if(quote == "") quote = "'";
							}
							break;
						*/
						case "\"":
							if(quote == "\"") quote = "";
							else{
								if(quote == "") quote = "\"";
							}
							break;
						
						case ",":
							if(quote == ""){
								li.push(st.substr(ist,i-ist));
								ist = i+1;
							}
							break;
						}
					}
					
					if(ist < st.length) li.push(st.substr(ist,st.length-ist));
				}
				break;

			default: error("Option not recognised 62"); break;
			}
			ele.push(li);
		}
	}
	
	// Removes spaces and quotation marks
	for(let r = 0; r < ele.length; r++){
		for(let c = 0; c < ele[r].length; c++){ 	
			let te = ele[r][c].trim();
			if(te.length >= 2){
				if((te.substr(0,1) == "'" && te.substr(te.length-1,1) == "'") ||
				   (te.substr(0,1) == "\"" && te.substr(te.length-1,1) == "\"")){
					te = te.substr(1,te.length-2).trim();
				}						 
			}
			ele[r][c] = te;
		}
	}
	
	if(ele.length == 0) return "No content was found";
		
	let ncol = ele[0].length;
	
	for(let r = 0; r < ele.length; r++){
		if(ele[r].length != ncol) return "Not all the columns have the same size";
	}
	
	let heading=[];
	
	if(head == true){
		heading = ele[0];
		ele.splice(0,1);
	}
	else{
		for(let c = 0; c < ncol; c++) heading.push("Column "+(c+1));
	}
	
	return { filename:cut_path(filename), heading:heading, ele:ele, ncol:ncol, nrow:ele.length};
}
 
 
/// Removes path from filename
function cut_path(filename)
{
	let j = filename.length-1;
	while(j >= 0 && filename.substr(j,1) != "/" && filename.substr(j,1) != "\\") j--;
	return filename.substr(j+1);
}


/// Gets a subtable based on a series of column heading names
function get_subtable(tab,col_name)
{
	let col = [];
	for(let i = 0; i < col_name.length; i++){
		let c = find_string_in(tab.heading,col_name[i]);
		if(c == undefined){
			return {error:"Cannot find heading '"+col_name[i]+"' in the file '"+tab.filename+"'"};
		}
		col.push(c);
	}
	
	let ele = [];
	for(let r = 0; r < tab.nrow; r++){
		ele[r]=[];
		for(let i = 0; i < col.length; i++){
			ele[r][i] = tab.ele[r][col[i]];
		}
	}		
	
	return {error:"", filename:tab.filename, heading:col_name, nrow:tab.nrow, ncol:col.length, ele:ele};
}
 

/// Saves an exported file
function save_file(filename,type)                               
{
	close_bubble(); generate_screen();
	
	let te;

	switch(type){
	case "BICI_file": 
		model.filename = filename;
		te = save_BICI_file();	
		break;
	
	case "Export script":
		create_output_file("siminf",filename);
		return;
		
	case "Export graph":
		inter.graph.export_image(filename);
		return;
		
	case "Export table":
		export_table(filename);
		return;
		
	case "Export matrix table":
		export_matrix_table(filename);
		return;
		
	default: error("Option not recognised 63"); break;
	}
	
	write_file(te,filename);
}


/// Prints an image
function print_image(outcan)
{
	let dataURL = outcan.toDataURL();

	let html  = '<html><head><title></title></head>';
	html += '<body style="width: 100%; padding: 0; margin: 0;">';
	html += '<img onload="window.print(); window.close();" ';
	html += 'style="width:95%"';
	html += ' src="' + dataURL + '"  /></body></html>';
	html += '</body></html>';
	
	let printWindow = window.open('', 'to_print', 'left=100, top=100, height=500,width=800');

	printWindow.document.open();
	printWindow.document.write(html);
}


/// Saves an image from a canvas
function save_image(outcan,filename)
{
	let dataURL = outcan.toDataURL();

	dataURL = dataURL.replace(/^data:image\/\w+;base64,/, "");
	
	let buf = new Buffer(dataURL, 'base64');
	
	write_file(buf,filename,"image");
}
	

/// Exports a table in csv format
function export_table(filename)
{
	let l = find(inter.layer,"name","TableContent");
	if(l == undefined){ error("Not layer found"); return;}
	
	let lay = inter.layer[l];
	
	let tab = lay.op.table;

	let te = "";
	for(let i = 0; i < tab.heading.length; i++){
		if(i != 0) te += ",";
		te += add_quote(tab.heading[i].name);
	}
	te += endl;
		
	for(let j = 0; j < tab.content.length; j++){
		let ele =  tab.content[j];
		for(let i = 0; i < ele.length; i++){
			if(i != 0) te += ",";
			let val = ele[i].te;
			if(typeof val == 'string'){
				if(val.includes(" ")) te += add_quote(val);
				else te += val;
			}
			else te += val;
		}
		te += endl;
	}
	
	write_file(te,filename,"export");
}


/// Exports the matrix table in csv format
function export_matrix_table(filename)
{
	if(!inter.graph){ error("Should be graph"); return;}
	
	let da = inter.graph.data[0];
	
	let te = "";
	for(let i = 0; i < da.xlab.length; i++) te += ","+da.xlab[i];
	te += endl;
	
	for(let j = 0; j < da.ylab.length; j++){
		te += da.ylab[j];
		for(let i = 0; i < da.xlab.length; i++) te += ","+da.mat[j][i];
		te += endl;
	}

	write_file(te,filename,"export");
}


/// Writes a file
function write_file(te,filename,op)
{
	if(filename.includes("NO_SAVE")) return;

	let fs = require('fs');
	// The "\ufeff" gets uft8 working in Excell
	//if(op == "export") te = "\ufeff"+te;
	if(op != "image") te = "\ufeff"+te;

	fs.writeFile(filename, te, function (err) {
		 if (err) {
			 alertp("There was an error attempting to save your data.");
			 return;
		 } else if (callback) { callback();}});
}
"use strict";
// Function related to input and output on the webpage sidebar

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
	let value = ById("fileToLoad").value;	
	let fileToLoad = ById("fileToLoad").files[0];

	if(!fileToLoad){ alertp("Could not load this file."); return;}
	
	data.filename = fileToLoad.name;

	ById("fileToLoad").value="";
	
	let op ={ file:fileToLoad, type:inter.file_type};
	if(inter.radio_heading) op.heading = inter.radio_heading.value;
	if(inter.radio_format) op.format = inter.radio_format.value;
	
	if(inter.file_type == "A matrix" || inter.file_type == "Ainv matrix"){
		op.p = edit_source.info.p;
		op.i = edit_source.info.i;
	}
	
	if(inter.file_type == "GEOJSON_annotation"){
		op.p = inter.p_cl_store.p;
		op.cl = inter.p_cl_store.cl;
	}
	
	op.path = fileToLoad.path;
	
	start_worker("Load File",op);
}


/// Saves an exported file
function save_file(filename,type)                               
{
	close_bubble(); generate_screen();
	
	switch(type){
	case "BICI_file":
		save_bici(filename);
		break;
	
	case "Export script":
		start_worker("Export BICI",{save_type:"save", map_store:map_store});
		return;
		
	case "Export graph":
		inter.graph.export_image(image_scale_factor,filename);
		return;
		
	case "Export video":
		inter.graph.export_video(filename);
		return;
		
	case "Export figure":
		inter.graph.export_figure(filename);
		return;
		
	case "Export table":
		export_table(filename);
		return;

	case "Export table content":
		export_table_content(filename);
		return;
		
	case "Export matrix table":
		export_matrix_table(filename);
		return;
		
	default: error("Option not recognised 63"); break;
	}
}


/// Saves BICI file
function save_bici(filename)
{
	let file_list = inter.file_store.file_list;
	
	let dir;
	//if(filename != undefined && file_list.length > 1 && inter.file_store.save_type == "export"){
	if(filename != undefined && file_list.length > 1){
		let root = find_root(filename);
		let file_local = find_file(filename);
		
		let k = file_local.length-1; while(k > 0 && file_local.substr(k,1) != ".") k--;
		if(k == 0) error("Cannot find '.'");
		else dir = file_local.substr(0,k)+"-data-files";
	
		dir = dir.replace(/\\/g,"/");
		
		let te = file_list[file_list.length-1].data;
		let na= 'data-dir folder="."';
		let i = te.indexOf(na);
		if(i == -1) error("cannot find 'data-dir'");
		else{
			i += na.length-2;
			
			file_list[file_list.length-1].data = te.substr(0,i)+dir+te.substr(i+1);
		}
		
		root = root.replace(/\\/g,"/");
		dir = root+dir;
	}
	
	inter.file_store.filename = filename;
	
	create_files(file_list,filename,dir,"save bici");
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
	
	write_file_async(buf,filename,"image");
}
	

/// Exports a table in csv format
function export_table(filename)
{
	let l = find(inter.layer,"name","TableContent");
	if(l == undefined){ error("No layer found"); return;}

	let tab = inter.layer[l].op.table;

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

	te = add_escape_char(te);
	te = te.replace(/—/g,"-");
	te = te.replace(/⟨/g,"<");
	te = te.replace(/⟩/g,">");
	
	write_file_async(te,filename,"export");
}


/// Exports a table (type used when click "view" in csv format
function export_table_content(filename)
{
	let pos = ["TableContent","CreateEditXvectorContent","CreateEditTableContent",
	           "CreateEditParamContent", "CreateEditAmatrixContent"];
	
	let loop = 0; 
	while(loop < pos.length && find(inter.layer,"name",pos[loop]) == undefined) loop++;
	
	if(loop == pos.length){ error("No layer found"); return;}
	
	let l = find(inter.layer,"name",pos[loop]);
	
	let te = "";

	switch(inter.layer[l].name){
	case "TableContent":
		export_table(filename);
		return;
		
	case "CreateEditTableContent":
		{
			let tab = inter.layer[l].op.table;

			for(let i = 0; i < tab.heading.length; i++){
				if(i != 0) te += ",";
				te += add_quote(tab.heading[i]);
			}
			te += endl;
				
			for(let j = 0; j < tab.ele.length; j++){
				let ele = tab.ele[j];
				for(let i = 0; i < ele.length; i++){
					if(i != 0) te += ",";
					let val = ele[i];
					if(typeof val == 'string'){
						if(val.includes(" ")) te += add_quote(val);
						else te += val;
					}
					else te += val;
				}
				te += endl;
			}
		}
		break;
		
	case "CreateEditXvectorContent":
		{
			let Xvec = inter.edit_Xvector;
			te += "Individual,Value"+endl;
			for(let j = 0; j < Xvec.ind_list.length; j++){
				te += '"'+Xvec.ind_list[j]+'",'+Xvec.X_value[j]+endl;
			}
		}
		break;

	case "CreateEditParamContent":
		{
			let i = inter.edit_param.i;
			
			let par = model.param[i];
			let dep = par.dep;
			let list = par_find_list(par);
			
			let value = inter.edit_param.value;
		
			let dim = get_dimensions(value);
			
			let ele_list = get_element_list(value,dim);
		
			for(let j = 0; j < dep.length; j++) te += dep[j] + ",";
			if(inter.edit_param.type == "PriorSplit") te += "prior";
			else te += "value";
			te += endl;
		
			for(let k = 0; k < ele_list.length; k++){
				for(let j = 0; j < dep.length; j++) te += list[j][ele_list[k][j]]+",";
				te += get_element(value,ele_list[k])+endl;
			}
		}
		break;
	
	case "CreateEditAmatrixContent":
		{
			let A = inter.edit_Amatrix;
			for(let i = 0; i < A.ind_list.length; i++){
				if(i != 0) te += ",";
				te += A.ind_list[i];
			}
			te += endl;
			
			for(let j = 0; j < A.ind_list.length; j++){
				for(let i = 0; i < A.ind_list.length; i++){
					if(i != 0) te += ",";
					te += A.A_value[j][i];
				}
				te += endl;
			}
		}
		break;
		
	default: 
		error("Export table problem");
		break;
	}

	te = add_escape_char(te);
	
	write_file_async(te,filename,"export");
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

	write_file_async(te,filename,"export");
}

/// Writes a file asynchronously
function write_file_async(te,filename,type)
{
	const fs = require('fs');
	
	fs.writeFile(filename, te, function (err) {
		if (err) {
			setTimeout(function(){ alertp("There was an error attempting to save your data.");}, 10);
			return;
		}
		else{
			if(type != "image"){
				inter.saving.num++;
				set_loading_percent(70+30*inter.saving.num/inter.saving.num_to_do);
				if(inter.saving.num == inter.saving.num_to_do){
					saving_done();
				}
			}
		}
	});
}


/// Gets a path
function get_path(file)
{
	let root = "Examples\\";
	if(begin(file,root)) file = "..\\"+file;
	
	if(ver=="windows") file = file.replace(/\//g,"\\");
	return file;
}


/// Loads any files within the BICI file
function load_BICI_files(ans,per_start,per_end)
{
	loading_symbol_message("Loading...");
	
	var fs = nw.require('fs');
	
	let dfl = ans.data_file_list;
	let num = 0, num_max = dfl.length;
	for(let i = 0; i < num_max; i++){
		let file = get_path(dfl[i].full_name);
		if(begin_str(file,"..\\")) file = file.substr(3);
	
		fs.readFile(file, 'utf8', function(err, txt) {
			if(err){ 
				setTimeout(function(){ 
					stop_loading_symbol();
					alert_help("Problem loading file","The following error occurred: "+err);
					generate_screen();
					}, 10);
	
				return;
			}
	
			dfl[i].te = txt;
			num++;
		
			set_loading_percent(per_start+((num+0.5)/num_max)*(per_end-per_start));
		
			if(num == num_max){
				setTimeout(function(){ start_worker("Import2",{ data_file_list:dfl});	}, 10);
			}
		});
	}
}

"use strict";
// Functions which generate the help pop-up box


/// Sets up the background buttons the the help page
function add_help_background_buts(lay)
{
	lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, type:"HelpBackground"});
	lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, ac:"Nothing", type:"Nothing"});
}


/// Close the help page
function close_help()
{
	inter.help = {};
}


/// Sets up the buttos for the help page
function add_help_buts(lay)
{
	let l_store = inter.layer.length;

	let help_type = "normal";
	if(inter.help.script != undefined) help_type = "import";
	if(inter.help.show_datatable == true) help_type = "datatable";
	if(inter.help.warn != undefined) help_type = "warning";
	if(test_comment) help_type = "testing";
	
	let marx = 1;
	let marx2 = 0; 
	switch(help_type){
	case "datatable":  marx2 = 0.5; break;
	case "import": marx2 = 0.8; break;
	case "normal": case "warning": case "testing": break;
	default: error("Option not recognised 42"+help_type); break;
	}
	
	let dx = lay.dx;

	let dx_cont = dx-2*marx-2*marx2-0.3;
	add_layer("HelpContent",lay.x+marx+marx2,0,dx_cont,0,{ymax:20.7, help_type:help_type});	
	
	let hei = inter.layer[l_store].dy;
	let sh = 0.3;
	let hei_import_text = 0;
	
	switch(help_type){
	case "import":
		if(inter.help.scroll_to_line == true){
			inter.help.scroll_to_line = false;
			for(let ll = l_store+1; ll < inter.layer.length; ll++){
				if(inter.layer[ll].name == "Yscroll"){
					let y = inter.help.liney - 1.4*5; if(y < 0) y = 0; 
					if(y > 0){
						change_scroll(y,inter.layer[ll].but[0],"page_set");
						generate_screen(); return;
					}
				}
			}
		}
	
		let text_anno = text_convert_annotation(inter.help.te,1,1.4,dx_cont,"",BLACK);
		
		hei_import_text = text_anno.height;
		hei += hei_import_text+1.5;
		sh = 0.6;
		break;
	
	case "datatable":
		hei += 7; sh = 0.6;
		break;

	case "normal":
		break;

	case "warning":
		break;
		
	case "testing":
		break;
		
	default: error("Option not recognised 43"+help_type); break;
	}
	
	let cy = lay.dy/2 - hei/2;

	let martop = 2.8, marbot = 2.8;
	lay.add_button({x:0, y:cy-martop, dx:dx, dy:hei+martop+marbot, type:"HelpBack"});

	let title = inter.help.title;
	if(test_comment){
		title = all_comments_name[inter.help.index]+" "+inter.help.index+" "+all_comments_name.length;
	}
	
	lay.add_title(title,1.3,cy-martop + 0.9,{col:WHITE, col_line:WHITE});
	
	let marrx = 1.3,  marry = 0.8, busi = 1.2;
	lay.add_button({x:dx-marx-marrx, y:cy-martop+marry, dx:busi, dy:busi, ac:"CloseHelp", type:"WhiteCross"});

	for(let l = l_store; l < inter.layer.length; l++){
		inter.layer[l].y += cy+sh;
		if(inter.layer[l].name == "Xscroll") inter.layer[l].y += inter.layer[l_store].dy;
	}

	switch(help_type){
	case "import":
		{
			let ma = 1.5;
			lay.add_button({x:ma, y:cy+0.4, dx:dx-2*ma, dy:inter.layer[l_store].dy+0.5, type:"Outline", back_col:WHITE, col:WHITE}); 
		
			if(inter.help.te != ""){
				lay.add_paragraph(inter.help.te,dx_cont,1.3,cy+hei+2.8-hei_import_text-2.8,WHITE,para_si,para_lh);
			}
		}
		break;

	case "datatable":
		{
			let ma = 1.5;
			lay.add_button({x:ma, y:cy+0.4, dx:dx-2*ma, dy:inter.layer[l_store].dy+0.5, type:"Outline", back_col:WHITE, col:WHITE}); 
			
			let y = cy+hei-5.5;
		
			lay.add_button({te:"Upload a new data table", x:1.3, y:y, dx:10, dy:1.3, col:WHITE, type:"SubTitle"});
			y += 1.8;
			
			lay.add_paragraph("Table heading:",dx_cont,1.3,y,WHITE,para_si,para_lh);
			y += 0.15;
			lay.add_radio_white(8.5,y,"Yes","Yes",inter.radio_heading,{back_col:BLUE_BACK});
			lay.add_radio_white(12.5,y,"No","No",inter.radio_heading,{back_col:BLUE_BACK});
			y += 1.8;
			lay.add_paragraph("File format:",dx_cont,1.3,y,WHITE,para_si,para_lh);
			y += 0.15;
			lay.add_radio_white(6.8,y,"csv","Comma separated (.csv)",inter.radio_format,{back_col:BLUE_BACK});
			lay.add_radio_white(18.3,y,"tsv","Tab separated (.txt,.tsv)",inter.radio_format,{back_col:BLUE_BACK});
		}
		break;

	case "normal":
		break;
		
	case "warning":
		break;

	case "testing":
		break;
		
	default: error("Option not recognised 44"+help_type); break;
	}
	
	let help = inter.help;

	let widbut = but_width, heibut = but_height;
	
	if(help_type == "warning" && help.ok){
		lay.add_corner_button([["Cancel","White","CloseHelp"],["Run","White",help.ok]],{x:dx-marx+0.2, y:cy+hei+2.5});
	}
	else{
		if(help.save){
			lay.add_corner_button([["Close","White","CloseHelp"],["Save","White",help.save]],{x:dx-marx+0.2, y:cy+hei+2.5});
		}
		else{
			if(help.ok){
				lay.add_corner_button([["Cancel","White","CloseHelp"],["OK","White",help.ok]],{x:dx-marx+0.2, y:cy+hei+2.5});
			}
			else{
				if(help.upload){
					lay.add_corner_button([["Cancel","White","CloseHelp"],["Upload","White",help.upload]],{x:dx-marx+0.2, y:cy+hei+2.5});
				}
				else{
					lay.add_corner_button([["Done","White","CloseHelp"]],{x:dx-marx+0.2, y:cy+hei+2.5});
				}
			}
		}
	}
}


/// Gets the width of each of the lines
function add_script_width(script)
{
	for(let li = 0; li < script.length; li++){
		let text_anno = text_convert_annotation(script[li].te,para_si,para_lh,LARGE,undefined,BLACK);
		script[li].w = text_anno.width_first+1;
	}	
}


/// Sets up the buttons in the scrolable help page
function add_help_content_buts(lay)
{	
	let cx = 0.5, cy = 0;
	
	let type = lay.op.help_type;
	switch(type){
	case "normal": case "warning": case "testing":
		{
			let dx = lay.inner_dx-2*cx;
			let te = inter.help.te;
			if(test_comment) te = all_comments[inter.help.index];

			if(type == "warning"){
				let warn = inter.help.warn;
				for(let j = 0; j < warn.length; j++){
					let wa = warn[j];
					lay.add_button({x:cx, y:cy, dx:1.3*1.5, dy:1.1*1.5, type:"WarnPicBlack"});
					cy = lay.add_paragraph(wa,dx-3,cx+2.3,cy+0.3,WHITE,para_si,para_lh);
					cy += 1;
				}			
			
				let sug = inter.help.sug;
				if(sug.length > 0){
					te = "This may be rectified by one of the following suggestions:\n";
				
					for(let k = 0; k < sug.length; k++) te += "• "+sug[k]+"\n";
				}
			}
			
			let para = te.split("\n");

			for(let j = 0; j < para.length; j++){
				let fir = para[j].substr(0,2);
				switch(fir){
				case "• ":
					lay.add_paragraph("•",dx,cx+0.5,cy,WHITE,para_si,para_lh);
					cy = lay.add_paragraph(para[j].substr(2),dx-1.4,cx+1.4,cy,WHITE,para_si,para_lh);
					break;
					
				case "* ":
					cy = lay.add_paragraph(para[j].substr(2),dx-2.4,cx+2.4,cy,WHITE,para_si,para_lh);
					break;
					
				case ">>":
					cy = lay.add_paragraph(para[j].substr(2),dx-1.4,cx+1.4,cy,WHITE,para_si,para_lh);
					break;
					
				case "]>":
					cy = lay.add_paragraph(para[j].substr(2),dx-2.4,cx+2.4,cy,WHITE,para_si,para_lh);
					break;
					
				default:
					cy = lay.add_paragraph(para[j],dx,cx,cy,WHITE,para_si,para_lh);
					break;
				}
				cy += 0.5;
			}
			cy -= 0.3;
		}
		break;
		
	case "import":
		{
			let script = inter.help.script;
			
			let fo = get_font(para_si);
			
			let tenum_st;
			let tab;
			
			cy = 0;
			for(let j = 0; j < script.length; j++){
				let tenum = String(script[j].line+1).length;
				if(tenum != tenum_st){
					let str = "000000000000".substr(0,tenum);
					tab = text_width(str,fo)+0.5;
					tenum_st = tenum;
				}
				
				lay.add_button({te:script[j].line+1, x:cx, y:cy, dx:tab, dy:para_lh, type:"Text", si:1, font:fo, col:RED}); 
				
				let back_col; if(script[j].line == inter.help.line){ back_col = LLBLUE; inter.help.liney = cy;}
				
				let text_anno = text_convert_annotation(script[j].te,para_si,para_lh,LARGE,undefined,BLACK);
			
				lay.add_button({text_anno:text_anno, x:cx+tab, y:cy, dx:text_anno.wmax, dy:para_lh, back_col:back_col, type:"Script Line"}); 
			
				cy += para_lh;
			}
		}
		break;
	
	case "datatable":
		{
			let table = { width:data.table_width, heading:[{name:""},{name:"File"},{name:"# cols"},{name:"# rows"},{name:""}], content:[]};
	
			if(data.table.length == 0) table = "There are currently no data tables loaded.";

			for(let i = 0; i < data.table.length; i++){
				let tab = data.table[i];
				table.content.push([{te:" "},{te:tab.filename, i:i, ac:"UsePreloadedTable"},{te:tab.ncol},{te:tab.nrow},{te:"Delete",ac:"DeleteSource3",source:data.table,i:i}]);
			}	

			add_table_content(lay,table,WHITE,1.5);
			let box = get_but_box(lay.but);

			cy = box.ymax;
		}
		break;
	
	default: error("Option not recognised 46"); break;
	}
	
	if(cy > lay.op.ymax) cy = lay.op.ymax;

	lay.dy = cy;
}


/// Outputs a help message
function alert_help(title,te)
{	
	if(te == undefined){ te = title; title = "Sorry an error occurred...";}

	if(te.length > 0 && te.substr(te.length-1,1) != "."){
		let ch = te.substr(te.length-1,1);
		if(ch != "." && ch != "\n") te += ".";
	}
	
	while(te.length > 0 && te.substr(te.length-1,1) == "\n"){
		te = te.substr(0,te.length-1);
	}
	
	inter.help = {title:title, te:te};
}

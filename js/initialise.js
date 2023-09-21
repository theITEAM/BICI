"use strict";

/// Intialises the software
function initialise_BICI()                                         
{
	document.onselectstart = function() { return false; };

	if(debug == true) require('nw.gui').Window.get().showDevTools();

	inter.graph = new Graph();

	initialise_pages();

	inter.canvas = ById("myCanvas");
	inter.canvas_cv = inter.canvas.getContext("2d");
	cv = inter.canvas_cv;

	add_listeners();
		
	ById("bod").style.visibility = "visible";	
	
	change_page({pa:"Home"});
	
	if(1 == 0){
		if(make_file == 1){
			load_import_test();
		
			//model.start_new("without results"); TURN ON TO SEE IF MODEL IS REFERNCED IN RESULTS
			
			if(inter.help.title == undefined){
				//change_page({pa:"Simulation", su:"Population"});
				//change_page({pa:"Model", su:"Parameters"});
	
				//change_page({pa:"Simulation", su:"Parameters"});
				//change_page({pa:"Model", su:"Compartments"});
				//change_page({pa:"Inference", su:"ResultsParameters"});
				change_page({pa:"Inference", su:"Run"});
				
			}
			
			//if(inter.help.title == undefined && sim_result.sample.length > 0){
		if(inter.help.title == undefined && (sim_result.sample.length > 0 )){
				//change_page({pa:"Inference", su:"Data"});
				//change_page({pa:"Model", su:"Compartments"});
				//change_page({pa:"Model", su:"Parameters"});
				//change_page({pa:"Simulation", su:"Run", susu:"Transitions"});
				//change_page({pa:"Simulation", su:"Results", susu:"Populations"});
				//change_page({pa:"Simulation", su:"Results", susu:"Splines"});
				//change_page({pa:"Simulation", su:"Run", susu:"Individuals"});
				//change_page({pa:"Simulation", su:"Population"});
					//inter.graph.export_image();
				//change_page({pa:"Simulation", su:"Generate Data"});
			}
			
			if(inter.help.title == undefined &&  inf_result.sample.length > 0){
				//	change_page({pa:"Inference", su:"Results", susu:"Individuals"});
				//change_page({pa:"Inference", su:"Results", susu:"Parameters"});
				//change_page({pa:"Inference", su:"Results", susu:"Populations"});
				//change_page({pa:"Inference", su:"Results", susu:"Populations"});
				//change_page({pa:"Inference", su:"Results", susu:"Transitions"});
			
			}
		}
		else{
			//load_example("temp2");
			load_example("temp");
			//change_page({pa:"Simulation", su:"Run"});
		
			//load_example("SIR");
			//load_example("marek");
		}
	}
}


/// Adds listeners
function add_listeners()
{
	let a = document;
	
	a.addEventListener('contextmenu', function(ev) {
    ev.preventDefault();
		// Right click
    //alert('success!');
    return false;
	}, false);

	a.addEventListener('paste', function(evt) {
		event.preventDefault();
    let paste = (event.clipboardData || window.clipboardData).getData('text');
		cursor_paste(paste);
	}, false);
	
	a.addEventListener('copy', function(evt) {
		event.preventDefault();
		cursor_copy();
	}, false);

	a.addEventListener('cut', function(evt) {
		event.preventDefault();
		cursor_copy(true);
	}, false);

	a.addEventListener('mousemove', function(evt) {
		inter.mouse_pos = mouse_position(inter.canvas, evt);
		
		if(inter.move_timeout != true){
			inter.move_timeout = true;
			inter.move_timeout_num = 0;
		
			setTimeout(function(){
				inter.move_timeout = false;
				mouse_move(inter.mouse_pos.x,inter.mouse_pos.y);
				},100);
		}
		inter.move_timeout_num++;
	}, false);

	a.addEventListener('mousedown', function(evt) {
		let mousePos = mouse_position(inter.canvas, evt);
		mouse_down(mousePos.x,mousePos.y,evt); 
	}, false);

	a.addEventListener ("mouseout", function(evt) {
		generate_screen();
	}, false);
	
	a.addEventListener('mouseup', function(evt) {
		let mousePos = mouse_position(inter.canvas, evt); 
		//ctrlkey = evt.ctrlKey;
		if(evt.altKey) location.reload(true);
		mouse_up(mousePos.x,mousePos.y);
	}, false);
	
	a.onkeydown = key_press;
	
	a.addEventListener('mousewheel', function(evt) {
		let mx = inter.mx, my = inter.my;
		let l;
		for(l = inter.layer.length-1; l >= 0; l--){
			let lay = inter.layer[l];
			
			if(lay.name == "HelpBackground") break;
		
			if(lay.name == "Yscroll"){
				let layc = lay.op.lay;
				if(mx >= layc.x && mx < layc.x+layc.dx && my >= layc.y && my < layc.y+layc.dy){
					change_scroll(evt.deltaY,lay.but[0],"wheel");
				
					let gen_screen = false;
				
					if(bubble_on() == true && find(inter.layer,"name",inter.bubble.lay_name) > l){
						generate_screen();
					}
					else{
						lay.initialise_buttons();
						lay.plot_buttons();
						layc.plot_buttons();
			
						plot_screen();
					}
					return;
				}
			}
		}
	}, false);
}


/// Loads all the pictures used in the graphical interface
function load_pictures()
{
	let file_name = ["logo_big","speech","newmod","manual","paper","warn","desc","SIR","SEIR","grid2","transexp","Simple population","SIbd","LV","gridSI","metapopUK","scotland_region","badger","zoomin","zoomout","desc","settings","SIRsex"];

	inter.npicloaded = 0;
	for(let i = 0; i < file_name.length; i++){
		inter.pics[i] = { name:file_name[i]};
		inter.pics[i].image = new Image; 
		inter.pics[i].image.src = "pics/"+file_name[i]+".png";
		inter.pics[i].image.onload = function(i){  
			inter.npicloaded++;
			if(inter.npicloaded == file_name.length) initialise_BICI(); 
		};
	}
}


/// Finds a pictures from the name
function find_pic(name)
{
	for(let i = 0; i < inter.pics.length; i++){
		if(name == inter.pics[i].name) return inter.pics[i].image;
	}
	
	error(name+" Cannot find picture");
	
	let im = new Image;
	im.src = "pics/warn.png";
	return im;
}




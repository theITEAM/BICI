"use strict";
/// Functions which deal with annotations on the model


/// Adds all the buttons associated with annotation
function add_annotation_buts(lay)
{
	let p = model.get_p();
	let cl = model.get_cl();
	if(model.species.length == 0 || cl == model.species[p].ncla) return;
	let claa = model.species[p].cla[cl];
	
	let cam = claa.camera;    

	for(let k = 0; k < claa.annotation.length; k++){
		let an = claa.annotation[k];
		switch(an.type){
			case "text":
				{
					let x = lay.dx/2 + (an.x-cam.x)*cam.scale; 
					let y = lay.dy/2 + (an.y-cam.y)*cam.scale;	

					let si = an.tesize*cam.scale*cam.ruler;
					let fo = get_font(si);
					let w = text_width(an.te,fo);
					lay.add_button({te:an.te, x:x-w/2, y:y-0.5*si, dx:w+1, dy:si, p:p, cl:cl, i:k, si:si, font:fo, type:"LabelText", col:an.color, ac:"LabelText"});
				}
				break;
			
			case "box":
				{
					let xmin = LARGE, xmax = -LARGE;
					let ymin = LARGE, ymax = -LARGE;
					
					for(let k = 0; k < an.comps.length; k++){
						let na = an.comps[k];
						let c = find(claa.comp,"name",na);
						let co = claa.comp[c];
					
						switch(co.type){
						case "box": 
							{
								let x = co.x-co.w/2; if(x < xmin) xmin = x;
								let y = co.y-co.h/2; if(y < ymin) ymin = y;
						
								x = co.x+co.w/2; if(x > xmax) xmax = x;
								y = co.y+co.h/2; if(y > ymax) ymax = y;
							}
							break;
							
						case "latlng": 
							{
								let x = co.x; if(x < xmin) xmin = x; if(x > xmax) xmax = x;
								let y = co.y; if(y < ymin) ymin = y; if(y > ymax) ymax = y;
							}
							break;
								
						case "boundary": 	
							{
								let x = co.xmid; if(x < xmin) xmin = x; if(x > xmax) xmax = x;
								let y = co.ymid; if(y < ymin) ymin = y; if(y > ymax) ymax = y;
							}
							break;

						default: error("Option not recognised 2"); break;
						}
					}
					
					let mar = cam.ruler;
					if(xmin == xmax || ymin == ymax) mar *= 2;
					
					let xc = xmin - mar;
					let yc = ymin - 1.5*mar;
					let wc = xmax-xmin+2*mar; 
					let hc = ymax-ymin+2.5*mar;
					
					let x = lay.dx/2 + (xc-cam.x)*cam.scale; 
					let y = lay.dy/2 + (yc-cam.y)*cam.scale;	
	
					let w = wc*cam.scale; 
					let h = hc*cam.scale; 
				
					lay.add_button({te:an.te, x:x, y:y, dx:w, dy:h, type:"Box", p:p, cl:cl, i:k, col:an.color, ac:"Box", si:an.tesize*cam.scale*cam.ruler});
				}
				break;
			
			case "map": break;

			default: error("Option not recognised 1"+an.type); break;
		}
	}
}


/// Adds all the buttons associated with annotation maps
function add_annotation_map_buts(lay)
{
	let p = model.get_p();
	let cl = model.get_cl();
	
	if(model.species.length == 0 || cl == model.species[p].ncla) return;

	let claa = model.species[p].cla[cl];
	
	//lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, ac:"ClassificationBack", type:"Nothing", p:p, cl:cl});
		
	let cam = claa.camera;    

	for(let k = 0; k < claa.annotation.length; k++){
		let an = claa.annotation[k];
		switch(an.type){
			case "map":
				{
					for(let i = 0; i < an.feature.length; i++){
						let fea = an.feature[i];
						let box = fea.box;
					
						let p1 = trans_point(box.xmin,box.ymin,cam,lay);
						let p2 = trans_point(box.xmax,box.ymax,cam,lay);
						
						lay.add_button({x:p1.x, y:p1.y, dx:p2.x-p1.x, dy:p2.y-p1.y, type:"Feature", polygon:fea.polygon});
					}
				}
				break;

			default: break;
		}
	}
	
	// Plots labels in latlng view
	for(let k = 0; k < claa.ncomp; k++){
		let c = claa.comp[k];
		if(c.type == "latlng"){
			let p = trans_point(c.x,c.y,cam,lay);
			
			let te = claa.comp[k].name
			let si = cam.scale*cam.ruler;
			let fo = get_font(si);
			let w = text_width(te,fo);
				
			lay.add_button({te:te, x:p.x-w/2, y:p.y+0.5*si, dx:w, dy:si, si:si, font:fo, type:"LatLngLabel", col:BLACK});
		}
	}
}


/// Transforms a lng lat position into an x,y position 
function transform_latlng(lng,lat)
{
	let x = lng*Math.PI/180;
	let sign = 1;
	let phi_p = lat*Math.PI/180; if(phi_p < 0){ phi_p = -phi_p; sign = -1;}
	let y = -sign*Math.log(Math.tan(0.5*phi_p + 0.25*Math.PI));
	
	return {x:x, y:y};
}


/// Transforms a x,y position into an lng lat position 
function transform_latlng_inv(x,y)
{
	let lng = x*180/Math.PI;
	
	let sign = 1;
	if(y < 0){ sign = -1; y = -y;}

	let phi_p = (Math.atan(Math.exp(y)) - 0.25*Math.PI)/0.5;
	
	let lat = -sign*phi_p*180/Math.PI;
	
	return {lng:lng, lat:lat};
}


/// Transforms a point based on camera
function trans_point(x,y,cam,lay)
{
	return {x:lay.dx/2 + (x-cam.x)*cam.scale, y:lay.dy/2 + (y-cam.y)*cam.scale};
}


/// Does the reverse function (i.e. gets position in model based on mouse location
function trans_point_rev(xx,yy,cam,lay)
{
	return {x:cam.x + (xx - lay.dx/2)/cam.scale, y:cam.y + (yy - lay.dy/2)/cam.scale};
}
	
	
/// Loads up a map into the annotations
function load_annotation_map(te,p,cl,def)
{			
	let data = JSON.parse(te);

	let feature = get_feature(data);
	if(feature.length == 0){
		inter.help = {title:"Error loading file", te:"This map has no features"};
		return;
	}	
	let box = get_map_bound_box(feature);
	
	let claa = model.species[p].cla[cl];
	
	claa.annotation.push({type:"map", feature:feature, box:box});
	
	set_camera(p,cl);
}


/// Loads up a map into compartments
function load_compartment_map(te,p,cl,def)
{
	if(p == undefined) p = model.get_p();
	if(cl == undefined) cl = model.get_cl();
				
	let da = JSON.parse(te);

	let feature = get_feature(da);
	if(feature.length == 0){
		inter.help = {title:"Error loading file", te:"This map has no features"};
		return;
	}	
		
	start_data_source("CompMap",{},{p:p, cl:cl});

	let tab = get_feature_table(da,feature);
	data.table.push(tab);
	
	let so = inter.edit_source;
	so.table_loaded = true;
	so.p = p; so.cl = cl;

	so.data_table_use = data.table.length-1;
			
	transfer_column(0);
}


/// Adds a compartment map
function add_compartment_map()
{
	let so = inter.edit_source;
	let p = so.p, cl = so.cl;
	let claa = model.species[p].cla[cl];
	
	let tab = so.table;
	let ele = so.table.ele;
	for(let r = 0; r < tab.nrow; r++){
		let fea = tab.ele[r][0];
		model.add_map_compartment(tab.ele[r][1],p,cl,fea,get_mask(fea),MAP_DEFAULT);
	}
	
	set_camera(p,cl);
	close_data_source();
}


/// Adds an individual compartmental bounardy based on file
function add_individual_compartment_boundary(name,p,cl,file,color)
{	
	let tab = import_geojson(file);
	
	if(tab == undefined) return;
	
	for(let c = 1; c < tab.ncol; c++){
		for(let r = 0; r < tab.nrow; r++){
			if(tab.ele[r][c] == name){
				let fea = tab.ele[r][0];
				model.add_map_compartment(name,p,cl,fea,get_mask(fea),color);
				return "success";
			}
		}
	}
 
	return "'"+name+"' could not be found in the file '"+file+"'";
}


/// Gets a mask of the feature so that mouse over works
function get_mask(fea)
{
	let can = document.createElement('canvas');
  let cv = can.getContext('2d');
	
	can.width = mask_size;
	can.height = mask_size;
	
	let poly = fea.polygon;
	for(let p = 0; p < poly.length; p++){
		let cor = poly[p].cor;
		
		cv.beginPath();
		for(let i = 0; i < cor.length; i++){
			let xx = mask_size*cor[i][0];
			let yy = mask_size*cor[i][1];
			if(i == 0) cv.moveTo(xx,yy);
			else cv.lineTo(xx,yy);
		}
		cv.closePath();
		cv.fillStyle = WHITE;
		cv.fill();
	}

	let data = cv.getImageData(0,0,mask_size,mask_size).data;
	let mask = [];
	for(let i = 0; i < mask_size*mask_size; i++){
		if(data[4*i+3] < 128) mask.push(0); else mask.push(1);
	}
	
	return mask;
}


/// Works out if mouse is over the compartment map
function mouse_over_comp_map(frx,fry,mask)
{
	if(mask[Math.floor(fry*mask_size)*mask_size+Math.floor(frx*mask_size)] == 1){
		return true;
	}		

	return false;	
}


/// Sets the camera to fit in all the objects 
function set_camera(p,cl)
{
	let claa = model.species[p].cla[cl];
	let cam = claa.camera;
	
	let box = get_model_box(p,cl); if(box.unset) return;
	
	cam.x = 0.5*(box.xmin+box.xmax); cam.y = 0.5*(box.ymin+box.ymax);
	cam.scale = box.scale;
	
	set_ruler(p,cl);
}


/// Sets the camera ruler (this is used to detemine text size in map view
function set_ruler(p,cl)
{
	let box = get_model_box(p,cl);
	
	let claa = model.species[p].cla[cl];
	let cam = claa.camera;
	switch(cam.coord){
	case "cartesian":
		cam.ruler = 1; 
		break;
		
	case "latlng":
		{
			cam.ruler = 1/box.scale; 
			
			let xmin = LARGE, xmax = -LARGE;
			let ymin = LARGE, ymax = -LARGE;
			for(let c = 0; c < claa.ncomp; c++){
				let co = claa.comp[c];
	
				if(co.type == "boundary"){
					if(co.xmid < xmin) xmin = co.xmid; 
					if(co.xmid > xmax) xmax = co.xmid; 
					if(co.ymid < ymin) ymin = co.ymid; 
					if(co.ymid > ymax) ymax = co.ymid; 
				}
				else{
					if(co.x < xmin) xmin = co.x; 
					if(co.x > xmax) xmax = co.x;
					if(co.y < ymin) ymin = co.y; 
					if(co.y > ymax) ymax = co.y;
				}
			}
			
			let N = 20;
			let grid=[];
			for(let k = 0; k < N*N; k++) grid[k]=[];
			
			for(let c = 0; c < claa.ncomp; c++){
				let co = claa.comp[c];
				let x, y;
				if(co.type == "boundary"){ x = co.xmid, y = co.ymid;}
				else{ x = co.x; y = co.y;}
				
				let i = Math.floor(ALMOST_ONE*N*(x-xmin)/(xmax-xmin+TINY));
				let j = Math.floor(ALMOST_ONE*N*(y-ymin)/(ymax-ymin+TINY));
				grid[j*N+i].push(c);
			}
			
			let av = 0, nav = 0;
			for(let k = 0; k < N*N; k++){
				let gr = grid[k];
				if(gr.length > 2){
					for(let i = 0; i < gr.length; i++){
						let min = LARGE;
						let coi = claa.comp[gr[i]];
						for(let j = 0; j < gr.length; j++){
							if(j != i){
								let coj = claa.comp[gr[j]];
								let d = (coi.x-coj.x)*(coi.x-coj.x) + (coi.y-coj.y)*(coi.y-coj.y);
								if(d < min) min = d;
							}
						}
						av += Math.sqrt(min); nav++;
					}
				}
			}
			
			if(nav > 0){
				av /= nav;
				let rule = 0.15*av;
				if(rule < cam.ruler) cam.ruler = rule;	
			}
			
			let dx = xmax-xmin, dy = ymax-ymin;
			let d = dx; if(dy > d) d = dy;
			
			if(d == 0){
				if(cam.ruler < 0.05) cam.ruler = 0.05;
			}
			else{
				if(cam.ruler > 0.1*d) cam.ruler = 0.1*d;
			}
		}
	
		break;

	default: error("Option not recognised 3"); break;
	}
}


// Gets a box surrounding the complete model
function get_model_box(p,cl)
{
	let claa = model.species[p].cla[cl];
	let cam = claa.camera;
	
	let xmin = LARGE, xmax = -LARGE;
	let ymin = LARGE, ymax = -LARGE;

	for(let c = 0; c < claa.ncomp; c++){
		let co = claa.comp[c];

		switch(co.type){
		case "boundary":
			let box = co.feature.box;
				
			if(box.xmin < xmin) xmin = box.xmin; 
			if(box.xmax > xmax) xmax = box.xmax; 
			if(box.ymin < ymin) ymin = box.ymin; 
			if(box.ymax > ymax) ymax = box.ymax; 
			break;
			
		case "box": case "latlng":
			if(co.x < xmin) xmin = co.x; 
			if(co.x > xmax) xmax = co.x;
			if(co.y < ymin) ymin = co.y; 
			if(co.y > ymax) ymax = co.y;
			break;

		default: error("Option not recognised 4"); break;
		}
	}
			
	for(let i = 0; i < claa.ntra; i++){
		let tr = claa.tra[i];
		for(let j = 0; j < tr.midp.length; j++){
			let po = tr.midp[j];
			if(po.x < xmin) xmin = po.x;
			if(po.x > xmax) xmax = po.x;
			if(po.y < ymin) ymin = po.y;
			if(po.y > ymax) ymax = po.y;			
		}
	}
	
	for(let i = 0; i < claa.annotation.length; i++){
		let anno = claa.annotation[i];

		if(!(anno.type == "map" && anno.feature.length == 248)){		
			switch(anno.type){
			case "label":
				if(anno.x < xmin) xmin = anno.x; 
				if(anno.x > xmax) xmax = anno.x;
				if(anno.y < ymin) ymin = anno.y; 
				if(anno.y > ymax) ymax = anno.y;
				break;
				
			case "map":
				let box = anno.box;
				if(box.xmin < xmin) xmin = box.xmin;
				if(box.xmax > xmax) xmax = box.xmax;
				if(box.ymin < ymin) ymin = box.ymin;
				if(box.ymax > ymax) ymax = box.ymax;
				break;

			default: error("Option not recognised 5"); break;
			}		
		}
	}
	
	let scalex = 0.7*(page_char_wid-menu_width)/(xmax-xmin); 
	if(cam.coord == "cartesian" && scalex > 1) scalex = 1;
	
	let scaley = 0.7*page_char_hei/(ymax-ymin);
	if(cam.coord == "cartesian" && scaley > 1) scaley = 1;

	let scale;
	if(scalex < scaley) scale = scalex;
	else scale = scaley;
	
	let unset = false; if(xmin == LARGE) unset = true;
	
	return {unset:unset, xmin:xmin, xmax:xmax, ymin:ymin, ymax:ymax, scale:scale};
}


/// Gets the bounding box for a series of polygons
function get_map_bound_box(feature)
{
	let xmin = LARGE, xmax = -LARGE;
	let ymin = LARGE, ymax = -LARGE;
	for(let i = 0; i < feature.length; i++){
		let box = feature[i].box;
		if(box.xmin < xmin) xmin = box.xmin;
		if(box.xmax > xmax) xmax = box.xmax;
		if(box.ymin < ymin) ymin = box.ymin;
		if(box.ymax > ymax) ymax = box.ymax;
	}
	return {xmin:xmin, xmax:xmax, ymin:ymin, ymax:ymax}; 
}

	
/// Gets polygons from the data
function get_feature(data) 
{
	let feature = [];
	
	for(let i = 0; i < data.features.length; i++){
		let geo = data.features[i].geometry;

		let polygon = [];
	
		let cor = geo.coordinates;
		
		switch(geo.type){
		case "Polygon":	
			for(let j = 0; j < cor.length; j++){
				polygon.push({cor:cor[j], foreground:true}); 
			}
			break;
			
		case "MultiPolygon":
			for(let pn = 0; pn < cor.length; pn++){
				let cor2 = cor[pn];
				for(let j = 0; j < cor2.length; j++){
					polygon.push({cor:cor2[j], foreground:true}); 
				}
			}
			break;
			
		default: error(geo[j].type+" Type not syported"); break;
		}
	
		/// Transforms lng lat to x y 
		for(let i = 0; i < polygon.length; i++){
			let points = polygon[i].cor;
			for(let j = 0; j < points.length; j++){
				let p = points[j];
				let pt = transform_latlng(p[0],p[1]);
				p[0] = pt.x; p[1] = pt.y; 
			}
		}

		// Gets bounding box
		let xmin = LARGE, xmax = -LARGE;
		let ymin = LARGE, ymax = -LARGE;
		
		for(let i = 0; i < polygon.length; i++){
			let points = polygon[i].cor;
			for(let j = 0; j < points.length; j++){
				let p = points[j];
				if(p[0] < xmin) xmin = p[0];
				if(p[0] > xmax) xmax = p[0];
				if(p[1] < ymin) ymin = p[1];
				if(p[1] > ymax) ymax = p[1];
			} 
		}
			
		// Scales points within the bounding box
		if(xmin != xmax && ymin != ymax){	
			for(let i = 0; i < polygon.length; i++){
				let points = polygon[i].cor;
				for(let j = 0; j < points.length; j++){
					let p = points[j];
					p[0] = (p[0] - xmin)/(xmax - xmin);
					p[1] = (p[1] - ymin)/(ymax - ymin);
				} 
			}
		}
		
		let box = {xmin:xmin, xmax:xmax, ymin:ymin, ymax:ymax};
		feature.push({ name:"", polygon:polygon, box:box});
	}
	
	return feature;
}


/// Gets names for all the features
function get_feature_table(da,feature) 
{
	let nrow = da.features.length;
	let geo = da.features[0].properties;
	
	let keys = Object.keys(geo);

	let ele=[];
	let heading=[];
	heading.push("Boundary");
	for(let r = 0; r < nrow; r++){		
		ele[r]=[];
		ele[r].push(feature[r]);
	}
		
	for(let i = 0; i < keys.length; i++){
		let prop = keys[i];
		let r = 0; while(r < nrow && da.features[r].properties[prop] != undefined) r++;
		if(r == nrow){
			heading.push(prop);
			for(let r = 0; r < nrow; r++){
				ele[r].push(da.features[r].properties[prop]);
			}
		}
	}

	return {filename:data.filename, heading:heading, col_used:[], ele:ele, ncol:heading.length, nrow:nrow, edit:false};
}

	
/// Sets up removal of annotation bubble
function remove_annotation_init()
{
	let p = model.get_p();
	let cl = model.get_cl();
	let anno = model.species[p].cla[cl].annotation;
	let list = [];
	let wmax = 0;
	
	for(let i = 0; i < anno.length; i++){			
		let an = anno[i];
		let te;
		switch(an.type){
		case "map": te = "Map with "+an.feature.length+" features";	break;
		case "text": te = "Label '"+an.te+"'"; break;
		case "box": te = "Box with text '"+an.te+"'"; break;
		default: error("Annotation type not recognised"); break;
		}
	
		let w = text_width(te,get_font(si_radio));
		if(w > wmax) wmax = w;
		list.push({te:te, box:{check:false}});
	}
	wmax += 3.5;
	if(wmax < 11.4) wmax = 11.4;
	
	inter.bubble.rem_annotation = { w:wmax, p:p, cl:cl, list:list};
	change_bubble_mode("remove");
}

	
	
/// Adds a scrollable list to allow for removal of annotations
function annotation_scrollable(lay)
{
	let list = lay.op.list;

	let cy = 0;
	for(let i = 0; i < list.length; i++){
		lay.add_checkbox(0,cy,list[i].te,list[i].te,list[i].box);
		cy += 1.2;
	}
	return cy;
}


/// Deletes selected annotations
function delete_annotation()
{
	let rem = inter.bubble.rem_annotation;
	let anno = model.species[rem.p].cla[rem.cl].annotation;
	
	for(let j = anno.length-1; j >= 0; j--){
		if(rem.list[j].box.check == true) anno.splice(j,1);
	}		
}


/// Zooms in and out
function zoom_factor(p,cl,fac)
{
	let show = model.get_show_model();
			
	let cam;

	if(show == true) cam = model.species[p].cla[cl].camera;
	else cam = inter.graph.op.species[p].cla[cl].camera;
	
	cam.scale *= fac;
	
	if(show == true) clone_camera(p,cl);
}


/// Zooms into a location selected by the mouse
function zoom_double_click(p,cl)
{
	let mx = inter.mx, my = inter.my;
	let lay = get_lay("Main");
	
	let show = model.get_show_model();
	
	let cam;
	if(show == true) cam = model.species[p].cla[cl].camera;
	else cam = inter.graph.op.species[p].cla[cl].camera;
	
	let x = mx-lay.x, y = my-lay.y;
	
	let f = 0.25;
	let xx = f*x + (1-f)*lay.dx/2;
	let yy = f*y + (1-f)*lay.dy/2;
	
	let rev = trans_point_rev(xx,yy,cam,lay);
	cam.x = rev.x; 
	cam.y = rev.y; 
	
	zoom_factor(p,cl,1.3);
}


/// Loads map
function load_map(val,p,cl)                                   
{	 
	start_loading_symbol(0);
	inter.load_map = {p:p, cl:cl};
	
	let file = "D:/BICI_nolongerused/Maps/"+val+".json"; // TO DO turn off	
	//let file = "Maps/"+val+".json";

	//let xmlhttp = new XMLHttpRequest();  
	xmlhttp.abort();
	
	xmlhttp.addEventListener("progress", function(ev) { 
		inter.loading_symbol.percent = Math.floor(100*ev.loaded/ev.total);
	});
	
	xmlhttp.addEventListener("loadend",function(){
		stop_loading_symbol();
		let p = inter.load_map.p;
		let sp = model.species[p];
		let cl = inter.load_map.cl;
	
		load_annotation_map(xmlhttp.responseText,p,cl);
		
		let cam = sp.cla[cl].camera;
		let po = transform_latlng(0,40);
		cam.x = po.x;
		cam.y = po.y;
		cam.scale = (page_char_wid-menu_width)/(2*Math.PI);
		generate_screen();
	});
	
	xmlhttp.open("GET",file,true); xmlhttp.send();
}


/// Adds a layer which allows the user to select compartments
function add_select_buts(lay)
{
	lay.add_button({lay_x:lay.x, lay_y:lay.y, x:0, y:0, dx:lay.dx, dy:lay.dy, ac:"Select", type:"Select"});
}


/// Adds a box to the model
function add_box()
{
	let p = model.get_p();
	let cl = model.get_cl();
	
	let claa = model.species[p].cla[cl];
	let comps = [];
	let lay = get_lay("Compartment");

	let mo = inter.mode;
	let x1 = mo.x1, x2 = mo.x2, y1 = mo.y1, y2 = mo.y2;
	if(x1 > x2){ let t = x1; x1 = x2; x2 = t;}
	if(y1 > y2){ let t = y1; y1 = y2; y2 = t;}
	
	for(let i = 0; i < lay.but.length; i++){
		let bu = lay.but[i];
	
		switch(bu.type){
		case "Compartment": case "CompLatLng":
			let x = bu.x+bu.dx; 
			let y = bu.y+bu.dy;
			if(x > x1 && x < x2 && y > y1 && y < y2){
				comps.push(bu.te);
			}
			break;
		}
	}

	claa.annotation.push({type:"box", te:"", tesize:si_annotation, comps:comps, color:BLACK});
}


/// Creates a GEOJSON feature from a feature (used to create files from loaded boundaries)
function generate_JSON_feature(name,feature)
{
	let cor = [];
	let poly = feature.polygon;
	let box = feature.box;
	for(let i = 0; i < poly.length; i++){
		cor[i] = [];
		cor[i].push(copy(poly[i].cor));
	}
	
	for(let i = 0; i < cor.length; i++){
		for(let j = 0; j < cor[i].length; j++){
			for(let k = 0; k < cor[i][j].length; k++){
				let x = box.xmin + cor[i][j][k][0]*(box.xmax-box.xmin);
				let y = box.ymin + cor[i][j][k][1]*(box.ymax-box.ymin);
				let p = transform_latlng_inv(x,y);
				cor[i][j][k][0] = p.lng;
				cor[i][j][k][1] = p.lat;
			}
		}
	}
	
	let fea = { type:"Feature", geometry:{type:"MultiPolygon", coordinates:cor}, properties:{name:name}};

	return fea;
}
